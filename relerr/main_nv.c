/*
 * Copyright 1991 Michael S. Warren and John K. Salmon.  All Rights Reserved.
 */

#include <stdio.h>		/* only use sprintf */
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stddef.h>
#include <unistd.h>
#include <errno.h>
#include "fastflpt.h"
#include "Assert.h"
#include "SDF.h"
#include "protos.h"
#include "macr.h"
#include "malloc.h"
#include "bigmalloc.h"
#include "SDFwrite.h"
#include "SDFread.h"
#include "physics_n.h"
#include "vop.h"
#include "Msgs.h"
#include "tree.h"
#include "timers.h"
#include "pqsort.h"
#include "singlio.h"
#include "mpmy.h"
#include "mpmy_io.h"
#include "mpmy_abnormal.h"
#include "gc.h"
#include "files.h"
#include "getparam.h"
#include "verify.h"
#include "randoms.h"
#include "ring.h"
#include "decomp.h"
#include "image.h"
#include "memfile.h"

void update_point_mass(body *btab, int nobj, body *p, float eps, float newt);
static void read_point_mass(body *point_mass, SDF *csdfp);
static void tidal_init(body **btab, int *nobj, int *gnobj, SDF *csdfp);
static void collision_init(body **btab, int *nobj, int *gnobj, SDF *csdfp);

static void SanityCheck(bodyptr btab, int nobj, int gnobj, double *mtotp);
static void ConvertVPofT(body *xptr, int n, float dp_dt, int to_physical);
static void SetLBTarget(sortresult_t *decompp, int hetero_load_balance);
static void IntegratePofT(body *xptr, const int n, const float dp,
		   float *tpos, float *tvel,
		   double *kep, double *pep);
static void Integrate(body *xptr, const int n, const float dt,
		      float *tpos, float *tvel,
		      double *kep, double *pep);
static void IntegratePofT_out(const body *xptr, outbody *yptr, const int n, 
			      const float dp,  float *tpos, float *tvel,
			      double *kep, double *pep);
static void Integrate_out(const body *xptr, outbody *yptr, const int n, 
			  const float dt, float *tpos, float *tvel,
			  double *kep, double *pep);
static void set_vels(body *p, int n, float real_time);
static SDF *startup(int argc, char **argv);
static void Periodic(tree_t *tp, float size);
static void WrapPeriodic(body *bp, int n, float *rmin, float *rmax, float sz,
	    int cosmology, int log_time, float tpos);
static void FixCube(body *b, int nobj, float l, float gm);
static void FixGlobalForce(body *bp, int n);
static void acc_zero(body *btab, int nobj, float mtot);
static int maxmem(void);
static int maxheap(void);

#ifdef _AIX
static void PrintLoad(Timer_t *cpu, Timer_t *wc);
#endif

float Znow(float time);
float Hnow(float time);

Timer_t StepTot, StepTotWC, BuildTot;
Timer_t FindForcesTm;
Timer_t PeriodicForceTm, FixCubeTm;
Counter_t NbodyCnt;
Counter_t MemCnt;
Counter_t HeapCnt_;	/* HeapCnt is in the SunOS name space?! */
Counter_t KNtermsCnt;

Timer_t WITm,WTermTm,WNTTm,WITm,PerTm;

static int CaughtTerm;

#ifdef __PARAGON__
/* Is this useful anywhere else??  It should really be NX-specific... */
#include <signal.h>
static void term_handler(int sig)
{
    singlPrintf("Caught TERM signal.\n");
    CaughtTerm = 1;
    signal(SIGTERM, SIG_DFL);
}
#endif

/* Hide the cosmological parameters in here.
   Keep them self-consistent... */
struct cosmo_s{
    float t;
    float a;
    float H0;
    float Omega0;
    float Lambda;
    float GNewt;
    float Zel_f;		/* the 'f' factor for linearly growing modes,
				 used only in set_vel = 1/H*Ddot/D.  It's
				 very close to 1 (exactly?) for flat models. */
} cosmo;

void CosmoPush(struct cosmo_s *p, float time);

void
main(int argc, char *argv[])
{
    int gnobj, nobj;
    bodyptr btab;
    float eps;			/* Plummer smoothing length */
    float tol;			/* MAC tolerance */
		/* for big MAC, this is multiplied by M/(rsize*rsize) */
    int i;
    float rmin[NDIM], rmax[NDIM];
    float sysradius;
    float dt;
    int nsteps;
    int first_step = 1;
    int do_output;
    int output_freq;
    int timer_freq;
    int image_freq;
    int x_pixels;
    int y_pixels;
    int log_image;
    int hetero_load_balance;
    float sort_tol;
    int iter;
    bodyptr p;
    float this_tol, this_eps;
    float frac_tol;
    int log_time = 0;		/* if true, use dt \propto t */
    int comov_eps = 0;		/* if true, use comoving epsilon*/
    int setpvel = 0;
    char outnamebase[256];
    SDF *csdfp;			/* SDF pointer to control file */
    SDF *sdfp;
    float tpos;			/* time positions are at */
    float tvel;
    float tposlast;
    int cosmology = 0;
    int save_first;		/* save first step (for acc testing) */
    double force[NDIM];		/* Things accumulated over all particles */
    double com[NDIM], comv[NDIM];
    double acc2;		/* must be double precision */
    double pe, ke;
    double mtot;
    MPMY_Comm_request req;
    sortresult_t sortedbtab;
    tree_t thetree;
    int massconf, xconf, yconf, zconf;
    int vxconf, vyconf, vzconf;
    int identconf, idconf;
    char name[256];
    int do_restrictvol;
    int read_nfiles, write_nfiles, do_sortoutput;
    int do_NlgN, do_nsquared;
    int vollist[128];
    int nvol;
    double ntermslocal;
    body point_mass;
    int do_point_mass = 0;
    int setup_tidal = 0;
    int setup_collision;
#ifdef __PARAGON__
    int catch_term;
#endif
#if 0
    void *decomptab = 0;
#endif
    int decomp_freq, decomp_iter = 0;
    int do_periodic;
    inherit_t inherit;
    macv_t mac;
    float R0;
    int timeout;
    int set_id;

    MPMY_Init(&argc, &argv);
    singlPrintf("Welcome to the variable O() integrator running on %d procs\n",
		MPMY_Nproc());
#if 0 && defined(__DELTA__)
    free(malloc(11000000));
#endif
    csdfp = startup(argc, argv);
    SDFgetintOrDefault(csdfp, "timeout", &timeout, 600);
    if (timeout > 0) MPMY_TimeoutSet(timeout);
    SDFgetstringOrDefault(csdfp, "datafile", name, sizeof(name), "");
    SDFgetintOrDefault(csdfp, "do_periodic", &do_periodic, 0);
    SDFgetintOrDefault(csdfp, "cosmology", &cosmology, 0);
    SDFgetintOrDefault(csdfp, "set_id", &set_id, 0);
    SDFgetintOrDefault(csdfp, "read_nfiles", &read_nfiles, 0);
    SDFgetintOrDefault(csdfp, "setpvel", &setpvel, 0);
    /* SDFsetbufsz(65536); */
    if( strlen(name) != 0 ) {
	singlPrintf("Reading \"%s\"\n", name);
	if (read_nfiles) MPMY_Nfileio(1);
	sdfp = SDFread(csdfp, (void **)&btab, &gnobj, &nobj, sizeof(body),
		      "mass", offsetof(body, mass), &massconf,
		      "x", offsetof(body, pos[0]), &xconf,
		      "y", offsetof(body, pos[1]), &yconf,
		      "z", offsetof(body, pos[2]), &zconf,
		      "vx", offsetof(body, vel[0]), &vxconf,
		      "vy", offsetof(body, vel[1]), &vyconf,
		      "vz", offsetof(body, vel[2]), &vzconf,
		      "ident", offsetof(body, ident), &identconf,
		      "id", offsetof(body, ident), &idconf,
		      NULL);
	if (read_nfiles) MPMY_Nfileio(0);
	Msgf(("Data read, nobj=%d, gnobj=%d\n", nobj, gnobj));
	Msgf(("Nproc:%d, Procnum: %d, Doc: %d\n",
	       MPMY_Nproc(), MPMY_Procnum(), ilog2(MPMY_Nproc())));
	if( identconf && idconf ){
	    SinglError("You can't have both an 'id' and an 'ident' in the data!\n");
	}
	if( massconf==0 || xconf==0 || yconf==0 || zconf==0 ){
	    SinglError("Could not find %s %s %s %s in data file!\n",
		  (massconf==0)? "mass" : "",
		  (xconf==0)? "x" : "",
		  (yconf==0)? "y" : "",
		  (zconf==0)? "z" : "");
	}
	if( vxconf != vyconf || vxconf != vzconf ){
	    if (setpvel) SinglError("Missing velocity components!\n");
	}
	if( (identconf == 0 && idconf == 0) || set_id){
	    int ni;
	    SinglWarning("No \"ident\" in file, numbering sequentially\n");
	    FixId(btab, nobj, gnobj);
	    /* decomp.c is currently broken unless this is done */
	    ni = ilog2(gnobj);
	    for(i=0; i<nobj; i++)
	      btab[i].ident <<= 31-ni;
	}
	/* With relerr MAC acc initialziation, nterms from file is no help */
	FixNterms(btab, nobj);
	SDFgetfloatOrDefault(sdfp, "Gnewt", &cosmo.GNewt, (float)1.0);
	if( SDFhasname("time", sdfp) )
	  SDFgetfloatOrDefault(sdfp, "time",  &tpos, (float)0.0);
	else
	  SDFgetfloatOrDefault(sdfp, "tpos",  &tpos, (float)0.0);

	if (cosmology) {
	    float Z;
	    cosmo.t = tpos;
	    SDFgetfloatOrDefault(sdfp, "Omega0",  &cosmo.Omega0, (float)1.0);
	    SDFgetfloatOrDefault(sdfp, "Lambda_prime",  &cosmo.Lambda, (float)0.0);
	    /* default is for h_100 = 0.5 */
	    SDFgetfloatOrDefault(sdfp, "H0",  &cosmo.H0, (float)0.0511365);
	    if( SDFhasname("box_size", sdfp) ) {
		SDFgetfloatOrDie(sdfp, "box_size",  &R0);
		R0 /= 2.0;
	    } else
	      SDFgetfloatOrDie(sdfp, "R0",  &R0);

	    /* Now we need to get initial values for cosmo.a */
	    if( SDFhasname("redshift", sdfp) ){
		SDFgetfloat(sdfp, "redshift", &Z);
		cosmo.a = 1.F/(1.F + Z);
	    }else{
		if( cosmo.Omega0 == 1.0F ){
		    cosmo.a = pow( 1.5*cosmo.t*cosmo.H0, 2./3.);
		}else{
		    SinglError("Sorry.  Tell me the redshift in the data file\n");
		}
	    }
	    /* The Zel'dovich 'f' factor is only needed for setting initial
	       velocities.  At this point, we don't know if we will be asked
	       to do setpvel, though, so we read it anyway. */
	    if( SDFhasname("velocity_fac", sdfp) ){
		SDFgetfloatOrDie(sdfp, "velocity_fac", &cosmo.Zel_f);
	    }else{
		cosmo.Zel_f = 1.F;
	    }
	    singlPrintf("float H0 = %g\n", cosmo.H0);
	    singlPrintf("float Omega0 = %g\n", cosmo.Omega0);
	    singlPrintf("float redshift = %g\n", Znow(tpos));
	    singlPrintf("float Lambda = %g\n", cosmo.Lambda);
	    singlPrintf("float Zel_f = %g\n", cosmo.Zel_f);
	}

    /* This doesn't work if there is roundoff error in tvel */
    /* SDFgetfloatOrDefault(sdfp, "tvel",  &tvel, tpos);*/
    /* Fix the tpos == tvel line in Integrate() if you don't like next line */
	tvel = tpos;

	SDFgetintOrDefault  (sdfp, "iter",  &iter, 0);
	SDFgetintOrDefault(csdfp, "do_point_mass", &do_point_mass, 0);
	if (do_point_mass) {
	    if (iter && !tidal_init)
	      read_point_mass(&point_mass, sdfp);
	    else
	      read_point_mass(&point_mass, csdfp);
	}
	if(sdfp)
	  SDFclose(sdfp);
    } else {
	int seed, cencon, start;
	ran_state ranstate;
	float rsq;

	singlPrintf("Generating random dataset\n");
	if( SDFgetint(csdfp, "nobj", &gnobj) )
	  SinglError("Sorry, you've got to have an \"nobj\"\n");
	SDFgetintOrDefault(csdfp, "seed", &seed, 123);
	SDFgetintOrDefault(csdfp, "cencon", &cencon, 0);
	singlPrintf("int seed = %d;\n", seed);
	singlPrintf("int cencon = %d;\n", cencon);

	NobjInitial(gnobj, MPMY_Nproc(), MPMY_Procnum(), &nobj, &start);
	btab = Calloc(nobj, sizeof(body));
	FixId(btab, nobj, gnobj);
	FixNterms(btab, nobj);
	ran_init(seed+(MPMY_Procnum()+1), &ranstate);
	for (p = &btab[0]; p < &btab[nobj]; p++) {
#ifdef __PARAGON__
	    clear_tregs();
#endif
	    p->mass = 1.0 / gnobj;		 /*   set masses equal */
	    if (do_periodic)
	      rsq = cube_rand(&ranstate, NDIM, p->pos);
	    else
	      rsq = sphere_rand(&ranstate, NDIM, p->pos);
	    sphere_rand(&ranstate, NDIM, p->vel);
	    if( cencon ){
		float scale;
		scale = uniform_rand(&ranstate)*recipsqrtf(rsq);
		VS(p->pos, *= scale);
	    }
	}
	if (do_periodic) R0 = 1.0;
	cosmo.GNewt = (float)1.0;
	tvel = tpos = (float)0.0;
	iter = 0;
    }
    singlPrintf("Maxmem after data read is %d (%d)\n", maxmem(), maxheap());
    if( Msg_test("memleak") ){
	Msg_do("Memory map after data read\n");
	malloc_print();
    }

    SDFgetfloatOrDie(csdfp, "epsilon", &eps);
    SDFgetfloatOrDie(csdfp, "errtol", &tol);
    SDFgetfloatOrDefault(csdfp, "frac_tol", &frac_tol, 0.0);
    SDFgetfloatOrDie(csdfp, "dt", &dt);
    SDFgetintOrDie(csdfp, "nsteps", &nsteps);
    SDFgetintOrDefault(csdfp, "log_time", &log_time, 0);
    SDFgetintOrDefault(csdfp, "comov_eps", &comov_eps, 0);
    SDFgetintOrDefault(csdfp, "save_first", &save_first, 0);
    SDFgetintOrDefault(csdfp, "do_restrictvol", &do_restrictvol, 0);
    SDFgetintOrDefault(csdfp, "write_nfiles", &write_nfiles, 0);
    SDFgetintOrDefault(csdfp, "do_sortoutput", &do_sortoutput, !write_nfiles);
    SDFgetintOrDefault(csdfp, "do_NlgN", &do_NlgN, 0);
    SDFgetintOrDefault(csdfp, "do_nsquared", &do_nsquared, 0);
    SDFgetintOrDefault(csdfp, "setup_tidal", &setup_tidal, 0);
    SDFgetintOrDefault(csdfp, "setup_collision", &setup_collision, 0);
    if (do_restrictvol) {
	int ret;
	SDFgetintOrDie(csdfp, "nvol", &nvol);
	ret = SDFseekrdvecs(csdfp, "vollist", 0, nvol, vollist, sizeof(int),
			    NULL);
	if (ret) 
	  Error("SDFseekrdvecs failed, %s\n", SDFerrstring);
    }
    if (setup_tidal) {
	tidal_init(&btab, &nobj, &gnobj, csdfp);
    }
    if (setup_collision) {
	collision_init(&btab, &nobj, &gnobj, csdfp);
	tpos = tvel = 0.0;
	iter = 0;
    }

    if( SDFgetstring(csdfp, "outfile", outnamebase, sizeof(outnamebase)) == 0){
	do_output = ( strlen(outnamebase) > 0 );
    }else{
	do_output = 0;
    }
    if( do_output ){
	SDFgetintOrDefault(csdfp, "output_freq", &output_freq, nsteps);
    }else{
	output_freq = 1;
    }
    SDFgetintOrDefault(csdfp, "timer_freq", &timer_freq, output_freq);
    SDFgetintOrDefault(csdfp, "image_freq", &image_freq, 0);
    SDFgetintOrDefault(csdfp, "decomp_freq", &decomp_freq, 0);
    SDFgetintOrDefault(csdfp, "x_pixels", &x_pixels, 512);
    SDFgetintOrDefault(csdfp, "y_pixels", &y_pixels, 512);
    SDFgetintOrDefault(csdfp, "log_image", &log_image, 0);
    SDFgetintOrDefault(csdfp, "hetero_load_balance", &hetero_load_balance, 0);
    SDFgetfloatOrDefault(csdfp, "sort_tol", &sort_tol, 0.01);

    if(csdfp) 
	SDFclose(csdfp);

#ifdef __PARAGON__
    SDFgetintOrDefault(csdfp, "catch_term", &catch_term, 0);
    if (catch_term){
	singlPrintf("Catching SIGTERM (usually sent by NQS)\n");
	signal(SIGTERM, term_handler);
    }
#endif

    if (do_periodic) {
	EnableTimer(&PeriodicForceTm, "Periodic F");
	EnableTimer(&FixCubeTm, "Fix Cube");
    }

    singlPrintf("float errtol = %g;\n", tol);
    singlPrintf("float dt = %g;\n", dt);
    singlPrintf("float epsilon = %g;\n", eps);
    singlPrintf("int iter = %d;\n", iter);
    singlPrintf("int nsteps = %d;\n", nsteps);
    singlPrintf("int nproc = %d;\n", MPMY_Nproc());
    singlPrintf("int do_NlgN = %d;\n", do_NlgN);
    if( do_output ){
	singlPrintf("Output to %s.nnnn, every %d steps\n", 
	       outnamebase, output_freq);
    }else{
	singlPrintf("No output.\n");
    }
    singlPrintf("int timer_freq = %d;\n", timer_freq);
    singlPrintf("int image_freq = %d;\n", image_freq);
    singlPrintf("int decomp_freq = %d;\n", decomp_freq);
    singlPrintf("int hetero_load_balance = %d;\n", hetero_load_balance);
    singlPrintf("float sort_tol = %.4f;\n", sort_tol);
    singlPrintf("int do_periodic = %d;\n", do_periodic);
    if (cosmology) {
	singlPrintf("int cosmology = %d;\n", cosmology);
	singlPrintf("int log_time = %d;\n", log_time);
	singlPrintf("int comov_eps = %d;\n", comov_eps);
	singlPrintf("int setpvel = %d;\n", setpvel);
	singlPrintf("float R0 = %f;\n", R0);
    }
    if (do_point_mass) {
	singlPrintf("Point mass %g at (%g,%g,%g)\n", point_mass.mass,
	       point_mass.pos[0], point_mass.pos[1], point_mass.pos[2]);
	singlPrintf("Point mass vel (%g,%g,%g)\n", 
	       point_mass.vel[0], point_mass.vel[1], point_mass.vel[2]);
    }

    singlFflush();
    SanityCheck(btab, nobj, gnobj, &mtot);

    pqsortsetup(&sortedbtab, btab, nobj, sizeof(body), sort_tol, Realloc_f);
    
    if (log_time) {
	float dpdt;
	dpdt = 2./3. * pow((double)tpos, -1./3.);
	ConvertVPofT(btab, nobj, dpdt, 0);
    }

    SetupTree(&thetree, NDIM, 
	      sizeof(body), sizeof(cell), TBODYSZ, sizeof(cofmdata),
	      (pq_keyproto)GetKeyFromStruct, (pq_wgtproto)GetCost,
	      CofmFromDaugh, (cellfromcofm_t)CellFromCofm);

    if (frac_tol != (float)0.0) { /* init acc_last for first step */
	/* Reset timers and counters */
	ClearEnabledTimers();
	ClearEnabledCounters();
	StartTimer(&StepTotWC);
	StartTimer(&StepTot);
	if (do_periodic) {
	    if (cosmology)
	      sysradius = R0 * (1.0+1e-5) / (1.0 + Znow(tpos));
	    else
	      sysradius = R0;
	    VS(rmin, = -sysradius);
	    VS(rmax, = sysradius);
	    WrapPeriodic(btab, nobj, rmin, rmax, 2.0*sysradius, cosmology,
			 log_time, tpos);
	    FixRsizeExact(rmin, rmax);
	} else {
	    FindBbox(btab, nobj, rmin, rmax);
	    sysradius = 0.5*FixRsize(rmin, rmax);
	}
	this_tol = tol*mtot/(sysradius*sysradius);
	if (comov_eps) {
	     /* comoving smoothing */
	    this_eps = eps/(Znow(tpos)+(float)1.0);
	} else {
	    this_eps = eps;
	}
	SetTol(this_tol, frac_tol, cosmo.GNewt, this_eps, gnobj);
	FixKeys(btab, nobj, GETKEY);
	singlPrintf("BuildTree, most approx. mac\n");
	if( Msg_test("printtree") ){
	  int i;
	  Msg_do("Btab before BuildTree\n");
	  for(i=0; i<nobj; i++){
	    Msg_do("%s\n", PrintBodyContents(&btab[i]));
	  }
	}
	BuildTree(&thetree, &sortedbtab);
	if( Msg_test("printtree") ){
	  int i;
	  Msg_do("Btab after BuildTree\n");
	  for(i=0; i<nobj; i++){
	    Msg_do("%s\n", PrintBodyContents(&btab[i]));
	  }
	  Msg_do("Tree\n");
	  PrintTree(&thetree, PrintBodyContents, PrintCellContents);
	}
	btab = sortedbtab.data;
	nobj = sortedbtab.nobj;
	singlPrintf("BuildTree done\n");

	MPMY_Sync(); /* No sync might cause msg buffer overflow? */
	StartTimer(&FindForcesTm);
	WalkInit(&thetree, &thetree, sizeof(Sink), (macv_t) Lowestmacv,
		 (inherit_t)InheritSink);

	for (p = btab; p < btab+nobj; p++) {
	    VS(p->acc, = (float)0.0);
	    p->phi = (float)0.0;
	    p->nterms = 0;
	}
	if (do_periodic) {
	    singlPrintf("FindForces (periodic), this_eps=%g\n", this_eps);
	    StartTimer(&PeriodicForceTm);
	    Periodic(&thetree, 2.0*sysradius);
	    StopTimer(&PeriodicForceTm);
	    singlPrintf("FindForces (periodic) done\n");
	} 
	singlPrintf("FindForces, this_eps=%g\n", this_eps);
	WalkNT(&thetree);
	WalkTerminate();
	if (cosmology) FixGlobalForce(btab, nobj);
	StopTimer(&FindForcesTm);
	singlPrintf("FindForces done\n");

	FreeTree(&thetree);
	singlPrintf("FreeTree done\n");
	Msgf(("FreeTree done\n"));
	if (do_periodic) FixCube(btab, nobj, sysradius, cosmo.GNewt*mtot);
	StopTimer(&StepTot);
	StopTimer(&StepTotWC);
	OutputTimer(&StepTot, singlPrintf);
	OutputTimer(&StepTotWC, singlPrintf);

	for (p = btab; p < btab+nobj; p++) {
	    float sacc2;
	    sacc2 = Dot(p->acc, p->acc);
	    p->acc_last = sqrt(sacc2);
	}
    }
    if (do_NlgN) {
	inherit = (inherit_t)InheritSinkNlogN;
	mac = (macv_t)Nlogngate;
    } else {
	inherit = (inherit_t)InheritSink;
	mac = (frac_tol == 0.0) ? (macv_t)Unifiedmacv : (macv_t)Fracmacv;
    }

    for (nsteps += iter; iter <= nsteps; iter++) {
	if (timeout > 0) MPMY_TimeoutReset(timeout);
	/* If hetero_load_balance flag set, store load balance target value
	 * else default constant value from pqsortsetup() used */
	if (hetero_load_balance)
	    SetLBTarget(&sortedbtab, hetero_load_balance);
	/* Reset timers and counters */
	ClearEnabledTimers();
	ClearEnabledCounters();
	StartTimer(&StepTotWC);
	StartTimer(&StepTot);
	if (do_periodic) {
	    if (cosmology)
	      sysradius = R0 * (1.0+1e-5) / (1.0 + Znow(tpos));
	    else
	      sysradius = R0;
	    VS(rmin, = -sysradius);
	    VS(rmax, = sysradius);
	    FixRsizeExact(rmin, rmax);
	} else {
	    FindBbox(btab, nobj, rmin, rmax);
	    sysradius = 0.5*FixRsize(rmin, rmax);
	}
	Msgf(("rmin=(%g, %g, %g) rmax=(%g, %g, %g)\n",
	      rmin[0], rmin[1], rmin[2], 
	      rmax[0], rmax[1], rmax[2]));
	this_tol = tol*mtot/(sysradius*sysradius);
	if (comov_eps) {
	    /* comoving smoothing */
	    this_eps = eps/(Znow(tpos)+(float)1.0);
	} else {
	    this_eps = eps;
	}
	
        if (do_nsquared) {
            set_eps(this_eps);
	    for (p = btab; p < btab+nobj; p++) {
		VS(p->acc, = (float)0.0);
		p->phi = (float)0.0;
		p->nterms = 0;
	    }
	    StartTimer(&FindForcesTm);
            Ring(btab, sizeof(body), nobj, btab, sizeof(body), nobj, 
                 TBODYSZ, set_body, do_grav2);
            for (p = btab; p < btab+nobj; p++) {
                VS(p->acc, *= cosmo.GNewt);
                p->phi *= cosmo.GNewt;
            }
	    StopTimer(&FindForcesTm);
	} else {
	    SetTol(this_tol, frac_tol, cosmo.GNewt, this_eps, gnobj);
	    FixKeys(btab, nobj, GETKEY);
	    if (frac_tol == 0.0)
	      singlPrintf("BuildTree, this_tol=%g\n", this_tol);
	    else
	      singlPrintf("BuildTree, frac_tol=%g\n", frac_tol);
	    StartTimer(&BuildTot);
#if 0
	    if (decomp_iter % decomp_freq == 0 || decomp_iter < 2) {
		singlPrintf("Setting decomptab\n");
		if (decomptab) Free(decomptab);
		decomptab = SaveDecomp();
	    }
	    else SetDecomp(decomptab);
#endif
	    decomp_iter++;
	    if( Msg_test("printtree") ){
	      int i;
	      Msg_do("Btab before BuildTree\n");
	      for(i=0; i<nobj; i++){
		Msg_do("%s\n", PrintBodyContents(&btab[i]));
	      }
	    }
	    BuildTree(&thetree, &sortedbtab);
	    if( Msg_test("printtree") ){
	      int i;
	      Msg_do("Btab after BuildTree\n");
	      for(i=0; i<nobj; i++){
		Msg_do("%s\n", PrintBodyContents(&btab[i]));
	      }
	      Msg_do("Tree\n");
	      PrintTree(&thetree, PrintBodyContents, PrintCellContents);
	    }
	    btab = sortedbtab.data;
	    nobj = sortedbtab.nobj;
	    StopTimer(&BuildTot);
	    singlPrintf("BuildTree done %d (%d)\n", maxmem(), maxheap());
	    AddCounter(&NbodyCnt, nobj);

	    /* Periodic does multiple calls to Walk, so we must init here */
	    /* rather than in inherit */
	    for (p = btab; p < btab+nobj; p++) {
		VS(p->acc, = (float)0.0);
		p->phi = (float)0.0;
		p->nterms = 0;
	    }

	    MPMY_Sync(); /* No sync might cause msg buffer overflow? */
	    StartTimer(&FindForcesTm);
	    StartTimer(&WITm);
	    WalkInit(&thetree, &thetree, sizeof(Sink), mac, inherit);
	    StopTimer(&WITm);
	    StartTimer(&PerTm);
	    if (do_periodic) {
		singlPrintf("FindForces (periodic), this_eps=%g\n", this_eps);
		StartTimer(&PeriodicForceTm);
		Periodic(&thetree, 2.0*sysradius);
		StopTimer(&PeriodicForceTm);
		singlPrintf("FindForces (periodic) done\n");
		for (p = btab; p < btab+nobj; p++) { /* only fundamental phi */
		    p->phi = (float)0.0;
		}
	    }
	    StopTimer(&PerTm);
	    singlPrintf("FindForces, this_eps=%g\n", this_eps);
	    StartTimer(&WNTTm);
	    WalkNT(&thetree);
	    StopTimer(&WNTTm);
	    StartTimer(&WTermTm);
	    if (cosmology) FixGlobalForce(btab, nobj);
	    WalkTerminate();
	    StopTimer(&WTermTm);
	    StopTimer(&FindForcesTm);
	    singlPrintf("FindForces done %d (%d)\n", maxmem(), maxheap());

	    MPMY_Sync();
	    /* This should be the high-water mark for memory use */
	    AddCounter(&MemCnt, malloc_used()/1024);

	    FreeTree(&thetree);
	    singlPrintf("FreeTree done %d (%d)\n", maxmem(), maxheap());
	    Msgf(("FreeTree done\n"));
	}
	if (do_periodic) FixCube(btab, nobj, sysradius, cosmo.GNewt*mtot);

	if (setpvel) {
	    setpvel = 0;
	    set_vels(btab, nobj, tpos);
	    singlPrintf("Velocities adjusted to linear approximation.\n");
	    if (log_time) {
		float dpdt;
		dpdt = 2./3. * pow((double)tpos, -1./3.);
		ConvertVPofT(btab, nobj, dpdt, 0);
	    }
	}

	if (image_freq && iter%image_freq == 0) {
	    char name[256];
	    float sysr, image_rmin[3], image_rmax[3];

	    if (cosmology)
	      sysr = R0 * (1.0+1e-5) / (1.0 + Znow(tpos));
	    else
	      sysr = R0;
	    VS(image_rmin, = -sysr);
	    VS(image_rmax, = sysr);
	    FixRsizeExact(image_rmin, image_rmax);

	    sprintf(name, "%s_img.%04d", outnamebase, iter);
	    Image(btab[0].pos, btab[0].pos+1, &(btab[0].mass),
		  sizeof(body), nobj, image_rmin, image_rmax, 
		  x_pixels, y_pixels, 10, 250, log_image, name);
	}

	if (do_point_mass) {
	    p = &point_mass;
	    update_point_mass(btab, nobj, p, this_eps, cosmo.GNewt);
	    MPMY_ICombine_Init(&req);
	    MPMY_ICombine(&(p->phi), &(p->phi), 1, MPMY_FLOAT, MPMY_SUM, req);
	    MPMY_ICombine(&(p->acc), &(p->acc), NDIM, MPMY_FLOAT,MPMY_SUM,req);
	    MPMY_ICombine_Wait(req);
	}


	if (ForceOutput() || CaughtTerm 
	    || (do_output && !first_step
		&& ((iter+output_freq) % output_freq == 0))
	    || (save_first && first_step)) {
	    sortresult_t outputsort;
	    outbodyptr output_btab;
	    float output_R0, output_z, output_h;
	    char outname[256];
	    int output_nobj = nobj;
	    float tpos_out = tpos;
	    float tvel_out = tvel; /* changed in Integrate() */

	    Msgf(("Doing output\n"));
	    output_btab = Malloc(output_nobj * sizeof(outbody));
	    for(i=0; i<output_nobj; i++){
		output_btab[i].mass = btab[i].mass;
		VV(output_btab[i].pos, = btab[i].pos);
		VV(output_btab[i].vel, = btab[i].vel);
#ifdef SAVE_ACC
		VV(output_btab[i].acc, = btab[i].acc);
		output_btab[i].phi = btab[i].phi;
#endif
		output_btab[i].ident = btab[i].ident;
	    }
	    /* Don't sort before Integrate_out or btab and output_btab */
	    /* will not be in the same order */
	    if (log_time) {
		IntegratePofT_out(btab, output_btab, output_nobj, dt,
				  &tpos_out, &tvel_out, &ke, &pe);
	    } else {
		Integrate_out(btab, output_btab, output_nobj, dt,
			      &tpos_out, &tvel_out, &ke, &pe);
	    }
	    pqsortsetup_order(&outputsort, output_btab, output_nobj,
			sizeof(outbody), 0.1, 1, Realloc_f);
	    output_btab = pqsort(&outputsort,
				 (pq_wgtproto)UnityCost, 
				 (pq_keyproto)OutIdentKey);
	    output_nobj = outputsort.nobj;
	    if (cosmology) {
		output_z = Znow(tpos_out);
		output_h = Hnow(tpos_out);
		output_R0 = R0;
	    } else {
		output_z = 0.0;
		output_h = 0.0;
		output_R0 = sysradius;
	    }
	    Msgf(("After output pqsort, %d outbodies\n", output_nobj));
	    MPMY_ICombine_Init(&req);
	    MPMY_ICombine(&ke, &ke, 1, MPMY_DOUBLE, MPMY_SUM, req);
	    MPMY_ICombine(&pe, &pe, 1, MPMY_DOUBLE, MPMY_SUM, req);
	    MPMY_ICombine_Wait(req);
	    if (do_point_mass) {
		p = &point_mass;
		ke += (float)0.5 * p->mass * Dot(p->vel, p->vel);
		pe += (float)0.5 * p->mass * p->phi;
	    }
	    sprintf(outname, "%s.%04d", outnamebase, iter);
	    if (write_nfiles) MPMY_Nfileio(1);
	    SDFwrite(outname, gnobj, 
		      output_nobj, output_btab, sizeof(outbody), OUTBODYDESC,
		      "npart", SDF_INT, gnobj,
		      "eps", SDF_FLOAT, eps,
		      "Gnewt", SDF_FLOAT, cosmo.GNewt,
		      "tolerance", SDF_FLOAT, this_tol*cosmo.GNewt,
		      "iter", SDF_INT, iter,
		      "tpos", SDF_FLOAT, tpos_out,
		      "tvel", SDF_FLOAT, tvel_out,
		      "R0", SDF_FLOAT, output_R0,
		      "Omega0", SDF_FLOAT, cosmo.Omega0,
		      "H0", SDF_FLOAT, cosmo.H0,
		      "Lambda_prime", SDF_FLOAT, cosmo.Lambda,
		      "hubble", SDF_FLOAT, output_h,
		      "redshift", SDF_FLOAT, output_z,
		      "ke", SDF_DOUBLE, ke,
		      "pe", SDF_DOUBLE, pe,
		      NULL);
	    if (write_nfiles) MPMY_Nfileio(0);
	    Free(output_btab);
	    singlPrintf("\nOutput to %s done.\n", outname);
	    Msgf(("Output to %s done\n", outname));
#ifndef __DELTA__
	    if (MPMY_Procnum() == 0) {
		char name[256];
		sprintf(name, "%s.restart", outnamebase);
		if (unlink(name))
		  Shout("unlink of %s failed, errno=%d\n", name, errno);
		if (symlink(outname, name))
		  Shout("symlink of %s failed, errno=%d\n", outname, errno);
	    }
#endif
	}

	if (ForceStop() || CaughtTerm) {
	    singlPrintf("Stopping.\n");
	    break;
	}

	Msgf(("integrating positions\n"));
	tposlast = tpos;
	if (do_point_mass) {
	    float tpos_tmp = tpos;
	    float tvel_tmp = tvel;
	    Integrate(&point_mass, 1, dt, &tpos_tmp, &tvel_tmp, &ke, &pe);
	}
	if (log_time) {
	    IntegratePofT(btab, nobj, dt, &tpos, &tvel, &ke, &pe);
	} else {
	    Integrate(btab, nobj, dt, &tpos, &tvel, &ke, &pe);
	}
	if(cosmology){
	    CosmoPush(&cosmo, tpos);
	    Msgf(("Pushed cosmo params to tpos=%g, Z=%g\n",
		  tpos, Znow(tpos)));
	}

	if (do_periodic) {
	    if (cosmology)
	      sysradius = R0 * 1.0 / (1.0 + Znow(tpos));
	    else
	      sysradius = R0;
	    VS(rmin, = -sysradius);
	    VS(rmax, = sysradius);
	    WrapPeriodic(btab, nobj, rmin, rmax, 2.0*sysradius, cosmology,
			 log_time, tpos);
	}

	VS(force, = 0.0);
	VS(com, = 0.0);
	VS(comv, = 0.0);
	acc2 = 0.0;
	mtot = 0.0;
	ntermslocal = 0.0;
	for (p = btab; p < btab+nobj; p++) {
	    float sacc2;
	    VV(com, += p->mass*p->pos);
	    VV(comv, += p->mass*p->vel);
	    VV(force, += p->mass*p->acc);
	    sacc2 = Dot(p->acc, p->acc);
	    p->acc_last = sqrt(sacc2);
	    acc2 += sacc2;
	    mtot += p->mass;
	    ntermslocal += p->nterms;
	    if (p->nterms <= 0) SeriousWarning("nterms is %f\n", p->nterms);
	}
	AddCounter(&KNtermsCnt, (int)(ntermslocal/1000.)); /* might overflow?? */

	Msgf(("doing MPMY_combine\n"));
	MPMY_ICombine_Init(&req);
	MPMY_ICombine(force, force, NDIM, MPMY_DOUBLE, MPMY_SUM, req);
	MPMY_ICombine(com, com, NDIM, MPMY_DOUBLE, MPMY_SUM, req);
	MPMY_ICombine(comv, comv, NDIM, MPMY_DOUBLE, MPMY_SUM, req);
	MPMY_ICombine(&acc2, &acc2, 1, MPMY_DOUBLE, MPMY_SUM, req);
	MPMY_ICombine(&ke, &ke, 1, MPMY_DOUBLE, MPMY_SUM, req);
	MPMY_ICombine(&pe, &pe, 1, MPMY_DOUBLE, MPMY_SUM, req);
	MPMY_ICombine(&mtot, &mtot, 1, MPMY_DOUBLE, MPMY_SUM, req);
	MPMY_ICombine_Wait(req);
	Msgf(("done MPMY_combine\n"));
	if (do_point_mass) {
	    p = &point_mass;
	    ke += (float)0.5 * p->mass * Dot(p->vel, p->vel);
	    pe += (float)0.5 * p->mass * p->phi;
	}

	StopTimer(&StepTot);
	StopTimer(&StepTotWC);

	if (cosmology) 
	    singlPrintf("\ntpos = %g, znow = %g, iter = %d, size = %g, eps = %g\n", 
			tposlast, Znow(tposlast),
			iter, sysradius, this_eps);
	else
	    singlPrintf("\ntpos = %g, iter = %d, size = %g\n",
			tposlast, iter, sysradius);
	singlPrintf("ke = %g, pe = %g, energy = %g\n", ke, pe, ke+pe);
	VS(force, /= mtot);
	VS(com, /= mtot);
	VS(comv, /= mtot);
	singlPrintf("CM accel: (" Sinfix("%g", " ") "): %g\n",
		    Vinfix(force, COMMA), sqrt(Dot(force, force)));
	singlPrintf("rms accel: %g\n", sqrt(acc2/gnobj));
	if (do_point_mass) {
	    singlPrintf("Point mass %g at (%g,%g,%g) (%g,%g,%g)\n",
			point_mass.mass, point_mass.pos[0], point_mass.pos[1],
			point_mass.pos[2], point_mass.vel[0],
			point_mass.vel[1], point_mass.vel[2]);
	    singlPrintf("CofM at (%g,%g,%g) (%g,%g,%g)\n",
			com[0], com[1], com[2], comv[0], comv[1], comv[2]);
	}
	AddCounter(&HeapCnt_, malloc_heapsz()/1024);
		    
	if( timer_freq && iter%timer_freq == 0 ){
	    OutputTimers(singlPrintf);
	    OutputCounters(singlPrintf);
	    if( Msg_test("timers") ){
		/* This can be very tedious on a big machine. */
		OutputIndividualTimers(Msg_do);
		OutputIndividualCounters(Msg_do);
	    }
	} else {
	    OutputTimer(&StepTot, singlPrintf);
	    OutputTimer(&StepTotWC, singlPrintf);
	}
#ifdef _AIX
	PrintLoad(&StepTot, &StepTotWC);
#endif	
	singlFflush();

	first_step = 0;
	if( Msg_test("memleak") ){
	    Msg_do("Memory map after iteration %d\n", iter);
	    malloc_print();
	}
    }
    singlPrintf("Bye!\n");
    Msgf(("Bye!\n"));
    Msg_flush();
    MPMY_Finalize();
    exit(0);
}

static SDF *startup(int argc, char **argv){
    SDF *csdfp;
    char msg_turn_on[512];
    char msgdir[256];
    char tmp[256];
    char *msgbase, *lastslash;
    char cfile[256];
    int Msg_memfile;

    if (argc > 1)
 	strncpy(cfile, argv[1], sizeof(cfile));
    else
 	Getsparam("control file", cfile);
    if ((csdfp = SDFopen(NULL, cfile)) == NULL) {
 	SinglError("Sorry, couldn't SDFopen %s\n%s\n",
	      cfile, SDFerrstring);
    }
    singlPrintf("cfile \"%s\" opened\n", cfile);
    SDFgetintOrDefault(csdfp, "Msg_memfile", &Msg_memfile, 0);
    if (Msg_memfile) {
#if defined(__PARAGON__) || defined(_AIX) || defined(sparc)
	sigio_setup();
#endif
	memfile_init(Msg_memfile);
	Msg_addfile(0, (Msgvfprintf_t)memfile_vfprintf, 0);
	singlPrintf("Putting all Msgs in memfile\n");
    } else {
	/* Get the msgdir either from:
	   argv[2] 
	   "msgbase" in csdfp
	   misc.argv[0]/msg

	   We then append .<procnum> to the name
	   */
	if( argc > 2 ){
	    msgbase = argv[2];
	}else if( SDFgetstring(csdfp, "msgbase", tmp, sizeof(tmp))==0 ){
	    msgbase = tmp;
	}else{
	    lastslash = strrchr(argv[0], '/');
	    if( lastslash ){
		msgbase = lastslash+1;
	    }else{
		msgbase = argv[0];
	    }
	    sprintf(tmp, "misc.%s/msg", msgbase);
	    msgbase = tmp;
	}	
	sprintf(msgdir, "%s.%d", msgbase, MPMY_Procnum());
	MsgdirInit(msgdir);
    }
    SDFgetstringOrDefault(csdfp, "Msg_turn_on", msg_turn_on, 
			  sizeof(msg_turn_on), "");
    Msg_turnon(msg_turn_on);

    if( Msg_test("bigmalloc.c") ){
	malloc_debug(2);
	Msg_do("Malloc_debug(2), expect slow mallocs\n");
    }else{
	malloc_debug(1);
    }

    EnableTimer(&StepTot, "Step Total");
    EnableWCTimer(&StepTotWC, "Step Tot(WC)");
    EnableTimer(&BuildTot, "Build Total");
    EnableTimer(&DecompTm, "Decomp");
    EnableTimer(&DecompCommTm, "DecompComm");
    EnableTimer(&SortTm, "Sort");
    EnableTimer(&MakeTreeTm, "Make Tree");
    EnableTimer(&FindForcesTm, "Force Eval");
    EnableTimer(&GravTm, "Grav Time");
    EnableTimer(&WalkDeferTm, "Walk Defer");
    EnableTimer(&WITm, "WalkInit");
    EnableTimer(&WTermTm, "WalkTerm");
    EnableTimer(&WNTTm, "WalkNT");
    EnableTimer(&WITm, "WalkInit");
    EnableTimer(&PerTm, "Periodic");
    EnableCounter(&CCInt, "Cell-cell");
    EnableCounter(&BCInt, "Body-cell");
    EnableCounter(&CBInt, "Cell-body");
    EnableCounter(&BBInt, "Body-body");
    EnableCounter(&KNtermsCnt, "KNterms");
    EnableCounter(&NbodyCnt, "Nbody");
    EnableCounter(&CCIntRej, "MAC fail");
    EnableCounter(&SharedCnt, "Shared Cells");
    EnableCounter(&TranslateCnt, "Translate");
    EnableCounter(&DeferCnt, "Deferred");
    EnableCounter(&MemCnt, "Mem Used (K)");
    EnableCounter(&HeapCnt_, "Heap Sz (K)");
    return csdfp;
}

static void SanityCheck(bodyptr btab, int nobj, int gnobj, double *mtotp){
    double mtot;
    bodyptr p;
    int sumnobj;
    MPMY_Comm_request req;

    mtot = 0.0;
    for (p = btab; p < btab+nobj; p++) {
	mtot += p->mass;
	/* We really could use more checks here! */
    }
    sumnobj = nobj;
    MPMY_ICombine_Init(&req);
    MPMY_ICombine(&mtot, &mtot, 1, MPMY_DOUBLE, MPMY_SUM, req);
    MPMY_ICombine(&sumnobj, &sumnobj, 1, MPMY_INT, MPMY_SUM, req);
    MPMY_ICombine_Wait(req);
    assert(sumnobj == gnobj);
    Msgf(("Particle 0 (%d), %g, %g %g %g, %g %g %g\n",
	  btab->ident, btab->mass, 
	  btab->pos[0], btab->pos[1], btab->pos[2],
	  btab->vel[0], btab->vel[1], btab->vel[2]));
    Msgf(("Particle %d (%d), %g, %g %g %g, %g %g %g\n", nobj-1,
	  btab[nobj-1].ident, btab[nobj-1].mass, 
	  btab[nobj-1].pos[0], btab[nobj-1].pos[1], btab[nobj-1].pos[2],
	  btab[nobj-1].vel[0], btab[nobj-1].vel[1], btab[nobj-1].vel[2]));
    singlPrintf("gnobj = %d, mtot = %f\n", gnobj, mtot);
    *mtotp = mtot;
}

static void
SetLBTarget(sortresult_t *decompp, int hetero_load_balance) {

/* 
   These are suffering from code-rot. The counters they refer to have
   disappeared.  I'm not sure what they were counting anyway! */
#if 0
    switch (hetero_load_balance) {
    case 1:
	if ((ReadCounter(&GravCnt) > 0)
	    && (ReadCounter(&CCInt) > 0)
	    && (ReadCounter(&TranslateCnt) > 0)
	    && (ReadTimer(&GravTm) > 0.0))
	    decompp->loadbal_target = ReadTimer(&GravTm)
		/ (ReadCounter(&GravCnt)
		   + ReadCounter(&CCInt)
		   + ReadCounter(&TranslateCnt));
	break;


    case 2:
	if ((ReadCounter(&GravCnt) > 0)
	    && (ReadCounter(&CCInt) > 0)
	    && (ReadCounter(&TranslateCnt) > 0)
	    && (ReadTimer(&GravTm) > 0.0))
	    decompp->loadbal_target = ReadTimer(&GravTm)
		/ (ReadCounter(&GravCnt)
		   + 2.0*ReadCounter(&CCInt)
		   + 1.5*ReadCounter(&TranslateCnt));
	break;

    case 3:
	if ((ReadCounter(&KNtermsCnt) > 0)
	    && (ReadTimer(&GravTm) > 0.0))
	    decompp->loadbal_target = ReadTimer(&GravTm)
		/ ReadCounter(&KNtermsCnt);
	break;

    }
#endif

    return;
}

static void
FixGlobalForce(body *xptr, int n){
    /* Make whatever corrections are necessary to the acceleration, etc.
       based on values of GNewt, Lambda, etc., etc. */
    float G = cosmo.GNewt;
    float lambdafac;
    body *p;

    lambdafac = cosmo.Lambda*cosmo.H0*cosmo.H0;
    while(--n){
	p = xptr++;
#if 0				/* These are still in grav for msw version */
	VS(p->acc, *= G);
	p->phi *= G;
	p->errsum *= G;
	p->errsum2 *= G*G;
#endif
	VV(p->acc, += lambdafac*p->pos);
    }
}

/* Where the time variable p is a function of the real time t */
/* We need to know physics info, because the energy equation is complicated */

/* The initial conditions contain x(t) and v(t)
    we want to use x(p(t)) and v(p(t))

    Thus via the chain rule, dx/dt = dx/dp dp/dt = x' dp/dt
    and d2x/dt2 = d2x/dp2 (dp/dt)^2 + dx/dp d2p/dt2 
       = x'' (dp/dt)^2 + x' d2p/dt2

    To convert vel to x' we divide by dp/dt, and must convert back when
    finding the physical energy, or writing output 

    The acc field is d2x/dt2, so acc = x'' (dp/dt)^2 + x' d2p/dt2
    We use x'' to evolve x', and x'' = (acc - x' d2p/dt2)/(dp/dt)^2
    This assumes x' is at the same time as x, which is only true
    for the first step of the leapfrog integrator.  However, when we
    write the equation for advancing x'_i-1/2 to x'_i+1/2, we can
    replace x'_i with the average of x'_i-1/2 and x'_i+1/2, and solve
    for x'_i+1/2.  This results in

    x'_i+1/2 = x_i-1/2 * (1 - z)/(1 + z) + acc \delta p / ( (dp/dt)^2 * z)

    where z = 1 - \delta p * d2p/dt2 / (2 * (dp/dt)^2)
    
 */

static void
IntegratePofT(body *xptr, const int n, const float dp, float *tpos, 
	      float *tvel, double *kep, double *pep)
{
    body *end = xptr + n;
    const float dp_dt = 2./3.*pow(*tpos, -1./3.);
    const float d2p_dt2 = -2./9.*pow(*tpos, -4./3.);
    const float dp_dt2 = dp_dt * dp_dt;
    const float dp_on_dp_dt2 = (float)dp / dp_dt2;
    const float dphalf_on_dp_dt2 = (float)0.5 * dp / dp_dt2;
    const float a = dp * d2p_dt2 / ((float)2.0 * dp_dt2);
    const float a_hi = (float)1.0 - a;
    const float a_lo = (float)1.0 / ((float)1.0 + a);
    float vcentered[NDIM];
    double ke = 0.0;
    double pe = 0.0;
    
    if (*tvel < *tpos) {
	for (; xptr < end; xptr++) {
	    VVV(vcentered, = a_hi * xptr->vel, + dphalf_on_dp_dt2 * xptr->acc);
	    VS(vcentered, *= dp_dt); /* convert to physical vel */
	    ke += xptr->mass * Dot(vcentered, vcentered);
	    pe += xptr->mass * xptr->phi;
	    VVV(xptr->vel, = a_hi * xptr->vel, + dp_on_dp_dt2 * xptr->acc);
	    VS(xptr->vel, *= a_lo);
	    VV(xptr->pos, += dp * xptr->vel);
	}
	*tvel = pow(pow(*tvel, 2./3.) + dp, 3./2.);
	*tpos = pow(pow(*tpos, 2./3.) + dp, 3./2.);
    } else if (*tvel == *tpos) {
	for (; xptr < end; xptr++) {
	    VV(vcentered, = dp_dt * xptr->vel); /* convert to physical vel */
	    ke += xptr->mass * Dot(vcentered, vcentered);
 	    pe += xptr->mass * xptr->phi;
	    VVV(xptr->vel, = a_hi * xptr->vel, + dphalf_on_dp_dt2 * xptr->acc);
	    VV(xptr->pos, += dp * xptr->vel);
	}
	*tvel = pow(pow(*tvel, 2./3) + 0.5 * dp, 3./2.);
	*tpos = pow(pow(*tpos, 2./3) + dp, 3./2.);
    } else {
	Error("Bad state in IntegratePofT\n");
    }
    *kep = 0.5 * ke;
    *pep = 0.5 * pe;
}

static void
Integrate(body *xptr, const int n, const float dt, float *tpos, float *tvel,
	  double *kep, double *pep)
{
    body *end = xptr + n;
    float vcentered[NDIM];
    float dt_half = (float)0.5 * dt;
    double ke = 0.0;
    double pe = 0.0;

    if (*tvel < *tpos) {	/* leapfrog step */
	for (; xptr < end; xptr++) {
	    VVV(vcentered, = xptr->vel, + dt_half * xptr->acc);
	    ke += xptr->mass * Dot(vcentered, vcentered);
	    pe += xptr->mass * xptr->phi;
	    VV(xptr->vel, += dt * xptr->acc);
	    VV(xptr->pos, += dt * xptr->vel);
	}
	*tvel += dt;
	*tpos += dt;
    } else if (*tvel == *tpos) { /* first step */
	for (; xptr < end; xptr++) {
	    ke += xptr->mass * Dot(xptr->vel, xptr->vel);
 	    pe += xptr->mass * xptr->phi;
	    VV(xptr->vel, += dt_half * xptr->acc);
	    VV(xptr->pos, += dt * xptr->vel);
	}
	*tvel += dt_half;
	*tpos += dt;
    } else {
	Error("Bad state in Integrate\n");
    }
    *kep = 0.5 * ke;
    *pep = 0.5 * pe;
}


static void
IntegratePofT_out(const body *xptr, outbody *yptr, const int n, 
		  const float dp, float *tpos, float *tvel, 
		  double *kep, double *pep)
{
    const body *end = xptr + n;
    const float dp_dt = 2./3.*pow(*tpos, -1./3.);
    const float d2p_dt2 = -2./9.*pow(*tpos, -4./3.);
    const float dp_dt2 = dp_dt * dp_dt;
    const float dphalf_on_dp_dt2 = (float)0.5 * dp / dp_dt2;
    const float a = dp * d2p_dt2 / ((float)2.0 * dp_dt2);
    const float a_hi = (float)1.0 - a;
    double ke = 0.0;
    double pe = 0.0;
    
    if (*tvel == *tpos) {
	/* It must be the first step, so don't update anything */
	for (; xptr < end; xptr++, yptr++) {
	    VS(yptr->vel, *= dp_dt); /* convert to physical vel */
	    ke += yptr->mass * Dot(yptr->vel, yptr->vel);
	    pe += yptr->mass * xptr->phi;
	}
    } else {
	for (; xptr < end; xptr++, yptr++) {
	    VVV(yptr->vel, = a_hi * yptr->vel, 
		+ dphalf_on_dp_dt2 * xptr->acc);
	    VS(yptr->vel, *= dp_dt); /* convert to physical vel */
	    ke += yptr->mass * Dot(yptr->vel, yptr->vel);
	    pe += yptr->mass * xptr->phi;
	}
	*tvel = pow(pow(*tvel, 2./3) + 0.5 * dp, 3./2.);
    }
    *kep = 0.5 * ke;
    *pep = 0.5 * pe;
}

static void
Integrate_out(const body *xptr, outbody *yptr, const int n, const float dt, 
	      float *tpos, float *tvel, double *kep, double *pep)
{
    const body *end = xptr + n;
    float dt_half = (float)0.5 * dt;
    double ke = 0.0;
    double pe = 0.0;

    if (*tvel == *tpos) {
	/* It must be the first step, so don't update anything */
	for (; xptr < end; xptr++, yptr++) {
	    ke += yptr->mass * Dot(yptr->vel, yptr->vel);
	    pe += yptr->mass * xptr->phi;
	}
    } else {
	for (; xptr < end; xptr++, yptr++) {
	    VV(yptr->vel, += dt_half * xptr->acc);
	    ke += yptr->mass * Dot(yptr->vel, yptr->vel);
	    pe += yptr->mass * xptr->phi;
	} 
	*tvel += dt_half;
    }
    *kep = 0.5 * ke;
    *pep = 0.5 * pe;
}


static void
ConvertVPofT(body *xptr, int n, float dp_dt, int to_physical)
{
    body *end = xptr + n;
    float one_on_dp_dt = (float)1.0 / dp_dt;

    if (to_physical) {
	for (; xptr < end; xptr++) {
	    /* convert to physical vel */
	    VV(xptr->vel, = dp_dt * xptr->vel);
	}
    } else {
	for (; xptr < end; xptr++) {
	    /* convert from physical vel */
	    VV(xptr->vel, = one_on_dp_dt * xptr->vel);
	}
    }
}

#define one_kpc (3.08567802e16) /* km */
#define one_Gyr (3.1558149984e16) /* sec */

float 
Anow(float time)
{
    struct cosmo_s foo;

    foo = cosmo;
    CosmoPush(&foo, time);
    return foo.a;
}    

float
Znow(float time)
{
    return 1.F/Anow(time) - 1.F;
}

float
Hnow(float time)
{
    struct cosmo_s foo;
    float Omega0 = cosmo.Omega0;
    float Lambda = cosmo.Lambda;
    float H0 = cosmo.H0;
    float a2, a3;

    foo = cosmo;
    CosmoPush(&foo, time);
    
    a2 = foo.a*foo.a;
    a3 = foo.a*a2;
    return H0 * sqrt( Omega0/a3 + Lambda - (Omega0 + Lambda - 1.F)/a2);
}

/* This erases the velocities, and sets them */
/* according to linear theory */

static void
set_vels(body *p, int n, float real_time)
{
    float tmp[NDIM];
    body *end = p + n;
    float H;
    float acc_back;
    float vel_fac, pos_fac;
    float a;
    float asum1, asum2;
    
    a = Anow(real_time);
    H = Hnow(real_time);
    /* Note that the Lambda force has already been added in.  We must
       subtract it, along with the background Omega0 term to get the
       'peculiar' acceleration.  Yuck! */
    acc_back = (0.5*cosmo.Omega0/(a*a*a) - cosmo.Lambda)*cosmo.H0*cosmo.H0;

    vel_fac = (a*a*a* cosmo.Zel_f * H)/(1.5F * cosmo.Omega0 * cosmo.H0 * cosmo.H0);
    pos_fac = H;
    /* Velocities really store dx/dp = dx/dt * (dp/dt)^-1, but that
       correction is done elsewhere! */

    singlPrintf("set_vels, Zel_f = %g, vel_fac is %.2f * t, H is %f\n", 
		cosmo.Zel_f, vel_fac/real_time, H);

    asum1 = 0.F;
    asum2 = 0.F;
    for (; p < end; p++) {
	VVV(tmp, = p->acc, + acc_back*p->pos);    /* peculiar acc */
	asum1 += Dot(p->acc, p->pos); /* diagnostic */
	asum2 += Dot(tmp, p->pos); /* diagnostic */
	VVV(p->vel, = vel_fac*tmp, + H*p->pos); /* pec. vel + Hubble flow */
    }
    singlPrintf("Mean(proper acc dot position) = %g\n", asum1/n);
    singlPrintf("Mean(peculiar acc dot position) = %g\n", asum2/n);
}

#if 0
static void
set_vels(body *p, int n, float Omega0, float H0, float real_time)
{
    float tmp[NDIM];
    body *end = p + n;
    float H, z;
    float acc_back;
    float vel_fac;
    
    z = Znow(Omega0, H0, real_time);
    H = Hubble_param(Omega0, H0, z);
    acc_back = 0.5*Omega0*H0*H0*(1.0+z)*(1.0+z)*(1.0+z);
    /* See Peebles LSS 14.8 */
    /* We approximate the perturbation parameter (f) as pow(Omega0, 0.6) */
    /* vel_fac is equal to physical time for Omega=1 */
    vel_fac = 2.0/(3.0*pow(Omega0, 0.4)*H); 

    singlPrintf("set_vels, vel_fac is %.2f * t, H is %f\n", 
		vel_fac/real_time, H);

    for (; p < end; p++) {
	VVS(tmp, = p->pos, * acc_back);    /* background acc */
	VVV(p->vel, = p->acc, + tmp);   /* peculiar acc */
	VS(p->vel, *= vel_fac);
	VVS(tmp, = p->pos, * H); /* hubble flow */
	VV(p->vel, += tmp);
    }
}

static void
set_vels_old(body *p, int n, float real_time)
{
    float tmp[NDIM];
    double one_on_t = 1.0/real_time;
    double factor = (2.0/3.0)*one_on_t;
    double factor2 = (2.0/9.0)*one_on_t*one_on_t; 
    body *end = p + n;

    for (; p < end; p++) {
	VVS(tmp, = p->pos, * factor2);    /* background acc */
	VVV(p->vel, = p->acc, + tmp);   /* peculiar acc */
	VS(p->vel, *= real_time);
	VVS(tmp, = p->pos, * factor); /* hubble flow */
	VV(p->vel, += tmp);
    }
}
#endif

void
update_point_mass(body *btab, int nobj, body *p, float eps, float newt)
{
    body *r;
    float dr2, oneor, oneor2;
    float phii;
    float smooth2 = eps*eps;
    Vxd(float r);
    Vxd(float ppos);

    p->phi = (float)0.0;
    VS(p->acc, = (float)0.0);
    VxV(ppos, = p->pos);

    for (r = btab; r < btab+nobj; r++) {
	VxVVx(r, = r->pos, - ppos); /* 3 flops */
	
	dr2 = Dotx(r, r);	/* 5 flops */
	
	dr2 += smooth2;
	
	oneor = recipsqrtf(dr2);	/* 8 flops */
	
	oneor2 = oneor * oneor;	/* 17 flops */
	phii = oneor * p->mass * newt;
	r->phi -= phii;
	VVx(r->acc, -= oneor2 * phii * r);
	phii = oneor * r->mass;
	p->phi -= phii;
	VVx(p->acc, += oneor2 * phii * r);
    }
}

static void
read_point_mass(body *point_mass, SDF *sdfp)
{
    SDFgetfloatOrDie(sdfp, "pt_x", &(point_mass->pos[0]));
    SDFgetfloatOrDie(sdfp, "pt_y", &(point_mass->pos[1]));
    SDFgetfloatOrDie(sdfp, "pt_z", &(point_mass->pos[2]));
    SDFgetfloatOrDie(sdfp, "pt_vx", &(point_mass->vel[0]));
    SDFgetfloatOrDie(sdfp, "pt_vy", &(point_mass->vel[1]));
    SDFgetfloatOrDie(sdfp, "pt_vz", &(point_mass->vel[2]));
    SDFgetfloatOrDie(sdfp, "pt_mass", &(point_mass->mass));
}


static void
tidal_init(body **btab, int *nobj, int *gnobj, SDF *csdfp)
{
    body *gal1;
    body *p1;
    int n = *nobj;
    float off[NDIM], offv[NDIM];

    SDFgetfloatOrDie(csdfp, "offset_x", &off[0]);
    SDFgetfloatOrDie(csdfp, "offset_y", &off[1]);
    SDFgetfloatOrDie(csdfp, "offset_z", &off[2]);
    SDFgetfloatOrDie(csdfp, "offset_vx", &offv[0]);
    SDFgetfloatOrDie(csdfp, "offset_vy", &offv[1]);
    SDFgetfloatOrDie(csdfp, "offset_vz", &offv[2]);

    gal1 = *btab;

    for (p1 = gal1; p1 < gal1+n; p1++) {
	VV(p1->pos, += off);
	VV(p1->vel, += offv);
    }
    singlPrintf("Did gal offset +/- (%g,%g,%g), (%g,%g,%g)\n",
		off[0], off[1], off[2], offv[0], offv[1], offv[2]);
}

static void
collision_init(body **btab, int *nobj, int *gnobj, SDF *csdfp)
{
    body *star1;
    body *star2;
    body *p1, *p2;
    int n = *nobj;
    int gn = *gnobj;
    float off[NDIM], offv[NDIM];

    SDFgetfloatOrDie(csdfp, "offset_x", &off[0]);
    SDFgetfloatOrDie(csdfp, "offset_y", &off[1]);
    SDFgetfloatOrDie(csdfp, "offset_z", &off[2]);
    SDFgetfloatOrDie(csdfp, "offset_vx", &offv[0]);
    SDFgetfloatOrDie(csdfp, "offset_vy", &offv[1]);
    SDFgetfloatOrDie(csdfp, "offset_vz", &offv[2]);

    /* Double space for btab */
    *btab = Realloc(*btab, 2 * n * sizeof(body));
    *nobj *= 2;
    *gnobj *= 2;
    star1 = *btab;
    star2 = *btab + n;

    for (p1 = star1, p2 = star2; p1 < star2; p1++, p2++) {
	memcpy(p2, p1, sizeof(body));
	p2->ident += gn;
	VV(p2->pos, -= off);
	VV(p2->vel, -= offv);

	VV(p1->pos, += off);
	VV(p1->vel, += offv);
    }
    singlPrintf("Did collision pair offset +/- (%g,%g,%g), (%g,%g,%g)\n",
		off[0], off[1], off[2], offv[0], offv[1], offv[2]);
}

static void
Periodic(tree_t *tp, float width)
{
    int i, j, k;
    float offset[NDIM];
    
    for (i = -1; i <= 1; i++) {
        offset[0] = i * width;
        for (j = -1; j <= 1; j++) {
            offset[1] = j * width;
            for (k = -1; k <= 1; k++) {
                offset[2] = k * width;
		SetGravOffset(offset);
                if (i || j || k) WalkNT(tp);
            }
        }
    }
    UnSetGravOffset();
}

static void
WrapPeriodic(body *bp, int n, float *rmin, float *rmax, float sz, 
	     int cosmology, int log_time, float tvel)
{
    body *b;
    int flux[NDIM] = {0, 0, 0};
    float vsz; /* hubble flow */
    
    for(b=bp; b<&bp[n]; b++) {
	VVVS(if LPAREN b->pos, > rmax, RPAREN flux, += 1);
	VVVS(if LPAREN b->pos, < rmin, RPAREN flux, -= 1);
    }
    MPMY_Combine(flux, flux, NDIM, MPMY_INT, MPMY_SUM);
    singlPrintf("Flux %d %d %d\n", flux[0], flux[1], flux[2]);
    if (!cosmology) {
	for(b=bp; b<&bp[n]; b++) {
	    VVVS(if LPAREN b->pos, > rmax, RPAREN b->pos, -= sz);
	    VVVS(if LPAREN b->pos, < rmin, RPAREN b->pos, += sz);
	}
    } else {
	if (log_time) 
	  vsz = sz*Hnow(tvel)*1.5*pow((double)tvel, 1./3.); /* ?? */
	else 
	  vsz = sz*Hnow(tvel);
	for(b=bp; b<&bp[n]; b++) {
	    VVVS(if LPAREN b->pos, > rmax, RPAREN b->vel, -= vsz);
	    VVVS(if LPAREN b->pos, > rmax, RPAREN b->pos, -= sz);
	    VVVS(if LPAREN b->pos, < rmin, RPAREN b->vel, += vsz);
	    VVVS(if LPAREN b->pos, < rmin, RPAREN b->pos, += sz);
	}
    }
}    


float gammln(float xx);
void gcf(float *gammcf, float a, float x, float *gln);
void gser(float *gamser, float a, float x, float *gln);
float gammq(float a, float x);
float gammp(float a, float x);
float erffc(float x);
void ewald(float *x, float L, float *phi, float *f);

#if 0
static void 
FixCube(body *b, int nobj, float l, float gm)
{
    int i;
    float f[NDIM];
    float r;
    float phi;

    StartTimer(&FixCubeTm);

    /* We did the nearest boxes already */
    l *= 3.0;
    gm *= NDIM*NDIM*NDIM;
    for (i = 0; i < nobj; i++) {
	ewald(b[i].pos, l, &phi, f);
	r = sqrt(Dot(b[i].pos,b[i].pos));
	/* Take out central point which we already counted */
	phi += 1.0/r;
	VVS(f, -= b[i].pos, / (r*r*r));
	phi *= gm;
	b[i].phi += phi;
	VS(f, *= -gm);		/* depends on sign convention of ewald */
	VV(b[i].acc, += f);
    }
    StopTimer(&FixCubeTm);
}

#else

static void 
FixCube(body *b, int nobj, float l, float gm)
{
    int i;
    float fac;
    float f[NDIM];
    float x,y,z;

    StartTimer(&FixCubeTm);

    fac = gm/(8.0*l*l*l);
    l *= 3.0;
    for (i = 0; i < nobj; i++) {
	x = b[i].pos[0]/l;
	y = b[i].pos[1]/l;
	z = b[i].pos[2]/l;
	f[0] = x*(1.5396007178390*(y*y+z*z) 
		  - 0.64150029909958*y*y*z*z
		  + 0.1069167165166*x*x*(y*y+z*z)
		  - 1.0264004785593*x*x - 0.021383343303322*x*x*x*x
		  + 0.05345835825829837*(y*y*y*y+z*z*z*z));
	f[1] = y*(1.5396007178390*(x*x+z*z) 
		  - 0.64150029909958*x*x*z*z
		  + 0.1069167165166*y*y*(x*x+z*z)
		  - 1.0264004785593*y*y - 0.021383343303322*y*y*y*y
		  + 0.05345835825829837*(x*x*x*x+z*z*z*z));
	f[2] = z*(1.5396007178390*(x*x+y*y) 
		  - 0.64150029909958*x*x*y*y
		  + 0.1069167165166*z*z*(x*x+y*y)
		  - 1.0264004785593*z*z - 0.021383343303322*z*z*z*z
		  + 0.05345835825829837*(x*x*x*x+y*y*y*y));
	VS(f, *= l);
	VS(f, *= fac);
	VV(b[i].acc, -= f);
    }
    StopTimer(&FixCubeTm);
}
#endif

#define N 2
#ifndef M_PI
#define	M_PI	3.14159265358979323846
#endif

/* sign convention used here assumes attractive force */
void
ewald(float *x, float L, float *phip, float *f)
{
    int i[NDIM];
    double xx[NDIM];
    double ax2;
    double axx;
    const double alpha = 2.0/L;
    const double a1 = 2.0*alpha/sqrt(M_PI);
    const double a2 = -M_PI*M_PI/(alpha*alpha*L*L);
    const double a3 = 2.0*M_PI/L;
    const double a4 =  2.0/(L*L);
    const double a5 =  1.0/(M_PI*L);
    double hh;
    double t, t1, t2;
    double phi;

    VS(f, = 0);
    phi = M_PI/(alpha*alpha*L*L*L);
    for (i[0] = -N; i[0] <= N; i[0]++) {
	for (i[1] = -N; i[1] <= N; i[1]++) {
	    for (i[2] = -N; i[2] <= N; i[2]++) {
		VVV(xx, = x, - L * i);
		ax2 = Dot(xx,xx);
		axx = sqrt(ax2);
		if (axx > 2.6*L) continue;
		t1 = erfc(alpha*axx);
		t2 = a1*exp(-alpha*alpha*ax2);
		t = (t1/axx+t2)/ax2;
		phi -= t1/axx;
		VV(f, += t * xx);
	    }
	}
    }
    for (i[0] = -N; i[0] <= N; i[0]++) {
	for (i[1] = -N; i[1] <= N; i[1]++) {
	    for (i[2] = -N; i[2] <= N; i[2]++) {
		if (!i[0] && !i[1] && !i[2]) continue;
		hh = Dot(i,i);
		if (hh > 8) continue;
		t = a5*exp(a2*hh)*cos(a3*Dot(i,x))/hh;
		phi -= t;
		t = a4*exp(a2*hh)*sin(a3*Dot(i,x))/hh;
		VV(f, += t * i);
	    }
	}
    }
    *phip = phi;
}

#if 0
float erffc(float x)
{
    return x < 0.0 ? 1.0+gammp(0.5,x*x) : gammq(0.5,x*x);
}

float gammp(float a, float x)
{
    float gamser,gammcf,gln;
    
    if (x < 0.0 || a <= 0.0) Error("Invalid arguments in routine gammp");
    if (x < (a+1.0)) {
	gser(&gamser,a,x,&gln);
	return gamser;
    } else {
	gcf(&gammcf,a,x,&gln);
	return 1.0-gammcf;
    }
}


float gammq(float a, float x)
{
    float gamser,gammcf,gln;

    if (x < 0.0 || a <= 0.0) Error("Invalid arguments in routine gammq");
    if (x < (a+1.0)) {
	gser(&gamser,a,x,&gln);
	return 1.0-gamser;
    } else {
	gcf(&gammcf,a,x,&gln);
	return gammcf;
    }
}

#define ITMAX 100
#define EPS 3.0e-7

void gser(float *gamser, float a, float x, float *gln)
{
    int n;
    float sum,del,ap;

    *gln=gammln(a);
    if (x <= 0.0) {
	if (x < 0.0) Error("x less than 0 in routine gser");
	*gamser=0.0;
	return;
    } else {
	ap=a;
	del=sum=1.0/a;
	for (n=1;n<=ITMAX;n++) {
	    ++ap;
	    del *= x/ap;
	    sum += del;
	    if (fabs(del) < fabs(sum)*EPS) {
		*gamser=sum*exp(-x+a*log(x)-(*gln));
		return;
	    }
	}
	Error("a too large, ITMAX too small in routine gser");
	return;
    }
}


#define FPMIN 1.0e-30

void gcf(float *gammcf, float a, float x, float *gln)
{
    int i;
    float an,b,c,d,del,h;

    *gln=gammln(a);
    b=x+1.0-a;
    c=1.0/FPMIN;
    d=1.0/b;
    h=d;
    for (i=1;i<=ITMAX;i++) {
	an = -i*(i-a);
	b += 2.0;
	d=an*d+b;
	if (fabs(d) < FPMIN) d=FPMIN;
	c=b+an/c;
	if (fabs(c) < FPMIN) c=FPMIN;
	d=1.0/d;
	del=d*c;
	h *= del;
	if (fabs(del-1.0) < EPS) break;
    }
    if (i > ITMAX) Error("a too large, ITMAX too small in gcf");
    *gammcf=exp(-x+a*log(x)-(*gln))*h;
}


float gammln(float xx)
{
    double x,y,tmp,ser;
    static double cof[6]={
	76.18009172947146,-86.50532032941677,24.01409824083091,
	-1.231739572450155,0.1208650973866179e-2,-0.5395239384953e-5};
    int j;
    
    y=x=xx;
    tmp=x+5.5;
    tmp -= (x+0.5)*log(tmp);
    ser=1.000000000190015;
    for (j=0;j<=5;j++) ser += cof[j]/++y;
    return -tmp+log(2.5066282746310005*ser/x);
}
#endif

static void
acc_zero(body *btab, int nobj, float mtot) 
{
    body *p;
    double force[NDIM];

    VS(force, = 0.0);
    for (p = btab; p < btab+nobj; p++) {
	VV(force, += p->mass*p->acc);
    }
    MPMY_Combine(force, force, NDIM, MPMY_DOUBLE, MPMY_SUM);
    VS(force, /= mtot);
    singlPrintf("correcting Acc (%g,%g,%g)\n", force[0], force[1], force[2]);
    for (p = btab; p < btab+nobj; p++) {
	VV(p->acc, -= force);
    }
}

static int 
maxheap(void)
{
    int memused = malloc_heapsz()/1024;
    MPMY_Combine(&memused, &memused, 1, MPMY_INT, MPMY_MAX);
    return memused;
}

static int 
maxmem(void)
{
    int memused = malloc_used()/1024;
    MPMY_Combine(&memused, &memused, 1, MPMY_INT, MPMY_MAX);
    return memused;
}

void CosmoPush(struct cosmo_s *p, float time){
    float Omega0 = p->Omega0;
    float Lambda = p->Lambda;
    float H0 = p->H0;
    float H, a2, a3, aold, a2dot;
    int i;
    float deltat, dt;
    int nstep;

    /* The cosmo structure holds,H0, Omega0, Lambda' = Lambda/3H0^2, a
       and t.  We integrate (forward or backward) to the new 'time' */

    deltat = time - p->t;
    if( deltat == 0.F )
	return;

    /* Felten et al do all their integrals with dt=1/(400 H0).  We can
       do the same by choosing Nstep appropriately.  In fact, we can
       do a little better by ensuring dt < 1/(400 H). */
    aold = p->a;
    a2 = aold*aold;
    a3 = a2*aold;
    H = H0 * sqrt( Omega0/a3 + Lambda - (Omega0 + Lambda - 1.F)/a2);
    nstep = (int)(400.*H*fabs(deltat)) + 1;
    Msgf(("Cosmo push %d steps, deltat=%g, H*deltat=%g\n", 
	  nstep, deltat, deltat*H));
    dt = deltat/(float)nstep;
    
    for(i=0; i<nstep; i++){
	aold = p->a;
	a2 = aold*aold;
	a3 = a2*aold;
	H = H0 * sqrt( Omega0/a3 + Lambda - (Omega0 + Lambda - 1.F)/a2);
	/* Follow the advice of Felten et al.  Do this to second-order */
#if 1
	a2dot = aold*H0*H0*( -0.5F*Omega0/a3 + Lambda);
#else
	a2dot = 0.F;
#endif
	p->a = aold + dt*H*aold + 0.5F*dt*dt*a2dot;
    }
    Msgf(("After push Z=%g\n", 1./p->a - 1.));
    p->t = time;
}


#ifdef _AIX
static void
PrintLoad(Timer_t *cpu, Timer_t *wc)
{
    int i;
    float percent;
    float *table;

    percent = 100.0 * ReadTimer(cpu)/ReadTimer(wc);

    table = Malloc(MPMY_Nproc() * sizeof(float));
    MPMY_Gather(&percent, 1, MPMY_FLOAT, table, 0);

    singlPrintf("CPU Efficiency: ");
    for (i = 0; i < MPMY_Nproc(); i++)
      singlPrintf("%4.0f ", table[i]);
    singlPrintf("\n");
    Free(table);
}
#endif
