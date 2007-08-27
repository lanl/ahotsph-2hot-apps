/*
 * Copyright 1994 Michael S. Warren and John K. Salmon.  All Rights Reserved.
 */

/* Posix says I need these to call open...Life is so complicated.*/
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
/* End of requirements for open */
#include <stdio.h>		/* only use sprintf */
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stddef.h>
#include <unistd.h>
#include <errno.h>
#include "Assert.h"
#include "SDF.h"
#include "protos.h"
#include "macr.h"
#include "malloc.h"
#include "bigmalloc.h"
#include "SDFwrite.h"
#include "SDFread.h"
#include "physics_sph.h"
#include "vop.h"
#include "Msgs.h"
#include "tree.h"  /* includes timers.h and pqsort.h */
#include "singlio.h"
#include "mpmy.h"
#include "gc.h"
#include "files.h"
#include "getparam.h"
#include "verify.h"
#include "integrate.h"
#include "fastflpt.h"  /* Karen */

#define MACargs float x0, float x1, float x2, cell *, float , Stk *, int *

int MAC_n(MACargs);
int MACbc(MACargs);
int MAC_b(MACargs, int);
int MACsph_c(MACargs);
int MACsph_b(MACargs, int);
static void SanityCheck(bodyptr btab, int nobj, int gnobj, double *mtotp);
static void SetLBTarget(sortresult_t *decompp, int hetero_load_balance);
static void ConvertVPofT(body *xptr, int n, float dp_dt, int to_physical);
static void collision_init(body **btab, int *nobj, int *gnobj, SDF *csdfp);
static void point_mass_offset(body *btab, int nobj,SDF *csdfp);
static void read_point_mass(body *point_mass, SDF *csdfp);
static void hpsel(unsigned long m, unsigned long n, float arr[], float heap[]);


void update_intermediate(body *btab, int nobj, float dt_last, int flag);
void update_final(body *btab, int nobj);
void update_point_mass(body *btab, int nobj, body *p);

static void ConvertVPofT(body *xptr, int n, float dp_dt, int to_physical);

static float Znow(float time);

extern void malloc_print(void);
extern int malloc_debug(int);
void body_interaction(body **pp, body **end, const float *pos0, float *mass0,
		      float *phi0, float *acc0, moment *qpole0, 
		      const float *eps2p, int *ncut, int *tot_interact);

static float Znow(float time);
static void set_vels(body *p, int n, float real_time);
static SDF *startup(int argc, char **argv);

Timer_t StepTot, StepTotWC, BuildTot;
Timer_t GravTm;
Timer_t FindForcesTm;
extern Timer_t SDFreadTm, SDFwriteTm;	/* should be in SDFread.h? */
Counter_t NbodyCnt;
Counter_t NtermsCnt;

extern Counter_t CCInt, CBInt, BCInt, BBInt;
extern Counter_t CCIntRej;
extern Counter_t TranslateCnt;


Timer_t RhoSPH, ForceSPH;
Counter_t NbodCnt;
extern Counter_t SPHCnt, SPHrej, nbrMACCnt;

int
main(int argc, char *argv[])
{
    int gnobj, nobj;
    bodyptr btab;
    float newton_const;
    float eps;			/* Plummer smoothing length */
    float tol;			/* MAC tolerance */
    int absolute_tol;
		/* for big MAC, this is multiplied by M/(rsize*rsize) */
    int i;
    float courant_number;
    float rmin[NDIM], rmax[NDIM];
    float sysradius;
    float dt;
    float dt_last = (float)0.0;
    int nsteps;
    int first_step = 1;
    int do_output;
    int hetero_load_balance;
    int timer_freq;
    int output_freq;
    int iter = 0;
    bodyptr p = 0;
    int stride = sizeof(body)/sizeof(float);
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
    float hubble;
    int cosmology = 0;
    int do_gravity = 1;
    int no_hydro = 0;
    int save_first;		/* save first step (for acc testing) */
    int loadbalance = 1;	/* use nterms to load balance */
    int id_sort_output;
    float new_h;
    float new_u;
    int ndim = NDIM;
    int boundary = 0;
    float max_vsound;
    float max_rho;
    float min_dt;
    float dt_lower_limit, dt_upper_limit; /* Added upper limit - Don */
    int min_nbrs, max_nbrs;
    int do_point_mass = 0;
    int do_setup;
    body point_mass;
    float min_h, max_h;
    int atmin, atmax;
    int adaptive_dt;
    float visc_alpha, visc_beta;

    double pe, ke, te;
    double mtot = 0.0;
    sortresult_t sortedbtab;
    tree_t thetree;
    int massconf, xconf, yconf, zconf;
    int vxconf, vyconf, vzconf;
    int identconf, idconf, ntermsconf;
    int uconf, hconf;
    int rhoconf; /* GR */
    char name[256];
    int exact_rho = 0;
    double nterms_total;
    MPMY_Comm_request req;
    int setup_collision;
    double com[NDIM], comv[NDIM];
    float max_dacc, rms_dacc;
    float dacc_factor;
    float density_max, kelvin_max, press_max; /* GR */
    float xx0,yy0,zz0,vx0,vy0,vz0, bhmass, Gamma;
    int do_kerr;
    float kerr_ang_mom;
    char statsname[256];
    FILE *statsfile = 0;
    int nbrcut_max, nbrcut_min;
    float nbrcut_fac;
    float hole_mass;	/* Don - Set to zero to remove the black hole */
    float hp, hx;	/* Don - Polarization states of the grav_rad. */
    int remove_com_motion; /* Karen */
    float vt, vx, vy, vz, ut, ux, uy, uz, v2, uut, gamaold; /* Karen */
    float r_horizon;

    MPMY_Init(&argc, &argv);
    singlPrintf("Welcome to the O(N) sph integrator\n");

#if 0 && defined(__DELTA__)
    free(malloc(11000000));
#endif

    csdfp = startup(argc, argv);
    SDFgetstring(csdfp, "datafile", name, sizeof(name));
    if( !((strncmp(name, "test", 4) == 0)) ) {
	int junk;
	sdfp = SDFread(csdfp, (void **)&btab, &gnobj, &nobj, sizeof(body),
		      "mass", offsetof(body, mass), &massconf,
		      "x", offsetof(body, pos[0]), &xconf,
		      "y", offsetof(body, pos[1]), &yconf,
		      "z", offsetof(body, pos[2]), &zconf,
		      "vx", offsetof(body, vel[0]), &vxconf,
		      "vy", offsetof(body, vel[1]), &vyconf,
		      "vz", offsetof(body, vel[2]), &vzconf,
		      "u", offsetof(body, u), &uconf,
		      "h", offsetof(body, h), &hconf,
		      "ident", offsetof(body, ident), &identconf,
		      "id", offsetof(body, ident), &idconf,
		      "nterms", offsetof(body, nterms), &ntermsconf,
		      "rho", offsetof(body, rho), &rhoconf,
		      "gama", offsetof(body, gama), &junk,
		      "enthalpy", offsetof(body, enth), &junk,
		      "sx", offsetof(body, mom[0]), &junk,
		      "sy", offsetof(body, mom[1]), &junk,
		      "sz", offsetof(body, mom[2]), &junk,
		      NULL);

	Msg_do("Data read, nobj=%d, gnobj=%d\n", nobj, gnobj);
	Msg_do("Nproc:%d, Procnum: %d\n", MPMY_Nproc(), MPMY_Procnum());
	if( identconf && idconf ){
	    Error("You can't have both an 'id' and an 'ident' in the data!\n");
	}
	if( massconf==0 || xconf==0 || yconf==0 || zconf==0 ){
	    Error("Could not find %s %s %s %s in data file!\n",
		  (massconf==0)? "mass" : "",
		  (xconf==0)? "x" : "",
		  (yconf==0)? "y" : "",
		  (zconf==0)? "z" : "");
	}
	if( vxconf != vyconf || vxconf != vzconf ){
	    Error("Found only some velocity components!\n");
	}
	if( identconf == 0 && idconf == 0 ){
	    SinglWarning("No \"ident\" in file, numbering sequentially\n");
	    FixId(btab, nobj, gnobj);
	}
	if( ntermsconf == 0 ){
	    SinglWarning("No \"nterms\" in file, using 1\n");
	    FixNterms(btab, nobj);
	}
	if( rhoconf==0){
	    Error("Missing rho datafile\n");
	}
	SDFgetfloatOrDefault(sdfp, "Gnewt", &newton_const, (float)1.0);
	if( SDFhasname("time", sdfp) )
	  SDFgetfloatOrDefault(sdfp, "time",  &tpos, (float)0.0);
	else
	  SDFgetfloatOrDefault(sdfp, "tpos",  &tpos, (float)0.0);

    /* This doesn't work if there is roundoff error in tvel */
    /* SDFgetfloatOrDefault(sdfp, "tvel",  &tvel, tpos);*/
    /* Fix the tpos == tvel line in Integrate() if you don't like next line */
	tvel = tpos;

	SDFgetfloatOrDefault(sdfp, "hubble",  &hubble, (float)0.0);
	SDFgetfloatOrDefault(sdfp, "Gamma",  &Gamma, (float)(5.0/3.0));
	SDFgetintOrDefault  (sdfp, "iter",  &iter, 0);
	SDFgetintOrDefault  (sdfp, "ndim",  &ndim, 3);
	SDFgetintOrDefault  (sdfp, "boundary",  &boundary, 0);
	SDFgetintOrDefault  (csdfp, "do_point_mass", &do_point_mass, 0);
	SDFgetintOrDefault  (csdfp, "do_setup",  &do_setup, 0);
	SDFgetintOrDefault  (csdfp, "do_kerr",  &do_kerr, 0);
	SDFgetintOrDefault(csdfp, "remove_com_motion", &remove_com_motion, 0);
        SDFgetfloatOrDefault(csdfp, "hole_mass",  &hole_mass, (float)1.0); /* Don */
	SDFgetfloatOrDefault(csdfp, "kerr_ang_mom",  &kerr_ang_mom, 0.0);
	SDFgetfloatOrDefault(csdfp, "bhmass",  &bhmass, (float)1e3);
	    singlPrintf("float bhmass = %g;\n", bhmass);
	    singlPrintf("int do_kerr = %d;\n", do_kerr);
	    singlPrintf("float hole_mass = %f;\n", hole_mass);
	    singlPrintf("float kerr_ang_mom = %f;\n", kerr_ang_mom);
	    singlPrintf("int remove_com_motion = %d;\n", remove_com_motion);

	setup_metric(do_kerr, hole_mass, kerr_ang_mom);

	if (do_setup) {
	    SDFgetfloatOrDefault(csdfp, "xx0",  &xx0, (float)100.0);
	    SDFgetfloatOrDefault(csdfp, "yy0",  &yy0, (float)0.0);
	    SDFgetfloatOrDefault(csdfp, "zz0",  &zz0, (float)0.0);
	    SDFgetfloatOrDefault(csdfp, "vx0",  &vx0, (float)-0.1);
	    SDFgetfloatOrDefault(csdfp, "vy0",  &vy0, (float)4.0e-2);
	    SDFgetfloatOrDefault(csdfp, "vz0",  &vz0, (float)0.0);
	    singlPrintf("float xx0 = %f;\n", xx0);
	    singlPrintf("float yy0 = %f;\n", yy0);
	    singlPrintf("float zz0 = %f;\n", zz0);
	    singlPrintf("float vx0 = %f;\n", vx0);
	    singlPrintf("float vy0 = %f;\n", vy0);
	    singlPrintf("float vz0 = %f;\n", vz0);
	    tvel = tpos = 0.0;
	    initial_cond(btab,nobj,xx0,yy0,zz0,vx0,vy0,vz0,bhmass,Gamma);
	} 

	if (do_point_mass) {
	    if (iter)
	      read_point_mass(&point_mass, sdfp);
	    else
	      read_point_mass(&point_mass, csdfp);
	}
	if(sdfp)
	  SDFclose(sdfp);
    } else {
	int seed, cencon;
	singlPrintf("Generating random dataset\n");
	if( SDFgetint(csdfp, "nobj", &gnobj) )
	  Error("Sorry, you've got to have an \"nobj\"\n");
	SDFgetintOrDefault(csdfp, "seed", &seed, 123);
	SDFgetintOrDefault(csdfp, "cencon", &cencon, 0);
	RdTest(&btab, gnobj, &nobj, seed, cencon);
	singlPrintf("int seed = %d;\n", seed);
	singlPrintf("int cencon = %d;\n", cencon);
	hconf = uconf = 1;
	newton_const = (float)1.0;
	tvel = tpos = (float)0.0;
	iter = 0;
    }

    SDFgetfloatOrDie(csdfp, "epsilon", &eps);
    SDFgetfloatOrDie(csdfp, "errtol", &tol);
    SDFgetfloatOrDefault(csdfp, "frac_tol", &frac_tol, 0.0);
    SDFgetintOrDefault(csdfp, "absolute_tol", &absolute_tol, 0);
    SDFgetfloatOrDie(csdfp, "dt", &dt);
    SDFgetfloatOrDefault(csdfp, "new_h", &new_h, (float)0.0);
    SDFgetfloatOrDefault(csdfp, "new_u", &new_u, (float)0.0);
    SDFgetfloatOrDefault(csdfp, "min_h", &min_h, (float)0.0);
    SDFgetfloatOrDefault(csdfp, "max_h", &max_h, (float)1e30);
    SDFgetintOrDefault(csdfp, "nbrcut_max", &nbrcut_max, 400);
    SDFgetintOrDefault(csdfp, "nbrcut_min", &nbrcut_min, 8);
    SDFgetfloatOrDefault(csdfp, "nbrcut_fac", &nbrcut_fac, (float)0.1);
    SDFgetfloatOrDefault(csdfp, "courant_number", &courant_number, (float)0.4);
    SDFgetfloatOrDefault(csdfp, "visc_alpha", &visc_alpha, (float)1.0);
    SDFgetfloatOrDefault(csdfp, "visc_beta", &visc_beta, (float)2.0);
    SDFgetfloatOrDefault(csdfp, "dt_lower_limit", &dt_lower_limit, (float)0.0);
    SDFgetfloatOrDefault(csdfp, "dt_upper_limit", &dt_upper_limit, (float)1000.0); /* Don */
    SDFgetfloatOrDefault(csdfp, "dacc_factor", &dacc_factor, (float)1e-10);
    SDFgetintOrDie(csdfp, "nsteps", &nsteps);
    SDFgetintOrDefault(csdfp, "cosmology", &cosmology, 0);
    SDFgetintOrDefault(csdfp, "do_gravity", &do_gravity, 1);
    SDFgetintOrDefault(csdfp, "no_hydro", &no_hydro, 0);
    SDFgetintOrDefault(csdfp, "log_time", &log_time, 0);
    SDFgetintOrDefault(csdfp, "comov_eps", &comov_eps, 0);
    SDFgetintOrDefault(csdfp, "setpvel", &setpvel, 0);
    SDFgetintOrDefault(csdfp, "save_first", &save_first, 0);
    SDFgetintOrDefault(csdfp, "exact_rho", &exact_rho, 0);
    SDFgetintOrDefault(csdfp, "loadbalance", &loadbalance, 1);
    SDFgetintOrDefault(csdfp, "id_sort_output", &id_sort_output, 1);
    SDFgetintOrDefault(csdfp, "adaptive_dt", &adaptive_dt, 1);
    SDFgetintOrDefault(csdfp, "setup_collision", &setup_collision, 0);
    SDFgetfloatOrDefault(csdfp, "r_horizon", &r_horizon, -1e30);
    if (setup_collision) {
	if (iter)
	  Error("You are trying to init collision with iter != 0\n");
	collision_init(&btab, &nobj, &gnobj, csdfp);
	}
    if (do_point_mass && iter == 0) {
	point_mass_offset(btab, nobj, csdfp);
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
    SDFgetintOrDefault(csdfp, "hetero_load_balance", &hetero_load_balance, 0);
    SDFgetstring(csdfp, "statsfile", statsname, sizeof(statsname));
    if (MPMY_Procnum() == 0) Fopen(statsfile, statsname, "a");
    singlPrintf("Appending stats to %s\n", statsname);
    if(csdfp) 
	SDFclose(csdfp);

    if( Msg_test("memleak") ){
	malloc_debug(2);
	singlPrintf("Malloc_debug(2), expect slow mallocs\n");
    }else{
	malloc_debug(1);
    }
    singlPrintf("float errtol = %g;\n", tol);
    singlPrintf("float dt = %g;\n", dt);
    singlPrintf("float courant_number = %g;\n", courant_number);
    singlPrintf("float visc_alpha = %g;\n", visc_alpha);
    singlPrintf("float visc_beta = %g;\n", visc_beta);
    singlPrintf("float epsilon = %g;\n", eps);
    singlPrintf("int nsteps = %d;\n", nsteps);
    singlPrintf("int do_gravity = %d;\n", do_gravity);
    singlPrintf("int loadbalance = %d;\n", loadbalance);
    singlPrintf("int nproc = %d;\n", MPMY_Nproc());
    singlPrintf("float min_h = %g;\n", min_h);
    singlPrintf("float max_h = %g;\n", max_h);
    singlPrintf("int exact_rho = %d;\n", exact_rho);
    singlPrintf("int adaptive_dt = %d;\n", adaptive_dt);
    if (do_point_mass) {
	singlPrintf("Point mass %g at (%g,%g,%g)\n", point_mass.mass,
	       point_mass.pos[0], point_mass.pos[1], point_mass.pos[2]);
	singlPrintf("Point mass vel (%g,%g,%g)\n", 
	       point_mass.vel[0], point_mass.vel[1], point_mass.vel[2]);
    }
    
    if( do_output ){
	singlPrintf("Output to %s.nnn, every %d steps\n", 
	       outnamebase, output_freq);
    }else{
	singlPrintf("No output.\n");
    }
    if (cosmology || log_time) {
	Error("flag not supported yet\n");
    }

/* OJO: For test runs with just a few collisionless particles, 
   un-comment this line: 
    nobj = gnobj = 10;              
*/

    SanityCheck(btab, nobj, gnobj, &mtot);

    pqsortsetup(&sortedbtab, btab, nobj, sizeof(body), 0.01F, Realloc_f);

    SetupTree(&thetree, NDIM, 
	      sizeof(body), sizeof(cell), TBODYSZ, sizeof(cofmdata),
	      (pq_keyproto)GetKeyFromStruct, (pq_wgtproto)GetCost,
	      CofmFromDaugh, (cellfromcofm_t)CellFromCofm);

    if (log_time) {
	float dpdt;
	dpdt = 2./3. * pow((double)tpos, -1./3.);
	ConvertVPofT(btab, nobj, dpdt, 0);
    }

    SPH_setup(ndim);

    if (new_h != (float)0.0) {
	singlPrintf("Setting h to %f\n", new_h);
	for (p = btab; p < btab+nobj; p++) {
	    p->h = new_h;
	}
    } else if (hconf == 0) Error("No h in data file\n");

    if (new_h != (float)0.0) {
	singlPrintf("Setting h to isothermal %f\n", new_h);
	for (p = btab; p < btab+nobj; p++) {
	    float r = sqrt(Dot(p->pos, p->pos));
	    p->h = new_h*pow(r,2.0/3.0);
	    VS(p->vel, = (float)0.0);
	}
    } else if (hconf == 0) Error("No h in data file\n");

    if (new_u != (float)0.0) {
	singlPrintf("Setting u to %f\n", new_u);
	for (p = btab; p < btab+nobj; p++) {
	    p->u = new_u;
	}
    } else if (uconf == 0) Error("No u in data file\n");
    dt_last = dt;
    
    /* This only works for the 1d shock tube */
    if (boundary) {
	for (p = btab; p < btab+nobj; p++) {
	    int id_x = p->ident & ((1 << 10)-1);
	    /* these are not updated */
	    if ((id_x < boundary) || (id_x >= gnobj - boundary)) 
	      /* use bit 30 as flag, must match in integrate.c */
	      p->ident |= (1 << 30);
	    else if (p->ident & (1 << 30))
	      Error("Set bit conficts with boundary flag\n");
	}
    }

    if (frac_tol != (float)0.0) { /* init acc_last for first step */
	/* Reset timers and counters */
	ClearEnabledTimers();
	ClearEnabledCounters();
	StartTimer(&StepTot);
	StartTimer(&StepTotWC);
	FindBbox(btab, nobj, rmin, rmax);
	sysradius = 0.5*FixRsize(rmin, rmax);
	this_tol = tol*mtot/(sysradius*sysradius);
	if (comov_eps) {
	    this_eps = eps/(Znow(tpos)+(float)1.0); /* comoving smoothing */
	} else {
	    this_eps = eps;
	}
	SetTol(this_tol, frac_tol, newton_const, this_eps, gnobj);
	FixKeys(btab, nobj, GETKEY);
	singlPrintf("BuildTree, most approx. mac\n");
	BuildTree(&thetree, &sortedbtab);
	btab = sortedbtab.data;
	nobj = sortedbtab.nobj;
	singlPrintf("BuildTree done\n");

	MPMY_Sync(); /* No sync might cause msg buffer overflow? */
	singlPrintf("FindForces, this_eps=%g\n", this_eps);
	Walk(&thetree, &thetree, sizeof(Sink), (macv_t) Lowestmacv,
	     (inherit_t)InheritSink);
	singlPrintf("FindForces done\n");

	FreeTree(&thetree);
	singlPrintf("FreeTree done\n");
	Msgf(("FreeTree done\n"));
	StopTimer(&StepTot);
	StopTimer(&StepTotWC);
	OutputTimer(&StepTot, singlPrintf);
	OutputTimer(&StepTotWC, singlPrintf);

	for (p = btab; p < btab+nobj; p++) {
	    float fx;
	    p->acc_last = fabs(p->acc[0]);
	    fx = fabs(p->acc[1]);
	    if (fx > p->acc_last) p->acc_last = fx;
	    fx = fabs(p->acc[2]);
	    if (fx > p->acc_last) p->acc_last = fx;
	}
    }

    if (first_step) {
	for (p = btab; p < btab+nobj; p++) {
	    p->gama_last = p->gama;
	    if (p->ident & (1 << 30)) 
	      Error("Set bit conficts with boundary flag\n");
	}
    }


    for (nsteps += iter; iter <= nsteps; iter++) {
	/* If hetero_load_balance flag set, store load balance target value
	 * else default constant value from pqsortsetup() used */
	if (hetero_load_balance)
	    SetLBTarget(&sortedbtab, hetero_load_balance);
	/* Reset timers and counters */
	ClearEnabledTimers();
	ClearEnabledCounters();
	StartTimer(&StepTot);
	StartTimer(&StepTotWC);
	ShrinkBtab(&btab, &nobj, r_horizon);
	MPMY_Combine(&nobj, &gnobj, 1, MPMY_FLOAT, MPMY_SUM);
	pqsortsetup(&sortedbtab, btab, nobj, sizeof(body), 0.01F, Realloc_f);
	singlPrintf("Subtracted particles inside horizon, gnobj now %d\n", gnobj);
	FindBbox(btab, nobj, rmin, rmax);
	sysradius = 0.5*FixRsize(rmin, rmax);
	Msgf(("rmin=(%g, %g, %g) rmax=(%g, %g, %g)\n",
	      rmin[0], rmin[1], rmin[2], 
	      rmax[0], rmax[1], rmax[2]));
	if (absolute_tol)
	  this_tol = tol;
	else 
	  this_tol = tol*mtot/(sysradius*sysradius);
	if (comov_eps) {
	    this_eps = eps/(Znow(tpos)+(float)1.0); /* comoving smoothing */
	} else {
	    this_eps = eps;
	}

	SetTol(this_tol, frac_tol, newton_const, this_eps, gnobj);
	FixKeys(btab, nobj, GETKEY);
	/* Must be done before BuildTree */
	if (!(first_step || exact_rho)) {
	    get_metric(btab, nobj); /* GR */ 
	    update_intermediate(btab, nobj, dt_last, 1);
	} 
	if (frac_tol == 0.0)
	  singlPrintf("BuildTree, this_tol=%g\n", this_tol);
	else
	  singlPrintf("BuildTree, frac_tol=%g\n", frac_tol);
	StartTimer(&BuildTot);
	BuildTree(&thetree, &sortedbtab);
	btab = sortedbtab.data;
	nobj = sortedbtab.nobj;
	StopTimer(&BuildTot);
	singlPrintf("BuildTree done\n");
	AddCounter(&NbodyCnt, nobj);

	MPMY_Sync(); /* No sync might cause msg buffer overflow? */
	singlPrintf("FindForces\n");
	StartTimer(&FindForcesTm);
	if (do_gravity) {
	    Walk(&thetree, &thetree, sizeof(Sink), 
		 (macv_t)((frac_tol == 0.0) ? Unifiedmacv : Fracmacv),
		 (inherit_t)InheritSink);
	} else {
	    for (p = btab; p < btab+nobj; p++) {
		VS(p->acc, = (float)0.0);
		p->phi = (float)0.0;
	    }
	}
	StopTimer(&FindForcesTm);
	singlPrintf("FindForces done\n");

	if (first_step || exact_rho) {
	    MPMY_Sync(); /* No sync might cause msg buffer overflow? */
	    singlPrintf("FindRho\n");
	    StartTimer(&RhoSPH);
	    get_metric(btab, nobj); /* GR - so update_final has gr_mass - Don */
	    SetSPH(visc_alpha, visc_beta, Gamma, gnobj, macRho, nbrMAC);
	    if (!no_hydro) {
		Walk(&thetree, &thetree, sizeof(SinkSPH), 
		     (macv_t)SPHgate, (inherit_t)InheritSPH);
	    }
	    update_final(btab, nobj);
	    /* This sets rho_est and pr for communication during BuildTree */
	    /* get_metric(btab, nobj); GR - Removed - Don */ 
	    update_intermediate(btab, nobj, dt_last, 0);
	    StopTimer(&RhoSPH);

	    FreeTree(&thetree);
	    singlPrintf("FreeTree done\n");

	    singlPrintf("BuildTree\n");
	    StartTimer(&BuildTot);
	    BuildTree(&thetree, &sortedbtab);
	    btab = sortedbtab.data;
	    nobj = sortedbtab.nobj;
	    StopTimer(&BuildTot);
	}

	MPMY_Sync(); /* No sync might cause msg buffer overflow? */
	singlPrintf("ForceSPH\n");
	StartTimer(&ForceSPH);
	SetSPH(visc_alpha, visc_beta, Gamma, gnobj, macSPH, nbrMAC);
	if (!no_hydro) {
	    Walk(&thetree, &thetree, sizeof(SinkSPH), 
		 (macv_t)SPHgate, (inherit_t)InheritSPH);
	}
	update_final(btab, nobj);
	StopTimer(&ForceSPH);
	FreeTree(&thetree);

	add_gr(btab, nobj);	/* GR */

	if (setpvel) {
	    setpvel = 0;
	    set_vels(btab, nobj, tpos);
	    singlPrintf("Velocities adjusted to Zel'dovich approximation.\n");
	    if (log_time) {
		float dpdt;
		dpdt = 2./3. * pow((double)tpos, -1./3.);
		ConvertVPofT(btab, nobj, dpdt, 0);
	    }
	}

	if (do_point_mass) {
	    p = &point_mass;
	    update_point_mass(btab, nobj, p);
	    MPMY_ICombine_Init(&req);
	    MPMY_ICombine(&(p->phi), &(p->phi), 1, MPMY_FLOAT, MPMY_SUM, req);
	    MPMY_ICombine(&(p->acc), &(p->acc), NDIM, MPMY_FLOAT,MPMY_SUM,req);
	    MPMY_ICombine_Wait(req);
	}

	if (ForceOutput()
	    || (do_output && !first_step
		&& ((iter+output_freq) % output_freq == 0))
	    || (save_first && first_step)) {
	    sortresult_t outputsort;
	    outbodyptr output_btab;
	    char outname[256];
	    int output_nobj = nobj;
	    float tpos_out = tpos;
	    float tvel_out = tvel; /* changed in Integrate() */

	    pe = ke = te = 0.0;
	    for (p = btab; p < btab+nobj; p++) {
		/* ke += (float)0.5 * p->mass * Dot(p->vel, p->vel); */
		ke += p->mass * ( (p->gama*p->alfa) * (p->gama*p->alfa) 
                    - (float)1.0 ) / (float)2.0;
		te += p->mass * p->u;
		pe += (float)0.5 * p->mass * p->phi;
	    }
	    output_btab = Malloc(output_nobj * sizeof(outbody));
	    for(i=0; i<output_nobj; i++){
		output_btab[i].mass = btab[i].mass;
		VV(output_btab[i].pos, = btab[i].pos);
		VV(output_btab[i].vel, = btab[i].vel);
		output_btab[i].u = btab[i].u;
		output_btab[i].h = btab[i].h;
		output_btab[i].rho = btab[i].rho;
		output_btab[i].phi = btab[i].phi;
		output_btab[i].nbrs = btab[i].nbrs;
		if (boundary)
		  output_btab[i].ident = btab[i].ident & ~(1<<30);
		else
		  output_btab[i].ident = btab[i].ident;
		/* GR */
		VV(output_btab[i].mom, = btab[i].mom);
		output_btab[i].mom[3] = btab[i].mom[3];
		output_btab[i].gama = btab[i].gama;
		output_btab[i].enth = btab[i].enth;
		/* end GR */
	    }
	    Msg("output", ("Doing output of %d bodies\n", output_nobj));
	    if (id_sort_output) {
		singlPrintf("Trying to sort output\n");
		pqsortsetup_order(&outputsort, output_btab, output_nobj,
				  sizeof(outbody), 0.1F, 1, Realloc_f);
		output_btab = pqsort(&outputsort,
				     (pq_wgtproto)UnityCost, 
				     (pq_keyproto)OutIdentKey);
		output_nobj = outputsort.nobj;
	    } else {
		singlPrintf("Output not sorted!\n");
	    }
	    Msg("output", ("After pqsort, %d outbodies\n", output_nobj));
	    MPMY_ICombine_Init(&req);
	    MPMY_ICombine(&ke, &ke, 1, MPMY_DOUBLE, MPMY_SUM, req);
	    MPMY_ICombine(&pe, &pe, 1, MPMY_DOUBLE, MPMY_SUM, req);
	    MPMY_ICombine(&te, &te, 1, MPMY_DOUBLE, MPMY_SUM, req);
	    MPMY_ICombine_Wait(req);
	    if (do_point_mass) {
		p = &point_mass;
		ke += (float)0.5 * p->mass * Dot(p->vel, p->vel);
		pe += (float)0.5 * p->mass * p->phi;
	    }
	    sprintf(outname, "%s.%04d", outnamebase, iter);
	    if (!do_point_mass) {
		SDFwrite(outname, gnobj, 
		       output_nobj, output_btab, sizeof(outbody),
		       OUTBODYDESC,
		       "npart", SDF_INT, gnobj,
		       "eps", SDF_FLOAT, eps,
		       "Gnewt", SDF_FLOAT, newton_const,
		       "Gamma", SDF_FLOAT, Gamma,
		       "tolerance", SDF_FLOAT, this_tol*newton_const,
		       "iter", SDF_INT, iter,
		       "ndim", SDF_INT, ndim,
		       "boundary", SDF_INT, boundary,
		       "tpos", SDF_FLOAT, tpos_out,
		       "tvel", SDF_FLOAT, tvel_out,
		       "ke", SDF_DOUBLE, ke,
		       "pe", SDF_DOUBLE, pe,
		       "te", SDF_DOUBLE, te,
		       NULL);
	    } else {
		SDFwrite(outname, gnobj, 
		       output_nobj, output_btab, sizeof(outbody),
		       OUTBODYDESC,
		       "npart", SDF_INT, gnobj,
		       "eps", SDF_FLOAT, eps,
		       "Gnewt", SDF_FLOAT, newton_const,
		       "tolerance", SDF_FLOAT, this_tol*newton_const,
		       "iter", SDF_INT, iter,
		       "ndim", SDF_INT, ndim,
		       "boundary", SDF_INT, boundary,
		       "tpos", SDF_FLOAT, tpos_out,
		       "tvel", SDF_FLOAT, tvel_out,
		       "ke", SDF_DOUBLE, ke,
		       "pe", SDF_DOUBLE, pe,
		       "te", SDF_DOUBLE, te,
		       "do_point_mass", SDF_INT, do_point_mass,
		       "pt_x", SDF_FLOAT, point_mass.pos[0],
		       "pt_y", SDF_FLOAT, point_mass.pos[1],
		       "pt_z", SDF_FLOAT, point_mass.pos[2],
		       "pt_vx", SDF_FLOAT, point_mass.vel[0],
		       "pt_vy", SDF_FLOAT, point_mass.vel[1],
		       "pt_vz", SDF_FLOAT, point_mass.vel[2],
		       "pt_mass", SDF_FLOAT, point_mass.mass,
		       "pt_h", SDF_FLOAT, point_mass.h,
		       NULL);
	    }
	    Free(output_btab);
	    singlPrintf("\nOutput done.\n");
	}

	if (ForceStop()) {
	    singlPrintf("Stopping.\n");
	    break;
	}

	if (first_step) {
	    for (i = 0; i < nobj; i++) {
		VV(btab[i].vel_last, = btab[i].vel);
		VV(btab[i].force_last, = btab[i].acc);
		btab[i].udot_last = btab[i].udot;
	    }
#if 0
	    if (do_point_mass) {
		VVV(point_mass.pos_last, 
		    = point_mass.pos, - dt * point_mass.vel);
	    }
#endif
	}

	max_dacc = rms_dacc = 0.0;
	for (p = btab; p < btab+nobj; p++) {
	    float fx, dacc, sacc;
	    fx = fabs(p->acc[0]-p->force_last[0]);
	    sacc = fx;
	    fx = fabs(p->acc[1]-p->force_last[1]);
	    if (fx > sacc) sacc = fx;
	    fx = fabs(p->acc[2]-p->force_last[2]);
	    if (fx > sacc) sacc = fx;
	    dacc = sacc/(dacc_factor);
	    if (dacc > max_dacc)
	      max_dacc = dacc;
	    rms_dacc += dacc*dacc;
	    p->acc_last = sacc;
	}

	ABUpdateXs(&btab[0].u, stride, &btab[0].udot, stride, 
		  &btab[0].udot_last, stride, 
		  &btab[0].ident, stride, nobj, dt, dt_last);
	ABUpdateVs(btab[0].mom, stride, btab[0].acc, stride, 
		   btab[0].force_last, stride, 
		   &btab[0].ident, stride, nobj, dt, dt_last);
	ABUpdateVs(btab[0].pos, stride, btab[0].vel, stride, 
		   btab[0].vel_last, stride, 
		   &btab[0].ident, stride, nobj, dt, dt_last);
	UpdateSXs(&btab[0].h, stride, &btab[0].hdot, stride,
		 &btab[0].ident, stride, nobj, dt, dt_last);
#if 0
	if (do_point_mass) {
	    PUpdateV(point_mass.vel, stride, point_mass.pos, stride, 
		     point_mass.pos_last, stride, point_mass.acc, stride, 
		     1, dt, dt_last);
	    PUpdateX(point_mass.pos, stride, point_mass.pos_last, stride,
		     point_mass.acc, stride, 1, dt, dt_last);
	}
#endif
	tpos += dt;
	tvel += dt;
	dt_last = dt;

	for (p = btab; p < btab+nobj; p++) {
	    if (p->nbrs > 8*nbrcut_max) 
	      p->h *= 0.75;
	    else if (p->nbrs > 2*nbrcut_max) 
	      p->h *= pow((double)nbrcut_max/p->nbrs, (1./3.));
	    else if (p->nbrs > nbrcut_max) p->h -= nbrcut_fac * p->h;
	    else if (p->nbrs < nbrcut_min) p->h += nbrcut_fac * p->h;
	}

	ke = pe = te = 0.0;
	max_vsound = (float)0.0;
	max_rho = (float)0.0;
	min_dt = 1e30;
	max_nbrs = 0;
	min_nbrs = 10000;
	density_max = kelvin_max = press_max = (float)0.0; /* GR */
	nterms_total = 0.0;
	VS(com, = 0.0);
	VS(comv, = 0.0);
#if 0
	max_dacc = 0.0;
	rms_dacc = 0.0;
#endif
	for (p = btab; p < btab+nobj; p++) {
	    float dti;
	    const float gm1 = Gamma-(float)1.0;	/* GR */
	    float rho_lab;

	    atmax = 0;
	    atmin = 0;
	    /* ke += (float)0.5 * p->mass * Dot(p->vel, p->vel); */
	    ke += p->mass * ( (p->gama*p->alfa) * (p->gama*p->alfa) 
                    - (float)1.0 ) / (float)2.0;
	    te += p->mass * p->u;
	    pe += (float)0.5 * p->mass * p->phi;
	    if (p->vsound/p->h > max_vsound)
	      max_vsound = p->vsound/p->h;
	    dti = (float)courant_number
	      /(fabs(p->drho_dt)/p->rho + p->vsound/p->h);
	    if (dti < min_dt)
	      min_dt = dti;
	    if (p->nbrs > max_nbrs)
	      max_nbrs = p->nbrs;
	    if (p->nbrs < min_nbrs)
	      min_nbrs = p->nbrs;
	    if (p->rho > max_rho) max_rho = p->rho;
	    if (loadbalance) {
		p->nterms += 4 * p->nbrs;
	    }
	    nterms_total += p->nterms;
	    AddCounter(&SPHCnt, p->nbrs);
	    if (p->h > max_h) {
		p->h = max_h;
		atmax++;
	    }
	    if (p->h < min_h) {
		p->h = min_h;
		atmin++;
	    }
	    if (p->u < 0) {
		SeriousWarning("Negative u (%g) for particle %d\n",
			       p->u, p->ident);
		p->u = 0.0;
	    }

/* OJO: For center of mass info use this: */
	    VV(com, += p->mass*p->pos);
	    VV(comv, += p->mass*p->vel);

/* OJO: For first particle info use this: */
/*
	    if (p->ident == 0) {
		VV(com, = p->pos);
		VV(comv, = p->vel);
	    }
*/
	    /* GR */

	    rho_lab = p->rho * p->alfa / p->gama;
	    if (rho_lab > density_max) density_max = rho_lab;
	    if (p->u*gm1 > kelvin_max) kelvin_max = p->u * gm1;
	    if (p->pr >  press_max) press_max = p->pr;

	    /* end GR */
	}
        grav_rad(btab,nobj,&hp,&hx); /* grav_rad subroutine - Don */
	AddCounter(&NbodCnt, nobj);
	AddCounter(&NtermsCnt, nterms_total/1000.0); /* prevent overflow */
	MPMY_ICombine_Init(&req);
	MPMY_ICombine(&ke, &ke, 1, MPMY_DOUBLE, MPMY_SUM, req);
	MPMY_ICombine(&pe, &pe, 1, MPMY_DOUBLE, MPMY_SUM, req);
/* Collect grav_rad from all the processors - Don */
	MPMY_ICombine(&hp, &hp, 1, MPMY_FLOAT, MPMY_SUM, req); /* Don */
	MPMY_ICombine(&hx, &hx, 1, MPMY_FLOAT, MPMY_SUM, req); /* Don */
	MPMY_ICombine(com, com, NDIM, MPMY_DOUBLE, MPMY_SUM, req);
	MPMY_ICombine(comv, comv, NDIM, MPMY_DOUBLE, MPMY_SUM, req);
	MPMY_ICombine(&max_vsound, &max_vsound, 1, MPMY_FLOAT, MPMY_MAX, req);
	MPMY_ICombine(&max_rho, &max_rho, 1, MPMY_FLOAT, MPMY_MAX, req);
	MPMY_ICombine(&min_dt, &min_dt, 1, MPMY_FLOAT, MPMY_MIN, req);
	MPMY_ICombine(&max_dacc, &max_dacc, 1, MPMY_FLOAT, MPMY_MAX, req);
	MPMY_ICombine(&rms_dacc, &rms_dacc, 1, MPMY_FLOAT, MPMY_SUM, req);
	MPMY_ICombine(&max_nbrs, &max_nbrs, 1, MPMY_INT, MPMY_MAX, req);
	MPMY_ICombine(&min_nbrs, &min_nbrs, 1, MPMY_INT, MPMY_MIN, req);
	MPMY_ICombine(&atmax, &atmax, 1, MPMY_INT, MPMY_SUM, req);
	MPMY_ICombine(&atmin, &atmin, 1, MPMY_INT, MPMY_SUM, req);
	MPMY_ICombine(&press_max, &press_max, 1, MPMY_FLOAT, MPMY_MAX, req); /* SUM->MAX - Don */
	MPMY_ICombine(&density_max, &density_max, 1, MPMY_FLOAT, MPMY_MAX, req);
	MPMY_ICombine(&kelvin_max, &kelvin_max, 1, MPMY_FLOAT, MPMY_MAX, req);
	MPMY_ICombine_Wait(req);

/* OJO: For center of mass info use this: */
/*
*/
	VS(com, /= mtot);
	VS(comv, /= mtot);

/* Remove center of mass motion for relaxation: (Karen) */
	if (remove_com_motion) {
	    for (p = btab; p < btab+nobj; p++) {
		p->pos[0] += xx0 - com[0];
		p->pos[1] += yy0 - com[1];
		p->pos[2] += zz0 - com[2];
		p->vel[0] += vx0 - comv[0];
		p->vel[1] += vy0 - comv[1];
		p->vel[2] += vz0 - comv[2];
            }
            get_metric(btab, nobj);
            for (p = btab; p < btab+nobj; p++) {
                vx = p->vel[0];
                vy = p->vel[1];
                vz = p->vel[2];
                vt = 1.;
                v2 = vt*vt*p->gtt + vx*vx*p->gxx
                     + vy*vy*p->gyy + vz*vz*p->gzz
                     + 2.0
                     *(vt*vx*p->gxt + vt*vy*p->gyt
                     + vt*vz*p->gzt + vx*vy*p->gxy
                     + vx*vz*p->gxz + vy*vz*p->gyz);
                uut = sqrt(-1.0/v2);
                ut = vt*p->gtt + vx*p->gxt + vy*p->gyt + vz*p->gzt;
                ux = vt*p->gxt + vx*p->gxx + vy*p->gxy + vz*p->gxz;
                uy = vt*p->gyt + vx*p->gxy + vy*p->gyy + vz*p->gyz;
                uz = vt*p->gzt + vx*p->gxz + vy*p->gyz + vz*p->gzz;
                ut *= uut;
                ux *= uut;
                uy *= uut;
                uz *= uut;
                p->mom[0] = p->enth*ux;
                p->mom[1] = p->enth*uy;
                p->mom[2] = p->enth*uz;
                p->mom[3] = p->enth*ut;
                p->alfa = sqrtf_fast(-1.0/p->gutt);
                gamaold = p->gama;
                p->gama = uut * p->alfa;
                p->rho *= p->gama/gamaold;
	    }
	}

	rms_dacc = 1.0/sqrt(rms_dacc/gnobj);
	max_dacc = 1.0/max_dacc;
	if (do_point_mass) {
	    p = &point_mass;
	    ke += (float)0.5 * p->mass * Dot(p->vel, p->vel);
	    pe += (float)0.5 * p->mass * p->phi;
	}

	if (adaptive_dt) {
	    float min = min_dt;
	    if (rms_dacc < min_dt && !first_step)
	      min = rms_dacc;
	    if (dt > min && dt > dt_lower_limit) {
		dt *= (float)(3./4.);
		Msgf(("Adjusting dt down by factor of 3/4\n"));
	    }
	    while (dt > min) {
		singlPrintf("dt changing too fast\n");
		dt *= (float)(3./4.);
	    }
	    if (min > 1.6*dt) {
		dt *= (float)(4./3.);
		Msgf(("Adjusting dt up by factor of 4/3\n"));
	    }
            if (dt > dt_upper_limit) dt = dt_upper_limit; /* Don */
	}

   singlPrintf("\n++++++++++++++++++++++ QUEPASA +++++++++++++++++++++++\n");
	singlPrintf("GENRL: tpos = %g, iter = %d, size = %g, rr = %g\n", 
		    tpos, iter, sysradius, sqrt(Dot(com, com)));
	singlPrintf("       dt = %g, min dt = %g, max dacc = %g, rms dacc = %g\n", 
		    dt_last, min_dt, max_dacc, rms_dacc);
	singlPrintf("       min h/vsound = %g, max nbrs = %4d, min nbrs = %4d\n",
	       1.0/max_vsound, max_nbrs, min_nbrs);
	singlPrintf("       %d at max h, %d at min h\n", atmax, atmin);
	singlPrintf("ENER:  kin = %g, pot = %g, int = %g, tot = %g\n",
	       ke, pe, te, ke+pe+te);
	singlPrintf("MAX:   rho = %g, den = %g, temp = %g, press = %g\n",
	       max_rho,density_max, kelvin_max, press_max);
	singlPrintf("CofM:  pos = %g %g %g vel = %g %g %g\n",
		    com[0], com[1], com[2], comv[0], comv[1], comv[2]);
   singlPrintf("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n\n");

	if (MPMY_Procnum() == 0) {
#if 0
	 /* Pablo */
	    fprintf(enefile, "%12.4e %12.4e %12.4e %12.4e %12.4e\n",
		    tpos, ke, pe, te, ke+pe+te);
	    fprintf(maxfile, "%12.4e %12.4e %12.4e %12.4e %12.4e\n",
		    tpos, max_rho,density_max, kelvin_max, press_max); 
	    fprintf(cofmfile, "%12.4e %12.4e %12.4e %12.4e %12.4e %12.4e\n",
		    com[0], com[1], com[2], comv[0], comv[1], comv[2]);
	    fflush(enefile);
	    fflush(maxfile);
	    fflush(cofmfile);
         /* End Pablo */
#endif
	    fprintf(statsfile, "%15e %15e %15e %15e %15e %15e %15e %15e %15e %15e %15e %15e %15e\n", tpos, max_rho, density_max, kelvin_max, press_max, com[0], com[1], com[2], comv[0], comv[1], comv[2], hp, hx ); /* Added grav_rad output - Don */
	    fflush(statsfile);
	}
	first_step = 0;
	StopTimer(&StepTot);
	StopTimer(&StepTotWC);

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
	singlFflush();

	first_step = 0;
	if( Msg_test("memleak") ){
	    Msg_do("Memory map after iteration %d\n", iter);
	    malloc_print();
	}
    }
    if (MPMY_Procnum() == 0) Fclose(statsfile);
    Msg_flush();
    singlPrintf("Bye!\n");
    exit(0);
}

static SDF *startup(int argc, char **argv){
    SDF *csdfp;
    char msg_turn_on[512];
    char msgdir[256];
    char tmp[256];
    char *msgbase, *lastslash;
    char cfile[256];

    if (argc > 1)
 	strncpy(cfile, argv[1], sizeof(cfile));
    else
 	Getsparam("control file", cfile);
    if ((csdfp = SDFopen(NULL, cfile)) == NULL) {
 	Error("Sorry, couldn't SDFopen %s\n%s\n",
	      cfile, SDFerrstring);
    }
    /* Get the msgdir either from:
       argv[2] 
       csdfp
       argv[0]
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
    EnableTimer(&StepTotWC, "Step Wall");
    EnableTimer(&BuildTot, "Build Total");
    EnableTimer(&FindForcesTm, "Force Eval");
    EnableTimer(&ForceSPH, "Force (SPH)");
    EnableTimer(&RhoSPH, "Rho (SPH)");
    EnableTimer(&SDFreadTm, "SDFread"); SDFreadTm.enabled = 0;
    EnableTimer(&SDFwriteTm, "SDFwrite"); SDFwriteTm.enabled = 0;
    EnableCounter(&CCInt, "Cell-cell");
    EnableCounter(&BCInt, "Body-cell");
    EnableCounter(&CBInt, "Cell-body");
    EnableCounter(&BBInt, "Body-body");
    EnableCounter(&NtermsCnt, "Nterms");
    EnableCounter(&NbodyCnt, "Nbody");
    EnableCounter(&CCIntRej, "MAC fail");
    EnableCounter(&TranslateCnt, "Translate");
    EnableCounter(&SPHCnt, "SPH inter");
    EnableCounter(&SPHrej, "SPH rejects");
    EnableCounter(&nbrMACCnt, "nbr MAC");
    EnableCounter(&DeferCnt, "Deferred"); /* in newwalk.c */
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
	if ((ReadCounter(&NtermsCnt) > 0)
	    && (ReadTimer(&GravTm) > 0.0))
	    decompp->loadbal_target = ReadTimer(&GravTm)
		/ ReadCounter(&NtermsCnt);
	break;

    }
#endif

    return;
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

static float
Znow(float time)
{
    /* This assumes h = 0.5 */
    return(pow(1.5*1.023*.05*time, -2./3.)-1);
}


/* This erases the velocities, and sets them */
/* according to the Zel'dovich approximation */

static void
set_vels(body *p, int n, float real_time)
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

static int
cmpfloat(const void *a, const void *b)
{
    
    if (*(float *)a > *(float *)b)
      return 1;
    else 
      return -1;
}

/* From numerical recipes */
/* heap[0] returns the mth largest element */
static void
hpsel(unsigned long m, unsigned long n, float arr[], float heap[])
{
    unsigned long i,j,k;
    float swap;

    if (m > n/2 || m < 1) Error("probable misuse of hpsel\n");
    for (i=0;i<m;i++) heap[i]=arr[i];
    qsort(heap, m, sizeof(float), cmpfloat);
    for (i=m;i<n;i++) {
	if (arr[i] > heap[1]) {
	    heap[0]=arr[i];
	    for (j=0;;) {
		k=j << 1;
		if (k > m) break;
		if (k != m && heap[k] > heap[k+1]) k++;
		if (heap[j] <= heap[k]) break;
		swap=heap[k];
		heap[k]=heap[j];
		heap[j]=swap;
		j=k;
	    }
	}
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
    SDFgetfloatOrDie(sdfp, "pt_h", &(point_mass->h));
}


static void
point_mass_offset(body *btab, int nobj, SDF *csdfp)
{
    body *p1;
    float off[NDIM], offv[NDIM];

    SDFgetfloatOrDie(csdfp, "offset_x", &off[0]);
    SDFgetfloatOrDie(csdfp, "offset_y", &off[1]);
    SDFgetfloatOrDie(csdfp, "offset_z", &off[2]);
    SDFgetfloatOrDie(csdfp, "offset_vx", &offv[0]);
    SDFgetfloatOrDie(csdfp, "offset_vy", &offv[1]);
    SDFgetfloatOrDie(csdfp, "offset_vz", &offv[2]);

    for (p1 = btab; p1 < btab+nobj; p1++) {
	VV(p1->pos, += off);
	VV(p1->vel, += offv);
    }
    singlPrintf("Did star offset +/- (%g,%g,%g), (%g,%g,%g)\n",
		off[0], off[1], off[2], offv[0], offv[1], offv[2]);
}

