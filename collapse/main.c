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
#include "physics.h"
#include "physics_sph.h"
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
#include "integrate.h"

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

extern void sigio_setup(void);

static void SanityCheck(body *btab, int nobj, int gnobj, double *mtotp);
static void SPHSanityCheck(SPHbody *btab, int nobj, int gnobj, double *mtotp);
static void set_vels(body *p, int n, float real_time);
static void Output(body *btab, int nobj, const char *outname, int iter);
static void SPHOutput(SPHbody *btab, int nobj, const char *outname, int iter);
static void WindOutput(SPHbody *btab, int nobj, windbody *windbtab, 
		       int windnobj, const char *outname, int iter);
static SDF *startup(int argc, char **argv);
static void Periodic(tree_t *tp, float size);
static void PeriodicSPH(tree_t *tp, float size, float vsize);
static void WrapPeriodic(body *bp, int n, float *rmin, float *rmax, float sz, int cosmology, int log_time, float tpos, float dt);
static void SPHWrapPeriodic(SPHbody *bp, int n, float *rmin, float *rmax, float sz, int cosmology, int log_time, float tpos, float dt);
static void FixCube(body *b, int nobj, float l, float gm);
static void FixGlobalForce(body *bp, int n);
static void Fix_h(SPHbody *btab, int nobj, int nbrcut_max, int nbrcut_min, float nbrcut_fac, float max_h, float min_h);
static void Diags(body *btab, int nobj, double ke, double pe, double *etot, float dt_last, int iter, int gnobj);
static void SPHDiags(SPHbody *btab, int nobj, double ke, double pe, double te, double *etot, float dt_last, int iter, int gnobj, float *tmin, int *tbad);
static void Fix_dt(float *dt, float dt_max, int tlow_cut, float tmin, int tbad, int dtshort, int dtlong, int limit_high, int limit_low);
static void ReadCosmo(SDF *sdfp, struct cosmo_s *cosmo, float tpos, float *R0p);
static void CosmoPush(struct cosmo_s *p, float time);
static int dark_need_update(float dark_tacc, float dark_dt);
/*  static float IdtSPHGetCost(const SPHbody *ptr); */

/* In shrink.c */
/*  void ShrinkBtab(SPHbody **SPHbtabp, body *btabp, int *nobj, float r_limit); */
/* void ShrinkBtab2(SPHbody **SPHbtabp, int *nobj, float r_limit); */
void AdjustBtab(SPHbody **SPHbtabp, int *nobj, int gnobj, windbody *windbtab, 
		int windnobj, float r_limit, float dt, int iter, float tpos,
		int *added_particles);
void AdjustBtab2(SPHbody **SPHbtabp, int *nobj, int gnobj, windbody *windbtab, 
		int windnobj, float r_limit, float dt, int iter, float tpos,
		int *added_particles, float *newmass);

static int maxmem(void);
static int maxheap(void);

float Znow(float time);
float Hnow(float time);

Timer_t StepTot, StepTotWC, BuildTot;
Timer_t FindForcesTm;
Timer_t RhoSPH, ForceSPH, PerTmSPH;
Timer_t FixCubeTm;
Counter_t NbodyCnt;
Counter_t MemCnt;
Counter_t HeapCnt_;	/* HeapCnt is in the SunOS name space?! */
Counter_t NtermsCnt;
Counter_t SPHbodyCnt;

Timer_t WTermTm,WNTTm,PerTm;

static float dt;
static float dt_max;
static float sysradius;
static float tpos, tvel;
static int cosmology;
static float R0;
static float this_tol, this_eps;
static float frac_tol;
static float Gamma;		/* lowercase is a math.h function */
static int default_nterms;
static float centmass;

static double gnterms;
static float courant_number;
static int adaptive_dt;
static int independent_dt;
static int dark_independent_dt;

int do_diffusion;  /* used in main and in sph.c */

#ifdef __PARAGON__
void
chk_slow(int die)
{
    int i, k;
    Timer_t t;
    double tt;

    int *array = Calloc(200000, sizeof(int));
    EnableTimer(&t, "Slow");
    MPMY_Sync();
    StartTimer(&t);
    for (i = 0; i < 10; i++) {
	for (k = 0; k < 200000; k++) {
	    array[k] += i;
	}
    }
    StopTimer(&t);
    Free(array);
    tt = ReadTimer(&t);
    DisableTimer(&t);
    if (tt > 0.40) Shout("Node %d is slow (%.3f)\n", _myphysnode(), tt);
    MPMY_Combine(&tt, &tt, 1, MPMY_DOUBLE, MPMY_MAX);
    if ((tt > 0.40) && die) exit(1);
}
#endif

int
main(int argc, char *argv[])
{
    int gnobj, nobj;
    int SPHgnobj, SPHnobj, SPHoldnobj;
    int windgnobj, windnobj;
    float r_inner;
    int PMgnobj, PMnobj;
    int SPHsinkgnobj, SPHsinknobj;
    body *btab, *p;
    body *pmtab;
    SPHbody *SPHbtab, *SPHsinkbtab = NULL, *q;
    windbody *windbtab;
    float eps;			/* Plummer smoothing length */
    float tol;			/* MAC tolerance */
		/* for big MAC, this is multiplied by M/(rsize*rsize) */
    int i;
    float rmin[NDIM], rmax[NDIM];
    int nsteps;
    int first_step = 1;
    int added_particles = 0;
    int stride = sizeof(body)/sizeof(float);
    int SPHstride = sizeof(SPHbody)/sizeof(float);
    int SPHstride2 = sizeof(SPHbody)/sizeof(double);
    int SPHstride3 = sizeof(SPHbody)/sizeof(unsigned int);
    int do_output;
    int output_freq;
    int timer_freq;
    float sort_tol;
    int iter;
    float CWfac;
    int ntimer_detail;
    int log_time = 0;		/* if true, use dt \propto t */
    int comov_eps = 0;		/* if true, use comoving epsilon*/
    float comov_eps_epoch;
    int setpvel = 0;
    char outnamebase[256];
    SDF *csdfp;			/* SDF pointer to control file */
    SDF *sdfp = NULL;
    float tposlast;
    int save_first;		/* save first step (for acc testing) */
    double pe, ke, te;
    double dark_ke, dark_pe;
    double etot;
    double mtot, SPHmtot;
    sortresult_t sortedbtab, SPHsortedbtab, sortedatab;
    tree_t thetree, SPHtree, SPHsinktree, *sinkptr = NULL;
    char name[256];
    int do_BH, do_DL, do_Bmax, do_Arel;
    int MACtype = BMAX_MAC;
    int image_freq, x_pixels, y_pixels, log_image;
    int do_periodic;
    inherit_t inherit;
    macv_t mac;
    int timeout;
    int set_id;
    float dt_last;
    float new_h, new_u;
    int do_sph, do_grav, do_winds;
    int do_point_mass, do_point_mass2;
    float newmass = 0.0, totnewmass = 0.0;
    int exact_rho;
    float visc_alpha, visc_beta, visc_epsilon, heat_f1;
    int nbrcut_max, nbrcut_min;
    float nbrcut_fac;
    float min_h, max_h;
    int udot_limit[2];
    float vsz;
    float tmin;
    int tbad;
    int tlow_cut, dt_short, dt_long;
    accbody *SPHatab, *pa;
    int SPHanobj;
    void *decomp_info = NULL;
    int do_restart;
    float dark_tacc = -1e30;	/* initialize so dark_need_update is true */
    float dark_dt;
    int did_dark_update;
    int SPHnupdate;
    int make_sink_tree;
    int has_grav_data;

    MPMY_Init(&argc, &argv);
    singlPrintf("Welcome to the variable O() integrator running on %d procs\n",
		MPMY_Nproc());
    csdfp = startup(argc, argv);

    SetBoundary(30); /* From integrate.c; seems not to be working (?) */

    SDFgetintOrDefault(csdfp, "timeout", &timeout, 600);
    if (timeout > 0) MPMY_TimeoutSet(timeout);
#ifdef __PARAGON__
    {
	int fail_if_slow;
	SDFgetintOrDefault(csdfp, "fail_if_slow", &fail_if_slow, 0);
	chk_slow(fail_if_slow);
    }
#endif
    SDFgetstring(csdfp, "datafile", name, sizeof(name));
    SDFgetintOrDefault(csdfp, "do_restart", &do_restart, 0);
    SDFgetintOrDefault(csdfp, "do_periodic", &do_periodic, 0);
    SDFgetintOrDefault(csdfp, "cosmology", &cosmology, 0);
    SDFgetintOrDefault(csdfp, "set_id", &set_id, 0);
    SDFgetintOrDefault(csdfp, "setpvel", &setpvel, 0);
    SDFgetintOrDefault(csdfp, "do_sph", &do_sph, 0);
    SDFgetintOrDefault(csdfp, "do_diffusion", &do_diffusion, 0);
    SDFgetintOrDefault(csdfp, "do_grav", &do_grav, 1);
    SDFgetintOrDefault(csdfp, "do_winds", &do_winds, 0);
    SDFgetintOrDefault(csdfp, "do_point_mass", &do_point_mass, 0);
    SDFgetintOrDefault(csdfp, "do_point_mass2", &do_point_mass2, 0);
    SDFgetintOrDefault(csdfp, "has_grav_data", &has_grav_data, do_grav);
    if (do_sph || do_grav) {
	if (!((strncmp(name, "test", 4) == 0))) {
	    if (SDFhasname("SPHdatafile", csdfp) || do_restart) {
		char iname[256];

		if (do_sph) {
		    SDFgetfloatOrDefault(csdfp, "new_h", &new_h, 0.0);
		    SDFgetfloatOrDefault(csdfp, "new_u", &new_u, 0.0);
		    if (do_restart) sprintf(iname, "%s_sph.restart", name);
		    else SDFgetstring(csdfp, "SPHdatafile", iname, 
				      sizeof(iname));
		    sdfp = SPHRead(iname, csdfp, &SPHbtab, &SPHgnobj, &SPHnobj,
				   set_id, setpvel, new_h, new_u);
		} else SPHgnobj = SPHnobj = 0;

		if (has_grav_data) {
		    if (do_restart) sprintf(iname, "%s.restart", name);
		    else SDFgetstring(csdfp, "datafile", iname, sizeof(iname));
		    sdfp = DarkRead(iname, csdfp, (void **)&btab, &gnobj, 
				    &nobj, set_id, setpvel);
		} else {
		    gnobj = nobj = 0;
		    btab = Malloc(sizeof(body)); /* realloced later */
		}

		if (do_point_mass) {
		    SDFgetfloatOrDefault(csdfp, "r_inner", &r_inner, 0.05);
		    if (do_restart) sprintf(iname, "%s.restart", name);
		    else SDFgetstring(csdfp, "datafile", iname, sizeof(iname));
		    sdfp = DarkRead(iname, csdfp, (void **)&pmtab, &PMgnobj, 
				    &PMnobj, set_id, setpvel);
		    Msgf(("lx = %e; ly = %e; lz = %e; accmass = %e\n", 
			  pmtab->l[0], pmtab->l[1], pmtab->l[2], 
			  pmtab->accmass));
		} else if (do_point_mass2) {
		    SDFgetfloatOrDefault(csdfp, "r_inner", &r_inner, 0.05);
		    SDFgetfloatOrDie(sdfp, "centmass", &centmass);
		    PMgnobj = PMnobj = 0;
		    pmtab = Malloc(sizeof(body)); /* realloced later */
/*  		    SDFgetstring(csdfp, "windfile", iname, sizeof(iname)); */
/*   		    sdfp = SPHRead(iname, csdfp, &SPHwind, &windnobj, */
/* 				   &windnobj, set_id, setpvel, new_h,new_u); */
		} else {
		    PMgnobj = PMnobj = 0;
		    pmtab = Malloc(sizeof(body)); /* realloced later */
		}

		if (do_winds) {
		    sdfp = WindRead(iname, csdfp, &windbtab, &windgnobj, 
				    &windnobj);
		} else windgnobj = windnobj = 0;

		SDFgetfloatOrDefault(sdfp, "dt", &dt, 0.0);
		SDFgetfloatOrDefault(sdfp, "dark_dt", &dark_dt, dt);
	    } else {
		sdfp = InitRead(name, csdfp, (void **)&btab, &gnobj, &nobj, 
				&SPHbtab, &SPHgnobj, &SPHnobj, 
				set_id, setpvel, new_h, new_u);
	    }
	    FixNterms(btab, nobj);
	    SPHFixNterms(SPHbtab, SPHnobj);
	    SDFgetfloatOrDefault(sdfp, "Gnewt", &cosmo.GNewt, (float)1.0);
	    SDFgetfloatOrDefault(sdfp, "tpos",  &tpos, (float)0.0);
	    tvel = tpos;
	    SDFgetintOrDefault  (sdfp, "iter",  &iter, 0);
	    if (cosmology) ReadCosmo(sdfp, &cosmo, tpos, &R0);
	    if(sdfp) SDFclose(sdfp);
	} else {
	    SPHTestData(csdfp, &SPHbtab, &SPHgnobj, &SPHnobj, do_periodic);
	    nobj = gnobj = 0;
	    cosmo.GNewt = (float)1.0;
	    tvel = tpos = (float)0.0;
	    iter = 0;
	    if (do_periodic) R0 = 1.0;
	}
    }
    singlPrintf("Maxmem after data read is %d (%d)\n", maxmem(), maxheap());
    if( Msg_test("memleak") ){
	Msg_do("Memory map after data read\n");
	malloc_print();
    }

    SDFgetfloatOrDefault(csdfp, "epsilon", &eps, 0.0);
    if (do_grav) {
	SDFgetintOrDefault(csdfp, "do_DL", &do_DL, 0);
	SDFgetintOrDefault(csdfp, "do_BH", &do_BH, 0);
	SDFgetintOrDefault(csdfp, "do_Bmax", &do_Bmax, 0);
	SDFgetintOrDefault(csdfp, "do_Arel", &do_Arel, 0);
	if (do_BH || do_Bmax) 
	  SDFgetfloatOrDie(csdfp, "theta", &tol);
	else
	  SDFgetfloatOrDie(csdfp, "errtol", &tol);
	SDFgetfloatOrDefault(csdfp, "frac_tol", &frac_tol, 0.0);
    }
    SDFgetfloatOrDefault(csdfp, "CWfac", &CWfac, 0.0);
    if (!do_restart) {
	SDFgetfloatOrDie(csdfp, "dt", &dt);
/*  SDFgetfloatOrDefault(csdfp, "dark_dt", &dark_dt, do_grav ? dt : 1e30); */
    }
    SDFgetfloatOrDefault(csdfp, "dark_dt", &dark_dt, do_grav ? dt : 1e30);
    SDFgetintOrDie(csdfp, "nsteps", &nsteps);
    SDFgetintOrDefault(csdfp, "log_time", &log_time, 0);
    SDFgetintOrDefault(csdfp, "comov_eps", &comov_eps, 0);
    SDFgetfloatOrDefault(csdfp, "comov_eps_epoch", &comov_eps_epoch, 10.0);
    SDFgetintOrDefault(csdfp, "save_first", &save_first, 0);
    SDFgetintOrDefault(csdfp, "ntimer_detail", &ntimer_detail, 0);
    SDFgetintOrDefault(csdfp, "exact_rho", &exact_rho, 0);
    SDFgetfloatOrDefault(csdfp, "visc_alpha", &visc_alpha, (float)1.0);
    SDFgetfloatOrDefault(csdfp, "visc_beta", &visc_beta, (float)2.0);
    SDFgetfloatOrDefault(csdfp, "visc_epsilon", &visc_epsilon, (float)1e-2);
    SDFgetfloatOrDefault(csdfp, "heat_f1", &heat_f1, (float)1.0);
    SDFgetfloatOrDefault(csdfp, "gamma", &Gamma, (float)(5.0/3.0));
    SDFgetfloatOrDefault(csdfp, "courant_number", &courant_number, (float)0.4);
    SDFgetfloatOrDefault(csdfp, "min_h", &min_h, (float)0.0);
    SDFgetfloatOrDefault(csdfp, "max_h", &max_h, (float)1e30);
    SDFgetintOrDefault(csdfp, "nbrcut_max", &nbrcut_max, 500);
    SDFgetintOrDefault(csdfp, "nbrcut_min", &nbrcut_min, 10);
    SDFgetfloatOrDefault(csdfp, "nbrcut_fac", &nbrcut_fac, (float)0.1);
    SDFgetintOrDefault(csdfp, "adaptive_dt", &adaptive_dt, 1);
    SDFgetintOrDefault(csdfp, "independent_dt", &independent_dt, 0);
    SDFgetintOrDefault(csdfp, "dark_independent_dt", &dark_independent_dt, 0);
    SDFgetintOrDefault(csdfp, "default_nterms", &default_nterms, 100);
    if (adaptive_dt) {
	SDFgetintOrDefault(csdfp, "tlow_cut", &tlow_cut, 40);
	SDFgetintOrDefault(csdfp, "dt_short", &dt_short, 0);
	SDFgetintOrDefault(csdfp, "dt_long", &dt_long, 10);
	SDFgetfloatOrDefault(csdfp, "dt_max", &dt_max, 1e30);
    }
    if (do_Bmax) MACtype = BMAX_MAC;
    else if (do_BH) MACtype = BH_MAC;
    else if (do_Arel) MACtype = AREL_MAC;
    else if (do_grav) Error("No MAC specified\n");

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
    SDFgetfloatOrDefault(csdfp, "sort_tol", &sort_tol, 0.01);
    SDFgetintOrDefault(csdfp, "image_freq", &image_freq, 0);
    SDFgetintOrDefault(csdfp, "x_pixels", &x_pixels, 512);
    SDFgetintOrDefault(csdfp, "y_pixels", &y_pixels, 512);
    SDFgetintOrDefault(csdfp, "log_image", &log_image, 0);

    if(csdfp) 
	SDFclose(csdfp);

    if (do_periodic) {
	EnableTimer(&FixCubeTm, "Fix Cube");
    }

    singlPrintf("float errtol = %g;\n", tol);
    singlPrintf("float dt = %g;\n", dt);
    singlPrintf("float dark_dt = %g;\n", dark_dt);
    singlPrintf("float epsilon = %g;\n", eps);
    singlPrintf("int iter = %d;\n", iter);
    singlPrintf("int nsteps = %d;\n", nsteps);
    singlPrintf("int nproc = %d;\n", MPMY_Nproc());
    singlPrintf("int do_Bmax = %d;\n", do_Bmax);
    singlPrintf("int do_BH = %d;\n", do_BH);
    singlPrintf("int do_Arel = %d;\n", do_Arel);
    singlPrintf("int do_DL = %d;\n", do_DL);
    singlPrintf("int exact_rho = %d;\n", exact_rho);
    singlPrintf("float courant_number = %g;\n", courant_number);
    singlPrintf("float gamma = %f;\n", Gamma);
    singlPrintf("float visc_alpha = %g;\n", visc_alpha);
    singlPrintf("float visc_beta = %g;\n", visc_beta);
    singlPrintf("float visc_epsilon = %g;\n", visc_epsilon);
    singlPrintf("float heat_f1 = %g;\n", heat_f1);
    singlPrintf("float min_h = %g;\n", min_h);
    singlPrintf("float max_h = %g;\n", max_h);
    singlPrintf("int adaptive_dt = %d;\n", adaptive_dt);
    singlPrintf("int independent_dt = %d;\n", independent_dt);
    singlPrintf("int dark_independent_dt = %d;\n", dark_independent_dt);
    if (adaptive_dt) {
	singlPrintf("int tlow_cut = %d;\n", tlow_cut);
	singlPrintf("int dt_long = %d;\n", dt_long);
	singlPrintf("int dt_short = %d;\n", dt_short);
	singlPrintf("float dt_max = %g;\n", dt_max);
    }
    if (do_point_mass || do_point_mass2) {
        singlPrintf("float r_inner = %f;\n", r_inner);
	singlPrintf("float GNewt = %e;\n", cosmo.GNewt);
    }
    singlPrintf("int do_diffusion = %d;\n", do_diffusion);
    if( do_output ){
	singlPrintf("Output to %s.nnnn, every %d steps\n", 
	       outnamebase, output_freq);
    }else{
	singlPrintf("No output.\n");
    }
    singlPrintf("int timer_freq = %d;\n", timer_freq);
    singlPrintf("float sort_tol = %.4f;\n", sort_tol);
    singlPrintf("int do_periodic = %d;\n", do_periodic);
    if (log_time) Error("This code does not support log_time\n");
    if (cosmology) {
	singlPrintf("int cosmology = %d;\n", cosmology);
	singlPrintf("int comov_eps = %d;\n", comov_eps);
	singlPrintf("float comov_eps_epoch = %f;\n", comov_eps_epoch);
	singlPrintf("int setpvel = %d;\n", setpvel);
	singlPrintf("float R0 = %f;\n", R0);
    }

    singlFflush();
    if (do_sph) SPHSanityCheck(SPHbtab, SPHnobj, SPHgnobj, &SPHmtot);

    SetupTree(&thetree, NDIM,
	      sizeof(body), sizeof(cell), TBODYSZ, sizeof(cofmdata),
	      (pq_keyproto)GetKeyFromStruct, (pq_wgtproto)GetCost,
	      CofmFromDaugh, (cellfromcofm_t)CellFromCofm);

    SetupTree(&SPHtree, NDIM,
	      sizeof(SPHbody), sizeof(SPHcell), SPHTBODYSZ, sizeof(SPHcofmdata),
	      (pq_keyproto)SPHGetKeyFromStruct, (pq_wgtproto)SPHGetCost,
	      SPHCofmFromDaugh, (cellfromcofm_t)SPHCellFromCofm);

    SetupTree(&SPHsinktree, NDIM,
	      sizeof(SPHbody), sizeof(SPHcell), SPHTBODYSZ, sizeof(SPHcofmdata),
	      (pq_keyproto)SPHGetKeyFromStruct, (pq_wgtproto)SPHGetCost,
	      SPHCofmFromDaugh, (cellfromcofm_t)SPHCellFromCofm);


    for (q = SPHbtab; q < SPHbtab+SPHnobj; q++) {
	q->dt = q->dt_next = dt;
	q->tacc = -1e30;
    }

    dt_last = dt;
    SPH_setup(NDIM);
    inherit = (inherit_t)InheritSinkNlogN;

    if (do_DL)
      mac = (macv_t)DLRcritMAC;
    else
      mac = (macv_t)RcritMAC;

    this_eps = eps;
    this_tol = tol;

    /* Testing initialization */
    for (q = SPHbtab; q < SPHbtab+nobj; q++) {
	VS(q->acc, = 0.0);
	VS(q->acc_last, = 0.0);
	VS(q->grav_acc, = 0.0);
	q->phi = 0.0;
    }

    for (nsteps += iter; iter <= nsteps; iter++) {
	if (timeout > 0) MPMY_TimeoutReset(timeout);
	/* Reset timers and counters */
	ClearEnabledTimers();
	ClearEnabledCounters();
	StartTimer(&StepTotWC);
	StartTimer(&StepTot);

	if (do_point_mass || do_point_mass2) {
	  SPHoldnobj = SPHnobj;
/*  	  ShrinkBtab((SPHbody **)&SPHbtab, pmtab, &SPHnobj, r_inner); */
/*   	  ShrinkBtab2((SPHbody **)&SPHbtab, &SPHnobj, r_inner);  */
/*  	  AdjustBtab((SPHbody **)&SPHbtab, &SPHnobj, SPHgnobj, windbtab,  */
/* 		     windnobj, r_inner, dt_last, iter, tpos,  */
/* 		     &added_particles); */
 	  AdjustBtab2((SPHbody **)&SPHbtab, &SPHnobj, SPHgnobj, windbtab,
		      windnobj, r_inner, dt_last, iter, tpos,
		      &added_particles, &newmass);

	  MPMY_Combine(&SPHnobj, &SPHgnobj, 1, MPMY_FLOAT, MPMY_SUM);
	  MPMY_Combine(&newmass, &totnewmass, 1, MPMY_FLOAT, MPMY_SUM);
	  centmass += totnewmass;
	  Msgf(("Iter: %d: Added %d bodies to SPHbtab\nBH mass = %f\n", iter, 
		SPHnobj-SPHoldnobj, centmass));
/* 	  SPHFixId(SPHbtab, SPHnobj, SPHgnobj); */
	}

	/* comoving smoothing */
	/* Note: behavior changed Jan. 25, 1996. Beware of old ctl files */
	if (comov_eps && (Znow(tpos)+1.0 >= comov_eps_epoch)) 
	  this_eps = eps*comov_eps_epoch/(Znow(tpos)+(float)1.0);
	else this_eps = eps;

	/* Add sph particles to btab for gravity */
	if (do_grav) {
	    GravPlusSPH((void **)&btab, &nobj, SPHbtab, SPHnobj);
	    SanityCheck(btab, nobj, gnobj+SPHgnobj, &mtot); /* need mtot */

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
	} else {
	    SPHFindBbox(SPHbtab, SPHnobj, rmin, rmax);
	    sysradius = 0.5*FixRsize(rmin, rmax);
	}
	Msgf(("rmin=(%g, %g, %g) rmax=(%g, %g, %g)\n",
	      rmin[0], rmin[1], rmin[2], 
	      rmax[0], rmax[1], rmax[2]));

	if (do_grav && dark_need_update(dark_tacc, dark_dt)) {
	    /* We aren't using the first two params */
	    SetTol(0, 0, cosmo.GNewt, this_eps, gnobj+SPHgnobj);
	    FixKeys(btab, nobj, GETKEY);
	    
	    if (MACtype == AREL_MAC) this_tol = tol*mtot/(sysradius*sysradius);
	    SetupCofm(MACtype, this_tol, frac_tol);
	    singlPrintf("BuildTree, tol=%g, frac_tol=%g\n", this_tol, frac_tol);
	    
	    StartTimer(&BuildTot);
	    pqsortsetup(&sortedbtab, btab, nobj, sizeof(body), sort_tol, Realloc_f);
	    BuildTree(&thetree, &sortedbtab);
	    /* PrintTree(&thetree, PrintBodyContents, PrintCellContents); */
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
	    
	    StartTimer(&FindForcesTm);
	    WalkInit(&thetree, &thetree, sizeof(Sink), mac, inherit);
	    StartTimer(&PerTm);
	    if (do_periodic) {
		singlPrintf("FindForces (periodic), this_eps=%g\n", this_eps);
		Periodic(&thetree, 2.0*sysradius);
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
	    WalkTerminate();
	    StopTimer(&WTermTm);
	    if (cosmology) FixGlobalForce(btab, nobj);
	    if (do_periodic) FixCube(btab, nobj, sysradius, cosmo.GNewt*mtot);
	    StopTimer(&FindForcesTm);
	    singlPrintf("FindForces done %d (%d)\n", maxmem(), maxheap());

	    FreeTree(&thetree);
	    singlPrintf("FreeTree done %d (%d)\n", maxmem(), maxheap());
	    Msgf(("FreeTree done\n"));
	    dark_tacc = tpos;
	    did_dark_update = 1;
	    singlPrintf("Updated %d grav accs\n", gnobj+SPHgnobj);
	} else did_dark_update = 0;

	/* Must do this before SPH */
	if (setpvel) {
	    setpvel = 0;
	    set_vels(btab, nobj, tpos);
	    singlPrintf("Velocities adjusted to linear approximation.\n");
	}
	
	/* Do image before GravMinusSPH if you want to image all particles */
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

	if (do_grav) GravMinusSPH((void **)&btab, &nobj, &SPHatab, &SPHanobj);

	/* This should be the high-water mark for memory use */
	AddCounter(&MemCnt, malloc_used()/1024);

	if (do_sph && (first_step || exact_rho)) {
	    singlPrintf("BuildTree\n");
	    StartTimer(&BuildTot);
	    pqsortsetup(&SPHsortedbtab, SPHbtab, SPHnobj, sizeof(SPHbody), sort_tol, Realloc_f);
	    SPHFixKeys(SPHbtab, SPHnobj, SPHGetKey);
	    BuildTree(&SPHtree, &SPHsortedbtab);
	    SPHbtab = SPHsortedbtab.data;
	    SPHnobj = SPHsortedbtab.nobj;
	    StopTimer(&BuildTot);
	    StartTimer(&RhoSPH);
	    SetSPH(visc_alpha, visc_beta, visc_epsilon, heat_f1,
		   Gamma, SPHgnobj, macRho, nbrMAC);
	    /* Periodic does multiple calls to Walk, so we must init here */
	    /* rather than in inherit */
	    for (q = SPHbtab; q < SPHbtab+SPHnobj; q++) {
		if (SPH_need_update(q)) {
		    q->rho = (float)0.0;
		    q->drho_dt = (float)0.0;
		    q->udot = (float)0.0;
		    q->nbrs = 0;
		    q->nterms = 0;
		    VS(q->lvel, = 0.0);
		}
	    }
	    WalkInit(&SPHtree, &SPHtree, sizeof(SinkSPH), 
		     (macv_t)SPHgate, (inherit_t)InheritSPH);
	    if (do_periodic) {
		singlPrintf("FindRho (periodic)\n");
		vsz = (cosmology) ? 2.0*sysradius*Hnow(tvel) : 0.0;
		PeriodicSPH(&SPHtree, 2.0*sysradius, vsz);
		singlPrintf("FindRho (periodic) done\n");
	    }
	    singlPrintf("FindRho\n");
	    WalkNT(&SPHtree);
	    WalkTerminate();
	    update_final(SPHbtab, SPHnobj, dt, &udot_limit[0], &udot_limit[1]);
	    StopTimer(&RhoSPH);
	    FreeTree(&SPHtree);
	    singlPrintf("FreeTree done\n");
	}

	if (do_sph) {
	    SPHFixKeys(SPHbtab, SPHnobj, SPHGetKey);
	    /* This sets rho_est and pr for communication during BuildTree */
	    update_intermediate(SPHbtab, SPHnobj, dt_last, 
				!(first_step || exact_rho), 0);

	    SPHsinknobj = 0;
	    for (q = SPHbtab; q < SPHbtab+SPHnobj; q++) {
		if (SPH_need_update(q)) SPHsinknobj++;
	    }
	    MPMY_Combine(&SPHsinknobj, &SPHsinkgnobj, 1, MPMY_INT, MPMY_SUM);
	    if (SPHsinkgnobj > SPHgnobj/4 || !SPHsinkgnobj) make_sink_tree = 0;
	    else make_sink_tree = 1;

	    if (make_sink_tree) {
		int n = 0;
		SPHsinkbtab = Malloc(SPHsinknobj * sizeof(SPHbody));
		for (q = SPHbtab; q < SPHbtab+SPHnobj; q++) {
		    if (SPH_need_update(q)) {
			memcpy(SPHsinkbtab+n++, q, sizeof(SPHbody));
		    }
		}
		singlPrintf("Build Sink Tree\n");
		StartTimer(&BuildTot);
		pqsortsetup(&SPHsortedbtab, SPHsinkbtab, SPHsinknobj, sizeof(SPHbody), sort_tol, Realloc_f);
		Msgf(("Build Sink Tree, nobj = %d\n", SPHsinknobj));
		BuildTree(&SPHsinktree, &SPHsortedbtab);
		SPHsinkbtab = SPHsortedbtab.data;
		SPHsinknobj = SPHsortedbtab.nobj;
		Msgf(("Build Sink Tree done, nobj = %d\n", SPHsinknobj));
		StopTimer(&BuildTot);
		sinkptr = &SPHsinktree;
	    }
	    singlPrintf("BuildTree\n");
	    StartTimer(&BuildTot);
	    if (did_dark_update || make_sink_tree) decomp_info = SaveDecomp();
	    pqsortsetup(&SPHsortedbtab, SPHbtab, SPHnobj, sizeof(SPHbody), sort_tol, Realloc_f);
	    Msgf(("Build Tree, nobj = %d\n", SPHnobj));
	    BuildTree(&SPHtree, &SPHsortedbtab);
	    SPHbtab = SPHsortedbtab.data;
	    SPHnobj = SPHsortedbtab.nobj;
	    Msgf(("Build Tree done, nobj = %d\n", SPHnobj));
	    StopTimer(&BuildTot);
	    AddCounter(&SPHbodyCnt, SPHnobj);

	    if (did_dark_update) {
		/* Make accs computed above congruent to SPH particle decomp */
		SetDecomp(decomp_info);
		singlPrintf("Re-arrange acctab\n");
		pqsortsetup(&sortedatab, SPHatab, SPHanobj, sizeof(accbody), sort_tol, Realloc_f);
		SPHatab = pqsort(&sortedatab, UnityCost, accbodyGetKey);
		SPHanobj = sortedatab.nobj;
		if (!make_sink_tree) ClearDecomp(decomp_info);

		if (SPHanobj != SPHnobj)
		  Error("acc table has %d entries, should be %d\n", SPHanobj, SPHnobj);
	    
		/* Update gravitational forces from entire system */
		for (i = 0; i < SPHnobj; i++) {
		    q = SPHbtab+i;
		    pa = SPHatab+i;
		    if (KeyNEQ(q->key, pa->key)) Error("Key mismatch in acc table\n");
		    VV(q->grav_acc, = pa->grav_acc);
		    q->phi = pa->phi;
		    q->grav_nterms = pa->grav_nterms;
		}
		Free(SPHatab);
	    }

	    if (!make_sink_tree) {
		SPHsinkgnobj = SPHgnobj;
		SPHsinknobj = SPHnobj;
		SPHsinkbtab = SPHbtab;
		sinkptr = &SPHtree;
	    }

	    StartTimer(&ForceSPH);
	    SetSPH(visc_alpha, visc_beta, visc_epsilon, heat_f1, Gamma, SPHgnobj, 
		   macSPH, nbrMAC);
	    for (q = SPHsinkbtab; q < SPHsinkbtab+SPHsinknobj; q++) {
		if (SPH_need_update(q)) {
		    q->rho = (float)0.0;
		    q->drho_dt = (float)0.0;
		    q->udot = (float)0.0;
		    q->nbrs = 0;
		    q->nterms = 0;
		    VS(q->acc, = 0.0);
		    VS(q->lvel, = 0.0);
		}
	    }

	    WalkInit(&SPHtree, sinkptr, sizeof(SinkSPH), 
		     (macv_t)SPHgate, (inherit_t)InheritSPH);
	    StartTimer(&PerTmSPH);
	    if (do_periodic) {
		singlPrintf("ForceSPH (periodic)\n");
		vsz = (cosmology) ? 2.0*sysradius*Hnow(tvel) : 0.0;
		PeriodicSPH(sinkptr, 2.0*sysradius, vsz);
		singlPrintf("ForceSPH (periodic) done\n");
	    }
	    StopTimer(&PerTmSPH);
	    singlPrintf("ForceSPH\n");
	    WalkNT(sinkptr);
	    WalkTerminate();
	    singlPrintf("ForceSPH done\n");
	    udot_limit[0] = udot_limit[1]  = 0;
	    update_final(SPHsinkbtab, SPHsinknobj, dt, &udot_limit[0], &udot_limit[1]);
	    StopTimer(&ForceSPH);
	    /* This should be the high-water mark for memory use */
	    AddCounter(&MemCnt, malloc_used()/1024);

	    FreeTree(&SPHtree);
	    singlPrintf("FreeTree done %d (%d)\n", maxmem(), maxheap());

	    if (make_sink_tree) {
		SPHbody *r;
		FreeTree(&SPHsinktree);
		/* Make sinkbtab congruent with SPHbtab */
		SetDecomp(decomp_info);
		pqsortsetup(&SPHsortedbtab, SPHsinkbtab, SPHsinknobj, sizeof(SPHbody), sort_tol, Realloc_f);
		singlPrintf("Re-arrange sinks\n");
		Msgf(("Re-arrange sinks, nobj = %d\n", SPHsinknobj));
		SPHsinkbtab = pqsort(&SPHsortedbtab, UnityCost, SPHGetKey);
		SPHsinknobj = SPHsortedbtab.nobj;
		Msgf(("Re-arrange sinks done, nobj = %d\n", SPHsinknobj));
		ClearDecomp(decomp_info);

		/* Replace updated bodies */
		r = SPHbtab;
		for (q = SPHsinkbtab; q < SPHsinkbtab+SPHsinknobj; q++) {
		    for (; KeyNEQ(r->key, q->key) && r < SPHbtab+SPHnobj; r++)
		      /* NULL */;
		    if (r->ident != q->ident) {
			/* Keys can match if pos's are close enough */
			for (; KeyNEQ(r->key, q->key) 
				 && r->ident != q->ident 
				 && r < SPHbtab+SPHnobj; r++)
			    /* NULL */;
			if (r == SPHbtab+SPHnobj)
			    Error("Ident mismatch in sinktab\n");
		    }
		    memcpy(r, q, sizeof(SPHbody));
		}
		Free(SPHsinkbtab);
	    }
	}

	if (do_point_mass) {
	    /* Need to add code for parallel stuff here */
	    for (p = pmtab; p < pmtab+PMgnobj; p++) {
	      VS(p->acc, = 0.0);
	      p->phi = 0.0;
	    }
	    for (p = pmtab; p < pmtab+PMgnobj; p++) {
	      update_point_mass(pmtab, PMnobj, p, eps*eps, cosmo.GNewt);
	    }
	    for (p = pmtab; p < pmtab+PMgnobj; p++) {
	      update_point_SPHmass(SPHbtab, SPHnobj, p, eps*eps, cosmo.GNewt);
	    }
	    singlPrintf("Updated %d point-mass accs\n", PMgnobj);
	}

	if (do_point_mass2) {
	    update_point_SPHmass2(SPHbtab, SPHnobj, eps*eps, cosmo.GNewt, 
				  centmass);
	}

	MPMY_Sync();

	if (ForceOutput()
	    || (do_output && !first_step
		&& ((iter+output_freq) % output_freq == 0))
	    || (save_first && first_step)) {
	    if (do_sph) {
		if (do_winds) WindOutput(SPHbtab, SPHnobj, windbtab, 
					 windnobj, outnamebase, iter);
		else SPHOutput(SPHbtab, SPHnobj, outnamebase, iter);
	    }
	    if (has_grav_data) Output(btab, nobj, outnamebase, iter);
	    if (do_point_mass) Output(pmtab, PMnobj, outnamebase, iter);
	}


	if (ForceStop()) {
	    singlPrintf("Stopping.\n");
	    break;
	}

	Msgf(("integrating positions\n"));
	dark_ke = dark_pe = 0.0;
	if (has_grav_data) {
	  for (p = btab; p < btab+nobj; p++) {
	    dark_ke += 0.5 * p->mass * Dot(p->vel, p->vel);
	    dark_pe += 0.5 * p->mass * p->phi;
	  }
	}
	if (do_point_mass) {
	  for (p = pmtab; p < pmtab+PMnobj; p++) {
	    dark_ke += 0.5 * p->mass * Dot(p->vel, p->vel);
	    dark_pe += 0.5 * p->mass * p->phi;
	  }
	}
	Fix_h(SPHbtab, SPHnobj, nbrcut_max, nbrcut_min, nbrcut_fac, max_h, min_h);
	ke = pe = te = 0.0;
	SPHnupdate = 0;
	for (q = SPHbtab; q < SPHbtab+SPHnobj; q++) {
	    ke += 0.5 * q->mass * Dot(q->vel, q->vel);
	    pe += 0.5 * q->mass * q->phi;
	    te += q->mass * q->u;
	    if (SPH_need_update(q)) {
		q->tacc = tpos;
		SPHnupdate++;
	    }
	}
	tposlast = tpos;
	if (first_step) {
	    for (i = 0; i < nobj; i++) {
		VVV(btab[i].pos_last, = btab[i].pos, - dt * btab[i].vel);
	    }
	    for (i = 0; i < PMnobj; i++) {
		VVV(pmtab[i].pos_last, = pmtab[i].pos, - dt * pmtab[i].vel);
	    }
	    for (i = 0; i < SPHnobj; i++) {
	      /* Roundoff problems happen here if dt is small and */
	      /* pos is large and sigle precision */
		VVV(SPHbtab[i].pos_last, = SPHbtab[i].pos,- dt*SPHbtab[i].vel);
		SPHbtab[i].udot_last = SPHbtab[i].udot;
	    }
	}
	else if (added_particles) {
	    for (i = 0; i < SPHnobj; i++) {
	      /* Roundoff problems happen here if dt is small and */
	      /* pos is large and sigle precision */
		if (SPHbtab[i].windid < windnobj) {
		    SPHbtab[i].udot_last = SPHbtab[i].udot;
		}
	    }
	    added_particles = 0;
	}
	/* One must be careful with this integration scheme, since v */
	/* is a derived variable.  To really adjust v, change pos_last */
	PUpdateV(btab[0].vel, stride, btab[0].pos, stride, btab[0].pos_last, 
		 stride, btab[0].acc, stride, nobj, dt, dt_last);
	/* v must be done before x, since pos_last is changed in PUpX */
	PUpdateX(btab[0].pos, stride, btab[0].pos_last, stride,
		 btab[0].acc, stride, nobj, dt, dt_last);

	PUpdateV(pmtab[0].vel, stride, pmtab[0].pos, stride, 
		 pmtab[0].pos_last, stride, pmtab[0].acc, stride, PMnobj, 
		 dt, dt_last);
	PUpdateX(pmtab[0].pos, stride, pmtab[0].pos_last, stride,
		 pmtab[0].acc, stride, PMnobj, dt, dt_last);

	ABUpdateXs(&SPHbtab[0].u, SPHstride, &SPHbtab[0].udot, SPHstride, 
		   &SPHbtab[0].udot_last, SPHstride, &(SPHbtab[0].ident),
		   SPHstride3, SPHnobj, dt, dt_last);
	/* One must be careful with this integration scheme, since v */
	/* is a derived variable.  To really adjust v, change pos_last */
#ifdef POS_IS_DOUBLE
	PUpdateVsd(SPHbtab[0].vel, SPHstride, SPHbtab[0].pos, SPHstride2, 
		   SPHbtab[0].pos_last, SPHstride2, SPHbtab[0].acc, SPHstride, 
		   &(SPHbtab[0].ident), SPHstride3, SPHnobj, dt, dt_last);
	/* v must be done before x, since pos_last is changed in PUpX */
	PUpdateXsd(SPHbtab[0].pos, SPHstride2, SPHbtab[0].pos_last, SPHstride2,
		   SPHbtab[0].acc, SPHstride, &(SPHbtab[0].ident), SPHstride3, 
		   SPHnobj, dt, dt_last);
#else
	PUpdateVs(SPHbtab[0].vel, SPHstride, SPHbtab[0].pos, SPHstride, 
		 SPHbtab[0].pos_last, SPHstride, SPHbtab[0].acc, SPHstride, 
		 &(SPHbtab[0].ident), SPHstride3, SPHnobj, dt, dt_last);
	/* v must be done before x, since pos_last is changed in PUpX */
	PUpdateXs(SPHbtab[0].pos, SPHstride, SPHbtab[0].pos_last, SPHstride,
		 SPHbtab[0].acc, SPHstride, &(SPHbtab[0].ident), SPHstride3, 
		  SPHnobj, dt, dt_last);
#endif
	UpdateSXs(&SPHbtab[0].h, SPHstride, &SPHbtab[0].hdot, SPHstride, 
		  &(SPHbtab[0].ident), SPHstride3, SPHnobj, dt, dt_last);
	tpos += dt;
	tvel += dt;
	dt_last = dt;
	for (i = 0; i < SPHnobj; i++) {
	    VV(SPHbtab[i].acc_last, = SPHbtab[i].acc);
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
			 log_time, tpos, dt_last);
	    SPHWrapPeriodic(SPHbtab, SPHnobj, rmin, rmax, 2.0*sysradius, cosmology,
			 log_time, tpos, dt_last);
	}
	
	if (cosmology) 
	  singlPrintf("\ntpos: %g znow: %.3f iter: %d size: %.2f, eps: %.0f\n", 
		      tposlast, Znow(tposlast), iter, sysradius, this_eps);
	else
	  singlPrintf("\ntpos: %g iter: %d size: %f\n",
		      tposlast, iter, sysradius);

	etot = 0.0;
	if (has_grav_data) Diags(btab, nobj, dark_ke, dark_pe, &etot, dt_last, iter, gnobj);
	if (do_point_mass) Diags(pmtab, PMnobj, dark_ke, dark_pe, &etot, dt_last, iter, PMgnobj);
	if (do_sph) SPHDiags(SPHbtab, SPHnobj, ke, pe, te, &etot, dt_last, iter, SPHgnobj, 
	      &tmin, &tbad);

	MPMY_Combine(udot_limit, udot_limit, 2, MPMY_INT, MPMY_SUM);

	if (adaptive_dt) Fix_dt(&dt, (dark_dt < dt_max) ? dark_dt : dt_max, 
				tlow_cut, tmin, tbad, dt_short, dt_long,
				udot_limit[0], udot_limit[1]);

	singlPrintf("udot_limit high: %d low: %d\n", udot_limit[0], 
		    udot_limit[1]);
	singlPrintf("Total Energy: %g\n", etot);

	StopTimer(&StepTot);
	StopTimer(&StepTotWC);

	AddCounter(&HeapCnt_, malloc_heapsz()/1024);
		    
	if( timer_freq && iter%timer_freq == 0 ){
	    OutputTimers(singlPrintf);
	    OutputCounters(singlPrintf);
	    if( Msg_test("timers") ){
		/* This can be very tedious on a big machine. */
		OutputIndividualTimers(Msg_do);
		OutputIndividualCounters(Msg_do);
	    }
	    if (ntimer_detail) {
		struct {
		    int node;
		    float grav_tm;
		    float sph_tm;
		    float imbal_tm;
		    float per_tm;
		    int nterms;
		    int nbody;
		    int SPHnbody;
		    int SPHnupdate;
		} perf, *gp = 0;

		perf.node = MPMY_Procnum();
		perf.grav_tm = ReadTimer(&GravTm);
		perf.sph_tm = ReadTimer(&ForceSPH);
		perf.per_tm = ReadTimer(&PerTm);
		perf.nterms = ReadCounter(&NtermsCnt);
		perf.nbody = nobj;
		perf.SPHnbody = SPHnobj;
		perf.SPHnupdate = SPHnupdate;

		if (MPMY_Procnum() == 0) {
		    gp = Malloc(MPMY_Nproc() * sizeof(perf));
		    MPMY_Gather(&perf, sizeof(perf), MPMY_CHAR, gp, 0);
		    for (i = 0; i < MPMY_Nproc(); i++) 
		      singlPrintf("%3d %8.2f %8.2f %10d %6d %6d %6d\n",
				  gp[i].node, gp[i].sph_tm,  gp[i].imbal_tm, 
				  gp[i].nterms, gp[i].nbody, gp[i].SPHnbody,
				  gp[i].SPHnupdate);
		    Free(gp);
		} else {
		    MPMY_Gather(&perf, sizeof(perf), MPMY_CHAR, gp, 0);
		}
	    }
	} else {
	    OutputTimer(&StepTot, singlPrintf);
	    OutputTimer(&StepTotWC, singlPrintf);
	}
	MPMY_Combine(&SPHnupdate, &SPHnupdate, 1, MPMY_INT, MPMY_SUM);
	singlPrintf("Updated %d SPH accs\n", SPHnupdate);
	singlFflush();

	/* This can greatly improve the load balance */
	if (CWfac != 0.0) {	/* CWfac = 1 seems to work well for do_periodic */
	    gnterms /= gnobj;
	    singlPrintf("Avg nterms = %.0f, CWfac is %.2f\n", gnterms, CWfac);
	    gnterms *= CWfac;
	    /* Account for 'constant' work associated with each particle */
	    for (p = btab; p < btab+nobj; p++) 
	      p->nterms += gnterms;
	}

	first_step = 0;
	if( Msg_test("memleak") ){
	    Msg_do("Memory map after iteration %d\n", iter);
	    malloc_print();
	}
    }
    singlPrintf("Bye!\n");
    Msgf(("Bye!\n"));
    Msg_flush();
    exit(0);			/* trex seems to hang in __exit() */
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
    EnableTimer(&MACTm, "MAC Time");
    EnableTimer(&ForceSPH, "Force (SPH)");
    EnableTimer(&RhoSPH, "Rho (SPH)");
    EnableTimer(&WalkDeferTm, "Walk Defer");
    EnableTimer(&WTermTm, "WalkTerm");
    EnableTimer(&WNTTm, "WalkNT");
    EnableTimer(&PerTm, "Periodic");
    EnableTimer(&PerTmSPH, "Periodic SPH");
#if 0
    EnableCounter(&CCInt, "Cell-cell");
    EnableCounter(&BCInt, "Body-cell");
    EnableCounter(&CBInt, "Cell-body");
    EnableCounter(&BBInt, "Body-body");
#endif
    EnableCounter(&NtermsCnt, "Nterms");
    EnableCounter(&NbodyCnt, "Nbody");
    EnableCounter(&SPHbodyCnt, "SPHbody");
    EnableCounter(&CCIntRej, "MAC fail");
    EnableCounter(&SharedCnt, "Shared Cells");
    EnableCounter(&TranslateCnt, "Translate");
    EnableCounter(&DeferCnt, "Deferred");
    EnableCounter(&MemCnt, "Mem Used (K)");
    EnableCounter(&HeapCnt_, "Heap Sz (K)");
    return csdfp;
}

static void SanityCheck(body *btab, int nobj, int gnobj, double *mtotp){
    double mtot;
    body *p;
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

static void SPHSanityCheck(SPHbody *btab, int nobj, int gnobj, double *mtotp){
    double mtot;
    SPHbody *p;
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
    Msgf(("SPH Particle 0 (%d), %g, %g %g %g, %g %g %g\n",
	  btab->ident, btab->mass, 
	  btab->pos[0], btab->pos[1], btab->pos[2],
	  btab->vel[0], btab->vel[1], btab->vel[2]));
    Msgf(("SPH Particle %d (%d), %g, %g %g %g, %g %g %g\n", nobj-1,
	  btab[nobj-1].ident, btab[nobj-1].mass, 
	  btab[nobj-1].pos[0], btab[nobj-1].pos[1], btab[nobj-1].pos[2],
	  btab[nobj-1].vel[0], btab[nobj-1].vel[1], btab[nobj-1].vel[2]));
    singlPrintf("SPHgnobj = %d, SPHmtot = %f\n", gnobj, mtot);
    *mtotp = mtot;
}

static void
FixGlobalForce(body *xptr, int n){
    /* Make whatever corrections are necessary to the acceleration, etc.
       based on values of GNewt, Lambda, etc., etc. */
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
PeriodicSPH(tree_t *tp, float width, float vwidth)
{
    int i, j, k;
    float offset[NDIM];
    float voffset[NDIM];

    /* if doing cosmology, SPH needs velocities wrapped */
    
    for (i = -1; i <= 1; i++) {
        offset[0] = i * width;
        voffset[0] = i * vwidth;
        for (j = -1; j <= 1; j++) {
            offset[1] = j * width;
            voffset[1] = j * vwidth;
            for (k = -1; k <= 1; k++) {
                offset[2] = k * width;
                voffset[2] = k * vwidth;
		SetSPHOffset(offset, voffset);
                if (i || j || k) WalkNT(tp);
            }
        }
    }
    UnSetSPHOffset();
}

static void
WrapPeriodic(body *bp, int n, float *rmin, float *rmax, float sz, 
	     int cosmology, int log_time, float tvel, float dt)
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
    /* If using AB integrator, one must also wrap pos_last */
    if (!cosmology) {
	for(b=bp; b<&bp[n]; b++) {
	    VVVS(if LPAREN b->pos, > rmax, RPAREN b->pos_last, -= sz);
	    VVVS(if LPAREN b->pos, > rmax, RPAREN b->pos, -= sz);
	    VVVS(if LPAREN b->pos, < rmin, RPAREN b->pos_last, += sz);
	    VVVS(if LPAREN b->pos, < rmin, RPAREN b->pos, += sz);
	}
    } else {
	if (log_time) 
	  vsz = sz*Hnow(tvel)*1.5*pow((double)tvel, 1./3.); /* ?? */
	else 
	  vsz = sz*Hnow(tvel);
	for(b=bp; b<&bp[n]; b++) {
	    VVVS(if LPAREN b->pos, > rmax, RPAREN b->vel, -= vsz);
	    VVVS(if LPAREN b->pos, > rmax, RPAREN b->pos_last, -= sz-vsz*dt);
	    VVVS(if LPAREN b->pos, > rmax, RPAREN b->pos, -= sz);
	    VVVS(if LPAREN b->pos, < rmin, RPAREN b->vel, += vsz);
	    VVVS(if LPAREN b->pos, < rmin, RPAREN b->pos_last, += sz-vsz*dt);
	    VVVS(if LPAREN b->pos, < rmin, RPAREN b->pos, += sz);
	}
    }
}    

static void
SPHWrapPeriodic(SPHbody *bp, int n, float *rmin, float *rmax, float sz, 
	     int cosmology, int log_time, float tvel, float dt)
{
    SPHbody *b;
    int flux[NDIM] = {0, 0, 0};
    float vsz; /* hubble flow */
    
    for(b=bp; b<&bp[n]; b++) {
	VVVS(if LPAREN b->pos, > rmax, RPAREN flux, += 1);
	VVVS(if LPAREN b->pos, < rmin, RPAREN flux, -= 1);
    }
    MPMY_Combine(flux, flux, NDIM, MPMY_INT, MPMY_SUM);
    singlPrintf("SPHFlux %d %d %d\n", flux[0], flux[1], flux[2]);
    /* If using AB integrator, one must also wrap pos_last */
    if (!cosmology) {
	for(b=bp; b<&bp[n]; b++) {
	    VVVS(if LPAREN b->pos, > rmax, RPAREN b->pos_last, -= sz);
	    VVVS(if LPAREN b->pos, > rmax, RPAREN b->pos, -= sz);
	    VVVS(if LPAREN b->pos, < rmin, RPAREN b->pos_last, += sz);
	    VVVS(if LPAREN b->pos, < rmin, RPAREN b->pos, += sz);
	}
    } else {
	if (log_time) 
	  vsz = sz*Hnow(tvel)*1.5*pow((double)tvel, 1./3.); /* ?? */
	else 
	  vsz = sz*Hnow(tvel);
	for(b=bp; b<&bp[n]; b++) {
	    VVVS(if LPAREN b->pos, > rmax, RPAREN b->vel, -= vsz);
	    VVVS(if LPAREN b->pos, > rmax, RPAREN b->pos_last, -= sz-vsz*dt);
	    VVVS(if LPAREN b->pos, > rmax, RPAREN b->pos, -= sz);
	    VVVS(if LPAREN b->pos, < rmin, RPAREN b->vel, += vsz);
	    VVVS(if LPAREN b->pos, < rmin, RPAREN b->pos_last, += sz-vsz*dt);
	    VVVS(if LPAREN b->pos, < rmin, RPAREN b->pos, += sz);
	}
    }
}    

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

static void CosmoPush(struct cosmo_s *p, float time){
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


/* Write out nobj SPH particles; only one node should write wind particles */
static void WindOutput(SPHbody *btab, int nobj, windbody *windbtab, 
		       int windnobj, const char *outnamebase, int iter)
{
    SPHbody *p;
    int i;
    sortresult_t outputsort;
    SPHoutbody *output_btab;
    int output_nobj = nobj;
    float tpos_out = tpos;
    float tvel_out = tvel; /* changed in Integrate() */
    double ke, pe, te;
    MPMY_Comm_request req;
    int output_gnobj;
    float output_z, output_h, output_R0;
    char outname[256];

    sprintf(outname, "%s_sph.%04d", outnamebase, iter);
    pe = ke = te = 0.0;
    for (p = btab; p < btab+nobj; p++) {
	ke += (float)0.5 * p->mass * Dot(p->vel, p->vel);
	te += p->mass * p->u;
	pe += (float)0.5 * p->mass * p->phi;
    }
    output_btab = Malloc(output_nobj * sizeof(SPHoutbody));
    for(i=0; i<output_nobj; i++){
	output_btab[i].mass = btab[i].mass;
	VV(output_btab[i].pos, = btab[i].pos);
	VV(output_btab[i].vel, = btab[i].vel);
	output_btab[i].u = btab[i].u;
	output_btab[i].h = btab[i].h;
	output_btab[i].rho = btab[i].rho;
/*  	output_btab[i].drho_dt = btab[i].drho_dt; */
	output_btab[i].udot = btab[i].udot;
#ifdef SPH_SAVE_ACC
	VV(output_btab[i].acc, = btab[i].acc);
	VV(output_btab[i].acc_last, = btab[i].acc_last);
	output_btab[i].phi = btab[i].phi;
	output_btab[i].dt = btab[i].dt;
#endif
 	output_btab[i].nbrs = btab[i].nbrs;
	output_btab[i].ident = btab[i].ident;
	output_btab[i].windid = btab[i].windid;
    }
    Msgf(("Doing output of %d bodies\n", output_nobj));
    singlPrintf("Trying to sort output\n");
    pqsortsetup_order(&outputsort, output_btab, output_nobj,
		      sizeof(SPHoutbody), 0.1F, 1, Realloc_f);
    output_btab = pqsort(&outputsort, UnityCost, (pq_keyproto)SPHOutIdentKey);
    output_nobj = outputsort.nobj;
    Msgf(("After pqsort, %d outbodies\n", output_nobj));
    MPMY_ICombine_Init(&req);
    MPMY_ICombine(&ke, &ke, 1, MPMY_DOUBLE, MPMY_SUM, req);
    MPMY_ICombine(&pe, &pe, 1, MPMY_DOUBLE, MPMY_SUM, req);
    MPMY_ICombine(&te, &te, 1, MPMY_DOUBLE, MPMY_SUM, req);
    MPMY_ICombine(&output_nobj, &output_gnobj, 1, MPMY_INT, MPMY_SUM, req);
    MPMY_ICombine_Wait(req);
    if (cosmology) {
	output_z = Znow(tpos_out);
	output_h = Hnow(tpos_out);
	output_R0 = R0;
    } else {
	output_z = 0.0;
	output_h = 0.0;
	output_R0 = sysradius;
    }
    SDFwritewind(outname, output_gnobj, output_nobj, 
		 output_btab, windnobj, windbtab, sizeof(SPHoutbody), 
		 sizeof(windbody), WINDOUTBODYDESC, SPHOUTBODYDESC, 
		/* "npart", SDF_INT, output_gnobj, */
		"iter", SDF_INT, iter,
		"dt", SDF_FLOAT, dt,
		"eps", SDF_FLOAT, this_eps,
		"Gnewt", SDF_FLOAT, cosmo.GNewt,
		"tolerance", SDF_FLOAT, this_tol,
		"frac_tolerance", SDF_FLOAT, frac_tol,
		"ndim", SDF_INT, NDIM,
		"tpos", SDF_FLOAT, tpos_out,
		"tvel", SDF_FLOAT, tvel_out,
		"R0", SDF_FLOAT, output_R0,
		"Omega0", SDF_FLOAT, cosmo.Omega0,
		"H0", SDF_FLOAT, cosmo.H0,
		"Lambda_prime", SDF_FLOAT, cosmo.Lambda,
		"hubble", SDF_FLOAT, output_h,
		"redshift", SDF_FLOAT, output_z,
		"gamma", SDF_FLOAT, Gamma,
		 "centmass", SDF_FLOAT, centmass, 
		"ke", SDF_DOUBLE, ke,
		"pe", SDF_DOUBLE, pe,
		"te", SDF_DOUBLE, te,
		NULL);
    Free(output_btab);
    singlPrintf("\nOutput done.\n");
#ifndef __DELTA__
    if (MPMY_Procnum() == 0) {
	char name[256];
	sprintf(name, "%s_sph.restart", outnamebase);
	if (unlink(name))
	  Shout("unlink of %s failed, errno=%d\n", name, errno);
	if (symlink(outname, name))
	  Shout("symlink of %s failed, errno=%d\n", outname, errno);
    }
#endif
}


static void SPHOutput(SPHbody *btab, int nobj, const char *outnamebase, int iter)
{
    SPHbody *p;
    int i;
    sortresult_t outputsort;
    SPHoutbody *output_btab;
    int output_nobj = nobj;
    float tpos_out = tpos;
    float tvel_out = tvel; /* changed in Integrate() */
    double ke, pe, te;
    MPMY_Comm_request req;
    int output_gnobj;
    float output_z, output_h, output_R0;
    char outname[256];

    sprintf(outname, "%s_sph.%04d", outnamebase, iter);
    pe = ke = te = 0.0;
    for (p = btab; p < btab+nobj; p++) {
	ke += (float)0.5 * p->mass * Dot(p->vel, p->vel);
	te += p->mass * p->u;
	pe += (float)0.5 * p->mass * p->phi;
    }
    output_btab = Malloc(output_nobj * sizeof(SPHoutbody));
    for(i=0; i<output_nobj; i++){
	output_btab[i].mass = btab[i].mass;
	VV(output_btab[i].pos, = btab[i].pos);
	VV(output_btab[i].vel, = btab[i].vel);
	output_btab[i].u = btab[i].u;
	output_btab[i].h = btab[i].h;
	output_btab[i].rho = btab[i].rho;
/*  	output_btab[i].drho_dt = btab[i].drho_dt; */
	output_btab[i].udot = btab[i].udot;
#ifdef SPH_SAVE_ACC
	VV(output_btab[i].acc, = btab[i].acc);
	VV(output_btab[i].acc_last, = btab[i].acc_last);
	output_btab[i].phi = btab[i].phi;
	output_btab[i].dt = btab[i].dt;
#endif
 	output_btab[i].nbrs = btab[i].nbrs;
	output_btab[i].ident = btab[i].ident;
    }
/*     Msg("output", ("Doing output of %d bodies\n", output_nobj)); */
    Msgf(("Doing output of %d bodies\n", output_nobj));
    singlPrintf("Trying to sort output\n");
    pqsortsetup_order(&outputsort, output_btab, output_nobj,
		      sizeof(SPHoutbody), 0.1F, 1, Realloc_f);
    output_btab = pqsort(&outputsort, UnityCost, (pq_keyproto)SPHOutIdentKey);
    output_nobj = outputsort.nobj;
/*     Msg("output", ("After pqsort, %d outbodies\n", output_nobj)); */
    Msgf(("After pqsort, %d outbodies\n", output_nobj));
    MPMY_ICombine_Init(&req);
    MPMY_ICombine(&ke, &ke, 1, MPMY_DOUBLE, MPMY_SUM, req);
    MPMY_ICombine(&pe, &pe, 1, MPMY_DOUBLE, MPMY_SUM, req);
    MPMY_ICombine(&te, &te, 1, MPMY_DOUBLE, MPMY_SUM, req);
    MPMY_ICombine(&output_nobj, &output_gnobj, 1, MPMY_INT, MPMY_SUM, req);
    MPMY_ICombine_Wait(req);
    if (cosmology) {
	output_z = Znow(tpos_out);
	output_h = Hnow(tpos_out);
	output_R0 = R0;
    } else {
	output_z = 0.0;
	output_h = 0.0;
	output_R0 = sysradius;
    }
    SDFwrite(outname, output_gnobj, 
	     output_nobj, output_btab, sizeof(SPHoutbody),
	     SPHOUTBODYDESC,
	     "npart", SDF_INT, output_gnobj,
	     "iter", SDF_INT, iter,
	     "dt", SDF_FLOAT, dt,
	     "eps", SDF_FLOAT, this_eps,
	     "Gnewt", SDF_FLOAT, cosmo.GNewt,
	     "tolerance", SDF_FLOAT, this_tol,
	     "frac_tolerance", SDF_FLOAT, frac_tol,
	     "ndim", SDF_INT, NDIM,
	     "tpos", SDF_FLOAT, tpos_out,
	     "tvel", SDF_FLOAT, tvel_out,
	     "R0", SDF_FLOAT, output_R0,
	     "Omega0", SDF_FLOAT, cosmo.Omega0,
	     "H0", SDF_FLOAT, cosmo.H0,
	     "Lambda_prime", SDF_FLOAT, cosmo.Lambda,
	     "hubble", SDF_FLOAT, output_h,
	     "redshift", SDF_FLOAT, output_z,
	     "gamma", SDF_FLOAT, Gamma,
	     "centmass", SDF_FLOAT, centmass, 
	     "ke", SDF_DOUBLE, ke,
	     "pe", SDF_DOUBLE, pe,
	     "te", SDF_DOUBLE, te,
	     NULL);
    Free(output_btab);
    singlPrintf("\nOutput done.\n");
#ifndef __DELTA__
    if (MPMY_Procnum() == 0) {
	char name[256];
	sprintf(name, "%s_sph.restart", outnamebase);
	if (unlink(name))
	  Shout("unlink of %s failed, errno=%d\n", name, errno);
	if (symlink(outname, name))
	  Shout("symlink of %s failed, errno=%d\n", outname, errno);
    }
#endif
}

static void Output(body *btab, int nobj, const char *outnamebase, int iter)
{
    body *p;
    int i;
    sortresult_t outputsort;
    outbody *output_btab;
    int output_nobj = nobj;
    float tpos_out = tpos;
    float tvel_out = tvel; /* changed in Integrate() */
    double ke, pe;
    MPMY_Comm_request req;
    int output_gnobj;
    float output_z, output_h, output_R0;
    char outname[256];

    sprintf(outname, "%s.%04d", outnamebase, iter);
    pe = ke = 0.0;
    for (p = btab; p < btab+nobj; p++) {
	ke += (float)0.5 * p->mass * Dot(p->vel, p->vel);
	pe += (float)0.5 * p->mass * p->phi;
    }
    output_btab = Malloc(output_nobj * sizeof(outbody));
    for(i=0; i<output_nobj; i++){
	output_btab[i].mass = btab[i].mass;
	VV(output_btab[i].pos, = btab[i].pos);
	VV(output_btab[i].vel, = btab[i].vel);
#ifdef SAVE_ACC
	VV(output_btab[i].acc, = btab[i].acc);
	output_btab[i].phi = btab[i].phi;
#endif
	/* Added l and accmass; this is pretty specific to */
	/* point masses */
	VV(output_btab[i].l, = btab[i].l);
	output_btab[i].accmass = btab[i].accmass;	
	output_btab[i].ident = btab[i].ident;
    }
    Msg("output", ("Doing output of %d bodies\n", output_nobj));
    singlPrintf("Trying to sort output\n");
    pqsortsetup_order(&outputsort, output_btab, output_nobj,
		      sizeof(outbody), 0.1F, 1, Realloc_f);
    output_btab = pqsort(&outputsort,
			 (pq_wgtproto)UnityCost, 
			 (pq_keyproto)OutIdentKey);
    output_nobj = outputsort.nobj;
    Msg("output", ("After pqsort, %d outbodies\n", output_nobj));
    MPMY_ICombine_Init(&req);
    MPMY_ICombine(&ke, &ke, 1, MPMY_DOUBLE, MPMY_SUM, req);
    MPMY_ICombine(&pe, &pe, 1, MPMY_DOUBLE, MPMY_SUM, req);
    MPMY_ICombine(&output_nobj, &output_gnobj, 1, MPMY_INT, MPMY_SUM, req);
    MPMY_ICombine_Wait(req);
    if (cosmology) {
	output_z = Znow(tpos_out);
	output_h = Hnow(tpos_out);
	output_R0 = R0;
    } else {
	output_z = 0.0;
	output_h = 0.0;
	output_R0 = sysradius;
    }
    SDFwrite(outname, output_gnobj, 
	     output_nobj, output_btab, sizeof(outbody),
	     OUTBODYDESC,
	     "npart", SDF_INT, output_gnobj,
	     "iter", SDF_INT, iter,
	     "dt", SDF_FLOAT, dt,
	     "eps", SDF_FLOAT, this_eps,
	     "Gnewt", SDF_FLOAT, cosmo.GNewt,
	     "tolerance", SDF_FLOAT, this_tol,
	     "frac_tolerance", SDF_FLOAT, frac_tol,
	     "iter", SDF_INT, iter,
	     "ndim", SDF_INT, NDIM,
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
    Free(output_btab);
    singlPrintf("\nOutput done.\n");
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

static void 
Fix_h(SPHbody *btab, int nobj, int nbrcut_max, int nbrcut_min,
		  float nbrcut_fac, float max_h, float min_h)
{
    SPHbody *p;
    int nn[6];
    
    nn[0] = nn[1] = nn[2] = nn[3] = nn[4] = nn[5] = 0;
    for (p = btab; p < btab+nobj; p++) {
	if (!SPH_need_update(p)) continue;
	if (p->nbrs > 8*nbrcut_max) {
	    p->h *= 0.75;
	    nn[2]++;
	} else if (p->nbrs > 2*nbrcut_max) {
	    p->h *= pow((double)nbrcut_max/p->nbrs, (1./3.));
	    nn[3]++;
	} else if (p->nbrs > nbrcut_max) {
	    p->h -= nbrcut_fac * p->h;
	    nn[0]++;
	}
	
	if (p->nbrs < nbrcut_min) {
	    p->h += nbrcut_fac * p->h;
	    nn[1]++;
	}

	if (p->h >= max_h) {
	    p->h = max_h;
	    nn[4]++;
	}
	if (p->h <= min_h) {
	    p->h = min_h;
	    nn[5]++;
	}
    }
    MPMY_Combine(&nn, &nn, 6, MPMY_INT, MPMY_SUM);
    singlPrintf("Nbr_cuts: over: %d under: %d 2x_over: %d 8x_over: %d\n",
		nn[0]+nn[2]+nn[3], nn[1], nn[2]+nn[3], nn[2]);
    singlPrintf("Smoothing cuts: max: %d min: %d\n", nn[4], nn[5]);
}

static void
Diags(body *btab, int nobj, double ke, double pe, double *etot,
      float dt_last, int iter, int gnobj)
{
    double force[NDIM];		/* Things accumulated over all particles */
    double com[NDIM], comv[NDIM];
    double acc2;		/* must be double precision */
    double avg_nbrs;
    double mtot;
    float sacc2;
    int local_gnobj;
    body *p;
    MPMY_Comm_request req;

    VS(force, = 0.0);
    VS(com, = 0.0);
    VS(comv, = 0.0);
    acc2 = 0.0;
    mtot = 0.0;
    gnterms = 0.0;
    avg_nbrs = 0.0;
    for (p = btab; p < btab+nobj; p++) {
	VV(com, += p->mass*p->pos);
	VV(comv, += p->mass*p->vel);
	VV(force, += p->mass*p->acc);
	sacc2 = Dot(p->acc, p->acc);
	acc2 += sacc2;
	mtot += p->mass;
	gnterms += p->nterms;
	if (p->nterms <= 0) SeriousWarning("nterms is %f\n", p->nterms);
    }
    AddCounter(&NtermsCnt, (int)gnterms); /* might overflow?? */
    
    Msgf(("doing MPMY_combine\n"));
    MPMY_ICombine_Init(&req);
    MPMY_ICombine(force, force, NDIM, MPMY_DOUBLE, MPMY_SUM, req);
    MPMY_ICombine(com, com, NDIM, MPMY_DOUBLE, MPMY_SUM, req);
    MPMY_ICombine(comv, comv, NDIM, MPMY_DOUBLE, MPMY_SUM, req);
    MPMY_ICombine(&acc2, &acc2, 1, MPMY_DOUBLE, MPMY_SUM, req);
    MPMY_ICombine(&ke, &ke, 1, MPMY_DOUBLE, MPMY_SUM, req);
    MPMY_ICombine(&pe, &pe, 1, MPMY_DOUBLE, MPMY_SUM, req);
    MPMY_ICombine(&mtot, &mtot, 1, MPMY_DOUBLE, MPMY_SUM, req);
    MPMY_ICombine(&gnterms, &gnterms, 1, MPMY_DOUBLE, MPMY_SUM, req);
    MPMY_ICombine(&nobj, &local_gnobj, 1, MPMY_INT, MPMY_SUM, req);
    MPMY_ICombine_Wait(req);
    Msgf(("done MPMY_combine\n"));
    
    if (gnobj != local_gnobj) {
	Error("gnobj not consistent %d vs %d\n", gnobj, local_gnobj);
    }
    singlPrintf("dark_ke: %g dark_pe: %g dark_energy: %g\n", 
		ke, pe, ke+pe);
    avg_nbrs /= gnobj;
    VS(force, /= mtot);
    VS(com, /= mtot);
    VS(comv, /= mtot);
    singlPrintf("dark CM accel: (" Sinfix("%g", " ") "): %g\n",
		Vinfix(force, COMMA), sqrt(Dot(force, force)));
    singlPrintf("dark rms accel: %g\n", sqrt(acc2/gnobj));
    *etot += ke+pe;
}

static void
SPHDiags(SPHbody *btab, int nobj, double ke, double pe, double te, double *etot, 
	 float dt_last, int iter, int gnobj, float *tmin, int *tbad)
{
    double force[NDIM];		/* Things accumulated over all particles */
    double com[NDIM], comv[NDIM];
    double acc2;		/* must be double precision */
    double avg_nbrs;
    double mtot;
    float max_vsound;
    float max_h, min_h, max_rho, min_rho, max_u, min_u;
    float rho_err, rms_rho_err2, max_rho_err;
    float min_dt;
    float sacc2;
    float tx;
    float dti;
    int min_nbrs, max_nbrs;
    int tlow;
    int local_gnobj;
    SPHbody *p;
    MPMY_Comm_request req;

    max_vsound = (float)0.0;
    rms_rho_err2 = max_rho_err = (float)0.0;
    max_h = max_rho = max_u = (float)0.0;
    min_h = min_rho = min_u = 1e30;
    min_dt = 1e30;
    max_nbrs = 0;
    min_nbrs = 10000;
    VS(force, = 0.0);
    VS(com, = 0.0);
    VS(comv, = 0.0);
    acc2 = 0.0;
    mtot = 0.0;
    gnterms = 0.0;
    avg_nbrs = 0.0;
    tlow = 0;
    for (p = btab; p < btab+nobj; p++) {
	VV(com, += p->mass*p->pos);
	VV(comv, += p->mass*p->vel);
	VV(force, += p->mass*p->acc);
	sacc2 = Dot(p->acc, p->acc);
	acc2 += sacc2;
	mtot += p->mass;
	gnterms += p->nterms;
	avg_nbrs += p->nbrs;
	if (p->nterms <= 0) SeriousWarning("nterms is %f\n", p->nterms);
	if (p->h > max_h) max_h = p->h;
	if (p->h < min_h) min_h = p->h;
	if (p->vsound/p->h > max_vsound)
	  max_vsound = p->vsound/p->h;
	dti = p->h/p->vsound;
	tx = p->rho/fabs(p->drho_dt);
	if (tx < dti) dti = tx;
	tx = p->u/fabs(p->udot);
	if (tx < dti) dti = tx;
	dti *= courant_number;
	if (p->min_nbr_dt == 1e30) { 
	   /* This could happen if there are no nbrs. */
	    SeriousWarning("Ignoring min_nbr_dt of %g\n", p->min_nbr_dt);
	} else {
	    p->dt = p->min_nbr_dt;
	}
	if (dti > dt_max && 2.0*p->dt < dt_max) p->dt_next = 2.0*p->dt;
	else if (dti > 2.0*p->dt) p->dt_next = 2.0*p->dt;
	else if (dti < p->dt) {
	    p->dt *= 0.5;
	    p->dt_next = p->dt;
	}
	if (dti < min_dt)
	  min_dt = dti;
	if (dti < dt)
	  tlow++;
	if (p->dt_next > dt_max) p->dt_next = dt_max;
	if (p->nbrs > max_nbrs)
	  max_nbrs = p->nbrs;
	if (p->nbrs < min_nbrs)
	  min_nbrs = p->nbrs;
	if (p->rho > max_rho) max_rho = p->rho;
	if (p->rho < min_rho) min_rho = p->rho;
	if (p->u > max_u) max_u = p->u;
	if (p->u < min_u) min_u = p->u;
	if (p->u < 0.0) {
	    SeriousWarning("Iter %d: particle has negative energy\n%s\n", 
			   iter, PrintSPHBodyContents(p));
	    p->u = 0.0;
	}
	rho_err = fabs(fabs(p->rho_est/p->rho) - (float)1.0);
	rms_rho_err2 += rho_err*rho_err;
	if (rho_err > max_rho_err) max_rho_err = rho_err;
    }
    AddCounter(&NtermsCnt, (int)gnterms); /* might overflow?? */
    
    Msgf(("doing MPMY_combine\n"));
    MPMY_ICombine_Init(&req);
    MPMY_ICombine(force, force, NDIM, MPMY_DOUBLE, MPMY_SUM, req);
    MPMY_ICombine(com, com, NDIM, MPMY_DOUBLE, MPMY_SUM, req);
    MPMY_ICombine(comv, comv, NDIM, MPMY_DOUBLE, MPMY_SUM, req);
    MPMY_ICombine(&acc2, &acc2, 1, MPMY_DOUBLE, MPMY_SUM, req);
    MPMY_ICombine(&ke, &ke, 1, MPMY_DOUBLE, MPMY_SUM, req);
    MPMY_ICombine(&pe, &pe, 1, MPMY_DOUBLE, MPMY_SUM, req);
    MPMY_ICombine(&te, &te, 1, MPMY_DOUBLE, MPMY_SUM, req);
    MPMY_ICombine(&mtot, &mtot, 1, MPMY_DOUBLE, MPMY_SUM, req);
    MPMY_ICombine(&gnterms, &gnterms, 1, MPMY_DOUBLE, MPMY_SUM, req);
    MPMY_ICombine(&avg_nbrs, &avg_nbrs, 1, MPMY_DOUBLE, MPMY_SUM, req);
    MPMY_ICombine(&max_vsound, &max_vsound, 1, MPMY_FLOAT, MPMY_MAX, req);
    MPMY_ICombine(&max_h, &max_h, 1, MPMY_FLOAT, MPMY_MAX, req);
    MPMY_ICombine(&max_rho, &max_rho, 1, MPMY_FLOAT, MPMY_MAX, req);
    MPMY_ICombine(&max_u, &max_u, 1, MPMY_FLOAT, MPMY_MAX, req);
    MPMY_ICombine(&min_h, &min_h, 1, MPMY_FLOAT, MPMY_MIN, req);
    MPMY_ICombine(&min_rho, &min_rho, 1, MPMY_FLOAT, MPMY_MIN, req);
    MPMY_ICombine(&min_u, &min_u, 1, MPMY_FLOAT, MPMY_MIN, req);
    MPMY_ICombine(&min_dt, &min_dt, 1, MPMY_FLOAT, MPMY_MIN, req);
    MPMY_ICombine(&rms_rho_err2, &rms_rho_err2, 1, MPMY_FLOAT, MPMY_SUM, req);
    MPMY_ICombine(&max_rho_err, &max_rho_err, 1, MPMY_FLOAT, MPMY_MAX, req);
    MPMY_ICombine(&nobj, &local_gnobj, 1, MPMY_INT, MPMY_SUM, req);
    MPMY_ICombine(&max_nbrs, &max_nbrs, 1, MPMY_INT, MPMY_MAX, req);
    MPMY_ICombine(&min_nbrs, &min_nbrs, 1, MPMY_INT, MPMY_MIN, req);
    MPMY_ICombine(&tlow, &tlow, 1, MPMY_INT, MPMY_SUM, req);
    MPMY_ICombine_Wait(req);
    Msgf(("done MPMY_combine\n"));
    
    if (gnobj != local_gnobj) {
	Error("gnobj not consistent %d vs %d\n", gnobj, local_gnobj);
    }
    singlPrintf("ke: %g pe: %g te: %g energy: %g\n", 
		ke, pe, te, ke+pe+te);
    avg_nbrs /= gnobj;
    VS(force, /= mtot);
    VS(com, /= mtot);
    VS(comv, /= mtot);
    singlPrintf("CM accel: (" Sinfix("%g", " ") "): %g\n",
		Vinfix(force, COMMA), sqrt(Dot(force, force)));
    singlPrintf("centmass: %g\n", centmass);
    singlPrintf("rms accel: %g\n", sqrt(acc2/gnobj));
    singlPrintf("dt: %g min_dt: %g h_vsound: %g bad_courant: %d\n",
		dt_last, min_dt, 1.0/max_vsound, tlow);
    singlPrintf("max_h: %g min_h: %g max_rho: %g min_rho: %g\n",
		max_h, min_h, max_rho, min_rho);
    singlPrintf("max_u: %g min_u: %g rms_rho_err:%g max_rho_err: %g\n", 
		max_u, min_u, sqrt(rms_rho_err2/gnobj), max_rho_err);
    singlPrintf("max_nbrs: %d min_nbrs: %d avg_nbrs: %.0f\n", 
		max_nbrs, min_nbrs, avg_nbrs);
    *etot += ke+pe+te;
    *tmin = min_dt;
    *tbad = tlow;
}

static void
Fix_dt(float *dt, float dt_max, int tlow_cut, float tmin, int tbad, 
       int dtshort, int dtlong, int limit_high, int limit_low)
{
    static int dtlongvote;
    static int dtshortvote;

    if (tmin > 4**dt) dtlongvote += 5;
    if (tmin > 2.0**dt) {
	dtshortvote = 0;
	dtlongvote++;
    }
    if (tbad < tlow_cut/10 || tbad == 0) {  /* Integer division, buddy */
	dtlongvote++;
    } else if (tbad > tlow_cut/2) {
	dtlongvote--;
    }
    if (tbad >= tlow_cut || limit_high > 0 || limit_low > 0 ) {
	dtlongvote = 0;
	dtshortvote++;
    }

    /* singlPrintf("Votes: short: %d; long %d\n", dtshortvote, dtlongvote); */
    if (dtshortvote > dtshort) {
	singlPrintf(("Adjusting dt down by factor of 1/2\n"));
	*dt *= (float)(1./2.);
	dtshortvote = dtlongvote = 0;
    } else if (dtlongvote > dtlong && (4./3.)**dt <= dt_max) {
	singlPrintf(("Adjusting dt up by factor of 2\n"));
	*dt *= (float)(2./1.);
	dtshortvote = dtlongvote = 0;
    }
}

static void
ReadCosmo(SDF *sdfp, struct cosmo_s *cosmo, float tpos, float *R0p)
{
    float Z;

    cosmo->t = tpos;
    SDFgetfloatOrDefault(sdfp, "Omega0",  &cosmo->Omega0, (float)1.0);
    SDFgetfloatOrDefault(sdfp, "Lambda_prime",  &cosmo->Lambda, (float)0.0);
    /* default is for h_100 = 0.5 */
    SDFgetfloatOrDefault(sdfp, "H0",  &cosmo->H0, (float)0.0511365);
    if( SDFhasname("box_size", sdfp) ) {
	SDFgetfloatOrDie(sdfp, "box_size",  R0p);
	*R0p /= 2.0;
    } else
      SDFgetfloatOrDie(sdfp, "R0",  &R0);
    
    /* Now we need to get initial values for cosmo->a */
    if( SDFhasname("redshift", sdfp) ){
	SDFgetfloat(sdfp, "redshift", &Z);
	cosmo->a = 1.F/(1.F + Z);
    }else{
	if( cosmo->Omega0 == 1.0F ){
	    cosmo->a = pow( 1.5*cosmo->t*cosmo->H0, 2./3.);
	}else{
	    SinglError("Sorry.  Tell me the redshift in the data file\n");
	}
    }
    /* The Zel'dovich 'f' factor is only needed for setting initial
       velocities.  At this point, we don't know if we will be asked
       to do setpvel, though, so we read it anyway. */
    if( SDFhasname("velocity_fac", sdfp) ){
	SDFgetfloatOrDie(sdfp, "velocity_fac", &cosmo->Zel_f);
    }else{
	cosmo->Zel_f = 1.F;
    }
    singlPrintf("float H0 = %g\n", cosmo->H0);
    singlPrintf("float Omega0 = %g\n", cosmo->Omega0);
    singlPrintf("float redshift = %g\n", Znow(tpos));
    singlPrintf("float Lambda = %g\n", cosmo->Lambda);
    singlPrintf("float Zel_f = %g\n", cosmo->Zel_f);
}

/* This is a mess. physics_generic.c needs a better abstraction */
Key_t SPHGetKey(const void *p)
{
    body t;
    VV(t.pos, = ((SPHbody *)p)->pos);
    return GETKEY(&t);
}


static int
dark_need_update(float dark_tacc, float dark_dt)
{
    if (!dark_independent_dt) return 1;
    return (dark_tacc + dark_dt <= tpos + dt * 1.00001);
}

/*  static float  */
/*  IdtSPHGetCost(const SPHbody *ptr) */
/*  { */
/*      if (SPH_need_update(ptr)) */
/*        return (float) ptr->nterms; */
/*      else */
/*        return (float) default_nterms; */
/*  } */

int
SPH_need_update(const SPHbody *p)
{
    if (!independent_dt) return 1;
    return (p->tacc + p->dt <= tpos + dt * 1.00001);
}

