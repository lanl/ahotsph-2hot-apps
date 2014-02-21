/*
 * Copyright 2012-2013 Michael S. Warren. All Rights Reserved.
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
#include "Malloc.h"
#include "SDFwrite.h"
#include "SDFread.h"
#include "order.h"
#include "physics.h"
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
#include "decomp.h"
#include "image.h"
#include "memfile.h"
#include "ewald_le.h"
#include "cosmo.h"
#include "integrate.h"
#include "output.h"
#include "cpu.h"
#include "version.h"

static void SanityCheck(bodyptr btab, int nobj, int64_t gnobj, double *mtotp);
static void set_vels(body *p, int n, float dp, cosmology *c, 
		     double *tpos, double *tvel, int do_periodic);
static void set_displacement(body *p, int n, double dt, cosmology *c,
			     double *tpos, double *tvel, float noise, int do_periodic, ran_state *ranstate);
static SDF *startup(int argc, char **argv);
static void WrapPeriodic(body *bp, int n, float *rmin, float *rmax);
static void FixCubeEwaldLE(body *b, int nobj, const float *l, float gm, int nimage,
			   int subtract_background);
static void FixGlobalForce(body *bp, int n, cosmology *c, float real_time, int do_periodic);
static void WriteLightCone(body *xptr, const int n, const double dt, const double dtv, 
			   cosmology *c, double tpos, double tvel, char *name, char *tag,
			   int iter, float *origin, float *lc_rmax);
int maxmem(void);
int maxheap(void);

struct CPU_s CPU;
Timer_t StepTot, StepTotWC, BuildTot;
Timer_t FindForcesTm;
Timer_t FixCubeTm, LightConeTm, LightConeOpenTm, LightConeWriteTm;
Counter_t NbodyCnt;
Counter_t MemCnt, Ncell, Nquadcell, Nhexacell;
Counter_t Ncell_local, Nquadcell_local, Nhexacell_local, Ntbody, Nhcell;
Counter_t HeapCnt_;	/* HeapCnt is in the SunOS name space?! */
Counter_t KNtermsCnt, Scycles, Mcycles, Qcycles, FQcycles, Hcycles, FHcycles;

Timer_t WITm,WTermTm,WNTTm,WITm;
extern Timer_t SharedCellsWaitTm; /* in tree.c */
extern Timer_t PQSortCommTm, PQSortWaitTm; /* in pqsort.c */
extern Timer_t PQSortAtoaTm, PQSortAtoavTm;
extern Counter_t PQSortSends, PQSortRecvs, PQSortMaxn;

extern Timer_t GravQFTm, GravHFTm;;

#ifdef CUDA
static int has_cuda = 1;
#else
static int has_cuda = 0;
#endif

int
main(int argc, char *argv[])
{
    int64_t gnobj;
    int nobj;
    bodyptr btab;
    float eps;			/* smoothing scale */
    int i;
    float rmin[NDIM], rmax[NDIM];
    float sysradius[NDIM];
    float dt, dt_base, dt_hiz_tol;
    int nsteps;
    int first_step = 1;
    int do_output, do_checkpoint;
    int checkpoint_steps_interval, checkpoint_wallclock_interval;
    double subsample_fraction;
    int subsample_steps_interval, subsample_random_seed;
    int job_max_wallclock, step_wallclock_estimate, output_wallclock_estimate;
    int steps_to_next_output;
    int timer_freq;
    float sort_tol;
    int iter, start_iter;
    bodyptr p;
    int force_smoothing_type;
    float comov_eps, comov_eps_Zmin;
    float this_eps, this_eps_scaled;
    float CWfac;
    int ntimer_detail;
    int light_cone = 0;
    float lc_origin[NDIM] = {0.0f, 0.0f, 0.0f};
    char lc_name_suffix[256];
    int Ztol = 0;		/* if true, divide errtol by linear growth factor ^ Ztol */
    int setpvel = 0;
    int remove_hubble_flow = 0;
    int setdisplacement = 0;
    float unit_mass = 0.0;
    float noise = 0.0;
    char outnamebase[256];
    SDF *csdfp;			/* SDF pointer to control file */
    SDF *sdfp;
    double tpos;			/* time positions are at */
    double tvel;
    double tposlast;
    double dtout, dtvout;
    int do_cosmology = 0;
    int save_first;		/* save first step (for acc testing) */
    double force[NDIM];		/* Things accumulated over all particles */
    double com[NDIM], comv[NDIM];
    double acc2;		/* must be double precision */
    double pe, ke;
    double mtot;
    MPMY_Comm_request req;
    sortresult_t sortedbtab;
    int static_decomp = 0;
    tree_t thetree;
    int massconf, xconf, yconf, zconf;
    int vxconf, vyconf, vzconf;
    int cmassconf;
    int identconf;
    char name[256];
    int read_nfiles, write_nfiles;
    int do_BH, do_DL, do_Bmax, do_Arel, do_n2_ewald;
    float n2_sample_frac;
    int image_freq, x_pixels, y_pixels, log_image;
    float image_size;
    double gnterms;
    int do_periodic, nimage;
    int fix_cube_ewald, fix_cube_ewald_le;
    walkinit_t walk_init_src;
    inherit_t walk_inherit;
#ifdef BODY_HAS_KEY
    pq_keyproto *getkey = (pq_keyproto *)GetKeyFromStruct;
#else
    pq_keyproto *getkey = (pq_keyproto *)GetKey;
#endif
    float R[NDIM];
    int N[NDIM];
    int timeout;
    int set_id;
    int seed;
    ran_state ranstate;
    int memory_tune_MB;
    int identsort_output = 0;
    int do_test = 0;
    int do_profile = 0;
    int do_profile_proc = 0;
    cosmology cosmo;
    int n_outlist = 12;
    double *a_outlist;
    double *t_outlist;
    double default_a_outlist[12] = 
	{0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95,1.0};
    int ic_Nmesh;
    double ic_growthfac;
    mac_s mac = {.type = AREL_MAC, .tol = 1e-3, .tree = &thetree};

    MPMY_Init(&argc, &argv);
    csdfp = startup(argc, argv);
    mxn_s mxn = {.hblock=4096, .min_qsink=64, .min_hsink=64, .min_qsrc=512, .min_hsrc=512};
    mxn_s mxn_cuda = {.hblock=16*1024*1024, .min_qsink=64, .min_hsink=64,
		      .min_qsrc=64, .min_hsrc=64, .do_pQ=1, .do_pH=1};
    if (has_cuda) mxn = mxn_cuda;
    CUDA_Init();
    /* Attempt to get a contiguous chunk of heap */
    SDFgetintOrDefault( csdfp, "memory_tune_MB", &memory_tune_MB, 0);
    if (memory_tune_MB) {
	char *p = Malloc(memory_tune_MB*1024LL*1024LL);
	Free(p);
    }

    SDFgetintOrDefault(csdfp, "timeout", &timeout, 600);
    SDFgetintOrDefault(csdfp, "job_max_wallclock", &job_max_wallclock, MPMY_JobRemaining());
    if (timeout > 0) MPMY_TimeoutSet(timeout);
#ifdef __PARAGON__
    SDFgetintOrDefault(csdfp, "fail_if_slow", &fail_if_slow, 0);
    chk_slow(fail_if_slow);
#endif
    SDFgetstring(csdfp, "datafile", name, sizeof(name));
    SDFgetintOrDefault(csdfp, "do_periodic", &do_periodic, 0);
    SDFgetintOrDefault(csdfp, "nimage", &nimage, 2);
    SDFgetintOrDefault(csdfp, "do_test", &do_test, 0);
    SDFgetintOrDefault(csdfp, "do_profile", &do_profile, 0);
    SDFgetintOrDefault(csdfp, "do_profile_proc", &do_profile_proc, MPMY_Nproc()/3);
    SDFgetintOrDefault(csdfp, "fix_cube_ewald", &fix_cube_ewald, 0);
    SDFgetintOrDefault(csdfp, "fix_cube_ewald_le", &fix_cube_ewald_le, 0);
    SDFgetintOrDefault(csdfp, "cosmology", &do_cosmology, 0);
    SDFgetintOrDefault(csdfp, "set_id", &set_id, 0);
    SDFgetintOrDefault(csdfp, "setpvel", &setpvel, 0);
    SDFgetintOrDefault(csdfp, "remove_hubble_flow", &remove_hubble_flow, 0);
    SDFgetintOrDefault(csdfp, "setdisplacement", &setdisplacement, 0);
    SDFgetintOrDefault(csdfp, "read_nfiles", &read_nfiles, 0);
    SDFgetintOrDefault(csdfp, "write_nfiles", &write_nfiles, 0);
    SDFgetintOrDefault(csdfp, "seed", &seed, 123);
    SDFgetintOrDefault(csdfp, "identsort_output", &identsort_output, 1);
    SDFgetintOrDefault(csdfp, "static_decomp", &static_decomp, 0);
    SDFgetfloatOrDefault(csdfp, "noise", &noise, 0.0);
    SDFgetintOrDefault(csdfp, "Nx",  &N[0], 0);
    SDFgetintOrDefault(csdfp, "Ny",  &N[1], 0);
    SDFgetintOrDefault(csdfp, "Nz",  &N[2], 0);
    ran_init(seed+(MPMY_Procnum()+1), &ranstate);

    if (!do_test) {
	singlPrintf("Reading \"%s\"\n", name);
	if (read_nfiles) MPMY_Nfileio(1);
	/* Use vz for cmass from initial conditions */
	/* Only used when setdisplacement is true */
	sdfp = SDFread64(csdfp, (void *)&btab, &gnobj, &nobj, sizeof(body),
		      "mass", offsetof(body, mass), &massconf,
		      "x", offsetof(body, pos[0]), &xconf,
		      "y", offsetof(body, pos[1]), &yconf,
		      "z", offsetof(body, pos[2]), &zconf,
		      "vx", offsetof(body, vel[0]), &vxconf,
		      "vy", offsetof(body, vel[1]), &vyconf,
		      "vz", offsetof(body, vel[2]), &vzconf,
		      "cmass", offsetof(body, vel[2]), &cmassconf,
		      "ident", offsetof(body, ident), &identconf,
		      NULL);
	if (read_nfiles) MPMY_Nfileio(0);
	Msgf(("Data read, nobj=%d, gnobj=%ld\n", nobj, gnobj));
	Msgf(("Nproc:%d, Procnum: %d, Doc: %d\n",
	       MPMY_Nproc(), MPMY_Procnum(), ilog2(MPMY_Nproc())));
	if (cmassconf && vzconf) {
	    SinglError("You can't have both a 'cmass' and a 'vz' in the data!\n");
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
	if( (identconf == 0) || set_id){
	    SinglWarning("No \"ident\" in file, numbering sequentially\n");
	    FixId(btab, nobj, gnobj);
	}
	/* With relerr MAC acc initialziation, nterms from file is no help */
	FixNterms(btab, nobj);
	SDFgetdoubleOrDefault(sdfp, "Gnewt", &cosmo.Gnewt, 1.0);
	if( SDFhasname("time", sdfp) )
	  SDFgetdoubleOrDefault(sdfp, "time",  &tpos, 0.0);
	else
	  SDFgetdoubleOrDefault(sdfp, "tpos",  &tpos, 0.0);
	SDFgetdoubleOrDefault(sdfp, "tvel",  &tvel, tpos);

	if (do_periodic) {
	    if (SDFhasname("Rx", sdfp)) {
		SDFgetfloatOrDie(sdfp, "Rx",  &R[0]);
		SDFgetfloatOrDie(sdfp, "Ry",  &R[1]);
		SDFgetfloatOrDie(sdfp, "Rz",  &R[2]);
		SDFgetintOrDefault(sdfp, "Nx",  &N[0], N[0]);
		SDFgetintOrDefault(sdfp, "Ny",  &N[1], N[1]);
		SDFgetintOrDefault(sdfp, "Nz",  &N[2], N[2]);
	    } else if (SDFhasname("box_size", sdfp)) {
		SDFgetfloatOrDie(sdfp, "box_size",  &R[0]);
		R[0] /= 2.0;
		R[1] = R[2] = R[0];
		N[0] = N[1] = N[2] = pow(gnobj, (1./3.));
	    } else {
		SDFgetfloatOrDie(sdfp, "R0",  &R[0]);
		R[1] = R[2] = R[0];
		N[0] = N[1] = N[2] = pow(gnobj, (1./3.));
	    }
	} else {
	    SDFgetfloatOrDie(sdfp, "Rx",  &R[0]);
	    SDFgetfloatOrDie(sdfp, "Ry",  &R[1]);
	    SDFgetfloatOrDie(sdfp, "Rz",  &R[2]);
	    /* integer factor to keep error behavior identical to periodic */
	    Error("Need to fix this, using cosmo before it is set\n");
	    VV(sysradius, = (2.0/(1.0 + cosmo.z_at_t(&cosmo, tpos)))*R);
	    VS(rmin, = -sysradius[NDIM-1]);
	    VS(rmax, = sysradius[NDIM-1]);
	    FixRsizeExact(rmin, rmax);
	    Msg_do("new_rmin=(%g, %g, %g) new_rmax=(%g, %g, %g)\n",
		   rmin[0], rmin[1], rmin[2], 
		   rmax[0], rmax[1], rmax[2]);
	}

	if (do_cosmology) {
	    int version;
	    SDFgetintOrDefault(sdfp, "version",  &version, 1);
	    memset(&cosmo, 0, sizeof(cosmo));
	    cosmo.t = tpos;
	    SDFgetdoubleOrDefault(sdfp, "H0",  &cosmo.H0, 0.0511365);
	    if (version == 1) {
		SDFgetdoubleOrDefault(sdfp, "Omega0",  &cosmo.Omega0, 1.0);
		SDFgetdoubleOrDefault(sdfp, "Omega_r",  &cosmo.Omega0_r, 0.0);
		SDFgetdoubleOrDefault(sdfp, "Omega_m",  &cosmo.Omega0_m, cosmo.Omega0-cosmo.Omega0_r);
		SDFgetdoubleOrDefault(sdfp, "Omega_de",  &cosmo.Omega0_fld, 0.0);
		SDFgetdoubleOrDefault(sdfp, "w0",  &cosmo.w0_fld, 0.0);
		SDFgetdoubleOrDefault(sdfp, "wa",  &cosmo.wa_fld, 0.0);
		SDFgetdoubleOrDefault(sdfp, "Lambda_prime",  &cosmo.Omega0_lambda, 0.0);
	    } else if (version == 2) {
		SDFgetdoubleOrDie(sdfp, "h_100",  &cosmo.h_100);
		SDFgetdoubleOrDie(sdfp, "Omega0",  &cosmo.Omega0);
		SDFgetdoubleOrDie(sdfp, "Omega0_r",  &cosmo.Omega0_r);
		SDFgetdoubleOrDie(sdfp, "Omega0_m",  &cosmo.Omega0_m);
		SDFgetdoubleOrDie(sdfp, "Omega0_lambda",  &cosmo.Omega0_lambda);
		SDFgetdoubleOrDie(sdfp, "Omega0_cdm",  &cosmo.Omega0_cdm);
		SDFgetdoubleOrDefault(sdfp, "Omega0_ncdm_tot",  &cosmo.Omega0_ncdm_tot, 0.0);
		SDFgetdoubleOrDie(sdfp, "Omega0_b",  &cosmo.Omega0_b);
		SDFgetdoubleOrDie(sdfp, "Omega0_g",  &cosmo.Omega0_g);
		SDFgetdoubleOrDefault(sdfp, "Omega0_ur",  &cosmo.Omega0_ur, 0.0);
		SDFgetdoubleOrDefault(sdfp, "Omega0_fld",  &cosmo.Omega0_fld, 0.0);
		SDFgetdoubleOrDefault(sdfp, "w0_fld",  &cosmo.w0_fld, 0.0);
		SDFgetdoubleOrDefault(sdfp, "wa_fld",  &cosmo.wa_fld, 0.0);
	    } else Error("Unsupported file version\n");
	    SDFgetdoubleOrDefault(sdfp, "Gnewt", &cosmo.Gnewt, 1.0);
	    /* Now we need to get initial values for cosmo.a */
	    if (SDFhasname("redshift", sdfp)) {
		double Z;
		SDFgetdouble(sdfp, "redshift", &Z);
		cosmo.a = 1.0/(1.0 + Z);
	    } else {
		SinglError("Sorry.  Tell me the redshift in the data file\n");
	    }
	    /* The Zel'dovich 'f' factor is only needed for setting initial
	       velocities.  At this point, we don't know if we will be asked
	       to do setpvel, though, so we read it anyway. */
	    if (SDFhasname("velocity_fac", sdfp)) {
		SDFgetdoubleOrDie(sdfp, "velocity_fac", &cosmo.velfac);
	    } else {
		cosmo.velfac = 1.0;
	    }
	    if (SDFhasname("cosmo_tbl", csdfp)) {
		char cosmo_tbl[256];
		SDFgetstring(csdfp, "cosmo_tbl", cosmo_tbl, sizeof(cosmo_tbl));
		tbl_init(&cosmo, cosmo_tbl);
	    } else {
		cosmo1_init(&cosmo);
	    }
	}
	if (setdisplacement) SDFgetfloatOrDie(sdfp, "unit_mass", &unit_mass);
	SDFgetintOrDefault(sdfp, "iter",  &iter, 0);
	start_iter = iter;
	SDFgetintOrDefault(sdfp, "ic_Nmesh",  &ic_Nmesh, 0);
	if (SDFhasname("ic_growthfac", sdfp)) {
	    SDFgetdoubleOrDie(sdfp, "ic_growthfac",  &ic_growthfac);
	} else if (iter == 0) {
	    ic_growthfac = cosmo.growthfac_at_z(&cosmo, 1.0/cosmo.a-1.0);
	}
	if (sdfp) SDFclose(sdfp);
    } else {
	int cencon;
	int64_t start;
	float rsq;

	singlPrintf("Generating random dataset\n");
	if( SDFgetint64(csdfp, "nobj", &gnobj) )
	  SinglError("Sorry, you've got to have an \"nobj\"\n");
	SDFgetintOrDefault(csdfp, "cencon", &cencon, 0);
	singlPrintf("int seed = %d;\n", seed);
	singlPrintf("int cencon = %d;\n", cencon);

	NobjInitial64(gnobj, MPMY_Nproc(), MPMY_Procnum(), &nobj, &start);
	btab = Calloc(nobj, sizeof(body));
	FixId(btab, nobj, gnobj);
	FixNterms(btab, nobj);
	for (p = &btab[0]; p < &btab[nobj]; p++) {
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

	iter = start_iter = 0;
	if (do_cosmology) {
	    tvel = tpos = 0.05;
	    memset(&cosmo, 0, sizeof(cosmo));
	    cosmo.t = tpos;
	    cosmo.H0 = 0.0511365;
	    cosmo.a = pow(1.5*cosmo.t*cosmo.H0, 2./3.);
	    cosmo.Omega0 = 1.0;
	    cosmo.Omega0_m = cosmo.Omega0-cosmo.Omega0_r;
	    cosmo.Gnewt = GNEWT;
	    R[0] = R[1] = R[2] = 100000.0;
	    for (p = &btab[0]; p < &btab[nobj]; p++) {
		p->mass *= gnobj;
		VV(p->pos, *= cosmo.a*R);
	    }

	} else {
	    tvel = tpos = 0.0;
	    cosmo.Gnewt = 1.0;
	    R[0] = R[1] = R[2] = 1.0;
	}
    }
    if (do_cosmology) {
	singlPrintf("double H0 = %g\n", cosmo.H0);
	singlPrintf("double Omega0 = %g\n", cosmo.Omega0);
	singlPrintf("double Omega0_m = %g\n", cosmo.Omega0_m);
	singlPrintf("double Omega0_r = %g\n", cosmo.Omega0_r);
	if (cosmo.Omega0_fld != 0.0) {
	    singlPrintf("double Omega0_fld = %g\n", cosmo.Omega0_fld);
	    singlPrintf("double w0_fld = %g\n", cosmo.w0_fld);
	    singlPrintf("double wa_fld = %g\n", cosmo.wa_fld);
	}
	singlPrintf("double redshift = %g\n", cosmo.z_at_t(&cosmo, tpos));
	singlPrintf("double Omega0_lambda = %g\n", cosmo.Omega0_lambda);
	singlPrintf("double velfac = %g\n", cosmo.velfac);
    }

    singlPrintf("Maxmem after data read is %d (%d)\n", maxmem(), maxheap());
    if( Msg_test("memleak") ){
	Msg_do("Memory map after data read\n");
	malloc_print();
    }

    /* epsilon_mscale is epsilon for a unit mass particle */
    /* actual epsilon is scaled by mass ^ 1/3 */
    SDFgetfloatOrDie(csdfp, "epsilon_mscale", &eps);
    SDFgetfloatOrDefault(csdfp, "comov_epsilon", &comov_eps, eps);
    SDFgetfloatOrDefault(csdfp, "comov_epsilon_Zmin", &comov_eps_Zmin, 0.0f);
    SDFgetintOrDefault(csdfp, "force_smoothing_type", &force_smoothing_type, 0);
    SDFgetintOrDefault(csdfp, "do_n2_ewald", &do_n2_ewald, 0);
    SDFgetfloatOrDefault(csdfp, "n2_sample_frac", &n2_sample_frac, 1.0);
    SDFgetintOrDefault(csdfp, "do_DL", &do_DL, 0);
    SDFgetintOrDefault(csdfp, "do_BH", &do_BH, 0);
    SDFgetintOrDefault(csdfp, "do_Bmax", &do_Bmax, 0);
    SDFgetintOrDefault(csdfp, "do_Arel", &do_Arel, 0);
    if (do_BH || do_Bmax) 
      SDFgetfloatOrDie(csdfp, "theta", &mac.tol);
    else
      SDFgetfloatOrDie(csdfp, "errtol", &mac.tol);
    SDFgetfloatOrDefault(csdfp, "frac_tol", &mac.rel_tol, 0.0);
    SDFgetfloatOrDefault(csdfp, "frac_tol0", &mac.rel_tol0, 0.0);
    SDFgetfloatOrDefault(csdfp, "CWfac", &CWfac, 0.0);
    SDFgetfloatOrDefault(csdfp, "DLfac", &mac.dlfac, 6.0);
    SDFgetfloatOrDefault(csdfp, "DLmax", &mac.dlmax, 10000.0);
    SDFgetfloatOrDefault(csdfp, "ptol_boost", &mac.ptol_boost, 0.0);
    SDFgetfloatOrDefault(csdfp, "stol_max", &mac.stol_max, 0.0);
    SDFgetintOrDefault(csdfp, "quad_ncut", &mac.p2cut, 7);
    SDFgetintOrDefault(csdfp, "hexa_ncut", &mac.p4cut, 20);
    SDFgetintOrDefault(csdfp, "geometric_center", &mac.geometric_center, 0);
    SDFgetintOrDefault(csdfp, "subtract_background", &mac.subtract_background, 0);
    SDFgetint(csdfp, "mxn_hblock", &mxn.hblock);
    SDFgetint(csdfp, "mxn_min_qsink", &mxn.min_qsink);
    SDFgetint(csdfp, "mxn_min_hsink", &mxn.min_hsink);
    SDFgetint(csdfp, "mxn_min_qsrc", &mxn.min_qsrc);
    SDFgetint(csdfp, "mxn_min_hsrc", &mxn.min_hsrc);
    SDFgetint(csdfp, "mxn_do_pQ", &mxn.do_pQ);
    SDFgetint(csdfp, "mxn_do_pH", &mxn.do_pH);
    SDFgetfloatOrDie(csdfp, "dt", &dt_base); dt = dt_base;
    SDFgetintOrDie(csdfp, "nsteps", &nsteps);
    SDFgetintOrDefault(csdfp, "Ztol", &Ztol, 0);
    SDFgetfloatOrDefault(csdfp, "dt_hiz_tol", &dt_hiz_tol, 0.04);
    SDFgetintOrDefault(csdfp, "save_first", &save_first, 0);
    SDFgetintOrDefault(csdfp, "ntimer_detail", &ntimer_detail, 0);
    SDFgetintOrDefault(csdfp, "light_cone", &light_cone, 0);
    if (light_cone) {
	SDFgetfloatOrDefault(csdfp, "lc_origin_x", &lc_origin[0], 0.0f);
	SDFgetfloatOrDefault(csdfp, "lc_origin_y", &lc_origin[1], 0.0f);
	SDFgetfloatOrDefault(csdfp, "lc_origin_z", &lc_origin[2], 0.0f);
	SDFgetstringOrDefault(csdfp, "lc_name_suffix",lc_name_suffix, 
			      sizeof(lc_name_suffix), "lc000");
    }
    if (do_Bmax) mac.type = BMAX_MAC;
    else if (do_BH) mac.type = BH_MAC;
    else if (do_Arel) mac.type = AREL_MAC;
    else Error("No MAC specified\n");

    if (SDFhasname("a_outlist", csdfp)) {
	n_outlist = SDFarrcnt("a_outlist", csdfp);
	a_outlist = Calloc(n_outlist, sizeof(double));
	SDFrdvecs(csdfp, "a_outlist", 1, a_outlist, n_outlist*sizeof(double), NULL);
    } else {
	a_outlist = default_a_outlist;
    }

    t_outlist = Calloc(n_outlist, sizeof(double));
    for (i = 0; i < n_outlist; i++) {
	t_outlist[i] = cosmo.t_at_z(&cosmo, 1.0/a_outlist[i]-1.0);
    }
    if (!save_first && do_cosmology && tpos*1.0001 > t_outlist[n_outlist-1]
	&& tvel*1.0001 > t_outlist[n_outlist-1]) {
	singlPrintf("Already done.\n");
	exit(0);
    } 
    if (Msg_test("hostname")) {
      char hostname[128];
      gethostname(hostname, sizeof(hostname));
      Msg_do("hostname %s\n", hostname);
    }

    SDFgetstringOrDefault(csdfp, "outfile", outnamebase, sizeof(outnamebase), "output");
    SDFgetintOrDefault(csdfp, "checkpoint_steps_interval", &checkpoint_steps_interval, 2000);
    SDFgetintOrDefault(csdfp, "checkpoint_wallclock_interval", &checkpoint_wallclock_interval, 7200);
    SDFgetintOrDefault(csdfp, "step_wallclock_estimate", &step_wallclock_estimate, 1800);
    SDFgetintOrDefault(csdfp, "output_wallclock_estimate", &output_wallclock_estimate, 1800);
    SDFgetintOrDefault(csdfp, "subsample_steps_interval", &subsample_steps_interval, 2000);
    SDFgetdoubleOrDefault(csdfp, "subsample_fraction", &subsample_fraction, 0.0);
    SDFgetintOrDefault(csdfp, "subsample_random_seed", &subsample_random_seed, 123);
    SDFgetintOrDefault(csdfp, "timer_freq", &timer_freq, 10);
    SDFgetfloatOrDefault(csdfp, "sort_tol", &sort_tol, 0.001);
    SDFgetintOrDefault(csdfp, "image_freq", &image_freq, 0);
    SDFgetfloatOrDefault(csdfp, "image_size", &image_size, 0.0);
    SDFgetintOrDefault(csdfp, "x_pixels", &x_pixels, 512);
    SDFgetintOrDefault(csdfp, "y_pixels", &y_pixels, 512);
    SDFgetintOrDefault(csdfp, "log_image", &log_image, 0);

    if(csdfp) 
	SDFclose(csdfp);

    MPMY_CheckpointSetup(job_max_wallclock, checkpoint_wallclock_interval, 
			 step_wallclock_estimate+output_wallclock_estimate);

    if (do_periodic) {
	EnableTimer(&FixCubeTm, "Fix Cube");
    }
    if (light_cone) {
	EnableTimer(&LightConeTm, "Light Cone");
	/* EnableTimer(&LightConeOpenTm, "Light Cone Open"); */
	/* EnableTimer(&LightConeWriteTm, "Light Cone Write"); */
    }

    singlPrintf("float errtol = %g;\n", mac.tol);
    singlPrintf("float stol_max = %g;\n", mac.stol_max);
    singlPrintf("float ptol_boost = %g;\n", mac.ptol_boost);
    singlPrintf("float dt = %g;\n", dt);
    singlPrintf("double tpos = %g;\n", tpos);
    singlPrintf("double tvel = %g;\n", tvel);
    singlPrintf("int force_smoothing_type = %d;\n", force_smoothing_type);
    singlPrintf("float epsilon_mscale = %g;\n", eps);
    singlPrintf("int start_iter = %d;\n", start_iter);
    singlPrintf("int nsteps = %d;\n", nsteps);
    singlPrintf("int nproc = %d;\n", MPMY_Nproc());
    singlPrintf("int do_Bmax = %d;\n", do_Bmax);
    singlPrintf("int do_BH = %d;\n", do_BH);
    singlPrintf("int do_Arel = %d;\n", do_Arel);
    singlPrintf("int do_DL = %d;\n", do_DL);
    singlPrintf("float CWfac = %.2f;\n", CWfac);
    singlPrintf("float DLfac = %.2f;\n", mac.dlfac);
    singlPrintf("float DLmax = %g;\n", mac.dlmax);
    singlPrintf("float frac_tol = %g;\n", mac.rel_tol);
    singlPrintf("float frac_tol0 = %g;\n", mac.rel_tol0);
    singlPrintf("int quad_ncut = %d;\n", mac.p2cut);
    singlPrintf("int hexa_ncut = %d;\n", mac.p4cut);
    singlPrintf("int geometric_center = %d;\n", mac.geometric_center);
    singlPrintf("int subtract_background = %d;\n", mac.subtract_background);
    singlPrintf("int mxn_hblock = %d;\n", mxn.hblock);
    singlPrintf("int mxn_min_qsink = %d;\n", mxn.min_qsink);
    singlPrintf("int mxn_min_hsink = %d;\n", mxn.min_hsink);
    singlPrintf("int mxn_min_qsrc = %d;\n", mxn.min_qsrc);
    singlPrintf("int mxn_min_hsrc = %d;\n", mxn.min_hsrc);
    singlPrintf("int mxn_do_pQ = %d;\n", mxn.do_pQ);
    singlPrintf("int mxn_do_pH = %d;\n", mxn.do_pH);
    singlPrintf("int body_size = %d;\n", sizeof(body));
    singlPrintf("int outbody_size = %d;\n", sizeof(outbody));
    singlPrintf("int cell_size = %d;\n", sizeof(cell));
    singlPrintf("int quadcell_size = %d;\n", sizeof(quadcell));
    singlPrintf("int hexacell_size = %d;\n", sizeof(hexacell));
    singlPrintf("int cofm_size = %d;\n", sizeof(cofmdata));
    singlPrintf("int tbody_size = %d;\n", TBODYSZ);
    singlPrintf("int hcell_size = %d;\n", sizeof(hcell));
    singlPrintf("int hash_bits = %d;\n", HASH_BITS);
    singlPrintf("int memory_tune_MB = %d;\n", memory_tune_MB);
    singlPrintf("int identsort_output = %d;\n", identsort_output);
    if (static_decomp) singlPrintf("int static_decomp = %d;\n", static_decomp);
    singlPrintf("Checkpoint to %s.nnnn, every %d steps or %d seconds\n", 
		outnamebase, checkpoint_steps_interval, checkpoint_wallclock_interval);
    if (job_max_wallclock > 0) {
	singlPrintf("Final checkpoint in %d seconds (%.2f hours)\n", 
		    job_max_wallclock, job_max_wallclock/3600.0);
    }
    singlPrintf("int timer_freq = %d;\n", timer_freq);
    if (subsample_fraction > 0.0) {
	singlPrintf("int subsample_steps_interval = %d;\n", subsample_steps_interval);
	singlPrintf("int subsample_random_seed = %d;\n", subsample_random_seed);
	singlPrintf("double subsample_fraction = %.7f;\n", subsample_fraction);
    }
    singlPrintf("float sort_tol = %.4f;\n", sort_tol);
    singlPrintf("int do_periodic = %d, nimage is %d;\n", do_periodic, nimage);
    singlPrintf("int fix_cube_ewald_le = %d;\n", fix_cube_ewald_le);

    if (do_cosmology) {
	singlPrintf("int cosmology = %d;\n", do_cosmology);
	singlPrintf("int light_cone = %d;\n", light_cone);
	singlPrintf("int Ztol = %d;\n", Ztol);
	singlPrintf("float comov_epsilon = %f;\n", comov_eps);
	singlPrintf("float comov_epsilon_Zmin = %f;\n", comov_eps_Zmin);
	singlPrintf("float dt_hiz_tol = %f;\n", dt_hiz_tol);
	if (remove_hubble_flow) singlPrintf("int remove_hubble_flow = %d;\n", remove_hubble_flow);
	singlPrintf("int setpvel = %d;\n", setpvel);
	singlPrintf("int setdisplacement = %d;\n", setdisplacement);
	singlPrintf("float noise = %f;\n", noise);
	singlPrintf("float R[0] = %f;\n", R[0]);
	singlPrintf("float R[1] = %f;\n", R[1]);
	singlPrintf("float R[2] = %f;\n", R[2]);
	singlPrintf("Outputs at\n  a      t\n");
	for (i = 0; i < n_outlist; i++) {
	    if (t_outlist[i]*1.0001 > tpos)
		singlPrintf("%.3f %6.3f\n", a_outlist[i], t_outlist[i]);
	}
    }

    singlFflush();
    SanityCheck(btab, nobj, gnobj, &mtot);

    pqsortsetup(&sortedbtab, btab, nobj, sizeof(body), sort_tol, Realloc_f);

    if (do_cosmology && remove_hubble_flow) {
	double H = cosmo.H_at_t(&cosmo, tpos);
	for (p = btab; p < btab+nobj; p++) {
	    VV(p->vel, -= H*p->pos);
	}
	singlPrintf("Hubble flow removed.\n");
    }
    if (do_cosmology && do_periodic && !setdisplacement) {
	/* set_displacement uses union of vel[2], so don't erase it */
	ConvertV(&btab[0].vel[0], sizeof(body), nobj, cosmo.a_at_t(&cosmo, tvel), 0);
    }

    SetupTree(&thetree, NDIM, 
	      sizeof(body), sizeof(cell), sizeof(quadcell), sizeof(hexacell),
	      TBODYSZ, sizeof(cofmdata), (pq_keyproto)getkey, 
	      (static_decomp) ? NULL : (pq_wgtproto)GetCost, CofmFromDaugh, (cellfromcofm_t)CellFromCofm,
	      CellSz);

    if (do_periodic) {
	walk_init_src = (walkinit_t)WalkInitSrcPeriodic;
    } else {
	walk_init_src = (walkinit_t)WalkInitSrc;
    }
    walk_inherit = (inherit_t)InheritSinkNlogN;


    if (do_DL)
	mac.rcrit_func = mac.subtract_background ? (macv_t)DLRcritMACsb : (macv_t)DLRcritMAC;
    else
	mac.rcrit_func = (macv_t)RcritMAC;

    mac.this_tol = mac.tol;

    /* Initial conditions converted from other formats may not be wrapped */
    if (do_periodic) {
      if (do_cosmology)
	  VV(sysradius, = (1.0 / (1.0 + cosmo.z_at_t(&cosmo, tpos)))*R);
      else
	  VV(sysradius,  = R);
      VV(rmin, = -sysradius);
      VV(rmax, = sysradius);
      WrapPeriodic(btab, nobj, rmin, rmax);
    }

    for (iter = start_iter; iter <= start_iter + nsteps; iter++) {
	if (timeout > 0) MPMY_TimeoutReset(timeout);
	/* Reset timers and counters */
	ClearEnabledTimers();
	ClearEnabledCounters();
	StartTimer(&StepTotWC);
	StartTimer(&StepTot);
#ifdef GPERF
        #include <gperftools/profiler.h>
	if (do_profile && (do_profile_proc == -1 || do_profile_proc == MPMY_Procnum())) {
	    char gperffile[256];
	    sprintf(gperffile, "nlnprof_%04d.%04d", iter, MPMY_Procnum());
	    ProfilerStart(gperffile);
	}
#endif

	if (do_cosmology) for (dt = dt_base; dt/tpos > dt_hiz_tol; dt *= 0.5) /* NULL */;

	if (do_periodic) {
	    if (do_cosmology) {
		VV(sysradius,  = (1.0/(1.0 + cosmo.z_at_t(&cosmo, tpos)))*R);
	    } else {
		VV(sysradius, = R);
	    }
	    VS(rmin, = -sysradius[NDIM-1]); /* keep tree cells cubical */
	    VS(rmax, = sysradius[NDIM-1]);
	    FixRsizeExact(rmin, rmax);
	} else {
	    FindBbox(btab, nobj, rmin, rmax);
	    FixRsize(rmin, rmax);
	    VVV(sysradius, = rmax, - rmin);
	    VS(sysradius, *= 0.5);
	}
	Msg_do("rmin=(%g, %g, %g) rmax=(%g, %g, %g)\n",
	      rmin[0], rmin[1], rmin[2], 
	      rmax[0], rmax[1], rmax[2]);


	if (mac.type == AREL_MAC) {
	    /* Perhaps not what you expect for high-aspect volumes */
	    mac.this_tol = mac.tol*mtot*5e5/((1.0+cosmo.z_at_t(&cosmo, tpos))*sysradius[0]*sysradius[0]*sysradius[0]);
	}
	if (do_cosmology && Ztol) {
	    double fac = cosmo.growthfac_at_t(&cosmo, tpos);
	    mac.this_tol *= fac;
	}
	mac.r0 = sysradius[NDIM-1];
	mac.m0 = mtot;
	mac.nx = cbrt(gnobj);
	if (setdisplacement) {
	    mac.this_tol /= 3.0;
	}
	SetupCofm(&mac);

	if (do_cosmology && !do_periodic) {
	    /* integer factor to keep error behavior identical to periodic */
	    VV(sysradius, = (2.0/(1.0 + cosmo.z_at_t(&cosmo, tpos)))*R);
	    VS(rmin, = -sysradius[NDIM-1]);
	    VS(rmax, = sysradius[NDIM-1]);
	    FixRsizeExact(rmin, rmax);
	    Msg_do("new_rmin=(%g, %g, %g) new_rmax=(%g, %g, %g)\n",
		   rmin[0], rmin[1], rmin[2], 
		   rmax[0], rmax[1], rmax[2]);
	}
	/* comoving smoothing */
	if (do_cosmology && (comov_eps > 0.0f) && (cosmo.z_at_t(&cosmo, tpos) > comov_eps_Zmin)) {
	    this_eps = eps*cosmo.a_at_t(&cosmo, tpos)*(1.0+comov_eps_Zmin);
	} else {
	    this_eps = eps;
	}
	if (setdisplacement) this_eps /= 100.0;
	this_eps_scaled = this_eps*pow(btab[0].mass, (float)(1./3.));

	SetupGrav(cosmo.Gnewt, this_eps, gnobj, &mac,
		  btab[0].mass, force_smoothing_type);
#ifdef BODY_HAS_KEY
	FixKeys(btab, nobj, GetKey);
	getkey = (pq_keyproto *)GetKeyFromStruct;
#endif

	if (do_periodic) {
	    float offset[NDIM];
	    VV(offset, = 2.0f*sysradius);
	    SetGravOffset(offset, nimage);
	}

	singlPrintf("BuildTree %d (%d), tol=%g\n", maxmem(), maxheap(), mac.this_tol);

	StartTimer(&BuildTot);
        if (iter == 0 && !static_decomp) {	/* avoid poor load balance on first step */
	    singlPrintf("preliminary BuildTree\n");
            BuildTree(&thetree, &sortedbtab);
            FreeTree(&thetree);
        }
	BuildTree(&thetree, &sortedbtab);
	btab = sortedbtab.data;
	nobj = sortedbtab.nobj;
	StopTimer(&BuildTot);
	singlPrintf("BuildTree done %d (%d)\n", maxmem(), maxheap());
	if( Msg_test("memleak") ){
	    Msg_do("Memory map after BuildTree\n");
	    malloc_print();
	}
	AddCounter(&NbodyCnt, nobj);
	AddCounter(&Ncell_local, ChnAllocCnt(&thetree.cellchn));
	AddCounter(&Nquadcell_local, ChnAllocCnt(&thetree.cell2chn));
	AddCounter(&Nhexacell_local, ChnAllocCnt(&thetree.cell4chn));

	/* Could init in Inherit? */
	for (p = btab; p < btab+nobj; p++) {
	    VS(p->acc, = (float)0.0);
	    p->phi = (float)0.0;
	    p->nterms = 0;
	}

	MPMY_Sync(); /* No sync might cause msg buffer overflow? */
	StartTimer(&FindForcesTm);
	if (do_n2_ewald) {
	    EwaldSetup(2.0*R[2]*cosmo.a_at_t(&cosmo, tpos), cosmo.Gnewt);
	    EwaldForces(btab, nobj, n2_sample_frac, &ranstate);
	} else {
	    StartTimer(&WITm);
	    WalkInitSink(&thetree, btab, nobj, &mxn);
	    WalkInit(&thetree, &thetree, sizeof(Sink), walk_init_src, mac.rcrit_func, 
		     walk_inherit);
	    StopTimer(&WITm);

	    singlPrintf("FindForces\n");
	    StartTimer(&WNTTm);
	    WalkNT(&thetree);
	    StopTimer(&WNTTm);

	    StartTimer(&WTermTm);
	    WalkTerminate();
	    WalkTerminateSink(&thetree, btab, nobj);
	    StopTimer(&WTermTm);
	}
	StopTimer(&FindForcesTm);
	singlPrintf("FindForces done %d (%d)\n", maxmem(), maxheap());

	MPMY_Sync();
	/* This should be the high-water mark for memory use */
	AddCounter(&MemCnt, malloc_used()/(1024*1024));
	AddCounter(&Ncell, ChnAllocCnt(&thetree.cellchn));
	AddCounter(&Nquadcell, ChnAllocCnt(&thetree.cell2chn));
	AddCounter(&Nhexacell, ChnAllocCnt(&thetree.cell4chn));
	AddCounter(&Ntbody, ChnAllocCnt(&thetree.tbodychn));
	AddCounter(&Nhcell, ChnAllocCnt(&thetree.hcellchn));

	if( Msg_test("memleak") ){
	    Msg_do("Memory map before FreeTree\n");
	    malloc_print();
	}
	FreeTree(&thetree);
	singlPrintf("FreeTree done %d (%d)\n", maxmem(), maxheap());
	Msgf(("FreeTree done\n"));
	if( Msg_test("memleak") ){
	    Msg_do("Memory map after FreeTree\n");
	    malloc_print();
	}

	if (do_periodic && !do_n2_ewald) {
	    if (fix_cube_ewald_le) {
		FixCubeEwaldLE(btab, nobj, sysradius, cosmo.Gnewt*mtot, nimage,
			       mac.subtract_background);
	    }
	}

	/* Get force statistics while it is still peculiar acc */
	VS(force, = 0.0);
	acc2 = 0.0;
	for (p = btab; p < btab+nobj; p++) {
	    float sacc2;
	    VV(force, += p->mass*p->acc);
	    sacc2 = Dot(p->acc, p->acc);
	    acc2 += sacc2;
	}
	if (do_cosmology) FixGlobalForce(btab, nobj, &cosmo, tpos, do_periodic);
	
	if (setpvel) {
	    if (setpvel == 2) {
		for (p = btab; p < btab+nobj; p++) {
		    VS(p->vel, = 0.0f);
		}
		singlPrintf("Velocities erased.\n");
	    } else {
		set_vels(btab, nobj, dt, &cosmo, &tpos, &tvel, do_periodic);
		singlPrintf("Velocities adjusted to linear approximation.\n");
	    }
	    setpvel = 0;
	    if (do_cosmology && do_periodic) {
		ConvertV(&btab[0].vel[0], sizeof(body), nobj, cosmo.a_at_t(&cosmo, tvel), 0);
	    }
	} else if (setdisplacement) {
	    setdisplacement = 0;
	    set_displacement(btab, nobj, dt, &cosmo, &tpos, &tvel, noise, do_periodic, &ranstate);
	    singlPrintf("Displacements and velocities adjusted to linear approximation.\n");
	    if (do_cosmology && do_periodic) {
		ConvertV(&btab[0].vel[0], sizeof(body), nobj, cosmo.a_at_t(&cosmo, tvel), 0);
	    }
	    if (do_periodic) {
		if (do_cosmology) {
		    VV(sysradius, = (1.0/(1.0 + cosmo.z_at_t(&cosmo, tpos)))*R);
		} else {
		    VV(sysradius, = R);
		}
		VV(rmin, = -sysradius);
		VV(rmax, = sysradius);
		WrapPeriodic(btab, nobj, rmin, rmax);
	    }
	}

	do_output = do_checkpoint = 0;
	dtout = dt; /* time for checkpoint, reset below for do_output */
	dtvout = tpos + 0.5 * dt - tvel;
	for (i = 0; i < n_outlist; i++) {
	    if ((tpos < 1.0001*t_outlist[i]) && (tpos + dt >= t_outlist[i])) {
		if (!(first_step && tpos == tvel)) { /* Don't rewrite output */
		    do_output = 1;
		    dtout = t_outlist[i] - tpos;
		    dtvout = t_outlist[i] - tvel;
		    break;
		}
	    }
	}
	if ((first_step && save_first) || ForceOutput()) {
	    do_output = 1;
	    dtout = 0.0;
	    dtvout = tpos - tvel;
	    save_first = 0;
	}
	steps_to_next_output = nsteps; /* init to large number */
	for (i = 0; i < n_outlist; i++) {
	    int steps = (t_outlist[i] - tpos)/dt;
	    if (steps > 0 && steps < steps_to_next_output)
		steps_to_next_output = steps;
	}
	if (MPMY_CheckpointDue(steps_to_next_output*step_wallclock_estimate))
	    do_checkpoint = 1;
	if ((iter-start_iter+1+checkpoint_steps_interval) % checkpoint_steps_interval == 0) do_checkpoint = 1;
	if (ForceCheckpoint()) do_checkpoint = 1;

	if (0 == iter % subsample_steps_interval && iter && subsample_fraction > 0.0) {
	    output(outnamebase, gnobj, nobj, btab, iter, dtout, dtvout,
		   &cosmo, tpos, tvel, do_cosmology, do_periodic, 
		   eps, this_eps_scaled, force_smoothing_type, mac.this_tol, 
		   mac.rel_tol, mac.rel_tol0,
		   R, N, write_nfiles, &ke, &pe, 1, identsort_output,
		   ic_Nmesh, ic_growthfac, subsample_fraction, subsample_random_seed);
	}

	if (do_output || do_checkpoint) {
	    output(outnamebase, gnobj, nobj, btab, iter, dtout, dtvout,
		   &cosmo, tpos, tvel, do_cosmology, do_periodic, 
		   eps, this_eps_scaled, force_smoothing_type, mac.this_tol, 
		   mac.rel_tol, mac.rel_tol0,
		   R, N, write_nfiles, &ke, &pe, do_output, identsort_output,
		   ic_Nmesh, ic_growthfac, 0.0, 0);
	    if (!do_output) MPMY_CheckpointFinished();
	}

	if (light_cone && do_output) {
	    /* We need these in case of a restart from a non-checkpoint file */
	    VS(lc_origin, = 0.0f);
	    WriteLightCone(btab, nobj, dtout, dtvout, &cosmo, tpos, tvel, 
			   outnamebase, "lc000o", iter, lc_origin, R);
	    VS(lc_origin, = -0.9*R[0]);
	    WriteLightCone(btab, nobj, dtout, dtvout, &cosmo, tpos, tvel, 
			   outnamebase, "lc111o", iter, lc_origin, R);
	}

	if (light_cone) {
	    VS(lc_origin, = 0.0f);
	    WriteLightCone(btab, nobj, dt, tpos+0.5*dt-tvel, &cosmo, tpos, tvel, 
			   outnamebase, "lc000", iter, lc_origin, R);
	    VS(lc_origin, = -0.9*R[0]);
	    WriteLightCone(btab, nobj, dt, tpos+0.5*dt-tvel, &cosmo, tpos, tvel, 
			   outnamebase, "lc111", iter, lc_origin, R);
	}
	
	if ((do_cosmology && tpos + dt >= t_outlist[n_outlist-1]) || (do_checkpoint && ForceStop()) || MPMY_JobDone()) {
	    singlPrintf("Stopping.\n");
	    break;
	} 

	Msgf(("integrating positions and velocities\n"));
	tposlast = tpos;
	if (do_cosmology && do_periodic) {
	    CosmoIntegrate(&btab[0].mass, &btab[0].pos[0], &btab[0].vel[0],
			   &btab[0].acc[0], &btab[0].phi, sizeof(body),
			   &btab[0].pos[0], &btab[0].vel[0], sizeof(body),
			   nobj, dt, tpos+0.5*dt-tvel, &cosmo, &tpos, &tvel, &ke, &pe);
	} else {
	    Integrate(&btab[0].mass, &btab[0].pos[0], &btab[0].vel[0],
		      &btab[0].acc[0], &btab[0].phi, sizeof(body),
		      &btab[0].pos[0], &btab[0].vel[0], sizeof(body),
		      nobj, dt, tpos+0.5*dt-tvel, &tpos, &tvel, &ke, &pe);
	}
	if (do_periodic) {
	    if (do_cosmology) {
		VV(sysradius, = (1.0 / (1.0 + cosmo.z_at_t(&cosmo, tpos)))*R);
	    } else {
		VV(sysradius, = R);
	    }
	    VV(rmin, = -sysradius);
	    VV(rmax, = sysradius);
	    WrapPeriodic(btab, nobj, rmin, rmax);
	}

	VS(com, = 0.0);
	VS(comv, = 0.0);
	mtot = 0.0;
	gnterms = 0.0;
	for (p = btab; p < btab+nobj; p++) {
	    VV(com, += p->mass*p->pos);
	    VV(comv, += p->mass*p->vel);
	    mtot += p->mass;
	    gnterms += p->nterms;
	    if (p->nterms <= 0) SeriousWarning("nterms is %f\n", p->nterms);
	}
	AddCounter(&KNtermsCnt, gnterms/1000.0);
	if (ReadCounter(&BSInt))
	    AddCounter(&Scycles, 10.0*CPU.Hz*ReadTimer(&GravSTm)/ReadCounter(&BSInt));
	if (ReadCounter(&BCInt))
	    AddCounter(&Mcycles, 10.0*CPU.Hz*ReadTimer(&GravMTm)/ReadCounter(&BCInt));
	if (ReadCounter(&BC2Int))
	    AddCounter(&Qcycles, 10.0*CPU.Hz*ReadTimer(&GravQTm)/ReadCounter(&BC2Int));
	if (ReadCounter(&FBC2Int))
	    AddCounter(&FQcycles, 10.0*CPU.Hz*ReadTimer(&GravQFTm)/ReadCounter(&FBC2Int));
	if (ReadCounter(&BC4Int))
	    AddCounter(&Hcycles, 10.0*CPU.Hz*ReadTimer(&GravHTm)/ReadCounter(&BC4Int));
	if (ReadCounter(&FBC4Int))
	    AddCounter(&FHcycles, 10.0*CPU.Hz*ReadTimer(&GravHFTm)/ReadCounter(&FBC4Int));
	if (ReadCounter(&Qcycles) > 230) {
	    char hostname[128];
	    gethostname(hostname, sizeof(hostname));
	    SeriousWarning("Proc %d %s is slow, Qcycles is %ld\n", 
			   MPMY_Procnum(), hostname, ReadCounter(&Qcycles));
	}
	   
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
	MPMY_ICombine_Wait(req);
	Msgf(("done MPMY_combine\n"));
	StopTimer(&StepTot);
	StopTimer(&StepTotWC);

	if (do_cosmology) 
	    singlPrintf("\ntpos = %g, znow = %g, iter = %d, size = %g, eps = %g\n", 
			tposlast, cosmo.z_at_t(&cosmo, tposlast),
			iter, (1.0 / (1.0 + cosmo.z_at_t(&cosmo, tposlast)))*R[0], this_eps_scaled);
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
	AddCounter(&HeapCnt_, malloc_heapsz()/(1024*1024));
		    
	if (timer_freq && iter%timer_freq == 0) {
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
		    float mac_tm;
		    float imbal_tm;
		    float per_tm;
		    int nterms;
		    int nbody;
		} perf, *gp = NULL;

		perf.node = MPMY_Procnum();
		perf.grav_tm = ReadTimer(&GravTm);
		perf.mac_tm = ReadTimer(&MACTm);
		perf.imbal_tm = ReadTimer(&WTermTm);
		perf.nterms = ReadCounter(&KNtermsCnt);
		perf.nbody = nobj;

		if (MPMY_Procnum() == 0) {
		    gp = Malloc(MPMY_Nproc() * sizeof(perf));
		    MPMY_Gather(&perf, sizeof(perf), MPMY_CHAR, gp, 0);
		    for (i = 0; i < MPMY_Nproc(); i++) 
		      singlPrintf("%3d %8.2f %8.2f %8.2f %8.2f %10d %6d\n",
				  gp[i].node, gp[i].grav_tm, gp[i].mac_tm, 
				  gp[i].imbal_tm, gp[i].per_tm, 
				  gp[i].nterms, gp[i].nbody);
		    Free(gp);
		} else {
		    MPMY_Gather(&perf, sizeof(perf), MPMY_CHAR, gp, 0);
		}
	    }
	    if ((GravTm.mean > 10.0) && (WTermTm.mean > 0.5*GravTm.mean) && (ReadTimer(&WTermTm) < 5.0)) {
		char hostname[128];
		gethostname(hostname, sizeof(hostname));
		SeriousWarning("High load imbalance.  Check proc %d at %s\n", 
			       MPMY_Procnum(), hostname);
	    }
	} else {
	    OutputTimer(&StepTot, singlPrintf);
	    OutputTimer(&StepTotWC, singlPrintf);
	}
	singlFflush();

	/* This can greatly improve the load balance */
	if (CWfac != 0.0) {	/* CWfac = 1 seems to work well for do_periodic */
	    gnterms /= gnobj;
	    singlPrintf("Avg nterms = %.0f\n", gnterms);
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
#ifdef GPERF
	if (do_profile && (do_profile_proc == -1 || do_profile_proc == MPMY_Procnum())) ProfilerStop();
#endif
    }
    singlPrintf("Bye!\n");
    MPMY_Finalize();
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
    char hostname[128];

    if (argc > 1)
 	strncpy(cfile, argv[1], sizeof(cfile));
    else
 	Getsparam("control file", cfile);
    if ((csdfp = SDFopen(NULL, cfile)) == NULL) {
 	SinglError("Sorry, couldn't SDFopen %s\n%s\n",
	      cfile, SDFerrstring);
    }

    singlPrintf("Welcome to the tree19 %sN-body integrator running on %d procs\n",
#ifdef HEXA
		"hexadecapole ",
#else
#ifdef QUAD
		"quadrupole ",
#else
		"monopole ",
#endif
#endif
		MPMY_Nproc());
    singlPrintf("Version %s\n", Version);
    singlPrintf("Compiled %s %s\n", Compiled_date, Compiled_time);
    singlPrintf("Compiler %s\n", Compiler);
    singlPrintf("Arch %s\n", Arch);
    CPU.Hz = clockspeedHz(&CPU.ncores, &CPU.cpuname);
    singlPrintf("CPU %s, %d cores, %.3f GHz\n", CPU.cpuname, CPU.ncores, 1e-9*CPU.Hz);

#ifdef SAVE_ACC
    singlPrintf("Saving phi and acc\n");
#endif
#ifdef BODY_HAS_KEY
    singlPrintf("Bodies cache their key\n");
#endif
    gethostname(hostname, sizeof(hostname));
    singlPrintf("Node 0 is %s\n", hostname);

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

    if( Msg_test("Malloc.c") ){
	malloc_debug(2);
	Msg_do("Malloc_debug(2), expect slow mallocs\n");
    }else{
	malloc_debug(1);
    }

    EnableWCTimer(&StepTotWC, "Step (Wall)");
    EnableCPUTimer(&StepTot, "Step (CPU)");
    EnableTimer(&BuildTot, "Build Total");
    EnableTimer(&PQSortTm, "PQSort");
    EnableTimer(&PQSortCommTm, "PQSortComm");
    EnableTimer(&PQSortWaitTm, "PQSortWait");
    EnableTimer(&PQSortAtoaTm, "Alltoall");
    EnableTimer(&PQSortAtoavTm, "Alltoallv");
    EnableTimer(&DecompTm, "Decomp");
    EnableTimer(&DecompCommTm, "DecompComm");
    EnableTimer(&DecompWaitTm, "DecompWait");
    EnableTimer(&SortTm, "Sort");
    EnableTimer(&MakeTreeTm, "MakeTree");
    EnableTimer(&SharedCellsTm, "SharedCells");
    EnableTimer(&SharedCellsCommTm, "SharedComm");
    EnableTimer(&SharedCellsWaitTm, "SharedWait");
    EnableTimer(&FindForcesTm, "Force Eval");
    EnableCPUTimer(&GravSTm, "Smth Time");
    EnableCPUTimer(&GravMTm, "Mono Time");
    EnableCPUTimer(&GravQTm, "Quad Time");
    EnableCPUTimer(&GravQFTm, "QuadF Time");
    EnableCPUTimer(&GravHTm, "Hexa Time");
    EnableCPUTimer(&GravHFTm, "HexaF Time");
    EnableCPUTimer(&GravTm, "Grav Time");
    EnableCPUTimer(&MACTm, "MAC Time");
    EnableTimer(&CUDAWtTm, "CUDA Wait");
    EnableTimer(&WalkDeferTm, "Walk Defer");
    EnableTimer(&WTermTm, "WalkTerm");
    EnableTimer(&WNTTm, "WalkNT");
    EnableCounter(&PQSortSends, "PQSort Sends");
    EnableCounter(&PQSortRecvs, "PQSort Recvs");
    EnableCounter(&PQSortMaxn, "PQSort Maxn");
    EnableCounter(&BSInt, "Body-smth");
    EnableCounter(&BCInt, "Body-mono");
    EnableCounter(&BC2Int, "Body-quad");
    EnableCounter(&FBC2Int, "Body-quadF");
    if (has_cuda) EnableCounter(&FBC2FInt, "Body-quadFF");
    EnableCounter(&BC4Int, "Body-hexa");
    EnableCounter(&FBC4Int, "Body-hexaF");
    if (has_cuda) EnableCounter(&FBC4FInt, "Body-hexaFF");
    EnableCounter(&MACcnt, "MAC");
    EnableCounter(&BBMACcnt, "BB MAC");
    EnableCounter(&CEmpty, "Cube Empty");
    EnableCounter(&MCAnti, "Cube Anti");
    EnableCounter(&MCCorr, "Cube Corr");
    EnableCounter(&BSMax, "Max smth");
    EnableCounter(&Scycles, "Smth cycles");
    EnableCounter(&Mcycles, "Mono cycles");
    EnableCounter(&Qcycles, "Quad cycles");
    EnableCounter(&FQcycles, "QuadF cycles");
    EnableCounter(&Hcycles, "Hexa cycles");
    EnableCounter(&FHcycles, "HexaF cycles");
    EnableCounter(&KNtermsCnt, "KNterms");
    EnableCounter(&NbodyCnt, "Bodies");
    EnableCounter(&Ncell_local, "Local Cells");
    EnableCounter(&Nquadcell_local, "Local Quad");
    EnableCounter(&Nhexacell_local, "Local Hexa");
    EnableCounter(&Ncell, "Cells");
    EnableCounter(&Nquadcell, "Quad Cells");
    EnableCounter(&Nhexacell, "Hexa Cells");
    EnableCounter(&Ntbody, "TBodies");
    EnableCounter(&Nhcell, "Hcells");
    EnableCounter(&SharedCnt, "Shared Cells");
    EnableCounter(&RequestCnt, "Requests");
    EnableCounter(&DeferCnt, "Deferred");
    EnableCounter(&MemCnt, "Memory (MB)");
    EnableCounter(&HeapCnt_, "Heap (MB)");
    return csdfp;
}

static void SanityCheck(bodyptr btab, int nobj, int64_t gnobj, double *mtotp){
    double mtot;
    bodyptr p;
    int64_t sumnobj;
    MPMY_Comm_request req;

    mtot = 0.0;
    for (p = btab; p < btab+nobj; p++) {
	mtot += p->mass;
	/* We really could use more checks here! */
    }
    sumnobj = nobj;
    MPMY_ICombine_Init(&req);
    MPMY_ICombine(&mtot, &mtot, 1, MPMY_DOUBLE, MPMY_SUM, req);
    MPMY_ICombine(&sumnobj, &sumnobj, 1, MPMY_INT64, MPMY_SUM, req);
    MPMY_ICombine_Wait(req);
    assert(sumnobj == gnobj);
    Msgf(("Particle 0 (%ld), %g, %g %g %g, %g %g %g\n",
	  btab->ident, btab->mass, 
	  btab->pos[0], btab->pos[1], btab->pos[2],
	  btab->vel[0], btab->vel[1], btab->vel[2]));
    Msgf(("Particle %d (%ld), %g, %g %g %g, %g %g %g\n", nobj-1,
	  btab[nobj-1].ident, btab[nobj-1].mass, 
	  btab[nobj-1].pos[0], btab[nobj-1].pos[1], btab[nobj-1].pos[2],
	  btab[nobj-1].vel[0], btab[nobj-1].vel[1], btab[nobj-1].vel[2]));
    singlPrintf("gnobj = %ld\n", gnobj);
    singlPrintf("mtot = %lf\n", mtot);
    *mtotp = mtot;
}

static void
FixGlobalForce(body *xptr, int n, cosmology *c, float real_time, int do_periodic){
    /* Make whatever corrections are necessary to the acceleration, etc.
       based on values of GNewt, Lambda, etc., etc. */
    float acc_back;
    body *p;

    if (do_periodic) return;	/* now handled by CosmoIntegrate */

    acc_back = c->Omega0_lambda*c->H0*c->H0;
    while(n--){
	p = xptr++;
	VV(p->acc, += acc_back*p->pos);
    }
}

/* This erases the velocities, and sets them */
/* according to linear theory */

static void
set_vels(body *p, int n, float dt, cosmology *c, double *tpos, double *tvel, int do_periodic)
{
    float peculiar_acc[NDIM], proper_acc[NDIM];
    body *end = p + n;
    double acc_back;
    double vel_fac, pos_fac;
    double a;
    double asum1, asum2;

    a = c->a_at_t(c, *tpos);

    /* mean was subtraced in FixCube for do_periodic */
    /* Lambda was added in FixGlobalForce for !do_periodic */
    acc_back = (0.5*c->Omega0/(a*a*a) - c->Omega0_lambda)*c->H0*c->H0;

    /* density */
    pos_fac = (a*a*a)/(1.5 * c->Omega0_m * c->H0 * c->H0);
    vel_fac = pos_fac * c->velfac * c->H;

    *tvel = *tpos;

    singlPrintf("set_vels: tpos = %g, tvel = %g\n", *tpos, *tvel);
    singlPrintf("set_vels: velfac = %g, vel_fac/t is %.4f * t, H is %f\n", 
		c->velfac, vel_fac/(*tpos), c->H);

    asum1 = 0.;
    asum2 = 0.;
    for (; p < end; p++) {
	if (do_periodic) {
	    VVV(proper_acc, = p->acc, + acc_back*p->pos);
	    VV(peculiar_acc, = p->acc);
	} else {
	    VV(proper_acc, = p->acc);
	    VVV(peculiar_acc, = p->acc, + acc_back*p->pos);
	}
	asum1 += Dot(proper_acc, p->pos); /* diagnostic */
	asum2 += Dot(peculiar_acc, p->pos); /* diagnostic */
	if (do_periodic) {
	    VV(p->vel, = vel_fac*p->acc);
	} else {
	    VVV(p->vel, = vel_fac*p->acc, + c->H*p->pos); /* pec. vel + Hubble flow */
	}
    }
    singlPrintf("Mean(proper acc dot position) = %g\n", asum1/n);
    singlPrintf("Mean(peculiar acc dot position) = %g\n", asum2/n);
}

static void
set_displacement(body *p, int n, double dt, cosmology *c, double *tpos, double *tvel, 
		 float noise, int do_periodic, ran_state *ranstate)
{
    float peculiar_acc[NDIM], proper_acc[NDIM];
    body *end = p + n;
    double acc_back;
    double vel_fac, pos_fac;
    double a;
    float mass;
    double asum1, asum2;
    
    a = c->a_at_t(c, *tpos);
 
    /* mean was subtraced in FixCube for do_periodic */
    /* Lambda was added in FixGlobalForce for !do_periodic */
    acc_back = (0.5*c->Omega0/(a*a*a) - c->Omega0_lambda)*c->H0*c->H0;

    /* density */
    pos_fac = (a*a*a)/(1.5 * c->Omega0_m * c->H0 * c->H0);
    vel_fac = pos_fac * c->velfac * c->H;

    *tvel = *tpos;

    singlPrintf("set_displacement: tpos = %g, tvel = %g\n", *tpos, *tvel);
    singlPrintf("set_displacement: velfac = %g, vel_fac/t = %.4f * t, H = %f\n", 
		c->velfac, vel_fac/(*tpos), c->H);

    asum1 = 0.;
    asum2 = 0.;
    for (; p < end; p++) {
	if (do_periodic) {
	    VVV(proper_acc, = p->acc, + acc_back*p->pos);
	    VV(peculiar_acc, = p->acc);
	} else {
	    VV(proper_acc, = p->acc);
	    VVV(peculiar_acc, = p->acc, + acc_back*p->pos);
	}
	asum1 += Dot(proper_acc, p->pos); /* diagnostic */
	asum2 += Dot(peculiar_acc, p->pos); /* diagnostic */
	mass = p->vel[2];	/* we stuck cmass into vel[2]. could be union? */
	if (p->mass/mass < 0.1 || p->mass/mass > 10.0) {
	    SeriousWarning("High mass ratio in set_displacement %f %f\n",
			   p->mass, mass);
	}
	p->mass = mass;
	VVV(p->pos, = pos_fac*peculiar_acc, + p->pos);
	if (do_periodic) {
	    VV(p->vel, = vel_fac*peculiar_acc);
	} else {
	    VVV(p->vel, = vel_fac*peculiar_acc, + c->H*p->pos); /* pec. vel + Hubble flow */
	}
	if (noise != 0.0) {
	  float vec[NDIM];
	  sphere_rand(ranstate, NDIM, vec);
	  VV(p->pos, += noise*a*vec);
	}
    }
    singlPrintf("Mean(proper acc dot position) = %lg\n", asum1/n);
    singlPrintf("Mean(peculiar acc dot position) = %lg\n", asum2/n);
}


static void
WrapPeriodic(body *bp, int n, float *rmin, float *rmax)
{
    body *b;
    int fluxp[NDIM] = {0, 0, 0};
    int fluxm[NDIM] = {0, 0, 0};
    float sz[NDIM];
    float center[NDIM];
    float exact_rmin[NDIM]; /* match IntPos Key logic exactly */
    float exact_rmax[NDIM];

    VVV(sz, = rmax, - rmin);
    VVVS(center, = LPAREN rmax, + rmin, RPAREN*0.5f);
    VVV(exact_rmin, = -0.5f*sz, + center);
    VVV(exact_rmax, = exact_rmin, + sz);

    for(b=bp; b<&bp[n]; b++) {
	VVVS(if LPAREN b->pos, >= exact_rmax, RPAREN fluxp, += 1);
	VVVV(if LPAREN b->pos, >= exact_rmax, RPAREN b->pos, -= sz);
	VVVS(if LPAREN b->pos, < exact_rmin, RPAREN fluxm, += 1);
	VVVV(if LPAREN b->pos, < exact_rmin, RPAREN b->pos, += sz);
    }
    MPMY_Combine(fluxp, fluxp, NDIM, MPMY_INT, MPMY_SUM);
    MPMY_Combine(fluxm, fluxm, NDIM, MPMY_INT, MPMY_SUM);
    singlPrintf("edge flux %d %d %d plus %d %d %d minus %d %d %d\n", 
		fluxp[0]-fluxm[0],
		fluxp[1]-fluxm[1],
		fluxp[2]-fluxm[2],
		fluxp[0], fluxp[1], fluxp[2], 
		fluxm[0], fluxm[1], fluxm[2]);
}

static float Q[1716];

static void 
FixCubeEwaldLE(body *btab, int nobj, const float *l, float gm, int nimage,
	       int subtract_background)

{
    float r[NDIM];
    body *p;
    double f[NDIM];
    double phi;
    float grho;
    int i[NDIM], ii;
    double phicorr = 0.0;

    StartTimer(&FixCubeTm);
    grho = gm/(8.*l[0]*l[1]*l[2]); /* G rho */

    for (i[0] = -nimage; i[0] <= nimage; i[0]++) {
	for (i[1] = -nimage; i[1] <= nimage; i[1]++) {
	    for (i[2] = -nimage; i[2] <= nimage; i[2]++) {
		ii = Dot(i, i);
		if (ii) phicorr += 1./sqrt((double)ii);
	    }
	}
    }

    calculate_cartesian_moments(btab, nobj, l[0]*2.0, Q, subtract_background);

    Msgf(("grho %g 2*L*grho %g\n", grho, l[0]*2.0f*grho));

    if (!subtract_background) {
	for (p = btab; p < btab+nobj; p++) {
	    VV(r, = p->pos);
	    VS(r, /= l[0]*2.0f);
	    ewald_background(r, p->mass, phicorr, f, &phi);
	    VS(f, *= l[0]*2.0f*grho);
	    VV(p->acc, += f);
	    phi *= l[0]*l[0]*4.0f*grho;
	    p->phi += phi;
	}
    }

    for (p = btab; p < btab+nobj; p++) {
	VV(r, = p->pos);
	VS(r, /= l[0]*2.0f);
	ewald_le(r, f, &phi, Q, nimage);
	VS(f, *= l[0]*2.0f*grho);
	phi *= l[0]*l[0]*4.0f*grho;
	VV(p->acc, += f);
	p->phi += phi;
    }
    StopTimer(&FixCubeTm);
}

int 
maxheap(void)
{
    int memused;
    size_t s = malloc_heapsz();
   
    memused = s/(1024LL*1024LL);
    MPMY_Combine(&memused, &memused, 1, MPMY_INT, MPMY_MAX);
    return memused;
}

int 
maxmem(void)
{
    int memused;
    size_t s = malloc_used();
   
    memused = s/(1024LL*1024LL);
    MPMY_Combine(&memused, &memused, 1, MPMY_INT, MPMY_MAX);
    return memused;
}

static void
WriteLightCone(body *xptr, const int n, const double dt, const double dtv, 
	       cosmology *c, double tpos, double tvel,
	       char *name, char *tag, int iter, float *lc_origin, float *rmax)
{
    body *end = xptr + n;
    double tpos0, tpos1;
    double a0, a1;
    double hubble0, hubble1;
    double dt_kick, dt_kick0, dt_kick1, dt_drift1;
    double alpha;
    double r0, r1, r02, r12;
    float pos[NDIM];
    float vel[NDIM];
    float p0[NDIM], p1[NDIM], v0[NDIM], v1[NDIM];
    Stk outstk;
    int nout, gnout;
    int outsize = 36;

    /* Write particles between tpos and tpos+dt */
    StartTimer(&LightConeTm);
    StkInitEz(&outstk);
    a0 = c->a_at_t(c, tpos);
    hubble0 = c->H_at_z(c, 1.0/a0-1.0);
    r0 = c->conformal_distance_at_z(c, 1.0/a0-1.0);
    tpos0 = tpos;
    tpos1 = tpos+dt;
    dt_kick = a0*a0*c->kick_t0_t1(c, tvel, tvel + dtv);
    dt_kick0 = a0*a0*c->kick_t0_t1(c, tvel, tpos0);
    dt_kick1 = a0*a0*c->kick_t0_t1(c, tvel, tpos1);
    dt_drift1 = c->drift_t0_t1(c, tpos, tpos1);
    a1 = c->a_at_t(c, tpos1);
    hubble1 = c->H_at_z(c, 1.0/a1-1.0);
    r1 = c->conformal_distance_at_z(c, 1.0/a1-1.0);
    if (r1 < 0.0) r1 = 0.0;

    for (; xptr < end; xptr++) {
	VV(p0, = (1.0/a0) * xptr->pos); /* to comoving */
	VV(p0, -= lc_origin);
	VV(p1, = p0);
	VVV(vel, = xptr->vel, + dt_kick * xptr->acc);
	VV(p1, += dt_drift1 * vel);
	r02 = Dot(p0, p0);
	r12 = Dot(p1, p1);
	if (r02 <= r0*r0 && r12 > r1*r1) {
	    VVV(v0, = xptr->vel, + dt_kick0 * xptr->acc);
	    VVV(v1, = xptr->vel, + dt_kick1 * xptr->acc);
	    /* to physical vel */
	    VS(v0, /= a0);
	    VS(v1, /= a1);
	    /* add hubble flow */
	    VV(v0, += hubble0 * a0 * p0);
	    VV(v1, += hubble1 * a1 * p1);
	    
	    r02 = sqrt(r02);
	    r12 = sqrt(r12);
	    alpha = (r12-r1)/(r0-r1+r12-r02);
	    VVV(pos, = alpha*p0, + (1.0-alpha)*p1);
	    VVV(vel, = alpha*v0, + (1.0-alpha)*v1);
	    assert(outsize == 7*sizeof(float) + sizeof(int64_t));
	    StkPushData(&outstk, &(xptr->mass), sizeof(float));
	    StkPushData(&outstk, pos , NDIM*sizeof(float));
	    StkPushData(&outstk, vel , NDIM*sizeof(float));
	    StkPushData(&outstk, &(xptr->ident), sizeof(int64_t));
	}
    }
    char outname[256];
    nout = StkSz(&outstk)/outsize;
    MPMY_Combine(&nout, &gnout, 1, MPMY_INT, MPMY_SUM); /* could be eliminated */
    sprintf(outname, "lc/%s_%s.%04d.%04d", name, tag, iter, MPMY_Procnum()/PROCS_PER_NODE);
    if (nout > 0) {
	StartTimer(&LightConeOpenTm);
	FILE *outfp = fopen(outname, "a");
	StopTimer(&LightConeOpenTm);
	StartTimer(&LightConeWriteTm);
	fwrite(StkBase(&outstk), nout, outsize, outfp);
	StopTimer(&LightConeWriteTm);
	fclose(outfp);
    }
    if (gnout > 0)
	singlPrintf("fwrite of %d to %s at time = %6.3f r0 = %8.2f z = %6.4f origin %g %g %g\n",
		    gnout, outname, tpos, r0, c->z_at_t(c, tpos), 
		    lc_origin[0], lc_origin[1], lc_origin[2]);
    StkTerminate(&outstk);
    StopTimer(&LightConeTm);
}


