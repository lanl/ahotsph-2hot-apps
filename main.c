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
#include "nrutil.h"
#include "units.h"
#include "cool.h"

#define MAXCOEF 16

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
static void SPHOutputA(SPHbody *btab, int nobj, const char *outname, int iter);
static void WindOutput(SPHbody *btab, int nobj, windbody *windbtab, 
		       int windnobj, const char *outname, int iter);
static void ShortWindOutput(SPHbody *btab, int nobj, windbody *windbtab, 
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
int make_spec_names(char ***chararr, char spec, int num);

/* In shrink.c */
/*  void ShrinkBtab(SPHbody **SPHbtabp, body *btabp, int *nobj, float r_limit); */
/* void ShrinkBtab2(SPHbody **SPHbtabp, int *nobj, float r_limit); */
void AdjustBtab(SPHbody **SPHbtabp, int *nobj, int gnobj, windbody *windbtab, 
		int windnobj, int windpartpershell, float r_limit, float dt, 
		int iter, float tpos, int *added_particles);
void AdjustBtab2(SPHbody **SPHbtabp, int *nobj, int gnobj, windbody *windbtab, 
		int windnobj, float r_limit, float dt, int iter, float tpos,
		int *added_particles, float *newmass);
void AdjustBtab3(SPHbody **SPHbtabp, int *nobj, int gnobj, float r_limit,
		 float r_outer);
void AdjustBtab4(SPHbody ** SPHbtabp, int *nobj, bndry_t b, float *newmass,
                 float *newr, float *newp, float *newl, float newt, float tpos);
void AddWinds(SPHbody **SPHbtabp, int *nobj, template_t *tempbtab, 
	      int windpartpershell, float r_wind, float v_wind, 
	      float mdot_wind, float u_wind, float *t_wind, float tpos, 
	      float dt, float *dt_next, float openangle_wind);
void AddNonconstWinds(SPHbody **SPHbtabp, int *nobj, template_t *temptab, 
		      int windpartpershell, winddata_t *wdata, int wnobj, 
		      float r_wind, float r_outer, float *t_wind, float tpos, 
		      float dt, float *dt_next, float openangle_wind);
void AddAccreting(SPHbody **SPHbtabp, int *nobj, template_t *tempbtab, 
		  int windpartpershell, float r_wind, float v_wind, 
		  float mdot_wind, float u_wind, float *t_wind, float tpos, 
		  float dt, float *dt_next, float omega);
void ReadTemplate(char *filename, template_t **temptab, int *tempnobj);
void ReadWindData(char *filename, winddata_t **wdata, int *wnobj);

static int maxmem(void);
static int maxheap(void);

float Znow(float time);
float Hnow(float time);

/* can some of these go into the function where they're solely used? */
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
static float dt_wind;
static float R0;
static float this_tol, this_eps;

static double gnterms;
/* static double ggravnterms; */

static bndry_t bndry;
setup_params_t params;

double massCF;
double lenCF;
double timeCF;

double ivmassCF, ivtimeCF, ivlenCF;
double timeCF2, ivtimeCF2;
double lenCF2, ivlenCF2,  ivlenCF3;
double ldivtCF, tdivlCF;

double grav_c, c_light;

int NNW;	/* number of isotopes in network */
int **inNW;
int nparr[NISO], nnarr[NISO];
float **tablep; /*array to hold cooling curve table values*/
float **ionfracp; /*array to hold ionfraction table values*/

/*
int **inNW;
int nparr[NISO], nnarr[NISO];
*/

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
    FILE *fp = NULL;
    int gnobj, nobj;
    int SPHgnobj, SPHnobj, SPHoldnobj;
    int windgnobj, windnobj;
    int PMgnobj, PMnobj;
    int SPHsinkgnobj, SPHsinknobj;
    body *btab, *p;
    body *pmtab;
    SPHbody *SPHbtab, *SPHsinkbtab = NULL, *q;
    windbody *windbtab;
    int i;
    float rmin[NDIM], rmax[NDIM];
    int first_step = 1;
    int added_particles = 0;
    int stride = sizeof(body)/sizeof(float);
    int SPHstride = sizeof(SPHbody)/sizeof(float);
    int SPHstride2 = sizeof(SPHbody)/sizeof(double);
    int SPHstride3 = sizeof(SPHbody)/sizeof(unsigned int);
    int iter;
    SDF *csdfp;			/* SDF pointer to control file */
    SDF *sdfp = NULL;
    float tposlast;
    double pe, ke, te;
    double dark_ke, dark_pe;
    double etot;
    double mtot, SPHmtot;
    sortresult_t sortedbtab, SPHsortedbtab, sortedatab;
    tree_t thetree, SPHtree, SPHsinktree, *sinkptr = NULL;
    int MACtype = BMAX_MAC;
    inherit_t inherit;
    macv_t mac;
    float dt_last;
    template_t *tempbtab;
    winddata_t *wdata;
    int tempnobj;
    int wnobj;
    float newp[NDIM], newl[NDIM];
    float newr = 0.0;
    float newmass = 0.0, totnewmass = 0.0, invmass = 0.0;
    float drag_coeff;
/*     float newmass = 0.0, totnewmass = 0.0; */
    int udot_limit[2];
    float vsz;
    float tmin;
    int tbad;
    accbody *SPHatab, *pa;
    int SPHanobj;
    void *decomp_info = NULL;
    float dark_tacc = -1e30;	/* initialize so dark_need_update is true */
    int did_dark_update;
    int SPHnupdate;
    int make_sink_tree;
    int Gridpts = 0, Nel = 0; 	/* for cooling tables */
    int status, done,rank,idbug;
	char netrcfn[20] = "                    "; /* build.f reserves 20 chars for net.rc filename */
    char **pnames, **nnames;
    int calc_gamma = 0;
    float tot_u, tot_pv;

/*
    argv[1]="/scratch/cellinge/runsnsph/casa16run4.ctl";
*/

    MPMY_Init(&argc, &argv);
    singlPrintf("Welcome to the variable O() integrator running on %d procs\n",
		MPMY_Nproc());

    csdfp = startup(argc, argv);

    SetBoundary(30); /* From integrate.c; if you change this, also change
			it in SPHDiags */

	read_initial_ctl(csdfp, &params);
    if (params.timeout > 0) MPMY_TimeoutSet(params.timeout);

	if (params.do_cooling || params.do_burning) {
/* read in Z and N for abundances. at some later point, populate nparr/nnarr ~CIE*/
	    if( (fp = fopen("networklist","r"))==NULL ) printf("error opening networklist\n");
	    fscanf(fp, "%d", &NNW);
	    inNW = (int **)malloc( 2 * sizeof(int *) );
	    inNW[0] = (int *)malloc( (NNW+1) * sizeof(int) );
	    inNW[1] = (int *)malloc( (NNW+1) * sizeof(int) );
	    for (i = 0; i < NNW; i++) fscanf(fp, "%d\t%d", &inNW[0][i], &inNW[1][i]);
	    inNW[0][NNW] = 0; /*for electron fraction */
	    inNW[1][NNW] = 0;
	    fclose(fp);
	    fp = NULL;
	 
	    make_spec_names(&pnames, 'p', NISO);
	    make_spec_names(&nnames, 'n', NISO);
	}

    if (params.do_sph || params.do_grav) {
		if (!((strncmp(params.name, "test", 4) == 0))) {
			if (strlen(params.SPHdatafile) > 0 || params.do_restart) {
				char iname[256];
	
				if (params.do_sph) {
				    if (params.do_restart) sprintf(iname, "%s_sph.restart", params.name);
					else sprintf(iname, "%s", params.SPHdatafile);
	/* this is where the SDF file is read in -CIE */
					if (params.do_cooling || params.do_burning) {
					    sdfp = SPHReadA(iname, csdfp, &SPHbtab, &SPHgnobj, &SPHnobj,
							   params.set_id, params.setpvel, params.new_h, params.new_u);
			            for ( i=0; i<NISO; i++ ) {
			            /* need to create the 'p1/n1' specifiers */
			                SDFgetintOrDie(sdfp, pnames[i], &nparr[i]);
			                SDFgetintOrDie(sdfp, nnames[i], &nnarr[i]);
			            }
					} else {
						sdfp = SPHRead(iname, csdfp, &SPHbtab, &SPHgnobj, &SPHnobj, params.set_id, params.setpvel, params.new_h, params.new_u);
					}
				} else SPHgnobj = SPHnobj = 0;
	
				if (params.has_grav_data) {
				    if (params.do_restart) sprintf(iname, "%s.restart", params.name);
					else sprintf(iname, "%s", params.name);
				    sdfp = DarkRead(iname, csdfp, (void **)&btab, &gnobj, 
						    &nobj, params.set_id, params.setpvel);
				} else {
				    gnobj = nobj = 0;
				    btab = Malloc(sizeof(body)); /* realloced later */
				}
	
				if (params.do_point_mass) {
				    if (params.do_restart) sprintf(iname, "%s.restart", params.name);
					else sprintf(iname, "%s", params.name);
				    sdfp = DarkRead(iname, csdfp, (void **)&pmtab, &PMgnobj, 
						    &PMnobj, params.set_id, params.setpvel);
				    Msgf(("lx = %e; ly = %e; lz = %e; accmass = %e\n", 
					  pmtab->l[0], pmtab->l[1], pmtab->l[2], 
					  pmtab->accmass));
				} else if (params.do_point_mass2) {
				    PMgnobj = PMnobj = 0;
				    pmtab = Malloc(sizeof(body)); /* realloced later */
					if (params.timeout > 0) MPMY_TimeoutSet(params.timeout);
				} else {
				    PMgnobj = PMnobj = 0;
				    pmtab = Malloc(sizeof(body)); /* realloced later */
				}
	
		        /* first iteration, get central particle data from ctl file */
		        if (params.do_absorbing_bndry) {
		 			read_absorb_bndry(csdfp, &bndry);
		        }
	
		        if (params.do_winds) {
		            if (params.old_winds) {
		                sdfp = WindRead(iname, csdfp, &windbtab, &windgnobj, 
		                        &windnobj);
		            } else windgnobj = windnobj = 0;
		
		            if (params.const_winds || params.nonconst_winds || params.accreting_winds) {
		                if (MPMY_Procnum() == 0) {
		                    ReadTemplate(params.template_name, 
		                            (template_t **)&tempbtab, &tempnobj);
		                    if (tempnobj != params.windpartpershell) 
		                        Error("tempnobj != windpartpershell");
		                }
		            }
		            if (params.nonconst_winds) {
		                if (MPMY_Procnum() == 0)
		                    ReadWindData(params.winddata_name,
		                            (winddata_t **)&wdata, &wnobj);
		            }
		        } 
	
				/* on restart, get dt from sdf file */
		        SDFgetfloatOrDefault(sdfp, "dt", &dt, 0.0);
				if (params.do_restart)
					params.dt = dt;
	        } else {
	            sdfp = InitRead(params.name, csdfp, (void **)&btab, &gnobj, &nobj, 
	                    &SPHbtab, &SPHgnobj, &SPHnobj, 
	                    params.set_id, params.setpvel, params.new_h, params.new_u);
	        }
	        FixNterms(btab, nobj);
	        SPHFixNterms(SPHbtab, SPHnobj);
	        /*SDFgetfloatOrDie(csdfp, "Gnewt", &cosmo.GNewt);*/
	        SDFgetfloatOrDefault(sdfp, "tpos",  &tpos, (float)0.0);
	        tvel = tpos;
	        SDFgetintOrDefault  (sdfp, "iter",  &iter, 0);

	        /* not first iteration, get central data particle from sdf file */
	        if (params.do_absorbing_bndry && (iter > 0)) {
				read_absorb_bndry(sdfp, &bndry);
	        }
	        if (params.cosmology) ReadCosmo(sdfp, &cosmo, tpos, &R0);
	        if(sdfp) SDFclose(sdfp);
	    } else {
	        SPHTestData(csdfp, &SPHbtab, &SPHgnobj, &SPHnobj, params.do_periodic);
	        nobj = gnobj = 0;
	        cosmo.GNewt = (float)1.0;
	        tvel = tpos = (float)0.0;
	        iter = 0;
	        if (params.do_periodic) R0 = 1.0;
	    }
    }
    singlPrintf("Maxmem after data read is %d (%d)\n", maxmem(), maxheap());
    if( Msg_test("memleak") ){
        Msg_do("Memory map after data read\n");
        malloc_print();
    }

    if (!params.do_restart) {
		dt = params.dt;
        /*  SDFgetfloatOrDefault(csdfp, "dark_dt", &dark_dt, params.do_grav ? dt : 1e30); */
    }
    massCF= (double)params.fmassCF;
    lenCF= (double)params.flenCF;
    timeCF= (double)params.ftimeCF;

    if (params.do_Bmax) MACtype = BMAX_MAC;
    else if (params.do_BH) MACtype = BH_MAC;
    else if (params.do_Arel) MACtype = AREL_MAC;
    else if (params.do_grav) Error("No MAC specified\n");

    if(csdfp) 
        SDFclose(csdfp);

    if (params.do_periodic) {
        EnableTimer(&FixCubeTm, "Fix Cube");
    }

    /* at this point, calculate some useful constants, factors, 
       so these don't have to be computed every time they're needed?
       (i.e. for every particle) */
    ivlenCF = 1./lenCF;
    ivmassCF = 1./massCF;
    ivtimeCF = 1./timeCF;

    ldivtCF = lenCF*ivtimeCF;
    tdivlCF = ivlenCF*timeCF;
    timeCF2 = timeCF*timeCF;
    ivtimeCF2 = ivtimeCF*ivtimeCF;
    lenCF2 = lenCF*lenCF;
    ivlenCF2 = ivlenCF*ivlenCF;
    ivlenCF3 = ivlenCF*ivlenCF*ivlenCF;

    cosmo.GNewt = GRAV_C *((double)massCF*ivlenCF *tdivlCF*tdivlCF);

    grav_c = cosmo.GNewt;
    c_light = C_LIGHT * tdivlCF;

    singlPrintf("float params.dt = %g;\n", params.dt);
    singlPrintf("int iter = %d;\n", iter);
    singlPrintf("int nproc = %d;\n", MPMY_Nproc());
    singlPrintf("float Gnewt = %g;\n", cosmo.GNewt);

	print_initial_ctl(params);

    if (params.do_winds) {
        if (params.old_winds) {
            singlPrintf("int windgnobj = %d;\n", windgnobj);
        }
    }
    if (params.do_absorbing_bndry) {
		print_absorb_bndry (bndry);
    }
    if (params.cosmology) {
        singlPrintf("float R0 = %f;\n", R0);
    }

    singlFflush();
    if (params.do_sph) SPHSanityCheck(SPHbtab, SPHnobj, SPHgnobj, &SPHmtot);

    /* read in necessary files on all processors */
    /*read in cooling curves and ion fraction tables*/
    /*
       if(params.do_cooling || params.do_burning) {
       }
       read in cooling curves and ion. tables in any case. need for temp calculation
       */
	if (params.do_cooling || params.do_burning) {
		singlPrintf("reading in cooling tables .... ");
	    init_CoolTable(&Gridpts, &Nel);
	    singlPrintf("success!\n");
	}

    rank = MPMY_Procnum();

    /*set up network for burn code. do this AFTER do_burning is set!!*/
    if(params.do_burning) {
        /* each processor needs its own 'net.rc' file */
		sprintf(netrcfn, "net.rc.%-d",rank);
        singlPrintf("building network library .... ");
        build_(&rank,&idbug,netrcfn);
        singlPrintf("success!\n");
    }

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
        q->dt = q->dt_next = params.dt;
        q->tacc = -1e30;
    }

    dt_last = params.dt;
    SPH_setup(NDIM, params.kernel_ncoef1, params.kernel_coef1, params.kernel_ncoef2, params.kernel_coef2);
    inherit = (inherit_t)InheritSinkNlogN;

    if (params.do_DL)
        mac = (macv_t)DLRcritMAC;
    else
        mac = (macv_t)RcritMAC;

#ifdef SPH_GRAV
    mac = (macv_t)SPHDLRcritMAC;
#endif

    this_eps = params.eps;
    this_tol = params.tol;

    /* Testing initialization */
    for (q = SPHbtab; q < SPHbtab+SPHnobj; q++) {
        VS(q->acc, = 0.0);
        VS(q->acc_last, = 0.0);
        VS(q->grav_acc, = 0.0);
        q->phi = 0.0;
    }

    for (params.nsteps += iter; iter <= params.nsteps; iter++) {
        if (params.timeout > 0) MPMY_TimeoutReset(params.timeout);
        /* Reset timers and counters */
        ClearEnabledTimers();
        ClearEnabledCounters();
        StartTimer(&StepTotWC);
        StartTimer(&StepTot);

        tot_u = 0.0;
        tot_pv = 0.0;

        /* calculating gamma */
        if(calc_gamma) {
            for( q = SPHbtab; q < SPHbtab+SPHnobj; q++) {
                tot_u += q->u*q->mass;
                tot_pv += q->pr * q->h*q->h*q->h;
            }

            /* must first MPI_COMBINE these things!! */
            MPMY_Combine(&tot_u, &tot_u, 1, MPMY_FLOAT, MPMY_SUM);
            MPMY_Combine(&tot_pv, &tot_pv, 1, MPMY_FLOAT, MPMY_SUM);

        }

        if (params.do_point_mass2) {
            for(q = SPHbtab; q < SPHbtab+SPHnobj; q++) {
                /* Even if grav isn't on, I still want to zero phi at the
                   beginning of each iteration */
                q->phi = 0.0;
            }
        }

        if (params.do_point_mass || params.do_point_mass2) {
            SPHoldnobj = SPHnobj;
            /*  	  ShrinkBtab((SPHbody **)&SPHbtab, pmtab, &SPHnobj, params.r_inner); */
            /*   	  ShrinkBtab2((SPHbody **)&SPHbtab, &SPHnobj, params.r_inner);  */
            AdjustBtab((SPHbody **)&SPHbtab, &SPHnobj, SPHgnobj, windbtab, 
                    windnobj, params.windpartpershell, params.r_inner, dt_last, iter, tpos, 
                    &added_particles);
            /*  	  AdjustBtab2((SPHbody **)&SPHbtab, &SPHnobj, SPHgnobj, windbtab,  */
            /* 		      windnobj, params.r_inner, dt_last, iter, tpos,  */
            /* 		      &added_particles, &newmass); */

            MPMY_Combine(&SPHnobj, &SPHgnobj, 1, MPMY_INT, MPMY_SUM);
            /* 	  MPMY_Combine(&newmass, &totnewmass, 1, MPMY_FLOAT, MPMY_SUM); */
            /* 	  params.centmass += totnewmass; */
            Msgf(("Iter: %d: Added %d bodies to SPHbtab\nBH mass = %f\n", iter, 
                        SPHnobj-SPHoldnobj, params.centmass));
            /* SPHFixId(SPHbtab, SPHnobj, SPHgnobj); */
        }

        if (params.do_winds && params.const_winds) {
            SPHoldnobj = SPHnobj;
            /* Only proc 0 adds wind particles; would be better to add
               particles only to node that included particles close to
               wind in the first place */
            dt_wind = 1e30;
            if (MPMY_Procnum() == 0) {
                AddWinds((SPHbody **)&SPHbtab, &SPHnobj, tempbtab, 
                        params.windpartpershell, params.r_wind, params.v_wind, params.mdot_wind, 
                        params.u_wind, &(params.t_wind), tpos, dt_last, &dt_wind, 
                        params.openangle_wind);
            }
            MPMY_Combine(&SPHnobj, &SPHgnobj, 1, MPMY_INT, MPMY_SUM);
            MPMY_Combine(&dt_wind, &dt_wind, 1, MPMY_FLOAT, MPMY_MIN);
            Msgf(("Iter: %d: Added %d bodies to SPHbtab\n", iter,
                        SPHnobj-SPHoldnobj));
        } else if (params.do_winds && params.nonconst_winds) {
            SPHoldnobj = SPHnobj;
            /* Only proc 0 adds wind particles; would be better to add
               particles only to node that included particles close to
               wind in the first place */
            dt_wind = 1e30;
            AddNonconstWinds((SPHbody **)&SPHbtab, &SPHnobj, tempbtab, 
                    params.windpartpershell, wdata, wnobj, params.r_wind, 
                    params.r_outer, &(params.t_wind), tpos, dt_last, &dt_wind,
                    params.openangle_wind);
            MPMY_Combine(&SPHnobj, &SPHgnobj, 1, MPMY_INT, MPMY_SUM);
            MPMY_Combine(&dt_wind, &dt_wind, 1, MPMY_FLOAT, MPMY_MIN);
            Msgf(("Iter: %d: Added %d bodies to SPHbtab\n", iter,
                        SPHnobj-SPHoldnobj));
        } else if (params.do_winds && params.accreting_winds) {
            SPHoldnobj = SPHnobj;
            dt_wind = 1e30;
            if (MPMY_Procnum() == 0) {
                AddAccreting((SPHbody **)&SPHbtab, &SPHnobj, tempbtab, 
                        params.windpartpershell, params.r_wind, params.v_wind, params.mdot_wind, 
                        params.u_wind, &(params.t_wind), tpos, dt_last, &dt_wind, 
                        params.omega_wind);
            }
            MPMY_Combine(&SPHnobj, &SPHgnobj, 1, MPMY_INT, MPMY_SUM);
            MPMY_Combine(&dt_wind, &dt_wind, 1, MPMY_FLOAT, MPMY_MIN);
            Msgf(("Iter: %d: Added %d bodies to SPHbtab\n", iter,
                        SPHnobj-SPHoldnobj));
        } else dt_wind = 1e30;

        if (params.do_boundary) {
            SPHoldnobj = SPHnobj;
            AdjustBtab3((SPHbody **)&SPHbtab, &SPHnobj, SPHgnobj, params.r_inner, 
                    params.r_outer);
            MPMY_Combine(&SPHnobj, &SPHgnobj, 1, MPMY_INT, MPMY_SUM);
            Msgf(("Iter: %d: Removed %d bodies from SPHbtab\n", iter, 
                        SPHoldnobj-SPHnobj));
        }

        if (params.do_absorbing_bndry) {
            SPHoldnobj = SPHnobj;
            AdjustBtab4((SPHbody **)&SPHbtab, &SPHnobj, bndry, &newmass, &newr,
                    newp, newl, cosmo.GNewt, dt);

            totnewmass = 0.0;
            MPMY_Combine(&SPHnobj, &SPHgnobj, 1, MPMY_INT, MPMY_SUM);
            MPMY_Combine(&newmass, &totnewmass, 1, MPMY_FLOAT, MPMY_SUM);
            MPMY_Combine(&newr, &newr, 1, MPMY_FLOAT, MPMY_MIN);
            MPMY_Combine(newp, newp, 3, MPMY_FLOAT, MPMY_SUM);
            MPMY_Combine(newl, newl, 3, MPMY_FLOAT, MPMY_SUM);

            bndry.mass += totnewmass;
            bndry.r = newr;
            VV(bndry.p, += newp);
            VV(bndry.l, += newl);
            /* what's the difference between using 'VV' and 'VVS' for this?
             * also, no VV ever divides, but always multiplies, and scalar is 
             * always first ~CIE */
            invmass = 1.0/bndry.mass;
            VVS(bndry.vel, = bndry.p, * invmass ); /* seems safer ~CIE */
            /* now, the call later to UpdateX will move central particle accordingly? */
            /* how do I prevent artificial added momentum due to rounding errors? */
            Msgf(("Iter %d: removed %d bodies from SPHbtab\nBndry mass = %g\n", iter, SPHoldnobj-SPHnobj, bndry.mass));
            Msgf(("p-vec = ( %g, %g, %g )\nl-vec = ( %g, %g, %g )\n", bndry.p[0],bndry.p[1],bndry.p[2],bndry.l[0],bndry.l[1],bndry.l[2]));
        }

        /* comoving smoothing */
        /* Note: behavior changed Jan. 25, 1996. Beware of old ctl files */
        if (params.comov_eps && (Znow(tpos)+1.0 >= params.comov_eps_epoch)) 
            this_eps = params.eps*params.comov_eps_epoch/(Znow(tpos)+(float)1.0);
        else this_eps = params.eps;

        /* Add sph particles to btab for gravity */
        if (params.do_grav) {
            GravPlusSPH((void **)&btab, &nobj, SPHbtab, SPHnobj);
            SanityCheck(btab, nobj, gnobj+SPHgnobj, &mtot); /* need mtot */

            if (params.do_periodic) {
                if (params.cosmology)
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

        if (params.do_grav && dark_need_update(dark_tacc, params.dark_dt)) {
            /* We aren't using the first two params */
            SetTol(0, 0, cosmo.GNewt, this_eps, gnobj+SPHgnobj);
            FixKeys(btab, nobj, GETKEY);

            if (MACtype == AREL_MAC) this_tol = params.tol*mtot/(sysradius*sysradius);
            SetupCofm(MACtype, this_tol, params.frac_tol);
            singlPrintf("BuildTree, tol=%g, frac_tol=%g, sysrad.=%.3g\n", this_tol, params.frac_tol,sysradius);

            StartTimer(&BuildTot);
            pqsortsetup(&sortedbtab, btab, nobj, sizeof(body), params.sort_tol, Realloc_f);
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
            if (params.do_periodic) {
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
            if (params.cosmology) FixGlobalForce(btab, nobj);
            if (params.do_periodic) FixCube(btab, nobj, sysradius, cosmo.GNewt*mtot);
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
        if (params.setpvel) {
            params.setpvel = 0;
            set_vels(btab, nobj, tpos);
            singlPrintf("Velocities adjusted to linear approximation.\n");
        }

        /* Do image before GravMinusSPH if you want to image all particles */
        if (params.image_freq && iter%params.image_freq == 0) {
            char name[256];
            float sysr, image_rmin[3], image_rmax[3];

            if (params.cosmology)
                sysr = R0 * (1.0+1e-5) / (1.0 + Znow(tpos));
            else
                sysr = R0;
            VS(image_rmin, = -sysr);
            VS(image_rmax, = sysr);
            FixRsizeExact(image_rmin, image_rmax);

            sprintf(name, "%s_img.%04d", params.outnamebase, iter);
            Image(btab[0].pos, btab[0].pos+1, &(btab[0].mass),
                    sizeof(body), nobj, image_rmin, image_rmax, 
                    params.x_pixels, params.y_pixels, 10, 250, params.log_image, name);
        }

        if (params.do_grav) GravMinusSPH((void **)&btab, &nobj, &SPHatab, &SPHanobj);

        /* This should be the high-water mark for memory use */
        AddCounter(&MemCnt, malloc_used()/1024);

        if (params.do_sph && (first_step || params.exact_rho)) {
            singlPrintf("BuildTree\n");
            StartTimer(&BuildTot);
            pqsortsetup(&SPHsortedbtab, SPHbtab, SPHnobj, sizeof(SPHbody), params.sort_tol, Realloc_f);
            SPHFixKeys(SPHbtab, SPHnobj, SPHGetKey);
            BuildTree(&SPHtree, &SPHsortedbtab);
            SPHbtab = SPHsortedbtab.data;
            SPHnobj = SPHsortedbtab.nobj;
            StopTimer(&BuildTot);
            StartTimer(&RhoSPH);
            SetSPH(params.visc_alpha, params.visc_beta, params.visc_epsilon, params.heat_f1,
                    params.Gamma, SPHgnobj, macRho, nbrMAC);
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
            if (params.do_periodic) {
                singlPrintf("FindRho (periodic)\n");
                vsz = (params.cosmology) ? 2.0*sysradius*Hnow(tvel) : 0.0;
                PeriodicSPH(&SPHtree, 2.0*sysradius, vsz);
                singlPrintf("FindRho (periodic) done\n");
            }
            singlPrintf("FindRho\n");
            WalkNT(&SPHtree);
            WalkTerminate();
            singlPrintf("updating final ....");
            update_final(SPHbtab, SPHnobj, Gridpts, Nel, dt, &udot_limit[0], &udot_limit[1],rank,tpos, R0);
            singlPrintf("updated final\n");
            /*update_final(SPHbtab, SPHnobj, dt, &udot_limit[0], &udot_limit[1]);*/
            StopTimer(&RhoSPH);
            FreeTree(&SPHtree);
            singlPrintf("FreeTree done\n");
        }

        if (params.do_sph) {
            SPHFixKeys(SPHbtab, SPHnobj, SPHGetKey);
            /* This sets rho_est and pr for communication during BuildTree */
            update_intermediate(SPHbtab, SPHnobj, Gridpts, Nel, dt_last, 
                    !(first_step || params.exact_rho), 0, sysradius);

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
                pqsortsetup(&SPHsortedbtab, SPHsinkbtab, SPHsinknobj, sizeof(SPHbody), params.sort_tol, Realloc_f);
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
            pqsortsetup(&SPHsortedbtab, SPHbtab, SPHnobj, sizeof(SPHbody), params.sort_tol, Realloc_f);
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
                pqsortsetup(&sortedatab, SPHatab, SPHanobj, sizeof(accbody), params.sort_tol, Realloc_f);
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
            SetSPH(params.visc_alpha, params.visc_beta, params.visc_epsilon, params.heat_f1, params.Gamma, SPHgnobj, 
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
            if (params.do_periodic) {
                singlPrintf("ForceSPH (periodic)\n");
                vsz = (params.cosmology) ? 2.0*sysradius*Hnow(tvel) : 0.0;
                PeriodicSPH(sinkptr, 2.0*sysradius, vsz);
                singlPrintf("ForceSPH (periodic) done\n");
            }
            StopTimer(&PerTmSPH);
            singlPrintf("ForceSPH\n");
            WalkNT(sinkptr);
            WalkTerminate();
            singlPrintf("ForceSPH done\n");
            udot_limit[0] = udot_limit[1]  = 0;
            update_final(SPHsinkbtab, SPHsinknobj, Gridpts, Nel, dt, &udot_limit[0], &udot_limit[1],rank,tpos, R0);
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
                pqsortsetup(&SPHsortedbtab, SPHsinkbtab, SPHsinknobj, sizeof(SPHbody), params.sort_tol, Realloc_f);
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

        if (params.do_point_mass) {
            /* Need to add code for parallel stuff here */
            for (p = pmtab; p < pmtab+PMgnobj; p++) {
                VS(p->acc, = 0.0);
                p->phi = 0.0;
            }
            for (p = pmtab; p < pmtab+PMgnobj; p++) {
                update_point_mass(pmtab, PMnobj, p, params.eps*params.eps, cosmo.GNewt);
            }
            for (p = pmtab; p < pmtab+PMgnobj; p++) {
                update_point_SPHmass(SPHbtab, SPHnobj, p, params.eps*params.eps, cosmo.GNewt);
            }
            singlPrintf("Updated %d point-mass accs\n", PMgnobj);
        }

        if (params.do_point_mass2 || params.do_boundary) {
            update_point_SPHmass2(SPHbtab, SPHnobj, params.eps*params.eps, cosmo.GNewt, 
                    params.centmass);
        }

        if (params.do_absorbing_bndry) {
            /*updates acc of all particles due to central particle (?) ~CIE */
            update_point_SPHmass_bndry(SPHbtab, SPHnobj, cosmo.GNewt, bndry);
        }

        if (params.do_drag) {
            for (q = SPHbtab; q < SPHbtab+SPHnobj; q++) {
                /* Just throws energy away */
                VV(q->acc, -= drag_coeff*q->vel);
            }
        }

        MPMY_Sync();

        if (ForceOutput()
                || (params.do_output && !first_step
                    && ((iter+params.output_freq) % params.output_freq == 0))
                || (params.save_first && first_step)) {
            if (params.do_sph) {
                if (params.do_winds) { 
                    if (params.old_winds) {
                        if (params.short_output) 
                            ShortWindOutput(SPHbtab, SPHnobj, windbtab, 
                                    windnobj, params.outnamebase, iter);
                        else
                            WindOutput(SPHbtab, SPHnobj, windbtab, 
                                    windnobj, params.outnamebase, iter);
                    } 
                    else
                        SPHOutput(SPHbtab, SPHnobj, params.outnamebase, iter);
                }
                else {
					if (params.do_cooling || params.do_burning)
					    SPHOutputA(SPHbtab, SPHnobj, params.outnamebase, iter);
					else
					    SPHOutput(SPHbtab, SPHnobj, params.outnamebase, iter);
				}
            }
            if (params.has_grav_data) Output(btab, nobj, params.outnamebase, iter);
            if (params.do_point_mass) Output(pmtab, PMnobj, params.outnamebase, iter);
        }


        if (ForceStop()) {
            singlPrintf("Stopping.\n");
            break;
        }

        Msgf(("integrating positions\n"));
        dark_ke = dark_pe = 0.0;
        if (params.has_grav_data) {
            for (p = btab; p < btab+nobj; p++) {
                dark_ke += 0.5 * p->mass * Dot(p->vel, p->vel);
                dark_pe += 0.5 * p->mass * p->phi;
            }
        }
        if (params.do_point_mass) {
            for (p = pmtab; p < pmtab+PMnobj; p++) {
                dark_ke += 0.5 * p->mass * Dot(p->vel, p->vel);
                dark_pe += 0.5 * p->mass * p->phi;
            }
        }
        Fix_h(SPHbtab, SPHnobj, params.nbrcut_max, params.nbrcut_min, params.nbrcut_fac, params.max_h, params.min_h);
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
        if (params.do_absorbing_bndry) {
            /* to get central velocity from particle absorption, need to Update 
             * bndry.vel under consideration from absorbed (linear, +angular?) momentum?
             */
            UpdateX(bndry.pos, sizeof(bndry_t), bndry.vel, sizeof(bndry_t), 1, dt, dt_last);
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
                &SPHbtab[0].udot_last, SPHstride, &SPHbtab[0].ident,
                SPHstride3, SPHnobj, dt, dt_last);
        /* One must be careful with this integration scheme, since v */
        /* is a derived variable.  To really adjust v, change pos_last */
#ifdef POS_IS_DOUBLE
        PUpdateVsd(SPHbtab[0].vel, SPHstride, SPHbtab[0].pos, SPHstride2, 
                SPHbtab[0].pos_last, SPHstride2, SPHbtab[0].acc, SPHstride, 
                &SPHbtab[0].ident, SPHstride3, SPHnobj, dt, dt_last);
        /* v must be done before x, since pos_last is changed in PUpX */
        PUpdateXsd(SPHbtab[0].pos, SPHstride2, SPHbtab[0].pos_last, SPHstride2,
                SPHbtab[0].acc, SPHstride, &SPHbtab[0].ident, SPHstride3, 
                SPHnobj, dt, dt_last);
#else
        PUpdateVs(SPHbtab[0].vel, SPHstride, SPHbtab[0].pos, SPHstride, 
                SPHbtab[0].pos_last, SPHstride, SPHbtab[0].acc, SPHstride, 
                &SPHbtab[0].ident, SPHstride3, SPHnobj, dt, dt_last);
        /* v must be done before x, since pos_last is changed in PUpX */
        PUpdateXs(SPHbtab[0].pos, SPHstride, SPHbtab[0].pos_last, SPHstride,
                SPHbtab[0].acc, SPHstride, &SPHbtab[0].ident, SPHstride3, 
                SPHnobj, dt, dt_last);
#endif
        UpdateSXs(&SPHbtab[0].h, SPHstride, &SPHbtab[0].hdot, SPHstride, 
                &SPHbtab[0].ident, SPHstride3, SPHnobj, dt, dt_last);
        tpos += dt;
        tvel += dt;
        dt_last = dt;
        for (i = 0; i < SPHnobj; i++) {
            VV(SPHbtab[i].acc_last, = SPHbtab[i].acc);
        }
        if(params.cosmology){
            CosmoPush(&cosmo, tpos);
            Msgf(("Pushed cosmo params to tpos=%g, Z=%g\n",
                        tpos, Znow(tpos)));
        }

        if (params.do_periodic) {
            if (params.cosmology)
                sysradius = R0 * 1.0 / (1.0 + Znow(tpos));
            else
                sysradius = R0;
            VS(rmin, = -sysradius);
            VS(rmax, = sysradius);
            WrapPeriodic(btab, nobj, rmin, rmax, 2.0*sysradius, params.cosmology,
                    params.log_time, tpos, dt_last);
            SPHWrapPeriodic(SPHbtab, SPHnobj, rmin, rmax, 2.0*sysradius, params.cosmology,
                    params.log_time, tpos, dt_last);
        }

        if (params.cosmology) 
            singlPrintf("\ntpos: %g znow: %.3f iter: %d size: %.2f, eps: %.0f\n", 
                    tposlast, Znow(tposlast), iter, sysradius, this_eps);
        else
            singlPrintf("\ntpos: %g iter: %d size: %f\n",
                    tposlast, iter, sysradius);

        etot = 0.0;
        if (params.has_grav_data) Diags(btab, nobj, dark_ke, dark_pe, &etot, dt_last, iter, gnobj);
        if (params.do_point_mass) Diags(pmtab, PMnobj, dark_ke, dark_pe, &etot, dt_last, iter, PMgnobj);
        if (params.do_sph) SPHDiags(SPHbtab, SPHnobj, ke, pe, te, &etot, dt_last, iter, SPHgnobj, 
                &tmin, &tbad);

        MPMY_Combine(udot_limit, udot_limit, 2, MPMY_INT, MPMY_SUM);

        /* Fear my nested ternary operators! */
        if (params.adaptive_dt) Fix_dt(&dt, 
                (((params.dark_dt<params.dt_max) ? params.dark_dt:params.dt_max) < dt_wind)
                ? ((params.dark_dt<params.dt_max) ? params.dark_dt:params.dt_max):dt_wind,
                params.tlow_cut, tmin, tbad, params.dt_short, params.dt_long,
                udot_limit[0], udot_limit[1]);

        singlPrintf("udot_limit high: %d low: %d\n", udot_limit[0], 
                udot_limit[1]);
        singlPrintf("Total Energy: %g\n", etot);
        singlPrintf("Central mass: %g; boundary radius: %g\n", 
                bndry.mass, bndry.r);

        StopTimer(&StepTot);
        StopTimer(&StepTotWC);

        AddCounter(&HeapCnt_, malloc_heapsz()/1024);

        if( params.timer_freq && iter%params.timer_freq == 0 ){
            OutputTimers(singlPrintf);
            OutputCounters(singlPrintf);
            if( Msg_test("timers") ){
                /* This can be very tedious on a big machine. */
                OutputIndividualTimers(Msg_do);
                OutputIndividualCounters(Msg_do);
            }
            if (params.ntimer_detail) {
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
        if (params.CWfac != 0.0) {	/* CWfac = 1 seems to work well for do_periodic */
            gnterms /= gnobj;
            singlPrintf("Avg nterms = %.0f, CWfac is %.2f\n", gnterms, params.CWfac);
            gnterms *= params.CWfac;
            /* Account for 'constant' work associated with each particle */
            for (p = btab; p < btab+nobj; p++) 
                p->nterms += gnterms;
        }

        if (params.SPHCWfac != 0.0) {
            gnterms /= SPHgnobj;
            singlPrintf("Avg nterms = %.0f, SPHCWfac is %.2f\n", gnterms,
                    params.SPHCWfac);
            gnterms *= params.SPHCWfac;
            for (q = SPHbtab; q < SPHbtab + SPHnobj; ++q)
                q->nterms += gnterms;

            /* 	    if (params.do_grav) { */
            /* 		ggravnterms /= SPHgnobj; */
            /* 		singlPrintf("Avg grav_nterms = %.0f, SPHCWfac is %.2f\n", ggravnterms, */
            /* 			    params.SPHCWfac); */
            /* 		ggravnterms *= params.SPHCWfac; */
            /* 		for (q = SPHbtab; q < SPHbtab + SPHnobj; ++q) */
            /* 		    q->grav_nterms += ggravnterms; */
            /* 	    } */
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
    /*free(tablep);*/ /*CE: necessary?*/
    /*free(ionfracp);*/
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
        output_btab[i].drho_dt = btab[i].drho_dt;
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
    if (params.cosmology) {
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
            "frac_tolerance", SDF_FLOAT, params.frac_tol,
            "ndim", SDF_INT, NDIM,
            "tpos", SDF_FLOAT, tpos_out,
            "tvel", SDF_FLOAT, tvel_out,
            "R0", SDF_FLOAT, output_R0,
            "Omega0", SDF_FLOAT, cosmo.Omega0,
            "H0", SDF_FLOAT, cosmo.H0,
            "Lambda_prime", SDF_FLOAT, cosmo.Lambda,
            "hubble", SDF_FLOAT, output_h,
            "redshift", SDF_FLOAT, output_z,
            "gamma", SDF_FLOAT, params.Gamma,
            "centmass", SDF_FLOAT, params.centmass, 
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


/* Write out nobj SPH particles; only one node should write wind particles */
static void ShortWindOutput(SPHbody *btab, int nobj, windbody *windbtab, 
        int windnobj, const char *outnamebase, int iter)
{
    int i;
    sortresult_t outputsort;
    SPHshortoutbody *output_btab;
    int output_nobj = nobj;
    float tpos_out = tpos;
    float tvel_out = tvel; /* changed in Integrate() */
    MPMY_Comm_request req;
    int output_gnobj;
    float output_z, output_h, output_R0;
    char outname[256];

    sprintf(outname, "%s_sph.%04d", outnamebase, iter);
    output_btab = Malloc(output_nobj * sizeof(SPHshortoutbody));
    for(i=0; i<output_nobj; i++){
        output_btab[i].mass = btab[i].mass;
        VV(output_btab[i].pos, = btab[i].pos);
        VV(output_btab[i].vel, = btab[i].vel);
        output_btab[i].u = btab[i].u;
        output_btab[i].h = btab[i].h;
        output_btab[i].rho = btab[i].rho;
        output_btab[i].nbrs = btab[i].nbrs;
        output_btab[i].ident = btab[i].ident;
        output_btab[i].windid = btab[i].windid;
    }
    Msgf(("Doing output of %d bodies\n", output_nobj));
    singlPrintf("Trying to sort output\n");
    pqsortsetup_order(&outputsort, output_btab, output_nobj,
            sizeof(SPHshortoutbody), 0.1F, 1, Realloc_f);
    output_btab = pqsort(&outputsort, UnityCost, (pq_keyproto)SPHShortOutIdentKey);
    output_nobj = outputsort.nobj;
    Msgf(("After pqsort, %d outbodies\n", output_nobj));
    MPMY_ICombine_Init(&req);
    MPMY_ICombine(&output_nobj, &output_gnobj, 1, MPMY_INT, MPMY_SUM, req);
    MPMY_ICombine_Wait(req);
    if (params.cosmology) {
        output_z = Znow(tpos_out);
        output_h = Hnow(tpos_out);
        output_R0 = R0;
    } else {
        output_z = 0.0;
        output_h = 0.0;
        output_R0 = sysradius;
    }
    SDFwritewind(outname, output_gnobj, output_nobj, 
            output_btab, windnobj, windbtab, sizeof(SPHshortoutbody), 
            sizeof(windbody), WINDOUTBODYDESC, SPHSHORTOUTBODYDESC, 
            /* "npart", SDF_INT, output_gnobj, */
            "iter", SDF_INT, iter,
            "dt", SDF_FLOAT, dt,
            "eps", SDF_FLOAT, this_eps,
            "Gnewt", SDF_FLOAT, cosmo.GNewt,
            "tolerance", SDF_FLOAT, this_tol,
            "frac_tolerance", SDF_FLOAT, params.frac_tol,
            "ndim", SDF_INT, NDIM,
            "tpos", SDF_FLOAT, tpos_out,
            "tvel", SDF_FLOAT, tvel_out,
            "R0", SDF_FLOAT, output_R0,
            "gamma", SDF_FLOAT, params.Gamma,
            "centmass", SDF_FLOAT, params.centmass, 
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
    int i,j;
    sortresult_t outputsort;
    SPHoutbody *output_btab;
    int output_nobj = nobj;
    float tpos_out = tpos;
    float tvel_out = tvel; /* changed in Integrate() */
    float twind_out = params.t_wind;
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
        output_btab[i].drho_dt = btab[i].drho_dt;
        output_btab[i].udot = btab[i].udot;
#ifdef SPH_SAVE_ACC
        VV(output_btab[i].acc, = btab[i].acc);
        VV(output_btab[i].acc_last, = btab[i].acc_last);
        output_btab[i].phi = btab[i].phi;
        output_btab[i].dt = btab[i].dt;
#endif
        output_btab[i].pr = btab[i].pr;
        output_btab[i].nbrs = btab[i].nbrs;
        output_btab[i].ident = btab[i].ident;
        output_btab[i].temp = btab[i].temp;
        output_btab[i].Y_el = btab[i].Y_el;
        output_btab[i].mfp = btab[i].mfp;
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
    if (params.cosmology) {
        output_z = Znow(tpos_out);
        output_h = Hnow(tpos_out);
        output_R0 = R0;
    } else {
        output_z = 0.0;
        output_h = 0.0;
        output_R0 = sysradius;
    }
    /* I'm guessing this writes whatever is in output_btab, matched to SPHOUTBODYDESC -CIE */
    SDFwrite(outname, output_gnobj, 
            output_nobj, output_btab, sizeof(SPHoutbody),
            SPHOUTBODYDESC,
            "npart", SDF_INT, output_gnobj,
            "iter", SDF_INT, iter,
            "dt", SDF_FLOAT, dt,
            "eps", SDF_FLOAT, this_eps,
            "Gnewt", SDF_FLOAT, cosmo.GNewt,
            "tolerance", SDF_FLOAT, this_tol,
            "frac_tolerance", SDF_FLOAT, params.frac_tol,
            "ndim", SDF_INT, NDIM,
            "tpos", SDF_FLOAT, tpos_out,
            "tvel", SDF_FLOAT, tvel_out,
            "t_wind", SDF_FLOAT, twind_out,
            "R0", SDF_FLOAT, output_R0,
            "Omega0", SDF_FLOAT, cosmo.Omega0,
            "H0", SDF_FLOAT, cosmo.H0,
            "Lambda_prime", SDF_FLOAT, cosmo.Lambda,
            "hubble", SDF_FLOAT, output_h,
            "redshift", SDF_FLOAT, output_z,
            "gamma", SDF_FLOAT, params.Gamma,
            "massCF", SDF_FLOAT, massCF,
            "lenCF", SDF_FLOAT, lenCF,
            "timeCF", SDF_FLOAT, timeCF,
            "centmass", SDF_FLOAT, params.centmass, 
            "bndry_x", SDF_FLOAT, bndry.pos[0],
            "bndry_y", SDF_FLOAT, bndry.pos[1],
            "bndry_z", SDF_FLOAT, bndry.pos[2],
            "bndry_vx", SDF_FLOAT, bndry.vel[0],
            "bndry_vy", SDF_FLOAT, bndry.vel[1],
            "bndry_vz", SDF_FLOAT, bndry.vel[2],
            "bndry_px", SDF_FLOAT, bndry.p[0],
            "bndry_py", SDF_FLOAT, bndry.p[1],
            "bndry_pz", SDF_FLOAT, bndry.p[2],
            "bndry_lx", SDF_FLOAT, bndry.l[0],
            "bndry_ly", SDF_FLOAT, bndry.l[1],
            "bndry_lz", SDF_FLOAT, bndry.l[2],
            "bndry_mass", SDF_FLOAT, bndry.mass,
            "bndry_r", SDF_FLOAT, bndry.r,
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

static void SPHOutputA(SPHbody *btab, int nobj, const char *outnamebase, int iter)
{
    SPHbody *p;
    int i,j;
    sortresult_t outputsort;
    SPHoutbody_NW *output_btab;
    int output_nobj = nobj;
    float tpos_out = tpos;
    float tvel_out = tvel; /* changed in Integrate() */
    float twind_out = params.t_wind;
    double ke, pe, te;
    MPMY_Comm_request req;
    int output_gnobj;
    float output_z, output_h, output_R0;
    char outname[256];
    char **pnames, **nnames;

    make_spec_names(&pnames, 'p', NISO);
    make_spec_names(&nnames, 'n', NISO);

    sprintf(outname, "%s_sph.%04d", outnamebase, iter);
    pe = ke = te = 0.0;
    for (p = btab; p < btab+nobj; p++) {
        ke += (float)0.5 * p->mass * Dot(p->vel, p->vel);
        te += p->mass * p->u;
        pe += (float)0.5 * p->mass * p->phi;
    }
    output_btab = Malloc(output_nobj * sizeof(SPHoutbody_NW));
    for(i=0; i<output_nobj; i++){
        output_btab[i].mass = btab[i].mass;
        VV(output_btab[i].pos, = btab[i].pos);
        VV(output_btab[i].vel, = btab[i].vel);
        output_btab[i].u = btab[i].u;
        output_btab[i].h = btab[i].h;
        output_btab[i].rho = btab[i].rho;
        output_btab[i].drho_dt = btab[i].drho_dt;
        output_btab[i].udot = btab[i].udot;
#ifdef SPH_SAVE_ACC
        VV(output_btab[i].acc, = btab[i].acc);
        VV(output_btab[i].acc_last, = btab[i].acc_last);
        output_btab[i].phi = btab[i].phi;
        output_btab[i].dt = btab[i].dt;
#endif
        output_btab[i].pr = btab[i].pr;
        output_btab[i].nbrs = btab[i].nbrs;
        output_btab[i].ident = btab[i].ident;
        output_btab[i].temp = btab[i].temp;
        output_btab[i].Y_el = btab[i].Y_el;
        output_btab[i].mfp = btab[i].mfp;
        for(j=0;j<NISO;j++){ /*will this work? -CIE: so far, it compiled and runs*/
            output_btab[i].abund[j] = btab[i].abund[j];
            /*
               output_btab[i].np[j] = btab[i].np[j];
               output_btab[i].nn[j] = btab[i].nn[j];
               */
        }
    }
    /*     Msg("output", ("Doing output of %d bodies\n", output_nobj)); */
    Msgf(("Doing output of %d bodies\n", output_nobj));
    singlPrintf("Trying to sort output\n");
    pqsortsetup_order(&outputsort, output_btab, output_nobj,
            sizeof(SPHoutbody_NW), 0.1F, 1, Realloc_f);
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
    if (params.cosmology) {
        output_z = Znow(tpos_out);
        output_h = Hnow(tpos_out);
        output_R0 = R0;
    } else {
        output_z = 0.0;
        output_h = 0.0;
        output_R0 = sysradius;
    }
    /* I'm guessing this writes whatever is in output_btab, matched to SPHOUTBODYDESC -CIE */
    SDFwrite(outname, output_gnobj, 
            output_nobj, output_btab, sizeof(SPHoutbody_NW),
            SPHOUTBODYDESC_NW,
            "npart", SDF_INT, output_gnobj,
            "iter", SDF_INT, iter,
            "dt", SDF_FLOAT, dt,
            "eps", SDF_FLOAT, this_eps,
            "Gnewt", SDF_FLOAT, cosmo.GNewt,
            "tolerance", SDF_FLOAT, this_tol,
            "frac_tolerance", SDF_FLOAT, params.frac_tol,
            "ndim", SDF_INT, NDIM,
            "tpos", SDF_FLOAT, tpos_out,
            "tvel", SDF_FLOAT, tvel_out,
            "t_wind", SDF_FLOAT, twind_out,
            "R0", SDF_FLOAT, output_R0,
            "Omega0", SDF_FLOAT, cosmo.Omega0,
            "H0", SDF_FLOAT, cosmo.H0,
            "Lambda_prime", SDF_FLOAT, cosmo.Lambda,
            "hubble", SDF_FLOAT, output_h,
            "redshift", SDF_FLOAT, output_z,
            "gamma", SDF_FLOAT, params.Gamma,
            "massCF", SDF_FLOAT, massCF,
            "lenCF", SDF_FLOAT, lenCF,
            "timeCF", SDF_FLOAT, timeCF,
            "centmass", SDF_FLOAT, params.centmass, 
            "bndry_x", SDF_FLOAT, bndry.pos[0],
            "bndry_y", SDF_FLOAT, bndry.pos[1],
            "bndry_z", SDF_FLOAT, bndry.pos[2],
            "bndry_vx", SDF_FLOAT, bndry.vel[0],
            "bndry_vy", SDF_FLOAT, bndry.vel[1],
            "bndry_vz", SDF_FLOAT, bndry.vel[2],
            "bndry_px", SDF_FLOAT, bndry.p[0],
            "bndry_py", SDF_FLOAT, bndry.p[1],
            "bndry_pz", SDF_FLOAT, bndry.p[2],
            "bndry_lx", SDF_FLOAT, bndry.l[0],
            "bndry_ly", SDF_FLOAT, bndry.l[1],
            "bndry_lz", SDF_FLOAT, bndry.l[2],
            "bndry_mass", SDF_FLOAT, bndry.mass,
            "bndry_r", SDF_FLOAT, bndry.r,
            "ke", SDF_DOUBLE, ke,
            "pe", SDF_DOUBLE, pe,
            "te", SDF_DOUBLE, te,
            pnames[0], SDF_INT, nparr[0],
            nnames[0], SDF_INT, nnarr[0],
            pnames[1], SDF_INT, nparr[1],
            nnames[1], SDF_INT, nnarr[1],
            pnames[2], SDF_INT, nparr[2],
            nnames[2], SDF_INT, nnarr[2],
            pnames[3], SDF_INT, nparr[3],
            nnames[3], SDF_INT, nnarr[3],
            pnames[4], SDF_INT, nparr[4],
            nnames[4], SDF_INT, nnarr[4],
            pnames[5], SDF_INT, nparr[5],
            nnames[5], SDF_INT, nnarr[5],
            pnames[6], SDF_INT, nparr[6],
            nnames[6], SDF_INT, nnarr[6],
            pnames[7], SDF_INT, nparr[7],
            nnames[7], SDF_INT, nnarr[7],
            pnames[8], SDF_INT, nparr[8],
            nnames[8], SDF_INT, nnarr[8],
            pnames[9], SDF_INT, nparr[9],
            nnames[9], SDF_INT, nnarr[9],
            pnames[10], SDF_INT, nparr[10],
            nnames[10], SDF_INT, nnarr[10],
            pnames[11], SDF_INT, nparr[11],
            nnames[11], SDF_INT, nnarr[11],
            pnames[12], SDF_INT, nparr[12],
            nnames[12], SDF_INT, nnarr[12],
            pnames[13], SDF_INT, nparr[13],
            nnames[13], SDF_INT, nnarr[13],
            pnames[14], SDF_INT, nparr[14],
            nnames[14], SDF_INT, nnarr[14],
            pnames[15], SDF_INT, nparr[15],
            nnames[15], SDF_INT, nnarr[15],
            pnames[16], SDF_INT, nparr[16],
            nnames[16], SDF_INT, nnarr[16],
            pnames[17], SDF_INT, nparr[17],
            nnames[17], SDF_INT, nnarr[17],
            pnames[18], SDF_INT, nparr[18],
            nnames[18], SDF_INT, nnarr[18],
            pnames[19], SDF_INT, nparr[19],
            nnames[19], SDF_INT, nnarr[19],
            /*
               pnames[20], SDF_INT, nparr[20];
               nnames[20], SDF_INT, nnarr[20];
               pnames[21], SDF_INT, nparr[21];
               nnames[21], SDF_INT, nnarr[21];
               */
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
    if (params.cosmology) {
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
            "frac_tolerance", SDF_FLOAT, params.frac_tol,
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
    /*     ggravnterms = 0.0; */
    avg_nbrs = 0.0;
    tlow = 0;
    for (p = btab; p < btab+nobj; p++) {
        avg_nbrs += p->nbrs;
        if (!(p->ident & (1<<30))) {  /* If not a boundary particle... */
            VV(com, += p->mass*p->pos);
            VV(comv, += p->mass*p->vel);
            VV(force, += p->mass*p->acc);
            sacc2 = Dot(p->acc, p->acc);
            acc2 += sacc2;
            mtot += p->mass;
            gnterms += p->nterms;
            /* 	    ggravnterms += p->grav_nterms; */
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
            dti *= params.courant_number;
            if (p->min_nbr_dt == 1e30) { 
                /* This could happen if there are no nbrs. */
                SeriousWarning("Ignoring min_nbr_dt of %g\n", p->min_nbr_dt);
            } else {
                p->dt = p->min_nbr_dt;
            }
            if (dti > params.dt_max && 2.0*p->dt < params.dt_max) p->dt_next = 2.0*p->dt;
            else if (dti > 2.0*p->dt) p->dt_next = 2.0*p->dt;
            else if (dti < p->dt) {
                p->dt *= 0.5;
                p->dt_next = p->dt;
            }
            if (dti < min_dt)
                min_dt = dti;
            if (dti < dt)
                tlow++;
            if (p->dt_next > params.dt_max) p->dt_next = params.dt_max;
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
        if (p->nbrs > max_nbrs)
            max_nbrs = p->nbrs;
        if (p->nbrs < min_nbrs)
            min_nbrs = p->nbrs;	    
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
    /*     MPMY_ICombine(&ggravnterms, &ggravnterms, 1, MPMY_DOUBLE, MPMY_SUM, req); */
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
    if (tbad >= tlow_cut || limit_high >= tlow_cut || limit_low >= tlow_cut ) {
        dtlongvote = 0;
        dtshortvote++;
    }

    /* singlPrintf("Votes: short: %d; long %d\n", dtshortvote, dtlongvote); */
    if (dtshortvote > dtshort) {
        singlPrintf(("Adjusting dt down by factor of 1/2\n"));
        *dt *= (float)(1./2.);
        dtshortvote = dtlongvote = 0;
    } else if (dtlongvote > dtlong && (2.0)**dt <= dt_max) {
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
    if (!params.dark_independent_dt) return 1;
    return (dark_tacc + dark_dt <= tpos + dt * 1.00001);
}

/*  static float  */
/*  IdtSPHGetCost(const SPHbody *ptr) */
/*  { */
/*      if (SPH_need_update(ptr)) */
/*        return (float) ptr->nterms; */
/*      else */
/*        return (float) params.default_nterms; */
/*  } */

    int
SPH_need_update(const SPHbody *p)
{
    if (!params.independent_dt) return 1;
    return (p->tacc + p->dt <= tpos + dt * 1.00001);
}

int make_spec_names(char ***specarr, char spec, int num)
{
    int i;
    char tmpchr[20];

    *specarr = (char **)malloc(num * sizeof(char *) );

    for( i = 0; i < num; i++ ){

        sprintf( tmpchr, "%c%-d\0", spec, (i+1));
        *(*specarr + i) = (char *)malloc( ( strlen(tmpchr) + 1) * sizeof(char) );
        sprintf((*specarr)[i], "%s",tmpchr); 
        /*printf("isotope specifier: %s  %d\n", specarr[i],(int)strlen(specarr[i]));*/

    }
    return i;
}
