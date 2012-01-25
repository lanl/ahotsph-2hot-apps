/* 
   Parallel SPH interpolation + some physics on the mesh
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <Assert.h>

#include "error.h"
#include "mpmy.h"
#include "Msgs.h"
#include "physics.h"
#include "physics_sph.h"
#include "wvt.h"
#include "ghosts.h"

#include "SDF.h"
#include "SDFwrite.h"
#include "singlio.h"
#include "sphinit.h"
#include "vop.h"
#include "decomp.h"
#include "consts.h"
#include "math.h"
#include <errno.h>


#define MAXCOEF 16
#define POSFIXED_FLAG (1<<30)
#define SPHFIXED_FLAG (1<<29)
/*1<<28 is already DUMMYSINK_FLAG */
#define DUAL_FLAG (1<<27)
#define SPECIAL0_FLAG (1<<26)
#define SPECIAL1_FLAG (1<<25)
#define ALL_FLAGS (1<<31 | 1<<30 | 1<<29 | 1<<28 | 1<<27 | 1<<26 | 1<<25)
/* Hmm, we're down to 16 Million with all these tags now. We should think about
   having a separate integer carried around with only "flag" content. This may
   be particularly useful if we want to have more information about the history
   of a particle. For example, we could check during runtime if a particle has
   fulfilled some criterion that makes it "interesting" for analysis later. */

static SDF *initfiles(int argc, char *argv[]);
static void SPHSanityCheck(SPHbody *btab, int nobj, int gnobj, double *mtotp);
static void AdjustBtab (SPHbody **SPHbtabp, int *nobj, int gnobj, 
			double *rmin, double *rmax);

static void AdjustBtab_Spherical (SPHbody **SPHbtabp, int *nobj, int gnobj, 
				  double innerbound, double outerbound);
static void AdjustBtab_Rho (SPHbody **SPHbtabp, int *nobj, int gnobj, 
			    double rhocut);

void MySPHFixId(SPHbody *btab, int nobj, int gnobj);


static void SPHOutput(SPHbody *btab, int nobj, const char *outname, int iter, int do_floatoutput);

Timer_t StepTot, StepTotWC, BuildTot;
Timer_t FindForcesTm;
Timer_t RhoSPH, ForceSPH, PerTmSPH;
Timer_t FixCubeTm;
Timer_t WTermTm, WNTTm, PerTm;
Timer_t EosTm;
Timer_t SDFreadTm;
Counter_t NbodyCnt;
Counter_t MemCnt;
Counter_t HeapCnt_;	/* HeapCnt is in the SunOS name space?! */
Counter_t NtermsCnt;
Counter_t SPHbodyCnt;

int do_diffusion = 0;
double tvel=0.;
double tpos=0.;
double this_eps=0.;
double this_tol=0.;
double frac_tol=0.;
double Gamma=1.6666666666666666;
struct cosmo_s{
    double t;
    double a;
    double H0;
    double Omega0;
    double Lambda;
    double GNewt;
    double b;  /* Cluster core radius for Plummer model */
    double Zel_f;		/* the 'f' factor for linearly growing modes,
				 used only in set_vel = 1/H*Ddot/D.  It's
				 very close to 1 (exactly?) for flat models. */
} cosmo;
static double dt=0.;
static double sysradius=0.;

int main(int argc, char *argv[]) 
{
    int iter;
    SDF *csdfp, *sdfp;
    SPHbody *btab, *SPHbtab, *p, **btabp;
/*     void *decomp_info = NULL; */
    sortresult_t sortedbtab;
    tree_t SPHtree;
    double mtot;
    int num[NDIM];  /* uniform mesh for now */
    double rmin[NDIM], rmax[NDIM];
    double outrmin[NDIM], outrmax[NDIM];
    double sysradius;
    double sort_tol = 0.01;
    int i, j;
    int gnobj, nobj, targetnobj;
    int kernel_ncoef1, kernel_ncoef2;
    double kernel_coef1[MAXCOEF], kernel_coef2[MAXCOEF];
    char outnamebase[256];
    char outdir[256];
    double totvol;
    int do_externalstart;
    double tothvol;
    int nghosts, gnghosts;
    double outerbound=320., innerbound=-1.;
    int nloop, nhloop, nmassloop;
    double targetneighbors;
    char startfile[256];
    int ngood;
    double npervol;
    int do_floatoutput;
    int keepcenterfixed;
    int do_hydrostatic; 
    int do_center, center_dual, center_sphfixed, center_posfixed;
    int special0_sphpoint, special1_sphpoint;
    double center_h, center_grav_mass;
    int do_eospolytrope;
    double kpolytrope;

    MPMY_Init(&argc, &argv);
    singlPrintf("Welcome to the SPH interpolator running on %d procs\n", 
		MPMY_Nproc());

    csdfp = initfiles(argc, argv);

    SDFgetintOrDefault(csdfp, "do_floatoutput", &do_floatoutput, 0);

    if (SDFhasname("kernel_ncoef1", csdfp)) {
      SDFgetintOrDie(csdfp, "kernel_ncoef1", &kernel_ncoef1);
      if (kernel_ncoef1 >= MAXCOEF) Error("Increase MAXCOEF\n");
      SDFgetintOrDie(csdfp, "kernel_ncoef2", &kernel_ncoef2);
      if (kernel_ncoef2 >= MAXCOEF) Error("Increase MAXCOEF\n");
      if (SDFseekrdvecs(csdfp, "kernel_coef1", 0, kernel_ncoef1, 
                       kernel_coef1, 0, NULL))
         Error("SDFread kernel_coef1 failed\n");
      if (SDFseekrdvecs(csdfp, "kernel_coef2", 0, kernel_ncoef2, 
                       kernel_coef2, 0, NULL))
         Error("SDFread kernel_coef2 failed\n");
    } else {
      /* Monaghan spline kernel is default */
      kernel_ncoef1 = kernel_ncoef2 = 4;
      kernel_coef1[0] = 1.0;           kernel_coef2[0] = 2.0;
      kernel_coef1[1] = 0.0;           kernel_coef2[1] = -3.0;
      kernel_coef1[2] = -3.0/2.0;      kernel_coef2[2] = 3.0/2.0;
      kernel_coef1[3] = 3.0/4.0;       kernel_coef2[3] = -1.0/4.0;
    }
    singlPrintf("Kernel: %g %g %g %g %g %g %g %g \n", kernel_coef1[0], 
		kernel_coef1[1], kernel_coef1[2], kernel_coef1[3], 
		kernel_coef2[0], kernel_coef2[1], kernel_coef2[2], 
		kernel_coef2[3]);

    SDFclose(csdfp);
    singlPrintf("number of arguments: %d %s", argc, argv[2]);

    /* What we need now is the function that loads the grid data*/
    double minr=0., maxr=1., minz=-1., maxz=1.;
    double mintheta=0., maxtheta=2*3.14159265;
    int dimr=128, dimz=50, dimtheta=256;
    double rhocut=0.;

    init_cylindricalgrid(&btab, &gnobj, &nobj, dimr, dimz, dimtheta, 
			 minr, maxr, minz, maxz, mintheta, maxtheta);

    singlPrintf("Nobj: %d, Gnobj: %d\n", nobj, gnobj);
    
    AdjustBtab_Rho (&btab, &nobj, &gnobj, rhocut);

    singlPrintf("Nobj: %d, Gnobj: %d\n", nobj, gnobj);
 /*    for (p=btab; p<btab+nobj; p++) */
/*       if (p->rho != 0.) */
/* 	singlPrintf("x %g, y %g, z %g, rho %g\n", p->pos[0], p->pos[1], */
/* 		  p->pos[2], p->rho); */

    sprintf(outdir, "");
    sprintf(outnamebase, "grid.sdf");
    iter=0;
    SPHOutput(btab, nobj, outnamebase, iter, do_floatoutput);

    free(btab);

}

double ***allocate3d(int xdim, int ydim, int zdim)
{
     int i, j;
     double ***a3d = (double ***)malloc(xdim * sizeof(double **));
     for(i = 0; i < xdim; i++) {
       a3d[i] = (double **)malloc(ydim * sizeof(double *));
       for(j = 0; j < ydim; j++)
	 a3d[i][j] = (double *)malloc(zdim * sizeof(double));
     }
     return a3d;
}

void free3d(double ***a3d, int xdim, int ydim, int zdim)
{
     int i, j;
     for(i = 0; i < xdim; i++) {
       for(j = 0; j < ydim; j++)
	 free(a3d[i][j]);
       free(a3d[i]);
     }
}


void init_cylindricalgrid(SPHbody **btabp, int *gnobj, int *nobj, int dimr, 
			  int dimz, int dimtheta, double minr, double maxr, 
			  double minz, double maxz, 
			  double mintheta, double maxtheta)
{
  SPHbody *btab, *q;
  int j, k, l;
  double r, z, theta, ddata;
  FILE *infile;

  *gnobj=dimr*dimz*dimtheta;
  *nobj=*gnobj;

  btab = Malloc((*nobj)*sizeof(SPHbody));
  /* remember to free again! */

  infile=fopen("dens_Q0.4.cdat","r");

  q=btab;
  for (l=0; l<dimtheta; l++) {
    theta=(double) mintheta+l/(dimtheta-1.)*(maxtheta-mintheta);
    for (k=0; k<dimz; k++) {
      z=(double) minz+k/(dimz-1.)*(maxz-minz);
      for (j=0; j<dimr; j++) {
	r=(double) minr+j/(dimr-1.)*(maxr-minr);
 	fread(&ddata, sizeof(ddata),1,infile); 
	if (l+k+j < 5) singlPrintf("r %g, z %g, theta %g, rho %g\n", 
				   r, z, theta, ddata);
	q->pos[0]=r*cos(theta);
	q->pos[1]=r*sin(theta);
	q->pos[2]=z;
	q->rho=ddata;
/* 	if (l+k+j < 5)  */
/* 	if (ddata != 0)  */
/* 	  singlPrintf("x %g, y %g, z %g, rho %g\n",  */
/* 				   q->pos[0], q->pos[1], q->pos[2], q->rho); */
	q++;
      }      
    }    
  }
  fclose(infile);	

  *btabp=btab;
  
}


static SDF *initfiles(int argc, char *argv[]) 
{
    SDF *csdfp;
    char msg_turn_on[512];
    char tmp[256];
    char msgdir[256];
    char *msgbase, *lastslash;

/*     if (argc < 3) */
/* 	SinglError("Usage: %s control-file data-file(s)\n", argv[0]); */

    if ( (csdfp = SDFopen(NULL, argv[1])) == NULL )
	SinglError("%s: couldn't SDFopen %s: %s\n", argv[0], argv[1], 
		   SDFerrstring);

    /* Set up message directory */
    if( SDFgetstring(csdfp, "msgbase", tmp, sizeof(tmp))==0 )
	msgbase = tmp;
    else {
	lastslash = strrchr(argv[0], '/');
	if (lastslash) 
	    msgbase = lastslash + 1;
	else 
	    msgbase = argv[0];
	sprintf(tmp, "misc.%s/msg", msgbase);
	msgbase = tmp;
    }
    sprintf(msgdir, "%s.%d", msgbase, MPMY_Procnum());
    MsgdirInit(msgdir);

    SDFgetstringOrDefault(csdfp, "Msg_turn_on", msg_turn_on, 
			  sizeof(msg_turn_on), "");
    Msg_turnon(msg_turn_on);

    EnableTimer(&StepTot, "Step Total");
    EnableWCTimer(&StepTotWC, "Step Tot(WC)");
    EnableTimer(&BuildTot, "Build Total");
    EnableTimer(&DecompTm, "Decomp");
    EnableTimer(&DecompCommTm, "DecompComm");
    EnableTimer(&SortTm, "Sort");
    EnableTimer(&MakeTreeTm, "Make Tree");
/*     EnableTimer(&FindForcesTm, "Force Eval"); */
/*     EnableTimer(&GravTm, "Grav Time"); */
/*     EnableTimer(&MACTm, "MAC Time"); */
/*     EnableTimer(&ForceSPH, "Force (SPH)"); */
/*     EnableTimer(&RhoSPH, "Rho (SPH)"); */
    EnableTimer(&WalkDeferTm, "Walk Defer");
    EnableTimer(&WTermTm, "WalkTerm");
    EnableTimer(&WNTTm, "WalkNT");
    EnableTimer(&EosTm, "Eos");
/*     EnableTimer(&PerTm, "Periodic"); */
/*     EnableTimer(&PerTmSPH, "Periodic SPH"); */

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
    singlPrintf("Sanity check: gnobj = %d, mtot = %f\n", gnobj, mtot);
    *mtotp = mtot;
}


int SPH_need_update(const SPHbody *p)
{
    return 1;
}


Key_t SPHGetKey(const void *p)
{
    body t;
    VV(t.pos, = ((SPHbody *)p)->pos);
    return GETKEY(&t);
}


static void
AdjustBtab (SPHbody **SPHbtabp, int *nobj, int gnobj, double *outrmin, double *outrmax)
{
    SPHbody *btab = *SPHbtabp;
    SPHbody *p;
    Stk s;
    SPHbody *q;

    StkInitEz(&s);
    singlPrintf("Removing particles that do not overlap the output region: \n"); 
    singlPrintf("%g <= x <= %g, %g <= y <=%g, %g <= z <= %g \n",outrmin[0],outrmax[0],outrmin[1], outrmax[1], outrmin[2], outrmax[2]);
   

    for (p = btab; p < btab+*nobj; p++) {
	/* keep all particles inside reasonable volume of solution */

	if ( ( p->pos[0]+2*p->h >= outrmin[0] ) &&
	     ( p->pos[1]+2*p->h >= outrmin[1] ) &&
	     ( p->pos[2]+2*p->h >= outrmin[2] ) &&
	     ( p->pos[0]-2*p->h <= outrmax[0] ) &&
	     ( p->pos[1]-2*p->h <= outrmax[1] ) &&
	     ( p->pos[2]-2*p->h <= outrmax[2] ) ) { 

	    q = StkPush(&s, sizeof(SPHbody));
	    *q = *p;

	}
    }

    Free(btab);
    StkCrunch(&s);
    *nobj = StkSz(&s)/sizeof(SPHbody);
    btab = StkBase(&s);
    *SPHbtabp = Realloc(btab, *nobj * sizeof(SPHbody));
}



static void
AdjustBtab_Spherical (SPHbody **SPHbtabp, int *nobj, int gnobj, double innerbound, double outerbound)
{
    SPHbody *btab = *SPHbtabp;
    SPHbody *p;
    Stk s;
    SPHbody *q;
    double r;

    StkInitEz(&s);
    singlPrintf("Removing particles that do not overlap the output region: \n"); 
    singlPrintf("Inner bound: %g, Outer bound: %g", innerbound, outerbound);
    for (p = btab; p < btab+*nobj; p++) {
	/* keep all particles inside reasonable volume of solution */
      r=sqrt(p->pos[0]*p->pos[0]+p->pos[1]*p->pos[1]+p->pos[2]*p->pos[2]);
/* 	  singlPrintf("r=%g ",r); */
      if ( r > innerbound && r < outerbound)
	{
	  q = StkPush(&s, sizeof(SPHbody));
	  *q = *p;
	}/*  else { */
/* 	  singlPrintf("Particle %d discarded \n",p->ident); */
/* 	} */
    }

    Free(btab);
    StkCrunch(&s);
    *nobj = StkSz(&s)/sizeof(SPHbody);
    btab = StkBase(&s);
    *SPHbtabp = Realloc(btab, *nobj * sizeof(SPHbody));
/*     for (p = btab; p < btab+*nobj; p++) { */
/*       singlPrintf("kept: %g %g %g \n", p->pos[0], p->pos[1], p->pos[2]); */
/*     } */

}

static void
AdjustBtab_Rho (SPHbody **SPHbtabp, int *nobj, int gnobj, double rhocut)
{
    SPHbody *btab = *SPHbtabp;
    SPHbody *p;
    Stk s;
    SPHbody *q;
    double r;
    int discard;

    StkInitEz(&s);
    singlPrintf("Removing particles with rho < %g:\n", rhocut); 
    discard=0;
    for (p = btab; p < btab+*nobj; p++) {
      /* keep all particles with rho > rhocut */
      if ( p->rho > rhocut)
	{
	  q = StkPush(&s, sizeof(SPHbody));
	  *q = *p;
	} else {
	  discard++;
	  /* 	  singlPrintf("Particle %d discarded \n",p->ident); */
	}
    }

    Free(btab);
    StkCrunch(&s);
    *nobj = StkSz(&s)/sizeof(SPHbody);
    btab = StkBase(&s);
    *SPHbtabp = Realloc(btab, *nobj * sizeof(SPHbody));

}


static void SPHOutput(SPHbody *btab, int nobj, const char *outnamebase, int iter, int do_outputfloat)
{
    SPHbody *p;
    int i;
    sortresult_t outputsort;
    SPHoutbody *output_btab;
    SPHfloatoutbody *output_fbtab;
    int output_nobj = nobj;
    double tpos_out = tpos;
    double tvel_out = tvel; /* changed in Integrate() */
    double ke, pe, te;
    MPMY_Comm_request req;
    int output_gnobj;
    double output_z, output_h, output_R0;
    char outname[256];
    sprintf(outname, "%s", outnamebase);
    singlPrintf("Saving to %s\n", outname);
    pe = ke = te = 0.0;
    for (p = btab; p < btab+nobj; p++) {
	ke += (double)0.5 * p->mass * Dot(p->vel, p->vel);
	te += p->mass * p->u;
	pe += (double)0.5 * p->mass * p->phi;
    }
    output_btab = Malloc(output_nobj * sizeof(SPHoutbody));
    for(i=0; i<output_nobj; i++){
	output_btab[i].mass = btab[i].mass;
	VV(output_btab[i].pos, =  btab[i].pos);
	VV(output_btab[i].vel, =  btab[i].vel);
	output_btab[i].u =  btab[i].u;
	output_btab[i].h =  btab[i].h;
	output_btab[i].rho =  btab[i].rho;
  	output_btab[i].drho_dt = btab[i].drho_dt; 
	output_btab[i].udot =  btab[i].udot;
	output_btab[i].temp =  btab[i].temp;
#ifdef SPH_SAVE_ACC
	VV(output_btab[i].acc, =  btab[i].acc);
	VV(output_btab[i].acc_last, =  btab[i].acc_last);
	VV(output_btab[i].grav_acc, =  btab[i].grav_acc);
	output_btab[i].grav_mass =  btab[i].grav_mass;
	output_btab[i].phi =  btab[i].phi;
	output_btab[i].dt =  btab[i].dt;
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
    output_z = 0.0;
    output_h = 0.0;
    output_R0 = sysradius;


    if ( do_outputfloat == 0 || do_outputfloat == 2 ) {

      SDFwrite(outname, output_gnobj, 
	     output_nobj, output_btab, sizeof(SPHoutbody),
	     SPHOUTBODYDESC,
	     "npart", SDF_INT, output_gnobj,
	     "iter", SDF_INT, iter,
	     "dt", SDF_DOUBLE, dt,
	     "eps", SDF_DOUBLE, this_eps,
	     "Gnewt", SDF_DOUBLE, cosmo.GNewt,
	     "tolerance", SDF_DOUBLE, this_tol,
	     "frac_tolerance", SDF_DOUBLE, frac_tol,
	     "ndim", SDF_INT, NDIM,
	     "tpos", SDF_DOUBLE, tpos_out,
	     "tvel", SDF_DOUBLE, tvel_out,
	     "gamma", SDF_DOUBLE, Gamma,
	     "ke", SDF_DOUBLE, ke,
	     "pe", SDF_DOUBLE, pe,
	     "te", SDF_DOUBLE, te,
	     NULL);
    }

    /*NOTE: the float output file does not have grav_acc and grav_mass, 
     since the codes using float don't have this capability anyway.*/

    if ( do_outputfloat == 1 || do_outputfloat == 2 ) {
      output_fbtab = Malloc(output_nobj * sizeof(SPHfloatoutbody));
      for(i=0; i<output_nobj; i++){
	output_fbtab[i].mass = (float) output_btab[i].mass;
	VV(output_fbtab[i].pos, =  (float) output_btab[i].pos);
	VV(output_fbtab[i].vel, =  (float) output_btab[i].vel);
	output_fbtab[i].u =  (float) output_btab[i].u;
	output_fbtab[i].h =  (float) output_btab[i].h;
	output_fbtab[i].rho =  (float) output_btab[i].rho;
/*  	output_btab[i].drho_dt = btab[i].drho_dt; */
	output_fbtab[i].udot =  (float) output_btab[i].udot;
	output_fbtab[i].temp =  (float) output_btab[i].temp;
#ifdef SPH_SAVE_ACC
	VV(output_fbtab[i].acc, =  (float) output_btab[i].acc);
	VV(output_fbtab[i].acc_last, =  (float) output_btab[i].acc_last);
	output_fbtab[i].phi =  (float) output_btab[i].phi;
	output_fbtab[i].dt =  (float) output_btab[i].dt;
#endif
 	output_fbtab[i].nbrs = output_btab[i].nbrs;
	output_fbtab[i].ident = output_btab[i].ident;
      }
      singlPrintf("i made it.");
      sprintf(outname, "%s_float", outnamebase);
      SDFwrite(outname, output_gnobj, 
	     output_nobj, output_fbtab, sizeof(SPHfloatoutbody),
	     SPHFLOATOUTBODYDESC,
	     "npart", SDF_INT, output_gnobj,
	     "iter", SDF_INT, iter,
	     "dt", SDF_FLOAT, (float) dt,
	     "eps", SDF_FLOAT, (float) this_eps,
	     "Gnewt", SDF_FLOAT, (float) cosmo.GNewt,
	     "tolerance", SDF_FLOAT, (float) this_tol,
	     "frac_tolerance", SDF_FLOAT, (float) frac_tol,
	     "ndim", SDF_INT, NDIM,
	     "tpos", SDF_FLOAT, (float) tpos_out,
	     "tvel", SDF_FLOAT, (float) tvel_out,
	     "gamma", SDF_FLOAT, (float) Gamma,
	     "ke", SDF_FLOAT, (float) ke,
	     "pe", SDF_FLOAT, (float) pe,
	     "te", SDF_FLOAT, (float) te,
	     NULL);
      Free(output_fbtab);
    }
    Free(output_btab);
    singlPrintf("\nOutput done.\n");


/* #ifndef __DELTA__ */
/*     if (MPMY_Procnum() == 0) { */
/* 	char name[256]; */
/* 	sprintf(name, "%s.restart", outnamebase); */
/* 	if (unlink(name)) */
/* 	  Shout("unlink of %s failed, errno=%d\n", name, errno); */
/* 	if (symlink(outname, name)) */
/* 	  Shout("symlink of %s failed, errno=%d\n", outname, errno); */
/*     } */
/* #endif */
}


static void SPHOutputold(SPHbody *btab, int nobj, const char *outnamebase, int iter)
{
    SPHbody *p;
    int i;
    sortresult_t outputsort;
    SPHoutbody *output_btab;
    int output_nobj = nobj;
    double tpos_out = tpos;
    double tvel_out = tvel; /* changed in Integrate() */
    double ke, pe, te;
    MPMY_Comm_request req;
    int output_gnobj;
    double output_z, output_h, output_R0;
    char outname[256];
    sprintf(outname, "%s", outnamebase);
    singlPrintf("Saving to %s\n", outname);
    pe = ke = te = 0.0;
    for (p = btab; p < btab+nobj; p++) {
	ke += (double)0.5 * p->mass * Dot(p->vel, p->vel);
	te += p->mass * p->u;
	pe += (double)0.5 * p->mass * p->phi;
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
		      sizeof(SPHoutbody), (double) 0.1, 1, Realloc_f);
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
    output_z = 0.0;
    output_h = 0.0;
    output_R0 = sysradius;

    SDFwrite(outname, output_gnobj, 
	     output_nobj, output_btab, sizeof(SPHoutbody),
	     SPHOUTBODYDESC,
	     "npart", SDF_INT, output_gnobj,
	     "iter", SDF_INT, iter,
	     "dt", SDF_DOUBLE, dt,
	     "eps", SDF_DOUBLE, this_eps,
	     "Gnewt", SDF_DOUBLE, cosmo.GNewt,
	     "tolerance", SDF_DOUBLE, this_tol,
	     "frac_tolerance", SDF_DOUBLE, frac_tol,
	     "ndim", SDF_INT, NDIM,
	     "tpos", SDF_DOUBLE, tpos_out,
	     "tvel", SDF_DOUBLE, tvel_out,
	     "gamma", SDF_DOUBLE, Gamma,
	     "ke", SDF_DOUBLE, ke,
	     "pe", SDF_DOUBLE, pe,
	     "te", SDF_DOUBLE, te,
	     NULL);
    Free(output_btab);
    singlPrintf("\nOutput done.\n");
/* #ifndef __DELTA__ */
/*     if (MPMY_Procnum() == 0) { */
/* 	char name[256]; */
/* 	sprintf(name, "%s.restart", outnamebase); */
/* 	if (unlink(name)) */
/* 	  Shout("unlink of %s failed, errno=%d\n", name, errno); */
/* 	if (symlink(outname, name)) */
/* 	  Shout("symlink of %s failed, errno=%d\n", outname, errno); */
/*     } */
/* #endif */
}

void MySPHFixId(SPHbody *btab, int nobj, int gnobj){
    int start;
    int mynobj;
    int i;
    int *nobjproc;
    
    
    nobjproc=Malloc(sizeof(int)*MPMY_Nproc());
    for (i=0; i<MPMY_Nproc(); i++)
      {
	if (i == MPMY_Procnum()) 
	  {
	    nobjproc[i]=nobj;
	  } else {
	    nobjproc[i]=0;
	  }
      }
    MPMY_Combine(nobjproc, nobjproc, MPMY_Nproc(), MPMY_INT, MPMY_SUM);
    for (i=0; i<MPMY_Nproc(); i++) {
      singlPrintf("%d ", nobjproc[i]);
    }

    start=0;
    for (i=0; i<MPMY_Procnum(); i++) {
      start+=nobjproc[i];
    }

    for(i=0; i<nobj; i++){
      btab[i].ident = ((start+i) & ~ALL_FLAGS) | (btab[i].ident & ALL_FLAGS);
    }
}



