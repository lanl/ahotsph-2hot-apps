/* 
   Parallel WVT SPH Initial Conditions Setup
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
#include "initial.h"

#include "SDF.h"
#include "SDFwrite.h"
#include "singlio.h"
#include "sphinit.h"
#include "vop.h"
#include "decomp.h"
#include "consts.h"
#include "math.h"
#include <errno.h>

#define SQUARE(x) ((x)*(x))

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
static void AdjustBtab_Spherical (SPHbody **SPHbtabp, int *nobj, int gnobj, 
                                  double innerbound, double outerbound);
static void AdjustBtab_Rhozero (SPHbody **SPHbtabp, int *nobj, int gnobj, 
                                double rhozero);
void MySPHFixId(SPHbody *btab, int nobj, int gnobj);


static void SPHOutput(SPHbody *btab, int nobj, const char *outname, int iter, 
                      int do_floatoutput);


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
    double Zel_f; /* the 'f' factor for linearly growing modes,
		     used only in set_vel = 1/H*Ddot/D.  It's
		     very close to 1 (exactly?) for flat models. */
} cosmo;
static double dt=0.;

int main(int argc, char *argv[]) 
{
    int iter;
    SDF *csdfp, *sdfp;
    SPHbody *btab, *p;
    /*     void *decomp_info = NULL; */
    sortresult_t sortedbtab;
    tree_t SPHtree;
    double mtot;
    int num[NDIM];  /* uniform mesh for now */
    float rmin[NDIM], rmax[NDIM];
    double outrmin[NDIM], outrmax[NDIM];
    double sort_tol = 0.01;
    int i;
    int gnobj, nobj, targetnobj;
    int kernel_ncoef1, kernel_ncoef2;
    double kernel_coef1[MAXCOEF], kernel_coef2[MAXCOEF];
    char outnamebase[256];
    char outdir[256];
    double totvol;
    int do_externalstart;
    double tothvol;
    double outerbound=320., innerbound=-1.;
    int nloop, nhloop, nmassloop;
    double targetneighbors, targetmtot, nneighbors;
    char startfile[256];
    int ngood, ngood2;
    double npervol;
    int do_floatoutput;
    int keepcenterfixed;
    int do_hydrostatic; 
    int do_center, center_dual, center_sphfixed, center_posfixed;
    int special0_sphpoint, special1_sphpoint;
    double center_h, center_grav_mass;
    int do_eospolytrope;
    double kpolytrope, rhomin;

    int dogrid, do2d, dim;

    int inputoption, dimr, dimz, dimtheta;
    double cyl_minr, cyl_maxr, cyl_minz, cyl_maxz, cyl_mintheta, cyl_maxtheta;
    double cyl_center[3];
    
    int dimx, dimy;
    double cart_minx, cart_maxx, cart_miny, cart_maxy, cart_minz, cart_maxz;
    double cart_center[3];
    char cartfile_rho[256], cartfile_h[256];

    /*argv[1] = "/home/cellinge/WVT.dir/WVT-example.dir/initial.ctl";*/
    
    MPMY_Init(&argc, &argv);
    singlPrintf("Welcome to the SPH interpolator running on %d procs\n", 
		MPMY_Nproc());
    
    csdfp = initfiles(argc, argv);
    
    SDFgetintOrDefault(csdfp, "do_externalstart", &do_externalstart, 0);
    if (do_externalstart) 
	SDFgetstring(csdfp, "startfile", startfile, sizeof(startfile)); 
    SDFgetintOrDefault(csdfp, "nloop", &nloop, 50);
    SDFgetintOrDefault(csdfp, "nhloop", &nhloop, 50);
    SDFgetintOrDefault(csdfp, "nmassloop", &nmassloop, 50);
    SDFgetdoubleOrDefault(csdfp, "targetmtot", &targetmtot, -1.);
    SDFgetintOrDefault(csdfp, "keepcenterfixed", &keepcenterfixed, 0);
    SDFgetintOrDie(csdfp, "n_x", &(num[0]));
    SDFgetintOrDie(csdfp, "n_y", &(num[1]));
    SDFgetstring(csdfp, "outdir", outdir, sizeof(outdir));
    SDFgetintOrDie(csdfp, "n_z", &(num[2]));
    SDFgetintOrDefault  (csdfp, "targetnobj",  &targetnobj, 1000);
    SDFgetdoubleOrDefault(csdfp, "outerbound", &outerbound, 1e30);
    SDFgetdoubleOrDefault(csdfp, "innerbound", &innerbound, -1);
    SDFgetdoubleOrDefault(csdfp, "targetneighbors", &targetneighbors, 100.);
    SDFgetintOrDefault(csdfp, "do_floatoutput", &do_floatoutput, 0);
    SDFgetintOrDefault(csdfp, "do_hydrostatic", &do_hydrostatic, 0);
    SDFgetintOrDefault(csdfp, "do_eospolytrope", &do_eospolytrope, 0);
    if (do_eospolytrope) 
	SDFgetdoubleOrDie(csdfp, "kpolytrope", &kpolytrope);
    SDFgetdoubleOrDefault(csdfp, "rhomin", &rhomin, 1e-20);
    SDFgetdoubleOrDefault(csdfp, "Gamma", &Gamma, 1.66666666666667);

    SDFgetintOrDefault(csdfp, "dogrid", &dogrid, 0);
    SDFgetintOrDefault(csdfp, "do2d", &do2d, 0);
    if (do2d) dim=2; else dim=3;

    /* Input Method */
    SDFgetintOrDefault(csdfp, "inputoption", &inputoption, 1);
    WVT_setinputoption(inputoption);

    if (inputoption == 2) { /* Cylindrical Grid */
	SDFgetintOrDie(csdfp, "dimr", &dimr);
	SDFgetintOrDie(csdfp, "dimz", &dimz);
	SDFgetintOrDie(csdfp, "dimtheta", &dimtheta);
	
	SDFgetdoubleOrDefault(csdfp, "cyl_minr", &cyl_minr, 0.);
	SDFgetdoubleOrDefault(csdfp, "cyl_maxr", &cyl_maxr, 1.);
	SDFgetdoubleOrDefault(csdfp, "cyl_minz", &cyl_minz, 0.);
	SDFgetdoubleOrDefault(csdfp, "cyl_maxz", &cyl_maxz, 1.);
	SDFgetdoubleOrDefault(csdfp, "cyl_mintheta", &cyl_mintheta, 0.);
	SDFgetdoubleOrDefault(csdfp, "cyl_maxtheta", &cyl_maxtheta, 6.27);
	
	SDFgetdoubleOrDefault(csdfp, "cyl_xcenter", &(cyl_center[0]), 0.);
	SDFgetdoubleOrDefault(csdfp, "cyl_ycenter", &(cyl_center[1]), 0.);
	SDFgetdoubleOrDefault(csdfp, "cyl_zcenter", &(cyl_center[2]), 0.);
	
	init_cylindricalgrid(dimr, dimz, dimtheta, cyl_minr, cyl_maxr, 
			     cyl_minz, cyl_maxz,
			     cyl_mintheta, cyl_maxtheta, cyl_center);
    }
    
    if (inputoption == 4) { /* Cartesion Grid */
	SDFgetintOrDie(csdfp, "cart_dimx", &dimx);
	SDFgetintOrDie(csdfp, "cart_dimy", &dimy);
	SDFgetintOrDie(csdfp, "cart_dimz", &dimz);
	
	SDFgetdoubleOrDefault(csdfp, "cart_minx", &cart_minx, -1.);
	SDFgetdoubleOrDefault(csdfp, "cart_maxx", &cart_maxx, 1.0);
	SDFgetdoubleOrDefault(csdfp, "cart_miny", &cart_miny, -1.);
	SDFgetdoubleOrDefault(csdfp, "cart_maxy", &cart_maxy, 1.0);
	SDFgetdoubleOrDefault(csdfp, "cart_minz", &cart_minz, -1.);
	SDFgetdoubleOrDefault(csdfp, "cart_maxz", &cart_maxz, 1.0);
	
	SDFgetdoubleOrDefault(csdfp, "cart_xcenter", &(cart_center[0]), 0.);
	SDFgetdoubleOrDefault(csdfp, "cart_ycenter", &(cart_center[1]), 0.);
	SDFgetdoubleOrDefault(csdfp, "cart_zcenter", &(cart_center[2]), 0.);
	
	SDFgetstring(csdfp,"cartfile_rho",cartfile_rho, sizeof(cartfile_rho)); 
	SDFgetstring(csdfp, "cartfile_h", cartfile_h, sizeof(cartfile_h)); 
	
	init_cartesiangrid(dimx, dimy, dimz, cart_minx, cart_maxx, 
			   cart_miny, cart_maxy,
			   cart_minz, cart_maxz, cart_center, 
			   cartfile_rho, cartfile_h);
    }
    
    /* Center */
    SDFgetintOrDefault(csdfp, "do_center", &do_center, 0);
    SDFgetintOrDefault(csdfp, "center_dual", &center_dual, 0);
    SDFgetintOrDefault(csdfp, "center_sphfixed", &center_sphfixed, 0);
    SDFgetintOrDefault(csdfp, "center_posfixed", &center_posfixed, 0);
    SDFgetintOrDefault(csdfp, "special0_sphpoint", &special0_sphpoint, 0);
    SDFgetintOrDefault(csdfp, "special1_sphpoint", &special1_sphpoint, 0);
    SDFgetdoubleOrDefault(csdfp, "center_h", &center_h, 0.1);
    SDFgetdoubleOrDefault(csdfp, "center_grav_mass", &center_grav_mass, 
			  0.391973);
    if (do_center) keepcenterfixed=1;
    
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
    singlPrintf("Kernel: %g %g %g %g %g %g %g %g \n", 
		kernel_coef1[0], kernel_coef1[1], kernel_coef1[2], 
		kernel_coef1[3], kernel_coef2[0], kernel_coef2[1], 
		kernel_coef2[2], kernel_coef2[3]);

    if (SDFhasname("outrmin", csdfp)) {
	if (SDFseekrdvecs(csdfp, "outrmin", 0, 3, 
			  outrmin, 0, NULL))
	    Error("SDFread outrmin failed\n");
    } else {
	outrmin[0]=-1e30 ;
	outrmin[1]=-1e30 ;
	outrmin[2]=-1e30 ;
    }
    if (SDFhasname("outrmax", csdfp)) {
	if (SDFseekrdvecs(csdfp, "outrmax", 0, 3, 
			  outrmax, 0, NULL))
	    Error("SDFread outrmin failed\n");
    } else {
	outrmax[0]=1e30 ;
	outrmax[1]=1e30 ;
	outrmax[2]=1e30 ;
    }
    
    singlPrintf("outrmin: %g %g %g \n",outrmin[0],outrmin[1],outrmin[2]);
    singlPrintf("outrmax: %g %g %g \n",outrmax[0],outrmax[1],outrmax[2]);
    SDFclose(csdfp);
    singlPrintf("number of arguments: %d %s", argc, argv[2]);
    
    ClearEnabledTimers();
    ClearEnabledCounters();
    StartTimer(&StepTotWC);
    StartTimer(&StepTot);
    
    totvol=pow(outerbound, dim);
    if (do_externalstart) {
	sdfp = SPHReadf(startfile, &btab, &gnobj, &nobj); 
	SDFgetdoubleOrDefault(sdfp, "tpos",  &tpos, (double)0.0);
	SDFgetintOrDefault  (sdfp, "iter",  &iter, 0);
	SDFclose(sdfp);
	singlPrintf("SPHReadf done");
	iter=0;
	
	AdjustBtab_Spherical((SPHbody **)&btab, &nobj, gnobj, innerbound, 
			     outerbound);    
	MPMY_Combine(&nobj, &gnobj, 1, MPMY_INT, MPMY_SUM); 
	p = Malloc(sizeof(SPHoutbody));
	WVT_hofpos(p,-1,totvol, &tothvol, dim);
	WVT_hofpos(p,1,totvol, &tothvol, dim);
	for (p = btab; p < btab+nobj; p++)
	    if (p->ident & DUAL_FLAG) 
		p->ident=p->ident | POSFIXED_FLAG;
   	MPMY_Combine(&nobj, &gnobj, 1, MPMY_INT, MPMY_SUM); 
	MySPHFixId(btab, nobj, gnobj); 
        sprintf(outnamebase, "%swvt.0000", outdir);
	iter=0;
	SPHOutput(btab, nobj, outnamebase, iter, do_floatoutput);
    } else { 
	/* initialize hofpos first, so you only read in data once */
	WVT_hofpos(btab,-1,totvol, &tothvol, dim);
	if (do2d) {
	    num[2]=1;
	    outrmax[2]=1e-10;
	    outrmin[2]=-1e-10;
	}
	singlPrintf("Cube num: %d %d %d \n", num[0], num[1], num[2]);
	 
	if (dogrid == 0) 
	    WVTInitProbdist(&btab, &gnobj, &nobj, outrmin, outrmax, 
			    targetnobj, totvol, outerbound, innerbound, num,
			    dim); 
	else if (dogrid == 1) 	    
	    WVTInitCube(&btab, &gnobj, &nobj, outrmin, outrmax, num, dim);   
 	else if (dogrid == 2) 
	    WVTInitHCP(&btab, &gnobj, &nobj, outrmin, outrmax, num, dim);  
 	else if (dogrid == 3) 
	    WVTInitCCP(&btab, &gnobj, &nobj, outrmin, outrmax, num, dim);  
 	else if (dogrid == 4) 
	    WVTInitProbdistlr(&btab, &gnobj, &nobj, outrmin, outrmax, 
			    targetnobj, totvol, outerbound, innerbound, num,
			    dim); 
	
	singlPrintf("done");	
	AdjustBtab_Spherical((SPHbody **)&btab, &nobj, gnobj, innerbound, 
			     outerbound);    
	
	if (do2d) for (p=btab; p<btab+nobj; p++) p->pos[2]=0.;
	
   	MPMY_Combine(&nobj, &gnobj, 1, MPMY_INT, MPMY_SUM); 
	MySPHFixId(btab, nobj, gnobj); 
        sprintf(outnamebase, "%swvt.0000", outdir);
	iter=0;
	SPHOutput(btab, nobj, outnamebase, iter, do_floatoutput);
    } 
    
    SPHSanityCheck(btab, nobj, gnobj, &mtot);
    SPHFindBbox(btab, nobj, rmin, rmax);
    singlPrintf("rmin: %g %g %g \n",rmin[0],rmin[1],rmin[2]);
    singlPrintf("rmax: %g %g %g \n",rmax[0],rmax[1],rmax[2]);
    SPHFixNterms(btab, nobj);
    
    if (keepcenterfixed && MPMY_Procnum() == 0 && !do_externalstart) {
	btab[0].pos[0]=0.;
	btab[0].pos[1]=0.;
	btab[0].pos[2]=0.;
        btab[0].vel[0]=0.;
        btab[0].vel[1]=0.;
        btab[0].vel[2]=0.;
	btab[0].h=0.;    /* h=0 should be equivalent of keeping pos fixed */
	btab[0].ident=btab[0].ident | POSFIXED_FLAG;
	btab[0].ident=btab[0].ident | DUAL_FLAG;
    }
    
    
    SetupTree(&SPHtree, NDIM, sizeof(SPHbody), sizeof(SPHcell), 0, 0,
	      SPHTBODYSZ, sizeof(SPHcofmdata), 
	      (pq_keyproto)SPHGetKeyFromStruct, (pq_wgtproto)SPHGetCost,
	      SPHCofmFromDaugh, (cellfromcofm_t)SPHCellFromCofm, SPHCellSz);
    
    SPH_setup(NDIM, kernel_ncoef1, kernel_coef1, kernel_ncoef2, kernel_coef2);
    WVT_setup(NDIM, kernel_ncoef1, kernel_coef1, kernel_ncoef2, kernel_coef2);

    SPHSetupCofm(&SPHtree);

    for(i = 1; i < nloop; ++i) {

	singlPrintf("ITER %d", i);
	
	MPMY_Combine(&nobj, &gnobj, 1, MPMY_INT, MPMY_SUM);
	singlPrintf("rmin: %g %g %g \n",rmin[0],rmin[1],rmin[2]);
	singlPrintf("rmax: %g %g %g \n",rmax[0],rmax[1],rmax[2]);
	
	singlPrintf("BuildTree\n");
	StartTimer(&BuildTot);
	
	nneighbors=targetneighbors*(1.+0.*(0.5*pow(2,dim)-1.)*
				    (1.-((double) i)/((double) nloop-1.)));

	/* Assign h according to the desired distribution, see wvt.c */
	totvol=pow(outerbound, dim)*nneighbors;
	if (innerbound > 0) 
	    totvol-=pow(innerbound,dim)*nneighbors;
 	for (p=btab; p<btab+nobj; p++) {
	    p->acc[0]=0.;
	    p->acc[1]=0.;
	    p->acc[2]=0.;
	}
	
	tothvol=-1.; /* Set to -1 so it is computed again */
	WVT_hofpos(btab,nobj,totvol, &tothvol, dim);
	
	SphericalGhosts(&btab, &nobj, outerbound, innerbound,
			&tothvol, totvol, nneighbors, dim);
	/* 	BoxGhosts(&btab, &nobj, outrmin, outrmax, &tothvol, totvol); */

	SPHFindBbox(btab, nobj, rmin, rmax);
	
	/* Initialize these variables, since they store the particle 
	   separations in this code, instead of physical properties */
	for (p=btab; p<btab+nobj; p++)
	    {
		p->udot=1e30;
		p->rho_est=1e30;
		p->vsound=1e30;
		p->temp=1e30;
		p->drho_dt=1e30;
		p->rho=0;
		p->acc[0]=0.;
		p->acc[1]=0.;
		p->acc[2]=0.;
		p->nbrs=0;
	  }

	singlPrintf("nobj:%d\n", nobj);
	MPMY_Combine(&nobj, &gnobj, 1, MPMY_INT, MPMY_SUM);
	SPHFixNterms(btab, nobj);
	pqsortsetup(&sortedbtab, btab, nobj, sizeof(SPHbody), sort_tol,
		    Realloc_f);
        sortedbtab.method = 0;
        sortedbtab.decomp_cycles = 1;
        sortedbtab.first_cycle = 0;
	singlPrintf("nobj:%d\n", nobj);
	SPHFixKeys(btab, nobj, SPHGetKey);
	singlPrintf("nobj:%d\n", nobj);
 	Warning("after fixkeys Proc:%d Nproc:%d Nobj:%d\n", 
		MPMY_Procnum(), MPMY_Nproc(), nobj );  
	
	BuildTree(&SPHtree, &sortedbtab);
	singlPrintf("nobj:%d\n", nobj);
	btab = sortedbtab.data;
	nobj = sortedbtab.nobj;
	/* 	decomp_info = SaveDecomp(); */
	/* 	SetDecomp(decomp_info); */
	StopTimer(&BuildTot);
	
        SetSPH(0., 0., 0., 0., Gamma, gnobj,
	       macWVT, nbrMAC);
	SetWVT(0., 0., 0., 0., Gamma, gnobj,
	       macWVT, nbrMAC);
	
	singlPrintf("Walk\n");
	
	
	WalkInit(&SPHtree, &SPHtree, sizeof(SinkSPH), WalkInitSrc,
                 (macv_t)WVTgate, (inherit_t)InheritWVT);
	singlPrintf("Walkinit done\n");
	WalkNT(&SPHtree);
	singlPrintf("WalkNT done\n");
	WalkTerminate();
	singlPrintf("WalkTerminate done\n");
	
	FreeTree(&SPHtree);
	singlPrintf("FreeTree done\n");
	
	StopTimer(&StepTot);
	StopTimer(&StepTotWC);
	
	OutputTimers(singlPrintf);

	/* Advance to the next "time" step */
	if (keepcenterfixed)
	    for (p = btab; p < btab+nobj; p++)
		if (p->ident & POSFIXED_FLAG)
		    p->h=0.;
	
	WVTupdate(btab, nobj, i, nloop, dim, nneighbors);	
	if (do2d) for (p=btab; p<btab+nobj; p++) p->pos[2]=0.;
	tothvol=-1.; /* Set to -1 so it is computed again */
	WVT_hofpos(btab,nobj,totvol, &tothvol, dim);
	
        if (i < 10) sprintf(outnamebase, "%swvt.000%d", outdir, i);
        else if (i < 100) sprintf(outnamebase, "%swvt.00%d", outdir, i);
        else if (i < 1000) sprintf(outnamebase, "%swvt.0%d", outdir, i);
        else sprintf(outnamebase, "%swvt.%d", outdir, i);
	SPHOutput(btab, nobj, outnamebase, iter, do_floatoutput);
	
/* 	for (p=btab; p<btab+nobj; p++)  */
/* 	    if (!(p->ident & (1<<28))) singlPrintf("%d ", p->nbrs); */
	
	singlPrintf("Nobj with ghosts: %d", nobj);
	RemoveGhosts(&btab, &nobj);
	singlPrintf("After ghosts removed:nobj=%d\n", nobj);
    }
    
    tothvol=-1.;
    WVT_hofpos(btab,nobj,totvol, &tothvol, dim);      
    WVT_hofpos(btab,nobj,totvol, &tothvol, dim);      
    
    
    if (nhloop > 0) {
	tothvol=-1.;
	WVT_hofpos(btab,nobj,totvol, &tothvol, dim);
    }


    /* Compute h to ensure a constant number of neighbors */
    for(i = 1; i <= nhloop; i++) {
	for (p=btab; p<btab+nobj; p++) p->nbrs=0;
	
	SphericalGhosts(&btab, &nobj, outerbound, innerbound,
			&tothvol, totvol, targetneighbors, dim);
	/* 	BoxGhosts(&btab, &nobj, outrmin, outrmax, &tothvol, totvol); */

	SPHFindBbox(btab, nobj, rmin, rmax);
	
	MPMY_Combine(&nobj, &gnobj, 1, MPMY_INT, MPMY_SUM);
	SPHFixNterms(btab, nobj);	      	
	pqsortsetup(&sortedbtab, btab, nobj, sizeof(SPHbody), sort_tol,
		    Realloc_f);
        sortedbtab.method = 0;
        sortedbtab.decomp_cycles = 1;
        sortedbtab.first_cycle = 0;
	SPHFixKeys(btab, nobj, SPHGetKey);
	BuildTree(&SPHtree, &sortedbtab);
	btab = sortedbtab.data;
	nobj = sortedbtab.nobj;
	
        SetSPH(0., 0., 0., 0., Gamma, gnobj,
	       macConstNeigh, nbrMAC);
	
	singlPrintf("Walk\n");	
	WalkInit(&SPHtree, &SPHtree, sizeof(SinkSPH), WalkInitSrc,
                 (macv_t)SPHgate, (inherit_t)InheritSPH);
	singlPrintf("Walkinit done\n");
	WalkNT(&SPHtree);
	singlPrintf("WalkNT done\n");
	WalkTerminate();
	singlPrintf("WalkTerminate done\n");
	
	FreeTree(&SPHtree);
	singlPrintf("FreeTree done\n");
	
        if (i < 10) sprintf(outnamebase, "%swvtneigh.000%d", outdir, i);
        else if (i < 100) sprintf(outnamebase, "%swvtneigh.00%d", outdir, i);
        else if (i < 1000) sprintf(outnamebase, "%swvtneigh.0%d", outdir, i);
        else sprintf(outnamebase, "%swvtneigh.%d", outdir, i);
	SPHOutput(btab, nobj, outnamebase, iter, do_floatoutput);
	
	/* Calculate the fraction of particles that meets the neighbor 
	   criterion +- 5%. If the fraction is > 95%, stop.*/
	ngood=0;
	for (p=btab; p<btab+nobj; p++) {
	    if (p->nbrs >(int) (0.95*targetneighbors)
		&& p->nbrs <(int) (1.05*targetneighbors))
		ngood++;
	}
	MPMY_Combine(&ngood, &ngood, 1, MPMY_INT, MPMY_SUM);
	if ( (double) ngood/gnobj > .999) {
	    singlPrintf("Smoothing lengths are good enough.");
	    break;
	}

	for (p=btab; p<btab+nobj; p++) {
	    if (p->nbrs <(int) (0.975*targetneighbors) 
		|| p->nbrs >(int) (1.025*targetneighbors))
		p->h*=pow((targetneighbors+1.)/
			  ((double) p->nbrs+1.),.33333333);
	}
	singlPrintf("Nobj with ghosts: %d", nobj);
	RemoveGhosts(&btab, &nobj);
	singlPrintf("After ghosts removed:nobj=%d\n", nobj);
    }


    
    /* Compute the mass per particle, so that interpolation yields rho */
    /* THIS WON'T WORK VERY WELL IN THE CENTER IF YOU ARE RESOLUION LIMITED  
       AND MAY BE A RUNAWAY PROCESS. TO AVOID THIS YOU WOULD HAVE TO 
       WRITE A NEW MAC THAT UPDATES ALL MASSES WITHIN 2H OF A RHO THAT HAS TO 
       BE CORRECTED. TOO MUCH WORK FOR NOW...*/
    if (inputoption == 1) SPHofpos(btab,nobj); 
    
    AdjustBtab_Rhozero((SPHbody **)&btab, &nobj, gnobj, 1e-20);

    if (nmassloop > 0) { 
	mtot=0.;
	for (p=btab; p<btab+nobj; p++) {
	    if (do2d) npervol=((double)targetneighbors)/(3.1459*p->h*p->h); 
	    else npervol=((double)targetneighbors) / 
		     (4./3.*3.1459*p->h*p->h*p->h*8.);
	    p->mass=p->rho/npervol;
	    p->du=p->rho;
	    mtot+=p->mass;
	}
	MPMY_Combine(&mtot, &mtot, 1, MPMY_DOUBLE, MPMY_SUM);
	if (targetmtot > 0) {
	    for (p=btab; p<btab+nobj; p++) {
		p->mass*=targetmtot/mtot;
	    }
	}
    }
    
    for(i = 1; i < nmassloop; ++i) {
	for (p=btab; p<btab+nobj; p++) p->nbrs=0; 
        for (p=btab; p<btab+nobj; p++) p->rho=0.;
	
	SPHFindBbox(btab, nobj, rmin, rmax);
	
	MPMY_Combine(&nobj, &gnobj, 1, MPMY_INT, MPMY_SUM);
	
	SPHFixNterms(btab, nobj);
	pqsortsetup(&sortedbtab, btab, nobj, sizeof(SPHbody), sort_tol,
		    Realloc_f);
        sortedbtab.method = 0;
        sortedbtab.decomp_cycles = 1;
        sortedbtab.first_cycle = 0;
	SPHFixKeys(btab, nobj, SPHGetKey);
	BuildTree(&SPHtree, &sortedbtab);
	btab = sortedbtab.data;
	nobj = sortedbtab.nobj;
	
        SetSPH(0., 0., 0., 0., Gamma, gnobj,
	       macRho, nbrMAC);
	
	singlPrintf("Walk\n");	
	WalkInit(&SPHtree, &SPHtree, sizeof(SinkSPH), WalkInitSrc,
                 (macv_t)SPHgate, (inherit_t)InheritSPH);
	singlPrintf("Walkinit done\n");
	WalkNT(&SPHtree);
	singlPrintf("WalkNT done\n");
	WalkTerminate();
	singlPrintf("WalkTerminate done\n");
	
	FreeTree(&SPHtree);
	singlPrintf("FreeTree done\n");
	
        for (p=btab; p<btab+nobj; p++) {
	    p->rho += 1./3.141592653589793238462 *p->mass/(p->h * p->h * p->h);
	}
	
        if (i < 10) sprintf(outnamebase, "%swvtsph.000%d", outdir, i);
        else if (i < 100) sprintf(outnamebase, "%swvtsph.00%d", outdir, i);
        else if (i < 1000) sprintf(outnamebase, "%swvtsph.0%d", outdir, i);
        else sprintf(outnamebase, "%swvtsph.%d", outdir, i);
	SPHOutput(btab, nobj, outnamebase, iter, do_floatoutput);
	
	/* Calculate the fraction of particles that meets reproduce rho within 
	   +- 1%. If the fraction is > 95%, stop.*/
	ngood=0;
	ngood2=0;
	mtot=0.;
	for (p=btab; p<btab+nobj; p++) {
	    mtot+=p->mass;
	    if (p->rho > 0.99*p->du && p->rho < 1.01*p->du)
		ngood++;
	    if (p->rho > 0.95*p->du && p->rho < 1.05*p->du)
		ngood2++;
	}
	MPMY_Combine(&ngood, &ngood, 1, MPMY_INT, MPMY_SUM);
	MPMY_Combine(&ngood2, &ngood2, 1, MPMY_INT, MPMY_SUM);
	MPMY_Combine(&mtot, &mtot, 1, MPMY_DOUBLE, MPMY_SUM);
	singlPrintf("%d of %d densities within 0.01 of target.\n",
		    ngood, gnobj);
	singlPrintf("%d of %d densities within 0.05 of target.\n", 
		    ngood2, gnobj);
	singlPrintf("Total Mass: %g.\n", mtot);

	if ( (double) ngood/gnobj > .999 || i == nmassloop ) {
	    singlPrintf("Masses are good enough.");
	    break;
	}

	for (p=btab; p<btab+nobj; p++) {
	    p->mass=p->mass*(1+p->du/p->rho)/2.;
	}
    }

    /* Fill the grav_mass column */
    for (p=btab; p<btab+nobj; p++) {
	p->grav_mass=p->mass;
    }

    if (do_eospolytrope) 
	for (p=btab; p<btab+nobj; p++) 
	    p->temp=kpolytrope;
    
    /* Make sure the pressure gradient and rho is correct, i.e. derive 
       internal energy from pr, rho and equation of state for ideal gas */
    if (do_hydrostatic)
	for (p=btab; p<btab+nobj; p++) p->u=p->pr/(Gamma-1.)/p->rho;
    

    /*     And now set the properties of the central particle */
    if (do_center) {
	singlPrintf("Setting central particle properties.\n");
	for (p = btab; p < btab+nobj; p++)
	    if (p->ident & DUAL_FLAG) {
		p->h=center_h;
		p->grav_mass=center_grav_mass;
		/*p->mass=1e-5;  */
		/*p->rho=p->mass*kernel_coef1[0]/3.14159/p->h/p->h/p->h; */
		p->rho_est=p->rho;
		/*p->u=7e-07; */
		if (center_posfixed) p->ident=p->ident | POSFIXED_FLAG; 
		else p->ident=p->ident & ~POSFIXED_FLAG; 
		
		if (center_sphfixed) p->ident=p->ident | SPHFIXED_FLAG; 
		else p->ident=p->ident & ~SPHFIXED_FLAG; 
		
		if (center_dual) p->ident=p->ident | DUAL_FLAG;
		else p->ident=p->ident & ~DUAL_FLAG;
		
		if (special0_sphpoint) p->ident=p->ident | SPECIAL0_FLAG;
		else p->ident=p->ident & ~SPECIAL0_FLAG;
		
		if (special1_sphpoint) p->ident=p->ident | SPECIAL1_FLAG;
		else p->ident=p->ident & ~SPECIAL1_FLAG;
		
		singlPrintf("\n");
	    }
    } else {
	for (p = btab; p < btab+nobj; p++) {
	    p->ident=p->ident & ~DUAL_FLAG;
	    p->ident=p->ident & ~POSFIXED_FLAG;
	}
    }
    
    /* Compute SPH properties (rho, u, press, vel) */
    /* I have to have similar options for the other input methods!!*/
    /*     if (inputoption == 1) SPHofpos(btab,nobj);  */
    
    /*interp_cylindricalgrid(btab, nobj, dimr, dimz, dimtheta, minr, maxr, */
    /*                            minz, maxz, mintheta, maxtheta); */
    AdjustBtab_Spherical((SPHbody **)&btab, &nobj, gnobj, innerbound, 
			 outerbound*1.05);    

    AdjustBtab_Rhozero((SPHbody **)&btab, &nobj, gnobj, rhomin);
    
    sprintf(outnamebase, "%swvtfinal.sdf", outdir);
    SPHOutput(btab, nobj, outnamebase, iter, do_floatoutput);
    
    singlPrintf("Bye!\n");

    MPMY_Finalize();
    
    return 0;
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
AdjustBtab_Spherical (SPHbody **SPHbtabp, int *nobj, int gnobj, 
		      double innerbound, double outerbound)
{
    SPHbody *btab = *SPHbtabp;
    SPHbody *p;
    Stk s;
    SPHbody *q;
    double r;
    
    StkInitEz(&s);
    singlPrintf("Removing particles that do not overlap the output region:\n");
    singlPrintf("Inner bound: %g, Outer bound: %g", innerbound, outerbound);
    for (p = btab; p < btab+*nobj; p++) {
	/* keep all particles inside reasonable volume of solution */
	r=sqrt(p->pos[0]*p->pos[0]+p->pos[1]*p->pos[1]+p->pos[2]*p->pos[2]);
	/* 	  singlPrintf("r=%g ",r); */
	if ( r > innerbound && r < outerbound) {
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
AdjustBtab_Rhozero (SPHbody **SPHbtabp, int *nobj, int gnobj, double rhozero)
{
    SPHbody *btab = *SPHbtabp;
    SPHbody *p;
    Stk s;
    SPHbody *q;
    
    StkInitEz(&s);
    singlPrintf("Removing particles with rho < %g \n", rhozero);
    for (p = btab; p < btab+*nobj; p++) {
	/* keep all particles inside the computational domain (rho>rhozero) */
	if ( p->rho > rhozero)
	    {
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


static void SPHOutput(SPHbody *btab, int nobj, const char *outnamebase, 
		      int iter, int do_outputfloat)
{
    int i;
    sortresult_t outputsort;
    SPHoutbody *output_btab;
    SPHfloatoutbody *output_fbtab;
    int output_nobj = nobj;
    double tpos_out = tpos;
    double tvel_out = tvel; /* changed in Integrate() */
    MPMY_Comm_request req;
    int output_gnobj;
    char outname[256];
    sprintf(outname, "%s", outnamebase);
    singlPrintf("Saving to %s\n", outname);
    
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
/* 	VV(output_btab[i].acc_last, =  btab[i].acc_last); */
/* 	VV(output_btab[i].grav_acc, =  btab[i].grav_acc); */
 	output_btab[i].grav_mass =  btab[i].grav_mass; 
/* 	output_btab[i].phi =  btab[i].phi; */
/* 	output_btab[i].dt =  btab[i].dt; */
#endif
 	output_btab[i].nbrs = btab[i].nbrs;
	output_btab[i].ident = btab[i].ident;
    }
    /*     Msg("output", ("Doing output of %d bodies\n", output_nobj)); */
    Msgf(("Doing output of %d bodies\n", output_nobj));
    singlPrintf("Trying to sort output\n");
    pqsortsetup_order(&outputsort, output_btab, output_nobj,
		      sizeof(SPHoutbody), 0.1F, 1, Realloc_f);
    outputsort.method = 1;
    outputsort.decomp_cycles = 1;
    outputsort.first_cycle = 0;
    output_btab = pqsort(&outputsort, UnityCost, (pq_keyproto)SPHOutIdentKey);
    output_nobj = outputsort.nobj;
    
    
    /*     Msg("output", ("After pqsort, %d outbodies\n", output_nobj)); */
    Msgf(("After pqsort, %d outbodies\n", output_nobj));
    MPMY_ICombine_Init(&req);
    MPMY_ICombine(&output_nobj, &output_gnobj, 1, MPMY_INT, MPMY_SUM, req);
    MPMY_ICombine_Wait(req);

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



void MySPHFixId(SPHbody *btab, int nobj, int gnobj){
    int start;
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

    free(nobjproc);
}




