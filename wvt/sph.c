/*
 * © 2026. Triad National Security, LLC. All rights reserved.
 * This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration. All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration. The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare. derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.
 */

#include <math.h>
#include <stdlib.h>
#include "stk.h"
#include "Malloc.h"


#include "physics_sph.h"
#include "Msgs.h"
#include "timers.h"
#include "vop.h"
#include "fastflpt.h"
#include "mpmy.h"
#include "gc.h"
#include "singlio.h"

#define PI 3.141592653589793238462
#ifndef M_1_PI
#define	M_1_PI 0.31830988618379067154
#endif
/*#define NKERNEL_TABLE 40000*/
#define NKERNEL_TABLE 40000
#define MAX_INDEX (NKERNEL_TABLE+2)

#define NO_UPDATE 2

static double dvtable; /* == 0.0001 ... */
static double invdvtable; /* == 10000.0 ... */
static double cnormk;
static double wij[MAX_INDEX];
static double grwij[MAX_INDEX];
static double fmass[MAX_INDEX];
static double fpoten[MAX_INDEX];
static double Gamma = (double)(4.0/3.0);  /* Yikes; remember this */
static double alpha = (double)1.0;
static double beta = (double)2.5;
static double epsil = (double)1e-2;
static double heatf1 = (double)1.0;
static int ndim;
static int Nobj;
static int add_offset;
static double offset[NDIM];
static double voffset[NDIM];
static void (*bodyfunc)(SinkSPH *sink, hcell **src_vec, int *res, int n);
static void (*cellfunc)(SinkSPH *sink, hcell **src_vec, int *res, int n);

extern int do_diffusion;

void
SetSPHOffset(double *off, double *voff)
{
    VV(offset, = off);
    VV(voffset, = voff);
    add_offset = 1;
}

void
UnSetSPHOffset(void)
{
    VS(offset, = 0.0);
    VS(voffset, = 0.0);
    add_offset = 0;
}

/* Note that InheritSPH in WVT is different from the normal one!!!! */
void 
InheritSPH(const SinkSPH *from, SinkSPH *to, hcell *pp)
{
    if( to == NULL ){
	SPHbody *bp = pp->ptr;
	if (from->isbody == NO_UPDATE)
	  return;
	/* Must accumulate for periodic BC to work */
	/* Must initialize to zero appropriately */
	bp->rho += from->rho;
	bp->du = from->du;       /* Why am I doing this? */
	bp->du_r += from->du_r;   /* Do I need to do anything else? */
	bp->drho_dt += from->drho_dt;
	bp->udot += from->udot;
	bp->nbrs += from->nbrs;
	VV(bp->acc, += from->M1);
	VV(bp->lvel, += from->lvel);
	bp->nterms += from->nterms;
	bp->min_nbr_dt = from->min_nbr_dt;   /* set dt to min of nbrs dt */
	if (from->interactions != Nobj)
	    Error("Ninteract is %d, should be %d\n", from->interactions, Nobj);
	return;
    }

    if (Sub_Flags(pp)) {
	/* Stuff needed for cell-cell neighbor evaluation */
	SPHcell *cp = pp->ptr;
	VV(to->pos, = cp->pos);
	to->extent = cp->bmax + cp->lap;
	to->isbody = 0;
    } else {
	/* Stuff needed for above, plus physics info */
	SPHbody *bp = pp->ptr;
	if (!SPH_need_update(bp)) {
	    to->isbody = NO_UPDATE;
	    return;
	}
	VV(to->pos, = bp->pos);
	to->extent = bp->h;
	to->h = bp->h;
	to->isbody = 1;
	VV(to->vel, = bp->vel);
	to->pr = bp->pr;
	to->rho_est = bp->rho_est;
	to->mass = bp->mass;
	to->vsound = bp->vsound;
	to->u = bp->u;
	to->rho = (double)0.0;
	to->drho_dt = (double)0.0;
	to->udot = (double)0.0;
	VS(to->lvel, = (double)0.0);
	to->nterms = 1;
	to->nbrs = 0;
	VS(to->M1, = (double)0.0);
	to->min_nbr_dt = 1e30;
	/* Diffusion quantities */
	/* Can some of these be set to zero here?  Like udot above? */
	to->temp = bp->temp;
	to->du = bp->du;
	to->u_r = bp->u_r;
	to->du_r = bp->du_r;  /* Or = (double)0.0; ? */
/* 	to->D = bp->D; */
    }

    if (add_offset) {
	VV(to->pos, += offset);
	VV(to->vel, += voffset);
    }

    if (from) {
	to->interactions = from->interactions;
    } else {
	to->interactions = 0;
    }

}

void
SPHgate(SinkSPH *sink, hcell **src_vec, int *flags, int *result, int n)
{
    int i;

    if (sink->isbody == NO_UPDATE) {
	for (i = 0; i < n; i++) result[i] = MAC_ACCEPT;
    } else if (sink->isbody)
      bodyfunc(sink, src_vec, result, n);
    else
      cellfunc(sink, src_vec, result, n);
}

void
nbrMAC(SinkSPH *sink, hcell **source_vec, int *result, int n)
{
    int i;

    for (i = 0; i < n; i++) {	/* This is about 20% less efficient */
	result[i] = MAC_SPLIT_SINK;
    }
}


void
macRho(SinkSPH *sink, hcell **source_vec, int *result, int n)
{
    const double extent_sink = sink->extent;
    VxdV(const double pos_sink, = sink->pos);
    const double h = sink->h;
    Vxd(double r);
    Vxd(double dv);
    double extent_src;
    int daughters;
    int i;
    SPHbody *bp = 0;
    double projv;
    double v2;
    double rij;
    double wtij, grwtij;
    double hmean11, hmean21;
    double dxx, dwdx, dgrwdx;
    int index;
    int nbrs = 0;
    double rhoi = (double)0.0;
    double divvi = (double)0.0;
    double dr2;
    int interactions = 0;
      
    for (i = 0; i < n; i++) {
	const hcellptr source = source_vec[i];

	if (Sub_Flags(source)) {
	    const SPHcell *cp = source->ptr;
	    VxV(r, = cp->pos);	
	    extent_src = cp->bmax + cp->lap;
	    daughters = cp->daughters;
	} else {
	    bp = source->ptr;
	    VxV(r, = bp->pos);
	    extent_src = bp->h;
	    daughters = 1;
	}

	VxVx(r, -= pos_sink);	/* 8 flops */
	dr2 = Dotx(r, r); /* == 400.0000000032623 */

	if ((rij = sqrtf_fast(dr2)) > extent_src + extent_sink
	    || dr2 == (double)0.0) {
	    goto accept;
	} else if (daughters != 1) {
	    goto failed;
	}

	hmean11 = (double)2.0 / (h + bp->h);
	hmean21 = hmean11 * hmean11; /* == 0.01 */
	    
	v2 = dr2 * hmean21;
	index = v2 * invdvtable;
	if (index >= MAX_INDEX) Error("Index too large\n");
	dxx = v2 - index * dvtable;
	dwdx = (wij[index+1] - wij[index]) * invdvtable;
	wtij = (wij[index] + dwdx * dxx ) * hmean21 * hmean11;
	if (wtij < (double)0.0) Error("Negative wtij (macRho) = %g\n", wtij);
	dgrwdx = (grwij[index+1] - grwij[index]) * invdvtable;
	grwtij = (grwij[index] + dgrwdx * dxx) * hmean21 * hmean21;

	rhoi += bp->mass * wtij;
/* 	singlPrintf("M%g %g ", rhoi, bp->mass); */
	/* velocity divergence times density */
	VxVV(dv, = bp->vel, - sink->vel);
	projv = grwtij * Dotx(dv, r) * recipsqrtf(dr2);
	divvi -= bp->mass * projv;

	nbrs++;
      accept:
	interactions += daughters;
	result[i] = MAC_ACCEPT;
	continue;
      failed:
	result[i] = MAC_SPLIT_SRC;
    }

    sink->interactions += interactions;
    sink->rho += rhoi;
/*     singlPrintf("S%g %g ", sink->rho, sink->mass); */
    sink->nbrs += nbrs;
    sink->drho_dt -= divvi;
}


void
SetSPH(double visc_alpha, double visc_beta, double visc_epsilon, double heat_f1,
       double eos_gamma, int gnobj,  void bfunc(), void cfunc())
{
    Nobj = gnobj;
    alpha = visc_alpha;
    beta = visc_beta;
    epsil = visc_epsilon;
    heatf1 = heat_f1;
    Gamma = eos_gamma;
    bodyfunc = bfunc;
    cellfunc = cfunc;
      
}

void
/*SPH_setup(int dim)*/
SPH_setup(int dim, int ncoef1, double *wcoef1, int ncoef2, double *wcoef2)
{
    double v2max;
    double v, v2 /* , v3 */;
    double w, dw, dm, dm1;
    /* double dif2;  */
    /* double sum;  */  /* sum needs to be double */
    double ddvtable;
    int i, i1, j;
  

    ndim = dim;


    /* Do any diffusion-specific initialization here: */
    /* Set rmax, luminosity, etc. */


    /* maximum interaction length and step size */
    v2max = (double)4.0;
    ddvtable = v2max / NKERNEL_TABLE;
    /*i1 = (double)1.0 / dvtable;*/
    invdvtable = (double)1.0 / ddvtable;
    dvtable = ddvtable;

    /* normalization constant */
    if (ndim == 3)
      cnormk = (double)M_1_PI; /* 1./pi */
    else if (ndim == 2)
      cnormk = (double)(M_1_PI*10.0/7.0);
    else if (ndim == 1)
      cnormk = (double)(2.0/3.0);
    else
      Error("Bad ndim in sph_ktable\n");

    /* build tables */
    /* a) v less than 1 */

    /*    for (i = 0; i < i1; i++) {
	v2 = i * dvtable;
	v = sqrtf_fast(v2);
	v3 = v * v2;
	sum = 1.0 - 1.5 * v2 + 0.75 * v3;
	wij[i] = cnormk * sum;
	sum = -3.0 * v + 2.25 * v2;
	grwij[i] = cnormk * sum;
	fmass[i] = (4.0/3.0)*v3 - 1.2*v2*v3 + 0.5*v3*v3;
	fpoten[i] = (2.0/3.0)*v2 - 0.3*v2*v2 + 0.1*v2*v3 - 1.4;
        singlPrintf("%g ", grwij[i]);
	}*/
    wij[0] = cnormk*wcoef1[0];
    grwij[0] = 0.0;            /* Is this right? */
    i1 = 1.0 / ddvtable;
/*     singlPrintf("i         v          Wij           GrWij         fMass");      */
    for (i = 1; i <= i1; i++) {
       v2 = i * ddvtable;
       v = sqrt(v2);
       w = wcoef1[ncoef1-1];
       for (j = ncoef1-1; j >= 1 ; j--)
         w = w * v + wcoef1[j-1];
       dw = (ncoef1 - 1.0) * wcoef1[ncoef1-1];
       for (j = ncoef1-2; j >= 1; j--)
          dw = dw * v + j * wcoef1[j];
       wij[i] = cnormk * w;
/*        grwij[i] = cnormk * dw / v; /\*OLD VERSION WRONG!*\/ */
       grwij[i] = cnormk * dw;

       /* Enclosed mass now for m=1, h=1*/
       dm = wcoef1[ncoef1-1]/(ncoef1-1+3);
       for (j = ncoef1-2; j >= 0; j--)
          dm = dm * v + wcoef1[j]/(j+3.);
       fmass[i]=4*PI*cnormk*v*v*v*dm;
/*        singlPrintf("%d %g %g %g %g \n",  */
/* 		   i, v, wij[i], grwij[i], fmass[i]);      */

       /* Potential now: fpoten is potential for G=1, m=1, h=1*/
       fpoten[i]=(double) 0.;
       for (j = 0;j <= ncoef1-1; j++) 
	 fpoten[i]+=wcoef1[j]/(j+3.)*pow(v,(j+2.))+
	   wcoef1[j]/(j+2.)*(1-pow(v,(j+2.)))+
	   wcoef2[j]/(j+2.)*(pow(2.,(j+2.))-1.);
       fpoten[i]*=4.*PI*cnormk;

     }
   
    v=1.;
    dm1=wcoef2[ncoef2-1]/(ncoef2-1+3.);
    for (j=ncoef2-2; j>=0;j--)  
	   dm1=dm1*v + wcoef2[j]/(j+3.);      
    dm1=fmass[i1]-4*PI*dm1*cnormk*v*v*v;


     /*  b) v greater than 1 */
    /*    for (i = i1; i <= NKERNEL_TABLE; i++) {
	v2 = i * dvtable;
	v = sqrtf_fast(v2);
	v3 = v * v2;
	dif2 = (double)2.0 - v;
	sum = 0.25 * dif2 * dif2 * dif2;
	wij[i] = cnormk * sum;
	sum = -0.75 * v2 + 3.0 * v - 3.0;
	grwij[i] = cnormk * sum;
	fmass[i] = (-1.0/6.0)*v3*v3 + 1.2*v2*v3 - 3.0*v2*v2 + (8.0/3.0)*v3
	  - (1.0/15.0);
	fpoten[i] = (-1.0/30.0)*v2*v3 + 0.3*v2*v2 - v3 + (4.0/3.0)*v2 - 1.6;
          singlPrintf("%g ", grwij[i]);
  }*/

    for (i = i1+1; i <= NKERNEL_TABLE; i++) {
       v2 = i * ddvtable;
       v = sqrt(v2);
       w = wcoef2[ncoef2-1];
       for (j = ncoef2-1; j >= 1 ; j--)
         w = w * v + wcoef2[j-1];
       dw = (ncoef2 - 1.0) * wcoef2[ncoef2-1];
       for (j = ncoef2-2; j >= 1; j--)
          dw = dw * v + j * wcoef2[j];
       wij[i] = cnormk * w;
/*        grwij[i] = cnormk * dw / v; /\*OLD VERSION WRONG?!?*\/ */
       grwij[i] = cnormk * dw;

       /* Enclosed mass now */
       dm = wcoef2[ncoef2-1]/(ncoef2-1+3);
       for (j = ncoef2-2; j >= 0; j--)
          dm = dm * v + wcoef2[j]/(j+3.);
       fmass[i]=4*PI*cnormk*v*v*v*dm+dm1;
/*        singlPrintf("%d %g %g %g %g \n",  */
/* 		   i, v, wij[i], grwij[i], fmass[i]);      */

       /* Potential now: fpoten is potential for G=1, m=1, h=1*/
       fpoten[i]=(double) 0.;
       for (j=0; j <= ncoef2-1; j++)
	 fpoten[i]+=1./v*wcoef1[j]/(j+3.)+
	   1./v*wcoef2[j]/(j+3.)*(pow(v,(j+3.))-1.)+
	   wcoef2[j]/(j+2.)*(pow(2.,(j+2.))-pow(v,(j+2.)));
       fpoten[i]*=4.*PI*cnormk;
     }

    /* Make sure the mass is perfectly normalized, don't generate mass */
    for (i = 0; i <= NKERNEL_TABLE; i++) {
      fmass[i]/=fmass[NKERNEL_TABLE];
    }

    /* For interpolation of last table entry */
    wij[NKERNEL_TABLE+1] = grwij[NKERNEL_TABLE+1] = (double) 0.0;
    fpoten[NKERNEL_TABLE+1]=(double) 0.5;
    fmass[NKERNEL_TABLE+1]=(double) 1.;
}

void
SPH_oldsetup(int dim)
{
    double v2max;
    double v, v2, v3;
    double dif2;
    double sum;   /* sum needs to be double */
    int i, i1;

    ndim = dim;


    /* Do any diffusion-specific initialization here: */
    /* Set rmax, luminosity, etc. */


    /* maximum interaction length and step size */
    v2max = (double)4.0;
    dvtable = v2max / NKERNEL_TABLE;
    i1 = (double)1.0 / dvtable;
    invdvtable = i1;

    /* normalisation constant */
    if (ndim == 3)
      cnormk = (double)M_1_PI;
    else if (ndim == 2)
      cnormk = (double)(M_1_PI*10.0/7.0);
    else if (ndim == 1)
      cnormk = (double)(2.0/3.0);
    else
      Error("Bad ndim in sph_ktable\n");

    /* build tables */
    /* a) v less than 1 */

    for (i = 0; i < i1; i++) {
	v2 = i * dvtable;
	v = sqrtf_fast(v2);
	v3 = v * v2;
	sum = 1.0 - 1.5 * v2 + 0.75 * v3;
	wij[i] = cnormk * sum;
	sum = -3.0 * v + 2.25 * v2;
	grwij[i] = cnormk * sum;
	fmass[i] = (4.0/3.0)*v3 - 1.2*v2*v3 + 0.5*v3*v3;
	fpoten[i] = (2.0/3.0)*v2 - 0.3*v2*v2 + 0.1*v2*v3 - 1.4;
	if (wij[i] < 0) singlPrintf("%g ", wij[i]);
    }

    /*  b) v greater than 1 */
    for (i = i1; i <= NKERNEL_TABLE; i++) {
	v2 = i * dvtable;
	v = sqrtf_fast(v2);
	v3 = v * v2;
	dif2 = (double)2.0 - v;
	sum = 0.25 * dif2 * dif2 * dif2;
	wij[i] = cnormk * sum;
	sum = -0.75 * v2 + 3.0 * v - 3.0;
	grwij[i] = cnormk * sum;
	fmass[i] = (-1.0/6.0)*v3*v3 + 1.2*v2*v3 - 3.0*v2*v2 + (8.0/3.0)*v3
	  - (1.0/15.0);
	fpoten[i] = (-1.0/30.0)*v2*v3 + 0.3*v2*v2 - v3 + (4.0/3.0)*v2 - 1.6;
	if (wij[i] < 0) singlPrintf("%g ", wij[i]);
    }
    /* For interpolation of last table entry */
    wij[NKERNEL_TABLE+1] = grwij[NKERNEL_TABLE+1] = (double) 0.0;
}

double eos_n, eos_u;

#include "physics.h"

