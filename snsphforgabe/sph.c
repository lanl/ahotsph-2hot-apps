#include <math.h>
#include <stdlib.h>
#include "physics_sph.h"
#include "vop.h"
#include "fastflpt.h"
#include "timers.h"
#include "error.h"
#include "singlio.h"
#include "Msgs.h"
#include "mpmy.h"

#ifndef M_1_PI
#define	M_1_PI 0.31830988618379067154
#endif
#define NKERNEL_TABLE 80000
#define MAX_INDEX (NKERNEL_TABLE+2)

#define NO_UPDATE 2
#define NUMRMAXTAB 1	/* Use nearest trapped particle as defining r */

extern Timer_t sphEOS;
Counter_t SPHCnt, SPHrej, nbrMACCnt;

float ftrape = 1.0;
float ftrapb = 1.0;
float ftrapx = 1.0;

konst_s *Konst;
output_s *Outputf;
neut_out_s *Neut_out;
nu_out_s *Nu_out;
nu_lums_s *Nu_lums;
units_s *Units;
unit2_s *Unit2;
beta_s *Nubeta;
nutrap_s *Nutrap;
static nu_lums_s old_nu_lums;

static float dvtable;
static float invdvtable;
static float cnormk;
static float wij[MAX_INDEX];
static float grwij[MAX_INDEX];
static float fmass[MAX_INDEX];
static float fpoten[MAX_INDEX];
static float Gamma = (5.0/3.0);
static float alpha = 1.0;
static float Beta = 2.5;
static float epsil = 1e-2;
static float heatf1 = 1.0;
static float jtrape = 0.;
static float jtrapb = 0.;
static float jtrapx = 0.;

static float rmaxnue=0.0;
static float rmaxnueb=0.0;
static float rmaxnux=0.0;
static float rmaxnu[3];
static float rinner;
static float *grmaxnuet_tab, *grmaxnuebt_tab, *grmaxnuxt_tab;
static int grmaxtabsize;

static int ndim;
static int Nobj;
static int add_offset;
static int add_rotate;
static float offset[NDIM];
static float voffset[NDIM];
static float cosa, sina;
static void (*bodyfunc)(SinkSPH *sink, hcell **src_vec, int *res, int n);
static void (*cellfunc)(SinkSPH *sink, hcell **src_vec, int *res, int n);



#ifdef __GNUC__
#define INLINE inline
#else
#undef INLINE
#define INLINE
#endif

INLINE float MAX(float a, float b)
{
  if (a > b) 
    return a;
  else
    return b;
}

INLINE float MIN(float a, float b)
{
  if (a < b) 
    return a;
  else
    return b;
}

static int
fltcompar(const void *x, const void *y)
{
  if (*(const float *)x < *(const float *)y)
    return (1);
  else if (*(const float *)x > *(const float *)y)
    return (-1);
  else return (0);
}

void
SetSPHOffset(float *off, float *voff)
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

void
SetSPHRotate(float angle)
{
    cosa = cos(angle);
    sina = sin(angle);
    add_rotate = 1;
}

void
UnSetSPHRotate(void)
{
    cosa = sina = 0.0;
    add_rotate = 0;
}

void 
InheritSPH(const SinkSPH *from, SinkSPH *to, hcell *pp)
{
    if( to == NULL ){
	SPHbody *bp = pp->ptr;
	if (from->isbody == NO_UPDATE)
	  return;
	/* Self contribution to density */
	if (add_offset == 0 && add_rotate == 0) {
#if NDIM==3
	  bp->rho += wij[0] * from->mass / (from->h * from->h * from->h);
#else
	  bp->rho += from->xfac * wij[0] * from->mass / (from->h * from->h);
#endif
	}
	/* Must accumulate for periodic BC to work */
	/* Must initialize to zero appropriately */
	bp->dynue += from->dynue;
	bp->dynueb += from->dynueb;
	bp->dynux += from->dynux;
	bp->dunue += from->dunue;
	bp->dunueb += from->dunueb;
	bp->dunux += from->dunux;
	bp->rho += from->rho;
	bp->drho_dt += from->drho_dt;
	bp->udot += from->udot;
	bp->udot2 += from->udot2;
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
	if (bp->bghost) {
	  to->isbody = NO_UPDATE;
	  return; /* Is this right? */
	}
	VV(to->pos, = bp->pos);
	to->extent = bp->h;
	to->h = bp->h;
	to->isbody = 1;
	VV(to->vel, = bp->vel);
	to->pr = bp->pr;
	to->prnu = bp->prnu;
	to->rho_est = bp->rho_est;
	to->mass = bp->mass;
	to->vsound = bp->vsound;
	to->u = bp->u;
	to->rho = (float)0.0;
	to->drho_dt = (float)0.0;
	to->udot = (float)0.0;
	to->udot2 = (float)0.0;
	VS(to->lvel, = (float)0.0);
	to->nterms = 1;
	to->nbrs = 0;
	VS(to->M1, = (float)0.0);
	to->min_nbr_dt = 1e30;
#if NDIM==3
	to->xfac = 1.0;
#else	/* 2-d geometrical factor */
	to->xfac = 0.5*M_1_PI/
	  ((fabs(bp->pos[0]) > bp->h*0.125) ? fabs(bp->pos[0]) : bp->h*0.125);
#endif
	to->r = bp->r;
	to->dt = bp->dt;
	to->dynue = bp->dynue;
	to->dynueb = bp->dynueb;
	to->dynux = bp->dynux;
	to->dunue = bp->dunue;
	to->dunueb = bp->dunueb;
	to->dunux = bp->dunux;
	to->ynue = bp->ynue;
	to->ynueb = bp->ynueb;
	to->ynux = bp->ynux;
	to->unue = bp->unue;
	to->unueb = bp->unueb;
	to->unux = bp->unux;
	to->enuet = bp->enuet;
	to->enuebt = bp->enuebt;
	to->enuxt = bp->enuxt;
	to->dnue = bp->dnue;
	to->dnueb = bp->dnueb;
	to->dnux = bp->dnux;
	to->tempnue = bp->tempnue;
	to->tempnueb = bp->tempnueb;
	to->tempnux = bp->tempnux;
	to->etanue = bp->etanue;
	to->etanueb = bp->etanueb;
	to->etanux = bp->etanux;
	to->gshift = bp->gshift;
	to->ident = bp->ident;
    }

    if (add_offset) {
	VV(to->pos, += offset);
	VV(to->vel, += voffset);
    }
    if (add_rotate) {
        float tx, ty;
	/* When rotation approaches 90 degrees, we could spuriously get */
	/* neighbors both ways without this test */
	if (to->isbody && to->pos[1] > 0.0 && sina > 0.0) {
	  to->isbody = NO_UPDATE;
	  to->interactions = Nobj;
	  return;
	}
	if (to->isbody && to->pos[1] < 0.0 && sina < 0.0) {
	  to->isbody = NO_UPDATE;
	  to->interactions = Nobj;
	  return;
	}
	tx = to->pos[0];
	ty = to->pos[1];
        to->pos[0] = tx * cosa - ty * sina;
        to->pos[1] = tx * sina + ty * cosa;
	tx = to->vel[0];
	ty = to->vel[1];
        to->vel[0] = tx * cosa - ty * sina;
        to->vel[1] = tx * sina + ty * cosa;
    }

    if (from) {
	to->interactions = from->interactions;
    } else {
	to->interactions = 0;
    }
}

void
SPHgate(SinkSPH *sink, hcell **src_vec, int *result, int n)
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
    const float extent_sink = sink->extent;
    VxdV(const float pos_sink, = sink->pos);
    const float h = sink->h;
    Vxd(float r);
    Vxd(float dv);
    float extent_src;
    int daughters;
    int i;
    SPHbody *bp = 0;
    float projv;
    float v2;
    float rij;
    float wtij, grwtij;
    float hmean11, hmean21;
    float dxx, dwdx, dgrwdx;
    int index;
    int nbrs = 0;
    float rhoi = (float)0.0;
    float divvi = (float)0.0;
    float dr2;
    float xfac;
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
	dr2 = Dotx(r, r);

	if ((rij = sqrtf_fast(dr2)) > extent_src + extent_sink
	    || dr2 == (float)0.0) {
	    goto accept;
	} else if (daughters != 1) {
	    goto failed;
	}

	hmean11 = (float)2.0 / (h + bp->h);
	hmean21 = hmean11 * hmean11;
	    
	v2 = dr2 * hmean21;
	index = v2 * invdvtable;
	if (index >= MAX_INDEX) Error("Index too large\n");
	dxx = v2 - index * dvtable;
	dwdx = (wij[index+1] - wij[index]) * invdvtable;
	dgrwdx = (grwij[index+1] - grwij[index]) * invdvtable;
#if NDIM==3
	wtij = (wij[index] + dwdx * dxx ) * hmean21 * hmean11;
	grwtij = (grwij[index] + dgrwdx * dxx) * hmean21 * hmean21;
#else
	wtij = (wij[index] + dwdx * dxx ) * hmean21;
	/* Note: we're using the 3d expression for grwtij */
	/* grwtij = (grwij[index] + dgrwdx * dxx) * hmean21 * hmean11; */
	grwtij = (grwij[index] + dgrwdx * dxx) * hmean21 * hmean21;
#endif

	/* if (bp->mass * wtij < 0.0) Error("rhoi < 0.0\n"); */
	rhoi += bp->mass * wtij;

	/* velocity divergence times density */
	VxVV(dv, = bp->vel, - sink->vel);
#if NDIM==3
	projv = grwtij * Dotx(dv, r) * recipsqrtf(dr2);
	divvi -= bp->mass * projv;
#else
	xfac = 0.5*M_1_PI/
	  ((fabs(bp->pos[0]) > bp->h*0.125) ? fabs(bp->pos[0]) : bp->h*0.125);
	projv = grwtij *
	  ((bp->vel[0]*sink->xfac-sink->vel[0]*xfac)*r0+
	   (bp->vel[1]*sink->xfac-sink->vel[1]*xfac)*r1);
	divvi -= bp->mass * projv;
#endif

	nbrs++;
      accept:
	interactions += daughters;
	result[i] = MAC_ACCEPT;
	continue;
      failed:
	result[i] = MAC_SPLIT_SRC;
    }
    sink->interactions += interactions;
    sink->rho += rhoi * sink->xfac;
    sink->nbrs += nbrs;
    sink->drho_dt -= divvi;
    if (fabs(sink->drho_dt > 1e14)) {
	Error("bad drho %g, %d\n", sink->drho_dt, sink->ident);
    }
}


void
macSPH(SinkSPH *sink, hcell **source_vec, int *result, int n)
{
    const float extent_sink = sink->extent;
    VxdV(const float pos_sink, = sink->pos);
    VxdV(const float v, = sink->vel);
    const float h = sink->h;
    const double pro2 = (sink->pr+sink->prnu) / sink->rho_est / sink->rho_est;
    const float mass = sink->mass;
    const float rho_est = sink->rho_est;
    const float vsound = sink->vsound;
    const float u = sink->u;
    /* all other variables necessary for neutrino diffusion */
    float cthird;
    float cnue, cnueb, cnux;
    float dcnue, dcnueb, dcnux, denue, denueb, denux;
    float dgfac;
    /* end of additions */

    Vxd(float r);
    Vxd(double f);
    Vxd(float dv);
    Vxd(float smv);
    float min_nbr_dt = sink->min_nbr_dt;
    float extent_src;
    int daughters;
    int i;
    SPHbody *bp = 0;
    float hmean11, hmean21;
    int index;
    int nbrs = 0;
    double rhoi = (double)0.0;
    double divvi = (double)0.0;
    float dr2;
    Vxd(float runi);
    double dq = (double)0.0;
    float vv, vv2;
    float dxx, dwdx, dgrwdx;
    double wtij, grwtij;
    float rapm, robar1;
    double grpm, wpm;
    double poro2;
    double projv, vsbar, est_divv, t12, projv2d;
    float rij, rij1;
    float xfac;
    int interactions = 0;

    Konst = (void *)&Fortran(konst);
    Unit2 = (void *)&Fortran(unit2);
#if NDIM==2
    cthird = Konst->clight*(2./3.);
#else
    cthird = Konst->clight;
#endif
    cnue = sink->rho_est*MAX(sink->ynue,0.0);
    cnueb = sink->rho_est*MAX(sink->ynueb,0.0);
    cnux = sink->rho_est*MAX(sink->ynux,0.0);
    dcnue = dcnueb = dcnux = denue = denueb = denux = 0.0;

    VxS(f, = (float)0.0);
    VxS(smv, = (float)0.0);
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
	    if (bp->rho_est <= (float)0.0)
	      Error("RhoEst is <= 0 for %s\n", hcellPrint(source));
	}

	VxVx(r, -= pos_sink);	/* 8 flops */
	dr2 = Dotx(r, r);
	
	sink->nterms += 1;

	if ((rij = sqrtf_fast(dr2)) > extent_src + extent_sink
	    || dr2 == (float)0.0) {
	    goto accept;
	} else if (daughters != 1) {
	    goto failed;
	}

	if (bp->dt_next < min_nbr_dt) min_nbr_dt = bp->dt_next;

	hmean11 = (float)2.0 / (h + bp->h);
	hmean21 = hmean11 * hmean11;

	vv2 = dr2 * hmean21;	/* v2 and v renamed to avoid conflict */
	index = vv2 * invdvtable;
	if (index >= MAX_INDEX) Error("Index too large\n");
	vv = rij * hmean11;
	dxx = vv2 - index * dvtable;
	dwdx = (wij[index+1] - wij[index]) * invdvtable;
	dgrwdx = (grwij[index+1] - grwij[index]) * invdvtable;
#if NDIM==3
	wtij = (wij[index] + dwdx * dxx) * hmean21 * hmean11;
	grwtij = (grwij[index] + dgrwdx * dxx) * hmean21 * hmean21;
	xfac = 1.0;
#else
	wtij = (wij[index] + dwdx * dxx) * hmean21;
	/* Note: we're using the 3d expression for grwtij */
	/* grwtij = (grwij[index] + dgrwdx * dxx) * hmean21 * hmean11; */
	grwtij = (grwij[index] + dgrwdx * dxx) * hmean21 * hmean21;
	xfac = 0.5*M_1_PI/
	  ((fabs(bp->pos[0]) > bp->h*0.125) ? fabs(bp->pos[0]) : bp->h*0.125);
#endif

	rapm = mass / bp->mass;
	robar1 = (float)2.0 / (rho_est + bp->rho_est);	
	grpm = bp->mass * grwtij;
	wpm = bp->mass * wtij;

	poro2 = grpm * (pro2 + (bp->pr + bp->prnu)/bp->rho_est/bp->rho_est);
	rij1 = (float)1.0 / rij;
	VxVVx(dv, = bp->vel, - v);
	VxVx(smv, += robar1 * wpm * dv);
	rhoi += wpm;
#if NDIM==3
	VxVx(runi, = rij1 * r);
	projv = Dotx(dv, runi);
#else
	VxVx(runi, = xfac * r);
	projv = ((bp->vel[0]*sink->xfac/xfac-v0)*runi0+
		 (bp->vel[1]*sink->xfac/xfac-v1)*runi1);
#endif
	VxVx(f, += poro2 * runi);
	divvi -= grpm * projv;

	/* artificial viscosity and energy dissipation */
	vsbar = (float)0.5 * (vsound + bp->vsound);
#if NDIM==3
	if (projv < (float)0.0 && alpha != (float)0.0) {
	    est_divv = projv * vv / (vv2 + epsil);
	    t12 = grpm * est_divv * (Beta * est_divv - alpha * vsbar) * robar1;
	    VxVx(f, += t12 * runi);
	    dq += t12 * projv;
	}
#else
	projv2d = Dotx(dv, r);
	if (projv2d < (float)0.0 && alpha != (float)0.0) {
	    est_divv = projv2d * hmean11 / (vv2 + epsil);
	    t12 = grpm * est_divv * (Beta * est_divv - alpha * vsbar) * robar1;
	    VxVx(f, += t12 * runi);
	    dq += t12 * projv;	/* why not projv2d? */
	}
#endif
	/* artificial heat conduction */
	if (heatf1 != (float)0.0) {
	    float qmean = heatf1 * vsbar / hmean11;
	    t12 = grpm * qmean * (bp->u - u) * robar1;
	    dq -= t12;
	}

	/* mike...  it seems that neutrino diffusion should be here
	   Neutrino diffusion for the "trapped" particles */
	/* rinner is the inner boundary */
#if NDIM==3
	if (grpm < 0.) {
#else
	if (grpm < 0. && sink->r > rinner) {
#endif
	  float gij, gji, wnuij, dnueij, tnueij, cnuej;
	  float enuej, dcin, dcout;
	  float dnuebij, tnuebij, cnuebj;
	  float enuebj;
	  float tnuxij, dnuxij, cnuxj;
	  float enuxj;
	  /* additions to do neutrino blocking */
	  float blockij,blockji,expij,expji;
#if NDIM==3
	  dgfac=-grpm/bp->rho_est*rij1;
	  wnuij=cthird*dgfac*rij;
#else
	  dgfac=-grpm/bp->rho_est*(sink->xfac+xfac);
	  wnuij=cthird*dgfac*rij;
#endif
	  /* Flux limiter removed 7/3/2001 */
	  /*  gravitational redshift */
	  gij=sink->gshift/bp->gshift;
	  gji=1./gij;
	  /* wave coefficient cthird=speed of light/3 */
	  if (bp->r < rmaxnue) {
	    /* average diffusion coefficients */
	    dnueij=dgfac*2.0*cthird*sink->dnue*bp->dnue
		/(sink->dnue+bp->dnue);
	    /* this is the flux limiter: */
	    /* I have just written max,min in fortran */
	    tnueij=MIN(dnueij,wnuij);
	    cnuej=bp->rho_est*MAX(bp->ynue,0.);
	    if (sink->r < rmaxnue ) {
	      /* effective energies: */
	      enuej=gji*bp->enuet;
	      /* changes to add blocking */
	      /* note that we are now using tempnu and etanu */
	      expij=exp(MIN(gij*sink->enuet/bp->tempnue-bp->etanue,50.));
	      expji=exp(MIN(enuej/sink->tempnue-sink->etanue,50.));
	      blockij=expij/(1+expij);
	      blockji=expji/(1+expji);
	      dcin=tnueij*cnuej*blockji;
	      dcout=tnueij*cnue*blockij;
	      /* end of changes for electron neutrinos */
	      dcnue += (dcin-dcout);
	      /*  neutrino e-density changes 
		  (units aren't correct - get fixed at the end) */
	      denue+=(dcin*enuej-dcout*sink->enuet);
	    } else if (cnue < cnuej) {
	      dcin=tnueij*cnuej;
	      dcnue += dcin;
	      enuej=gji*bp->enuet;
	      denue+=dcin*enuej;
	    } 
	  } else if (sink->r < rmaxnue) {
	    cnuej=bp->rho_est*MAX(bp->ynue,0.);
            if (cnuej < cnue) {
	      dnueij=dgfac*2.0*cthird*sink->dnue*bp->dnue
		/(sink->dnue+bp->dnue);
	      tnueij=MIN(dnueij,wnuij);
	      dcin=tnueij*cnue;
	      dcnue -= dcin;
	      denue -= dcin*sink->enuet;
            }
	  }
	  if (bp->r < rmaxnueb) {
	    /* average diffusion coefficients */
	    dnuebij=dgfac*2.0*cthird*sink->dnueb*bp->dnueb
		/(sink->dnueb+bp->dnueb);
	    /* this is the flux limiter: */
	    /* I have just written max,min in fortran */
	    tnuebij=MIN(dnuebij,wnuij);
	    cnuebj=bp->rho_est*MAX(bp->ynueb,0.);
	    if (sink->r < rmaxnueb) {
	      /* effective energies: */
	      enuebj=gji*bp->enuebt;
	      /* changes to add blocking */
	      expij=exp(MIN(gij*sink->enuebt/bp->tempnueb-bp->etanueb,50.));
	      expji=exp(MIN(enuebj/sink->tempnueb-sink->etanueb,50.));
	      blockij=expij/(1+expij);
	      blockji=expji/(1+expji);
	      dcin=tnuebij*cnuebj*blockji;
	      dcout=tnuebij*cnueb*blockij;
	      /* end of changes for anti-electron neutrinos */
	      dcnueb += (dcin-dcout);
	      /*  neutrino e-density changes 
		  (units aren't correct - get fixed at the end) */
	      denueb += (dcin*enuebj-dcout*sink->enuebt);
	    } else if (cnueb < cnuebj) {
	      dcin=tnuebij*cnuebj;
	      dcnueb+=dcin;
	      enuebj=gji*bp->enuebt;
	      denueb += dcin*enuebj;
	    }
	  } else if (sink->r < rmaxnueb) {
	    cnuebj=bp->rho_est*MAX(bp->ynueb,0.);
            if (cnuebj < cnueb) {
	      dnuebij=dgfac*2.0*cthird*sink->dnueb*bp->dnueb
		/(sink->dnueb+bp->dnueb);
	      tnuebij=MIN(dnuebij,wnuij);
	      dcin=tnuebij*cnueb;
	      dcnueb -= dcin;
	      denueb -= dcin*sink->enuebt;
            }
	  }
	  if (bp->r < rmaxnux) {
	    /* average diffusion coefficients */
	    dnuxij=dgfac*2.0*cthird*sink->dnux*bp->dnux
		/(sink->dnux+bp->dnux);
	    /* this is the flux limiter: */
	    /* I have just written max,min in fortran */
	    tnuxij=MIN(dnuxij,wnuij);
	    cnuxj=bp->rho_est*MAX(bp->ynux,0.);
	    if (sink->r < rmaxnux) {
	      /* effective energies: */
	      enuxj=gji*bp->enuxt;
	      /* changes to add blocking */
	      expij=exp(MIN(gij*sink->enuxt/bp->tempnux-bp->etanux,50.));
	      expji=exp(MIN(enuxj/sink->tempnux-sink->etanux,50.));
	      blockij=expij/(1+expij);
	      blockji=expji/(1+expji);
	      dcin=tnuxij*cnuxj*blockji;
	      dcout=tnuxij*cnux*blockij;
	      /* end of changes for mu/tau neutrinos */
	      dcnux += (dcin-dcout);
	      /*  neutrino e-density changes 
		  (units aren't correct - get fixed at the end) */
	      denux += (dcin*enuxj-dcout*sink->enuxt);
	    } else if (cnux < cnuxj) {
	      /* undo symmetry, if neighbor trapped but sink is not */
	      dcin=tnuxij*cnuxj;
	      dcnux += dcin;
	      enuxj=gji*bp->enuxt;
	      denux += dcin*enuxj;
	    }
	  } else if (sink->r < rmaxnux) {
	    cnuxj=bp->rho_est*MAX(bp->ynux,0.);
            if (cnuxj < cnux) {
	      dnuxij=dgfac*2.0*cthird*sink->dnux*bp->dnux
		/(sink->dnux+bp->dnux);
	      tnuxij=MIN(dnuxij,wnuij);
	      dcin=tnuxij*cnux;
	      dcnux -= dcin;
	      denux -= dcin*sink->enuxt;
            }
	  }
	}
	nbrs++;
      accept:
	IncrCounter(&SPHrej);
	interactions += daughters;
	result[i] = MAC_ACCEPT;
	continue;
      failed:
	IncrCounter(&SPHrej);
	result[i] = MAC_SPLIT_SRC;
    }
    sink->dynue += dcnue/sink->rho_est;
    sink->dynueb += dcnueb/sink->rho_est;
    sink->dynux += dcnux/sink->rho_est;
    sink->dunue += denue*Unit2->umevnuc/sink->rho_est;
    sink->dunueb += denueb*Unit2->umevnuc/sink->rho_est;
    sink->dunux += denux*Unit2->umevnuc/sink->rho_est;
    sink->interactions += interactions;
    sink->rho += rhoi * sink->xfac;
    sink->drho_dt -= divvi;
    if (fabs(sink->drho_dt > 1e14)) {
	Error("bad drho %g, %d\n", sink->drho_dt, sink->ident);
    }
    sink->udot += (float)0.5 * dq;
    sink->udot2 += (float)0.5 * dq;
    VVx(sink->M1, += f);
    VVx(sink->lvel, += smv);
    sink->nbrs += nbrs;
    sink->nterms += nbrs*8;
    sink->min_nbr_dt = min_nbr_dt;
}

void
SetSPH(float visc_alpha, float visc_beta, float visc_epsilon, 
       float heat_f1,float eos_gamma, int gnobj,  void bfunc(), void cfunc())
{
    Nobj = gnobj;
    alpha = visc_alpha;
    Beta = visc_beta;
    epsil = visc_epsilon;
    heatf1 = heat_f1;
    Gamma = eos_gamma;
    bodyfunc = bfunc;
    cellfunc = cfunc;
}

void
SPHaux(float rb)
{
  float tmp;
  rinner = rb;
  singlPrintf("rmaxnuet: %12g rmaxnuebt: %12g rmaxnuxt: %12g\n",
	      rmaxnu[0], rmaxnu[1], rmaxnu[2]);

  /* Limit rate at which rmaxnus can change */
  if (rmaxnue == 0.0) rmaxnue = rmaxnu[0]*0.9;
  else {
    tmp = rmaxnu[0]/rmaxnue;
    if (tmp > 1.1) tmp = 1.1;
    else if (tmp < 0.95) tmp = 0.95;
    rmaxnue *= tmp;
  }
  if (rmaxnueb == 0.0) rmaxnueb = rmaxnu[1]*0.9;
  else {
    tmp = rmaxnu[1]/rmaxnueb;
    if (tmp > 1.1) tmp = 1.1;
    else if (tmp < 0.95) tmp = 0.95;
    rmaxnueb *= tmp;
  }
  if (rmaxnux == 0.0) rmaxnux = rmaxnu[2]*0.9;
  else {
    tmp = rmaxnu[2]/rmaxnux;
    if (tmp > 1.1) tmp = 1.1;
    else if (tmp < 0.95) tmp = 0.95;
    rmaxnux *= tmp;
  }
  singlPrintf("rmaxnue:  %12g rmaxnueb:  %12g rmaxnux:  %12g\n",
	      rmaxnue, rmaxnueb, rmaxnux);
}

void
Getrmax(float *maxnue, float *maxnueb, float *maxnux, 
	float *enue, float *enueb, float *enux, 
	float *e2nue, float *e2nueb, float *e2nux)
{
  *maxnue = rmaxnue;
  *maxnueb = rmaxnueb;
  *maxnux = rmaxnux;
  *enue = old_nu_lums.enue;
  *enueb = old_nu_lums.enueb;
  *enux = old_nu_lums.enux;
  *e2nue = old_nu_lums.e2nue;
  *e2nueb = old_nu_lums.e2nueb;
  *e2nux = old_nu_lums.e2nux;
}

void
SPH_setup(int dim, int ncoef1, double *wcoef1, int ncoef2, double *wcoef2,
	  float maxnue, float maxnueb, float maxnux,
	  float enue, float enueb, float enux,
	  float e2nue, float e2nueb, float e2nux, 
	  float ftrapex, float ftrapbx, float ftrapxx)
{
    double v2max;
    double v, v2;
    double w, dw;
    double ddvtable;
    double dm, dm1;
    int i, i1, j;

    rmaxnue = rmaxnu[0] = maxnue;
    rmaxnueb = rmaxnu[1] = maxnueb;
    rmaxnux = rmaxnu[2] = maxnux;
    old_nu_lums.enue = enue;
    old_nu_lums.enueb = enueb;
    old_nu_lums.enux = enux;
    old_nu_lums.e2nue = e2nue;
    old_nu_lums.e2nueb = e2nueb;
    old_nu_lums.e2nux = e2nux;
    ftrape = ftrapex;
    ftrapb = ftrapbx;
    ftrapx = ftrapxx;
    ndim = dim;

    singlPrintf("Using %dth particle to set rmaxnus\n", NUMRMAXTAB);
    grmaxtabsize = MPMY_Nproc()*NUMRMAXTAB;
    if (MPMY_Procnum() == 0) {
      grmaxnuet_tab = Malloc(grmaxtabsize*sizeof(float));
      grmaxnuebt_tab = Malloc(grmaxtabsize*sizeof(float));
      grmaxnuxt_tab = Malloc(grmaxtabsize*sizeof(float));
    }

    /* maximum interaction length and step size */
    v2max = 4.0;
    ddvtable = v2max / NKERNEL_TABLE;
    invdvtable = 1.0 / ddvtable;
    dvtable = ddvtable;

    /* normalization constant */
    if (ndim == 3) 
      cnormk = M_1_PI;
    else if (ndim == 2)
      cnormk = M_1_PI*10.0/7.0;
    else if (ndim == 1)
      cnormk = 2.0/3.0;
    else
      Error("Bad ndim in sph_ktable\n");

    /* build tables */
    /* a) v less than 1 */

    wij[0] = cnormk*wcoef1[0];
    grwij[0] = 0.0;
    fmass[0] = 0.0;
    fpoten[0] = 0.0;
    for (j = 0; j <= ncoef1-1; j++) 
	fpoten[0] += wcoef1[j]/(j+2.0) +
	    wcoef2[j]/(j+2.) * (pow(2.0, j+2.0)-1.0);
    fpoten[0]*=4.*M_PI*cnormk;

    i1 = 1.0 / ddvtable;
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
	grwij[i] = cnormk * dw;

       /* Enclosed mass now for m=1, h=1*/
       dm = wcoef1[ncoef1-1]/(ncoef1-1+3);
       for (j = ncoef1-2; j >= 0; j--)
	   dm = dm * v + wcoef1[j]/(j+3.);
       fmass[i]=4*M_PI*cnormk*v*v*v*dm;

       /* Potential now: fpoten is potential for G=1, m=1, h=1*/
       fpoten[i]=(double)0.;
       for (j = 0;j <= ncoef1-1; j++) 
	   fpoten[i]+=wcoef1[j]/(j+3.)*pow(v,(j+2.))+
	       wcoef1[j]/(j+2.)*(1-pow(v,(j+2.)))+
	       wcoef2[j]/(j+2.)*(pow(2.,(j+2.))-1.);
       fpoten[i]*=4.*M_PI*cnormk;
    }

    v=1.;
    dm1=wcoef2[ncoef2-1]/(ncoef2-1+3.);
    for (j=ncoef2-2; j>=0;j--)  
	dm1=dm1*v + wcoef2[j]/(j+3.);      
    dm1=fmass[i1]-4*M_PI*dm1*cnormk*v*v*v;

    /*  b) v greater than 1 */
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
	grwij[i] = cnormk * dw;

	/* Enclosed mass now */
	dm = wcoef2[ncoef2-1]/(ncoef2-1+3);
	for (j = ncoef2-2; j >= 0; j--)
	    dm = dm * v + wcoef2[j]/(j+3.);
	fmass[i]=4*M_PI*cnormk*v*v*v*dm+dm1;

	/* Potential now: fpoten is potential for G=1, m=1, h=1*/
	fpoten[i]=(double) 0.;
	for (j=0; j <= ncoef2-1; j++)
	    fpoten[i]+=1./v*wcoef1[j]/(j+3.)+
		1./v*wcoef2[j]/(j+3.)*(pow(v,(j+3.))-1.)+
		wcoef2[j]/(j+2.)*(pow(2.,(j+2.))-pow(v,(j+2.)));
	fpoten[i]*=4.*M_PI*cnormk;
    }

    /* Make sure the mass is perfectly normalized, don't generate mass */
    for (i = 0; i <= NKERNEL_TABLE; i++) {
	fmass[i]/=fmass[NKERNEL_TABLE];
    }

    /* For interpolation of last table entry */
    wij[NKERNEL_TABLE+1] = grwij[NKERNEL_TABLE+1] = 0.0;
    fpoten[NKERNEL_TABLE+1]=(double)0.5;
    fmass[NKERNEL_TABLE+1]=(double)1.0;

#ifdef WRITE_WIJ
    for (i = 0; i <= NKERNEL_TABLE; i++) 
      singlPrintf("%5d %8g %8g\n", i, wij[i], grwij[i]);
#endif
    Nu_lums = (void *)&Fortran(nulums);
    Nu_lums->enue = 3.0;
    Nu_lums->enueb = 3.0;
    Nu_lums->enux = 3.0;
    Nu_lums->e2nue = 3.0;
    Nu_lums->e2nueb = 3.0;
    Nu_lums->e2nux = 3.0;
}


void
update_final(SPHbody *btab, int nobj, float dt, int *limit_high, int *limit_low, float dttol)
{
    SPHbody *p;
    float avokb, sfac;
    int nwarn = 0;
    
    Unit2 = (void *)&Fortran(unit2);
    Units = (void *)&Fortran(units);
    avokb = 6.02e23*1.381e-16;
    sfac = avokb*Unit2->utemp/Units->uergg;      

    for (p = btab; p < btab+nobj; p++) {
	if (!SPH_need_update(p)) continue;
	if (p->bghost) continue; /* Is this right? */
	VV(p->acc, += p->grav_acc);
#if NDIM==2
	/* p->hdot = (float)(-1.0/2.0) * p->h * p->drho_dt / p->rho; */
	/* for SN code, rho and drho/dt are scaled to 3d quantities */
	/* so we use 1/3 instead of 1/2 */
	p->hdot = (float)(-1.0/3.0) * p->h * p->drho_dt / p->rho;
#else
	p->hdot = (float)(-1.0/3.0) * p->h * p->drho_dt / p->rho;
#endif
	if (p->hdot * dt > dttol*p->h) {
	    SeriousWarning("Hdot limit (high)\n%s\n", PrintSPHBodyContents(p));
	    p->hdot = dttol*p->h/dt;
	    ++*limit_high;
	}
	if (p->hdot * dt < -dttol*p->h) {
	    SeriousWarning("Hdot limit (low)\n%s\n", PrintSPHBodyContents(p));
	    p->hdot = -dttol*p->h/dt;
	    ++*limit_low;
	}
	{
	  float tmp;
	  tmp = p->udot;
	  p->dq = p->udot;

	  p->udot -= p->dunu;
	  /* xxx check these changes to udot2 */
	  p->udot2 -= p->dunu;
	  if (p->ebeta < 1.0) 
	    p->udot2 += sfac*p->dye*(p->xmuhat-p->xmue);
#if NDIM==3
	  p->udot += (p->drho_dt / p->rho) * ((p->pr + p->prnu) / p->rho);
#else
				/* WTF is vx/x here? */
	  p->udot += p->pr/p->rho * (p->drho_dt/p->rho - p->vel[0]/p->pos[0]);
	  /* p->udot += p->drho_dt * p->pr / (p->rho * p->rho);  */

#endif
	}
	if (!finite(p->udot)) {
	    SeriousWarning("Bad value for udot\n%s\n", PrintSPHBodyContents(p));
	    p->udot = 0.0;
	}
	p->udot2 /= p->temp;
	if (!finite(p->udot2)) {
	    Error("Bad value for udot2\n%s\n", PrintSPHBodyContents(p));
	}

	if (p->temp/p->temprev > (1.0+dttol)) {
	    SeriousWarning("Tempdot limit (high)\n%s\n", 
		PrintSPHBodyContents(p));
/* 	    p->temp = p->temprev * 1.0+dttol; */
	    *limit_high += 20;
	}
	if (p->temp/p->temprev < (1.0-dttol)) {
	    SeriousWarning("Tempdot limit (low)\n%s\n", 
		PrintSPHBodyContents(p));
/* 	    p->temp = p->temprev * 1.0-dttol; */
	    ++*limit_low;
	}

	if (p->u > 20.0 && p->udot * dt > dttol * p->u) {
	    Msg_do("Udot limit (high)\n%s\n", PrintSPHBodyContents(p));
	    p->udot = dttol*p->u/dt;
	    ++*limit_high;
	} else if (p->u < -20.0 && p->udot * dt < dttol * p->u) {
	    Msg_do("Udot limit (high)\n%s\n", PrintSPHBodyContents(p));
	    p->udot = dttol*p->u/dt;
	    ++*limit_high;
	}

	if (p->ye > 1e-4 && p->dye * dt > dttol * p->ye) {
              if (++nwarn < 3) Msg_do("ye limit (high)\n%s\n",
				      PrintSPHBodyContents(p));
              p->dye = dttol*p->ye/dt;
              ++*limit_high;
	}
	if (p->dye * dt < -0.9 * p->ye) {
              if (++nwarn < 3) Msg_do("ye limit (low)\n%s\n",
				      PrintSPHBodyContents(p));
              p->dye = -0.9*p->ye/dt;
              ++*limit_low;
	}

	if (p->ynue > 1e-4 && p->dynue * dt > dttol * p->ynue) {
              if (++nwarn < 3) Msg_do("ynue limit (high) %.2f\n%s\n",
				      p->dynue*dt/p->ynue,
				      PrintSPHBodyContents(p));
              p->dynue = dttol*p->ynue/dt;
	      p->dunue = dttol*p->unue/dt;
              ++*limit_high;
	}
	if (p->dynue * dt < -0.9 * p->ynue) {
              if (++nwarn < 3) Msg_do("ynue limit (low) %.2f\n%s\n",
				      p->dynue*dt/p->ynue,
				      PrintSPHBodyContents(p));
              p->dynue = -0.9*p->ynue/dt;
	      p->dunue = -0.9*p->unue/dt;
              ++*limit_low;
	}
	if (p->unue > 20.0 && p->dunue * dt > dttol * p->unue) {
	    if (++nwarn < 3) Msg_do("unue limit (high) %.2f\n%s\n", 
				    p->dunue*dt/p->unue,
				    PrintSPHBodyContents(p));
            p->dynue = dttol*p->ynue/dt;
	    p->dunue = dttol*p->unue/dt;
	    ++*limit_high;
	}
	if (p->dunue * dt < -0.9 * p->unue) {
	    if (++nwarn < 3) Msg_do("unue limit (low) %.2f\n%s\n", 
				    p->dunue*dt/p->unue,
				    PrintSPHBodyContents(p));
            p->dynue = -0.9*p->ynue/dt;
	    p->dunue = -0.9*p->unue/dt;
	    ++*limit_low;
	}

	if (p->ynueb > 1e-4 && p->dynueb * dt > dttol * p->ynueb) {
              if (++nwarn < 3) Msg_do("ynueb limit (high)\n%s\n",
				      PrintSPHBodyContents(p));
              p->dynueb = dttol*p->ynueb/dt;
	      p->dunueb = dttol*p->unueb/dt;
              ++*limit_high;
	}
	if (p->dynueb * dt < -0.9 * p->ynueb) {
              if (++nwarn < 3) Msg_do("ynueb limit (low)\n%s\n",
				      PrintSPHBodyContents(p));
              p->dynueb = -0.9*p->ynueb/dt;
	      p->dunueb = -0.9*p->unueb/dt;
              ++*limit_low;
	}
	if (p->unueb > 20.0 && p->dunueb * dt > dttol * p->unueb) {
	    if (++nwarn < 3) Msg_do("unueb limit (high)\n%s\n", 
				    PrintSPHBodyContents(p));
            p->dynueb = dttol*p->ynueb/dt;
	    p->dunueb = dttol*p->unueb/dt;
	    ++*limit_high;
	}
	if (p->dunueb * dt < -0.9 * p->unueb) {
	    if (++nwarn < 3) Msg_do("unueb limit (low)\n%s\n", 
				    PrintSPHBodyContents(p));
            p->dynueb = -0.9*p->ynueb/dt;
	    p->dunueb = -0.9*p->unueb/dt;
	    ++*limit_low;
	}

	if (p->ynux > 1e-4 && p->dynux * dt > dttol * p->ynux) {
              if (++nwarn < 3) Msg_do("ynux limit (high)\n%s\n",
				      PrintSPHBodyContents(p));
              p->dynux = dttol*p->ynux/dt;
	      p->dunux = dttol*p->unux/dt;
              ++*limit_high;
	}
	if (p->dynux * dt < -0.9 * p->ynux) {
              if (++nwarn < 3) Msg_do("ynux limit (low)\n%s\n",
				      PrintSPHBodyContents(p));
              p->dynux = -0.9*p->ynux/dt;
	      p->dunux = -0.9*p->unux/dt;
              ++*limit_low;
	}
	if (p->unux > 20.0 && p->dunux * dt > dttol * p->unux) {
	    if (++nwarn < 3) Msg_do("unux limit (high)\n%s\n", 
				    PrintSPHBodyContents(p));
            p->dynux = dttol*p->ynux/dt;
	    p->dunux = dttol*p->unux/dt;
	    ++*limit_high;
	}
	if (p->dunux * dt < -0.9 * p->unux) {
	    if (++nwarn < 3) Msg_do("unux limit (low)\n%s\n", 
				    PrintSPHBodyContents(p));
            p->dynux = -0.9*p->ynux/dt;
	    p->dunux = -0.9*p->unux/dt;
	    ++*limit_low;
	}
    }
}

void
update_intermediate(SPHbody *btab, int nobj, float dt_last, int flag, int *limit, float xmtheo)
{
    SPHbody *p;
    double rho, u, u2, ye, temp, abar, xp, xn;
    double xpf, p2, p3, p4, temprev, rhoprev, xpprev, xnprev, yeprev, ufreez;
    int ifleos;
    int iident;
    int i;
    float gshift, steps, h, mass;
    float enuef, enuebf, enuxf, e2nuef, e2nuebf, e2nuxf;
    float enues, enuebs, enuxs, dee, deeb, dex;
    float rlumnuef,rlumnuebf,rlumnuxf;
    float dlumnu;
    float xlumnu,ylumnu,zlumnu;
    double hgw1,hgw2,hgw3,hgw4,hgw5,hgw6;
    double ixx, iyy, izz, ixy, ixz, iyz;
    float dtrapnue, dtrapnueb;
    MPMY_Comm_request req;
    float rmaxnuet_tab[NUMRMAXTAB];
    float rmaxnuebt_tab[NUMRMAXTAB];
    float rmaxnuxt_tab[NUMRMAXTAB];
    int procnum = MPMY_Procnum();

    StartTimer(&sphEOS);
    Outputf = (void *)&Fortran(output);
    Neut_out = (void *)&Fortran(neutout);
    Nu_out = (void *)&Fortran(nuout);
    Nu_lums = (void *)&Fortran(nulums);
    Units = (void *)&Fortran(units);
    Unit2 = (void *)&Fortran(unit2);
    Nubeta = (void *)&Fortran(beta);
    Nutrap = (void *)&Fortran(nutrap);
    rlumnuef = rlumnuebf = rlumnuxf = 0.0;
    dlumnu = 0.0;
    xlumnu = ylumnu = zlumnu = 0.0;
    hgw1 = hgw2 = hgw3 = hgw4 = hgw5 = hgw6 = 0.0;
    ixx = iyy = izz = ixy = ixz = iyz = 0.0;
    enuef = enuebf = enuxf = e2nuef = e2nuebf = e2nuxf = 0.0;
    enues = enuebs = enuxs = dee = deeb = dex = 0.0;
    xmtheo = 1;
    steps = dt_last;
    rmaxnu[0] = rmaxnu[1] = rmaxnu[2] = 0.0;
    dtrapnue = dtrapnueb = 0.0;
    for (i = 0; i < NUMRMAXTAB; i++) {
      rmaxnuet_tab[i] = rmaxnuebt_tab[i] = rmaxnuxt_tab[i] = 0.0;
    }

    for (p = btab; p < btab+nobj; p++) {
	if (!SPH_need_update(p)) continue;
	if (p->bghost) continue; /* Is this right? */
	if (flag) p->rho_est = p->rho + p->drho_dt * dt_last;
	else p->rho_est = p->rho;
	if (p->rho_est <= (float)0.0) 
	  Error("Rho_est is 0\n%s\n", PrintSPHBodyContents(p));
	rho = p->rho_est;
	u = p->u;
	u2 = p->u2;
	ye = p->ye;
	temp = p->temp;
	ifleos = p->ifleos;
	iident = p->ident;
	abar = p->abar;
	xp = p->xp;
	xn = p->xn;
	xpf = p->xpf;
	p2 = p->p2;
	p3 = p->p3;
	p4 = p->p4;
	temprev = p->temprev;
	p->temprev = temp;	/* temp gets updated in eos3 */
	rhoprev = p->rhoprev;
	xpprev = p->xpprev;
	p->xpprev = xp;
	xnprev = p->xnprev;
	p->xnprev = xn;
	yeprev = p->yeprev;
	ufreez = p->ufreez;

#if 0
	if (p->ident == 78627) {
	  Msg_do("id %d rho %g u %g u2 %g ye %g temp %g ifleos %d abar %g xp %g xn %g\n", p->ident, rho, u, u2, ye, temp, ifleos, abar, xp, xn);
	}
#endif

	Fortran(eos3)(&rho, &u, &u2, &ye, &temp, &ifleos, &abar, &xp, &xn, 
		      &xpf, &p2, &p3, &p4, &temprev, &rhoprev, 
		      &xpprev, &xnprev, &yeprev, &ufreez, &iident, &procnum);
	p->u = u;
	p->ifleos = ifleos;
	p->u2 = u2;
	p->temp = temp;		/* updated in eos3 */
	p->abar = abar;
	p->xp = xp;
	p->xn = xn;
	p->vsound = Outputf->vsound;
	p->pr = Outputf->pr;
	p->eta = Outputf->eta;
	p->xmuhat = Outputf->xmuhat;
	p->xmue = Outputf->xmue;
	p->xpf = xpf;
	p->p2 = p2;
	p->p3 = p3;
	p->p4 = p4;
	p->rhoprev = rhoprev;
	p->yeprev = yeprev;
#if 0				/* don't update! */
	p->temprev = temprev;
	p->xpprev = xpprev;
	p->xnprev = xnprev;
#endif
	p->ufreez = ufreez;
	h=p->h;
	mass=p->mass;
	gshift = p->gshift;
	Nu_lums->enue = old_nu_lums.enue;
	Nu_lums->enueb = old_nu_lums.enueb;
	Nu_lums->enux = old_nu_lums.enux;
	Nu_lums->e2nue = old_nu_lums.e2nue;
	Nu_lums->e2nueb = old_nu_lums.e2nueb;
	Nu_lums->e2nux = old_nu_lums.e2nux;
	Neut_out->dye = p->dye;
	Neut_out->dynue = p->dynue;
	Neut_out->dynueb = p->dynueb;
	Neut_out->dynux = p->dynux;
	Neut_out->tempnue = p->tempnue;
	Neut_out->tempnueb = p->tempnueb;
	Neut_out->tempnux = p->tempnux;
	Neut_out->enuet = p->enuet;
	Neut_out->enuebt = p->enuebt;
	Neut_out->enuxt = p->enuxt;
	Neut_out->dnue = p->dnue;
	Neut_out->dnueb = p->dnueb;
	Neut_out->dnux = p->dnux;
	Neut_out->dunue = p->dunue;
	Neut_out->dunueb = p->dunueb;
	Neut_out->dunux = p->dunux;
	Neut_out->dunu = p->dunu;
	Neut_out->etanue = p->etanue;
	Neut_out->etanueb = p->etanueb;
	Neut_out->etanux = p->etanux;
	Neut_out->prnu = p->prnu;
	Fortran(neutrino)(&steps, &p->rho_est, &p->ye, &p->xp, &p->xn, &p->h, 
			  &Outputf->xheavy, &Outputf->xalpha, &Outputf->yeh, 
			  &p->eta, &p->temp, &p->abar, &p->gshift, 
			  &p->r, &p->mass, &p->vsound, 
			  &p->xmuhat, &p->ynue, &p->ynueb, &p->ynux,
			  &p->unue, &p->unueb, &p->unux, 
			  &rmaxnue,& rmaxnueb, &rmaxnux,
			  &ftrape, &ftrapb, &ftrapx,
                          &jtrape, &jtrapb, &jtrapx, &p->ident);
	p->dye = Neut_out->dye;
	p->dynue = Neut_out->dynue;
	p->dynueb = Neut_out->dynueb;
	p->dynux = Neut_out->dynux;
	p->tempnue = Neut_out->tempnue;
	p->tempnueb = Neut_out->tempnueb;
	p->tempnux = Neut_out->tempnux;
	p->enuet = Neut_out->enuet;
	p->enuebt = Neut_out->enuebt;
	p->enuxt = Neut_out->enuxt;
	p->dnue = Neut_out->dnue;
	p->dnueb = Neut_out->dnueb;
	p->dnux = Neut_out->dnux;
	p->dunue = Neut_out->dunue;
	p->dunueb = Neut_out->dunueb;
	p->dunux = Neut_out->dunux;
	p->dunu = Neut_out->dunu;
	p->etanue = Neut_out->etanue;
	p->etanueb = Neut_out->etanueb;
	p->etanux = Neut_out->etanux;
	p->prnu = Neut_out->prnu;
	p->ebeta = Nubeta->ebetaeq;
	p->pbeta = Nubeta->pbetaeq;
	for (i = 0; i < NUMRMAXTAB; i++) {
	    if (Nu_out->rmxnue > rmaxnuet_tab[i]) {
	      rmaxnuet_tab[i] = Nu_out->rmxnue;
	      break;
	    }
	}
	for (i = 0; i < NUMRMAXTAB; i++) {
	    if (Nu_out->rmxnueb > rmaxnuebt_tab[i]) {
	      rmaxnuebt_tab[i] = Nu_out->rmxnueb;
	      break;
	    }
	}
	for (i = 0; i < NUMRMAXTAB; i++) {
	    if (Nu_out->rmxnux > rmaxnuxt_tab[i]) {
	      rmaxnuxt_tab[i] = Nu_out->rmxnux;
	      break;
	    }
	}
	rlumnuef += Nu_lums->rlumnue;
	rlumnuebf += Nu_lums->rlumnueb;
	rlumnuxf += Nu_lums->rlumnux;
	xlumnu += Nu_lums->dlumnu*p->pos[0]/p->r;
	ylumnu += Nu_lums->dlumnu*p->pos[1]/p->r;
	zlumnu += Nu_lums->dlumnu*p->pos[2]/p->r;
	hgw1 += Nu_lums->dlumnu*(1.0+p->pos[2]/p->r)*(2.0*p->pos[0]*p->pos[0]
						 /p->r/p->r - 1.0);
	hgw2 += Nu_lums->dlumnu*(1.0-p->pos[2]/p->r)*(2.0*p->pos[0]*p->pos[0]
						 /p->r/p->r - 1.0);
	hgw3 += Nu_lums->dlumnu*(1.0+
		 sqrtf_fast(1.0-p->pos[2]/p->r*p->pos[2]/p->r)*p->pos[0]/p->r)*
				 (p->pos[2]/p->r*p->pos[2]/p->r - 
				  sqrtf_fast(1.0-p->pos[2]/p->r*p->pos[2]/p->r)*
				  (1.0-p->pos[0]*p->pos[0]/p->r/p->r))/
				 (p->pos[2]/p->r*p->pos[2]/p->r + 
				  sqrtf_fast(1.0-p->pos[2]/p->r*p->pos[2]/p->r)*
				  (1.0-p->pos[0]*p->pos[0]/p->r/p->r));
	hgw4 += Nu_lums->dlumnu*(1.0-
		 sqrtf_fast(1.0-p->pos[2]/p->r*p->pos[2]/p->r)*p->pos[0]/p->r)*
				 (p->pos[2]/p->r*p->pos[2]/p->r - 
				  sqrtf_fast(1.0-p->pos[2]/p->r*p->pos[2]/p->r)*
				  (1.0-p->pos[0]*p->pos[0]/p->r/p->r))/
				 (p->pos[2]/p->r*p->pos[2]/p->r + 
				  sqrtf_fast(1.0-p->pos[2]/p->r*p->pos[2]/p->r)*
				  (1.0-p->pos[0]*p->pos[0]/p->r/p->r));
	hgw5 += Nu_lums->dlumnu*(1.0+
		 sqrtf_fast(1.0-p->pos[2]/p->r*p->pos[2]/p->r)*p->pos[1]/p->r)*
				 (p->pos[2]/p->r*p->pos[2]/p->r - 
				  sqrtf_fast(1.0-p->pos[2]/p->r*p->pos[2]/p->r)*
				  (1.0-p->pos[1]*p->pos[1]/p->r/p->r))/
				 (p->pos[2]/p->r*p->pos[2]/p->r + 
				  sqrtf_fast(1.0-p->pos[2]/p->r*p->pos[2]/p->r)*
				  (1.0-p->pos[1]*p->pos[1]/p->r/p->r));
	hgw6 += Nu_lums->dlumnu*(1.0-
		 sqrtf_fast(1.0-p->pos[2]/p->r*p->pos[2]/p->r)*p->pos[1]/p->r)*
				 (p->pos[2]/p->r*p->pos[2]/p->r - 
				  sqrtf_fast(1.0-p->pos[2]/p->r*p->pos[2]/p->r)*
				  (1.0-p->pos[1]*p->pos[1]/p->r/p->r))/
				 (p->pos[2]/p->r*p->pos[2]/p->r + 
				  sqrtf_fast(1.0-p->pos[2]/p->r*p->pos[2]/p->r)*
				  (1.0-p->pos[1]*p->pos[1]/p->r/p->r));
	/* From Centrella & McMillan (1993) via Fryer, Holz & Hughes (2004) */
	ixx += 2.0/3.0*p->mass*(2.0*p->pos[0]*p->acc[0] - p->pos[1]*p->acc[1] - 
			p->pos[2]*p->acc[2] + 2.0*p->vel[0]*p->vel[0] - 
			p->vel[1]*p->vel[1] - p->vel[2]*p->vel[2]);

	iyy += 2.0/3.0*p->mass*(2.0*p->pos[1]*p->acc[1] - p->pos[2]*p->acc[2] - 
			p->pos[0]*p->acc[0] + 2.0*p->vel[1]*p->vel[1] - 
			p->vel[2]*p->vel[2] - p->vel[0]*p->vel[0]);

	izz += 2.0/3.0*p->mass*(2.0*p->pos[2]*p->acc[2] - p->pos[0]*p->acc[0] - 
			p->pos[1]*p->acc[1] + 2.0*p->vel[2]*p->vel[2] - 
			p->vel[0]*p->vel[0] - p->vel[1]*p->vel[1]);

	ixy += p->mass*(p->pos[0]*p->acc[0] + p->pos[1]*p->acc[1] + 
			2.0*p->vel[0]*p->vel[1]);

	ixz += p->mass*(p->pos[0]*p->acc[0] + p->pos[2]*p->acc[2] + 
			2.0*p->vel[0]*p->vel[2]);

	iyz += p->mass*(p->pos[1]*p->acc[1] + p->pos[2]*p->acc[2] + 
			2.0*p->vel[1]*p->vel[2]);

	p->acc[0] -= Nu_lums->dlumnu*p->pos[0]/p->r/Konst->clight/mass;
	p->acc[1] -= Nu_lums->dlumnu*p->pos[1]/p->r/Konst->clight/mass;
	p->acc[2] -= Nu_lums->dlumnu*p->pos[2]/p->r/Konst->clight/mass;
	enuef += Nu_lums->enue;
	enuebf += Nu_lums->enueb;
	enuxf += Nu_lums->enux;
	e2nuef += Nu_lums->e2nue;
	e2nuebf += Nu_lums->e2nueb;
	e2nuxf += Nu_lums->e2nux;
	enues += Nu_lums->enues;
	enuebs += Nu_lums->enuebs;
	enuxs += Nu_lums->enuxs;
	dee += Nu_lums->dee;
	deeb += Nu_lums->deeb;
	dex += Nu_lums->dex;
	if (Nu_lums->enue < 0.0) Error("Bad value for enue (%g)\n", 
		Nu_lums->enue);
	if (Nu_lums->rlumnue < 0.0) Error("Bad value for enue (%g)\n",
		Nu_lums->rlumnue);
    }
    StopTimer(&sphEOS);
    MPMY_ICombine_Init(&req);
    MPMY_ICombine(&rlumnuef, &rlumnuef, 1, MPMY_FLOAT, MPMY_SUM, req);
    MPMY_ICombine(&rlumnuebf, &rlumnuebf, 1, MPMY_FLOAT, MPMY_SUM, req);
    MPMY_ICombine(&rlumnuxf, &rlumnuxf, 1, MPMY_FLOAT, MPMY_SUM, req);
    MPMY_ICombine(&xlumnu, &xlumnu, 1, MPMY_FLOAT, MPMY_SUM, req);
    MPMY_ICombine(&ylumnu, &ylumnu, 1, MPMY_FLOAT, MPMY_SUM, req);
    MPMY_ICombine(&zlumnu, &zlumnu, 1, MPMY_FLOAT, MPMY_SUM, req);
    MPMY_ICombine(&hgw1, &hgw1, 1, MPMY_DOUBLE, MPMY_SUM, req);
    MPMY_ICombine(&hgw2, &hgw2, 1, MPMY_DOUBLE, MPMY_SUM, req);
    MPMY_ICombine(&hgw3, &hgw3, 1, MPMY_DOUBLE, MPMY_SUM, req);
    MPMY_ICombine(&hgw4, &hgw4, 1, MPMY_DOUBLE, MPMY_SUM, req);
    MPMY_ICombine(&hgw5, &hgw5, 1, MPMY_DOUBLE, MPMY_SUM, req);
    MPMY_ICombine(&hgw6, &hgw6, 1, MPMY_DOUBLE, MPMY_SUM, req);
    MPMY_ICombine(&ixx, &ixx, 1, MPMY_DOUBLE, MPMY_SUM, req);
    MPMY_ICombine(&iyy, &iyy, 1, MPMY_DOUBLE, MPMY_SUM, req);
    MPMY_ICombine(&izz, &izz, 1, MPMY_DOUBLE, MPMY_SUM, req);
    MPMY_ICombine(&ixy, &ixy, 1, MPMY_DOUBLE, MPMY_SUM, req);
    MPMY_ICombine(&ixz, &ixz, 1, MPMY_DOUBLE, MPMY_SUM, req);
    MPMY_ICombine(&iyz, &iyz, 1, MPMY_DOUBLE, MPMY_SUM, req);
    MPMY_ICombine(&enuef, &enuef, 1, MPMY_FLOAT, MPMY_SUM, req);
    MPMY_ICombine(&enuebf, &enuebf, 1, MPMY_FLOAT, MPMY_SUM, req);
    MPMY_ICombine(&enuxf, &enuxf, 1, MPMY_FLOAT, MPMY_SUM, req);
    MPMY_ICombine(&e2nuef, &e2nuef, 1, MPMY_FLOAT, MPMY_SUM, req);
    MPMY_ICombine(&e2nuebf, &e2nuebf, 1, MPMY_FLOAT, MPMY_SUM, req);
    MPMY_ICombine(&e2nuxf, &e2nuxf, 1, MPMY_FLOAT, MPMY_SUM, req);
    MPMY_ICombine(&enues, &enues, 1, MPMY_FLOAT, MPMY_SUM, req);
    MPMY_ICombine(&enuebs, &enuebs, 1, MPMY_FLOAT, MPMY_SUM, req);
    MPMY_ICombine(&enuxs, &enuxs, 1, MPMY_FLOAT, MPMY_SUM, req);
    MPMY_ICombine(&dee, &dee, 1, MPMY_FLOAT, MPMY_SUM, req);
    MPMY_ICombine(&deeb, &deeb, 1, MPMY_FLOAT, MPMY_SUM, req);
    MPMY_ICombine(&dex, &dex, 1, MPMY_FLOAT, MPMY_SUM, req);
    MPMY_ICombine_Wait(req);

    MPMY_Gather(rmaxnuet_tab, NUMRMAXTAB, MPMY_FLOAT, grmaxnuet_tab, 0);
    MPMY_Gather(rmaxnuebt_tab, NUMRMAXTAB, MPMY_FLOAT, grmaxnuebt_tab, 0);
    MPMY_Gather(rmaxnuxt_tab, NUMRMAXTAB, MPMY_FLOAT, grmaxnuxt_tab, 0);
    if (MPMY_Procnum() == 0) {
      qsort(grmaxnuet_tab, grmaxtabsize, sizeof(float), fltcompar);
      qsort(grmaxnuebt_tab, grmaxtabsize, sizeof(float), fltcompar);
      qsort(grmaxnuxt_tab, grmaxtabsize, sizeof(float), fltcompar);
      rmaxnu[0] = grmaxnuet_tab[NUMRMAXTAB-1];
      rmaxnu[1] = grmaxnuebt_tab[NUMRMAXTAB-1];
      rmaxnu[2] = grmaxnuxt_tab[NUMRMAXTAB-1];
    }
    MPMY_Bcast(rmaxnu, 3, MPMY_FLOAT, 0);

    enuef = enuef/(rlumnuef+1e-20);
    e2nuef = e2nuef/(rlumnuef+1e-20);
    enuebf = enuebf/(rlumnuebf+1e-20);
    e2nuebf = e2nuebf/(rlumnuebf+1e-20);
    enuxf = enuxf/(rlumnuxf+1e-20);
    e2nuxf = e2nuxf/(rlumnuxf+1e-20);
    enues = enues/(dee+1e-20);
    enuebs = enuebs/(deeb+1e-20);
    enuxs = enuxs/(dex+1e-20);
    old_nu_lums.enue = enuef;
    old_nu_lums.enueb = enuebf;
    old_nu_lums.enux = enuxf;
    old_nu_lums.e2nue = e2nuef;
    old_nu_lums.e2nueb = e2nuebf;
    old_nu_lums.e2nux = e2nuxf;
    rlumnuef = rlumnuef/xmtheo;
    rlumnuebf = rlumnuebf/xmtheo;
    rlumnuxf = rlumnuxf/xmtheo;
    xlumnu = xlumnu/xmtheo;
    ylumnu = ylumnu/xmtheo;
    zlumnu = zlumnu/xmtheo;
    hgw1 = hgw1/xmtheo;
    hgw2 = hgw2/xmtheo;
    hgw3 = hgw3/xmtheo;
    hgw4 = hgw4/xmtheo;
    hgw5 = hgw5/xmtheo;
    hgw6 = hgw6/xmtheo;
    ixx = ixx/xmtheo;
    iyy = iyy/xmtheo;
    izz = izz/xmtheo;
    ixy = ixy/xmtheo;
    ixz = ixz/xmtheo;
    iyz = iyz/xmtheo;

    {
      double f = Unit2->ufoe/Units->utime;
      singlPrintf(" nue loss: %12g foes/s at <E> (MeV) %g\n", 
		  rlumnuef*f, enuef);
      singlPrintf("nueb loss: %12g foes/s at <E> (MeV) %g\n", 
		  rlumnuebf*f, enuebf);
      singlPrintf(" nux loss: %12g foes/s at <E> (MeV) %g\n", 
		  rlumnuxf*f, enuxf);
      singlPrintf("nusphere nue  losses: %12g foes/s at %g\n", dee*f, enues);
      singlPrintf("nusphere nueb losses: %12g foes/s at %g\n", deeb*f, enuebs);
      singlPrintf("nusphere nux  losses: %12g foes/s at %g\n", dex*f, enuxs);
      singlPrintf("x-dir losses: %12g\n", xlumnu*f);
      singlPrintf("y-dir losses: %12g\n", ylumnu*f);
      singlPrintf("z-dir losses: %12g\n", zlumnu*f);
      singlPrintf("gw emission: %12g,%12g,%12g,%12g,%12g,%12g\n",hgw1*f,hgw2*f,
             hgw3*f, hgw4*f,hgw5*f,hgw6*f);
      singlPrintf("mass motion: %12g,%12g,%12g,%12g,%12g,%12g\n", ixx, iyy,
             izz, ixy, ixz, iyz);
    }

    StartTimer(&sphEOS);
    for (p = btab; p < btab+nobj; p++) {
	if (!SPH_need_update(p)) continue;
	if (p->bghost) continue; /* Is this right? */
	Neut_out->dye = p->dye;
	Neut_out->dynue = p->dynue;
	Neut_out->dynueb = p->dynueb;
	Neut_out->dynux = p->dynux;
	Neut_out->tempnue = p->tempnue;
	Neut_out->tempnueb = p->tempnueb;
	Neut_out->tempnux = p->tempnux;
	Neut_out->enuet = p->enuet;
	Neut_out->enuebt = p->enuebt;
	Neut_out->enuxt = p->enuxt;
	Neut_out->dnue = p->dnue;
	Neut_out->dnueb = p->dnueb;
	Neut_out->dnux = p->dnux;
	Neut_out->dunue =  p->dunue;
	Neut_out->dunueb = p->dunueb;
	Neut_out->dunux = p->dunux;
	Neut_out->dunu = p->dunu;
	Neut_out->etanue = p->etanue;
	Neut_out->etanueb = p->etanueb;
	Neut_out->etanux = p->etanux;
	Neut_out->prnu = p->prnu;
	Nubeta->ebetaeq = p->ebeta;
	Nubeta->pbetaeq = p->pbeta;

	Fortran(neutrino2)(&steps, &p->rho_est, &p->xp, &p->xn, &p->eta,
	   &p->temp, &p->r, &p->mass, &p->vsound, &p->xmuhat, 
           &p->ynue, &p->ynueb, &p->ynux, &p->unue, &p->unueb, &rlumnuef, 
           &rlumnuebf, &rlumnuxf, &enuef,
           &enuebf, &enuxf, &e2nuef, &e2nuebf, &e2nuxf,
           &p->gshift, &rmaxnue, &rmaxnueb, &rmaxnux, &p->h);
	p->dunue = Neut_out->dunue;
	p->dunueb = Neut_out->dunueb;
	p->dunux = Neut_out->dunux;
	p->dynue = Neut_out->dynue;
	p->dynueb = Neut_out->dynueb;
	p->dynux = Neut_out->dynux;
	p->dye = Neut_out->dye;
	p->dunu = Neut_out->dunu;
	dtrapnue += Nutrap->dtrapnue;
	dtrapnueb += Nutrap->dtrapnueb;
	if (!finite(p->dunue)) {
	    SeriousWarning("Bad dunue\n%s\n", PrintSPHBodyContents(p));
	    p->dunue = 0.0;
	}
	if (!finite(p->dynue)) {
	    SeriousWarning("Bad dynue\n%s\n", PrintSPHBodyContents(p));
	    p->dynue = 0.0;
	}
	if (!finite(p->dunu)) {
	    SeriousWarning("Bad dunu\n%s\n", PrintSPHBodyContents(p));
	    p->dunu = 0.0;
	}
	if (!finite(p->dye)) {
	    SeriousWarning("Bad dye\n%s\n", PrintSPHBodyContents(p));
	    p->dye = 0.0;
	}
    }
    StopTimer(&sphEOS);
    MPMY_ICombine_Init(&req);
    MPMY_ICombine(&dtrapnue, &dtrapnue, 1, MPMY_FLOAT, MPMY_SUM, req);
    MPMY_ICombine(&dtrapnueb, &dtrapnueb, 1, MPMY_FLOAT, MPMY_SUM, req);
    MPMY_ICombine_Wait(req);
    dtrapnue=dtrapnue/xmtheo;
    dtrapnueb=dtrapnueb/xmtheo;
    if (dtrapnue > 0.1*rlumnuef) ftrape *= 1.05;
    if (dtrapnue < 0.03*rlumnuef && ftrape > 1.05) ftrape *= 0.97;
    if (dtrapnueb > 0.1*rlumnuebf) ftrapb *= 1.05;
    if (dtrapnueb < 0.03*rlumnuebf && ftrapb > 1.05) ftrapb *= 0.97;
    singlPrintf("dtrapnue: %12g dtrapnueb: %12g ftrape: %f ftrapb: %f ftrapx: %f\n", dtrapnue, dtrapnueb, ftrape, ftrapb, ftrapx);
}

#include "physics.h"

void
update_point_SPHmass(SPHbody *btab, int SPHnobj, 
		  void *pp, float smooth2, float newt)
{
    SPHbody *r;
    body *p = pp;
    float dr2, oneor, oneor2;
    float phii;
    Vxd(float r);
    Vxd(float ppos);

    VxV(ppos, = p->pos);

    for (r = btab; r < btab+SPHnobj; r++) {
	VxVVx(r, = r->pos, - ppos); /* 3 flops */
	
	dr2 = Dotx(r, r);	/* 5 flops */
	if (dr2 != (float)0.0) {	
	  dr2 += smooth2;
	
	  oneor = recipsqrtf(dr2);	/* 8 flops */
	
	  oneor2 = oneor * oneor;	/* 17 flops */
	  phii = newt * oneor * p->mass;
	  r->phi -= phii;
	  VVx(r->acc, -= oneor2 * phii * r);
	  phii = newt * oneor * r->mass;
	  p->phi -= phii;
	  VVx(p->acc, += oneor2 * phii * r);
	}
    }
}


void
update_bardeen(SPHbody *btab, int nobj, float G, float c, bndry_t b)
{
    /* From Nelson & Papaloizou (2000), eqs. 8 and 14 */
    SPHbody *p;
    Vxd(float r);
    Vxd(float v);
    Vxd(float S);
    Vxd(float ppos);
    Vxd(float pvel);
    float dr2, oneor, A, B, C;
    float rplus = G*b.mass/(c*c);

    VxS(ppos, = (float)0.0);  /* BH fixed at origin */
    VxS(pvel, = (float)0.0);
    VxV(S, = G/(c*c)*b.j);

    for (p = btab; p < btab+nobj; p++) {

	VxVVx(v, = p->vel, - pvel);
	VxVVx(r, = p->pos, - ppos);
	dr2 = Dotx(r, r);
	oneor = recipsqrtf(dr2);

 	A = -G*b.mass*oneor*oneor*oneor;

 	p->acc[0] += A*r0;
 	p->acc[1] += A*r1;
 	p->acc[2] += A*r2;

	p->phi -= G*b.mass*oneor;

/*  	A = -G*b.mass*oneor*oneor*oneor*(1.0 + 6.0*rplus*oneor); */
/*  	B = 2.0*oneor*oneor*oneor; */
/*  	C = 6.0*Dotx(S, r)*oneor*oneor*oneor*oneor*oneor; */

/*  	p->acc[0] += A*r0 + B*(v1*S2 - v2*S1) + C*(r1*v2 - r2*v1); */
/*  	p->acc[1] += A*r1 + B*(v2*S0 - v0*S2) + C*(r2*v0 - r0*v2); */
/*  	p->acc[2] += A*r2 + B*(v0*S1 - v1*S0) + C*(r0*v1 - r1*v0); */
    }
}


void 
do_SPHgrav(const float *p, const float *end, const float *pos0, float *mass0, 
	   float *acc0, float *phi0, const float *eps2p, int *ncut)
{
    float dr2, h, h2, dxx, dphidx, dmassdx, v2;
    Vxd(float r);
    float phii, mor3, mass;
    VxdV(const float ppos, = pos0);
    Vxd(float a);
    float phi = *phi0;
    float total_mass = *mass0;
    const float eps2 = *eps2p; /* Nuts! Eps2 is not eps^2, but hsink!*/
    int index;

    VxV(a, = acc0);

    while (p < end) {
	mass = *p++;
	r0 = *p++;
	r1 = *p++;
	r2 = *p++;
	h = *p++;

	h=(h+eps2)/2.0;

	h2= h * h;
	VxVx(r, -= ppos);	    /* 3 flops */
	dr2 = Dotx(r, r);	    /* 5 flops */
	total_mass += mass; /* Hmm, do I need total mass or enclosed
			       here?  If enclosed, put this statement
			       after if clause*/
	
	if (dr2 >= 4.*h2) { /* Beyond 2h, point source for phi and acc! */
	    phii = recipsqrtf(dr2); /* 8 flops */
	    mor3 = phii * phii;	    /* 5 flops */
	    phii *= mass;
	    mor3 *= phii;
	} 
	else if (dr2 > (float) 0.) { /* Within 2h, use SPH particle
					smoothing for gravity */
	    v2 = dr2 / h2;
	    index = v2 * invdvtable;
	    phii=fpoten[index]*mass/h;
	    mor3=mass*fmass[index]/dr2*recipsqrtf(dr2);
	}
	phi -= phii;
	VxVx(a, += mor3 * r); /* 6 flops */
    }
    VVx(acc0, = a);
    *mass0 = total_mass;
    *phi0 = phi;
}
