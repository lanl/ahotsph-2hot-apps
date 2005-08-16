#include <math.h>
#include "physics_sph.h"
#include "vop.h"
#include "fastflpt.h"
#include "timers.h"
#include "error.h"
#include "singlio.h"
#include "mpmy.h"

#ifndef M_1_PI
#define	M_1_PI 0.31830988618379067154
#endif
#define NKERNEL_TABLE 80000
#define MAX_INDEX (NKERNEL_TABLE+2)

#define NO_UPDATE 2

extern Timer_t sphEOS;
Counter_t SPHCnt, SPHrej, nbrMACCnt;

konst_s *konst;
output_s *output;
units_s *units;
unit2_s *unit2;

static float dvtable;
static float invdvtable;
static float cnormk;
static float wij[MAX_INDEX];
static float grwij[MAX_INDEX];
static float Gamma = (5.0/3.0);
static float alpha = 1.0;
static float beta = 2.5;
static float epsil = 1e-2;
static float heatf1 = 1.0;
static float rinner;

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
	to->gshift = bp->gshift;
    }

    if (add_offset) {
	VV(to->pos, += offset);
	VV(to->vel, += voffset);
    }
    if (add_rotate) {
        to->pos[0] = to->pos[0] * cosa - to->pos[1] * sina;
        to->pos[1] = to->pos[0] * sina + to->pos[1] * cosa;
        to->vel[0] = to->vel[0] * cosa - to->vel[1] * sina;
        to->vel[1] = to->vel[0] * sina + to->vel[1] * cosa;
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
	grwtij = (grwij[index] + dgrwdx * dxx) * hmean21 * hmean11;
#endif

	/* if (bp->mass * wtij < 0.0) Error("rhoi < 0.0\n"); */
	rhoi += bp->mass * wtij;

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
    sink->rho += rhoi * sink->xfac;
    sink->nbrs += nbrs;
    sink->drho_dt -= divvi * sink->xfac;
}


void
macSPH(SinkSPH *sink, hcell **source_vec, int *result, int n)
{
    const float extent_sink = sink->extent;
    VxdV(const float pos_sink, = sink->pos);
    VxdV(const float v, = sink->vel);
    const float h = sink->h;
    const double pro2 = (sink->pr / sink->rho_est) / sink->rho_est;
    const float mass = sink->mass;
    const float rho_est = sink->rho_est;
    const float vsound = sink->vsound;
    const float u = sink->u;

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
    double dq = (float)0.0;
    float vv, vv2;
    float dxx, dwdx, dgrwdx;
    double wtij, grwtij;
    float rapm, robar1;
    double wpm, grpm;
    double poro2;
    double projv, vsbar, est_divv, t12;
    float rij, rij1;
    int interactions = 0;

    konst = (void *)&Fortran(konst);
    unit2 = (void *)&Fortran(unit2);

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
#else
	wtij = (wij[index] + dwdx * dxx) * hmean21;
	grwtij = (grwij[index] + dgrwdx * dxx) * hmean21 * hmean11;
#endif

	rapm = mass / bp->mass;
	robar1 = (float)2.0 / (rho_est + bp->rho_est);	
	grpm = bp->mass * grwtij;
	wpm = bp->mass * wtij;

	poro2 = grpm * (pro2 + (bp->pr / bp->rho_est) / bp->rho_est);
	rij1 = (float)1.0 / rij;
	VxVx(runi, = rij1 * r);
	VxVx(f, += poro2 * runi);
	VxVVx(dv, = bp->vel, - v);
	VxVx(smv, += robar1 * wpm * dv);
	projv = Dotx(dv, runi);

	rhoi += bp->mass * wtij;
	divvi -= bp->mass * grwtij * projv;

	/* artificial viscosity and energy dissipation */
	vsbar = (float)0.5 * (vsound + bp->vsound);
	if (projv < (float)0.0 && alpha != (float)0.0) {
	    est_divv = projv * vv / (vv2 + epsil);
	    t12 = grpm * est_divv * (beta * est_divv - alpha * vsbar) * robar1;
	    VxVx(f, += t12 * runi);
	    dq += t12 * projv;
	}
	/* artificial heat conduction */
	if (heatf1 != (float)0.0) {
	    float qmean = heatf1 * vsbar / hmean11;
	    t12 = grpm * qmean * (bp->u - u) * robar1;
	    dq -= t12;
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
    sink->interactions += interactions;
    sink->rho += rhoi * sink->xfac;
    sink->drho_dt -= divvi * sink->xfac;
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
    beta = visc_beta;
    epsil = visc_epsilon;
    heatf1 = heat_f1;
    Gamma = eos_gamma;
    bodyfunc = bfunc;
    cellfunc = cfunc;
}

void
SPHaux(float rb)
{
  rinner = rb;
}
void
SPH_setup(int dim, int ncoef1, double *wcoef1, int ncoef2, double *wcoef2)
{
    double v2max;
    double v, v2;
    double w, dw;
    double ddvtable;
    int i, i1, j;

    ndim = dim;

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
    grwij[0] = 0.0;		/* Is this right? */
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
	grwij[i] = cnormk * dw / v;
    }

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
	grwij[i] = cnormk * dw / v;
    }
    /* For interpolation of last table entry */
    wij[NKERNEL_TABLE+1] = grwij[NKERNEL_TABLE+1] = 0.0;

#ifdef WRITE_WIJ
    for (i = 0; i <= NKERNEL_TABLE; i++) 
      singlPrintf("%5d %8g %8g\n", i, wij[i], grwij[i]);
#endif
}

void
update_final(SPHbody *btab, int nobj, float dt, int *limit_high, int *limit_low)
{
    SPHbody *p;

    for (p = btab; p < btab+nobj; p++) {
	if (!SPH_need_update(p)) continue;
	VV(p->acc, += p->grav_acc);
	/* more changes here -1.0 -> -100.0 */
	p->hdot = (float)(-1.0/3.0) * p->h * p->drho_dt / p->rho;
	/* more changes here 1. -> 200. */
	if (p->hdot * dt > p->h) {
	    SeriousWarning("Hdot limit (high)\n%s\n", PrintSPHBodyContents(p));
	    p->hdot = p->h/dt;
	    ++*limit_high;
	}
	/* more changes here -0.5 -> -100.0 */
	if (p->hdot * dt < -0.5*p->h) {
	    SeriousWarning("Hdot limit (low)\n%s\n", PrintSPHBodyContents(p));
	    p->hdot = -0.5*p->h/dt;
	    ++*limit_low;
	}
	p->udot += (p->drho_dt / p->rho) * (p->pr / p->rho);
	if (!finite(p->udot)) {
	    Error("Bad value for udot: #%d: drho_dt: %e; pr: %e; rho; %e; nbrs: %d\n", p->ident, p->drho_dt, p->pr, p->rho, p->nbrs);
	}
	p->udot2 /= p->temp;
	if (!finite(p->udot2)) Error("Bad value for udot2\n");
	if (p->temp/p->temprev > 1.05) {
	    SeriousWarning("Tempdot limit (high)\n%s\n", 
			   PrintSPHBodyContents(p));
	    *limit_high += 4;
	}
#if 0
	if (p->udot * dt > p->u) {
	    SeriousWarning("Udot limit (high)\n%s\n", PrintSPHBodyContents(p));
	    p->udot = p->u/dt;
	    ++*limit_high;
	}
	if (p->udot * dt < -0.5*p->u) {
	    SeriousWarning("Udot limit (low)\n%s\n", PrintSPHBodyContents(p));
	    p->udot = -0.5*p->u/dt;
	    ++*limit_low;
	}
#endif
	if (p->temp > 500.) {
	  singlPrintf("udot: %12g udot2: %12g\n",p->udot,p->udot2);
	  singlPrintf("drho_dt: %12g rho: %12g\n",p->drho_dt,p->rho);
	}

    }
}

void
update_intermediate(SPHbody *btab, int nobj, float dt_last, int flag, int *limit)
{
    SPHbody *p;
    double rho, u, udot, u2, ye, temp, abar, xp, xn, ufreez;
    double xpf, p2, p3, p4, temprev, rhoprev, xpprev, xnprev, yeprev;
    double dtbrn,xhe4n,xc12n,xo16n,xne20n,xmg24n,xsi28n;
    double xni56n,xco56n,xfe56n;
    int ifleos;
    
    float xmtheo, steps, mass;

    StartTimer(&sphEOS);
    output = (void *)&Fortran(output);
    units = (void *)&Fortran(units);
    unit2 = (void *)&Fortran(unit2);
    xmtheo = 1;
    steps = dt_last;
    
    for (p = btab; p < btab+nobj; p++) {
	if (!SPH_need_update(p)) continue;
	if (flag) p->rho_est = p->rho + p->drho_dt * dt_last;
	else p->rho_est = p->rho;
	if (p->rho_est <= (float)0.0) 
	  Error("Rho_est is 0\n%s\n", PrintSPHBodyContents(p));
	rho = p->rho_est;
	u = p->u;
	udot = p->udot;
	u2 = p->u2;
	ye = p->ye;
	temp = p->temp;
	ifleos = p->ifleos;
	abar = p->abar;
	xp = p->xp;
	xn = p->xn;
	xpf = p->xpf;
	p2 = p->p2;
	p3 = p->p3;
	p4 = p->p4;
	mass=p->mass;
	temprev = p->temprev;
	rhoprev = p->rhoprev;
	xpprev = p->xpprev;
	xnprev = p->xnprev;
	yeprev = p->yeprev;
	ufreez = p->ufreez;
	dtbrn = p->dtbrn;
	xhe4n = p->xhe4n;
	xc12n = p->xc12n;
	xo16n = p->xo16n;
	xne20n = p->xne20n;
	xmg24n = p->xmg24n;
	xsi28n = p->xsi28n;
	xni56n = p->xni56n;
	xco56n = p->xco56n;
	xfe56n = p->xfe56n;
	if (xhe4n < 1.e-10) {
	  xhe4n = 0;
	}
	if (xc12n < 1.e-10) {
	  xc12n = 0;
	}
	if (xo16n < 1.e-10) {
	  xo16n = 0;
	}
	if (xne20n < 1.e-10) {
	  xne20n = 0;
	}
	if (xmg24n < 1.e-10) {
	  xmg24n = 0;
	}
	if (xsi28n < 1.e-10) {
	  xsi28n = 0;
	}
	if (xni56n < 1.e-10) {
	  xni56n = 0;
	}
	if (xco56n < 1.e-15) {
	  xco56n = 0;
	}
	if (xfe56n < 1.e-15) {
	  xfe56n = 0;
	}
	Fortran(eos3)(&steps, &mass, &rho, &u, &udot, &u2, &ye, &temp, 
		      &ifleos, &abar, &xp, &xn, &xpf, &p2, &p3, &p4, 
		      &temprev, &rhoprev, &xpprev, &xnprev, &yeprev, 
		      &ufreez, &xhe4n, &xc12n, &xo16n, &xne20n, &xmg24n, 
		      &xsi28n, &xni56n, &xco56n, &xfe56n, &dtbrn);  
	if (rho != p->rho_est) Error("rho mismatch\n");
	if (ye != p->ye) Error("ye mismatch\n");
	if (mass != p->mass) Error("mass mismatch\n");
	p->u=u;
	p->ifleos = ifleos;
	p->u2 = u2;
	temp = (float)temp;
	p->temp = temp;
	abar = (float)abar;
	p->abar = abar;
	xp = (float)xp;
	p->xp = xp;
	xn = (float)xn;
	p->xn = xn;
	xhe4n = (float)xhe4n;
	p->xhe4n = xhe4n;
	xc12n = (float)xc12n;
	p->xc12n = xc12n;
	xo16n = (float)xo16n;
	p->xo16n = xo16n;
	xne20n = (float)xne20n;
	p->xne20n = xne20n;
	xmg24n = (float)xmg24n;
	p->xmg24n = xmg24n;
	xsi28n = (float)xsi28n;
	p->xsi28n = xsi28n;
	xni56n = (float)xni56n;
	p->xni56n = xni56n;
	xco56n = (float)xco56n;
	p->xco56n = xco56n;
	xfe56n = (float)xfe56n;
	p->xfe56n = xfe56n;
	dtbrn = (float)dtbrn;
	p->dtbrn = dtbrn;
	output->vsound = (float)output->vsound;
	p->vsound = output->vsound;
	output->pr = (float)output->pr;
	p->pr = output->pr;
	output->eta = (float)output->eta;
	p->eta = output->eta;
	output->xmuhat = (float)output->xmuhat;
	p->xmuhat = output->xmuhat;
	output->xmue = (float)output->xmue;
	p->xmue = output->xmue;
	p->xpf = xpf;
	p->p2 = p2;
	p->p3 = p3;
	p->p4 = p4;
	p->temprev = temprev;
	p->rhoprev = rhoprev;
	p->xpprev = xpprev;
	p->xnprev = xnprev;
	p->yeprev = yeprev;
	p->ufreez = ufreez;
	}
    StopTimer(&sphEOS);
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
update_point_SPHmass_bndry(SPHbody *btab, int SPHnobj, float newt, 
			   bndry_t bndry)
{
    SPHbody *r;
    float dr2, oneor, oneor2;
    float phii;
    Vxd(float r);

    for (r = btab; r < btab+SPHnobj; r++) {

	VxVV(r, = r->pos, - bndry.pos);
	
	dr2 = Dotx(r, r);

	oneor = recipsqrtf(dr2);
	
	oneor2 = oneor * oneor;
	phii = newt * oneor * bndry.mass;
	r->phi -= phii;
	VVx(r->acc, -= oneor2 * phii * r);
    }
}
