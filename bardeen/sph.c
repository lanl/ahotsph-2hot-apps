#include <math.h>
#include <stdlib.h>
#include "physics_sph.h"
#include "vop.h"
#include "fastflpt.h"
#include "timers.h"
#include "error.h"

#ifndef M_1_PI
#define	M_1_PI 0.31830988618379067154
#endif
#define NKERNEL_TABLE 40000
#define MAX_INDEX (NKERNEL_TABLE+2)

#define NO_UPDATE 2

Counter_t SPHCnt, SPHrej, nbrMACCnt;

static float dvtable; /* == 0.0001 ... */
static float invdvtable; /* == 10000.0 ... */
static float cnormk;
static float wij[MAX_INDEX];
static float grwij[MAX_INDEX];
static float fmass[MAX_INDEX];
static float fpoten[MAX_INDEX];
static float Gamma = (float)(4.0/3.0);  /* Yikes; remember this */
static float alpha = (float)1.0;
static float beta = (float)2.5;
static float epsil = (float)1e-2;
static float heatf1 = (float)1.0;
static int ndim;
static int Nobj;
static int add_offset;
static float offset[NDIM];
static float voffset[NDIM];
static void (*bodyfunc)(SinkSPH *sink, hcell **src_vec, int *res, int n);
static void (*cellfunc)(SinkSPH *sink, hcell **src_vec, int *res, int n);

extern int do_diffusion;
extern int do_adiabatic;
extern float eos_K;

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
InheritSPH(const SinkSPH *from, SinkSPH *to, hcell *pp)
{
    if( to == NULL ){
	SPHbody *bp = pp->ptr;
	if (from->isbody == NO_UPDATE)
	  return;
	/* Must accumulate for periodic BC to work */
	/* Must initialize to zero appropriately */
	bp->rho += from->rho;
	bp->du += from->du;       /* Why am I doing this? */
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
	to->rho = (float)0.0;
	to->drho_dt = (float)0.0;
	to->udot = (float)0.0;
	VS(to->lvel, = (float)0.0);
	to->nterms = 1;
	to->nbrs = 0;
	VS(to->M1, = (float)0.0);
	to->min_nbr_dt = 1e30;
	/* Diffusion quantities */
	/* Can some of these be set to zero here?  Like udot above? */
	to->temp = bp->temp;
	to->du = bp->du;
	to->u_r = bp->u_r;
	to->du_r = bp->du_r;  /* Or = (float)0.0; ? */
	to->D = bp->D;
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
	dr2 = Dotx(r, r); /* == 400.0000000032623 */

	if ((rij = sqrtf_fast(dr2)) > extent_src + extent_sink
	    || dr2 == (float)0.0) {
	    goto accept;
	} else if (daughters != 1) {
	    goto failed;
	}

	hmean11 = (float)2.0 / (h + bp->h);
	hmean21 = hmean11 * hmean11; /* == 0.01 */
	    
	v2 = dr2 * hmean21;
	index = v2 * invdvtable;
	if (index >= MAX_INDEX) Error("Index too large\n");
	dxx = v2 - index * dvtable;
	dwdx = (wij[index+1] - wij[index]) * invdvtable;
	wtij = (wij[index] + dwdx * dxx ) * hmean21 * hmean11;
	if (wtij < (float)0.0) Error("Negative wtij (macRho) = %g\n", wtij);
	dgrwdx = (grwij[index+1] - grwij[index]) * invdvtable;
	grwtij = (grwij[index] + dgrwdx * dxx) * hmean21 * hmean21;

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
    sink->rho += rhoi;
    sink->nbrs += nbrs;
    sink->drho_dt -= divvi;
}


void
macSPH(SinkSPH *sink, hcell **source_vec, int *result, int n)
{
    const float extent_sink = sink->extent;
    VxdV(const float pos_sink, = sink->pos);
    VxdV(const float v, = sink->vel);
    const float h = sink->h;
    const float pro2 = sink->pr / (sink->rho_est * sink->rho_est);
    const float mass = sink->mass;
    const float rho_est = sink->rho_est;
    const float vsound = sink->vsound;
    const float u = sink->u;
    Vxd(float r);
    Vxd(float f);
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
    float rhoi = (float)0.0;
    float divvi = (float)0.0;
    float dr2;
    Vxd(float runi);
    float dq = (float)0.0;
    float vv, vv2;
    float dxx, dwdx, wtij, dgrwdx, grwtij;
    float rapm, robar1, grpm, wpm;
    float poro2;
    float projv, vsbar, est_divv, t12;
    float rij, rij1;
    int interactions = 0;

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
	wtij = (wij[index] + dwdx * dxx) * hmean21 * hmean11;
	if (wtij < (float)0.0) Error("Negative wtij (macSPH) = %g\n", wtij);
	dgrwdx = (grwij[index+1] - grwij[index]) * invdvtable;
	grwtij = (grwij[index] + dgrwdx * dxx) * hmean21 * hmean21;

	rapm = mass / bp->mass;
	robar1 = (float)2.0 / (rho_est + bp->rho_est);
	grpm = bp->mass * grwtij;
	wpm = bp->mass * wtij;

	poro2 = grpm * (pro2 + bp->pr / (bp->rho_est * bp->rho_est));
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
	/* flux-limited diffusion */
	if (do_diffusion)
	    if (grpm < 0.0) {  /* What does this condition really mean? */
		float Dmeanr = 2.0*rij1*sink->D/(sink->D+bp->D)*bp->D;

		sink->du_r += ( (C_LIGHT < Dmeanr) ? C_LIGHT : Dmeanr ) *
		    (sink->u_r - bp->u_r) * grpm / bp->rho_est;
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
    sink->rho += rhoi;
    sink->drho_dt -= divvi;
    sink->udot += (float)0.5 * dq;
    VVx(sink->M1, += f);
    VVx(sink->lvel, += smv);
    sink->nbrs += nbrs;
    sink->nterms += nbrs*8;
    sink->min_nbr_dt = min_nbr_dt;
}

void
SetSPH(float visc_alpha, float visc_beta, float visc_epsilon, float heat_f1,
       float eos_gamma, int gnobj,  void bfunc(), void cfunc())
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
SPH_setup(int dim)
{
    float v2max;
    float v, v2, v3;
    float dif2;
    double sum;   /* sum needs to be double */
    int i, i1;

    ndim = dim;


    /* Do any diffusion-specific initialization here: */
    /* Set rmax, luminosity, etc. */


    /* maximum interaction length and step size */
    v2max = (float)4.0;
    dvtable = v2max / NKERNEL_TABLE;
    i1 = (float)1.0 / dvtable;
    invdvtable = i1;

    /* normalisation constant */
    if (ndim == 3)
      cnormk = (float)M_1_PI;
    else if (ndim == 2)
      cnormk = (float)(M_1_PI*10.0/7.0);
    else if (ndim == 1)
      cnormk = (float)(2.0/3.0);
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
    }

    /*  b) v greater than 1 */
    for (i = i1; i <= NKERNEL_TABLE; i++) {
	v2 = i * dvtable;
	v = sqrtf_fast(v2);
	v3 = v * v2;
	dif2 = (float)2.0 - v;
	sum = 0.25 * dif2 * dif2 * dif2;
	wij[i] = cnormk * sum;
	sum = -0.75 * v2 + 3.0 * v - 3.0;
	grwij[i] = cnormk * sum;
	fmass[i] = (-1.0/6.0)*v3*v3 + 1.2*v2*v3 - 3.0*v2*v2 + (8.0/3.0)*v3
	  - (1.0/15.0);
	fpoten[i] = (-1.0/30.0)*v2*v3 + 0.3*v2*v2 - v3 + (4.0/3.0)*v2 - 1.6;
    }
    /* For interpolation of last table entry */
    wij[NKERNEL_TABLE+1] = grwij[NKERNEL_TABLE+1] = (float) 0.0;
}

void
update_final(SPHbody *btab, int nobj, float dt, int *limit_high, int *limit_low)
{
    SPHbody *p;

    for (p = btab; p < btab+nobj; p++) {
	if (!SPH_need_update(p)) continue;
	VV(p->acc, += p->grav_acc);
	p->rho += cnormk * p->mass / (p->h * p->h * p->h);
	p->hdot = (float)(-1.0/3.0) * p->h * p->drho_dt / p->rho;
	if (p->hdot * dt > p->h) {
	    SeriousWarning("Hdot limit (high)\n%s\n", PrintSPHBodyContents(p));
	    p->hdot = p->h/dt;
	    ++*limit_high;
	}
	if (p->hdot * dt < -0.5*p->h) {
	    SeriousWarning("Hdot limit (low)\n%s\n", PrintSPHBodyContents(p));
	    p->hdot = -0.5*p->h/dt;
	    ++*limit_low;
	}

	if (!do_adiabatic) {
	    p->udot += p->drho_dt * p->pr / (p->rho * p->rho) + 
		( (do_diffusion) ? (p->du_r/p->rho) /* Diffusion */
		  : 0.0 );

	    if (!finite(p->udot)) 
		Error("Bad value for udot\n");

	    /* Are these limits appropriate? */
	    /* Does this enforce the Courant limit correctly with diffusion? */
	    if ( (p->udot * dt > p->u) && !(p->ident & (1<<30)) ) {
		p->udot = p->u/dt;
		++*limit_high;
	    }
	    if ( (p->udot * dt < -0.333*p->u) && !(p->ident & (1<<30)) ) {
		p->udot = -0.333*p->u/dt;
		++*limit_low;
	    }
	}
    }
}


double eos_n, eos_u;

void
update_intermediate(SPHbody *btab, int nobj, float dt_last, int flag, int *limit)
{
    float kes, kff;  /* Opacities (Thomson, free-free) */
    SPHbody *p;

    for (p = btab; p < btab+nobj; p++) {
	if (!SPH_need_update(p)) continue;
	if (flag) p->rho_est = p->rho + p->drho_dt * dt_last;
	else p->rho_est = p->rho;
	if (p->rho_est <= (float)0.0) 
	  Error("Rho_est is 0\n%s\n", PrintSPHBodyContents(p));

	if (do_adiabatic) 
	    p->pr = eos_K * pow(p->rho_est, Gamma);
	else 
	    p->pr = p->u * (Gamma - (float)1.0) * p->rho_est;

	p->vsound = sqrtf_fast(Gamma * p->pr / p->rho_est);

	if (do_diffusion) {
	    /* Set constants in physics_sph.h for now */
	    /* Or read in from the control file (global)? */

	    /* Calculate temperature from u, then "create" photons (a*T^4) */
	    eos_n = ((double)(p->rho_est))/((double)(MH));
	    eos_u = ((double)(p->u))*((double)(p->rho_est));

	    /* Figure out good upper and lower limits for temp */
	    p->temp = newtraph(4.0e4, 1.5e7, eos_u*1.0e-6, uvst, duvst);
	    p->u_r = A_COEFF*p->temp*p->temp*p->temp*p->temp;
	    p->du_r = 0.0;

	    /* Calculate diffusion coefficient */
	    kes = KES_COEFF;
	    kff = (KFF_COEFF)*p->rho_est*pow(p->temp, -3.5);
	    p->D = C_LIGHT / ( 3.0*(kes+kff)*p->rho_est );

	    /* Also, eventually, handle lightbulb approximation here */
	}
    }
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
update_point_SPHmass2(SPHbody *btab, int SPHnobj, float smooth2, float newt, 
		      float mass)
{
    SPHbody *r;
    float dr2, oneor, oneor2;
    float phii;
    Vxd(float r);
    Vxd(float ppos);

    VxS(ppos, = (float)0.0);  /* Body fixed at origin */

    for (r = btab; r < btab+SPHnobj; r++) {
	VxVVx(r, = r->pos, - ppos); /* 3 flops */
	
	dr2 = Dotx(r, r);	/* 5 flops */
	if (dr2 != (float)0.0) {	
	  dr2 += smooth2;
	
	  oneor = recipsqrtf(dr2);	/* 8 flops */
	
	  oneor2 = oneor * oneor;	/* 17 flops */
	  phii = newt * oneor * mass;
	  r->phi -= phii;
	  VVx(r->acc, -= oneor2 * phii * r);
	}
    }
}


void
update_point_SPHmass3(SPHbody *btab, int SPHnobj, float smooth2, float newt, 
		      float mass, float b)
{
    /* Plummer model */
    SPHbody *r;
    float dr2b2, oneor, oneor2;
    float phii;
    Vxd(float r);
    Vxd(float ppos);

    VxS(ppos, = (float)0.0);  /* Body fixed at origin */

    for (r = btab; r < btab+SPHnobj; r++) {

	VxVVx(r, = r->pos, - ppos);
	
	dr2b2 = Dotx(r, r) + b*b;

	oneor = recipsqrtf(dr2b2);
	
	oneor2 = oneor * oneor;
	phii = newt * oneor * mass;
	r->phi -= phii;
	VVx(r->acc, -= oneor2 * phii * r);

    }
}


void
update_point_mass_bardeen(SPHbody *btab, int nobj, float G, 
			  float mass, float rplus, float S)
{
    /* From Nelson & Papaloizou (2000), eqs. 8 and 14 */
    SPHbody *p;
    Vxd(float r);
    Vxd(float v);
    Vxd(float ppos);
    Vxd(float pvel);
    float dr2, oneor, A, B, C;

    VxS(ppos, = (float)0.0);  /* BH fixed at origin */
    VxS(pvel, = (float)0.0);

    for (p = btab; p < btab+nobj; p++) {

	VxVVx(v, = p->vel, - pvel);
	VxVVx(r, = p->pos, - ppos);
	dr2 = Dotx(r, r);
	oneor = recipsqrtf(dr2);

 	A = -G*mass*oneor*oneor*oneor*(1.0 + 6.0*rplus*oneor);
 	B = 2.0*S*oneor*oneor*oneor;
 	C = 6.0*S*r2*oneor*oneor*oneor*oneor*oneor;

 	p->acc[0] += A*r0 + B*v1 + C*(r1*v2 - r2*v1);
 	p->acc[1] += A*r1 - B*v0 + C*(r2*v0 - r0*v2);
 	p->acc[2] += A*r2 + C*(r0*v1 - r1*v0);

/*  	A = -G*mass*oneor*oneor*oneor; */

/*  	p->acc[0] += A*r0; */
/*  	p->acc[1] += A*r1; */
/*  	p->acc[2] += A*r2; */
    }
}
