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

Counter_t SPHCnt, SPHrej, nbrMACCnt;

static float dvtable;
static float invdvtable;
static float cnormk;
static float wij[MAX_INDEX];
static float grwij[MAX_INDEX];
static float Gamma = (float)(5.0/3.0);
static float alpha = (float)1.0;
static float beta = (float)2.5;
static float epsil = (float)1e-2;
static int ndim;
static int Nobj;
static void (*bodyfunc)(SinkSPH *sink, hcell **src_vec, int *res, int n);
static void (*cellfunc)(SinkSPH *sink, hcell **src_vec, int *res, int n);

void 
InheritSPH(const SinkSPH *from, SinkSPH *to, hcell *pp)
{
    if( to == NULL ){
	body *bp = pp->ptr;
	bp->rho = from->rho;
	bp->drho_dt = from->drho_dt;
	bp->udot = from->udot;
	bp->nbrs = from->nbrs;
	VV(bp->acc, += from->M1);
	/* bp->nterms = from->nterms; */
	if (from->interactions != Nobj)
	    Error("Ninteract is %d, should be %d\n", from->interactions, Nobj);
	return;
    }

    if (Sub_Flags(pp)) {
	/* Stuff needed for cell-cell neighbor evaluation */
	cell *cp = pp->ptr;
	VV(to->pos, = cp->pos);
	to->extent = cp->bmax + cp->lap;
	to->isbody = 0;
    } else {
	/* Stuff needed for above, plus physics info */
	body *bp = pp->ptr;
	VV(to->pos, = bp->pos);
	to->extent = bp->h;
	to->h = bp->h;
	to->isbody = 1;
	VV(to->vel, = bp->vel);
	to->pr = bp->pr;
	to->rho_est = bp->rho_est;
	to->mass = bp->mass;
	to->vsound = bp->vsound;
	to->rho = (float)0.0;
	to->drho_dt = (float)0.0;
	to->udot = (float)0.0;
	/* to->nterms = 0; */
	to->nbrs = 0;
	VS(to->M1, = (float)0.0);
   to->alfa = bp->alfa;
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
    if (sink->isbody)
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
    body *bp = 0;
    float projv;
    float v2;
    float rij;
    float wtij, grwtij;
    float hmean11, hmean21, hmean31, hmean41;
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
	    const cell *cp = source->ptr;
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
	hmean31 = hmean21 * hmean11;
	hmean41 = hmean21 * hmean21;
	    
	v2 = dr2 * hmean21;
	index = v2 * invdvtable;
	if (index >= MAX_INDEX) Error("Index too large\n");
	dxx = v2 - index * dvtable;
	dwdx = (wij[index+1] - wij[index]) * invdvtable;
	wtij = (wij[index] + dwdx * dxx ) * hmean31;
	if (wtij < (float)0.0) Error("Negative wtij = %g\n", wtij);
	dgrwdx = (grwij[index+1] - grwij[index]) * invdvtable;
	grwtij = (grwij[index] + dgrwdx * dxx) * hmean41;

	rhoi += bp->gr_mass * wtij;    

	/* velocity divergence times density */
	VxVV(dv, = bp->vel, - sink->vel);
	projv = grwtij * Dotx(dv, r) * recipsqrtf(dr2);
	divvi -= bp->gr_mass * projv;  

	nbrs++;
      accept:
	interactions += daughters;
	result[i] = MAC_ACCEPT;
	continue;
      failed:
	result[i] = MAC_SPLIT_SRC;
    }

    sink->interactions += interactions;
    /* sink->nterms += (float)nterms; */
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
    Vxd(float r);
    Vxd(float f);
    Vxd(float dv);
    float extent_src;
    int daughters;
    int i;
    body *bp = 0;
    float hmean11, hmean21, hmean31, hmean41;
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
    for (i = 0; i < n; i++) {
	const hcellptr source = source_vec[i];

	if (Sub_Flags(source)) {
	    const cell *cp = source->ptr;
	    VxV(r, = cp->pos);	
	    extent_src = cp->bmax + cp->lap;
	    daughters = cp->daughters;
	} else {
	    bp = source->ptr;
	    VxV(r, = bp->pos);
	    extent_src = bp->h;
	    daughters = 1;
	    if (bp->rho_est == (float)0.0)
	      Error("RhoEst is zero for %s\n", hcellPrint(source));
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
	hmean31 = hmean21 * hmean11;
	hmean41 = hmean21 * hmean21;

	vv2 = dr2 * hmean21;	/* v2 and v renamed to avoid conflict */
	index = vv2 * invdvtable;
	if (index >= MAX_INDEX) Error("Index too large\n");
	vv = rij * hmean11;
	dxx = vv2 - index * dvtable;
	dwdx = (wij[index+1] - wij[index]) * invdvtable;
	wtij = (wij[index] + dwdx * dxx) * hmean31;
	if (wtij < (float)0.0) Error("Negative wtij = %g\n", wtij);
	dgrwdx = (grwij[index+1] - grwij[index]) * invdvtable;
	grwtij = (grwij[index] + dgrwdx * dxx) * hmean41;

	rapm = mass / bp->mass;
	robar1 = (float)2.0 / (rho_est + bp->rho_est);
	grpm = bp->gr_mass * grwtij;
	wpm = bp->gr_mass * wtij;

	poro2 = grpm * (pro2 + bp->pr / (bp->rho_est * bp->rho_est));
	rij1 = (float)1.0 / rij;
	VxVx(runi, = rij1 * r);
	VxVx(f, += poro2 * runi);
	VxVVx(dv, = bp->vel, - v);
	projv = Dotx(dv, runi);

	rhoi += bp->gr_mass * wtij;
	divvi -= bp->gr_mass * grwtij * projv;

	/* artificial viscosity and energy dissipation */
	if (projv < (float)0.0 && alpha != (float)0.0) {
	    vsbar = (float)0.5 * (vsound + bp->vsound);
	    est_divv = projv * vv / (vv2 + epsil);
	    t12 = grpm * est_divv * (beta * est_divv - alpha * vsbar) * robar1;
	    VxVx(f, += t12 * runi);
	    dq += t12 * projv;
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
    VVx(sink->M1, += sink->alfa * f);
    sink->nbrs += nbrs;
}

void
SetSPH(float visc_alpha, float visc_beta, float eos_gamma, int gnobj, 
       void bfunc(), void cfunc())
{
    Nobj = gnobj;
    alpha = visc_alpha;
    Gamma = eos_gamma;
    beta = visc_beta;
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

    /* maximum interaction length and step size */
    v2max = (float)4.0;
    dvtable = v2max / NKERNEL_TABLE;
    i1 = (float)1.0 / dvtable;
    invdvtable = i1;

    /* normalisation constant */
    if (ndim == 3)
      cnormk = (float)M_1_PI;
    else if (ndim == 2)
      Error("Melvyn needs to find this constant\n");
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
    }

    /*  b) v greater than 1 */
    for (i = i1; i <= NKERNEL_TABLE; i++) {
	v2 = i * dvtable;
	v = sqrtf_fast(v2);
	dif2 = (float)2.0 - v;
	sum = 0.25 * dif2 * dif2 * dif2;
	wij[i] = cnormk * sum;
	sum = -0.75 * v2 + 3.0 * v - 3.0;
	grwij[i] = cnormk * sum;
    }
    /* For interpolation of last table entry */
    wij[NKERNEL_TABLE+1] = grwij[NKERNEL_TABLE+1] = (float) 0.0;
}

void
update_final(body *btab, int nobj)
{
    body *p;

    for (p = btab; p < btab+nobj; p++) {
  	p->rho += cnormk * p->gr_mass / (p->h * p->h * p->h);   
	p->udot += p->drho_dt * p->pr * p->gama/ (p->rho * p->rho);
	p->hdot = (float)(-1.0/3.0) * p->h * p->drho_dt / p->rho;
    }
}

void
update_intermediate(body *btab, int nobj, float dt_last, int flag)
{
    body *p;
    float vx,vy,vz,ux,uy,uz,ut,uut,c1,c2,c3;
	 float gamadot;

    for (p = btab; p < btab+nobj; p++) {

     	/* GR */

      	p->enth = (float)1.0+p->u*Gamma;
       	ux = p->mom[0]/p->enth;
       	uy = p->mom[1]/p->enth;
       	uz = p->mom[2]/p->enth;

         c1 = p->gutt;
         c2 = ux*p->guxt+uy*p->guyt+uz*p->guzt;
         c3 = (float)1.0
              +ux*ux*p->guxx
              +uy*uy*p->guyy
              +uz*uz*p->guzz
              +2.0*
              (ux*uy*p->guxy
              +ux*uz*p->guxz
              +uy*uz*p->guyz);
         ut = (-c2+sqrtf_fast(c2*c2-c1*c3))/c1;

         uut
           =ut*p->gutt
           +ux*p->guxt
           +uy*p->guyt
           +uz*p->guzt;
         vx
           =ut*p->guxt
           +ux*p->guxx
           +uy*p->guxy
           +uz*p->guxz;
         vy
           =ut*p->guyt
           +ux*p->guxy
           +uy*p->guyy
           +uz*p->guyz;
         vz
           =ut*p->guzt
           +ux*p->guxz
           +uy*p->guyz
           +uz*p->guzz;

         p->vel[0] = vx/uut;
         p->vel[1] = vy/uut;
         p->vel[2] = vz/uut;

     	   if (flag)
	         p->rho_est = p->rho + p->drho_dt * dt_last;
     	   else
            p->rho_est = p->rho;

         p->mom[3]    = ut;
         p->alfa      = sqrtf_fast(-(float)1.0/p->gutt);
         p->gama_last = p->gama;
       	p->gama      = uut*p->alfa;
     	   p->pr        = p->u*p->rho_est*(Gamma-(float)1.0)/(p->gama);

      	if (p->u <= (float)0.0) {
	    Shout("u is less than zero %f\n", p->u);
	    p->u = 0.0;
	}


         p->vsound    = sqrtf_fast(p->u*Gamma*(Gamma-(float)1.0)
 			                 / (p->u*Gamma+(float)1.0));

     	   gamadot = (p->gama - p->gama_last)/dt_last;
      	p->udot = -p->pr * gamadot / p->rho_est  
      	        +  p->drho_dt * p->pr * p->gama 
	               / (p->rho_est * p->rho_est);

     	/* end GR */
    }
}

void
update_point_mass(body *btab, int nobj, body *p)
{
    body *r;
    float dr2, oneor, oneor2;
    float phii;
    float smooth2 = p->h*p->h;
    Vxd(float r);
    Vxd(float ppos);

    p->phi = (float)0.0;
    VS(p->acc, = (float)0.0);
    VxV(ppos, = p->pos);

    for (r = btab; r < btab+nobj; r++) {
	VxVVx(r, = r->pos, - ppos); /* 3 flops */
	
	dr2 = Dotx(r, r);	/* 5 flops */
	
	dr2 += smooth2;
	
	oneor = recipsqrtf(dr2);	/* 8 flops */
	
	oneor2 = oneor * oneor;	/* 17 flops */
	phii = oneor * p->gr_mass;
                        /* OJO */
	r->phi -= phii;
	VVx(r->acc, -= oneor2 * phii * r);
	phii = oneor * r->gr_mass;
	p->phi -= phii;
	VVx(p->acc, += oneor2 * phii * r);
    }
}
