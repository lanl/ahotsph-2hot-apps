/* #define NO_MSGS */
#include "physics_sph.h"
#include "vop.h"
#include "tensop.h"
#include "fastflpt.h"
#include "Msgs.h"
#include "timers.h"
#include "stk.h"

Counter_t CCInt, CBInt, BCInt, BBInt;
Counter_t CCIntRej;
Counter_t TranslateCnt;

Timer_t Imbal;
Timer_t GravTm;

static float acc_tolerance;
static float frac_tolerance;
static float eps2;
static float GNewt;
static int Nobj;

void
SetTol(float tol, float frac_tol, float newton_const, float eps, int gnobj)
{
    acc_tolerance = tol;
    cofm_setup(tol);
    frac_tolerance = frac_tol/newton_const;
    eps2 = eps*eps;
    GNewt = newton_const;
    Nobj = gnobj;
}

void 
InheritSink(const Sink *from, Sink *to, hcell *pp)
{
    float xtau[NDIM];
    float dot;

    if( to == NULL ){
	body *bp = pp->ptr;
	bp->phi = GNewt*from->M0;
	VV(bp->acc, = -GNewt*from->M1);
	bp->nterms = from->nterms;
	if (from->interactions != Nobj)
	    Error("Ninteract is %d, should be %d\n", from->interactions, Nobj);
	return;
    }

    IncrCounter(&TranslateCnt);
    if (Sub_Flags(pp)) {
	cell *cp = pp->ptr;

	VV(to->pos, = cp->pos);
	to->bmax = cp->bmax;
	to->acc_last_max = cp->acc_last_max;
	to->daughters = cp->daughters;
	to->isbody = 0;
    } else {
	body *bp = pp->ptr;
	VV(to->pos, = bp->pos);
	to->bmax = (float)0.0;
	to->acc_last_max = bp->acc_last;
	to->daughters = 1.F;
	to->isbody = 1;
    }

    if (from) {
	to->interactions = from->interactions;
	to->nterms = from->nterms * ((float)to->daughters/(float)from->daughters);
	VVV(xtau, = from->pos, - to->pos); /* 8 flops */
	dot = Dot(xtau, from->M1);

	to->M0 = from->M0 - dot;
	VV(to->M1, = from->M1);

	if (!to->isbody) {
	    to->M2 = from->M2;
	}

	to->M1[0] += xtau[0]*from->M2.xx; /* 18 flops */
	to->M1[0] += xtau[1]*from->M2.xy;
	to->M1[0] += xtau[2]*from->M2.xz;
	to->M1[1] += xtau[0]*from->M2.xy;
	to->M1[1] += xtau[1]*from->M2.yy;
	to->M1[1] += xtau[2]*from->M2.yz;
	to->M1[2] += xtau[0]*from->M2.xz;
	to->M1[2] += xtau[1]*from->M2.yz;
	to->M1[2] += xtau[2]*from->M2.zz;
	Msgf(("inherit %f to %f, %f to %f\n", from->pos[0], to->pos[0],
	      from->M1[0], to->M1[0]));
    } else {
	to->interactions = 0;
	to->nterms = 0;
	to->M0 = (float)0.0;
	VS(to->M1, = (float)0.0);
	TS(to->M2, = (float)0.0);
    }
}


/* 74 flops for cell-cell interactions */
/* 51 flops for body-cell interactions */
/* 32 flops for failed MAC */

void
Unifiedmacv(Sink *sink, const hcell **source_vec, int *result, int n)
{
    const float bmaxsink = sink->bmax;
    VxdV(const float pos_sink, = sink->pos);
    float phi = 0.F;
    Vxd(float a);
    int interactions = 0;
    int nterms = 0;
    float dr2;
    Vxd(float r);
    float mass;
    float mor3;
    float er2, rinv, rinv2;
    float bmaxsrc;
    float B2, B3;
    int daughters;
    int i;

    VxS(a, = 0.F);
    StartTimer(&GravTm);
    for (i = 0; i < n; i++) {
	const hcell *source = source_vec[i];

	/* Be aware of memory access patterns */
	if (Sub_Flags(source)) {
	    /* non-terminal source */
	    const cell *cp = source->ptr;
	    mass = cp->mass;	/* Access in same order as cell struct */
	    VxV(r, = cp->pos);	
	    B2 = cp->B2;
	    B3 = cp->B3;
	    bmaxsrc = cp->bmax;
	    daughters = cp->daughters;
	} else {
	    const body *bp = source->ptr;
	    mass = bp->mass;
	    VxV(r, = bp->pos);
	    bmaxsrc = B2 = B3 = (float)0.0;
	    daughters = 1;
	}
	
	VxVx(r, -= pos_sink);	/* 8 flops */
	dr2 = Dotx(r, r);

	if (dr2 == (float)0.0) {
	    if (bmaxsink == (float)0.0)
	      goto accept;
	    else
	      goto failed;
	}

	rinv = recipsqrt8bit(dr2);		/* 3 or 10 flops */
	er2 = (bmaxsink + bmaxsrc) * rinv;

	if (er2 >= (float)1.0)  		/* 14 flops */
	  goto failed;
	er2 = (float)1.0 - er2;
	er2 *= er2;
	rinv2 = rinv*rinv;
	if (er2 * acc_tolerance < rinv2*rinv2 * 
	    (B2 + (float)3.0*mass*bmaxsink*bmaxsink - rinv*B3))
	  goto failed;

	dr2 += eps2;			/* 10 flops */
	rinv = recipsqrtf(dr2);
	rinv2 = rinv*rinv;

	phi -= mass * rinv;		/* 9 flops */
	mor3 = mass * rinv * rinv2;

	VxVx(a, -= mor3 * r);
	Msgf(("%f %f %f\n", pos_sink0, r0, mor3 * r0));
	nterms++;		/* one 'term' */

	if (!sink->isbody) {		/* 23 flops */
	    float mor5 = (float)3.0 * mor3 * rinv2;
	    moment *qpole0 = &sink->M2;
	    qpole0->xx += r0*r0*mor5 - mor3;
	    qpole0->yy += r1*r1*mor5 - mor3;
	    qpole0->zz += r2*r2*mor5 - mor3;
	    qpole0->xy += r0*r1*mor5;
	    qpole0->xz += r0*r2*mor5;
	    qpole0->yz += r1*r2*mor5;
	    nterms++;		/* add another 'term' for your trouble */
	    if (daughters > 1)
	      IncrCounter(&CCInt);
	    else 
	      IncrCounter(&CBInt);
	} else {
	    if (daughters > 1)
	      IncrCounter(&BCInt);
	    else 
	      IncrCounter(&BBInt);
	}
      accept:
	interactions += daughters;
	result[i] = MAC_ACCEPT;
	Msgf(("a %s %d %.2f %.2f\n", (sink->isbody) ? "b" : "c", daughters,
	      bmaxsink, bmaxsrc));
	continue;
      failed:
	IncrCounter(&CCIntRej);
	Msgf(("r %s %d %.2f %.2f\n", (sink->isbody) ? "b" : "c", daughters,
	      bmaxsink, bmaxsrc));
	result[i] = (bmaxsink > bmaxsrc) ? MAC_SPLIT_SINK : MAC_SPLIT_SRC;
    }
    sink->M0 += phi;
    VVx(sink->M1, += a);
    sink->interactions += interactions;
    sink->nterms += (float)nterms;
    StopTimer(&GravTm);
}

void
Lowestmacv(Sink *sink, const hcell **source_vec, int *result, int n)
{
    const float bmaxsink = sink->bmax;
    VxdV(const float pos_sink, = sink->pos);
    float phi = 0.F;
    Vxd(float a);
    int interactions = 0;
    int nterms = 0;
    float dr2;
    Vxd(float r);
    float mass;
    float mor3;
    float er2, rinv, rinv2;
    float bmaxsrc;
    float B2, B3;
    int daughters;
    int i;

    VxS(a, = 0.F);
    StartTimer(&GravTm);
    for (i = 0; i < n; i++) {
	const hcell *source = source_vec[i];

	/* Be aware of memory access patterns */
	if (Sub_Flags(source)) {
	    /* non-terminal source */
	    const cell *cp = source->ptr;
	    mass = cp->mass;	/* Access in same order as cell struct */
	    VxV(r, = cp->pos);	
	    B2 = cp->B2;
	    B3 = cp->B3;
	    bmaxsrc = cp->bmax;
	    daughters = cp->daughters;
	} else {
	    const body *bp = source->ptr;
	    mass = bp->mass;
	    VxV(r, = bp->pos);
	    bmaxsrc = B2 = B3 = (float)0.0;
	    daughters = 1;
	}
	
	VxVx(r, -= pos_sink);	/* 8 flops */
	dr2 = Dotx(r, r);

	if (dr2 == (float)0.0) {
	    if (bmaxsink == (float)0.0)
	      goto accept;
	    else
	      goto failed;
	}

	rinv = recipsqrt8bit(dr2);		/* 3 or 10 flops */
	er2 = (bmaxsink + bmaxsrc) * rinv;

	if (er2 >= (float)1.0)  		/* 14 flops */
	  goto failed;

	dr2 += eps2;			/* 10 flops */
	rinv = recipsqrtf(dr2);
	rinv2 = rinv*rinv;

	phi -= mass * rinv;		/* 9 flops */
	mor3 = mass * rinv * rinv2;

	VxVx(a, -= mor3 * r);
	Msgf(("%f %f %f\n", pos_sink0, r0, mor3 * r0));
	nterms++;		/* one 'term' */

	if (!sink->isbody) {		/* 23 flops */
	    float mor5 = (float)3.0 * mor3 * rinv2;
	    moment *qpole0 = &sink->M2;
	    qpole0->xx += r0*r0*mor5 - mor3;
	    qpole0->yy += r1*r1*mor5 - mor3;
	    qpole0->zz += r2*r2*mor5 - mor3;
	    qpole0->xy += r0*r1*mor5;
	    qpole0->xz += r0*r2*mor5;
	    qpole0->yz += r1*r2*mor5;
	    nterms++;		/* add another 'term' for your trouble */
	    if (daughters > 1)
	      IncrCounter(&CCInt);
	    else 
	      IncrCounter(&CBInt);
	} else {
	    if (daughters > 1)
	      IncrCounter(&BCInt);
	    else 
	      IncrCounter(&BBInt);
	}
      accept:
	interactions += daughters;
	result[i] = MAC_ACCEPT;
	Msgf(("a %s %d %.2f %.2f\n", (sink->isbody) ? "b" : "c", daughters,
	      bmaxsink, bmaxsrc));
	continue;
      failed:
	IncrCounter(&CCIntRej);
	Msgf(("r %s %d %.2f %.2f\n", (sink->isbody) ? "b" : "c", daughters,
	      bmaxsink, bmaxsrc));
	result[i] = (bmaxsink > bmaxsrc) ? MAC_SPLIT_SINK : MAC_SPLIT_SRC;
    }
    sink->M0 += phi;
    VVx(sink->M1, += a);
    sink->interactions += interactions;
    sink->nterms += (float)nterms;
    StopTimer(&GravTm);
}

void
Fracmacv(Sink *sink, const hcell **source_vec, int *result, int n)
{
    const float bmaxsink = sink->bmax;
    float acc_last_tol;
    VxdV(const float pos_sink, = sink->pos);
    float phi = 0.F;
    Vxd(float a);
    int interactions = 0;
    int nterms = 0;
    float dr2;
    Vxd(float r);
    float mass;
    float mor3;
    float er2, rinv, rinv2;
    float bmaxsrc;
    float B2, B3;
    int daughters;
    int i;

    VxS(a, = 0.F);
    StartTimer(&GravTm);
    acc_last_tol = sink->acc_last_max*frac_tolerance;
    if (acc_last_tol < acc_tolerance)
      acc_last_tol = acc_tolerance;
    for (i = 0; i < n; i++) {
	const hcell *source = source_vec[i];

	/* Be aware of memory access patterns */
	if (Sub_Flags(source)) {
	    /* non-terminal source */
	    const cell *cp = source->ptr;
	    mass = cp->mass;	/* Access in same order as cell struct */
	    VxV(r, = cp->pos);	
	    B2 = cp->B2;
	    B3 = cp->B3;
	    bmaxsrc = cp->bmax;
	    daughters = cp->daughters;
	} else {
	    const body *bp = source->ptr;
	    mass = bp->mass;
	    VxV(r, = bp->pos);
	    bmaxsrc = B2 = B3 = (float)0.0;
	    daughters = 1;
	}
	
	VxVx(r, -= pos_sink);	/* 8 flops */
	dr2 = Dotx(r, r);

	if (dr2 == (float)0.0) {
	    if (bmaxsink == (float)0.0)
	      goto accept;
	    else
	      goto failed;
	}

	rinv = recipsqrt8bit(dr2);		/* 3 or 10 flops */
	er2 = (bmaxsink + bmaxsrc) * rinv;

	if (er2 >= (float)1.0)  		/* 14 flops */
	  goto failed;
	er2 = (float)1.0 - er2;
	er2 *= er2;
	rinv2 = rinv*rinv;
	if (er2 * acc_last_tol < rinv2*rinv2 * 
	    (B2 + (float)3.0*mass*bmaxsink*bmaxsink - rinv*B3))
	  goto failed;

	dr2 += eps2;			/* 10 flops */
	rinv = recipsqrtf(dr2);
	rinv2 = rinv*rinv;

	phi -= mass * rinv;		/* 9 flops */
	mor3 = mass * rinv * rinv2;

	VxVx(a, -= mor3 * r);
	Msgf(("%f %f %f\n", pos_sink0, r0, mor3 * r0));
	nterms++;		/* one 'term' */

	if (!sink->isbody) {		/* 23 flops */
	    float mor5 = (float)3.0 * mor3 * rinv2;
	    moment *qpole0 = &sink->M2;
	    qpole0->xx += r0*r0*mor5 - mor3;
	    qpole0->yy += r1*r1*mor5 - mor3;
	    qpole0->zz += r2*r2*mor5 - mor3;
	    qpole0->xy += r0*r1*mor5;
	    qpole0->xz += r0*r2*mor5;
	    qpole0->yz += r1*r2*mor5;
	    nterms++;		/* add another 'term' for your trouble */
	    if (daughters > 1)
	      IncrCounter(&CCInt);
	    else 
	      IncrCounter(&CBInt);
	} else {
	    if (daughters > 1)
	      IncrCounter(&BCInt);
	    else 
	      IncrCounter(&BBInt);
	}
      accept:
	interactions += daughters;
	result[i] = MAC_ACCEPT;
	Msgf(("a %s %d %.2f %.2f\n", (sink->isbody) ? "b" : "c", daughters,
	      bmaxsink, bmaxsrc));
	continue;
      failed:
	IncrCounter(&CCIntRej);
	Msgf(("r %s %d %.2f %.2f\n", (sink->isbody) ? "b" : "c", daughters,
	      bmaxsink, bmaxsrc));
	result[i] = (bmaxsink > bmaxsrc) ? MAC_SPLIT_SINK : MAC_SPLIT_SRC;
    }
    sink->M0 += phi;
    VVx(sink->M1, += a);
    sink->interactions += interactions;
    sink->nterms += (float)nterms;
    StopTimer(&GravTm);
}


void
Nlognmacv(Sink *sink, const hcell **source_vec, int *result, int n)
{
    VxdV(float pos_sink, = sink->pos);
    float phi = 0.F;
    Vxd(float a);
    int interactions = 0;
    int nterms = 0;
    float dr2;
    Vxd(float r);
    float mass;
    float mor3;
    float rinv, rinv2;
    float rcrit;
    int daughters;
    int i;

    VxS(a, = 0.F);
    if (!sink->isbody) {
	for (i = 0; i < n; i++) {
	    result[i] = MAC_SPLIT_SINK;
	}
	return;
    }

    StartTimer(&GravTm);
    for (i = 0; i < n; i++) {
	const hcell *source = source_vec[i];
	if (Sub_Flags(source)) {
	    const cell *cp = source->ptr;
	    mass = cp->mass;	/* Access in same order as cell struct */
	    VxV(r, = cp->pos);	
	    rcrit = cp->rcrit;
	    daughters = cp->daughters;
	} else {
	    const body *bp = source->ptr;
	    mass = bp->mass;
	    VxV(r, = bp->pos);
	    rcrit = (float)0.0;
	    daughters = 1;
	}

	VxVx(r, -= pos_sink);
	dr2 = Dotx(r, r);

	if (dr2 == (float)0.0)
	  goto accept;

	if (dr2 < rcrit*rcrit)
	  goto failed;

	dr2 += eps2;
	rinv = recipsqrtf(dr2);
	rinv2 = rinv*rinv;

	phi -= mass * rinv;
	mor3 = mass * rinv * rinv2;

	VxVx(a, -= mor3 * r);
	nterms++;

      accept:
	if (daughters > 1)
	  IncrCounter(&BCInt);
	else 
	  IncrCounter(&BBInt);
	interactions += daughters;
	result[i] = MAC_ACCEPT;
	continue;
      failed:
	IncrCounter(&CCIntRej);
	result[i] = MAC_SPLIT_SRC;
    }
    sink->M0 += phi;
    VVx(sink->M1, += a);
    sink->interactions += interactions;
    sink->nterms += (float)nterms;
    StopTimer(&GravTm);
}
