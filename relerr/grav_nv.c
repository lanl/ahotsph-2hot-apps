#define NO_MSGS
#define NOTIMERS  /* Timers are a major performance hit on the delta */
#include "physics_n.h"
#include "vop.h"
#include "tensop.h"
#include "fastflpt.h"
#include "Msgs.h"
#include "timers.h"
#include "stk.h"


#define SINK_BIAS ((float)2.1)

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
static int add_offset;
static float offset[NDIM];
int last_icnt;

void
SetGravOffset(float *off)
{
    VV(offset, = off);
    add_offset = 1;
}

void
UnSetGravOffset(void)
{
    VS(offset, = 0.0);
    add_offset = 0;
}


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
	/* Make sure these are initialized to zero externally */
	bp->phi += GNewt*from->M0;
	VV(bp->acc, += -GNewt*from->M1);
	bp->nterms += from->nterms;
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
    if (add_offset) {
	VV(to->pos, += offset);
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
	    (B2 + (float)6.0*mass*bmaxsrc*bmaxsink +
	     (float)3.0*mass*bmaxsink*bmaxsink - rinv*B3))
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
	result[i] = (SINK_BIAS*bmaxsink >= bmaxsrc) ? MAC_SPLIT_SINK 
	  : MAC_SPLIT_SRC;
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
	result[i] = (SINK_BIAS*bmaxsink >= bmaxsrc) ? MAC_SPLIT_SINK 
	  : MAC_SPLIT_SRC;
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
	    (B2 + (float)6.0*mass*bmaxsrc*bmaxsink + 
	     (float)3.0*mass*bmaxsink*bmaxsink - rinv*B3))
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
	result[i] = (SINK_BIAS*bmaxsink >= bmaxsrc) ? MAC_SPLIT_SINK 
	  : MAC_SPLIT_SRC;
    }
    sink->M0 += phi;
    VVx(sink->M1, += a);
    sink->interactions += interactions;
    sink->nterms += (float)nterms;
    StopTimer(&GravTm);
}

/* The fast N log N stuff is below here */

#define UNROLL 2		/* unroll by 2 by default */
#define INTERACTF do_grav_u2

#if defined(__INTEL_SSD__) && !defined(NO_ASM)
#undef UNROLL
#undef INTERACTF
#define UNROLL 3
#define INTERACTF do_grav_fast
#endif 

#if defined(__T3D__) || defined (_IBMR2)
#define USE_CHEB_RSQRT
#endif

#ifdef USE_CHEB_RSQRT
#undef INTERACTF
#define INTERACTF do_cheb_u2
#endif

#ifdef __ncube__
#undef UNROLL
#undef INTERACTF
#endif

#ifdef __CM5VU__
#undef UNROLL
#undef INTERACTF
#define UNROLL 32
#define INTERACTF do_grav_cm5
#endif

#if 0
/* If you want Mflops, this is the number */
/* It roughly doubles nterms over using 3.0-4.0 */
#define RcritFac ((float)2.0)	/* should be >= 2.0 */
#else
/* If you want bodies/sec, this is close to the right number */
#define RcritFac ((float)4.0)	/* should be >= 2.0 */
#endif

/* This should be dynamically extensible */
#define IVECSZ 40960
struct {
    float mass;
    float pos[NDIM];
} Ivec[IVECSZ];

void DLmacv(Sink *sink, const hcell **source_vec, int *result, int n);
void Nlognmacv(Sink *sink, const hcell **source_vec, int *result, int);
void INTERACTF(const float *p, const float *end, const float *pos0, 
	       float *mass0, float *acc0, float *phi0, const float *eps2p, 
	       int *ncut);

void 
InheritSinkNlogN(const Sink *from, Sink *to, hcell *pp)
{
    if( to == NULL ){
	body *bp = pp->ptr;
	/* must init mtot or else you get quiet exceptions in asm code */
	float mtot = (float)0.0;
	int ijunk = 0, nn;
	float acc[NDIM];
	float phi;

	VS(acc, = (float)0.0);
	phi = (float)0.0;
	/* putting a getrusage based timer here can slow things down a lot */
	StartTimer(&GravTm);
#ifdef UNROLL
	nn = from->icnt - (from->icnt % UNROLL);
	/* Msg_do("icnt is %d\n", from->icnt); */
	/* Use the interface to the fast assembly code */
#ifdef __INTEL_SSD__
	if ((int)&Ivec[0] & 07 || (int)&Ivec[1] & 07)
	  Error("Ivec not aligned for asm code\n");
#endif
	INTERACTF((float *)&Ivec[0], (float *)&Ivec[nn], 
		  from->pos, &mtot, acc, &phi, &eps2, &ijunk);
	if (from->icnt % UNROLL)
	  do_grav((float *)&Ivec[nn], (float *)&Ivec[from->icnt], 
		  from->pos, &mtot, acc, &phi, &eps2, &ijunk);
#else
	do_grav((float *)&Ivec[0], (float *)&Ivec[from->icnt], 
		from->pos, &mtot, acc, &phi, &eps2, &ijunk);
#endif
	StopTimer(&GravTm);

#if 0	/* This is only needed for testing */
	if (from->m+mtot+bp->mass < .999 || from->m+mtot+bp->mass > 1.001) {
	    Error("bad m (%f)\n", from->m+mtot+bp->mass);
	}
#endif
	/* Make sure these are initialized to zero externally */
	bp->phi += GNewt*from->M0;
	bp->phi += GNewt*phi;
	VV(bp->acc, += -GNewt*from->M1);
	VV(bp->acc, += GNewt*acc);
	/* guess that icnt interactions are 3x faster than nterms */
	bp->nterms += from->nterms + from->icnt;
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
    if (add_offset) {
	VV(to->pos, += offset);
    }

    if (from) {
	to->interactions = from->interactions;
	to->nterms = from->nterms;
	to->m = from->m;
	to->M0 = from->M0;
	VV(to->M1, = from->M1);
	to->icnt = from->icnt;
	if (to->icnt < last_icnt) last_icnt = (to->icnt & ~31);
	if (to->icnt >= IVECSZ) Error("ivec overflow\n");

	Msgf(("inherit %f to %f, %f to %f\n", from->pos[0], to->pos[0],
	      from->M1[0], to->M1[0]));
    } else {
	to->interactions = 0;
	to->nterms = 0;
	to->m = (float)0.0;
	to->M0 = (float)0.0;
	VS(to->M1, = (float)0.0);
	to->icnt = 0;
	last_icnt = 0;		/* use to optimize cm5 VU copy */
    }
}

void
Nlogngate(Sink *sink, const hcell **source_vec, int *result, int n)
{
    if (sink->isbody) {
	Nlognmacv(sink, source_vec, result, n);
    } else {
	DLmacv(sink, source_vec, result, n);
    }
}

void
DLmacv(Sink *sink, const hcell **source_vec, int *result, int n)
{
    VxdV(float pos_sink, = sink->pos);
    float bmax = sink->bmax;
    int icnt = sink->icnt;
    int interactions = 0;
    float dr2;
    Vxd(float r);
    Vxd(float dx);
    float rcrit_bmax;
    int i;

    for (i = 0; i < n; i++) {
	if (Sub_Flags(source_vec[i])) {
	    const cell *cp = source_vec[i]->ptr;
	    VxV(r, = cp->pos);	
	    rcrit_bmax = cp->rcrit + bmax;

	    VxVxVx(dx, = r, - pos_sink);
	    dr2 = Dotx(dx, dx);

	    if (dr2 > rcrit_bmax*rcrit_bmax) {
		Ivec[icnt].mass = cp->mass;
		VV(Ivec[icnt].pos, = cp->pos);
		icnt++;
		interactions += cp->daughters;
		result[i] = MAC_ACCEPT;
		IncrCounter(&CCInt);
	    } else if (RcritFac * bmax > rcrit_bmax || dr2 == (float)0.0) { 
		result[i] = MAC_SPLIT_SINK;
		IncrCounter(&CCIntRej);
	    } else {
		result[i] = MAC_SPLIT_SRC;
		IncrCounter(&CCIntRej);
	    }
	} else {
	    const body *bp = source_vec[i]->ptr;
	    Ivec[icnt].mass = bp->mass;
	    VV(Ivec[icnt].pos, = bp->pos);
	    icnt++;
	    interactions += 1;
	    result[i] = MAC_ACCEPT;
	    IncrCounter(&CBInt);
	}
    }
    sink->interactions += interactions;
    sink->icnt = icnt;
}


void
Nlognmacv(Sink *sink, const hcell **source_vec, int *result, int n)
{
    VxdV(float pos_sink, = sink->pos);
    float phi = 0.F;
    float m = (float)0.0;
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

	m += mass;
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
    sink->m += m;
    StopTimer(&GravTm);
}

void 
do_grav(const float *p, const float *end, const float *pos0, float *mass0, 
	float *acc0, float *phi0, const float *eps2p, int *ncut)
{
    float dr2;
    Vxd(float r);
    float phii, mor3, mass;
    VxdV(const float ppos, = pos0);
    Vxd(float a);
    float phi = *phi0;
    float total_mass = *mass0;
    const float eps2 = *eps2p;

    VxV(a, = acc0);

    while (p < end) {
	mass = *p++;
	r0 = *p++;
	r1 = *p++;
	r2 = *p++;
	VxVx(r, -= ppos);	/* 3 flops */

	dr2 = Dotx(r, r);	/* 5 flops */
	dr2 += eps2;

	phii = recipsqrtf(dr2);	/* 8 flops */
	
	mor3 = phii * phii;	/* 5 flops */
	phii *= mass;
	total_mass += mass;
	mor3 *= phii;
	phi -= phii;

	VxVx(a, += mor3 * r);	/* 6 flops */
    }
    VVx(acc0, = a);
    *mass0 = total_mass;
    *phi0 = phi;
}

void 
do_grav_u2(const float *p, const float *end, const float *pos0, float *mass0,
	   float *acc0,	float *phi0, const float *eps2p, int *ncut)
{
    float dr2a, dr2b;
    Vxd(float ra);
    Vxd(float rb);
    float phiia, phiib;
    float mor3a, mor3b;
    float massa, massb;
    VxdV(const float ppos, = pos0);
    Vxd(float a);
    float phi = *phi0;
    float total_mass = *mass0;
    const float eps2 = *eps2p;
    int n = end-p;
    int i;

    VxV(a, = acc0);

    /* This dummy loop index improves R10000 IRIX64 performance by 50% */
    for (i = 0; i < n/8; i++) {
	massa = *p++;
	ra0 = *p++;
	ra1 = *p++;
	ra2 = *p++;

	massb = *p++;
	rb0 = *p++;
	rb1 = *p++;
	rb2 = *p++;

	VxVx(ra, -= ppos);
	VxVx(rb, -= ppos);

	dr2a = Dotx(ra, ra);
	dr2b = Dotx(rb, rb);

	dr2a += eps2;
	dr2b += eps2;

	phiia = recipsqrtf(dr2a);
	phiib = recipsqrtf(dr2b);
	
	mor3a = phiia * phiia;
	mor3b = phiib * phiib;
	phiia *= massa;
	phiib *= massb;
	total_mass += massa + massb;
	mor3a *= phiia;
	mor3b *= phiib;
	phi -= phiia;
	phi -= phiib;

	VxVx(a, += mor3a * ra);
	VxVx(a, += mor3b * rb);
    }
    VVx(acc0, = a);
    *mass0 = total_mass;
    *phi0 = phi;
}

#ifdef USE_CHEB_RSQRT
#include "karp.h"

/* This does 38 actual flops per interaction */
/* The approximate rsqrt uses 13 multiplies and 6 adds */

void 
do_cheb_u2(const float *p, const float *end, const float *pos0, float *mass0,
	   float *acc0,	float *phi0, const float *eps2p, int *ncut)
{
    float dr2a, dr2b;
    Vxd(float ra);
    Vxd(float rb);
    float mor3a, mor3b;
    float massa, massb;
    VxdV(const float ppos, = pos0);
    Vxd(float a);
    float phi = *phi0;
    float total_mass = *mass0;
    const float eps2 = *eps2p;
    float xa, xb;
    unsigned int ita, itb;
    const int mbits = 23;

    VxV(a, = acc0);

    while (p < end) {
	massa = *p++;
	ra0 = *p++;
	ra1 = *p++;
	ra2 = *p++;

	massb = *p++;
	rb0 = *p++;
	rb1 = *p++;
	rb2 = *p++;

	VxVx(ra, -= ppos);
	VxVx(rb, -= ppos);

	dr2a = Dotx(ra, ra);
	dr2b = Dotx(rb, rb);
	dr2a += eps2;
	dr2b += eps2;

	ita = (*(FLOAT_PUN *)&dr2a) >> mbits;
	itb = (*(FLOAT_PUN *)&dr2b) >> mbits; 
	xa = dr2a*u[ita];
	xb = dr2b*u[itb];
	mor3a = g0 + xa*(g1 + xa*(g2 + xa*(g3 + xa*(g4 + xa*g5))));
	mor3b = g0 + xb*(g1 + xb*(g2 + xb*(g3 + xb*(g4 + xb*g5))));
	mor3a *= (float)1.5 - (float)0.5*xa*xa*xa*mor3a*mor3a;
	mor3b *= (float)1.5 - (float)0.5*xb*xb*xb*mor3b*mor3b;
	mor3a *= t[ita];
	mor3b *= t[itb];

	mor3a *= massa;
	mor3b *= massb;
	total_mass += massa + massb;
	phi -= dr2a*mor3a;
	phi -= dr2b*mor3b;

	VxVx(a, += mor3a * ra);
	VxVx(a, += mor3b * rb);
    }
    VVx(acc0, = a);
    *mass0 = total_mass;
    *phi0 = phi;
}
#endif

