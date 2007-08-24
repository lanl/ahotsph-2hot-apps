#define NO_MSGS
#define NOTIMERS  /* Timers are a major performance hit on the delta */
#include "physics.h"
#include "vop.h"
#include "tensop.h"
#include "fastflpt.h"
#include "Msgs.h"
#include "timers.h"
#include "stk.h"


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
update_point_mass(body *btab, int nobj, 
		  body *p, float smooth2, float newt)
{
    body *r;
    float dr2, oneor, oneor2;
    float phii;
    Vxd(float r);
    Vxd(float ppos);

    VxV(ppos, = p->pos);

    for (r = btab; r < btab+nobj; r++) {
	VxVVx(r, = r->pos, - ppos); /* 3 flops */
	
	dr2 = Dotx(r, r);	/* 5 flops */
	
	if (dr2 != (float)0.0) {
	
	  dr2 += smooth2;
	
	  oneor = recipsqrtf(dr2);	/* 8 flops */
	
	  oneor2 = oneor * oneor;	/* 17 flops */
	  phii = newt * oneor * r->mass;
	  p->phi -= phii;
	  VVx(p->acc, += oneor2 * phii * r);
	}
    }
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

#if defined(__T3D__) || defined (_IBMR2)
#define USE_CHEB_RSQRT
#endif

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

