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
#if NDIM > 2
	r2 = *p++;
#endif
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
#if NDIM > 2
	ra2 = *p++;
#endif

	massb = *p++;
	rb0 = *p++;
	rb1 = *p++;
#if NDIM > 2
	rb2 = *p++;
#endif

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
