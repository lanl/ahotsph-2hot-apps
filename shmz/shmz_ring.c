#include <math.h>
#include "physics.h"
#include "vop.h"
#include "fastflpt.h"

static float k;

void
set_body(void *o, void *p)
{
    memcpy(o, p, TBODYSZ);
}

void
set_k(float lambda)
{
    k = 2.0*M_PI/lambda;
}

void 
do_shmz(void *p0, void *list, int bsize, int n)
{
    Vxd(float r);
    float stren;
    double dr2;
    double kx, kx_inv;
    double phi_r, phi_i;
    double s, c;
    VxdV(const float ppos, = ((const body *)p0)->pos);
    float *p = list;
    float *end = ((float *)list) + n * bsize/sizeof(float);

    phi_r = phi_i = 0.0;

    while (p < end) {
	stren = *p++;
	r0 = *p++;
	r1 = *p++;
	r2 = *p++;
	VxVx(r, -= ppos);
	dr2 = Dotx(r, r)*k*k;
	if (dr2 == 0.0) continue;
	kx_inv = 1.0/sqrt(dr2);
	kx = kx_inv * dr2;
#ifdef HAS_SINCOS
	sincos(kx, &s, &c);
#else
	s = sin(kx); c = cos(kx); 
#endif
	phi_r += s*kx_inv*stren;
	phi_i += c*kx_inv*stren;
    }
    ((body *)p0)->phi_r += phi_r;
    ((body *)p0)->phi_i += phi_i;
}
	  
