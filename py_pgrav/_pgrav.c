#include <stdio.h>
#include <math.h>

#define MSZ 4

void
_pgravm_sK1(const float *p, float *accp, const int m, const int stride, 
	    const float *f, const int n, const float eps_inv, int *nsmoothed)
{
    double eps_inv2 = eps_inv*eps_inv;
    double eps2 = 1.0/eps_inv2;

    for (int i = 0; i < m; i++) {
	double acc[4] = {};
	for (int j = 0; j < n; j++) {
	    double dx = p[i*stride+1] - f[j*MSZ+1];
	    double dy = p[i*stride+2] - f[j*MSZ+2];
	    double dz = p[i*stride+3] - f[j*MSZ+3];
	    double r2 = dx*dx + dy*dy + dz*dz;
	    if (r2 > eps2) {
		double rinv = -1.0/sqrt(r2);
		double t = rinv*rinv;
		rinv *= f[j*MSZ+0];
		t *= rinv;
		acc[0] += dx*t;
		acc[1] += dy*t;
		acc[2] += dz*t;
		acc[3] += rinv;
	    } else {
		double u2 = r2 * eps_inv2 - 1.0;
		double t = (((45.0/32.0)*u2 + (-3.0/8.0)) * u2 + (1.0/2.0)) * u2 - 1.0;
		t *= f[j*MSZ+0] * eps_inv;
		acc[3] += t;
		t = ((-135.0/16.0) * u2 + (3.0/2.0)) * u2 - 1.0;
		t *= f[j*MSZ+0] * eps_inv * eps_inv2;
		acc[0] += dx*t;
		acc[1] += dy*t;
		acc[2] += dz*t;
	    }
	    /* printf("%2d %2d %12g %12g %12g %12g\n", 
	       i, j, p[i*stride+1], f[j*MSZ+1], f[j*4+0], t); */
	}
	accp[i*stride+0] += acc[0];
	accp[i*stride+1] += acc[1];
	accp[i*stride+2] += acc[2];
	accp[i*stride+3] += acc[3];
    }
}
