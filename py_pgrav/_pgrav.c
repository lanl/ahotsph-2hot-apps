#include <stdio.h>
#include <math.h>

#define MSZ 4

void
_pgravm(const float *p, float *accp, const int m, const int stride, 
	const float *f, const int n)
{
    for (int i = 0; i < m; i++) {
	double acc[4] = {};
	for (int j = 0; j < n; j++) {
	    double dx = p[i*stride+1] - f[j*MSZ+1];
	    double dy = p[i*stride+2] - f[j*MSZ+2];
	    double dz = p[i*stride+3] - f[j*MSZ+3];
	    double r2 = dx*dx + dy*dy + dz*dz;
	    double rinv = -1.0/sqrt(r2);
	    double t = rinv*rinv;
	    rinv *= f[j*4+0];
	    t *= rinv;
	    acc[0] += dx*t;
	    acc[1] += dy*t;
	    acc[2] += dz*t;
	    acc[3] += rinv;
	    printf("%2d %2d %12g %12g %12g %12g\n", 
		   i, j, p[i*stride+1], f[j*MSZ+1], f[j*4+0], t);
	}
	accp[i*stride+0] += acc[0];
	accp[i*stride+1] += acc[1];
	accp[i*stride+2] += acc[2];
	accp[i*stride+3] += acc[3];
    }
}
