/* Driver for routine fourn */

#include <stdio.h>
#include <stdlib.h>
#include "protos.h"
#include "randoms.h"
#include "bigmalloc.h"
#define NRANSI
#include "nr.h"

#define NDIM 3

#define MAXNPROC 64
void ranp_reset(ran_state *rs);

#define Index(i,j,k) (2*(((i)*nn[1]+(j))*nn[2]+(k)))

void
ft(unsigned long *nn, float *data, ran_state *rs, 
   void spectrum(int, int, int, float *, float *))
{
    int isign;
    unsigned long i,j,k,l;
    float *data1;
    int ndat2;
    int ii, ij, ik, is, js;
    float real, imag;
    float x;
    float max = 0.0;

    printf("doing %ldx%ldx%ld transform\n", nn[0],nn[1],nn[2]);
    
    ndat2 = nn[0]*nn[1]*nn[2]*2;
    
    data1=Calloc(ndat2, sizeof(float));
    
    for (i=0;i<nn[0];i++) {
	/* This scheme doesn't work if nn[0] is less than MAXNPROC */
	if (i && (i % (nn[0]/MAXNPROC) == 0)) ranp_reset(rs);
	is = (i < nn[0]/2) ? i : nn[0]-i;
	ii = (i) ? nn[0]-i: i;
	for (j=0;j<nn[1];j++) {
	    js = (j < nn[1]/2) ? j : nn[1]-j;
	    ij = (j) ? nn[1] - j : j;
	    for (k=0;k<=nn[2]/2;k++) {
		ik = (k) ? nn[2] - k : k;
		spectrum(is, js, k, &real, &imag);
		l = Index(i,j,k);
		data1[l] += real;
		data1[l+1] += imag;
		l = Index(ii, ij, ik);
		data1[l] += real;
		data1[l+1] -= imag;
	    }
	}
    }
    isign = -1;
    
    fourn(data1-1,nn-1,NDIM,isign);

    for (i=0;i<nn[0];i++) {
	for (j=0;j<nn[1];j++) {
	    for (k=0;k<nn[2];k++) {
		l = Index(i, j, k);
		data[l/2] = data1[l];
		x = data1[l+1];
		if (x > max) max = x;
	    }
	}
    }
    Free(data1);
    printf("max imaginary value is %f\n", max);
}
