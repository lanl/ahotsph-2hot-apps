/* Driver for routine fourn */

#include <stdio.h>
#include <stdlib.h>
#define NRANSI
#include "nr.h"
#include "nrutil.h"

#define NDIM 2

#define Index(i,j) (2*((i)*nn[0]+(j)))

int main(int argc, char *argv[])
{
    int isign;
    long idum=(-23);
    unsigned long i,j,k,l,ndum=2,*nn;
    float *data1,*data2,*data3;
    int l2n, ndat2;

    if (argc >= 2)
      l2n = atoi(argv[1]);
    else
      l2n = 4;

    ndat2 = (1 << (NDIM*l2n+1));

    nn=lvector(0,NDIM-1);

    data1=vector(0,ndat2-1);
    data2=vector(0,ndat2-1);
    data3=vector(0,ndat2-1);
    
    for (i=0;i<NDIM;i++)
      nn[i] = (1<<l2n);

    printf ("doing %dx%d transform\n", nn[0],nn[1]);

    for (i=0;i<nn[1];i++) {
	for (j=0;j<nn[0];j++) {
	    l = Index(i,j);
	    data1[l]=data2[l]=data3[l]=2*ran1(&idum)-1;
	    l++;
	    data1[l]=data2[l]=data3[l]=2*ran1(&idum)-1;
	}
    }
    isign = 1;
    fourn(data2-1,nn-1,NDIM,isign);
    fourn(data3-1,nn-1,NDIM,isign);
    isign = -1;
    fourn(data3-1,nn-1,NDIM,isign);

    if (argc > 2) {
	for (i = 0; i < ndat2; i++) {
	    printf("%12f %12f %12f\n", data1[i], data2[i], data3[i]);
	}
    }
    exit(0);
}

