/* Driver for routine fourn */

#include <stdio.h>
#include <stdlib.h>
#define NRANSI
#include "nr.h"
#include "nrutil.h"

#define NDIM 3

#define Index(i,j,k) (2*(((i)*nn[1]+(j))*nn[2]+(k)))

int main(int argc, char *argv[])
{
    int isign;
    long idum=(-23);
    unsigned long i,j,k,l,ndum=2,*nn;
    float *data1,*data2,*data3;
    int ndat2;

    nn=lvector(0,NDIM-1);

    nn[0] = atoi(argv[1]);
    nn[1] = atoi(argv[2]);
    nn[2] = atoi(argv[3]);

    printf ("doing %dx%dx%d transform\n", nn[0],nn[1],nn[2]);

    ndat2 = nn[0]*nn[1]*nn[2]*2;

    data1=vector(0,ndat2-1);
    data2=vector(0,ndat2-1);
    data3=vector(0,ndat2-1);

    for (i=0;i<nn[0];i++) {
	for (j=0;j<nn[1];j++) {
	    for (k=0;k<nn[2];k++) {
		l = Index(i,j,k);
		data1[l]=data2[l]=data3[l]=2*ran1(&idum)-1;
		l++;
		data1[l]=data2[l]=data3[l]=2*ran1(&idum)-1;
	    }
	}
    }
    isign = 1;
    fourn(data2-1,nn-1,NDIM,isign);
    fourn(data3-1,nn-1,NDIM,isign);
    isign = -1;
    fourn(data3-1,nn-1,NDIM,isign);

        if (argc > 4) {
	for (i = 0; i < ndat2; i++) {
	    printf("%12f %12f %12f\n", 
		   data1[i], data2[i], data3[i]/(1e-20+data1[i]));
	}
    }
    exit(0);
}
