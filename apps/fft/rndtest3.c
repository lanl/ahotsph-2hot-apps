#include <stdio.h>
#include <stdlib.h>
#include "protos.h"

float ran1(long *idum);
void fftndrf(int ndim, int npts[], float *data);
void fftndrb(int ndim, int npts[], float *data);
void fftndcf(int ndim, int npts[], float *data);
void fftndcb(int ndim, int npts[], float *data);

#define NDIM 3

#define Index(i,j,k) ((((i)*nn[1]+(j))*nn[2]+(k)))

int main(int argc, char *argv[])
{
    int isign;
    long idum=(-23);
    unsigned long i,j,k,l;
    int *nn;
    float *data1;
    int ndat2;
    FILE *fp;

    nn=malloc(3*sizeof(int));

    nn[0] = atoi(argv[1]);
    nn[1] = atoi(argv[2]);
    nn[2] = atoi(argv[3]);

    printf ("doing %dx%dx%d real transform\n", nn[0],nn[1],nn[2]);

    ndat2 = nn[0]*nn[1]*nn[2];

    data1=malloc(ndat2*sizeof(float));
    memset(data1, ndat2*sizeof(float), 0);

    for (i=0;i<nn[0];i++) {
	float real, imag;
	for (j=0;j<nn[1];j++) {
	    for (k=0;k<=nn[2];k+=2) {
		real = 2*ran1(&idum)-1;
		imag = 2*ran1(&idum)-1;
		l = Index(i,j,k);
		if (k == 0) {
		    int ii, ij, ll;
		    ii = i;
		    ij = j;
		    if (ii) ii = nn[0]-ii;
		    if (ij) ij = nn[1]-ij;
		    ll = Index(ii, ij, 0);
		    data1[l] += real;
		    data1[l+1] += imag;
		    data1[ll] += real;
		    data1[ll+1] -= imag;
		} else if (k == nn[2]) {
		    int ii, ij, ll;
		    l = Index(i,j,0);
		    ii = i;
		    ij = j;
		    if (ii) ii = nn[0]-ii;
		    if (ij) ij = nn[1]-ij;
		    ll = Index(ii, ij, 0);
		    data1[l] -= imag;
		    data1[l+1] += real;
		    data1[ll] += imag;
		    data1[ll+1] += real;
		} else {
		    data1[l] = real;
		    data1[l+1] = imag;
		}
	    }
	}
    }

    isign = -1;

    fp = fopen("rndtest3in", "w");
    for (i=0;i<nn[0];i++) {
	for (j=0;j<nn[1];j++) {
	    for (k=0;k<nn[2];k++) {
		l = Index(i, j, k);
		fwrite(data1+l, 1, sizeof(float), fp);
	    }
	}
    }
    fclose(fp);

    if (isign == 1)
      fftndrf(NDIM, nn, data1);
    else
      fftndrb(NDIM, nn, data1);

    fp = fopen("rndtest3out", "w");

    for (i=0;i<nn[0];i++) {
	for (j=0;j<nn[1];j++) {
	    for (k=0;k<nn[2];k++) {
		l = Index(i,j,k);
		fwrite(&data1[l], 1, sizeof(float), fp);
	    }
	}
    }
    fclose(fp);

    exit(0);
}
