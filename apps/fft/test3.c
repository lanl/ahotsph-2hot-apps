/* Driver for routine fourn */

#include <stdio.h>
#include <stdlib.h>
#include "protos.h"
#define NRANSI
#include "nr.h"

float *vector(long nl, long nh);
void free_vector(float *v, long nl, long nh);
unsigned long *lvector(long nl, long nh);
void free_lvector(unsigned long *v, long nl, long nh);
float ran1(long *idum);

#define NDIM 3

#define Index(i,j,k) (2*(((i)*nn[1]+(j))*nn[2]+(k)))

int main(int argc, char *argv[])
{
    int isign;
    long idum=(-23);
    unsigned long i,j,k,l,*nn;
    float *data1;
    int ndat2;
    FILE *fp;

    nn=lvector(0,NDIM-1);

    nn[0] = atoi(argv[1]);
    nn[1] = atoi(argv[2]);
    nn[2] = atoi(argv[3]);

    printf ("doing %ldx%ldx%ld transform\n", nn[0],nn[1],nn[2]);

    ndat2 = nn[0]*nn[1]*nn[2]*2;

    data1=vector(0,ndat2-1);

    memset(data1, ndat2, 0);

    for (i=0;i<nn[0];i++) {
	float real, imag;
	int ii, ij, ik;
	ii = i;
	if (ii) ii = nn[0] - ii;
	for (j=0;j<nn[1];j++) {
	    ij = j;
	    if (ij) ij = nn[1] - ij;
	    for (k=0;k<=nn[2]/2;k++) {
		ik = k;
		if (ik) ik = nn[2] - ik;
		real = 2*ran1(&idum)-1;
		imag = 2*ran1(&idum)-1;
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

    fp = fopen("test3in", "w");
    for (i=0;i<nn[0];i++) {
	for (j=0;j<nn[1];j++) {
	    for (k=0;k<nn[2]/2;k++) {
		l = Index(i, j, k);
		fwrite(data1+l, 2, sizeof(float), fp);
	    }
	}
    }
    fclose(fp);

    fourn(data1-1,nn-1,NDIM,isign);

    fp = fopen("test3out", "w");

    for (i=0;i<nn[0];i++) {
	for (j=0;j<nn[1];j++) {
	    for (k=0;k<nn[2];k++) {
		l = Index(i, j, k);
		fwrite(data1+l, 1, sizeof(float), fp);
	    }
	}
    }
    fclose(fp);

    exit(0);
}
