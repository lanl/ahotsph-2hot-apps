#include <stdio.h>
#include <stdlib.h>

float ran1(long *idum);
void fftndrf(int ndim, int npts[], float *data);
void fftndrb(int ndim, int npts[], float *data);
void fftndcf(int ndim, int npts[], float *data);
void fftndcb(int ndim, int npts[], float *data);

#define NDIM 3

#define Index(i,j,k) (2*(((i)*nn[1]+(j))*nn[2]+(k)))

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

    printf ("doing %dx%dx%d transform\n", nn[0],nn[1],nn[2]);

    ndat2 = nn[0]*nn[1]*nn[2]*2;

    data1=malloc(ndat2*sizeof(float));

    for (i=0;i<nn[0];i++) {
	for (j=0;j<nn[1];j++) {
	    for (k=0;k<nn[2];k++) {
		l = Index(i,j,k);
		data1[l]=2*ran1(&idum)-1;
		l++;
		data1[l]=2*ran1(&idum)-1;
	    }
	}
    }
    isign = 1;
    if (isign == 1)
      fftndcb(NDIM, nn, data1);
    else
      fftndcf(NDIM, nn, data1);

    fp = fopen("ndtest3out", "w");

    for (i=0;i<nn[0];i++) {
	for (j=0;j<nn[1];j++) {
	    for (k=0;k<nn[2];k++) {
		l = Index(i,j,k);
		fwrite(&data1[l], 1, sizeof(float), fp);
		fwrite(&data1[l+1], 1, sizeof(float), fp);
	    }
	}
    }
    fclose(fp);

    exit(0);
}
