/* Driver for routine fourn */

#include <stdio.h>
#include <stdlib.h>
#include "protos.h"
#define NRANSI
#include "nr.h"

float ***f3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh);
float **matrix(long nrl, long nrh, long ncl, long nch);
float *vector(long nl, long nh);
void free_vector(float *v, long nl, long nh);
unsigned long *lvector(long nl, long nh);
void free_lvector(unsigned long *v, long nl, long nh);
float ran1(long *idum);

#define NDIM 3

int main(int argc, char *argv[])
{
    int isign;
    long idum=(-23);
    unsigned long i,j,k,*nn;
    float ***data, **speq;
    FILE *fp;
    int i1, i2, i3;

    nn=lvector(0,NDIM-1);

    nn[0] = atoi(argv[1]);
    nn[1] = atoi(argv[2]);
    nn[2] = atoi(argv[3]);

    printf ("doing %ldx%ldx%ld real transform\n", nn[0],nn[1],nn[2]);

    data=f3tensor(1,nn[0],1,nn[1],1,nn[2]);
    speq=matrix(1,nn[0],1,2*nn[1]);

    for (i=1;i<=nn[0];i++) {
	for (j=1;j<=nn[1];j++) {
	    speq[i][2*j-1] = 0;
	    speq[i][2*j] = 0;
	    for (k=1;k<=nn[2];k++) {
		data[i][j][k] = 0;
	    }
	}
    }
    
    for (i=1;i<=nn[0];i++) {
	float real, imag;
	for (j=1;j<=nn[1];j++) {
	    for (k=1;k<nn[2];k+=2) {
		real = 2*ran1(&idum)-1;
		imag = 2*ran1(&idum)-1;
		if (k == nn[2]+1) {
		    speq[i][2*j-1] = real;
		    speq[i][2*j] = imag;
		} else {
		    data[i][j][k] = real;
		    data[i][j][k+1] = imag;
		}
	    }
	}
    }
    isign = -1;

    fp = fopen("rtest3in", "w");
    for (i=1;i<=nn[0];i++) {
	for (j=1;j<=nn[1];j++) {
	    for (k=1;k<=nn[2];k++) {
		fwrite(&data[i][j][k], 1, sizeof(float), fp);
	    }
	}
    }
    fclose(fp);

    rlft3(data, speq, nn[0], nn[1], nn[2], isign);

    fp = fopen("rtest3out", "w");

    for (i=1;i<=nn[0];i++)
      for (j=1;j<=nn[1];j++)
	for (k=1;k<=nn[2];k++)
	  fwrite(&data[i][j][k], 1, sizeof(float), fp);

    fclose(fp);

    fp = fopen("rtest3out2", "w");

    for (i=1;i<=nn[0];i++)
      for (j=1;j<=2*nn[1];j++)
	  fwrite(&speq[i][j], 1, sizeof(float), fp);

    fclose(fp);

    exit(0);
}
