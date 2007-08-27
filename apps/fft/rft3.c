#include <stdio.h>
#include <stdlib.h>
#include "bigmalloc.h"
#include "Msgs.h"
#include "timers.h"
#include "protos.h"
#include "singlio.h"
#include "mpmy.h"
#include "mpmy_io.h"
#define NRANSI
#include "complex.h"

/* From nrutil.h and nr.h */
float *vector(long nl, long nh);
void free_vector(float *v, long nl, long nh);
unsigned long *lvector(long nl, long nh);
void free_lvector(unsigned long *v, long nl, long nh);
float ran1(long *idum);
void four1(float data[], unsigned long nn, int isign);
void realft(float data[], unsigned long nn, int isign);

/* local */
static void fourP(fcomplex data[], unsigned long nn[], int ndim, int isign);
static void vcopy(fcomplex *y, fcomplex *x, unsigned int n, 
      unsigned int ystride, unsigned int xstride);
static void vvcopy(fcomplex *y, fcomplex *x, unsigned int n, unsigned int size,
      unsigned int ystride, unsigned int xstride);

static fcomplex *direct_untranspose(void);
static fcomplex *direct_transpose(void);


#define NC (sizeof(fcomplex)/sizeof(float))

#define NDIM 3

#define Index(i,j,k) ((((i)*nn[1]+(j))*nn[2]+(k)))

static int nproc;
static int procnum;
static fcomplex *tptr;
static unsigned int xn, xp, xnum;

Timer_t Tot, TotWC, Four1, FourP;
Timer_t CommTm;

void
main(int argc, char *argv[])
{
    int isign;
    long idum=(-23);
    unsigned long i,j,k,l,*nn;
    fcomplex *data1;
    int ndat2;
    int ii;
    MPMYFile *fp;

    MPMY_Init(&argc, &argv);
    MsgdirInit("msgs");
    nproc = MPMY_Nproc();
    procnum = MPMY_Procnum();

    EnableTimer(&Tot, "Total");
    EnableWCTimer(&TotWC, "Total(WC)");
    EnableTimer(&Four1, "Four1");
    EnableTimer(&FourP, "FourP");
    EnableTimer(&CommTm, "Comm");
    ClearEnabledTimers();
    ClearEnabledCounters();

    StartTimer(&Tot);
    StartTimer(&TotWC);

    nn=lvector(0,NDIM-1);

    nn[0] = atoi(argv[1]);
    nn[1] = atoi(argv[2]);
    nn[2] = atoi(argv[3]);


    nn[2] /= 2;			/* real transform */
    ndat2 = nn[0]*nn[1]*nn[2]/nproc;
    idum -= procnum;


    singlPrintf ("doing %dx%dx%d real transform\n", nn[0],nn[1],nn[2]*2);
    singlPrintf ("int nproc = %d;\n", nproc);
    singlPrintf ("Estimate %dk memory use per node\n", 
		 (2*ndat2*sizeof(fcomplex))/1024);

    data1=(fcomplex *)vector(0,NC*ndat2-1);

    if (nn[0] < nproc)
      Error("1st dimension smaller than Nproc\n");
    if (nn[0] % nproc)
      Error("1st dimension not divisible by nproc\n");

    memset(data1, ndat2*sizeof(fcomplex), 0);

    for (ii=procnum*nn[0]/nproc; ii<(procnum+1)*nn[0]/nproc; ii++) {
	float real, imag;
	i = ii % (nn[0]/nproc);
	for (j=0;j<nn[1];j++) {
	    for (k=0;k<=nn[2];k++) {
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
		    data1[l].r += real;
		    data1[l].i += imag;
		    data1[ll].r += real;
		    data1[ll].i -= imag;
		} else if (k == nn[2]) {
		    int ii, ij, ll;
		    l = Index(i,j,0);
		    ii = i;
		    ij = j;
		    if (ii) ii = nn[0]-ii;
		    if (ij) ij = nn[1]-ij;
		    ll = Index(ii, ij, 0);
		    data1[l].r -= imag;
		    data1[l].i += real;
		    data1[ll].r += imag;
		    data1[ll].i += real;
		} else {
		    data1[l].r = real;
		    data1[l].i = imag;
		}
	    }
	}
    }
    fp = MPMY_Fopen("rft3in", 
		    MPMY_WRONLY | MPMY_CREAT | MPMY_TRUNC | MPMY_MULTI);

    MPMY_Fwrite(data1, ndat2, sizeof(fcomplex), fp);
    MPMY_Fclose(fp);

    isign = -1;
    fourP(data1-1,nn-1,NDIM,isign);

    StopTimer(&Tot);
    StopTimer(&TotWC);

    fp = MPMY_Fopen("rft3out", 
		    MPMY_WRONLY | MPMY_CREAT | MPMY_TRUNC | MPMY_MULTI);

    MPMY_Fwrite(data1, ndat2, sizeof(fcomplex), fp);
    MPMY_Fclose(fp);

    Msg_flush();
    OutputTimers(singlPrintf);
    OutputCounters(singlPrintf);
    exit(0);
}

static void
fourP(fcomplex data[], unsigned long nn[], int ndim, int isign)
{
    int nmax;
    fcomplex *ds;
    fcomplex *p, *q;
    fcomplex *buf;
    int sendproc;

    nn++;			/* for compatibility with fourn */
    data++;
    nmax = nn[0]*nn[1]*nn[2]/nproc;

    StartTimer(&FourP);

    if (isign == 1) {
	StartTimer(&Four1);
	for (q = data; q < data + nmax; q += nn[2]) {
	    realft((float *)q-1,2*nn[2],isign);
	}
	StopTimer(&Four1);
    }

    q = (fcomplex *)vector(0,NC*nn[1]-1);
    for (ds = data; ds < data + nmax; ds += nn[1] * nn[2]) {
	for (p = ds; p < ds + nn[2]; p++) {
	    vcopy(q, p, nn[1], 1, nn[2]);
	    StartTimer(&Four1);
	    four1((float *)q-1,nn[1],isign);
	    StopTimer(&Four1);
	    vcopy(p, q, nn[1], nn[2], 1);
	}
    }
    free_vector((float *)q,0,NC*nn[1]-1);


    q = (fcomplex *)vector(0,NC*nn[0]-1);
    if (nproc == 1) {
	for (p = data; p < data + nn[1] * nn[2] / nproc; p++) {
	    vcopy(q, p, nn[0], 1, nn[1]*nn[2]);

	    StartTimer(&Four1);
	    four1((float *)q-1,nn[0],isign);
	    StopTimer(&Four1);

	    vcopy(p, q, nn[0], nn[1]*nn[2], 1);
	}
    } else {
	tptr = data;
	xp = nn[0]/nproc;
	xn = nn[1]*nn[2]/nproc;
	xnum = nn[1]*nn[2];
	buf = direct_transpose();
	memcpy(data, buf, nmax * sizeof(fcomplex));
	Free(buf);

	for (p = data; p < data + nn[1] * nn[2] / nproc; p++) {
	    vcopy(q, p, nn[0], 1, xn);

	    StartTimer(&Four1);
	    four1((float *)q-1,nn[0],isign);
	    StopTimer(&Four1);

	    vcopy(p, q, nn[0], xn, 1);
	}

	tptr = data;
	direct_untranspose();
	for (sendproc = 0; sendproc < nproc; sendproc++) {
	    vvcopy(data+xn*sendproc, buf+xp*xn*sendproc, xp, 
		   xn, xnum, xn);
	}
	Free(buf);

    }
    free_vector((float *)q,0,NC*nn[0]-1);

    if (isign == -1) {
	StartTimer(&Four1);
	for (q = data; q < data + nmax; q += nn[2]) {
	    realft((float *)q-1,2*nn[2],isign);
	}
	StopTimer(&Four1);
    }

    StopTimer(&FourP);
}

#define ATYPE 2000

static fcomplex *
direct_transpose(void) 
{
    int i;
    int xsendproc, xrecvproc;
    int sz = xn * xp * sizeof(fcomplex);
    fcomplex *buf = Malloc(nproc*sz);
    fcomplex *small_buf = Malloc(sz);
    MPMY_Comm_request req, inreq;

    for (i = 0; i < nproc; i++) {
	xsendproc = (i+procnum) % nproc;
	xrecvproc = (nproc+procnum-i) % nproc;
	vvcopy(small_buf, tptr+xn*xsendproc, xp, xn, xn, xnum);
	StartTimer(&CommTm);
	MPMY_Irecv(buf+xp*xn*xrecvproc, sz, xrecvproc, ATYPE, &inreq);
	MPMY_Isend(small_buf, sz, xsendproc, ATYPE, &req);
	MPMY_Wait(inreq, 0);
	MPMY_Wait(req, 0);
	StopTimer(&CommTm);
    }
    Free(small_buf);
    return(buf);
}

static fcomplex *
direct_untranspose(void)
{
    int i;
    int xsendproc, xrecvproc;
    int sz = xn * xp * sizeof(fcomplex);
    fcomplex *buf = Malloc(nproc*sz);
    MPMY_Comm_request req, inreq;

    for (i = 0; i < nproc; i++) {
	xsendproc = (i+procnum) % nproc;
	xrecvproc = (nproc+procnum-i) % nproc;
	StartTimer(&CommTm);
	MPMY_Irecv(buf+xp*xn*xrecvproc, sz, xsendproc, ATYPE, &inreq);
	MPMY_Isend(tptr+xp*xn*xsendproc, sz, xsendproc, ATYPE, &req);
	MPMY_Wait(req, 0);
	MPMY_Wait(inreq, 0);
	StopTimer(&CommTm);
    }
    return(buf);
}


static void
vcopy(fcomplex *y, fcomplex *x, unsigned int n, 
      unsigned int ystride, unsigned int xstride)
{
    fcomplex *ylast = y + n * ystride;

    while (y < ylast) {
	*y = *x;
	x += xstride;
	y += ystride;
    }
}


static void
vvcopy(fcomplex *y, fcomplex *x, unsigned int n, unsigned int size,
      unsigned int ystride, unsigned int xstride)
{
    fcomplex *ylast = y + n * ystride;

    while (y < ylast) {
	memcpy(y, x, size * sizeof(fcomplex));
	x += xstride;
	y += ystride;
    }
}

