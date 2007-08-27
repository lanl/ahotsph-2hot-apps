#include <stdio.h>
#include <stdlib.h>
#include "collective.h"
#include "bigmalloc.h"
#include "Msgs.h"
#include "timers.h"
#include "protos.h"
#include "sysdep.h"
#define NRANSI
#include "complex.h"

/* From nrutil.h and nr.h */
float *vector(long nl, long nh);
void free_vector(float *v, long nl, long nh);
unsigned long *lvector(long nl, long nh);
void free_lvector(unsigned long *v, long nl, long nh);
float ran1(long *idum);
void four1(float data[], unsigned long nn, int isign);

/* local */
static void fourt(fcomplex data[], unsigned long nn[], int ndim, int isign);
static void vcopy(fcomplex *y, fcomplex *x, unsigned int n, 
      unsigned int ystride, unsigned int xstride);
static void vvcopy(fcomplex *y, fcomplex *x, unsigned int n, unsigned int size,
      unsigned int ystride, unsigned int xstride);
static void transpose();
static void untranspose();
static void transpose_sort(int proc, char *buf, int size);

static fcomplex *direct_untranspose(void);
static fcomplex *direct_transpose(void);


#define NC (sizeof(fcomplex)/sizeof(float))

#define NDIM 3

#define Index(i,j,k) ((((i)*nn[1]+(j))*nn[2]+(k)))

static fcomplex *tptr;
static unsigned int xn, xp, xnum;

Timer_t Tot, Four1, FourT;
Timer_t Wait, ShiftTm, CommTm;
Counter_t ShiftCalls, BytesRecvd,  BytesSent;
Counter_t ArecvCalls, AsendCalls;

void
main(int argc, char *argv[])
{
    int isign;
    long idum=(-23);
    unsigned long i,j,k,l,*nn;
    fcomplex *data1,*data2,*data3;
    int ndat2;
    int offset;

    SysdepInit(&argc, &argv);
    IOinit("msgs");

    EnableTimer(&Tot, "Total");
    EnableTimer(&Four1, "Four1");
    EnableTimer(&FourT, "FourT");
    EnableTimer(&ShiftTm, "Shift");
    EnableTimer(&CommTm, "Comm");
    EnableTimer(&Wait, "Wait");
    EnableCounter(&BytesSent, "BytesSent");
    EnableCounter(&BytesRecvd, "BytesRecv");
    EnableCounter(&ShiftCalls, "Shift");
    EnableCounter(&ArecvCalls, "Arecv");
    EnableCounter(&AsendCalls, "Asend");
    ClearEnabledTimers();
    ClearEnabledCounters();

    StartTimer(&Tot);

    nn=lvector(0,NDIM-1);

    nn[0] = atoi(argv[1]);
    nn[1] = atoi(argv[2]);
    nn[2] = atoi(argv[3]);

    singlPrintf ("doing %dx%dx%d transform\n", nn[0],nn[1],nn[2]);
    singlPrintf ("int nproc = %d\n", Nproc());

    ndat2 = nn[0]*nn[1]*nn[2]/Nproc();
    idum -= Procnum();

    data1=(fcomplex *)vector(0,NC*ndat2-1);
    data2=(fcomplex *)vector(0,NC*ndat2-1);
    data3=(fcomplex *)vector(0,NC*ndat2-1);

    if (nn[0] < Nproc())
      Error("1st dimension smaller than Nproc()\n");

    offset = Procnum()*nn[0]/Nproc();

    for (i=0;i<nn[0]/Nproc();i++) {
	for (j=0;j<nn[1];j++) {
	    for (k=0;k<nn[2];k++) {
		l = Index(i,j,k);
#if 1
		data1[l].r=data2[l].r=data3[l].r=2*ran1(&idum)-1;
		data1[l].i=data2[l].i=data3[l].i=2*ran1(&idum)-1;
#else
		data1[l].r=data2[l].r=data3[l].r
		  =(i+offset)*nn[1]*nn[2]+j*nn[2]+k;
		data1[l].i=data2[l].i=data3[l].i
		  = -1.0*((i+offset)*nn[1]*nn[2]+j*nn[2]+k);
#endif
	    }
	}
    }
    isign = 1;
    fourt(data2-1,nn-1,NDIM,isign);
    fourt(data3-1,nn-1,NDIM,isign);
    isign = -1;
    fourt(data3-1,nn-1,NDIM,isign);

    if (argc > 4) {
	for (i = 0; i < ndat2; i++) {
	    Msg_do("%12f %12f %12f\n", 
		   data1[i].r, data2[i].r, data3[i].r/(1e-20+data1[i].r));
	    Msg_do("%12f %12f %12f\n", 
		   data1[i].i, data2[i].i, data3[i].i/(1e-20+data1[i].i));
	}
    }
    StopTimer(&Tot);
    Msg_flush();
    OutputTimers(singlPrintf);
    OutputCounters(singlPrintf);
    exit(0);
}

/* #define four1 (void) */

static void
fourt(fcomplex data[], unsigned long nn[], int ndim, int isign)
{
    int nmax;
    fcomplex *ds;
    fcomplex *p, *q;
    fcomplex *buf;
    int sendproc;

    nn++;			/* for compatibility with fourn */
    data++;
    nmax = nn[0]*nn[1]*nn[2]/Nproc();

    StartTimer(&FourT);

    StartTimer(&Four1);
    for (q = data; q < data + nmax; q += nn[2]) {
	four1((float *)q-1,nn[2],isign);
    }
    StopTimer(&Four1);

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
    if (Nproc() == 1) {
	for (p = data; p < data + nn[1] * nn[2] / Nproc(); p++) {
	    vcopy(q, p, nn[0], 1, nn[1]*nn[2]);

	    StartTimer(&Four1);
	    four1((float *)q-1,nn[0],isign);
	    StopTimer(&Four1);

	    vcopy(p, q, nn[0], nn[1]*nn[2], 1);
	}
    } else {
	tptr = data;
	xp = nn[0]/Nproc();
	xn = nn[1]*nn[2]/Nproc();
	xnum = nn[1]*nn[2];
#ifdef USE_ROUTER
	buf = xcom(ROUTER, ~0, transpose);
	transpose_sort(Procnum(), buf, xp*xn*sizeof(fcomplex));
#else
	buf = direct_transpose();
#endif
	memcpy(data, buf, nmax * sizeof(fcomplex));
	Free(buf);

	for (p = data; p < data + nn[1] * nn[2] / Nproc(); p++) {
	    vcopy(q, p, nn[0], 1, xn);

	    StartTimer(&Four1);
	    four1((float *)q-1,nn[0],isign);
	    StopTimer(&Four1);

	    vcopy(p, q, nn[0], xn, 1);
	}

	tptr = data;
#ifdef USE_ROUTER
	buf = xcom(ROUTER, ~0, untranspose);
	transpose_sort(Procnum(), buf, xp*xn*sizeof(fcomplex));
#else
	direct_untranspose();
#endif
	for (sendproc = 0; sendproc < Nproc(); sendproc++) {
	    vvcopy(data+xn*sendproc, buf+xp*xn*sendproc, xp, 
		   xn, xnum, xn);
	}
	Free(buf);

    }
    free_vector((float *)q,0,NC*nn[0]-1);

    StopTimer(&FourT);
}

int arecv (void *inb, int size, int type);
int arecvfrom(void *inb, int size, int type, int *proc);
void asend (void *outb, int outcnt, int dest, int type);

#define ATYPE 2000

static fcomplex *
direct_transpose(void) 
{
    int i;
    int xsendproc, xrecvproc;
    int procnum = Procnum();
    int nproc = Nproc();
    int sz = xn * xp * sizeof(fcomplex);
    fcomplex *buf = Malloc(nproc*sz);
    fcomplex *small_buf = Malloc(sz);

    for (i = 0; i < nproc; i++) {
	xsendproc = (i+procnum) % nproc;
	xrecvproc = (nproc+procnum-i) % nproc;
	vvcopy(small_buf, tptr+xn*xsendproc, xp, xn, xn, xnum);
	StartTimer(&CommTm);
	asend(small_buf, sz, xsendproc, ATYPE+procnum);
	arecvfrom(buf+xp*xn*xrecvproc, sz, ATYPE+xrecvproc, &xrecvproc);
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
    int procnum = Procnum();
    int nproc = Nproc();
    int sz = xn * xp * sizeof(fcomplex);
    fcomplex *buf = Malloc(nproc*sz);

    for (i = 0; i < nproc; i++) {
	xsendproc = (i+procnum) % nproc;
	xrecvproc = (nproc+procnum-i) % nproc;
	StartTimer(&CommTm);
	asend(tptr+xp*xn*xsendproc, sz, xsendproc, ATYPE+procnum);
	arecvfrom(buf+xp*xn*xrecvproc, sz, ATYPE+xrecvproc, &xrecvproc);
	StopTimer(&CommTm);
    }
    return(buf);
}

#ifdef USE_ROUTER

static void
transpose(struct combuf *xout, unsigned int xsendproc)
{
    int sz = xn * xp * sizeof(fcomplex);

    BufExtend(xout, xout->used + sz);
    vvcopy((fcomplex *)(xout->buf+xout->ptr), tptr+xn*xsendproc, xp, 
	   xn, xn, xnum);
    xout->ptr += sz;
}

static void
untranspose(struct combuf *xout, unsigned int xsendproc)
{
    int sz = xn * xp * sizeof(fcomplex);

    BufExtend(xout, xout->used + sz);
    memcpy(xout->buf+xout->ptr, tptr+xp*xn*xsendproc, sz);
    xout->ptr += sz;
}

static void
transpose_sort(int proc, char *buf, int size)
{
    int high_bit = 1, high_mask, low_mask;
    int from, to, temp, i;
    char *swap_buf;

    if (proc == 0)		/* proc == 0 is already sorted */
      return;
    swap_buf = (char *) Malloc(size);
    temp = proc;
    while(temp >>= 1)
      high_bit <<= 1;
    low_mask = high_bit - 1;
    high_mask = ~low_mask << 1;
    
    for (i = 0; i < Nproc()/2; i++) {
	temp = (i & low_mask) | ((i << 1) & high_mask);
	from = ~high_bit & temp;
	to = high_bit | ((temp ^ proc) & (temp | proc));
	memcpy(swap_buf, buf+size*from, size);
	memcpy(buf+size*from, buf+size*to, size);
	memcpy(buf+size*to, swap_buf, size);
    }
    Free(swap_buf);
}

#endif

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

