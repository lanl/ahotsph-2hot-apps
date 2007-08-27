#include <stdio.h>
#include <stdlib.h>
#include "bigmalloc.h"
#include "Msgs.h"
#include "timers.h"
#include "protos.h"
#include "singlio.h"
#include "mpmy.h"
#include "Assert.h"
#include "randoms.h"
#define NRANSI
#include "complex.h"

/* From nrutil.h and nr.h */
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
static void exchange_cdata(fcomplex *cdata, unsigned long *nn);

void ranp_reset(int i, int n, ran_state *rs);

#define MAXNPROC 64

#define NC (sizeof(fcomplex)/sizeof(float))

#define NDIM 3

#define Index(i,j,k) ((((i)*nn[1]+(j))*half_nn2+(k)))
#define ConjIndex(i,j) \
    (Index((i) ? nn[0]-i : 0, (j) ? nn[1]-j : 0, 0)/half_nn2)
#define ConjIndex2(i,j) \
    (Index(i, (j) ? nn[1]-j : 0, 0)/half_nn2)

static int nproc;
static int procnum;
static fcomplex *tptr;
static unsigned int xn, xp, xnum;

Timer_t Four1, FourP, CommTm;

static int 
maxheap(void)
{
    int memused = malloc_heapsz()/1024;
    MPMY_Combine(&memused, &memused, 1, MPMY_INT, MPMY_MAX);
    return memused;
}

static int 
maxmem(void)
{
    int memused = malloc_used()/1024;
    MPMY_Combine(&memused, &memused, 1, MPMY_INT, MPMY_MAX);
    return memused;
}


/* Note: This transform multiplies output data by 2.0 */
/* It is then entirely consistent with the complex fft in ft.c */

void
prft(unsigned long *nn, float *data, ran_state *rs, 
     void spectrum(int, int, int, float *, float *))
{
    int isign;
    unsigned long i,j,k,l;
    fcomplex *data1;
    fcomplex *cdata;
    int ndat2;
    int ip;
    int ll;
    int is, js;
    float real, imag;
    long half_nn2;

    nproc = MPMY_Nproc();
    procnum = MPMY_Procnum();

    EnableTimer(&Four1, "Four1");
    EnableTimer(&FourP, "FourP");
    EnableTimer(&CommTm, "Comm");

    /* There was a very bad bug in previous versions of this code */
    /* Where I divided nn[2] by 2 here, thus modifying it non-locally */
    half_nn2 = nn[2]/2;			/* real transform */
    ndat2 = nn[0]*nn[1]*half_nn2/nproc;

    singlPrintf ("doing %dx%dx%d real transform\n", nn[0],nn[1],nn[2]);
    singlPrintf ("int nproc = %d;\n", nproc);
    singlPrintf ("Estimate 2x%dk memory use per node for fft\n", 
		 (ndat2*sizeof(fcomplex))/1024);

    data1=(fcomplex *)data;
    cdata=Calloc(ndat2/half_nn2, sizeof(fcomplex));

    if (nn[0] < nproc)
      Error("1st dimension smaller than Nproc\n");
    if (nn[0] % nproc)
      Error("1st dimension not divisible by nproc\n");

    for (ip=procnum*nn[0]/nproc; ip<(procnum+1)*nn[0]/nproc; ip++) {
	i = ip % (nn[0]/nproc);
	is = (ip < nn[0]/2) ? ip : nn[0]-ip;
	/* This scheme doesn't work if nn[0] is less than MAXNPROC */
	ranp_reset(i, nn[0], rs);
	for (j=0;j<nn[1];j++) {
	    js = (j < nn[1]/2) ? j : nn[1]-j;
	    spectrum(is, js, 0, &real, &imag);
	    l = Index(i,j,0);	/* k = 0 mode */
	    ll = l/half_nn2;
	    data1[l].r = real;
	    data1[l].i = imag;
	    cdata[ll].r += real;
	    cdata[ll].i -= imag;
	
	    for (k=1;k<half_nn2;k++) {
		spectrum(is, js, k, &real, &imag);
		l = Index(i,j,k);
		data1[l].r = real;
		data1[l].i = imag;
	    }

	    spectrum(is, js, half_nn2, &real, &imag);
	    l = Index(i,j,0);	/* k = half_nn2 mode */
	    ll = l/half_nn2;
	    data1[l].r -= imag;
	    data1[l].i += real;
	    cdata[ll].r += imag;
	    cdata[ll].i += real;
	}
    }

    /* Assure symmetry relation for k = 0 and half_nn2 modes */
    exchange_cdata(cdata, nn);
    for (ip=procnum*nn[0]/nproc; ip<(procnum+1)*nn[0]/nproc; ip++) {
	i = ip % (nn[0]/nproc);
	for (j=0;j<nn[1];j++) {
	    l = Index(i,j,0);
	    ll = ConjIndex2(i,j);
	    data1[l].r += cdata[ll].r;
	    data1[l].i += cdata[ll].i;
	}
    }
    Free(cdata);

    isign = -1;
    fourP(data1-1,nn-1,NDIM,isign);

    for (i = 0; i < ndat2; i++) {
	data1[i].r *= (float)2.0;
	data1[i].i *= (float)2.0;
    }
}

static void
fourP(fcomplex data[], unsigned long nn[], int ndim, int isign)
{
    int nmax;
    fcomplex *ds;
    fcomplex *p, *q;
    fcomplex *buf;
    int sendproc;
    int half_nn2;

    nn++;			/* for compatibility with fourn */
    data++;
    half_nn2 = nn[2]/2;
    nmax = nn[0]*nn[1]*half_nn2/nproc;

    StartTimer(&FourP);

    if (isign == 1) {
	StartTimer(&Four1);
	for (q = data; q < data + nmax; q += half_nn2) {
	    realft((float *)q-1,nn[2],isign);
	}
	StopTimer(&Four1);
    }

    q = Malloc(nn[1]*sizeof(fcomplex));
    for (ds = data; ds < data + nmax; ds += nn[1] * half_nn2) {
	for (p = ds; p < ds + half_nn2; p++) {
	    vcopy(q, p, nn[1], 1, half_nn2);
	    StartTimer(&Four1);
	    four1((float *)q-1,nn[1],isign);
	    StopTimer(&Four1);
	    vcopy(p, q, nn[1], half_nn2, 1);
	}
    }
    Free(q);

    q = Malloc(nn[0]*sizeof(fcomplex));
    if (nproc == 1) {
	for (p = data; p < data + nn[1] * half_nn2 / nproc; p++) {
	    vcopy(q, p, nn[0], 1, nn[1]*half_nn2);

	    StartTimer(&Four1);
	    four1((float *)q-1,nn[0],isign);
	    StopTimer(&Four1);

	    vcopy(p, q, nn[0], nn[1]*half_nn2, 1);
	}
    } else {
	tptr = data;
	xp = nn[0]/nproc;
	xn = nn[1]*half_nn2/nproc;
	xnum = nn[1]*half_nn2;

	/* Theses transposes cost a factor of two in memory */
	buf = direct_transpose();
	memcpy(data, buf, nmax * sizeof(fcomplex));
	Free(buf);

	for (p = data; p < data + nn[1] * half_nn2 / nproc; p++) {
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
    Free(q);

    if (isign == -1) {
	StartTimer(&Four1);
	for (q = data; q < data + nmax; q += half_nn2) {
	    realft((float *)q-1,nn[2],isign);
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
#ifndef __DELTA
	MPMY_Irecv(buf+xp*xn*xrecvproc, sz, xrecvproc, ATYPE, &inreq);
	MPMY_Isend(small_buf, sz, xsendproc, ATYPE, &req);
#else
	MPMY_Irecv(buf+xp*xn*xrecvproc, sz, xrecvproc, ATYPE+xrecvproc,&inreq);
	MPMY_Isend(small_buf, sz, xsendproc, ATYPE+procnum, &req);
#endif
	MPMY_Wait2(inreq, 0, req, 0);
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
	/* Could this be made into a shift? */
#ifndef __DELTA__
	MPMY_Irecv(buf+xp*xn*xrecvproc, sz, xrecvproc, ATYPE, &inreq);
	MPMY_Isend(tptr+xp*xn*xsendproc, sz, xsendproc, ATYPE, &req);
#else
	MPMY_Irecv(buf+xp*xn*xrecvproc, sz, xrecvproc, ATYPE+xrecvproc,&inreq);
	MPMY_Isend(tptr+xp*xn*xsendproc, sz, xsendproc, ATYPE+procnum, &req);
#endif
	MPMY_Wait2(req, 0, inreq, 0);
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

static void
exchange_cdata(fcomplex *cdata, unsigned long *nn)
{
    fcomplex *buf;
    int i, ip, dest;
    int sz = nn[1] * sizeof(fcomplex);

    buf = Malloc(sz);
    if (nproc == 1) {	/* Need to reverse the nn[0] co-ordinate anyway */
	for (i=1; i<nn[0]/2; i++) {
	    dest = nn[0]-i;
	    memcpy(buf, cdata+i*nn[1], sz);
	    memcpy(cdata+i*nn[1], cdata+dest*nn[1], sz);
	    memcpy(cdata+dest*nn[1], buf, sz);
	}
    } else if (procnum < nproc/2) {
	for (ip=procnum*nn[0]/nproc; ip<(procnum+1)*nn[0]/nproc; ip++) {
	    i = ip % (nn[0]/nproc);
	    dest = (ip) ? nn[0]-ip : 0;
	    dest /= nn[0]/nproc;
	    if (dest != procnum) {
		MPMY_Shift(dest, buf, sz, cdata+i*nn[1], sz, 0);
		memcpy(cdata+i*nn[1], buf, sz);
	    }
	} 
    } else {
	for (ip=(procnum+1)*nn[0]/nproc-1; ip>=procnum*nn[0]/nproc; ip--) {
	    i = ip % (nn[0]/nproc);
	    dest = (ip) ? nn[0]-ip : 0;
	    dest /= nn[0]/nproc;
	    if (dest != procnum) {
		MPMY_Shift(dest, buf, sz, cdata+i*nn[1], sz, 0);
		memcpy(cdata+i*nn[1], buf, sz);
	    }
	} 
    }
    Free(buf);
}

