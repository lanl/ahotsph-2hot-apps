/* In this version, I decide that the idea of putting all */
/* possible variations in a single routine is misguided.  */
/* A set of general purpose N-dimensional fft routine. */
/* The array of ints, NPTS[] tells how many points are in the data */
/* array, DATA, in each of NDIM dimensions. */
/* The last two letters, [rc][fb] indicate real/complex and */
/* forward/backward transform, respectively. */
/* The convention is that a forward */
/* FFT takes one from "physical" space to "Fourier" space, while */
/* an backward FFT takes one from Fourier space to physical space. */
/* Neither direction does the division by Npts. */
/* The complex forward transform does :

x[a][b]...[c] = \sum_{i,j,...,k} 
    exp(-2pi\sqrt{-1} ( ai/N[0] + bj/N[1] +... ck/N[ndim-1] )) x[i][j]...[k]
*/
/* The complex backward transform does:
x[a][b]...[c] = \sum_{i,j,...,k} 
    exp(2pi\sqrt{-1} ( ai/N[0] + bj/N[1] +... ck/N[ndim-1] )) x[i][j]...[k]
*/
/* The real forward transform takes a real array and converts it */
/* into a set of complex values.  There are half as many complex values. */
/* Complex fourier space values are stored for all values of the first */
/* ndim-1 indices, and only half of the values of the last index.  */
/* Symmetry relations must be used to obtain the remaining values. */
/* The "signs" in the exponent which define the FT are as in the */
/* complex case shown above. */
/* Denote f_{a,b,...,c} as the actual complex value of the fourier component */
/* with index a,b,...,c.  Let N_c be the number of points in the last */
/* dimension. Let cc(x) be the complex conjugate of x.*/
/* Then the values in a real fourier transform are stored as follows: */
/*
 The values for 1<=c<N_c/2 are stored in the "obvious way":
  f_{a,b,...,c} = x[a][b]...[c]  for all a,b
 For values of c in the upper half, N_c/2 < c < N_c, we use the
"usual symmetry relation" that holds for real Fourier transforms:
  f_{a, b, ..., c} = cc(x[N_a-a][N_b-b]...[N_c-c]) for all a,b

The c=0 and c=N_c/2 components are folded together 
in a slightly complicated way:
 f_{a,b,...,0} = (x[a][b]...[0] + cc(x[N_a-a][N_b-b]...[0]))/2.
 f_{a,b,...,N_c/2} = (x[a][b]...[0] - cc(x[N_a-a][N_b-b]...[0]))/2.

  To go the other way, the complex value actually stored in memory is
 x[a][b]...[0] = f_{a,b,...,0} + i f_{a,b,..., N_c/2}
*/

/*
It is unlikely that any of this works at all if N_c is odd.
Otherwise, there is no restriction on the values of N_i.  Of course,
FFT's always go MUCH faster when N_i are products of small integers.
In fact, the integers should be 2,3,5.  The Fortran routines don't
do use any special tricks to do N=7 or higher butterflies.
*/

/* EXPORTS */
void fft1drf(int npts, float *data);
void fft1drb(int npts, float *data);
void fft1dcf(int npts, float *data);
void fft1dcb(int npts, float *data);
void fftndrf(int ndim, int npts[], float *data);
void fftndrb(int ndim, int npts[], float *data);
void fftndcf(int ndim, int npts[], float *data);
void fftndcb(int ndim, int npts[], float *data);
int offsetnd(int ndim, int npts[], int index[]);
/* ENDEXPORTS */

static void extractr(int n, int skip, float *start, float *dest);
static void replacer(int n, int skip, float *start, float *from);
static void extractc(int n, int skip, float *start, float *dest);
static void replacec(int n, int skip, float *start, float *from);

static void setup_fft1dr(int npts);
static void setup_fft1dc(int npts);

void fftndcf(int ndim, int npts[], float *data){
    int d, n1, n2;
    int nfft, fftpts, fftpts2, nptstot, npts_max;
    int skip, skip2, i1, i2, i;
    float *tmp, *start, *start2;
    
    nptstot = 1;
    npts_max = 0;
    for(d=0; d<ndim-1; d++){	/* Don't do last dimension! */
	nptstot *= npts[d];
	if(npts_max < npts[d]) {
	    npts_max = npts[d];
	}
    }
    tmp = (float *)malloc(2*npts_max*sizeof(float));
    nfft = nptstot;
    nptstot *= npts[ndim-1];

    /* Do the last dimension outside the loop. */
    /* These transforms can be done in place. */
    start = data;
    fftpts = npts[d];
    fftpts2 = 2*fftpts;
    for(i=0; i<nfft; i++){
	fft1dcf(fftpts, start);
	start += fftpts2;
    }

    /* Now do the remaining dimensions. */
    skip = fftpts;
    for(d=ndim-2; d >= 0; d--){
	n1 = skip;
	fftpts = npts[d];
	n2 = nptstot/(fftpts*skip);
	/* This could be optimized a bit more by eliminating the */
	/* integer counters, i1, i2 and just comparing the pointers */
	/* start2, start against their terminal values. */
	/* A smart compiler will probably do this anyway... */
	start2 = data;
	skip2 = 2*fftpts*n1;
	for(i2=0; i2<n2; i2++, start2 += skip2){
	    start = start2;
	    for(i1=0; i1<n1; i1++, start += 2){
		/* The increments on start and start2 are equivalent to: */
		/* start = data + 2*(fftpts*n1*i2 + i1); */
		extractc(fftpts, skip, start, tmp);
		fft1dcf(fftpts, tmp);
		replacec(fftpts, skip, start, tmp);
	    }
	}
	skip *= fftpts;
    }
}

void fftndcb(int ndim, int npts[], float *data){
    int d, n1, n2;
    int nfft, fftpts, fftpts2, nptstot, npts_max;
    int skip, skip2, i1, i2, i;
    float *tmp, *start, *start2;
    
    nptstot = 1;
    npts_max = 0;
    for(d=0; d<ndim-1; d++){	/* Don't do last dimension! */
	nptstot *= npts[d];
	if(npts_max < npts[d]) {
	    npts_max = npts[d];
	}
    }
    tmp = (float *)malloc(2*npts_max*sizeof(float));
    nfft = nptstot;		/* used in final loop.  Do not change! */
    nptstot *= npts[ndim-1];

    skip = nptstot;
    n2 = 1;
    for(d=0; d < ndim-1; d++){
	fftpts = npts[d];
	skip /= fftpts;
	n1 = skip;
	/* This could be optimized a bit more by eliminating the */
	/* integer counters, i1, i2 and just comparing the pointers */
	/* start2, start against their terminal values. */
	/* A smart compiler will probably do this anyway... */
	start2 = data;
	skip2 = 2*fftpts*n1;
	for(i2=0; i2<n2; i2++, start2 += skip2){
	    start = start2;
	    for(i1=0; i1<n1; i1++, start += 2){
		/* The increments on start and start2 are equivalent to: */
		/* start = data + 2*(fftpts*n1*i2 + i1); */
		extractc(fftpts, skip, start, tmp);
		fft1dcb(fftpts, tmp);
		replacec(fftpts, skip, start, tmp);
	    }
	}
	n2 *= fftpts;
    }

    /* Do the last dimension outside the loop. */
    /* These transforms can be done in place. */
    start = data;
    fftpts = npts[ndim-1];
    fftpts2 = 2*fftpts;
    for(i=0; i<nfft; i++){
	fft1dcb(fftpts, start);
	start += fftpts2;
    }

}

void fftndrf(int ndim, int npts[], float *data){
    int d, n1, n2;
    int nfft, fftpts, nptstot, npts_max;
    int skip, skip2, i1, i2, i;
    float *tmp, *start, *start2;
    
    nptstot = 1;
    npts_max = 0;
    for(d=0; d<ndim-1; d++){	/* Don't do last dimension! */
	nptstot *= npts[d];
	if(npts_max < npts[d]) {
	    npts_max = npts[d];
	}
    }
    tmp = (float *)malloc(2*npts_max*sizeof(float));
    nfft = nptstot;
    nptstot *= npts[ndim-1];

    /* Do the last dimension outside the loop. */
    /* These transforms can be done in place. */
    start = data;
    fftpts = npts[d];
    for(i=0; i<nfft; i++){
	fft1drf(fftpts, start);
	start += fftpts;
    }

    /* Now do the remaining dimensions. */
    /* Rely on the last dimension having an even number of points! */
    skip = fftpts/2;
    for(d=ndim-2; d >= 0; d--){
	n1 = skip;
	fftpts = npts[d];
	n2 = nptstot/(2*fftpts*skip);
	/* This could be optimized a bit more by eliminating the */
	/* integer counters, i1, i2 and just comparing the pointers */
	/* start2, start against their terminal values. */
	/* A smart compiler will probably do this anyway... */
	start2 = data;
	skip2 = 2*fftpts*n1;
	for(i2=0; i2<n2; i2++, start2 += skip2){
	    start = start2;
	    for(i1=0; i1<n1; i1++, start+=2){
		/* The increments on start and start2 are equivalent to: */
		/* start = data + (fftpts*n1*i2 + i1); */
		extractc(fftpts, skip, start, tmp);
		fft1dcf(fftpts, tmp);
		replacec(fftpts, skip, start, tmp);
	    }
	}
	skip *= fftpts;
    }
}

void fftndrb(int ndim, int npts[], float *data){
    int d, n1, n2;
    int nfft, fftpts, nptstot, npts_max;
    int skip, skip2, i1, i2, i;
    float *tmp, *start, *start2;
    
    nptstot = 1;
    npts_max = 0;
    for(d=0; d<ndim-1; d++){	/* Don't do last dimension! */
	nptstot *= npts[d];
	if(npts_max < npts[d]) {
	    npts_max = npts[d];
	}
    }
    tmp = (float *)malloc(2*npts_max*sizeof(float));
    nfft = nptstot;		/* used in final loop.  Do not change! */
    nptstot *= npts[ndim-1];

    skip = nptstot/2;
    n2 = 1;
    for(d=0; d < ndim-1; d++){
	fftpts = npts[d];
	skip /= fftpts;
	n1 = skip;
	/* This could be optimized a bit more by eliminating the */
	/* integer counters, i1, i2 and just comparing the pointers */
	/* start2, start against their terminal values. */
	/* A smart compiler will probably do this anyway... */
	start2 = data;
	skip2 = 2*fftpts*n1;
	for(i2=0; i2<n2; i2++, start2 += skip2){
	    start = start2;
	    for(i1=0; i1<n1; i1++, start += 2){
		/* The increments on start and start2 are equivalent to: */
		/* start = data + 2*(fftpts*n1*i2 + i1); */
		extractc(fftpts, skip, start, tmp);
		fft1dcb(fftpts, tmp);
		replacec(fftpts, skip, start, tmp);
	    }
	}
	n2 *= fftpts;
    }

    /* Do the last dimension outside the loop. */
    /* These transforms can be done in place. */
    start = data;
    fftpts = npts[ndim-1];
    for(i=0; i<nfft; i++){
	fft1drb(fftpts, start);
	start += fftpts;
    }
}

static void extractr(int n, int skip, float *start, float *dest){
    while(n--){
	*dest++ = *start;
	start += skip;
    }
}

static void extractc(int n, int skip, float *start, float *dest){
    skip *= 2;
    while(n--){
	*dest++ = *start;
	*dest++ = start[1];
	start += skip;
    }
}

static void replacer(int n, int skip, float *start, float *from){
    while(n--){
	*start = *from++;
	start += skip;
    }
}

static void replacec(int n, int skip, float *start, float *from){	
    skip *= 2;
    while(n--){
	*start = *from++;
	start[1] = *from++;
	start += skip;
    }
}


int offsetnd(int ndim, int npts[], int index[]){
    int offset, d;

    offset = 0;
    for(d=0; d<ndim; d++){
	offset *= npts[d];
	offset += index[d];
    }
    return offset;
}

/* For our inner loop, we call the fortran routines obtained */
/* from the fftpack library on netlib@ornl.  See the documentation */
/* in fftpack.index.  rfftf does forward, real xform.  rfftb does */
/* backward real xform.  rffti initializes.  Similarly for the complex */
/* routines. */

#ifdef risc6000
#define call_fortran(func, args) func args
#define decl_fortran(func, args) func args
#endif

#ifdef alliant
/* It doesn't seem to grok ANSI ## */
#define call_fortran(func, args) func\
_ args
#define decl_fortran(func, args) func\
_ args
#endif

#ifndef call_fortran
/* Hopefully, we've got an ANSI pre-processor */
#define call_fortran(func, args) func##_ args
#define decl_fortran(func, args) func##_ args
#endif

void decl_fortran( rfftf, (int *npts, float *data, float *work));
void decl_fortran( rfftb, (int *npts, float *data, float *work));
void decl_fortran( rffti, (int *npts, float *work));
void decl_fortran( cfftf, (int *npts, float *data, float *work));
void decl_fortran( cfftb, (int *npts, float *data, float *work));
void decl_fortran( cffti, (int *npts, float *work));
static void rfftpack_to_c(int n, float *data);
static void c_to_rfftpack(int n, float *data);

static float *rwork;
static float *cwork;
static int last_nptsr;
static int last_nptsc;

void fft1drf(int npts, float *data){
    if(npts != last_nptsr){
	setup_fft1dr(npts);
	last_nptsr= npts;
    }
    call_fortran( rfftf, (&npts, data, rwork));
    rfftpack_to_c(npts, data);
}

void fft1drb(int npts, float *data){
    if(npts != last_nptsr){
	setup_fft1dr(npts);
	last_nptsr= npts;
    }
    c_to_rfftpack(npts, data);
    call_fortran(rfftb, (&npts, data, rwork));
}

static void setup_fft1dr(int npts){
    if(rwork)
	free((void *)rwork);
    rwork = (float *)malloc((2*npts+15)*sizeof(float));
    call_fortran(rffti, (&npts, rwork));
}

/* The fortran routines expects the transormed array to be an array like */
/* r[0], r[1], c[1], r[2], c[2], ..., r[n/2-2], c[n/2-2], r[n/2-1] */
/* We would rather stuff 0'th and n/2'th element into a single complex like */
/* r[0], r[n/2-1], r[1], c[1], r[2], c[2], ..., r[n/2-2], c[n/2-2] */
/* The question is:  Which is faster: two calls to memcpy (guaranteed */
/* non-overlapping case) or one call to memmove (overlapping by one word) */
/* or a roll-your-own memmove.  (See below) */
static void rfftpack_to_c(int n, float *data)
{
    float last = data[n-1];
#ifdef HAS_MEMMOVE
    (void)memmove((void *)(data+2), (void *)(data+1), (n-2)*sizeof(float)); 
#else
    float *p = data+n-1;
    n -= 2;
    while(n--){
	p[0] = p[-1];
	p--;
    }
#endif
    data[1] = last;
}

static void c_to_rfftpack(int n, float *data)
{
    int nn;
    float last = data[1];
#ifdef HAS_MEMMOVE
    (void)memmove((void *)(data+1), (void *)(data+2), (n-2)*sizeof(float));
#else
    float *p = data+1;
    nn = n-2;
    while(nn--){
	p[0] = p[1];
	p++;
    }
#endif
    data[n-1] = last;
}


void fft1dcf(int npts, float *data){
    if(npts != last_nptsc){
	setup_fft1dc(npts);
	last_nptsc= npts;
    }
    call_fortran(cfftf, (&npts, data, cwork));
}

void fft1dcb(int npts, float *data){
    if(npts != last_nptsc){
	setup_fft1dc(npts);
	last_nptsc= npts;
    }
    call_fortran(cfftb, (&npts, data, cwork));
}

static void setup_fft1dc(int npts){
    if(cwork)
	free((void *)cwork);
    cwork = (float *)malloc((4*npts+15)*sizeof(float));
    call_fortran(cffti, (&npts, cwork));
}

