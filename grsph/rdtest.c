#include <math.h>
#include "randoms.h"
#include "physics_sph.h"
#include "bigmalloc.h"
#include "vop.h"
#include "fastflpt.h"
#include "mpmy.h"

static ran_state ranstate;
static void testdata (bodyptr btab, int nobj, int gnobj, 
			 int seed, int cencon);
static void pickvec (float *x, float *y, int cf);
static double xrand(double hi, double lo);

void
RdTest(bodyptr *btabp, int gnobj, int *nobjp, int seed, int cencon)
{
    int start, leftover;
    int nobj;
    bodyptr btab;
    int i;

    nobj = gnobj / MPMY_Nproc();
    start = nobj*MPMY_Procnum();
    leftover = gnobj - nobj*MPMY_Nproc();
    if (MPMY_Procnum() < leftover){
	++nobj;
	start += MPMY_Procnum();
    }else{
	start += leftover;
    }
    btab = (bodyptr) Malloc(nobj*sizeof(body));
    testdata(btab, nobj, gnobj, seed, cencon);
    /* Fill in idents. */
    for(i=0; i<nobj; i++){
	btab[i].ident = start+i;
	btab[i].nterms = 1;
    }
    *nobjp = nobj;
    *btabp = btab;
}

static void
testdata(bodyptr btab, int nobj, int gnobj, int seed, int cencon)
{
    bodyptr p;
    int cnt = 0;
    float h = .47 * pow((float)8.5/gnobj, .333333);
    float sys_rho = 1e-6/(.47*.47*.47);
    float sys_u = .1 * pow(sys_rho, .6666);
 
    ran_init(seed*(MPMY_Procnum()+1), &ranstate);
    for (p = &btab[0]; p < &btab[nobj]; cnt++, p++) {
	p->mass = 1e-6 / gnobj;		 /*   set masses equal */
	p->u = sys_u;
	p->h = h;
	p->gama = 1.0;
	p->enth = 1.0 + (5./3.) * p->u;
	pickvec(p->pos, p->vel, 0); /*   pick position, velocity */
	VS(p->pos, *= .47);
	VS(p->vel, = (float)0.0);
	VS(p->mom, = (float)0.0);
	p->pos[0] += 47.0;
    }
}

static void
pickvec(float *x, float *y, int cf)
{
    int k;
    float rsqx, rsqy, xsc;
 
    do {
	rsqx = (float)0.0;
	rsqy = (float)0.0;
	for (k = 0; k < NDIM; k++) {
            x[k] = xrand(-1.0, 1.0);    /* pick a point at random */
            y[k] = xrand(-1.0, 1.0);    /* pick a point at random */
            rsqx += x[k] * x[k];
            rsqy += y[k] * y[k];
	}
        xsc = xrand(0.0, 1.0);
    } while (rsqx > (float)1.0 || rsqy > (float)1.0);
    if (cf)
      VS(x, *= xsc/sqrtf_fast(rsqx));        /* M(r) proportional to r */
}

/*
 * XRAND: generate floating-point random number.
 */
 
static double
xrand(double xl, double xh)
{
    return (xl + (xh - xl) * uniform_rand(&ranstate));
}
