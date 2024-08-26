#include "Malloc.h"
#include "Msgs.h"
#include "error.h"
#include "mpmy.h"
#include "physics.h"
#include "randoms.h"
#include "ring.h"
#include "singlio.h"
#include "vop.h"

/* accumulate acc and phi in double precision */
typedef struct {
    float mass;      /* mass of body */
    float pos[NDIM]; /* position of body */
    double acc[NDIM];
    double phi;
    int64_t ident; /* Not strictly necessary, but could catch errors */
} qbody;

static double L;
static double G;

void EwaldSetup(double box_length, double Gnewt) {
    L = box_length;
    G = Gnewt;
}

#if 0
static void
set_body(void *o, void *p)
{
    qbody *q = o;
    body *b = p;

    q->mass = b->mass;
    VV(q->pos, = b->pos);
    q->ident = b->ident;
    VS(q->acc, = 0.0);
    q->phi = 0.0;
}

static void
finish_body(void *o, void *p)
{
    qbody *q = o;
    body *b = p;

    assert(b->ident == q->ident);
    VV(b->acc, = G*q->acc);	/* double => float */
    b->phi = G*q->phi;
}

static void
grav_ewald(void *p0, void *list, int bsize, int n)
{
    qbody *q = p0;
    body *b, *btab = list;
    double x[NDIM];
    double phi, acc[NDIM];
    double pot, f[NDIM];

    if (bsize != sizeof(body)) Error("Sizes inconsistent\n");

    VS(acc, = 0.0);
    phi = 0.0;
    for (b = btab; b < btab + n; b++) {
	VVV(x, = b->pos, - q->pos);
	ewald(x, L, f, &pot);
	phi += b->mass*pot;
	VV(acc, += b->mass*f);
    }
    VV(q->acc, += acc);
    q->phi += phi;
}
#endif

static int Qgnobj;

void sum_qbody(const void *in0, const void *in1, void *out) {
    int i;
    const qbody *q0 = in0;
    const qbody *q1 = in1;
    qbody *qout = out;

    for (i = 0; i < Qgnobj; i++) {
        VVV(qout[i].acc, = q0[i].acc, +q1[i].acc);
        qout[i].phi = q0[i].phi + q1[i].phi;
    }
}


void EwaldForces(body *btab, int nobj, float sample_frac, void *rs) {
    int i, k, qnobj;
    int *rcount, *roffset;
    float p;
    qbody *q, *qtab, *myqtab;
    int my_offset = 0;
    body *b;
    double x[NDIM];
    double pot, f[NDIM];
    MPMY_Comm_request req;
    ran_state *ranstate = rs;

#if 0
    if (sample_frac == 1.0 || sample_frac == 0.0) {
	/* Ring decomposition, send qbodies around */
	Ring2(btab, sizeof(body), nobj, btab, sizeof(body), nobj,
	      sizeof(qbody), set_body, grav_ewald, finish_body);
	return;
    }
#endif
    /* Keep the global qtab on each processor. */
    /* Since this is an N^2 calculation, I can't imagine a case where it won't fit */
    myqtab = Malloc(nobj * sizeof(qbody));
    p = 2.0 * 1.0 / sample_frac;
    qnobj = 0;
    k = 0;
    while (1) {
        if (sample_frac < 1.0) {
            k += (int)(p * uniform_rand(ranstate));
        }
        if (k < nobj) {
            myqtab[qnobj].mass = btab[k].mass;
            VV(myqtab[qnobj].pos, = btab[k].pos);
            VS(myqtab[qnobj].acc, = 0.0);
            myqtab[qnobj].phi = 0.0;
            myqtab[qnobj].ident = btab[k].ident;
            k++;
            qnobj++;
        } else {
            break;
        }
    }
    myqtab = Realloc(myqtab, qnobj * sizeof(qbody));
    rcount = Malloc(MPMY_Nproc() * sizeof(int));
    roffset = Malloc(MPMY_Nproc() * sizeof(int));
    Native_MPMY_Allgather(&qnobj, 1, MPMY_INT, rcount);

    Qgnobj = rcount[0];
    rcount[0] *= sizeof(qbody);
    roffset[0] = 0;
    for (i = 1; i < MPMY_Nproc(); i++) {
        Qgnobj += rcount[i];
        rcount[i] *= sizeof(qbody);
        roffset[i] = roffset[i - 1] + rcount[i - 1];
    }
    my_offset = roffset[MPMY_Procnum()] / sizeof(qbody);
    Msgf(("qtab my_offset %d, Qgnobj = %d\n", my_offset, Qgnobj));
    singlPrintf("Doing n2 gravity on reduced sample of %d particles\n", Qgnobj);

    qtab = Malloc(Qgnobj * sizeof(qbody));
    Native_MPMY_Allgatherv(myqtab, qnobj * sizeof(qbody), MPMY_CHAR, qtab, rcount, roffset);
    Free(roffset);
    Free(rcount);
    Free(myqtab);

    for (q = qtab; q < qtab + Qgnobj; q++) {
        if ((q - qtab) % 100 == 0) {
            singlPrintf("cycle %5ld of %d in EwaldForces\n", q - qtab, Qgnobj);
            /* Msgf(("cycle %5ld of %d in EwaldForces\n", q-qtab, Qgnobj)); */
        }
        for (b = btab; b < btab + nobj; b++) {
            VVV(x, = b->pos, -q->pos);
            ewald(x, L, f, &pot);
            q->phi += b->mass * pot;
            VV(q->acc, += b->mass * f);
        }
    }
    /* sum acc and phi across procs */
    MPMY_ICombine_Init(&req);
    MPMY_ICombine_func(qtab, qtab, Qgnobj * sizeof(qbody), sum_qbody, req);
    MPMY_ICombine_Wait(req);

    b = btab;
    for (i = 0; i < qnobj; i++) {
        q = qtab + my_offset + i;
        while ((b->ident != q->ident) && (b < btab + nobj)) {
            b->nterms = 1; /* suppress nterms is 0 warning */
            b++;
        }
        if (b >= btab + nobj)
            Error("ident mismatch\n");
        VV(b->acc, = G * q->acc); /* double => float */
        b->phi = G * q->phi;
        b->nterms = Qgnobj;
    }
    for (; b < btab + nobj; b++) { b->nterms = 1; }

    Free(qtab);
    Msgf(("EwaldForces Done\n"));
}
