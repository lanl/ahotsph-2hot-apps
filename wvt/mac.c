#ifdef __DELTA__
#define NOTIMERS
#endif

#include "physics.h"
#include "vop.h"
#include "fastflpt.h"
#include "Msgs.h"
#include "timers.h"
#include "stk.h"


Counter_t CCInt, CBInt, BCInt, BBInt;
Counter_t CCIntRej;
Counter_t TranslateCnt;

Timer_t GravTm;
Timer_t MACTm;

static double eps2;
static double GNewt;
static int Nobj;
static int add_offset;
static double offset[NDIM];

void
SetGravOffset(double *off)
{
    VV(offset, = off);
    add_offset = 1;
}

void
UnSetGravOffset(void)
{
    VS(offset, = 0.0);
    add_offset = 0;
}


void
SetTol(double tol, double frac_tol, double newton_const, double eps, int gnobj)
{
    eps2 = eps*eps;
    GNewt = newton_const;
    Nobj = gnobj;
}

/* The fast N log N stuff is below here */

#define UNROLL 2		/* unroll by 2 by default */
#define INTERACTF do_grav_u2

#if defined(__INTEL_SSD__) && !defined(NO_ASM)
#undef UNROLL
#undef INTERACTF
#define UNROLL 3
#define INTERACTF do_grav_fast
#endif 

#if defined(__T3D__) || defined (_IBMR2)
#define USE_CHEB_RSQRT
#endif

#ifdef USE_CHEB_RSQRT
#undef INTERACTF
#define INTERACTF do_cheb_u2
#endif

/* This should be dynamically extensible */
#define IVECSZ 40960
struct {
    double mass;
    double pos[NDIM];
} Ivec[IVECSZ];

void DLmacv(Sink *sink, const hcell **source_vec, int *result, int n);
void Nlognmacv(Sink *sink, const hcell **source_vec, int *result, int);
void INTERACTF(const double *p, const double *end, const double *pos0, 
	       double *mass0, double *acc0, double *phi0, const double *eps2p, 
	       int *ncut);

void 
InheritSinkNlogN(const Sink *from, Sink *to, hcell *pp)
{
    if (to == NULL) {
	body *bp = pp->ptr;
	/* must init mtot or else you get quiet exceptions in asm code */
	double mtot = (double)0.0;
	int ijunk = 0, nn;
	double acc[NDIM];
	double phi;

	VS(acc, = (double)0.0);
	phi = (double)0.0;
	/* putting a getrusage based timer here can slow things down a lot */
	StartTimer(&GravTm);
#ifdef UNROLL
#ifdef __DELTA__
	while (from->icnt % UNROLL) {
	    Ivec[from->icnt].mass = (double)0.0;
	    VS(Ivec[from->icnt].pos, = (double)0.0);
	    from->icnt++;
	}
#endif
	nn = from->icnt - (from->icnt % UNROLL);

	/* Use the interface to the fast assembly code */
#ifdef __INTEL_SSD__
	if ((int)&Ivec[0] & 07 || (int)&Ivec[1] & 07)
	  Error("Ivec not aligned for asm code\n");
#endif
	INTERACTF((double *)&Ivec[0], (double *)&Ivec[nn], 
		  from->pos, &mtot, acc, &phi, &eps2, &ijunk);
	if (from->icnt % UNROLL)
	  do_grav((double *)&Ivec[nn], (double *)&Ivec[from->icnt], 
		  from->pos, &mtot, acc, &phi, &eps2, &ijunk);
#else
	do_grav((double *)&Ivec[0], (double *)&Ivec[from->icnt], 
		from->pos, &mtot, acc, &phi, &eps2, &ijunk);
#endif
	StopTimer(&GravTm);

	/* Make sure these are initialized to zero externally */
	bp->phi += GNewt*from->M0;
	bp->phi += GNewt*phi;
	VV(bp->acc, += -GNewt*from->M1);
	VV(bp->acc, += GNewt*acc);
	bp->nterms += from->nterms + from->icnt;
	if (from->interactions != Nobj)
	  Error("Ninteract is %d, should be %d\n", from->interactions, Nobj);
	return;
    }

    if (Sub_Flags(pp)) {
	cell *cp = pp->ptr;

	VV(to->pos, = cp->pos);
	to->bmax = cp->bmax;
	to->daughters = cp->daughters;
	to->isbody = 0;
    } else {
	body *bp = pp->ptr;
	VV(to->pos, = bp->pos);
	to->bmax = (double)0.0;
	to->daughters = 1.F;
	to->isbody = 1;
    }
    if (add_offset) {
	VV(to->pos, += offset);
    }

    if (from) {
	to->interactions = from->interactions;
	to->nterms = from->nterms;
	to->M0 = from->M0;
	VV(to->M1, = from->M1);
	to->icnt = from->icnt;
	if (to->icnt >= IVECSZ) Error("ivec overflow\n");
    } else {
	to->interactions = 0;
	to->nterms = 0;
	to->M0 = (double)0.0;
	VS(to->M1, = (double)0.0);
	to->icnt = 0;
    }
}

void
RcritMAC(Sink *sink, const hcell **source_vec, int *result, int n)
{
    VxdV(double pos_sink, = sink->pos);
    Vxd(double a);
    int icnt = sink->icnt;
    int interactions = 0;
    double dr2;
    Vxd(double pos);
    Vxd(double r);
    double mass;
    double rcrit;
    int daughters;
    int i;

    StartTimer(&MACTm);
    if (!sink->isbody) {
	for (i = 0; i < n; i++) result[i] = MAC_SPLIT_SINK;
	return;
    }

    VxS(a, = 0.F);

    for (i = 0; i < n; i++) {
	const hcell *source = source_vec[i];
	if (Sub_Flags(source)) {
	    const cell *cp = source->ptr;
	    mass = cp->mass;	/* Access in same order as cell struct */
	    VxV(pos, = cp->pos);	
	    daughters = cp->daughters;
	    rcrit = cp->rcrit;
	} else {
	    const body *bp = source->ptr;
	    Ivec[icnt].mass = bp->mass;
	    VV(Ivec[icnt].pos, = bp->pos);
	    icnt++;
	    interactions++;
	    result[i] = MAC_ACCEPT;
	    continue;
	}

	VxVxVx(r, = pos, - pos_sink);
	dr2 = Dotx(r, r);

	if (dr2 < rcrit*rcrit) {
	    result[i] = MAC_SPLIT_SRC;
	} else {
	    Ivec[icnt].mass = mass;
	    VVx(Ivec[icnt].pos, = pos);
	    icnt++;
	    interactions += daughters;
	    result[i] = MAC_ACCEPT;
	}
    }
    sink->interactions += interactions;
    sink->icnt = icnt;
    StopTimer(&MACTm);
}

#define RcritFac ((double)4.0)	/* should be >= 2.0 */

/* RcritMAC with Don't Laugh-like traversal */
void
DLRcritMAC(Sink *sink, const hcell **source_vec, int *result, int n)
{
    VxdV(double pos_sink, = sink->pos);
    double bmax = sink->bmax;
    int icnt = sink->icnt;
    int interactions = 0;
    int daughters;
    double dr2;
    double mass;
    Vxd(double r);
    Vxd(double dx);
    double rcrit_bmax;
    int i;

    StartTimer(&MACTm);
    for (i = 0; i < n; i++) {
	if (Sub_Flags(source_vec[i])) {
	    const cell *cp = source_vec[i]->ptr;
	    mass = cp->mass;
	    VxV(r, = cp->pos);
	    /* bmax is 0 if sink is a body */
	    rcrit_bmax = cp->rcrit + bmax;
	    daughters = cp->daughters;
	} else {
	    /* cell-body or body-body */
	    const body *bp = source_vec[i]->ptr;
	    Ivec[icnt].mass = bp->mass;
	    VV(Ivec[icnt].pos, = bp->pos);
	    icnt++;
	    interactions += 1;
	    result[i] = MAC_ACCEPT;
	    continue;
	}
	    
	VxVxVx(dx, = r, - pos_sink);
	dr2 = Dotx(dx, dx);
		
	if (dr2 > rcrit_bmax*rcrit_bmax) {
	    /* cell-cell or body-cell */
	    Ivec[icnt].mass = mass;
	    VVx(Ivec[icnt].pos, = r);
	    icnt++;
	    interactions += daughters;
	    result[i] = MAC_ACCEPT;
	} else if (RcritFac * bmax > rcrit_bmax) {
	    result[i] = MAC_SPLIT_SINK;
	    if (sink->isbody) Error("Trying to split body\n");
	} else {
	    result[i] = MAC_SPLIT_SRC;
	}
    }
    sink->interactions += interactions;
    sink->icnt = icnt;
    StopTimer(&MACTm);
}



