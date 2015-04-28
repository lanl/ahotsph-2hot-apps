#define NO_MSGS

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

static float eps2;
static float GNewt;
static int Nobj;
static int add_offset;
static float offset[NDIM];

void
SetGravOffset(float *off)
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
SetTol(float tol, float frac_tol, float newton_const, float eps, int gnobj)
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

#if defined(__T3D__)
#undef INTERACTF
#define INTERACTF do_cheb_u2
#endif

#if defined(__alpha)
#undef INTERACTF
#undef UNROLL
#define INTERACTF do_grav_cheb_unroll3
#define UNROLL 3
#endif

/* This should be dynamically extensible */
#define IVECSZ 81920
struct {
    float mass;
    float pos[NDIM];
} Ivec[IVECSZ];

void DLmacv(Sink *sink, const hcell **source_vec, int *result, int n);
void Nlognmacv(Sink *sink, const hcell **source_vec, int *result, int);
void INTERACTF(const float *p, const float *end, const float *pos0, 
	       float *mass0, float *acc0, float *phi0, const float *eps2p, 
	       int *ncut);

void 
InheritSinkNlogN(const Sink *from, Sink *to, hcell *pp)
{
    if (to == NULL) {
	body *bp = pp->ptr;
	/* must init mtot or else you get quiet exceptions in asm code */
	float mtot = (float)0.0;
	int ijunk = 0, nn;
	float acc[NDIM];
	float phi;

	VS(acc, = (float)0.0);
	phi = (float)0.0;
	/* putting a getrusage based timer here can slow things down a lot */
	StartTimer(&GravTm);
#ifdef UNROLL
#ifdef __DELTA__
	while (from->icnt % UNROLL) {
	    Ivec[from->icnt].mass = (float)0.0;
	    VS(Ivec[from->icnt].pos, = (float)0.0);
	    from->icnt++;
	}
#endif
	nn = from->icnt - (from->icnt % UNROLL);

#ifdef __DELTA__
	/* Use the interface to the fast assembly code */
	if ((int)&Ivec[0] & 07 || (int)&Ivec[1] & 07)
	  Error("Ivec not aligned for asm code\n");
#endif
	INTERACTF((float *)&Ivec[0], (float *)&Ivec[nn], 
		  from->pos, &mtot, acc, &phi, &eps2, &ijunk);
	if (from->icnt % UNROLL)
	  do_grav((float *)&Ivec[nn], (float *)&Ivec[from->icnt], 
		  from->pos, &mtot, acc, &phi, &eps2, &ijunk);
#else
	do_grav((float *)&Ivec[0], (float *)&Ivec[from->icnt], 
		from->pos, &mtot, acc, &phi, &eps2, &ijunk);
#endif
	Msgf(("%4d interactions for %s\n", from->icnt, PrintKey(pp->key)));
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
	to->bmax = (float)0.0;
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
	to->M0 = (float)0.0;
	VS(to->M1, = (float)0.0);
	to->icnt = 0;
    }
}

void
RcritMAC(Sink *sink, const hcell **source_vec, int *result, int n)
{
    VxdV(float pos_sink, = sink->pos);
    Vxd(float a);
    int icnt = sink->icnt;
    int interactions = 0;
    float dr2;
    Vxd(float pos);
    Vxd(float r);
    float mass;
    float rcrit;
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

#define RcritFac ((float)4.0)	/* should be >= 2.0 */

/* RcritMAC with Don't Laugh-like traversal */
void
DLRcritMAC(Sink *sink, const hcell **source_vec, int *result, int n)
{
    VxdV(float pos_sink, = sink->pos);
    float bmax = sink->bmax;
    int icnt = sink->icnt;
    int interactions = 0;
    int daughters;
    float dr2;
    float mass;
    Vxd(float r);
    Vxd(float dx);
    float rcrit_bmax;
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
	    Msgf(("B %4d %s\n", icnt, PrintKey(source_vec[i]->key)));
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
	    Msgf(("A %4d %s\n", icnt, PrintKey(source_vec[i]->key)));
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



