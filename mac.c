/* #define NOTIMERS */

#include <math.h>
#include "order.h"
#include "physics.h"
#include "vop.h"
#include "fastflpt.h"
#include "Msgs.h"
#include "timers.h"
#include "stk.h"
#include "mpmy.h"

Counter_t CCInt, CBInt, BSInt, BSMax, BCInt, BC2Int, BC4Int, BBInt;
Counter_t FBC2Int, FBC4Int;
Counter_t CCIntRej;
Counter_t TranslateCnt;

Timer_t GravTm, PGravTm, GravSTm, GravMTm, GravQTm, GravHTm, GravQFTm, GravHFTm;
Timer_t MACTm;

static int64_t GNobj, Nobj;
static float Eps2, Eps;
static int Smooth_type;
static float GNewt;
static float DLfac, DLmax;
static int Quad_Ncut = 7;
static int Hexa_Ncut = 20;
#define MAX_IMAGE 125
static int Nimage = 1;
static float offset_array[MAX_IMAGE][NDIM];
static tree_t *SinkTree;
static body *Btab;

static void mxn_hexa(Sink *to, hcell *pp);
static int MxN_hblock = 4*1024;
static int MxN_min_sink = 256;
static int MxN_min_hsrc = 512;
static int MxN_do_pH = 0;

static grav_t Sinteract;
static grav_t Minteract;
static grav_t Qinteract;
static grav_t Hinteract;
#ifdef AMD6100
static int amd6100 = 1;
#else
static int amd6100 = 0;
#endif

void
SetGravOffset(float *off, int n)
{
    int i, j, k;
    int idx;
    int f = 2*n+1;
    
    Nimage = f*f*f;
    if (Nimage > MAX_IMAGE) Error("MAX_IMAGE too small\n");
    for (i = -n; i <= n; i++) {
	for (j = -n; j <= n; j++) {
	    for (k = -n; k <= n; k++) {
		idx = (i+n)*f*f+(j+n)*f+(k+n);
		if (idx >= Nimage) Error("Bad index\n");
		offset_array[idx][0] = i*off[0];
		offset_array[idx][1] = j*off[1];
		offset_array[idx][2] = k*off[2];
	    }
	}
    }
}

#ifdef __AVX__
#define NSSE 8 /* Number of floats in an SSE register */
#define Arch(a) a##_avx8
#else
#define NSSE 4 /* Number of floats in an SSE register */
#define Arch(a) a##_sse4
#endif
/* costs of interaction relative to monopole */
#define QUAD_COST 3
#define HEXA_COST 9

#define SVECSZ (5623) /* This should be dynamically extensible */
struct {
    float mass[NSSE];
    float pos[NDIM][NSSE];
} Svec[SVECSZ];

#define MVECSZ (180073) /* This should be dynamically extensible */
struct {
    float mass[NSSE];
    float pos[NDIM][NSSE];
} Mvec[MVECSZ];

#ifdef QUAD
#define QVECSZ (7919)
struct {
    float mass[NSSE];
    float pos[NDIM][NSSE];
    float R[NSSE];
    float qxx[NSSE];
    float qxy[NSSE];
    float qyy[NSSE];
    float qxz[NSSE];
    float qyz[NSSE];
} Qvec[QVECSZ];
#endif

#ifdef HEXA
#define HVECSZ (24576*2)
struct {
    float mass[NSSE];
    float pos[NDIM][NSSE];
    float R[NSSE];
    float qxx[NSSE];
    float qxy[NSSE];
    float qyy[NSSE];
    float qxz[NSSE];
    float qyz[NSSE];
    float qxxx[NSSE];
    float qxxy[NSSE];
    float qxyy[NSSE];
    float qyyy[NSSE];
    float qxxz[NSSE];
    float qxyz[NSSE];
    float qyyz[NSSE];
    float qxxxx[NSSE];
    float qxxxy[NSSE];
    float qxxyy[NSSE];
    float qxyyy[NSSE];
    float qyyyy[NSSE];
    float qxxxz[NSSE];
    float qxxyz[NSSE];
    float qxyyz[NSSE];
    float qyyyz[NSSE];
} Hvec[HVECSZ];
#endif


void
SetupGrav(float newton_const, float e, int64_t gnobj, float dl_fac, float dl_max,
	  int qcut, int hcut, float particle_mass, int smooth_type)
{
    GNewt = newton_const;
    Smooth_type = smooth_type;
    GNobj = gnobj;
    DLfac = dl_fac;
    DLmax = dl_max;
    Quad_Ncut = qcut;
    Hexa_Ncut = hcut;
    if (smooth_type == 0) {
	Eps = e;
	Sinteract = Minteract = Arch(do_gravp);
	Qinteract = Arch(do_gravpq);
	Hinteract = (amd6100) ? Arch(do_gravph_amd6100) : Arch(do_gravph);
    } else {
	if (smooth_type == 1) {
	    Eps = 1.6f*e;
	    Sinteract = Arch(do_gravsU);
	} else if (smooth_type == 2) {
	    Eps = 1.95f*e;
	    Sinteract = Arch(do_gravsF1);
	} else if (smooth_type == 3) {
	    Eps = 2.2f*e;
	    Sinteract = Arch(do_gravsF2);
	} else if (smooth_type == 4) {
	    Eps = 3.35f*e;
	    Sinteract = Arch(do_gravsK1);
	} else if (smooth_type == 5) {
	    Eps = 2.8f*e;
	    Sinteract = Arch(do_gravsS);
	} else if (smooth_type == 6) {
	    Eps = 3.0f*e;
	    Sinteract = Arch(do_gravsCP);
	} else {
	    Error("Unknown smoothing type\n");
	}
	Minteract = Arch(do_grav);
	Qinteract = Arch(do_gravq);
	Hinteract = (amd6100) ? Arch(do_gravh_amd6100) : Arch(do_gravh);
    }
    Eps2 = Eps*Eps*pow(particle_mass, (float)(2./3.));
}

void Nlognmacv(Sink *sink, const hcell **source_vec, int *result, int);

void
WalkInitSink(tree_t *tp, body *btab, int64_t nobj, int mxn_hblock)
{
    SinkTree = tp;
    Btab = btab;
    Nobj = nobj;
    MxN_hblock = mxn_hblock;
}

void
WalkInitSrc(Stk *kstk, Stk *ostk)
{
    StkPushType(kstk, KeyInt(1), Key_t);
    StkPushType(ostk, offset_array[Nimage/2], float *);
}

void
WalkInitSrcPeriodic(Stk *kstk, Stk *ostk)
{
    int i;

    for (i = 0; i < Nimage; i++) {
	StkPushType(kstk, KeyInt(1), Key_t);
	StkPushType(ostk, offset_array[i], float *);
    }
}

body *
FirstBody(hcell *pp)
{
    Key_t k;
    hcell *p;
    tree_t *tp = SinkTree;
    int nsub = 1<<tp->ndim;
    int sub_flags, i;

    if (pp->type & SHARED) return NULL;
    while ((sub_flags = Sub_Flags(pp))) {
	k = KeyLshift(pp->key, tp->ndim);
	for (i = 0; i < nsub; i++) {
	    if (sub_flags & (1 << i)) {
		p = Find(tp, KeyOrInt(k, i));
		if (p == NULL) Error("FirstBody failed\n");
		pp = p;
		break;
	    }
	}
    }
    return pp->ptr;
}

body *
LastBody(hcell *pp)
{
    Key_t k;
    hcell *p;
    tree_t *tp = SinkTree;
    int nsub = 1<<tp->ndim;
    int sub_flags, i;

    if (pp->type & SHARED) return NULL;
    while ((sub_flags = Sub_Flags(pp))) {
	k = KeyLshift(pp->key, tp->ndim);
	for (i = nsub-1; i >= 0; i--) {
	    if (sub_flags & (1 << i)) {
		p = Find(tp, KeyOrInt(k, i));
		if (p == NULL) Error("LastBody failed\n");
		pp = p;
		break;
	    }
	}
    }
    return pp->ptr;
}

void 
InheritSinkNlogN(const Sink *from, Sink *to, hcell *pp)
{
    if (to == NULL) {
	body *bp = pp->ptr;
	/* must init mtot or else you get quiet exceptions in asm code */
	float mtot = 0.0f;
	int ijunk = 0, nn;
	float acc[NDIM];
	float phi;
	float e;

	VS(acc, = 0.0f);
	phi = 0.0f;
	/* Not symmetric if eps varies, so will not be precisely momentum conserving */
	if (Smooth_type == 0 || Smooth_type == 6) {
	    /* plummer, e ~ (rho*eps)^2 */
	    e = Eps*Eps*pow(bp->mass, 2./3.);
	} else {
	    /* polynomial, e ~ 1/(rho*eps) */
	    e = pow(bp->mass, -1./3.)/Eps;
	}
	/* putting a getrusage based timer here can slow things down a lot */
	StartTimer(&GravTm);
#ifdef HEXA
	StartTimer(&GravHTm);
	nn = from->hcnt;
	while (nn % NSSE) {
	    Hvec[nn/NSSE].mass[nn%NSSE] = 0.0f;
	    VS(Hvec[nn/NSSE].pos,[nn%NSSE] = 0.0f);
	    Hvec[nn/NSSE].R[nn%NSSE] = 0.0f;
	    Hvec[nn/NSSE].qxx[nn%NSSE] = 0.0f;
	    Hvec[nn/NSSE].qxy[nn%NSSE] = 0.0f;
	    Hvec[nn/NSSE].qyy[nn%NSSE] = 0.0f;
	    Hvec[nn/NSSE].qxz[nn%NSSE] = 0.0f;
	    Hvec[nn/NSSE].qyz[nn%NSSE] = 0.0f;
	    Hvec[nn/NSSE].qxxx[nn%NSSE] = 0.0f;
	    Hvec[nn/NSSE].qxxy[nn%NSSE] = 0.0f;
	    Hvec[nn/NSSE].qxyy[nn%NSSE] = 0.0f;
	    Hvec[nn/NSSE].qyyy[nn%NSSE] = 0.0f;
	    Hvec[nn/NSSE].qxxz[nn%NSSE] = 0.0f;
	    Hvec[nn/NSSE].qxyz[nn%NSSE] = 0.0f;
	    Hvec[nn/NSSE].qyyz[nn%NSSE] = 0.0f;
	    Hvec[nn/NSSE].qxxxx[nn%NSSE] = 0.0f;
	    Hvec[nn/NSSE].qxxxy[nn%NSSE] = 0.0f;
	    Hvec[nn/NSSE].qxxyy[nn%NSSE] = 0.0f;
	    Hvec[nn/NSSE].qxyyy[nn%NSSE] = 0.0f;
	    Hvec[nn/NSSE].qyyyy[nn%NSSE] = 0.0f;
	    Hvec[nn/NSSE].qxxxz[nn%NSSE] = 0.0f;
	    Hvec[nn/NSSE].qxxyz[nn%NSSE] = 0.0f;
	    Hvec[nn/NSSE].qxyyz[nn%NSSE] = 0.0f;
	    Hvec[nn/NSSE].qyyyz[nn%NSSE] = 0.0f;
	    nn++;
	}
	if ((long long)&Hvec[0] & 0xF || (long long)&Hvec[1] & 0xF)
	  Error("Hvec not aligned for asm code\n");
	if (nn) Hinteract((float *)&Hvec[from->hcnt_done/NSSE], (float *)&Hvec[nn/NSSE], 
			  from->pos, &mtot, acc, &phi, &e, &ijunk);
	AddCounter(&BC4Int, from->hcnt-from->hcnt_done);
	StopTimer(&GravHTm);
#endif
#ifdef QUAD
	StartTimer(&GravQTm);
	nn = from->qcnt;
	while (nn % NSSE) {
	    Qvec[nn/NSSE].mass[nn%NSSE] = 0.0f;
	    VS(Qvec[nn/NSSE].pos,[nn%NSSE] = 0.0f);
	    Qvec[nn/NSSE].R[nn%NSSE] = 0.0f;
	    Qvec[nn/NSSE].qxx[nn%NSSE] = 0.0f;
	    Qvec[nn/NSSE].qxy[nn%NSSE] = 0.0f;
	    Qvec[nn/NSSE].qyy[nn%NSSE] = 0.0f;
	    Qvec[nn/NSSE].qxz[nn%NSSE] = 0.0f;
	    Qvec[nn/NSSE].qyz[nn%NSSE] = 0.0f;
	    nn++;
	}
	if ((long long)&Qvec[0] & 0xF || (long long)&Qvec[1] & 0xF)
	  Error("Qvec not aligned for asm code\n");
	if (nn) Qinteract((float *)&Qvec[from->qcnt_done/NSSE], (float *)&Qvec[nn/NSSE], 
			  from->pos, &mtot, acc, &phi, &e, &ijunk);
	AddCounter(&BC2Int, from->qcnt-from->qcnt_done);
	StopTimer(&GravQTm);
#endif
	StartTimer(&GravMTm);
	nn = from->mcnt;
	while (nn % NSSE) {
	    Mvec[nn/NSSE].mass[nn%NSSE] = 0.0f;
	    VS(Mvec[nn/NSSE].pos,[nn%NSSE] = 0.0f);
	    nn++;
	}
	if ((long long)&Mvec[0] & 0xF || (long long)&Mvec[1] & 0xF)
	  Error("Mvec not aligned for asm code\n");
	if (nn) Minteract((float *)&Mvec[0], (float *)&Mvec[nn/NSSE], 
			  from->pos, &mtot, acc, &phi, &e, &ijunk);
	AddCounter(&BCInt, from->mcnt);
	StopTimer(&GravMTm);
	StartTimer(&GravSTm);
	nn = from->scnt;
	if (from->scnt > BSMax.counter) BSMax.counter = from->scnt;
	while (nn % NSSE) {
	    Svec[nn/NSSE].mass[nn%NSSE] = 0.0f;
	    VS(Svec[nn/NSSE].pos,[nn%NSSE] = 0.0f);
	    nn++;
	}
	if ((long long)&Svec[0] & 0xF || (long long)&Svec[1] & 0xF)
	  Error("Svec not aligned for asm code\n");
	if (nn) Sinteract((float *)&Svec[0], (float *)&Svec[nn/NSSE], 
			  from->pos, &mtot, acc, &phi, &e, &ijunk);
	AddCounter(&BSInt, from->scnt);
	StopTimer(&GravSTm);
	StopTimer(&GravTm);
	if (!finite(acc[0]) || !finite(acc[1]) || !finite(acc[2]) || !finite(phi)) {
	    Error("bad results from do_grav for (%g,%g,%g), ax=%g ay=%g az=%g phi=%g\n", from->pos[0], from->pos[1], from->pos[2],
		  acc[0], acc[1], acc[2], phi);
	}
	
	if (Minteract == do_grav_sse16_ivec_asm) {
	    /* Fix self-interaction for phi from fast asm code */
	    /* If eps is very small, roundoff becomes a problem */
	    phi += bp->mass*recipsqrtf(e);
	}

	/* Make sure these are initialized to zero externally */
	bp->phi += from->M0;
	bp->phi += phi;
	bp->phi *= GNewt;
	VV(bp->acc, -= from->M1);
	VV(bp->acc, += acc);
	VS(bp->acc, *= GNewt);
	bp->nterms += from->nterms + from->scnt + from->mcnt 
#ifdef QUAD
	    + QUAD_COST*from->qcnt 
#endif
#ifdef HEXA
	    + HEXA_COST*from->hcnt
#endif
	    ;
	if (from->interactions != Nimage*GNobj)
	  Error("Ninteract is %ld, should be %ld\n", from->interactions, Nimage*GNobj);
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
	to->bmax = 0.0f;
	to->daughters = 1;
	to->isbody = 1;
    }

    if (from) {
	to->interactions = from->interactions;
	to->nterms = from->nterms;
	to->M0 = from->M0;
	VV(to->M1, = from->M1);
	to->scnt = from->scnt;
	if (to->scnt >= NSSE*SVECSZ) Error("svec overflow\n");
	to->mcnt = from->mcnt;
	if (to->mcnt >= NSSE*MVECSZ) Error("mvec overflow\n");
#ifdef QUAD
	to->qcnt = from->qcnt;
	to->qcnt_done = from->qcnt_done;
	if (to->qcnt >= NSSE*QVECSZ) Error("qvec overflow\n");
#endif
#ifdef HEXA
	to->hcnt = from->hcnt;
	to->hcnt_done = from->hcnt_done;
	if (to->hcnt >= NSSE*HVECSZ) Error("hvec overflow\n");
#endif
	if (MxN_hblock && to->daughters >= MxN_min_sink && 
	    to->hcnt-to->hcnt_done >= MxN_min_hsrc) {
	    mxn_hexa(to, pp);
	}
    } else {
	to->interactions = 0;
	to->nterms = 0;
	to->M0 = 0.0f;
	VS(to->M1, = 0.0f);
	to->mcnt = 0;
	to->scnt = 0;
#ifdef QUAD
	to->qcnt = to->qcnt_done = 0;
#endif
#ifdef HEXA
	to->hcnt = to->hcnt_done = 0;
#endif
    }
}

static void
mxn_hexa(Sink *to, hcell *pp)
{
    body *p;
    body *first = FirstBody(pp);
    body *last = LastBody(pp)+1;
    int i, n0, n1, m_block, block;

    if (!first || !last) return;
    if (first < Btab || last > first+Nobj) Error("first/last out of range\n");
    StartTimer(&GravTm);
    StartTimer(&GravHFTm);
    n0 = to->hcnt_done/NSSE;
    n1 = to->hcnt/NSSE;
    /* Size MxN_hblock for appropriate WalkPoll() latency */
    /* If Walk Defer timer is large, make MxN_hblock smaller */
    if (n1-n0 > MxN_hblock) m_block = 1;
    else if (n1 == n0) Error("mxn_hexa called with n == 0\n");
    else m_block = MxN_hblock / (n1-n0);
    block = m_block;
    for (p = first; p < last; p += m_block) {
	if (p + block > last) block = last-p;
	if (MxN_do_pH) {
	    pHinteract(&p->mass, p->acc, block, sizeof(body)/sizeof(float),
		       (float *)&Hvec[n0], n1-n0);
	} else {
	    float mtot = 0.0f;
	    float e = 0.0f;
	    int ijunk = 0;
	    for (i = 0; i < block; i++) {
		Hinteract((float *)&Hvec[n0], (float *)&Hvec[n1],
			  (p+i)->pos, &mtot, (p+i)->acc, &(p+i)->phi, &e, &ijunk);
	    }
	}
	WalkPoll();
    }
    AddCounter(&FBC4Int, (last-first)*(n1-n0)*NSSE);
    to->hcnt_done = n1*NSSE;

    if (!finite(first->acc[0]) || !finite(first->acc[1]) || !finite(first->acc[2]) || !finite(first->phi)) {
	Error("bad results from do_grav for (%g,%g,%g), ax=%g ay=%g az=%g phi=%g\n", 
	      first->pos[0], first->pos[1], first->pos[2],
	      first->acc[0], first->acc[1], first->acc[2], first->phi);
    }
    StopTimer(&GravHFTm);
    StopTimer(&GravTm);
}

void
RcritMAC(Sink *sink, const hcell **source_vec, const float **offset_vec, int *result, int n)
{
    VxdV(float pos_sink, = sink->pos);
    int mcnt = sink->mcnt;
    int scnt = sink->scnt;
#ifdef QUAD
    int qcnt = sink->qcnt;
#endif
#ifdef HEXA
    int hcnt = sink->hcnt;
#endif
    int64_t interactions = 0;
    float dr2;
    Vxd(float r);
    Vxd(float dx);
    int i;
    const cell *cp = NULL;
    const quadcell *qcp = NULL;
    const hexacell *hcp = NULL;

    StartTimer(&MACTm);
    if (!sink->isbody) {
	for (i = 0; i < n; i++) result[i] = MAC_SPLIT_SINK;
	return;
    }

    for (i = 0; i < n; i++) {
	const hcell *source = source_vec[i];
	if (Sub_Flags(source)) {
	    cp = source->ptr;
	    qcp = source->ptr;
	    hcp = source->ptr;
	} else {
	    const body *bp = source->ptr;
	    VxVV(r, = bp->pos, + offset_vec[i]);
	    VxVxVx(dx, = r, - pos_sink);
	    dr2 = Dotx(dx, dx);
	    if (dr2 > Eps2) {
		Mvec[mcnt/NSSE].mass[mcnt%NSSE] = bp->mass;
		VVx(Mvec[mcnt/NSSE].pos,[mcnt%NSSE] = r);
		mcnt++;
		if (mcnt/NSSE >= MVECSZ) Error("mvec overflow\n");
	    } else if (dr2 > 0.0f) {
		Svec[scnt/NSSE].mass[scnt%NSSE] = bp->mass;
		VVx(Svec[scnt/NSSE].pos,[scnt%NSSE] = r);
		scnt++;
		if (scnt/NSSE >= SVECSZ) Error("svec overflow\n");
	    }
	    interactions++;
	    result[i] = MAC_ACCEPT;
	    continue;
	}

	VxVV(r, = cp->pos, + offset_vec[i]);
	VxVxVx(dx, = r, - pos_sink);
	dr2 = Dotx(dx, dx);

	if (dr2 >= cp->rcrit*cp->rcrit) {
	    Mvec[mcnt/NSSE].mass[mcnt%NSSE] = cp->mass;
	    VVx(Mvec[mcnt/NSSE].pos,[mcnt%NSSE] = r);
	    mcnt++;
	    if (mcnt/NSSE >= MVECSZ) Error("mvec overflow\n");
	    interactions += cp->daughters;
	    result[i] = MAC_ACCEPT;
#ifdef QUAD
	} else if (Quad_Ncut && (cp->daughters >= Quad_Ncut) &&
		   (dr2 > qcp->rcrit_q*qcp->rcrit_q)) {
	    Qvec[qcnt/NSSE].mass[qcnt%NSSE] = qcp->mass;
	    VVx(Qvec[qcnt/NSSE].pos,[qcnt%NSSE] = r);
	    Qvec[qcnt/NSSE].R[qcnt%NSSE] = qcp->R;
	    Qvec[qcnt/NSSE].qxx[qcnt%NSSE] = qcp->qxx;
	    Qvec[qcnt/NSSE].qxy[qcnt%NSSE] = qcp->qxy;
	    Qvec[qcnt/NSSE].qyy[qcnt%NSSE] = qcp->qyy;
	    Qvec[qcnt/NSSE].qxz[qcnt%NSSE] = qcp->qxz;
	    Qvec[qcnt/NSSE].qyz[qcnt%NSSE] = qcp->qyz;
	    qcnt++;
	    interactions += qcp->daughters;
	    result[i] = MAC_ACCEPT;
#endif
#ifdef HEXA
	} else if (Hexa_Ncut && (cp->daughters >= Hexa_Ncut)
		   && (dr2 > hcp->rcrit_h*hcp->rcrit_h)) {
	    Hvec[hcnt/NSSE].mass[hcnt%NSSE] = hcp->mass;
	    VVx(Hvec[hcnt/NSSE].pos,[hcnt%NSSE] = r);
	    Hvec[hcnt/NSSE].R[hcnt%NSSE] = hcp->R;
	    Hvec[hcnt/NSSE].qxx[hcnt%NSSE] = hcp->qxx;
	    Hvec[hcnt/NSSE].qxy[hcnt%NSSE] = hcp->qxy;
	    Hvec[hcnt/NSSE].qyy[hcnt%NSSE] = hcp->qyy;
	    Hvec[hcnt/NSSE].qxz[hcnt%NSSE] = hcp->qxz;
	    Hvec[hcnt/NSSE].qyz[hcnt%NSSE] = hcp->qyz;
	    Hvec[hcnt/NSSE].qxxx[hcnt%NSSE] = hcp->qxxx;
	    Hvec[hcnt/NSSE].qxxy[hcnt%NSSE] = hcp->qxxy;
	    Hvec[hcnt/NSSE].qxyy[hcnt%NSSE] = hcp->qxyy;
	    Hvec[hcnt/NSSE].qyyy[hcnt%NSSE] = hcp->qyyy;
	    Hvec[hcnt/NSSE].qxxz[hcnt%NSSE] = hcp->qxxz;
	    Hvec[hcnt/NSSE].qxyz[hcnt%NSSE] = hcp->qxyz;
	    Hvec[hcnt/NSSE].qyyz[hcnt%NSSE] = hcp->qyyz;
	    Hvec[hcnt/NSSE].qxxxx[hcnt%NSSE] = hcp->qxxxx;
	    Hvec[hcnt/NSSE].qxxxy[hcnt%NSSE] = hcp->qxxxy;
	    Hvec[hcnt/NSSE].qxxyy[hcnt%NSSE] = hcp->qxxyy;
	    Hvec[hcnt/NSSE].qxyyy[hcnt%NSSE] = hcp->qxyyy;
	    Hvec[hcnt/NSSE].qyyyy[hcnt%NSSE] = hcp->qyyyy;
	    Hvec[hcnt/NSSE].qxxxz[hcnt%NSSE] = hcp->qxxxz;
	    Hvec[hcnt/NSSE].qxxyz[hcnt%NSSE] = hcp->qxxyz;
	    Hvec[hcnt/NSSE].qxyyz[hcnt%NSSE] = hcp->qxyyz;
	    Hvec[hcnt/NSSE].qyyyz[hcnt%NSSE] = hcp->qyyyz;
	    hcnt++;
	    interactions += hcp->daughters;
	    result[i] = MAC_ACCEPT;
#endif
	} else {
	    result[i] = MAC_SPLIT_SRC;
	}
    }
    sink->interactions += interactions;
    sink->scnt = scnt;
    sink->mcnt = mcnt;
#ifdef QUAD
    sink->qcnt = qcnt;
#endif
#ifdef HEXA
    sink->hcnt = hcnt;
#endif
    StopTimer(&MACTm);
}


/* RcritMAC with Don't Laugh-like traversal */
void
DLRcritMAC(Sink *sink, const hcell **source_vec, const float **offset_vec, int *result, int n)
{
    VxdV(float pos_sink, = sink->pos);
    float bmax = sink->bmax;
    int scnt = sink->scnt;
    int mcnt = sink->mcnt;
#ifdef QUAD
    int qcnt = sink->qcnt;
#endif
#ifdef HEXA
    int hcnt = sink->hcnt;
#endif
    int64_t interactions = 0;
    int64_t daughters;
    float dr2;
    Vxd(float r);
    Vxd(float dx);
    int i;
    float rcrit;

    StartTimer(&MACTm);
    for (i = 0; i < n; i++) {
	if (Sub_Flags(source_vec[i])) {
	    const cell *cp = source_vec[i]->ptr;
	    daughters = cp->daughters;
	    VxVV(r, = cp->pos, + offset_vec[i]);
	    VxVxVx(dx, = r, - pos_sink);
	    rcrit = cp->rcrit;
	    dr2 = Dotx(dx, dx);
		
	    /* bmax is 0 if sink is a body */
	    if (dr2 > (rcrit + bmax)*(rcrit + bmax)) {
		/* cell-cell or body-cell */
		Mvec[mcnt/NSSE].mass[mcnt%NSSE] = cp->mass;
		VVx(Mvec[mcnt/NSSE].pos,[mcnt%NSSE] = r);
		mcnt++;
		if (mcnt/NSSE >= MVECSZ) Error("mvec overflow\n");
		interactions += daughters;
		result[i] = MAC_ACCEPT;
		continue;
	    }
#ifdef QUAD
	    if (Quad_Ncut && (daughters >= Quad_Ncut)) {
		const quadcell *qcp = source_vec[i]->ptr;
		rcrit = qcp->rcrit_q;
		if (dr2 > (rcrit + bmax)*(rcrit + bmax)) {
		    Qvec[qcnt/NSSE].mass[qcnt%NSSE] = qcp->mass;
		    VVx(Qvec[qcnt/NSSE].pos,[qcnt%NSSE] = r);
		    Qvec[qcnt/NSSE].R[qcnt%NSSE] = qcp->R;
		    Qvec[qcnt/NSSE].qxx[qcnt%NSSE] = qcp->qxx;
		    Qvec[qcnt/NSSE].qxy[qcnt%NSSE] = qcp->qxy;
		    Qvec[qcnt/NSSE].qyy[qcnt%NSSE] = qcp->qyy;
		    Qvec[qcnt/NSSE].qxz[qcnt%NSSE] = qcp->qxz;
		    Qvec[qcnt/NSSE].qyz[qcnt%NSSE] = qcp->qyz;
		    qcnt++;
		    if (qcnt/NSSE >= QVECSZ) Error("qvec overflow\n");
		    interactions += daughters;
		    result[i] = MAC_ACCEPT;
		    continue;
		}
	    }
#endif
#ifdef HEXA
	    if (Hexa_Ncut && (daughters >= Hexa_Ncut)) {
		const hexacell *hcp = source_vec[i]->ptr;
		rcrit = hcp->rcrit_h;
		if (dr2 > (rcrit + bmax)*(rcrit + bmax)) {
		    Hvec[hcnt/NSSE].mass[hcnt%NSSE] = hcp->mass;
		    VVx(Hvec[hcnt/NSSE].pos,[hcnt%NSSE] = r);
		    Hvec[hcnt/NSSE].R[hcnt%NSSE] = hcp->R;
		    Hvec[hcnt/NSSE].qxx[hcnt%NSSE] = hcp->qxx;
		    Hvec[hcnt/NSSE].qxy[hcnt%NSSE] = hcp->qxy;
		    Hvec[hcnt/NSSE].qyy[hcnt%NSSE] = hcp->qyy;
		    Hvec[hcnt/NSSE].qxz[hcnt%NSSE] = hcp->qxz;
		    Hvec[hcnt/NSSE].qyz[hcnt%NSSE] = hcp->qyz;
		    Hvec[hcnt/NSSE].qxxx[hcnt%NSSE] = hcp->qxxx;
		    Hvec[hcnt/NSSE].qxxy[hcnt%NSSE] = hcp->qxxy;
		    Hvec[hcnt/NSSE].qxyy[hcnt%NSSE] = hcp->qxyy;
		    Hvec[hcnt/NSSE].qyyy[hcnt%NSSE] = hcp->qyyy;
		    Hvec[hcnt/NSSE].qxxz[hcnt%NSSE] = hcp->qxxz;
		    Hvec[hcnt/NSSE].qxyz[hcnt%NSSE] = hcp->qxyz;
		    Hvec[hcnt/NSSE].qyyz[hcnt%NSSE] = hcp->qyyz;
		    Hvec[hcnt/NSSE].qxxxx[hcnt%NSSE] = hcp->qxxxx;
		    Hvec[hcnt/NSSE].qxxxy[hcnt%NSSE] = hcp->qxxxy;
		    Hvec[hcnt/NSSE].qxxyy[hcnt%NSSE] = hcp->qxxyy;
		    Hvec[hcnt/NSSE].qxyyy[hcnt%NSSE] = hcp->qxyyy;
		    Hvec[hcnt/NSSE].qyyyy[hcnt%NSSE] = hcp->qyyyy;
		    Hvec[hcnt/NSSE].qxxxz[hcnt%NSSE] = hcp->qxxxz;
		    Hvec[hcnt/NSSE].qxxyz[hcnt%NSSE] = hcp->qxxyz;
		    Hvec[hcnt/NSSE].qxyyz[hcnt%NSSE] = hcp->qxyyz;
		    Hvec[hcnt/NSSE].qyyyz[hcnt%NSSE] = hcp->qyyyz;
		    hcnt++;
		    if (hcnt/NSSE >= HVECSZ) Error("hvec overflow\n");
		    interactions += daughters;
		    result[i] = MAC_ACCEPT;
		    continue;
		}
	    }
#endif
	    if ((bmax > DLmax) || (DLfac * bmax > rcrit)) {
		result[i] = MAC_SPLIT_SINK;
		if (sink->isbody) Error("Trying to split body\n");
	    } else {
		result[i] = MAC_SPLIT_SRC;
	    }
	} else if (sink->isbody) {
	    /* body-body */
	    const body *bp = source_vec[i]->ptr;
	    VxVV(r, = bp->pos, + offset_vec[i]);
	    VxVxVx(dx, = r, - pos_sink);
	    dr2 = Dotx(dx, dx);
	    if (dr2 > Eps2) {
		Mvec[mcnt/NSSE].mass[mcnt%NSSE] = bp->mass;
		VVx(Mvec[mcnt/NSSE].pos,[mcnt%NSSE] = r);
		mcnt++;
		if (mcnt/NSSE >= MVECSZ) Error("mvec overflow\n");
	    } else if (dr2 > 0.0f) {
		Svec[scnt/NSSE].mass[scnt%NSSE] = bp->mass;
		VVx(Svec[scnt/NSSE].pos,[scnt%NSSE] = r);
		scnt++;
		if (scnt/NSSE >= SVECSZ) Error("svec overflow\n");
	    }
	    interactions++;
	    result[i] = MAC_ACCEPT;
	    continue;
	} else {
	    /* cell-body */
	    /* Comparing dr2 with Eps2 only works for body-body interactions */
	    result[i] = MAC_SPLIT_SINK;
	}
    }
    sink->interactions += interactions;
    sink->scnt = scnt;
    sink->mcnt = mcnt;
#ifdef QUAD
    sink->qcnt = qcnt;
#endif
#ifdef HEXA
    sink->hcnt = hcnt;
#endif
    StopTimer(&MACTm);
}
