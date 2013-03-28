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
#include "ewald_le.h"

Counter_t CCInt, CBInt, BSInt, BSMax, BCInt, BC2Int, BC4Int, BBInt;
Counter_t CEmpty, MCAnti, MCCorr;
Counter_t FBC2Int, FBC4Int;
Counter_t CCIntRej;
Counter_t TranslateCnt;

Timer_t GravTm, PGravTm, GravSTm, GravMTm, GravQTm, GravHTm, GravQFTm, GravHFTm;
Timer_t MACTm;

#if 0
int64_t WatchId = 7310803995;
Key_t WatchKey = {.k = {781836800001, 0}};
#include <stdio.h>
/* #define Msg_do printf */
#define DebugWatchId(a, ...) \
   if (bp->ident == WatchId) \
       Msg_do(a, __VA_ARGS__)
#define DebugWatchKey(a, ...) \
   if (KeyContained(sink->key, WatchKey, NDIM)) \
       Msg_do(a, __VA_ARGS__)
#else
#define DebugWatchId(a, ...)
#define DebugWatchKey(a, ...)
int WatchId = -1;
Key_t WatchKey = {.k = {-1, -1}};
#endif
int Nwarn;

typedef float v4sf __attribute__ ((vector_size (16)));

#define BACKGROUND_FLAG (1<<10)
#define OFFSET_MASK (BACKGROUND_FLAG-1) /* needs to hold MAX_IMAGE */
#define MAX_IMAGE 125
#define offset_index(f) ((f) & OFFSET_MASK)

static v4sf Eps2v;
static int64_t GNobj, Nobj;
static float Eps2, Eps;
static int Smooth_type;
static float GNewt;
static const mac_s *mac;
static int Nimage = 1;
static float offset_array[MAX_IMAGE][NDIM];
static tree_t *SinkTree;
static body *Btab;

static void mxn_hexa(Sink *to, hcell *pp);
static int MxN_hblock = 4*1024;
static int MxN_min_sink = 256;
static int MxN_min_hsrc = 512;
static int MxN_do_pH = 0;

static grav_f Sinteract;
static grav_f Minteract;
static grav_f Qinteract;
static grav_f Hinteract;
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
static struct {
    v4sf mass, x, y, z;
} Svec[SVECSZ];

#define MVECSZ (180073) /* This should be dynamically extensible */
static struct {
    v4sf mass, x, y, z;
} Mvec[MVECSZ];

#ifdef QUAD
#define QVECSZ (7919)
static struct Qvec {
    v4sf mass, x, y, z;
    v4sf R;
#ifdef DIPOLE
    v4sf qx, qy, qz;
#endif
    v4sf qxx, qxy, qyy, qxz, qyz;
} Qvec[QVECSZ];
#endif

#ifdef HEXA
#define HVECSZ (24576*2)
static struct Hvec {
    v4sf mass, x, y, z;
    v4sf R;
#ifdef DIPOLE
    v4sf qx, qy, qz;
#endif
    v4sf qxx, qxy, qyy, qxz, qyz;
    v4sf qxxx, qxxy, qxyy, qyyy, qxxz, qxyz, qyyz;
    v4sf qxxxx, qxxxy, qxxyy, qxyyy, qyyyy, qxxxz, qxxyz, qxyyz, qyyyz;
} Hvec[HVECSZ];
#endif

static struct ucell_s {
    float halfsz, mass, x4, x2y2;
} ucell[CHUBITS];

void
SetupGrav(float newton_const, float e, int64_t gnobj, mac_s *m,
	  float particle_mass, int smooth_type)
{
    Nwarn = 0;
    GNewt = newton_const;
    Smooth_type = smooth_type;
    GNobj = gnobj;
    mac = m;
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
	if (mac->geometric_center) {
	    Qinteract = Arch(do_gravdq);
	    Hinteract = (amd6100) ? Arch(do_gravdh_amd6100) : Arch(do_gravdh);
	} else {
	    Qinteract = Arch(do_gravq);
	    Hinteract = (amd6100) ? Arch(do_gravh_amd6100) : Arch(do_gravh);
	}
    }
    Eps2 = Eps*Eps*pow(particle_mass, (float)(2./3.));
    Eps2v = (v4sf){Eps2, Eps2, Eps2, Eps2};
    float a1 = mac->r0;
    float a3 = mac->rho0*pow3(2.0f*mac->r0);
    float a5 = a3*pow2(2.0f*mac->r0);
    float a7 = a5*pow2(2.0f*mac->r0);
    for (int i = 0; i < CHUBITS; i++) {
	ucell[i].halfsz = a1;
	ucell[i].mass = a3;
	ucell[i].x4 = -a7/300.;
	ucell[i].x2y2 = a7/600.;
	a1 *= 0.5f;
	a3 *= pow3(0.5f);
	a5 *= pow3(0.5f)*pow2(0.5f);
	a7 *= pow3(0.5f)*pow2(0.5f)*pow2(0.5f);
    }
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
    StkPushType(ostk, Nimage/2, int);
}

void
WalkInitSrcPeriodic(Stk *kstk, Stk *ostk)
{
    int i;

    for (i = 0; i < Nimage; i++) {
	StkPushType(kstk, KeyInt(1), Key_t);
	StkPushType(ostk, i, int);
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
	float acc[NDIM], phi;
	float accd[NDIM], phid;

	float e;

	DebugWatchId("   %12g %12g %12g %ld %ld\n", bp->pos[0], bp->pos[1], bp->pos[2], 
		     pp->key.k[0], pp->key.k[1]);
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
	    Hvec[nn/NSSE].x[nn%NSSE] = 0.0f;
	    Hvec[nn/NSSE].y[nn%NSSE] = 0.0f;
	    Hvec[nn/NSSE].z[nn%NSSE] = 0.0f;
	    Hvec[nn/NSSE].R[nn%NSSE] = 0.0f;
#ifdef DIPOLE
	    Hvec[nn/NSSE].qx[nn%NSSE] = 0.0f;
	    Hvec[nn/NSSE].qy[nn%NSSE] = 0.0f;
	    Hvec[nn/NSSE].qz[nn%NSSE] = 0.0f;
#endif
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
	VS(accd, = 0.0); phid = 0.0;
	if (nn) Hinteract((float *)&Hvec[from->hcnt_done/NSSE], (float *)&Hvec[nn/NSSE], 
			  from->pos, &mtot, accd, &phid, &e, &ijunk);
	AddCounter(&BC4Int, from->hcnt-from->hcnt_done);
	StopTimer(&GravHTm);
	VV(acc, += accd); phi += phid;
	DebugWatchId("p4 %12g %12g %12g %d\n", acc[0], acc[1], acc[2], from->hcnt);
#endif
#ifdef QUAD
	StartTimer(&GravQTm);
	nn = from->qcnt;
	while (nn % NSSE) {
	    Qvec[nn/NSSE].mass[nn%NSSE] = 0.0f;
	    Qvec[nn/NSSE].x[nn%NSSE] = 0.0f;
	    Qvec[nn/NSSE].y[nn%NSSE] = 0.0f;
	    Qvec[nn/NSSE].z[nn%NSSE] = 0.0f;
	    Qvec[nn/NSSE].R[nn%NSSE] = 0.0f;
#ifdef DIPOLE
	    Qvec[nn/NSSE].qx[nn%NSSE] = 0.0f;
	    Qvec[nn/NSSE].qy[nn%NSSE] = 0.0f;
	    Qvec[nn/NSSE].qz[nn%NSSE] = 0.0f;
#endif
	    Qvec[nn/NSSE].qxx[nn%NSSE] = 0.0f;
	    Qvec[nn/NSSE].qxy[nn%NSSE] = 0.0f;
	    Qvec[nn/NSSE].qyy[nn%NSSE] = 0.0f;
	    Qvec[nn/NSSE].qxz[nn%NSSE] = 0.0f;
	    Qvec[nn/NSSE].qyz[nn%NSSE] = 0.0f;
	    nn++;
	}
	if ((long long)&Qvec[0] & 0xF || (long long)&Qvec[1] & 0xF)
	  Error("Qvec not aligned for asm code\n");
	VS(accd, = 0.0); phid = 0.0;
	if (nn) Qinteract((float *)&Qvec[from->qcnt_done/NSSE], (float *)&Qvec[nn/NSSE], 
			  from->pos, &mtot, accd, &phid, &e, &ijunk);
	DebugWatchId("p2 %12g %12g %12g %d\n", accd[0], accd[1], accd[2], from->qcnt);
	AddCounter(&BC2Int, from->qcnt-from->qcnt_done);
	StopTimer(&GravQTm);
	VV(acc, += accd); phi += phid;
#endif
	StartTimer(&GravMTm);
	nn = from->mcnt;
	while (nn % NSSE) {
	    Mvec[nn/NSSE].mass[nn%NSSE] = 0.0f;
	    Mvec[nn/NSSE].x[nn%NSSE] = 0.0f;
	    Mvec[nn/NSSE].y[nn%NSSE] = 0.0f;
	    Mvec[nn/NSSE].z[nn%NSSE] = 0.0f;
	    nn++;
	}
	if ((long long)&Mvec[0] & 0xF || (long long)&Mvec[1] & 0xF)
	  Error("Mvec not aligned for asm code\n");
	VS(accd, = 0.0); phid = 0.0;
	if (nn) Minteract((float *)&Mvec[0], (float *)&Mvec[nn/NSSE], 
			  from->pos, &mtot, accd, &phid, &e, &ijunk);
	DebugWatchId("p1 %12g %12g %12g %d\n", accd[0], accd[1], accd[2], from->mcnt);
	AddCounter(&BCInt, from->mcnt);
	StopTimer(&GravMTm);
	VV(acc, += accd); phi += phid;
	StartTimer(&GravSTm);
	nn = from->scnt;
	if (from->scnt > BSMax.counter) BSMax.counter = from->scnt;
	while (nn % NSSE) {
	    Svec[nn/NSSE].mass[nn%NSSE] = 0.0f;
	    Svec[nn/NSSE].x[nn%NSSE] = 0.0f;
	    Svec[nn/NSSE].y[nn%NSSE] = 0.0f;
	    Svec[nn/NSSE].z[nn%NSSE] = 0.0f;
	    nn++;
	}
	if ((long long)&Svec[0] & 0xF || (long long)&Svec[1] & 0xF)
	  Error("Svec not aligned for asm code\n");
	VS(accd, = 0.0); phid = 0.0;
	if (nn) Sinteract((float *)&Svec[0], (float *)&Svec[nn/NSSE], 
			  from->pos, &mtot, accd, &phid, &e, &ijunk);
	AddCounter(&BSInt, from->scnt);
	StopTimer(&GravSTm);
	StopTimer(&GravTm);
	VV(acc, += accd); phi += phid;
	if (!isfinite(acc[0]) || !isfinite(acc[1]) || !isfinite(acc[2]) || !isfinite(phi)) {
	    Error("bad results from do_grav for id %ld (%g,%g,%g), ax=%g ay=%g az=%g phi=%g\n", bp->ident, from->pos[0], from->pos[1], from->pos[2],
		  acc[0], acc[1], acc[2], phi);
	}
	
	if (Minteract == do_grav_sse16_ivec_asm) {
	    /* Fix self-interaction for phi from fast asm code */
	    /* If eps is very small, roundoff becomes a problem */
	    phi += bp->mass*recipsqrtf(e);
	}
	if (mac->subtract_background) {
	    float r[NDIM];
	    VVV(r, = from->pos, - from->cen);
	    VS(accd, = 0.0);
	    cubic_acc(r, from->cr, accd);
	    VS(accd, *= -mac->rho0);
	    VV(acc, += accd);
	    DebugWatchId("ca %12g - %g %g %g\n", from->cr, r[0], r[1], r[2]);
	    DebugWatchId("c %12g %12g %12g\n", accd[0], accd[1], accd[2]);
	    double cmass = -mac->rho0*pow3(2.0f*from->cr);
	    DebugWatchId("fmass %12g cmass %12g\n", from->fmass, cmass);
#if 0
	    if (fabs(from->fmass+cmass) > mac->m0*1e-7 || -cmass < bp->mass) 
		if (Nwarn++ < 2) SeriousWarning("Background subtraction off by %.2f for id %ld key %s {%lu,%lu} fmass %g cmass %g near %d\n", fabs(from->fmass+cmass)/bp->mass, bp->ident, PrintKey(pp->key), pp->key.k[0], pp->key.k[1], from->fmass/bp->mass, cmass/bp->mass, from->near);
	    if (from->near != 27) {
		Error("Bad Near for id %ld key %s {%lu,%lu} fmass %g cmass %g near %d\n", bp->ident, PrintKey(pp->key), pp->key.k[0], pp->key.k[1], from->fmass/bp->mass, cmass/bp->mass, from->near);
	    }
#endif
	}

	/* Make sure these are initialized to zero externally */
	bp->phi += from->M0;
	bp->phi += phi;
	bp->phi *= GNewt;
	VV(bp->acc, -= from->M1);
	VV(bp->acc, += acc);
	DebugWatchId("a %12g %12g %12g\n", bp->acc[0], bp->acc[1], bp->acc[2]);
	VS(bp->acc, *= GNewt);
	DebugWatchId("g %12g %12g %12g\n", bp->acc[0], bp->acc[1], bp->acc[2]);
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

    cell *ccp = NULL;
    if (Sub_Flags(pp)) {
	cell *cp = pp->ptr;

	VV(to->pos, = cp->pos);
	to->bmax = cp->bmax;
	to->daughters = cp->daughters;
	to->isbody = 0;
	if (cp->daughters <= 32 || 4.0f*mac->r0/mac->nx == ucell[cp->level].halfsz) ccp = cp;
    } else {
	body *bp = pp->ptr;
	VV(to->pos, = bp->pos);
	to->bmax = 0.0f;
	to->daughters = 1;
	to->isbody = 1;
    }

    if (from) {
	if (ccp && from->clevel == CHUBITS) {
	    to->clevel = ccp->level;
	    VV(to->cen, = ccp->pos);
	    to->cr = 3.0f*ucell[ccp->level].halfsz;
	    to->cr2 = pow2(3.99f*ucell[ccp->level].halfsz);
	} else {
	    to->clevel = from->clevel;
	    VV(to->cen, = from->cen);
	    to->cr = from->cr;
	    to->cr2 = from->cr2;
	} 
	if (to->isbody && from->clevel == CHUBITS) {
	    float pos[NDIM], cellsz;
	    CellCorner(pp->key, pos, &cellsz);
	    VS(pos, += 0.5f*cellsz);
	    to->clevel = TreeLevel(pp->key, NDIM);
	    VV(to->cen, = pos);
	    to->cr = 1.5f*cellsz;
	    to->cr2 = pow2(1.995f*cellsz);
	}
	to->near = from->near;
	to->fmass = from->fmass;
	if ((TreeLevel(pp->key, NDIM) >= 1) && KeyContained(pp->key, WatchKey, NDIM)) 
	    Msg_do("Inherit to %s, cr %g\n", PrintKey(pp->key), to->cr);
	to->interactions = from->interactions;
	to->key = pp->key;
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
	to->fmass = 0.0;
	to->near = 0;
	to->clevel = CHUBITS;
	VS(to->cen, = 0.0f);
	to->cr = 0.0f;
	to->cr2 = 0.0f;
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

    if (!isfinite(first->acc[0]) || !isfinite(first->acc[1]) || !isfinite(first->acc[2]) || !isfinite(first->phi)) {
	Error("bad results from do_grav for %ld (%g,%g,%g), ax=%g ay=%g az=%g phi=%g\n", 
	      first->ident, first->pos[0], first->pos[1], first->pos[2],
	      first->acc[0], first->acc[1], first->acc[2], first->phi);
    }
    StopTimer(&GravHFTm);
    StopTimer(&GravTm);
}

void
RcritMAC(Sink *sink, const hcell **source_vec, int *flags, int *result, int n)
{
#if 0
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
	    VxVV(r, = bp->pos, + offset_array[offset_index(flags[i])]);
	    VxVxVx(dx, = r, - pos_sink);
	    dr2 = Dotx(dx, dx);
	    if (dr2 > Eps2) {
		Mvec[mcnt/NSSE].mass[mcnt%NSSE] = bp->mass;
		Mvec[mcnt/NSSE].x[mcnt%NSSE] = r0;
		Mvec[mcnt/NSSE].y[mcnt%NSSE] = r1;
		Mvec[mcnt/NSSE].z[mcnt%NSSE] = r2;
		mcnt++;
		if (mcnt/NSSE >= MVECSZ) Error("mvec overflow\n");
	    } else if (dr2 > 0.0f) {
		Svec[scnt/NSSE].mass[scnt%NSSE] = bp->mass;
		Svec[scnt/NSSE].x[scnt%NSSE] = r0;
		Svec[scnt/NSSE].y[scnt%NSSE] = r1;
		Svec[scnt/NSSE].z[scnt%NSSE] = r2;
		scnt++;
		if (scnt/NSSE >= SVECSZ) Error("svec overflow\n");
	    }
	    interactions++;
	    result[i] = MAC_ACCEPT;
	    continue;
	}

	VxVV(r, = cp->pos, + offset_array[offset_index(flags[i])]);
	VxVxVx(dx, = r, - pos_sink);
	dr2 = Dotx(dx, dx);

	if (dr2 >= cp->rcrit*cp->rcrit) {
	    Mvec[mcnt/NSSE].mass[mcnt%NSSE] = cp->mass;
	    Mvec[mcnt/NSSE].x[mcnt%NSSE] = r0;
	    Mvec[mcnt/NSSE].y[mcnt%NSSE] = r1;
	    Mvec[mcnt/NSSE].z[mcnt%NSSE] = r2;
	    mcnt++;
	    if (mcnt/NSSE >= MVECSZ) Error("mvec overflow\n");
	    interactions += cp->daughters;
	    result[i] = MAC_ACCEPT;
#ifdef QUAD
	} else if (Quad_Ncut && (cp->daughters >= Quad_Ncut) &&
		   (dr2 > qcp->rcrit_q*qcp->rcrit_q)) {
	    Qvec[qcnt/NSSE].mass[qcnt%NSSE] = qcp->mass;
	    Qvec[qcnt/NSSE].x[qcnt%NSSE] = r0;
	    Qvec[qcnt/NSSE].y[qcnt%NSSE] = r1;
	    Qvec[qcnt/NSSE].z[qcnt%NSSE] = r2;
	    Qvec[qcnt/NSSE].R[qcnt%NSSE] = qcp->bmax;
#ifdef DIPOLE
	    Qvec[qcnt/NSSE].qx[qcnt%NSSE] = qcp->qx;
	    Qvec[qcnt/NSSE].qy[qcnt%NSSE] = qcp->qy;
	    Qvec[qcnt/NSSE].qz[qcnt%NSSE] = qcp->qz;
#endif
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
	    Hvec[hcnt/NSSE].x[hcnt%NSSE] = r0;
	    Hvec[hcnt/NSSE].y[hcnt%NSSE] = r1;
	    Hvec[hcnt/NSSE].z[hcnt%NSSE] = r2;
	    Hvec[hcnt/NSSE].R[hcnt%NSSE] = hcp->bmax;
#ifdef DIPOLE
	    Hvec[hcnt/NSSE].qx[hcnt%NSSE] = hcp->qx;
	    Hvec[hcnt/NSSE].qy[hcnt%NSSE] = hcp->qy;
	    Hvec[hcnt/NSSE].qz[hcnt%NSSE] = hcp->qz;
#endif
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
#endif
}

#if 0

/* Transpose the 4x4 matrix composed of row[0-3].  */
#define _MM_TRANSPOSE4_PS(row0, row1, row2, row3)			\
do {									\
  v4sf __r0 = (row0), __r1 = (row1), __r2 = (row2), __r3 = (row3);	\
  v4sf __t0 = __builtin_ia32_unpcklps (__r0, __r1);			\
  v4sf __t1 = __builtin_ia32_unpcklps (__r2, __r3);			\
  v4sf __t2 = __builtin_ia32_unpckhps (__r0, __r1);			\
  v4sf __t3 = __builtin_ia32_unpckhps (__r2, __r3);			\
  (row0) = __builtin_ia32_movlhps (__t0, __t1);				\
  (row1) = __builtin_ia32_movhlps (__t1, __t0);				\
  (row2) = __builtin_ia32_movlhps (__t2, __t3);				\
  (row3) = __builtin_ia32_movhlps (__t3, __t2);				\
} while (0)


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
		Mvec[mcnt/NSSE].x[mcnt%NSSE] = r0;
		Mvec[mcnt/NSSE].y[mcnt%NSSE] = r1;
		Mvec[mcnt/NSSE].z[mcnt%NSSE] = r2;
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
		    Qvec[qcnt/NSSE].x[qcnt%NSSE] = r0;
		    Qvec[qcnt/NSSE].y[qcnt%NSSE] = r1;
		    Qvec[qcnt/NSSE].z[qcnt%NSSE] = r2;
		    Qvec[qcnt/NSSE].R[qcnt%NSSE] = qcp->bmax;
#ifdef DIPOLE
		    Qvec[qcnt/NSSE].qx[qcnt%NSSE] = qcp->qx;
		    Qvec[qcnt/NSSE].qy[qcnt%NSSE] = qcp->qy;
		    Qvec[qcnt/NSSE].qz[qcnt%NSSE] = qcp->qz;
#endif
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
		    Hvec[hcnt/NSSE].x[hcnt%NSSE] = r0;
		    Hvec[hcnt/NSSE].y[hcnt%NSSE] = r1;
		    Hvec[hcnt/NSSE].z[hcnt%NSSE] = r2;
		    Hvec[hcnt/NSSE].R[hcnt%NSSE] = hcp->bmax;
#ifdef DIPOLE
		    Hvec[hcnt/NSSE].qx[hcnt%NSSE] = hcp->qx;
		    Hvec[hcnt/NSSE].qy[hcnt%NSSE] = hcp->qy;
		    Hvec[hcnt/NSSE].qz[hcnt%NSSE] = hcp->qz;
#endif
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
	    if ((bmax > mac->dlmax) || (mac->dlfac * bmax > rcrit)) {
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
		Mvec[mcnt/NSSE].x[mcnt%NSSE] = r0;
		Mvec[mcnt/NSSE].y[mcnt%NSSE] = r1;
		Mvec[mcnt/NSSE].z[mcnt%NSSE] = r2;
		mcnt++;
		if (mcnt/NSSE >= MVECSZ) Error("mvec overflow\n");
	    } else if (dr2 > 0.0f) {
		Svec[scnt/NSSE].mass[scnt%NSSE] = bp->mass;
		Svec[scnt/NSSE].x[scnt%NSSE] = r0;
		Svec[scnt/NSSE].y[scnt%NSSE] = r1;
		Svec[scnt/NSSE].z[scnt%NSSE] = r2;
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

/* RcritMAC with Don't Laugh-like traversal */
void
DLRcritMAC(Sink *sink, const hcell **source_vec, const float **offset_vec, int *result, int n)
{
    float dr2;
    Vxd(float r);
    Vxd(float dx);
    int nh = 0, nq = 0, nm = 0, ns = 0;
    const hexacell *hptr[n+NSSE]; /* could use struct, but would have padding */
    const quadcell *qptr[n+NSSE];
    v4sf hbuf[n+NSSE];
    v4sf qbuf[n+NSSE];
    v4sf mbuf[n+NSSE];
    v4sf sbuf[n+NSSE];

    StartTimer(&MACTm);
    for (int i = 0; i < n; i++) {
	if (Sub_Flags(source_vec[i])) {
	    const cell *cp = source_vec[i]->ptr;
	    const quadcell *qcp = source_vec[i]->ptr;
	    const hexacell *hcp = source_vec[i]->ptr;
	    VxVV(r, = cp->pos, + offset_vec[i]);
	    VxVxV(dx, = r, - sink->pos);
	    dr2 = Dotx(dx, dx);

	    int isquad = Quad_Ncut && (cp->daughters >= Quad_Ncut);
	    int ishexa = Hexa_Ncut && (cp->daughters >= Hexa_Ncut);
	    float smallest_rcrit = ishexa ? hcp->rcrit_h : (isquad ? qcp->rcrit_q : cp->rcrit);
	    
	    /* bmax is 0 if sink is a body */
	    if (dr2 > Square(cp->rcrit + sink->bmax)) {
		/* cell-cell or body-cell */
		mbuf[nm++] = (v4sf){cp->mass, r0, r1, r2};
		sink->interactions += cp->daughters;
		result[i] = MAC_ACCEPT;
	    } else if (isquad && dr2 > Square(qcp->rcrit_q + sink->bmax)) {
		qptr[nq] = qcp;
		qbuf[nq++] = (v4sf){cp->mass, r0, r1, r2};
		sink->interactions += cp->daughters;
		result[i] = MAC_ACCEPT;
	    } else if (ishexa && dr2 > Square(hcp->rcrit_h + sink->bmax)) {
		hptr[nh] = hcp;
		hbuf[nh++] = (v4sf){cp->mass, r0, r1, r2};
		sink->interactions += cp->daughters;
		result[i] = MAC_ACCEPT;
	    } else if ((sink->bmax > mac->dlmax) || (mac->dlfac * sink->bmax > smallest_rcrit)) {
		result[i] = sink->isbody ? MAC_ERROR : MAC_SPLIT_SINK;
	    } else {
		result[i] = MAC_SPLIT_SRC;
	    }
	} else if (sink->isbody) {
	    /* body-body */
	    const body *bp = source_vec[i]->ptr;
	    VxVV(r, = bp->pos, + offset_vec[i]);
	    VxVxV(dx, = r, - sink->pos);
	    dr2 = Dotx(dx, dx);

	    if (dr2 > Eps2) mbuf[nm++] = (v4sf){bp->mass, r0, r1, r2};
	    else if (dr2 > 0.0f) sbuf[ns++] = (v4sf){bp->mass, r0, r1, r2};
	    sink->interactions++;
	    result[i] = MAC_ACCEPT;
	} else {
	    /* cell-body */
	    /* Comparing dr2 with Eps2 only works for body-body interactions */
	    result[i] = MAC_SPLIT_SINK;
	}
    }

    hexacell hzero = {};
    quadcell qzero = {};
    int i, j, k, m;

    j = sink->hcnt / NSSE;
    k = sink->hcnt % NSSE;
    m = k ? Min(NSSE-k, nh) : 0;
    sink->hcnt += nh;
    if (sink->hcnt/NSSE >= HVECSZ) Error("hvec overflow\n"); 
    if (m) {
	for (i = 0; i < m; i++) Hvec[j].mass[i+k] = hbuf[i][0];
	for (i = 0; i < m; i++) Hvec[j].x[i+k] = hbuf[i][1];
	for (i = 0; i < m; i++) Hvec[j].y[i+k] = hbuf[i][2];
	for (i = 0; i < m; i++) Hvec[j].z[i+k] = hbuf[i][3];
	for (i = 0; i < m; i++) Hvec[j].R[i+k] = hptr[i]->bmax;
	for (i = 0; i < m; i++) Hvec[j].qx[i+k] = hptr[i]->qx;
	for (i = 0; i < m; i++) Hvec[j].qy[i+k] = hptr[i]->qy;
	for (i = 0; i < m; i++) Hvec[j].qz[i+k] = hptr[i]->qz;
	for (i = 0; i < m; i++) Hvec[j].qxx[i+k] = hptr[i]->qxx;
	for (i = 0; i < m; i++) Hvec[j].qxy[i+k] = hptr[i]->qxy;
	for (i = 0; i < m; i++) Hvec[j].qyy[i+k] = hptr[i]->qyy;
	for (i = 0; i < m; i++) Hvec[j].qxz[i+k] = hptr[i]->qxz;
	for (i = 0; i < m; i++) Hvec[j].qyz[i+k] = hptr[i]->qyz;
	for (i = 0; i < m; i++) Hvec[j].qxxx[i+k] = hptr[i]->qxxx;
	for (i = 0; i < m; i++) Hvec[j].qxxy[i+k] = hptr[i]->qxxy;
	for (i = 0; i < m; i++) Hvec[j].qxyy[i+k] = hptr[i]->qxyy;
	for (i = 0; i < m; i++) Hvec[j].qyyy[i+k] = hptr[i]->qyyy;
	for (i = 0; i < m; i++) Hvec[j].qxxz[i+k] = hptr[i]->qxxz;
	for (i = 0; i < m; i++) Hvec[j].qxyz[i+k] = hptr[i]->qxyz;
	for (i = 0; i < m; i++) Hvec[j].qyyz[i+k] = hptr[i]->qyyz;
	for (i = 0; i < m; i++) Hvec[j].qxxxx[i+k] = hptr[i]->qxxxx;
	for (i = 0; i < m; i++) Hvec[j].qxxxy[i+k] = hptr[i]->qxxxy;
	for (i = 0; i < m; i++) Hvec[j].qxxyy[i+k] = hptr[i]->qxxyy;
	for (i = 0; i < m; i++) Hvec[j].qxyyy[i+k] = hptr[i]->qxyyy;
	for (i = 0; i < m; i++) Hvec[j].qyyyy[i+k] = hptr[i]->qyyyy;
	for (i = 0; i < m; i++) Hvec[j].qxxxz[i+k] = hptr[i]->qxxxz;
	for (i = 0; i < m; i++) Hvec[j].qxxyz[i+k] = hptr[i]->qxxyz;
	for (i = 0; i < m; i++) Hvec[j].qxyyz[i+k] = hptr[i]->qxyyz;
	for (i = 0; i < m; i++) Hvec[j].qyyyz[i+k] = hptr[i]->qyyyz;
	j++;
    }
    for (i = nh; i < nh + NSSE; i++) {
	hptr[i] = &hzero;
	hbuf[i] = (v4sf){0.0f, 0.0f, 0.0f, 0.0f};
    }
    for (i = m; i < nh; i += NSSE, j++) {
	v4sf r0 = hbuf[i+0];
	v4sf r1 = hbuf[i+1];
	v4sf r2 = hbuf[i+2];
	v4sf r3 = hbuf[i+3];
	_MM_TRANSPOSE4_PS(r0, r1, r2, r3);
	Hvec[j].mass = r0;
	Hvec[j].x = r1;
	Hvec[j].y = r2;
	Hvec[j].z = r3;
	r0 = __builtin_ia32_loadups(&hptr[i+0]->bmax);
	r1 = __builtin_ia32_loadups(&hptr[i+1]->bmax);
	r2 = __builtin_ia32_loadups(&hptr[i+2]->bmax);
	r3 = __builtin_ia32_loadups(&hptr[i+3]->bmax);
	_MM_TRANSPOSE4_PS(r0, r1, r2, r3);
	Hvec[j].R = r0;
	Hvec[j].qx = r1;
	Hvec[j].qy = r2;
	Hvec[j].qz = r3;
	r0 = __builtin_ia32_loadups(&hptr[i+0]->qxx);
	r1 = __builtin_ia32_loadups(&hptr[i+1]->qxx);
	r2 = __builtin_ia32_loadups(&hptr[i+2]->qxx);
	r3 = __builtin_ia32_loadups(&hptr[i+3]->qxx);
	_MM_TRANSPOSE4_PS(r0, r1, r2, r3);
	Hvec[j].qxx = r0;
	Hvec[j].qxy = r1;
	Hvec[j].qyy = r2;
	Hvec[j].qxz = r3;
	r0 = __builtin_ia32_loadups(&hptr[i+0]->qyz);
	r1 = __builtin_ia32_loadups(&hptr[i+1]->qyz);
	r2 = __builtin_ia32_loadups(&hptr[i+2]->qyz);
	r3 = __builtin_ia32_loadups(&hptr[i+3]->qyz);
	_MM_TRANSPOSE4_PS(r0, r1, r2, r3);
	Hvec[j].qyz = r0;
	Hvec[j].qxxx = r1;
	Hvec[j].qxxy = r2;
	Hvec[j].qxyy = r3;
	r0 = __builtin_ia32_loadups(&hptr[i+0]->qyyy);
	r1 = __builtin_ia32_loadups(&hptr[i+1]->qyyy);
	r2 = __builtin_ia32_loadups(&hptr[i+2]->qyyy);
	r3 = __builtin_ia32_loadups(&hptr[i+3]->qyyy);
	_MM_TRANSPOSE4_PS(r0, r1, r2, r3);
	Hvec[j].qyyy = r0;
	Hvec[j].qxxz = r1;
	Hvec[j].qxyz = r2;
	Hvec[j].qyyz = r3;
	r0 = __builtin_ia32_loadups(&hptr[i+0]->qxxxx);
	r1 = __builtin_ia32_loadups(&hptr[i+1]->qxxxx);
	r2 = __builtin_ia32_loadups(&hptr[i+2]->qxxxx);
	r3 = __builtin_ia32_loadups(&hptr[i+3]->qxxxx);
	_MM_TRANSPOSE4_PS(r0, r1, r2, r3);
	Hvec[j].qxxxx = r0;
	Hvec[j].qxxxy = r1;
	Hvec[j].qxxyy = r2;
	Hvec[j].qxyyy = r3;
	r0 = __builtin_ia32_loadups(&hptr[i+0]->qyyyy);
	r1 = __builtin_ia32_loadups(&hptr[i+1]->qyyyy);
	r2 = __builtin_ia32_loadups(&hptr[i+2]->qyyyy);
	r3 = __builtin_ia32_loadups(&hptr[i+3]->qyyyy);
	_MM_TRANSPOSE4_PS(r0, r1, r2, r3);
	Hvec[j].qyyyy = r0;
	Hvec[j].qxxxz = r1;
	Hvec[j].qxxyz = r2;
	Hvec[j].qxyyz = r3;
	Hvec[j].qyyyz = (v4sf){hptr[i]->qyyyz, hptr[i+1]->qyyyz, hptr[i+2]->qyyyz, hptr[i+3]->qyyyz};
    }

    j = sink->qcnt / NSSE;
    k = sink->qcnt % NSSE;
    m = k ? Min(NSSE-k, nq) : 0;
    sink->qcnt += nq;
    if (sink->qcnt/NSSE >= QVECSZ) Error("qvec overflow\n"); 
    if (m) {
	for (i = 0; i < m; i++) Qvec[j].mass[i+k] = qbuf[i][0];
	for (i = 0; i < m; i++) Qvec[j].x[i+k] = qbuf[i][1];
	for (i = 0; i < m; i++) Qvec[j].y[i+k] = qbuf[i][2];
	for (i = 0; i < m; i++) Qvec[j].z[i+k] = qbuf[i][3];
	for (i = 0; i < m; i++) Qvec[j].R[i+k] = qptr[i]->bmax;
	for (i = 0; i < m; i++) Qvec[j].qx[i+k] = qptr[i]->qx;
	for (i = 0; i < m; i++) Qvec[j].qy[i+k] = qptr[i]->qy;
	for (i = 0; i < m; i++) Qvec[j].qz[i+k] = qptr[i]->qz;
	for (i = 0; i < m; i++) Qvec[j].qxx[i+k] = qptr[i]->qxx;
	for (i = 0; i < m; i++) Qvec[j].qxy[i+k] = qptr[i]->qxy;
	for (i = 0; i < m; i++) Qvec[j].qyy[i+k] = qptr[i]->qyy;
	for (i = 0; i < m; i++) Qvec[j].qxz[i+k] = qptr[i]->qxz;
	for (i = 0; i < m; i++) Qvec[j].qyz[i+k] = qptr[i]->qyz;
	j++;
    }
    for (i = nq; i < nq + NSSE; i++) {
	qptr[i] = &qzero;
	qbuf[i] = (v4sf){0.0f, 0.0f, 0.0f, 0.0f};
    }
    for (i = m; i < nq; i += NSSE, j++) {
	v4sf r0 = qbuf[i+0];
	v4sf r1 = qbuf[i+1];
	v4sf r2 = qbuf[i+2];
	v4sf r3 = qbuf[i+3];
	_MM_TRANSPOSE4_PS(r0, r1, r2, r3);
	Qvec[j].mass = r0;
	Qvec[j].x = r1;
	Qvec[j].y = r2;
	Qvec[j].z = r3;
	r0 = __builtin_ia32_loadups(&qptr[i+0]->bmax);
	r1 = __builtin_ia32_loadups(&qptr[i+1]->bmax);
	r2 = __builtin_ia32_loadups(&qptr[i+2]->bmax);
	r3 = __builtin_ia32_loadups(&qptr[i+3]->bmax);
	_MM_TRANSPOSE4_PS(r0, r1, r2, r3);
	Qvec[j].R = r0;
	Qvec[j].qx = r1;
	Qvec[j].qy = r2;
	Qvec[j].qz = r3;
	r0 = __builtin_ia32_loadups(&qptr[i+0]->qxx);
	r1 = __builtin_ia32_loadups(&qptr[i+1]->qxx);
	r2 = __builtin_ia32_loadups(&qptr[i+2]->qxx);
	r3 = __builtin_ia32_loadups(&qptr[i+3]->qxx);
	_MM_TRANSPOSE4_PS(r0, r1, r2, r3);
	Qvec[j].qxx = r0;
	Qvec[j].qxy = r1;
	Qvec[j].qyy = r2;
	Qvec[j].qxz = r3;
	Qvec[j].qyz = (v4sf){qptr[i]->qyz, qptr[i+1]->qyz, qptr[i+2]->qyz, qptr[i+3]->qyz};
    }

    j = sink->mcnt / NSSE;
    k = sink->mcnt % NSSE;
    m = k ? Min(NSSE-k, nm) : 0;
    sink->mcnt += nm;
    if (sink->mcnt/NSSE >= MVECSZ) Error("mvec overflow\n"); 
    if (m) {
	for (i = 0; i < m; i++) Mvec[j].mass[i+k] = mbuf[i][0];
	for (i = 0; i < m; i++) Mvec[j].x[i+k] = mbuf[i][1];
	for (i = 0; i < m; i++) Mvec[j].y[i+k] = mbuf[i][2];
	for (i = 0; i < m; i++) Mvec[j].z[i+k] = mbuf[i][3];
	j++;
    }
    for (i = nm; i < nm + NSSE; i++) {
	mbuf[i] = (v4sf){0.0f, 0.0f, 0.0f, 0.0f};
    }
    for (i = m; i < nm; i += NSSE, j++) {
	v4sf r0 = mbuf[i+0];
	v4sf r1 = mbuf[i+1];
	v4sf r2 = mbuf[i+2];
	v4sf r3 = mbuf[i+3];
	_MM_TRANSPOSE4_PS(r0, r1, r2, r3);
	Mvec[j].mass = r0;
	Mvec[j].x = r1;
	Mvec[j].y = r2;
	Mvec[j].z = r3;
    }

    j = sink->scnt / NSSE;
    k = sink->scnt % NSSE;
    m = k ? Min(NSSE-k, ns) : 0;
    sink->scnt += ns;
    if (sink->scnt/NSSE >= SVECSZ) Error("svec overflow\n"); 
    if (m) {
	for (i = 0; i < m; i++) Svec[j].mass[i+k] = sbuf[i][0];
	for (i = 0; i < m; i++) Svec[j].x[i+k] = sbuf[i][1];
	for (i = 0; i < m; i++) Svec[j].y[i+k] = sbuf[i][2];
	for (i = 0; i < m; i++) Svec[j].z[i+k] = sbuf[i][3];
	j++;
    }
    for (i = ns; i < ns + NSSE; i++) {
	sbuf[i] = (v4sf){0.0f, 0.0f, 0.0f, 0.0f};
    }
    for (i = m; i < ns; i += NSSE, j++) {
	v4sf r0 = sbuf[i+0];
	v4sf r1 = sbuf[i+1];
	v4sf r2 = sbuf[i+2];
	v4sf r3 = sbuf[i+3];
	_MM_TRANSPOSE4_PS(r0, r1, r2, r3);
	Svec[j].mass = r0;
	Svec[j].x = r1;
	Svec[j].y = r2;
	Svec[j].z = r3;
    }
    StopTimer(&MACTm);
}

#define appendMvec(mxyz, p) \
    do { \
      int _i = sink->mcnt/NSSE; \
      int _j = sink->mcnt%NSSE;	 \
      Mvec[_i].mass[_j] = mxyz[0]; \
      Mvec[_i].x[_j] = mxyz[1]; \
      Mvec[_i].y[_j] = mxyz[2]; \
      Mvec[_i].z[_j] = mxyz[3]; \
      sink->mcnt++;		\
    } while(0)

#define appendSvec(mxyz, p) \
    do { \
      int _i = sink->scnt/NSSE; \
      int _j = sink->scnt%NSSE;	 \
      Svec[_i].mass[_j] = mxyz[0]; \
      Svec[_i].x[_j] = mxyz[1]; \
      Svec[_i].y[_j] = mxyz[2]; \
      Svec[_i].z[_j] = mxyz[3]; \
      sink->scnt++;		\
    } while(0)

#define appendQvec(mxyz, p) \
    do { \
      int _i = sink->qcnt/NSSE; \
      int _j = sink->qcnt%NSSE;	 \
      Qvec[_i].mass[_j] = mxyz[0]; \
      Qvec[_i].x[_j] = mxyz[1]; \
      Qvec[_i].y[_j] = mxyz[2]; \
      Qvec[_i].z[_j] = mxyz[3]; \
      Qvec[_i].R[_j] = (p)->bmax; \
      Qvec[_i].qx[_j] = (p)->qx; \
      Qvec[_i].qy[_j] = (p)->qy; \
      Qvec[_i].qz[_j] = (p)->qz; \
      Qvec[_i].qxx[_j] = (p)->qxx; \
      Qvec[_i].qxy[_j] = (p)->qxy; \
      Qvec[_i].qyy[_j] = (p)->qyy; \
      Qvec[_i].qxz[_j] = (p)->qxz; \
      Qvec[_i].qyz[_j] = (p)->qyz; \
      sink->qcnt++;		\
    } while(0)

#define appendHvec(mxyz, p) \
    do { \
      int _i = sink->hcnt/NSSE; \
      int _j = sink->hcnt%NSSE;	 \
      Hvec[_i].mass[_j] = mxyz[0]; \
      Hvec[_i].x[_j] = mxyz[1]; \
      Hvec[_i].y[_j] = mxyz[2]; \
      Hvec[_i].z[_j] = mxyz[3]; \
      Hvec[_i].R[_j] = (p)->bmax; \
      Hvec[_i].qx[_j] = (p)->qx; \
      Hvec[_i].qy[_j] = (p)->qy; \
      Hvec[_i].qz[_j] = (p)->qz; \
      Hvec[_i].qxx[_j] = (p)->qxx; \
      Hvec[_i].qxy[_j] = (p)->qxy; \
      Hvec[_i].qyy[_j] = (p)->qyy; \
      Hvec[_i].qxz[_j] = (p)->qxz; \
      Hvec[_i].qyz[_j] = (p)->qyz; \
      Hvec[_i].qxxx[_j] = (p)->qxxx; \
      Hvec[_i].qxxy[_j] = (p)->qxxy; \
      Hvec[_i].qxyy[_j] = (p)->qxyy; \
      Hvec[_i].qyyy[_j] = (p)->qyyy; \
      Hvec[_i].qxxz[_j] = (p)->qxxz; \
      Hvec[_i].qxyz[_j] = (p)->qxyz; \
      Hvec[_i].qyyz[_j] = (p)->qyyz; \
      Hvec[_i].qxxxx[_j] = (p)->qxxxx; \
      Hvec[_i].qxxxy[_j] = (p)->qxxxy; \
      Hvec[_i].qxxyy[_j] = (p)->qxxyy; \
      Hvec[_i].qxyyy[_j] = (p)->qxyyy; \
      Hvec[_i].qyyyy[_j] = (p)->qyyyy; \
      Hvec[_i].qxxxz[_j] = (p)->qxxxz; \
      Hvec[_i].qxxyz[_j] = (p)->qxxyz; \
      Hvec[_i].qxyyz[_j] = (p)->qxyyz; \
      Hvec[_i].qyyyz[_j] = (p)->qyyyz; \
      sink->hcnt++;		\
    } while(0)

/* RcritMAC with Don't Laugh-like traversal */
void
DLRcritMAC(Sink *sink, const hcell **source_vec, const float **offset_vec, int *result, int n)
{
    const v4sf sink_pos = {0.0f, sink->pos[0], sink->pos[1], sink->pos[2]};
    const v4sf sink_bmax = {sink->bmax, sink->bmax, sink->bmax, sink->bmax};
    const v4sf DLfac_bmax = {mac->dlfac*sink->bmax, mac->dlfac*sink->bmax, mac->dlfac*sink->bmax, mac->dlfac*sink->bmax};
    const v4sf zero = {0.0f, 0.0f, 0.0f, 0.0f};

    StartTimer(&MACTm);
    for (int i = 0; i < n; i += 4) {
	const cell *cp[4] = {source_vec[i+0]->ptr, 
			     source_vec[i+1 < n ? i+1 : i]->ptr,
			     source_vec[i+2 < n ? i+2 : i]->ptr,
			     source_vec[i+3 < n ? i+3 : i]->ptr};
	int iscell = 
	    (Sub_Flags(source_vec[i+0]) ? 1 : 0) | 
	    (Sub_Flags(source_vec[i+1 < n ? i+1 : i]) ? 2 : 0) | 
	    (Sub_Flags(source_vec[i+2 < n ? i+2 : i]) ? 4 : 0) | 
	    (Sub_Flags(source_vec[i+3 < n ? i+3 : i]) ? 8 : 0);
	v4sf r[4];
	r[0] = __builtin_ia32_loadups((const float *)cp[0]);
	r[1] = __builtin_ia32_loadups((const float *)cp[1]);
	r[2] = __builtin_ia32_loadups((const float *)cp[2]);
	r[3] = __builtin_ia32_loadups((const float *)cp[3]);
	v4sf rc = {iscell & 1 ? cp[0]->rcrit : 0.0f,
		   iscell & 2 ? cp[1]->rcrit : 0.0f,
		   iscell & 4 ? cp[2]->rcrit : 0.0f,
		   iscell & 8 ? cp[3]->rcrit : 0.0f};
	rc += sink_bmax;
	int splitsink = __builtin_ia32_movmskps(__builtin_ia32_cmpleps(rc, DLfac_bmax));
	v4sf rc2 = rc*rc;
	int hashexa =
	    (iscell & 1 && cp[0]->daughters >= Hexa_Ncut ? 1 : 0) | 
	    (iscell & 2 && cp[1]->daughters >= Hexa_Ncut ? 2 : 0) | 
	    (iscell & 4 && cp[2]->daughters >= Hexa_Ncut ? 4 : 0) | 
	    (iscell & 8 && cp[3]->daughters >= Hexa_Ncut ? 8 : 0);
	int hasquad =
	    (iscell & 1 && cp[0]->daughters >= Quad_Ncut ? 1 : 0) | 
	    (iscell & 2 && cp[1]->daughters >= Quad_Ncut ? 2 : 0) | 
	    (iscell & 4 && cp[2]->daughters >= Quad_Ncut ? 4 : 0) | 
	    (iscell & 8 && cp[3]->daughters >= Quad_Ncut ? 8 : 0);
	r[0] += *(v4sf *)(offset_vec[i+0]);
	r[1] += *(v4sf *)(offset_vec[i+1 < n ? i+1 : i]);
	r[2] += *(v4sf *)(offset_vec[i+2 < n ? i+2 : i]);
	r[3] += *(v4sf *)(offset_vec[i+3 < n ? i+3 : i]);
	v4sf dr0 = r[0]-sink_pos;
	v4sf dr1 = r[1]-sink_pos;
	v4sf dr2 = r[2]-sink_pos;
	v4sf dr3 = r[3]-sink_pos;
	_MM_TRANSPOSE4_PS(dr0, dr1, dr2, dr3);
	v4sf r2 = dr1*dr1;
	r2 += dr2*dr2;
	r2 += dr3*dr3;
	int mactrue = __builtin_ia32_movmskps(__builtin_ia32_cmpleps(rc2,r2));
	int notsmooth = __builtin_ia32_movmskps(__builtin_ia32_cmpleps(Eps2v,r2));
	int notself = __builtin_ia32_movmskps(__builtin_ia32_cmpneqps(r2,zero));

	for (int j = 0; j < 4 && i+j < n; j++) {
	    int b = 1<<j;
	    if (iscell & b) {
		if (mactrue & b) {
		    result[i+j] = MAC_ACCEPT;
		    sink->interactions += cp[j]->daughters;
		    if (hashexa & b) appendHvec(r[j], (hexacell *)cp[j]);
		    else if (hasquad & b) appendQvec(r[j], (quadcell *)cp[j]);
		    else appendMvec(r[j], cp[j]);
		} else result[i+j] = splitsink & b ? MAC_SPLIT_SINK : MAC_SPLIT_SRC;
	    } else if (sink->isbody) {
		result[i+j] = MAC_ACCEPT;
		sink->interactions++;
		if (notsmooth & b) appendMvec(r[j], cp[j]);
		else if (notself & b) appendSvec(r[j], cp[j]);
	    } else result[i+j] = MAC_SPLIT_SINK;
	} 
    }
    StopTimer(&MACTm);
    
    if (sink->scnt/NSSE >= SVECSZ) Error("svec overflow\n");
    if (sink->mcnt/NSSE >= MVECSZ) Error("mvec overflow\n");
    if (sink->qcnt/NSSE >= QVECSZ) Error("qvec overflow\n");
    if (sink->hcnt/NSSE >= HVECSZ) Error("hvec overflow\n");
}

#endif

#define appendMvec(p, m)				\
    do { \
      int _i = sink->mcnt/NSSE; \
      int _j = sink->mcnt%NSSE;	 \
      Mvec[_i].mass[_j] = m; \
      Mvec[_i].x[_j] = r0; \
      Mvec[_i].y[_j] = r1; \
      Mvec[_i].z[_j] = r2; \
      sink->mcnt++; \
    } while(0)

#define appendSvec(p) \
    do { \
      int _i = sink->scnt/NSSE; \
      int _j = sink->scnt%NSSE;	 \
      Svec[_i].mass[_j] = (p)->mass; \
      Svec[_i].x[_j] = r0; \
      Svec[_i].y[_j] = r1; \
      Svec[_i].z[_j] = r2; \
      sink->scnt++; \
    } while(0)

#define appendQdvec(p)							\
    float _m = bs ? (p)->mass + ucell[(p)->level].mass : (p)->mass;	\
    do {								\
	int _i = sink->qcnt/NSSE;					\
	int _j = sink->qcnt%NSSE;					\
	sink->fmass += _m;						\
	Qvec[_i].mass[_j] = _m;						\
	Qvec[_i].x[_j] = r0;						\
	Qvec[_i].y[_j] = r1;						\
	Qvec[_i].z[_j] = r2;						\
	Qvec[_i].R[_j] = (p)->bmax;					\
	Qvec[_i].qx[_j] = (p)->qx;					\
	Qvec[_i].qy[_j] = (p)->qy;					\
	Qvec[_i].qz[_j] = (p)->qz;					\
	Qvec[_i].qxx[_j] = (p)->qxx;					\
	Qvec[_i].qxy[_j] = (p)->qxy;					\
	Qvec[_i].qyy[_j] = (p)->qyy;					\
	Qvec[_i].qxz[_j] = (p)->qxz;					\
	Qvec[_i].qyz[_j] = (p)->qyz;					\
	sink->qcnt++;							\
    } while(0)

#define appendHdvec(p)						\
    float _m;							\
    do {							\
	int _i = sink->hcnt/NSSE;				\
	int _j = sink->hcnt%NSSE;				\
	if (bs) {						\
	    _m = (p)->umass;					\
	    sink->fmass += _m;					\
	    Hvec[_i].mass[_j] = _m;				\
	    Hvec[_i].qxxxx[_j] = (p)->uxxxx;			\
	    Hvec[_i].qyyyy[_j] = (p)->uyyyy;			\
	    Hvec[_i].qxxyy[_j] = (p)->uxxyy;			\
	} else {						\
	    _m = (p)->mass;					\
	    sink->fmass += _m;					\
	    Hvec[_i].mass[_j] = _m;				\
	    Hvec[_i].qxxxx[_j] = (p)->qxxxx;			\
	    Hvec[_i].qyyyy[_j] = (p)->qyyyy;			\
	    Hvec[_i].qxxyy[_j] = (p)->qxxyy;			\
	}							\
	Hvec[_i].x[_j] = r0;					\
	Hvec[_i].y[_j] = r1;					\
	Hvec[_i].z[_j] = r2;					\
	Hvec[_i].R[_j] = (p)->bmax;				\
	Hvec[_i].qx[_j] = (p)->qx;				\
	Hvec[_i].qy[_j] = (p)->qy;				\
	Hvec[_i].qz[_j] = (p)->qz;				\
	Hvec[_i].qxx[_j] = (p)->qxx;				\
	Hvec[_i].qxy[_j] = (p)->qxy;				\
	Hvec[_i].qyy[_j] = (p)->qyy;				\
	Hvec[_i].qxz[_j] = (p)->qxz;				\
	Hvec[_i].qyz[_j] = (p)->qyz;				\
	Hvec[_i].qxxx[_j] = (p)->qxxx;				\
	Hvec[_i].qxxy[_j] = (p)->qxxy;				\
	Hvec[_i].qxyy[_j] = (p)->qxyy;				\
	Hvec[_i].qyyy[_j] = (p)->qyyy;				\
	Hvec[_i].qxxz[_j] = (p)->qxxz;				\
	Hvec[_i].qxyz[_j] = (p)->qxyz;				\
	Hvec[_i].qyyz[_j] = (p)->qyyz;				\
	Hvec[_i].qxxxy[_j] = (p)->qxxxy;			\
	Hvec[_i].qxyyy[_j] = (p)->qxyyy;			\
	Hvec[_i].qxxxz[_j] = (p)->qxxxz;			\
	Hvec[_i].qxxyz[_j] = (p)->qxxyz;			\
	Hvec[_i].qxyyz[_j] = (p)->qxyyz;			\
	Hvec[_i].qyyyz[_j] = (p)->qyyyz;			\
	sink->hcnt++;						\
    } while(0)

#define appendQvec(p)				\
    do {					\
	int _i = sink->qcnt/NSSE;					\
	int _j = sink->qcnt%NSSE;					\
	Qvec[_i].mass[_j] = (p)->mass;					\
	Qvec[_i].x[_j] = r0;						\
	Qvec[_i].y[_j] = r1;						\
	Qvec[_i].z[_j] = r2;						\
	Qvec[_i].R[_j] = (p)->bmax;					\
	Qvec[_i].qxx[_j] = (p)->qxx;					\
	Qvec[_i].qxy[_j] = (p)->qxy;					\
	Qvec[_i].qyy[_j] = (p)->qyy;					\
	Qvec[_i].qxz[_j] = (p)->qxz;					\
	Qvec[_i].qyz[_j] = (p)->qyz;					\
	sink->qcnt++;							\
    } while(0)

#define appendHvec(p) \
    do { \
      int _i = sink->hcnt/NSSE; \
      int _j = sink->hcnt%NSSE;	 \
      Hvec[_i].mass[_j] = (p)->mass; \
      Hvec[_i].x[_j] = r0; \
      Hvec[_i].y[_j] = r1; \
      Hvec[_i].z[_j] = r2; \
      Hvec[_i].R[_j] = (p)->bmax; \
      Hvec[_i].qxx[_j] = (p)->qxx; \
      Hvec[_i].qxy[_j] = (p)->qxy; \
      Hvec[_i].qyy[_j] = (p)->qyy; \
      Hvec[_i].qxz[_j] = (p)->qxz; \
      Hvec[_i].qyz[_j] = (p)->qyz; \
      Hvec[_i].qxxx[_j] = (p)->qxxx; \
      Hvec[_i].qxxy[_j] = (p)->qxxy; \
      Hvec[_i].qxyy[_j] = (p)->qxyy; \
      Hvec[_i].qyyy[_j] = (p)->qyyy; \
      Hvec[_i].qxxz[_j] = (p)->qxxz; \
      Hvec[_i].qxyz[_j] = (p)->qxyz; \
      Hvec[_i].qyyz[_j] = (p)->qyyz; \
      Hvec[_i].qxxxx[_j] = (p)->qxxxx; \
      Hvec[_i].qxxxy[_j] = (p)->qxxxy; \
      Hvec[_i].qxxyy[_j] = (p)->qxxyy; \
      Hvec[_i].qxyyy[_j] = (p)->qxyyy; \
      Hvec[_i].qyyyy[_j] = (p)->qyyyy; \
      Hvec[_i].qxxxz[_j] = (p)->qxxxz; \
      Hvec[_i].qxxyz[_j] = (p)->qxxyz; \
      Hvec[_i].qxyyz[_j] = (p)->qxyyz; \
      Hvec[_i].qyyyz[_j] = (p)->qyyyz; \
      sink->hcnt++;		\
    } while(0)

#define ddot(a, b) (pow2(a##0-b[0])+pow2(a##1-b[1])+pow2(a##2-b[2]))

/* RcritMAC with Don't Laugh-like traversal */
void
DLRcritMACsb(Sink *sink, const hcell **source_vec, int *restrict flags_vec, int *restrict result, int n)
{
    float dr2;
    Vxd(float r);

    StartTimer(&MACTm);
    for (int i = 0; i < n; i++) {
	int sf = Sub_Flags(source_vec[i]);
	if (!sf && !sink->isbody) {
	    result[i] = MAC_SPLIT_SINK;
	    continue;
	}
	const cell *cp = source_vec[i]->ptr;
	const quadcell *qcp = source_vec[i]->ptr;
	const hexacell *hcp = source_vec[i]->ptr;
	int bs = flags_vec[i] & BACKGROUND_FLAG;
	VxVV(r, = cp->pos, + offset_array[offset_index(flags_vec[i])]);
	if (sf) {
	    float bmax;
	    int isquad = mac->qcut && cp->daughters >= mac->qcut;
	    int ishexa = mac->hcut && cp->daughters >= mac->hcut;
	    float smallest_rcrit = ishexa ? hcp->rcrit_h : (isquad ? qcp->rcrit_q : cp->rcrit);
	    /* Mvec does not compute dipole, so don't test cp->rcrit here if doing geometric_center */
	    if (sink->clevel != CHUBITS && cp->level < sink->clevel) {
		dr2 = ddot(r, sink->cen);
		bmax = 1.2f*sink->cr;
	    } else {
		dr2 = ddot(r, sink->pos);
		bmax = sink->bmax;
	    }
	    if (!bs && cp->level == sink->clevel && ddot(r, sink->cen) < sink->cr2) {
		bs = BACKGROUND_FLAG;
		flags_vec[i] |= bs;
		sink->near++;
		DebugWatchKey("b %12g %12g Near %ld %s\n", 0.0, sqrt(ddot(r, sink->cen)), (long int)cp->daughters, PrintKey(source_vec[i]->key));
	    } 
	    if (isquad && dr2 > Square(qcp->rcrit_q + bmax)) {
		appendQdvec(qcp);
		DebugWatchKey("%s %12g %12g Qvec %ld %s\n", bs ? "b" : " ", _m, sqrt(dr2), (long int)qcp->daughters, PrintKey(source_vec[i]->key));
		sink->interactions += qcp->daughters;
		result[i] = MAC_ACCEPT;
	    } else if (ishexa && dr2 > Square(hcp->rcrit_h + bmax)) {
		appendHdvec(hcp);
		DebugWatchKey("%s %12g %12g Hvec %ld %s %g %g\n", bs ? "b" : " ", _m, sqrt(dr2), (long int)hcp->daughters, PrintKey(source_vec[i]->key), hcp->rcrit_h, bmax);
		sink->interactions += hcp->daughters;
		result[i] = MAC_ACCEPT;
	    } else {
		result[i] = mac->dlfac * sink->bmax > smallest_rcrit ? MAC_SPLIT_SINK : MAC_SPLIT_SRC;
	    }
	} else {
	    /* body-body */
	    dr2 = ddot(r, sink->pos);
	    if (dr2 > Eps2) {
		appendMvec(cp, cp->mass);
	    } else if (dr2 > 0.0f) appendSvec(cp);
	    sink->fmass += cp->mass;
	    DebugWatchKey("  %12g %12g Mvec 1 %s\n", cp->mass, sqrt(dr2), PrintKey(source_vec[i]->key));
	    sink->interactions++;
	    result[i] = MAC_ACCEPT;
	    if (!bs) {
		float mxyz[4], cellsz;
		CellCorner(source_vec[i]->key, mxyz+1, &cellsz);
		VV((mxyz+1), += 0.5f*cellsz + offset_array[offset_index(flags_vec[i])]);
		if (fabsf(mxyz[1]-sink->cen[0]) > sink->cr ||
		    fabsf(mxyz[2]-sink->cen[1]) > sink->cr ||
		    fabsf(mxyz[3]-sink->cen[2]) > sink->cr) {
		    mxyz[0] = -mac->rho0*cellsz*cellsz*cellsz;
		    appendMvec(mxyz, mxyz[0]);
		    sink->fmass += mxyz[0];
		    AddCounter(&MCAnti, 1);
		    DebugWatchKey("b %12g %12g Mvec anti\n", mxyz[0], sqrt(ddot(r, sink->cen)));
		} else {
		    DebugWatchKey("  %12g %12g Mvec anti Near %s \n", 0.0, sqrt(ddot(r, sink->cen)), PrintKey(source_vec[i]->key));
		    sink->near++;
		}
	    }
	}
	if (result[i] == MAC_SPLIT_SRC && !bs) {
	    if (cp->level >= sink->clevel) {
		float um = -ucell[cp->level].mass;
		appendMvec(cp, um);
		sink->fmass += um;
		AddCounter(&MCCorr, 1);
		flags_vec[i] |= BACKGROUND_FLAG;
		DebugWatchKey("  %12g %12g Cvec %ld %s\n", um, sqrt(dr2), (long int)cp->daughters, PrintKey(source_vec[i]->key));
	    } else if (sf != (1<<MAXNSUB) - 1) {
		Key_t k = KeyLshift(source_vec[i]->key, NDIM);
		for (int j = 0; j < MAXNSUB; sf >>= 1, k.k[0]++, j++) {
		    if ((sf & 1) == 0) {
			float mxyz[4], cellsz;
			CellCorner(k, mxyz+1, &cellsz);
			VV((mxyz+1), += 0.5f*cellsz + offset_array[offset_index(flags_vec[i])]);
			if (fabsf(mxyz[1]-sink->cen[0]) > sink->cr ||
			    fabsf(mxyz[2]-sink->cen[1]) > sink->cr ||
			    fabsf(mxyz[3]-sink->cen[2]) > sink->cr) {
			    mxyz[0] = -ucell[cp->level+1].mass;
			    appendMvec(mxyz, mxyz[0]);
			    sink->fmass += mxyz[0];
			    AddCounter(&CEmpty, 1);
			    DebugWatchKey("  %12g %12g %12g %12g Mvec empty %s %d\n", mxyz[0], mxyz[1]-sink->cen[0], mxyz[2]-sink->cen[1], mxyz[3]-sink->cen[2], PrintKey(k), offset_index(flags_vec[i]));
			} else {
			    DebugWatchKey("  %12g %12g %12g %12g Mvec empty Near %s %d\n", mxyz[0], mxyz[1]-sink->cen[0], mxyz[2]-sink->cen[1], mxyz[3]-sink->cen[2], PrintKey(k), offset_index(flags_vec[i]));
			    sink->near++;
			}
		    }
		}
	    }
	}
    }
    StopTimer(&MACTm);
    
    if (sink->scnt/NSSE >= SVECSZ) Error("svec overflow\n");
    if (sink->mcnt/NSSE >= MVECSZ) Error("mvec overflow\n");
    if (sink->qcnt/NSSE >= QVECSZ) Error("qvec overflow\n");
    if (sink->hcnt/NSSE >= HVECSZ) Error("hvec overflow\n");
}

/* RcritMAC with Don't Laugh-like traversal */
void
DLRcritMAC(Sink *sink, const hcell **source_vec, int *restrict flags_vec, int *restrict result, int n)
{
    float dr2;
    Vxd(float r);

    StartTimer(&MACTm);
    for (int i = 0; i < n; i++) {
	int sf = Sub_Flags(source_vec[i]);
	if (!sf && !sink->isbody) {
	    result[i] = MAC_SPLIT_SINK;
	    continue;
	}
	const cell *cp = source_vec[i]->ptr;
	const quadcell *qcp = source_vec[i]->ptr;
	const hexacell *hcp = source_vec[i]->ptr;
	VxVV(r, = cp->pos, + offset_array[offset_index(flags_vec[i])]);
	dr2 = ddot(r, sink->pos);
	if (sf) {
	    int isquad = mac->qcut && cp->daughters >= mac->qcut;
	    int ishexa = mac->hcut && cp->daughters >= mac->hcut;
	    float smallest_rcrit = ishexa ? hcp->rcrit_h : (isquad ? qcp->rcrit_q : cp->rcrit);
	    if (dr2 > Square(cp->rcrit + cp->bmax)) {
		appendMvec(cp, cp->mass);
		sink->interactions += cp->daughters;
		result[i] = MAC_ACCEPT;
	    } else if (isquad && dr2 > Square(qcp->rcrit_q + qcp->bmax)) {
		appendQvec(qcp);
		sink->interactions += qcp->daughters;
		result[i] = MAC_ACCEPT;
	    } else if (ishexa && dr2 > Square(hcp->rcrit_h + hcp->bmax)) {
		appendHvec(hcp);
		sink->interactions += hcp->daughters;
		result[i] = MAC_ACCEPT;
	    } else {
		result[i] = mac->dlfac * sink->bmax > smallest_rcrit ? MAC_SPLIT_SINK : MAC_SPLIT_SRC;
	    }
	} else {
	    /* body-body */
	    if (dr2 > Eps2) appendMvec(cp, cp->mass);
	    else if (dr2 > 0.0f) appendSvec(cp);
	    sink->interactions++;
	    result[i] = MAC_ACCEPT;
	}
    }
    StopTimer(&MACTm);
    
    if (sink->scnt/NSSE >= SVECSZ) Error("svec overflow\n");
    if (sink->mcnt/NSSE >= MVECSZ) Error("mvec overflow\n");
    if (sink->qcnt/NSSE >= QVECSZ) Error("qvec overflow\n");
    if (sink->hcnt/NSSE >= HVECSZ) Error("hvec overflow\n");
}
