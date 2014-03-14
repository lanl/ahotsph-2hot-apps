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
#include "segment.h"
#include "gravcuda.h"

Counter_t CCInt, CBInt, BSInt, BSMax, BCInt, BC2Int, BC4Int, BBInt;
Counter_t CEmpty, MCAnti, MCCorr;
Counter_t FBCInt, FBC2Int, FBC4Int, FBCFInt, FBC2FInt, FBC4FInt, LBC2Int;
Counter_t MACcnt, BBMACcnt, EmptyMACcnt, MACcnt0, MACcnt1, MACcnt2, MACcnt3;

Timer_t GravTm, PGravTm, GravSTm, GravMTm, GravMMTm, GravQTm, GravHTm, GravMFTm, GravQFTm, GravHFTm, GravQLTm;
Timer_t MACTm, CUDAWtTm;

#if 0
int64_t WatchId = 500000000;
// Key_t WatchKey = {.k = {13229148200, 0}};
//Key_t WatchKey = {.k = {30745751, 0}};
// Key_t WatchKey = {.k = {2129464, 0}};
Key_t WatchKey = {.k = {2330112, 0}};
/* Key_t WatchKey = {.k = {3772552, 0}}; */
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

#define BACKGROUND_FLAG (1<<10)
#define OFFSET_MASK (BACKGROUND_FLAG-1) /* needs to hold MAX_IMAGE */
#define MAX_IMAGE 125
#define offset_index(f) ((f) & OFFSET_MASK)

static vsf Eps2v;
static int64_t GNobj, Nobj;
static float Eps2, Eps;
static int Smooth_type;
static float GNewt;
#define offset_index(f) ((f) & OFFSET_MASK)
static const mac_s *mac;
static const mxn_s *MxN;
static int Nimage = 1;
static float offset_array[MAX_IMAGE][NDIM];
static tree_t *SinkTree;
static body *Btab;
static hexacell *Htab;

static void mxn_hexa(Sink *s, const hcell *pp);
static void mxn_quad(Sink *s, const hcell *pp);
static void mxn_quad2(Sink *s, const hcell *pp);
static void mxn_qquad(Sink *s, const hcell *pp);
static void mxn_mono(Sink *s, const hcell *pp);

static grav_f Sinteract;
static grav_f Minteract;
static grav_ff MMinteract;
static grav_f Qinteract;
static grav_qf QQinteract;
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

/* costs of interaction relative to monopole */
#define QUAD_COST 3
#define HEXA_COST 9

#define SVECSZ (5623) /* This should be dynamically extensible */
static struct {
    vsf mass, x, y, z;
} Svec[SVECSZ];

#define MVECSZ (8*1873) /* This should be dynamically extensible */
static struct {
    vsf mass, x, y, z;
} Mvec[MVECSZ];

#define MMVECSZ (1021) /* This should be dynamically extensible */
segment MMvec[MMVECSZ];

#ifdef QUAD
#define QVECSZ (7919)
typedef struct Q_s {
    vsf mass, x, y, z;
    vsf R;
#ifdef DIPOLE
    vsf qx, qy, qz;
#endif
    vsf qxx, qxy, qyy, qxz, qyz;
} Q_s;
static Q_s Qvec[QVECSZ];
#define QQVECSZ (7919)
int QQvec[QQVECSZ];
#endif

#ifdef HEXA
#define HVECSZ (24576*2)
static struct Hvec {
    vsf mass, x, y, z;
    vsf R;
#ifdef DIPOLE
    vsf qx, qy, qz;
#endif
    vsf qxx, qxy, qyy, qxz, qyz;
    vsf qxxx, qxxy, qxyy, qyyy, qxxz, qxyz, qyyz;
    vsf qxxxx, qxxxy, qxxyy, qxyyy, qyyyy, qxxxz, qxxyz, qxyyz, qyyyz;
} Hvec[HVECSZ];
#endif

static struct ucell_s {
    float halfsz, mass, x4, x2y2;
} ucell[CHUBITS];

#define MNSQ_MAX 65536
#define MNSQ_MAX_M 65536

typedef struct mnsq_s {
    int inuse;
    int sink_base;
    int m;
    int ss_len;
    int seg_n;
    int source_base;
    int source_n[MNSQ_MAX/64];
    segment ss_seg[MNSQ_MAX/64];
    uint16_t ss_index[MNSQ_MAX_M];
    segment seg[MNSQ_MAX];
} mnsq_s;

mnsq_s mnsq;

#define QNSQ_MAX 20164
#define QNSQ_MAX_M (512+256)

typedef struct qnsq_s {
    int inuse;
    int sink_base;
    int m;
    int source_len;
    segment ss_seg[QNSQ_MAX_M];
    Q_s source[QNSQ_MAX/NSSE];
} qnsq_s;

qnsq_s qnsq;

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
	Minteract = Arch(do_grav);
	if (mac->geometric_center) {
	    Qinteract = Arch(do_gravdq);
	    QQinteract = Arch(do_gravdqq);
	    Hinteract = (amd6100) ? Arch(do_gravdh_amd6100) : Arch(do_gravdh);
	} else {
	    Qinteract = Arch(do_gravq);
	    Hinteract = (amd6100) ? Arch(do_gravh_amd6100) : Arch(do_gravh);
	}
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
	    Minteract = Arch(do_grav_sK1);
	    MMinteract = Arch(do_gravmm_sK1);
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
    }
    Eps2 = Eps*Eps*pow(particle_mass, (float)(2./3.));
    Eps2v = (vsf)vsf_scalar(Eps2);
    float a1 = mac->r0 * (1.0f + mac->expand_root);
    float a3 = mac->rho0 * pow3(2.0f*a1);
    float a5 = a3 * pow2(2.0f*a1);
    float a7 = a5 * pow2(2.0f*a1);
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
WalkInitHSrc(hexacell *htab, int ncells)
{
    Htab = htab;
    WalkInitSrcCUDA((float *)htab, sizeof(hexacell), ncells);
}

void
WalkInitSink(tree_t *tp, body *btab, int64_t nobj, mxn_s *mxn)
{
    SinkTree = tp;
    Btab = btab;
    Nobj = nobj;
    MxN = mxn;
    WalkInitSinkCUDA((float *)btab, sizeof(body)/sizeof(float), nobj);
    memset(&mnsq, 0, sizeof(mnsq));
}

void
WalkTerminateSink(tree_t *tp, body *btab, int64_t nobj)
{
    WalkTerminateSinkCUDA((float *)btab, sizeof(body)/sizeof(float), nobj);

    /* scale by GNewt */
    for (int64_t i = 0; i < nobj; i++) {
	btab[i].phi *= GNewt;
	VS(btab[i].acc, *= GNewt);
    }
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
FirstBody(const hcell *pp)
{
    Key_t k = pp->key;
    const hcell *p;
    const tree_t *tp = SinkTree;
    int nsub = 1<<tp->ndim;
    int sub_flags, i;

    while ((sub_flags = Sub_Flags(pp))) {
	if (pp->type & (SHARED|NONLOCAL)) return NULL;
	k = KeyLshift(pp->key, tp->ndim);
	for (i = 0; i < nsub; i++) {
	    if (sub_flags & (1 << i)) {
		p = Find(tp, KeyOrInt(k, i));
		if (p == NULL) Error("FirstBody failed %s %d\n", PrintKey(k), i);
		pp = p;
		break;
	    }
	}
    }
    return pp->ptr;
}

body *
LastBody(const hcell *pp)
{
    Key_t k = pp->key;
    const hcell *p;
    const tree_t *tp = SinkTree;
    int nsub = 1<<tp->ndim;
    int sub_flags, i;

    while ((sub_flags = Sub_Flags(pp))) {
	if (pp->type & (SHARED|NONLOCAL)) return NULL;
	k = KeyLshift(pp->key, tp->ndim);
	for (i = nsub-1; i >= 0; i--) {
	    if (sub_flags & (1 << i)) {
		p = Find(tp, KeyOrInt(k, i));
		if (p == NULL) Error("LastBody failed %s %d %d\n", PrintKey(k), i, pp->type);
		pp = p;
		break;
	    }
	}
    }
    return pp->ptr;
}

int
mxn_poll(double *last)
{
    double now = MPMY_Wtime();
    if (now - *last > 0.0005) {	/* should be an MxN parameter */
	WalkPoll();
	*last = now;
	return 1;
    } return 0;
}

/* Queue an MxN interaction as an array of segments */
void
grav_mns_queue(int this_base, int this_m, const segment *this_seg, int this_seg_n, int this_source_n)
{
    if (!mnsq.inuse) {
	mnsq.sink_base = this_base;
	mnsq.inuse = 1;
    }
    int base = this_base - mnsq.sink_base;

    Msgf(("queue %3d base %5d m %5d seg_n %5d -- ",
	  mnsq.ss_len, mnsq.sink_base, mnsq.m, mnsq.seg_n));
    Msgf(("base %5d m %5d seg_n %5d source_n %5d\n",
	  this_base, this_m, this_seg_n, this_source_n));

    if (this_base < mnsq.sink_base + mnsq.m) {
	/* could fix this by updating entries already in queue */
	Error("repeated or out-of-sequence sink in grav_mns_queue\n");
    } else if (this_base > mnsq.sink_base + mnsq.m) {
	int gap_m = this_base - (mnsq.sink_base + mnsq.m);
	Msgf(("Filling gap of %d\n", gap_m));
	mnsq.ss_seg[mnsq.ss_len].base = mnsq.seg_n;
	mnsq.ss_seg[mnsq.ss_len].length = 0; /* fill gap with zero work placeholders */
	mnsq.source_n[mnsq.ss_len] = 0;
	for (int i = 0; i < gap_m; i++) {
	    mnsq.ss_index[mnsq.m + i] = mnsq.ss_len;
	}
	mnsq.m += gap_m;
	mnsq.ss_len++;
    }
    mnsq.m += this_m;
    memcpy(mnsq.seg + mnsq.seg_n, this_seg, this_seg_n * sizeof(segment));
    
    mnsq.ss_seg[mnsq.ss_len].base = mnsq.seg_n;
    mnsq.ss_seg[mnsq.ss_len].length = this_seg_n;
    mnsq.source_n[mnsq.ss_len] = this_source_n;
    if (base + this_m >= MNSQ_MAX_M) Error("mnsq.ss_index overflow\n");
    for (int i = 0; i < this_m; i++) {
	mnsq.ss_index[base+i] = mnsq.ss_len;
    }
    mnsq.seg_n += this_seg_n;
    mnsq.ss_len++;
    if (mnsq.seg_n >= MNSQ_MAX) Error("mnsq.seg_n overflow\n");
    if (mnsq.ss_len >= MNSQ_MAX/64) Error("mnsq.ss_len overflow\n");
}

void
grav_qns_queue(int this_base, int this_m, const Q_s *this_source, int this_source_n)
{
    if (!qnsq.inuse) {
	qnsq.sink_base = this_base;
	qnsq.inuse = 1;
    }
    int base = this_base - qnsq.sink_base;

    Msgf(("queue base %5d m %5d source_len %5d -- ",
	  qnsq.sink_base, qnsq.m, qnsq.source_len));
    Msgf(("base %5d m %5d source_len %5d\n",
	  this_base, this_m, this_source_n));

    if (this_base < qnsq.sink_base + qnsq.m) {
	Error("repeated or out-of-sequence sink in grav_qns_queue\n");
    } else if (this_base > qnsq.sink_base + qnsq.m) {
	int gap_m = this_base - (qnsq.sink_base + qnsq.m);
	Msgf(("Filling gap of %d\n", gap_m));
	for (int i = 0; i < gap_m; i++) {
	    qnsq.ss_seg[qnsq.m + i].base = qnsq.source_len;
	    qnsq.ss_seg[qnsq.m + i].length = 0;
	}
	qnsq.m += gap_m;
    }
    qnsq.m += this_m;
    int this_source_bytes = this_source_n * sizeof(Q_s)/NSSE;
    if (qnsq.source_len + this_source_n >= QNSQ_MAX) 
	Error("qnsq.source overflow\n");
    memcpy(qnsq.source + qnsq.source_len/NSSE, this_source, this_source_bytes);
    
    if (base + this_m >= QNSQ_MAX_M) Error("qnsq overflow\n");
    for (int i = 0; i < this_m; i++) {
	qnsq.ss_seg[base+i].base = qnsq.source_len;
	qnsq.ss_seg[base+i].length = this_source_n;
    }
    qnsq.source_len += this_source_n;
}

void
grav_mns_flush(float mass, float e, int *ncut)
{
    double last_poll = MPMY_Wtime();
    int q;

    StartTimer(&CUDAWtTm);
    while ((q = qallocCUDA(NULL)) < 0) mxn_poll(&last_poll);
    StopTimer(&CUDAWtTm);
    grav_mnss_CUDA("pMM_mnss_sK1", mnsq.sink_base, mnsq.m, 
		   mnsq.ss_index, mnsq.ss_seg, mnsq.ss_len,
		   mnsq.seg, mnsq.seg_n, mnsq.source_base, mnsq.source_n,
		   mass, e, ncut, q);
    
    mnsq.sink_base += mnsq.m;
    mnsq.inuse = 0;
    mnsq.m = 0;
    mnsq.ss_len = 0;
    mnsq.seg_n = 0;
}

void
grav_qns_flush(void)
{
    double last_poll = MPMY_Wtime();
    int q;

    StartTimer(&CUDAWtTm);
    while ((q = qallocCUDA(NULL)) < 0) mxn_poll(&last_poll);
    StopTimer(&CUDAWtTm);
    StartTimer(&GravQLTm);
    grav_qnss_CUDA("pQ1_qnss", qnsq.sink_base, qnsq.m, 
		   qnsq.ss_seg, (float *)qnsq.source, qnsq.source_len, q);
    
    qnsq.sink_base += qnsq.m;
    qnsq.inuse = 0;
    qnsq.m = 0;
    qnsq.source_len = 0;
    StopTimer(&GravQLTm);
}

void 
InheritSinkNlogN(const Sink *from, Sink *to, hcell *pp)
{
    if (to == NULL) {
	body *bp = pp->ptr;
	/* must init mtot or else you get quiet exceptions in asm code */
	float mtot = 0.0f;
	int nn, smooth_cnt = from->smooth_cnt;
	float acc[NDIM], phi;
	float accd[NDIM], phid;

	float e;

	DebugWatchId("---%12g %12g %12g %ld %ld\n", bp->pos[0], bp->pos[1], bp->pos[2], 
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
	if (bp-Btab >= Nobj-1) grav_mns_flush(Btab[0].mass, e, &smooth_cnt);
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
			  from->pos, &mtot, accd, &phid, &e, NULL);
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
			  from->pos, &mtot, accd, &phid, &e, NULL);
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
	if (nn) Minteract((float *)&Mvec[from->mcnt_done/NSSE], (float *)&Mvec[nn/NSSE], 
			  from->pos, &mtot, accd, &phid, &e, &smooth_cnt);
	DebugWatchId("p1 %12g %12g %12g %d\n", accd[0], accd[1], accd[2], from->mcnt);
	AddCounter(&BCInt, from->mcnt-from->mcnt_done);
	StopTimer(&GravMTm);
	VV(acc, += accd); phi += phid;

	StartTimer(&GravMTm);
	nn = from->mmcnt-from->mmcnt_done;
	VS(accd, = 0.0); phid = 0.0;
	if (nn) MMinteract(Btab[0].pos, sizeof(Btab[0])/sizeof(float), Btab[0].mass, MMvec+from->mmcnt_done, nn, 
			   from->mmterms-from->mmterms_done, from->pos, &mtot, accd, &phid, &e, &smooth_cnt);
	DebugWatchId("p1 %12g %12g %12g %d\n", accd[0], accd[1], accd[2], from->mcnt);
	AddCounter(&BCInt, from->mmterms-from->mmterms_done);
	VV(acc, += accd); phi += phid;
	StopTimer(&GravMTm);

	StartTimer(&GravSTm);
	nn = from->scnt;
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
			  from->pos, &mtot, accd, &phid, &e, &smooth_cnt);
	if (from->scnt > BSMax.counter) BSMax.counter = from->scnt;
	AddCounter(&BSInt, from->scnt);
	if (smooth_cnt > BSMax.counter) BSMax.counter = smooth_cnt;
	AddCounter(&BSInt, smooth_cnt); /* Includes self */
	StopTimer(&GravSTm);
	StopTimer(&GravTm);
	VV(acc, += accd); phi += phid;
	
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
#if 0
	    double cmass = -mac->rho0*pow3(2.0f*from->cr);
	    DebugWatchId("fmass %12g cmass %12g\n", from->fmass, cmass);
	    if (fabs(from->fmass+cmass) > mac->m0*1e-7 || -cmass < bp->mass) 
		if (Nwarn++ < 2) SeriousWarning("Background subtraction off by %.2f for id %ld key %s {%lu,%lu} fmass %g cmass %g near %d\n", fabs(from->fmass+cmass)/bp->mass, bp->ident, PrintKey(pp->key), pp->key.k[0], pp->key.k[1], from->fmass/bp->mass, cmass/bp->mass, from->near);
	    if (from->near != 27) {
		Error("Bad Near for id %ld key %s {%lu,%lu} fmass %g cmass %g near %d\n", bp->ident, PrintKey(pp->key), pp->key.k[0], pp->key.k[1], from->fmass/bp->mass, cmass/bp->mass, from->near);
	    }
#endif
	}
	if (!isfinite(acc[0]) || !isfinite(acc[1]) || !isfinite(acc[2]) || !isfinite(phi)) {
	    Error("bad results from do_grav for id %ld (%g,%g,%g), ax=%g ay=%g az=%g phi=%g\n", bp->ident, from->pos[0], from->pos[1], from->pos[2],
		  acc[0], acc[1], acc[2], phi);
	}

	/* Make sure these are initialized to zero externally */
	bp->phi += from->M0;
	bp->phi += phi;
	VV(bp->acc, -= from->M1);
	VV(bp->acc, += acc);
	DebugWatchId("a %12g %12g %12g\n", bp->acc[0], bp->acc[1], bp->acc[2]);
	// Msgf(("a %5ld %12g %12g %12g\n", bp-Btab, bp->acc[0], bp->acc[1], bp->acc[2]));
	bp->nterms += from->nterms + from->scnt + from->mcnt + from->mmterms 
#ifdef QUAD
	    + QUAD_COST*from->qcnt + QUAD_COST*from->qqcnt 
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
	if (cp->daughters <= 32 || 4.0001f*mac->r0/mac->nx >= ucell[cp->level].halfsz) ccp = cp;
    } else {
	body *bp = pp->ptr;
	VV(to->pos, = bp->pos);
	to->bmax = 0.0f;
	to->daughters = 1;
	to->isbody = 1;
    }

    if (from) {
	if (!from->processed /* && MPMY_Procnum() == MPMY_Nproc()/3 */) {
	    int msrc = from->mcnt-from->mcnt_done;
	    int mmsrc = from->mmterms-from->mmterms_done;
	    int qsrc = from->qcnt-from->qcnt_done;
	    int qqsrc = from->qqcnt-from->qqcnt_done;
	    int hsrc = from->hcnt-from->hcnt_done;
	    Msgf(("%ld m %d/%d mm %d/%d q %d/%d qq %d/%d h %d/%d\n", from->daughters, 
		  msrc, from->mcnt, mmsrc, from->mmterms, qsrc, from->qcnt, qqsrc, from->qqcnt, hsrc, from->hcnt));
	}
#ifdef HEXA
	if (from->hcnt >= NSSE*HVECSZ) Error("hvec overflow\n");
	if (!from->processed && mac->p4cut && MxN->hblock && from->daughters >= MxN->min_hsink && 
	    from->hcnt-from->hcnt_done >= MxN->min_hsrc) {
	    mxn_hexa((Sink *)from, from->pp);
	}
	to->hcnt = from->hcnt;
	to->hcnt_done = from->hcnt_done;
#endif
#ifdef QUAD
	if (from->qcnt >= NSSE*QVECSZ) Error("qvec overflow\n");
	if (!from->processed && mac->p2cut && MxN->hblock && from->daughters >= MxN->min_qsink && 
	    from->qcnt-from->qcnt_done >= MxN->min_qsrc) {
	    mxn_quad2((Sink *)from, from->pp);
	}
	to->qcnt = from->qcnt;
	to->qcnt_done = from->qcnt_done;
	if (from->qqcnt >= QQVECSZ) Error("qqvec overflow\n");
	if (!from->processed && mac->p2cut && MxN->hblock && from->daughters >= MxN->min_qsink && 
	    from->qqcnt-from->qqcnt_done >= MxN->min_qsrc) {
	    mxn_qquad((Sink *)from, from->pp);
	}
	to->qqcnt = from->qqcnt;
	to->qqcnt_done = from->qqcnt_done;
#endif
	if (from->mmcnt >= MMVECSZ) Error("mmvec overflow\n");
	if (!from->processed && MxN->hblock && from->daughters >= MxN->min_msink && 
	    from->mmterms-from->mmterms_done >= MxN->min_msrc) {
	    mxn_mono((Sink *)from, from->pp);
	}
	to->mmterms = from->mmterms;
	to->mmterms_done = from->mmterms_done;
	to->mmcnt = from->mmcnt;
	to->mmcnt_done = from->mmcnt_done;
	to->mcnt = from->mcnt;
	to->mcnt_done = from->mcnt_done;
	if (!from->processed) ((Sink *)from)->processed = 1;
	to->scnt = from->scnt;
	if (to->scnt >= NSSE*SVECSZ) Error("svec overflow\n");
	    
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
	    Msgf(("Inherit to %s, cr %g\n", PrintKey(pp->key), to->cr));
	to->interactions = from->interactions;
	if (to->interactions == Nimage*GNobj) 
	    to->done = 1;
	else
	    to->done = from->done;
	to->processed = 0;
	to->pp = pp;
	to->key = pp->key;
	to->nterms = from->nterms;
	to->M0 = from->M0;
	VV(to->M1, = from->M1);
	to->smooth_cnt = from->smooth_cnt;
	to->smooth_len = from->smooth_len;
    } else {
	memset(to, 0, sizeof(Sink));
	to->daughters = GNobj;
	to->pp = pp;
	to->clevel = CHUBITS;
	/* Assumes constant particle masses */
	if (Smooth_type == 0 || Smooth_type == 6) {
	    /* plummer, e ~ (rho*eps)^2 */
	    to->smooth_len = Eps*Eps*pow(Btab[0].mass, 2./3.);
	} else {
	    /* polynomial, e ~ 1/(rho*eps) */
	    to->smooth_len = pow(Btab[0].mass, -1./3.)/Eps;
	}
    }
}


void
mxn_mono(Sink *s, const hcell *pp)
{
    body *first = FirstBody(pp);
    body *last = LastBody(pp);

    if (!first || !last) return;
    last++;
    if (first < Btab || last > first+Nobj) Error("first/last out of range\n");
    StartTimer(&GravTm);
    StartTimer(&GravMFTm);

    int seg_n = s->mmcnt - s->mmcnt_done;
    int source_n = s->mmterms-s->mmterms_done;
#ifdef CUDA
    if (MxN->do_pM && first-Btab >= mnsq.sink_base + mnsq.m) { /* Don't do more than one level */
	grav_mns_queue(first-Btab, last-first, MMvec, seg_n, source_n);
	/* need 100k transfer to get 4 GB/sec on Titan (PCI-Express 2.0) */
	if (mnsq.seg_n >= 32768) grav_mns_flush(Btab[0].mass, s->smooth_len, &s->smooth_cnt);
    } else {
	float mtot = 0.0;
	for (body *p = first; p < last; p++) {
	    MMinteract(Btab[0].pos, sizeof(Btab[0])/sizeof(float), Btab[0].mass, MMvec+s->mmcnt_done, seg_n, 
		       s->mmterms-s->mmterms_done, p->pos, &mtot, p->acc, &p->phi, &s->smooth_len, &s->smooth_cnt);
	}
    }
#else
    {
	float mtot = 0.0;
	for (body *p = first; p < last; p++) {
	    MMinteract(Btab[0].pos, sizeof(Btab[0])/sizeof(float), Btab[0].mass, MMvec+s->mmcnt_done, seg_n, 
		       s->mmterms-s->mmterms_done, p->pos, &mtot, p->acc, &p->phi, &s->smooth_len, &s->smooth_cnt);
	}
    }
#endif
    AddCounter(&FBCInt, (last-first)*source_n);
    s->mmcnt_done = s->mmcnt;
    s->mmterms_done = s->mmterms;

    StopTimer(&GravMFTm);
    StopTimer(&GravTm);
}

void
mxn_quad2(Sink *s, const hcell *pp)
{
    body *first = FirstBody(pp);
    body *last = LastBody(pp);

    if (!first || !last) return;
    last++;
    if (first < Btab || last > first+Nobj) Error("first/last out of range\n");
    StartTimer(&GravTm);

    int source_n = (s->qcnt/NSSE-s->qcnt_done/NSSE)*NSSE;
    if (MxN->do_pQ) {
	if (MxN->do_pQL && last-first < 256 && first-Btab >= qnsq.sink_base + qnsq.m) { /* Don't do more than one level */
	    grav_qns_queue(first-Btab, last-first, &Qvec[s->qcnt_done/NSSE], source_n);
	    if (qnsq.m >= 512) grav_qns_flush();
	} else {
	    double last_poll = MPMY_Wtime();
	    int q, qinuse;;
	    StartTimer(&CUDAWtTm);
	    while ((q = qallocCUDA(&qinuse)) < 0) {
		if (!mxn_poll(&last_poll)) break;
	    }
	    StopTimer(&CUDAWtTm);
	    if (q >= 0 && qinuse < 24) {
		StartTimer(&GravQFTm);
		grav_mn_CUDA("pQ", &first->mass, first->acc, last-first, sizeof(body)/sizeof(float),
			     (float *)&Qvec[s->qcnt_done/NSSE], source_n, sizeof(Qvec[0])/(NSSE*sizeof(float)),
			     s->smooth_len, &s->smooth_cnt, q);
		StopTimer(&GravQFTm);
	    } else {
		StartTimer(&GravQTm);
		double last_poll = MPMY_Wtime();
		for (body *p = first; p < last; p++) {
		    float mtot = 0.0;
		    Qinteract((float *)&Qvec[s->qcnt_done/NSSE], (float *)&Qvec[s->qcnt/NSSE], p->pos, &mtot, p->acc, &p->phi, 
			      &s->smooth_len, &s->smooth_cnt);
		    if ((p-first+1) % 16 == 0) mxn_poll(&last_poll);
		}
		AddCounter(&FBC2FInt, (last-first)*source_n);
		StopTimer(&GravQTm);
	    }
	}
    } else {
	StartTimer(&GravQTm);
	double last_poll = MPMY_Wtime();
	for (body *p = first; p < last; p++) {
	    float mtot = 0.0;
	    Qinteract((float *)&Qvec[s->qcnt_done/NSSE], (float *)&Qvec[s->qcnt/NSSE], p->pos, &mtot, p->acc, &p->phi, 
		      &s->smooth_len, &s->smooth_cnt);
	    if ((p-first+1) % 16 == 0) mxn_poll(&last_poll);
	}
	StopTimer(&GravQTm);
    }
    AddCounter(&FBC2Int, (last-first)*source_n);
    s->qcnt_done = (s->qcnt/NSSE)*NSSE;

    StopTimer(&GravTm);
}

void
mxn_quad(Sink *s, const hcell *pp)
{
    body *p;
    body *first = FirstBody(pp);
    body *last = LastBody(pp);
    int i, n0, n1, m_block, block;
    double last_poll = MPMY_Wtime();

    if (!first || !last) return;
    last++;
    if (first < Btab || last > first+Nobj) Error("first/last out of range\n");
    StartTimer(&GravTm);
    StartTimer(&GravQFTm);
    n0 = s->qcnt_done/NSSE;
    n1 = s->qcnt/NSSE;
    /* Size MxN->hblock for appropriate WalkPoll() latency */
    /* If Walk Defer timer is large, make MxN->hblock smaller */
    if (n1-n0 > MxN->hblock) m_block = 1;
    else if (n1 == n0) Error("mxn_quad called with n == 0\n");
    else m_block = MxN->hblock / (n1-n0);
    block = m_block;
    for (p = first; p < last; p += m_block) {
	if (p + block > last) block = last-p;
	if (MxN->do_pQ) {
#ifdef CUDA
	    int q, qinuse;;
	    StartTimer(&CUDAWtTm);
	    while ((q = qallocCUDA(&qinuse)) < 0) {
		if (!mxn_poll(&last_poll)) break;
	    }
	    StopTimer(&CUDAWtTm);
	    if (q >= 0 && qinuse < 24) {
		grav_mn_CUDA("pQ", &p->mass, p->acc, block, sizeof(body)/sizeof(float),
			     (float *)&Qvec[n0], (n1-n0)*NSSE, sizeof(Qvec[0])/(NSSE*sizeof(float)),
			     s->smooth_len, &s->smooth_cnt, q);
	    } else {
		float mtot = 0.0f;
		for (i = 0; i < block; i++) {
		    Qinteract((float *)&Qvec[n0], (float *)&Qvec[n1],
			      (p+i)->pos, &mtot, (p+i)->acc, &(p+i)->phi, &s->smooth_len, &s->smooth_cnt);
		    if ((i+1) % 16 == 0) mxn_poll(&last_poll); /* ~ ratio in speed between GPU core/ CPU core */
		}
		AddCounter(&FBC2FInt, block*(n1-n0)*NSSE);
	    }
#else
	    pQinteract(&p->mass, p->acc, block, sizeof(body)/sizeof(float),
		       (float *)&Qvec[n0], n1-n0);
#endif
	} else {
	    float mtot = 0.0f;
	    float e = 0.0f;
	    int ijunk = 0;
	    for (i = 0; i < block; i++) {
		Qinteract((float *)&Qvec[n0], (float *)&Qvec[n1],
			  (p+i)->pos, &mtot, (p+i)->acc, &(p+i)->phi, &e, &ijunk);
	    }
	}
	mxn_poll(&last_poll); /* ~ ratio in speed between GPU core/ CPU core */
    }
    // Msg_do("%ld %g %g %g\n", first->ident, first->acc[0], first->acc[1], first->acc[2]);
    AddCounter(&FBC2Int, (last-first)*(n1-n0)*NSSE);
    s->qcnt_done = n1*NSSE;

    StopTimer(&GravQFTm);
    StopTimer(&GravTm);
}

void
mxn_qquad(Sink *s, const hcell *pp)
{
    body *first = FirstBody(pp);
    body *last = LastBody(pp);

    if (!first || !last) return;
    last++;
    if (first < Btab || last > first+Nobj) Error("first/last out of range\n");
    StartTimer(&GravTm);
    StartTimer(&GravQLTm);

    int source_n = s->qqcnt - s->qqcnt_done;
#ifdef CUDA
    if (MxN->do_pQ) {
	double last_poll = MPMY_Wtime();
	int q;

	StartTimer(&CUDAWtTm);
	while ((q = qallocCUDA(NULL)) < 0) mxn_poll(&last_poll);
	StopTimer(&CUDAWtTm);
	grav_qns_CUDA("pQ1", first-Btab, last-first, QQvec + s->qqcnt_done, source_n, q);
    } else {
	float mtot = 0.0;
	for (body *p = first; p < last; p++) {
	    QQinteract((float *)Htab, sizeof(hexacell), QQvec+s->qqcnt_done, source_n, p->pos, &mtot, p->acc, &p->phi);
	}
    }
#else
    {
	float mtot = 0.0;
	for (body *p = first; p < last; p++) {
	    QQinteract((float *)Htab, sizeof(hexacell), QQvec+s->qqcnt_done, source_n, p->pos, &mtot, p->acc, &p->phi);
	}
    }
#endif
    AddCounter(&LBC2Int, source_n);
    s->qqcnt_done = s->qqcnt;

    StopTimer(&GravQLTm);
    StopTimer(&GravTm);
}


void
mxn_hexa(Sink *s, const hcell *pp)
{
    body *p;
    body *first = FirstBody(pp);
    body *last = LastBody(pp);
    int i, n0, n1, m_block, block;
    double last_poll = MPMY_Wtime();

    if (!first || !last) return;
    last++;
    if (first < Btab || last > first+Nobj) Error("first/last out of range\n");
    StartTimer(&GravTm);
    StartTimer(&GravHFTm);
    n0 = s->hcnt_done/NSSE;
    n1 = s->hcnt/NSSE;
    /* Size MxN->hblock for appropriate WalkPoll() latency */
    /* If Walk Defer timer is large, make MxN->hblock smaller */
    if (n1-n0 > MxN->hblock) m_block = 1;
    else if (n1 == n0) Error("mxn_hexa called with n == 0\n");
    else m_block = MxN->hblock / (n1-n0);
    block = m_block;
    for (p = first; p < last; p += m_block) {
	if (p + block > last) block = last-p;
	if (MxN->do_pH) {
#ifdef CUDA
	    int q;
	    StartTimer(&CUDAWtTm);
	    while ((q = qallocCUDA(NULL)) < 0) {
		if (!mxn_poll(&last_poll)) break;
	    }
	    StopTimer(&CUDAWtTm);
	    if (q >= 0) {
		grav_mn_CUDA("pH", &p->mass, p->acc, block, sizeof(body)/sizeof(float),
			     (float *)&Hvec[n0], (n1-n0)*NSSE, sizeof(Hvec[0])/(NSSE*sizeof(float)),
			     s->smooth_len, &s->smooth_cnt, q);
	    } else {
		float mtot = 0.0f;
		for (i = 0; i < block; i++) {
		    Hinteract((float *)&Hvec[n0], (float *)&Hvec[n1],
			      (p+i)->pos, &mtot, (p+i)->acc, &(p+i)->phi, &s->smooth_len, &s->smooth_cnt);
		    if ((i+1) % 16 == 0) mxn_poll(&last_poll); /* ~ ratio in speed between GPU core/ CPU core */
		}
		AddCounter(&FBC4FInt, block*(n1-n0)*NSSE);
	    }
#else
	    pHinteract(&p->mass, p->acc, block, sizeof(body)/sizeof(float),
		       (float *)&Hvec[n0], n1-n0);
#endif
	} else {
	    float mtot = 0.0f;
	    float e = 0.0f;
	    int nsmoothed = 0;
	    for (i = 0; i < block; i++) {
		Hinteract((float *)&Hvec[n0], (float *)&Hvec[n1],
			  (p+i)->pos, &mtot, (p+i)->acc, &(p+i)->phi, &e, &nsmoothed);
	    }
	}
	mxn_poll(&last_poll);
    }
    // Msg_do("%ld %g %g %g\n", first->ident, first->acc[0], first->acc[1], first->acc[2]);
    AddCounter(&FBC4Int, (last-first)*(n1-n0)*NSSE);
    s->hcnt_done = n1*NSSE;

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

    if (sink->done) {
	for (int i = 0; i < n; i++) result[i] = MAC_SPLIT_SINK;
	return;
    }

    StartTimer(&MACTm);
    AddCounter(&MACcnt, n);
    if (sink->daughters < 32) 
	AddCounter(&MACcnt0, n);
    else if (sink->daughters < 8*32) 
	AddCounter(&MACcnt1, n);
    else if (sink->daughters < 8*8*32) 
	AddCounter(&MACcnt2, n);
    else if (sink->daughters < 8*8*8*32) 
	AddCounter(&MACcnt3, n);

    for (int i = 0; i < n; i++) {
	int sf = Sub_Flags(source_vec[i]);
	const cell *cp = source_vec[i]->ptr;
	const quadcell *qcp = source_vec[i]->ptr;
	const hexacell *hcp = source_vec[i]->ptr;
	int bs = flags_vec[i] & BACKGROUND_FLAG;
	VxVV(r, = cp->pos, + offset_array[offset_index(flags_vec[i])]);
	if (sf) {
	    float bmax;
	    int isquad = mac->p2cut && cp->daughters >= mac->p2cut;
	    int ishexa = mac->p4cut && cp->daughters >= mac->p4cut;
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
		if (0 && !bs && MxN->do_pQL && offset_index(flags_vec[i]) == Nimage >> 1 && !(source_vec[i]->type & NONLOCAL)) {
		    QQvec[sink->qqcnt++] = hcp-Htab;
		} else {
		    appendQdvec(qcp);
		}
		DebugWatchKey("%s %12g %12g Qvec %ld %s\n", bs ? "b" : " ", _m, sqrt(dr2), (long int)qcp->daughters, PrintKey(source_vec[i]->key));
		sink->interactions += qcp->daughters;
		result[i] = MAC_ACCEPT;
	    } else if (ishexa && dr2 > Square(hcp->rcrit_h + bmax)) {
		appendHdvec(hcp);
		DebugWatchKey("%s %12g %12g Hvec %ld %s %g %g\n", bs ? "b" : " ", _m, sqrt(dr2), (long int)hcp->daughters, PrintKey(source_vec[i]->key), hcp->rcrit_h, bmax);
		sink->interactions += hcp->daughters;
		result[i] = MAC_ACCEPT;
	    } else {
		result[i] = (mac->dlfac * sink->bmax > smallest_rcrit && sink->daughters > mac->leaf_max_n) ? MAC_SPLIT_SINK : MAC_SPLIT_SRC;
	    }
	    DebugWatchKey("%10g %10g %10ld %s\n", sink->bmax, smallest_rcrit, (long int)cp->daughters, (result[i] == MAC_ACCEPT) ? "Accept" : 
			  (result[i] == MAC_SPLIT_SINK) ? "Split Sink" : "Split Source");
	} else {
	    /* body-body */
	    AddCounter(&BBMACcnt, 1);
	    appendMvec(cp, cp->mass);
	    sink->fmass += cp->mass;
	    DebugWatchKey("  %12g %12g Mvec 1 %s\n", cp->mass, sqrt(ddot(r, sink->pos)), PrintKey(source_vec[i]->key));
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
		DebugWatchKey("  %12g %12g Cvec %ld %s\n", um, sqrt(ddot(r, sink->cen)), (long int)cp->daughters, PrintKey(source_vec[i]->key));
	    } else if (sf != (1<<MAXNSUB) - 1 && cp->level > 4) { /* assumes empty cells near root are from expand_root */
		Key_t k = KeyLshift(source_vec[i]->key, NDIM);
		for (int j = 0; j < MAXNSUB; sf >>= 1, k.k[0]++, j++) {
		    if ((sf & 1) == 0) {
			float mxyz[4], cellsz;
			AddCounter(&EmptyMACcnt, 1);
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
	/* Optimization of terminal traversal */
	if (result[i] != MAC_ACCEPT && cp->daughters <= mac->leaf_max_n && !(source_vec[i]->type & (SHARED|NONLOCAL))) {
	    int idx = offset_index(flags_vec[i]);
	    if (MxN->do_pM && idx == Nimage >> 1) {
		MMvec[sink->mmcnt].base = cp->bptr-Btab;
		MMvec[sink->mmcnt++].length = cp->daughters;
		sink->mmterms += cp->daughters;
	    } else {
		for (body *b = cp->bptr; b < cp->bptr + cp->daughters; b++) {
		    VxVV(r, = b->pos, + offset_array[idx]);
		    appendMvec(b, b->mass);
		}
	    }
	    result[i] = MAC_ACCEPT;
	    sink->interactions += cp->daughters;
	}
    }
    StopTimer(&MACTm);
    
    if (sink->scnt/NSSE >= SVECSZ) Error("svec overflow\n");
    if (sink->mcnt/NSSE >= MVECSZ) Error("mvec overflow\n");
    if (sink->mmcnt >= MMVECSZ) Error("mmvec overflow\n");
    if (sink->qcnt/NSSE >= QVECSZ) Error("qvec overflow\n");
    if (sink->qqcnt >= QQVECSZ) Error("qqvec overflow\n");
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
	if (sf) {
	    dr2 = ddot(r, sink->pos);
	    int isquad = mac->p2cut && cp->daughters >= mac->p2cut;
	    int ishexa = mac->p4cut && cp->daughters >= mac->p4cut;
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
	    appendMvec(cp, cp->mass);
	    sink->interactions++;
	    result[i] = MAC_ACCEPT;
	}
        /* Optimization of terminal traversal */
        if (result[i] == MAC_SPLIT_SRC && cp->bptr && !(source_vec[i]->type & (SHARED|NONLOCAL))) {
            result[i] = MAC_ACCEPT;
            sink->interactions += cp->daughters;
            for (body *b = cp->bptr; b < cp->bptr + cp->daughters; b++) {
                VxVV(r, = b->pos, + offset_array[offset_index(flags_vec[i])]);
                if (b->mass != Btab[0].mass) {
                    Error("Bad particle mass %f type %d %s\n",
                          b->mass, source_vec[i]->type, PrintType(source_vec[i]->type));
                }
                appendMvec(b, b->mass);
	    }
	}
    }
    StopTimer(&MACTm);
    
    if (sink->scnt/NSSE >= SVECSZ) Error("svec overflow\n");
    if (sink->mcnt/NSSE >= MVECSZ) Error("mvec overflow\n");
    if (sink->qcnt/NSSE >= QVECSZ) Error("qvec overflow\n");
    if (sink->hcnt/NSSE >= HVECSZ) Error("hvec overflow\n");
}
