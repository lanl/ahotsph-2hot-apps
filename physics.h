/*
 * Copyright 1991 Michael S. Warren and John K. Salmon.  All Rights Reserved.
 */

#ifndef physics_NdotH
#define physics_NdotH

#include <sys/types.h>
#include "tree.h"
#include "key.h"
#include "timers.h"
#include "vec.h"

#define SAVE_ACC
/* #define BODY_HAS_KEY */

#ifdef USE_PH
/* An ugly hack! */
#define CELLCORNER CellCornerPH
#define GETKEY GetKeyPH
#else
#define CELLCORNER CellCorner
#define GETKEY GetKey
#endif

#ifndef NDIM
#define NDIM 3
#endif

typedef struct {
    float mass;			/* mass of body */
    float pos[NDIM];		/* position of body */
    float vel[NDIM];		/* velocity of body */
    float acc[NDIM];
    float phi;
    float nterms;
    int64_t ident;
#ifdef BODY_HAS_KEY
    Key_t key;
#endif
} __attribute__ ((packed)) body, *bodyptr;

/* When we send a body from node to node, how much must we send??? */
#define TBODYSZ (1+NDIM)*sizeof(float)

typedef struct  {
    float mass;			/* mass of body */
    float pos[NDIM];		/* position of body */
    float vel[NDIM];		/* velocity of body */
#ifdef SAVE_ACC
    float acc[NDIM];
    float phi;
#endif
    int64_t ident;		/* unique? identifier */
} __attribute__ ((packed)) outbody, *outbodyptr;

/* This is the descriptor that goes into the SDF header. */

#ifdef SAVE_ACC
#if NDIM==3
#define OUTBODYDESC \
"struct {\n\
    float mass;			/* mass of body */\n\
    float x, y, z;		/* position of body */\n\
    float vx, vy, vz;		/* velocity of body */\n\
    float ax, ay, az;		/* acceleration */\n\
    float phi;			/* potential */\n\
    int64_t ident;		/* unique? identifier */\n\
}"
#else
#if NDIM==2
#define OUTBODYDESC \
"struct {\n\
    float mass;			/* mass of body */\n\
    float x, y;			/* position of body */\n\
    float vx, vy;		/* velocity of body */\n\
    float ax, ay;		/* acceleration */\n\
    float phi;			/* potential */\n\
    int64_t ident;		/* unique? identifier */\n\
}"
#else
 # error No case for NDIM
#endif /* NDIM==2 */
#endif /* NDIM==3 */
#else
#if NDIM==3
#define OUTBODYDESC \
"struct {\n\
    float mass;			/* mass of body */\n\
    float x, y, z;		/* position of body */\n\
    float vx, vy, vz;		/* velocity of body */\n\
    int64_t ident;		/* unique? identifier */\n\
}"
#else
#if NDIM==2
#define OUTBODYDESC \
"struct {\n\
    float mass;			/* mass of body */\n\
    float x, y;			/* position of body */\n\
    float vx, vy;		/* velocity of body */\n\
    int64_t ident;		/* unique? identifier */\n\
}"
#else
 # error No case for NDIM
#endif /* NDIM==2 */
#endif /* NDIM==3 */
#endif /* SAVE_ACC */


/* The cell, quadcell and hexacells are type-punned in cofm, so don't */
/* independently mess with the common sections at the start */

typedef struct {
    float mass;
    float pos[NDIM];
    float rcrit, bmax;
    int64_t level : 16;
    int64_t daughters : 48;
} cell, *cellptr;

typedef struct {
    float mass;
    float pos[NDIM];
    float rcrit_m, bmax;
    int64_t level : 16;
    int64_t daughters : 48;
    float rcrit_q, pad;
#ifdef DIPOLE
    float qx, qy, qz;
#else
    int padding;
#endif
    float qxx, qxy, qyy, qxz, qyz;
} quadcell, *quadcellptr;

typedef struct {
    float mass;
    float pos[NDIM];
    float rcrit_m, bmax;
    int64_t level : 16;
    int64_t daughters : 48;
    float rcrit_q, rcrit_h;
#ifdef DIPOLE
    float qx, qy, qz;
#else
    int padding;
#endif
    float qxx, qxy, qyy, qxz, qyz;
    float qxxx, qxxy, qxyy, qyyy, qxxz, qxyz, qyyz;
    float qxxxx, qxxxy, qxxyy, qxyyy, qyyyy, qxxxz, qxxyz, qxyyz, qyyyz;
    float umass, uxxxx, uyyyy, uxxyy;
} hexacell, *hexacellptr;

/* This is the intermediate data structure used to construct cofm */
typedef struct{
    double m;
    double center[NDIM];
    double bmax;
    double B2;
    double sz;
    int64_t level : 16;
    int64_t ndaughters : 48;
#ifdef QUAD
    double x, y, z;
    double x2, xy, y2, xz, yz, z2;
    double x3, x2y, xy2, y3, x2z, xyz, y2z, xz2, yz2, z3;
    double x4, x3y, x2y2, xy3, y4, x3z, x2yz, xy2z, y3z, x2z2, xyz2, y2z2, xz3, yz3, z4;
#endif
#ifdef HEXA
    double x5, x4y, x3y2, x2y3, xy4, y5, x4z, x3yz, x2y2z, xy3z, y4z, x3z2, x2yz2, xy2z2, y3z2, x2z3, xyz3, y2z3, xz4, yz4, z5;
    double x6, x4y2, x2y4, y6, x4z2, x2y2z2, y4z2, x2z4, y2z4, z6; /* just the symmetric components */
#endif
} cofmdata;

typedef struct{
    float bmax;
    float pos[NDIM];
    int isbody;
    float nterms;
    float M0;
    float M1[NDIM];
    int64_t daughters;
    int64_t interactions;
    Key_t key;
    double fmass;
    int clevel;
    int near;
    float cen[NDIM];
    float cr, cr2;
    int scnt;
    int mcnt;
#ifdef QUAD
    int qcnt;
    int qcnt_done;
#endif
#ifdef HEXA
    int hcnt;
    int hcnt_done;
#endif
} Sink;

/* Tell physics.c that we have nterms in the body struct */
#define HAS_NTERMS
#define HAS_IDENT

enum mac_type {BMAX_MAC, BH_MAC, AREL_MAC};

typedef struct mac_s {
    /* specified */
    macv_t rcrit_func;
    tree_t *tree;
    enum mac_type type;
    float tol;		/* absolute acc */
    float this_tol;	/* absolute acc for this step */
    float rel_tol;
    float rel_tol0;
    float r0;
    float dlfac;
    float dlmax;
    float nx;
    float ptol_boost;
    float stol_max;
    union {
	int p2cut;
	int qcut;
    };
    union {
	int p4cut;
	int hcut;
    };
    int p8cut;
    int geometric_center;
    int subtract_background;
    double m0;
    /* derived */
    double rho0;
    float inv_tol;
    float inv_rel_tol;
    float inv_rel_tol0;
    float bmax0;
    float cr;
} mac_s;


/* Prototypes for all the functions which are "friends" of physics.h */
#include "physics_generic.h"

/* In main_n.c */
extern Timer_t StepTot;
extern Timer_t BuildTot;
extern Timer_t FindForcesTm;
extern Counter_t NbodyCnt;

/* In cofm.c */
void SetupCofm(mac_s *mac);
void CofmFromDaugh(hcellptr hptr, hcellptr daughters[]);
int CellSz(void *p);
void *CellFromCofm(cofmdata *cmp);

/* In print.c */
char *PrintCellContents(const cell *cp);
char *PrintCellContents4(const hexacell *cp);
char *PrintBodyContents(const body *bp);
char *PrintBodyContentsLong(const body *vp);
char *PrintBranch(const cofmdata *cmp);

/* In mac.c */
extern Timer_t GravTm, GravSTm, GravMTm, GravQTm, GravHTm, EwaldTm, MACTm, MACswzlTm;
extern Counter_t CCInt, BSInt, BSMax, CBInt, BCInt, BC2Int, BC4Int, BBInt;
extern Counter_t CEmpty, MCCorr, MCAnti;
extern Counter_t FBC2Int, FBC4Int;
extern Counter_t MACcnt, BBMACcnt;

void SetupGrav(float newton_const, float eps, int64_t gnobj, mac_s *mac,
	       float particle_mass, int smoothing_type);
void InheritSink(const Sink *from, Sink *to, hcell *pp);
void DLRcritMAC(Sink *sink, const hcell **source, int *flags, int *result, int n);
void DLRcritMACsb(Sink *sink, const hcell **source, int *flags, int *result, int n);
void RcritMAC(Sink *sink, const hcell **source, int *flags, int *result, int n);
void SetGravOffset(float *off, int nimage);
void WalkInitSink(tree_t *tp, body *btab, int64_t nobj, int mxn_hblock);
void WalkInitSrc(Stk *kstk, Stk *ostk);
void WalkInitSrcPeriodic(Stk *kstk, Stk *ostk);
void InheritSinkNlogN(const Sink *from, Sink *to, hcell *pp);

/* In grav.c */
typedef void (*grav_f)(const float *p, const float *end, const float *pos0, float *mass0,
		       float *acc0, float *phi0, const float *eps2p, int *ncut);

#define grav_decl(f) void Arch(f)(const float *f, const float *fend, const float *pos0, float *mass0, float *acc, float *phi0, const float *e, int *ncut)

grav_decl(do_gravdq);
grav_decl(do_gravdh);
grav_decl(do_gravdh_amd6100);

void Arch(do_gravh)(const float *f, const float *fend, const float *pos0, float *mass0, float *acc, float *phi0, const float *e, int *ncut);
void Arch(do_gravh_amd6100)(const float *f, const float *fend, const float *pos0, float *mass0, float *acc, float *phi0, const float *e, int *ncut);
void Arch(do_gravq)(const float *f, const float *fend, const float *pos0, float *mass0, float *acc, float *phi0, const float *e, int *ncut);
void Arch(do_grav)(const float *f, const float *fend, const float *pos0, float *mass0, float *acc, float *phi0, const float *e, int *ncut);
void do_grav_sse16_ivec_asm(const float *f, const float *fend, const float *pos0, float *mass0, float *acc, float *phi0, const float *e, int *ncut);
void Arch(do_gravsU)(const float *f, const float *fend, const float *pos0, float *mass0, float *acc, float *phi0, const float *e, int *ncut);
void Arch(do_gravsF1)(const float *f, const float *fend, const float *pos0, float *mass0, float *acc, float *phi0, const float *e, int *ncut);
void Arch(do_gravsF2)(const float *f, const float *fend, const float *pos0, float *mass0, float *acc, float *phi0, const float *e, int *ncut);
void Arch(do_gravsK1)(const float *f, const float *fend, const float *pos0, float *mass0, float *acc, float *phi0, const float *e, int *ncut);
void Arch(do_gravsS)(const float *f, const float *fend, const float *pos0, float *mass0, float *acc, float *phi0, const float *e, int *ncut);
void Arch(do_gravsCP)(const float *f, const float *fend, const float *pos0, float *mass0, float *acc, float *phi0, const float *e, int *ncut);

void Arch(do_gravph)(const float *f, const float *fend, const float *pos0, float *mass0, float *acc, float *phi0, const float *e, int *ncut);
void Arch(do_gravph_amd6100)(const float *f, const float *fend, const float *pos0, float *mass0, float *acc, float *phi0, const float *e, int *ncut);
void Arch(do_gravpq)(const float *f, const float *fend, const float *pos0, float *mass0, float *acc, float *phi0, const float *e, int *ncut);
void Arch(do_gravp)(const float *f, const float *fend, const float *pos0, float *mass0, float *acc, float *phi0, const float *e, int *ncut);

void pHinteract(const float *p, float *accp, const int n, const int stride, 
		const float *f, const int source_n);
void pQinteract(const float *p, float *accp, const int n, const int stride, 
		const float *f, const int source_n);

#ifdef CUDA
void CUDA_Init(void);
void pinteractCUDA(const float *p, float *accp, const int n, const int stride, 
		   const float *f, const int source_n, const int sz);
#endif


/* In ewald.c */
void ewald(double *x, double L, double *f, double *phi);

/* In grav_n2.c */
void EwaldSetup(double box_length, double Gnewt);
void EwaldForces(body *btab, int nobj, float sample_frac, void *ranstate);
#endif
