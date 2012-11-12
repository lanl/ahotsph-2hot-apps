/*
 * Copyright 1991 Michael S. Warren and John K. Salmon.  All Rights Reserved.
 */

#ifndef physics_NdotH
#define physics_NdotH

#include <sys/types.h>
#include "tree.h"
#include "key.h"
#include "timers.h"

/* #define SAVE_ACC */
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
    float bmax, R, rcrit;
    int padding;
    int64_t daughters;
} cell, *cellptr;

typedef struct {
    float mass;
    float pos[NDIM];
    float bmax, R, rcrit;
    int padding;
    int64_t daughters;
    float rcrit_q;
    float qxx, qxy, qyy, qxz, qyz;
} quadcell, *quadcellptr;

typedef struct {
    float mass;
    float pos[NDIM];
    float bmax, R, rcrit;
    int padding;
    int64_t daughters;
    float rcrit_q;
    float qxx, qxy, qyy, qxz, qyz;
    float rcrit_h;
    float qxxx, qxxy, qxyy, qyyy, qxxz, qxyz, qyyz;
    float qxxxx, qxxxy, qxxyy, qxyyy, qyyyy, qxxxz, qxxyz, qxyyz, qyyyz;
    float padding_junk;
} hexacell, *hexacellptr;

/* This is the intermediate data structure used to construct cofm */
typedef struct{
    float mass;
    float pos[NDIM];
    float massinv;
    float bmax;
    float B2;
    float sz;
    int64_t ndaughters;
#if defined(QUAD) || defined(HEXA)
    double x2, xy, y2, xz, yz, z2;
    double x3, x2y, xy2, y3, x2z, xyz, y2z, xz2, yz2, z3;
    double x4, x3y, x2y2, xy3, y4, x3z, x2yz, xy2z, y3z, x2z2, xyz2, y2z2, xz3, yz3, z4;
    double x5, x4y, x3y2, x2y3, xy4, y5, x4z, x3yz, x2y2z, xy3z, y4z, x3z2, x2yz2, xy2z2, y3z2, x2z3, xyz3, y2z3, xz4, yz4, z5;
    double x6, x4y2, x2y4, y6, x4z2, x2y2z2, y4z2, x2z4, y2z4, z6;
#endif
} cofmdata;

typedef struct{
    float bmax;
    float pos[NDIM];
    float M0;
    float M1[NDIM];
    int isbody;
    int64_t daughters;
    float nterms;
    int64_t interactions;
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

#define Mass(x) ((x)->mass)
#define Pos(x)  ((x)->pos)

#define BMAX_MAC 1
#define BH_MAC  2
#define AREL_MAC 3

/* Prototypes for all the functions which are "friends" of physics.h */
#include "physics_generic.h"

/* In main_n.c */
extern Timer_t StepTot;
extern Timer_t BuildTot;
extern Timer_t FindForcesTm;
extern Counter_t NbodyCnt;

/* In cofm.c */
void SetupCofm(int MACtype, float tol, float rel_tol, float rel_tol0, float R0, float ptol_boost,
	       float stol_max, int qcut, int hcut,  tree_t *t);
void CofmFromDaugh(hcellptr hptr, hcellptr daughters[]);
int CellSz(void *p);
void *CellFromCofm(cofmdata *cmp);

/* In print.c */
char *PrintCellContents(const hexacell *cp);
char *PrintCellContents4(const hexacell *cp);
char *PrintBodyContents(const body *bp);
char *PrintBodyContentsLong(const body *vp);
char *PrintBranch(const cofmdata *cmp);

/* In mac.c */
extern Timer_t GravTm, GravSTm, GravMTm, GravQTm, GravHTm, EwaldTm, MACTm;
extern Counter_t CCInt, BSInt, BSMax, CBInt, BCInt, BC2Int, BC4Int, BBInt;
extern Counter_t FBC2Int, FBC4Int;
extern Counter_t CCIntRej;
extern Counter_t TranslateCnt;

void SetupGrav(float newton_const, float eps, int64_t gnobj, float dl_fac, float dl_max,
	       int qcut, int hcut, float particle_mass, int smoothing_type);
void InheritSink(const Sink *from, Sink *to, hcell *pp);
void DLRcritMAC(Sink *sink, const hcell **source, const float **offset, int *result, int n);
void RcritMAC(Sink *sink, const hcell **source, const float **offset, int *result, int n);
void SetGravOffset(float *off, int nimage);
void WalkInitSink(tree_t *tp, body *btab, int64_t nobj, int mxn_hblock);
void WalkInitSrc(Stk *kstk, Stk *ostk);
void WalkInitSrcPeriodic(Stk *kstk, Stk *ostk);
void InheritSinkNlogN(const Sink *from, Sink *to, hcell *pp);

/* In grav.c */
typedef void (*grav_t)(const float *p, const float *end, const float *pos0, float *mass0,
		       float *acc0, float *phi0, const float *eps2p, int *ncut);

void do_gravh_sse4(const float *f, const float *fend, const float *pos0, float *mass0, float *acc, float *phi0, const float *e, int *ncut);
void do_gravh_amd6100_sse4(const float *f, const float *fend, const float *pos0, float *mass0, float *acc, float *phi0, const float *e, int *ncut);
void do_gravq_sse4(const float *f, const float *fend, const float *pos0, float *mass0, float *acc, float *phi0, const float *e, int *ncut);
void do_grav_sse4(const float *f, const float *fend, const float *pos0, float *mass0, float *acc, float *phi0, const float *e, int *ncut);
void do_grav_sse16_ivec_asm(const float *f, const float *fend, const float *pos0, float *mass0, float *acc, float *phi0, const float *e, int *ncut);
void do_gravsU_sse4(const float *f, const float *fend, const float *pos0, float *mass0, float *acc, float *phi0, const float *e, int *ncut);
void do_gravsF1_sse4(const float *f, const float *fend, const float *pos0, float *mass0, float *acc, float *phi0, const float *e, int *ncut);
void do_gravsF2_sse4(const float *f, const float *fend, const float *pos0, float *mass0, float *acc, float *phi0, const float *e, int *ncut);
void do_gravsK1_sse4(const float *f, const float *fend, const float *pos0, float *mass0, float *acc, float *phi0, const float *e, int *ncut);
void do_gravsS_sse4(const float *f, const float *fend, const float *pos0, float *mass0, float *acc, float *phi0, const float *e, int *ncut);
void do_gravsCP_sse4(const float *f, const float *fend, const float *pos0, float *mass0, float *acc, float *phi0, const float *e, int *ncut);

void do_gravph_sse4(const float *f, const float *fend, const float *pos0, float *mass0, float *acc, float *phi0, const float *e, int *ncut);
void do_gravph_amd6100_sse4(const float *f, const float *fend, const float *pos0, float *mass0, float *acc, float *phi0, const float *e, int *ncut);
void do_gravpq_sse4(const float *f, const float *fend, const float *pos0, float *mass0, float *acc, float *phi0, const float *e, int *ncut);
void do_gravp_sse4(const float *f, const float *fend, const float *pos0, float *mass0, float *acc, float *phi0, const float *e, int *ncut);
void pHinteract(const float *p, float *accp, const int n, const int stride, 
		const float *f, const int source_n);
void pQinteract(const float *p, float *accp, const int n, const int stride, 
		const float *f, const int source_n);

/* In ewald.c */
void ewald(double *x, double L, double *f, double *phi);

/* In grav_n2.c */
void EwaldSetup(double box_length, double Gnewt);
void EwaldForces(body *btab, int nobj, float sample_frac, void *ranstate);
#endif
