/*
 * Copyright 1991 Michael S. Warren and John K. Salmon.  All Rights Reserved.
 */

#ifndef physics_NdotH
#define physics_NdotH

#include "tree.h"
#include "key.h"
#include "timers.h"

#ifdef USE_PH
/* An ugly hack! */
#define CELLCORNER CellCornerPH
#define GETKEY GetKeyPH
#else
#define CELLCORNER CellCorner
#define GETKEY GetKey
#endif

#define NDIM 3

typedef struct {
    float mass;			/* mass of body */
    float pos[NDIM];		/* position of body */
    float vel[NDIM];		/* velocity of body */
    float acc[NDIM];
    float phi;
    Key_t key;
    unsigned int ident;
    float nterms;
    float acc_last;
} body, *bodyptr;

/* When we send a body from node to node, how much must we send??? */
#define TBODYSZ (1+NDIM)*sizeof(float)

typedef struct {
    float mass;			/* mass of body */
    float pos[NDIM];		/* position of body */
    float vel[NDIM];		/* velocity of body */
#ifdef SAVE_ACC
    float acc[NDIM];
    float phi;
#endif
    unsigned int ident;		/* unique? identifier */
} outbody, *outbodyptr;

/* This is the descriptor that goes into the SDF header. */

#ifdef SAVE_ACC
#if NDIM==3
#define OUTBODYDESC \
"struct {\n\
    float mass;			/* mass of body */\n\
    float x, y, z;			/* position of body */\n\
    float vx, vy, vz;		/* velocity of body */\n\
    float ax, ay, az;		/* acceleration */\n\
    float phi;			/* potential */\n\
    unsigned int ident;		/* unique? identifier */\n\
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
    unsigned int ident;		/* unique? identifier */\n\
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
    float x, y, z;			/* position of body */\n\
    float vx, vy, vz;		/* velocity of body */\n\
    unsigned int ident;		/* unique? identifier */\n\
}"
#else
#if NDIM==2
#define OUTBODYDESC \
"struct {\n\
    float mass;			/* mass of body */\n\
    float x, y;			/* position of body */\n\
    float vx, vy;		/* velocity of body */\n\
    unsigned int ident;		/* unique? identifier */\n\
}"
#else
 # error No case for NDIM
#endif /* NDIM==2 */
#endif /* NDIM==3 */
#endif /* SAVE_ACC */

typedef struct {
    float mass;
    float pos[NDIM];
    float rcrit;
    int daughters;
    float B2, B3, bmax;
    float acc_last_max;
} cell, *cellptr;

/* Now for Order N */

#if (NDIM==3)
typedef struct {
    float xx;
    float yy;
    float zz;
    float xy;
    float xz;
    float yz;
} moment;
#else
typedef struct {
    float xx;
    float yy;
    float xy;
} moment;
#endif

/* This is the intermediate data structure used to construct cofm */
typedef struct{
    float mass;
    float pos[NDIM];
    float B2;
    float bmax;
    float acc_last_max;
    float massinv;
    int ndaughters;
} cofmdata;

typedef struct{
    float bmax;
    float pos[NDIM];
    float M0;
    float M1[NDIM];
    float m;		/* This is only needed for testing */
    float acc_last_max;
    int isbody;
    moment M2;
    float daughters;
    float nterms;
    int interactions;
    int icnt;
} Sink;

/* Tell physics.c that we have nterms in the body struct */
#define HAS_NTERMS
#define HAS_IDENT
#define HAS_KEY

#define Mass(x) ((x)->mass)
#define Pos(x)  ((x)->pos)

/* Prototypes for all the functions which are "friends" of physics.h */
#include "physics_generic.h"

/* In main_n.c */
extern Timer_t StepTot;
extern Timer_t BuildTot;
extern Timer_t FindForcesTm;
extern Counter_t NbodyCnt;
Key_t getkey(const body *);

/* In cofm_n.c */
void cofm_setup(float tol);
void CofmFromDaugh(hcellptr hptr, hcellptr daughters[]);
void CellFromCofm(cell *cp, cofmdata *cmp);

/* In print.c */
char *PrintCellContents(const cell *cp);
char *PrintBodyContents(const body *bp);
char *PrintBodyContentsLong(const body *vp);
char *PrintBranch(const cofmdata *cmp);

/* In grav_nv.c */
extern Counter_t CCInt, CBInt, BCInt, BBInt;
extern Counter_t CCIntRej;
extern Counter_t TranslateCnt;

extern Timer_t GravTm;

void SetTol(float tol, float frac_tol, float newton_const, float eps, int gnobj);
void Unifiedmacv(Sink *sink, const hcell **source, int *result, int n);
void Fracmacv(Sink *sink, const hcell **source, int *result, int n);
void Lowestmacv(Sink *sink, const hcell **source, int *result, int n);
void InheritSink(const Sink *from, Sink *to, hcell *pp);
void SetGravOffset(float *off);
void UnSetGravOffset(void);
void InheritSinkNlogN(const Sink *from, Sink *to, hcell *pp);
void Nlogngate(Sink *sink, const hcell **source_vec, int *result, int n);
void do_grav(const float *p, const float *end, const float *pos0, float *mass0,
	     float *acc0, float *phi0, const float *eps2p, int *ncut);


/* In grav_ring.c */
void set_body(void *o, void *p);
void do_grav2(void *p0, void *list, int bsize, int n);
void set_eps(float eps);

#endif
