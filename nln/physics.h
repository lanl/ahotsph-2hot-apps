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

#ifndef NDIM
#define NDIM 3
#endif
/* #define SAVE_ACC */

typedef struct {
    float mass;			/* mass of body */
    float pos[NDIM];		/* position of body */
    float vel[NDIM];		/* velocity of body */
    float acc[NDIM];
    float phi;
    Key_t key;
    unsigned int ident;
    float nterms;
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
    float bmax, rcrit;
    int daughters;
    float padding_junk;
} cell, *cellptr;


/* This is the intermediate data structure used to construct cofm */
typedef struct{
    float mass;
    float pos[NDIM];
    float massinv;
    float bmax;
    float B2;
    float sz;
    int ndaughters;
} cofmdata;

typedef struct{
    float bmax;
    float pos[NDIM];
    float M0;
    float M1[NDIM];
    int isbody;
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
Key_t getkey(const body *);

/* In cofm.c */
void SetupCofm(int MACtype, float tol, float rel_tol);
void CofmFromDaugh(hcellptr hptr, hcellptr daughters[]);
void CellFromCofm(cell *cp, cofmdata *cmp);

/* In print.c */
char *PrintCellContents(const cell *cp);
char *PrintBodyContents(const body *bp);
char *PrintBodyContentsLong(const body *vp);
char *PrintBranch(const cofmdata *cmp);



/* In mac.c */
extern Timer_t GravTm, MACTm;
extern Counter_t CCInt, CBInt, BCInt, BBInt;
extern Counter_t CCIntRej;
extern Counter_t TranslateCnt;

void SetTol(float tol, float frac_tol, float newton_const, float eps, int gnobj);
void InheritSink(const Sink *from, Sink *to, hcell *pp);
void DLRcritMAC(Sink *sink, const hcell **source, int *result, int n);
void RcritMAC(Sink *sink, const hcell **source, int *result, int n);
void SetGravOffset(float *off);
void UnSetGravOffset(void);
void InheritSinkNlogN(const Sink *from, Sink *to, hcell *pp);

/* In grav.v */
void do_grav(const float *p, const float *end, const float *pos0, float *mass0,
	     float *acc0, float *phi0, const float *eps2p, int *ncut);
#endif
