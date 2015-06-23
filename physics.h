/*
 * Copyright 1991 Michael S. Warren and John K. Salmon.  All Rights Reserved.
 */

#ifndef physics_NdotH
#define physics_NdotH

#include <tree.h>
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

#include "ndim.h"
#define SAVE_ACC

typedef struct {
    double mass;			/* mass of body */
    double pos[NDIM];		/* position of body */
    double vel[NDIM];		/* velocity of body */
    double h;
    double acc[NDIM];       /* EVERYTHING ABOVE ACC WILL BE COMMUNICATED! */
    double phi;
    Key_t key;
    unsigned int ident;
    unsigned int type;
    double nterms;
    double pos_last[NDIM];	/* position of body */
} body, *bodyptr;

/* When we send a body from node to node, how much must we send??? */
/*#define TBODYSZ (1+NDIM)*sizeof(double)*/
#define TBODYSZ (offsetof(body, acc)) /*SD--08/25/2005*/

typedef struct {
    double mass;			/* mass of body */
    double pos[NDIM];		/* position of body */
    double vel[NDIM];		/* velocity of body */
    double h;
#ifdef SAVE_ACC
    double acc[NDIM];
    double phi;
#endif
    unsigned int ident;		/* unique? identifier */
    unsigned int type;	/* fill the gap for SDF file */
} outbody, *outbodyptr;

/* This is the descriptor that goes into the SDF header. */

#ifdef SAVE_ACC
#if NDIM==3
#define OUTBODYDESC \
"struct {\n\
    double mass;			/* mass of body */\n\
    double x, y, z;			/* position of body */\n\
    double vx, vy, vz;		/* velocity of body */\n\
    double h;                   /* smoothing length for gravity */\n\
    double ax, ay, az;		/* acceleration */\n\
    double phi;			/* potential */\n\
    unsigned int ident;		/* unique? identifier */\n\
    unsigned int type;       /* fill the gap for SDF file */\n\
}"
#else
#if NDIM==2
#define OUTBODYDESC \
"struct {\n\
    double mass;			/* mass of body */\n\
    double x, y;			/* position of body */\n\
    double vx, vy;		/* velocity of body */\n\
    double h;                   /* smoothing length for gravity */\n\
    double ax, ay;		/* acceleration */\n\
    double phi;			/* potential */\n\
    unsigned int ident;		/* unique? identifier */\n\
    unsigned int type;       /* fill the gap for SDF file */\n\
}"
#else
 # error No case for NDIM
#endif /* NDIM==2 */
#endif /* NDIM==3 */
#else
#if NDIM==3
#define OUTBODYDESC \
"struct {\n\
    double mass;			/* mass of body */\n\
    double x, y, z;			/* position of body */\n\
    double vx, vy, vz;		/* velocity of body */\n\
    double h;                   /* smoothing length for gravity */\n\
    unsigned int ident;		/* unique? identifier */\n\
    unsigned int type;       /* fill the gap for SDF file */\n\
}"
#else
#if NDIM==2
#define OUTBODYDESC \
"struct {\n\
    double mass;			/* mass of body */\n\
    double x, y;			/* position of body */\n\
    double vx, vy;		/* velocity of body */\n\
    double h;                   /* smoothing length for gravity */\n\
    unsigned int ident;		/* unique? identifier */\n\
    unsigned int type;       /* fill the gap for SDF file */\n\
}"
#else
 # error No case for NDIM
#endif /* NDIM==2 */
#endif /* NDIM==3 */
#endif /* SAVE_ACC */

typedef struct {
    double mass;
    double pos[NDIM];
    double bmax, rcrit;
    int daughters;
    double padding_junk;
} cell, *cellptr;


/* This is the intermediate data structure used to construct cofm */
typedef struct{
    double mass;
    double pos[NDIM];
    double massinv;
    double bmax;
    double B2;
    double sz;
    int ndaughters;
} cofmdata;

typedef struct{
    double bmax;
    double pos[NDIM];
    double h;
    double M0;
    double M1[NDIM];
    int isbody;
    double daughters;
    double nterms;
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

/* In main.c */
extern Timer_t StepTot;
extern Timer_t BuildTot;
extern Timer_t FindForcesTm;
extern Counter_t NbodyCnt;
Key_t getkey(const body *);

/* In cofm.c */
void SetupCofm(int MACtype, double tol, double rel_tol);
void CofmFromDaugh(hcellptr hptr, hcellptr daughters[]);
void CellFromCofm(cell *cp, cofmdata *cmp);

/* In print.c */
char *PrintCellContents(const cell *cp);
char *PrintBodyContents(const body *bp);
char *PrintBodyContentsLong(const body *vp);
char *PrintBranch(const cofmdata *cmp);

/* In sphplus.c */
void NbodyPlusNbody(void **btabp, int *nobj, body *btab2, int nobj2);

/* In mac.c */
extern Timer_t GravTm, MACTm;
extern Counter_t CCInt, CBInt, BCInt, BBInt;
extern Counter_t CCIntRej;
extern Counter_t TranslateCnt;

void SetTol(double tol, double frac_tol, double newton_const, double eps, 
	    int gnobj);
void InheritSink(const Sink *from, Sink *to, hcell *pp);
void DLRcritMAC(Sink *sink, const hcell **source, int *result, int n);
void RcritMAC(Sink *sink, const hcell **source, int *result, int n);
void SetGravOffset(double *off);
void UnSetGravOffset(void);
void InheritSinkNlogN(const Sink *from, Sink *to, hcell *pp);

/* In grav.c */
void do_grav(const double *p, const double *end, const double *pos0, 
	     double *mass0, double *acc0, double *phi0, 
	     const double *eps2p, int *ncut);
void update_point_mass(body *btab, int nobj, body *p, double smooth2, 
		       double newt);

/* In sph.c */
void do_SPHgrav(const double *p, const double *end, const double *pos0, 
		double *mass0, double *acc0, double *phi0, 
		const double *eps2p, int *ncut);

#endif


