/*
 * Copyright 1993 Michael S. Warren and John K. Salmon.  All Rights Reserved.
 */

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
    float mass;			        /* mass of body */
    float pos[NDIM];		        /* position of body */
    float vel[NDIM];		        /* velocity of body */
    float h;			           /* smoothing length */
    float rho;			           /* density */
    float pr;			           /* pressure */
    float vsound;		           /* sound speed */
    float rho_est;	           /* estimated density */
    float gr_mass;              /* mass/sqrt(hdet) */
    Key_t key;
    /* Things declared above this line are communicated between processors */
    /* so they can be used in in the loop over nbrs in FindRho and ForceSPH */
    /* Don't add anything above this line unless you fix TBODYSZ */
    float acc[NDIM];            /* forces */
    float phi;                  /* self-gravity potential */
    float u;		   	        /* internal energy */
    float udot;
    float vel_last[NDIM];       
    float force_last[NDIM];       
    float udot_last;            
    float drho_dt;               
    float hdot;
    float acc_last;
    unsigned int ident;
    unsigned int nterms;
    unsigned int nbrs;
    /* GR */
    float vflow[NDIM];          /* flow velocity */
    float mom[NDIM+1];		     /* momentum */
    float enth;                 /* enthalpy */
    float gama;                 /* Lorentz factor */
    float gama_last;            /* Old gama */
    float alfa;                 /* redshift */
    float gxx;                  /* metric */
    float gyy;
    float gzz;
    float gxy;
    float gxz;
    float gyz;
    float gxt;
    float gyt;
    float gzt;
    float gtt;
    float guxx;                 /* metric inverse */
    float guyy;
    float guzz;
    float guxy;
    float guxz;
    float guyz;
    float guxt;
    float guyt;
    float guzt;
    float gutt;
    float hdetx;
    float hdety;
    float hdetz;
    /* end GR */
} body, *bodyptr;

/* When we send a body from node to node, how much must we send??? */
#define TBODYSZ ((7+2*NDIM)*sizeof(float)+sizeof(Key_t))

typedef struct {
    float mass;			/* mass of body */
    float pos[NDIM];		/* position of body */
    float vel[NDIM];		/* velocity of body */
    float u;			/* internal energy */
    float h;
    float rho;
    float phi;
    unsigned int nbrs;
    unsigned int ident;		/* unique? identifier */
    /* GR */
    float mom[NDIM+1];
    float gama;
    float enth;
    /* end GR */
} outbody, *outbodyptr;

/* This is the descriptor that goes into the SDF header. */

#if NDIM ==3
#define OUTBODYDESC \
"struct {\n\
    float mass;			/* mass of body */\n\
    float x, y, z;			/* position of body */\n\
    float vx, vy, vz;		/* velocity of body */\n\
    float u;			/* internal energy */\n\
    float h;			/* smoothing length */\n\
    float rho;			/* density */\n\
    float phi;			/* grav potential */\n\
    unsigned int nbrs;		/* number of neighbors */\n\
    unsigned int ident;		/* unique? identifier */\n\
    float sx, sy, sz, st; \n\
    float gama; \n\
    float enthalpy; \n\
}"
#else
#if NDIM==2
#define OUTBODYDESC \
"struct {\n\
    float mass;			/* mass of body */\n\
    float x, y;			/* position of body */\n\
    float vx, vy;		/* velocity of body */\n\
    float u;			/* internal energy */\n\
    float h;			/* smoothing length */\n\
    float rho;			/* density */\n\
    float phi;			/* grav potential */\n\
    unsigned int nbrs;		/* number of neighbors */\n\
    unsigned int ident;		/* unique? identifier */\n\
    float sx, sy, st; \n\
    float gama; \n\
    float enthalpy; \n\
}"
#else
#error No case for NDIM
#endif
#endif

typedef struct {
    float mass;
    float pos[NDIM];
    float rcrit, B2, B3, bmax;
    float lap;
    unsigned int daughters;
    float acc_last_max;
} cell, *cellptr;

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

typedef struct{
    float bmax;
    float pos[NDIM];
    float M0;
    float M1[NDIM];
    float acc_last_max;
    int isbody;
    moment M2;
    int daughters;
    unsigned int nterms;
    int interactions;
} Sink;

typedef struct{
    float extent;
    float pos[NDIM];
    float vel[NDIM];
    float rho;
    float pr;
    float rho_est;
    float vsound;
    float mass;
    float drho_dt;
    float udot;
    float M1[NDIM];
    float h;
	 float alfa;
    int isbody;
    int nbrs;
    unsigned int nterms;
    int interactions;
} SinkSPH;

/* Now for Order N */

/* This is the intermediate data structure used to construct cofm */
typedef struct{
    float mass;
    float pos[NDIM];
    float B2;
    float bmax;
    float acc_last_max;
    float lap;
    float massinv;
    int ndaughters;
} cofmdata;

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

/* In walk_n.c */
extern Counter_t MACccCnt;
extern Counter_t MACccPass;
extern Counter_t MACcbCnt;
extern Counter_t MACcbPass;
extern Counter_t MACbcCnt;
extern Counter_t MACbcPass;
extern Counter_t WalkCnt;
extern Counter_t DeferCnt;
extern Counter_t NobjCnt;
extern Counter_t NtermsCnt;
extern Timer_t ImbalTm;
extern Timer_t GravTm;
void
FindForces(tree_t *tp, float GNewt, float eps,
	int check_parents,  void (*init_physdata)(body *, cell *),
	void (*CCinteract)(body **pp, body **end, const float *pos0, 
		      float *mass0, float *phi0, float *acc0,
		      moment *qpole0, const float *eps2p, int *ncut,
		      int *tot_interact),
	void (*BCinteract)(body *p, body **pp, body **end),
	float (*bmaxf)(cell *), 
      int (*CCmac)(float x0, float x1, float x2, cell *, float, Stk *, int *),
      int (*BCmac)(float x0, float x1, float x2, cell *, float, Stk *, int *),
      int (*CBmac)(float x0, float x1, float x2, cell *, float, Stk *, int *, int));

/* In grav_n.c */
extern Counter_t GravCnt;
extern Counter_t CCInt;
extern Counter_t TranslateCnt;

void do_body(body *b, body **pp, body **end);
void FindRho(body *p, body **nbr_list, body **end_list);
void forceSPH(body *p, body **nbr_list, body **end_list);

void Translate(const float *new_pos,  const float *old_pos,
	  float phi, const float *acc, const moment *qpole, 
	  float phi_old, const float *acc_old,
	  float M0_old, const float *M1_old, const moment *M2_old,
	  float *M0_new, float *M1_new, moment *M2_new,
	  float *phi_new, float *acc_new);

/* In physics_n.c */
char *PrintCellContents(const cell *cp);
char *PrintBodyContents(const body *bp);
char *PrintBodyContentsLong(const body *vp);
char *PrintBranch(const cofmdata *cmp);

/* In rdtest.c */
void RdTest(body **btabp, int gnobj, int *nobjp, int seed, int cencon);

void SetTol(float tol, float frac_tol, float newton_const, float eps, int gnobj);
void Unifiedmacv(Sink *sink, const hcell **source, int *result, int n);
void Fracmacv(Sink *sink, const hcell **source, int *result, int n);
void Lowestmacv(Sink *sink, const hcell **source, int *result, int n);
void Nlognmacv(Sink *sink, const hcell **src_vec, int *result, int n);
void InheritSink(const Sink *from, Sink *to, hcell *pp);
void SPHgate(SinkSPH *sink, hcell **src_vec, int *result, int n);
void nbrMAC(SinkSPH *sink, hcell **src_vec, int *result, int n);
void macRho(SinkSPH *sink, hcell **source, int *result, int n);
void macSPH(SinkSPH *sink, hcell **source, int *result, int n);
void InheritSPH(const SinkSPH *from, SinkSPH *to, hcell *pp);

void SPH_setup(int dim);
void SetSPH(float visc_alpha, float visc_beta, float eos_gamma, int gnobj, 
       void bfunc(SinkSPH *sink, hcell **source, int *result, int n), 
       void cfunc(SinkSPH *sink, hcell **source, int *result, int n));
/* GR */
void setup_metric(int kerr_flag, float hole_mass, float kerr_ang_mom);
void get_metric(body *btab, int nobj);
void add_gr(body *btab, int nobj);
void initial_cond(body *btab, int nobj,
       float xx0, float yy0, float zz0,
       float vx0, float vy0, float vz0,
       float bhmass, float Gamma);
/* Don */
void grav_rad(body *btab, int nobj, float *hp, float *hx);
int df(int x, int y);
/* end Don */
/* end GR */
/* In shrink.c */
void ShrinkBtab(body **btabp, int *nobj, float r_limit);

