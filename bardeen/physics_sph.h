/*
 * Copyright 1996 Michael S. Warren and John K. Salmon.  All Rights Reserved.
 */

#include "tree.h"
#include "key.h"
#include "timers.h"

#define NDIM 3
#define SPH_SAVE_ACC
#define POS_IS_DOUBLE

/* Some physical constants, in cgs units */
#define A_COEFF (1.043565e-17)
#define C_LIGHT (3.424758e+02)
#define KES_COEFF (1.043946e+02)
#define KFF_COEFF (1.591470e+12)
#define K_BOLTZ (9.059183e-66)
#define MH (8.411685e-58)

typedef struct {
#ifdef POS_IS_DOUBLE
  /* double first for alignment */
    double pos[NDIM];		/* position of body */
    float mass;			/* mass of body */
#else
    float mass;			/* mass of body */
    float pos[NDIM];		/* position of body */
#endif
    float vel[NDIM];		/* velocity of body */
    float h;			/* smoothing length */
    float rho;			/* density */
    float pr;			/* pressure */
    float vsound;		/* sound speed */
    float rho_est;		/* estimated density */
    float u;			/* internal energy */
    float temp;                 /* temperature, used to enforce LTE */
    float du;                   /* change in internal energy this timestep */
    float dt_next;
    /* Things declared above this line are communicated between processors */
    /* so they can be used in in the loop over nbrs in FindRho and ForceSPH */
    /* Don't add anything above this line unless you fix TBODYSZ */
    float acc[NDIM];
    float grav_acc[NDIM];
    float acc_last[NDIM];
    /* Do these need to go between nodes?  Can things above come down here? */
    float u_r;                  /* radiation energy density */
    float du_r;                 /* change in u_r this timestep */
    float D;                    /* Diffusion coefficient */
    float phi;
    Key_t key;
    unsigned int ident;
    float nterms;
    float grav_nterms;
    float lvel[NDIM];
    float drho_dt;
#ifdef POS_IS_DOUBLE
    double pos_last[NDIM];
#else
    float pos_last[NDIM];
#endif
    float hdot;
    float udot;
    float udot_last;
    unsigned int nbrs;
    float tacc;
    float dt;
    float min_nbr_dt;
    unsigned int windid;
} SPHbody;


/* windbody and WINDOUTBODYDESC need to be padded to a double boundary for
   correct alignment in memory and on disk */
typedef struct {
#ifdef POS_IS_DOUBLE
    double pos[NDIM];
#else
    float pos[NDIM];
#endif
    float vel[NDIM];
    float rhowind;
    float vwind;
    float uwind;
    unsigned int ident;
    int dummy;
} windbody;

typedef struct {		/* don't need all of this info */
    float grav_acc[NDIM];
    float phi;
    int grav_nterms;
    int ident;
    Key_t key;
} accbody;

/* When we send a body from node to node, how much must we send??? */
/*  #define SPHTBODYSZ (8+2*NDIM)*sizeof(float) */
#define SPHTBODYSZ offsetof(SPHbody, acc)

/* If you add anything to the outbody structure, make sure to add an */
/* assignment to the Output routine */
typedef struct {
#ifdef POS_IS_DOUBLE
    double pos[NDIM];		/* position of body */
#else
    float pos[NDIM];		/* position of body */
#endif
    float mass;			/* mass of body */
    float vel[NDIM];		/* velocity of body */
    float u;
    float h;
    float rho;
    float drho_dt;
    float udot;
#ifdef SPH_SAVE_ACC
    float acc[NDIM];
    float acc_last[NDIM];
    float phi;
    float dt;
#endif
    unsigned int nbrs; 
    unsigned int ident;		/* unique? identifier */
    unsigned int windid;
} SPHoutbody;

typedef struct {
#ifdef POS_IS_DOUBLE
    double pos[NDIM];		/* position of body */
#else
    float pos[NDIM];		/* position of body */
#endif
    float mass;			/* mass of body */
    float vel[NDIM];		/* velocity of body */
    float u;
    float h;
    float rho;
    unsigned int nbrs; 
    unsigned int ident;		/* unique? identifier */
    unsigned int windid;
} SPHshortoutbody;

/* This is the descriptor that goes into the SDF header. */

#ifdef SPH_SAVE_ACC
#define SPHOUTBODYDESC \
"struct {\n\
    double x, y, z;		/* position of body */\n\
    float mass;			/* mass of body */\n\
    float vx, vy, vz;		/* velocity of body */\n\
    float u;			/* internal energy */\n\
    float h;			/* smoothing length */\n\
    float rho;			/* density */\n\
    float drho_dt;              /* time derivative of rho */\n\
    float udot;			/* time derivative of u */\n\
    float ax, ay, az;		/* acceleration */\n\
    float lax, lay, laz;	/* acceleration at tpos-dt */\n\
    float phi;			/* potential */\n\
    float idt;			/* timestep */\n\
    unsigned int nbrs;          /* number of neighbors */\n\
    unsigned int ident;		/* unique identifier */\n\
    unsigned int windid;        /* wind id */\n\
}"
#define SPHSHORTOUTBODYDESC \
"struct {\n\
    double x, y, z;           /* position of body */\n\
    float mass;                       /* mass of body */\n\
    float vx, vy, vz;         /* velocity of body */\n\
    float u;                  /* internal energy */\n\
    float h;                  /* smoothing length */\n\
    float rho;                        /* density */\n\
    unsigned int nbrs;          /* number of neighbors */\n\
    unsigned int ident;               /* unique identifier */\n\
    unsigned int windid;        /* wind id */\n\
}"
#else
#define SPHOUTBODYDESC \
"struct {\n\
    double x, y, z;		/* position of body */\n\
    float mass;			/* mass of body */\n\
    float vx, vy, vz;		/* velocity of body */\n\
    float u;			/* internal energy */\n\
    float h;			/* smoothing length */\n\
    float rho;			/* density */\n\
    float drho_dt;              /* time derivative of rho */\n\
    float udot;			/* time derivative of u */\n\
    unsigned int nbrs;          /* number of neighbors */\n\
    unsigned int ident;		/* unique identifier */\n\
    unsigned int windid;        /* wind id */\n\
}"
#define SPHSHORTOUTBODYDESC \
"struct {\n\
    double x, y, z;           /* position of body */\n\
    float mass;                       /* mass of body */\n\
    float vx, vy, vz;         /* velocity of body */\n\
    float u;                  /* internal energy */\n\
    float h;                  /* smoothing length */\n\
    float rho;                        /* density */\n\
    unsigned int nbrs;          /* number of neighbors */\n\
    unsigned int ident;               /* unique identifier */\n\
    unsigned int windid;        /* wind id */\n\
}"
#endif /* SPH_SAVE_ACC */

#define WINDOUTBODYDESC \
"struct {\n\
	double xwind, ywind, zwind;\n\
	float vxwind, vywind, vzwind;\n\
	float rhowind;\n\
	float vwind;\n\
	float uwind;\n\
	unsigned int identwind;\n\
        int dummy;\n\
}"

typedef struct {
    float mass;
    float pos[NDIM];
    float bmax, rcrit;
    int daughters;
    float lap;
} SPHcell;


/* This is the intermediate data structure used to construct cofm */
typedef struct{
    float mass;
    float pos[NDIM];
    float massinv;
    float bmax;
    float B2;
    float sz;
    float lap;
    int ndaughters;
} SPHcofmdata;

typedef struct{
    float extent;
    float pos[NDIM];
    float vel[NDIM];
    float rho;
    float pr;
    float rho_est;
    float vsound;
    float u;
    float temp;
    float du;
    float u_r;
    float du_r;
    float D;
    float mass;
    float drho_dt;
    float udot;
    float M1[NDIM];
    float lvel[NDIM];
    float h;
    int isbody;
    int nbrs;
    unsigned int nterms;
    int interactions;
    float min_nbr_dt;
} SinkSPH;

/* In main.c */
int SPH_need_update(const SPHbody *p);

/* In physics_generic.c */
void CellCorner(Key_t key, float *corner, float *size);

/* In physics_sph.c */
/* There are various void * decls here, since we don't want to have body *s */
void SPHFindBbox(SPHbody *bp, int n, float *rmin, float *rmax);
void SPHFixKeys(SPHbody *btab, int nobj, Key_t (*func)(const void *));
Key_t SPHGetKey(const void *p);
float SPHGetCost(const SPHbody *p);
Key_t SPHGetKeyFromStruct(const SPHbody *p);
void SPHFixId(SPHbody *btab, int nobj, int gnobj);
void SPHFixNterms(SPHbody *btab, int nobj);
Key_t accbodyGetKey(const void *ptr);
Key_t SPHOutIdentKey(const SPHoutbody *bp);
Key_t SPHShortOutIdentKey(const SPHshortoutbody *bp);

/* In sphcofm.c */
void SPHSetupCofm(int MACtype, float tol, float rel_tol);
void SPHCofmFromDaugh(hcellptr hptr, hcellptr daughters[]);
void SPHCellFromCofm(SPHcell *cp, SPHcofmdata *cmp);

/* In sphprint.c */
char *PrintSPHCellContents(const SPHcell *cp);
char *PrintSPHBodyContents(const SPHbody *bp);
char *PrintSPHBodyContentsLong(const SPHbody *vp);
char *PrintSPHBranch(const SPHcofmdata *cmp);

/* In sph.c */
void SetSPH(float visc_alpha, float visc_beta, float visc_epsilon, 
	    float heat_f1, float eos_gamma, int gnobj,  
	    void bfunc(), void cfunc());
void SPHgate(SinkSPH *sink, hcell **src_vec, int *result, int n);
void nbrMAC(SinkSPH *sink, hcell **src_vec, int *result, int n);
void macRho(SinkSPH *sink, hcell **source, int *result, int n);
void macSPH(SinkSPH *sink, hcell **source, int *result, int n);
void InheritSPH(const SinkSPH *from, SinkSPH *to, hcell *pp);
void update_final(SPHbody *btab, int nobj, float dt, int *limit_high, int *limit_low);
void update_intermediate(SPHbody *btab, int nobj, float dt_last, int flag, int *limit);
void SPH_setup(int dim);
void SetSPHOffset(float *off, float *voff);
void UnSetSPHOffset(void);
void update_point_SPHmass(SPHbody *btab, int nobj, void *p, float smooth2, float newt);
void update_point_SPHmass2(SPHbody *btab, int nobj, float smooth2, float newt, float mass);
void update_point_SPHmass3(SPHbody *btab, int nobj, float smooth2, float newt, float mass, float b);
void update_point_mass_bardeen(SPHbody *btab, int nobj, float G, float mass, float rplus, float S);


/* In sphinit.c */
void *DarkRead(char *name, void *csdfp, void **btabp, int *gnobjp, int *nobjp, int set_id, int setpvel);
void *SPHRead(char *name, void *csdfp, SPHbody **btabp, int *gnobjp, int *nobjp, int set_id, int setpvel, float new_h, float new_u);
void SPHTestData(void *csdfp, SPHbody **btabp, int *gnobjp, int *nobjp, int periodic);
void *InitRead(char *name, void *csdfp, void **btabp, int *gnobjp, int *nobjp, 
	 SPHbody **SPHbtabp, int *SPHgnobjp, int *SPHnobjp, 
	 int set_id, int setpvel, float new_h, float new_u);
void DarkSPHTestData(void *csdfp, void **btabp, int *gnobjp, int *nobjp, 
		     SPHbody **SPHbtabp, int *SPHgnobjp, int *SPHnobjp, 
		     int periodic);
void *WindRead(char *name, void *csdfp, windbody **btabp, int *gnobjp, 
	       int *nobjp);

/* In sphplus.c */
void GravPlusSPH(void **btab, int *nobj, SPHbody *SPHbtab, int SPHnobj);
void GravMinusSPH(void **btab, int *nobj, accbody **atab, int *anobj);

/* In eos.c */
double uvst(double t);
double duvst(double t);

/* In newtraph.c */
float newtraph(double xl, double xr, double prec, double (*f)(double x), 
	     double (*df)(double x));

