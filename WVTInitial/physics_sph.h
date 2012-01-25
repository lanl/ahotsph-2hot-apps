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
    double pos[NDIM];		/* position of body */
    double mass;		/* mass of body */
    double vel[NDIM];		/* velocity of body */
    double h;			/* smoothing length */
    double rho;			/* density */
    double pr;			/* pressure */
    double vsound;		/* sound speed */
    double rho_est;		/* estimated density */
    double u;			/* internal energy */
    double temp;                /* temperature, used to enforce LTE */
    double du;                  /* change in internal energy this timestep */
    double dt_next;
    /* Things declared above this line are communicated between processors */
    /* so they can be used in the loop over nbrs in FindRho and ForceSPH */
    /* Don't add anything above this line unless you fix TBODYSZ */
    double acc[NDIM];
/*     double grav_acc[NDIM]; */
/*     double acc_last[NDIM]; */
    /* Do these need to go between nodes?  Can things above come down here? */
    double u_r;                  /* radiation energy density */
    double du_r;                 /* change in u_r this timestep */
/*     double D;                    /\* Diffusion coefficient *\/ */
/*     double phi; */
    double grav_mass;      /* normally = mass, different for dual particles*/
    Key_t key;
    unsigned int ident;
    double nterms;
    double grav_nterms;
    double lvel[NDIM];
    double drho_dt;
    double pos_last[NDIM];
    double hdot;
    double udot;
    double udot_last;
    unsigned int nbrs;
/*     double tacc; */
/*     double dt; */
    double min_nbr_dt;
    unsigned int windid;
    unsigned int type;
} SPHbody;

/* typedef struct { */
/*     double pos[NDIM];		/\* position of body *\/ */
/*     double mass;		/\* mass of body *\/ */
/*     double vel[NDIM];		/\* velocity of body *\/ */
/*     double h;			/\* smoothing length *\/ */
/*     double rho;			/\* density *\/ */
/*     double pr;			/\* pressure *\/ */
/*     double vsound;		/\* sound speed *\/ */
/*     double rho_est;		/\* estimated density *\/ */
/*     double u;			/\* internal energy *\/ */
/*     double temp;                /\* temperature, used to enforce LTE *\/ */
/*     double du;                  /\* change in internal energy this timestep *\/ */
/*     double dt_next; */
/* } MYTREESPH; */

/* windbody and WINDOUTBODYDESC need to be padded to a double boundary for
   correct alignment in memory and on disk */
/* typedef struct { */
/* #ifdef POS_IS_DOUBLE */
/*     double pos[NDIM]; */
/* #else */
/*     double pos[NDIM]; */
/* #endif */
/*     double vel[NDIM]; */
/*     double rhowind; */
/*     double vwind; */
/*     double uwind; */
/*     unsigned int ident; */
/*     int dummy; */
/* } windbody; */

typedef struct {		/* don't need all of this info */
    double grav_acc[NDIM];
    double phi;
    int grav_nterms;
    int ident;
    Key_t key;
} accbody;

/* When we send a body from node to node, how much must we send??? */
/*  #define SPHTBODYSZ (8+2*NDIM)*sizeof(double) */
#define SPHTBODYSZ offsetof(SPHbody, acc)

/* If you add anything to the outbody structure, make sure to add an */
/* assignment to the Output routine */

typedef struct {
    double pos[NDIM];		/* position of body */
    double mass;	   /* mass of body */
    double vel[NDIM];     /* velocity of body */
    double u;      	   /* specific energy of body*/
    double h;      	   /* smoothing length of body */
    double rho;            /* density of body */
    double drho_dt;        /* drho/dt of body */
    double udot;           /* du/dt of body */
    double temp;           /* temperature of body */
#ifdef SPH_SAVE_ACC
    double acc[NDIM];     /* acceleration of body */
    double acc_last[NDIM];  /* last acceleration of body */
    double grav_acc[NDIM];  /* grav acceleration of body */
    double grav_mass;      /* normally = mass, for dual particles different */
    double phi;            /* potential at body location */
    double dt;             /* timestep of body */
#endif
    unsigned int nbrs;     /* number of neighbors */
    unsigned int ident;	   /* unique identifier */
    unsigned int windid;   /* wind id */
    unsigned int type;     /* to fill the double block (4int=1double)*/
} SPHoutbody;

typedef struct {
    double pos[NDIM];		/* position of body */
    float mass;	   /* mass of body */
    float vel[NDIM];     /* velocity of body */
    float u;      	   /* specific energy of body*/
    float h;      	   /* smoothing length of body */
    float rho;            /* density of body */
    float drho_dt;        /* drho/dt of body */
    float udot;           /* du/dt of body */
    float temp;           /* temperature of body */
#ifdef SPH_SAVE_ACC
    float acc[NDIM];     /* acceleration of body */
    float acc_last[NDIM];  /* last acceleration of body */
    float phi;            /* potential at body location */
    float dt;             /* timestep of body */
#endif
    unsigned int nbrs;     /* number of neighbors */
    unsigned int ident;	   /* unique identifier */
    unsigned int windid;   /* wind id */
    unsigned int padding;  /* padding for structure alignment */
} SPHfloatoutbody;

typedef struct {
    double pos[NDIM];		/* position of body */
    double mass;			/* mass of body */
    double vel[NDIM];		/* velocity of body */
    double u;
    double h;
    double rho;
    unsigned int nbrs; 
    unsigned int ident;		/* unique? identifier */
    unsigned int windid;
    unsigned int type;  /* to fill the double block (4int=1double)*/
} SPHshortoutbody;

typedef struct {
    float pos[NDIM];		/* position of body */
    float mass;			/* mass of body */
    float vel[NDIM];		/* velocity of body */
    float u;
    float h;
    float rho;
    unsigned int nbrs; 
    unsigned int ident;		/* unique? identifier */
    unsigned int windid;
    unsigned int type;  /* to fill the double block (4int=1double)*/
} SPHshortfloatoutbody;

/* This is the descriptor that goes into the SDF header. */

#ifdef SPH_SAVE_ACC
#define SPHOUTBODYDESC \
"struct {\n\
    double x, y, z;		/* position of body */\n\
    double mass;	   /* mass of body */\n\
    double vx, vy, vz;     /* velocity of body */\n\
    double u;      	   /* specific energy of body*/\n\
    double h;      	   /* smoothing length of body */\n\
    double rho;            /* density of body */\n\
    double drho_dt;        /* drho/dt of body */\n\
    double udot;           /* du/dt of body */\n\
    double temp;           /* temperature of body */\n\
    double ax, ay, az;     /* acceleration of body */\n\
    double lax, lay, laz;  /* last acceleration of body */\n\
    double gax, gay, gaz;  /* gravity acceleration of body */\n\
    double grav_mass;      /* gravitational mass of body */\n\
    double phi;            /* potential at body location */\n\
    double idt;             /* timestep of body */\n\
    unsigned int nbrs;     /* number of neighbors */\n\
    unsigned int ident;	   /* unique identifier */\n\
    unsigned int windid;   /* wind id */\n\
    unsigned int useless;  /* to fill the double block (4int=1double)*/\n\
}"

#define SPHFLOATOUTBODYDESC \
"struct {\n\
    double x, y, z;		/* position of body */\n\
    float mass;	   /* mass of body */\n\
    float vx, vy, vz;     /* velocity of body */\n\
    float u;      	   /* specific energy of body*/\n\
    float h;      	   /* smoothing length of body */\n\
    float rho;            /* density of body */\n\
    float drho_dt;        /* drho/dt of body */\n\
    float udot;           /* du/dt of body */\n\
    float temp;           /* temperature of body */\n\
    float ax, ay, az;     /* acceleration of body */\n\
    float lax, lay, laz;  /* last acceleration of body */\n\
    float phi;            /* potential at body location */\n\
    float idt;             /* timestep of body */\n\
    unsigned int nbrs;     /* number of neighbors */\n\
    unsigned int ident;	   /* unique identifier */\n\
    unsigned int windid;   /* wind id */\n\
    unsigned int padding;  /* padding for structure alighnment */\n\
}"


#define SPHSHORTOUTBODYDESC \
"struct {\n\
    double x, y, z;           /* position of body */\n\
    double mass;                       /* mass of body */\n\
    double vx, vy, vz;         /* velocity of body */\n\
    double u;                  /* internal energy */\n\
    double h;                  /* smoothing length */\n\
    double rho;                        /* density */\n\
    unsigned int nbrs;          /* number of neighbors */\n\
    unsigned int ident;               /* unique identifier */\n\
    unsigned int windid;        /* wind id */\n\
    unsigned int useless;  /* to fill the double block (4int=1double)*/\n\
}"

#define SPHSHORTFLOATOUTBODYDESC \
"struct {\n\
    float x, y, z;           /* position of body */\n\
    float mass;                       /* mass of body */\n\
    float vx, vy, vz;         /* velocity of body */\n\
    float u;                  /* internal energy */\n\
    float h;                  /* smoothing length */\n\
    float rho;                        /* density */\n\
    unsigned int nbrs;          /* number of neighbors */\n\
    unsigned int ident;               /* unique identifier */\n\
    unsigned int windid;        /* wind id */\n\
    unsigned int useless;  /* to fill the double block (4int=1double)*/\n\
}"

#else
#define SPHOUTBODYDESC \
"struct {\n\
    double x, y, z;	   /* position of body */\n\
    double mass;	   /* mass of body */\n\
    double vx, vy, vz;     /* velocity of body */\n\
    double u;      	   /* specific energy of body*/\n\
    double h;      	   /* smoothing length of body */\n\
    double rho;            /* density of body */\n\
    double drho_dt;              /* time derivative of rho */\n\
    double udot;			/* time derivative of u */\n\
    unsigned int nbrs;     /* number of neighbors */\n\
    unsigned int ident;	   /* unique identifier */\n\
    unsigned int windid;   /* wind id */\n\
    unsigned int useless;  /* to fill the double block (4int=1double)*/\n\
}"

#define SPHFLOATOUTBODYDESC \
"struct {\n\
    double x, y, z;	   /* position of body */\n\
    float mass;	   /* mass of body */\n\
    float vx, vy, vz;     /* velocity of body */\n\
    float u;      	   /* specific energy of body*/\n\
    float h;      	   /* smoothing length of body */\n\
    float rho;            /* density of body */\n\
    float drho_dt;              /* time derivative of rho */\n\
    float udot;			/* time derivative of u */\n\
    unsigned int nbrs;     /* number of neighbors */\n\
    unsigned int ident;	   /* unique identifier */\n\
    unsigned int windid;   /* wind id */\n\
}"


#define SPHSHORTOUTBODYDESC \
"struct {\n\
    double x, y, z;           /* position of body */\n\
    double mass;                       /* mass of body */\n\
    double vx, vy, vz;         /* velocity of body */\n\
    double u;                  /* internal energy */\n\
    double h;                  /* smoothing length */\n\
    double rho;                        /* density */\n\
    unsigned int nbrs;          /* number of neighbors */\n\
    unsigned int ident;               /* unique identifier */\n\
    unsigned int windid;        /* wind id */\n\
    unsigned int useless;  /* to fill the double block (4int=1double)*/\n\
}"

#define SPHSHORTFLOATOUTBODYDESC \
"struct {\n\
    float x, y, z;           /* position of body */\n\
    float mass;                       /* mass of body */\n\
    float vx, vy, vz;         /* velocity of body */\n\
    float u;                  /* internal energy */\n\
    float h;                  /* smoothing length */\n\
    float rho;                        /* density */\n\
    unsigned int nbrs;          /* number of neighbors */\n\
    unsigned int ident;               /* unique identifier */\n\
    unsigned int windid;        /* wind id */\n\
    unsigned int useless;  /* to fill the double block (4int=1double)*/\n\
}"

#endif /* SPH_SAVE_ACC */

/* #define WINDOUTBODYDESC \ */
/* "struct {\n\ */
/* 	double xwind, ywind, zwind;\n\ */
/* 	double vxwind, vywind, vzwind;\n\ */
/* 	double rhowind;\n\ */
/* 	double vwind;\n\ */
/* 	double uwind;\n\ */
/* 	unsigned int identwind;\n\ */
/*         int dummy;\n\ */
/* }" */

typedef struct {
    double mass;
    double pos[NDIM];
    double bmax, rcrit;
    int daughters;
    double lap;
} SPHcell;


/* This is the intermediate data structure used to construct cofm */
typedef struct{
    double mass;
    double pos[NDIM];
    double massinv;
    double bmax;
    double B2;
    double sz;
    double lap;
    int ndaughters;
} SPHcofmdata;

typedef struct{
    double extent;
    double pos[NDIM];
    double vel[NDIM];
    double rho;
    double pr;
    double rho_est;
    double vsound;
    double u;
    double temp;
    double du;
    double u_r;
    double du_r;
    double D;
    double mass;
    double drho_dt;
    double udot;
    double M1[NDIM];
    double lvel[NDIM];
    double h;
    int isbody;
    int nbrs;
    unsigned int nterms;
    int interactions;
    double min_nbr_dt;
  unsigned int ident; /*SD, 08-18-2005*/
} SinkSPH;

/* In main.c */
int SPH_need_update(const SPHbody *p);

/* In physics_generic.c */
void CellCorner(Key_t key, double *corner, double *size);

/* In physics_sph.c */
/* There are various void * decls here, since we don't want to have body *s */
void SPHFindBbox(SPHbody *bp, int n, double *rmin, double *rmax);
void SPHFixKeys(SPHbody *btab, int nobj, Key_t (*func)(const void *));
Key_t SPHGetKey(const void *p);
double SPHGetCost(const SPHbody *p);
Key_t SPHGetKeyFromStruct(const SPHbody *p);
void SPHFixId(SPHbody *btab, int nobj, int gnobj);
void SPHFixNterms(SPHbody *btab, int nobj);
Key_t accbodyGetKey(const void *ptr);
Key_t SPHOutIdentKey(const SPHoutbody *bp);
/* Key_t TESTOutIdentKey(const TESToutbody *bp); */
Key_t SPHShortOutIdentKey(const SPHshortoutbody *bp);


/* In sphcofm.c */
void SPHSetupCofm(int MACtype, double tol, double rel_tol);
void SPHCofmFromDaugh(hcellptr hptr, hcellptr daughters[]);
void SPHCellFromCofm(SPHcell *cp, SPHcofmdata *cmp);

/* In sphprint.c */
char *PrintSPHCellContents(const SPHcell *cp);
char *PrintSPHBodyContents(const SPHbody *bp);
char *PrintSPHBodyContentsLong(const SPHbody *vp);
char *PrintSPHBranch(const SPHcofmdata *cmp);

/* In sph.c */
void SetSPH(double visc_alpha, double visc_beta, double visc_epsilon, 
	    double heat_f1, double eos_gamma, int gnobj,  
	    void bfunc(), void cfunc());
void SPHgate(SinkSPH *sink, hcell **src_vec, int *result, int n);
void nbrMAC(SinkSPH *sink, hcell **src_vec, int *result, int n);
void macRho(SinkSPH *sink, hcell **source, int *result, int n);
void macSPH(SinkSPH *sink, hcell **source, int *result, int n);
void InheritSPH(const SinkSPH *from, SinkSPH *to, hcell *pp);
void update_final(SPHbody *btab, int nobj, double dt, int *limit_high, 
		  int *limit_low);
void update_intermediate(SPHbody *btab, int nobj, double dt_last, int flag, 
			 int *limit);
/*void SPH_setup(int dim);*/
void SPH_setup(int dim, int ncoef1, double *wcoef1, int ncoef2, 
	       double *wcoef2);
void SetSPHOffset(double *off, double *voff);
void UnSetSPHOffset(void);
void update_point_SPHmass(SPHbody *btab, int nobj, void *p, double smooth2, 
			  double newt);
void update_point_SPHmass2(SPHbody *btab, int nobj, double smooth2, 
			   double newt, double mass);
void update_point_SPHmass3(SPHbody *btab, int nobj, double smooth2, 
			   double newt, double mass, double b);

/* In sphinit.c */
void *DarkRead(char *name, void *csdfp, void **btabp, int *gnobjp, 
	       int *nobjp, int set_id, int setpvel);
void *SPHRead(char *name, void *csdfp, SPHbody **btabp, int *gnobjp, 
	      int *nobjp, int set_id, int setpvel, double new_h, double new_u);
void SPHTestData(void *csdfp, SPHbody **btabp, int *gnobjp, int *nobjp, 
		 int periodic);
void *InitRead(char *name, void *csdfp, void **btabp, int *gnobjp, int *nobjp, 
	       SPHbody **SPHbtabp, int *SPHgnobjp, int *SPHnobjp, 
	       int set_id, int setpvel, double new_h, double new_u);
void DarkSPHTestData(void *csdfp, void **btabp, int *gnobjp, int *nobjp, 
		     SPHbody **SPHbtabp, int *SPHgnobjp, int *SPHnobjp, 
		     int periodic);
/* void *WindRead(char *name, void *csdfp, windbody **btabp, int *gnobjp,  */
/* 	       int *nobjp); */

/* In sphplus.c */
void GravPlusSPH(void **btab, int *nobj, SPHbody *SPHbtab, int SPHnobj);
void GravMinusSPH(void **btab, int *nobj, accbody **atab, int *anobj);
void SPHPlusSPH(void **btabp, int *nobj, SPHbody *SPHbtab, int SPHnobj);


/* In eos.c */
double uvst(double t);
double duvst(double t);

/* In newtraph.c */
double newtraph(double xl, double xr, double prec, double (*f)(double x), 
		double (*df)(double x));

