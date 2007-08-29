/*
 * Copyright 1996 Michael S. Warren and John K. Salmon.  All Rights Reserved.
 */

#include "tree.h"
#include "key.h"
#include "timers.h"
#include "ndim.h"

#define SPH_SAVE_ACC
/* #define POS_IS_DOUBLE */

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
    float abar;
    float temp;
    float ye;
    float xp;
    float xn;
    float xmu;
    int ifleos;
    float dt_next;
    float gshift;
    float r;
    /* Things declared above this line are communicated between processors */
    /* so they can be used in in the loop over nbrs in FindRho and ForceSPH */
    /* Don't add anything above this line unless you fix TBODYSZ */
    float acc[NDIM];
    float grav_acc[NDIM];
    float acc_last[NDIM];
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
    float udot2;
    float udot2_last;
    unsigned int nbrs;
    float tacc;
    float dt;
    float min_nbr_dt;
    float u2;
    float xpf;
    float p2;
    float p3;
    float p4;
    float dye;
    float temprev;
    float rhoprev;
    float xpprev;
    float xnprev;
    float yeprev;
    float ufreez;
    float dtbrn;
    float xhe4n;
    float xc12n;
    float xo16n;
    float xne20n;
    float xmg24n;
    float xsi28n;
    float xni56n;
    float xco56n;
    float xfe56n;
    float eta;
    float xmuhat;
    float xmue;
    float prg;
    int bghost;
    void *ireal;
    short ebeta, pbeta;
} SPHbody;

typedef struct {		/* don't need all of this info */
    float grav_acc[NDIM];
    float phi;
    int grav_nterms;
    int ident;
    Key_t key;
} accbody;

/* When we send a body from node to node, how much must we send??? */
#define SPHTBODYSZ offsetof(SPHbody, acc)

/* If you add anything to the outbody structure, make sure to add an */
/* assignment to the Output routine */
typedef struct {
#ifdef POS_IS_DOUBLE
    double pos[NDIM];		/* position of body */
    float mass;			/* mass of body */
#else
    float mass;			/* mass of body */
    float pos[NDIM];		/* position of body */
#endif
    float vel[NDIM];		/* velocity of body */
    float u;
    float h;
    float rho;
    float drho_dt;
    float udot;
    float pr;
    float vsound;
    float temp;
    float ye;
    float xp;
    float xn;
    float u2;
    float abar;
    float ufreez;
    float dtbrn;
    float xhe4n;
    float xc12n;
    float xo16n;
    float xne20n;
    float xmg24n;
    float xsi28n;
    float xni56n;
    float xco56n;
    float xfe56n;
    int  ifleos;
#ifdef SPH_SAVE_ACC
    float acc[NDIM];
    float acc_last[NDIM];
    float phi;
    float dt;
#endif
    unsigned int nbrs;
    unsigned int ident;		/* unique? identifier */
    unsigned int dummy;
} SPHoutbody;

/* This is the descriptor that goes into the SDF header. */

#if NDIM==3
#ifdef SPH_SAVE_ACC
#define SPHOUTBODYDESC \
"struct {\n\
    float mass;			/* mass of body */\n\
    float x, y, z;			/* position of body */\n\
    float vx, vy, vz;		/* velocity of body */\n\
    float u;			/* internal energy */\n\
    float h;			/* smoothing length */\n\
    float rho;			/* density */\n\
    float drho_dt;		/* time derivative of density */\n\
    float udot;			/* time derivative of u */\n\
    float pr;			\n\
    float vsound;		\n\
    float temp;			\n\
    float ye;			\n\
    float xp;			\n\
    float xn;			\n\
    float u2;			\n\
    float abar;			\n\
    float ufreez;		\n\
    float dtbrn;		\n\
    float xhe4n;		\n\
    float xc12n;		\n\
    float xo16n;		\n\
    float xne20n;		\n\
    float xmg24n;		\n\
    float xsi28n;		\n\
    float xni56n;		\n\
    float xco56n;		\n\
    float xfe56n;		\n\
    int  ifleos;		\n\
    float ax, ay, az;		/* acceleration */\n\
    float lax, lay, laz;	/* acceleration at tpos-dt */\n\
    float phi;			/* potential */\n\
    float idt;			/* timestep */\n\
    unsigned int nbrs;		/* number of neighbors */\n\
    unsigned int ident;		/* unique identifier */\n\
    unsigned int dummy;\n\
}"
#else
#define SPHOUTBODYDESC \
"struct {\n\
    float mass;			/* mass of body */\n\
    float x, y, z;		/* position of body */\n\
    float vx, vy, vz;		/* velocity of body */\n\
    float u;			/* internal energy */\n\
    float h;			/* smoothing length */\n\
    float rho;			/* density */\n\
    float drho_dt;		/* time derivative of density */\n\
    float udot;			/* time derivative of u */\n\
    float pr;			\n\
    float temp;			\n\
    float ye;			\n\
    float xp;			\n\
    float xn;			\n\
    float u2;			\n\
    float abar;			\n\
    float ufreez;		\n\
    float dtbrn;		\n\
    float xhe4n;		\n\
    float xc12n;		\n\
    float xo16n;		\n\
    float xne20n;		\n\
    float xmg24n;		\n\
    float xsi28n;		\n\
    float xni56n;		\n\
    float xco56n;		\n\
    float xfe56n;		\n\
    int  ifleos;		\n\
    float vsound;		\n\
    unsigned int nbrs;		/* number of neighbors */\n\
    unsigned int ident;		/* unique identifier */\n\
    unsigned int dummy;\n\
}"
#endif /* SPH_SAVE_ACC */
#else  /* NDIM==2 */
#ifdef SPH_SAVE_ACC
#define SPHOUTBODYDESC \
"struct {\n\
    float mass;			/* mass of body */\n\
    float x, y;			/* position of body */\n\
    float vx, vy;		/* velocity of body */\n\
    float u;			/* internal energy */\n\
    float h;			/* smoothing length */\n\
    float rho;			/* density */\n\
    float drho_dt;		/* time derivative of density */\n\
    float udot;			/* time derivative of u */\n\
    float pr;			\n\
    float vsound;		\n\
    float temp;			\n\
    float ye;			\n\
    float xp;			\n\
    float xn;			\n\
    float u2;			\n\
    float abar;			\n\
    float ufreez;		\n\
    float dtbrn;		\n\
    float xhe4n;		\n\
    float xc12n;		\n\
    float xo16n;		\n\
    float xne20n;		\n\
    float xmg24n;		\n\
    float xsi28n;		\n\
    float xni56n;		\n\
    float xco56n;		\n\
    float xfe56n;		\n\
    int  ifleos;		\n\
    float ax, ay;		/* acceleration */\n\
    float lax, lay;		/* acceleration at tpos-dt */\n\
    float phi;			/* potential */\n\
    float idt;			/* timestep */\n\
    unsigned int nbrs;		/* number of neighbors */\n\
    unsigned int ident;		/* unique identifier */\n\
}"
#else
#define SPHOUTBODYDESC \
"struct {\n\
    float mass;			/* mass of body */\n\
    float x, y;			/* position of body */\n\
    float vx, vy;		/* velocity of body */\n\
    float u;			/* internal energy */\n\
    float h;			/* smoothing length */\n\
    float rho;			/* density */\n\
    float drho_dt;		/* time derivative of density */\n\
    float udot;			/* time derivative of u */\n\
    float pr;			\n\
    float vsound;		\n\
    float temp;			\n\
    float ye;			\n\
    float xp;			\n\
    float xn;			\n\
    float u2;			\n\
    float abar;			\n\
    float ufreez;		\n\
    float dtbrn;		\n\
    float xhe4n;		\n\
    float xc12n;		\n\
    float xo16n;		\n\
    float xne20n;		\n\
    float xmg24n;		\n\
    float xsi28n;		\n\
    float xni56n;		\n\
    float xco56n;		\n\
    float xfe56n;		\n\
    int  ifleos;		\n\
    unsigned int nbrs;		/* number of neighbors */\n\
    unsigned int ident;		/* unique identifier */\n\
}"
#endif
#endif /* NDIM==2 */

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
    float mass;
    float drho_dt;
    float udot;
    float udot2;
    float M1[NDIM];
    float lvel[NDIM];
    float h;
    int isbody;
    int nbrs;
    unsigned int nterms;
    int interactions;
    float min_nbr_dt;
    float dt;
    float r;
    float gshift;
    float xfac;			/* geometrical factor */
} SinkSPH;

typedef struct {
    float pos[NDIM];
    float vel[NDIM];
    float mass;
    float r;
} bndry_t;

/* External Fortran linkage */
#define Fortran(x) x##_
/* GNU Fortran adds two underscores if there is an underscore in the name */
#ifdef __GNUC__
#define Fortran2(x) x##__
#else
#define Fortran2(x) x##_
#endif

/* In main.c */
int SPH_need_update(const SPHbody *p);

/* In physics_generic.c */
void CellCorner(Key_t key, float *corner, float *size);

/* In physics_sph.c */
/* There are various void * decls here, since we don't want to have body *s */
void SPHFindBbox(SPHbody *bp, int n, float *rmin, float *rmax);
void SPHFixKeys(SPHbody *btab, int nobj, Key_t (*func)(const SPHbody *));
Key_t SPHGetKey(const SPHbody *p);
float SPHGetCost(const SPHbody *p);
Key_t SPHGetKeyFromStruct(const SPHbody *p);
void SPHFixId(SPHbody *btab, int nobj, int gnobj);
void SPHFixNterms(SPHbody *btab, int nobj);
Key_t accbodyGetKey(const void *ptr);
Key_t SPHOutIdentKey(const SPHoutbody *bp);

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
void SPHaux(float rinner);
void SPHgate(SinkSPH *sink, hcell **src_vec, int *result, int n);
void nbrMAC(SinkSPH *sink, hcell **src_vec, int *result, int n);
void macRho(SinkSPH *sink, hcell **source, int *result, int n);
void macSPH(SinkSPH *sink, hcell **source, int *result, int n);
void InheritSPH(const SinkSPH *from, SinkSPH *to, hcell *pp);
void update_final(SPHbody *btab, int nobj, float dt, int *limit_high, int *limit_low);
void update_intermediate(SPHbody *btab, int nobj, float dt_last, int flag, int *limit);
void SPH_setup(int dim, int ncoef1, double *coef1, int ncoef2, double *coef2);

void SetSPHOffset(float *off, float *voff);
void UnSetSPHOffset(void);
void SetSPHRotate(float angle);
void UnSetSPHRotate(void);
void update_point_SPHmass(SPHbody *btab, int nobj, void *p, float smooth2, float newt);
void update_point_SPHmass_bndry(SPHbody *btab, int SPHnobj, float newt, 
				bndry_t bndry);

/* In sphinit.c */
void *DarkRead(char *name, void *csdfp, void **btabp, int *gnobjp, int *nobjp, int set_id, int setpvel);
void *SPHRead(char *name, void *csdfp, SPHbody **btabp, int *gnobjp, int *nobjp, int set_id, int setpvel, float new_h, float new_u);
void SPHTestData(void *csdfp, SPHbody **btabp, int *gnobjp, int *nobjp, int periodic);
void *InitRead(char *name, void *csdfp, void **btabp, int *gnobjp, int *nobjp, 
	 SPHbody **SPHbtabp, int *SPHgnobjp, int *SPHnobjp, 
	 int set_id, int setpvel, float new_h, float new_u);
void DarkSPHTestData(void *csdfp, void **btabp, int *gnobjp, int *nobjp, 
		SPHbody **SPHbtabp, int *SPHgnobjp, int *SPHnobjp, int periodic);

/* In sphplus.c */
void GravPlusSPH(void **btab, int *nobj, SPHbody *SPHbtab, int SPHnobj);
void GravMinusSPH(void **btab, int *nobj, accbody **atab, int *anobj);

/* In sn.c */
void mmw(SPHbody *btab, int nobj);
void eosaux_setup(SPHbody *btab, int nobj);
void eos_prev(SPHbody *btab, int nobj);
void movebound(SPHbody *btab, int nobj, float t, float rb, float *vb, 
int *icore);
void pghost(SPHbody *btab, int nobj, int *nghost, Stk *ghosts, 
       float rb, float vb, float rbout, int iextf, int icore, float gg, 
       float xmcore, float aleph);
void remove_ghosts(SPHbody **btabp, int *nobjp);
void sn_gravity(SPHbody *btab, int nobj, float xmcore, float xmtheo, float gg, 
	   float clight, int icore, float rmin, float rmax);

/* In eos3.f */ 
void Fortran(eossetup)(void);
void Fortran(eos3)(float *steps, float *mass, double *rhoi, double *ui, 
                   double *udoti, 
		   double *u2i, double *yei, 
		   double *tempi, int *ifleosi, double *abari, double *xpi, 
		   double *xni, double *xpfi, double *p2i, double *p3i, 
		   double *p4i, double *temprev, double*rhoprev, 
		   double *xpprev, double *xnprev, double *yeprev, 
		   double *ufreez, double *xhe4n, double *xc12n,
		   double *xo16n, double *xne20n, double *xmg24n,
		   double *xsi28n, double *xni56n, double *xco56n, 
		   double *xfe56n, double *dtbrn); 
extern void *Fortran(output);
extern void *Fortran(konst);
extern void *Fortran(units);
extern void *Fortran(unit2);

/* common /konst/ gg, clight, arad, bigr, xsecnn, xsecne */
typedef struct {
  float gg;
  float clight;
  float arad;
  float bigr;
  float xsecnn;
  float xsecne;
} konst_s;

/* common /output/ vsoundi,pri,etai,yehi,xmuei,xmuhati,xalphai,
   $     xheavyi */
typedef struct {
  double vsound;
  double pr;
  double eta;
  double yeh;
  double xmue;
  double xmuhat;
  double xalpha;
  double xheavy;
} output_s;

/* common /units/ umass, udist, udens, utime, uergg, uergcc */
typedef struct {
  double umass;
  float udist;
  float udens;
  float utime;
  float uergg;
  float uergcc;
} units_s;

/* common /unit2/ utemp, utmev, ufoe, umevnuc, umeverg */
typedef struct {
  float utemp;
  float utemv;
  float ufoe;
  float umevnuc;
  float umeverg;
} unit2_s;

extern konst_s *konst;
extern output_s *output;
extern units_s *units;
extern unit2_s *unit2;
