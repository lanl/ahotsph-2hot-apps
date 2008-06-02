/*
 * Copyright 1996 Michael S. Warren and John K. Salmon.  All Rights Reserved.
 */

#include "tree.h"
#include "key.h"
#include "timers.h"
#include "ndim.h"

#define SPH_SAVE_ACC
#define POS_IS_DOUBLE
#define SPH_GRAV

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
    float prnu;
    int ifleos;
    float dt_next;
    float gshift;
    float r;
    float dnue, dnueb, dnux;
    float ynue, ynueb, ynux;
    float enuet, enuebt, enuxt;
    float tempnue, tempnueb, tempnux;
    float etanue, etanueb, etanux;
    /* These are needed for making cells for SPHreduce */
    float unue, unueb, unux;
    float u2;
    float ufreez;
    /* Things declared above this line are communicated between processors */
    /* so they can be used in in the loop over nbrs in FindRho and ForceSPH */
    /* Don't add anything above this line unless you fix TBODYSZ */
    float acc[NDIM];
    float grav_acc[NDIM];
    float acc_last[NDIM];
    float phi;
#ifdef POS_IS_DOUBLE
    double pos_last[NDIM];
#else
    float pos_last[NDIM];
#endif
    Key_t key;
    unsigned int ident;
    float nterms;
    float grav_nterms;
    float lvel[NDIM];
    float drho_dt;
    float hdot;
    float udot;
    float udot_last;
    float udot2;
    float udot2_last;
    unsigned int nbrs;
    float tacc;
    float dt;
    float min_nbr_dt;
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
    float dynue, dynueb, dynux;
    float dunue, dunueb, dunux;
    float eta;
    float xmuhat;
    float xmue;
    float dunu;
    float prg;
    float taccreted;
    int iteraccreted;
    int bghost;
    void *ireal;
    short ebeta, pbeta;
    float dq;                   /* viscous heating */
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
    float mass;		/* mass of body */
    float pos[NDIM];		/* position of body */
    float vel[NDIM];		/* velocity of body */
    float u;
    float h;
    float rho;
    float drho_dt;
    float udot;
    float pr;
    float prnu;
    float vsound;
    float temp;
    float ye;
    float xp;
    float xn;
    float u2;
    float abar;
    float ynue;
    float ynueb;
    float ynux;
    float unue;
    float unueb;
    float unux;
    float ufreez;
    float dnue;
    float dnueb;
    float dnux;
    float enuet;
    float enuebt;
    float enuxt;
    float dye;
    float dunu;
    float dynue;
    float dunue;
    float eta;
    float tempnue;
    float etanue;
    float xpf;
    float p2;
    float p3;
    float p4;
    int  ifleos;
#ifdef SPH_SAVE_ACC
    float acc[NDIM];
    float grav_acc[NDIM];
 /* float acc_last[NDIM]; */
    float phi;
    float dt;
#endif
    float taccreted;
    int iteraccreted;
    unsigned int nbrs;
    unsigned int ident;		/* unique? identifier */
    float dq;
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
    float prnu;			\n\
    float vsound;		\n\
    float temp;			\n\
    float ye;			\n\
    float xp;			\n\
    float xn;			\n\
    float u2;			\n\
    float abar;			\n\
    float ynue;			\n\
    float ynueb;		\n\
    float ynux;			\n\
    float unue;			\n\
    float unueb;		\n\
    float unux;			\n\
    float ufreez;		\n\
    float dnue;			\n\
    float dnueb;		\n\
    float dnux;			\n\
    float enuet;		\n\
    float enuebt;		\n\
    float enuxt;		\n\
    float dye;			\n\
    float dunu;			\n\
    float dynue;		\n\
    float dunue;		\n\
    float eta;			\n\
    float tempnue;		\n\
    float etanue;		\n\
    float xpf;			\n\
    float p2;			\n\
    float p3;			\n\
    float p4;			\n\
    int  ifleos;		\n\
    float ax, ay, az;		/* acceleration */\n\
    float gax, gay, gaz;	/* acceleration due to gravity */\n\
    float phi;			/* potential */\n\
    float idt;			/* timestep */\n\
    float taccreted;            /* time of accretion */\n\
    int iteraccreted;           /* iteration of accretion */\n\
    unsigned int nbrs;		/* number of neighbors */\n\
    unsigned int ident;		/* unique identifier */\n\
    float dq;                   /* viscous heating */\n\
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
    float prnu;			\n\
    float vsound;		\n\
    float temp;			\n\
    float ye;			\n\
    float xp;			\n\
    float xn;			\n\
    float u2;			\n\
    float abar;			\n\
    float ynue;			\n\
    float ynueb;		\n\
    float ynux;			\n\
    float unue;			\n\
    float unueb;		\n\
    float unux;			\n\
    float ufreez;		\n\
    float dnue;			\n\
    float dnueb;		\n\
    float dnux;			\n\
    float enuet;		\n\
    float enuebt;		\n\
    float enuxt;		\n\
    float dye;			\n\
    float dunu;			\n\
    float dynue;		\n\
    float dunue;		\n\
    float eta;			\n\
    float tempnue;		\n\
    float etanue;		\n\
    float xpf;			\n\
    float p2;			\n\
    float p3;			\n\
    float p4;			\n\
    int  ifleos;		\n\
    float taccreted;            /* time of accretion */\n\
    int iteraccreted;           /* iteration of accretion */\n\
    unsigned int nbrs;		/* number of neighbors */\n\
    unsigned int ident;		/* unique identifier */\n\
    float dq;                   /* viscous heating */\n\
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
    float prnu;			\n\
    float vsound;		\n\
    float temp;			\n\
    float ye;			\n\
    float xp;			\n\
    float xn;			\n\
    float u2;			\n\
    float abar;			\n\
    float ynue;			\n\
    float ynueb;		\n\
    float ynux;			\n\
    float unue;			\n\
    float unueb;		\n\
    float unux;			\n\
    float ufreez;		\n\
    float dnue;			\n\
    float dnueb;		\n\
    float dnux;			\n\
    float enuet;		\n\
    float enuebt;		\n\
    float enuxt;		\n\
    float dye;			\n\
    float dunu;			\n\
    float dynue;		\n\
    float dunue;		\n\
    float eta;			\n\
    float tempnue;		\n\
    float etanue;		\n\
    float xpf;			\n\
    float p2;			\n\
    float p3;			\n\
    float p4;			\n\
    int  ifleos;		\n\
    float ax, ay;		/* acceleration */\n\
    float gax, gay;		/* acceleration due to gravity */\n\
    float phi;			/* potential */\n\
    float idt;			/* timestep */\n\
    float taccreted;            /* time of accretion */\n\
    int iteraccreted;           /* iteration of accretion */\n\
    unsigned int nbrs;		/* number of neighbors */\n\
    unsigned int ident;		/* unique identifier */\n\
    float dq;                   /* viscous heating */\n\
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
    float prnu;			\n\
    float vsound;		\n\
    float temp;			\n\
    float ye;			\n\
    float xp;			\n\
    float xn;			\n\
    float u2;			\n\
    float abar;			\n\
    float ynue;			\n\
    float ynueb;		\n\
    float ynux;			\n\
    float unue;			\n\
    float unueb;		\n\
    float unux;			\n\
    float ufreez;		\n\
    float dnue;			\n\
    float dnueb;		\n\
    float dnux;			\n\
    float enuet;		\n\
    float enuebt;		\n\
    float enuxt;		\n\
    float dye;			\n\
    float dunu;			\n\
    float dynue;		\n\
    float dunue;		\n\
    float eta;			\n\
    float tempnue;		\n\
    float etanue;		\n\
    float xpf;			\n\
    float p2;			\n\
    float p3;			\n\
    float p4;			\n\
    int  ifleos;		\n\
    float taccreted;            /* time of accretion */\n\
    int iteraccreted;           /* iteration of accretion */\n\
    unsigned int nbrs;		/* number of neighbors */\n\
    unsigned int ident;		/* unique identifier */\n\
    float dq;                   /* viscous heating */\n\
}"
#endif
#endif /* NDIM==2 */

typedef struct {
    float mass;
    float pos[NDIM];
    float vel[NDIM];
    float bmax;
    float lap;
    int daughters;
    int ident;
    float u;
    float ifleos;
    float abar;
    float temp;
    float ye;
    float xp;
    float xn;
    float u2;
    float ynue;
    float ynueb;
    float ynux;
    float unue;
    float unueb;
    float unux;
    float ufreez;
} SPHcell;


/* This is the intermediate data structure used to construct cofm */
typedef struct{
    float mass;
    float pos[NDIM];
    float vel[NDIM];
    float massinv;
    float bmax;
    float B2;
    float sz;
    float lap;
    int ndaughters;
    int ident;
    float u;
    float ifleos;
    float abar;
    float temp;
    float ye;
    float xp;
    float xn;
    float u2;
    float ynue;
    float ynueb;
    float ynux;
    float unue;
    float unueb;
    float unux;
    float ufreez;
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
    float prnu;
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
    float dynue, dynueb, dynux;
    float dunue, dunueb, dunux;
    float ynue, ynueb, ynux;
    float unue, unueb, unux;
    float enuet, enuebt, enuxt;
    float dnue, dnueb, dnux;
    float tempnue, tempnueb, tempnux;
    float etanue, etanueb, etanux;
    float gshift;
    float xfac;			/* geometrical factor */
    int ident;			/* not necessary, but useful for debugging */
} SinkSPH;

/* bndry_t needs double-precision members to store the accumulating
   sums of small quantities from accreted particles */
typedef struct {
    double pos[NDIM];
    double vel[NDIM];
    double j[NDIM];
    double mass;
    double a;
    double r;
    double force_r;
} bndry_t;

/* External Fortran linkage */
#define Fortran(x) x##_
#define Fortran2(x) x##_

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
void update_final(SPHbody *btab, int nobj, float dt, int *limit_high, int *limit_low, float dttol);
void update_intermediate(SPHbody *btab, int nobj, float dt_last, int flag, int *limit, float xmtheo);
void SPH_setup(int dim, int ncoef1, double *coef1, int ncoef2, double *coef2,
	       float maxnue, float maxnueb, float maxnux,
	       float enue, float enueb, float enux,
	       float e2nue, float e2nueb, float e2nux,
	       float ftrape, float ftrapb, float ftrapx);
void Getrmax(float *maxnue, float *maxnueb, float *maxnux, 
	     float *enue, float *enueb, float *enux, 
	     float *e2nue, float *e2nueb, float *e2nux);

void SetSPHOffset(float *off, float *voff);
void UnSetSPHOffset(void);
void SetSPHRotate(float angle);
void UnSetSPHRotate(void);
void update_point_SPHmass(SPHbody *btab, int nobj, void *p, float smooth2, float newt);
void update_bardeen(SPHbody *btab, int nobj, float G, float c, bndry_t b);

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
void SPHreduce(tree_t *sphtree, float bmax, float rmax, 
	       char *outnamebase, int iter, 
	       float gnewt, float dt, float tpos, float tvel,
	       float rmaxnue, float rmaxnueb, float rmaxnux, 
	       float enue, float enueb, float enux, 
	       float e2nue, float e2nueb, float e2nux,
	       float ftrape, float ftrapb, float ftrapx);

/* In sn.c */
void mmw(SPHbody *btab, int nobj);
void eosaux_setup(SPHbody *btab, int nobj);
void eos_prev(SPHbody *btab, int nobj);
void eosgen_setup(SPHbody *btab, int nobj);
void movebound(SPHbody *btab, int nobj, float t, double rb, double *vb, 
int *icore);
void pghost(SPHbody **btab, int *nobj, int *gnobj, 
	    double rb, double vb, double rbout, int iextf, int icore, 
	    float gg, float xmcore, float aleph, int do_ghosts);
void remove_ghosts(SPHbody **btabp, int *nobjp, int *gnobjp);
void sn_gravity(SPHbody *btab, int nobj, float xmcore, float xmtheo, float gg, 
	   float clight, int icore, float rmin, float rmax);

/* In eos3.f */ 
void Fortran(eossetup)(void);
void Fortran(eos3)(double *rhoi, double *ui, double *u2i, double *yei, 
		   double *tempi, int *ifleosi, double *abari, double *xpi, 
		   double *xni, double *xpfi, double *p2i, double *p3i, 
		   double *p4i, double *temprev, double*rhoprev, 
		   double *xpprev, double *xnprev, double *yeprev, 
		   double *ufreez, int *iident, int *iprocnum);
void Fortran(neutrino)(float *steps, float *rhok, float *yek, float *xpk,
		       float *xnk, float *hi, 
		       double *xheavyk, double *xalphak, double *yehk, 
		       float *etak,
		       float *tempk, float *abark, float *gshifti, float *ri,
		       float *pmassi, float *vsoundk, float *xmuhatk,
		       float *ynuei, float *ynuebi, float *ynuxi, float *unuei,
		       float *unuebi, float *unuxi, float *rmaxnue,
		       float *rmaxnueb, float *rmaxnux, float *ftrape,
		       float *ftrapb, float *ftrapx, float *jtrape,
		       float *jtrapb, float *jtrapx, int *ident);
void Fortran(neutrino2)(float *steps, float *rhok, float *xpk, float *xnk,
			float *etak, float *tempk, float *ri, float *pmassi,
			float *vsoundi, float *xmuhati, float *ynuei,
			float *ynuebi, float *ynuxi, 
			float *unuei, float *unuebi, float *rlumnue,
			float *rlumnueb, float *rlumnux, float *enue,
			float *enueb, float *enux, float *e2nue, float *e2nueb,
			float *e2nux, float *gshifti, float *rmaxnue, 
			float *rmaxnueb, float *rmaxnux, float *hi);
extern void *Fortran(output);
extern void *Fortran(neutout);
extern void *Fortran(nuout);
extern void *Fortran(nulums);
extern void *Fortran(konst);
extern void *Fortran(units);
extern void *Fortran(unit2);
extern void *Fortran(beta);
extern void *Fortran(nutrap);

/* In eosgen.f */ 
void Fortran(eosgen)(double *rhoi, double *tempi, double *yei, 
		     double *abari, double *ui, double *u2i, 
		     double *pri, double *xpi, double *xni, 
		     double *ufreez, int *ifleosi, int *iident, int *iprocnum);

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

/* common /neutout/ dyei,dynuei,dynuebi,dynuxi,
   $     tempnuei,tempnuebi,tempnuxi,enueti,enuebti,enuxti,
   $     dnuei,dnuebi,dnuxi,dunuei,dunuebi,dunuxi,dunui,
   $     etanuei,etanuebi,etanuxi,prnui */
typedef struct {
  float dye;
  float dynue;
  float dynueb;
  float dynux;
  float tempnue;
  float tempnueb;
  float tempnux;

  float enuet;
  float enuebt;
  float enuxt;
  float dnue;
  float dnueb;
  float dnux;
  float dunue;
  float dunueb;
  float dunux;
  float dunu;
  float etanue;
  float etanueb;
  float etanux;
  float prnu;
} neut_out_s;

/* common /nuout/ rmxnue,rmxnueb,rmxnux */
typedef struct {
  float rmxnue;
  float rmxnueb;
  float rmxnux;
} nu_out_s;

/* common /nulums/ rlumnue, rlumnueb, rlumnux,
   $     enue, enueb, enux, e2nue, e2nueb, e2nux */
typedef struct {
  float rlumnue;
  float rlumnueb;
  float rlumnux;
  float dlumnu;
  float enue;
  float enueb;
  float enux;
  float e2nue;
  float e2nueb;
  float e2nux;
  float enues;
  float enuebs;
  float enuxs;
  float dee;
  float deeb;
  float dex;
} nu_lums_s;

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

typedef struct {
  int ebetaeq;
  int pbetaeq;
} beta_s;

typedef struct {
  float dtrapnue;
  float dtrapnueb;
} nutrap_s;

extern konst_s *Konst;
extern output_s *Outputf;
extern neut_out_s *Neut_out;
extern nu_out_s *Nu_out;
extern nu_lums_s *Nu_lums;
extern units_s *Units;
extern unit2_s *Unit2;
extern beta_s *Nubeta;
extern nutrap_s *Nutrap;
extern float ftrape;
extern float ftrapb;
extern float ftrapx;

