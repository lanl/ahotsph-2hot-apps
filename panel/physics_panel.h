/* Compute u = curl(psi) + grad(phi) */

/* Panel method for only a scalar charge (sigma) */

/* Error estimates analogous to gravitational problem (grad of potential) */

/* Multipole expansions centered at the centroid of charge abs value  */

#define NDIM 3
#include "tree.h"



typedef struct {
 /* this part is used for exact integration formula: */
    float pos1[3];            /* position of corner 1 */
    float pos2[3];            /* position of corner 2 */
    float pos3[3];            /* position of corner 3 */
    float pos[3];             /* position of centroid */
    float ex[3];              /* ex= unit vector parallel to 1-3, */
    float ey[3];              /*   ez= outward normal when traveling as 1-2-3, */
    float ez[3];              /*   ey= ez X ex       */
    float x1l, y1l;
    float x2l, y2l;
    float x3l;                /* y3l=y1l, hence s31=0 */
    float d12, c12, s12;
    float d23, c23, s23;
    float d31, c31;           /* s31=0 */
    float sphi, dxlsphi, dylsphi; /* self-potential phi at panel centroid (xl=yl=0,
                                     zl=0_) and its local derivatives, for a unit 
                                     panel scalar strength. dzlsphi=2*pi... */ 

 /* this part is used for panel multipole expansion formula: */

    float ip;                /* panel area = panel monopole term */ 
               
    float ixxlp;             /* panel multipole expansion coefficients with respect */
    float iyylp;              /*   to panel centroid and in local coodinate system */
    float ixylp;             /*    (panel of unit strength). */

    float ixxp;              /* panel multipole expansion coefficients with respect */
    float iyyp;              /*   to panel centroid and in absolute coodinate system */
    float izzp;              /*   (panel of unit strength).  */
    float ixyp;
    float ixzp;
    float iyzp;

    float size, dist2crit; 
                             /* size=max(d1c, d2c, d3c) where c=centroid
                                dist2crit = (dist2)_crit for switching from
                                full formula to multipole expansion formula
                                in mono_panel.c                              */

    float sigma;		/* panel scalar charge. */

/* The transmitted data stops here. */
/* BELOW THIS LINE is information used only by the iteration.  It is not necessary */
/* transmit this information between processors for the multipole method. */

    float phi;			/* scalar potential */
    float vel[3];             /* velocity at panel centroid = curl(psi) + grad(phi) */
    float dsigma;		/* change in scalar charge */
    float uext[3];		/* the rhs of the equations */
    float uexact[3];		/* the 'exact' analytic solution, if known */

    float errsum, errsum2;	/* sums of error bounds */
    unsigned int ident;       /* panel identification number */
    int nterms;	              /* number of terms in field eval */
} body, *bodyptr;

/* How much to send when we communicate a body */
#define TBODYSZ offsetof(body, phi)

typedef struct {
  float pos[3];
  float sigma;
  float vel[3];
/*
  float uext[3];
  float gamma[3];
  float vnorm, vtan;
  float errsum, errsum2;
*/
  int ident;
} outbody;
  
#define OUTBODYDESC \
"struct {\n\
   float x, y, z;\n\
   float sigma;\n\
   float velx, vely, velz;\n\
/*\n\
   float uextx, uexty, uextz;\n\
   float gammax, gammay, gammaz;\n\
   float vnorm, vtan;\n\
   float errsum, errsum2;\n\
*/\n\
   int ident;\n\
}"

#define WHOLEBODYDESC \
"struct {\n\
    float x1,y1,z1;            \n\
    float x2,y2,z2;            \n\
    float x3,y3,z3;            \n\
    float x,y,z;               \n\
    float exx, exy, exz;       \n\
    float eyx, eyy, eyz;       \n\
    float ezx, ezy, ezz;       \n\
    float x1l, y1l;            \n\
    float x2l, y2l;            \n\
    float x3l;                 \n\
    float d12, c12, s12;       \n\
    float d23, c23, s23;       \n\
    float d31, c31;            \n\
    float sphi, dxlsphi, dylsphi; \n\
    float ip;     \n\
    float ixxlp;   \n\
    float iyylp;   \n\
    float ixylp;   \n\
    float ixxp;   \n\
    float iyyp;   \n\
    float izzp;   \n\
    float ixyp;   \n\
    float ixzp;   \n\
    float iyzp;   \n\
    float size, dist2crit; \n\
    float sigma;	\n\
    float phi;		\n\
    float vx, vy, vz;          \n\
    float dsigma;		\n\
    float uextx, uexty, uextz;  \n\
    float uexactx, uexacty, uexactz; \n\
    float errsum, errsum2;	\n\
    unsigned int ident;        \n\
    int nterms;		       \n\
}"

struct dmoment{      /* dipole terms of multipole expansion for a group of panels */
    /* The scalar components */
    float x, y, z;
};

struct qmoment{    /* quadrupole terms of multipole expansion for a group of panels */
    float xx;             
    float yy;
    float zz;
    float xy;
    float xz;
    float yz;
};

typedef struct {
    float pos[3];             /* expansion center (= geometric cell center)*/
    float sigma;		/* scalar charge */
    struct dmoment dpole;     /* dipole terms */
    struct qmoment qpole;     /* quadrapole terms */
    float bmax;	              /* bmax of a cell of cell */
    float b3;                 /* estimate of b3, with b3 <= bmax*b2       */ 
    float b4;                 /* estimate of b4, with b4 >= b2^2/b0       */
    float rcrit2;	      /* minimum squared interaction distance */
    int daughters;            /* how many particles in that cell */
    unsigned int sub_flags;
} cell, *cellptr;

typedef struct {
    body *bp;
} Sink;

typedef  struct {
    float pos[3];
    float sigma;
    struct dmoment dpole;
    struct qmoment qpole;
    float b0;
    float b2;
    float bmax;
    float maxsize;
    int daughters;
} cofm_data;


   
/* Tell physics.c that we have nterms and ident */
#define HAS_NTERMS
#define HAS_IDENT

#define Pos(x) 		((x)->pos)
#define Pos1(x)         ((x)->pos1)
#define Pos2(x)         ((x)->pos2)
#define Pos3(x)         ((x)->pos3)
#define Sigma(x)	((x)->sigma)
#define Phi(x)		((x)->phi)
#define Vel(x)         ((x)->vel)
#define Dpole(x)       ((x)->dpole)
#define Qpole(x)       ((x)->qpole)
#define Bmax(x)        ((x)->bmax)
#define B3(x)          ((x)->b3)
#define B4(x)          ((x)->b4)
#define Daughters(x)   ((x)->daughters)
#define Ip(x)          ((x)->ip)
#define Ixxp(x)        ((x)->ixxp)
#define Iyyp(x)        ((x)->iyyp)
#define Izzp(x)        ((x)->izzp)
#define Ixyp(x)        ((x)->ixyp)
#define Ixzp(x)        ((x)->ixzp)
#define Iyzp(x)        ((x)->iyzp)
#define Ixxlp(x)       ((x)->ixxlp)
#define Iyylp(x)       ((x)->iyylp)
#define Ixylp(x)       ((x)->ixylp)
#define Size(x)        ((x)->size)
#define Ex(x)          ((x)->ex)
#define Ey(x)          ((x)->ey)
#define Ez(x)          ((x)->ez)
#define Uext(x)		((x)->uext)
#define Errsum(x)      ((x)->errsum)
#define Errsum2(x)     ((x)->errsum2)

#define B0(x)          ((x)->b0)
#define B2(x)          ((x)->b2)

/* Prototypes for all the functions which are "friends" of physics.h */
#include "physics_generic.h"


/* From cofm_panel.c */
extern float errtol; 
void CofmFromDaugh(hcell *, hcell **);
void CellFromCofm(cell *, cofm_data *);

/* From mac.c */
void NlgNInherit(const Sink *from, Sink *to, hcell *pp);
void NlgNMACv(Sink *sink, const hcell **source_vec, int *result, int n);

/* From prepare_panel.c */
void PreparePanel(bodyptr btab, int n,
		  void (*Uexternal)(body *bp));

/* From integrate_panel.c */
void Update(bodyptr btab, int n, float relax, float *residual);

/* From quads_panel.c */
extern Counter_t CellInt;
void Cinter(body *me, cell *src);

/* From mono_panel.c */
extern Counter_t BodyFullCnt;
extern Counter_t BodyQuadCnt;
void Binter(body *me, body *src);
