#include "ndim.h"
#include "tree.h"
#include "key.h"
#include "complex.h"

typedef struct {
    float strength;		/* strength of body */
    float pos[NDIM];		/* position of body */
    float phi_r;		/* real part of scalar field */
    float phi_i;		/* imaginary part of scalar field */
    Key_t key;
    unsigned int ident;		/* identifier */
    float nterms;
} body, *bodyptr;

/* When we send a body from node to node, how much must we send??? */
#define TBODYSZ (1+NDIM)*sizeof(float)

typedef struct {
    float strength;		/* strength of body */
    float pos[NDIM];		/* position of body */
    float phi_r;		/* real part of scalar field */
    float phi_i;		/* imaginary part of scalar field */
    unsigned int ident;		/* identifier */
} outbody, *outbodyptr;

/* This is the descriptor that goes into the SDF header. */

#define OUTBODYDESC \
"struct {\n\
    float strength;		/* strength of body */\n\
    float x, y, z;		/* position of body */\n\
    float phi_r, phi_i;		/* complex scalar field */\n\
    unsigned int ident;		/* unique identifier */\n\
}"

typedef struct {
    float strength;
    float pos[NDIM];
    float sz;
    int nu, nv, ndeg;
    complex *ffsf;
    int daughters;
} cell, *cellptr;


/* This is the intermediate data structure used to construct cofm */
typedef struct {
    float strength;
    float pos[NDIM];
    float sz;
    int nu, nv, ndeg;
    complex *ffsf;
    int ndaughters;
} cofmdata;

typedef struct{
    Key_t key;
    float pos[NDIM];
    float sz;
    complex phi;
    int nu;
    int nv;
    int ndeg;
    complex *ffsf;
    int isbody;
    float daughters;
    float nterms;
    int interactions;
} Sink;

#define HAS_NTERMS
#define HAS_IDENT
#define HAS_KEY

/* Prototypes for all the functions which are "friends" of physics.h */
#include "physics_generic.h"

/* In shmz_ring.c */
void set_body(void *o, void *p);
void set_k(float lambda);
void do_shmz(void *p0, void *list, int bsize, int n);

/* In cofm.c */
void SetupCofm(int order, double lambda);
void CofmFromDaugh(hcellptr hptr, hcellptr daughters[]);
void CellFromCofm(cell *cp, cofmdata *cmp);

void SetTol(int gnobj);
void InheritSink(const Sink *from, Sink *to, hcell *pp);
void OutToIn(Sink *sink, const hcell **source, int *result, int n);


/* In print.c */
char *PrintCellContents(const cell *p);
char *PrintBodyContents(const body *p);
