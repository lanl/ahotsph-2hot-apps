#ifndef PHYSICS_vrtxDOTh
#define PHYSICS_vrtxDOTh

/* Version for the gaussian smoothing 
   Also, error estimates are on the norm of vorticity vector (march 93)
*/
/* Jan 1995,
     added optional (#ifdef) double precision for certain accumulators. */

#define NDIM 3
#include "tree.h"
#include "key.h"

extern float epsilon;
extern float kc;		/* kernel cutoff */
extern float kc2;		/* kernel cutoff squared */
extern float epsinv;
extern float errtol;
extern float nu;

/* The transmitted bodies must cary pos[3], str[3] and vol. */
#define TBODYSZ (7*sizeof(float))


#ifdef DOUBLE_ACCUMULATORS
/* Use ACCUM and PAD_DECL in the struct declaration, and use
   ACCUM_S and PAD_DECL_S in the BODY_DESC string.  Ugly! */
#define ACCUM double
#define ACCUM_S "double"
#define PAD_DECL float _junk_pad1;
#define PAD_DECL_S "float _junk_pad1;"
#else
#define ACCUM float
#define ACCUM_S "float"
#define PAD_DECL 
#define PAD_DECL_S "/* No padding */"
#endif

/* structures.h contains the particular structure definitions appropriate for
   mono, di or quadrupole interactions.  It is in the private directories. */
#include "structures.h"

/* Tell physics.c that we have nterms and ident */
#define HAS_NTERMS
#define HAS_IDENT
#define HAS_KEY

#define Errsum(x)      ((x)->errsum)
#define Errsum2(x)     ((x)->errsum2)
#define Pos(x)         ((x)->pos)
#define Strength(x)    ((x)->strength)
#define Vol(x)         ((x)->vol)
#define Psi(x)         ((x)->psi)
#define Vel(x)         ((x)->vel)
#define Gradvel(x)     ((x)->gradvel)
#define Dstr(x)        ((x)->dstr)
#define Vel_old(x)     ((x)->vel_old)
#define Dstr_old(x)    ((x)->dstr_old)
#define Dpole(x)       ((x)->dpole)
#define Qpole(x)       ((x)->qpole)
#define Bmax(x)        ((x)->bmax)
#define Daughters(x)   ((x)->daughters)
#define B0(x) ((x)->b0)
#define B1(x) ((x)->b1)
#define B2(x) ((x)->b2)
#define B3(x) ((x)->b3)
#define B4(x) ((x)->b4)
#define Omegat(x) ((x)->omegat)
    
/* Prototypes for all the functions which are "friends" of physics.h */
#include "physics_generic.h"

#ifdef __cplusplus
extern "C"{
#endif /* __cplusplus */
/* From cofm_vrtx.c */
void CofmFromDaugh(hcell *, hcell **);
void CellFromCofm(cell *, cofm_data *);

/* From mac.c */
void NlgNInherit(const Sink *from, Sink *to, hcell *pp);
int NlgNMAC(Sink *sink, const hcell *source);
void NlgNMACv(Sink *sink, const hcell **sources, int *results, int nsrc);
extern Timer_t VrtxTm;

/* From integrate_vrtx_rk.c */
void Update(bodyptr btab, int n, float dta, float dtb, float dtrel, int iflag);

/* from quads_vrtx.c */
void InteractCell(cellptr cp, bodyptr me);
extern Counter_t CellCnt;

/* from mono_vrtx.c */
void InteractBody(body *bp, bodyptr me, 
		    float eps2inv12, float nu);
extern Counter_t BodyCnt;
extern Counter_t FullKernelCnt;
extern Counter_t TaylorKernelCnt;

/* from diags.c */
void GlobalDiags(body *btab, int nobj);
extern Counter_t NtermsCnt;

/* from relaxomega.c */
void RelaxOmega(body *btab, int nobj, float relaxw);

/* from omegat.c */
void OmegatBody(body *bp, body *me, float eps2inv12);

/* from remesh.c */
void remesh(body **btabp, int *nobjp, float remesh_h, float remesh_min_str);

/* from fixomegatot.c */
void FixOmegaTot(bodyptr btab, int nobj, int gnobj);

/* from print.c */
char *PrintBodyContents(const body *p);
char *PrintCellContents(const cell *p);
char *PrintBranch(const cofm_data *cmp);

#ifdef __cplusplus
}
#endif

#endif
