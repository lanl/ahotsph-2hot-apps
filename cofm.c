#include <math.h>
#include <stdio.h>
#include "tree.h"
#include "order.h"
#include "physics.h"
#include "vop.h"
#include "chn.h"
#include "Malloc.h"
#include "Msgs.h"
#include "error.h"
#include "Assert.h"
#include "fastflpt.h"
#include "protos.h"

static int MACtype = BMAX_MAC;		/* default */
static float Tol;
static float invTol;
static float RelTol;
static float invRelTol;
static float RelTol0;
static float invRelTol0;
static float Bmax0;
static float Ptol_boost;
static float Stol_max;
static int Quad_Ncut = 7;
static int Hexa_Ncut = 20;
static tree_t *Tree;

static float min_sigma_m = 1.0f;
static float max_sigma_m = 7.2e+10; /* exp(25.0) */
static float sigma_m[] = {
    1.0,			/* 0, 10^10 Msol */
    1.0,
    1.0,
    1.0,
    1.0,
    1.0,
    1.0,
    1.0,
    1.0,			/* 8, 3 x 10^13 Msol, r = 7.8/h Mpc */
    1.5,
    2.0,
    3.0,
    6.0,
    10.88,
    15.18,
    21.34,
    31.58,
    46.73,
    74.63,
    120.3,
    199.9,
    353.1,
    611.2,
    1117,			/* 23, r = 1000/h Mpc */
    1983};			/* 24, 2.65 x 10^20 Msol */

static void mpole_add_mono(cofmdata *cmp, float m, float x, float y, float z, float Rinv);
static void mpole_add_shift(cofmdata *cmp, float m, float x, float y, float z, float Rinv, 
			    cofmdata *dp);

void SetupCofm(int type, float tol, float rel_tol, float rel_tol0, float r0, float ptol_boost, 
	       float stol_max, int qcut, int hcut,  tree_t *t)
{
    MACtype = type;
    Tol = tol;
    invTol = 1.0/tol;
    RelTol = rel_tol;
    invRelTol = 1.0/rel_tol;
    RelTol0 = rel_tol0;
    invRelTol0 = 1.0/rel_tol0;
    Bmax0 = r0;
    Ptol_boost = ptol_boost;
    Stol_max = stol_max;
    Quad_Ncut = qcut;
    Hexa_Ncut = hcut;
    Tree = t;
}

void CofmFromDaugh(hcellptr hptr, hcellptr daughters[]){
    int i;
    cofmdata *dp;
    cofmdata *cmp;
    body *bp = NULL;
    float dmass;
    float newbmax;
    float center[NDIM], cellsz;
    Vxd(float dx);
#if defined(QUAD) || defined(HEXA)
    float Rinv;
#endif

    assert(Sub_Flags(hptr));

    cmp = hptr->ptr;
    assert(cmp);
    memset(cmp, 0, sizeof(cofmdata));
    CELLCORNER(hptr->key, center, &cellsz);
#if defined(QUAD) || defined(HEXA)
    Rinv = recipf(0.5f*cellsz);
#endif

    /* First get the cm of the new cell. */
    for(i=0; i<(1<<NDIM); i++){
	if (daughters[i] == NULL)
	  continue;
	if (Sub_Flags(daughters[i]) == 0) {
	    bp = daughters[i]->ptr;
	    dmass = bp->mass;
	    cmp->mass += dmass;
	    VV(cmp->pos, += dmass * bp->pos);
	    cmp->ndaughters++;
	} else {
	    dp = daughters[i]->ptr;
	    dmass = dp->mass;
	    cmp->mass += dmass;
	    VV(cmp->pos, += dmass * dp->pos);
	    cmp->ndaughters += dp->ndaughters;
	} 
    }
    /* Divide out the total mass */
    if (cmp->mass != (float)0.) {
	cmp->massinv = recipf(cmp->mass);
	VS(cmp->pos, *= cmp->massinv);
    } else {
	Error("Zero mass in BranchFromDaughters!\n");
    }
    /* Now loop again to pick up B2, etc.  */
    for (i=0; i<(1<<NDIM); i++) {
	float tmp[NDIM];
	float tmpsq;
	
	if(daughters[i] == NULL)
	    continue;
	if (Sub_Flags(daughters[i]) == 0) {
	    bp = daughters[i]->ptr;
	    VVV(tmp, = bp->pos, - cmp->pos);
	    newbmax = (float)0.;
	    tmpsq = Dot(tmp, tmp);
	    if (tmpsq != 0.F) {
		cmp->B2 += bp->mass * tmpsq;
		newbmax += sqrtf_fast(tmpsq);
#if defined(QUAD) || defined(HEXA)
		mpole_add_mono(cmp, bp->mass, tmp[0], tmp[1], tmp[2], Rinv);
#endif
	    }
	} else {
	    dp = daughters[i]->ptr;
	    VVV(tmp, = dp->pos, -cmp->pos);
	    newbmax = dp->bmax;
	    tmpsq = Dot(tmp, tmp);
	    cmp->B2 += dp->B2;
#if defined(QUAD) || defined(HEXA)
	    mpole_add_shift(cmp, dp->mass, tmp[0], tmp[1], tmp[2], Rinv, dp);
#endif
	    if (tmpsq != 0.F) {
		cmp->B2 += dp->mass * tmpsq;
		newbmax += sqrtf_fast(tmpsq);
	    }
	}
	if (newbmax > cmp->bmax)
	  cmp->bmax = newbmax;
    }
    if (cmp->B2 == 0.0) {
      for(i=0; i<(1<<NDIM); i++){
	if (daughters[i] == NULL)
	  continue;
	if (Sub_Flags(daughters[i]) == 0) {
	  bp = daughters[i]->ptr;
	  Msg_do("B2 body %s\n", PrintBodyContents(bp));
	} else {
	  dp = daughters[i]->ptr;
	  Msg_do("B2 cell\n");
	} 
      }
      cmp->B2 = fabs(bp->pos[0] * 1e-7); /* HACK */
      Warning("cmp->B2 is zero, suspect identical particle positions\n");
    }

    /* This is an alternative bound on bmax, which is sometimes tighter */
    /* than the cumulative bound computed above. */
    cmp->sz = cellsz;		/* for pure Barnes-But */
    cellsz *= (float)0.5;
    VS(center, += cellsz);
    VxVVS(dx, = cellsz+ fabs LPAREN cmp->pos, - center,  RPAREN);
    newbmax = sqrtf_fast(Dotx(dx, dx));
    cmp->bmax = (newbmax < cmp->bmax) ? newbmax : cmp->bmax;
    hptr->ptr = cmp;
}

static double a[9];		/* coef of error poly */
static void rcrit_poly(int n, double r, double *value, double *deriv);
static double rtnewt(int n, void (*funcd)(int, double, double *, double *),  
		    double x1, double xacc);

/* Tell the tree internals how big an object to copy */
int 
CellSz(void *p)
{
    const cell *cp = p;

    if (Hexa_Ncut && (cp->daughters >= Hexa_Ncut)) {
	return sizeof(hexacell);
    } else if (Quad_Ncut && (cp->daughters >= Quad_Ncut)) {
	return sizeof(quadcell);
    } else {
	return sizeof(cell);
    }
}

/* Turn the ptr from a cofmdata to a cell. */
void *CellFromCofm(cofmdata *cmp)
{
    hexacell *cp;

    if (Hexa_Ncut && (cmp->ndaughters >= Hexa_Ncut)) {
	cp = ChnAlloc(&Tree->hexacellchn);
    } else if (Quad_Ncut && (cmp->ndaughters >= Quad_Ncut)) {
	cp = ChnAlloc(&Tree->quadcellchn);
    } else {
	cp = ChnAlloc(&Tree->cellchn);
    }
    cp->mass = cmp->mass;
    VV(cp->pos, = cmp->pos);
    cp->bmax = cmp->bmax;
    cp->R = 0.5f*cmp->sz;
    if (MACtype == AREL_MAC) {
	float abs_rcrit;
	float rel_rcrit, rel_rcrit0;
	float B3, bmaxhalf, rcritmax;
	float B2 = cmp->B2;
#if defined(QUAD) || defined(HEXA)
	double R2, R4, B4, B5;
#endif
	float ptol;

	if (Stol_max > 1.0f) {
	    if (cmp->mass <= min_sigma_m) {
		ptol = Tol;
	    } else if (cmp->mass > max_sigma_m) {
		ptol = Stol_max*Tol;
	    } else {
		int idx = logf(cmp->mass);
		float s = sigma_m[idx];
		if (s < Stol_max) {
		    ptol = s*Tol;
		} else {
		    ptol = Stol_max*Tol;
		}
	    }
	} else {
	    ptol = Tol/(1.+Ptol_boost*cmp->bmax/Bmax0);
	}
	bmaxhalf = cmp->bmax * (float)0.5;
	rcritmax = bmaxhalf + sqrtf_fast(bmaxhalf*bmaxhalf
					 + sqrtf_fast((float)3.*B2/ptol));
	if (!finite(rcritmax))
	  Error("Bad rcritmax, q->bmax = %g, B2 = %g\n", cmp->bmax, B2);
	if (B2 == (float)0.0) Error("B2 is zero\n");
	B3 = B2 * sqrtf_fast(B2*cmp->massinv); /* lower bound */
	if (!finite(B2) || !finite(B3) || !finite(cmp->bmax))
	  Error("Bad value B2 = %g, B3 = %g, bmax = %g\n", B2, B3, cmp->bmax);
	a[0] = 2.*B3;
	a[1] = -3. * B2;
	a[2] = 0.;
	a[3] = ptol*cmp->bmax*cmp->bmax;
	a[4] = -2. * ptol * cmp->bmax;
	a[5] = ptol;
	abs_rcrit =  rtnewt(5, rcrit_poly, rcritmax, .01*rcritmax);
	abs_rcrit += 0.01*rcritmax;

	rcritmax = cmp->bmax + sqrtf_fast((float)3.*invRelTol0*B2*cmp->massinv);
	a[0] = 2.*B3;
	a[1] = -3. * B2 + RelTol0 * cmp->bmax*cmp->bmax * cmp->mass;
	a[2] = -2. * RelTol0 * cmp->bmax * cmp->mass;
	a[3] = RelTol0 * cmp->mass;
	rel_rcrit0 =  rtnewt(3, rcrit_poly, rcritmax, .01*rcritmax);
	rel_rcrit0 += 0.01*rcritmax;

	rcritmax = cmp->bmax + sqrtf_fast((float)3.*invRelTol*B2*cmp->massinv);
	a[0] = 2.*B3;
	a[1] = -3. * B2 + RelTol * cmp->bmax*cmp->bmax * cmp->mass;
	a[2] = -2. * RelTol * cmp->bmax * cmp->mass;
	a[3] = RelTol * cmp->mass;
	rel_rcrit =  rtnewt(3, rcrit_poly, rcritmax, .01*rcritmax);
	rel_rcrit += 0.01*rcritmax;

	/* rcrit is the least accurate of abs_rcrit or rel_rcrit */
	cp->rcrit = (abs_rcrit < rel_rcrit) ? abs_rcrit : rel_rcrit;
	/* if rel_rcrit0 is more accurate, then use it */
	if (rel_rcrit0 > cp->rcrit) cp->rcrit = rel_rcrit0;
#ifdef QUAD
	if (Quad_Ncut && (cmp->ndaughters >= Quad_Ncut)) {
	    R2 = 0.25*cmp->sz*cmp->sz;
	    R4 = R2*R2;
	    B3 = cmp->bmax*(cmp->x2+cmp->y2+cmp->z2)*R2; 	/* upper bound */
	    B4 = (cmp->x4 + cmp->y4 + cmp->z4)*R4;
	    if (!finite(B4) || !finite(B3))
		Error("Bad B3 or B4, B2 = %g, massinv = %g\n", B2, cmp->massinv);
	    rcritmax = abs_rcrit;
	    a[0] = 3.*B4;
	    a[1] = -4. * B3;
	    a[2] = 0.;
	    a[3] = 0.0;
	    a[4] = ptol*cmp->bmax*cmp->bmax;
	    a[5] = -2. * ptol * cmp->bmax;
	    a[6] = ptol;
	    abs_rcrit =  rtnewt(6, rcrit_poly, rcritmax, .01*rcritmax);
	    abs_rcrit += 0.01*rcritmax;

	    rcritmax = rel_rcrit0;
	    a[0] = 3.*B4;
	    a[1] = -4. * B3;
	    a[2] = RelTol0 * cmp->bmax*cmp->bmax * cmp->mass;
	    a[3] = -2. * RelTol0 * cmp->bmax * cmp->mass;
	    a[4] = RelTol0 * cmp->mass;
	    rel_rcrit0 =  rtnewt(4, rcrit_poly, rcritmax, .01*rcritmax);
	    rel_rcrit0 += 0.01*rcritmax;

	    rcritmax = rel_rcrit;
	    a[0] = 3.*B4;
	    a[1] = -4. * B3;
	    a[2] = RelTol * cmp->bmax*cmp->bmax * cmp->mass;
	    a[3] = -2. * RelTol * cmp->bmax * cmp->mass;
	    a[4] = RelTol * cmp->mass;
	    rel_rcrit =  rtnewt(4, rcrit_poly, rcritmax, .01*rcritmax);
	    rel_rcrit += 0.01*rcritmax;

	    cp->rcrit_q = (abs_rcrit < rel_rcrit) ? abs_rcrit : rel_rcrit;
	    if (rel_rcrit0 > cp->rcrit_q) cp->rcrit_q = rel_rcrit0;
	    {
		double t;
		t = (cmp->x2 + cmp->y2 + cmp->z2)/3.0;
		cp->qxx = cmp->x2 - t;
		cp->qxy = cmp->xy;
		cp->qyy = cmp->y2 - t;
		cp->qxz = cmp->xz;
		cp->qyz = cmp->yz;
	    }
	}
#endif

#ifdef HEXA
	if (Hexa_Ncut && (cmp->ndaughters >= Hexa_Ncut)) {
	    R2 = 0.25*cmp->sz*cmp->sz;
	    R4 = R2*R2;
	    B5 = cmp->bmax*(cmp->x4+cmp->y4+cmp->z4)*R4; /* upper bound */
	    if (!finite(B5) || !finite(cmp->B6))
		Error("Bad B5 or B6, B2 = %g, massinv = %g\n", B2, cmp->massinv);
	    rcritmax = abs_rcrit;
	    a[0] = 5.*cmp->B6*R4*R2;
	    a[1] = -6.*B5;
	    a[2] = 0.;
	    a[3] = 0.;
	    a[4] = 0.;
	    a[5] = 0.;
	    a[6] = ptol*cmp->bmax*cmp->bmax;
	    a[7] = -2. * ptol * cmp->bmax;
	    a[8] = ptol;
	    abs_rcrit =  rtnewt(8, rcrit_poly, rcritmax, .001*rcritmax);
	    abs_rcrit += 0.001*rcritmax;

	    rcritmax = rel_rcrit0;
	    a[0] = 5.*cmp->B6*R4*R2;
	    a[1] = -6.*B5;
	    a[2] = 0.;
	    a[3] = 0.;
	    a[4] = RelTol0 * cmp->bmax*cmp->bmax * cmp->mass;
	    a[5] = -2. * RelTol0 * cmp->bmax * cmp->mass;
	    a[6] = RelTol0 * cmp->mass;
	    rel_rcrit0 =  rtnewt(6, rcrit_poly, rcritmax, .001*rcritmax);
	    rel_rcrit0 += 0.001*rcritmax;

	    rcritmax = rel_rcrit;
	    a[0] = 5.*cmp->B6*R4*R2;
	    a[1] = -6.*B5;
	    a[2] = 0.;
	    a[3] = 0.;
	    a[4] = RelTol * cmp->bmax*cmp->bmax * cmp->mass;
	    a[5] = -2. * RelTol * cmp->bmax * cmp->mass;
	    a[6] = RelTol * cmp->mass;
	    rel_rcrit =  rtnewt(6, rcrit_poly, rcritmax, .001*rcritmax);
	    rel_rcrit += 0.001*rcritmax;

	    cp->rcrit_h = (abs_rcrit < rel_rcrit) ? abs_rcrit : rel_rcrit;
	    if (rel_rcrit0 > cp->rcrit_h) cp->rcrit_h = rel_rcrit0;
	    {
		double t, tx, ty, tz, txx, txy, tyy, txz, tyz, tzz;
		t = (cmp->x2 + cmp->y2 + cmp->z2)/3.0;
		cp->qxx = cmp->x2 - t;
		cp->qxy = cmp->xy;
		cp->qyy = cmp->y2 - t;
		cp->qxz = cmp->xz;
		cp->qyz = cmp->yz;
		tx = (cmp->x3 + cmp->xy2 + cmp->xz2)/5.0;
		ty = (cmp->x2y + cmp->y3 + cmp->yz2)/5.0;
		tz = (cmp->x2z + cmp->y2z + cmp->z3)/5.0;
		cp->qxxx = cmp->x3 - 3.0*tx;
		cp->qxxy = cmp->x2y - ty;
		cp->qxyy = cmp->xy2 - tx;
		cp->qyyy = cmp->y3 - 3.0*ty;
		cp->qxxz = cmp->x2z - tz;
		cp->qxyz = cmp->xyz;
		cp->qyyz = cmp->y2z - tz;
		txx = (cmp->x4 + cmp->x2y2 + cmp->x2z2)/7.0;
		txy = (cmp->x3y + cmp->xy3 + cmp->xyz2)/7.0;
		txz = (cmp->x3z + cmp->xy2z + cmp->xz3)/7.0;
		tyy = (cmp->x2y2 + cmp->y4 + cmp->y2z2)/7.0;
		tyz = (cmp->x2yz + cmp->y3z + cmp->yz3)/7.0;
		tzz = (cmp->x2z2 + cmp->y2z2 + cmp->z4)/7.0;
		t = 0.1*(txx + tyy + tzz);
		cp->qxxxx = cmp->x4 - 6.0*(txx - t);
		cp->qxxxy = cmp->x3y - 3.0*txy;
		cp->qxxyy = cmp->x2y2 - (txx + tyy - 2.0*t);
		cp->qxyyy = cmp->xy3 - 3.0*txy;
		cp->qyyyy = cmp->y4 - 6.0*(tyy - t);
		cp->qxxxz = cmp->x3z - 3.0*txz;
		cp->qxxyz = cmp->x2yz - tyz;
		cp->qxyyz = cmp->xy2z - txz;
		cp->qyyyz = cmp->y3z - 3.0*tyz;
	    }
	}
#endif
    } else if(MACtype == BH_MAC) cp->rcrit = cmp->sz*invTol;
    else if (MACtype == BMAX_MAC) cp->rcrit = cmp->bmax*invTol;
    else Error("Bad MAC type (%d)\n", MACtype);

    cp->daughters = cmp->ndaughters;
    Msgf(("Cell: %s\n", PrintCellContents(cp)));
    return cp;
}

/* Use doubles here to avoid catastrophe from roundoff. */
static void rcrit_poly(int n, double r, double *value, double *deriv)
{
    /* Do we care that a[2] and da[1] are zero? */
    double dp = 0.;
    double p;

    p = a[n];
    /* See pg. 149 of Numerical Rec. */
    /* We could unroll it... */
    while(n>0) { dp = dp*r + p; p = p*r + a[--n]; }
    if (!finite(p) || !finite(dp))
      Error("Bad p or dp, p = %g, dp = %g, n = %d\n", p, dp, n);
    *value = p;
    *deriv = dp;
}

/* From Numerical Recipes, rtnewt.c (modified) */
/* It's even more dangerous than the version that NR says is too */
/* dangerous to use...No checking of bounds.  We might just run off */
/* to infinity... */
#define JMAX 30

static double
rtnewt(int n, void (*funcd)(int, double, double *, double *), double x1, double xacc)
{
    int j;
    double df,dx,f,rtn;
    
    rtn=x1;
    for (j=1;j<=JMAX;j++) {
	(*funcd)(n,rtn,&f,&df);
	dx=f/df;
	Msg("rcrit", ("f(%g)=%g, dx=%g\n", rtn, f, -dx));
	rtn -= dx;
	if (fabs(dx) < xacc) return rtn;
    }
    SeriousWarning("Maximum number of iterations exceeded in RTNEWT, %lf %lf %lf %lf, %d\n", 
		   x1, xacc, rtn, dx, n);
    return x1;
}

#if defined(QUAD) || defined(HEXA)
static void
mpole_add_mono(cofmdata *cmp, float m, float x, float y, float z, float Rinv)
{
    float x2, y2, z2;
    float x3, y3, z3;
    float x4, y4, z4;
    float x5, y5, z5;
    float x6, y6, z6;

    x *= Rinv;
    y *= Rinv;
    z *= Rinv;

    x2 = x*x; y2 = y*y; z2 = z*z;
    cmp->x2 += m*x2;
    cmp->xy += m*x*y;
    cmp->y2 += m*y2;
    cmp->xz += m*x*z;
    cmp->yz += m*y*z;
    cmp->z2 += m*z2;

    x3 = x2*x; y3 = y2*y; z3 = z2*z;
    cmp->x3 += m*x3;
    cmp->x2y += m*x2*y;
    cmp->xy2 += m*x*y2;
    cmp->y3 += m*y3;
    cmp->x2z += m*x2*z;
    cmp->xyz += m*x*y*z;
    cmp->y2z += m*y2*z;
    cmp->xz2 += m*x*z2;
    cmp->yz2 += m*y*z2;
    cmp->z3 += m*z3;

    x4 = x3*x; y4 = y3*y; z4 = z3*z;
    cmp->x4 += m*x4;
    cmp->x3y += m*x3*y;
    cmp->x2y2 += m*x2*y2;
    cmp->xy3 += m*x*y3;
    cmp->y4 += m*y4;
    cmp->x3z += m*x3*z;
    cmp->x2yz += m*x2*y*z;
    cmp->xy2z += m*x*y2*z;
    cmp->y3z += m*y3*z;
    cmp->x2z2 += m*x2*z2;
    cmp->xyz2 += m*x*y*z2;
    cmp->y2z2 += m*y2*z2;
    cmp->xz3 += m*x*z3;
    cmp->yz3 += m*y*z3;
    cmp->z4 += m*z4;

    x5 = x4*x; y5 = y4*y; z5 = z4*z;

    cmp->x5 += m*x5;
    cmp->y5 += m*y5;
    cmp->z5 += m*z5;

    x6 = x5*x; y6 = y5*y; z6 = z5*z;

    cmp->B6 += m * (x6 + y6 + z6);
}

static void
mpole_add_shift(cofmdata *cmp, float m, float x, float y, float z, float Rinv, cofmdata *dp)
{
    float x2, y2, z2;
    float x3, y3, z3;
    float x4, y4, z4;
    float x5, y5, z5;
    float x6, y6, z6;

    x *= Rinv;
    y *= Rinv;
    z *= Rinv;

    x2 = x*x; y2 = y*y; z2 = z*z;
    cmp->x2 += 0.25f*dp->x2 + m*x2;
    cmp->xy += 0.25f*dp->xy + m*x*y;
    cmp->y2 += 0.25f*dp->y2 + m*y2;
    cmp->xz += 0.25f*dp->xz + m*x*z;
    cmp->yz += 0.25f*dp->yz + m*y*z;
    cmp->z2 += 0.25f*dp->z2 + m*z2;

    x3 = x2*x; y3 = y2*y; z3 = z2*z;
    cmp->x3 += 0.125f*dp->x3 + 0.75f*dp->x2*x + m*x3;
    cmp->x2y += 0.125f*dp->x2y + 0.5f*dp->xy*x + 0.25f*dp->x2*y + m*x2*y;
    cmp->xy2 += 0.125f*dp->xy2 + 0.25f*dp->y2*x + 0.5f*dp->xy*y + m*x*y2;
    cmp->y3 += 0.125f*dp->y3 + 0.75f*dp->y2*y + m*y3;
    cmp->x2z += 0.125f*dp->x2z + 0.5f*dp->xz*x + 0.25f*dp->x2*z + m*x2*z;
    cmp->xyz += 0.125f*dp->xyz + 0.25f*dp->yz*x + 0.25f*dp->xz*y + 0.25f*dp->xy*z + m*x*y*z;
    cmp->y2z += 0.125f*dp->y2z + 0.5f*dp->yz*y + 0.25f*dp->y2*z + m*y2*z;
    cmp->xz2 += 0.125f*dp->xz2 + 0.25f*dp->z2*x + 0.5f*dp->xz*z + m*x*z2;
    cmp->yz2 += 0.125f*dp->yz2 + 0.25f*dp->z2*y + 0.5f*dp->yz*z + m*y*z2;
    cmp->z3 += 0.125f*dp->z3 + 0.75f*dp->z2*z + m*z3;
	
    x4 = x3*x; y4 = y3*y; z4 = z3*z;
    cmp->x4 += 0.0625f*dp->x4 + 0.5f*dp->x3*x + 1.5f*dp->x2*x2 + m*x4;
    cmp->x3y += 0.0625f*dp->x3y + 0.375f*dp->x2y*x + 0.75f*dp->xy*x2 + 0.125f*dp->x3*y + 0.75f*dp->x2*x*y + m*x3*y;
    cmp->x2y2 += 0.0625f*dp->x2y2 + 0.25f*dp->xy2*x + 0.25f*dp->y2*x2 + 0.25f*dp->x2y*y + dp->xy*x*y + 0.25f*dp->x2*y2 + m*x2*y2;
    cmp->xy3 += 0.0625f*dp->xy3 + 0.125f*dp->y3*x + 0.375f*dp->xy2*y + 0.75f*dp->y2*x*y + 0.75f*dp->xy*y2 + m*x*y3;
    cmp->y4 += 0.0625f*dp->y4 + 0.5f*dp->y3*y + 1.5f*dp->y2*y2 + m*y4;
    cmp->x3z += 0.0625f*dp->x3z + 0.375f*dp->x2z*x + 0.75f*dp->xz*x2 + 0.125f*dp->x3*z + 0.75f*dp->x2*x*z + m*x3*z;
    cmp->x2yz += 0.0625f*dp->x2yz + 0.25f*dp->xyz*x + 0.25f*dp->yz*x2 + 0.125f*dp->x2z*y + 0.5f*dp->xz*x*y + 0.125f*dp->x2y*z + 0.5f*dp->xy*x*z + 0.25f*dp->x2*y*z + m*x2*y*z;
    cmp->xy2z += 0.0625f*dp->xy2z + 0.125f*dp->y2z*x + 0.25f*dp->xyz*y + 0.5f*dp->yz*x*y + 0.25f*dp->xz*y2 + 0.125f*dp->xy2*z + 0.25f*dp->y2*x*z + 0.5f*dp->xy*y*z + m*x*y2*z;
    cmp->y3z += 0.0625f*dp->y3z + 0.375f*dp->y2z*y + 0.75f*dp->yz*y2 + 0.125f*dp->y3*z + 0.75f*dp->y2*y*z + m*y3*z;
    cmp->x2z2 += 0.0625f*dp->x2z2 + 0.25f*dp->xz2*x + 0.25f*dp->z2*x2 + 0.25f*dp->x2z*z + dp->xz*x*z + 0.25f*dp->x2*z2 + m*x2*z2;
    cmp->xyz2 += 0.0625f*dp->xyz2 + 0.125f*dp->yz2*x + 0.125f*dp->xz2*y + 0.25f*dp->z2*x*y + 0.25f*dp->xyz*z + 0.5f*dp->yz*x*z + 0.5f*dp->xz*y*z + 0.25f*dp->xy*z2 + m*x*y*z2;
    cmp->y2z2 += 0.0625f*dp->y2z2 + 0.25f*dp->yz2*y + 0.25f*dp->z2*y2 + 0.25f*dp->y2z*z + dp->yz*y*z + 0.25f*dp->y2*z2 + m*y2*z2;
    cmp->xz3 += 0.0625f*dp->xz3 + 0.125f*dp->z3*x + 0.375f*dp->xz2*z + 0.75f*dp->z2*x*z + 0.75f*dp->xz*z2 + m*x*z3;
    cmp->yz3 += 0.0625f*dp->yz3 + 0.125f*dp->z3*y + 0.375f*dp->yz2*z + 0.75f*dp->z2*y*z + 0.75f*dp->yz*z2 + m*y*z3;
    cmp->z4 += 0.0625f*dp->z4 + 0.5f*dp->z3*z + 1.5f*dp->z2*z2 + m*z4;
    
    x5 = x4*x; y5 = y4*y; z5 = z4*z;
    cmp->x5 += 0.03125f*dp->x5 + 0.3125f*dp->x4*x + 1.25f*dp->x3*x2 + 2.5f*dp->x2*x3 + m*x5;
    cmp->y5 += 0.03125f*dp->y5 + 0.3125f*dp->y4*y + 1.25f*dp->y3*y2 + 2.5f*dp->y2*y3 + m*y5;
    cmp->z5 += 0.03125f*dp->z5 + 0.3125f*dp->z4*z + 1.25f*dp->z3*z2 + 2.5f*dp->z2*z3 + m*z5;

    x6 = x5*x; y6 = y5*y; z6 = z5*z;
    cmp->B6 += 0.015625f*dp->B6;
    cmp->B6 += 0.1875f*dp->x5*x + 0.9375f*dp->x4*x2 + 2.5f*dp->x3*x3 + 3.75f*dp->x2*x4 + m*x6;
    cmp->B6 += 0.1875f*dp->y5*y + 0.9375f*dp->y4*y2 + 2.5f*dp->y3*y3 + 3.75f*dp->y2*y4 + m*y6;
    cmp->B6 += 0.1875f*dp->z5*z + 0.9375f*dp->z4*z2 + 2.5f*dp->z3*z3 + 3.75f*dp->z2*z4 + m*z6;
}
#endif
