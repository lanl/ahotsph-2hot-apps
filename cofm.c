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
	cp = ChnAlloc(&Tree->cell4chn);
    } else if (Quad_Ncut && (cmp->ndaughters >= Quad_Ncut)) {
	cp = ChnAlloc(&Tree->cell2chn);
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
	double R2, R4, B4, B5, B6;
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
	    B4 = (cmp->x4 + cmp->y4 + cmp->z4 + 2*cmp->x2y2 + 2*cmp->x2z2 + 2*cmp->y2z2)*R4;
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
	    B4 = (cmp->x4 + cmp->y4 + cmp->z4 + 2.0f*cmp->x2y2 + 2.0f*cmp->x2z2 + 2.0f*cmp->y2z2)*R4;
	    B5 = cmp->bmax*B4; /* upper bound */
	    B6 = (cmp->x6+cmp->y6+cmp->z6+3.0f*cmp->x4y2+3.0f*cmp->x2y4+3.0f*cmp->x4z2+6.0f*cmp->x2y2z2+3.0f*cmp->y4z2+3.0f*cmp->x2z4+3.0f*cmp->y2z4)*R2*R4;
	    if (!finite(B5) || !finite(B6))
		Error("Bad B5 or B6, B2 = %g, massinv = %g\n", B2, cmp->massinv);
	    rcritmax = abs_rcrit;
	    a[0] = 5.*B6;
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
	    a[0] = 5.*B6;
	    a[1] = -6.*B5;
	    a[2] = 0.;
	    a[3] = 0.;
	    a[4] = RelTol0 * cmp->bmax*cmp->bmax * cmp->mass;
	    a[5] = -2. * RelTol0 * cmp->bmax * cmp->mass;
	    a[6] = RelTol0 * cmp->mass;
	    rel_rcrit0 =  rtnewt(6, rcrit_poly, rcritmax, .001*rcritmax);
	    rel_rcrit0 += 0.001*rcritmax;

	    rcritmax = rel_rcrit;
	    a[0] = 5.*B6;
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
	    Msgf(("Cell4: %s\n", PrintCellContents4(cp)));
	}
#endif
    } else if(MACtype == BH_MAC) cp->rcrit = cmp->sz*invTol;
    else if (MACtype == BMAX_MAC) cp->rcrit = cmp->bmax*invTol;
    else Error("Bad MAC type (%d)\n", MACtype);

    cp->daughters = cmp->ndaughters;
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
    cmp->x4y += m*x4*y;
    cmp->x4z += m*x4*z;
    cmp->x3y2 += m*x3*y2;
    cmp->x3yz += m*x3*y*z;
    cmp->x3z2 += m*x3*z2;
    cmp->x2y3 += m*x2*y3;
    cmp->x2y2z += m*x2*y2*z;
    cmp->x2yz2 += m*x2*y*z2;
    cmp->x2z3 += m*x2*z3;
    cmp->xy4 += m*x*y4;
    cmp->xy3z += m*x*y3*z;
    cmp->xy2z2 += m*x*y2*z2;
    cmp->xyz3 += m*x*y*z3;
    cmp->xz4 += m*x*z4;
    cmp->y5 += m*y5;
    cmp->y4z += m*y4*z;
    cmp->y3z2 += m*y3*z2;
    cmp->y2z3 += m*y2*z3;
    cmp->yz4 += m*y*z4;
    cmp->z5 += m*z5;

    x6 = x5*x; y6 = y5*y; z6 = z5*z;
    cmp->x6 += m*x6;
    cmp->x4y2 += m*x4*y2;
    cmp->x4z2 += m*x4*z2;
    cmp->x2y4 += m*x2*y4;
    cmp->x2y2z2 += m*x2*y2*z2;
    cmp->x2z4 += m*x2*z4;
    cmp->y6 += m*y6;
    cmp->y4z2 += m*y4*z2;
    cmp->y2z4 += m*y2*z4;
    cmp->z6 += m*z6;
}

static void
mpole_add_shift(cofmdata *cmp, float m, float x, float y, float z, float Rinv, cofmdata *dp)
{
    float x2, y2, z2;
    float x3, y3, z3;
    float x4, y4, z4;
    float x5, y5, z5;
    float x6, y6, z6;

    x *= 2.0f*Rinv;
    y *= 2.0f*Rinv;
    z *= 2.0f*Rinv;

    x2 = x*x; y2 = y*y; z2 = z*z;

    cmp->x2 += (dp->x2 + m*x2)/4.0f;
    cmp->xy += (dp->xy + m*x*y)/4.0f;
    cmp->xz += (dp->xz + m*x*z)/4.0f;
    cmp->y2 += (dp->y2 + m*y2)/4.0f;
    cmp->yz += (dp->yz + m*y*z)/4.0f;
    cmp->z2 += (dp->z2 + m*z2)/4.0f;

    x3 = x2*x; y3 = y2*y; z3 = z2*z;

    cmp->x3 += (dp->x3 + 3*dp->x2*x + m*x3)/8.0f;
    cmp->x2y += (dp->x2y + 2*dp->xy*x + dp->x2*y + m*x2*y)/8.0f;
    cmp->x2z += (dp->x2z + 2*dp->xz*x + dp->x2*z + m*x2*z)/8.0f;
    cmp->xy2 += (dp->xy2 + dp->y2*x + 2*dp->xy*y + m*x*y2)/8.0f;
    cmp->xyz += (dp->xyz + dp->yz*x + dp->xz*y + dp->xy*z + m*x*y*z)/8.0f;
    cmp->xz2 += (dp->xz2 + dp->z2*x + 2*dp->xz*z + m*x*z2)/8.0f;
    cmp->y3 += (dp->y3 + 3*dp->y2*y + m*y3)/8.0f;
    cmp->y2z += (dp->y2z + 2*dp->yz*y + dp->y2*z + m*y2*z)/8.0f;
    cmp->yz2 += (dp->yz2 + dp->z2*y + 2*dp->yz*z + m*y*z2)/8.0f;
    cmp->z3 += (dp->z3 + 3*dp->z2*z + m*z3)/8.0f;

    x4 = x3*x; y4 = y3*y; z4 = z3*z;

    cmp->x4 += (dp->x4 + 4*dp->x3*x + 6*dp->x2*x2 + m*x4)/16.0f;
    cmp->x3y += (dp->x3y + 3*dp->x2y*x + 3*dp->xy*x2 + dp->x3*y + 3*dp->x2*x*y + m*x3*y)/16.0f;
    cmp->x3z += (dp->x3z + 3*dp->x2z*x + 3*dp->xz*x2 + dp->x3*z + 3*dp->x2*x*z + m*x3*z)/16.0f;
    cmp->x2y2 += (dp->x2y2 + 2*dp->xy2*x + dp->y2*x2 + 2*dp->x2y*y + 4*dp->xy*x*y + dp->x2*y2 + m*x2*y2)/16.0f;
    cmp->x2yz += (dp->x2yz + 2*dp->xyz*x + dp->yz*x2 + dp->x2z*y + 2*dp->xz*x*y + dp->x2y*z + 2*dp->xy*x*z + dp->x2*y*z  + m*x2*y*z)/16.0f;
    cmp->x2z2 += (dp->x2z2 + 2*dp->xz2*x + dp->z2*x2 + 2*dp->x2z*z + 4*dp->xz*x*z + dp->x2*z2 + m*x2*z2)/16.0f;
    cmp->xy3 += (dp->xy3 + dp->y3*x + 3*dp->xy2*y + 3*dp->y2*x*y + 3*dp->xy*y2 + m*x*y3)/16.0f;
    cmp->xy2z += (dp->xy2z + dp->y2z*x + 2*dp->xyz*y + 2*dp->yz*x*y + dp->xz*y2 + dp->xy2*z + dp->y2*x*z + 2*dp->xy*y*z + m*x*y2*z)/16.0f;
    cmp->xyz2 += (dp->xyz2 + dp->yz2*x + dp->xz2*y + dp->z2*x*y + 2*dp->xyz*z + 2*dp->yz*x*z + 2*dp->xz*y*z + dp->xy*z2 + m*x*y*z2)/16.0f;
    cmp->xz3 += (dp->xz3 + dp->z3*x + 3*dp->xz2*z + 3*dp->z2*x*z + 3*dp->xz*z2 + m*x*z3)/16.0f;
    cmp->y4 += (dp->y4 + 4*dp->y3*y + 6*dp->y2*y2 + m*y4)/16.0f;
    cmp->y3z += (dp->y3z + 3*dp->y2z*y + 3*dp->yz*y2 + dp->y3*z + 3*dp->y2*y*z + m*y3*z)/16.0f;
    cmp->y2z2 += (dp->y2z2 + 2*dp->yz2*y + dp->z2*y2 + 2*dp->y2z*z + 4*dp->yz*y*z + dp->y2*z2 + m*y2*z2)/16.0f;
    cmp->yz3 += (dp->yz3 + dp->z3*y + 3*dp->yz2*z + 3*dp->z2*y*z + 3*dp->yz*z2 + m*y*z3)/16.0f;
    cmp->z4 += (dp->z4 + 4*dp->z3*z + 6*dp->z2*z2 + m*z4)/16.0f;

    x5 = x4*x; y5 = y4*y; z5 = z4*z;

    cmp->x5 += (dp->x5 + 5*dp->x4*x + 10*dp->x3*x2 + 10*dp->x2*x3 + m*x5)/32.0f;
    cmp->x4y += (dp->x4y + 4*dp->x3y*x + 6*dp->x2y*x2 + 4*dp->xy*x3 + dp->x4*y + 4*dp->x3*x*y + 6*dp->x2*x2*y + m*x4*y)/32.0f;
    cmp->x4z += (dp->x4z + 4*dp->x3z*x + 6*dp->x2z*x2 + 4*dp->xz*x3 + dp->x4*z + 4*dp->x3*x*z + 6*dp->x2*x2*z + m*x4*z)/32.0f;
    cmp->x3y2 += (dp->x3y2 + 3*dp->x2y2*x + 3*dp->xy2*x2 + dp->y2*x3 + 2*dp->x3y*y + 6*dp->x2y*x*y + 6*dp->xy*x2*y + dp->x3*y2 + 3*dp->x2*x*y2 + m*x3*y2)/32.0f;
    cmp->x3yz += (dp->x3yz + 3*dp->x2yz*x + 3*dp->xyz*x2 + dp->yz*x3 + dp->x3z*y + 3*dp->x2z*x*y + 3*dp->xz*x2*y + dp->x3y*z + 3*dp->x2y*x*z + 3*dp->xy*x2*z + dp->x3*y*z + 3*dp->x2*x*y*z + m*x3*y*z)/32.0f;
    cmp->x3z2 += (dp->x3z2 + 3*dp->x2z2*x + 3*dp->xz2*x2 + dp->z2*x3 + 2*dp->x3z*z + 6*dp->x2z*x*z + 6*dp->xz*x2*z + dp->x3*z2 + 3*dp->x2*x*z2 + m*x3*z2)/32.0f;
    cmp->x2y3 += (dp->x2y3 + 2*dp->xy3*x + dp->y3*x2 + 3*dp->x2y2*y + 6*dp->xy2*x*y + 3*dp->y2*x2*y + 3*dp->x2y*y2 + 6*dp->xy*x*y2 + dp->x2*y3 + m*x2*y3)/32.0f;
    cmp->x2y2z += (dp->x2y2z + 2*dp->xy2z*x + dp->y2z*x2 + 2*dp->x2yz*y + 4*dp->xyz*x*y + 2*dp->yz*x2*y + dp->x2z*y2 + 2*dp->xz*x*y2 + dp->x2y2*z + 2*dp->xy2*x*z + dp->y2*x2*z + 2*dp->x2y*y*z + 4*dp->xy*x*y*z + dp->x2*y2*z + m*x2*y2*z)/32.0f;
    cmp->x2yz2 += (dp->x2yz2 + 2*dp->xyz2*x + dp->yz2*x2 + dp->x2z2*y + 2*dp->xz2*x*y + dp->z2*x2*y + 2*dp->x2yz*z + 4*dp->xyz*x*z + 2*dp->yz*x2*z + 2*dp->x2z*y*z + 4*dp->xz*x*y*z + dp->x2y*z2 + 2*dp->xy*x*z2 + dp->x2*y*z2 + m*x2*y*z2)/32.0f;
    cmp->x2z3 += (dp->x2z3 + 2*dp->xz3*x + dp->z3*x2 + 3*dp->x2z2*z + 6*dp->xz2*x*z + 3*dp->z2*x2*z + 3*dp->x2z*z2 + 6*dp->xz*x*z2 + dp->x2*z3 + m*x2*z3)/32.0f;
    cmp->xy4 += (dp->xy4 + dp->y4*x + 4*dp->xy3*y + 4*dp->y3*x*y + 6*dp->xy2*y2 + 6*dp->y2*x*y2 + 4*dp->xy*y3 + m*x*y4)/32.0f;
    cmp->xy3z += (dp->xy3z + dp->y3z*x + 3*dp->xy2z*y + 3*dp->y2z*x*y + 3*dp->xyz*y2 + 3*dp->yz*x*y2 + dp->xz*y3 + dp->xy3*z + dp->y3*x*z + 3*dp->xy2*y*z + 3*dp->y2*x*y*z + 3*dp->xy*y2*z + m*x*y3*z)/32.0f;
    cmp->xy2z2 += (dp->xy2z2 + dp->y2z2*x + 2*dp->xyz2*y + 2*dp->yz2*x*y + dp->xz2*y2 + dp->z2*x*y2 + 2*dp->xy2z*z + 2*dp->y2z*x*z + 4*dp->xyz*y*z + 4*dp->yz*x*y*z + 2*dp->xz*y2*z + dp->xy2*z2 + dp->y2*x*z2 + 2*dp->xy*y*z2 + m*x*y2*z2)/32.0f;
    cmp->xyz3 += (dp->xyz3 + dp->yz3*x + dp->xz3*y + dp->z3*x*y + 3*dp->xyz2*z + 3*dp->yz2*x*z + 3*dp->xz2*y*z + 3*dp->z2*x*y*z + 3*dp->xyz*z2 + 3*dp->yz*x*z2 + 3*dp->xz*y*z2 + dp->xy*z3 + m*x*y*z3)/32.0f;
    cmp->xz4 += (dp->xz4 + dp->z4*x + 4*dp->xz3*z + 4*dp->z3*x*z + 6*dp->xz2*z2 + 6*dp->z2*x*z2 + 4*dp->xz*z3 + m*x*z4)/32.0f;
    cmp->y5 += (dp->y5 + 5*dp->y4*y + 10*dp->y3*y2 + 10*dp->y2*y3 + m*y5)/32.0f;
    cmp->y4z += (dp->y4z + 4*dp->y3z*y + 6*dp->y2z*y2 + 4*dp->yz*y3 + dp->y4*z + 4*dp->y3*y*z + 6*dp->y2*y2*z + m*y4*z)/32.0f;
    cmp->y3z2 += (dp->y3z2 + 3*dp->y2z2*y + 3*dp->yz2*y2 + dp->z2*y3 + 2*dp->y3z*z + 6*dp->y2z*y*z + 6*dp->yz*y2*z + dp->y3*z2 + 3*dp->y2*y*z2 + m*y3*z2)/32.0f;
    cmp->y2z3 += (dp->y2z3 + 2*dp->yz3*y + dp->z3*y2 + 3*dp->y2z2*z + 6*dp->yz2*y*z + 3*dp->z2*y2*z + 3*dp->y2z*z2 + 6*dp->yz*y*z2 + dp->y2*z3 + m*y2*z3)/32.0f;
    cmp->yz4 += (dp->yz4 + dp->z4*y + 4*dp->yz3*z + 4*dp->z3*y*z + 6*dp->yz2*z2 + 6*dp->z2*y*z2 + 4*dp->yz*z3 + m*y*z4)/32.0f;
    cmp->z5 += (dp->z5 + 5*dp->z4*z + 10*dp->z3*z2 + 10*dp->z2*z3 + m*z5)/32.0f;

    x6 = x5*x; y6 = y5*y; z6 = z5*z;

    cmp->x6 += (dp->x6 + 6*dp->x5*x + 15*dp->x4*x2 + 20*dp->x3*x3 + 15*dp->x2*x4 + m*x6)/64.0f;
    cmp->x4y2 += (dp->x4y2 + 4*dp->x3y2*x + 6*dp->x2y2*x2 + 4*dp->xy2*x3 + dp->y2*x4 + 2*dp->x4y*y + 8*dp->x3y*x*y + 12*dp->x2y*x2*y + 8*dp->xy*x3*y + dp->x4*y2 + 4*dp->x3*x*y2 + 6*dp->x2*x2*y2 + m*x4*y2)/64.0f;
    cmp->x4z2 += (dp->x4z2 + 4*dp->x3z2*x + 6*dp->x2z2*x2 + 4*dp->xz2*x3 + dp->z2*x4 + 2*dp->x4z*z + 8*dp->x3z*x*z + 12*dp->x2z*x2*z + 8*dp->xz*x3*z + dp->x4*z2 + 4*dp->x3*x*z2 + 6*dp->x2*x2*z2 + m*x4*z2)/64.0f;
    cmp->x2y4 += (dp->x2y4 + 2*dp->xy4*x + dp->y4*x2 + 4*dp->x2y3*y + 8*dp->xy3*x*y + 4*dp->y3*x2*y + 6*dp->x2y2*y2 + 12*dp->xy2*x*y2 + 6*dp->y2*x2*y2 + 4*dp->x2y*y3 + 8*dp->xy*x*y3 + dp->x2*y4 + m*x2*y4)/64.0f;
    cmp->x2y2z2 += (dp->x2y2z2 + 2*dp->xy2z2*x + dp->y2z2*x2 + 2*dp->x2yz2*y + 4*dp->xyz2*x*y + 2*dp->yz2*x2*y + dp->x2z2*y2 + 2*dp->xz2*x*y2 + dp->z2*x2*y2 + 2*dp->x2y2z*z + 4*dp->xy2z*x*z + 2*dp->y2z*x2*z + 4*dp->x2yz*y*z + 8*dp->xyz*x*y*z + 4*dp->yz*x2*y*z + 2*dp->x2z*y2*z + 4*dp->xz*x*y2*z + dp->x2y2*z2 + 2*dp->xy2*x*z2 + dp->y2*x2*z2 + 2*dp->x2y*y*z2 + 4*dp->xy*x*y*z2 + dp->x2*y2*z2 + m*x2*y2*z2)/64.0f;
    cmp->x2z4 += (dp->x2z4 + 2*dp->xz4*x + dp->z4*x2 + 4*dp->x2z3*z + 8*dp->xz3*x*z + 4*dp->z3*x2*z + 6*dp->x2z2*z2 + 12*dp->xz2*x*z2 + 6*dp->z2*x2*z2 + 4*dp->x2z*z3 + 8*dp->xz*x*z3 + dp->x2*z4 + m*x2*z4)/64.0f;
    cmp->y6 += (dp->y6 + 6*dp->y5*y + 15*dp->y4*y2 + 20*dp->y3*y3 + 15*dp->y2*y4 + m*y6)/64.0f;
    cmp->y4z2 += (dp->y4z2 + 4*dp->y3z2*y + 6*dp->y2z2*y2 + 4*dp->yz2*y3 + dp->z2*y4 + 2*dp->y4z*z + 8*dp->y3z*y*z + 12*dp->y2z*y2*z + 8*dp->yz*y3*z + dp->y4*z2 + 4*dp->y3*y*z2 + 6*dp->y2*y2*z2 + m*y4*z2)/64.0f;
    cmp->y2z4 += (dp->y2z4 + 2*dp->yz4*y + dp->z4*y2 + 4*dp->y2z3*z + 8*dp->yz3*y*z + 4*dp->z3*y2*z + 6*dp->y2z2*z2 + 12*dp->yz2*y*z2 + 6*dp->z2*y2*z2 + 4*dp->y2z*z3 + 8*dp->yz*y*z3 + dp->y2*z4 + m*y2*z4)/64.0f;
    cmp->z6 += (dp->z6 + 6*dp->z5*z + 15*dp->z4*z2 + 20*dp->z3*z3 + 15*dp->z2*z4 + m*z6)/64.0f;
}
#endif
