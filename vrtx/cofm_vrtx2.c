/*
 * Copyright 1993 Michael S. Warren, John K. Salmon, and
 * Gregoire S. Winckelmans. All Rights Reserved.
 */

/* Also, error estimates are based on norm of vorticity vector (March 93) */
/* Also put the cell at the center of mass of the vorticity norm. */

#include <math.h>
#include <stddef.h>
#include "physics_vrtx.h"
#include "Assert.h"
#include "vop.h"
#include "timers.h"
#include "Msgs.h"
#include "fastflpt.h"
#include "error.h"
#include "protos.h"

#define SQRT3 (1.73205080757)

static void zero_moment(cofm_data *q);
static void add_cell_moment(cofm_data *q, const cofm_data *c);
static void add_body_moment(cofm_data *q, const body *b);

#define NCOEF 7
#define NEWTON_RAPHSON_TOL 0.01

#ifndef HUGE
#define HUGE 3.e30F
#endif

static double rcrit_floor;
static double rcrit_coef[NCOEF];		/* coef of error poly */
static void rcrit_poly(double r, double *value, double *deriv);
static double rtnewt(void (*funcd)(double, double *, double *),  
		    double x1, double xacc);

extern float errtol;
float kc;
float kc2;

void CofmFromDaugh(hcell *hptr, hcell **daughters)
{
    unsigned int i;
    cofm_data *cmp;
    float dxfromcenter[NDIM];
    float xsize, newbmax;
    float *dpos;
    float moment[NDIM];
    float wgtinv, sumwgt, wgt;

    cmp = hptr->ptr;
    assert(cmp);
    assert(Sub_Flags(hptr));
    zero_moment(cmp);
    
    /* First figure out where the new cell's center is. */
    VS(moment, = 0.F);
    sumwgt = 0.F;
    Msgf(("CofmFromDaugh(hp=%s)\n", hcellPrint(hptr)));
    for(i=0; i< (1<<NDIM); i++) {
	if( !daughters[i] )
	    continue;
	if( Sub_Flags(daughters[i]) ){
	    cofm_data *dcmp = (cofm_data *)daughters[i]->ptr;
	    wgt = B0(dcmp);
	    dpos = Pos(dcmp);
	}else{
	    body *bp = (body *)daughters[i]->ptr;
	    float str2;	
	    /* sqrtf_fast may bomb if given 0 as input! */
	    str2 = Dot(Strength(bp), Strength(bp));
	    Msgf(("bp=%p, id=%d, Str[0]=%g\n",
		  bp, bp->ident, Strength(bp)[0]));
	    if( !finite(str2) ){
	      Error("non-finite str2: %g, looking at daughter %d of %s\n"
		    "bp=%p, btab[%d], Pos=%g %g %g, Strength=%g %g %g\n",
		    str2, i, hcellPrint(hptr), 
		    bp, bp->ident, Pos(bp)[0], Pos(bp)[1], Pos(bp)[2],
		    Strength(bp)[0], Strength(bp)[1], Strength(bp)[2]);
	    }
	    wgt = (str2 > 0.F) ? sqrtf_fast(str2) : 0.F;
	    dpos = Pos(bp);
	}

	VV(moment, += wgt*dpos);
	sumwgt += wgt;
    }	
    if( sumwgt > 0.F ){
	wgtinv = recipf(sumwgt);
	VV(Pos(cmp), = wgtinv*moment);
    }else{
	float corner[NDIM];
	float sz;
	CellCorner(hptr->key, corner, &sz);
	sz *= 0.5F;
	VV(Pos(cmp), = sz + corner);
    }
    
    /* Then compute the moments induced by the daughters around the center */
    Msgf(("New cell created with hcell=%s, pos=(%g %g %g)\n",
	  hcellPrint(hptr), cmp->pos[0], cmp->pos[1], cmp->pos[2]));

    for (i = 0; i < (1<<NDIM); i++) {
	if( !daughters[i] )
	    continue;
	if( Sub_Flags(daughters[i]) )
	    add_cell_moment(cmp, (cofm_data *)(daughters[i]->ptr));
	else
	    add_body_moment(cmp, (body *)(daughters[i]->ptr));
    }

    /* Now compute the alternative bound on bmax based on how far */
    /* we have moved from the center? */
    CellCorner(hptr->key, dxfromcenter, &xsize);
    xsize *= 0.5F;
    VV(dxfromcenter, += xsize - Pos(cmp));
    VVS(dxfromcenter, = fabs LPAREN dxfromcenter, RPAREN + xsize);
    newbmax = sqrtf_fast(Dot(dxfromcenter, dxfromcenter));
    if (newbmax < Bmax(cmp))
	Bmax(cmp) = newbmax;
	
}

void CellFromCofm(cell *cell, cofm_data *cofm)
{
    float tmpa[NDIM];
    static float three_tol;
    static float errtol_last;
    float b3, b4;
    float bmaxhalf, rcritmax, rcrit;

    if( errtol_last != errtol && errtol > 0. ){
	three_tol = 3.F/errtol;
	errtol_last = errtol;
    }

    VV(Pos(cell), = Pos(cofm));
    Daughters(cell) = Daughters(cofm);
    
    Bmax(cell) = Bmax(cofm);

    b3 = Bmax(cofm)*B2(cofm);

    /* Avoid div-by-zero by using ? :  */
    b4 = (0.F != B0(cofm)) ? B2(cofm)*B2(cofm)/B0(cofm) : 0.F ; 
    B3(cell) = b3;
    B4(cell) = b4;
    
    if( errtol < 0.0F || b3 < 0.0F || b4 < 0.0F ){
	Shout("Impossible value: errtol=%g, daughters=%d, B0=%g, B2=%g, Bmax=%g, b3=%g, b4=%g\n",
	      errtol, Daughters(cofm), B0(cofm), B2(cofm), Bmax(cofm), b3, b4);
	Shout("Cell position: %g %g %g\n",
	      cofm->pos[0], cofm->pos[1], cofm->pos[2]);
	Error("bye!\n");
    }else if( B0(cofm) == 0.F ){
      /* This case really does correspond to a cell full of zero-strength
	 bodies.  We can "evaluate" this interaction anywhere.  Including
	 inside the cell.  */
      cell->rcrit = 0.F;
    }else if( b3 == 0.0F && b4 == 0.0F ){
      /* How can b3 and b4 both vanish?
	  If the number of non-zero strengths is <= one.  Note that one
	  non-zero str mixed up with a bunch of zeros will have a cm
	  directly on top of the non-zero-str, so b2, b3, etc will vanish
	  even though b0 is non-zero.

	  We used to have 0.0F here, but that allows a bizarre kind of
	  self-interaction between a body and a cell containing that
	  body plus a bunch of zero-strength bodies.  This is more
	  conservative, excluding all self-interactions but possibly
	  excluding a bunch of others too... */
      cell->rcrit = 1.0001F*Bmax(cofm);
    }else if( errtol == 0.0F ){
	cell->rcrit = HUGE;
    }else{
	/* The usual case! */
	rcrit_coef[0] = 3.F*b4;
	rcrit_coef[1] = -4.F * b3;
	rcrit_coef[2] = 0.;
	rcrit_coef[3] = 0.;
	rcrit_coef[4] = errtol*cofm->bmax*cofm->bmax;
	rcrit_coef[5] = -2. * errtol *cofm->bmax;
	rcrit_coef[6] = errtol;
	bmaxhalf = cofm->bmax * 0.5F;
	rcrit_floor = cofm->bmax;
	/* 
	   rcritmax follows from the dipole error estimate, neglecting the 
	   -4*B3 term, and solving exactly for the resulting rcrit.
	   It MUST be an over-estimate of the rcrit for the quad case,
	   which is just what we want.
	   */
	rcritmax = bmaxhalf + sqrtf_fast(bmaxhalf*bmaxhalf + sqrtf_fast(three_tol*B2(cofm))); /* ??? */

	rcrit = rtnewt(rcrit_poly, rcritmax, NEWTON_RAPHSON_TOL*rcritmax);
	rcrit += NEWTON_RAPHSON_TOL*rcritmax;
	/* If B3 and B4 are very small, rcrit approaches bmax. */
	/* allow for roundoff in the Newton-Raphson approx. */
#if 0
	{
	    float err = (1./((rcrit-cofm->bmax)*(rcrit-cofm->bmax))) * 
		(4.*b3/(rcrit*rcrit*rcrit) - 
		 3.*b4/(rcrit*rcrit*rcrit*rcrit));
	    assert(err < errtol);
	}
#endif
	assert(rcrit > cofm->bmax);
	cell->rcrit =  rcrit;
	if( cell->rcrit < kc )
	    cell->rcrit = kc;
    }

    VV(Strength(cell), = Strength(cofm));
    
    VV(Dpole(cell).x, = Dpole(cofm).x);
    VV(Dpole(cell).y, = Dpole(cofm).y);
    VV(Dpole(cell).z, = Dpole(cofm).z);


    VVVV(tmpa, = Qpole(cofm).xx, +Qpole(cofm).yy, +Qpole(cofm).zz);
    VS(tmpa, *= 0.3333333333333F);
    
    VVV(Qpole(cell).xx, = Qpole(cofm).xx, -tmpa);
    VVV(Qpole(cell).yy, = Qpole(cofm).yy, -tmpa);
    VVV(Qpole(cell).zz, = Qpole(cofm).zz, -tmpa);

/*
    VV(Qpole(cell).xx, = Qpole(cofm).xx);
    VV(Qpole(cell).yy, = Qpole(cofm).yy);
    VV(Qpole(cell).zz, = Qpole(cofm).zz);
*/

    VV(Qpole(cell).xy, = Qpole(cofm).xy);
    VV(Qpole(cell).xz, = Qpole(cofm).xz);
    VV(Qpole(cell).yz, = Qpole(cofm).yz);

    
    if( Msg_test(__FILE__) ){
	Msg_do("\nFinishing:\n");
	Msg_do("pos=(%g %g %g)\n",
	       Pos(cofm)[0], Pos(cofm)[1], Pos(cofm)[2]);
	Msg_do("Str=(%g %g %g)\n", 
	       Strength(cofm)[0], Strength(cofm)[1], Strength(cofm)[2]);
	Msg_do("Daughters=%d\n", Daughters(cofm));
	Msg_do("B0=%g\n", B0(cofm));
	Msg_do("B2=%g\n", B2(cofm));
	Msg_do("Bmax=%g\n", Bmax(cofm));
	Msg_do("B3=%g\n", b3);
	Msg_do("B4=%g\n", b4);
	Msg_do("rcrit=%g\n", cell->rcrit);
    }
}

static void zero_moment(cofm_data *q)
{
    
    Daughters(q) = 0;
    
    VS(Strength(q), = (float)0.0);
    
    VS(Dpole(q).x, = (float)0.0);
    VS(Dpole(q).y, = (float)0.0);
    VS(Dpole(q).z, = (float)0.0);
    
    VS(Qpole(q).xx, = (float)0.0);
    VS(Qpole(q).yy, = (float)0.0);
    VS(Qpole(q).zz, = (float)0.0);
    VS(Qpole(q).xy, = (float)0.0);
    VS(Qpole(q).xz, = (float)0.0);
    VS(Qpole(q).yz, = (float)0.0);
    
    
    B0(q) = (float)0.0;
    
    B2(q) = (float)0.0;

    Bmax(q) = (float)0.0;
}

static void add_body_moment(cofm_data *q, const body *b) {
    float B0b;
    float rad2;
    float radius[NDIM];
    float newbmax;

    VVV(radius, = Pos(b), - Pos(q));

    
    VV(Strength(q), += Strength(b));
    
    
    VV(Dpole(q).x, += radius[0]*Strength(b));
    VV(Dpole(q).y, += radius[1]*Strength(b));
    VV(Dpole(q).z, += radius[2]*Strength(b));
    
    
    VV(Qpole(q).xx, += radius[0]*radius[0]*Strength(b));
    VV(Qpole(q).yy, += radius[1]*radius[1]*Strength(b));
    VV(Qpole(q).zz, += radius[2]*radius[2]*Strength(b));
    VV(Qpole(q).xy, += radius[0]*radius[1]*Strength(b));
    VV(Qpole(q).xz, += radius[0]*radius[2]*Strength(b));
    VV(Qpole(q).yz, += radius[1]*radius[2]*Strength(b));
    
    
    B0b = sqrtf_fast(Dot(Strength(b), Strength(b)));

    B0(q) += B0b;
    

    rad2 = Dot(radius, radius);    
    newbmax = sqrtf_fast(rad2);
    if( newbmax > Bmax(q) )
	Bmax(q) = newbmax;

    B2(q) += rad2*B0b;
    
    Daughters(q) ++;
}

static void 
    add_cell_moment(cofm_data *q, const cofm_data *c)
{
    float rad2, newbmax;
    float radius[NDIM];
    
    VVV(radius, = Pos(c), - Pos(q));

    VV(Strength(q), += Strength(c));
    
    
    VVV(Dpole(q).x, += Dpole(c).x, + radius[0]*Strength(c));
    VVV(Dpole(q).y, += Dpole(c).y, + radius[1]*Strength(c));
    VVV(Dpole(q).z, += Dpole(c).z, + radius[2]*Strength(c));
    
    
    VVVV(Qpole(q).xx, += Qpole(c).xx
	 , +  (float)2.0*radius[0]*Dpole(c).x
	 , +  radius[0]*radius[0]*Strength(c));
    
    VVVV(Qpole(q).yy, += Qpole(c).yy
	 , +  (float)2.0*radius[1]*Dpole(c).y
	 , +  radius[1]*radius[1]*Strength(c));
    
    VVVV(Qpole(q).zz, += Qpole(c).zz
	 , +  (float)2.0*radius[2]*Dpole(c).z
	 , +  radius[2]*radius[2]*Strength(c));
    
    
    VVVVV(Qpole(q).xy, += Qpole(c).xy
	  , +  radius[0]*Dpole(c).y
	  , +  radius[1]*Dpole(c).x
	  , +  radius[0]*radius[1]*Strength(c));
    
    VVVVV(Qpole(q).xz, += Qpole(c).xz
	  , +  radius[0]*Dpole(c).z
	  , +  radius[2]*Dpole(c).x
	  , +  radius[0]*radius[2]*Strength(c));
    
    VVVVV(Qpole(q).yz, += Qpole(c).yz
	  , +  radius[1]*Dpole(c).z
	  , +  radius[2]*Dpole(c).y
	  , +  radius[1]*radius[2]*Strength(c));
    
    
    
    
    B0(q) += B0(c);
    
    
    rad2 = Dot(radius, radius);    
    newbmax = sqrtf_fast(rad2) + Bmax(c);
    if( newbmax > Bmax(q) )
	Bmax(q) = newbmax;

    B2(q) += B2(c)
	   +rad2*B0(c);
    
    Daughters(q) += Daughters(c);
    
}

/* Use doubles here to avoid catastrophe from roundoff. */
static void rcrit_poly(double r, double *value, double *deriv){
    /* Do we care that rcrit_coef[2] and drcrit_coef[1] are zero? */
    double p = rcrit_coef[NCOEF-1];
    double dp = 0.;
    int n = NCOEF-1;

    /* See pg. 149 of Numerical Rec. */
    /* We could unroll it... */
    while(n>0) { dp = dp*r + p; p = p*r + rcrit_coef[--n]; }
    if (!finite(p) || !finite(dp))
      Error("Bad p or dp, p = %g, dp = %g, n = %d\n", p, dp, n);
    *value = p;
    *deriv = dp;
}

/* From Numerical Recipes, rtnewt.c (modified) */
/* It's even more dangerous than the version that NR says is too */
/* dangerous to use...No checking of bounds.  We might just run off */
/* to infinity... */
#define JMAX 20

static double
    rtnewt(void (*funcd)(double, double *, double *),
	   double x1, double xacc)
{
    int j;
    double df,dx,f,rtn;
    
    rtn=x1;
    for (j=1;j<=JMAX;j++) {
	(*funcd)(rtn,&f,&df);
	dx=f/df;
	Msg("rcrit", ("f(%g)=%g, dx=%g\n", rtn, f, -dx));
	rtn -= dx;
	if( rtn < rcrit_floor ){
	    rtn += dx;
	    dx = 0.5 * (rtn - rcrit_floor);
	    rtn -= dx;
	    Warning("Newton-Raphson iteration overshoot!\n");
	}
	if (fabs(dx) < xacc) return rtn;
    }
    Error("Maximum number of iterations exceeded in RTNEWT\n");
}

#undef JMAX
