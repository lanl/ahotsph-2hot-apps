/*
 * Copyright 1993 Michael S. Warren, John K. Salmon, and
 * Gregoire S. Winckelmans. All Rights Reserved.
 */

/* Also, error estimates are based on norm of vorticity vector (March 93) */
/* Also put the cell at the center of mass of the vorticity norm. */

#include <math.h>
#include <stddef.h>
#include "physics_panel.h"
#include "physics_generic.h"
#include "Assert.h"
#include "vop.h"
#include "timers.h"
#include "Msgs.h"
#include "fastflpt.h"

static void zero_moment(cofm_data *q);
static void add_cell_moment(cofm_data *q, const cofm_data *c);
static void add_body_moment(cofm_data *q, const body *b);

#define NCOEF 7
#define NEWTON_RAPHSON_TOL 0.01
static double rcrit_coef[NCOEF];		/* coef of error poly */
static void rcrit_poly(double r, double *value, double *deriv);
static double rtnewt(void (*funcd)(double, double *, double *),  
		    double x1, double xacc);

static float three_tol;

float errtol;


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
    
    /* First figure out where the new cell's centroid is. */
    VS(moment, = 0.F);
    sumwgt = 0.F;
    for(i=0; i< (1<<NDIM); i++) {
	if( !daughters[i] )
	    continue;
	if( Sub_Flags(daughters[i]) ){
	    cofm_data *dcmp = (cofm_data *)daughters[i]->ptr;
	    wgt = B0(dcmp);
	    dpos = Pos(dcmp);
	}else{
	    body *bp = (body *)daughters[i]->ptr;
	    wgt = fabs(Sigma(bp));
	    dpos = Pos(bp);
	}

	VV(moment, += wgt*dpos);
	sumwgt += wgt;
    }	

    if( sumwgt == 0. ){
	/* All the individual weights are identically zero!! */
	/* It doesn't matter what we do, but to be safe, we */
	/* put the position at the center of the cell */
	float sz;
	CellCorner(hptr->key, Pos(cmp), &sz);
	sz *= 0.5F;
	VS(Pos(cmp), += sz);
    }else{
	wgtinv = recipf(sumwgt);
	VV(Pos(cmp), = wgtinv*moment);
    }
    
    /* Then compute the moments induced by the daughters around the center */
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
    newbmax = sqrtf_fast(Dot(dxfromcenter, dxfromcenter)) + cmp->maxsize;
    if (newbmax < Bmax(cmp))
	Bmax(cmp) = newbmax;
	
}

void CellFromCofm(cell *cell, cofm_data *cofm)
{
    float trace;
    float bmaxhalf, rcrit, rcritmax;
    static int first_time = 1;

    if( first_time ){
	three_tol = 3.F/errtol;
	first_time = 0;
    }

    VV(Pos(cell), = Pos(cofm));
    Daughters(cell) = Daughters(cofm);
    
    Bmax(cell) = Bmax(cofm);

    B3(cell) = Bmax(cofm)*B2(cofm);

    /* Avoid div-by-zero by using ? :  */
    B4(cell) = (0.F != B0(cofm)) ? B2(cofm)*B2(cofm)/B0(cofm) : 0.F ; 
    
    if( errtol > 0.F ){
	rcrit_coef[0] = 3.*B4(cell);
	rcrit_coef[1] = -4. * B3(cell);
	rcrit_coef[2] = 0.;
	rcrit_coef[3] = 0.;
	rcrit_coef[4] = errtol*cofm->bmax*cofm->bmax;
	rcrit_coef[5] = -2. * errtol *cofm->bmax;
	rcrit_coef[6] = errtol;
	bmaxhalf = cofm->bmax * 0.5F;
	rcritmax = bmaxhalf + sqrtf_fast(bmaxhalf*bmaxhalf + sqrtf_fast(three_tol*B2(cofm))); /* ??? */
	rcrit = rtnewt(rcrit_poly, rcritmax, NEWTON_RAPHSON_TOL*rcritmax);
	rcrit += NEWTON_RAPHSON_TOL*rcritmax;
	/* If B3 and B4 are very small, rcrit approaches bmax. */
	/* allow for roundoff in the Newton-Raphson approx. */
#if 0
	{
	    float err = (1./((rcrit-cofm->bmax)*(rcrit-cofm->bmax))) * 
		(4.*B3(cell)/(rcrit*rcrit*rcrit) - 
		 3.*B4(cell)/(rcrit*rcrit*rcrit*rcrit));
	    assert(err < errtol);
	}
#endif
	assert(rcrit > cofm->bmax);
	cell->rcrit2 =  rcrit*rcrit;
    }else{
	cell->rcrit2 = HUGE;
    }

    Sigma(cell) = Sigma(cofm);
    
    Dpole(cell).x = Dpole(cofm).x;
    Dpole(cell).y = Dpole(cofm).y;
    Dpole(cell).z = Dpole(cofm).z;

    trace = (1.F/3.F) * (Qpole(cofm).xx + Qpole(cofm).yy + Qpole(cofm).zz);
    Qpole(cell).xx = Qpole(cofm).xx - trace;
    Qpole(cell).yy = Qpole(cofm).yy - trace;
    Qpole(cell).zz = Qpole(cofm).zz - trace;
    Qpole(cell).xy = Qpole(cofm).xy;
    Qpole(cell).xz = Qpole(cofm).xz;
    Qpole(cell).yz = Qpole(cofm).yz;
    
    if( Msg_test(__FILE__) ){
	Msg_do("\nFinishing:\n");
	Msg_do("pos=(%g %g %g)\n",
	       Pos(cofm)[0], Pos(cofm)[1], Pos(cofm)[2]);
	Msg_do("Daughters=%d\n", Daughters(cofm));
	Msg_do("B0=%g\n", B0(cofm));
	Msg_do("B2=%g\n", B2(cofm));
	Msg_do("Bmax=%g\n", Bmax(cofm));
        }
}

static void zero_moment(cofm_data *q)
{
    
    Daughters(q) = 0;
    
    Sigma(q) = (float)0.0;
    
    Dpole(q).x = Dpole(q).y = Dpole(q).z = (float)0.0;

    Qpole(q).xx = Qpole(q).yy = Qpole(q).zz 
	= Qpole(q).xy = Qpole(q).xz = Qpole(q).yz = (float)0.0;
    
    B0(q) = (float)0.0;
    
    B2(q) = (float)0.0;

    q->maxsize = (float)0.0;
    Bmax(q) = (float)0.0;
}

static void add_body_moment(cofm_data *q, const body *b) {
/* A panel already has a multipole expansion about its centroid: non-zero
   monopole and quadrapole terms; zero dipole term since panel position
   is located at panel centroid */
    float B0b;
    float rad2;
    float radius[NDIM];
    float newbmax;
    float magstr;

    VVV(radius, = Pos(b), - Pos(q));

    Sigma(q) += Ip(b)*Sigma(b);
    
    Dpole(q).x += radius[0]*Ip(b)*Sigma(b);
    Dpole(q).y += radius[1]*Ip(b)*Sigma(b);
    Dpole(q).z += radius[2]*Ip(b)*Sigma(b);
    
    Qpole(q).xx += (Ixxp(b)+radius[0]*radius[0]*Ip(b))*Sigma(b);    
    Qpole(q).yy += (Iyyp(b)+radius[1]*radius[1]*Ip(b))*Sigma(b);
    Qpole(q).zz += (Izzp(b)+radius[2]*radius[2]*Ip(b))*Sigma(b);
    Qpole(q).xy += (Ixyp(b)+radius[0]*radius[1]*Ip(b))*Sigma(b);
    Qpole(q).xz += (Ixzp(b)+radius[0]*radius[2]*Ip(b))*Sigma(b);
    Qpole(q).yz += (Iyzp(b)+radius[1]*radius[2]*Ip(b))*Sigma(b);
    
    
    magstr = fabs(Sigma(b));
    B0b = Ip(b) * magstr;

    B0(q) += B0b;
    

    rad2 = Dot(radius, radius);    
    newbmax = sqrtf_fast(rad2) + Size(b); /* a panel has a size too ! */
    if( newbmax > Bmax(q) )
	Bmax(q) = newbmax;

    B2(q) += rad2*B0b + (Ixxp(b) + Iyyp(b) + Izzp(b))*magstr;

    if( q->maxsize < Size(b) )
      q->maxsize = Size(b);
    
    Daughters(q) ++;
}

static void 
    add_cell_moment(cofm_data *q, const cofm_data *c)
{
    float rad2, newbmax;
    float radius[NDIM];
    
    VVV(radius, = Pos(c), - Pos(q));

    Sigma(q) += Sigma(c);
    
    
    Dpole(q).x += Dpole(c).x + radius[0]*Sigma(c);
    Dpole(q).y += Dpole(c).y + radius[1]*Sigma(c);
    Dpole(q).z += Dpole(c).z + radius[2]*Sigma(c);

    
    
    Qpole(q).xx += Qpole(c).xx
	  +  (float)2.0*radius[0]*Dpole(c).x
	  +  radius[0]*radius[0]*Sigma(c);
    
    Qpole(q).yy += Qpole(c).yy
	  +  (float)2.0*radius[1]*Dpole(c).y
	  +  radius[1]*radius[1]*Sigma(c);
    
    Qpole(q).zz += Qpole(c).zz
	  +  (float)2.0*radius[2]*Dpole(c).z
	  +  radius[2]*radius[2]*Sigma(c);
    
    
    Qpole(q).xy += Qpole(c).xy
	  +  radius[0]*Dpole(c).y
	  +  radius[1]*Dpole(c).x
	  +  radius[0]*radius[1]*Sigma(c);
    
    Qpole(q).xz += Qpole(c).xz
	  +  radius[0]*Dpole(c).z
	  +  radius[2]*Dpole(c).x
	  +  radius[0]*radius[2]*Sigma(c);
    
    Qpole(q).yz += Qpole(c).yz
	  +  radius[1]*Dpole(c).z
	  +  radius[2]*Dpole(c).y
	  +  radius[1]*radius[2]*Sigma(c);
    
    
    B0(q) += B0(c);
    
    
    rad2 = Dot(radius, radius);    
    newbmax = sqrtf_fast(rad2) + Bmax(c);
    if( newbmax > Bmax(q) )
	Bmax(q) = newbmax;

    B2(q) += B2(c) +rad2*B0(c);

    if( q->maxsize < c->maxsize )
        q->maxsize = c->maxsize;
    
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
	if (fabs(dx) < xacc) return rtn;
    }
    Error("Maximum number of iterations exceeded in RTNEWT\n");
}

#undef JMAX
