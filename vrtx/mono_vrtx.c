/* Monopole interaction for smooth kernel: Gaussian smoothing
   exp(x) computed using shift of x together with a polynomial fit, 
   Revised on Nov 27 1995...
   erf(x) computed using a polynomial fit together with exp(x).
   */

#include <math.h>
#include "Msgs.h"
#include "physics_vrtx.h"
#include "vop.h"
#include "fastflpt.h"

Counter_t BodyCnt;
Counter_t TaylorKernelCnt;
Counter_t FullKernelCnt;

#define OO1  1.0000000000000e+00F
#define OO2  5.0000000000000e-01F
#define OO3  3.3333333333333e-01F
#define OO4  2.5000000000000e-01F
#define OO5  2.0000000000000e-01F
#define OO6  1.6666666666667e-01F
#define OO7  1.4285714285714e-01F
#define OO9  1.1111111111111e-01F
#define OO11 9.0909090909091e-02F
#define OO13 7.6923076923077e-02F
#define OO15 6.6666666666667e-02F
#define OO17 5.8823529411765e-02F

#define SQRT2OPI 7.9788456080287e-01F

#define AP 2.3164188e-01F
#define A1 2.54829592e-01F
#define A2 2.84496736e-01F
#define A3 1.421413741F
#define A4 1.453152027F
#define A5 1.061405429F

#define BB1 9.999999995e-01F
#define BB2 4.999999206e-01F
#define BB3 1.666653019e-01F
#define BB4 4.16573475e-02F
#define BB5 8.3013598e-03F
#define BB6 1.3298820e-03F
#define BB7 1.413161e-04F


/* exp(-0.0), exp(-0.5), exp(-1.0), exp(-1.5), ..., exp(-18.0) are stored
in aexp[j] for j=0 to j=36.
These will be used to find exp(-r2o2), with r2o2 <= 18, i.e., r2 <= 36,
i.e., dist2 <= 36 *eps2, i.e., dist <= 6 * eps. Hence kernel_cutoff should
always be less or equal to 6.0 !! */

static float aexp[37]={ 
    1.000000e+00,
    6.065307e-01,
    3.678794e-01,
    2.231302e-01,
    1.353353e-01,
    8.208500e-02,
    4.978707e-02,
    3.019738e-02,
    1.831564e-02,
    1.110900e-02,
    6.737947e-03,
    4.086771e-03,
    2.478752e-03,
    1.503439e-03,
    9.118820e-04,
    5.530844e-04,
    3.354626e-04,
    2.034684e-04,
    1.234098e-04,
    7.485183e-05,
    4.539993e-05,
    2.753645e-05,
    1.670170e-05,
    1.013009e-05,
    6.144212e-06,
    3.726653e-06,
    2.260329e-06,
    1.370959e-06,
    8.315287e-07,
    5.043477e-07,
    3.059023e-07,
    1.855391e-07,
    1.125352e-07,
    6.825603e-08,
    4.139938e-08,
    2.510999e-08,
    1.522998e-08 };

void InteractBody(body *bp, bodyptr me, 
		  float eps2inv12, float nu)
{

    float radius[3];
    float eps2inv, eps2inv32, eps2inv52;
    float dist2, r2, r2o2, r2inv, r2inv12;
    float xi, cf0, cf1, cf3;
    float tmp, z;
    int j;
    float dxpsi[3], dypsi[3], dzpsi[3];
    float dxxpsi[3], dyypsi[3], dzzpsi[3], dxypsi[3], dxzpsi[3], 
    dyzpsi[3];

    IncrCounter(&BodyCnt);
    me->nterms +=1;

    eps2inv=eps2inv12*eps2inv12;

    eps2inv32=eps2inv*eps2inv12;
    eps2inv52=eps2inv*eps2inv32;    
    
    VVV(radius, = Pos(me), - Pos(bp));
    
    dist2 = Dot(radius, radius);

    r2 = eps2inv*dist2;
    r2o2=.5F*r2;


    if(r2o2 <= 0.25F) {	/* Taylor series and no need for a square root */
	IncrCounter(&TaylorKernelCnt);
        me->nterms +=1; /* ???? */

	cf0=eps2inv12*SQRT2OPI*(OO1-r2o2*(OO3-OO2*r2o2*(OO5-OO3*r2o2*
			(OO7-OO4*r2o2*(OO9-OO5*r2o2*(OO11-OO6*r2o2*OO13))))));
	cf1=eps2inv32*SQRT2OPI*(OO3-r2o2*(OO5-OO2*r2o2*(OO7-OO3*r2o2*
			(OO9-OO4*r2o2*(OO11-OO5*r2o2*(OO13-OO6*r2o2*OO15))))));
	cf3=eps2inv52*SQRT2OPI*(OO5-r2o2*(OO7-OO2*r2o2*(OO9-OO3*r2o2*
			(OO11-OO4*r2o2*(OO13-OO5*r2o2*(OO15-OO6*r2o2*OO17))))));

	xi =eps2inv32*SQRT2OPI*(OO1-r2o2*(OO1-OO2*r2o2*(OO1-OO3*r2o2*
			(OO1-OO4*r2o2*(OO1-OO5*r2o2*(OO1-OO6*r2o2*OO1))))));
	/* xi=3.F*cf1-dist2*cf3;  this is true as well! */
        /* xi=eps2inv*(cf0-dist2*cf1);  and so is this! */


    }else{  /* need for a square root for either of the two following cases */
	/*    r2inv=1./r2;
	      r2inv12=sqrtf(r2inv);  */

	r2inv12=recipsqrtf(r2);
	r2inv=r2inv12*r2inv12;

	if(dist2 <= kc2){	/* full glory formulas */

	    IncrCounter(&FullKernelCnt);
            me->nterms +=1; /* ???? */

	    /* first evaluate tmp=exp(-r2o2) 
                0 <= z < 0.5 Hence can use Abramowitz and Stegun
	        pp.71,  art. 4.2.45 which requires 0<= z <= ln2  */

            j=(int)r2;
            z=r2o2-.5F*(float)j;
 
	    tmp=aexp[j]*(1.-z*(BB1-z*(BB2-z*(BB3-z*(BB4-z*(BB5-z*(BB6-z*BB7)))))));
	
	    xi=SQRT2OPI*tmp;

	 /* then evaluate tmp=erf(r/sqrt(2)) Abramowitz and Stegun pp. 299, art.7.1.26  */


	    z=r2inv12/(r2inv12+AP);

	    tmp=1.-tmp*z*(A1-z*(A2-z*(A3-z*(A4-z*A5))));
	
	    cf0=r2inv12*tmp;
	    cf1=r2inv  *(   cf0-xi);
	    cf3=r2inv  *(3.*cf1-xi);
	
	    xi =eps2inv32*xi;
	    cf0=eps2inv12*cf0;
	    cf1=eps2inv32*cf1;
	    cf3=eps2inv52*cf3;

     
	}else{		/* singular kernel */
        
	    xi=(float)0.;
	    cf0=eps2inv12*r2inv12;
	    cf1=eps2inv*r2inv*   cf0;
	    cf3=eps2inv*r2inv*3.F*cf1;

	}

    }

	
    VV(Psi(me), += cf0*Strength(bp));
	
	
    tmp = -cf1*radius[0];
    dxpsi[1] = tmp*Strength(bp)[1];
    dxpsi[2] = tmp*Strength(bp)[2];
	
    tmp = -cf1*radius[1];
    dypsi[0] = tmp*Strength(bp)[0];
    dypsi[2] = tmp*Strength(bp)[2];
	
    tmp = -cf1*radius[2];
    dzpsi[0] = tmp*Strength(bp)[0];
    dzpsi[1] = tmp*Strength(bp)[1];
	

	
    tmp = cf3*radius[0]*radius[0]-cf1;
    dxxpsi[1] = tmp*Strength(bp)[1]; 
    dxxpsi[2] = tmp*Strength(bp)[2]; 
	
    tmp = cf3*radius[1]*radius[1]-cf1;
    dyypsi[0] = tmp*Strength(bp)[0];
    dyypsi[2] = tmp*Strength(bp)[2];
	
    tmp = cf3*radius[2]*radius[2]-cf1;
    dzzpsi[0] = tmp*Strength(bp)[0];   
    dzzpsi[1] = tmp*Strength(bp)[1];   

	
    tmp = cf3*radius[0]*radius[1];
    VV(dxypsi, = tmp*Strength(bp)); 
	
    tmp = cf3*radius[1]*radius[2];
    VV(dyzpsi, = tmp*Strength(bp));

    tmp = cf3*radius[0]*radius[2];
    VV(dxzpsi, = tmp*Strength(bp));
	

	
    Vel(me)[0] += dypsi[2]-dzpsi[1]; 
    Vel(me)[1] += dzpsi[0]-dxpsi[2];
    Vel(me)[2] += dxpsi[1]-dypsi[0];
	
	
    Gradvel(me)[0][0] += dxypsi[2]-dxzpsi[1]; 
    Gradvel(me)[0][1] += dyypsi[2]-dyzpsi[1]; 
    Gradvel(me)[0][2] += dyzpsi[2]-dzzpsi[1]; 
	
    Gradvel(me)[1][0] += dxzpsi[0]-dxxpsi[2];
    Gradvel(me)[1][1] += dyzpsi[0]-dxypsi[2];
    Gradvel(me)[1][2] += dzzpsi[0]-dxzpsi[2];
	
    Gradvel(me)[2][0] += dxxpsi[1]-dxypsi[0];
    Gradvel(me)[2][1] += dxypsi[1]-dyypsi[0];
    Gradvel(me)[2][2] += dxzpsi[1]-dyzpsi[0];



    /* viscous interaction */
	   
    tmp=(2.F*nu*eps2inv)*xi;

    Dstr(me)[0] += tmp*( Vol(me)*Strength(bp)[0]
			-Vol(bp)*Strength(me)[0]);
    Dstr(me)[1] += tmp*( Vol(me)*Strength(bp)[1]
			-Vol(bp)*Strength(me)[1]);
    Dstr(me)[2] += tmp*( Vol(me)*Strength(bp)[2]
			-Vol(bp)*Strength(me)[2]);

}
