/* omegat.c */

/* May 95: this function computes the contribution of neighbor particle ``bp"
           to the particle representation of the vorticity field at
           particle ``me". 
   Revised Nov 27 1995 to make it faster.
*/

#include <math.h>
#include "protos.h"
#include "Msgs.h"
#include "physics_vrtx.h"
#include "vop.h"
#include "fastflpt.h"
#include "timers.h"
#include "neigh.h"

#define SQRT2OPI 7.9788456080287e-01F

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


void OmegatBody(body *bp, body *me, float eps2inv12){

    float radius[3];
    float eps2inv, eps2inv32;
    float dist2, r2, r2o2;
    float z, xi;
    int j;


    eps2inv=eps2inv12*eps2inv12;
    eps2inv32=eps2inv*eps2inv12;
    
    VVV(radius, = Pos(me), - Pos(bp));
    
    dist2 = Dot(radius, radius);

/* The tree-search algorithm may give some bodies that are beyond
   kernel_cutoff.  Is it better to throw them away */
    if( dist2 > kc2 )
	return;

    Msgf(("me->id=%d, omegat dist = %g\n", me->ident, sqrt(dist2)));
    r2 = eps2inv*dist2;
    r2o2=.5F*r2;

    IncrCounter(&NfindAcceptsCnt);

/* Full glory formula for exp(-r2o2) */

    j=(int)r2;
    z=r2o2-.5F*(float)j;

    xi=aexp[j]*(1.-z*(BB1-z*(BB2-z*(BB3-z*(BB4-z*(BB5-z*(BB6-z*BB7)))))));

    xi*=eps2inv32*SQRT2OPI;


    VV(Omegat(me), += xi*Strength(bp));
	
}
