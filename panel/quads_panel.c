

/* Multipole expansion of singular kernel (used for a group of panels)    */


#include <stddef.h>
#include "protos.h"
#include "physics_panel.h"
#include "vop.h"
#include "Assert.h"
#include "Msgs.h"
#include "fastflpt.h"

Counter_t CellInt;

void Cinter(body *me, cell *cp)
{    
    float radius[3];
    float x, y, z;
    float dist2, d2inv, d2inv12;
    float cf0, cf1, cf3, cf5;
    float tmpa, XdotD, XdotQdotX;
    float dxphi, dyphi, dzphi;
    struct vect{float x,y,z;} QdotX;
    float bmaxfac, err;

    IncrCounter(&CellInt);

    VVV(radius, = Pos(me), -Pos(cp));

    x=radius[0];
    y=radius[1];
    z=radius[2];

		
    dist2=x*x+y*y+z*z;

    d2inv12=recipsqrtf(dist2);
    d2inv=d2inv12*d2inv12;


    cf0=    d2inv12;		/* 1/d */
    cf1=    d2inv*cf0;		/* 1/d^3 */
    cf3=3.F*d2inv*cf1;		/* 3/d^5 */
    cf5=5.F*d2inv*cf3;		/* 15/d^7 */

    bmaxfac = recipf(1.F - cf0*Bmax(cp));
    err = (1.F/3.F)*cf3*bmaxfac*bmaxfac*(4.F*B3(cp) - 3.F*B4(cp)*cf0);

    XdotD = x*Dpole(cp).x + y*Dpole(cp).y + z*Dpole(cp).z;

    QdotX.x = x*Qpole(cp).xx + y*Qpole(cp).xy + z*Qpole(cp).xz;
    QdotX.y = x*Qpole(cp).xy + y*Qpole(cp).yy + z*Qpole(cp).yz;
    QdotX.z = x*Qpole(cp).xz + y*Qpole(cp).yz + z*Qpole(cp).zz;

    XdotQdotX = 0.5F*(x*QdotX.x + y*QdotX.y + z*QdotX.z);

	
    Phi(me) += cf0*Sigma(cp) + cf1*XdotD + cf3*XdotQdotX;

    tmpa = cf1*Sigma(cp) + cf3*XdotD + cf5*XdotQdotX;
    dxphi = -x*tmpa + (cf1*Dpole(cp).x + cf3*QdotX.x);
    dyphi = -y*tmpa + (cf1*Dpole(cp).y + cf3*QdotX.y);
    dzphi = -z*tmpa + (cf1*Dpole(cp).z + cf3*QdotX.z);
	

    Vel(me)[0] += dxphi; 
    Vel(me)[1] += dyphi;
    Vel(me)[2] += dzphi;

    Errsum(me) += err;
    Errsum2(me) += err*err;
}    
