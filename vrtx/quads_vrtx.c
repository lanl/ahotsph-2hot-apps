/* Multipole expansion of singular kernel  */

#include <stddef.h>
#include "physics_vrtx.h"
#include "protos.h"
#include "vop.h"
#include "Assert.h"
#include "Msgs.h"
#include "fastflpt.h"

Counter_t CellCnt;

void InteractCell(cellptr cp, bodyptr me)
{    
    float radius[3];
    float x, y, z, xx, xy, xz, yy, yz, zz;
    float dist2, d2inv, d2inv12;
    float cf0, cf1, cf3, cf5, cf7;
    float tmpa[3], tmpb[3], XdotD[3], XdotQdotX[3];
    float err, bmaxfac;

    struct vect{
	float x[3];
	float y[3];
	float z[3];
    }QdotX, wrk;


    float dxpsi[3], dypsi[3], dzpsi[3];
    float dxxpsi[3], dyypsi[3], dzzpsi[3], dxypsi[3], dxzpsi[3], dyzpsi[3];

    IncrCounter(&CellCnt);
    me->nterms += 2;		/* ??? */

    VVV(radius, = Pos(me), -Pos(cp));
	
    x=radius[0];
    y=radius[1];
    z=radius[2];
	
    xx=x*x;
    xy=x*y;
    xz=x*z;
    yy=y*y;
    yz=y*z;
    zz=z*z;
	
    dist2=xx+yy+zz;

    if( dist2 == 0.F ){
      /* This IS remotely possible, if the cell contains only zero
	 strength bodies it will have an rcrit2 of zero, so we might
	 evaluate the field inside the cell.  We might even evaluate
	 it right on top of the cell's center, which gets us here. */
      assert(Strength(cp)[0] == 0.F && Strength(cp)[1] == 0.F && Strength(cp)[2] == 0.F);
      return;
    }
    /*      d2inv=1.0F/dist2;
	    d2inv12=sqrtf(d2inv);    */

    d2inv12=recipsqrtf(dist2);
    d2inv=d2inv12*d2inv12;


    cf0=    d2inv12;
    cf1=    d2inv*cf0;
    cf3=3.F*d2inv*cf1;
    cf5=5.F*d2inv*cf3;
    cf7=7.F*d2inv*cf5;

    bmaxfac = recipf(1.F - cf0*Bmax(cp));
    err = (1.F/3.F)*cf3*bmaxfac*bmaxfac*(4.F*B3(cp) - 3.F*B4(cp)*cf0);
    Errsum(me) += err;
    Errsum2(me) += err*err;
	
    VVVV(XdotD, = x*Dpole(cp).x, + y*Dpole(cp).y, + z*Dpole(cp).z);

    VVVV(QdotX.x, = x*Qpole(cp).xx, + y*Qpole(cp).xy, + z*Qpole(cp).xz);
    VVVV(QdotX.y, = x*Qpole(cp).xy, + y*Qpole(cp).yy, + z*Qpole(cp).yz);
    VVVV(QdotX.z, = x*Qpole(cp).xz, + y*Qpole(cp).yz, + z*Qpole(cp).zz);

    VVVV(XdotQdotX, = x*QdotX.x, + y*QdotX.y, + z*QdotX.z);
    VS(XdotQdotX, *= 0.5F);


	
    VVVV(Psi(me),  += cf0*Strength(cp), + cf1*XdotD, + cf3*XdotQdotX); 


    VVVV(tmpa, = cf1*Strength(cp), + cf3*XdotD, + cf5*XdotQdotX);

    /*      dxpsi[0] = -x*tmpa[0] + (cf1*Dpole(cp).x[0]+cf3*QdotX.x[0]); */	
    dxpsi[1] = -x*tmpa[1] + (cf1*Dpole(cp).x[1]+cf3*QdotX.x[1]);	
    dxpsi[2] = -x*tmpa[2] + (cf1*Dpole(cp).x[2]+cf3*QdotX.x[2]);
	
    dypsi[0] = -y*tmpa[0] + (cf1*Dpole(cp).y[0]+cf3*QdotX.y[0]);
    /*      dypsi[1] = -y*tmpa[1] + (cf1*Dpole(cp).y[1]+cf3*QdotX.y[1]); */
    dypsi[2] = -y*tmpa[2] + (cf1*Dpole(cp).y[2]+cf3*QdotX.y[2]);	
	
    dzpsi[0] = -z*tmpa[0] + (cf1*Dpole(cp).z[0]+cf3*QdotX.z[0]);	
    dzpsi[1] = -z*tmpa[1] + (cf1*Dpole(cp).z[1]+cf3*QdotX.z[1]);
    /*      dzpsi[2] = -z*tmpa[2] + (cf1*Dpole(cp).z[2]+cf3*QdotX.z[2]); */

	
	
    VVVV(tmpb, = cf3*Strength(cp), + cf5*XdotD, + cf7*XdotQdotX);
        
    /* store what will be used more than once: */
    VVV(wrk.x, = cf3*Dpole(cp).x, + cf5*QdotX.x);
    VVV(wrk.y, = cf3*Dpole(cp).y, + cf5*QdotX.y);
    VVV(wrk.z, = cf3*Dpole(cp).z, + cf5*QdotX.z);

    /*	dxxpsi[0] = xx*tmpb[0]-2.0F*x*wrk.x[0]+cf3*Qpole(cp).xx[0]
	-tmpa[0]; */	
    dxxpsi[1] = xx*tmpb[1]-2.0F*x*wrk.x[1]+cf3*Qpole(cp).xx[1]
	-tmpa[1];	
    dxxpsi[2] = xx*tmpb[2]-2.0F*x*wrk.x[2]+cf3*Qpole(cp).xx[2]
	-tmpa[2];	
	
    dyypsi[0] = yy*tmpb[0]-2.0F*y*wrk.y[0]+cf3*Qpole(cp).yy[0]
	-tmpa[0];	
    /*	dyypsi[1] = yy*tmpb[1]-2.0F*y*wrk.y[1]+cf3*Qpole(cp).yy[1]
	-tmpa[1]; */	
    dyypsi[2] = yy*tmpb[2]-2.0F*y*wrk.y[2]+cf3*Qpole(cp).yy[2]
	-tmpa[2];	
	
    dzzpsi[0] = zz*tmpb[0]-2.0F*z*wrk.z[0]+cf3*Qpole(cp).zz[0]
	-tmpa[0];	
    dzzpsi[1] = zz*tmpb[1]-2.0F*z*wrk.z[1]+cf3*Qpole(cp).zz[1]
	-tmpa[1];	
    /*	dzzpsi[2] = zz*tmpb[2]-2.0F*z*wrk.z[2]+cf3*Qpole(cp).zz[2]
	-tmpa[2]; */	
	
 
    VVVVV(dxypsi, = xy*tmpb, -x*wrk.y, -y*wrk.x, +cf3*Qpole(cp).xy); 
    VVVVV(dxzpsi, = xz*tmpb, -x*wrk.z, -z*wrk.x, +cf3*Qpole(cp).xz); 
    VVVVV(dyzpsi, = yz*tmpb, -y*wrk.z, -z*wrk.y, +cf3*Qpole(cp).yz); 
	
	

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
	
}    

