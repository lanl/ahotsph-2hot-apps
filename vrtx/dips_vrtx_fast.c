/* Multipole expansion of singular kernel  */

#include <stddef.h>
#include "physics_vrtx.h"
#include "protos.h"
#include "vop.h"
#include "Assert.h"
#include "Msgs.h"
#include "fastflpt.h"

Counter_t CellCnt;

/* 218 flops */
/* 38 loads, 17 stores */
void InteractCell(cellptr cp, bodyptr me)
{    
    float radius[3];
    float x, y, z, xx, xy, xz, yy, yz, zz;
    float dist2, d2inv, d2inv12;
    float cf0, cf1, cf3, cf5, cf7;
    float tmpa[3], tmpb[3], XdotD[3];
    float err, bmaxfac;

    struct vect{
    float x[3];
    float y[3];
    float z[3];
    } wrk;


    float dxpsi[3], dypsi[3], dzpsi[3];
    float dxxpsi[3], dyypsi[3], dzzpsi[3], dxypsi[3], dxzpsi[3], dyzpsi[3];

    IncrCounter(&CellCnt);

    /* 3 flops, 6 loads */
    VVV(radius, = Pos(me), -Pos(cp));
	
    x=radius[0];
    y=radius[1];
    z=radius[2];
	
    /* 6 flops */
    xx=x*x;
    xy=x*y;
    xz=x*z;
    yy=y*y;
    yz=y*z;
    zz=z*z;

    /* 11 flops */
    dist2=xx+yy+zz;
    d2inv12=recipsqrtf(dist2);
    d2inv=d2inv12*d2inv12;

    /* 7 flops */
    cf0=    d2inv12;
    cf1=    d2inv*cf0;
    cf3=3.F*d2inv*cf1;
    cf5=5.F*d2inv*cf3;
    cf7=7.F*d2inv*cf5;

    /* 17 flops, 5 loads, 2 stores */
    bmaxfac = recipf(1.F - cf0*Bmax(cp));
    err = d2inv*d2inv*bmaxfac*bmaxfac*(3.F*B2(cp) - 2.F*B3(cp)*cf0);
    Errsum(me) += err;
    Errsum2(me) += err*err;

    /* 15 flops, 9 loads */
    VVVV(XdotD, = x*Dpole(cp).x, + y*Dpole(cp).y, + z*Dpole(cp).z);

    /* 15 flops, 3 loads, 3 dloads, 3 dstores */
    VVV(Psi(me),  += cf0*Strength(cp), + cf1*XdotD); 

    /* 9 flops */
    VVV(tmpa, = cf1*Strength(cp), + cf3*XdotD);

    /* 18 flops */
    dxpsi[1] = -x*tmpa[1] + (cf1*Dpole(cp).x[1]);	
    dxpsi[2] = -x*tmpa[2] + (cf1*Dpole(cp).x[2]);
    dypsi[0] = -y*tmpa[0] + (cf1*Dpole(cp).y[0]);
    dypsi[2] = -y*tmpa[2] + (cf1*Dpole(cp).y[2]);	
    dzpsi[0] = -z*tmpa[0] + (cf1*Dpole(cp).z[0]);	
    dzpsi[1] = -z*tmpa[1] + (cf1*Dpole(cp).z[1]);
	
    /* 9 flops */
    VVV(tmpb, = cf3*Strength(cp), + cf5*XdotD);
    
    /* store what will be used more than once: */
    /* 9 flops */
    VV(wrk.x, = cf3*Dpole(cp).x);
    VV(wrk.y, = cf3*Dpole(cp).y);
    VV(wrk.z, = cf3*Dpole(cp).z);

    /* 30 flops */
    dxxpsi[1] = xx*tmpb[1]-2.0F*x*wrk.x[1]-tmpa[1];		      
    dxxpsi[2] = xx*tmpb[2]-2.0F*x*wrk.x[2]-tmpa[2];		      
    dyypsi[0] = yy*tmpb[0]-2.0F*y*wrk.y[0]-tmpa[0];		      
    dyypsi[2] = yy*tmpb[2]-2.0F*y*wrk.y[2]-tmpa[2];		      
    dzzpsi[0] = zz*tmpb[0]-2.0F*z*wrk.z[0]-tmpa[0];		      
    dzzpsi[1] = zz*tmpb[1]-2.0F*z*wrk.z[1]-tmpa[1];		      
 
    /* 45 flops */
    VVVV(dxypsi, = xy*tmpb, -x*wrk.y, -y*wrk.x); 
    VVVV(dxzpsi, = xz*tmpb, -x*wrk.z, -z*wrk.x); 
    VVVV(dyzpsi, = yz*tmpb, -y*wrk.z, -z*wrk.y); 
	
    /* 6 flops, 3 dloads, 3 dstores */
    Vel(me)[0] += dypsi[2]-dzpsi[1]; 
    Vel(me)[1] += dzpsi[0]-dxpsi[2];
    Vel(me)[2] += dxpsi[1]-dypsi[0];
	
    /* 6 flops, 3 dloads, 3 dstores */
    Gradvel(me)[0][0] += dxypsi[2]-dxzpsi[1]; 
    Gradvel(me)[0][1] += dyypsi[2]-dyzpsi[1]; 
    Gradvel(me)[0][2] += dyzpsi[2]-dzzpsi[1]; 
	
    /* 6 flops, 3 dloads, 3 dstores */
    Gradvel(me)[1][0] += dxzpsi[0]-dxxpsi[2];
    Gradvel(me)[1][1] += dyzpsi[0]-dxypsi[2];
    Gradvel(me)[1][2] += dzzpsi[0]-dxzpsi[2];
	
    /* 6 flops, 3 dloads, 3 dstores */
    Gradvel(me)[2][0] += dxxpsi[1]-dxypsi[0];
    Gradvel(me)[2][1] += dxypsi[1]-dyypsi[0];
    Gradvel(me)[2][2] += dxzpsi[1]-dyzpsi[0];
}    
    
