/* Multipole expansion of singular kernel  */

#include <stddef.h>
#include "physics_vrtx.h"
#include "protos.h"
#include "Assert.h"
#include "Msgs.h"
#include "vop.h"
#include "fastflpt.h"

Counter_t CellCnt;

void Ddot(float *z, float r0, float r1, float r2, 
	  const float *y0, const float *y1, const float *y2);

void Ddot_asm(float *z, float r0, float r1, float r2, 
	  const float *y0, const float *y1, const float *y2);

#if !defined(__INTEL_SSD__) || defined(NO_ASM)
void Ddot(float *z, float r0, float r1, float r2, 
	  const float *y0, const float *y1, const float *y2)
{
    z[0] = r0*y0[0];
    z[1] = r0*y0[1];
    z[2] = r0*y0[2];

    z[0] += r1*y1[0];
    z[1] += r1*y1[1];
    z[2] += r1*y1[2];

    z[0] += r2*y2[0];
    z[1] += r2*y2[1];
    z[2] += r2*y2[2];
}
#else
#define Ddot Ddot_asm
#endif

/* inputs: me.pos (3), cp.pos (3), cp.dpole (9), cp.qpole (18) */
/* outputs: me.vel (3), me.psi (3), me.gradvel (9) */

void InteractCell(cellptr cp, bodyptr me)
{    
    float x, y, z;
    float d2inv, d2inv12;
    float cf0, cf1, cf3, cf5, cf7;
    float tmpa[3], tmpb[3], tmppsi[3], XdotD[3], XdotQdotX[3];

    struct vect{
	float x[3];
	float y[3];
	float z[3];
    }QdotX, wrk;

    float dxpsi[3], dypsi[3], dzpsi[3];
    float dxxpsi[3], dyypsi[3], dzzpsi[3], dxypsi[3], dxzpsi[3], dyzpsi[3];

    IncrCounter(&CellCnt);
    me->nterms += 2;		/* ??? */

#ifdef __INTEL_SSD__
    clear_tregs();
#endif

    /* 3 flops, 6 loads */
    x = Pos(me)[0] - Pos(cp)[0];
    y = Pos(me)[1] - Pos(cp)[1];
    z = Pos(me)[2] - Pos(cp)[2];

    /* 21 flops */
    cf0 = recipsqrtf(x*x+y*y+z*z);
    d2inv = cf0*cf0;
    cf1 = d2inv*cf0;
    cf3 = 3.F*d2inv*cf1;
    cf5 = 5.F*d2inv*cf3;
    cf7 = 7.F*d2inv*cf5;

    /* 15 flops, 9 loads */
    Ddot(XdotD, x, y, z, Dpole(cp).x, Dpole(cp).y, Dpole(cp).z);

    /* 45 flops, 18 loads */
    Ddot(QdotX.x, x, y, z, Qpole(cp).xx, Qpole(cp).xy, Qpole(cp).xz);
    Ddot(QdotX.y, x, y, z, Qpole(cp).xy, Qpole(cp).yy, Qpole(cp).yz);
    Ddot(QdotX.z, x, y, z, Qpole(cp).xz, Qpole(cp).yz, Qpole(cp).zz);

    /* 18 flops */
    Ddot(XdotQdotX, x, y, z, QdotX.x, QdotX.y, QdotX.z);
    VS(XdotQdotX, *= 0.5F);
	
    /* 15 flops */
    Ddot(tmppsi, cf0, cf1, cf3, Strength(cp), XdotD, XdotQdotX);
    /* 3 flops, 3 dloads, 3 dstores */
    VV(Psi(me), += tmppsi);

    /* 15 flops */
    Ddot(tmpa, cf1, cf3, cf5, Strength(cp), XdotD, XdotQdotX);

    /* 15 flops */
    Ddot(tmpb, cf3, cf5, cf7, Strength(cp), XdotD, XdotQdotX);

    /* 30 useful flops (Don't need diagonal components) */
    Ddot(dxpsi, -x, cf1, cf3, tmpa, Dpole(cp).x, QdotX.x);
    Ddot(dypsi, -y, cf1, cf3, tmpa, Dpole(cp).y, QdotX.y);
    Ddot(dzpsi, -z, cf1, cf3, tmpa, Dpole(cp).z, QdotX.z);

    /* 6 flops, 3 dloads, 3 dstores */
    Vel(me)[0] += dypsi[2]-dzpsi[1]; 
    Vel(me)[1] += dzpsi[0]-dxpsi[2];
    Vel(me)[2] += dxpsi[1]-dypsi[0];
	
    /* store what will be used more than once: */
    /* 27 useful flops */
    Ddot(wrk.x, cf3, cf5, 0.0F, Dpole(cp).x, QdotX.x, QdotX.x /* dummy */);
    Ddot(wrk.y, cf3, cf5, 0.0F, Dpole(cp).y, QdotX.y, QdotX.y /* dummy */);
    Ddot(wrk.z, cf3, cf5, 0.0F, Dpole(cp).z, QdotX.z, QdotX.z /* dummy */);

    /* 42 useful flops (Don't need diagonal components) */
    Ddot(dxxpsi, x*x, -2.0F*x, cf3, tmpb, wrk.x, Qpole(cp).xx);
    Ddot(dyypsi, y*y, -2.0F*y, cf3, tmpb, wrk.y, Qpole(cp).yy);
    Ddot(dzzpsi, z*z, -2.0F*z, cf3, tmpb, wrk.z, Qpole(cp).zz);
    VV(dxxpsi, -= tmpa);
    VV(dyypsi, -= tmpa);
    VV(dzzpsi, -= tmpa);

    /* 63 flops */
    Ddot(dxypsi, x*y, -x, -y, tmpb, wrk.y, wrk.x);
    Ddot(dxzpsi, x*z, -x, -z, tmpb, wrk.z, wrk.x);
    Ddot(dyzpsi, y*z, -y, -z, tmpb, wrk.z, wrk.y);
    VV(dxypsi, += cf3*Qpole(cp).xy);
    VV(dxzpsi, += cf3*Qpole(cp).xz);
    VV(dyzpsi, += cf3*Qpole(cp).yz);

	
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

