
/* Monopole interaction for the panel (evaluated in local coordinate system): 
   full formula if too close,
   otherwise use the panel multipole expansion centered at panel centroid */



#include <stddef.h>
#include <math.h>
#include "protos.h"
#include "physics_panel.h"
#include "Assert.h"
#include "Msgs.h"
#include "fastflpt.h"
#include "vop.h"

Counter_t BodyFullCnt;
Counter_t BodyQuadCnt;

void Binter(body *me, body *bp){


    float radius[3];

    float xl, yl, zl, zzl;
    float tmpa, Ixxlpp, Iyylpp, XdotQdotX;

    struct vect{
    float xl;
    float yl;
    }QdotX;


    float xmx1l, ymy1l, xmx2l, ymy2l, xmx3l;
    float r1, r2, r3;
    float q12, q23, q31;
    float a12_1, a12_2, a23_2, a23_3, a31_3, a31_1;
    float r12, r23, r31;
    float azl;
    float anum12, anum23, anum31;
    float dnom12, dnom23, dnom31;
    float aj12, aj23, aj31;
    float deltheta; 

    float dist2, d2inv12, d2inv;
    float cf0, cf1, cf3, cf5;

    float phi, dxlphi,   dylphi,   dzlphi;
    float      dxphi,    dyphi,    dzphi;

    VVV(radius, = Pos(me), - Pos(bp));
    
    xl=Dot(Ex(bp), radius);
    yl=Dot(Ey(bp), radius);
    zl=Dot(Ez(bp), radius);
    
    zzl=zl*zl;
    
    dist2=xl*xl+yl*yl+zzl;
    
    /* Full formula if we are too close to the panel:
       This has to be done in the local panel coordinate system.
       For computational savings, we have stored (x1l, y1l), (x2l, y2l)
       and x3l: (Recall that y3l=y1l and hence that s31=0.)             */
    
    if( dist2 == 0.F ){
	
	/* self interaction must be handled as a special case!! 
	   */
	
        deltheta=6.283185307179586F;
	
	dzlphi=-deltheta;   /* that's for z=0+  */
	/* dzlphi=deltheta;  that's for z=0- */
	
        phi=bp->sphi;
	
        dxlphi=bp->dxlsphi;
        dylphi=bp->dylsphi;
	
    }
    
    else {
	float boverd, distinv, err;
	int fullformula;

	distinv = recipsqrtf(dist2);
	boverd = Size(bp)*distinv;
	if( boverd >= 1.0 ){
	    fullformula = 1;
	}else{
	    float dist2inv, dist3inv, dist4inv;
	    float factor, B2, B3, B4, B0;
	    float absSigma;

	    absSigma = fabs(Sigma(bp));
	    if( absSigma != 0.F ){
		dist2inv = distinv*distinv;
		dist3inv = distinv*dist2inv;
		dist4inv = dist2inv*dist2inv;
		factor = recipf(1.F - boverd);
	
		B2 = absSigma*(Ixxlp(bp) + Iyylp(bp));
		B3 = B2*Size(bp);
		B0 = Ip(bp)*absSigma;
		B4 = B2*B2/B0;
		err = factor*factor * dist2inv * 
		    (4.F*B3*dist3inv - 3.F*B4*dist4inv);
		fullformula = (err > errtol);
	    }else{
		/* In fact, we could just return! */
		err = 0.F;
		fullformula = 0;
	    }
	}

	if( fullformula ){
	
	IncrCounter(&BodyFullCnt);

        xmx1l=xl-(bp->x1l);
        ymy1l=yl-(bp->y1l);
	
        xmx2l=xl-(bp->x2l);
        ymy2l=yl-(bp->y2l);
	
        xmx3l=xl-(bp->x3l);
	
	
        r1=sqrtf_fast(xmx1l*xmx1l+ymy1l*ymy1l+zzl);
        r2=sqrtf_fast(xmx2l*xmx2l+ymy2l*ymy2l+zzl);
        r3=sqrtf_fast(xmx3l*xmx3l+ymy1l*ymy1l+zzl);
	
        q12=log((r1+r2+(bp->d12))/(r1+r2-(bp->d12)));
        q23=log((r2+r3+(bp->d23))/(r2+r3-(bp->d23)));
        q31=log((r3+r1+(bp->d31))/(r3+r1-(bp->d31)));
	
        a12_1=-(xmx1l*(bp->c12)+ymy1l*(bp->s12));
        a12_2=-(xmx2l*(bp->c12)+ymy2l*(bp->s12));
	
        a23_2=-(xmx2l*(bp->c23)+ymy2l*(bp->s23));
        a23_3=-(xmx3l*(bp->c23)+ymy1l*(bp->s23));
	
        a31_3=-(xmx3l*(bp->c31));
        a31_1=-(xmx1l*(bp->c31));
	
        r12=xmx1l*(bp->s12)-ymy1l*(bp->c12);  /* signed distance to segment 1-2 */
        r23=xmx2l*(bp->s23)-ymy2l*(bp->c23);
        r31=               -ymy1l*(bp->c31);
	
        azl=fabs(zl);
	
        anum12=r12*azl*(r1*a12_2-r2*a12_1);
        anum23=r23*azl*(r2*a23_3-r3*a23_2);
        anum31=r31*azl*(r3*a31_1-r1*a31_3);
	
        dnom12=r1*r2*r12*r12+zzl*a12_1*a12_2;
        dnom23=r2*r3*r23*r23+zzl*a23_2*a23_3;
        dnom31=r3*r1*r31*r31+zzl*a31_3*a31_1;
	
        aj12=atan2(anum12,dnom12);
        aj23=atan2(anum23,dnom23);
        aj31=atan2(anum31,dnom31);
	
        
        deltheta=0.F;
        if(r12<0.F && r23<0.F && r31< 0.F) deltheta=6.283185307179586F;
	
	dzlphi=-((aj12+aj23+aj31)+deltheta);
        if(zl<0.F) dzlphi=-dzlphi;
	
        phi=-(r12*q12+r23*q23+r31*q31)+zl*dzlphi;
	
        dxlphi=-((bp->s12)*q12+(bp->s23)*q23              );
        dylphi= ((bp->c12)*q12+(bp->c23)*q23+(bp->c31)*q31);
	
	
    }
    
    
    else{
	
	/* We are far enough from the panel that we can use its multipole
	   expansion about its centroid. (Recall that no dipole term here !!):     */
	
	IncrCounter(&BodyQuadCnt);

        tmpa = .333333333333F*(Ixxlp(bp)+Iyylp(bp));
        Ixxlpp = Ixxlp(bp)-tmpa;
        Iyylpp = Iyylp(bp)-tmpa;
	
	
	d2inv12=recipsqrtf(dist2);
        d2inv=d2inv12*d2inv12;
	
	
        cf0=   d2inv12;
        cf1=   d2inv*cf0;
        cf3=3.F*d2inv*cf1;
        cf5=5.F*d2inv*cf3;
	
        QdotX.xl = xl*Ixxlpp    + yl*Ixylp(bp);
        QdotX.yl = xl*Ixylp(bp) + yl*Iyylpp   ;
	
        XdotQdotX = xl*QdotX.xl + yl*QdotX.yl;
	
        phi = cf0*Ip(bp) + cf3*XdotQdotX;
	
	
        tmpa = cf1*Ip(bp) + cf5*XdotQdotX;
	
        dxlphi = -xl*tmpa + cf3*QdotX.xl; 
        dylphi = -yl*tmpa + cf3*QdotX.yl; 
        dzlphi = -zl*tmpa               ; 

	me->errsum += err;
	me->errsum2 += err*err;
    }
    }
    /* derivatives induced by a panel of unit scalar strength, but in absolute
       coordinate system: */
    
    dxphi=Ex(bp)[0]*dxlphi+Ey(bp)[0]*dylphi+Ez(bp)[0]*dzlphi; 
    dyphi=Ex(bp)[1]*dxlphi+Ey(bp)[1]*dylphi+Ez(bp)[1]*dzlphi; 
    dzphi=Ex(bp)[2]*dxlphi+Ey(bp)[2]*dylphi+Ez(bp)[2]*dzlphi; 
    
    /* Panel contribution to velocity field at observation point: */
    /* grad(phi) */
    
    Vel(me)[0] += dxphi*Sigma(bp); 
    Vel(me)[1] += dyphi*Sigma(bp);
    Vel(me)[2] += dzphi*Sigma(bp);
}
