
/* This function ``prepares" the panels, i.e., precomputes all the stuff one might
   need many times and does not want to compute more than once. */

/* It also initializes the panels strength vector for the iterative solver. */ 

#include <stdlib.h>
#include <math.h>
#include "physics_panel.h"
#include "vop.h"

double drand48(void);


void PreparePanel(bodyptr btab, int n, 
		  void (*Uexternal)(body *bp))
{

      float oo2, oo3, oo4;
      float xp1, yp1, zp1;
      float xp2, yp2, zp2;
      float xp3, yp3, zp3;
      float xpc, ypc, zpc;
      float xp13, yp13, zp13;
      float xp21, yp21, zp21;
      float anorminv;
      float exx, exy, exz;
      float eyx, eyy, eyz;
      float ezx, ezy, ezz;
      float xpt1, ypt1, zpt1;
      float xpt2, ypt2, zpt2;
      float xpt3, ypt3, zpt3;
      float x1, y1, x2, y2, x3;
      float d12, d23, d31;
      float c12, c23, c31;
      float s12, s23;
      float r1, r2, r3;
      float q12, q23, q31;
      float r12, r23, r31;
      float deltheta; 
      float phi, dxphi, dyphi, dzphi;
      float xlc, ylc, xl2, yl2, xl3, xl23;
      float aI, aIxx, aIyy, aIxy;
      float aIxxp, aIyyp, aIzzp, aIxyp, aIxzp, aIyzp;
      int i;
      bodyptr bp;
      float size;



      oo2=(float)0.5;
      oo3=(float)0.333333333333333333;
      oo4=(float)0.25;


      for(i=0; i<n; i++)
      {

      bp=btab+i;


/* Panel corner coordinates in absolute coordinate system: */

      xp1=Pos1(bp)[0];
      yp1=Pos1(bp)[1];
      zp1=Pos1(bp)[2];

      xp2=Pos2(bp)[0];
      yp2=Pos2(bp)[1];
      zp2=Pos2(bp)[2];

      xp3=Pos3(bp)[0];
      yp3=Pos3(bp)[1];
      zp3=Pos3(bp)[2];

/* panel centroid: */

      xpc=oo3*(xp1+xp2+xp3);
      ypc=oo3*(yp1+yp2+yp3);
      zpc=oo3*(zp1+zp2+zp3);


/* segment 1-3: */

      xp13=xp3-xp1;
      yp13=yp3-yp1;
      zp13=zp3-zp1;

/* segment 2-1: */

      xp21=xp1-xp2;
      yp21=yp1-yp2;
      zp21=zp1-zp2;

/* ex= unit vector parallel to 1-3: */

      exx=xp13;
      exy=yp13;
      exz=zp13;
      anorminv=(float)1.0/sqrt(exx*exx+exy*exy+exz*exz);
      exx*=anorminv;
      exy*=anorminv;
      exz*=anorminv;

/* ez= outward unit vector: */

      ezx=yp13*zp21-zp13*yp21;
      ezy=zp13*xp21-xp13*zp21;
      ezz=xp13*yp21-yp13*xp21;
      anorminv=(float)1.0/sqrt(ezx*ezx+ezy*ezy+ezz*ezz);
      ezx*=anorminv;
      ezy*=anorminv;
      ezz*=anorminv;

/* ey = ez X ex: */

      eyx=ezy*exz-ezz*exy;
      eyy=ezz*exx-ezx*exz;
      eyz=ezx*exy-ezy*exx;


/* panel corner points in local coordinate system centered at panel centroid: */

/* first translate: */

      xpt1=xp1-xpc;
      ypt1=yp1-ypc;
      zpt1=zp1-zpc;

      xpt2=xp2-xpc;
      ypt2=yp2-ypc;
      zpt2=zp2-zpc;

      xpt3=xp3-xpc;
      ypt3=yp3-ypc;
      zpt3=zp3-zpc;


/* then rotate: */

      x1=exx*xpt1+exy*ypt1+exz*zpt1;
      y1=eyx*xpt1+eyy*ypt1+eyz*zpt1;

      x2=exx*xpt2+exy*ypt2+exz*zpt2;
      y2=eyx*xpt2+eyy*ypt2+eyz*zpt2;

      x3=exx*xpt3+exy*ypt3+exz*zpt3;    /* y3=y1 */


/* Lengths of the panel sides: */

      d12=sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
      d23=sqrt((x3-x2)*(x3-x2)+(y1-y2)*(y1-y2));
      d31=sqrt((x1-x3)*(x1-x3)                );

/* Sines and Cosines: */

      c12=(x2-x1)/d12;
      s12=(y2-y1)/d12;
      c23=(x3-x2)/d23;
      s23=(y1-y2)/d23;
      c31=(x1-x3)/d31;    /* s31=0 */


/* 
 Self-Interaction potential and its derivatives (In local coordinate system!!)
 for panel of unit scalar strength: 
*/

        r1=sqrt(x1*x1+y1*y1);
        r2=sqrt(x2*x2+y2*y2);
        r3=sqrt(x3*x3+y1*y1);

        q12=log((r1+r2+d12)/(r1+r2-d12));
        q23=log((r2+r3+d23)/(r2+r3-d23));
        q31=log((r3+r1+d31)/(r3+r1-d31));

        r12=-x1*s12+y1*c12;  /* signed distance to segment 1-2 */
        r23=-x2*s23+y2*c23;
        r31=        y1*c31;

        deltheta=6.283185307179586;

	dzphi=deltheta;      /* that's for z=0- */

        phi=-(r12*q12+r23*q23+r31*q31);

        dxphi=-(s12*q12+s23*q23        );
        dyphi= (c12*q12+c23*q23+c31*q31);



/* Panel multipole expansion coefficients with respect to panel centroid, 
   and in local coordinate system centered at panel centroid:
   (note: the formulas were derived in a local coordinate system centered at
   point 1. Hence the shift: */

      xlc=-x1;
      ylc=-y1;

      xl2=x2+xlc;
      yl2=y2+ylc;

      xl3=x3+xlc;

/* also xl1=yl1=yl3=0
   also xlc=(xl2+xl3)/3, ylc=yl2/3  */

      xl23=xl3-xl2;
      
      aI=-yl2*xl3*oo2;

      aIxx=-yl2*( xl2*xl2*xl2*oo4
                 +xl23*(xl3*xl3*oo2-xl23*(2.*xl3*oo3-xl23*oo4)) );
      aIxx-=xlc*xlc*aI;

      aIyy=-yl2*yl2*yl2*xl3*oo3*oo4;
      aIyy-=ylc*ylc*aI;

      aIxy=-yl2*yl2*oo2*( xl2*xl2*oo4
                         +xl23*(xl3*oo3-xl23*oo4) );
      aIxy-=xlc*ylc*aI;
      

/* Panel multipole expansion coefficients with respect to panel centroid,
   but in absolute coordinate system: */

      aIxxp=exx*exx*aIxx+eyx*eyx*aIyy+2.*exx*eyx*aIxy;
      aIyyp=exy*exy*aIxx+eyy*eyy*aIyy+2.*exy*eyy*aIxy;
      aIzzp=exz*exz*aIxx+eyz*eyz*aIyy+2.*exz*eyz*aIxy;
      aIxyp=exx*exy*aIxx+eyx*eyy*aIyy+(exx*eyy+eyx*exy)*aIxy;
      aIxzp=exx*exz*aIxx+eyx*eyz*aIyy+(exx*eyz+eyx*exz)*aIxy;
      aIyzp=exy*exz*aIxx+eyy*eyz*aIyy+(exy*eyz+eyy*exz)*aIxy;




/* Multipole expansion error stuff: */
/* Squared distances from corner points to centroid and their max: */
      
      size= ( r1 > r2 )?  ( ( r1 > r3 )? r1 : r3 ) 
                          : ( ( r2 > r3 )? r2 : r3 ) ;


#if 0				/* not used any more */
      {float aItrace, beta,
            x_new, x, xx, xxx, dist2crit;
      aItrace=aIxx+aIyy;
      beta=aItrace/(aI*size*size);      /* shape factor: 0 < beta < 1 */

/* find the zeroes of: (1-x)^2 rel_tol = 4 beta x^3 - 3 beta^2 x^4, where rel_tol is
   the relative error tolerance acceptable for using the multipole expansion
   formula instead of the full glory formula in mono_panel.c, and x=(size/dist)
   acceptable */
/* first guess: */

      x_new= pow((rel_tol/(4.*beta)), .3333);

/* then N-R iteration: */

      do {
       x=x_new;
       xx=x*x;
       xxx=xx*x;
       x_new=x-( ( beta*xxx*(3.*beta*x-4.)+rel_tol*(1.-2.*x+xx) )/
                 ( 12.*beta*xx*(beta*x-1.)-2.*rel_tol*(1-x)     )  );
      } while( fabs((x_new-x)/x_new) > 1.e-5 );
      
      dist2crit=(size*size)/(x_new*x_new);
/*      bp->shapefact=beta; */
      bp->dist2crit=dist2crit;
   }
#endif

/* Now we can transfer the information to the tree code */

      Pos(bp)[0]=xpc;
      Pos(bp)[1]=ypc;
      Pos(bp)[2]=zpc;

      Ex(bp)[0]=exx; 
      Ex(bp)[1]=exy;    
      Ex(bp)[2]=exz;

      Ey(bp)[0]=eyx; 
      Ey(bp)[1]=eyy;    
      Ey(bp)[2]=eyz;

      Ez(bp)[0]=ezx; 
      Ez(bp)[1]=ezy;    
      Ez(bp)[2]=ezz;

      Uexternal(bp);

      bp->x1l=x1;
      bp->y1l=y1;
      bp->x2l=x2;
      bp->y2l=y2;
      bp->x3l=x3;

      bp->d12=d12;
      bp->d23=d23;
      bp->d31=d31;

      bp->c12=c12;
      bp->c23=c23;
      bp->c31=c31;

      bp->s12=s12;
      bp->s23=s23;

      bp->sphi=phi;
      bp->dxlsphi=dxphi;
      bp->dylsphi=dyphi;

      Ip(bp)=aI;

      Ixxlp(bp)=aIxx;
      Iyylp(bp)=aIyy;
      Ixylp(bp)=aIxy;

      Ixxp(bp)=aIxxp;
      Iyyp(bp)=aIyyp;
      Izzp(bp)=aIzzp;
      Ixyp(bp)=aIxyp;
      Ixzp(bp)=aIxzp;
      Iyzp(bp)=aIyzp;

      Size(bp)=size;

      /* Set sigma and vel self-consistently to zero. */
      /* This way we can start the iteration with an Jacobi Update */
      /* Which gives us an initial guess equal to uext.ez/4pi */
      Sigma(bp) = 0.F;
      VS(Vel(bp), = 0.F);
      }

}


