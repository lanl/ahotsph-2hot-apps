#include <math.h> /* for pow() prototype */
#include "physics_sph.h"
#include "fastflpt.h"
#include "vop.h"

void get_metric(body *btab, int nobj);
void add_gr(body *btab, int nobj);

static float hole_m;
static float ang_mom;
static int do_kerr;
static int did_setup;

void
setup_metric(int kerr_flag, float hole_mass, float kerr_ang_mom)
{
/* Add hole_mass and kerr_ang_mom to setup - Don */
    hole_m = hole_mass; /* Don */
    ang_mom = kerr_ang_mom; /* Don */
    do_kerr = kerr_flag;
    did_setup = 1;
}

void
get_metric(body *btab, int nobj)
{
   body *p;

   float tiny = 1.0e-20;
   float x, y, z;
   float r, r2, r4, ct, st, ct2, st2, c2t, s2t, c2t2, s2t2, cp, sp;
   float cp2, sp2, c2p, s2p, m, a, a2, a3, a4, a5, a6;
   float rho, rho2, rho3, rho4, rho5, rho6, rho7;

   for (p = btab; p < btab+nobj; p++) {

	/* Auxiliary junk */

         x = p->pos[0];
         y = p->pos[1];
         z = p->pos[2];

         if ((x==0.0)&&(y==0.0)) {
            x += tiny;
            y += tiny;
         }

         r       = sqrtf_fast(Dot(p->pos,p->pos));
         r2      =  r*r;
         r4      =  r2*r2;
         ct      =  z/r;
         st      =  sqrtf_fast(x*x+y*y)/r;
         ct2     =  ct*ct;
         st2     =  st*st;
         c2t     =  ct2-st2;
         s2t     =  2.*st*ct;
         c2t2    =  c2t*c2t;
         s2t2    =  s2t*s2t;
         cp      =  x/r/st;
         sp      =  y/r/st;
         cp2     =  cp*cp;
         sp2     =  sp*sp;
         c2p     =  cp2-sp2;
         s2p     =  2.*sp*cp;
	 m = hole_m;		/* Don */
         a = ang_mom;   	/* Don */
         a2      =  a*a;
         a3      =  a*a2;
         a4      =  a*a3;
	 a5      =  a*a4;
	 a6      =  a*a5;
	 rho2    =  0.5*(r2-a2)+0.5*sqrtf_fast(r4+a4+2.*a2*r2*c2t);
         rho     =  sqrtf_fast(rho2);
	 rho3    =  rho*rho2;
         rho4    =  rho*rho3;
         rho5    =  rho*rho4;
         rho6    =  rho*rho5;
         rho7    =  rho*rho6;

	/* contravariant components of the Kerr metric */

         p->gtt=-1+(2*m*rho)/(a2*ct2+rho2);
         p->gxt=(2*m*r*rho*(cp*rho+a*sp)*st)/((a2+rho2)*(a2*ct2+rho2));
         p->gyt=(-2*m*r*rho*(a*cp-rho*sp)*st)/((a2+rho2)*(a2*ct2+rho2));
         p->gzt=(2*m*r*ct)/(a2*ct2+rho2);
         p->gxx=1+(2*m*r2*pow(a2+rho2,-2)*pow(cp*rho+a*sp,2)*st2*rho)/
            (a2*ct2+rho2);
         p->gxy=(-2*m*r2*pow(a2+rho2,-2)*st2*rho*(cp*rho+a*sp)*(a*cp-rho*sp))/
            (a2*ct2+rho2);
         p->gxz=(2*m*ct*r2*(cp*rho+a*sp)*st)/((a2+rho2)*(a2*ct2+rho2));
         p->gyy=1+(2*m*r2*pow(a2+rho2,-2)*pow(a*cp-rho*sp,2)*st2*rho)/
            (a2*ct2+rho2);
         p->gyz=(-2*m*ct*r2*(a*cp-rho*sp)*st)/((a2+rho2)*(a2*ct2+rho2));
         p->gzz=1+(2*m*r2*ct2)/((a2*ct2+rho2)*rho);

	/* covariant components of the Kerr metric */

         p->gutt=(-2*(a2*ct2+rho2)*(2*m*r2+rho3+a2*rho))/(4*m*a2*r2*ct2+4*m*
            (-a2+r2)*rho2+(3+c2t)*a2*rho3-4*m*rho4+2*rho5+2*a4*ct2*rho);
         p->guxt=(4*m*r*rho2*(cp*rho+a*sp)*st)/(4*m*a2*r2*ct2+4*m*(-a2+r2)*
            rho2+(3+c2t)*a2*rho3-4*m*rho4+2*rho5+2*a4*ct2*rho);
         p->guyt=(4*m*r*rho2*(-(a*cp)+rho*sp)*st)/(4*m*a2*r2*ct2+4*m*(-a2+r2)*
            rho2+(3+c2t)*a2*rho3-4*m*rho4+2*rho5+2*a4*ct2*rho);
         p->guzt=(4*m*r*ct*(a2+rho2)*rho)/(4*m*a2*r2*ct2+4*m*(-a2+r2)*rho2+
            (3+c2t)*a2*rho3-4*m*rho4+2*rho5+2*a4*ct2*rho);
         p->guxx=(2*(2*m*a4*r2*ct2+((5+c2t)*a2*rho5)/2.-2*m*rho6+rho7+2*m*a2*
            rho2*(-a2+2*r2*ct2+r2*cp2*st2)+2*m*rho4*(-2*a2+r2*ct2+r2*sp2*st2)+
            a6*ct2*rho+rho3*(a4+2*a4*ct2-2*a*m*r2*st2*s2p)))/((a2+rho2)*
            (4*m*a2*r2*ct2+4*m*(-a2+r2)*rho2+(3+c2t)*a2*rho3-4*m*rho4+2*rho5+
            2*a4*ct2*rho));
         p->guxy=(-2*m*r2*rho2*st2*(-2*a*c2p*rho-a2*s2p+rho2*s2p))/((a2+rho2)*
            (4*m*a2*r2*ct2+4*m*(-a2+r2)*rho2+(3+c2t)*a2*rho3-4*m*rho4+2*rho5+
            2*a4*ct2*rho));
         p->guxz=(-2*m*r2*rho*(cp*rho+a*sp)*s2t)/(4*m*a2*r2*ct2+4*m*(-a2+r2)*
            rho2+(3+c2t)*a2*rho3-4*m*rho4+2*rho5+2*a4*ct2*rho);
         p->guyy=(2*(2*m*a4*r2*ct2+((5+c2t)*a2*rho5)/2.-2*m*rho6+rho7+2*m*
            rho4*(-2*a2+r2*ct2+r2*cp2*st2)+2*m*a2*rho2*(-a2+2*r2*ct2+r2*sp2*
            st2)+a6*ct2*rho+rho3*(a4+2*a4*ct2+2*a*m*r2*st2*s2p)))/((a2+rho2)*
            (4*m*a2*r2*ct2+4*m*(-a2+r2)*rho2+(3+c2t)*a2*rho3-4*m*rho4+2*rho5+
            2*a4*ct2*rho));
         p->guyz=(-4*m*ct*r2*rho*(-(a*cp)+rho*sp)*st)/(4*m*a2*r2*ct2+4*m*(-a2+
            r2)*rho2+(3+c2t)*a2*rho3-4*m*rho4+2*rho5+2*a4*ct2*rho);
         p->guzz=(rho*(2*a4*ct2+(3+c2t)*a2*rho2-4*m*rho3+2*rho4-2*m*(2*a2-r2+
            c2t*r2)*rho))/(4*m*a2*r2*ct2+4*m*(-a2+r2)*rho2+(3+c2t)*a2*rho3-4*
            m*rho4+2*rho5+2*a4*ct2*rho);

   /* Compute the redshifted mass */

         p->gr_mass = p->mass / sqrt(
                      p->gxx * (p->gyy*p->gzz - p->gyz*p->gyz)
		    - p->gxy * (p->gxy*p->gzz - p->gyz*p->gxz)
		    + p->gxz * (p->gxy*p->gyz - p->gyy*p->gxz));

   }
}

void
add_gr(body *btab, int nobj)
{
    body *p;

    float tiny = 1.0e-20;
    float x, y, z, vx, vy, vz;
    float r, r2, r4;
    float ct, st, ct2, st2, ct3, st3, c2t, s2t, c2t2, s2t2, c3t, s3t, s4t;
    float cp, sp, cp2, sp2, c2p, s2p;
    float m, a, a2, a3, a4, a5, a6;
    float rho, rho2, rho3, rho4, rho5, rho6, rho7, rhor, rhot;

   float muk, q, tosc, trelax; /* Karen: For Roche stuff */
    
    float uut, ut, ux, uy, uz;
    float gttx, gtty, gttz;
    float gxxx, gxxy, gxxz;
    float gxyx, gxyy, gxyz;
    float gxzx, gxzy, gxzz;
    float gyyx, gyyy, gyyz;
    float gyzx, gyzy, gyzz;
    float gzzx, gzzy, gzzz;

    float gxtx, gxty, gxtz;
    float gytx, gyty, gytz;
    float gztx, gzty, gztz;

    float curx, cury, curz;

    float gxxr, gxxt, gxxp;
    float gxyr, gxyt, gxyp;
    float gxzr, gxzt, gxzp;
    float gxtr, gxtt, gxtp;
    float gyyr, gyyt, gyyp;
    float gyzr, gyzt, gyzp;
    float gytr, gytt, gytp;
    float gzzr, gzzt, gzzp;
    float gztr, gztt, gztp;
    float gttr, gttt, gttp;

    float drdx, drdy, drdz, dtdx, dtdy, dtdz, dpdx, dpdy, dpdz;
    float dhdr, dhdt, dhdp, deth, hdetx, hdety, hdetz;

    /* begin el loopo grande */

    for (p = btab; p < btab+nobj; p++) {

	/* Auxiliary junk */

         x = p->pos[0];
         y = p->pos[1];
         z = p->pos[2];

         vx = p->vel[0];
         vy = p->vel[1];
         vz = p->vel[2];

         if ((x==0.0)&&(y==0.0)) {
            x += tiny;
            y += tiny;
         }

         r       = sqrtf_fast(Dot(p->pos,p->pos));
         r2      =  r*r;
         r4      =  r2*r2;
         ct      =  z/r;
	 st      =  sqrtf_fast(x*x+y*y)/r;
         ct2     =  ct*ct;
         st2     =  st*st;
	 ct3     =  ct*ct2;
	 st3     =  st*st2;
         c2t     =  ct2-st2;
         s2t     =  2.*st*ct;
         c2t2    =  c2t*c2t;
         s2t2    =  s2t*s2t;
	 c3t     =  ct*c2t-st*s2t;
	 s3t     =  st*c2t+s2t*ct;
	 s4t     =  2.*s2t*c2t;
         cp      =  x/r/st;
         sp      =  y/r/st;
         cp2     =  cp*cp;
         sp2     =  sp*sp;
         c2p     =  cp2-sp2;
         s2p     =  2.*sp*cp;
	 m = hole_m;		/* Don */
         a = ang_mom;  	 	/* Don */
         a2      =  a*a;
         a3      =  a*a2;
         a4      =  a*a3;
	 a5      =  a*a4;
	 a6      =  a*a5;
	 rho2    =  0.5*(r2-a2)+0.5*sqrtf_fast(r4+a4+2.*a2*r2*c2t);
         rho     =  sqrtf_fast(rho2);
	 rho3    =  rho*rho2;
         rho4    =  rho*rho3;
         rho5    =  rho*rho4;
         rho6    =  rho*rho5;
         rho7    =  rho*rho6;
	 rhor = 0.5*r/rho*(1.+(r2+a2*c2t)/sqrtf_fast(r4+a4+2.*a2*r2*c2t));
	 rhot = -0.5*a2*r2*s2t/(rho*sqrtf_fast(r4+a4+2.*a2*r2*c2t));

/* Derivatives with respect to (r,theta,phi) of the inverse metric */

         gttr=(-4*rho*(2*m*r2+rho3+a2*rho)*rhor)/(4*m*a2*r2*ct2+4*m*(-a2+r2)*
            rho2+(3+c2t)*a2*rho3-4*m*rho4+2*rho5+2*a4*ct2*rho)-(2*(a2*ct2+
            rho2)*(4*m*r+a2*rhor+3*rho2*rhor))/(4*m*a2*r2*ct2+4*m*(-a2+r2)*
            rho2+(3+c2t)*a2*rho3-4*m*rho4+2*rho5+2*a4*ct2*rho)+2*(a2*ct2+
            rho2)*pow(4*m*a2*r2*ct2+4*m*(-a2+r2)*rho2+(3+c2t)*a2*rho3-4*m*
            rho4+2*rho5+2*a4*ct2*rho,-2)*(2*m*r2+rho3+a2*rho)*(8*m*r*a2*ct2+
            8*m*r*rho2+2*a4*ct2*rhor+3*(3+c2t)*a2*rho2*rhor-16*m*rho3*rhor+
            10*rho4*rhor+8*m*(-a2+r2)*rho*rhor);
         gxtr=(4*m*rho2*(cp*rho+a*sp)*st)/(4*m*a2*r2*ct2+4*m*(-a2+r2)*rho2+
            (3+c2t)*a2*rho3-4*m*rho4+2*rho5+2*a4*ct2*rho)+(4*m*r*cp*rho2*st*
            rhor)/(4*m*a2*r2*ct2+4*m*(-a2+r2)*rho2+(3+c2t)*a2*rho3-4*m*rho4+
            2*rho5+2*a4*ct2*rho)+(8*m*r*rho*(cp*rho+a*sp)*st*rhor)/(4*m*a2*
            r2*ct2+4*m*(-a2+r2)*rho2+(3+c2t)*a2*rho3-4*m*rho4+2*rho5+2*a4*
            ct2*rho)-4*m*r*rho2*pow(4*m*a2*r2*ct2+4*m*(-a2+r2)*rho2+(3+c2t)*
            a2*rho3-4*m*rho4+2*rho5+2*a4*ct2*rho,-2)*(cp*rho+a*sp)*st*(8*m*r*
            a2*ct2+8*m*r*rho2+2*a4*ct2*rhor+3*(3+c2t)*a2*rho2*rhor-16*m*rho3*
            rhor+10*rho4*rhor+8*m*(-a2+r2)*rho*rhor);
         gytr=(4*m*rho2*(-(a*cp)+rho*sp)*st)/(4*m*a2*r2*ct2+4*m*(-a2+r2)*
            rho2+(3+c2t)*a2*rho3-4*m*rho4+2*rho5+2*a4*ct2*rho)+(4*m*r*rho2*
            sp*st*rhor)/(4*m*a2*r2*ct2+4*m*(-a2+r2)*rho2+(3+c2t)*a2*rho3-4*m*
            rho4+2*rho5+2*a4*ct2*rho)+(8*m*r*rho*(-(a*cp)+rho*sp)*st*rhor)/
            (4*m*a2*r2*ct2+4*m*(-a2+r2)*rho2+(3+c2t)*a2*rho3-4*m*rho4+2*rho5+
            2*a4*ct2*rho)-4*m*r*rho2*pow(4*m*a2*r2*ct2+4*m*(-a2+r2)*rho2+(3+
            c2t)*a2*rho3-4*m*rho4+2*rho5+2*a4*ct2*rho,-2)*(-(a*cp)+rho*sp)*
            st*(8*m*r*a2*ct2+8*m*r*rho2+2*a4*ct2*rhor+3*(3+c2t)*a2*rho2*rhor-
            16*m*rho3*rhor+10*rho4*rhor+8*m*(-a2+r2)*rho*rhor);
         gztr=(4*m*ct*(a2+rho2)*rho)/(4*m*a2*r2*ct2+4*m*(-a2+r2)*rho2+(3+
            c2t)*a2*rho3-4*m*rho4+2*rho5+2*a4*ct2*rho)+(8*m*r*ct*rho2*rhor)/
            (4*m*a2*r2*ct2+4*m*(-a2+r2)*rho2+(3+c2t)*a2*rho3-4*m*rho4+2*rho5+
            2*a4*ct2*rho)+(4*m*r*ct*(a2+rho2)*rhor)/(4*m*a2*r2*ct2+4*m*(-a2+
            r2)*rho2+(3+c2t)*a2*rho3-4*m*rho4+2*rho5+2*a4*ct2*rho)-4*m*r*ct*
            (a2+rho2)*pow(4*m*a2*r2*ct2+4*m*(-a2+r2)*rho2+(3+c2t)*a2*rho3-4*
            m*rho4+2*rho5+2*a4*ct2*rho,-2)*rho*(8*m*r*a2*ct2+8*m*r*rho2+2*a4*
            ct2*rhor+3*(3+c2t)*a2*rho2*rhor-16*m*rho3*rhor+10*rho4*rhor+8*m*
            (-a2+r2)*rho*rhor);
         gxxr=(-4*pow(a2+rho2,-2)*rho*(2*m*a4*r2*ct2+((5+c2t)*a2*rho5)/2.-2*
            m*rho6+rho7+2*m*a2*rho2*(-a2+2*r2*ct2+r2*cp2*st2)+2*m*rho4*(-2*a2+
            r2*ct2+r2*sp2*st2)+a6*ct2*rho+rho3*(a4+2*a4*ct2-2*a*m*r2*st2*
            s2p))*rhor)/(4*m*a2*r2*ct2+4*m*(-a2+r2)*rho2+(3+c2t)*a2*rho3-4*
            m*rho4+2*rho5+2*a4*ct2*rho)-(2*pow(4*m*a2*r2*ct2+4*m*(-a2+r2)*
            rho2+(3+c2t)*a2*rho3-4*m*rho4+2*rho5+2*a4*ct2*rho,-2)*(2*m*a4*r2*
            ct2+((5+c2t)*a2*rho5)/2.-2*m*rho6+rho7+2*m*a2*rho2*(-a2+2*r2*ct2+
            r2*cp2*st2)+2*m*rho4*(-2*a2+r2*ct2+r2*sp2*st2)+a6*ct2*rho+rho3*
            (a4+2*a4*ct2-2*a*m*r2*st2*s2p))*(8*m*r*a2*ct2+8*m*r*rho2+2*a4*ct2*
            rhor+3*(3+c2t)*a2*rho2*rhor-16*m*rho3*rhor+10*rho4*rhor+8*m*(-a2+
            r2)*rho*rhor))/(a2+rho2)+(2*(4*m*r*a4*ct2+2*m*a2*rho2*(4*r*ct2+2*
            r*cp2*st2)+2*m*rho4*(2*r*ct2+2*r*sp2*st2)-4*a*m*r*rho3*st2*s2p+
            a6*ct2*rhor+(5*(5+c2t)*a2*rho4*rhor)/2.-12*m*rho5*rhor+7*rho6*
            rhor+8*m*rho3*(-2*a2+r2*ct2+r2*sp2*st2)*rhor+4*m*a2*(-a2+2*r2*
            ct2+r2*cp2*st2)*rho*rhor+3*rho2*(a4+2*a4*ct2-2*a*m*r2*st2*s2p)*
            rhor))/((a2+rho2)*(4*m*a2*r2*ct2+4*m*(-a2+r2)*rho2+(3+c2t)*a2*
            rho3-4*m*rho4+2*rho5+2*a4*ct2*rho));
         gxyr=(-4*m*r*rho2*st2*(-2*a*c2p*rho-a2*s2p+rho2*s2p))/((a2+rho2)*(4*
            m*a2*r2*ct2+4*m*(-a2+r2)*rho2+(3+c2t)*a2*rho3-4*m*rho4+2*rho5+2*
            a4*ct2*rho))+(4*m*r2*pow(a2+rho2,-2)*rho3*st2*(-2*a*c2p*rho-a2*
            s2p+rho2*s2p)*rhor)/(4*m*a2*r2*ct2+4*m*(-a2+r2)*rho2+(3+c2t)*a2*
            rho3-4*m*rho4+2*rho5+2*a4*ct2*rho)-(4*m*r2*st2*rho*(-2*a*c2p*rho-
            a2*s2p+rho2*s2p)*rhor)/((a2+rho2)*(4*m*a2*r2*ct2+4*m*(-a2+r2)*
            rho2+(3+c2t)*a2*rho3-4*m*rho4+2*rho5+2*a4*ct2*rho))+(2*m*r2*rho2*
            pow(4*m*a2*r2*ct2+4*m*(-a2+r2)*rho2+(3+c2t)*a2*rho3-4*m*rho4+2*
            rho5+2*a4*ct2*rho,-2)*st2*(-2*a*c2p*rho-a2*s2p+rho2*s2p)*(8*m*r*
            a2*ct2+8*m*r*rho2+2*a4*ct2*rhor+3*(3+c2t)*a2*rho2*rhor-16*m*rho3*
            rhor+10*rho4*rhor+8*m*(-a2+r2)*rho*rhor))/(a2+rho2)-(2*m*r2*rho2*
            st2*(-2*a*c2p*rhor+2*rho*s2p*rhor))/((a2+rho2)*(4*m*a2*r2*ct2+4*m*
            (-a2+r2)*rho2+(3+c2t)*a2*rho3-4*m*rho4+2*rho5+2*a4*ct2*rho));
         gxzr=(-4*m*r*rho*(cp*rho+a*sp)*s2t)/(4*m*a2*r2*ct2+4*m*(-a2+r2)*rho2+
            (3+c2t)*a2*rho3-4*m*rho4+2*rho5+2*a4*ct2*rho)-(2*m*cp*r2*rho*s2t*
            rhor)/(4*m*a2*r2*ct2+4*m*(-a2+r2)*rho2+(3+c2t)*a2*rho3-4*m*rho4+2*
            rho5+2*a4*ct2*rho)-(2*m*r2*(cp*rho+a*sp)*s2t*rhor)/(4*m*a2*r2*ct2+
            4*m*(-a2+r2)*rho2+(3+c2t)*a2*rho3-4*m*rho4+2*rho5+2*a4*ct2*rho)+2*
            m*r2*pow(4*m*a2*r2*ct2+4*m*(-a2+r2)*rho2+(3+c2t)*a2*rho3-4*m*rho4+
            2*rho5+2*a4*ct2*rho,-2)*rho*(cp*rho+a*sp)*s2t*(8*m*r*a2*ct2+8*m*r*
            rho2+2*a4*ct2*rhor+3*(3+c2t)*a2*rho2*rhor-16*m*rho3*rhor+10*rho4*
            rhor+8*m*(-a2+r2)*rho*rhor);
         gyyr=(-4*pow(a2+rho2,-2)*rho*(2*m*a4*r2*ct2+((5+c2t)*a2*rho5)/2.-2*
            m*rho6+rho7+2*m*rho4*(-2*a2+r2*ct2+r2*cp2*st2)+2*m*a2*rho2*(-a2+2*
            r2*ct2+r2*sp2*st2)+a6*ct2*rho+rho3*(a4+2*a4*ct2+2*a*m*r2*st2*
            s2p))*rhor)/(4*m*a2*r2*ct2+4*m*(-a2+r2)*rho2+(3+c2t)*a2*rho3-4*m*
            rho4+2*rho5+2*a4*ct2*rho)-(2*pow(4*m*a2*r2*ct2+4*m*(-a2+r2)*rho2+
            (3+c2t)*a2*rho3-4*m*rho4+2*rho5+2*a4*ct2*rho,-2)*(2*m*a4*r2*ct2+
            ((5+c2t)*a2*rho5)/2.-2*m*rho6+rho7+2*m*rho4*(-2*a2+r2*ct2+r2*cp2*
            st2)+2*m*a2*rho2*(-a2+2*r2*ct2+r2*sp2*st2)+a6*ct2*rho+rho3*(a4+2*
            a4*ct2+2*a*m*r2*st2*s2p))*(8*m*r*a2*ct2+8*m*r*rho2+2*a4*ct2*rhor+
            3*(3+c2t)*a2*rho2*rhor-16*m*rho3*rhor+10*rho4*rhor+8*m*(-a2+r2)*
            rho*rhor))/(a2+rho2)+(2*(4*m*r*a4*ct2+2*m*rho4*(2*r*ct2+2*r*cp2*
            st2)+2*m*a2*rho2*(4*r*ct2+2*r*sp2*st2)+4*a*m*r*rho3*st2*s2p+a6*
            ct2*rhor+(5*(5+c2t)*a2*rho4*rhor)/2.-12*m*rho5*rhor+7*rho6*rhor+
            8*m*rho3*(-2*a2+r2*ct2+r2*cp2*st2)*rhor+4*m*a2*(-a2+2*r2*ct2+r2*
            sp2*st2)*rho*rhor+3*rho2*(a4+2*a4*ct2+2*a*m*r2*st2*s2p)*rhor))/
            ((a2+rho2)*(4*m*a2*r2*ct2+4*m*(-a2+r2)*rho2+(3+c2t)*a2*rho3-4*m*
            rho4+2*rho5+2*a4*ct2*rho));
         gyzr=(-8*m*r*ct*rho*(-(a*cp)+rho*sp)*st)/(4*m*a2*r2*ct2+4*m*(-a2+r2)*
            rho2+(3+c2t)*a2*rho3-4*m*rho4+2*rho5+2*a4*ct2*rho)-(4*m*ct*r2*rho*
            sp*st*rhor)/(4*m*a2*r2*ct2+4*m*(-a2+r2)*rho2+(3+c2t)*a2*rho3-4*m*
            rho4+2*rho5+2*a4*ct2*rho)-(4*m*ct*r2*(-(a*cp)+rho*sp)*st*rhor)/
            (4*m*a2*r2*ct2+4*m*(-a2+r2)*rho2+(3+c2t)*a2*rho3-4*m*rho4+2*rho5+
            2*a4*ct2*rho)+4*m*ct*r2*pow(4*m*a2*r2*ct2+4*m*(-a2+r2)*rho2+(3+
            c2t)*a2*rho3-4*m*rho4+2*rho5+2*a4*ct2*rho,-2)*rho*(-(a*cp)+rho*
            sp)*st*(8*m*r*a2*ct2+8*m*r*rho2+2*a4*ct2*rhor+3*(3+c2t)*a2*rho2*
            rhor-16*m*rho3*rhor+10*rho4*rhor+8*m*(-a2+r2)*rho*rhor);
         gzzr=((2*a4*ct2+(3+c2t)*a2*rho2-4*m*rho3+2*rho4-2*m*(2*a2-r2+c2t*r2)*
            rho)*rhor)/(4*m*a2*r2*ct2+4*m*(-a2+r2)*rho2+(3+c2t)*a2*rho3-4*m*
            rho4+2*rho5+2*a4*ct2*rho)+(rho*(-2*m*(-2*r+2*r*c2t)*rho-2*m*(2*a2-
            r2+c2t*r2)*rhor-12*m*rho2*rhor+8*rho3*rhor+2*(3+c2t)*a2*rho*
            rhor))/(4*m*a2*r2*ct2+4*m*(-a2+r2)*rho2+(3+c2t)*a2*rho3-4*m*rho4+
            2*rho5+2*a4*ct2*rho)-pow(4*m*a2*r2*ct2+4*m*(-a2+r2)*rho2+(3+c2t)*
            a2*rho3-4*m*rho4+2*rho5+2*a4*ct2*rho,-2)*rho*(2*a4*ct2+(3+c2t)*a2*
            rho2-4*m*rho3+2*rho4-2*m*(2*a2-r2+c2t*r2)*rho)*(8*m*r*a2*ct2+8*m*
            r*rho2+2*a4*ct2*rhor+3*(3+c2t)*a2*rho2*rhor-16*m*rho3*rhor+10*
            rho4*rhor+8*m*(-a2+r2)*rho*rhor);
         gttt=(-2*(a2*ct2+rho2)*(a2*rhot+3*rho2*rhot))/(4*m*a2*r2*ct2+4*m*
            (-a2+r2)*rho2+(3+c2t)*a2*rho3-4*m*rho4+2*rho5+2*a4*ct2*rho)-(2*
            (2*m*r2+rho3+a2*rho)*(-2*ct*a2*st+2*rho*rhot))/(4*m*a2*r2*ct2+4*
            m*(-a2+r2)*rho2+(3+c2t)*a2*rho3-4*m*rho4+2*rho5+2*a4*ct2*rho)+2*
            (a2*ct2+rho2)*pow(4*m*a2*r2*ct2+4*m*(-a2+r2)*rho2+(3+c2t)*a2*rho3-
            4*m*rho4+2*rho5+2*a4*ct2*rho,-2)*(2*m*r2+rho3+a2*rho)*(-8*m*ct*a2*
            r2*st-4*ct*a4*rho*st-2*a2*rho3*s2t+2*a4*ct2*rhot+3*(3+c2t)*a2*
            rho2*rhot-16*m*rho3*rhot+10*rho4*rhot+8*m*(-a2+r2)*rho*rhot);
         gxtt=(4*m*r*ct*rho2*(cp*rho+a*sp))/(4*m*a2*r2*ct2+4*m*(-a2+r2)*rho2+
            (3+c2t)*a2*rho3-4*m*rho4+2*rho5+2*a4*ct2*rho)+(4*m*r*cp*rho2*st*
            rhot)/(4*m*a2*r2*ct2+4*m*(-a2+r2)*rho2+(3+c2t)*a2*rho3-4*m*rho4+2*
            rho5+2*a4*ct2*rho)+(8*m*r*rho*(cp*rho+a*sp)*st*rhot)/(4*m*a2*r2*
            ct2+4*m*(-a2+r2)*rho2+(3+c2t)*a2*rho3-4*m*rho4+2*rho5+2*a4*ct2*
            rho)-4*m*r*rho2*pow(4*m*a2*r2*ct2+4*m*(-a2+r2)*rho2+(3+c2t)*a2*
            rho3-4*m*rho4+2*rho5+2*a4*ct2*rho,-2)*(cp*rho+a*sp)*st*(-8*m*ct*
            a2*r2*st-4*ct*a4*rho*st-2*a2*rho3*s2t+2*a4*ct2*rhot+3*(3+c2t)*a2*
            rho2*rhot-16*m*rho3*rhot+10*rho4*rhot+8*m*(-a2+r2)*rho*rhot);
         gytt=(4*m*r*ct*rho2*(-(a*cp)+rho*sp))/(4*m*a2*r2*ct2+4*m*(-a2+r2)*
            rho2+(3+c2t)*a2*rho3-4*m*rho4+2*rho5+2*a4*ct2*rho)+(4*m*r*rho2*
            sp*st*rhot)/(4*m*a2*r2*ct2+4*m*(-a2+r2)*rho2+(3+c2t)*a2*rho3-4*m*
            rho4+2*rho5+2*a4*ct2*rho)+(8*m*r*rho*(-(a*cp)+rho*sp)*st*rhot)/(4*
            m*a2*r2*ct2+4*m*(-a2+r2)*rho2+(3+c2t)*a2*rho3-4*m*rho4+2*rho5+2*
            a4*ct2*rho)-4*m*r*rho2*pow(4*m*a2*r2*ct2+4*m*(-a2+r2)*rho2+(3+
            c2t)*a2*rho3-4*m*rho4+2*rho5+2*a4*ct2*rho,-2)*(-(a*cp)+rho*sp)*st*
            (-8*m*ct*a2*r2*st-4*ct*a4*rho*st-2*a2*rho3*s2t+2*a4*ct2*rhot+3*(3+
            c2t)*a2*rho2*rhot-16*m*rho3*rhot+10*rho4*rhot+8*m*(-a2+r2)*rho*
            rhot);
         gztt=(-4*m*r*(a2+rho2)*rho*st)/(4*m*a2*r2*ct2+4*m*(-a2+r2)*rho2+(3+
            c2t)*a2*rho3-4*m*rho4+2*rho5+2*a4*ct2*rho)+(8*m*r*ct*rho2*rhot)/
            (4*m*a2*r2*ct2+4*m*(-a2+r2)*rho2+(3+c2t)*a2*rho3-4*m*rho4+2*rho5+
            2*a4*ct2*rho)+(4*m*r*ct*(a2+rho2)*rhot)/(4*m*a2*r2*ct2+4*m*(-a2+
            r2)*rho2+(3+c2t)*a2*rho3-4*m*rho4+2*rho5+2*a4*ct2*rho)-4*m*r*ct*
            (a2+rho2)*pow(4*m*a2*r2*ct2+4*m*(-a2+r2)*rho2+(3+c2t)*a2*rho3-4*m*
            rho4+2*rho5+2*a4*ct2*rho,-2)*rho*(-8*m*ct*a2*r2*st-4*ct*a4*rho*st-
            2*a2*rho3*s2t+2*a4*ct2*rhot+3*(3+c2t)*a2*rho2*rhot-16*m*rho3*rhot+
            10*rho4*rhot+8*m*(-a2+r2)*rho*rhot);
         gxxt=(-4*pow(a2+rho2,-2)*rho*(2*m*a4*r2*ct2+((5+c2t)*a2*rho5)/2.-2*m*
            rho6+rho7+2*m*a2*rho2*(-a2+2*r2*ct2+r2*cp2*st2)+2*m*rho4*(-2*a2+
            r2*ct2+r2*sp2*st2)+a6*ct2*rho+rho3*(a4+2*a4*ct2-2*a*m*r2*st2*
            s2p))*rhot)/(4*m*a2*r2*ct2+4*m*(-a2+r2)*rho2+(3+c2t)*a2*rho3-4*m*
            rho4+2*rho5+2*a4*ct2*rho)-(2*pow(4*m*a2*r2*ct2+4*m*(-a2+r2)*rho2+
            (3+c2t)*a2*rho3-4*m*rho4+2*rho5+2*a4*ct2*rho,-2)*(2*m*a4*r2*ct2+
            ((5+c2t)*a2*rho5)/2.-2*m*rho6+rho7+2*m*a2*rho2*(-a2+2*r2*ct2+r2*
            cp2*st2)+2*m*rho4*(-2*a2+r2*ct2+r2*sp2*st2)+a6*ct2*rho+rho3*(a4+
            2*a4*ct2-2*a*m*r2*st2*s2p))*(-8*m*ct*a2*r2*st-4*ct*a4*rho*st-2*a2*
            rho3*s2t+2*a4*ct2*rhot+3*(3+c2t)*a2*rho2*rhot-16*m*rho3*rhot+10*
            rho4*rhot+8*m*(-a2+r2)*rho*rhot))/(a2+rho2)+(2*(-4*m*ct*a4*r2*st-
            2*ct*a6*rho*st+2*m*a2*rho2*(-4*ct*r2*st+2*ct*r2*cp2*st)+2*m*rho4*
            (-2*ct*r2*st+2*ct*r2*sp2*st)+rho3*(-4*ct*a4*st-4*a*m*ct*r2*s2p*
            st)-a2*rho5*s2t+a6*ct2*rhot+(5*(5+c2t)*a2*rho4*rhot)/2.-12*m*rho5*
            rhot+7*rho6*rhot+8*m*rho3*(-2*a2+r2*ct2+r2*sp2*st2)*rhot+4*m*a2*
            (-a2+2*r2*ct2+r2*cp2*st2)*rho*rhot+3*rho2*(a4+2*a4*ct2-2*a*m*r2*
            st2*s2p)*rhot))/((a2+rho2)*(4*m*a2*r2*ct2+4*m*(-a2+r2)*rho2+(3+
            c2t)*a2*rho3-4*m*rho4+2*rho5+2*a4*ct2*rho));
         gxyt=(-4*m*ct*r2*rho2*(-2*a*c2p*rho-a2*s2p+rho2*s2p)*st)/((a2+rho2)*
            (4*m*a2*r2*ct2+4*m*(-a2+r2)*rho2+(3+c2t)*a2*rho3-4*m*rho4+2*rho5+
            2*a4*ct2*rho))+(4*m*r2*pow(a2+rho2,-2)*rho3*st2*(-2*a*c2p*rho-a2*
            s2p+rho2*s2p)*rhot)/(4*m*a2*r2*ct2+4*m*(-a2+r2)*rho2+(3+c2t)*a2*
            rho3-4*m*rho4+2*rho5+2*a4*ct2*rho)-(4*m*r2*st2*rho*(-2*a*c2p*rho-
            a2*s2p+rho2*s2p)*rhot)/((a2+rho2)*(4*m*a2*r2*ct2+4*m*(-a2+r2)*
            rho2+(3+c2t)*a2*rho3-4*m*rho4+2*rho5+2*a4*ct2*rho))+(2*m*r2*rho2*
            pow(4*m*a2*r2*ct2+4*m*(-a2+r2)*rho2+(3+c2t)*a2*rho3-4*m*rho4+2*
            rho5+2*a4*ct2*rho,-2)*st2*(-2*a*c2p*rho-a2*s2p+rho2*s2p)*(-8*m*ct*
            a2*r2*st-4*ct*a4*rho*st-2*a2*rho3*s2t+2*a4*ct2*rhot+3*(3+c2t)*a2*
            rho2*rhot-16*m*rho3*rhot+10*rho4*rhot+8*m*(-a2+r2)*rho*rhot))/
            (a2+rho2)-(2*m*r2*rho2*st2*(-2*a*c2p*rhot+2*rho*s2p*rhot))/
            ((a2+rho2)*(4*m*a2*r2*ct2+4*m*(-a2+r2)*rho2+(3+c2t)*a2*rho3-4*m*
            rho4+2*rho5+2*a4*ct2*rho));
         gxzt=(-4*m*c2t*r2*rho*(cp*rho+a*sp))/(4*m*a2*r2*ct2+4*m*(-a2+r2)*
            rho2+(3+c2t)*a2*rho3-4*m*rho4+2*rho5+2*a4*ct2*rho)-(2*m*cp*r2*rho*
            s2t*rhot)/(4*m*a2*r2*ct2+4*m*(-a2+r2)*rho2+(3+c2t)*a2*rho3-4*m*
            rho4+2*rho5+2*a4*ct2*rho)-(2*m*r2*(cp*rho+a*sp)*s2t*rhot)/(4*m*a2*
            r2*ct2+4*m*(-a2+r2)*rho2+(3+c2t)*a2*rho3-4*m*rho4+2*rho5+2*a4*ct2*
            rho)+2*m*r2*pow(4*m*a2*r2*ct2+4*m*(-a2+r2)*rho2+(3+c2t)*a2*rho3-4*
            m*rho4+2*rho5+2*a4*ct2*rho,-2)*rho*(cp*rho+a*sp)*s2t*(-8*m*ct*a2*
            r2*st-4*ct*a4*rho*st-2*a2*rho3*s2t+2*a4*ct2*rhot+3*(3+c2t)*a2*
            rho2*rhot-16*m*rho3*rhot+10*rho4*rhot+8*m*(-a2+r2)*rho*rhot);
         gyyt=(-4*pow(a2+rho2,-2)*rho*(2*m*a4*r2*ct2+((5+c2t)*a2*rho5)/2.-2*m*
            rho6+rho7+2*m*rho4*(-2*a2+r2*ct2+r2*cp2*st2)+2*m*a2*rho2*(-a2+2*
            r2*ct2+r2*sp2*st2)+a6*ct2*rho+rho3*(a4+2*a4*ct2+2*a*m*r2*st2*
            s2p))*rhot)/(4*m*a2*r2*ct2+4*m*(-a2+r2)*rho2+(3+c2t)*a2*rho3-4*m*
            rho4+2*rho5+2*a4*ct2*rho)-(2*pow(4*m*a2*r2*ct2+4*m*(-a2+r2)*rho2+
            (3+c2t)*a2*rho3-4*m*rho4+2*rho5+2*a4*ct2*rho,-2)*(2*m*a4*r2*ct2+
            ((5+c2t)*a2*rho5)/2.-2*m*rho6+rho7+2*m*rho4*(-2*a2+r2*ct2+r2*cp2*
            st2)+2*m*a2*rho2*(-a2+2*r2*ct2+r2*sp2*st2)+a6*ct2*rho+rho3*(a4+2*
            a4*ct2+2*a*m*r2*st2*s2p))*(-8*m*ct*a2*r2*st-4*ct*a4*rho*st-2*a2*
            rho3*s2t+2*a4*ct2*rhot+3*(3+c2t)*a2*rho2*rhot-16*m*rho3*rhot+10*
            rho4*rhot+8*m*(-a2+r2)*rho*rhot))/(a2+rho2)+(2*(-4*m*ct*a4*r2*st-
            2*ct*a6*rho*st+2*m*rho4*(-2*ct*r2*st+2*ct*r2*cp2*st)+2*m*a2*rho2*
            (-4*ct*r2*st+2*ct*r2*sp2*st)+rho3*(-4*ct*a4*st+4*a*m*ct*r2*s2p*
            st)-a2*rho5*s2t+a6*ct2*rhot+(5*(5+c2t)*a2*rho4*rhot)/2.-12*m*rho5*
            rhot+7*rho6*rhot+8*m*rho3*(-2*a2+r2*ct2+r2*cp2*st2)*rhot+4*m*a2*
            (-a2+2*r2*ct2+r2*sp2*st2)*rho*rhot+3*rho2*(a4+2*a4*ct2+2*a*m*r2*
            st2*s2p)*rhot))/((a2+rho2)*(4*m*a2*r2*ct2+4*m*(-a2+r2)*rho2+(3+
            c2t)*a2*rho3-4*m*rho4+2*rho5+2*a4*ct2*rho));
         gyzt=(-4*m*r2*ct2*rho*(-(a*cp)+rho*sp))/(4*m*a2*r2*ct2+4*m*(-a2+r2)*
            rho2+(3+c2t)*a2*rho3-4*m*rho4+2*rho5+2*a4*ct2*rho)+(4*m*r2*st2*
            rho*(-(a*cp)+rho*sp))/(4*m*a2*r2*ct2+4*m*(-a2+r2)*rho2+(3+c2t)*a2*
            rho3-4*m*rho4+2*rho5+2*a4*ct2*rho)-(4*m*ct*r2*rho*sp*st*rhot)/
            (4*m*a2*r2*ct2+4*m*(-a2+r2)*rho2+(3+c2t)*a2*rho3-4*m*rho4+2*rho5+
            2*a4*ct2*rho)-(4*m*ct*r2*(-(a*cp)+rho*sp)*st*rhot)/(4*m*a2*r2*ct2+
            4*m*(-a2+r2)*rho2+(3+c2t)*a2*rho3-4*m*rho4+2*rho5+2*a4*ct2*rho)+4*
            m*ct*r2*pow(4*m*a2*r2*ct2+4*m*(-a2+r2)*rho2+(3+c2t)*a2*rho3-4*m*
            rho4+2*rho5+2*a4*ct2*rho,-2)*rho*(-(a*cp)+rho*sp)*st*(-8*m*ct*a2*
            r2*st-4*ct*a4*rho*st-2*a2*rho3*s2t+2*a4*ct2*rhot+3*(3+c2t)*a2*
            rho2*rhot-16*m*rho3*rhot+10*rho4*rhot+8*m*(-a2+r2)*rho*rhot);
         gzzt=((2*a4*ct2+(3+c2t)*a2*rho2-4*m*rho3+2*rho4-2*m*(2*a2-r2+c2t*r2)*
            rho)*rhot)/(4*m*a2*r2*ct2+4*m*(-a2+r2)*rho2+(3+c2t)*a2*rho3-4*m*
            rho4+2*rho5+2*a4*ct2*rho)+(rho*(-4*ct*a4*st-2*a2*rho2*s2t+4*m*r2*
            rho*s2t-2*m*(2*a2-r2+c2t*r2)*rhot-12*m*rho2*rhot+8*rho3*rhot+2*(3+
            c2t)*a2*rho*rhot))/(4*m*a2*r2*ct2+4*m*(-a2+r2)*rho2+(3+c2t)*a2*
            rho3-4*m*rho4+2*rho5+2*a4*ct2*rho)-pow(4*m*a2*r2*ct2+4*m*(-a2+r2)*
            rho2+(3+c2t)*a2*rho3-4*m*rho4+2*rho5+2*a4*ct2*rho,-2)*rho*(2*a4*
            ct2+(3+c2t)*a2*rho2-4*m*rho3+2*rho4-2*m*(2*a2-r2+c2t*r2)*rho)*(-8*
            m*ct*a2*r2*st-4*ct*a4*rho*st-2*a2*rho3*s2t+2*a4*ct2*rhot+3*(3+
            c2t)*a2*rho2*rhot-16*m*rho3*rhot+10*rho4*rhot+8*m*(-a2+r2)*rho*
            rhot);
         gttp=0;
         gxtp=(4*m*r*rho2*(a*cp-rho*sp)*st)/(4*m*a2*r2*ct2+4*m*(-a2+r2)*rho2+
            (3+c2t)*a2*rho3-4*m*rho4+2*rho5+2*a4*ct2*rho);
         gytp=(4*m*r*rho2*(cp*rho+a*sp)*st)/(4*m*a2*r2*ct2+4*m*(-a2+r2)*rho2+
            (3+c2t)*a2*rho3-4*m*rho4+2*rho5+2*a4*ct2*rho);
         gztp=0;
         gxxp=(2*(-4*a*m*c2p*r2*rho3*st2-4*m*cp*a2*r2*rho2*st2*sp+4*m*cp*r2*
            rho4*st2*sp))/((a2+rho2)*(4*m*a2*r2*ct2+4*m*(-a2+r2)*rho2+(3+c2t)*
            a2*rho3-4*m*rho4+2*rho5+2*a4*ct2*rho));
         gxyp=(-2*m*r2*rho2*st2*(-2*c2p*a2+2*c2p*rho2+4*a*rho*s2p))/((a2+
            rho2)*(4*m*a2*r2*ct2+4*m*(-a2+r2)*rho2+(3+c2t)*a2*rho3-4*m*rho4+2*
            rho5+2*a4*ct2*rho));
         gxzp=(-2*m*r2*rho*(a*cp-rho*sp)*s2t)/(4*m*a2*r2*ct2+4*m*(-a2+r2)*
            rho2+(3+c2t)*a2*rho3-4*m*rho4+2*rho5+2*a4*ct2*rho);
         gyyp=(2*(4*a*m*c2p*r2*rho3*st2+4*m*cp*a2*r2*rho2*st2*sp-4*m*cp*r2*
            rho4*st2*sp))/((a2+rho2)*(4*m*a2*r2*ct2+4*m*(-a2+r2)*rho2+(3+c2t)*
            a2*rho3-4*m*rho4+2*rho5+2*a4*ct2*rho));
         gyzp=(-4*m*ct*r2*rho*(cp*rho+a*sp)*st)/(4*m*a2*r2*ct2+4*m*(-a2+r2)*
            rho2+(3+c2t)*a2*rho3-4*m*rho4+2*rho5+2*a4*ct2*rho);
         gzzp=0;

   /* Derivative of the spatial determinant */

      deth=-(4*m*a2*r2*ct2+4*m*(-a2+r2)*rho2+(3+c2t)*a2*rho3-4*m*rho4+2*rho5+
         2*a4*ct2*rho)/(2.*(a2+rho2)*(a2*ct2+rho2)*rho);

      dhdr=m*pow(a2+rho2,-2)*pow(a2*ct2+rho2,-2)*pow(rho,-2)*(-2*r*(5+c2t)*a4*
         ct2*rho3-4*r*(2+c2t)*a2*rho5-4*r*rho7-4*r*a6*pow(ct,4)*rho+2*a6*r2*
         pow(ct,4)*rhor+a4*(2*a2+7*r2+3*c2t*r2)*ct2*rho2*rhor+2*a2*(4*r2+c2t*
         (a2+3*r2))*rho4*rhor+(-3*a2+c2t*a2+6*r2)*rho6*rhor-2*pow(rho,8)*
         rhor);
      dhdt=m*pow(a2+rho2,-2)*pow(a2*ct2+rho2,-2)*pow(rho,-2)*(2*a6*rho3*s2t+4*
         a4*rho5*s2t+2*a2*rho7*s2t+2*a6*r2*pow(ct,4)*rhot+a4*(2*a2+7*r2+3*c2t*
         r2)*ct2*rho2*rhot+2*a2*(4*r2+c2t*(a2+3*r2))*rho4*rhot+(-3*a2+c2t*a2+
         6*r2)*rho6*rhot-2*pow(rho,8)*rhot);
      dhdp=0;

	/* The Jacobian */

	drdx  = st*cp;
	drdy  = st*sp;
	drdz  = ct;

	dtdx  = ct*cp/r;
	dtdy  = ct*sp/r;
	dtdz  = -st/r;
	
	dpdx  = -sp/st/r;
	dpdy  = cp/st/r;
	dpdz  = 0.;


   /* Partial derivative of the logarithm of the spatial determinant */

   hdetx = (drdx*dhdr + dtdx*dhdt + dpdx*dhdp)/deth;
   hdety = (drdy*dhdr + dtdy*dhdt + dpdy*dhdp)/deth;
   hdetz = (drdz*dhdr + dtdz*dhdt + dpdz*dhdp)/deth;

	/* Partial derivatives (x,y,z) of the inverse 
	   cartesian-like metric */

	gxxx  = gxxr*drdx + gxxt*dtdx + gxxp*dpdx;
	gxxy  = gxxr*drdy + gxxt*dtdy + gxxp*dpdy;
	gxxz  = gxxr*drdz + gxxt*dtdz + gxxp*dpdz;
	
	gxyx  = gxyr*drdx + gxyt*dtdx + gxyp*dpdx;
	gxyy  = gxyr*drdy + gxyt*dtdy + gxyp*dpdy;
	gxyz  = gxyr*drdz + gxyt*dtdz + gxyp*dpdz;
	
	gxzx  = gxzr*drdx + gxzt*dtdx + gxzp*dpdx;
	gxzy  = gxzr*drdy + gxzt*dtdy + gxzp*dpdy;
	gxzz  = gxzr*drdz + gxzt*dtdz + gxzp*dpdz;
	
	gxtx  = gxtr*drdx + gxtt*dtdx + gxtp*dpdx;
	gxty  = gxtr*drdy + gxtt*dtdy + gxtp*dpdy;
	gxtz  = gxtr*drdz + gxtt*dtdz + gxtp*dpdz;
	
	gyyx  = gyyr*drdx + gyyt*dtdx + gyyp*dpdx;
	gyyy  = gyyr*drdy + gyyt*dtdy + gyyp*dpdy;
	gyyz  = gyyr*drdz + gyyt*dtdz + gyyp*dpdz;
	
	gyzx  = gyzr*drdx + gyzt*dtdx + gyzp*dpdx;
	gyzy  = gyzr*drdy + gyzt*dtdy + gyzp*dpdy;
	gyzz  = gyzr*drdz + gyzt*dtdz + gyzp*dpdz;
	
	gytx  = gytr*drdx + gytt*dtdx + gytp*dpdx;
	gyty  = gytr*drdy + gytt*dtdy + gytp*dpdy;
	gytz  = gytr*drdz + gytt*dtdz + gytp*dpdz;
	
	gzzx  = gzzr*drdx + gzzt*dtdx + gzzp*dpdx;
	gzzy  = gzzr*drdy + gzzt*dtdy + gzzp*dpdy;
	gzzz  = gzzr*drdz + gzzt*dtdz + gzzp*dpdz;
	
	gztx  = gztr*drdx + gztt*dtdx + gztp*dpdx;
	gzty  = gztr*drdy + gztt*dtdy + gztp*dpdy;
	gztz  = gztr*drdz + gztt*dtdz + gztp*dpdz;
	
	gttx  = gttr*drdx + gttt*dtdx + gttp*dpdx;
	gtty  = gttr*drdy + gttt*dtdy + gttp*dpdy;
	gttz  = gttr*drdz + gttt*dtdz + gttp*dpdz;
	
	ux = p->mom[0];
	uy = p->mom[1];
	uz = p->mom[2];
        ut = p->mom[3];
        uut =  ut*p->gutt+ux*p->guxt+uy*p->guyt+uz*p->guzt;

	curx = ut*ut*gttx+ux*ux*gxxx+uy*uy*gyyx+uz*uz*gzzx
	  +(float)2.0*(ux*uy*gxyx+ux*uz*gxzx+ux*ut*gxtx+
                       uy*uz*gyzx+uy*ut*gytx+uz*ut*gztx);
	cury = ut*ut*gtty+ux*ux*gxxy+uy*uy*gyyy+uz*uz*gzzy
	  +(float)2.0*(ux*uy*gxyy+ux*uz*gxzy+ux*ut*gxty+
                       uy*uz*gyzy+uy*ut*gyty+uz*ut*gzty);
	curz = ut*ut*gttz+ux*ux*gxxz+uy*uy*gyyz+uz*uz*gzzz
	  +(float)2.0*(ux*uy*gxyz+ux*uz*gxzz+ux*ut*gxtz+
                       uy*uz*gyzz+uy*ut*gytz+uz*ut*gztz);

	p->acc[0] += -0.5/uut*curx - p->alfa*p->pr/p->rho*hdetx;
	p->acc[1] += -0.5/uut*cury - p->alfa*p->pr/p->rho*hdety;
	p->acc[2] += -0.5/uut*curz - p->alfa*p->pr/p->rho*hdetz;

 /* Roche stuff - Karen */
/*        muk = 1./(8.*8.*8.);
        q = 1.4/4.5;
        tosc = 3.e5*q;
        trelax = 0.9*tosc;

        p->acc[0] += muk*(3.+q)*(x-8.);
        p->acc[1] += muk*q*y;
        p->acc[2] += -muk*z; */

/* Damping for relaxation - note that this is not general. */
/*        p->acc[0] += -vx/100.;
        p->acc[1] += -vy/100.;
        p->acc[2] += -vz/100.; */

    }
}

