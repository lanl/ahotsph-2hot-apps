#include "singlio.h"
#include <math.h>
#include "physics_sph.h"
#include "vop.h"

void
initial_cond(body *btab, int nobj, 
        float xx0, float yy0, float zz0,
        float vx0, float vy0, float vz0, 
        float bhmass, float Gamma)
{
    body *p;

    float vt, vx, vy, vz;
    float ut, ux, uy, uz;
    float v2, uut;

    float umass, udist, utime, uvel, ubden, ueden;
    float hp, hx;

/* Convert from solar units to black hole units - Don */

    umass = 1.0/bhmass;
    udist = 1.0/bhmass;
    utime = sqrt(udist*udist*udist/umass);
    uvel  = udist/utime;
    ubden = umass/udist/udist/udist; /* Corrected ubden - Don */
    ueden = uvel*uvel;

    for (p = btab; p < btab+nobj; p++) {
         VS(p->pos, *= udist);
         p->h *= udist;
         VS(p->vel, *= uvel);
         p->mass *= umass;
         p->rho *= ubden;
         p->u *= ueden;
    }

    /* Move star to initial position and initial velocity */
      
    for (p = btab; p < btab+nobj; p++) {

        p->pos[0] += xx0;
        p->pos[1] += yy0;
        p->pos[2] += zz0;
        p->vel[0] += vx0;
        p->vel[1] += vy0;
        p->vel[2] += vz0;

    }

    /* Compute initial metric */

    get_metric(btab, nobj);

    /* Compute four-velocity and other stuff */

    for (p = btab; p < btab+nobj; p++) {

        vx = p->vel[0];
        vy = p->vel[1];
        vz = p->vel[2];
        vt = (float)1.0;

        v2 = vt*vt*p->gtt + vx*vx*p->gxx 
           + vy*vy*p->gyy + vz*vz*p->gzz
           + 2.0
           *(vt*vx*p->gxt + vt*vy*p->gyt 
           + vt*vz*p->gzt + vx*vy*p->gxy
           + vx*vz*p->gxz + vy*vz*p->gyz);

        uut = sqrt(-(float)1.0/v2);

        ut = vt*p->gtt + vx*p->gxt + vy*p->gyt + vz*p->gzt;
        ux = vt*p->gxt + vx*p->gxx + vy*p->gxy + vz*p->gxz;
        uy = vt*p->gyt + vx*p->gxy + vy*p->gyy + vz*p->gyz;
        uz = vt*p->gzt + vx*p->gxz + vy*p->gyz + vz*p->gzz;

        ut *= uut;
        ux *= uut;
        uy *= uut;
        uz *= uut;

        p->alfa = sqrt(-(float)1.0/p->gutt);
        p->gama = uut*p->alfa;
        p->pr   = p->u*p->rho*(Gamma-1.0);
        p->enth = 1.0+p->u*Gamma;
        p->rho  = p->rho*p->gama;

        p->mom[0] = p->enth*ux;
        p->mom[1] = p->enth*uy;
        p->mom[2] = p->enth*uz;
        p->mom[3] = p->enth*ut;
	
    }
}
