From pablo Fri Feb  4 12:59:20 1994
Return-Path: <pablo>
Date: Fri, 4 Feb 94 11:46:12 MST
From: pablo (pablo laguna)
To: msw

#include "physics_sph.h"
#include "fastflpt.h"

static float gamma = (float)(5.0/3.0);

void
quepasa(body *btab, int nobj)
{
    body *p;
    float rho_lab, masa;
    float star_mass = 0.0;
    float xcm = 0.0;
    float ycm = 0.0;
    float zcm = 0.0;
    float rcm;
    float gm1 = gamma-(float)1.0;
    float density_max = 0.0;                /* Maximum Density */
    float kelvin_max = 0.0;                 /* Maximum temperature */
    float entropy_tot = 0.0;                /* Total entropy */
    float energy_int = 0.0;                 /* Internal energy */
    float energy_sg = 0.0;                  /* Self-gravity energy */
    float energy_gr  = 0.0;                 /* GR (kinetic+BH potential) energy */
    float energy_tot = 0.0;                 /* Total energy */
    

    for (p = btab; p < btab+nobj; p++) {

         masa = p->mass;
         rho_lab = p->rho * p->alfa / p->gama;
         density_max = max(density_max,rho_lab);
         kelvin_max = max(kelvin_max,p->u * gm1);
/*         entropy_tot += masa * (float)1.0 / gm1 * alog( p->pr / rho_lab**gamma ); Karen */
           entropy_tot += masa * (float)1.0 / gm1 * alog( p->pr / pow(rho_lab,gamma) ); 
         energy_int += masa * p->u;
         energy_sg  += masa * p->phi;
         energy_gr  += masa * ( (p->gama*p->alfa) * (p->gama*p->alfa) 
                    - (float)1.0 ) / (float)2.0;
         xcm += masa * p->pos[0];
         ycm += masa * p->pos[1];
         zcm += masa * p->pos[2];
         star_mass += masa;

    }
    energy_tot = energy_gr+energy_int+energy_sg;
    rcm = sqrtf(xcm*xcm + ycm*ycm + zcm*zcm)/star_mass;

/* Mike: print energy_tot, entrophy_tot, density_max, kelvin_max, and rcm */

}




From pablo Fri Feb  4 12:59:29 1994
Return-Path: <pablo>
Date: Fri, 4 Feb 94 11:46:12 MST
From: pablo (pablo laguna)
To: msw

#include "physics_sph.h"
#include "fastflpt.h"
#include <math.h>

static float gamma = (float)(5.0/3.0);

void
quepasa(body *btab, int nobj)
{
    body *p;
    float rho_lab, masa;
    float star_mass = 0.0;
    float xcm = 0.0;
    float ycm = 0.0;
    float zcm = 0.0;
    float rcm;
    float gm1 = gamma-(float)1.0;
    float density_max = 0.0;                /* Maximum Density */
    float kelvin_max = 0.0;                 /* Maximum temperature */
    float entropy_tot = 0.0;                /* Total entropy */
    float energy_int = 0.0;                 /* Internal energy */
    float energy_sg = 0.0;                  /* Self-gravity energy */
    float energy_gr  = 0.0;                 /* GR (kinetic+BH potential) energy */
    float energy_tot = 0.0;                 /* Total energy */
    

    for (p = btab; p < btab+nobj; p++) {

         masa = p->mass;
         rho_lab = p->rho * p->alfa / p->gama;
         density_max = max(density_max,rho_lab);
         kelvin_max = max(kelvin_max,p->u * gm1);
         entropy_tot += masa * (float)1.0 / gm1 * alog( p->pr / rho_lab**gamma );
         energy_int += masa * p->u;
         energy_sg  += masa * p->phi;
         energy_gr  += masa * ( (p->gama*p->alfa) * (p->gama*p->alfa) 
                    - (float)1.0 ) / (float)2.0;
         xcm += masa * p->pos[0];
         ycm += masa * p->pos[1];
         zcm += masa * p->pos[2];
         star_mass += masa;

    }
    energy_tot = energy_gr+energy_int+energy_sg;
    rcm = sqrtf(xcm*xcm + ycm*ycm + zcm*zcm)/star_mass;

/* Mike: print energy_tot, entrophy_tot, density_max, kelvin_max, and rcm */

}




