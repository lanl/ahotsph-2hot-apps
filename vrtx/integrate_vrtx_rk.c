/*
 * Copyright 1992 Michael S. Warren, John K. Salmon, and
 * Gregoire S. Winckelmans. All Rights Reserved.
 */

/* Nov. 94 Add Pedrizzetti's divergence filtering method */


#include "physics_vrtx.h"
#include "vop.h"
#include "fastflpt.h"


void
Update(bodyptr btab, int n, float dta, float dtb, float dtrel, int iflag)
{
    int i;
    bodyptr bp;
    float omega[3], omega2, str2, term;

    for(i=0; i<n; i++)
    {
	bp = btab+i;


	/* Vel_old could be NaN */
	if (dtb != 0.0)
	  VVV(Pos(bp), += dta * Vel(bp), + dtb * Vel_old(bp));
	else
	  VV(Pos(bp), += dta * Vel(bp));

	if (!finite(bp->pos[0])) 
	  Error("pos0 is infinite! particle %d of %d\n%g %g %g\n",
		i, n, bp->pos[0], bp->vel[0], bp->vel_old[0]);
	if (!finite(bp->pos[1])) 
	  Error("pos1 is infinite! particle %d of %d\n%g %g %g\n",
		i, n, bp->pos[1], bp->vel[1], bp->vel_old[1]);
	if (!finite(bp->pos[2])) 
	  Error("pos2 is infinite! particle %d of %d\n%g %g %g\n",
		i, n, bp->pos[2], bp->vel[2], bp->vel_old[2]);

        if(iflag){
	  VV(Vel_old(bp), = Vel(bp));
        }

#if 0
/* Add the stretching term to the viscous interaction. The symmetric
   scheme is used */

	Dstr(bp)[0] += 0.5F *( (Gradvel(bp)[0][0]+Gradvel(bp)[0][0])
	         		             *Strength(bp)[0]
	               	      +(Gradvel(bp)[0][1]+Gradvel(bp)[1][0])
				             *Strength(bp)[1]
			      +(Gradvel(bp)[0][2]+Gradvel(bp)[2][0])
				             *Strength(bp)[2] );

	Dstr(bp)[1] += 0.5F *( (Gradvel(bp)[1][0]+Gradvel(bp)[0][1])
		             		     *Strength(bp)[0]
			      +(Gradvel(bp)[1][1]+Gradvel(bp)[1][1])
				             *Strength(bp)[1]
			      +(Gradvel(bp)[1][2]+Gradvel(bp)[2][1])
				             *Strength(bp)[2] );

	Dstr(bp)[2] += 0.5F *( (Gradvel(bp)[2][0]+Gradvel(bp)[0][2])
				             *Strength(bp)[0]
			      +(Gradvel(bp)[2][1]+Gradvel(bp)[1][2])
				             *Strength(bp)[1]
			      +(Gradvel(bp)[2][2]+Gradvel(bp)[2][2])
				             *Strength(bp)[2] );
#endif

/* Add the stretching term to the viscous interaction. The transpose 
   scheme is used */

        Dstr(bp)[0] +=   Gradvel(bp)[0][0]*Strength(bp)[0]
                        +Gradvel(bp)[1][0]*Strength(bp)[1]
                        +Gradvel(bp)[2][0]*Strength(bp)[2];

        Dstr(bp)[1] +=   Gradvel(bp)[0][1]*Strength(bp)[0]
                        +Gradvel(bp)[1][1]*Strength(bp)[1]
                        +Gradvel(bp)[2][1]*Strength(bp)[2];

        Dstr(bp)[2] +=   Gradvel(bp)[0][2]*Strength(bp)[0]
                        +Gradvel(bp)[1][2]*Strength(bp)[1]
                        +Gradvel(bp)[2][2]*Strength(bp)[2];


/* Relax the particle vorticity divergence using Predizzetti's scheme, 
  Fluid Dyn. Res. 10, 101 (1992) */

        omega[0] = Gradvel(bp)[2][1]-Gradvel(bp)[1][2];
        omega[1] = Gradvel(bp)[0][2]-Gradvel(bp)[2][0];
        omega[2] = Gradvel(bp)[1][0]-Gradvel(bp)[0][1];

        omega2 = Dot(omega, omega);
        if(omega2 != 0.F){ 

            str2 = Dot(Strength(bp), Strength(bp));
            term= sqrt(str2/omega2);
           
            VVV(Strength(bp), = (1.F-dtrel)*Strength(bp), + dtrel*term*omega); 
        } 

	if (dtb != 0.0)
	  VVV(Strength(bp), += dta * Dstr(bp), + dtb * Dstr_old(bp));
	else
	  VV(Strength(bp), += dta * Dstr(bp));

        if(iflag){
	  VV(Dstr_old(bp), = Dstr(bp));
        }
        
    }

}

