#include "physics_sph.h"

void grav_rad(body *btab, int nobj, float *hx, float *hp);
int df(int x, int y);

void grav_rad(body *btab, int nobj, float *hx, float *hp)
{
   body *p;
	float x[3];	/* Position */
	float v[3];	/* Velocity */
	float g[3];	/* Acceleration */
	float ms,P,dns;	/* Mass, pressure, and density */ 
	int l, m;		/* Tensor indices */
	float Itt[NDIM][NDIM];	/* Components of the second time derviative
				 of the quadrupole moment */
/*	float hp, hx;*/		/* The two polarazation states of the outgoing
				 gravitational radiation */

	for(l=0;l<NDIM;l++)
		for(m=0;m<NDIM;m++)
			Itt[l][m] = 0.0;

   for (p = btab; p < btab+nobj; p++) {

        x[0] = p->pos[0];
        x[1] = p->pos[1];
        x[2] = p->pos[2];

        v[0] = p->vel[0];
        v[1] = p->vel[1];
        v[2] = p->vel[2];

	g[0] = p->acc[0];
	g[1] = p->acc[1];
	g[2] = p->acc[2]; 

	ms = p->mass;
	P = p->pr;
	dns = p->rho;

/* Determine the quadrupole moments */
	for(l=0;l<NDIM;l++)
		for(m=0;m<NDIM;m++)
/* Karen: fixing what I think is a bug... */
/*			Itt[l][m] += m*(2*v[l]*v[m]
				+ 2*P/dns*df(l,m)
				+ (x[m]*g[l]+x[l]*g[m]) ); */
			Itt[l][m] += ms*(2*v[l]*v[m]
				+ 2*P/dns*df(l,m)
				+ (x[m]*g[l]+x[l]*g[m]) );
   }
   (*hp) = 2.0*(Itt[0][0]-Itt[1][1]);	 /* h+ polarization */
   (*hx) = 4.0*Itt[0][1];		 /* hx polarization */
}

/* Delta function */
int df(int x, int y)
{
	if( x == y ) return 1;
		else return 0;
}
