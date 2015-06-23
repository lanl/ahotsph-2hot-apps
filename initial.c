#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "stk.h"
#include "bigmalloc.h"
#include "randoms.h"
#include "singlio.h"
#include "ranlib.h"


#include "physics_sph.h"
#include "Msgs.h"
#include "timers.h"
#include "vop.h"
#include "fastflpt.h"
#include "mpmy.h"
#include "gc.h"
#include "wvt.h"
#include "polint.h"
#include "initial.h"

#define PI 3.141592653589793238462
#ifndef M_1_PI
#define	M_1_PI 0.31830988618379067154
#endif
/*#define NKERNEL_TABLE 40000*/
#define NKERNEL_TABLE 40000
#define MAX_INDEX (NKERNEL_TABLE+2)

#define NO_UPDATE 2

void locate(double xx[], unsigned long n, double x, int *j);


/* From numerical recipes */
void locate(double xx[], unsigned long n, double x, int *j)
{
    unsigned long ju,jm,jl;
    int ascnd;
    
    jl=0;
    ju=n+1;
    ascnd=(xx[n] >= xx[1]);
    while (ju-jl > 1) {
	jm=(ju+jl) >> 1;
	if ((x >= xx[jm]) == ascnd)
	    jl=jm;
	else
	    ju=jm;
    }
    if (x == xx[1]) *j=1;
    else if(x == xx[n]) *j=n-1;
    else *j=jl;
}



void SPHofpos(SPHbody *btab, int nobj)
{
    SPHbody *p;
    double rad,vr;
    FILE *fp;
    double rin, rhoin, uin, pressin, velin;
    double *r, *rho, *u, *press, *vel;
    int i, index;

    fp = fopen("inputmodel.dat", "r");
    i=0;
    while ( !feof(fp) )	{
	fscanf(fp, "%lf %lf %lf %lf %lf\n", &rin, &rhoin, &uin, 
		       &pressin, &velin);
	i++;
    }
    fclose(fp);
    
    r=Malloc(i*sizeof(double));
    rho=Malloc(i*sizeof(double));
    u=Malloc(i*sizeof(double));
    press=Malloc(i*sizeof(double));
    vel=Malloc(i*sizeof(double));
    
    fp = fopen("inputmodel.dat", "r");
    i=0;
    while ( !feof(fp) ) {
	fscanf(fp, "%lf %lf %lf %lf %lf\n", &rin, &rhoin, &uin, 
		       &pressin, &velin);
	r[i]=rin;
	rho[i]=rhoin;
	u[i]=uin;
	press[i]=pressin;
	vel[i]=velin;
	i++;
    }
    fclose(fp);

    for (p=btab; p<btab+nobj; p++) {
	rad=sqrt( p->pos[0]*p->pos[0] + p->pos[1]*p->pos[1] + 
		  p->pos[2]*p->pos[2]);
	if (rad < r[0]) rad=r[0];
	if (rad > r[i-1]) rad=r[i-1];
		     
	locate(r-1, i, rad, &index);   /* REMEMBER! Numerical recipes is */
	index--;	   /* unitordered, i.e. arrays start at 1! */
	
	/*  polint(r-1,u-1,i+1,rad,&uin,&temp); */
	p->u=u[index]+(rad-r[index])/(r[index+1]-r[index])*
	    (u[index+1]-u[index]);
	p->rho=rho[index]+(rad-r[index])/(r[index+1]-r[index])*
	    (rho[index+1]-rho[index]);
	p->pr=press[index]+(rad-r[index])/(r[index+1]-r[index])*
	    (press[index+1]-press[index]);
        vr=vel[index]+(rad-r[index])/(r[index+1]-r[index])*
            (vel[index+1]-vel[index]);
        if (rad > 0.0) {
            p->vel[0]=vr*p->pos[0]/rad;
            p->vel[1]=vr*p->pos[1]/rad;
            p->vel[2]=vr*p->pos[2]/rad;
        }
        else {
            p->vel[0] = p->vel[1] = p->vel[2] = 0.0;
        }
    }
}

