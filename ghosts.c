#include <math.h>
#include <stdlib.h>
#include "stk.h"
#include "Malloc.h"
#include "randoms.h"


#include "physics_sph.h"
#include "Msgs.h"
#include "timers.h"
#include "vop.h"
#include "fastflpt.h"
#include "mpmy.h"
#include "gc.h"
#include "singlio.h"
#include "wvt.h"

#define PI 3.141592653589793238462
#define GHOST_FLAG (1<<28)



void
SPHPlusSPH(void **btabp, int *nobj, SPHbody *SPHbtab, int SPHnobj);

void
SphericalGhosts(SPHbody **btab, int *nobj, double rmax, double rmin, 
		double *tothvol, double totvol, double targetneighbors, 
		int dim)
{
    /* If you want to save memory+computing time when you don't have rmin
       or rmax, set them to something below 0 and over 1e29   */
    
    SPHbody *p, *ghost, *currentghost;
    double r, rhat[NDIM], rghost, h;
    int nghosts=0, maxnghosts=0;
    double q=2.; /* Make sure this is the same as in the wvt.c */
    int do_rmin=0, do_rmax=0;
    double delr, scalfac;
    
    if (rmin > 0.) {
	maxnghosts+=*nobj;
	do_rmin=1;
    }
    if (rmax < 1e29) {
	maxnghosts+=*nobj;
	do_rmax=1;
    }
    singlPrintf("boundaries:%g %g, maxnghosts:%d do_rmin:%d do_rmax:%d\n",
		rmin, rmax, maxnghosts, do_rmin, do_rmax);
    /*   ghost = Malloc(ghostp, maxnghosts * sizeof(SPHbody));  */
    ghost = Malloc(maxnghosts*sizeof(SPHbody));

    /*   ghost = Realloc(*ghostp, maxnghosts * sizeof(SPHbody)); */
    currentghost=ghost;
    singlPrintf("boundaries:%g %g\n",rmin, rmax);
  
    /* avg distance to neighbor with targetneighbors within 4pi/3*(2)^3 */
/*     scalfac=pow(64./targetneighbors,0.333333333) ;  /\* in 3d *\/ */
    scalfac=4*pow(targetneighbors,-1./((double) dim));
    for (p=*btab; p<*btab+*nobj; p++) {
	r=sqrt(p->pos[0]*p->pos[0]+p->pos[1]*p->pos[1]+p->pos[2]*p->pos[2]);
	if (r < 1e-29) r=1e-29;

	delr=p->h*scalfac;
	
	VV(rhat, =1./r*p->pos); /* Unit radial vector */
	h=p->h;
	
	/* Check if the particle is beyond the outer boundary */
	if (r > rmax && do_rmax) {
	    VV(p->pos, =rmax*rhat);
	    r=rmax;
	}
	
	/* Check if the particle is beyond the inner boundary */
	if (r < rmin && do_rmin) {
	    VV(p->pos, =rmin*rhat);
	    r=rmin;
	}
	
	/* Add ghost particles if particles are close to the outer boundary */
	if (r+q*h > rmax && do_rmax && r != rmax) {	
	    rghost=2*rmax-r+0.5*delr;
	    VV(currentghost->pos, =rghost*rhat);
	    WVT_hofpos(currentghost,1,totvol,tothvol, dim);
	    currentghost->rho=p->rho;
	    currentghost->h=p->h;
	    currentghost->u=p->u;
	    currentghost->pr=p->pr;
	    currentghost->mass=p->mass;
	    currentghost->ident=nghosts | GHOST_FLAG;
	    currentghost++;
	    nghosts++;
	}
	
	/* Add ghost particles if particles are close to the inner boundary */
	if (r-q*h < rmin && r-q*h>0 && do_rmin && r !=rmin) {	
	    rghost=2*rmax-r;
	    VV(currentghost->pos, =rghost*rhat);
	    WVT_hofpos(currentghost,1,totvol,tothvol, dim);
	    currentghost->rho=p->rho;
	    currentghost->h=p->h;
	    currentghost->u=p->u;
	    currentghost->pr=p->pr;
	    currentghost->mass=p->mass;
	    currentghost->ident=nghosts | GHOST_FLAG;
	    currentghost++;
	    nghosts++;
	}
	
	
	/* Add ghosts close to the outer boundary if they sit on top of it */
	if (r == rmax && do_rmax) {	
	    rghost=rmax+0.5*delr;
	    VV(currentghost->pos, =rghost*rhat);
	    WVT_hofpos(currentghost,1,totvol,tothvol, dim);
	    currentghost->rho=p->rho;
	    currentghost->h=p->h;
	    currentghost->u=p->u;
	    currentghost->pr=p->pr;
	    currentghost->mass=p->mass;
	    currentghost->ident=nghosts | GHOST_FLAG;
	    currentghost++;
	    nghosts++;
	}
	
	/* Add ghosts close to the inner boundary if they sit on top of it */
	if (r == rmin && do_rmin) {	
	    rghost=rmin-1e-20;
	    VV(currentghost->pos, =rghost*rhat);
	    WVT_hofpos(currentghost,1,totvol,tothvol, dim);
	    currentghost->rho=p->rho;
	    currentghost->h=p->h;
	    currentghost->u=p->u;
	    currentghost->pr=p->pr;
	    currentghost->mass=p->mass;
	    currentghost->ident=nghosts | GHOST_FLAG;
	    currentghost++;
	    nghosts++;
	}
    }


    singlPrintf("nghosts=%d ", nghosts);
    SPHPlusSPH((void **)btab, nobj, ghost, nghosts);
    Free(ghost);
}



void
RemoveGhosts(SPHbody **btabp, int *nobj)
{
    SPHbody *btab = *btabp;
    SPHbody *p, *next;

    /* Shrink btab, taking out SPH particles and copying them to atab */
    for (p = next = btab; p < btab+*nobj; p++) {
	if (!(p->ident & GHOST_FLAG)) *next++ = *p;
    }
    *nobj = next-btab;
    *btabp = Realloc(btab, *nobj * sizeof(SPHbody));
}


void
RemoveGhosts2(SPHbody **btabp, int *nobj)
{
    SPHbody *btab = *btabp;
    SPHbody *p, *next, *good;
/*     Stk s; */
/*     SPHbody *q; */
    int ngood=0;

    /*     StkInitEz(&s); */

    next=btab;
    good=Malloc(*nobj*sizeof(SPHbody));
    next=good;

    for (p=btab; p<btab+*nobj; p++)
      {
	if (!(p->ident & GHOST_FLAG)) { 
	  p->h=3.3333333;
	  *good=*p;
	  good++;
	  ngood++;
	}
      }


/*     /\* Shrink btab, taking out SPH particles and copying them to atab *\/ */
/*     for (p = next = btab; p < btab+*nobj; p++) { */
/* /\*       singlPrintf("%d ", p->ident); *\/ */
/* 	if (p->ident & GHOST_FLAG) { */
/* 	    q = StkPush(&s, sizeof(SPHbody)); */
/* 	    VV(q->grav_acc, = p->acc); */
/* 	    q->phi = p->phi; */
/* 	    q->grav_nterms = p->nterms; */
/* /\* 	    q->ident = p->ident & ~GHOST_FLAG; *\/ */
/* 	    q->key = p->key; */
/* 	} else *next++ = *p; */
/*     } */
/*     StkCrunch(&s); */
/*     *nobj = next-btab; */
/*     *btabp = Realloc(btab, *nobj * sizeof(SPHbody)); */

    *btabp=next;
	*nobj=ngood;
/* 	Free(good); */
}







