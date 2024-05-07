#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "stk.h"
#include "Malloc.h"
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
#include "initial.h"

#define PI 3.141592653589793238462
#ifndef M_1_PI
#define	M_1_PI 0.31830988618379067154
#endif
/*#define NKERNEL_TABLE 40000*/
#define NKERNEL_TABLE 40000
#define MAX_INDEX (NKERNEL_TABLE+2)

#define NO_UPDATE 2
#define GHOST_FLAG (1<<28)

Counter_t SPHCnt, SPHrej, nbrMACCnt;

static double dvtable; /* == 0.0001 ... */
static double invdvtable; /* == 10000.0 ... */
static double cnormk;
static double wij[MAX_INDEX];
static double grwij[MAX_INDEX];
static double fmass[MAX_INDEX];
static double fpoten[MAX_INDEX];
static double Gamma = (double)(4.0/3.0);  /* Yikes; remember this */
static double alpha = (double)1.0;
static double beta = (double)2.5;
static double epsil = (double)1e-2;
static double heatf1 = (double)1.0;
static int ndim;
static int Nobj;
static int add_offset;
static double offset[NDIM];
static double voffset[NDIM];
static void (*bodyfunc)(SinkSPH *sink, hcell **src_vec, int *res, int n);
static void (*cellfunc)(SinkSPH *sink, hcell **src_vec, int *res, int n);
static double *rglob, *hglob;
static int iglob;

int inputoption=1;



/* Cylindrical Grid global variables (Initialized in InitCylGrid)*/
double ***cylgrid;
double ***cylgrid_h;
double minr=0.0039525692, maxr=1., minz=-1, maxz=1;
double mintheta=0., maxtheta=6.2586416;
int dimr=127, dimz=49, dimtheta=256;
double cylcenter[3];

/* Cartesian Grid global variables (Initialized in InitCartGrid*/
double ***cartgrid;
double ***cartgrid_h;
double minx=-1, maxx=1., miny=-1, maxy=1;
/* double minz=0., maxz=6.2586416; */
int dimx=127, dimy=49;/*  dimz=256 */
double cartcenter[3];


extern int do_diffusion;


void
SetWVTOffset(double *off, double *voff)
{
    VV(offset, = off);
    VV(voffset, = voff);
    add_offset = 1;
}

void
UnSetWVTOffset(void)
{
    VS(offset, = 0.0);
    VS(voffset, = 0.0);
    add_offset = 0;
}

void 
InheritWVT(const SinkSPH *from, SinkSPH *to, hcell *pp)
{
    if( to == NULL ){
	SPHbody *bp = pp->ptr;
	if (from->isbody == NO_UPDATE)
	  return;
	/* Must accumulate for periodic BC to work */
	/* Must initialize to zero appropriately */
	bp->rho = from->rho;
	bp->du += from->du;       /* Why am I doing this? */
	bp->du_r += from->du_r;   /* Do I need to do anything else? */
	bp->drho_dt += from->drho_dt;
	bp->udot += from->udot;
	bp->nbrs += from->nbrs;
	VV(bp->acc, += from->M1);
	VV(bp->lvel, += from->lvel);
	bp->nterms += from->nterms;
	bp->min_nbr_dt = from->min_nbr_dt;   /* set dt to min of nbrs dt */
	if (from->interactions != Nobj)
	    Error("Ninteract is %d, should be %d\n", from->interactions, Nobj);
	return;
    }

    if (Sub_Flags(pp)) {
	/* Stuff needed for cell-cell neighbor evaluation */
	SPHcell *cp = pp->ptr;
	VV(to->pos, = cp->pos);
	to->extent = cp->bmax + cp->lap;
	to->isbody = 0;
    } else {
	/* Stuff needed for above, plus physics info */
	SPHbody *bp = pp->ptr;
	if (!SPH_need_update(bp)) {
	    to->isbody = NO_UPDATE;
	    return;
	}
	VV(to->pos, = bp->pos);
	to->extent = bp->h;
	to->h = bp->h;
	to->isbody = 1;
	VV(to->vel, = bp->vel);
	to->pr = bp->pr;
	to->rho_est = (double)1e30;  /* Store resolution for WVT */
	to->mass = bp->mass;
	to->vsound = bp->vsound;
	to->u = bp->u;
	to->rho = (double)1e30; /* Store resolution for WVT */
	to->drho_dt = (double)1e30; /* Store resolution for WVT */
	to->udot = (double)1e30; /* Store resolution for WVT */
	VS(to->lvel, = (double)0.0);
	to->nterms = 1;
	to->nbrs = 0;
	VS(to->M1, = (double)0.0);
	to->min_nbr_dt = 1e30;
	to->temp = (double)1e30; /* Store resolution for WVT */
	to->du = bp->du;
	to->u_r = bp->u_r;
	to->du_r = bp->du_r;  /* Or = (double)0.0; ? */
/* 	to->D = bp->D; */
    }

    if (add_offset) {
	VV(to->pos, += offset);
	VV(to->vel, += voffset);
    }

    if (from) {
	to->interactions = from->interactions;
    } else {
	to->interactions = 0;
    }

}


void WVTInitCube(SPHbody **btabp, int *gnobj, int *nobj, double min[NDIM], 
		 double max[NDIM], int num[NDIM], int dim)
{
    Stk s;
    SPHbody *q;
    double pos[NDIM], delta[NDIM];
    int start = 0, i = 0;
    float randoffset[NDIM];
    ran_state ranstate;
    int seed=1;

    delta[0] = (max[0]-min[0])/((double)(num[0]-0.9));
    delta[1] = (max[1]-min[1])/((double)(num[1]-0.9));
    delta[2] = (max[2]-min[2])/((double)(num[2]-0.9));

    *gnobj = ((int)ceil((max[0]-min[0])/delta[0])) * 
	((int)ceil((max[1]-min[1])/delta[1])) * 
	((int)ceil((max[2]-min[2])/delta[2]));
    
    NobjInitial(*gnobj, MPMY_Nproc(), MPMY_Procnum(), nobj, &start);

    (*btabp) = Malloc((*nobj)*sizeof(SPHbody));
    
    StkInitEz(&s);
    ran_init(seed*(MPMY_Procnum()+1), &ranstate);

    for(pos[0] = min[0]; pos[0] <= max[0]+delta[0]/3.0; 
	pos[0] += delta[0]) {

	for(pos[1] = min[1]; pos[1] <= max[1]+delta[1]/3.0; 
	    pos[1] += delta[1]) {

	    for(pos[2] = min[2]; pos[2] <= max[2]+delta[2]/3.0; 
		pos[2] += delta[2]) {
	      
		if ( (i >= start) && (i < start + *nobj) ) {
		    q = StkPush(&s, sizeof(SPHbody));
		    /* Take the 1e-20 out to "randomize" the grid */
 		    VVVV(q->pos, =randoffset,*1e-20*delta, + pos);
		    q->rho = 1.0; /* Can be anything but 0...*/
		    q->u = 1.0; /* Can be anything but 0...*/
		    q->mass=1.0; /* Can be anything but 0...*/
		    VS(q->vel,=0.); /* Can be anything but 0...*/
		    q->h=0.;
		    q->ident=i;
		    /* q->nterms = 1; */
		}
		i++;
	    }

	}

    }

    StkCrunch(&s);
    *nobj = StkSz(&s)/sizeof(SPHbody);
    q = StkBase(&s);
    *btabp = Realloc(q, *nobj * sizeof(SPHbody));
    
    Msgf(("gnobj = %d; nobj = %d; start = %d\n", *gnobj, *nobj, start));
}


void WVTInitHex(SPHbody **btabp, int *gnobj, int *nobj, double min[NDIM], 
		double max[NDIM], int num[NDIM], int dim)
{
    /* INCORRECT VERSION !!! USE WVTINITHCP INSTEAD */
    Stk s;
    SPHbody *q;
    double pos[NDIM], delta[NDIM];
    int start = 0, i = 0;
    float randoffset[NDIM];
    ran_state ranstate;
    int seed=1;
    double xoff, xmin, yoff, ymin, zmin;

    delta[0] = (max[0]-min[0])/((double)(num[0]-1.));
    delta[1] = delta[0]*0.5*sqrt(3.0);
    delta[2] = delta[0]*sqrt(6.0)/3.0;

    zmin= (int)trunc( min[2]/delta[2] )*delta[2];

    *gnobj = ((int)ceil((max[0]-min[0])/delta[0])) * 
	((int)ceil((max[1]-min[1])/delta[1])) * 
	((int)ceil((max[2]-min[2])/delta[2]));
    
    NobjInitial(*gnobj, MPMY_Nproc(), MPMY_Procnum(), nobj, &start);

    (*btabp) = Malloc((*nobj)*sizeof(SPHbody));
    
    StkInitEz(&s);
    ran_init(seed*(MPMY_Procnum()+1), &ranstate);

    for(pos[2] = zmin; pos[2] <= max[2]; pos[2] += delta[2]) {

        yoff = ((int)round(fabs(pos[2])/delta[2]) % 2) * delta[1] 
	    * (1.0 - 1.0/sqrt(3.0));
	ymin=min[1]+yoff;

	for(pos[1] = ymin; pos[1] <= max[1]; pos[1] += delta[1]) {

            xoff = ((int)round((pos[1] + yoff)/delta[1]) % 2) * 0.5*delta[0];
            xmin = min[0] + xoff;

	    for(pos[0] = xmin; pos[0] <= max[0]; pos[0] += delta[0]) {
	      
		/* Make sure I am really in my box */
		if ( pos[0] >= min[0] && pos[0] <=max[0] &&
		     pos[1] >= min[1] && pos[1] <=max[1] &&
		     pos[2] >= min[2] && pos[2] <=max[2] ) {
		    if ( (i >= start) && (i < start + *nobj) ) {
			q = StkPush(&s, sizeof(SPHbody));
			/* Take the 1e-20 out to "randomize" the grid */
			VVVV(q->pos, =randoffset,*1e-20*delta, + pos);
			q->rho = 1.0; /* Can be anything but 0...*/
			q->u = 1.0; /* Can be anything but 0...*/
			q->mass=1.0; /* Can be anything but 0...*/
			VS(q->vel,=0.); /* Can be anything but 0...*/
			q->h=0.;
			q->ident=i;
			/* q->nterms = 1; */
		    }
		    i++;
		}
	    }

	}

    }

    StkCrunch(&s);
    *nobj = StkSz(&s)/sizeof(SPHbody);
    q = StkBase(&s);
    *btabp = Realloc(q, *nobj * sizeof(SPHbody));
    
    Msgf(("gnobj = %d; nobj = %d; start = %d\n", *gnobj, *nobj, start));
}


void WVTInitCCP(SPHbody **btabp, int *gnobj, int *nobj, double min[NDIM], 
		double max[NDIM], int num[NDIM], int dim)
{
    Stk s;
    SPHbody *q;
    double pos[NDIM], delta[NDIM];
    int start = 0, i = 0;
    float randoffset[NDIM];
    ran_state ranstate;
    int seed=1;
    double xoff, xmin, yoff, ymin, zmin;
    int zlayer, ylayer;

    delta[0] = (max[0]-min[0])/((double)(num[0]-1.));
    delta[1] = delta[0]*0.5*sqrt(3.0);
    delta[2] = delta[0]*sqrt(6.0)/3.0;

    zmin= (int)trunc( min[2]/delta[2] )*delta[2];

    *gnobj = ((int)ceil((max[0]-min[0])/delta[0])) * 
	((int)ceil((max[1]-min[1])/delta[1])) * 
	((int)ceil((max[2]-min[2])/delta[2]));
    
    NobjInitial(*gnobj, MPMY_Nproc(), MPMY_Procnum(), nobj, &start);

    (*btabp) = Malloc((*nobj)*sizeof(SPHbody));
    
    StkInitEz(&s);
    ran_init(seed*(MPMY_Procnum()+1), &ranstate);

    zlayer=0; /* 0=A, 1=B, 2=C, 3=A, 4=B, ... */
    for(pos[2] = zmin; pos[2] <= max[2]; pos[2] += delta[2]) {
	zlayer=zlayer % 3;
	singlPrintf("zlayer: %d\n", zlayer);

/*         yoff = ((double) zlayer) * delta[1]  */
/* 	    * (1.0 - 1.0/sqrt(3.0)); */

	if (zlayer == 0) yoff = 0.;
       	if (zlayer == 1) yoff = 2./3*delta[1];
       	if (zlayer == 2) yoff = -2./3*delta[1];

/*        	if (zlayer == 1) yoff = delta[1] * (1.0 - 1.0/sqrt(3.0)); */
/* 	if (zlayer == 2) yoff = -delta[1] * (1.0 - 1.0/sqrt(3.0));	 */

	ymin=min[1]+yoff;

	ylayer=0;
       	for(pos[1] = ymin; pos[1] <= max[1]; pos[1] += delta[1]) {
	    ylayer=ylayer % 2;
	    
	    xoff=0.;
	    if (zlayer == 0) xoff += 0.;
	    if (zlayer == 1) xoff += 0.;
	    if (zlayer == 2) xoff += 0.;
	
	    xoff += ((double) ylayer) * 0.5*delta[0]; 
/*             xoff = ((double) ylayer) * 0.5*delta[0]; */
            xmin = min[0] + xoff;

	    for(pos[0] = xmin; pos[0] <= max[0]; pos[0] += delta[0]) {
	      
		/* Make sure I am really in my box */
		if ( pos[0] >= min[0] && pos[0] <=max[0] &&
		     pos[1] >= min[1] && pos[1] <=max[1] &&
		     pos[2] >= min[2] && pos[2] <=max[2] ) {
		    if ( (i >= start) && (i < start + *nobj) ) {
			q = StkPush(&s, sizeof(SPHbody));
			/* Take the 1e-20 out to "randomize" the grid */
			VVVV(q->pos, =randoffset,*1e-20*delta, + pos);
			q->rho = 1.0; /* Can be anything but 0...*/
			q->u = 1.0; /* Can be anything but 0...*/
			q->mass=1.0; /* Can be anything but 0...*/
			VS(q->vel,=0.); /* Can be anything but 0...*/
			q->h=0.;
			q->ident=i;
			/* q->nterms = 1; */
		    }
		    i++;
		}
	    }

	    ylayer++;
	}

	zlayer++;
    }

    StkCrunch(&s);
    *nobj = StkSz(&s)/sizeof(SPHbody);
    q = StkBase(&s);
    *btabp = Realloc(q, *nobj * sizeof(SPHbody));
    
    Msgf(("gnobj = %d; nobj = %d; start = %d\n", *gnobj, *nobj, start));
}

void WVTInitHCP(SPHbody **btabp, int *gnobj, int *nobj, double min[NDIM], 
		double max[NDIM], int num[NDIM], int dim)
{
    Stk s;
    SPHbody *q;
    double pos[NDIM], delta[NDIM];
    int start = 0, i = 0;
    float randoffset[NDIM];
    ran_state ranstate;
    int seed=1;
    double xoff, xmin, yoff, ymin, zmin;
    int zlayer, ylayer;

    delta[0] = (max[0]-min[0])/((double)(num[0]-1.));
    delta[1] = delta[0]*0.5*sqrt(3.0);
    delta[2] = delta[0]*sqrt(6.0)/3.0;

    zmin= (int)trunc( min[2]/delta[2] )*delta[2];

    *gnobj = ((int)ceil((max[0]-min[0])/delta[0])) * 
	((int)ceil((max[1]-min[1])/delta[1])) * 
	((int)ceil((max[2]-min[2])/delta[2]));
    
    NobjInitial(*gnobj, MPMY_Nproc(), MPMY_Procnum(), nobj, &start);

    (*btabp) = Malloc((*nobj)*sizeof(SPHbody));
    
    StkInitEz(&s);
    ran_init(seed*(MPMY_Procnum()+1), &ranstate);

    zlayer=0; /* 0=A, 1=B, 2=C, 3=A, 4=B, ... */
    for(pos[2] = zmin; pos[2] <= max[2]; pos[2] += delta[2]) {
	zlayer=zlayer % 2;
	singlPrintf("zlayer: %d\n", zlayer);

/*         yoff = ((double) zlayer) * delta[1]  */
/* 	    * (1.0 - 1.0/sqrt(3.0)); */

	if (zlayer == 0) yoff = 0.;
       	if (zlayer == 1) yoff = 2./3*delta[1];
       	if (zlayer == 2) yoff = -2./3*delta[1];

/*        	if (zlayer == 1) yoff = delta[1] * (1.0 - 1.0/sqrt(3.0)); */
/* 	if (zlayer == 2) yoff = -delta[1] * (1.0 - 1.0/sqrt(3.0));	 */

	ymin=min[1]+yoff;

	ylayer=0;
       	for(pos[1] = ymin; pos[1] <= max[1]; pos[1] += delta[1]) {
	    ylayer=ylayer % 2;
	    
	    xoff=0.;
	    if (zlayer == 0) xoff += 0.;
	    if (zlayer == 1) xoff += 0.;
	    if (zlayer == 2) xoff += 0.;
	
	    xoff += ((double) ylayer) * 0.5*delta[0]; 
/*             xoff = ((double) ylayer) * 0.5*delta[0]; */
            xmin = min[0] + xoff;

	    for(pos[0] = xmin; pos[0] <= max[0]; pos[0] += delta[0]) {
	      
		/* Make sure I am really in my box */
		if ( pos[0] >= min[0] && pos[0] <=max[0] &&
		     pos[1] >= min[1] && pos[1] <=max[1] &&
		     pos[2] >= min[2] && pos[2] <=max[2] ) {
		    if ( (i >= start) && (i < start + *nobj) ) {
			q = StkPush(&s, sizeof(SPHbody));
			/* Take the 1e-20 out to "randomize" the grid */
			VVVV(q->pos, =randoffset,*1e-20*delta, + pos);
			q->rho = 1.0; /* Can be anything but 0...*/
			q->u = 1.0; /* Can be anything but 0...*/
			q->mass=1.0; /* Can be anything but 0...*/
			VS(q->vel,=0.); /* Can be anything but 0...*/
			q->h=0.;
			q->ident=i;
			/* q->nterms = 1; */
		    }
		    i++;
		}
	    }

	    ylayer++;
	}

	zlayer++;
    }

    StkCrunch(&s);
    *nobj = StkSz(&s)/sizeof(SPHbody);
    q = StkBase(&s);
    *btabp = Realloc(q, *nobj * sizeof(SPHbody));
    
    Msgf(("gnobj = %d; nobj = %d; start = %d\n", *gnobj, *nobj, start));
}


void WVTInitHex2(SPHbody **btabp, int *gnobj, int *nobj, double min[NDIM], 
		 double max[NDIM], int num[NDIM], int dim)
{
    /* INCORRECT VERSION !!! USE WVTINITHCP INSTEAD */
    double pos[NDIM];
    int start = 0, i = 0;
    float randoffset[NDIM];
    ran_state ranstate;
    int seed=1;

    SPHbody *p = NULL;
    double dz, dy, dx, dr;
    double z, y, x, h, r_outer;
    double yoff, xoff;
    double xmin, xmax, ymin, ymax, zmin, zmax;

    singlPrintf("xmin=%g, ymin=%g, zmin=%g \n", min[0], min[1], min[2]);
    singlPrintf("xmax=%g, ymax=%g, zmax=%g \n", max[0], max[1], max[2]);    

    xmin=min[0];
    ymin=min[1];
    zmin=min[2];

    xmax=max[0];
    ymax=max[1];
    zmax=max[2];

    singlPrintf("xmin=%g, ymin=%g, zmin=%g \n", xmin, ymin, zmin);
    singlPrintf("xmax=%g, ymax=%g, zmax=%g \n", xmax, ymax, zmax);    

    dr=(xmax-xmin)/((double)(num[0]-0.9));
    dz = dr*sqrt(6.0)/3.0;
    dy = dr*0.5*sqrt(3.0);
    dx = dr;
    r_outer=xmax;
    singlPrintf("dx=%g, dy=%g, dz=%g \n", dx, dy, dz);

    *gnobj = ((int)ceil((max[0]-min[0])/dx)) * 
	((int)ceil((max[1]-min[1])/dy)) * 
	((int)ceil((max[2]-min[2])/dz));
    
    NobjInitial(*gnobj, MPMY_Nproc(), MPMY_Procnum(), nobj, &start);

    (*btabp) = Malloc((*nobj)*sizeof(SPHbody));
    ran_init(seed*(MPMY_Procnum()+1), &ranstate);

    for (z = zmin; z <= zmax; z += dz) {
        h = sqrt(r_outer*r_outer - pos[2]*pos[2]);

        yoff = ((int)round(fabs(z)/dz) % 2) * dy * (1.0 - 1.0/sqrt(3.0));
        ymin = (int)trunc( (-h - yoff)/dy )*dy + yoff;

        if (ymin < -h) continue;

        for (y = ymin; y <= ymax; y += dy) {
            xoff = ((int)round((y + yoff)/dy) % 2) * 0.5*dx;
            xmin = (int)trunc( (-sqrt(h*h - y*y) - xoff)/dx )*dx + xoff;
            
            if (xmin < -sqrt(h*h - y*y)) continue;

            for (x = xmin; x <= sqrt(h*h - y*y); x += dx) {
		if ( (i >= start) && (i < start + *nobj) ) {
		    p = (SPHbody *)realloc(p, ++(*nobj) * sizeof(SPHbody));
		    /* Take the 1e-20 out to "randomize" the grid */
		    /* You need this factor for the tree not to fail due to 
		     having "perfect planes" */
		    p[*nobj-1].pos[0] = x+randoffset[0]*1e-20*dx;
		    p[*nobj-1].pos[1] = y+randoffset[1]*1e-20*dy;
		    p[*nobj-1].pos[2] = z+randoffset[2]*1e-20*dz;
		    p[*nobj-1].rho = 1;		    
		    p[*nobj-1].u = 1;		    
		    p[*nobj-1].mass = 1;		    
		    p[*nobj-1].h = 1;		    
		    p[*nobj-1].ident = i;		    
		    VS(p[*nobj-1].vel, = 0.);		    
		}
	    i++;	    
	    singlPrintf("x=%g, y=%g, z=%g \n", x, y, z);
	    }
	}

    }

    *btabp = p;
    singlPrintf("gnobj = %d; nobj = %d; start = %d\n", *gnobj, *nobj, start);
}


/* From Gabe */
void hexcp(double dr, double r_outer, SPHbody **btabp, int *nobj)
{
    /* INCORRECT VERSION !!! USE WVTINITHCP INSTEAD */
    SPHbody *p = NULL;
    double x, y, z, h;
    double dz = dr*sqrt(6.0)/3.0, dy = dr*0.5*sqrt(3.0), dx = dr;
    double yoff, ymin;
    double xoff, xmin;
    double zmin = (int)trunc( -r_outer/dz )*dz;

    *nobj = 0;

    for (z = zmin; z <= r_outer; z += dz) {
        h = sqrt(r_outer*r_outer - z*z);

        yoff = ((int)round(fabs(z)/dz) % 2) * dy * (1.0 - 1.0/sqrt(3.0));
	/* that factor should be 2/3 instead of 1-1./sqrt(3) */

        ymin = (int)trunc( (-h - yoff)/dy )*dy + yoff;

        if (ymin < -h) continue;

        for (y = ymin; y <= h; y += dy) {
            xoff = ((int)round((y + yoff)/dy) % 2) * 0.5*dx;
            xmin = (int)trunc( (-sqrt(h*h - y*y) - xoff)/dx )*dx + xoff;
            
            if (xmin < -sqrt(h*h - y*y)) continue;

            for (x = xmin; x <= sqrt(h*h - y*y); x += dx) {
                p = (SPHbody *)realloc(p, ++(*nobj) * sizeof(SPHbody));
                p[*nobj-1].pos[0] = x;
                p[*nobj-1].pos[1] = y;
                p[*nobj-1].pos[2] = z;
            }
        }
    }

    *btabp = p;
}


void WVTInitProbdist(SPHbody **btabp, int *gnobj, int *nobj, double min[NDIM], 
		     double max[NDIM], int targetnobj, double totvol, 
		     double outerbound, double innerbound, int num[NDIM], 
		     int dim)
{
    Stk s, ss;
    SPHbody *q, *qq;
    double pos[NDIM], delta[NDIM];
    int start = 0, i = 0, j=0;
    float randoffset[NDIM];
    ran_state ranstate;
    int seed=1;
    double h, rho, r;
    int nadded;
    double tothvol=1.;
    double subtotrho, totrho=0.;

    if (num[0] < 0) num[0]=1024;
    if (num[1] < 0) num[1]=1024;
    if (num[2] < 0) num[2]=1024;

    delta[0] = (max[0]-min[0])/((double)(num[0]));
    delta[1] = (max[1]-min[1])/((double)(num[1]));
    delta[2] = (max[2]-min[2])/((double)(num[2]));

    *gnobj=targetnobj;
    NobjInitial(*gnobj, MPMY_Nproc(), MPMY_Procnum(), nobj, &start);

    (*btabp) = Malloc((*nobj)*sizeof(SPHbody));
    q=Malloc(sizeof(SPHbody));

    StkInitEz(&s);
    ran_init(seed*(MPMY_Procnum()+1), &ranstate);
    
    i=0;
    for(pos[0] = min[0]+delta[0]/2.; pos[0] <= max[0]; 
	pos[0] += delta[0]) {
	
	for(pos[1] = min[1]+delta[1]/2.; pos[1] <= max[1]; 
	    pos[1] += delta[1]) {
	    
	    for(pos[2] = min[2]+delta[2]/2.; pos[2] <= max[2]; 
		pos[2] += delta[2]) {
		
		r=sqrt(pos[0]*pos[0]+pos[1]*pos[1]+pos[2]*pos[2]); 
		if (r < outerbound && r > innerbound) 
		    {
			VV(q->pos, =pos);
			WVT_hofpos(q,1,totvol, &tothvol, dim);
			h=q->h;
			rho=pow(h,-dim);
			totrho+=rho;
			i++;
		    }
	    }
	    
	}
	
    }

    StkInitEz(&ss);
    i=start;   
    subtotrho=0.;
    for(pos[0] = min[0]+delta[0]/2.; pos[0] <= max[0]; 
	pos[0] += delta[0]) {
	
	for(pos[1] = min[1]+delta[1]/2.; pos[1] <= max[1]; 
	    pos[1] += delta[1]) {
	    
	    for(pos[2] = min[2]+delta[2]/2.; pos[2] <= max[2]; 
		pos[2] += delta[2]) {
		
		r=sqrt(pos[0]*pos[0]+pos[1]*pos[1]+pos[2]*pos[2]); 
		if (r < outerbound && r > innerbound) {
		    VV(q->pos, =pos);
		    WVT_hofpos(q,1,totvol, &tothvol, dim);
		    h=q->h;
		    rho=pow(h,-dim);
		    subtotrho+=rho;
		    rho*=((double) *gnobj)/totrho;
		    nadded=ignpoi((float)rho);
		    
		    if ( subtotrho/totrho >=
			 (double)MPMY_Procnum()/((double) MPMY_Nproc()) &&
			 subtotrho/totrho <
			 ((double)MPMY_Procnum()+1)/((double) MPMY_Nproc())) {
			for (j=0; j<nadded; j++) {
			    /* 	    VVVV(pos, =randoffset,*delta, + p->pos); */
			    /* 	    VV(pos, =p->pos); */
			    VVVV(q->pos, =randoffset,*delta,/2. + pos); 
			    r=sqrt(q->pos[0]*q->pos[0]+q->pos[1]*q->pos[1]+q->pos[2]*q->pos[2]);  
			    /* 	    r=1;  */
			    
			    if ( r < outerbound && r > innerbound ) {
				qq = StkPush(&ss, sizeof(SPHbody));
				VV(qq->pos, =q->pos); 
				qq->rho = 1.0; /* anything but 0...*/
				qq->u = 1.0; /* anything but 0...*/
				qq->mass=1.0; /* anything but 0...*/
				VS(qq->vel,=0.); /* anything but 0...*/
				qq->h=0.;
				qq->ident=i;
				/* q->nterms = 1; */
				i++;
			    }
			}
		    }
		}
	    }
	}
    }
    singlPrintf("totrho:%g subtotrho:%g\n", totrho, subtotrho);
    /*    MPMY_Combine(&totrho, &totrho, 1, MPMY_DOUBLE, MPMY_SUM);*/

    StkCrunch(&ss);
    *nobj = StkSz(&ss)/sizeof(SPHbody);
    qq = StkBase(&ss);
    *btabp = Realloc(qq, *nobj * sizeof(SPHbody));
    Free(q);
    Msgf(("gnobj = %d; nobj = %d; start = %d\n", *gnobj, *nobj, start));
}


void WVTInitProbdistlr(SPHbody **btabp, int *gnobj, int *nobj, double min[NDIM], 
		     double max[NDIM], int targetnobj, double totvol, 
		     double outerbound, double innerbound, int num[NDIM], 
		     int dim)
{
    Stk s, ss;
    SPHbody *q, *qq;
    double delta[NDIM];
    int start = 0, i = 0, j=0;
    ran_state ranstate;
    int seed=1;
    double h, rho, r, rad, r1, r2;
    int nadded;
    double tothvol=1.;
    double subtotrho, totrho=0.;
    double testran, testran2;
    int id;

    if (num[0] < 0) num[0]=1024;
    if (num[1] < 0) num[1]=1024;
    if (num[2] < 0) num[2]=1024;

    delta[0] = (log(outerbound)-log(innerbound))/((double)(num[0]));

    *gnobj=targetnobj;
    NobjInitial(*gnobj, MPMY_Nproc(), MPMY_Procnum(), nobj, &start);

    (*btabp) = Malloc((*nobj)*sizeof(SPHbody));
    q=Malloc(sizeof(SPHbody));

    StkInitEz(&s);
    ran_init(seed*(MPMY_Procnum()+1), &ranstate);
    
    for(i=0; i < num[0]; ++i) {
      r = exp(log(innerbound)+((double)i+0.5)*delta[0]);
      if (r < outerbound && r > innerbound) 
	{
	  q->pos[0] = r;
	  q->pos[1] = 0.0;
	  q->pos[2] = 0.0;
	  WVT_hofpos(q,1,totvol, &tothvol, dim);
	  h=q->h;
	  rho=pow(h,-dim);
          r2=exp(log(r)+0.5*delta[0]);
          r1=exp(log(r)-0.5*delta[0]);
	  totrho+=rho*((r2*r2*r2-r1*r1*r1)/(outerbound*outerbound*outerbound));
	}
	
    }

    StkInitEz(&ss);
    id=start;   
    subtotrho=0.;
    for(i=0; i < num[0]; ++i) {
      r = exp(log(innerbound)+((double)i+0.5)*delta[0]);
      if (r < outerbound && r > innerbound) 
	{
	  q->pos[0] = r;
	  q->pos[1] = 0.0;
	  q->pos[2] = 0.0;
	  WVT_hofpos(q,1,totvol, &tothvol, dim);
	  h=q->h;
	  rho=pow(h,-dim);
          r2=exp(log(r)+0.5*delta[0]);
          r1=exp(log(r)-0.5*delta[0]);
	  subtotrho+=rho*((r2*r2*r2-r1*r1*r1)/(outerbound*outerbound*outerbound));
	  rho*=((double) *gnobj)/totrho*((r2*r2*r2-r1*r1*r1)/(outerbound*outerbound*outerbound));
	  nadded=ignpoi((float)rho);
	  
	  if ( ( subtotrho/totrho >=
                 (double)MPMY_Procnum()/((double) MPMY_Nproc()) &&
                 subtotrho/totrho <
                 ((double)MPMY_Procnum()+1)/((double) MPMY_Nproc()) ) ||
               ( i == num[0]-1 &&
                 (MPMY_Procnum()+1) == MPMY_Nproc() ) ) {
	    for (j=0; j<nadded; j++) {
	      testran=uniform_rand(&ranstate)-0.5;
	      rad=exp(log(r)+testran*delta[0]);
	      testran=2.0*(uniform_rand(&ranstate)-0.5);
	      testran2=2.0*PI*(uniform_rand(&ranstate));
	      q->pos[0] = rad*sqrt(1-testran*testran)*cos(testran2);
	      q->pos[1] = rad*sqrt(1-testran*testran)*sin(testran2);
	      q->pos[2] = rad*testran;
              /*WVT_hofpos(q,1,totvol,&tothvol,dim);*/
	      if ( rad < outerbound && rad > innerbound ) {
		qq = StkPush(&ss, sizeof(SPHbody));
		VV(qq->pos, =q->pos); 
		qq->rho = 1.0; /* anything but 0...*/
		qq->u = 1.0; /* anything but 0...*/
		qq->mass=1.0; /* anything but 0...*/
		VS(qq->vel,=0.); /* anything but 0...*/
                qq->h=0.;
		//qq->h=q->h;
		qq->ident=id;
		/* q->nterms = 1; */
		id++;
	      }
	    }
	    Msgf(("r1 = %g, r2 = %g, nadded = %d\n", r1, r2, nadded));
	  }
	}
    }
    singlPrintf("totrho:%g subtotrho:%g\n", totrho, subtotrho);
    /*    MPMY_Combine(&totrho, &totrho, 1, MPMY_DOUBLE, MPMY_SUM);*/

    StkCrunch(&ss);
    *nobj = StkSz(&ss)/sizeof(SPHbody);
    qq = StkBase(&ss);
    *btabp = Realloc(qq, *nobj * sizeof(SPHbody));
    Free(q);
    Msgf(("gnobj = %d; nobj = %d; start = %d\n", *gnobj, *nobj, start));
}



void WVTInitProbdist_lotsmem(SPHbody **btabp, int *gnobj, int *nobj, 
			     double min[NDIM], double max[NDIM], 
			     int targetnobj, double totvol, 
			     double outerbound, double innerbound, int dim)
{
    Stk s, ss;
    SPHbody *q, *p, *qq;
    double pos[NDIM], delta[NDIM];
    int start = 0, i = 0, j=0;
    float randoffset[NDIM];
    ran_state ranstate;
    int seed=1;
    double h, rho, r;
    int nadded;
    int num[NDIM];
    double tothvol=1.;
    double totrho=0.;

    long int nbox;

    num[0]=128;
    num[1]=128;
    num[2]=128;

    delta[0] = (max[0]-min[0])/((double)(num[0]-0.9));
    delta[1] = (max[1]-min[1])/((double)(num[1]-0.9));
    delta[2] = (max[2]-min[2])/((double)(num[2]-0.9));

    *gnobj=targetnobj;
    NobjInitial(*gnobj, MPMY_Nproc(), MPMY_Procnum(), nobj, &start);

    (*btabp) = Malloc((*nobj)*sizeof(SPHbody));

    StkInitEz(&s);
    ran_init(seed*(MPMY_Procnum()+1), &ranstate);

    for(pos[0] = min[0]; pos[0] <= max[0]+delta[0]/3.0; 
	pos[0] += delta[0]) {
	
	for(pos[1] = min[1]; pos[1] <= max[1]+delta[1]/3.0; 
	    pos[1] += delta[1]) {
	    
	    for(pos[2] = min[2]; pos[2] <= max[2]+delta[2]/3.0; 
		pos[2] += delta[2]) {
		
		r=sqrt(pos[0]*pos[0]+pos[1]*pos[1]+pos[2]*pos[2]); 
		if (r < outerbound && r > innerbound) {
		    q = StkPush(&s, sizeof(SPHbody));
		    VV(q->pos, =pos);
		    WVT_hofpos(q,1,totvol, &tothvol, dim);
		    h=q->h;
		    
		    q->u=r; 
		    q->rho=1./(h*h*h);
		    totrho+=q->rho;
		    
		    /* Since we do the probability distribution now, we need 
		       to know the whole grid on each node. */
		    /* 	      q->rho = 1.0; /\* Can be anything but 0...*\/ */
		    /* 	      q->u = 1.0; /\* Can be anything but 0...*\/ */
		    q->mass=totrho; /* Can be anything but 0...*/
		    VS(q->vel,=0.); /* Can be anything but 0...*/
		    /* 	      q->h=0.; */
		    q->ident=i;
		    /* q->nterms = 1; */
		    i++;
		}
	    }
	}
    }

    /*     MPMY_Combine(&totrho, &totrho, 1, MPMY_DOUBLE, MPMY_SUM); */
    StkCrunch(&s);
    nbox = StkSz(&s)/sizeof(SPHbody); 
    q = StkBase(&s);
    /*     *btabp = Realloc(q, *nobj * sizeof(SPHbody)); */
    StkInitEz(&ss);
    
    i=start;
    /*     while (i <= start+*nobj) {*/
    for (p=q; p<q+nbox; p++) 
	{
	    if (p->mass/totrho>=(double)MPMY_Procnum()/((double) MPMY_Nproc()) 
		&& p->mass/totrho<((double)MPMY_Procnum()+1)/
		((double) MPMY_Nproc())) {
		
		rho=p->rho/totrho* ((double) *gnobj);
		/* 	  singlPrintf("%g ", p->rho); */
		nadded=ignpoi((float)rho);
		
		for (j=0; j<nadded; j++) {
		    /*testran=uniform_rand(&ranstate);*/
		    /* 	    VVVV(pos, =randoffset,*delta, + p->pos); */
		    /* 	    VV(pos, =p->pos); */
		    VVVV(pos, =randoffset,*delta,/2. + p->pos); 
		    r=sqrt(pos[0]*pos[0]+pos[1]*pos[1]+pos[2]*pos[2]);  
		    /* 	    r=1;  */
		    
		    if ( r < outerbound && r > innerbound) {
			qq = StkPush(&ss, sizeof(SPHbody));
			VV(qq->pos, =pos); 
			qq->rho = 1.0; /* Can be anything but 0...*/
			qq->u = 1.0; /* Can be anything but 0...*/
			qq->mass=1.0; /* Can be anything but 0...*/
			VS(qq->vel,=0.); /* Can be anything but 0...*/
			qq->h=0.;
			qq->ident=i;
			/* q->nterms = 1; */
			i++;
		    }
		}
	    }
	}	        
    /*        }  */

    StkCrunch(&ss);
    *nobj = StkSz(&ss)/sizeof(SPHbody);
    qq = StkBase(&ss);
    *btabp = Realloc(qq, *nobj * sizeof(SPHbody));
    Free(q);

    Warning("Processor Id:%d Nproc:%d\n", MPMY_Procnum(), MPMY_Nproc() ); 
    Msgf(("gnobj = %d; nobj = %d; start = %d\n", *gnobj, *nobj, start));
}



void WVTInitCube2(SPHbody **btabp, int *gnobj, int *nobj, double min[NDIM], 
		  double max[NDIM], int num[NDIM], int dim)
{
    Stk s;
    SPHbody *q;
    double pos[NDIM], delta[NDIM];
    int start = 0, i = 0;
    float randoffset[NDIM];
    ran_state ranstate;
    int seed=1;
    
    delta[0] = (max[0]-min[0])/((double)(num[0]-0.9));
    delta[1] = (max[1]-min[1])/((double)(num[1]-0.9));
    delta[2] = (max[2]-min[2])/((double)(num[2]-0.9));
    
    *gnobj = ((int)ceil((max[0]-min[0])/delta[0])) * 
	((int)ceil((max[1]-min[1])/delta[1])) * 
	((int)ceil((max[2]-min[2])/delta[2]));
    
    NobjInitial(*gnobj, MPMY_Nproc(), MPMY_Procnum(), nobj, &start);
    
    (*btabp) = Malloc((*nobj)*sizeof(SPHbody));
    
    StkInitEz(&s);
    ran_init(seed*(MPMY_Procnum()+1), &ranstate);
    
    for(pos[0] = min[0]; pos[0] <= max[0]+delta[0]/3.0; 
	pos[0] += delta[0]) {
	
	for(pos[1] = min[1]; pos[1] <= max[1]+delta[1]/3.0; 
	    pos[1] += delta[1]) {
	    
	    for(pos[2] = min[2]; pos[2] <= max[2]+delta[2]/3.0; 
		pos[2] += delta[2]) {
		
		if ( (i >= start) && (i < start + *nobj) ) {
		    q = StkPush(&s, sizeof(SPHbody));
		    /* 		    VV(q->pos, = pos); */
 		    VVVV(q->pos, =randoffset,*delta, + pos);
		    q->rho = 1.0; /* Can be anything but 0...*/
		    q->u = 1.0; /* Can be anything but 0...*/
		    q->mass=1.0; /* Can be anything but 0...*/
		    VS(q->vel,=0.); /* Can be anything but 0...*/
		    q->h=0.;
		    q->ident=i;
		    /* q->nterms = 1; */
		}		
		i++;
	    }
	}
    }
    
    StkCrunch(&s);
    *nobj = StkSz(&s)/sizeof(SPHbody);
    q = StkBase(&s);
    *btabp = Realloc(q, *nobj * sizeof(SPHbody));
    
    Msgf(("gnobj = %d; nobj = %d; start = %d\n", *gnobj, *nobj, start));
}


void WVTInitProbdist2(SPHbody **btabp, int *gnobj, int *nobj, 
		      double min[NDIM], double max[NDIM], int targetnobj, 
		      double totvol, double outerbound, double innerbound, 
		      int dim)
{
    Stk s, ss;
    SPHbody *q, *p, *qq;
    double pos[NDIM], delta[NDIM];
    int start = 0, i = 0;
    float randoffset[NDIM];
    ran_state ranstate;
    int seed=1;
    double h, rho, r;
    int nadded;
    int num[NDIM];
    double tothvol=1.;
    double totrho=0.;
    double testran;
    long int nbox;

    num[0]=64;
    num[1]=64;
    num[2]=64;

    delta[0] = (max[0]-min[0])/((double)(num[0]-0.9));
    delta[1] = (max[1]-min[1])/((double)(num[1]-0.9));
    delta[2] = (max[2]-min[2])/((double)(num[2]-0.9));

    *gnobj=targetnobj;
    NobjInitial(*gnobj, MPMY_Nproc(), MPMY_Procnum(), nobj, &start);

    (*btabp) = Malloc((*nobj)*sizeof(SPHbody));

    StkInitEz(&s);
    ran_init(seed*(MPMY_Procnum()+1), &ranstate);

    for(pos[0] = min[0]; pos[0] <= max[0]+delta[0]/3.0; 
	pos[0] += delta[0]) {

	for(pos[1] = min[1]; pos[1] <= max[1]+delta[1]/3.0; 
	    pos[1] += delta[1]) {

	    for(pos[2] = min[2]; pos[2] <= max[2]+delta[2]/3.0; 
		pos[2] += delta[2]) {
	      
		q = StkPush(&s, sizeof(SPHbody));
		VV(q->pos, =pos);
		WVT_hofpos(q,1,totvol, &tothvol, dim);
		h=q->h;
		r=sqrt(q->pos[0]*q->pos[0]+q->pos[1]*q->pos[1]+
		       q->pos[2]*q->pos[2]); 
		q->u=r; 
		q->rho=1./(h*h*h);
		totrho+=q->rho;
		/* Since we do the probability distribution now, we need to 
		   know the whole grid on each node. */
		q->mass=1.0; /* Can be anything but 0...*/
		VS(q->vel,=0.); /* Can be anything but 0...*/
		/* 	      q->h=0.; */
		q->ident=i;
		/* q->nterms = 1; */
		i++;
	    }
	}
    }

    MPMY_Combine(&totrho, &totrho, 1, MPMY_DOUBLE, MPMY_SUM);

    StkCrunch(&s);
    nbox = StkSz(&s)/sizeof(SPHbody); 

    q = StkBase(&s);
    StkInitEz(&ss);

    i=start;
    /*     while (i <= start+*nobj)  */
    /*        {  */
    for (p=q; p<q+nbox; p++) {
	rho=p->rho/totrho* ((double) *nobj);
	/* 	  singlPrintf("%g ", p->rho); */
	nadded=0;
	testran=uniform_rand(&ranstate);
	/* 	    VVVV(pos, =randoffset,*delta, + p->pos); */
	/* 	    VV(pos, =p->pos); */
	r=sqrt(p->pos[0]*p->pos[0]+p->pos[1]*p->pos[1]+p->pos[2]*p->pos[2]);  
	/* 	    r=1;  */
	if ( testran < rho && 
	     r < outerbound && r > innerbound) 
	    {
		qq = StkPush(&ss, sizeof(SPHbody));
		VVVV(qq->pos, =randoffset,*delta, + p->pos); 
		qq->rho = 1.0; /* Can be anything but 0...*/
		qq->u = 1.0; /* Can be anything but 0...*/
		qq->mass=1.0; /* Can be anything but 0...*/
		VS(qq->vel,=0.); /* Can be anything but 0...*/
		qq->h=0.;
		qq->ident=i;
		/* q->nterms = 1; */
		i++;
		nadded++;
	    }
	
    }
    /*       }  */


    StkCrunch(&ss);
    *nobj = StkSz(&ss)/sizeof(SPHbody);
    qq = StkBase(&ss);
    *btabp = Realloc(qq, *nobj * sizeof(SPHbody));



    Msgf(("gnobj = %d; nobj = %d; start = %d\n", *gnobj, *nobj, start));
}




void
WVTgate(SinkSPH *sink, hcell **src_vec, int *flags, int *result, int n)
{
    int i;
    
    if (sink->isbody == NO_UPDATE) {
	for (i = 0; i < n; i++) result[i] = MAC_ACCEPT;
    } else if (sink->isbody)
	bodyfunc(sink, src_vec, result, n);
    else
	cellfunc(sink, src_vec, result, n);
}

void
macWVT(SinkSPH *sink, hcell **source_vec, int *result, int n)
{
    VxdV(const double pos_sink, = sink->pos);
    Vxd(double r);
    Vxd(double f);
    double min_nbr_dt = sink->min_nbr_dt;
    double extent_src;
    int daughters;
    int i;
    SPHbody *bp = 0;
    int nbrs = 0;
    double dr2;
    Vxd(double runi);
    double rij;
    int interactions = 0;
    double deltai, deltaj;
    double eps=1e-2*sink->h; /* REMEMBER THIS !!! */
    double norm;
    double fudge=1.; /* Fudge factor, should be somewhat close to 1 */
    double q=2.; /* Should definitely be >= 2! Extent to which neighbors are
		    enclosed */
    double res1=sink->udot, res2=sink->rho_est, res3=sink->vsound, 
	res4=sink->temp; 
    double hmean11, hmean21, v2, h;
    int index;
    /* Store the resolution in some other variables you don't need now */
    
    /* In the following comments, sink is j, sources are i's */
    
    deltaj=1./pow(sink->h,fudge); 
    VxS(f, = (double)0.0); /* Initialize "force" to 0 */
    for (i = 0; i < n; i++) {
	const hcellptr source = source_vec[i];
	
	if (Sub_Flags(source)) {
	    const SPHcell *cp = source->ptr;
	    VxV(r, = cp->pos);	
	    extent_src = cp->bmax + cp->lap;
	    daughters = cp->daughters;
	} else {
	    bp = source->ptr;
	    VxV(r, = bp->pos);
	    extent_src = bp->h;
	    daughters = 1;
	}
	
	VxVx(r, -= pos_sink);	/* r012 are the comp. of rj->ri vector now */
	dr2 = Dotx(r, r);
	
	sink->nterms += 1;
	
	/* rij= distance ri to rj */
	rij = sqrtf_fast(dr2);
	h=(extent_src+sink->h)/2.;
	if (rij > q*h
	    || dr2 == (double) 0.) {
	    goto accept; /* skip if too far away */
	} else if (daughters != 1) {
	    goto failed;
	}
	
	if (bp->dt_next < min_nbr_dt) min_nbr_dt = bp->dt_next;
	
	VxVx(runi, = 1./rij*r);  /*Unity vector from rj to ri */
	deltai=1./pow(bp->h, fudge);

 	hmean11 = (double)2.0 / (sink->h + bp->h); 
 	hmean21 = hmean11 * hmean11; 
	
 	v2 = dr2 * hmean21; 
 	index = v2 * invdvtable; 
 	if (index >= MAX_INDEX) Error("Index too large\n"); 

/*  	norm=-deltaj/deltai * wtij;  */
/*  	norm=-deltaj/deltai/((rij+eps)*(rij+eps));  */
  	norm=-deltaj/deltai*(1./(v2+eps)-1./(4.+eps));  
	VxVx(f, += norm * runi); /* f is one piece of the sum that 
				    makes up the "acceleration term */
	
	/* Now find the resolution by determining the distance to the 
	   4 closest neighbors */
	if (rij < res1) {
	    res4=res3;
	    res3=res2;
	    res2=res1;
	    res1=rij;
	} else if (rij < res2) {
	    res4=res3;
	    res3=res2;
	    res2=rij;
	} else if (rij < res3) {
	    res4=res3;
	    res3=rij;
	} else if (rij < res4) {
	    res4=rij;
	}

	sink->udot=res1;
	sink->rho_est=res2;
	sink->vsound=res3;
	sink->temp=res3;
	sink->drho_dt=(res1+res2+res3+res4)/4.;

	sink->rho=(res1+res2+res3+res4)/4.;

	
	nbrs++;
    accept:
	IncrCounter(&SPHrej);
	interactions += daughters;
	result[i] = MAC_ACCEPT;
	continue;
    failed:
	IncrCounter(&SPHrej);
	result[i] = MAC_SPLIT_SRC;
    }
    sink->interactions += interactions;
    VVx(sink->M1, += f);
    sink->nbrs += nbrs;
    
}



void WVTupdate(SPHbody *btab, int nobj, int loop, int nloop, int dim, 
	       int nneighbors)
{
    SPHbody *p;
    double dt, delx,dely,delz, h;
    double q=0.5, oneoverdim=1./((double) dim);

    for (p=btab; p<btab+nobj; p++) {
	h=p->h/pow( (double) nneighbors, oneoverdim);
	/* h/N^.33 is the actual separation */
 	if (dim == 3) dt=1./((double) nneighbors)*10*h/
			  (1.+0*3*( (double) loop/( (double) nloop))); 
 	else /* dim == 2 */
            dt=h*h*h/32./(1.+2*( (double) loop/( (double) nloop))); 
/*  	dt=h*h*h/10.;  */
	if (!(p->ident & GHOST_FLAG)) {
	    delx=dt*p->acc[0];
	    dely=dt*p->acc[1];
	    delz=dt*p->acc[2];
	    
	    if (abs(delx) <= q*h) p->pos[0]+=delx; 
	    else if (delx > q*h) p->pos[0]+=q*h;
	    else if (delx < -q*h) p->pos[0]+=-q*h;
	    
	    if (abs(dely) <= q*h) p->pos[1]+=dely; 
	    else if (dely > q*h) p->pos[1]+=q*h;
	    else if (dely < -q*h) p->pos[1]+=-q*h;
	    
	    if (abs(delz) <= q*h) p->pos[2]+=delz; 
	    else if (delz > q*h) p->pos[2]+=q*h;
	    else if (delz < -q*h) p->pos[2]+=-q*h;
	}
    }  
}

void WVT_setinputoption(int option)
{
    inputoption=option;
}

void WVT_hofpos(SPHbody *btab, int nobj, double totvol, double *tothvol, 
		int dim)
{
    
    if (inputoption == 1) {
	WVT_hofpos_inputh(btab,nobj,totvol, tothvol, dim);
    } else   if (inputoption == 2) {
	WVT_hofpos_cylgrid(btab,nobj,totvol, tothvol, dim);
    } else   if (inputoption == 3) {
	WVT_hofpos_pwl(btab,nobj,totvol, tothvol, dim);
    } else   if (inputoption == 4) {
	WVT_hofpos_cartgrid(btab,nobj,totvol, tothvol, dim);
    }     
}

void WVT_hofpos_pwl(SPHbody *btab, int nobj, double totvol, double *tothvol,
		    int dim)
{
    SPHbody *p;
    double rad;
    double q=1.;
    int redo_tothvol=0;
    if (*tothvol < 0) redo_tothvol=1;
    
    if (redo_tothvol)  *tothvol=0.;
    
    for (p=btab; p<btab+nobj; p++) {
	rad=sqrt( p->pos[0]*p->pos[0] + p->pos[1]*p->pos[1] + 
		  p->pos[2]*p->pos[2]);
	
	/* Maybe I should use input parameters here. Right now this is 
	   hardwired. */
	p->h=pow(rad,1.25)+0.75; 
	if (redo_tothvol) *tothvol+=p->h*p->h*p->h;
    }
    
    if (redo_tothvol)  *tothvol *=4.*sqrt(3.);
    
    if (redo_tothvol) MPMY_Combine(tothvol, tothvol, 1, MPMY_DOUBLE, MPMY_SUM);
    
    for (p=btab; p<btab+nobj; p++) {
	p->h*=q*pow( (totvol/ *tothvol), 0.333333333 );
    }
}



void
SetWVT(double visc_alpha, double visc_beta, double visc_epsilon, 
       double heat_f1, double eos_gamma, int gnobj,  void bfunc(), 
       void cfunc())
{
    Nobj = gnobj;
    alpha = visc_alpha;
    beta = visc_beta;
    epsil = visc_epsilon;
    heatf1 = heat_f1;
    Gamma = eos_gamma;
    bodyfunc = bfunc;
    cellfunc = cfunc;    
}


void
WVT_setup(int dim, int ncoef1, double *wcoef1, int ncoef2, double *wcoef2)
{
    double v2max;
    double v, v2;
    double w, dw, dm, dm1;
    double ddvtable;
    int i, i1, j;
    
    
    ndim = dim;
    
    /* maximum interaction length and step size */
    v2max = (double)4.0;
    ddvtable = v2max / NKERNEL_TABLE;
    /*i1 = (double)1.0 / dvtable;*/
    invdvtable = (double)1.0 / ddvtable;
    dvtable = ddvtable;
    
    /* normalization constant */
    if (ndim == 3)
	cnormk = (double)M_1_PI; /* 1./pi */
    else if (ndim == 2)
	cnormk = (double)(M_1_PI*10.0/7.0);
    else if (ndim == 1)
	cnormk = (double)(2.0/3.0);
    else
	Error("Bad ndim in sph_ktable\n");

    /* build tables */
    /* a) v less than 1 */
    wij[0] = cnormk*wcoef1[0];
    grwij[0] = 0.0;     
    i1 = 1.0 / ddvtable;

    for (i = 1; i <= i1; i++) {
	v2 = i * ddvtable;
	v = sqrt(v2);
	w = wcoef1[ncoef1-1];
	for (j = ncoef1-1; j >= 1 ; j--)
	    w = w * v + wcoef1[j-1];
	dw = (ncoef1 - 1.0) * wcoef1[ncoef1-1];
	for (j = ncoef1-2; j >= 1; j--)
	    dw = dw * v + j * wcoef1[j];
	wij[i] = cnormk * w;
	grwij[i] = cnormk * dw;

	/* Enclosed mass now for m=1, h=1*/
	dm = wcoef1[ncoef1-1]/(ncoef1-1+3);
	for (j = ncoef1-2; j >= 0; j--)
	    dm = dm * v + wcoef1[j]/(j+3.);
	fmass[i]=4*PI*cnormk*v*v*v*dm;

	/* Potential now: fpoten is potential for G=1, m=1, h=1*/
	fpoten[i]=(double) 0.;
	for (j = 0;j <= ncoef1-1; j++) 
	    fpoten[i]+=wcoef1[j]/(j+3.)*pow(v,(j+2.))+
		wcoef1[j]/(j+2.)*(1-pow(v,(j+2.)))+
		wcoef2[j]/(j+2.)*(pow(2.,(j+2.))-1.);
	fpoten[i]*=4.*PI*cnormk;
    }
    
    v=1.;
    dm1=wcoef2[ncoef2-1]/(ncoef2-1+3.);
    for (j=ncoef2-2; j>=0;j--)  
	dm1=dm1*v + wcoef2[j]/(j+3.);      
    dm1=fmass[i1]-4*PI*dm1*cnormk*v*v*v;
    

     /*  b) v greater than 1 */
    for (i = i1+1; i <= NKERNEL_TABLE; i++) {
	v2 = i * ddvtable;
	v = sqrt(v2);
	w = wcoef2[ncoef2-1];
	for (j = ncoef2-1; j >= 1 ; j--)
	    w = w * v + wcoef2[j-1];
	dw = (ncoef2 - 1.0) * wcoef2[ncoef2-1];
	for (j = ncoef2-2; j >= 1; j--)
	    dw = dw * v + j * wcoef2[j];
	wij[i] = cnormk * w;
	grwij[i] = cnormk * dw;
	
	/* Enclosed mass now */
	dm = wcoef2[ncoef2-1]/(ncoef2-1+3);
	for (j = ncoef2-2; j >= 0; j--)
	    dm = dm * v + wcoef2[j]/(j+3.);
	fmass[i]=4*PI*cnormk*v*v*v*dm+dm1;

	/* Potential now: fpoten is potential for G=1, m=1, h=1*/
	fpoten[i]=(double) 0.;
	for (j=0; j <= ncoef2-1; j++)
	    fpoten[i]+=1./v*wcoef1[j]/(j+3.)+
		1./v*wcoef2[j]/(j+3.)*(pow(v,(j+3.))-1.)+
		wcoef2[j]/(j+2.)*(pow(2.,(j+2.))-pow(v,(j+2.)));
	fpoten[i]*=4.*PI*cnormk;
    }
    
    /* Make sure the mass is perfectly normalized, don't generate mass */
    for (i = 0; i <= NKERNEL_TABLE; i++) {
	fmass[i]/=fmass[NKERNEL_TABLE];
    }

    /* For interpolation of last table entry */
    wij[NKERNEL_TABLE+1] = grwij[NKERNEL_TABLE+1] = (double) 0.0;
    fpoten[NKERNEL_TABLE+1]=(double) 0.5;
    fmass[NKERNEL_TABLE+1]=(double) 1.;
}

void
update_WVT(SPHbody *btab, int nobj, double dt, int *limit_high, int *limit_low)
{
    SPHbody *p;
    
    for (p = btab; p < btab+nobj; p++) {
	if (!SPH_need_update(p)) continue;
/* 	VV(p->acc, += p->grav_acc); */
	p->rho += cnormk * p->mass / (p->h * p->h * p->h);
	p->hdot = (double)(-1.0/3.0) * p->h * p->drho_dt / p->rho;
	if (p->hdot * dt > p->h) {
	    SeriousWarning("Hdot limit (high)\n%s\n", PrintSPHBodyContents(p));
	    p->hdot = p->h/dt;
	    ++*limit_high;
	}
	if (p->hdot * dt < -0.5*p->h) {
	    SeriousWarning("Hdot limit (low)\n%s\n", PrintSPHBodyContents(p));
	    p->hdot = -0.5*p->h/dt;
	    ++*limit_low;
	}
	
	p->udot += p->drho_dt * p->pr / (p->rho * p->rho);
	
	if (!finite(p->udot)) 
	    Error("Bad value for udot\n");
	
	/* Are these limits appropriate? */
	/* Does this enforce the Courant limit correctly with diffusion? */
	if (p->udot * dt > p->u) {
	    p->udot = p->u/dt;
	    ++*limit_high;
	}
	if (p->udot * dt < -0.333*p->u) {
	    p->udot = -0.333*p->u/dt;
	    ++*limit_low;
	}
    }
}



void
macConstNeigh(SinkSPH *sink, hcell **source_vec, int *result, int n)
{
    const double extent_sink = sink->extent;
    VxdV(const double pos_sink, = sink->pos);
    const double h = sink->h;
    Vxd(double r);
    Vxd(double dv);
    double extent_src;
    int daughters;
    int i;
    SPHbody *bp = 0;
    double projv;
    double v2;
    double rij;
    double wtij, grwtij;
    double hmean11, hmean21;
    double dxx, dwdx, dgrwdx;
    int index;
    int nbrs = 0;
    double rhoi = (double)0.0;
    double divvi = (double)0.0;
    double dr2;
    int interactions = 0;
    
    for (i = 0; i < n; i++) {
	const hcellptr source = source_vec[i];
	
	if (Sub_Flags(source)) {
	    const SPHcell *cp = source->ptr;
	    VxV(r, = cp->pos);	
	    extent_src = cp->bmax + cp->lap;
	    daughters = cp->daughters;
	} else {
	    bp = source->ptr;
	    VxV(r, = bp->pos);
	    extent_src = bp->h;
	    daughters = 1;
	}
	
	VxVx(r, -= pos_sink);	/* 8 flops */
	dr2 = Dotx(r, r); /* == 400.0000000032623 */
	
	if ((rij = sqrtf_fast(dr2)) > extent_sink+extent_src
	    || dr2 == (double)0.0) {
	    goto accept;
	} else if (daughters != 1) {
	    goto failed;
	}
	
 	hmean11 = (double)2.0 / (h + bp->h); 
 	hmean21 = hmean11 * hmean11; 
	
 	v2 = dr2 * hmean21; 
 	index = v2 * invdvtable; 
 	if (index >= MAX_INDEX) Error("Index too large\n"); 
 	dxx = v2 - index * dvtable; 
 	dwdx = (wij[index+1] - wij[index]) * invdvtable; 
 	wtij = (wij[index] + dwdx * dxx ) * hmean21 * hmean11; 
 	if (wtij < (double)0.0) Error("Negative wtij (macRho) = %g\n", wtij); 
 	dgrwdx = (grwij[index+1] - grwij[index]) * invdvtable; 
 	grwtij = (grwij[index] + dgrwdx * dxx) * hmean21 * hmean21; 
	
 	rhoi += bp->mass * wtij; 

 	/* velocity divergence times density */
 	VxVV(dv, = bp->vel, - sink->vel); 
 	projv = grwtij * Dotx(dv, r) * recipsqrtf(dr2); 
 	divvi -= bp->mass * projv; 
	
	nbrs++;
    accept:
	interactions += daughters;
	result[i] = MAC_ACCEPT;
	continue;
    failed:
	result[i] = MAC_SPLIT_SRC;
    }
    
    sink->interactions += interactions; 
    sink->rho += rhoi; 
    sink->nbrs += nbrs;
    sink->drho_dt -= divvi; 
}


void WVT_hofpos_inputh(SPHbody *btab, int nobj, double totvol, 
		       double *tothvol, int dim)
{
    SPHbody *p;
    double rad;
    double q=0.5;
    FILE *fp;
    double rin, hin;
    /*   double *r, *h; */
    int i, index;
    int redo_tothvol=0;
    
    if (*tothvol < 0) redo_tothvol=1;
    if (redo_tothvol)  *tothvol=0.;
    
    
    /*IF nobj=-1 initialize rglob and hglob only */
    if (nobj == -1 ) {
	
	fp = fopen("inputh.dat", "r");
	i=0;
	while ( !feof(fp) ) {
	    fscanf(fp, "%lf %lf\n", &rin, &hin);
	    i++;
	}
	fclose(fp);
	
	rglob=Malloc(i*sizeof(double));
	hglob=Malloc(i*sizeof(double));
	
	fp = fopen("inputh.dat", "r");
	i=0;
	while ( !feof(fp) ) {
	    fscanf(fp, "%lf %lf\n", &rin, &hin);
	    rglob[i]=rin;
	    hglob[i]=hin;
	    i++;
	}
	iglob=i;
	fclose(fp);
    } 
    else {	
	i=iglob;
	for (p=btab; p<btab+nobj; p++) {
	    rad=sqrt( p->pos[0]*p->pos[0] + p->pos[1]*p->pos[1] + 
		      p->pos[2]*p->pos[2]);
	    if (rad < rglob[0]) rad=rglob[0];
	    if (rad > rglob[i-1]) rad=rglob[i-1];
	    
	    locate(rglob-1, i, rad, &index);/* REMEMBER! Numerical recipes */
	    index--;	   /* is unitordered, i.e. arrays start at 1! */
	    
	    p->h=hglob[index]+(rad-rglob[index])/(rglob[index+1]-rglob[index])*
		(hglob[index+1]-hglob[index]);
	    
	    if (redo_tothvol) *tothvol+=pow(p->h,dim);
	}
	if (redo_tothvol) 
	    MPMY_Combine(tothvol, tothvol, 1, MPMY_DOUBLE, MPMY_SUM);
	
	for (p=btab; p<btab+nobj; p++) {
	    p->h*=q*pow( (totvol/ *tothvol), 1./dim);
	}
    }
}




void WVT_hofpos_cylgrid(SPHbody *btab, int nobj, double totvol, 
			double *tothvol, int dim)
{
    SPHbody *p;
    double q=1.;
    int redo_tothvol=0;
    double bgrho=1e-5;
    
    if (nobj < 0) singlPrintf("Using Cylindrical Grid data.\n");
    
    if (*tothvol < 0) redo_tothvol=1;
    if (redo_tothvol)  *tothvol=0.;
    
    /* What we need now is the function that loads the grid data*/

    
    interp_cylindricalgrid(btab, nobj, dimr, dimz, dimtheta, minr, maxr,
			   minz, maxz, mintheta, maxtheta, bgrho);
    
    
    if (redo_tothvol) {
	for (p=btab; p<btab+nobj; p++) {
	    *tothvol+=p->h*p->h*p->h;
	}
	*tothvol *=4.*sqrt(3.);
	MPMY_Combine(tothvol, tothvol, 1, MPMY_DOUBLE, MPMY_SUM);
    }
    
    if (nobj > 0) {
	for (p=btab; p<btab+nobj; p++) {
	    p->h*=q*pow( (totvol/ *tothvol), 0.333333333 );
	}
    }
}


double ***allocate3d(int xdim, int ydim, int zdim)
{
    int i, j;
    double ***a3d = (double ***)malloc(xdim * sizeof(double **));
    for(i = 0; i < xdim; i++) {
	a3d[i] = (double **)malloc(ydim * sizeof(double *));
	for(j = 0; j < ydim; j++)
	    a3d[i][j] = (double *)malloc(zdim * sizeof(double));
    }
    return a3d;
}

void free3d(double ***a3d, int xdim, int ydim, int zdim)
{
    int i, j;
    for(i = 0; i < xdim; i++) {
	for(j = 0; j < ydim; j++)
	    free(a3d[i][j]);
	free(a3d[i]);
    }
}




void interp_cylindricalgrid(SPHbody *SPHbtab, int nobj, int dimr, int dimz, 
			    int dimtheta, double minr, double maxr, 
			    double minz, double maxz, 
			    double mintheta, double maxtheta, double bgrho)
{
    int j, k, l, tj, tk, tl;
    double jj, kk, ll;
    double r, z, theta, dr, dz, dtheta;
    double w, sumw, sumprop, sumprop_h;
    SPHbody *q;
    
    
    dr=(double) (maxr-minr)/(dimr-1.);
    dz=(double) (maxz-minz)/(dimz-1.);
    dtheta=(double) (maxtheta-mintheta)/(dimtheta-1.);
    

    if (nobj != -1) {
	for (q=SPHbtab; q<SPHbtab+nobj; q++) {
	    z=q->pos[2]+cylcenter[2];
	    r=sqrt( pow(q->pos[0]+cylcenter[0], 2)+pow(q->pos[1]+
						       cylcenter[1],2) );
	    theta=atan2pi(q->pos[0]+cylcenter[0], q->pos[1]+cylcenter[1]);
	    
	    ll=(double) (theta-mintheta)/dtheta;
	    kk=(double) (fabs(z)-minz)/dz;
	    jj=(double) (r-minr)/dr;
	
	    l=(int) ll;
	    k=(int) kk;
	    j=(int) jj;
	    
	    ll=(double) ll-l;
	    kk=(double) kk-k;
	    jj=(double) jj-j;
	    
	    /* The 8 closest neighbors are the ones with l,k,j and l+1,k+1,j+1
	       indeces. */
	    sumw=0.;
	    sumprop=0.;
	    sumprop_h=0.;
	    for (tl=0; tl<2; tl++)	{
		for (tk=0; tk<2; tk++)	{
		    for (tj=0; tj<2; tj++) {
			if (k+tk >= dimz || j+tj >= dimr) {
			    w=0.;
			} else {
			    /* weight is the inverse squared of the distance */
			    w=(double) 1./(pow(ll-tl,2)+pow(kk-tk,2)+
					   pow(jj-tj,2)+1e-10);
			    sumprop+=w*cylgrid[j+tj][k+tk][(l+tl) % dimtheta];
			    sumprop_h+=w*cylgrid_h[j+tj][k+tk][(l+tl) 
							       % dimtheta]; 
			    sumw+=w;	    
			}
		    }
		}
	    }
	    if (sumw == 0.) {
		q->rho=0.;
		q->h=0.1; /* STANDARDIZE THIS!!!*/
	    } else {
		q->rho=sumprop/sumw;
		q->h=sumprop_h/sumw;
	    }	
	}
    }

    /* Not good programming, but hey, who cares for now?*/ 
    /*   free(cylgrid); */

}


double atan2pi(double x, double y)
     /* Returns the arctangens of y/x, from 0 to 2pi, taking signs of x and
	y into account */
{
    double theta=atan(y/x);
    double pi=3.14159265;
    
    if (x >= 0 && y >= 0) {
	return theta;
    } else if (x < 0 && y >= 0) {
	return pi+theta;
    } else if (x < 0 && y < 0) {
	return pi+theta;
    } else if (x >= 0 && y < 0) {
	return 2*pi+theta;
    }    
    return theta;
}


void init_cylindricalgrid(int dimr, int dimz, 
			  int dimtheta, double minr, double maxr, 
			  double minz, double maxz, 
			  double mintheta, double maxtheta, double center[3])
{
    int j, k, l;
    double ddata, dr, dz, dtheta;
    FILE *infile, *infile2;
    
    dr=(double) (maxr-minr)/(dimr-1.);
    dz=(double) (maxz-minz)/(dimz-1.);
    dtheta=(double) (maxtheta-mintheta)/(dimtheta-1.);
    
    /* Set global variables */
    cylcenter[0]=center[0];
    cylcenter[1]=center[1];
    cylcenter[2]=center[2];
    
    singlPrintf("Interpolating Grid\n");
    cylgrid=allocate3d(dimr, dimz, dimtheta);
    cylgrid_h=allocate3d(dimr, dimz, dimtheta);
    /* remember to free again! */
    
    singlPrintf("Dim: %d %d %d, \n", dimr, dimz, dimtheta);
    singlPrintf("minr %g maxr %g minz %g maxz %g mintheta %g maxtheta %g\n", 
		minr, maxr, minz, maxz, mintheta, maxtheta);
    singlPrintf("dr: %g dz: %g dtheta %g, \n", dr, dz, dtheta);
    
    
    /* cylindrical grid from Patrick Motl */
    /*   infile=fopen("dens_Q0.4","r"); */
    /*   float junk; */
    /*   fread(&junk, sizeof(junk),1,infile); */
    infile=fopen("dens_Q0.4_cutout.cdat","r");
    infile2=fopen("outh_Q0.4_cutout.cdat","r");
    
    /* When you read unformatted fortran data, it is padded by 4 bytes at 
       beginning and end. Remember, data is little endian */
    
    for (l=0; l<dimtheta; l++) {
	for (k=0; k<dimz; k++) {
	    for (j=0; j<dimr; j++) {
		fread(&ddata, sizeof(ddata),1,infile); 
		cylgrid[j][k][l]=ddata;
		fread(&ddata, sizeof(ddata),1,infile2); 
		cylgrid_h[j][k][l]=ddata;
	    }      
	}    
    }
    fclose(infile);	
    fclose(infile2);	
    
    singlPrintf("%g, %g, %g, %g\n", cylgrid[0][0][0], cylgrid[63][0][0], 
		cylgrid[63][0][128], cylgrid[37][7][137]);    
}


void init_cartesiangrid(int tdimx, int tdimy, int tdimz, double tminx, 
			double tmaxx, double tminy, double tmaxy, 
			double tminz, double tmaxz, double center[3], 
			char *cartfile_rho, char *cartfile_h)
{
    int j, k, l;
    double ddata, dx, dy, dz;
    FILE *infile, *infile2;
    
    minx=tminx;
    maxx=tmaxx;
    miny=tminy;
    maxy=tmaxy;
    minz=tminz;
    maxz=tmaxz;
    
    dimx=tdimx;
    dimy=tdimy;
    dimz=tdimz;
    
    dx=(double) (maxx-minx)/(dimx-1.);
    dy=(double) (maxy-miny)/(dimy-1.);
    dz=(double) (maxz-minz)/(dimz-1.);
    
    singlPrintf("spacing: %g %g %g\n", dx, dy, dz);
    
    /* Set global variables */
    cartcenter[0]=center[0];
    cartcenter[1]=center[1];
    cartcenter[2]=center[2];
    
    singlPrintf("Interpolating Cartesian Grid\n");
    cartgrid=allocate3d(dimx, dimy, dimz);
    cartgrid_h=allocate3d(dimx, dimy, dimz);
    /* remember to free again! */
    
    singlPrintf("Dim: %d %d %d, \n", dimx, dimy, dimz);
    
    /* cylindrical grid from Patrick Motl */
    /*   infile=fopen("dens_Q0.4","r"); */
    /*   float junk; */
    /*   fread(&junk, sizeof(junk),1,infile); */
    infile=fopen(cartfile_rho,"r");
    infile2=fopen(cartfile_h,"r");
    
    /* When you read unformatted fortran data, it is padded by 4 bytes at 
       beginning and end. Remember, data is little endian */
    
    for (l=0; l<dimz; l++) {
	for (k=0; k<dimy; k++) {
	    for (j=0; j<dimx; j++) {
		fread(&ddata, sizeof(ddata),1,infile); 
		cartgrid[j][k][l]=ddata;
		fread(&ddata, sizeof(ddata),1,infile2); 
		cartgrid_h[j][k][l]=ddata;
/* 		singlPrintf("%g ", cartgrid_h[j][k][l]); */
	    }      
	}    
    }
    fclose(infile);	
    fclose(infile2);	
    
}



void interp_cartesiangrid(SPHbody *SPHbtab, int nobj, double bgrho)
{
    int j, k, l, tj, tk, tl;
    double jj, kk, ll;
    double x, y, z, dx, dy, dz;
    double w, sumw, sumprop, sumprop_h;
    SPHbody *q;
    
    
    dx=(double) (maxx-minx)/(dimx-1.);
    dy=(double) (maxy-miny)/(dimy-1.);
    dz=(double) (maxz-minz)/(dimz-1.);
    
    /*   singlPrintf("spacing: %g %g %g", dx, dy, dz); */
    
    if (nobj != -1) {
	for (q=SPHbtab; q<SPHbtab+nobj; q++) {
	    x=q->pos[0]+cartcenter[0];
	    y=q->pos[1]+cartcenter[1];
	    z=q->pos[2]+cartcenter[2];
	    
	    ll=(double) (z-minz)/dz;
	    kk=(double) (y-miny)/dy;
	    jj=(double) (x-minx)/dx;
	    
	    l=(int) ll;
	    k=(int) kk;
	    j=(int) jj;
	    
	    ll=(double) ll-l;
	    kk=(double) kk-k;
	    jj=(double) jj-j;

	    /* Make sure you are treating the boundary */
	    if (l > dimz) {
		l=dimz;
		ll=0.;
	    }
	    if (l < 0) {
		l=0;
		ll=0.;
	    }
	    if (k > dimy) {
		k=dimy;
		kk=0.;
	    }
	    if (k < 0) {
		k=0;
		kk=0.;
	    }
	    if (j > dimx) {
		j=dimz;
		jj=0.;
	    }
	    if (j < 0) {
		j=0;
		jj=0.;
	    }

	    /* The 8 closest neighbors are the ones with l,k,j and l+1,k+1,j+1
	       indeces. */
	    sumw=0.;
	    sumprop=0.;
	    sumprop_h=0.;
	    for (tl=0; tl<2; tl++)	{
		for (tk=0; tk<2; tk++)	{
		    for (tj=0; tj<2; tj++) {
			if ( j+tj < 0 || 
			     j+tj >= dimx || 
			     k+tk < 0 || 
			     k+tk >= dimy 
			     || l+tl < 0 || 
			     l+tl >= dimz) {
			    w=0.;
			} else {
			    /* weight is the inverse squared of the distance */
			    w=(double) 1./(pow(ll-tl,2)+pow(kk-tk,2)+
					   pow(jj-tj,2)+1e-10);
			    sumprop+=w*cartgrid[j+tj][k+tk][(l+tl)];
			    sumprop_h+=w*cartgrid_h[j+tj][k+tk][(l+tl)]; 
			    sumw+=w;	    
			}
		    }
		}
	    }
	    if (sumw == 0.) {
		q->rho=0.;
		q->h=1.; /* STANDARDIZE THIS!!!*/
	    } else {
		q->rho=sumprop/sumw;
		q->h=sumprop_h/sumw;
	    }	
	}
    }

    /* Not good programming, but hey, who cares for now?*/ 
    /*   free(cylgrid); */

}


void WVT_hofpos_cartgrid(SPHbody *btab, int nobj, double totvol, 
			 double *tothvol, int dim)
{
    SPHbody *p;
    double q=1.;
    int redo_tothvol=0;
    double bgrho=1e-5;

    if (nobj < 0) singlPrintf("Using Cartesian Grid data.\n");
    
    if (*tothvol < 0) redo_tothvol=1;
    if (redo_tothvol)  *tothvol=0.;
    
    /* What we need now is the function that loads the grid data*/

    
    interp_cartesiangrid(btab, nobj, bgrho);
    
    
    if (redo_tothvol) {
	for (p=btab; p<btab+nobj; p++) {
	    *tothvol+=p->h*p->h*p->h;
	}
	*tothvol *=4.*sqrt(3.);
	MPMY_Combine(tothvol, tothvol, 1, MPMY_DOUBLE, MPMY_SUM);
    }
    
    if (nobj > 0) {
	for (p=btab; p<btab+nobj; p++) {
	    p->h*=q*pow( (totvol/ *tothvol), 0.333333333 );
	}
    }
}
