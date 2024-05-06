#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <strings.h>
#include <math.h>
#include <assert.h>
#include "fastflpt.h"
#include "Msgs.h"
#include "physics.h"
#include "physics_sph.h"
#include "stk.h"
#include "vop.h"
#include "singlio.h"
#include "error.h"
#include "mpmy.h"
#include "units.h"

void hunt(winddata_t *w, int wnobj, float t, float *v, float *mdot, float *u);


void
ShrinkBtab (SPHbody **SPHbtabp, body *btabp, int *nobj, float r_limit)
{
    SPHbody *btab = *SPHbtabp;
    SPHbody *p;
    Stk s;
    SPHbody *q;
    float r2;

    StkInitEz(&s);
    r2 = r_limit*r_limit;

    for (p = btab; p < btab+*nobj; p++) {
      if (Dot(p->pos, p->pos) >= r2) { /* acceptable */
	q = StkPush(&s, sizeof(SPHbody));
	*q = *p;
      } else {
	btabp->accmass += p->mass;
	btabp->l[0] += p->pos[1]*p->vel[2] - p->pos[2]*p->vel[1];
	btabp->l[1] += p->pos[2]*p->vel[0] - p->pos[0]*p->vel[2];
	btabp->l[2] += p->pos[0]*p->vel[1] - p->pos[1]*p->vel[0];
	Msgf(("Point mass gobbled m = %e; total = %e\nAccreted ang momentum = (%e, %e, %e)\n", 
	      p->mass, btabp->accmass, 
	      btabp->l[0], btabp->l[1], btabp->l[2]));
      }
    }
    Free(btab);
    StkCrunch(&s);
    *nobj = StkSz(&s)/sizeof(SPHbody);
    btab = StkBase(&s);
    *SPHbtabp = Realloc(btab, *nobj * sizeof(SPHbody));
}

void
ShrinkBtab2 (SPHbody **SPHbtabp, int *nobj, float r_limit)
{
    SPHbody *btab = *SPHbtabp;
    SPHbody *p;
    Stk s;
    SPHbody *q;
    float r2;

    StkInitEz(&s);
    r2 = r_limit*r_limit;

    for (p = btab; p < btab+*nobj; p++) {
      if (Dot(p->pos, p->pos) >= r2) { /* acceptable */
	q = StkPush(&s, sizeof(SPHbody));
	*q = *p;
      } 
/*        else { */
/*  	btabp->accmass += p->mass; */
/*  	btabp->l[0] += p->pos[1]*p->vel[2] - p->pos[2]*p->vel[1]; */
/*  	btabp->l[1] += p->pos[2]*p->vel[0] - p->pos[0]*p->vel[2]; */
/*  	btabp->l[2] += p->pos[0]*p->vel[1] - p->pos[1]*p->vel[0]; */
/*  	Msgf(("Point mass gobbled m = %e; total = %e\nAccreted ang momentum = (%e, %e, %e)\n",  */
/*  	      p->mass, btabp->accmass,  */
/*  	      btabp->l[0], btabp->l[1], btabp->l[2])); */
/*        } */
    }
    Free(btab);
    StkCrunch(&s);
    *nobj = StkSz(&s)/sizeof(SPHbody);
    btab = StkBase(&s);
    *SPHbtabp = Realloc(btab, *nobj * sizeof(SPHbody));
}

void
AdjustBtab (SPHbody **SPHbtabp, int *nobj, int gnobj, windbody *windbtab, 
	    int windnobj, int windpartpershell, float r_limit, float dt, 
	    int iter, float tpos, int *added_particles)
{
    SPHbody *btab = *SPHbtabp;
    SPHbody *p;
    Stk s;
    SPHbody *q;
    unsigned int id;
    float r2 = r_limit*r_limit;
    float r_wind = 5.0; /* Decoupled from inner boundary size */
    float wr;
    double wpos[NDIM];

    float d = r_wind * 
	sqrt( 4.0-1.0 / (pow( sin( M_PI*(windpartpershell)/
				   (6.0*((windpartpershell)-2)) ), 2.0 )) );

    StkInitEz(&s);

    for (p = btab; p < btab+*nobj; p++) {
	/* Keep all particles outside of BH at origin, and
	   keep all particles inside reasonable volume of solution */

	if ( (Dot(p->pos, p->pos) >= r2)
	     && (fabs(p->pos[0]) <= 3400.0) 
	     && (fabs(p->pos[1]) <= 3400.0) 
	     && (fabs(p->pos[2]) <= 3400.0) 
	     && (p->u <= 1e5) ) { 

	    q = StkPush(&s, sizeof(SPHbody));
	    *q = *p;

	    if ( p->windid < windnobj ) {  /* Particle on inner shell? */
		VVV(wpos, = p->pos, - windbtab[p->windid].pos);
		wr = sqrt(Dot(wpos, wpos));

		if (wr > r_wind + 0.8*d){ /* Particle far from source? */
		    *added_particles = 1;  /* Indicate particle addition */

		    id = q->windid;
		    q->windid += windnobj;  /* Turn off addition for
					       recently pushed particle */
		    q = StkPush(&s, sizeof(SPHbody));

		    /* Be aware that some quantities not set here are set
		       only when exact_rho = 1 */

		    q->mass = p->mass;

		    VVS(q->pos, = wpos, * r_wind / wr);
		    VV(q->vel, = windbtab[id].vwind/r_wind*q->pos);
		    VV(q->pos, += windbtab[id].pos);

		    VVV(q->pos_last, = q->pos, - dt*q->vel);

		    q->h = 1.8*d;  /* Match calculation in writewind.c */

		    q->u = (p->u + windbtab[id].uwind) / 2.0;
		    q->udot = 0.0;
		    q->udot_last = 0.0;  /* Just in case */
		    q->pr = 0.0;  /* Fixed in update_intermediate */

		    VS(q->acc, = 0.0);
		    VS(q->acc_last, = 0.0);
		    VS(q->grav_acc, = 0.0);

		    q->nterms = 1;  /* Equivalent to SPHFixNterms */

		    q->tacc = -1e30;

		    /* Lots of possibly-unnecessary initializations */
		    /* Without diffusion, these should all stay 0 */
		    q->dt = q->dt_next = dt;  /* CORRECT?? */
		    q->min_nbr_dt = 1e30;  /* Just testing */
		    q->du = 0.0;
		    q->du_r = 0.0;
		    q->u_r = 0.0;
		    q->phi = 0.0;  /* Set this correctly? */

		    q->windid = id;
		    q->ident = 100000000;  /* Fix in call to SPHFixId? No */

/* 		    Msgf(("p->pos: %f %f %f; windid: %d; u: %e\n",  */
/* 			  q->pos[0], q->pos[1], q->pos[2], q->windid,  */
/* 			  q->u)); */
		}
	    }
	}
	/* Else track accreted/ejected material; do this right sometime */
	else {
	    Msgf(("%d: %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %u %u\n", iter, tpos, p->pos[0], p->pos[1], p->pos[2], p->vel[0], p->vel[1], p->vel[2], p->mass, p->rho, p->u, p->h, p->windid, p->ident));
	}

    }

    Free(btab);
    StkCrunch(&s);
    *nobj = StkSz(&s)/sizeof(SPHbody);
    btab = StkBase(&s);
    *SPHbtabp = Realloc(btab, *nobj * sizeof(SPHbody));
}


void
AdjustBtab2 (SPHbody **SPHbtabp, int *nobj, int gnobj, windbody *windbtab, 
	    int windnobj, float r_limit, float dt, int iter, float tpos,
	    int *added_particles, float *newmass)
{
    SPHbody *btab = *SPHbtabp;
    SPHbody *p;
    Stk s;
    SPHbody *q;
    float r2 = r_limit*r_limit;

    StkInitEz(&s);

    for (*newmass = 0.0, p = btab; p < btab+*nobj; p++) {
	/* Keep all particles outside of BH at origin, and
	   keep all particles inside reasonable volume of solution */

	if (Dot(p->pos, p->pos) >= r2) {
	    q = StkPush(&s, sizeof(SPHbody));
	    *q = *p;
	}
	/* Else add accreted material */
	else {
	    *newmass += p->mass;
	}
    }

    Free(btab);
    StkCrunch(&s);
    *nobj = StkSz(&s)/sizeof(SPHbody);
    btab = StkBase(&s);
    *SPHbtabp = Realloc(btab, *nobj * sizeof(SPHbody));
}


void
AdjustBtab3(SPHbody **SPHbtabp, int *nobj, int gnobj, float r_limit, 
	    float r_outer)
{
    SPHbody *btab = *SPHbtabp;
    SPHbody *p;
    Stk s;
    SPHbody *q;
    float r2 = r_limit*r_limit;
    float or2 = r_outer*r_outer;

    StkInitEz(&s);

    for (p = btab; p < btab+*nobj; p++) {
	/* Keep all particles outside of BH at origin and inside r_outer */

	if ( (Dot(p->pos, p->pos) >= r2) && (Dot(p->pos, p->pos) <= or2) ) {
	    q = StkPush(&s, sizeof(SPHbody));
	    *q = *p;
	}
    }

    Free(btab);
    StkCrunch(&s);
    *nobj = StkSz(&s)/sizeof(SPHbody);
    btab = StkBase(&s);
    *SPHbtabp = Realloc(btab, *nobj * sizeof(SPHbody));
}

/* adjust btab from snevolbrna to let the central particle acrete mass ~CIE */
void
AdjustBtab4(SPHbody **SPHbtabp, int *nobj, bndry_t b, float *newmass,
            float *newr, float *newp, float *newl, float G, float dt)
{
    SPHbody *btab = *SPHbtabp;
    SPHbody *p, *q;
    Stk s;
    float r1, r2, v2, b2, minb2 = 1e30;
    float j[NDIM], jhat[NDIM], r_vec[NDIM], v_vec[NDIM];
    float jm, jmax;
    float small = 1.e-12;
    float r_ns2, r_sw2, vel_i;
    float v_max, tff, m_accret;

    v_max = -0.033*C_LIGHT*tdivlCF;

    StkInitEz(&s);

    VS(newp, = 0.0);
    VS(newl, = 0.0);

    for (*newmass = 0.0, p = btab; p < btab+*nobj; p++) {

        /* r and v w.r.t. to bndry */
        VVV(v_vec, = p->vel, - b.vel);
        v2 = Dot(v_vec, v_vec);
        
        VVV(r_vec, = p->pos, - b.pos);
        r2 = Dot(r_vec, r_vec);
        r1 = sqrt(r2);

        /* One option: adjust r2 based on particle velocities to
           simulate capture-radius behavior */
        /* r2 = 4.0*newt*newt*b.mass*b.mass / (v2 * v2); */

        /* Another option: start small and move r2 out after eating
           all particles to 10% of the radius of the next-nearest
           particle */

	b2 = b.r;

        /* the purpose of this is to impose a sort of minimum bndry radius */
        /* miminum bndry_r: Neutron star radius of 10km */
        r_ns2 = 1.e6*ivlenCF; 

        /* "Schwarzschild radius"; i.e. where v_esc==0.1c */
        r_sw2 = 2.*GRAV_C*(b.mass*massCF)/(C_LIGHT*C_LIGHT*0.05); 
        r_sw2 = r_sw2*ivlenCF; /* convert to code-units */

        b2 = (r_ns2 > b2 ? r_ns2 : b2); /* pick the bigger one */ 
        b2 = (r_sw2 > b2 ? r_sw2 : b2); /* pick the bigger one */ 


        vel_i = Dot(v_vec, r_vec); /* vel of particle in bndry frame of ref. */
        vel_i = vel_i/r1;

        /* if within bndry.r || falling in at greater than 14000km/s and 
           within 5*bndry.r then eat particle */
        /*if ( (b2+p->h) >= r1 || (vel_i < v_max && r1 <= 2.*b2)) {*/
        if ( (b2+p->h) >= r1 ) {

            /* eat partial particle */
            /* outside bndry_r, but overlapping and not moving too fast */
            if( (fabs(b2-r1) < p->h) && (vel_i > v_max) ) {
                tff= sqrt( 2.*b.r*b.r*b.r/ (G*b.mass));
                m_accret = p->mass * dt/tff;

                /* partially eat */
                if( p->mass > m_accret) {
                    q = StkPush(&s, sizeof(SPHbody));
                    *q = *p;
                    q->mass = p->mass - m_accret;
                    p->mass = m_accret;
                } 

            }
            /* eat whole particle */

            *newmass += p->mass;

            VVS(newp, += v_vec, * p->mass);

            /* this assumes that the central particle is at the origin? ~CIE */
/*
            j[0] = p->mass * (p->pos[1]*p->vel[2] - p->pos[2]*p->vel[1]);
            j[1] = p->mass * (p->pos[2]*p->vel[0] - p->pos[0]*p->vel[2]);
            j[2] = p->mass * (p->pos[0]*p->vel[1] - p->pos[1]*p->vel[0]);
*/

            j[0] = p->mass * ( r_vec[1]*v_vec[2] - r_vec[2]*v_vec[1] );
            j[1] = p->mass * ( r_vec[2]*v_vec[0] - r_vec[0]*v_vec[2] );
            j[2] = p->mass * ( r_vec[0]*v_vec[1] - r_vec[1]*v_vec[0] );

            jm = sqrt( j[0]*j[0] + j[1]*j[1] + j[2]*j[2] );
            jhat[0] = j[0]/(jm + small);
            jhat[1] = j[1]/(jm + small);
            jhat[2] = j[2]/(jm + small);

            jmax = sqrt( G * b.mass * b.r ) / p->mass;

            jm = ( jm < jmax ? jm : jmax );

            VVS(j, = jhat, * jm);

            VV(newl, += j);

            Msgf(("dt: %g: #%d: m: %g; x: %g; y: %g; z: %g; vx: %g; vy: %g; vz: %g\n", dt, p->ident, p->mass, p->pos[0], p->pos[1], p->pos[2], p->vel[0], p->vel[1], p->vel[2]));
        } else { /* dont eat particle*/
            q = StkPush(&s, sizeof(SPHbody));
            *q = *p;
            if (r1 < minb2) minb2 = r1;
        }
    }

    Free(btab);
    StkCrunch(&s);
    *nobj = StkSz(&s)/sizeof(SPHbody);
    btab = StkBase(&s);
    *SPHbtabp = Realloc(btab, *nobj * sizeof(SPHbody));

    *newr = 0.8*minb2;  /* Candidate new boundary radius =
                                 innermost particle's
                                 distance-to-boundary * 25% */
    if (*newr < b.r) *newr = b.r;  /* Never shrink boundary */

}


void 
AddWinds(SPHbody **SPHbtabp, int *nobj, template_t *temptab, 
	 int windpartpershell, float r_wind, float v_wind, float mdot_wind, 
	 float u_wind, float *t_wind, float tpos, float dt, float *dt_next,
	 float openangle_wind)
{
    SPHbody *btab = *SPHbtabp;
    SPHbody *p;
    Stk s;
    SPHbody *q;
    template_t *t;
    float theta, phi;
    double x, y, z;
    float d = r_wind * 
	sqrt( 4.0-1.0 / (pow( sin( M_PI*(windpartpershell)/
				   (6.0*((windpartpershell)-2)) ), 2.0 )) );

    StkInitEz(&s);

    /* Push all existing particles onto stack */
    for (p = btab; p < btab + *nobj; ++p) {
	q = StkPush(&s, sizeof(SPHbody));
	*q = *p;
    }

    /* If enough time has passed since last addition, add another
       shell of wind particles */
    if (tpos - *t_wind >= 1.8*d/v_wind) {  

	/* Pick random theta, phi */
	theta = acos(1.0 - 2.0*drand48());  /* so that p(y) = sin(y) */
	phi = 2.0*M_PI*drand48();

	/* Loop over shell template, add particles */
	for (t = temptab; t < temptab + windpartpershell; ++t) {
	    q = StkPush(&s, sizeof(SPHbody));
	
	    q->mass = mdot_wind/windpartpershell * (tpos - *t_wind);

	    /* Rotate by (theta, phi) */	    
	    /*comment out if no density variation, also comment out 
	      openangle_wind initialization upstairs. -CE*/
 	    q->pos[0] = t->pos[0]*sin(theta)*cos(phi); 
 	    q->pos[1] = t->pos[1]*sin(theta)*sin(phi); 
 	    q->pos[2] = t->pos[2]*cos(theta); 

	    VV(q->pos, = t->pos);

	    x = cos(theta)*q->pos[0] - sin(theta)*q->pos[2];
	    z = sin(theta)*q->pos[0] + cos(theta)*q->pos[2];

	    q->pos[0] = x;
	    q->pos[2] = z;

	    x = cos(phi)*q->pos[0] - sin(phi)*q->pos[1];
	    y = sin(phi)*q->pos[0] + cos(phi)*q->pos[1];

	    q->pos[0] = x;
	    q->pos[1] = y;

	    if ( cos(openangle_wind/180.0*M_PI) > 
		 fabs(q->pos[2])/sqrtf_fast(Dot(q->pos, q->pos)) ) {
		StkPop(&s, sizeof(SPHbody));
		continue;
	    }

	    VS(q->pos, *= r_wind);

	    VV(q->vel, = v_wind/r_wind*q->pos);  /* Outward radial vel */
	    
	    VVV(q->pos_last, = q->pos, - dt*q->vel);
	    
	    q->h = 1.8*d;  /* Match calculation in writewind.c */
	    
	    q->u = u_wind;
	    q->udot = 0.0;
	    q->udot_last = 0.0;  /* Just in case */
	    q->pr = 0.0;  /* Fixed in update_intermediate */
	    
	    VS(q->acc, = 0.0);
	    VS(q->acc_last, = 0.0);
	    VS(q->grav_acc, = 0.0);
	    
	    q->nterms = 1;  /* Equivalent to SPHFixNterms */
	    
	    q->tacc = -1e30;
	    
	    /* Lots of possibly-unnecessary initializations */
	    /* Without diffusion, these should all stay 0 */
	    q->dt = q->dt_next = dt;  /* CORRECT?? */
	    q->min_nbr_dt = 1e30;  /* Just testing */
	    q->du = 0.0;
	    q->du_r = 0.0;
	    q->u_r = 0.0;
	    q->phi = 0.0;  /* Set this correctly? */
	    
	    q->ident = (*nobj)++;
	}
	
	*t_wind = tpos;
	*dt_next = 1.8*d/v_wind;
    }

    Free(btab);
    StkCrunch(&s);
    *nobj = StkSz(&s)/sizeof(SPHbody);
    btab = StkBase(&s);
    *SPHbtabp = Realloc(btab, *nobj * sizeof(SPHbody));
}


void ReadTemplate(char *filename, template_t **temptab, int *tempnobj) 
{
    FILE *fp;
    int i;
    double r;

    if ( (fp = fopen(filename, "r")) == NULL )
	Error("Can't open %s: %s\n", filename, strerror(errno));

    fscanf(fp, "%d %*g", tempnobj);  /* Read number of lines in file */

    *temptab = (template_t *)Malloc(*tempnobj * sizeof(template_t));

    for(i = 0; i < *tempnobj; ++i) {
	if (fscanf(fp, "%lg %lg %lg", &((*temptab)[i].pos[0]),
		   &((*temptab)[i].pos[1]), &((*temptab)[i].pos[2])) != 3)
	    Error("Error reading positions from %s\n", filename);
	r = sqrt(Dot((*temptab)[i].pos, (*temptab)[i].pos));
	(*temptab)[i].pos[0] /= r;
	(*temptab)[i].pos[1] /= r;
	(*temptab)[i].pos[2] /= r;
    }

    fclose(fp);
}


const double MAS = 1.9889e27;    /* g */
const double LEN = 1.0e14;       /* cm */
const double TIM = 3600*24*365;  /* s */


void ReadWindData(char *filename, winddata_t **wdata, int *wnobj)
{
    FILE *fp;
    char input[100];
    int i;
    float t, dt, dt2, mdot, v_inf;

    if ( (fp = fopen(filename, "r")) == NULL )
	Error("Can't open %s: %s\n", filename, strerror(errno));

    fgets(input, sizeof(input), fp);  /* Read and discard header line */
    if (ferror(fp))
	Error("%s: %s\n", filename, strerror(errno));

    *wdata = NULL;
    *wnobj = 0;
    /* BEWARE OF UNITS BELOW!!! */
    /* assume it's in cgs, convert to user-units */
    while (fscanf(fp, "%g %g %*g %*g %g %g", &t, &dt, &mdot, &v_inf) == 4) {
	(*wnobj)++;
	*wdata = (winddata_t *)Realloc(*wdata, *wnobj * sizeof(winddata_t));
	(*wdata)[*wnobj - 1].t = dt/timeCF;
	(*wdata)[*wnobj - 1].mdot = -mdot*timeCF/massCF;  /* Flip sign from input */
	(*wdata)[*wnobj - 1].v_inf = v_inf*tdivlCF;
	(*wdata)[*wnobj - 1].u = 0.102547;  /* From src/winds/proto */
    }

    if (ferror(fp))
	Error("%s: %s\n", filename, strerror(errno));

    Msgf(("Read %d winddata_t's from %s\n", *wnobj, filename));

    dt = (*wdata)[0].t;
    (*wdata)[0].t = 0.0;

    for(i = 1; i < *wnobj; ++i) {
	dt2 = (*wdata)[i].t;
	(*wdata)[i].t = (*wdata)[i-1].t + dt;
	dt = dt2;
    }

    fclose(fp);    
}


void AddNonconstWinds(SPHbody **SPHbtabp, int *nobj, template_t *temptab, 
		      int windpartpershell, winddata_t *wdata, int wnobj, 
		      float r_wind, float r_outer, float *t_wind, float tpos, 
		      float dt, float *dt_next, float openangle_wind)
{
    SPHbody *btab = *SPHbtabp;
    SPHbody *p;
    Stk s;
    SPHbody *q;
    template_t *t;
    float theta, phi;
    double x, y, z;
    float v_wind, mdot_wind, u_wind;
    float d = r_wind * 
	sqrt( 4.0-1.0 / (pow( sin( M_PI*(windpartpershell)/
				   (6.0*((windpartpershell)-2)) ), 2.0 )) );

    StkInitEz(&s);

    /* Push all existing particles inside r = r_outer onto stack */
    for (p = btab; p < btab + *nobj; ++p) {
	if (Dot(p->pos, p->pos) <= r_outer*r_outer) {
	    q = StkPush(&s, sizeof(SPHbody));
	    *q = *p;
	}
    }

    if (MPMY_Procnum() == 0) {
	hunt(wdata, wnobj, *t_wind, &v_wind, &mdot_wind, &u_wind);

	/* If enough time has passed since last addition, add another
	   shell of wind particles */
	if (tpos - *t_wind >= 1.8*d/v_wind) {  

	    hunt(wdata, wnobj, tpos, &v_wind, &mdot_wind, &u_wind);

	    Msgf(("t = %g; v = %g; mdot = %g; u = %g\n", tpos, v_wind,
		  mdot_wind, u_wind));

	    /* Pick random theta, phi */
	    theta = acos(1.0 - 2.0*drand48());  /* so that p(y) = sin(y) */
	    phi = 2.0*M_PI*drand48();

	    /* Loop over shell template, add particles */
	    for (t = temptab; t < temptab + windpartpershell; ++t) {
		q = StkPush(&s, sizeof(SPHbody));
	
		q->mass = mdot_wind/windpartpershell * (tpos - *t_wind);

		/* Rotate by (theta, phi) */	    
		VV(q->pos, = t->pos);

		x = cos(theta)*q->pos[0] - sin(theta)*q->pos[2];
		z = sin(theta)*q->pos[0] + cos(theta)*q->pos[2];
		q->pos[0] = x;
		q->pos[2] = z;

		x = cos(phi)*q->pos[0] - sin(phi)*q->pos[1];
		y = sin(phi)*q->pos[0] + cos(phi)*q->pos[1];
		q->pos[0] = x;
		q->pos[1] = y;

		if ( cos(openangle_wind/180.0*M_PI) > 
		     fabs(q->pos[2])/sqrtf_fast(Dot(q->pos, q->pos)) ) {
		    StkPop(&s, sizeof(SPHbody));
		    continue;
		}

		VS(q->pos, *= r_wind);

		VV(q->vel, = v_wind/r_wind*q->pos);  /* Outward radial vel */
		VVV(q->pos_last, = q->pos, - dt*q->vel);
	    
		q->h = 1.8*d;  /* Match calculation in writewind.c */
	    
		q->u = u_wind;
		q->udot = 0.0;
		q->udot_last = 0.0;  /* Just in case */
		q->pr = 0.0;  /* Fixed in update_intermediate */
	    
		VS(q->acc, = 0.0);
		VS(q->acc_last, = 0.0);
		VS(q->grav_acc, = 0.0);
	    
		q->nterms = 1;  /* Equivalent to SPHFixNterms */
	    
		q->tacc = -1e30;
	    
		/* Lots of possibly-unnecessary initializations */
		/* Without diffusion, these should all stay 0 */
		q->dt = q->dt_next = dt;  /* CORRECT?? */
		q->min_nbr_dt = 1e30;  /* Just testing */
		q->du = 0.0;
		q->du_r = 0.0;
		q->u_r = 0.0;
		q->phi = 0.0;  /* Set this correctly? */
	    
		q->ident = (*nobj)++;
	    }
	
	    *t_wind = tpos;
	    *dt_next = 1.8*d/v_wind;
	}
    }

    Free(btab);
    StkCrunch(&s);
    *nobj = StkSz(&s)/sizeof(SPHbody);
    btab = StkBase(&s);
    *SPHbtabp = Realloc(btab, *nobj * sizeof(SPHbody));
}


void hunt(winddata_t *w, int wnobj, float t, float *v, float *mdot, float *u)
{
    int i = wnobj/2;
    int il = 0, ih = wnobj-2;
    float s;

    while( !((w[i].t <= t) && (w[i+1].t >= t)) ) {
	if (w[i].t < t)
	    il = i;
	else if (w[i].t > t)
	    ih = i;
	i = (il + ih)/2;
    }

    /* Need to interpolate here */
    s = (w[i+1].v_inf - w[i].v_inf) / (w[i+1].t - w[i].t);
    *v = w[i].v_inf + s * (t - w[i].t);

    s = (w[i+1].mdot - w[i].mdot) / (w[i+1].t - w[i].t);
    *mdot = w[i].mdot + s * (t - w[i].t);

    s = (w[i+1].u - w[i].u) / (w[i+1].t - w[i].t);
    *u = w[i].u + s * (t - w[i].t);
}


void 
AddAccreting(SPHbody **SPHbtabp, int *nobj, template_t *temptab, 
	     int windpartpershell, float r_wind, float v_r, 
	     float mdot_wind, float u_wind, float *t_wind, float tpos, 
	     float dt, float *dt_next, float omega)
{
    SPHbody *btab = *SPHbtabp;
    SPHbody *p;
    Stk s;
    SPHbody *q;
    template_t *t;
    float theta, phi;
    double x, y, z;
    float d = r_wind * 
	sqrt( 4.0-1.0 / (pow( sin( M_PI*(windpartpershell)/
				   (6.0*((windpartpershell)-2)) ), 2.0 )) );

    StkInitEz(&s);

    /* Push all existing particles onto stack */
    for (p = btab; p < btab + *nobj; ++p) {
	q = StkPush(&s, sizeof(SPHbody));
	*q = *p;
    }

    /* If enough time has passed since last addition, add another
       shell of wind particles */
    if (tpos - *t_wind >= 0.9*d/fabs(v_r)) {  

	/* Pick random theta, phi */
	theta = acos(1.0 - 2.0*drand48());  /* so that p(y) = sin(y) */
	phi = 2.0*M_PI*drand48();

	/* Loop over shell template, add particles */
	for (t = temptab; t < temptab + windpartpershell; ++t) {
	    q = StkPush(&s, sizeof(SPHbody));
	
	    q->mass = mdot_wind/windpartpershell * (tpos - *t_wind);

	    /* Rotate by (theta, phi) */	    
/* 	    q->pos[0] = t->pos[0]*sin(theta)*cos(phi); */
/* 	    q->pos[1] = t->pos[1]*sin(theta)*sin(phi); */
/* 	    q->pos[2] = t->pos[2]*cos(theta); */

	    VV(q->pos, = t->pos);

	    x = cos(theta)*q->pos[0] - sin(theta)*q->pos[2];
	    z = sin(theta)*q->pos[0] + cos(theta)*q->pos[2];

	    q->pos[0] = x;
	    q->pos[2] = z;

	    x = cos(phi)*q->pos[0] - sin(phi)*q->pos[1];
	    y = sin(phi)*q->pos[0] + cos(phi)*q->pos[1];

	    q->pos[0] = x;
	    q->pos[1] = y;

	    VS(q->pos, *= r_wind);

	    VV(q->vel, = v_r/r_wind*q->pos);  /* Radial velocity */

	    /* Add v_phi = sqrt(x*x+y*y)*omega */
	    q->vel[0] -= q->pos[1]*omega;
	    q->vel[1] += q->pos[0]*omega;
	    
	    VVV(q->pos_last, = q->pos, - dt*q->vel);
	    
	    q->h = 1.7*d;  /* Match calculation in writewind.c */
	    
	    q->u = u_wind;
	    q->udot = 0.0;
	    q->udot_last = 0.0;  /* Just in case */
	    q->pr = 0.0;  /* Fixed in update_intermediate */
	    
	    VS(q->acc, = 0.0);
	    VS(q->acc_last, = 0.0);
	    VS(q->grav_acc, = 0.0);
	    
	    q->nterms = 1;  /* Equivalent to SPHFixNterms */
	    
	    q->tacc = -1e30;
	    
	    /* Lots of possibly-unnecessary initializations */
	    /* Without diffusion, these should all stay 0 */
	    q->dt = q->dt_next = dt;  /* CORRECT?? */
	    q->min_nbr_dt = 1e30;  /* Just testing */
	    q->du = 0.0;
	    q->du_r = 0.0;
	    q->u_r = 0.0;
	    q->phi = 0.0;  /* Set this correctly? */
	    
	    q->ident = (*nobj)++;
	}
	
	*t_wind = tpos;
	*dt_next = 0.9*d/v_r;
    }

    Free(btab);
    StkCrunch(&s);
    *nobj = StkSz(&s)/sizeof(SPHbody);
    btab = StkBase(&s);
    *SPHbtabp = Realloc(btab, *nobj * sizeof(SPHbody));
}
