#include <stddef.h>
#include <math.h>
#include "SDF.h"
#include "mpmy.h"
#define NDIM 3
#include "vop.h"
#include "fastflpt.h"
#include "singlio.h"
#include "SDFread.h"
#include "SDFwrite.h"
#include "error.h"
#include "ring.h"
#include "randoms.h"
#include "bigmalloc.h"

void SDFwrite(const char *filename, int gnobj, int nobj, const void *btab,
	    int bsize, const char *bodydesc, 
	      /* const char *name, SDF_type_enum type, <type> val */ ...);
void memfile_init(int sz);
void memfile_vfprintf(void *junk, const char *fmt, va_list args);

typedef struct {
    float mass;			/* mass of body */
    float pos[NDIM];		/* position of body */
    float vel[NDIM];		/* velocity of body */
    float rho;
} body;

#define HIST_BINS 25		/* remember to change OUTBODYDESC */
#define RHO_BINS 10

typedef struct {
    float mass;			/* mass of body */
    float pos[NDIM];		/* position of body */
    float vel[NDIM];		/* velocity of body */
    float mass_hist[RHO_BINS][HIST_BINS];
    int   nobj_hist[RHO_BINS][HIST_BINS];
    float rvel_hist[RHO_BINS][HIST_BINS];
    float rvel2_hist[RHO_BINS][HIST_BINS];
    float vel2_hist[RHO_BINS][HIST_BINS];
} cbody;

#define OUTBODYDESC \
"struct {\n\
    float mass;			/* mass of body */\n\
    float x, y, z;			/* position of body */\n\
    float mass_hist[25];	/* histogram */\n\
    int   nobj_hist[25];	/* count histogram */\n\
    float rvel_hist[25];	/* radial velocity histogram */\n\
    float rvel2_hist[25];	/* radial velocity^2 histogram */\n\
    float vel2_hist[25];	/* radial velocity^2 histogram */\n\
}"

static float Max_radius2, Min_radius2;
static float Physical_L0, LogFactor;
static float Hubble;

static void
setup_corr(float min, float max, float l0, float hubble)
{
    Physical_L0 = l0;
    Max_radius2 = max*max;
    Min_radius2 = min*min;
    LogFactor = HIST_BINS/log(Max_radius2/Min_radius2);
    Hubble = hubble;
}


void
init(void *o, void *p)
{
    memcpy(o, p, sizeof(body));
}

void 
do_corr(void *p0, void *list, int bsize, int n)
{
    unsigned int s;
    int i, j, k;
    cbody *p = p0;
    body *q;
    void *last = (char *)list + bsize*n;
    Vxd(float dr);
    float oneor, r2, partial_r2;
    Vxd(float dv);
    float vv;
    /* Try to get these into regisiters */
    VxdV(const float ppos, = p->pos);
    VxdV(const float pvel, = p->vel);
    const float lfac = LogFactor;
    const float imr2 = 1.0/Min_radius2;
    const float minr2 = Min_radius2;
    const float maxr2 = Max_radius2*(1.0-1e-6);	/* avoid roundoff errors */
    const float offset[NDIM] = {-Physical_L0, (float)0.0, Physical_L0};

    while (list < last) {
	q = list;
	for (i = 0; i < NDIM; i++) {
	    dr0 = ppos0 - q->pos[0] + offset[i];
	    partial_r2 = dr0*dr0;
	    if (partial_r2 >= maxr2) continue;
	    dv0 = pvel0 - q->vel[0] + Hubble*offset[i] - Hubble*dr0;
	    for (j = 0; j < NDIM; j++) {
		dr1 = ppos1 - q->pos[1] + offset[j];
		partial_r2 = dr0*dr0 + dr1*dr1;
		if (partial_r2 >= maxr2) continue;
		dv1 = pvel1 - q->vel[1] + Hubble*offset[j] - Hubble*dr1;
		for (k = 0; k < NDIM; k++) {
		    dr2 = ppos2 - q->pos[2] + offset[k];
		    dv2 = pvel2 - q->vel[2] + Hubble*offset[k] - Hubble*dr2;
		    r2 = partial_r2 + dr2*dr2;

		    if (r2 < maxr2) {
			if (r2 < minr2)
			  s = 0;
			else
#ifdef __PARAGON__
			  s = logf(r2*imr2) * lfac;
#else
			  s = log(r2*imr2) * lfac;
#endif
			if (s >= HIST_BINS) {
			    Shout("s is %d\n", s);
			    continue;
			}
			p->mass_hist[s] += q->mass;
			p->nobj_hist[s]++;

			vv = Dotx(dv, dr);
			oneor = recipsqrtf(r2);
			p->rvel_hist[s] += vv*oneor;
			vv *= vv;
			p->rvel2_hist[s] += vv*oneor*oneor;
			p->vel2_hist[s] += Dotx(dv, dv);
		    }
		}
	    }
	}
	list = (char *)list + bsize;
    }
}

void
main(int argc, char *argv[])
{
    int massconf, rhoconf;
    int xconf, yconf, zconf;
    int vxconf, vyconf, vzconf;
    char name[256], outname[256];
    SDF *csdfp, *sdfp;
    body *btab;
    cbody *ctab;
    int i, j;
    int nobj, gnobj;
    int cnobj, cgnobj;
    float tpos, redshift, R0;
    int iter;
    int seed;
    int samples, nsamples;
    float win_min, win_max;
    float physical_l0, physical_vol;
    float radius[HIST_BINS], mean_radius[HIST_BINS];
    float mass_corr[RHO_BINS][HIST_BINS];
    int nobj_corr[RHO_BINS][HIST_BINS];
    float nobj_corrf[RHO_BINS][HIST_BINS];
    float rvel_corr[RHO_BINS][HIST_BINS];
    float rvel2_corr[RHO_BINS][HIST_BINS];
    float vel2_corr[RHO_BINS][HIST_BINS];
    float massb[HIST_BINS], nobjb[HIST_BINS];
    float total_mass = 0.0;
    float shrink_fac;
    float p;
    float hubble;
    ran_state ranstate;
    char mass_name[64];
 
    MPMY_Init(&argc, &argv);
    if (argc != 2) {
	singlPrintf("usage: %s ctlfile\n", argv[0]);
	exit(1);
    }
    singlPrintf("Welcome to the machine\n");
    
    csdfp = SDFopen(0, argv[1]);
    if (csdfp == 0)
      SinglError("Sorry, couldn't SDFopen %s\n%s\n", argv[2], SDFerrstring);

    SDFgetstring(csdfp, "datafile", name, sizeof(name));
    SDFgetstring(csdfp, "outfile", outname, sizeof(outname));
    SDFgetstring(csdfp, "mass_name", mass_name, sizeof(mass_name));

    singlPrintf("Reading \"%s\"\n", name);
    sdfp = SDFread(csdfp, (void **)&btab, &gnobj, &nobj, sizeof(body),
		   mass_name, offsetof(body, mass), &massconf,
		   "x", offsetof(body, pos[0]), &xconf,
		   "y", offsetof(body, pos[1]), &yconf,
		   "z", offsetof(body, pos[2]), &zconf,
		   "vx", offsetof(body, vel[0]), &vxconf,
		   "vy", offsetof(body, vel[1]), &vyconf,
		   "vz", offsetof(body, vel[2]), &vzconf,
		   rho, offsetof(body, rho), &rhoconf,
		   NULL);
    if( massconf==0 || rhoconf == 0 || xconf==0 || yconf==0 || zconf==0 ){
	SinglError("Could not find %s %s %s %s in data file!\n",
		   (massconf==0)? "mass" : "",
		   (rhoconf==0)? "rho" : "",
		   (xconf==0)? "x" : "",
		   (yconf==0)? "y" : "",
		   (zconf==0)? "z" : "");
    }
    if( vxconf==0 || vyconf==0 || vzconf==0 ){
	SinglError("Could not find %s %s %s in data file!\n",
		   (vxconf==0)? "vx" : "",
		   (vyconf==0)? "vy" : "",
		   (vzconf==0)? "vz" : "");
    }
    SDFgetfloatOrDefault(sdfp, "tpos",  &tpos, (float)0.0);
    SDFgetfloatOrDefault(sdfp, "redshift",  &redshift, 0.0);
    SDFgetfloatOrDefault(sdfp, "hubble",  &hubble, 0.0511);
    SDFgetfloatOrDefault(sdfp, "R0",  &R0, 0.0);
    if (R0 == 0.0) 
        SDFgetfloatOrDefault(csdfp, "R0",  &R0, 0.0);
    SDFgetintOrDefault  (sdfp, "iter",  &iter, 0);

    SDFgetfloatOrDefault(csdfp, "win_max", &win_max, R0);
    SDFgetfloatOrDefault(csdfp, "win_min", &win_min, R0/100.0);
    SDFgetintOrDefault(csdfp, "samples", &samples, -1);
    SDFgetintOrDefault(csdfp, "seed", &seed, 123);
    SDFgetfloatOrDefault(csdfp, "hubble",  &hubble, hubble);

    SDFclose(csdfp);
    SDFclose(sdfp);

    for (i = 0; i < nobj; i++) {
	total_mass += btab[i].mass;
    }

    if (samples == -1) {
	ctab = Malloc(nobj*sizeof(cbody));
	p = shrink_fac = 1.0;
    } else {
	shrink_fac = (float) samples/gnobj;
	ran_init(seed*(MPMY_Procnum()+1), &ranstate);
	/* Add a bit more than we will need, hopefully */
	ctab = Malloc((shrink_fac*nobj+1000)*sizeof(cbody));
	p = 2.0 * 1.0/shrink_fac;
    }
    cnobj = 0;

    for (i = 0; i < nobj; i++, cnobj++) {
	if (cnobj >= shrink_fac*nobj+1000) Error("ctab overflowed\n");
	ctab[cnobj].mass = btab[i].mass;
	VV(ctab[cnobj].pos, = btab[i].pos);
	VV(ctab[cnobj].vel, = btab[i].vel);
	for (k = 0; k < RHO_BINS; k++) {
	  for (j = 0; j < HIST_BINS; j++) {
	    ctab[cnobj].mass_hist[k][j] = (float)0.0;
	    ctab[cnobj].nobj_hist[k][j] = 0;
	    ctab[cnobj].rvel_hist[k][j] = (float)0.0;
	    ctab[cnobj].rvel2_hist[k][j] = (float)0.0;
	    ctab[cnobj].vel2_hist[k][j] = (float)0.0;
	  }
	}
	if (samples != -1) i += (int)(p * uniform_rand(&ranstate));
    }
    ctab = Realloc(ctab, cnobj*sizeof(cbody));
    MPMY_Combine(&cnobj, &cgnobj, 1, MPMY_INT, MPMY_SUM);
    nsamples = cgnobj;
    singlPrintf("Doing %d centers\n", nsamples);

    physical_l0 = 2.0*R0;
    win_max /= 1.0+redshift;
    win_min /= 1.0+redshift;
    physical_l0 /= 1.0+redshift;
    physical_vol = physical_l0*physical_l0*physical_l0;
    
    setup_corr(win_min, win_max, physical_l0, hubble);

    Ring(ctab, sizeof(cbody), cnobj, btab, sizeof(body), nobj, sizeof(body),
	 init, do_corr);

    for (k = 0; k < RHO_BINS; k++) {
      for (j = 0; j < HIST_BINS; j++) {
	mass_corr[k][j] = 0.0;
	nobj_corr[k][j] = 0;
	rvel_corr[k][j] = 0.0;
	rvel2_corr[k][j] = 0.0;
	vel2_corr[k][j] = 0.0;
    }

    for (i = 0; i < cnobj; i++) {
	for (k = 0; k < RHO_BINS; k++) {
	  for (j = 0; j < HIST_BINS; j++) {
	    mass_corr[k][j] += ctab[i].mass_hist[k][j];
	    nobj_corr[k][j] += ctab[i].nobj_hist[k][j];
	    rvel_corr[k][j] += ctab[i].rvel_hist[k][j];
	    rvel2_corr[k][j] += ctab[i].rvel2_hist[k][j];
	    vel2_corr[k][j] += ctab[i].vel2_hist[k][j];
	}
    }
    MPMY_Combine(&total_mass, &total_mass, 1, MPMY_FLOAT, MPMY_SUM);

    for (k = 0; k < RHO_BINS; k++) {
      for (i = 0; i < HIST_BINS; i++) {
	nobj_corrf[k][i] = nobj_corr[k][i];
	mean_radius[i] = exp(0.5*(i+0.5)/LogFactor) * win_min;
	radius[i] = exp(0.5*(i+1.0)/LogFactor) * win_min;
	massb[k][i] = (4.0*M_PI/3.0) * (total_mass[k]/physical_vol)
	  * radius[i] * radius[i] * radius[i];
	nobjb[k][i] = (4.0*M_PI/3.0) * (gnobj[k]/physical_vol)
	  * radius[i] * radius[i] * radius[i];
      }
    }

    MPMY_Combine(mass_corr, mass_corr, RHO_BINS*HIST_BINS, MPMY_FLOAT, MPMY_SUM);
    MPMY_Combine(nobj_corr, nobj_corr, RHO_BINS*HIST_BINS, MPMY_INT, MPMY_SUM);
    MPMY_Combine(nobj_corrf, nobj_corrf, RHO_BINS*HIST_BINS, MPMY_FLOAT, MPMY_SUM);
    MPMY_Combine(rvel_corr, rvel_corr, RHO_BINS*HIST_BINS, MPMY_FLOAT, MPMY_SUM);
    MPMY_Combine(rvel2_corr, rvel2_corr, RHO_BINS*HIST_BINS, MPMY_FLOAT, MPMY_SUM);
    MPMY_Combine(vel2_corr, vel2_corr, RHO_BINS*HIST_BINS, MPMY_FLOAT, MPMY_SUM);


    for (k = 0; k < RHO_BINS; k++) {
      rvel_corr[k][0] /= (nobj_corrf[k][0]+1e-6);
      rvel2_corr[k][0] /= (nobj_corrf[k][0]+1e-6);
      vel2_corr[k][0] /= (nobj_corrf[k][0]+1e-6);
      mass_corr[k][0] /= nsamples*massb[k][0];
      nobj_corrf[k][0] /= nsamples*nobjb[k][0];
      for (i = 1; i < HIST_BINS; i++) {
	rvel_corr[k][i] /= (nobj_corrf[k][i]+1e-6);
	rvel2_corr[k][i] /= (nobj_corrf[k][i]+1e-6);
	vel2_corr[k][i] /= (nobj_corrf[k][i]+1e-6);
	mass_corr[k][i] /= nsamples*(massb[k][i]-massb[k][i-1]);
	nobj_corrf[k][i] /= nsamples*(nobjb[k][i]-nobjb[k][i-1]);
      }
    }


    singlPrintf("# input datfile %s\n", name);
    singlPrintf("# Mass variable is = %s\n", mass_name);
    singlPrintf("# float win_min = %g\n", win_min*(1.0+redshift));
    singlPrintf("# float win_max = %g\n", win_max*(1.0+redshift));
    singlPrintf("# float physical_l0 = %g\n", physical_l0);
    singlPrintf("# float redshift = %g\n", redshift);
    singlPrintf("# float hubble = %g\n", hubble);
    singlPrintf("# float total_mass = %g\n", total_mass);
    singlPrintf("# int gnobj = %d\n", gnobj);
    singlPrintf("# int nsamples = %d\n", nsamples);
    singlPrintf("# values for radius and velocities are comoving\n");
    singlPrintf("# bin_num radius nobj_corr mass_corr v_rad v_rad_rms v_rms n\n");

      for (j = 0; j < HIST_BINS; j++) {
	singlPrintf("%2d %10g %10g %10g %10g %10g %10g %d\n", 
		    j, mean_radius[j]*(1.0+redshift), nobj_corrf[j], 
		    mass_corr[j], rvel_corr[j]/(1.0+redshift), 
		    sqrt(rvel2_corr[j])/(1.0+redshift), 
		    sqrt(vel2_corr[j])/(1.0+redshift), 
		    nobj_corr[j]);
    }
    exit(0);
}
