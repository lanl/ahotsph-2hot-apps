#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "protos.h"
#include "timers.h"
#include "singlio.h"
#include "mpmy.h"
#include "mpmy_io.h"
#include "randoms.h"
#include "Assert.h"
#include "bigmalloc.h"
#include "Msgs.h"
#include "SDF.h"
#include "error.h"
#include "getparam.h"
#include "SDFwrite.h"

#define MAXNDIM 3
#define PKTSIZE (4*16384)

#define GNEWT 4.49865897e4 /* units of kpc^3 / 10^10Msolar-Gyr^2 */
#define one_kpc (3.08567802e16) /* km */
#define one_Gyr (3.1558149984e16) /* sec */
#ifndef M_PI
#define	M_PI	3.14159265358979323846
#endif

#define Rindex(i,j,k,n1,n2) ((((i)*(n1)+(j))*(n2)+(k)))

ran_state Ranstate;
Timer_t Tot, TotWC;
float ran1(long *idum);
void Holtzman_cdm(int i, int j, int k, float *real, float *imag);
void EBW_cdm(int i, int j, int k, float *real, float *imag);
void powlaw(int i, int j, int k, float *real, float *imag);
void ranp_setup(long seed, long n, long max, ran_state *rs);
void setup_Holtzman_cdm(float L0, float Omega0, float h, float T0, float Tquad,
	       float t2, float t3, float t4, float t5, ran_state *ranstate);
void setup_EBW_cdm(float L0, float Omega0, float h, float T0, float Tquad,
	       float nu, float a, float b, float c, ran_state *ranstate);
void setup_powlaw(float L0, float Omega0, float h, float T0, float Tquad,
		  float n, ran_state *rs);
void prft(unsigned long *nn, float *data, ran_state *rs, 
     void spectrum(int, int, int, float *, float *));
void ft(unsigned long *nn, float *data, ran_state *rs, 
     void spectrum(int, int, int, float *, float *));

float t_from_Z(float Omega0, float H0, float Lambda_prime, float Z);
float growthfac_from_Z(float Omega0, float H0, float Lambda_prime, float Z);
float velfac_from_Z(float Omega0, float H0, float Lambda_prime, float Z);
float hubble_from_Z(float Omega0, float H0, float Lambda_prime, float Z);
float Wsq(int i, int j, int k);
static void shrink(float **datahndl, int *gnpts, int *npts, unsigned long nn[3], int ir);
static void radius_cut(unsigned long nn[3], float *data, float rmin, 
		       float rmax);
static void center_data(unsigned long nn[3], float *data, int center[3]);
void maxidx_f(const void *a, const void *b, void *c);
static void sinc_setup(unsigned long npts[]);
static float *sinctbl[MAXNDIM];
static int Ndim;

struct max_st {
    float max;
    int index;
};

#define MASSDESC \
"struct {\n\
    float drho;		/* mass = unit_mass * (1.0 + drho) */\n\
}"

typedef struct {
    float mass;			/* mass of body */
    float pos[3];		/* position of body */
} body, *bodyptr;

#define BODYDESC \
"struct {\n\
    float mass;			/* mass of body */\n\
    float x, y, z;		/* position of body */\n\
}"

#define HDRPARAMS \
     "ICdesc", SDF_STRING, ICdesc, \
     "ICparams", SDF_STRING, ICparams, \
     "Gnewt", SDF_FLOAT, GNEWT, \
     "Omega0", SDF_FLOAT, Omega0, \
     "Lambda_prime", SDF_FLOAT, Lambda_prime, \
     "Zinitial", SDF_FLOAT, Zinitial, \
     "redshift", SDF_FLOAT, Zinitial, \
     "tpos", SDF_FLOAT, tpos, \
     "H0", SDF_FLOAT, H0, \
     "h_100", SDF_FLOAT, h_100, \
     "hubble", SDF_FLOAT, H, \
     "growth_fac", SDF_FLOAT, growth_fac, \
     "velocity_fac", SDF_FLOAT, velocity_fac, \
     "Tquad", SDF_FLOAT, Tquad, \
     "mtot", SDF_FLOAT, massout, \
     "unit_mass", SDF_FLOAT, unit_mass, \
     "Ndim", SDF_INT, Ndim, \
     "seed", SDF_INT, seed, \
     "Nmesh", SDF_INT, Nmesh, \
     "data_level", SDF_INT, data_level, \
     "boundary_difference", SDF_INT, boundary_difference, \
     "ic_center_x", SDF_INT, center[0], \
     "ic_center_y", SDF_INT, center[1], \
     "ic_center_z", SDF_INT, center[2], \
     "rcut_min", SDF_FLOAT, rcut_min, \
     "rcut_max", SDF_FLOAT, rcut_max, \
     "nreject", SDF_INT, nreject, \
     "R0", SDF_FLOAT, L0/2

#define TAG 0x147

void
main(int argc, char *argv[])
{
    char cfile[256];
    char outname[256];
    char msgfile[256];
    unsigned long nn[3];
    SDF *csdfp;
    float *data;
    int i;
    int npts, gnpts, initial_gnpts;
    MPMYFile *fp;
    int write_sdf, write_mass_only, write_image, write_header;
    int Nmesh;
    int seed;
    int do_complex_fft;
    float L0, Omega0, Lambda_prime, T_cmb, Tquad, h_100;
    float Zinitial;
    float t2, t3, t4, t5;
    float nu, a, b, c;
    float pt_target;
    float Z;
    int EBW_cdmspec;
    int do_powlaw;
    int negmass_ok;
    float powlaw_n, norm;
    void (*Pkcdm)(int i, int j, int k, float *real, float *imag);
    MPMY_Comm_request req;
    float sum, sumsq, min, avg, var, rms, range;
    float rcut_min, rcut_max;
    char ICdesc[256];
    char ICparams[256];
    int data_level, boundary_difference;
    int center[3] = {0.,0.,0.};
    struct max_st max_s;
    int do_center, do_offset;

    MPMY_Init(&argc, &argv);
#ifndef __DELTA__
    sprintf(msgfile, "msgs/msg.%d", MPMY_Procnum());
    MsgdirInit(msgfile);
    Msg_turnon("SDFwrite.c");
#endif

    EnableTimer(&Tot, "Total");
    EnableWCTimer(&TotWC, "Total(WC)");
    ClearEnabledTimers();
    ClearEnabledCounters();
    StartTimer(&Tot);
    StartTimer(&TotWC);

    if (argc > 1)
      strncpy(cfile, argv[1], sizeof(cfile));
    else
      Getsparam("control file", cfile);
    
    if ((csdfp = SDFopen(NULL, cfile)) == NULL) {
 	SinglError("Sorry, couldn't SDFopen %s\n%s\n",
	      cfile, SDFerrstring);
    }
    singlPrintf("cfile \"%s\" opened\n", cfile);

    SDFgetstringOrDie(csdfp, "outfile", outname, sizeof(outname));
    SDFgetintOrDefault(csdfp, "write_sdf", &write_sdf, 1);
    SDFgetintOrDefault(csdfp, "write_header", &write_header, 0);
    SDFgetintOrDefault(csdfp, "data_level", &data_level, 0);
    SDFgetintOrDefault(csdfp, "boundary_difference", &boundary_difference, 0);
    SDFgetintOrDefault(csdfp, "do_center", &do_center, 0);
    SDFgetintOrDefault(csdfp, "do_offset", &do_offset, 0);
    if (do_offset) {
	SDFgetintOrDie(csdfp, "ic_center_x", &center[0]); 
	SDFgetintOrDie(csdfp, "ic_center_y", &center[1]); 
	SDFgetintOrDie(csdfp, "ic_center_z", &center[2]); 
    }
    SDFgetintOrDefault(csdfp, "write_mass_only", &write_mass_only, 0);
    SDFgetintOrDefault(csdfp, "Nmesh", &Nmesh, 128);
    SDFgetintOrDefault(csdfp, "seed", &seed, 123);
    SDFgetintOrDefault(csdfp, "negmass_ok", &negmass_ok, 0);
    SDFgetintOrDefault(csdfp, "do_complex_fft", &do_complex_fft, 0);
    SDFgetintOrDefault(csdfp, "Ndim", &Ndim, 3);
    SDFgetintOrDefault(csdfp, "write_image", &write_image, (Ndim==2) ? 1 : 0);
    SDFgetfloatOrDefault(csdfp, "rcut_min", &rcut_min, 0.0);
    SDFgetfloatOrDefault(csdfp, "L0", &L0, 250000.0);
    SDFgetfloatOrDefault(csdfp, "rcut_max", &rcut_max, L0);
    SDFgetfloatOrDefault(csdfp, "Zinitial", &Zinitial, 0.0);
    SDFgetfloatOrDefault(csdfp, "Omega0", &Omega0, 1.0);
    SDFgetfloatOrDefault(csdfp, "Lambda_prime", &Lambda_prime, 0.0);
    SDFgetfloatOrDefault(csdfp, "h_100", &h_100, 0.5);
    SDFgetfloatOrDefault(csdfp, "T_cmb", &T_cmb, 2.735);
    SDFgetfloatOrDefault(csdfp, "Tquad", &Tquad, 15.3e-6);
    SDFgetfloatOrDefault(csdfp, "Z", &Z, -1.0);
    /* Defaults are  the h=0.5, Omega=1, Omega_baryon=0.05 Holtzman model. */
    SDFgetintOrDefault(csdfp, "EBW_cdmspec", &EBW_cdmspec, 0);
    SDFgetintOrDefault(csdfp, "do_powlaw", &do_powlaw, 0);
    if (EBW_cdmspec) {
	SDFgetfloatOrDefault(csdfp, "EBW_nu", &nu, 1.13);
	SDFgetfloatOrDefault(csdfp, "EBW_a", &a, 6.4);
	SDFgetfloatOrDefault(csdfp, "EBW_b", &b, 3.0);
	SDFgetfloatOrDefault(csdfp, "EBW_c", &c, 1.7);
	Pkcdm = EBW_cdm;
	sprintf(ICdesc, "EBW cdm");
	sprintf(ICparams, "nu=%.2f, a=%.2f, b=%.2f, c=%.2f", nu, a, b, c);
    } else if (do_powlaw) {
	SDFgetfloatOrDefault(csdfp, "powlaw_n", &powlaw_n, 1.0);
	SDFgetfloatOrDefault(csdfp, "norm", &norm, 1.0);
	Pkcdm = powlaw;
	sprintf(ICdesc, "Power Law");
	sprintf(ICparams, "n = %f, norm = %f", powlaw_n, norm);
    } else {
	SDFgetfloatOrDefault(csdfp, "Holtzman_t2", &t2, -.9876);
	SDFgetfloatOrDefault(csdfp, "Holtzman_t3", &t3, 26.27);
	SDFgetfloatOrDefault(csdfp, "Holtzman_t4", &t4, 43.51);
	SDFgetfloatOrDefault(csdfp, "Holtzman_t5", &t5, 50.45);
	Pkcdm = Holtzman_cdm;
	sprintf(ICdesc, "Holtzman cdm");
	sprintf(ICparams, "t2=%.4f, t3=%.4f, t4=%.4f, t5=%.4f",
		t2, t3, t4, t5);
    }
    
    if (Z != -1.0) 
      SDFgetfloatOrDefault(csdfp, "pt_target", &pt_target, 0.2);

    if (Ndim == 3) {
	nn[0] = nn[1] = nn[2] = Nmesh;
    } else if (Ndim == 2) {
	nn[1] = 1;
	nn[0] = nn[2] = Nmesh;
    }
    gnpts = initial_gnpts = nn[0]*nn[1]*nn[2];
    npts = gnpts/MPMY_Nproc();

    /* The idea is to give the same results on any nproc < nn[0] */
    ranp_setup(seed, 3 * 2*nn[0]*nn[1]*(nn[2]+1), nn[0], &Ranstate);
    if (EBW_cdmspec) 
      setup_EBW_cdm(L0, Omega0, h_100, T_cmb, Tquad, 
		    nu, a, b, c, &Ranstate);      
    else if (do_powlaw)
      setup_powlaw(L0, Omega0, h_100, norm, Tquad, powlaw_n, &Ranstate);      
    else 
      setup_Holtzman_cdm(L0, Omega0, h_100, T_cmb, Tquad, 
			 t2, t3, t4, t5, &Ranstate);

    data = Calloc(npts, sizeof(float));

    sinc_setup(nn);
    if (do_complex_fft)		
      ft(nn, data, &Ranstate, Pkcdm); /* This one does not work in parallel */
    else
      prft(nn, data, &Ranstate, Pkcdm);

    if (do_offset) {
	/* Use given offset */
	center_data(nn, data, center);
    }

    if (data_level) 
      shrink(&data, &gnpts, &npts, nn, data_level);
    
    sum = 0.;
    sumsq = 0.;
    max_s.max = min = 0.;
    for (i = 0; i < npts; i++) {
	if(max_s.max < data[i]) {
	    max_s.max = data[i];
	    max_s.index = i + npts*MPMY_Procnum();
	}
	if(min > data[i]) min = data[i];
	sum += data[i];
	sumsq += data[i]*data[i];
    }
    MPMY_ICombine_Init(&req);
    MPMY_ICombine(&sum, &sum, 1, MPMY_FLOAT, MPMY_SUM, req);
    MPMY_ICombine(&sumsq, &sumsq, 1, MPMY_FLOAT, MPMY_SUM, req);
    MPMY_ICombine(&min, &min, 1, MPMY_FLOAT, MPMY_MIN, req);
    MPMY_ICombine_func(&max_s, &max_s, sizeof(struct max_st), maxidx_f, req);
    MPMY_ICombine_Wait(req);

    /* We 'know' these are zero-mean variables, so we don't subtract */
    /* the average. */
    avg = sum/gnpts;
    var = sumsq/gnpts;
    rms = sqrt(var);
    range = max_s.max-min;

    singlPrintf("max %f, min %f, avg %f, rms %f\n", max_s.max, min, avg, rms);

    if (do_center) {
	/* Offset so highest point is at center */
	center[0] = (max_s.index / (nn[0]*nn[1]));
	center[1] = ((max_s.index / nn[2]) % nn[1]);
	center[2] = (max_s.index % nn[2]);
	center_data(nn, data, center);
    } 
    
    if (write_image) {
	char imgname[256];
	int i;
	float scale;
	char *image;
	sprintf(imgname, "%s_img", outname);
	image = Malloc(npts);
	scale = (max_s.max > -min) ? 127.0/max_s.max : -127.0/min;
	for (i = 0; i < npts; i++) {
	    image[i] = 128+scale*data[i];
	}
	fp = MPMY_Fopen(imgname, 
			MPMY_WRONLY | MPMY_CREAT | MPMY_TRUNC | MPMY_MULTI);
	MPMY_Fwrite(image, npts, 1, fp);
	MPMY_Fclose(fp);
	Free(image);
    }
    
    if (write_sdf) {
	float M0, H0, T0, H, tpos, unit_mass;

	H0 = h_100*0.1*(one_Gyr/one_kpc); /* in Gyr^-1 */
	/* H = H0*(1.0+Zinitial) * sqrt(1.0 + Omega0*Zinitial); */
	H = hubble_from_Z(Omega0, H0, Lambda_prime, Zinitial);
	M0 = Omega0*(3./(8.*M_PI*GNEWT))*H0*H0*L0*L0*L0;
	unit_mass = M0/gnpts;
	tpos = t_from_Z(Omega0, H0, Lambda_prime, Zinitial);
	T0 = t_from_Z(Omega0, H0, Lambda_prime, 0.0);

	if (write_mass_only) {
	/* We write only the delta-mass field to save memory and disk space */
	/* In this case, it is up to the N-body code to compute the rest */

	    SDFwrite(outname, gnpts, npts, data, sizeof(float), MASSDESC,
		     "ICdesc", SDF_STRING, ICdesc,
		     "ICparams", SDF_STRING, ICparams,
		     "npart", SDF_INT, gnpts,
		     "regular_mesh", SDF_INT, 1,
		     "Gnewt", SDF_FLOAT, GNEWT,
		     "Omega0", SDF_FLOAT, Omega0,
		     "Lambda_prime", SDF_FLOAT, Lambda_prime,
		     "Zinitial", SDF_FLOAT, Zinitial,
		     "redshift", SDF_FLOAT, Zinitial,
		     "tpos", SDF_FLOAT, tpos,
		     "H0", SDF_FLOAT, H0,
		     "hubble", SDF_FLOAT, H,
		     "h_100", SDF_FLOAT, h_100,
		     "Tquad", SDF_FLOAT, Tquad,
		     "unit_mass", SDF_FLOAT, unit_mass,
		     "Ndim", SDF_INT, Ndim,
		     "seed", SDF_INT, seed,
		     "Nmesh", SDF_INT, Nmesh,
		     "R0", SDF_FLOAT, L0/2,
		     NULL);
	} else {
	    float midpt[MAXNDIM], posfac[MAXNDIM];
	    float midpt_b[MAXNDIM], posfac_b[MAXNDIM];
	    int bfac;
	    int index[MAXNDIM];
	    float r, rsq;
	    int gnobj;
	    int nobj = 0;
	    int nreject = 0;
	    float massout = 0.0;
	    float rcut_max2, rcut_min2;
	    float scale_fac, growth_fac, velocity_fac;
	    int d;
	    float mass;
	    body *btab;

	    scale_fac = 1.0/(1.0+Zinitial);
	    growth_fac = growthfac_from_Z(Omega0, H0, Lambda_prime, Zinitial)
	      / growthfac_from_Z(Omega0, H0, Lambda_prime, 0.0);
	    velocity_fac = velfac_from_Z(Omega0, H0, Lambda_prime, Zinitial);
	    rcut_min2 = rcut_min*scale_fac; rcut_min2 *= rcut_min2;
	    rcut_max2 = rcut_max*scale_fac; rcut_max2 *= rcut_max2;
	    btab = Malloc(npts * sizeof(body));
	    bfac = 1 << boundary_difference;
	    /* nn has already been reduced by 1 << data_level */
	    assert(nn[0] % bfac == 0);
	    for(d=0; d<Ndim; d++){
		index[d] = 0;
		posfac[d] = scale_fac*L0/nn[d];
		posfac_b[d] = scale_fac*L0/(nn[d]/bfac);
		midpt[d] = (nn[d]-1.)/2.;
		midpt_b[d] = ((nn[d]/bfac)-1.)/2.;
	    }
	    index[0] = MPMY_Procnum()*nn[0]/MPMY_Nproc();
	    for (i = 0; i < npts; i++) {
		if (data[i] != (float)0.0) {
		    mass = unit_mass * (1.0 + data[i]*growth_fac);
		    btab[nobj].mass = mass;
		    rsq = (float)0.0;
		    for (d = 0; d < Ndim; d++) {
			r = posfac[d]*(index[d]-midpt[d]);
			btab[nobj].pos[d] = r;
    /* We need to match the boundary with a possibly coarser outer sphere */
			r = posfac_b[d]*((index[d]/bfac)-midpt_b[d]);
			rsq += r*r;
		    }
		    if (negmass_ok || mass > (float)0.0) {
			if (rsq >= rcut_min2 && rsq < rcut_max2) {
			    massout += mass;
			    nobj++;
			}
		    } else {
			nreject++;
		    }
		}
		d = Ndim-1;
		while (d >= 0){
		    if (++index[d] == nn[d]) {
			index[d--] = 0;
		    } else {
			break;
		    }
		}
	    }
	    MPMY_ICombine_Init(&req);
	    MPMY_ICombine(&nobj, &gnobj, 1, MPMY_INT, MPMY_SUM, req);
	    MPMY_ICombine(&massout, &massout, 1, MPMY_FLOAT, MPMY_SUM, req);
	    MPMY_ICombine(&nreject, &nreject, 1, MPMY_INT, MPMY_SUM, req);
	    MPMY_ICombine_Wait(req);
	    singlPrintf("Fluctuations will grow by %.2f\n", 1.0/growth_fac);
	    singlPrintf("Initial point fluctuations are %f\n", growth_fac*rms);
	    singlPrintf("Initial velocity factor (f) is %f\n", velocity_fac);
	    singlPrintf("T0 is %.2f Gyr\n", T0);
	    singlPrintf("Rejected %d negative mass particles\n", nreject);
	    singlPrintf("Writing %d particles\n", gnobj);

	    if (write_header) {
		char hdrname[256];
		sprintf(hdrname, "%s.hdr", outname);
		/* We don't want npart in the header */
		SDFwritehdr(hdrname, BODYDESC, HDRPARAMS, NULL);
	    }

	    SDFwrite(outname, gnobj, nobj, btab, sizeof(body), BODYDESC,
		     "npart", SDF_INT, gnobj,
		     HDRPARAMS,
		     NULL);

	    Free(btab);
	}
    }

    StopTimer(&Tot);
    StopTimer(&TotWC);

    Free(data);
    Msg_flush();
    OutputTimers(singlPrintf);
    OutputCounters(singlPrintf);
    exit(0);
}

#if 0
/* Crays don't have acosh */
static double Acosh(double x)
{
    return log(x + sqrt(x*x-1.0));
}

static float
growthfac_from_Z(float Omega0, float H0, float Z)
{
    /* This is just the growing mode */
    /* See Weinberg 15.9.27--15.9.31 or Peebles LSS 11.16 */
    double d, d0;

    if (Omega0 == 1.0) {
	d = 1.0/(1.0+Z);
	d0 = 1.0;
    } else if(Omega0 < 1.0) {
	/* Using floats can cause roundoff problems near Omega0=1 */
	double psi, coshpsi;
	coshpsi = 1.0 + 2.0*(1.0 - Omega0)/(Omega0*(1.0+Z));
	psi = Acosh(coshpsi);
	d = - 3.0 * psi * sinh(psi)/((coshpsi-1.0)*(coshpsi-1.0))
	  + (5.0+coshpsi)/(coshpsi-1.0);
	coshpsi = 1.0 + 2.0*(1.0 - Omega0)/Omega0;
	psi = Acosh(coshpsi);
	d0 = - 3.0 * psi * sinh(psi)/((coshpsi-1.0)*(coshpsi-1.0))
	  + (5.0+coshpsi)/(coshpsi-1.0);
    } else {
	double theta, costheta;
	costheta = 1.0 - 2.0*(Omega0-1.0)/(Omega0*(1.0+Z));
	theta = acos(costheta);
	d = - 3.0 * theta * sin(theta)/((1.0-costheta)*(1.0-costheta))
	  + (5.0+costheta)/(1.0-costheta);
	costheta = 1.0 - 2.0*(Omega0-1.0)/Omega0;
	theta = acos(costheta);
	d0 = - 3.0 * theta * sin(theta)/((1.0-costheta)*(1.0-costheta))
	  + (5.0+costheta)/(1.0-costheta);
    }
    return d/d0;
}

static float
t_from_Z(float Omega0, float H0, float Z)
{
    float t, theta, psi;

    if(Omega0 == 1.0){
	t = (2.0/3.0) * pow(1.0+Z, -1.5);
    }else if(Omega0 < 1.0){
	psi = Acosh( 1.0 + 2.0*(1.0 - Omega0)/(Omega0*(1.0+Z)) );
	t = (Omega0/2.0)*pow(1.0-Omega0, -1.5)*(sinh(psi) - psi) ;
    }else{
	theta = acos( 1.0 - 2.0*(Omega0-1.)/(Omega0*(1.0+Z)) );
	t = (Omega0/2.0)*pow(Omega0-1.0, -1.5)*(theta-sin(theta));
    }
    t /= H0;
    return t;
}
#endif

/* This should be the transform of the window function */
/* for the "kernel" that represents each mass point. */
float 
Wsq(int i, int j, int k)
{				/* This is wrong for 2d */
    return (sinctbl[0][i]*sinctbl[1][j]*sinctbl[2][k]);
}

static void 
sinc_setup(unsigned long npts[]){
    int i, j, halfn;
    float k2;

    for(i=0; i<Ndim; i++){
	if(i==0 || npts[i] != npts[i-1]){
	    halfn = (npts[i]/2)+1;
	    sinctbl[i] = Malloc(halfn*sizeof(float));
	    sinctbl[i][0] = 1.;
	    for(j=1; j<halfn; j++){
		k2 = M_PI*j/(npts[i]);
		sinctbl[i][j] = sin(k2)/k2;
	    }
	}else{ 
	    sinctbl[i] = sinctbl[i-1];
	}
    }
}

static long *ranp_seeds;
static int ranp_ntaps;
static int ranp_n;
static int ranp_available;
static int ranp_active;
static int ranp_maxnproc;

void
ranp_setup(long seed, long n, long maxnproc, ran_state *rs)
{
    int i;

    assert(maxnproc % MPMY_Nproc() == 0);
    assert(n % MPMY_Nproc() == 0);
    ranp_ntaps = maxnproc/MPMY_Nproc();
    assert(ranp_ntaps);
    ranp_seeds = Malloc(ranp_ntaps * sizeof(long));

    if (seed <= 0) Error("Bad seed (%ld)\n", seed);
    for (i = 0; i < ranp_ntaps; i++)
      ranp_seeds[i] = -seed - i - MPMY_Procnum()*ranp_ntaps;
    ranp_n = n/maxnproc;
    ranp_active = 0;
    ranp_available = ranp_n;
    ranp_maxnproc = maxnproc;
    
    rs->next_norml_ok = 0;
}

void
ranp_reset(int i, int n, ran_state *rs)
{
    assert(n >= ranp_maxnproc);
    if (i && (i % (n/ranp_maxnproc) == 0)) {
	ranp_available = ranp_n;
	ranp_active++;
	rs->next_norml_ok = 0;
	Msg_do("ranp_reset, active is now %d\n", ranp_active);
    }
}

float
ranp(ran_state *rs)
{
    float ret;
    if (ranp_active >= ranp_ntaps || ranp_active < 0)
      Error("ranp_active too large or no call to ranp_setup\n");
    ret = ran1(ranp_seeds+ranp_active);
    if (--ranp_available == 0) {
	Msg_do("ranp_available reset\n");
	ranp_available = ranp_n;
	ranp_active++;
	rs->next_norml_ok = 0;
    }
    return ret;
}


float normal_rand(ran_state *st)
/*
This is the Polar method for normal distributions, as described on or near
page 104 of Knuth, Semi-numerical Algorithms.  To quote Knuth, "The polar
method is quite slow, but it has essentially perfect accuracy, and it is very
easy to write a program for the polar method..."  'nuf said.  Algorithm due
to Box, Muller and Marsaglia.
*/
{
    float	v1, v2;	/* uniformly distributed on [-1, 1) */
    float s;	/* radius of a point pulled from a uniform circle */
    float	foo;	/* A useful intermediate value. */
    
    if(st->next_norml_ok){
	st->next_norml_ok = 0;
	return st->next_norml;
    }
    
    do{
	v1 = 2.0F * ranp(st) - 1.0F;
	v2 = 2.0F * ranp(st) - 1.0F;
	s = v1*v1 + v2*v2;
    } while(s >= 1.0F);
    foo = sqrt( -2.0F * log(s)/s);
    st->next_norml_ok = 1;
    st->next_norml = v1*foo;
    return v2*foo;
}


/* Shrink dataset by data_level in each dimension */
/* Replaces input data and modifies nn[] and npts */


static void
shrink(float **datahndl, int *gnpts, int *npts, unsigned long nn[3], int data_level)
{
    int procnum = MPMY_Procnum();
    int nproc = MPMY_Nproc();
    int ip, i, j, k;
    int proc_gather;
    unsigned long ns[3];
    float *data = *datahndl;
    float *data2;
    float mfac;
    
    int ifac = 1 << data_level;
    assert(nn[0] % ifac == 0);
    singlPrintf("Reducing by factor of %d in each dimension.\n", ifac);
    *gnpts /= (ifac*ifac*ifac);
    ns[0] = nn[0]/ifac;
    ns[1] = (nn[1] > 1) ? nn[1]/ifac : 1;
    ns[2] = nn[2]/ifac;
    /* If ns[0] is less than nproc, our data layout must be modified */
    proc_gather = nproc/ns[0];
    if (proc_gather > 1) {
	*npts = ns[1]*ns[2];
    } else {
	*npts = ns[0]*ns[1]*ns[2]/nproc;
    }
    data2 = Calloc(*npts, sizeof(float));
    for (ip=procnum*nn[0]/nproc; ip<(procnum+1)*nn[0]/nproc; ip++) {
	i = ip % (nn[0]/nproc);
	for (j=0;j<nn[1];j++) {
	    for (k=0;k<nn[2];k++) {
		data2[Rindex(i/ifac,j/ifac,k/ifac,ns[1],ns[2])] 
		  += data[Rindex(i,j,k,nn[1],nn[2])];
	    }
	}
    }
    Free(data);
    if (proc_gather > 1) {
	if (procnum % proc_gather == 0) {
	    float *tmpbuf = Malloc(*npts*sizeof(float));
	    for (i = procnum+1; i < procnum+proc_gather; i++) {
		MPMY_recvn(tmpbuf, *npts*sizeof(float), i, TAG);
		for (j = 0; j < *npts; j++) {
		    data2[j] += tmpbuf[j];
		}
	    }
	    Free(tmpbuf);
	} else {
	    int sendproc = procnum - procnum % proc_gather;
	    MPMY_send(data2, *npts*sizeof(float), sendproc, TAG);
	    *npts = 0;
	    Free(data2);
	    data2 = 0;
	}
    } 
    mfac = 1.0/(ifac*ifac*ifac);
    for (i = 0; i < *npts; i++)
      data2[i] *= mfac;

    *datahndl = data2;
    nn[0] = ns[0]; nn[1] = ns[1]; nn[2] = ns[2];
}

/* Take grid data and set mass to zero outside some domain */
/* We assume the output filter does not write zero mass particles */

static void
radius_cut(unsigned long nn[3], float *data, float rmin, float rmax)
{
    float midpt[MAXNDIM], posfac[MAXNDIM];
    int index[MAXNDIM];
    float r, rsq;
    int i, d;
    int npts = nn[0]*nn[1]*nn[2]/MPMY_Nproc();
    float min_rsq = rmin*rmin;
    float max_rsq = rmax*rmax;
    int nd=3;

    singlPrintf("Radius cut %f %f\n", rmin, rmax);
    for(d=0; d<nd; d++){
	posfac[d] = 1.0/nn[d];
	index[d] = 0;
	midpt[d] = (nn[d]-1.)/2.;
    }
    assert(nn[0] % MPMY_Nproc() == 0);
    index[0] = MPMY_Procnum()*nn[0]/MPMY_Nproc();
    for (i = 0; i < npts; i++) {
	rsq = (float)0.0;
	for (d = 0; d < nd; d++) {
	    r = posfac[d]*(index[d]-midpt[d]);
	    rsq += r*r;
	}

	if (rsq >= max_rsq || rsq < min_rsq)
	  data[i] = (float)0.0;

	d = nd-1;
	while (d >= 0){
	    if (++index[d] == nn[d]) {
		index[d--] = 0;
	    } else {
		break;
	    }
	}
    }
}

static void
center_data(unsigned long nn[3], float *data, int center[3])
{
    int i, j, k;
    int n, ni, idx, np, nq;
    int nsend1, nsend2;
    float *tmpbuf;
    int offset[3];
    int nproc = MPMY_Nproc();
    int procnum = MPMY_Procnum();
    MPMY_Comm_request req1, req2;

    singlPrintf("putting point (%d,%d,%d) at center\n", 
		center[0], center[1], center[2]);

    offset[0] = (nn[0] + nn[0]/2 - center[0]) % nn[0];
    offset[1] = (nn[1] + nn[1]/2 - center[1]) % nn[1];
    offset[2] = (nn[2] + nn[2]/2 - center[2]) % nn[2];

    singlPrintf("nn (%d,%d,%d) offset (%d,%d,%d)\n", 
		nn[0], nn[1], nn[2], offset[0], offset[1], offset[2]);


    ni = nn[0]/nproc;
    n = ni*nn[1]*nn[2];
    tmpbuf = Malloc(n * sizeof(float));

    /* Do last two on processor dimensions */
    for (i=0;i<nn[0]/nproc;i++) {
	for (j=0;j<nn[1];j++) {
	    for (k=0;k<nn[2];k++) {
		idx = Rindex(i,(j+offset[1])%nn[1],(k+offset[2])%nn[2],
			     nn[1],nn[2]);
		tmpbuf[idx] = data[Rindex(i,j,k,nn[1],nn[2])];
	    }
	}
    }
    for (i = 0; i < n; i++)
      data[i] = 0;
    
    /* Do the across processor dimension */
    /* This code may use send to self, and certainly does for nproc==1 */
    np = offset[0] / ni;
    nq = offset[0] % ni;
    nsend1 = (ni-nq) * nn[1] * nn[2];
    nsend2 = nq * nn[1] * nn[2];
    MPMY_Irecv(data+nsend2, nsend1*sizeof(float), (procnum+nproc-np)%nproc,
	       TAG, &req2);
    MPMY_Isend(tmpbuf, nsend1*sizeof(float), (procnum+np)%nproc, 
	       TAG, &req1);
    MPMY_Wait2(req1, 0, req2, 0);

    MPMY_Irecv(data, nsend2*sizeof(float), (procnum+nproc-np-1)%nproc, 
	       TAG+1, &req2);
    MPMY_Isend(tmpbuf+nsend1, nsend2*sizeof(float), (procnum+np+1)%nproc, 
	       TAG+1,  &req1);
    MPMY_Wait2(req1, 0, req2, 0);

    Free(tmpbuf);
}

/* We need the index that goes with max, and a userfunc is the only way */
void
maxidx_f(const void *a, const void *b, void *c)
{
    const struct max_st *aa, *bb;
    struct max_st *cc;

    aa = a; bb = b; cc = c;
    if (aa->max > bb->max)
      *cc = *aa;
    else
      *cc = *bb;
}
