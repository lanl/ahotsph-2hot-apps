#include <stdio.h>		/* only use sprintf */
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stddef.h>
#include <unistd.h>
#include <errno.h>
#include "SDF.h"
#include "malloc.h"
#include "Malloc.h"
#include "SDFwrite.h"
#include "physics.h"
#include "vop.h"
#include "Msgs.h"
#include "singlio.h"
#include "mpmy.h"
#include "mpmy_io.h"
#include "cosmo.h"
#include "integrate.h"
#include "output.h"

int maxmem(void);
int maxheap(void);

void
output(const char *outnamebase, int64_t gnobj, int nobj, const body *btab, int iter, 
       double dt, double dtv, struct cosmo_s *cosmo, 
       double tpos, double tvel, int cosmology, int do_periodic, 
       float eps, float this_eps_scaled, int force_smoothing_type,
       float this_tol, float frac_tol, float frac_tol0, 
       const float *R, const int *N, 
       int write_nfiles, double *ke, double *pe, 
       int do_output, int identsort_output)
{
    int i;
    sortresult_t outputsort;
    outbodyptr output_btab;
    double z, a, h, tacc;
    char outname[256];
    outbodyptr b;
    float rmin[NDIM], rmax[NDIM];
    float sysradius[NDIM];
    MPMY_Comm_request req;
    int checkpoint;

    if (do_output) checkpoint = 0;
    else checkpoint = 1;

    if( Msg_test("memleak") ){
	Msg_do("Memory map before output\n");
	malloc_print();
    }

    Msgf(("Doing output\n"));
    output_btab = Malloc(nobj * sizeof(outbody));
    Msgf(("output_btab Malloc done\n"));
    for (i=0; i<nobj; i++) {
	output_btab[i].mass = btab[i].mass;
	/* pos and vel set in Integrate() */
#ifdef SAVE_ACC
	VV(output_btab[i].acc, = btab[i].acc);
	output_btab[i].phi = btab[i].phi;
#endif
	output_btab[i].ident = btab[i].ident;
    }
    /* Don't sort before Integrate_out or btab and output_btab */
    /* will not be in the same order */
    tacc = tpos;	   /* tpos and tvel advanced in Integrate()  */
    if (cosmology && do_periodic) {
	CosmoIntegrate(&btab[0].mass, &btab[0].pos[0], &btab[0].vel[0],
		       &btab[0].acc[0], &btab[0].phi, sizeof(body),
		       &output_btab[0].pos[0], &output_btab[0].vel[0], sizeof(outbody),
		       nobj, dt, dtv, cosmo, &tpos, &tvel, ke, pe);
    } else {
	Integrate(&btab[0].mass, &btab[0].pos[0], &btab[0].vel[0],
		  &btab[0].acc[0], &btab[0].phi, sizeof(body),
		  &output_btab[0].pos[0], &output_btab[0].vel[0], sizeof(outbody),
		  nobj, dt, dtv, &tpos, &tvel, ke, pe);
    }
    if (do_output && identsort_output) {
	pqsortsetup_order(&outputsort, output_btab, nobj,
			  sizeof(outbody), 0.1, 1, Realloc_f);
	output_btab = pqsort(&outputsort,
			     (pq_wgtproto)UnityCost, 
			     (pq_keyproto)OutIdentKey);
	nobj = outputsort.nobj;
	Msgf(("After output pqsort, %d outbodies\n", nobj));
    }
    if (cosmology) {
	a = Anow(cosmo, tpos);
	z = Znow(cosmo, tpos);
	h = Hnow(cosmo, tpos);
	VV(sysradius, = a*R);
	VV(rmin, = -sysradius);
	VV(rmax, = sysradius);
	ConvertV(&output_btab[0].vel[0], sizeof(outbody), nobj, Anow(cosmo, tpos), 1);
	/* WrapPeriodic */
	for(b=output_btab; b<&output_btab[nobj]; b++) {
	    VVVV(if LPAREN b->pos, > rmax, RPAREN b->pos, -= 2.0*sysradius);
	    VVVV(if LPAREN b->pos, < rmin, RPAREN b->pos, += 2.0*sysradius);
	}
    } else {
	a = 0.0;
	z = 0.0;
	h = 0.0;
    }
    MPMY_ICombine_Init(&req);
    MPMY_ICombine(ke, ke, 1, MPMY_DOUBLE, MPMY_SUM, req);
    MPMY_ICombine(pe, pe, 1, MPMY_DOUBLE, MPMY_SUM, req);
    MPMY_ICombine_Wait(req);

    if (checkpoint) {
	iter++;
	sprintf(outname, "%s.%04d", outnamebase, iter);
    } else if (cosmology) {
	sprintf(outname, "%s_%.03f", outnamebase, a);
    } else {
	sprintf(outname, "%s_t%.03f", outnamebase, tpos);
    }
    singlPrintf("\n%s to %s %d (%d)\n", 
		(checkpoint) ? "Checkpoint" : "Output", outname, maxmem(), maxheap());

    if (write_nfiles) MPMY_Nfileio(1);
    SDFwrite64(outname, gnobj, 
	       nobj, output_btab, sizeof(outbody), OUTBODYDESC,
	       "npart", SDF_INT64, gnobj,
	       "iter", SDF_INT, iter,
	       "tpos", SDF_DOUBLE, tpos,
	       "tvel", SDF_DOUBLE, tvel,
	       "tacc", SDF_DOUBLE, tacc,
	       "R0", SDF_FLOAT, R[2],
	       "redshift", SDF_DOUBLE, z,
	       "a", SDF_DOUBLE, a,
	       "Omega0", SDF_DOUBLE, cosmo->Omega0,
	       "Omega_m", SDF_DOUBLE, cosmo->Omega_m,
	       "Omega_r", SDF_DOUBLE, cosmo->Omega_r,
	       "Omega_de", SDF_DOUBLE, cosmo->Omega_de,
	       "w0", SDF_DOUBLE, cosmo->w0,
	       "wa", SDF_DOUBLE, cosmo->wa,
	       "H0", SDF_DOUBLE, cosmo->H0,
	       "Lambda_prime", SDF_DOUBLE, cosmo->Lambda,
	       "hubble", SDF_DOUBLE, h,
	       "Gnewt", SDF_DOUBLE, cosmo->Gnewt,
	       "tolerance", SDF_FLOAT, this_tol,
	       "frac_tolerance", SDF_FLOAT, frac_tol,
	       "frac_tolerance0", SDF_FLOAT, frac_tol0,
	       "Rx", SDF_FLOAT, R[0],
	       "Ry", SDF_FLOAT, R[1],
	       "Rz", SDF_FLOAT, R[2],
	       "Nx", SDF_INT, N[0],
	       "Ny", SDF_INT, N[1],
	       "Nz", SDF_INT, N[2],
	       "epsilon_mscale", SDF_FLOAT, eps,
	       "epsilon_scaled", SDF_FLOAT, this_eps_scaled,
	       "epsilon0", SDF_FLOAT, this_eps_scaled/a,
	       "force_smoothing_type", SDF_INT, force_smoothing_type,
	       "checkpoint", SDF_INT, checkpoint,
	       "ke", SDF_DOUBLE, *ke,
	       "pe", SDF_DOUBLE, *pe,
	       NULL);
    if (write_nfiles) MPMY_Nfileio(0);
    singlPrintf("%s to %s done %d (%d)\n", 
		(checkpoint) ? "Checkpoint" : "Output", outname, maxmem(), maxheap());
    Msgf(("Output to %s done\n", outname));
    Free(output_btab);
    if (MPMY_Procnum() == 0) {
	char name[256];
	sprintf(name, "%s.restart", outnamebase);
	if (unlink(name))
	    Shout("unlink of %s failed, errno=%d\n", name, errno);
	if (symlink(outname, name))
	    Shout("symlink of %s failed, errno=%d\n", outname, errno);
    }
}
