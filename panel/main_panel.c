/*
 * Copyright 1992, 1993, 1994  Michael S. Warren, John K. Salmon, and
 * Gregoire S. Winckelmans.  All Rights Reserved.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "fastflpt.h"
#include "physics_panel.h"
#include "SDF.h"
#include "SDFwrite.h"
#include "Assert.h"
#include "protos.h"
#include "vop.h"
#include "timers.h"
#include "Msgs.h"
#include "tree.h"
#include "files.h"
#include "error.h"
#include "getparam.h"
#include "bigmalloc.h"
#include "malloc.h"
#include "mpmy.h"
#include "singlio.h"
#include "SDFread.h"

static SDF *startup(int argc, char **argv);
static SDF *csdfp_s;
Timer_t StepTot, BuildTot, FindFieldTm;


void Uinfty_f(body *bp);
static float U_infty[3];

void Udipole_f(body *bp);
static float dipole_str[3], dipole_pos[3];


void GlobalDiags(int nobj, int gnobj, body *btab);

int main(int argc, char *argv[])
{
    bodyptr btab;
    int nobj, gnobj;
    char outnamebase[256], restartnamebase[256];
    int do_output, output_freq, save_first;
    float relax;
    int do_restart, restart_freq;
    int timer_freq;
    SDF *sdfp, *csdfp;
    float rmin[NDIM], rmax[NDIM];
    int firststep = 1;
    int nsteps;
    int iter, laststep;
    int stopnow;
    int x1conf, y1conf, z1conf;
    int x2conf, y2conf, z2conf;
    int x3conf, y3conf, z3conf;
    int identconf, idconf, ntermsconf;
    tree_t thetree;
    sortresult_t sortedbtab;
    float residual;
    void (*Uexternal)(body *bp);

    singlPrintf("Welcome to the Panel Code\n");
    csdfp_s = csdfp = startup(argc, argv);
    sdfp = SDFread(csdfp, (void **)&btab, &gnobj, &nobj, sizeof(body),
		  "x1", offsetof(body, pos1[0]), &x1conf,
		  "y1", offsetof(body, pos1[1]), &y1conf,
		  "z1", offsetof(body, pos1[2]), &z1conf,
		  "x2", offsetof(body, pos2[0]), &x2conf,
		  "y2", offsetof(body, pos2[1]), &y2conf,
		  "z2", offsetof(body, pos2[2]), &z2conf,
		  "x3", offsetof(body, pos3[0]), &x3conf,
		  "y3", offsetof(body, pos3[1]), &y3conf,
		  "z3", offsetof(body, pos3[2]), &z3conf,
		  "ident", offsetof(body, ident), &identconf,
		  "id", offsetof(body, ident), &idconf,
		  "nterms", offsetof(body, nterms), &ntermsconf,
		  NULL);
	
    if( x1conf==0 || y1conf==0 || z1conf==0 
       || x2conf==0 || y2conf==0 || z2conf==0 
       || x3conf==0 || y3conf==0 || z3conf==0 )
	Error("Data doesn't have coordinates x,y,z!\n");

    if( identconf && idconf )
	Error("Data has both id and ident!\n");

    if( ntermsconf == 0 )
	FixNterms(btab, nobj);

    if( identconf == 0 && idconf == 0 )
	FixId(btab, nobj, gnobj);

    SDFgetintOrDefault  (sdfp, "iter",  &iter, 0);
    if(sdfp)
	SDFclose(sdfp);
    
    SDFgetfloatOrDie(csdfp, "errtol", &errtol);
    SDFgetfloatOrDie(csdfp, "Jacobi_relax", &relax);
    SDFgetintOrDie(csdfp, "nsteps", &nsteps);
    SDFgetintOrDefault(csdfp, "save_first", &save_first, 0);
    SDFgetstringOrDefault(csdfp, "outfile", 
			  outnamebase, sizeof(outnamebase), "");
    do_output = ( strlen(outnamebase) > 0 );
    if( do_output ){
	SDFgetintOrDefault(csdfp, "output_freq", &output_freq, nsteps);
    }
    if( output_freq == 0 ){
	do_output = 0;
    }
    SDFgetintOrDefault(csdfp, "timer_freq", &timer_freq, output_freq);

    SDFgetstringOrDefault(csdfp, "restartfile", 
			  restartnamebase, sizeof(restartnamebase), "");
    do_restart = ( strlen(restartnamebase) > 0 );
    if( do_restart ){
	SDFgetintOrDefault(csdfp, "restart_freq", &restart_freq, nsteps);
    }
    if( restart_freq == 0 ){
	do_restart = 0;
    }

    singlPrintf("errtol = %g;\n", errtol);
    singlPrintf("Jacobi_relax = %g;\n", relax);
    singlPrintf("nsteps = %d;\n", nsteps);
    singlPrintf("gnobj = %d;\n", gnobj);
    if( do_output ){
	singlPrintf("Output to %s.nnn, every %d steps\n", 
	       outnamebase, output_freq);
    }else{
	singlPrintf("No output.\n");
    }
    singlPrintf("int timer_freq = %d;\n", timer_freq);
    if( SDFhasname("U_infty_x", csdfp) ){
	Uexternal = Uinfty_f;
	SDFgetfloatOrDie(csdfp, "U_infty_x", &U_infty[0]);
	SDFgetfloatOrDie(csdfp, "U_infty_y", &U_infty[1]);
	SDFgetfloatOrDie(csdfp, "U_infty_z", &U_infty[2]);
	singlPrintf("Using free-stream velocity (%.2g, %.2g, %.2g)\n",
		    U_infty[0], U_infty[1], U_infty[2]);
    }else{
	Uexternal = Udipole_f;
	SDFgetfloatOrDie(csdfp, "Dipole_x", &dipole_str[0]);
	SDFgetfloatOrDie(csdfp, "Dipole_y", &dipole_str[1]);
	SDFgetfloatOrDie(csdfp, "Dipole_z", &dipole_str[2]);
	SDFgetfloatOrDie(csdfp, "Dipole_pos_x", &dipole_pos[0]);
	SDFgetfloatOrDie(csdfp, "Dipole_pos_y", &dipole_pos[1]);
	SDFgetfloatOrDie(csdfp, "Dipole_pos_z", &dipole_pos[2]);
	singlPrintf("Using point vortex dipole:\n");
	singlPrintf("\tstrength: (%.2g, %.2g, %.2g), location: (%.2g, %.2g, %.2g)\n",
		    dipole_str[0], dipole_str[1], dipole_str[2],
		    dipole_pos[0], dipole_pos[1], dipole_pos[2]);
    }
    
    /* Prepare the panels by precomputing as much as possible */
    PreparePanel(btab, nobj, Uexternal);
    if(csdfp) 
	SDFclose(csdfp);
    
    pqsortsetup(&sortedbtab, btab, nobj, sizeof(body), 0.01F, Realloc_f);
    SetupTree(&thetree, NDIM, 
	      sizeof(body), sizeof(cell), TBODYSZ, sizeof(cofm_data),
	      (pq_keyproto)GetKey, (pq_wgtproto)GetCost,
	      CofmFromDaugh, (cellfromcofm_t)CellFromCofm);

    /* Add one to nstep because of human psychology. */
    /* If you ask for 100 steps with outputs every 10, you want */
    /* outputs at 0, 10, 20, ... 100.  That's really 101 steps. */
    nsteps += 1;
    laststep = iter + nsteps;
    for (; iter < laststep; iter++){
	StartTimer(&StepTot);
	/* Perform Relaxed Jacobi iteration */
	Update(btab, nobj, relax, &residual);

	FindBbox(btab, nobj, rmin, rmax);
	FixRsize(rmin, rmax);
	singlPrintf("Building...");
	StartTimer(&BuildTot);
	BuildTree(&thetree, &sortedbtab);
	StopTimer(&BuildTot);
	singlPrintf("Done\n");
	btab = sortedbtab.data;
	nobj = sortedbtab.nobj;

	singlPrintf("Fields...");
	StartTimer(&FindFieldTm);
	WalkInit(&thetree, &thetree, sizeof(Sink), (macv_t) NlgNMACv, 
		 (inherit_t) NlgNInherit);
	WalkNT(&thetree);
	WalkTerminate();
	StopTimer(&FindFieldTm);
	singlPrintf("Done\n");

	FreeTree(&thetree);

	StopTimer(&StepTot);

	GlobalDiags(nobj, gnobj, btab);

	singlPrintf("iter=%d, CPU: %g, residual(before): %g\n", 
	       iter, ReadTimer(&StepTot), residual);
	    
	/* We now have velocities  */
	/* at the same "time".  We write an output file here, if it is */
	/* the right time. */
	
	if( (stopnow = ForceStop()) )
	    singlPrintf("Emergency STOP!\n");

	/* It's time to write a restart file if:
	   1 - Restart is enabled AND
	   2 - This is not the first step (or we explicitly 
	   want the first step) AND
	   3 - 
	     a) the restart_freq divides the iteration number OR
	     b) we have just been told to stop
	  Simple enough???
	  */
	if (do_restart && 
	    (iter % restart_freq == 0 || stopnow) && 
	    (save_first || !firststep)) {
	    char restartname[256];
	    
	    /* Don't try to sort a restart-file */
	    sprintf(restartname, "%s.%03d", restartnamebase, iter);
	    SDFwrite(restartname, gnobj, 
		   nobj, btab, sizeof(body), WHOLEBODYDESC,
		   "npart", SDF_INT, gnobj,
		   "iter", SDF_INT, iter,
		     "residual", SDF_FLOAT, residual,
		   NULL);
	    singlPrintf("\nRestart file %s created.\n", restartname);
	}
	if (do_output && 
	    (iter % output_freq == 0) &&
	    (save_first || !firststep) ){
	    char outname[256];
	    int output_nobj = nobj;
	    int i;
	    outbody *output_btab;
	    sortresult_t outputsort;

	    output_btab = Malloc(output_nobj * sizeof(outbody));
	    for(i=0; i<nobj; i++){
		bodyptr bp = &btab[i];
		outbody *output_bp = &output_btab[i];

		VVVV(output_bp->pos, = bp->pos1, 
		   + bp->pos2, + bp->pos3);
		VS(output_bp->pos, *= 0.33333333F);
		output_bp->sigma = bp->sigma;
		VVV(output_bp->vel, = bp->vel, + bp->uext);
		/* let's save some space in the output files... */
/*
		VVV(vtot, = bp->uext, + bp->vel);
		VV(output_bp->uext, = bp->uext);
		output_bp->vnorm = Dot(vtot, bp->ez);
		VVV(vtan, = vtot, - output_bp->vnorm*bp->ez);
		output_bp->vtan = sqrtf(Dot(vtan, vtan));
		output_bp->gamma[0] = vtot[1]*bp->ez[2] - vtot[2]*bp->ez[1];
		output_bp->gamma[1] = vtot[2]*bp->ez[0] - vtot[0]*bp->ez[2];
		output_bp->gamma[2] = vtot[0]*bp->ez[1] - vtot[1]*bp->ez[0];
		VS(output_bp->gamma, *= (-1.F/(4.F*3.1415926535F)));
		output_bp->errsum = bp->errsum;
		output_bp->errsum2 = sqrt(bp->errsum2);
*/
		output_bp->ident = bp->ident;
	    }
	    pqsortsetup(&outputsort, output_btab, output_nobj, 
			sizeof(outbody), 0.1F, Realloc_f);
	    output_btab = pqsort(&outputsort, 
				 (pq_wgtproto)UnityCost, 
				 (pq_keyproto)OutIdentKey);
	    output_nobj = outputsort.nobj;
	    Msg("output", ("After pqsort, %d outbodies\n", output_nobj));
	    sprintf(outname, "%s.%03d", outnamebase, iter);
	    SDFwrite(outname, gnobj, 
		   output_nobj, output_btab, sizeof(outbody), OUTBODYDESC,
		   "npart", SDF_INT, gnobj,
		   "iter", SDF_INT, iter,
		     "residual", SDF_FLOAT, residual,
		   NULL);
	    Free(output_btab);
	    singlPrintf("\nOutput file %s created.\n", outname);
	}
	if( timer_freq && iter%timer_freq == 0 ){
	    OutputTimers(singlPrintf);
	    OutputCounters(singlPrintf);
	}
	singlFflush();
	if( stopnow )
	    break;

	if( Msg_test("memleak") ){
	    Msg_do("Memory map after iteration %d\n", iter);
	    malloc_print();
	}
	ClearEnabledCounters();
	ClearEnabledTimers();
	firststep = 0;
    }
    Msg_flush();
    MPMY_Finalize();
    exit(0);
}

#ifdef sun
#include <floatingpoint.h>
#endif

/* Do some fairly generic startup stuff. */
SDF *startup(int argc, char **argv){
    SDF *csdfp;
    char cfile[256];
    char msgdir[256];
    char msg_turn_on[512];
    char tmp[256];
    char *msgbase, *lastslash;

#ifdef sun
    ieee_handler("set", "common", SIGFPE_ABORT);
#endif
    if (argc > 1)
 	strncpy(cfile, argv[1], sizeof(cfile));
    else
 	Getsparam("control file", cfile); /* disregard warning */
    if ((csdfp = SDFopen(NULL, cfile)) == NULL) {
 	Error("Sorry, couldn't SDFopen %s\n%s\n",
	      cfile, SDFerrstring);
    }
    /* Get the msgdir either from:
       argv[2] 
       "msgbase" in csdfp
       misc.argv[0]/msg
       
       We then append .<procnum> to the name
       */
    if( argc > 2 ){
	msgbase = argv[2];
    }else if( SDFgetstring(csdfp, "msgbase", tmp, sizeof(tmp))==0 ){
	msgbase = tmp;
    }else{
	lastslash = strrchr(argv[0], '/');
	if( lastslash ){
	    msgbase = lastslash+1;
	}else{
	    msgbase = argv[0];
	}
	sprintf(tmp, "misc.%s/msg", msgbase);
	msgbase = tmp;
    }	
    sprintf(msgdir, "%s.%d", msgbase, MPMY_Procnum());
    MsgdirInit(msgdir);
    SDFgetstringOrDefault(csdfp, "Msg_turn_on", msg_turn_on, 
			  sizeof(msg_turn_on), "");
    Msg_turnon(msg_turn_on);

    if( Msg_test("bigmalloc.c") ){
	malloc_debug(2);
	Msg_do("Malloc_debug(2), expect slow mallocs\n");
    }else{
	malloc_debug(1);
    }

    EnableTimer(&StepTot, "Step Total");
    EnableTimer(&BuildTot, "Build Total");
    EnableTimer(&FindFieldTm, "Field Comput.");
    EnableCounter(&BodyFullCnt, "Body Full Inter.");
    EnableCounter(&BodyQuadCnt, "Body Quad Inter.");
    EnableCounter(&CellInt, "Cell Inter.");

    return csdfp;
}

/* Compute the externally imposed velocity at the point pos. */
/* It may be as simple as VV(extvel, = Uinfty), or it may */
/* call another tree code to find the vel induced by a set of vortices */
void Uinfty_f(body *bp){
    float ndotUinf;

    VV(bp->uext, = U_infty);
    ndotUinf = Dot(bp->pos, U_infty);
    /* WARNING!! Uexact is only correct if the panels are a unit sphere!!! */
    VVV(bp->uexact, = U_infty, - ndotUinf*bp->pos);
    VS(bp->uexact, *= 1.5F);
}

#if 0
/* Compute uext and uexact at the point, pos, due to a dipole */
/* located anywhere outside the unit-sphere, but with strength vector purely radial. */
/* WARNING: uexact only works for the unit-sphere, and for a dipole strength vector */
/* purely radial, even though it doesn't look like it! .*/
void Udipole_f(body *bp){
    float r[3];
    float r2, r2inv12, r2inv, r2inv32, r2inv52, GammadotR;
    float image_pos[3], image_str[3], radial, transverse[3], 
          b2, b2inv12, b2inv, b2inv32;
    static int image_found = 0;

/*  uext due to the dipole. */
    VVV(r, = bp->pos, - dipole_pos);
    r2 = Dot(r, r);
/*    r2inv12 = recipsqrtf(r2); */
    r2inv12 = 1.F/sqrt(r2);
    r2inv = r2inv12*r2inv12;
    r2inv32 = r2inv*r2inv12;
    r2inv52 = r2inv*r2inv32;
    GammadotR = Dot(dipole_str, r);
    VVV(bp->uext, = -3.F*r2inv52*GammadotR*r, + r2inv32*dipole_str);

    /* The image dipole. */
    if( !image_found ){
	b2 = Dot(dipole_pos, dipole_pos);
/*        b2inv12 = recipsqrtf(b2); */
        b2inv12 = 1.F/sqrt(b2);
        b2inv = b2inv12*b2inv12;
        b2inv32 = b2inv*b2inv12;
	VV(image_pos, = b2inv*dipole_pos);
	/* Figure out the 'radial' and 'transverse' components of the */
	/* dipole strength */
	radial = Dot(dipole_pos, dipole_str)*b2inv;
	VVV(transverse, = dipole_str, - radial*dipole_pos);
	VVV(image_str, = transverse, - radial*dipole_pos);
	VS(image_str, *= b2inv32);
	image_found = 1;
    }

    /* contribution to uexact from the image dipole */
    VVV(r, = bp->pos, - image_pos);
    r2 = Dot(r, r);
/*    r2inv12 = recipsqrtf(r2); */
    r2inv12 = 1.F/sqrt(r2);
    r2inv = r2inv12*r2inv12;
    r2inv32 = r2inv*r2inv12;
    r2inv52 = r2inv*r2inv32;
    GammadotR = Dot(image_str, r);
    VVV(bp->uexact, = -3.F*r2inv52*GammadotR*r, + r2inv32*image_str);

    /* Plus the contribution to uexact from the original dipole... */
    VV(bp->uexact, += bp->uext);
}
#endif



/* Compute uext and uexact at the point, pos, due to a dipole */
/* located anywhere outside the unit-sphere, and with arbitrary strength vector */
/* WARNING: uexact only works for the unit-sphere .*/
void Udipole_f(body *bp){
    double r[3];
    double r2, r2inv12, r2inv, r2inv32, r2inv52, GammadotR;
    double radial, transverse[3];
    double b2, b2inv, b2inv32, gt2, gt2inv12;
    double d2, d2inv12, d2inv, d2inv32, x, y, z, d212, xpd212inv;
    double d2c, d2cinv12, d2cinv, d2cinv32, xc, d2c12, xcpd2c12inv;
    double stuff0, stuff1, stuff2, u[3];
    /* These statics are part of the 'exact' calculation.  We compute */
    /* them once and keep them around. */
    static int image_found = 0, line_dipole = 0 ;
    static double image_pos[3], image_str[3];
    static double b2inv12, gt212; /* used in line-dipole calculation */
    static double ex[3], ey[3], ez[3]; /* used in line-dipole calculation */


    /*  uext due to the dipole ( u=-grad(phi) ). */
    VVV(r, = bp->pos, - dipole_pos);
    r2 = Dot(r, r);
    /*    r2inv12 = recipsqrtf(r2); */
    r2inv12 = 1./sqrt(r2);
    r2inv = r2inv12*r2inv12;
    r2inv32 = r2inv*r2inv12;
    r2inv52 = r2inv*r2inv32;
    GammadotR = Dot(dipole_str, r);
    VVV(bp->uext, = -3.F*r2inv52*GammadotR*r, + r2inv32*dipole_str);

    /* The image system: image dipole + image line dipole if necessary */
    if( !image_found ){
	b2 = Dot(dipole_pos, dipole_pos);
        b2inv12 = 1./sqrt(b2);
        b2inv = b2inv12*b2inv12;
        b2inv32 = b2inv*b2inv12;
	VV(image_pos, = b2inv*dipole_pos);
	/* Figure out the 'radial' and 'transverse' components of the */
	/* dipole strength */
	radial = Dot(dipole_pos, dipole_str)*b2inv;
	VVV(transverse, = dipole_str, - radial*dipole_pos);
	VVV(image_str, = transverse, - radial*dipole_pos);
	VS(image_str, *= b2inv32);
	image_found = 1;

        gt2= Dot(transverse, transverse);
        if (gt2 >= 1.e-10){         
	    /* there is a transverse component to the dipole strength!
	       Hence there will also be a linear line dipole for the
	       image system!  Local coordinate system with ex aligned
	       with dipole_pos, ey aligned with transverse component
	       of dipole strength, and ez=ex X ey.  */
	    line_dipole = 1;
	    gt2inv12 = 1./sqrt(gt2);
	    VV(ex, = b2inv12*dipole_pos);
	    VV(ey, = gt2inv12*transverse);
	    ez[0] = ex[1]*ey[2]-ex[2]*ey[1];
	    ez[1] = ex[2]*ey[0]-ex[0]*ey[2];
	    ez[2] = ex[0]*ey[1]-ex[1]*ey[0];
	    gt212=gt2*gt2inv12;      
	    /* We need to remember gt212, ex, ey, ez and b2inv12 */
	}

    }

    /* contribution to uexact from the original dipole (u=-grad(phi) ). */
    VV(bp->uexact, = bp->uext);

    /* Plus the contribution to uexact from the image dipole (u=-grad(phi) ) */
    VVV(r, = bp->pos, - image_pos);
    r2 = Dot(r, r);
    /*    r2inv12 = recipsqrtf(r2); */
    r2inv12 = 1./sqrt(r2);
    r2inv = r2inv12*r2inv12;
    r2inv32 = r2inv*r2inv12;
    r2inv52 = r2inv*r2inv32;
    GammadotR = Dot(image_str, r);
    VVV(bp->uexact, += -3.*r2inv52*GammadotR*r, + r2inv32*image_str);

    if ( line_dipole ){
	/* Plus the contribution to uexact from the image line
	   dipole. This is first evaluated in the local coordinate
	   system. Then its is transformed back to the absolute
	   coordinate system.  */


	d2= Dot(bp->pos, bp->pos);
	d2inv12=1./sqrt(d2);
	d2inv=d2inv12*d2inv12;
	d2inv32=d2inv*d2inv12;

	x= Dot(bp->pos, ex);
	y= Dot(bp->pos, ey);
	z= Dot(bp->pos, ez);
	d212=d2*d2inv12;
	xpd212inv=1./(x+d212);



	d2c=r2;
	d2cinv12=r2inv12;
	d2cinv=r2inv;
	d2cinv32=r2inv32;

	xc= x-b2inv12;
	d2c12=d2c*d2cinv12;
	xcpd2c12inv=1.F/(xc+d2c12);

    

	stuff0 = b2inv12*( d2cinv12*(1.-(x*xcpd2c12inv))
			  -d2inv12 *(1.-(x*xpd212inv  ))   );

	stuff1 = b2inv12*( d2cinv32*(1.-(x*xcpd2c12inv))
			  -d2inv32 *(1.-(x*xpd212inv  ))
			  -x*d2cinv*(xcpd2c12inv*xcpd2c12inv)
			  +x*d2inv *(xpd212inv  *xpd212inv  )  );

	stuff2 = b2inv12*( xc*d2cinv32*(1.-(x*xcpd2c12inv))
			  -x *d2inv32 *(1.-(x*xpd212inv  ))
			  +(d2cinv12*xcpd2c12inv)*(1.-x*d2cinv12)
			  -(d2inv12 *xpd212inv  )*(1.-x*d2inv12 )  );

	/* phi=-y*stuff0, dphidy=y*y*stuff1-stuff0,
	   dphidz=y*z*stuff1, dphidx= y*stuff2, and u= -gt212*grad(phi) */

	u[0]= -gt212*(y*stuff2);
	u[1]= -gt212*(y*y*stuff1-stuff0);
	u[2]= -gt212*(y*z*stuff1);

	bp->uexact[0] +=u[0]*ex[0]+u[1]*ey[0]+u[2]*ez[0];
	bp->uexact[1] +=u[0]*ex[1]+u[1]*ey[1]+u[2]*ez[1];
	bp->uexact[2] +=u[0]*ex[2]+u[1]*ey[2]+u[2]*ez[2];
   
    

    }


}




void GlobalDiags(int nobj, int gnobj, body *btab){
    body *bp;
    int i;
    double errsum, errsum2;
    float maxerr, maxerr2;
    float vtot[NDIM], v2;
    double lift[NDIM];
    float exacterrmax;
    double exacterr, exacterrsum;
    float uerr[NDIM];
    float exactnormal, exactnormalmax;
    MPMY_Comm_request req;

    errsum = 0.;
    errsum2 = 0.;
    maxerr = 0.;
    maxerr2 = 0.;
    exacterrsum = 0.;
    exacterrmax = 0.;
    exactnormal = 0.;
    exactnormalmax = 0.;
    VS(lift, = 0.);
    for(i=0; i<nobj; i++){
	bp = &btab[i];
	errsum += bp->errsum;
	errsum2 += sqrt(bp->errsum2);
	if( maxerr < bp->errsum ) maxerr = bp->errsum;
	if( maxerr2 < bp->errsum2 ) maxerr2 = bp->errsum2;
	/* Omega?, ??? */
	VVV(vtot, = bp->vel, + bp->uext);
	v2 = 0.5F*Dot(vtot, vtot);
	/* assume that p = -0.5*v^2 based on Bernoulli */
	/* assume that all ez are outward normals. */
	VV(lift, -= v2*bp->ip*bp->ez);
	VVV(uerr, = vtot, - bp->uexact);
	exacterr = sqrt(Dot(uerr, uerr));
	exacterrsum += exacterr;
	if( exacterrmax < exacterr ) exacterrmax = exacterr;
	exactnormal = Dot(bp->uexact, bp->ez);
	if( exactnormalmax < exactnormal ) exactnormalmax = exactnormal;
    }

    MPMY_ICombine_Init(&req);
    MPMY_ICombine(&errsum, &errsum, 1, MPMY_DOUBLE, MPMY_SUM, req);
    MPMY_ICombine(&errsum2, &errsum2, 1, MPMY_DOUBLE, MPMY_SUM, req);
    MPMY_ICombine(&lift, &lift, NDIM, MPMY_DOUBLE, MPMY_SUM, req);
    MPMY_ICombine(&exacterrsum, &exacterrsum, 1, MPMY_DOUBLE, MPMY_SUM, req);
    MPMY_ICombine(&maxerr, &maxerr, 1, MPMY_FLOAT, MPMY_MAX, req);
    MPMY_ICombine(&maxerr2, &maxerr2, 1, MPMY_FLOAT, MPMY_MAX, req);
    MPMY_ICombine(&exacterrmax, &exacterrmax, 1, MPMY_FLOAT, MPMY_MAX, req);
    MPMY_ICombine(&exactnormalmax, &exactnormalmax, 1, MPMY_FLOAT, MPMY_MAX, req);
    MPMY_ICombine_Wait(req);

    singlPrintf("lift: (%g %g %g), |lift|=%g\n",
		lift[0], lift[1], lift[2], sqrt(Dot(lift, lift)));
    singlPrintf("velocity error bounds:\n");
    singlPrintf("\tsum error bounds: mean: %g, max: %g\n", errsum/gnobj, maxerr);
    singlPrintf("\tsqrt(sum sq err bounds): mean: %g, max: %g\n",
		errsum2/gnobj, sqrt(maxerr2));
    singlPrintf("\texact errors (%g): mean: %g, max: %g\n",
		exactnormalmax, exacterrsum/gnobj, exacterrmax);
    /* Print the same stuff again, but in a form suitable for cut-and-paste */
    singlPrintf("%d %.3g %.3g %.3g %.3g %.3g %.3g %.3g %.3g\n",
		gnobj, ReadTimer(&StepTot), errsum/gnobj, maxerr, 
		errsum2/gnobj, sqrt(maxerr2), 
		exacterrsum/gnobj, exacterrmax, exactnormalmax);
}
