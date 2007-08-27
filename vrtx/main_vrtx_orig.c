/*
Copyright 1992, 1993, 1994, 1995. All Rights Reserved.
Michael S. Warren, John K. Salmon, Gregoire S. Winckelmans
*/
/* 
April 95: modified to produce linear diagnostics and to use RK2 on
first step 
May 95: modified to include the relaxation scheme for div(omega) with
relaxed Jacobi scheme.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "physics_vrtx.h"
#include "SDF.h"
#include "SDFwrite.h"
#include "Assert.h"
#include "protos.h"
#include "vop.h"
#include "timers.h"
#include "Msgs.h"
#include "tree.h"
#include "singlio.h"
#include "files.h"
#include "SDFread.h"
#include "error.h"
#include "getparam.h"
#include "bigmalloc.h"
#include "malloc.h"
#include "mpmy.h"
#include "fastflpt.h"
#include "memfile.h"
#include "abm.h"
#include "neigh.h"

static SDF *startup(int argc, char **argv);
Timer_t StepTot, BuildTot, FindFieldTm;
Timer_t StepWCTm;

float epsilon;
float errtol;
float epsinv;
float nu;

double omega_tot[3], lin_impulse[3], ang_impulse[3], ke, en, he;
          /* shared with GlobalDiags */

static void FixStrengths(body *bp, int nobj);
static void FixVol(body *bp, int nobj);
static int maxmem(void);
static int maxheap(void);

int main(int argc, char *argv[])
{
    bodyptr btab;
    int nobj, gnobj;
    float time;
    char outnamebase[256], restartnamebase[256];
    int do_output, output_freq, save_first;
    int do_restart, restart_freq;
    int do_relaxom, relaxom_freq, nrelaxom, irelaxom;
    float kernel_cutoff;
    int timer_freq;
    SDF *sdfp, *csdfp;
    float rmin[NDIM], rmax[NDIM];
    float dt, dt12, dt32;
    int nsteps;
    int first_step = 1;
    int rk_step, rk_substep;
    int have_velold;
    int iter, laststep;
    int stopnow;
    int xconf, yconf, zconf, strxconf, stryconf, strzconf;
    int volconf, voxconf, voyconf, vozconf;
    int dstrxconf, dstryconf, dstrzconf, identconf, idconf, ntermsconf;
    int omegaxconf, omegayconf, omegazconf;
    tree_t thetree;
    tree_t neightree;
    sortresult_t sortedbtab;
    float sysradius;
    float relax, dtrel, relaxw;
    
    MPMY_Init(&argc, &argv);
    singlPrintf("Welcome to the Vortex integrator\n");
    Msgf(("sizeof(body): %d, offsetof(psix): %d, vx_old: %d, key: %d\n",
		(int)sizeof(body), (int)offsetof(body, psi[0]), 
		(int)offsetof(body,vel_old[0]), (int)offsetof(body, key)));
    csdfp = startup(argc, argv);
    sdfp = SDFread(csdfp, (void **)&btab, &gnobj, &nobj, sizeof(body),
		  "x", offsetof(body, pos[0]), &xconf,
		  "y", offsetof(body, pos[1]), &yconf,
		  "z", offsetof(body, pos[2]), &zconf,
		  "vol", offsetof(body, vol), &volconf,
		  "strx", offsetof(body, strength[0]), &strxconf,
		  "stry", offsetof(body, strength[1]), &stryconf,
		  "strz", offsetof(body, strength[2]), &strzconf,
		   "omegax", offsetof(body, strength[0]), &omegaxconf,
		   "omegay", offsetof(body, strength[1]), &omegayconf,
		   "omegaz", offsetof(body, strength[2]), &omegazconf,
		  "vx_old", offsetof(body, vel_old[0]), &voxconf,
		  "vy_old", offsetof(body, vel_old[1]), &voyconf,
		  "vz_old", offsetof(body, vel_old[2]), &vozconf,
		  "dstrx_old", offsetof(body, dstr_old[0]), &dstrxconf,
		  "dstry_old", offsetof(body, dstr_old[1]), &dstryconf,
		  "dstrz_old", offsetof(body, dstr_old[2]), &dstrzconf,
		  "ident", offsetof(body, ident), &identconf,
		  "id", offsetof(body, ident), &idconf,
		  "nterms", offsetof(body, nterms), &ntermsconf,
		  NULL);
	
    if( xconf==0 || yconf==0 || zconf==0 )
	SinglError("Data doesn't have coordinates x,y,z!\n");

    if( (strxconf && omegaxconf) || (stryconf && omegayconf) 
       || (strzconf && omegazconf) )
	SinglError("Data has both strength and omega.  Too confusing\n");

    if( strxconf != stryconf || strxconf != strzconf )
	SinglError("Inconsistency in finding all components of strength\n");

    if( omegaxconf != omegayconf || omegaxconf != omegazconf )
	SinglError("Inconsistency in finding all components of omega\n");

    if( strxconf==0 && omegaxconf==0 )
	SinglError("Data doesn't have strengths or omegas!\n");

    if( volconf == 0 )
	FixVol(btab, nobj);

    if( omegaxconf ){
	if( volconf == 0 )
	    SinglWarning("Data has Omega but no vol.  >>Guessing<< vol=1.\n");
	FixStrengths(btab, nobj);
    }

    if( identconf && idconf )
	SinglError("Data has both id and ident!\n");

    if( ntermsconf == 0 )
	FixNterms(btab, nobj);

    if( identconf == 0 && idconf == 0 )
	FixId(btab, nobj, gnobj);

    have_velold = (voxconf && voyconf && vozconf 
		   && dstrxconf && dstryconf && dstrzconf);

    singlPrintf("nproc = %d\n", MPMY_Nproc());
    singlPrintf("gnobj = %d\n", gnobj);
    SDFgetfloatOrDefault(sdfp, "time",  &time, (float)0.0);
    singlPrintf("initial_time = %g\n", time);
    SDFgetintOrDefault  (sdfp, "iter",  &iter, 0);
    singlPrintf("initial_iter = %d\n", iter);

    SDFgetfloatOrDefault(sdfp, "epsilon", &epsilon, -1.0);
    if(sdfp)
	SDFclose(sdfp);
    
    if( epsilon == -1.0 ){
	SDFgetfloatOrDie(csdfp, "epsilon", &epsilon);
    }else{
	singlPrintf("Using epsilon=%g from data file.\n"
	       "Disregarding  epsilon from ctrl file!\n", epsilon);
    }
    singlPrintf("epsilon = %g;\n", epsilon);
    epsinv = 1.0F/epsilon;

    SDFgetfloatOrDie(csdfp, "kernel_cutoff", &kernel_cutoff);
    kc = kernel_cutoff * epsilon;
    kc2 = kc*kc;
    singlPrintf("kernel_cutoff = %g * epsilon = %g\n", kernel_cutoff, kc);
    
    SDFgetfloatOrDie(csdfp, "nu", &nu);
    singlPrintf("nu = %g;\n", nu);
    SDFgetfloatOrDie(csdfp, "errtol", &errtol);
    singlPrintf("error_tolerance = %g;\n", errtol);
    SDFgetfloatOrDie(csdfp, "dt", &dt);
    singlPrintf("dt = %g;\n", dt);
    SDFgetfloatOrDie(csdfp, "relax", &relax); 
    singlPrintf("Pedrizzetti's relax = %g;\n", relax);
    dtrel = relax*dt;
    singlPrintf("dt*relax = %g;\n", dtrel);
    SDFgetfloatOrDie(csdfp, "relaxw", &relaxw);
    singlPrintf("Winck relax = %g;\n", relaxw);
    SDFgetintOrDie(csdfp, "nsteps", &nsteps);
    singlPrintf("nsteps = %d;\n", nsteps);
    SDFgetintOrDefault(csdfp, "save_first", &save_first, 0);
    singlPrintf("save first step: %d\n", save_first);
    SDFgetintOrDefault(csdfp, "do_relaxom", &do_relaxom, 0);
    singlPrintf("do_relaxom: %d\n", do_relaxom);
    if( do_relaxom ){
	SDFgetintOrDie(csdfp, "nrelaxom", &nrelaxom);
	singlPrintf("nrelaxom: %d\n", nrelaxom);
	SDFgetintOrDefault(csdfp, "relaxom_freq", &relaxom_freq, 0);
	singlPrintf("relaxom_freq: %d\n", relaxom_freq);
    }

    SDFgetstringOrDefault(csdfp, "outfile", 
			  outnamebase, sizeof(outnamebase), "");
    do_output = ( strlen(outnamebase) > 0 );
    if( do_output ){
	SDFgetintOrDefault(csdfp, "output_freq", &output_freq, nsteps);
    }
    if( output_freq == 0 ){
	do_output = 0;
    }
    if( do_output ){
	singlPrintf("Output to %s.nnn, every %d steps\n", 
	       outnamebase, output_freq);
    }else{
	singlPrintf("No output.\n");
    }
    SDFgetintOrDefault(csdfp, "timer_freq", &timer_freq, output_freq);
    singlPrintf("timer_freq = %d;\n", timer_freq);

    SDFgetstringOrDefault(csdfp, "restartfile", 
			  restartnamebase, sizeof(restartnamebase), "");
    do_restart = ( strlen(restartnamebase) > 0 );
    if( do_restart ){
	SDFgetintOrDefault(csdfp, "restart_freq", &restart_freq, nsteps);
    }
    if( restart_freq == 0 ){
	do_restart = 0;
    }
    if( do_restart ){
	singlPrintf("Restart to %s.nnn, every %d steps\n", 
	       restartnamebase, restart_freq);
    }else{
	singlPrintf("No restart files.\n");
    }

    if(csdfp) 
	SDFclose(csdfp);
    
    if( iter > 0 && have_velold ){
	rk_step = 0;
        rk_substep = 0;
	singlPrintf("Using Adams Bashforth 2  on first step\n");
    }else{
	rk_step = 1;
        rk_substep = 1;
	singlPrintf("Using Runge-Kutta 2  on first step\n");
    }
    singlPrintf("Ready to start. maxmem: %dK, maxheap: %dK\n", maxmem(), maxheap());
    singlFflush();
    
    dt12 = 0.5F*dt;
    dt32 = 1.5F*dt;

    pqsortsetup(&sortedbtab, btab, nobj, sizeof(body), 0.01F, Realloc_f);
    SetupTree(&thetree, NDIM, 
	      sizeof(body), sizeof(cell), TBODYSZ, sizeof(cofm_data),
	      (pq_keyproto)GetKeyFromStruct, (pq_wgtproto)GetCost,
	      CofmFromDaugh, (cellfromcofm_t)CellFromCofm);
    SetupTree(&neightree, NDIM,
	      sizeof(body), sizeof(ncell), TBODYSZ, sizeof(ncofm),
	      (pq_keyproto)GetKeyFromStruct, (pq_wgtproto)GetCost,
	      NeighCofmFromDaugh, (cellfromcofm_t)NeighCellFromCofm);

    singlPrintf("After SetupTree. maxmem: %dK, maxheap: %dK\n", maxmem(), maxheap());


    /* Add one to nsteps because of human psychology. */
    /* If you ask for 100 steps with outputs every 10, you want */
    /* outputs at 0, 10, 20, ... 100.  That's really 101 steps. */

    nsteps += 1;

    laststep = iter + nsteps;

    for (; iter < laststep; iter++){
	StartTimer(&StepTot);
	StartTimer(&StepWCTm);
	if(Msg_test(__FILE__)){
	    Msg_do("btab[0]: %s\n", PrintBodyContents(btab));
	    Msg_do("btab[%d]: %s\n", nobj-1, PrintBodyContents(btab+nobj-1));
	}
	FindBbox(btab, nobj, rmin, rmax);
	sysradius = 0.5*FixRsize(rmin, rmax);
	FixKeys(btab, nobj, GetKey);
	singlPrintf("system radius: %g\n", sysradius);
	singlPrintf("Building...");
	singlFflush();
	StartTimer(&BuildTot);
	BuildTree(&thetree, &sortedbtab);
	StopTimer(&BuildTot);
	singlPrintf("Done. maxmem: %dK, maxheap: %dK\n", maxmem(), maxheap());
	btab = sortedbtab.data;
	nobj = sortedbtab.nobj;

	singlPrintf("Fields...");
	singlFflush();
	StartTimer(&FindFieldTm);
	Walk(&thetree, &thetree, sizeof(Sink), (macv_t)NlgNMACv, (inherit_t)NlgNInherit);
	     
	StopTimer(&FindFieldTm);
	singlFflush();

	FreeTree(&thetree);
	singlPrintf("Done. maxmem: %dK, maxheap: %dK\n", maxmem(), maxheap());
	    
	StopTimer(&StepWCTm);
	StopTimer(&StepTot);

	singlPrintf("iter=%d, time=%g, WC: %g, CPU: %g\n", 
	    iter, time, ReadTimer(&StepWCTm), ReadTimer(&StepTot));

	GlobalDiags(btab, nobj);
	/* We now have positions, velocities and diagnostics           */
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
	*/
	if (do_restart && 
            (( rk_substep || (!rk_step && iter%restart_freq == 0) ) 
              || stopnow) && 
	    (save_first || !first_step)) {

	    char restartname[256];
	    /* Don't try to sort a restart-file */
	    sprintf(restartname, "%s.%03d", restartnamebase, iter);
	    SDFwrite(restartname, gnobj, 
		   nobj, btab, sizeof(body), WHOLEBODYDESC,
		   "npart", SDF_INT, gnobj,
		   "epsilon", SDF_FLOAT, epsilon,
		   "iter", SDF_INT, iter,
		   "time", SDF_FLOAT, time,
		   NULL);
	    singlPrintf("\nRestart file %s created.\n", restartname);
	}

	if (do_output && 
            ( rk_substep || (!rk_step && iter%output_freq == 0) ) ){

	    char outname[256];
	    int output_nobj = nobj;
	    int i;
	    outbody *output_btab;
	    sortresult_t outputsort;

	    output_btab = Malloc(output_nobj * sizeof(outbody));
	    for(i=0; i<nobj; i++){
		VV(output_btab[i].pos, = btab[i].pos);
		VVS(output_btab[i].omega, = btab[i].strength, / btab[i].vol);
		output_btab[i].vol = btab[i].vol;
		output_btab[i].ident = btab[i].ident;
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
		   "epsilon", SDF_FLOAT, epsilon,
		   "iter", SDF_INT, iter,
		   "time", SDF_FLOAT, time,
                   "omega_tot[0]", SDF_DOUBLE, omega_tot[0],
                   "omega_tot[1]", SDF_DOUBLE, omega_tot[1],
                   "omega_tot[2]", SDF_DOUBLE, omega_tot[2],
                   "lin_impulse[0]", SDF_DOUBLE, lin_impulse[0],
                   "lin_impulse[1]", SDF_DOUBLE, lin_impulse[1],
                   "lin_impulse[2]", SDF_DOUBLE, lin_impulse[2],
                   "ang_impulse[0]", SDF_DOUBLE, ang_impulse[0],
                   "ang_impulse[1]", SDF_DOUBLE, ang_impulse[1],
                   "ang_impulse[2]", SDF_DOUBLE, ang_impulse[2],
		   "ke", SDF_DOUBLE, ke,
		   "en", SDF_DOUBLE, en,
		   "he", SDF_DOUBLE, he,
		   NULL);
	    singlPrintf("\nOutput file %s created.\n", outname);
	}


	
	/* It's time to relax the particle vorticity field  */
	if (do_relaxom && 
            ( rk_substep || (!rk_step && iter%relaxom_freq == 0) ) ){
	    
	    /* 
	       First build up another tree so as to know, for each
	       particle, the list of particles that are within
	       kernel_cutoff*epsilon distance from that particle,
	       i.e., dist2 <=kc2 .  This tree is based on this
	       geometrical criterium only.  It is not based on error
	       estimates for multipole expansions.  Hence it is build
	       once and is then used for all iterations */
	    
	    BuildTree(&neightree, &sortedbtab);
	    btab = sortedbtab.data;
	    nobj = sortedbtab.nobj;
	    
            for(irelaxom=0; irelaxom<nrelaxom; irelaxom++){

        	/* Compute the particle vorticity field, Omegat(me)[i] 
        	   for that iteration using that neightree and omegat.c: */
  
		Walk(&neightree, &neightree, sizeof(nsink),
		     (macv_t)NeighMACv, (inherit_t)NeighInherit);

		/* change the particle strengths for that iteration: 
		   see code in relaxomega.c */

                RelaxOmega(btab, nobj, relaxw);

            }
	    singlPrintf("\n Relaxation of particle vorticity completed\n");

	    /* singlPrintf("Diagnostics after relaxation:\n"); 
	    GlobalDiags(btab, nobj); */
	}



	if (timer_freq && 
            ( rk_substep || (!rk_step && iter%timer_freq == 0) ) ){
	    OutputTimers(singlPrintf);
	    OutputCounters(singlPrintf);
	    singlFflush();
	}

	if( stopnow )
	    break;


	/* Perform time integration. */
	/* Runge-Kutta for the first step, Adams-Bashforth for the rest */

	if(rk_substep){
	    Update(btab, nobj, dt, (float)0.0, dtrel, 1); 
            iter -=1;
            rk_substep = 0;
	    first_step = 0;
        }else{
	    Update(btab, nobj, (rk_step)? dt12 : dt32, -dt12, 
                               (rk_step)? (float)0.0 : dtrel,
                               (rk_step)? 0 : 1);
	    time += dt;
	    rk_step = 0;
        }

	if( Msg_test("memleak") ){
	    Msg_do("Memory map after iteration %d\n", iter);
	    malloc_print();
	    MPMY_Diagnostic(Msg_do);
	}

	ClearEnabledCounters();
	ClearEnabledTimers();
    }
    Msg_flush();
    exit(0);
}

/* Do some fairly generic startup stuff. */
SDF *startup(int argc, char **argv){
    SDF *csdfp;
    char cfile[256];
    char msgdir[256];
    char tmp[256];
    char msg_turn_on[512];
    char *msgbase, *lastslash;
    int Msg_memfile;

    if (argc > 1)
 	strncpy(cfile, argv[1], sizeof(cfile));
    else
 	Getsparam("control file", cfile); /* disregard warning */
    if ((csdfp = SDFopen(NULL, cfile)) == NULL) {
 	SinglError("Sorry, couldn't SDFopen %s\n%s\n",
	      cfile, SDFerrstring);
    }
    singlPrintf("cfile \"%s\" opened\n", cfile);
    SDFgetintOrDefault(csdfp, "Msg_memfile", &Msg_memfile, 0);
    if (Msg_memfile) {
#ifdef __PARAGON__
	sigio_setup();
#endif
	memfile_init(Msg_memfile);
	Msg_addfile(0, (Msgvfprintf_t)memfile_vfprintf, 0);
	singlPrintf("Putting all Msgs in memfile\n");
    } else {
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
    }
    SDFgetstringOrDefault(csdfp, "Msg_turn_on", msg_turn_on, 
			  sizeof(msg_turn_on), "");
    Msg_turnon(msg_turn_on);

    if( Msg_test("bigmalloc.c") ){
	malloc_debug(2);
	Msg_do("Malloc_debug(2), expect slow mallocs\n");
    }else{
	malloc_debug(1);
    }

    Msgf(("Messages turned on.  Enabling timers\n")); Msg_flush();
    EnableTimer(&StepTot, "Step Total (CPU)");
    EnableWCTimer(&StepWCTm, "Step Total (WC)");
    EnableTimer(&BuildTot, "Build Total");
    EnableTimer(&FindFieldTm, "Field Comput.");
    EnableTimer(&WalkDeferTm, "Walk Defer");
    EnableCounter(&CellCnt, "Cell Int.");
    EnableCounter(&BodyCnt, "Body Int.");
    EnableCounter(&TaylorKernelCnt, "Kernel(Tayl)");
    EnableCounter(&FullKernelCnt, "Kernel(Full)");
    EnableCounter(&NtermsCnt, "Nterms");
    EnableCounter(&DeferCnt, "Defer");
    EnableCounter(&NfindAcceptsCnt, "Neighbors");
    EnableCounter(&NfindTestsCnt, "Nbr tests");
    if( MPMY_Nproc() > 1 ){
	EnableCounter(&MPMYSendCnt, "MPMY Sends");
	EnableCounter(&MPMYRecvCnt, "MPMY Recvs");
	EnableCounter(&MPMYDoneCnt, "MPMY Done");
	EnableCounter(&ABMByteCnt, "ABM Bytes");
	EnableCounter(&ABMPostCnt, "ABM Posts");
	EnableCounter(&ABMIsendCnt, "ABM Sends");
	ABMHistEnable(3, 13);
    }
    Msgf(("Timers enabled\n"));
    return csdfp;
}

static void 
FixStrengths(body *bp, int nobj){
    int i;

    singlPrintf("Read Omega from file, multiplying through by vol\n");
    for(i=0; i<nobj; i++){
	VS(bp->strength, *= bp->vol);
	bp++;
    }
}

static void 
FixVol(body *bp, int nobj){
    int i;

    singlPrintf("No vol in data file, setting to unity\n");
    for(i=0; i<nobj; i++){
	bp->vol = 1.0F;
	bp++;
    }
}

static int 
maxheap(void)
{
    int memused = malloc_heapsz()/1024;
    MPMY_Combine(&memused, &memused, 1, MPMY_INT, MPMY_MAX);
    return memused;
}

static int 
maxmem(void)
{
    int memused = malloc_used()/1024;
    MPMY_Combine(&memused, &memused, 1, MPMY_INT, MPMY_MAX);
    return memused;
}

/*
From gsw@jedi.gme.usherb.ca  Fri May 19 02:11:31 1995
Return-Path: gsw@jedi.gme.usherb.ca
Received: from jedi.gme.usherb.ca (jedi.gme.usherb.ca [132.210.90.95]) by mullet.anu.edu.au (8.6.11/8.6.9) with SMTP id CAA25456 for <johns@mullet.anu.edu.au>; Fri, 19 May 1995 02:11:27 +1000
Received: by jedi.gme.usherb.ca (AIX 3.2/UCB 5.64/4.03)
          id AA23807; Thu, 18 May 1995 12:08:21 -0400
Date: Thu, 18 May 1995 12:08:21 -0400
From: gsw@jedi.gme.usherb.ca (Gregoire Winckelmans)
Message-Id: <9505181608.AA23807@jedi.gme.usherb.ca>
To: johns@mullet.anu.edu.au
Status: R
*/

