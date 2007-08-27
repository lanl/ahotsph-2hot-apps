/* 
   Implementation of the mpmy interface in terms of CMMD primitives.
   DOES NOT USE send_async
*/

#include <stdio.h>
#include <stdlib.h>
#include <cm/cmmd.h>
#include "gc.h"
#include "Msgs.h"
#include "mpmy.h"
#include "mpmy_time.h"
#include "mpmy_abnormal.h"
#include "Assert.h"
#include "bigmalloc.h"
#include "chn.h"
#include "error.h"
#include "dll.h"

struct comm_s{
    CMMD_mcb mcb;
};

static Chn commchn;
static struct comm_s send_comm;

#define NCOMM 100

int 
MPMY_Isend(const void *buf, int cnt, int dest, int tag, 
	   MPMY_Comm_request *req)
{
    if (CMMD_send_noblock(dest, tag, buf, cnt)) {
	SeriousWarning("send_noblock failed, errno=%d\n", CMMD_get_errno());
	return MPMY_FAILED;
    }
    *req = &send_comm;

    Msgf(("send %d bytes to %d, tag=%d\n", cnt, dest, tag));
    if (CMMD_available_sends() < 10) {
	int i;
	for (i = 0; i < 10000000; i++) {
	    if (CMMD_available_sends() >= 10) return MPMY_SUCCESS;
	}
	SeriousWarning("Low on send rports (%d), delaying failed\n",
		       CMMD_available_sends());
    }
    return MPMY_SUCCESS;
}

int
MPMY_Irecv(void *buf, int cnt, int src, int tag, MPMY_Comm_request *reqp)
{
    struct comm_s *req;

    req = ChnAlloc(&commchn);
    if( req == NULL ){
	SeriousWarning("ChnAlloc (MPMY_Irecv) returns NULL\n");
	return MPMY_FAILED;
    }
    if (tag == MPMY_TAG_ANY) tag = CMMD_ANY_TAG;
    if (src == MPMY_SOURCE_ANY) src = CMMD_ANY_NODE;

    if (CMMD_available_receives() < 10) {
	SeriousWarning("Low on recv rports\n");
    }
    req->mcb = CMMD_receive_async(src, tag, buf, cnt, NULL, NULL);
    if( req->mcb < 0 ){
	SeriousWarning("CMMD_recv_async failed, CMMD_get_errno=%d, src=%d, tag=%d, buf=%lx, cnt=%d\n",
		       CMMD_get_errno(), src, tag, (unsigned long)buf, cnt);
	return MPMY_FAILED;
    }
    *reqp = req;
    Msgf(("Irecv: hndl=%p, src=%d, tag=%d\n", *reqp, src, tag));
    return MPMY_SUCCESS;
}

int
MPMY_Test(MPMY_Comm_request req, int *flag, MPMY_Status *stat)
{
    struct comm_s *comm = req;
    CMMD_mcb mcb;

    if( req == NULL ){
	SeriousWarning("MPMY_Test(NULL): illegal.\n");
	return MPMY_FAILED;
    }
    if( comm == &send_comm ){
	/* It's a send request.  It's done. */
	*flag = 1;
	return MPMY_SUCCESS;
    }
    mcb = comm->mcb;
    *flag = CMMD_msg_done(mcb);
    if( *flag < 0 ){
	SeriousWarning("CMMD_msg_done returns %d\n", *flag);
	return MPMY_FAILED;
    }
    if (*flag) {
	if (stat) {
	    stat->count = CMMD_mcb_bytes(mcb);
	    stat->src = CMMD_mcb_node(mcb);
	    stat->tag = CMMD_mcb_tag(mcb);
	}
	CMMD_free_mcb(mcb);
	ChnFree(&commchn, req);
    }
    return MPMY_SUCCESS;
}

#define HAVE_MPMY_WAIT
int 
MPMY_Wait(MPMY_Comm_request req, MPMY_Status *stat)
{
    struct comm_s *comm = req;
    CMMD_mcb mcb;

    if( req == NULL ){
	SeriousWarning("MPMY_Wait(NULL): illegal.\n");
	return MPMY_FAILED;
    }
    if( comm == &send_comm ){
	/* It's a send request.  It's done. */
	return MPMY_SUCCESS;
    }
    mcb = comm->mcb;
    CMMD_msg_wait(mcb);
    if (stat) {
	stat->count = CMMD_mcb_bytes(mcb);
	stat->src = CMMD_mcb_node(mcb);
	stat->tag = CMMD_mcb_tag(mcb);
    }
    CMMD_free_mcb(mcb);
    ChnFree(&commchn, req);
    return MPMY_SUCCESS;
}

#define HAVE_MPMY_SYNC
int
MPMY_Sync(void)
{
    CMMD_sync_with_nodes();
    return MPMY_SUCCESS;
}

#define HAVE_SYSTEM_ABORT
void MPMY_SystemAbort(void){
    CMMD_error("CMMD_error with signal %d, errno=%d\n", 
	       MPMY_Abnormal_signum, CMMD_get_errno());
    exit(100);
}

/* This should probably be handled via some MPMY construction... */
#include <signal.h>
#include <unistd.h>
static void 
sigint(int sig){
    Msg_do("sigint: Caught signal %d\n", sig);
    MPMY_Diagnostic(Msg_do);
    Msg_flush();
    exit(0);
}

static void 
sigterm(int sig){
    Msg_do("sigterm: Caught signal %d\n", sig);
    MPMY_Diagnostic(Msg_do);
    Msg_flush();
    exit(0);
}

int
MPMY_Init(int *argcp, char ***argvp)
{
    int argc, newargc, i;
    char **argv, **newargv;
#ifdef __CM5VU__
    extern void cm5_alloc_heap(int n);
    extern void __acmain(void);

    /* __acmain() allocates aux stack space */
    /* It gets called automatically if main() is compiled with ac, */
    /* but we don't want to compile main with ac */
    __acmain();
    
    /* aux heap space needs the native malloc, so we need to do it here */
    /* or else sbrk will fail later */

    /* This is defined in asm-cm5/do_grav_fast.c, or else there is an */
    /* empty definition in grav_nv.c if we don't need aux space */
    cm5_alloc_heap(16384);	/* allocates words, not bytes */
#endif

    CMMD_fset_io_mode(stdout,CMMD_independent);
    CMMD_fset_io_mode(stderr,CMMD_independent);
    /* Let's play with argc/argv.  It seems such a shame to have a triply 
       addressed pointer and not use it...
       If we pass --debug_cmmd or --nproc ### as arguments, they are picked up
       by the Init, and DELETED FROM ARGC/ARGV
       */
    argc = *argcp;
    argv = *argvp;
    newargv = Calloc(argc+1, sizeof(char *));
    newargc = 0;
    /* I should use gnu getopt... */
    for(i=0; i<argc; i++){
	if(strcmp(argv[i], "--nproc")==0){
	    /* This lets us set the effective partition size to be
	       smaller than the allocated one.  Possibly useful for
	       debugging, or benchmarking. */
	    _MPMY_nproc_ = atoi(argv[++i]);
	    CMMD_reset_partition_size(_MPMY_nproc_);
	    continue;
	}else if(strcmp(argv[i], "--debug_cmmd")==0){
	    CMMD_enable_safety();
	    continue;
	}
	newargv[newargc++] = argv[i];
    }
    *argvp = newargv;
    *argcp = newargc;
	    
    _MPMY_nproc_ = CMMD_partition_size();
    _MPMY_procnum_ = CMMD_self_address();
    _MPMY_initialized_ = 1;
    /* This can happen if we change the partition size based on
       argc/argv.  We aren't planning on restoring the full size! */
    
    if( _MPMY_procnum_ > _MPMY_nproc_ )
	MPMY_SystemExit();

    ChnInit(&commchn, sizeof(struct comm_s), _MPMY_nproc_, Realloc_f);
    _MPMY_setup_absigs();
    MPMY_OnAbnormal(MPMY_SystemAbort);
    MPMY_OnAbnormal((Abhndlr)Msg_flush);
    MPMY_OnAbnormal(MPMY_Abannounce);
    signal(SIGINT, sigint);
    signal(SIGTERM, sigterm);
    return MPMY_SUCCESS;
}

#include "timers_cm5.c"
#define HAVE_MPMY_TIMERS
#define __CM5__
#define CANT_USE_ALARM
#include "mpmy_io.c"
#include "mpmy_abnormal.c"
#include "mpmy_generic.c"
