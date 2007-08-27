/* Msg-users note: This file emits messages for both __FILE__ and (somewhat
   more selectively) for "rports" */

/* 
   Implementation of the mpmy interface in terms of CMMD primitives.
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

struct comm_s {
    int deferred;
    CMMD_mcb mcb;
} ;

static int *send_seq_no;
static int *recv_seq_no;
static int debug = 0;


/* Best guess to determine if a request is (or is not) a recv request
   This is only used for debugging.  It doesn't have to be (and is
   not) reliable!! */
#define ISRECV(req) ((int)(req) < 16384)

/* How often to spit out the contents of the queues when msgs are on. */
#define DIAGNOSTIC_FREQ 10000

/* We need to maintain queues of unsent outgoing messages for each
   processor.  WHY ISN'T THIS IN THE OS??? */

struct Qelmt_s{
    int cnt;
    const void *buf;
    int tag;
    struct comm_s *comm;
} ;

#define DELAY 1.0		/* how long (sec) to wait before retry */

static void *Timer;		/* pointer to the MPMY_Timer control */
static Chn QelmtChn;		/* Chain of data in the Dllarr */
static Chn destChn;		/* Chain of data in destDll */
static Dll destDll;		/* destinations that are jammed */
static Dll *Dllarr;		/* a separate Dll for each dest */
static double *failure_time;	/* the time of the most recent jam */
static int NQed;		/* how many messages are jammed in total */
static void InitQ(int nproc);
static int IsClearQ(int dest);
static int RetryQ(void);
static int PushQ(int dest, int tag, const void *buf, int cnt, struct comm_s *req);

static Chn reqchn;

static double Now(void){
    double ret;

    MPMY_StopTimer(Timer);
    ret = MPMY_ReadTimer(Timer);
    MPMY_StartTimer(Timer);
    return ret;
}

int 
MPMY_Isend(const void *buf, int cnt, int dest, int tag, 
	   MPMY_Comm_request *reqp)
{
    struct comm_s *req = ChnAlloc(&reqchn);

    if( req == NULL ){
	SeriousWarning("ChnAlloc returns NULL\n");
	return MPMY_FAILED;
    }

    if( IsClearQ(dest) ){
	req->mcb = CMMD_send_async(dest, tag, buf, cnt, NULL, NULL);
	if( req->mcb < 0 ){
	    if( CMMD_get_errno() == CMMD_ERR_NO_RPORTS ){
		req->deferred = 1;
		failure_time[dest] = Now();
		Msg("rports", ("No rports on %d.\n", dest));
	    }else{
		SeriousWarning("CMMD_send_async failed errno:%d\n",
			       CMMD_get_errno());
		ChnFree(&reqchn, req);
		return MPMY_FAILED;
	    }
	}else{
	    req->deferred = 0;
	}
    }else{
      Msg("rports", ("Q to %d jammed.\n", dest));
      req->deferred = 1;
    }
    if( req->deferred ){
	if( PushQ(dest, tag, buf, cnt, req) != MPMY_SUCCESS ){
	    return MPMY_FAILED;
	}
	req->mcb = -1;
    }
    *reqp = req;
    if( debug )
	Msgf(("Isend(%p,%d)cnt=%d, dest=%d, tag=%d, seq=%d, sum=%#x\n",
	       req, (int)req->mcb, cnt, dest, tag, send_seq_no[dest]++,
	       cksum(buf, cnt)));
    else
	Msgf(("Isend(%p,%d)cnt=%d, dest=%d, tag=%d\n",
	   req, (int)req->mcb, cnt, dest, tag));

    IncrCounter(&MPMYSendCnt);
    return MPMY_SUCCESS;
}

int
MPMY_Irecv(void *buf, int cnt, int src, int tag, MPMY_Comm_request *reqp)
{
    struct comm_s *req;

    req = ChnAlloc(&reqchn);
    if( req == NULL ){
	SeriousWarning("ChnAlloc (MPMY_Irecv) returns NULL\n");
	return MPMY_FAILED;
    }
    if (tag == MPMY_TAG_ANY) tag = CMMD_ANY_TAG;
    if (src == MPMY_SOURCE_ANY) src = CMMD_ANY_NODE;

    req->mcb = CMMD_receive_async(src, tag, buf, cnt, NULL, NULL);
    if( req->mcb < 0 ){
	SeriousWarning("CMMD_recv_async failed, CMMD_get_errno=%d, src=%d, tag=%d, buf=%lx, cnt=%d\n",
		       CMMD_get_errno(), src, tag, (unsigned long)buf, cnt);
	return MPMY_FAILED;
    }
    req->deferred = 0;

    *reqp = req;
    Msgf(("Irecv(%p,%d)src=%d, tag=%d\n", 
	  req, (int)req->mcb, src, tag));
    
    IncrCounter(&MPMYRecvCnt);
    return MPMY_SUCCESS;
}

int
MPMY_Test(MPMY_Comm_request reqq, int *flag, MPMY_Status *stat)
{
    struct comm_s *req = (struct comm_s *)reqq;
    CMMD_mcb mcb;
    static struct comm_s *lastfailedreq;

    if( req == NULL ){
	SeriousWarning("MPMY_Test(NULL): illegal.\n");
	return MPMY_FAILED;
    }
#ifdef AGGRESSIVE_TESTING
    /* NOTE: We have to do this for recv tests too.  Otherwise, we can
       deadlock with both ends spinning on receive Tests.  In fact, I
       think there is still a deadlock problem.  No time now, I've got
       to watch the final scene of Res Dogs again... */
    if( RetryQ() != MPMY_SUCCESS ){
	return MPMY_FAILED;
    }
#endif
    if( req->deferred ){
	Msg("rports-all", ("Test %p: (deferred)\n", req));
	*flag = 0;
	return MPMY_SUCCESS;
    }

    mcb = req->mcb;
    *flag = CMMD_msg_done(mcb);
    if( *flag < 0 ){
	SeriousWarning("CMMD_msg_done returns %d\n", *flag);
	return MPMY_FAILED;
    }
    if (*flag) {
	IncrCounter(&MPMYDoneCnt);
	if (stat) {
	    stat->count = CMMD_mcb_bytes(mcb);
	    stat->src = CMMD_mcb_node(mcb);
	    stat->tag = CMMD_mcb_tag(mcb);
	}
	if( _Msg_enabled && Msg_test(__FILE__)){
	    int cnt = CMMD_mcb_bytes(mcb);
	    int tag = CMMD_mcb_tag(mcb);
	    int srcdest = CMMD_mcb_node(mcb);
	    char info[64];
	    
	    if( debug && ISRECV(mcb) ){
		sprintf(info, ", seq=%d, sum: %#x\n",
			recv_seq_no[CMMD_mcb_node(mcb)]++, 
			cksum(CMMD_mcb_buffer(mcb), cnt));
	    }else{
		strcpy(info, "\n");
	    }
	    
	    if( ISRECV(mcb) ){
		Msg_do("Test(%p,%d) T R cnt=%d, src=%d, tag=%d%s", 
		       req, mcb, cnt, srcdest, tag, info);
	    }else{
		Msg_do("Test(%p,%d) T S cnt=%d, dest=%d, tag=%d%s",
		      req, mcb, cnt, srcdest, tag, info);
	    }
	}
	CMMD_free_mcb(mcb);
	ChnFree(&reqchn, req);
	lastfailedreq = 0;
    } else {
	if( req!=lastfailedreq ) {
	    Msgf(("Test(%p,%d) F\n", req, mcb));
	    lastfailedreq = req;
	}
    }
    return MPMY_SUCCESS;
}

#define HAVE_MPMY_WAIT
int 
MPMY_Wait(MPMY_Comm_request reqq, MPMY_Status *stat)
{
    struct comm_s *req = (struct comm_s *)reqq;
    CMMD_mcb mcb;
    int defcnt = 0;
    int done;

    if( req == NULL ){
	SeriousWarning("MPMY_Wait(NULL): illegal.\n");
	return MPMY_FAILED;
    }

    Msgf(("Wait"));
    Msg_flush();
    while( req->deferred ){
	if( defcnt == 0 )
	    Msg("rports", ("(deferred)"));
	if( RetryQ() == MPMY_FAILED )
	    return MPMY_FAILED;
	if( ++defcnt % 1000000 == 0 ){
	    Warning("MPMY_Wait spun %d times\n", defcnt);
	    MPMY_Diagnostic(Msg_do);
	    Msg_flush();
	}
    }
    mcb = req->mcb;
    Msgf(("(%p,%d), ", req, mcb));
    /* We can't just wait because we might have some outgoing messages that 
       haven't been sent.  */
    do{
	if(NQed){
	    /* Maybe it's done even though others are queued */
	    done = CMMD_msg_done(mcb);
	}else{
	    /* CMMD_msg_wait is a void func.  We can't test for an error... */
	    CMMD_msg_wait(mcb);
	    done = 1;
	}
	if( !done ){
	    /* Might reset NQed */
	    if( RetryQ() == MPMY_FAILED )
		return MPMY_FAILED;
	    /* Looping here does not necessarily indicate that
	       anything is wrong.  We may just be waiting for
	       something that's slow in coming, and we also,
	       coincidentally, have some jammed processors. */
	    if( ++defcnt % 1000000 == 0 ){
		/* Maybe this isn't an error?? */
		if( _Msg_enabled && Msg_test("rports") ){
		    Msg_do("MPMY_Wait spun %d times with NQed=%d\n", 
			   defcnt, NQed);
		    MPMY_Diagnostic(Msg_do);
		    Msg_flush();
		}
	    }
	}
    }while(!done);

    if( stat ){
	stat->count = CMMD_mcb_bytes(mcb);
	stat->src = CMMD_mcb_node(mcb);
	stat->tag = CMMD_mcb_tag(mcb);
    }
    
    if( _Msg_enabled  && Msg_test(__FILE__) ){
	int cnt = CMMD_mcb_bytes(mcb);
	int tag = CMMD_mcb_tag(mcb);
	int srcdest = CMMD_mcb_node(mcb);
	char info[64];
	    
	if( debug && ISRECV(mcb) ){
	    sprintf(info, ", seq=%d, sum: %#x\n",
		    recv_seq_no[CMMD_mcb_node(mcb)]++, 
		    cksum(CMMD_mcb_buffer(mcb), cnt));
	}else{
	    strcpy(info, "\n");
	}
	    
	if( ISRECV(mcb) ){
	    Msg_do("cnt=%d, src=%d, tag=%d%s", 
		  cnt, srcdest, tag, info);
	}else{
	    Msg_do("cnt=%d, dest=%d, tag=%d%s",
		  cnt, srcdest, tag, info);
	}
    }
    CMMD_free_mcb(mcb);
    ChnFree(&reqchn, req);
    IncrCounter(&MPMYDoneCnt);
    return MPMY_SUCCESS;
}

#define HAVE_MPMY_SYNC
int
MPMY_Sync(void)
{
    CMMD_sync_with_nodes();
    return MPMY_SUCCESS;
}

#define HAVE_MPMY_FLICK
int MPMY_Flick(void)
{
    if( RetryQ() != MPMY_SUCCESS )
	return MPMY_FAILED;
    /* I could call CMMD_poll(), but I don't think it would do anything */
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
    double now;

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
	    debug = 1;
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

#if 0 /* Only needed if there are problems REALLY EARLY */
    if( debug ){
	char name[128];
	FILE *fp;
	sprintf(name, "foo.%d", _MPMY_procnum_);
	fp = fopen(name, "w");
	setvbuf(fp, NULL, _IONBF, 0);
	Msg_addfile(fp, (Msgvfprintf_t)vfprintf, (Msgfflush_t)fflush);
	Msg_do("Hello world\n");Msg_flush();
	Msg_on(__FILE__);
	Msg_on("rports");
    }
#endif
    Timer = MPMY_CreateTimer(MPMY_WC_TIME);
    failure_time = Malloc(MPMY_Nproc()*sizeof(double));
    now = Now();
    for(i=0; i<MPMY_Nproc(); i++){
	failure_time[i] = now - (DELAY + 0.00001);
    }

    ChnInit(&reqchn, sizeof(struct comm_s), _MPMY_nproc_, Realloc_f);
    InitQ(_MPMY_nproc_);
    _MPMY_setup_absigs();
    MPMY_OnAbnormal(MPMY_SystemAbort);
    MPMY_OnAbnormal((Abhndlr)Msg_flush);
    MPMY_OnAbnormal(MPMY_Abannounce);
    if( debug ){
	send_seq_no = Calloc(_MPMY_nproc_, sizeof(int));
	recv_seq_no = Calloc(_MPMY_nproc_, sizeof(int));
    }
    /* We probably need another set of functions to deal with SIGINT... */
    signal(SIGINT, sigint);
    signal(SIGTERM, sigterm);
    return MPMY_SUCCESS;
}

#define HAVE_MPMY_DIAGNOSTIC
void MPMY_Diagnostic(void (*printflike)(const char *, ...)){
    struct Qelmt_s *Qelmt;
    Dll_elmt *p;
    int dest;
    Dll *Q;

    (*printflike)("%s list of jammed processors\n", __FILE__); 
    (*printflike)("req cnt\n");
    for(dest=0;dest<MPMY_Nproc(); dest++){
	Q = &Dllarr[dest];
	if( DllLength(Q) )
	    (*printflike)("dest=%d\n", dest);
	for(p = DllTop(Q); p!=DllInf(Q); p=DllDown(p)){
	    Qelmt = DllData(p);
	    (*printflike)("\t%p, %d\n",
			  Qelmt->comm, Qelmt->cnt);
	}
    }
}

static void
InitQ(int nproc){
    int i;

    DllCreateChn(&destChn, sizeof(int), nproc);
    DllCreate(&destDll, &destChn);
    Dllarr = Calloc(nproc, sizeof(Dll));
    /* We really need to create them all to use the same chain... */
    DllCreateChn(&QelmtChn, sizeof(struct Qelmt_s), 10);
    for(i=0; i<nproc; i++){
	DllCreate(&Dllarr[i], &QelmtChn);
    }
}

static int
IsClearQ(int dest){
    return !DllLength(&Dllarr[dest]);
}

/* WARNING!  There is evidence that talking to the I/O system can hang
   when we are low on rports!  Thus, turning on Msgs might be fatal, even
   though the program might extricate itself from its predicament if left
   to its own devices.  One solution might be to use the 'memfile' 
   technology... */
static int
RetryQ(void){
    struct Qelmt_s *Qelmt;
    Dll *Q;
    Dll_elmt *p, *q;
    CMMD_mcb mcb;
    int dest;
    static long int nretry;

    if( _Msg_enabled && Msg_test("rports") && DllLength(&destDll) && (++nretry)%DIAGNOSTIC_FREQ == 0 ){
	Msg_do("RetryQ counter=%ld\n", nretry);
	MPMY_Diagnostic(Msg_do);
    }
    if( DllLength(&destDll) )
	Msg("rports-all", ("Retry Q:\n"));
    /* Loop over the destinations that are pending */
    for(p = DllTop(&destDll); p!=DllInf(&destDll); p=DllDown(p)){
	dest = *(int *)DllData(p);
	if( failure_time[dest] + DELAY > Now() ){
	    Msgf(("Retry to %d canceled.  Need to wait %g more\n",
		   failure_time[dest]+DELAY - Now()));
	    continue;
	}
	Q = &Dllarr[dest];
	/* For each one, loop over all the entries. */
	for(q = DllTop(Q); q!=DllInf(Q); q=DllDown(q)){
	    Qelmt = (struct Qelmt_s *)DllData(q);
	    mcb = CMMD_send_async(dest, Qelmt->tag, Qelmt->buf, 
				  Qelmt->cnt, NULL, NULL);
	    if( mcb < 0 ){
		if( CMMD_get_errno() != CMMD_ERR_NO_RPORTS ){
		    SeriousWarning("Unexpected error in Retry: send(dest=%d, tag=%d, cnt=%d) errno: %d\n",
				   dest, Qelmt->tag, Qelmt->cnt, 
				   CMMD_get_errno());
		    return MPMY_FAILED;
		}
		failure_time[dest] = Now();
		break;
	    }
	    /* Copy the relevant info to the 'comm_s' structure */
	    Qelmt->comm->mcb = mcb;
	    Qelmt->comm->deferred = 0;
	    Msg("rports", ("\tSuccess! req=%p, mcb=%d, dest=%d, tag=%d, cnt=%d.  %d left\n", 
			   Qelmt->comm, mcb, dest, Qelmt->tag, 
			   Qelmt->cnt, DllLength(Q)-1));
	    /* And excise this node from the chain */
	    q = DllDeleteUp(Q, q);
	    NQed--;
	}
	/* If the Q is empty remove the destination from destDll. */
	if( DllLength(Q) == 0 ){
	    p = DllDeleteUp(&destDll, p);
	}
    }
    return MPMY_SUCCESS;
}

static int 
PushQ(int dest, int tag, const void *buf, int cnt, struct comm_s *req){
    struct Qelmt_s *Qelmt;
    Dll *Q = &Dllarr[dest];

    if( DllLength(Q) == 0 ){
	*(int *)DllData(DllInsertAtBottom(&destDll)) = dest;
    }
    Qelmt = DllData(DllInsertAtBottom(Q));
    if( Qelmt == NULL ){
      SeriousWarning("Retry failed to allocate a Qelmt!\n");
      return MPMY_FAILED;
    }
    Msg("rports", ("PushQ req=%p, tag=%d, cnt=%d, dest=%d. %d left\n", 
	  req, tag, cnt, dest, DllLength(Q)));
    Qelmt->tag = tag;
    Qelmt->buf = buf;
    Qelmt->cnt = cnt;
    Qelmt->comm = req;
    NQed++;
    return MPMY_SUCCESS;
}

#include "timers_cm5.c"
#define HAVE_MPMY_TIMERS
#define __CM5__
#define CANT_USE_ALARM
#include "mpmy_pario.c"
#include "mpmy_abnormal.c"
#include "mpmy_generic.c"
