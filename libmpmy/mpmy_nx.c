/* 
   Implementation of the mpmy interface in terms of nx primitives.
*/
#ifdef __DELTA__
#include <mesh.h>
#endif
#ifdef __GAMMA__
#include <cube.h>
#endif
#ifdef __PARAGON__
#include <nx.h>
#endif
#include <signal.h>
#include <stdlib.h>
#include <errno.h>
#include "bigmalloc.h"
#include "chn.h"
#include "protos.h"
#include "Msgs.h"
#include "mpmy.h"
#include "mpmy_io.h"
#include "mpmy_abnormal.h"
#include "timers.h"
#include "Assert.h"

#ifdef __DELTA__
/* why aren't these in mesh.h???  */
extern void msgwait(int);
extern void msgcancel(int);
extern void gsync(void);
#endif

#define NX_ANY (-1)
#define IN 1
#define OUT 2

/* The default 'process type' is 0 according to man setptype */
/* I have seen 'invalid ptype' errors in msgdone, but I don't know what
   to attribute them to.  - js, 17 dec 94. */
#define PTYPE 0

/* NCOMM: initial allocation of comm structures */
#define NCOMM 800 
struct comm_s{
    int hndl;
    int inout;
    void *buf;			/* not necessary except for debugging */
#ifdef __PARAGON__
    struct {
	long int type;
	long int length;
	long int src;
	long int ptype;
	long int sys[4];} info;
#else
    int src;
#endif
};

static Chn commchn;

int MPMY_Isend(const void *buf, int cnt, int dest, int tag, MPMY_Comm_request *reqp){
    struct comm_s *req = ChnAlloc(&commchn);
    IncrCounter(&MPMYSendCnt);

    if( req == NULL ){
	SeriousWarning("Out of buffers in Isend\n");
	return MPMY_FAILED;
    }
    if( tag == MPMY_TAG_ANY ) tag = -1;
    req->hndl = _isend(tag, (void *)buf, cnt, dest, PTYPE);
    if( req->hndl < 0 ){
	SeriousWarning("isend returns %d, trying to send %d to %d\n", 
		       req->hndl, cnt, dest);
	MPMY_Diagnostic(Msg_do);
	return MPMY_FAILED;
    }

    req->inout = OUT;
    Msgf(("Isend: hndl=%d, dest=%d, tag=%d, len=%d\n", req->hndl, dest, tag, cnt));
    if( Msg_test(__FILE__)){
	int i;
	int sum = 0;
	const char *cbuf = buf;
	for(i=0; i<cnt; i++){
	    sum ^= cbuf[i];
	}
	Msg_do("\tMPMY_Isend cksum: %d\n", sum);
    }
    *reqp = (MPMY_Comm_request)req;
    return MPMY_SUCCESS;
}

int MPMY_Irecv(void *buf, int cnt, int src, int tag, MPMY_Comm_request *reqp){
    struct comm_s *req = ChnAlloc(&commchn);
    IncrCounter(&MPMYRecvCnt);

    if( req == NULL ){
	Shout("Out of buffers in Irecv\n");
	return MPMY_FAILED;
    }

    if( tag == MPMY_TAG_ANY ) tag = NX_ANY;
    if( src == MPMY_SOURCE_ANY ) src = NX_ANY;
#ifdef __PARAGON__
    req->hndl = _irecvx(tag, buf, cnt, src, PTYPE, (long *)&req->info);
#else
    req->hndl = _irecv(tag, buf, cnt);
    req->src = src;
#endif
    if( req->hndl < 0 ){
	Shout("irecv returns %d\n", req->hndl);
	return MPMY_FAILED;
    }

    req->buf = buf;
    req->inout = IN;
    Msgf(("Irecv: hndl=%d, src=%d, tag=%d, len=%d\n", req->hndl, src, tag, cnt));
    *reqp = (MPMY_Comm_request)req;
    return MPMY_SUCCESS;
}

int MPMY_Test(MPMY_Comm_request reqq, int *flag, MPMY_Status *stat){
    struct comm_s *req = reqq;
    int Src, from;

    if( req == NULL ){
	Shout("MPMY_Test(NULL)\n");
	return MPMY_FAILED;
    }

    if( req->hndl < 0 ){
	Shout("Testing handle %d\n", req->hndl);
	return MPMY_FAILED;
    }

    *flag = _msgdone(req->hndl);
    if( *flag < 0 ){
	SeriousWarning("_msgdone(hndl=%d) returns %d in MPMY_Test, errno=%d\n", 
		       req->hndl, *flag, errno);
	MPMY_Diagnostic(Msg_do);
#if 0	/* Extremely wishful thinking ... */
	*flag = 0;
	return MPMY_SUCCESS;
#else
	return MPMY_FAILED;
#endif
    }

    if(*flag){
	IncrCounter(&MPMYDoneCnt);
	if(req->inout == IN){
#ifndef __PARAGON__
/* nx doesn't allow selecting by source.  Thus, we always have no choice
   but to accept data from any source and verify when it arrives that it
   really came from where we expected.  We could try to save it, but
   it would be a lot of work.  Do we need to? */
	    Src = req->src;
	    from = infonode();
	    Msgf(("Test: hndl=%d, Src=%d, from=%d, count=%d, type=%d\n",
		  req->hndl, Src, from, infocount(), infotype()));
	    if( Src!=NX_ANY && Src!=from  ){
		Shout("MPMY_Test.  Message arrived from %d, expecting %d!\n", 
		      from, Src);
		ChnFree(&commchn, req);
		return MPMY_FAILED;
	    }
	    if( stat ){
		stat->count = infocount();
		stat->src = infonode();
		stat->tag = infotype();
	    }
#else
	    if( stat ){
		stat->src = req->info.src;
		stat->tag = req->info.type;
		stat->count = req->info.length;
	    }
	    Msgf(("Recvd(T) hndl=%d from %d, tag %d, count: %d\n", req->hndl,
		  req->info.src, req->info.type, req->info.length));
#endif	
	}
	req->hndl = -req->hndl; /* prevent accidental re-use */
	ChnFree(&commchn, req);
    }
    return MPMY_SUCCESS;
}

#define HAVE_MPMY_WAIT
int MPMY_Wait(MPMY_Comm_request reqq, MPMY_Status *stat){
    struct comm_s *req = reqq;
    int Src;
    int from;

    if( req == NULL ){
	Shout("MPMY_Wait(NULL)\n");
	return MPMY_FAILED;
    }

    if( req->hndl < 0 ){
	Shout("Negative hndl (%d) in MPMY_Wait\n", req->hndl);
	return MPMY_FAILED;
    }

    Msgf(("Wait for %d...\n", req->hndl));
    msgwait(req->hndl);
    IncrCounter(&MPMYDoneCnt);
    if( req->inout == IN ){
#ifndef __PARAGON__
	Src = req->src;
	from = infonode();
	if( Src!=NX_ANY && Src!= from  ){
	    Shout("MPMY_Wait.  Message arrived from %d, expecting %d!\n",
		  from, Src);
	    ChnFree(&commchn, req);
	    return MPMY_FAILED;
	}
	if( stat ){
	    stat->count = infocount();
	    stat->src = infonode();
	    stat->tag = infotype();
	}
#else
	if( stat ){
	    stat->src = req->info.src;
	    stat->tag = req->info.type;
	    stat->count = req->info.length;
	}
	Msgf(("Recvd(W) hndl=%d from %d, tag %d, count: %d, buf: %p\n", req->hndl,
	      req->info.src, req->info.type, req->info.length, req->buf));
	if( Msg_test(__FILE__)){
	    int i;
	    int sum = 0;
	    char *cbuf = req->buf;
	    for(i=0; i<req->info.length; i++){
		sum ^= cbuf[i];
	    }
	    Msg_do("\tMPMY_Wait cksum: %d\n", sum);
	}
#endif
    }
    req->hndl = -req->hndl;
    ChnFree(&commchn, req);
    return MPMY_SUCCESS;
}

#define HAVE_MPMY_SYNC
int MPMY_Sync(void){
    gsync();
    return MPMY_SUCCESS;
}

#ifdef __PARAGON__		/* This MAY NOT BE NECESSARY!!  */
#define HAVE_MPMY_PHYSNODE
static char _Physnode[64];
const char *MPMY_Physnode(void){
  return _Physnode;
}

#define HAVE_MPMY_SHIFT

#include "error.h"
#include "verify.h"

#define SHIFT_TAG 0x1492
/* Because NX can't distinguish different sources when reading */
/* messages, we help it out by adding processor info to the tag. */
int MPMY_Shift(int proc, void *recvbuf, int recvcnt, 
	       const void *sendbuf, int sendcnt, MPMY_Status *stat){
    Msgf(("Starting MPMY_Shift(proc=%d, recvcnt=%d, sendcnt=%d:\n",
	  proc, recvcnt, sendcnt));
    if (proc > MPMY_Procnum()) {
	Verify0(_crecv(SHIFT_TAG+proc, recvbuf, recvcnt));
	Verify0(_csend(SHIFT_TAG+MPMY_Procnum(), (void *)sendbuf, sendcnt, proc, 0));
    } else if (proc == MPMY_Procnum()) {
	memcpy(recvbuf, sendbuf, sendcnt);
    } else {
	Verify0(_csend(SHIFT_TAG+MPMY_Procnum(), (void*)sendbuf, sendcnt, proc, 0));
	Verify0(_crecv(SHIFT_TAG+proc, recvbuf, recvcnt));
    }
    if (stat)
      stat->count = infocount();
    Msgf(("Finished MPMY_Shift\n"));
    
    return MPMY_SUCCESS;
}
#endif /* __PARAGON__ */

#define HAVE_MPMY_FLICK
int MPMY_Flick(void){
    /* Calling flick() seems like a good idea, but it's uses up msec
       (!) on delta and achieves absolutely nothing.  Since we aren't
       multi-threading, there's no value in it whatsoever.  Avoid! */
    /* flick(); */
    return MPMY_SUCCESS;
}

/* The last abnormal handler either falls through or calls exit, 
   depending on the signal being handled... */
static void NXFinalAbFunc(void){
    Msg_do("Node %d in NXFinalAbFunc()\n", MPMY_Procnum()); Msg_flush();
#ifdef __PARAGON__ /* Or anywhere else that implements PrintMemfileVia0 */
    if( MPMY_Procnum() == 0 ){
	/* On node 0, we hang on for 60 sec, hoping that somebody will
	   try to do a PrintMemfileVia0.  If node 0 dies then nobody
	   else can get the 'semaphore' in PrintMemfileVia0.  It's not
	   perfect, but failing to produce output is better than crashing
	   the machine via NORMA... */
	double give_up = dclock() + 60.0;

	/* Don't even think about trying to use signal/alarm! There's
	   nothing else running.  We might as well busy-wait... */
	while( dclock() < give_up ) flick(); /* a reason for flick()? */
    }
#endif

    /* Now what to do?  Shall we fall through, call exit or abort??
       It's extraordinarily hard to get facts about what goes on here
       in the final death-throes of a process.  It has been observed
       that exit(0) is sometimes not sufficient to cause the process
       to actually give up the machine.  The exact circumstances are
       murky.  It may be related to whether one or all processes try
       ot exit, or it may be dependent on what Node0 is doing.  Who
       knows? 

       On the other hand, in some cases, e.g., SEGV, SIGFPE, it might
       be really helpful to fall through and let the OS tell us whatever
       it feels like about the error.  This remains to be explored.  */

    /* As of 2 Jan 1995, it is no longer necessary to reset the
       SIGABRT signal because we don't trap it (and presumably never
       raies it) in mpmy_abnormal.c */
    /* signal(SIGABRT, SIG_DFL); */
    /* But on the otherr hand (Jan 1995), calling abort and dumping
       core is implicated in crashing the paragon.  Let's just try
       to exit as gracefully as we can... */
    /* abort();*/
    exit(1);
}

extern void PrintMemfileVia0(void);

#define HAVE_MPMY_DIAGNOSTIC
void MPMY_Diagnostic(int (*printflike)(const char *, ...)){
    printflike("NX commchn free_cnt = %d, nalloced = %d\n",
	       ChnFreeCnt(&commchn), ChnAllocCnt(&commchn));
}    
    
static void *
NXChnRealloc(void *p, size_t sz){
    Msg("memleak", ("CommChn realloc(%#lx, %ld)\n", 
		    (unsigned long)p, (long)sz));
    Msg("memleak", ("commchn.freecnt=%d, Sendcnt: %d, Recvcnt: %d, Donecnt: %d\n", 
	ChnFreeCnt(&commchn), ReadCounter(&MPMYSendCnt), ReadCounter(&MPMYRecvCnt),
	ReadCounter(&MPMYDoneCnt)));
    return Realloc(p, sz);
}

int MPMY_Init(int *argcp, char ***argvp){
    if( _MPMY_initialized_ ){
	SeriousWarning("MPMY_Init already called!  Not changing state!\n");
	return MPMY_FAILED;
    }
    ChnInit(&commchn, sizeof(struct comm_s), NCOMM, NXChnRealloc);
    /* call sigio_setup() from user code if desired */
    _MPMY_nproc_ = numnodes();
    _MPMY_procnum_ = mynode();
    sprintf(_Physnode, "%d", _myphysnode());
    _MPMY_initialized_ = 1;
    _MPMY_setup_absigs();
    MPMY_OnAbnormal(NXFinalAbFunc);
    MPMY_OnAbnormal(PrintMemfileVia0);
    MPMY_OnAbnormal(PrintMPMYDiags);
    MPMY_OnAbnormal(MPMY_Abannounce);
    return MPMY_SUCCESS;
}

#include "timers_nx.c"
#include "mpmy_io.c"
#include "mpmy_abnormal.c"
#include "mpmy_generic.c"
