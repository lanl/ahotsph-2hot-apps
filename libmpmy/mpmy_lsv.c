/* 
   Implementation of the mpmy interface in terms of lsv primitives.
*/
#include <unistd.h>		/* only used for getpid proto in Flick */
#include <signal.h>
#include <stdio.h>		/* only used for sprintf */
#include <string.h>
#include "Assert.h"
#include "Msgs.h"
#include "chn.h"
#include "bigmalloc.h"
#include "protos.h"
#include "mpmy.h"
#include "lsv.h"
#include "mpmy_time.h"
#include "mpmy_io.h"
#include "mpmy_abnormal.h"

#define NCOMM 20
static Chn commchn;

struct comm_s{
    void *buf;
    int cnt;
    int src;
    int tag;
};
static struct comm_s send_comm;

int MPMY_Isend(const void *buf, int cnt, int dest, int tag, MPMY_Comm_request *req){
    IncrCounter(&MPMYSendCnt);
    Ssend(buf, cnt, dest, tag);	/* blocking! */
    *req = &send_comm;
    return MPMY_SUCCESS;
}

int MPMY_Irecv(void *buf, int cnt, int src, int tag, MPMY_Comm_request *reqp){
    struct comm_s *req = ChnAlloc(&commchn);

    IncrCounter(&MPMYRecvCnt);
    if( req == NULL ){
	SeriousWarning("MPMY_Irecv: ChnAlloc failed\n");
	return MPMY_FAILED;
    }
    *reqp = req;
    req->buf = buf;
    req->cnt = cnt;
    req->tag = (tag==MPMY_TAG_ANY)? LSV_ANY : tag;
    req->src = (src==MPMY_SOURCE_ANY)? LSV_ANY : src;
    return MPMY_SUCCESS;
}

int MPMY_Test(MPMY_Comm_request req, int *flag, MPMY_Status *stat){
    struct comm_s *comm = req;
    int inbytes;

    if( comm == &send_comm ){
	/* It's a send request.  It's done. */
	*flag = 1;
	IncrCounter(&MPMYDoneCnt);
	return MPMY_SUCCESS;
    }

    inbytes = Srecv(comm->buf, comm->cnt, 
		    comm->tag, &comm->src);
    if( inbytes >= 0 ){
	if( stat ){
	    stat->src = comm->src;
	    stat->tag = comm->tag;
	    stat->count = inbytes;
	}
	IncrCounter(&MPMYDoneCnt);
	ChnFree(&commchn, req);
	*flag = 1;
    }else{
	*flag = 0;
    }
    return MPMY_SUCCESS;
}

#if 1 /* is this related to the seq-number hanging? */
#define HAVE_MPMY_WAIT
int MPMY_Wait(MPMY_Comm_request req, MPMY_Status *stat){
    struct comm_s *comm = req;
    int inbytes;

    if( comm == &send_comm ){
	/* It's a send request.  It's done. */
	IncrCounter(&MPMYDoneCnt);
	return MPMY_SUCCESS;
    }

    inbytes = Srecv_block(comm->buf, comm->cnt, 
			  comm->tag, &comm->src);
    if( inbytes >= 0 ){
	IncrCounter(&MPMYDoneCnt);
	if( stat ){
	    stat->src = comm->src;
	    stat->tag = comm->tag;
	    stat->count = inbytes;
	}
	ChnFree(&commchn, req);
    }else{
	SeriousWarning("Impossible return (%d) from Srecv_block\n", inbytes);
	return MPMY_FAILED;
    }
    return MPMY_SUCCESS;
}
#endif

#define HAVE_MPMY_FLICK
int MPMY_Flick(void){
    /* What to do?  We'd like to cause this process/thread to
       relinquish the cpu, but we don't want to make it wait for a
       long time.  Talking to the OS is a bad idea.  I'm guessing that
       under Unix, getpid has to query the kernel, and so might let
       some other process run.  Another possibility is to try
       sleep(0), or maybe set a signal handler on some unlikely
       signal.  I dunno... */
/*    (void)getpid(); */
    return MPMY_SUCCESS;
}

#define HAVE_MPMY_DIAGNOSTIC
void MPMY_Diagnostic(int (*printflike)(const char *, ...)){
    Sdiag(printflike);
}

#define HAVE_MPMY_PHYSNODE
static char _Physnode[512];
const char *MPMY_Physnode(void){
  return _Physnode;
}

#if defined(HAVE_UNAME) || 1 /* assume it's ok for now. */
#include <sys/utsname.h>
static void _MPMY_setupphysnode(void){
  struct utsname utsn;
  uname(&utsn);
  strncpy(_Physnode, utsn.nodename, sizeof(_Physnode));
}

#else
static void _MPMY_setupphysnode(void){
  strncpy(_Physnode, "?", sizeof(_Physnode));
}
#endif

static void *
LSVChnRealloc(void *p, size_t sz){
    Msg("memleak", ("CommChn realloc(%p, %ld)\n", p, (long)sz));
    Msg("memleak", ("commchn.freecnt=%d, Sendcnt: %d, Recvcnt: %d, Donecnt: %d\n", 
	commchn.free_cnt, ReadCounter(&MPMYSendCnt), ReadCounter(&MPMYRecvCnt),
	ReadCounter(&MPMYDoneCnt)));
    return Realloc(p, sz);
}

int MPMY_Init(int *argcp, char ***argvp){
  int debug = 0;
  int argc;
  char **argv;
  char **newargv;
  int newargc;
  int i;

  argc = *argcp;
  argv = *argvp;
  newargv = Calloc(argc+1, sizeof(char *));
  newargc = 0;
  /* I should use gnu getopt... */
  for(i=0; i<argc; i++){
    if(strcmp(argv[i], "--debug-lsv")==0){
      debug = 1;
      continue;
    }
    newargv[newargc++] = argv[i];
  }
  *argvp = newargv;
  *argcp = newargc;
	    
  if( debug ){
    char name[128];
    FILE *fp;
    sprintf(name, "lsv-dbg.%d", getpid());
    fp = fopen(name, "w");
    setvbuf(fp, NULL, _IONBF, 0);
    Msg_addfile(fp, (Msgvfprintf_t)vfprintf, (Msgfflush_t)fflush);
    Msg_do("Hello world\n");Msg_flush();
    Msg_on(__FILE__);
    Msg_on("lsv.c");
  }

  Sinit_elt();
  ChnInit(&commchn, sizeof(struct comm_s), NCOMM, LSVChnRealloc);
  _MPMY_nproc_ = LSV_nproc;
  _MPMY_procnum_ = LSV_procnum;
  _MPMY_initialized_ = 1;
  signal(SIGINT, SIG_DFL);	/* why do I need to do this?  If I don't, I don't see them??? */
  _MPMY_setupphysnode();
  _MPMY_setup_absigs();
  MPMY_OnAbnormal(MPMY_SystemAbort);
  /* Notice that we create core dirs in the current directory. */
  sprintf(MPMY_Abchdir_arg, "lsv.core.%d", MPMY_Procnum());
  MPMY_OnAbnormal(MPMY_Abchdir);
  MPMY_OnAbnormal(MPMY_Abannounce);
  return MPMY_SUCCESS;
}

/* #define CANT_USE_ALARM */
#include "mpmy_io.c"
#include "mpmy_abnormal.c"
#include "mpmy_generic.c"
