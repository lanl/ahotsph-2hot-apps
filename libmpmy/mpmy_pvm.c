/* 
   Implementation of the mpmy interface in terms of pvm (non-T3D) primitives.
   This code is based on the vertex implementation.
*/
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <unistd.h>
#include <pvm3.h>
#include "bigmalloc.h"
#include "chn.h"
#include "protos.h"
#include "Msgs.h"
#include "mpmy.h"
#include "mpmy_abnormal.h"
#include "mpmy_io.h"
#include "verify.h"

#define INIT_TYPE 0

static int tids[1024];            /* array of task id, should be dynamic */

struct comm_s{
    void *buf;
    int cnt;
    int src;
    int tag;
};
static Chn commchn;
static struct comm_s send_comm;
#define NCOMM 100

/* This is not a very fast way to do this */
int
pvm_get_PE(int tid)
{
    int i;
    for (i = 0; i < _MPMY_nproc_; i++) {
	if (tid == tids[i]) return i;
    }
    Error("Bad tid in pvm_get_PE\n");
}

int
MPMY_Isend(const void *buf, int cnt, int dest, int tag, 
	   MPMY_Comm_request *req)
{
    int ret;
    ret = pvm_initsend(PvmDataRaw);
    if (ret < 0) Error("pvm_initsend returned %d\n", ret);
    Verify0(pvm_pkbyte((void *)buf, cnt, 1));
    while (pvm_send(tids[dest], tag) < 0);
    Msgf(("pvm_send %d bytes to %d, tag=%d\n", cnt, dest, tag));
    *req = &send_comm;
    return MPMY_SUCCESS;
}

int 
MPMY_Irecv(void *buf, int cnt, int src, int tag, MPMY_Comm_request *reqp)
{
    struct comm_s *comm = ChnAlloc(&commchn);

    if( comm == NULL ){
	return MPMY_FAILED;
    }
    *reqp = comm;
    comm->buf = buf;
    comm->cnt = cnt;
    comm->tag = (tag==MPMY_TAG_ANY)? -1 : tag;
    comm->src = (src==MPMY_SOURCE_ANY)? -1 : src;
    Msgf(("Irecv %d bytes from %d, tag=%d\n", cnt, src, tag));
    return MPMY_SUCCESS;
}

int
MPMY_Test(MPMY_Comm_request req, int *flag, MPMY_Status *stat)
{
    int inbytes;
    int intid, intag, bufid;
    struct comm_s *comm = req;
    int ret;

    if( comm == &send_comm ){
	/* It's a send request.  It's done. */
	*flag = 1;
	return MPMY_SUCCESS;
    }

    if (comm->src < -1 || comm->src > _MPMY_nproc_)
      Error("Bad source in MPMY_Test (%d)\n", comm->src);
    if (comm->src == -1)
      bufid = pvm_nrecv(-1, comm->tag);
    else
      bufid = pvm_nrecv(tids[comm->src], comm->tag);
    if( bufid > 0 ){
	Verify0(pvm_bufinfo(bufid, &inbytes, &intag, &intid));
	Msgf(("pvm_bufinfo(bufid=%d, inbytes=%d, intag=%d, intid=%d)\n", 
	      bufid, inbytes, intag, intid));
	if( inbytes <= comm->cnt ) {
	    ret = pvm_upkbyte(comm->buf, inbytes, 1);
	    if (ret < 0) Error("pvm_upkbyte returned %d\n", ret);
	}		   
	else
	    Error("Incoming message too long, tag=%d, inbytes=%d\n",
		  intag, inbytes);
	ChnFree(&commchn, comm);
	if( stat ){
	    stat->src = pvm_get_PE(intid);
	    stat->tag = intag;
	    stat->count = inbytes;
	}
	*flag = 1;
    }else{
	Msgf(("Msg_Test(src=%d, tag=%d) not ready\n", comm->src, comm->tag));
	*flag = 0;
    }
    return MPMY_SUCCESS;
}

#define HAVE_MPMY_WAIT
int
MPMY_Wait(MPMY_Comm_request req, MPMY_Status *stat)
{
    int inbytes;
    int intid, intag, bufid;
    struct comm_s *comm = req;
    int ret;

    if( comm == &send_comm ){
	/* It's a send request.  It's done. */
	return MPMY_SUCCESS;
    }

    if (comm->src < -1 || comm->src > _MPMY_nproc_)
      Error("Bad source in MPMY_Wait (%d)\n", comm->src);
    if (comm->src == -1)
      bufid = pvm_recv(-1, comm->tag);
    else
      bufid = pvm_recv(tids[comm->src], comm->tag);
    if( bufid < 0 ){
	Error("pvm_recv failed\n");
    } else {
	Verify0(pvm_bufinfo(bufid, &inbytes, &intag, &intid));
	Msgf(("pvm_bufinfo(bufid=%d, inbytes=%d, intag=%d, intid=%d)\n", 
	      bufid, inbytes, intag, intid));
	if( inbytes <= comm->cnt ) {
	    ret = pvm_upkbyte(comm->buf, inbytes, 1);
	    if (ret < 0) Error("pvm_upkbyte returned %d\n", ret);
	}		   
	else
	    Error("Incoming message too long, tag=%d, inbytes=%d\n",
		  intag, inbytes);
	ChnFree(&commchn, comm);
	if( stat ){
	    stat->src = pvm_get_PE(intid);
	    stat->tag = intag;
	    stat->count = inbytes;
	}
    }
    return MPMY_SUCCESS;
}

#ifdef sun
extern int on_exit( void (*func)(int, void *), void *);
static void
on_exit_func(int status, void *arg)
{
    pvm_exit();
}
#else
static void
atexit_func(void)
{
    pvm_exit();
}
#endif

int MPMY_Init(int *argcp, char ***argvp)
{
    int i;
    int mytid;
    char **argv = *argvp;
    char *progname;
    char prog_path[256];
    char *envs;

    /* enroll in pvm */
    if ((mytid = pvm_mytid()) < 0) {
	pvm_perror(0);
	Error("MPMY_Init Failed\n");
    }
    
    /* find out if I am parent or child */
    tids[0] = pvm_parent();
    if (tids[0] == PvmNoParent) {      /* then I am the parent */
	tids[0] = mytid;
	_MPMY_procnum_ = 0;

	if ((envs = getenv("MPMY_DOC")))
	  _MPMY_nproc_ = 1 << atoi(envs);
	else if ((envs = getenv("MPMY_NPROC")))
	  _MPMY_nproc_ = atoi(envs);
	else
	  _MPMY_nproc_ = 1;

	if (_MPMY_nproc_ < 1 || _MPMY_nproc_ > 2048)
	  Error("Bad nproc, (%d)\n", _MPMY_nproc_);
	
	if( (progname = strrchr(argv[0], '/')) != NULL ){
	    progname++;
	} else {
	    progname = argv[0];
	}
	fprintf(stderr, "Starting %d instances of %s\n", _MPMY_nproc_, progname);
	sprintf(prog_path, "%s/%s", getenv("PWD"), progname);
	
	/* start up copies of myself */
	if (_MPMY_nproc_ > 1) {
	    int flags = 0;
	    if (getenv("MPMY_DEBUG"))
	      flags = PvmTaskDebug;
	    /* This will fail before pvm 3.3.3 */
	    i = pvm_spawn(prog_path, argv+1, flags, "", _MPMY_nproc_-1, &tids[1]);
	    if (i < 0) {
		pvm_perror(0);
		Error("pvm_spawn failed\n");
	    } else if (i+1 < _MPMY_nproc_) {
		Error("task %d failed to start, error code %d\n",
		      i+1, tids[i+1]);
	    }
	    /* multicast _MPMY_nproc_ to children */
	    pvm_initsend(PvmDataRaw);
	    pvm_pkint(&_MPMY_nproc_, 1, 1);
	    if (pvm_mcast(&tids[1], _MPMY_nproc_-1, INIT_TYPE)) {
		pvm_perror(0);
		Error("pvm_mcast failed\n");
	    }
	    /* multicast tids array to children */
	    pvm_initsend(PvmDataRaw);
	    pvm_pkint(tids, _MPMY_nproc_, 1);
	    if (pvm_mcast(&tids[1], _MPMY_nproc_-1, INIT_TYPE)) {
		pvm_perror(0);
		Error("pvm_mcast failed\n");
	    }
	}
    } else {   /* I am a child */
	/* receive _MPMY_nproc_ */
	pvm_recv(tids[0], INIT_TYPE);
	pvm_upkint(&_MPMY_nproc_, 1, 1);
	/* receive tids array */
	pvm_recv(tids[0], INIT_TYPE);
	pvm_upkint(tids, _MPMY_nproc_, 1);
	for (i=1; i<_MPMY_nproc_ ; i++) {
	    if (mytid == tids[i]) {
		_MPMY_procnum_ = i; 
		break; 
	    }
	}
	if (i == _MPMY_nproc_) Error("Fell through tid assignment\n");
	/* get working dir from argv[0] and cd there */
	strcpy(prog_path, argv[0]);
	envs = strrchr(prog_path, '/');
	*envs = 0;
	if (chdir(prog_path))
	  Error("Change dir to %s failed\n", prog_path);
    }
    _MPMY_initialized_ = 1;
    ChnInit(&commchn, sizeof(struct comm_s), NCOMM, Realloc_f);
    _MPMY_setup_absigs();
    MPMY_OnAbnormal(MPMY_Abannounce);
#ifdef sun
    on_exit(on_exit_func, NULL);
#else
    atexit(atexit_func);
#endif
    return MPMY_SUCCESS;
}

#include "mpmy_io.c"
#include "mpmy_abnormal.c"
#include "mpmy_generic.c"
