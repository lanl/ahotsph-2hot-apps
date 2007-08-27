/* 
   Implementation of the mpmy interface in terms of sunmos primitives.
   (This is the cleanest implementation yet!)
*/
#include <string.h>
#include <stdlib.h>
#include "chn.h"
#include "mpmy.h"
#include "mpmy_abnormal.h"
#include "mpmy_io.h"
#include "bigmalloc.h"
#include "Msgs.h"

/* To control the contents of pario.c, abnormal.c, etc. */
#define __SUNMOS__

/* possible values for inout */
#define IN 2
#define OUT 3

#define NCOMM 300		/* initial allocation of comm_s */
struct comm_s{
    int uflag;
    void *buf;
    void *tmp;
    int cnt;
    int src;
    int tag;
    int inout;
};

static Chn commchn;

/* The (Oct 93) man page for nsend/nrecv says non-blocking receives
   must be aligned on double boundaries.  I'm guessing that this 
   paragon-specific?  Maybe one day I'll get to test this code on an NCUBE*/

#ifdef __INTEL_SSD__
#define MisAligned(p) (((unsigned long int)(p)) & (sizeof(double)-1))
#else
#define MisAligned(p) 0
#endif

int MPMY_Isend(const void *buf, int cnt, int dest, int tag, MPMY_Comm_request *reqp){
    struct comm_s *comm = ChnAlloc(&commchn);
    int ret;

    Msgf(("Isend %d to %d, tag:%d\n", cnt, dest, tag)); Msg_flush();
    if( comm == 0 ){
	return MPMY_FAILED;
    }
    *reqp = comm;
    comm->uflag = 0;
    comm->tmp = NULL;
    ret = _nsend((void *)buf, cnt, dest, tag, &comm->uflag, 0);
    return (ret==0)?MPMY_SUCCESS:MPMY_FAILED;
}

static int junk;
int MPMY_Irecv(void *buf, int cnt, int src, int tag, MPMY_Comm_request *reqp){
    struct comm_s *comm = ChnAlloc(&commchn);
    int ret;

    Msgf(("Irecv %d from %d tag:%d\n", cnt, src, tag)); Msg_flush();
    if( comm == 0 ){
	return MPMY_FAILED;
    }
    *reqp = comm;
    comm->cnt = cnt;
    comm->tag = (tag==MPMY_TAG_ANY)? -1 : tag;
    comm->src = (src==MPMY_SOURCE_ANY)? -1 : src;
    /* According to the Oct '93 man page for nsend/nrecv
     "Unaligned buffers are not allowed with non-blocking receives." */
    if( MisAligned(buf) ){
	/* We have to malloc space and do a memcpy when it arrives! */
	comm->buf = buf;
	buf = comm->tmp = Malloc(cnt);
    }else{
	comm->tmp = NULL;
    }
    comm->uflag = 0;
    ret = _nrecv(buf, &comm->cnt, &comm->src, &comm->tag, &comm->uflag, &junk);
    return (ret==0)?MPMY_SUCCESS:MPMY_FAILED;
}

int MPMY_Test(MPMY_Comm_request req, int *flag, MPMY_Status *stat){
    struct comm_s *comm = req;

    if( comm->uflag ){
	*flag = 1;
	if( comm->inout == IN && stat ){
	    stat->src = comm->src;
	    stat->tag = comm->tag;
	    stat->count = comm->cnt;
	}
	if( comm->tmp ){
	    memcpy(comm->buf, comm->tmp, comm->cnt);
	    Free(comm->tmp);
	}
	ChnFree(&commchn, req);
    }else{
	*flag = 0;
    }
    return MPMY_SUCCESS;
}

/* Sunmos doesn't specify how to figure out who you are.  It's in the */
/* 'compatibilty' library, which depends on the underlying hardware.  Sigh. */
#ifdef __INTEL_SSD__
/* All this just for prototypes formynode and numnodes() ? */
#ifdef __DELTA__
#include <mesh.h>
#endif
#ifdef __GAMMA__
#include <cube.h>
#endif
#ifdef __PARAGON__
#include <nx.h>
#endif

int MPMY_Init(int *argcp, char ***argvp){
    int host, proc, dim;

    _MPMY_procnum_ = mynode();
    _MPMY_nproc_ = numnodes();
    _MPMY_initialized_ = 1;
    ChnInit(&commchn, sizeof(struct comm_s), NCOMM, Realloc_f);
    _MPMY_setup_absigs();
    MPMY_OnAbnormal(MPMY_SystemAbort);
    return MPMY_SUCCESS;
}

/* Now it gets wierd.  Sunmos doesn't have brk or sbrk, so we can't
   link against our home-grown malloc.c.  Cross your fingers: */
#undef malloc
#undef calloc
#undef realloc
#undef free
void *sw_malloc(size_t n){ return malloc(n); }
void *sw_calloc(size_t n, size_t m){ return calloc(n,m); }
void *sw_realloc(void *p, size_t n){ return realloc(p, n);}
void sw_free(void *p){ free(p); }
int malloc_debug(int i){ return -1; }
int malloc_verify(void){return 0;}
size_t malloc_heapsz(void){return -1;}
size_t malloc_avail(void){return -1;}
size_t malloc_used(void){return -1;}
void malloc_print(void){Msg_do("Can't print malloc structures for system malloc\n");}
/* We do the same thing to sigio_dump.c  In fact, sigio_dump should probably
   be part of mpmy and not part of libsw! */
void sigio_setup(void){}

#include "timers_nx.c"
#include "mpmy_pario.c"

#else
/* non-intel==ncube?? */
/* nnodes() and mynode() stopped working in vertex 3.2 ???!!! */

int MPMY_Init(int *argcp, char ***argvp){
    int host, proc, dim;

    whoami(&_MPMY_procnum_, &proc, &host, &dim);
    CommInit();
    _MPMY_nproc_ = 1<<dim;
    _MPMY_initialized_ = 1;
    return MPMY_SUCCESS;
}

#include "timers_posix.c"	/* isn't there an 'ntime'??  */
#include "mpmy_io.c"

#endif

#define HAVE_SYSTEM_ABORT
void MPMY_SystemAbort(void){
    abort();
}

#define NO_SIGNALS
#define CANT_USE_ALARM
#include "mpmy_abnormal.c"
#include "mpmy_generic.c"
