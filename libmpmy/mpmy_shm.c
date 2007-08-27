/* 
   Implementation of the mpmy interface in terms of pvm (non-T3D) primitives.
   This code is based on the vertex implementation.
*/
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/ipc.h>
#include <sys/shm.h>
#include "bigmalloc.h"
#include "chn.h"
#include "protos.h"
#include "Msgs.h"
#include "mpmy.h"
#include "mpmy_abnormal.h"
#include "mpmy_io.h"
#include "verify.h"

struct comm_s{
    void *buf;
    int cnt;
    int src;
    int tag;
};

struct shm_tag {
    void *buf;
    int cnt;
    int dest;
    int tag;
};

static Chn commchn;
static struct comm_s send_comm;

#define NCOMM 128
#define SHM_PKTSIZE 512
#define SHM_NCOMM 128
int Shmid;
static void *Shm_seg;
static struct shm_tag *Shm_ctl;

static void
init_shm()
{
    key_t key;
    int shmflag, size, shm_offset, slot;
    char *p;
    
    key = 6;
    shm_offset = _MPMY_nproc_ * SHM_NCOMM * sizeof(struct comm_s);
    size = _MPMY_nproc_ * SHM_NCOMM* SHM_PKTSIZE + shm_offset;
    shmflag = IPC_CREAT | 00600; /* read/write permission */
    Shmid = shmget(key, size, shmflag);
    if (Shmid == -1) Error("shmget failed, errno=%d\n", errno);
    Msgf(("shmid is %d\n", Shmid));
    
    p = shmat(Shmid, 0, 0);
    if ((int)p == -1) Error("shmat failed, errno=%d\n", errno);
    Msgf(("Shm_seg starts at %p\n", p));
    Shm_seg = p + shm_offset + _MPMY_procnum_ * SHM_NCOMM * SHM_PKTSIZE;
    Msgf(("my Shm_seg is at %p\n", Shm_seg));
    Shm_ctl = (struct shm_tag *)p;
    Msgf(("Shm_ctl array is at %p\n", Shm_ctl));
    for (slot = 0; slot < SHM_NCOMM; slot++)
      Shm_ctl[_MPMY_procnum_*SHM_NCOMM+slot].dest = -1; /* unlock initially */
}

int
MPMY_Isend(const void *buf, int cnt, int dest, int tag, 
	   MPMY_Comm_request *req)
{
    volatile int done = 0;
    int slot;

    if (cnt > SHM_PKTSIZE*SHM_NCOMM/2) 
      Error("Msg too large (%d), max %d\n", cnt, SHM_PKTSIZE*SHM_NCOMM/2);
    if (cnt > SHM_PKTSIZE) {
	slot = SHM_NCOMM/2;
	do {
	    done = Shm_ctl[_MPMY_procnum_*SHM_NCOMM+slot].dest;
	    if (done != -1) Msgf(("Big slot busy, trying again\n"));
	} while (done != -1);
    } else {
	while (done != -1) {
	    for (slot = 0; slot < SHM_NCOMM; slot++) {
		done = Shm_ctl[_MPMY_procnum_*SHM_NCOMM+slot].dest;
		if (done == -1) break;
	    }
	    if (done != -1) Msgf(("All slots busy, trying again\n"));
	}
    }
    Shm_ctl[_MPMY_procnum_*SHM_NCOMM+slot].buf = Shm_seg+slot*SHM_PKTSIZE;
    Shm_ctl[_MPMY_procnum_*SHM_NCOMM+slot].cnt = cnt;
    Shm_ctl[_MPMY_procnum_*SHM_NCOMM+slot].dest = dest;
    Shm_ctl[_MPMY_procnum_*SHM_NCOMM+slot].tag = tag;
    memcpy(Shm_seg+slot*SHM_PKTSIZE, buf, cnt);
    *req = &send_comm;
    Msgf(("send %d bytes to %d, tag=%d\n", cnt, dest, tag));
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
    comm->tag = tag;
    comm->src = src;
    Msgf(("Irecv %d bytes from %d, tag=%d\n", cnt, src, tag));
    return MPMY_SUCCESS;
}

int
MPMY_Test(MPMY_Comm_request req, int *flag, MPMY_Status *stat)
{
    struct comm_s *comm = req;
    int src, dest, tag, cnt;
    int slot;
    void *buf;

    if( comm == &send_comm ){
	/* It's a send request.  It's done. */
	*flag = 1;
	return MPMY_SUCCESS;
    }

    if (comm->src < MPMY_SOURCE_ANY || comm->src > _MPMY_nproc_) 
      Error("Bad source in MPMY_Test (%d)\n", comm->src);

    for (src = 0; src < _MPMY_nproc_; src++) {
	if (comm->src != MPMY_SOURCE_ANY && comm->src != src) continue;
	for (slot = 0; slot < SHM_NCOMM; slot++) {
	    dest = Shm_ctl[src*SHM_NCOMM+slot].dest;
	    if (dest != _MPMY_procnum_) continue;
	    tag = Shm_ctl[src*SHM_NCOMM+slot].tag;
	    if (tag != comm->tag && comm->tag != MPMY_TAG_ANY) continue;
	    cnt = Shm_ctl[src*SHM_NCOMM+slot].cnt;
	    Msgf(("Recvd(T) from %d, tag %d, count: %d\n", src, tag, cnt));
	    if (cnt > comm->cnt)
	      Error("Cnt too large (%d), max %d\n", cnt, comm->cnt);
	    buf = Shm_ctl[src*SHM_NCOMM+slot].buf;
	    memcpy(comm->buf, buf, cnt);
	    if (stat) {
		stat->src = src;
		stat->tag = tag;
		stat->count = cnt;
	    }
	    ChnFree(&commchn, req);
	    *flag = 1;
	    Shm_ctl[src*SHM_NCOMM+slot].dest = -1; /* unlock buffer */
	    return MPMY_SUCCESS;
	}
    }
    /* Msgf(("Tested out, 0\n")); */
    *flag = 0;
    return MPMY_SUCCESS;
}


#ifdef sun
extern int on_exit( void (*func)(int, void *), void *);
static void
on_exit_func(int status, void *arg)
{
    atexit_func();
}
#endif

static void
atexit_func(void)
{
    int ret;
    ret = shmdt((char *)Shm_ctl);
    Msgf(("shmdt ret is %d\n", ret));

    if (MPMY_Procnum() == 0) {
	ret = shmctl(Shmid, IPC_RMID, 0);
	Msgf(("shmctl ret is %d\n", ret));
    }
}

extern void PrintMemfile(void);

int MPMY_Init(int *argcp, char ***argvp)
{
    char *envs;

    if ((envs = getenv("MPMY_PROCNUM")))
      _MPMY_procnum_ = atoi(envs);
    else
      _MPMY_procnum_ = 0;
    if ((envs = getenv("MPMY_DOC")))
      _MPMY_nproc_ = 1 << atoi(envs);
    else if ((envs = getenv("MPMY_NPROC")))
      _MPMY_nproc_ = atoi(envs);
    else
      _MPMY_nproc_ = 1;

    if (_MPMY_nproc_ < 1 || _MPMY_nproc_ > 2048)
      Error("Bad nproc, (%d)\n", _MPMY_nproc_);

    _MPMY_initialized_ = 1;
    ChnInit(&commchn, sizeof(struct comm_s), NCOMM, Realloc_f);
    _MPMY_setup_absigs();
    MPMY_OnAbnormal(MPMY_SystemAbort);
    /* Notice that we create core dirs in the current directory. */
    sprintf(MPMY_Abchdir_arg, "shm.core.%d", MPMY_Procnum());
    MPMY_OnAbnormal(MPMY_Abchdir);
    MPMY_OnAbnormal(MPMY_Abannounce);
    MPMY_OnAbnormal(PrintMemfile);
#ifdef sun
    on_exit(on_exit_func, NULL);
#else
    atexit(atexit_func);
#endif
    init_shm();
    return MPMY_SUCCESS;
}

#include "mpmy_io.c"
#include "mpmy_abnormal.c"
#include "mpmy_generic.c"
