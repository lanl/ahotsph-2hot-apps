#define HAVE_MPMY_TIMERS
#include <sys/time.h>
#include "mpmy_time.h"
#include "chn.h"
#include "bigmalloc.h"

extern int readrtc(struct timestruc_t *);

static Chn timer_chn;
static int initialized;

typedef struct {
    int type;
    struct timestruc_t start, accum;
} MPMY_Timer;

void *MPMY_CreateTimer(int type){
    MPMY_Timer *ret;

    if( initialized == 0 ){
        ChnInit(&timer_chn, sizeof(MPMY_Timer), 40, Realloc_f);
        initialized = 1;
    }

    ret = ChnAlloc(&timer_chn);
    ret->type = type;
    return (void *)ret;
}

int MPMY_DestroyTimer(void *p){
    ChnFree(&timer_chn, p);
    return MPMY_SUCCESS;
}

int MPMY_ClearTimer(void *p)
{
    MPMY_Timer *t = p;
    ntimerclear(&(t->accum));
    return MPMY_SUCCESS;
}

double MPMY_ReadTimer(void *p)
{
    MPMY_Timer *t = p;
    return (t->accum.tv_sec + t->accum.tv_nsec/(double)NS_PER_SEC);
}

int MPMY_StartTimer(void *p)
{
    MPMY_Timer *t = p;

    if (readrtc(&(t->start)))
	Warning("failed readrtc\n");
    return MPMY_SUCCESS;
}

int MPMY_StopTimer(void *p)
{
    MPMY_Timer *t = p;
    struct timestruc_t T__now, T__diff;

    if (readrtc(&T__now))
	Warning("failed readrtc\n");

    /* If the RTC value did not monotonically increase because it was
       corrected backwards in between samples (e.g. xntp), ignore it
       altogether instead of bombing from assertion */
    if ((T__now.tv_sec > t->start.tv_sec)
	|| ((T__now.tv_sec == t->start.tv_sec)
	    && (T__now.tv_nsec >= t->start.tv_nsec))) {
	ntimersub(T__now, t->start, T__diff);
	ntimeradd(T__diff, t->accum, t->accum);
    } else
	Msgf(("MPMY_StopTimer:Timer inconsistent\n"));
    return MPMY_SUCCESS;
}

