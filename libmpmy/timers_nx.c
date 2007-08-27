#define HAVE_MPMY_TIMERS
/* Assume that nx.h/cube.h/mesh.h is already included. */
#include "mpmy_time.h"
#include "chn.h"
#include "bigmalloc.h"

#ifndef __PARAGON__
extern void hwclock(esize_t *);
#else
/* The value in nx.h is apparently incorrect */
#undef HWHZ
#define HWHZ 50000000
#endif

/* 
#if !( defined(__SUNMOS__) && defined(__PARAGON__) )
#define USE_E_ROUTINES
#endif
*/

static Chn timer_chn;
static int initialized;

typedef struct {
    int type;
    esize_t start, accum;
} MPMY_Timer;

void *MPMY_CreateTimer(int type){
    MPMY_Timer *ret;

    if( initialized == 0 ){
	ChnInit(&timer_chn, sizeof(MPMY_Timer), 40, Realloc_f);
	initialized = 1;
    }

    ret = ChnAlloc(&timer_chn);
    ret->type = type;		/* NOT USED! */
    MPMY_ClearTimer(ret);
    return (void *)ret;
}

int MPMY_DestroyTimer(void *p){
    ChnFree(&timer_chn, p);
    return MPMY_SUCCESS;
}

int MPMY_StartTimer(void *p)
{
    MPMY_Timer *t = p;
    hwclock(&(t->start));
    return MPMY_SUCCESS;
}

#if defined(USE_E_ROUTINES)
/* This version doesn't assume anything about the internals organization */
/* of an esize.  It uses the 'official' e-routines interface*/
#define LGBIG 30
#define BIG (1<<LGBIG)

int MPMY_ClearTimer(void *p)
{
    MPMY_Timer *t = p;
    static int ezero_init;
    static esize_t ezero;
    /* Who knows how slow stoe is... */
    if( !ezero_init ){
	ezero = stoe("0");
	ezero_init = 1;
    }
    t->accum = ezero;
    return MPMY_SUCCESS;
}

double MPMY_ReadTimer(void *p)
{
    MPMY_Timer *t = p;
    double div, rem;
    /* Sometimes hwclock returns garbage, which causes ediv to */
    /* get a range error.  We don't want this minor inconvenience to */
    /* blow away a long run. */
    div = (double)_ediv(t->accum, BIG);
    rem = (double)_emod(t->accum, BIG);
    return (BIG*div + rem)/(double)HWHZ;
}

int MPMY_StopTimer(void *p)
{
    MPMY_Timer *t = p;
    esize_t T__now;

    hwclock(&T__now);
    T__now = esub(T__now, t->start);
    t->accum = eadd(t->accum, T__now);
    return MPMY_SUCCESS;
}
#else
/* Here, we assume we know something about the internal organization */
/* of an esize */
#define TWO_TO_32 (65536.*65536.)

int MPMY_ClearTimer(void *p)
{
    MPMY_Timer *t = p;
    t->accum.shigh = 0;
    t->accum.slow = 0;
    return MPMY_SUCCESS;
}

double MPMY_ReadTimer(void *p)
{
    MPMY_Timer *t = p;
    return (((unsigned long)t->accum.shigh)*TWO_TO_32 
	    + (unsigned long)t->accum.slow)*(1./HWHZ);
}

#define MAX (0xffffffff)
int MPMY_StopTimer(void *p)
{
    MPMY_Timer *t = p;
    esize_t T__now;
    typedef unsigned long int u32;
    u32 high, low, tmp;

    hwclock(&T__now);
    /* (high,low) = esub(T__now, t->start) */
    if( (u32)T__now.slow > (u32)t->start.slow ){
	low = (u32)T__now.slow - (u32)t->start.slow;
	high = T__now.shigh - t->start.shigh;
    }else{
	low = (MAX - (u32)t->start.slow) + (u32)T__now.slow + 1;
	high = (T__now.shigh - 1) - t->start.shigh;
    }
    /* t->accum = eadd(t->accum, (high,low)); */
    tmp = (u32)t->accum.slow;
    t->accum.slow += low;	
    t->accum.shigh += high;
    /* detect overflow of low bits here... */
    if( (u32)t->accum.slow < low || (u32)t->accum.slow < tmp )
	t->accum.shigh++;
    return MPMY_SUCCESS;
}
#endif



