/*
 * Copyright 1991, 1992, 1993 Michael S. Warren and John K. Salmon.  
 *  All Rights Reserved.
 */

#ifndef _TimersDOTh
#define _TimersDOTh

/* This used to look more like a Counter_t, but it was just too much */
/* trouble to deal with the archtecture specific differences, so we */
/* hide all the ugliness behind another layer of indirection.  In */
/* fact, it would make more sense if, e.g., EnableTimer returned a void* */
/* but that would require too much changing of "application" code. */
typedef struct {
    int enabled;
    void *mpmy_tm;
    char *name;
    double min, max, mean;
} Timer_t;
    
typedef struct {
    int enabled;
    unsigned int counter;
    unsigned int max, min;
    double mean, sum;
    char *name;
} Counter_t;

#ifdef __cplusplus
extern "C"{
#endif /* __cplusplus */
void StartTimer(Timer_t *t);
void StopTimer(Timer_t *t);
void SumTimers(void);
void EnableWCTimer(Timer_t *t, char *name);
void EnableCPUTimer(Timer_t *t, char *name);
void DisableTimer(Timer_t *t);
void ClearTimer(Timer_t *t);
void ClearEnabledTimers(void);
void OutputTimers(int (*Printf_Like)(const char *, ...));
void OutputTimer(Timer_t *t, int (*Printf_Like)(const char *, ...));
void OutputIndividualTimers(int (*Printf_Like)(const char *, ...));
double ReadTimer(Timer_t *t);

void SumCounters(void);
void EnableCounter(Counter_t *t, char *name);
void DisableCounter(Counter_t *t);
void ClearCounter(Counter_t *c);
void ClearEnabledCounters(void);
void OutputCounters(int (*Printf_Like)(const char *, ...));
void OutputIndividualCounters(int (*Printf_Like)(const char *, ...));
void OutputOneCounter(Counter_t *c, int (*Printf_Like)(const char *, ...));
int ReadCounter(Counter_t *c);
#ifdef __cplusplus
}
#endif /* __cplusplus */

#define EnableTimer(t, name) EnableCPUTimer(t, name)

#ifdef NOTIMERS
#define StartTimer(x) /**/
#define StopTimer(x) /**/
#define StartWCTimer(x) /**/
#define StopWCTimer(x) /**/
#endif

#ifndef NOCOUNTERS
#define IncrCounter(c) ((c)->counter++)
#define AddCounter(c, add) ((c)->counter += (add))
#else
#define IncrCounter(c)
#define AddCounter(c, add)
#endif

#endif /* _TimersDOTh */
