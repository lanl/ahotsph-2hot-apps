#include "mpmy_time.h"

/* The CMMD timers have a really nice interface.
   We can use integers (of our own chosing) between 0-64 to identify timers. 
   So we play fast-and-loose and stick the type in the upper bits, and
   then assume the whole thing will fit in a void*. 
   I think the code here is ok.  It's completely lost if the caller
   frees the same timer twice, or frees an unused timer.  Caveat emptor.*/

static int freeid[64];
static int initialized = 0;
static int ntimer = 0;

void *MPMY_CreateTimer(int type){
    int id;
    if( !initialized ){
	int i;
	initialized = 1;
	for(i=0; i<64; i++){
	    freeid[i] = i;
	}
    }

    if( ntimer == 64 ){
	return (void *)-1;
    }
    id = freeid[ntimer++];
    CMMD_node_timer_clear(id);
    return (void *) ( id | (type<<6));
}

int MPMY_DestroyTimer(void *p){
    int id = ((int)p)&63;
    if( ntimer == 0 )
	return MPMY_FAILED;

    freeid[--ntimer] = id;
    return MPMY_SUCCESS;
}

int MPMY_StartTimer(void *p){
    int id = (int)p & 63;
    
    CMMD_node_timer_start(id);
    return MPMY_SUCCESS;
}

int MPMY_StopTimer(void *p){
    int id = (int)p & 63;
    
    CMMD_node_timer_stop(id);
    return MPMY_SUCCESS;
}

int MPMY_ClearTimer(void *p){
    int id = (int)p & 63;
    
    CMMD_node_timer_clear(id);
    return MPMY_SUCCESS;
}

double MPMY_ReadTimer(void *p){
    int id = (int)p & 63;
    int type = ((int)p) >> 6;
    
    switch(type){
    case MPMY_WC_TIME:
	return CMMD_node_timer_elapsed(id);
    case MPMY_CPU_TIME:
	return CMMD_node_timer_busy(id);
    }
    return -1.0;
}
