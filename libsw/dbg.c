#include <nx.h>

#define FORCE(t) ( (t) | (1<<30))
#define DBG_REQUEST_TYPE FORCE(0x1a2b3c4)

#define WHAT_DBGBUF 1
#define WHAT_MEMORY 2
#define WHAT_STACK 3
#define WHAT_CODE 4

typedef struct {
    int what;
} dbg_request_t;

void PrintMemfile();
static dbg_request_t dbg_request;

void
jab_dbg_handler(int proc)
{
    csend(DBG_REQUEST_TYPE, (char *)&dbg_request, sizeof(dbg_request), proc,0);
{

void
set_dbg_handler(void){
    hrecv(DBG_REQUEST_TYPE, (char *)&dbg_request,
	  sizeof(dbg_request), dbg_handler);
}

static void
dbg_handler(long type, long count, long node, long pid)
{
    PrintMemfile();
    set_dbg_handler();
    return;
}
