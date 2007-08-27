		/* The generic 'IH' program */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <errno.h>
#include "Assert.h"
#include "error.h"
#include "swampi.h"
#include "mpmy_abnormal.h"
#include "Msgs.h"

extern time_t time(time_t *tloc);
#ifndef __POWERPC__
extern char *cuserid(char *);
#else
char *cuserid(char *n) {
    static char buf[16];

    if (n == NULL) {
        strcpy(buf, "Who am I?");
        return buf;
    } else {
        return NULL;
    }
}
#endif

static void parse(int argc, char **argv);
static void usage(void);
static char *file_to_string(const char *name);
void log_start(void);
void log_end(void);
void run(void);

#ifndef LOGFILE
#define LOGFILE "/ss1-raid/gaber/swampi_log"
#endif
static time_t date;
static FILE *logfile;

int gargc, gdoc = 0;
char *machine_list = NULL;
char *message_turnon = NULL;
int timeout = 0;
int verbose = 0;
int inet_numeric = 0;
int suspend_proc = -2;
int gnproc = 1;
char **gargv, *progelt;
char *start_prog = "swampistart";

int
main(int argc, char **argv, char **envp)
{
    parse(argc, argv);
    log_start();

    run();

    log_end();
    exit(0);		/* not reached */
}

static void
parse(int argc, char **argv)
/*
	parse the command line.  One arg at a time, some, like -d
	cause the next (few) args to be read.  The first one without a
	leading '-' is the elt program, and all following ones are passed
	to it in gargc, gargv.  Note that the elt program name must not
	begin with a '-'.
*/
{
    char *progname;

    progname = argv[0];
    if ( argc < 2 ) {
	fprintf(stderr,"Usage: %s [options] progelt [eltargs]\n", 
		argv[0]);
	fprintf(stderr,"%s: not enough arguments\n", argv[0]);
	exit(1);
    }
    
    while ( --argc && (*++argv)[0] == '-' ) {
	switch( (*argv)[1] ) {
	case 'd':
	    gdoc = atoi(*++argv); /* Dimension of cube */
	    gnproc= 1<<gdoc;
	    argc--;
	    break;
	case 'n':
	    gnproc = atoi(*++argv);
	    argc--;
	    break;
	case 'm':
	    machine_list = *++argv;
	    argc--;
	    break;
	case 'S':
	    suspend_proc = atoi(*++argv);
	    argc--;
	    break;
        case 't':
	    timeout = atoi(*++argv);
	    argc--;
	    break;
	case 'M':
	    message_turnon = *++argv;
	    argc--;
	    break;
	case 'i':
	    inet_numeric = 1;
	    break;
	case 'v':
	    verbose = 1;
	    break;
	case 'x':
	    start_prog = *++argv;
	    argc--;
	    break;
	default:
	    fprintf(stderr, "Unknown argv: %s\n", *argv);
	    usage();
	}
    }
    
    if(message_turnon){
	extern int _MPMY_procnum_;
	MsgdirInit("msgs/msg.host");
	Msg_turnon(message_turnon);
	_MPMY_procnum_ = -1;
	_MPMY_setup_absigs();
	MPMY_OnAbnormal(MPMY_SystemAbort);
	MPMY_OnAbnormal(MPMY_Abannounce);
    }
    if( machine_list && machine_list[0] == '<' ){
      /* Perl does this really well... */
      machine_list = file_to_string(machine_list+1);
    }
      
    gargc = argc;
    gargv = argv;
    progelt = argv[0];
    if ( gargc < 1 ) {
	fprintf(stderr,"Usage: %s [options] progelt [eltargs]\n",
		progname);
	fprintf(stderr,"%s: missing progelt argv\n", progname);
	exit(1);
    }
}

static void
usage(void)
{
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "\t -d dim     \t\t set number of elt procs to 2**dim\n");
    fprintf(stderr, "\t -n nproc   \t\t set number of elt procs to nproc\n");
    fprintf(stderr, "\t -v         \t\t be verbose\n");
    fprintf(stderr, "\t -i         \t\t use numeric inet address for host\n");
    fprintf(stderr, "\t -t time    \t\t Set timeout alarm in seconds\n");
    fprintf(stderr, "\t -M file.c  \t\t Turn on msgs for file.c on nodes\n");
    fprintf(stderr, "\t -S proc    \t\t Suspend proc (-1 for all)\n");
    fprintf(stderr, "\t -m machine \t\t run node program on machine\n");
    fprintf(stderr, "\t -m \"machine list\" \t run node program on several machines\n");
    fprintf(stderr, "\t -x startscript \t\t use alternative startscript\n");
}

void
log_start(void) 
{
    logfile = fopen(LOGFILE, "a");
    date = time(0);
    if (logfile) {
	if( machine_list )
	    fprintf(logfile, "%10s %10s %2d %10s  %s", cuserid(NULL),
		    machine_list, gdoc, progelt, ctime(&date));
	else
	    fprintf(logfile, "%10s (local) %2d %10s  %s", cuserid(NULL),
		    gdoc, progelt, ctime(&date));
	    
	fflush(logfile);
	fclose(logfile);
    }
}

void
log_end(void)
{
    logfile = fopen(LOGFILE, "a");
    if (logfile) {
	fprintf(logfile, "%10s %5ld sec\n", cuserid(NULL), (long)(time(0)-date));
	fclose(logfile);
    }
}

void putenv_int(const char *name, int value){
    char *p = malloc(strlen(name) + 14);
    sprintf(p, "%s=%d", name, value);
    putenv(p);
}

void putenv_str(const char *name, const char *value){
    char *p = malloc(strlen(name) + strlen(value) + 2);
    sprintf(p, "%s=%s", name, value);
    putenv(p);
}

#define GROWBY 1024
static char *
file_to_string(const char *name){
  FILE *fp = fopen(name, "r");
  char *pstart, *pend, *p;
  int done, c, sz;
  
  if( fp == NULL )
    return NULL;
  sz = GROWBY;
  p = pstart = malloc(sz);
  pend = pstart + sz;
  done = 0;
  while( !done ){
    c=getc(fp);
    /* A convenient way to append the trailing NUL */
    if( c == EOF ){
      done = 1;
      c = '\0';
    }
    assert(p<=pend);
    if( p==pend ){
      int deltap = p - pstart;
      assert( deltap == sz );
      sz += GROWBY;
      pstart = realloc(pstart, sz);
      if( pstart == NULL ){
	fclose(fp);
	return NULL;
      }
      p = pstart + deltap;
      pend = pstart + sz;
    }
    *p++ = c;
  }
  return pstart;
}

void run(void)
{
    int host_pid;
    int port;
    char *myname;
    /* Very simple for now.  */

    _MPI_init_host1(&port, &myname, gnproc);
    if( (host_pid=fork()) == 0 ){
	if (timeout) MPMY_TimeoutSet(timeout);
	_MPI_init_host(gnproc);
	exit(0);
    }else{
	if (timeout) MPMY_TimeoutSet(timeout);
	putenv_int("MPI_NPROC", gnproc);
	putenv_int("MPI_HOSTPORT", port);
	putenv_str("MPI_HOST", myname);
	/* gargv[0] is progelt, but that is lost by the exec, so we */
	/* pass it in the environment */
	putenv_str("MPI_PROGELT", progelt);
	putenv_int("MPI_PID", host_pid);
	if( machine_list )
	    putenv_str("MPI_REMOTE", machine_list);
	if( message_turnon )
	    putenv_str("MPI_MESSAGE_TURNON", message_turnon);
	if( timeout )
	    putenv_int("MPI_TIMEOUT", timeout);
	if( suspend_proc != -2 )
	    putenv_int("MPI_SUSPEND", suspend_proc);
	if( verbose )
	    putenv_int("MPI_VERBOSE", 1);
	execvp(start_prog, gargv);
	if( errno == ENOENT ){
	  Error("Can't exec %s, file does not exist (ENOENT)\n", start_prog);
	}else{
	  Error("Returned from execvp.  Very Bad News.  errno=%d\n", errno);
	}
    }
}

