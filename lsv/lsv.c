		/* The generic 'IH' program */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include "lsv.h"
#include "Msgs.h"

extern time_t time(time_t *tloc);
extern char *cuserid(char *);

static void parse(int argc, char **argv);
static void usage(void);
void log_start(void);
void log_end(void);
void run(void);

#ifndef LOGFILE
#define LOGFILE "/tmp/lsv_log"
#endif
static time_t date;
static FILE *logfile;

int gargc, gdoc = 0;
char *machine_list = NULL;
int verbose = 0;
int inet_numeric = 0;
int messages = 0;
int suspend_proc = -2;
int gnproc = 1;
char **gargv, *progelt;
char *start_prog = "lsvstart";

void
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
	case 'M':
	    messages = 1;
	    break;
	case 'S':
	    suspend_proc = atoi(*++argv);
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
    
    if(messages){
	extern int vfprintf();
	extern int fflush();
	Msg_addfile(stdout, (Msgvfprintf_t)vfprintf, (Msgfflush_t)fflush);
	Msg_on("mpmy_lsv.c");
	Msg_on("lsv.c");
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
    fprintf(stderr, "\t -M         \t\t Turn on msgs\n");
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

void run(void)
{
    int host_pid;
    int port;
    char *myname;
    /* Very simple for now.  */

    Sinit_host1(&port, &myname);
    if( (host_pid=fork()) == 0 ){
	Sinit_host(gnproc);
	exit(0);
    }else{
	putenv_int("LSV_NPROC", gnproc);
	putenv_int("LSV_HOSTPORT", port);
	putenv_str("LSV_HOST", myname);
	/* gargv[0] is progelt, but that is lost by the exec, so we */
	/* pass it in the environment */
	putenv_str("LSV_PROGELT", progelt);
	putenv_int("LSV_PID", host_pid);
	if( machine_list )
	    putenv_str("LSV_REMOTE", machine_list);
	if( suspend_proc != -2 )
	    putenv_int("LSV_SUSPEND", suspend_proc);
	if( verbose )
	    putenv_int("LSV_VERBOSE", 1);
	execvp(start_prog, gargv);
    }
}

