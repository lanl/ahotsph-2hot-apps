#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <netinet/in.h>
#include <netdb.h>
#include <sys/socket.h>
#include <fcntl.h>
#include <signal.h>
#include <errno.h>
#include "error.h"
#include "singlio.h"
#include "Msgs.h"
#include "protos.h"
#include "mpmy.h"
#include "byteswap.h"

#if defined(__SUN4__) && !defined(_SUNOS4_PROTOTYPES_)
/* These should really be prototyped somewhere else, 
   but this will do for now */
int close(int);
int getpid(void);
int sigvec(int sig, struct sigvec *vec, struct sigvec *ovec);
int alarm(int);
int ioctl(int fd, int cmd, void *p);
int gethostname(char *name, int namelen);
int inet_addr(char *cp);
char *inet_ntoa(struct in_addr in);
int socket(int domain, int type, int protocol);
int bind(int s, struct sockaddr *name, int namelen);
int recv(int s, void *buf, int len, int flags);
int recvfrom(int s, void *buf, int len, int flags, 
	     struct sockaddr *from, int *fromlen);
int sendto(int s, const void *msg, int len, 
	   int flags, struct sockaddr *to, int tolen);
#endif

void PrintMemfile(void);
#ifdef __PARAGON__
static void jab_dbg_handler(int proc);
static void set_dbg_handler(void);
static void dbg_handler(long type, long count, long node, long pid);
#endif

static void sock_init(char *hostname, int *port, 
		      struct sockaddr_in *acc, int bind_flag);

static void setup_handler(void);
static void io_ready(int);

static int sock;		/* file descriptor for my UDP socket */

void
sigio_setup(void)
{
    int port = 4000;
    struct sockaddr_in my_addr;

#ifdef __PARAGON__
    set_dbg_handler();
    if (MPMY_Procnum())
      return;
#endif

    setup_handler();
    sock_init(NULL, &port, &my_addr, 1); /* get my sockaddr */
    
    if (fcntl(sock, F_SETOWN, getpid()) < 0) 
      Error("F_SETOWN error\n");

#ifdef FASYNC
    if (fcntl(sock, F_SETFL, FASYNC) < 0)
      Error("F_SETFL FASYNC error\n");
#endif

}


static void
setup_handler(void)
{
#if !defined(__iX86__) && !defined(_AIX) && !defined(__SUN5__)
    struct sigvec vec;

    vec.sv_handler = io_ready;
    vec.sv_flags = SV_INTERRUPT;
    vec.sv_mask = 0;
    sigvec(SIGIO, &vec, NULL);
#else
    signal(SIGIO, io_ready);
#endif
}

static void
io_ready(int sig)
{
    int node;


#ifndef __PARAGON__
    PrintMemfile();
#if defined(__iX86__) || defined(_AIX) || defined(__SUN5__)
    signal(SIGIO, io_ready);
#endif
    return;
#endif

    if (recv(sock, &node, sizeof(int), 0) != sizeof(int))
      Error("recv failed, errno=%d", errno);

    /* This is not very pretty, we probably should transmit a byteorder */

    if (node < 0 || node >= MPMY_Nproc())
      Byteswap(sizeof(int), 1, &node, &node);
    if (node < 0 || node >= MPMY_Nproc()) {
	fprintf(stderr, 
		"Node to jab (%d) is out of range, and swapping didn't help\n",
		node);
    } else {
#ifdef __PARAGON__
    jab_dbg_handler(node);
#endif
    }
}

/* This is virtually identical to the lsv code */

static void
sock_init(char *hostname, int *port, struct sockaddr_in *acc, int bind_flag)
{
    struct hostent *hp;
    char host_name[256];
    unsigned long inaddr;
    int tries = 0;
    
    if (hostname == NULL) {
	if( (hostname = getenv("LSV_HOSTNAME")) == NULL ){
	    if (gethostname(host_name, 256))
		Error("sock_create: gethostname failed\n");
	    hostname = host_name;
	}
    }
    memset(acc, 0, sizeof(struct sockaddr_in));
    acc->sin_family = htons(AF_INET);

    if ((inaddr = inet_addr(hostname)) != -1) /* it is numeric */
      acc->sin_addr.s_addr = inaddr;
    else if ((hp = gethostbyname(hostname)) != (struct hostent *)0)
      memcpy(&(acc->sin_addr), hp->h_addr, hp->h_length);
    else
      Error("gethostbyname failed\n");
    /* printf("addr for %s is %s\n", 
	hostname, inet_ntoa(acc->sin_addr));
    printf("htons(port) is %d\n", htons(*port)); */

    if (bind_flag) {
	sock = socket(AF_INET, SOCK_DGRAM, 0);
    }
  try_again:
    acc->sin_port = htons(*port);
    if (bind_flag) {
	int ret;
	ret = bind(sock,(struct sockaddr *)acc,sizeof(struct sockaddr_in));
	if (ret < 0 ) { 
	    if (tries < 100)  {
		/* printf("bind returns %d\n", ret); */
		(*port)++;
		tries++;
		goto try_again;
	    } else {
		Error("Can't bind socket. Tried %d, up to port %d\n", 
		      tries, *port);
		exit(1);
	    }
	}
	singlPrintf("sigio_dump at %s port %d\n", hostname, *port);
	errno = 0;		/* clear errors */
    }
}

#ifdef __PARAGON__
#include <nx.h>

#define FORCE(t) ( (t) | (1<<30))
#define DBG_REQUEST_TYPE FORCE(0x1a2b3c4)

typedef struct {
    int what;
} dbg_request_t;

static dbg_request_t dbg_request;

static void
jab_dbg_handler(int proc)
{
    fprintf(stderr, "jab node %d\n", proc);
    csend(DBG_REQUEST_TYPE, (char *)&dbg_request, sizeof(dbg_request), proc,0);
}

static void
set_dbg_handler(void)
{
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
#endif
