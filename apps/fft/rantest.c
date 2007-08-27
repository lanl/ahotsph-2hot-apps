#include <stdio.h>
#include <math.h>
#include "Assert.h"
#include "error.h"
#include "mpmy.h"
#include "randoms.h"
#include "bigmalloc.h"
#include "Msgs.h"

#define MAXNPROC 16

static long *ranp_seeds;
static int ranp_ntaps;
static int ranp_n;
static int ranp_available;
static int ranp_active;
float ran1(long *idum);
float ranp(ran_state *rs);
void ranp_setup(long seed, long n, ran_state *rs);
float normal_rand2(void);

ran_state rs;

void
main()
{
    int i;
    
    ranp_setup(123, 100000, &rs);

    for (i = 0; i < 10000; i++)
      printf("%f %f\n", normal_rand(&rs), normal_rand2());

    exit(0);
}

float
normal_rand2(void)
{
    float krlz, phase;

    krlz = sqrt(-log(ranp(&rs)));
    phase = 2.0*M_PI*ranp(&rs);
    return(krlz*cos(phase));
}

void
ranp_setup(long seed, long n, ran_state *rs)
{
    int i;

    assert(MAXNPROC % MPMY_Nproc() == 0);
    assert(n % MPMY_Nproc() == 0);
    ranp_ntaps = MAXNPROC/MPMY_Nproc();
    assert(ranp_ntaps);
    ranp_seeds = Malloc(ranp_ntaps * sizeof(long));

    if (seed <= 0) Error("Bad seed (%ld)\n", seed);
    for (i = 0; i < ranp_ntaps; i++)
      ranp_seeds[i] = -seed - i - MPMY_Procnum()*ranp_ntaps;
    ranp_n = n/MAXNPROC;
    ranp_active = 0;
    ranp_available = ranp_n;
    
    rs->next_norml_ok = 0;
}

void
ranp_reset(ran_state *rs)
{
    ranp_available = ranp_n;
    ranp_active++;
    rs->next_norml_ok = 0;
    Msg_do("ranp_reset, active is now %d\n", ranp_active);
}

float
ranp(ran_state *rs)
{
    float ret;
    if (ranp_active >= ranp_ntaps || ranp_active < 0)
      Error("ranp_active too large or no call to ranp_setup\n");
    ret = ran1(ranp_seeds+ranp_active);
    if (--ranp_available == 0) {
	Msg_do("ranp_available reset\n");
	ranp_available = ranp_n;
	ranp_active++;
	rs->next_norml_ok = 0;
    }
    return ret;
}


float normal_rand(ran_state *st)
/*
This is the Polar method for normal distributions, as described on or near
page 104 of Knuth, Semi-numerical Algorithms.  To quote Knuth, "The polar
method is quite slow, but it has essentially perfect accuracy, and it is very
easy to write a program for the polar method..."  'nuf said.  Algorithm due
to Box, Muller and Marsaglia.
*/
{
    float	v1, v2;	/* uniformly distributed on [-1, 1) */
    float s;	/* radius of a point pulled from a uniform circle */
    float	foo;	/* A useful intermediate value. */
    
    if(st->next_norml_ok){
	st->next_norml_ok = 0;
	return st->next_norml;
    }
    
    do{
	v1 = 2.0F * ranp(st) - 1.0F;
	v2 = 2.0F * ranp(st) - 1.0F;
	s = v1*v1 + v2*v2;
    } while(s >= 1.0F);
    foo = sqrt( -2.0F * log(s)/s);
    st->next_norml_ok = 1;
    st->next_norml = v1*foo;
    return v2*foo;
}
