# SDF Control file testing a 1M centralized data set
char datafile[]="";
int nobj=1000000;
int seed=3141590;
int cencon=1;
int do_NlgN = 0;
float epsilon=0.02;
float frac_tol=0.025;
float errtol=0.01;
float dt=0.001;
int nsteps=40;
int cosmology=0;
int log_time=0;
int comov_eps=0;
int save_first=0;
char outfile[]="testfile";
int output_freq=5;
int timer_freq = 1;
int decomp_freq = 5;
char Msg_turn_on[] = "pqsort.c,decomp.c,image.c";
#char Msg_turn_on[] = "nomsgs";
int Msg_memfile = 65536;
# SDF-EOH

