# SDF Control file testing a 4M centralized data set
char datafile[]="";
int nobj=4000000;
int seed=3141590;
int cencon=1;
int do_NlgN = 0;
float epsilon=0.02;
float frac_tol=0.025;
float errtol=0.01;
float dt=0.001;
int nsteps=100;
int cosmology=0;
int log_time=0;
int comov_eps=0;
int save_first=0;
char outfile[]="testfile";
int output_freq=100;
int timer_freq = 1;
float sort_tol = .001;
#char Msg_turn_on[] = "nomsgs";
char Msg_turn_on[] = "nlcomm.c";
int Msg_memfile = 65536;
# SDF-EOH

