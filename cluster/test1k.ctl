# SDF Control file testing a 1M data set
char datafile[]="test";
int nobj=1000;
int seed=3141590;
int cencon=0;
int do_sph = 1;
int exact_rho = 0;
float new_u = 0.25;
int nbrcut_max = 150;
int nbrcut_min = 30;
float nbrcut_fac = 0.1;
int do_periodic = 0;
int do_Arel = 1;
int do_DL = 1;
float epsilon=0.02;
float errtol=0.05;
float frac_tol=0.05;
float dt=0.01;
int nsteps=300;
int cosmology=0;
int log_time=0;
int comov_eps=0;
int save_first=1;
char outfile[]="testfile";
int output_freq=25;
int timer_freq = 10;
char Msg_turn_on[] = "main.c,tree.c,walk.c,cofm.c,sph.c";
#int Msg_memfile = 65536;
# SDF-EOH

