# SDF Control file testing a 55k merger of two isothermals.
# see mergecmd for the details of how it was made.
# The smaller one is a 
# modified isothermal with M=0.1, R=0.1, G=1, eps=0.01, gamma=0.01 halo.
# Vc is 1, so rotaton period is ~0.6.  We integrate with a timestep
# of 0.01.  Will we still see a big spray?
# an errtol of .02 may be ok.  Who knows??
char datafile[]="/scratch/sd3c/johns/mrgr2/mrgr2.1050";
char outfile[]="/scratch/sd3c/johns/mrgr2/mrgr2cm5";
#char hdrfile[]="/data/merlin/pjq/pjq.sdfh";
int do_NlgN = 0;
float epsilon=0.01;
float errtol=0.02;
float dt=0.01;
int nsteps=50;
int cosmology=0;
int log_time=0;
int comov_eps=0;
int save_first=0;
int output_freq=50;
int timer_freq = 5;
char Msg_turn_on[] = "nomsgs";
#int Msg_memfile = 65536;
# SDF-EOH

