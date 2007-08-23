# SDF
char SPHdatafile[] = "mss2.0";
char outfile[] = "msms7k";
int do_sph = 1;
int do_grav = 0;
int exact_rho = 0;
float epsilon = 0.01;
float errtol=0.05;
float frac_tol=0.05;
float dt = .03;
int do_Arel = 1;
int nbrcut_max = 150;
int nbrcut_min = 30;
float nbrcut_fac = 0.1;
float max_h = 0.2;
float min_h = .01;
int adaptive_dt = 1;
int independent_dt = 0;
int nsteps = 2000;
int output_freq = 50;
int timer_freq = 10;
int save_first = 0;
float visc_alpha = 1.0;
float visc_beta = 2.5;
#char Msg_turn_on[]="pqsort.c,main_sph.c";
# SDF-EOH

