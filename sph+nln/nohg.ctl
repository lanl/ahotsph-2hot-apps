# SDF Control file testing a 1M data set
char SPHdatafile[]="noh1k.0";
char outfile[]="noh1k";
int do_sph = 1;
int do_grav = 1;
int do_DL = 1;
int do_Arel = 1;
int do_periodic = 0;
float errtol = .01;
float frac_tol = .01;
float epsilon = 1e-3;
float new_u = 1e-4;
int exact_rho = 0;
int nbrcut_max = 150;
int nbrcut_min = 30;
float nbrcut_fac = 0.1;
int do_periodic = 0;
int adaptive_dt = 1;
float dt = 0.001;
float dt_max = 0.064;
float default_nterms = 10;
int nsteps=300;
int save_first=0;
int output_freq=25;
int timer_freq = 1;
int ntimer_detail = 1;
#char Msg_turn_on[] = "main.c,decomp.c,pqsort.c";
#int Msg_memfile = 65536;
# SDF-EOH

