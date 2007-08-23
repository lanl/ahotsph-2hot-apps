# SDF
char datafile[] = "cdm_ns.0650";
char SPHdatafile[] = "cdm_ns_sph.0650";
char outfile[] = "cdm_ns3";
int do_restart = 0;
int timeout = 1800;
float CWfac = 1.2;
int do_DL = 1;
int do_Arel = 1;
int do_periodic = 0;
float errtol = .5;
float frac_tol = .01;
float dt = 0.002;
int comov_eps = 1;
float comov_eps_epoch = 10.0;
float epsilon = 20;
int nsteps = 6;
int output_freq = 5;
int timer_freq = 5;
int cosmology = 1;
int save_first = 0;
#
int do_sph = 1;
float gamma = 1.66666666;
float visc_alpha = 0.0;
float visc_beta = 0.0;
float visc_epsilon = 0.0;
float heat_f1 = 0.0;
float courant_number = 0.25;
int exact_rho = 1;
int nbrcut_max = 150;
int nbrcut_min = 30;
float nbrcut_fac = 0.05;
int adaptive_dt = 1;
float dt_max = 0.02;
int tlow_cut = 200;
int dt_short = 0;
int dt_long = 10;
#
char Msg_turn_on[] = "main_nv.c,msgdone";
#int Msg_memfile = 65536;
# SDF-EOH

