# SDF
char datafile[] = "cdm_m";
char outfile[] = "cdm_test";
int timeout = 1800;
float CWfac = 1.2;
int do_DL = 1;
int do_Arel = 1;
int do_periodic = 1;
float errtol = .5;
float frac_tol = .01;
float dt = 0.002;
int comov_eps = 1;
float comov_eps_epoch = 10.0;
float epsilon = 20;
int nsteps = 650;
int output_freq = 10;
int timer_freq = 5;
int cosmology = 1;
int save_first = 1;
int setpvel = 1;
#
int do_sph = 1;
float new_h = 10.0;
float new_u = 2.0;
float gamma = 1.66666666;
float visc_alpha = 1.0;
float visc_beta = 2.0;
float visc_epsilon = 1e-4;
float heat_f1 = 1.0;
float courant_number = 0.25;
int exact_rho = 0;
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

