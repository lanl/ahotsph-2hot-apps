# SDF
char datafile[] = "data/panel_d4.0";
char outfile[] = "data/pd4e"; /* use empty string, or don't define for none */
char restartfile[] = "data/panel_d4.rst"; 
int nsteps = 10; 
int output_freq = 2;
int restart_freq = 0;
int timer_freq = 1;
int save_first = 0;
#float U_infty_x =0.;
#float U_infty_y =0.;
#float U_infty_z =1.;

float Dipole_x = 1.0;
float Dipole_y = 0.;
float Dipole_z = 0.;
float Dipole_pos_x = 1.5;
float Dipole_pos_y = 0.;
float Dipole_pos_z = 0.;
	
float errtol = 1.e-3;
float Jacobi_relax = 0.75;
int check_errs=0;
char Msg_turn_on[]="nomsgs";
# SDF-EOH
