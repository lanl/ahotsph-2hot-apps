int n_x = 60;
int n_y = 60;
int n_z = 60;
#char Msg_turn_on[] = "decomp.c,tree.c,sph.c";

# New Spline Kernel Coefficients from Chris
# if omitted, the default is the Monaghan kernel 
int kernel_ncoef1 = 5;
struct {double kernel_coef1;}[5] = {1.2798, 0.0, -3.25703558174757,
3.0425, -0.83492};
int kernel_ncoef2 = 5;
struct {double kernel_coef2;}[5] = {2.1164, -3.3596555, 1.797,
-0.32368, 0.000272};

# Define interpolation range
# if omitted the size of the simulation will be used
struct {double outrmin;}[3] = {-10.,-20., -30.};
struct {double outrmax;}[3] = {10.,20., 30.};
# SDF-EOF
