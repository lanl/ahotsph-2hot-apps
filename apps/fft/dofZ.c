#include <math.h>
#include <stdlib.h>

static float Omega0;
static float Lambda_prime;
static float H0;

float qromo(float (*func)(float), float a, float b,
	float (*choose)(float (*)(float), float, float, int));
float midpnt(float (*func)(float), float a, float b, int n);

float
adot(float a)
{
    return (H0/a)
      *sqrt(Omega0*a*(1.0-a) + a*a * (1.0 + Lambda_prime * (a*a - 1.0)));
}

float
addot(float a)
{
    return a * H0 * H0 * (Lambda_prime-0.5*Omega0/(a*a*a));
}

float
integrand(float a)
{
    float x;

    x = adot(a);
    return 1.0/(x*x*x);
    
}

float
t_integrand(float a)
{
    float x;

    x = adot(a);
    return 1.0/x;
    
}

float 
growthfac_from_Z(float omega0, float h0, float lambda_prime, float z)
{
    float a = 1.0/(1.0+z);
    Omega0 = omega0;
    H0 = h0;
    Lambda_prime = lambda_prime;
    return 2.5*H0*H0*adot(a)*qromo(integrand, 0.0, a, midpnt)/a;
}

float 
velfac_from_Z(float omega0, float h0, float lambda_prime, float z)
{
    float d, a_dot;
    float a = 1.0/(1.0+z);
    Omega0 = omega0;
    H0 = h0;
    Lambda_prime = lambda_prime;
    d = qromo(integrand, 0.0, a, midpnt);
    a_dot = adot(a);
    return addot(a)
      *a/(a_dot*a_dot) - 1.0 + a/(a_dot*a_dot*a_dot*d);
}

float
t_from_Z(float omega0, float h0, float lambda_prime, float z)
{
    float d;
    float a = 1.0/(1.0+z);
    Omega0 = omega0;
    H0 = h0;
    Lambda_prime = lambda_prime;
    d = qromo(t_integrand, 0.0, a, midpnt);
    return (d);
}

float
hubble_from_Z(float omega0, float h0, float lambda_prime, float z)
{
    float a = 1.0/(1.0+z);
    Omega0 = omega0;
    H0 = h0;
    Lambda_prime = lambda_prime;
    return adot(a)/a;
}

#ifdef STANDALONE

main(int argc, char *argv[])
{
    float z, norm;
    float f = sqrt(2.0);
    float omega0 = atof(argv[1]);
    float h0 = atof(argv[2]);
    float lp = atof(argv[3]);

    norm = growthfac_from_Z(1.0, h0, 0.0, 1e4)
      /growthfac_from_Z(omega0, h0, lp, 1e4);

    for (z = 0; z <= 1000; z++) {
	printf("%g %g %g %g\n", z, growthfac_from_Z(omega0, h0, lp, z)*norm,
	       velfac_from_Z(omega0, h0, lp, z), t_from_Z(omega0, h0, lp, z));
    }
}

#endif
