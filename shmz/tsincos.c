#include <math.h>
#include <stdio.h>
#include "fastflpt.h"

void
main(int argc, char **argv)
{
    double s, c, x;
    
    x = 1.0;

    sincos(x, &s, &c);
    printf("%f %f\n", s, c);

    s = sin(x);
    c = cos(x);
    printf("%f %f\n", s, c);


    exit(0);
}
