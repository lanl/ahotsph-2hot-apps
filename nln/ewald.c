/*
 * © 2026. Triad National Security, LLC. All rights reserved.
 * This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration. All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration. The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare. derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.
 */

#include <math.h>
#define NDIM 3
#include "error.h"
#include "vop.h"

#define N 4

void ewald(double *x, double L, double *f, double *phi) {
    int i[NDIM];
    double xx[NDIM];
    double ax2;
    double axx;
    double hh;
    double t, t1, t2;
    double alpha;
    double a1, a2, a3, a4, a5;

    alpha = 2.0 / L;
    VS(f, = 0);
    *phi = M_PI / (alpha * alpha * L * L * L);
    if (x[0] == 0.0 && x[1] == 0.0 && x[2] == 0.0)
        return;
    a1 = 2.0 * alpha / sqrt(M_PI);
    a2 = -M_PI * M_PI / (alpha * alpha * L * L);
    a3 = 2.0 * M_PI / L;
    a4 = 2.0 / (L * L);
    a5 = 1.0 / (M_PI * L);


    for (i[0] = -N; i[0] <= N; i[0]++) {
        for (i[1] = -N; i[1] <= N; i[1]++) {
            for (i[2] = -N; i[2] <= N; i[2]++) {
                VVV(xx, = x, -L * i);
                ax2 = Dot(xx, xx);
                axx = sqrt(ax2);
                if (axx > 4.0 * L)
                    continue;
                t1 = erfc(alpha * axx);
                t2 = a1 * exp(-alpha * alpha * ax2);
                t = (t1 / axx + t2) / ax2;
                VV(f, += t * xx);
                *phi -= t1 / axx;
            }
        }
    }
    for (i[0] = -N; i[0] <= N; i[0]++) {
        for (i[1] = -N; i[1] <= N; i[1]++) {
            for (i[2] = -N; i[2] <= N; i[2]++) {
                hh = Dot(i, i);
                if (hh > 16 || !hh)
                    continue;
                t = a4 * exp(a2 * hh) * sin(a3 * Dot(i, x)) / hh;
                VV(f, += t * i);
                t = a5 * exp(a2 * hh) * cos(a3 * Dot(i, x)) / hh;
                *phi -= t;
            }
        }
    }
}
