/*
 * © 2026. Triad National Security, LLC. All rights reserved.
 * This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration. All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration. The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare. derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.
 */

#define NOTIMERS /* Timers are a major performance hit on the delta */
#include "Msgs.h"
#include "fastflpt.h"
#include "physics.h"
#include "stk.h"
#include "tensop.h"
#include "timers.h"
#include "vop.h"

#if defined(__T3D__) || defined(_IBMR2)
#define USE_CHEB_RSQRT
#endif

void do_grav(const float *p,
             const float *end,
             const float *pos0,
             float *mass0,
             float *acc0,
             float *phi0,
             const float *eps2p,
             int *ncut) {
    float dr2;
    Vxd(float r);
    float phii, mor3, mass;
    VxdV(const float ppos, = pos0);
    Vxd(float a);
    float phi = *phi0;
    float total_mass = *mass0;
    const float eps2 = *eps2p;
#ifdef ERR_INFO
    static int select;
    float bmax;
    float rcrit;
    if (++select % 1001 == 0) {
        Msgf(("# Particle\n"));
    }
#endif

    VxV(a, = acc0);

    while (p < end) {
        mass = *p++;
        r0 = *p++;
        r1 = *p++;
        r2 = *p++;
#ifdef ERR_INFO
        bmax = *p++;
        rcrit = *p++;
#endif

        VxVx(r, -= ppos); /* 3 flops */

        dr2 = Dotx(r, r); /* 5 flops */
        dr2 += eps2;

#ifdef ERR_INFO
        if (select % 1001 == 0) {
            Msgf(("%g %g %g %g\n", sqrt(dr2), mass, bmax, rcrit));
        }
#endif

        phii = recipsqrtf(dr2); /* 8 flops */

        mor3 = phii * phii; /* 5 flops */
        phii *= mass;
        total_mass += mass;
        mor3 *= phii;
        phi -= phii;

        VxVx(a, += mor3 * r); /* 6 flops */
    }
    VVx(acc0, = a);
    *mass0 = total_mass;
    *phi0 = phi;
}

void do_grav_u2(const float *p,
                const float *end,
                const float *pos0,
                float *mass0,
                float *acc0,
                float *phi0,
                const float *eps2p,
                int *ncut) {
    float dr2a, dr2b;
    Vxd(float ra);
    Vxd(float rb);
    float phiia, phiib;
    float mor3a, mor3b;
    float massa, massb;
    VxdV(const float ppos, = pos0);
    Vxd(float a);
    float phi = *phi0;
    float total_mass = *mass0;
    const float eps2 = *eps2p;

    VxV(a, = acc0);

    while (p < end) {
        massa = *p++;
        ra0 = *p++;
        ra1 = *p++;
        ra2 = *p++;

        massb = *p++;
        rb0 = *p++;
        rb1 = *p++;
        rb2 = *p++;

        VxVx(ra, -= ppos);
        VxVx(rb, -= ppos);

        dr2a = Dotx(ra, ra);
        dr2b = Dotx(rb, rb);

        dr2a += eps2;
        dr2b += eps2;

        phiia = recipsqrtf(dr2a);
        phiib = recipsqrtf(dr2b);

        mor3a = phiia * phiia;
        mor3b = phiib * phiib;
        phiia *= massa;
        phiib *= massb;
        total_mass += massa + massb;
        mor3a *= phiia;
        mor3b *= phiib;
        phi -= phiia;
        phi -= phiib;

        VxVx(a, += mor3a * ra);
        VxVx(a, += mor3b * rb);
    }
    VVx(acc0, = a);
    *mass0 = total_mass;
    *phi0 = phi;
}
