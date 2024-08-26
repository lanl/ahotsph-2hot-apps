#include <math.h>
#include <stdio.h>

#include "Assert.h"
#include "Malloc.h"
#include "Msgs.h"
#include "chn.h"
#include "error.h"
#include "fastflpt.h"
#include "order.h"
#include "physics.h"
#include "protos.h"
#include "tree.h"
#include "vop.h"

static const mac_s *mac;
static membuf_s *cellbuf;
static float min_sigma_m = 1.0f;
static float max_sigma_m = 7.2e+10; /* exp(25.0) */
                                    // clang-format off
static float sigma_m[] = {
    1.0,			/* 0, 10^10 Msol */
    1.0,
    1.0,
    1.0,
    1.0,
    1.0,
    1.0,
    1.0,
    1.0,			/* 8, 3 x 10^13 Msol, r = 7.8/h Mpc */
    1.5,
    2.0,
    3.0,
    6.0,
    10.88,
    15.18,
    21.34,
    31.58,
    46.73,
    74.63,
    120.3,
    199.9,
    353.1,
    611.2,
    1117,			/* 23, r = 1000/h Mpc */
    1983};			/* 24, 2.65 x 10^20 Msol */
// clang-format on
static void mpole_add_cube(double a, double rho, cofmdata *q);
static void mpole_add_mono(cofmdata *cmp, body *bp);
static void mpole_add(cofmdata *cmp, cofmdata *dp);

static double a[14]; /* coef of error poly */
static void rcrit_poly(int n, double r, double *value, double *deriv);
static double rtnewt(int n, void (*funcd)(int, double, double *, double *), double x1, double xacc);

void SetupCofm(mac_s *m, membuf_s *c) {
    m->rho0 = m->m0 / pow3(2.0 * m->r0);
    m->inv_tol = 1.0 / m->this_tol;
    m->inv_rel_tol = 1.0 / m->rel_tol;
    m->inv_rel_tol0 = 1.0 / m->rel_tol0;
    m->bmax0 = m->r0;
    mac = m;
    cellbuf = c;
}

/* Make Cell bounds correct if system is smaller than root cell */
static void CellCenterFix(float c[NDIM], float *size) {
    float rmin, rmax, s[NDIM];

    for (int i = 0; i < NDIM; i++) {
        rmin = Max(c[i], -mac->r0);
        rmax = Min(c[i] + *size, mac->r0);
        c[i] = 0.5f * (rmin + rmax);
        s[i] = rmax - rmin;
    }

    float smin = mac->r0 * 100.0f;
    float smax = 0.0f;

    for (int i = 0; i < NDIM; i++) {
        smin = Min(smin, s[i]);
        smax = Max(smax, s[i]);
    }
    if (smin / smax < 0.99f || smin / smax > 1.01f) {
        *size = -1.0f; /* flag non-cubic cells with negative size */
    } else if (smax < *size * 0.99f) {
        *size = smax;
    }
}

void CofmFromDaugh(hcellptr hptr, hcellptr daughters[]) {
    int i;
    cofmdata *dp;
    cofmdata *cmp;
    body *bp = NULL;
    float newbmax;
    float center[NDIM], cellsz;
    double dmass;
    double cofm[NDIM];
    double tmp[NDIM];
    double tmpsq;
    Vxd(float dx);

    assert(Sub_Flags(hptr));

    cmp = hptr->ptr;
    assert(cmp);
    memset(cmp, 0, sizeof(cofmdata));
    CELLCORNER(hptr->key, center, &cellsz);
    cmp->level = TreeLevel(hptr->key, NDIM);
    if (mac->expand_root > 0.0) {
        CellCenterFix(center, &cellsz);
        cmp->sz = cellsz;
    } else {
        VS(center, += 0.5f * cellsz);
        cmp->sz = cellsz;
    }

    /* Count daughters. */
    /* Point bptr at first body in cell */
    for (i = 0; i < (1 << NDIM); i++) {
        if (daughters[i] == NULL)
            continue;
        if (Sub_Flags(daughters[i]) == 0) {
            bp = daughters[i]->ptr;
            if (cmp->ndaughters == 0)
                cmp->bptr = bp;
            cmp->ndaughters++;
        } else {
            dp = daughters[i]->ptr;
            if (cmp->ndaughters == 0)
                cmp->bptr = dp->bptr;
            cmp->ndaughters += dp->ndaughters;
        }
    }
    if (mac->geometric_center) {
#ifndef DIPOLE
        Error("Can't do geometric center expansion without DIPOLE\n");
#endif
        /* Use geometric center */
        VV(cmp->center, = center);
    } else {
        /* Use cofm for center */
        dmass = 0.0;
        VS(cofm, = 0.0);
        /* relative to cell center for smaller truncation error */
        for (i = 0; i < (1 << NDIM); i++) {
            if (daughters[i] == NULL)
                continue;
            if (Sub_Flags(daughters[i]) == 0) {
                bp = daughters[i]->ptr;
                dmass += bp->mass;
                VVV(tmp, = center, -bp->pos);
                VV(cofm, += bp->mass * tmp);
            } else {
                dp = daughters[i]->ptr;
                dmass += dp->m;
                VVV(tmp, = center, -dp->center);
                VV(cofm, += dp->m * tmp);
            }
        }
        if (dmass != 0.) {
            VS(cofm, *= 1.0 / dmass);
            VVV(cmp->center, = center, -cofm);
        } else {
            Error("Zero mass in BranchFromDaughters!\n");
        }
    }
    /* Now loop again to pick up B2, etc.  */
    for (i = 0; i < (1 << NDIM); i++) {
        if (daughters[i] == NULL)
            continue;
        if (Sub_Flags(daughters[i]) == 0) {
            bp = daughters[i]->ptr;
            VVV(tmp, = cmp->center, -bp->pos);
            newbmax = 0.;
            tmpsq = Dot(tmp, tmp);
            if (tmpsq != 0.) {
                cmp->B2 += bp->mass * tmpsq;
                newbmax += sqrt(tmpsq);
                mpole_add_mono(cmp, bp);
            } else {
                cmp->m += bp->mass;
            }
        } else {
            dp = daughters[i]->ptr;
            VVV(tmp, = cmp->center, -dp->center);
            newbmax = dp->bmax;
            tmpsq = Dot(tmp, tmp);
            cmp->B2 += dp->B2;
            mpole_add(cmp, dp);
            if (tmpsq != 0.) {
                cmp->B2 += dp->m * tmpsq;
                newbmax += sqrt(tmpsq);
            }
        }
        if (newbmax > cmp->bmax)
            cmp->bmax = newbmax;
    }
    if (cmp->B2 == 0.0) {
        for (i = 0; i < (1 << NDIM); i++) {
            if (daughters[i] == NULL)
                continue;
            if (Sub_Flags(daughters[i]) == 0) {
                bp = daughters[i]->ptr;
                Msg_do("B2 body %s\n", PrintBodyContents(bp));
            } else {
                dp = daughters[i]->ptr;
                Msg_do("B2 cell\n");
            }
        }
        cmp->bmax = fabs(bp->pos[0] * 1e-7); /* HACK */
        cmp->B2 = bp->mass * cmp->bmax * cmp->bmax;
        Warning("cmp->B2 is zero, suspect identical particle positions\n");
    }

    if (cellsz > 0.0f) {
        /* This is an alternative bound on bmax, which is sometimes tighter */
        /* than the cumulative bound computed above. */
        cellsz *= (float)0.5;
        VxVVS(dx, = cellsz + fabs LPAREN cmp->center, -center, RPAREN);
        newbmax = sqrtf_fast(Dotx(dx, dx));
        cmp->bmax = (newbmax < cmp->bmax) ? newbmax : cmp->bmax;
    }
    hptr->ptr = cmp;
}

/* Tell the tree internals how big an object to copy */
int CellSz(void *p) {
    const cell *cp = p;

    if (mac->p4cut && (cp->daughters >= mac->p4cut)) {
        return sizeof(hexacell);
    } else if (mac->p2cut && (cp->daughters >= mac->p2cut)) {
        return sizeof(quadcell);
    } else {
        return sizeof(cell);
    }
}

/* Turn the ptr from a cofmdata to a cell. */
void *CellFromCofm(cofmdata *cmp) {
    cell *cp = NULL;
    quadcell *qcp = NULL;
    hexacell *hcp = NULL;
    cofmdata u;

    if (mac->p4cut && (cmp->ndaughters >= mac->p4cut)) {
        if (cellbuf->p4store_enable && cellbuf->p4store_used < cellbuf->p4store_len) {
            hcp = (hexacell *)cellbuf->p4store + cellbuf->p4store_used++;
            qcp = (quadcell *)hcp;
            cp = (cell *)hcp;
        } else {
            cp = ChnAlloc(&mac->tree->cell4chn);
            qcp = (quadcell *)cp;
            hcp = (hexacell *)cp;
        }
    } else if (mac->p2cut && (cmp->ndaughters >= mac->p2cut)) {
        cp = ChnAlloc(&mac->tree->cell2chn);
        qcp = (quadcell *)cp;
    } else {
        cp = ChnAlloc(&mac->tree->cellchn);
    }
    cp->mass = cmp->m;
    VV(cp->pos, = cmp->center);
    cp->bmax = cmp->bmax;
    cp->level = cmp->level;
    cp->daughters = cmp->ndaughters;
    cp->bptr = cmp->bptr;
    /* cp->R = 0.5f*cmp->sz; */
    if (mac->type == AREL_MAC) {
        float abs_rcrit;
        float rel_rcrit, rel_rcrit0;
        float B3, bmaxhalf, rcritmax;
        float B2 = cmp->B2;
        float massinv = recipf(cmp->m);
#if defined(QUAD) || defined(HEXA)
        double B4, B5, B6, rinv = 1.0 / cp->bmax;
#endif
        float ptol;

        if (mac->stol_max > 1.0f) {
            if (cmp->m <= min_sigma_m) {
                ptol = mac->this_tol;
            } else if (cmp->m > max_sigma_m) {
                ptol = mac->stol_max * mac->this_tol;
            } else {
                int idx = logf(cmp->m);
                float s = sigma_m[idx];
                if (s < mac->stol_max) {
                    ptol = s * mac->this_tol;
                } else {
                    ptol = mac->stol_max * mac->this_tol;
                }
            }
        } else {
            ptol = mac->this_tol / (1. + mac->ptol_boost * cmp->bmax / mac->bmax0);
        }
        bmaxhalf = cmp->bmax * (float)0.5;
        rcritmax = bmaxhalf + sqrtf_fast(bmaxhalf * bmaxhalf + sqrtf_fast((float)3. * B2 / ptol));
        if (!isfinite(rcritmax))
            Error("Bad rcritmax, q->bmax = %g, B2 = %g\n", cmp->bmax, B2);
        if (B2 == (float)0.0)
            Error("B2 is zero\n");
        B3 = B2 * sqrtf_fast(B2 * massinv); /* lower bound */
        if (!isfinite(B2) || !isfinite(B3) || !isfinite(cmp->bmax))
            Error("Bad value B2 = %g, B3 = %g, bmax = %g\n", B2, B3, cmp->bmax);
        a[0] = 2. * B3;
        a[1] = -3. * B2;
        a[2] = 0.;
        a[3] = ptol * cmp->bmax * cmp->bmax;
        a[4] = -2. * ptol * cmp->bmax;
        a[5] = ptol;
        abs_rcrit = rtnewt(5, rcrit_poly, rcritmax, .01 * rcritmax);
        abs_rcrit += 0.01 * rcritmax;

        rcritmax = cmp->bmax + sqrtf_fast((float)3. * mac->inv_rel_tol0 * B2 * massinv);
        a[0] = 2. * B3;
        a[1] = -3. * B2 + mac->rel_tol0 * cmp->bmax * cmp->bmax * cmp->m;
        a[2] = -2. * mac->rel_tol0 * cmp->bmax * cmp->m;
        a[3] = mac->rel_tol0 * cmp->m;
        rel_rcrit0 = rtnewt(3, rcrit_poly, rcritmax, .01 * rcritmax);
        rel_rcrit0 += 0.01 * rcritmax;

        rcritmax = cmp->bmax + sqrtf_fast((float)3. * mac->inv_rel_tol * B2 * massinv);
        a[0] = 2. * B3;
        a[1] = -3. * B2 + mac->rel_tol * cmp->bmax * cmp->bmax * cmp->m;
        a[2] = -2. * mac->rel_tol * cmp->bmax * cmp->m;
        a[3] = mac->rel_tol * cmp->m;
        rel_rcrit = rtnewt(3, rcrit_poly, rcritmax, .01 * rcritmax);
        rel_rcrit += 0.01 * rcritmax;

        /* rcrit is the least accurate of abs_rcrit or rel_rcrit */
        cp->rcrit = (abs_rcrit < rel_rcrit) ? abs_rcrit : rel_rcrit;
        /* if rel_rcrit0 is more accurate, then use it */
        if (rel_rcrit0 > cp->rcrit)
            cp->rcrit = rel_rcrit0;
        if (cmp->sz < 0.0f)
            cp->rcrit = 5.0f * mac->r0;
#ifdef QUAD
        memcpy(&u, cmp, sizeof(u));
        if (mac->subtract_background)
            mpole_add_cube(cmp->sz, -mac->rho0, &u);

        if (mac->p2cut && (cmp->ndaughters >= mac->p2cut)) {
            if (mac->subtract_background && cmp->ndaughters > 8) {
                B3 = (fabs(cmp->x3) + 3 * fabs(cmp->x2y) + 3 * fabs(cmp->xy2) + fabs(cmp->y3)
                      + 3 * fabs(cmp->x2z) + 6 * fabs(cmp->xyz) + 3 * fabs(cmp->y2z)
                      + 3 * fabs(cmp->xz2) + 3 * fabs(cmp->yz2) + fabs(cmp->z3))
                     / 9.0;
                B4 = 0.0;
            } else {
                B3 = cmp->bmax * (cmp->x2 + cmp->y2 + cmp->z2); /* upper bound */
                B4 = cmp->x4 + cmp->y4 + cmp->z4 + 2 * cmp->x2y2 + 2 * cmp->x2z2 + 2 * cmp->y2z2;
            }
            if (!isfinite(B4) || !isfinite(B3))
                Error("Bad B3 or B4, B2 = %g, massinv = %g\n", B2, massinv);
            rcritmax = abs_rcrit;
            a[0] = 3. * B4;
            a[1] = -4. * B3;
            a[2] = 0.;
            a[3] = 0.0;
            /* add another factor of r here? */
            a[4] = ptol * cmp->bmax * cmp->bmax;
            a[5] = -2. * ptol * cmp->bmax;
            a[6] = ptol;
            abs_rcrit = rtnewt(6, rcrit_poly, rcritmax, .01 * rcritmax);
            abs_rcrit += 0.01 * rcritmax;

            rcritmax = rel_rcrit0;
            a[0] = 3. * B4;
            a[1] = -4. * B3;
            a[2] = mac->rel_tol0 * cmp->bmax * cmp->bmax * cmp->m;
            a[3] = -2. * mac->rel_tol0 * cmp->bmax * cmp->m;
            a[4] = mac->rel_tol0 * cmp->m;
            rel_rcrit0 = rtnewt(4, rcrit_poly, rcritmax, .01 * rcritmax);
            rel_rcrit0 += 0.01 * rcritmax;

            rcritmax = rel_rcrit;
            a[0] = 3. * B4;
            a[1] = -4. * B3;
            a[2] = mac->rel_tol * cmp->bmax * cmp->bmax * cmp->m;
            a[3] = -2. * mac->rel_tol * cmp->bmax * cmp->m;
            a[4] = mac->rel_tol * cmp->m;
            rel_rcrit = rtnewt(4, rcrit_poly, rcritmax, .01 * rcritmax);
            rel_rcrit += 0.01 * rcritmax;

            qcp->rcrit_q = (abs_rcrit < rel_rcrit) ? abs_rcrit : rel_rcrit;
            if (rel_rcrit0 > qcp->rcrit_q)
                qcp->rcrit_q = rel_rcrit0;
            if (cmp->sz < 0.0f)
                qcp->rcrit_q = 5.0f * mac->r0;
            {
                double t;
#ifdef DIPOLE

                qcp->mass = u.m;
                t = rinv;
                qcp->qx = t * u.x;
                qcp->qy = t * u.y;
                qcp->qz = t * u.z;
#endif
                t = (u.x2 + u.y2 + u.z2) / 3.0;
                qcp->qxx = (u.x2 - t);
                qcp->qxy = u.xy;
                qcp->qyy = (u.y2 - t);
                qcp->qxz = u.xz;
                qcp->qyz = u.yz;
                t = rinv * rinv;
                qcp->qxx *= t;
                qcp->qxy *= t;
                qcp->qyy *= t;
                qcp->qxz *= t;
                qcp->qyz *= t;
            }
        }
#endif

#ifdef HEXA
        if (mac->p4cut && (cmp->ndaughters >= mac->p4cut)) {
            if (mac->subtract_background && cmp->ndaughters > 8) {
                B5 = (fabs(cmp->x5) + 5 * fabs(cmp->x4y) + 10 * fabs(cmp->x3y2)
                      + 10 * fabs(cmp->x2y3) + 5 * fabs(cmp->xy4) + fabs(cmp->y5)
                      + 5 * fabs(cmp->x4z) + 20 * fabs(cmp->x3yz) + 30 * fabs(cmp->x2y2z)
                      + 20 * fabs(cmp->xy3z) + 5 * fabs(cmp->y4z) + 10 * fabs(cmp->x3z2)
                      + 30 * fabs(cmp->x2yz2) + 30 * fabs(cmp->xy2z2) + 10 * fabs(cmp->y3z2)
                      + 10 * fabs(cmp->x2z3) + 20 * fabs(cmp->xyz3) + 10 * fabs(cmp->y2z3)
                      + 5 * fabs(cmp->xz4) + 5 * fabs(cmp->yz4) + fabs(cmp->z5))
                     / 27.0;
                B6 = 0.0;
            } else {
                B4 = (cmp->x4 + cmp->y4 + cmp->z4 + 2.0f * cmp->x2y2 + 2.0f * cmp->x2z2
                      + 2.0f * cmp->y2z2);
                B5 = cmp->bmax * B4; /* upper bound */
                B6 = cmp->x6 + cmp->y6 + cmp->z6 + 3.0f * cmp->x4y2 + 3.0f * cmp->x2y4
                     + 3.0f * cmp->x4z2 + 6.0f * cmp->x2y2z2 + 3.0f * cmp->y4z2 + 3.0f * cmp->x2z4
                     + 3.0f * cmp->y2z4;
            }
            if (!isfinite(B5) || !isfinite(B6))
                Error("Bad B5 or B6, B2 = %g, massinv = %g\n", B2, massinv);
            rcritmax = abs_rcrit;
            a[0] = 5. * B6;
            a[1] = -6. * B5;
            a[2] = 0.;
            a[3] = 0.;
            a[4] = 0.;
            a[5] = 0.;
            a[6] = ptol * cmp->bmax * cmp->bmax;
            a[7] = -2. * ptol * cmp->bmax;
            a[8] = ptol;
            abs_rcrit = rtnewt(8, rcrit_poly, rcritmax, .001 * rcritmax);
            abs_rcrit += 0.001 * rcritmax;

            rcritmax = rel_rcrit0;
            a[0] = 5. * B6;
            a[1] = -6. * B5;
            a[2] = 0.;
            a[3] = 0.;
            a[4] = mac->rel_tol0 * cmp->bmax * cmp->bmax * cmp->m;
            a[5] = -2. * mac->rel_tol0 * cmp->bmax * cmp->m;
            a[6] = mac->rel_tol0 * cmp->m;
            rel_rcrit0 = rtnewt(6, rcrit_poly, rcritmax, .001 * rcritmax);
            rel_rcrit0 += 0.001 * rcritmax;

            rcritmax = rel_rcrit;
            a[0] = 5. * B6;
            a[1] = -6. * B5;
            a[2] = 0.;
            a[3] = 0.;
            a[4] = mac->rel_tol * cmp->bmax * cmp->bmax * cmp->m;
            a[5] = -2. * mac->rel_tol * cmp->bmax * cmp->m;
            a[6] = mac->rel_tol * cmp->m;
            rel_rcrit = rtnewt(6, rcrit_poly, rcritmax, .001 * rcritmax);
            rel_rcrit += 0.001 * rcritmax;

            hcp->rcrit_h = (abs_rcrit < rel_rcrit) ? abs_rcrit : rel_rcrit;
            if (rel_rcrit0 > hcp->rcrit_h)
                hcp->rcrit_h = rel_rcrit0;
            if (cmp->sz < 0.0f)
                hcp->rcrit_h = 5.0f * mac->r0;
            {
                double t, tx, ty, tz, txx, txy, tyy, txz, tyz, tzz;
#ifdef DIPOLE
                hcp->mass = u.m;
                t = rinv;
                hcp->qx = t * u.x;
                hcp->qy = t * u.y;
                hcp->qz = t * u.z;
#endif
                t = (u.x2 + u.y2 + u.z2) / 3.0;
                hcp->qxx = u.x2 - t;
                hcp->qxy = u.xy;
                hcp->qyy = u.y2 - t;
                hcp->qxz = u.xz;
                hcp->qyz = u.yz;
                t = rinv * rinv;
                hcp->qxx *= t;
                hcp->qxy *= t;
                hcp->qyy *= t;
                hcp->qxz *= t;
                hcp->qyz *= t;
                tx = (u.x3 + u.xy2 + u.xz2) / 5.0;
                ty = (u.x2y + u.y3 + u.yz2) / 5.0;
                tz = (u.x2z + u.y2z + u.z3) / 5.0;
                hcp->qxxx = u.x3 - 3.0 * tx;
                hcp->qxxy = u.x2y - ty;
                hcp->qxyy = u.xy2 - tx;
                hcp->qyyy = u.y3 - 3.0 * ty;
                hcp->qxxz = u.x2z - tz;
                hcp->qxyz = u.xyz;
                hcp->qyyz = u.y2z - tz;
                t = rinv * rinv * rinv;
                hcp->qxxx *= t;
                hcp->qxxy *= t;
                hcp->qxyy *= t;
                hcp->qyyy *= t;
                hcp->qxxz *= t;
                hcp->qxyz *= t;
                hcp->qyyz *= t;
                txx = (u.x4 + u.x2y2 + u.x2z2) / 7.0;
                txy = (u.x3y + u.xy3 + u.xyz2) / 7.0;
                txz = (u.x3z + u.xy2z + u.xz3) / 7.0;
                tyy = (u.x2y2 + u.y4 + u.y2z2) / 7.0;
                tyz = (u.x2yz + u.y3z + u.yz3) / 7.0;
                tzz = (u.x2z2 + u.y2z2 + u.z4) / 7.0;
                t = 0.1 * (txx + tyy + tzz);
                hcp->qxxxx = u.x4 - 6.0 * (txx - t);
                hcp->qxxxy = u.x3y - 3.0 * txy;
                hcp->qxxyy = u.x2y2 - (txx + tyy - 2.0 * t);
                hcp->qxyyy = u.xy3 - 3.0 * txy;
                hcp->qyyyy = u.y4 - 6.0 * (tyy - t);
                hcp->qxxxz = u.x3z - 3.0 * txz;
                hcp->qxxyz = u.x2yz - tyz;
                hcp->qxyyz = u.xy2z - txz;
                hcp->qyyyz = u.y3z - 3.0 * tyz;
                txx = (cmp->x4 + cmp->x2y2 + cmp->x2z2) / 7.0;
                tyy = (cmp->x2y2 + cmp->y4 + cmp->y2z2) / 7.0;
                tzz = (cmp->x2z2 + cmp->y2z2 + cmp->z4) / 7.0;
                t = 0.1 * (txx + tyy + tzz);
                hcp->umass = cmp->m;
                hcp->uxxxx = cmp->x4 - 6.0 * (txx - t);
                hcp->uxxyy = cmp->x2y2 - (txx + tyy - 2.0 * t);
                hcp->uyyyy = cmp->y4 - 6.0 * (tyy - t);
                t = rinv * rinv * rinv * rinv;
                hcp->qxxxx *= t;
                hcp->qxxxy *= t;
                hcp->qxxyy *= t;
                hcp->qxyyy *= t;
                hcp->qyyyy *= t;
                hcp->qxxxz *= t;
                hcp->qxxyz *= t;
                hcp->qxyyz *= t;
                hcp->qyyyz *= t;
                hcp->uxxxx *= t;
                hcp->uxxyy *= t;
                hcp->uyyyy *= t;
            }
            Msgf(("Cell4: %s\n", PrintCellContents4(hcp)));
        }
#endif
    } else if (mac->type == BH_MAC)
        cp->rcrit = cmp->sz * mac->inv_tol;
    else if (mac->type == BMAX_MAC)
        cp->rcrit = cmp->bmax * mac->inv_tol;
    else
        Error("Bad MAC type (%d)\n", mac->type);

    return cp;
}

/* Use doubles here to avoid catastrophe from roundoff. */
static void rcrit_poly(int n, double r, double *value, double *deriv) {
    /* Do we care that a[2] and da[1] are zero? */
    double dp = 0.;
    double p;

    p = a[n];
    /* See pg. 149 of Numerical Rec. */
    /* We could unroll it... */
    while (n > 0) {
        dp = dp * r + p;
        p = p * r + a[--n];
    }
    if (!isfinite(p) || !isfinite(dp))
        Error("Bad p or dp, p = %g, dp = %g, n = %d\n", p, dp, n);
    *value = p;
    *deriv = dp;
}

/* From Numerical Recipes, rtnewt.c (modified) */
/* It's even more dangerous than the version that NR says is too */
/* dangerous to use...No checking of bounds.  We might just run off */
/* to infinity... */
#define JMAX 30

static double rtnewt(int n,
                     void (*funcd)(int, double, double *, double *),
                     double x1,
                     double xacc) {
    int j;
    double df, dx, f, rtn;

    rtn = x1;
    for (j = 1; j <= JMAX; j++) {
        (*funcd)(n, rtn, &f, &df);
        dx = f / df;
        Msg("rcrit", ("f(%g)=%g, dx=%g\n", rtn, f, -dx));
        rtn -= dx;
        if (fabs(dx) < xacc)
            return rtn;
    }
    SeriousWarning("Maximum number of iterations exceeded in RTNEWT, %lf %lf %lf %lf, %d\n",
                   x1,
                   xacc,
                   rtn,
                   dx,
                   n);
    return x1;
}

static void mpole_add_cube(double a, double rho, cofmdata *q) {
    double a3, a5, a7;

    a3 = a * a * a * rho;
    q->m += a3;

#ifdef QUAD
    a5 = a * a * a3;
    q->x2 += a5 / 12.;
    q->y2 += a5 / 12.;
    q->z2 += a5 / 12.;
#endif

#ifdef HEXA
    a7 = a * a * a5;
    q->x4 += a7 / 80.;
    q->x2y2 += a7 / 144.;
    q->y4 += a7 / 80.;
    q->x2z2 += a7 / 144.;
    q->y2z2 += a7 / 144.;
    q->z4 += a7 / 80.;
#endif
}

static void mpole_add_mono(cofmdata *p, body *q) {
    double x, y, z;
    double x2, y2, z2;
    double x3, y3, z3;
    double x4, y4, z4;
    double x5, y5, z5;
    double x6, y6, z6;
    double m;

    m = q->mass;

    p->m += m;

    x = -p->center[0] + q->pos[0];
    y = -p->center[1] + q->pos[1];
    z = -p->center[2] + q->pos[2];

    if (fabs(x) > mac->bmax0 || fabs(y) > mac->bmax0 || fabs(z) > mac->bmax0)
        Error("Bad wrap %.12g %.12g %.12g, R0 = %.12g\n",
              q->pos[0],
              q->pos[1],
              q->pos[2],
              mac->bmax0);

    p->x += m * x;
    p->y += m * y;
    p->z += m * z;

#ifdef QUAD
    x2 = x * x;
    y2 = y * y;
    z2 = z * z;

    p->x2 += m * x2;
    p->xy += m * x * y;
    p->xz += m * x * z;
    p->y2 += m * y2;
    p->yz += m * y * z;
    p->z2 += m * z2;

    x3 = x2 * x;
    y3 = y2 * y;
    z3 = z2 * z;

    p->x3 += m * x3;
    p->x2y += m * x2 * y;
    p->x2z += m * x2 * z;
    p->xy2 += m * x * y2;
    p->xyz += m * x * y * z;
    p->xz2 += m * x * z2;
    p->y3 += m * y3;
    p->y2z += m * y2 * z;
    p->yz2 += m * y * z2;
    p->z3 += m * z3;

    x4 = x3 * x;
    y4 = y3 * y;
    z4 = z3 * z;

    p->x4 += m * x4;
    p->x3y += m * x3 * y;
    p->x3z += m * x3 * z;
    p->x2y2 += m * x2 * y2;
    p->x2yz += m * x2 * y * z;
    p->x2z2 += m * x2 * z2;
    p->xy3 += m * x * y3;
    p->xy2z += m * x * y2 * z;
    p->xyz2 += m * x * y * z2;
    p->xz3 += m * x * z3;
    p->y4 += m * y4;
    p->y3z += m * y3 * z;
    p->y2z2 += m * y2 * z2;
    p->yz3 += m * y * z3;
    p->z4 += m * z4;
#endif

#ifdef HEXA
    x5 = x4 * x;
    y5 = y4 * y;
    z5 = z4 * z;

    p->x5 += m * x5;
    p->x4y += m * x4 * y;
    p->x4z += m * x4 * z;
    p->x3y2 += m * x3 * y2;
    p->x3yz += m * x3 * y * z;
    p->x3z2 += m * x3 * z2;
    p->x2y3 += m * x2 * y3;
    p->x2y2z += m * x2 * y2 * z;
    p->x2yz2 += m * x2 * y * z2;
    p->x2z3 += m * x2 * z3;
    p->xy4 += m * x * y4;
    p->xy3z += m * x * y3 * z;
    p->xy2z2 += m * x * y2 * z2;
    p->xyz3 += m * x * y * z3;
    p->xz4 += m * x * z4;
    p->y5 += m * y5;
    p->y4z += m * y4 * z;
    p->y3z2 += m * y3 * z2;
    p->y2z3 += m * y2 * z3;
    p->yz4 += m * y * z4;
    p->z5 += m * z5;

    x6 = x5 * x;
    y6 = y5 * y;
    z6 = z5 * z;

    p->x6 += m * x6;
    p->x4y2 += m * x4 * y2;
    p->x4z2 += m * x4 * z2;
    p->x2y4 += m * x2 * y4;
    p->x2y2z2 += m * x2 * y2 * z2;
    p->x2z4 += m * x2 * z4;
    p->y6 += m * y6;
    p->y4z2 += m * y4 * z2;
    p->y2z4 += m * y2 * z4;
    p->z6 += m * z6;
#endif
}

static void mpole_add(cofmdata *p, cofmdata *q) {
    double x, y, z;
    double x2, y2, z2;
    double x3, y3, z3;
    double x4, y4, z4;
    double x5, y5, z5;
    double x6, y6, z6;
    double m;

    m = q->m;

    p->m += m;

    x = q->center[0] - p->center[0];
    y = q->center[1] - p->center[1];
    z = q->center[2] - p->center[2];

    p->x += q->x + m * x;
    p->y += q->y + m * y;
    p->z += q->z + m * z;

#ifdef QUAD
    x2 = x * x;
    y2 = y * y;
    z2 = z * z;

    p->x2 += q->x2 + 2 * q->x * x + m * x2;
    p->xy += q->xy + q->y * x + q->x * y + m * x * y;
    p->xz += q->xz + q->z * x + q->x * z + m * x * z;
    p->y2 += q->y2 + 2 * q->y * y + m * y2;
    p->yz += q->yz + q->z * y + q->y * z + m * y * z;
    p->z2 += q->z2 + 2 * q->z * z + m * z2;

    x3 = x2 * x;
    y3 = y2 * y;
    z3 = z2 * z;

    p->x3 += q->x3 + 3 * q->x2 * x + 3 * q->x * x2 + m * x3;
    p->x2y += q->x2y + 2 * q->xy * x + q->y * x2 + q->x2 * y + 2 * q->x * x * y + m * x2 * y;
    p->x2z += q->x2z + 2 * q->xz * x + q->z * x2 + q->x2 * z + 2 * q->x * x * z + m * x2 * z;
    p->xy2 += q->xy2 + q->y2 * x + 2 * q->xy * y + 2 * q->y * x * y + q->x * y2 + m * x * y2;
    p->xyz += q->xyz + q->yz * x + q->xz * y + q->z * x * y + q->xy * z + q->y * x * z
              + q->x * y * z + m * x * y * z;
    p->xz2 += q->xz2 + q->z2 * x + 2 * q->xz * z + 2 * q->z * x * z + q->x * z2 + m * x * z2;
    p->y3 += q->y3 + 3 * q->y2 * y + 3 * q->y * y2 + m * y3;
    p->y2z += q->y2z + 2 * q->yz * y + q->z * y2 + q->y2 * z + 2 * q->y * y * z + m * y2 * z;
    p->yz2 += q->yz2 + q->z2 * y + 2 * q->yz * z + 2 * q->z * y * z + q->y * z2 + m * y * z2;
    p->z3 += q->z3 + 3 * q->z2 * z + 3 * q->z * z2 + m * z3;

    x4 = x3 * x;
    y4 = y3 * y;
    z4 = z3 * z;

    p->x4 += q->x4 + 4 * q->x3 * x + 6 * q->x2 * x2 + 4 * q->x * x3 + m * x4;
    p->x3y += q->x3y + 3 * q->x2y * x + 3 * q->xy * x2 + q->y * x3 + q->x3 * y + 3 * q->x2 * x * y
              + 3 * q->x * x2 * y + m * x3 * y;
    p->x3z += q->x3z + 3 * q->x2z * x + 3 * q->xz * x2 + q->z * x3 + q->x3 * z + 3 * q->x2 * x * z
              + 3 * q->x * x2 * z + m * x3 * z;
    p->x2y2 += q->x2y2 + 2 * q->xy2 * x + q->y2 * x2 + 2 * q->x2y * y + 4 * q->xy * x * y
               + 2 * q->y * x2 * y + q->x2 * y2 + 2 * q->x * x * y2 + m * x2 * y2;
    p->x2yz += q->x2yz + 2 * q->xyz * x + q->yz * x2 + q->x2z * y + 2 * q->xz * x * y
               + q->z * x2 * y + q->x2y * z + 2 * q->xy * x * z + q->y * x2 * z + q->x2 * y * z
               + 2 * q->x * x * y * z + m * x2 * y * z;
    p->x2z2 += q->x2z2 + 2 * q->xz2 * x + q->z2 * x2 + 2 * q->x2z * z + 4 * q->xz * x * z
               + 2 * q->z * x2 * z + q->x2 * z2 + 2 * q->x * x * z2 + m * x2 * z2;
    p->xy3 += q->xy3 + q->y3 * x + 3 * q->xy2 * y + 3 * q->y2 * x * y + 3 * q->xy * y2
              + 3 * q->y * x * y2 + q->x * y3 + m * x * y3;
    p->xy2z += q->xy2z + q->y2z * x + 2 * q->xyz * y + 2 * q->yz * x * y + q->xz * y2
               + q->z * x * y2 + q->xy2 * z + q->y2 * x * z + 2 * q->xy * y * z
               + 2 * q->y * x * y * z + q->x * y2 * z + m * x * y2 * z;
    p->xyz2 += q->xyz2 + q->yz2 * x + q->xz2 * y + q->z2 * x * y + 2 * q->xyz * z
               + 2 * q->yz * x * z + 2 * q->xz * y * z + 2 * q->z * x * y * z + q->xy * z2
               + q->y * x * z2 + q->x * y * z2 + m * x * y * z2;
    p->xz3 += q->xz3 + q->z3 * x + 3 * q->xz2 * z + 3 * q->z2 * x * z + 3 * q->xz * z2
              + 3 * q->z * x * z2 + q->x * z3 + m * x * z3;
    p->y4 += q->y4 + 4 * q->y3 * y + 6 * q->y2 * y2 + 4 * q->y * y3 + m * y4;
    p->y3z += q->y3z + 3 * q->y2z * y + 3 * q->yz * y2 + q->z * y3 + q->y3 * z + 3 * q->y2 * y * z
              + 3 * q->y * y2 * z + m * y3 * z;
    p->y2z2 += q->y2z2 + 2 * q->yz2 * y + q->z2 * y2 + 2 * q->y2z * z + 4 * q->yz * y * z
               + 2 * q->z * y2 * z + q->y2 * z2 + 2 * q->y * y * z2 + m * y2 * z2;
    p->yz3 += q->yz3 + q->z3 * y + 3 * q->yz2 * z + 3 * q->z2 * y * z + 3 * q->yz * z2
              + 3 * q->z * y * z2 + q->y * z3 + m * y * z3;
    p->z4 += q->z4 + 4 * q->z3 * z + 6 * q->z2 * z2 + 4 * q->z * z3 + m * z4;
#endif

#ifdef HEXA
    x5 = x4 * x;
    y5 = y4 * y;
    z5 = z4 * z;

    p->x5 += q->x5 + 5 * q->x4 * x + 10 * q->x3 * x2 + 10 * q->x2 * x3 + 5 * q->x * x4 + m * x5;
    p->x4y += q->x4y + 4 * q->x3y * x + 6 * q->x2y * x2 + 4 * q->xy * x3 + q->y * x4 + q->x4 * y
              + 4 * q->x3 * x * y + 6 * q->x2 * x2 * y + 4 * q->x * x3 * y + m * x4 * y;
    p->x4z += q->x4z + 4 * q->x3z * x + 6 * q->x2z * x2 + 4 * q->xz * x3 + q->z * x4 + q->x4 * z
              + 4 * q->x3 * x * z + 6 * q->x2 * x2 * z + 4 * q->x * x3 * z + m * x4 * z;
    p->x3y2 += q->x3y2 + 3 * q->x2y2 * x + 3 * q->xy2 * x2 + q->y2 * x3 + 2 * q->x3y * y
               + 6 * q->x2y * x * y + 6 * q->xy * x2 * y + 2 * q->y * x3 * y + q->x3 * y2
               + 3 * q->x2 * x * y2 + 3 * q->x * x2 * y2 + m * x3 * y2;
    p->x3yz += q->x3yz + 3 * q->x2yz * x + 3 * q->xyz * x2 + q->yz * x3 + q->x3z * y
               + 3 * q->x2z * x * y + 3 * q->xz * x2 * y + q->z * x3 * y + q->x3y * z
               + 3 * q->x2y * x * z + 3 * q->xy * x2 * z + q->y * x3 * z + q->x3 * y * z
               + 3 * q->x2 * x * y * z + 3 * q->x * x2 * y * z + m * x3 * y * z;
    p->x3z2 += q->x3z2 + 3 * q->x2z2 * x + 3 * q->xz2 * x2 + q->z2 * x3 + 2 * q->x3z * z
               + 6 * q->x2z * x * z + 6 * q->xz * x2 * z + 2 * q->z * x3 * z + q->x3 * z2
               + 3 * q->x2 * x * z2 + 3 * q->x * x2 * z2 + m * x3 * z2;
    p->x2y3 += q->x2y3 + 2 * q->xy3 * x + q->y3 * x2 + 3 * q->x2y2 * y + 6 * q->xy2 * x * y
               + 3 * q->y2 * x2 * y + 3 * q->x2y * y2 + 6 * q->xy * x * y2 + 3 * q->y * x2 * y2
               + q->x2 * y3 + 2 * q->x * x * y3 + m * x2 * y3;
    p->x2y2z += q->x2y2z + 2 * q->xy2z * x + q->y2z * x2 + 2 * q->x2yz * y + 4 * q->xyz * x * y
                + 2 * q->yz * x2 * y + q->x2z * y2 + 2 * q->xz * x * y2 + q->z * x2 * y2
                + q->x2y2 * z + 2 * q->xy2 * x * z + q->y2 * x2 * z + 2 * q->x2y * y * z
                + 4 * q->xy * x * y * z + 2 * q->y * x2 * y * z + q->x2 * y2 * z
                + 2 * q->x * x * y2 * z + m * x2 * y2 * z;
    p->x2yz2 += q->x2yz2 + 2 * q->xyz2 * x + q->yz2 * x2 + q->x2z2 * y + 2 * q->xz2 * x * y
                + q->z2 * x2 * y + 2 * q->x2yz * z + 4 * q->xyz * x * z + 2 * q->yz * x2 * z
                + 2 * q->x2z * y * z + 4 * q->xz * x * y * z + 2 * q->z * x2 * y * z + q->x2y * z2
                + 2 * q->xy * x * z2 + q->y * x2 * z2 + q->x2 * y * z2 + 2 * q->x * x * y * z2
                + m * x2 * y * z2;
    p->x2z3 += q->x2z3 + 2 * q->xz3 * x + q->z3 * x2 + 3 * q->x2z2 * z + 6 * q->xz2 * x * z
               + 3 * q->z2 * x2 * z + 3 * q->x2z * z2 + 6 * q->xz * x * z2 + 3 * q->z * x2 * z2
               + q->x2 * z3 + 2 * q->x * x * z3 + m * x2 * z3;
    p->xy4 += q->xy4 + q->y4 * x + 4 * q->xy3 * y + 4 * q->y3 * x * y + 6 * q->xy2 * y2
              + 6 * q->y2 * x * y2 + 4 * q->xy * y3 + 4 * q->y * x * y3 + q->x * y4 + m * x * y4;
    p->xy3z += q->xy3z + q->y3z * x + 3 * q->xy2z * y + 3 * q->y2z * x * y + 3 * q->xyz * y2
               + 3 * q->yz * x * y2 + q->xz * y3 + q->z * x * y3 + q->xy3 * z + q->y3 * x * z
               + 3 * q->xy2 * y * z + 3 * q->y2 * x * y * z + 3 * q->xy * y2 * z
               + 3 * q->y * x * y2 * z + q->x * y3 * z + m * x * y3 * z;
    p->xy2z2 += q->xy2z2 + q->y2z2 * x + 2 * q->xyz2 * y + 2 * q->yz2 * x * y + q->xz2 * y2
                + q->z2 * x * y2 + 2 * q->xy2z * z + 2 * q->y2z * x * z + 4 * q->xyz * y * z
                + 4 * q->yz * x * y * z + 2 * q->xz * y2 * z + 2 * q->z * x * y2 * z + q->xy2 * z2
                + q->y2 * x * z2 + 2 * q->xy * y * z2 + 2 * q->y * x * y * z2 + q->x * y2 * z2
                + m * x * y2 * z2;
    p->xyz3 += q->xyz3 + q->yz3 * x + q->xz3 * y + q->z3 * x * y + 3 * q->xyz2 * z
               + 3 * q->yz2 * x * z + 3 * q->xz2 * y * z + 3 * q->z2 * x * y * z + 3 * q->xyz * z2
               + 3 * q->yz * x * z2 + 3 * q->xz * y * z2 + 3 * q->z * x * y * z2 + q->xy * z3
               + q->y * x * z3 + q->x * y * z3 + m * x * y * z3;
    p->xz4 += q->xz4 + q->z4 * x + 4 * q->xz3 * z + 4 * q->z3 * x * z + 6 * q->xz2 * z2
              + 6 * q->z2 * x * z2 + 4 * q->xz * z3 + 4 * q->z * x * z3 + q->x * z4 + m * x * z4;
    p->y5 += q->y5 + 5 * q->y4 * y + 10 * q->y3 * y2 + 10 * q->y2 * y3 + 5 * q->y * y4 + m * y5;
    p->y4z += q->y4z + 4 * q->y3z * y + 6 * q->y2z * y2 + 4 * q->yz * y3 + q->z * y4 + q->y4 * z
              + 4 * q->y3 * y * z + 6 * q->y2 * y2 * z + 4 * q->y * y3 * z + m * y4 * z;
    p->y3z2 += q->y3z2 + 3 * q->y2z2 * y + 3 * q->yz2 * y2 + q->z2 * y3 + 2 * q->y3z * z
               + 6 * q->y2z * y * z + 6 * q->yz * y2 * z + 2 * q->z * y3 * z + q->y3 * z2
               + 3 * q->y2 * y * z2 + 3 * q->y * y2 * z2 + m * y3 * z2;
    p->y2z3 += q->y2z3 + 2 * q->yz3 * y + q->z3 * y2 + 3 * q->y2z2 * z + 6 * q->yz2 * y * z
               + 3 * q->z2 * y2 * z + 3 * q->y2z * z2 + 6 * q->yz * y * z2 + 3 * q->z * y2 * z2
               + q->y2 * z3 + 2 * q->y * y * z3 + m * y2 * z3;
    p->yz4 += q->yz4 + q->z4 * y + 4 * q->yz3 * z + 4 * q->z3 * y * z + 6 * q->yz2 * z2
              + 6 * q->z2 * y * z2 + 4 * q->yz * z3 + 4 * q->z * y * z3 + q->y * z4 + m * y * z4;
    p->z5 += q->z5 + 5 * q->z4 * z + 10 * q->z3 * z2 + 10 * q->z2 * z3 + 5 * q->z * z4 + m * z5;

    x6 = x5 * x;
    y6 = y5 * y;
    z6 = z5 * z;

    p->x6 += q->x6 + 6 * q->x5 * x + 15 * q->x4 * x2 + 20 * q->x3 * x3 + 15 * q->x2 * x4
             + 6 * q->x * x5 + m * x6;
    p->x4y2 += q->x4y2 + 4 * q->x3y2 * x + 6 * q->x2y2 * x2 + 4 * q->xy2 * x3 + q->y2 * x4
               + 2 * q->x4y * y + 8 * q->x3y * x * y + 12 * q->x2y * x2 * y + 8 * q->xy * x3 * y
               + 2 * q->y * x4 * y + q->x4 * y2 + 4 * q->x3 * x * y2 + 6 * q->x2 * x2 * y2
               + 4 * q->x * x3 * y2 + m * x4 * y2;
    p->x4z2 += q->x4z2 + 4 * q->x3z2 * x + 6 * q->x2z2 * x2 + 4 * q->xz2 * x3 + q->z2 * x4
               + 2 * q->x4z * z + 8 * q->x3z * x * z + 12 * q->x2z * x2 * z + 8 * q->xz * x3 * z
               + 2 * q->z * x4 * z + q->x4 * z2 + 4 * q->x3 * x * z2 + 6 * q->x2 * x2 * z2
               + 4 * q->x * x3 * z2 + m * x4 * z2;
    p->x2y4 += q->x2y4 + 2 * q->xy4 * x + q->y4 * x2 + 4 * q->x2y3 * y + 8 * q->xy3 * x * y
               + 4 * q->y3 * x2 * y + 6 * q->x2y2 * y2 + 12 * q->xy2 * x * y2 + 6 * q->y2 * x2 * y2
               + 4 * q->x2y * y3 + 8 * q->xy * x * y3 + 4 * q->y * x2 * y3 + q->x2 * y4
               + 2 * q->x * x * y4 + m * x2 * y4;
    p->x2y2z2 += q->x2y2z2 + 2 * q->xy2z2 * x + q->y2z2 * x2 + 2 * q->x2yz2 * y
                 + 4 * q->xyz2 * x * y + 2 * q->yz2 * x2 * y + q->x2z2 * y2 + 2 * q->xz2 * x * y2
                 + q->z2 * x2 * y2 + 2 * q->x2y2z * z + 4 * q->xy2z * x * z + 2 * q->y2z * x2 * z
                 + 4 * q->x2yz * y * z + 8 * q->xyz * x * y * z + 4 * q->yz * x2 * y * z
                 + 2 * q->x2z * y2 * z + 4 * q->xz * x * y2 * z + 2 * q->z * x2 * y2 * z
                 + q->x2y2 * z2 + 2 * q->xy2 * x * z2 + q->y2 * x2 * z2 + 2 * q->x2y * y * z2
                 + 4 * q->xy * x * y * z2 + 2 * q->y * x2 * y * z2 + q->x2 * y2 * z2
                 + 2 * q->x * x * y2 * z2 + m * x2 * y2 * z2;
    p->x2z4 += q->x2z4 + 2 * q->xz4 * x + q->z4 * x2 + 4 * q->x2z3 * z + 8 * q->xz3 * x * z
               + 4 * q->z3 * x2 * z + 6 * q->x2z2 * z2 + 12 * q->xz2 * x * z2 + 6 * q->z2 * x2 * z2
               + 4 * q->x2z * z3 + 8 * q->xz * x * z3 + 4 * q->z * x2 * z3 + q->x2 * z4
               + 2 * q->x * x * z4 + m * x2 * z4;
    p->y6 += q->y6 + 6 * q->y5 * y + 15 * q->y4 * y2 + 20 * q->y3 * y3 + 15 * q->y2 * y4
             + 6 * q->y * y5 + m * y6;
    p->y4z2 += q->y4z2 + 4 * q->y3z2 * y + 6 * q->y2z2 * y2 + 4 * q->yz2 * y3 + q->z2 * y4
               + 2 * q->y4z * z + 8 * q->y3z * y * z + 12 * q->y2z * y2 * z + 8 * q->yz * y3 * z
               + 2 * q->z * y4 * z + q->y4 * z2 + 4 * q->y3 * y * z2 + 6 * q->y2 * y2 * z2
               + 4 * q->y * y3 * z2 + m * y4 * z2;
    p->y2z4 += q->y2z4 + 2 * q->yz4 * y + q->z4 * y2 + 4 * q->y2z3 * z + 8 * q->yz3 * y * z
               + 4 * q->z3 * y2 * z + 6 * q->y2z2 * z2 + 12 * q->yz2 * y * z2 + 6 * q->z2 * y2 * z2
               + 4 * q->y2z * z3 + 8 * q->yz * y * z3 + 4 * q->z * y2 * z3 + q->y2 * z4
               + 2 * q->y * y * z4 + m * y2 * z4;
    p->z6 += q->z6 + 6 * q->z5 * z + 15 * q->z4 * z2 + 20 * q->z3 * z3 + 15 * q->z2 * z4
             + 6 * q->z * z5 + m * z6;
#endif
}
