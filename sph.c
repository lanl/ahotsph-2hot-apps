#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "physics_sph.h"
#include "vop.h"
#include "fastflpt.h"
#include "timers.h"
#include "error.h"
#include "singlio.h"
#include "nrutil.h"
#include "units.h"
#include "cool.h"

#ifndef M_1_PI
#define	M_1_PI 0.31830988618379067154
#endif
#define NKERNEL_TABLE 80000
#define MAX_INDEX (NKERNEL_TABLE+2)

#define NO_UPDATE 2

Counter_t SPHCnt, SPHrej, nbrMACCnt;

static float dvtable; /* == 0.0001 ... */
static float invdvtable; /* == 10000.0 ... */
static float cnormk;
static float wij[MAX_INDEX];
static float grwij[MAX_INDEX];
static float fmass[MAX_INDEX];
static float fpoten[MAX_INDEX];
static float Gamma;/* = (float)(5.0/3.0);*/
static float alpha;/* = (float)1.0;*/
static float beta;/* = (float)2.5;*/
static float epsil;/* = (float)1e-2;*/
static float heatf1;/* = (float)1.0;*/
static int ndim;
static int Nobj;
static int add_offset;
static float offset[NDIM];
static float voffset[NDIM];
static void (*bodyfunc)(SinkSPH *sink, hcell **src_vec, int *res, int n);
static void (*cellfunc)(SinkSPH *sink, hcell **src_vec, int *res, int n);

//extern int do_diffusion;
//extern int do_cooling;
//extern int do_burning;

    void
SetSPHOffset(float *off, float *voff)
{
    VV(offset, = off);
    VV(voffset, = voff);
    add_offset = 1;
}

    void
UnSetSPHOffset(void)
{
    VS(offset, = 0.0);
    VS(voffset, = 0.0);
    add_offset = 0;
}

    void 
InheritSPH(const SinkSPH *from, SinkSPH *to, hcell *pp)
{
    if( to == NULL ){
        SPHbody *bp = pp->ptr;
        if (from->isbody == NO_UPDATE)
            return;
        /* Must accumulate for periodic BC to work */
        /* Must initialize to zero appropriately */
        bp->rho += from->rho;
        bp->du += from->du;       /* Why am I doing this? */
        bp->du_r += from->du_r;   /* Do I need to do anything else? */
        bp->drho_dt += from->drho_dt;
        bp->udot += from->udot;
        bp->nbrs += from->nbrs;
        VV(bp->acc, += from->M1);
        VV(bp->lvel, += from->lvel);
        bp->nterms += from->nterms;
        bp->min_nbr_dt = from->min_nbr_dt;   /* set dt to min of nbrs dt */
        if (from->interactions != Nobj)
            Error("Ninteract is %d, should be %d\n", from->interactions, Nobj);
        return;
    }

    if (Sub_Flags(pp)) {
        /* Stuff needed for cell-cell neighbor evaluation */
        SPHcell *cp = pp->ptr;
        VV(to->pos, = cp->pos);
        to->extent = cp->bmax + cp->lap;
        to->isbody = 0;
    } else {
        /* Stuff needed for above, plus physics info */
        SPHbody *bp = pp->ptr;
        if (!SPH_need_update(bp)) {
            to->isbody = NO_UPDATE;
            return;
        }
        VV(to->pos, = bp->pos);
        to->extent = bp->h;
        to->h = bp->h;
        to->isbody = 1;
        VV(to->vel, = bp->vel);
        to->pr = bp->pr;
        to->rho_est = bp->rho_est;
        to->mass = bp->mass;
        to->vsound = bp->vsound;
        to->u = bp->u;
        to->rho = (float)0.0;
        to->drho_dt = (float)0.0;
        to->udot = (float)0.0;
        VS(to->lvel, = (float)0.0);
        to->nterms = 1;
        to->nbrs = 0;
        VS(to->M1, = (float)0.0);
        to->min_nbr_dt = 1e30;
        /* Diffusion quantities */
        /* Can some of these be set to zero here?  Like udot above? */
        to->temp = bp->temp;
        to->du = bp->du;
        to->u_r = bp->u_r;
        to->du_r = bp->du_r;  /* Or = (float)0.0; ? */
        to->D = bp->D;
    }

    if (add_offset) {
        VV(to->pos, += offset);
        VV(to->vel, += voffset);
    }

    if (from) {
        to->interactions = from->interactions;
    } else {
        to->interactions = 0;
    }

}

    void
SPHgate(SinkSPH *sink, hcell **src_vec, int *result, int n)
{
    int i;

    if (sink->isbody == NO_UPDATE) {
        for (i = 0; i < n; i++) result[i] = MAC_ACCEPT;
    } else if (sink->isbody)
        bodyfunc(sink, src_vec, result, n);
    else
        cellfunc(sink, src_vec, result, n);
}

    void
nbrMAC(SinkSPH *sink, hcell **source_vec, int *result, int n)
{
    int i;

    for (i = 0; i < n; i++) {	/* This is about 20% less efficient */
        result[i] = MAC_SPLIT_SINK;
    }
}


    void
macRho(SinkSPH *sink, hcell **source_vec, int *result, int n)
{
    const float extent_sink = sink->extent;
    VxdV(const float pos_sink, = sink->pos);
    const float h = sink->h;
    Vxd(float r);
    Vxd(float dv);
    float extent_src;
    int daughters;
    int i;
    SPHbody *bp = 0;
    float projv;
    float v2;
    float rij;
    float wtij, grwtij;
    float hmean11, hmean21;
    float dxx, dwdx, dgrwdx;
    int index;
    int nbrs = 0;
    float rhoi = (float)0.0;
    float divvi = (float)0.0;
    float dr2;
    int interactions = 0;

    for (i = 0; i < n; i++) {
        const hcellptr source = source_vec[i];

        if (Sub_Flags(source)) {
            const SPHcell *cp = source->ptr;
            VxV(r, = cp->pos);	
            extent_src = cp->bmax + cp->lap;
            daughters = cp->daughters;
        } else {
            bp = source->ptr;
            VxV(r, = bp->pos);
            extent_src = bp->h;
            daughters = 1;
        }

        VxVx(r, -= pos_sink);	/* 8 flops */
        dr2 = Dotx(r, r); /* == 400.0000000032623 */

        if ((rij = sqrtf_fast(dr2)) > extent_src + extent_sink
                || dr2 == (float)0.0) {
            goto accept;
        } else if (daughters != 1) {
            goto failed;
        }

        hmean11 = (float)2.0 / (h + bp->h);
        hmean21 = hmean11 * hmean11; /* == 0.01 */

        v2 = dr2 * hmean21;
        index = v2 * invdvtable;
        if (index >= MAX_INDEX) Error("Index too large\n");
        dxx = v2 - index * dvtable;
        dwdx = (wij[index+1] - wij[index]) * invdvtable;
        wtij = (wij[index] + dwdx * dxx ) * hmean21 * hmean11;
        /* if (wtij < (float)0.0) Error("Negative wtij (macRho) = %g\n", wtij); */
        dgrwdx = (grwij[index+1] - grwij[index]) * invdvtable;
        grwtij = (grwij[index] + dgrwdx * dxx) * hmean21 * hmean21;

        rhoi += bp->mass * wtij;

        /* velocity divergence times density */
        VxVV(dv, = bp->vel, - sink->vel);
        projv = grwtij * Dotx(dv, r) * recipsqrtf(dr2);
        divvi -= bp->mass * projv;

        nbrs++;
accept:
        interactions += daughters;
        result[i] = MAC_ACCEPT;
        continue;
failed:
        result[i] = MAC_SPLIT_SRC;
    }

    sink->interactions += interactions;
    sink->rho += rhoi;
    sink->nbrs += nbrs;
    sink->drho_dt -= divvi;
}


    void
macSPH(SinkSPH *sink, hcell **source_vec, int *result, int n)
{
    const float extent_sink = sink->extent;
    VxdV(const float pos_sink, = sink->pos);
    VxdV(const float v, = sink->vel);
    const float h = sink->h;
    const float pro2 = sink->pr / (sink->rho_est * sink->rho_est);
    const float mass = sink->mass;
    const float rho_est = sink->rho_est;
    const float vsound = sink->vsound;
    const float u = sink->u;
    Vxd(float r);
    Vxd(float f);
    Vxd(float dv);
    Vxd(float smv);
    float min_nbr_dt = sink->min_nbr_dt;
    float extent_src;
    int daughters;
    int i;
    SPHbody *bp = 0;
    float hmean11, hmean21;
    int index;
    int nbrs = 0;
    float rhoi = (float)0.0;
    float divvi = (float)0.0;
    float dr2;
    Vxd(float runi);
    float dq = (float)0.0;
    float vv, vv2;
    float dxx, dwdx, wtij, dgrwdx, grwtij;
    float rapm, robar1, grpm, wpm;
    float poro2;
    float projv, vsbar, est_divv, t12;
    float rij, rij1;
    float clight;
    int interactions = 0;

    VxS(f, = (float)0.0);
    VxS(smv, = (float)0.0);
    for (i = 0; i < n; i++) {
        const hcellptr source = source_vec[i];

        if (Sub_Flags(source)) {
            const SPHcell *cp = source->ptr;
            VxV(r, = cp->pos);	
            extent_src = cp->bmax + cp->lap;
            daughters = cp->daughters;
        } else {
            bp = source->ptr;
            VxV(r, = bp->pos);
            extent_src = bp->h;
            daughters = 1;
            if (bp->rho_est <= (float)0.0)
                Error("RhoEst is <= 0 for %s\n", hcellPrint(source));
        }

        VxVx(r, -= pos_sink);	/* 8 flops */
        dr2 = Dotx(r, r);

        sink->nterms += 1;

        if ((rij = sqrtf_fast(dr2)) > extent_src + extent_sink
                || dr2 == (float)0.0) {
            goto accept;
        } else if (daughters != 1) {
            goto failed;
        }

        if (bp->dt_next < min_nbr_dt) min_nbr_dt = bp->dt_next;

        hmean11 = (float)2.0 / (h + bp->h);
        hmean21 = hmean11 * hmean11;

        vv2 = dr2 * hmean21;	/* v2 and v renamed to avoid conflict */
        index = vv2 * invdvtable;
        if (index >= MAX_INDEX) Error("Index too large\n");
        vv = rij * hmean11;
        dxx = vv2 - index * dvtable;
        dwdx = (wij[index+1] - wij[index]) * invdvtable;
        wtij = (wij[index] + dwdx * dxx) * hmean21 * hmean11;
        /* if (wtij < (float)0.0) Error("Negative wtij (macSPH) = %g\n", wtij); */
        dgrwdx = (grwij[index+1] - grwij[index]) * invdvtable;
        grwtij = (grwij[index] + dgrwdx * dxx) * hmean21 * hmean21;

        rapm = mass / bp->mass;
        robar1 = (float)2.0 / (rho_est + bp->rho_est);
        grpm = bp->mass * grwtij;
        wpm = bp->mass * wtij;

        poro2 = grpm * (pro2 + bp->pr / (bp->rho_est * bp->rho_est));
        rij1 = (float)1.0 / rij;
        VxVx(runi, = rij1 * r);
        VxVx(f, += poro2 * runi);
        VxVVx(dv, = bp->vel, - v);
        VxVx(smv, += robar1 * wpm * dv);
        projv = Dotx(dv, runi);

        rhoi += bp->mass * wtij;
        divvi -= bp->mass * grwtij * projv;

        /* artificial viscosity and energy dissipation */
        vsbar = (float)0.5 * (vsound + bp->vsound);
        if (projv < (float)0.0 && alpha != (float)0.0) {
            est_divv = projv * vv / (vv2 + epsil);
            t12 = grpm * est_divv * (beta * est_divv - alpha * vsbar) * robar1;
            VxVx(f, += t12 * runi);
            dq += t12 * projv;
        }
        /* artificial heat conduction */
        if (heatf1 != (float)0.0) {
            float qmean = heatf1 * vsbar / hmean11;
            t12 = grpm * qmean * (bp->u - u) * robar1;
            dq -= t12;
        }
        /* flux-limited diffusion */
        if (params.do_diffusion)
            if (grpm < 0.0) {  /* What does this condition really mean? */
                float Dmeanr = 2.0*rij1*sink->D/(sink->D+bp->D)*bp->D;
                clight = C_LIGHT*tdivlCF* rij1;

                sink->du_r += ( (clight < Dmeanr) ? clight : Dmeanr ) *
                    (sink->u_r - bp->u_r) * grpm / bp->rho_est;
            }

        nbrs++;
accept:
        IncrCounter(&SPHrej);
        interactions += daughters;
        result[i] = MAC_ACCEPT;
        continue;
failed:
        IncrCounter(&SPHrej);
        result[i] = MAC_SPLIT_SRC;
    }
    sink->interactions += interactions;
    sink->rho += rhoi;
    sink->drho_dt -= divvi;
    sink->udot += (float)0.5 * dq;
    VVx(sink->M1, += f);
    VVx(sink->lvel, += smv);
    sink->nbrs += nbrs;
    sink->nterms += nbrs*8;
    sink->min_nbr_dt = min_nbr_dt;
}

    void
SetSPH(float visc_alpha, float visc_beta, float visc_epsilon, float heat_f1,
        float eos_gamma, int gnobj,  void bfunc(), void cfunc())
{
    Nobj = gnobj;
    alpha = visc_alpha;
    beta = visc_beta;
    epsil = visc_epsilon;
    heatf1 = heat_f1;
    Gamma = eos_gamma;
    bodyfunc = bfunc;
    cellfunc = cfunc;

}

    void
SPH_setup(int dim, int ncoef1, double *wcoef1, int ncoef2, double *wcoef2)
{
    double v2max;
    double v, v2;
    double w, dw;
    double ddvtable;
    double dm, dm1;
    int i, i1, j;

    ndim = dim;


    /* Do any diffusion-specific initialization here: */
    /* Set rmax, luminosity, etc. */


    /* maximum interaction length and step size */
    v2max = 4.0;
    ddvtable = v2max / NKERNEL_TABLE;
    invdvtable = 1.0 / ddvtable;
    dvtable = ddvtable;

    /* normalization constant */
    if (ndim == 3) 
        cnormk = M_1_PI;
    else if (ndim == 2)
        cnormk = M_1_PI*10.0/7.0;
    else if (ndim == 1)
        cnormk = 2.0/3.0;
    else
        Error("Bad ndim in sph_ktable\n");

    /* build tables */
    /* a) v less than 1 */

    wij[0] = cnormk*wcoef1[0];
    grwij[0] = 0.0;
    fmass[0] = 0.0;
    fpoten[0] = 0.0;
    for (j = 0; j <= ncoef1-1; j++) 
        fpoten[0] += wcoef1[j]/(j+2.0) +
            wcoef2[j]/(j+2.) * (pow(2.0, j+2.0)-1.0);
    fpoten[0]*=4.*M_PI*cnormk;

    i1 = 1.0 / ddvtable;
    for (i = 1; i <= i1; i++) {
        v2 = i * ddvtable;
        v = sqrt(v2);
        w = wcoef1[ncoef1-1];
        for (j = ncoef1-1; j >= 1 ; j--)
            w = w * v + wcoef1[j-1];
        dw = (ncoef1 - 1.0) * wcoef1[ncoef1-1];
        for (j = ncoef1-2; j >= 1; j--)
            dw = dw * v + j * wcoef1[j];
        wij[i] = cnormk * w;
        grwij[i] = cnormk * dw;

        /* Enclosed mass now for m=1, h=1*/
        dm = wcoef1[ncoef1-1]/(ncoef1-1+3);
        for (j = ncoef1-2; j >= 0; j--)
            dm = dm * v + wcoef1[j]/(j+3.);
        fmass[i]=4*M_PI*cnormk*v*v*v*dm;

        /* Potential now: fpoten is potential for G=1, m=1, h=1*/
        fpoten[i]=(double)0.;
        for (j = 0;j <= ncoef1-1; j++) 
            fpoten[i]+=wcoef1[j]/(j+3.)*pow(v,(j+2.))+
                wcoef1[j]/(j+2.)*(1-pow(v,(j+2.)))+
                wcoef2[j]/(j+2.)*(pow(2.,(j+2.))-1.);
        fpoten[i]*=4.*M_PI*cnormk;
    }

    v=1.;
    dm1=wcoef2[ncoef2-1]/(ncoef2-1+3.);
    for (j=ncoef2-2; j>=0;j--)  
        dm1=dm1*v + wcoef2[j]/(j+3.);      
    dm1=fmass[i1]-4*M_PI*dm1*cnormk*v*v*v;

    /*  b) v greater than 1 */
    for (i = i1+1; i <= NKERNEL_TABLE; i++) {
        v2 = i * ddvtable;
        v = sqrt(v2);
        w = wcoef2[ncoef2-1];
        for (j = ncoef2-1; j >= 1 ; j--)
            w = w * v + wcoef2[j-1];
        dw = (ncoef2 - 1.0) * wcoef2[ncoef2-1];
        for (j = ncoef2-2; j >= 1; j--)
            dw = dw * v + j * wcoef2[j];
        wij[i] = cnormk * w;
        grwij[i] = cnormk * dw;

        /* Enclosed mass now */
        dm = wcoef2[ncoef2-1]/(ncoef2-1+3);
        for (j = ncoef2-2; j >= 0; j--)
            dm = dm * v + wcoef2[j]/(j+3.);
        fmass[i]=4*M_PI*cnormk*v*v*v*dm+dm1;

        /* Potential now: fpoten is potential for G=1, m=1, h=1*/
        fpoten[i]=(double) 0.;
        for (j=0; j <= ncoef2-1; j++)
            fpoten[i]+=1./v*wcoef1[j]/(j+3.)+
                1./v*wcoef2[j]/(j+3.)*(pow(v,(j+3.))-1.)+
                wcoef2[j]/(j+2.)*(pow(2.,(j+2.))-pow(v,(j+2.)));
        fpoten[i]*=4.*M_PI*cnormk;
    }

    /* Make sure the mass is perfectly normalized, don't generate mass */
    for (i = 0; i <= NKERNEL_TABLE; i++) {
        fmass[i]/=fmass[NKERNEL_TABLE];
    }

    /* For interpolation of last table entry */
    wij[NKERNEL_TABLE+1] = grwij[NKERNEL_TABLE+1] = 0.0;
    fpoten[NKERNEL_TABLE+1]=(double)0.5;
    fmass[NKERNEL_TABLE+1]=(double)1.0;
}

#include "Msgs.h"
double eos_n, eos_u;


    void
update_final(SPHbody *btab, int nobj, int Gridpts, const int Nel, float dt, int *limit_high, int *limit_low, int rank, float tpos, float R0)
{
    SPHbody *p;
    int i,j,k; /*coupla indices for loops*/
    double u, n, lcool, deltah, dt_cool;
    double *molfrac; /*float or double?? */
    double dt_tot,udot_tot,dt_sub,dt_save,udot,frac=0.05, minfrac=0.001;
    double temp,rho, dt_cgs, ndens = 0., ne=0.;
    double tlo, tup, mfp, radius;
    int decr,notprinted;
    long cycles=0, countc; 
    static long cycled = 0;
    int temp_ok, partid;

    /*molfrac = (double *)malloc( (NNW+1) * sizeof(double) );*/

    tlo = 1.0e0;
    tup = 2.5e11;

    notprinted = 1;

    for (p = btab; p < btab+nobj; p++) {
        /********* checking udot and hdot ***********/
        if (!SPH_need_update(p)) continue;
        VV(p->acc, += p->grav_acc);
        /* Changed cnormk to wij[0] to allow for non-standard kernels; thanks Steven */
        p->rho += wij[0] * p->mass / (p->h * p->h * p->h); 
        p->hdot = (float)(-1.0/3.0) * p->h * p->drho_dt / p->rho;
        if (p->hdot * dt > p->h) {
            SeriousWarning("Hdot limit (high)\n%s\n", PrintSPHBodyContents(p));
            p->hdot = p->h/dt;
            ++*limit_high;
        }
        if (p->hdot * dt < -0.5*p->h) {
            SeriousWarning("Hdot limit (low)\n%s\n", PrintSPHBodyContents(p));
            p->hdot = -0.5*p->h/dt;
            ++*limit_low;
        }

        p->udot += p->drho_dt * p->pr / (p->rho * p->rho) + 
            ( (params.do_diffusion) ? (p->du_r/p->rho) /* Diffusion */
              : 0.0 );

        if (!isfinite(p->udot)) 
            Error("%d Bad value for udot: drho/dt=%.4E pr=%.4E rho=%.4E\n",
                    p->ident,p->drho_dt,p->pr,p->rho);

        /* Are these limits appropriate? */
        /* Does this enforce the Courant limit correctly with diffusion? */
        if ( (p->udot * dt > p->u) && !(p->ident & (1<<30)) ) {
            p->udot = p->u/dt;
            ++*limit_high;
        }
        if ( (p->udot * dt < -0.333*p->u) && !(p->ident & (1<<30)) ) {
            p->udot = -0.333*p->u/dt;
            ++*limit_low;
        }

		if (!params.poly_eos) {
        /************* update T, eos_n, eos_u *************/
		if (params.do_cooling || params.do_burning) {
	        temp_ok = prep_cool_burn(p, 1.e1, 2.5e11, Gridpts, Nel, 0);
		} else {
		    /* Calculate temperature from u, then "create" photons (a*T^4) */
		    eos_n = (double)(p->rho) * ivlenCF3 * massCF;
			eos_n /= (double) 28 * MH;
		    eos_u = ((double)(p->u))*((double)(p->rho));
			eos_u *= massCF * ivlenCF * ivtimeCF2; 
	
		    /* Figure out good upper and lower limits for temp */
		    p->temp = newtraph(4.0e3, 1.5e11, eos_u*1.0e-6, uvst, duvst);
		}

        /********** do the burning ***********/
        if(params.do_burning && temp_ok)
            deltah = burning(p, dt, rank);

        /********** do the cooling ***********/
        if (params.do_cooling && (tpos >= 2.62980e6*ivtimeCF) ) {/*(tstar > 0.5 * M_PI*1.e7 / t)) */
            dt_cool = cooling(p, dt, frac, Gridpts, Nel, &notprinted);
            /*
               if ( iter%interval == 0) {
               print_ionization_state();
               }
               */
        }
		}
    }
}


    void
update_intermediate(SPHbody *btab, int nobj, int Gridpts, const int Nel, float dt_last, int flag, int *limit, float R0)
{
    float kes, kff;  /* Opacities (Thomson, free-free) */
    float pgas, prad;
    float ne, rho, max_rad, r2;
    double P_ratio, Gammai; 
    double tlo = 1.e0, tup = 2.5e11;
    int j, temp_ok;
    SPHbody *p;

    max_rad = 0.95*R0*R0;

    for (p = btab; p < btab+nobj; p++) {
        if (!SPH_need_update(p)) continue;
        if (flag) p->rho_est = p->rho + p->drho_dt * dt_last;
        else p->rho_est = p->rho;
        if (p->rho_est <= (float)0.0) 
            Error("Rho_est is 0\n%s\n", PrintSPHBodyContents(p));

        /* keep these in cgs-units */
		if (params.do_cooling || params.do_burning) 
	        temp_ok = prep_cool_burn(p, 1.e1, 1.0e11, Gridpts, Nel, 1);
		else if (params.poly_eos) {
			p->pr = POLY_EOS_K * pow (p->rho_est, POLY_EOS_GAMMA);
			p->vsound = sqrtf_fast (p->pr / p->rho_est);
		}
		else {
		    /* Calculate temperature from u, then "create" photons (a*T^4) */
		    eos_n = (double)(p->rho_est) * ivlenCF3 * massCF;
			eos_n /= (double) 28 * MH;
		    eos_u = ((double)(p->u))*((double)(p->rho_est));
			eos_u *= massCF * ivlenCF * ivtimeCF2; 
	
		    /* Figure out good upper and lower limits for temp */
		    p->temp = newtraph(4.0e3, 1.5e11, eos_u*1.0e-6, uvst, duvst);
			p->mfp = p->h;
		}

        /* calculate the total pressure by calculating the respective 
           contributions of gas and radiation pressure */
        pgas = prad = 0.0;
        pgas = (double)(eos_n * K_BOLTZ * p->temp);
        if( 20.*p->h > p->mfp) 
            prad = (double)(0.33333333333*A_RAD * p->temp*p->temp*p->temp*p->temp);
        p->pr = (float)(pgas + prad)*lenCF*timeCF2*ivmassCF;

        /* this is still using Gamma. Perhaps could calculate gamma 
           at this point? */
        P_ratio = (double)pgas / (double)(pgas+prad);
        /* from D. Clayton's Stellar Evolution book, p.119 */
        Gammai = (double)(32. - 24.*P_ratio - 3.*P_ratio*P_ratio) / 
            (double)(24. - 18.*P_ratio - 3.*P_ratio*P_ratio);
        /*p->vsound = sqrtf_fast(Gammai * p->pr*P_ratio / p->rho_est);*/ /*code-units*/
		p->vsound = sqrtf_fast(Gamma * p->pr / p->rho_est);

        if (params.do_diffusion) {
            /* NOTE: MOST LIKELY VERY BROKEN!!!! DON'T DIFFUSE!! */

            /* Calculate temperature from u, then "create" photons (a*T^4) */
            /* Figure out good upper and lower limits for temp */
            p->u_r = A_RAD*p->temp*p->temp*p->temp*p->temp*
                lenCF*timeCF2*ivmassCF;
            p->du_r = 0.0;

            /* Calculate diffusion coefficient in user-units */
            kes = KES_COEFF/MH *(massCF * ivlenCF2);
            kff = (KFF_COEFF) * p->rho_est*pow(p->temp, -3.5)*
                (massCF * ivlenCF); 
            p->D = C_LIGHT * tdivlCF /
                ( 3.0*(kes+kff)*p->rho_est );

            /* Also, eventually, handle lightbulb approximation here */
        }
    }
}

#include "physics.h"

    void
update_point_SPHmass(SPHbody *btab, int SPHnobj, 
        void *pp, float smooth2, float newt)
{
    SPHbody *r;
    body *p = pp;
    float dr2, oneor, oneor2;
    float phii;
    Vxd(float r);
    Vxd(float ppos);

    VxV(ppos, = p->pos);

    for (r = btab; r < btab+SPHnobj; r++) {
        VxVVx(r, = r->pos, - ppos); /* 3 flops */

        dr2 = Dotx(r, r);	/* 5 flops */
        if (dr2 != (float)0.0) {	
            dr2 += smooth2;

            oneor = recipsqrtf(dr2);	/* 8 flops */

            oneor2 = oneor * oneor;	/* 17 flops */
            phii = newt * oneor * p->mass;
            r->phi -= phii;
            VVx(r->acc, -= oneor2 * phii * r);
            phii = newt * oneor * r->mass;
            p->phi -= phii;
            VVx(p->acc, += oneor2 * phii * r);
        }
    }
}


    void
update_point_SPHmass2(SPHbody *btab, int SPHnobj, float smooth2, float newt, 
        float mass)
{
    SPHbody *r;
    float dr2, oneor, oneor2;
    float phii;
    Vxd(float r);
    Vxd(float ppos);

    VxS(ppos, = (float)0.0);  /* Body fixed at origin */

    for (r = btab; r < btab+SPHnobj; r++) {
        VxVVx(r, = r->pos, - ppos); /* 3 flops */

        dr2 = Dotx(r, r);	/* 5 flops */
        if (dr2 != (float)0.0) {	
            dr2 += smooth2;

            oneor = recipsqrtf(dr2);	/* 8 flops */

            oneor2 = oneor * oneor;	/* 17 flops */
            phii = newt * oneor * mass;
            r->phi -= phii;
            VVx(r->acc, -= oneor2 * phii * r);
        }
    }
}


    void
update_point_SPHmass_bndry(SPHbody *btab, int SPHnobj, float newt, 
        bndry_t bndry)
{
    SPHbody *r;
    float dr2, oneor, oneor2;
    float phii;
    Vxd(float r);

    for (r = btab; r < btab+SPHnobj; r++) {

        VxVV(r, = r->pos, - bndry.pos); /* 3 flops */

        dr2 = Dotx(r, r);	/* 5 flops */

        oneor = recipsqrtf(dr2);	/* 8 flops */

        oneor2 = oneor * oneor;	/* 17 flops */
        phii = newt * oneor * bndry.mass;
        r->phi -= phii;
        VVx(r->acc, -= oneor2 * phii * r);
    }
}


    void 
do_SPHgrav(const float *p, const float *end, const float *pos0, float *mass0, 
        float *acc0, float *phi0, const float *eps2p, int *ncut)
{
    float dr2, h, h2, dxx, dphidx, dmassdx, v2;
    Vxd(float r);
    float phii, mor3, mass;
    VxdV(const float ppos, = pos0);
    Vxd(float a);
    float phi = *phi0;
    float total_mass = *mass0;
    const float eps2 = *eps2p; /* Nuts! Eps2 is not eps^2, but hsink!*/
    int index;

    VxV(a, = acc0);

    while (p < end) {
        mass = *p++;
        r0 = *p++;
        r1 = *p++;
        r2 = *p++;
        h = *p++;

        h=(h+eps2)/2.0;

        h2= h * h;
        VxVx(r, -= ppos);	    /* 3 flops */
        dr2 = Dotx(r, r);	    /* 5 flops */
        total_mass += mass; /* Hmm, do I need total mass or enclosed
                               here?  If enclosed, put this statement
                               after if clause*/

        if (dr2 >= 4.*h2) { /* Beyond 2h, point source for phi and acc! */
            phii = recipsqrtf(dr2); /* 8 flops */
            mor3 = phii * phii;	    /* 5 flops */
            phii *= mass;
            mor3 *= phii;
        } 
        else if (dr2 > (float) 0.) { /* Within 2h, use SPH particle
                                        smoothing for gravity */
            v2 = dr2 / h2;
            index = v2 * invdvtable;
            phii=fpoten[index]*mass/h;
            mor3=mass*fmass[index]/dr2*recipsqrtf(dr2);
        }
        phi -= phii;
        VxVx(a, += mor3 * r); /* 6 flops */
    }
    VVx(acc0, = a);
    *mass0 = total_mass;
    *phi0 = phi;
}
