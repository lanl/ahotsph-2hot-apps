#include <math.h>
#include "mpmy.h"
#include "physics_vrtx.h"
#include "vop.h"
#include "fastflpt.h"
#include "singlio.h"
#include "protos.h"

#define FPI 1.256637062e+01F

extern double omega_tot[3], lin_impulse[3], ang_impulse[3], ke, en, he;
Counter_t NtermsCnt;

void GlobalDiags(body *btab, int nobj){
    double pcs[3];
    double esum, esum2;
    double ntermslocal;
    float maxerr, maxerr2;
    int gnobj;
    body *p;
    MPMY_Comm_request req;

    /* Evaluate total vorticity, linear impulse, angular impulse,
       kinetic energy, enstrophy and helicity  */
	    
    VS(omega_tot, = (float)0.0);
    VS(lin_impulse, = (float)0.0);
    VS(ang_impulse, = (float)0.0);
    ke = (float)0.0;
    en = (float)0.0;
    he = (float)0.0;
    esum = (float)0.0;
    esum2 = (float)0.0;
    maxerr = (float)0.0;
    maxerr2 = (float)0.0;
    ntermslocal = 0.0;

    for (p = btab; p < btab+nobj; p++){
	float str2;

        VV(omega_tot, += Strength(p));
        pcs[0] = Pos(p)[1]*Strength(p)[2]-Pos(p)[2]*Strength(p)[1]; 
        pcs[1] = Pos(p)[2]*Strength(p)[0]-Pos(p)[0]*Strength(p)[2]; 
        pcs[2] = Pos(p)[0]*Strength(p)[1]-Pos(p)[1]*Strength(p)[0]; 
        VV(lin_impulse, += 0.5F * pcs); 
        ang_impulse[0] += (Pos(p)[1]*pcs[2]-Pos(p)[2]*pcs[1])/3.F;
        ang_impulse[1] += (Pos(p)[2]*pcs[0]-Pos(p)[0]*pcs[2])/3.F;
        ang_impulse[2] += (Pos(p)[0]*pcs[1]-Pos(p)[1]*pcs[0])/3.F;

	ke += 0.5 * Dot(Strength(p), Psi(p));
	en += Strength(p)[0]*(Gradvel(p)[2][1]-Gradvel(p)[1][2])
	    +Strength(p)[1]*(Gradvel(p)[0][2]-Gradvel(p)[2][0])
		+Strength(p)[2]*(Gradvel(p)[1][0]-Gradvel(p)[0][1]);
	he +=              Dot(Strength(p), Vel(p));
	esum += Errsum(p);
	esum2 += sqrtf_fast( Errsum2(p) );
	str2 = Dot(Strength(p), Strength(p));
	if( !finite(str2) || str2 > 1.e20 ){
	    SeriousWarning("Particle %d looks bad.\n", p - btab);
	    Shout("id=%d, Pos=%g %g %g, Str=%g %g %g\n",
		  p->ident, 
		  Pos(p)[0], Pos(p)[1], Pos(p)[2],
		  Strength(p)[0], Strength(p)[1], Strength(p)[2]);
	    Shout("dstr = %g %g %g\n", 
		  (p->dstr)[0], (p->dstr)[1], (p->dstr)[2]);
	}
		  
	if( maxerr < Errsum(p) ) maxerr = Errsum(p);
	if( maxerr2 < Errsum2(p) ) maxerr2 = Errsum2(p);
	ntermslocal += p->nterms;
    }
    AddCounter(&NtermsCnt, (int)ntermslocal);  /* might overflow? */

    gnobj = nobj;
    MPMY_ICombine_Init(&req);
    MPMY_ICombine(&omega_tot, &omega_tot, 3, MPMY_DOUBLE, MPMY_SUM, req);
    MPMY_ICombine(&lin_impulse, &lin_impulse, 3, MPMY_DOUBLE, MPMY_SUM, req);
    MPMY_ICombine(&ang_impulse, &ang_impulse, 3, MPMY_DOUBLE, MPMY_SUM, req);
    MPMY_ICombine(&ke, &ke, 1, MPMY_DOUBLE, MPMY_SUM, req);
    MPMY_ICombine(&en, &en, 1, MPMY_DOUBLE, MPMY_SUM, req);
    MPMY_ICombine(&he, &he, 1, MPMY_DOUBLE, MPMY_SUM, req);
    MPMY_ICombine(&esum, &esum, 1, MPMY_DOUBLE, MPMY_SUM, req);
    MPMY_ICombine(&esum2, &esum2, 1, MPMY_DOUBLE, MPMY_SUM, req);
    MPMY_ICombine(&maxerr, &maxerr, 1, MPMY_FLOAT, MPMY_MAX, req);
    MPMY_ICombine(&maxerr2, &maxerr2, 1, MPMY_FLOAT, MPMY_MAX, req);
    MPMY_ICombine(&gnobj, &gnobj, 1, MPMY_INT, MPMY_SUM, req);
    MPMY_ICombine_Wait(req);

    singlPrintf("omega_tot= %g %g %g\n", FPI*omega_tot[0],
                            FPI*omega_tot[1], FPI*omega_tot[2]); 
    singlPrintf("lin_impulse= %g %g %g\n", FPI*lin_impulse[0],
                            FPI*lin_impulse[1], FPI*lin_impulse[2]);
    singlPrintf("ang_impulse= %g %g %g\n", FPI*ang_impulse[0],
                            FPI*ang_impulse[1], FPI*ang_impulse[2]);

    singlPrintf("ke = %g, en = %g, he = %g\n", FPI*ke, FPI*en, FPI*he); 
    esum /= gnobj;
    esum2 /=  gnobj;
    singlPrintf("sum error bounds: mean: %g, max: %g\n", esum, maxerr);
    singlPrintf("sqrt(sum sq err bounds): mean: %g, max: %g\n",
		esum2, sqrt(maxerr2));
    singlPrintf("Errs: %g %g %g %g\n", esum, maxerr, esum2, sqrt(maxerr2));
    singlFflush();
}
	    
