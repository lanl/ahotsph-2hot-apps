#include <math.h>
#include <stdlib.h>
#include "physics_sph.h"
#include "stk.h"
#include "vop.h"
#include "singlio.h"
#include "mpmy.h"
#include "fastflpt.h"
#include "randoms.h"

/*  subroutine sets the mean molecular weight of the gas */
/*  assuming complete ionization. */
void
mmw(SPHbody *btab, int nobj)
{
  SPHbody *p;

  for (p = btab; p < btab+nobj; p++) {
    p->xmu = p->abar/(p->abar*p->ye+1.0);
  }
}

void
eosaux_setup(SPHbody *btab, int nobj)
{
  SPHbody *p;

  for (p = btab; p < btab+nobj; p++) {
      p->temprev = p->temp;
      p->rhoprev = p->rho;
      p->xpprev = p->xp;
      p->xnprev = p->xn;
      p->yeprev = p->ye;
  }
}

void
eos_prev(SPHbody *btab, int nobj)
{
  SPHbody *p;

  for (p = btab; p < btab+nobj; p++) {
      p->rhoprev = p->rho_est;
      p->yeprev = p->ye;
#if 0				/* don't update! */
      p->temprev = p->temp;
      p->xpprev = p->xp;
      p->xnprev = p->xn;
#endif
  }
}

void
eosgen_setup(SPHbody *btab, int nobj)
{
  SPHbody *p;
  double rho;
  double temp;
  double ye;
  double abar;
  double u;
  double u2;
  double pr;
  double xp;
  double xn;
  double ufreez;
  int ifleos;
  int ident;
  int procnum = MPMY_Procnum();

  for (p = btab; p < btab+nobj; p++) {
    rho = p->rho;
    temp = p->temp;
    ye = p->ye;
    abar = p->abar;
    u = p->u;
    u2 = p->u2;
    pr = p->pr;
    xp = p->xp;
    xn = p->xn;
    ufreez = p->ufreez;
    ifleos = p->ifleos;
    ident = p->ident;

    Fortran(eosgen)(&rho, &temp, &ye, &abar, &u, &u2, &pr, &xp,  &xn, 
		    &ufreez, &ifleos, &ident, &procnum);
    p->u = u;
    p->u2 = u2;
    p->pr = pr;
    p->xp = xp;
    p->xn = xn;
    p->ufreez = ufreez;
    p->ifleos = ifleos;
    /* set rhoprev to rho, since initial conditions may be significantly off */
    p->rhoprev = p->rho;
    
  }
}


void
movebound(SPHbody *btab, int nobj, float t, double rb, double *vb, int *icore)
{
  /* subroutine moving the inner boundary. */
  float sumrv, sumv, sumr, sumr2;
  float radmin, homfac;
  float ri2, ri;
  float vri, vrmin;
  int nsum;
  SPHbody *p;
  MPMY_Comm_request req;

  if (rb > 2.5e-4) {
    sumrv=0.0;
    sumv=0.0;
    sumr=0.0;
    sumr2=0.0;
    nsum=0;
    radmin=1e9;
    for (p = btab; p < btab+nobj; p++) {
      ri2=Dot(p->pos, p->pos);
      ri=sqrtf_fast(ri2);
      p->r=ri;
      if (ri < radmin) {
	vri=Dot(p->vel, p->pos)/ri;
	vrmin=vri;
	radmin=ri;
      }
      if (ri <= 5.0*rb && ri >= 2.0*rb) {
	vri=Dot(p->vel, p->pos)/ri;
	sumrv += vri*ri;
	sumv += vri;
	sumr += ri;
	sumr2 += ri2;
	nsum++;
      }
    }
    MPMY_ICombine_Init(&req);
    MPMY_ICombine(&sumrv, &sumrv, 1, MPMY_FLOAT, MPMY_SUM, req);
    MPMY_ICombine(&sumv, &sumv, 1, MPMY_FLOAT, MPMY_SUM, req);
    MPMY_ICombine(&sumr, &sumr, 1, MPMY_FLOAT, MPMY_SUM, req);
    MPMY_ICombine(&sumr2, &sumr2, 1, MPMY_FLOAT, MPMY_SUM, req);
    MPMY_ICombine(&vrmin, &vrmin, 1, MPMY_FLOAT, MPMY_MIN, req);
    MPMY_ICombine(&radmin, &radmin, 1, MPMY_FLOAT, MPMY_MIN, req);
    MPMY_ICombine(&nsum, &nsum, 1, MPMY_INT, MPMY_SUM, req);
    MPMY_ICombine_Wait(req);

    homfac=sumrv/sumr2;
    *vb=homfac*rb;
    if (*vb > 0.0) *vb = 0.0;
    singlPrintf("moveb: rb %f, vb %f, nsum %d\n", rb, *vb, nsum);
  } else {
    /* We still need to compute the radii!!! */
    for (p = btab; p < btab+nobj; p++) {
      p->r=sqrtf_fast(Dot(p->pos, p->pos));
    }
    *vb=0.0;
  }
  if (*vb >= 0.0) *icore=2;
}

void
pghost(SPHbody **btabp, int *nobjp, int *gnobjp, double rb, double vb, 
       double rbout, int iextf, int icore, float gg, float xmcore, 
       float aleph, int do_ghosts)
{
  float rb2 = rb*rb;
  float gcore;
  float vbout;
  float ri;
  float distbound, ratio;
  float vri, vrg, delta;
  float sina, cosa, sin2a, cos2a;
  SPHbody *btab = *btabp;
  SPHbody *p, *q;
  int i, rnobj;
  int n = 0;
  int idoffset;
  Stk s, *ghosts;
  int nghost;
  int nobj = *nobjp;

  /* do_ghosts == 1, do both inner and side ghosts */
  /* do_ghosts == 2, do only inner ghosts */
  /* do_ghosts == 2 is broken due to xfac problems */

  if (do_ghosts == 0) {
    return;
  }

  ghosts = &s;
  StkInitEz(ghosts);
  cosa = cos(aleph);
  sina = sin(aleph);
  cos2a = cos(2.0*aleph);
  sin2a = sin(2.0*aleph);

  if (iextf==1 && icore!=0) {
    gcore=gg*xmcore/rb2;
  } else {
    gcore=0.0;
  }
  vbout=0.0;

  rnobj = 0;
  idoffset = *gnobjp * (MPMY_Procnum()+1);
  for (p = btab; p < btab + nobj; p++) {
    ri = p->r;
    p->bghost = 0;
    if (ri-2.0*p->h < rb) {
      /* inner boundary ghosts */
      q = StkPush(ghosts, sizeof(SPHbody));
      *q = *p;
      q->ireal = p;
      distbound = ri-rb;
      ratio = rb/(ri+distbound);
      VV(q->pos, = ratio*p->pos);
      q->r = sqrtf_fast(Dot(q->pos, q->pos));
      vri = Dot(p->pos, p->vel)/ri;
      vrg = ((rb*vri+ri*vb)*(2.0*ri-rb)-(2.0*vri-vb)*ri*rb)/
	((2.0*ri-rb)*(2.0*ri-rb));
      q->vel[0] = vrg*p->pos[0]/ri+(p->vel[0]-p->pos[0]/ri*vri)*ratio;
      q->vel[1] = vrg*p->pos[1]/ri+(p->vel[1]-p->pos[1]/ri*vri)*ratio;
      q->mass = p->mass*q->pos[0]/p->pos[0];
      if (q->mass < 0.0) Error("negative mass\n");
      q->u = p->u;
      q->h = p->h*ratio;
      q->prg=2.0*distbound*gcore;
      q->rho = q->rho_est = p->rho * p->r * p->r / (q->r * q->r);
      q->pr = q->prg * q->rho + p->pr;
      q->bghost=1;
      q->ident = idoffset + n++;
    } else if (ri+2.0*p->h > rbout) {
      /* outer boundary ghosts */
      q = StkPush(ghosts, sizeof(SPHbody));
      *q = *p;
      q->ireal = p;
      distbound = rbout-ri;
      ratio = (2.0*rbout-ri)/ri;
      VV(q->pos, = ratio*p->pos);
      vri = Dot(p->pos, p->vel)/ri;
      q->vel[0] = p->vel[0]*(2.0*rbout/ri-1.0)+
	p->pos[0]*(2.0*vbout/ri-2.0*rbout*vri/(ri*ri));
      q->vel[1] = p->vel[1]*(2.0*rbout/ri-1.0)+
	p->pos[1]*(2.0*vbout/ri-2.0*rbout*vri/(ri*ri));
      q->mass = p->mass*q->pos[0]/p->pos[0];
      if (q->mass < 0.0) Error("negative mass\n");
      q->u = p->u;
      q->h = p->h*ratio;
      q->prg = 0.0;
      q->rho_est = p->rho;
      q->bghost = 1;
      q->ident = idoffset + n++;
    }
  }
  if (do_ghosts == 1) {
    /* side ghosts */
    rnobj = StkSz(ghosts)/sizeof(SPHbody);
    for (p = btab; p < btab + nobj; p++) {
      delta = sina*p->pos[0]-cosa*p->pos[1];
      if (delta <= 2.0*p->h) {
	q = StkPush(ghosts, sizeof(SPHbody));
	*q = *p;
	q->ireal = p;
	q->pos[0] =  p->pos[0]*cos2a+p->pos[1]*sin2a;
	q->pos[1] = -p->pos[0]*sin2a+p->pos[1]*cos2a;
	q->vel[0] =  p->vel[0]*cos2a+p->vel[1]*sin2a;
	q->vel[1] = -p->vel[0]*sin2a+p->vel[1]*cos2a;
	q->mass = fabs(p->mass*q->pos[0]/p->pos[0]);
	q->u = p->u;
	q->h = p->h;
	q->prg = 0.0;
	q->rho_est = p->rho;
	q->bghost = 1;
	q->ident = idoffset + n++;
      }
      delta = sina*p->pos[0]+cosa*p->pos[1];
      if (delta <= 2.0*p->h) {
	q = StkPush(ghosts, sizeof(SPHbody));
	*q = *p;
	q->ireal = p;
	q->pos[0] = p->pos[0]*cos2a-p->pos[1]*sin2a;
	q->pos[1] = p->pos[0]*sin2a+p->pos[1]*cos2a;
	q->vel[0] = p->vel[0]*cos2a-p->vel[1]*sin2a;
	q->vel[1] = p->vel[0]*sin2a+p->vel[1]*cos2a;
	q->mass = fabs(p->mass*q->pos[0]/p->pos[0]);
	q->u = p->u;
	q->h = p->h;
	q->prg = 0.0;
	q->rho_est = p->rho;
	q->bghost = 1;
	q->ident = idoffset + n++;
      }
    }
    /* side ghosts of ghosts */
    for (i = 0; i < rnobj; i++) {
      p = ((SPHbody *)StkBase(ghosts))+i;
      delta = sina*p->pos[0]-cosa*p->pos[1];
      if (delta <= 2.0*p->h) {
	q = StkPush(ghosts, sizeof(SPHbody));
	*q = *p;
	q->ireal = p;
	q->pos[0] =  p->pos[0]*cos2a+p->pos[1]*sin2a;
	q->pos[1] = -p->pos[0]*sin2a+p->pos[1]*cos2a;
	q->vel[0] =  p->vel[0]*cos2a+p->vel[1]*sin2a;
	q->vel[1] = -p->vel[0]*sin2a+p->vel[1]*cos2a;
	q->mass = fabs(p->mass*q->pos[0]/p->pos[0]);
	q->u = p->u;
	q->h = p->h;
	q->prg = 0.0;
	q->rho_est = p->rho;
	q->bghost = 1;
	q->ident = idoffset + n++;
      }
      delta = sina*p->pos[0]+cosa*p->pos[1];
      if (delta <= 2.0*p->h) {
	q = StkPush(ghosts, sizeof(SPHbody));
	*q = *p;
	q->ireal = p;
	q->pos[0] = p->pos[0]*cos2a-p->pos[1]*sin2a;
	q->pos[1] = p->pos[0]*sin2a+p->pos[1]*cos2a;
	q->vel[0] = p->vel[0]*cos2a-p->vel[1]*sin2a;
	q->vel[1] = p->vel[0]*sin2a+p->vel[1]*cos2a;
	q->mass = fabs(p->mass*q->pos[0]/p->pos[0]);
	q->u = p->u;
	q->h = p->h;
	q->prg = 0.0;
	q->rho_est = p->rho;
	q->bghost = 1;
	q->ident = idoffset + n++;
      }
    }
  }
  nghost = StkSz(ghosts)/sizeof(SPHbody);

  btab = Realloc(btab, (nobj+nghost) * sizeof(SPHbody));
  memcpy(btab+nobj, StkBase(ghosts), nghost * sizeof(SPHbody));
  *btabp = btab;
  *nobjp += nghost;
  MPMY_Combine(&nghost, &nghost, 1, MPMY_INT, MPMY_SUM);
  *gnobjp += nghost;
  StkTerminate(ghosts);
  singlPrintf("Created %d ghosts\n", nghost);
}

void
remove_ghosts(SPHbody **btabp, int *nobjp, int *gnobjp)
{
    SPHbody *btab = *btabp;
    SPHbody *p, *next, *q;
    Stk s;
    int removed = 0;

    StkInitEz(&s);
    /* Shrink btab, taking out ghost particles  */
    for (p = next = btab; p < btab+*nobjp; p++) {
	if (p->ireal == NULL) {
	    q = StkPush(&s, sizeof(SPHbody));
	    memcpy(q, p, sizeof(SPHbody));
	} else {
	  *next++ = *p;
	  removed++;
	}
    }
    StkCrunch(&s);
    memcpy(btab, StkBase(&s), StkSz(&s));
    *btabp = Realloc(btab, StkSz(&s));
    *nobjp = StkSz(&s)/sizeof(SPHbody);
    StkTerminate(&s);
    MPMY_Combine(&removed, &removed, 1, MPMY_INT, MPMY_SUM);
    *gnobjp -= removed;
    singlPrintf("Removed %d ghosts\n", removed);
}

typedef struct {
  float r;
  float mass;
  SPHbody *ptr;
  float mofr;
} sortbody;

int 
rcompare(const void *a1, const void *a2)
{
  const sortbody *b1 = a1;
  const sortbody *b2 = a2;

  if (b1->r < b2->r) return -1;
  else if (b1->r > b2->r) return 1;
  else return 0;
}

#define NBINS 10000
void
sn_gravity(SPHbody *btab, int nobj, float xmcore, float xmtheo, float gg, 
	   float clight, int icore, float rmin, float rmax)
{
  int i;
  SPHbody *p;
  float *hist;
  float gconst, pconst;
  float xmr;
  float rinv;
  float fac, r, h, massr;
  int bin;
  static int init = 0;
  static ran_state ranstate;

#define EPS (1e-4)

  if (!init) {
	ran_init(MPMY_Procnum()+1, &ranstate);
	init = 1;
  }

  hist = Calloc(NBINS, sizeof(float));
  fac = (NBINS-1)/log(rmax/rmin);
  for (i = 0; i < nobj; i++) {
    r = sqrtf_fast(Dot(btab[i].pos, btab[i].pos));
    btab[i].r = r;
    if (btab[i].bghost) continue; /* don't use ghosts ??? */
    h = btab[i].h;
    /* Fuzz the mass distributon over 0.4*h in each direction */
    if (r+0.4*h < rmax && r-0.4*h > rmin) {
        r += (uniform_rand(&ranstate)*2.0F-1.0F)*h*0.4;
    }
    if (r > rmax) Error("rmax is too small\n");
    if (r < rmin) Error("rmin is too big\n");
    bin = log(r/rmin)*fac;
    hist[bin] += btab[i].mass;
  }
  MPMY_Combine(hist, hist, NBINS, MPMY_FLOAT, MPMY_SUM);
  /* compute internal mass for all particles */
  if (icore != 0) {
    hist[0] += xmcore*xmtheo;
  } 
  for (i = 1; i < NBINS; i++) {
    hist[i] += hist[i-1];
  }

  gconst=gg/xmtheo;
  pconst=gg/(clight*clight)/xmtheo;
  for (i = 0; i < nobj; i++) {
    p = btab+i;
    if (p->bghost) continue; /* don't use ghosts */
    bin = log(p->r/rmin)*fac;
    if (bin <= 0 || bin >= NBINS) {
	Error("mass bin (%d) out of range\n", bin);
    }
    massr = hist[bin-1];	/* Shell does not feel itself */
    /* gravitational force */
    /* smooth forces near origin. Is Plummer model smoothing adequate? */
    rinv = 1.0/sqrtf_fast(EPS*EPS+p->r*p->r);
    xmr = gconst * massr * rinv * rinv * rinv;
    p->phi = -pconst * massr * rinv;
    VV(p->grav_acc, = -xmr * p->pos);
    /* gravitational redshift (w.r.t. r=infinity) */
    p->gshift = 1.0/sqrtf_fast(1.0-2.0*p->phi);
  }
  Free(hist);
}

void
swerror(void)
{
  Error("error in Fortran space\n");
}

void
swerror_(void)
{
  Error("error in Fortran space\n");
}

#include <stdio.h>

void
swflush(void)
{
  fflush(stdout);
  fflush(stderr);
}
