#include <math.h>
#include <stdlib.h>
#include "physics_sph.h"
#include "stk.h"
#include "vop.h"
#include "singlio.h"
#include "mpmy.h"

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
      p->temprev = p->temp;
      p->rhoprev = p->rho;
      p->xpprev = p->xp;
      p->xnprev = p->xn;
      p->yeprev = p->ye;
  }
}

void
movebound(SPHbody *btab, int nobj, float t, float rb, float *vb, int *icore)
{
  /* subroutine moving the inner boundary. */
  float sumrv, sumv, sumr, sumr2;
  float radmin, homfac;
  float ri2, ri;
  float vri, vrmin;
  int nsum;
  SPHbody *p;

  if (rb > 2.5e-4) {
    if (*icore==1) {
      sumrv=0.0;
      sumv=0.0;
      sumr=0.0;
      sumr2=0.0;
      nsum=0;
      radmin=1e9;
      for (p = btab; p < btab+nobj; p++) {
	ri2=Dot(p->pos, p->pos);
	ri=sqrt(ri2);
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
      homfac=sumrv/sumr2;
      *vb=homfac*rb;
      singlPrintf("moveb: rb %f, vb %f, nsum %d\n", rb, *vb, nsum);
    } else {
      /* We still need to compute the radii!!! */
      for (p = btab; p < btab+nobj; p++) {
	p->r=sqrt(Dot(p->pos, p->pos));
      }
      *vb=0.0;
    }
    if (*vb >= 0.0) *icore=2;
  }
}

void
pghost(SPHbody *btab, int nobj, int *nghost, Stk *ghosts, 
       float rb, float vb, float rbout, int iextf, int icore, float gg, 
       float xmcore, float aleph)
{
  float rb2 = rb*rb;
  float gcore;
  float vbout;
  float ri;
  float distbound, ratio;
  float vri, vrg, delta;
  float sina, cosa, sin2a, cos2a;
  SPHbody *p, *q;
  int i, rnobj;

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
      q->bghost=1;
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
      q->bghost = 1;
    }
  }
  rnobj = StkSz(ghosts)/sizeof(SPHbody);
  /* side ghosts */
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
      q->mass = p->mass*q->pos[0]/p->pos[0];
      if (q->mass < 0.0) Error("negative mass\n");
      q->u = p->u;
      q->h = p->h;
      q->prg = 0.0;
      q->bghost = p->bghost;
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
      q->mass = p->mass*q->pos[0]/p->pos[0];
      if (q->mass < 0.0) Error("negative mass\n");
      q->u = p->u;
      q->h = p->h;
      q->prg = 0.0;
      q->bghost = p->bghost;
    }
  }
  /* side ghosts of ghosts */
  for (i = 0; i < rnobj; i++) {
    delta = sina*p->pos[0]-cosa*p->pos[1];
    if (delta <= 2.0*p->h) {
      q = StkPush(ghosts, sizeof(SPHbody));
      p = ((SPHbody *)StkBase(ghosts))+i;
      *q = *p;
      q->ireal = p;
      q->pos[0] =  p->pos[0]*cos2a+p->pos[1]*sin2a;
      q->pos[1] = -p->pos[0]*sin2a+p->pos[1]*cos2a;
      q->vel[0] =  p->vel[0]*cos2a+p->vel[1]*sin2a;
      q->vel[1] = -p->vel[0]*sin2a+p->vel[1]*cos2a;
      q->mass = fabs(p->mass*q->pos[0]/p->pos[0]);
      if (q->mass < 0.0) Error("negative mass\n");
      q->u = p->u;
      q->h = p->h;
      q->prg = 0.0;
      q->bghost = p->bghost;
    }
    delta = sina*p->pos[0]+cosa*p->pos[1];
    if (delta <= 2.0*p->h) {
      q = StkPush(ghosts, sizeof(SPHbody));
      p = ((SPHbody *)StkBase(ghosts))+i;
      *q = *p;
      q->ireal = p;
      q->pos[0] = p->pos[0]*cos2a-p->pos[1]*sin2a;
      q->pos[1] = p->pos[0]*sin2a+p->pos[1]*cos2a;
      q->vel[0] = p->vel[0]*cos2a-p->vel[1]*sin2a;
      q->vel[1] = p->vel[0]*sin2a+p->vel[1]*cos2a;
      q->mass = fabs(p->mass*q->pos[0]/p->pos[0]);
      if (q->mass < 0.0) Error("negative mass\n");
      q->u = p->u;
      q->h = p->h;
      q->prg = 0.0;
      q->bghost = p->bghost;
    }
  }
  *nghost = StkSz(ghosts)/sizeof(SPHbody);
}

void
remove_ghosts(SPHbody **btabp, int *nobjp)
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
    singlPrintf("Removed %d ghost particles\n", removed);
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

#define NBINS 100000
void
sn_gravity(SPHbody *btab, int nobj, float xmcore, float xmtheo, float gg, 
	   float clight, int icore, float rmin, float rmax)
{
  int i;
  SPHbody *p;
  float *hist;
  float gconst, pconst;
  float xmr;
  float fac, r, massr;
  int bin;

  hist = Calloc(NBINS+1, sizeof(float));
  fac = NBINS/log(rmax/rmin);
  for (i = 0; i < nobj; i++) {
    r = btab[i].r;
    if (r > rmax) Error("rmax is too small\n");
    if (r < rmin) Error("rmin is too big\n");
    bin = log(btab[i].r/rmin)*fac;
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
    bin = log(btab[i].r/rmin)*fac;
    massr = hist[bin];		/* Shell feels itself */
    p = btab+i;
    /* gravitational force */
#if 0
    if (p->r < 2.0 * p->h) {
      /* smooth forces near origin. Is Plummer model smoothing adequate? */
      xmr = gconst * massr/(p->r * p->r * p->r + p->h * p->h);
      p->phi = -pconst * massr / (p->r + p->h);
    } else {
      xmr = gconst * massr/(p->r * p->r * p->r);
      p->phi = -pconst * massr / p->r;
    }
#else
#define EPS (5e-4)
    /* smooth forces near origin. Is Plummer model smoothing adequate? */
    xmr = gconst * massr/(EPS*EPS*EPS+p->r*p->r*p->r);
    p->phi = -pconst * massr / (EPS+p->r);
#endif
    VV(p->grav_acc, = -xmr * p->pos);
    /* gravitational redshift (w.r.t. r=infinity) */
    p->gshift = 1.0/sqrt(1.0-2.0*p->phi);
  }
  Free(hist);
}




