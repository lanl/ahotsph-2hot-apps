
static void sfzero(int nu, int nv, complex *f, int nd)
{
    /* zero signature functions */
    int q;

    if (nd != 1) Error("Assumes nd = 1\n");

    for (q = 0; q < nu*nv; q++) {
	f[q].r = 0.0;
	f[q].i = 0.0;
    }
}

static void sfadd(int nu, int nv, complex *f, complex *g, complex *h, int nd)
{
    /* add signature functions */
    int q;

    if (nd != 1) Error("Assumes nd = 1\n");

    for (q = 0; q < nu*nv; q++) {
	h[q].r = f[q].r + g[q].r;
	h[q].i = f[q].i + g[q].i;
    }
}

static void sfmult(int nu, int nv, complex *f, complex *g, complex *h, int nd)
{
    /* multiply signature functions */
    /* treat g as translation operator */
    int q;
    double r, i;	  /* source and dest might be identical */

    for (q = 0; q < nu*nv; q++) {
	r = f[q].r * g[q].r - f[q].i * g[q].i;
	i = f[q].r * g[q].i + f[q].i * g[q].r;
	h[q].r = r;
	h[q].i = i;
    }
}


/*     generate theta and phi grids for the signature function tabulation */
/*     The theta grid is an equispaced 1/2 shifted grid on [0,pi]         */
/*     The phi grid is equispaced with nph intervals & points on [0,2*pi) */
static void gnthph(int nth, int nph, double *csth, double *csph)
{
    int ith, iph;
    double arg, dth, dph;

    dth = M_PI/nth;
    dph = (2.*M_PI)/nph;

    if (2*nth > N_u0 || 2*nph > N_v0) Error("Bad args\n");

    for (ith = 0; ith < nth; ith++) {
	arg = (ith+.5)*dth;
	csth[2*ith] = cos(arg);
	csth[2*ith+1] = sin(arg);
    }

    for (iph = 0; iph < nph; iph++) {
	arg = iph*dph;
	csph[2*iph] = cos(arg);
	csph[2*iph+1] = sin(arg);
    }
}

/* 
   generate the diagonal form for the outer to outer translation
   or
   generate the diagonal form for the inner to inner translation
  
   d(j,i) = exp(-ik*s(i,j).r)  where r=b-a with a=(old center of
                               expansion) and b=(new center of
                               expansion)
   (i,j) represent grid points on the unit sphere:
         (theta_i, phi_j)
         -->this is swapped for d(j,i)
 */
static void gendto(int nx, int mx, complex *d, double *csu, double *csv, 
		   float *r, double k)
{
    int i, j;
    double rdotcs, arg;

    /* loop over all points on the sphere */
    for (j = 0; j < mx; j++) {
	rdotcs = csv[2*j]*r[0] + csv[2*j+1]*r[1];
	for (i = 0; i < nx; i++) {
	    arg = -k*(csu[2*i+1]*rdotcs + csu[2*i]*r[2]);
	    /* This loop is inside out, and should be fixed */
	    d[i*mx+j].r = cos(arg);
	    d[i*mx+j].i = sin(arg);
	}
    }
}

static void genffsf(double wt, int nu, int nv, complex *ffsf, complex *c, 
		    int nd)
{
    /* generate signature functions for scalar Helmholtz */
    int q;

    if (nd != 1) Error("Assumes nd = 1\n");
    
    for (q = 0; q < nu*nv; q++) {
	ffsf[q].r = wt*c[q].r;
	ffsf[q].i = wt*c[q].i;
    }
}



	    sortresult_t fakesortedbtab;
	    tree_t faketree;
	    body faketab[2];
	    body *fakebtab = faketab;


	{
	    faketab[0].strength = 1.;
	    faketab[0].pos[0] = .105;
	    faketab[0].pos[1] = -.025;
	    faketab[0].pos[2] = 1.015;
	    faketab[0].phi_i = faketab[0].phi_r = 0;
	    faketab[1].strength = 1.;
	    faketab[1].pos[0] = .235;
	    faketab[1].pos[1] = -.025;
	    faketab[1].pos[2] = 1.015;
	    faketab[1].phi_i = faketab[1].phi_r = 0;
	    FixKeys(faketab, 2, GetKey);
	    pqsortsetup(&fakesortedbtab, fakebtab, 2, sizeof(body), 
			0.01, Realloc_f);
	    SetupTree(&faketree, NDIM, 
		      sizeof(body), sizeof(cell), TBODYSZ, sizeof(cofmdata),
		      (pq_keyproto)GetKeyFromStruct, (pq_wgtproto)GetCost,
		      CofmFromDaugh, (cellfromcofm_t)CellFromCofm);
	    BuildTree(&faketree, &fakesortedbtab);
	    fakebtab = fakesortedbtab.data;
	    Walk(&thetree, &faketree, sizeof(Sink), mac, inherit);
	}
