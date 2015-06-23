/* In wvt.c */
void SetWVT(double visc_alpha, double visc_beta, double visc_epsilon, 
	    double heat_f1, double eos_gamma, int gnobj,  
	    void bfunc(), void cfunc());
void WVTgate(SinkSPH *sink, hcell **src_vec, int *result, int n);
void macWVT(SinkSPH *sink, hcell **source, int *result, int n);
void macConstNeigh(SinkSPH *sink, hcell **source_vec, int *result, int n);
void InheritWVT(const SinkSPH *from, SinkSPH *to, hcell *pp);
void WVT_setup(int dim, int ncoef1, double *wcoef1, int ncoef2, 
	       double *wcoef2);
void SetWVTOffset(double *off, double *voff);
void UnSetWVTOffset(void);
void update_WVT(SPHbody *btab, int nobj, double dt, int *limit_high, 
		int *limit_low);

void WVTInitCube(SPHbody **btabp, int *gnobj, int *nobj, double min[NDIM], 
		 double max[NDIM], int num[NDIM], int dim);
void WVTInitHex(SPHbody **btabp, int *gnobj, int *nobj, double min[NDIM], 
		double max[NDIM], int num[NDIM], int dim);
void WVTInitCCP(SPHbody **btabp, int *gnobj, int *nobj, double min[NDIM], 
		double max[NDIM], int num[NDIM], int dim);
void WVTInitHCP(SPHbody **btabp, int *gnobj, int *nobj, double min[NDIM], 
		double max[NDIM], int num[NDIM], int dim);
void WVTupdate(SPHbody *btab, int nobj, int loop, int nloop, int dim, 
	       int nneighbors);
void WVTInitProbdist(SPHbody **btabp, int *gnobj, int *nobj, double min[NDIM], 
		     double max[NDIM], int targetnobj, double totvol, 
		     double outerbound, double innerbound, int num[NDIM], 
		     int dim);
void WVTInitProbdistlr(SPHbody **btabp, int *gnobj, int *nobj, double min[NDIM], 
		     double max[NDIM], int targetnobj, double totvol, 
		     double outerbound, double innerbound, int num[NDIM], 
		     int dim);
void interp_cylindricalgrid(SPHbody *SPHbtab, int nobj, int dimr, int dimz, 
			    int dimtheta, double minr, double maxr, 
			    double minz, double maxz, 
			    double mintheta, double maxtheta, double bgrho);
double atan2pi(double x, double y);
void init_cylindricalgrid(int dimr, int dimz, 
			  int dimtheta, double minr, double maxr, 
			  double minz, double maxz, 
			  double mintheta, double maxtheta, double center[3]);

void init_cartesiangrid(int dimx, int dimy, int dimz, double minx, 
			double maxx, double miny, double maxy, 
			double minz, double maxz, double center[3], 
			char *cartfile_rho, char *cartfile_h);
void interp_cartesiangrid(SPHbody *SPHbtab, int nobj, double bgrho);


void WVT_setinputoption(int option);
void WVT_hofpos(SPHbody *btab, int nobj, double totvol, double *tothvol, 
		int dim);
void WVT_hofpos_pwl(SPHbody *btab, int nobj, double totvol, double *tothvol,
		    int dim);
void WVT_hofpos_inputh(SPHbody *btab, int nobj, double totvol, 
		       double *tothvol, int dim);
void WVT_hofpos_cylgrid(SPHbody *btab, int nobj, double totvol, 
			double *tothvol, int dim);
void WVT_hofpos_cartgrid(SPHbody *btab, int nobj, double totvol, 
			 double *tothvol, int dim);
