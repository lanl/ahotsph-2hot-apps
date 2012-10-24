void output(const char *outnamebase, int64_t gnobj, int nobj, const body *btab, int iter, 
	    double dt, double dtv, 
	    struct cosmo_s *cosmo, double tpos, double tvel, 
	    int cosmology, int do_periodic, 
	    float eps, float this_eps_scaled, int force_smoothing_type,
	    float this_tol, float frac_tol, float frac_tol0, 
	    const float *R, const int *N, int write_nfiles, double *ke, double *pe, 
	    int do_output, int identsort_output);
