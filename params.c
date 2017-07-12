#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "SDF.h"
#include "singlio.h"
#include "error.h"
#include "params.h"
#include "ndim.h"

void read_initial_ctl(SDF *sdfp, setup_params_t *params) {

	/* default values as found in main.c */
	params->setpvel = 0;

	/* read in settings */
    SDFgetintOrDefault(sdfp, "timeout", &(params->timeout), 600);
//    if (params->timeout > 0) MPMY_TimeoutSet(params->timeout);
#ifdef __PARAGON__
    {
	int fail_if_slow;
	SDFgetintOrDefault(sdfp, "fail_if_slow", &fail_if_slow, 0);
	chk_slow(fail_if_slow);
    }
#endif

    SDFgetstring(sdfp, "datafile", params->name, sizeof(params->name));
    SDFgetintOrDefault(sdfp, "do_restart", &(params->do_restart), 0);
    SDFgetintOrDefault(sdfp, "do_periodic", &(params->do_periodic), 0);
    SDFgetintOrDefault(sdfp, "cosmology", &(params->cosmology), 0);
    SDFgetintOrDefault(sdfp, "set_id", &(params->set_id), 0);
    SDFgetintOrDefault(sdfp, "setpvel", &(params->setpvel), 0);
    SDFgetintOrDefault(sdfp, "do_sph", &(params->do_sph), 1);
    SDFgetintOrDefault(sdfp, "do_diffusion", &(params->do_diffusion), 0);
    SDFgetintOrDefault(sdfp, "do_cooling", &(params->do_cooling), 0);
    SDFgetintOrDefault(sdfp, "do_burning", &(params->do_burning), 0);
    SDFgetintOrDefault(sdfp, "do_grav", &(params->do_grav), 0);
    SDFgetintOrDefault(sdfp, "do_winds", &(params->do_winds), 0);
    SDFgetintOrDefault(sdfp, "do_point_mass", &(params->do_point_mass), 0);
    SDFgetintOrDefault(sdfp, "do_point_mass2", &(params->do_point_mass2), 0);
    SDFgetintOrDefault(sdfp, "do_boundary", &(params->do_boundary), 0);
    SDFgetintOrDefault(sdfp, "do_absorbing_bndry", &(params->do_absorbing_bndry), 0);
    SDFgetintOrDefault(sdfp, "do_drag", &(params->do_drag), 0);
    SDFgetintOrDefault(sdfp, "has_grav_data", &(params->has_grav_data), params->do_grav);

	SDFgetintOrDefault (sdfp, "poly_eos", &(params->poly_eos), 0);

	if (params->do_sph || params->do_grav) {
	    if (SDFhasname("SPHdatafile", sdfp))
	        SDFgetstring(sdfp, "SPHdatafile", params->SPHdatafile, sizeof(params->SPHdatafile));
		/* read in parameters if not doing a test case or not a restart */
		if (strlen (params->SPHdatafile) > 0 || params->do_restart) {
		
			if (params->do_sph) {
				SDFgetfloatOrDefault(sdfp, "new_h", &(params->new_h), 0.0);
				SDFgetfloatOrDefault(sdfp, "new_u", &(params->new_u), 0.0);
			}
		
		    if (params->do_point_mass) 
		        SDFgetfloatOrDie(sdfp, "r_inner", &(params->r_inner));
		
		    if (params->do_point_mass2) {
		        SDFgetfloatOrDie(sdfp, "r_inner", &(params->r_inner));
		        SDFgetfloatOrDie(sdfp, "centmass", &(params->centmass));
			}
		
			if (params->do_boundary) {
			    SDFgetfloatOrDie(sdfp, "r_inner", &(params->r_inner));
			    SDFgetfloatOrDie(sdfp, "r_outer", &(params->r_outer));
			    SDFgetfloatOrDie(sdfp, "centmass", &(params->centmass));
			}
		
		    if (params->do_winds) {
		        SDFgetintOrDie(sdfp, "windpart_per_shell", &(params->windpartpershell));
		        SDFgetintOrDefault(sdfp, "old_winds", &(params->old_winds), 1);
		        SDFgetintOrDefault(sdfp, "const_winds", &(params->const_winds), 0);
		        SDFgetintOrDefault(sdfp, "nonconst_winds", &(params->nonconst_winds), 0);
		        SDFgetintOrDefault(sdfp, "accreting_winds", &(params->accreting_winds), 0);
		        if ( (params->const_winds && params->nonconst_winds) || 
		                (params->const_winds && params->accreting_winds) || 
		                (params->nonconst_winds && params->accreting_winds) )
		            Error("const_winds && nonconst_winds && accreting_winds; pick one\n");
			
			    if (params->const_winds || params->nonconst_winds || params->accreting_winds) {
			        SDFgetfloatOrDie(sdfp, "r_wind", &(params->r_wind));
			        SDFgetfloatOrDefault(sdfp, "t_wind", &(params->t_wind), 0.0);
			        if (params->const_winds || params->nonconst_winds)
			            SDFgetfloatOrDefault(sdfp, "openangle_wind", 
			                    &(params->openangle_wind), 180.0);
			        if (params->accreting_winds)
			            SDFgetfloatOrDie(sdfp, "omega_wind", &(params->omega_wind));
			        SDFgetstring(sdfp, "template_name", params->template_name,
			                sizeof(params->template_name));
			    }
			    
			    if (params->const_winds || params->accreting_winds) {
			        SDFgetfloatOrDie(sdfp, "v_wind", &(params->v_wind));
			        SDFgetfloatOrDie(sdfp, "mdot_wind", &(params->mdot_wind));
			        SDFgetfloatOrDie(sdfp, "u_wind", &(params->u_wind));
			    }
			    if (params->nonconst_winds) {
			        SDFgetstring (sdfp, "winddata_name", params->winddata_name, sizeof(params->winddata_name));
			        SDFgetfloatOrDefault (sdfp, "r_outer", &(params->r_outer), 1e30);
			    }
			}
		}
	}
    
    SDFgetfloatOrDefault(sdfp, "epsilon", &(params->eps), 0.0);
    if (params->do_grav) {
        SDFgetintOrDefault(sdfp, "do_DL", &(params->do_DL), 0);
        SDFgetintOrDefault(sdfp, "do_BH", &(params->do_BH), 0);
        SDFgetintOrDefault(sdfp, "do_Bmax", &(params->do_Bmax), 0);
        SDFgetintOrDefault(sdfp, "do_Arel", &(params->do_Arel), 0);
        if (params->do_BH || params->do_Bmax) 
            SDFgetfloatOrDie(sdfp, "theta", &(params->tol));
        else
            SDFgetfloatOrDie(sdfp, "errtol", &(params->tol));
        SDFgetfloatOrDefault(sdfp, "frac_tol", &(params->frac_tol), 0.0);
    }

    SDFgetfloatOrDefault(sdfp, "CWfac", &(params->CWfac), 0.0);
    SDFgetfloatOrDefault(sdfp, "SPHCWfac", &(params->SPHCWfac), 0.0);

	/* get time step from ctl file on new starts */
    if (!params->do_restart) 
        SDFgetfloatOrDie(sdfp, "dt", &(params->dt));

    SDFgetfloatOrDefault(sdfp, "dark_dt", &(params->dark_dt), (params->do_grav ? params->dt : 1e30));
    SDFgetintOrDie(sdfp, "nsteps", &(params->nsteps));
    SDFgetintOrDefault(sdfp, "log_time", &(params->log_time), 0);
    SDFgetintOrDefault(sdfp, "comov_eps", &(params->comov_eps), 0);
    SDFgetfloatOrDefault(sdfp, "comov_eps_epoch", &(params->comov_eps_epoch), 10.0);
    SDFgetintOrDefault(sdfp, "save_first", &(params->save_first), 0);
    SDFgetintOrDefault(sdfp, "ntimer_detail", &(params->ntimer_detail), 0);
    SDFgetintOrDefault(sdfp, "exact_rho", &(params->exact_rho), 0);
    SDFgetfloatOrDefault(sdfp, "visc_alpha", &(params->visc_alpha), (float)1.0);
    SDFgetfloatOrDefault(sdfp, "visc_beta", &(params->visc_beta), (float)2.0);
    SDFgetfloatOrDefault(sdfp, "visc_epsilon", &(params->visc_epsilon), (float)1e-2);
    SDFgetfloatOrDefault(sdfp, "heat_f1", &(params->heat_f1), (float)0.0);
    SDFgetfloatOrDie(sdfp, "gamma", &(params->Gamma));
    SDFgetfloatOrDefault(sdfp, "courant_number", &(params->courant_number), (float)0.4);
    SDFgetfloatOrDefault(sdfp, "min_h", &(params->min_h), (float)0.0);
    SDFgetfloatOrDefault(sdfp, "max_h", &(params->max_h), (float)1e30);
    SDFgetintOrDefault(sdfp, "nbrcut_max", &(params->nbrcut_max), 500);
    SDFgetintOrDefault(sdfp, "nbrcut_min", &(params->nbrcut_min), 10);
    SDFgetfloatOrDefault(sdfp, "nbrcut_fac", &(params->nbrcut_fac), (float)0.1);
    SDFgetintOrDefault(sdfp, "adaptive_dt", &(params->adaptive_dt), 1);
    SDFgetintOrDefault(sdfp, "independent_dt", &(params->independent_dt), 0);
    SDFgetintOrDefault(sdfp, "dark_independent_dt", &(params->dark_independent_dt), 0);
    SDFgetintOrDefault(sdfp, "default_nterms", &(params->default_nterms), 100);

    SDFgetfloatOrDefault(sdfp, "massCF", &(params->fmassCF), 1.0);/*mass conversion factor; CE*/
    SDFgetfloatOrDefault(sdfp, "lengthCF", &(params->flenCF), 1.0);/*length conversion factor; CE*/
    SDFgetfloatOrDefault(sdfp, "timeCF", &(params->ftimeCF), 1.0);/*time conversion factor; CE*/

    if (params->adaptive_dt) {
        SDFgetintOrDefault(sdfp, "tlow_cut", &(params->tlow_cut), 40);
        SDFgetintOrDefault(sdfp, "dt_short", &(params->dt_short), 0);
        SDFgetintOrDefault(sdfp, "dt_long", &(params->dt_long), 10);
        SDFgetfloatOrDefault(sdfp, "dt_max", &(params->dt_max), 1e30);
    }

	SDFgetstringOrDefault(sdfp, "outfile", params->outnamebase, sizeof (params->outnamebase), "");
	/*
    if (SDFhasname ("outfile", sdfp)) 
        SDFgetstring (sdfp, "outfile", params->outnamebase, sizeof (params->outnamebase));
    else
        sprintf(params->outnamebase,"%s", "");
		*/
	if (strlen (params->outnamebase) > 0) {
		params->do_output = 1;
	}

	if (params->do_output) {
        SDFgetintOrDefault (sdfp, "output_freq", &(params->output_freq), params->nsteps);
        SDFgetintOrDefault (sdfp, "short_output", &(params->short_output), 0);
    } else {
        params->output_freq = 1;
	}

    SDFgetintOrDefault(sdfp, "timer_freq", &(params->timer_freq), params->output_freq);
    SDFgetfloatOrDefault(sdfp, "sort_tol", &(params->sort_tol), 0.01);
    SDFgetintOrDefault(sdfp, "image_freq", &(params->image_freq), 0);
    SDFgetintOrDefault(sdfp, "x_pixels", &(params->x_pixels), 512);
    SDFgetintOrDefault(sdfp, "y_pixels", &(params->y_pixels), 512);
    SDFgetintOrDefault(sdfp, "log_image", &(params->log_image), 0);

    /* read in the kernel coefficients */
    if (SDFhasname("kernel_ncoef1", sdfp)) {
        SDFgetintOrDie(sdfp, "kernel_ncoef1", &(params->kernel_ncoef1));
        if (params->kernel_ncoef1 >= MAXCOEF) Error("Increase MAXCOEF\n");
        SDFgetintOrDie(sdfp, "kernel_ncoef2", &(params->kernel_ncoef2));
        if (params->kernel_ncoef2 >= MAXCOEF) Error("Increase MAXCOEF\n");
        if (SDFseekrdvecs(sdfp, "kernel_coef1", 0, params->kernel_ncoef1, 
                    &(params->kernel_coef1), 0, NULL))
            Error("SDFread kernel_coef1 failed\n");
        if (SDFseekrdvecs(sdfp, "kernel_coef2", 0, params->kernel_ncoef2, 
                    &(params->kernel_coef2), 0, NULL))
            Error("SDFread kernel_coef2 failed\n");
    } else {
        /* Monaghan spline kernel is default */
        params->kernel_ncoef1 = params->kernel_ncoef2 = 4;
        params->kernel_coef1[0] = 1.0;		params->kernel_coef2[0] = 2.0;
        params->kernel_coef1[1] = 0.0;		params->kernel_coef2[1] = -3.0;
        params->kernel_coef1[2] = -3.0/2.0;	params->kernel_coef2[2] = 3.0/2.0;
        params->kernel_coef1[3] = 3.0/4.0;	params->kernel_coef2[3] = -1.0/4.0;
    }

    if (params->do_drag) {
        SDFgetfloatOrDie(sdfp, "drag_coeff", &(params->drag_coeff));
    }
}

void print_initial_ctl(setup_params_t params) {
    singlPrintf("printing params structure\n");

	singlPrintf("float errtol = %g;\n", params.tol);
    singlPrintf("float dark_dt = %g;\n", params.dark_dt);
	singlPrintf("float eps = %g;\n", params.eps);
    singlPrintf("int nsteps = %d;\n", params.nsteps);
	singlPrintf("int do_Bmax = %d;\n", params.do_Bmax);
	singlPrintf("int do_BH = %d;\n", params.do_BH);
	singlPrintf("int do_Arel = %d;\n", params.do_Arel);
	singlPrintf("int do_DL = %d;\n", params.do_DL);
    singlPrintf("int exact_rho = %d;\n", params.exact_rho);
    singlPrintf("float courant_number = %g;\n", params.courant_number);
    singlPrintf("float gamma = %f;\n", params.Gamma);
    singlPrintf("float massCF = %g;\n", params.fmassCF);/*added by CE*/
    singlPrintf("float lenCF = %g;\n", params.flenCF);/*added by CE*/
    singlPrintf("float timeCF = %g;\n", params.ftimeCF);/*added by CE*/
    singlPrintf("float visc_alpha = %g;\n", params.visc_alpha);
    singlPrintf("float visc_beta = %g;\n", params.visc_beta);
    singlPrintf("float visc_epsilon = %g;\n", params.visc_epsilon);
    singlPrintf("float heat_f1 = %g;\n", params.heat_f1);
    singlPrintf("float min_h = %g;\n", params.min_h);
    singlPrintf("float max_h = %g;\n", params.max_h);
    singlPrintf("int adaptive_dt = %d;\n", params.adaptive_dt);
    singlPrintf("int independent_dt = %d;\n", params.independent_dt);
    singlPrintf("int dark_independent_dt = %d;\n", params.dark_independent_dt);
	if (params.adaptive_dt) {
		singlPrintf("int tlow_cut = %d;\n", params.tlow_cut);
		singlPrintf("int dt_short = %d;\n", params.dt_short);
		singlPrintf("int dt_long = %d;\n", params.dt_long);
		singlPrintf("float dt_max = %f;\n", params.dt_max);
	}
    if (params.do_winds) {
        singlPrintf("int do_winds = %d;\n", params.do_winds);
        singlPrintf("int windpartpershell = %d;\n", params.windpartpershell);
        singlPrintf("int old_winds = %d;\n", params.old_winds);
        if (params.old_winds) {
            //singlPrintf("int windgnobj = %d;\n", windgnobj);
        }
        singlPrintf("int const_winds = %d;\n", params.const_winds);
        singlPrintf("int nonconst_winds = %d;\n", params.nonconst_winds);
        singlPrintf("int accreting_winds = %d;\n", params.accreting_winds);
        if (params.const_winds || params.accreting_winds) {
            singlPrintf("float v_wind = %g;\n", params.v_wind);
            singlPrintf("float mdot_wind = %g;\n", params.mdot_wind);
            singlPrintf("float u_wind = %g;\n", params.u_wind);
        }
        if (params.const_winds || params.nonconst_winds || params.accreting_winds) {
            singlPrintf("float r_wind = %g;\n", params.r_wind);
            singlPrintf("float t_wind = %g;\n", params.t_wind);
            if (params.const_winds || params.nonconst_winds)
                singlPrintf("float openangle_wind = %g;\n", params.openangle_wind);
            singlPrintf("char template_name[] = \"%s\"\n", params.template_name);
            srand48(192837465);
        }
        if (params.accreting_winds)
            singlPrintf("float omega_wind = %g;\n", params.omega_wind);
        if (params.nonconst_winds) {
            singlPrintf("char winddata_name[] = \"%s\"\n", params.winddata_name);
            singlPrintf("float r_outer = %g;\n", params.r_outer);
        }
    }
    if (params.do_point_mass || params.do_point_mass2) {
        singlPrintf("float r_inner = %f;\n", params.r_inner);
        singlPrintf("float centmass = %e;\n", params.centmass);
    }
    if (params.do_boundary) {
        singlPrintf("float r_inner = %f;\n", params.r_inner);
        singlPrintf("float r_outer = %f;\n", params.r_outer);
        singlPrintf("float centmass = %e;\n", params.centmass);
    }
	if (params.do_drag) {
		singlPrintf("int do_drag = %d;\n", params.do_drag);
		singlPrintf("float drag_coeff = %g;\n", params.drag_coeff);
	}
    singlPrintf("int do_cooling = %d;\n", params.do_cooling);
    singlPrintf("int do_diffusion = %d;\n", params.do_diffusion);
    singlPrintf("int do_burning = %d;\n", params.do_burning);

	if (params.do_output) {
		if (params.short_output) 
			singlPrintf("Output to %s.nnnn, every %d steps\n", 
					params.outnamebase, params.output_freq);
	} else {
		singlPrintf("No output.\n");
	}
	singlPrintf("int timer_freq = %d;\n", params.timer_freq);
	singlPrintf("float sort_tol = %g;", params.sort_tol);
    singlPrintf("int do_periodic = %d;\n", params.do_periodic);

    singlPrintf("kernel coefficients:\n\t");
    for (int i = 0; i < params.kernel_ncoef1; i++)
        singlPrintf("%12.9f ", params.kernel_coef1[i]);
    singlPrintf("\n\t");
    for (int i = 0; i < params.kernel_ncoef2; i++)
        singlPrintf("%12.9f ", params.kernel_coef2[i]);
    singlPrintf("\n");
    if (params.log_time) Error("This code does not support log_time\n");
	if (params.cosmology) {
		singlPrintf("int cosmology = %d;\n", params.cosmology);
		singlPrintf("int comov_eps = %d;\n", params.comov_eps);
		singlPrintf("float comov_eps_epoch = %f;\n", params.comov_eps_epoch);
		singlPrintf("int setpvel = %d;\n", params.setpvel);
	}

	/* these are not printed
	singlPrintf("int timeout = %d;\n", params.timeout);
	singlPrintf("string datafile[] = %s;\n", params.name);
	singlPrintf("int do_restart = %d;\n", params.do_restart);
	singlPrintf("int set_id = %d;\n", params.set_id);
	singlPrintf("int do_sph = %d;\n", params.do_sph);
	singlPrintf("int do_grav = %d;\n", params.do_grav);
    singlPrintf("int do_winds = %d;\n", params.do_winds);
	singlPrintf("int do_point_mass = %d;\n", params.do_point_mass);
	singlPrintf("int do_point_mass2 = %d;\n", params.do_point_mass2);
	singlPrintf("int do_boundary = %d;\n", params.do_boundary);
	singlPrintf("int do_absorbing_bndry = %d;\n", params.do_absorbing_bndry);
	singlPrintf("int has_grav_data = %d;\n", params.has_grav_data);

	singlPrintf("string SPHdatafile[] = %s;\n", params.SPHdatafile);

	singlPrintf("float CWfac = %g;\n", params.CWfac);
	singlPrintf("float SPHCWfac = %g;\n", params.SPHCWfac);

    singlPrintf("float epsilon = %g;\n", params.eps);

	singlPrintf("int log_time = %d;\n", params.log_time);
    singlPrintf("int save_first = %d;\n", params.save_first);
    singlPrintf("int ntimer_detail = %d;\n", params.ntimer_detail);
    singlPrintf("int nbrcut_max = %d;\n", params.nbrcut_max);
    singlPrintf("int nbrcut_min = %d;\n", params.nbrcut_min);
    singlPrintf("float nbrcut_fac = %g;\n", params.nbrcut_fac);
    singlPrintf("int default_nterms = %d;\n", params.default_nterms);

	singlPrintf("string outnamebase[] = %s;\n", params.outnamebase);
	singlPrintf("int do_output = %d;\n", params.do_output);
	singlPrintf("int output_freq = %d;\n", params.output_freq);
	singlPrintf("int short_output = %d;\n", params.short_output);
	singlPrintf("int image_freq = %d;\n", params.image_freq);
	singlPrintf("int x_pixels = %d;\n", params.x_pixels);
	singlPrintf("int y_pixels = %d;\n", params.y_pixels);
	singlPrintf("int log_image = %d;\n", params.log_image);
	*/

    singlPrintf("end printing params structure\n");
}

void print_absorb_bndry(bndry_t bndry) {
    if (params.do_absorbing_bndry) {
		/* print position */
		singlPrintf("float bndry[] = [ %g", bndry.pos[0]);
#if NDIM>=2
		singlPrintf(", %g", bndry.pos[1]);
#if NDIM>=3
		singlPrintf(", %g", bndry.pos[2]);
#endif
#endif
		singlPrintf(" ];\n");

		/* print velocity */
		singlPrintf("float bndry_vel[] = [ %g", bndry.vel[0]);
#if NDIM>=2
		singlPrintf(", %g", bndry.vel[1]);
#if NDIM>=3
		singlPrintf(", %g", bndry.vel[2]);
#endif
#endif
		singlPrintf(" ];\n");

		/* print linear momentum */
		singlPrintf("float bndry_p[] = [ %g", bndry.p[0]);
#if NDIM>=2
		singlPrintf(", %g", bndry.p[1]);
#if NDIM>=3
		singlPrintf(", %g", bndry.p[2]);
#endif
#endif
		singlPrintf(" ];\n");

		/* print angular momentum */
		singlPrintf("float bndry_l[] = [ %g", bndry.l[0]);
#if NDIM>=2
		singlPrintf(", %g", bndry.l[1]);
#if NDIM>=3
		singlPrintf(", %g", bndry.l[2]);
#endif
#endif
		singlPrintf(" ];\n");

        singlPrintf("float bndry_mass = %g;\n", bndry.mass);
        singlPrintf("float bndry_r = %g;\n", bndry.r);
    }
}

void read_absorb_bndry(SDF *sdfp, bndry_t *bndry) {
    SDFgetfloatOrDie(sdfp, "bndry_x", &(bndry->pos[0]));
    SDFgetfloatOrDie(sdfp, "bndry_vx", &(bndry->vel[0]));
    SDFgetfloatOrDefault(sdfp, "bndry_lx", &(bndry->l[0]), 0.0);
    SDFgetfloatOrDefault(sdfp, "bndry_px", &(bndry->p[0]), 0.0);
#if NDIM>=2        
    SDFgetfloatOrDie(sdfp, "bndry_y", &(bndry->pos[1]));
    SDFgetfloatOrDie(sdfp, "bndry_vy", &(bndry->vel[1]));
    SDFgetfloatOrDefault(sdfp, "bndry_ly", &(bndry->l[1]), 0.0);
    SDFgetfloatOrDefault(sdfp, "bndry_py", &(bndry->p[1]), 0.0);
#if NDIM>=3
    SDFgetfloatOrDie(sdfp, "bndry_z", &(bndry->pos[2]));
    SDFgetfloatOrDie(sdfp, "bndry_vz", &(bndry->vel[2]));
    SDFgetfloatOrDefault(sdfp, "bndry_lz", &(bndry->l[2]), 0.0);
    SDFgetfloatOrDefault(sdfp, "bndry_pz", &(bndry->p[2]), 0.0);
#endif
#endif
    SDFgetfloatOrDie(sdfp, "bndry_mass", &(bndry->mass));
    SDFgetfloatOrDie(sdfp, "bndry_r", &(bndry->r));

}
