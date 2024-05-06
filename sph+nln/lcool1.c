/************************************************************************
 * PURPOSE:								*
 *   to calculate the cooling term for a given temperature and 		*
 *   composition, and return the cooling value to the calling routine.	*
 *   Returns total cooling term in erg*cm^3/s.				*
 *   The composition is taken directly from SPHbody->abund[i] and	* 
 *   converted from mass fraction to number density.			*
 *   									*
 * NOTE:								*
 *   expects table* and ionfrac* to have been initialized by 		*
 *   readintables.c							*
 *   the analytic cooling routine is still in here as analytic_cool	*
 *									*
 * DOES:								*
 *   - uses NR locate to find the index corresponding to the temperature*
 *   - extrapolate or return analytic cooling term if outside table	*
 *   - calculate cooling term for each element/ion and sum over those	*
 *   - uses NR polint to interpolate over the tables. 			*
 *   - it populates X_el from SPHbody->abund[i] and puts elements in 	*
 *     ascending order of Z, and also sums over isotopes. Thus the order*
 *     of isotopes/elements in SPHbody->abund does not matter. 		*
 ************************************************************************/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <Msgs.h>
#include <singlio.h>
#include "physics_sph.h"
#include "units.h"
#include "cool.h"
#include "nrutil.h" /*ok to have this in*/

#ifndef M_1_PI
#define	M_1_PI 0.31830988618379067154
#endif


/*
 * take all command line outputs like error or status messages out for 
 * now. eventually should have a flag in the .ctl files to turn 
 * debugging on (i.e. with output messages) or off
 */

/*arrays in C: array[row-index][col-index]*/

double calc_lcool1(float abundarr[], int nparr[], int nnarr[], double temp, double rho, int Gridpts,int Nel,int extrapolate)
{
    //extern float **tablep;	/*global array with cooling values*/
    //extern float **ionfracp;	/*global array with ion fractions*/
    double lcool;	/*holds final cooling term;returned*/
    double *lcoolp;	/*pointer to lcool*/
    double dy,df;	/*measure of accuracy returned from interp.*/
    double *dyp;	/*pointer to dy*/
    double *dfp;	/*pointer to df*/
    double ioncool=0.0;		/*hold intermediate cooling term*/
    double fracn;	/*holds ion fraction returned by interp. */
    double *fracnp;	/*pointer to frac*/
    double fracneu;	/*fraction of neutral atoms per element*/
    double rowarr[2],interp[2],fracns[2];	/*temporary arrays for interp.*/
    double logtemp; 	/*log(temp) for ionfracp interpolation*/
    float temps[Gridpts];
    float X_el[Nel];		/*element fraction, i.e. number density*/
    float ndens = 0.;	/* total number density */
    float n_e = 0.;	/* total electron number density */
    float n_ion = 0.;	/* total ion number density */
    long j;	/*holds index returned by locate routine*/
    long *jp;	/*pointer to j*/
    long j_prev;
    //int nparr[Nel], nnarr[Nel];	/*array for A,N numbers for each isotope*/
    int n,m,N,index; 	/*some indices for looping and arrays*/
	
    /*initialize things so I don't get stupid warnings all the time:*/
    j=-999;
    dy=-999;
    jp=&j;
    dyp=&dy;
    lcoolp=&lcool;
    dfp=&df;
    fracnp=&fracn;

    /* we're doing temperature in logspace */
    logtemp=log10(temp);

    /*locate the indices of the table;
      same for both tables as they go over the same range/grid points 
    */
    for( n = 0; n < Gridpts; n++)
	temps[n] = ionfracp[0][n];

    locate(&temps[0], Gridpts, logtemp, jp); 
    j_prev=j;
	
    /*if we're outside the table, do analytic cooling if extrapolate=0
      or extrapolate if extrapolate=1:
      extrapolation is still a little funky - CIE
      above table= j=-2, below table= j=-99
    */
    if ((j==-2) || (j==-99))
    {
	if (extrapolate)
	{
	    /*printf("extrapolating......\n");*/
	    /*reset j for extrapolation*/
	    if (j==-2) j=0;
	    /*now do interpolation, should automatically 
	     * extrapolate as set up below*/

	    /*reset j for extrapolation*/
	    if (j==-99) j=Gridpts-1;
	    /*now do interpolation, should automatically 
	     * extrapolate as set up below*/
	}
	else
	    return analytic_cool(temp); /*we're done here*/
    } 

    /*technically, need to interpolate over every ion in every element
      then weigh by equil fraction at current T, and mix according to 
      composition.
      So, create array with element fractions from the SPHbody array of
      isotope fractions (mole fractions?)
    */

    for( n = 0; n < Nel; n++) X_el[n] = 0.; /*initialize all to zero*/

    for( n = 0; n < NISO; n++) {/*sum individual isotopes*/
        /*exclude bare neutrons;in tablep/ionfracp index=0 is H*/
	/* this also puts X_el in order of ascending Z, if an element
	 * does not exist in abund[i], it is just zero in X_el */
        if(nparr[n] > 0) {
           X_el[ nparr[n]-1 ] += abundarr[n] * rho *
               N_AVOG / ((float)(nnarr[n] + nparr[n]));
	   ndens += X_el[ nparr[n]-1]; /*this double-counts isotopes*/
/*
	   n_e += (double)(rho * massCF /(lengthCF*lengthCF*lengthCF) * 
	       N_AVOG / (double)(nparr[j] + nnarr[j]) * abundarr[j]) * 
               (nparr[j]); 
*/
        }
    }

    for ( n = 0; n < Nel; n++) {
	    N = n+1;
	    /*loop over ions for each element,dont skip any*/
        for ( m = 0; m < (N+1); m++) {
            index = (int)(N * (float)( (N+1)/2 ) ); /* to preempt int division probs */
	    /*re-assign table values so interpolation can be done in 
	      double precision (table is floats). */
	    rowarr[0] = ionfracp[0][ j ];
	    rowarr[1] = ionfracp[0][ j + 1 ];
	    interp[0] = tablep[ index-1+m ][ j ];
	    interp[1] = tablep[ index-1+m ][ j + 1 ];

	    /*interpolate*/
	    polint(rowarr,interp,2,logtemp,lcoolp,dyp);

	    /*reset value if extrapolated to negative value*/
	    if (lcool<0.0)
	    {
		if ((interp[0]-interp[1]) <0) lcool=0.0;
		/*else fracn=1.0;*//*this case shouldn't happen*/
	    }
		
	    fracns[0] = ionfracp[ index+m ][ j ];
	    fracns[1] = ionfracp[ index+m ][ j + 1 ];

	    /*interpolate*/
	    polint(rowarr,fracns,2,logtemp,fracnp,dfp);

        if(m == 0) fracneu = fracn;

	    /*reset value if extrapolated to unphysical value*/
	    if (fracn<0.0)
	    {
	        if ((fracns[0]-fracns[1]) <0) fracn=1.0e-40;
		else fracn=1.0;
	    }

            n_e += X_el[n] * (double)(m) * fracn;

            /* fracn weighs it by fraction in this ionization state */
            /* X_el[n]/ndens weighs it by fractional abundance of this element */
	    ioncool += lcool * fracn * X_el[n] / ndens; 
	}
        /* sum total ion number density, excluding neutral ions */
        n_ion += X_el[n] * ((double)1.0 - (double)fracneu);
    }
    /*lcool is in erg*cm^3/s, X_el/ndens is n_i/n_tot */
    /*multiplying by n_e * n_totalion hopefully converts to erg/cm^3/s */
    ioncool = ioncool * n_e * n_ion;

    return ioncool;
} 

double find_ne(float abundarr[], int nparr[], int nnarr[] ,double temp, double rho, int Gridpts,int Nel)
{
/* could also calculate n_ion and set a flag that determines which one is returned */
    extern float **ionfracp;	/*global array with ion fractions*/
    double dy,df;	/*measure of accuracy returned from interp.*/
    double *dyp;	/*pointer to dy*/
    double *dfp;	/*pointer to df*/
    double fracn;	/*holds ion fraction returned by interp. */
    double *fracnp;	/*pointer to frac*/
    double fracneu;	/*fraction of neutral atoms per element*/
    double rowarr[2],fracns[2];	/*temporary arrays for interp.*/
    double logtemp; 	/*log(temp) for ionfracp interpolation*/
    float temps[Gridpts];
    float X_el[Nel];		/*element fraction, i.e. number density*/
    double nelectron = 0.;	/* total electron number density */
    long j;	/*holds index returned by locate routine*/
    long *jp;	/*pointer to j*/
    long j_prev;
    //int nparr[Nel], nnarr[Nel];	/*array for A,N numbers for each isotope*/
    int n,m,N,index; 	/*some indices for looping and arrays*/
	
    /*initialize things so I don't get stupid warnings all the time:*/
    j=-999;
    dy=-999;
    jp=&j;
    dyp=&dy;
    dfp=&df;
    fracnp=&fracn;

    /* we're doing temperature in logspace */
    logtemp=log10(temp);
    if( logtemp < ionfracp[0][0] ) logtemp = 4.0;
    if( logtemp > ionfracp[0][Gridpts-1] ) logtemp = 9.0;

    /*locate the indices of the table;
      same for both tables as they go over the same range/grid points 
    */
    for( n = 0; n < Gridpts; n++) {
        temps[n] = ionfracp[0][n];
    }

    locate(&temps[0], Gridpts, logtemp, jp); 
    j_prev=j;
	
    /*if we're outside the table, do analytic cooling if extrapolate=0
      or extrapolate if extrapolate=1:
      extrapolation is still a little funky - CE
      above table= j=-2, below table= j=-99
      edit: should the above line be vice versa?? -CIE
    */
    if (j==-2) j=0;
    if (j == -99) j=Gridpts-1;

    for( n = 0; n < Nel; n++) X_el[n] = 0.; /*initialize all to zero*/

    for( n = 0; n < NISO; n++) {/*sum individual isotopes*/
        /* exclude bare neutrons;in tablep/ionfracp index=0 is H*/
	/* this also puts X_el in order of ascending Z, if an element
	 * does not exist in abund[i], it is just zero in X_el */
        if(nparr[n] > 0) {
           X_el[ nparr[n]-1 ] += abundarr[n] * rho *
               N_AVOG / ((float)(nnarr[n] + nparr[n]));
        }
    }

    nelectron = 0.0;

    for ( n = 0; n < Nel; n++)	
    {
	N = n+1;
	/*loop over ions for each element,dont skip any*/
        for ( m = 0; m < (N+1); m++) {
            index = (int)(N * (float)( (N+1)/2 ) ); /* to preempt int division probs */
	        /* re-assign table values so interpolation can be done in 
	           double precision (table is floats). */
	        rowarr[0] = ionfracp[0][ j ];
	        rowarr[1] = ionfracp[0][ j + 1 ];
		
	        fracns[0] = ionfracp[ index+m ][ j ];
	        fracns[1] = ionfracp[ index+m ][ j + 1 ];

	        /*interpolate*/
	        polint(rowarr,fracns,2,logtemp,fracnp,dfp);

	        /*reset value if extrapolated to unphysical value*/
	        if (fracn<0.0)
	        {
	            if ((fracns[0]-fracns[1]) <0) fracn=1.0e-40;
	    	    else fracn=1.0;
	        }

            nelectron += X_el[n] * (double)(m) * fracn;
	    }
    }
    if(nelectron != nelectron) nelectron = 0.0;
    if(nelectron < 0.0) nelectron = 0.0;
    return nelectron;
} /*end find_ne*/



double analytic_cool(double temp)
{
    /* From Chris's email; fit to Dalgarno and McCray (ARA&A
	 1972, 10, 375) and Sutherland and Dopita (ApJS, 88,
	 253) 
	 */
	double lcool;
	
	if (temp < 1.0e4)
		lcool = 1.0e-26 * exp(-1.0e5/temp) * sqrt(temp);//guessed; for O - CE
		//lcool = 1.0e-27 * exp(-1.0e2/temp) * sqrt(temp);
	else if (temp < 3.0e5)
		lcool=1.0e-21;
	else
		lcool=1.0e-21/(3.0*(log10(temp)-5.5)+1.0);
		//lcool=1.0e-21/(3.0*(log10(temp)-5.5)+1.0);
	
	return lcool;
}/*end analytic_cool*/



/* this bisection routine is from NR for C: */
/*"Given an array xx[1..n], and given an value x, returns a value
j such that x is between xx[j] and xx[j+1]. xx must be monotonic,
either increasing or decreasing. j=-2 or j=-99 is returned to indicate
that x is out of range.*/
/*#includes: none*/
/*call syntax: locate(&xx, N, x, j) */
/*~~~~~~~~~~~~~~~~WORKS!!! 04/23/2009~~~~~~~~~~~~~*/
void 
locate(float xx[], long Nel, float x, long *j)
{
   /*floats should be enough for finding correct indices*/
   long ju, jm, jl;
   int ascnd,sign;

   jl=-1;
   ju=Nel;
   /*check whether the table is in increasing (ascnd=1) or 
     decreasing (ascnd=0) order
   */
   ascnd=(xx[Nel-1] >= xx[0]);   /*what does this line do - check whether
                                   xx[N]>xx[1] or not -CE*/
   if (ascnd) sign=1;
   if (!ascnd) sign=-1;

   while (ju-jl > 1)
   {
      jm=(ju+jl) >> 1;
      if ((x >= xx[jm]) == ascnd)   
         jl=jm;
      else
         ju=jm;
   }
   /*if (true && true) && true (1 is true) then below (j<0) table*/
   if ( (((x - xx[0])*sign <0) && ((x-xx[Nel-1])*sign <0)) && 1)
	   *j=-2; 
   /*if (not false && not false) && true (1 is true) then above (j>Nel) table*/
   else if ( (!((x - xx[0])*sign <0) && !((x-xx[Nel-1])*sign <0)) && 1)
	   *j=-99; 
   else *j=jl;
} /*end locate*/



/*this is a polynomial interpolation routine from NR for C (S3.1):
"Given arrays xa[1..n] and and ya[1..n], and given a value x, this 
routine returns a value y, and an error estimate dy. If P(xP) is the 
polynomial of degree N-1 such that P(xa_i)=ya_i, i=1,...,n, then 
the returned value y=P(x).*/
/*#includes: math.h, "nrutil.h"*/
/*call syntax: polint(&xx[14],&yy[14],4,x,yp,dyp) for 4-point 
interpolation between tabulated points [14..17]*/
/*~~~~~~~~~~~~~~~~~ WORKS!!!! 04/23/2009 ~~~~~~~~~~~~~~~~~*/
void polint(double xa[], double ya[], int n, double x, double *y, double *dy)
{
   //ZERO!!! offset is assumed in all indices
   int i,m,ns=0,size;
   double den,dif,dift,ho,hp,w;
   double *c, *d;

   //printf("interpolating between %E and %E .... get ", ya[0], ya[1]);
   //printf("xa= {%.1E, %.1E} \n", xa[0], xa[1]);
   //printf("ya= {%.1E, %.1E} \n", ya[0], ya[1]);

   //if we're exactly at one grid point, just return that value*/
   if ((x-xa[0]) == 0.0) 
	   *y=ya[0];
   else if ((x-xa[1]) == 0.0) 
	   *y=ya[1];
   else	//else interpolate*/
   {
   dif=fabs(x-xa[0]);
   c=dvector(0,n-1);
   d=dvector(0,n-1);
   //find the index ns of the closest table entry
   for (i=0;i<n;i++)
   {
      if ((dift=fabs(x-xa[i])) < dif)
      {
         ns=i;
         dif=dift;
      }
      //initialize the tableau of c's and d's
      c[i]=ya[i];
      d[i]=ya[i];
   }
   //initial approximation to y
   *y=ya[ns--];
   //for each column of the tableau ...
   for (m=0;m<n-1;m++)
   {
      //... loop over current c's and d's and update them
      for (i=0;i<n-m-1;i++)
      {
         ho=xa[i]-x;
         hp=xa[i+m+1]-x;
         w=c[i+1]-d[i];
	 den=ho-hp;
         //error message: two input xa's are identical to within roundoff
         //if (den==0.0) nerror("Error in routine polint");
         den=w/den;
         //update c's and d's
         d[i]=hp*den;
         c[i]=ho*den;
      }
      /*After each column in the tableau is completed, we decide wich 
        correction, c or d, we want to add to our accumulating value of
        y, i.e., which path to take through the tableau - forking up or 
        down. We do this in such a way as to take the most "straight 
        line" route through the tableau to its apex, updating ns 
        accordingly to keep track of where we are. This route keeps the 
        partial approximations centered (insofar as possible) on the 
        target x. the last dy added is thus the error indication. */
      *y += (*dy=(2*(ns+1)<(n-(m+1)) ? c[ns+1] : d[ns--]));
      //printf("%E\n",*y);
   }
   //printf("%E\n", *y);
   free_dvector(d,0,n-1);
   free_dvector(c,0,n-1);
   }
} //end polint



/*bilinear polynomial interpolation in 2 dimensions from NR for C: */
/*Given arrays x1a[1..m] and x2a[1..n] of independent variables, and a 
submatrix of function values ya[1..m][1..n], tabulated at the grid points
defined by x1a and x2a; and given values x1 and x2 of the independent 
variables; this routine teturns as interpolated function value y, and an 
accuracy idication dy (based only on the interpolation in the x1 direction)*/
/*#includes: "nrutil.h" */
/* call syntax: polin2d(&x1a[jj],&x2a[kk],&yap[0][0],2,2,x1,x2,yp,dyp) 
 * where x1a,x1, refers to rows and jj is the row number from which to 
 * start the interpolation, x2a,x2 refers to columns and kk is the column
 * number from which to start the interpolation, 
 * and yap is an array of pointers (which point to the values
 * to interpolate) like so:
 * yap[0][0]=&table[jj][kk]
 * yap[1][0]=&table[jj][kk+1]
 * yap[0][1]=&table[jj+1][kk]
 * yap[1][1]=&table[jj+1][kk+1]
 * (really, whatever the order, polin2d treats the first index of yap as 
 * rows and interpolates over those first)
 */
/*~~~~~~~~~~~~~~~~ WORKS!!!! 04/23/2009 ~~~~~~~~~~~~~~*/
void
polin2d(double x1a[], double x2a[], double **ya, int m, int n, double x1, 
        double x2, double *y, double *dy)
{
   //ZERO!!! offset is assumed in all indices
   void polint(double xa[], double ya[], int n, double x, double *y, double *dy);

   int j;
   double *ymtmp;
 
   ymtmp=dvector(0,m-1);
   //loop over rows
   for (j=0;j<m;j++)
   {
      //interpolate over the 'rows' of the grid square containing the point
      //(x1,x2), i.e. over (x1a[j],x2). Put answer into temporary storage
      polint(x2a,ya[j],n,x2,&ymtmp[j],dy);
   }
   //do the final interpolation, i.e. sort of interpolate over the 
   //interpolated 'rows'
   polint(x1a,ymtmp,m,x1,y,dy);
   free_dvector(ymtmp,0,m-1);
} //end polin2d

