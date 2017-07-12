#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<strings.h>
#include"singlio.h"
/*
#include<assert.h>
#include<time.h>
#include<errno.h>
#include"error.h"
#include"fastflpt.h"
#include"Msgs.h"
#include"physics.h"
#include"stk.h"
#include"vop.h"
#include"mpmy.h"
#include"timers.h"
*/
#include"physics_sph.h"
#include"cool.h"
#include"nrutil.h"

/******************************************************************************
 * the two arrays that hold the table values for the ion fractions (ionfracp) *
 * and cooling terms (tablep) are now 2D arrays, with the first index going   *
 * over the ions (order: H (Z=1) first, Zn (Z=30) last; for each element      *
 * neutral atom first, bare ion last), and the second index going over the    *
 * the data points for each temperature gridpoint. The correct row (ion) can  *
 * be found with: [ Nel*(Nel+1)/2 -1 + ionstate] ( -1 because arrays in C     *
 * start at 0), where ionstate is the ionization state (0 for neutral, Z+1    *
 * for bare ion), and Nel is the number of the element (=Z).                  *
 ******************************************************************************/
/*
extern float **tablep;
extern float **ionfracp;
*/

/* return Gridpts and Nel to calling function*/
void init_CoolTable(int *Gridpts, int *Nel)
{
    FILE *File1p;	/*pointer to file with cooling curves*/
    FILE *file2p;	/*pointer to file with ion fractions*/
    long lSize;		/*holds file size (number of characters in file)*/
    int i,j,k;		/*indices for looping through arrays*/
    int counter1 = 0, counter2 = 0;
    int index;		/*index to access correct array element*/
    int myint;		/*holds un-needed integers read in from files*/
    int tot_ion;	/* total number of ions in database (w/ bare ions)*/
    char mychar,myline[50];	/*holds new-lines and text read in from files*/

    /*open table with cooling curves*/
    File1p = fopen(
           "CHIANTI-COOLING.dat", "r");
           /*"/home/cellinge/SNSPH.dir/tree16/sph+nln/CHIANTI-COOLING.dat", "r");*/
    if (File1p == NULL) 
	singlPrintf("error opening cooling curves: \n");

    /*open table with ion fractions*/
    file2p = fopen(
           "mazzotta_etal_9.ioneq","r");
           /*"/home/cellinge/SNSPH.dir/tree16/sph+nln/mazzotta_etal_9.ioneq","r");*/
    if (file2p == NULL) 
	singlPrintf("error opening ion fractions: \n");

    fscanf(file2p, "%i %i", Gridpts,Nel);
     
    /* (Nel*(Nel+1)/2 + Nel is total number of ions, incl. bare ions */
    tot_ion = (*Nel) * ( (*Nel)+3 ) / 2;
    /*singlPrintf("total ions: %d\n", tot_ion);*/


    /* for ionfracp, need (Nel*(Nel+1)/2 + Nel+1) by Gridpts array */
    ionfracp = (float**)malloc( (tot_ion + 1) * sizeof(float *) );

    for ( i = 0; i < (tot_ion + 1); i++) {
	ionfracp[i] = (float *)malloc( (*Gridpts) * sizeof(float) );
        for ( j = 0; j < (*Gridpts); j++)
            ionfracp[i][j] = 0.0;
    }


    /* for tablep, need (Nel*(Nel+1)/2) by Gridpts array */ 
    tablep = (float **)malloc( (tot_ion) * sizeof(float *) );

    for ( i = 0; i < (tot_ion); i++) {
	tablep[i] = (float *)malloc( (*Gridpts) * sizeof(float *) );
        for ( j = 0; j < (*Gridpts); j++)
            tablep[i][j] = 0.0;
    }

    fgets(myline,50,File1p);/*read in first line of text in cooling curves*/
    fgets(myline,50,File1p);/*read in second line of text in cooling curves*/

    /* since the ions are from H to Zn in ascending order, don't need
     * Z of element (?) */
    for( i = 0; i < (*Nel); i++) 
	    fscanf(File1p, "%*i");


    mychar=fgetc(File1p);/*read in extra new-line in cooling curves*/
    fgets(myline,50,File1p);/*read in line "temperatures...." in cooling curves*/


    /* we're getting log(T) from ionfractions, so skip T from cooling table*/
    for ( i = 0; i < (*Gridpts); i++) 
	    fscanf(File1p, "%*g");

    /*read in log(temperatures) for ion fractions*/
    for ( i = 0; i < (*Gridpts); i++) 
	    fscanf(file2p,"%4g",&ionfracp[0][i]);
	

    fgets(myline,50,File1p);/*get trailing new-line in cooling curves*/
    fgets(myline,50,File1p);/*get trailing new-line in cooling curves*/


    /*loop over elements*/
    for ( j = 0; j < (*Nel); j++) { 
	/*get line with element name*/
        fgets(myline,50,File1p);

        /*loop over ions in current element*/
	for ( k = 0; k <= (j+1); k++) {
	    fscanf(file2p,"%3i",&myint);
	    fscanf(file2p,"%3i",&myint);

            for ( i = 0; i < (*Gridpts); i++) {

    	        fscanf(File1p,"%13E",&tablep[ counter1 ][i]);

		/*add one to counter2 since the 1st contains log(T) */
	        fscanf(file2p,"%10E",&ionfracp[ counter2+1 ][i]);

	/* in CHIANTI-file, (neutrals?) bare ions seem to be missing, but are
	 * present in Mazzotta et al. ion fractions file. In the above set-up
	 * the bare ions from the ion-fractions file are read in, and the ions
	 * in tablep and ionfracp are lined up (i.e. each row of ionfrac
	 * and tablep contain the same ion of the same element), and the 
	 * bare ions are read in as zeros in tablep).
	 */
	    }

            /*printf("counter: %4d, j=%2d, k=%2d, tab= %.6E  ion=%.6E\n",counter1,j,k,tablep[counter1][12],ionfracp[counter1][12]);*/
	    counter1++;
	    counter2++;
	}
	mychar=fgetc(File1p);/*get first trailing new-line*/
	mychar=fgetc(File1p);/*get second trailing new-line*/
    }

    /*close file*/
    fclose(File1p);
    fclose(file2p);
    /*DO NOT free(tablep); UNTIL THE END OF THE WHOLE PROGRAM!!!!*/
}

