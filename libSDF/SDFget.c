#include <stddef.h>
/* stdio only needed for sprintf proto */
#include <stdio.h>
#include "error.h"
#include "SDF.h"
/* 

These routines return 0 on success and 1 on failure.  They set
SDFerrstring if they think something strange is happening.  They
don't change the value under the VALUE argument until they are sure
that everything is OK, so it's possible (but what would Henry Spencer
say?) to use them by setting a default value, then calling the routine,
and not worrying about the error condition.

*/

/* An abbreviation used many times... */

#define SDFget(sdfp, name, addr) \
    if(SDFseekrdvecs(sdfp, name, 0, 1, addr, 0, NULL)) return -1

/* I'll resist the temptation to make one macro to cover all these cases */
/* Should we check for loss of precision in float too? */
int 
SDFgetfloat(SDF *sdfp, char *name, float *value)
{
    double double_value;
    int int_value;
    float float_value;
    
    if( sdfp == NULL || !SDFhasname(name, sdfp) ){
	return -1;
    }
    switch(SDFtype(name, sdfp)){
    case SDF_FLOAT:
	SDFget(sdfp, name, &float_value);
	*value = float_value;
	return 0;
    case SDF_DOUBLE:
	SDFget(sdfp, name, &double_value);
	*value = (float) double_value;
	return 0;
    case SDF_INT:
	SDFget(sdfp, name, &int_value);
	*value = (float) int_value;
	return 0;
    default:
	sprintf(SDFerrstring, 
		"SDFgetfloat: \"%s\" must be either float or int.\n", name);
	return -1;
    }
}

int
SDFgetdouble(SDF *sdfp, char *name, double *value)
{
    double double_value;
    float float_value;
    int int_value;

    if( sdfp == NULL || !SDFhasname(name, sdfp) ){
	return -1;
    }
    switch(SDFtype(name, sdfp)){
      case SDF_DOUBLE:
	SDFget(sdfp, name, &double_value);
	*value = double_value;
	return 0;
      case SDF_FLOAT:
	SDFget(sdfp, name, &float_value);
	*value = (double) float_value;
	return 0;
      case SDF_INT:
	SDFget(sdfp, name, &int_value);
	*value = (double) int_value;
	return 0;
      default:
	sprintf(SDFerrstring, 
		"SDFgetdouble: \"%s\" must be either float or int.\n", name);
	return -1;
    }
}

int
SDFgetint(SDF *sdfp, char *name, int *value)
{
    int int_value;
    float float_value;
    double double_value;

    if( sdfp == NULL || !SDFhasname(name, sdfp) ){
	return -1;
    }
    switch(SDFtype(name, sdfp)){
      case SDF_INT:
	SDFget(sdfp, name, &int_value);
	*value = int_value;
	return 0;
      case SDF_DOUBLE:
	SDFget(sdfp, name, &double_value);
	int_value = (int) double_value;
	if ((double)int_value != double_value){
	    sprintf(SDFerrstring, 
		    "SDFgetint: \"%s\" has lost precision\n", name);
	    return -1;
	}
	*value = int_value;
	return 0;
      case SDF_FLOAT:
	SDFget(sdfp, name, &float_value);
	int_value = (int) float_value;
	if ((float)int_value != float_value){
	    Warning("SDFgetint: \"%s\" has lost precision\n", name);
	    return -1;
	}
	*value = int_value;
	return 0;
      default:
	sprintf(SDFerrstring, "SDFgetint: \"%s\" must be either float or int.\n", name);
	return -1;
    }
}

int
SDFgetstring(SDF *sdfp, const char *name, char *string, int size){
    int len;
    int ret;

    if( sdfp == NULL ){
	sprintf(SDFerrstring, "SDFgetint: NULL sdfp\n");
	return -1;
    }
    if( !SDFhasname(name, sdfp) ){
	return -1;
    }
    if( SDFtype(name, sdfp) != SDF_CHAR ){
	return -1;
    }
    if( (len=SDFarrcnt(name, sdfp)) >= size ){
	return -1;
    }

    ret = SDFseekrdvecs(sdfp, name, 0, 1, string, 0, NULL);
    if( string[len-1] != '\0' ){
	sprintf(SDFerrstring, 
		"Read non-terminated string in SDFgetstring(\"%s\")\n", name);
	string[len] = '\0';
    }
    return ret;
}

    
