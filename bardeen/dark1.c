#include <stdio.h>
#include <stdlib.h>
#include "physics.h"
#include "SDF.h"
#include "SDFwrite.h"

void
main(int argc, char *argv[])
{
    outbody b;

    b.mass = 1.0;
    
    b.pos[0] = 1.0;
    b.pos[1] = 1.0;
    b.pos[2] = 0.0;
    
    b.vel[0] = 0.0;
    b.vel[1] = -5.0;
    b.vel[2] = 0.0;

    b.ident = 11569;

    SDFwrite("dark1", 1, 1, &b, sizeof(outbody),
	     OUTBODYDESC,
	     "npart", SDF_INT, 1,
	     NULL);

    exit(0);
}
