#include <stdlib.h>
#include "gc.h"
#include "mpmy.h"
#include "Msgs.h"
#include "physics_sph.h"
#include "SDF.h"
#include "SDFreadf.h"
#include "singlio.h"


/* Added ability to read in rho from SDF file (to avoid calculating it
   myself */

SDF *SPHReadf(char *name, SPHbody **btabp, int *gnobjp, int *nobjp)
{
    SDF *sdfp;
    int massconf, gravmassconf, xconf, yconf, zconf;
    int vxconf, vyconf, vzconf;
    int hconf, uconf, identconf, windidconf;
    int rhoconf;
    SPHbody *btab; 
    int nobj, gnobj, i;
    
    singlPrintf("Reading \"%s\"\n", name);
    sdfp = SDFreadf(name, (void **)btabp, gnobjp, nobjp, sizeof(SPHbody),
		    "mass", offsetof(SPHbody, mass), &massconf,
		    "grav_mass", offsetof(SPHbody, grav_mass), &gravmassconf,
		    "rho", offsetof(SPHbody, rho), &rhoconf,
		    "x", offsetof(SPHbody, pos[0]), &xconf,
		    "y", offsetof(SPHbody, pos[1]), &yconf,
		    "z", offsetof(SPHbody, pos[2]), &zconf,
		    "vx", offsetof(SPHbody, vel[0]), &vxconf,
		    "vy", offsetof(SPHbody, vel[1]), &vyconf,
		    "vz", offsetof(SPHbody, vel[2]), &vzconf,
		    "u", offsetof(SPHbody, u), &uconf,
		    "h", offsetof(SPHbody, h), &hconf,
		    "ident", offsetof(SPHbody, ident), &identconf,
		    "windid", offsetof(SPHbody, windid), &windidconf,
		    NULL);
    nobj = *nobjp;
    gnobj = *gnobjp;
    btab = *btabp;
    
    Msgf(("Data read, SPHnobj=%d, SPHgnobj=%d\n", *nobjp, *gnobjp));
    Msgf(("Nproc:%d, Procnum: %d, Doc: %d\n",
	  MPMY_Nproc(), MPMY_Procnum(), ilog2(MPMY_Nproc())));
    
    if (massconf==0 || xconf==0 || yconf==0 || zconf==0) {
	SinglError("Could not find %s %s %s %s in data file!\n",
		   (massconf==0)? "mass" : "",
		   (xconf==0)? "x" : "",
		   (yconf==0)? "y" : "",
		   (zconf==0)? "z" : "");
    }
    
    if (gravmassconf==0) {
	singlPrintf("Grav_mass not in file, using sphmass instead.");
	for(i=0; i<nobj; i++){
	    btab[i].grav_mass = btab[i].mass;
	}
    }
    
    if (identconf == 0)
	SinglWarning("No \"ident\" in file\n");
    
    if (windidconf == 0)
	SinglWarning("No \"windid\" in file; are you using wind source?\n");
    
    return sdfp;
}
