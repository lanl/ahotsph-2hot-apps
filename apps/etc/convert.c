#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <SDF.h>
#include "bigmalloc.h"
#include "macr.h"
#include "protos.h"

int
main(int argc, char *argv[])
{
   SDF *sdfp;
   FILE *outfp;
   char *infile;
   short image;
   int nx, ny;
   int ret;

   infile = argv[1];

   Fopen(outfp, "topo_bs.img", "w");
   sdfp = SDFopen("img.hdr", infile);
   
   if( SDFgetint(sdfp, "nx", &nx) )
     Error("Sorry, you've got to have an \"nx\"\n");
   
   if( SDFgetint(sdfp, "ny", &ny) )
     Error("Sorry, you've got to have an \"ny\"\n");
   
   for (i = 0; i < nx*ny; i++) {
     ret = SDFseekrdvecs(sdfp, "image", i, 1, &image, NULL);
     if (ret) 
       Error("SDFseekrdvecs failed, %s\n", SDFerrstring);
     Fwrite(image, sizeof(short), 1, outfp);
   }

   Fclose(outfp);

   exit(0);
}
