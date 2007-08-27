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
   char *basename;
   char outfile[256];
   char infile[256];
   int ret;
   int nx, ny;
   unsigned char *image;
   unsigned char *outimage;
   int max_pixel, min_pixel;
   float max_value, min_value;
   int i, j, k, start, end, stride;
   float minc, maxc;

   if (argc == 7) {
       basename = argv[1];
       start = atoi(argv[2]);
       end = atoi(argv[3]);
       stride = atoi(argv[4]);
       minc = atoi(argv[5]);
       maxc = atoi(argv[6]);
   } else {
       fprintf(stderr, "usage: %s basename start end stride minc maxc\n", 
	       argv[0]);
       exit(1);
   }

   for (k = start; k <= end; k += stride) {
       sprintf(infile, "%s.%04d", basename, k);
   
       sdfp = SDFopen(NULL, infile);
       if (sdfp == NULL) {
	   fprintf(stderr, "Can't SDFopen(%s)\n%s\n", infile, SDFerrstring);
	   exit(1);
       }
       
       if( SDFgetint(sdfp, "nx", &nx) )
	 Error("Sorry, you've got to have an \"nx\"\n");
       if( SDFgetint(sdfp, "ny", &ny) )
	 Error("Sorry, you've got to have an \"ny\"\n");
       if( SDFgetint(sdfp, "max_pixel", &max_pixel) )
	 Error("Sorry, you've got to have an \"max_pixel\"\n");
       if( SDFgetint(sdfp, "min_pixel", &min_pixel) )
	 Error("Sorry, you've got to have an \"min_pixel\"\n");
       if( SDFgetfloat(sdfp, "max_value", &max_value) )
	 Error("Sorry, you've got to have a \"max_value\"\n");
       if( SDFgetfloat(sdfp, "min_value", &min_value) )
	 Error("Sorry, you've got to have an \"min_value\"\n");

       fprintf(stderr, "%s %10.2f %10.2f\n", infile, min_value, max_value);

       image = Malloc(nx * ny);
       outimage = Malloc(nx * ny);

       ret = SDFseekrdvecs(sdfp, "value", 0, nx*ny, image, 1, NULL);

       if (ret) 
	 Error("SDFseekrdvecs failed, %s\n", SDFerrstring);
       SDFclose(sdfp);

       for (i = 0; i < nx; i++) {
	   for (j = 0; j < ny; j++) {
	       float x = image[i*ny+j];
	       x -= min_pixel;
	       x *= (max_value-min_value)/(max_pixel-min_pixel);
	       x += min_value-minc;
	       x *= (max_pixel-min_pixel)/(maxc-minc);
	       x += min_pixel;
	       if (x > max_pixel) x = max_pixel+1;
	       if (x < min_pixel) x = min_pixel-1;
	       /* invert top-to-bottom */
	       outimage[(nx-1-i)*ny+j] = (char)x;
	   }
       }
       Free(image);
       sprintf(outfile, "%s.%04d.raw", basename, k);
       Fopen(outfp, outfile, "w");
       
       Fwrite(outimage, 1, nx*ny, outfp);
       Fclose(outfp);
       Free(outimage);
   }
   exit(0);
}
