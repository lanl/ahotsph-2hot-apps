#include <math.h>
#include <stdlib.h>
#include "SDF.h"
#include "bigmalloc.h"
#include "macr.h"

typedef struct {
  float strength;
  float pos[3];
  int ident;
} outbody;

void
main(int argc, char **argv)
{
  outbody *btab, *p;
  int ntheta = 32;
  int nphi = 64;
  int nobj;
  float dtheta, dphi;
  float *costh, *sinth, *cosph, *sinph;
  int i, j;
  float x;
  FILE *fp;
  char outname[64] = "fmmsrc.1";

  if (argc == 4) {
      strcpy(outname, argv[1]);
      ntheta = atoi(argv[2]);
      nphi = atoi(argv[3]);
  }
  nobj = ntheta * nphi;

  btab = Malloc(nobj * sizeof(outbody));
  costh = Malloc(ntheta * sizeof(float));
  sinth = Malloc(ntheta * sizeof(float));
  cosph = Malloc(nphi * sizeof(float));
  sinph = Malloc(nphi * sizeof(float));
  
  dtheta = M_PI/ntheta;
  for (i = 0; i < ntheta; i++) {
    x = (0.5 + i) * dtheta;
    costh[i] = cos(x);
    sinth[i] = sin(x);
  }

  dphi = 2.0*M_PI/nphi;
  for (i = 0; i < nphi; i++) {
    x = i * dphi;
    cosph[i] = cos(x);
    sinph[i] = sin(x);
  }

  p = btab;
  for (i = 0; i < nphi; i++) {
    for (j = 0; j < ntheta; j++) {
      p->strength = 1.0;
      p->pos[0] = sinth[j]*cosph[i];
      p->pos[1] = sinth[j]*sinph[i];
      p->pos[2] = costh[j];
      p->ident = p-btab;
      p++;
    }
  }

  if (p - btab != nobj) Error("nobj incorrect\n");

  Fopen(fp, outname, "w");
  fprintf(fp, 
	  "# SDF\n"
	  "parameter byteorder = 0x%x;\n"
	  "int npart = %d;\n"
	  "struct {\n"
	  "\tfloat strength;             /* strength of body */\n"
	  "\tfloat x, y, z;              /* position of body */\n"
	  "\tunsigned int ident;         /* unique identifier */\n"
	  "}[%d];\n"
          "#\f\n"
	  "# SDF-EOH\n",
	  SDFcpubyteorder(), nobj, nobj);
  Fwrite(btab, sizeof(outbody), nobj, fp);
  Fclose(fp);
}
