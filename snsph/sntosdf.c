#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include "bigmalloc.h"
#include "macr.h"
#include "error.h"
#include "byteswap.h"

static unsigned char junk[4] = {0x12, 0x34, 0x56, 0x78};
static unsigned int *cpubyteorder = (unsigned int *)&junk[0];

#define Fread_swap(a, b, c, d) \
  Fread(a, b, c, d); \
  if (swap) Byteswap(b, c, a, a)

void
main(int argc, char *argv[])
{
    FILE *infp;
    FILE *fp;
    char *in_name, *out_name;
    int npart, this_npart;
    float t, gamma, tkin, tterm;
    float *h, *x, *y, *vx, *vy, *u, *mass, *rho;
    float *abar, *temp, *ye, *xp, *xn;
    int *ifleos;
    float *ynue, *ynueb, *ynux, *unue, *unueb, *unux, *ufreez;
    float *pr, *u2, *cent2, *te, *teb, *tx;
    float xmcore, rb, ftrape, ftrapb, ftrapx;
    float *kepcel, *xpf, *pvar2, *pvar3, *pvar4;
    int fortran_crap;
    int i, j;
    int ndumps;
    char tmp[256];
    int swap = 0;

    if (argc == 3) {
	in_name = argv[1];
	out_name = argv[2];
	ndumps = 1;
    } else if (argc == 4) {
	in_name = argv[1];
	out_name = argv[2];
	ndumps = atoi(argv[3]);
    } else {
	fprintf(stderr, "Usage: %s infile outfile [ndumps]\n", argv[0]);
	exit(1);
    }

    Fopen(infp, in_name, "r");

    for (j = 0; j < ndumps; j++) {
      if (ndumps == 1) {
	strcpy(tmp, out_name);
      } else {
	sprintf(tmp, "%s.%d", out_name, j);
      }
      Fopen(fp, tmp, "w");
      Fread(&fortran_crap, sizeof(int), 1, infp);
      Fread(&this_npart, sizeof(int), 1, infp);
      if (swap == 0 && (this_npart < 0 || this_npart > 1000000)) {
	swap = 1;
	fprintf(stderr, "Npart is %d, trying byte-swap\n", this_npart);
      }
      if (swap) Byteswap(sizeof(int), 1, &this_npart, &this_npart);
      if (this_npart < 0 || this_npart > 1000000) {
	fprintf(stderr, "Bad Npart value\n");
	exit(1);
      }
      if (this_npart > 0) {
	npart = this_npart;	/* rewritten files get npart set to zero */
      }
      Fread_swap(&t, sizeof(float), 1, infp);
      Fread_swap(&gamma, sizeof(float), 1, infp);
      fprintf(stderr, "npart %d, t %f, gamma %f\n",
	      npart, t, gamma);
      
      h = Malloc(npart * sizeof(float));
      x = Malloc(npart * sizeof(float));
      y = Malloc(npart * sizeof(float));
      vx = Malloc(npart * sizeof(float));
      vy = Malloc(npart * sizeof(float));
      u = Malloc(npart * sizeof(float));
      mass = Malloc(npart * sizeof(float));
      abar = Malloc(npart * sizeof(float));
      rho = Malloc(npart * sizeof(float));
      temp = Malloc(npart * sizeof(float));
      ye = Malloc(npart * sizeof(float));
      xp = Malloc(npart * sizeof(float));
      xn = Malloc(npart * sizeof(float));
      ifleos = Malloc(npart * sizeof(int));
      ynue = Malloc(npart * sizeof(float));
      ynueb = Malloc(npart * sizeof(float));
      ynux = Malloc(npart * sizeof(float));
      unue = Malloc(npart * sizeof(float));
      unueb = Malloc(npart * sizeof(float));
      unux = Malloc(npart * sizeof(float));
      ufreez = Malloc(npart * sizeof(float));
      pr = Malloc(npart * sizeof(float));
      u2 = Malloc(npart * sizeof(float));
      cent2 = Malloc(npart * sizeof(float));
      te = Malloc(npart * sizeof(float));
      teb = Malloc(npart * sizeof(float));
      tx = Malloc(npart * sizeof(float));
      kepcel = Malloc(npart * sizeof(float));
      xpf = Malloc(npart * sizeof(float));
      pvar2 = Malloc(npart * sizeof(float));
      pvar3 = Malloc(npart * sizeof(float));
      pvar4 = Malloc(npart * sizeof(float));
      
      Fread_swap(h, sizeof(float), npart, infp);
      Fread_swap(&tkin, sizeof(float), 1, infp);
      Fread_swap(&tterm, sizeof(float), 1, infp);
      Fread_swap(x, sizeof(float), npart, infp);
      Fread_swap(y, sizeof(float), npart, infp);
      Fread_swap(vx, sizeof(float), npart, infp);
      Fread_swap(vy, sizeof(float), npart, infp);
      Fread_swap(u, sizeof(float), npart, infp);
      Fread_swap(mass, sizeof(float), npart, infp);
      Fread_swap(abar, sizeof(float), npart, infp);
      Fread_swap(rho, sizeof(float), npart, infp);
      Fread_swap(temp, sizeof(float), npart, infp);
      Fread_swap(ye, sizeof(float), npart, infp);
      Fread_swap(xp, sizeof(float), npart, infp);
      Fread_swap(xn, sizeof(float), npart, infp);
      Fread_swap(ifleos, sizeof(int), npart, infp);
      Fread_swap(ynue, sizeof(float), npart, infp);
      Fread_swap(ynueb, sizeof(float), npart, infp);
      Fread_swap(ynux, sizeof(float), npart, infp);
      Fread_swap(unue, sizeof(float), npart, infp);
      Fread_swap(unueb, sizeof(float), npart, infp);
      Fread_swap(unux, sizeof(float), npart, infp);
      Fread_swap(ufreez, sizeof(float), npart, infp);
      Fread_swap(pr, sizeof(float), npart, infp);
      Fread_swap(u2, sizeof(float), npart, infp);
      Fread_swap(cent2, sizeof(float), npart, infp);
      Fread_swap(te, sizeof(float), npart, infp);
      Fread_swap(teb, sizeof(float), npart, infp);
      Fread_swap(tx, sizeof(float), npart, infp);
      Fread_swap(&xmcore, sizeof(float), 1, infp);
      Fread_swap(&rb, sizeof(float), 1, infp);
      Fread_swap(&ftrape, sizeof(float), 1, infp);
      Fread_swap(&ftrapb, sizeof(float), 1, infp);
      Fread_swap(&ftrapx, sizeof(float), 1, infp);
      Fread_swap(kepcel, sizeof(float), npart, infp);
      Fread_swap(xpf, sizeof(float), npart, infp);
      Fread_swap(pvar2, sizeof(float), npart, infp);
      Fread_swap(pvar3, sizeof(float), npart, infp);
      Fread_swap(pvar4, sizeof(float), npart, infp);
      Fread_swap(&fortran_crap, sizeof(int), 1, infp);
      
      fprintf(fp, "# SDF\n");
      fprintf(fp, "char name[] = \"sntosdf on file %s\";\n", in_name);
      fprintf(fp, "int npart = %d;\n", npart);
      fprintf(fp, "float tpos = %g;\n", t);
      fprintf(fp, "float Gamma = %g;\n", gamma);
      fprintf(fp, "float tkin = %g;\n", tkin);
      fprintf(fp, "float tterm = %g;\n", tterm);
      fprintf(fp, "float xmcore = %g;\n", xmcore);
      fprintf(fp, "float rb = %f;\n", rb);
      fprintf(fp, "float ftrape = %g;\n", ftrape);
      fprintf(fp, "float ftrapb = %g;\n", ftrapb);
      fprintf(fp, "float ftrapx = %g;\n", ftrapx);
      fprintf(fp, "int ndim = %d;\n", 2);
      fprintf(fp, "parameter byteorder = 0x%x;\n", *cpubyteorder);
      fputs  ("struct {\n", fp);
      fputs  ("\tfloat mass;\n", fp);
      fputs  ("\tfloat x, y;\n", fp);
      fputs  ("\tfloat vx, vy;\n", fp);
      fputs  ("\tfloat u;\n", fp);
      fputs  ("\tfloat h;\n", fp);
      fputs  ("\tfloat rho;\n", fp);
      fputs  ("\tint ident;\n", fp);
      fputs  ("\tfloat abar;\n", fp);
      fputs  ("\tfloat temp;\n", fp);
      fputs  ("\tfloat ye;\n", fp);
      fputs  ("\tfloat xp;\n", fp);
      fputs  ("\tfloat xn;\n", fp);
      fputs  ("\tint ifleos;\n", fp);
      fputs  ("\tfloat ynue;\n", fp);
      fputs  ("\tfloat ynueb;\n", fp);
      fputs  ("\tfloat ynux;\n", fp);
      fputs  ("\tfloat unue;\n", fp);
      fputs  ("\tfloat unueb;\n", fp);
      fputs  ("\tfloat unux;\n", fp);
      fputs  ("\tfloat ufreez;\n", fp);
      fputs  ("\tfloat pr;\n", fp);
      fputs  ("\tfloat u2;\n", fp);
      fputs  ("\tfloat p2;\n", fp);
      fputs  ("\tfloat p3;\n", fp);
      fputs  ("\tfloat p4;\n", fp);
      fputs  ("}", fp);
      fprintf(fp, "[%d];\n", npart);
      fputs  ("#\f\n", fp);
      fputs  ("# SDF-EOH\n", fp);
      
      for (i = 0; i < npart; i++) {
	Fwrite(mass+i, sizeof(float), 1, fp);
	Fwrite(x+i, sizeof(float), 1, fp);
	Fwrite(y+i, sizeof(float), 1, fp);
	Fwrite(vx+i, sizeof(float), 1, fp);
	Fwrite(vy+i, sizeof(float), 1, fp);
	Fwrite(u+i, sizeof(float), 1, fp);
	Fwrite(h+i, sizeof(float), 1, fp);
	Fwrite(rho+i, sizeof(float), 1, fp);
	Fwrite(&i, sizeof(int), 1, fp);
	Fwrite(abar+i, sizeof(float), 1, fp);
	Fwrite(temp+i, sizeof(float), 1, fp);
	Fwrite(ye+i, sizeof(float), 1, fp);
	Fwrite(xp+i, sizeof(float), 1, fp);
	Fwrite(xn+i, sizeof(float), 1, fp);
	Fwrite(ifleos+i, sizeof(int), 1, fp);
	Fwrite(ynue+i, sizeof(float), 1, fp);
	Fwrite(ynueb+i, sizeof(float), 1, fp);
	Fwrite(ynux+i, sizeof(float), 1, fp);
	Fwrite(unue+i, sizeof(float), 1, fp);
	Fwrite(unueb+i, sizeof(float), 1, fp);
	Fwrite(unux+i, sizeof(float), 1, fp);
	Fwrite(ufreez+i, sizeof(float), 1, fp);
	Fwrite(pr+i, sizeof(float), 1, fp);
	Fwrite(u2+i, sizeof(float), 1, fp);
	Fwrite(pvar2+i, sizeof(float), 1, fp);
	Fwrite(pvar3+i, sizeof(float), 1, fp);
	Fwrite(pvar4+i, sizeof(float), 1, fp);
      }
      Fclose(fp);
    }
    Fclose(infp);
    
    /* forego the frees */
    exit(0);
}
