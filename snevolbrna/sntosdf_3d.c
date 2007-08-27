#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include "bigmalloc.h"
#include "macr.h"
#include "error.h"

static unsigned char junk[4] = {0x12, 0x34, 0x56, 0x78};
static unsigned int *cpubyteorder = (unsigned int *)&junk[0];

void
main(int argc, char *argv[])
{
    FILE *infp;
    FILE *fp;
    char *in_name, *out_name;
    int npart;
    float t, gamma, tkin, tterm;
    float *h, *x, *y, *z, *vx, *vy, *vz, *u, *mass, *rho;
    float *abar, *temp, *ye, *xp, *xn;
    int *ifleos;
    float *ynue, *ynueb, *ynux, *unue, *unueb, *unux, *ufreez;
    float *pr, *u2, *te, *teb, *tx;
    float xmcore, rb, ftrape, ftrapb, ftrapx;
    float *kepcel, *xpf, *pvar2, *pvar3, *pvar4;
    int fortran_crap;
    int i;
    char msgfile[256];

    if (argc == 3) {
	in_name = argv[1];
	out_name = argv[2];
    } else {
	fprintf(stderr, "Usage: %s infile outfile\n", argv[0]);
	exit(1);
    }

    Fopen(infp, in_name, "r");
    Fopen(fp, out_name, "w");

    Fread(&fortran_crap, sizeof(int), 1, infp);
    Fread(&npart, sizeof(int), 1, infp);
    if (npart < 0 || npart > 1000000) {
	fprintf(stderr, "Npart is %d, suspect byte order problem\n", npart);
	exit(1);
    }
    Fread(&t, sizeof(float), 1, infp);
    Fread(&gamma, sizeof(float), 1, infp);
    fprintf(stderr, "npart %d, t %f, gamma %f\n",
	    npart, t, gamma);

    h = malloc(npart * sizeof(float));
    x = malloc(npart * sizeof(float));
    y = malloc(npart * sizeof(float));
    z = malloc(npart * sizeof(float));
    vx = malloc(npart * sizeof(float));
    vy = malloc(npart * sizeof(float));
    vz = malloc(npart * sizeof(float));
    u = malloc(npart * sizeof(float));
    mass = malloc(npart * sizeof(float));
    abar = malloc(npart * sizeof(float));
    rho = malloc(npart * sizeof(float));
    temp = malloc(npart * sizeof(float));
    ye = malloc(npart * sizeof(float));
    xp = malloc(npart * sizeof(float));
    xn = malloc(npart * sizeof(float));
    ifleos = malloc(npart * sizeof(int));
    ynue = malloc(npart * sizeof(float));
    ynueb = malloc(npart * sizeof(float));
    ynux = malloc(npart * sizeof(float));
    unue = malloc(npart * sizeof(float));
    unueb = malloc(npart * sizeof(float));
    unux = malloc(npart * sizeof(float));
    ufreez = malloc(npart * sizeof(float));
    pr = malloc(npart * sizeof(float));
    u2 = malloc(npart * sizeof(float));
    te = malloc(npart * sizeof(float));
    teb = malloc(npart * sizeof(float));
    tx = malloc(npart * sizeof(float));
    kepcel = malloc(npart * sizeof(float));
    xpf = malloc(npart * sizeof(float));
    pvar2 = malloc(npart * sizeof(float));
    pvar3 = malloc(npart * sizeof(float));
    pvar4 = malloc(npart * sizeof(float));

    Fread(h, sizeof(float), npart, infp);
    Fread(&tkin, sizeof(float), 1, infp);
    Fread(&tterm, sizeof(float), 1, infp);
    Fread(x, sizeof(float), npart, infp);
    Fread(y, sizeof(float), npart, infp);
    Fread(z, sizeof(float), npart, infp);
    Fread(vx, sizeof(float), npart, infp);
    Fread(vy, sizeof(float), npart, infp);
    Fread(vz, sizeof(float), npart, infp);
    Fread(u, sizeof(float), npart, infp);
    Fread(mass, sizeof(float), npart, infp);
    Fread(abar, sizeof(float), npart, infp);
    Fread(rho, sizeof(float), npart, infp);
    Fread(temp, sizeof(float), npart, infp);
    Fread(ye, sizeof(float), npart, infp);
    Fread(xp, sizeof(float), npart, infp);
    Fread(xn, sizeof(float), npart, infp);
    Fread(ifleos, sizeof(int), npart, infp);
    Fread(ynue, sizeof(float), npart, infp);
    Fread(ynueb, sizeof(float), npart, infp);
    Fread(ynux, sizeof(float), npart, infp);
    Fread(unue, sizeof(float), npart, infp);
    Fread(unueb, sizeof(float), npart, infp);
    Fread(unux, sizeof(float), npart, infp);
    Fread(ufreez, sizeof(float), npart, infp);
    Fread(pr, sizeof(float), npart, infp);
    Fread(u2, sizeof(float), npart, infp);
    Fread(te, sizeof(float), npart, infp);
    Fread(teb, sizeof(float), npart, infp);
    Fread(tx, sizeof(float), npart, infp);
    Fread(&xmcore, sizeof(float), 1, infp);
    Fread(&rb, sizeof(float), 1, infp);
    Fread(&ftrape, sizeof(float), 1, infp);
    Fread(&ftrapb, sizeof(float), 1, infp);
    Fread(&ftrapx, sizeof(float), 1, infp);
    Fread(kepcel, sizeof(float), 1, infp);
    Fread(xpf, sizeof(float), 1, infp);
    Fread(pvar2, sizeof(float), 1, infp);
    Fread(pvar3, sizeof(float), 1, infp);
    Fread(pvar4, sizeof(float), 1, infp);

    Fclose(infp);

    fprintf(fp, "# SDF\n");
    fprintf(fp, "char name[] = \"sntosdf on file %s\";\n", in_name);
    fprintf(fp, "int npart = %d;\n", npart);
    fprintf(fp, "float tpos = %f;\n", t);
    fprintf(fp, "float Gamma = %f;\n", gamma);
    fprintf(fp, "float tkin = %f;\n", tkin);
    fprintf(fp, "float tterm = %f;\n", tterm);
    fprintf(fp, "float xmcore = %f;\n", xmcore);
    fprintf(fp, "float rb = %f;\n", rb);
    fprintf(fp, "float ftrape = %f;\n", ftrape);
    fprintf(fp, "float ftrapb = %f;\n", ftrapb);
    fprintf(fp, "float ftrapx = %f;\n", ftrapx);
    fprintf(fp, "int ndim = %d;\n", 3);
    fprintf(fp, "parameter byteorder = 0x%x;\n", *cpubyteorder);
    fputs  ("struct {\n", fp);
    fputs  ("\tfloat mass;\n", fp);
    fputs  ("\tfloat x, y, z;\n", fp);
    fputs  ("\tfloat vx, vy, vz;\n", fp);
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
    fputs  ("\tfloat te;\n", fp);
    fputs  ("\tfloat teb;\n", fp);
    fputs  ("\tfloat tx;\n", fp);
    fputs  ("}", fp);
    fprintf(fp, "[%d];\n", npart);
    fputs  ("#\f\n", fp);
    fputs  ("# SDF-EOH\n", fp);

    for (i = 0; i < npart; i++) {
	Fwrite(mass+i, sizeof(float), 1, fp);
	Fwrite(x+i, sizeof(float), 1, fp);
	Fwrite(y+i, sizeof(float), 1, fp);
	Fwrite(z+i, sizeof(float), 1, fp);
	Fwrite(vx+i, sizeof(float), 1, fp);
	Fwrite(vy+i, sizeof(float), 1, fp);
	Fwrite(vz+i, sizeof(float), 1, fp);
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
	Fwrite(te+i, sizeof(float), 1, fp);
	Fwrite(teb+i, sizeof(float), 1, fp);
	Fwrite(tx+i, sizeof(float), 1, fp);
    }
    Fclose(fp);
    
    /* forego the frees */
    exit(0);
}
