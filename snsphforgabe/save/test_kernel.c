#include <stdio.h>

/* External Fortran linkage */
#define Fortran(x) x##_
/* GNU Fortran adds two underscores if there is an underscore in the name */
#ifdef __GNUC__
#define Fortran2(x) x##__
#else
#define Fortran2(x) x##_
#endif


#define ITABLE 80000
extern void Fortran(ktable)(void);
/* common /table/ wij(0:itable), grwij(0:itable), dvtable */
extern float Fortran(table)[ITABLE*2];

main(int argc, char *argv)
{
  int i;
  float *wij, *grwij;

  Fortran(ktable)();
  wij = Fortran(table);
  grwij = Fortran(table)+ITABLE+1;
  for (i = 0; i < ITABLE; i++) {
    printf("%5d %8g %8g\n", i, wij[i], grwij[i]);
  }
}
