
/* Code to produce triangular panels on the sphere, and to output all the
   necessary stuff for the panel tree code     */

/* This version uses an isocahedron (20 faces, each being an equilateral triangle;
   12 vertices, each sitting on the unit sphere) as a first cut. Then each face
   of the isocahedron is broken up into 4 equilateral triangles. Each new triangle
   is itself decomposed into 4 equilateral triangles, etc... When the decomposition
   has been done up to level, l, = level_deep, d, then the small triangular 
   panels are projected onto the unit sphere.
*/

#include <stdio.h>
#include <math.h>
#include "SDF.h"
#include "getparam.h"

/* This isocahedron has its vertices on the unit sphere!
    ratio R/S is (1 + sqrt(5))/2 = 1.618033989  */
#define R 0.850650808  /* sqrt( (5 + sqrt(5))/10 ) */
#define S 0.525731112  /* sqrt( (5 - sqrt(5))/10 ) */

/* Face f, vertex i has coordinates (X(f,i), Y(f,i), Z(f,i))
   with f between 0 and 19, and i between 0 and 2 */
#define X(f,i)        icoscoord[icosptr[f][i]][0]
#define Y(f,i)        icoscoord[icosptr[f][i]][1]
#define Z(f,i)        icoscoord[icosptr[f][i]][2]

#define Output(x) fwrite(&x, sizeof(x), 1, stdout);


void CutTriangle(float pos0[3], float pos1[3], float pos2[3], int d, 
                 int l, int b[50]){

  float pos01[3], pos12[3], pos20[3];
  float x, y, z, cc;
  int lc, ident;

  l+=1;
  
  pos01[0]=.5*(pos0[0]+pos1[0]);
  pos01[1]=.5*(pos0[1]+pos1[1]);
  pos01[2]=.5*(pos0[2]+pos1[2]);

  pos12[0]=.5*(pos1[0]+pos2[0]);
  pos12[1]=.5*(pos1[1]+pos2[1]);
  pos12[2]=.5*(pos1[2]+pos2[2]);
  
  pos20[0]=.5*(pos2[0]+pos0[0]);
  pos20[1]=.5*(pos2[1]+pos0[1]);
  pos20[2]=.5*(pos2[2]+pos0[2]);

  if(d == l){          /* We are at the finest level! Do output */


/* ident of first triangle: */

      b[l]=0;

/* compute sum (lc=1,l) b[lc] * 4^(l-lc): */
      ident=b[1];
      if(l!=1){
        for(lc=2; lc<=l; lc++) ident=b[lc]+4*ident;
      }

      ident+=b[0];


      x=pos0[0];
      y=pos0[1];
      z=pos0[2];
      cc=(float)1./sqrt(x*x+y*y+z*z);
      x*=cc;
      y*=cc;
      z*=cc;
      Output( x );
      Output( y );
      Output( z );

      x=pos01[0];
      y=pos01[1];
      z=pos01[2];
      cc=(float)1./sqrt(x*x+y*y+z*z);
      x*=cc;
      y*=cc;
      z*=cc;
      Output( x );
      Output( y );
      Output( z );

      x=pos20[0];
      y=pos20[1];
      z=pos20[2];
      cc=(float)1./sqrt(x*x+y*y+z*z);
      x*=cc;
      y*=cc;
      z*=cc;
      Output( x );
      Output( y );
      Output( z );


      Output( ident );


/* Note: 
  ident of second triangle could be computed as: 
      b[l]=1;

      ident=b[1];
      if(l!=1) for(lc=2; lc<=l; lc++){
        ident=b[lc]+4*ident;
      }

      ident+=b[0];

but instead it is also  simply = ident of first triangle + 1 !!! 
*/

      ident+=1;             /* ident of second triangle */

      x=pos1[0];
      y=pos1[1];
      z=pos1[2];
      cc=(float)1./sqrt(x*x+y*y+z*z);
      x*=cc;
      y*=cc;
      z*=cc;
      Output( x );
      Output( y );
      Output( z );

      x=pos12[0];
      y=pos12[1];
      z=pos12[2];
      cc=(float)1./sqrt(x*x+y*y+z*z);
      x*=cc;
      y*=cc;
      z*=cc;
      Output( x );
      Output( y );
      Output( z );

      x=pos01[0];
      y=pos01[1];
      z=pos01[2];
      cc=(float)1./sqrt(x*x+y*y+z*z);
      x*=cc;
      y*=cc;
      z*=cc;
      Output( x );
      Output( y );
      Output( z );

      Output( ident );


      ident+=1;            /* ident of third triangle */

      x=pos2[0];
      y=pos2[1];
      z=pos2[2];
      cc=(float)1./sqrt(x*x+y*y+z*z);
      x*=cc;
      y*=cc;
      z*=cc;
      Output( x );
      Output( y );
      Output( z );

      x=pos20[0];
      y=pos20[1];
      z=pos20[2];
      cc=(float)1./sqrt(x*x+y*y+z*z);
      x*=cc;
      y*=cc;
      z*=cc;
      Output( x );
      Output( y );
      Output( z );

      x=pos12[0];
      y=pos12[1];
      z=pos12[2];
      cc=(float)1./sqrt(x*x+y*y+z*z);
      x*=cc;
      y*=cc;
      z*=cc;
      Output( x );
      Output( y );
      Output( z );


      Output( ident );


      ident+=1;            /* ident of fourth triangle */

      x=pos01[0];
      y=pos01[1];
      z=pos01[2];
      cc=(float)1./sqrt(x*x+y*y+z*z);
      x*=cc;
      y*=cc;
      z*=cc;
      Output( x );
      Output( y );
      Output( z );

      x=pos12[0];
      y=pos12[1];
      z=pos12[2];
      cc=(float)1./sqrt(x*x+y*y+z*z);
      x*=cc;
      y*=cc;
      z*=cc;
      Output( x );
      Output( y );
      Output( z );

      x=pos20[0];
      y=pos20[1];
      z=pos20[2];
      cc=(float)1./sqrt(x*x+y*y+z*z);
      x*=cc;
      y*=cc;
      z*=cc;
      Output( x );
      Output( y );
      Output( z );


      Output( ident );

  }

  else{               /* Cut the triangle into four triangles (next level) */

      b[l]=0;
      CutTriangle(pos0, pos01, pos20, d, l, b);

      b[l]=1;
      CutTriangle(pos1, pos12, pos01, d, l, b);

      b[l]=2;
      CutTriangle(pos2, pos20, pos12, d, l, b);

      b[l]=3; 
      CutTriangle(pos01, pos12, pos20, d, l, b);
 
  }

}



main(int argc, char *argv){

  float pos0[3], pos1[3], pos2[3];
  float x, y, z, cc;

  int f, l, d, dc, fd, b[50], Ntot, ident;

  double icoscoord[12][3] = {
   R, S, 0.,   -R, S, 0.,   -R, -S, 0.,   R, -S, 0.,
   0., R, S,   0., -R, S,   0., -R, -S,   0., R, -S,
   S, 0., R,   S, 0., -R,   -S, 0., -R,   -S, 0., R
   };

  int icosptr[20][3] = {
   8,4,11,     11,5,8,     8,0,4,     4,1,11,    11,2,5,
   5,3,8,      8,3,0,      11,1,2,    0,7,4,     4,7,1,
   2,6,5,      5,6,3,      9,0,3,     10,2,1,    9,3,6,
   6,2,10,     10,1,7,     7,0,9,     10,7,9,    9,6,10
   };


    Getiparam("how many times to cut each face of the isocahedron?", d);

/* compute 4^d: */
    fd=1;
    if(d!=0){
     for(dc=0; dc<d; dc++) fd*=4;
    }

    Ntot=20 * fd;
    fprintf(stderr, "number of panels per face= %d\n", fd);
    fprintf(stderr, "Ntot= %d\n", Ntot);


    printf("# SDF\n");
    printf("parameter byteorder = 0x%x;\n", SDFcpubyteorder());
    printf("int npart = %d;\n", Ntot);
    printf("int iter=0;\n");
/*    printf("float tpos=0;\n"); */ 
    printf("struct{\nfloat x1,y1,z1;\nfloat x2,y2,z2;\nfloat x3,y3,z3;\nint ident;}[%d];\n", Ntot);
    printf("# \f\n");
    printf("# SDF-EOH \n");


/* *****************************************************************************
First cut: an isocahedron:   */
 
    for(f=0; f<20; f++){

      l=0;
      b[0]=f*fd;

      pos0[0]=X(f,0);
      pos0[1]=Y(f,0);
      pos0[2]=Z(f,0);

      pos1[0]=X(f,1);
      pos1[1]=Y(f,1);
      pos1[2]=Z(f,1);

      pos2[0]=X(f,2);
      pos2[1]=Y(f,2);
      pos2[2]=Z(f,2);

      if(d == l){            /* then output the deformed isocahedron only */

      ident=b[0];

      x=pos0[0];
      y=pos0[1];
      z=pos0[2];
      cc=(float)1./sqrt(x*x+y*y+z*z);
      x*=cc;
      y*=cc;
      z*=cc;
      Output( x );
      Output( y );
      Output( z );

      x=pos1[0];
      y=pos1[1];
      z=pos1[2];
      cc=(float)1./sqrt(x*x+y*y+z*z);
      x*=cc;
      y*=cc;
      z*=cc;
      Output( x );
      Output( y );
      Output( z );

      x=pos2[0];
      y=pos2[1];
      z=pos2[2];
      cc=(float)1./sqrt(x*x+y*y+z*z);
      x*=cc;
      y*=cc;
      z*=cc;
      Output( x );
      Output( y );
      Output( z );

      Output( ident );

     }

     else{                   /* Cut each face of the isocahedron into 4 triangles */

     CutTriangle(pos0, pos1, pos2, d, l, b);

     }

   }

}


