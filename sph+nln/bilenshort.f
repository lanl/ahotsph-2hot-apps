c © 2026. Triad National Security, LLC. All rights reserved.
c This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration. All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration. The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare. derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.

      subroutine bilenshort(ax,ay,aa,x,y,a,lx,ly,lz,i,j,k)
c..bilinear interpolation with derivatives
      implicit none
      integer*4 lx, ly, lz
      integer*4 i, j, k
      real*8 aa(lz,lx,ly), ax(lx), ay(ly)
      real*8 x, y, a, t, u
      real*8 dtdx, dudy, dadx, dady
c---------------------------------------------------------------------
c..interpolation factors
      dtdx = 1.0d0/( ax(i+1) - ax(i) )
      t    = (x - ax(i) )*dtdx

      dudy = 1.0d0/( ay(j+1) - ay(j) )
      u    = (y - ay(j) )*dudy

c..bilinear interpolation function
      a = (1.0d0 - t)*(1.0d0 - u)*aa(k,i  ,j)
     1  +          t *(1.0d0 - u)*aa(k,i+1,j)
     2  +          t *         u *aa(k,i+1,j+1)
     3  + (1.0d0 - t)*         u *aa(k,i  ,j+1)

c..analytic derivatives of approximation function a
c      dadx = dtdx*( (1.0d0 - u)*( -aa(  i,j  ) + aa(i+1,j  ) )
c     1             +         u *(  aa(i+1,j+1) - aa(i  ,j+1) ) )
c      dady = dudy*( (1.0d0 - t)*( -aa(  i,j  ) + aa(i  ,j+1) )
c     1             +         t *( -aa(i+1,j  ) + aa(i+1,j+1) ) )

      return

      end

