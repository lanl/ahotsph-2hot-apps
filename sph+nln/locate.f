c © 2026. Triad National Security, LLC. All rights reserved.
c This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration. All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration. The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare. derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.

      subroutine locate(xx,n,x,j)
      implicit none
c..   numerical receipes bisection algorithm
c..   locate index j for variable x lying between j, j+1
c..   in array xx of dimension n
c..   j  = desired index
c..   xx = array of tabulated values
c..   n  = dimension of xx
c..   x  = variable value to be interpolated

      integer*4 j,n
      real*8    x, xx(n)
      integer*4 jl,jm,ju
c---------------------------------------------------------------------
      jl = 0
      ju = n+1
 10   if( ju-jl .gt. 1 )then
         jm = (ju+jl)/2
         if( (xx(n) .gt. xx(1)) .eqv. (x .gt. xx(jm)) )then
            jl = jm
         else
            ju = jm
         endif
         goto 10
      endif

      j = jl
      
      return
      end
