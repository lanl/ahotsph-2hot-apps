c © 2026. Triad National Security, LLC. All rights reserved.
c This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration. All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration. The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare. derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.

      subroutine leqs (a,b,n,idim)
      implicit none

      include 'dimen'

c..   external
      real*8 a,b
      integer*4 n, idim

c..   internal
      integer*4 n1,i,j,k,l,jj,imax
      real*8 r,c,zero

      dimension  a(ndim,ndim),b(ndim)

      parameter( zero = 0.0d0 )

c-----------------------------------------------------
c.....leqs performs matrix inversion
c     n-dimension of a

      n1 = n-1

c..   Find maximum element in each row, and divide

      do 19 i = 1, n
         r = dabs( a(i,1) )
         do 16 j = 2, n
            c = dabs( a(i,j) )
            r = dmax1( r, c )
 16      continue
         do 17 j = 1, n
            a(i,j) = a(i,j) / r
 17      continue
         b(i) = b(i) / r
 19   continue

      do 9 j = 1, n1
         l = j + 1
c..   no pivoting
         do 9 i = l, n
            r = -a(i,j)
            if ( r .ne. zero )then
               r = r / a(j,j)
               do 7 k = l, n
                  a(i,k) = a(i,k) + r*a(j,k)
 7             continue
               b(i) = b(i) + r*b(j)
            endif
 9       continue

c..   The matrix is now in triangular form
c..   start the back substitution and find solution

         b(n) = b(n) / a(n,n)
         do 13 l = 1, n1
            i = n-l
            r = zero
            imax = i+1
            do 12 j = imax, n
               jj  = i + n+1 - j
               r   = r + a(i,jj)*b(jj)
 12         continue
            b(i) = (b(i)-r)/a(i,i)
 13      continue

         return
         end

