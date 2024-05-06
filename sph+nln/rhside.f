      subroutine rhside(idebug)
c..   defines right hand side for instantaneous rate
c     and for solution
c..   specify b array elements for ax = b
c..   derivatives of rhs with respect to T9, V (wda 5/9/99)
c..   revised 1/8/00

      implicit none

      include 'dimen'
      include 'crate'
      include 'comcsolve'

      integer*4 idebug, i, k
      integer*4 k1,k2, i1,i2,i3,i4,i5,i6
      real*8    yysig, ysig, y3sig
c---------------------------------------------------------------
      do i = 1, itot
         b(i)  = 0.0d0
         bt(i) = 0.0d0
         bv(i) = 0.0d0
      enddo

c..   deck 1: i ---> j
c..   beta decays, positron decays,
c..   i + e ---> j for electron captures (see rate.f)
      weaksum = 0.0d0
      k1 = k1deck(1)
      k2 = k2deck(1)
      if( k1 .ne. 0 .and. k2 .ne. 0 )then
         do k = k1,k2
            i1 = nrr(1,k)
            i2 = nrr(2,k)
            b(i1)    = b(i1)    - y(i1)*sig(k)
            b(i2)    = b(i2)    + y(i1)*sig(k)
            weaksum = weaksum + ( nz(i2)-nz(i1) )*y(i1)*sig(k)
c..   dYe/dt is weaksum
            bt(i1)    = bt(i1)    - y(i1)*sigt(k)
            bv(i1)    = bv(i1)    - y(i1)*sigv(k)
            bt(i2)    = bt(i2)    + y(i1)*sigt(k)
            bv(i2)    = bv(i2)    + y(i1)*sigv(k)
         enddo
      endif

      if( idebug .ne. 0 )then
         do k = k1,k2
            i1 = nrr(1,k)
            i2 = nrr(2,k)
            write(*,'(2i5,a5,2x,a5,2x,a5,5x,a4,2a1,1p6e12.3)')
     1           k,ideck(k),
     1           xid(i1), " --> ", xid(i2),rlkh(k),rnr(k),rvw(k),
     1           y(i1)*sig(k),y(i1),sig(k)
         enddo
      endif
     
c----------------------------------------------------
c..   deck 2: dissociations, i ---> j + k
      k1 = k1deck(2)
      k2 = k2deck(2)
      if( k1 .ne. 0 .and. k2 .ne. 0 )then
         do k = k1,k2
            i1 = nrr(1,k)
            i2 = nrr(2,k)
            i3 = nrr(3,k)
            ysig = y(i1)*sig(k)
            b(i1)    = b(i1)    - ysig
            b(i2)    = b(i2)    + ysig
            b(i3)    = b(i3)    + ysig
            bt(i1)   = bt(i1)   - y(i1)*sigt(k)
            bv(i1)   = bv(i1)   - y(i1)*sigv(k)
            bt(i2)   = bt(i2)   + y(i1)*sigt(k)
            bv(i2)   = bv(i2)   + y(i1)*sigv(k)
            bt(i3)   = bt(i3)   + y(i1)*sigt(k)
            bv(i3)   = bv(i3)   + y(i1)*sigv(k)
         enddo
      endif

      if( idebug .ne. 0 )then
         do k = k1,k2
            i1 = nrr(1,k)
            i2 = nrr(2,k)
            i3 = nrr(3,k)
            ysig = y(i1)*sig(k)
            
            write(*,'(i4,i2,3(a5,1x),a5,2x,a4,2a1,1p7e10.2)')
     1           k,ideck(k),
     1           xid(i1), " --> ", xid(i2), xid(i3),rlkh(k),rnr(k),
     1           rvw(k),ysig,y(i1),sig(k)
         enddo
      endif
     
c----------------------------------------------------
c..   deck 3: dissociations, i ---> j + k + l
      k1 = k1deck(3)
      k2 = k2deck(3)
      if( k1 .ne. 0 .and. k2 .ne. 0 )then
         do k = k1,k2
            i1 = nrr(1,k)
            i2 = nrr(2,k)
            i3 = nrr(3,k)
            i4 = nrr(4,k)
            ysig = y(i1)*sig(k)
            b(i1)    = b(i1)    - ysig
            b(i2)    = b(i2)    + ysig
            b(i3)    = b(i3)    + ysig
            b(i4)    = b(i4)    + ysig

            bt(i1)    = bt(i1)    - y(i1)*sigt(k)
            bv(i1)    = bv(i1)    - y(i1)*sigv(k)
            bt(i2)    = bt(i2)    + y(i1)*sigt(k)
            bv(i2)    = bv(i2)    + y(i1)*sigv(k)
            bt(i3)    = bt(i3)    + y(i1)*sigt(k)
            bv(i3)    = bv(i3)    + y(i1)*sigv(k)
            bt(i4)    = bt(i4)    + y(i1)*sigt(k)
            bv(i4)    = bv(i4)    + y(i1)*sigv(k)
         enddo
      endif

      if( idebug .ne. 0 )then
         do k = k1,k2
            i1 = nrr(1,k)
            i2 = nrr(2,k)
            i3 = nrr(3,k)
            i4 = nrr(4,k)
            ysig = y(i1)*sig(k)

            write(*,'(i4,i2,4(a5,1x),a5,2x,a4,2a1,1p7e10.2)')
     1           k,ideck(k),
     1           xid(i1), " --> ", xid(i2), xid(i3), xid(i4),
     1           rlkh(k),rnr(k),rvw(k),
     1           ysig,y(i1),sig(k)
         enddo
      endif
      
c----------------------------------------------------
c..   deck 4: captures, i + j ---> k
      k1 = k1deck(4)
      k2 = k2deck(4)
      if( k1 .ne. 0 .and. k2 .ne. 0 )then
         do k = k1,k2
            i1 = nrr(1,k)
            i2 = nrr(2,k)
            i3 = nrr(3,k)
            yysig = y(i1)*y(i2)*sig(k)
            b(i1)    = b(i1)    - yysig
            b(i2)    = b(i2)    - yysig
            b(i3)    = b(i3)    + yysig

            bt(i1)    = bt(i1)    - y(i1)*y(i2)*sigt(k)
            bv(i1)    = bv(i1)    - y(i1)*y(i2)*sigv(k)
            bt(i2)    = bt(i2)    - y(i1)*y(i2)*sigt(k)
            bv(i2)    = bv(i2)    - y(i1)*y(i2)*sigv(k)
            bt(i3)    = bt(i3)    + y(i1)*y(i2)*sigt(k)
            bv(i3)    = bv(i3)    + y(i1)*y(i2)*sigv(k)
         enddo
      endif

      if( idebug .ne. 0 )then
         do k = k1,k2
            i1 = nrr(1,k)
            i2 = nrr(2,k)
            i3 = nrr(3,k)
            yysig = y(i1)*y(i2)*sig(k)

            write(*,'(i4,i2,3(a5,1x),a5,2x,a4,2a1,1p7e10.2)')
     1           k,ideck(k),
     1           xid(i1), xid(i2), " --> ", xid(i3),rlkh(k),rnr(k),
     1           rvw(k),yysig,y(i1),y(i2),sig(k)
         enddo
      endif
      
c----------------------------------------------------
c..   deck 5: exchange, i + j ---> k + l
      k1 = k1deck(5)
      k2 = k2deck(5)
      if( k1 .ne. 0 .and. k2 .ne. 0 )then
         do k = k1,k2
            i1 = nrr(1,k)
            i2 = nrr(2,k)
            i3 = nrr(3,k)
            i4 = nrr(4,k)
            yysig = y(i1)*y(i2)*sig(k)
            b(i1)    = b(i1)    - yysig
            b(i2)    = b(i2)    - yysig
            b(i3)    = b(i3)    + yysig
            b(i4)    = b(i4)    + yysig

            bt(i1)    = bt(i1)    - y(i1)*y(i2)*sigt(k)
            bv(i1)    = bv(i1)    - y(i1)*y(i2)*sigv(k)
            bt(i2)    = bt(i2)    - y(i1)*y(i2)*sigt(k)
            bv(i2)    = bv(i2)    - y(i1)*y(i2)*sigv(k)
            bt(i3)    = bt(i3)    + y(i1)*y(i2)*sigt(k)
            bv(i3)    = bv(i3)    + y(i1)*y(i2)*sigv(k)
            bt(i4)    = bt(i4)    + y(i1)*y(i2)*sigt(k)
            bv(i4)    = bv(i4)    + y(i1)*y(i2)*sigv(k)
         enddo
      endif

      if( idebug .ne. 0 )then
         do k = k1,k2
            i1 = nrr(1,k)
            i2 = nrr(2,k)
            i3 = nrr(3,k)
            i4 = nrr(4,k)
            yysig = y(i1)*y(i2)*sig(k)
           
            write(*,'(i4,i2,4(a5,1x),a5,2x,a4,2a1,1p5e10.2)')
     1           k,ideck(k),
     1           xid(i1), xid(i2), " --> ", xid(i3),xid(i4),
     1           rlkh(k),rnr(k),rvw(k),
     1           yysig,y(i1),y(i2),sig(k)
         enddo
      endif
      
c--------------------------------------------------------------
c..   deck 6: exchange, i + j ---> k + l + m
      k1 = k1deck(6)
      k2 = k2deck(6)
      if( k1 .ne. 0 .and. k2 .ne. 0 )then
         do k = k1,k2
            i1 = nrr(1,k)
            i2 = nrr(2,k)
            i3 = nrr(3,k)
            i4 = nrr(4,k)
            i5 = nrr(5,k)
            yysig = y(i1)*y(i2)*sig(k)
            b(i1)    = b(i1)    - yysig
            b(i2)    = b(i2)    - yysig
            b(i3)    = b(i3)    + yysig
            b(i4)    = b(i4)    + yysig
            b(i5)    = b(i5)    + yysig

            bt(i1)    = bt(i1)    - y(i1)*y(i2)*sigt(k)
            bt(i2)    = bt(i2)    - y(i1)*y(i2)*sigt(k)
            bt(i3)    = bt(i3)    + y(i1)*y(i2)*sigt(k)
            bt(i4)    = bt(i4)    + y(i1)*y(i2)*sigt(k)
            bt(i5)    = bt(i5)    + y(i1)*y(i2)*sigt(k)
            bv(i1)    = bv(i1)    - y(i1)*y(i2)*sigv(k)
            bv(i2)    = bv(i2)    - y(i1)*y(i2)*sigv(k)
            bv(i3)    = bv(i3)    + y(i1)*y(i2)*sigv(k)
            bv(i4)    = bv(i4)    + y(i1)*y(i2)*sigv(k)
            bv(i5)    = bv(i5)    + y(i1)*y(i2)*sigv(k)
         enddo
      endif

      if( idebug .ne. 0 )then
         do k = k1,k2
            i1 = nrr(1,k)
            i2 = nrr(2,k)
            i3 = nrr(3,k)
            i4 = nrr(4,k)
            i5 = nrr(5,k)
            yysig = y(i1)*y(i2)*sig(k)

            write(*,'(i4,i2,5(a5,1x),a5,2x,a4,2a1,1p6e10.2)')
     1           k,ideck(k),
     1           xid(i1), xid(i2), " --> ", xid(i3),xid(i4),xid(i5),
     1           rlkh(k),rnr(k),rvw(k),
     1           yysig,y(i1),y(i2),sig(k)
         enddo
      endif
      
c-------------------------------------------------------------------
c..   deck 7: exchange, i + j ---> k + l + m + n
      k1 = k1deck(7)
      k2 = k2deck(7)
      if( k1 .ne. 0 .and. k2 .ne. 0 )then
         do k = k1,k2
            i1 = nrr(1,k)
            i2 = nrr(2,k)
            i3 = nrr(3,k)
            i4 = nrr(4,k)
            i5 = nrr(5,k)
            i6 = nrr(6,k)
            yysig = y(i1)*y(i2)*sig(k)
            b(i1)    = b(i1)    - yysig
            b(i2)    = b(i2)    - yysig
            b(i3)    = b(i3)    + yysig
            b(i4)    = b(i4)    + yysig
            b(i5)    = b(i5)    + yysig
            b(i6)    = b(i6)    + yysig


            bt(i1)    = bt(i1)    - y(i1)*y(i2)*sigt(k)
            bt(i2)    = bt(i2)    - y(i1)*y(i2)*sigt(k)
            bt(i3)    = bt(i3)    + y(i1)*y(i2)*sigt(k)
            bt(i4)    = bt(i4)    + y(i1)*y(i2)*sigt(k)
            bt(i5)    = bt(i5)    + y(i1)*y(i2)*sigt(k)
            bt(i6)    = bt(i6)    + y(i1)*y(i2)*sigt(k)
            bv(i1)    = bv(i1)    - y(i1)*y(i2)*sigv(k)
            bv(i2)    = bv(i2)    - y(i1)*y(i2)*sigv(k)
            bv(i3)    = bv(i3)    + y(i1)*y(i2)*sigv(k)
            bv(i4)    = bv(i4)    + y(i1)*y(i2)*sigv(k)
            bv(i5)    = bv(i5)    + y(i1)*y(i2)*sigv(k)
            bv(i6)    = bv(i6)    + y(i1)*y(i2)*sigv(k)
         enddo
      endif

      if( idebug .ne. 0 )then
         do k = k1,k2
            i1 = nrr(1,k)
            i2 = nrr(2,k)
            i3 = nrr(3,k)
            i4 = nrr(4,k)
            i5 = nrr(5,k)
            i6 = nrr(6,k)
            yysig = y(i1)*y(i2)*sig(k)

            write(*,'(i4,i2,6(a5,1x),a5,2x,a4,2a1,1p7e9.2)')
     1           k,ideck(k),
     1           xid(i1), xid(i2), " --> ", xid(i3),xid(i4),xid(i5),
     1           xid(i6),rlkh(k),rnr(k),rvw(k),
     1           yysig,y(i1),y(i2),sig(k)
         enddo
      endif
      
c----------------------------------------------------------
c..   deck 8: exchange, i + j + k ---> l
c..   and  i + j + k ---> l + m
      k1 = k1deck(8)
      k2 = k2deck(8)
      if( k1 .ne. 0 .and. k2 .ne. 0 )then
         do k = k1,k2
            i1 = nrr(1,k)
            i2 = nrr(2,k)
            i3 = nrr(3,k)
            i4 = nrr(4,k)
            i5 = nrr(5,k)

c..   triple alpha checked
               y3sig = y(i1)*y(i2)*y(i3)*sig(k)
               b(i1)    = b(i1)    - y3sig
               b(i2)    = b(i2)    - y3sig
               b(i3)    = b(i3)    - y3sig
               b(i4)    = b(i4)    + y3sig
             
               bt(i1)    = bt(i1)    - y(i1)*y(i2)*y(i3)*sigt(k)
               bt(i2)    = bt(i2)    - y(i1)*y(i2)*y(i3)*sigt(k)
               bt(i3)    = bt(i3)    - y(i1)*y(i2)*y(i3)*sigt(k)
               bt(i4)    = bt(i4)    + y(i1)*y(i2)*y(i3)*sigt(k)

               bv(i1)    = bv(i1)    - y(i1)*y(i2)*y(i3)*sigv(k)
               bv(i2)    = bv(i2)    - y(i1)*y(i2)*y(i3)*sigv(k)
               bv(i3)    = bv(i3)    - y(i1)*y(i2)*y(i3)*sigv(k)
               bv(i4)    = bv(i4)    + y(i1)*y(i2)*y(i3)*sigv(k)
            if( i5 .ne. 0 )then
               b(i5)     = b(i5)     + y(i1)*y(i2)*y(i3)*sig(k)
               bt(i5)    = bt(i5)    + y(i1)*y(i2)*y(i3)*sigt(k)
               bv(i5)    = bv(i5)    + y(i1)*y(i2)*y(i3)*sigv(k)
            endif
         enddo
      endif
      
c----------------------------------------------------
      if( idebug .ne. 0 )then
         do k = k1,k2
            i1 = nrr(1,k)
            i2 = nrr(2,k)
            i3 = nrr(3,k)
            i4 = nrr(4,k)
            i5 = nrr(5,k)
            if( i5 .eq. 0 )then
               y3sig = y(i1)*y(i2)*y(i3)*sig(k)

               write(*,'(i4,i2,4(a5,1x),a5,2x,a4,2a1,1p5e10.2)')
     1              k,ideck(k),
     1              xid(i1), xid(i2), xid(i3), " --> ",xid(i4),
     1              rlkh(k),rnr(k),rvw(k),
     1              y3sig,y(i1),y(i2),y(i3),sig(k)
            else
               y3sig = y(i1)*y(i2)*y(i3)*sig(k)

               write(*,'(i4,i2,5(a5,1x),a5,1x,a4,2a1,1p6e9.2)')
     1              k,ideck(k),
     1              xid(i1), xid(i2), xid(i3), " --> ",xid(i4),xid(i5), 
     1              rlkh(k),rnr(k),rvw(k),
     1              y3sig,y(i1),y(i2),y(i3),sig(k)
            endif
         enddo
c..   final idebug command
c         stop'rhside'
      endif

      return
      end


