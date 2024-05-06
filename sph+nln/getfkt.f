      subroutine getfkt(ibug,irank)

c..reads Thielemann file format
c..identifies symbol in Z,N
c..and echo to standard i/o
c---------------------------------------------------------------
c..to do: read netwinv to get symbols, Z, N, q, w and part. func.
c..and avoid a lot of this searching
c---------------------------------------------------------------

      implicit real*8 (a-h,o-z)

      include 'dimen'

      integer*4 irank

c..these values apply to fkt reaclib 1988
c      parameter(niso = 806)
c      parameter( nzmax=36,nnmax=73)

c..these values apply to fkt reaclib.nosmo
      parameter(niso = 5410)
      parameter( nzmax=85,nnmax=168)

      dimension p(7),wkp(7)
      character*4 lkh,wklkh
      character*5 nam(6),inam(niso),blank,rnam(6),wknam(6)
      character*1 vw,nr,wkvw,wknr

      dimension isoz(niso), ison(niso)
      dimension izr(6),inr(6),id(6)
      dimension nnn(nzmax,nnmax)

      logical tobe

      include 'crate'
      include 'comcsolve'

c===============================================================
      data blank/'     '/
c............................................
      do j=1,nnmax
        do i=1,nzmax
          nnn(i,j) = 0
        enddo
      enddo
c............................................
c..get isotope names, Z and N
      inquire(file='isotope.lib',exist=tobe)
      if( .not. tobe )then
         if(irank.eq.0) write(*,*)'No file isotope.lib found'
         stop'getfkt.f'
      endif
      open(12,file='isotope.lib')
      i = 0
10    continue
      i = i + 1
      read(12,'(a5,2i4)',end=11)inam(i),isoz(i),ison(i)
c      write(*,'(a5,2i4)')inam(i),isoz(i),ison(i)
      goto 10
11    continue
      i = i-1
      if(irank.eq.0) write(*,*)'isotope.lib read,',i,' isotopes'
c      write(3,*)'isotope.lib read,',i,' isotopes'
      close(12)

      
c..read thielemann coef.s

      open(5,file='netsu')
c..friedel-style format
300   FORMAT(i1,4x,6a5,8x,a4,a1,a1,3x,1pe12.5,i5)
310   FORMAT(4e13.6)

      kk = 0
      ktype = 0
20    continue
      READ(5,300,END=90) k,(nam(j),j=1,6),lkh,nr,vw,qqf,iffnu
      READ(5,310) (p(j),j=1,7)

c      write(*,300) k,(nam(j),j=1,6),lkh,nr,vw,qqf,iffnu
cccccccccccccccccccccccccc

      if( k .ne. 0 )then
c..flag for type of reaction (collected in "decks")
        ktype = ktype + 1
c        write(*,*)'deck is ',ktype
        do j =1,6
         nam(j) = blank
        enddo
      else
c..store values in arrays for rate subroutine
        kk = kk+1
        do j = 1,6
          rname(j,kk) = nam(j)
        enddo
        do j = 1,7
          rcoef(j,kk) = p(j)
c          if(ktype .eq. 1)then
c             write(*,*)j,rcoef(j,kk)
c          endif
        enddo
        qval(kk) = qqf
        rlkh(kk) = lkh
        rnr(kk)  = nr
        rvw(kk)  = vw
        ideck(kk) = ktype
        iffn(kk) = iffnu
c..identify all participating nuclei
        call find(nam(1),id(1),inam,blank)
        call find(nam(2),id(2),inam,blank)
        call find(nam(3),id(3),inam,blank)
        call find(nam(4),id(4),inam,blank)
        call find(nam(5),id(5),inam,blank)
        call find(nam(6),id(6),inam,blank)
      endif

c      write(*,*)'before decks '

c--deck 1--------------------------------------------
      if( ktype .eq. 1 .and. nam(1) .ne. blank )then
         if( id(1) .gt. 0 .and. id(1) .le. niso )then
           rnam(1) = inam(id(1))
           izr(1)  = isoz(id(1))
           inr(1)  = ison(id(1))
         else
           rnam(1) = blank
           izr(1) = 0
           inr(1) = 0
         endif
         if( id(2) .gt. 0 .and. id(2) .le. niso )then
           rnam(2) = inam(id(2))
           izr(2) = isoz(id(2))
           inr(2) = ison(id(2))
         else
           rnam(2) = blank
           izr(2) = 0
           inr(2) = 0
         endif
         ncharge = izr(2) - izr(1)
         nnucleon = ncharge +  inr(2) - inr(1)
         if( ibug .ne. 0 )then
          if(irank.eq.0) write(*,
     1        '(a5,3i5,10x,a5,3i5,10x,2i5,2x,a4,2a1)')
     1        rnam(1),izr(1),inr(1),id(1),
     1        rnam(2),izr(2),inr(2),id(2),
     1        ncharge,nnucleon,lkh,nr,vw
          write(3,'(a5,3i5,10x,a5,3i5,10x,2i5,2x,a4,2a1)')
     1        rnam(1),izr(1),inr(1),id(1),
     1        rnam(2),izr(2),inr(2),id(2),
     1        ncharge,nnucleon,lkh,nr,vw
        endif
c..check reaction identities
         if( nnucleon .ne. 0 )then
           if(irank.eq.0) write(*,*)'nucleon nonconservation',ktype,kk
           stop
         endif
         if( ncharge .eq. 0 )then
           if(irank.eq.0) write(*,*)'charge nonconservation',ktype,kk
           stop
         endif
         if( ncharge .eq. 1 )then
           if( lkh .ne. 'bet-' )then
            if( lkh .ne. ' ffn' )then
             if( lkh .ne. 'bkmo' )then
              if( lkh .ne. 'mo92' )then
              if(irank.eq.0) write(*,*)'error in beta decay'
              if(irank.eq.0) write(*,
     1            '(a5,3i5,10x,a5,3i5,10x,2i5,2x,a4,2a1)')
     1            rnam(1),izr(1),inr(1),id(1),
     1            rnam(2),izr(2),inr(2),id(2),
     1            ncharge,nnucleon,lkh,nr,vw
             
c             stop
              endif
             endif
            endif
           endif
         endif
         if( ncharge .eq. -1 )then
           if( lkh .ne. 'bet+')then
             if( lkh .ne. ' bec')then
               if( lkh .ne. '  ec')then
                if( lkh .ne. ' ffn')then
                 if( lkh .ne. 'bkmo' )then
                  if( lkh .ne. 'btyk' )then
                  if(irank.eq.0) write(*,*)'error in positron decay'
                  if(irank.eq.0) write(*,
     1            '(a5,3i5,10x,a5,3i5,10x,2i5,2x,a4,2a1)')
     1            rnam(1),izr(1),inr(1),id(1),
     1            rnam(2),izr(2),inr(2),id(2),
     1            ncharge,nnucleon,lkh,nr,vw
c                 stop
                  endif
                 endif
                endif
               endif
             endif
           endif
         endif
      endif

c--deck 2--------------------------------------------
      if( ktype .eq. 2 .and. nam(1) .ne. blank )then
         if( id(1) .gt. 0 .and. id(1) .le. niso )then
           rnam(1) = inam(id(1))
           izr(1) = isoz(id(1))
           inr(1) = ison(id(1))
         else
           rnam(1) = blank
           izr(1) = 0
           inr(1) = 0
         endif
         if( id(2) .gt. 0 .and. id(2) .le. niso )then
           rnam(2) = inam(id(2))
           izr(2) = isoz(id(2))
           inr(2) = ison(id(2))
         else
           rnam(2) = blank
           izr(2) = 0
           inr(2) = 0
         endif
         if( id(3) .gt. 0 .and. id(3) .le. niso )then
           rnam(3) = inam(id(3))
           izr(3) = isoz(id(3))
           inr(3) = ison(id(3))
         else
           rnam(3) = blank
           izr(3) = 0
           inr(3) = 0
         endif
         ncharge = izr(3) + izr(2) - izr(1)
         nnucleon = ncharge + inr(3) + inr(2) - inr(1)
         if( ibug .ne. 0 )then
          if(irank.eq.0) write(*,'(3(a5,3i5,5x),2i5,2x,a4,2a1)')
     1        rnam(1),izr(1),inr(1),id(1),
     1        rnam(2),izr(2),inr(2),id(2),
     1        rnam(3),izr(3),inr(3),id(3),
     1        ncharge,nnucleon,lkh,nr,vw
          write(3,'(3(a5,3i5,5x),2i5,2x,a4,2a1)')
     1        rnam(1),izr(1),inr(1),id(1),
     1        rnam(2),izr(2),inr(2),id(2),
     1        rnam(3),izr(3),inr(3),id(3),
     1        ncharge,nnucleon,lkh,nr,vw
         endif
c..check reaction identities
         if( nnucleon .ne. 0 )then
           if(irank.eq.0) write(*,*)'nucleon nonconservation',ktype
           stop
         endif
         if( ncharge .gt. 1 .or. ncharge .lt. -1 )then
           if(irank.eq.0) write(*,*)
     1       '2error in beta decay, ncharge=',ncharge
         endif
         if( ncharge .eq. 1 )then
           if( lkh .ne. 'bet-' )then
            if( lkh .ne. ' ffn')then
             if( lkh .ne. 'bkmo' )then
              if( lkh .ne. 'mo92' )then
               if(irank.eq.0) write(*,*)'2error in beta decay'
               if(irank.eq.0) write(*,
     1            '(a5,3i5,10x,a5,3i5,10x,2i5,2x,a4,2a1)')
     1            rnam(1),izr(1),inr(1),id(1),
     1            rnam(2),izr(2),inr(2),id(2),
     1            ncharge,nnucleon,lkh,nr,vw
c             stop
              endif
             endif
            endif
           endif
         endif
         if( ncharge .eq. -1 )then
           if( lkh .ne. 'bet+')then
             if( lkh .ne. ' bec')then
               if( lkh .ne. '  ec')then
                if( lkh .ne. ' ffn')then
                 if( lkh .ne. 'bkmo' )then
                  if( lkh .ne. 'btyk' )then
                 if(irank.eq.0) write(*,*)'2error in positron decay'
c                 stop
                  endif
                 endif
                endif
               endif
             endif
           endif
         endif
      endif
c
c--deck 3--------------------------------------------
      if( ktype .eq. 3 .and. nam(1) .ne. blank )then
         if( id(1) .gt. 0 .and. id(1) .le. niso )then
           rnam(1) = inam(id(1))
           izr(1) = isoz(id(1))
           inr(1) = ison(id(1))
         else
           rnam(1) = blank
           izr(1) = 0
           inr(1) = 0
         endif
         if( id(2) .gt. 0 .and. id(2) .le. niso )then
           rnam(2) = inam(id(2))
           izr(2) = isoz(id(2))
           inr(2) = ison(id(2))
         else
           rnam(2) = blank
           izr(2) = 0
           inr(2) = 0
         endif
         if( id(3) .gt. 0 .and. id(3) .le. niso )then
           rnam(3) = inam(id(3))
           izr(3) = isoz(id(3))
           inr(3) = ison(id(3))
         else
           rnam(3) = blank
           izr(3) = 0
           inr(3) = 0
         endif
         if( id(4) .gt. 0 .and. id(4) .le. niso )then
           rnam(4) = inam(id(4))
           izr(4) = isoz(id(4))
           inr(4) = ison(id(4))
         else
           rnam(4) = blank
           izr(4) = 0
           inr(4) = 0
         endif
         ncharge  = izr(4) + izr(3) + izr(2) - izr(1)
         nnucleon = ncharge
     1            + inr(4) + inr(3) + inr(2) - inr(1)
         if( ibug .ne. 0 )then
          if(irank.eq.0) write(*,'(4(a5,2i3,i4,3x),2i5,2x,a4,2a1)')
     1        rnam(1),izr(1),inr(1),id(1),
     1        rnam(2),izr(2),inr(2),id(2),
     1        rnam(3),izr(3),inr(3),id(3),
     1        rnam(4),izr(4),inr(4),id(4),
     1        ncharge,nnucleon,lkh,nr,vw
          write(3,'(4(a5,2i3,i4,3x),2i5,2x,a4,2a1)')
     1        rnam(1),izr(1),inr(1),id(1),
     1        rnam(2),izr(2),inr(2),id(2),
     1        rnam(3),izr(3),inr(3),id(3),
     1        rnam(4),izr(4),inr(4),id(4),
     1        ncharge,nnucleon,lkh,nr,vw
         endif
c..check reaction identities
         if( nnucleon .ne. 0 )then
           if(irank.eq.0) write(*,*)'nucleon nonconservation',ktype
           stop
         endif
         if( ncharge .gt. 1 .or. ncharge .lt. -1 )then
           if(irank.eq.0) write(*,*)
     1        '3error in beta decay, ncharge=',ncharge
         endif
         if( ncharge .eq. 1 )then
           if( lkh .ne. 'bet-' )then
            if( lkh .ne. ' ffn')then
             if( lkh .ne. 'bkmo' )then
              if( lkh .ne. 'mo92' )then
               if(irank.eq.0) write(*,*)'3error in beta decay'
c             stop
              endif
             endif
            endif
           endif
         endif
         if( ncharge .eq. -1 )then
           if( lkh .ne. 'bet+')then
             if( lkh .ne. ' bec')then
               if( lkh .ne. '  ec')then
                if( lkh .ne. ' ffn')then
                 if( lkh .ne. 'bkmo' )then
                  if( lkh .ne. 'btyk' )then
                 if(irank.eq.0) write(*,*)'3error in positron decay'
c                 stop
                  endif
                 endif
                endif
               endif
             endif
           endif
         endif
      endif
c
c--deck 4--------------------------------------------
      if( ktype .eq. 4 .and. nam(1) .ne. blank )then
         if( id(1) .gt. 0 .and. id(1) .le. niso )then
           rnam(1) = inam(id(1))
           izr(1) = isoz(id(1))
           inr(1) = ison(id(1))
         else
           rnam(1) = blank
           izr(1) = 0
           inr(1) = 0
         endif
         if( id(2) .gt. 0 .and. id(2) .le. niso )then
           rnam(2) = inam(id(2))
           izr(2) = isoz(id(2))
           inr(2) = ison(id(2))
         else
           rnam(2) = blank
           izr(2) = 0
           inr(2) = 0
         endif
         if( id(3) .gt. 0 .and. id(3) .le. niso )then
           rnam(3) = inam(id(3))
           izr(3) = isoz(id(3))
           inr(3) = ison(id(3))
         else
           rnam(3) = blank
           izr(3) = 0
           inr(3) = 0
         endif
         ncharge = izr(3) - izr(2) - izr(1)
         nnucleon = ncharge + inr(3) - inr(2) - inr(1)
         if( ibug .ne. 0 )then
          if(irank.eq.0) write(*,'(3(a5,3i5,5x),2i5,2x,a4,2a1)')
     1        rnam(1),izr(1),inr(1),id(1),
     1        rnam(2),izr(2),inr(2),id(2),
     1        rnam(3),izr(3),inr(3),id(3),
     1        ncharge,nnucleon,lkh,nr,vw
          write(3,'(3(a5,3i5,5x),2i5,2x,a4,2a1)')
     1        rnam(1),izr(1),inr(1),id(1),
     1        rnam(2),izr(2),inr(2),id(2),
     1        rnam(3),izr(3),inr(3),id(3),
     1        ncharge,nnucleon,lkh,nr,vw
         endif
c..check reaction identities
         if( nnucleon .ne. 0 )then
          if(irank.eq.0) write(*,*)'nucleon nonconservation',ktype
           stop
         endif
         if( ncharge .gt. 1 .or. ncharge .lt. -1 )then
           if(irank.eq.0) write(*,*)
     1        '4error in beta decay, ncharge=',ncharge
         endif
         if( ncharge .eq. 1 )then
           if( lkh .ne. 'bet-' )then
             if(irank.eq.0) write(*,*)'4error in beta decay'
             stop
           endif
         endif
         if( ncharge .eq. -1 )then
           if( lkh .ne. 'bet+')then
             if( lkh .ne. ' bec')then
               if( lkh .ne. '  ec')then
                 if(irank.eq.0) write(*,*)'4error in positron decay'
                 stop
               endif
             endif
           endif
         endif
      endif
c
c--deck 5--------------------------------------------
      if( ktype .eq. 5 .and. nam(1) .ne. blank )then
         if( id(1) .gt. 0 .and. id(1) .le. niso )then
           rnam(1) = inam(id(1))
           izr(1) = isoz(id(1))
           inr(1) = ison(id(1))
         else
           rnam(1) = blank
           izr(1) = 0
           inr(1) = 0
         endif
         if( id(2) .gt. 0 .and. id(2) .le. niso )then
           rnam(2) = inam(id(2))
           izr(2) = isoz(id(2))
           inr(2) = ison(id(2))
         else
           rnam(2) = blank
           izr(2) = 0
           inr(2) = 0
         endif
         if( id(3) .gt. 0 .and. id(3) .le. niso )then
           rnam(3) = inam(id(3))
           izr(3) = isoz(id(3))
           inr(3) = ison(id(3))
         else
           rnam(3) = blank
           izr(3) = 0
           inr(3) = 0
         endif
         if( id(4) .gt. 0 .and. id(4) .le. niso )then
           rnam(4) = inam(id(4))
           izr(4) = isoz(id(4))
           inr(4) = ison(id(4))
         else
           rnam(4) = blank
           izr(4) = 0
           inr(4) = 0
         endif
         ncharge  = izr(4) + izr(3) - izr(2) - izr(1)
         nnucleon = ncharge
     1            + inr(4) + inr(3) - inr(2) - inr(1)
         if( ibug .ne. 0 )then
         if(irank.eq.0) write(*,'(4(a5,2i3,i4,3x),2i5,2x,a4,2a1)')
     1        rnam(1),izr(1),inr(1),id(1),
     1        rnam(2),izr(2),inr(2),id(2),
     1        rnam(3),izr(3),inr(3),id(3),
     1        rnam(4),izr(4),inr(4),id(4),
     1        ncharge,nnucleon,lkh,nr,vw
         write(3,'(4(a5,2i3,i4,3x),2i5,2x,a4,2a1)')
     1        rnam(1),izr(1),inr(1),id(1),
     1        rnam(2),izr(2),inr(2),id(2),
     1        rnam(3),izr(3),inr(3),id(3),
     1        rnam(4),izr(4),inr(4),id(4),
     1        ncharge,nnucleon,lkh,nr,vw
         endif
c..check reaction identities
         if( nnucleon .ne. 0 )then
           if(irank.eq.0) write(*,*)'nucleon nonconservation',ktype
           stop
         endif
         if( ncharge .gt. 1 .or. ncharge .lt. -1 )then
           if(irank.eq.0) write(*,*)
     1       '5error in beta decay, ncharge=',ncharge
         endif
         if( ncharge .eq. 1 )then
           if( lkh .ne. 'bet-' )then
             if(irank.eq.0) write(*,*)'5error in beta decay'
             stop
           endif
         endif
         if( ncharge .eq. -1 )then
           if( lkh .ne. 'bet+')then
             if( lkh .ne. ' bec')then
               if( lkh .ne. '  ec')then
                 if(irank.eq.0) write(*,*)'5error in positron decay'
                 stop
               endif
             endif
           endif
         endif
      endif
c
c--deck 6--------------------------------------------
      if( ktype .eq. 6 .and. nam(1) .ne. blank )then
         if( id(1) .gt. 0 .and. id(1) .le. niso )then
           rnam(1) = inam(id(1))
           izr(1) = isoz(id(1))
           inr(1) = ison(id(1))
         else
           rnam(1) = blank
           izr(1) = 0
           inr(1) = 0
         endif
         if( id(2) .gt. 0 .and. id(2) .le. niso )then
           rnam(2) = inam(id(2))
           izr(2) = isoz(id(2))
           inr(2) = ison(id(2))
         else
           rnam(2) = blank
           izr(2) = 0
           inr(2) = 0
         endif
         if( id(3) .gt. 0 .and. id(3) .le. niso )then
           rnam(3) = inam(id(3))
           izr(3) = isoz(id(3))
           inr(3) = ison(id(3))
         else
           rnam(3) = blank
           izr(3) = 0
           inr(3) = 0
         endif
         if( id(4) .gt. 0 .and. id(4) .le. niso )then
           rnam(4) = inam(id(4))
           izr(4) = isoz(id(4))
           inr(4) = ison(id(4))
         else
           rnam(4) = blank
           izr(4) = 0
           inr(4) = 0
         endif
         if( id(5) .gt. 0 .and. id(5) .le. niso )then
           rnam(5) = inam(id(5))
           izr(5) = isoz(id(5))
           inr(5) = ison(id(5))
         else
           rnam(5) = blank
           izr(5) = 0
           inr(5) = 0
         endif
         ncharge  = izr(5) + izr(4) + izr(3) - izr(2) - izr(1)
         nnucleon = ncharge
     1            + inr(5) + inr(4) + inr(3) - inr(2) - inr(1)
         if( ibug .ne. 0 )then
          if(irank.eq.0) write(*,'(5(a5,2i3,i4,2x),2i5,1x,a4,2a1)')
     1        rnam(1),izr(1),inr(1),id(1),
     1        rnam(2),izr(2),inr(2),id(2),
     1        rnam(3),izr(3),inr(3),id(3),
     1        rnam(4),izr(4),inr(4),id(4),
     1        rnam(5),izr(5),inr(5),id(5),
     1        ncharge,nnucleon,lkh,nr,vw
          write(3,'(5(a5,2i3,i4,2x),2i5,1x,a4,2a1)')
     1        rnam(1),izr(1),inr(1),id(1),
     1        rnam(2),izr(2),inr(2),id(2),
     1        rnam(3),izr(3),inr(3),id(3),
     1        rnam(4),izr(4),inr(4),id(4),
     1        rnam(5),izr(5),inr(5),id(5),
     1        ncharge,nnucleon,lkh,nr,vw
         endif
c..check reaction identities
         if( nnucleon .ne. 0 )then
           if(irank.eq.0) write(*,*)'nucleon nonconservation',ktype
           stop
         endif
         if( ncharge .gt. 1 .or. ncharge .lt. -1 )then
           if(irank.eq.0) write(*,*)
     1        '6error in beta decay, ncharge=',ncharge
         endif
         if( ncharge .eq. 1 )then
           if( lkh .ne. 'bet-' )then
             if(irank.eq.0) write(*,*)'6error in beta decay'
             stop
           endif
         endif
         if( ncharge .eq. -1 )then
           if( lkh .ne. 'bet+')then
             if( lkh .ne. ' bec')then
               if( lkh .ne. '  ec')then
                 if(irank.eq.0) write(*,*)'6error in positron decay'
                 stop
               endif
             endif
           endif
         endif
      endif
c--deck 7--------------------------------------------
      if( ktype .eq. 7 .and. nam(1) .ne. blank )then
         if( id(1) .gt. 0 .and. id(1) .le. niso )then
           rnam(1) = inam(id(1))
           izr(1) = isoz(id(1))
           inr(1) = ison(id(1))
         else
           rnam(1) = blank
           izr(1) = 0
           inr(1) = 0
         endif
         if( id(2) .gt. 0 .and. id(2) .le. niso )then
           rnam(2) = inam(id(2))
           izr(2) = isoz(id(2))
           inr(2) = ison(id(2))
         else
           rnam(2) = blank
           izr(2) = 0
           inr(2) = 0
         endif
         if( id(3) .gt. 0 .and. id(3) .le. niso )then
           rnam(3) = inam(id(3))
           izr(3) = isoz(id(3))
           inr(3) = ison(id(3))
         else
           rnam(3) = blank
           izr(3) = 0
           inr(3) = 0
         endif
         if( id(4) .gt. 0 .and. id(4) .le. niso )then
           rnam(4) = inam(id(4))
           izr(4) = isoz(id(4))
           inr(4) = ison(id(4))
         else
           rnam(4) = blank
           izr(4) = 0
           inr(4) = 0
         endif
         if( id(5) .gt. 0 .and. id(5) .le. niso )then
           rnam(5) = inam(id(5))
           izr(5) = isoz(id(5))
           inr(5) = ison(id(5))
         else
           rnam(5) = blank
           izr(5) = 0
           inr(5) = 0
         endif
         if( id(6) .gt. 0 .and. id(6) .le. niso )then
           rnam(6) = inam(id(6))
           izr(6) = isoz(id(6))
           inr(6) = ison(id(6))
         else
           rnam(6) = blank
           izr(6) = 0
           inr(6) = 0
         endif
         ncharge  = izr(6) + izr(5) + izr(4)
     1            + izr(3) - izr(2) - izr(1)
         nnucleon = ncharge
     1            + inr(6) + inr(5) + inr(4)
     2            + inr(3) - inr(2) - inr(1)
         if( ibug .ne. 0 )then
          if(irank.eq.0) write(*,'(6(a5,2i3,i4,2x),2i5,1x,a4,2a1)')
     1        rnam(1),izr(1),inr(1),id(1),
     1        rnam(2),izr(2),inr(2),id(2),
     1        rnam(3),izr(3),inr(3),id(3),
     1        rnam(4),izr(4),inr(4),id(4),
     1        rnam(5),izr(5),inr(5),id(5),
     1        rnam(6),izr(6),inr(6),id(6),
     1        ncharge,nnucleon,lkh,nr,vw
          write(3,'(6(a5,2i3,i4,2x),2i5,1x,a4,2a1)')
     1        rnam(1),izr(1),inr(1),id(1),
     1        rnam(2),izr(2),inr(2),id(2),
     1        rnam(3),izr(3),inr(3),id(3),
     1        rnam(4),izr(4),inr(4),id(4),
     1        rnam(5),izr(5),inr(5),id(5),
     1        rnam(6),izr(6),inr(6),id(6),
     1        ncharge,nnucleon,lkh,nr,vw
         endif
c..check reaction identities
         if( nnucleon .ne. 0 )then
           if(irank.eq.0) write(*,*)'nucleon nonconservation',ktype
           stop
         endif
         if( ncharge .gt. 1 .or. ncharge .lt. -1 )then
           if(irank.eq.0) write(*,*)
     1        '7error in beta decay, ncharge=',ncharge
         endif
         if( ncharge .eq. 1 )then
           if( lkh .ne. 'bet-' )then
             if(irank.eq.0) write(*,*)'7error in beta decay'
             stop
           endif
         endif
         if( ncharge .eq. -1 )then
           if( lkh .ne. 'bet+')then
             if( lkh .ne. ' bec')then
               if( lkh .ne. '  ec')then
                 if(irank.eq.0) write(*,*)'7error in positron decay'
                 stop
               endif
             endif
           endif
         endif
      endif

c--deck 8--------------------------------------------
      if( ktype .eq. 8 .and. nam(1) .ne. blank )then
         if( id(1) .gt. 0 .and. id(1) .le. niso )then
           rnam(1) = inam(id(1))
           izr(1) = isoz(id(1))
           inr(1) = ison(id(1))
         else
           rnam(1) = blank
           izr(1) = 0
           inr(1) = 0
         endif
         if( id(2) .gt. 0 .and. id(2) .le. niso )then
           rnam(2) = inam(id(2))
           izr(2) = isoz(id(2))
           inr(2) = ison(id(2))
         else
           rnam(2) = blank
           izr(2) = 0
           inr(2) = 0
         endif
         if( id(3) .gt. 0 .and. id(3) .le. niso )then
           rnam(3) = inam(id(3))
           izr(3) = isoz(id(3))
           inr(3) = ison(id(3))
         else
           rnam(3) = blank
           izr(3) = 0
           inr(3) = 0
         endif
         if( id(4) .gt. 0 .and. id(4) .le. niso )then
           rnam(4) = inam(id(4))
           izr(4) = isoz(id(4))
           inr(4) = ison(id(4))
         else
           rnam(4) = blank
           izr(4) = 0
           inr(4) = 0
         endif
         if( id(5) .gt. 0 .and. id(5) .le. niso )then
           rnam(5) = inam(id(5))
           izr(5) = isoz(id(5))
           inr(5) = ison(id(5))
         else
           rnam(5) = blank
           izr(5) = 0
           inr(5) = 0
         endif
         ncharge  = izr(5) + izr(4) - izr(3) - izr(2) - izr(1)
         nnucleon = ncharge
     1            + inr(5) + inr(4) - inr(3) - inr(2) - inr(1)
         if( ibug .ne. 0 )then
          if(irank.eq.0) write(*,'(5(a5,2i3,i4,3x),2i5,2x,a4,2a1)')
     1        rnam(1),izr(1),inr(1),id(1),
     1        rnam(2),izr(2),inr(2),id(2),
     1        rnam(3),izr(3),inr(3),id(3),
     1        rnam(4),izr(4),inr(4),id(4),
     1        rnam(5),izr(5),inr(5),id(5),
     1        ncharge,nnucleon,lkh,nr,vw
          write(3,'(5(a5,2i3,i4,3x),2i5,2x,a4,2a1)')
     1        rnam(1),izr(1),inr(1),id(1),
     1        rnam(2),izr(2),inr(2),id(2),
     1        rnam(3),izr(3),inr(3),id(3),
     1        rnam(4),izr(4),inr(4),id(4),
     1        rnam(5),izr(5),inr(5),id(5),
     1        ncharge,nnucleon,lkh,nr,vw
         endif
c..check reaction identities
         if( nnucleon .ne. 0 )then
           if(irank.eq.0) write(*,*)'nucleon nonconservation',ktype
           stop
         endif
         if( ncharge .gt. 1 .or. ncharge .lt. -1 )then
           if(irank.eq.0) write(*,*)
     1        '8error in beta decay, ncharge=',ncharge
         endif
         if( ncharge .eq. 1 )then
           if( lkh .ne. 'bet-' )then
             if(irank.eq.0) write(*,*)'8error in beta decay'
             stop
           endif
         endif
         if( ncharge .eq. -1 )then
           if( lkh .ne. 'bet+')then
             if( lkh .ne. ' bec')then
               if( lkh .ne. '  ec')then
                 if(irank.eq.0) write(*,*)'8error in positron decay'
                 stop
               endif
             endif
           endif
         endif
      endif
      do nc = 1, 6
      if( nam(nc) .ne. blank )then

c        if(  nnn(izr(nc),inr(nc))  .eq. 0 )then
c          if( izr(nc) .ne. 0 .and. inr(nc) .ne. 0)then

        if( izr(nc) .ne. 0 .and. inr(nc) .ne. 0)then
          if(  nnn(izr(nc),inr(nc))  .eq. 0 )then
            nnn(izr(nc),inr(nc)) = id(nc)
          endif
        endif
      endif

      enddo

c-------------------------------------------------------

      goto 20

90    continue
      close(5)
c      close(6)
      if(irank.eq.0) write(*,*)' getfkt has finished reading netsu',
     1      kk,' reaction entries'

c      write(3,*)' getfkt has finished reading netsu',
c     1      kk,' reaction entries'

c--------------------------------------------------------------------
c--------------------------------------------------------------------

      number = 0
      do j=1,nnmax
        do i=1,nzmax
c..do not count alpha as nucleus
          if(nnn(i,j) .ne. 0 .and. nnn(i,j) .ne. 6 )then
            number = number + 1
            if(irank.eq.0) write(*,'(4i5,a6)')
     1         i,j,nnn(i,j),number,inam(nnn(i,j))
cccccccccccccccccccccc


          endif
        enddo
      enddo
      itot = number + 3

      if(irank.eq.0) write(*,*)itot,' nuclei',nzmax,nnmax
c      stop'sssssssssssssss'

      if( itot+1 .gt. ndim )then
        if(irank.eq.0) write(*,'(a20,2(a10,i5))')
     1     'getfkt:','itot+1',itot+1,'ndim',ndim
        if(irank.eq.0) write(*,*)'network too large for dimensions'
        if(irank.eq.0) write(*,*)'in dimenfile,'
        if(irank.eq.0) write(*,*)'change ndim to .ge. ',itot+1
        if(irank.eq.0) write(*,*)'change nreac to .ge. ',kk

        stop'getfkt 1'
      endif
      inuc = number
      ireac = kk

c..setup nz,nn,xid vectors for burn
      l = 0
      do k = 1, niso
        do j=1,nnmax
          do i=1,nzmax
c..do not count alpha as nucleus
            if( k .eq. nnn(i,j) .and. k .ne. 6 )then
              l = l+1
              nz(l) = isoz(k)
              nn(l) = ison(k)
              xid(l) = inam(k)
              if( ibug .ne. 0 )then
                if(irank.eq.0) write(*,'(3i5,2x,a5)')
     1             l,nz(l),nn(l),xid(l)
c                write(3,'(3i5,2x,a5)')l,nz(l),nn(l),xid(l)
              endif
            endif
          enddo
        enddo
      enddo
      
c==============================================
c 1,2, and 6 are n,p, and alphas in isotope.lib
c==============================================
      xid(number+1) = inam(1)
      nz (number+1) = isoz(1)
      nn (number+1) = ison(1)
      xid(number+2) = inam(2)
      nz (number+2) = isoz(2)
      nn (number+2) = ison(2)
      xid(number+3) = inam(6)
      nz (number+3) = isoz(6)
      nn (number+3) = ison(6)
c electrons are in number + 4 = ndim usually
      xid(number+4) = '   Ye'

      if( ibug .ne. 0 )then
        do l = number+1, number+3
              if(irank.eq.0) write(*,'(3i5,2x,a5)')l,nz(l),nn(l),xid(l)
c              write(3,'(3i5,2x,a5)')l,nz(l),nn(l),xid(l)
        enddo
      endif
c
c..count number of reactions of each type (deck)
      do j = 1, 8
        ndeck(j) = 0
      enddo
      do i = 1, ireac
        do j = 1,8
        if( ideck(i) .eq. j )ndeck(j) = ndeck(j) + 1
        enddo
      enddo

      if(irank.eq.0) write(*,*)
     1 'number of reactions ndeck(n) in each deck (n=1 to 8)'
      if(irank.eq.0) write(*,'(8i5)')ndeck
      if( ndeck(2) .ne. ndeck(4) )then
        if(irank.eq.0) write(*,'(a50)')
     1 'WARNING: 1-->2+3 does not match inverse 2+3-->1'
        if(irank.eq.0) write(*,'(a50)')
     1     '         because ndeck(2) .ne. ndeck(4)'
      endif
      if( mod( ndeck(5),2 ) .ne. 0 )then
        if(irank.eq.0) write(*,'(a50)')
     1 'WARNING: 1-->2+3 does not match inverse 2+3-->1'
        if(irank.eq.0) write(*,'(a50)')
     1     '         because ndeck(2) .ne. ndeck(4)'
      endif

c      write(3,*)'number of reactions in each deck (1 to 8)'
c      write(3,'(8i5)')ndeck
c      if( ndeck(2) .ne. ndeck(4) )then
c        write(3,'(a50)')
c     1  'WARNING: 1-->2+3 does not match inverse 2+3-->1'
c        write(3,'(a50)')'         because ndeck(2) .ne. ndeck(4)'
c      endif
c      if( mod( ndeck(5),2 ) .ne. 0 )then
c         write(3,'(a50)')'WARNING: ndeck(5) is not even, '
c         write(3,'(a50)')'         so inverses may not balance'
c      endif

      return
      end



