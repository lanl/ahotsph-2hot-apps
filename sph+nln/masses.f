      subroutine masses(qq,w,itot,irank)
      implicit none

c..sets up the binding energies using Friedel's NETWInV file.
c..searches for n,p,alpha and puts at end of list

      include 'dimen'

      integer*4 irank

      character*5  nuc(nnuc)
      character*5  cdummy
      character*72 line
      real*8       mxcess(nnuc)

      integer*4 itot, i, j, inuc, ih1, in1, ia4
      integer*4 nt(24), nz(nnuc), n(nnuc)
      real*8    p(nnuc), a(nnuc), qq(nnuc+1),w(nnuc+1)

      real*8     dum(nnuc),dumw(nnuc)
      real*8     p1,p2,p3, w1,w2,w3

      logical tobe
c--------------------------------------------------------------

      inquire(file='netwinv',exist=tobe)
      if( .not. tobe )then
         if(irank.eq.0) write(*,*)
     1      'masses.f: no netwinv file in this directory'
      endif

      open (12, file='netwinv', status='old')

c..     Read in the temperature array on the second line of NETWInV:
c..     These are the T's at which Friedel's partition functions
c..     are defined.

      read (12,'(a72)')line
      read (12,1) (nt(i),i=1,24)
    1 format (24i3)

c..skip the nuclear symbols

      
cccccccccccccccccccccccccccccccccccccccccc?????????????????????????
      
      do i=1,itot
         read (12,'(a5)')cdummy
      enddo

c..     Start reading

      i = 1
 1000 continue
      read (12,*,end=2000) nuc(i), a(i), nz(i), n(i), w(i), mxcess(i)
c      write(*,*) nuc(i), a(i),nz(i), n(i), w(i), mxcess(i)
      read (12,4) (p(j),j=1,24)
    2 format (a5,4x,f9.3,i3,1x,i3,2x,f5.1,2x,f7.3)
    4 format(8f9.2)

c cie: something seems to be stomping on the memory of some of these
c      arrays as or some time after they are read in (precise timing
c      varies). So this is observed in the first couple of elements 
c      of either 'n' or 'nz'. Don't know what's causing this. But,
c      this results in a bug in the reordering done below that results
c      in 'w' containing zeros in its last three elements. 

c..     update count for reading next line
      i = i+1
      goto 1000

c..     All nuclei read
 2000 inuc = i-1

      if(irank.eq.0) print *,'MASSES:', inuc, ' nuclei found'

c..construct binding energy from mass excesses
c..use mass excesses
      do i = 1, itot
         qq(i) = mxcess(i)
      enddo

    3 format (i4,2x,a5,4x,f9.3,i3,1x,i3,2x,f5.1,2x,f7.3,2x,f8.3)

c..change J's to 2J+1
      do i = 1, itot
        w(i) = 2.0d0*w(i) + 1.0d0
      enddo

c..     Rearrange with n, p, and alpha at the end (this improves
c..efficiency of most linear equation solvers)

      do i=1,itot
        dum(i) = qq(i)
        dumw(i) = w(i)
      enddo
c..find them
      ih1 = 0
      in1 = 0
      ia4 = 0
      do i = 1, itot
        if( nz(i) .eq. 0 .and. n(i) .eq. 1 )then
           in1 = i
        endif
        if( nz(i) .eq. 1 .and. n(i) .eq. 0 )then
           ih1 = i
        endif
        if( nz(i) .eq. 2 .and. n(i) .eq. 2 )then
           ia4 = i
        endif
      enddo

      p1 = qq(in1)
      p2 = qq(ih1)
      p3 = qq(ia4)
      w1 =  w(in1)
      w2 =  w(ih1)
      w3 =  w(ia4)

      j = 1
      do i=1,itot
         if( i .ne. in1 .and. i .ne. ih1 .and. i .ne. ia4 )then
           qq(j) = dum(i)
           w(j)  = dumw(i)
           j = j+1
         endif
      enddo

      qq(itot-2) = p1
      qq(itot-1) = p2
      qq(itot  ) = p3
      w(itot-2)  = w1
      w(itot-1)  = w2
      w(itot  )  = w3

      close(12)

      return
      end



