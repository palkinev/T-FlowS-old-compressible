C@a====================================================================*
      PROGRAM profile
*----------------------------------------------------------------------*
C!  Outputs the velocity profile to be used as input for T-Rex         *
C@e--------------------------------------------------------------------*
      IMPLICIT NONE
*----------------------------------------------------------------------*
      REAL*8  L, x, v
      INTEGER N, i
c#rcs-----------------------
      character*80 rcs1,rcs2
      data rcs1/
     >'$Id: profile.f,v 1.3 2002/06/01 17:17:40 bojan Exp $'/
      data rcs2/
     >'$Source: /home/IUS/mhadziabdic/.CVSROOT/T-Rex/Test/profile.f,v $'/
C@d====================================================================*

* Number of points
      N=200

      write(*,'(A20)') '#==================#'
      write(*,'(A20)') '# Number of points #'
      write(*,'(A20)') '#==================#'
      write(*,*) N      
      write(*,'(A13)') '#===========#'
      write(*,'(A13)') '# Direction #'
      write(*,'(A13)') '#===========#'
      write(*,*) 'z'     
      write(*,'(A20)') '#==================#'
      write(*,'(A20)') '# coordinate U V W #'
      write(*,'(A20)') '#==================#'
* Parabolic profile in range 0:2, extrema=1 *
c      L = 2.d0
c      do i=0,N
c        x = real(i) * L/real(n)
c        v = 1.d0 - (x-L/2.d0)**2
c        write(*,'(4F18.6)') x, v, 0.d0, 0.d0
c      end do

* Parabolic profile in range 0:3, extrema=1 *
!      L = 3.d0
!      do i=0,N
!        x = real(i) * L/real(N)
!        v = (2.25-(x-1.5)*(x-1.5))/2.25
!        write(*,'(4F18.6)') x, v, 0.d0, 0.d0
!      end do

* Parabolic profile in range 0:1, extrema=1 *
      L = 1.d0
      do i=0,N
        x = real(i) * L/real(n)
        v = 6.d0 * x * (1.d0 - x)
        write(*,'(5F18.6)') x, v, 0.d0, 0.d0, 0.d0
      end do

      END
