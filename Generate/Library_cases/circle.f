      PROGRAM circle

      IMPLICIT NONE

      REAL x0, y0, z0, r0, PI, PHI, DelPHI
      REAL x, y, z
      INTEGER  i,N

      PI=2.0*asin(1.0)

      write(*,*) PI 

      write(*,*) 'Input center: '
      read(*,*) x0, y0, z0


      write(*,*) 'Input radius: '
      read(*,*) r0 

      write(*,*) 'Input number of segments: '
      read(*,*) N

      DelPHI = 2.0*PI/(real(N))
 
      do i=0,N-1
        PHI = DelPHI * i
        x = x0 + r0 * cos(PHI)
        z = z0 + r0 * sin(PHI)
        y = y0
        write(*,'(I5,3F12.6)') i+1, x, y, z
      end do

      END 



      
