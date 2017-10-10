!======================================================================!
  subroutine Prec_Form(N, A, D, prec) 
!----------------------------------------------------------------------!
! Forms preconditioning matrix "D" from provided matrix "A".           !
!----------------------------------------------------------------------!
!------------------------------[Modules]-------------------------------!
  use allt_mod, only: Matrix
  use par_mod
!----------------------------------------------------------------------!
  implicit none
!-----------------------------[Parameters]-----------------------------!
  integer      :: N
  type(Matrix) :: A
  type(Matrix) :: D    
  integer      :: prec
!-------------------------------[Locals]-------------------------------!
  real     :: sum1
  integer  :: i, j, k
!======================================================================!
                 
!+++++++++++++++++++++++++!
!     preconditioning     !
!+++++++++++++++++++++++++!

!----- 1) diagonal preconditioning -----!
!     => Parallelized                   ! 
  if(prec == 1) then        
    do i=1,N                     
      D % val(D % dia(i)) = A % val(A % dia(i))           
    end do                      

!----- 2) incomplete cholesky preconditioning -----!
  else if(prec == 2) then   
    do i = 1,N
      sum1 = A % val(A % dia(i))       ! take diaginal entry   
      do j = A % col(i), A % dia(i)-1  ! only lower traingular
        k = A % row(j)                    
        sum1 = sum1 - D % val(D % dia(k)) * A % val(j) * A % val(j)  
      end do
      D % val(D % dia(i)) = 1.0 / sum1
    end do

!----- .) no preconditioning -----!
!     => Parallelized             ! 
  else                          
    do i=1,N
      D % val(D % dia(i)) = 1.0
    end do
  end if 

  end subroutine Prec_Form