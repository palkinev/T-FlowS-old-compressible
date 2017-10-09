!======================================================================!
  subroutine GloSum(phi) 
!----------------------------------------------------------------------!
!   Estimates global summ among all processors.                        !
!----------------------------------------------------------------------!
  implicit none
!------------------------------[Include]-------------------------------!
  include 'mpif.h'
!-----------------------------[Parameters]-----------------------------!
  real    :: phi
!-------------------------------[Locals]-------------------------------!
  real    :: phi_new
  integer :: error
!======================================================================!

!================================================
      call MPI_ALLREDUCE      &               
!-----------------------------------+------------
             (phi,            & ! send buffer
              phi_new,        & ! recv buffer 
              1,              & ! length     
              MPI_REAL,       & ! datatype  
              MPI_SUM,        & ! operation 
              MPI_COMM_WORLD, &             
              error) 
!================================================

  phi = phi_new

  end subroutine GloSum
