!======================================================================!
  SUBROUTINE EndPar
!----------------------------------------------------------------------!
!   Ends parallel execution.                                           !
!----------------------------------------------------------------------!
  IMPLICIT NONE
!------------------------------[Include]-------------------------------!
  INCLUDE 'mpif.h'
!-------------------------------[Locals]-------------------------------!
  INTEGER :: ERROR
!======================================================================!

  call MPI_FINALIZE(ERROR)

  END SUBROUTINE EndPar
