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
!--------------------------------[CVS]---------------------------------!
!  $Id: EndPar.f90,v 1.1 2002/11/01 15:12:12 niceno Exp $  
!  $Source: /home/muhamed/.CVSROOT/T-Rex/Parallel/Double/EndPar.f90,v $  
!======================================================================!

  call MPI_FINALIZE(ERROR)

  END SUBROUTINE EndPar
