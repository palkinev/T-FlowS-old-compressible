!======================================================================!
  SUBROUTINE Wait 
!----------------------------------------------------------------------!
!   Puts barrier for parallel execution.                               !
!----------------------------------------------------------------------!
  IMPLICIT NONE
!------------------------------[Include]-------------------------------!
  INCLUDE 'mpif.h'
!-------------------------------[Locals]-------------------------------!
  INTEGER :: error
!--------------------------------[CVS]---------------------------------!
!  $Id: Wait.f90,v 1.1 2014/11/24 11:39:07 muhamed Exp $  
!  $Source: /home/mhadziabdic/Dropbox/cvsroot/T-FlowS-CVS/Parallel/Single/Wait.f90,v $  
!======================================================================!

!==================================
      call MPI_BARRIER        &
!----------------------------------
	     (MPI_COMM_WORLD, &
	      error) 
!==================================

  END SUBROUTINE Wait
