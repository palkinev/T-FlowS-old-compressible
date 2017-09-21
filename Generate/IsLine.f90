!======================================================================!
  INTEGER FUNCTION IsLine(n1, n2, b) 
!----------------------------------------------------------------------!
!   Checks if the line defined n1 and n2 is inside the block b.        !
!----------------------------------------------------------------------!
!------------------------------[Modules]-------------------------------!
  USE all_mod
  USE gen_mod
!----------------------------------------------------------------------! 
  IMPLICIT NONE
!-----------------------------[Parameters]-----------------------------!
  INTEGER :: b, n1, n2
!-------------------------------[Locals]-------------------------------!
  INTEGER :: l1, l2
!--------------------------------[CVS]---------------------------------!
!  $Id: IsLine.f90,v 1.1 2014/11/24 11:31:30 muhamed Exp $  
!  $Source: /home/mhadziabdic/Dropbox/cvsroot/T-FlowS-CVS/Generate/IsLine.f90,v $  
!======================================================================!

  do l1=1,8
    do l2=1,8
      if( (BlkPnt(b,l1) == n1) .and.                                &
	  (BlkPnt(b,l2) == n2) ) then
	   goto 1
      end if 
    end do
  end do     

  IsLine=0
  return

1 IsLine=b
  return

  END FUNCTION IsLine
