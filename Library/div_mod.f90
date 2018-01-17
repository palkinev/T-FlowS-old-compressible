!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
!                                 !                                    !
!                                 !   Bojan Niceno                     !
!   Global variable definitions   !   Delft University of Technology   !
!         for all modules         !   Section Heat Transfer            !
!                                 !   niceno@duttwta.wt.tn.tudelft.nl  !
!                                 !                                    !
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
!..RCS/CVS ident
! $Id: div_mod.h90,v 1.4 2001/02/19 15:20:03 niceno Exp $
! $Source: /home/muhamed/.CVSROOT/T-Rex/Library/div_mod.h90,v $
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
MODULE div_mod

    IMPLICIT NONE

    INTEGER,ALLOCATABLE, PUBLIC :: ix(:), iy(:), iz(:), iin(:)
    REAL,ALLOCATABLE, PUBLIC    :: criter(:)

    !------------!
    ! Parameters !
    !------------!
    INTEGER, PUBLIC :: ALGOR, COORDINATE, INERTIAL

END MODULE
