!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
!                                 !                                    !
!                                 !   Bojan Niceno                     !
!   RANS models  variable         !   Delft University of Technology   !
!   definitions for the processor !   Section Heat Transfer             !
!                                 !   niceno@duttwta.wt.tn.tudelft.nl  !
!                                 !                                    !
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
!..RCS/CVS ident
! $Id: rans_mod.h90,v 1.5 2005/01/25 12:38:53 muhamed Exp $
! $Source: /home/muhamed/.CVSROOT/T-Rex/Library/rans_mod.h90,v $
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!

MODULE rans_mod

    USE allp_mod

    IMPLICIT NONE
 
    !----- Constants for the k-eps model:
    REAL, PUBLIC :: Ce1, Ce2, Ce3, Cmu, Cmu25, Cmu75, kappa, Elog
 
    !----- Constants for the k-eps-v2f model:
    REAL, PUBLIC :: CmuD, Cl, Ct, alpha, Cni, cf1, cf2, cf3, Cf_1, Cf_2
    REAL, PUBLIC :: Lim
    REAL, PUBLIC :: g1, g1_star, g2, g3, g3_star, g4, g5

    !----- Constants for the Spalart-Allmaras model:
    REAL, PUBLIC :: Cb1, Cb2, SIGMAv, Cw1, Cw2, Cw3, Cvis1

    !----- Vorticity
    REAL,ALLOCATABLE, PUBLIC :: Vort(:), VortMean(:)

    !----- Turbulent viscosity
    REAL,ALLOCATABLE, PUBLIC :: VISt(:), CmuS(:)
 
    !----- Turbulent conductivity
    REAL,ALLOCATABLE, PUBLIC :: CONt(:)
 
    !----- Lenght and Time Scales
    REAL,ALLOCATABLE, PUBLIC :: Lsc(:)
    REAL,ALLOCATABLE, PUBLIC :: Tsc(:)

    !----- Production of turbulent kinetic energy
    REAL,ALLOCATABLE, PUBLIC :: Pk(:)
 
    !----- Non-dimensional distance
    REAL,ALLOCATABLE, PUBLIC :: Ynd(:)
 
    !----- Friction velocity
    REAL,ALLOCATABLE, PUBLIC :: Uf(:)
    REAL,ALLOCATABLE, PUBLIC :: Ufmean(:)

    !----- Wall viscosity (wall function approuch)
    REAL,ALLOCATABLE, PUBLIC :: VISwall(:)

    !  REAL,ALLOCATABLE, PUBLIC :: EE(:)
    REAL,ALLOCATABLE, PUBLIC :: Fs(:)
    !  REAL,ALLOCATABLE, PUBLIC :: Feps(:)

    REAL,ALLOCATABLE, PUBLIC :: nn1(:)
    REAL,ALLOCATABLE, PUBLIC :: nn2(:)
    REAL,ALLOCATABLE, PUBLIC :: nn3(:)

    REAL,ALLOCATABLE, PUBLIC :: Bud1(:)
    REAL,ALLOCATABLE, PUBLIC :: Bud2(:)
    REAL,ALLOCATABLE, PUBLIC :: Bud3(:)
    REAL,ALLOCATABLE, PUBLIC :: Bud4(:)
    REAL,ALLOCATABLE, PUBLIC :: Bud5(:)
    REAL,ALLOCATABLE, PUBLIC :: Bud6(:)
    REAL,ALLOCATABLE, PUBLIC :: Bud7(:)
    REAL,ALLOCATABLE, PUBLIC :: Bud8(:)
    REAL,ALLOCATABLE, PUBLIC :: Bud9(:)
 
    REAL,ALLOCATABLE, PUBLIC :: uu_star(:)
    REAL,ALLOCATABLE, PUBLIC :: vv_star(:)
    REAL,ALLOCATABLE, PUBLIC :: ww_star(:)
    REAL,ALLOCATABLE, PUBLIC :: uv_star(:)
    REAL,ALLOCATABLE, PUBLIC :: uw_star(:)
    REAL,ALLOCATABLE, PUBLIC :: vw_star(:)
END MODULE 
