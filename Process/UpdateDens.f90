!======================================================================!
  SUBROUTINE UpdateDens(DENS, SCAL,n,sigma)
!----------------------------------------------------------------------!
!   Discretizes and solves momentum conservation equations             !
!----------------------------------------------------------------------!
!------------------------------[Modules]-------------------------------!
  USE all_mod
  USE pro_mod
  USE les_mod
  USE rans_mod
  USE par_mod
!----------------------------------------------------------------------!
  IMPLICIT NONE
!-----------------------------[Parameters]-----------------------------!
  TYPE(Unknown) :: DENS, SCAL
  REAL          :: temper(NC), alpha_tem, SCAL_st, Rho_0, sigma
  INTEGER       :: n
!-------------------------------[Locals]-------------------------------!
  INTEGER :: s, c, c1, c2
!======================================================================!

  do c=1,NC

!    DENS % n(c) =  (SCAL % n(c) / 1.0 + (1 - SCAL % n(c)) / 20.0)**(-1.0)

  end do

  call Exchng(DENS % n)
  call Exchng(Zmix % n)

  END SUBROUTINE UpdateDens
