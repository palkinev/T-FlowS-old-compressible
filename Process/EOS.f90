!==============================================================================!
  subroutine eos(dens_0, dens_1)
!------------------------------------------------------------------------------!
! Applies equation of state for Moin model problems                            !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use all_mod
  use pro_mod
  use les_mod
  use rans_mod
  use par_mod
  use Moin_Problem_mod !  module to control conflicting Moin modules
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Parameters]---------------------------------!
  real,          intent (in) :: dens_0, dens_1 ! dens_0 >= dens_1
!-----------------------------------[Locals]-----------------------------------!
  integer :: c
!==============================================================================!

  !------------------!
  !   applying EOS   !
  !------------------!

  select case (Moin_Problem())
  case (4) ! RT
    do c = 1, NC
      rho % n(c) =  1.0 / ( 9.0*zmix % n(c) + 1.0)
    end do
  case (10)
    do c = 1, NC
      rho % n(c) = Moin_Problem_EOS (zmix % n(c))
    end do
  case (30)
    do c = 1, NC
      rho % n(c) = Moin_Problem_EOS (zmix % n(c))
    end do
  case default
    do c = 1, NC
      rho % n(c) = 1.0 / ( (zmix % n(c)/dens_1) + (1.0-zmix % n(c))/dens_0 )
    end do
  end select

  !--------------------------------------------------!
  !   applying BC and updating dependent variables   !
  !--------------------------------------------------!

  call CalcDens(rho,.FALSE.,dens_0,dens_1,Moin_Problem())  ! update other variables and apply BC

  !updating dynamic viscosity
  select case (Moin_Problem())
!    case (10)
!        do c = 1, NC
!            visc_dyn(c)     = Moin_Problem_visc_dyn_of_Z(zmix % n(c)) / 224.28074617774
!            visczmix_dyn(c) = visc_dyn(c) / 0.7
!        end do
    !VISc                = minval(visc_dyn)
    !visczmix            = minval(visczmix_dyn)
  end select


  !-----------------------------!
  !   Applying EOS to borders   !
  !-----------------------------!

  select case (Moin_Problem())
  case (4) ! RT
    do c = -NbC,-1
      rho % n(c) =  1.0/( 9.0*zmix % n(c) + 1.0)
    end do
  case (10)
    do c = -NbC,-1
      rho % n(c) = Moin_Problem_EOS (zmix % n(c))
    end do
  case (30)
    do c = -NbC,-1
      rho % n(c) = Moin_Problem_EOS (zmix % n(c))
    end do
  case default
    do c = -NbC,-1
      rho % n(c) = 1.0/( (zmix % n(c)/dens_1) + (1.0-zmix % n(c))/dens_0 )
    end do
  end select

  select case (Moin_Problem())
!    case (10)
!        do c = -NbC, 1
!            visc_dyn(c)     = Moin_Problem_visc_dyn_of_Z(zmix % n(c)) / 224.28074617774
!            visczmix_dyn(c) = visc_dyn(c) / 0.7
!        end do
    !VISc                = minval(visc_dyn)
    !visczmix            = minval(visczmix_dyn)

    !write (*,*) "MIN(visc_dyn)    =", minval(visc_dyn), "visc_dyn(1)=", visc_dyn(1)
    !write (*,*) "MIN(visczmix_dYn)=", minval(visczmix_dyn)
  end select

  !--------------------------------------------------!
  !   applying BC and updating dependent variables   !
  !--------------------------------------------------!
  call wait
  call Exchng(rho % n) !main source for rho presicion
  call Exchng(visc_dyn)
  call Exchng(visczmix_dyn)

end subroutine eos

! if(this < 2) open(991, FILE="out.dat", status='replace', access="sequential", position='append', action='write')
! z_test = 0.0
! do
!     write(991,'(2ES26.16E3)') z_test, Moin_Problem_EOS(z_test)
!     if ( abs (z_test - 1.0) < 1.D-15 ) exit
!     z_test = z_test + 0.001
! end do
! if(this < 2) close(991)