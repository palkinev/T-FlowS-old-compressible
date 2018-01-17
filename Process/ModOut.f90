!==============================================================================!
  subroutine ModOut
!------------------------------------------------------------------------------!
!   Modifies the fluxes at outflow boundaries to conserve the mass             !
!   in the whole domain                                                        !
!----------------------------------[Modules]-----------------------------------!
  use allp_mod, only: convect, inflow, outflow, pressure, tiny
  use all_mod
  use pro_mod
  use Moin_Problem_Mod
!------------------------------------------------------------------------------!
  implicit none
!-----------------------------------[Locals]-----------------------------------!
  integer :: m, s, c1, c2
  real    :: fac(256)
  real    :: delta_mass, RhoS
!---------------------------------[Interface]----------------------------------!
  ! these function must be pushed into par_mod, then interface is not needed
  interface
    subroutine glosum(phi) 
      implicit none
      include 'mpif.h'
      real    :: phi
      real    :: phinew
      integer :: error

    end subroutine glosum
  end interface
!==============================================================================!


  !--------------------------------------------------------!
  !   Calculate density change in domain = d rho / dt*dV   !
  !--------------------------------------------------------!

  ! - d rho / dt dv
  delta_mass = 0.0
  do m = 1, Nmat
    delta_mass = sum(volume(1: NC)/dt* &
    (1.5*Rho % n(1: NC) - 2.0*Rho % o(1: NC) + 0.5*Rho % oo(1: NC)))

    call glosum(delta_mass)
  end do

  !--------------------------------------!
  !   Calculate the inflow mass fluxes   !
  !--------------------------------------!

  do m = 1, Nmat
    massin(m) = 0.0
    do s = 1, NS
      c1 = SideC(1,s)
      c2 = SideC(2,s)
      if(c2 < 0) then

        if(TypeBC(c2)  ==  INFLOW) then
          if(material(c1) == m) massin(m) = massin(m) + Flux(s)
        endif
        if(TypeBC(c2)  ==  PRESSURE .and. Flux(s) < 0.0) then
          if(material(c1) == m) massin(m) = massin(m) + Flux(s)
        end if
        if(TypeBC(c2)  ==  CONVECT .and. Flux(s) < 0.0) then
          if(material(c1) == m) massin(m) = massin(m) + Flux(s)
        end if
      end if
    end do
    call glosum(massin(m))
  end do

  !---------------------------------------!
  !   Calculate the outflow mass fluxes   !
  !---------------------------------------!

  do m = 1, Nmat
    masout(m) = 0.0
    do s = 1, NS
      c1 = SideC(1,s)
      c2 = SideC(2,s)
      if(c2  < 0) then
        if(TypeBC(c2) == OUTFLOW) then
          if(material(c1) == m) masout(m) = masout(m) + Flux(s)
        endif
        if(TypeBC(c2) == CONVECT .and. Flux(s) > 0.0) then
          if(material(c1) == m) masout(m) = masout(m) + Flux(s)
        endif
        if(TypeBC(c2) == PRESSURE .and. Flux(s) > 0.0) then
          if(material(c1) == m) masout(m) = masout(m) + Flux(s)
        endif

      endif
    end do
    call glosum(masout(m))
  end do

  !---------------------------------------!
  !   Correct it via fac to satisfy the   !
  !   overall mass balance                !
  !---------------------------------------!
  do m = 1, Nmat
    fac(m) = -(massin(m)+delta_mass)/(masout(m)+TINY)
  end do

  do m = 1, Nmat
    do s = 1, NS
      c1 = SideC(1,s)
      c2 = SideC(2,s)
      if(c2  < 0) then
        if(TypeBC(c2) == OUTFLOW .or. TypeBC(c2) == CONVECT .or. TypeBC(c2) == PRESSURE) then
          if(material(c1) == m) then
      
      
            U % n(c2) = rho % n(c1) * U % n(c1) / rho % n(c2) * fac(m)
            V % n(c2) = rho % n(c1) * V % n(c1) / rho % n(c2) * fac(m)
            W % n(c2) = rho % n(c1) * W % n(c1) / rho % n(c2) * fac(m)

            RhoS = f(s) * rho % n(c1) + ( 1.0-f(s) ) * rho % n(c2)

            Flux(s) =   f(s) * ( rho % n(c1) * U % n(c1) * Sx(s) +   &
                                 rho % n(c1) * V % n(c1) * Sy(s) +   &
                                 rho % n(c1) * W % n(c1) * Sz(s) ) + &
                ( 1.0-f(s) ) * ( rho % n(c2) * U % n(c2) * Sx(s) +   &
                                 rho % n(c2) * V % n(c2) * Sy(s) +   &
                                 rho % n(c2) * W % n(c2) * Sz(s) )
            Flux_u(s) = Flux(s) / RhoS

          endif
        endif
      endif
    end do
  end do

end subroutine ModOut