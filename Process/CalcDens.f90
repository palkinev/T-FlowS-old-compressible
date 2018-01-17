!==============================================================================!
subroutine CalcDens(dens, update_rho, dens_0, dens_1)
!------------------------------------------------------------------------------!
!   Discretizes and solves continuity equation                                 !
!   to produce an estimation for rho^(n+1) before CalcZmix on ini=1            !
!   Apply B.C. otherwise                                                       !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use allp_mod, only : unknown, buffer, inflow, outflow, pressure, wall
  use all_mod
  use pro_mod
  use Moin_Problem_mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Parameters]---------------------------------!
  type(unknown)              :: dens
  real,          intent (in) :: dens_0, dens_1 ! dens_0 >= dens_1
  logical,       intent (in) :: update_rho     ! update rho from continuity eq.?
!-----------------------------------[Locals]-----------------------------------!
  integer :: s, c1, c2
  real    :: rhos
!---------------------------------[Interface]----------------------------------!
  ! these function must be pushed into par_mod, then interface is not needed
  interface
    subroutine exchng(phi)
      use allp_mod, only : maxpro
      use all_mod
      use par_mod
      use pro_mod
      implicit none
      include 'mpif.h'
      real    :: phi(-nbc:nc)
      integer :: c1, c2, sub, rtag, stag, length, error
      integer :: status(mpi_status_size)
    end subroutine exchng
  end interface
!==============================================================================!

  !------------------------------------------------!
  !   Estimate rho on n+1 timestep only in ini=1   !
  !   from continuity equation                     !
  !------------------------------------------------!

UpdateDensityBlock: &
  if ( ini == 1 .and. update_rho ) then

    ! b = - flux
    b(:) = 0.0

    !----- old values (o and oo)
    dens % oo(1:NC) = dens % o(1:NC)
    dens % o (1:NC) = dens % n(1:NC)

    do s = 1, NS

      c1 = SideC(1,s)
      c2 = SideC(2,s)

      ! volume integral of div(rho*u)
      ! c1->c2 is considered as a positive direction

      if ( c2 > 0 .or. (c2 < 0 .and. TypeBC(c2) == BUFFER) ) then
        b(c1) = b(c1) - Flux(s)
        if (c2 > 0) b(c2) = b(c2) + Flux(s)
        !elseif(c2 < 0 .and. TypeBC(c2) == BUFFER) then
        else
          b(c1) = b(c1) - Flux(s)
      endif

    end do

    ! Temporal discretization
    if ( INERT == LIN ) then !----- Two time levels; Linear interpolation
      dens % n(1: NC) = dens % o(1: NC) + b(1: NC)*dt/volume(1: NC)
    elseif ( INERT == PAR ) then !----- Three time levels; parabolic interpolation
      dens % n(1: NC) = 4.0/3.0*dens % o(1: NC) - 1.0/3.0*dens % oo(1: NC) &
                     + 2.0/3.0*b(1: NC)*dt/volume(1: NC)
    end if

    !---- Density limiter
    dens  % n(1: NC) = max(dens % n(1: NC) , dens_1)
    dens  % n(1: NC) = min(dens % n(1: NC) , dens_0)

    U     % n(1: NC) = Rho % n(1: NC) * U    % n(1: NC) / dens % n(1: NC)
    V     % n(1: NC) = Rho % n(1: NC) * V    % n(1: NC) / dens % n(1: NC)
    W     % n(1: NC) = Rho % n(1: NC) * W    % n(1: NC) / dens % n(1: NC)
    Zmix  % n(1: NC) = Rho % n(1: NC) * Zmix % n(1: NC) / dens % n(1: NC)

    do s = 1, NS

      c1 = SideC(1,s)
      c2 = SideC(2,s)

      Flux(s) =   f(s) * ( dens % n(c1) * U % n(c1) * Sx(s) +   &
                           dens % n(c1) * V % n(c1) * Sy(s) +   &
                           dens % n(c1) * W % n(c1) * Sz(s) ) + &
          ( 1.0-f(s) ) * ( dens % n(c2) * U % n(c2) * Sx(s) +   &
                           dens % n(c2) * V % n(c2) * Sy(s) +   &
                           dens % n(c2) * W % n(c2) * Sz(s) )

      Flux_u(s) = Flux(s) / rhoS

    end do

  end if UpdateDensityBlock

  !-------------------------!
  !   boundary conditions   !
  !-------------------------!

  ! Since Rho, Rho*Phi - fundamental variables,
  ! then new Phi = (Rho*Phi)_old / Rho_new

BoundaryConditions: &
  do s = 1, NS
    c1 = SideC(1,s)
    c2 = SideC(2,s)

    if ( c2 < 0 .and. TypeBC(c2) == INFLOW ) then
      select case (Moin_Problem())
      case (1)
        ! Dirichlet b.c. for problem 1 (as in Moin)
        Zmix  % n(c2)  = Zmix % n(c1)
        dens  % n(c2)  = dens % n(c1)

        U     %  n(c2) = 0.0
        V     %  n(c2) = 0.0
        W     %  n(c2) = 0.0
  ! this B.C. is correct, but Moin uses another one
  !                case (10)
  !                    ! Dirichlet b.c. for problem 1 (as in Moin)
  !                    dens  % n(c2)    = dens_0
  !
  !                    U     %  n(c2)   = 1.0
  !                    V     %  n(c2)   = 0.0
  !                    W     %  n(c2)   = 0.0
  !                    Zmix  %  n(c2)   = 0.0
      case (10)
        ! Dirichlet b.c. for problem 1 (as in Moin)
        dens  % n(c2)  = dens   % n(c1)

        U     %  n(c2) = 0.0
        V     %  n(c2) = 0.0
        W     %  n(c2) = 0.0
        Zmix  %  n(c2) = Zmix  % n(c1)
      case (2)
        ! Dirichlet b.c. for problem 2 (as in Moin)
        dens  %  n(c2) = 1.0

        U     %  n(c2) = 0.0
        V     %  n(c2) = 0.5
        W     %  n(c2) = 0.0
        Zmix  %  n(c2) = 1.0
      case (4)
        ! Dirichlet b.c. for problem 4_RT
        ! slip condition (face-normal is 0, face-tangential der is zero)
        dens  %  n(c2) = dens % n(c1)
        Zmix  % n(c2)  = Zmix % n(c1) * dens % n(c1) / dens % n(c2)
        U     %  n(c2) = U    %  n(c1)
        V     %  n(c2) = 0.0
        W     %  n(c2) = W    %  n(c1)
      end select

    elseif ( c2 < 0.and.TypeBC(c2) == OUTFLOW ) then
      dens  % n(c2) = dens   % n(c1)
      Zmix  % n(c2) = Rho % n(c1) * Zmix  % n(c1) / dens % n(c2)

      U     % n(c2) = Rho % n(c1) * U     % n(c1) / dens % n(c2)
      V     % n(c2) = Rho % n(c1) * V     % n(c1) / dens % n(c2)
      W     % n(c2) = Rho % n(c1) * W     % n(c1) / dens % n(c2)

    elseif (c2  < 0.and.TypeBC(c2) == WALL ) then
      dens  % n(c2)   = dens  % n(c1)

      U     % n(c2)   = 0.0
      V     % n(c2)   = 0.0
      W     % n(c2)   = 0.0
      Zmix  % n(c2)   = Rho % n(c1) * Zmix  % n(c1) / dens % n(c1)
    elseif (c2  < 0.and.TypeBC(c2) == PRESSURE ) then
      dens  % n(c2)   = dens  % n(c1)

      U     % n(c2)   = Rho % n(c1) * U  % n(c1) / dens % n(c1)
      V     % n(c2)   = 0.0
      W     % n(c2)   = Rho % n(c1) * W  % n(c1) / dens % n(c1)
      Zmix  % n(c2)   = Zmix   % n(c1) * dens % n(c1) / dens % n(c2)
    end if
  end do BoundaryConditions

  !--------------------------------!
  !   Update dependent variables   !
  !--------------------------------!

  ! vectorize this, by setting initial rho to 1.
  U     % n(-NbC: -1) = Rho % n(-NbC: -1) * U    % n(-NbC: -1) / dens % n(-NbC: -1)
  V     % n(-NbC: -1) = Rho % n(-NbC: -1) * V    % n(-NbC: -1) / dens % n(-NbC: -1)
  W     % n(-NbC: -1) = Rho % n(-NbC: -1) * W    % n(-NbC: -1) / dens % n(-NbC: -1)
  Zmix  % n(-NbC: -1) = Rho % n(-NbC: -1) * Zmix % n(-NbC: -1) / dens % n(-NbC: -1)

  !-------------------!
  !   Update Fluxes   !
  !-------------------!

  do s = 1, NS
    c1 = SideC(1,s)
    c2 = SideC(2,s)

    rhoS = f(s) * dens % n(c1) + ( 1.0-f(s) ) * dens % n(c2)

    if ( c2 < 0 .and. TypeBC(c2).ne.BUFFER ) then
      Flux(s) =   f(s) * ( Rho % n(c1) * U % n(c1) * Sx(s) +   &
                           Rho % n(c1) * V % n(c1) * Sy(s) +   &
                           Rho % n(c1) * W % n(c1) * Sz(s) ) + &
          ( 1.0-f(s) ) * ( Rho % n(c2) * U % n(c2) * Sx(s) +   &
                           Rho % n(c2) * V % n(c2) * Sy(s) +   &
                           Rho % n(c2) * W % n(c2) * Sz(s) )
    end if

    Flux_u(s) = Flux(s) / rhoS
  end do

  !---------------------------!
  !   Update Rho at bufffer   !
  !---------------------------!

  call Exchng(dens % n)

  call Exchng(Zmix % n)
  call Exchng(U    % n)
  call Exchng(V    % n)
  call Exchng(W    % n)

end subroutine CalcDens
