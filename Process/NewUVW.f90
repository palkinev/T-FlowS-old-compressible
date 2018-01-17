!==============================================================================!
subroutine newuvw( var,                 & ! 1 - x, 2 - y, 3 - z
                   Ui,                  & ! U_i
                   dUidi, dUidj, dUidk, &
                   dUjdi, dUkdi,        &
                   dUjdj, dUkdk,        & !Sij conponents
                   Si,    Sj,    Sk,    &
                   Di,    Dj,    Dk,    &
                   Pi )
                   ! i - main component, j and k are secondary
                   ! this means:
                   ! var = 1 -> Ui = U -> i = 1, j = 2, k = 3
                   ! var = 2 -> Ui = V -> i = 2, j = 1, k = 3
                   ! var = 3 -> Ui = W -> i = 3, j = 1, k = 2
!------------------------------------------------------------------------------!
!   Discretizes and solves momentum conservation equations for (rho*u_i)       !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use allp_mod, only: buffer, convect, inflow, outflow, wall, wallfl, unknown
  use all_mod
  use pro_mod
  use les_mod
  use rans_mod
  use par_mod
  use Moin_Problem_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Parameters]---------------------------------!
  type(unknown)             :: Ui
  integer,       intent(in) :: var
  real,          intent(in) :: dUidi(-NbC: NC), &
                               dUidj(-NbC: NC), &
                               dUidk(-NbC: NC), &
                               dUjdi(-NbC: NC), &
                               dUkdi(-NbC: NC), &
                               dUjdj(-NbC: NC), &
                               dUkdk(-NbC: NC)

  real,          intent(in) :: Si(NS), Sj(NS), Sk(NS)
  real,          intent(in) :: Di(NS), Dj(NS), Dk(NS)
  real,          intent(in) :: Pi(-NbC:NC)
!-----------------------------------[Locals]-----------------------------------!
  integer :: c, c1, c2, s, niter, miter, mat
  real    :: Fex, Fim
  real    :: rhoUiS
  real    :: A0, A12, A21, D12, D21
  real    :: error
  real    :: VISeff
  real    :: dUidiS,dUidjS,dUidkS,dUjdiS,dUkdiS,dUjdjS,dUkdkS
  real    :: RhoS
  real, allocatable :: rhoUi(:)
!------------------------------------------------------------------------------!
!
!     [A]{u} = {b}   [kgm/s^2]   [N]
!
!  Dimensions of certain variables
!
!     A              [kg/s]
!     U, V, W        [m/s]
!     bU, bV, bW     [kgm/s^2]   [N]
!     P, PP          [kg/ms^2]   [N/m^2]
!     Flux           [kg/s]
!     CU*, CV*, CW*  [kgm/s^2]   [N]
!     DU*, DV*, DW*  [kgm/s^2]   [N]
!     XU*, XV*, XW*  [kgm/s^2]   [N]
!
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

    subroutine wait
      implicit none
      include 'mpif.h'
      integer :: error
    end subroutine wait

    subroutine bicg(N, NB, NONZ, A, Acol, Arow, Adia, Ab, x, r1, &
      prec, niter, tol, IniRes, FinRes)
      use allp_mod, only: tiny
      use sol_mod
      use par_mod
      implicit none
      integer  :: n, nb, nonz

      real     :: a(nonz),ab(-nb:-1)
      integer  :: acol(n+1),adia(n)
      integer  :: arow(nonz)
      real     :: x(-nb:n), r1(n)

      integer  :: prec,  niter
      real     :: tol
      real     :: inires, finres
    end subroutine bicg

    subroutine calminmax(phi)
      use allp_mod, only : buffer
      use all_mod
      use pro_mod
      implicit none
      real    :: phi(-nbc:nc)
      integer :: c1, c2, s
    end subroutine calminmax

    subroutine convscheme(phis, s,                            &
                          phi,                                &
                          dphidi, dphidj, dphidk, di, dj, dk, &
                          blenda)
      use all_mod
      use pro_mod
      implicit none
      real          :: phis
      integer       :: s
      real          :: phi(-nbc:nc)
      real          :: dphidi(-nbc:nc), dphidj(-nbc:nc), dphidk(-nbc:nc)
      real          :: di(ns),          dj(ns),          dk(ns)
      integer       :: blenda
      integer       :: c1, c2, c, d
      real          :: fj
      real          :: gd, gu, alfa, beta1, beta2
      real          :: phij, phiu, phiustar, rj, sign, gammac, beta
    end subroutine convscheme
  end interface
!==============================================================================!

  ! equation is solved for this variable
  allocate(rhoUi (-NbC:NC) );
  rhoUi = Rho % n * Ui % n

! currently only works with
!    CROSS=FI
!    DIFFUS=FI
!    CONVEC=FI

  !----- This is important for "copy" boundary conditions. Find out why !
  Abou(-NbC: -1) = 0.0

  Aval(:) = 0.0
  b(:) = 0.0
  !-----------------------------------------!
  !     Initialize variables and fluxes     !
  !-----------------------------------------!

  !----- old values (o and oo)
  if(ini == 1) then
    Ui % oo(1: NC) = Ui % o(1: NC)
    Ui % o (1: NC) = Ui % n(1: NC)
  end if

  !====================!
  !                    !
  !     Convection     !
  !                    !
  !====================!

  !----- Compute PHImax and PHImin
  do mat = 1, Nmat
    if(BLEND(mat) /= NO) then
      call CalMinMax(rhoUi)  ! or Ui % o ???
      goto 1
    end if
  end do

  !----- new values
  1 Ui % C = 0.0
    Ui % X = 0.0


  !----- move this in some othe function
  select case (Moin_Problem())
    case (4) ! RT
      select case (var)
        case (2)
          do c = 1, NC
            Ui % source(c) = Rho % n(c)*Moin_Problem_v_src (xc(c), yc(c), Moin_Time())
          end do
      end select
    case default
      select case (var)
        case (1)
          do c = 1, NC
            Ui % source(c) = Moin_Problem_u_src (xc(c), yc(c), Moin_Time())
          end do
        case (2)
          do c = 1, NC
            Ui % source(c) = Moin_Problem_v_src (xc(c), yc(c), Moin_Time())
          end do
      end select
  end select

  !--------------------------------!
  !     Spatial Discretization     !
  !--------------------------------!

  do s = 1, NS

    c1 = SideC(1,s)
    c2 = SideC(2,s)

    !---- Central differencing
    rhoUiS = f(s)*rhoUi(c1) + (1.0-f(s))*rhoUi(c2)

    if(BLEND(material(c1)) /= NO .or. BLEND(material(c2)) /= NO) then
      call ConvScheme(rhoUiS, s, rhoUi, dUidi, dUidj, dUidk, Di, Dj, Dk, &
        max(BLEND(material(c1)),BLEND(material(c2))) )
    end if

    !---- Central differencing for convection
    if(c2 > 0) then
      Ui % C(c1) = Ui % C(c1) - Flux_u(s)*rhoUiS
      Ui % C(c2) = Ui % C(c2) + Flux_u(s)*rhoUiS
    else
      Ui % C(c1) = Ui % C(c1) - Flux_u(s)*rhoUiS
    endif

    !---- Upwind
    if(BLEND_TEM(material(c1)) /= NO .or. BLEND_TEM(material(c2)) /= NO) then
      if(Flux_u(s) < 0) then   ! from c2 to c1
        Ui % X(c1)   = Ui % X(c1) - Flux_u(s)*rhoUi(c2)
        if(c2 > 0) then
          Ui % X(c2) = Ui % X(c2) + Flux_u(s)*rhoUi(c2)
        endif
      else !from c1 to c2
        Ui % X(c1)   = Ui % X(c1) - Flux_u(s)*rhoUi(c1)
        if(c2 > 0) then
          Ui % X(c2) = Ui % X(c2) + Flux_u(s)*rhoUi(c1)
        endif
      end if
    end if   ! BLEND_TEM
  end do  ! through sides

  
  !---------------------------------!
  !     Temporal discretization     !
  !---------------------------------!

  !----- Adams-Bashforth scheeme for convective fluxes

  !----- Crank-Nicholson scheeme for convective fluxes

  !----- Fully implicit treatment of convective fluxes
  if(CONVEC == FI) then
    do c = 1, NC
      b(c) = b(c) + URFC(material(c)) * &
        (Ui % C(c) - Ui % X(c))
    end do
  end if
  
  !----------------------------------------------------!
  !     Browse through all the faces, where else ?     !
  !----------------------------------------------------!

  !----- new values
  do c = 1, NC
    Ui % X(c) = 0.0
  end do

  !==================!
  !                  !
  !     Diffusion    !
  !                  !
  !==================!

  !--------------------------------!
  !     Spatial Discretization     !
  !--------------------------------!
  do s = 1, NS
    A12 = 0.0
    A21 = 0.0
    D12 = 0.0
    D21 = 0.0
    c1 = SideC(1,s)
    c2 = SideC(2,s)
  
    RhoS = fF(s)*Rho % n(c1) + (1.0-fF(s))*Rho % n(c2) !check this !!!!

    select case (Moin_Problem())
      case (4) ! RT
        VISeff = VISc * RhoS
!            case (10)
!                VISeff = (fF(s)*VISc_Dyn(c1) + (1.0-fF(s))*VISc_Dyn(c2))
!                !VISeff = (fF(s)*VISc_Dyn(c1)*Rho % n(c1) + (1.0-fF(s))*VISc_Dyn(c2)*Rho % n(c2))
      case default
        VISeff = VISc
    end select


    !-------------CROSS diffusive flux-------------------------------
    !Sij 1
    dUidiS = fF(s)*dUidi(c1) + (1.0-fF(s))*dUidi(c2)
    dUidjS = fF(s)*dUidj(c1) + (1.0-fF(s))*dUidj(c2)
    dUidkS = fF(s)*dUidk(c1) + (1.0-fF(s))*dUidk(c2)
    !Sij 2
    dUjdiS = fF(s)*dUjdi(c1) + (1.0-fF(s))*dUjdi(c2)
    dUkdiS = fF(s)*dUkdi(c1) + (1.0-fF(s))*dUkdi(c2)
    !Sij 3
    dUjdjS = fF(s)*dUjdj(c1) + (1.0-fF(s))*dUjdj(c2)
    dUkdkS = fF(s)*dUkdk(c1) + (1.0-fF(s))*dUkdk(c2)

    !---- total exact diffusive flux ( Sij = 1/2 dUi/dxj + 1/2 dUj/dxi - 1/3 deltaij dUk/dxk)
    !----                                        (1)            (2)              (3)
    Fex = VISeff*(4.0/3.0*  dUidiS          *Si(s) & ! 2*Sij (1)
                         + (dUidjS + dUjdiS)*Sj(s) & ! 2*Sij (2)
                         + (dUidkS + dUkdiS)*Sk(s) & ! 2*Sij (2)
               - 2.0/3.0 * (dUjdjS + dUkdkS)*Si(s) ) ! 2*Sij (3)

    ! this total exact is without (3)
    !Fex=VISeff*(     dUidiS *Si(s)  & ! 2*Sij 1
    !             +   dUjdiS *Sj(s)  & ! 2*Sij 2
    !             +   dUkdiS *Sk(s)  )

    !---- implicit diffusive flux of Sij (1)
    A0 = VISeff * Scoef(s)

    Fim = A0*( dUidiS * Di(s) &
             + dUidjS * Dj(s) &
             + dUidkS * Dk(s) )

    !Fim=A0*    (4.0/3.0* dUidiS        *Di(s)  & ! 2*Sij 1
    !                 +  (dUidjS+dUjdiS)*Dj(s)  & ! 2*Sij 2
    !                 +  (dUidkS+dUkdiS)*Dk(s)  &
    !          - 2.0/3.0*(dUjdjS+dUkdkS)*Di(s)  ) ! 2*Sij 3


    !-------------cross diffusive flux-------------------------------

    !---- cross diffusion part
    Ui % X(c1) = Ui % X(c1) + Fex - Fim
    if(c2  > 0) then
      Ui % X(c2) = Ui % X(c2) - Fex + Fim
    end if

    !----- calculate the coefficients for the system matrix
    if( (DIFFUS == CN) .or. (DIFFUS == FI) ) then

      !if(DIFFUS  ==  CN) then       ! Crank Nicholson
      !  A12 = 0.5 * A0
      !  A21 = 0.5 * A0
      !end if

      if(DIFFUS  ==  FI) then       ! Fully implicit
        D12 = D12 + A0 / Rho % n(c2)
        D21 = D21 + A0 / Rho % n(c1)
      end if

      if(BLEND(material(c1)) /= NO .or. BLEND(material(c2)) /= NO) then
        A12 = A12 - min(Flux_u(s), real(0.0))
        A21 = A21 + max(Flux_u(s), real(0.0))
      endif

      !----- fill the system matrix
      if(c2  > 0) then
        !m1 convective tests - yes, diffusive tests - yes
        Aval(Adia(c1))    = Aval(Adia(c1))    + A21 + D21
        Aval(SidAij(2,s)) = Aval(SidAij(2,s)) - A21 - D21
        Aval(Adia(c2))    = Aval(Adia(c2))    + A12 + D12
        Aval(SidAij(1,s)) = Aval(SidAij(1,s)) - A12 - D12

        !m2 (old code) convective tests - yes, diffusive tests - no
        !Aval(SidAij(1,s)) = Aval(SidAij(1,s)) - A12 - D12
        !Aval(Adia(c1))    = Aval(Adia(c1))    + A12 + D12
        !Aval(SidAij(2,s)) = Aval(SidAij(2,s)) - A21 - D21
        !Aval(Adia(c2))    = Aval(Adia(c2))    + A21 + D21

        ! it seems m2 has wrong diffusive part, but acceptable convective parts, therefore
        ! m3 convective tests - 1st order - bad
        !Aval(SidAij(1,s)) = Aval(SidAij(1,s)) - A12 - D12
        !Aval(Adia(c1))    = Aval(Adia(c1))    + A12 + D21
        !Aval(SidAij(2,s)) = Aval(SidAij(2,s)) - A21 - D21
        !Aval(Adia(c2))    = Aval(Adia(c2))    + A21 + D12
      else if(c2  < 0) then
        !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
        ! Outflow is not included because it was causing problems     !
        !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
        if( (TypeBC(c2) == WALL)    .or.     &
          (TypeBC(c2) == CONVECT) .or.     &
          (TypeBC(c2) == INFLOW)  .or.     &
          !(TypeBC(c2) == OUTFLOW) .or.     &
          (TypeBC(c2) == WALLFL) ) then
          !m1
          Aval(Adia(c1)) = Aval(Adia(c1)) + A21 + D21
          b(c1) = b(c1) + ( A12 + D12 ) * rhoUi(c2)
          !m3
          !Aval(Adia(c1)) = Aval(Adia(c1)) + A12 + D21
          !b(c1) = b(c1) + ( A12 + D21 ) * rhoUi(c2)
        elseif( TypeBC(c2) == OUTFLOW)  then
          !m1
          Aval(Adia(c1)) = Aval(Adia(c1)) + A21
          b(c1) = b(c1) + ( A12 ) * rhoUi(c2)
          !m3
          !Aval(Adia(c1)) = Aval(Adia(c1)) + A12
          !b(c1) = b(c1) + ( A12 ) * rhoUi(c2)
        else if( TypeBC(c2) == BUFFER ) then
          !m1
          Aval(Adia(c1)) = Aval(Adia(c1)) + A21 + D21
          Abou(c2) = - A12 - D12  ! cool parallel stuff

          !m3
          !Aval(Adia(c1)) = Aval(Adia(c1)) + A12 + D21
          !Abou(c2) = - A21 - D12  ! cool parallel stuff
        endif
      end if

    end if

  end do  ! through sides

  !---------------------------------!
  !     Temporal discretization     !
  !---------------------------------!

  !----- Adams-Bashfort scheeme for diffusion fluxes

  !----- Crank-Nicholson scheme for difusive terms

  !----- Fully implicit treatment for diffusive terms
  !      is handled via the linear system of equations

  !----- Adams-Bashfort scheeme for cross diffusion

  !----- Crank-Nicholson scheme for cross diffusive terms

  !----- Fully implicit treatment for cross diffusive terms
  if(CROSS == FI) then
    do c = 1, NC
      b(c) = b(c) + Ui % X(c)
    end do
  end if

  !========================!
  !                        !
  !     Inertial terms     !
  !                        !
  !========================!

  !----- Two time levels; Linear interpolation
  if(INERT.eq.LIN) then
    Aval(Adia(1: NC)) = Aval(Adia(1: NC)) + volume(1: NC)/dt
    b(1: NC) = b(1: NC) + volume(1: NC)/dt * Rho % o(1: NC)*Ui % o(1: NC) + Ui % source(1: NC)*volume(1: NC)!*vd2d_mms_rho (zc(1: NC), dt)  ! remove the source when tested
  end if

  !----- Three time levels; parabolic interpolation
  if(INERT.eq.PAR) then
    Aval(Adia(1: NC)) = Aval(Adia(1: NC)) + 1.5 * volume(1: NC)/dt
    b(1: NC) = b(1: NC) + 2.0 * volume(1: NC)/dt * Rho % o(1: NC)*Ui % o(1: NC) - 0.5 * volume(1: NC)/dt * Rho % oo(1: NC)*Ui % oo(1: NC) + Ui % source(1: NC)*volume(1: NC) ! remove the source when tested
  end if

  !=====================================!
  !                                     !
  !     Pressure term contributions     !
  !                                     !
  !=====================================!

  !------------------------------!
  !     Global pressure drop     !
  !------------------------------!
  !  if(var == 1) then
  !    do c = 1, NC
  !      b(c) = b(c)  + PdropX(material(c)) * volume(c)
  !    end do
  !  else if(var == 2) then
  !    do c = 1, NC
  !      b(c) = b(c)  + PdropY(material(c)) * volume(c)
  !    end do
  !  else if(var == 3) then
  !    do c = 1, NC
  !      b(c) = b(c)  + PdropZ(material(c)) * volume(c)
  !    end do
  !  end if

  !-------------------------------------!
  !     Local pressure distribution     !
  !-------------------------------------!
  b(1: NC) = b(1: NC) - Pi(1: NC)*volume(1: NC)
  !--------------------------------------------!
  !     All other terms defined by the user    !
  !--------------------------------------------!
  !!!  if(HOT == YES) call UserForce(var)


  !=======================================!
  !                                       !
  !     Solve the equations for U,V,W     !
  !                                       !
  !=======================================!

  if (var == 1) Asave(1:NC) = Aval(Adia(1:NC)) ! diagonal values before relaxation, used in CalcPS
  b(1:NC) = b(1:NC) + Aval(Adia(1:NC)) * (1.0 - U % URF)*rhoUi(1:NC) / U % URF
  Aval(Adia(1:NC)) = Aval(Adia(1:NC)) / U % URF

  if(ALGOR == SIMPLE)   miter = 10
  if(ALGOR == FRACT)    miter = 5

  niter = 200

  call bicg(NC, Nbc, NONZERO, Aval,Acol,Arow,Adia,Abou, &
            rhoUi, b, PREC,                             &
            niter, Ui % STol, res(var), error)


  if(var == 1) then
    write(LineRes(17:28), '(1PE12.3)')  res(var)
    write(LineRes(77:80), '(I4)')       niter
    write(LineRes1(17:28), '(1PE12.3)') error
  end if
  if(var == 2) then
    write(LineRes(29:40), '(1PE12.3)')  res(var)
    write(LineRes(81:84), '(I4)')       niter
    write(LineRes1(29:40), '(1PE12.3)') error
  end if
  if(var == 3) then
    write(LineRes(41:52), '(1PE12.3)')  res(var)
    write(LineRes(85:88), '(I4)')       niter
    write(LineRes1(41:52), '(1PE12.3)') error
  end if

  Ui  % n(1: NC) = rhoUi(1: NC) / Rho % n(1: NC)

  call Exchng(Ui  % n)
  call wait

  if (allocated(rhoUi) ) deallocate (rhoUi)


end subroutine newuvw
