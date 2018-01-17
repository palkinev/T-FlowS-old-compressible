!==============================================================================!
subroutine CalcZmix( var,  & ! = 7
                     dZdx, & ! dz/dx
                     dZdy, & ! dz/dy
                     dZdz  ) ! dz/dz
!------------------------------------------------------------------------------!
!  Purpose: Solve transport equation for Rho * Z                               !
!  Equation is described in latex manual (to do)                               !
!  main variables are rho, rho * Zmix                                          !
!------------------------------[Modules]---------------------------------------!
  use allp_mod, only: Unknown, buffer, convect, inflow, outflow, wall, wallfl
  use all_mod
  use pro_mod
  use rans_mod
  use par_mod

  use moin_problem_mod !  module to control conflicting Moin modules
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Parameters]---------------------------------!
  integer,       intent(in) :: var
  real,          intent(in) :: dZdx(-NbC:NC),  dZdy(-NbC:NC),  dZdz(-NbC:NC)
  !-------------------------------[Locals]-------------------------------!
  real    :: rhoZS, dZdxS,  dZdyS,  dZdzS
  integer :: c,s,c1,c2,niter,miter,mat
  real    :: A0,A12,A21,D12,D21,error
  real    :: Fex, Fim
  real    :: RhoS, VISeff
  real, allocatable :: rhoZ(:)
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

    subroutine cg(N, NB, NONZ, A, Acol, Arow, Adia, Ab, x, r1, &
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
    end subroutine cg

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
  allocate(rhoZ (-NbC:NC) )
  rhoZ = Rho % n * Zmix % n


  !----- This is important for "copy" boundary conditions. Find out why !
  Abou(-NbC: -1) = 0.0

  Aval(:) = 0.0
  b(:) = 0.0
  !-----------------------------------------!
  !     Initialize variables and fluxes     !
  !-----------------------------------------!

  !----- old values (o and oo)
  if (ini == 1) then
    Zmix % oo(1: NC) = Zmix % o(1: NC)
    Zmix % o (1: NC) = Zmix % n(1: NC)
  end if


  !====================!
  !                    !
  !     Convection     !
  !                    !
  !====================!

  !----- Compute Zmax and Zmin
  do mat = 1, Nmat
    if(BLEND_TEM(mat) /= NO) then
      call CalMinMax(rhoZ)  ! or Z % o ???
      goto 1
    end if
  end do

  !----- new values
  1 Zmix % C = 0.0
    Zmix % X = 0.0

  !--------------------- calculating Z_SRC coefficient and defining Z_SRC
  ! Z_SRC_of_Z fixed at current timestep
  select case (Moin_Problem())
!        case (10)
!            if (ini == 1) then
!                do c = 1, NC
!                    Zmix % source(c) = Z_SRC_INT * moin_problem_z_src (Zmix % n(c), 0.0, 0.0) ! Z_SRC_of_Z
!                end do
!            end if
    case default
      do c = 1, NC
        Zmix % source(c) = moin_problem_z_src (xc(c), yc(c), Moin_Time())     ! Z_SRC_of_x_y_t
      end do
  end select

  !--------------------------------!
  !     Spatial Discretization     !
  !--------------------------------!
  do s = 1, NS

    c1 = SideC(1,s)
    c2 = SideC(2,s)

    !---- Central differencing
    rhoZS = f(s)*rhoZ(c1) + (1.0-f(s))*rhoZ(c2)

    if(BLEND_TEM(material(c1)) /= NO .or. BLEND_TEM(material(c2)) /= NO) then
      call ConvScheme(rhoZS, s, rhoZ, dZdx, dZdy, dZdz, Dx, Dy, Dz, &
        max(BLEND_TEM(material(c1)),BLEND_TEM(material(c2))) )
    end if

    if(c2 > 0) then
      Zmix % C(c1) = Zmix % C(c1) - Flux_u(s)*rhoZS
      Zmix % C(c2) = Zmix % C(c2) + Flux_u(s)*rhoZS
    else
      Zmix % C(c1) = Zmix % C(c1) - Flux_u(s)*rhoZS
    endif

    !---- Upwind
    if(BLEND_TEM(material(c1)) /= NO .or. BLEND_TEM(material(c2)) /= NO) then
      if(Flux_u(s)<0) then   ! from c2 to c1
        Zmix % X(c1)   = Zmix % X(c1) - Flux_u(s)*rhoZ(c2)
        if(c2 > 0) then
          Zmix % X(c2) = Zmix % X(c2) + Flux_u(s)*rhoZ(c2)
        endif
      else !from c1 to c2
        Zmix % X(c1)   = Zmix % X(c1) - Flux_u(s)*rhoZ(c1)
        if(c2 > 0) then
          Zmix % X(c2) = Zmix % X(c2) + Flux_u(s)*rhoZ(c1)
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
  if(CONVEC.eq.FI) then
    do c = 1, NC
      b(c) = b(c) + URFC_Tem(material(c)) * &
        (Zmix % C(c) - Zmix % X(c))
            !( Zmix % C(c)) !CDS right side
            !(Zmix % X(c)) !UDS right side
    end do
  end if

  !==================!
  !                  !
  !     Diffusion     !
  !                  !
  !==================!

  !----- Set Z % X back to zero
  do c = 1, NC
    Zmix % X(c) = 0.0
  end do

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

    RhoS = fF(s)*Rho % n(c1) + (1.0-fF(s))*Rho % n(c2)

    select case (Moin_Problem())
      case (4) ! RT
        VISeff = VIScZmix * RhoS
!            case (10)
!                VISeff = fF(s)*VIScZmix_Dyn(c1) + (1.0-fF(s))*VIScZmix_Dyn(c2)
!                !VISeff = (fF(s)*VIScZmix_Dyn(c1)*Rho % n(c1) + (1.0-fF(s))*VIScZmix_Dyn(c2)*Rho % n(c2))
      case default
        VISeff = VIScZmix
    end select

    !----- gradients on the cell fFace
    dZdxS = fF(s)*dZdx(c1) + (1.0-fF(s))*dZdx(c2)
    dZdyS = fF(s)*dZdy(c1) + (1.0-fF(s))*dZdy(c2)
    dZdzS = fF(s)*dZdz(c1) + (1.0-fF(s))*dZdz(c2)

    !---- total (exact) diffusive flux
    Fex = VISeff*( dZdxS*Sx(s)+dZdyS*Sy(s)+dZdzS*Sz(s) )

    !---- implicit diffusive flux
    A0 = VISeff * Scoef(s)

    Fim = A0 *        &
      ( dZdxS*Dx(s)   &
      + dZdyS*Dy(s)   &
      + dZdzS*Dz(s) )

    !---- cross diffusion part
      Zmix % X(c1) = Zmix % X(c1) + Fex - Fim
    if(c2 > 0) then
      Zmix % X(c2) = Zmix % X(c2) - Fex + Fim
    end if


    !----- calculate the coefficients for the sysytem matrix
    if( (DIFFUS.eq.CN) .or. (DIFFUS.eq.FI) ) then

      ! Crank Nicholson

      if(DIFFUS .eq. FI) then       ! Fully implicit
        D12 = D12 + A0 / Rho % n(c2)
        D21 = D21 + A0 / Rho % n(c1)
      end if

      !convective term UDS
      if(BLEND_TEM(material(c1)) /= NO .or. BLEND_TEM(material(c2)) /= NO) then
        A12 = A12 - min(Flux_u(s), 0.0)
        A21 = A21 + max(Flux_u(s), 0.0)
      endif

      !----- fill the system matrix
      if(c2 > 0) then
        !m1
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
        if( &
          (TypeBC(c2) == INFLOW)  .or.     &
          (TypeBC(c2) == WALL)    .or.     &
          (TypeBC(c2) == CONVECT) .or.     &
          !(TypeBC(c2) == OUTFLOW) .or.     &
          (TypeBC(c2) == WALLFL) ) then
          !m1
          Aval(Adia(c1)) = Aval(Adia(c1)) + A21 + D21
          b(c1) = b(c1) + ( A12 + D12 ) * rhoZ(c2)
          !b(c1) = b(c1) + ( A21 + D21 ) * rhoZ(c2)
          !m3
          !Aval(Adia(c1)) = Aval(Adia(c1)) + A12 + D21
          !b(c1) = b(c1) + ( A12 + D21 ) * rhoZ(c2)
        elseif( TypeBC(c2) == OUTFLOW)  then
          !m1
          Aval(Adia(c1)) = Aval(Adia(c1)) + A21
          b(c1) = b(c1) + ( A12 ) * rhoZ(c2)
          !m3
          !Aval(Adia(c1)) = Aval(Adia(c1)) + A12
          !b(c1) = b(c1) + ( A12 ) * Uri % n(c2)
        else if( TypeBC(c2) == BUFFER ) then
          !m1
          Aval(Adia(c1)) = Aval(Adia(c1)) + A21 + D21
          Abou(c2) = - A12 - D12  ! cool parallel stuff
          !Abou(c2) = - A21 - D21  ! cool parallel stuff
          !m3
          !Aval(Adia(c1)) = Aval(Adia(c1)) + A12 + D21
          !Abou(c2) = - A12 - D21  ! cool parallel stuff
        endif
      end if

    end if

  end do  ! through sides

  !---------------------------------*
  !     Temporal discretization     *
  !---------------------------------*

  !----- Adams-Bashfort scheeme for diffusion fluxes

  !----- Crank-Nicholson scheme for difusive terms

  !----- Fully implicit treatment for diffusive terms
  !      is handled via the linear system of equations

  !----- Adams-Bashfort scheeme for cross diffusion

  !----- Crank-Nicholson scheme for cross difusive terms

  !----- Fully implicit treatment for cross difusive terms
  if(CROSS.eq.FI) then
    b(1: NC) = b(1: NC) + Zmix % X(1: NC)
  end if
  !========================*
  !                        *
  !     Inertial terms     *
  !                        *
  !========================*

  !----- Two time levels; Linear interpolation
  if(INERT.eq.LIN) then
    Aval(Adia(1: NC)) = Aval(Adia(1: NC)) + volume(1: NC)/dt
    b(1: NC)  = b(1: NC) + volume(1: NC)/dt * Zmix % o(1: NC)*Rho % o(1: NC) + Zmix % source(1: NC)*volume(1: NC)!*vd2d_mms_rho (zc(1: NC), dt)  ! remove the source when tested
  end if

  !----- Three time levels; parabolic interpolation
  if(INERT.eq.PAR) then
    Aval(Adia(1: NC)) = Aval(Adia(1: NC)) + 1.5 * volume(1: NC)/dt
    b(1: NC)  = b(1: NC) + 2.0 * volume(1: NC)/dt * Rho % o(1: NC)*Zmix % o(1: NC) - 0.5 * volume(1: NC)/dt * Rho % oo(1: NC)*Zmix % oo(1: NC) + Zmix % source(1: NC)*volume(1: NC) ! remove the source when tested
  end if

  !===================================!
  !                                   !
  !     Solve the equations for Z     !
  !                                   !
  !===================================!

  b(1:NC) = b(1:NC) + Aval(Adia(1:NC)) * (1.0 - Zmix % URF)*rhoZ(1:NC) / Zmix % URF !?? * Z % n(1:NC)
  Aval(Adia(1:NC)) = Aval(Adia(1:NC)) / Zmix % URF

  if(ALGOR == SIMPLE)   miter = 20
  if(ALGOR == FRACT)    miter = 5

  niter = miter

  call cg( NC, Nbc, NONZERO, Aval, Acol, Arow, Adia, Abou, &
       rhoZ, b, PREC, niter, Zmix % STol, res(var), error )

  write(LineRes(65:76),  '(1PE12.3)') res(var)
  write(LineRes(93:96),  '(I4)')      niter
  write(LineRes1(65:76), '(1PE12.3)') error

  Zmix  % n(1: NC) = rhoZ(1: NC) / Rho % n(1: NC)

  ! limiter for Zmix < 0 and > 1 after CalcZmix
  Zmix  % n(1: NC) = max(Zmix % n(1: NC) , 0.0)
  Zmix  % n(1: NC) = min(Zmix % n(1: NC) , 1.0)

  call Exchng(Zmix  % n)
  call wait


  if (allocated(rhoZ) ) deallocate (rhoZ)

end subroutine CalcZmix