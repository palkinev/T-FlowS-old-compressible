!==============================================================================!
  subroutine Calc_Sgs_Coefficients_Dynamic_Smagorinsky
!------------------------------------------------------------------------------!
!   Calculates coefficients in Smagorinsky model                               !
!   with dynamic procedure for compressible case.                              !
!   This function follows the manual in latex called "???"                     !
!------------------------------------------------------------------------------!
!   Attention: LaTex formalism is used in this function                        !
!------------------------------------------------------------------------------!
!   To do: check on inflow and outflow                                         !
!------------------------------------------------------------------------------!
!   All variables considered in LES code are spatially filtered by cells       !
!   themselves and denoted as vars with bar above them:                        !
!   and denoted as vars with bar:                                              !
!   \bar{rho}, \bar{U}, \bar{V},_\bar{W}, \bar{Zmix}.                          !
!   This procedure uses spectral information on the flow, thus larger test     !
!   filter is introduced and denoted as \hat{} (^), it is applied after cell   !
!   filter \bar{}.                                                             !
!----------------------------------[Modules]-----------------------------------!
  use allp_mod, only : TINY, BUFFER
  use all_mod
  use pro_mod
  use les_mod
  use rans_mod
  use par_mod
!------------------------------------------------------------------------------!
  implicit none
!-----------------------------------[Locals]-----------------------------------!
  integer           :: s, c1, c2
  real, allocatable :: hat_volume(:), hat_delta(:), bar_delta(:)
  real, allocatable :: hat_rho(:)
  real, allocatable :: tilde_u_x(:), tilde_u_y(:), tilde_u_z(:), &
                       tilde_v_x(:), tilde_v_y(:), tilde_v_z(:), &
                       tilde_w_x(:), tilde_w_y(:), tilde_w_z(:)
  real, allocatable :: tilde_sij_11(:), tilde_sij_22(:), tilde_sij_33(:), &
                       tilde_sij_12(:), tilde_sij_13(:), tilde_sij_23(:), &
                       tilde_sij_kk(:)
  real, allocatable :: pi_g(:), pi_t(:)
  real, allocatable :: hat_rui_1(:), hat_rui_2(:), hat_rui_3(:)
  real, allocatable :: hat_r_tilde_sij_11(:), hat_r_tilde_sij_22(:), &
                       hat_r_tilde_sij_33(:), hat_r_tilde_sij_12(:), &
                       hat_r_tilde_sij_13(:), hat_r_tilde_sij_23(:), &
                       hat_r_tilde_sij_kk(:)
  real, allocatable :: hat_ruiuj_11(:), hat_ruiuj_22(:), hat_ruiuj_33(:), &
                       hat_ruiuj_12(:), hat_ruiuj_13(:), hat_ruiuj_23(:)
  real, allocatable :: hat_p_term(:)
  real, allocatable :: hat_mij_term_11(:), hat_mij_term_22(:), &
                       hat_mij_term_33(:), &
                       hat_mij_term_12(:), hat_mij_term_13(:), &
                       hat_mij_term_23(:)
  real, allocatable :: lij_11(:), lij_22(:), lij_33(:), &
                       lij_12(:), lij_13(:), lij_23(:), &
                       lij_kk(:)
  real, allocatable :: mij_11(:), mij_22(:), mij_33(:), &
                       mij_12(:), mij_13(:), mij_23(:), &
                       mij_kk(:), mijmij(:)
  real, allocatable :: mijld_ij(:)
  real, allocatable :: c_i(:), c_r(:)
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

    subroutine graphi(phi, i, phii, boundary)
      use allp_mod, only: symmetry, buffer
      use all_mod
      use pro_mod
      implicit none
      integer :: i
      real    :: phi(-nbc:nc), phii(-nbc:nc)
      logical :: boundary
      integer :: s, c, c1, c2
      real    :: dphi1, dphi2, dxc1, dyc1, dzc1, dxc2, dyc2, dzc2
    end subroutine graphi
  end interface
!==============================================================================!

  ! I do not yet know why we do this
  call Exchng( U % n  )
  call Exchng( V % n  )
  call Exchng( W % n  )

  !---------------------!
  !   allocating vars   !
  !---------------------!

  allocate( tilde_u_x(-NbC:NC)       ); tilde_u_x          = 0.
  allocate( tilde_u_y(-NbC:NC)       ); tilde_u_y          = 0.
  allocate( tilde_u_z(-NbC:NC)       ); tilde_u_z          = 0.
  allocate( tilde_v_x(-NbC:NC)       ); tilde_v_x          = 0.
  allocate( tilde_v_y(-NbC:NC)       ); tilde_v_y          = 0.
  allocate( tilde_v_z(-NbC:NC)       ); tilde_v_z          = 0.
  allocate( tilde_w_x(-NbC:NC)       ); tilde_w_x          = 0.
  allocate( tilde_w_y(-NbC:NC)       ); tilde_w_y          = 0.
  allocate( tilde_w_z(-NbC:NC)       ); tilde_w_z          = 0.
  allocate( bar_delta(-NbC:NC)       ); bar_delta          = 0.
  allocate( hat_volume(1:NC)         ); hat_volume         = 0.
  allocate( hat_delta(1:NC)          ); hat_delta          = 0.
  allocate( hat_rho(1:NC)            ); hat_rho            = 0.
  allocate( tilde_sij_11(-NbC:NC)    ); tilde_sij_11       = 0.
  allocate( tilde_sij_22(-NbC:NC)    ); tilde_sij_22       = 0.
  allocate( tilde_sij_33(-NbC:NC)    ); tilde_sij_33       = 0.
  allocate( tilde_sij_12(-NbC:NC)    ); tilde_sij_12       = 0.
  allocate( tilde_sij_13(-NbC:NC)    ); tilde_sij_13       = 0.
  allocate( tilde_sij_23(-NbC:NC)    ); tilde_sij_23       = 0.
  allocate( tilde_sij_kk(-NbC:NC)    ); tilde_sij_kk       = 0.
  allocate( pi_g(-NbC:NC)            ); pi_g               = 0.
  allocate( pi_t(1:NC)               ); pi_t               = 0.
  allocate( hat_rui_1(1:NC)          ); hat_rui_1          = 0.
  allocate( hat_rui_2(1:NC)          ); hat_rui_2          = 0.
  allocate( hat_rui_3(1:NC)          ); hat_rui_3          = 0.
  allocate( hat_r_tilde_sij_11(1:NC) ); hat_r_tilde_sij_11 = 0.
  allocate( hat_r_tilde_sij_22(1:NC) ); hat_r_tilde_sij_22 = 0.
  allocate( hat_r_tilde_sij_33(1:NC) ); hat_r_tilde_sij_33 = 0.
  allocate( hat_r_tilde_sij_12(1:NC) ); hat_r_tilde_sij_12 = 0.
  allocate( hat_r_tilde_sij_13(1:NC) ); hat_r_tilde_sij_13 = 0.
  allocate( hat_r_tilde_sij_23(1:NC) ); hat_r_tilde_sij_23 = 0.
  allocate( hat_r_tilde_sij_kk(1:NC) ); hat_r_tilde_sij_kk = 0.
  allocate( hat_ruiuj_11(1:NC)       ); hat_ruiuj_11       = 0.
  allocate( hat_ruiuj_22(1:NC)       ); hat_ruiuj_22       = 0.
  allocate( hat_ruiuj_33(1:NC)       ); hat_ruiuj_33       = 0.
  allocate( hat_ruiuj_12(1:NC)       ); hat_ruiuj_12       = 0.
  allocate( hat_ruiuj_13(1:NC)       ); hat_ruiuj_13       = 0.
  allocate( hat_ruiuj_23(1:NC)       ); hat_ruiuj_23       = 0.
  allocate( hat_p_term(1:NC)         ); hat_p_term         = 0.
  allocate( hat_mij_term_11(1:NC)    ); hat_mij_term_11    = 0.
  allocate( hat_mij_term_22(1:NC)    ); hat_mij_term_22    = 0.
  allocate( hat_mij_term_33(1:NC)    ); hat_mij_term_33    = 0.
  allocate( hat_mij_term_12(1:NC)    ); hat_mij_term_12    = 0.
  allocate( hat_mij_term_13(1:NC)    ); hat_mij_term_13    = 0.
  allocate( hat_mij_term_23(1:NC)    ); hat_mij_term_23    = 0.
  allocate( lij_11(1:NC)             ); lij_11             = 0.
  allocate( lij_22(1:NC)             ); lij_22             = 0.
  allocate( lij_33(1:NC)             ); lij_33             = 0.
  allocate( lij_12(1:NC)             ); lij_12             = 0.
  allocate( lij_13(1:NC)             ); lij_13             = 0.
  allocate( lij_23(1:NC)             ); lij_23             = 0.
  allocate( lij_kk(1:NC)             ); lij_kk             = 0.
  allocate( mij_11(1:NC)             ); mij_11             = 0.
  allocate( mij_22(1:NC)             ); mij_22             = 0.
  allocate( mij_33(1:NC)             ); mij_33             = 0.
  allocate( mij_12(1:NC)             ); mij_12             = 0.
  allocate( mij_13(1:NC)             ); mij_13             = 0.
  allocate( mij_23(1:NC)             ); mij_23             = 0.
  allocate( mij_kk(1:NC)             ); mij_kk             = 0.
  allocate( mijmij(1:NC)             ); mijmij             = 0.
  allocate( mijld_ij(1:NC)           ); mijld_ij           = 0.
  allocate( c_i(1:NC)                ); c_i                = 0.
  allocate( c_r(1:NC)                ); c_r                = 0.

  !-----------------------------------------!
  !   Favre-averaged velocity derevatives   !
  !-----------------------------------------!
  ! \bar{\Delta} (characteristic filter width)
  bar_delta(-NbC:NC) = volume(-NbC:NC)**(1./3.)

  ! \frac{\partial \tilde{u_i}}{\partial x}
  call GraPhi(u % n, 1, tilde_u_x,.TRUE.)
  call GraPhi(u % n, 2, tilde_u_y,.TRUE.)
  call GraPhi(u % n, 3, tilde_u_z,.TRUE.)
  call GraPhi(v % n, 1, tilde_v_x,.TRUE.)
  call GraPhi(v % n, 2, tilde_v_y,.TRUE.)
  call GraPhi(v % n, 3, tilde_v_z,.TRUE.)
  call GraPhi(w % n, 1, tilde_w_x,.TRUE.)
  call GraPhi(w % n, 2, tilde_w_y,.TRUE.)
  call GraPhi(w % n, 3, tilde_w_z,.TRUE.)

  ! \tilde{S}_{ij} =                                           &
  ! \frac{1}{2} * (\frac{\partial \tilde{u}_i}{\partial x_j} + &
  !                \frac{\partial \tilde{u}_j}{\partial x_i})
  tilde_sij_11(-NbC:NC) =       tilde_u_x(-NbC:NC)
  tilde_sij_22(-NbC:NC) =       tilde_v_y(-NbC:NC)
  tilde_sij_33(-NbC:NC) =       tilde_w_z(-NbC:NC)
  tilde_sij_12(-NbC:NC) = 0.5*( tilde_u_y(-NbC:NC) + tilde_v_x(-NbC:NC) )
  tilde_sij_13(-NbC:NC) = 0.5*( tilde_u_z(-NbC:NC) + tilde_w_x(-NbC:NC) )
  tilde_sij_23(-NbC:NC) = 0.5*( tilde_v_z(-NbC:NC) + tilde_w_y(-NbC:NC) )

  tilde_sij_kk = tilde_sij_11 + tilde_sij_22 + tilde_sij_33

  ! \tilde{S_{ij}} \tilde{S_{ij}} convolution
  pi_g =        tilde_sij_11**2 + tilde_sij_22**2 + tilde_sij_33**2 + &
         2. * ( tilde_sij_12**2 + tilde_sij_13**2 + tilde_sij_23**2 )

  !-----------------------------!
  !   test filter ^ variables   !
  !-----------------------------!

  ! volume of test filter cell
  hat_volume = volume(1:NC)

  ! \hat{\bar{\rho}}
  hat_rho = volume(1:NC) * rho % n(1:NC)

  ! \hat{\bar{\rho}\tilde{u_i}} or \hat{\bar{\rho u_i}}
  hat_rui_1 = volume(1:NC) * rho % n(1:NC) * U % n(1:NC)
  hat_rui_2 = volume(1:NC) * rho % n(1:NC) * V % n(1:NC)
  hat_rui_3 = volume(1:NC) * rho % n(1:NC) * W % n(1:NC)

  ! \hat{\bar{\rho}\tilde{S}_{ij}}
  hat_r_tilde_sij_11 = volume(1:NC) * rho % n(1:NC) * tilde_sij_11(1:NC)
  hat_r_tilde_sij_22 = volume(1:NC) * rho % n(1:NC) * tilde_sij_22(1:NC)
  hat_r_tilde_sij_33 = volume(1:NC) * rho % n(1:NC) * tilde_sij_33(1:NC)
  hat_r_tilde_sij_12 = volume(1:NC) * rho % n(1:NC) * tilde_sij_12(1:NC)
  hat_r_tilde_sij_13 = volume(1:NC) * rho % n(1:NC) * tilde_sij_13(1:NC)
  hat_r_tilde_sij_23 = volume(1:NC) * rho % n(1:NC) * tilde_sij_23(1:NC)

  ! \hat{\bar{\rho}\tilde{u_i}\tilde{u_j}} (part of Germano identity)
  hat_ruiuj_11 = volume(1:NC) * rho % n(1:NC) * U % n(1:NC) * U % n(1:NC)
  hat_ruiuj_22 = volume(1:NC) * rho % n(1:NC) * V % n(1:NC) * V % n(1:NC)
  hat_ruiuj_33 = volume(1:NC) * rho % n(1:NC) * W % n(1:NC) * W % n(1:NC)
  hat_ruiuj_12 = volume(1:NC) * rho % n(1:NC) * U % n(1:NC) * V % n(1:NC)
  hat_ruiuj_13 = volume(1:NC) * rho % n(1:NC) * U % n(1:NC) * W % n(1:NC)
  hat_ruiuj_23 = volume(1:NC) * rho % n(1:NC) * V % n(1:NC) * W % n(1:NC)

  ! \hat{2\bar{\Delta}^2\bar{rho}pi_g} (term in c_i equation)
  hat_p_term = volume(1:NC) * 2. * bar_delta(1:NC)**2 * rho % n(1:NC) * pi_g(1:NC)

  ! \hat{2\bar{\Delta}^2\bar{rho}pi_g^0.5 &
  ! (\tilde{S}_ij - \frac{1}{3}\tilde{S}_kk\delta_{ij})}
  hat_mij_term_11 = volume(1:NC) * 2. * bar_delta(1:NC)**2 * &
                               rho % n(1:NC) * pi_g(1:NC)**0.5 * &
                               ( tilde_sij_11(1:NC) - tilde_sij_kk(1:NC) / 3.)
  hat_mij_term_22 = volume(1:NC) * 2. * bar_delta(1:NC)**2 * &
                               rho % n(1:NC) * pi_g(1:NC)**0.5 * &
                               ( tilde_sij_22(1:NC) - tilde_sij_kk(1:NC) / 3.)
  hat_mij_term_33 = volume(1:NC) * 2. * bar_delta(1:NC)**2 * &
                               rho % n(1:NC) * pi_g(1:NC)**0.5 * &
                               ( tilde_sij_33(1:NC) - tilde_sij_kk(1:NC) / 3.)
  hat_mij_term_12 = volume(1:NC) * 2. * bar_delta(1:NC)**2 * &
                               rho % n(1:NC) * pi_g(1:NC)**0.5 * &
                               ( tilde_sij_12(1:NC) )
  hat_mij_term_13 = volume(1:NC) * 2. * bar_delta(1:NC)**2 * &
                               rho % n(1:NC) * pi_g(1:NC)**0.5 * &
                               ( tilde_sij_13(1:NC) )
  hat_mij_term_23 = volume(1:NC) * 2. * bar_delta(1:NC)**2 * &
                               rho % n(1:NC) * pi_g(1:NC)**0.5 * &
                               ( tilde_sij_23(1:NC) )

  do s = 1, NS ! one day this loop will be from 1 to NSinside
      c1 = SideC(1, s)
      c2 = SideC(2, s)

      if (c2 > 0 .or. typebc(c2) == BUFFER) then
        ! volume of test filter cell
        hat_volume(c1) = hat_volume(c1) + volume(c2)

        ! \hat{\bar{\rho}}
        hat_rho(c1) = hat_rho(c1) + volume(c2) * rho % n(c2)

        ! \hat{\bar{\rho}\tilde{u_i}} or \hat{\bar{\rho u_i}}
        hat_rui_1(c1) = hat_rui_1(c1) + &
                        volume(c2) * rho % n(c2) * U % n(c2)
        hat_rui_2(c1) = hat_rui_2(c1) + &
                        volume(c2) * rho % n(c2) * V % n(c2)
        hat_rui_3(c1) = hat_rui_3(c1) + &
                        volume(c2) * rho % n(c2) * W % n(c2)

        ! \hat{\bar{\rho}\tilde{S}_{ij}}
        hat_r_tilde_sij_11(c1) = hat_r_tilde_sij_11(c1) + &
                                 volume(c2) * rho % n(c2) * tilde_sij_11(c2)
        hat_r_tilde_sij_22(c1) = hat_r_tilde_sij_22(c1) + &
                                 volume(c2) * rho % n(c2) * tilde_sij_22(c2)
        hat_r_tilde_sij_33(c1) = hat_r_tilde_sij_33(c1) + &
                                 volume(c2) * rho % n(c2) * tilde_sij_33(c2)
        hat_r_tilde_sij_12(c1) = hat_r_tilde_sij_12(c1) + &
                                 volume(c2) * rho % n(c2) * tilde_sij_12(c2)
        hat_r_tilde_sij_13(c1) = hat_r_tilde_sij_13(c1) + &
                                 volume(c2) * rho % n(c2) * tilde_sij_13(c2)
        hat_r_tilde_sij_23(c1) = hat_r_tilde_sij_23(c1) + &
                                 volume(c2) * rho % n(c2) * tilde_sij_23(c2)

        ! \hat{\bar{\rho}\tilde{u_i}\tilde{u_j}} (part of Germano identity)
        hat_ruiuj_11(c1) = hat_ruiuj_11(c1) + &
                           volume(c2) * rho % n(c2) * U % n(c2) * U % n(c2)
        hat_ruiuj_22(c1) = hat_ruiuj_22(c1) + &
                           volume(c2) * rho % n(c2) * V % n(c2) * V % n(c2)
        hat_ruiuj_33(c1) = hat_ruiuj_33(c1) + &
                           volume(c2) * rho % n(c2) * W % n(c2) * W % n(c2)
        hat_ruiuj_12(c1) = hat_ruiuj_12(c1) + &
                           volume(c2) * rho % n(c2) * U % n(c2) * V % n(c2)
        hat_ruiuj_13(c1) = hat_ruiuj_13(c1) + &
                           volume(c2) * rho % n(c2) * U % n(c2) * W % n(c2)
        hat_ruiuj_23(c1) = hat_ruiuj_23(c1) + &
                           volume(c2) * rho % n(c2) * V % n(c2) * W % n(c2)

        ! \hat{2\bar{\Delta}^2\bar{rho}pi_g} (term in c_i equation)
        hat_p_term(c1) = hat_p_term(c1) + &
          volume(c2) * 2. * bar_delta(c2)**2 * rho % n(c2) * pi_g(c2)

        ! \hat{2\bar{\Delta}^2\bar{rho}pi_g^0.5 &
        ! (\tilde{S}_ij - \frac{1}{3}\tilde{S}_kk\delta_{ij})}
        hat_mij_term_11(c1) = hat_mij_term_11(c1) + &
                              volume(c2) * 2. * bar_delta(c2)**2 * &
                              rho % n(c2) * pi_g(c2)**0.5 * &
                              ( tilde_sij_11(c2) - tilde_sij_kk(c2) / 3.)
        hat_mij_term_22(c1) = hat_mij_term_22(c1) + &
                              volume(c2) * 2. * bar_delta(c2)**2 * &
                              rho % n(c2) * pi_g(c2)**0.5 * &
                              ( tilde_sij_22(c2) - tilde_sij_kk(c2) / 3.)
        hat_mij_term_33(c1) = hat_mij_term_33(c1) + &
                              volume(c2) * 2. * bar_delta(c2)**2 * &
                              rho % n(c2) * pi_g(c2)**0.5 * &
                              ( tilde_sij_33(c2) - tilde_sij_kk(c2) / 3.)
        hat_mij_term_12(c1) = hat_mij_term_12(c1) + &
                              volume(c2) * 2. * bar_delta(c2)**2 * &
                              rho % n(c2) * pi_g(c2)**0.5 * &
                              ( tilde_sij_12(c2) )
        hat_mij_term_13(c1) = hat_mij_term_13(c1) + &
                              volume(c2) * 2. * bar_delta(c2)**2 * &
                              rho % n(c2) * pi_g(c2)**0.5 * &
                              ( tilde_sij_13(c2) )
        hat_mij_term_23(c1) = hat_mij_term_23(c1) + &
                              volume(c2) * 2. * bar_delta(c2)**2 * &
                              rho % n(c2) * pi_g(c2)**0.5 * &
                              ( tilde_sij_23(c2) )
    end if
  end do

  do s = 1, NS ! one day this loop will be from 1 to NSboundary
    c1 = SideC(1, s)
    c2 = SideC(2, s) ! < 0

    if (c2 > 0) then
      ! volume of test filter cell
      hat_volume(c2) = hat_volume(c2) + volume(c1)

      ! \hat{\bar{\rho}}
      hat_rho(c2) = hat_rho(c2) + volume(c1) * rho % n(c1)

      ! \hat{\bar{\rho}\tilde{u_i}} or \hat{\bar{\rho u_i}}
      hat_rui_1(c2) = hat_rui_1(c2) + &
                      volume(c1) * rho % n(c1) * U % n(c1)

      hat_rui_2(c2) = hat_rui_2(c2) + &
                      volume(c1) * rho % n(c1) * V % n(c1)

      hat_rui_3(c2) = hat_rui_3(c2) + &
                      volume(c1) * rho % n(c1) * W % n(c1)

      ! \hat{\bar{\rho}\tilde{S}_{ij}}

      hat_r_tilde_sij_11(c2) = hat_r_tilde_sij_11(c2) + &
                               volume(c1) * rho % n(c1) * tilde_sij_11(c1)

      hat_r_tilde_sij_22(c2) = hat_r_tilde_sij_22(c2) + &
                               volume(c1) * rho % n(c1) * tilde_sij_22(c1)

      hat_r_tilde_sij_33(c2) = hat_r_tilde_sij_33(c2) + &
                               volume(c1) * rho % n(c1) * tilde_sij_33(c1)

      hat_r_tilde_sij_12(c2) = hat_r_tilde_sij_12(c2) + &
                               volume(c1) * rho % n(c1) * tilde_sij_12(c1)

      hat_r_tilde_sij_13(c2) = hat_r_tilde_sij_13(c2) + &
                               volume(c1) * rho % n(c1) * tilde_sij_13(c1)

      hat_r_tilde_sij_23(c2) = hat_r_tilde_sij_23(c2) + &
                               volume(c1) * rho % n(c1) * tilde_sij_23(c1)

      ! \hat{\bar{\rho}\tilde{u_i}\tilde{u_j}} (part of Germano identity)
      hat_ruiuj_11(c2) = hat_ruiuj_11(c2) + &
                         volume(c1) * rho % n(c1) * U % n(c1) * U % n(c1)

      hat_ruiuj_22(c2) = hat_ruiuj_22(c2) + &
                         volume(c1) * rho % n(c1) * V % n(c1) * V % n(c1)

      hat_ruiuj_33(c2) = hat_ruiuj_33(c2) + &
                         volume(c1) * rho % n(c1) * W % n(c1) * W % n(c1)

      hat_ruiuj_12(c2) = hat_ruiuj_12(c2) + &
                         volume(c1) * rho % n(c1) * U % n(c1) * V % n(c1)

      hat_ruiuj_13(c2) = hat_ruiuj_13(c2) + &
                         volume(c1) * rho % n(c1) * U % n(c1) * W % n(c1)

      hat_ruiuj_23(c2) = hat_ruiuj_23(c2) + &
                         volume(c1) * rho % n(c1) * V % n(c1) * W % n(c1)

      ! \hat{2\bar{\Delta}^2\bar{rho}pi_g} (term in c_i equation)
      hat_p_term(c2) = hat_p_term(c2) + &
        volume(c1) * 2. * bar_delta(c1)**2 * rho % n(c1) * pi_g(c1)


      ! \hat{2\bar{\Delta}^2\bar{rho}pi_g^0.5 &
      ! (\tilde{S}_ij - \frac{1}{3}\tilde{S}_kk\delta_{ij})}
      hat_mij_term_11(c2) = hat_mij_term_11(c2) + &
                            volume(c1) * 2. * bar_delta(c1)**2 * &
                            rho % n(c1) * pi_g(c1)**0.5 * &
                            ( tilde_sij_11(c1) - tilde_sij_kk(c1) / 3.)

      hat_mij_term_22(c2) = hat_mij_term_22(c2) + &
                            volume(c1) * 2. * bar_delta(c1)**2 * &
                            rho % n(c1) * pi_g(c1)**0.5 * &
                            ( tilde_sij_22(c1) - tilde_sij_kk(c1) / 3.)

      hat_mij_term_33(c2) = hat_mij_term_33(c2) + &
                            volume(c1) * 2. * bar_delta(c1)**2 * &
                            rho % n(c1) * pi_g(c1)**0.5 * &
                            ( tilde_sij_33(c1) - tilde_sij_kk(c1) / 3.)

      hat_mij_term_12(c2) = hat_mij_term_12(c2) + &
                            volume(c1) * 2. * bar_delta(c1)**2 * &
                            rho % n(c1) * pi_g(c1)**0.5 * &
                            ( tilde_sij_12(c1) )

      hat_mij_term_13(c2) = hat_mij_term_13(c2) + &
                            volume(c1) * 2. * bar_delta(c1)**2 * &
                            rho % n(c1) * pi_g(c1)**0.5 * &
                            ( tilde_sij_13(c1) )

      hat_mij_term_23(c2) = hat_mij_term_23(c2) + &
                            volume(c1) * 2. * bar_delta(c1)**2 * &
                            rho % n(c1) * pi_g(c1)**0.5 * &
                            ( tilde_sij_23(c1) )
    end if
  end do

  !------------------------------------------!
  !   Averaging through test filter volume   !
  !------------------------------------------!
  hat_rho = hat_rho / hat_volume

  hat_rui_1 = hat_rui_1 / hat_volume
  hat_rui_2 = hat_rui_2 / hat_volume
  hat_rui_3 = hat_rui_3 / hat_volume

  hat_r_tilde_sij_11 = hat_r_tilde_sij_11 / hat_volume
  hat_r_tilde_sij_22 = hat_r_tilde_sij_22 / hat_volume
  hat_r_tilde_sij_33 = hat_r_tilde_sij_33 / hat_volume
  hat_r_tilde_sij_12 = hat_r_tilde_sij_12 / hat_volume
  hat_r_tilde_sij_13 = hat_r_tilde_sij_13 / hat_volume
  hat_r_tilde_sij_23 = hat_r_tilde_sij_23 / hat_volume

  hat_r_tilde_sij_kk = hat_r_tilde_sij_11 + hat_r_tilde_sij_22 + &
                       hat_r_tilde_sij_33

  hat_ruiuj_11 = hat_ruiuj_11 / hat_volume
  hat_ruiuj_22 = hat_ruiuj_22 / hat_volume
  hat_ruiuj_33 = hat_ruiuj_33 / hat_volume
  hat_ruiuj_12 = hat_ruiuj_12 / hat_volume
  hat_ruiuj_13 = hat_ruiuj_13 / hat_volume
  hat_ruiuj_23 = hat_ruiuj_23 / hat_volume

  hat_p_term = hat_p_term / hat_volume

  hat_mij_term_11 = hat_mij_term_11 / hat_volume
  hat_mij_term_22 = hat_mij_term_22 / hat_volume
  hat_mij_term_33 = hat_mij_term_33 / hat_volume
  hat_mij_term_12 = hat_mij_term_12 / hat_volume
  hat_mij_term_13 = hat_mij_term_13 / hat_volume
  hat_mij_term_23 = hat_mij_term_23 / hat_volume

  !--------------------------------------------!
  !   Constructing relations for c_i and c_r   !
  !--------------------------------------------!

  ! \hat{\bar{\rho}\tilde{S}_{ij}} \hat{\bar{\rho}\tilde{S}_{ij}} / &
  ! \hat{\bar{\rho}}^2 convolution
  pi_t = (        hat_r_tilde_sij_11**2 + hat_r_tilde_sij_22**2 + &
                  hat_r_tilde_sij_33**2 + &
           2. * ( hat_r_tilde_sij_12**2 + hat_r_tilde_sij_13**2 + &
                  hat_r_tilde_sij_23**2 ) ) / hat_rho**2

  ! L_{ij} (Germano identity)
  lij_11 = hat_ruiuj_11 - hat_rui_1 * hat_rui_1 / hat_rho
  lij_22 = hat_ruiuj_22 - hat_rui_2 * hat_rui_2 / hat_rho
  lij_33 = hat_ruiuj_33 - hat_rui_3 * hat_rui_3 / hat_rho
  lij_12 = hat_ruiuj_12 - hat_rui_1 * hat_rui_2 / hat_rho
  lij_13 = hat_ruiuj_13 - hat_rui_1 * hat_rui_3 / hat_rho
  lij_23 = hat_ruiuj_23 - hat_rui_2 * hat_rui_3 / hat_rho

  lij_kk = lij_11 + lij_22 + lij_33

  ! \hat{\bar{\Delta}}
  hat_delta  = hat_volume**(1./3.)

  ! M_{ij}
  mij_11 = hat_mij_term_11 - 2. * hat_delta**2 * pi_t**0.5 * &
           ( hat_r_tilde_sij_11 - hat_r_tilde_sij_kk / 3. )
  mij_22 = hat_mij_term_22 - 2. * hat_delta**2 * pi_t**0.5 * &
           ( hat_r_tilde_sij_22 - hat_r_tilde_sij_kk / 3. )
  mij_33 = hat_mij_term_33 - 2. * hat_delta**2 * pi_t**0.5 * &
           ( hat_r_tilde_sij_33 - hat_r_tilde_sij_kk / 3. )
  mij_12 = hat_mij_term_12 - 2. * hat_delta**2 * pi_t**0.5 * &
           ( hat_r_tilde_sij_12 )
  mij_13 = hat_mij_term_13 - 2. * hat_delta**2 * pi_t**0.5 * &
           ( hat_r_tilde_sij_13 )
  mij_23 = hat_mij_term_23 - 2. * hat_delta**2 * pi_t**0.5 * &
           ( hat_r_tilde_sij_23 )

  mij_kk = mij_11 + mij_22 + mij_33

  ! M_{ij} M_{ij} convolution
  mijmij =        mij_11**2 + mij_22**2 + mij_33**2 + &
           2. * ( mij_12**2 + mij_13**2 + mij_23**2 )

  ! M_{ij} (L_{ij}-\frac{1}{3}L_{kk}\delta_{ij}) convolution
  mijld_ij =      mij_11 * lij_11 + mij_22 * lij_22 + mij_33 * lij_33 + &
           2. * ( mij_12 * lij_12 + mij_13 * lij_13 + mij_23 * lij_23 ) - &
           mij_kk * lij_kk / 3.

  !------------------------------------!
  !   SGS stress tensor coefficients   !
  !------------------------------------!

  ! coefficient in SGS energy
  ! this c_i is only needed when energy equation is considered, otherwise
  ! it is part of the pressure
  c_i = lij_kk / ( hat_p_term - 2. * hat_delta**2 * hat_rho * pi_t )

  ! coefficient in SGS stresses
  c_r = mijld_ij / ( mijmij + TINY )

  !----------------------------------!
  !   Updating turbulent viscosity   !
  !----------------------------------!
  vist = c_r * bar_delta**2 * rho % n * pi_g**(0.5)

  call Exchng(vist)

  !-----------------------!
  !   deallocating vars   !
  !-----------------------!

  if ( allocated( hat_volume         ) ) deallocate( hat_volume         )
  if ( allocated( hat_delta          ) ) deallocate( hat_delta          )
  if ( allocated( bar_delta          ) ) deallocate( bar_delta          )
  if ( allocated( hat_rho            ) ) deallocate( hat_rho            )
  if ( allocated( tilde_u_x          ) ) deallocate( tilde_u_x          )
  if ( allocated( tilde_u_y          ) ) deallocate( tilde_u_y          )
  if ( allocated( tilde_u_z          ) ) deallocate( tilde_u_z          )
  if ( allocated( tilde_v_x          ) ) deallocate( tilde_v_x          )
  if ( allocated( tilde_v_y          ) ) deallocate( tilde_v_y          )
  if ( allocated( tilde_v_z          ) ) deallocate( tilde_v_z          )
  if ( allocated( tilde_w_x          ) ) deallocate( tilde_w_x          )
  if ( allocated( tilde_w_y          ) ) deallocate( tilde_w_y          )
  if ( allocated( tilde_w_z          ) ) deallocate( tilde_w_z          )
  if ( allocated( tilde_sij_11       ) ) deallocate( tilde_sij_11       )
  if ( allocated( tilde_sij_22       ) ) deallocate( tilde_sij_22       )
  if ( allocated( tilde_sij_33       ) ) deallocate( tilde_sij_33       )
  if ( allocated( tilde_sij_12       ) ) deallocate( tilde_sij_12       )
  if ( allocated( tilde_sij_13       ) ) deallocate( tilde_sij_13       )
  if ( allocated( tilde_sij_23       ) ) deallocate( tilde_sij_23       )
  if ( allocated( pi_g               ) ) deallocate( pi_g               )
  if ( allocated( pi_t               ) ) deallocate( pi_t               )
  if ( allocated( hat_rui_1          ) ) deallocate( hat_rui_1          )
  if ( allocated( hat_rui_2          ) ) deallocate( hat_rui_2          )
  if ( allocated( hat_rui_3          ) ) deallocate( hat_rui_3          )
  if ( allocated( hat_r_tilde_sij_11 ) ) deallocate( hat_r_tilde_sij_11 )
  if ( allocated( hat_r_tilde_sij_22 ) ) deallocate( hat_r_tilde_sij_22 )
  if ( allocated( hat_r_tilde_sij_33 ) ) deallocate( hat_r_tilde_sij_33 )
  if ( allocated( hat_r_tilde_sij_12 ) ) deallocate( hat_r_tilde_sij_12 )
  if ( allocated( hat_r_tilde_sij_13 ) ) deallocate( hat_r_tilde_sij_13 )
  if ( allocated( hat_r_tilde_sij_23 ) ) deallocate( hat_r_tilde_sij_23 )
  if ( allocated( hat_r_tilde_sij_kk ) ) deallocate( hat_r_tilde_sij_kk )
  if ( allocated( hat_ruiuj_11       ) ) deallocate( hat_ruiuj_11       )
  if ( allocated( hat_ruiuj_22       ) ) deallocate( hat_ruiuj_22       )
  if ( allocated( hat_ruiuj_33       ) ) deallocate( hat_ruiuj_33       )
  if ( allocated( hat_ruiuj_12       ) ) deallocate( hat_ruiuj_12       )
  if ( allocated( hat_ruiuj_13       ) ) deallocate( hat_ruiuj_13       )
  if ( allocated( hat_ruiuj_23       ) ) deallocate( hat_ruiuj_23       )
  if ( allocated( hat_p_term         ) ) deallocate( hat_p_term         )
  if ( allocated( hat_mij_term_11    ) ) deallocate( hat_mij_term_11    )
  if ( allocated( hat_mij_term_22    ) ) deallocate( hat_mij_term_22    )
  if ( allocated( hat_mij_term_33    ) ) deallocate( hat_mij_term_33    )
  if ( allocated( hat_mij_term_12    ) ) deallocate( hat_mij_term_12    )
  if ( allocated( hat_mij_term_13    ) ) deallocate( hat_mij_term_13    )
  if ( allocated( hat_mij_term_23    ) ) deallocate( hat_mij_term_23    )
  if ( allocated( lij_11             ) ) deallocate( lij_11             )
  if ( allocated( lij_22             ) ) deallocate( lij_22             )
  if ( allocated( lij_33             ) ) deallocate( lij_33             )
  if ( allocated( lij_12             ) ) deallocate( lij_12             )
  if ( allocated( lij_13             ) ) deallocate( lij_13             )
  if ( allocated( lij_23             ) ) deallocate( lij_23             )
  if ( allocated( lij_kk             ) ) deallocate( lij_kk             )
  if ( allocated( mij_11             ) ) deallocate( mij_11             )
  if ( allocated( mij_22             ) ) deallocate( mij_22             )
  if ( allocated( mij_33             ) ) deallocate( mij_33             )
  if ( allocated( mij_12             ) ) deallocate( mij_12             )
  if ( allocated( mij_13             ) ) deallocate( mij_13             )
  if ( allocated( mij_23             ) ) deallocate( mij_23             )
  if ( allocated( mij_kk             ) ) deallocate( mij_kk             )
  if ( allocated( mijmij             ) ) deallocate( mijmij             )
  if ( allocated( mijld_ij           ) ) deallocate( mijld_ij           )
  if ( allocated( c_i                ) ) deallocate( c_i                )
  if ( allocated( c_r                ) ) deallocate( c_r                )

  return

end subroutine Calc_Sgs_Coefficients_Dynamic_Smagorinsky

!------------------------------------------------------------------------------!
!-------export along the line
!  filename = "hat-ppp.dat"
!  write(filename(5:7),'(I3.3)') this
!
!  open(991, FILE=adjustl(trim(filename)), status='replace', access="sequential", position='append', action='write')
!  do c=1,NC
!    if ( c .ne. 0 .and. zc(c)>=-1.0 .and. zc(c)<=1.0 .and.  yc(c)<=0.03125+0.0001 .and. yc(c)>=0.03125-0.0001 ) then !problem 3
!        write(991,'(66ES26.16E3)') xc(c),& ! 1
!                                   yc(c),& ! 2
!                                   zc(c),& ! 3
!                                   U % n(c),& ! 4
!                                   Moin_Problem_u (xc(c), yc(c), 0.125),& ! 5
!                                   V % n(c),& ! 6
!                                   Moin_Problem_v (xc(c), yc(c), 0.125),& ! 7
!                                   W % n(c),& ! 8
!                                   Zmix % n(c),& ! 9
!                                   Moin_Problem_z (xc(c), yc(c), 0.125),& ! 10
!                                   Rho % n(c),& ! 11
!                                   Moin_Problem_rho (xc(c), yc(c), 0.125),& ! 12
!                                   hat_volume(c),& ! 13
!                                   hat_rho(c),& ! 14
!                                   hat_rui_1(c),& ! 15
!                                   hat_rui_2(c),& ! 16
!                                   hat_rui_3(c),& ! 17
!                                   hat_r_tilde_sij_11(c),& ! 18
!                                   hat_r_tilde_sij_22(c),& ! 19
!                                   hat_r_tilde_sij_33(c),& ! 20
!                                   hat_r_tilde_sij_12(c),& ! 21
!                                   hat_r_tilde_sij_13(c),& ! 22
!                                   hat_r_tilde_sij_23(c),& ! 23
!                                   hat_ruiuj_11(c),& ! 24
!                                   hat_ruiuj_22(c),& ! 25
!                                   hat_ruiuj_33(c),& ! 26
!                                   hat_ruiuj_12(c),& ! 27
!                                   hat_ruiuj_13(c),& ! 28
!                                   hat_ruiuj_23(c),& ! 29
!                                   hat_p_term(c),& ! 30
!                                   hat_mij_term_11(c),& ! 31
!                                   hat_mij_term_22(c),& ! 32
!                                   hat_mij_term_33(c),& ! 33
!                                   hat_mij_term_12(c),& ! 34
!                                   hat_mij_term_13(c),& ! 35
!                                   hat_mij_term_23(c),& ! 36
!                                   pi_t(c),& ! 37
!                                   lij_11(c),& ! 38
!                                   lij_22(c),& ! 39
!                                   lij_33(c),& ! 40
!                                   lij_12(c),& ! 41
!                                   lij_13(c),& ! 42
!                                   lij_23(c),& ! 43
!                                   hat_delta(c),& ! 44
!                                   mij_11(c),& ! 45
!                                   mij_22(c),& ! 46
!                                   mij_33(c),& ! 47
!                                   mij_12(c),& ! 48
!                                   mij_13(c),& ! 49
!                                   mij_23(c),& ! 50
!                                   mijmij(c),& ! 51
!                                   mijld_ij(c),& ! 52
!                                   c_i(c),& ! 53
!                                   c_r(c),& ! 54
!                                   tilde_u_x(c),& ! 55
!                                   tilde_u_y(c),& ! 56
!                                   tilde_u_z(c),& ! 57
!                                   tilde_v_x(c),& ! 58
!                                   tilde_v_y(c),& ! 59
!                                   tilde_v_z(c),& ! 60
!                                   tilde_w_x(c),& ! 61
!                                   tilde_w_y(c),& ! 62
!                                   tilde_w_z(c),& ! 63
!                                   volume(c),& ! 64
!                                   pi_g(c),& ! 65
!                                   bar_delta(c) ! 66
!
!  end if
!end do
!close(991)
!
!call wait
!call exit(1)
!------------------------------------------------------------------------------!