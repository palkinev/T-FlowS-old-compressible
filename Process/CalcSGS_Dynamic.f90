!==============================================================================!
  subroutine Calc_Sgs_Coefficients_Dynamic_Smagorinsky
!------------------------------------------------------------------------------!
!   Calculates coefficients in Smagorinsky model                               !
!   with dynamic procedure for compressible case.                              !
!   This function follows the manual in latex called "???"                     !
!------------------------------------------------------------------------------!
!   Attention: LaTex formalism is used in this function                        !
!------------------------------------------------------------------------------!
!   All variables considered in LES code are spatially filtered by cells       !
!   themselves and denoted as vars with bar above them:                        !
!   and denoted as vars with bar:                                              !
!   \bar{rho}, \bar{U}, \bar{V},_\bar{W}, \bar{Zmix}.                          !
!   This procedure uses spectral information on the flow, thus larger test     !
!   filter is introduced and denoted as \hat{} (^), it is applied after cell   !
!   filter \bar{}.                                                             !
!----------------------------------[Modules]-----------------------------------!
  use allp_mod, only : tiny
  use all_mod
  use pro_mod
  use les_mod
  use rans_mod
!------------------------------------------------------------------------------!
  implicit none
!-----------------------------------[Locals]-----------------------------------!
  integer           :: c, j, cj
  real, allocatable :: hat_volume(:), hat_delta(:), bar_delta(:)
  real, allocatable :: hat_Rho(:)
  real, allocatable :: tilde_U_x(:), tilde_U_y(:), tilde_U_z(:), &
                       tilde_V_x(:), tilde_V_y(:), tilde_V_z(:), &
                       tilde_W_x(:), tilde_W_y(:), tilde_W_z(:)
  real, allocatable :: tilde_Sij_11(:), tilde_Sij_22(:), tilde_Sij_33(:), &
                       tilde_Sij_12(:), tilde_Sij_13(:), tilde_Sij_23(:), &
                       tilde_Sij_kk(:)
  real, allocatable :: PI_G(:), PI_T(:)
  real, allocatable :: hat_rUi_1(:), hat_rUi_2(:), hat_rUi_3(:)
  real, allocatable :: hat_r_tilde_Sij_11(:), hat_r_tilde_Sij_22(:), hat_r_tilde_Sij_33(:), &
                       hat_r_tilde_Sij_12(:), hat_r_tilde_Sij_13(:), hat_r_tilde_Sij_23(:), &
                       hat_r_tilde_Sij_kk(:)
  real, allocatable :: hat_rUiUj_11(:), hat_rUiUj_22(:), hat_rUiUj_33(:), &
                       hat_rUiUj_12(:), hat_rUiUj_13(:), hat_rUiUj_23(:)
  real, allocatable :: hat_P_term(:)
  real, allocatable :: hat_Mij_term_11(:), hat_Mij_term_22(:), hat_Mij_term_33(:), &
                       hat_Mij_term_12(:), hat_Mij_term_13(:), hat_Mij_term_23(:)
  real, allocatable :: Lij_11(:), Lij_22(:), Lij_33(:), &
                       Lij_12(:), Lij_13(:), Lij_23(:), &
                       Lij_kk(:)
  real, allocatable :: Mij_11(:), Mij_22(:), Mij_33(:), &
                       Mij_12(:), Mij_13(:), Mij_23(:), &
                       Mij_kk(:), MijMij(:)
  real, allocatable :: MijLd_ij(:)
  real, allocatable :: C_I(:), C_R(:)
                       
  character(len=80) :: name_out
!==============================================================================!

  ! I do not yet know why we are doing this
  call Exchng(Ur % n)
  call Exchng(Vr % n)
  call Exchng(Wr % n)
  call Exchng(U % n)
  call Exchng(V % n)
  call Exchng(W % n)

  !---------------------!
  !   allocating vars   !
  !---------------------!
  allocate( hat_volume(1:NC) ); hat_volume = 0.
  allocate( hat_delta(1:NC) ); hat_delta = 0.
  allocate( bar_delta(1:NC) ); bar_delta = 0.
  allocate( hat_Rho(1:NC) ); hat_Rho = 0.
  allocate( tilde_U_x(1:NC) ); tilde_U_x = 0.
  allocate( tilde_U_y(1:NC) ); tilde_U_y = 0.
  allocate( tilde_U_z(1:NC) ); tilde_U_z = 0.
  allocate( tilde_V_x(1:NC) ); tilde_V_x = 0.
  allocate( tilde_V_y(1:NC) ); tilde_V_y = 0.
  allocate( tilde_V_z(1:NC) ); tilde_V_z = 0.
  allocate( tilde_W_x(1:NC) ); tilde_W_x = 0.
  allocate( tilde_W_y(1:NC) ); tilde_W_y = 0.
  allocate( tilde_W_z(1:NC) ); tilde_W_z = 0.
  allocate( tilde_Sij_11(1:NC) ); tilde_Sij_11 = 0.
  allocate( tilde_Sij_22(1:NC) ); tilde_Sij_22 = 0.
  allocate( tilde_Sij_33(1:NC) ); tilde_Sij_33 = 0.
  allocate( tilde_Sij_12(1:NC) ); tilde_Sij_12 = 0.
  allocate( tilde_Sij_13(1:NC) ); tilde_Sij_13 = 0.
  allocate( tilde_Sij_23(1:NC) ); tilde_Sij_23 = 0.
  allocate( PI_G(1:NC) ); PI_G = 0.
  allocate( PI_T(1:NC) ); PI_T = 0.
  allocate( hat_rUi_1(1:NC) ); hat_rUi_1 = 0.
  allocate( hat_rUi_2(1:NC) ); hat_rUi_2 = 0.
  allocate( hat_rUi_3(1:NC) ); hat_rUi_3 = 0.
  allocate( hat_r_tilde_Sij_11(1:NC) ); hat_r_tilde_Sij_11 = 0.
  allocate( hat_r_tilde_Sij_22(1:NC) ); hat_r_tilde_Sij_22 = 0.
  allocate( hat_r_tilde_Sij_33(1:NC) ); hat_r_tilde_Sij_33 = 0.
  allocate( hat_r_tilde_Sij_12(1:NC) ); hat_r_tilde_Sij_12 = 0.
  allocate( hat_r_tilde_Sij_13(1:NC) ); hat_r_tilde_Sij_13 = 0.
  allocate( hat_r_tilde_Sij_23(1:NC) ); hat_r_tilde_Sij_23 = 0.
  allocate( hat_r_tilde_Sij_kk(1:NC) ); hat_r_tilde_Sij_kk = 0.
  allocate( hat_rUiUj_11(1:NC) ); hat_rUiUj_11 = 0.
  allocate( hat_rUiUj_22(1:NC) ); hat_rUiUj_22 = 0.
  allocate( hat_rUiUj_33(1:NC) ); hat_rUiUj_33 = 0.
  allocate( hat_rUiUj_12(1:NC) ); hat_rUiUj_12 = 0.
  allocate( hat_rUiUj_13(1:NC) ); hat_rUiUj_13 = 0.
  allocate( hat_rUiUj_23(1:NC) ); hat_rUiUj_23 = 0.
  allocate( hat_P_term(1:NC) ); hat_P_term = 0.
  allocate( hat_Mij_term_11(1:NC) ); hat_Mij_term_11 = 0.
  allocate( hat_Mij_term_22(1:NC) ); hat_Mij_term_22 = 0.
  allocate( hat_Mij_term_33(1:NC) ); hat_Mij_term_33 = 0.
  allocate( hat_Mij_term_12(1:NC) ); hat_Mij_term_12 = 0.
  allocate( hat_Mij_term_13(1:NC) ); hat_Mij_term_13 = 0.
  allocate( hat_Mij_term_23(1:NC) ); hat_Mij_term_23 = 0.
  allocate( Lij_11(1:NC) ); Lij_11 = 0.
  allocate( Lij_22(1:NC) ); Lij_22 = 0.
  allocate( Lij_33(1:NC) ); Lij_33 = 0.
  allocate( Lij_12(1:NC) ); Lij_12 = 0.
  allocate( Lij_13(1:NC) ); Lij_13 = 0.
  allocate( Lij_23(1:NC) ); Lij_23 = 0.
  allocate( Lij_kk(1:NC) ); Lij_kk = 0.
  allocate( Mij_11(1:NC) ); Mij_11 = 0.
  allocate( Mij_22(1:NC) ); Mij_22 = 0.
  allocate( Mij_33(1:NC) ); Mij_33 = 0.
  allocate( Mij_12(1:NC) ); Mij_12 = 0.
  allocate( Mij_13(1:NC) ); Mij_13 = 0.
  allocate( Mij_23(1:NC) ); Mij_23 = 0.
  allocate( Mij_kk(1:NC) ); Mij_kk = 0.
  allocate( MijMij(1:NC) ); MijMij = 0.
  allocate( MijLd_ij(1:NC) ); MijLd_ij = 0.
  allocate( C_I(1:NC) ); C_I = 0.
  allocate( C_R(1:NC) ); C_R = 0.


  ! \hat{\bar{\rho u_i}} or \hat{\bar{\rho}\tilde{u_i}}
  hat_rU = 0.0
  hat_rV = 0.0
  hat_rW = 0.0

  ! \hat{\bar{\Delta}}
  hat_delta = 0.0

  ! volume of test filter
  hat_volume = 0.0






  ! \bar{\Delta} (characteristic filter width)
  bar_delta = volume(c)**(1./3.) 

  !-----------------------------------------!
  !   Favre-averaged velocity derevatives   !
  !-----------------------------------------!
  ! \frac{\partial \tilde{u_i}}{\partial x}
  call GraPhi(U, 1, tilde_U_x,.TRUE.)
  call GraPhi(U, 2, tilde_U_y,.TRUE.)
  call GraPhi(U, 3, tilde_U_z,.TRUE.)
  call GraPhi(V, 1, tilde_V_x,.TRUE.)
  call GraPhi(V, 2, tilde_V_y,.TRUE.)
  call GraPhi(V, 3, tilde_V_z,.TRUE.)
  call GraPhi(W, 1, tilde_W_x,.TRUE.)
  call GraPhi(W, 2, tilde_W_y,.TRUE.)
  call GraPhi(W, 3, tilde_W_z,.TRUE.)

  ! \tilde{S}_{ij} =                                           &
  ! \frac{1}{2} * (\frac{\partial \tilde{u}_i}{\partial x_j} + &
  !                \frac{\partial \tilde{u}_j}{\partial x_i})
  tilde_Sij_11 =       tilde_U_x 
  tilde_Sij_22 =       tilde_V_y
  tilde_Sij_33 =       tilde_W_z
  tilde_Sij_12 = 0.5*( tilde_U_y + tilde_V_x )
  tilde_Sij_13 = 0.5*( tilde_U_z + tilde_W_x )
  tilde_Sij_23 = 0.5*( tilde_V_z + tilde_W_y )

  tilde_Sij_kk = tilde_Sij_11 + tilde_Sij_22 + tilde_Sij_33

  ! \tilde{S_{ij}} \tilde{S_{ij}} convolution
  PI_G =        tilde_Sij_11**2. + tilde_Sij_22**2. + tilde_Sij_33**2. + &
         2. * ( tilde_Sij_12**2. + tilde_Sij_13**2. + tilde_Sij_23**2. )

  !---------------------------!
  !   test filter variables   !
  !---------------------------!
  do c = 1, NC
    do j = Acol(c), Acol(c + 1) - 1 ! going through all cell neighbours (and cell itself)
      cj = Arow(j) ! neighbour cell ID

        ! \hat{\bar{dV}}
        hat_volume(c) = hat_volume(c) + volume(cj)

        ! \hat{\bar{\rho}}
        hat_Rho(c) = hat_Rho(c) + volume(cj) * Rho % n(cj)

        ! \hat{\bar{\rho}\tilde{u_i}}
        hat_rUi_1(c) = hat_rUi_1(c) + volume(cj) * Rho % n(cj) * U % n(cj)
        hat_rUi_2(c) = hat_rUi_2(c) + volume(cj) * Rho % n(cj) * V % n(cj)
        hat_rUi_3(c) = hat_rUi_3(c) + volume(cj) * Rho % n(cj) * W % n(cj)

        ! \hat{\bar{\rho}\tilde{S}_{ij}} (part of Germano identity)
        hat_r_tilde_Sij_11(c) = hat_r_tilde_Sij_11(c) + volume(cj) * Rho % n(cj) * tilde_Sij_11(cj)
        hat_r_tilde_Sij_22(c) = hat_r_tilde_Sij_22(c) + volume(cj) * Rho % n(cj) * tilde_Sij_22(cj)
        hat_r_tilde_Sij_33(c) = hat_r_tilde_Sij_33(c) + volume(cj) * Rho % n(cj) * tilde_Sij_33(cj)
        hat_r_tilde_Sij_12(c) = hat_r_tilde_Sij_12(c) + volume(cj) * Rho % n(cj) * tilde_Sij_12(cj)
        hat_r_tilde_Sij_13(c) = hat_r_tilde_Sij_13(c) + volume(cj) * Rho % n(cj) * tilde_Sij_13(cj)
        hat_r_tilde_Sij_23(c) = hat_r_tilde_Sij_23(c) + volume(cj) * Rho % n(cj) * tilde_Sij_23(cj)
        
        ! \hat{\bar{\rho}\tilde{u_i}\tilde{u_j}} (part of Germano identity)
        hat_rUiUj_11(c) = hat_rUiUj_11(c) + volume(cj) * Rho % n(cj) * U % n(cj) * U % n(cj)
        hat_rUiUj_22(c) = hat_rUiUj_22(c) + volume(cj) * Rho % n(cj) * V % n(cj) * V % n(cj)
        hat_rUiUj_33(c) = hat_rUiUj_33(c) + volume(cj) * Rho % n(cj) * W % n(cj) * W % n(cj)
        hat_rUiUj_12(c) = hat_rUiUj_12(c) + volume(cj) * Rho % n(cj) * U % n(cj) * V % n(cj)
        hat_rUiUj_13(c) = hat_rUiUj_13(c) + volume(cj) * Rho % n(cj) * U % n(cj) * W % n(cj)
        hat_rUiUj_23(c) = hat_rUiUj_23(c) + volume(cj) * Rho % n(cj) * V % n(cj) * W % n(cj)

        ! \hat{2\bar{\Delta}^2\bar{rho}PI_G} (term in C_I equation)
        hat_P_term(c) = hat_P_term(c) + volume(cj) * 2. * bar_delta(cj)**2. * Rho % n(cj) * PI_G(cj)

        ! \hat{2\bar{\Delta}^2\bar{rho}PI_G^0.5 (\tilde{S}_ij - \frac{1}{3}\tilde{S}_kk\delta_{ij})}
        hat_Mij_term_11(c) = hat_Mij_term_11(c) + volume(cj) * 2. * bar_delta(cj)**2. * Rho % n(cj) * PI_G(cj)**0.5 * ( tilde_Sij_11(cj) - tilde_Sij_kk(cj) / 3.)
        hat_Mij_term_22(c) = hat_Mij_term_22(c) + volume(cj) * 2. * bar_delta(cj)**2. * Rho % n(cj) * PI_G(cj)**0.5 * ( tilde_Sij_22(cj) - tilde_Sij_kk(cj) / 3.)
        hat_Mij_term_33(c) = hat_Mij_term_33(c) + volume(cj) * 2. * bar_delta(cj)**2. * Rho % n(cj) * PI_G(cj)**0.5 * ( tilde_Sij_33(cj) - tilde_Sij_kk(cj) / 3.)
        hat_Mij_term_12(c) = hat_Mij_term_12(c) + volume(cj) * 2. * bar_delta(cj)**2. * Rho % n(cj) * PI_G(cj)**0.5 * ( tilde_Sij_12(cj) )
        hat_Mij_term_13(c) = hat_Mij_term_13(c) + volume(cj) * 2. * bar_delta(cj)**2. * Rho % n(cj) * PI_G(cj)**0.5 * ( tilde_Sij_13(cj) )
        hat_Mij_term_23(c) = hat_Mij_term_23(c) + volume(cj) * 2. * bar_delta(cj)**2. * Rho % n(cj) * PI_G(cj)**0.5 * ( tilde_Sij_23(cj) )

    end do
  end do
  
  !------------------------------------------!
  !   averaging through test filter volume   !
  !------------------------------------------!
  hat_Rho = hat_Rho / hat_volume

  hat_rUi_1 = hat_rUi_1 / hat_volume
  hat_rUi_2 = hat_rUi_2 / hat_volume
  hat_rUi_3 = hat_rUi_3 / hat_volume

  hat_r_tilde_Sij_11 = hat_r_tilde_Sij_11 / hat_volume
  hat_r_tilde_Sij_22 = hat_r_tilde_Sij_22 / hat_volume
  hat_r_tilde_Sij_33 = hat_r_tilde_Sij_33 / hat_volume
  hat_r_tilde_Sij_12 = hat_r_tilde_Sij_12 / hat_volume
  hat_r_tilde_Sij_13 = hat_r_tilde_Sij_13 / hat_volume
  hat_r_tilde_Sij_23 = hat_r_tilde_Sij_23 / hat_volume

  hat_r_tilde_Sij_kk = hat_r_tilde_Sij_11 + hat_r_tilde_Sij_22 + hat_r_tilde_Sij_33

  hat_rUiUj_11 = hat_rUiUj_11 / hat_volume
  hat_rUiUj_22 = hat_rUiUj_22 / hat_volume
  hat_rUiUj_33 = hat_rUiUj_33 / hat_volume
  hat_rUiUj_12 = hat_rUiUj_12 / hat_volume
  hat_rUiUj_13 = hat_rUiUj_13 / hat_volume
  hat_rUiUj_23 = hat_rUiUj_23 / hat_volume

  hat_P_term = hat_P_term / hat_volume

  hat_Mij_term_11 = hat_Mij_term_11 / hat_volume
  hat_Mij_term_22 = hat_Mij_term_22 / hat_volume
  hat_Mij_term_33 = hat_Mij_term_33 / hat_volume
  hat_Mij_term_12 = hat_Mij_term_12 / hat_volume
  hat_Mij_term_13 = hat_Mij_term_13 / hat_volume
  hat_Mij_term_23 = hat_Mij_term_23 / hat_volume

  !--------------------------------------------!
  !   constructing relations for C_I and C_R   !
  !--------------------------------------------!

  ! \hat{\bar{\rho}\tilde{S}_{ij}} \hat{\bar{\rho}\tilde{S}_{ij}} / \hat{\bar{\rho}}^2 convolution
  PI_T = (        hat_r_tilde_Sij_11**2. + hat_r_tilde_Sij_22**2. + hat_r_tilde_Sij_33**2. + &
           2. * ( hat_r_tilde_Sij_12**2. + hat_r_tilde_Sij_13**2. + hat_r_tilde_Sij_23**2. ) ) / hat_Rho**2.

  ! L_{ij} (Germano identity)
  Lij_11 = hat_rUiUj_11 - hat_rUi_1 * hat_rUi_1 / hat_Rho
  Lij_22 = hat_rUiUj_22 - hat_rUi_2 * hat_rUi_2 / hat_Rho
  Lij_33 = hat_rUiUj_33 - hat_rUi_3 * hat_rUi_3 / hat_Rho
  Lij_12 = hat_rUiUj_12 - hat_rUi_1 * hat_rUi_2 / hat_Rho
  Lij_13 = hat_rUiUj_13 - hat_rUi_1 * hat_rUi_3 / hat_Rho
  Lij_23 = hat_rUiUj_23 - hat_rUi_2 * hat_rUi_3 / hat_Rho

  Lij_kk = Lij_11 + Lij_22 + Lij_33

  ! \hat{\bar{\Delta}}
  hat_delta  = hat_volume(c)**(1./3.)

  ! M_{ij}
  Mij_11 = hat_Mij_term_11 - 2. * hat_delta**2. * PI_T**0.5 * ( hat_r_tilde_Sij_11 - hat_r_tilde_Sij_kk / 3. )
  Mij_22 = hat_Mij_term_22 - 2. * hat_delta**2. * PI_T**0.5 * ( hat_r_tilde_Sij_22 - hat_r_tilde_Sij_kk / 3. )
  Mij_33 = hat_Mij_term_33 - 2. * hat_delta**2. * PI_T**0.5 * ( hat_r_tilde_Sij_33 - hat_r_tilde_Sij_kk / 3. )
  Mij_12 = hat_Mij_term_12 - 2. * hat_delta**2. * PI_T**0.5 * ( hat_r_tilde_Sij_12 )
  Mij_13 = hat_Mij_term_13 - 2. * hat_delta**2. * PI_T**0.5 * ( hat_r_tilde_Sij_13 )
  Mij_23 = hat_Mij_term_23 - 2. * hat_delta**2. * PI_T**0.5 * ( hat_r_tilde_Sij_23 )

  Mij_kk = Mij_11 + Mij_22 + Mij_33

  ! M_{ij} M_{ij} convolution
  MijMij =        Mij_11**2. + Mij_22**2. + Mij_33**2. + &
           2. * ( Mij_12**2. + Mij_13**2. + Mij_23**2. )

  ! M_{ij} (L_{ij}-\frac{1}{3}L_{kk}\delta_{ij}) convolution
  MijLd_ij =      Mij_11 * Lij_11 + Mij_22 * Lij_22 + Mij_33 * Lij_33 + &
           2. * ( Mij_12 * Lij_12 + Mij_13 * Lij_13 + Mij_23 * Lij_23 ) - &
           Mij_kk * Lij_kk / 3.


  ! coefficient in SGS energy
  C_I = Lij_kk / ( hat_P_term - 2. * hat_delta**2. * hat_Rho * PI_T )

  ! coefficient in SGS stresses
  C_R = MijLd_ij / ( MijMij + tiny )

  do c=1,NC

    if(C_R(c) < 0.0) then
      C_R(c) = 0.0
    else if(C_R(c) > 0.04) then
      C_R(c) = 0.04
    end if

  end do

  return
END subroutine Calc_Sgs_Coefficients_Dynamic_Smagorinsky
