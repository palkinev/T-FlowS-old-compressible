!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
!                                 !                                    !
!                                 !   Bojan Niceno                     !
!      Parameter definitions      !   Delft University of Technology   !
!         for all modules         !   Section Heat Transfer            !
!                                 !   niceno@duttwta.wt.tn.tudelft.nl  !
!                                 !                                    !
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
module allp_mod

  implicit none

  integer, parameter, public :: &
    MAXP =   200, MAXL = 1000,  &
    MAXPRO   = 1024,            & ! max. n. of processors    
    INFLOW   =    1,            & ! boundary condition       
    WALL     =    2,            & ! boundary condition
    OUTFLOW  =    3,            & ! boundary condition
    SYMMETRY =    4,            & ! boundary condition
    CONVECT  =    5,            & ! boundary condition
    PRESSURE =   12,            & ! boundary condition
    PERIODIC =   13,            & ! boundary condition
    BUFFER   =   11,            & ! boundary condition
    WALLFL   =    6               ! boundary condition

  integer, parameter, public :: &
    FLUID    =    7,            & ! material state: fluid
    SOLID    =    8               ! material state: solid

  real, parameter, public ::  &
    HUGE=1.e+30, TINY=1.e-64

  !----- Unknown type (rename to FIELD)
  type, public :: Unknown !this
    real, pointer :: n(:)         ! new value for this
    real, pointer :: o(:),  oo(:) ! old and older then old for this
    real, pointer :: C(:)         ! convective fluxes for this
    real, pointer :: X(:)         ! surfce sources for this
    real, pointer :: mean(:)      ! long time average for this
    real, pointer :: filt(:)      ! LES spatial filter
    real, pointer :: q(:)         ! flux for a variable for this
    real, pointer :: fluc(:)
    real, pointer :: source(:)    ! source for this
    real          :: URF          ! under relaxation factor for this
    real          :: Stol         ! solver tolerance for this
    real          :: bound(128)   ! boundary values
    real          :: init(128)    ! initial values
    real          :: pro(11024)   ! inlfow profile
    real          :: Sigma        ! sigma
  end type Unknown
  !————————————————————————————————————————————————————————————————————————————————————————
end module

! precision parameters (you can just use compiler precision flags)
!integer, parameter ::                              &
!    sp = kind(1.0),                                &
!    dp = selected_real_kind(2*precision(1.0_sp)),  &
!    qp = selected_real_kind(2*precision(1.0_dp))