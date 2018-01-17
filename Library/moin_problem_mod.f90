!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
!                              !                                       !
!                              !   egor palkin                         !
!   module that controls       !   IT SB RAS                           !
!   conflicting sub-modules    !                                       !
!   built by moin for various  !   palkinev89@gmail.com                !
!   laminar and turbulent      !                                       !
!   test cases to model        !                                       !
!   reacting flows             !                                       !
!                              !                                       !
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
module Moin_Problem_Mod

    implicit none

    integer, private  :: Moin_Problem_ID
    real,    private  :: rho_0, rho_1, VIS, VISZmix
    real,    private  :: time
    ! rho_0 > rho_1 (maybe reverse it?)

    ! functions
    public :: Moin_Problem
    public :: Rho0
    public :: Rho1
    public :: Moin_Set_Time
    public :: Moin_Time
    public :: Moin_print_error
    public :: Moin_Problem_Init
    public :: Moin_Problem_1_Init
    public :: Moin_Problem_2_Init
    public :: Moin_Problem_3_Init
    public :: Moin_Problem_4_Init
    public :: Moin_Problem_Table_Init
    public :: Moin_Problem_u
    public :: Moin_Problem_u_src
    public :: Moin_Problem_v
    public :: Moin_Problem_v_src
    public :: Moin_Problem_p
    public :: Moin_Problem_rho
    public :: Moin_Problem_z
    public :: Moin_Problem_z_src
    public :: Moin_Problem_visc_dyn_of_x
    public :: Moin_Problem_visc_dyn_of_Z
    public :: Moin_Problem_EOS
    public :: Moin_Problem_1_u
    public :: Moin_Problem_1_u_src
    public :: Moin_Problem_1_v
    public :: Moin_Problem_1_v_src
    public :: Moin_Problem_1_p
    public :: Moin_Problem_1_rho
    public :: Moin_Problem_1_z
    public :: Moin_Problem_1_z_src
    public :: Moin_Problem_2_u
    public :: Moin_Problem_2_u_src
    public :: Moin_Problem_2_v
    public :: Moin_Problem_2_v_src
    public :: Moin_Problem_2_p
    public :: Moin_Problem_2_rho
    public :: Moin_Problem_2_z
    public :: Moin_Problem_2_z_src
    public :: Moin_Problem_3_u
    public :: Moin_Problem_3_u_src
    public :: Moin_Problem_3_v
    public :: Moin_Problem_3_v_src
    public :: Moin_Problem_3_p
    public :: Moin_Problem_3_rho
    public :: Moin_Problem_3_z
    public :: Moin_Problem_3_z_src
    public :: Moin_Problem_4_u
    public :: Moin_Problem_4_u_src
    public :: Moin_Problem_4_v
    public :: Moin_Problem_4_v_src
    public :: Moin_Problem_4_p
    public :: Moin_Problem_4_rho
    public :: Moin_Problem_4_z
    public :: Moin_Problem_4_z_src
    public :: Moin_Problem_Table_u
    public :: Moin_Problem_Table_u_src
    public :: Moin_Problem_Table_v
    public :: Moin_Problem_Table_v_src
    public :: Moin_Problem_Table_p
    public :: Moin_Problem_Table_rho
    public :: Moin_Problem_Table_z
    public :: Moin_Problem_Table_z_src
    public :: Moin_Problem_Table_EOS
    public :: Moin_Problem_Table_visc_dyn_of_Z
    public :: Moin_Problem_Table_visc_dyn_of_x
!————————————————————————————————————————————————————————————————————————————————————————
contains
    !————————————————————————————————————————————————————————————————————————————————————————

    !————————————————————————————————————————————————————————————————————————————————————————
    !---- CONTROL SUBMODULE ( OR REDIRECT)
    !————————————————————————————————————————————————————————————————————————————————————————
    subroutine Moin_Problem_Init ( &
        Moin_Problem_ID_in, &
        rho_0_in, &
        rho_1_in, &
        vis_in, &
        visZmix_in, &
        EOS_of_Zmix_in, &
        ZmixVisc_Dyn_of_Zmix_in, &
        Zmix_SRC_of_Zmix_in, &
        vel_Init_in, &
        rho_Init_in, &
        Zmix_Init_in, &
        Visc_dyn_Init_in, &
        pressure_Init_in)

        use par_Mod
        implicit none
        integer,                 intent (in) :: Moin_Problem_ID_in
        real,        intent (in) :: rho_0_in, rho_1_in, vis_in, visZmix_in
        character (*), optional, intent (in) :: EOS_of_Zmix_in, ZmixVisc_Dyn_of_Zmix_in, Zmix_SRC_of_Zmix_in
        character (*), optional, intent (in) :: vel_Init_in, rho_Init_in, Zmix_Init_in, Visc_dyn_Init_in, pressure_Init_in

225     format (" # Moin module: ", A) ! output style

        if (this < 2) write(*,225) "this module was meant to test this code in various Moin test problems"

        Moin_Problem_ID  = Moin_Problem_ID_in
        rho_0            = rho_0_in
        rho_1            = rho_1_in
        VIS              = VIS_in
        VISZmix          = VISZmix_in


        select case (Moin_Problem_ID)
            case (1)
                if (this < 2) write(*,225) "1d problem 1 is called (isothermal mixing)"
                call Moin_Problem_1_Init
            case (2)
                if (this < 2) write(*,225) "2d problem 2 is called (flame front)"
                call Moin_Problem_2_Init
            case (3)
                if (this < 2) write(*,225) "2d problem 3 is called (periodic mixing)"
                call Moin_Problem_3_Init
            case (4)
                if (this < 2) write(*,225) "2d problem 4 is called (Railegh-Tailor)"
                call Moin_Problem_4_Init
            case default

                if (this < 2 .and. Moin_Problem_ID.eq.10) write(*,225) "1d problem 1 with table variables is called"
                if (this < 2 .and. Moin_Problem_ID.eq.30) write(*,225) "2d problem 3 with table variables is called (periodic mixing)"

                call Moin_Problem_Table_Init( &
                    EOS_of_Zmix_in, &
                    ZmixVisc_Dyn_of_Zmix_in, &
                    Zmix_SRC_of_Zmix_in, &
                    vel_Init_in, &
                    rho_Init_in, &
                    Zmix_Init_in, &
                    Visc_dyn_Init_in, &
                    pressure_Init_in &
                    )                

        end select


    end subroutine Moin_Problem_Init
    !————————————————————————————————————————————————————————————————————————————————————————
    function Moin_Problem () !just returns Moin_Problem_ID
        implicit none
        integer :: Moin_Problem

        Moin_Problem = Moin_Problem_ID

    end function Moin_Problem
    !————————————————————————————————————————————————————————————————————————————————————————
    function Rho0 () !just returns rho_0
        implicit none
        real :: Rho0

        Rho0 = rho_0

    end function Rho0
    !————————————————————————————————————————————————————————————————————————————————————————
    function Rho1 () !just returns rho_0
        implicit none
        real :: Rho1

        Rho1 = rho_1

    end function Rho1
    !————————————————————————————————————————————————————————————————————————————————————————
    subroutine Moin_Set_Time (time_in) ! sets time for moin functions
        implicit none
        real :: time_in

        time = time_in

    end subroutine Moin_Set_Time
    !————————————————————————————————————————————————————————————————————————————————————————
    function Moin_Time () !just returns current time
        implicit none
        real :: Moin_Time

        Moin_Time = time

    end function Moin_Time
    !————————————————————————————————————————————————————————————————————————————————————————
    subroutine Moin_print_error (time_cur) ! print L2 norm errors
        use allp_Mod, only : huge
        use all_Mod
        use pro_Mod
        use par_Mod

        implicit none
        real, intent (in) :: time_cur
        real :: L2_U, L2_V, L2_P, L2_Z, L2_RHO, vol, p_sum, pn_min, pex_min
        integer :: c

        L2_U   = 0.0
        L2_V   = 0.0
        L2_P   = 0.0
        L2_Z   = 0.0
        L2_RHO = 0.0

        vol    = 0.0
        p_sum  = 0.0

        pn_min = HUGE
        pex_min = HUGE

215     format ("# Error (L2 norm)--- ") ! output style
225     format ("# ",A7,5A26)               ! output style
235     format ("# ",F7.5,5ES26.16)         ! output style

        !do c=-NbC,NC
        !    if ( c > 0) then
        !        p_sum = p_sum + P % n(c)*volume(c)
        !        vol  = vol  + volume(c)
        !    end if
        !end do
        !call GloSum(vol)
        !call GloSum(p_sum)

        pn_min = HUGE
        pex_min = HUGE

        do c = 1, NC
            vol     = vol + volume(c)
            pn_min  = min(P % n(c), pn_min)
            pex_min = min(Moin_Problem_p (xc(c), yc(c), time_cur), pex_min)
        end do
        call GloSum(vol)
        call GloMin(pn_min)
        call GloMin(pex_min)


        do c = 1, NC
            L2_U = L2_U + volume(c) * ( U % n(c) - Moin_Problem_u(xc(c), yc(c), time_cur) )**2
            L2_V = L2_V + volume(c) * ( V % n(c) - Moin_Problem_v(xc(c), yc(c), time_cur) )**2
            if ( Moin_Problem_ID.eq.1 .and. abs(time_cur-1.) >= 1.0D-8) then ! to save time
                L2_P = 0.0
            else
                !L2_P = L2_P   + volume(c)*( P % n(c)-p_sum/vol - Moin_Problem_p   (xc(c), yc(c), time_cur) )**2
                L2_P = L2_P + volume(c)*( P % n(c) - pn_min - Moin_Problem_p   (xc(c), yc(c), time_cur) - pex_min )**2
            end if
            L2_Z = L2_Z + volume(c)*( Zmix % n(c) - Moin_Problem_z (xc(c), yc(c), time_cur) )**2
            L2_RHO = L2_RHO + volume(c)*( Rho % n(c) - Moin_Problem_rho (xc(c), yc(c), time_cur) )**2
        end do

        ! averaging over processors

        call GloSum(L2_U)
        call GloSum(L2_V)
        call GloSum(L2_P)
        call GloSum(L2_Z)
        call GloSum(L2_RHO)

        L2_U   = sqrt(L2_U/vol)
        L2_V   = sqrt(L2_V/vol)
        L2_P   = sqrt(L2_P/vol)
        L2_Z   = sqrt(L2_Z/vol)
        L2_RHO = sqrt(L2_RHO/vol)

        if (this < 2) then
            write(*,215)
            write(*,225) "T", "U", "V", "P", "Z", "RHO"
            write(*,235) time_cur, L2_U, L2_V, L2_P, L2_Z, L2_RHO
        end if

    end subroutine Moin_print_error
    !————————————————————————————————————————————————————————————————————————————————————————
    function Moin_Problem_u (x, y, t)
        implicit none

        real,  intent (in)    :: x, y, t
        real :: Moin_Problem_u

        Moin_Problem_u = 0.0

        select case (Moin_Problem_ID)
            case (1)
                Moin_Problem_u = Moin_Problem_1_u(x,y,t)
            case (2)
                Moin_Problem_u = Moin_Problem_2_u(x,y,t)
            case (3)
                Moin_Problem_u = Moin_Problem_3_u(x,y,t)
            case (4)
                Moin_Problem_u = Moin_Problem_4_u(x,y,t)
            case default
                Moin_Problem_u = Moin_Problem_Table_u(x,y,t)
        end select

    end function Moin_Problem_u
    !————————————————————————————————————————————————————————————————————————————————————————
    function Moin_Problem_u_src (x, y, t)
        implicit none

        real,  intent (in)    :: x, y, t
        real :: Moin_Problem_u_src

        Moin_Problem_u_src = 0.0

        select case (Moin_Problem_ID)
            case (1)
                Moin_Problem_u_src = Moin_Problem_1_u_src(x,y,t)
            case (2)
                Moin_Problem_u_src = Moin_Problem_2_u_src(x,y,t)
            case (3)
                Moin_Problem_u_src = Moin_Problem_3_u_src(x,y,t)
            case (4)
                Moin_Problem_u_src = Moin_Problem_4_u_src(x,y,t)
            case default
                Moin_Problem_u_src = Moin_Problem_Table_u_src(x,y,t)
        end select

    end function Moin_Problem_u_src
    !————————————————————————————————————————————————————————————————————————————————————————
    function Moin_Problem_v (x, y, t)
        implicit none

        real,  intent (in)    :: x, y, t
        real :: Moin_Problem_v

        Moin_Problem_v = 0.0

        select case (Moin_Problem_ID)
            case (1)
                Moin_Problem_v = Moin_Problem_1_v(x,y,t)
            case (2)
                Moin_Problem_v = Moin_Problem_2_v(x,y,t)
            case (3)
                Moin_Problem_v = Moin_Problem_3_v(x,y,t)
            case (4)
                Moin_Problem_v = Moin_Problem_4_v(x,y,t)
            case default
                Moin_Problem_v = Moin_Problem_Table_v(x,y,t)
        end select

    end function Moin_Problem_v
    !————————————————————————————————————————————————————————————————————————————————————————
    function Moin_Problem_v_src (x, y, t)
        implicit none

        real,  intent (in)    :: x, y, t
        real :: Moin_Problem_v_src

        Moin_Problem_v_src = 0.0

        select case (Moin_Problem_ID)
            case (1)
                Moin_Problem_v_src = Moin_Problem_1_v_src(x,y,t)
            case (2)
                Moin_Problem_v_src = Moin_Problem_2_v_src(x,y,t)
            case (3)
                Moin_Problem_v_src = Moin_Problem_3_v_src(x,y,t)
            case (4)
                Moin_Problem_v_src = Moin_Problem_4_v_src(x,y,t)
            case default
                Moin_Problem_v_src = Moin_Problem_Table_v_src(x,y,t)
        end select

    end function Moin_Problem_v_src
    !————————————————————————————————————————————————————————————————————————————————————————
    function Moin_Problem_p (x, y, t)
        implicit none

        real,  intent (in)    :: x, y, t
        real :: Moin_Problem_p

        Moin_Problem_p = 0.0

        select case (Moin_Problem_ID)
            case (1)
                Moin_Problem_p = Moin_Problem_1_p(x,y,t)
            case (2)
                Moin_Problem_p = Moin_Problem_2_p(x,y,t)
            case (3)
                Moin_Problem_p = Moin_Problem_3_p(x,y,t)
            case (4)
                Moin_Problem_p = Moin_Problem_4_p(x,y,t)
            case default
                Moin_Problem_p = Moin_Problem_Table_p(x,y,t)
        end select

    end function Moin_Problem_p
    !————————————————————————————————————————————————————————————————————————————————————————
    function Moin_Problem_rho (x, y, t)
        implicit none

        real,  intent (in)    :: x, y, t
        real :: Moin_Problem_rho

        Moin_Problem_rho = 0.0

        select case (Moin_Problem_ID)
            case (1)
                Moin_Problem_rho = Moin_Problem_1_rho(x,y,t)
            case (2)
                Moin_Problem_rho = Moin_Problem_2_rho(x,y,t)
            case (3)
                Moin_Problem_rho = Moin_Problem_3_rho(x,y,t)
            case (4)
                Moin_Problem_rho = Moin_Problem_4_rho(x,y,t)
            case default
                Moin_Problem_rho = Moin_Problem_Table_rho(x,y,t)
        end select

    end function Moin_Problem_rho
    !————————————————————————————————————————————————————————————————————————————————————————
    function Moin_Problem_z (x, y, t)
        implicit none

        real,  intent (in)    :: x, y, t
        real :: Moin_Problem_z

        Moin_Problem_z = 0.0

        select case (Moin_Problem_ID)
            case (1)
                Moin_Problem_z = Moin_Problem_1_z(x,y,t)
            case (2)
                Moin_Problem_z = Moin_Problem_2_z(x,y,t)
            case (3)
                Moin_Problem_z = Moin_Problem_3_z(x,y,t)
            case (4)
                Moin_Problem_z = Moin_Problem_4_z(x,y,t)
            case default
                Moin_Problem_z = Moin_Problem_Table_z(x,y,t)
        end select

    end function Moin_Problem_z
    !————————————————————————————————————————————————————————————————————————————————————————
    function Moin_Problem_z_src (x, y, t)
        implicit none

        real,  intent (in)    :: x, y, t
        real :: Moin_Problem_z_src

        Moin_Problem_z_src = 0.0

        select case (Moin_Problem_ID)
            case (1)
                Moin_Problem_z_src = Moin_Problem_1_z_src(x,y,t)
            case (2)
                Moin_Problem_z_src = Moin_Problem_2_z_src(x,y,t)
            case (3)
                Moin_Problem_z_src = Moin_Problem_3_z_src(x,y,t)
            case (4)
                Moin_Problem_z_src = Moin_Problem_4_z_src(x,y,t)
            case default
                Moin_Problem_z_src = Moin_Problem_Table_z_src(x,y,t)
        end select

    end function Moin_Problem_z_src
    !————————————————————————————————————————————————————————————————————————————————————————
    function Moin_Problem_visc_dyn_of_x (x, y, t) ! it is table only
        implicit none

        real,  intent (in)    :: x, y, t
        real :: Moin_Problem_visc_dyn_of_x
        Moin_Problem_visc_dyn_of_x = 0.0

        select case (Moin_Problem_ID)
            case default
                Moin_Problem_visc_dyn_of_x = Moin_Problem_Table_visc_dyn_of_x(x,y,t) + y*0.0 + t*0.0
        end select

    end function Moin_Problem_visc_dyn_of_x
    !————————————————————————————————————————————————————————————————————————————————————————
    function Moin_Problem_visc_dyn_of_Z (Z) ! it is table only
        implicit none

        real,  intent (in)    :: Z
        real :: Moin_Problem_visc_dyn_of_Z
        Moin_Problem_visc_dyn_of_Z = 0.0

        select case (Moin_Problem_ID)
            case default
                Moin_Problem_visc_dyn_of_Z = Moin_Problem_Table_visc_dyn_of_Z(Z)
        end select

    end function Moin_Problem_visc_dyn_of_Z
    !————————————————————————————————————————————————————————————————————————————————————————
    function Moin_Problem_EOS (Z) ! it is table only yet
        implicit none
        real,  intent (in)    :: Z
        real :: Moin_Problem_EOS
        Moin_Problem_EOS = 0.0

        select case (Moin_Problem_ID)
            case default
                Moin_Problem_EOS = Moin_Problem_Table_EOS(Z)
        end select

    end function Moin_Problem_EOS
    !————————————————————————————————————————————————————————————————————————————————————————
    !---- MOIN SUBMODULE FOR TEST 1
    !————————————————————————————————————————————————————————————————————————————————————————
    subroutine Moin_Problem_1_Init
        use par_Mod
        use vd1d_mms_Mod

        implicit none

225     format (" # Moin problem 1: ", A) ! output style
        if (this < 2) write(*,225) "VIScZmix input value is ignored. "
        if (this < 2) write(*,225) "VIScZmix = VISc"


        call vd1d_mms_init (&
            rho_0, &         !rho_0
            rho_1, &         !rho_1
            4.0, &           !k1
            2.0, &           !k2
            50.0, &          !width
            VIS)             !diff
        if (this < 2) write(*,225) "domain is x:[0 , 2], y:[-3 , 3], z:[-3 , 3] "
        if (this < 2) write(*,225) "successfully initiallized"

    end subroutine Moin_Problem_1_Init
    !————————————————————————————————————————————————————————————————————————————————————————
    function Moin_Problem_1_u(x, y, t)
        use vd1d_mms_Mod

        implicit none

        real,  intent (in)    :: x, y, t
        real :: Moin_Problem_1_u

        Moin_Problem_1_u = vd2d_mms_u (x,y,t)

    end function Moin_Problem_1_u
    !————————————————————————————————————————————————————————————————————————————————————————
    function Moin_Problem_1_u_src(x, y, t)
        use vd1d_mms_Mod

        implicit none

        real,  intent (in)    :: x, y, t
        real :: Moin_Problem_1_u_src

        Moin_Problem_1_u_src = vd2d_mms_u_src (x,y,t)

    end function Moin_Problem_1_u_src
    !————————————————————————————————————————————————————————————————————————————————————————
    function Moin_Problem_1_v(x, y, t)
        use vd1d_mms_Mod

        implicit none

        real,  intent (in)    :: x, y, t
        real :: Moin_Problem_1_v

        Moin_Problem_1_v = vd2d_mms_v (x,y,t)

    end function Moin_Problem_1_v
    !————————————————————————————————————————————————————————————————————————————————————————
    function Moin_Problem_1_v_src(x, y, t)
        use vd1d_mms_Mod

        implicit none

        real,  intent (in)    :: x, y, t
        real :: Moin_Problem_1_v_src

        Moin_Problem_1_v_src = vd2d_mms_v_src (x,y,t)

    end function Moin_Problem_1_v_src
    !————————————————————————————————————————————————————————————————————————————————————————
    function Moin_Problem_1_p(x, y, t)
        !use vd1d_mms_Mod
        use vd1d_mms_Mod_table ! force to use two_column_Module
        implicit none

        real,  intent (in)    :: x, y, t
        real :: Moin_Problem_1_p

        Moin_Problem_1_p = vd2d_mms_p (x,y,t)
    !Moin_Problem_1_p = 0.0

    end function Moin_Problem_1_p
    !————————————————————————————————————————————————————————————————————————————————————————
    function Moin_Problem_1_rho(x, y, t)
        use vd1d_mms_Mod

        implicit none

        real,  intent (in)    :: x, y, t
        real :: Moin_Problem_1_rho

        Moin_Problem_1_rho = vd2d_mms_rho (x,y,t)

    end function Moin_Problem_1_rho
    !————————————————————————————————————————————————————————————————————————————————————————
    function Moin_Problem_1_z(x, y, t)
        use vd1d_mms_Mod

        implicit none

        real,  intent (in)    :: x, y, t
        real :: Moin_Problem_1_z

        Moin_Problem_1_z = vd2d_mms_z (x,y,t)

    end function Moin_Problem_1_z
    !————————————————————————————————————————————————————————————————————————————————————————
    function Moin_Problem_1_z_src(x, y, t)
        use vd1d_mms_Mod

        implicit none

        real,  intent (in)    :: x, y, t
        real :: Moin_Problem_1_z_src

        Moin_Problem_1_z_src = vd2d_mms_z_src (x,y,t)

    end function Moin_Problem_1_z_src
    !————————————————————————————————————————————————————————————————————————————————————————
    !---- MOIN SUBMODULE FOR TEST TABLE MODIFICATION
    !————————————————————————————————————————————————————————————————————————————————————————
    subroutine Moin_Problem_Table_Init( &
        EOS_of_Zmix_in, &
        ZmixVisc_Dyn_of_Zmix_in, &
        Zmix_SRC_of_Zmix_in, &
        vel_Init_in, &
        rho_Init_in, &
        Zmix_Init_in, &
        Visc_dyn_Init_in, &
        pressure_Init_in &
        )

        use par_Mod
        use vd1d_mms_Mod_table
        implicit none
        character (*), optional, intent (in) :: EOS_of_Zmix_in, ZmixVisc_Dyn_of_Zmix_in, Zmix_SRC_of_Zmix_in
        character (*), optional, intent (in) :: vel_Init_in, rho_Init_in, Zmix_Init_in, Visc_dyn_Init_in, pressure_Init_in

225     format (" # Moin problem - Table Modification: ", A) ! output style

        call vd1d_mms_init ( &
            Moin_Problem_ID, &         ! MoinProblemID
            rho_0, &                   ! rho_0
            rho_1, &                   ! rho_1
            VIS, &                     ! diff
            VISZmix, &                 ! visc
            EOS_of_Zmix_in, &          ! EOS file
            ZmixVisc_Dyn_of_Zmix_in, & ! ViscZmix_Dyn file
            Zmix_SRC_of_Zmix_in, &     ! ViscZmix_SRC file
            vel_Init_in, &             ! Vel_ini  file
            rho_Init_in, &             ! rho_ini  file
            Zmix_Init_in, &            ! Zmix_ini file
            Visc_dyn_Init_in, &        ! Visc_Dyn_ini file
            pressure_Init_in &         ! pressure_ini
            )

        if (this < 2) write(*,225) "successfully initiallized"


    end subroutine Moin_Problem_Table_Init
    !————————————————————————————————————————————————————————————————————————————————————————
    function Moin_Problem_Table_u(x, y, t)
        use vd1d_mms_Mod_table

        implicit none

        real,  intent (in)    :: x, y, t
        real :: Moin_Problem_Table_u

        Moin_Problem_Table_u = vd2d_mms_u (x,y,t)

    end function Moin_Problem_Table_u
    !————————————————————————————————————————————————————————————————————————————————————————
    function Moin_Problem_Table_u_src(x, y, t)
        use vd1d_mms_Mod_table

        implicit none

        real,  intent (in)    :: x, y, t
        real :: Moin_Problem_Table_u_src

        Moin_Problem_Table_u_src = vd2d_mms_u_src (x,y,t)

    end function Moin_Problem_Table_u_src
    !————————————————————————————————————————————————————————————————————————————————————————
    function Moin_Problem_Table_v(x, y, t)
        use vd1d_mms_Mod_table

        implicit none

        real,  intent (in)    :: x, y, t
        real :: Moin_Problem_Table_v

        Moin_Problem_Table_v = vd2d_mms_v (x,y,t)

    end function Moin_Problem_Table_v
    !————————————————————————————————————————————————————————————————————————————————————————
    function Moin_Problem_Table_v_src(x, y, t)
        use vd1d_mms_Mod_table

        implicit none

        real,  intent (in)    :: x, y, t
        real :: Moin_Problem_Table_v_src

        Moin_Problem_Table_v_src = vd2d_mms_v_src (x,y,t)

    end function Moin_Problem_Table_v_src
    !————————————————————————————————————————————————————————————————————————————————————————
    function Moin_Problem_Table_p(x, y, t)
        use vd1d_mms_Mod_table

        implicit none

        real,  intent (in)    :: x, y, t
        real :: Moin_Problem_Table_p

        Moin_Problem_Table_p = vd2d_mms_p (x,y,t)

    end function Moin_Problem_Table_p
    !————————————————————————————————————————————————————————————————————————————————————————
    function Moin_Problem_Table_rho(x, y, t)
        use vd1d_mms_Mod_table

        implicit none

        real,  intent (in)    :: x, y, t
        real :: Moin_Problem_Table_rho

        Moin_Problem_Table_rho = vd2d_mms_rho (x,y,t)

    end function Moin_Problem_Table_rho
    !————————————————————————————————————————————————————————————————————————————————————————
    function Moin_Problem_Table_z(x, y, t)
        use vd1d_mms_Mod_table

        implicit none

        real,  intent (in)    :: x, y, t
        real :: Moin_Problem_Table_z

        Moin_Problem_Table_z = vd2d_mms_z (x,y,t)

    end function Moin_Problem_Table_z
    !————————————————————————————————————————————————————————————————————————————————————————
    function Moin_Problem_Table_z_src(x, y, t)
        use vd1d_mms_Mod_table

        implicit none

        real,  intent (in)    :: x, y, t
        real :: Moin_Problem_Table_z_src

        Moin_Problem_Table_z_src = vd2d_mms_z_src (x,y,t)

    end function Moin_Problem_Table_z_src
    !————————————————————————————————————————————————————————————————————————————————————————
    function Moin_Problem_Table_EOS (Z)
        use vd1d_mms_Mod_table

        implicit none

        real,  intent (in)    :: Z
        real :: Moin_Problem_Table_EOS

        Moin_Problem_Table_EOS = moin_table_EOS (Z)

    end function Moin_Problem_Table_EOS
    !————————————————————————————————————————————————————————————————————————————————————————
    function Moin_Problem_Table_visc_dyn_of_Z (Z)
        use vd1d_mms_Mod_table

        implicit none

        real,  intent (in)    :: Z
        real :: Moin_Problem_Table_visc_dyn_of_Z

        Moin_Problem_Table_visc_dyn_of_Z = moin_table_visc_dyn_of_Z (Z)

    end function Moin_Problem_Table_visc_dyn_of_Z
    !————————————————————————————————————————————————————————————————————————————————————————
    function Moin_Problem_Table_visc_dyn_of_x (x, y, t)
        use vd1d_mms_Mod_table

        implicit none

        real,  intent (in)    :: x, y, t
        real :: Moin_Problem_Table_visc_dyn_of_x

        Moin_Problem_Table_visc_dyn_of_x = moin_table_visc_dyn_of_x (x,y,t)

    end function Moin_Problem_Table_visc_dyn_of_x
    !————————————————————————————————————————————————————————————————————————————————————————
    !---- MOIN SUBMODULE FOR TEST 2
    !————————————————————————————————————————————————————————————————————————————————————————
    subroutine Moin_Problem_2_Init
        use par_Mod
        use vd2d_mms_Mod_2

        implicit none

225     format (" # Moin problem 2: ", A) ! output style
        if (this < 2) write(*,225) "VIScZmix input value is ignored. "
        if (this < 2) write(*,225) "VIScZmix = VISc"

        call vd2d_mms_init_2(rho_0, & !rho_0
            rho_1, & !rho_1
            VIS,   & !diff
            VIS,   & !visc
            1.0, & !xvel
            0.5, & !yvel
            0.2, & !ampl
            4.0, & !wnum
            20.0,& !delx
            1.5)   !freq
        if (this < 2) write(*,225) "domain is x:[0 , 3], y:[-0.5 , 0.5], z:[-3 , 3] "
        if (this < 2) write(*,225) "successfully initiallized"


    end subroutine Moin_Problem_2_Init
    !————————————————————————————————————————————————————————————————————————————————————————
    function Moin_Problem_2_u(x, y, t)
        use vd2d_mms_Mod_2

        implicit none

        real,  intent (in)    :: x, y, t
        real :: Moin_Problem_2_u

        Moin_Problem_2_u = vd2d_mms_u (x,y,t)

    end function Moin_Problem_2_u
    !————————————————————————————————————————————————————————————————————————————————————————
    function Moin_Problem_2_u_src(x, y, t)
        use vd2d_mms_Mod_2

        implicit none

        real,  intent (in)    :: x, y, t
        real :: Moin_Problem_2_u_src

        Moin_Problem_2_u_src = vd2d_mms_u_src (x,y,t)

    end function Moin_Problem_2_u_src
    !————————————————————————————————————————————————————————————————————————————————————————
    function Moin_Problem_2_v(x, y, t)
        use vd2d_mms_Mod_2

        implicit none

        real,  intent (in)    :: x, y, t
        real :: Moin_Problem_2_v

        Moin_Problem_2_v = vd2d_mms_v (x,y,t)

    end function Moin_Problem_2_v
    !————————————————————————————————————————————————————————————————————————————————————————
    function Moin_Problem_2_v_src(x, y, t)
        use vd2d_mms_Mod_2

        implicit none

        real,  intent (in)    :: x, y, t
        real :: Moin_Problem_2_v_src

        Moin_Problem_2_v_src = vd2d_mms_v_src (x,y,t)

    end function Moin_Problem_2_v_src
    !————————————————————————————————————————————————————————————————————————————————————————
    function Moin_Problem_2_p(x, y, t)
        use vd2d_mms_Mod_2

        implicit none

        real,  intent (in)    :: x, y, t
        real :: Moin_Problem_2_p

        Moin_Problem_2_p = vd2d_mms_p (x,y,t)

    end function Moin_Problem_2_p
    !————————————————————————————————————————————————————————————————————————————————————————
    function Moin_Problem_2_rho(x, y, t)
        use vd2d_mms_Mod_2

        implicit none

        real,  intent (in)    :: x, y, t
        real :: Moin_Problem_2_rho

        Moin_Problem_2_rho = vd2d_mms_rho (x,y,t)

    end function Moin_Problem_2_rho
    !————————————————————————————————————————————————————————————————————————————————————————
    function Moin_Problem_2_z(x, y, t)
        use vd2d_mms_Mod_2

        implicit none

        real,  intent (in)    :: x, y, t
        real :: Moin_Problem_2_z

        Moin_Problem_2_z = vd2d_mms_z (x,y,t)

    end function Moin_Problem_2_z
    !————————————————————————————————————————————————————————————————————————————————————————
    function Moin_Problem_2_z_src(x, y, t)
        use vd2d_mms_Mod_2

        implicit none

        real,  intent (in)    :: x, y, t
        real :: Moin_Problem_2_z_src

        Moin_Problem_2_z_src = vd2d_mms_z_src (x,y,t)

    end function Moin_Problem_2_z_src
    !————————————————————————————————————————————————————————————————————————————————————————
    !---- MOIN SUBMODULE FOR TEST 3
    !————————————————————————————————————————————————————————————————————————————————————————
    subroutine Moin_Problem_3_Init
        use par_Mod
        use vd2d_mms_Mod_3

        implicit none

225     format (" # Moin problem 3: ", A) ! output style
        if (this < 2) write(*,225) "VIScZmix input value is ignored. "
        if (this < 2) write(*,225) "VIScZmix = VISc"

        call vd2d_mms_init_3 (2.0, & ! wave_num
            2.0, & ! freq
            0.0, & ! uf=vf wave_vel
            rho_0, & ! rho__0
            rho_1, & ! rho__1
            VIS,   & ! visc
            VIS)     ! diff
        if (this < 2) write(*,225) "domain is x:[-1 , 1], y:[-1 , 1], z:[-3 , 3] "
        if (this < 2) write(*,225) "successfully initiallized"

    end subroutine Moin_Problem_3_Init
    !————————————————————————————————————————————————————————————————————————————————————————
    function Moin_Problem_3_u(x, y, t)
        use vd2d_mms_Mod_3

        implicit none

        real,  intent (in)    :: x, y, t
        real :: Moin_Problem_3_u

        Moin_Problem_3_u = vd2d_mms_u (x,y,t)

    end function Moin_Problem_3_u
    !————————————————————————————————————————————————————————————————————————————————————————
    function Moin_Problem_3_u_src(x, y, t)
        use vd2d_mms_Mod_3

        implicit none

        real,  intent (in)    :: x, y, t
        real :: Moin_Problem_3_u_src

        Moin_Problem_3_u_src = vd2d_mms_u_src (x,y,t)

    end function Moin_Problem_3_u_src
    !————————————————————————————————————————————————————————————————————————————————————————
    function Moin_Problem_3_v(x, y, t)
        use vd2d_mms_Mod_3

        implicit none

        real,  intent (in)    :: x, y, t
        real :: Moin_Problem_3_v

        Moin_Problem_3_v = vd2d_mms_v (x,y,t)

    end function Moin_Problem_3_v
    !————————————————————————————————————————————————————————————————————————————————————————
    function Moin_Problem_3_v_src(x, y, t)
        use vd2d_mms_Mod_3

        implicit none

        real,  intent (in)    :: x, y, t
        real :: Moin_Problem_3_v_src

        Moin_Problem_3_v_src = vd2d_mms_v_src (x,y,t)

    end function Moin_Problem_3_v_src
    !————————————————————————————————————————————————————————————————————————————————————————
    function Moin_Problem_3_p(x, y, t)
        use vd2d_mms_Mod_3

        implicit none

        real,  intent (in)    :: x, y, t
        real :: Moin_Problem_3_p

        Moin_Problem_3_p = vd2d_mms_p (x,y,t)

    end function Moin_Problem_3_p
    !————————————————————————————————————————————————————————————————————————————————————————
    function Moin_Problem_3_rho(x, y, t)
        use vd2d_mms_Mod_3

        implicit none

        real,  intent (in)    :: x, y, t
        real :: Moin_Problem_3_rho

        Moin_Problem_3_rho = vd2d_mms_rho (x,y,t)

    end function Moin_Problem_3_rho
    !————————————————————————————————————————————————————————————————————————————————————————
    function Moin_Problem_3_z(x, y, t)
        use vd2d_mms_Mod_3

        implicit none

        real,  intent (in)    :: x, y, t
        real :: Moin_Problem_3_z

        Moin_Problem_3_z = vd2d_mms_z (x,y,t)

    end function Moin_Problem_3_z
    !————————————————————————————————————————————————————————————————————————————————————————
    function Moin_Problem_3_z_src(x, y, t)
        use vd2d_mms_Mod_3

        implicit none

        real,  intent (in)    :: x, y, t
        real :: Moin_Problem_3_z_src

        Moin_Problem_3_z_src = vd2d_mms_z_src (x,y,t)

    end function Moin_Problem_3_z_src
    !————————————————————————————————————————————————————————————————————————————————————————
    !---- MOIN SUBMODULE FOR TEST 4
    !————————————————————————————————————————————————————————————————————————————————————————
    subroutine Moin_Problem_4_Init
        use par_Mod
        use   vd2d_mms_Mod_4_RT

        implicit none

225     format (" # Moin problem 4: ", A) ! output style

        !if (this < 2) write(*,225) "VIScZmix input value is ignored. "
        !if (this < 2) write(*,225) "VIScZmix = VISc"

        call vd2d_mms_init_4 (0.001, & ! magnitude
            0.002, & ! y-width
            9.0,   & !acceleration
            4.0,   & ! freqencies
            14.0,  &
            23.0,  &
            28.0,  &
            33.0,  &
            42.0,  &
            51.0,  &
            59.0   )
        if (this < 2) write(*,225) "domain is x:[-0.5 , 0.5], y:[-0.5 , 0.5], y:[-1.5 , 1.5] "
        if (this < 2) write(*,225) "successfully initiallized"


    end subroutine Moin_Problem_4_Init
    !————————————————————————————————————————————————————————————————————————————————————————
    function Moin_Problem_4_u(x, y, t)
        use   vd2d_mms_Mod_4_RT

        implicit none

        real,  intent (in)    :: x, y, t
        real :: Moin_Problem_4_u

        Moin_Problem_4_u = vd2d_mms_u (x,y,t)

    end function Moin_Problem_4_u
    !————————————————————————————————————————————————————————————————————————————————————————
    function Moin_Problem_4_u_src(x, y, t)
        use   vd2d_mms_Mod_4_RT

        implicit none

        real,  intent (in)    :: x, y, t
        real :: Moin_Problem_4_u_src

        Moin_Problem_4_u_src = vd2d_mms_u_src (x,y,t)

    end function Moin_Problem_4_u_src
    !————————————————————————————————————————————————————————————————————————————————————————
    function Moin_Problem_4_v(x, y, t)
        use   vd2d_mms_Mod_4_RT

        implicit none

        real,  intent (in)    :: x, y, t
        real :: Moin_Problem_4_v

        Moin_Problem_4_v = vd2d_mms_v (x,y,t)

    end function Moin_Problem_4_v
    !————————————————————————————————————————————————————————————————————————————————————————
    function Moin_Problem_4_v_src(x, y, t)
        use   vd2d_mms_Mod_4_RT

        implicit none

        real,  intent (in)    :: x, y, t
        real :: Moin_Problem_4_v_src

        Moin_Problem_4_v_src = vd2d_mms_v_src (x,y,t)

    end function Moin_Problem_4_v_src
    !————————————————————————————————————————————————————————————————————————————————————————
    function Moin_Problem_4_p(x, y, t)
        use   vd2d_mms_Mod_4_RT

        implicit none

        real,  intent (in)    :: x, y, t
        real :: Moin_Problem_4_p

        Moin_Problem_4_p = vd2d_mms_p (x,y,t)

    end function Moin_Problem_4_p
    !————————————————————————————————————————————————————————————————————————————————————————
    function Moin_Problem_4_rho(x, y, t)
        use   vd2d_mms_Mod_4_RT

        implicit none

        real,  intent (in)    :: x, y, t
        real :: Moin_Problem_4_rho

        Moin_Problem_4_rho = vd2d_mms_rho (x,y,t)

    end function Moin_Problem_4_rho
    !————————————————————————————————————————————————————————————————————————————————————————
    function Moin_Problem_4_z(x, y, t)
        use   vd2d_mms_Mod_4_RT

        implicit none

        real,  intent (in)    :: x, y, t
        real :: Moin_Problem_4_z

        Moin_Problem_4_z = vd2d_mms_z (x,y,t)

    end function Moin_Problem_4_z
    !————————————————————————————————————————————————————————————————————————————————————————
    function Moin_Problem_4_z_src(x, y, t)
        use   vd2d_mms_Mod_4_RT

        implicit none

        real,  intent (in)    :: x, y, t
        real :: Moin_Problem_4_z_src

        Moin_Problem_4_z_src = vd2d_mms_z_src (x,y,t)

    end function Moin_Problem_4_z_src
    
end module Moin_Problem_Mod
