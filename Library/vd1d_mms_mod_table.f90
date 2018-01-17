module vd1d_mms_mod_table
    !====================================================================
    ! This module made for Moin problem1 like vd1d_mms_mod,
    ! but for table values: u, rho, viscZmix are in files
    ! It should be read interpolated
    ! via linear interpolation OR spline method
    !
    ! use moin analitical function is file is not set
    !
    !====================================================================
    use two_column_file_with_x_and_y

    implicit none

    type(file), private :: EOS_of_Zmix
    type(file), private :: Visc_Dyn_of_Zmix
    type(file), private :: Zmix_SRC_of_Zmix 
    type(file), private :: rho
    type(file), private :: velocity
    type(file), private :: Zmix
    type(file), private :: VISc_Dyn
    type(file), private :: press

    ! moin 1 problem analitycal parameters
    real, parameter     :: pi = 3.1415926535897932
    integer, private    :: Moin_Problem_ID
    real, private       :: rho0, rho1, k1, k2, w0, D, r01, k12, k21
    real, private       :: k, w, mu, uf

    !functions
    public :: vd1d_mms_init
    public :: moin_table_EOS
    public :: vd2d_mms_u
    public :: vd2d_mms_rho
    public :: vd2d_mms_z
    public :: moin_table_visc_dyn_of_x
    public :: moin_table_visc_dyn_of_Z
    public :: vd2d_mms_z_src
    public :: vd2d_mms_v
    public :: vd2d_mms_p
    public :: vd2d_mms_u_src
    public :: vd2d_mms_v_src
!————————————————————————————————————————————————————————————————————————————————————————
contains
!————————————————————————————————————————————————————————————————————————————————————————
    subroutine vd1d_mms_init ( &
        Moin_Problem_ID_in, &      ! MoinProblemID
        rho_0, &                ! rho_0
        rho_1, &                ! rho_1
        VIS, &                  ! diff
        VISZmix, &              ! visc
        ! candera tables
        EOS_of_Zmix_in, &       ! rho = EOS(Z)
        Visc_Dyn_of_Zmix_in, &  ! visc_dyn(Z)
        Zmix_SRC_of_Zmix_in, &  ! Z_Source(Z)
        ! initial field files
        vel_init_in, &          ! initial u(x)
        rho_init_in, &          ! initial rho(x)
        Zmix_init_in, &         ! initial Z(x)
        Visc_dyn_init_in, &     ! initial VISc_Dyn(x)
        pressure_init_in &     ! initial pressure(x)
        )

        use par_mod
        implicit none

        integer :: Moin_Problem_ID_in
        real    :: rho_0, rho_1, VIS, VISZmix
        character (*), optional, intent (in) :: EOS_of_Zmix_in, Visc_Dyn_of_Zmix_in, Zmix_SRC_of_Zmix_in
        character (*), optional, intent (in) :: vel_init_in, rho_init_in, Zmix_init_in, Visc_dyn_init_in, pressure_init_in

225     format (" # Moin problem table mod: Initialization: ", A)

        Moin_Problem_ID = Moin_Problem_ID_in
        
        ! init parameters like in Moin paper
        select case( Moin_Problem_ID )
        case (10)
            k1 = 4.0
            k2 = 2.0
            k12 = k1-k2
            k21 = k2-k1
            rho0 = rho_0
            rho1 = rho_1
            r01 = rho0-rho1
            w0 = 50.0
            D = VIS 
        case (30)
            k = 2.0
            w = 2.0
            uf = 0.5
            rho0 = rho_0
            rho1 = rho_1
            mu = VIS
            D = VIS
        end select

        !—————————— EOS_of_Zmix_in
        if ( present(EOS_of_Zmix_in) ) then
            call set_filename(EOS_of_Zmix, EOS_of_Zmix_in)

            if ( (this < 2) .and. file_is_set(EOS_of_Zmix) ) then
                write(*,225) "custom EOS is set by file"
                !call print_info(EOS_of_Zmix)
            end if
        end if

        !—————————— Visc_Dyn_of_Zmix_in
        if ( present(Visc_Dyn_of_Zmix_in) ) then
            call set_filename(Visc_Dyn_of_Zmix, Visc_Dyn_of_Zmix_in)

            if ( (this < 2) .and. file_is_set(Visc_Dyn_of_Zmix) ) then
                write(*,225) "custom molecular viscosity of Zmix is set by file"
                !call print_info(Visc_Dyn_of_Zmix)
            end if
        end if

        !—————————— Zmix_SRC_of_Zmix_in
        if ( present(Zmix_SRC_of_Zmix_in) ) then
            call set_filename(Zmix_SRC_of_Zmix, Zmix_SRC_of_Zmix_in)

            if ( (this < 2) .and. file_is_set(Zmix_SRC_of_Zmix) ) then
                write(*,225) "custom Zmix source is set by file"
                !call print_info(Zmix_SRC_of_Zmix)
            end if
        end if

        !——————————Initial velocity
        if ( present(vel_init_in) ) then
            call set_filename(velocity, vel_init_in)

            if ( (this < 2) .and. file_is_set(velocity) ) then
                write(*,225) "custom initial velocity is set by file"
                !call print_info(velocity)
            end if
        end if

        !——————————Initial density
        if ( present(rho_init_in) ) then
            call set_filename(rho, rho_init_in)

            if ( (this < 2) .and. file_is_set(Visc_Dyn_of_Zmix) ) then
                write(*,225) "custom initial density is rho by file"
                !call print_info(rho)
            end if
        end if

        !——————————Initial Zmix
        if ( present(Zmix_init_in) ) then
            call set_filename(Zmix, Zmix_init_in)

            if ( (this < 2) .and. file_is_set(Zmix) ) then
                write(*,225) "custom initial Zmix is set by file"
                !call print_info(Zmix)
            end if
        end if

        !——————————Initial VISc_Dyn
        if ( present(Visc_dyn_init_in) ) then
            call set_filename(VISc_Dyn, Visc_dyn_init_in)

            if ( (this < 2) .and. file_is_set(VISc_Dyn) ) then
                write(*,225) "custom initial dinamic viscosity is set by file"
                !call print_info(VISc_Dyn)
            end if
        end if

        !——————————Initial Pressure
        if ( present(pressure_init_in) ) then
            call set_filename(press, pressure_init_in)

            if ( (this < 2) .and. file_is_set(press) ) then
                write(*,225) "custom initial pressure is set by file"
                !call print_info(press)
            end if
        end if    

    end subroutine vd1d_mms_init
    !————————————————————————————————————————————————————————————————————————————————————————
    function moin_table_EOS (Z)
        implicit none
        real, intent (in) :: Z
        real              :: moin_table_EOS

        if ( file_is_set(EOS_of_Zmix) ) then
            moin_table_EOS = linear_interp(EOS_of_Zmix, Z)
        else
            moin_table_EOS = 1.0 / ( (Z /rho1) + (1.0-Z)/rho0 )
        end if

    end function moin_table_EOS
    !————————————————————————————————————————————————————————————————————————————————————————
    function vd2d_mms_u (x,y,t)
        implicit none
        real, intent (in) :: x,y,t
        real              :: vd2d_mms_u
        real              :: u1
        real  :: xt,yt

        if ( file_is_set(velocity) ) then
            vd2d_mms_u = linear_interp(velocity, x) + y*0.0 + t*0.0
        else
            select case( Moin_Problem_ID )
            case(10)
                u1 = exp(w0 * exp (-k2 * t) * abs(x)) + 0.0*y
                vd2d_mms_u = (2.0 * k2 * abs(x) * r01 * exp(- k1 * t) * u1/(u1**2 + 1.0)+r01 * k12 * &
                    exp(-k12 * t)/w0 * (2.0 * atan(u1) - pi/2.0)) / vd2d_mms_rho(x,y,t)
            case(30)
                xt = x - uf*t
                yt = y - uf*t
                vd2d_mms_u=-w/k/4.0 * cos (k * pi * xt) * sin (k * pi * yt) * sin (w * pi * t)&
                    * (rho1-rho0)/vd2d_mms_rho (x,y,t)
            end select
        end if

    end function vd2d_mms_u
    !————————————————————————————————————————————————————————————————————————————————————————
    function vd2d_mms_rho (x,y,t)
        implicit none
        real, intent (in) :: x,y,t
        real              :: vd2d_mms_rho
        real              :: z

        if ( file_is_set(rho) ) then
          vd2d_mms_rho = linear_interp(rho, x) + y*0.0 + t*0.0
        else
            select case( Moin_Problem_ID )
            case (10)
                z = vd2d_mms_z (abs(x),y,t) + 0.0*y
                vd2d_mms_rho = 1.0/(z/rho1+(1.0-z)/rho0)

            case (30)
                vd2d_mms_rho = 1.0/(vd2d_mms_z (x,y,t)/rho1+(1.0-vd2d_mms_z (x,y,t))/rho0) 
            end select
        end if

    end function vd2d_mms_rho
    !————————————————————————————————————————————————————————————————————————————————————————
    function vd2d_mms_z (x,y,t)
        implicit none
        real, intent (in) :: x,y,t
        real              :: vd2d_mms_z
        real              :: z1,z2
        real  :: xt,yt

        if ( file_is_set(Zmix) ) then
            vd2d_mms_z = linear_interp(Zmix, x) + y*0.0 + t*0.0
        else
            select case( Moin_Problem_ID )
            case(10)
                z1 = exp (-k1 * t) + 0.0*y
                z2 = cosh (w0 * exp (-k2 * t) * abs(x))
                vd2d_mms_z=(z1-z2)/(z1 * (1.0-rho0/rho1)-z2)
            case(30)
                xt = x - uf*t
                yt = y - uf*t

        vd2d_mms_z = &
            (1.0 + sin (k * pi * xt) * sin (k * pi * yt) * cos (w * pi * t))/ &
            (1.0 + rho0/rho1+(1.0-rho0/rho1)&
            * sin (k * pi * xt) * sin (k * pi * yt) * cos (w * pi * t))
            end select
        end if

    end function vd2d_mms_z
    !————————————————————————————————————————————————————————————————————————————————————————
    function moin_table_visc_dyn_of_x (x,y,t)
        implicit none
        real, intent (in) :: x,y,t
        real              :: moin_table_visc_dyn_of_x

        moin_table_visc_dyn_of_x = linear_interp(VISc_Dyn, x) + y*0.0 + t*0.0
    end function moin_table_visc_dyn_of_x
    !————————————————————————————————————————————————————————————————————————————————————————
    function moin_table_visc_dyn_of_Z (Z)
        implicit none
        real, intent (in) :: Z
        real              :: moin_table_visc_dyn_of_Z

        moin_table_visc_dyn_of_Z = linear_interp(Visc_Dyn_of_Zmix, Z)
    end function moin_table_visc_dyn_of_Z
    !————————————————————————————————————————————————————————————————————————————————————————
    function vd2d_mms_z_src (x,y,t)
        implicit none
        real, intent (in) :: x,y,t
        real              :: vd2d_mms_z_src
        real :: s0,s1,s2,s3,s4,s5,s6,s7,s8,s9,s10, &
            s11,s12,s13,s14,s15,s16,s17,s18,s19,s20

        if ( file_is_set(Zmix_SRC_of_Zmix) ) then
            vd2d_mms_z_src = linear_interp(Zmix_SRC_of_Zmix, x) + y*0.0 + t*0.0
        else
            select case( Moin_Problem_ID )
            case(10)
                s0 = exp (-k1 * t) + 0.0*y
                s1 = exp (-k2 * t)
                s2 = sech (w0 * s1 * abs(x))
                s3 = tanh (w0 * s1 * abs(x))
                s4 = exp (w0 * s1 * abs(x))
                s5 = s0/s1
                s6 = s0 * s2
                s7 = rho1 + r01 * s6! rho
                s8 = 2.0 * atan (s4)-1.0/2.0 * pi
                s9 = s6 * s3 * w0
                s10=-s9 * s1 * rho1 + s9 * s1 * rho0
                s11 = 1.0-s3 ** 2
                s12 = s4 ** 2 + 1.0
                s13=-1.0 + s6
                s14=-k1 * s6 + s9 * k2 * s1 * abs(x)
                s15 = s6 * s3 ** 2 * w0 ** 2 * s1 ** 2
                s16 = 2.0 * k2 * r01 * s0 * s4/s12
                s17 = s6 * s11 * w0 ** 2 * s1 ** 2
                s18 = s16 * abs(x) + r01 * k12 * s5/w0 * s8
                s19 = s9 * s1 * s7-s13 * s10
                s20 = 2.0 * s4 * s1 * (s16 * abs(x) * s4 * w0-r01 * k12 * s5)/s12
                vd2d_mms_z_src = &
                    -s14 * rho1-(s16 + s16 * abs(x) * w0 * s1-s20) * s13 * rho1/s7 + s18 * rho1 * s19/s7 ** 2 &
                    -D * rho1 * (-s15 * s7 ** 2 + s17 * s7 ** 2 + 2.0 * s9 * s1 * s10 * s7-2.0 * s13 * s10 ** 2 &
                    -s13 * r01 * s7 * s17 + s13 * r01 * s7 * s15)/s7 ** 3
            case(30)
                s0 = cos (w * pi * t)
                s1 = sin (w * pi * t)
                s2 = cos (k * pi * x)
                s3 = sin (k * pi * x)
                s4 = cos (k * pi * y)
                s5 = sin (k * pi * y)
                s6 = rho1-rho0
                s7 = 1.0 + s3 * s5 * s0
                s8 = (0.5 + 0.5 * s3 * s5 * s0) * s6 + rho0! rho
                s9 = rho1 * s3 * s5 * s0 + rho1-s3 * s5 * s0 * rho0 + rho0
                s10 = rho1 * s2 * k * pi * s5 * s0-s2 * k * pi * s5 * s0 * rho0
                s11 = rho1 * s3 * s4 * k * pi * s0-s3 * s4 * k * pi * s0 * rho0
                s12 = rho0 * s3 * s5 * s1 * w * pi-s3 * s5 * s1 * w * pi * rho1
                s13 = s3 * k ** 2 * pi ** 2 * s5 * s0 * rho0-rho1 * s3 * k ** 2 * pi ** 2 * s5 * s0
                s14 = 1.0/4.0
                s15 = 2.0 * D
                s16 = s15 * pi * pi
                vd2d_mms_z_src=&
                    -rho1 * (s8 * s3 * s5 * s1 * w * pi * s9 ** 2 * k + s8 * s7 * s12 * s9 * k + s14 * s2 ** 2 * pi * s5 ** 2 &
                    * s0 * s1 * w * s6 * s9 ** 2 * k-s14 * s7 * s5 * s1 * w * s6 * s2 * s10 * s9 + s14 * s3 ** 2 * s4 ** 2 * pi &
                    * s0 * s1 * w * s6 * s9 ** 2 * k-s14 * s7 * s3 * s1 * w * s6 * s4 * s11 * s9-s16 * k ** 3 * s3 * s5 * s0 &
                    * s9 ** 2-s15 * k ** 2 * s2 * pi * s5 * s0 * s10 * s9 + s15 * k * s7 * s10 ** 2-s15 * k * s7 * s13 * s9 &
                    -s15 * k ** 2 * s3 * s4 * pi * s0 * s11 * s9 + s15 * k * s7 * s11 ** 2)/s9 ** 3/k
            end select

        end if
    end function vd2d_mms_z_src
    !————————————————————————————————————————————————————————————————————————————————————————
    function vd2d_mms_v (x,y,t)
        implicit none
        real, intent (in) :: x,y,t
        real              :: vd2d_mms_v

        select case( Moin_Problem_ID )
        case(30)
            vd2d_mms_v = vd2d_mms_u (y,x,t)
            !! NOTE: v obtained from u by transposing x <–> y
        case default
            vd2d_mms_v=x*0.0 + y*0.0 + t*0.0
        end select

    end function vd2d_mms_v
    !!————————————————————————————————————————————————————————————————————————————————————————
    function vd2d_mms_p (x,y,t)
        implicit none
        real, intent (in) :: x,y,t
        real              :: vd2d_mms_p
        
        if ( file_is_set(press) ) then
            vd2d_mms_p = linear_interp(press, x) + y*0.0 + t*0.0
        else
            select case( Moin_Problem_ID )
            case (30)
                vd2d_mms_p = vd2d_mms_rho (x,y,t) * vd2d_mms_u (x,y,t) * vd2d_mms_v (x,y,t)/2.0
                !! NOTE: pressure can vary by a constant and still satisfy the mms
            case default
                vd2d_mms_p = x*0.0 + y*0.0 + t*0.0
            end select

        end if

    end function vd2d_mms_p
    !————————————————————————————————————————————————————————————————————————————————————————
    function vd2d_mms_u_src (x,y,t)
        implicit none
        real, intent (in) :: x,y,t
        real              :: vd2d_mms_u_src
        real  :: s0,s1,s2,s3,s4,s5,s6,s7,s8,&
            s9,s10,s11,s12,s13,s14,s15

        select case( Moin_Problem_ID )
        case (30)
            s0 = cos (w * pi * t)
            s1 = sin (w * pi * t)
            s2 = cos (k * pi * x)
            s3 = sin (k * pi * x)
            s4 = cos (k * pi * y)
            s5 = sin (k * pi * y)
            s6 = rho1-rho0
            s7=(0.5 + 0.5 * s3 * s5 * s0) * s6 + rho0! rho
            s8 = 1.0/64.0
            s9 = s8 * 2.0
            s10 = s9 * 2.0
            s11 = s10 * 3.0
            s12 = s10 * 4.0
            s13 = pi * pi/6.0
            s14 = s13 * 2.0
            s15 = s14 * 2.0
            vd2d_mms_u_src=&
                w * s6 * (-s12 * s5 * s0 * w * pi * s2 * s7 ** 3-s11 * s5 ** 2 * s1 ** 2 * w * s6 * s2 * s3 * pi * s7 ** 2 &
                -s9 * s5 ** 3 * s1 ** 2 * w * s6 ** 2 * s2 ** 3 * pi * s0 * s7 + s10 * s4 ** 2 * pi * s1 ** 2 * w * s6 * s2 * s3 &
                * s7 ** 2-s9 * s5 * s1 ** 2 * w * s6 ** 2 * s2 * s3 ** 2 * s4 ** 2 * pi * s0 * s7-s15 * mu * s1 * s2 * k ** 2 &
                * s5 * s7 ** 2 + s15 * mu * s1 * s2 * k ** 2 * s5 ** 2 * s6 * s3 * s0 * s7 + s13 * mu * s1 * s2 ** 3 * k ** 2 &
                * s5 ** 3 * s6 ** 2 * s0 ** 2-s14 * mu * s1 * s2 * k ** 2 * s3 * s6 * s4 ** 2 * s0 * s7 + s13 * mu * s1 * s2 &
                * k ** 2 * s3 ** 2 * s6 ** 2 * s4 ** 2 * s0 ** 2 * s5-s9 * s5 * s1 ** 2 * w * s6 * s3 ** 2 * pi * s4 * s7 ** 2 &
                + s9 * s5 * s1 ** 2 * w * s6 * s2 ** 2 * pi * s4 * s7 ** 2-s8 * s5 ** 2 * s1 ** 2 * w * s6 ** 2 * s2 ** 2 * s3 &
                * s4 * pi * s0 * s7)/k/s7 ** 3
        case default
            vd2d_mms_u_src = x*0.0 + y*0.0 + t*0.0
        end select

    end function vd2d_mms_u_src
    !————————————————————————————————————————————————————————————————————————————————————————
    function vd2d_mms_v_src (x,y,t)
        implicit none
        real, intent (in) :: x,y,t
        real              :: vd2d_mms_v_src
        
        select case( Moin_Problem_ID )
        case(30)
            vd2d_mms_v_src = vd2d_mms_u_src (y,x,t)
            !! NOTE: v obtained from u by transposing x <–> y
        case default
            vd2d_mms_v_src = x*0.0 + y*0.0 + t*0.0
        end select
    end function vd2d_mms_v_src
    !————————————————————————————————————————————————————————————————————————————————————————
    function sech (x)
        implicit none
        real,intent (in) :: x
        real :: sech

        sech = 1.0/cosh (abs(x))
    end function sech
    !————————————————————————————————————————————————————————————————————————————————————————
end module vd1d_mms_mod_table

