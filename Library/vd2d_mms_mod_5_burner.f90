module vd2d_mms_mod_5_burner

    implicit none

    real,  parameter :: pi = 3.1415926535897932
    real,  private :: dens0,dens1,alpha,delta,&
        omega1,omega2,omega3,omega4,omega5,omega6,omega7,omega8

    !functions
    PUBLIC :: vd2d_mms_init_5
    PUBLIC :: vd2d_mms_u
    PUBLIC :: vd2d_mms_rho
    PUBLIC :: vd2d_mms_z
    PUBLIC :: vd2d_mms_z_src
    PUBLIC :: vd2d_mms_v
    PUBLIC :: vd2d_mms_p
    PUBLIC :: vd2d_mms_u_src
    PUBLIC :: vd2d_mms_v_src
    PUBLIC :: x_int
!————————————————————————————————————————————————————————————————————————————————————————
contains
    !————————————————————————————————————————————————————————————————————————————————————————
    subroutine vd2d_mms_init_5 (den0,den1,a1,d1,o1,o2,o3,o4,o5,o6,o7,o8)
        implicit none
        real, intent (in) :: den0,den1,a1,d1,o1,o2,o3,o4,o5,o6,o7,o8
        dens0 = den0 ! high-density liquid
        dens1 = den1 ! low-density liquid
        alpha = a1 ! magnitude
        delta = d1 ! x-width
        !frequencies
        omega1 = o1
        omega2 = o2
        omega3 = o3
        omega4 = o4
        omega5 = o5
        omega6 = o6
        omega7 = o7
        omega8 = o8
    end subroutine vd2d_mms_init_5
    !————————————————————————————————————————————————————————————————————————————————————————
    function vd2d_mms_u (x,y,t) !t here is z !!!
        implicit none
        real, intent (in) :: x,y,t
        real  :: vd2d_mms_u

        !vd2d_mms_u = 0.5
        !if (sqrt(x**2 + y**2) < 1.0) then
        !    vd2d_mms_u = 0.75 + 0.25*tanh( (x_int(y)-x+0.7)/(2*delta) ) +t*0.0+x*0.0
        !end if

        vd2d_mms_u = 0.5
        if (sqrt(x**2 + y**2) < 1.0 .and. t < 0.0) then
            vd2d_mms_u = 1.0
        end if


    end function vd2d_mms_u
    !————————————————————————————————————————————————————————————————————————————————————————
    function vd2d_mms_rho (x,y,t)
        implicit none
        real, intent (in) :: x,y,t
        real  :: vd2d_mms_rho

        vd2d_mms_rho =  1.0/( (vd2d_mms_z(x,y,t)/dens1) + (1.0-vd2d_mms_z(x,y,t))/dens0 )
    end function vd2d_mms_rho
    !————————————————————————————————————————————————————————————————————————————————————————
    function vd2d_mms_z (x,y,t) !t here is z !!!
        implicit none
        real, intent (in) :: x,y,t
        real  :: vd2d_mms_z

        !vd2d_mms_z = 1.0
        !if (sqrt(x**2 + y**2) < 1.0) then
        !    vd2d_mms_z = 0.5*(1.0 + tanh( (x-x_int(y)-0.7)/(2*delta) )) +t*0.0+x*0.0
        !end if

        vd2d_mms_z = 1.0
        if (sqrt(x**2 + y**2) < 1.0 .and. t < 0.0) then
            vd2d_mms_z = 0.0
        end if


    end function vd2d_mms_z
    !————————————————————————————————————————————————————————————————————————————————————————
    function vd2d_mms_z_src (x,y,t)
        implicit none
        real, intent (in) :: x,y,t
        real  :: vd2d_mms_z_src

        vd2d_mms_z_src=x*0.0 + y*0.0 + t*0.0

    end function vd2d_mms_z_src
    !————————————————————————————————————————————————————————————————————————————————————————
    function vd2d_mms_v (x,y,t)
        implicit none
        real, intent (in) :: x,y,t
        real  :: vd2d_mms_v

        vd2d_mms_v=x*0.0 + y*0.0 + t*0.0
    end function vd2d_mms_v
    !!————————————————————————————————————————————————————————————————————————————————————————
    function vd2d_mms_p (x,y,t)
        implicit none
        real, intent (in) :: x,y,t
        real  :: vd2d_mms_p

        vd2d_mms_p=x*0.0 + y*0.0 + t*0.0

    end function vd2d_mms_p
    !————————————————————————————————————————————————————————————————————————————————————————

    function vd2d_mms_u_src (x,y,t)
        implicit none
        real, intent (in) :: x,y,t
        real  :: vd2d_mms_u_src

        vd2d_mms_u_src = x*0.0 + y*0.0 + t*0.0
    end function vd2d_mms_u_src
    !————————————————————————————————————————————————————————————————————————————————————————
    function vd2d_mms_v_src (x,y,t)
        implicit none
        real, intent (in) :: x,y,t
        real  :: vd2d_mms_v_src

        vd2d_mms_v_src = x*0.0 + t*0.0 + y*0.0
    end function vd2d_mms_v_src
    !————————————————————————————————————————————————————————————————————————————————————————

    function x_int (y)
        implicit none
        real, intent (in) :: y
        real  :: x_int

        x_int = -alpha*( cos(omega1*pi*y) + &
            cos(omega2*pi*y) + cos(omega3*pi*y) + cos(omega4*pi*y) + &
            cos(omega5*pi*y) + cos(omega6*pi*y) + cos(omega7*pi*y) + &
            cos(omega8*pi*y)  ) !maybe reference paper has a mistake in formulation and do not have 8th mode
    end function x_int

end module vd2d_mms_mod_5_burner

