module vd2d_mms_mod_4_RT

    implicit none

    real,  parameter :: pi = 3.141592653589793238462
    real,  private :: alpha,delta,g,omega1,omega2,&
        omega3,omega4,omega5,omega6,omega7,omega8

    !functions
    PUBLIC :: vd2d_mms_init_4
    PUBLIC :: vd2d_mms_u
    PUBLIC :: vd2d_mms_rho
    PUBLIC :: vd2d_mms_z
    PUBLIC :: vd2d_mms_z_src
    PUBLIC :: vd2d_mms_v
    PUBLIC :: vd2d_mms_p
    PUBLIC :: vd2d_mms_u_src
    PUBLIC :: vd2d_mms_v_src
    PUBLIC :: y_int

!————————————————————————————————————————————————————————————————————————————————————————
contains
    !————————————————————————————————————————————————————————————————————————————————————————
    subroutine vd2d_mms_init_4 (a1,d1,g1,o1,o2,o3,o4,o5,o6,o7,o8)
        implicit none
        real, intent (in) :: a1,d1,g1,o1,o2,o3,o4,o5,o6,o7,o8
        alpha = a1 ! magnitude
        delta = d1 ! y-width
        g = g1     !acceleration
        !frequencies
        omega1 = o1
        omega2 = o2
        omega3 = o3
        omega4 = o4
        omega5 = o5
        omega6 = o6
        omega7 = o7
        omega8 = o8
    end subroutine vd2d_mms_init_4
    !————————————————————————————————————————————————————————————————————————————————————————
    function vd2d_mms_u (x,y,t)
        implicit none
        real, intent (in) :: x,y,t
        real  :: vd2d_mms_u

        vd2d_mms_u = x*0.0 + y*0.0 + t*0.0
    end function vd2d_mms_u
    !————————————————————————————————————————————————————————————————————————————————————————
    function vd2d_mms_rho (x,y,t)
        implicit none
        real, intent (in) :: x,y,t
        real  :: vd2d_mms_rho

        vd2d_mms_rho =  1.0/( 9.0*vd2d_mms_z(x,y,t) + 1.0)
    end function vd2d_mms_rho
    !————————————————————————————————————————————————————————————————————————————————————————
    function vd2d_mms_z (x,y,t)
        implicit none
        real, intent (in) :: x,y,t
        real  :: vd2d_mms_z

        vd2d_mms_z= 0.5*(1.0 + tanh( (y_int(x)-y)/(2*delta) )) +t*0.0
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

        vd2d_mms_v_src = x*0.0 - g + t*0.0 + y*0.0
    end function vd2d_mms_v_src
    !————————————————————————————————————————————————————————————————————————————————————————

    function y_int (x)
        implicit none
        real, intent (in) :: x
        real  :: y_int

        y_int = -alpha*( cos(omega1*pi*x) + &
            cos(omega2*pi*x) + cos(omega3*pi*x) + cos(omega4*pi*x) + &
            cos(omega5*pi*x) + cos(omega6*pi*x) + cos(omega7*pi*x) + &
            cos(omega8*pi*x)  )
    end function y_int


end module vd2d_mms_mod_4_RT

