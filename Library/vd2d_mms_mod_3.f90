module vd2d_mms_mod_3
    implicit none
    real, parameter :: pi = 3.1415926535897932
    real, private :: k,w,rho0,rho1,D,mu, uf

    !functions
    PUBLIC :: vd2d_mms_init_3
    PUBLIC :: vd2d_mms_rho_of_z
    PUBLIC :: vd2d_mms_z
    PUBLIC :: vd2d_mms_rho
    PUBLIC :: vd2d_mms_u
    PUBLIC :: vd2d_mms_v
    PUBLIC :: vd2d_mms_p
    PUBLIC :: vd2d_mms_z_src
    PUBLIC :: vd2d_mms_u_src
    PUBLIC :: vd2d_mms_v_src
!————————————————————————————————————————————————————————————————————————————————————————
contains
!————————————————————————————————————————————————————————————————————————————————————————
    subroutine vd2d_mms_init_3 (wave_num,freq,wave_vel,rho_0,rho_1,visc,diff)
        implicit none
        real, intent (in) :: wave_num, freq, wave_vel, rho_0, rho_1, visc, diff
        k = wave_num
        w = freq
        rho0 = rho_0
        rho1 = rho_1
        mu = visc
        D = diff
        uf = wave_vel
    end subroutine vd2d_mms_init_3
    !————————————————————————————————————————————————————————————————————————————————————————
    function vd2d_mms_rho_of_z (z)
        implicit none
        real, intent (in) :: z
        real  :: vd2d_mms_rho_of_z
        vd2d_mms_rho_of_z = 1.0/(z/rho1+(1.0-z)/rho0)
    end function vd2d_mms_rho_of_z
    !————————————————————————————————————————————————————————————————————————————————————————
    function vd2d_mms_z (x,y,t)
        implicit none
        real, intent (in) :: x,y,t
        real  :: vd2d_mms_z
        real  :: xt,yt
        xt = x - uf*t
        yt = y - uf*t

        vd2d_mms_z = &
            (1.0 + sin (k * pi * xt) * sin (k * pi * yt) * cos (w * pi * t))/ &
            (1.0 + rho0/rho1+(1.0-rho0/rho1)&
            * sin (k * pi * xt) * sin (k * pi * yt) * cos (w * pi * t))
    end function vd2d_mms_z
    !————————————————————————————————————————————————————————————————————————————————————————
    function vd2d_mms_rho (x,y,t)
        implicit none
        real, intent (in) :: x,y,t
        real  :: vd2d_mms_rho
        vd2d_mms_rho = vd2d_mms_rho_of_z (vd2d_mms_z (x,y,t))
    end function vd2d_mms_rho
    !————————————————————————————————————————————————————————————————————————————————————————
    function vd2d_mms_u (x,y,t)
        implicit none
        real, intent (in) :: x,y,t
        real  :: vd2d_mms_u
        real  :: xt,yt
        xt = x - uf*t
        yt = y - uf*t

        vd2d_mms_u=-w/k/4.0 * cos (k * pi * xt) * sin (k * pi * yt) * sin (w * pi * t)&
            * (rho1-rho0)/vd2d_mms_rho (x,y,t)
    end function vd2d_mms_u
    !————————————————————————————————————————————————————————————————————————————————————————
    function vd2d_mms_v (x,y,t)
        implicit none
        real, intent (in) :: x,y,t
        real  :: vd2d_mms_v
        vd2d_mms_v = vd2d_mms_u (y,x,t)
    !! NOTE: v obtained from u by transposing x <–> y
    end function vd2d_mms_v
    !————————————————————————————————————————————————————————————————————————————————————————
    function vd2d_mms_p (x,y,t)
        implicit none
        real, intent (in) :: x,y,t
        real  :: vd2d_mms_p
        vd2d_mms_p = &
            vd2d_mms_rho (x,y,t) * vd2d_mms_u (x,y,t) * vd2d_mms_v (x,y,t)/2.0
    !! NOTE: pressure can vary by a constant and still satisfy the mms
    end function vd2d_mms_p
    !————————————————————————————————————————————————————————————————————————————————————————
    function vd2d_mms_z_src (x,y,t)
        implicit none
        real, intent (in) :: x,y,t
        real  :: vd2d_mms_z_src
        real  :: s0,s1,s2,s3,s4,s5,s6,s7,s8,&
            s9,s10,s11,s12,s13,s14,s15,s16
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
    end function vd2d_mms_z_src
    !————————————————————————————————————————————————————————————————————————————————————————
    function vd2d_mms_u_src (x,y,t)
        implicit none
        real, intent (in) :: x,y,t
        real  :: vd2d_mms_u_src
        real  :: s0,s1,s2,s3,s4,s5,s6,s7,s8,&
            s9,s10,s11,s12,s13,s14,s15
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
    end function vd2d_mms_u_src
    !————————————————————————————————————————————————————————————————————————————————————————
    function vd2d_mms_v_src (x,y,t)
        implicit none
        real, intent (in) :: x,y,t
        real  :: vd2d_mms_v_src
        vd2d_mms_v_src = vd2d_mms_u_src (y,x,t)
    !! NOTE: v obtained from u by transposing x <–> y
    end function vd2d_mms_v_src
end module vd2d_mms_mod_3
