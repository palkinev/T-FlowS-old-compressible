module vd2d_mms_mod_2
    implicit none
    real, parameter :: pi = 3.1415926535897932
    real, private :: rho0,rho1,D,mu,uf,vf,a,k,b0,w,r10

    !functions
    PUBLIC :: vd2d_mms_init_2
    PUBLIC :: vd2d_mms_rho
    PUBLIC :: vd2d_mms_rho_of_z
    PUBLIC :: vd2d_mms_z
    PUBLIC :: vd2d_mms_z_src
    PUBLIC :: vd2d_mms_u
    PUBLIC :: vd2d_mms_u_src
    PUBLIC :: vd2d_mms_v
    PUBLIC :: vd2d_mms_v_src
    PUBLIC :: vd2d_mms_p
!————————————————————————————————————————————————————————————————————————————————————————
contains
!————————————————————————————————————————————————————————————————————————————————————————
    subroutine vd2d_mms_init_2 (den0,den1,diff,visc,xvel,yvel,ampl,wnum,delx,freq)
        implicit none
        real, intent (in) :: den0,den1,diff,visc
        real, intent (in) :: xvel,yvel,ampl,wnum,delx,freq
        rho0 = den0! density at z = 0
        rho1 = den1! density at z = 1
        mu = visc ! dynamic viscosity
        D = diff ! diffusion coefficient: rho * alpha_z
        uf = xvel ! x convective velocity
        vf = yvel ! y convective velocity
        a = ampl ! amplitude of rho perturbation
        k = wnum * pi! wave number of rho perturbation
        b0 = delx ! parameter controlling sharpness of transition
        w = freq ! inverse time-scale for ‘‘diffusion"
        r10 = rho1-rho0
    end subroutine vd2d_mms_init_2
    !————————————————————————————————————————————————————————————————————————————————————————
    function vd2d_mms_rho (x,y,t)
        implicit none
        real, intent (in) :: x,y,t
        real  :: vd2d_mms_rho
        real  :: z
        z = vd2d_mms_z (x,y,t)
        vd2d_mms_rho = vd2d_mms_rho_of_z (z)
    end function vd2d_mms_rho
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
        real  :: xt,yt,xi,b,c
        xt = uf * t-x
        yt = vf * t-y
        b = b0 * exp (-w * t)
        xi = a * cos (k * yt)
        c = tanh (b * (xi + xt))
        vd2d_mms_z=(1.0 + c)/(1.0 + c+rho0/rho1 * (1.0-c))
    end function vd2d_mms_z
    !————————————————————————————————————————————————————————————————————————————————————————
    function vd2d_mms_z_src (x,y,t)
        implicit none
        real, intent (in) :: x,y,t
        real  :: vd2d_mms_z_src
        real  :: s0,s1,s2,s3,s4,s5,s6,s7,s8,s9,s10, &
            s11,s12,s13,s14,s15,s16,s17,s18,s19,s20
        s0 = uf * t-x
        s1 = vf * t-y
        s2 = exp (-w * t)
        s3 = 1.0/s2
        s4 = cos (k * s1)
        s5 = sin (-k * s1)
        s6 = a * s4 + s0
        s7 = tanh (b0 * s2 * s6)
        s8 = exp (2.0 * b0 * s2 * s6)
        s9 = 1.0/2.0 + 1.0/2.0 * s7
        s10 = 1.0-s7 ** 2
        s11 = s8 + 1.0
        s12 = 2.0 * (w * a * s4 + w * uf * t-w * x-uf)
        s13 = s9 * r10 + rho0
        s14 = r10 * s7 + rho1 + rho0
        s15 = -s10 * b0 * s2 * r10
        s16 = -s10 * b0 * s2 * a * s5 * k * r10
        s17 = -b0 * w * s2 * s6 + b0 * s2 * (a * s5 * k * vf + uf)
        s18 = (1.0 + s7) * rho1
        s19 = 2.0 * (a * w * s4 * b0 + w * uf * t * b0-x * w * b0)-1.0/s11 * s12 * b0-w * s3 * log (s11)
        s20 = (-s11 ** 2 + s11 + s3 * s2 * s11 * s8) * w-s12 * b0 * s2 * s8
        vd2d_mms_z_src = &
            1.0/2.0 * s10 * s17 * r10 * s18/s14 + s13 * s10 * s17 * rho1/s14-s13 * s18/s14 ** 2 &
            * s10 * b0 * s2 * (-w * s6 + a * s5 * k * vf + uf) * r10 + 1.0/2.0 * s10 * s2 * rho1/s14 * r10 &
            * s19 + 1.0/2.0 * s18/s14 ** 2 * r10 * s19/b0 * s15-s18/s14 * r10 * s20/s11 ** 2 &
            -1.0/2.0 * s10 * b0 * s2 * a * s5 * k * r10 * s18/s14 * vf-s13 * s10 * b0 * s2 * a * s5 * k &
            * rho1/s14 * vf-s13 * s18/s14 ** 2 * vf * s16-D * (-2.0 * s7 * s10 * b0 ** 2 * s2 ** 2 * rho1 &
            /s14 + 2.0 * s10 * b0 * s2 * rho1/s14 ** 2 * s15 + 2.0 * s18/s14 ** 3 * s15 ** 2 + 2.0 &
            * s18/s14 ** 2 * s7 * s10 * b0 ** 2 * s2 ** 2 * r10-2.0 * s7 * s10 * b0 ** 2 * s2 ** 2 * a ** 2 * s5 &
            ** 2 * k ** 2 * rho1/s14-s10 * b0 * s2 * a * s4 * k ** 2 * rho1/s14 + 2.0 * s10 * b0 * s2 * a * s5 &
            * k * rho1/s14 ** 2 * s16 + 2.0 * s18/s14 ** 3 * s16 ** 2 + s18/s14 ** 2 * s10 * b0 * s2 * a * k &
            ** 2 * r10 * (2.0 * a * s7 * s5 ** 2 * s2 * b0 + s4))
    end function vd2d_mms_z_src
    !————————————————————————————————————————————————————————————————————————————————————————
    function vd2d_mms_u (x,y,t)
        implicit none
        real, intent (in) :: x,y,t
        real  :: vd2d_mms_u
        real  :: xt,yt,xi,b
        xt = uf * t-x
        yt = vf * t-y
        b = b0 * exp (-w * t)
        xi = a * cos (k * yt)
        vd2d_mms_u=&
            (rho1-rho0) * (-w * (xi + xt)+(w * (xi + xt)-uf)/(exp (2.0 * b * (xi + xt))+1.0) &
            + 1.0/2.0 * w/b * log (exp (2.0 * b * (xi + xt))+1.0)) &
            /vd2d_mms_rho (x,y,t)
    end function vd2d_mms_u
    !————————————————————————————————————————————————————————————————————————————————————————
    function vd2d_mms_u_src (x,y,t)
        implicit none
        real, intent (in) :: x,y,t
        real  :: vd2d_mms_u_src
        real  :: s0,s1,s2,s3,s4,s5,s6,s7,s8,s9,s10, &
            s11,s12,s13,s14,s15,s16,s17,s18,s19,s20
        s0 = uf * t-x
        s1 = vf * t-y
        s2 = exp (-w * t)
        s3 = 1.0/s2
        s4 = cos (k * s1)
        s5 = sin (-k * s1)
        s6 = a * s4 + s0
        s7 = tanh (b0 * s2 * s6)
        s8 = exp (2.0 * b0 * s2 * s6)
        s9 = 1.0/2.0 + 1.0/2.0 * s7
        s10 = 1.0-s7 ** 2
        s11 = s8 + 1.0
        s12 = 2.0 * (w * a * s4 + w * uf * t-w * x-uf)
        s13 = s9 * r10 + rho0
        s14 = a * s5 * k * vf + uf
        s15 = log (s11) * s11
        s16 = r10 * (-w * a * s4-w * uf * t + w * x)+1.0/2.0 * r10/s11 * (s12 + w * s3 * s15/b0)
        s17=-k * a * s5 * r10 * (w * s11-s12 * b0 * s2 * s8 + s3 * w * s2 * s8 * s11-w * s11 ** 2)/s11 ** 2
        s18 = r10 * (-w * s11 ** 2-s12 * b0 * s2 * s8 + w * s11 + s3 * w * s2 * s8 * s11)/s11 ** 2
        s19 = s15 * w + 2.0 * (-b0 * s2 * s8 * s6 * w + b0 * s2 * s8 * a * s5 * k * vf + b0 * s2 * s8 * uf)
        s20=(s12 * b0 * s2 * s4 * s8 * s11 + w * s4 * s11 ** 3-w * s4 * s11 ** 2-2.0 * s12 * b0 ** 2 * s2 ** 2 &
            * a * s5 ** 2 * s8 * s11 + 2.0 * w * s3 * s2 ** 2 * a * s5 ** 2 * b0 * s8 * s11 ** 2-2.0 * w * s3 * s2 &
            ** 2 * a * s5 ** 2 * s8 ** 2 * b0 * s11 + 4.0 * s12 * b0 ** 2 * s2 ** 2 * a * s5 ** 2 * s8 ** 2-w * s2 * s3 &
            * s4 * s8 * s11 ** 2-4.0 * a * w * s5 ** 2 * b0 * s2 * s8 * s11) * a * k ** 2 * r10/s11 ** 3
        vd2d_mms_u_src=&
            -w * r10 * s14-1.0/2.0/s11 ** 2 * r10 * s12 * (-2.0 * b0 * w * s2 * s6 + 2.0 * b0 &
            * s2 * s14) * s8 + 1.0/2.0/s11 * r10 * (2.0 * w * a * s5 * k * vf + 2.0 * w * uf) &
            + 1.0/2.0 * s3 * w * r10 * s19/s11/b0-2.0 * s16/s13 * s18 + 1.0/2.0 &
            * s16 ** 2/s13 ** 2 * s10 * b0 * s2 * r10 + s17 * vf-mu * ((16.0/3.0/s11 ** 3 * r10 &
            * s12 * b0 ** 2 * s2 ** 2 * s8 ** 2-16.0/3.0/s11 ** 2 * r10 * w * b0 * s2 * s8-8.0 &
            /3.0/s11 ** 2 * r10 * s12 * b0 ** 2 * s2 ** 2 * s8-8.0/3.0 * s3 * w * s2 ** 2 * b0 * s8 &
            * r10 * (-s11 + s8)/s11 ** 2)/s13-4.0/3.0 * s18/s13 ** 2 * s10 * b0 * s2 * r10 &
            + 2.0/3.0 * s16/s13 ** 3 * s10 ** 2 * b0 ** 2 * s2 ** 2 * r10 ** 2 + 4.0/3.0 * s16 &
            /s13 ** 2 * s7 * s10 * b0 ** 2 * s2 ** 2 * r10 + s20/s13 + s17/s13 ** 2 * s10 * b0 * s2 * a * s5 * k &
            * r10 + 1.0/2.0 * s16/s13 ** 3 * s10 ** 2 * b0 ** 2 * s2 ** 2 * a ** 2 * s5 ** 2 * k ** 2 * r10 &
            ** 2 + s16/s13 ** 2 * s7 * s10 * b0 ** 2 * s2 ** 2 * a ** 2 * s5 ** 2 * k ** 2 * r10 + 1.0/2.0 &
            * s16/s13 ** 2 * s10 * b0 * s2 * a * s4 * k ** 2 * r10)
    end function vd2d_mms_u_src
    !————————————————————————————————————————————————————————————————————————————————————————
    function vd2d_mms_v (x,y,t)
        implicit none
        real, intent (in) :: x,y,t
        real  :: vd2d_mms_v
        !useless
        vd2d_mms_v = vf+x-x+y-y+t-t
    end function vd2d_mms_v
    !————————————————————————————————————————————————————————————————————————————————————————
    function vd2d_mms_v_src (x,y,t)
        implicit none
        real, intent (in) :: x,y,t
        real  :: vd2d_mms_v_src
        real  :: s0,s1,s2,s3,s4,s5,s6,s7,s8,s9,s10, &
            s11,s12,s13,s14,s15,s16,s17,s18,s19
        s0 = uf * t-x
        s1 = vf * t-y
        s2 = exp (-w * t)
        s3 = 1.0/s2
        s4 = cos (k * s1)
        s5 = sin (-k * s1)
        s6 = a * s4 + s0
        s7 = tanh (b0 * s2 * s6)
        s8 = exp (2.0 * b0 * s2 * s6)
        s9 = 1.0/2.0 + 1.0/2.0 * s7
        s10 = 1.0-s7 ** 2
        s11 = s8 + 1.0
        s12 = 2.0 * (w * a * s4 + w * uf * t-w * x-uf)
        s13 = s9 * r10 + rho0
        s14=-r10 * (-w * s11 ** 2-s12 * b0 * s2 * s8 + w * s11 + s3 * w * s2 * s11 * s8)/s11 ** 2
        s15 = 2.0 * (-s12 * b0 * s2 * s8 + s11 * w)+s2 * s11 * (s12 * b0-s3 * w * s11 + s3 * w * s8)
        s16 = 2.0 * b0 * s11 * (w * a * s4 + w * uf * t-w * x)-s12 * b0-s3 * w * log (s11) * s11
        s17=-w * s11 ** 2-s12 * b0 * s2 * s8 + s11 * w + s3 * w * s2 * s11 * s8
        s18=-b0 * w * s2 * s6 + b0 * s2 * (a * s5 * k * vf + uf)
        s19 = 1.0/2.0 * (s10 * s18 * r10 * vf-s10 * b0 * s2 * a * s5 * k * r10 * vf ** 2)+s14 * vf
        vd2d_mms_v_src=&
            1.0/12.0 * mu * r10 ** 3 * s16 * b0/s11/s13 ** 3 * s10 ** 2 * s2 ** 2 * a * s5 * k &
            -mu * (1.0/6.0 * s14/s13 ** 2 * b0 * s2 * a * s5 * k * r10-1.0/6.0 * a &
            * s5 * k * r10 ** 2 * s17/s11 ** 2/s13 ** 2 * b0 * s2-1.0/6.0 * r10 ** 2 * s16 &
            * b0/s11/s13 ** 2 * s7 * s2 ** 2 * a * s5 * k) * s10 + s19 + 2.0/3.0 * mu * b0 * s2 &
            * s8 * a * s5 * k * r10 * s15/s11 ** 3/s13
    end function vd2d_mms_v_src
    !————————————————————————————————————————————————————————————————————————————————————————
    function vd2d_mms_p (x,y,t)
        implicit none
        real, intent (in) :: x,y,t
        real  :: vd2d_mms_p
        vd2d_mms_p = 0.0*x*y*t
    !! NOTE: pressure can vary by a constant and still satisfy the mms
    end function vd2d_mms_p

!useless functions here
!subroutine vd1d_mms_init (rho_0,rho_1,k_1,k_2,width,diff)
!implicit none
!real, intent (in) :: rho_0,rho_1,k_1,k_2,width,diff
!k1 = k_1 ! inverse time-scale for amplitude decay
!k2 = k_2 ! inverse time-scale for width decay
!rho0 = rho_0! density at z = 0
!rho1 = rho_1! density at z = 1
!w0 = width ! initial width
!D = diff ! z-diffusivity
!r01 = rho0-rho1
!k12 = k1-k2
!k21 = k2-k1   !diff
!
end module vd2d_mms_mod_2
