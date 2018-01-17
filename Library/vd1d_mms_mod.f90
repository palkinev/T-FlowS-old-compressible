module vd1d_mms_mod

    implicit none

    real, parameter :: pi = 3.1415926535897932
    real, private :: rho0,rho1,k1,k2,w0,D,r01,k12,k21

    !functions
    PUBLIC :: vd1d_mms_init
    PUBLIC :: vd2d_mms_u
    PUBLIC :: vd2d_mms_rho
    PUBLIC :: vd2d_mms_z
    PUBLIC :: vd2d_mms_z_src
    PUBLIC :: vd2d_mms_v
    PUBLIC :: vd2d_mms_p
    PUBLIC :: vd2d_mms_u_src
    PUBLIC :: vd2d_mms_v_src
    PUBLIC :: sech
!————————————————————————————————————————————————————————————————————————————————————————
contains
    !————————————————————————————————————————————————————————————————————————————————————————
    subroutine vd1d_mms_init (rho_0,rho_1,k_1,k_2,width,diff)
        implicit none
        real,intent (in) :: rho_0,rho_1,k_1,k_2,width,diff
        k1 = k_1 ! inverse time-scale for amplitude decay
        k2 = k_2 ! inverse time-scale for width decay
        rho0 = rho_0! density at z = 0
        rho1 = rho_1! density at z = 1
        w0 = width ! initial width
        D = diff ! z-diffusivity
        r01 = rho0-rho1
        k12 = k1-k2
        k21 = k2-k1
    end subroutine vd1d_mms_init
    !————————————————————————————————————————————————————————————————————————————————————————
    function vd2d_mms_u (x,y,t)
        implicit none
        real,intent (in) :: x,y,t
        real :: vd2d_mms_u
        real :: u1

        u1 = exp(w0 * exp (-k2 * t) * abs(x)) + 0.0*y
        vd2d_mms_u=(2.0 * k2 * abs(x) * r01 * exp(- k1 * t) * u1/(u1**2 + 1.0)+r01 * k12 * &
            exp(-k12 * t)/w0 * (2.0 * atan(u1) - pi/2.0)) / vd2d_mms_rho(x,y,t)
    end function vd2d_mms_u
    !————————————————————————————————————————————————————————————————————————————————————————
    function vd2d_mms_rho (x,y,t)
        implicit none
        real,intent (in) :: x,y,t
        real :: vd2d_mms_rho
        real :: z

        z = vd2d_mms_z (abs(x),y,t) + 0.0*y
        vd2d_mms_rho = 1.0/(z/rho1+(1.0-z)/rho0)
    end function vd2d_mms_rho
    !————————————————————————————————————————————————————————————————————————————————————————
    function vd2d_mms_z (x,y,t)
        implicit none
        real,intent (in) :: x,y,t
        real :: vd2d_mms_z
        real :: z1,z2

        z1 = exp (-k1 * t) + 0.0*y
        z2 = cosh (w0 * exp (-k2 * t) * abs(x))
        vd2d_mms_z=(z1-z2)/(z1 * (1.0-rho0/rho1)-z2)
    end function vd2d_mms_z
    !————————————————————————————————————————————————————————————————————————————————————————
    function vd2d_mms_z_src (x,y,t)
        implicit none
        real,intent (in) :: x,y,t
        real :: vd2d_mms_z_src
        real :: s0,s1,s2,s3,s4,s5,s6,s7,s8,s9,s10, &
            s11,s12,s13,s14,s15,s16,s17,s18,s19,s20

        !x = abs(x) ! for double mesh

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
        vd2d_mms_z_src=&
            -s14 * rho1-(s16 + s16 * abs(x) * w0 * s1-s20) * s13 * rho1/s7 + s18 * rho1 * s19/s7 ** 2 &
            -D * rho1 * (-s15 * s7 ** 2 + s17 * s7 ** 2 + 2.0 * s9 * s1 * s10 * s7-2.0 * s13 * s10 ** 2 &
            -s13 * r01 * s7 * s17 + s13 * r01 * s7 * s15)/s7 ** 3

    end function vd2d_mms_z_src
    !————————————————————————————————————————————————————————————————————————————————————————
    function vd2d_mms_v (x,y,t)
        implicit none
        real,intent (in) :: x,y,t
        real :: vd2d_mms_v
        vd2d_mms_v=x*0.0 + y*0.0 + t*0.0
    end function vd2d_mms_v
    !!————————————————————————————————————————————————————————————————————————————————————————
    function vd2d_mms_p (x,y,t)
        implicit none
        real,intent (in) :: x,y,t
        real :: vd2d_mms_p
        real :: xl, Pl !, xr, dummy
        integer          :: i, ReadWriteError
        character        :: PressureFile*250, CurrentPath*250,dummmy_str*10

        xl = 0.0
        Pl = 0.0 + 0.0*y
        vd2d_mms_p = 0.0

        PressureFile="pressure/xxxx.dat"
        PressureFile=trim(adjustl(PressureFile))
        write(PressureFile(10:13),'(I4.4)') int(t*800)
        !call getcwd(CurrentPath) !current path
        CurrentPath="/home/palkine/Fortran/Projects/combustion/T-FlowS-comb-Problem-3-CalcPS/Problem1/"

        open(unit=998, file=trim(adjustl(CurrentPath))//PressureFile, &
            iostat=ReadWriteError, access="sequential", status="old", action="read")
        !write(*,*) "Opened file is:", trim(adjustl(CurrentPath))//PressureFile
        if ( ReadWriteError /= 0 ) stop "Error opening pressure file"

        if ( abs(x) > 0.0 ) then
            read(998,*) dummmy_str ! dummy
  
            !if ( x > 0 ) then
            do i = 1, nint(abs(x)*4096/2)+1
                read(998,*) xl,Pl
            end do
            !else
            !  do i = 1, int(-x*4096/2)
            !    read(998,*) xl,Pl
            !  end do
            !end if
            !read(998,*) xr,dummy
            close(998)

            !write(*,*) xl-2/4096, x, xr-2/4096

            vd2d_mms_p = Pl
        end if

        if (isnan(Pl)) then
            write (*,*) "x=", x, "xl=", x, "P=", Pl
            stop '"Pl" is a NaN'
        end if

    end function vd2d_mms_p
    !————————————————————————————————————————————————————————————————————————————————————————

    function vd2d_mms_u_src (x,y,t)
        implicit none
        real,intent (in) :: x,y,t
        real :: vd2d_mms_u_src
        vd2d_mms_u_src = x*0.0 + y*0.0 + t*0.0
    end function vd2d_mms_u_src
    !————————————————————————————————————————————————————————————————————————————————————————
    function vd2d_mms_v_src (x,y,t)
        implicit none
        real,intent (in) :: x,y,t
        real :: vd2d_mms_v_src
        vd2d_mms_v_src = x*0.0 + y*0.0 + t*0.0
    end function vd2d_mms_v_src
    !————————————————————————————————————————————————————————————————————————————————————————

    function sech (x)
        implicit none
        real,intent (in) :: x
        real :: sech

        sech = 1.0/cosh (abs(x))
    end function sech
end module vd1d_mms_mod


!function vd2d_mms_p_old (x,y,t)
!implicit none
!real,intent (in) :: x,y,t
!real :: vd2d_mms_p
!real,allocatable :: xf(:), Pf(:)
!integer          :: i, ReadWriteError
!character        :: PressureFile*250, CurrentPath*250,dummmy_str*10

!allocate(xf(4098)); xf = 0.0
!allocate(Pf(4098)); Pf = 0.0
!vd2d_mms_p = 0.0

!PressureFile="pressure/xxxx.dat"
!PressureFile=trim(adjustl(PressureFile))
!write(PressureFile(10:13),'(I4.4)') int(t*800)
!!call getcwd(CurrentPath) !current path
!CurrentPath="/home/palkine/Fortran/Projects/combustion/T-FlowS-comb-Problem-3-CalcPS/Problem1/"

!open(unit=998, file=trim(adjustl(CurrentPath))//PressureFile, iostat=ReadWriteError, status="old", action="read")
!!write(*,*) "Opened file is:", trim(adjustl(CurrentPath))//PressureFile
!if ( ReadWriteError /= 0 ) stop "Error opening pressure file"


!  read(998,*) dummmy_str ! dummy

!  do i = 1, 4098 !int(x*4096/2)+1
!    read(998,*) xf(i),Pf(i)
!  end do
!  close(998)

!if ( x == 0 ) then
!  vd2d_mms_p = 0.0*y !+ xr*0  + dummy*0
!  return
!else if ( x == 2 ) then
!  vd2d_mms_p = Pf(4098)
!  return
!end if

!do i = 2, 4097
!  if ( xf(i)-2/4096 <= x .and. x <= xf(i+1)-2/4096 ) then !caught the cell
!    vd2d_mms_p = Pf(i)
!    return
!  end if
!end do

!deallocate(xf)
!deallocate(Pf)

!end function vd2d_mms_p_old
