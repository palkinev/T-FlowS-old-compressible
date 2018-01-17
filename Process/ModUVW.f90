!======================================================================!
SUBROUTINE ModUVW() !must be placed after NewUVW.f90 and before CalcPS.f90
    !----------------------------------------------------------------------!
    !   Modifies the fluxes and Ur,Vr,Wr by multiplying with a factor
    !   at whole domain in order to conserve mass.
    !----------------------------------------------------------------------!
    !------------------------------[Modules]-------------------------------!
    USE allp_mod, only : inflow, outflow
    USE all_mod
    USE pro_mod
    !----------------------------------------------------------------------!
    IMPLICIT NONE
    !-------------------------------[Locals]-------------------------------!
    INTEGER :: m, s, c, c1, c2
    REAL    :: fac(256) , FL
    REAL    :: dMass, Urs, Vrs, Wrs, Flux_total(1:NC)
    REAL    :: RHOs
    !--------------------------------[CVS]---------------------------------!
    !  $Id: ModOut.f90,v 1.14 2008/12/10 14:46:52 IUS\mhadziabdic Exp $
    !  $Source: /home/IUS/mhadziabdic/.CVSROOT/T-Rex/Process/ModOut.f90,v $
    !======================================================================!

    !------------------------------------------!
    !     Calculate density change in domain   !
    !------------------------------------------!

    do s=1, NS
        c1=SideC(1,s)
        c2=SideC(2,s)
        if(c2  < 0) then
            if(TypeBC(c2) == OUTFLOW) then
                Rho  % n(c2)   = Rho   % n(c1)

                Ur    % n(c2)   = Ur     % n(c1)
                Vr    % n(c2)   = Vr     % n(c1)
                Wr    % n(c2)   = Wr     % n(c1)
                Zmixr % n(c2)   = Zmixr  % n(c1)
          
                U     % n(c2)   = Ur     % n(c1) / Rho % n(c2)
                V     % n(c2)   = Vr     % n(c1) / Rho % n(c2)
                W     % n(c2)   = Wr     % n(c1) / Rho % n(c2)
                Zmix  % n(c2)   = Zmixr  % n(c1) / Rho % n(c2)
            end if
        end if
    end do


    Flux_total = 0.0
    do m=1,Nmat
        do s=1,NS
            c1=SideC(1,s)
            c2=SideC(2,s)
            Urs  = f(s) * Ur  % n(c1) + (1.0-f(s)) * Ur  % n(c2)
            Vrs  = f(s) * Vr  % n(c1) + (1.0-f(s)) * Vr  % n(c2)
            Wrs  = f(s) * Wr  % n(c1) + (1.0-f(s)) * Wr  % n(c2)

            Flux_total(c1) = Flux_total(c1) + ( Urs*Sx(s) + Vrs*Sy(s) + Wrs*Sz(s) )
            if(c2  > 0) Flux_total(c2) = Flux_total(c2) - ( Urs*Sx(s) + Vrs*Sy(s) + Wrs*Sz(s) )
        end do
    end do
  
    FL    = 0.0
    dMass = 0.0
    do m=1,Nmat
        do c=1,NC
            dMass = dMass + volume(c)/dt*(1.5*Rho % n(c) - 2.0*Rho % o(c) + 0.5*Rho % oo(c))

            FL    = FL + Flux_total(c)
        end do
    end do

    fac(1) = - dMass / FL

    write(*,*) "dMass = ", dMass, " FL = ", FL , "fac=", fac(1)

  
        !P %mean(:) = 0.0
        !do s=1,NS
        !  c1=SideC(1,s)
        !  c2=SideC(2,s)
        !  if(c2  > 0 .or. c2  < 0 .and. TypeBC(c2) == BUFFER) then
        !                P %mean(c1) = P %mean(c1)-Flux(s)
        !    if(c2  > 0) P %mean(c2) = P %mean(c2)+Flux(s)
        !  else
        !                P %mean(c1) = P %mean(c1)-Flux(s)
        !  end if
        !end do
        !do c=1,NC     ! - d (rho^n+1) /dt
        !  ! 1)
        !  P %mean(c) = - P %mean(c) + volume(c)/dt*(1.5*Rho % n(c) - 2.0*Rho % o(c) + 0.5*Rho % oo(c) )
        !  ! 2)
        !  !P % mean(c) = P % mean(c) - volume(c)/dt*(1.5*Rho % n(c) - 2.0*Rho % o(c) )
        !  !P % mean(c) = P % mean(c) + volume(c)/dt*( 0.5*Rho % oo(c) )
        !end do
        !  open(996, FILE='moduvw_1.dat')
        !  write(*,*) 'Saving to moduvw_1.dat'
        !  do c=-NbC,NC
        !    if ( c .ne. 0 .AND. zc(c)>=-1.0 .AND. zc(c)<=1.0 .AND.  yc(c)<=1.0/75.0 .AND. yc(c)>=0.0 ) then
        !    write(996,'(5ES26.16E3)') xc(c), yc(c), &
        !    Ur % n(c), volume(c)/dt*(1.5*Rho % n(c) - 2.0*Rho % o(c) + 0.5*Rho % oo(c) ), P % mean(c)
        !    end if
        !  end do
        !  call wait
        !  close(996)
        !  !call exit(1)


    do c=-NbC,NC
        if ( TypeBC(c)/=INFLOW) then
            Ur % n(c) = Ur % n(c) * fac(1)
            Vr % n(c) = Vr % n(c) * fac(1)
            Wr % n(c) = Wr % n(c) * fac(1)
        end if
    end do



    !---- Update Flux & Flux_u
    do s=1,NS
        c1=SideC(1,s)
        c2=SideC(2,s)

        if  ( TypeBC(c2) /=INFLOW) then
            RhoS = f(s) * Rho  % n(c1) + (1.0-f(s)) * Rho % n(c2)
            Urs  = f(s) * Ur   % n(c1) + (1.0-f(s)) * Ur  % n(c2)
            Vrs  = f(s) * Vr   % n(c1) + (1.0-f(s)) * Vr  % n(c2)
            Wrs  = f(s) * Wr   % n(c1) + (1.0-f(s)) * Wr  % n(c2)

            Flux  (s) = Urs*Sx(s) + Vrs*Sy(s) + Wrs*Sz(s)
            Flux_u(s) = Flux(s)/RhoS

        end if
    end do


    !P %mean(:) = 0.0
    !do s=1,NS
    !  c1=SideC(1,s)
    !  c2=SideC(2,s)
    !  if(c2  > 0 .or. c2  < 0 .and. TypeBC(c2) == BUFFER) then
    !                P %mean(c1) = P %mean(c1)-Flux(s)
    !    if(c2  > 0) P %mean(c2) = P %mean(c2)+Flux(s)
    !  else
    !                P %mean(c1) = P %mean(c1)-Flux(s)
    !  end if
    !end do
    !do c=1,NC     ! - d (rho^n+1) /dt
    !  ! 1)
    !  P %mean(c) = - P %mean(c) + volume(c)/dt*(1.5*Rho % n(c) - 2.0*Rho % o(c) + 0.5*Rho % oo(c) )
    !  ! 2)
    !  !P % mean(c) = P % mean(c) - volume(c)/dt*(1.5*Rho % n(c) - 2.0*Rho % o(c) )
    !  !P % mean(c) = P % mean(c) + volume(c)/dt*( 0.5*Rho % oo(c) )
    !end do
    !  open(996, FILE='moduvw_2.dat')
    !  write(*,*) 'Saving to moduvw_2.dat'
    !  do c=-NbC,NC
    !    if ( c .ne. 0 .AND. zc(c)>=-1.0 .AND. zc(c)<=1.0 .AND.  yc(c)<=1.0/75.0 .AND. yc(c)>=0.0 ) then
    !    write(996,'(5ES26.16E3)') xc(c), yc(c), &
    !    Ur % n(c), volume(c)/dt*(1.5*Rho % n(c) - 2.0*Rho % o(c) + 0.5*Rho % oo(c) ), P % mean(c)
    !    end if
    !  end do
    !  call wait
    !  close(996)
    !  !call exit(1)

END SUBROUTINE ModUVW
