!======================================================================!
SUBROUTINE ModFlux(Scale_UVW) !must be placed after NewUVW.f90 and before CalcPS.f90
    !----------------------------------------------------------------------!
    !   Modifies the fluxes
    !   in order to conserve mass inside CalcPS.f90, but not actually
    !   interferes in Ur,Vr,Wr and flux, flux_u directly
    !   The idea behind this is following: int_ right hand side of pressure
    !   equation (delta mass conservation eq.) must be 0.
    !   But between newUVW and CalcPS matrix constucted from
    !   fluxes (flux_u, Ur,Vr,Wr) is saved, so no interfering with those is
    !   allowed. Therefore only impression of real fluxes in calcPS.f90
    !  is modified
    !----------------------------------------------------------------------!
    !------------------------------[Modules]-------------------------------!
    USE allp_mod, only : outflow
    USE all_mod
    USE par_mod
    USE pro_mod
    !----------------------------------------------------------------------!
    IMPLICIT NONE
    !-------------------------------[Locals]-------------------------------!
    INTEGER :: m, s, c, c1, c2
    REAL    :: FL
    REAL    :: dMass, Urs, Vrs, Wrs, Flux_total(1:NC)
    REAL    :: Scale_UVW(Nmat) !why 256
    !--------------------------------[CVS]---------------------------------!
    !  $Id: ModOut.f90,v 1.14 2008/12/10 14:46:52 IUS\mhadziabdic Exp $
    !  $Source: /home/IUS/mhadziabdic/.CVSROOT/T-Rex/Process/ModOut.f90,v $
    !======================================================================!



    !do m=1,Nmat
    !  MassIn(m) = 0.0
    !  MASOUT(m) = 0.0
    !  do s=1,NS
    !    c1=SideC(1,s)
    !    c2=SideC(2,s)
    !    if(c2 < 0) then
    !      if(TypeBC(c2) == OUTFLOW) then !!b.c. on outflow
    !        Rho   % n(c2)   = Rho   % n(c1)
    !        Ur    % n(c2)   = Ur    % n(c1)
    !        Vr    % n(c2)   = Vr    % n(c1)
    !        Wr    % n(c2)   = Wr    % n(c1)
    !        Zmixr % n(c2)   = Zmixr % n(c1)
    !      ! renew flux at boundaries
    !      Urs  = f(s) * Ur  % n(c1) + (1.0-f(s)) * Ur  % n(c2)
    !      Vrs  = f(s) * Vr  % n(c1) + (1.0-f(s)) * Vr  % n(c2)
    !      Wrs  = f(s) * Wr  % n(c1) + (1.0-f(s)) * Wr  % n(c2)
    !      Flux(s) = ( Urs*Sx(s) + Vrs*Sy(s) + Wrs*Sz(s) )
    !      end if
    !      if(TypeBC(c2)  ==  INFLOW   .OR. &
    !        (TypeBC(c2)  ==  PRESSURE .AND. Flux(s) < 0.0) .OR. & ! acts like an inflow
    !        (TypeBC(c2)  ==  CONVECT  .AND. Flux(s) < 0.0) ) then ! acts like an inflow
    !         if(material(c1) == m) MassIn(m) = MassIn(m) + Flux(s)
    !      elseif(TypeBC(c2) == OUTFLOW) then !!b.c. on outflow) then
    !         if(material(c1) == m) MASOUT(m) = MASOUT(m) + Flux(s)
    !      end if
    !    end if
    !  end do
    !  call glosum(MASSIN(m))
    !  call glosum(MASOUT(m))  ! not checked
    !end do
    !dMass = 0.0
    !do m=1,Nmat
    !  do c=1,NC
    !      dMass = dMass + volume(c)/dt*(1.5*Rho % n(c) - 2.0*Rho % o(c) + 0.5*Rho % oo(c))
    !  end do
    !end do
    !call glosum(dMass)
    !do m=1,Nmat
    !  Scale_UVW(m) = -(MASSIN(m)+dMass)/(MASOUT(m)+TINY)
    !end do
    !
    !
    !
    !!! ---- apply b.c.
    !    do s=1, NS
    !      c1=SideC(1,s)
    !      c2=SideC(2,s)
    !      if(c2  < 0) then
    !        if(TypeBC(c2) == OUTFLOW) then
    !          Rho  % n(c2)    = Rho   % n(c1)
    !
    !          Ur    % n(c2)   = Ur     % n(c1)
    !          Vr    % n(c2)   = Vr     % n(c1)
    !          Wr    % n(c2)   = Wr     % n(c1)
    !          Zmixr % n(c2)   = Zmixr  % n(c1)
    !
    !          U     % n(c2)   = Ur     % n(c1) / Rho % n(c2)
    !          V     % n(c2)   = Vr     % n(c1) / Rho % n(c2)
    !          W     % n(c2)   = Wr     % n(c1) / Rho % n(c2)
    !          Zmix  % n(c2)   = Zmixr  % n(c1) / Rho % n(c2)
    !        end if
    !      end if
    !    end do


    !------------------------------------------!
    !     Calculate density change in domain   !
    !------------------------------------------!
    Flux_total = 0.0
    do m=1,Nmat
        do s=1,NS
            c1=SideC(1,s)
            c2=SideC(2,s)
            if(TypeBC(c2) == OUTFLOW) then !!b.c. on outflow
                Rho   % n(c2)   = Rho   % n(c1)
                Ur    % n(c2)   = Ur    % n(c1)
                Vr    % n(c2)   = Vr    % n(c1)
                Wr    % n(c2)   = Wr    % n(c1)
                Zmixr % n(c2)   = Zmixr % n(c1)

                U     % n(c2)   = Ur    % n(c1) / Rho % n(c2)
                V     % n(c2)   = Vr    % n(c1) / Rho % n(c2)
                W     % n(c2)   = Wr    % n(c1) / Rho % n(c2)
                Zmix  % n(c2)   = Zmixr % n(c1) / Rho % n(c2)

            end if
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

    !call glosum(dMass)
    !call glosum(FL)

    do m=1,Nmat
        Scale_UVW(m) = - dMass / FL
    end do
  
    do m=1,Nmat
        if (this < 2) write(*,*) "dMass = ", dMass, " FL = ", FL , " Scale_UVW= ", Scale_UVW(m)
    end do




END SUBROUTINE ModFlux
