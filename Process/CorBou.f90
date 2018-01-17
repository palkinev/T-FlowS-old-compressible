!======================================================================!
SUBROUTINE CorBou(n)
    !----------------------------------------------------------------------!
    !   Extrapoloate variables on the boundaries where needed              !
    !----------------------------------------------------------------------!
    !------------------------------[Modules]-------------------------------!
    USE allp_mod, only : inflow, outflow, pressure, wallfl
    USE all_mod
    USE pro_mod
    USE rans_mod
    USE vd1d_mms_mod   !# problem_1
    !USE vd2d_mms_mod_2 !# problem_2
    !USE vd2d_mms_mod_3 !# problem_3
    !----------------------------------------------------------------------!
    IMPLICIT NONE
    !-----------------------------[Parameters]-----------------------------!
    INTEGER       :: n
    !-------------------------------[Locals]-------------------------------!
    INTEGER :: c1, c2, s
    !--------------------------------[CVS]---------------------------------!
    !  $Id: CorBou.f90,v 1.20 2008/12/11 09:38:58 IUS\mhadziabdic Exp $
    !  $Source: /home/IUS/mhadziabdic/.CVSROOT/T-Rex/Process/CalBou.f90,v $
    !======================================================================!

    do s=1,NS
        c1=SideC(1,s)
        c2=SideC(2,s)

        !----------------------------------------------------!
        !     Outflow (and inflow, if needed) boundaries     !
        !----------------------------------------------------!

        !---- on the boundary perform the extrapolation
        if(c2  < 0) then

            if( TypeBC(c2) == INFLOW ) then

                !Ur    %  n(c2) = vd2d_mms_u   (xc(c2), yc(c2), (n  )*dt) * vd2d_mms_rho (xc(c2), yc(c2), (n  )*dt)
                !Vr    %  n(c2) = 0.0
                !Wr    %  n(c2) = 0.0
                !Zmixr %  n(c2) = vd2d_mms_z   (xc(c2), yc(c2), (n  )*dt) * vd2d_mms_rho (xc(c2), yc(c2), (n  )*dt)

                Rho   %  n(c2) = vd2d_mms_rho (xc(c2), yc(c2), (n  )*dt)
                U     %  n(c2) = vd2d_mms_u   (xc(c2), yc(c2), (n  )*dt)
                V     %  n(c2) = 0.0
                W     %  n(c2) = 0.0
                Zmix  %  n(c2) = vd2d_mms_z   (xc(c2), yc(c2), (n  )*dt)

            elseif( TypeBC(c2) == OUTFLOW.or.TypeBC(c2) == PRESSURE.or.TypeBC(c2) == WALLFL ) then

                !Ur      % n(c2)   = Ur    % n(c1)
                !Vr      % n(c2)   = Vr    % n(c1)
                !Wr      % n(c2)   = Wr    % n(c1)
                !Zmixr   % n(c2)   = Zmixr % n(c1)


                Rho   % n(c2)   = Rho   % n(c1)
        
                !U     % n(c2)   = Ur     % n(c1) / Rho % n(c2)
                !V     % n(c2)   = Vr     % n(c1) / Rho % n(c2)
                !W     % n(c2)   = Wr     % n(c1) / Rho % n(c2)
                !Zmix  % n(c2)   = Zmixr  % n(c1) / Rho % n(c2)

            !P     % n(c2)   = P     % n(c1)
            !PP    % n(c2)   = PP    % n(c1)

            end if
        end if
    end do

END SUBROUTINE CorBou
