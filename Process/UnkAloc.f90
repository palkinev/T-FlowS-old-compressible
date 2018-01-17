!======================================================================!
SUBROUTINE UnkAloc
    !----------------------------------------------------------------------!
    ! Allocates the memory for unknowns. It is called either from LoaRes   !
    ! or from Processor.                                                   !
    !------------------------------[Modules]-------------------------------!
    USE all_mod
    USE pro_mod
    USE les_mod
    USE par_mod
    USE rans_mod
    !----------------------------------------------------------------------!
    IMPLICIT NONE
    !--------------------------------[CVS]---------------------------------!
    !  $Id: UnkAloc.f90,v 1.28 2008/11/19 14:55:56 IUS\mhadziabdic Exp $
    !  $Source: /home/IUS/mhadziabdic/.CVSROOT/T-Rex/Process/UnkAloc.f90,v $
    !======================================================================!

    ! rho
    allocate (Rho  % n(-NbC:NC));    Rho % n    = 1.0
    allocate (Rho  % mean(-NbC:NC)); Rho % mean = 0.0
    allocate (Rho  % o     (   1:NC) )
    allocate (Rho  % oo    (   1:NC) )

    !Zmix
    allocate (Zmix  % n(-NbC:NC));    Zmix % n  = 0.0
    allocate (Zmix  % mean(-NbC:NC)); Zmix % mean = 0.0
    allocate (Zmix  % o     (   1:NC) )
    allocate (Zmix  % oo    (   1:NC) )
    allocate (Zmix  % C     (   1:NC) )
    allocate (Zmix  % X     (   1:NC) )
    allocate (Zmix  % source(   1:NC) )

    ! U
    allocate (U  % n(-NbC:NC));    U % n    = 0.0
    allocate (U  % mean(-NbC:NC)); U % mean = 0.0
    allocate (U  % o     (   1:NC) )
    allocate (U  % oo    (   1:NC) )
    allocate (U  % C     (   1:NC) )
    allocate (U  % X     (   1:NC) )
    allocate (U  % source(   1:NC) )

    ! V
    allocate (V  % n(-NbC:NC));    V % n    = 0.0
    allocate (V  % mean(-NbC:NC)); V % mean = 0.0
    allocate (V  % o     (   1:NC) )
    allocate (V  % oo    (   1:NC) )
    allocate (V  % C     (   1:NC) )
    allocate (V  % X     (   1:NC) )
    allocate (V  % source(   1:NC) )

    ! W
    allocate (W  % n(-NbC:NC));    W % n    = 0.0
    allocate (W  % mean(-NbC:NC)); W % mean = 0.0
    allocate (W  % o     (   1:NC) )
    allocate (W  % oo    (   1:NC) )
    allocate (W  % C     (   1:NC) )
    allocate (W  % X     (   1:NC) )
    allocate (W  % source(   1:NC) )

    !dynamic viscocity
    allocate (VISc_Dyn(-NbC:NC)); VISc_Dyn = 0.0

    !velocity gradients
    allocate (Ux(-NbC:NC)); Ux=0.0
    allocate (Uy(-NbC:NC)); Uy=0.0
    allocate (Uz(-NbC:NC)); Uz=0.0

    allocate (Vx(-NbC:NC)); Vx=0.0
    allocate (Vy(-NbC:NC)); Vy=0.0
    allocate (Vz(-NbC:NC)); Vz=0.0

    allocate (Wx(-NbC:NC)); Wx=0.0
    allocate (Wy(-NbC:NC)); Wy=0.0
    allocate (Wz(-NbC:NC)); Wz=0.0

    !dynamic VIScZmix
    allocate (VIScZmix_Dyn(-NbC:NC)); VIScZmix_Dyn = 0.0

    ! derivative of Zmix for diffusion term in CalcZmix
    allocate (Zmixx(-NbC:NC)); Zmixx=0.0
    allocate (Zmixy(-NbC:NC)); Zmixy=0.0
    allocate (Zmixz(-NbC:NC)); Zmixz=0.0

    !pressure
    allocate (P  % n(-NbC:NC));    P % n    = 0.0
    allocate (P  % mean(-NbC:NC)); P % mean = 0.0
    allocate (PP % n(-NbC:NC));    PP % n   = 0.0

    !pressure gradients
    allocate (Px(-NbC:NC)); Px=0.
    allocate (Py(-NbC:NC)); Py=0.
    allocate (Pz(-NbC:NC)); Pz=0.

    allocate (Ps(NS)); Ps=0.;

    allocate (PHIx(-NbC:NC)); PHIx=0.
    allocate (PHIy(-NbC:NC)); PHIy=0.
    allocate (PHIz(-NbC:NC)); PHIz=0.
    allocate (PHIside(NS)); PHIside=0.

    allocate (PHImax(-NbC:NC)); PHImax=0.
    allocate (PHImin(-NbC:NC)); PHImin=0.

    allocate (G(6,NC)); G=0


    allocate (Flux(NS));     Flux=0.
    allocate (Flux_u(NS));   Flux_u=0.

    allocate (PdropX(Nmat)); PdropX=0.0
    allocate (PdropY(Nmat)); PdropY=0.0
    allocate (PdropZ(Nmat)); PdropZ=0.0

    allocate (Utau(Nmat));   Utau=0.0
    allocate (Vtau(Nmat));   Vtau=0.0
    allocate (Wtau(Nmat));   Wtau=0.0

    allocate (FLUXx(Nmat));  FLUXx=0.0
    allocate (FLUXy(Nmat));  FLUXy=0.0
    allocate (FLUXz(Nmat));  FLUXz=0.0

    allocate (FLUXoX(Nmat)); FLUXoX=0.0
    allocate (FLUXoY(Nmat)); FLUXoY=0.0
    allocate (FLUXoZ(Nmat)); FLUXoZ=0.0

    allocate (Ubulk(Nmat));  Ubulk=0.0
    allocate (Vbulk(Nmat));  Vbulk=0.0
    allocate (Wbulk(Nmat));  Wbulk=0.0

    allocate (MassIn(Nmat)); MassIn=0.0
    allocate (MasOut(Nmat)); MasOut=0.0

    allocate (BadForG(NC));  BadForG = .FALSE.
    allocate (NumGood(NC));  NumGood = 0
    allocate (NumNeig(NC));  NumNeig = 0

    allocate (near(-NbC:NC));  near  = 0
    allocate (VISwall(-NbC:NC)); VISwall =0.0

    !---- variables for temperature
    if(HOT==YES) then
        allocate (T % n(-NbC:NC)); T % n=0.
        allocate (T % o(NC));      T % o=0.
        allocate (T % oo(NC));     T % oo=0.
        allocate (T % C(NC));      T % C=0.
        allocate (T % X(NC));      T % X=0.
        allocate (T % q(-NbC:-1)); T % q=0.
    end if

    !---- variables defined in les_mod.h90:
    if(SIMULA == LES) then
        if(MODE == DYN.or.MODE==MIX) then

            allocate (Cdyn(-NbC:NC)); Cdyn = 0.0

            allocate (Rho_filt(-NbC:NC));    Rho_filt = 0.

            allocate (U_filt(-NbC:NC));      U_filt  = 0.
            allocate (V_filt(-NbC:NC));      V_filt  = 0.
            allocate (W_filt(-NbC:NC));      W_filt  = 0.

            allocate (hat_Ux(-NbC:NC));    hat_Ux = 0.
            allocate (hat_Uy(-NbC:NC));    hat_Uy = 0.
            allocate (hat_Uz(-NbC:NC));    hat_Uz = 0.
            allocate (hat_Vx(-NbC:NC));    hat_Vx = 0.
            allocate (hat_Vy(-NbC:NC));    hat_Vy = 0.
            allocate (hat_Vz(-NbC:NC));    hat_Vz = 0.
            allocate (hat_Wx(-NbC:NC));    hat_Wx = 0.
            allocate (hat_Wy(-NbC:NC));    hat_Wy = 0.
            allocate (hat_Wz(-NbC:NC));    hat_Wz = 0.

            allocate (Ur_filt(1:NC));     Ur_filt = 0.
            allocate (Vr_filt(1:NC));     Vr_filt = 0.
            allocate (Wr_filt(1:NC));     Wr_filt = 0.

            allocate (UUr_filt(1:NC));    UUr_filt = 0.
            allocate (VVr_filt(1:NC));    VVr_filt = 0.
            allocate (WWr_filt(1:NC));    WWr_filt = 0.
            allocate (UVr_filt(1:NC));    UVr_filt = 0.
            allocate (UWr_filt(1:NC));    UWr_filt = 0.
            allocate (VWr_filt(1:NC));    VWr_filt = 0.

            allocate (Sr11_filt(1:NC));  Sr11_filt = 0.
            allocate (Sr22_filt(1:NC));  Sr22_filt = 0.
            allocate (Sr33_filt(1:NC));  Sr33_filt = 0.
            allocate (Sr12_filt(1:NC));  Sr12_filt = 0.
            allocate (Sr13_filt(1:NC));  Sr13_filt = 0.
            allocate (Sr23_filt(1:NC));  Sr23_filt = 0.

            allocate (Srkk_filt(1:NC));  Srkk_filt = 0.
   
        end if

        allocate(ShearTest(-NbC:NC));  ShearTest = 0.0
        allocate (Ksgs(-NbC:NC));      Ksgs=0.
        allocate (Cdyn_mean(-NbC:NC)); Cdyn_mean = 0
    end if

    if(SIMULA == LES.or.SIMULA==DNS) then
        allocate (uu % mean(-NbC:NC)); uu % mean=0.
        allocate (vv % mean(-NbC:NC)); vv % mean=0.
        allocate (ww % mean(-NbC:NC)); ww % mean=0.
        allocate (uv % mean(-NbC:NC)); uv % mean=0.
        allocate (uw % mean(-NbC:NC)); uw % mean=0.
        allocate (vw % mean(-NbC:NC)); vw % mean=0.

        if(BUDG==YES) then
            allocate (uu % n(-NbC:NC)); uu % n=0.
            allocate (vv % n(-NbC:NC)); vv % n=0.
            allocate (ww % n(-NbC:NC)); ww % n=0.
            allocate (uv % n(-NbC:NC)); uv % n=0.
            allocate (uw % n(-NbC:NC)); uw % n=0.
            allocate (vw % n(-NbC:NC)); vw % n=0.

            allocate (U % fluc(-NbC:NC));  U % fluc =0.
            allocate (V % fluc(-NbC:NC));  V % fluc =0.
            allocate (W % fluc(-NbC:NC));  W % fluc =0.
            allocate (P % fluc(-NbC:NC));  P % fluc =0.

            allocate (Kx(-NbC:NC)); Kx=0.

            if(HOT==YES) then
                allocate (T % fluc(-NbC:NC));  T % fluc=0.
            end if

        end if

        allocate(VISt_mean(NC)); VISt_mean = 0.0
        allocate (ShearMean(NC));  ShearMean=0.
        if(HOT==YES) then
            allocate (T % mean(-NbC:NC));  T % mean=0.
            allocate (TT % mean(-NbC:NC)); TT % mean=0.
            allocate (uT % mean(-NbC:NC)); uT % mean=0.
            allocate (vT % mean(-NbC:NC)); vT % mean=0.
            allocate (wT % mean(-NbC:NC)); wT % mean=0.
        end if
    end if

 
    allocate (VISt(-NbC:NC)); VISt=0
    allocate (IsNearWall(NC)); IsNearWall = .FALSE.

    allocate (Vort(-NbC:NC));  Vort=0.
    allocate (Shear(-NbC:NC)); Shear=0.
    allocate (TauWall(NC));    TauWall=0.

!??????????????????????????????????????????!
!     Is there enough allocated memory     !
!??????????????????????????????????????????!
! Do something !  

END SUBROUTINE UnkAloc
