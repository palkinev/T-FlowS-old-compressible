!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
!                                 !                                    !
!                                 !   Bojan Niceno                     !
!   Global variable definitions   !   Delft University of Technology   !
!        for the processor        !   Section Heat Transfer            !
!                                 !   niceno@duttwta.wt.tn.tudelft.nl  !
!                                 !                                    !
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
!..RCS/CVS ident
! $Id: pro_mod.h90,v 1.43 2005/01/25 12:38:34 muhamed Exp $
! $Source: /home/muhamed/.CVSROOT/T-Rex/Library/pro_mod.h90,v $
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
MODULE pro_mod

    USE allp_mod

    IMPLICIT NONE

    !----- Right hand side for velocity and pressure equations 
    REAL,ALLOCATABLE, PUBLIC :: b(:)

    !----- System matrix for velocity and pressure correction
    REAL,ALLOCATABLE, PUBLIC    :: Aval(:)
    REAL,ALLOCATABLE, PUBLIC    :: Asave(:)
    INTEGER,ALLOCATABLE, PUBLIC :: Arow(:),  Acol(:),  Adia(:)  
    INTEGER,ALLOCATABLE, PUBLIC :: SidAij(:,:)     

    !----- Used in Dynamic Smgaorinsky model ----------------------------!
    REAL,ALLOCATABLE, PUBLIC    :: Aval_dif(:)
    !--------------------------------------------------------------------!

    !----- Parts of the matrix for boundary conditions. 
    REAL,ALLOCATABLE, PUBLIC :: Abou(:)  

    !----- Correlation points
    REAL, PUBLIC :: R11_1, R11_2, R11_3, R11_4, R11_5
    REAL, PUBLIC :: R11_6, R11_7, R11_8, R11_9, R11_10
    REAL, PUBLIC :: A11_1, A11_2, A11_3, A11_4, A11_5
    REAL, PUBLIC :: A11_6, A11_7, A11_8, A11_9, A11_10

    !----- Diffusivity coefficients 
    REAL,ALLOCATABLE, PUBLIC :: DIFz(:)                                        ! NEW

    !----- Velocity derivativeses 
    REAL,ALLOCATABLE, PUBLIC :: Ux(:), Uy(:), Uz(:)
    REAL,ALLOCATABLE, PUBLIC :: Vx(:), Vy(:), Vz(:)
    REAL,ALLOCATABLE, PUBLIC :: Wx(:), Wy(:), Wz(:)

    !----- Zmix derivativeses
    REAL,ALLOCATABLE, PUBLIC :: Zmixx(:), Zmixy(:), Zmixz(:)

    !----- Pressure derivativeses dP/dx .... 
    REAL,ALLOCATABLE, PUBLIC :: Px(:), Py(:), Pz(:)

    REAL,ALLOCATABLE, PUBLIC :: Kx(:)

    !----- Pressure at the cell faces  
    REAL,ALLOCATABLE, PUBLIC :: Ps(:)

!    REAL,ALLOCATABLE, PUBLIC :: VAR1x(:),   VAR1y(:),   VAR1z(:)
!    REAL,ALLOCATABLE, PUBLIC :: VAR2x(:),   VAR2y(:),   VAR2z(:)
!    REAL,ALLOCATABLE, PUBLIC :: VAR3x(:),   VAR3y(:),   VAR3z(:)
!    REAL,ALLOCATABLE, PUBLIC :: PHI1x(:),   PHI1y(:),   PHI1z(:)
!    REAL,ALLOCATABLE, PUBLIC :: PHI2x(:),   PHI2y(:),   PHI2z(:)
!    REAL,ALLOCATABLE, PUBLIC :: PHI3x(:),   PHI3y(:),   PHI3z(:)
!    REAL,ALLOCATABLE, PUBLIC :: PHI4x(:),   PHI4y(:),   PHI4z(:)
!    REAL,ALLOCATABLE, PUBLIC :: PHI5x(:),   PHI5y(:),   PHI5z(:)
!    REAL,ALLOCATABLE, PUBLIC :: PHI6x(:),   PHI6y(:),   PHI6z(:)
!    REAL,ALLOCATABLE, PUBLIC :: PHI7x(:),   PHI7y(:),   PHI7z(:)
!    REAL,ALLOCATABLE, PUBLIC :: PHI8x(:),   PHI8y(:),   PHI8z(:)
!    REAL,ALLOCATABLE, PUBLIC :: PHI9x(:),   PHI9y(:),   PHI9z(:)
!    REAL,ALLOCATABLE, PUBLIC :: PHI10x(:),  PHI10y(:),  PHI10z(:)
!    REAL,ALLOCATABLE, PUBLIC :: PHI11x(:),  PHI11y(:),  PHI11z(:)

    !----- For convective schemes
    REAL,ALLOCATABLE, PUBLIC :: PHImax(:), PHImin(:) 

    !----- Scalar
    REAL,ALLOCATABLE, PUBLIC :: PHIx(:),   PHIy(:),   PHIz(:)
    REAL,ALLOCATABLE, PUBLIC :: PHIrx(:),  PHIry(:),  PHIrz(:)
    REAL,ALLOCATABLE, PUBLIC :: PHIside(:)


    !----- Velocity components
    TYPE(Unknown), PUBLIC :: U
    TYPE(Unknown), PUBLIC :: V
    TYPE(Unknown), PUBLIC :: W
    REAL,POINTER, PUBLIC  :: VISc_Dyn(:) ! in case of heterogeneous viscosity profile

    !----- "Mass flux" components = rho*u_i
    !TYPE(Unknown), PUBLIC :: Ur
    !TYPE(Unknown), PUBLIC :: Vr
    !TYPE(Unknown), PUBLIC :: Wr

    !----- Temperature
    TYPE(Unknown), PUBLIC :: T

    !----- Density
    TYPE(Unknown), PUBLIC :: Rho

    !----- Z scalar
    TYPE(Unknown), PUBLIC :: Zmix
    !TYPE(Unknown), PUBLIC :: Zmixr
    REAL,POINTER, PUBLIC  :: VIScZmix_Dyn(:) ! in case of heterogeneous Diffusivity profile
    REAL, PUBLIC          :: Z_SRC_INT

    !----- Pressure
    TYPE(Unknown), PUBLIC :: P
    TYPE(Unknown), PUBLIC :: PP

    !----- Turbulence models variables
    TYPE(Unknown), PUBLIC :: KIN
    TYPE(Unknown), PUBLIC :: EPS
    TYPE(Unknown), PUBLIC :: V_2
    TYPE(Unknown), PUBLIC :: F22
    TYPE(Unknown), PUBLIC :: VIS

    !----- Stresses
    TYPE(Unknown), PUBLIC :: uu
    TYPE(Unknown), PUBLIC :: vv
    TYPE(Unknown), PUBLIC :: ww
    TYPE(Unknown), PUBLIC :: uv
    TYPE(Unknown), PUBLIC :: uw
    TYPE(Unknown), PUBLIC :: vw

    TYPE(Unknown), PUBLIC :: TT, uT, vT, wT

    TYPE(Unknown), PUBLIC :: uuu, uuv, uuw
    TYPE(Unknown), PUBLIC :: vvu, vvv, vvw
    TYPE(Unknown), PUBLIC :: wwu, wwv, www
    TYPE(Unknown), PUBLIC :: uwu, uwv, uww


    !=====================================================================!
    !        Hybrid apriori variables
    !=====================================================================!

    !----- Turbulent viscosity
    REAL,ALLOCATABLE, PUBLIC :: VISt_sgs(:)
    REAL,ALLOCATABLE, PUBLIC :: VISt_eff(:)

    !----- Mass fluxes throught cell faces
    REAL,ALLOCATABLE, PUBLIC :: Flux_r(:), Alfa_lim(:)

    !----- Mass fluxes throught the whole domain
    REAL,ALLOCATABLE, PUBLIC :: FLUXx_r(:),  FLUXy_r(:),  FLUXz_r(:)
    !=====================================================================!

    !------------------------------!
    !     Algorythm parameters     !
    !------------------------------!
!   INTEGER, PUBLIC :: K_EPS
!   INTEGER, PUBLIC :: K_EPS_VV
    INTEGER, PUBLIC :: HRe
    INTEGER, PUBLIC :: MODE   
    INTEGER, PUBLIC :: LRe
    !INTEGER, PUBLIC :: SPA_ALL
    !INTEGER, PUBLIC :: DES_SPA
    INTEGER, PUBLIC :: J_L    
    INTEGER, PUBLIC :: NAG     
    INTEGER, PUBLIC :: S_L_Y   
    INTEGER, PUBLIC :: WOLF   
!   INTEGER, PUBLIC :: ZETA
!   INTEGER, PUBLIC :: HYB_ZETA
!   INTEGER, PUBLIC :: HYB_PITM
    INTEGER, PUBLIC :: RNG   
    INTEGER, PUBLIC :: SMAG
    INTEGER, PUBLIC :: DYN 
!    INTEGER, PUBLIC :: WALE
    INTEGER, PUBLIC :: MIX  
!    INTEGER, PUBLIC :: ZPANS
    !INTEGER, PUBLIC :: ZETAM
    !INTEGER, PUBLIC :: EBM
    INTEGER, PUBLIC :: HYB
    !INTEGER, PUBLIC :: HJ
    INTEGER, PUBLIC :: WF
    INTEGER, PUBLIC :: STAN

    !----- Mass fluxes throught cell faces (rho*U*S; U*S)
    REAL,ALLOCATABLE, PUBLIC :: Flux(:)                               
    REAL,ALLOCATABLE, PUBLIC :: Flux_u(:)                                      ! Flux_u - NEW

    !---- Geometrical staff 
    REAL,ALLOCATABLE, PUBLIC :: Scoef(:)
    REAL,ALLOCATABLE, PUBLIC :: G(:,:) 
    REAL,ALLOCATABLE, PUBLIC :: fF(:)   ! weight factors for the fluid phase

    LOGICAL,ALLOCATABLE, PUBLIC :: IsNearWall(:)
    LOGICAL,ALLOCATABLE, PUBLIC :: IsNearPeri(:)
!    LOGICAL,ALLOCATABLE, PUBLIC :: IsNearWall_2(:)
    !LOGICAL,ALLOCATABLE, PUBLIC :: IsNearWall_3(:)
    LOGICAL,ALLOCATABLE, PUBLIC :: IsNearInflow(:)
    LOGICAL,ALLOCATABLE, PUBLIC :: ConvZone1(:)

    !---- Cells which are bad for calculation of gradients
    LOGICAL,ALLOCATABLE, PUBLIC :: BadForG(:)
    INTEGER,ALLOCATABLE, PUBLIC :: NumGood(:), NumNeig(:)

    !----- Mass fluxes throught the whole domain
    REAL,ALLOCATABLE, PUBLIC :: MassIn(:), MasOut(:) 
    REAL,ALLOCATABLE, PUBLIC :: FLUXx(:),  FLUXy(:),  FLUXz(:)
    REAL,ALLOCATABLE, PUBLIC :: FLUXoX(:), FLUXoY(:), FLUXoZ(:) 
    REAL,ALLOCATABLE, PUBLIC :: Ubulk(:),  Vbulk(:),  Wbulk(:)

    !----- Viscosity, Density, Conductivity
    INTEGER, PUBLIC :: StateMat(100)
    INTEGER, PUBLIC :: SimulMat(100)
    REAL, PUBLIC    :: VISc, VIScZmix, DENc(100), CONc(100), CAPc(100)

    !---- angular velocity 
    REAL, PUBLIC    :: omegaX, omegaY, omegaZ, omega

    !---- Time step and total time
    REAL, PUBLIC    :: dt, Time

    !----- Constants needed for UserProbe2d (cut lines)
    REAL, PUBLIC      :: x_o, y_o, Href
    INTEGER, PUBLIC   :: Ncuts
    CHARACTER, PUBLIC :: namCut*80

    !---- Integer variable needed for interpolation of
    !---- results between different meshes tranfer (LoaIni)
    INTEGER, PUBLIC          :: NClast, N_sign, eqn
    INTEGER,ALLOCATABLE, PUBLIC :: near(:)
    INTEGER,ALLOCATABLE, PUBLIC :: near_2(:)
    INTEGER,ALLOCATABLE, PUBLIC :: near_3(:)
    !INTEGER,ALLOCATABLE, PUBLIC :: connect(:)
    INTEGER,ALLOCATABLE, PUBLIC :: connect2(:)

    !----- Residuals
    REAL, PUBLIC    :: errmax, res(100)  

    !----- Monitoring planes for each material (domain)
    REAL,ALLOCATABLE, PUBLIC :: xp(:), yp(:), zp(:)

    !---------------------------!
    !     Solver parameters     !
    !---------------------------!
    REAL, PUBLIC    :: URFC(100), SIMTol, URFC_Tur(100), URFC_Tem(100)
    REAL, PUBLIC    :: TRFC(100)

    !----- Under-relaxation parameter for turbulent quantity
    REAL, PUBLIC    :: URFT, Alfa_fin1, Alfa_fin2

    !-----------------------------------!
    !     Area of the cross section     !
    !-----------------------------------!
    REAL,ALLOCATABLE, PUBLIC :: AreaX(:), AreaY(:), AreaZ(:)           
    REAL, PUBLIC :: Area, Tflux, Qflux, Xmax, Ymax, Zmax           

    !------------------------------!
    !     Algorythm parameters     !
    !------------------------------!
    INTEGER, PUBLIC :: INERT,    CONVEC,    CROSS,    DIFFUS 
    INTEGER, PUBLIC :: LIN,      PAR,       AB,       CN,       FI
    INTEGER, PUBLIC :: ALGOR,    SIMPLE,    FRACT
    INTEGER, PUBLIC :: SIMULA,   DNS,       LES 
    INTEGER, PUBLIC :: POSPRO,   AVS,       GMV
    INTEGER, PUBLIC :: CHANNEL,  TEST,      OTHER,    HOT, HOTini, PIPE, JET, ROT, BUDG, COHERENT
    INTEGER, PUBLIC :: SHAKE(100),    BLEND(100),BLEND_TUR(100), BLEND_TEM(100), YES, NO
    INTEGER, PUBLIC :: SHAKE_PER(100),SHAKE_INT(100)
    INTEGER, PUBLIC :: PREC 
    INTEGER, PUBLIC :: CDS,      QUICK,    LUDS,     MINMOD,   SMART,    AVL_SMART, SUPERBEE, GAMMA
    INTEGER, PUBLIC :: XHOM,     YHOM,     ZHOM

    INTEGER, PUBLIC ::  MoinID

    INTEGER,PARAMETER, PUBLIC :: MAXM=100 
    INTEGER, PUBLIC   :: Cm(MAXM), Nmon
    REAL, PUBLIC      :: NOM(MAXM), DEN(MAXM), R11(MAXM), U_f(MAXM)

    INTEGER, PUBLIC :: Ndt, Ndtt, Nstat, Nini, ini, Ndyn, Nstat2, NewSta, NK, Nbudg 

    INTEGER, PUBLIC :: NONZERO

    !---------------------------------------------------------------------------!
    ! LineMon:   1:  6 -> Time step 
    ! ~~~~~~~~   7: 18 -> Time                
    !           19: 66 -> U,V,W,P monitoring
    !           67: 78 -> T  monitoring 
    !           79: 90 -> FLUXx
    !           91:102 -> P drop
    !          103:114 -> CFL
    !          115:126 -> Pe
    !          127:138 -> Kin.en.           
    !---------------------------------------------------------------------------!
    CHARACTER(138), PUBLIC :: LinMon0 ! everything that goes on the screen
    CHARACTER(138), PUBLIC :: LinMon1 ! everything that goes on the screen
    CHARACTER(138), PUBLIC :: LinMon2 ! everything that goes on the screen
    CHARACTER(138), PUBLIC :: LinMon3 ! everything that goes on the screen
    CHARACTER(138), PUBLIC :: LinMon4 ! everything that goes on the screen
    CHARACTER(138), PUBLIC :: LinMon5 ! everything that goes on the screen
    CHARACTER(138), PUBLIC :: LinMon6 ! everything that goes on the screen
    CHARACTER(138), PUBLIC :: LinMon7 ! everything that goes on the screen
    CHARACTER(138), PUBLIC :: LinMon8 ! everything that goes on the screen
    CHARACTER(138), PUBLIC :: LinMon9 ! everything that goes on the screen
    !---------------------------------------------------------------------------!
    ! LineRes:   1:  1 -> #
    ! ~~~~~~~~   2:  4 -> ini
    !            5: 16 -> errmax 
    !           17: 28 -> res U
    !           29: 40 -> res V
    !           41: 52 -> res W
    !           53: 64 -> res PP
    !           65: 76 -> res T
    !           77: 80 -> iter U
    !           81: 84 -> iter V
    !           85: 88 -> iter W
    !           89: 92 -> iter P
    !           93: 96 -> iter T
    !---------------------------------------------------------------------------!
    CHARACTER(100), PUBLIC :: LineRes              ! everything that goes on the screen
    CHARACTER(100), PUBLIC :: LineRes1             ! everything that goes on the screen
    CHARACTER, PUBLIC      :: namIni(128)*80
END MODULE
