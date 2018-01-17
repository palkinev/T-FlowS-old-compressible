!======================================================================!
  SUBROUTINE LoaRes(restart)
!----------------------------------------------------------------------!
! Reads: NAME.restart                                                  !
! ~~~~~~                                                               ! 
!------------------------------[Modules]-------------------------------!
  USE all_mod
  USE pro_mod
  USE les_mod
  USE par_mod
  USE rans_mod

  USE Moin_Problem_mod
!----------------------------------------------------------------------!
  IMPLICIT NONE
!-----------------------------[Parameters]-----------------------------!
  LOGICAL   :: restart 
!-------------------------------[Locals]-------------------------------!
  INTEGER   :: c, s, m
  INTEGER   :: i_1, i_2, i_3, i_4, i_5, i_6
  CHARACTER :: nameIn*80, answer*80
  REAL      :: version
  REAL      :: r_1, r_2, r_3, r_4, r_5, r_6
!--------------------------------[CVS]---------------------------------!
!  $Id: LoaRes.f90,v 1.32 2008/12/10 14:45:17 IUS\mhadziabdic Exp $  
!  $Source: /home/IUS/mhadziabdic/.CVSROOT/T-Rex/Process/LoaRes.f90,v $  
!======================================================================!

  if(this  < 2) &              
    write(*,*) '# Input restart file name [skip cancels]:'
  call ReadC(7,inp,tn,ts,te)
  read(inp(ts(1):te(1)), '(A80)')  nameIn
  answer=nameIn
  call ToUppr(answer) 

  if(answer == 'SKIP') then
    restart = .false.
    return 
  end if

!---- save the name
  answer = name
  name = nameIn

!---------------------------!
!     Read restart file     !
!---------------------------!
  call NamFil(this, nameIn, '.restart', len_trim('.restart') )
  open(9, FILE=nameIn, FORM='UNFORMATTED')
  write(6, *) '# Now reading the file:', nameIn

!---- version
  read(9) version ! version

!---- 60 INTEGER parameters ----------------------------------------
  read(9)      i_1,      NbC,       NC,       NS,     Ndtt,    Nstat
  read(9)       Cm,      i_2,      i_3,      i_4,      i_5,      ini
  read(9)    ALGOR,    INERT,   CONVEC,    CROSS,   DIFFUS,   SIMULA
  read(9)   POSPRO,  CHANNEL,     TEST,    OTHER,      HOT,      i_6
  read(9)    BLEND,      i_2,      i_3,     MODE,     PIPE,      i_6
  read(9)      i_1,      i_2,      i_3,      i_4,      i_5,      i_6
  read(9)      i_1,      i_2,      i_3,      i_4,      i_5,      i_6
  read(9)      i_1,      i_2,      i_3,      i_4,      i_5,      i_6
  read(9)      i_1,      i_2,      i_3,      i_4,      i_5,      i_6
  read(9)      i_1,      i_2,      i_3,      i_4,      i_5,      i_6
!-------------------------------------------------------------------

!---- 60 REAL parameters --------------------------------------
  read(9)     r_1,        r_2,        r_3,     xp,      yp,     zp  
  read(9)     r_1,        r_2,        r_3,    r_4,     r_4,    r_6             
  read(9)   ReTau,        r_2,        Cs0,    r_4,     r_4,    r_6 
  read(9)      dt,       Time,      Kflow,    r_4,     r_5,    r_6 
  read(9)   U%URF,      P%URF,       URFC, SIMTol,  U%Stol,    r_6 
  read(9) PP%Stol,   Zmix%URF,  Zmix%STol,    r_4,     r_5,    r_6   
  read(9)     r_1,        r_2,        r_3,    r_4,     r_5,    r_6
  read(9)     r_1,        r_2,        r_3,    r_4,     r_5,    r_6 
  read(9)     r_1,        r_2,        r_3,    r_4,     r_5,    r_6   
  read(9)     r_1,        r_2,        r_3,    r_4,     r_5,    r_6   
!--------------------------------------------------------------

  call UnkAloc

!----density
  read(9) (Rho % n(c),   c=-NbC,NC)
  read(9) (Rho % o(c),   c=-NbC,NC)
  read(9) (Rho % oo(c),  c=-NbC,NC)
!-----U
  read(9) (U  % n(c),   c=-NbC,NC)
!  read(9) (Ur % n(c),   c=-NbC,NC)
!  read(9) (Ur % o(c),   c=1,NC)
!  read(9) (Ur %oo(c),   c=1,NC)

  !-----V
  read(9) (V  % n(c),   c=-NbC,NC)
!  read(9) (Vr % n(c),   c=-NbC,NC)
!  read(9) (Vr % o(c),   c=1,NC)
!  read(9) (Vr %oo(c),   c=1,NC)

!-----W
  read(9) (W  % n(c),   c=-NbC,NC)
!  read(9) (Wr % n(c),   c=-NbC,NC)
!  read(9) (Wr % o(c),   c=1,NC)
!  read(9) (Wr %oo(c),   c=1,NC)

  read(9) (P % n(c),   c=-NbC,NC)
  read(9) (PP % n(c),  c=-NbC,NC)

  read(9) (Px(c),   c=-NbC,NC)
  read(9) (Py(c),   c=-NbC,NC)
  read(9) (Pz(c),   c=-NbC,NC)


  read(9) (VISt(c), c=-NbC,NC)

!---- Fluxes 
  read(9) (Flux(s),   s=1,NS)
  read(9) (Flux_u(s), s=1,NS)

!---- LES and DNS
  if(SIMULA == LES) then
    read(9) (U % mean(c),  c=-NbC,NC)
    read(9) (V % mean(c),  c=-NbC,NC)
    read(9) (W % mean(c),  c=-NbC,NC)
    read(9) (uu % mean(c), c=-NbC,NC)
    read(9) (vv % mean(c), c=-NbC,NC)
    read(9) (ww % mean(c), c=-NbC,NC)
    read(9) (uv % mean(c), c=-NbC,NC)
    read(9) (uw % mean(c), c=-NbC,NC)
    read(9) (vw % mean(c), c=-NbC,NC)

    read(9) (P % mean(c),  c=1,NC)
    read(9) (VISt_mean(c), c=1,NC)
    read(9) (near(c),   c=-NbC,NC)
    if(MODE == DYN) read(9) (Cdyn_mean(c), c=1,NC)
  end if

!---- Temperature


  if(HOT == YES) then
    read(9) (Zmix  % n(c),   c=-NbC,NC)
    !read(9) (Zmixr % q(c),   c=-NbC,-1)
!    read(9) (Zmixr % n(c),   c=-NbC,NC)
!    read(9) (Zmixr % o(c),   c=1,NC)
!    read(9) (Zmixr %oo(c),   c=1,NC)

    !only for moin problem 10
    if (Moin_Problem().eq.10) then
      read(9) (VISc_Dyn(c),    c=1,NC)
      read(9) (VIScZmix_Dyn(c),c=1,NC)
    end if 

!    if(SIMULA == LES .or. SIMULA == DNS) then
!      read(9) (Zmixr % mean(c), c=-NbC,NC)
!      read(9) (TT % mean(c), c=-NbC,NC)
!      read(9) (uT % mean(c), c=-NbC,NC)
!      read(9) (vT % mean(c), c=-NbC,NC)
!      read(9) (wT % mean(c), c=-NbC,NC)
!    end if
  end if

!---- Pressure drops in each material (domain)
  do m=1,Nmat
    read(9) PdropX(m), PdropY(m), PdropZ(m)
    read(9) FLUXoX(m), FLUXoY(m), FLUXoZ(m)
    read(9) FLUXx(m),  FLUXy(m),  FLUXz(m)
    read(9) AreaX(m),  AreaY(m),  AreaZ(m)
    read(9) Ubulk(m),  Vbulk(m),  Wbulk(m)
  end do

  close(9)

  restart = .true.

!---- restore the name
  name = answer 

  END SUBROUTINE LoaRes 
