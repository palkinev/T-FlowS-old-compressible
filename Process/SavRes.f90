!======================================================================!
  SUBROUTINE SavRes(namAut)
!----------------------------------------------------------------------!
! Writes: NAME.restart                                                 !
! ~~~~~~~                                                              !
!------------------------------[Modules]-------------------------------!
  USE all_mod
  USE pro_mod
  USE les_mod
  USE par_mod
  USE rans_mod

  USE Moin_Problem_mod 
!----------------------------------------------------------------------!
  IMPLICIT NONE
!-------------------------------[Locals]-------------------------------!
  INTEGER             :: c, s, m
  CHARACTER           :: namOut*80, answer*80
  CHARACTER, OPTIONAL :: namAut*(*)
!--------------------------------[CVS]---------------------------------!
!  $Id: SavRes.f90,v 1.33 2008/12/10 14:58:16 IUS\mhadziabdic Exp $  
!  $Source: /home/IUS/mhadziabdic/.CVSROOT/T-Rex/Process/SavRes.f90,v $  
!======================================================================!

  if(PRESENT(namAut)) then
!---- save the name
    answer = name
    name = namAut
  else
    if(this  < 2)                                                     &
      write(*,*) '# Output restart file name [skip cancels]:'
    call ReadC(7,inp,tn,ts,te)
!->>> write(*,*) inp(1:300)
    read(inp(ts(1):te(1)), '(A80)')  namOut
    answer=namOut
    call ToUppr(answer) 

    if(answer == 'SKIP') return 

!---- save the name
    answer = name
    name = namOut
  end if

!-----------------------------!
!     Create restart file     !
!-----------------------------!
  call NamFil(this, namOut, '.restart', len_trim('.restart') )
  open(9, FILE=namOut, FORM='UNFORMATTED')
  if(this  < 2) write(6, *) '# Now creating the file:', namOut

!---- version
  write(9) 0.0  ! version

!---- 60 INTEGER parameters -----------------------------------------
  write(9)        0,      NbC,       NC,       NS,     Ndtt,    Nstat
  write(9)       Cm,        0,        0,        0,        0,    ini+1
  write(9)    ALGOR,    INERT,   CONVEC,    CROSS,   DIFFUS,   SIMULA
  write(9)   POSPRO,  CHANNEL,     TEST,    OTHER,      HOT,        0
  write(9)    BLEND,        0,        0,     MODE,     PIPE,        0
  write(9)        0,        0,        0,        0,        0,        0
  write(9)        0,        0,        0,        0,        0,        0
  write(9)        0,        0,        0,        0,        0,        0
  write(9)        0,        0,        0,        0,        0,        0
  write(9)        0,        0,        0,        0,        0,        0
!--------------------------------------------------------------------

!---- 60 REAL parameters --------------------------------------
  write(9)     0.0,     0.0,        0.0,     xp,     yp,     zp  
  write(9)     0.0,     0.0,        0.0,    0.0,     0.0,    0.0     
  write(9)   ReTau,     0.0,        Cs0,    0.0,     0.0,    0.0
  write(9)      dt,     Time,     Kflow,    0.0,     0.0,    0.0
  write(9)   U%URF,    P%URF,      URFC, SIMTol,  U%Stol,    0.0 
  write(9) PP%Stol, Zmix%URF,Zmix %STol,    0.0,     0.0,    0.0
  write(9)     0.0,      0.0,       0.0,    0.0,     0.0,    0.0
  write(9)     0.0,      0.0,       0.0,    0.0,     0.0,    0.0
  write(9)     0.0,      0.0,       0.0,    0.0,     0.0,    0.0
  write(9)     0.0,      0.0,       0.0,    0.0,     0.0,    0.0
!--------------------------------------------------------------

!----density
  write(9) (Rho % n(c),   c=-NbC,NC)
  write(9) (Rho % o(c),   c=-NbC,NC)
  write(9) (Rho % oo(c),  c=-NbC,NC)
!-----U
  write(9) (U  % n(c),   c=-NbC,NC)
!  write(9) (Ur % n(c),   c=-NbC,NC)
!  write(9) (Ur % o(c),   c=1,NC)
!  write(9) (Ur %oo(c),   c=1,NC)
!  write(9) (Ur % C(c),   c=1,NC)
!  write(9) (Ur % Co(c),  c=1,NC)
!  write(9) (Ur %Coo(c),  c=1,NC)
!  write(9) (Ur % Do(c),  c=1,NC)
!  write(9) (Ur %Doo(c),  c=1,NC)
!  write(9) (Ur % X(c),   c=1,NC)
!  write(9) (Ur % Xo(c),  c=1,NC)
!  write(9) (Ur %Xoo(c),  c=1,NC)
!-----V
  write(9) (V  % n(c),   c=-NbC,NC)
!  write(9) (Vr % n(c),   c=-NbC,NC)
!  write(9) (Vr % o(c),   c=1,NC)
!  write(9) (Vr %oo(c),   c=1,NC)
!  write(9) (Vr % C(c),   c=1,NC)
!  write(9) (Vr % Co(c),  c=1,NC)
!  write(9) (Vr %Coo(c),  c=1,NC)
!  write(9) (Vr % Do(c),  c=1,NC)
!  write(9) (Vr %Doo(c),  c=1,NC)
!  write(9) (Vr % X(c),   c=1,NC)
!  write(9) (Vr % Xo(c),  c=1,NC)
!  write(9) (Vr %Xoo(c),  c=1,NC)
!-----W
  write(9) (W  % n(c),   c=-NbC,NC)
!  write(9) (Wr % n(c),   c=-NbC,NC)
!  write(9) (Wr % o(c),   c=1,NC)
!  write(9) (Wr %oo(c),   c=1,NC)
!  write(9) (Wr % C(c),   c=1,NC)
!  write(9) (Wr % Co(c),  c=1,NC)
!  write(9) (Wr %Coo(c),  c=1,NC)
!  write(9) (Wr % Do(c),  c=1,NC)
!  write(9) (Wr %Doo(c),  c=1,NC)
!  write(9) (Wr % X(c),   c=1,NC)
!  write(9) (Wr % Xo(c),  c=1,NC)
!  write(9) (Wr %Xoo(c),  c=1,NC)

  write(9) (P % n(c),   c=-NbC,NC)
  write(9) (PP % n(c),  c=-NbC,NC)

  write(9) (Px(c),   c=-NbC,NC)
  write(9) (Py(c),   c=-NbC,NC)
  write(9) (Pz(c),   c=-NbC,NC)


  write(9) (VISt(c), c=-NbC,NC)

!---- Fluxes 
  write(9) (Flux(s),   s=1,NS)
  write(9) (Flux_u(s), s=1,NS)

!---- LES and DNS
  if(SIMULA == LES) then
    write(9) (U % mean(c),  c=-NbC,NC)
    write(9) (V % mean(c),  c=-NbC,NC)
    write(9) (W % mean(c),  c=-NbC,NC)
    write(9) (uu % mean(c), c=-NbC,NC)
    write(9) (vv % mean(c), c=-NbC,NC)
    write(9) (ww % mean(c), c=-NbC,NC)
    write(9) (uv % mean(c), c=-NbC,NC)
    write(9) (uw % mean(c), c=-NbC,NC)
    write(9) (vw % mean(c), c=-NbC,NC)

    write(9) (P % mean(c),  c=1,NC)
    write(9) (VISt_mean(c), c=1,NC)
    write(9) (near(c),   c=-NbC,NC)
    if(MODE == DYN) write(9) (Cdyn_mean(c), c=1,NC)
  end if

!---- Temperature


  if(HOT == YES) then
    write(9) (Zmix  % n(c),   c=-NbC,NC)
    !write(9) (Zmixr % q(c),   c=-NbC,-1)
!    write(9) (Zmixr % n(c),   c=-NbC,NC)
!    write(9) (Zmixr % o(c),   c=1,NC)
!    write(9) (Zmixr %oo(c),   c=1,NC)
!    write(9) (Zmixr % C(c),   c=1,NC)
!    write(9) (Zmixr % Co(c),  c=1,NC)
!    write(9) (Zmixr %Coo(c),  c=1,NC)
!    write(9) (Zmixr % Do(c),  c=1,NC)
!    write(9) (Zmixr %Doo(c),  c=1,NC)
!    write(9) (Zmixr % X(c),   c=1,NC)
!    write(9) (Zmixr % Xo(c),  c=1,NC)
!    write(9) (Zmixr %Xoo(c),  c=1,NC)

    !only for moin problem 10
    if (Moin_Problem().eq.10) then
      write(9) (VISc_Dyn(c),    c=1,NC)
      write(9) (VIScZmix_Dyn(c),c=1,NC)
    end if 

!    if(SIMULA == LES .or. SIMULA == DNS) then
!      write(9) (Zmixr % mean(c), c=-NbC,NC)
!      write(9) (TT % mean(c), c=-NbC,NC)
!      write(9) (uT % mean(c), c=-NbC,NC)
!      write(9) (vT % mean(c), c=-NbC,NC)
!      write(9) (wT % mean(c), c=-NbC,NC)
!    end if
  end if

!---- Pressure drops in each material (domain)
  do m=1,Nmat
    write(9) PdropX(m), PdropY(m), PdropZ(m)
    write(9) FLUXoX(m), FLUXoY(m), FLUXoZ(m)
    write(9) FLUXx(m),  FLUXy(m),  FLUXz(m)
    write(9) AreaX(m),  AreaY(m),  AreaZ(m)
    write(9) Ubulk(m),  Vbulk(m),  Wbulk(m)
  end do

  close(9)


!---- restore the name
  name = answer 

  END SUBROUTINE SavRes
