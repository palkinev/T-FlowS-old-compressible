!======================================================================!
SUBROUTINE LoaIni()
    !----------------------------------------------------------------------!
    ! This version of LoaIni is optimised for very large meshes
    ! Program SUB_INI needs to be used to create files needed by this
    ! subroutine
    !------------------------------[Modules]-------------------------------!
    USE allp_mod, only : huge
    USE all_mod
    USE pro_mod
    USE les_mod
    USE par_mod
    USE rans_mod
    !----------------------------------------------------------------------!
    IMPLICIT NONE
    !-------------------------------[Locals]-------------------------------!
    INTEGER          :: j, k,  c, nearest, c1, c2, s
    INTEGER          :: NCold
    !======================================================================*

    REAL,ALLOCATABLE :: Xold(:),Yold(:),Zold(:)

    REAL,ALLOCATABLE :: Uold(:),Vold(:),Wold(:)
    REAL,ALLOCATABLE :: UCold(:),VCold(:),WCold(:)
    REAL,ALLOCATABLE :: UCoold(:),VCoold(:),WCoold(:)
    REAL,ALLOCATABLE :: Uoold(:),Voold(:),Woold(:)
    REAL,ALLOCATABLE :: UDoold(:),VDoold(:),WDoold(:)
    REAL,ALLOCATABLE :: UXold(:),VXold(:),WXold(:)
    REAL,ALLOCATABLE :: UXoold(:),VXoold(:),WXoold(:)
    !stresses----------------------------------------
    !uu--------------------------------------------
    REAL,ALLOCATABLE :: uuold(:),uuCold(:),uuCoold(:),uuoold(:),uuDoold(:),uuXold(:),uuXoold(:)
    REAL,ALLOCATABLE :: vvold(:),vvCold(:),vvCoold(:),vvoold(:),vvDoold(:),vvXold(:),vvXoold(:)
    REAL,ALLOCATABLE :: wwold(:),wwCold(:),wwCoold(:),wwoold(:),wwDoold(:),wwXold(:),wwXoold(:)
    REAL,ALLOCATABLE :: uvold(:),uvCold(:),uvCoold(:),uvoold(:),uvDoold(:),uvXold(:),uvXoold(:)
    REAL,ALLOCATABLE :: uwold(:),uwCold(:),uwCoold(:),uwoold(:),uwDoold(:),uwXold(:),uwXoold(:)
    REAL,ALLOCATABLE :: vwold(:),vwCold(:),vwCoold(:),vwoold(:),vwDoold(:),vwXold(:),vwXoold(:)
    REAL,ALLOCATABLE :: f22old(:),f22oold(:),f22Doold(:),f22Xold(:),f22Xoold(:)
    REAL,ALLOCATABLE :: Epsold(:),EpsCold(:),EpsCoold(:),Epsoold(:),EpsDoold(:),EpsXold(:),EpsXoold(:)


    REAL,ALLOCATABLE :: Pold(:)
    REAL,ALLOCATABLE :: PPold(:)
    REAL,ALLOCATABLE :: Pxold(:),Pyold(:),Pzold(:)

    REAL          :: DISTold
    REAL          :: Dist, Us, Ws, Vs

    !---- Variables for ReadC:
    CHARACTER     :: answer*80
    CHARACTER(80) :: namAut
    !--------------------------------[CVS]---------------------------------!
    character(80) :: rcs1,rcs2
    data rcs1/                                                        &
        '$Id: LoaIni.f90,v 1.4 2008/12/10 14:44:15 IUS\mhadziabdic Exp $'/
    data rcs2/                                                        &
        '$Source: /home/IUS/mhadziabdic/.CVSROOT/T-Rex/Process/LoaIni.f90,v $'/
    !======================================================================!

    call ReadC(7,inp,tn,ts,te)
    !->>> write(*,*) inp(1:300)
    read(inp(ts(1):te(1)), '(A80)') namAut
    answer=namAut
    call ToUppr(answer)
    if(answer == 'SKIP') return
    !---- save the name
    answer = name
    name = namAut


    !  if(tn==2) read(inp(ts(2):te(2)),'(A8)') answer
    !  call ToUppr(answer)
    !  if(answer == 'HOT') then
    !    HOTini = YES
    !  end if

    call NamFil(THIS, namAut, '.ini', len_trim('.ini'))

    if(this < 2) write(*,*)'now reading file:', namAut

    open(5, FILE=namAut)
    read(5,*) NCold

    write(*,*) 'Old cells count =', NCold
    !xyz--------------------------
    allocate (Xold(NCold)); Xold = 0.0
    allocate (Yold(NCold)); Yold = 0.0
    allocate (Zold(NCold)); Zold = 0.0

    !UVW--------------------------
    allocate (Uold(NCold)); Uold = 0.0
    allocate (Vold(NCold)); Vold = 0.0
    allocate (Wold(NCold)); Wold = 0.0

    allocate (Uoold(NCold)); Uoold = 0.0
    allocate (Voold(NCold)); Voold = 0.0
    allocate (Woold(NCold)); Woold = 0.0

    allocate (UCold(NCold)); UCold = 0.0
    allocate (VCold(NCold)); VCold = 0.0
    allocate (WCold(NCold)); WCold = 0.0

    allocate (UCoold(NCold)); UCoold = 0.0
    allocate (VCoold(NCold)); VCoold = 0.0
    allocate (WCoold(NCold)); WCoold = 0.0

    allocate (UDoold(NCold)); UDoold = 0.0
    allocate (VDoold(NCold)); VDoold = 0.0
    allocate (WDoold(NCold)); WDoold = 0.0

    allocate (UXold(NCold)); UXold = 0.0
    allocate (VXold(NCold)); VXold = 0.0
    allocate (WXold(NCold)); WXold = 0.0

    allocate (UXoold(NCold)); UXoold = 0.0
    allocate (VXoold(NCold)); VXoold = 0.0
    allocate (WXoold(NCold)); WXoold = 0.0

    !P-----------------------------
    allocate (Pold(NCold));  Pold  = 0.0
    allocate (PPold(NCold)); PPold = 0.0
    allocate (Pxold(NCold)); Pxold = 0.0
    allocate (Pyold(NCold)); Pyold = 0.0
    allocate (Pzold(NCold)); Pzold = 0.0

    !stresses---------------------------
    !uu-------------------------------
    allocate (uuold(NCold));   uuold   = 0.0
    allocate (uuoold(NCold));  uuoold  = 0.0
    allocate (uuCold(NCold));  uuCold  = 0.0
    allocate (uuCoold(NCold)); uuCoold = 0.0
    allocate (uuDoold(NCold)); uuDoold = 0.0
    allocate (uuXold(NCold));  uuXold  = 0.0
    allocate (uuXoold(NCold)); uuXoold = 0.0

    !vv-------------------------------
    allocate (vvold(NCold));   vvold   = 0.0
    allocate (vvoold(NCold));  vvoold  = 0.0
    allocate (vvCold(NCold));  vvCold  = 0.0
    allocate (vvCoold(NCold)); vvCoold = 0.0
    allocate (vvDoold(NCold)); vvDoold = 0.0
    allocate (vvXold(NCold));  vvXold  = 0.0
    allocate (vvXoold(NCold)); vvXoold = 0.0

    !ww-------------------------------
    allocate (wwold(NCold));   wwold   = 0.0
    allocate (wwoold(NCold));  wwoold  = 0.0
    allocate (wwCold(NCold));  wwCold  = 0.0
    allocate (wwCoold(NCold)); wwCoold = 0.0
    allocate (wwDoold(NCold)); wwDoold = 0.0
    allocate (wwXold(NCold));  wwXold  = 0.0
    allocate (wwXoold(NCold)); wwXoold = 0.0

    !uv-------------------------------
    allocate (uvold(NCold));   uvold   = 0.0
    allocate (uvoold(NCold));  uvoold  = 0.0
    allocate (uvCold(NCold));  uvCold  = 0.0
    allocate (uvCoold(NCold)); uvCoold = 0.0
    allocate (uvDoold(NCold)); uvDoold = 0.0
    allocate (uvXold(NCold));  uvXold  = 0.0
    allocate (uvXoold(NCold)); uvXoold = 0.0

    !uw-------------------------------
    allocate (uwold(NCold));   uwold   = 0.0
    allocate (uwoold(NCold));  uwoold  = 0.0
    allocate (uwCold(NCold));  uwCold  = 0.0
    allocate (uwCoold(NCold)); uwCoold = 0.0
    allocate (uwDoold(NCold)); uwDoold = 0.0
    allocate (uwXold(NCold));  uwXold  = 0.0
    allocate (uwXoold(NCold)); uwXoold = 0.0

    !vw-------------------------------
    allocate (vwold(NCold));   vwold   = 0.0
    allocate (vwoold(NCold));  vwoold  = 0.0
    allocate (vwCold(NCold));  vwCold  = 0.0
    allocate (vwCoold(NCold)); vwCoold = 0.0
    allocate (vwDoold(NCold)); vwDoold = 0.0
    allocate (vwXold(NCold));  vwXold  = 0.0
    allocate (vwXoold(NCold)); vwXoold = 0.0
  
    !Eps-----------------------------
    allocate (Epsold(NCold));   Epsold   = 0.0
    allocate (Epsoold(NCold));  Epsoold  = 0.0
    allocate (EpsCold(NCold));  EpsCold  = 0.0
    allocate (EpsCoold(NCold)); EpsCoold = 0.0
    allocate (EpsDoold(NCold)); EpsDoold = 0.0
    allocate (EpsXold(NCold));  EpsXold  = 0.0
    allocate (EpsXoold(NCold)); EpsXoold = 0.0

    !f22-----------------------------
    allocate (f22old(NCold));   f22old   = 0.0
    allocate (f22oold(NCold));  f22oold  = 0.0
    allocate (f22Doold(NCold)); f22Doold = 0.0
    allocate (f22Xold(NCold));  f22Xold  = 0.0
    allocate (f22Xoold(NCold)); f22Xoold = 0.0

    j = NCold
    do k = 1, j
        if(this < 2) then
            if(mod(k,20000) == 0) write(*,*) (100.*k/(1.*j)), '% complete...'
        end if
        read(5,'(78ES18.8)')  Xold(k), Yold(k), Zold(k), &
            Uold(k), Uoold(k), UCold(k), UCoold(k), UDoold(k), UXold(k), UXoold(k), &
            Vold(k), Voold(k), VCold(k), VCoold(k), VDoold(k), VXold(k), VXoold(k), &
            Wold(k), Woold(k), WCold(k), WCoold(k), WDoold(k), WXold(k), WXoold(k), &
            Pold(k), PPold(k), Pxold(k), Pyold(k), Pzold(k), &
            uuold(k), uuoold(k), uuCold(k), uuCoold(k), uuDoold(k), uuXold(k), uuXoold(k), &
            vvold(k), vvoold(k), vvCold(k), vvCoold(k), vvDoold(k), vvXold(k), vvXoold(k), &
            wwold(k), wwoold(k), wwCold(k), wwCoold(k), wwDoold(k), wwXold(k), wwXoold(k), &
            uvold(k), uvoold(k), uvCold(k), uvCoold(k), uvDoold(k), uvXold(k), uvXoold(k), &
            uwold(k), uwoold(k), uwCold(k), uwCoold(k), uwDoold(k), uwXold(k), uwXoold(k), &
            vwold(k), vwoold(k), vwCold(k), vwCoold(k), vwDoold(k), vwXold(k), vwXoold(k), &
            Epsold(k), Epsoold(k), EpsCold(k), EpsCoold(k), EpsDoold(k), EpsXold(k), EpsXoold(k)!, &
                            !f22old(k), f22oold(k), f22Doold(k), f22Xold(k), f22Xoold(k)
    end do
    close(5)
    if(this < 2) write(*,*) 'LoaInI: finished with reading the files'

    nearest = 0
    near = 0
    DISTold = HUGE
    do c = 1, NC
        if(this < 2) then
            if(mod(c,20000) == 0) write(*,*) (100.*c/(1.*NC)), '% complete...'
        end if
        DISTold = HUGE
        do k = 1, j
            if(Dist(Xold(k),Yold(k),Zold(k),xc(c),yc(c),zc(c)) < DISTold) then
                DISTold = Dist(Xold(k),Yold(k),Zold(k),xc(c),yc(c),zc(c))
                nearest =  k
            end if
        end do
        !U---------------------------
        U % n(c)  = Uold(nearest)
        U % o(c)  = Uoold(nearest)
        U % C(c)  = UCold(nearest)
        U % X(c)  = UXold(nearest)
        !V---------------------------
        V % n(c)  = Vold(nearest)
        V % o(c)  = Voold(nearest)
        V % C(c)  = VCold(nearest)
        V % X(c)  = VXold(nearest)
        !W---------------------------
        W % n(c)  = Wold(nearest)
        W % o(c)  = Woold(nearest)
        W % C(c)  = WCold(nearest)
        W % X(c)  = WXold(nearest)
        !P---------------------------
        P % n(c)  = Pold(nearest)
        PP % n(c) = PPold(nearest)
        Px(c)     = Pxold(nearest)
        Py(c)     = Pyold(nearest)
        Pz(c)     = Pzold(nearest)

        !uu--------------------------
        uu % n(c) = uuold(nearest)
        uu % o(c) = uuoold(nearest)
        uu % C(c) = uuCold(nearest)
        uu % X(c) = uuXold(nearest)
        !vv--------------------------
        vv % n(c) = vvold(nearest)
        vv % o(c) = vvoold(nearest)
        vv % C(c) = vvCold(nearest)
        vv % X(c) = vvXold(nearest)
        !ww--------------------------
        ww % n(c) = wwold(nearest)
        ww % o(c) = wwoold(nearest)
        ww % C(c) = wwCold(nearest)
        ww % X(c) = wwXold(nearest)
        !uv--------------------------
        uv % n(c) = uvold(nearest)
        uv % o(c) = uvoold(nearest)
        uv % C(c) = uvCold(nearest)
        uv % X(c) = uvXold(nearest)
        !uw--------------------------
        uw % n(c) = uwold(nearest)
        uw % o(c) = uwoold(nearest)
        uw % C(c) = uwCold(nearest)
        uw % X(c) = uwXold(nearest)
        !vw--------------------------
        vw % n(c) = vwold(nearest)
        vw % o(c) = vwoold(nearest)
        vw % C(c) = vwCold(nearest)
        vw % X(c) = vwXold(nearest)
    end do
    do s=1, NS
        c1=SideC(1,s)
        c2=SideC(2,s)

        !---- interpoliraj gustocu i brzine
        Us = f(s) * U % n(c1) + (1.0-f(s)) * U % n(c2)
        Vs = f(s) * V % n(c1) + (1.0-f(s)) * V % n(c2)
        Ws = f(s) * W % n(c1) + (1.0-f(s)) * W % n(c2)
        Flux(s) = ( Us*Sx(s) + Vs*Sy(s) + Ws*Sz(s) )
    end do

    name = answer

    !xyz--------------------------
    deallocate(Xold);
    deallocate(Yold);
    deallocate(Zold);
    !UVW--------------------------
    deallocate(Uold);
    deallocate(Vold);
    deallocate(Wold);
    deallocate(Uoold);
    deallocate(Voold);
    deallocate(Woold);
    deallocate(UDoold);
    deallocate(VDoold);
    deallocate(WDoold);
    deallocate(UCold);
    deallocate(VCold);
    deallocate(WCold);
    deallocate(UCoold);
    deallocate(VCoold);
    deallocate(WCoold);
    deallocate(UXold);
    deallocate(VXold);
    deallocate(WXold);
    deallocate(UXoold);
    deallocate(VXoold);
    deallocate(WXoold);
    !P-----------------------------
    deallocate(Pold);
    deallocate(PPold);
    deallocate(Pxold);
    deallocate(Pyold);
    deallocate(Pzold);
    !stresses---------------------------
    !uu-------------------------------
    deallocate(uuold);
    deallocate(uuoold);
    deallocate(uuDoold);
    deallocate(uuCold);
    deallocate(uuCoold);
    deallocate(uuXold);
    deallocate(uuXoold);
    !vv-------------------------------
    deallocate(vvold);
    deallocate(vvoold);
    deallocate(vvDoold);
    deallocate(vvCold);
    deallocate(vvCoold);
    deallocate(vvXold);
    deallocate(vvXoold);
    !ww-------------------------------
    deallocate(wwold);
    deallocate(wwoold);
    deallocate(wwDoold);
    deallocate(wwCold);
    deallocate(wwCoold);
    deallocate(wwXold);
    deallocate(wwXoold);
    !uv-------------------------------
    deallocate(uvold);
    deallocate(uvoold);
    deallocate(uvDoold);
    deallocate(uvCold);
    deallocate(uvCoold);
    deallocate(uvXold);
    deallocate(uvXoold);
    !uw-------------------------------
    deallocate(uwold);
    deallocate(uwoold);
    deallocate(uwDoold);
    deallocate(uwCold);
    deallocate(uwCoold);
    deallocate(uwXold);
    deallocate(uwXoold);
    !vw-------------------------------
    deallocate(vwold);
    deallocate(vwoold);
    deallocate(vwDoold);
    deallocate(vwCold);
    deallocate(vwCoold);
    deallocate(vwXold);
    deallocate(vwXoold);
  
    !Eps-----------------------------
    deallocate(Epsold);
    deallocate(Epsoold);
    deallocate(EpsDoold);
    deallocate(EpsCold);
    deallocate(EpsCoold);
    deallocate(EpsXold);
    deallocate(EpsXoold);
    !f22-----------------------------
    deallocate(f22old);
    deallocate(f22oold);
    deallocate(f22Doold);
    deallocate(f22Xold);
    deallocate(f22Xoold);

    write(*,*) 'Finished with LoaIni  Processor: ', this

    !---- restore the name
    name = answer

END SUBROUTINE LoaIni

