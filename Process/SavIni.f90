!======================================================================!
SUBROUTINE SavIni()
    !----------------------------------------------------------------------!
    !
    !----------------------------------------------------------------------!
    !------------------------------[Modules]-------------------------------!
    USE all_mod
    USE pro_mod
    USE les_mod
    USE par_mod
    USE rans_mod
    !----------------------------------------------------------------------!
    IMPLICIT NONE
    !-------------------------------[Locals]-------------------------------!
    INTEGER   ::  c, nn
    CHARACTER :: namOut*80, answer*80, ext*4
    !--------------------------------[CVS]---------------------------------!
    character(80) :: rcs1,rcs2
    data rcs1/                                                        &
        '$Id: SavIni.f90,v 1.4 2008/12/10 14:57:50 IUS\mhadziabdic Exp $'/
    data rcs2/                                                        &
        '$Source: /home/IUS/mhadziabdic/.CVSROOT/T-Rex/Process/SavIni.f90,v $'/
    !======================================================================!

    !---- store the name
    if(this  < 2)                                                     &
        write(*,*) '# Now saving initial files [skip cancels]:'
    call ReadC(7,inp,tn,ts,te)
  
    read(inp(ts(1):te(1)), '(A80)')  namOut
    answer=namOut
    call ToUppr(answer)
  
    if(answer == 'SKIP') return
  
    !---- save the name
    answer = name
    name = namOut
    nn = 0
    ext = '.xyz'
    call NamFil(THIS, namOut, ext, len_trim(ext))
    open(9,FILE=namOut)
    do c= 1, NC
        nn = nn + 1
    end do    ! through centers 
    write(9,'(I10)') nn
    do c= 1, NC
        write(9,'(3E25.8)') xc(c),yc(c),zc(c)
    end do    ! through centers 
    close(9)

    ext = '.U__'
    call NamFil(THIS, namOut, ext, len_trim(ext))
    open(9,FILE=namOut)
    do c= 1, NC
        write(9,'(7E18.8)') U % n(c), U % o(c), U % C(c), &
            U % X(c)
    end do    ! through centers 
    close(9)

    ext = '.V__'
    call NamFil(THIS, namOut, ext, len_trim(ext))
    open(9,FILE=namOut)
    do c= 1, NC
        write(9,'(7E18.8)') V % n(c), V % o(c), V % C(c), &
            V % X(c)
    end do    ! through centers 
    close(9)

    ext = '.W__'
    call NamFil(THIS, namOut, ext, len_trim(ext))
    open(9,FILE=namOut)
    do c= 1, NC
        write(9,'(7E18.8)') W % n(c), W % o(c), W % C(c), &
            W % X(c)
    end do    ! through centers 
    close(9)

    ext = '.P__'
    call NamFil(THIS, namOut, ext, len_trim(ext))
    open(9,FILE=namOut)
    do c= 1, NC
        write(9,'(5E18.8)') P % n(c), PP % n(c), Px(c), Py(c),  &
            Pz(c)
    end do    ! through centers 
    close(9)
 
    if(HOT == YES) then 
        ext = '.T__'
        call NamFil(THIS, namOut, ext, len_trim(ext))
        open(9,FILE=namOut)
        do c= 1, NC
            write(9,'(7E18.8)') T % n(c), T % o(c), T % C(c), &
                T % X(c)
        end do    ! through centers
        close(9)
    end if 

    name = answer

END SUBROUTINE SavIni
