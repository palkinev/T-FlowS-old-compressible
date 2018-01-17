!======================================================================!
SUBROUTINE BouLoa(OVERWRITE_BC)
    !----------------------------------------------------------------------!
    ! Reads: NAME.b                                                        !
    ! ~~~~~~                                                               !
    !------------------------------[Modules]-------------------------------!
    USE allp_mod, only: buffer, convect, fluid, huge, inflow, outflow, pressure, solid, symmetry, wall, wallfl
    USE all_mod
    USE pro_mod
    USE rans_mod
    USE par_mod
    !----------------------------------------------------------------------!
    IMPLICIT NONE
    !-----------------------------[Parameters]-----------------------------!
    LOGICAL       :: OVERWRITE_BC
    !-------------------------------[Locals]-------------------------------!
    INTEGER       :: c, n, dum1, NB, NP, Ninit, m
    CHARACTER(80) :: namBou, namPro(128), dir
    INTEGER       :: typBou(128)
    REAL          :: xyz(1024)
    REAL          :: wi
    LOGICAL       :: here
    !--------------------------------[CVS]---------------------------------!
    !  $Id: BouLoa.f90,v 1.38 2009/06/30 11:42:00 IUS\mhadziabdic Exp $
    !  $Source: /home/IUS/mhadziabdic/.CVSROOT/T-Rex/Process/BouLoa.f90,v $
    !======================================================================!

    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
    !     Read the file with boundary conditions     !
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
    namBou = name
    namBou(len_trim(name)+1:len_trim(name)+2) = '.b'
    open(9, FILE=namBou)
    if(this < 2) write(*,*) '# Now reading the file:', namBou

    !---------------------!
    ! Phisical properties !
    !---------------------!
    call ReadC(9,inp,tn,ts,te)
    read(inp,*) Nmat
    do n=1,Nmat
        call ReadC(9,inp,tn,ts,te)
        call ToUppr(  inp(ts(2):te(2))  )
        if( inp(ts(2):te(2))  ==  'FLUID') then
            StateMat(n)=FLUID
        else if( inp(ts(2):te(2))  ==  'SOLID') then
            StateMat(n)=SOLID
        else
            if(this < 2) write(*,*) 'BouLoa: Unknown material state'
            stop
        end if
        read(inp(ts(3):te(3)),*) VISc
        read(inp(ts(4):te(4)),*) DENc(n)
        if(HOT==YES) read(inp(ts(5):te(5)),*) CONc(n)
        if(HOT==YES) read(inp(ts(6):te(6)),*) CAPc(n)
    end do
    !---------------------!
    ! Boundary conditions !
    !---------------------!
    call ReadC(9,inp,tn,ts,te)
    read(inp,*) NB
    do n=1,NB
        call ReadC(9,inp,tn,ts,te)
        call ToUppr(  inp(ts(2):te(2))  )
        call ToUppr(  inp(ts(3):te(3))  )
        read(inp(ts(1):te(1)),*) dum1
        if( inp(ts(2):te(2)) == 'INFLOW') then
            typBou(n)=INFLOW
            if (this < 2) write(*,*) '# INFLOW:'
        else if( inp(ts(2):te(2)) == 'WALL') then
            typBou(n)=WALL
            if (this < 2) write(*,*) '# WALL:'
        else if( inp(ts(2):te(2)) == 'OUTFLOW') then
            typBou(n)=OUTFLOW
            if (this < 2) write(*,*) '# OUTFLOW:'
        else if( inp(ts(2):te(2)) == 'SYMMETRY') then
            typBou(n)=SYMMETRY
            if (this < 2) write(*,*) '# SYMMETRY:'
        else if( inp(ts(2):te(2)) == 'WALLFLUX') then
            typBou(n)=WALLFL
            if (this < 2) write(*,*) '# WALLFLUX:'
        else if( inp(ts(2):te(2)) == 'CONVECTIVE') then
            typBou(n)=CONVECT
            if (this < 2) write(*,*) '# CONVECTIVE:'
        else if( inp(ts(2):te(2)) == 'PRESSURE') then
            typBou(n)=PRESSURE
            if (this < 2) write(*,*) '# PRESSURE:'
        else
            if(this < 2) write(*,*) 'BouLoa: Unknown boundary condition type: ', inp(ts(2):te(2))
            stop
        end if
        if( inp(ts(3):te(3))  ==  'FILE') then
            read(inp(ts(4):te(4)),'(A80)') namPro(n)
        else
            if (this<2) write(*,*) '# Reading boundary condition : ', inp(ts(2):te(2)), ' :'
            read(inp(ts(3):te(3)),*) U % bound(n)
            if (this < 2) write(*,*) '# U     :', U % bound(n)
            read(inp(ts(4):te(4)),*) V % bound(n)
            if (this < 2) write(*,*) '# V     :', V % bound(n)
            read(inp(ts(5):te(5)),*) W % bound(n)
            if (this < 2) write(*,*) '# W     :', W % bound(n)
            if(typBou(n)==PRESSURE) then
                read(inp(ts(6):te(6)),*) P % bound(n)
                if (this < 2) write(*,*) '# P     :', P % bound(n)
                if(HOT==YES) then
                    !          read(inp(ts(7):te(7)),*) T % bound(n)
                    read(inp(ts(7):te(7)),*) Zmix % bound(n)
                    if (this < 2) write(*,*) '# Zmix  :', Zmix % bound(n)
                else
                end if
                namPro(n)=''
            else
                if(HOT==YES) then
                    !          read(inp(ts(6):te(6)),*) T % bound(n)
                    read(inp(ts(6):te(6)),*) Zmix % bound(n)
                    if (this < 2) write(*,*) '# Zmix  :', Zmix % bound(n)
                else
                end if
                namPro(n)=''
            end if
        end if
    end do

    !  write(*,*) 'before init??'

    !--------------------!
    ! Initial conditions !
    !--------------------!
    call ReadC(9,inp,tn,ts,te)
    read(inp,*) Ninit
    if(Ninit > Nmat) then
        if(this < 2) write(*,*) 'Warning: there are more initial conditions then materials'
    end if

    do n=1,Ninit
        call ReadC(9,inp,tn,ts,te)
        call ToUppr(inp(ts(2):te(2)))

        !-----Initial conditions given in GMV file
        if(inp(ts(2):te(2)) == 'FILE') then
            read(inp(ts(3):te(3)),'(A80)') namIni(n)
            write(*,*) '@BouLoa: material ', n, '; init. cond. given by file: ', namIni(n)
        else
            namIni(n) = ''

            !-----Initial conditions given by constant
            if (this < 2) write(*,*) '# Reading initial conditions'
            read(inp(ts(2):te(2)),*) U % init(n)
            if (this < 2) write(*,*) '# U     :', U % init(n)
            read(inp(ts(3):te(3)),*) V % init(n)
            if (this < 2) write(*,*) '# V     :', V % init(n)
            read(inp(ts(4):te(4)),*) W % init(n)
            if (this < 2) write(*,*) '# W     :', W % init(n)
 
            if(HOT==YES) then
                !      read(inp(ts(5):te(5)),*) T % init(n)
                read(inp(ts(5):te(5)),*) Zmix % init(n)
                if (this < 2) write(*,*) '# Zmix  :', Zmix % init(n)
            else ! HOT /= YES
            end if
        end if
    end do

    close(9)

    !----------------------------------!
    !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
    !----------------------------------!
    do n=1,NB

        !---- Boundary condition is given by a single constant
        if(namPro(n) == '') then
            do c=-1,-NbC,-1
                if(bcmark(c) == n) then
                    TypeBC(c) = typBou(n)

                    !===== if OVERWRITE_BC is set to true, set boundary values,
                    !===== otherwise, just the TypeBC remains set.

                    if(OVERWRITE_BC) then
                        U % n(c) = U % bound(n)
                        V % n(c) = V % bound(n)
                        W % n(c) = W % bound(n)
                        !Ur % n(c)= U % n(c) * Rho % n(c)
                        !Vr % n(c)= V % n(c) * Rho % n(c)
                        !Wr % n(c)= W % n(c) * Rho % n(c)
                        P % n(c) = P % bound(n)
                        if(HOT == YES) then
                            if(TypeBC(c).eq.WALLFL) then
                                !                T % q(c) =  T % bound(n)
                                Zmix % q(c) =  Zmix % bound(n)         !??
                              !Zmixr % q(c) =  Zmix % q(c) * Rho % n(c) !??
                            else
                                !                T % n(c) =  T % bound(n)
                                Zmix % n(c) =  Zmix % bound(n)
                                !Zmixr % n(c) =  Zmix % n(c) * Rho % n(c)
                            endif
                        end if  ! for HOT==YES
                    end if
                end if
            end do
        !---- Boundary condition is prescribed in a file
        else
            open(9, FILE=namPro(n))
            if(this < 2) write(*,*) '# Now reading the file:', namPro(n)
            call ReadC(9,inp,tn,ts,te)
            read(inp(ts(1):te(1)),*) NP                  ! number of points
            if(NP  > 1000) then
                if(this < 2) write(*,*) 'BouLoa: Too much points in a profile !'
                stop
            end if
            call ReadC(9,inp,tn,ts,te)
            read(inp(ts(1):te(1)),*) dir  ! direction

            do m=1,NP
                call ReadC(9,inp,tn,ts,te)
                read(inp(ts(1):te(1)),*) xyz(m)
                read(inp(ts(2):te(2)),*) U % pro(m)
                read(inp(ts(3):te(3)),*) V % pro(m)
                read(inp(ts(4):te(4)),*) W % pro(m)
                read(inp(ts(5):te(5)),*) Zmix % pro(m)
            end do

            do c=-1,-NbC,-1
                if(bcmark(c) == n) then
                    TypeBC(c) = typBou(n)
          
                    !===== if OVERWRITE_BC is set to true, set boundary values,
                    !===== otherwise, just the TypeBC remains set.

                    if(OVERWRITE_BC) then
                        do m=1,NP-1
                            here = .FALSE.
                            !----- compute the weight factors
                            if( (dir == 'X' .or. dir == 'x') .and.                  &
                                xc(c) >= xyz(m) .and. xc(c) <= xyz(m+1) ) then
                                wi = ( xyz(m+1)-xc(c) ) / ( xyz(m+1) - xyz(m) )
                                here = .TRUE.
                            else if( (dir == 'Y' .or. dir == 'y') .and.             &
                                yc(c) >= xyz(m) .and. yc(c) <= xyz(m+1) ) then
                                wi = ( xyz(m+1)-yc(c) ) / ( xyz(m+1) - xyz(m) )
                                here = .TRUE.
                            else if( (dir == 'Z' .or. dir == 'z') .and.             &
                                zc(c) >= xyz(m) .and. zc(c) <= xyz(m+1) ) then
                                wi = ( xyz(m+1)-zc(c) ) / ( xyz(m+1) - xyz(m) )
                                here = .TRUE.
                            else if( (dir == 'RX' .or. dir == 'rx') .and.           &
                                sqrt(yc(c)*yc(c)+zc(c)*zc(c)) >= xyz(m) .and.      &
                                sqrt(yc(c)*yc(c)+zc(c)*zc(c)) <= xyz(m+1) ) then
                                wi = ( xyz(m+1) - sqrt(yc(c)*yc(c)+zc(c)*zc(c)) )     &
                                    / ( xyz(m+1) - xyz(m) )
                                here = .TRUE.
                            else if( (dir == 'RY' .or. dir == 'ry') .and.           &
                                sqrt(xc(c)*xc(c)+zc(c)*zc(c)) >= xyz(m) .and.      &
                                sqrt(xc(c)*xc(c)+zc(c)*zc(c)) <= xyz(m+1) ) then
                                wi = ( xyz(m+1) - sqrt(xc(c)*xc(c)+zc(c)*zc(c)) )     &
                                    / ( xyz(m+1) - xyz(m) )
                                here = .TRUE.
                            else if( (dir == 'RZ' .or. dir == 'rz') .and.           &
                                sqrt(xc(c)*xc(c)+yc(c)*yc(c)) <= xyz(m) .and.      &
                                sqrt(xc(c)*xc(c)+yc(c)*yc(c)) >= xyz(m+1) ) then
                                wi = ( xyz(m+1) - sqrt(xc(c)*xc(c)+yc(c)*yc(c)) )     &
                                    / ( xyz(m+1) - xyz(m) )
                                here = .TRUE.
                            end if

                            !----- interpolate the profiles
                            if(here) then
                                U % n(c) = wi*U % pro(m) + (1.-wi)*U % pro(m+1)
                                V % n(c) = wi*V % pro(m) + (1.-wi)*V % pro(m+1)
                                W % n(c) = wi*W % pro(m) + (1.-wi)*W % pro(m+1)
                                !Ur % n(c) = U % n(c) * Rho % n(c)
                                !Vr % n(c) = V % n(c) * Rho % n(c)
                                !Wr % n(c) = W % n(c) * Rho % n(c)
                                Zmix % n(c) = wi*Zmix % pro(m) + (1.-wi)*Zmix % pro(m+1)
                                !Zmixr % n(c) = Zmix % n(c) * Rho % n(c)
                                !                if(HOT==YES) &
                                !                  T % n(c) = wi*T % pro(m) + (1.-wi)*T % pro(m+1)
                                !                  Zmix % n(c) = wi*Zmix % pro(m) + (1.-wi)*Zmix % pro(m+1)
                                !                  Zmixr % n(c) = Zmix % n(c) * Rho % n(c)
                            end if
                        end do
                    end if  ! if(OVERWRITE_BC)
                end if
            end do
            close(9)
        end if
    end do

    !---------------------------------!
    ! Finally handle the buffer cells !
    !---------------------------------!
    do c=-1,-NbC,-1
        if(bcmark(c) == BUFFER) TypeBC(c)=BUFFER
    end do

    RETURN

END SUBROUTINE BouLoa
