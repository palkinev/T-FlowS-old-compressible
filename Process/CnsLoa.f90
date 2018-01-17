!======================================================================!
SUBROUTINE CnsLoa
    !----------------------------------------------------------------------!
    ! Reads:  NAME.cns                                                     !
    ! ~~~~~~                                                               !
    !------------------------------[Modules]-------------------------------!
    USE all_mod
    USE pro_mod
    USE sol_mod
    USE les_mod
    USE par_mod

    !----------------------------------------------------------------------!
    IMPLICIT NONE
    !-------------------------------[Locals]-------------------------------!
    INTEGER       :: c, s, it
    CHARACTER(80) :: nameIn
    CHARACTER(8)  :: answer
    !--------------------------------[CVS]---------------------------------!
    !  $Id: CnsLoa.f90,v 1.7 2009/06/30 12:01:16 IUS\mhadziabdic Exp $
    !  $Source: /home/IUS/mhadziabdic/.CVSROOT/T-Rex/Process/CnsLoa.f90,v $
    !======================================================================!

    if(this < 2) write(*,*) '# Input problem name:'
    call ReadC(7,inp,tn,ts,te)
    read(inp(ts(1):te(1)), '(A80)')  name

    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
    !       Read the file with the      !
    !     connections between cells     !
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
    call NamFil(THIS, nameIn, '.cns', len_trim('.cns'))

    open(9, FILE=nameIn,FORM='UNFORMATTED')
    if(this < 2) write(*,*) '# Now reading the file:', nameIn

    !///// number of cells, boundary cells and sides
    read(9) NC
    read(9) NbC
    read(9) NS
    read(9) c   ! NSsh is not used in Processor.
    read(9) Nmat

    !///// cell materials
    allocate (material(-NbC:NC))
    read(9) (material(c), c=1,NC)
    read(9) (material(c), c=-1,-NBC,-1)

    !///// sides
    allocate (SideC(0:2,NS))
    read(9) (SideC(0,s), s=1,NS)
    read(9) (SideC(1,s), s=1,NS)
    read(9) (SideC(2,s), s=1,NS)

    !///// boundary cells
  
    !allocate (TypeBC(-NbC:-1)); TypeBC = 0
    allocate (TypeBC(-NbC:NC)); TypeBC = 0
    allocate (bcmark(-NbC:-1))
    read(9) (bcmark(c), c=-1,-NbC, -1)

    !///// boundary copy cells
    allocate (CopyC(-NbC:-1))
    read(9) (CopyC(c), c=-1,-NbC, -1)

    close(9)

    !----- Type of the problem
    if(this  < 2) then
        write(*,*) '# Type of problem: '
        write(*,*) '# CHANNEL          -> Channel flow'
        write(*,*) '# PIPE             -> Pipe flow'
        write(*,*) '# JET              -> Impinging jet flow'
        write(*,*) '# TEST             -> Test Laplacian equation'
        write(*,*) '# OTHER            -> All the other problems'
        write(*,*) '# HOT              -> Problems with temperature'
        write(*,*) '# XHOM, YHOM, ZHOM -> Homogeneous directions'
        write(*,*) '# TGV -> Taylor-Green Vortex test case' ! actually used for combustion
        write(*,*) '# COHERENT -> Gather coherent stats (check CalcMn)'
    endif
    call ReadC(7,inp,tn,ts,te)
    do it=1,tn
        read(inp(ts(it):te(it)),'(A8)')  answer
        call ToUppr(answer)
        if(answer == 'CHANNEL') then
            CHANNEL = YES
        else if(answer == 'PIPE') then
            PIPE = YES
        else if(answer == 'JET') then
            JET = YES
        else if(answer == 'TEST') then
            TEST = YES
        else if(answer == 'OTHER') then
            OTHER = YES
        else if(answer == 'HOT') then
            HOT = YES
        else if(answer == 'XHOM') then
            XHOM = YES
        else if(answer == 'YHOM') then
            YHOM = YES
        else if(answer == 'ZHOM') then
            ZHOM = YES
        else if(answer == 'MOIN_1') then ! moin problem 1
            MoinID = 1
        else if(answer == 'MOIN_2') then ! moin problem 2
            MoinID = 2
        else if(answer == 'MOIN_3') then ! moin problem 3
            MoinID = 3
        else if(answer == 'MOIN_4') then ! moin problem 4 (RT)
            MoinID = 4
        else if(answer == 'MOIN_1_T') then ! moin problem 1 with table values
            MoinID = 10
        else if(answer == 'MOIN_3_T') then ! moin problem 3 with table values
            MoinID = 30            
        !    else if(answer == 'TGV') then
        !      TGV = YES
        else if(answer == 'COHERENT') then
            COHERENT = YES
        else if(answer == 'ROT') then
            ROT = YES
        else if(answer == 'BUDG') then
            BUDG = YES
        else
            write(*,*) 'Error in input ! Exiting'
            stop
        endif
    end do

END SUBROUTINE CnsLoa
