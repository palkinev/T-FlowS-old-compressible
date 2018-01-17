SUBROUTINE SavParView_xml_bin_ascii_old(sub, NCsub, nameIN)
    !----------------------------------------------------------------------!
    ! Reads: NAME.gmv and generates NAME.vti Paraview XML output file      !
    ! now supports ascii or binary format                                  !
    !------------------------------[Modules]-------------------------------!
    USE all_mod
    USE par_mod
    USE allp_mod
    USE pro_mod
    USE les_mod
    USE rans_mod

    !USE vtk_fortran ! VTKFortran-current and latest
    USE LIB_VTK_IO ! VTKFortran-1.1.1

    IMPLICIT NONE

    INCLUDE 'mpif.h'

    integer             :: error_MPI

    CHARACTER           :: nameIN*(*)
    INTEGER             :: NCsub, sub, sub_local

    CHARACTER           :: stringadummy*100, nameTMP*80, nameTMP2*80
    INTEGER             :: NNsub, NCsub_new
    INTEGER             :: off_set_connection
    INTEGER,ALLOCATABLE :: connessione(:,:) ! connection
    INTEGER             :: celleconnessione
    INTEGER             :: c, n, s, i
    REAL,ALLOCATABLE    :: x(:), y(:), z(:)    ! self evident
    ! mesh 2 VTK_IO vars
    integer*1, dimension(1:NCsub) :: cell_type_vtk_io
    integer*4, dimension(1:NCsub) :: offset_vtk_io, PID
    integer*4, allocatable        :: connect_vtk_io(:)
    integer*4                     :: E_IO
    ! auxilary vars
    real*8,   allocatable         :: stresses_2_vtk_io(:,:)
    real*8,   allocatable         :: Q_2_vtk_io(:)
    real*8,   allocatable         :: Q_LES_filt_x1_2_vtk_io(:)
    real*8,   allocatable         :: Q_LES_filt_x2_2_vtk_io(:)
    real*8,   allocatable         :: Q_LES_filt_x3_2_vtk_io(:)
    real*8,   allocatable         :: Q_LES_filt_x4_2_vtk_io(:)
    real                          :: Ua,Va,Wa,DE
    integer                       :: cj, j
    real*8,   allocatable         :: TT_2_vtk_io(:)
    real*8,   allocatable         :: pp_2_vtk_io(:)
    ! to cut new geometry
    integer                       :: NC_new
    REAL                          :: x1, x2, x3, x4
    REAL                          :: y1, y2, y3, y4
    REAL                          :: z1, z2
    logical                       :: FileExist
    integer*1, allocatable        :: SaveAreaIsEmpty(:)
    integer*1, allocatable        :: SaveAreaIsEmptySum(:)

    nameTMP2 = name ! task name

    ! feature :
    ! draw only zones, which are inside area specified in Paraview.dat
    ! if file does not not exists -> draw full zone as usual
    sub_local = 1
    if (sub .ne. 0) then ! multi processor case
        allocate(SaveAreaIsEmpty(1:Npro));
        allocate(SaveAreaIsEmptySum(1:Npro));
    else                 ! one processor case
        sub_local = 0
        sub = 1
        allocate(SaveAreaIsEmpty(1:1));
        allocate(SaveAreaIsEmptySum(1:1));
    end if

    SaveAreaIsEmpty = 0 ! 0 -> not empty
    SaveAreaIsEmptySum = 0 ! 0 -> not empty
    NC_new = 0

    inquire (file="ParaView.dat",exist=FileExist)
    if (FileExist) then

        open(119, FILE='ParaView.dat', status="old", action="read")
        call wait

        !     x3,y3           x2,y2             z1         z2
        !       |---------------\                |---------|
        !       |      ^ y       \               |         |
        !       |      |          \              |         | ---> z-plane
        !       |      ---> x      \             |         |
        !       |-------------------|            |---------|
        !     x1,y1               x4,y4         z1         z2

        if(this < 2) write(*,*) '# Now reading the file: Paraview.dat '
        call ReadC(119,inp,tn,ts,te)
        read(inp(ts(1):te(1)),*) x1
        read(inp(ts(2):te(2)),*) x2
        x3 = x1
        x4 = x2
        !if(tn==3 .or. tn==5)   c_start=-NbC !paraview do not store boundary cells

        if(tn==4)   read(inp(ts(3):te(3)),*) x3
        if(tn==4)   read(inp(ts(4):te(4)),*) x4

        call ReadC(119,inp,tn,ts,te)
        read(inp(ts(1):te(1)),*) y1
        read(inp(ts(2):te(2)),*) y2
        y4 = y1
        y3 = y2
        if(tn==4)   read(inp(ts(3):te(3)),*) y3
        if(tn==4)   read(inp(ts(4):te(4)),*) y4

        call ReadC(119,inp,tn,ts,te)
        read(inp(ts(1):te(1)),*) z1
        read(inp(ts(2):te(2)),*) z2

        if(this < 2) write(*,*) '#  x:(', x1, ';', x2, ';', x3, ';', x4, ')'
        if(this < 2) write(*,*) '#  x:(', y1, ';', y2, ';', y3, ';', y4, ')'
        if(this < 2) write(*,*) '#  z:(', z1, ';', z2, ')'

        call wait

        ! filtering cells
        do c = 1,NCsub
            if (  xc(c) > x1                              .and. xc(c) < x2                              &
                .and. xc(c) > x3                              .and. xc(c) < x4                              &
                .and. yc(c) > y1 + (y4-y1)*(xc(c)-x1)/(x4-x1) .and. yc(c) < y2 + (y3-y2)*(xc(c)-x2)/(x3-x2) &
                .and. zc(c) > z1                              .and. zc(c) < z2                              ) then

                NC_new = NC_new + 1

            end if
        end do

        !write(*,*) 'PID =', this,' . ', NC_new, ' cells were found'

        if (NC_new == 0) then
            SaveAreaIsEmpty(sub) = 1
        end if

    end if ! else draw full area

    ! in case of single processor
    if (sub_local == 0)  sub = 0

    !<<<<<<<<<<<<<<<<<<<<<<<<<!
    !                         !
    !     reads GMV file      !
    !                         !
    !<<<<<<<<<<<<<<<<<<<<<<<<<!
    call NamFil(sub, nameTMP, '.gmv', len_trim('.gmv'))
    open(9, FILE = nameTMP)
    if (this <2) write(*,*) 'Now reading the file: ', nameTMP

    !---------------!
    !     start     !
    !---------------!
    call ReadC(9,inp,tn,ts,te)  !read 'gmvinput ascii' line

    !---------------!
    !     nodes     !
    !---------------!
    call ReadC(9,inp,tn,ts,te)
    read(inp(ts(1):te(1)),*) stringadummy
    read(inp(ts(2):te(2)),*) NNsub
    allocate(x(NNsub)); x=0.
    allocate(y(NNsub)); y=0.
    allocate(z(NNsub)); z=0.
  
    do n=1,NNsub
        call ReadC(9,inp,tn,ts,te)
        read(inp(ts(1):te(1)),*) x(n)
    end do
    do n=1,NNsub
        call ReadC(9,inp,tn,ts,te)
        read(inp(ts(1):te(1)),*) y(n)
    end do
    do n=1,NNsub
        call ReadC(9,inp,tn,ts,te)
        read(inp(ts(1):te(1)),*) z(n)
    end do

    !----------------------!
    !     cell section     !
    !----------------------!

    celleconnessione = 0
    off_set_connection = 0

    call ReadC(9,inp,tn,ts,te)
    read(inp(ts(1):te(1)),*) stringadummy
    read(inp(ts(2):te(2)),*) NCsub_new

    if (NCsub_new.ne.NCsub) then
        write(*,*) 'number of cells read and processed is different, exiting!'
        stop
    end if
  
    allocate(connessione(NCsub,9)); connessione=0

    do n=1,NCsub
        call ReadC(9,inp,tn,ts,te)
        read(inp(ts(1):te(1)),*) stringadummy
        read(inp(ts(2):te(2)),*) off_set_connection
        if (n==1) then
            connessione(n,1) = off_set_connection
        else
            connessione(n,1) = connessione (n-1,1) + off_set_connection
        end if
    
        celleconnessione = celleconnessione + off_set_connection

        call ReadC(9,inp,tn,ts,te)
        do c=1,off_set_connection
            read(inp(ts(c):te(c)),*) connessione(n,c+1)
        end do
    end do

    close(9)

    ! BEGIN - FOR VTK_IO mesh construction

    ! filling nodes connectivity
    s = 1
    allocate(connect_vtk_io(1:celleconnessione)); connect_vtk_io(:) = 0

    do c=1,NCsub
        !hex
        if (connessione(c,9) /= 0) then
            do i=2,9
                connect_vtk_io(s) = connessione(c,i)-1; s = s + 1
              !write(*,*) 'connect_vtk_io(',s-1,')=',connect_vtk_io(s-1), 'connessione(c,i)-1=', connessione(c,i)-1
            end do
        else
            !prism
            if ((connessione(c,8) == 0).and.(connessione(c,7) /= 0)) then
                connect_vtk_io(s) = connessione(c,2)-1; s = s + 1
                connect_vtk_io(s) = connessione(c,4)-1; s = s + 1
                connect_vtk_io(s) = connessione(c,3)-1; s = s + 1
                connect_vtk_io(s) = connessione(c,5)-1; s = s + 1
                connect_vtk_io(s) = connessione(c,7)-1; s = s + 1
                connect_vtk_io(s) = connessione(c,6)-1; s = s + 1
            else
                if ((connessione(c,7) == 0).and.(connessione(c,6) /= 0)) then ! pir ?
                    connect_vtk_io(s) = connessione(c,5)-1; s = s + 1
                    connect_vtk_io(s) = connessione(c,4)-1; s = s + 1
                    connect_vtk_io(s) = connessione(c,3)-1; s = s + 1
                    connect_vtk_io(s) = connessione(c,6)-1; s = s + 1
                    connect_vtk_io(s) = connessione(c,2)-1; s = s + 1
                else                                                         ! tetra ?
                    connect_vtk_io(s) = connessione(c,5)-1; s = s + 1
                    connect_vtk_io(s) = connessione(c,4)-1; s = s + 1
                    connect_vtk_io(s) = connessione(c,3)-1; s = s + 1
                    connect_vtk_io(s) = connessione(c,2)-1; s = s + 1
                end if
            end if
        end if

    end do

    ! filling cells offset_vtk_io
    do c=1,NCsub
        offset_vtk_io(c) = connessione(c,1)
    end do

    ! filling cell types
    do c=1,NCsub
        if ((connessione(c,6)==0)) then !tetra
            cell_type_vtk_io(c) = 10
        else
            if ((connessione(c,7)==0)) then
                cell_type_vtk_io(c) = 14          !pyramid
            else
                if ((connessione(c,8)==0)) then
                    cell_type_vtk_io(c) = 13      !prism
                else
                    cell_type_vtk_io(c) = 12      !hex
                end if
            end if
        end if
    end do

    ! -- Processors ID
    PID = 0
    PID(1:NCsub) = sub

    ! END - FOR VTK_IO mesh construction

    ! -- auxilary block ----------------

    call GraPhi(U % mean, 1, Ux,.TRUE.)    ! dU/dx
    call GraPhi(U % mean, 2, Uy,.TRUE.)    ! dU/dy
    call GraPhi(U % mean, 3, Uz,.TRUE.)    ! dU/dz
    call GraPhi(V % mean, 1, Vx,.TRUE.)    ! dV/dx
    call GraPhi(V % mean, 2, Vy,.TRUE.)    ! dV/dy
    call GraPhi(V % mean, 3, Vz,.TRUE.)    ! dV/dz
    call GraPhi(W % mean, 1, Wx,.TRUE.)    ! dW/dx
    call GraPhi(W % mean, 2, Wy,.TRUE.)    ! dW/dy
    call GraPhi(W % mean, 3, Wz,.TRUE.)    ! dW/dz


    !-- Q-criterium common
    call CalcShear(U % n, V % n, W % n, Shear)

    allocate(Q_2_vtk_io(1:NCsub)); Q_2_vtk_io(:) = 0. !< [1:NC_or_NN].

    Q_2_vtk_io(1:NCsub) = (Vort(1:NCsub)**2.0 - Shear(1:NCsub)**2.0)/4.0

    !-- Q-criterium LES for filtered velocity

    if (SIMULA == LES ) then

        call CalcShear(U % filt, V % filt, W % filt, Shear)

        !allocate(Q_LES_filt_x1_2_vtk_io(1:NCsub)); Q_LES_filt_x1_2_vtk_io(:) = 0. !< [1:NC_or_NN].

        !Q_LES_filt_x1_2_vtk_io(1:NCsub) = (Vort(1:NCsub)**2.0 - Shear(1:NCsub)**2.0)/4.0
  
        do c =1, NC
            Ua   = 0.0
            Va   = 0.0
            Wa   = 0.0
            DE   = 0.0
  
            do j = Acol(c), Acol(c + 1) - 1
                cj = Arow(j)
                if(cj /= c) then
  
                    !--- Test velocitys
                    Ua = Ua + volume(cj) * U % filt(cj)
                    Va = Va + volume(cj) * V % filt(cj)
                    Wa = Wa + volume(cj) * W % filt(cj)
  
                    !--- Test volume
                    DE = DE + volume(cj)
                end if
            end do
  
            !--- Now it's taking into account influence of central cell
            !--- within test molecule
      
            DE = DE + volume(c)
  
            Ua = Ua + volume(c) * U % filt(c)
            Va = Va + volume(c) * V % filt(c)
            Wa = Wa + volume(c) * W % filt(c)
  
            !--- Now calculating test values
            U % filt(c) = Ua / ( DE)
            V % filt(c) = Va / ( DE)
            W % filt(c) = Wa / ( DE)
        end do
  
        call CalcShear(U % filt, V % filt, W % filt, Shear)
  
        !allocate(Q_LES_filt_x2_2_vtk_io(1:NCsub)); Q_LES_filt_x2_2_vtk_io(:) = 0. !< [1:NC_or_NN].

        !Q_LES_filt_x2_2_vtk_io(1:NCsub) = (Vort(1:NCsub)**2.0 - Shear(1:NCsub)**2.0)/4.0
  
        do c =1, NC
            Ua   = 0.0
            Va   = 0.0
            Wa   = 0.0
            DE   = 0.0
  
            do j = Acol(c), Acol(c + 1) - 1
                cj = Arow(j)
                if(cj /= c) then
  
                    !--- Test velocitys
                    Ua = Ua + volume(cj) * U % filt(cj)
                    Va = Va + volume(cj) * V % filt(cj)
                    Wa = Wa + volume(cj) * W % filt(cj)
  
                    !--- Test volume
                    DE = DE + volume(cj)
                end if
            end do
  
            !--- Now it's taking into account influence of central cell
            !--- within test molecule
  
            DE = DE + volume(c)
  
            Ua = Ua + volume(c) * U % filt(c)
            Va = Va + volume(c) * V % filt(c)
            Wa = Wa + volume(c) * W % filt(c)
  
            !--- Now calculating test values
            U % filt(c) = Ua / ( DE)
            V % filt(c) = Va / ( DE)
            W % filt(c) = Wa / ( DE)
        end do
  
        call CalcShear(U % filt, V % filt, W % filt, Shear)
  
        !allocate(Q_LES_filt_x3_2_vtk_io(1:NCsub)); Q_LES_filt_x3_2_vtk_io(:) = 0. !< [1:NC_or_NN].

        !Q_LES_filt_x3_2_vtk_io(1:NCsub) = (Vort(1:NCsub)**2.0 - Shear(1:NCsub)**2.0)/4.0
  
        do c =1, NC
            Ua   = 0.0
            Va   = 0.0
            Wa   = 0.0
            DE   = 0.0
  
            do j = Acol(c), Acol(c + 1) - 1
                cj = Arow(j)
                if(cj /= c) then
  
                    !--- Test velocitys
                    Ua = Ua + volume(cj) * U % filt(cj)
                    Va = Va + volume(cj) * V % filt(cj)
                    Wa = Wa + volume(cj) * W % filt(cj)
  
                    !--- Test volume
                    DE = DE + volume(cj)
                end if
            end do
  
            !--- Now it's taking into account influence of central cell
            !--- within test molecule
  
            DE = DE + volume(c)
  
            Ua = Ua + volume(c) * U % filt(c)
            Va = Va + volume(c) * V % filt(c)
            Wa = Wa + volume(c) * W % filt(c)
  
            !--- Now calculating test values
            U % filt(c) = Ua / ( DE)
            V % filt(c) = Va / ( DE)
            W % filt(c) = Wa / ( DE)
        end do
  
        call CalcShear(U % filt, V % filt, W % filt, Shear)
  
        allocate(Q_LES_filt_x4_2_vtk_io(1:NCsub)); Q_LES_filt_x4_2_vtk_io(:) = 0. !< [1:NC_or_NN].

        Q_LES_filt_x4_2_vtk_io(1:NCsub) = (Vort(1:NCsub)**2.0 - Shear(1:NCsub)**2.0)/4.0
  
    end if ! Q_LES_filt

    !-- mean uu,vv,ww,uv,uw,vw stresses-------------------------------------------------------------------------------------------

    if (SIMULA == LES) then
        allocate(stresses_2_vtk_io(1:6,1:NCsub)); stresses_2_vtk_io(:,:) = 0. !< [1:N_COL,1:NC_or_NN].
                                      !|------------------------resolved-------------------------|
        stresses_2_vtk_io(1,1:NCsub) = UU  % mean(1:NCsub) - U % mean(1:NCsub) **2
        stresses_2_vtk_io(2,1:NCsub) = VV  % mean(1:NCsub) - V % mean(1:NCsub) **2
        stresses_2_vtk_io(3,1:NCsub) = WW  % mean(1:NCsub) - W % mean(1:NCsub) **2
        stresses_2_vtk_io(4,1:NCsub) = UV  % mean(1:NCsub) - U % mean(1:NCsub) * V % mean(1:NCsub)
        stresses_2_vtk_io(5,1:NCsub) = UW  % mean(1:NCsub) - U % mean(1:NCsub) * W % mean(1:NCsub)
        stresses_2_vtk_io(6,1:NCsub) = VW  % mean(1:NCsub) - V % mean(1:NCsub) * W % mean(1:NCsub)

    end if
    !-- mean tt stresses-------------------------------------------------------------------------------------------
    if (HOT == YES .and. MoinID .eq. 0) then
        allocate(TT_2_vtk_io(1:NCsub)); TT_2_vtk_io(:) = 0. !< [1:NC_or_NN].

        TT_2_vtk_io(1:NCsub) = TT % mean(1:NCsub) - T % mean(1:NCsub) **2
    end if
    !-- mean pp stresses-------------------------------------------------------------------------------------------

    if     (SIMULA == LES)  then
        allocate(pp_2_vtk_io(1:NCsub)); pp_2_vtk_io(:) = 0. !< [1:NC_or_NN].

        pp_2_vtk_io(1:NCsub) = PP % mean(1:NCsub) - P % mean(1:NCsub) **2
    end if

    call wait
    !write(6, *) 'Here is an error.. '
  
    ! -- MAIN BLOCK ----------------
  
    ! vtu - xml - parts
    ! ascii and binary are ok, raw is not
    ! use only float64 (base64) format, because vtk supports only base64 binary and int32
    name = nameIN
    call NamFil(sub, nameTMP, '.vtu', len_trim('.vtu'))
    if (this < 2 ) write(6, *) 'Now writing the file: ', nameTMP

    ! in case of single processor
    if (sub_local == 0)  sub = 1

    if (SaveAreaIsEmpty(sub) == 0) then ! do not write zones that do not cross area specified in Paraview.dat file
        E_IO = VTK_INI_XML(output_format = 'ascii', filename = nameTMP, mesh_topology = 'UnstructuredGrid')
        E_IO = VTK_GEO_XML(NN = NNsub, NC = NCsub, X = x, Y = y, Z = z)
        E_IO = VTK_CON_XML(NC = NCsub, connect = connect_vtk_io, offset = offset_vtk_io, cell_type = cell_type_vtk_io )
        E_IO = VTK_DAT_XML(var_location = 'cell', var_block_action = 'OPEN')
        E_IO = VTK_VAR_XML(NC_NN = NCsub           , varname = 'density'      , var = Rho % n(1:NCsub))
        E_IO = VTK_VAR_XML(NC_NN = NCsub           , varname = 'velocity'     , varX = U % n(1:NCsub)   , varY = V % n(1:NCsub)   , varZ = W % n(1:NCsub))
          !E_IO = VTK_VAR_XML(NC_NN = NCsub           , varname = 'velocity_mean', varX = U % mean(1:NCsub), varY = V % mean(1:NCsub), varZ = W % mean(1:NCsub))
        if (SIMULA /= DNS) then
          !E_IO = VTK_VAR_XML(NC_NN = NCsub, N_COL = 6, varname = 'stresses_mean', var  = stresses_2_vtk_io(:,:))
        end if
        E_IO = VTK_VAR_XML(NC_NN = NCsub           , varname = 'pressure'     , var  = P % n(1:NCsub))
        !E_IO = VTK_VAR_XML(NC_NN = NCsub           , varname = 'pressure_mean', var  = P % mean(1:NCsub))
        E_IO = VTK_VAR_XML(NC_NN = NCsub           , varname = 'PID'          , var  = PID(1:NCsub))
        E_IO = VTK_VAR_XML(NC_NN = NCsub           , varname = 'Q'            , var  = Q_2_vtk_io(1:NCsub))
        if (SIMULA == LES) then
            !E_IO = VTK_VAR_XML(NC_NN = NCsub           , varname = 'pp_mean'      , var  = pp_2_vtk_io(1:NCsub))
            E_IO = VTK_VAR_XML(NC_NN = NCsub           , varname = 'Q_LES_filt_x4', var  = Q_LES_filt_x4_2_vtk_io(1:NCsub))
        end if
        if (HOT == YES .and. MoinID .eq. 0) then
            E_IO = VTK_VAR_XML(NC_NN = NCsub           , varname = 'T'            , var  = T % n(1:NCsub))
          !E_IO = VTK_VAR_XML(NC_NN = NCsub           , varname = 'T_mean'       , var  = T % mean(1:NCsub))
          !E_IO = VTK_VAR_XML(NC_NN = NCsub           , varname = 'TT_mean'      , var  = TT_2_vtk_io(1:NCsub))
        end if

        E_IO = VTK_DAT_XML(var_location = 'cell', var_block_action = 'CLOSE')
        E_IO = VTK_GEO_XML()
        E_IO = VTK_END_XML()
    end if !SaveAreaIsEmpty

    ! in case of single processor
    if (sub_local == 0)  sub = 0
    ! ----------------------

    ! pvtu - xml - header
    ! header must be written by first processor with not empty area
    !
    call wait
    !write(*,*) 'PID =', this,' . ', SaveAreaIsEmpty

    if (FileExist) then ! if there are instructions to cut the area
        call MPI_ALLREDUCE                     &
            (SaveAreaIsEmpty,    & ! send buffer
            SaveAreaIsEmptySum,  & ! recv buffer
            Npro,                & ! length
            MPI_INTEGER,         & ! datatype
            MPI_SUM,             & ! operation
            MPI_COMM_WORLD,      &
            error_MPI)
    end if

    call wait
    !write(*,*) 'PID =', this,' . ', SaveAreaIsEmptySum

    if (this < 2) then
        name = nameIN
        call NamFil(0, nameTMP, '.pvtu', len_trim('.pvtu'))
        write(6, *) 'Now writing the file: ', nameTMP

        E_IO = PVTK_INI_XML(filename = nameTMP, mesh_topology = 'PUnstructuredGrid', tp='Float64')

        name = nameIN
        if (Npro > 0) then
            do c=1,Npro
                call NamFil(c, nameTMP, '.vtu', len_trim('.vtu'))
                if ( SaveAreaIsEmptySum(c) == 0 ) E_IO = PVTK_GEO_XML(source = nameTMP)
            end do
        else
            call NamFil(0, nameTMP, '.vtu', len_trim('.vtu'))
            E_IO = PVTK_GEO_XML(source = nameTMP)
        end if


        E_IO = PVTK_DAT_XML(var_location = 'cell', var_block_action = 'OPEN')

        E_IO = PVTK_VAR_XML(varname = 'density'      , tp='Float64')
        E_IO = PVTK_VAR_XML(Nc = 3, varname = 'velocity'     , tp='Float64' )
          !E_IO = PVTK_VAR_XML(Nc = 3, varname = 'velocity_mean', tp='Float64' )
        if (SIMULA /= DNS) then
          !E_IO = PVTK_VAR_XML(Nc = 6, varname = 'stresses_mean', tp='Float64' )
        end if
        E_IO = PVTK_VAR_XML(        varname = 'pressure'     , tp='Float64' )
        !E_IO = PVTK_VAR_XML(        varname = 'pressure_mean', tp='Float64' )
        E_IO = PVTK_VAR_XML(        varname = 'PID'          , tp='Int32'   )
        E_IO = PVTK_VAR_XML(        varname = 'Q'            , tp='Float64' )
        if (SIMULA == LES) then
            !E_IO = PVTK_VAR_XML(        varname = 'pp_mean'      , tp='Float64' )
            E_IO = PVTK_VAR_XML(        varname = 'Q_LES_filt_x4', tp='Float64' )
        end if
        if (HOT == YES .and. MoinID .eq. 0) then
            E_IO = PVTK_VAR_XML(        varname = 'T'            , tp='Float64' )
          !E_IO = PVTK_VAR_XML(        varname = 'T_mean'      , tp='Float64' )
          !E_IO = PVTK_VAR_XML(        varname = 'TT_mean'      , tp='Float64' )
        end if
        E_IO = PVTK_DAT_XML(var_location = 'cell', var_block_action = 'CLOSE')
        E_IO = PVTK_END_XML()

    end if ! if (this < 2) then

    call wait

    name = nameTMP2 ! return task name

    if ( allocated( connect_vtk_io         ) ) deallocate( connect_vtk_io         )
    if ( allocated( stresses_2_vtk_io      ) ) deallocate( stresses_2_vtk_io      )

    if ( allocated( Q_2_vtk_io             ) ) deallocate( Q_2_vtk_io             )
    if ( allocated( Q_LES_filt_x1_2_vtk_io ) ) deallocate( Q_LES_filt_x1_2_vtk_io )
    if ( allocated( Q_LES_filt_x2_2_vtk_io ) ) deallocate( Q_LES_filt_x2_2_vtk_io )
    if ( allocated( Q_LES_filt_x3_2_vtk_io ) ) deallocate( Q_LES_filt_x3_2_vtk_io )
    if ( allocated( Q_LES_filt_x4_2_vtk_io ) ) deallocate( Q_LES_filt_x4_2_vtk_io )

    if ( allocated( TT_2_vtk_io            ) ) deallocate( TT_2_vtk_io            )
    if ( allocated( pp_2_vtk_io            ) ) deallocate( pp_2_vtk_io            )

    if ( allocated( SaveAreaIsEmpty        ) ) deallocate( SaveAreaIsEmpty        )
    if ( allocated( SaveAreaIsEmptySum     ) ) deallocate( SaveAreaIsEmptySum     )

END SUBROUTINE SavParView_xml_bin_ascii_old
