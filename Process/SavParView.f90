!======================================================================!
subroutine SavParView(sub, NCsub, namAut)
    !----------------------------------------------------------------------!
    ! Reads: NAME.gmv and generates NAME.vti Paraview XML output file      !
    ! ~~~~~~~                                                              !
    !------------------------------[Modules]-------------------------------!
    use all_mod
    use par_mod
    use allp_mod
    use les_mod
    use pro_mod
    use rans_mod
    !----------------------------------------------------------------------!
    implicit none
    
    !-----------------------------[Parameters]-----------------------------!
    integer ::  NCsub, sub
    character(len=80) :: storename, namTem, namXML
    character(len=10) :: namAut
    !-------------------------------[Locals]-------------------------------!
    integer   :: c,  c1,  c2,  n, s
    character(len=80) :: namOut, nameIn
    character(len=100):: stringadummy
    real,allocatable :: x(:), y(:), z(:)    ! self evident
    integer,allocatable :: connessione(:,:) ! connection
    integer :: celleconnessione
    integer :: NNsub, NCsub_new
    integer :: off_set_connection
    integer :: i
    integer :: numprocessi
    real    :: Unor
    real    :: p_sum, pressure_err(1:NC),vol


    !======================================================================!

    !  write(*,*) 'writing paraview XML data file: ', namAut

    !  if (this < 2) then
    !      write(*,*) 'writing paraview XML data file: ', namAut
    !  end if

    pressure_err(:) = 0.0

    namTem = name
    storename = namAut
    !<<<<<<<<<<<<<<<<<<<<<<<<<!
    !                         !
    !     reads GMV file      !
    !                         !
    !<<<<<<<<<<<<<<<<<<<<<<<<<!
    call NamFil(sub, namOut, '.gmv', len_trim('.gmv'))
    open(9, FILE=namOut)
    if (this <2) then
        write(*,*) 'Now reading the file: ', namOut
    end if

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
            connessione(n,1) = connessione (n-1,1)+ off_set_connection
        end if
    
        celleconnessione = celleconnessione + off_set_connection

        call ReadC(9,inp,tn,ts,te)
        do c=1,off_set_connection
            read(inp(ts(c):te(c)),*) connessione(n,c+1)
        end do
    end do

    close(9)

    name = namAut
    call NamFil(sub, namXML, '.vtu', len_trim('.vtu'))

    open(9, FILE=namXML)
    if (this <2) then
        write(6, *) 'Now writing the file:', namXML
    end if


    write(9,'(A)') '<?xml version="1.0"?>'
    write(9,'(A)') '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">'
    write(9,'(A)') '  <UnstructuredGrid>'
    write(9,*    ) '    <Piece NumberOfPoints="', NNsub,'" NumberOfCells="', NCsub, '" >'
    write(9,'(A)') '      <CellData>'

    !--vector: velocity
    write(9,'(A)') '        <DataArray type="Float32" Name="velocity" NumberOfComponents="3" format="ascii">'
    do c=1,NCsub
        write(9,*) U % n(c), V % n(c), W % n(c)
    end do
    write(9,'(A)') '        </DataArray>'


    call GraPhi(U % n , 1, Ux ,.TRUE.)    ! dU/dx
    call GraPhi(U % n , 2, Uy ,.TRUE.)    ! dU/dy
    call GraPhi(U % n , 3, Uz ,.TRUE.)    ! dU/dz
  
    write(9,'(A)') '        <DataArray type="Float32" Name="dv_dx_i" NumberOfComponents="3" format="ascii">'
    do c=1,NCsub
        write(9,*) Ux(c), Uy(c), Uz(c)
    end do
    write(9,'(A)') '        </DataArray>'


    !--Processor ID
    write(9,'(A)') '        <DataArray type="Float32" Name="ProcessorID" format="ascii">'
    do c=1,NCsub
        write(9,*) sub
    end do
    write(9,'(A)') '        </DataArray>'

    !--scalar: pressure
    write(9,'(A)') '        <DataArray type="Float32" Name="pressure" format="ascii">'
    do c=1,NCsub
        write(9,*) P % n(c)
    end do
    write(9,'(A)') '        </DataArray>'

    call GradP (P % n,Px,Py,Pz) ! dP/dx, dP/dy, dP/dz
    !call GradP3(P % n,Px,Py,Pz) ! dP/dx, dP/dy, dP/dz

    !--scalar: pressure
    write(9,'(A)') '        <DataArray type="Float32" Name="dp_dx_i" NumberOfComponents="3" format="ascii">'
    do c=1,NCsub
        write(9,*) Px(c),Py(c),Pz(c)
    end do
    write(9,'(A)') '        </DataArray>'


    write(9,'(A)') '        <DataArray type="Float32" Name="pres.correction" format="ascii">'
    do c=1,NCsub
        write(9,*) PP % n(c)
    end do
    write(9,'(A)') '        </DataArray>'


    p_sum  = 0.0
    vol  = 0.0

    do c=1,NC
        if ( c .ne. 0) then
            p_sum = p_sum + P % n(c)*volume(c)
            vol  = vol  + volume(c)
        end if
    end do

    call GloSum(vol)
    call GloSum(p_sum)

    do c=1,NC
        !--------------L^2 error
        pressure_err(c) = pressure_err(c) + volume(c)*(P % n(c) - p_sum/vol )**2
    end do


    write(9,'(A)') '        <DataArray type="Float32" Name="pressure_err" format="ascii">'
    do c=1,NCsub
        write(9,*) pressure_err(c)
    end do
    write(9,'(A)') '        </DataArray>'

    !--vector: velocity
    write(9,'(A)') '        <DataArray type="Float32" Name="density" NumberOfComponents="3" format="ascii">'
    do c=1,NCsub
        write(9,*) Rho % n(c), Rho % o(c), Rho % oo(c)
    end do
    write(9,'(A)') '        </DataArray>'

    !------------- Calculate the max mass error (INCOMMENT ORIGINAL IN CorUVW.f90)
    b(:) = 0.0

    do s=1,NS
        c1=SideC(1,s)
        c2=SideC(2,s)
        if(c2  > 0 .or. c2  < 0 .and. TypeBC(c2) == BUFFER) then
            b(c1) = b(c1)-Flux(s)
            if(c2  > 0) b(c2) = b(c2)+Flux(s)
        else
            b(c1) = b(c1)-Flux(s)
        end if
    end do

    !--scalar: wall distance
    write(9,'(A)') '        <DataArray type="Float32" Name="flux" format="ascii">'
    do c=1,NCsub
        write(9,*) - b(c)
    end do
    write(9,'(A)') '        </DataArray>'

    Unor = 0.0

    do c=1,NC     ! - d (rho^n+1) /dt
        ! 1)
        b(c) = - b(c) + volume(c)/dt*(1.5*Rho % n(c) - 2.0*Rho % o(c) + 0.5*Rho % oo(c) )
      ! 2)
      !b(c) = b(c) - volume(c)/dt*(1.5*Rho % n(c) - 2.0*Rho % o(c) )
      !b(c) = b(c) + volume(c)/dt*( 0.5*Rho % oo(c) )
      ! 3)
      !Unor = (1.5*Rho % n(c) - 2.0*Rho % o(c) + 0.5*Rho % oo(c) )
      !b(c) = b(c) - volume(c)/dt*Unor
      ! 4)
      !b(c) = b(c) - (volume(c)/dt)*(1.5*Rho % n(c) - 2.0*Rho % o(c) + 0.5*Rho % oo(c) )
    end do


    write(9,'(A)') '        <DataArray type="Float32" Name="mass_err" format="ascii">'
    do c=1,NCsub
        write(9,*) b(c) / volume(c)
    end do
    write(9,'(A)') '        </DataArray>'



    !------------- Calculate the max mass error  - END
  
    !--scalar: wall distance
    write(9,'(A65)') '<DataArray type="Float32" Name="wall distance" format="ascii">'
    do c=1,NCsub
        if(IsNearWall(c)) then
            write(9,*) sqrt(WallDs(c)*sqrt(U%n(c)**2+V%n(c)**2+W%n(c)**2)/VISc)
        else
            write(9,*) 0.0
        end if
    end do
    write(9,'(A20)') '        </DataArray>'



    !------------------------------------------------------------------------------------------------------------
    ! T H E R M A L  F I E L D  R E L A T E D   V A R I A B L E S (HOT == YES)
    !------------------------------------------------------------------------------------------------------------

    if(HOT == YES) then

        ! scalar: temperature
        write(9,'(A)') '        <DataArray type="Float32" Name="Zmix" format="ascii">'
        do c=1,NCsub
            write(9,*) Zmix % n(c)
        end do
        write(9,'(A)') '        </DataArray>'

    end if !HOT == YES




    write(9,'(A)') '      </CellData>'
  
    write(9,'(A)') '      <Points>'
    write(9,*) '         <DataArray type="Float32" NumberOfComponents="3" format="ascii">'
    do c=1,NNsub
        write(9,*) x(c), y(c), z(c)
    end do
    write(9,'(A)') '        </DataArray>'
    write(9,'(A)') '      </Points>'
    write(9,'(A)') '      <Cells>'
    write(9,'(A)') '        <DataArray type="Int32" Name="connectivity" format="ascii">'
    do n=1,NCsub
        !hexa
        if (connessione(n,9) /= 0) then
            do i=2,9
                write(9,*) connessione(n,i)-1
            end do
        else
            !prism
            if ((connessione(n,8) == 0).and.(connessione(n,7) /= 0)) then
                write(9,*) connessione(n,2)-1, connessione(n,4)-1, &
                    connessione(n,3)-1, connessione(n,5)-1, connessione(n,7)-1, connessione(n,6)-1
            else
                if ((connessione(n,7) == 0).and.(connessione(n,6) /= 0)) then
                    write(9,*) connessione(n,5)-1, connessione(n,4)-1,&
                        connessione(n,3)-1, connessione(n,6)-1, connessione(n,2)-1
                else
                    write(9,*) connessione(n,5)-1, connessione(n,4)-1,&
                        connessione(n,3)-1, connessione(n,2)-1
                end if
            end if
        end if

    end do
    write(9,'(A)') '        </DataArray>'
    write(9,'(A)') '        <DataArray type="Int32" Name="offsets" format="ascii">'
    do n=1,NCsub
        write(9,*) connessione(n,1)
    end do
    write(9,'(A)') '        </DataArray>'
    write(9,'(A)') '        <DataArray type="UInt8" Name="types" format="ascii">'
    do n=1,NCsub
        if ((connessione(n,6)==0)) then !thetra
            write(9,*) '10'
        else
            if ((connessione(n,7)==0)) then
                write(9,*) '14'          !pyr
            else
                if ((connessione(n,8)==0)) then
                    write(9,*) '13'      !prism
                else
                    write(9,*) '12'      !hexa
                end if
            end if
        end if
    end do
    write(9,'(A)') '        </DataArray>'
    write(9,'(A)') '      </Cells>'
    write(9,'(A)') '</Piece>'
    write(9,'(A)') '</UnstructuredGrid>'
    write(9,'(A)') '</VTKFile>'

    ! write down the master file for parallel jobs
    ! actually it writes it down also for sequential process but it's not
    ! necessary for paraview
    numprocessi = Npro
    call GloMax(real(numprocessi))
    if (numprocessi==0) then
        numprocessi=1
    end if

    if (this < 2) then
        name = storename
        call NamFil(0, name, '.pvtu', len_trim('.pvtu'))
        open(12, FILE=name)
        !     write(*,*) 'writing parallel master file: ', name
        write(12,'(A)') '<?xml version="1.0"?>'
        write(12,'(A)') '<VTKFile type="PUnstructuredGrid" version="0.1" byte_order="LittleEndian">'
        write(12,'(A)') '  <PUnstructuredGrid GhostLevel="0">'
        write(12,'(A)') '    <PCellData>'

        write(12,'(A)') '      <PDataArray type="Float32" Name="velocity" NumberOfComponents="3"/>'
        write(12,'(A)') '      <PDataArray type="Float32" Name="dv_dx_i" NumberOfComponents="3"/>'
        write(12,'(A)') '      <PDataArray type="Float32" Name="ProcessorID"/>'
        write(12,'(A)') '      <PDataArray type="Float32" Name="pressure"/>'
        write(12,'(A)') '      <PDataArray type="Float32" Name="dp_dx_i"  NumberOfComponents="3"/>'
        write(12,'(A)') '      <PDataArray type="Float32" Name="pres.correction"/>'
        write(12,'(A)') '      <PDataArray type="Float32" Name="pressure_err"/>'
        write(12,'(A)') '      <PDataArray type="Float32" Name="density"  NumberOfComponents="3"/>'
        write(12,'(A)') '      <PDataArray type="Float32" Name="flux"/>'
        write(12,'(A)') '      <PDataArray type="Float32" Name="mass_err"/>'
        write(12,'(A)') '      <PDataArray type="Float32" Name="wall distance"/>'
        write(12,'(A)') '      <PDataArray type="Float32" Name="Zmix"/>'
     
        write(12,'(A)') '    </PCellData>'
        write(12,'(A)') '    <PPoints>'
        write(12,'(A)') '      <PDataArray type="Float32" NumberOfComponents="3"/>' !x,y,z,
        write(12,'(A)') '    </PPoints>'
        ! qui ci devo mettere la stampa dei vari files di sottoprocesso  <Piece Source="chan1_180.vtu"/>
        name = namAut
        do i=1,numprocessi
            call NamFil(i, nameIn, '.vtu"/>', len_trim('.vtu"/>'))
            !           write(12,*) '<Piece Source="',nameIn, '"/>'
            write(12,'(A)',advance="no") '    <Piece Source="'
            write(12,'(A)') nameIn
        end do
     
        write(12,*) '  </PUnstructuredGrid>'
        write(12,*) '</VTKFile>'
        close(12)
    end if !(this <2)

    close(9)
    name = namTem

    call wait

    deallocate(x)
    deallocate(y)
    deallocate(z)
    deallocate(connessione)

    return

end subroutine SavParView
