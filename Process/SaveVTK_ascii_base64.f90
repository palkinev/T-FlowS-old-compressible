subroutine SaveVTK_ascii_base64(sub, NCsub, namAut)
  !----------------------------------------------------------------------!
  ! Reads: NAME.gmv and generates NAME.vti Paraview XML output file      !
  ! now supports ascii or binary format                                  !
  !------------------------------[Modules]-------------------------------!
  use all_mod
  use par_mod
  use allp_mod
  use les_mod
  use pro_mod
  use rans_mod
  use Base64_Mod

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

  !---- base64 mod tests
  integer(kind=4)              :: i_kind4(1)
  real(kind=8)                 :: r_kind8(1)
  integer(kind=1), allocatable :: pack(:)
  character(len=:), allocatable  :: code64
  i_kind4(1) = 0
  r_kind8(1) = 1

  call B64_Mod_Init()

  ! base 64 tests
  call pack_data_I4_R8(i_kind4, r_kind8, packed=pack)
  print *, pack(size(pack, dim=1))
  ! => 63 <<<
  print *, Get_Bit_Size_Real_8(1._8)
  ! => 64 <<<
  call b64_encode_I1_a(n=[120_1,-1_1], code=code64)
  print "(A)", code64
  !=> eP8= <<<
  call Encode_B64_R8_a(n=[1._8,2._8], code=code64)
  !=> AAAAAAAA8D8AAAAAAAAAQA== <<<
  print "(A)", code64
  

  print *, byte_size_int_1
  call Encode_B64_R8_a_xyz(x=[1._8,2._8], y=[1._8,2._8], z=[1._8,2._8], code=code64)
  print "(A)", code64


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


  !  <VTKFile type=”UnstructuredGrid” version=”0.1” byte_order=”LittleEndian”>
  !  <UnstructuredGrid>
  !  <!-- Actual data block  -->
  !  <Piece NumberOfPoints=”#” NumberOfCells=”#”>
  !
  !    <!-- Data defined in nodes  -->
  !    <PointData>
  !      ...
  !      <DataArray NumberOfComponents=”3” .../>
  !    </PointData>
  !
  !    <!-- Data defined in cell centers  -->
  !    <CellData>
  !      ...
  !      <DataArray NumberOfComponents=”3” .../>
  !    </CellData>
  !
  !    <!-- Geometry block. Can be stored once, then just copied  -->
  !    <Points>
  !      <DataArray type="Float32" NumberOfComponents="3" format="ascii">
  !    </Points>
  !    <Cells>
  !      <DataArray type=”Int32” Name=”connectivity” .../>
  !      <DataArray type=”Int32” Name=”offsets” .../>
  !      <DataArray type=”UInt8” Name=”types” .../>
  !    </Cells>
  !  </Piece>
  !  </UnstructuredGrid>
  !  </VTKFile>


  write(9,'(A)') '<?xml version="1.0"?>'
  write(9,'(A)') '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">'
  write(9,'(A)') '  <UnstructuredGrid>'
  write(9,*    ) '    <Piece NumberOfPoints="', NNsub,'" NumberOfCells="', NCsub, '" >'
  write(9,'(A)') '      <CellData>'

  !--vector: velocity
  write(9,'(A)') '        <DataArray type="Float64" Name="velocity" NumberOfComponents="3" format="binary">'
  call Encode_B64_R8_a_xyz(x=U % n(1: NC), y=V % n(1: NC), z=W % n(1: NC), code=code64)
  write(9,*) code64
  !deallocate(code64)
  write(9,'(A)') '        </DataArray>'


  write(9,'(A)') '      </CellData>'

  write(9,'(A)') '      <Points>'
  write(9,*) '         <DataArray type="Float64" NumberOfComponents="3" format="binary">'
  call Encode_B64_R8_a_xyz(x=x(1:NNsub), y=x(1:NNsub), z=x(1:NNsub), code=code64)
  write(9,*) code64
  !do c=1,NNsub
  !  write(9,*) x(c), y(c), z(c)
  !end do
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

end subroutine SaveVTK_ascii_base64
