subroutine Save_Mesh_CGNS_Lib(sub)                                             !
!   Writes in parallel 3-D unstructured grid and to CGNS grid file 'grid.cgns' !
!------------------------------------------------------------------------------!
!   All conventions for this lib are here:                                     !
!   https://cgns.github.io/CGNS_docs_current/sids/conv.html                    !
!------------------------------------------------------------------------------!
!   Parallel functions are describes here:                                     !
!   https://cgns.github.io/CGNS_docs_current/pcgns/library.html                !
!------------------------------------------------------------------------------!
!   this function inherited syntax mostly from:                                !
!   CGNS/src/Test_UserGuideCode/Fortran_code/write_grid_unst.F90               !
!------------------------------------------------------------------------------!
!   Array structures in current function are strictly followings:              !
!   Processor:    |        P_1        |               P_2               | ...  !
!   x,y,z:        |      (1 : NN_1)   |       NN_1 + 1 : NN_1 + NN_2    | ...  ! 
!   Connections:  |   (8, 1 : NC_1)   |   (8, NC_1 + 1 : NC_1 + NC_2)   | ...  ! 
!                                                                              !
!   Array:     ...| 1 + sum_(i=1)^(N-1) (NC_i) | 1 + sum_(i=1)^(N) (NC_i) |    ! 
!   Connectivity:     ...| 1 + sum_(i=1)^(N-1) (NC_i) | 1 + sum_(i=1)^(N) (NC_i) | 
!
!------------------------------------------------------------------------------!
!   TO DO:                                                                     !
!   parallel writing  - OK                                                     !
!   tetra cells                                                                !
!   partial geometry                                                           !
!----------------------------------[Modules]-----------------------------------!
  use all_mod
  use par_mod

  use cgns !on  Windows machines use "include 'cgnswin_f.h'"
!------------------------------------------------------------------------------!
  implicit none
  include 'mpif.h'
!-----------------------------------[Locals]-----------------------------------!
  integer(C_INT)     :: commsize, commrank

  ! CGNS variables
  integer            :: boundry_elems
  integer            :: file_idx ! CGNS file index number
  integer            :: section_idx ! Element section index, where 1 ≤ S ≤ nsections
  integer            :: base_idx ! Base index number, where 1 ≤ B ≤ nbases
  integer            :: zone_idx ! Zone index number, where 1 ≤ Z ≤ nzones
  integer            :: coord_idx ! Coordinate array index number, where 1 ≤ C ≤ ncoords
  integer            :: ier, iset

  integer(cgsize_t)  :: mesh_info(3)
  integer(cgsize_t),allocatable  :: cell_connections(:, :)
  integer(cgsize_t)  :: first_cell, last_cell
  integer(cgsize_t)  :: first_node, last_node
  character          :: base_name*32, zone_name*32

  ! Other variables
  integer             :: sub, sub_local

  character           :: dummy*100, nameTMP*80
  integer             :: n_nodes
  integer             :: c, n, i
  integer             :: nodes_in_current_cell ! 8 -> hex, 4-> tetra

  real, allocatable   :: x(:), y(:), z(:) ! cell node coordinates

  ! sub zone variables
  integer                :: NC_In_Sub_Zone
  REAL                   :: x1, x2, x3, x4
  REAL                   :: y1, y2, y3, y4
  REAL                   :: z1, z2
  logical                :: FileExist
  integer*1, allocatable :: SaveAreaIsEmpty(:)
  integer*1, allocatable :: SaveAreaIsEmptySum(:)
  integer*4              :: Array_at_root(1:Npro)
!---------------------------------[Interface]----------------------------------!
!==============================================================================!
  !-------------------------------------------------------------------!
  !   Parallel commands
  call cgp_mpi_comm_f(mpi_comm_world, ier) ! Set the MPI communicator for CGNS
  !call cgp_mpi_info_f(integer info, ier) ! Set the MPI info object for CGNS
  call cgp_mpi_info_f(0, ier) ! Set the MPI info object for CGNS
  call cgp_pio_mode_f(CGP_INDEPENDENT, ier) ! Set the parallel IO mode for CGNS
  !call cgp_pio_mode_f(CGP_COLLECTIVE, ier) ! Set the parallel IO mode for CGNS
  !-------------------------------------------------------------------!

  if (this < 2) write(6,'('' Program write_grid_unst'')')
  if (CG_BUILD_64BIT) then
    if (this < 2) write(6,'('' ...using 64-bit mode for particular integers'')')
  end if



  !----------------------------------------------------------------------!
  !   feature:                                                           !
  !   draw only zones, which are inside area specified in Paraview.dat   !
  !   if file does not not exists -> draw full zone as usual             !
  !----------------------------------------------------------------------!
    
!  sub_local = 1
!  if (sub .ne. 0) then ! multi processor case
!    allocate(SaveAreaIsEmpty(1:Npro));
!    allocate(SaveAreaIsEmptySum(1:Npro));
!  else                 ! one processor case
!    sub_local = 0
!    sub = 1
!    allocate(SaveAreaIsEmpty(1:1));
!    allocate(SaveAreaIsEmptySum(1:1));
!  end if
!
!  SaveAreaIsEmpty = 0    ! 0 -> not empty
!  SaveAreaIsEmptySum = 0 ! 0 -> not empty
!  NC_In_Sub_Zone = 0
!
!  inquire (file="ParaView.dat",exist=FileExist)
!  if (FileExist) then
!
!    open(119, FILE='ParaView.dat', status="old", action="read")
!    call wait
!
!    !     x3,y3           x2,y2             z1         z2
!    !       |---------------\                |---------|
!    !       |      ^ y       \               |         |
!    !       |      |          \              |         | ---> z-plane
!    !       |      ---> x      \             |         |
!    !       |-------------------|            |---------|
!    !     x1,y1               x4,y4         z1         z2
!
!    if(this < 2) write(*,*) '# Now reading the file: Paraview.dat '
!    call ReadC(119,inp,tn,ts,te)
!    read(inp(ts(1):te(1)),*) x1
!    read(inp(ts(2):te(2)),*) x2
!    x3 = x1
!    x4 = x2
!
!    if(tn==4)   read(inp(ts(3):te(3)),*) x3
!    if(tn==4)   read(inp(ts(4):te(4)),*) x4
!
!    call ReadC(119,inp,tn,ts,te)
!    read(inp(ts(1):te(1)),*) y1
!    read(inp(ts(2):te(2)),*) y2
!    y4 = y1
!    y3 = y2
!    if(tn==4)   read(inp(ts(3):te(3)),*) y3
!    if(tn==4)   read(inp(ts(4):te(4)),*) y4
!
!    call ReadC(119,inp,tn,ts,te)
!    read(inp(ts(1):te(1)),*) z1
!    read(inp(ts(2):te(2)),*) z2
!
!    if(this < 2) write(*,*) '#  x:(', x1, ';', x2, ';', x3, ';', x4, ')'
!    if(this < 2) write(*,*) '#  x:(', y1, ';', y2, ';', y3, ';', y4, ')'
!    if(this < 2) write(*,*) '#  z:(', z1, ';', z2, ')'
!
!    call wait
!
!    ! filtering cells
!    do c = 1, NC
!      if (    xc(c) > x1 &
!        .and. xc(c) < x2 &
!        .and. xc(c) > x3 &
!        .and. xc(c) < x4 &
!        .and. yc(c) > y1 + (y4-y1)*(xc(c)-x1)/(x4-x1) &
!        .and. yc(c) < y2 + (y3-y2)*(xc(c)-x2)/(x3-x2) &
!        .and. zc(c) > z1 &
!        .and. zc(c) < z2 ) then
!
!        NC_In_Sub_Zone = NC_In_Sub_Zone + 1
!
!      end if
!    end do
!
!    if (NC_In_Sub_Zone == 0) then
!      SaveAreaIsEmpty(sub) = 1
!    end if
!
!  end if ! else draw full area
!
!  if (FileExist) then ! if there are instructions to cut the area
!    call MPI_ALLREDUCE(      &
!        SaveAreaIsEmpty,     & ! send buffer
!        SaveAreaIsEmptySum,  & ! recv buffer
!        Npro,                & ! length
!        MPI_INTEGER,         & ! datatype
!        MPI_SUM,             & ! operation
!        MPI_COMM_WORLD,      &
!        ier )
!    end if
!
!  ! in case of single processor
!  if (sub_local == 0)  sub = 0

  !-------------------------------------!
  ! Read .gmv block for additional data !
  !-------------------------------------!

  call NamFil(sub, nameTMP, '.gmv', len_trim('.gmv'))
  open(9, FILE = nameTMP)
  if (this <2 ) write(*,*) 'Now reading the file: ', nameTMP

  call ReadC(9,inp,tn,ts,te)  ! read line "gmvinput ascii"
  call ReadC(9,inp,tn,ts,te)  ! read "nodes  $number_of_nodes"
  read(inp(ts(1):te(1)),*) dummy ! "nodes" word
  read(inp(ts(2):te(2)),*) n_nodes ! nodes count


  !---------- x, y, z arrays structure
  ! |      (1 : NN_1)   |       NN_1 + 1 : NN_1 + NN_2    | ... 
    Array_at_root = 0
    write(*,*) "PID = ", this, " Nodes = ", n_nodes
    call wait
    if(this < 2 ) write(*,*) "----------------------------"

  call wait
  call mpi_gather(  &
    n_nodes,         & ! send number of nodes, which is
    1,              & ! 1 element
    MPI_INTEGER4,   & ! of 32-bit integer type
    Array_at_root,  & ! to array "Nodes"
    1,              & ! which is 1 per processor long.
    MPI_INTEGER4,   & ! Recieved data type is 32-bit integer
    0,              & ! at the proc 1
    mpi_comm_world, & ! communicator
    ier)              ! mpi_error
    
    write(*,*) "PID = ", this, " Nodes at root after gather = ", Array_at_root
    
    if(this < 2 ) write(*,*) "----------------------------"

    if (this < 2) then

      Array_at_root(1) = 1
      do i = 2, Npro
        Array_at_root(i) = Array_at_root(i) + Array_at_root(i-1)
      end do
      if(this < 2 ) write(*,*) "PID = ", this, " Target_array = ", Array_at_root
      if(this < 2 ) write(*,*) "----------------------------"

    end if

  call wait
  call mpi_scatter( &
    Array_at_root,  & ! send number of nodes, which is
    1,              & ! 1 elements per processor
    MPI_INTEGER4,   & ! of 32-bit integer type
    first_node,     & ! to integer "first_node"
    Npro,           & ! which is Npro - 1 elements long.
    MPI_INTEGER4,   & ! Recieved data type is 32-bit integer
    0,              & ! from the proc 1
    mpi_comm_world, & ! communicator
    ier)              ! mpi_error
    write(*,*) "PID = ", this, " fist node_after_scatter = ", first_node

  !---------- read x, y, z of cell nodes consequently
  allocate(x(first_node:first_node+n_nodes-1), stat = ier); x = 0.
    if (ier.ne.0)then
       print*, '*FAILED* allocation of x'
       call cgp_error_exit_f()
    endif
  allocate(y(first_node:first_node+n_nodes-1), stat = ier); y = 0.
    if (ier.ne.0)then
       print*, '*FAILED* allocation of y'
       call cgp_error_exit_f()
    endif
  allocate(z(first_node:first_node+n_nodes-1), stat = ier); z = 0.
    if (ier.ne.0)then
       print*, '*FAILED* allocation of z'
       call cgp_error_exit_f()
    endif

  !write(*,*) "pid=", this, "NC=", NC, "NN=", n_nodes

  do n = lbound(x,dim=1), ubound(x,dim=1)
    call ReadC(9,inp,tn,ts,te)
    read(inp(ts(1):te(1)),*) x(n)
  end do
  do n = lbound(y,dim=1), ubound(y,dim=1)
    call ReadC(9,inp,tn,ts,te)
    read(inp(ts(1):te(1)),*) y(n)
  end do
  do n = lbound(z,dim=1), ubound(z,dim=1)
    call ReadC(9,inp,tn,ts,te)
    read(inp(ts(1):te(1)),*) z(n)
  end do

  call ReadC(9,inp,tn,ts,te) !  read line "cells  $number_o_cells"
  read(inp(ts(1):te(1)),*) dummy ! "cells" word
  read(inp(ts(2):te(2)),*) NC_In_Sub_Zone ! cells count

  if (NC_In_Sub_Zone .ne. NC) then
    write(*,*) 'number of cells read and processed is different, exiting!'
    stop
  end if

  call wait
  !---------- cell_connections array structure
  !   Connections:  |   (8, 1 : NC_1)   |   (8, NC_1 + 1 : NC_1 + NC_2)   | ...  ! 
    Array_at_root = 0
!    write(*,*) "PID = ", this, " NC = ", NC
    call wait
    if(this < 2 ) write(*,*) "----------------------------"

  call wait
  call mpi_gather(  &
    NC,             & ! send number of nodes, which is
    1,              & ! 1 element
    MPI_INTEGER4,   & ! of 32-bit integer type
    Array_at_root,  & ! to array "Nodes"
    1,              & ! which is 1 per processor long.
    MPI_INTEGER4,   & ! Recieved data type is 32-bit integer
    0,              & ! at the proc 1
    mpi_comm_world, & ! communicator
    ier)              ! mpi_error

    write(*,*) "PID = ", this
    
!    write(*,*) "PID = ", this, " Cells at root after gather = ", Array_at_root
!    
!    if(this < 2 ) write(*,*) "----------------------------"


    stop
  !---------- read nodes connection for cells
  ! gmv syntax ex:
  ! hex 8
  ! 1        2       67       66      261      262      327      326
  ! 

  if (this == 1) then
    allocate(cell_connections(8,1:288)); cell_connections = 0.
  end if

  if (this == 2) then
    allocate(cell_connections(8,289:576)); cell_connections = 0.
  end if

  do c = lbound(cell_connections,dim=2),  ubound(cell_connections,dim=2)
    nodes_in_current_cell = 0

    call ReadC(9,inp,tn,ts,te)
    read(inp(ts(1):te(1)),*) dummy ! "hex"/"tetra"/"pyramid"/ word
    read(inp(ts(2):te(2)),*) nodes_in_current_cell ! nodes count to make a cell

    call ReadC(9,inp,tn,ts,te)
    do n = 1, nodes_in_current_cell
      read(inp(ts(n):te(n)),*) cell_connections(n, c) ! nodes index is the same as in gmv
    end do
  end do

  close(9) ! close .gmv file

  if (this == 2) then
    cell_connections(:,:) = cell_connections(:,:) + 684
    !write(*,"(A,8I7)") "cc=", cell_connections(:,:)
  end if
  call wait
  !-------------------------!
  ! CGNS mesh writing block !
  !-------------------------!
  if (this < 2) write(6,'('' created simple 3-D grid points'')')

  !----------   open CGNS file for write
  
  call cgp_open_f('grid.cgns',CG_MODE_WRITE,file_idx,ier)
    if (ier.ne.CG_OK) then
       print*,'*FAILED* cgp_open_f'
       call cgp_error_exit_f()
    endif

  !----------   write x, y, z grid points to CGNS file (use SIDS-standard)
  !if (this == 1) base_name = "Base 1" !define zone name (user can give any name)
  !if (this == 2) base_name = "Base 2" !define zone name (user can give any name)

  base_name = "Base 1"
  call cg_base_write_f(file_idx,base_name,3,3,base_idx,ier)
    if (ier.ne.CG_OK) then
       print*,'*FAILED* cg_base_write_f'
       call cgp_error_exit_f()
    endif

  zone_name = "Zone 1" !define zone name (user can give any name)

  !if (this == 1) zone_name = "Zone 1" !define zone name (user can give any name)
  !if (this == 2) zone_name = "Zone 2" !define zone name (user can give any name)


  !---------- sum of nodes through processors
  i = n_nodes
  call IglSum(i)
  !---------- total cell count
  mesh_info(1) = i
  !---------- sum of cells through processors
  i = NC_In_Sub_Zone
  call IglSum(i)
  mesh_info(2) = i
  mesh_info(3) = 0 ! boundary vertex size (zero if elements not sorted)

  call cg_zone_write_f(file_idx,base_idx,zone_name,mesh_info, &
    Unstructured,zone_idx,ier) ! create zone
    if (ier.ne.CG_OK) then
     print*,'*FAILED* cg_zone_write_f'
     call cgp_error_exit_f()
  endif

  first_node = lbound(x,dim=1)
  last_node = ubound(x,dim=1)
  first_cell = lbound(cell_connections,dim=2) ! first cell
  last_cell = ubound(cell_connections,dim=2) ! last cell

  !---------- create empty "CoordinateX" node in DB
  call cgp_coord_write_f(file_idx,base_idx,zone_idx,RealDouble, &
    'CoordinateX',coord_idx,ier)
    if (ier.ne.CG_OK) then
      print*,'*FAILED* cgp_coord_write_f (Coord_X)'
      call cgp_error_exit_f()
    endif

  !---------- fill "CoordinateX" node in DB
  call cgp_coord_write_data_f(file_idx,base_idx,zone_idx,coord_idx, &
    first_node,last_node,x,ier)
    if (ier.ne.CG_OK) then
      print*,'*FAILED* cgp_coord_write_data_f (Coord_X)'
      call cgp_error_exit_f()
    endif

  !---------- create empty "CoordinateY" node in DB
  call cgp_coord_write_f(file_idx,base_idx,zone_idx,RealDouble, &
    'CoordinateY',coord_idx,ier)
    if (ier.ne.CG_OK) then
      print*,'*FAILED* cgp_coord_write_f (Coord_Y)'
      call cgp_error_exit_f()
    endif

    !---------- fill "CoordinateY" node in DB
  call cgp_coord_write_data_f(file_idx,base_idx,zone_idx,coord_idx, &
      first_node,last_node,y,ier)
    if (ier.ne.CG_OK) then
      print*,'*FAILED* cgp_coord_write_data_f (Coord_Y)'
      call cgp_error_exit_f()
    endif

  !---------- create empty "CoordinateZ" node in DB
  call cgp_coord_write_f(file_idx,base_idx,zone_idx,RealDouble, &
    'CoordinateZ',coord_idx,ier)
    if (ier.ne.CG_OK) then
      print*,'*FAILED* cgp_coord_write_f (Coord_Z)'
      call cgp_error_exit_f()
    endif

  !---------- fill "CoordinateZ" node in DB
  call cgp_coord_write_data_f(file_idx,base_idx,zone_idx,coord_idx, &
    first_node,last_node,z,ier) ! write z grid coordinates
    if (ier.ne.CG_OK) then
      print*,'*FAILED* cgp_coord_write_data_f (Coord_Z)'
      call cgp_error_exit_f()
    endif

  boundry_elems = 0 ! unsorted boundary elements

  !---------- create empty "Elem" node in DB with nodes connectivity
  call cgp_section_write_f(file_idx,base_idx,zone_idx, &
    'Elem',HEXA_8,int(1,8),int(i,8),boundry_elems, &
    section_idx,ier) ! prepare to write connectivity data
    if (ier.ne.CG_OK) then
      print*,'*FAILED* cgp_section_write_f'
      call cgp_error_exit_f()
    endif
  !---------- fill "Elem" node in DB with nodes connectivity
  call cgp_elements_write_data_f(file_idx,base_idx,zone_idx, &
    section_idx,first_cell,last_cell,cell_connections,ier) ! write HEXA_8 element connectivity
      if (ier.ne.CG_OK) then
      print*,'*FAILED* cgp_elements_write_data_f (elements)'
      call cgp_error_exit_f()
    endif

  !---------- close DB
  call cgp_close_f(file_idx,ier)
  if (this < 2) write(6,'('' Successfully wrote unstructured grid to file'', &
    '' grid.cgns'')') ! close CGNS file

end subroutine Save_Mesh_CGNS_Lib

!----------------- sequential

! call cg_coord_write_f(file_idx,base_idx,zone_idx,RealDouble, &
!     'CoordinateX',x,coord_idx,ier) ! write x grid coordinates 
! call cg_coord_write_f(file_idx,base_idx,zone_idx,RealDouble, &
!     'CoordinateY',y,coord_idx,ier) ! write y grid coordinates 
! call cg_coord_write_f(file_idx,base_idx,zone_idx,RealDouble, &
!     'CoordinateZ',z,coord_idx,ier) ! write z grid coordinates 

!  boundry_elems = 0 ! unsorted boundary elements
!  call cg_section_write_f(file_idx,base_idx,zone_idx, &
!       'Elem',HEXA_8,first_cell,last_cell,boundry_elems,cell_connections, &
!       section_idx,ier) ! write HEXA_8 element connectivity
