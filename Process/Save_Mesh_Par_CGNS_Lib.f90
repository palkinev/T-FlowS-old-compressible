subroutine Save_Mesh_Par_CGNS_Lib(sub)                                         !
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
!------------------------------------------------------------------------------!
!   TO DO:                                                                     !
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
  ! CGNS variables
  integer :: boundry_elems
  integer :: file_idx ! CGNS file index number
  integer :: section_idx ! Element section index, where 1 ≤ S ≤ nsections
  integer :: base_idx ! Base index number, where 1 ≤ B ≤ nbases
  integer :: zone_idx ! Zone index number, where 1 ≤ Z ≤ nzones
  integer :: coord_idx ! Coordinate array index number, where 1 ≤ C ≤ ncoords
  integer :: ier

  integer(cgsize_t)              :: mesh_info(3)
  integer(cgsize_t), allocatable :: cell_connections(:, :)
  character                      :: base_name*32, zone_name*32
  integer(cgsize_t)              :: first_cell_cgns, last_cell_cgns
  integer(cgsize_t)              :: first_node_cgns, last_node_cgns

  ! Other variables
  integer*8           :: first_cell
  integer*8           :: first_node
  integer             :: sub

  character           :: dummy*100, nameTMP*80
  integer*8           :: n_nodes
  integer             :: c, n, i
  integer             :: nodes_in_current_cell ! 8 -> hex, 4-> tetra
  real, allocatable   :: x(:), y(:), z(:) ! cell node coordinates

  ! sub zone variables
  integer*8           :: NC_In_Sub_Zone
!---------------------------------[Interface]----------------------------------!
!==============================================================================!
  !---------- cgns parallel mode
  call cgp_pio_mode_f(CGP_INDEPENDENT, ier) ! Set the parallel IO mode for CGNS

  if (this < 2) write(6,'('' Program write_grid_unst'')')
  if (CG_BUILD_64BIT) then
    if (this < 2) write(6,'('' ...using 64-bit mode for particular integers'')')
  end if

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

  !---------- fetch x, y, z arrays structure
  call Fetch_Arrays_Dimensions_Par_CGNS_Lib(first_node, n_nodes)

  !---------- allocate x, y, z
  allocate(x(first_node: first_node + n_nodes - 1), stat = ier); x = 0.
    if (ier .ne. 0) then
       print*, '*FAILED* allocation of x'
       call cgp_error_exit_f()
    endif
  allocate(y(first_node: first_node + n_nodes - 1), stat = ier); y = 0.
    if (ier .ne. 0) then
       print*, '*FAILED* allocation of y'
       call cgp_error_exit_f()
    endif
  allocate(z(first_node: first_node + n_nodes - 1), stat = ier); z = 0.
    if (ier .ne. 0) then
       print*, '*FAILED* allocation of z'
       call cgp_error_exit_f()
    endif

  !---------- read x, y, z of cell nodes consequently
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

  !---------- fetch cell_connections array structure
  call Fetch_Arrays_Dimensions_Par_CGNS_Lib(first_cell, NC_In_Sub_Zone)

  !---------- read nodes connection for cells
  ! gmv syntax ex:
  ! hex 8
  ! 1        2       67       66      261      262      327      326
  !

  !---------- allocate cell_connections
  allocate(cell_connections(8,first_cell: first_cell + NC_In_Sub_Zone - 1), stat = ier)
  cell_connections = 0.
    if (ier.ne.0)then
       print*, '*FAILED* allocation of cell_connections'
       call cgp_error_exit_f()
    endif

  !---------- read cell_connections
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

  call wait
  close(9) ! close .gmv file

  ! shift nodes id according to a number of nodes
  cell_connections = cell_connections + first_node - 1
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

  !----------   create a base node
  base_name = "Base 1"
  call cg_base_write_f(file_idx,base_name,3,3,base_idx,ier)
    if (ier.ne.CG_OK) then
       print*,'*FAILED* cg_base_write_f'
       call cgp_error_exit_f()
    endif

  zone_name = "Zone 1" !define zone name (user can give any name)

  !---------- sum of nodes through processors
  i = n_nodes
  call IglSum(i)
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

  first_node_cgns = lbound(x,dim=1)
  last_node_cgns  = ubound(x,dim=1)

  call wait

  !---------- create empty "CoordinateX" node in DB
  call cgp_coord_write_f(file_idx,base_idx,zone_idx,RealDouble, &
    'CoordinateX',coord_idx,ier)
    if (ier.ne.CG_OK) then
      print*,'*FAILED* cgp_coord_write_f (Coord_X)'
      call cgp_error_exit_f()
    endif

  !---------- fill "CoordinateX" node in DB
  call cgp_coord_write_data_f(file_idx,base_idx,zone_idx,coord_idx, &
    first_node_cgns,last_node_cgns,x,ier)
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
      first_node_cgns,last_node_cgns,y,ier)
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
    first_node_cgns,last_node_cgns,z,ier) ! write z grid coordinates
    if (ier.ne.CG_OK) then
      print*,'*FAILED* cgp_coord_write_data_f (Coord_Z)'
      call cgp_error_exit_f()
    endif

  !---------- create empty "Elem" node in DB with nodes connectivity
  boundry_elems = 0 ! unsorted boundary elements
  first_cell_cgns = lbound(cell_connections,dim=2)
  last_cell_cgns  = ubound(cell_connections,dim=2)
  call cgp_section_write_f(file_idx,base_idx,zone_idx, &
    'Elem',HEXA_8,int(1,8),int(i,8),boundry_elems, &
    section_idx,ier) ! prepare to write connectivity data
    if (ier.ne.CG_OK) then
      print*,'*FAILED* cgp_section_write_f'
      call cgp_error_exit_f()
    endif

  !---------- fill "Elem" node in DB with nodes connectivity
  call cgp_elements_write_data_f(file_idx,base_idx,zone_idx, &
    section_idx,first_cell_cgns,last_cell_cgns,cell_connections,ier) ! write HEXA_8 element connectivity
      if (ier.ne.CG_OK) then
      print*,'*FAILED* cgp_elements_write_data_f (elements)'
      call cgp_error_exit_f()
    endif

  !---------- close DB
  call cgp_close_f(file_idx,ier)
  if (this < 2) write(6,'('' Successfully wrote unstructured grid to file'', &
    '' grid.cgns'')') ! close CGNS file

end subroutine Save_Mesh_Par_CGNS_Lib