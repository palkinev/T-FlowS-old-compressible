subroutine Save_Mesh_Seq_CGNS_Lib(sub)                                         !
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
!   TO DO:                                                                     !
!   tetra cells                                                                !
!   partial geometry                                                           !
!----------------------------------[Modules]-----------------------------------!
  use all_mod
  use cgns !on  Windows machines use "include 'cgnswin_f.h'"
!------------------------------------------------------------------------------!
  implicit none
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

  ! Other variables
  integer             :: sub

  character           :: dummy*100, nameTMP*80
  integer*8           :: n_nodes
  integer             :: c, n, i
  integer             :: nodes_in_current_cell ! 8 -> hex, 4-> tetra

  real, allocatable   :: x(:), y(:), z(:) ! cell node coordinates

  ! sub zone variables
  integer*8              :: NC_In_Sub_Zone
!---------------------------------[Interface]----------------------------------!
!==============================================================================!

  write(6,'('' Program write_grid_unst'')')
  if (CG_BUILD_64BIT) then
    write(6,'('' ...using 64-bit mode for particular integers'')')
  end if

  !-------------------------------------!
  ! Read .gmv block for additional data !
  !-------------------------------------!

  call NamFil(sub, nameTMP, '.gmv', len_trim('.gmv'))
  open(9, FILE = nameTMP)
  write(*,*) 'Now reading the file: ', nameTMP

  call ReadC(9,inp,tn,ts,te)  ! read line "gmvinput ascii"
  call ReadC(9,inp,tn,ts,te)  ! read "nodes  $number_of_nodes"
  read(inp(ts(1):te(1)),*) dummy ! "nodes" word
  read(inp(ts(2):te(2)),*) n_nodes ! nodes count

  !---------- allocate x, y, z
  allocate(x(1: n_nodes), stat = ier); x = 0.
    if (ier.ne.0)then
       print*, '*FAILED* allocation of x'
       call cg_error_exit_f()
    endif
  allocate(y(1: n_nodes), stat = ier); y = 0.
    if (ier.ne.0)then
       print*, '*FAILED* allocation of y'
       call cg_error_exit_f()
    endif
  allocate(z(1: n_nodes), stat = ier); z = 0.
    if (ier.ne.0)then
       print*, '*FAILED* allocation of z'
       call cg_error_exit_f()
    endif

  !---------- read x, y, z of cell nodes consequently
  do n = 1, n_nodes
    call ReadC(9,inp,tn,ts,te)
    read(inp(ts(1):te(1)),*) x(n)
  end do
  do n = 1, n_nodes
    call ReadC(9,inp,tn,ts,te)
    read(inp(ts(1):te(1)),*) y(n)
  end do
  do n = 1, n_nodes
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

  !---------- read nodes connection for cells
  ! gmv syntax ex:
  ! hex 8
  ! 1        2       67       66      261      262      327      326
  !

  !---------- allocate cell_connections
  allocate(cell_connections(8, 1: NC_In_Sub_Zone), stat = ier)
  cell_connections = 0.
    if (ier.ne.0)then
       print*, '*FAILED* allocation of cell_connections'
       call cg_error_exit_f()
    endif

  !---------- read cell_connections
  do c = 1, NC_In_Sub_Zone
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

  !-------------------------!
  ! CGNS mesh writing block !
  !-------------------------!
  write(6,'('' created simple 3-D grid points'')')

  !----------   open CGNS file for write
  
  call cg_open_f('grid.cgns',CG_MODE_WRITE,file_idx,ier)
    if (ier.ne.CG_OK) then
       print*,'*FAILED* cg_open_f'
       call cg_error_exit_f()
    endif

  !----------   create a base node
  base_name = "Base 1"
  call cg_base_write_f(file_idx,base_name,3,3,base_idx,ier)
    if (ier.ne.CG_OK) then
       print*,'*FAILED* cg_base_write_f'
       call cg_error_exit_f()
    endif

  zone_name = "Zone 1" !define zone name (user can give any name)

  !---------- sum of nodes
  mesh_info(1) = n_nodes
  !---------- total cell count
  mesh_info(2) = NC_In_Sub_Zone
  mesh_info(3) = 0 ! boundary vertex size (zero if elements not sorted)

  call cg_zone_write_f(file_idx,base_idx,zone_name,mesh_info, &
    Unstructured,zone_idx,ier) ! create zone
    if (ier.ne.CG_OK) then
     print*,'*FAILED* cg_zone_write_f'
     call cg_error_exit_f()
  endif

!---------- create "CoordinateX" node in DB and fill it
  call cg_coord_write_f(file_idx,base_idx,zone_idx,RealDouble, &
    'CoordinateX',x,coord_idx,ier)
!---------- create "CoordinateY" node in DB and fill it
  call cg_coord_write_f(file_idx,base_idx,zone_idx,RealDouble, &
    'CoordinateY',y,coord_idx,ier)
!---------- create "CoordinateZ" node in DB and fill it
  call cg_coord_write_f(file_idx,base_idx,zone_idx,RealDouble, &
    'CoordinateZ',z,coord_idx,ier)

  !---------- create "Elem" node in DB with nodes connectivity and fill it
  boundry_elems = 0 ! unsorted boundary elements
  first_cell_cgns = 1
  last_cell_cgns = NC_In_Sub_Zone
  call cg_section_write_f(file_idx,base_idx,zone_idx, &
    'Elem',HEXA_8,first_cell_cgns,last_cell_cgns,boundry_elems,cell_connections, &
    section_idx,ier) ! write HEXA_8 element connectivity

  !---------- close DB
  call cg_close_f(file_idx,ier)
  write(6,'('' Successfully wrote unstructured grid to file'', &
    '' grid.cgns'')') ! close CGNS file

end subroutine Save_Mesh_Seq_CGNS_Lib