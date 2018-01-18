subroutine Save_Mesh_CGNS_Lib(sub)                                             !
!   Writes 3-D unstructured grid and to CGNS grid file 'grid.cgns'             !
!------------------------------------------------------------------------------!
!   All conventions for this lib are here:                                     !
!   https://cgns.github.io/CGNS_docs_current/sids/conv.html                    !
!------------------------------------------------------------------------------!
!   Parallel functions are describes hereЖ                                     !
!   https://cgns.github.io/CGNS_docs_current/pcgns/library.html                !
!------------------------------------------------------------------------------!
!   this function inherited syntax from                                        !
!   CGNS/src/Test_UserGuideCode/Fortran_code/write_grid_unst.F90               !
!------------------------------------------------------------------------------!
!   TO DO:                                                                     !
!   parallel writing                                                           !
!   tetra cells                                                                !
!   partial geometry                                                           !
!----------------------------------[Modules]-----------------------------------!
  use all_mod
  use par_mod
  !use allp_mod
  !use pro_mod
  !use les_mod
  !use rans_mod

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
  integer            :: ielem_no ! 
  integer            :: ier, iset, iphysdim, icelldim

  integer(cgsize_t)  :: mesh_info(1, 3), cell_connections(8, NC)
  integer(cgsize_t)  :: first_cell, last_cell
  integer(cgsize_t)  :: firt_node, last_node
  character          :: basename*32, zonename*32

  ! Other variables
  integer             :: sub, sub_local

  character           :: dummy*100, nameTMP*80
  integer             :: Nnodes
  integer             :: c, n, i
  integer             :: nodes_in_current_cell ! 8 -> hex, 4-> tetra

  real, allocatable   :: x(:), y(:), z(:) ! nodes coordinates

  ! sub zone variables
  integer                :: NC_In_Sub_Zone
  REAL                   :: x1, x2, x3, x4
  REAL                   :: y1, y2, y3, y4
  REAL                   :: z1, z2
  logical                :: FileExist
  integer*1, allocatable :: SaveAreaIsEmpty(:)
  integer*1, allocatable :: SaveAreaIsEmptySum(:)
!---------------------------------[Interface]----------------------------------!
!==============================================================================!
  !-------------------------------------------------------------------!
  !   Parallel commands
  call cgp_mpi_comm_f(MPI_COMM_WORLD, ier) ! Set the MPI communicator for CGNS
  !call cgp_mpi_info_f(integer info, ier) ! Set the MPI info object for CGNS
  call cgp_mpi_info_f(0, ier) ! Set the MPI info object for CGNS
  call cgp_pio_mode_f(CGP_INDEPENDENT, ier) ! Set the parallel IO mode for CGNS

  !call cgp_error_exit_f()  ! Exit with error message 
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

  SaveAreaIsEmpty = 0    ! 0 -> not empty
  SaveAreaIsEmptySum = 0 ! 0 -> not empty
  NC_In_Sub_Zone = 0

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
    do c = 1, NC
      if (    xc(c) > x1 &
        .and. xc(c) < x2 &
        .and. xc(c) > x3 &
        .and. xc(c) < x4 &
        .and. yc(c) > y1 + (y4-y1)*(xc(c)-x1)/(x4-x1) &
        .and. yc(c) < y2 + (y3-y2)*(xc(c)-x2)/(x3-x2) &
        .and. zc(c) > z1 &
        .and. zc(c) < z2 ) then

        NC_In_Sub_Zone = NC_In_Sub_Zone + 1

      end if
    end do

    if (NC_In_Sub_Zone == 0) then
      SaveAreaIsEmpty(sub) = 1
    end if

  end if ! else draw full area

  if (FileExist) then ! if there are instructions to cut the area
    call MPI_ALLREDUCE(      &
        SaveAreaIsEmpty,     & ! send buffer
        SaveAreaIsEmptySum,  & ! recv buffer
        Npro,                & ! length
        MPI_INTEGER,         & ! datatype
        MPI_SUM,             & ! operation
        MPI_COMM_WORLD,      &
        ier )
    end if

  ! in case of single processor
  if (sub_local == 0)  sub = 0

  !-------------------------------------!
  ! Read .gmv block for additional data !
  !-------------------------------------!

  call NamFil(sub, nameTMP, '.gmv', len_trim('.gmv'))
  open(9, FILE = nameTMP)
  if (this <2 ) write(*,*) 'Now reading the file: ', nameTMP

  call ReadC(9,inp,tn,ts,te)  ! read line "gmvinput ascii"
  call ReadC(9,inp,tn,ts,te)  ! read "nodes  $number_of_nodes"
  read(inp(ts(1):te(1)),*) dummy ! "nodes" word
  read(inp(ts(2):te(2)),*) Nnodes ! nodes count

  !---------- read x, y, z of cell nodes consequently
  allocate(x(Nnodes)); x=0.
  allocate(y(Nnodes)); y=0.
  allocate(z(Nnodes)); z=0.
  do n = 1, Nnodes
    call ReadC(9,inp,tn,ts,te)
    read(inp(ts(1):te(1)),*) x(n)
  end do
  do n = 1, Nnodes
    call ReadC(9,inp,tn,ts,te)
    read(inp(ts(1):te(1)),*) y(n)
  end do
  do n = 1, Nnodes
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

  !allocate(cell_connections(8, NC)); cell_connections = 0
  nodes_in_current_cell = 0
  ielem_no = 0 ! ??? index of first element in current cell
  
  do c = 1, NC
    call ReadC(9,inp,tn,ts,te)
    read(inp(ts(1):te(1)),*) dummy ! "hex"/"tetra"/"pyramid"/ word
    read(inp(ts(2):te(2)),*) nodes_in_current_cell ! nodes count to make a cell

    call ReadC(9,inp,tn,ts,te)
    do n = 1, nodes_in_current_cell
      read(inp(ts(n):te(n)),*) cell_connections(n, c) ! nodes index is the same as in gmv
    end do
    ielem_no = ielem_no + 1
  end do

  close(9) ! close .gmv file

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
  basename = "Base 1"
  icelldim = 3
  iphysdim = 3
  call cg_base_write_f(file_idx,basename,icelldim,iphysdim,base_idx,ier)
    if (ier.ne.CG_OK) then
       print*,'*FAILED* cg_base_write_f'
       call cgp_error_exit_f()
    endif
  zonename = "Zone  1" !define zone name (user can give any name)
  mesh_info(1,1) = Nnodes ! vertex count
  mesh_info(1,2) = NC_In_Sub_Zone ! cell count
  mesh_info(1,3) = 0 ! boundary vertex size (zero if elements not sorted)
  call cg_zone_write_f(file_idx,base_idx,zonename,mesh_info, &
    Unstructured,zone_idx,ier) ! create zone
    if (ier.ne.CG_OK) then
     print*,'*FAILED* cg_zone_write_f'
     call cgp_error_exit_f()
  endif

!--------- parallel

  firt_node = 1
  last_node = Nnodes
  first_cell = 1 ! ??? start element
  last_cell = NC ! ??? start element

  !---------- create empty "CoordinateX" node in DB
  call cgp_coord_write_f(file_idx,base_idx,zone_idx,RealDouble, &
    'CoordinateX',coord_idx,ier)
    if (ier.ne.CG_OK) then
      print*,'*FAILED* cgp_coord_write_f (Coord_X)'
      call cgp_error_exit_f()
    endif
  !---------- fill "CoordinateX" node in DB
  call cgp_coord_write_data_f(file_idx,base_idx,zone_idx,coord_idx, &
    firt_node,last_node,x,ier)
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
      firt_node,last_node,y,ier)
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
    firt_node,last_node,z,ier) ! write z grid coordinates
    if (ier.ne.CG_OK) then
      print*,'*FAILED* cgp_coord_write_data_f (Coord_Z)'
      call cgp_error_exit_f()
    endif
!----------------- sequential

! call cg_coord_write_f(file_idx,base_idx,zone_idx,RealDouble, &
!     'CoordinateX',x,coord_idx,ier) ! write x grid coordinates 
! call cg_coord_write_f(file_idx,base_idx,zone_idx,RealDouble, &
!     'CoordinateY',y,coord_idx,ier) ! write y grid coordinates 
! call cg_coord_write_f(file_idx,base_idx,zone_idx,RealDouble, &
!     'CoordinateZ',z,coord_idx,ier) ! write z grid coordinates 

!----------------- sequential
!  boundry_elems = 0 ! unsorted boundary elements
!  call cg_section_write_f(file_idx,base_idx,zone_idx, &
!       'Elem',HEXA_8,first_cell,last_cell,boundry_elems,cell_connections, &
!       section_idx,ier) ! write HEXA_8 element connectivity


 ! CALL cgp_section_write_f(F, B, Z, 'Hex', HEXA_8, start_1, totelems, 0, E, ierr)
 ! CALL cgp_elements_write_data_f(F, B, Z, E, start, END, ie, ierr)

!--------- parallel
  boundry_elems = 0 ! unsorted boundary elements

  !---------- create empty "Elem" node in DB with nodes connectivity
  call cgp_section_write_f(file_idx,base_idx,zone_idx, &
    'Elem',HEXA_8,first_cell,last_cell,boundry_elems, &
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