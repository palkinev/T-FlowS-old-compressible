subroutine Add_BC_To_Mesh_CGNS_Lib(sub, NCsub, nameIN)
!   Writes 3-D unstructured grid and to CGNS grid file 'grid.cgns'
!------------------------------------------------------------------------------!
!   All conventions for this lib are here:                                     !
!   https://cgns.github.io/CGNS_docs_current/sids/conv.html                    !
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
  use allp_mod
  !use pro_mod
  use les_mod
  use rans_mod

  use cgns !on  Windows machines use "include 'cgnswin_f.h'"
!------------------------------------------------------------------------------!
  implicit none
  include 'mpif.h'
!-----------------------------------[Locals]-----------------------------------!
  integer            :: commsize, commrank, ierr
  integer            :: maxelemi
  integer            :: maxelemj

  ! CGNS variables
  integer ifirstnode, nbdyelem
  integer index_file, index_section, ielem_no
  integer index_base, index_zone, index_coord
  integer ier, iset, iphysdim, icelldim

  integer(cgsize_t) isize(1, 3), ielem(8, NC)
  integer(cgsize_t) nelem_start, nelem_end
  character basename*32, zonename*32


  integer             :: error_MPI

  character           :: nameIN*(*)
  integer             :: NCsub, sub, sub_local

  character           :: dummy*100, nameTMP*80, nameTMP2*80
  integer             :: NNsub, NCsub_new
  integer             :: c, n, s, i, c1, c2
  integer             :: nodes_in_current_cell ! 8 -> hex, 4-> tetra

  real, allocatable    :: x(:), y(:), z(:) ! nodes coordinates

  ! mesh 2 CGNS_IO vars
  integer                       :: NC_new
  REAL                          :: x1, x2, x3, x4
  REAL                          :: y1, y2, y3, y4
  REAL                          :: z1, z2
  logical                       :: FileExist
  integer*1, allocatable        :: SaveAreaIsEmpty(:)
  integer*1, allocatable        :: SaveAreaIsEmptySum(:)
!---------------------------------[Interface]----------------------------------!
!==============================================================================!

end subroutine Add_BC_To_Mesh_CGNS_Lib