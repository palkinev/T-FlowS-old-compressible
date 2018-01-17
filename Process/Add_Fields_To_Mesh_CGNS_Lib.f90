subroutine Add_Fields_To_Mesh_CGNS_Lib(sub)
!   Opens an existing CGNS file 'grid.cgns' that contains a simple 3-D         !
!   unstructured grid, and adds a flow solution (at cell centers) to it        !
!------------------------------------------------------------------------------!
!   The CGNS grid file 'grid.cgns' must already exist
!------------------------------------------------------------------------------!
!   All conventions for this lib are here:                                     !
!   https://cgns.github.io/CGNS_docs_current/sids/conv.html                    !
!------------------------------------------------------------------------------!
!   this function inherited syntax from                                        !
!   CGNS/src/Test_UserGuideCode/Fortran_code/write_flowvert_unst.F90           !
!   except that function writes to cell vertexes dimensional data              !
!------------------------------------------------------------------------------!
!   TO DO:                                                                     !
!   parallel writing                                                           !
!   tetra cells                                                                !
!   partial geometry                                                           !
!----------------------------------[Modules]-----------------------------------!
  use all_mod
  use par_mod
  use allp_mod
  use pro_mod
  use les_mod
  use rans_mod

  use cgns !on  Windows machines use "include 'cgnswin_f.h'"
!------------------------------------------------------------------------------!
  implicit none
  include 'mpif.h'
!---------------------------------[Parameters]----------------------------------!
  integer             :: sub
!-----------------------------------[Locals]-----------------------------------!
  ! CGNS variables
  integer             :: index_field,index_flow,index_zone,index_base,ier
  integer             :: index_file,iset
  character           :: solname*32

  ! Other variables
  real*8, allocatable :: flux_2_cgns(:)
  integer             :: c1, c2, s
!---------------------------------[Interface]----------------------------------!


  if (this < 2) write(6,'('' Program write_flowvert_unst'')')
  if (CG_BUILD_64BIT) then
    if (this < 2) write(6,'('' ...compiled in 64-bit mode, but not needed'')')
  end if

  !---------- open CGNS file for modify
  call cg_open_f('grid.cgns',CG_MODE_MODIFY,index_file,ier)
  if (ier .ne. CG_OK) call cg_error_exit_f

  if (this < 2) write(6,'('' created simple 3-D rho and p flow solution'')')

  !---------- write flow solution to existing cgns file

  index_base = 1 ! we know there is only one base (real working code would check!)

  index_zone = 1 ! we know there is only one zone (real working code would check!)

  solname = 'FlowSolution' ! define flow solution node name (user can give any name)
  call cg_sol_write_f(index_file,index_base,index_zone,solname, &
    CellCenter,index_flow,ier) ! create flow solution node with CellCentered data

  !---------- density
  call cg_field_write_f(index_file,index_base,index_zone,index_flow, &
    RealDouble,'Density',Rho % n(1:NC),index_field,ier)

  !---------- velocity
  call cg_field_write_f(index_file,index_base,index_zone,index_flow, &
    RealDouble,'Velocity_X',U % n(1:NC),index_field,ier) 
  call cg_field_write_f(index_file,index_base,index_zone,index_flow, &
    RealDouble,'Velocity_Y',V % n(1:NC),index_field,ier) 
  call cg_field_write_f(index_file,index_base,index_zone,index_flow, &
    RealDouble,'Velocity_W',W % n(1:NC),index_field,ier)

  !---------- presssure
  call cg_field_write_f(index_file,index_base,index_zone,index_flow, &
    RealDouble,'Pressure',P % n(1:NC),index_field,ier)

  !---------- flux
  allocate(flux_2_cgns(1:NC)); flux_2_cgns(:) = 0. !< [1:NC_or_NN]

  do s = 1, NS
      c1 = SideC(1,s)
      c2 = SideC(2,s)
      if(c2  > 0 .or. c2  < 0 .and. TypeBC(c2) == BUFFER) then
          flux_2_cgns(c1) = flux_2_cgns(c1) - Flux(s)
          if(c2  > 0) flux_2_cgns(c2) = flux_2_cgns(c2)+Flux(s)
      else
          flux_2_cgns(c1) = flux_2_cgns(c1)-Flux(s)
      end if
  end do

  call cg_field_write_f(index_file,index_base,index_zone,index_flow, &
    RealDouble,'Flux',flux_2_cgns(1:NC),index_field,ier)


  call cg_close_f(index_file,ier) ! close CGNS file

  if (this < 2) write(6,'('' Successfully added flow solution data to file'', &
   '' grid.cgns (unstructured)'')')
  if (this < 2) write(6,'(''   Note:  if the original CGNS file already had'', &
   '' a FlowSolution_t node,'')')
  if (this < 2) write(6,'(''          it has been overwritten'')')

  if ( allocated( flux_2_cgns ) ) deallocate( flux_2_cgns )
!==============================================================================!

end subroutine Add_Fields_To_Mesh_CGNS_Lib