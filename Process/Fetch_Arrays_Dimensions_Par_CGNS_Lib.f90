subroutine Fetch_Arrays_Dimensions_Par_CGNS_Lib(idx, NN_or_NC)
!   Fetches correct dimensions for arrays in CGNS lib dependent functions      !
!------------------------------------------------------------------------------!
!   Arrays structure in CGNS parallel functions are strictly followings:       !
!   Processor:    |        P_1        |               P_2               | ...  !
!   x,y,z:        |      (1 : NN_1)   |       NN_1 + 1 : NN_1 + NN_2    | ...  ! 
!   Connections:  |   (8, 1 : NC_1)   |   (8, NC_1 + 1 : NC_1 + NC_2)   | ...  ! 
!----------------------------------[Modules]-----------------------------------!
  use par_mod
!------------------------------------------------------------------------------!
  implicit none
  include 'mpif.h'
!-----------------------------------[Locals]-----------------------------------!
  integer*8   :: idx, NN_or_NC
  integer*8 :: Array_at_root(1:Npro)
  integer   :: i, ier
!------------------------------------------------------------------------------!

  Array_at_root = 0

  call wait
  call mpi_gather(  &
    NN_or_NC,       & ! send number of nodes/cells, which is
    1,              & ! 1 element
    MPI_INTEGER8,   & ! of 32-bit integer type
    Array_at_root,  & ! to the root array
    1,              & ! which is 1 per processor long message.
    MPI_INTEGER8,   & ! Received data type is 32-bit integer
    0,              & ! at root processor
    mpi_comm_world, & ! communicator
    ier)              ! mpi_error

  if (this < 2) then

    do i = 2, Npro
      Array_at_root(i) = 1 + Array_at_root(i-1)
    end do
    Array_at_root(1) = 1

  end if

  call wait
  call mpi_scatter( &
    Array_at_root,  & ! send number of nodes/cells, which is
    1,              & ! 1 elements per processor
    MPI_INTEGER8,   & ! of 32-bit integer type
    idx,            & ! to integer "idx"
    Npro,           & ! which is Npro elements long message.
    MPI_INTEGER8,   & ! Received data type is 32-bit integer
    0,              & ! from the proc 1
    mpi_comm_world, & ! communicator
    ier)              ! mpi_error

end subroutine Fetch_Arrays_Dimensions_Par_CGNS_Lib