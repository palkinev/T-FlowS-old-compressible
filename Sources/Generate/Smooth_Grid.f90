!==============================================================================!
  subroutine Smooth_Grid(grid)
!------------------------------------------------------------------------------!
!   Smooths the grid lines by a Laplacian-like algorythm.                      !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use all_mod
  use gen_mod
  use Grid_Mod
!------------------------------------------------------------------------------! 
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!-----------------------------------[Locals]-----------------------------------!
  integer              :: c, n, i, j, k, m
  real                 :: x_new_tmp, y_new_tmp, z_new_tmp    
  real                 :: x_max, y_max, z_max, x_min, y_min, z_min 
  real                 :: x1, y1, z1, x8, y8, z8 
  integer              :: reg
  real, allocatable    :: x_node_new(:), y_node_new(:), z_node_new(:) 
  integer, allocatable :: node_to_nodes(:,:)
!==============================================================================!

  ! Allocate memory for additional arrays
  allocate(x_node_new(grid % max_n_nodes)); x_node_new=0
  allocate(y_node_new(grid % max_n_nodes)); y_node_new=0
  allocate(z_node_new(grid % max_n_nodes)); z_node_new=0
  allocate(node_to_nodes(grid % max_n_nodes,0:40)); node_to_nodes=0

  write(*,*) '# Now smoothing the cells. This may take a while !' 

  do n=1,Nn
    node_to_nodes(n,0) = 0
  end do 

  x_max=-HUGE
  y_max=-HUGE
  z_max=-HUGE
  x_min=+HUGE
  y_min=+HUGE
  z_min=+HUGE
  do n=1,Nn
    x_max=max(grid % xn(n), x_max) 
    y_max=max(grid % yn(n), y_max) 
    z_max=max(grid % zn(n), z_max) 
    x_min=min(grid % xn(n), x_min) 
    y_min=min(grid % yn(n), y_min) 
    z_min=min(grid % zn(n), z_min) 
  end do 

  !-----------------------!
  !   Connect the nodes   !
  !-----------------------!
  do c=1,Nc                        ! through cells
    do i=1,8                       ! through nodes of a cell 
      n = grid % cells_n(i,c)      ! first cell
      do j=1,8                     ! through nodes of a cell 
        m = grid % cells_n(j,c)    ! second cell 
        if(n /=  m) then 
          do k=1,node_to_nodes(n,0)
            if(node_to_nodes(n,k) == m) goto 10            
          end do
          node_to_nodes(n,0)=node_to_nodes(n,0)+1
          node_to_nodes(n,node_to_nodes(n,0)) = m
        end if
10      end do
    end do
  end do 

  !----------------------------!
  !   Browse through regions   !
  !----------------------------!
  do reg=1,n_smoothing_regions
    if( ( .not. smooth_in_x(reg) ) .and.  &
        ( .not. smooth_in_y(reg) ) .and.  &
        ( .not. smooth_in_z(reg) ) ) then
      do n=1,Nn
        x1=smooth_regions(reg,1)
        y1=smooth_regions(reg,2)
        z1=smooth_regions(reg,3)
        x8=smooth_regions(reg,4)
        y8=smooth_regions(reg,5)
        z8=smooth_regions(reg,6)
        if( (x1<=grid % xn(n)) .and. (grid % xn(n)<=x8) .and. &
            (y1<=grid % yn(n)) .and. (grid % yn(n)<=y8) .and. &
            (z1<=grid % zn(n)) .and. (grid % zn(n)<=z8) ) then
          node_to_nodes(n,0) = 0
        endif
      end do
    end if
  end do

  !---------------------!
  !   Smooth the grid   !
  !---------------------!
  do reg=1,n_smoothing_regions
    write(*,*) '# Now smoothing region ',reg,' with:',              &
                smooth_iters(reg), ' iterations.'

    do j=1,smooth_iters(reg)         

      ! Calculate new coordinates using the old values 
      do n=1,Nn
        if(node_to_nodes(n,0)   >  0) then
          x_new_tmp=0.0
          y_new_tmp=0.0
          z_new_tmp=0.0
          do i=1,node_to_nodes(n,0)
            x_new_tmp = x_new_tmp + grid % xn(node_to_nodes(n,i))
            y_new_tmp = y_new_tmp + grid % yn(node_to_nodes(n,i))
            z_new_tmp = z_new_tmp + grid % zn(node_to_nodes(n,i))
          end do
          x_new_tmp = x_new_tmp / (1.0*node_to_nodes(n,0)) 
          y_new_tmp = y_new_tmp / (1.0*node_to_nodes(n,0)) 
          z_new_tmp = z_new_tmp / (1.0*node_to_nodes(n,0)) 
          if(grid % xn(n) > 0.001*x_min .and.  &
             grid % xn(n) < 0.999*x_max)       &
          x_node_new(n) = (1.0-smooth_relax(reg)*walln(n))*grid % xn(n) &
                        +      smooth_relax(reg)*walln(n) *x_new_tmp
          if(grid % yn(n) > 0.001*y_min .and.  &
             grid % yn(n) < 0.999*y_max)       &
          y_node_new(n) = (1.0-smooth_relax(reg)*walln(n))*grid % yn(n) &
                        +      smooth_relax(reg)*walln(n) *y_new_tmp
          if(grid % zn(n) > 0.001*z_min .and.  &
             grid % zn(n) < 0.999*z_max)       &
          z_node_new(n) = (1.0-smooth_relax(reg)*walln(n))*grid % zn(n) &
                        +      smooth_relax(reg)*walln(n)* z_new_tmp
        end if
      end do

      ! Update coordinates
      do n=1,Nn
        if(node_to_nodes(n,0)   >  0) then

          x1=smooth_regions(reg,1)
          y1=smooth_regions(reg,2)
          z1=smooth_regions(reg,3)
          x8=smooth_regions(reg,4)
          y8=smooth_regions(reg,5)
          z8=smooth_regions(reg,6)

          if( (x1 <= grid % xn(n)) .and.  &
              (grid % xn(n) <= x8) .and.  &
              (y1 <= grid % yn(n)) .and.  &
              (grid % yn(n) <= y8) .and.  &
              (z1 <= grid % zn(n)) .and.  &
              (grid % zn(n) <= z8) ) then

            if(smooth_in_x(reg)) then
              if(grid % xn(n) > 0.001*x_min .and. &
                 grid % xn(n) < 0.999*x_max)  &
              grid % xn(n)=x_node_new(n)
            end if

            if(smooth_in_y(reg)) then
              if(grid % yn(n) > 0.001*y_min .and.  &
                 grid % yn(n) < 0.999*y_max)  &
              grid % yn(n)=y_node_new(n)
            end if 

            if(smooth_in_z(reg)) then
              if(grid % zn(n) > 0.001*z_min .and.  &
                 grid % zn(n) < 0.999*z_max)  &
              grid % zn(n)=z_node_new(n)
            end if 

          end if  ! if the point belongs to region
        end if    ! if the point is not excluded from smoothing
      end do      ! through nodes
    end do
  end do

  deallocate(x_node_new)
  deallocate(y_node_new)
  deallocate(z_node_new)
  deallocate(node_to_nodes)

  end subroutine
