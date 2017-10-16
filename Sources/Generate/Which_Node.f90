!==============================================================================!
  integer function Which_Node(c, n) 
!------------------------------------------------------------------------------!
!   Returns the local number (1-8) of node n in cell c.                        !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use all_mod
  use gen_mod
!------------------------------------------------------------------------------! 
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer :: c, n
!-----------------------------------[Locals]-----------------------------------!
  integer :: i, j
!==============================================================================!

  Which_Node=0
  if (c  < 0) then 
    write(*,*) 'Which node: Cell non existent !'
    return
  endif

  ! Try the node only 
  do i=1,8
    if( CellN(c,i)  ==  n) then
      goto 1
    end if 
  end do     

  ! If it failed try his twins also
  do j=1,TwinN(n,0)
    do i=1,8
      if( CellN(c,i)  ==  TwinN(n,j)) then
        goto 1
      end if 
    end do     
  end do

  Which_Node = 0
  write(*,*) 'Which node: Trouble, node not found !'
  write(*,*) 'x, y, z = ', x_node(n), y_node(n), z_node(n)
  write(*,*) 'cell    = ', c, level(c)
  return

1 Which_Node = i
  return

  end function Which_Node
