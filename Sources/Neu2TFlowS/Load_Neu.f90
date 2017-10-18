!==============================================================================!
  subroutine Load_Neu() 
!------------------------------------------------------------------------------!
!   Reads the Fluents (Gambits) neutral file format.                           !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use all_mod 
  use gen_mod 
  use neu_mod 
  use par_mod 
!------------------------------------------------------------------------------!
  implicit none
!-----------------------------------[Locals]-----------------------------------!
  character(len=300)  :: Line
  character(len=130)  :: name_in
  integer             :: i, j, dum1, dum2, dum3
  integer,allocatable :: temp(:)
  integer             :: c, dir, type
!==============================================================================!

  name_in = name
  name_in(len_trim(name)+1:len_trim(name)+4) = '.neu'

  open(9,file=name_in)
  write(*,*) '# Now reading the file: ', name_in

  ! Skip first 6 lines
  do i=1,6
    call Read_Line(9,Line,tn,ts,te)
  end do 

  ! Read the line which contains usefull information  
  call Read_Line(9,Line,tn,ts,te)
  read(Line(ts(1):te(1)),*) NN  
  read(Line(ts(2):te(2)),*) NC
  read(Line(ts(3):te(3)),*) NBloc
  read(Line(ts(4):te(4)),*) NBS

  write(*,*) '# Total number of nodes:  ', NN
  write(*,*) '# Total number of cells:  ', NC
  write(*,*) '# Total number of blocks: ', NBloc
  write(*,*) '# Total number of boundary sections: ', NBS

  ! Count the boundary cells
  NbC = 0
  do 
    call Read_Line(9,Line,tn,ts,te)
    if( Line(ts(1):te(1)) == 'BOUNDARY' ) then
      do j=1,NBS
        if(j>1) call Read_Line(9,Line,tn,ts,te) ! BOUNDARY CONDITIONS
        call Read_Line(9,Line,tn,ts,te)
        read(Line(ts(3):te(3)),*) dum1  
        NBc = NBc + dum1 
        do i=1,dum1
          read(9,*) c, dum2, dir
        end do
        call Read_Line(9,Line,tn,ts,te)         ! ENDOFSECTION
      end do
      write(*,*) '# Total number of boundary cells: ', NbC
      go to 1
    end if
  end do 

1 rewind(9)

  ! Skip first 7 lines
  do i=1,7
    call Read_Line(9,Line,tn,ts,te)
  end do 

  ! Allocate memory =--> carefull, there is no checking!
  allocate(x_node(NN)); x_node=0
  allocate(y_node(NN)); y_node=0
  allocate(z_node(NN)); z_node=0

  allocate(material(-NbC:NC));     material=0 
  allocate(BCtype(NC,6));          BCtype=0
  allocate(BCmark(-NbC-1:-1));     BCmark=0
  allocate(CellN(-NbC-1:NC,-1:8)); CellN=0
  allocate(SideC(0:2,NC*5));       SideC=0    
  allocate(SideN(NC*5,0:4));       SideN=0
  n_copy = 1000000 
  allocate(CopyC(-NbC:-1));  CopyC=0
  allocate(CopyS(2,n_copy));  CopyS=0

  allocate(NewN(NN));        NewN=0  
  allocate(NewC(-NbC-1:NC)); NewC=0  
  allocate(NewS(NC*5));      NewS=0  

  allocate(proces(0:NC)); proces=0

  allocate(temp(NC)); temp=0

  ! Skip one line 
  call Read_Line(9,Line,tn,ts,te)

  !--------------------------------!
  !   Read the nodal coordinates   !
  !--------------------------------!
  call Read_Line(9,Line,tn,ts,te)          ! NODAL COORDINATES
  do i=1,NN
    call Read_Line(9,Line,tn,ts,te)
    read(Line(ts(2):te(2)),*) x_node(i)  
    read(Line(ts(3):te(3)),*) y_node(i)
    read(Line(ts(4):te(4)),*) z_node(i) 
  end do
  call Read_Line(9,Line,tn,ts,te)          ! ENDOFSECTION

  !-----------------------------!
  !   Read nodes of each cell   !
  !-----------------------------!
  call Read_Line(9,Line,tn,ts,te)          ! ELEMENTS/CELLS
  do i=1,NC
    read(9,'(I8,1X,I2,1X,I2,1X,7I8:/(15X,7I8:))') dum1, dum2, &
           CellN(i,0), (CellN(i,j), j=1,CellN(i,0))
  end do
  call Read_Line(9,Line,tn,ts,te)          ! ENDOFSECTION

  ! Read block data
  do j=1,NBloc
    call Read_Line(9,Line,tn,ts,te)        ! ELEMENT GROUP
    call Read_Line(9,Line,tn,ts,te)
    read(Line(ts(4):te(4)),'(I10)') dum1  
    call Read_Line(9,Line,tn,ts,te)        ! block*
    call Read_Line(9,Line,tn,ts,te)        ! 0
    read(9,'(10I8)') (temp(i), i=1,dum1)
    do i=1,dum1
      material(temp(i)) = j
    end do
    call Read_Line(9,Line,tn,ts,te)        ! ENDOFSECTION
  end do

  !-------------------------!
  !   Boundary conditions   !
  !-------------------------!
  write(*,*) '# NBS=', NBS 
  do j=1,NBS
    call Read_Line(9,Line,tn,ts,te)        ! BOUNDARY CONDITIONS
    call Read_Line(9,Line,tn,ts,te)
    if( Line(ts(1):te(1)) == 'one')      type =  1
    if( Line(ts(1):te(1)) == 'two')      type =  2
    if( Line(ts(1):te(1)) == 'three')    type =  3
    if( Line(ts(1):te(1)) == 'four')     type =  4
    if( Line(ts(1):te(1)) == 'five')     type =  5
    if( Line(ts(1):te(1)) == 'six')      type =  6
    if( Line(ts(1):te(1)) == 'seven')    type =  7
    if( Line(ts(1):te(1)) == 'eight')    type =  8
    if( Line(ts(1):te(1)) == 'nine')     type =  9
    if( Line(ts(1):te(1)) == 'ten')      type = 10
    if( Line(ts(1):te(1)) == 'eleven')   type = 11
    if( Line(ts(1):te(1)) == 'onetwo')   type = 12
    if( Line(ts(1):te(1)) == 'onethree') type = 13
    if( Line(ts(1):te(1)) == 'onefour')  type = 14
    if( Line(ts(1):te(1)) == 'onefive')  type = 15
    if( Line(ts(1):te(1)) == 'onesix')   type = 16
    if( Line(ts(1):te(1)) == 'oneseven') type = 17
    if( Line(ts(1):te(1)) == 'oneeight') type = 18
    if( Line(ts(1):te(1)) == 'onenine')  type = 19
    write(*,*)  Line(ts(1):te(1)), type 
    read(Line(ts(3):te(3)),'(I8)') dum1  
    do i=1,dum1
      read(9,*) c, dum2, dir
      BCtype(c,dir) = type 
    end do
    call Read_Line(9,Line,tn,ts,te)        ! ENDOFSECTION
  end do

  close(9)

  end subroutine Load_Neu