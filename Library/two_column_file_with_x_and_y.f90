!————————————————————————————————————————————————————————————————————————————————————————
!---- SUB MODULE TO HANDLE FILES
!————————————————————————————————————————————————————————————————————————————————————————
module two_column_file_with_x_and_y
! known bugs: not working if total point count < 3
    implicit none
    !private

    !————————————————————————————————————————————————————————————————————————————————————————
    type, public :: file     ! type for two column file
        character(len=100) :: filename ! it has filename and must be given from outside
        logical  :: isset   ! = true if file is not "", else false
        integer  :: file_id ! it is assigned file ID - it is free
        integer  :: lines   ! it is file line count
        real, allocatable  :: x(:), y(:)   ! first and second column [has to be sorted on x]
        logical            :: homogeneous_mesh ! if mesh is homogeneous
        real               :: dx ! dx, if mesh is homogeneous
    end type file

    ! functions
    public :: set_filename
    public :: file_is_set
    public :: assign_file_id
    public :: get_lines_count
    public :: assign_x_and_y
    public :: assign_dx
    public :: print_info
    public :: linear_interp
        
!————————————————————————————————————————————————————————————————————————————————————————
contains
    !————————————————————————————————————————————————————————————————————————————————————————
    subroutine set_filename(this,input)
        implicit none
        type(file)               :: this
        character(*), intent(in) :: input
        !integer,parameter           :: seed = 61
900     format (" # WARNING: ZERO VALID LINES FOUND. ZERO VALUES WILL BE USED FOR: ", A)

        !call srand(seed)

        this % homogeneous_mesh = .false.

        this % filename = input
        this % isset = .true.

        if (input == "") then
            this % isset = .false.
            return
        end if

        call assign_file_id (this)
        call get_lines_count(this)
        if (this % lines == 0) this % isset = .false.
        if (this % isset) call assign_x_and_y(this)
        if (this % isset) call assign_dx(this)

        if (.not.this % isset) write(*,900) trim(this % filename)

    end subroutine set_filename
    !————————————————————————————————————————————————————————————————————————————————————————
    function file_is_set(this)
        implicit none
        type(file)  :: this
        logical     :: file_is_set

        file_is_set = this % isset

    end function file_is_set
    !————————————————————————————————————————————————————————————————————————————————————————
    subroutine assign_file_id (this)
        implicit none
        !type(file), intent(in) :: this
        type(file) :: this
        logical itsopen
        integer :: FID

        itsopen = .true. !initial value to trigger check below

        do while ( itsopen ) ! it will stop when FID is unique
            FID = 115 + mod(irand(0),1500)
            inquire(unit=FID, opened=itsopen)
        end do
800     format (" # TWO COLUMN MODULE: FILE OPENED: ", A)
        write(*,800) trim(this % filename)
        this % file_id = FID
    end subroutine assign_file_id
    !————————————————————————————————————————————————————————————————————————————————————————
    subroutine get_lines_count(this)
        implicit none
        type(file)    :: this
        logical       :: exist
        integer       :: file_lines_count, EOF
        character(10) :: dummy

        !            800 format (" # FILE OPENED: ",A)
        !            810 format (" # FILE CLOSED: ",A)
        !           900 format (" # FILE READING ERROR: ",A)
910     format (" # FILE DOES NOT EXIST: ",A)

        inquire(file = trim(this % filename), exist=exist)
        if (exist) then
            open(this % file_id, file=trim(this % filename), status="old", action="read")
            !write(*,800) trim(this % filename)
        else
            write(*,910) trim(this % filename)
        end if

        file_lines_count = -1
        EOF = 0

        do while ( EOF == 0 ) !until end of file
            read(this % file_id,*,iostat=EOF) dummy
            !write(*,*) "string= ", dummy, "EOF= ", EOF
            if (EOF > 0)  then ! some reading error
                !write(*,900) this % file_id
            elseif (EOF < 0) then ! end of file or file is closed
                !rewind(this % file_ID)
                close(this % file_id)
                !write(*,810) trim(this % filename)
            end if

            if ( dummy == '') then
                ! empty string
            elseif(                    &
                dummy(1:1) == '!' .or. &
                dummy(1:1) == '#' .or. &
                dummy(1:1) == '$' .or. &
                dummy(1:1) == '%' .or. &
                dummy(1:1) == '/' .or. &
                dummy(1:1) == '&' .or. &
                dummy(1:1) == '[' .or. &
                dummy(1:1) == ']' .or. &
                dummy(1:1) == '{' .or. &
                dummy(1:1) == '}' .or. &
                dummy(1:1) == ':' .or. &
                dummy(1:1) == ';' .or. &
                dummy(1:1) == '@' .or. &
                dummy(1:1) == '<' .or. &
                dummy(1:1) == '>' .or. &
                dummy(1:1) == ' ' .or. &
                dummy(1:1) == '*' ) then
            else
                file_lines_count = file_lines_count + 1
                !write(*,*) "line ccount = ", file_lines_count
            endif
        end do

        this % lines = file_lines_count

    end subroutine get_lines_count
    !————————————————————————————————————————————————————————————————————————————————————————
    subroutine assign_x_and_y (this)
        !use all_mod
        implicit none
        type(file)  :: this
        logical     :: exist
        integer     :: i
        integer     :: tn, ts(300), te(300)
        character   :: inp*300

        allocate (this % x(1:this % lines)); this % x = 0.0
        allocate (this % y(1:this % lines)); this % y = 0.0


        inquire(file=trim(this % filename), exist=exist)
        if (exist) then
            open(this % file_id, file=trim(this % filename), status="old", action="read")
        end if

        do i=1, this % lines
            call ReadC(this % file_id,inp,tn,ts,te)
            read(inp(ts(1):te(1)),*) this % x(i)
            read(inp(ts(2):te(2)),*) this % y(i)
            !write(*,*) "x=", this % x(i), "f= ",this % y(i)
        end do

        close(this % file_id)
        !write(*,810) this % filename

    end subroutine assign_x_and_y
    !————————————————————————————————————————————————————————————————————————————————————————
    subroutine assign_dx (this)
        implicit none
        type(file) :: this
        real :: dx_0
        integer :: i

275     format (" # File", A, " has not homogeneous mesh ")          ! output style
285     format (" # File", A, " has homogeneous mesh ")              ! output style
        this % dx = 0.
        dx_0 = abs( this % x(2) - this % x(1) )


        do i = 2, this % lines - 1
            this % dx = this % x(i+1) - this % x(i)

            if ( abs( this % dx - dx_0 ) > 1.0D-16  ) then
                write(*,275) trim(this % filename)
                return
            end if
        end do
        
        this % homogeneous_mesh = .true.
        write(*,285) trim(this % filename)
        return

    end subroutine assign_dx
    !————————————————————————————————————————————————————————————————————————————————————————
    subroutine print_info (this)
        implicit none
        type(file) :: this
        integer    :: i

225     format (" # 2-column filetype: ")          ! output style
235     format (" # name    : ", A)                ! output style
240     format (" # status  : ", A)                ! output style
245     format (" # file_ID : ", I6)               ! output style
255     format (" # ", A6, A12,   A12)   ! output style
265     format (" # ", I6, F12.8, F12.8) ! output style

        write(*,225)
        write(*,235) trim(this % filename)
        write(*,240) this % isset
        write(*,245) this % file_id

        write(*,255) "line", "x", "y"
        do i = 1, this % lines
            write(*,265) i, this % x(i), this % y(i)
        end do

    end subroutine print_info
    !————————————————————————————————————————————————————————————————————————————————————————
    !        function linear_interp (this, x) ! x is an array of real numbers
    !            implicit none
    !            type(file)       :: this
    !            real, intent (in) :: x(:)
    !            !real, optional, intent (in) :: x
    !            real, allocatable :: linear_interp(:)
    !            integer           :: i, j
    !
    !            allocate(linear_interp(1:size(x)));
    !            linear_interp = 0.0
    !            ! we need to find interval, where x has falled
    !            ! but first, we need to make sure it is not extrapolation
    !            do j = 1, size(x)
    !                if ( ( x(j) > this % x(this % lines) ) .or. ( x(j) < this % x(1) ) ) then
    !                      return ! go out of function
    !                end if
    !                do i = 1, this % lines - 1 ! check all inner intervals
    !                    if ( ( x(j) <= this % x(i+1) ) .and. ( x(j) >= this % x(i) ) ) then
    !                        linear_interp = this % y(i) + (this % y(i+1) - this % y(i))*(x(j) - this % x(i))/(this % x(i+1) - this % x(i))
    !                        return ! go out of function
    !                    end if
    !                end do
    !            end do
    !
    !        end function linear_interp

    !————————————————————————————————————————————————————————————————————————————————————————
    function linear_interp (this, x) ! x is a real number
        implicit none
        type(file)        :: this
        real, intent (in) :: x
        real              :: linear_interp
        integer           :: i

        linear_interp = 0.0
        if (.not. this % isset) return

        ! we need to find interval, where x has fallen
        ! but first, we need to make sure it is not extrapolation

        !if ( ( x > this % x(this % lines) ) .or. ( x < this % x(1) ) ) then
        !    return ! go out of function
        !end if

        i = int (x / this % dx) + 1

        if ( x > this % x(this % lines) .or. i == this % lines )  then
            linear_interp = this % y(this % lines)
            return ! go out of function
        end if
            
        if ( x < this % x(1) )  then
            linear_interp = this % y(1)
            return ! go out of function
        end if
        
        if ( this % homogeneous_mesh ) then ! if mesh if homogeneous

            linear_interp = this % y(i) + (this % y(i+1) - this % y(i))*(x - this % x(i))/(this % x(i+1) - this % x(i))
            return ! go out of function
        else
            do i = 1, this % lines - 1 ! check all inner intervals
                if ( ( x <= this % x(i+1) ) .and. ( x >= this % x(i) ) ) then
                    linear_interp = this % y(i) + (this % y(i+1) - this % y(i))*(x - this % x(i))/(this % x(i+1) - this % x(i))
                    return ! go out of function
                end if
            end do
        end if

    end function linear_interp

!        function linear_interp (this, x) ! x is a poiner to array or number
!            implicit none
!            type(file)       :: this
!            real, pointer, intent (in) :: x
!            real, allocatable :: linear_interp
!            integer           :: i
!
!            linear_interp = 0.0
!            ! we need to find interval, where x has falled
!            ! but first, we need to make sure it is not extrapolation
!            
!            if ( ( x > this % x(this % lines) ) .or. ( x < this % x(1) ) ) then
!                  return ! go out of function
!            end if
!            do i = 1, this % lines - 1 ! check all inner intervals
!                if ( ( x <= this % x(i+1) ) .and. ( x >= this % x(i) ) ) then
!                    linear_interp = this % y(i) + (this % y(i+1) - this % y(i))*(x - this % x(i))/(this % x(i+1) - this % x(i))
!                    return ! go out of function
!                end if
!            end do
!        end function linear_interp
!————————————————————————————————————————————————————————————————————————————————————————
end module two_column_file_with_x_and_y
!————————————————————————————————————————————————————————————————————————————————————————


!program module_test
!    use two_column_file_with_x_and_y
!    implicit none
!    type(file) :: a , b
!
!    call set_filename(a,"/home/palkine/Fortran/Projects/combustion/T-FlowS-comb-Problem-3-CalcPS/Problem1/tables/density_of_x.dat")
!call print_info(a)
!
!call set_filename(b,"/home/palkine/Fortran/Projects/combustion/T-FlowS-comb-Problem-3-CalcPS/Problem1/tables/viscosity_of_x.dat")
!call print_info(b)
!
!
!
!end program
!
!!————————————————————————————————————————————————————————————————————————————————————————
!!---- Just working example
!!————————————————————————————————————————————————————————————————————————————————————————
!module class_Circle
!implicit none
!real :: pi = 3.1415926535897931d0 ! class-wide private constant
!
!class, public :: Circle
!real :: radius
!    contains
!        procedure :: area => circle_area
!        procedure :: print => circle_print
!end class Circle
!contains
!    function circle_area(this) result(area)
!    class(Circle), intent(in) :: this
!    real :: area
!        area = pi * this%radius**2
!    end function circle_area
!
!    subroutine circle_print(this)
!    class(Circle), intent(in) :: this
!    real :: area
!        area = this%area()  ! Call the class-bound function
!        print *, 'Circle: r = ', this%radius, ' area = ', area
!    end subroutine circle_print
!end module class_Circle

!program circle_test
!  use class_Circle
!  implicit none
!
!  class(Circle) :: c     ! Declare a variable of type Circle.
!  c = Circle(1.5)       ! Use the implicit constructor, radius = 1.5.
!  call c%print          ! Call the type-bound subroutine
!end program circle_test
