module test_tools
    use fruit
    implicit none

contains
    subroutine test_real_number(var)
        real :: var
        var = -1.1
        var = var + 1.0
        call assert_true ( abs(var + 0.1) < 1.0D-16)
    end subroutine test_real_number

    subroutine test_real_pointer(var)
        real, POINTER :: var(:)

        call assert_true (.not.associated (var))

        allocate (var (-10:20))
        call assert_true (associated (var))
        call assert_true (size(var) == 31)

        var = 0.0
        var(-10) = -1.1
        var(-10) = var(-10) + 1.0
        call assert_true ( abs(var(-10) + 0.1) < 1.0D-16)
        
        deallocate (var)
        call assert_true (.not.associated (var))

    end subroutine test_real_pointer

    subroutine test_real_allocatable_array(var)
        real, allocatable :: var(:)

        allocate (var (-10:20))
        call assert_true (size(var) == 31)

        var = 0.0
        var(-10) = -1.1
        var(-10) = var(-10) + 1.0
        call assert_true ( abs(var(-10) + 0.1) < 1.0D-16)
        
        deallocate (var)

    end subroutine test_real_allocatable_array

end module test_tools
