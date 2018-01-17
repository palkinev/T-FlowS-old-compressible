module all_mod_test
    use fruit
    use test_tools
    implicit none

    contains
    subroutine test_all_mod_xc
    use all_mod, only : xc

    call test_real_allocatable_array(xc)
    end subroutine test_all_mod_xc

end module all_mod_test
