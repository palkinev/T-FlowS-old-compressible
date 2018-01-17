module allp_mod_test
    use fruit
    use test_tools
    implicit none

contains
    subroutine test_allp_mod_maxp
    use allp_mod, only : maxp

    call assert_equals (maxp, 200)

    end subroutine test_allp_mod_maxp

    subroutine test_allp_mod_maxl
    use allp_mod, only : maxl

    call assert_equals (maxl, 1000)

    end subroutine test_allp_mod_maxl

    subroutine test_allp_mod_maxpro
    use allp_mod, only : maxpro

    call assert_equals (maxpro, 1024)

    end subroutine test_allp_mod_maxpro

    subroutine test_allp_mod_inflow
    use allp_mod, only : inflow

    call assert_equals (inflow, 1)

    end subroutine test_allp_mod_inflow

    subroutine test_allp_mod_wall
    use allp_mod, only : wall

    call assert_equals (wall, 2)

    end subroutine test_allp_mod_wall

    subroutine test_allp_mod_outflow
    use allp_mod, only : outflow

    call assert_equals (outflow, 3)

    end subroutine test_allp_mod_outflow

    subroutine test_allp_mod_symmetry
    use allp_mod, only : symmetry

    call assert_equals (symmetry, 4)

    end subroutine test_allp_mod_symmetry

    subroutine test_allp_mod_convect
    use allp_mod, only : convect

    call assert_equals (convect, 5)

    end subroutine test_allp_mod_convect

    subroutine test_allp_mod_pressure
    use allp_mod, only : pressure

    call assert_equals (pressure, 12)

    end subroutine test_allp_mod_pressure

    subroutine test_allp_mod_periodic
    use allp_mod, only : periodic

    call assert_equals (periodic, 13)

    end subroutine test_allp_mod_periodic

    subroutine test_allp_mod_buffer
    use allp_mod, only : buffer

    call assert_equals (buffer, 11)

    end subroutine test_allp_mod_buffer

    subroutine test_allp_mod_wallfl
    use allp_mod, only : wallfl

    call assert_equals (wallfl, 6)

    end subroutine test_allp_mod_wallfl

    subroutine test_allp_mod_fluid
    use allp_mod, only : fluid

    call assert_equals (fluid, 7)

    end subroutine test_allp_mod_fluid

    subroutine test_allp_mod_solid
    use allp_mod, only : solid

    call assert_equals (solid, 8)

    end subroutine test_allp_mod_solid

    subroutine test_allp_mod_huge
    use allp_mod, only : huge

    call assert_equals (huge, 1.e+30)

    end subroutine test_allp_mod_huge

    subroutine test_allp_mod_tiny
    use allp_mod, only : tiny

    call assert_equals (tiny, 1.e-64)

    end subroutine test_allp_mod_tiny

    subroutine test_allp_mod_type_unknown_field_n
    use allp_mod, only : unknown
    type(Unknown) :: var

    call test_real_pointer(var % n)

    end subroutine test_allp_mod_type_unknown_field_n

    subroutine test_allp_mod_type_unknown_field_o
    use allp_mod, only : unknown
    type(Unknown) :: var

    call test_real_pointer(var % o)

    end subroutine test_allp_mod_type_unknown_field_o

    subroutine test_allp_mod_type_unknown_field_oo
    use allp_mod, only : unknown
    type(Unknown) :: var

    call test_real_pointer(var % oo)

    end subroutine test_allp_mod_type_unknown_field_oo

    subroutine test_allp_mod_type_unknown_field_X
    use allp_mod, only : unknown
    type(Unknown) :: var

    call test_real_pointer(var % X)

    end subroutine test_allp_mod_type_unknown_field_X

    subroutine test_allp_mod_type_unknown_field_C
    use allp_mod, only : unknown
    type(Unknown) :: var

    call test_real_pointer(var % C)

    end subroutine test_allp_mod_type_unknown_field_C

    subroutine test_allp_mod_type_unknown_field_mean
    use allp_mod, only : unknown
    type(Unknown) :: var

    call test_real_pointer(var % mean)

    end subroutine test_allp_mod_type_unknown_field_mean

    subroutine test_allp_mod_type_unknown_field_fluc
    use allp_mod, only : unknown
    type(Unknown) :: var

    call test_real_pointer(var % fluc)

    end subroutine test_allp_mod_type_unknown_field_fluc

    subroutine test_allp_mod_type_unknown_field_source
    use allp_mod, only : unknown
    type(Unknown) :: var

    call test_real_pointer(var % source)

    end subroutine test_allp_mod_type_unknown_field_source

    subroutine test_allp_mod_type_unknown_field_URF
    use allp_mod, only : unknown
    type(Unknown) :: var

    var % URF = -1.1
    var % URF = var % URF + 1.0
    call assert_true ( abs(var % URF + 0.1) < 1.0D-16)

    end subroutine test_allp_mod_type_unknown_field_URF

    subroutine test_allp_mod_type_unknown_field_Stol
    use allp_mod, only : unknown
    type(Unknown) :: var

    var % Stol = -1.1
    var % Stol = var % Stol + 1.0
    call assert_true ( abs(var % Stol + 0.1) < 1.0D-16)

    end subroutine test_allp_mod_type_unknown_field_Stol

    subroutine test_allp_mod_type_unknown_field_bound
    use allp_mod, only : unknown
    type(Unknown) :: var

    call assert_true (size(var % bound) == 128)
    var % bound(128) = -1.3
    var % bound(128) = var % bound(128) + 1.0
    call assert_true ( abs(var % bound(128) + 0.3) < 1.0D-16)


    end subroutine test_allp_mod_type_unknown_field_bound

    subroutine test_allp_mod_type_unknown_field_init
    use allp_mod, only : unknown
    type(Unknown) :: var

    call assert_true (size(var % init) == 128)
    var % init(128) = -1.3
    var % init(128) = var % init(128) + 1.0
    call assert_true ( abs(var % init(128) + 0.3) < 1.0D-16)

    end subroutine test_allp_mod_type_unknown_field_init

    subroutine test_allp_mod_type_unknown_field_pro
    use allp_mod, only : unknown
    type(Unknown) :: var

    call assert_true (size(var % pro) == 11024)
    var % pro(11024) = -1.3
    var % pro(11024) = var % pro(11024) + 1.0
    call assert_true ( abs(var % pro(11024) + 0.3) < 1.0D-16)

    end subroutine test_allp_mod_type_unknown_field_pro

    subroutine test_allp_mod_type_unknown_field_Sigma
    use allp_mod, only : unknown
    type(Unknown) :: var

    call test_real_number(var % Sigma)

    end subroutine test_allp_mod_type_unknown_field_Sigma

end module allp_mod_test
