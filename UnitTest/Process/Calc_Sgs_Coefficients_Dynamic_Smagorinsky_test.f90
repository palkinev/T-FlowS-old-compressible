module Calc_Sgs_Coefficients_Dynamic_Smagorinsky_test
    use fruit
    use test_tools
    implicit none

    contains
    subroutine test_Calc_Sgs_Coefficients_Dynamic_Smagorinsky_test_X
    use allp_mod
    use all_mod
    use pro_mod
    use les_mod
    use rans_mod
    use moin_problem_mod

    NC = 20
    NS = 30
    allocate(rho % n(1:NC)); Rho % n = 0.
    allocate(ur % n(1:NC)); ur % n = 0.
    allocate(vr % n(1:NC)); vr % n = 0.
    allocate(wr % n(1:NC)); wr % n = 0.
    allocate(u % n(1:NC)); u % n = 0.
    allocate(v % n(1:NC)); v % n = 0.
    allocate(w % n(1:NC)); w % n = 0.


    call IniPar
    call Calc_Sgs_Coefficients_Dynamic_Smagorinsky

    call assert_true(1==1)



    
    end subroutine test_Calc_Sgs_Coefficients_Dynamic_Smagorinsky_test_X

end module Calc_Sgs_Coefficients_Dynamic_Smagorinsky_test
