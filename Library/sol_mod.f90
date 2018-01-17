MODULE sol_mod 

    USE allp_mod

    IMPLICIT NONE

    REAL,ALLOCATABLE, PUBLIC :: D(:)
    REAL,ALLOCATABLE, PUBLIC :: p1(:), p2(:)
    REAL,ALLOCATABLE, PUBLIC :: q1(:), q2(:), r2(:)
    REAL,ALLOCATABLE, PUBLIC :: u1(:), u2(:), v1(:), v2(:)
    REAL, PUBLIC :: alfa, beta, rho, rhoold, bnrm2, sum1, sum2, error

    INTEGER, PUBLIC :: i, j, k, iter, sub
END MODULE 
