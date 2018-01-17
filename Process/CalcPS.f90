!======================================================================!
subroutine CalcPS()
    !----------------------------------------------------------------------!
    !   Forms and solves pressure equation for the S.I.M.P.L.E. method.    !
    !   par. 8.8. Ferziger, Peric                                          !
    !----------------------------------------------------------------------!
    !------------------------------[Modules]-------------------------------!
    use allp_mod, only: buffer, convect, inflow, outflow, pressure
    use all_mod
    use pro_mod
    use par_mod
    use Moin_Problem_mod
    !----------------------------------------------------------------------!
    implicit none
    !-------------------------------[Parameters]---------------------------!

    !-------------------------------[Locals]-------------------------------!
    integer :: s, c, c1, c2, niter, m
    real    :: Urs, Vrs, Wrs, A12
    real    :: error
    real    :: SMDPN
    real    :: dPxi, dPyi, dPzi, bc_copy(1:NC)
    real    :: bc_sum, vol

    !--------------------------------[CVS]---------------------------------!
    !  $Id: CalcPS.f90,v 1.23 2008/12/05 13:25:42 IUS\mhadziabdic Exp $
    !  $Source: /home/IUS/mhadziabdic/.CVSROOT/T-Rex/Process/CalcPS.f90,v $
    !----------------------------------------------------------------------!
    !
    !  The form of equations which I am solving:
    !
    !     /               /
    !    |               |
    !    | rho u dS = dt | GRAD pp dS
    !    |               |
    !   /               /
    !
    !  Dimension of the system under consideration
    !
    !     [App] {pp} = {bpp}               [kg/s]
    !
    !  Dimensions of certain variables
    !
    !     Flux           [kg/s]
    !     A_(rhou)_i     [m^3/s]
    !     P = PP         [kg/ms^2] = [N]
    !     A12 = A_PP     [ms]
    !     b              [kg/s] = [flux] = [d rho / dt *dV]
    !
    !======================================================================!
    Aval(:) = 0.0
    !------------------------------------------!
    !      Initialize the source term for      !
    !     the pressure correction equation     !
    !------------------------------------------!
    b(:) = 0.0
    bc_copy(:) = 0.0


    !---------------------------------------------!
    !     Initialize the pressure corrections     !
    !---------------------------------------------!
    PP % n = 0.0

    !-----------------------------------------------------!
    !     Calculate the mass fluxes on the cell faces     !
    !-----------------------------------------------------!

    do m = 1, Nmat
        do s = 1, NS
            c1 = SideC(1,s)
            c2 = SideC(2,s)

            !---- side is inside the domain
            if(c2  > 0 .or. (c2  < 0 .and. TypeBC(c2) == BUFFER) ) then

                SMDPN = ( Sx(s)*Sx(s) + Sy(s)*Sy(s) + Sz(s)*Sz(s) ) &
                    / ( Sx(s)*Dx(s) + Sy(s)*Dy(s) + Sz(s)*Dz(s) )

                !---- interpoliraj gustocu i brzine
                Urs  = f(s) * Rho % n(c1)*U % n(c1) + (1.0-f(s)) * Rho % n(c2)*U % n(c2)
                Vrs  = f(s) * Rho % n(c1)*V % n(c1) + (1.0-f(s)) * Rho % n(c2)*V % n(c2)
                Wrs  = f(s) * Rho % n(c1)*W % n(c1) + (1.0-f(s)) * Rho % n(c2)*W % n(c2)

                !---- calculate coeficients for the system matrix
                !---- equation (8.59) first part before ~~ in Ferziger, Peric
                if(c2  > 0) then
                    A12 = 0.5 * SMDPN *           &
                        (  volume(c1) / Asave(c1) &
                        + volume(c2) / Asave(c2) )
                    Aval(SidAij(1,s))  = - A12 !+Aval(SidAij(1,s))
                    Aval(SidAij(2,s))  = - A12 !+Aval(SidAij(2,s))
                    Aval(Adia(c1)) = Aval(Adia(c1)) + A12
                    Aval(Adia(c2)) = Aval(Adia(c2)) + A12

                else ! BUFFER, volume(c2<0)>0
                    A12 = 0.5 * SMDPN *          &
                        ( volume(c1) / Asave(c1) &
                        + volume(c2) / Asave(c2) )
                    Aval(Adia(c1)) = Aval(Adia(c1)) +  A12
                    Abou(c2)  = - A12
                end if

                !---- interpoliraj razliku pritiska
                dPxi=.5*( Px(c1) + Px(c2) )*Dx(s)
                dPyi=.5*( Py(c1) + Py(c2) )*Dy(s)
                dPzi=.5*( Pz(c1) + Pz(c2) )*Dz(s)

                !---- now calculate the flux through cell face
                !---- equation (8.56)*(rhoui)*S in Ferziger, Peric
                Flux(s) = ( Urs*Sx(s) + Vrs*Sy(s) + Wrs*Sz(s) ) &
                    + A12 * (P % n(c1) - P % n(c2))             &
                    + A12 * (dPxi + dPyi + dPzi)
                !UDS flux
                b(c1) = b(c1) - Flux(s)
                if(c2  > 0) b(c2) = b(c2) + Flux(s)

            !----- side is on the boundary
            else ! (c2 < 0)

                if(TypeBC(c2) == INFLOW) then
                    !Aval(Adia(c1)) = Aval(Adia(c1)) -(1.0 * SMDPN) !* dt
        
                    Urs  = f(s) * Rho % n(c1)*U % n(c1) + (1.0-f(s)) * Rho % n(c2)*U % n(c2)
                    Vrs  = f(s) * Rho % n(c1)*V % n(c1) + (1.0-f(s)) * Rho % n(c2)*V % n(c2)
                    Wrs  = f(s) * Rho % n(c1)*W % n(c1) + (1.0-f(s)) * Rho % n(c2)*W % n(c2)

                    !! pin P to 0 on INFLOW. If commented, then P'|inf = 0
                    SMDPN = ( Sx(s)*Sx(s) + Sy(s)*Sy(s) + Sz(s)*Sz(s) ) &
                          / ( Sx(s)*Dx(s) + Sy(s)*Dy(s) + Sz(s)*Dz(s) )
                    A12 = SMDPN * volume(c1) / Asave(c1)
                    !Aval(Adia(c1)) = Aval(Adia(c1)) +  A12

                    Flux(s) = ( Urs*Sx(s) + Vrs*Sy(s) + Wrs*Sz(s) )
                    b(c1) = b(c1) - Flux(s) ! + A12 * PP % n(c2)
                else if(TypeBC(c2) == OUTFLOW .or.  TypeBC(c2) == CONVECT) then
                    !Aval(Adia(c1)) = Aval(Adia(c1)) - (1.0 * SMDPN) !* dt

                    Urs  = f(s) * Rho % n(c1)*U % n(c1) + (1.0-f(s)) * Rho % n(c2)*U % n(c2)
                    Vrs  = f(s) * Rho % n(c1)*V % n(c1) + (1.0-f(s)) * Rho % n(c2)*V % n(c2)
                    Wrs  = f(s) * Rho % n(c1)*W % n(c1) + (1.0-f(s)) * Rho % n(c2)*W % n(c2)

                    !! pin P to 0 on OUTFLOW. If commented, then P'|out = 0
                    SMDPN = ( Sx(s)*Sx(s) + Sy(s)*Sy(s) + Sz(s)*Sz(s) ) &
                        / ( Sx(s)*Dx(s) + Sy(s)*Dy(s) + Sz(s)*Dz(s) )
                    A12 = SMDPN * volume(c1) / Asave(c1)
                    Aval(Adia(c1)) = Aval(Adia(c1)) +  A12

                    Flux(s) = ( Urs*Sx(s) + Vrs*Sy(s) + Wrs*Sz(s) )
                    b(c1)   = b(c1) - Flux(s) ! + A12 * PP % n(c2)
                else if(TypeBC(c2) == PRESSURE) then
                    Urs  = Rho % n(c1)*U % n(c1)
                    Vrs  = Rho % n(c1)*V % n(c1)
                    Wrs  = Rho % n(c1)*W % n(c1)
                    Flux(s) = ( Urs*Sx(s) + Vrs*Sy(s) + Wrs*Sz(s) )
                    b(c1) = b(c1) - Flux(s)
                    SMDPN = ( Sx(s)*Sx(s) + Sy(s)*Sy(s) + Sz(s)*Sz(s) ) &
                        / ( Sx(s)*Dx(s) + Sy(s)*Dy(s) + Sz(s)*Dz(s) )
                    A12 = SMDPN * volume(c1) / Asave(c1)
                    Aval(Adia(c1)) = Aval(Adia(c1)) +  A12
                else  ! it is SYMMETRY
                    Flux(s) = 0.0
                end if
            end if

        end do
    end do


    ! d (rho^n+1) /dt
    errmax = 0.0
    do c = 1, NC
        errmax = max(errmax, abs(b(c)))
    end do

    !-----------------------------------------------------------!
    !     Do not solve the pressure corection too accurate.      !
    !     Value 1.e-18 blows the solution.                      !
    !     Value 1.e-12 keeps the solution stable                !
    !-----------------------------------------------------------!


    bc_sum = 0.0 
    vol = 0.0 

    do m = 1, Nmat
        do c = 1, NC
            b(c) = b(c) - volume(c)/dt*(1.5*Rho % n(c) - 2.0*Rho % o(c) + 0.5*Rho % oo(c))

            bc_sum = bc_sum + b(c)*volume(c)
            vol    = vol    + volume(c) 
        end do
    end do
    call glosum(bc_sum)
    call glosum(vol)

    !    if(this < 2) write(*,*) "cont = ", bc_sum/vol

    do c = 1, NC
        b(c) = b(c) - bc_sum/vol ! earlier it was just a sum(b)/NC - showed same results in hex
        !b(c) = b(c) - bc_sum/NC ! but in wedge configuration sum(b)/NC is the only right choice (why??)
    end do

    !do s = 1, NS
    !  c1 = SideC(1,s)
    !  c2 = SideC(2,s)
    !  if (TypeBC(c2)/=BUFFER) then
    !    b(c1) = 0.0
    !  end if
    !end do

    bc_copy = b

    niter=800


    call cg(NC, NbC, NONZERO, Aval, Acol,Arow,Adia,Abou,  &
        PP % n, b, PREC, niter, PP % STol,            &
        res(4), error)

    write(LineRes(53:64),  '(1PE12.3)') res(4)
    write(LineRes(89:92),  '(I4)')      niter
    write(LineRes1(53:64), '(1PE12.3)') error 


    !-----------------------------------!
    !     Update the pressure field     !
    !----------------------------------`-!
    do c = 1, NC
        P % n(c)  =  P % n(c)  +  P  % URF * PP % n(c)
    end do
    !--------------------------------------!


    !     Normalize the pressure field     !
    !--------------------------------------!

    !call Exchng(P % n)
    !call Wait
    call Exchng(PP % n)


    return

end subroutine CalcPS

    !!-------matrix visualization
    !do c = 1, NC
    !  if ( c.ne.0 .AND. zc(c)<=1.0 .AND.  zc(c)>=-1.0 ) then
    !    write(*,('(A,I10,A,F7.3,A,F7.3,A,I2)')) 'Cell:', c, ' x=', xc(c), ' y=', yc(c), &
    !    ' Width: ', Acol(c+1)-Acol(c)
    !    write(*,'(A,I10)'),'diag:',  Adia(c)
    !    write(*,'(A,25I10)') 'Acol:', ( s, s=Acol(c),Acol(c+1)-1 )
    !    !write(*,'(25I10)')   ( Arow(s), s=Acol(c),Acol(c+1)-1 )
    !    write(*,'(A,25F10.4)') '   x:', ( xc(Arow(s)), s=Acol(c),Acol(c+1)-1 )
    !    write(*,'(A,25F10.4)') '   y:', ( yc(Arow(s)), s=Acol(c),Acol(c+1)-1 )
    !    write(*,'(A,25F10.4)') '   z:', ( zc(Arow(s)), s=Acol(c),Acol(c+1)-1 )
    !    write(*,'(A,25F10.4)') 'Aval:', ( Aval(s), s=Acol(c),Acol(c+1)-1 )
    !    write(*,*) '- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -'
    !  end if
    !end do
