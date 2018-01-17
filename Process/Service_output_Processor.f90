SUBROUTINE Service_output_Processor(n,str_in)  
    USE all_mod
    USE pro_mod
    USE par_mod
    USE Moin_Problem_mod

                !------------------------ Service output ------------------
                !if ( .not.restar  .and. ini==1 .and.  ( n==1  .or. n==Ndt .or. mod(n,800)==0) ) then
                !    call Service_output_Processor(n,"pp2_1.dat")
                !end if

                !if ( .not.restar  .and. ini==2 .and.  ( n==1  .or. n==Ndt .or. mod(n,800)==0) ) then
                !    call Service_output_Processor(n,"pp2_2.dat")
                !end if
                !------------------------ Service output ------------------

    IMPLICIT NONE

    CHARACTER*(*) :: str_in
    INTEGER :: n, c
    !INTEGER :: s, c1, c2
    CHARACTER :: str*10

    str = adjustl(trim(str_in))


    if(this < 2) open(991, FILE=adjustl(trim(str_in)), status='replace', access="sequential", position='append', action='write')

    !P %mean(:) = 0.0 

    !do s=1,NS
    !  c1=SideC(1,s)
    !  c2=SideC(2,s)
    !  if(c2  > 0 .or. c2  < 0 .and. TypeBC(c2) == BUFFER) then
    !                P %mean(c1) = P %mean(c1)-Flux(s)
    !    if(c2  > 0) P %mean(c2) = P %mean(c2)+Flux(s)
    !  else
    !                P %mean(c1) = P %mean(c1)-Flux(s)
    !  end if
    !end do

    !do c=1,NC      !!- d (rho^n+1) /dt
    !  !! 1)
    !  P %mean(c) = - P %mean(c) + volume(c)/dt*(1.5*Rho % n(c) - 2.0*Rho % o(c) + 0.5*Rho % oo(c) )
    !  !! 2)
    !  !P % mean(c) = P % mean(c) - volume(c)/dt*(1.5*Rho % n(c) - 2.0*Rho % o(c) )
    !  !P % mean(c) = P % mean(c) + volume(c)/dt*( 0.5*Rho % oo(c) )
    !end do

    ! P % mean is Z_SOURCE now
    do c=1,NC
        P % mean(c) = Zmix % source(c)
    end do

    ! P % mean is VISC_DYN
    do c=-NbC,NC
        P % mean(c) = VISc_Dyn(c) * 224.28074617774
    end do

    do c=-NbC,-1
        !if ( c .ne. 0 .AND. zc(c)>=-1.0 .AND. zc(c)<=1.0 .AND.  yc(c)<=0.125+1.0/200 .AND. yc(c)>=0.125-1.0/200 ) then
        if ( c .ne. 0 .and. xc(c) > 1 .AND. zc(c)>=-1.0 .AND. zc(c)<=1.0 .AND.  yc(c)<=1.0 .AND. yc(c)>=-1.0 ) then !problem 1
            !if ( c .ne. 0 .AND. zc(c)>=-0.5 .AND. zc(c)<=0.5) then
            !write(991,'(13ES26.16E3)') xc(c), yc(c), U % n(c), V % n(c), W % n(c), Zmix % n(c), &
            !Rho % n(c), Ur % n(c), Vr % n(c), Wr % n(c), Zmixr % n(c), P % n(c), P % mean(c)
            write(991,'(15ES26.16E3)') xc(c), yc(c), U % n(c), V % n(c), W % n(c), & 
                Moin_Problem_u (xc(c), yc(c), (n  )*dt),  Moin_Problem_u (xc(c), yc(c), (n-1  )*dt), &
                P % n(c), Moin_Problem_p (xc(c), yc(c), (n  )*dt), Px(c), &
                Zmix % n(c), Moin_Problem_z (xc(c), yc(c), (n  )*dt), &
                Rho % n(c), Moin_Problem_rho (xc(c), yc(c), (n  )*dt), P % mean(c)
        end if
    end do


    do c=1,NC
        !if ( c .ne. 0 .AND. zc(c)>=-1.0 .AND. zc(c)<=1.0 .AND.  yc(c)<=0.125+1.0/200 .AND. yc(c)>=0.125-1.0/200 ) then
        if ( c .ne. 0 .AND. zc(c)>=-1.0 .AND. zc(c)<=1.0 .AND.  yc(c)<=1.0 .AND. yc(c)>=-1.0 ) then !problem 1
            !if ( c .ne. 0 .AND. zc(c)>=-0.5 .AND. zc(c)<=0.5) then
            !write(991,'(13ES26.16E3)') xc(c), yc(c), U % n(c), V % n(c), W % n(c), Zmix % n(c), &
            !Rho % n(c), Ur % n(c), Vr % n(c), Wr % n(c), Zmixr % n(c), P % n(c), P % mean(c)
            write(991,'(15ES26.16E3)') xc(c), yc(c), U % n(c), V % n(c), W % n(c), & 
                Moin_Problem_u (xc(c), yc(c), (n  )*dt),  Moin_Problem_u (xc(c), yc(c), (n-1  )*dt), &
                P % n(c), Moin_Problem_p (xc(c), yc(c), (n  )*dt), Px(c), &
                Zmix % n(c), Moin_Problem_z (xc(c), yc(c), (n  )*dt), &
                Rho % n(c), Moin_Problem_rho (xc(c), yc(c), (n  )*dt), P % mean(c)
        end if
    end do

    if(this<2) write(*,*) 'Saving to ', adjustl(trim(str_in))
    do c=-NbC,-1
        !if ( c .ne. 0 .AND. zc(c)>=-1.0 .AND. zc(c)<=1.0 .AND.  yc(c)<=0.125+1.0/200 .AND. yc(c)>=0.125-1.0/200 ) then
        if ( c .ne. 0 .and. xc(c) < 1.AND. zc(c)>=-1.0 .AND. zc(c)<=1.0 .AND.  yc(c)<=1.0 .AND. yc(c)>=-1.0 ) then !problem 1
            !if ( c .ne. 0 .AND. zc(c)>=-0.5 .AND. zc(c)<=0.5) then
            !write(991,'(13ES26.16E3)') xc(c), yc(c), U % n(c), V % n(c), W % n(c), Zmix % n(c), &
            !Rho % n(c), Ur % n(c), Vr % n(c), Wr % n(c), Zmixr % n(c), P % n(c), P % mean(c)
            write(991,'(15ES26.16E3)') xc(c), yc(c), U % n(c), V % n(c), W % n(c), & 
                Moin_Problem_u (xc(c), yc(c), (n  )*dt),  Moin_Problem_u (xc(c), yc(c), (n-1  )*dt), &
                P % n(c), Moin_Problem_p (xc(c), yc(c), (n  )*dt), Px(c), &
                Zmix % n(c), Moin_Problem_z (xc(c), yc(c), (n  )*dt), &
                Rho % n(c), Moin_Problem_rho (xc(c), yc(c), (n  )*dt), P % mean(c)
        end if
    end do

    call wait 
    !call exit(1)
    if(this < 2) close(991)

    RETURN

END SUBROUTINE Service_output_Processor
