!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!                    ________  ______                                  !
!                   |        ||      \                                 !
!                   `--.  .--'|  ,-.  \________  ___                   !
!                      |  |___|  |_/  /  _  \  \/  /                   !
!                      |  |___|      /  _____\    /                    !
!                      |  |   |  |\  \  \____/    \                    !
!                      |__|   |__| \__\_____/__/\__\                   !
!                                                                      !
!                   UNSTRUCTURED PRAGMATIC LES SOLVER                  !
!                                                                      !
!                   COMBUSTION CODE                                    !
!                                                                      !
!----------------------------------------------------------------------!
!                                                                      !
!   Pragmatic means that it can use both FRACTIONAL STEP and SIMPLE    !
!                                                                      !
!                            version: Chatou                           !
!                                                                      !
!----------------------------------------------------------------------!
!                                                                      !
!                                     Bojan NICENO                     !
!                                     Delft University of Technology   !
!                                     Faculty of Applied Sciences      !
!                                     Section Thermofluids             !
!                                     niceno@ws.tn.tudelft.nl          !
!                                                                      !
!======================================================================!
PROGRAM Processor
    !----------------------------------------------------------------------!
    !   Unstructured Finite Volume solver.                                 !
    !----------------------------------------------------------------------!
    !------------------------------[Modules]-------------------------------!
    use all_mod
    use pro_mod
    use les_mod
    use par_mod
    use rans_mod

    use Moin_Problem_mod !  module to control conflicting Moin modules
    !----------------------------------------------------------------------!
    implicit none
    !------------------------------[Calling]-------------------------------!
    real :: CorUVW
    !-------------------------------[Locals]-------------------------------!
    integer          :: i, m, n, Ndtt_temp, c
    integer          :: ini_start
    real             :: Mres, CPUtim
    real             :: start, finish
    character        :: namSav*17
    logical          :: restar, multiple, Dismiss_Regular_Save
    real             :: vol
    real,allocatable :: Rho_ini_1(:)
    real             :: Rho_error
    !------------------------------[Interface]-------------------------------!
!    interface operator (.X.) ! defines * for two solver as call to function
!      procedure multiply_unknowns_for_solver
!    end interface
    !--------------------------------[CVS]---------------------------------!
    !  $Id: Processor.f90,v 1.92 2009/01/21 10:10:51 IUS\mhadziabdic Exp $
    !  $Source: /home/IUS/mhadziabdic/.CVSROOT/T-Rex/Process/Processor.f90,v $
    !======================================================================!
 
    !-----------------------------[Interface]------------------------------!


    call cpu_time(start)
    !---- Test the precision
    open(90,FORM='UNFORMATTED',FILE='Processor.real');
    write(90) 3.1451592
    close(90)
               
    !//////////////////////////////////!
    !     Start parallel execution     !
    !//////////////////////////////////!

    call StaPar
    call Timex(CPUtim)

    !///////////////////////////////!
    !     Open the command file     !
    !///////////////////////////////!

    open(7, FILE='T-FlowS.cmn')

    if(this  < 2) call logo
    call wait

    !---- initialize parameters
    call IniPar

    !---- load the finite volume grid
    call CnsLoa
    call GeoAloc
    call GeoLoa
    call BufLoa
    call Exchng(volume(-NbC))

    call wait
    call TopolM
    call wait

    !>>>>>>>>>>>>>>>>>>>!
    !                   !
    !     TIME LOOP     !
    !                   !
    !<<<<<<<<<<<<<<<<<<<!

    Ndt  = 0
    Ndtt = 0

    call LoaRes(restar)

    if(restar) then
      call BouLoa(.false.)
    end if

    !----- Read command file (T-FlowS.cmn)
    call ReaCom(restar)

    !----- Initialize variables
    if(.not. restar) call BouLoa(.true.)

    if(.not. restar) call IniVar ! call initial variables and moin functions
    if (Moin_Problem().eq.1 .or. Moin_Problem().eq.10) then
        allocate (Rho_ini_1(-NbC:NC)); Rho_ini_1 = 0.0 ! Moin exit criterium
    end if

    call Wait

    write(*,*) "Ndtt=", Ndtt
    if(restar) then
        if (ini.ge.Nini) then ! loaded saves already finished outer iterations
            ini = 1
            Ndtt = Ndtt + 1 
        end if 
        !continue from loaded ini
        Ndtt = Ndtt - 1 
        ini_start = ini ! value from LoaRes
        Dismiss_Regular_Save=.true.
    else
        ini_start = 1
    end if

    !----- Interpolate between diff. meshes
    call LoaIni()

    !----- Check if there are more materials
    multiple = .false.
    i = StateMat(1)
    do m = 1, Nmat
        if(StateMat(m) /= i) multiple = .true.
    end do

    if(this  < 2)             &
        write(*,'(A18,15A11)')  &
        '#        N        ',   &
        ' Mass inb. ',          &
        '    U:     ',          &
        '    V:     ',          &
        '    W:     ',          &
        '    P:     '

    if(this  < 2)             &
        write(*,'(A18,15A11)')  &
        '  CFL max: ',          &
        '  Pe max:  ',          &
        '  Flux x:  ',          &
        '  Flux y:  ',          &
        '  Flux z:  ',          &
        'iterations ',          &
        '  dp/dx:   ',          &
        '    K:     '
    call wait

    !----- Prepare ...
    call Calc3()
    call FindBad()
    if(SIMULA==LES.and.MODE==SMAG.and..NOT.restar) call NearWallCell()

    !----- Prepare the gradient matrix for velocities
    call CalcG(.true.)

    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
    !        LET THE TIME        !
    !     INTEGRATION  BEGIN     !
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<!


    LinMon0 = ' '
    LinMon1 = '# [1]'
    LinMon2 = '# [2]'
    LinMon3 = '# [3]'
    LineRes = '#'
    LineRes1= ' '

    !---- Print the areas of monitoring planes
    if(This < 2) then
        do m = 1, Nmat
            write(*,'(A5,I2,A2,1PE12.3)') '# Ax(',m,')=', AreaX(m)
            write(*,'(A5,I2,A2,1PE12.3)') '# Ay(',m,')=', AreaY(m)
            write(*,'(A5,I2,A2,1PE12.3)') '# Az(',m,')=', AreaZ(m)
        end do
    end if

    if (this < 2) write(*,*) '# VISc = ', VISc, ' VIScZmix = ', VIScZmix

    if(Ndt == 0)  goto 6

    !---------------------------------------------------
    !       START TIME SERIES
    !---------------------------------------------------
    do n = Ndtt+1, Ndt + 1 !!! should be without 1

        Time = Time + dt


        if(Cm(1) > 0) then
            write(LinMon0(1: 6), '(I6)')      n;
            write(LinMon0(7:18), '(1PE12.3)') Time;
        end if

        if(SIMULA == LES) then
            call CalcShear(U % n, V % n, W % n, Shear)
            if(MODE == DYN)  call Calc_Sgs_Coefficients_Dynamic_Smagorinsky()
        end if

        !call CalcConvect

        do ini = ini_start, Nini ! start outer iterations for FRACTION & SIMPLE

            if ( ini_start > 1 ) ini_start = 1

            !---Moin Error
            if (ini == 1 .and. n > 1 .and. MoinID .ne. 0 ) call Moin_Print_Error(Moin_Time())
            
            !---Moin Time
            if (ini == 1) call Moin_Set_Time(Time)

            if (ini==1) then
                if (n.eq.Ndt+1) then ! EXIT
                    goto 3
                end if
            end if ! ini == 1

            if (Moin_Problem().eq.1 .or. Moin_Problem().eq.10) Rho_ini_1 = Rho % n ! these 2 problems have a special exit option

!            !--------------------- calculating Z_SRC coefficient
!            ! Z_SRC_of_Z fixed at current timestep
!            select case (Moin_Problem())
!                case (10)
!                    if (ini == 1) then
!                        Z_SRC_INT   = 0.0
!                        vol         = 0.0
!
!                        do c=1,NC
!                            Z_SRC_INT = Z_SRC_INT + Moin_Problem_z_src (Zmix % n(c), 0.0, 0.0)*volume(c)
!                            vol  = vol  + volume(c)
!                        end do
!
!                        call GloSum(vol)
!                        call GloSum(Z_SRC_INT)
!
!                        Z_SRC_INT= vol/Z_SRC_INT/2.0
!                        !Z_SRC_INT= 60.563188008109336 !exact value
!
!                        if (this < 2) write(*,*) "Z_coef=", Z_SRC_INT
!
!                        !do c=1,NC
!                        !    Zmixr % source(c) = Z_SRC_INT * Moin_Problem_z_src (Zmix % n(c), 0.0, 0.0) ! Z_SRC_of_Z
!                        !end do
!                    end if
!            end select

            call CalcDens(Rho, .true., Rho0(), Rho1(), Moin_Problem()) ! estimate Rho on ini==1, update other variables, apply b.c.

            ! get new source term for Mass Fraction Equation
            if ( SIMULA == LES ) call Mass_Fraction_Equation_Source(Zmix % n)

            !calculate Zmix gradients
            call GraPhi(Zmix % n, 1, Zmixx,.true.) ! dZ/dx
            call GraPhi(Zmix % n, 2, Zmixy,.true.) ! dZ/dy
            call GraPhi(Zmix % n, 3, Zmixz,.true.) ! dZ/dz

            call CalcZmix( 7    , & ! = 7
                           Zmixx, & ! dz/dx
                           Zmixy, & ! dz/dy
                           Zmixz  ) ! dz/dz

            !----------------- equation of state
            call EOS(Rho0(), Rho1())

            !!calculate pressure gradients
            if(.NOT. multiple) then
                call GradP (P % n,Px,Py,Pz) ! dP/dx, dP/dy, dP/dz
            else
                call GradP3(P % n,Px,Py,Pz) ! dP/dx, dP/dy, dP/dz
            end if

            !calculate velocity gradients
            call GraPhi(U % n , 1, Ux ,.true.)    ! dU/dx
            call GraPhi(U % n , 2, Uy ,.true.)    ! dU/dy
            call GraPhi(U % n , 3, Uz ,.true.)    ! dU/dz

            call GraPhi(V % n , 1, Vx ,.true.)    ! dV/dx
            call GraPhi(V % n , 2, Vy ,.true.)    ! dV/dy
            call GraPhi(V % n , 3, Vz ,.true.)    ! dV/dz

            call GraPhi(W % n , 1, Wx ,.true.)    ! dW/dx
            call GraPhi(W % n , 2, Wy ,.true.)    ! dW/dy
            call GraPhi(W % n , 3, Wz ,.true.)    ! dW/dz

            !---- Ur momentum component -----------------------!
            call NewUVW( 1,              & !
                         U,              & !
                         Ux,   Uy,   Uz, & !
                         Vx,   Wx,       & !
                         Vy,   Wz,       & !
                         Sx,   Sy,   Sz, & !
                         Dx,   Dy,   Dz, & !
                         Px )              !
            !---- Vr momentum component -----------------------!
            call NewUVW( 2,              & !
                         V,              & !
                         Vy,   Vx,   Vz, & !
                         Uy,   Wy,       & !
                         Ux,   Wz,       & !
                         Sy,   Sx,   Sz, & !
                         Dy,   Dx,   Dz, & !
                         Py )              !
            !---- Wr momentum component -----------------------!
            call NewUVW( 3,              & !
                         W,              & !
                         Wz,   Wx,   Wy, & !
                         Uz,   Vz,       & !
                         Ux,   Vy,       & !
                         Sz,   Sx,   Sy, & !
                         Dz,   Dx,   Dy, & !
                         Pz )              !

            !seem useless because of ModOut, but not for ModUVW
            call CalcDens(Rho, .false., Rho0(), Rho1(), Moin_Problem())  ! update other variables

            if(ALGOR == SIMPLE) then
                call Exchng(Asave)
                if (Moin_Problem() .ne. 3)  call ModOut()
                call CalcPS()
            end if


            ! grad pp is required in CorUVW according to Ferziger
            if(.NOT. multiple) then
                call GradP(PP % n,Px,Py,Pz)
            else
                call GradP3(PP % n,Px,Py,Pz)
            end if



            !call CalcFlux               !  NOT RESTED
            Mres = CorUVW()              !  project the velocities

            call CalcDens(Rho, .false., Rho0(), Rho1(), Moin_Problem())  ! update other variables

            write(LineRes(2:6),'(I5)') ini

            !----- Errors output (LinRes - input error, LinRes1 - output error)
            if(Cm(1) /= 0) write(*,'(A100)') LineRes
            !if(Cm(1) /= 0) write(*,'(A100)') LineRes1



            !---- TOLERANCE CHECK
            if(ALGOR == SIMPLE) then
              if ( Moin_Problem() == 1 .or. &
                   Moin_Problem() == 10 ) then
                Rho_error = 0.0
                call wait
                Rho_error = max(maxval( abs( Rho % n - Rho_ini_1 )), Rho_error)
                if (Rho_error < 1.e-8) then
                  if (this < 2) write(*,*) "End of outer iteraions: |rho(ini)-rho(ini-1)|=",  Rho_error,  "< 10^-8"
                  exit
                end if
              elseif( res(1) <= SIMTol .and. & ! Ur   
                      res(2) <= SIMTol .and. & ! Vr
                      res(3) <= SIMTol .and. & ! Wr
                      res(4) <= SIMTol .and. & ! P
                      res(7) <= SIMTol &       ! Zmixr
                ) then
                exit ! end outer iterations
              end if
            end if ! SIMPLE

            !----- End of the current iteration
        end do ! do ini =1,Nini


        !----- End of the current time step n
        if(Cm(1) /= 0 ) then
            write(*,'(A138)') LinMon0
            do m = 1, Nmat
                if(m .eq.1) write(*,'(A138)') LinMon1
                if(m .eq.2) write(*,'(A138)') LinMon2
                if(m .eq.3) write(*,'(A138)') LinMon3
            end do
        else
            write(*,*) "No Monitoring points found, though Nmon=", Nmon
        end if

        !----- Write the values in monitoring points
        do i = 1, Nmon
            if(Cm(i) > 0) then
                if(HOT == NO) then
                    write(10+i,'(I9,5E16.6)')                    &
                        n, U % n(Cm(i)), V%n(Cm(i)), W%n(Cm(i)), P%n(Cm(i)), Cdyn(Cm(i))
                else
                    write(10+i,'(I9,5E16.6)')                    &
                        n, U % n(Cm(i)), V%n(Cm(i)), W%n(Cm(i)), P%n(Cm(i)), Zmix%n(Cm(i))
                end if
            end if
        end do


        call CalcMn(Nstat, n)  !  calculate mean values
        call wait

        !-----------------------------------------------------!
        !                                                     !
        !   Recalculate the pressure drop                     !
        !   to keep the constant mass flux                    !
        !                                                     !
        !   First Newtons law:                                !
        !   ~~~~~~~~~~~~~~~~~~                                !
        !   F = m * a                                         !
        !                                                     !
        !   where:                                            !
        !   ~~~~~~                                            !
        !   a = dv / dt = dFlux / dt * 1 / (A * rho)          !
        !   m = rho * V                                       !
        !   F = Pdrop * l * A = Pdrop * V                     !
        !                                                     !
        !   finally:                                          !
        !   ~~~~~~~~                                          !
        !   Pdrop * V = rho * V * dFlux / dt * 1 / (A * rho)  !
        !                                                     !
        !   after cancelling: V and rho, it yields:           !
        !                                                     !
        !   Pdrop = dFlux/dt/A                                !
        !                                                     !
        !-----------------------------------------------------!
        !    do m = 1, Nmat
        !      if(SHAKE(m) == YES) then
        !        if( n  < SHAKE_PER(m) .and. MOD(n+1,SHAKE_INT(m)) == 0 ) then
        !          call UserPerturb2(5.0,n,m)
        !        endif
        !      endif

        !      if( FLUXoX(m)  /=  0.0 ) then                                  ! WRONG??
        !!        PdropX(m) = (FLUXoX(m)-FLUXx(m)) / (dt*AreaX(m)+TINY)
        !        PdropX(m) = (2.0-FLUXx(m)) / (dt*AreaX(m)+TINY)
        !      end if
        !      if( FLUXoY(m)  /=  0.0 ) then
        !!        PdropY(m) = (FLUXoY(m)-FLUXy(m)) / (dt*AreaY(m)+TINY)
        !        PdropY(m) = (0.0-FLUXy(m)) / (dt*AreaY(m)+TINY)
        !      end if
        !      if( FLUXoZ(m)  /=  0.0 ) then
        !!        PdropZ(m) = (FLUXoZ(m)-FLUXz(m)) / (dt*AreaZ(m)+TINY)
        !        PdropZ(m) = (0.0-FLUXz(m)) / (dt*AreaZ(m)+TINY)
        !      end if
        !    end do

        !---- Regular savings, each 1000 time steps
        if(mod(n,5000) == 0 .and. .not.Dismiss_Regular_Save) then
            Ndtt_temp = Ndtt
            Ndtt      = n
            namSav = 'SAVExxxxxx'
            write(namSav(5:10),'(I6.6)') n
            call SavParView(this,NC,namSav)
            call SavRes(namSav)
            Ndtt = Ndtt_temp
        end if
        Dismiss_Regular_Save = .false.

    ! !---- Is user forcing it to stop ?
    !     open(1, file='exit_now', err=7, status='old')
    !       call Wait
    !       if(This < 2) close(1, status='delete')
    !       Ndtt = n
    !       goto 6
    continue 

    !---- Is user forcing it to save ?
    open(1, file='save_now', err=8, status='old') 
    call Wait
    if(This < 2) close(1, status='delete')
    Ndtt_temp = Ndtt
    Ndtt      = n
    namSav = 'SAVExxxxxx'
    write(namSav(5:10),'(I6.6)') n
    call SavRes(namSav)
    !call DatSav(namSav)
    !call ProSav(namSav)
    call SavParView(this,NC,namSav)
    !if(CHANNEL == YES) then
    !  call UserCutLines_channel(zc)
    !else if(PIPE == YES) then
    !  call UserCutLines_pipe
    !end if
    Ndtt = Ndtt_temp
8 continue

  end do                    ! n, number of time steps  

  if(This < 2) then
      open(9,FILE='stop')
      close(9)
  end if

  Ndtt = n - 1                                 

  !-----------------------------------------!
  !     Save the results after finishing    !
  !-----------------------------------------!

  !6 call UserCutLines_Horiz_LES(namSav)

6 call SavRes

  !  call SavIni
  !  call ProSav
  !  call DatSav
  !namSav = 'SAVExxxxxx'
  !write(namSav(5:10),'(I6.6)') n
  !call SavParView(this,NC,namSav)
  !call Wait
  
 
  !  call UserDiffuser

3 if(this  < 2) write(*,*) '# Exiting !'

  !////////////////////////////////!
  !     Close the command file     !
  !////////////////////////////////!
  close(7)                


  !////////////////////////////////////!
  !     Close the monitoring files     !
  !////////////////////////////////////!
  do n=1,Nmon
      if(Cm(n)  > 0) close(10+n)
  end do

  !////////////////////////////////!
  !     End parallel execution     !
  !////////////////////////////////!
  !call SaveVTK_VTKFortran_lib(this,NC,'SAVEffffff')

  call endpar            
  call cpu_time(finish)
  if(this < 2) print '("Time = ",f14.3," seconds.")',finish-start

  END PROGRAM Processor  
