!======================================================================!
subroutine IniVar
	!----------------------------------------------------------------------!
	!   Initialize dependent variables                                     !
	!----------------------------------------------------------------------!
	!------------------------------[Modules]-------------------------------!
	use allp_mod
	use all_mod
	use pro_mod
	use les_mod
	use par_mod
	use rans_mod

	use Moin_Problem_mod
	!----------------------------------------------------------------------!
  implicit none
	!------------------------------[Calling]-------------------------------!
	!-------------------------------[Locals]-------------------------------!
	integer :: c, c1, c2, m, s, n
	integer :: N1, N2, N3, N4, N5, N6
	real    :: RhoS
	real             :: Rho_0, Rho_1
	character        :: path_to_tables*100
	!--------------------------------[CVS]---------------------------------!
	!  $Id: IniVar.f90,v 1.28 2009/06/30 12:09:24 IUS\mhadziabdic Exp $
	!  $Source: /home/IUS/mhadziabdic/.CVSROOT/T-Rex/Process/IniVar.f90,v $
	!----------------------------------------------------------------------!
	!                                                                      !
	!                                                                      !
	!                 +-------------+                                      !
	!                /|            /|                                      !
	!               /             / |                                      !
	!              /  |          /  |                                      !
	!        +--- +-------------+   |                                      !
	!        |    |   |         |   |                                      !
	!        |    |   +- - - - -|- -+  ----+                               !
	!        H    |  /          |  /      /                                !
	!        |    | /           | /      W                                 !
	!        |    |/            |/      /                                  !
	!        +--- +-------------+  ----+                                   !
	!                                                                      !
	!             |      L      |                                          !
	!             +-------------+                                          !
	!                                                                      !
	!                                                                      !
	!       2 * Tau * L * W = dp/dx * H * W * L  -> law                    !
	!                                                                      !
	!       2 * Tau = dp/dx * H                                            !
	!                                                                      !
	!       H = 2 * delta                                                  !
	!                                                                      !
	!       Tau = dp/dx                                                    !
	!       ===========                                                    !
	!                                                                      !
	!       uTau = (Tau/rho)^0.5       -> definition                       !
	!                                                                      !
	!       uTau = ReTau * nu * delta  -> definition                       !
	!                                                                      !
	!       rho   = 1                                                      !
	!       nu    = 1                                                      !
	!       delta = 1                                                      !
	!                                                                      !
	!       Tau = dp/dx = uTau^2.0 = ReTau^2.0                             !
	!       ==================================                             !
	!                                                                      !
	!        +                                                             !
	!       y = y * uTau / nu                                              !
	!                                                                      !
	!                                                                      !
	!    +--------------+                                                  !
	!    |   +          |                                                  !
	!    |  y    <= 11  |  ?                                               !
	!    |   min        |                                                  !
	!    +--------------+                                                  !
	!                                                                      !
	!======================================================================!

	!----- Initialize model functions (Moin problem module)
	!---- overwrite values from T-Flows.cmn and *.b


	! these are defaults for Moin problems
	! after this, user is not allowed to alter Moin_Problem_ID
	! if value is interested, Moin_Problem() function should be called

	select case (MoinID)
		case (1)
			Nini = 1000 ! max outer iterations
			VIScZmix = VISc
			Rho_0 = 20.0
			Rho_1 = 1.0

			path_to_tables = "../"
			call Moin_Problem_init( &
				MoinID, &
				Rho_0, &
				Rho_1, &
				VISc, &
				VIScZmix, &
				"", &  ! rho = EOS(Z)
				"", &  ! visc_dyn(Z)
				"", &  ! Z_Source(Z)
				"", &  ! initial u(x)
				"", &  ! initial rho(x)
				"", &  ! initial Z(x)
				"", &  ! initial VISc_Dyn(x)
				"" )
			call Moin_Problem_init(MoinID, Rho_0, Rho_1, VISc, VIScZmix)
!        case (10)
!            P    % URF = 0.4
!            Nini = 200 ! Moin parameter
!
!            VISc = 1.0 / 224.28074617774
!            VIScZmix = VISc / 0.7
!            Rho_0 = 1.0
!            Rho_1 = 0.1558951
!
!            path_to_tables = "/home/palkine/Fortran/Projects/combustion/T-FlowS-comb-Problem-3-CalcPS/Problem1/tables/"
!            call Moin_Problem_init( MoinID, Rho_0, Rho_1, VISc, VIScZmix, &
!                trim(path_to_tables) // "EOS_of_Z.dat", &                  ! rho = EOS(Z)
!                trim(path_to_tables) // "viscosity_of_Z.dat", &            ! visc_dyn(Z)
!                trim(path_to_tables) // "Z_source_of_Z.dat", &             ! Z_Source(Z)
!                trim(path_to_tables) // "velocity_of_x.dat", &             ! initial u(x)
!                trim(path_to_tables) // "density_of_x.dat", &              ! initial rho(x)
!                trim(path_to_tables) // "Z_of_x.dat", &                    ! initial Z(x)
!                trim(path_to_tables) // "viscosity_of_x.dat", &            ! initial VISc_Dyn(x)
!                trim(path_to_tables) // "pressure_of_x.dat")               ! initial pressure(x)
!            do c = -NbC, NC
!                VISc_Dyn(c)     = Moin_Problem_visc_dyn_of_x(xc(c),0.0,0.0) / 224.28074617774
!                VIScZmix_Dyn(c) = VISc_Dyn(c) / 0.7
!            end do
		case (10)
			!Ur   % URF = 0.6
			!Vr   % URF = 0.6
			!Wr   % URF = 0.6
			!Zmixr% URF = 0.6
			P    % URF = 0.1
			!Ur % STol    = 1.e-10
			!Zmixr % STol = 1.e-10
			P % STol     = 1.e-10

			Nini = 1000 ! max outer iterations
			Ndt = 800
			dt = 1.0/Ndt
			VISc = 0.03
			VIScZmix = VISc
			Rho_0 = 20.0
			Rho_1 = 1.0

			path_to_tables = "../"
			call Moin_Problem_init( &
				MoinID, &
				Rho_0, &
				Rho_1, &
				VISc, &
				VIScZmix, &
				trim(path_to_tables) // "Rho_of_Zmix.dat", &  ! rho = EOS(Z)
				!"", &  ! rho = EOS(Z)
				"", &  ! visc_dyn(Z)
				"", &  ! Z_Source(Z)
				"", &  ! initial u(x)
				"", &  ! initial rho(x)
				"", &  ! initial Z(x)
				"", &  ! initial VISc_Dyn(x)
				"" )
			call Moin_Problem_init(MoinID, Rho_0, Rho_1, VISc, VIScZmix)
		case (2)
			Nini = 20 ! Moin parameter

			VISc = 1.e-3
			VIScZmix = VISc
			Rho_0 = 20.0
			Rho_1 = 1.0

			call Moin_Problem_init(MoinID, Rho_0, Rho_1, VISc, VIScZmix)
		case (3)
			Nini = 20 ! Moin parameter

			VISc = 1.e-3
			VIScZmix = VISc
			Rho_0 = 5.0
			Rho_1 = 1.0

			call Moin_Problem_init(MoinID, Rho_0, Rho_1, VISc, VIScZmix)
		case (30)
			Nini = 20 ! Moin parameter

			VISc = 1.e-3
			VIScZmix = VISc
			Rho_0 = 5.0
			Rho_1 = 1.0
			
			path_to_tables = "../"

			call Moin_Problem_init( &
				MoinID, &
				Rho_0, &
				Rho_1, &
				VISc, &
				VIScZmix, &
				trim(path_to_tables) // "Rho_of_Zmix.dat", &  ! rho = EOS(Z)
				!"", &  ! rho = EOS(Z)
				"", &  ! visc_dyn(Z)
				"", &  ! Z_Source(Z)
				"", &  ! initial u(x)
				"", &  ! initial rho(x)
				"", &  ! initial Z(x)
				"", &  ! initial VISc_Dyn(x)
				"" )
		case (4)
			VISc = 1.e-3
			VIScZmix = 0.5e-3
			Rho_0 = 1.0
			Rho_1 = 0.1

			call Moin_Problem_init(MoinID, Rho_0, Rho_1, VISc, VIScZmix)
	end select
	!convection test
	!VISc = 1.e-6
	!diffusion test
	!VISc = 3.0
	!----- Moin module part


	do n=1,Nmat
		if(namIni(n) == '') then

		!-----Initialize variables with GMV file
		else
			!if(namIni(n)(len_trim(namIni(n))-3:len_trim(namIni(n))) == ".gmv") then
			!  call ReadGMV(namIni(n))
			!else if(namIni(n)(len_trim(namIni(n))-3:len_trim(namIni(n))) == ".neu") then
			!  call ReadFLUENT(namIni(n))
			!else
			write(*,*) '@IniVar: Unknown file type for initialization!!'
			STOP
		  !end if
		end if  !end if namIni=''
	end do   !end do n=1,Nmat


	if(MoinID .ne. 0) then
		n = 1
		if (this < 2) write(*,*) 'dt = ', dt
		!-------------Vector fields initialization (for problem 1-4)
		do c = -NbC, NC
			if (c < 0 ) then ! .AND. n < 2 ???
				!density
				Rho % n(c) = Moin_Problem_rho (xc(c), yc(c), (n-1)*dt)
				!velocity
				U  %  n(c) = Moin_Problem_u   (xc(c), yc(c), (n-1)*dt)
				V  %  n(c) = Moin_Problem_v   (xc(c), yc(c), (n-1)*dt)
				W  %  n(c) = 0.0
				!pressure
				P %  n(c)  = Moin_Problem_p   (xc(c), yc(c), (n-1)*dt)
				!Zmix
				Zmix %n(c) = Moin_Problem_z   (xc(c), yc(c), (n-1)*dt)

			elseif(c .ne. 0 ) then !c > 0 .and. n < 2
				!--------------initial parameters
				!rho
				Rho % n(c) = Moin_Problem_rho (xc(c), yc(c), (n-1)*dt)
				Rho % o(c) = Moin_Problem_rho (xc(c), yc(c), (n-2)*dt)
				Rho %oo(c) = Moin_Problem_rho (xc(c), yc(c), (n-3)*dt)
				!Zmix
				Zmix %n(c) = Moin_Problem_z (xc(c), yc(c), (n-1)*dt)
				Zmix %o(c) = Moin_Problem_z (xc(c), yc(c), (n-2)*dt)
				Zmix %oo(c)= Moin_Problem_z (xc(c), yc(c), (n-3)*dt)
				!U
				U  % n(c)  = Moin_Problem_u (xc(c), yc(c), (n-1)*dt)
				U  % o(c)  = Moin_Problem_u (xc(c), yc(c), (n-2)*dt)
				U  % oo(c) = Moin_Problem_u (xc(c), yc(c), (n-3)*dt)
				!V
				V  % n(c)  = Moin_Problem_v (xc(c), yc(c), (n-1)*dt)
				V  % o(c)  = Moin_Problem_v (xc(c), yc(c), (n-2)*dt)
				V  % oo(c) = Moin_Problem_v (xc(c), yc(c), (n-3)*dt)
				!W
				W  % n(c)   = 0.0
				W  % o(c)   = 0.0
				W  % oo(c)  = 0.0
				!P
				P % n(c) = Moin_Problem_p (xc(c), yc(c), (n-1)*dt)

			end if ! initialize the parameters
		
		end do


	end if

	!=====================================!
	!        Calculate the inflow         !
	!     and initializes the Flux(s)     !
	!     at both inflow and outflow      !
	!=====================================!
	N1 = 0
	N2 = 0
	N3 = 0
	N4 = 0
	N5 = 0
	N6 = 0
	do m=1,Nmat
		MassIn(m) = 0.0
		do s=1,NS
			c1=SideC(1,s)
			c2=SideC(2,s)


			if(c2  < 0) then
				if(TypeBC(c2)  ==  INFLOW) then
					if(material(c1) == m) MassIn(m) = MassIn(m) - Flux(s)
				endif
				if(TypeBC(c2)  ==  WALL)     N1 = N1 + 1
				if(TypeBC(c2)  ==  INFLOW)   N2 = N2 + 1
				if(TypeBC(c2)  ==  OUTFLOW)  N3 = N3 + 1
				if(TypeBC(c2)  ==  SYMMETRY) N4 = N4 + 1
				if(TypeBC(c2)  ==  WALLFL)   N5 = N5 + 1
				if(TypeBC(c2)  ==  CONVECT)  N6 = N6 + 1
			end if

			RhoS = f(s) * Rho % n(c1) + ( 1.0-f(s) ) * Rho % n(c2)
      Flux(s) =   f(s) * ( Rho % n(c1) * U % n(c1) * Sx(s) +   &
                           Rho % n(c1) * V % n(c1) * Sy(s) +   &
                           Rho % n(c1) * W % n(c1) * Sz(s) ) + &
          ( 1.0-f(s) ) * ( Rho % n(c2) * U % n(c2) * Sx(s) +   &
                           Rho % n(c2) * V % n(c2) * Sy(s) +   &
                           Rho % n(c2) * W % n(c2) * Sz(s) )

			Flux_u(s) = Flux(s) / RhoS

		end do
		call iglsum(N1)
		call iglsum(N2)
		call iglsum(N3)
		call iglsum(N4)
		call iglsum(N5)
		call iglsum(N6)
		call glosum(MassIn(m))
	end do

	!---------------EMERGENCY SAVING
	call wait
	! call SaveVTK_VTKFortran_lib(this,NC,'SAVEcccccc')
	call SaveVTK_ascii_base64(this,NC,'SAVEtttttt')
	call SavParView(this,NC,'SAVEcccccc')

	call Save_Mesh_CGNS_Lib(this) ! -> grid.cgns
	!call Save_Mesh_CGNS_Lib(this) ! -> grid.cgns
	!call Add_Fields_To_Mesh_CGNS_Lib(this) ! -> grid.cgns

	!==========================!
	!     Initializes Time     !
	!==========================!
	Time = 0.0

	!>>> This is very handy test to perform
	if(this  < 2) then
		write(*,*) '# MassIn=', MassIn
		write(*,*) '# Number of faces on the wall        : ',N1
		write(*,*) '# Number of inflow faces             : ',N2
		write(*,*) '# Number of outflow faces            : ',N3
		write(*,*) '# Number of symetry faces            : ',N4
		write(*,*) '# Number of faces on the heated wall : ',N5
		write(*,*) '# Number of convective outflow faces : ',N6

		write(*,*) '# Variables initialized !'
	end if

	return

end subroutine IniVar
