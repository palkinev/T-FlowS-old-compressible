!======================================================================!
subroutine Mass_Fraction_Equation_Source(Mass_Fraction)
!------------------------------------------------------------------------------!
!   Calculates "recovered" Mass fraction field by applying inverse spatial     !
!   filter on Mass_Fraction. Then use this value to get new Source term for    !
!   Mass Fraction Equation from Candera tables.                                !
!------------------------------------------------------------------------------!
!   Only for LES method                                                        !
!----------------------------------[Modules]-----------------------------------!
  use all_mod
  use pro_mod
  use les_mod
  use rans_mod
  use par_mod
!------------------------------------------------------------------------------!
  implicit none
!-----------------------------[Parameters]-------------------------------------!
  real :: Mass_Fraction(-NbC: NC)
!-----------------------------------[Locals]-----------------------------------!
  real, allocatable :: zmix_xx(:), zmix_yy(:), zmix_zz(:)

!------------------------------[Interface]-------------------------------!
  ! these function must be pushed into par_mod, then interface is not needed
  interface 
    subroutine exchng(phi)
      use allp_mod, only : maxpro
      use all_mod
      use par_mod
      use pro_mod
      implicit none
      include 'mpif.h'
      real    :: phi(-nbc: nc)
      integer :: c1, c2, sub, rtag, stag, length, error
      integer :: status(mpi_status_size)
    end subroutine exchng

    subroutine graphi(phi, i, phii, boundary)
      use allp_mod, only: symmetry, buffer
      use all_mod
      use pro_mod
      implicit none
      integer :: i
      real    :: phi(-nbc: nc), phii(-nbc: nc)
      logical :: boundary
      integer :: s, c, c1, c2
      real    :: dphi1, dphi2, dxc1, dyc1, dzc1, dxc2, dyc2, dzc2
    end subroutine graphi
  end interface
!==============================================================================!

  ! I do not yet know why we are doing this
  call Exchng( Mass_Fraction )

  !---------------------!
  !   allocating vars   !
  !---------------------!

  allocate( zmix_xx(-NbC: NC) ); zmix_xx = 0.
  allocate( zmix_yy(-NbC: NC) ); zmix_yy = 0.
  allocate( zmix_zz(-NbC: NC) ); zmix_zz = 0.

  !---------------------------------------!
  !   calculate Mass fraction gradients   !
  !---------------------------------------!

  call GraPhi(Mass_Fraction, 1, Zmixx,.TRUE.) ! dZ/dx
  call GraPhi(Mass_Fraction, 2, Zmixy,.TRUE.) ! dZ/dy
  call GraPhi(Mass_Fraction, 3, Zmixz,.TRUE.) ! dZ/dz

  !-----------------------------------------------!
  !   calculate Mass fraction laplas components   !
  !-----------------------------------------------!
  call GraPhi(Zmixx, 1, zmix_xx,.TRUE.) ! dZ/dxx
  call GraPhi(Zmixy, 2, zmix_yy,.TRUE.) ! dZ/dyy
  call GraPhi(Zmixz, 3, zmix_zz,.TRUE.) ! dZ/dzz

  ! check b. c. ???
  ! call Get_Source_Term_From_Candera(Mass_Fraction - &
  !   volume(-NbC: NC)**(1./3.) * (zmix_xx**.2 -zmix_yy**2. -zmix_zz**2.)/24)


  !-----------------------!
  !   deallocating vars   !
  !-----------------------!
  
  if ( allocated( zmix_xx ) ) deallocate( zmix_xx )
  if ( allocated( zmix_yy ) ) deallocate( zmix_yy )
  if ( allocated( zmix_zz ) ) deallocate( zmix_zz )

end subroutine Mass_Fraction_Equation_Source
!======================================================================!

!!------------------------------------------------------------------------------!
!!-------export along the line
!  use Moin_Problem_mod
!  character(len=100) :: filename
!  integer :: c
!
!  filename = "lap-ppp.dat"
!  write(filename(5:7),'(I3.3)') this
!    
!  open(991, FILE=adjustl(trim(filename)), status='replace', access="sequential", position='append', action='write')
!  do c=1,NC
!    if ( c .ne. 0 .and. zc(c)>=-1.0 .and. zc(c)<=1.0 .and.  yc(c)<=0.03125+0.0001 .and. yc(c)>=0.03125-0.0001 ) then !problem 3
!        write(991,'(18ES26.16E3)') xc(c),& ! 1
!                                   yc(c),& ! 2
!                                   zc(c),& ! 3
!                                   U % n(c),& ! 4
!                                   Moin_Problem_u (xc(c), yc(c), 0.125),& ! 5
!                                   V % n(c),& ! 6
!                                   Moin_Problem_v (xc(c), yc(c), 0.125),& ! 7
!                                   W % n(c),& ! 8
!                                   Mass_Fraction(c),& ! 9
!                                   Moin_Problem_z (xc(c), yc(c), 0.125),& ! 10
!                                   Rho % n(c),& ! 11
!                                   Moin_Problem_rho (xc(c), yc(c), 0.125),& ! 12,
!                                   Mass_Fraction(c) - volume(c)**(2./3.) * (zmix_xx(c)**2. -zmix_yy(c)**2. -zmix_zz(c)**2.)/24.,& ! 13
!                                   Zmixx(c),& ! 14
!                                   zmix_xx(c),& ! 15
!                                   zmix_yy(c),& ! 16
!                                   zmix_zz(c),& ! 17
!                                   volume(c)**(2./3.) ! 18
!
!
!  end if
!end do
!close(991)
!call wait
!call exit(1)
!!------------------------------------------------------------------------------!