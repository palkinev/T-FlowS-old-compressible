!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
!                                 !                                    !
!                                 !   Bojan Niceno                     !
!   Global variable definitions   !   Delft University of Technology   !
!         for all modules         !   Section Heat Transfer            !
!                                 !   niceno@duttwta.wt.tn.tudelft.nl  !
!                                 !                                    !
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
!..RCS/CVS ident
! $Id: all_mod.h90,v 1.8 2002/10/31 11:26:48 niceno Exp $
! $Source: /home/muhamed/.CVSROOT/T-Rex/Library/all_mod.h90,v $
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
!   Note: cell_n, parent, A_row, A_col, A_dia, side_c, side_cc, 
!         sideAij, are for all grids
!======================================================================!
module all_mod

  implicit none

  real, allocatable, public :: xc(:),yc(:),zc(:)
  real, allocatable, public :: Sx(:),Sy(:),Sz(:)
  real, allocatable, public :: volume(:)            ! cell's volume
  real, allocatable, public :: delta(:)             ! delta (max(dx,dy,dz))
  real, allocatable, public :: a1(:)                ! scotti
  real, allocatable, public :: a2(:)                ! scotti
  real, allocatable, public :: Dx(:),Dy(:),Dz(:)
  real, allocatable, public :: xsp(:),ysp(:),zsp(:) ! face coordinates
  real, allocatable, public :: WallDs(:), f(:)

  character, public :: name*80
  character, public :: inp*300
  integer, public   :: tn, ts(300), te(300)

  integer, public   :: NC, NS                    ! num. of nodes and cells
  integer, public   :: NbC
  integer, public   :: MNBS
  integer, public   :: NRL
  integer, public   :: Ncopy
  integer, public   :: Nmat                      ! number of materials
  LOGICAL, public   :: Mater(1024)               ! is the material present ?

  integer, allocatable, public :: material(:)     ! material markers
  integer, allocatable, public :: SideC(:,:)      !  c0, c1, c2

  integer, allocatable, public :: TypeBC(:)       ! type of boundary condition
  integer, allocatable, public :: bcmark(:)

  integer, allocatable, public :: CopyC(:)        !  might be shorter
  integer, allocatable, public :: CopyS(:,:)      !  similar to SideC

  integer, allocatable, public :: SideC1C2(:,:)   !  similar to SideC

  real, allocatable, public   :: Dxsp(:,:)       !  similar to SideC
  real, allocatable, public   :: Dysp(:,:)       !  similar to SideC
  real, allocatable, public   :: Dzsp(:,:)       !  similar to SideC

!  ! functions
!  public :: multiply_unknowns_for_solver
!  public :: deallocate_unknown
!!————————————————————————————————————————————————————————————————————————————————————————
!  contains
!  !————————————————————————————————————————————————————————————————————————————————————————
!  function multiply_unknowns_for_solver( struct_1, struct_2 ) result ( struct_result )
!    use allp_mod, only: Unknown
!    implicit none
!    
!    type(Unknown), intent(in) :: struct_1, struct_2
!    type(Unknown)             :: struct_result
!
!    allocate(struct_result % n     (-NbC:NC) )
!    !allocate(struct_result % o     (   1:NC) )
!    !allocate(struct_result % oo    (   1:NC) )
!    !allocate(struct_result % C     (   1:NC) )
!    !allocate(struct_result % X     (   1:NC) )
!    !allocate(struct_result % mean  (-NbC:NC) )
!    !allocate(struct_result % filt  (-NbC:NC) )
!    !allocate(struct_result % q     (-NbC:NC) )
!    !allocate(struct_result % fluc  (-NbC:NC) )
!    !allocate(struct_result % source(   1:NC) )
!
!    struct_result % n       = struct_1 % n       * struct_2 % n
!    !struct_result % o       = struct_1 % o       * struct_2 % o
!    !struct_result % oo      = struct_1 % oo      * struct_2 % oo
!    !struct_result % C       = struct_1 % C       * struct_2 % C
!    !struct_result % X       = struct_1 % X       * struct_2 % X
!    !struct_result % mean    = struct_1 % mean    * struct_2 % mean
!    !struct_result % filt    = struct_1 % filt    * struct_2 % filt
!    !struct_result % q       = struct_1 % q       * struct_2 % q
!    !struct_result % fluc    = struct_1 % fluc    * struct_2 % fluc
!    !struct_result % source  = struct_1 % source  * struct_2 % source
!
!    struct_result % URF   = struct_2 % URF
!    struct_result % Stol  = struct_2 % Stol
!    struct_result % bound = struct_2 % bound
!    struct_result % init  = struct_2 % init
!    struct_result % pro   = struct_2 % pro
!    struct_result % Sigma = struct_2 % Sigma
!
!  end function multiply_unknowns_for_solver
!  !————————————————————————————————————————————————————————————————————————————————————————
!  subroutine deallocate_unknown( struct )
!    use allp_mod, only: Unknown
!    implicit none
!    
!    type(Unknown) :: struct
!
!    if ( associated( struct % n      ) ) deallocate( struct % n      );
!    if ( associated( struct % o      ) ) deallocate( struct % o      );
!    if ( associated( struct % oo     ) ) deallocate( struct % oo     );
!    if ( associated( struct % C      ) ) deallocate( struct % C      );
!    if ( associated( struct % X      ) ) deallocate( struct % X      );
!    if ( associated( struct % mean   ) ) deallocate( struct % mean   );
!    if ( associated( struct % filt   ) ) deallocate( struct % filt   );
!    if ( associated( struct % q      ) ) deallocate( struct % q      );
!    if ( associated( struct % fluc   ) ) deallocate( struct % fluc   );
!    if ( associated( struct % source ) ) deallocate( struct % source );
!
!  end subroutine deallocate_unknown

end module