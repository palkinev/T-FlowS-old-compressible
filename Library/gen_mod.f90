!+++++++++++++++++++++++++++++++++++++!
!                                     !
!     Global variable definitions     !
!       for the mesh generator        !
!                                     !
!+++++++++++++++++++++++++++++++++++++!
!..RCS/CVS ident
! $Id: gen_mod.h90,v 1.8 2000/03/22 21:10:47 bojan Exp $
! $Source: /home/muhamed/.CVSROOT/T-Rex/Library/gen_mod.h90,v $ 
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
MODULE gen_mod

  USE allp_mod

  IMPLICIT NONE

  !INTEGER,PARAMETER :: MAXP   = 200

  REAL,ALLOCATABLE,PUBLIC    :: x(:),  y(:),  z(:)   ! node coordinates
  REAL,ALLOCATABLE,PUBLIC    :: walln(:)             ! node distance from the wall
  INTEGER,ALLOCATABLE,PUBLIC :: SideN(:,:)           ! numb, n1, n2, n3, n4
  INTEGER,ALLOCATABLE,PUBLIC :: SideCc(:,:)
						
  INTEGER,ALLOCATABLE,PUBLIC :: CellC(:,:)           ! cell's neighbours
  INTEGER,ALLOCATABLE,PUBLIC :: CellN(:,:)           ! cell nodes

  INTEGER,ALLOCATABLE,PUBLIC :: TwinN(:,:)

  INTEGER,ALLOCATABLE,PUBLIC :: NewN(:)    ! new number for the nodes and cells
  INTEGER,ALLOCATABLE,PUBLIC :: NewC(:)    ! new number for cells
  INTEGER,ALLOCATABLE,PUBLIC :: NewS(:)    ! new number for sides
  INTEGER,ALLOCATABLE,PUBLIC :: CelMar(:)  ! cell marker

  INTEGER,ALLOCATABLE,PUBLIC :: NodeN2(:,:)
  INTEGER,ALLOCATABLE,PUBLIC :: NodeN4(:,:)
  INTEGER,ALLOCATABLE,PUBLIC :: NodeN8(:,:)

  INTEGER,ALLOCATABLE,PUBLIC :: level(:)   ! refinement level

  INTEGER,PUBLIC :: MAXN, MAXB, MAXS

  INTEGER,PUBLIC :: NR(MAXP)               ! refin. levels, refin. regions

  REAL,PUBLIC    :: xp(MAXP), yp(MAXP), zp(MAXP)       ! point coordinates
  REAL,PUBLIC    :: xl(MAXP,MAXL),yl(MAXP,MAXL),zl(MAXP,MAXL),LinWgt(MAXP)
  REAL,PUBLIC    :: BlkWgt(MAXL,3), BlFaWt(MAXL,3)     ! leave this
  REAL,PUBLIC    :: FRegio(MAXP,MAXP,0:6)              ! levels, regions

  REAL,PUBLIC    :: SRegio(MAXP,0:6), Srelax(MAXP)  ! levels, regions
  LOGICAL,PUBLIC :: SdirX(MAXP), SdirY(MAXP), SdirZ(MAXP)
  INTEGER,PUBLIC :: Siter(MAXP)

  INTEGER,PUBLIC :: BlkPnt(MAXP,0:8),  & ! 0 for orientation
	     BlkRes(MAXP,6),    & ! NI,NJ,NK,NI*NJ*NK,NNo,NVo       
	     BlkFac(MAXP,6,4),  &                                    
	     BlFaLa(MAXP),      &                                   
	     Bound(MAXP,8),     &                                  
	     Period(MAXP,8),    &                                 
	     Copy(MAXP,0:8)

  INTEGER,PUBLIC :: LinPnt(MAXL,2), LinRes(MAXL)
  INTEGER,PUBLIC :: Nbloc, NP, Nline, Nsurf, Nboun, Nperi
  INTEGER,PUBLIC :: NN, NN2, NN4, NN8
  INTEGER,PUBLIC:: NSR                   ! smoothing regions
  INTEGER,PUBLIC :: NSsh                  ! number of shadow faces

  INTEGER,PUBLIC :: WallFacFst, WallFacLst

  INTEGER,PUBLIC :: ELIPSO, RECTAN, PLANE,YES,NO

  CHARACTER*4,PUBLIC :: BndFac(MAXP)

END MODULE
