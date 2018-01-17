!======================================================================!
  SUBROUTINE UserPerturb2(n, Dom)
!----------------------------------------------------------------------!
!   Perturbs the flow field for any flow.                              !
!----------------------------------------------------------------------!
!------------------------------[Modules]-------------------------------!
  USE all_mod
  USE pro_mod
  USE par_mod
  USE les_mod
!----------------------------------------------------------------------!
  IMPLICIT NONE
!------------------------------[Calling]-------------------------------!
  INTEGER :: n, Dom
!-------------------------------[Locals]-------------------------------!
  INTEGER :: c, q
  REAL    :: randn, RRR
  integer, allocatable :: a(:) 
!--------------------------------[CVS]---------------------------------!
!  $Id: UserPerturb2.f90,v 1.11 2002/10/30 16:30:03 niceno Exp $  
!  $Source: /home/muhamed/.CVSROOT/T-Rex/User/UserPerturb2.f90,v $  
!======================================================================!

  call random_seed(size=q)
  allocate(a(q)) 
  a = This*100 + n
  call random_seed(PUT = a)    ! Set user seed

!----------------------------------!
!      add fluctuating values      !
!----------------------------------!
  do c=-NbC,NC
    if(material(c) == Dom) then
    call random_number(randn)

!---- 10 % of the maximum values
!    U % n(c)  = U % n(c) + .1*fac*(0.5-randn) & 
!              * max(abs(U % n(c)), abs(V % n(c)), abs(W % n(c))) 
!    U % o(c)  = U % n(c)
!    U % oo(c) = U % n(c)

!---- 10 % of the maximum values
!    V % n(c)  = V % n(c)   + .1*fac*(0.5-randn) &
!              * max(abs(U % n(c)), abs(V % n(c)), abs(W % n(c))) 
!    V % o(c)  = V % n(c)
!    V % oo(c) = V % n(c) 

!---- 10 % of the maximum values
!    W % n(c)  = W % n(c)   + .1*fac*(0.5-randn) &
!              * max(abs(U % n(c)), abs(V % n(c)), abs(W % n(c))) 
!    W % o(c)  = W % n(c)
!    W % oo(c) = W % n(c)


!     if(sqrt((xc(c)-0.0)**2.0+(yc(c)-0.0)**2.0+(zc(c)-1.1)**2.0)<0.1) then
!     Zmixr % n(c)  = 1.0
!     Zmixr % o(c)  = Zmixr % n(c)
!     Zmixr % oo(c) = Zmixr % n(c)
!
!     Zmix % n(c)  = Zmixr % n(c) / Rho % n(c)
!     endif

!     if(zc(c) > 0.0) then
!     sigma = 5.0
!     Zmix % n(c)   = 1.0
!     Rho % n(c)    = 1.0 / (1.0 + sigma) 
!     Rho % o(c)    = Rho % n(c)
!     Rho % oo(c)   = Rho % n(c)
!
!     Ur  % n(c)    = Ur  % n(c) / (1.0 + sigma) 
!     Ur  % o(c)    = Ur  % o(c) / (1.0 + sigma) 
!     Ur  % oo(c)   = Ur  % oo(c) / (1.0 + sigma) 
!     Vr  % n(c)    = Vr  % n(c) / (1.0 + sigma) 
!     Vr  % o(c)    = Vr  % o(c) / (1.0 + sigma) 
!     Vr  % oo(c)   = Vr  % oo(c) / (1.0 + sigma) 
!     Wr  % n(c)    = Wr  % n(c) / (1.0 + sigma) 
!     Wr  % o(c)    = Wr  % o(c) / (1.0 + sigma) 
!     Wr  % oo(c)   = Wr  % oo(c) / (1.0 + sigma) 

!     Ur  % C(c)    = Ur  % C(c) / (1.0 + sigma)
!     Ur  % Co(c)   = Ur  % Co(c) / (1.0 + sigma)
!     Vr  % C(c)    = Vr  % C(c) / (1.0 + sigma)
!     Vr  % Co(c)   = Vr  % Co(c) / (1.0 + sigma)
!     Wr  % C(c)    = Wr  % C(c) / (1.0 + sigma)
!     Wr  % Co(c)   = Wr  % Co(c) / (1.0 + sigma)

!     Ur  % Do(c)   = Ur  % Do(c) / (1.0 + sigma)
!     Vr  % Do(c)   = Vr  % Do(c) / (1.0 + sigma)
!     Wr  % Do(c)   = Wr  % Do(c) / (1.0 + sigma)

!     Ur  % X(c)    = Ur  % X(c) / (1.0 + sigma)
!     Vr  % X(c)    = Vr  % X(c) / (1.0 + sigma)
!     Wr  % X(c)    = Wr  % X(c) / (1.0 + sigma)
!     Ur  % Xo(c)   = Ur  % Xo(c) / (1.0 + sigma)
!     Vr  % Xo(c)   = Vr  % Xo(c) / (1.0 + sigma)
!     Wr  % Xo(c)   = Wr  % Xo(c) / (1.0 + sigma)

!     P % n(c)      = P  % n(c) / (1.0 + sigma)
!     PP % n(c)     = PP  % n(c) / (1.0 + sigma)

!     Px(c)         = Px(c) / (1.0 + sigma)
!     Py(c)         = Py(c) / (1.0 + sigma)
!     Pz(c)         = Pz(c) / (1.0 + sigma)


!     if(this < 2) then
!     write(*,*) Zmix % n(1), Zmixr % n(1)
!     end if
!     Zmixr % n(c)  = Zmix % n(c) * Rho % n(c) 
!     Zmixr % o(c)  = Zmixr % n(c)
!     Zmixr % oo(c) = Zmixr % n(c)
!     if(this < 2) then
!     write(*,*) Zmix % n(1), Zmixr % n(1)
!     end if
!    endif

    RRR = sqrt(xc(c)**2.0 + yc(c)**2.0)

    if(zc(c) < -0.0 .and. RRR > 0.73 .and. RRR < 1.0) then
     if (zc(c) < -0.1) then
     !Wr  % n(c)  = 1.0
     !Wr  % o(c)  = 1.0
     !Wr  % oo(c) = 1.0
     W   % n(c)  = 1.0
     endif
     Zmix%n(c)   = 0.0
     !Zmixr%n(c)  = 0.0
     !Zmixr%o(c)  = 0.0
     !Zmixr%oo(c) = 0.0
     Rho % n(c)  = 1.0
    endif

    if(zc(c) < -0.0 .and. RRR > 0.365 .and. RRR < 0.678) then
     if (zc(c) < -0.1) then
     !Wr  % n(c)  = 0.444385
     !Wr  % o(c)  = 0.444385
     !Wr  % oo(c) = 0.444385
     W   % n(c)  = 0.444385
     endif
     Zmix%n(c)   = 0.0
     !Zmixr%n(c)  = 0.0
     !Zmixr%o(c)  = 0.0
     !Zmixr%oo(c) = 0.0
     Rho % n(c)  = 1.0
    endif

    end if
  end do

!   do s=1,NS
!    c1=SideC(1,s)
!    c2=SideC(2,s)
!
!    RHOs      = fF(s)*Rho%n(c1)+(1.0-fF(s))*Rho%n(c2)
!    Flux(s)   = Flux_u(s) * RHOs
!   enddo

  END SUBROUTINE UserPerturb2
