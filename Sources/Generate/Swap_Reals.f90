!==============================================================================!
  subroutine Swap_Reals(a, b) 
!------------------------------------------------------------------------------!
!   Swaps two double precision reals.                                          !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  real :: a, b
!-----------------------------------[Locals]-----------------------------------!
  real :: t
!==============================================================================!

  t = a
  a = b
  b = t

  end subroutine Swap_Reals