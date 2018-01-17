!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
!                              !                                       !
!                              !                                       !
!   module that can convert    !                                       !
!   ASCII data to base64 MIME  !                                       !
!   ASCII representation of    !                                       !
!   binary format and reverse  !                                       !
!   Disadvantages: slower than !                                       !
!   binary, weights 1/3 times  !                                       !
!   more then binary           !
!   Advantages: cross-platform !                                       !
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
module Base64_Mod
! requires gfortran 4.8 ("Deferred-length character strings")

  implicit none

  character(len=64), parameter, private :: &
    base64_chars  = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/"
            !(A-Z, a-z, 0-9, +, /) - 64 symbols
            ! padd symbol "=" as 65-symbol
            ! A is 0, / is 63
  character(len=1),  parameter, private :: &
    padd_character = "="
  real(kind=8),    parameter, public :: max_real_8 = huge(1._8 )
  integer(kind=1), parameter, public :: max_int_1  = huge(1_1)
  ! defined in init
  integer(kind=1), private :: bit_size_real_8 
  integer(kind=1), private :: byte_size_real_8
  integer(kind=1), public  :: byte_size_int_1 


  logical, private :: is_initialized = .false.
  ! personal notes:

  ! _8 or kind=8 is equivalent of Float64 (16 digits after point)
  ! _1 if binary 8-bit
  ! real   (kind=8) require 64 bits or 8 bytes
  ! integer(kind=1) require 8  bits or 1 byte

  ! |7 6 5 4 3 2 1 0|7 6 5 4 3 2 1 0|7 6 5 4 3 2 1 0| <- little-endian bit order
  ! |0 1 2 3 4 5 6 7|0 1 2 3 4 5 6 7|0 1 2 3 4 5 6 7| <-    big-endian bit order

  ! mvbits is intrinsic function to COPY bits
  ! mvbits (FROM(in), FROMPOS(in), LEN(in), TO(in-out), TOPOS(in))
  ! all arguments have type integer(kind=1)
  ! example for little-endian bit order:
  ! integer(kind=1) :: from = 13        ! 00001101
  ! integer(kind=1) :: to = 6           ! 00000110
  ! call mvbits(from, 2, 2, to, 0) ! returns to = 00000111

  !character(len=:), allocatable, intent(out) :: code    ! valid since Fortran 2003

  ! functions
  public :: B64_Mod_Init
  public :: Encode_B64
  public :: Decode_B64
  public :: Get_Bit_Size_Real_8
  public :: Encode_B64_R8_a
!————————————————————————————————————————————————————————————————————————————————————————
contains
!————————————————————————————————————————————————————————————————————————————————————————

  !————————————————————————————————————————————————————————————————————————————————————————
  !---- Init module
  !————————————————————————————————————————————————————————————————————————————————————————
  subroutine B64_Mod_Init()
    implicit none

    bit_size_real_8  = Get_Bit_Size_Real_8(max_real_8)
    byte_size_real_8 = bit_size_real_8  / 8_1

    byte_size_int_1 = bit_size(max_int_1)/8_1 ! Number of bytes of kind=I1P integer

    is_initialized = .true.
  end subroutine B64_Mod_Init


  !————————————————————————————————————————————————————————————————————————————————————————
  !---- Encode bits stream to base 64 representation (ASCII) string
  !————————————————————————————————————————————————————————————————————————————————————————
  subroutine Encode_B64 (bits, padd, code)
    !|7 6 5 4 3 2 1 0|7 6 5 4 3 2 1 0|7 6 5 4 3 2 1 0|<- little-endian bit order
    !|0 1 0 0 0 0 0 1|0 1 0 0 0 0 1 0|0 1 0 0 0 0 1 1|<- 8-bits representation
    !|---8 bits(1)---|---8 bits(2)---|---8 bits(3)---|<- 3-bytes (octets)
    !|-----------|---|-------|-------|---|-----------|
    !|0 1 0 0 0 0|0 1 0 1 0 0|0 0 1 0 0 1|0 0 0 0 1 1|<- 6-bits representation
    !|5 4 3 2 1 0|5 4 3 2 1 0|5 4 3 2 1 0|5 4 3 2 1 0|<- little-endian bit order
    !|-6 bits(1)-|-6 bits(2)-|-6 bits(3)-|-6 bits(4)-|
    !|      16   |      20   |      9    |      3    |<- decimal representation
    !|-b64 index-|-b64 index-|-b64 index-|-b64 index-|
    !|      Q    |      U    |      J    |      D    |<- 4-streak base64 code
    !|-b64 char -|-b64 char -|-b64 char -|-b64 char -|
    implicit none
    integer(kind=1), intent(in)  :: bits(1:)            ! 3x8-bits streak to be encoded
    integer(kind=8)              :: bits_len            ! Length of bits array
    integer(kind=4), intent(in)  :: padd                ! Number of padding characters ('=')
    character(len=*),intent(out) :: code                ! 4-streak characters b64_code
    integer(kind=1)              :: six_bits_group(1:4) ! 4x6-bits slices
    integer(kind=8)              :: ib6                 ! 6-bit index
    integer(kind=8)              :: ib8                 ! 8-bit index

    bits_len = size( bits, dim = 1, kind = selected_int_kind(18) )

    ib6 = 1_8
    do ib8 = 1_8, bits_len, 3_8
    ! loop over array elements: 3 bytes (24 bits) scanning
      six_bits_group = 0_1
         call mvbits( bits(ib8  ), 2, 6, six_bits_group(1), 0 ) ! copy (1) 2:7 to 0:5 (1)
         call mvbits( bits(ib8  ), 0, 2, six_bits_group(2), 4 ) ! copy (1) 0:1 to 4:5 (2)
      if (ib8+1<=bits_len) then
         call mvbits( bits(ib8+1), 4, 4, six_bits_group(2), 0 ) ! copy (2) 4:7 to 0:3 (2)
         call mvbits( bits(ib8+1), 0, 4, six_bits_group(3), 2 ) ! copy (2) 0:3 to 2:5 (3)
      endif
      if (ib8+2<=bits_len) then
         call mvbits( bits(ib8+2), 6, 2, six_bits_group(3), 0 ) ! copy (3) 6:7 to 0:1 (3)
         call mvbits( bits(ib8+2), 0, 6, six_bits_group(4), 0 ) ! copy (3) 0:5 to 0:5 (4)
      endif
      six_bits_group = six_bits_group + 1_1
      code(ib6  :ib6  ) = base64_chars( six_bits_group(1): six_bits_group(1) )
      code(ib6+1:ib6+1) = base64_chars( six_bits_group(2): six_bits_group(2) )
      code(ib6+2:ib6+2) = base64_chars( six_bits_group(3): six_bits_group(3) )
      code(ib6+3:ib6+3) = base64_chars( six_bits_group(4): six_bits_group(4) )
      ib6 = ib6 + 4_8
    enddo
    if (padd>0) code(len(code)-padd+1:)=repeat(padd_character,padd)
  end subroutine Encode_B64
  !————————————————————————————————————————————————————————————————————————————————————————
  !---- Decode base64 representation(ASCII) strings of 4 characters into bits stream
  !————————————————————————————————————————————————————————————————————————————————————————
  subroutine Decode_B64 ( b64_code, bits )
    !|-b64 char -|-b64 char -|-b64 char -|-b64 char -|
    !|      Q    |      U    |      J    |      D    |<- 4-streak base64 code 
    !|-b64 index-|-b64 index-|-b64 index-|-b64 index-|
    !|      16   |      20   |      9    |      3    |<- decimal representation
    !|-6 bits(1)-|-6 bits(2)-|-6 bits(3)-|-6 bits(4)-|
    !|0 1 0 0 0 0|0 1 0 1 0 0|0 0 1 0 0 1|0 0 0 0 1 1|<- 6-bits representation
    !|5 4 3 2 1 0|5 4 3 2 1 0|5 4 3 2 1 0|5 4 3 2 1 0|<- little-endian bit order
    !|-----------|---|-------|-------|---|-----------|
    !|---8 bits(1)---|---8 bits(2)---|---8 bits(3)---|<- 3-bytes (octets)
    !|0 1 0 0 0 0 0 1|0 1 0 0 0 0 1 0|0 1 0 0 0 0 1 1|<- 8-bits representation
    !|7 6 5 4 3 2 1 0|7 6 5 4 3 2 1 0|7 6 5 4 3 2 1 0|<- little-endian bit order
    implicit none

    character(len=*),intent(in)  :: b64_code            ! 4-streak characters b64_code
    integer(kind=1), intent(out) :: bits(1:)            ! (1-3)x8-bits decoded streak
    integer(kind=8)              :: bits_len            ! Length of bits array
    integer(kind=1)              :: six_bits_group(1:4) ! 4x6-bits slices
    integer(kind=8)              :: ib6                 ! 6-bit index
    integer(kind=8)              :: ib8                 ! 8-bit index

    bits_len = size( bits, dim = 1, kind = selected_int_kind(18) )

    ib8 = 1_8
    ! loop over b64_code characters: 3 bytes (24 bits) scanning
    do ib6 = 1_8, len(b64_code), 4_8
      six_bits_group = 0_1
      ! b64 -> binary 6-bit representation
      six_bits_group(1) = index( base64_chars, b64_code(ib6  : ib6  ) ) - 1
      six_bits_group(2) = index( base64_chars, b64_code(ib6+1: ib6+1) ) - 1
      six_bits_group(3) = index( base64_chars, b64_code(ib6+2: ib6+2) ) - 1
      six_bits_group(4) = index( base64_chars, b64_code(ib6+3: ib6+3) ) - 1
        ! filling 8-bits array from left to right
        call mvbits( six_bits_group(1), 0, 6, bits(ib8  ), 2 ) ! copy (1) 0:5 to 2:7 (1)
        call mvbits( six_bits_group(2), 4, 2, bits(ib8  ), 0 ) ! copy (2) 4:5 to 0:1 (1)
      if (ib8+1<=bits_len) then
        call mvbits( six_bits_group(2), 0, 4, bits(ib8+1), 4 ) ! copy (2) 0:3 to 4:7 (2)
        call mvbits( six_bits_group(3), 2, 4, bits(ib8+1), 0 ) ! copy (3) 2:5 to 0:3 (2)
      endif
      if (ib8+2<=bits_len) then
        call mvbits( six_bits_group(3), 0, 2, bits(ib8+2), 6 ) ! copy (3) 0:1 to 6:7 (3)
        call mvbits( six_bits_group(4), 0, 6, bits(ib8+2), 0 ) ! copy (4) 0:5 to 0:5 (3)
      endif
      ib8 = ib8 + 3_8
    enddo
  end subroutine Decode_B64
  !————————————————————————————————————————————————————————————————————————————————————————
  !---- Encode real(kind=8) array to base64 and append padd character where needed
  !————————————————————————————————————————————————————————————————————————————————————————
  subroutine Encode_B64_R8_a(n, code)
    !  character(len=:), allocatable :: code64
    !  call Encode_B64_R8_a(n=[1._8,2._8], code=code64)
    !  print "(A)", code64
    !=> AAAAAAAA8D8AAAAAAAAAQA== <<<

    ! personal note: 
    ! 3-byte -> 4 characters
    ! but padding character introduces exceptions:
    ! '='  means that last 4-characters group contain only 2 byte
    ! '==' means that last 4-characters group contain only 1 byte

    ! see example from subroutine Encode_B64 without last byte
    ! imagine this is last 3-byte group we have in array:

    !|7 6 5 4 3 2 1 0|7 6 5 4 3 2 1 0|7 6 5 4 3 2 1 0|<- little-endian bit order
    !|0 1 0 0 0 0 0 1|0 1 0 0 0 0 1 0|- - - - - - - -|<- 8-bits representation
    !|---8 bits(1)---|---8 bits(2)---|---8 bits(3)---|<- 3-bytes (octets)
    !|-----------|---|-------|-------|---|-----------|
    !|0 1 0 0 0 0|0 1 0 1 0 0|0 0 1 0 0 0|- - - - - -|<- 6-bits representation
    !|5 4 3 2 1 0|5 4 3 2 1 0|5 4 3 2 1 0|5 4 3 2 1 0|<- little-endian bit order
    !|-6 bits(1)-|-6 bits(2)-|-6 bits(3)-|-6 bits(4)-|
    !|      16   |      20   |      8    |      -    |<- decimal representation
    !|-b64 index-|-b64 index-|-b64 index-|-b64 index-|
    !|      Q    |      U    |      I    |      =    |<- 4-streak base64 code
    !|-b64 char -|-b64 char -|-b64 char -|-b64 char -|


    implicit none
    real(kind=8),       intent(in)             :: n(1:)   ! Array of real*8 to be encoded
    character(len=:), allocatable, intent(out) :: code    ! Encoded array
    integer(kind=1), allocatable               :: n1(:)   ! One byte integer array containing n
    integer(kind=4)                            :: padd    ! Number of padding characters ('=')
    integer(kind=8)                            :: ns      ! Size of n

    ! byte_size_real_8 = 8
    ! ((ns*byte_size_real_8+2)/3)*3) is number of bytes
    ! ((ns*byte_size_real_8+2)/3)*4) is number of base64 characters
    ns = size(n, dim=1)
    allocate(n1(1:((ns*byte_size_real_8+2)/3)*3))
    n1 = 0_1

    if (.not. allocated(code)) code = repeat(' ',((ns*byte_size_real_8+2)/3)*4)
    n1 = transfer(n, n1)
    padd = mod((ns*byte_size_real_8),3_8)
    if (padd>0_4) padd = 3_4 - padd

    call Encode_B64(bits=n1, padd=padd, code=code)
  end subroutine Encode_B64_R8_a
  !————————————————————————————————————————————————————————————————————————————————————————
  !---- Encode integer(kind=1) array to base64.
  !————————————————————————————————————————————————————————————————————————————————————————
  subroutine b64_encode_I1_a(n, code)
    ! Encode array numbers to base64 (I1P)
    !
    ! use befor64
    ! use penf
    ! character(len=:), allocatable :: code64
    ! call b64_encode(n=[120_I1P,-1_I1P], code=code64)
    ! print "(A)", code64
    !=> eP8= <<<
    integer(kind=1),               intent(in)  :: n(1:) ! Array of numbers to be encoded
    character(len=:), allocatable, intent(out) :: code  ! Encoded array
    integer(kind=1),     allocatable           :: n1(:) ! One byte integer array containing n
    integer(kind=4)                            :: padd  ! Number of padding characters ('=')
    integer(kind=8)                            :: ns    ! Size of n

    ns = size(n,dim=1)
    allocate(n1(1:((ns*byte_size_int_1+2)/3)*3))
    n1 = 0_1
    if (.not. allocated(code)) code = repeat(' ',((ns*byte_size_int_1+2)/3)*4)
    n1 = transfer(n,n1)
    padd = mod((ns*byte_size_int_1),3_8)
    if (padd>0_4) padd = 3_4 - padd
    call Encode_B64(bits=n1,padd=padd,code=code)
  end subroutine b64_encode_I1_a

  !————————————————————————————————————————————————————————————————————————————————————————
  !---- Encode (Base64) a dataarray with 3 components of rank 1 (kind=8).
  !————————————————————————————————————————————————————————————————————————————————————————
  subroutine Encode_B64_R8_a_xyz(x, y, z, code)
    real(kind=8),    intent(in)      :: x(1:)  ! X component
    real(kind=8),    intent(in)      :: y(1:)  ! Y component
    real(kind=8),    intent(in)      :: z(1:)  ! Z component
    character(len=:), allocatable    :: code   ! Encoded base64 dataarray
    integer(kind=1),  allocatable    :: xyz(:) ! Packed data
    integer(kind=4)                  :: nn     ! Number of elements
    integer(kind=4)                  :: n      ! Counter

    nn = size(x, dim=1)
    call pack_data_I4_R8( a1=[int(3*nn*byte_size_real_8, kind=4)], &
                          a2=[(x(n), y(n), z(n), n=1, nn)],        &
                          packed=xyz )
    call b64_encode_I1_a(n=xyz, code=code)
  end subroutine Encode_B64_R8_a_xyz

  subroutine pack_data_I4_R8(a1, a2, packed)
  !————————————————————————————————————————————————————————————————————————————————————————
  !---- Pack different kinds of data into single kind=1 array
  !————————————————————————————————————————————————————————————————————————————————————————
    !
    ! integer(kind=4)              :: a1(1)
    ! real(kind=8)                 :: a2(1)
    ! integer(kind=1), allocatable :: pack(:)
    ! a1(1) = 0
    ! a2(1) = 1
    ! call pack_data(a1=a1, a2=a2, packed=pack)
    ! print *, pack(size(pack, dim=1))
    !=> 63 <<<
    integer(kind=4),              intent(in)    :: a1(1:)    ! First data stream
    real(kind=8),                 intent(in)    :: a2(1:)    ! Second data stream
    integer(kind=1), allocatable, intent(inout) :: packed(:) ! Packed data into kind=1 array
    integer(kind=1), allocatable                :: p1(:)     ! Temporary packed data of first stream
    integer(kind=1), allocatable                :: p2(:)     ! Temporary packed data of second stream

    p1 = transfer(a1,p1)
    p2 = transfer(a2,p2)
    packed = [p1,p2]
  end subroutine pack_data_I4_R8

  !————————————————————————————————————————————————————————————————————————————————————————
  !---- Compute the number of bits of a real*8 variable
  !————————————————————————————————————————————————————————————————————————————————————————
  function Get_Bit_Size_Real_8(i) result(bits)
    ! print Get_Bit_Size_Real_8(1._8)
    ! => 64 <<<
    real(kind=8), intent(in) :: i       ! Real variable whose number of bits must be computed
    integer(kind=1)          :: bits    ! Number of bits of r
    integer(kind=1)          :: mold(1) ! "Molding" dummy variable for bits counting

    bits = size(transfer(i, mold), dim=1, kind=1) * 8_1
    !bits = 64
  end function Get_Bit_Size_Real_8

end module Base64_Mod
