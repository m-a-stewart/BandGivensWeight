module nested_types
use prec
implicit none

type d_ub
   integer(kind=int32) :: lbw, ubw
   integer(kind=int32), private :: n, lbwmax, ubwmax
   real(kind=dp), dimension(:,:), allocatable :: b
   integer(kind=int32), dimension(:), allocatable :: numrotsu
   integer(kind=int32), dimension(:,:), allocatable :: jsu
   real(kind=dp), dimension(:,:), allocatable :: csu, ssu
end type d_ub

type c_ub
   integer(kind=int32) :: lbw, ubw
   integer(kind=int32), private :: n, lbwmax, ubwmax
   complex(kind=dp), dimension(:,:), allocatable :: b
   integer(kind=int32), dimension(:), allocatable :: numrotsu
   integer(kind=int32), dimension(:,:), allocatable :: jsu
   complex(kind=dp), dimension(:,:), allocatable :: csu, ssu
end type c_ub

type d_bv
   integer(kind=int32) :: lbw, ubw
   integer(kind=int32), private :: n, lbwmax, ubwmax
   real(kind=dp), dimension(:,:), allocatable :: b
   integer(kind=int32), dimension(:), allocatable :: numrotsv
   integer(kind=int32), dimension(:,:), allocatable :: ksv
   real(kind=dp), dimension(:,:), allocatable :: csv, ssv
end type d_bv

type c_bv
   integer(kind=int32) :: lbw, ubw
   integer(kind=int32), private :: n, lbwmax, ubwmax
   complex(kind=dp), dimension(:,:), allocatable :: b
   integer(kind=int32), dimension(:), allocatable :: numrotsv
   integer(kind=int32), dimension(:,:), allocatable :: ksv
   complex(kind=dp), dimension(:,:), allocatable :: csv, ssv
end type c_bv

interface deallocate_ub
   module procedure d_deallocate_ub, c_deallocate_ub
end interface deallocate_ub

interface deallocate_bv
   module procedure d_deallocate_bv, c_deallocate_bv
end interface deallocate_bv

interface get_n
   module procedure d_ub_get_n, d_bv_get_n, c_ub_get_n, c_bv_get_n
end interface get_n

interface get_lbwmax
   module procedure d_ub_get_lbwmax, d_bv_get_lbwmax, c_ub_get_lbwmax, c_bv_get_lbwmax
end interface get_lbwmax

interface get_ubwmax
   module procedure d_ub_get_ubwmax, d_bv_get_ubwmax, c_ub_get_ubwmax, c_bv_get_ubwmax
end interface get_ubwmax

contains

type(d_ub) function d_new_ub(n,lbwmax,ubwmax) result(ub)
  integer(kind=int32), intent(in) :: n, lbwmax, ubwmax
  ub%n=n; ub%lbwmax=lbwmax; ub%ubwmax=ubwmax
  ub%ubw=0; ub%lbw=0
  allocate(ub%b(lbwmax+ubwmax+1,n))
  allocate(ub%csu(ubwmax,n))
  allocate(ub%ssu(ubwmax,n))
  allocate(ub%jsu(ubwmax,n))
  allocate(ub%numrotsu(n))
  ub%b=0.0_dp; ub%csu=0.0_dp; ub%ssu=0.0_dp
  ub%jsu=0; ub%numrotsu=0
end function d_new_ub

subroutine d_deallocate_ub(ub)
  type(d_ub), intent(inout) :: ub
  deallocate(ub%b)
  deallocate(ub%csu)
  deallocate(ub%ssu)
  deallocate(ub%jsu)
  deallocate(ub%numrotsu)
  ub%n=0; ub%lbwmax=0; ub%ubwmax=0
end subroutine d_deallocate_ub

type(c_ub) function c_new_ub(n,lbwmax,ubwmax) result(ub)
  integer(kind=int32), intent(in) :: n, lbwmax, ubwmax
  ub%n=n; ub%lbwmax=lbwmax; ub%ubwmax=ubwmax
  ub%ubw=0; ub%lbw=0
  allocate(ub%b(lbwmax+ubwmax+1,n))
  allocate(ub%csu(ubwmax,n))
  allocate(ub%ssu(ubwmax,n))
  allocate(ub%jsu(ubwmax,n))
  allocate(ub%numrotsu(n))
  ub%b=(0.0_dp,0.0_dp); ub%csu=(0.0_dp,0.0_dp); ub%ssu=(0.0_dp,0.0_dp)
  ub%jsu=0; ub%numrotsu=0
end function c_new_ub

subroutine c_deallocate_ub(ub)
  type(c_ub), intent(inout) :: ub
  deallocate(ub%b)
  deallocate(ub%csu)
  deallocate(ub%ssu)
  deallocate(ub%jsu)
  deallocate(ub%numrotsu)
  ub%n=0; ub%lbwmax=0; ub%ubwmax=0
end subroutine c_deallocate_ub

! BV

type(d_bv) function d_new_bv(n,lbwmax,ubwmax) result(bv)
  integer(kind=int32), intent(in) :: n, lbwmax, ubwmax
  bv%n=n; bv%lbwmax=lbwmax; bv%ubwmax=ubwmax
  bv%ubw=0; bv%lbw=0
  allocate(bv%b(n,lbwmax+ubwmax+1))
  allocate(bv%csv(n,ubwmax))
  allocate(bv%ssv(n,ubwmax))
  allocate(bv%ksv(n,ubwmax))
  allocate(bv%numrotsv(n))
  bv%b=0.0_dp; bv%csv=0.0_dp; bv%ssv=0.0_dp
  bv%ksv=0; bv%numrotsv=0
end function d_new_bv

subroutine d_deallocate_bv(bv)
  type(d_bv), intent(inout) :: bv
  deallocate(bv%b)
  deallocate(bv%csv)
  deallocate(bv%ssv)
  deallocate(bv%ksv)
  deallocate(bv%numrotsv)
  bv%n=0; bv%lbwmax=0; bv%ubwmax=0
end subroutine d_deallocate_bv

type(c_bv) function c_new_bv(n,lbwmax,ubwmax) result(bv)
  integer(kind=int32), intent(in) :: n, lbwmax, ubwmax
  bv%n=n; bv%lbwmax=lbwmax; bv%ubwmax=ubwmax
  bv%ubw=0; bv%lbw=0
  allocate(bv%b(n,lbwmax+ubwmax+1))
  allocate(bv%csv(n,ubwmax))
  allocate(bv%ssv(n,ubwmax))
  allocate(bv%ksv(n,ubwmax))
  allocate(bv%numrotsv(n))
  bv%b=(0.0_dp, 0.0_dp); bv%csv=(0.0_dp, 0.0_dp); bv%ssv=(0.0_dp, 0.0_dp)
  bv%ksv=0; bv%numrotsv=0
end function c_new_bv

subroutine c_deallocate_bv(bv)
  type(c_bv), intent(inout) :: bv
  deallocate(bv%b)
  deallocate(bv%csv)
  deallocate(bv%ssv)
  deallocate(bv%ksv)
  deallocate(bv%numrotsv)
  bv%n=0; bv%lbwmax=0; bv%ubwmax=0
end subroutine c_deallocate_bv

integer(kind=int32) function d_ub_get_n(ub) result(n)
  type(d_ub) :: ub
  n=ub%n
end function d_ub_get_n

integer(kind=int32) function d_bv_get_n(bv) result(n)
  type(d_bv) :: bv
  n=bv%n
end function d_bv_get_n

integer(kind=int32) function c_ub_get_n(ub) result(n)
  type(c_ub) :: ub
  n=ub%n
end function c_ub_get_n

integer(kind=int32) function c_bv_get_n(bv) result(n)
  type(c_bv) :: bv
  n=bv%n
end function c_bv_get_n

! lbwmax

integer(kind=int32) function d_ub_get_lbwmax(ub) result(n)
  type(d_ub) :: ub
  n=ub%lbwmax
end function d_ub_get_lbwmax

integer(kind=int32) function d_bv_get_lbwmax(bv) result(n)
  type(d_bv) :: bv
  n=bv%lbwmax
end function d_bv_get_lbwmax

integer(kind=int32) function c_ub_get_lbwmax(ub) result(n)
  type(c_ub) :: ub
  n=ub%lbwmax
end function c_ub_get_lbwmax

integer(kind=int32) function c_bv_get_lbwmax(bv) result(n)
  type(c_bv) :: bv
  n=bv%lbwmax
end function c_bv_get_lbwmax

! ubwmax

integer(kind=int32) function d_ub_get_ubwmax(ub) result(n)
  type(d_ub) :: ub
  n=ub%ubwmax
end function d_ub_get_ubwmax

integer(kind=int32) function d_bv_get_ubwmax(bv) result(n)
  type(d_bv) :: bv
  n=bv%ubwmax
end function d_bv_get_ubwmax

integer(kind=int32) function c_ub_get_ubwmax(ub) result(n)
  type(c_ub) :: ub
  n=ub%ubwmax
end function c_ub_get_ubwmax

integer(kind=int32) function c_bv_get_ubwmax(bv) result(n)
  type(c_bv) :: bv
  n=bv%ubwmax
end function c_bv_get_ubwmax





end module nested_types
