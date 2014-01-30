module nested_types
use prec
implicit none

type d_ub
   integer(kind=int32) :: n, lbw, ubw, lbwmax, ubwmax
   real(kind=dp), dimension(:,:), allocatable :: b
   integer(kind=int32), dimension(:), allocatable :: numrotsu
   integer(kind=int32), dimension(:,:), allocatable :: jsu
   real(kind=dp), dimension(:,:), allocatable :: csu, ssu
end type d_ub

type c_ub
   integer(kind=int32) :: n, lbw, ubw, lbwmax, ubwmax
   complex(kind=dp), dimension(:,:), allocatable :: b
   integer(kind=int32), dimension(:), allocatable :: numrotsu
   integer(kind=int32), dimension(:,:), allocatable :: jsu
   complex(kind=dp), dimension(:,:), allocatable :: csu, ssu
end type c_ub

type d_bv
   integer(kind=int32) :: n, lbw, ubw, lbwmax, ubwmax
   real(kind=dp), dimension(:,:), allocatable :: b
   integer(kind=int32), dimension(:), allocatable :: numrotsv
   integer(kind=int32), dimension(:,:), allocatable :: ksv
   real(kind=dp), dimension(:,:), allocatable :: csv, ssv
end type d_bv

type c_bv
   integer(kind=int32) :: n, lbw, ubw, lbwmax, ubwmax
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
end subroutine c_deallocate_bv

end module nested_types
