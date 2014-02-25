module sweeps
use prec
use rotation
implicit none

!
! These types represent a linear transformation
!
! Q = Q_1 Q_2 ... Q_{numsweeps}
!
! where
!
! Q_k = G_{k,1} ... G_{k,n-1}


type d_sweeps
   integer(kind=int32), private :: maxsweeps, n
   logical, private :: transposed
   integer(kind=int32) :: numsweeps
   real(kind=dp), allocatable, dimension(:,:) :: cs, ss
end type d_sweeps

type c_sweeps
   integer(kind=int32), private :: maxsweeps, n
   logical, private :: transposed
   integer(kind=int32) :: numsweeps
   complex(kind=dp), allocatable, dimension(:,:) :: cs, ss
end type c_sweeps

interface deallocate_sweeps
   module procedure d_deallocate_sweeps, c_deallocate_sweeps
end interface

interface sweeps_times_general
   module procedure d_sweeps_times_general, c_sweeps_times_general
end interface

contains

type(d_sweeps) function d_new_sweeps(n,maxsweeps) result(sw)
  integer(kind=int32), intent(in) :: n, maxsweeps
  sw%n=n; sw%maxsweeps=maxsweeps
  allocate(sw%cs(n,maxsweeps), sw%ss(n,maxsweeps))
  sw%transposed=.false.
  sw%numsweeps=0
  sw%cs=1.0_dp; sw%ss=0.0_dp
end function d_new_sweeps

type(c_sweeps) function c_new_sweeps(n,maxsweeps) result(sw)
  integer(kind=int32), intent(in) :: n, maxsweeps
  sw%n=n; sw%maxsweeps=maxsweeps
  allocate(sw%cs(n,maxsweeps), sw%ss(n,maxsweeps))
  sw%transposed=.false.
  sw%numsweeps=0
  sw%cs=(1.0_dp,0.0_dp); sw%ss=(0.0_dp,0.0_dp)
end function c_new_sweeps

subroutine d_deallocate_sweeps(sw)
  type(d_sweeps), intent(inout) :: sw
  deallocate(sw%cs, sw%ss)
end subroutine d_deallocate_sweeps

subroutine c_deallocate_sweeps(sw)
  type(c_sweeps), intent(inout) :: sw
  deallocate(sw%cs, sw%ss)
end subroutine c_deallocate_sweeps

subroutine d_sweeps_times_general(sw,a)
  type(d_sweeps) :: sw
  real(kind=dp), dimension(:,:), intent(inout) :: a
  type(d_rotation) :: rot
  integer(kind=int32) :: j, n, l
  if (sw%transposed) then
     do l=1,sw%numsweeps
        n=size(a,1)
        do j=1,n-1
           rot%cosine=sw%cs(j,l); rot%sine=sw%ss(j,l)
           call d_rotation_times_general(trp_rot(rot),a,j,j+1)
        end do
     end do
  else
     do l=sw%numsweeps,1,-1
        n=size(a,1)
        do j=n-1,1,-1
           rot%cosine=sw%cs(j,l); rot%sine=sw%ss(j,l)
           call d_rotation_times_general(rot,a,j,j+1)
        end do
     end do
  end if
end subroutine d_sweeps_times_general

subroutine c_sweeps_times_general(sw,a)
  type(c_sweeps) :: sw
  complex(kind=dp), dimension(:,:), intent(inout) :: a
  type(c_rotation) :: rot
  integer(kind=int32) :: j, n, l
  if (sw%transposed) then
     do l=1,sw%numsweeps
        n=size(a,1)
        do j=1,n-1
           rot%cosine=sw%cs(j,l); rot%sine=sw%ss(j,l)
           call c_rotation_times_general(trp_rot(rot),a,j,j+1)
        end do
     end do
  else
     do l=sw%numsweeps,1,-1
        n=size(a,1)
        do j=n-1,1,-1
           rot%cosine=sw%cs(j,l); rot%sine=sw%ss(j,l)
           call c_rotation_times_general(rot,a,j,j+1)
        end do
     end do
  end if
end subroutine c_sweeps_times_general

end module sweeps
