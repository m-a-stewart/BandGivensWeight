module band
use prec
use rotation
use utility
implicit none

interface get_el_bc
   module procedure d_get_el_bc, c_get_el_bc
end interface get_el_bc

interface set_el_bc
   module procedure d_set_el_bc, c_set_el_bc
end interface set_el_bc

interface get_el_br
   module procedure d_get_el_br, c_get_el_br
end interface get_el_br

interface set_el_br
   module procedure d_set_el_br, c_set_el_br
end interface set_el_br

! to general

interface bc_to_general
   module procedure d_bc_to_general, c_bc_to_general
end interface bc_to_general

interface br_to_general
   module procedure d_br_to_general, c_br_to_general
end interface br_to_general

! printing

interface print_bc
   module procedure d_print_bc, c_print_bc
end interface print_bc

interface print_br
   module procedure d_print_br, c_print_br
end interface print_br

! rotations

interface rotation_times_tbc
   module procedure d_rotation_times_tbc, c_rotation_times_tbc
end interface rotation_times_tbc

interface tbc_times_rotation
   module procedure d_tbc_times_rotation, c_tbc_times_rotation
end interface tbc_times_rotation

interface tbr_times_rotation
   module procedure d_tbr_times_rotation, c_tbr_times_rotation
end interface tbr_times_rotation

interface rotation_times_tbr
   module procedure d_rotation_times_tbr, c_rotation_times_tbr
end interface rotation_times_tbr

contains

! Get and set functions for band row and band column.

real(kind=dp) function d_get_el_bc(b,ubw,j,k) result(x)
  real(kind=dp), dimension(:,:), intent(in) :: b
  integer(kind=int32), intent(in) :: ubw,j,k
  !
  x=b(j-k+ubw+1,k)
end function d_get_el_bc

complex(kind=dp) function c_get_el_bc(b,ubw,j,k) result(x)
  complex(kind=dp), dimension(:,:), intent(in) :: b
  integer(kind=int32), intent(in) :: ubw,j,k
  !
  x=b(j-k+ubw+1,k)
end function c_get_el_bc

subroutine d_set_el_bc(b,ubw,j,k,x)
  real(kind=dp), dimension(:,:), intent(inout) :: b
  integer(kind=int32), intent(in) :: ubw,j,k
  real(kind=dp), intent(in) :: x
  !
  b(j-k+ubw+1,k) = x
end subroutine d_set_el_bc

subroutine c_set_el_bc(b,ubw,j,k,x)
  complex(kind=dp), dimension(:,:), intent(inout) :: b
  integer(kind=int32), intent(in) :: ubw,j,k
  complex(kind=dp), intent(in) :: x
  !
  b(j-k+ubw+1,k) = x
end subroutine c_set_el_bc

real(kind=dp) function d_get_el_br(b,lbw,j,k) result(x)
  real(kind=dp), dimension(:,:), intent(in) :: b
  integer(kind=int32), intent(in) :: lbw,j,k
  !
  x=b(j,k-j+lbw+1)
end function d_get_el_br

complex(kind=dp) function c_get_el_br(b,lbw,j,k) result(x)
  complex(kind=dp), dimension(:,:), intent(in) :: b
  integer(kind=int32), intent(in) :: lbw,j,k
  !
  x=b(j,k-j+lbw+1)
end function c_get_el_br

subroutine d_set_el_br(b,lbw,j,k,x)
  real(kind=dp), dimension(:,:), intent(inout) :: b
  integer(kind=int32), intent(in) :: lbw,j,k
  real(kind=dp), intent(in) :: x
  !
  b(j,k-j+lbw+1)=x
end subroutine d_set_el_br

subroutine c_set_el_br(b,lbw,j,k,x)
  complex(kind=dp), dimension(:,:), intent(inout) :: b
  integer(kind=int32), intent(in) :: lbw,j,k
  complex(kind=dp), intent(in) :: x
  !
  b(j,k-j+lbw+1)=x
end subroutine c_set_el_br

! band to general

subroutine d_bc_to_general(bc,lbw,ubw,a)
  real(kind=dp), dimension(:,:), intent(in) :: bc
  integer(kind=int32), intent(in) :: ubw, lbw
  real(kind=dp), dimension(:,:), intent(out) :: a
  ! 
  integer(kind=int32) :: d, n, k
  a=0.0_dp
  n=size(bc,2)
  do d=1,ubw+1
     do k=ubw-d+2,n
        a(d+k-ubw-1,k)=bc(d,k)
     end do
  end do
  do d=ubw+2, ubw+lbw+1
     do k=1,n-d+ubw+1
        a(d+k-ubw-1,k)=bc(d,k)
     end do
  end do
end subroutine d_bc_to_general

subroutine c_bc_to_general(bc,lbw,ubw,a)
  complex(kind=dp), dimension(:,:), intent(in) :: bc
  integer(kind=int32), intent(in) :: ubw, lbw
  complex(kind=dp), dimension(:,:), intent(out) :: a
  ! 
  integer(kind=int32) :: d, n, k
  a=(0.0_dp,0.0_dp)
  n=size(bc,2)
  do d=1,ubw+1
     do k=ubw-d+2,n
        a(d+k-ubw-1,k)=bc(d,k)
     end do
  end do
  do d=ubw+2, ubw+lbw+1
     do k=1,n-d+ubw+1
        a(d+k-ubw-1,k)=bc(d,k)
     end do
  end do
end subroutine c_bc_to_general

subroutine d_br_to_general(br,lbw,ubw,a)
  real(kind=dp), dimension(:,:), intent(in) :: br
  integer(kind=int32), intent(in) :: ubw, lbw
  real(kind=dp), dimension(:,:), intent(out) :: a
  !
  integer(kind=int32) :: d, n, j
  a=0.0_dp
  n=size(br,1)
  do d=1,lbw+1
     do j=lbw-d+2, n
        a(j,d+j-lbw-1)=br(j,d)
     end do
  end do
  do d=lbw+2,ubw+lbw+1
     do j=1, n-d+lbw+1
        a(j,d+j-lbw-1)=br(j,d)
     end do
  end do
end subroutine d_br_to_general

subroutine c_br_to_general(br,lbw,ubw,a)
  complex(kind=dp), dimension(:,:), intent(in) :: br
  integer(kind=int32), intent(in) :: ubw, lbw
  complex(kind=dp), dimension(:,:), intent(out) :: a
  !
  integer(kind=int32) :: d, n, j
  a=(0.0_dp,0.0_dp)
  n=size(br,1)
  do d=1,lbw+1
     do j=lbw-d+2, n
        a(j,d+j-lbw-1)=br(j,d)
     end do
  end do
  do d=lbw+2,ubw+lbw+1
     do j=1, n-d+lbw+1
        a(j,d+j-lbw-1)=br(j,d)
     end do
  end do
end subroutine c_br_to_general

! printing

subroutine d_print_bc(bc,lbw,ubw)
  real(kind=dp), dimension(:,:), intent(in) :: bc
  integer(kind=int32), intent(in) :: ubw, lbw
  !
  real(kind=dp), dimension(size(bc,2),size(bc,2)) :: a
  call d_bc_to_general(bc,lbw,ubw,a)
  call print_matrix(a)
end subroutine d_print_bc

subroutine c_print_bc(bc,lbw,ubw)
  complex(kind=dp), dimension(:,:), intent(in) :: bc
  integer(kind=int32), intent(in) :: ubw, lbw
  !
  complex(kind=dp), dimension(size(bc,2),size(bc,2)) :: a
  call c_bc_to_general(bc,lbw,ubw,a)
  call print_matrix(a)
end subroutine c_print_bc

subroutine d_print_br(br,lbw,ubw)
  real(kind=dp), dimension(:,:), intent(in) :: br
  integer(kind=int32), intent(in) :: ubw, lbw
  !
  real(kind=dp), dimension(size(br,1),size(br,1)) :: a
  call d_br_to_general(br,lbw,ubw,a)
  call print_matrix(a)
end subroutine d_print_br

subroutine c_print_br(br,lbw,ubw)
  complex(kind=dp), dimension(:,:), intent(in) :: br
  integer(kind=int32), intent(in) :: ubw, lbw
  !
  complex(kind=dp), dimension(size(br,1),size(br,1)) :: a
  call c_br_to_general(br,lbw,ubw,a)
  call print_matrix(a)
end subroutine c_print_br


! Rotations for band matrices with aligned columns

! Rotation times truncated band matrix with aligned columns.
! First l columns are not modified.
subroutine d_rotation_times_tbc(r,b,n,lbw,ubw,l,j)
  real(kind=dp), dimension(:,:), intent(inout) :: b
  type(d_rotation), intent(in) :: r
  integer(kind=int32), intent(in) :: j,n, ubw,lbw, l
  !
  real(kind=dp) :: c, s, tmp
  integer(kind=int32) :: k, k0, k1, d
  c=r%cosine; s=r%sine
  ! first and last column operated on
  k0=max(j-lbw+1,l+1)
  k1=max(min(j+ubw,n),l+1)
  ! loop over relevant columns in b.
  do k=k0,k1
     d=j-k+ubw+1
     tmp = b(d,k)
     b(d,k)=c*tmp - s*b(d+1,k)
     b(d+1,k) = s*tmp + c*b(d+1,k)
  end do
end subroutine d_rotation_times_tbc

! Rotation times truncated band matrix with aligned columns.
! First l columns are not modified.
subroutine c_rotation_times_tbc(r,b,n,lbw,ubw,l,j)
  complex(kind=dp), dimension(:,:), intent(inout) :: b
  type(c_rotation), intent(in) :: r
  integer(kind=int32), intent(in) :: j,n, ubw,lbw, l
  !
  complex(kind=dp) :: c, s, tmp
  integer(kind=int32) :: k, k0, k1, d
  c=r%cosine; s=r%sine
  ! first and last column operated on
  k0=max(j-lbw+1,l+1)
  k1=max(min(j+ubw,n),l+1)
  ! loop over relevant columns in b.
  do k=k0,k1
     d=j-k+ubw+1
     tmp = b(d,k)
     b(d,k)=c*tmp - conjg(s)*b(d+1,k)
     b(d+1,k) = s*tmp + c*b(d+1,k)
  end do
end subroutine c_rotation_times_tbc

! apply a rotation to columns k and k+1 of a truncated band matrix with aligned columns.
! last l rows are not modified.
subroutine d_tbc_times_rotation(b,n,lbw,ubw,l,r,k)
  real(kind=dp), dimension(:,:), intent(inout) :: b
  type(d_rotation), intent(in) :: r
  integer(kind=int32), intent(in) :: k,n, ubw,lbw, l
  !
  real(kind=dp) :: c, s, tmp
  integer(kind=int32) :: d, j0, j1, j
  c=r%cosine; s=r%sine
  j0=min(n-l,(max(1,k-ubw+1)))
  j1=min(n-l,k+lbw)
  do j=j0,j1
     d=j-k+ubw+1
     tmp=b(d,k)
     b(d,k)=c*tmp+s*b(d-1,k+1)
     b(d-1,k+1)=-s*tmp+c*b(d-1,k+1)
  end do
end subroutine d_tbc_times_rotation

subroutine c_tbc_times_rotation(b,n,lbw,ubw,l,r,k)
  complex(kind=dp), dimension(:,:), intent(inout) :: b
  type(d_rotation), intent(in) :: r
  integer(kind=int32), intent(in) :: k,n, ubw,lbw, l
  !
  complex(kind=dp) :: c, s, tmp
  integer(kind=int32) :: d, j0, j1, j
  c=r%cosine; s=r%sine
  j0=min(n-l,(max(1,k-ubw+1)))
  j1=min(n-l,k+lbw)
  do j=j0,j1
     d=j-k+ubw+1
     tmp=b(d,k)
     b(d,k)=c*tmp+s*b(d-1,k+1)
     b(d-1,k+1)=-conjg(s)*tmp+c*b(d-1,k+1)
  end do
end subroutine c_tbc_times_rotation

! rotations for band matrices with aligned rows.

! apply a rotation to columns k and k+1 of a truncated band matrix with aligned rows.
! The rotation is not applied to the last l rows.
subroutine d_tbr_times_rotation(b,m,lbw,ubw,l,r,k)
  real(kind=dp), dimension(:, :), intent(inout) :: b
  type(d_rotation), intent(in) :: r
  integer(kind=int32), intent(in) :: k,m, ubw,lbw, l
  !
  real(kind=dp) :: c, s, tmp
  integer(kind=int32) :: j, j0, j1, d
  c=r%cosine; s=r%sine
  ! first and last row operated on
  j0=min(m-l,max(1,k-ubw+1))
  j1=min(m-l,k+lbw)
  do j=j0,j1
     d=k-j+lbw+1
     tmp=b(j,d)
     b(j,d)=c*tmp+s*b(j,d+1)
     b(j,d+1)=-s*tmp+c*b(j,d+1)
  end do
end subroutine d_tbr_times_rotation

subroutine c_tbr_times_rotation(b,m,lbw,ubw,l,r,k)
  complex(kind=dp), dimension(:,:), intent(inout) :: b
  type(c_rotation), intent(in) :: r
  integer(kind=int32), intent(in) :: k,m, ubw,lbw, l
  !
  complex(kind=dp) :: c, s, tmp
  integer(kind=int32) :: j, j0, j1, d
  c=r%cosine; s=r%sine
  ! first and last row operated on
  j0=min(m-l,max(1,k-ubw+1))
  j1=min(m-l,k+lbw)
  do j=j0,j1
     d=k-j+lbw+1
     tmp=b(j,d)
     b(j,d)=c*tmp+s*b(j,d+1)
     b(j,d+1)=-conjg(s)*tmp+c*b(j,d+1)
  end do
end subroutine c_tbr_times_rotation

! A rotation times a truncated band matrix with aligned rows.  The
! rotation is not applied to the first l columns.
subroutine d_rotation_times_tbr(r,b,m,lbw,ubw,l,j)
  real(kind=dp), dimension(:,:), intent(inout) :: b
  type(d_rotation), intent(in) :: r
  integer(kind=int32), intent(in) :: j,m, ubw,lbw, l
  !
  real(kind=dp) :: c, s, tmp
  integer(kind=int32) :: d, k0,k1,k
  c=r%cosine; s=r%sine
  k0=max(l+1,j-lbw+1)
  k1=max(l+1,min(m,j+ubw))
  do k=k0,k1
     d=k-j+lbw+1
     tmp=b(j,d)
     b(j,d)=c*tmp-s*b(j+1,d-1)
     b(j+1,d-1)=s*tmp+c*b(j+1,d-1)
  end do
end subroutine d_rotation_times_tbr

subroutine c_rotation_times_tbr(r,b,m,lbw,ubw,l,j)
  complex(kind=dp), dimension(:,:), intent(inout) :: b
  type(d_rotation), intent(in) :: r
  integer(kind=int32), intent(in) :: j,m, ubw,lbw, l
  !
  complex(kind=dp) :: c, s, tmp
  integer(kind=int32) :: d, k0,k1,k
  c=r%cosine; s=r%sine
  k0=max(l+1,j-lbw+1)
  k1=max(l+1,min(m,j+ubw))
  do k=k0,k1
     d=k-j+lbw+1
     tmp=b(j,d)
     b(j,d)=c*tmp-conjg(s)*b(j+1,d-1)
     b(j+1,d-1)=s*tmp+c*b(j+1,d-1)
  end do
end subroutine c_rotation_times_tbr

end module band
