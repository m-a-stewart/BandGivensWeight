module mod_band_types
  use mod_prec
  use mod_rotation
  use mod_utility
  implicit none
  ! Routines for working with compact representations of banded
  ! matrices.  A banded matrix $A=[a_{jk}]$ with $a_{jk}=0$ for
  ! $k-j>s$ and $k-j < -r$ is stored in one of two ways: In the first
  ! format, elements of $A$ are stored in a matrix $B$ with columns
  ! aligned with those of $A$ and $a_{j,k}=b_{j-k+s+1,k}$. Thus each
  ! row of $B$ corresponds to a diagonal of $A$ and $B$ can be chosen
  ! to be $(r+s+1)\times n$.  Routines working with this format are
  ! given the tag `_bc`. In the second format, $B$ has rows aligned
  ! with those of $A$ and $a_{j,k}=b_{j,k-j+r+1}$.  In this case, each
  ! column of $B$ corresponds to a diagonal of $A$ and $B$ can be
  ! chosen to be $n\times (r+s+1)$.  Routines working with this format
  ! are given the tag `_br`.  The module provides routines for getting
  ! and setting elements $a_{jk}$ of a matrix $A$ stored in this
  ! format; for applying applying rotations in a truncated manner
  ! corresponding to $j$-leading or $j$-trailing transformations; for
  ! printing band matrices in a readable format; and for converting
  ! between different storage formats.

  private

  public :: get_el_bc, d_get_el_bc, z_get_el_bc, &
       set_el_bc, d_set_el_bc, z_set_el_bc, &
       get_el_br, d_get_el_br, z_get_el_br, &
       set_el_br, d_set_el_br, z_set_el_br

  public :: bc_to_general, d_bc_to_general, z_bc_to_general, &
       br_to_general, d_br_to_general, z_br_to_general

  public :: print_bc, d_print_bc, z_print_bc, &
       print_br, d_print_br, z_print_br, &
       print_abs_bc, d_print_abs_bc, z_print_abs_bc, &
       print_abs_br, d_print_abs_br, z_print_abs_br

  public :: f_rotation_times_tbc, f_d_rotation_times_tbc, &
       f_z_rotation_times_tbc, f_tbc_times_rotation, f_d_tbc_times_rotation, &
       f_z_tbc_times_rotation, f_tbr_times_rotation, f_d_tbr_times_rotation, &
       f_z_tbr_times_rotation, f_rotation_times_tbr, f_d_rotation_times_tbr, &
       f_z_rotation_times_tbr

  public :: rotation_times_tbc, d_rotation_times_tbc, z_rotation_times_tbc, &
       tbc_times_rotation, d_tbc_times_rotation, z_tbc_times_rotation, &
       tbr_times_rotation, d_tbr_times_rotation, z_tbr_times_rotation, &
       rotation_times_tbr, d_rotation_times_tbr, z_rotation_times_tbr

  public :: bc_to_br, d_bc_to_br, z_bc_to_br, &
       br_to_bc, d_br_to_bc, z_br_to_bc

  public :: extract_diagonals_bc, d_extract_diagonals_bc, &
       z_extract_diagonals_bc, extract_diagonals_br, d_extract_diagonals_br, &
       z_extract_diagonals_br

  public :: truncate_profile_bc, d_truncate_profile_bc, z_truncate_profile_bc, &
       truncate_profile_br, d_truncate_profile_br, z_truncate_profile_br

  public :: d_band_norm_bc, z_band_norm_bc, band_norm_bc, &
       d_band_norm_br, z_band_norm_br, band_norm_br

  interface get_el_bc
     module procedure d_get_el_bc, z_get_el_bc
  end interface get_el_bc

  interface set_el_bc
     module procedure d_set_el_bc, z_set_el_bc
  end interface set_el_bc

  interface get_el_br
     module procedure d_get_el_br, z_get_el_br
  end interface get_el_br

  interface set_el_br
     module procedure d_set_el_br, z_set_el_br
  end interface set_el_br

  ! to general

  interface bc_to_general
     module procedure d_bc_to_general, z_bc_to_general
  end interface bc_to_general

  interface br_to_general
     module procedure d_br_to_general, z_br_to_general
  end interface br_to_general

  ! printing

  interface print_bc
     module procedure d_print_bc, z_print_bc
  end interface print_bc

  interface print_br
     module procedure d_print_br, z_print_br
  end interface print_br

  interface print_abs_bc
     module procedure d_print_abs_bc, z_print_abs_bc
  end interface print_abs_bc

  interface print_abs_br
     module procedure d_print_abs_br, z_print_abs_br
  end interface print_abs_br


  ! rotations

  interface f_rotation_times_tbc
     module procedure f_d_rotation_times_tbc, f_z_rotation_times_tbc
  end interface f_rotation_times_tbc

  interface f_tbc_times_rotation
     module procedure f_d_tbc_times_rotation, f_z_tbc_times_rotation
  end interface f_tbc_times_rotation

  interface f_tbr_times_rotation
     module procedure f_d_tbr_times_rotation, f_z_tbr_times_rotation
  end interface f_tbr_times_rotation

  interface f_rotation_times_tbr
     module procedure f_d_rotation_times_tbr, f_z_rotation_times_tbr
  end interface f_rotation_times_tbr

  interface rotation_times_tbc
     module procedure d_rotation_times_tbc, z_rotation_times_tbc
  end interface rotation_times_tbc

  interface tbc_times_rotation
     module procedure d_tbc_times_rotation, z_tbc_times_rotation
  end interface tbc_times_rotation

  interface tbr_times_rotation
     module procedure d_tbr_times_rotation, z_tbr_times_rotation
  end interface tbr_times_rotation

  interface rotation_times_tbr
     module procedure d_rotation_times_tbr, z_rotation_times_tbr
  end interface rotation_times_tbr

  ! conversions

  interface bc_to_br
     module procedure d_bc_to_br, z_bc_to_br
  end interface bc_to_br

  interface br_to_bc
     module procedure d_br_to_bc, z_br_to_bc
  end interface br_to_bc

  interface extract_diagonals_br
     module procedure d_extract_diagonals_br, z_extract_diagonals_br
  end interface extract_diagonals_br

  interface extract_diagonals_bc
     module procedure d_extract_diagonals_bc, z_extract_diagonals_bc
  end interface extract_diagonals_bc

  interface truncate_profile_bc
     module procedure d_truncate_profile_bc, z_truncate_profile_bc
  end interface truncate_profile_bc

  interface truncate_profile_br
     module procedure d_truncate_profile_br, z_truncate_profile_br
  end interface truncate_profile_br

  interface band_norm_bc
     module procedure d_band_norm_bc, z_band_norm_bc
  end interface band_norm_bc

  interface band_norm_br
     module procedure d_band_norm_br, z_band_norm_br
  end interface band_norm_br

contains

  ! Get and set functions for band row and band column.

  real(kind=dp) function d_get_el_bc(bc,ubw,j,k) result(x)
    real(kind=dp), dimension(:,:), intent(in) :: bc
    integer(kind=int32), intent(in) :: ubw,j,k
    !
    x=bc(j-k+ubw+1,k)
  end function d_get_el_bc

  complex(kind=dp) function z_get_el_bc(bc,ubw,j,k) result(x)
    complex(kind=dp), dimension(:,:), intent(in) :: bc
    integer(kind=int32), intent(in) :: ubw,j,k
    !
    x=bc(j-k+ubw+1,k)
  end function z_get_el_bc

  subroutine d_set_el_bc(bc,ubw,j,k,x)
    real(kind=dp), dimension(:,:), intent(inout) :: bc
    integer(kind=int32), intent(in) :: ubw,j,k
    real(kind=dp), intent(in) :: x
    !
    bc(j-k+ubw+1,k) = x
  end subroutine d_set_el_bc

  subroutine z_set_el_bc(bc,ubw,j,k,x)
    complex(kind=dp), dimension(:,:), intent(inout) :: bc
    integer(kind=int32), intent(in) :: ubw,j,k
    complex(kind=dp), intent(in) :: x
    !
    bc(j-k+ubw+1,k) = x
  end subroutine z_set_el_bc

  real(kind=dp) function d_get_el_br(br,lbw,j,k) result(x)
    real(kind=dp), dimension(:,:), intent(in) :: br
    integer(kind=int32), intent(in) :: lbw,j,k
    !
    x=br(j,k-j+lbw+1)
  end function d_get_el_br

  complex(kind=dp) function z_get_el_br(br,lbw,j,k) result(x)
    complex(kind=dp), dimension(:,:), intent(in) :: br
    integer(kind=int32), intent(in) :: lbw,j,k
    !
    x=br(j,k-j+lbw+1)
  end function z_get_el_br

  subroutine d_set_el_br(br,lbw,j,k,x)
    real(kind=dp), dimension(:,:), intent(inout) :: br
    integer(kind=int32), intent(in) :: lbw,j,k
    real(kind=dp), intent(in) :: x
    !
    br(j,k-j+lbw+1)=x
  end subroutine d_set_el_br

  subroutine z_set_el_br(br,lbw,j,k,x)
    complex(kind=dp), dimension(:,:), intent(inout) :: br
    integer(kind=int32), intent(in) :: lbw,j,k
    complex(kind=dp), intent(in) :: x
    !
    br(j,k-j+lbw+1)=x
  end subroutine z_set_el_br

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

  subroutine z_bc_to_general(bc,lbw,ubw,a)
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
  end subroutine z_bc_to_general

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

  subroutine z_br_to_general(br,lbw,ubw,a)
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
  end subroutine z_br_to_general

  ! printing

  subroutine d_print_bc(bc,lbw,ubw)
    real(kind=dp), dimension(:,:), intent(in) :: bc
    integer(kind=int32), intent(in) :: ubw, lbw
    !
    real(kind=dp), dimension(size(bc,2),size(bc,2)) :: a
    call d_bc_to_general(bc,lbw,ubw,a)
    call print_matrix(a)
  end subroutine d_print_bc

  subroutine z_print_bc(bc,lbw,ubw)
    complex(kind=dp), dimension(:,:), intent(in) :: bc
    integer(kind=int32), intent(in) :: ubw, lbw
    !
    complex(kind=dp), dimension(size(bc,2),size(bc,2)) :: a
    call z_bc_to_general(bc,lbw,ubw,a)
    call print_matrix(a)
  end subroutine z_print_bc

  subroutine d_print_br(br,lbw,ubw)
    real(kind=dp), dimension(:,:), intent(in) :: br
    integer(kind=int32), intent(in) :: ubw, lbw
    !
    real(kind=dp), dimension(size(br,1),size(br,1)) :: a
    call d_br_to_general(br,lbw,ubw,a)
    call print_matrix(a)
  end subroutine d_print_br

  subroutine z_print_br(br,lbw,ubw)
    complex(kind=dp), dimension(:,:), intent(in) :: br
    integer(kind=int32), intent(in) :: ubw, lbw
    !
    complex(kind=dp), dimension(size(br,1),size(br,1)) :: a
    call z_br_to_general(br,lbw,ubw,a)
    call print_matrix(a)
  end subroutine z_print_br

  ! abs print

  subroutine d_print_abs_bc(bc,lbw,ubw)
    real(kind=dp), dimension(:,:), intent(in) :: bc
    integer(kind=int32), intent(in) :: ubw, lbw
    !
    real(kind=dp), dimension(size(bc,2),size(bc,2)) :: a
    call d_bc_to_general(bc,lbw,ubw,a)
    call print_matrix(abs(a))
  end subroutine d_print_abs_bc

  subroutine z_print_abs_bc(bc,lbw,ubw)
    complex(kind=dp), dimension(:,:), intent(in) :: bc
    integer(kind=int32), intent(in) :: ubw, lbw
    !
    complex(kind=dp), dimension(size(bc,2),size(bc,2)) :: a
    call z_bc_to_general(bc,lbw,ubw,a)
    call print_matrix(abs(a))
  end subroutine z_print_abs_bc

  subroutine d_print_abs_br(br,lbw,ubw)
    real(kind=dp), dimension(:,:), intent(in) :: br
    integer(kind=int32), intent(in) :: ubw, lbw
    !
    real(kind=dp), dimension(size(br,1),size(br,1)) :: a
    call d_br_to_general(br,lbw,ubw,a)
    call print_matrix(abs(a))
  end subroutine d_print_abs_br

  subroutine z_print_abs_br(br,lbw,ubw)
    complex(kind=dp), dimension(:,:), intent(in) :: br
    integer(kind=int32), intent(in) :: ubw, lbw
    !
    complex(kind=dp), dimension(size(br,1),size(br,1)) :: a
    call z_br_to_general(br,lbw,ubw,a)
    call print_matrix(abs(a))
  end subroutine z_print_abs_br

  ! Rotations applied to band matrices.
  !
  ! In applying rotations to a band matrix:
  !
  ! Column k has nonzeros in rows [max(k-ubw,1),min(k+lbw,m)]
  ! Column k+1 has nonzeros in rows [max(k-ubw+1,1),min(k+1+lbw,m)]
  ! A rotation acting from the right in columns k and k+1
  ! acts on the intersection of these row sets: [max(k-ubw+1,1),min(k+1+lbw,m)].
  !
  ! Row j has nonzeros in columns [max(1,j-lbw), min(j+ubw,n)]
  ! Row j+1 has nonzeros in columns [max(1,j+1-lbw), min(j+1+ubw,n)]
  ! A rotation acting on rows j and j+1 acts on the intersection:
  ! [max(1,j+1-lbw),min(n,j+ubw)]

  ! Rotations for band matrices with aligned columns
  !
  ! Rotation times truncated band matrix with aligned columns.
  ! First l and last p columns are not modified.

  subroutine f_d_rotation_times_tbc(c,s,b,n,lbw,ubw,l,p,j)
    real(kind=dp), dimension(:,:), intent(inout) :: b
    real(kind=dp), intent(in) :: c, s
    integer(kind=int32), intent(in) :: j,n, ubw,lbw, l, p
    !
    real(kind=dp) :: tmp
    integer(kind=int32) :: k, k0, k1, d
    ! first and last columns operated on, given by the
    ! intersection of [l+1,n-p] and [max(1,j+1-lbw),min(n,j+ubw)]
    ! which is [maxl+1,max(1,j+1-lbw)), min(n-p,min(n,j+ubw))]
    k0=max(l+1,j+1-lbw)
    k1=min(n-p,j+ubw)
    ! loop over relevant columns in b.
    do k=k0,k1
       d=j-k+ubw+1
       tmp = b(d,k)
       b(d,k)=c*tmp - s*b(d+1,k)
       b(d+1,k) = s*tmp + c*b(d+1,k)
    end do
  end subroutine f_d_rotation_times_tbc

  subroutine d_rotation_times_tbc(r,b,n,lbw,ubw,l,p,j)
    real(kind=dp), dimension(:,:), intent(inout) :: b
    type(d_rotation), intent(in) :: r
    integer(kind=int32), intent(in) :: j,n, ubw,lbw, l, p
    !
    real(kind=dp) :: c, s, tmp
    integer(kind=int32) :: k, k0, k1, d
    c=r%cosine; s=r%sine
    ! first and last columns operated on, given by the
    ! intersection of [l+1,n-p] and [max(1,j+1-lbw),min(n,j+ubw)]
    ! which is [maxl+1,max(1,j+1-lbw)), min(n-p,min(n,j+ubw))]
    k0=max(l+1,j+1-lbw)
    k1=min(n-p,j+ubw)
    ! loop over relevant columns in b.
    do k=k0,k1
       d=j-k+ubw+1
       tmp = b(d,k)
       b(d,k)=c*tmp - s*b(d+1,k)
       b(d+1,k) = s*tmp + c*b(d+1,k)
    end do
  end subroutine d_rotation_times_tbc

  ! Rotation times truncated band matrix with aligned columns.
  ! First l and last p columns are not modified.
  subroutine f_z_rotation_times_tbc(c,s,b,n,lbw,ubw,l,p,j)
    complex(kind=dp), dimension(:,:), intent(inout) :: b
    complex(kind=dp), intent(in) :: s
    real(kind=dp), intent(in) :: c
    integer(kind=int32), intent(in) :: j,n, ubw,lbw, l, p
    !
    complex(kind=dp) :: tmp
    integer(kind=int32) :: k, k0, k1, d
    ! first and last columns operated on, given by the
    ! intersection of [l+1,n-p] and [max(1,j+1-lbw),min(n,j+ubw)]
    ! which is [maxl+1,max(1,j+1-lbw)), min(n-p,min(n,j+ubw))]
    k0=max(l+1,j+1-lbw)
    k1=min(n-p,j+ubw)
    ! loop over relevant columns in b.
    do k=k0,k1
       d=j-k+ubw+1
       tmp = b(d,k)
       b(d,k)=c*tmp - conjg(s)*b(d+1,k)
       b(d+1,k) = s*tmp + c*b(d+1,k)
    end do
  end subroutine f_z_rotation_times_tbc

  ! Rotation times truncated band matrix with aligned columns.
  ! First l and last p columns are not modified.
  subroutine z_rotation_times_tbc(r,b,n,lbw,ubw,l,p,j)
    complex(kind=dp), dimension(:,:), intent(inout) :: b
    type(z_rotation), intent(in) :: r
    integer(kind=int32), intent(in) :: j,n, ubw,lbw, l, p
    !
    complex(kind=dp) :: s, tmp
    real(kind=dp) :: c
    integer(kind=int32) :: k, k0, k1, d
    c=r%cosine; s=r%sine
    ! first and last columns operated on, given by the
    ! intersection of [l+1,n-p] and [max(1,j+1-lbw),min(n,j+ubw)]
    ! which is [maxl+1,max(1,j+1-lbw)), min(n-p,min(n,j+ubw))]
    k0=max(l+1,j+1-lbw)
    k1=min(n-p,j+ubw)
    ! loop over relevant columns in b.
    do k=k0,k1
       d=j-k+ubw+1
       tmp = b(d,k)
       b(d,k)=c*tmp - conjg(s)*b(d+1,k)
       b(d+1,k) = s*tmp + c*b(d+1,k)
    end do
  end subroutine z_rotation_times_tbc

  ! apply a rotation to columns k and k+1 of a truncated band matrix with
  ! aligned columns. first l and last p rows are not modified.
  subroutine f_d_tbc_times_rotation(b,n,lbw,ubw,l,p,c,s,k)
    real(kind=dp), dimension(:,:), intent(inout) :: b
    real(kind=dp), intent(in) :: c, s
    integer(kind=int32), intent(in) :: k,n, ubw,lbw, l, p
    !
    real(kind=dp) :: tmp
    integer(kind=int32) :: d, j0, j1, j
    ! first and last rows operated on are given by the intersection
    ! of [l+1,n-p] with [max(k-ubw+1,1),min(k+lbw,n)] which is
    ! [max(l+1,k-ubw+1),min(k+lbw,n-p)]
    j0=max(l+1,k-ubw+1)
    j1=min(k+lbw,n-p)
    do j=j0,j1
       d=j-k+ubw+1
       tmp=b(d,k)
       b(d,k)=c*tmp+s*b(d-1,k+1)
       b(d-1,k+1)=-s*tmp+c*b(d-1,k+1)
    end do
  end subroutine f_d_tbc_times_rotation

  ! apply a rotation to columns k and k+1 of a truncated band matrix with
  ! aligned columns. first l and last p rows are not modified.
  subroutine d_tbc_times_rotation(b,n,lbw,ubw,l,p,r,k)
    real(kind=dp), dimension(:,:), intent(inout) :: b
    type(d_rotation), intent(in) :: r
    integer(kind=int32), intent(in) :: k,n, ubw,lbw, l, p
    !
    real(kind=dp) :: c, s, tmp
    integer(kind=int32) :: d, j0, j1, j
    c=r%cosine; s=r%sine
    ! first and last rows operated on are given by the intersection
    ! of [l+1,n-p] with [max(k-ubw+1,1),min(k+lbw,n)] which is
    ! [max(l+1,k-ubw+1),min(k+lbw,n-p)]
    j0=max(l+1,k-ubw+1)
    j1=min(k+lbw,n-p)
    do j=j0,j1
       d=j-k+ubw+1
       tmp=b(d,k)
       b(d,k)=c*tmp+s*b(d-1,k+1)
       b(d-1,k+1)=-s*tmp+c*b(d-1,k+1)
    end do
  end subroutine d_tbc_times_rotation

  subroutine f_z_tbc_times_rotation(b,n,lbw,ubw,l,p,c,s,k)
    complex(kind=dp), dimension(:,:), intent(inout) :: b
    complex(kind=dp), intent(in) :: s
    real(kind=dp), intent(in) :: c
    integer(kind=int32), intent(in) :: k,n, ubw,lbw, l, p
    !
    complex(kind=dp) :: tmp
    integer(kind=int32) :: d, j0, j1, j
    ! first and last rows operated on are given by the intersection
    ! of [l+1,n-p] with [max(k-ubw+1,1),min(k+lbw,n)] which is
    ! [max(l+1,k-ubw+1),min(k+lbw,n-p)]
    j0=max(l+1,k-ubw+1)
    j1=min(k+lbw,n-p)
    do j=j0,j1
       d=j-k+ubw+1
       tmp=b(d,k)
       b(d,k)=c*tmp+s*b(d-1,k+1)
       b(d-1,k+1)=-conjg(s)*tmp+c*b(d-1,k+1)
    end do
  end subroutine f_z_tbc_times_rotation

  subroutine z_tbc_times_rotation(b,n,lbw,ubw,l,p,r,k)
    complex(kind=dp), dimension(:,:), intent(inout) :: b
    type(z_rotation), intent(in) :: r
    integer(kind=int32), intent(in) :: k,n, ubw,lbw, l, p
    !
    complex(kind=dp) :: s, tmp
    real(kind=dp) :: c
    integer(kind=int32) :: d, j0, j1, j
    c=r%cosine; s=r%sine
    ! first and last rows operated on are given by the intersection
    ! of [l+1,n-p] with [max(k-ubw+1,1),min(k+lbw,n)] which is
    ! [max(l+1,k-ubw+1),min(k+lbw,n-p)]
    j0=max(l+1,k-ubw+1)
    j1=min(k+lbw,n-p)
    do j=j0,j1
       d=j-k+ubw+1
       tmp=b(d,k)
       b(d,k)=c*tmp+s*b(d-1,k+1)
       b(d-1,k+1)=-conjg(s)*tmp+c*b(d-1,k+1)
    end do
  end subroutine z_tbc_times_rotation

  ! rotations for band matrices with aligned rows.
  !
  ! apply a rotation to columns k and k+1 of a truncated band matrix with
  ! aligned rows. The rotation is not applied to the first l or last p rows.
  subroutine f_d_tbr_times_rotation(b,m,lbw,ubw,l,p,c,s,k)
    real(kind=dp), dimension(:, :), intent(inout) :: b
    real(kind=dp), intent(in) :: c, s
    integer(kind=int32), intent(in) :: k,m, ubw,lbw, l,p
    !
    real(kind=dp) :: tmp
    integer(kind=int32) :: j, j0, j1, d
    ! first and last rows operated on are given by
    ! the intersection of [l+1,m-p] and [max(k-ubw+1,1),min(k+lbw,m)]
    ! which is [max(l+1,k-ubw+1), min(k+lbw,m-p)]
    j0=max(l+1,k-ubw+1)
    j1=min(k+lbw,m-p)
    do j=j0,j1
       d=k-j+lbw+1
       tmp=b(j,d)
       b(j,d)=c*tmp+s*b(j,d+1)
       b(j,d+1)=-s*tmp+c*b(j,d+1)
    end do
  end subroutine f_d_tbr_times_rotation

  ! apply a rotation to columns k and k+1 of a truncated band matrix with
  ! aligned rows. The rotation is not applied to the first l or last p rows.
  subroutine d_tbr_times_rotation(b,m,lbw,ubw,l,p,r,k)
    real(kind=dp), dimension(:, :), intent(inout) :: b
    type(d_rotation), intent(in) :: r
    integer(kind=int32), intent(in) :: k,m, ubw,lbw, l,p
    !
    real(kind=dp) :: c, s, tmp
    integer(kind=int32) :: j, j0, j1, d
    c=r%cosine; s=r%sine
    ! first and last rows operated on are given by
    ! the intersection of [l+1,m-p] and [max(k-ubw+1,1),min(k+lbw,m)]
    ! which is [max(l+1,k-ubw+1), min(k+lbw,m-p)]
    j0=max(l+1,k-ubw+1)
    j1=min(k+lbw,m-p)
    do j=j0,j1
       d=k-j+lbw+1
       tmp=b(j,d)
       b(j,d)=c*tmp+s*b(j,d+1)
       b(j,d+1)=-s*tmp+c*b(j,d+1)
    end do
  end subroutine d_tbr_times_rotation

  subroutine f_z_tbr_times_rotation(b,m,lbw,ubw,l,p,c,s,k)
    complex(kind=dp), dimension(:,:), intent(inout) :: b
    complex(kind=dp), intent(in) :: s
    real(kind=dp), intent(in) :: c
    integer(kind=int32), intent(in) :: k,m, ubw,lbw, l, p
    !
    complex(kind=dp) :: tmp
    integer(kind=int32) :: j, j0, j1, d
    ! first and last rows operated on are given by
    ! the intersection of [l+1,m-p] and [max(k-ubw+1,1),min(k+lbw,m)]
    ! which is [max(l+1,k-ubw+1), min(k+lbw,m-p)]
    j0=max(l+1,k-ubw+1)
    j1=min(k+lbw,m-p)
    do j=j0,j1
       d=k-j+lbw+1
       tmp=b(j,d)
       b(j,d)=c*tmp+s*b(j,d+1)
       b(j,d+1)=-conjg(s)*tmp+c*b(j,d+1)
    end do
  end subroutine f_z_tbr_times_rotation

  subroutine z_tbr_times_rotation(b,m,lbw,ubw,l,p,r,k)
    complex(kind=dp), dimension(:,:), intent(inout) :: b
    type(z_rotation), intent(in) :: r
    integer(kind=int32), intent(in) :: k,m, ubw,lbw, l, p
    !
    complex(kind=dp) :: s, tmp
    real(kind=dp) :: c
    integer(kind=int32) :: j, j0, j1, d
    c=r%cosine; s=r%sine
    ! first and last rows operated on are given by
    ! the intersection of [l+1,m-p] and [max(k-ubw+1,1),min(k+lbw,m)]
    ! which is [max(l+1,k-ubw+1), min(k+lbw,m-p)]
    j0=max(l+1,k-ubw+1)
    j1=min(k+lbw,m-p)
    do j=j0,j1
       d=k-j+lbw+1
       tmp=b(j,d)
       b(j,d)=c*tmp+s*b(j,d+1)
       b(j,d+1)=-conjg(s)*tmp+c*b(j,d+1)
    end do
  end subroutine z_tbr_times_rotation

  ! A rotation times a truncated band matrix with aligned rows.  The
  ! rotation is not applied to the first l or last p columns.
  subroutine f_d_rotation_times_tbr(c,s,b,m,lbw,ubw,l,p,j)
    real(kind=dp), dimension(:,:), intent(inout) :: b
    real(kind=dp), intent(in) :: c, s
    integer(kind=int32), intent(in) :: j,m, ubw,lbw, l, p
    !
    real(kind=dp) :: tmp
    integer(kind=int32) :: d, k0,k1,k
    ! first and last columns operated on, given by the
    ! intersection of [l+1,m-p] and [max(1,j+1-lbw),min(m,j+ubw)]
    ! which is [maxl+1,max(1,j+1-lbw)), min(m-p,min(n,j+ubw))]
    k0=max(l+1,j+1-lbw)
    k1=min(m-p,j+ubw)
    do k=k0,k1
       d=k-j+lbw+1
       tmp=b(j,d)
       b(j,d)=c*tmp-s*b(j+1,d-1)
       b(j+1,d-1)=s*tmp+c*b(j+1,d-1)
    end do
  end subroutine f_d_rotation_times_tbr

  ! A rotation times a truncated band matrix with aligned rows.  The
  ! rotation is not applied to the first l or last p columns.
  subroutine d_rotation_times_tbr(r,b,m,lbw,ubw,l,p,j)
    real(kind=dp), dimension(:,:), intent(inout) :: b
    type(d_rotation), intent(in) :: r
    integer(kind=int32), intent(in) :: j,m, ubw,lbw, l, p
    !
    real(kind=dp) :: c, s, tmp
    integer(kind=int32) :: d, k0,k1,k
    c=r%cosine; s=r%sine
    ! first and last columns operated on, given by the
    ! intersection of [l+1,m-p] and [max(1,j+1-lbw),min(m,j+ubw)]
    ! which is [maxl+1,max(1,j+1-lbw)), min(m-p,min(n,j+ubw))]
    k0=max(l+1,j+1-lbw)
    k1=min(m-p,j+ubw)
    do k=k0,k1
       d=k-j+lbw+1
       tmp=b(j,d)
       b(j,d)=c*tmp-s*b(j+1,d-1)
       b(j+1,d-1)=s*tmp+c*b(j+1,d-1)
    end do
  end subroutine d_rotation_times_tbr

  subroutine f_z_rotation_times_tbr(c,s,b,m,lbw,ubw,l,p,j)
    complex(kind=dp), dimension(:,:), intent(inout) :: b
    complex(kind=dp), intent(in) :: s
    real(kind=dp), intent(in) :: c
    integer(kind=int32), intent(in) :: j,m, ubw,lbw, l,p
    !
    complex(kind=dp) :: tmp
    integer(kind=int32) :: d, k0,k1,k
    ! first and last columns operated on, given by the
    ! intersection of [l+1,m-p] and [max(1,j+1-lbw),min(m,j+ubw)]
    ! which is [maxl+1,max(1,j+1-lbw)), min(m-p,min(n,j+ubw))]
    k0=max(l+1,j+1-lbw)
    k1=min(m-p,j+ubw)
    do k=k0,k1
       d=k-j+lbw+1
       tmp=b(j,d)
       b(j,d)=c*tmp-conjg(s)*b(j+1,d-1)
       b(j+1,d-1)=s*tmp+c*b(j+1,d-1)
    end do
  end subroutine f_z_rotation_times_tbr

  subroutine z_rotation_times_tbr(r,b,m,lbw,ubw,l,p,j)
    complex(kind=dp), dimension(:,:), intent(inout) :: b
    type(z_rotation), intent(in) :: r
    integer(kind=int32), intent(in) :: j,m, ubw,lbw, l,p
    !
    complex(kind=dp) :: s, tmp
    real(kind=dp) :: c
    integer(kind=int32) :: d, k0,k1,k
    c=r%cosine; s=r%sine
    ! first and last columns operated on, given by the
    ! intersection of [l+1,m-p] and [max(1,j+1-lbw),min(m,j+ubw)]
    ! which is [maxl+1,max(1,j+1-lbw)), min(m-p,min(n,j+ubw))]
    k0=max(l+1,j+1-lbw)
    k1=min(m-p,j+ubw)
    do k=k0,k1
       d=k-j+lbw+1
       tmp=b(j,d)
       b(j,d)=c*tmp-conjg(s)*b(j+1,d-1)
       b(j+1,d-1)=s*tmp+c*b(j+1,d-1)
    end do
  end subroutine z_rotation_times_tbr

  subroutine d_bc_to_br(bc,br,lbw,ubw)
    real(kind=dp), dimension(:,:), intent(in) :: bc
    real(kind=dp), dimension(:,:), intent(out) :: br
    integer(kind=int32), intent(in) :: lbw, ubw

    integer(kind=int32) :: d, j, n
    n=size(bc,2)
    do d=1,lbw+1
       do j=1,n-lbw+d-1
          br(lbw-d+j+1,d)=bc(ubw+lbw+2-d,j)
       end do
    end do
    do d=lbw+2,lbw+ubw+1
       do j=1,n-d+lbw+1
          br(j,d)=bc(ubw+lbw+2-d,j+d-lbw-1)
       end do
    end do
  end subroutine d_bc_to_br

  subroutine z_bc_to_br(bc,br,lbw,ubw)
    complex(kind=dp), dimension(:,:), intent(in) :: bc
    complex(kind=dp), dimension(:,:), intent(out) :: br
    integer(kind=int32), intent(in) :: lbw, ubw

    integer(kind=int32) :: d, j, n
    n=size(bc,2)
    do d=1,lbw+1
       do j=1,n-lbw+d-1
          br(lbw-d+j+1,d)=bc(ubw+lbw+2-d,j)
       end do
    end do
    do d=lbw+2,lbw+ubw+1
       do j=1,n-d+lbw+1
          br(j,d)=bc(ubw+lbw+2-d,j+d-lbw-1)
       end do
    end do
  end subroutine z_bc_to_br

  subroutine d_br_to_bc(br,bc,lbw,ubw)
    real(kind=dp), dimension(:,:), intent(out) :: bc
    real(kind=dp), dimension(:,:), intent(in) :: br
    integer(kind=int32), intent(in) :: lbw, ubw

    integer(kind=int32) :: d, j, n

    n=size(bc,2)
    do d=1,lbw+1
       do j=1,n-lbw+d-1
          bc(ubw+lbw+2-d,j)=br(lbw-d+j+1,d)
       end do
    end do
    do d=lbw+2,lbw+ubw+1
       do j=1,n-d+lbw+1
          bc(ubw+lbw+2-d,j+d-lbw-1)=br(j,d)
       end do
    end do
  end subroutine d_br_to_bc

  subroutine z_br_to_bc(br,bc,lbw,ubw)
    complex(kind=dp), dimension(:,:), intent(out) :: bc
    complex(kind=dp), dimension(:,:), intent(in) :: br
    integer(kind=int32), intent(in) :: lbw, ubw

    integer(kind=int32) :: d, j, n
    n=size(bc,2)
    do d=1,lbw+1
       do j=1,n-lbw+d-1
          bc(ubw+lbw+2-d,j)=br(lbw-d+j+1,d)
       end do
    end do
    do d=lbw+2,lbw+ubw+1
       do j=1,n-d+lbw+1
          bc(ubw+lbw+2-d,j+d-lbw-1)=br(j,d)
       end do
    end do
  end subroutine z_br_to_bc

  subroutine d_extract_diagonals_bc(a, n, b, lbw, ubw, lbwmax, ubwmax)
    real(kind=dp), target, dimension(n,n), intent(in) :: a
    real(kind=dp), dimension(lbwmax+ubwmax+1,n), intent(out) :: b
    integer(kind=int32), intent(in) :: n, lbw, ubw, lbwmax, ubwmax
    !
    integer(kind=int32) :: k, d
    ! put diagonals in b
    b=0.0_dp
    do d=1,ubw+1
       do k=ubw-d+2,n
          b(d,k) = a(d+k-ubw-1,k)
       end do
    end do
    do d=ubw+2, ubw+lbw+1
       do k=1,n-d+ubw+1
          b(d,k) = a(d+k-ubw-1,k)
       end do
    end do
  end subroutine d_extract_diagonals_bc

  subroutine z_extract_diagonals_bc(a, n, b, lbw, ubw, lbwmax, ubwmax)
    complex(kind=dp), target, dimension(n,n), intent(in) :: a
    complex(kind=dp), dimension(lbwmax+ubwmax+1,n), intent(out) :: b
    integer(kind=int32), intent(in) :: n, lbw, ubw, lbwmax, ubwmax
    !
    integer(kind=int32) :: k, d
    ! put diagonals in b
    b=0.0_dp
    do d=1,ubw+1
       do k=ubw-d+2,n
          b(d,k) = a(d+k-ubw-1,k)
       end do
    end do
    do d=ubw+2, ubw+lbw+1
       do k=1,n-d+ubw+1
          b(d,k) = a(d+k-ubw-1,k)
       end do
    end do
  end subroutine z_extract_diagonals_bc

  subroutine d_extract_diagonals_br(a, n, b, lbw, ubw, lbwmax, ubwmax)
    real(kind=dp), target, dimension(n,n), intent(in) :: a
    real(kind=dp), dimension(n,lbwmax+ubwmax+1), intent(out) :: b
    integer(kind=int32), intent(in) :: n, lbw, ubw, lbwmax, ubwmax
    !
    integer(kind=int32) :: j, d
    ! put diagonals in b
    b=0.0_dp
    do d=1,lbw+1
       do j=lbw-d+2, n
          b(j,d)=a(j,d+j-lbw-1)
       end do
    end do
    do d=lbw+2,ubw+lbw+1
       do j=1, n-d+lbw+1
          b(j,d)=a(j,d+j-lbw-1)
       end do
    end do
  end subroutine d_extract_diagonals_br

  subroutine z_extract_diagonals_br(a, n, b, lbw, ubw, lbwmax, ubwmax)
    complex(kind=dp), target, dimension(n,n), intent(in) :: a
    complex(kind=dp), dimension(n,lbwmax+ubwmax+1), intent(out) :: b
    integer(kind=int32), intent(in) :: n, lbw, ubw, lbwmax, ubwmax
    !
    integer(kind=int32) :: j, d
    ! put diagonals in b
    b=0.0_dp
    do d=1,lbw+1
       do j=lbw-d+2, n
          b(j,d)=a(j,d+j-lbw-1)
       end do
    end do
    do d=lbw+2,ubw+lbw+1
       do j=1, n-d+lbw+1
          b(j,d)=a(j,d+j-lbw-1)
       end do
    end do
  end subroutine z_extract_diagonals_br

  ! Impose a particular upper and lower variable bandwidth on a band
  ! matrix.  Used for generating random band matrices for tests.
  subroutine d_truncate_profile_br(br,n,lbw,lbwmax,ubwmax,lower,upper)
    real(kind=dp), dimension(n,lbwmax+ubwmax+1), intent(inout) :: br
    integer(kind=int32), intent(in) :: n, lbw, lbwmax, ubwmax
    integer(kind=int32), dimension(n), intent(in), optional :: lower, upper

    integer(kind=int32) :: j

    if (present(upper)) then
       do j=1,n
          br(j,lbw+2+upper(j):lbwmax+ubwmax+1) = 0.0_dp
       end do
    end if

    if (present(lower)) then
       do j=1,n
          br(j,1:lbw-lower(j)) = 0.0_dp
       end do
    end if

  end subroutine d_truncate_profile_br

  subroutine z_truncate_profile_br(br,n,lbw,lbwmax,ubwmax,lower,upper)
    complex(kind=dp), dimension(n,lbwmax+ubwmax+1), intent(inout) :: br
    integer(kind=int32), intent(in) :: n, lbw, lbwmax, ubwmax
    integer(kind=int32), dimension(n), intent(in), optional :: lower, upper

    integer(kind=int32) :: j

    if (present(upper)) then
       do j=1,n
          br(j,lbw+2+upper(j):lbwmax+ubwmax+1) = (0.0_dp,0.0_dp)
       end do
    end if

    if (present(lower)) then
       do j=1,n
          br(j,1:lbw-lower(j)) = (0.0_dp,0.0_dp)
       end do
    end if
  end subroutine z_truncate_profile_br

  subroutine d_truncate_profile_bc(bc,n,ubw,lbwmax,ubwmax,lower,upper)
    real(kind=dp), dimension(lbwmax+ubwmax+1,n), intent(inout) :: bc
    integer(kind=int32), intent(in) :: n, ubw, lbwmax, ubwmax
    integer(kind=int32), dimension(n), intent(in), optional :: lower, upper

    integer(kind=int32) :: j

    if (present(upper)) then
       do j=1,n
          bc(1:ubw-upper(j),j)=0.0_dp
       end do
    end if

    if (present(lower)) then
       do j=1,n
          bc(ubw+2+lower(j):lbwmax+ubwmax+1,j)=0.0_dp
       end do
    end if

  end subroutine d_truncate_profile_bc

  subroutine z_truncate_profile_bc(bc,n,ubw,lbwmax,ubwmax,lower,upper)
    complex(kind=dp), dimension(lbwmax+ubwmax+1,n), intent(inout) :: bc
    integer(kind=int32), intent(in) :: n, ubw, lbwmax, ubwmax
    integer(kind=int32), dimension(n), intent(in), optional :: lower, upper

    integer(kind=int32) :: j

    if (present(upper)) then
       do j=1,n
          bc(1:ubw-upper(j),j)=(0.0_dp,0.0_dp)
       end do
    end if

    if (present(lower)) then
       do j=1,n
          bc(ubw+2+lower(j):lbwmax+ubwmax+1,j)=(0.0_dp,0.0_dp)
       end do
    end if

  end subroutine z_truncate_profile_bc

  real(kind=dp) function d_band_norm_bc(bc,n,lbw,ubw,lbwmax, &
       ubwmax,d0,d1) result(x)
    real(kind=dp), dimension(lbwmax+ubwmax+1,n), intent(in) :: bc
    integer(kind=int32), intent(in) :: n, lbw, ubw, lbwmax, ubwmax, d0, d1

    integer(kind=int32) :: k,l,l0,l1
    real(kind=dp) :: y, tmp

    l0 = ubw+1-d1
    l1 = ubw+1-d0

    y=0.0_dp
    ! maxabs for the specified band.
    do l=max(1,l0),min(ubw+1,l1)
       do k=ubw-l+2,n
          tmp=abs(bc(l,k))
          if (tmp > y) y=tmp
       end do
    end do
    do l=max(l0,ubw+2), min(ubw+lbw+1,l1)
       do k=1,n-l+ubw+1
          tmp=abs(bc(l,k))
          if (tmp > y) y=tmp
       end do
    end do
    x=0.0_dp;
    if (y==0.0_dp) return
    y=(2.0_dp)**(floor(log(y)/log(2.0_dp)))
    y=1.0_dp/y
    if (y > 1e150_dp .or. y < 1e-150_dp) then
       do l=max(1,l0),min(ubw+1,l1)
          do k=ubw-l+2,n
             tmp=bc(l,k)*y
             x=x+tmp*tmp
          end do
       end do
       do l=max(l0,ubw+2), min(ubw+lbw+1,l1)
          do k=1,n-l+ubw+1
             tmp=bc(l,k)*y
             x=x+tmp*tmp
          end do
       end do
       x=sqrt(x)/y
    else
       do l=max(1,l0),min(ubw+1,l1)
          do k=ubw-l+2,n
             tmp=bc(l,k)
             x=x+tmp*tmp
          end do
       end do
       do l=max(l0,ubw+2), min(ubw+lbw+1,l1)
          do k=1,n-l+ubw+1
             tmp=bc(l,k)
             x=x+tmp*tmp
          end do
       end do
       x=sqrt(x)
    end if
  end function d_band_norm_bc

  real(kind=dp) function z_band_norm_bc(bc,n,lbw,ubw,lbwmax, &
       ubwmax,d0,d1) result(x)
    complex(kind=dp), dimension(lbwmax+ubwmax+1,n), intent(in) :: bc
    integer(kind=int32), intent(in) :: n, lbw, ubw, lbwmax, ubwmax, d0, d1

    integer(kind=int32) :: k,l,l0,l1
    real(kind=dp) :: y, tmp

    l0 = ubw+1-d1
    l1 = ubw+1-d0

    y=0.0_dp
    ! maxabs for the specified band.
    do l=max(1,l0),min(ubw+1,l1)
       do k=ubw-l+2,n
          tmp=abs(bc(l,k))
          if (tmp > y) y=tmp
       end do
    end do
    do l=max(l0,ubw+2), min(ubw+lbw+1,l1)
       do k=1,n-l+ubw+1
          tmp=abs(bc(l,k))
          if (tmp > y) y=tmp
       end do
    end do

    x=0.0_dp;
    if (y==0.0_dp) return
    y=(2.0_dp)**(floor(log(y)/log(2.0_dp)))
    y=1.0_dp/y
    if (y > 1e150_dp .or. y < 1e-150_dp) then
       do l=max(1,l0),min(ubw+1,l1)
          do k=ubw-l+2,n
             tmp=abs(bc(l,k)*y)
             x=x+tmp*tmp
          end do
       end do
       do l=max(l0,ubw+2), min(ubw+lbw+1,l1)
          do k=1,n-l+ubw+1
             tmp=abs(bc(l,k)*y)
             x=x+tmp*tmp
          end do
       end do
       x=sqrt(x)/y
    else
       do l=max(1,l0),min(ubw+1,l1)
          do k=ubw-l+2,n
             tmp=abs(bc(l,k))
             x=x+tmp*tmp
          end do
       end do
       do l=max(l0,ubw+2), min(ubw+lbw+1,l1)
          do k=1,n-l+ubw+1
             tmp=abs(bc(l,k))
             x=x+tmp*tmp
          end do
       end do
       x=sqrt(x)
    end if
  end function z_band_norm_bc

  real(kind=dp) function d_band_norm_br(br,n,lbw,ubw,lbwmax, &
       ubwmax,d0,d1) result(x)
    real(kind=dp), dimension(n,lbwmax+ubwmax+1), intent(in) :: br
    integer(kind=int32), intent(in) :: n, lbw, ubw, lbwmax, ubwmax, d0, d1

    integer(kind=int32) :: j,l,l0,l1
    real(kind=dp) :: y, tmp

    l0 = d0+lbw+1
    l1 = d1+lbw+1

    y=0.0_dp
    ! maxabs for the specified band.

    do l=max(1,l0),min(lbw+1,l1)
       do j=lbw-l+2,n
          tmp=abs(br(j,l))
          if (tmp > y) y=tmp
       end do
    end do
    do l=max(l0,lbw+2), min(ubw+lbw+1,l1)
       do j=1,n-l+lbw+1
          tmp=abs(br(j,l))
          if (tmp > y) y=tmp
       end do
    end do

    x=0.0_dp;
    if (y==0.0_dp) return
    y=(2.0_dp)**(floor(log(y)/log(2.0_dp)))
    y=1.0_dp/y
    if (y > 1e150_dp .or. y < 1e-150_dp) then
       do l=max(1,l0),min(lbw+1,l1)
          do j=lbw-l+2,n
             tmp=br(j,l)*y
             x=x+tmp*tmp
          end do
       end do
       do l=max(l0,lbw+2), min(ubw+lbw+1,l1)
          do j=1,n-l+lbw+1
             tmp=br(j,l)*y
             x=x+tmp*tmp
          end do
       end do
       x=sqrt(x)/y
    else
       do l=max(1,l0),min(lbw+1,l1)
          do j=lbw-l+2,n
             tmp=br(j,l)
             x=x+tmp*tmp
          end do
       end do
       do l=max(l0,lbw+2), min(ubw+lbw+1,l1)
          do j=1,n-l+lbw+1
             tmp=br(j,l)
             x=x+tmp*tmp
          end do
       end do
       x=sqrt(x)
    end if
  end function d_band_norm_br

  real(kind=dp) function z_band_norm_br(br,n,lbw,ubw,lbwmax, &
       ubwmax,d0,d1) result(x)
    complex(kind=dp), dimension(n,lbwmax+ubwmax+1), intent(in) :: br
    integer(kind=int32), intent(in) :: n, lbw, ubw, lbwmax, ubwmax, d0, d1

    integer(kind=int32) :: j,l,l0,l1
    real(kind=dp) :: y, tmp

    l0 = d0+lbw+1
    l1 = d1+lbw+1

    y=0.0_dp
    ! maxabs for the specified band.
    do l=max(1,l0),min(lbw+1,l1)
       do j=lbw-l+2,n
          tmp=abs(br(j,l))
          if (tmp > y) y=tmp
       end do
    end do
    do l=max(l0,lbw+2), min(ubw+lbw+1,l1)
       do j=1,n-l+lbw+1
          tmp=abs(br(j,l))
          if (tmp > y) y=tmp
       end do
    end do

    x=0.0_dp;
    if (y==0.0_dp) return
    y=(2.0_dp)**(floor(log(y)/log(2.0_dp)))
    y=1.0_dp/y
    if (y > 1e150_dp .or. y < 1e-150_dp) then
       do l=max(1,l0),min(lbw+1,l1)
          do j=lbw-l+2,n
             tmp=abs(br(j,l)*y)
             x=x+tmp*tmp
          end do
       end do
       do l=max(l0,lbw+2), min(ubw+lbw+1,l1)
          do j=1,n-l+lbw+1
             tmp=abs(br(j,l)*y)
             x=x+tmp*tmp
          end do
       end do
       x=sqrt(x)/y
    else
       do l=max(1,l0),min(lbw+1,l1)
          do j=lbw-l+2,n
             tmp=abs(br(j,l))
             x=x+tmp*tmp
          end do
       end do
       do l=max(l0,lbw+2), min(ubw+lbw+1,l1)
          do j=1,n-l+lbw+1
             tmp=abs(br(j,l))
             x=x+tmp*tmp
          end do
       end do
       x=sqrt(x)
    end if
  end function z_band_norm_br

end module mod_band_types
