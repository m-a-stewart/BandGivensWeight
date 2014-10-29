module mod_band_types
  use mod_prec
  use mod_rotation
  use mod_shift
  use mod_utility
  implicit none

  private

  public :: get_el_bc, d_get_el_bc, c_get_el_bc, &
       set_el_bc, d_set_el_bc, c_set_el_bc, &
       get_el_br, d_get_el_br, c_get_el_br, &
       set_el_br, d_set_el_br, c_set_el_br

  public :: bc_to_general, d_bc_to_general, c_bc_to_general, &
       br_to_general, d_br_to_general, c_br_to_general

  public :: submatrix_bc, d_submatrix_bc, c_submatrix_bc, &
       submatrix_br, d_submatrix_br, c_submatrix_br

  public :: print_bc, d_print_bc, c_print_bc, &
       print_br, d_print_br, c_print_br, &
       print_abs_bc, d_print_abs_bc, c_print_abs_bc, &
       print_abs_br, d_print_abs_br, c_print_abs_br

  public :: rotation_times_tbc, d_rotation_times_tbc, c_rotation_times_tbc, &
       tbc_times_rotation, d_tbc_times_rotation, c_tbc_times_rotation, &
       tbr_times_rotation, d_tbr_times_rotation, c_tbr_times_rotation, &
       rotation_times_tbr, d_rotation_times_tbr, c_rotation_times_tbr

  public :: bc_to_br, d_bc_to_br, c_bc_to_br, &
       br_to_bc, d_br_to_bc, c_br_to_bc, &
       f_bw_expand_bc, f_d_bw_expand_bc, f_c_bw_expand_bc, &
       f_bw_expand_br, f_d_bw_expand_br, f_c_bw_expand_br, &
       f_bw_contract_bc, f_d_bw_contract_bc, f_c_bw_contract_bc, &
       f_bw_contract_br, f_d_bw_contract_br, f_c_bw_contract_br

  public :: extract_diagonals_bc, d_extract_diagonals_bc, c_extract_diagonals_bc, &
       extract_diagonals_br, d_extract_diagonals_br, c_extract_diagonals_br

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

  ! submatrices

  interface submatrix_bc
     module procedure d_submatrix_bc, c_submatrix_bc
  end interface submatrix_bc

  interface submatrix_br
     module procedure d_submatrix_br, c_submatrix_br
  end interface submatrix_br

  ! printing

  interface print_bc
     module procedure d_print_bc, c_print_bc
  end interface print_bc

  interface print_br
     module procedure d_print_br, c_print_br
  end interface print_br

  interface print_abs_bc
     module procedure d_print_abs_bc, c_print_abs_bc
  end interface print_abs_bc

  interface print_abs_br
     module procedure d_print_abs_br, c_print_abs_br
  end interface print_abs_br


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

  ! conversions

  interface bc_to_br
     module procedure d_bc_to_br, c_bc_to_br
  end interface bc_to_br

  interface br_to_bc
     module procedure d_br_to_bc, c_br_to_bc
  end interface br_to_bc

  interface f_bw_expand_bc
     module procedure f_d_bw_expand_bc, f_c_bw_expand_bc
  end interface f_bw_expand_bc

  interface f_bw_expand_br
     module procedure f_d_bw_expand_br, f_c_bw_expand_br
  end interface f_bw_expand_br

  interface f_bw_contract_bc
     module procedure f_d_bw_contract_bc, f_c_bw_contract_bc
  end interface f_bw_contract_bc

  interface f_bw_contract_br
     module procedure f_d_bw_contract_br, f_c_bw_contract_br
  end interface f_bw_contract_br

  interface extract_diagonals_br
     module procedure d_extract_diagonals_br, c_extract_diagonals_br 
  end interface extract_diagonals_br

  interface extract_diagonals_bc
     module procedure d_extract_diagonals_bc, c_extract_diagonals_bc 
  end interface extract_diagonals_bc

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

  ! Submatrix

  subroutine d_submatrix_bc(bc,lbw,ubw,j1,j2,k1,k2,a)
    real(kind=dp), dimension(:,:), intent(in) :: bc
    integer(kind=int32), intent(in) :: ubw, lbw,j1,j2,k1,k2
    real(kind=dp), dimension(:,:), intent(out) :: a
    ! 
    integer(kind=int32) :: n, k, j
    a=0.0_dp
    n=size(bc,2)
    do k=k1,k2
       do j=max(k-ubw,j1),min(k+lbw,j2)
          a(j-j1+1,k-k1+1)=bc(j-k+ubw+1,k)
       end do
    end do
  end subroutine d_submatrix_bc

  subroutine c_submatrix_bc(bc,lbw,ubw,j1,j2,k1,k2,a)
    complex(kind=dp), dimension(:,:), intent(in) :: bc
    integer(kind=int32), intent(in) :: ubw, lbw,j1,j2,k1,k2
    complex(kind=dp), dimension(:,:), intent(out) :: a
    ! 
    integer(kind=int32) :: n, k, j
    a=(0.0_dp,0.0_dp)
    n=size(bc,2)
    do k=k1,k2
       do j=max(k-ubw,j1),min(k+lbw,j2)
          a(j-j1+1,k-k1+1)=bc(j-k+ubw+1,k)
       end do
    end do
  end subroutine c_submatrix_bc

  subroutine d_submatrix_br(br,lbw,ubw,j1,j2,k1,k2,a)
    real(kind=dp), dimension(:,:), intent(in) :: br
    integer(kind=int32), intent(in) :: ubw, lbw,j1,j2,k1,k2
    real(kind=dp), dimension(:,:), intent(out) :: a
    ! 
    integer(kind=int32) :: n, k, j
    a=0.0_dp
    n=size(br,2)
    do j=j1,j2
       do k=max(j-lbw,k1),min(j+ubw,k2)
          a(j-j1+1,k-k1+1)=br(j,k-j+lbw+1)
       end do
    end do
  end subroutine d_submatrix_br

  subroutine c_submatrix_br(br,lbw,ubw,j1,j2,k1,k2,a)
    complex(kind=dp), dimension(:,:), intent(in) :: br
    integer(kind=int32), intent(in) :: ubw, lbw,j1,j2,k1,k2
    complex(kind=dp), dimension(:,:), intent(out) :: a
    ! 
    integer(kind=int32) :: n, k, j
    a=(0.0_dp,0.0_dp)
    n=size(br,2)
    do j=j1,j2
       do k=max(j-lbw,k1),min(j+ubw,k2)
          a(j-j1+1,k-k1+1)=br(j,k-j+lbw+1)
       end do
    end do
  end subroutine c_submatrix_br

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

  ! abs print

  subroutine d_print_abs_bc(bc,lbw,ubw)
    real(kind=dp), dimension(:,:), intent(in) :: bc
    integer(kind=int32), intent(in) :: ubw, lbw
    !
    real(kind=dp), dimension(size(bc,2),size(bc,2)) :: a
    call d_bc_to_general(bc,lbw,ubw,a)
    call print_matrix(abs(a))
  end subroutine d_print_abs_bc

  subroutine c_print_abs_bc(bc,lbw,ubw)
    complex(kind=dp), dimension(:,:), intent(in) :: bc
    integer(kind=int32), intent(in) :: ubw, lbw
    !
    complex(kind=dp), dimension(size(bc,2),size(bc,2)) :: a
    call c_bc_to_general(bc,lbw,ubw,a)
    call print_matrix(abs(a))
  end subroutine c_print_abs_bc

  subroutine d_print_abs_br(br,lbw,ubw)
    real(kind=dp), dimension(:,:), intent(in) :: br
    integer(kind=int32), intent(in) :: ubw, lbw
    !
    real(kind=dp), dimension(size(br,1),size(br,1)) :: a
    call d_br_to_general(br,lbw,ubw,a)
    call print_matrix(abs(a))
  end subroutine d_print_abs_br

  subroutine c_print_abs_br(br,lbw,ubw)
    complex(kind=dp), dimension(:,:), intent(in) :: br
    integer(kind=int32), intent(in) :: ubw, lbw
    !
    complex(kind=dp), dimension(size(br,1),size(br,1)) :: a
    call c_br_to_general(br,lbw,ubw,a)
    call print_matrix(abs(a))
  end subroutine c_print_abs_br

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
  ! First l columns are not modified.
  subroutine c_rotation_times_tbc(r,b,n,lbw,ubw,l,p,j)
    complex(kind=dp), dimension(:,:), intent(inout) :: b
    type(c_rotation), intent(in) :: r
    integer(kind=int32), intent(in) :: j,n, ubw,lbw, l, p
    !
    complex(kind=dp) :: c, s, tmp
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
  end subroutine c_rotation_times_tbc

  ! apply a rotation to columns k and k+1 of a truncated band matrix with aligned columns.
  ! first l and last p rows are not modified.
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

  subroutine c_tbc_times_rotation(b,n,lbw,ubw,l,p,r,k)
    complex(kind=dp), dimension(:,:), intent(inout) :: b
    type(c_rotation), intent(in) :: r
    integer(kind=int32), intent(in) :: k,n, ubw,lbw, l, p
    !
    complex(kind=dp) :: c, s, tmp
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
  end subroutine c_tbc_times_rotation

  ! rotations for band matrices with aligned rows.
  !
  ! apply a rotation to columns k and k+1 of a truncated band matrix with aligned rows.
  ! The rotation is not applied to the first l or last p rows.
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

  subroutine c_tbr_times_rotation(b,m,lbw,ubw,l,p,r,k)
    complex(kind=dp), dimension(:,:), intent(inout) :: b
    type(c_rotation), intent(in) :: r
    integer(kind=int32), intent(in) :: k,m, ubw,lbw, l, p
    !
    complex(kind=dp) :: c, s, tmp
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
  end subroutine c_tbr_times_rotation

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

  subroutine c_rotation_times_tbr(r,b,m,lbw,ubw,l,p,j)
    complex(kind=dp), dimension(:,:), intent(inout) :: b
    type(c_rotation), intent(in) :: r
    integer(kind=int32), intent(in) :: j,m, ubw,lbw, l,p
    !
    complex(kind=dp) :: c, s, tmp
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
  end subroutine c_rotation_times_tbr

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

  subroutine c_bc_to_br(bc,br,lbw,ubw)
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
  end subroutine c_bc_to_br

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

  subroutine c_br_to_bc(br,bc,lbw,ubw)
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
  end subroutine c_br_to_bc

  ! expand lbw and ubw with permanent and temporary extra diagonals.
  ! dubw, dubw_tmp, dlbw, dlbw_tmp give requested values.
  ! On output they contain actual values, constrained by the size of the matrix.
  subroutine f_d_bw_expand_bc(b_bc, n, lbw_bc, ubw_bc, lbwmax_bc, ubwmax_bc, &
       dlbw, dlbw_tmp, dubw, dubw_tmp, lbw_br, ubw_br)
    integer(kind=int32), intent(in) :: n, lbwmax_bc, ubwmax_bc
    integer(kind=int32), intent(inout) :: lbw_bc, ubw_bc, dubw, dlbw, dubw_tmp, dlbw_tmp
    real(kind=dp), dimension(lbwmax_bc+ubwmax_bc+1,n), intent(out) :: b_bc
    integer(kind=int32), intent(out) :: lbw_br, ubw_br


    integer(kind=int32) :: j

    ! add extra subdiagonals and put in zeros as needed.
    if (lbw_bc+dlbw >= n-1) then
       dlbw=n-1 - lbw_bc
       dlbw_tmp=0
    else if (lbw_bc+dlbw+dlbw_tmp >= n-1) then
       dlbw_tmp = n-1 - lbw_bc - dlbw
    end if
    if (dlbw+dlbw_tmp > 0) then
       b_bc(lbw_bc+ubw_bc+2:lbw_bc+ubw_bc+1+dlbw + dlbw_tmp,:)=0.0_dp
    end if
    lbw_br=lbw_bc+dlbw
    lbw_bc=lbw_bc+dlbw+dlbw_tmp

    ! add extra superdiagonals and shift as needed.
    if (ubw_bc+dubw >= n-1) then
       dubw=n-1 - ubw_bc
       dubw_tmp=0
    else if (ubw_bc+dubw+dubw_tmp >= n-1) then
       dubw_tmp = n-1 - ubw_bc - dubw
    end if
    do j=1,dubw+dubw_tmp
       call shift2(b_bc,1,0)
    end do
    ubw_br=ubw_bc+dubw
    ubw_bc=ubw_bc+dubw+dubw_tmp

  end subroutine f_d_bw_expand_bc

  subroutine f_c_bw_expand_bc(b_bc, n, lbw_bc, ubw_bc, lbwmax_bc, ubwmax_bc, &
       dlbw, dlbw_tmp, dubw, dubw_tmp, lbw_br, ubw_br)
    integer(kind=int32), intent(in) :: n, lbwmax_bc, ubwmax_bc
    integer(kind=int32), intent(inout) :: lbw_bc, ubw_bc, dubw, dlbw, dubw_tmp, dlbw_tmp
    complex(kind=dp), dimension(lbwmax_bc+ubwmax_bc+1,n), intent(out) :: b_bc
    integer(kind=int32), intent(out) :: lbw_br, ubw_br


    integer(kind=int32) :: j

    ! add extra subdiagonals and put in zeros as needed.
    if (lbw_bc+dlbw >= n-1) then
       dlbw=n-1 - lbw_bc
       dlbw_tmp=0
    else if (lbw_bc+dlbw+dlbw_tmp >= n-1) then
       dlbw_tmp = n-1 - lbw_bc - dlbw
    end if
    if (dlbw+dlbw_tmp > 0) then
       b_bc(lbw_bc+ubw_bc+2:lbw_bc+ubw_bc+1+dlbw + dlbw_tmp,:)=(0.0_dp, 0.0_dp)
    end if
    lbw_br=lbw_bc+dlbw
    lbw_bc=lbw_bc+dlbw+dlbw_tmp

    ! add extra superdiagonals and shift as needed.
    if (ubw_bc+dubw >= n-1) then
       dubw=n-1 - ubw_bc
       dubw_tmp=0
    else if (ubw_bc+dubw+dubw_tmp >= n-1) then
       dubw_tmp = n-1 - ubw_bc - dubw
    end if
    do j=1,dubw+dubw_tmp
       call shift2(b_bc,1,0)
    end do
    ubw_br=ubw_bc+dubw
    ubw_bc=ubw_bc+dubw+dubw_tmp

  end subroutine f_c_bw_expand_bc

  subroutine f_d_bw_contract_bc(b_bc, n, lbw_bc, ubw_bc, lbwmax_bc, ubwmax_bc, &
       dlbw, dlbw_tmp, dubw, dubw_tmp)
    integer(kind=int32), intent(in) :: n, lbwmax_bc, ubwmax_bc
    integer(kind=int32), intent(inout) :: lbw_bc, ubw_bc, dubw, dlbw, dubw_tmp, dlbw_tmp
    real(kind=dp), dimension(lbwmax_bc+ubwmax_bc+1,n), intent(out) :: b_bc

    integer(kind=int32) :: j

    lbw_bc=lbw_bc - dlbw_tmp
    ubw_bc=ubw_bc-dubw_tmp
    do j=1,dubw_tmp
       call shift2(b_bc,-1,0)
    end do

  end subroutine f_d_bw_contract_bc

  subroutine f_c_bw_contract_bc(b_bc, n, lbw_bc, ubw_bc, lbwmax_bc, ubwmax_bc, &
       dlbw, dlbw_tmp, dubw, dubw_tmp)
    integer(kind=int32), intent(in) :: n, lbwmax_bc, ubwmax_bc
    integer(kind=int32), intent(inout) :: lbw_bc, ubw_bc, dubw, dlbw, dubw_tmp, dlbw_tmp
    complex(kind=dp), dimension(lbwmax_bc+ubwmax_bc+1,n), intent(out) :: b_bc

    integer(kind=int32) :: j

    lbw_bc=lbw_bc - dlbw_tmp
    ubw_bc=ubw_bc-dubw_tmp
    do j=1,dubw_tmp
       call shift2(b_bc,-1,0)
    end do

  end subroutine f_c_bw_contract_bc

  ! expand lbw and ubw with permanent and temporary extra diagonals.
  ! dubw, dubw_tmp, dlbw, dlbw_tmp give requested values.
  ! On output they contain actual values, constrained by the size of the matrix.
  subroutine f_d_bw_expand_br(b_br, n, lbw_br, ubw_br, lbwmax_br, ubwmax_br, &
       dlbw, dlbw_tmp, dubw, dubw_tmp, lbw_bc, ubw_bc)
    integer(kind=int32), intent(in) :: n, lbwmax_br, ubwmax_br
    integer(kind=int32), intent(inout) :: lbw_br, ubw_br, dubw, dlbw, dubw_tmp, dlbw_tmp
    real(kind=dp), dimension(n,lbwmax_br+ubwmax_br+1), intent(out) :: b_br
    integer(kind=int32), intent(out) :: lbw_bc, ubw_bc

    integer(kind=int32) :: j

    ! add extra subdiagonals and shift as needed
    if (lbw_br+dlbw >= n-1) then
       dlbw=n-1 - lbw_br
       dlbw_tmp=0
    else if (lbw_br+dlbw+dlbw_tmp >= n-1) then
       dlbw_tmp = n-1 - lbw_br - dlbw
    end if
    lbw_bc=lbw_br+dlbw
    lbw_br=lbw_br+dlbw+dlbw_tmp
    do j=1,dlbw+dlbw_tmp
       call shift2(b_br,0,1)
    end do

    ! add extra superdiagonals.
    if (ubw_br+dubw >= n-1) then
       dubw=n-1 - ubw_br
       dubw_tmp=0
    else if (ubw_br+dubw+dubw_tmp >= n-1) then
       dubw_tmp = n-1 - ubw_br - dubw
    end if
    if (dubw+dubw_tmp > 0) then
       b_br(:,lbw_br+ubw_br+2:lbw_br+ubw_br+1+dubw + dubw_tmp)=0.0_dp
    end if
    ubw_bc=ubw_br+dubw
    ubw_br=ubw_br+dubw+dubw_tmp

  end subroutine f_d_bw_expand_br

  subroutine f_c_bw_expand_br(b_br, n, lbw_br, ubw_br, lbwmax_br, ubwmax_br, &
       dlbw, dlbw_tmp, dubw, dubw_tmp, lbw_bc, ubw_bc)
    integer(kind=int32), intent(in) :: n, lbwmax_br, ubwmax_br
    integer(kind=int32), intent(inout) :: lbw_br, ubw_br, dubw, dlbw, dubw_tmp, dlbw_tmp
    complex(kind=dp), dimension(n,lbwmax_br+ubwmax_br+1), intent(out) :: b_br
    integer(kind=int32), intent(out) :: lbw_bc, ubw_bc

    integer(kind=int32) :: j

    ! add extra subdiagonals and shift as needed
    if (lbw_br+dlbw >= n-1) then
       dlbw=n-1 - lbw_br
       dlbw_tmp=0
    else if (lbw_br+dlbw+dlbw_tmp >= n-1) then
       dlbw_tmp = n-1 - lbw_br - dlbw
    end if
    lbw_bc=lbw_br+dlbw
    lbw_br=lbw_br+dlbw+dlbw_tmp
    do j=1,dlbw+dlbw_tmp
       call shift2(b_br,0,1)
    end do

    ! add extra superdiagonals.
    if (ubw_br+dubw >= n-1) then
       dubw=n-1 - ubw_br
       dubw_tmp=0
    else if (ubw_br+dubw+dubw_tmp >= n-1) then
       dubw_tmp = n-1 - ubw_br - dubw
    end if
    if (dubw+dubw_tmp > 0) then
       b_br(:,lbw_br+ubw_br+2:lbw_br+ubw_br+1+dubw + dubw_tmp)=(0.0_dp,0.0_dp)
    end if
    ubw_bc=ubw_br+dubw
    ubw_br=ubw_br+dubw+dubw_tmp

  end subroutine f_c_bw_expand_br

  subroutine f_d_bw_contract_br(b_br, n, lbw_br, ubw_br, lbwmax_br, ubwmax_br, &
       dlbw, dlbw_tmp, dubw, dubw_tmp)
    integer(kind=int32), intent(in) :: n, lbwmax_br, ubwmax_br
    integer(kind=int32), intent(inout) :: lbw_br, ubw_br, dubw, dlbw, dubw_tmp, dlbw_tmp
    real(kind=dp), dimension(n,lbwmax_br+ubwmax_br+1), intent(out) :: b_br

    integer(kind=int32) :: j

    lbw_br=lbw_br - dlbw_tmp
    ubw_br=ubw_br - dubw_tmp
    do j=1,dlbw_tmp
       call shift2(b_br,0,-1)
    end do

  end subroutine f_d_bw_contract_br

  subroutine f_c_bw_contract_br(b_br, n, lbw_br, ubw_br, lbwmax_br, ubwmax_br, &
       dlbw, dlbw_tmp, dubw, dubw_tmp)
    integer(kind=int32), intent(in) :: n, lbwmax_br, ubwmax_br
    integer(kind=int32), intent(inout) :: lbw_br, ubw_br, dubw, dlbw, dubw_tmp, dlbw_tmp
    complex(kind=dp), dimension(n,lbwmax_br+ubwmax_br+1), intent(out) :: b_br

    integer(kind=int32) :: j

    lbw_br=lbw_br - dlbw_tmp
    ubw_br=ubw_br - dubw_tmp
    do j=1,dlbw_tmp
       call shift2(b_br,0,-1)
    end do

  end subroutine f_c_bw_contract_br

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

  subroutine c_extract_diagonals_bc(a, n, b, lbw, ubw, lbwmax, ubwmax)
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
  end subroutine c_extract_diagonals_bc

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

subroutine c_extract_diagonals_br(a, n, b, lbw, ubw, lbwmax, ubwmax)
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
  end subroutine c_extract_diagonals_br

end module mod_band_types
