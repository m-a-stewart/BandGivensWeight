module triangular
  use prec
  implicit none

  private
  public :: lower_right_invert, f_c_lower_right_invert, f_d_lower_right_invert, &
       lower_tr_right_invert, f_c_lower_tr_right_invert, f_d_lower_tr_right_invert, &
       lower_left_invert, f_c_lower_left_invert, f_d_lower_left_invert, &
       lower_tr_left_invert, f_c_lower_tr_left_invert, f_d_lower_tr_left_invert, &
       upper_right_invert, f_c_upper_right_invert, f_d_upper_right_invert

  public :: first_zero_diagonal, c_first_zero_diagonal, d_first_zero_diagonal, &
       last_zero_diagonal, c_last_zero_diagonal, d_last_zero_diagonal

  interface lower_right_invert
     module procedure f_c_lower_right_invert, f_d_lower_right_invert
  end interface lower_right_invert

  interface lower_tr_right_invert
     module procedure f_c_lower_tr_right_invert, f_d_lower_tr_right_invert
  end interface lower_tr_right_invert

  interface lower_left_invert
     module procedure f_c_lower_left_invert, f_d_lower_left_invert
  end interface lower_left_invert

  interface lower_tr_left_invert
     module procedure f_c_lower_tr_left_invert, f_d_lower_tr_left_invert
  end interface lower_tr_left_invert

  interface upper_right_invert
     module procedure f_c_upper_right_invert, f_d_upper_right_invert
  end interface upper_right_invert

  interface first_zero_diagonal
     module procedure c_first_zero_diagonal, d_first_zero_diagonal
  end interface first_zero_diagonal

  interface last_zero_diagonal
     module procedure c_last_zero_diagonal, d_last_zero_diagonal
  end interface last_zero_diagonal

contains


  ! Lower triangular

  ! x^T L = b^T
  subroutine f_c_lower_right_invert(x,l,b)
    complex(kind=dp), dimension(:,:), intent(in) :: l
    complex(kind=dp), dimension(:), intent(in) :: b
    complex(kind=dp), dimension(:), intent(out) :: x
    complex(kind=dp) :: tmp
    !
    integer(kind=int32) :: m,k,j
    !
    m=size(l,1)
    x(m)=b(m)/l(m,m)
    do k=m-1,1,-1
       tmp = b(k)
       do j=k+1,m
          tmp = tmp - x(j)*l(j,k)
       end do
       x(k) = tmp/l(k,k)
    end do
  end subroutine f_c_lower_right_invert

  ! x^T L = b^T
  subroutine f_d_lower_right_invert(x,l,b)
    real(kind=dp), dimension(:,:), intent(in) :: l
    real(kind=dp), dimension(:), intent(in) :: b
    real(kind=dp), dimension(:), intent(out) :: x
    real(kind=dp) :: tmp
    !
    integer(kind=int32) :: m,k,j
    !
    m=size(l,1)
    x(m)=b(m)/l(m,m)
    do k=m-1,1,-1
       tmp = b(k)
       do j=k+1,m
          tmp = tmp - x(j)*l(j,k)
       end do
       x(k) = tmp/l(k,k)
    end do
  end subroutine f_d_lower_right_invert

  ! x^T L^T = b^T
  subroutine f_d_lower_tr_right_invert(x,l,b)
    real(kind=dp), dimension(:,:), intent(in) :: l
    real(kind=dp), dimension(:), intent(in) :: b
    real(kind=dp), dimension(:), intent(out) :: x
    real(kind=dp) :: tmp
    !
    integer(kind=int32) :: m,k,j
    !
    m=size(l,1)
    x(1)=b(1)/l(1,1)
    do k=2,m
       tmp = b(k)
       do j=1,k-1
          tmp = tmp - x(j)*l(k,j)
       end do
       x(k) = tmp/l(k,k)
    end do
  end subroutine f_d_lower_tr_right_invert

  ! x^T L^H = b^T
  subroutine f_c_lower_tr_right_invert(x,l,b)
    complex(kind=dp), dimension(:,:), intent(in) :: l
    complex(kind=dp), dimension(:), intent(in) :: b
    complex(kind=dp), dimension(:), intent(out) :: x
    complex(kind=dp) :: tmp
    !
    integer(kind=int32) :: m,k,j
    !
    m=size(l,1)
    x(1)=b(1)/conjg(l(1,1))
    do k=2,m
       tmp = b(k)
       do j=1,k-1
          tmp = tmp - x(j)*conjg(l(k,j))
       end do
       x(k) = tmp/conjg(l(k,k))
    end do
  end subroutine f_c_lower_tr_right_invert

  !
  ! applying inverses from the left.
  !

  ! L x = b
  subroutine f_c_lower_left_invert(l,x,b)
    complex(kind=dp), dimension(:,:), intent(in) :: l
    complex(kind=dp), dimension(:), intent(in) :: b
    complex(kind=dp), dimension(:), intent(out) :: x
    complex(kind=dp) :: tmp
    !
    integer(kind=int32) :: m,k,j
    !
    m=size(l,1)
    x(1)=b(1)/l(1,1)
    do k=2,m
       tmp=b(k)
       do j=1,k-1
          tmp = tmp - l(k,j)*x(j)
       end do
       x(k)=tmp/l(k,k)
    end do
  end subroutine f_c_lower_left_invert

  ! L x = b
  subroutine f_d_lower_left_invert(l,x,b)
    real(kind=dp), dimension(:,:), intent(in) :: l
    real(kind=dp), dimension(:), intent(in) :: b
    real(kind=dp), dimension(:), intent(out) :: x
    real(kind=dp) :: tmp
    !
    integer(kind=int32) :: m,k,j
    !
    m=size(l,1)
    x(1)=b(1)/l(1,1)
    do k=2,m
       tmp=b(k)
       do j=1,k-1
          tmp = tmp - l(k,j)*x(j)
       end do
       x(k)=tmp/l(k,k)
    end do
  end subroutine f_d_lower_left_invert

  ! L^T x = b
  subroutine f_d_lower_tr_left_invert(l,x,b)
    real(kind=dp), dimension(:,:), intent(in) :: l
    real(kind=dp), dimension(:), intent(in) :: b
    real(kind=dp), dimension(:), intent(out) :: x
    real(kind=dp) :: tmp
    !
    integer(kind=int32) :: m,k,j
    !
    m=size(l,1)
    x(m)=b(m)/l(m,m)
    do k=m-1,1,-1
       tmp=b(k)
       do j=k+1,m
          tmp=tmp-l(j,k)*x(j)
       end do
       x(k)=tmp/l(k,k)
    end do
  end subroutine f_d_lower_tr_left_invert

  ! L^H x = b
  subroutine f_c_lower_tr_left_invert(l,x,b)
    complex(kind=dp), dimension(:,:), intent(in) :: l
    complex(kind=dp), dimension(:), intent(in) :: b
    complex(kind=dp), dimension(:), intent(out) :: x
    complex(kind=dp) :: tmp
    !
    integer(kind=int32) :: m,k,j
    !
    m=size(l,1)
    x(m)=b(m)/conjg(l(m,m))
    do k=m-1,1,-1
       tmp=b(k)
       do j=k+1,m
          tmp=tmp-conjg(l(j,k))*x(j)
       end do
       x(k)=tmp/conjg(l(k,k))
    end do
  end subroutine f_c_lower_tr_left_invert

  ! Upper triangular

  ! x^T R = b^T
  subroutine f_d_upper_right_invert(x,r,b)
    real(kind=dp), dimension(:,:), intent(in) :: r
    real(kind=dp), dimension(:), intent(in) :: b
    real(kind=dp), dimension(:), intent(out) :: x
    real(kind=dp) :: tmp
    !
    integer(kind=int32) :: m,k,j
    !
    m=size(r,1)
    x(1)=b(1)/r(1,1)
    do k=2,m
       tmp=b(k)
       do j=1,k-1
          tmp = tmp - x(j)*r(j,k)
       end do
       x(k)=tmp/r(k,k)
    end do
  end subroutine f_d_upper_right_invert

  ! x^T R = b^T
  subroutine f_c_upper_right_invert(x,r,b)
    complex(kind=dp), dimension(:,:), intent(in) :: r
    complex(kind=dp), dimension(:), intent(in) :: b
    complex(kind=dp), dimension(:), intent(out) :: x
    complex(kind=dp) :: tmp
    !
    integer(kind=int32) :: m,k,j
    !
    m=size(r,1)
    x(1)=b(1)/r(1,1)
    do k=2,m
       tmp=b(k)
       do j=1,k-1
          tmp = tmp - x(j)*r(j,k)
       end do
       x(k)=tmp/r(k,k)
    end do
  end subroutine f_c_upper_right_invert

  ! x^T R^T = b^T
  subroutine f_d_upper_tr_right_invert(x,r,b)
    real(kind=dp), dimension(:,:), intent(in) :: r
    real(kind=dp), dimension(:), intent(in) :: b
    real(kind=dp), dimension(:), intent(out) :: x
    real(kind=dp) :: tmp
    !
    integer(kind=int32) :: m,k,j
    !
    m=size(r,1)
    x(m)=b(m)/r(m,m)
    do k=m-1,1,-1
       tmp=b(k)
       do j=k+1,m
          tmp = tmp - x(j)*r(k,j)
       end do
       x(k) = tmp/r(k,k)
    end do
  end subroutine f_d_upper_tr_right_invert

! x^T R^H = b^T
  subroutine f_c_upper_tr_right_invert(x,r,b)
    complex(kind=dp), dimension(:,:), intent(in) :: r
    complex(kind=dp), dimension(:), intent(in) :: b
    complex(kind=dp), dimension(:), intent(out) :: x
    complex(kind=dp) :: tmp
    !
    integer(kind=int32) :: m,k,j
    !
    m=size(r,1)
    x(m)=b(m)/conjg(r(m,m))
    do k=m-1,1,-1
       tmp=b(k)
       do j=k+1,m
          tmp = tmp - x(j)*conjg(r(k,j))
       end do
       x(k) = tmp/conjg(r(k,k))
    end do
  end subroutine f_c_upper_tr_right_invert

  ! R x = b
  subroutine f_d_upper_left_invert(r,x,b)
    real(kind=dp), dimension(:,:), intent(in) :: r
    real(kind=dp), dimension(:), intent(in) :: b
    real(kind=dp), dimension(:), intent(out) :: x
    real(kind=dp) :: tmp
    !
    integer(kind=int32) :: m,k,j
    !
    m=size(r,1)
    x(m)=b(m)/r(m,m)
    do k=m-1,1,-1
       tmp=b(k)
       do j=k+1,m
          tmp = tmp - r(k,j)*x(j)
       end do
       x(k)=tmp/r(k,k)
    end do
  end subroutine f_d_upper_left_invert

  ! R x = b
  subroutine f_c_upper_left_invert(r,x,b)
    complex(kind=dp), dimension(:,:), intent(in) :: r
    complex(kind=dp), dimension(:), intent(in) :: b
    complex(kind=dp), dimension(:), intent(out) :: x
    complex(kind=dp) :: tmp
    !
    integer(kind=int32) :: m,k,j
    !
    m=size(r,1)
    x(m)=b(m)/r(m,m)
    do k=m-1,1,-1
       tmp=b(k)
       do j=k+1,m
          tmp = tmp - r(k,j)*x(j)
       end do
       x(k)=tmp/r(k,k)
    end do
  end subroutine f_c_upper_left_invert

  ! R^T x = b
  subroutine f_d_upper_tr_left_invert(r,x,b)
    real(kind=dp), dimension(:,:), intent(in) :: r
    real(kind=dp), dimension(:), intent(in) :: b
    real(kind=dp), dimension(:), intent(out) :: x
    real(kind=dp) :: tmp
    !
    integer(kind=int32) :: m,k,j
    !
    m=size(r,1)
    x(1)=b(1)/r(1,1)
    do k=2,m
       tmp=b(k)
       do j=1,k-1
          tmp=tmp-r(j,k)*x(j)
       end do
       x(k)=tmp/r(k,k)
    end do
  end subroutine f_d_upper_tr_left_invert

  ! R^H x = b
  subroutine f_c_upper_tr_left_invert(r,x,b)
    complex(kind=dp), dimension(:,:), intent(in) :: r
    complex(kind=dp), dimension(:), intent(in) :: b
    complex(kind=dp), dimension(:), intent(out) :: x
    complex(kind=dp) :: tmp
    !
    integer(kind=int32) :: m,k,j
    !
    m=size(r,1)
    x(1)=b(1)/conjg(r(1,1))
    do k=2,m
       tmp=b(k)
       do j=1,k-1
          tmp=tmp-conjg(r(j,k))*x(j)
       end do
       x(k)=tmp/conjg(r(k,k))
    end do
  end subroutine f_c_upper_tr_left_invert

  ! finding zeros on the diagonal

  integer(kind=int32) function d_first_zero_diagonal(a,tol) result(d)
    real(kind=dp), dimension(:,:), intent(in) :: a
    real(kind=dp), intent(in) :: tol
    !
    integer(kind=int32) :: j, n
    n=min(size(a,1), size(a,2))
    d=1
    do j=1,n
       if (abs(a(j,j)) <= tol) then
          return
       end if
       d=d+1
    end do
  end function d_first_zero_diagonal

  integer(kind=int32) function c_first_zero_diagonal(a,tol) result(d)
    complex(kind=dp), dimension(:,:), intent(in) :: a
    real(kind=dp), intent(in) :: tol
    !
    integer(kind=int32) :: j, n
    n=min(size(a,1), size(a,2))
    d=1
    do j=1,n
       if (abs(a(j,j)) <= tol) then
          return
       end if
       d=d+1
    end do
  end function c_first_zero_diagonal

  integer(kind=int32) function d_last_zero_diagonal(a,tol) result(d)
    real(kind=dp), dimension(:,:), intent(in) :: a
    real(kind=dp), intent(in) :: tol
    !
    integer(kind=int32) :: j, n
    n=min(size(a,1), size(a,2))
    d=n
    do j=n,1,-1
       if (abs(a(j,j)) <= tol) then
          return
       end if
       d=d-1
    end do
  end function d_last_zero_diagonal

  integer(kind=int32) function c_last_zero_diagonal(a,tol) result(d)
    complex(kind=dp), dimension(:,:), intent(in) :: a
    real(kind=dp), intent(in) :: tol
    !
    integer(kind=int32) :: j, n
    n=min(size(a,1), size(a,2))
    d=n
    do j=n,1,-1
       if (abs(a(j,j)) <= tol) then
          return
       end if
       d=d-1
    end do
  end function c_last_zero_diagonal

end module triangular
