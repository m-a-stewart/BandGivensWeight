module mod_triangular
  use mod_prec
  implicit none

  private

  public :: lower_right_multiply, f_z_lower_right_multiply, f_d_lower_right_multiply, &
       lower_tr_right_multiply, f_z_lower_tr_right_multiply, f_d_lower_tr_right_multiply, &
       lower_left_multiply, f_z_lower_left_multiply, f_d_lower_left_multiply, &
       lower_tr_left_multiply, f_z_lower_tr_left_multiply, f_d_lower_tr_left_multiply, &
       upper_right_multiply, f_z_upper_right_multiply, f_d_upper_right_multiply, &
       upper_tr_right_multiply, f_z_upper_tr_right_multiply, f_d_upper_tr_right_multiply, &
       upper_left_multiply, f_z_upper_left_multiply, f_d_upper_left_multiply, &
       upper_tr_left_multiply, f_z_upper_tr_left_multiply, f_d_upper_tr_left_multiply

  public :: lower_right_multiply_linpack, f_z_lower_right_multiply_linpack, &
       f_d_lower_right_multiply_linpack, lower_tr_right_multiply_linpack, &
       f_z_lower_tr_right_multiply_linpack, f_d_lower_tr_right_multiply_linpack, &
       lower_left_multiply_linpack, f_z_lower_left_multiply_linpack, &
       f_d_lower_left_multiply_linpack, lower_tr_left_multiply_linpack, &
       f_z_lower_tr_left_multiply_linpack, f_d_lower_tr_left_multiply_linpack, &
       upper_right_multiply_linpack, f_z_upper_right_multiply_linpack, &
       f_d_upper_right_multiply_linpack, upper_tr_right_multiply_linpack, &
       f_z_upper_tr_right_multiply_linpack, f_d_upper_tr_right_multiply_linpack, &
       upper_left_multiply_linpack, f_z_upper_left_multiply_linpack, f_d_upper_left_multiply_linpack, &
       upper_tr_left_multiply_linpack, f_z_upper_tr_left_multiply_linpack, &
       f_d_upper_tr_left_multiply_linpack


  public :: lower_right_invert, f_z_lower_right_invert, f_d_lower_right_invert, &
       lower_tr_right_invert, f_z_lower_tr_right_invert, f_d_lower_tr_right_invert, &
       lower_left_invert, f_z_lower_left_invert, f_d_lower_left_invert, &
       lower_tr_left_invert, f_z_lower_tr_left_invert, f_d_lower_tr_left_invert, &
       upper_right_invert, f_z_upper_right_invert, f_d_upper_right_invert, &
       upper_tr_right_invert, f_z_upper_tr_right_invert, f_d_upper_tr_right_invert, &
       upper_left_invert, f_z_upper_left_invert, f_d_upper_left_invert, &
       upper_tr_left_invert, f_z_upper_tr_left_invert, f_d_upper_tr_left_invert
  
  public :: lower_right_invert_linpack, f_z_lower_right_invert_linpack, &
       f_d_lower_right_invert_linpack, lower_tr_right_invert_linpack, &
       f_z_lower_tr_right_invert_linpack, f_d_lower_tr_right_invert_linpack, &
       lower_left_invert_linpack, f_z_lower_left_invert_linpack, &
       f_d_lower_left_invert_linpack, lower_tr_left_invert_linpack, &
       f_z_lower_tr_left_invert_linpack, f_d_lower_tr_left_invert_linpack, &
       upper_right_invert_linpack, f_z_upper_right_invert_linpack, &
       f_d_upper_right_invert_linpack, upper_tr_right_invert_linpack, &
       f_z_upper_tr_right_invert_linpack, f_d_upper_tr_right_invert_linpack, &
       upper_left_invert_linpack, f_z_upper_left_invert_linpack, f_d_upper_left_invert_linpack, &
       upper_tr_left_invert_linpack, f_z_upper_tr_left_invert_linpack, &
       f_d_upper_tr_left_invert_linpack

  public :: first_zero_diagonal, z_first_zero_diagonal, d_first_zero_diagonal, &
       reverse_first_zero_diagonal, z_reverse_first_zero_diagonal, d_reverse_first_zero_diagonal

  ! Multiplication
  
  interface lower_right_multiply
     module procedure f_z_lower_right_multiply, f_d_lower_right_multiply
  end interface lower_right_multiply

  interface lower_tr_right_multiply
     module procedure f_z_lower_tr_right_multiply, f_d_lower_tr_right_multiply
  end interface lower_tr_right_multiply

  interface lower_left_multiply
     module procedure f_z_lower_left_multiply, f_d_lower_left_multiply
  end interface lower_left_multiply

  interface lower_tr_left_multiply
     module procedure f_z_lower_tr_left_multiply, f_d_lower_tr_left_multiply
  end interface lower_tr_left_multiply

  interface upper_right_multiply
     module procedure f_z_upper_right_multiply, f_d_upper_right_multiply
  end interface upper_right_multiply

  interface upper_tr_right_multiply
     module procedure f_z_upper_tr_right_multiply, f_d_upper_tr_right_multiply
  end interface upper_tr_right_multiply

  interface upper_left_multiply
     module procedure f_z_upper_left_multiply, f_d_upper_left_multiply
  end interface upper_left_multiply

  interface upper_tr_left_multiply
     module procedure f_z_upper_tr_left_multiply, f_d_upper_tr_left_multiply
  end interface upper_tr_left_multiply

  ! Multiplication, Linpack
  
  interface lower_right_multiply_linpack
     module procedure f_z_lower_right_multiply_linpack, f_d_lower_right_multiply_linpack
  end interface lower_right_multiply_linpack

  interface lower_tr_right_multiply_linpack
     module procedure f_z_lower_tr_right_multiply_linpack, f_d_lower_tr_right_multiply_linpack
  end interface lower_tr_right_multiply_linpack

  interface lower_left_multiply_linpack
     module procedure f_z_lower_left_multiply_linpack, f_d_lower_left_multiply_linpack
  end interface lower_left_multiply_linpack

  interface lower_tr_left_multiply_linpack
     module procedure f_z_lower_tr_left_multiply_linpack, f_d_lower_tr_left_multiply_linpack
  end interface lower_tr_left_multiply_linpack

  interface upper_right_multiply_linpack
     module procedure f_z_upper_right_multiply_linpack, f_d_upper_right_multiply_linpack
  end interface upper_right_multiply_linpack

  interface upper_tr_right_multiply_linpack
     module procedure f_z_upper_tr_right_multiply_linpack, f_d_upper_tr_right_multiply_linpack
  end interface upper_tr_right_multiply_linpack

  interface upper_left_multiply_linpack
     module procedure f_z_upper_left_multiply_linpack, f_d_upper_left_multiply_linpack
  end interface upper_left_multiply_linpack

  interface upper_tr_left_multiply_linpack
     module procedure f_z_upper_tr_left_multiply_linpack, f_d_upper_tr_left_multiply_linpack
  end interface upper_tr_left_multiply_linpack
  
  ! Inversion
  
  interface lower_right_invert
     module procedure f_z_lower_right_invert, f_d_lower_right_invert
  end interface lower_right_invert

  interface lower_tr_right_invert
     module procedure f_z_lower_tr_right_invert, f_d_lower_tr_right_invert
  end interface lower_tr_right_invert

  interface lower_left_invert
     module procedure f_z_lower_left_invert, f_d_lower_left_invert
  end interface lower_left_invert

  interface lower_tr_left_invert
     module procedure f_z_lower_tr_left_invert, f_d_lower_tr_left_invert
  end interface lower_tr_left_invert

  interface upper_right_invert
     module procedure f_z_upper_right_invert, f_d_upper_right_invert
  end interface upper_right_invert

  interface upper_tr_right_invert
     module procedure f_z_upper_tr_right_invert, f_d_upper_tr_right_invert
  end interface upper_tr_right_invert

  interface upper_left_invert
     module procedure f_z_upper_left_invert, f_d_upper_left_invert
  end interface upper_left_invert

  interface upper_tr_left_invert
     module procedure f_z_upper_tr_left_invert, f_d_upper_tr_left_invert
  end interface upper_tr_left_invert

  ! Linpack inversion.
  
  interface lower_right_invert_linpack
     module procedure f_z_lower_right_invert_linpack, f_d_lower_right_invert_linpack
  end interface lower_right_invert_linpack

  interface lower_tr_right_invert_linpack
     module procedure f_z_lower_tr_right_invert_linpack, f_d_lower_tr_right_invert_linpack
  end interface lower_tr_right_invert_linpack

  interface lower_left_invert_linpack
     module procedure f_z_lower_left_invert_linpack, f_d_lower_left_invert_linpack
  end interface lower_left_invert_linpack

  interface lower_tr_left_invert_linpack
     module procedure f_z_lower_tr_left_invert_linpack, f_d_lower_tr_left_invert_linpack
  end interface lower_tr_left_invert_linpack

  interface upper_right_invert_linpack
     module procedure f_z_upper_right_invert_linpack, f_d_upper_right_invert_linpack
  end interface upper_right_invert_linpack

  interface upper_tr_right_invert_linpack
     module procedure f_z_upper_tr_right_invert_linpack, f_d_upper_tr_right_invert_linpack
  end interface upper_tr_right_invert_linpack

  interface upper_left_invert_linpack
     module procedure f_z_upper_left_invert_linpack, f_d_upper_left_invert_linpack
  end interface upper_left_invert_linpack

  interface upper_tr_left_invert_linpack
     module procedure f_z_upper_tr_left_invert_linpack, f_d_upper_tr_left_invert_linpack
  end interface upper_tr_left_invert_linpack
  
  
  interface first_zero_diagonal
     module procedure z_first_zero_diagonal, d_first_zero_diagonal
  end interface first_zero_diagonal

  interface reverse_first_zero_diagonal
     module procedure z_reverse_first_zero_diagonal, d_reverse_first_zero_diagonal
  end interface reverse_first_zero_diagonal

contains

  ! Multiply
  
  ! x^T L = b^T
  subroutine f_d_lower_right_multiply(x,l,b)
    real(kind=dp), dimension(:,:), intent(in) :: l
    real(kind=dp), dimension(:), intent(out) :: b
    real(kind=dp), dimension(:), intent(in) :: x
    !
    integer(kind=int32) :: m,k,j
    !
    m=size(l,1)
    do k=1,m
       b(k)=0.0_dp
       do j=k,m
          b(k)=b(k)+x(j)*l(j,k)
       end do
    end do
  end subroutine f_d_lower_right_multiply

  ! x^T L = b^T
  subroutine f_z_lower_right_multiply(x,l,b)
    complex(kind=dp), dimension(:,:), intent(in) :: l
    complex(kind=dp), dimension(:), intent(out) :: b
    complex(kind=dp), dimension(:), intent(in) :: x
    !
    integer(kind=int32) :: m,k,j
    !
    m=size(l,1)
    do k=1,m
       b(k)=(0.0_dp,0.0_dp)
       do j=k,m
          b(k)=b(k)+x(j)*l(j,k)
       end do
    end do
  end subroutine f_z_lower_right_multiply

  ! x^T L^T = b^T
  subroutine f_d_lower_tr_right_multiply(x,l,b)
    real(kind=dp), dimension(:,:), intent(in) :: l
    real(kind=dp), dimension(:), intent(out) :: b
    real(kind=dp), dimension(:), intent(in) :: x
    !
    integer(kind=int32) :: m,k,j
    !
    m=size(l,1)
    do k=1,m
       b(k)=0.0_dp
       do j=1,k
          b(k)=b(k)+x(j)*l(k,j)
       end do
    end do
  end subroutine f_d_lower_tr_right_multiply

  ! x^T L^H = b^T
  subroutine f_z_lower_tr_right_multiply(x,l,b)
    complex(kind=dp), dimension(:,:), intent(in) :: l
    complex(kind=dp), dimension(:), intent(out) :: b
    complex(kind=dp), dimension(:), intent(in) :: x
    !
    integer(kind=int32) :: m,k,j
    !
    m=size(l,1)
    do k=1,m
       b(k)=(0.0_dp,0.0_dp)
       do j=1,k
          b(k)=b(k)+x(j)*conjg(l(k,j))
       end do
    end do
  end subroutine f_z_lower_tr_right_multiply
  
  ! L x = b
  subroutine f_d_lower_left_multiply(l,x,b)
    real(kind=dp), dimension(:,:), intent(in) :: l
    real(kind=dp), dimension(:), intent(out) :: b
    real(kind=dp), dimension(:), intent(in) :: x
    !
    integer(kind=int32) :: m,k,j
    !
    m=size(l,1)
    do j=1,m
       b(j)=0.0_dp
       do k=1,j
          b(j)=b(j)+l(j,k)*x(k)
       end do
    end do
  end subroutine f_d_lower_left_multiply

  ! L x = b
  subroutine f_z_lower_left_multiply(l,x,b)
    complex(kind=dp), dimension(:,:), intent(in) :: l
    complex(kind=dp), dimension(:), intent(out) :: b
    complex(kind=dp), dimension(:), intent(in) :: x
    !
    integer(kind=int32) :: m,k,j
    !
    m=size(l,1)
    do j=1,m
       b(j)=(0.0_dp,0.0_dp)
       do k=1,j
          b(j)=b(j)+l(j,k)*x(k)
       end do
    end do
  end subroutine f_z_lower_left_multiply

  ! L^T x = b
  subroutine f_d_lower_tr_left_multiply(l,x,b)
    real(kind=dp), dimension(:,:), intent(in) :: l
    real(kind=dp), dimension(:), intent(out) :: b
    real(kind=dp), dimension(:), intent(in) :: x
    !
    integer(kind=int32) :: m,k,j
    !
    m=size(l,1)
    do j=1,m
       b(j)=0.0_dp
       do k=j,m
          b(j)=b(j)+l(k,j)*x(k)
       end do
    end do
  end subroutine f_d_lower_tr_left_multiply

  ! L^H x = b
  subroutine f_z_lower_tr_left_multiply(l,x,b)
    complex(kind=dp), dimension(:,:), intent(in) :: l
    complex(kind=dp), dimension(:), intent(out) :: b
    complex(kind=dp), dimension(:), intent(in) :: x
    !
    integer(kind=int32) :: m,k,j
    !
    m=size(l,1)
    do j=1,m
       b(j)=(0.0_dp,0.0_dp)
       do k=j,m
          b(j)=b(j)+conjg(l(k,j))*x(k)
       end do
    end do
  end subroutine f_z_lower_tr_left_multiply

  ! x^T U = b^T
  subroutine f_d_upper_right_multiply(x,u,b)
    real(kind=dp), dimension(:,:), intent(in) :: u
    real(kind=dp), dimension(:), intent(out) :: b
    real(kind=dp), dimension(:), intent(in) :: x
    !
    integer(kind=int32) :: m,k,j
    !
    m=size(u,1)
    do k=1,m
       b(k)=0.0_dp
       do j=1,k
          b(k)=b(k)+x(j)*u(j,k)
       end do
    end do
  end subroutine f_d_upper_right_multiply

  ! x^T U = b^T
  subroutine f_z_upper_right_multiply(x,u,b)
    complex(kind=dp), dimension(:,:), intent(in) :: u
    complex(kind=dp), dimension(:), intent(out) :: b
    complex(kind=dp), dimension(:), intent(in) :: x
    !
    integer(kind=int32) :: m,k,j
    !
    m=size(u,1)
    do k=1,m
       b(k)=(0.0_dp,0.0_dp)
       do j=1,k
          b(k)=b(k)+x(j)*u(j,k)
       end do
    end do
  end subroutine f_z_upper_right_multiply

  ! x^T U^T = b^T
  subroutine f_d_upper_tr_right_multiply(x,u,b)
    real(kind=dp), dimension(:,:), intent(in) :: u
    real(kind=dp), dimension(:), intent(out) :: b
    real(kind=dp), dimension(:), intent(in) :: x
    !
    integer(kind=int32) :: m,k,j
    !
    m=size(u,1)
    do k=1,m
       b(k)=0.0_dp
       do j=k,m
          b(k)=b(k)+x(j)*u(k,j)
       end do
    end do
  end subroutine f_d_upper_tr_right_multiply

  ! x^T U^H = b^T
  subroutine f_z_upper_tr_right_multiply(x,u,b)
    complex(kind=dp), dimension(:,:), intent(in) :: u
    complex(kind=dp), dimension(:), intent(out) :: b
    complex(kind=dp), dimension(:), intent(in) :: x
    !
    integer(kind=int32) :: m,k,j
    !
    m=size(u,1)
    do k=1,m
       b(k)=(0.0_dp,0.0_dp)
       do j=k,m
          b(k)=b(k)+x(j)*conjg(u(k,j))
       end do
    end do
  end subroutine f_z_upper_tr_right_multiply
  
  ! U x = b
  subroutine f_d_upper_left_multiply(u,x,b)
    real(kind=dp), dimension(:,:), intent(in) :: u
    real(kind=dp), dimension(:), intent(out) :: b
    real(kind=dp), dimension(:), intent(in) :: x
    !
    integer(kind=int32) :: m,k,j
    !
    m=size(u,1)
    do j=1,m
       b(j)=0.0_dp
       do k=j,m
          b(j)=b(j)+u(j,k)*x(k)
       end do
    end do
  end subroutine f_d_upper_left_multiply

  ! U x = b
  subroutine f_z_upper_left_multiply(u,x,b)
    complex(kind=dp), dimension(:,:), intent(in) :: u
    complex(kind=dp), dimension(:), intent(out) :: b
    complex(kind=dp), dimension(:), intent(in) :: x
    !
    integer(kind=int32) :: m,k,j
    !
    m=size(u,1)
    do j=1,m
       b(j)=(0.0_dp,0.0_dp)
       do k=j,m
          b(j)=b(j)+u(j,k)*x(k)
       end do
    end do
  end subroutine f_z_upper_left_multiply

  ! U^T x = b
  subroutine f_d_upper_tr_left_multiply(u,x,b)
    real(kind=dp), dimension(:,:), intent(in) :: u
    real(kind=dp), dimension(:), intent(out) :: b
    real(kind=dp), dimension(:), intent(in) :: x
    !
    integer(kind=int32) :: m,k,j
    !
    m=size(u,1)
    do j=1,m
       b(j)=0.0_dp
       do k=1,j
          b(j)=b(j)+u(k,j)*x(k)
       end do
    end do
  end subroutine f_d_upper_tr_left_multiply

  ! U^* x = b
  subroutine f_z_upper_tr_left_multiply(u,x,b)
    complex(kind=dp), dimension(:,:), intent(in) :: u
    complex(kind=dp), dimension(:), intent(out) :: b
    complex(kind=dp), dimension(:), intent(in) :: x
    !
    integer(kind=int32) :: m,k,j
    !
    m=size(u,1)
    do j=1,m
       b(j)=(0.0_dp,0.0_dp)
       do k=1,j
          b(j)=b(j)+conjg(u(k,j))*x(k)
       end do
    end do
  end subroutine f_z_upper_tr_left_multiply

  ! Multiply Linpack

  ! x^T L = b^T
  subroutine f_d_lower_right_multiply_linpack(x,l,b)
    real(kind=dp), dimension(:,:), intent(in) :: l
    real(kind=dp), dimension(:), intent(out) :: b
    real(kind=dp), dimension(:), intent(out) :: x
    !
    integer(kind=int32) :: m,k,j
    !
    m=size(l,1)
    x=1.0_dp
    b=0.0_dp
    do k=m,1,-1
       do j=k+1,m
          b(k)=b(k)+x(j)*l(j,k)
       end do
       if (b(k) /= 0.0_dp .and. l(k,k) /= 0.0_dp) then
          x(k)=b(k)/l(k,k); x(k)=x(k)/abs(x(k))
       end if
       b(k)=b(k)+x(k)*l(k,k)
    end do
  end subroutine f_d_lower_right_multiply_linpack

  ! x^T L = b^T
  subroutine f_z_lower_right_multiply_linpack(x,l,b)
    complex(kind=dp), dimension(:,:), intent(in) :: l
    complex(kind=dp), dimension(:), intent(out) :: b
    complex(kind=dp), dimension(:), intent(out) :: x
    !
    integer(kind=int32) :: m,k,j
    !
    m=size(l,1)
    x=(1.0_dp,0.0_dp)
    b=(0.0_dp,0.0_dp)
    do k=m,1,-1
       do j=k+1,m
          b(k)=b(k)+x(j)*l(j,k)
       end do
       if (b(k) /= (0.0_dp,0.0_dp) .and. l(k,k) /= (0.0_dp,0.0_dp)) then
          x(k)=b(k)/l(k,k); x(k)=x(k)/abs(x(k))
       end if
       b(k)=b(k)+x(k)*l(k,k)
    end do
  end subroutine f_z_lower_right_multiply_linpack

  ! x^T L^T = b^T
  subroutine f_d_lower_tr_right_multiply_linpack(x,l,b)
    real(kind=dp), dimension(:,:), intent(in) :: l
    real(kind=dp), dimension(:), intent(out) :: b
    real(kind=dp), dimension(:), intent(out) :: x
    !
    integer(kind=int32) :: m,k,j
    !
    m=size(l,1)
    x=1.0_dp
    b=0.0_dp
    do k=1,m
       do j=1,k-1
          b(k)=b(k)+x(j)*l(k,j)
       end do
       if (b(k) /= 0.0_dp .and. l(k,k) /= 0.0_dp) then
          x(k)=b(k)/l(k,k); x(k)=x(k)/abs(x(k))
       end if
       b(k)=b(k)+x(k)*l(k,k)
    end do
  end subroutine f_d_lower_tr_right_multiply_linpack

  ! x^T L^T = b^T
  subroutine f_z_lower_tr_right_multiply_linpack(x,l,b)
    complex(kind=dp), dimension(:,:), intent(in) :: l
    complex(kind=dp), dimension(:), intent(out) :: b
    complex(kind=dp), dimension(:), intent(out) :: x
    !
    integer(kind=int32) :: m,k,j
    !
    m=size(l,1)
    x=(1.0_dp,0.0_dp)
    b=(0.0_dp,0.0_dp)
    do k=1,m
       do j=1,k-1
          b(k)=b(k)+x(j)*conjg(l(k,j))
       end do
       if (b(k) /= (0.0_dp,0.0_dp) .and. l(k,k) /= (0.0_dp,0.0_dp)) then
          x(k)=b(k)/conjg(l(k,k)); x(k)=x(k)/abs(x(k))
       end if
       b(k)=b(k)+x(k)*conjg(l(k,k))
    end do
  end subroutine f_z_lower_tr_right_multiply_linpack

  ! L x = b
  subroutine f_d_lower_left_multiply_linpack(l,x,b)
    real(kind=dp), dimension(:,:), intent(in) :: l
    real(kind=dp), dimension(:), intent(out) :: b
    real(kind=dp), dimension(:), intent(out) :: x
    !
    integer(kind=int32) :: m,k,j
    !
    m=size(l,1)
    x=1.0_dp
    b=0.0_dp
    do j=1,m
       do k=1,j-1
          b(j)=b(j)+l(j,k)*x(k)
       end do
       if (b(j) /= 0.0_dp .and. l(j,j) /= 0.0_dp) then
          x(j)=b(j)/l(j,j); x(j)=x(j)/abs(x(j))
       end if
       b(j)=b(j)+l(j,j)*x(j)
    end do
  end subroutine f_d_lower_left_multiply_linpack

  ! L x = b
  subroutine f_z_lower_left_multiply_linpack(l,x,b)
    complex(kind=dp), dimension(:,:), intent(in) :: l
    complex(kind=dp), dimension(:), intent(out) :: b
    complex(kind=dp), dimension(:), intent(out) :: x
    !
    integer(kind=int32) :: m,k,j
    !
    m=size(l,1)
    x=(1.0_dp,0.0_dp)
    b=(0.0_dp,0.0_dp)
    do j=1,m
       do k=1,j-1
          b(j)=b(j)+l(j,k)*x(k)
       end do
       if (b(j) /= (0.0_dp,0.0_dp) .and. l(j,j) /= (0.0_dp,0.0_dp)) then
          x(j)=b(j)/l(j,j); x(j)=x(j)/abs(x(j))
       end if
       b(j)=b(j)+l(j,j)*x(j)
    end do
  end subroutine f_z_lower_left_multiply_linpack

  ! L^T x = b
  subroutine f_d_lower_tr_left_multiply_linpack(l,x,b)
    real(kind=dp), dimension(:,:), intent(in) :: l
    real(kind=dp), dimension(:), intent(out) :: b
    real(kind=dp), dimension(:), intent(out) :: x
    !
    integer(kind=int32) :: m,k,j
    !
    m=size(l,1)
    x=1.0_dp
    b=0.0_dp
    do j=m,1,-1
       do k=j+1,m
          b(j)=b(j)+l(k,j)*x(k)
       end do
       if (b(j) /= 0.0_dp .and. l(j,j) /= 0.0_dp) then
          x(j)=b(j)/l(j,j); x(j)=x(j)/abs(x(j))
       end if
       b(j)=b(j)+l(j,j)*x(j)
    end do
  end subroutine f_d_lower_tr_left_multiply_linpack

  ! L^H x = b
  subroutine f_z_lower_tr_left_multiply_linpack(l,x,b)
    complex(kind=dp), dimension(:,:), intent(in) :: l
    complex(kind=dp), dimension(:), intent(out) :: b
    complex(kind=dp), dimension(:), intent(out) :: x
    !
    integer(kind=int32) :: m,k,j
    !
    m=size(l,1)
    x=(1.0_dp,0.0_dp)
    b=(0.0_dp,0.0_dp)
    do j=m,1,-1
       do k=j+1,m
          b(j)=b(j)+conjg(l(k,j))*x(k)
       end do
       if (b(j) /= (0.0_dp,0.0_dp) .and. l(j,j) /= (0.0_dp,0.0_dp)) then
          x(j)=b(j)/conjg(l(j,j)); x(j)=x(j)/abs(x(j))
       end if
       b(j)=b(j)+conjg(l(j,j))*x(j)
    end do
  end subroutine f_z_lower_tr_left_multiply_linpack

  ! x^T U = b^T
  subroutine f_d_upper_right_multiply_linpack(x,u,b)
    real(kind=dp), dimension(:,:), intent(in) :: u
    real(kind=dp), dimension(:), intent(out) :: b
    real(kind=dp), dimension(:), intent(out) :: x
    !
    integer(kind=int32) :: m,k,j
    !
    m=size(u,1)
    x=1.0_dp
    b=0.0_dp
    do k=1,m
       do j=1,k-1
          b(k)=b(k)+x(j)*u(j,k)
       end do
       if (b(k) /= 0.0_dp .and. u(k,k) /= 0.0_dp) then
          x(k)=b(k)/u(k,k); x(k)=x(k)/abs(x(k))
       end if
       b(k)=b(k)+x(k)*u(k,k)
    end do
  end subroutine f_d_upper_right_multiply_linpack

  ! x^T U = b^T
  subroutine f_z_upper_right_multiply_linpack(x,u,b)
    complex(kind=dp), dimension(:,:), intent(in) :: u
    complex(kind=dp), dimension(:), intent(out) :: b
    complex(kind=dp), dimension(:), intent(out) :: x
    !
    integer(kind=int32) :: m,k,j
    !
    m=size(u,1)
    x=(1.0_dp,0.0_dp)
    b=(0.0_dp,0.0_dp)
    do k=1,m
       do j=1,k-1
          b(k)=b(k)+x(j)*u(j,k)
       end do
       if (b(k) /= (0.0_dp,0.0_dp) .and. u(k,k) /= (0.0_dp,0.0_dp)) then
          x(k)=b(k)/u(k,k); x(k)=x(k)/abs(x(k))
       end if
       b(k)=b(k)+x(k)*u(k,k)
    end do
  end subroutine f_z_upper_right_multiply_linpack

  ! x^T U^T = b^T
  subroutine f_d_upper_tr_right_multiply_linpack(x,u,b)
    real(kind=dp), dimension(:,:), intent(in) :: u
    real(kind=dp), dimension(:), intent(out) :: b
    real(kind=dp), dimension(:), intent(out) :: x
    !
    integer(kind=int32) :: m,k,j
    !
    m=size(u,1)
    x=1.0_dp
    b=0.0_dp
    do k=m,1,-1
       do j=k+1,m
          b(k)=b(k)+x(j)*u(k,j)
       end do
       if (b(k) /= 0.0_dp .and. u(k,k) /= 0.0_dp) then
          x(k)=b(k)/u(k,k); x(k)=x(k)/abs(x(k))
       end if
       b(k)=b(k)+x(k)*u(k,k)
    end do
  end subroutine f_d_upper_tr_right_multiply_linpack

  ! x^T U^H = b^T
  subroutine f_z_upper_tr_right_multiply_linpack(x,u,b)
    complex(kind=dp), dimension(:,:), intent(in) :: u
    complex(kind=dp), dimension(:), intent(out) :: b
    complex(kind=dp), dimension(:), intent(out) :: x
    !
    integer(kind=int32) :: m,k,j
    !
    m=size(u,1)
    x=(1.0_dp,0.0_dp)
    b=(0.0_dp,0.0_dp)
    do k=m,1,-1
       do j=k+1,m
          b(k)=b(k)+x(j)*conjg(u(k,j))
       end do
       if (b(k) /= (0.0_dp,0.0_dp) .and. u(k,k) /= (0.0_dp,0.0_dp)) then
          x(k)=b(k)/conjg(u(k,k)); x(k)=x(k)/abs(x(k))
       end if
       b(k)=b(k)+x(k)*conjg(u(k,k))
    end do
  end subroutine f_z_upper_tr_right_multiply_linpack
  
  ! U x = b
  subroutine f_d_upper_left_multiply_linpack(u,x,b)
    real(kind=dp), dimension(:,:), intent(in) :: u
    real(kind=dp), dimension(:), intent(out) :: b
    real(kind=dp), dimension(:), intent(out) :: x
    !
    integer(kind=int32) :: m,k,j
    !
    m=size(u,1)
    x=1.0_dp
    b=0.0_dp
    do j=m,1,-1
       do k=j+1,m
          b(j)=b(j)+u(j,k)*x(k)
       end do
       if (b(j) /= 0.0_dp .and. u(j,j) /= 0.0_dp) then
          x(j)=b(j)/u(j,j); x(j)=x(j)/abs(x(j))
       end if
       b(j)=b(j)+u(j,j)*x(j)
    end do
  end subroutine f_d_upper_left_multiply_linpack

  ! U x = b
  subroutine f_z_upper_left_multiply_linpack(u,x,b)
    complex(kind=dp), dimension(:,:), intent(in) :: u
    complex(kind=dp), dimension(:), intent(out) :: b
    complex(kind=dp), dimension(:), intent(out) :: x
    !
    integer(kind=int32) :: m,k,j
    !
    m=size(u,1)
    x=(1.0_dp,0.0_dp)
    b=(0.0_dp,0.0_dp)
    do j=m,1,-1
       do k=j+1,m
          b(j)=b(j)+u(j,k)*x(k)
       end do
       if (b(j) /= (0.0_dp,0.0_dp) .and. u(j,j) /= (0.0_dp,0.0_dp)) then
          x(j)=b(j)/u(j,j); x(j)=x(j)/abs(x(j))
       end if
       b(j)=b(j)+u(j,j)*x(j)
    end do
  end subroutine f_z_upper_left_multiply_linpack
  
  ! U^T x = b
  subroutine f_d_upper_tr_left_multiply_linpack(u,x,b)
    real(kind=dp), dimension(:,:), intent(in) :: u
    real(kind=dp), dimension(:), intent(out) :: b
    real(kind=dp), dimension(:), intent(out) :: x
    !
    integer(kind=int32) :: m,k,j
    !
    m=size(u,1)
    x=1.0_dp
    b=0.0_dp
    do j=1,m
       do k=1,j-1
          b(j)=b(j)+u(k,j)*x(k)
       end do
       if (b(j) /= 0.0_dp .and. u(j,j) /= 0.0_dp) then
          x(j)=b(j)/u(j,j); x(j)=x(j)/abs(x(j))
       end if
       b(j)=b(j)+u(j,j)*x(j)
    end do
  end subroutine f_d_upper_tr_left_multiply_linpack

  ! U^H x = b
  subroutine f_z_upper_tr_left_multiply_linpack(u,x,b)
    complex(kind=dp), dimension(:,:), intent(in) :: u
    complex(kind=dp), dimension(:), intent(out) :: b
    complex(kind=dp), dimension(:), intent(out) :: x
    !
    integer(kind=int32) :: m,k,j
    !
    m=size(u,1)
    x=(1.0_dp,0.0_dp)
    b=(0.0_dp,0.0_dp)
    do j=1,m
       do k=1,j-1
          b(j)=b(j)+conjg(u(k,j))*x(k)
       end do
       if (b(j) /= (0.0_dp,0.0_dp) .and. u(j,j) /= (0.0_dp,0.0_dp)) then
          x(j)=b(j)/conjg(u(j,j)); x(j)=x(j)/abs(x(j))
       end if
       b(j)=b(j)+conjg(u(j,j))*x(j)
    end do
  end subroutine f_z_upper_tr_left_multiply_linpack
  
  
  
  ! Invert

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

  ! x^T L = b^T
  subroutine f_z_lower_right_invert(x,l,b)
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
  end subroutine f_z_lower_right_invert

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
  subroutine f_z_lower_tr_right_invert(x,l,b)
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
  end subroutine f_z_lower_tr_right_invert

  !
  ! applying inverses from the left.
  !

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

  ! L x = b
  subroutine f_z_lower_left_invert(l,x,b)
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
  end subroutine f_z_lower_left_invert

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
  subroutine f_z_lower_tr_left_invert(l,x,b)
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
  end subroutine f_z_lower_tr_left_invert

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
  subroutine f_z_upper_right_invert(x,r,b)
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
  end subroutine f_z_upper_right_invert

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
  subroutine f_z_upper_tr_right_invert(x,r,b)
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
  end subroutine f_z_upper_tr_right_invert

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
  subroutine f_z_upper_left_invert(r,x,b)
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
  end subroutine f_z_upper_left_invert

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
  subroutine f_z_upper_tr_left_invert(r,x,b)
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
  end subroutine f_z_upper_tr_left_invert

  ! Linpack inversion for condition estimation.

  ! x^T L = b^T
  subroutine f_d_lower_right_invert_linpack(x,l,b)
    real(kind=dp), dimension(:,:), intent(in) :: l
    real(kind=dp), dimension(:), intent(out) :: b
    real(kind=dp), dimension(:), intent(out) :: x

    real(kind=dp) :: tmp
    integer(kind=int32) :: m,k,j
    !
    m=size(l,1)
    b=1.0_dp
    x(m)=b(m)/l(m,m)
    do k=m-1,1,-1
       tmp=0.0_dp
       do j=k+1,m
          tmp = tmp + x(j)*l(j,k)
       end do
       if (tmp /= 0.0_dp) b(k)=-tmp/abs(tmp)
       x(k) = (b(k)-tmp)/l(k,k)
    end do
  end subroutine f_d_lower_right_invert_linpack

  ! x^T L = b^T
  subroutine f_z_lower_right_invert_linpack(x,l,b)
    complex(kind=dp), dimension(:,:), intent(in) :: l
    complex(kind=dp), dimension(:), intent(out) :: b
    complex(kind=dp), dimension(:), intent(out) :: x

    complex(kind=dp) :: tmp
    integer(kind=int32) :: m,k,j
    !
    m=size(l,1)
    b=(1.0_dp,0.0_dp)
    x(m)=b(m)/l(m,m)
    do k=m-1,1,-1
       tmp=(0.0_dp,0.0_dp)
       do j=k+1,m
          tmp = tmp + x(j)*l(j,k)
       end do
       if (tmp /= (0.0_dp,0.0_dp)) b(k)=-tmp/abs(tmp)
       x(k) = (b(k)-tmp)/l(k,k)
    end do
  end subroutine f_z_lower_right_invert_linpack

  ! x^T L^T = b^T
  subroutine f_d_lower_tr_right_invert_linpack(x,l,b)
    real(kind=dp), dimension(:,:), intent(in) :: l
    real(kind=dp), dimension(:), intent(out) :: b
    real(kind=dp), dimension(:), intent(out) :: x
    real(kind=dp) :: tmp
    !
    integer(kind=int32) :: m,k,j
    !
    m=size(l,1)
    b=1.0_dp
    x(1)=b(1)/l(1,1)
    do k=2,m
       tmp = 0.0_dp
       do j=1,k-1
          tmp = tmp + x(j)*l(k,j)
       end do
       if (tmp /= 0.0_dp) b(k)=-tmp/abs(tmp)
       x(k) = (b(k)-tmp)/l(k,k)
    end do
  end subroutine f_d_lower_tr_right_invert_linpack

  ! x^T L^H = b^T
  subroutine f_z_lower_tr_right_invert_linpack(x,l,b)
    complex(kind=dp), dimension(:,:), intent(in) :: l
    complex(kind=dp), dimension(:), intent(out) :: b
    complex(kind=dp), dimension(:), intent(out) :: x
    complex(kind=dp) :: tmp
    !
    integer(kind=int32) :: m,k,j
    !
    m=size(l,1)
    b=(1.0_dp,0.0_dp)
    x(1)=b(1)/conjg(l(1,1))
    do k=2,m
       tmp = (0.0_dp,0.0_dp)
       do j=1,k-1
          tmp = tmp + x(j)*conjg(l(k,j))
       end do
       if (tmp /= (0.0_dp,0.0_dp)) b(k)=-tmp/abs(tmp)
       x(k) = (b(k)-tmp)/conjg(l(k,k))
    end do
  end subroutine f_z_lower_tr_right_invert_linpack

  ! L x = b
  subroutine f_d_lower_left_invert_linpack(l,x,b)
    real(kind=dp), dimension(:,:), intent(in) :: l
    real(kind=dp), dimension(:), intent(out) :: b
    real(kind=dp), dimension(:), intent(out) :: x
    real(kind=dp) :: tmp
    !
    integer(kind=int32) :: m,k,j
    !
    m=size(l,1)
    b=1.0_dp
    x(1)=b(1)/l(1,1)
    do k=2,m
       tmp=0.0_dp
       do j=1,k-1
          tmp = tmp + l(k,j)*x(j)
       end do
       if (tmp /= 0.0_dp) b(k)=-tmp/abs(tmp)
       x(k)=(b(k)-tmp)/l(k,k)
    end do
  end subroutine f_d_lower_left_invert_linpack

  ! L x = b
  subroutine f_z_lower_left_invert_linpack(l,x,b)
    complex(kind=dp), dimension(:,:), intent(in) :: l
    complex(kind=dp), dimension(:), intent(out) :: b
    complex(kind=dp), dimension(:), intent(out) :: x
    complex(kind=dp) :: tmp
    !
    integer(kind=int32) :: m,k,j
    !
    m=size(l,1)
    b=(1.0_dp,0.0_dp)
    x(1)=b(1)/l(1,1)
    do k=2,m
       tmp=(0.0_dp,0.0_dp)
       do j=1,k-1
          tmp = tmp + l(k,j)*x(j)
       end do
       if (tmp /= (0.0_dp,0.0_dp)) b(k)=-tmp/abs(tmp)
       x(k)=(b(k)-tmp)/l(k,k)
    end do
  end subroutine f_z_lower_left_invert_linpack

  ! L^T x = b
  subroutine f_d_lower_tr_left_invert_linpack(l,x,b)
    real(kind=dp), dimension(:,:), intent(in) :: l
    real(kind=dp), dimension(:), intent(out) :: b
    real(kind=dp), dimension(:), intent(out) :: x
    real(kind=dp) :: tmp
    !
    integer(kind=int32) :: m,k,j
    !
    m=size(l,1)
    b=1.0_dp
    x(m)=b(m)/l(m,m)
    do k=m-1,1,-1
       tmp=0.0_dp
       do j=k+1,m
          tmp=tmp+l(j,k)*x(j)
       end do
       if (tmp /= 0.0_dp) b(k)=-tmp/abs(tmp)
       x(k)=(b(k)-tmp)/l(k,k)
    end do
  end subroutine f_d_lower_tr_left_invert_linpack

  ! L^H x = b
  subroutine f_z_lower_tr_left_invert_linpack(l,x,b)
    complex(kind=dp), dimension(:,:), intent(in) :: l
    complex(kind=dp), dimension(:), intent(out) :: b
    complex(kind=dp), dimension(:), intent(out) :: x
    complex(kind=dp) :: tmp
    !
    integer(kind=int32) :: m,k,j
    !
    m=size(l,1)
    b=(1.0_dp,0.0_dp)
    x(m)=b(m)/conjg(l(m,m))
    do k=m-1,1,-1
       tmp=(0.0_dp,0.0_dp)
       do j=k+1,m
          tmp=tmp+conjg(l(j,k))*x(j)
       end do
       if (tmp /= (0.0_dp,0.0_dp)) b(k)=-tmp/abs(tmp)
       x(k)=(b(k)-tmp)/conjg(l(k,k))
    end do
  end subroutine f_z_lower_tr_left_invert_linpack

  ! x^T U = b^T
  subroutine f_d_upper_right_invert_linpack(x,u,b)
    real(kind=dp), dimension(:,:), intent(in) :: u
    real(kind=dp), dimension(:), intent(out) :: b
    real(kind=dp), dimension(:), intent(out) :: x
    real(kind=dp) :: tmp
    !
    integer(kind=int32) :: m,k,j
    !
    m=size(u,1)
    b=1.0_dp
    x(1)=b(1)/u(1,1)
    do k=2,m
       tmp=0.0_dp
       do j=1,k-1
          tmp = tmp + x(j)*u(j,k)
       end do
       if (tmp /= 0.0_dp) b(k)=-tmp/abs(tmp)
       x(k)=(b(k)-tmp)/u(k,k)
    end do
  end subroutine f_d_upper_right_invert_linpack

  ! x^T U = b^T
  subroutine f_z_upper_right_invert_linpack(x,u,b)
    complex(kind=dp), dimension(:,:), intent(in) :: u
    complex(kind=dp), dimension(:), intent(out) :: b
    complex(kind=dp), dimension(:), intent(out) :: x
    complex(kind=dp) :: tmp
    !
    integer(kind=int32) :: m,k,j
    !
    m=size(u,1)
    b=(1.0_dp,0.0_dp)
    x(1)=b(1)/u(1,1)
    do k=2,m
       tmp=(0.0_dp,0.0_dp)
       do j=1,k-1
          tmp = tmp + x(j)*u(j,k)
       end do
       if (tmp /= (0.0_dp,0.0_dp)) b(k)=-tmp/abs(tmp)
       x(k)=(b(k)-tmp)/u(k,k)
    end do
  end subroutine f_z_upper_right_invert_linpack

  ! x^T U^T = b^T
  subroutine f_d_upper_tr_right_invert_linpack(x,u,b)
    real(kind=dp), dimension(:,:), intent(in) :: u
    real(kind=dp), dimension(:), intent(out) :: b
    real(kind=dp), dimension(:), intent(out) :: x
    real(kind=dp) :: tmp
    !
    integer(kind=int32) :: m,k,j
    !
    m=size(u,1)
    b=1.0_dp
    x(m)=b(m)/u(m,m)
    do k=m-1,1,-1
       tmp=0.0_dp
       do j=k+1,m
          tmp = tmp + x(j)*u(k,j)
       end do
       if (tmp /= 0.0_dp) b(k)=-tmp/abs(tmp)
       x(k) = (b(k)-tmp)/u(k,k)
    end do
  end subroutine f_d_upper_tr_right_invert_linpack

  ! x^T U^H = b^T
  subroutine f_z_upper_tr_right_invert_linpack(x,u,b)
    complex(kind=dp), dimension(:,:), intent(in) :: u
    complex(kind=dp), dimension(:), intent(out) :: b
    complex(kind=dp), dimension(:), intent(out) :: x
    complex(kind=dp) :: tmp
    !
    integer(kind=int32) :: m,k,j
    !
    m=size(u,1)
    b=(1.0_dp,0.0_dp)
    x(m)=b(m)/conjg(u(m,m))
    do k=m-1,1,-1
       tmp=(0.0_dp,0.0_dp)
       do j=k+1,m
          tmp = tmp + x(j)*conjg(u(k,j))
       end do
       if (tmp /= (0.0_dp,0.0_dp)) b(k)=-tmp/abs(tmp)
       x(k) = (b(k)-tmp)/conjg(u(k,k))
    end do
  end subroutine f_z_upper_tr_right_invert_linpack

  ! U x = b
  subroutine f_d_upper_left_invert_linpack(u,x,b)
    real(kind=dp), dimension(:,:), intent(in) :: u
    real(kind=dp), dimension(:), intent(out) :: b
    real(kind=dp), dimension(:), intent(out) :: x
    real(kind=dp) :: tmp
    !
    integer(kind=int32) :: m,k,j
    !
    m=size(u,1)
    b=1.0_dp
    x(m)=b(m)/u(m,m)
    do k=m-1,1,-1
       tmp=0.0_dp
       do j=k+1,m
          tmp = tmp + u(k,j)*x(j)
       end do
       if (tmp /= 0.0_dp) b(k)=-tmp/abs(tmp)
       x(k)=(b(k)-tmp)/u(k,k)
    end do
  end subroutine f_d_upper_left_invert_linpack

  subroutine f_z_upper_left_invert_linpack(u,x,b)
    complex(kind=dp), dimension(:,:), intent(in) :: u
    complex(kind=dp), dimension(:), intent(out) :: b
    complex(kind=dp), dimension(:), intent(out) :: x
    complex(kind=dp) :: tmp
    !
    integer(kind=int32) :: m,k,j
    !
    m=size(u,1)
    b=(1.0_dp,0.0_dp)
    x(m)=b(m)/u(m,m)
    do k=m-1,1,-1
       tmp=(0.0_dp,0.0_dp)
       do j=k+1,m
          tmp = tmp + u(k,j)*x(j)
       end do
       if (tmp /= (0.0_dp,0.0_dp)) b(k)=-tmp/abs(tmp)
       x(k)=(b(k)-tmp)/u(k,k)
    end do
  end subroutine f_z_upper_left_invert_linpack

  ! U^T x = b
  subroutine f_d_upper_tr_left_invert_linpack(u,x,b)
    real(kind=dp), dimension(:,:), intent(in) :: u
    real(kind=dp), dimension(:), intent(out) :: b
    real(kind=dp), dimension(:), intent(out) :: x
    real(kind=dp) :: tmp
    !
    integer(kind=int32) :: m,k,j
    !
    m=size(u,1)
    b=1.0_dp
    x(1)=b(1)/u(1,1)
    do k=2,m
       tmp=0.0_dp
       do j=1,k-1
          tmp=tmp+u(j,k)*x(j)
       end do
       if (tmp /= 0.0_dp) b(k)=-tmp/abs(tmp)
       x(k)=(b(k)-tmp)/u(k,k)
    end do
  end subroutine f_d_upper_tr_left_invert_linpack

  ! U^H x = b
  subroutine f_z_upper_tr_left_invert_linpack(u,x,b)
    complex(kind=dp), dimension(:,:), intent(in) :: u
    complex(kind=dp), dimension(:), intent(out) :: b
    complex(kind=dp), dimension(:), intent(out) :: x
    complex(kind=dp) :: tmp
    !
    integer(kind=int32) :: m,k,j
    !
    m=size(u,1)
    b=(1.0_dp,0.0_dp)
    x(1)=b(1)/conjg(u(1,1))
    do k=2,m
       tmp=(0.0_dp,0.0_dp)
       do j=1,k-1
          tmp=tmp+conjg(u(j,k))*x(j)
       end do
       if (tmp /= (0.0_dp,0.0_dp)) b(k)=-tmp/abs(tmp)
       x(k)=(b(k)-tmp)/conjg(u(k,k))
    end do
  end subroutine f_z_upper_tr_left_invert_linpack
  

  ! finding zeros on the diagonal

  integer(kind=int32) function d_first_zero_diagonal(a,tol) result(d)
    real(kind=dp), dimension(:,:), intent(in) :: a
    real(kind=dp), intent(in) :: tol
    !
    integer(kind=int32) :: j, n
    n=min(size(a,1), size(a,2))
    d=1
    do j=1,n
       if (abs(a(j,j)) <= tol) return
       d=d+1
    end do
  end function d_first_zero_diagonal

  integer(kind=int32) function z_first_zero_diagonal(a,tol) result(d)
    complex(kind=dp), dimension(:,:), intent(in) :: a
    real(kind=dp), intent(in) :: tol
    !
    integer(kind=int32) :: j, n
    n=min(size(a,1), size(a,2))
    d=1
    do j=1,n
       if (abs(a(j,j)) <= tol) return
       d=d+1
    end do
  end function z_first_zero_diagonal

  ! Return the first zero diagonal from the bottom right.
  integer(kind=int32) function d_reverse_first_zero_diagonal(a,tol) result(d)
    real(kind=dp), dimension(:,:), intent(in) :: a
    real(kind=dp), intent(in) :: tol
    !
    integer(kind=int32) :: j, n
    n=min(size(a,1), size(a,2))
    d=n
    do j=n,1,-1
       if (abs(a(j,j)) <= tol) return
       d=d-1
    end do
  end function d_reverse_first_zero_diagonal

  integer(kind=int32) function z_reverse_first_zero_diagonal(a,tol) result(d)
    complex(kind=dp), dimension(:,:), intent(in) :: a
    real(kind=dp), intent(in) :: tol
    !
    integer(kind=int32) :: j, n
    n=min(size(a,1), size(a,2))
    d=n
    do j=n,1,-1
       if (abs(a(j,j)) <= tol) return
       d=d-1
    end do
  end function z_reverse_first_zero_diagonal

end module mod_triangular
