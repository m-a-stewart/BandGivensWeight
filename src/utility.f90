module utility
use prec
implicit none

interface lower_right_invert
   module procedure f_c_lower_right_invert, f_d_lower_right_invert
end interface lower_right_invert

interface lower_tr_right_invert
   module procedure f_c_lower_tr_right_invert, f_d_lower_tr_right_invert
end interface lower_tr_right_invert

interface first_zero_diagonal
   module procedure c_first_zero_diagonal, d_first_zero_diagonal
end interface first_zero_diagonal

interface maxabs
   module procedure c_maxabs, d_maxabs, c_maxabs_v, d_maxabs_v
end interface maxabs

interface norm2
   module procedure c_norm2_v, d_norm2_v
end interface norm2

interface normf
   module procedure c_normf, d_normf
end interface normf

interface print_matrix
   module procedure d_print_matrix, c_print_matrix, i_print_matrix
end interface print_matrix

interface random_complex
   module procedure random_complex_m, random_complex_v, random_complex_s
end interface random_complex

public :: f_c_lower_right_invert, f_d_lower_right_invert, f_c_lower_tr_right_invert, &
          f_d_lower_tr_right_invert, c_first_zero_diagonal, &
          d_first_zero_diagonal, c_maxabs, d_maxabs, d_print_matrix, c_print_matrix

contains

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

  integer(kind=int32) function c_first_zero_diagonal(a,tol) result(d)
    complex(kind=dp), dimension(:,:), intent(in) :: a
    real(kind=dp), intent(in) :: tol
    !
    integer(kind=int32) :: j, n
    n=min(size(a,1), size(a,2))
    d=1
    do j=1,n
       if (abs(a(j,j)) < tol) then
          return
       end if
       d=d+1
    end do
  end function c_first_zero_diagonal

  integer(kind=int32) function d_first_zero_diagonal(a,tol) result(d)
    real(kind=dp), dimension(:,:), intent(in) :: a
    real(kind=dp), intent(in) :: tol
    !
    integer(kind=int32) :: j, n
    n=min(size(a,1), size(a,2))
    d=1
    do j=1,n
       if (abs(a(j,j)) < tol) then
          return
       end if
       d=d+1
    end do
  end function d_first_zero_diagonal

  real(kind=dp) function c_maxabs(a) result(x)
    complex(kind=dp), dimension(:,:), intent(in) :: a
    !
    integer(kind=int32) :: j, k, m, n
    real(kind=dp) :: y
    m=size(a,1)
    n=size(a,2)
    x=0.0_dp
    do j=1,m
       do k=1,n
          y=abs(a(j,k))
          if (y > x) then
             x=y
          end if
       end do
    end do
  end function c_maxabs

  real(kind=dp) function d_maxabs(a) result(x)
    real(kind=dp), dimension(:,:), intent(in) :: a
    !
    integer(kind=int32) :: j, k, m, n
    real(kind=dp) :: y
    m=size(a,1)
    n=size(a,2)
    x=0.0_dp
    do j=1,m
       do k=1,n
          y=abs(a(j,k))
          if (y > x) then
             x=y
          end if
       end do
    end do
  end function d_maxabs

  real(kind=dp) function c_maxabs_v(a) result(x)
    complex(kind=dp), dimension(:), intent(in) :: a
    !
    integer(kind=int32) :: j, m
    real(kind=dp) :: y
    m=size(a)
    x=0.0_dp
    do j=1,m
       y=abs(a(j))
       if (y > x) then
          x=y
       end if
    end do
  end function c_maxabs_v

  real(kind=dp) function d_maxabs_v(a) result(x)
    real(kind=dp), dimension(:), intent(in) :: a
    !
    integer(kind=int32) :: j, m
    real(kind=dp) :: y
    m=size(a)
    x=0.0_dp
    do j=1,m
       y=abs(a(j))
       if (y > x) then
          x=y
       end if
    end do
  end function d_maxabs_v

  real(kind=dp) function c_norm2_v(a) result(x)
    complex(kind=dp), dimension(:), intent(in) :: a
    !
    integer(kind=int32) :: j, m
    real(kind=dp) :: y
    m=size(a)
    y=c_maxabs_v(a)
    if (y==0.0_dp) then
       x=0.0_dp
       return
    end if
    y=(2.0_dp)**(floor(log(y)/log(2.0_dp)))
    x=0.0_dp
    do j=1,m
       x=x+abs(a(j)/y)**2
    end do
    x=sqrt(x)*y
  end function c_norm2_v

  real(kind=dp) function d_norm2_v(a) result(x)
    real(kind=dp), dimension(:), intent(in) :: a
    !
    integer(kind=int32) :: j, m
    real(kind=dp) :: y, tmp
    m=size(a)
    y=d_maxabs_v(a)
    if (y==0.0_dp) then
       x=0.0_dp
       return
    end if
    y=(2.0_dp)**(floor(log(y)/log(2.0_dp)))
    x=0.0_dp
    if (y > 1e150_dp .or. y < 1e-150_dp) then
       do j=1,m
          tmp=a(j)/y
          x=x+tmp*tmp
       end do
       x=sqrt(x)*y
    else
       do j=1,m
          tmp=a(j)
          x=x+tmp*tmp
       end do
       x=sqrt(x)
    end if
  end function d_norm2_v

  real(kind=dp) function c_normf(a) result(x)
    complex(kind=dp), dimension(:,:), intent(in) :: a
    !
    integer(kind=int32) :: j, k, m, n
    real(kind=dp) :: y, tmp
    m=size(a,1)
    n=size(a,2)
    y=c_maxabs(a)
    if (y==0.0_dp) then
       x=0.0_dp
       return
    end if
    y=(2.0_dp)**(floor(log(y)/log(2.0_dp)))
    x=0.0_dp
    if (y > 1e150_dp .or. y < 1e-150_dp) then
       do j=1,m
          do k=1,n
             tmp=abs(a(j,k)/y)
             x=x+tmp*tmp
          end do
       end do
       x=sqrt(x)*y
    else
       do j=1,m
          do k=1,n
             tmp=abs(a(j,k))
             x=x+tmp*tmp
          end do
       end do
       x=sqrt(x)
    end if
  end function c_normf

  real(kind=dp) function d_normf(a) result(x)
    real(kind=dp), dimension(:,:), intent(in) :: a
    !
    integer(kind=int32) :: j, k, m, n
    real(kind=dp) :: y, tmp
    m=size(a,1)
    n=size(a,2)
    y=d_maxabs(a)
    if (y==0.0_dp) then
       x=0.0_dp
       return
    end if
    y=(2.0_dp)**(floor(log(y)/log(2.0_dp)))
    x=0.0_dp
    if (y > 1e150_dp .or. y < 1e-150_dp) then
       do j=1,m
          do k=1,n
             tmp = a(j,k)/y
             x=x+tmp*tmp
          end do
       end do
       x=sqrt(x)*y
    else
       do j=1,m
          do k=1,n
             tmp = a(j,k)
             x=x+tmp*tmp
          end do
       end do
       x=sqrt(x)
    end if
  end function d_normf

  subroutine d_print_matrix(a)
    real(kind=dp), dimension(:,:), intent(in) :: a
    integer(kind=int32) :: j,k,m,n
    m=size(a,1)
    n=size(a,2)
    do j=1,m
       do k=1,n
          write(*,'(ES12.4, ",  ")',advance='no') a(j,k)
       end do
       write(*,*)
    end do
  end subroutine d_print_matrix

  subroutine c_print_matrix(a)
    complex(kind=dp), dimension(:,:), intent(in) :: a
    integer(kind=int32) :: j,k,m,n
    m=size(a,1)
    n=size(a,2)
    do j=1,m
       do k=1,n
          write(*,'(ES12.4, A, ES11.4,A)',advance='no') real(a(j,k)), ' +i*',aimag(a(j,k)), '   '
       end do
       write(*,*)
    end do
  end subroutine c_print_matrix

  subroutine i_print_matrix(a)
    integer(kind=int32), dimension(:,:), intent(in) :: a
    integer(kind=int32) :: j,k,m,n
    m=size(a,1)
    n=size(a,2)
    do j=1,m
       do k=1,n
          write(*,'(i4, "  ")',advance='no') a(j,k)
       end do
       write(*,*)
    end do
  end subroutine i_print_matrix

  subroutine random_complex_m(a)
    complex(kind=dp), dimension(:,:), intent(out) :: a
    real(kind=dp) :: x,y
    integer(kind=int32) :: j, k, m, n
    m=size(a,1)
    n=size(a,2)
    do j=1,m
       do k=1,n
          call random_number(x)
          call random_number(y)
          a(j,k)=cmplx(x,y)
       end do
    end do
  end subroutine random_complex_m

  subroutine random_complex_v(a)
    complex(kind=dp), dimension(:), intent(out) :: a
    real(kind=dp) :: x,y
    integer(kind=int32) :: j, m
    m=size(a)
    do j=1,m
       call random_number(x)
       call random_number(y)
       a(j)=cmplx(x,y)
    end do
  end subroutine random_complex_v

  subroutine random_complex_s(a)
    complex(kind=dp), intent(out) :: a
    real(kind=dp) :: x,y
    call random_number(x)
    call random_number(y)
    a=cmplx(x,y)
  end subroutine random_complex_s
    
end module utility
