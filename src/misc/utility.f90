module mod_utility
  use mod_prec
  implicit none

  private

  public :: maxabs, c_maxabs, d_maxabs, c_maxabs_v, d_maxabs_v, &
       norm2, c_norm2_v, d_norm2_v, &
       normf, c_normf, d_normf

  public :: print_matrix, d_print_matrix, c_print_matrix, i_print_matrix

  public :: random_matrix, d_random_matrix, d_v_random_matrix, d_s_random_matrix, &
       c_random_matrix, c_v_random_matrix, c_s_random_matrix

  public :: ip_transpose, d_ip_transpose, c_ip_transpose

  public :: includes_nan, d_includes_nan, d_v_includes_nan, c_includes_nan, c_v_includes_nan

  public :: c_delta, d_delta

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

  interface random_matrix
     module procedure d_random_matrix, d_v_random_matrix, d_s_random_matrix, &
          c_random_matrix, c_v_random_matrix, c_s_random_matrix
  end interface random_matrix

  interface ip_transpose
     module procedure d_ip_transpose, c_ip_transpose
  end interface ip_transpose

  interface includes_nan
     module procedure d_includes_nan, d_v_includes_nan, c_includes_nan, c_v_includes_nan
  end interface includes_nan

  real(kind=dp), parameter :: d_random_shift=-1.0_dp, d_random_scale=2.0_dp
  complex(kind=dp), parameter :: c_random_shift=(-1.0_dp,-1.0_dp), c_random_scale=(2.0_dp,0.0_dp)

contains

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
          if (y > x .or. y /= y) then
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
          if (y > x .or. y /= y) then
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
       if (y > x .or. y /= y) then
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
       if (y > x .or. y /= y) then
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
          !          write(*,'(ES12.4, A, ES11.4,A)',advance='no') real(a(j,k)), ' +i*',aimag(a(j,k)), '   '
          write(*,'(A, ES12.4, A, ES11.4,A)',advance='no') '(', real(a(j,k)), ', ',aimag(a(j,k)), ' )    '
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

  subroutine d_random_matrix(a)
    real(kind=dp), dimension(:,:), intent(out) :: a
    call random_number(a)
    a=d_random_scale*a+d_random_shift
  end subroutine d_random_matrix

  subroutine d_v_random_matrix(a)
    real(kind=dp), dimension(:), intent(out) :: a
    call random_number(a)
    a=d_random_scale*a+d_random_shift
  end subroutine d_v_random_matrix

  subroutine d_s_random_matrix(a)
    real(kind=dp), intent(out) :: a
    call random_number(a)
    a=d_random_scale*a+d_random_shift
  end subroutine d_s_random_matrix


  subroutine c_random_matrix(a)
    complex(kind=dp), dimension(:,:), intent(out) :: a
    real(kind=dp) :: x,y
    integer(kind=int32) :: j, k
    do j=1,size(a,1)
       do k=1,size(a,2)
          call random_number(x)
          call random_number(y)
          a(j,k)=c_random_scale*cmplx(x,y)+c_random_shift
       end do
    end do
  end subroutine c_random_matrix

  subroutine c_v_random_matrix(a)
    complex(kind=dp), dimension(:), intent(out) :: a
    real(kind=dp) :: x,y
    integer(kind=int32) :: j
    do j=1,size(a)
       call random_number(x)
       call random_number(y)
       a(j)= c_random_scale*cmplx(x,y)+c_random_shift
    end do
  end subroutine c_v_random_matrix

  subroutine c_s_random_matrix(a)
    complex(kind=dp), intent(out) :: a
    real(kind=dp) :: x,y
    call random_number(x)
    call random_number(y)
    a=c_random_scale*cmplx(x,y)+c_random_shift
  end subroutine c_s_random_matrix

  real(kind=dp) function d_delta(j,k)
    integer(kind=int32), intent(in) :: j,k
    if (j==k) then
       d_delta=1.0_dp
    else
       d_delta=0.0_dp
    end if
  end function d_delta

  complex(kind=dp) function c_delta(j,k)
    integer(kind=int32), intent(in) :: j,k
    if (j==k) then
       c_delta=(1.0_dp,0.0_dp)
    else
       c_delta=(0.0_dp,0.0_dp)
    end if
  end function c_delta

  type(logical) function d_includes_nan(a) result(nanflag)
    real(kind=dp), dimension(:,:) :: a
    integer(kind=int32) :: m,n,j,k
    nanflag=.false.
    m=size(a,1)
    n=size(a,2)
    do j=1,m
       do k=1,n
          if (a(j,k) /= a(j,k)) then
             nanflag=.true.
          end if
       end do
    end do
  end function d_includes_nan

  type(logical) function d_v_includes_nan(a) result(nanflag)
    real(kind=dp), dimension(:) :: a
    integer(kind=int32) :: n,j
    nanflag=.false.
    n=size(a)
    do j=1,n
       if (a(j) /= a(j)) then
          nanflag=.true.
       end if
    end do
  end function d_v_includes_nan

  type(logical) function c_includes_nan(a) result(nanflag)
    complex(kind=dp), dimension(:,:) :: a
    integer(kind=int32) :: m,n,j,k
    nanflag=.false.
    m=size(a,1)
    n=size(a,2)
    do j=1,m
       do k=1,n
          if (a(j,k) /= a(j,k)) then
             nanflag=.true.
          end if
       end do
    end do
  end function c_includes_nan

  type(logical) function c_v_includes_nan(a) result(nanflag)
    complex(kind=dp), dimension(:) :: a
    integer(kind=int32) :: n,j
    nanflag=.false.
    n=size(a)
    do j=1,n
       if (a(j) /= a(j)) then
          nanflag=.true.
       end if
    end do
  end function c_v_includes_nan

  subroutine d_ip_transpose(a)
    real(kind=dp), dimension(:,:), intent(inout) :: a
    real(kind=dp) :: tmp
    integer(kind=int32) :: n, j,k
    n=size(a,1)
    do j=2,n
       do k=1,j-1
          tmp=a(j,k)
          a(j,k)=a(k,j)
          a(k,j)=tmp
       end do
    end do
  end subroutine d_ip_transpose

  subroutine c_ip_transpose(a)
    complex(kind=dp), dimension(:,:), intent(inout) :: a
    complex(kind=dp) :: tmp
    integer(kind=int32) :: n, j,k
    n=size(a,1)
    do j=2,n
       do k=1,j-1
          tmp=a(j,k)
          a(j,k)=conjg(a(k,j))
          a(k,j)=conjg(tmp)
       end do
    end do
    do j=1,n
       a(j,j)=conjg(a(j,j))
    end do
  end subroutine c_ip_transpose

end module mod_utility
