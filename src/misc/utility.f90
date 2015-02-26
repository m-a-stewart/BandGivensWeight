module mod_utility
  use mod_prec
  implicit none

  private

  public :: maybe_deallocate, not_all_allocated

  public :: maxabs, c_maxabs, d_maxabs, c_maxabs_v, d_maxabs_v, &
       norm2, c_norm2_v, d_norm2_v, &
       normf, c_normf, d_normf

  public :: print_matrix, d_print_matrix, c_print_matrix, i_print_matrix

  public :: random_matrix_to, d_random_matrix_to, d_v_random_matrix_to, d_s_random_matrix_to, &
       c_random_matrix_to, c_v_random_matrix_to, c_s_random_matrix_to, &
       d_random_matrix, d_random_vector, c_random_matrix, c_random_vector

  public :: ip_transpose, d_ip_transpose, c_ip_transpose

  public :: c_delta, d_delta

  public :: equals_option, i_equals_option

  interface maybe_deallocate
     module procedure d_maybe_deallocate2, c_maybe_deallocate2, i_maybe_deallocate2, &
          d_maybe_deallocate1, c_maybe_deallocate1, i_maybe_deallocate1
  end interface maybe_deallocate

  interface not_all_allocated
     module procedure d_not_all_allocated2, c_not_all_allocated2, i_not_all_allocated2, &
          d_not_all_allocated1, c_not_all_allocated1, i_not_all_allocated1
  end interface not_all_allocated

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

  interface random_matrix_to
     module procedure d_random_matrix_to, d_v_random_matrix_to, d_s_random_matrix_to, &
          c_random_matrix_to, c_v_random_matrix_to, c_s_random_matrix_to
  end interface random_matrix_to

  interface ip_transpose
     module procedure d_ip_transpose, c_ip_transpose
  end interface ip_transpose

  interface equals_option
     module procedure i_equals_option
  end interface equals_option

  real(kind=dp), parameter :: d_random_shift=-1.0_dp, d_random_scale=2.0_dp
  complex(kind=dp), parameter :: c_random_shift=(-1.0_dp,-1.0_dp), c_random_scale=(2.0_dp,0.0_dp)

contains

  logical function d_not_all_allocated2(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10) result(x)
    real(kind=dp), dimension(:,:), allocatable :: a1
    real(kind=dp), dimension(:,:), optional, allocatable :: a2, a3, a4, a5, &
         a6, a7, a8, a9, a10
    x=.false.
    if (.not. allocated(a1)) goto 1
    if (.not. present(a2)) return
    if (.not. allocated(a2)) goto 1
    if (.not. present(a3)) return
    if (.not. allocated(a3)) goto 1
    if (.not. present(a4)) return
    if (.not. allocated(a4)) goto 1
    if (.not. present(a5)) return
    if (.not. allocated(a5)) goto 1
    if (.not. present(a6)) return
    if (.not. allocated(a6)) goto 1
    if (.not. present(a7)) return
    if (.not. allocated(a7)) goto 1
    if (.not. present(a8)) return
    if (.not. allocated(a8)) goto 1
    if (.not. present(a9)) return
    if (.not. allocated(a9)) goto 1
    if (.not. present(a10)) return
    if (allocated(a10)) return
1   x=.true.
  end function d_not_all_allocated2

  logical function c_not_all_allocated2(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10) result(x)
    complex(kind=dp), dimension(:,:), allocatable :: a1
    complex(kind=dp), dimension(:,:), optional, allocatable :: a2, a3, a4, a5, &
         a6, a7, a8, a9, a10
    x=.false.
    if (.not. allocated(a1)) goto 1
    if (.not. present(a2)) return
    if (.not. allocated(a2)) goto 1
    if (.not. present(a3)) return
    if (.not. allocated(a3)) goto 1
    if (.not. present(a4)) return
    if (.not. allocated(a4)) goto 1
    if (.not. present(a5)) return
    if (.not. allocated(a5)) goto 1
    if (.not. present(a6)) return
    if (.not. allocated(a6)) goto 1
    if (.not. present(a7)) return
    if (.not. allocated(a7)) goto 1
    if (.not. present(a8)) return
    if (.not. allocated(a8)) goto 1
    if (.not. present(a9)) return
    if (.not. allocated(a9)) goto 1
    if (.not. present(a10)) return
    if (allocated(a10)) return
1   x=.true.
  end function c_not_all_allocated2

  logical function i_not_all_allocated2(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10) result(x)
    integer(kind=int32), dimension(:,:), allocatable :: a1
    integer(kind=int32), dimension(:,:), optional, allocatable :: a2, a3, a4, a5, &
         a6, a7, a8, a9, a10
    x=.false.
    if (.not. allocated(a1)) goto 1
    if (.not. present(a2)) return
    if (.not. allocated(a2)) goto 1
    if (.not. present(a3)) return
    if (.not. allocated(a3)) goto 1
    if (.not. present(a4)) return
    if (.not. allocated(a4)) goto 1
    if (.not. present(a5)) return
    if (.not. allocated(a5)) goto 1
    if (.not. present(a6)) return
    if (.not. allocated(a6)) goto 1
    if (.not. present(a7)) return
    if (.not. allocated(a7)) goto 1
    if (.not. present(a8)) return
    if (.not. allocated(a8)) goto 1
    if (.not. present(a9)) return
    if (.not. allocated(a9)) goto 1
    if (.not. present(a10)) return
    if (allocated(a10)) return
1   x=.true.
  end function i_not_all_allocated2

  logical function d_not_all_allocated1(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10) result(x)
    real(kind=dp), dimension(:), allocatable :: a1
    real(kind=dp), dimension(:), optional, allocatable :: a2, a3, a4, a5, &
         a6, a7, a8, a9, a10
    x=.false.
    if (.not. allocated(a1)) goto 1
    if (.not. present(a2)) return
    if (.not. allocated(a2)) goto 1
    if (.not. present(a3)) return
    if (.not. allocated(a3)) goto 1
    if (.not. present(a4)) return
    if (.not. allocated(a4)) goto 1
    if (.not. present(a5)) return
    if (.not. allocated(a5)) goto 1
    if (.not. present(a6)) return
    if (.not. allocated(a6)) goto 1
    if (.not. present(a7)) return
    if (.not. allocated(a7)) goto 1
    if (.not. present(a8)) return
    if (.not. allocated(a8)) goto 1
    if (.not. present(a9)) return
    if (.not. allocated(a9)) goto 1
    if (.not. present(a10)) return
    if (allocated(a10)) return
1   x=.true.
  end function d_not_all_allocated1

  logical function c_not_all_allocated1(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10) result(x)
    complex(kind=dp), dimension(:), allocatable :: a1
    complex(kind=dp), dimension(:), optional, allocatable :: a2, a3, a4, a5, &
         a6, a7, a8, a9, a10
    x=.false.
    if (.not. allocated(a1)) goto 1
    if (.not. present(a2)) return
    if (.not. allocated(a2)) goto 1
    if (.not. present(a3)) return
    if (.not. allocated(a3)) goto 1
    if (.not. present(a4)) return
    if (.not. allocated(a4)) goto 1
    if (.not. present(a5)) return
    if (.not. allocated(a5)) goto 1
    if (.not. present(a6)) return
    if (.not. allocated(a6)) goto 1
    if (.not. present(a7)) return
    if (.not. allocated(a7)) goto 1
    if (.not. present(a8)) return
    if (.not. allocated(a8)) goto 1
    if (.not. present(a9)) return
    if (.not. allocated(a9)) goto 1
    if (.not. present(a10)) return
    if (allocated(a10)) return
1   x=.true.
  end function c_not_all_allocated1

  logical function i_not_all_allocated1(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10) result(x)
    integer(kind=int32), dimension(:), allocatable :: a1
    integer(kind=int32), dimension(:), optional, allocatable :: a2, a3, a4, a5, &
         a6, a7, a8, a9, a10
    x=.false.
    if (.not. allocated(a1)) goto 1
    if (.not. present(a2)) return
    if (.not. allocated(a2)) goto 1
    if (.not. present(a3)) return
    if (.not. allocated(a3)) goto 1
    if (.not. present(a4)) return
    if (.not. allocated(a4)) goto 1
    if (.not. present(a5)) return
    if (.not. allocated(a5)) goto 1
    if (.not. present(a6)) return
    if (.not. allocated(a6)) goto 1
    if (.not. present(a7)) return
    if (.not. allocated(a7)) goto 1
    if (.not. present(a8)) return
    if (.not. allocated(a8)) goto 1
    if (.not. present(a9)) return
    if (.not. allocated(a9)) goto 1
    if (.not. present(a10)) return
    if (allocated(a10)) return
1   x=.true.
  end function i_not_all_allocated1
  

  subroutine d_maybe_deallocate2(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10)
    real(kind=dp), dimension(:,:), allocatable :: a1
    real(kind=dp), dimension(:,:), optional, allocatable :: a2, a3, a4, a5, &
         a6, a7, a8, a9, a10
    if (allocated(a1)) deallocate(a1)
    if (.not. present(a2)) go to 1
    if (allocated(a2)) deallocate(a2)
    if (.not. present(a3)) go to 1
    if (allocated(a3)) deallocate(a3)
    if (.not. present(a4)) go to 1
    if (allocated(a4)) deallocate(a4)
    if (.not. present(a5)) go to 1
    if (allocated(a5)) deallocate(a5)
    if (.not. present(a6)) go to 1
    if (allocated(a6)) deallocate(a6)
    if (.not. present(a7)) go to 1
    if (allocated(a7)) deallocate(a7)
    if (.not. present(a8)) go to 1
    if (allocated(a8)) deallocate(a8)
    if (.not. present(a9)) go to 1
    if (allocated(a9)) deallocate(a9)
    if (.not. present(a10)) go to 1
    if (allocated(a10)) deallocate(a10)
1   return
  end subroutine d_maybe_deallocate2

  subroutine c_maybe_deallocate2(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10)
    complex(kind=dp), dimension(:,:), allocatable :: a1
    complex(kind=dp), dimension(:,:), optional, allocatable :: a2, a3, a4, a5, &
         a6, a7, a8, a9, a10
    if (allocated(a1)) deallocate(a1)
    if (.not. present(a2)) go to 1
    if (allocated(a2)) deallocate(a2)
    if (.not. present(a3)) go to 1
    if (allocated(a3)) deallocate(a3)
    if (.not. present(a4)) go to 1
    if (allocated(a4)) deallocate(a4)
    if (.not. present(a5)) go to 1
    if (allocated(a5)) deallocate(a5)
    if (.not. present(a6)) go to 1
    if (allocated(a6)) deallocate(a6)
    if (.not. present(a7)) go to 1
    if (allocated(a7)) deallocate(a7)
    if (.not. present(a8)) go to 1
    if (allocated(a8)) deallocate(a8)
    if (.not. present(a9)) go to 1
    if (allocated(a9)) deallocate(a9)
    if (.not. present(a10)) go to 1
    if (allocated(a10)) deallocate(a10)
1   return
  end subroutine c_maybe_deallocate2

  subroutine i_maybe_deallocate2(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10)
    integer(kind=int32), dimension(:,:), allocatable :: a1
    integer(kind=int32), dimension(:,:), optional, allocatable :: a2, a3, a4, a5, &
         a6, a7, a8, a9, a10
    if (allocated(a1)) deallocate(a1)
    if (.not. present(a2)) go to 1
    if (allocated(a2)) deallocate(a2)
    if (.not. present(a3)) go to 1
    if (allocated(a3)) deallocate(a3)
    if (.not. present(a4)) go to 1
    if (allocated(a4)) deallocate(a4)
    if (.not. present(a5)) go to 1
    if (allocated(a5)) deallocate(a5)
    if (.not. present(a6)) go to 1
    if (allocated(a6)) deallocate(a6)
    if (.not. present(a7)) go to 1
    if (allocated(a7)) deallocate(a7)
    if (.not. present(a8)) go to 1
    if (allocated(a8)) deallocate(a8)
    if (.not. present(a9)) go to 1
    if (allocated(a9)) deallocate(a9)
    if (.not. present(a10)) go to 1
    if (allocated(a10)) deallocate(a10)
1   return
  end subroutine i_maybe_deallocate2

  subroutine d_maybe_deallocate1(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10)
    real(kind=dp), dimension(:), allocatable :: a1
    real(kind=dp), dimension(:), optional, allocatable :: a2, a3, a4, a5, &
         a6, a7, a8, a9, a10
    if (allocated(a1)) deallocate(a1)
    if (.not. present(a2)) go to 1
    if (allocated(a2)) deallocate(a2)
    if (.not. present(a3)) go to 1
    if (allocated(a3)) deallocate(a3)
    if (.not. present(a4)) go to 1
    if (allocated(a4)) deallocate(a4)
    if (.not. present(a5)) go to 1
    if (allocated(a5)) deallocate(a5)
    if (.not. present(a6)) go to 1
    if (allocated(a6)) deallocate(a6)
    if (.not. present(a7)) go to 1
    if (allocated(a7)) deallocate(a7)
    if (.not. present(a8)) go to 1
    if (allocated(a8)) deallocate(a8)
    if (.not. present(a9)) go to 1
    if (allocated(a9)) deallocate(a9)
    if (.not. present(a10)) go to 1
    if (allocated(a10)) deallocate(a10)
1   return
  end subroutine d_maybe_deallocate1

  subroutine c_maybe_deallocate1(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10)
    complex(kind=dp), dimension(:), allocatable :: a1
    complex(kind=dp), dimension(:), optional, allocatable :: a2, a3, a4, a5, &
         a6, a7, a8, a9, a10
    if (allocated(a1)) deallocate(a1)
    if (.not. present(a2)) go to 1
    if (allocated(a2)) deallocate(a2)
    if (.not. present(a3)) go to 1
    if (allocated(a3)) deallocate(a3)
    if (.not. present(a4)) go to 1
    if (allocated(a4)) deallocate(a4)
    if (.not. present(a5)) go to 1
    if (allocated(a5)) deallocate(a5)
    if (.not. present(a6)) go to 1
    if (allocated(a6)) deallocate(a6)
    if (.not. present(a7)) go to 1
    if (allocated(a7)) deallocate(a7)
    if (.not. present(a8)) go to 1
    if (allocated(a8)) deallocate(a8)
    if (.not. present(a9)) go to 1
    if (allocated(a9)) deallocate(a9)
    if (.not. present(a10)) go to 1
    if (allocated(a10)) deallocate(a10)
1   return
  end subroutine c_maybe_deallocate1

  subroutine i_maybe_deallocate1(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10)
    integer(kind=int32), dimension(:), allocatable :: a1
    integer(kind=int32), dimension(:), optional, allocatable :: a2, a3, a4, a5, &
         a6, a7, a8, a9, a10
    if (allocated(a1)) deallocate(a1)
    if (.not. present(a2)) go to 1
    if (allocated(a2)) deallocate(a2)
    if (.not. present(a3)) go to 1
    if (allocated(a3)) deallocate(a3)
    if (.not. present(a4)) go to 1
    if (allocated(a4)) deallocate(a4)
    if (.not. present(a5)) go to 1
    if (allocated(a5)) deallocate(a5)
    if (.not. present(a6)) go to 1
    if (allocated(a6)) deallocate(a6)
    if (.not. present(a7)) go to 1
    if (allocated(a7)) deallocate(a7)
    if (.not. present(a8)) go to 1
    if (allocated(a8)) deallocate(a8)
    if (.not. present(a9)) go to 1
    if (allocated(a9)) deallocate(a9)
    if (.not. present(a10)) go to 1
    if (allocated(a10)) deallocate(a10)
1   return
  end subroutine i_maybe_deallocate1

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

  function d_random_matrix(m,n) result(a)
    real(kind=dp), dimension(m,n) :: a
    integer(kind=int32), intent(in) :: m,n
    call d_random_matrix_to(a)
  end function d_random_matrix

  function d_random_vector(n) result(a)
    real(kind=dp), dimension(n) :: a
    integer(kind=int32), intent(in) :: n
    call d_v_random_matrix_to(a)
  end function d_random_vector

  function c_random_matrix(m,n) result(a)
    complex(kind=dp), dimension(m,n) :: a
    integer(kind=int32), intent(in) :: m,n
    call c_random_matrix_to(a)
  end function c_random_matrix

  function c_random_vector(n) result(a)
    complex(kind=dp), dimension(n) :: a
    integer(kind=int32), intent(in) :: n
    call c_v_random_matrix_to(a)
  end function c_random_vector
  
  subroutine d_random_matrix_to(a)
    real(kind=dp), dimension(:,:), intent(out) :: a
    call random_number(a)
    a=d_random_scale*a+d_random_shift
  end subroutine d_random_matrix_to

  subroutine d_v_random_matrix_to(a)
    real(kind=dp), dimension(:), intent(out) :: a
    call random_number(a)
    a=d_random_scale*a+d_random_shift
  end subroutine d_v_random_matrix_to

  subroutine d_s_random_matrix_to(a)
    real(kind=dp), intent(out) :: a
    call random_number(a)
    a=d_random_scale*a+d_random_shift
  end subroutine d_s_random_matrix_to

  subroutine c_random_matrix_to(a)
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
  end subroutine c_random_matrix_to

  subroutine c_v_random_matrix_to(a)
    complex(kind=dp), dimension(:), intent(out) :: a
    real(kind=dp) :: x,y
    integer(kind=int32) :: j
    do j=1,size(a)
       call random_number(x)
       call random_number(y)
       a(j)= c_random_scale*cmplx(x,y)+c_random_shift
    end do
  end subroutine c_v_random_matrix_to

  subroutine c_s_random_matrix_to(a)
    complex(kind=dp), intent(out) :: a
    real(kind=dp) :: x,y
    call random_number(x)
    call random_number(y)
    a=c_random_scale*cmplx(x,y)+c_random_shift
  end subroutine c_s_random_matrix_to

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

  integer(kind=int32) function i_equals_option(j,k) result(i)
    integer(kind=int32) :: j
    integer(kind=int32), optional :: k

    if(present(k)) then
       i=k
    else
       i=j
    end if
  end function i_equals_option

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
