module mod_utility
  use mod_prec
  implicit none

  private

  public :: maybe_deallocate, not_all_allocated

  public :: maxabs, z_maxabs, d_maxabs, z_maxabs_v, d_maxabs_v, &
       norm2, z_norm2_v, d_norm2_v, &
       normf, z_normf, d_normf

  public :: print_matrix, d_print_matrix, z_print_matrix, i_print_matrix

  public :: random_matrix_to, d_random_matrix_to, d_v_random_matrix_to, d_s_random_matrix_to, &
       z_random_matrix_to, z_v_random_matrix_to, z_s_random_matrix_to, &
       d_random_matrix, d_random_vector, z_random_matrix, z_random_vector

  public :: ip_transpose, d_ip_transpose, z_ip_transpose

  public :: equals_option, i_equals_option, d_equals_option, z_equals_option

  public :: dot_no_conjg

  public :: first_zero_diagonal, z_first_zero_diagonal, d_first_zero_diagonal, &
       reverse_first_zero_diagonal, z_reverse_first_zero_diagonal, d_reverse_first_zero_diagonal

  public :: first_zero, z_first_zero, d_first_zero, &
       reverse_first_zero, z_reverse_first_zero, d_reverse_first_zero
  
  interface maybe_deallocate
     module procedure d_maybe_deallocate2, z_maybe_deallocate2, i_maybe_deallocate2, &
          d_maybe_deallocate1, z_maybe_deallocate1, i_maybe_deallocate1
  end interface maybe_deallocate

  interface not_all_allocated
     module procedure d_not_all_allocated2, z_not_all_allocated2, &
          i_not_all_allocated2, &
          d_not_all_allocated1, z_not_all_allocated1, i_not_all_allocated1
  end interface not_all_allocated

  interface maxabs
     module procedure z_maxabs, d_maxabs, z_maxabs_v, d_maxabs_v
  end interface maxabs

  interface norm2
     module procedure z_norm2_v, d_norm2_v
  end interface norm2

  interface normf
     module procedure z_normf, d_normf
  end interface normf

  interface print_matrix
     module procedure d_print_matrix, z_print_matrix, i_print_matrix
  end interface print_matrix

  interface random_matrix_to
     module procedure d_random_matrix_to, d_v_random_matrix_to, d_s_random_matrix_to, &
          z_random_matrix_to, z_v_random_matrix_to, z_s_random_matrix_to
  end interface random_matrix_to

  interface ip_transpose
     module procedure d_ip_transpose, z_ip_transpose
  end interface ip_transpose

  interface equals_option
     module procedure i_equals_option, d_equals_option, z_equals_option
  end interface equals_option
  
  interface first_zero_diagonal
     module procedure z_first_zero_diagonal, d_first_zero_diagonal
  end interface first_zero_diagonal

  interface reverse_first_zero_diagonal
     module procedure z_reverse_first_zero_diagonal, d_reverse_first_zero_diagonal
  end interface reverse_first_zero_diagonal

  interface first_zero
     module procedure z_first_zero, d_first_zero
  end interface first_zero

  interface reverse_first_zero
     module procedure z_reverse_first_zero, d_reverse_first_zero
  end interface reverse_first_zero

  
  real(kind=dp), parameter :: d_random_shift=-1.0_dp, d_random_scale=2.0_dp
  complex(kind=dp), parameter :: z_random_shift=(-1.0_dp,-1.0_dp), z_random_scale=(2.0_dp,0.0_dp)

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

  logical function z_not_all_allocated2(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10) result(x)
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
  end function z_not_all_allocated2

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

  logical function z_not_all_allocated1(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10) result(x)
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
  end function z_not_all_allocated1

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

  subroutine z_maybe_deallocate2(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10)
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
  end subroutine z_maybe_deallocate2

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

  subroutine z_maybe_deallocate1(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10)
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
  end subroutine z_maybe_deallocate1

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

  real(kind=dp) function z_maxabs(a) result(x)
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
  end function z_maxabs

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

  real(kind=dp) function z_maxabs_v(a) result(x)
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
  end function z_maxabs_v

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

  real(kind=dp) function z_norm2_v(a) result(x)
    complex(kind=dp), dimension(:), intent(in) :: a
    !
    integer(kind=int32) :: j, m
    real(kind=dp) :: y
    m=size(a)
    y=z_maxabs_v(a)
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
  end function z_norm2_v

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

  real(kind=dp) function z_normf(a) result(x)
    complex(kind=dp), dimension(:,:), intent(in) :: a
    !
    integer(kind=int32) :: j, k, m, n
    real(kind=dp) :: y, tmp
    m=size(a,1)
    n=size(a,2)
    y=z_maxabs(a)
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
  end function z_normf

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

  subroutine z_print_matrix(a)
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
  end subroutine z_print_matrix

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

  function z_random_matrix(m,n) result(a)
    complex(kind=dp), dimension(m,n) :: a
    integer(kind=int32), intent(in) :: m,n
    call z_random_matrix_to(a)
  end function z_random_matrix

  function z_random_vector(n) result(a)
    complex(kind=dp), dimension(n) :: a
    integer(kind=int32), intent(in) :: n
    call z_v_random_matrix_to(a)
  end function z_random_vector
  
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

  subroutine z_random_matrix_to(a)
    complex(kind=dp), dimension(:,:), intent(out) :: a
    real(kind=dp) :: x,y
    integer(kind=int32) :: j, k
    do j=1,size(a,1)
       do k=1,size(a,2)
          call random_number(x)
          call random_number(y)
          a(j,k)=z_random_scale*cmplx(x,y)+z_random_shift
       end do
    end do
  end subroutine z_random_matrix_to

  subroutine z_v_random_matrix_to(a)
    complex(kind=dp), dimension(:), intent(out) :: a
    real(kind=dp) :: x,y
    integer(kind=int32) :: j
    do j=1,size(a)
       call random_number(x)
       call random_number(y)
       a(j)= z_random_scale*cmplx(x,y)+z_random_shift
    end do
  end subroutine z_v_random_matrix_to

  subroutine z_s_random_matrix_to(a)
    complex(kind=dp), intent(out) :: a
    real(kind=dp) :: x,y
    call random_number(x)
    call random_number(y)
    a=z_random_scale*cmplx(x,y)+z_random_shift
  end subroutine z_s_random_matrix_to

  integer(kind=int32) function i_equals_option(j,k) result(i)
    integer(kind=int32) :: j
    integer(kind=int32), optional :: k

    if(present(k)) then
       i=k
    else
       i=j
    end if
  end function i_equals_option

  real(kind=dp) function d_equals_option(x,y) result(z)
    real(kind=dp) :: x
    real(kind=dp), optional :: y

    if(present(y)) then
       z=y
    else
       z=x
    end if
  end function d_equals_option

  complex(kind=dp) function z_equals_option(x,y) result(z)
    complex(kind=dp) :: x
    complex(kind=dp), optional :: y

    if(present(y)) then
       z=y
    else
       z=x
    end if
  end function z_equals_option
  
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

  subroutine z_ip_transpose(a)
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
  end subroutine z_ip_transpose

  function dot_no_conjg(x,y) result(xy)
    complex(kind=dp), dimension(:), intent(in) :: x, y
    complex(kind=dp) :: xy
    integer(kind=int32) :: j
    xy=(0.0_dp,0.0_dp)
    do j=1,size(x)
       xy=xy+x(j)*y(j)
    end do
  end function dot_no_conjg

  ! finding zeros in a vector

  integer(kind=int32) function d_first_zero(a,tol) result(d)
    real(kind=dp), dimension(:), intent(in) :: a
    real(kind=dp), intent(in) :: tol
    !
    integer(kind=int32) :: j
    d=1
    do j=1,size(a)
       if (abs(a(j)) <= tol) return
       d=d+1
    end do
  end function d_first_zero

  integer(kind=int32) function z_first_zero(a,tol) result(d)
    complex(kind=dp), dimension(:), intent(in) :: a
    real(kind=dp), intent(in) :: tol
    !
    integer(kind=int32) :: j
    d=1
    do j=1,size(a)
       if (abs(a(j)) <= tol) return
       d=d+1
    end do
  end function z_first_zero
  
  integer(kind=int32) function d_reverse_first_zero(a,tol) result(d)
    real(kind=dp), dimension(:), intent(in) :: a
    real(kind=dp), intent(in) :: tol
    !
    integer(kind=int32) :: j
    d=size(a)
    do j=d,1,-1
       if (abs(a(j)) <= tol) return
       d=d-1
    end do
  end function d_reverse_first_zero

  integer(kind=int32) function z_reverse_first_zero(a,tol) result(d)
    complex(kind=dp), dimension(:), intent(in) :: a
    real(kind=dp), intent(in) :: tol
    !
    integer(kind=int32) :: j
    d=size(a)
    do j=d,1,-1
       if (abs(a(j)) <= tol) return
       d=d-1
    end do
  end function z_reverse_first_zero
  

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

end module mod_utility
