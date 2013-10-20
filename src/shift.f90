module shift
use prec

interface upper_left_shift
   module procedure d_upper_left_shift
end interface upper_left_shift

interface right_shift
   module procedure d_right_shift
end interface right_shift

interface down_shift
   module procedure d_down_shift
end interface down_shift

interface left_shift
   module procedure d_left_shift
end interface left_shift

contains

  subroutine d_upper_left_shift(a)
    real(kind=dp), dimension(:,:), intent(inout) :: a
    integer(kind=int32) :: m,n,j,k
    m=size(a,1)
    n=size(a,2)
    do j=1,m-1
       do k=1,n-1
          a(j,k) = a(j+1,k+1)
       end do
    end do
    a(:,n)=0.0_dp
    a(m,:)=0.0_dp
  end subroutine d_upper_left_shift

  subroutine d_left_shift(a)
    real(kind=dp), dimension(:,:), intent(inout) :: a
    integer(kind=int32) :: m,n,j,k
    m=size(a,1)
    n=size(a,2)
    do j=1,m
       do k=1,n-1
          a(j,k) = a(j,k+1)
       end do
    end do
    a(:,n)=0.0_dp
  end subroutine d_left_shift

  subroutine d_right_shift(a)
    real(kind=dp), dimension(:,:), intent(inout) :: a
    integer(kind=int32) :: m,n,j,k
    m=size(a,1)
    n=size(a,2)
    do j=1,m
       do k=n,2,-1
          a(j,k) = a(j,k-1)
       end do
    end do
    a(:,1)=0.0_dp
  end subroutine d_right_shift

  subroutine d_down_shift(a)
    real(kind=dp), dimension(:,:), intent(inout) :: a
    integer(kind=int32) :: m,n,j,k
    m=size(a,1)
    n=size(a,2)
    do j=m,2,-1
       do k=1,n
          a(j,k) = a(j-1,k)
       end do
    end do
    a(1,:)=0.0_dp
  end subroutine d_down_shift

end module shift
