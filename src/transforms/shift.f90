module mod_shift
  use mod_prec
  implicit none
  ! Simple routines for applying shifts to vectors and matrices.
  private
  public :: shift, d_shift2, z_shift2, d_shift1, z_shift1, i_shift1

  interface shift
     module procedure d_shift2, z_shift2, d_shift1, z_shift1, i_shift1
  end interface shift

contains

  subroutine d_shift1(a,d)
    real(kind=dp), dimension(:), intent(inout) :: a
    integer(kind=int32), intent(in) :: d
    integer(kind=int32) :: n,j

    n=size(a)
    if (d > 0) then
       do j=n,d+1,-1
          a(j)=a(j-d)
       end do
       a(1:d)=0.0_dp
    else if (d < 0) then
       do j=1,n+d
          a(j)=a(j-d)
       end do
       a(n+d+1:n)=0.0_dp
    end if
  end subroutine d_shift1

  subroutine z_shift1(a,d)
    complex(kind=dp), dimension(:), intent(inout) :: a
    integer(kind=int32), intent(in) :: d
    integer(kind=int32) :: n,j

    n=size(a)
    if (d > 0) then
       do j=n,d+1,-1
          a(j)=a(j-d)
       end do
       a(1:d)=(0.0_dp,0.0_dp)
    else if (d < 0) then
       do j=1,n+d
          a(j)=a(j-d)
       end do
       a(n+d+1:n)=(0.0_dp,0.0_dp)
    end if
  end subroutine z_shift1

  subroutine i_shift1(a,d)
    integer(kind=int32), dimension(:), intent(inout) :: a
    integer(kind=int32), intent(in) :: d
    integer(kind=int32) :: n,j

    n=size(a)
    if (d > 0) then
       do j=n,d+1,-1
          a(j)=a(j-d)
       end do
       a(1:d)=0
    else if (d < 0) then
       do j=1,n+d
          a(j)=a(j-d)
       end do
       a(n+d+1:n)=0
    end if
  end subroutine i_shift1

  subroutine d_shift2(a,dj,dk)
    real(kind=dp), dimension(:,:), intent(inout) :: a
    integer(kind=int32), intent(in) :: dj, dk
    integer(kind=int32) :: m,n,j,k

    m=size(a,1)
    n=size(a,2)
    if (dj >= 0 .and. dk >= 0) then
       do j=m,dj+1,-1
          do k=n,dk+1,-1
             a(j,k)=a(j-dj,k-dk)
          end do
       end do
       a(1:dj,:)=0.0_dp
       a(:,1:dk)=0.0_dp
    else if (dj >=0 .and. dk < 0) then
       do j=m,dj+1,-1
          do k=1,n+dk
             a(j,k)=a(j-dj,k-dk)
          end do
       end do
       a(1:dj,:)=0.0_dp
       a(:,n+dk+1:n)=0.0_dp
    else if (dj < 0 .and. dk >= 0) then
       do j=1,m+dj
          do k=n,dk+1,-1
             a(j,k)=a(j-dj,k-dk)
          end do
       end do
       a(m+dj+1:m,:)=0.0_dp
       a(:,1:dk)=0.0_dp
    else
       do j=1,m+dj
          do k=1,n+dk
             a(j,k)=a(j-dj,k-dk)
          end do
       end do
       a(m+dj+1:m,:)=0.0_dp
       a(:,n+dk+1:n)=0.0_dp
    end if
  end subroutine d_shift2

  subroutine z_shift2(a,dj,dk)
    complex(kind=dp), dimension(:,:), intent(inout) :: a
    integer(kind=int32), intent(in) :: dj, dk
    integer(kind=int32) :: m,n,j,k

    m=size(a,1)
    n=size(a,2)
    if (dj >= 0 .and. dk >= 0) then
       do j=m,dj+1,-1
          do k=n,dk+1,-1
             a(j,k)=a(j-dj,k-dk)
          end do
       end do
       a(1:dj,:)=(0.0_dp,0.0_dp)
       a(:,1:dk)=(0.0_dp,0.0_dp)
    else if (dj >=0 .and. dk < 0) then
       do j=m,dj+1,-1
          do k=1,n+dk
             a(j,k)=a(j-dj,k-dk)
          end do
       end do
       a(1:dj,:)=(0.0_dp,0.0_dp)
       a(:,n+dk+1:n)=(0.0_dp,0.0_dp)
    else if (dj < 0 .and. dk >= 0) then
       do j=1,m+dj
          do k=n,dk+1,-1
             a(j,k)=a(j-dj,k-dk)
          end do
       end do
       a(m+dj+1:m,:)=(0.0_dp,0.0_dp)
       a(:,1:dk)=(0.0_dp,0.0_dp)
    else
       do j=1,m+dj
          do k=1,n+dk
             a(j,k)=a(j-dj,k-dk)
          end do
       end do
       a(m+dj+1:m,:)=(0.0_dp,0.0_dp)
       a(:,n+dk+1:n)=(0.0_dp,0.0_dp)
    end if
  end subroutine z_shift2

end module mod_shift
