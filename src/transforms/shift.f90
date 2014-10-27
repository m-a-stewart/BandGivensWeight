module shift
  use prec

  interface shift2
     module procedure d_shift2, c_shift2
  end interface shift2

contains

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

  subroutine c_shift2(a,dj,dk)
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
  end subroutine c_shift2

end module shift
