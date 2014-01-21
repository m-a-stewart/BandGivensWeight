module test_data
use prec
implicit none

integer, parameter :: n=11, rmax=3, ubwmax=rmax+1, lbw=1, lbwmax=3
real(kind=dp), parameter :: tol=1e-14, tol1=1e-14, tol2=1e-10

contains
  
  subroutine d_assemble_a(a,u,v,d,lbw)
    real(kind=dp), dimension(:,:), intent(out) :: a
    real(kind=dp), dimension(:,:), intent(in) :: u, v
    real(kind=dp), dimension(:), intent(in) :: d
    integer(kind=int32), intent(in) :: lbw

    integer(kind=int32) :: j,k,n
    n=size(a,1)
    a=matmul(u,v)
    do k=1,n-lbw-1
       do j=k+lbw+1,n
          a(j,k)=0.0_dp
       end do
    end do
    do j=1,n
       a(j,j)=d(j)
    end do
  end subroutine d_assemble_a

  subroutine c_assemble_a(a,u,v,d,lbw)
    complex(kind=dp), dimension(:,:), intent(out) :: a
    complex(kind=dp), dimension(:,:), intent(in) :: u, v
    complex(kind=dp), dimension(:), intent(in) :: d
    integer(kind=int32), intent(in) :: lbw

    integer(kind=int32) :: j,k,n
    n=size(a,1)
    a=matmul(u,v)
    do k=1,n-lbw-1
       do j=k+lbw+1,n
          a(j,k)=(0.0_dp, 0.0_dp)
       end do
    end do
    do j=1,n
       a(j,j)=d(j)
    end do
  end subroutine c_assemble_a


end module test_data
