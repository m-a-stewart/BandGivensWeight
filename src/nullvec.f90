module nullvec
use prec
use utility
use triangular
implicit none

interface lower_left_nullvec
   module procedure f_c_lower_left_nullvec, f_d_lower_left_nullvec
end interface lower_left_nullvec

interface lower_right_nullvec
   module procedure f_c_lower_right_nullvec, f_d_lower_right_nullvec
end interface lower_right_nullvec

contains

  ! null vector of a lower triangular matrix.
  ! error: 0 no error
  !        1 no null vector within tolerance
  subroutine f_d_lower_left_nullvec(x,l,tol,maxit, error)
    real(kind=dp), dimension(:,:), intent(in) :: l
    real(kind=dp), dimension(:), intent(out) :: x
    real(kind=dp), intent(in) :: tol
    integer(kind=int32), intent(out) :: error
    integer(kind=int32), intent(in) :: maxit
    !
    integer(kind=int32) :: n, j, k, p
    real(kind=dp) :: nrmx
    real(kind=dp) :: d, tmp
    real(kind=dp), dimension(size(x)) :: y
    !
    n=size(l,1)
    x=0.0_dp
    y=0.0_dp
    p=first_zero_diagonal(l, 1.0e-1_dp*tol)
    if (p==1) then
       x(1)=1.0_dp
       error=1
    else if (p <= n) then
       ! if a small enough diagonal element is found, solve
       ! for a null vector directly.
       x(p)=1.0_dp
       call lower_right_invert(x(1:p-1),l(1:p-1, 1:p-1), -l(p,1:p-1))
       error = p
    else
       ! Otherwise use the Linpack condition estimator approach.
       x(n)=1.0_dp/l(n,n)
       do k=n-1,1,-1
          tmp = 0.0_dp
          do j=k+1,n
             tmp = tmp + x(j)*l(j,k)
          end do
          if (tmp==0) then
             d=1.0_dp
          else
             d=-tmp/abs(tmp)
          end if
          x(k) = (d-tmp)/l(k,k)
       end do
       x=x/maxabs(x)
       ! solve y^T L^T = x^T
       call lower_tr_right_invert(y,l,x)
       y=y/norm2(y)
       call lower_right_invert(x,l,y)
       nrmx=norm2(x)
       x=x/nrmx
       k=1
       error=-1
       do while (k < maxit)
          if (1/nrmx < tol) then
             error=0
             return
          end if
          call lower_tr_right_invert(y,l,x)
          y=y/norm2(y)
          call lower_right_invert(x,l,y)
          nrmx=norm2(x)
          x=x/nrmx
          k=k+1
       end do
    end if
  end subroutine f_d_lower_left_nullvec

  ! null vector of a lower triangular matrix.
  subroutine f_c_lower_left_nullvec(x,l,tol,maxit,error)
    complex(kind=dp), dimension(:,:), intent(in) :: l
    complex(kind=dp), dimension(:), intent(out) :: x
    real(kind=dp), intent(in) :: tol
    integer(kind=int32), intent(out) :: error
    integer(kind=int32), intent(in) :: maxit
    !
    integer(kind=int32) :: n, j, k, p
    real(kind=dp) :: nrmx
    complex(kind=dp) :: d, tmp
    complex(kind=dp), dimension(size(x)) :: y
    !
    n=size(l,1)
    x=(0.0_dp, 0.0_dp)
    y=(0.0_dp, 0.0_dp)
    p=first_zero_diagonal(l, 1.0e-1_dp*tol)
    if (p==1) then
       x(p)=(1.0_dp, 0.0_dp)
       error=1
    else if (p <= n) then
       ! if a small enough diagonal element is found, solve
       ! for a null vector directly.
       x(p)=(1.0_dp, 0.0_dp)
       call lower_right_invert(x(1:p-1),l(1:p-1, 1:p-1), -l(p,1:p-1))
       error=p
    else
       ! Otherwise use the Linpack condition estimator approach.
       x(n)=(1.0_dp,0.0_dp)/l(n,n)
       do k=n-1,1,-1
          tmp = (0.0_dp, 0.0_dp)
          do j=k+1,n
             tmp = tmp + x(j)*l(j,k)
          end do
          if (tmp==0) then
             d=(1.0_dp,0.0_dp)
          else
             d=-tmp/abs(tmp)
          end if
          x(k) = (d-tmp)/l(k,k)
       end do
       x=x/maxabs(x)
       ! solve y^T L^H = x^T
       call lower_tr_right_invert(y,l,x)
       y=y/norm2(y)
       call lower_right_invert(x,l,y)
       nrmx=norm2(x)
       x=x/nrmx
       k=1
       error=-1
       do while (k < maxit)
          if (1/nrmx < tol) then
             error=0
             return
          end if
          call lower_tr_right_invert(y,l,x)
          y=y/norm2(y)
          call lower_right_invert(x,l,y)
          nrmx=norm2(x)
          x=x/nrmx
          k=k+1
       end do
    end if
  end subroutine f_c_lower_left_nullvec


  ! right null vector of a lower triangular matrix.
  subroutine f_d_lower_right_nullvec(x,l,tol,maxit, error)
    real(kind=dp), dimension(:,:), intent(in) :: l
    real(kind=dp), dimension(:), intent(out) :: x
    real(kind=dp), intent(in) :: tol
    integer(kind=int32), intent(out) :: error
    integer(kind=int32), intent(in) :: maxit
    !
    integer(kind=int32) :: n, j, k, p
    real(kind=dp) :: nrmx
    real(kind=dp) :: d, tmp
    real(kind=dp), dimension(size(x)) :: y
    !
    n=size(l,1)
    x=0.0_dp
    y=0.0_dp
    p=last_zero_diagonal(l, 1.0e-1_dp*tol)
    if (p==n) then
       x(n)=1.0_dp
       error=n
    else if (p >= 1) then
       ! if a small enough diagonal element is found, solve
       ! for a null vector directly.
       x(p)=1.0_dp
       call lower_left_invert(l(p+1:n,p+1:n),x(p+1:n),-l(p+1:n,p))
       error = p
    else
       ! Otherwise use the Linpack condition estimator approach.
       x(1)=1.0_dp/l(1,1)
       do k=2,n
          tmp = 0.0_dp
          do j=1,k
             tmp=tmp+l(k,j)*x(j)
          end do
          if (tmp==0) then
             d=1.0_dp
          else
             d=-tmp/abs(tmp)
          end if
          x(k)=(d-tmp)/l(k,k)
       end do
       x=x/maxabs(x)
       ! solve L^T y = x
       call lower_tr_left_invert(l,y,x)
       y=y/norm2(y)
       ! and L x = y
       call lower_left_invert(l,x,y)
       nrmx=norm2(x)
       x=x/nrmx
       k=1
       error=-1
       do while (k < maxit)
          if (1/nrmx < tol) then
             error=0
             return
          end if
          call lower_tr_left_invert(l,y,x)
          y=y/norm2(y)
          call lower_left_invert(l,x,y)
          nrmx=norm2(x)
          x=x/nrmx
          k=k+1
       end do
    end if
  end subroutine f_d_lower_right_nullvec

  ! right null vector of a lower triangular matrix.
  subroutine f_c_lower_right_nullvec(x,l,tol,maxit, error)
    complex(kind=dp), dimension(:,:), intent(in) :: l
    complex(kind=dp), dimension(:), intent(out) :: x
    real(kind=dp), intent(in) :: tol
    integer(kind=int32), intent(out) :: error
    integer(kind=int32), intent(in) :: maxit
    !
    integer(kind=int32) :: n, j, k, p
    real(kind=dp) :: nrmx
    complex(kind=dp) :: d, tmp
    complex(kind=dp), dimension(size(x)) :: y
    !
    n=size(l,1)
    x=(0.0_dp, 0.0_dp)
    y=(0.0_dp, 0.0_dp)
    p=last_zero_diagonal(l, 1.0e-1_dp*tol)
    if (p==n) then
       x(n)=(1.0_dp, 0.0_dp)
       error=n
    else if (p >= 1) then
       ! if a small enough diagonal element is found, solve
       ! for a null vector directly.
       x(p)=(1.0_dp, 0.0_dp)
       call lower_left_invert(l(p+1:n,p+1:n),x(p+1:n),-l(p+1:n,p))
       error = p
    else
       ! Otherwise use the Linpack condition estimator approach.
       x(1)=(1.0_dp,0.0_dp)/l(1,1)
       do k=2,n
          tmp = (0.0_dp, 0.0_dp)
          do j=1,k
             tmp=tmp+l(k,j)*x(j)
          end do
          if (tmp==0) then
             d=(1.0_dp, 0.0_dp)
          else
             d=-tmp/abs(tmp)
          end if
          x(k)=(d-tmp)/l(k,k)
       end do
       x=x/maxabs(x)
       ! solve L^T y = x
       call lower_tr_left_invert(l,y,x)
       y=y/norm2(y)
       ! and L x = y
       call lower_left_invert(l,x,y)
       nrmx=norm2(x)
       x=x/nrmx
       k=1
       error=-1
       do while (k < maxit)
          if (1/nrmx < tol) then
             error=0
             return
          end if
          call lower_tr_left_invert(l,y,x)
          y=y/norm2(y)
          call lower_left_invert(l,x,y)
          nrmx=norm2(x)
          x=x/nrmx
          k=k+1
       end do
    end if
  end subroutine f_c_lower_right_nullvec

end module nullvec
