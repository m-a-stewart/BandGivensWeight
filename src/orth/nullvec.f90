module mod_nullvec
  use mod_triangular
  use mod_prec
  use mod_error_id
  use mod_utility
  implicit none

  private

  public :: lower_left_nullvec, z_lower_left_nullvec, d_lower_left_nullvec, &
       lower_right_nullvec, z_lower_right_nullvec, d_lower_right_nullvec

  interface lower_left_nullvec
     module procedure z_lower_left_nullvec, d_lower_left_nullvec
  end interface lower_left_nullvec

  interface lower_right_nullvec
     module procedure z_lower_right_nullvec, d_lower_right_nullvec
  end interface lower_right_nullvec

contains

  ! null vector of a lower triangular matrix.
  ! error: 1 no null vector within tolerance
  ! p returns a value such that x(p+1:n)=0
  ! If tol=0.0_dp, just do maxits iterations.
  subroutine d_lower_left_nullvec(x,l,tol,maxit, p, error)
    real(kind=dp), dimension(:,:), intent(in) :: l
    real(kind=dp), dimension(:), intent(out) :: x
    real(kind=dp), intent(in) :: tol
    type(error_info), intent(inout), optional :: error
    integer(kind=int32), intent(in) :: maxit
    integer(kind=int32), intent(out) :: p

    integer(kind=int32) :: n, j, k
    real(kind=dp) :: nrmx
    real(kind=dp) :: d, tmp
    real(kind=dp), dimension(size(x)) :: y
    type(routine_info), parameter :: info=info_d_lower_left_nullvec
    !
    if (failure(error)) return
    call push_id(info, error)
    
    n=size(l,1)
    x=0.0_dp
    y=0.0_dp
    p=first_zero_diagonal(l, tol)
    If (p==1) then
       x(1)=1.0_dp
       call pop_id(error); return
    else if (p <= n) then
       ! if a small enough diagonal element is found, solve
       ! for a null vector directly.
       x(p)=1.0_dp
       call lower_right_invert(x(1:p-1),l(1:p-1, 1:p-1), -l(p,1:p-1))
       call pop_id(error); return
    else
       ! Otherwise use the Linpack condition estimator approach.
       p=n
       call lower_right_invert_linpack(x,l,y)
       nrmx=norm2(x)
       x=x/nrmx
       nrmx=nrmx/norm2(y)
       k=1
       do while (k < maxit)
          if (1/nrmx < tol) then
             call pop_id(error)
             return
          end if
          call lower_tr_right_invert(y,l,x)
          y=y/norm2(y)
          call lower_right_invert(x,l,y)
          nrmx=norm2(x)
          x=x/nrmx
          k=k+1
       end do
       if (tol==0.0_dp) then
          call pop_id(error)
       else
          call set_error(1, info, error); return
       end if
    end if
  end subroutine d_lower_left_nullvec

  ! null vector of a lower triangular matrix.
  subroutine z_lower_left_nullvec(x,l,tol,maxit,p, error)
    complex(kind=dp), dimension(:,:), intent(in) :: l
    complex(kind=dp), dimension(:), intent(out) :: x
    real(kind=dp), intent(in) :: tol
    type(error_info), intent(inout), optional :: error
    integer(kind=int32), intent(in) :: maxit
    integer(kind=int32), intent(out) :: p
    !
    integer(kind=int32) :: n, j, k
    real(kind=dp) :: nrmx
    complex(kind=dp) :: d, tmp
    complex(kind=dp), dimension(size(x)) :: y
    type(routine_info), parameter :: info=info_z_lower_left_nullvec
    !
    n=size(l,1)
    if (failure(error)) return
    call push_id(info, error)

    x=(0.0_dp, 0.0_dp)
    y=(0.0_dp, 0.0_dp)
    p=first_zero_diagonal(l, tol)
    if (p==1) then
       x(p)=(1.0_dp, 0.0_dp)
       call pop_id(error); return
    else if (p <= n) then
       ! if a small enough diagonal element is found, solve
       ! for a null vector directly.
       x(p)=(1.0_dp, 0.0_dp)
       call lower_right_invert(x(1:p-1),l(1:p-1, 1:p-1), -l(p,1:p-1))
       call pop_id(error); return
    else
       ! Otherwise use the Linpack condition estimator.
       p=n
       call lower_right_invert_linpack(x,l,y)
       nrmx=norm2(x)
       x=x/nrmx
       nrmx=nrmx/norm2(y)
       k=1
       do while (k < maxit)
          if (1/nrmx < tol) then
             call pop_id(error)
             return
          end if
          call lower_tr_right_invert(y,l,x)
          y=y/norm2(y)
          call lower_right_invert(x,l,y)
          nrmx=norm2(x)
          x=x/nrmx
          k=k+1
       end do
       if (tol==0.0_dp) then
          call pop_id(error)
       else
          call set_error(1, info, error); return
       end if
    end if
  end subroutine z_lower_left_nullvec

  ! right null vector of a lower triangular matrix.
  ! Returns p such that x(1:p-1)=0
  subroutine d_lower_right_nullvec(x,l,tol,maxit, p, error)
    real(kind=dp), dimension(:,:), intent(in) :: l
    real(kind=dp), dimension(:), intent(out) :: x
    real(kind=dp), intent(in) :: tol
    type(error_info), intent(inout), optional :: error
    integer(kind=int32), intent(in) :: maxit
    integer(kind=int32), intent(out) :: p
    !
    integer(kind=int32) :: n, j, k
    real(kind=dp) :: nrmx
    real(kind=dp) :: d, tmp
    real(kind=dp), dimension(size(x)) :: y
    type(routine_info), parameter :: info=info_d_lower_right_nullvec
    !
    n=size(l,1)
    if (failure(error)) return
    call push_id(info, error)
    x=0.0_dp
    y=0.0_dp
    p=reverse_first_zero_diagonal(l, tol)
    if (p==n) then
       x(n)=1.0_dp
       call pop_id(error); return
    else if (p >= 1) then
       ! if a small enough diagonal element is found, solve
       ! for a null vector directly.
       x(p)=1.0_dp
       call lower_left_invert(l(p+1:n,p+1:n),x(p+1:n),-l(p+1:n,p))
       call pop_id(error); return
    else
       ! Otherwise use the Linpack condition estimator approach.
       p=1
       call lower_left_invert_linpack(l,x,y)
       nrmx=norm2(x)
       x=x/nrmx
       nrmx=nrmx/norm2(y)
       k=1
       do while (k < maxit)
          if (1/nrmx < tol) then
             call pop_id(error); return
          end if
          call lower_tr_left_invert(l,y,x)
          y=y/norm2(y)
          call lower_left_invert(l,x,y)
          nrmx=norm2(x)
          x=x/nrmx
          k=k+1
       end do
       if (tol==0.0_dp) then
          call pop_id(error)
       else
          call set_error(1, info, error); return
       end if
    end if
  end subroutine d_lower_right_nullvec

  ! right null vector of a lower triangular matrix.
  subroutine z_lower_right_nullvec(x, l, tol, maxit, p, error)
    complex(kind=dp), dimension(:,:), intent(in) :: l
    complex(kind=dp), dimension(:), intent(out) :: x
    real(kind=dp), intent(in) :: tol
    type(error_info), intent(inout), optional :: error
    integer(kind=int32), intent(in) :: maxit
    integer(kind=int32), intent(out) :: p
    !
    integer(kind=int32) :: n, j, k
    real(kind=dp) :: nrmx
    complex(kind=dp) :: d, tmp
    complex(kind=dp), dimension(size(x)) :: y
    type(routine_info), parameter :: info=info_z_lower_right_nullvec
    !
    n=size(l,1)
    if (failure(error)) return
    call push_id(info, error)
    x=(0.0_dp, 0.0_dp)
    y=(0.0_dp, 0.0_dp)
    p=reverse_first_zero_diagonal(l, tol)
    if (p==n) then
       x(n)=(1.0_dp, 0.0_dp)
       call pop_id(error); return
    else if (p >= 1) then
       ! if a small enough diagonal element is found, solve
       ! for a null vector directly.
       x(p)=(1.0_dp, 0.0_dp)
       call lower_left_invert(l(p+1:n,p+1:n),x(p+1:n),-l(p+1:n,p))
       call pop_id(error); return
    else
       ! Otherwise use the Linpack condition estimator approach.
       p=1
       call lower_left_invert_linpack(l,x,y)
       nrmx=norm2(x)
       x=x/nrmx
       nrmx=nrmx/norm2(y)
       k=1
       do while (k < maxit)
          if (1/nrmx < tol) then
             call pop_id(error)
             return
          end if
          call lower_tr_left_invert(l,y,x)
          y=y/norm2(y)
          call lower_left_invert(l,x,y)
          nrmx=norm2(x)
          x=x/nrmx
          k=k+1
       end do
       if (tol==0.0_dp) then
          call pop_id(error)
       else
          call set_error(1, info, error); return
       end if
    end if
  end subroutine z_lower_right_nullvec

end module mod_nullvec
