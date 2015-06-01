module mod_cond
  use mod_triangular
  use mod_prec
  use mod_error_id
  use mod_utility
  implicit none
 
  private

  public :: lower_left_nullvec, z_lower_left_nullvec, d_lower_left_nullvec, &
       lower_right_nullvec, z_lower_right_nullvec, d_lower_right_nullvec, &
       lower_min_sv, d_lower_min_sv, z_lower_min_sv

  interface lower_left_nullvec
     module procedure z_lower_left_nullvec, d_lower_left_nullvec
  end interface lower_left_nullvec

  interface lower_right_nullvec
     module procedure z_lower_right_nullvec, d_lower_right_nullvec
  end interface lower_right_nullvec

  interface lower_min_sv
     module procedure  d_lower_min_sv, z_lower_min_sv
  end interface lower_min_sv

  integer(kind=int32), parameter :: default_maxit=20
  real(kind=dp), parameter :: default_tol=4*eps, default_tolres=4*eps

contains

  real(kind=dp) function d_lower_min_sv(l,un,vn,res,tolres0, &
       maxit0,error) result(sigmau)
    real(kind=dp), dimension(:,:), intent(in) :: l
    real(kind=dp), dimension(:), intent(out) :: un,vn,res
    real(kind=dp), intent(in), optional :: tolres0
    type(error_info), intent(inout), optional :: error
    integer(kind=int32), intent(in), optional :: maxit0

    integer(kind=int32) :: n, k, maxit, p, q
    real(kind=dp) :: tolres, sigmav, maxl
    type(routine_info), parameter :: info=info_d_lower_min_sv
    !
    if (failure(error)) return
    call push_id(info, error)

    n=size(l,1)
    maxl=maxabs(l)
    maxit=equals_option(default_maxit,maxit0)
    tolres=equals_option(maxl*default_tolres,tolres0)

    un=0.0_dp
    vn=0.0_dp
    if (maxl==0.0_dp) then
       un(1)=1.0_dp; vn(1)=1.0_dp
       sigmau=0.0_dp
       res=0.0_dp
       call pop_id(error)
       return
    else
       ! treat very small diagonal elements as hard zeros.
       p=first_zero_diagonal(l, maxl*eps*eps)
       q=reverse_first_zero_diagonal(l, maxl*eps*eps)
       if (p<=n .or. q >= 1) then
          un(p)=1.0_dp
          vn(q)=1.0_dp
          if (p > 1) call lower_right_invert(un(1:p-1),l(1:p-1,1:p-1),-l(p,1:p-1))
          if (q < n) call lower_left_invert(l(q+1:n, q+1:n), vn(q+1:n), -l(q+1:n,q))
          un=un/norm2(un); vn=vn/norm2(vn)
          sigmau=0.0_dp
          ! Residual is L*vn-sigma*un or just L*vn
          res=0.0_dp
          res(q)=l(q,q)*vn(q)
          call lower_left_multiply(l(q+1:n,q+1:n),vn(q+1:n),res(q+1:n))
          res(q+1:n)=res(q+1:n)+l(q+1:n,q)*vn(q)
          call pop_id(error); return
       else
          ! Get initial vectors.
          call lower_right_invert_linpack(un,l,vn)
          sigmau=1/norm2(un)
          un=sigmau*un
          sigmau=sigmau*norm2(vn)
          k=1
          do while (k < maxit)
             ! Update v
             call lower_left_invert(l,vn,un)
             sigmav=1/norm2(vn)
             vn=sigmav*vn
             ! Update un
             res=un
             call lower_tr_left_invert(l,un,vn)
             sigmau=1/norm2(un)
             un=sigmau*un
             res=res*sigmav - un*sigmau
             if (norm2(res) < tolres) then
                call pop_id(error)
                return
             end if
             k=k+1
          end do
          if (tolres==0.0_dp) then
             ! maxits is OK if tolerance is zero.
             call pop_id(error)
          else
             call set_error(1, info, error); return
          end if
       end if
    end if
  end function d_lower_min_sv


  real(kind=dp) function z_lower_min_sv(l,un,vn,res,tolres0, &
       maxit0,error) result(sigmau)
    complex(kind=dp), dimension(:,:), intent(in) :: l
    complex(kind=dp), dimension(:), intent(out) :: un,vn,res
    real(kind=dp), intent(in), optional :: tolres0
    type(error_info), intent(inout), optional :: error
    integer(kind=int32), intent(in), optional :: maxit0

    integer(kind=int32) :: n, k, maxit, p, q
    real(kind=dp) :: tolres, sigmav, maxl
    type(routine_info), parameter :: info=info_z_lower_min_sv
    !
    if (failure(error)) return
    call push_id(info, error)

    n=size(l,1)
    maxl=maxabs(l)
    maxit=equals_option(default_maxit,maxit0)
    tolres=equals_option(maxl*default_tolres,tolres0)

    un=(0.0_dp,0.0_dp)
    vn=(0.0_dp,0.0_dp)
    if (maxl==0.0_dp) then
       un(1)=(1.0_dp,0.0_dp); vn(1)=(1.0_dp,0.0_dp)
       sigmau=0.0_dp
       res=(0.0_dp,0.0_dp)
       call pop_id(error)
       return
    else
       ! treat very small diagonal elements as hard zeros.
       p=first_zero_diagonal(l, maxl*eps*eps)
       q=reverse_first_zero_diagonal(l, maxl*eps*eps)
       if (p<=n .or. q >= 1) then
          un(p)=(1.0_dp,0.0_dp)
          vn(q)=(1.0_dp,0.0_dp)
          if (p > 1) call lower_right_invert(un(1:p-1),l(1:p-1,1:p-1),-l(p,1:p-1))
          if (q < n) call lower_left_invert(l(q+1:n, q+1:n), vn(q+1:n), -l(q+1:n,q))
          un=conjg(un)
          un=un/norm2(un); vn=vn/norm2(vn)
          sigmau=0.0_dp
          ! Residual is L*v-sigma*u or just L*v
          res=(0.0_dp,0.0_dp)
          res(q)=l(q,q)*vn(q)
          call lower_left_multiply(l(q+1:n,q+1:n),vn(q+1:n),res(q+1:n))
          res(q+1:n)=res(q+1:n)+l(q+1:n,q)*vn(q)
          call pop_id(error); return
       else
          ! Get initial vectors.
          call lower_right_invert_linpack(un,l,vn)
          sigmau=1/norm2(un)
          un=sigmau*un
          sigmau=sigmau*norm2(vn)
          k=1
          do while (k < maxit)
             ! Update vn
             call lower_left_invert(l,vn,un)
             sigmav=1/norm2(vn)
             vn=sigmav*vn
             ! Update un
             res=un
             call lower_tr_left_invert(l,un,vn)
             sigmau=1/norm2(un)
             un=sigmau*un
             res=res*sigmav - un*sigmau
             if (norm2(res) < tolres) then
                call pop_id(error)
                return
             end if
             k=k+1
          end do
          if (tolres==0.0_dp) then
             ! maxits shouldn't raise an error if the tolerance is zero.
             call pop_id(error)
          else
             call set_error(1, info, error); return
          end if
       end if
    end if
  end function z_lower_min_sv
   
  ! errors:
  ! 1: A is not square.
  ! 2: failure to converge in computing norm2(A)
  ! 3: failure to converge in computing norm2(inverse(A))
  ! real(kind=dp) function d_cond2_upper(a,tol,maxit,error) result(k2)
  !   real(kind=dp), dimension(:,:), intent(in) :: a
  !   real(kind=dp), intent(in) :: tol
  !   integer(kind=int32), intent(in) :: maxit
  !   type(error_info), intent(inout), optional :: error
    
  !   integer(kind=int32) :: n, j, k
  !   real(kind=dp) :: maxa, nrma, nrmainv
  !   real(kind=dp) :: d, tmp
  !   real(kind=dp), dimension(size(a,1)) :: y, x
  !   type(routine_info), parameter :: info=info_d_cond2_upper

  !   if (failure(error)) return
  !   call push_id(info, error)
  !   n=size(a,1)
  !   if (size(a,2)/=n) then
  !      call set_error(1, info, error); return       
  !   end if
    
  !   x=0.0_dp
  !   y=0.0_dp
  !   maxa=maxabs(a)

  !   k2=0.0_dp
  ! end function d_cond2_upper

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

    integer(kind=int32) :: n, k
    real(kind=dp) :: nrmx
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
    integer(kind=int32) :: n, k
    real(kind=dp) :: nrmx
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
    integer(kind=int32) :: n, k
    real(kind=dp) :: nrmx
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
    integer(kind=int32) :: n, k
    real(kind=dp) :: nrmx
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

end module mod_cond