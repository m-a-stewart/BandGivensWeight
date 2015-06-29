module mod_cond_orth_band
  use mod_prec
  use mod_error_id
  use mod_rotation
  use mod_orth_band_types
  use mod_band_types
  use mod_solve
  use mod_products
  use mod_utility
  use mod_submatrix
  use mod_convert

  implicit none
 
  private

  public :: ub_min_sv, d_ub_min_sv, z_ub_min_sv, &
       ub_max_sv, d_ub_max_sv, z_ub_max_sv

  interface ub_min_sv
     module procedure d_ub_min_sv, z_ub_min_sv
  end interface ub_min_sv

  interface ub_max_sv
     module procedure d_ub_max_sv, z_ub_max_sv
  end interface ub_max_sv
  
  integer(kind=int32), parameter :: default_maxit=20
  real(kind=dp), parameter :: default_tol=4*eps, default_tolres=4*eps

contains

  ! minimum singular values/vectors for upper triangular stored as a
  ! UB decomp.
  ! errors:
  ! 1: n<1
  ! 2: size error
  ! 3: ub is not upper triangular
  ! 4: numerically singular; no null vector computed
  ! 5: failed to converge
  real(kind=dp) function d_ub_min_sv(ub,un,vn,res,tolres0, &
       maxit0,error) result(sigmau)
    type(d_ub), intent(in) :: ub
    real(kind=dp), dimension(:), intent(out) :: un,vn,res
    real(kind=dp), intent(in), optional :: tolres0
    type(error_info), intent(inout), optional :: error
    integer(kind=int32), intent(in), optional :: maxit0

    integer(kind=int32) :: n, k, maxit, p, q
    real(kind=dp) :: tolres, sigmav, maxub
    real(kind=dp), dimension(size(un)) :: tmpv
    real(kind=dp), dimension(size(un),1) :: tmpc
    real(kind=dp), dimension(1,size(un)) :: tmpr    
    type(routine_info), parameter :: info=info_d_ub_min_sv
    type(d_bv), allocatable :: bv, bvt
    type(d_ub), allocatable :: ubl

    sigmau=0.0_dp
    n=get_n(ub)
    if (failure(error)) return
    call push_id(info, error)
    if (n<1) then
       call set_error(1, info, error); return
    else if (n /= size(un) .or. n /= size(vn) .or. n /= size(res)) then
       call set_error(2, info, error); return
    else if (ub%lbw /=0) then
       call set_error(3, info, error); return
    end if
 
    maxub=maxabs(ub%bc(1:ub%lbw+ub%ubw+1,:))
    maxit=equals_option(default_maxit,maxit0)
    tolres=equals_option(maxub*default_tolres,tolres0)

    un=0.0_dp
    vn=0.0_dp
    if (maxub==0.0_dp) then
       un(1)=1.0_dp; vn(1)=1.0_dp
       sigmau=0.0_dp
       res=0.0_dp
       call pop_id(error)
       return
    else
       ! treat very small diagonal elements as hard zeros.
       p=first_zero(ub%bc(ub%ubw+1,:), maxub*eps*eps)
       q=reverse_first_zero(ub%bc(ub%ubw+1,:), maxub*eps*eps)
       bv=bv_of(ub,error)
       if (p<=n .or. q >= 1) then
          vn(p)=1.0_dp
          un(q)=1.0_dp
          if (p > 1) then
             ubl=leading(ub,p-1)
             call ub_to_columns(p,p,ub,tmpc,error)
             tmpv(1:p-1)=-tmpc(1:p-1,1)
             call back_solve_ub(ubl,vn(1:p-1),tmpv(1:p-1),error)
          end if
          if (q < n) then
             bvt=trailing(bv,n-q)
             call bv_to_rows(q,q,bv,tmpr,error)
             tmpv(q+1:n)=-tmpr(1,q+1:n)
             call forward_solve_bv(un(q+1:n),bvt,tmpv(q+1:n))
          end if
          if (failure(error)) return
          if (.not. is_number(un(q+1:n)) .or. .not. is_number(vn(1:p-1))) then
             call set_error(4,info,error)
             return
          end if
          ! un=conjg(un)
          un=un/norm2(un); vn=vn/norm2(vn)
          sigmau=0.0_dp
          ! Residual is R*vn-sigma*un or just R*vn
          res=0.0_dp
          ubl=leading(ub,p)
          tmpc(1:p,1)=vn(1:p)
          call ub_times_general(ubl,reshape(vn(1:p),[p,1]),tmpc(1:p,:),error)          
          res(1:p)=tmpc(1:p,1)
          if (failure(error)) return          
          call pop_id(error); return
       else
          ! start with a random vector
          call random_matrix_to(un)
          un=un/norm2(un)
          k=1
          do while (k <= maxit)
             ! Update v
             tmpv=un
             call back_solve_ub(ub,vn,tmpv,error)
             if (failure(error)) return
             if (.not. is_number(vn)) then
                call set_error(4,info,error)
                return
             end if
             sigmav=1/norm2(vn)
             vn=sigmav*vn
             ! Update un
             res=un
             tmpv=vn
             call forward_solve_bv(un,bv,tmpv,error)
             if (failure(error)) return
             if (.not. is_number(un)) then
                call set_error(4,info,error)
                return
             end if
             ! un=conjg(un) and conjg(vn) above
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
             call set_error(5, info, error); return
          end if
       end if
    end if
  end function d_ub_min_sv

  ! minimum singular values/vectors for upper triangular stored as a
  ! UB decomp.
  real(kind=dp) function z_ub_min_sv(ub,un,vn,res,tolres0, &
       maxit0,error) result(sigmau)
    type(z_ub), intent(in) :: ub
    complex(kind=dp), dimension(:), intent(out) :: un,vn,res
    real(kind=dp), intent(in), optional :: tolres0
    type(error_info), intent(inout), optional :: error
    integer(kind=int32), intent(in), optional :: maxit0

    integer(kind=int32) :: n, k, maxit, p, q
    real(kind=dp) :: tolres, sigmav, maxub
    complex(kind=dp), dimension(size(un)) :: tmpv
    complex(kind=dp), dimension(size(un),1) :: tmpc
    complex(kind=dp), dimension(1,size(un)) :: tmpr    
    type(routine_info), parameter :: info=info_z_ub_min_sv
    type(z_bv), allocatable :: bv, bvt
    type(z_ub), allocatable :: ubl
    

    sigmau=0.0_dp
    n=get_n(ub)
    if (failure(error)) return
    call push_id(info, error)
    if (n<1) then
       call set_error(1, info, error); return
    else if (n /= size(un) .or. n /= size(vn) .or. n /= size(res)) then
       call set_error(2, info, error); return
    else if (ub%lbw /=0) then
       call set_error(3, info, error); return
    end if

    maxub=maxabs(ub%bc(1:ub%lbw+ub%ubw+1,:))
    maxit=equals_option(default_maxit,maxit0)
    tolres=equals_option(maxub*default_tolres,tolres0)

    un=(0.0_dp,0.0_dp)
    vn=(0.0_dp,0.0_dp)
    if (maxub==0.0_dp) then
       un(1)=(1.0_dp,0.0_dp); vn(1)=(1.0_dp,0.0_dp)
       sigmau=0.0_dp
       res=(0.0_dp,0.0_dp)
       call pop_id(error)
       return
    else
       ! treat very small diagonal elements as hard zeros.
       p=first_zero(ub%bc(ub%ubw+1,:), maxub*eps*eps)
       q=reverse_first_zero(ub%bc(ub%ubw+1,:), maxub*eps*eps)
       bv=bv_of(ub,error)
       if (p<=n .or. q >= 1) then
          vn(p)=(1.0_dp,0.0_dp)
          un(q)=(1.0_dp,0.0_dp)
          if (p > 1) then
             ubl=leading(ub,p-1)
             call ub_to_columns(p,p,ub,tmpc,error)
             tmpv(1:p-1)=-tmpc(1:p-1,1)
             call back_solve_ub(ubl,vn(1:p-1),tmpv(1:p-1),error)
          end if
          if (q < n) then
             bvt=trailing(bv,n-q)
             call bv_to_rows(q,q,bv,tmpr,error)
             tmpv(q+1:n)=-tmpr(1,q+1:n)
             call forward_solve_bv(un(q+1:n),bvt,tmpv(q+1:n))
          end if
          if (failure(error)) return
          if (.not. is_number(un(q+1:n)) .or. .not. is_number(vn(1:p-1))) then
             call set_error(4,info,error)
             return
          end if
          un=conjg(un)
          un=un/norm2(un); vn=vn/norm2(vn)
          sigmau=0.0_dp
          ! Residual is R*vn-sigma*un or just R*vn
          res=(0.0_dp,0.0_dp)
          ubl=leading(ub,p)
          tmpc(1:p,1)=vn(1:p)
          call ub_times_general(ubl,reshape(vn(1:p),[p,1]),tmpc(1:p,:),error)          
          res(1:p)=tmpc(1:p,1)
          if (failure(error)) return          
          call pop_id(error); return
       else
          ! start with a random vector
          call random_matrix_to(un)
          un=un/norm2(un)
          k=1
          do while (k <= maxit)
             ! Update v
             tmpv=un
             call back_solve_ub(ub,vn,tmpv,error)
             if (failure(error)) return
             if (.not. is_number(vn)) then
                call set_error(4,info,error)
                return
             end if
             sigmav=1/norm2(vn)
             vn=sigmav*vn
             ! Update un
             res=un
             tmpv=conjg(vn)
             call forward_solve_bv(un,bv,tmpv,error)
             if (failure(error)) return
             if (.not. is_number(un)) then
                call set_error(4,info,error)
                return
             end if
             un=conjg(un)
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
             call set_error(4, info, error); return
          end if
       end if
    end if
  end function z_ub_min_sv
  
  ! maximum singular values/vectors for upper triangular stored as a
  ! UB decomp.
  ! errors:
  ! 1: n<1
  ! 2: size error
  ! 3: failed to converge
  real(kind=dp) function d_ub_max_sv(ub,u1,v1,res,tolres0, &
       maxit0,error) result(sigmau)
    type(d_ub), intent(in) :: ub
    real(kind=dp), dimension(:), intent(out) :: u1,v1,res
    real(kind=dp), intent(in), optional :: tolres0
    type(error_info), intent(inout), optional :: error
    integer(kind=int32), intent(in), optional :: maxit0

    integer(kind=int32) :: n, k, maxit
    real(kind=dp) :: tolres, sigmav, maxub, sigmau0
    real(kind=dp), dimension(size(u1),1) :: tmpc
    real(kind=dp), dimension(1,size(u1)) :: tmpr    
    type(routine_info), parameter :: info=info_d_ub_max_sv

    sigmau=0.0_dp
    n=get_n(ub)
    if (failure(error)) return
    call push_id(info, error)
    if (n<1) then
       call set_error(1, info, error); return
    else if (n /= size(u1) .or. n /= size(v1) .or. n /= size(res)) then
       call set_error(2, info, error); return
    end if
 
    maxub=maxabs(ub%bc(1:ub%lbw+ub%ubw+1,:))
    maxit=equals_option(default_maxit,maxit0)
    tolres=equals_option(maxub*default_tolres,tolres0)

    u1=0.0_dp
    v1=0.0_dp
    if (maxub==0.0_dp) then
       u1(1)=1.0_dp; v1(1)=1.0_dp
       sigmau=0.0_dp
       res=0.0_dp
       call pop_id(error)
       return
    else
       ! start with a random vector
       call random_matrix_to(v1)
       v1=v1/norm2(v1)
       k=1
       do while (k <= maxit)
          ! Update u1
          res=u1; sigmau0=sigmau
          call ub_times_general(ub,reshape(v1,[n,1]),tmpc,error)
          u1=tmpc(:,1)
          if (failure(error)) return
          sigmav=norm2(u1)
          u1=u1/sigmav
          ! update v1
          call general_times_ub(reshape(u1,[1,n]),ub,tmpr)
          v1=tmpr(1,:)
          sigmau=norm2(v1)
          v1=v1/sigmau
          res = u1*sigmav-res*sigmau0
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
          call set_error(3, info, error); return
       end if
    end if
  end function d_ub_max_sv


  ! maximum singular values/vectors for upper triangular stored as a
  ! UB decomp.
  ! errors:
  ! 1: n<1
  ! 2: size error
  ! 3: failed to converge
  real(kind=dp) function z_ub_max_sv(ub,u1,v1,res,tolres0, &
       maxit0,error) result(sigmau)
    type(z_ub), intent(in) :: ub
    complex(kind=dp), dimension(:), intent(out) :: u1,v1,res
    real(kind=dp), intent(in), optional :: tolres0
    type(error_info), intent(inout), optional :: error
    integer(kind=int32), intent(in), optional :: maxit0

    integer(kind=int32) :: n, k, maxit
    real(kind=dp) :: tolres, sigmav, maxub, sigmau0
    complex(kind=dp), dimension(size(u1),1) :: tmpc
    complex(kind=dp), dimension(1,size(u1)) :: tmpr    
    type(routine_info), parameter :: info=info_z_ub_max_sv

    sigmau=0.0_dp
    n=get_n(ub)
    if (failure(error)) return
    call push_id(info, error)
    if (n<1) then
       call set_error(1, info, error); return
    else if (n /= size(u1) .or. n /= size(v1) .or. n /= size(res)) then
       call set_error(2, info, error); return
    end if
 
    maxub=maxabs(ub%bc(1:ub%lbw+ub%ubw+1,:))
    maxit=equals_option(default_maxit,maxit0)
    tolres=equals_option(maxub*default_tolres,tolres0)

    u1=(0.0_dp,0.0_dp)
    v1=(0.0_dp,0.0_dp)
    if (maxub==0.0_dp) then
       u1(1)=1.0_dp; v1(1)=1.0_dp
       sigmau=0.0_dp
       res=(0.0_dp,0.0_dp)
       call pop_id(error)
       return
    else
       ! start with a random vector
       call random_matrix_to(v1)
       v1=v1/norm2(v1)
       k=1
       do while (k <= maxit)
          ! Update u1
          res=u1; sigmau0=sigmau
          call ub_times_general(ub,reshape(v1,[n,1]),tmpc,error)
          u1=tmpc(:,1)
          if (failure(error)) return
          sigmav=norm2(u1)
          u1=u1/sigmav
          ! update v1
          call general_times_ub(reshape(conjg(u1),[1,n]),ub,tmpr)
          v1=conjg(tmpr(1,:))
          sigmau=norm2(v1)
          v1=v1/sigmau
          res = u1*sigmav-res*sigmau0
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
          call set_error(3, info, error); return
       end if
    end if
  end function z_ub_max_sv
  
end module mod_cond_orth_band
