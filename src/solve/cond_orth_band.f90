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

  implicit none
 
  private

  ! public :: ub_min_sv, d_ub_min_sv, z_ub_min_sv, &
  !      d_ub_max_sv, z_ub_max_sv
  
  integer(kind=int32), parameter :: default_maxit=20
  real(kind=dp), parameter :: default_tol=4*eps, default_tolres=4*eps

contains

  ! minimum singular values/vectors for upper triangular stored as a
  ! UB decomp.
  real(kind=dp) function d_ub_min_sv(ub,un,vn,res,tolres0, &
       maxit0,error) result(sigmau)
    type(d_ub), intent(in) :: ub
    real(kind=dp), dimension(:), intent(out) :: un,vn,res
    real(kind=dp), intent(in), optional :: tolres0
    type(error_info), intent(inout), optional :: error
    integer(kind=int32), intent(in), optional :: maxit0

    integer(kind=int32) :: n, k, maxit, p, q
    real(kind=dp) :: tolres, sigmav, maxub
    type(routine_info), parameter :: info=info_d_ub_min_sv

    sigmau=0.0_dp
    n=get_n(ub)
    if (failure(error)) return
    call push_id(info, error)
    if (n<1) then
       call set_error(1, info, error); return
    else if (n /= size(un) .or. n /= size(vn) .or. n /= size(res)) then
       call set_error(2, info, error); return
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
    end if
  end function d_ub_min_sv
  
end module mod_cond_orth_band
