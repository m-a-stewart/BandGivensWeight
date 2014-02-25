module qr_factorization
use prec
use assemble
use shift
use rotation
use band_types
use nested_types
use nullvec
use compressions_ub_to_bv
use conversions_ub_to_bv

implicit none

contains

!
! Errors
! 0: no error
! 1: ub%n /= bv%n
! 2: bv%lbw <= 0
! 3: dimension of cs or ss /= n
!
subroutine d_reduce_lbw_bv_to_ub(bv,ub,cs,ss,error)
    type(d_ub) :: ub
    type(d_bv) :: bv
    type(error_info), intent(out) :: error
    real(kind=dp), dimension(:) :: cs, ss
    if (get_n(ub) /= get_n(bv)) then
       call set_error(error, 1, id_d_reduce_lbw_bv_to_ub); return
    end if
    if (bv%lbw <= 0) then
       call set_error(error, 2, id_d_reduce_lbw_bv_to_ub); return
    end if
    if (size(cs) /= get_n(ub) .or. size(ss) /= get_n(ub)) then
       call set_error(error, 3, id_d_reduce_lbw_bv_to_ub); return
    end if
    call f_d_reduce_lbw_bv_to_ub(bv%b, get_n(bv), bv%lbw, bv%ubw, get_lbwmax(bv), &
         get_ubwmax(bv), bv%numrotsv, bv%ksv, bv%csv, bv%ssv, & 
         ub%b, ub%lbw, ub%ubw, get_lbwmax(ub), get_ubwmax(ub), ub%numrotsu, ub%jsu, ub%csu, ub%ssu, &
         cs, ss, error)
  end subroutine d_reduce_lbw_bv_to_ub

  ! 1: n < 1
  ! 2: Insufficient temp storage in b_bv
  ! 3: Insufficient storage in b_ub

  subroutine f_d_reduce_lbw_bv_to_ub(b_bv, n, lbw_bv, ubw_bv, lbwmax_bv, ubwmax_bv, numrots_bv, &
       ks_bv, cs_bv, ss_bv, &
       b_ub, lbw_ub, ubw_ub, lbwmax_ub, ubwmax_ub, numrots_ub, js_ub, cs_ub, ss_ub, cs, ss, error)
    integer(kind=int32), intent(in) :: n, lbw_bv, ubw_bv, lbwmax_ub, ubwmax_ub, lbwmax_bv, ubwmax_bv
    real(kind=dp), dimension(n,lbwmax_bv+ubwmax_bv+1), intent(inout) :: b_bv
    integer(kind=int32), dimension(n), intent(in) :: numrots_bv
    integer(kind=int32), dimension(n,ubwmax_bv), intent(in) :: ks_bv
    real(kind=dp), dimension(n,ubwmax_bv), intent(in) :: cs_bv, ss_bv

    real(kind=dp), dimension(lbwmax_ub+ubwmax_ub+1,n), intent(out) :: b_ub
    integer(kind=int32), dimension(n), intent(out) :: numrots_ub
    integer(kind=int32), dimension(ubwmax_ub,n), intent(out) :: js_ub
    real(kind=dp), dimension(ubwmax_ub,n), intent(out) :: cs_ub, ss_ub

    real(kind=dp), dimension(n) :: cs, ss
    integer(kind=int32), intent(out) :: lbw_ub, ubw_ub
    type(error_info), intent(out) :: error

    integer(kind=int32) :: j, k, d, k0,k1, lbw, ubw
    type(d_rotation) :: rot
    
    call clear_error(error)
    numrots_ub=0
    ss_ub=0.0_dp; cs_ub=0.0_dp
    ss=0.0_dp; cs=1.0_dp
    js_ub=0
    b_ub=0.0_dp
    ubw=ubw_bv+2
    lbw=lbw_bv
    
    if (n < 1) then
       call set_error(error, 1, id_f_d_reduce_lbw_bv_to_ub); return
    end if
    if (n==1) then
       b_ub(1,1)=b_bv(1,1)
    end if

    if (lbwmax_bv+ubwmax_bv+1<ubw+lbw+1) then
       call set_error(error, 2, id_f_d_reduce_lbw_bv_to_ub); return
    end if
    if (lbwmax_ub < lbw_bv-1 .or. ubwmax_ub < ubw_bv+1) then
       call set_error(error, 3, id_f_d_reduce_lbw_bv_to_ub); return
    end if
    do k=1,n-1
       ! apply v_{n-k}
       do j=1,numrots_bv(n-k)
          rot%cosine=cs_bv(n-k,j); rot%sine=ss_bv(n-k,j)
          call tbr_times_rotation(b_bv,n,lbw,ubw,0,n-k,trp_rot(rot),ks_bv(n-k,j))
       end do
       ! eliminate subdiagonal in row k+1 and column k-lbw+1
       if (k-lbw+1 >= 1 .and. k+1 <= n) then
          rot=lgivens(get_el_br(b_bv,lbw,k,k-lbw+1),get_el_br(b_bv,lbw,k+1,k-lbw+1))
          call rotation_times_tbr(trp_rot(rot),b_bv,n,lbw,ubw,0,0,k)
          cs(k)=rot%cosine; ss(k)=rot%sine
       end if
       ! Eliminate superdiagonal ubw.
       ! columns that have a nonzero in superdiagonal ubw
       ! Three cases:
       ! 1. For k <= ubw-1, columns ubw+1, ..., k+ubw-1 might
       !    have a nonzero in superdiagonal ubw.
       ! 2. For ubw <= k <= n-ubw+1, columns k+1, .... k+ubw-1
       !    might have a nonzero in superdiagonal ubw.
       ! 3. For n-ubw+2 <= k <= n-1, columns k+1, ..., n
       !    might have a nonzero in superdiagonal ubw.
       k0=max(k+1,ubw+1)
       k1=min(k+ubw-1,n)
       numrots_ub(k)=k1-k0+1
       do j=k1,k0,-1
          rot=lgivens2(get_el_br(b_bv,lbw,j-ubw,j), get_el_br(b_bv,lbw,j-ubw+1,j))
          js_ub(j-k0+1,k)=j-ubw
          cs_ub(j-k0+1,k)=rot%cosine; ss_ub(j-k0+1,k)=rot%sine
          call rotation_times_tbr(trp_rot(rot),b_bv,n,lbw,ubw,k,0,j-ubw)
       end do
    end do
    ubw_ub=ubw_bv+1
    lbw_ub=lbw_bv-1
    ! put diagonals in b
    do d=1,ubw_ub+1
       do k=ubw_ub-d+2,n
          b_ub(d,k) = get_el_br(b_bv,lbw_bv,d+k-ubw_ub-1,k)
       end do
    end do
    do d=ubw_ub+2, ubw_ub+lbw_ub+1
       do k=1,n-d+ubw_ub+1
          b_ub(d,k) = get_el_br(b_bv,lbw_bv,d+k-ubw_ub-1,k)
       end do
    end do
  end subroutine f_d_reduce_lbw_bv_to_ub


  subroutine c_reduce_lbw_bv_to_ub(bv,ub,cs,ss,error)
    type(c_ub) :: ub
    type(c_bv) :: bv
    type(error_info), intent(out) :: error
    complex(kind=dp), dimension(:) :: cs, ss
    if (get_n(ub) /= get_n(bv)) then
       call set_error(error, 1, id_d_reduce_lbw_bv_to_ub); return
    end if
    if (bv%lbw <= 0) then
       call set_error(error, 2, id_d_reduce_lbw_bv_to_ub); return
    end if
    if (size(cs) /= get_n(ub) .or. size(ss) /= get_n(ub)) then
       call set_error(error, 3, id_d_reduce_lbw_bv_to_ub); return
    end if
    call f_c_reduce_lbw_bv_to_ub(bv%b, get_n(bv), bv%lbw, bv%ubw, get_lbwmax(bv), &
         get_ubwmax(bv), bv%numrotsv, bv%ksv, bv%csv, bv%ssv, & 
         ub%b, ub%lbw, ub%ubw, get_lbwmax(ub), get_ubwmax(ub), ub%numrotsu, ub%jsu, ub%csu, ub%ssu, &
         cs, ss, error)
  end subroutine c_reduce_lbw_bv_to_ub

  ! 1: n < 1
  ! 2: Insufficient temp storage in b_bv
  ! 3: Insufficient storage in b_ub

  subroutine f_c_reduce_lbw_bv_to_ub(b_bv, n, lbw_bv, ubw_bv, lbwmax_bv, ubwmax_bv, numrots_bv, &
       ks_bv, cs_bv, ss_bv, &
       b_ub, lbw_ub, ubw_ub, lbwmax_ub, ubwmax_ub, numrots_ub, js_ub, cs_ub, ss_ub, cs, ss, error)
    integer(kind=int32), intent(in) :: n, lbw_bv, ubw_bv, lbwmax_ub, ubwmax_ub, lbwmax_bv, ubwmax_bv
    complex(kind=dp), dimension(n,lbwmax_bv+ubwmax_bv+1), intent(inout) :: b_bv
    integer(kind=int32), dimension(n), intent(in) :: numrots_bv
    integer(kind=int32), dimension(n,ubwmax_bv), intent(in) :: ks_bv
    complex(kind=dp), dimension(n,ubwmax_bv), intent(in) :: cs_bv, ss_bv

    complex(kind=dp), dimension(lbwmax_ub+ubwmax_ub+1,n), intent(out) :: b_ub
    integer(kind=int32), dimension(n), intent(out) :: numrots_ub
    integer(kind=int32), dimension(ubwmax_ub,n), intent(out) :: js_ub
    complex(kind=dp), dimension(ubwmax_ub,n), intent(out) :: cs_ub, ss_ub

    complex(kind=dp), dimension(n) :: cs, ss
    integer(kind=int32), intent(out) :: lbw_ub, ubw_ub
    type(error_info), intent(out) :: error

    integer(kind=int32) :: j, k, d, k0,k1, lbw, ubw
    type(c_rotation) :: rot
    
    call clear_error(error)
    numrots_ub=0
    ss_ub=(0.0_dp, 0.0_dp); cs_ub=(0.0_dp, 0.0_dp)
    ss=(0.0_dp, 0.0_dp); cs=(1.0_dp, 0.0_dp)
    js_ub=0
    b_ub=(0.0_dp, 0.0_dp)
    ubw=ubw_bv+2
    lbw=lbw_bv
    
    if (n < 1) then
       call set_error(error, 1, id_f_d_reduce_lbw_bv_to_ub); return
    end if
    if (n==1) then
       b_ub(1,1)=b_bv(1,1)
    end if

    if (lbwmax_bv+ubwmax_bv+1<ubw+lbw+1) then
       call set_error(error, 2, id_f_d_reduce_lbw_bv_to_ub); return
    end if
    if (lbwmax_ub < lbw_bv-1 .or. ubwmax_ub < ubw_bv+1) then
       call set_error(error, 3, id_f_d_reduce_lbw_bv_to_ub); return
    end if
    do k=1,n-1
       ! apply v_{n-k}
       do j=1,numrots_bv(n-k)
          rot%cosine=cs_bv(n-k,j); rot%sine=ss_bv(n-k,j)
          call tbr_times_rotation(b_bv,n,lbw,ubw,0,n-k,trp_rot(rot),ks_bv(n-k,j))
       end do
       ! eliminate subdiagonal in row k+1 and column k-lbw+1
       if (k-lbw+1 >= 1 .and. k+1 <= n) then
          rot=lgivens(get_el_br(b_bv,lbw,k,k-lbw+1),get_el_br(b_bv,lbw,k+1,k-lbw+1))
          call rotation_times_tbr(trp_rot(rot),b_bv,n,lbw,ubw,0,0,k)
          cs(k)=rot%cosine; ss(k)=rot%sine
       end if
       ! Eliminate superdiagonal ubw.
       ! columns that have a nonzero in superdiagonal ubw
       ! Three cases:
       ! 1. For k <= ubw-1, columns ubw+1, ..., k+ubw-1 might
       !    have a nonzero in superdiagonal ubw.
       ! 2. For ubw <= k <= n-ubw+1, columns k+1, .... k+ubw-1
       !    might have a nonzero in superdiagonal ubw.
       ! 3. For n-ubw+2 <= k <= n-1, columns k+1, ..., n
       !    might have a nonzero in superdiagonal ubw.
       k0=max(k+1,ubw+1)
       k1=min(k+ubw-1,n)
       numrots_ub(k)=k1-k0+1
       do j=k1,k0,-1
          rot=lgivens2(get_el_br(b_bv,lbw,j-ubw,j), get_el_br(b_bv,lbw,j-ubw+1,j))
          js_ub(j-k0+1,k)=j-ubw
          cs_ub(j-k0+1,k)=rot%cosine; ss_ub(j-k0+1,k)=rot%sine
          call rotation_times_tbr(trp_rot(rot),b_bv,n,lbw,ubw,k,0,j-ubw)
       end do
    end do
    ubw_ub=ubw_bv+1
    lbw_ub=lbw_bv-1
    ! put diagonals in b
    do d=1,ubw_ub+1
       do k=ubw_ub-d+2,n
          b_ub(d,k) = get_el_br(b_bv,lbw_bv,d+k-ubw_ub-1,k)
       end do
    end do
    do d=ubw_ub+2, ubw_ub+lbw_ub+1
       do k=1,n-d+ubw_ub+1
          b_ub(d,k) = get_el_br(b_bv,lbw_bv,d+k-ubw_ub-1,k)
       end do
    end do
  end subroutine f_c_reduce_lbw_bv_to_ub


end module qr_factorization

