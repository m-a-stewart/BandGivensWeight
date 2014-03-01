module conversions_bv_to_ub
  use misc
  use transforms
  use types
  implicit none

  interface convert_bv_to_ub
     module procedure d_convert_bv_to_ub, c_convert_bv_to_ub
  end interface convert_bv_to_ub

  interface f_convert_bv_to_ub
     module procedure f_d_convert_bv_to_ub, f_c_convert_bv_to_ub
  end interface f_convert_bv_to_ub

  type(routine_info), parameter :: info_d_convert_bv_to_ub=routine_info(id_d_convert_bv_to_ub, &
       'd_convert_bv_to_ub', &
       [ character(len=error_message_length) :: 'n<1', 'Insufficient storage in bv.', &
       'Insufficient Storage in ub.', 'ub%n /= bv%n' ] )

  type(routine_info), parameter :: info_c_convert_bv_to_ub=routine_info(id_c_convert_bv_to_ub, &
       'c_convert_bv_to_ub', &
       [ character(len=error_message_length) :: 'n<1', 'Insufficient storage in bv.', &
       'Insufficient Storage in ub.', 'ub%n /= bv%n' ] )

contains

  ! Errors:
  ! 0: no error
  ! 1: n<1
  ! 2: Insufficient storage in bv
  ! 3: Insufficient stroage in ub
  ! 4: ub%n /= bv%n
  subroutine d_convert_bv_to_ub(bv, ub, error)
    type(d_bv) :: bv
    type(d_ub) :: ub
    type(error_info), intent(out) :: error
    call clear_error(error)
    if (get_n(ub) < 1) then
       call set_error(error, 1, id_d_convert_bv_to_ub); return
    end if
    ! must allow for temporary fill-in
    if (get_lbwmax(bv)+get_ubwmax(bv)+1<bv%ubw+bv%lbw+2) then
       call set_error(error, 2, id_d_convert_bv_to_ub); return
    end if
    if (get_lbwmax(ub) < bv%lbw .or. get_ubwmax(ub) < bv%ubw) then
       call set_error(error, 3, id_d_convert_bv_to_ub); return
    end if
    if (get_n(ub) /= get_n(bv)) then
       call set_error(error, 4, id_d_convert_bv_to_ub); return
    end if

    call f_d_convert_bv_to_ub(bv%b, get_n(bv), bv%lbw, bv%ubw, get_lbwmax(bv), &
         get_ubwmax(bv), bv%numrotsv, bv%ksv, bv%csv, bv%ssv, ub%b,  ub%lbw, ub%ubw, get_lbwmax(ub), &
         get_ubwmax(ub), ub%numrotsu, ub%jsu, ub%csu, ub%ssu, error)
  end subroutine d_convert_bv_to_ub

  subroutine f_d_convert_bv_to_ub(b_bv, n, lbw, ubw, lbwmax_bv, ubwmax_bv, &
       numrots_bv, ks_bv, cs_bv, ss_bv, b_ub, lbw_ub, ubw_ub, lbwmax_ub, ubwmax_ub, numrots_ub, js_ub, &
       cs_ub, ss_ub, error)
    real(kind=dp), dimension(n,lbwmax_bv+ubwmax_bv+1), intent(inout) :: b_bv
    integer(kind=int32), intent(in) :: n, lbw, ubw, lbwmax_ub, ubwmax_ub, lbwmax_bv, ubwmax_bv
    integer(kind=int32), dimension(n), intent(in) :: numrots_bv
    integer(kind=int32), dimension(n,ubwmax_bv), intent(in) :: ks_bv
    real(kind=dp), dimension(n,ubwmax_bv), intent(in) :: cs_bv, ss_bv

    real(kind=dp), dimension(lbwmax_ub+ubwmax_ub+1,n), intent(out) :: b_ub
    integer(kind=int32), dimension(n), intent(out) :: numrots_ub
    integer(kind=int32), dimension(ubwmax_ub,n), intent(out) :: js_ub
    integer(kind=int32), intent(out) :: lbw_ub, ubw_ub
    real(kind=dp), dimension(ubwmax_ub,n), intent(out) :: cs_ub, ss_ub
    type(error_info), intent(out) :: error

    integer(kind=int32) :: j, k, ubw1, lbw1, d, k0, k1
    type(d_rotation) :: rot

    call clear_error(error)
    b_ub=0.0_dp; numrots_ub=0
    ss_ub=0.0_dp; cs_ub=0.0_dp
    js_ub=0
    lbw_ub=lbw; ubw_ub=ubw
    ubw1=ubw+1; lbw1=lbw
    if (n == 1) then
       b_ub(1,1)=b_bv(1,1)
       return
    end if
    b_bv(:,ubw1+lbw1+1)=0.0_dp
    do k=1,n-1
       do j=1,numrots_bv(n-k)
          rot%cosine=cs_bv(n-k,j); rot%sine=ss_bv(n-k,j)
          call tbr_times_rotation(b_bv,n,lbw1,ubw1,0,n-k,trp_rot(rot),ks_bv(n-k,j))
       end do
       !
       k0=max(k+2,ubw1+1)
       k1=min(k+ubw1,n)
       numrots_ub(k+1)=k1-k0+1
       do j=k1,k0,-1
          rot=lgivens2(get_el_br(b_bv,lbw1,j-ubw1,j), get_el_br(b_bv,lbw1,j-ubw1+1,j))
          js_ub(j-k0+1,k+1)=j-ubw1
          cs_ub(j-k0+1,k+1)=rot%cosine; ss_ub(j-k0+1,k+1)=rot%sine
          call rotation_times_tbr(trp_rot(rot),b_bv,n,lbw1,ubw1,k+1,0,j-ubw1)
       end do
    end do
    do d=1,lbw+1
       do j=1,n-lbw+d-1
          b_ub(ubw+lbw+2-d,j)=b_bv(lbw-d+j+1,d)
       end do
    end do
    do d=lbw+2,lbw+ubw+1
       do j=1,n-d+lbw+1
          b_ub(ubw+lbw+2-d,j+d-lbw-1)=b_bv(j,d)
       end do
    end do
  end subroutine f_d_convert_bv_to_ub

  subroutine c_convert_bv_to_ub(bv, ub, error)
    type(c_bv) :: bv
    type(c_ub) :: ub
    type(error_info), intent(out) :: error
    call clear_error(error)
    if (get_n(ub) < 1) then
       call set_error(error, 1, id_d_convert_bv_to_ub); return
    end if
    ! must allow for temporary fill-in
    if (get_lbwmax(bv)+get_ubwmax(bv)+1<bv%ubw+bv%lbw+2) then
       call set_error(error, 2, id_d_convert_bv_to_ub); return
    end if
    if (get_lbwmax(ub) < bv%lbw .or. get_ubwmax(ub) < bv%ubw) then
       call set_error(error, 3, id_d_convert_bv_to_ub); return
    end if
    if (get_n(ub) /= get_n(bv)) then
       call set_error(error, 4, id_d_convert_bv_to_ub); return
    end if
    call f_c_convert_bv_to_ub(bv%b, get_n(bv), bv%lbw, bv%ubw, get_lbwmax(bv), &
         get_ubwmax(bv), bv%numrotsv, bv%ksv, &
         bv%csv, bv%ssv, ub%b,  ub%lbw, ub%ubw, get_lbwmax(ub), get_ubwmax(ub), ub%numrotsu, ub%jsu, &
         ub%csu, ub%ssu, error)
  end subroutine c_convert_bv_to_ub

  subroutine f_c_convert_bv_to_ub(b_bv, n, lbw, ubw, lbwmax_bv, ubwmax_bv, numrots_bv, &
       ks_bv, cs_bv, ss_bv, b_ub, lbw_ub, ubw_ub, lbwmax_ub, ubwmax_ub, numrots_ub, js_ub, cs_ub, ss_ub, error)
    complex(kind=dp), dimension(n,lbwmax_bv+ubwmax_bv+1), intent(inout) :: b_bv
    integer(kind=int32), intent(in) :: n, lbw, ubw, lbwmax_bv, ubwmax_bv, lbwmax_ub, ubwmax_ub
    integer(kind=int32), dimension(n), intent(in) :: numrots_bv
    integer(kind=int32), dimension(n,ubwmax_bv), intent(in) :: ks_bv
    complex(kind=dp), dimension(n,ubwmax_bv), intent(in) :: cs_bv, ss_bv

    complex(kind=dp), dimension(lbwmax_ub+ubwmax_ub+1,n), intent(out) :: b_ub
    integer(kind=int32), dimension(n), intent(out) :: numrots_ub
    integer(kind=int32), dimension(ubwmax_ub,n), intent(out) :: js_ub
    integer(kind=int32), intent(out) :: lbw_ub, ubw_ub
    complex(kind=dp), dimension(ubwmax_ub,n), intent(out) :: cs_ub, ss_ub
    type(error_info), intent(out) :: error

    integer(kind=int32) :: j, k, ubw1, lbw1, d, k0, k1
    type(c_rotation) :: rot

    call clear_error(error)
    b_ub=(0.0_dp, 0.0_dp); numrots_ub=0
    ss_ub=(0.0_dp, 0.0_dp); cs_ub=(0.0_dp, 0.0_dp)
    js_ub=0
    lbw_ub=lbw; ubw_ub=ubw
    ubw1=ubw+1; lbw1=lbw

    if (n == 1) then
       b_ub(1,1)=b_bv(1,1)
       return
    end if
    b_bv(:,ubw1+lbw1+1)=(0.0_dp, 0.0_dp)
    ! apply v_{n-1}
    do k=1,n-1
       do j=1,numrots_bv(n-k)
          rot%cosine=cs_bv(n-k,j); rot%sine=ss_bv(n-k,j)
          call tbr_times_rotation(b_bv,n,lbw1,ubw1,0,n-k,trp_rot(rot),ks_bv(n-k,j))
       end do
       !
       k0=max(k+2,ubw1+1)
       k1=min(k+ubw1,n)
       numrots_ub(k+1)=k1-k0+1
       do j=k1,k0,-1
          rot=lgivens2(get_el_br(b_bv,lbw1,j-ubw1,j), get_el_br(b_bv,lbw1,j-ubw1+1,j))
          js_ub(j-k0+1,k+1)=j-ubw1
          cs_ub(j-k0+1,k+1)=rot%cosine; ss_ub(j-k0+1,k+1)=rot%sine
          call rotation_times_tbr(trp_rot(rot),b_bv,n,lbw1,ubw1,k+1,0,j-ubw1)
       end do
    end do
    ! store the results in b_ub
    do d=1,lbw+1
       do j=1,n-lbw+d-1
          b_ub(ubw+lbw+2-d,j)=b_bv(lbw-d+j+1,d)
       end do
    end do
    do d=lbw+2,lbw+ubw+1
       do j=1,n-d+lbw+1
          b_ub(ubw+lbw+2-d,j+d-lbw-1)=b_bv(j,d)
       end do
    end do
  end subroutine f_c_convert_bv_to_ub

end module conversions_bv_to_ub
