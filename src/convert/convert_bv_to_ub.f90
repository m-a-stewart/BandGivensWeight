module mod_convert_bv_to_ub
  use mod_prec
  use mod_error_id
  use mod_rotation
  use mod_nested_types
  use mod_band_types
  implicit none

  private
  public :: convert_bv_to_ub, d_convert_bv_to_ub, c_convert_bv_to_ub, &
       f_convert_bv_to_ub, f_d_convert_bv_to_ub, f_c_convert_bv_to_ub

  public :: info_d_convert_bv_to_ub, info_c_convert_bv_to_ub

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
    if (get_ubwmax(bv) < bv%ubw+1 .and. bv%ubw < get_n(bv)-1) then
       call set_error(error, 2, id_d_convert_bv_to_ub); return
    end if
    if (get_lbwmax(ub) < bv%lbw .or. get_ubwmax(ub) < bv%ubw) then
       call set_error(error, 3, id_d_convert_bv_to_ub); return
    end if
    if (get_n(ub) /= get_n(bv)) then
       call set_error(error, 4, id_d_convert_bv_to_ub); return
    end if

    call f_d_convert_bv_to_ub(bv%br, get_n(bv), bv%lbw, bv%ubw, get_lbwmax(bv), &
         get_ubwmax(bv), bv%numrotsv, bv%ksv, bv%csv, bv%ssv, ub%bc,  ub%lbw, ub%ubw, get_lbwmax(ub), &
         get_ubwmax(ub), ub%numrotsu, ub%jsu, ub%csu, ub%ssu, error)
  end subroutine d_convert_bv_to_ub

  subroutine f_d_convert_bv_to_ub(b_bv, n, lbw, ubw, lbwmax_bv, ubwmax_bv, &
       numrotsv, ksv, csv, ssv, b_ub, lbw_ub, ubw_ub, lbwmax_ub, ubwmax_ub, numrotsu, jsu, &
       csu, ssu, error)
    real(kind=dp), dimension(n,lbwmax_bv+ubwmax_bv+1), intent(inout) :: b_bv
    integer(kind=int32), intent(in) :: n, lbw, ubw, lbwmax_ub, ubwmax_ub, lbwmax_bv, ubwmax_bv
    integer(kind=int32), dimension(n), intent(in) :: numrotsv
    integer(kind=int32), dimension(n,ubwmax_bv), intent(in) :: ksv
    real(kind=dp), dimension(n,ubwmax_bv), intent(in) :: csv, ssv

    real(kind=dp), dimension(lbwmax_ub+ubwmax_ub+1,n), intent(out) :: b_ub
    integer(kind=int32), dimension(n), intent(out) :: numrotsu
    integer(kind=int32), dimension(ubwmax_ub,n), intent(out) :: jsu
    integer(kind=int32), intent(out) :: lbw_ub, ubw_ub
    real(kind=dp), dimension(ubwmax_ub,n), intent(out) :: csu, ssu
    type(error_info), intent(out) :: error

    integer(kind=int32) :: j, k, ubw1, lbw1, k0, k1
    type(d_rotation) :: rot

    call clear_error(error)
    b_ub(1:lbw+ubw+1,:)=0.0_dp; numrotsu=0
    ssu(1:ubw,:)=0.0_dp; csu(1:ubw,:)=0.0_dp
    jsu(1:ubw,:)=0
    lbw_ub=lbw; ubw_ub=ubw
    if (n == 1) then
       b_ub(1,1)=b_bv(1,1);
       lbw_ub=0; ubw_ub=0; numrotsu=0; return
    end if
    lbw1=lbw
    if (ubw < n-1) then
       ubw1=ubw+1
       b_bv(:,ubw1+lbw1+1)=0.0_dp
    else
       ubw1=ubw
    end if
    ! k is the size of the leading principal submatrix
    do k=1,n-2
       ! Apply V_k
       do j=1,numrotsv(k)
          rot%cosine=csv(k,j); rot%sine=ssv(k,j)
          call tbr_times_rotation(b_bv,n,lbw1,ubw1,0,n-k,trp_rot(rot),ksv(k,j))
       end do
       !
       k0=max(k+2,ubw1+1)
       k1=min(k+ubw1,n)
       numrotsu(k+1)=max(k1-k0+1,0)
       do j=k1,k0,-1
          rot=lgivens2(get_el_br(b_bv,lbw1,j-ubw1,j), get_el_br(b_bv,lbw1,j-ubw1+1,j))
          jsu(j-k0+1,k+1)=j-ubw1
          csu(j-k0+1,k+1)=rot%cosine; ssu(j-k0+1,k+1)=rot%sine
          call rotation_times_tbr(trp_rot(rot),b_bv,n,lbw1,ubw1,k+1,0,j-ubw1)
       end do
    end do
    call br_to_bc(b_bv,b_ub,lbw,ubw)
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
    if (get_ubwmax(bv) < bv%ubw+1 .and. bv%ubw < get_n(bv)-1) then
       call set_error(error, 2, id_d_convert_bv_to_ub); return
    end if
    if (get_lbwmax(ub) < bv%lbw .or. get_ubwmax(ub) < bv%ubw) then
       call set_error(error, 3, id_d_convert_bv_to_ub); return
    end if
    if (get_n(ub) /= get_n(bv)) then
       call set_error(error, 4, id_d_convert_bv_to_ub); return
    end if
    call f_c_convert_bv_to_ub(bv%br, get_n(bv), bv%lbw, bv%ubw, get_lbwmax(bv), &
         get_ubwmax(bv), bv%numrotsv, bv%ksv, &
         bv%csv, bv%ssv, ub%bc,  ub%lbw, ub%ubw, get_lbwmax(ub), get_ubwmax(ub), ub%numrotsu, ub%jsu, &
         ub%csu, ub%ssu, error)
  end subroutine c_convert_bv_to_ub

  subroutine f_c_convert_bv_to_ub(b_bv, n, lbw, ubw, lbwmax_bv, ubwmax_bv, numrotsv, &
       ksv, csv, ssv, b_ub, lbw_ub, ubw_ub, lbwmax_ub, ubwmax_ub, numrotsu, jsu, csu, ssu, error)
    complex(kind=dp), dimension(n,lbwmax_bv+ubwmax_bv+1), intent(inout) :: b_bv
    integer(kind=int32), intent(in) :: n, lbw, ubw, lbwmax_bv, ubwmax_bv, lbwmax_ub, ubwmax_ub
    integer(kind=int32), dimension(n), intent(in) :: numrotsv
    integer(kind=int32), dimension(n,ubwmax_bv), intent(in) :: ksv
    real(kind=dp), dimension(n,ubwmax_bv), intent(in) :: csv
    complex(kind=dp), dimension(n,ubwmax_bv), intent(in) :: ssv

    complex(kind=dp), dimension(lbwmax_ub+ubwmax_ub+1,n), intent(out) :: b_ub
    integer(kind=int32), dimension(n), intent(out) :: numrotsu
    integer(kind=int32), dimension(ubwmax_ub,n), intent(out) :: jsu
    integer(kind=int32), intent(out) :: lbw_ub, ubw_ub
    real(kind=dp), dimension(ubwmax_ub,n), intent(out) :: csu
    complex(kind=dp), dimension(ubwmax_ub,n), intent(out) :: ssu
    type(error_info), intent(out) :: error

    integer(kind=int32) :: j, k, ubw1, lbw1, k0, k1
    type(c_rotation) :: rot

    call clear_error(error)
    b_ub(1:lbw+ubw+1,:)=(0.0_dp, 0.0_dp); numrotsu=0
    ssu(1:ubw,:)=(0.0_dp, 0.0_dp); csu(1:ubw,:)=0.0_dp
    jsu(1:ubw,:)=0
    lbw_ub=lbw; ubw_ub=ubw
    if (n == 1) then
       b_ub(1,1)=b_bv(1,1);
       lbw_ub=0; ubw_ub=0; numrotsu=0; return
    end if
    lbw1=lbw
    if (ubw < n-1) then
       ubw1=ubw+1
       b_bv(:,ubw1+lbw1+1)=0.0_dp
    else
       ubw1=ubw
    end if
    ! apply v_{n-1}
    do k=1,n-2
       do j=1,numrotsv(k)
          rot%cosine=csv(k,j); rot%sine=ssv(k,j)
          call tbr_times_rotation(b_bv,n,lbw1,ubw1,0,n-k,trp_rot(rot),ksv(k,j))
       end do
       !
       k0=max(k+2,ubw1+1)
       k1=min(k+ubw1,n)
       numrotsu(k+1)=max(k1-k0+1,0)
       do j=k1,k0,-1
          rot=lgivens2(get_el_br(b_bv,lbw1,j-ubw1,j), get_el_br(b_bv,lbw1,j-ubw1+1,j))
          jsu(j-k0+1,k+1)=j-ubw1
          csu(j-k0+1,k+1)=rot%cosine; ssu(j-k0+1,k+1)=rot%sine
          call rotation_times_tbr(trp_rot(rot),b_bv,n,lbw1,ubw1,k+1,0,j-ubw1)
       end do
    end do
    call br_to_bc(b_bv,b_ub,lbw,ubw)
  end subroutine f_c_convert_bv_to_ub

end module mod_convert_bv_to_ub
