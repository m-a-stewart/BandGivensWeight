module row_compression
  use misc
  use sweeps
  use rotation
  use shift
  implicit none

  ! interface row_compress
  !    module procedure d_row_compress, c_row_compress
  ! end interface row_compress

  ! interface f_row_compress
  !    module procedure f_d_row_compress, f_c_row_compress
  ! end interface f_row_compress

  type(routine_info), parameter :: info_d_row_compress=routine_info(id_d_row_compress, &
       'd_row_compress', &
       [ character(len=error_message_length) :: 'ub%n /= bv%n or sw%n /= ub%n', &
       'Not enough stroage for each sweeps', 'Not enough storage for number of sweeps', &
       'Insufficient storage in ubt.', 'Insufficient storage in bv.'])

  type(routine_info), parameter :: info_c_row_compress=routine_info(id_c_row_compress, &
       'c_row_compress', &
       [ character(len=error_message_length) :: 'ub%n /= bv%n or sw%n /= ub%n', &
       'Not enough stroage for each sweeps', 'Not enough storage for number of sweeps', &
       'Insufficient storage in ubt.', 'Insufficient storage in bv.'])

contains


  ! Errors
  ! 0: no error
  ! 1: ub%n /= bv%n .or. sw%n /= ub%n
  ! 2: Not enough storage for each sweep.
  ! 3: Not enough storage for the number of sweeps.
  subroutine d_row_compress(ubt,bv,sw,error)
    type(d_ubt) :: ubt
    type(d_bv) :: bv
    type(d_sweeps) :: sw
    type(error_info), intent(out) :: error

    integer(kind=int32) :: n, lbw
    lbw=ubt%lbw
    n=get_n(ubt)
    call clear_error(error)
    if (n /= get_n(bv) .or. n /= get_n(sw)) then
       call set_error(error, 1, id_d_row_compress); return
    end if
    if (get_maxord(sw) < ubt%lbw) then
       call set_error(error, 2, id_d_row_compress); return
    end if
    if (get_maxnum(sw) < n-2) then
       call set_error(error, 3, id_d_row_compress); return
    end if
    if (get_ubwmax(ubt) < min(get_lbwmax(ubt)+get_ubwmax(ubt)+1,n-1)) then
       call set_error(error, 4, id_d_row_compress); return
    end if
    if (get_lbwmax(bv) < min(get_lbwmax(bv),n-1) .or. &
         get_ubwmax(bv) < min(get_lbwmax(ubt)+get_ubwmax(ubt),n-1)) then
       call set_error(error, 5, id_d_row_compress); return
    end if
    call f_d_row_compress(ubt%bc, get_n(ubt), ubt%lbw, ubt%ubw, get_lbwmax(ubt), &
         get_ubwmax(ubt), ubt%numrotsu, ubt%jsu, ubt%csu, ubt%ssu, & 
         ubt%numrotst, ubt%kst, ubt%cst, ubt%sst, & 
         bv%br, bv%lbw, bv%ubw, get_lbwmax(bv), get_ubwmax(bv), bv%numrotsv, bv%ksv, bv%csv, bv%ssv, &
         sw%num, get_maxnum(sw), get_maxord(sw), sw%numrots, sw%js, &
         sw%cs(:,1:lbw), sw%ss(:,1:lbw), error)
  end subroutine d_row_compress

  subroutine f_d_row_compress(b_ubt, n, lbw_ubt, ubw_ubt, lbwmax_ubt, ubwmax_ubt, numrotsu, &
       jsu, csu, ssu, numrotst, kst, cst, sst, &
       b_bv, lbw_bv, ubw_bv, lbwmax_bv, ubwmax_bv, numrotsv, ksv, csv, ssv, &
       num, maxnum, maxord, numrotsq, jsq, csq, ssq, error)
    integer(kind=int32), intent(in) :: n, lbw_ubt, ubw_ubt, lbwmax_ubt, ubwmax_ubt, lbwmax_bv, ubwmax_bv
    real(kind=dp), dimension(lbwmax_ubt+ubwmax_ubt+1,n), intent(inout) :: b_ubt
    integer(kind=int32), dimension(n), intent(in) :: numrotsu, numrotst
    integer(kind=int32), dimension(ubwmax_ubt,n), intent(in) :: jsu
    real(kind=dp), dimension(ubwmax_ubt,n), intent(in) :: csu, ssu
    integer(kind=int32), dimension(n,lbwmax_ubt), intent(in) :: kst
    real(kind=dp), dimension(n,lbwmax_ubt), intent(in) :: cst, sst

    real(kind=dp), dimension(n,lbwmax_bv+ubwmax_bv+1), intent(out) :: b_bv
    integer(kind=int32), dimension(n), intent(out) :: numrotsv
    integer(kind=int32), dimension(n,ubwmax_bv), intent(out) :: ksv
    real(kind=dp), dimension(n,ubwmax_bv), intent(out) :: csv, ssv
    integer(kind=int32), intent(out) :: lbw_bv, ubw_bv
    
    integer(kind=int32), intent(in) :: maxnum, maxord
    real(kind=dp), dimension(maxord,maxnum) :: csq, ssq
    integer(kind=int32), dimension(n), intent(out) :: numrotsq
    integer(kind=int32), dimension(maxord,maxnum), intent(out) :: jsq
    integer(kind=int32), intent(out) :: num
    type(error_info), intent(out) :: error
    integer(kind=int32) :: j, lbw, ubw
    logical :: full_ubw, full_lbw


    call clear_error(error)
    if (n == 1) then
       b_bv(1,1)=b_ubt(1,1);
       lbw_bv=0; ubw_bv=0; numrotsv=0; return
    end if

    ubw=min(n-1,ubw_ubt+lbw_ubt)
    call shift2(b_ubt,ubw-ubw_ubt,0)

    if (ubw < n-1) then
       ubw=ubw+1
       full_ubw=.false.
       call shift2(b_ubt,1,0)
    else
       full_ubw=.true.
    end if
    
    lbw=lbw_ubt
    if (lbw < n-1) then
       lbw=lbw+1
       full_lbw=.false.
       b_ubt(lbw+ubw+1,:) = 0.0_dp
    else
       full_lbw=.true.
    end if

    ! if (lbw <= 0) then
    !    call f_convert_bv_to_ub(b_bv, n, lbw, ubw, lbwmax_bv, ubwmax_bv, &
    !         numrotsv, ksv, csv, ssv, b_ub, lbw_ub, ubw_ub, &
    !         lbwmax_ub, ubwmax_ub, numrotsu, jsu, &
    !         csu, ssu, error)
    !    return
    ! end if
    ! do j=1,lbw-1
    !    call f_d_reduce_lbw_bv_to_ub(b_bv, n, lbw, ubw, lbwmax_bv, ubwmax_bv, numrotsv, &
    !         ksv, csv, ssv, &
    !         b_ub, lbw_ub, ubw_ub, lbwmax_ub, ubwmax_ub, numrotsu, jsu, &
    !         csu, ssu, cs(:,j), ss(:,j), error)
    !    call f_d_convert_ub_to_bv(b_ub, n, lbw_ub, ubw_ub, lbwmax_ub, ubwmax_ub, numrotsu, &
    !         jsu, csu, ssu, b_bv, lbw, ubw, lbwmax_bv, ubwmax_bv, numrotsv, ksv, &
    !         csv, ssv, error)
    ! end do
    ! call f_d_reduce_lbw_bv_to_ub(b_bv, n, lbw, ubw, lbwmax_bv, ubwmax_bv, numrotsv, &
    !      ksv, csv, ssv, &
    !      b_ub, lbw_ub, ubw_ub, lbwmax_ub, ubwmax_ub, numrotsu, jsu, &
    !      csu, ssu, cs(:,j), ss(:,j), error)
  end subroutine f_d_row_compress

end module row_compression

