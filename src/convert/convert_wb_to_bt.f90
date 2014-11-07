module mod_convert_wb_to_bt
  use mod_prec
  use mod_error_id
  use mod_rotation
  use mod_nested_types
  use mod_band_types
  implicit none

  private

  public :: convert_wb_to_bt, d_convert_wb_to_bt, c_convert_wb_to_bt, &
       f_convert_wb_to_bt, f_d_convert_wb_to_bt, f_c_convert_wb_to_bt

  public :: info_d_convert_wb_to_bt, info_c_convert_wb_to_bt

  interface convert_wb_to_bt
     module procedure d_convert_wb_to_bt, c_convert_wb_to_bt
  end interface convert_wb_to_bt

  interface f_convert_wb_to_bt
     module procedure f_d_convert_wb_to_bt, f_c_convert_wb_to_bt
  end interface f_convert_wb_to_bt

  type(routine_info), parameter :: info_d_convert_wb_to_bt=routine_info(id_d_convert_wb_to_bt, &
       'd_convert_wb_to_bt', &
       [ character(len=error_message_length) :: 'n<1', 'Insufficient storage in wb.', &
       'Insufficient Storage in bt.', 'bt%n /= wb%n' ] )

  type(routine_info), parameter :: info_c_convert_wb_to_bt=routine_info(id_c_convert_wb_to_bt, &
       'c_convert_wb_to_bt', &
       [ character(len=error_message_length) :: 'n<1', 'Insufficient storage in wb.', &
       'Insufficient Storage in bt.', 'bt%n /= wb%n' ] )

contains

  ! Errors:
  ! 0: no error
  ! 1: n<1
  ! 2: Insufficient storage in wb
  ! 3: Insufficient stroage in bt
  ! 4: bt%n /= wb%n
  subroutine d_convert_wb_to_bt(wb, bt, error)
    type(d_wb) :: wb
    type(d_bt) :: bt
    type(error_info), intent(out) :: error
    call clear_error(error)
    if (get_n(bt) < 1) then
       call set_error(error, 1, id_d_convert_wb_to_bt); return
    end if
    ! must allow for temporary fill-in
    if (get_lbwmax(wb) < wb%lbw+1 .and. wb%lbw < get_n(wb)-1) then
       call set_error(error, 2, id_d_convert_wb_to_bt); return
    end if
    if (get_ubwmax(bt) < wb%ubw .or. get_lbwmax(bt) < wb%lbw) then
       call set_error(error, 3, id_d_convert_wb_to_bt); return
    end if
    if (get_n(bt) /= get_n(wb)) then
       call set_error(error, 4, id_d_convert_wb_to_bt); return
    end if

    call f_d_convert_wb_to_bt(wb%bc, get_n(wb), wb%lbw, wb%ubw, get_lbwmax(wb), &
         get_ubwmax(wb), wb%numrotsw, wb%jsw, wb%csw, wb%ssw, bt%br,  bt%lbw, bt%ubw, get_lbwmax(bt), &
         get_ubwmax(bt), bt%numrotst, bt%kst, bt%cst, bt%sst, error)
  end subroutine d_convert_wb_to_bt

  subroutine f_d_convert_wb_to_bt(b_wb, n, lbw, ubw, lbwmax_wb, ubwmax_wb, &
       numrotsw, jsw, csw, ssw, b_bt, lbw_bt, ubw_bt, lbwmax_bt, ubwmax_bt, numrotst, kst, &
       cst, sst, error)
    real(kind=dp), dimension(lbwmax_wb+ubwmax_wb+1,n), intent(inout) :: b_wb
    integer(kind=int32), intent(in) :: n, lbw, ubw, lbwmax_bt, ubwmax_bt, lbwmax_wb, ubwmax_wb
    integer(kind=int32), dimension(n), intent(in) :: numrotsw
    integer(kind=int32), dimension(lbwmax_wb,n), intent(in) :: jsw
    real(kind=dp), dimension(lbwmax_wb,n), intent(in) :: csw, ssw

    real(kind=dp), dimension(n,lbwmax_bt+ubwmax_bt+1), intent(out) :: b_bt
    integer(kind=int32), dimension(n), intent(out) :: numrotst
    integer(kind=int32), dimension(n,lbwmax_bt), intent(out) :: kst
    integer(kind=int32), intent(out) :: lbw_bt, ubw_bt
    real(kind=dp), dimension(n,lbwmax_bt), intent(out) :: cst, sst
    type(error_info), intent(out) :: error

    integer(kind=int32) :: j, k, ubw1, lbw1, k0, k1
    type(d_rotation) :: rot

    call clear_error(error)
    b_bt(:,1:lbw+ubw+1)=0.0_dp; numrotst=0
    sst(:,1:lbw)=0.0_dp; cst(:,1:lbw)=0.0_dp
    kst(:,1:lbw)=0
    lbw_bt=lbw; ubw_bt=ubw
    if (n == 1) then
       b_bt(1,1)=b_wb(1,1);
       lbw_bt=0; ubw_bt=0; numrotst=0; return
    end if
    ubw1=ubw
    if (lbw < n-1) then
       lbw1=lbw+1
       b_wb(lbw1+ubw1+1,:)=0.0_dp
    else
       lbw1=lbw
    end if
    ! k is the size of the leading principal submatrix
    do k=1,n-2
       ! Apply W_k
       do j=1,numrotsw(k)
          rot%cosine=csw(j,k); rot%sine=ssw(j,k)
          call rotation_times_tbc(rot,b_wb,n,lbw1,ubw1,0,n-k,jsw(j,k))
       end do
       ! columns in which nonzeros have been introduced into the extra subdiagonal
       k0=max(k-lbw+1,1)
       k1=min(k,n-lbw1)
       ! Apply T_{k+1}
       numrotst(k+1)=max(k1-k0+1,0)
       do j=k1,k0,-1
          rot=rgivens2(get_el_bc(b_wb,ubw1,j+lbw1,j), get_el_bc(b_wb,ubw1,j+lbw1,j+1))
          kst(k+1,j-k0+1)=j
          cst(k+1,j-k0+1)=rot%cosine; sst(k+1,j-k0+1)=rot%sine
          call tbc_times_rotation(b_wb,n,lbw1,ubw1,k+1,0,rot,j)
       end do
    end do
    call bc_to_br(b_wb,b_bt,lbw,ubw)
  end subroutine f_d_convert_wb_to_bt


  ! Errors:
  ! 0: no error
  ! 1: n<1
  ! 2: Insufficient storage in wb
  ! 3: Insufficient stroage in bt
  ! 4: bt%n /= wb%n
  subroutine c_convert_wb_to_bt(wb, bt, error)
    type(c_wb) :: wb
    type(c_bt) :: bt
    type(error_info), intent(out) :: error
    call clear_error(error)
    if (get_n(bt) < 1) then
       call set_error(error, 1, id_c_convert_wb_to_bt); return
    end if
    ! must allow for temporary fill-in
    if (get_lbwmax(wb) < wb%lbw+1 .and. wb%lbw < get_n(wb)-1) then
       call set_error(error, 2, id_c_convert_wb_to_bt); return
    end if
    if (get_ubwmax(bt) < wb%ubw .or. get_lbwmax(bt) < wb%lbw) then
       call set_error(error, 3, id_c_convert_wb_to_bt); return
    end if
    if (get_n(bt) /= get_n(wb)) then
       call set_error(error, 4, id_c_convert_wb_to_bt); return
    end if

    call f_c_convert_wb_to_bt(wb%bc, get_n(wb), wb%lbw, wb%ubw, get_lbwmax(wb), &
         get_ubwmax(wb), wb%numrotsw, wb%jsw, wb%csw, wb%ssw, bt%br,  bt%lbw, bt%ubw, get_lbwmax(bt), &
         get_ubwmax(bt), bt%numrotst, bt%kst, bt%cst, bt%sst, error)
  end subroutine c_convert_wb_to_bt

  subroutine f_c_convert_wb_to_bt(b_wb, n, lbw, ubw, lbwmax_wb, ubwmax_wb, &
       numrotsw, jsw, csw, ssw, b_bt, lbw_bt, ubw_bt, lbwmax_bt, ubwmax_bt, numrotst, kst, &
       cst, sst, error)
    complex(kind=dp), dimension(lbwmax_wb+ubwmax_wb+1,n), intent(inout) :: b_wb
    integer(kind=int32), intent(in) :: n, lbw, ubw, lbwmax_bt, ubwmax_bt, lbwmax_wb, ubwmax_wb
    integer(kind=int32), dimension(n), intent(in) :: numrotsw
    integer(kind=int32), dimension(lbwmax_wb,n), intent(in) :: jsw
    real(kind=dp), dimension(lbwmax_wb,n), intent(in) :: csw
    complex(kind=dp), dimension(lbwmax_wb,n), intent(in) :: ssw

    complex(kind=dp), dimension(n,lbwmax_bt+ubwmax_bt+1), intent(out) :: b_bt
    integer(kind=int32), dimension(n), intent(out) :: numrotst
    integer(kind=int32), dimension(n,lbwmax_bt), intent(out) :: kst
    integer(kind=int32), intent(out) :: lbw_bt, ubw_bt
    real(kind=dp), dimension(n,lbwmax_bt), intent(out) :: cst
    complex(kind=dp), dimension(n,lbwmax_bt), intent(out) :: sst
    type(error_info), intent(out) :: error

    integer(kind=int32) :: j, k, ubw1, lbw1, k0, k1
    type(c_rotation) :: rot

    call clear_error(error)
    b_bt(:,1:lbw+ubw+1)=(0.0_dp,0.0_dp); numrotst=0
    sst(:,1:lbw)=(0.0_dp,0.0_dp); cst(:,1:lbw)=0.0_dp
    kst(:,1:lbw)=0
    lbw_bt=lbw; ubw_bt=ubw
    if (n == 1) then
       b_bt(1,1)=b_wb(1,1);
       lbw_bt=0; ubw_bt=0; numrotst=0; return
    end if
    ubw1=ubw
    if (lbw < n-1) then
       lbw1=lbw+1
       b_wb(lbw1+ubw1+1,:)=(0.0_dp,0.0_dp)
    else
       lbw1=lbw
    end if
    ! k is the size of the leading principal submatrix
    do k=1,n-2
       ! Apply W_k
       do j=1,numrotsw(k)
          rot%cosine=csw(j,k); rot%sine=ssw(j,k)
          call rotation_times_tbc(rot,b_wb,n,lbw1,ubw1,0,n-k,jsw(j,k))
       end do
       ! columns in which nonzeros have been introduced into the extra subdiagonal
       k0=max(k-lbw+1,1)
       k1=min(k,n-lbw1)
       ! Apply T_{k+1}
       numrotst(k+1)=max(k1-k0+1,0)
       do j=k1,k0,-1
          rot=rgivens2(get_el_bc(b_wb,ubw1,j+lbw1,j), get_el_bc(b_wb,ubw1,j+lbw1,j+1))
          kst(k+1,j-k0+1)=j
          cst(k+1,j-k0+1)=rot%cosine; sst(k+1,j-k0+1)=rot%sine
          call tbc_times_rotation(b_wb,n,lbw1,ubw1,k+1,0,rot,j)
       end do
    end do
    call bc_to_br(b_wb,b_bt,lbw,ubw)
  end subroutine f_c_convert_wb_to_bt


end module mod_convert_wb_to_bt
