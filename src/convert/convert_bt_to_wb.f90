module mod_convert_bt_to_wb
  use mod_prec
  use mod_error_id
  use mod_shift
  use mod_rotation
  use mod_orth_band_types
  use mod_band_types
  implicit none

  private

  public :: convert_bt_to_wb, d_convert_bt_to_wb, z_convert_bt_to_wb, &
       f_convert_bt_to_wb, f_d_convert_bt_to_wb, f_z_convert_bt_to_wb, &
       d_wb_of_bt, z_wb_of_bt, wb_of

  interface convert_bt_to_wb
     module procedure d_convert_bt_to_wb, z_convert_bt_to_wb
  end interface convert_bt_to_wb

  interface f_convert_bt_to_wb
     module procedure f_d_convert_bt_to_wb, f_z_convert_bt_to_wb
  end interface f_convert_bt_to_wb

  interface wb_of
     module procedure d_wb_of_bt, z_wb_of_bt
  end interface wb_of

contains

  function d_wb_of_bt(bt,error) result(wb)
    type(d_wb) :: wb
    type(d_bt), intent(in) :: bt
    type(error_info), intent(inout), optional :: error

    type(d_bt), allocatable :: bt1
    type(routine_info), parameter :: info=info_d_wb_of_bt

    if (failure(error)) return
    call push_id(info, error)
    
    wb=d_new_wb(get_n(bt), bt%lbw, bt%ubw)
    bt1=d_new_bt(get_n(bt),bt%lbw+1,  bt%ubw)
    call copy(bt1,bt)
    call d_convert_bt_to_wb(bt1,wb,error)

    call pop_id(error)
  end function d_wb_of_bt


  ! Errors:
  ! 0: no error
  ! 1: n<1
  ! 2: Insufficient storage in bt
  ! 3: Insufficient stroage in wb
  ! 4: wb%n /= bt%n
  subroutine d_convert_bt_to_wb(bt, wb, error)
    type(d_bt) :: bt
    type(d_wb) :: wb
    type(error_info), intent(inout), optional :: error
    type(routine_info), parameter :: info=info_d_convert_bt_to_wb
    
    if (failure(error)) return
    call push_id(info, error)

    if (get_n(wb) < 1) then
       call set_error(1, info, error); return
    end if
    ! must allow for temporary fill-in
    if (get_lbwmax(bt) < bt%lbw+1 .and. bt%lbw < get_n(bt)-1) then
       call set_error(2, info, error); return
    end if
    if (get_ubwmax(wb) < bt%ubw .or. get_lbwmax(wb) < bt%lbw) then
       call set_error(3, info, error); return
    end if
    if (get_n(wb) /= get_n(bt)) then
       call set_error(4, info, error); return
    end if

    call f_d_convert_bt_to_wb(bt%br, get_n(bt), bt%lbw, bt%ubw, get_lbwmax(bt), &
         get_ubwmax(bt), bt%numrotst, bt%kst, bt%cst, bt%sst, wb%bc,  wb%lbw, &
         wb%ubw, get_lbwmax(wb), &
         get_ubwmax(wb), wb%numrotsw, wb%jsw, wb%csw, wb%ssw)
    call pop_id(error)
  end subroutine d_convert_bt_to_wb

  subroutine f_d_convert_bt_to_wb(b_bt, n, lbw, ubw, lbwmax_bt, ubwmax_bt, &
       numrotst, kst, cst, sst, b_wb, lbw_wb, ubw_wb, lbwmax_wb, ubwmax_wb, numrotsw, jsw, &
       csw, ssw)
    real(kind=dp), dimension(n,lbwmax_bt+ubwmax_bt+1), intent(inout) :: b_bt
    integer(kind=int32), intent(in) :: n, lbw, ubw, lbwmax_wb, ubwmax_wb, lbwmax_bt, ubwmax_bt
    integer(kind=int32), dimension(n), intent(in) :: numrotst
    integer(kind=int32), dimension(n,lbwmax_bt), intent(in) :: kst
    real(kind=dp), dimension(n,lbwmax_bt), intent(in) :: cst, sst

    real(kind=dp), dimension(lbwmax_wb+ubwmax_wb+1,n), intent(out) :: b_wb
    integer(kind=int32), dimension(n), intent(out) :: numrotsw
    integer(kind=int32), dimension(lbwmax_wb,n), intent(out) :: jsw
    real(kind=dp), dimension(lbwmax_wb,n), intent(out) :: csw, ssw
    integer(kind=int32), intent(out) :: lbw_wb, ubw_wb

    integer(kind=int32) :: j, k, ubw1, lbw1, j0, j1
    type(d_rotation) :: rot
    logical :: full_lbw

    b_wb(1:lbw+ubw+1,:)=0.0_dp; numrotsw=0
    ssw(1:lbw,:)=0.0_dp; csw(:,1:lbw)=0.0_dp
    jsw(1:lbw,:)=0
    lbw_wb=lbw; ubw_wb=ubw
    if (n == 1) then
       b_wb(1,1)=b_bt(1,1);
       lbw_wb=0; ubw_wb=0; return
    end if
    ubw1=ubw
    if (lbw < n-1) then
       lbw1=lbw+1
       call shift(b_bt,0,1)
       full_lbw=.false.
    else
       lbw1=lbw
       full_lbw=.true.
    end if

    ! k is the size of the leading principal submatrix
    do k=n-1,2,-1
       ! Apply T_k
       do j=1,numrotst(k)
          rot%cosine=cst(k,j); rot%sine=sst(k,j)
          call tbr_times_rotation(b_bt,n,lbw1,ubw1,k,0,trp_rot(rot),kst(k,j))
       end do
       ! Rows in which nonzeros have been introduced into the extra subdiagonal.
       j0=max(k+1,lbw+2)
       j1=min(k+lbw,n)
       ! Apply W_{k-1}
       numrotsw(k-1)=max(j1-j0+1,0)
       do j=j0,j1
          rot=lgivens(get_el_br(b_bt,lbw1,j-1,j-lbw1),get_el_br(b_bt,lbw1,j,j-lbw1))
          jsw(j1-j+1,k-1)=j-1
          csw(j1-j+1,k-1)=rot%cosine; ssw(j1-j+1,k-1)=rot%sine
          call rotation_times_tbr(trp_rot(rot),b_bt,n,lbw1,ubw1,0,n-k+1,j-1)
       end do
    end do
    if (.not. full_lbw) then
       call shift(b_bt,0,-1)
    end if
    call br_to_bc(b_bt,b_wb,lbw,ubw)
  end subroutine f_d_convert_bt_to_wb

  function z_wb_of_bt(bt,error) result(wb)
    type(z_wb) :: wb
    type(z_bt), intent(in) :: bt
    type(error_info), intent(inout), optional :: error

    type(z_bt), allocatable :: bt1
    type(routine_info), parameter :: info=info_z_wb_of_bt

    if (failure(error)) return
    call push_id(info, error)
    
    wb=z_new_wb(get_n(bt), bt%lbw, bt%ubw)
    bt1=z_new_bt(get_n(bt),bt%lbw+1,  bt%ubw)
    call copy(bt1,bt)
    call z_convert_bt_to_wb(bt1,wb,error)

    call pop_id(error)
  end function z_wb_of_bt
  

  ! Errors:
  ! 0: no error
  ! 1: n<1
  ! 2: Insufficient storage in bt
  ! 3: Insufficient stroage in wb
  ! 4: wb%n /= bt%n
  subroutine z_convert_bt_to_wb(bt, wb, error)
    type(z_bt) :: bt
    type(z_wb) :: wb
    type(error_info), intent(inout), optional :: error
    type(routine_info), parameter :: info=info_z_convert_bt_to_wb
    
    if (failure(error)) return
    call push_id(info, error)
    
    if (get_n(wb) < 1) then
       call set_error(1, info, error); return
    end if
    ! must allow for temporary fill-in
    if (get_lbwmax(bt) < bt%lbw+1 .and. bt%lbw < get_n(bt)-1) then
       call set_error(2, info, error); return
    end if
    if (get_ubwmax(wb) < bt%ubw .or. get_lbwmax(wb) < bt%lbw) then
       call set_error(3, info, error); return
    end if
    if (get_n(wb) /= get_n(bt)) then
       call set_error(4, info, error); return
    end if

    call f_z_convert_bt_to_wb(bt%br, get_n(bt), bt%lbw, bt%ubw, get_lbwmax(bt), &
         get_ubwmax(bt), bt%numrotst, bt%kst, bt%cst, bt%sst, wb%bc,  wb%lbw, &
         wb%ubw, get_lbwmax(wb), &
         get_ubwmax(wb), wb%numrotsw, wb%jsw, wb%csw, wb%ssw)
    call pop_id(error)
  end subroutine z_convert_bt_to_wb

  subroutine f_z_convert_bt_to_wb(b_bt, n, lbw, ubw, lbwmax_bt, ubwmax_bt, &
       numrotst, kst, cst, sst, b_wb, lbw_wb, ubw_wb, lbwmax_wb, ubwmax_wb, numrotsw, jsw, &
       csw, ssw)
    complex(kind=dp), dimension(n,lbwmax_bt+ubwmax_bt+1), intent(inout) :: b_bt
    integer(kind=int32), intent(in) :: n, lbw, ubw, lbwmax_wb, ubwmax_wb, lbwmax_bt, ubwmax_bt
    integer(kind=int32), dimension(n), intent(in) :: numrotst
    integer(kind=int32), dimension(n,lbwmax_bt), intent(in) :: kst
    real(kind=dp), dimension(n,lbwmax_bt), intent(in) :: cst
    complex(kind=dp), dimension(n,lbwmax_bt), intent(in) :: sst

    complex(kind=dp), dimension(lbwmax_wb+ubwmax_wb+1,n), intent(out) :: b_wb
    integer(kind=int32), dimension(n), intent(out) :: numrotsw
    integer(kind=int32), dimension(lbwmax_wb,n), intent(out) :: jsw
    real(kind=dp), dimension(lbwmax_wb,n), intent(out) :: csw
    complex(kind=dp), dimension(lbwmax_wb,n), intent(out) :: ssw
    integer(kind=int32), intent(out) :: lbw_wb, ubw_wb

    integer(kind=int32) :: j, k, ubw1, lbw1, j0, j1
    type(z_rotation) :: rot
    logical :: full_lbw

    b_wb(1:lbw+ubw+1,:)=(0.0_dp,0.0_dp); numrotsw=0
    ssw(1:lbw,:)=(0.0_dp,0.0_dp); csw(:,1:lbw)=0.0_dp
    jsw(1:lbw,:)=0
    lbw_wb=lbw; ubw_wb=ubw
    if (n == 1) then
       b_wb(1,1)=b_bt(1,1);
       lbw_wb=0; ubw_wb=0; numrotsw=0; return
    end if
    ubw1=ubw
    if (lbw < n-1) then
       lbw1=lbw+1
       call shift(b_bt,0,1)
       full_lbw=.false.
    else
       lbw1=lbw
       full_lbw=.true.
    end if

    ! k is the size of the leading principal submatrix
    do k=n-1,2,-1
       ! Apply T_k
       do j=1,numrotst(k)
          rot%cosine=cst(k,j); rot%sine=sst(k,j)
          call tbr_times_rotation(b_bt,n,lbw1,ubw1,k,0,trp_rot(rot),kst(k,j))
       end do
       ! Rows in which nonzeros have been introduced into the extra subdiagonal.
       j0=max(k+1,lbw+2)
       j1=min(k+lbw,n)
       ! Apply W_{k-1}
       numrotsw(k-1)=max(j1-j0+1,0)
       do j=j0,j1
          rot=lgivens(get_el_br(b_bt,lbw1,j-1,j-lbw1),get_el_br(b_bt,lbw1,j,j-lbw1))
          jsw(j1-j+1,k-1)=j-1
          csw(j1-j+1,k-1)=rot%cosine; ssw(j1-j+1,k-1)=rot%sine
          call rotation_times_tbr(trp_rot(rot),b_bt,n,lbw1,ubw1,0,n-k+1,j-1)
       end do
    end do
    if (.not. full_lbw) then
       call shift(b_bt,0,-1)
    end if
    call br_to_bc(b_bt,b_wb,lbw,ubw)
  end subroutine f_z_convert_bt_to_wb

end module mod_convert_bt_to_wb
