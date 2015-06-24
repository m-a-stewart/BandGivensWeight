module mod_submatrix
  use mod_prec
  use mod_error_id
  use mod_band_types
  use mod_orth_band_types
  use mod_rotation
  use mod_utility
  implicit none
  ! This module contains routines for selecting
  ! structured and unstructured submatrices of
  ! various Givens-weight representations.

  private

  public :: leading, d_ub_leading, z_ub_leading, d_bt_leading, z_bt_leading, &
       d_ubt_leading, z_ubt_leading, &
       trailing, d_bv_trailing, z_bv_trailing, d_wb_trailing, z_wb_trailing, &
       d_wbv_trailing, z_wbv_trailing

  public :: columns_of_ub, ub_to_columns, &
       d_columns_of_ub, d_ub_to_columns, f_d_ub_to_columns, &
       z_columns_of_ub, z_ub_to_columns, f_z_ub_to_columns       

  public :: columns_of_wb, wb_to_columns, &
       d_columns_of_wb, d_wb_to_columns, f_d_wb_to_columns, &
       z_columns_of_wb, z_wb_to_columns, f_z_wb_to_columns

  public :: rows_of_bv, bv_to_rows, &
       d_rows_of_bv, d_bv_to_rows, f_d_bv_to_rows, &
       z_rows_of_bv, z_bv_to_rows, f_z_bv_to_rows

  public :: rows_of_bt, bt_to_rows, &
       d_rows_of_bt, d_bt_to_rows, f_d_bt_to_rows, &
       z_rows_of_bt, z_bt_to_rows, f_z_bt_to_rows
  
  interface leading
     module procedure d_ub_leading, z_ub_leading, d_bt_leading, z_bt_leading, &
          d_ubt_leading, z_ubt_leading
  end interface leading

  interface trailing
     module procedure d_bv_trailing, z_bv_trailing, d_wb_trailing, z_wb_trailing, &
          d_wbv_trailing, z_wbv_trailing
  end interface trailing

  interface columns_of_ub
     module procedure d_columns_of_ub, z_columns_of_ub
  end interface columns_of_ub

  interface ub_to_columns
     module procedure d_ub_to_columns, z_ub_to_columns
  end interface ub_to_columns

  interface rows_of_bv
     module procedure d_rows_of_bv, z_rows_of_bv
  end interface rows_of_bv

  interface bv_to_rows
     module procedure d_bv_to_rows, z_bv_to_rows
  end interface bv_to_rows

  interface columns_of_wb
     module procedure d_columns_of_wb, z_columns_of_wb
  end interface columns_of_wb

  interface wb_to_columns
     module procedure d_wb_to_columns, z_wb_to_columns
  end interface wb_to_columns

  interface rows_of_bt
     module procedure d_rows_of_bt, z_rows_of_bt
  end interface rows_of_bt

  interface bt_to_rows
     module procedure d_bt_to_rows, z_bt_to_rows
  end interface bt_to_rows
  
contains

  ! Structured leading and trailing principal submatrices.

  type(d_ub) function d_ub_leading(ub,l,lbwmaxl0,ubwmaxl0,error) result(ubl)
    type(d_ub), intent(in) :: ub
    integer(kind=int32), intent(in) :: l
    integer(kind=int32), intent(in), optional :: lbwmaxl0,ubwmaxl0
    type(error_info), intent(inout), optional :: error
    
    integer(kind=int32) :: lbwmax, ubwmax, lbwmaxl, ubwmaxl, n, lbw, ubw
    type(routine_info), parameter :: info=info_d_ub_leading

    if (failure(error)) return
    call push_id(info,error)

    lbwmax=get_lbwmax(ub)
    ubwmax=get_ubwmax(ub)
    n=get_n(ub)
    lbw=ub%lbw
    ubw=ub%ubw
    if (n < 1) then
       call set_error(1, info, error); return
    end if

    lbwmaxl=equals_option(lbw, lbwmaxl0)
    ubwmaxl=equals_option(ubw, ubwmaxl0)
    if (lbwmaxl < lbw .or. ubwmaxl < ubw) then
       call set_error(2, info, error); return
    end if

    ubl = d_new_ub(l,lbwmaxl,ubwmaxl)

    ubl%bc(1:lbw+ubw+1,:)=ub%bc(1:lbw+ubw+1,1:l)
    ubl%lbw=lbw
    ubl%ubw=ubw
    ubl%jsu(1:ubw,1:l-1)=ub%jsu(1:ubw,1:l-1)
    ubl%csu(1:ubw,1:l-1)=ub%csu(1:ubw,1:l-1)
    ubl%ssu(1:ubw,1:l-1)=ub%ssu(1:ubw,1:l-1)
    ubl%numrotsu=ub%numrotsu(1:l-1)
    call pop_id(error)
  end function d_ub_leading

  type(z_ub) function z_ub_leading(ub,l,lbwmaxl0,ubwmaxl0,error) result(ubl)
    type(z_ub), intent(in) :: ub
    integer(kind=int32), intent(in) :: l
    integer(kind=int32), intent(in), optional :: lbwmaxl0,ubwmaxl0
    type(error_info), intent(inout), optional :: error
    
    integer(kind=int32) :: lbwmax, ubwmax, lbwmaxl, ubwmaxl, n, lbw, ubw
    type(routine_info), parameter :: info=info_z_ub_leading

    if (failure(error)) return
    call push_id(info,error)

    lbwmax=get_lbwmax(ub)
    ubwmax=get_ubwmax(ub)
    n=get_n(ub)
    lbw=ub%lbw
    ubw=ub%ubw
    if (n < 1) then
       call set_error(1, info, error); return
    end if

    lbwmaxl=equals_option(lbw, lbwmaxl0)
    ubwmaxl=equals_option(ubw, ubwmaxl0)
    if (lbwmaxl < lbw .or. ubwmaxl < ubw) then
       call set_error(2, info, error); return
    end if

    ubl = z_new_ub(l,lbwmaxl,ubwmaxl)

    ubl%bc(1:lbw+ubw+1,:)=ub%bc(1:lbw+ubw+1,1:l)
    ubl%lbw=lbw
    ubl%ubw=ubw
    ubl%jsu(1:ubw,1:l-1)=ub%jsu(1:ubw,1:l-1)
    ubl%csu(1:ubw,1:l-1)=ub%csu(1:ubw,1:l-1)
    ubl%ssu(1:ubw,1:l-1)=ub%ssu(1:ubw,1:l-1)
    ubl%numrotsu=ub%numrotsu(1:l-1)
    call pop_id(error)
  end function z_ub_leading

  type(d_bv) function d_bv_trailing(bv,l,lbwmaxl0,ubwmaxl0,error) result(bvl)
    type(d_bv), intent(in) :: bv
    integer(kind=int32), intent(in) :: l
    integer(kind=int32), intent(in), optional :: lbwmaxl0,ubwmaxl0
    type(error_info), intent(inout), optional :: error
    
    integer(kind=int32) :: lbwmax, ubwmax, lbwmaxl, ubwmaxl, n, lbw, ubw, j, k
    type(routine_info), parameter :: info=info_d_bv_trailing

    if (failure(error)) return
    call push_id(info,error)

    lbwmax=get_lbwmax(bv)
    ubwmax=get_ubwmax(bv)
    n=get_n(bv)
    lbw=bv%lbw
    ubw=bv%ubw
    if (n < 1) then
       call set_error(1, info, error); return
    end if

    lbwmaxl=equals_option(lbw, lbwmaxl0)
    ubwmaxl=equals_option(ubw, ubwmaxl0)
    if (lbwmaxl < lbw .or. ubwmaxl < ubw) then
       call set_error(2, info, error); return
    end if

    bvl = d_new_bv(l,lbwmaxl,ubwmaxl)

    bvl%br(:,1:lbw+ubw+1)=bv%br(n-l+1:n,1:lbw+ubw+1)
    bvl%lbw=lbw
    bvl%ubw=ubw
    bvl%ksv(:,1:ubw)=bv%ksv(n-l+1:n,1:ubw)
    bvl%csv(:,1:ubw)=bv%csv(n-l+1:n,1:ubw)
    bvl%ssv(:,1:ubw)=bv%ssv(n-l+1:n,1:ubw)
    bvl%numrotsv=bv%numrotsv(n-l+1:n)
    do j=1,l
       do k=1,bvl%numrotsv(j)
          bvl%ksv(j,k)=bvl%ksv(j,k)-(n-l)
       end do
    end do
    call pop_id(error)
  end function d_bv_trailing  

  type(z_bv) function z_bv_trailing(bv,l,lbwmaxl0,ubwmaxl0,error) result(bvl)
    type(z_bv), intent(in) :: bv
    integer(kind=int32), intent(in) :: l
    integer(kind=int32), intent(in), optional :: lbwmaxl0,ubwmaxl0
    type(error_info), intent(inout), optional :: error
    
    integer(kind=int32) :: lbwmax, ubwmax, lbwmaxl, ubwmaxl, n, lbw, ubw, j, k
    type(routine_info), parameter :: info=info_z_bv_trailing

    if (failure(error)) return
    call push_id(info,error)

    lbwmax=get_lbwmax(bv)
    ubwmax=get_ubwmax(bv)
    n=get_n(bv)
    lbw=bv%lbw
    ubw=bv%ubw
    if (n < 1) then
       call set_error(1, info, error); return
    end if

    lbwmaxl=equals_option(lbw, lbwmaxl0)
    ubwmaxl=equals_option(ubw, ubwmaxl0)
    if (lbwmaxl < lbw .or. ubwmaxl < ubw) then
       call set_error(2, info, error); return
    end if

    bvl = z_new_bv(l,lbwmaxl,ubwmaxl)

    bvl%br(:,1:lbw+ubw+1)=bv%br(n-l+1:n,1:lbw+ubw+1)
    bvl%lbw=lbw
    bvl%ubw=ubw
    bvl%ksv(:,1:ubw)=bv%ksv(n-l+1:n,1:ubw)
    bvl%csv(:,1:ubw)=bv%csv(n-l+1:n,1:ubw)
    bvl%ssv(:,1:ubw)=bv%ssv(n-l+1:n,1:ubw)
    bvl%numrotsv=bv%numrotsv(n-l+1:n)
    do j=1,l
       do k=1,bvl%numrotsv(j)
          bvl%ksv(j,k)=bvl%ksv(j,k)-(n-l)
       end do
    end do
    call pop_id(error)
  end function z_bv_trailing

  type(d_bt) function d_bt_leading(bt,l,lbwmaxl0,ubwmaxl0,error) result(btl)
    type(d_bt), intent(in) :: bt
    integer(kind=int32), intent(in) :: l
    integer(kind=int32), intent(in), optional :: lbwmaxl0,ubwmaxl0
    type(error_info), intent(inout), optional :: error
    
    integer(kind=int32) :: lbwmax, ubwmax, lbwmaxl, ubwmaxl, n, lbw, ubw
    type(routine_info), parameter :: info=info_d_bt_leading

    if (failure(error)) return
    call push_id(info,error)

    lbwmax=get_lbwmax(bt)
    ubwmax=get_ubwmax(bt)
    n=get_n(bt)
    lbw=bt%lbw
    ubw=bt%ubw
    if (n < 1) then
       call set_error(1, info, error); return
    end if

    lbwmaxl=equals_option(lbw, lbwmaxl0)
    ubwmaxl=equals_option(ubw, ubwmaxl0)
    if (lbwmaxl < lbw .or. ubwmaxl < ubw) then
       call set_error(2, info, error); return
    end if

    btl = d_new_bt(l,lbwmaxl,ubwmaxl)

    btl%br(:,1:lbw+ubw+1)=bt%br(1:l,1:lbw+ubw+1)
    btl%lbw=lbw
    btl%ubw=ubw
    btl%kst(1:l-1,1:lbw)=bt%kst(1:l-1,1:lbw)
    btl%cst(1:l-1,1:lbw)=bt%cst(1:l-1,1:lbw)
    btl%sst(1:l-1,1:lbw)=bt%sst(1:l-1,1:lbw)
    btl%numrotst=bt%numrotst(1:l-1)
    call pop_id(error)
  end function d_bt_leading

  type(z_bt) function z_bt_leading(bt,l,lbwmaxl0,ubwmaxl0,error) result(btl)
    type(z_bt), intent(in) :: bt
    integer(kind=int32), intent(in) :: l
    integer(kind=int32), intent(in), optional :: lbwmaxl0,ubwmaxl0
    type(error_info), intent(inout), optional :: error
    
    integer(kind=int32) :: lbwmax, ubwmax, lbwmaxl, ubwmaxl, n, lbw, ubw
    type(routine_info), parameter :: info=info_z_bt_leading

    if (failure(error)) return
    call push_id(info,error)

    lbwmax=get_lbwmax(bt)
    ubwmax=get_ubwmax(bt)
    n=get_n(bt)
    lbw=bt%lbw
    ubw=bt%ubw
    if (n < 1) then
       call set_error(1, info, error); return
    end if

    lbwmaxl=equals_option(lbw, lbwmaxl0)
    ubwmaxl=equals_option(ubw, ubwmaxl0)
    if (lbwmaxl < lbw .or. ubwmaxl < ubw) then
       call set_error(2, info, error); return
    end if

    btl = z_new_bt(l,lbwmaxl,ubwmaxl)

    btl%br(:,1:lbw+ubw+1)=bt%br(1:l,1:lbw+ubw+1)
    btl%lbw=lbw
    btl%ubw=ubw
    btl%kst(1:l-1,1:lbw)=bt%kst(1:l-1,1:lbw)
    btl%cst(1:l-1,1:lbw)=bt%cst(1:l-1,1:lbw)
    btl%sst(1:l-1,1:lbw)=bt%sst(1:l-1,1:lbw)
    btl%numrotst=bt%numrotst(1:l-1)
    call pop_id(error)
  end function z_bt_leading

  type(d_wb) function d_wb_trailing(wb,l,lbwmaxl0,ubwmaxl0,error) result(wbl)
    type(d_wb), intent(in) :: wb
    integer(kind=int32), intent(in) :: l
    integer(kind=int32), intent(in), optional :: lbwmaxl0,ubwmaxl0
    type(error_info), intent(inout), optional :: error
    
    integer(kind=int32) :: lbwmax, ubwmax, lbwmaxl, ubwmaxl, n, lbw, ubw, j, k
    type(routine_info), parameter :: info=info_d_wb_trailing

    if (failure(error)) return
    call push_id(info,error)

    lbwmax=get_lbwmax(wb)
    ubwmax=get_ubwmax(wb)
    n=get_n(wb)
    lbw=wb%lbw
    ubw=wb%ubw
    if (n < 1) then
       call set_error(1, info, error); return
    end if

    lbwmaxl=equals_option(lbw, lbwmaxl0)
    ubwmaxl=equals_option(ubw, ubwmaxl0)
    if (lbwmaxl < lbw .or. ubwmaxl < ubw) then
       call set_error(2, info, error); return
    end if

    wbl = d_new_wb(l,lbwmaxl,ubwmaxl)

    wbl%bc(1:lbw+ubw+1,:)=wb%bc(1:lbw+ubw+1,n-l+1:n)
    wbl%lbw=lbw
    wbl%ubw=ubw
    wbl%jsw(1:lbw,:)=wb%jsw(1:lbw,n-l+1:n)
    wbl%csw(1:lbw,:)=wb%csw(1:lbw,n-l+1:n)
    wbl%ssw(1:lbw,:)=wb%ssw(1:lbw,n-l+1:n)
    wbl%numrotsw=wb%numrotsw(n-l+1:n)
    do k=1,l
       do j=1,wbl%numrotsw(k)
          wbl%jsw(j,k)=wbl%jsw(j,k)-(n-l)
       end do
    end do
    call pop_id(error)
  end function d_wb_trailing  

  type(z_wb) function z_wb_trailing(wb,l,lbwmaxl0,ubwmaxl0,error) result(wbl)
    type(z_wb), intent(in) :: wb
    integer(kind=int32), intent(in) :: l
    integer(kind=int32), intent(in), optional :: lbwmaxl0,ubwmaxl0
    type(error_info), intent(inout), optional :: error
    
    integer(kind=int32) :: lbwmax, ubwmax, lbwmaxl, ubwmaxl, n, lbw, ubw, j, k
    type(routine_info), parameter :: info=info_z_wb_trailing

    if (failure(error)) return
    call push_id(info,error)

    lbwmax=get_lbwmax(wb)
    ubwmax=get_ubwmax(wb)
    n=get_n(wb)
    lbw=wb%lbw
    ubw=wb%ubw
    if (n < 1) then
       call set_error(1, info, error); return
    end if

    lbwmaxl=equals_option(lbw, lbwmaxl0)
    ubwmaxl=equals_option(ubw, ubwmaxl0)
    if (lbwmaxl < lbw .or. ubwmaxl < ubw) then
       call set_error(2, info, error); return
    end if

    wbl = z_new_wb(l,lbwmaxl,ubwmaxl)

    wbl%bc(1:lbw+ubw+1,:)=wb%bc(1:lbw+ubw+1,n-l+1:n)
    wbl%lbw=lbw
    wbl%ubw=ubw
    wbl%jsw(1:lbw,:)=wb%jsw(1:lbw,n-l+1:n)
    wbl%csw(1:lbw,:)=wb%csw(1:lbw,n-l+1:n)
    wbl%ssw(1:lbw,:)=wb%ssw(1:lbw,n-l+1:n)
    wbl%numrotsw=wb%numrotsw(n-l+1:n)
    do k=1,l
       do j=1,wbl%numrotsw(k)
          wbl%jsw(j,k)=wbl%jsw(j,k)-(n-l)
       end do
    end do
    call pop_id(error)
  end function z_wb_trailing

  type(d_ubt) function d_ubt_leading(ubt,l,lbwmaxl0,ubwmaxl0,error) result(ubtl)
    type(d_ubt), intent(in) :: ubt
    integer(kind=int32), intent(in) :: l
    integer(kind=int32), intent(in), optional :: lbwmaxl0,ubwmaxl0
    type(error_info), intent(inout), optional :: error
    
    integer(kind=int32) :: lbwmax, ubwmax, lbwmaxl, ubwmaxl, n, lbw, ubw
    type(routine_info), parameter :: info=info_d_ubt_leading

    if (failure(error)) return
    call push_id(info,error)

    lbwmax=get_lbwmax(ubt)
    ubwmax=get_ubwmax(ubt)
    n=get_n(ubt)
    lbw=ubt%lbw
    ubw=ubt%ubw
    if (n < 1) then
       call set_error(1, info, error); return
    end if

    lbwmaxl=equals_option(lbw, lbwmaxl0)
    ubwmaxl=equals_option(ubw, ubwmaxl0)
    if (lbwmaxl < lbw .or. ubwmaxl < ubw) then
       call set_error(2, info, error); return
    end if

    ubtl = d_new_ubt(l,lbwmaxl,ubwmaxl)

    ubtl%bc(1:lbw+ubw+1,:)=ubt%bc(1:lbw+ubw+1,1:l)
    ubtl%lbw=lbw
    ubtl%ubw=ubw
    ubtl%jsu(1:ubw,1:l-1)=ubt%jsu(1:ubw,1:l-1)
    ubtl%csu(1:ubw,1:l-1)=ubt%csu(1:ubw,1:l-1)
    ubtl%ssu(1:ubw,1:l-1)=ubt%ssu(1:ubw,1:l-1)
    ubtl%numrotsu=ubt%numrotsu(1:l-1)

    ubtl%kst(1:l-1,1:lbw)=ubt%kst(1:l-1,1:lbw)
    ubtl%cst(1:l-1,1:lbw)=ubt%cst(1:l-1,1:lbw)
    ubtl%sst(1:l-1,1:lbw)=ubt%sst(1:l-1,1:lbw)
    ubtl%numrotst=ubt%numrotst(1:l-1)

    call pop_id(error)
  end function d_ubt_leading

  type(z_ubt) function z_ubt_leading(ubt,l,lbwmaxl0,ubwmaxl0,error) result(ubtl)
    type(z_ubt), intent(in) :: ubt
    integer(kind=int32), intent(in) :: l
    integer(kind=int32), intent(in), optional :: lbwmaxl0,ubwmaxl0
    type(error_info), intent(inout), optional :: error
    
    integer(kind=int32) :: lbwmax, ubwmax, lbwmaxl, ubwmaxl, n, lbw, ubw
    type(routine_info), parameter :: info=info_z_ubt_leading

    if (failure(error)) return
    call push_id(info,error)

    lbwmax=get_lbwmax(ubt)
    ubwmax=get_ubwmax(ubt)
    n=get_n(ubt)
    lbw=ubt%lbw
    ubw=ubt%ubw
    if (n < 1) then
       call set_error(1, info, error); return
    end if

    lbwmaxl=equals_option(lbw, lbwmaxl0)
    ubwmaxl=equals_option(ubw, ubwmaxl0)
    if (lbwmaxl < lbw .or. ubwmaxl < ubw) then
       call set_error(2, info, error); return
    end if

    ubtl = z_new_ubt(l,lbwmaxl,ubwmaxl)

    ubtl%bc(1:lbw+ubw+1,:)=ubt%bc(1:lbw+ubw+1,1:l)
    ubtl%lbw=lbw
    ubtl%ubw=ubw
    ubtl%jsu(1:ubw,1:l-1)=ubt%jsu(1:ubw,1:l-1)
    ubtl%csu(1:ubw,1:l-1)=ubt%csu(1:ubw,1:l-1)
    ubtl%ssu(1:ubw,1:l-1)=ubt%ssu(1:ubw,1:l-1)
    ubtl%numrotsu=ubt%numrotsu(1:l-1)

    ubtl%kst(1:l-1,1:lbw)=ubt%kst(1:l-1,1:lbw)
    ubtl%cst(1:l-1,1:lbw)=ubt%cst(1:l-1,1:lbw)
    ubtl%sst(1:l-1,1:lbw)=ubt%sst(1:l-1,1:lbw)
    ubtl%numrotst=ubt%numrotst(1:l-1)

    call pop_id(error)
  end function z_ubt_leading

  type(d_wbv) function d_wbv_trailing(wbv,l,lbwmaxl0,ubwmaxl0,error) result(wbvl)
    type(d_wbv), intent(in) :: wbv
    integer(kind=int32), intent(in) :: l
    integer(kind=int32), intent(in), optional :: lbwmaxl0,ubwmaxl0
    type(error_info), intent(inout), optional :: error
    
    integer(kind=int32) :: lbwmax, ubwmax, lbwmaxl, ubwmaxl, n, lbw, ubw, j, k
    type(routine_info), parameter :: info=info_d_wbv_trailing

    if (failure(error)) return
    call push_id(info,error)

    lbwmax=get_lbwmax(wbv)
    ubwmax=get_ubwmax(wbv)
    n=get_n(wbv)
    lbw=wbv%lbw
    ubw=wbv%ubw
    if (n < 1) then
       call set_error(1, info, error); return
    end if

    lbwmaxl=equals_option(lbw, lbwmaxl0)
    ubwmaxl=equals_option(ubw, ubwmaxl0)
    if (lbwmaxl < lbw .or. ubwmaxl < ubw) then
       call set_error(2, info, error); return
    end if

    wbvl = d_new_wbv(l,lbwmaxl,ubwmaxl)

    wbvl%br(:,1:lbw+ubw+1)=wbv%br(n-l+1:n,1:lbw+ubw+1)
    wbvl%lbw=lbw
    wbvl%ubw=ubw

    wbvl%ksv(:,1:ubw)=wbv%ksv(n-l+1:n,1:ubw)
    wbvl%csv(:,1:ubw)=wbv%csv(n-l+1:n,1:ubw)
    wbvl%ssv(:,1:ubw)=wbv%ssv(n-l+1:n,1:ubw)
    wbvl%numrotsv=wbv%numrotsv(n-l+1:n)

    do j=1,l
       do k=1,wbvl%numrotsv(j)
          wbvl%ksv(j,k)=wbvl%ksv(j,k)-(n-l)
       end do
    end do
    call pop_id(error)

    wbvl%jsw(1:lbw,:)=wbv%jsw(1:lbw,n-l+1:n)
    wbvl%csw(1:lbw,:)=wbv%csw(1:lbw,n-l+1:n)
    wbvl%ssw(1:lbw,:)=wbv%ssw(1:lbw,n-l+1:n)
    wbvl%numrotsw=wbv%numrotsw(n-l+1:n)
    do k=1,l
       do j=1,wbvl%numrotsw(k)
          wbvl%jsw(j,k)=wbvl%jsw(j,k)-(n-l)
       end do
    end do
    
  end function d_wbv_trailing  

  type(z_wbv) function z_wbv_trailing(wbv,l,lbwmaxl0,ubwmaxl0,error) result(wbvl)
    type(z_wbv), intent(in) :: wbv
    integer(kind=int32), intent(in) :: l
    integer(kind=int32), intent(in), optional :: lbwmaxl0,ubwmaxl0
    type(error_info), intent(inout), optional :: error
    
    integer(kind=int32) :: lbwmax, ubwmax, lbwmaxl, ubwmaxl, n, lbw, ubw, j, k
    type(routine_info), parameter :: info=info_z_wbv_trailing

    if (failure(error)) return
    call push_id(info,error)

    lbwmax=get_lbwmax(wbv)
    ubwmax=get_ubwmax(wbv)
    n=get_n(wbv)
    lbw=wbv%lbw
    ubw=wbv%ubw
    if (n < 1) then
       call set_error(1, info, error); return
    end if

    lbwmaxl=equals_option(lbw, lbwmaxl0)
    ubwmaxl=equals_option(ubw, ubwmaxl0)
    if (lbwmaxl < lbw .or. ubwmaxl < ubw) then
       call set_error(2, info, error); return
    end if

    wbvl = z_new_wbv(l,lbwmaxl,ubwmaxl)

    wbvl%br(:,1:lbw+ubw+1)=wbv%br(n-l+1:n,1:lbw+ubw+1)
    wbvl%lbw=lbw
    wbvl%ubw=ubw

    wbvl%ksv(:,1:ubw)=wbv%ksv(n-l+1:n,1:ubw)
    wbvl%csv(:,1:ubw)=wbv%csv(n-l+1:n,1:ubw)
    wbvl%ssv(:,1:ubw)=wbv%ssv(n-l+1:n,1:ubw)
    wbvl%numrotsv=wbv%numrotsv(n-l+1:n)

    do j=1,l
       do k=1,wbvl%numrotsv(j)
          wbvl%ksv(j,k)=wbvl%ksv(j,k)-(n-l)
       end do
    end do
    call pop_id(error)

    wbvl%jsw(1:lbw,:)=wbv%jsw(1:lbw,n-l+1:n)
    wbvl%csw(1:lbw,:)=wbv%csw(1:lbw,n-l+1:n)
    wbvl%ssw(1:lbw,:)=wbv%ssw(1:lbw,n-l+1:n)
    wbvl%numrotsw=wbv%numrotsw(n-l+1:n)
    do k=1,l
       do j=1,wbvl%numrotsw(k)
          wbvl%jsw(j,k)=wbvl%jsw(j,k)-(n-l)
       end do
    end do
    
  end function z_wbv_trailing

  ! Unstructured submatrices.
  
  ! general submatrices

  function d_columns_of_ub(k0,k1,ub, error) result(a)
    real(kind=dp), dimension(:,:), allocatable :: a
    type(d_ub), intent(in) :: ub
    type(error_info), intent(inout), optional :: error
    integer(kind=int32), intent(in) :: k0,k1
    type(routine_info), parameter :: info=info_d_columns_of_ub
    integer(kind=int32) :: n

    if (failure(error)) return
    call push_id(info, error)
    
    n=get_n(ub)
    allocate(a(n,k1-k0+1))
    call d_ub_to_columns(k0,k1,ub,a,error)
    call pop_id(error)
  end function d_columns_of_ub

  ! Errors
  ! 0: no error
  ! 1: Size error.
  subroutine d_ub_to_columns(k0,k1,ub,a,error)
    type(d_ub), intent(in) :: ub
    real(kind=dp), dimension(:,:), intent(out) :: a
    integer(kind=int32), intent(in) :: k0,k1
    type(error_info), intent(inout), optional :: error
    type(routine_info), parameter :: info=info_d_ub_to_columns

    if (failure(error)) return
    call push_id(info, error)

    if (get_n(ub) /= size(a,1) .or. (k1-k0+1) /= size(a,2) .or. k1<k0 .or. &
         k1 > get_n(ub) .or. k0 < 1) then
       call set_error(1, info, error); return
    end if
    call f_d_ub_to_columns(k0,k1,ub%bc, get_n(ub), ub%lbw, ub%ubw, get_lbwmax(ub), &
         get_ubwmax(ub), ub%numrotsu, ub%jsu, ub%csu, ub%ssu, a)
    call pop_id(error)
  end subroutine d_ub_to_columns

  subroutine f_d_ub_to_columns(k0,k1, bc, n, lbw, ubw, lbwmax, ubwmax, numrotsu, jsu, csu, ssu, a)
    real(kind=dp), target, dimension(n,k1-k0+1), intent(out) :: a
    integer(kind=int32), dimension(ubwmax,n), intent(in) :: jsu
    real(kind=dp), dimension(ubwmax,n), intent(in) :: csu, ssu
    real(kind=dp), dimension(lbwmax+ubwmax+1,n), intent(in) :: bc
    integer(kind=int32), dimension(n), intent(in) :: numrotsu
    integer(kind=int32), intent(in) :: k0, k1, ubw, lbw, n, lbwmax, ubwmax
    !
    integer(kind=int32) :: j,k,d,ka0,ka1
    type(d_rotation) :: rot

    a=0.0_dp
    do d=1,ubw+1
       do k=max(k0,ubw-d+2),k1          
          a(d+k-ubw-1,k-k0+1)=bc(d,k)
       end do
    end do
    do d=ubw+2, ubw+lbw+1
       do k=k0,min(k1,n-d+ubw+1)
          a(d+k-ubw-1,k-k0+1)=bc(d,k)
       end do
    end do

    if (n==1) return

    ka1=k1-k0+1
    do k=k1-1,2,-1
       ka0=max(1,k+1-k0+1)
       do j=1,numrotsu(k)
          rot%cosine=csu(j,k); rot%sine=ssu(j,k)
          call rotation_times_general(rot,a(:,ka0:ka1),jsu(j,k),jsu(j,k)+1)
       end do
    end do
  end subroutine f_d_ub_to_columns


  function z_columns_of_ub(k0,k1,ub, error) result(a)
    complex(kind=dp), dimension(:,:), allocatable :: a
    type(z_ub), intent(in) :: ub
    type(error_info), intent(inout), optional :: error
    integer(kind=int32), intent(in) :: k0,k1
    type(routine_info), parameter :: info=info_z_columns_of_ub
    integer(kind=int32) :: n

    if (failure(error)) return
    call push_id(info, error)
    
    n=get_n(ub)
    allocate(a(n,k1-k0+1))
    call z_ub_to_columns(k0,k1,ub,a,error)
    call pop_id(error)
  end function z_columns_of_ub

  ! Errors
  ! 0: no error
  ! 1: Size error.
  subroutine z_ub_to_columns(k0,k1,ub,a,error)
    type(z_ub), intent(in) :: ub
    complex(kind=dp), dimension(:,:), intent(out) :: a
    integer(kind=int32), intent(in) :: k0,k1
    type(error_info), intent(inout), optional :: error
    type(routine_info), parameter :: info=info_z_ub_to_columns

    if (failure(error)) return
    call push_id(info, error)

    if (get_n(ub) /= size(a,1) .or. (k1-k0+1) /= size(a,2) .or. k1<k0 .or. &
         k1 > get_n(ub) .or. k0 < 1) then
       call set_error(1, info, error); return
    end if
    call f_z_ub_to_columns(k0,k1,ub%bc, get_n(ub), ub%lbw, ub%ubw, get_lbwmax(ub), &
         get_ubwmax(ub), ub%numrotsu, ub%jsu, ub%csu, ub%ssu, a)
    call pop_id(error)
  end subroutine z_ub_to_columns

  subroutine f_z_ub_to_columns(k0,k1, bc, n, lbw, ubw, lbwmax, ubwmax, numrotsu, jsu, csu, ssu, a)
    complex(kind=dp), target, dimension(n,k1-k0+1), intent(out) :: a
    integer(kind=int32), dimension(ubwmax,n), intent(in) :: jsu
    real(kind=dp), dimension(ubwmax,n), intent(in) :: csu
    complex(kind=dp), dimension(ubwmax,n), intent(in) :: ssu    
    complex(kind=dp), dimension(lbwmax+ubwmax+1,n), intent(in) :: bc
    integer(kind=int32), dimension(n), intent(in) :: numrotsu
    integer(kind=int32), intent(in) :: k0, k1, ubw, lbw, n, lbwmax, ubwmax
    !
    integer(kind=int32) :: j,k,d,ka0,ka1
    type(z_rotation) :: rot

    a=(0.0_dp,0.0_dp)
    do d=1,ubw+1
       do k=max(k0,ubw-d+2),k1          
          a(d+k-ubw-1,k-k0+1)=bc(d,k)
       end do
    end do
    do d=ubw+2, ubw+lbw+1
       do k=k0,min(k1,n-d+ubw+1)
          a(d+k-ubw-1,k-k0+1)=bc(d,k)
       end do
    end do

    if (n==1) return

    ka1=k1-k0+1
    do k=k1-1,2,-1
       ka0=max(1,k+1-k0+1)
       do j=1,numrotsu(k)
          rot%cosine=csu(j,k); rot%sine=ssu(j,k)
          call rotation_times_general(rot,a(:,ka0:ka1),jsu(j,k),jsu(j,k)+1)
       end do
    end do
  end subroutine f_z_ub_to_columns

  function d_rows_of_bv(j0,j1,bv, error) result(a)
    real(kind=dp), dimension(:,:), allocatable :: a
    type(d_bv), intent(in) :: bv
    integer(kind=int32), intent(in) :: j0,j1
    type(error_info), intent(inout), optional :: error
    type(routine_info), parameter :: info=info_d_rows_of_bv
    integer(kind=int32) :: n

    if (failure(error)) return
    call push_id(info, error)
    
    n=get_n(bv)
    allocate(a(j1-j0+1,n))
    call d_bv_to_rows(j0,j1,bv,a,error)
    call pop_id(error)
  end function d_rows_of_bv


  subroutine d_bv_to_rows(j0,j1,bv,a,error)
    type(d_bv), intent(in) :: bv
    real(kind=dp), dimension(:,:), intent(out) :: a
    integer(kind=int32), intent(in) :: j0,j1
    type(error_info), intent(inout), optional :: error
    type(routine_info), parameter :: info=info_d_bv_to_rows
    
    if (failure(error)) return
    call push_id(info, error)
    if (get_n(bv) /= size(a,2) .or. (j1-j0+1) /= size(a,1) .or. j0 > j1 .or. &
         j0 < 1 .or. j1 > get_n(bv)) then
       call set_error(1, info, error); return
    end if
    call f_d_bv_to_rows(j0,j1,bv%br, get_n(bv), bv%lbw, bv%ubw, get_lbwmax(bv), get_ubwmax(bv), &
         bv%numrotsv, bv%ksv, bv%csv, bv%ssv, a)
    call pop_id(error)
  end subroutine d_bv_to_rows

  subroutine f_d_bv_to_rows(j0,j1,br, n, lbw, ubw, lbwmax, ubwmax, numrotsv, ksv, csv, ssv, a)
    real(kind=dp), target, dimension(j1-j0+1,n), intent(out) :: a
    integer(kind=int32), dimension(n,ubwmax), intent(in) :: ksv
    real(kind=dp), dimension(n, ubwmax), intent(in) :: csv, ssv
    real(kind=dp), dimension(n,lbwmax+ubwmax+1), intent(in) :: br
    integer(kind=int32), dimension(n), intent(in) :: numrotsv
    integer(kind=int32), intent(in) :: ubw, lbw, n, lbwmax, ubwmax, j0, j1
    !
    integer(kind=int32) :: j,k,d,ja1
    type(d_rotation) :: rot

    a=0.0_dp
    do d=1,lbw+1
       do j=max(j0,lbw-d+2), j1       
          a(j-j0+1,d+j-lbw-1)=br(j,d)
       end do
    end do
    do d=lbw+2,ubw+lbw+1
       do j=j0,min(j1,n-d+lbw+1)
          a(j-j0+1,d+j-lbw-1)=br(j,d)
       end do
    end do

    if (n==1) return

    do j=j0,n-2
       ja1=min(j1-j0+1,j-j0+1)
       do k=1,numrotsv(j)
          rot%cosine=csv(j,k); rot%sine=ssv(j,k)
          call general_times_rotation(a(1:ja1,:),trp_rot(rot),ksv(j,k), ksv(j,k)+1)
       end do
    end do

  end subroutine f_d_bv_to_rows

  function z_rows_of_bv(j0,j1,bv, error) result(a)
    complex(kind=dp), dimension(:,:), allocatable :: a
    type(z_bv), intent(in) :: bv
    integer(kind=int32), intent(in) :: j0,j1
    type(error_info), intent(inout), optional :: error
    type(routine_info), parameter :: info=info_z_rows_of_bv
    integer(kind=int32) :: n

    if (failure(error)) return
    call push_id(info, error)
    
    n=get_n(bv)
    allocate(a(j1-j0+1,n))
    call z_bv_to_rows(j0,j1,bv,a,error)
    call pop_id(error)
  end function z_rows_of_bv


  subroutine z_bv_to_rows(j0,j1,bv,a,error)
    type(z_bv), intent(in) :: bv
    complex(kind=dp), dimension(:,:), intent(out) :: a
    integer(kind=int32), intent(in) :: j0,j1
    type(error_info), intent(inout), optional :: error
    type(routine_info), parameter :: info=info_z_bv_to_rows
    
    if (failure(error)) return
    call push_id(info, error)
    if (get_n(bv) /= size(a,2) .or. (j1-j0+1) /= size(a,1) .or. j0 > j1 .or. &
         j0 < 1 .or. j1 > get_n(bv)) then
       call set_error(1, info, error); return
    end if
    call f_z_bv_to_rows(j0,j1,bv%br, get_n(bv), bv%lbw, bv%ubw, get_lbwmax(bv), get_ubwmax(bv), &
         bv%numrotsv, bv%ksv, bv%csv, bv%ssv, a)
    call pop_id(error)
  end subroutine z_bv_to_rows

  subroutine f_z_bv_to_rows(j0,j1,br, n, lbw, ubw, lbwmax, ubwmax, numrotsv, ksv, csv, ssv, a)
    complex(kind=dp), target, dimension(j1-j0+1,n), intent(out) :: a
    integer(kind=int32), dimension(n,ubwmax), intent(in) :: ksv
    real(kind=dp), dimension(n, ubwmax), intent(in) :: csv
    complex(kind=dp), dimension(n, ubwmax), intent(in) :: ssv    
    complex(kind=dp), dimension(n,lbwmax+ubwmax+1), intent(in) :: br
    integer(kind=int32), dimension(n), intent(in) :: numrotsv
    integer(kind=int32), intent(in) :: ubw, lbw, n, lbwmax, ubwmax, j0, j1
    !
    integer(kind=int32) :: j,k,d,ja1
    type(z_rotation) :: rot

    a=(0.0_dp,0.0_dp)
    do d=1,lbw+1
       do j=max(j0,lbw-d+2), j1       
          a(j-j0+1,d+j-lbw-1)=br(j,d)
       end do
    end do
    do d=lbw+2,ubw+lbw+1
       do j=j0,min(j1,n-d+lbw+1)
          a(j-j0+1,d+j-lbw-1)=br(j,d)
       end do
    end do

    if (n==1) return

    do j=j0,n-2
       ja1=min(j1-j0+1,j-j0+1)
       do k=1,numrotsv(j)
          rot%cosine=csv(j,k); rot%sine=ssv(j,k)
          call general_times_rotation(a(1:ja1,:),trp_rot(rot),ksv(j,k), ksv(j,k)+1)
       end do
    end do

  end subroutine f_z_bv_to_rows

  function d_columns_of_wb(k0,k1,wb, error) result(a)
    real(kind=dp), dimension(:,:), allocatable :: a
    type(d_wb), intent(in) :: wb
    integer(kind=int32), intent(in) :: k0,k1
    type(error_info), intent(inout), optional :: error
    type(routine_info), parameter :: info=info_d_columns_of_wb
    integer(kind=int32) :: n

    if (failure(error)) return
    call push_id(info, error)
    
    n=get_n(wb)
    allocate(a(n,k1-k0+1))
    call d_wb_to_columns(k0,k1,wb,a,error)
    call pop_id(error)
  end function d_columns_of_wb

  subroutine d_wb_to_columns(k0,k1,wb,a,error)
    type(d_wb), intent(in) :: wb
    real(kind=dp), dimension(:,:), intent(out) :: a
    integer(kind=int32), intent(in) :: k0,k1
    type(error_info), intent(inout), optional :: error
    type(routine_info), parameter :: info=info_d_wb_to_columns

    if (failure(error)) return
    call push_id(info, error)
    if (get_n(wb) /= size(a,1) .or. (k1-k0+1) /= size(a,2) .or. k1<k0 .or. &
         k1 > get_n(wb) .or. k0 < 1) then
       call set_error(1, info, error); return
    end if
    call f_d_wb_to_columns(k0,k1,wb%bc, get_n(wb), wb%lbw, wb%ubw, get_lbwmax(wb), get_ubwmax(wb), &
         wb%numrotsw, wb%jsw, wb%csw, wb%ssw, a)
    call pop_id(error)
  end subroutine d_wb_to_columns

  subroutine f_d_wb_to_columns(k0,k1,bc, n, lbw, ubw, lbwmax, ubwmax, numrotsw, jsw, csw, ssw, a)
    real(kind=dp), target, dimension(n,k1-k0+1), intent(out) :: a
    integer(kind=int32), dimension(lbwmax,n), intent(in) :: jsw
    real(kind=dp), dimension(lbwmax,n), intent(in) :: csw, ssw
    real(kind=dp), dimension(lbwmax+ubwmax+1,n), intent(in) :: bc
    integer(kind=int32), dimension(n), intent(in) :: numrotsw
    integer(kind=int32), intent(in) :: ubw, lbw, n, lbwmax, ubwmax, k0, k1
    !
    integer(kind=int32) :: j,k,d,ka1
    type(d_rotation) :: rot

    a=0.0_dp
    do d=1,ubw+1
       do k=max(k0,ubw-d+2),k1          
          a(d+k-ubw-1,k-k0+1)=bc(d,k)
       end do
    end do
    do d=ubw+2, ubw+lbw+1
       do k=k0,min(k1,n-d+ubw+1)
          a(d+k-ubw-1,k-k0+1)=bc(d,k)
       end do
    end do

    if (n==1) return

    do k=k0,n-2
       ka1=min(k-k0+1,k1-k0+1)
       do j=1,numrotsw(k)
          rot%cosine=csw(j,k); rot%sine=ssw(j,k)
          call rotation_times_general(rot,a(:,1:ka1),jsw(j,k),jsw(j,k)+1)
       end do
    end do
  end subroutine f_d_wb_to_columns

  function z_columns_of_wb(k0,k1,wb, error) result(a)
    complex(kind=dp), dimension(:,:), allocatable :: a
    type(z_wb), intent(in) :: wb
    integer(kind=int32), intent(in) :: k0,k1
    type(error_info), intent(inout), optional :: error
    type(routine_info), parameter :: info=info_z_columns_of_wb
    integer(kind=int32) :: n

    if (failure(error)) return
    call push_id(info, error)
    
    n=get_n(wb)
    allocate(a(n,k1-k0+1))
    call z_wb_to_columns(k0,k1,wb,a,error)
    call pop_id(error)
  end function z_columns_of_wb

  subroutine z_wb_to_columns(k0,k1,wb,a,error)
    type(z_wb), intent(in) :: wb
    complex(kind=dp), dimension(:,:), intent(out) :: a
    integer(kind=int32), intent(in) :: k0,k1
    type(error_info), intent(inout), optional :: error
    type(routine_info), parameter :: info=info_z_wb_to_columns

    if (failure(error)) return
    call push_id(info, error)
    if (get_n(wb) /= size(a,1) .or. (k1-k0+1) /= size(a,2) .or. k1<k0 .or. &
         k1 > get_n(wb) .or. k0 < 1) then
       call set_error(1, info, error); return
    end if
    call f_z_wb_to_columns(k0,k1,wb%bc, get_n(wb), wb%lbw, wb%ubw, get_lbwmax(wb), get_ubwmax(wb), &
         wb%numrotsw, wb%jsw, wb%csw, wb%ssw, a)
    call pop_id(error)
  end subroutine z_wb_to_columns

  subroutine f_z_wb_to_columns(k0,k1,bc, n, lbw, ubw, lbwmax, ubwmax, numrotsw, jsw, csw, ssw, a)
    complex(kind=dp), target, dimension(n,k1-k0+1), intent(out) :: a
    integer(kind=int32), dimension(lbwmax,n), intent(in) :: jsw
    real(kind=dp), dimension(lbwmax,n), intent(in) :: csw
    complex(kind=dp), dimension(lbwmax,n), intent(in) :: ssw    
    complex(kind=dp), dimension(lbwmax+ubwmax+1,n), intent(in) :: bc
    integer(kind=int32), dimension(n), intent(in) :: numrotsw
    integer(kind=int32), intent(in) :: ubw, lbw, n, lbwmax, ubwmax, k0, k1
    !
    integer(kind=int32) :: j,k,d,ka1
    type(z_rotation) :: rot

    a=(0.0_dp,0.0_dp)
    do d=1,ubw+1
       do k=max(k0,ubw-d+2),k1          
          a(d+k-ubw-1,k-k0+1)=bc(d,k)
       end do
    end do
    do d=ubw+2, ubw+lbw+1
       do k=k0,min(k1,n-d+ubw+1)
          a(d+k-ubw-1,k-k0+1)=bc(d,k)
       end do
    end do

    if (n==1) return

    do k=k0,n-2
       ka1=min(k-k0+1,k1-k0+1)
       do j=1,numrotsw(k)
          rot%cosine=csw(j,k); rot%sine=ssw(j,k)
          call rotation_times_general(rot,a(:,1:ka1),jsw(j,k),jsw(j,k)+1)
       end do
    end do
  end subroutine f_z_wb_to_columns

  function d_rows_of_bt(j0,j1,bt, error) result(a)
    real(kind=dp), dimension(:,:), allocatable :: a
    type(d_bt), intent(in) :: bt
    integer(kind=int32), intent(in) :: j0,j1
    type(error_info), intent(inout), optional :: error
    type(routine_info), parameter :: info=info_d_rows_of_bt
    integer(kind=int32) :: n

    if (failure(error)) return
    call push_id(info, error)
    
    n=get_n(bt)
    allocate(a(j1-j0+1,n))
    call d_bt_to_rows(j0,j1,bt,a,error)
    call pop_id(error)
  end function d_rows_of_bt

  ! Errors
  ! 0: no error
  ! 1: bt%n /= n
  subroutine d_bt_to_rows(j0,j1,bt,a,error)
    type(d_bt), intent(in) :: bt
    real(kind=dp), dimension(:,:), intent(out) :: a
    integer(kind=int32), intent(in) :: j0,j1
    type(error_info), intent(inout), optional :: error
    type(routine_info), parameter :: info=info_d_bt_to_rows

    if (failure(error)) return
    call push_id(info, error)
    
    if (get_n(bt) /= size(a,2) .or. j1-j0+1 /= size(a,1) .or. j1<j0 .or. &
         j0<1 .or. j1>get_n(bt)) then
       call set_error(1, info, error); return
    end if
    call f_d_bt_to_rows(j0,j1,bt%br, get_n(bt), bt%lbw, bt%ubw, get_lbwmax(bt), get_ubwmax(bt), &
         bt%numrotst, bt%kst, bt%cst, bt%sst, a)
    call pop_id(error)
  end subroutine d_bt_to_rows

  subroutine f_d_bt_to_rows(j0,j1,br, n, lbw, ubw, lbwmax, ubwmax, numrotst, kst, cst, sst, a)
    real(kind=dp), target, dimension(j1-j0+1,n), intent(out) :: a
    integer(kind=int32), dimension(n,lbwmax), intent(in) :: kst
    real(kind=dp), dimension(n,lbwmax), intent(in) :: cst, sst
    real(kind=dp), dimension(n,lbwmax+ubwmax+1), intent(in) :: br
    integer(kind=int32), dimension(n), intent(in) :: numrotst
    integer(kind=int32), intent(in) :: j0,j1    
    integer(kind=int32), intent(in) :: ubw, lbw, n, lbwmax, ubwmax
    !
    integer(kind=int32) :: j,k, d, ja0, ja1
    type(d_rotation) :: rot

    a=0.0_dp
    do d=1,lbw+1
       do j=max(j0,lbw-d+2), j1       
          a(j-j0+1,d+j-lbw-1)=br(j,d)
       end do
    end do
    do d=lbw+2,ubw+lbw+1
       do j=j0,min(j1,n-d+lbw+1)
          a(j-j0+1,d+j-lbw-1)=br(j,d)
       end do
    end do

    if (n==1) return

    ja1=j1-j0+1
    do k=j1-1,2,-1
       ja0=max(k+1-j0+1,1)
       do j=1,numrotst(k)
          rot%cosine=cst(k,j); rot%sine=sst(k,j)
          call general_times_rotation(a(ja0:ja1,:), trp_rot(rot), kst(k,j),kst(k,j)+1)
       end do
    end do

  end subroutine f_d_bt_to_rows

  function z_rows_of_bt(j0,j1,bt, error) result(a)
    complex(kind=dp), dimension(:,:), allocatable :: a
    type(z_bt), intent(in) :: bt
    integer(kind=int32), intent(in) :: j0,j1
    type(error_info), intent(inout), optional :: error
    type(routine_info), parameter :: info=info_z_rows_of_bt
    integer(kind=int32) :: n

    if (failure(error)) return
    call push_id(info, error)
    
    n=get_n(bt)
    allocate(a(j1-j0+1,n))
    call z_bt_to_rows(j0,j1,bt,a,error)
    call pop_id(error)
  end function z_rows_of_bt

  ! Errors
  ! 0: no error
  ! 1: bt%n /= n
  subroutine z_bt_to_rows(j0,j1,bt,a,error)
    type(z_bt), intent(in) :: bt
    complex(kind=dp), dimension(:,:), intent(out) :: a
    integer(kind=int32), intent(in) :: j0,j1
    type(error_info), intent(inout), optional :: error
    type(routine_info), parameter :: info=info_z_bt_to_rows

    if (failure(error)) return
    call push_id(info, error)
    
    if (get_n(bt) /= size(a,2) .or. j1-j0+1 /= size(a,1) .or. j1<j0 .or. &
         j0<1 .or. j1>get_n(bt)) then
       call set_error(1, info, error); return
    end if
    call f_z_bt_to_rows(j0,j1,bt%br, get_n(bt), bt%lbw, bt%ubw, get_lbwmax(bt), get_ubwmax(bt), &
         bt%numrotst, bt%kst, bt%cst, bt%sst, a)
    call pop_id(error)
  end subroutine z_bt_to_rows

  subroutine f_z_bt_to_rows(j0,j1,br, n, lbw, ubw, lbwmax, ubwmax, numrotst, kst, cst, sst, a)
    complex(kind=dp), target, dimension(j1-j0+1,n), intent(out) :: a
    integer(kind=int32), dimension(n,lbwmax), intent(in) :: kst
    real(kind=dp), dimension(n,lbwmax), intent(in) :: cst
    complex(kind=dp), dimension(n,lbwmax), intent(in) :: sst
    complex(kind=dp), dimension(n,lbwmax+ubwmax+1), intent(in) :: br
    integer(kind=int32), dimension(n), intent(in) :: numrotst
    integer(kind=int32), intent(in) :: j0,j1    
    integer(kind=int32), intent(in) :: ubw, lbw, n, lbwmax, ubwmax
    !
    integer(kind=int32) :: j,k, d, ja0, ja1
    type(z_rotation) :: rot

    a=(0.0_dp,0.0_dp)
    do d=1,lbw+1
       do j=max(j0,lbw-d+2), j1       
          a(j-j0+1,d+j-lbw-1)=br(j,d)
       end do
    end do
    do d=lbw+2,ubw+lbw+1
       do j=j0,min(j1,n-d+lbw+1)
          a(j-j0+1,d+j-lbw-1)=br(j,d)
       end do
    end do

    if (n==1) return

    ja1=j1-j0+1
    do k=j1-1,2,-1
       ja0=max(k+1-j0+1,1)
       do j=1,numrotst(k)
          rot%cosine=cst(k,j); rot%sine=sst(k,j)
          call general_times_rotation(a(ja0:ja1,:), trp_rot(rot), kst(k,j),kst(k,j)+1)
       end do
    end do

  end subroutine f_z_bt_to_rows
  
  
  
end module mod_submatrix
