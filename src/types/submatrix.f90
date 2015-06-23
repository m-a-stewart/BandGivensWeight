module mod_submatrix
  use mod_prec
  use mod_error_id
  use mod_band_types
  use mod_orth_band_types
  use mod_rotation
  use mod_utility
  implicit none
  ! This module contains routines for assembling a Givens weight
  ! parameterized matrix into an unstructured representation in
  ! general $n\times n$ array.

  private

  public :: leading, d_ub_leading, z_ub_leading, d_bt_leading, z_bt_leading, &
       d_ubt_leading, z_ubt_leading, &
       trailing, d_bv_trailing, z_bv_trailing, d_wb_trailing, z_wb_trailing, &
       d_wbv_trailing, z_wbv_trailing

  interface leading
     module procedure d_ub_leading, z_ub_leading, d_bt_leading, z_bt_leading, &
          d_ubt_leading, z_ubt_leading
  end interface leading

  interface trailing
     module procedure d_bv_trailing, z_bv_trailing, d_wb_trailing, z_wb_trailing, &
          d_wbv_trailing, z_wbv_trailing
  end interface trailing
  
contains

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
  
  
  
end module mod_submatrix
