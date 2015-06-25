module mod_convert_wbv_to_ubt
  use mod_prec
  use mod_error_id
  use mod_shift
  use mod_rotation
  use mod_orth_band_types
  use mod_band_types
  implicit none

  private

  public :: convert_wbv_to_ubt, d_convert_wbv_to_ubt, z_convert_wbv_to_ubt, &
       f_convert_wbv_to_ubt, f_d_convert_wbv_to_ubt, f_z_convert_wbv_to_ubt, &
       d_ubt_of_wbv, z_ubt_of_wbv, ubt_of

  interface convert_wbv_to_ubt
     module procedure d_convert_wbv_to_ubt, z_convert_wbv_to_ubt
  end interface convert_wbv_to_ubt

  interface f_convert_wbv_to_ubt
     module procedure f_d_convert_wbv_to_ubt, f_z_convert_wbv_to_ubt
  end interface f_convert_wbv_to_ubt

  interface ubt_of
     module procedure d_ubt_of_wbv, z_ubt_of_wbv
  end interface ubt_of

contains

  function d_ubt_of_wbv(wbv,error) result(ubt)
    type(d_ubt) :: ubt
    type(d_wbv), intent(in) :: wbv
    type(error_info), intent(inout), optional :: error

    type(d_wbv), allocatable :: wbv1
    type(routine_info), parameter :: info=info_d_ubt_of_wbv

    if (failure(error)) return
    call push_id(info, error)
    
    ubt=d_new_ubt(get_n(wbv), wbv%lbw, wbv%ubw)
    wbv1=d_new_wbv(get_n(wbv),wbv%lbw+1,  wbv%ubw+1)
    call copy(wbv1,wbv)
    call d_convert_wbv_to_ubt(wbv1,ubt,error)

    call pop_id(error)
  end function d_ubt_of_wbv


  ! Errors:
  ! 0: no error
  ! 1: n<1
  ! 2: Insufficient storage in wbv
  ! 3: Insufficient stroage in ubt
  ! 4: wbv%n /= ubt%n

  subroutine d_convert_wbv_to_ubt(wbv, ubt, error)
    type(d_wbv) :: wbv
    type(d_ubt) :: ubt
    type(error_info), intent(inout), optional :: error
    type(routine_info), parameter :: info=info_d_convert_wbv_to_ubt

    if (failure(error)) return
    call push_id(info, error)
    
    if (get_n(wbv) < 1) then
       call set_error(1, info, error); return
    end if
    if ((get_ubwmax(wbv) < wbv%ubw+1 .and. wbv%ubw < get_n(wbv)-1) .or. &
         (get_lbwmax(wbv) < wbv%lbw+1 .and. wbv%lbw < get_n(wbv)-1)) then
       call set_error(2, info, error); return
    end if
    if (get_lbwmax(ubt) < wbv%lbw .or. get_ubwmax(ubt) < wbv%ubw) then
       call set_error(3, info, error); return
    end if
    if (get_n(wbv) /= get_n(ubt)) then
       call set_error(4, info, error); return
    end if
    call f_d_convert_wbv_to_ubt(wbv%br, get_n(wbv), wbv%lbw, wbv%ubw, get_lbwmax(wbv), &
         get_ubwmax(wbv), wbv%numrotsw, wbv%jsw, wbv%csw, wbv%ssw, &
         wbv%numrotsv, wbv%ksv, wbv%csv, wbv%ssv, &
         ubt%bc, ubt%lbw, ubt%ubw, get_lbwmax(ubt), get_ubwmax(ubt), &
         ubt%numrotsu, ubt%jsu, ubt%csu, ubt%ssu, &
         ubt%numrotst, ubt%kst, ubt%cst, ubt%sst)
    call pop_id(error)
  end subroutine d_convert_wbv_to_ubt

  subroutine f_d_convert_wbv_to_ubt(b_wbv, n, lbw, ubw, lbwmax_wbv, ubwmax_wbv, numrotsw, &
       jsw, csw, ssw, numrotsv, ksv, csv, ssv, &
       b_ubt, lbw_ubt, ubw_ubt, lbwmax_ubt, ubwmax_ubt, &
       numrotsu, jsu, csu, ssu, numrotst, kst, cst, sst)
    real(kind=dp), dimension(n,lbwmax_wbv+ubwmax_wbv+1), intent(inout) :: b_wbv
    integer(kind=int32), intent(in) :: n, lbw, ubw, lbwmax_wbv, ubwmax_wbv, lbwmax_ubt, ubwmax_ubt
    integer(kind=int32), dimension(n), intent(in) :: numrotsw
    integer(kind=int32), dimension(lbwmax_wbv,n), intent(in) :: jsw
    real(kind=dp), dimension(lbwmax_wbv,n), intent(in) :: csw, ssw
    integer(kind=int32), dimension(n), intent(in) :: numrotsv
    integer(kind=int32), dimension(n,ubwmax_wbv), intent(in) :: ksv
    real(kind=dp), dimension(n,ubwmax_wbv), intent(in) :: csv, ssv

    real(kind=dp), dimension(lbwmax_ubt+ubwmax_ubt+1,n), intent(out) :: b_ubt
    integer(kind=int32), dimension(n), intent(out) :: numrotsu
    integer(kind=int32), dimension(ubwmax_ubt,n), intent(out) :: jsu
    real(kind=dp), dimension(ubwmax_ubt,n), intent(out) :: csu, ssu
    integer(kind=int32), dimension(n), intent(out) :: numrotst
    integer(kind=int32), dimension(n,lbwmax_ubt), intent(out) :: kst
    real(kind=dp), dimension(n,lbwmax_ubt), intent(out) :: cst, sst
    integer(kind=int32), intent(out) :: lbw_ubt, ubw_ubt

    integer(kind=int32) :: j, k, k0, k1, ubw1, lbw1
    type(d_rotation) :: rot
    logical :: full_lbw

    b_ubt(1:lbw+ubw+1,:)=0.0_dp; numrotsu=0
    ssu(1:ubw,:)=0.0_dp; csu(1:ubw,:)=0.0_dp
    jsu(1:ubw,:)=0
    numrotst=0
    sst(:,1:lbw)=0.0_dp; cst(:,1:lbw)=0.0_dp
    kst(:,1:lbw)=0
    lbw_ubt=lbw; ubw_ubt=ubw
    if (n == 1) then
       b_ubt(1,1)=b_wbv(1,1);
       lbw_ubt=0; ubw_ubt=0; return
    end if
    ! must allow for temporary fill-in
    if (lbw < n-1) then
       lbw1=lbw+1
       call shift(b_wbv,0,1)
       full_lbw=.false.
    else
       lbw1=lbw
       full_lbw=.true.
    end if
    if (ubw < n-1) then
       ubw1=ubw+1
       b_wbv(:,lbw1+ubw1+1)=0.0_dp
    else
       lbw1=lbw
    end if
    ! Do the lower triangular part.
    ! k is the size of the leading principal submatrix
    do k=1,n-2
       ! Apply W_k
       do j=1,numrotsw(k)
          rot%cosine=csw(j,k); rot%sine=ssw(j,k)
          call rotation_times_tbr(rot,b_wbv,n,lbw1,ubw1,0,n-k,jsw(j,k))
       end do
       ! columns in which nonzeros have been introduced into the extra subdiagonal
       k0=max(k-lbw+1,1)
       k1=min(k,n-lbw1)
       ! Apply T_{k+1}
       numrotst(k+1)=max(k1-k0+1,0)
       do j=k1,k0,-1
          rot=rgivens2(get_el_br(b_wbv,lbw1,j+lbw1,j), get_el_br(b_wbv,lbw1,j+lbw1,j+1))
          kst(k+1,j-k0+1)=j
          cst(k+1,j-k0+1)=rot%cosine; sst(k+1,j-k0+1)=rot%sine
          call tbr_times_rotation(b_wbv,n,lbw1,ubw1,k+1,0,rot,j)
       end do
    end do
    ! upper
    do k=1,n-2
       ! Apply V_k
       do j=1,numrotsv(k)
          rot%cosine=csv(k,j); rot%sine=ssv(k,j)
          call tbr_times_rotation(b_wbv,n,lbw1,ubw1,0,n-k,trp_rot(rot),ksv(k,j))
       end do
       !
       k0=max(k+2,ubw1+1)
       k1=min(k+ubw1,n)
       numrotsu(k+1)=max(k1-k0+1,0)
       do j=k1,k0,-1
          rot=lgivens2(get_el_br(b_wbv,lbw1,j-ubw1,j), get_el_br(b_wbv,lbw1,j-ubw1+1,j))
          jsu(j-k0+1,k+1)=j-ubw1
          csu(j-k0+1,k+1)=rot%cosine; ssu(j-k0+1,k+1)=rot%sine
          call rotation_times_tbr(trp_rot(rot),b_wbv,n,lbw1,ubw1,k+1,0,j-ubw1)
       end do
    end do
    if (.not. full_lbw) then
       call shift(b_wbv,0,-1)
    end if
    call br_to_bc(b_wbv,b_ubt,lbw,ubw)
  end subroutine f_d_convert_wbv_to_ubt

  function z_ubt_of_wbv(wbv,error) result(ubt)
    type(z_ubt) :: ubt
    type(z_wbv), intent(in) :: wbv
    type(error_info), intent(inout), optional :: error

    type(z_wbv), allocatable :: wbv1
    type(routine_info), parameter :: info=info_z_ubt_of_wbv

    if (failure(error)) return
    call push_id(info, error)
    
    ubt=z_new_ubt(get_n(wbv), wbv%lbw, wbv%ubw)
    wbv1=z_new_wbv(get_n(wbv),wbv%lbw+1,  wbv%ubw+1)
    call copy(wbv1,wbv)
    call z_convert_wbv_to_ubt(wbv1,ubt,error)

    call pop_id(error)
  end function z_ubt_of_wbv
  

  ! Errors:
  ! 0: no error
  ! 1: n<1
  ! 2: Insufficient storage in wbv
  ! 3: Insufficient stroage in ubt
  ! 4: wbv%n /= ubt%n

  subroutine z_convert_wbv_to_ubt(wbv, ubt, error)
    type(z_wbv) :: wbv
    type(z_ubt) :: ubt
    type(error_info), intent(inout), optional :: error
    type(routine_info), parameter :: info=info_z_convert_wbv_to_ubt

    if (failure(error)) return
    call push_id(info, error)

    if (get_n(wbv) < 1) then
       call set_error(1, info, error); return
    end if
    if ((get_ubwmax(wbv) < wbv%ubw+1 .and. wbv%ubw < get_n(wbv)-1) .or. &
         (get_lbwmax(wbv) < wbv%lbw+1 .and. wbv%lbw < get_n(wbv)-1)) then
       call set_error(2, info, error); return
    end if
    if (get_lbwmax(ubt) < wbv%lbw .or. get_ubwmax(ubt) < wbv%ubw) then
       call set_error(3, info, error); return
    end if
    if (get_n(wbv) /= get_n(ubt)) then
       call set_error(4, info, error); return
    end if
    call f_z_convert_wbv_to_ubt(wbv%br, get_n(wbv), wbv%lbw, wbv%ubw, get_lbwmax(wbv), &
         get_ubwmax(wbv), wbv%numrotsw, wbv%jsw, wbv%csw, wbv%ssw, &
         wbv%numrotsv, wbv%ksv, wbv%csv, wbv%ssv, &
         ubt%bc, ubt%lbw, ubt%ubw, get_lbwmax(ubt), get_ubwmax(ubt), &
         ubt%numrotsu, ubt%jsu, ubt%csu, ubt%ssu, &
         ubt%numrotst, ubt%kst, ubt%cst, ubt%sst)
    call pop_id(error)
  end subroutine z_convert_wbv_to_ubt

  subroutine f_z_convert_wbv_to_ubt(b_wbv, n, lbw, ubw, lbwmax_wbv, ubwmax_wbv, numrotsw, &
       jsw, csw, ssw, numrotsv, ksv, csv, ssv, &
       b_ubt, lbw_ubt, ubw_ubt, lbwmax_ubt, ubwmax_ubt, &
       numrotsu, jsu, csu, ssu, numrotst, kst, cst, sst)
    complex(kind=dp), dimension(n,lbwmax_wbv+ubwmax_wbv+1), intent(inout) :: b_wbv
    integer(kind=int32), intent(in) :: n, lbw, ubw, lbwmax_wbv, ubwmax_wbv, lbwmax_ubt, ubwmax_ubt
    integer(kind=int32), dimension(n), intent(in) :: numrotsw
    integer(kind=int32), dimension(lbwmax_wbv,n), intent(in) :: jsw
    real(kind=dp), dimension(lbwmax_wbv,n), intent(in) :: csw
    complex(kind=dp), dimension(lbwmax_wbv,n), intent(in) :: ssw
    integer(kind=int32), dimension(n), intent(in) :: numrotsv
    integer(kind=int32), dimension(n,ubwmax_wbv), intent(in) :: ksv
    real(kind=dp), dimension(n,ubwmax_wbv), intent(in) :: csv
    complex(kind=dp), dimension(n,ubwmax_wbv), intent(in) :: ssv

    complex(kind=dp), dimension(lbwmax_ubt+ubwmax_ubt+1,n), intent(out) :: b_ubt
    integer(kind=int32), dimension(n), intent(out) :: numrotsu
    integer(kind=int32), dimension(ubwmax_ubt,n), intent(out) :: jsu
    real(kind=dp), dimension(ubwmax_ubt,n), intent(out) :: csu
    complex(kind=dp), dimension(ubwmax_ubt,n), intent(out) :: ssu
    integer(kind=int32), dimension(n), intent(out) :: numrotst
    integer(kind=int32), dimension(n,lbwmax_ubt), intent(out) :: kst
    real(kind=dp), dimension(n,lbwmax_ubt), intent(out) :: cst
    complex(kind=dp), dimension(n,lbwmax_ubt), intent(out) :: sst
    integer(kind=int32), intent(out) :: lbw_ubt, ubw_ubt

    integer(kind=int32) :: j, k, k0, k1, ubw1, lbw1
    type(z_rotation) :: rot
    logical :: full_lbw

    b_ubt(1:lbw+ubw+1,:)=(0.0_dp,0.0_dp); numrotsu=0
    ssu(1:ubw,:)=(0.0_dp,0.0_dp); csu(1:ubw,:)=0.0_dp
    jsu(1:ubw,:)=0
    numrotst=0
    sst(:,1:lbw)=(0.0_dp,0.0_dp); cst(:,1:lbw)=0.0_dp
    kst(:,1:lbw)=0
    lbw_ubt=lbw; ubw_ubt=ubw
    if (n == 1) then
       b_ubt(1,1)=b_wbv(1,1);
       lbw_ubt=0; ubw_ubt=0; return
    end if
    ! must allow for temporary fill-in
    if (lbw < n-1) then
       lbw1=lbw+1
       call shift(b_wbv,0,1)
       full_lbw=.false.
    else
       lbw1=lbw
       full_lbw=.true.
    end if
    if (ubw < n-1) then
       ubw1=ubw+1
       b_wbv(:,lbw1+ubw1+1)=(0.0_dp,0.0_dp)
    else
       lbw1=lbw
    end if
    ! Do the lower triangular part.
    ! k is the size of the leading principal submatrix
    do k=1,n-2
       ! Apply W_k
       do j=1,numrotsw(k)
          rot%cosine=csw(j,k); rot%sine=ssw(j,k)
          call rotation_times_tbr(rot,b_wbv,n,lbw1,ubw1,0,n-k,jsw(j,k))
       end do
       ! columns in which nonzeros have been introduced into the extra subdiagonal
       k0=max(k-lbw+1,1)
       k1=min(k,n-lbw1)
       ! Apply T_{k+1}
       numrotst(k+1)=max(k1-k0+1,0)
       do j=k1,k0,-1
          rot=rgivens2(get_el_br(b_wbv,lbw1,j+lbw1,j), get_el_br(b_wbv,lbw1,j+lbw1,j+1))
          kst(k+1,j-k0+1)=j
          cst(k+1,j-k0+1)=rot%cosine; sst(k+1,j-k0+1)=rot%sine
          call tbr_times_rotation(b_wbv,n,lbw1,ubw1,k+1,0,rot,j)
       end do
    end do
    ! upper
    do k=1,n-2
       ! Apply V_k
       do j=1,numrotsv(k)
          rot%cosine=csv(k,j); rot%sine=ssv(k,j)
          call tbr_times_rotation(b_wbv,n,lbw1,ubw1,0,n-k,trp_rot(rot),ksv(k,j))
       end do
       !
       k0=max(k+2,ubw1+1)
       k1=min(k+ubw1,n)
       numrotsu(k+1)=max(k1-k0+1,0)
       do j=k1,k0,-1
          rot=lgivens2(get_el_br(b_wbv,lbw1,j-ubw1,j), get_el_br(b_wbv,lbw1,j-ubw1+1,j))
          jsu(j-k0+1,k+1)=j-ubw1
          csu(j-k0+1,k+1)=rot%cosine; ssu(j-k0+1,k+1)=rot%sine
          call rotation_times_tbr(trp_rot(rot),b_wbv,n,lbw1,ubw1,k+1,0,j-ubw1)
       end do
    end do
    if (.not. full_lbw) then
       call shift(b_wbv,0,-1)
    end if
    call br_to_bc(b_wbv,b_ubt,lbw,ubw)
  end subroutine f_z_convert_wbv_to_ubt


end module mod_convert_wbv_to_ubt
