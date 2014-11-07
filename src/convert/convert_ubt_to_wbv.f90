module mod_convert_ubt_to_wbv
  use mod_prec
  use mod_error_id
  use mod_shift
  use mod_rotation
  use mod_nested_types
  use mod_band_types
  implicit none

  private

  public :: convert_ubt_to_wbv, d_convert_ubt_to_wbv, c_convert_ubt_to_wbv, &
       f_convert_ubt_to_wbv, f_d_convert_ubt_to_wbv, f_c_convert_ubt_to_wbv

  public :: info_d_convert_ubt_to_wbv, info_c_convert_ubt_to_wbv

  interface convert_ubt_to_wbv
     module procedure d_convert_ubt_to_wbv, c_convert_ubt_to_wbv
  end interface convert_ubt_to_wbv

  interface f_convert_ubt_to_wbv
     module procedure f_d_convert_ubt_to_wbv, f_c_convert_ubt_to_wbv
  end interface f_convert_ubt_to_wbv

  type(routine_info), parameter :: info_d_convert_ubt_to_wbv=routine_info(id_d_convert_ubt_to_wbv, &
       'd_convert_ubt_to_wbv', &
       [ character(len=error_message_length) :: 'n<1', 'Insufficient storage in ubt', &
       'Insufficient storage in wbv.', 'ubt%n /= wbv%n' ] )

  type(routine_info), parameter :: info_c_convert_ubt_to_wbv=routine_info(id_c_convert_ubt_to_wbv, &
       'c_convert_ubt_to_wbv', &
       [ character(len=error_message_length) :: 'Insufficient storage in ubt', &
       'Insufficient storage in wbv.', 'ubt%n /= wbv%n' ] )

contains

  ! Errors:
  ! 0: no error
  ! 1: n<1
  ! 2: Insufficient storage in ubt
  ! 3: Insufficient stroage in wbv
  ! 4: ubt%n /= wbv%n

  subroutine d_convert_ubt_to_wbv(ubt, wbv, error)
    type(d_ubt) :: ubt
    type(d_wbv) :: wbv
    type(error_info), intent(out) :: error
    call clear_error(error)
    if (get_n(ubt) < 1) then
       call set_error(error, 1, id_d_convert_ubt_to_wbv); return
    end if
    if ((get_ubwmax(ubt) < ubt%ubw+1 .and. ubt%ubw < get_n(ubt)-1) .or. &
         (get_lbwmax(ubt) < ubt%lbw+1 .and. ubt%lbw < get_n(ubt)-1)) then
       call set_error(error, 2, id_d_convert_ubt_to_wbv); return
    end if
    if (get_lbwmax(wbv) < ubt%lbw .or. get_ubwmax(wbv) < ubt%ubw) then
       call set_error(error, 3, id_d_convert_ubt_to_wbv); return
    end if
    if (get_n(ubt) /= get_n(wbv)) then
       call set_error(error, 4, id_d_convert_ubt_to_wbv); return
    end if
    call f_d_convert_ubt_to_wbv(ubt%bc, get_n(ubt), ubt%lbw, ubt%ubw, get_lbwmax(ubt), &
         get_ubwmax(ubt), ubt%numrotsu, ubt%jsu, ubt%csu, ubt%ssu, &
         ubt%numrotst, ubt%kst, ubt%cst, ubt%sst, &
         wbv%br, wbv%lbw, wbv%ubw, get_lbwmax(wbv), get_ubwmax(wbv), &
         wbv%numrotsw, wbv%jsw, wbv%csw, wbv%ssw, &
         wbv%numrotsv, wbv%ksv, wbv%csv, wbv%ssv, error)
  end subroutine d_convert_ubt_to_wbv

  subroutine f_d_convert_ubt_to_wbv(b_ubt, n, lbw, ubw, lbwmax_ubt, ubwmax_ubt, numrotsu, &
       jsu, csu, ssu, numrotst, kst, cst, sst, &
       b_wbv, lbw_wbv, ubw_wbv, lbwmax_wbv, ubwmax_wbv, &
       numrotsw, jsw, csw, ssw, numrotsv, ksv, csv, ssv, error)
    real(kind=dp), dimension(lbwmax_ubt+ubwmax_ubt+1,n), intent(inout) :: b_ubt
    integer(kind=int32), intent(in) :: n, lbw, ubw, lbwmax_ubt, ubwmax_ubt, lbwmax_wbv, ubwmax_wbv
    integer(kind=int32), dimension(n), intent(in) :: numrotsu
    integer(kind=int32), dimension(ubwmax_ubt,n), intent(in) :: jsu
    real(kind=dp), dimension(ubwmax_ubt,n), intent(in) :: csu, ssu
    integer(kind=int32), dimension(n), intent(in) :: numrotst
    integer(kind=int32), dimension(n,lbwmax_ubt), intent(in) :: kst
    real(kind=dp), dimension(n,lbwmax_ubt), intent(in) :: cst, sst

    real(kind=dp), dimension(n,lbwmax_wbv+ubwmax_wbv+1), intent(out) :: b_wbv
    integer(kind=int32), dimension(n), intent(out) :: numrotsw
    integer(kind=int32), dimension(lbwmax_wbv,n), intent(out) :: jsw
    real(kind=dp), dimension(lbwmax_wbv,n), intent(out) :: csw, ssw
    integer(kind=int32), dimension(n), intent(out) :: numrotsv
    integer(kind=int32), dimension(n,ubwmax_wbv), intent(out) :: ksv
    real(kind=dp), dimension(n,ubwmax_wbv), intent(out) :: csv, ssv
    integer(kind=int32), intent(out) :: lbw_wbv, ubw_wbv
    type(error_info), intent(out) :: error

    integer(kind=int32) :: j, k, k0, k1, ubw1, lbw1, j0, j1
    type(d_rotation) :: rot
    logical :: full_ubw

    call clear_error(error)
    b_wbv(:,1:lbw+ubw+1)=0.0_dp; numrotsv=0
    ssv(:,1:ubw)=0.0_dp; csv(:,1:ubw)=0.0_dp
    ksv(:,1:ubw)=0
    numrotsw=0
    ssw(1:lbw,:)=0.0_dp; csw(1:lbw,:)=0.0_dp
    jsw(1:lbw,:)=0
    lbw_wbv=lbw; ubw_wbv=ubw
    if (n == 1) then
       b_wbv(1,1)=b_ubt(1,1);
       lbw_wbv=0; ubw_wbv=0; return
    end if
    ! must allow for temporary fill-in
    if (ubw < n-1) then
       ubw1=ubw+1
       call shift2(b_ubt,1,0)
       full_ubw=.false.
    else
       ubw1=ubw
       full_ubw=.true.
    end if
    if (lbw < n-1) then
       lbw1=lbw+1
       b_ubt(lbw1+ubw1+1,:)=0.0_dp
    else
       lbw1=lbw
    end if
    ! Do the lower triangular part first
    ! k is the size of the leading principal submatrix
    do k=n-1,2,-1
       ! Apply T_k
       do j=1,numrotst(k)
          rot%cosine=cst(k,j); rot%sine=sst(k,j)
          call tbc_times_rotation(b_ubt,n,lbw1,ubw1,k,0,trp_rot(rot),kst(k,j))
       end do
       ! Rows in which nonzeros have been introduced into the extra subdiagonal.
       j0=max(k+1,lbw+2)
       j1=min(k+lbw,n)
       ! Apply W_{k-1}
       numrotsw(k-1)=max(j1-j0+1,0)
       do j=j0,j1
          rot=lgivens(get_el_bc(b_ubt,ubw1,j-1,j-lbw1),get_el_bc(b_ubt,ubw1,j,j-lbw1))
          jsw(j1-j+1,k-1)=j-1
          csw(j1-j+1,k-1)=rot%cosine; ssw(j1-j+1,k-1)=rot%sine
          call rotation_times_tbc(trp_rot(rot),b_ubt,n,lbw1,ubw1,0,n-k+1,j-1)
       end do
    end do
    ! Upper triangular part
    do k=1,n-2 ! size of trailing principal submatrix
       do j=1,numrotsu(n-k)
          rot%cosine=csu(j,n-k); rot%sine=ssu(j,n-k)
          call rotation_times_tbc(rot,b_ubt,n,lbw1,ubw1,n-k,0,jsu(j,n-k))
       end do
       k0=max(n-k+1,ubw+2)
       k1=min(n-k+ubw,n)
       numrotsv(n-(k+1))=max(k1-k0+1,0)
       do j=k0,k1
          rot=rgivens(get_el_bc(b_ubt,ubw1,j-ubw1,j-1), get_el_bc(b_ubt,ubw1,j-ubw1,j))
          ksv(n-(k+1),k1-j+1)=j-1
          csv(n-(k+1),k1-j+1)=rot%cosine; ssv(n-(k+1),k1-j+1)=rot%sine
          call tbc_times_rotation(b_ubt,n,lbw1,ubw1,0,k+1,rot,j-1)
       end do
    end do
    ! Store the results in b_wbv
    if (.not. full_ubw) then
       call shift2(b_ubt,-1,0)
    end if
    call bc_to_br(b_ubt, b_wbv, lbw, ubw)
  end subroutine f_d_convert_ubt_to_wbv


  ! Errors:
  ! 0: no error
  ! 1: n<1
  ! 2: Insufficient storage in ubt
  ! 3: Insufficient stroage in wbv
  ! 4: ubt%n /= wbv%n

  subroutine c_convert_ubt_to_wbv(ubt, wbv, error)
    type(c_ubt) :: ubt
    type(c_wbv) :: wbv
    type(error_info), intent(out) :: error
    call clear_error(error)
    if (get_n(ubt) < 1) then
       call set_error(error, 1, id_c_convert_ubt_to_wbv); return
    end if
    if ((get_ubwmax(ubt) < ubt%ubw+1 .and. ubt%ubw < get_n(ubt)-1) .or. &
         (get_lbwmax(ubt) < ubt%lbw+1 .and. ubt%lbw < get_n(ubt)-1)) then
       call set_error(error, 2, id_c_convert_ubt_to_wbv); return
    end if
    if (get_lbwmax(wbv) < ubt%lbw .or. get_ubwmax(wbv) < ubt%ubw) then
       call set_error(error, 3, id_c_convert_ubt_to_wbv); return
    end if
    if (get_n(ubt) /= get_n(wbv)) then
       call set_error(error, 4, id_c_convert_ubt_to_wbv); return
    end if
    call f_c_convert_ubt_to_wbv(ubt%bc, get_n(ubt), ubt%lbw, ubt%ubw, get_lbwmax(ubt), &
         get_ubwmax(ubt), ubt%numrotsu, ubt%jsu, ubt%csu, ubt%ssu, &
         ubt%numrotst, ubt%kst, ubt%cst, ubt%sst, &
         wbv%br, wbv%lbw, wbv%ubw, get_lbwmax(wbv), get_ubwmax(wbv), &
         wbv%numrotsw, wbv%jsw, wbv%csw, wbv%ssw, &
         wbv%numrotsv, wbv%ksv, wbv%csv, wbv%ssv, error)
  end subroutine c_convert_ubt_to_wbv

  subroutine f_c_convert_ubt_to_wbv(b_ubt, n, lbw, ubw, lbwmax_ubt, ubwmax_ubt, numrotsu, &
       jsu, csu, ssu, numrotst, kst, cst, sst, &
       b_wbv, lbw_wbv, ubw_wbv, lbwmax_wbv, ubwmax_wbv, &
       numrotsw, jsw, csw, ssw, numrotsv, ksv, csv, ssv, error)
    complex(kind=dp), dimension(lbwmax_ubt+ubwmax_ubt+1,n), intent(inout) :: b_ubt
    integer(kind=int32), intent(in) :: n, lbw, ubw, lbwmax_ubt, ubwmax_ubt, lbwmax_wbv, ubwmax_wbv
    integer(kind=int32), dimension(n), intent(in) :: numrotsu
    integer(kind=int32), dimension(ubwmax_ubt,n), intent(in) :: jsu
    real(kind=dp), dimension(ubwmax_ubt,n), intent(in) :: csu
    complex(kind=dp), dimension(ubwmax_ubt,n), intent(in) :: ssu
    integer(kind=int32), dimension(n), intent(in) :: numrotst
    integer(kind=int32), dimension(n,lbwmax_ubt), intent(in) :: kst
    real(kind=dp), dimension(n,lbwmax_ubt), intent(in) :: cst
    complex(kind=dp), dimension(n,lbwmax_ubt), intent(in) :: sst

    complex(kind=dp), dimension(n,lbwmax_wbv+ubwmax_wbv+1), intent(out) :: b_wbv
    integer(kind=int32), dimension(n), intent(out) :: numrotsw
    integer(kind=int32), dimension(lbwmax_wbv,n), intent(out) :: jsw
    real(kind=dp), dimension(lbwmax_wbv,n), intent(out) :: csw
    complex(kind=dp), dimension(lbwmax_wbv,n), intent(out) :: ssw
    integer(kind=int32), dimension(n), intent(out) :: numrotsv
    integer(kind=int32), dimension(n,ubwmax_wbv), intent(out) :: ksv
    real(kind=dp), dimension(n,ubwmax_wbv), intent(out) :: csv
    complex(kind=dp), dimension(n,ubwmax_wbv), intent(out) :: ssv
    integer(kind=int32), intent(out) :: lbw_wbv, ubw_wbv
    type(error_info), intent(out) :: error

    integer(kind=int32) :: j, k, k0, k1, ubw1, lbw1, j0, j1
    type(c_rotation) :: rot
    logical :: full_ubw

    call clear_error(error)
    b_wbv(:,1:lbw+ubw+1)=(0.0_dp, 0.0_dp); numrotsv=0
    ssv(:,1:ubw)=(0.0_dp, 0.0_dp); csv(:,1:ubw)=0.0_dp
    ksv(:,1:ubw)=0
    numrotsw=0
    ssw(1:lbw,:)=(0.0_dp, 0.0_dp); csw(1:lbw,:)=0.0_dp
    jsw(1:lbw,:)=0
    lbw_wbv=lbw; ubw_wbv=ubw
    if (n == 1) then
       b_wbv(1,1)=b_ubt(1,1);
       lbw_wbv=0; ubw_wbv=0; return
    end if
    ! must allow for temporary fill-in
    if (ubw < n-1) then
       ubw1=ubw+1
       call shift2(b_ubt,1,0)
       full_ubw=.false.
    else
       ubw1=ubw
       full_ubw=.true.
    end if
    if (lbw < n-1) then
       lbw1=lbw+1
       b_ubt(lbw1+ubw1+1,:)=(0.0_dp, 0.0_dp)
    else
       lbw1=lbw
    end if
    ! Do the lower triangular part first
    ! k is the size of the leading principal submatrix
    do k=n-1,2,-1
       ! Apply T_k
       do j=1,numrotst(k)
          rot%cosine=cst(k,j); rot%sine=sst(k,j)
          call tbc_times_rotation(b_ubt,n,lbw1,ubw1,k,0,trp_rot(rot),kst(k,j))
       end do
       ! Rows in which nonzeros have been introduced into the extra subdiagonal.
       j0=max(k+1,lbw+2)
       j1=min(k+lbw,n)
       ! Apply W_{k-1}
       numrotsw(k-1)=max(j1-j0+1,0)
       do j=j0,j1
          rot=lgivens(get_el_bc(b_ubt,ubw1,j-1,j-lbw1),get_el_bc(b_ubt,ubw1,j,j-lbw1))
          jsw(j1-j+1,k-1)=j-1
          csw(j1-j+1,k-1)=rot%cosine; ssw(j1-j+1,k-1)=rot%sine
          call rotation_times_tbc(trp_rot(rot),b_ubt,n,lbw1,ubw1,0,n-k+1,j-1)
       end do
    end do
    ! Upper triangular part
    do k=1,n-2 ! size of trailing principal submatrix
       do j=1,numrotsu(n-k)
          rot%cosine=csu(j,n-k); rot%sine=ssu(j,n-k)
          call rotation_times_tbc(rot,b_ubt,n,lbw1,ubw1,n-k,0,jsu(j,n-k))
       end do
       k0=max(n-k+1,ubw+2)
       k1=min(n-k+ubw,n)
       numrotsv(n-(k+1))=max(k1-k0+1,0)
       do j=k0,k1
          rot=rgivens(get_el_bc(b_ubt,ubw1,j-ubw1,j-1), get_el_bc(b_ubt,ubw1,j-ubw1,j))
          ksv(n-(k+1),k1-j+1)=j-1
          csv(n-(k+1),k1-j+1)=rot%cosine; ssv(n-(k+1),k1-j+1)=rot%sine
          call tbc_times_rotation(b_ubt,n,lbw1,ubw1,0,k+1,rot,j-1)
       end do
    end do
    ! Store the results in b_wbv
    if (.not. full_ubw) then
       call shift2(b_ubt,-1,0)
    end if
    call bc_to_br(b_ubt, b_wbv, lbw, ubw)
  end subroutine f_c_convert_ubt_to_wbv

end module mod_convert_ubt_to_wbv
