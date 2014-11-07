module mod_convert_ub_to_bv
  use mod_prec
  use mod_error_id
  use mod_shift
  use mod_rotation
  use mod_nested_types
  use mod_band_types
  implicit none

  private

  public :: convert_ub_to_bv, d_convert_ub_to_bv, c_convert_ub_to_bv, &
       f_convert_ub_to_bv, f_d_convert_ub_to_bv, f_c_convert_ub_to_bv

  public :: info_d_convert_ub_to_bv, info_c_convert_ub_to_bv

  interface convert_ub_to_bv
     module procedure d_convert_ub_to_bv, c_convert_ub_to_bv
  end interface convert_ub_to_bv

  interface f_convert_ub_to_bv
     module procedure f_d_convert_ub_to_bv, f_c_convert_ub_to_bv
  end interface f_convert_ub_to_bv

  type(routine_info), parameter :: info_d_convert_ub_to_bv=routine_info(id_d_convert_ub_to_bv, &
       'd_convert_ub_to_bv', &
       [ character(len=error_message_length) :: 'n<1', 'Insufficient storage in ub', &
       'Insufficient storage in bv.', 'ub%n /= bv%n' ] )

  type(routine_info), parameter :: info_c_convert_ub_to_bv=routine_info(id_c_convert_ub_to_bv, &
       'c_convert_ub_to_bv', &
       [ character(len=error_message_length) :: 'Insufficient storage in ub', &
       'Insufficient storage in bv.', 'ub%n /= bv%n' ] )

contains

  ! Errors:
  ! 0: no error
  ! 1: n<1
  ! 2: Insufficient storage in ub
  ! 3: Insufficient stroage in bv
  ! 4: ub%n /= bv%n

  subroutine d_convert_ub_to_bv(ub, bv, error)
    type(d_ub) :: ub
    type(d_bv) :: bv
    type(error_info), intent(out) :: error
    call clear_error(error)
    if (get_n(ub) < 1) then
       call set_error(error, 1, id_d_convert_ub_to_bv); return
    end if
    if (get_ubwmax(ub) < ub%ubw+1 .and. ub%ubw < get_n(ub)-1) then
       call set_error(error, 2, id_d_convert_ub_to_bv); return
    end if
    if (get_lbwmax(bv) < ub%lbw .or. get_ubwmax(bv) < ub%ubw) then
       call set_error(error, 3, id_d_convert_ub_to_bv); return
    end if
    if (get_n(ub) /= get_n(bv)) then
       call set_error(error, 4, id_d_convert_ub_to_bv); return
    end if
    call f_d_convert_ub_to_bv(ub%bc, get_n(ub), ub%lbw, ub%ubw, get_lbwmax(ub), &
         get_ubwmax(ub), ub%numrotsu, ub%jsu, ub%csu, ub%ssu, bv%br, bv%lbw, bv%ubw, get_lbwmax(bv), &
         get_ubwmax(bv), bv%numrotsv, bv%ksv, bv%csv, bv%ssv, error)
  end subroutine d_convert_ub_to_bv

  subroutine f_d_convert_ub_to_bv(b_ub, n, lbw, ubw, lbwmax_ub, ubwmax_ub, numrotsu, &
       jsu, csu, ssu, b_bv, lbw_bv, ubw_bv, lbwmax_bv, ubwmax_bv, numrotsv, ksv, csv, ssv, error)
    real(kind=dp), dimension(lbwmax_ub+ubwmax_ub+1,n), intent(inout) :: b_ub
    integer(kind=int32), intent(in) :: n, lbw, ubw, lbwmax_ub, ubwmax_ub, lbwmax_bv, ubwmax_bv
    integer(kind=int32), dimension(n), intent(in) :: numrotsu
    integer(kind=int32), dimension(ubwmax_ub,n), intent(in) :: jsu
    real(kind=dp), dimension(ubwmax_ub,n), intent(in) :: csu, ssu

    real(kind=dp), dimension(n,lbwmax_bv+ubwmax_bv+1), intent(out) :: b_bv
    integer(kind=int32), dimension(n), intent(out) :: numrotsv
    integer(kind=int32), dimension(n,ubwmax_bv), intent(out) :: ksv
    real(kind=dp), dimension(n,ubwmax_bv), intent(out) :: csv, ssv
    integer(kind=int32), intent(out) :: lbw_bv, ubw_bv
    type(error_info), intent(out) :: error

    integer(kind=int32) :: j, k, k0, k1, ubw1, lbw1
    type(d_rotation) :: rot
    logical :: full_ubw

    call clear_error(error)
    b_bv(:,1:lbw+ubw+1)=0.0_dp; numrotsv=0
    ssv(:,1:ubw)=0.0_dp; csv(:,1:ubw)=0.0_dp
    ksv(:,1:ubw)=0
    lbw_bv=lbw; ubw_bv=ubw
    if (n == 1) then
       b_bv(1,1)=b_ub(1,1);
       lbw_bv=0; ubw_bv=0; numrotsv=0; return
    end if
    ! must allow for temporary fill-in
    lbw1=lbw
    if (ubw < n-1) then
       ubw1=ubw+1
       call shift2(b_ub,1,0)
       full_ubw=.false.
    else
       ubw1=ubw
       full_ubw=.true.
    end if
    do k=1,n-2 ! size of trailing principal submatrix
       do j=1,numrotsu(n-k)
          rot%cosine=csu(j,n-k); rot%sine=ssu(j,n-k)
          call rotation_times_tbc(rot,b_ub,n,lbw1,ubw1,n-k,0,jsu(j,n-k))
       end do
       k0=max(n-k+1,ubw+2)
       k1=min(n-k+ubw,n)
       numrotsv(n-(k+1))=max(k1-k0+1,0)
       do j=k0,k1
          rot=rgivens(get_el_bc(b_ub,ubw1,j-ubw1,j-1), get_el_bc(b_ub,ubw1,j-ubw1,j))
          ksv(n-(k+1),k1-j+1)=j-1
          csv(n-(k+1),k1-j+1)=rot%cosine; ssv(n-(k+1),k1-j+1)=rot%sine
          call tbc_times_rotation(b_ub,n,lbw1,ubw1,0,k+1,rot,j-1)
       end do
    end do
    ! Store the results in b_bv
    if (.not. full_ubw) then
       call shift2(b_ub,-1,0)
    end if
    call bc_to_br(b_ub, b_bv, lbw, ubw)
  end subroutine f_d_convert_ub_to_bv

  subroutine c_convert_ub_to_bv(ub, bv, error)
    type(c_ub) :: ub
    type(c_bv) :: bv
    type(error_info), intent(out) :: error
    call clear_error(error)
    if (get_n(ub) < 1) then
       call set_error(error, 1, id_c_convert_ub_to_bv); return
    end if
    if (get_ubwmax(ub) < ub%ubw+1 .and. ub%ubw < get_n(ub)-1) then
       call set_error(error, 2, id_c_convert_ub_to_bv); return
    end if
    if (get_lbwmax(bv) < ub%lbw .or. get_ubwmax(bv) < ub%ubw) then
       call set_error(error, 3, id_c_convert_ub_to_bv); return
    end if
    if (get_n(ub) /= get_n(bv)) then
       call set_error(error, 4, id_c_convert_ub_to_bv); return
    end if
    call f_c_convert_ub_to_bv(ub%bc, get_n(ub), ub%lbw, ub%ubw, get_lbwmax(ub), &
         get_ubwmax(ub), ub%numrotsu, ub%jsu, ub%csu, ub%ssu, bv%br, bv%lbw, bv%ubw, get_lbwmax(bv), &
         get_ubwmax(bv), bv%numrotsv, bv%ksv, bv%csv, bv%ssv, error)
  end subroutine c_convert_ub_to_bv

  subroutine f_c_convert_ub_to_bv(b_ub, n, lbw, ubw, lbwmax_ub, ubwmax_ub, numrotsu, &
       jsu, csu, ssu, b_bv, lbw_bv, ubw_bv, lbwmax_bv, ubwmax_bv, numrotsv, ksv, csv, &
       ssv, error)
    complex(kind=dp), dimension(lbwmax_ub+ubwmax_bv+1,n), intent(inout) :: b_ub
    integer(kind=int32), intent(in) :: n, lbw, ubw, lbwmax_ub, ubwmax_ub, lbwmax_bv, ubwmax_bv
    integer(kind=int32), dimension(n), intent(in) :: numrotsu
    integer(kind=int32), dimension(ubwmax_ub,n), intent(in) :: jsu
    real(kind=dp), dimension(ubwmax_ub,n), intent(in) :: csu
    complex(kind=dp), dimension(ubwmax_ub,n), intent(in) :: ssu

    complex(kind=dp), dimension(n,lbwmax_bv+ubwmax_bv+1), intent(out) :: b_bv
    integer(kind=int32), dimension(n), intent(out) :: numrotsv
    integer(kind=int32), dimension(n,ubwmax_bv), intent(out) :: ksv
    real(kind=dp), dimension(n,ubwmax_bv), intent(out) :: csv
    complex(kind=dp), dimension(n,ubwmax_bv), intent(out) :: ssv
    integer(kind=int32), intent(out) :: lbw_bv, ubw_bv
    type(error_info), intent(out) :: error

    integer(kind=int32) :: j, k, k0, k1, ubw1, lbw1
    type(c_rotation) :: rot
    logical :: full_ubw

    call clear_error(error)
    b_bv(:,1:lbw+ubw+1)=(0.0_dp,0.0_dp); numrotsv=0
    ssv(:,1:ubw)=(0.0_dp, 0.0_dp); csv(:,1:ubw)=0.0_dp
    ksv(:,1:ubw)=0
    lbw_bv=lbw; ubw_bv=ubw

    if (n == 1) then
       b_bv(1,1)=b_ub(1,1);
       lbw_bv=0; ubw_bv=0; numrotsv=0; return
    end if
    ! must allow for temporary fill-in
    ubw1=ubw+1; lbw1=lbw
    if (ubw < n-1) then
       ubw1=ubw+1
       call shift2(b_ub,1,0)
       full_ubw=.false.
    else
       ubw1=ubw
       full_ubw=.true.
    end if
    do k=1,n-2 ! trailing principal submatrix
       do j=1,numrotsu(n-k)
          rot%cosine=csu(j,n-k); rot%sine=ssu(j,n-k)
          call rotation_times_tbc(rot,b_ub,n,lbw1,ubw1,n-k,0,jsu(j,n-k))
       end do
       k0=max(n-k+1,ubw+2)
       k1=min(n-k+ubw,n)
       numrotsv(n-(k+1))=max(k1-k0+1,0)
       do j=k0,k1
          rot=rgivens(get_el_bc(b_ub,ubw1,j-ubw1,j-1), get_el_bc(b_ub,ubw1,j-ubw1,j))
          ksv(n-(k+1),k1-j+1)=j-1
          csv(n-(k+1),k1-j+1)=rot%cosine; ssv(n-(k+1),k1-j+1)=rot%sine
          call tbc_times_rotation(b_ub,n,lbw1,ubw1,0,k+1,rot,j-1)
       end do
    end do
    if (.not. full_ubw) then
       call shift2(b_ub,-1,0)
    end if
    call bc_to_br(b_ub, b_bv, lbw, ubw)
  end subroutine f_c_convert_ub_to_bv

end module mod_convert_ub_to_bv
