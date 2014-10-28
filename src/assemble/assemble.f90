module assemble
  use prec
  use error_id
  use band_types
  use nested_types
  use rotation
  implicit none

  private

  public :: ub_to_upper, d_ub_to_upper, c_ub_to_upper, &
       f_ub_to_upper, f_d_ub_to_upper, f_c_ub_to_upper, &
       bt_to_lower, d_bt_to_lower, c_bt_to_lower, &
       f_bt_to_lower, f_d_bt_to_lower, f_c_bt_to_lower, &
       ubt_to_general, d_ubt_to_general, c_ubt_to_general, &
       f_ubt_to_general, f_d_ubt_to_general, f_c_ubt_to_general, &
       bv_to_upper, d_bv_to_upper, c_bv_to_upper, &
       f_bv_to_upper, f_d_bv_to_upper, f_c_bv_to_upper, &
       wb_to_lower, d_wb_to_lower, c_wb_to_lower, &
       f_wb_to_lower, f_d_wb_to_lower, f_c_wb_to_lower, &
       wbv_to_general, d_wbv_to_general, c_wbv_to_general, &
       f_wbv_to_general, f_d_wbv_to_general, f_c_wbv_to_general

  public :: info_d_ub_to_upper, info_c_ub_to_upper, info_d_bt_to_lower, info_c_bt_to_lower, &
       info_d_ubt_to_general, info_c_ubt_to_general, info_d_bv_to_upper, &
       info_c_bv_to_upper, info_d_wb_to_lower, info_c_wb_to_lower, info_d_wbv_to_general, &
       info_c_wbv_to_general

  interface ub_to_upper
     module procedure d_ub_to_upper, c_ub_to_upper
  end interface ub_to_upper

  interface f_ub_to_upper
     module procedure f_d_ub_to_upper, f_c_ub_to_upper
  end interface f_ub_to_upper

  interface bt_to_lower
     module procedure d_bt_to_lower, c_bt_to_lower
  end interface bt_to_lower

  interface f_bt_to_lower
     module procedure f_d_bt_to_lower, f_c_bt_to_lower
  end interface f_bt_to_lower

  interface ubt_to_general
     module procedure d_ubt_to_general, c_ubt_to_general
  end interface ubt_to_general

  interface f_ubt_to_general
     module procedure f_d_ubt_to_general, f_c_ubt_to_general
  end interface f_ubt_to_general

  interface bv_to_upper
     module procedure d_bv_to_upper, c_bv_to_upper
  end interface bv_to_upper

  interface f_bv_to_upper
     module procedure f_d_bv_to_upper, f_c_bv_to_upper
  end interface f_bv_to_upper

  interface wb_to_lower
     module procedure d_wb_to_lower, c_wb_to_lower
  end interface wb_to_lower

  interface f_wb_to_lower
     module procedure f_d_wb_to_lower, f_c_wb_to_lower
  end interface f_wb_to_lower

  interface wbv_to_general
     module procedure d_wbv_to_general, c_wbv_to_general
  end interface wbv_to_general

  interface f_wbv_to_general
     module procedure f_d_wbv_to_general, f_c_wbv_to_general
  end interface f_wbv_to_general

  type(routine_info), parameter :: info_d_ub_to_upper=routine_info(id_d_ub_to_upper, 'd_ub_to_upper', &
       [ character(len=error_message_length) :: 'Size error in A.' ] )
  type(routine_info), parameter :: info_c_ub_to_upper=routine_info(id_c_ub_to_upper, 'c_ub_to_upper', &
       [ character(len=error_message_length) :: 'Size error in A.' ] )

  type(routine_info), parameter :: info_d_bt_to_lower=routine_info(id_d_bt_to_lower, 'd_bt_to_lower', &
       [ character(len=error_message_length) :: 'Size error in A.' ] )
  type(routine_info), parameter :: info_c_bt_to_lower=routine_info(id_c_bt_to_lower, 'c_bt_to_lower', &
       [ character(len=error_message_length) :: 'Size error in A.' ] )

  type(routine_info), parameter :: info_d_ubt_to_general=routine_info(id_d_ubt_to_general, 'd_ubt_to_general', &
       [ character(len=error_message_length) :: 'Size error in A.' ] )
  type(routine_info), parameter :: info_c_ubt_to_general=routine_info(id_c_ubt_to_general, 'c_ubt_to_general', &
       [ character(len=error_message_length) :: 'Size error in A.' ] )

  type(routine_info), parameter :: info_d_bv_to_upper=routine_info(id_d_bv_to_upper, 'd_bv_to_upper', &
       [ character(len=error_message_length) :: 'Size error in A.' ] )
  type(routine_info), parameter :: info_c_bv_to_upper=routine_info(id_c_bv_to_upper, 'c_bv_to_upper', &
       [ character(len=error_message_length) :: 'Size error in A.' ] )

  type(routine_info), parameter :: info_d_wb_to_lower=routine_info(id_d_wb_to_lower, 'd_wb_to_lower', &
       [ character(len=error_message_length) :: 'Size error in A.' ] )
  type(routine_info), parameter :: info_c_wb_to_lower=routine_info(id_c_wb_to_lower, 'c_wb_to_lower', &
       [ character(len=error_message_length) :: 'Size error in A.' ] )

  type(routine_info), parameter :: info_d_wbv_to_general=routine_info(id_d_wbv_to_general, 'd_wbv_to_general', &
       [ character(len=error_message_length) :: 'Size error in A.' ] )
  type(routine_info), parameter :: info_c_wbv_to_general=routine_info(id_c_wbv_to_general, 'c_wbv_to_general', &
       [ character(len=error_message_length) :: 'Size error in A.' ] )

contains

  ! Errors
  ! 0: no error
  ! 1: ub%n /= n
  subroutine d_ub_to_upper(ub,a,error)
    type(d_ub), intent(in) :: ub
    real(kind=dp), dimension(:,:), intent(out) :: a
    type(error_info), intent(out) :: error
    call clear_error(error)
    if (get_n(ub) /= size(a,1) .or. get_n(ub) /= size(a,2)) then
       call set_error(error, 1, id_d_ub_to_upper); return
    end if
    call f_d_ub_to_upper(ub%bc, get_n(ub), ub%lbw, ub%ubw, get_lbwmax(ub), get_ubwmax(ub), &
         ub%numrotsu, ub%jsu, ub%csu, ub%ssu, a)
  end subroutine d_ub_to_upper

  subroutine f_d_ub_to_upper(bc, n, lbw, ubw, lbwmax, ubwmax, numrotsu, jsu, csu, ssu, a)
    real(kind=dp), target, dimension(n,n), intent(out) :: a
    integer(kind=int32), dimension(ubwmax,n), intent(in) :: jsu
    real(kind=dp), dimension(ubwmax,n), intent(in) :: csu, ssu
    real(kind=dp), dimension(lbwmax+ubwmax+1,n), intent(in) :: bc
    integer(kind=int32), dimension(n), intent(in) :: numrotsu
    integer(kind=int32), intent(in) :: ubw, lbw, n, lbwmax, ubwmax
    !
    integer(kind=int32) :: j,k
    type(d_rotation) :: rot
    call bc_to_general(bc,lbw,ubw,a)
    if (n==1) then
       return
    end if
    do k=n-1,2,-1
       do j=1,numrotsu(k)
          rot%cosine=csu(j,k); rot%sine=ssu(j,k)
          call rotation_times_general(rot,a(:,k+1:n),jsu(j,k),jsu(j,k)+1)
       end do
    end do
  end subroutine f_d_ub_to_upper


  subroutine c_ub_to_upper(ub,a,error)
    type(c_ub), intent(in) :: ub
    complex(kind=dp), dimension(:,:), intent(out) :: a
    type(error_info), intent(out) :: error
    call clear_error(error)
    if (get_n(ub) /= size(a,1) .or. get_n(ub) /= size(a,2)) then
       call set_error(error, 1, id_c_ub_to_upper); return
    end if
    call f_c_ub_to_upper(ub%bc, get_n(ub), ub%lbw, ub%ubw, get_lbwmax(ub), get_ubwmax(ub), &
         ub%numrotsu, ub%jsu, ub%csu, ub%ssu, a)
  end subroutine c_ub_to_upper

  subroutine f_c_ub_to_upper(bc, n, lbw, ubw, lbwmax, ubwmax, numrotsu, jsu, csu, ssu, a)
    complex(kind=dp), target, dimension(n,n), intent(out) :: a
    integer(kind=int32), dimension(ubwmax,n), intent(in) :: jsu
    complex(kind=dp), dimension(ubwmax,n), intent(in) :: csu, ssu
    complex(kind=dp), dimension(lbwmax+ubwmax+1,n), intent(in) :: bc
    integer(kind=int32), dimension(n), intent(in) :: numrotsu
    integer(kind=int32), intent(in) :: ubw, lbw, n, lbwmax, ubwmax
    !
    integer(kind=int32) :: j,k
    type(c_rotation) :: rot

    call bc_to_general(bc,lbw,ubw,a)
    if (n==1) then
       return
    end if
    do k=n-1,2,-1
       do j=1,numrotsu(k)
          rot%cosine=csu(j,k); rot%sine=ssu(j,k)
          call rotation_times_general(rot,a(:,k+1:n),jsu(j,k),jsu(j,k)+1)
       end do
    end do
  end subroutine f_c_ub_to_upper

  ! Errors
  ! 0: no error
  ! 1: bt%n /= n
  subroutine d_bt_to_lower(bt,a,error)
    type(d_bt), intent(in) :: bt
    real(kind=dp), dimension(:,:), intent(out) :: a
    type(error_info), intent(out) :: error
    call clear_error(error)
    if (get_n(bt) /= size(a,1) .or. get_n(bt) /= size(a,2)) then
       call set_error(error, 1, id_d_bt_to_lower); return
    end if
    call f_d_bt_to_lower(bt%br, get_n(bt), bt%lbw, bt%ubw, get_lbwmax(bt), get_ubwmax(bt), &
         bt%numrotst, bt%kst, bt%cst, bt%sst, a)
  end subroutine d_bt_to_lower

  subroutine f_d_bt_to_lower(br, n, lbw, ubw, lbwmax, ubwmax, numrotst, kst, cst, sst, a)
    real(kind=dp), target, dimension(n,n), intent(out) :: a
    integer(kind=int32), dimension(n,lbwmax), intent(in) :: kst
    real(kind=dp), dimension(n,lbwmax), intent(in) :: cst, sst
    real(kind=dp), dimension(n,lbwmax+ubwmax+1), intent(in) :: br
    integer(kind=int32), dimension(n), intent(in) :: numrotst
    integer(kind=int32), intent(in) :: ubw, lbw, n, lbwmax, ubwmax
    !
    integer(kind=int32) :: j,k
    type(d_rotation) :: rot
    call br_to_general(br,lbw,ubw,a)
    if (n==1) then
       return
    end if
    do k=n-1,2,-1
       do j=1,numrotst(k)
          rot%cosine=cst(k,j); rot%sine=sst(k,j)
          call general_times_rotation(a(k+1:n,:), trp_rot(rot), kst(k,j),kst(k,j)+1)
       end do
    end do
  end subroutine f_d_bt_to_lower

  ! Errors
  ! 0: no error
  ! 1: bt%n /= n
  subroutine c_bt_to_lower(bt,a,error)
    type(c_bt), intent(in) :: bt
    complex(kind=dp), dimension(:,:), intent(out) :: a
    type(error_info), intent(out) :: error
    call clear_error(error)
    if (get_n(bt) /= size(a,1) .or. get_n(bt) /= size(a,2)) then
       call set_error(error, 1, id_c_bt_to_lower); return
    end if
    call f_c_bt_to_lower(bt%br, get_n(bt), bt%lbw, bt%ubw, get_lbwmax(bt), get_ubwmax(bt), &
         bt%numrotst, bt%kst, bt%cst, bt%sst, a)
  end subroutine c_bt_to_lower

  subroutine f_c_bt_to_lower(br, n, lbw, ubw, lbwmax, ubwmax, numrotst, kst, cst, sst, a)
    complex(kind=dp), target, dimension(n,n), intent(out) :: a
    integer(kind=int32), dimension(n,lbwmax), intent(in) :: kst
    complex(kind=dp), dimension(n,lbwmax), intent(in) :: cst, sst
    complex(kind=dp), dimension(n,lbwmax+ubwmax+1), intent(in) :: br
    integer(kind=int32), dimension(n), intent(in) :: numrotst
    integer(kind=int32), intent(in) :: ubw, lbw, n, lbwmax, ubwmax
    !
    integer(kind=int32) :: j,k
    type(c_rotation) :: rot
    call br_to_general(br,lbw,ubw,a)
    if (n==1) then
       return
    end if
    do k=n-1,2,-1
       do j=1,numrotst(k)
          rot%cosine=cst(k,j); rot%sine=sst(k,j)
          call general_times_rotation(a(k+1:n,:), trp_rot(rot), kst(k,j),kst(k,j)+1)
       end do
    end do
  end subroutine f_c_bt_to_lower

  ! Errors
  ! 0: no error
  ! 1: ubt%n /= n
  subroutine d_ubt_to_general(ubt,a,error)
    type(d_ubt), intent(in) :: ubt
    real(kind=dp), dimension(:,:), intent(out) :: a
    type(error_info), intent(out) :: error
    call clear_error(error)
    if (get_n(ubt) /= size(a,1) .or. get_n(ubt) /= size(a,2)) then
       call set_error(error, 1, id_d_ubt_to_general); return
    end if
    call f_d_ubt_to_general(ubt%bc, get_n(ubt), ubt%lbw, ubt%ubw, get_lbwmax(ubt), get_ubwmax(ubt), &
         ubt%numrotsu, ubt%jsu, ubt%csu, ubt%ssu, &
         ubt%numrotst, ubt%kst, ubt%cst, ubt%sst, a)
  end subroutine d_ubt_to_general

  subroutine f_d_ubt_to_general(bc, n, lbw, ubw, lbwmax, ubwmax, numrotsu, jsu, csu, ssu, &
       numrotst, kst, cst, sst, a)
    real(kind=dp), target, dimension(n,n), intent(out) :: a
    integer(kind=int32), dimension(ubwmax,n), intent(in) :: jsu
    real(kind=dp), dimension(ubwmax,n), intent(in) :: csu, ssu
    integer(kind=int32), dimension(n,lbwmax), intent(in) :: kst
    real(kind=dp), dimension(n,lbwmax), intent(in) :: cst, sst
    real(kind=dp), dimension(lbwmax+ubwmax+1,n), intent(in) :: bc
    integer(kind=int32), dimension(n), intent(in) :: numrotsu, numrotst
    integer(kind=int32), intent(in) :: ubw, lbw, n, lbwmax, ubwmax
    !
    integer(kind=int32) :: j,k
    type(d_rotation) :: rot
    call bc_to_general(bc,lbw,ubw,a)
    if (n==1) then
       return
    end if
    do k=n-1,2,-1
       do j=1,numrotsu(k)
          rot%cosine=csu(j,k); rot%sine=ssu(j,k)
          call rotation_times_general(rot,a(:,k+1:n),jsu(j,k),jsu(j,k)+1)
       end do
    end do

    do k=n-1,2,-1
       do j=1,numrotst(k)
          rot%cosine=cst(k,j); rot%sine=sst(k,j)
          call general_times_rotation(a(k+1:n,:), trp_rot(rot), kst(k,j),kst(k,j)+1)
       end do
    end do
  end subroutine f_d_ubt_to_general

  ! Errors
  ! 0: no error
  ! 1: ubt%n /= n
  subroutine c_ubt_to_general(ubt,a,error)
    type(c_ubt), intent(in) :: ubt
    complex(kind=dp), dimension(:,:), intent(out) :: a
    type(error_info), intent(out) :: error
    call clear_error(error)
    if (get_n(ubt) /= size(a,1) .or. get_n(ubt) /= size(a,2)) then
       call set_error(error, 1, id_c_ubt_to_general); return
    end if
    call f_c_ubt_to_general(ubt%bc, get_n(ubt), ubt%lbw, ubt%ubw, get_lbwmax(ubt), get_ubwmax(ubt), &
         ubt%numrotsu, ubt%jsu, ubt%csu, ubt%ssu, &
         ubt%numrotst, ubt%kst, ubt%cst, ubt%sst, a)
  end subroutine c_ubt_to_general

  subroutine f_c_ubt_to_general(bc, n, lbw, ubw, lbwmax, ubwmax, numrotsu, jsu, csu, ssu, &
       numrotst, kst, cst, sst, a)
    complex(kind=dp), target, dimension(n,n), intent(out) :: a
    integer(kind=int32), dimension(ubwmax,n), intent(in) :: jsu
    complex(kind=dp), dimension(ubwmax,n), intent(in) :: csu, ssu
    integer(kind=int32), dimension(n,lbwmax), intent(in) :: kst
    complex(kind=dp), dimension(n,lbwmax), intent(in) :: cst, sst
    complex(kind=dp), dimension(lbwmax+ubwmax+1,n), intent(in) :: bc
    integer(kind=int32), dimension(n), intent(in) :: numrotsu, numrotst
    integer(kind=int32), intent(in) :: ubw, lbw, n, lbwmax, ubwmax
    !
    integer(kind=int32) :: j,k
    type(c_rotation) :: rot
    call bc_to_general(bc,lbw,ubw,a)
    if (n==1) then
       return
    end if
    do k=n-1,2,-1
       do j=1,numrotsu(k)
          rot%cosine=csu(j,k); rot%sine=ssu(j,k)
          call rotation_times_general(rot,a(:,k+1:n),jsu(j,k),jsu(j,k)+1)
       end do
    end do

    do k=n-1,2,-1
       do j=1,numrotst(k)
          rot%cosine=cst(k,j); rot%sine=sst(k,j)
          call general_times_rotation(a(k+1:n,:), trp_rot(rot), kst(k,j),kst(k,j)+1)
       end do
    end do
  end subroutine f_c_ubt_to_general

  subroutine d_bv_to_upper(bv,a,error)
    type(d_bv), intent(in) :: bv
    real(kind=dp), dimension(:,:), intent(out) :: a
    type(error_info), intent(out) :: error
    call clear_error(error)
    if (get_n(bv) /= size(a,1) .or. get_n(bv) /= size(a,2)) then
       call set_error(error, 1, id_d_bv_to_upper); return
    end if
    call f_d_bv_to_upper(bv%br, get_n(bv), bv%lbw, bv%ubw, get_lbwmax(bv), get_ubwmax(bv), &
         bv%numrotsv, bv%ksv, bv%csv, bv%ssv, a)
  end subroutine d_bv_to_upper

  subroutine f_d_bv_to_upper(br, n, lbw, ubw, lbwmax, ubwmax, numrotsv, ksv, csv, ssv, a)
    real(kind=dp), target, dimension(n,n), intent(out) :: a
    integer(kind=int32), dimension(n,ubwmax), intent(in) :: ksv
    real(kind=dp), dimension(n, ubwmax), intent(in) :: csv, ssv
    real(kind=dp), dimension(n,lbwmax+ubwmax+1), intent(in) :: br
    integer(kind=int32), dimension(n), intent(in) :: numrotsv
    integer(kind=int32), intent(in) :: ubw, lbw, n, lbwmax, ubwmax
    !
    integer(kind=int32) :: j,k
    type(d_rotation) :: rot
    call br_to_general(br,lbw,ubw,a)
    if (n==1) then
       return
    end if
    do j=1,n-2
       do k=1,numrotsv(j)
          rot%cosine=csv(j,k); rot%sine=ssv(j,k)
          call general_times_rotation(a(1:j,:),trp_rot(rot),ksv(j,k), ksv(j,k)+1)
       end do
    end do
  end subroutine f_d_bv_to_upper

  subroutine c_bv_to_upper(bv,a,error)
    type(c_bv), intent(in) :: bv
    complex(kind=dp), dimension(:,:), intent(out) :: a
    type(error_info), intent(out) :: error
    call clear_error(error)
    if (get_n(bv) /= size(a,1) .or. get_n(bv) /= size(a,2)) then
       call set_error(error, 1, id_c_bv_to_upper); return
    end if
    call f_c_bv_to_upper(bv%br, get_n(bv), bv%lbw, bv%ubw, get_lbwmax(bv), get_ubwmax(bv), &
         bv%numrotsv, bv%ksv, bv%csv, bv%ssv, a)
  end subroutine c_bv_to_upper

  subroutine f_c_bv_to_upper(br, n, lbw, ubw, lbwmax, ubwmax, numrotsv, ksv, csv, ssv, a)
    complex(kind=dp), target, dimension(n,n), intent(out) :: a
    integer(kind=int32), dimension(n,ubwmax), intent(in) :: ksv
    complex(kind=dp), dimension(n, ubwmax), intent(in) :: csv, ssv
    complex(kind=dp), dimension(n,lbwmax+ubwmax+1), intent(in) :: br
    integer(kind=int32), dimension(n), intent(in) :: numrotsv
    integer(kind=int32), intent(in) :: ubw, lbw, n, lbwmax, ubwmax
    !
    integer(kind=int32) :: j,k
    type(c_rotation) :: rot
    call br_to_general(br,lbw,ubw,a)
    if (n==1) then
       return
    end if
    do j=1,n-2
       do k=1,numrotsv(j)
          rot%cosine=csv(j,k); rot%sine=ssv(j,k)
          call general_times_rotation(a(1:j,:),trp_rot(rot),ksv(j,k), ksv(j,k)+1)
       end do
    end do
  end subroutine f_c_bv_to_upper

  subroutine d_wb_to_lower(wb,a,error)
    type(d_wb), intent(in) :: wb
    real(kind=dp), dimension(:,:), intent(out) :: a
    type(error_info), intent(out) :: error
    call clear_error(error)
    if (get_n(wb) /= size(a,1) .or. get_n(wb) /= size(a,2)) then
       call set_error(error, 1, id_d_wb_to_lower); return
    end if
    call f_d_wb_to_lower(wb%bc, get_n(wb), wb%lbw, wb%ubw, get_lbwmax(wb), get_ubwmax(wb), &
         wb%numrotsw, wb%jsw, wb%csw, wb%ssw, a)
  end subroutine d_wb_to_lower

  subroutine f_d_wb_to_lower(bc, n, lbw, ubw, lbwmax, ubwmax, numrotsw, jsw, csw, ssw, a)
    real(kind=dp), target, dimension(n,n), intent(out) :: a
    integer(kind=int32), dimension(lbwmax,n), intent(in) :: jsw
    real(kind=dp), dimension(lbwmax,n), intent(in) :: csw, ssw
    real(kind=dp), dimension(lbwmax+ubwmax+1,n), intent(in) :: bc
    integer(kind=int32), dimension(n), intent(in) :: numrotsw
    integer(kind=int32), intent(in) :: ubw, lbw, n, lbwmax, ubwmax
    !
    integer(kind=int32) :: j,k
    type(d_rotation) :: rot
    call bc_to_general(bc,lbw,ubw,a)
    if (n==1) then
       return
    end if

    do k=1,n-2
       do j=1,numrotsw(k)
          rot%cosine=csw(j,k); rot%sine=ssw(j,k)
          call rotation_times_general(rot,a(:,1:k),jsw(j,k),jsw(j,k)+1)
       end do
    end do
  end subroutine f_d_wb_to_lower

  subroutine c_wb_to_lower(wb,a,error)
    type(c_wb), intent(in) :: wb
    complex(kind=dp), dimension(:,:), intent(out) :: a
    type(error_info), intent(out) :: error
    call clear_error(error)
    if (get_n(wb) /= size(a,1) .or. get_n(wb) /= size(a,2)) then
       call set_error(error, 1, id_c_wb_to_lower); return
    end if
    call f_c_wb_to_lower(wb%bc, get_n(wb), wb%lbw, wb%ubw, get_lbwmax(wb), get_ubwmax(wb), &
         wb%numrotsw, wb%jsw, wb%csw, wb%ssw, a)
  end subroutine c_wb_to_lower

  subroutine f_c_wb_to_lower(bc, n, lbw, ubw, lbwmax, ubwmax, numrotsw, jsw, csw, ssw, a)
    complex(kind=dp), target, dimension(n,n), intent(out) :: a
    integer(kind=int32), dimension(lbwmax,n), intent(in) :: jsw
    complex(kind=dp), dimension(lbwmax,n), intent(in) :: csw, ssw
    complex(kind=dp), dimension(lbwmax+ubwmax+1,n), intent(in) :: bc
    integer(kind=int32), dimension(n), intent(in) :: numrotsw
    integer(kind=int32), intent(in) :: ubw, lbw, n, lbwmax, ubwmax
    !
    integer(kind=int32) :: j,k
    type(c_rotation) :: rot
    call bc_to_general(bc,lbw,ubw,a)
    if (n==1) then
       return
    end if

    do k=1,n-2
       do j=1,numrotsw(k)
          rot%cosine=csw(j,k); rot%sine=ssw(j,k)
          call rotation_times_general(rot,a(:,1:k),jsw(j,k),jsw(j,k)+1)
       end do
    end do
  end subroutine f_c_wb_to_lower

  subroutine d_wbv_to_general(wbv,a,error)
    type(d_wbv), intent(in) :: wbv
    real(kind=dp), dimension(:,:), intent(out) :: a
    type(error_info), intent(out) :: error
    call clear_error(error)
    if (get_n(wbv) /= size(a,1) .or. get_n(wbv) /= size(a,2)) then
       call set_error(error, 1, id_d_wbv_to_general); return
    end if
    call f_d_wbv_to_general(wbv%br, get_n(wbv), wbv%lbw, wbv%ubw, get_lbwmax(wbv), get_ubwmax(wbv), &
         wbv%numrotsw, wbv%jsw, wbv%csw, wbv%ssw, &
         wbv%numrotsv, wbv%ksv, wbv%csv, wbv%ssv, a)
  end subroutine d_wbv_to_general

  subroutine f_d_wbv_to_general(br, n, lbw, ubw, lbwmax, ubwmax, numrotsw, jsw, csw, ssw, &
       numrotsv, ksv, csv, ssv, a)
    real(kind=dp), target, dimension(n,n), intent(out) :: a
    integer(kind=int32), dimension(lbwmax,n), intent(in) :: jsw
    real(kind=dp), dimension(lbwmax,n), intent(in) :: csw, ssw
    integer(kind=int32), dimension(n,ubwmax), intent(in) :: ksv
    real(kind=dp), dimension(n,ubwmax), intent(in) :: csv, ssv
    real(kind=dp), dimension(n,ubwmax+lbwmax+1), intent(in) :: br
    integer(kind=int32), dimension(n), intent(in) :: numrotsw, numrotsv
    integer(kind=int32), intent(in) :: ubw, lbw, n, lbwmax, ubwmax
    !
    integer(kind=int32) :: j,k
    type(d_rotation) :: rot
    call br_to_general(br,lbw,ubw,a)
    if (n==1) then
       return
    end if

    do j=1,n-2
       do k=1,numrotsv(j)
          rot%cosine=csv(j,k); rot%sine=ssv(j,k)
          call general_times_rotation(a(1:j,:),trp_rot(rot),ksv(j,k), ksv(j,k)+1)
       end do
    end do

    do k=1,n-2
       do j=1,numrotsw(k)
          rot%cosine=csw(j,k); rot%sine=ssw(j,k)
          call rotation_times_general(rot,a(:,1:k),jsw(j,k),jsw(j,k)+1)
       end do
    end do
  end subroutine f_d_wbv_to_general

  subroutine c_wbv_to_general(wbv,a,error)
    type(c_wbv), intent(in) :: wbv
    complex(kind=dp), dimension(:,:), intent(out) :: a
    type(error_info), intent(out) :: error
    call clear_error(error)
    if (get_n(wbv) /= size(a,1) .or. get_n(wbv) /= size(a,2)) then
       call set_error(error, 1, id_c_wbv_to_general); return
    end if
    call f_c_wbv_to_general(wbv%br, get_n(wbv), wbv%lbw, wbv%ubw, get_lbwmax(wbv), get_ubwmax(wbv), &
         wbv%numrotsw, wbv%jsw, wbv%csw, wbv%ssw, &
         wbv%numrotsv, wbv%ksv, wbv%csv, wbv%ssv, a)
  end subroutine c_wbv_to_general

  subroutine f_c_wbv_to_general(br, n, lbw, ubw, lbwmax, ubwmax, numrotsw, jsw, csw, ssw, &
       numrotsv, ksv, csv, ssv, a)
    complex(kind=dp), target, dimension(n,n), intent(out) :: a
    integer(kind=int32), dimension(lbwmax,n), intent(in) :: jsw
    complex(kind=dp), dimension(lbwmax,n), intent(in) :: csw, ssw
    integer(kind=int32), dimension(n,ubwmax), intent(in) :: ksv
    complex(kind=dp), dimension(n,ubwmax), intent(in) :: csv, ssv
    complex(kind=dp), dimension(n,ubwmax+lbwmax+1), intent(in) :: br
    integer(kind=int32), dimension(n), intent(in) :: numrotsw, numrotsv
    integer(kind=int32), intent(in) :: ubw, lbw, n, lbwmax, ubwmax
    !
    integer(kind=int32) :: j,k
    type(c_rotation) :: rot
    call br_to_general(br,lbw,ubw,a)
    if (n==1) then
       return
    end if

    do j=1,n-2
       do k=1,numrotsv(j)
          rot%cosine=csv(j,k); rot%sine=ssv(j,k)
          call general_times_rotation(a(1:j,:),trp_rot(rot),ksv(j,k), ksv(j,k)+1)
       end do
    end do

    do k=1,n-2
       do j=1,numrotsw(k)
          rot%cosine=csw(j,k); rot%sine=ssw(j,k)
          call rotation_times_general(rot,a(:,1:k),jsw(j,k),jsw(j,k)+1)
       end do
    end do
  end subroutine f_c_wbv_to_general

end module assemble
