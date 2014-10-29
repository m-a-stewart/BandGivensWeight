module mod_general_ubt
  use mod_prec
  use mod_error_id
  use mod_utility
  use mod_nested_types
  use mod_band_types
  use mod_general_ub
  implicit none

  private

  public :: f_general_ubt, f_d_general_ubt, f_c_general_ubt, &
       f_general_to_ubt, f_d_general_to_ubt, f_c_general_to_ubt, &
       general_to_ubt, d_general_to_ubt, c_general_to_ubt

  public :: info_d_general_to_ubt, info_f_d_general_to_ubt, info_c_general_to_ubt, &
       info_f_c_general_to_ubt

  interface f_general_ubt
     module procedure f_d_general_ubt, f_c_general_ubt
  end interface f_general_ubt

  interface f_general_to_ubt
     module procedure f_d_general_to_ubt, f_c_general_to_ubt
  end interface f_general_to_ubt

  interface general_to_ubt
     module procedure d_general_to_ubt, c_general_to_ubt
  end interface general_to_ubt

  type(routine_info), parameter :: info_d_general_to_ubt=routine_info(id_d_general_to_ubt, &
       'd_general_to_ubt', &
       [ character(len=error_message_length) :: 'n<1', &
       'Size of a and ubt not the same.' ] )

  type(routine_info), parameter :: info_f_d_general_to_ubt=routine_info(id_f_d_general_to_ubt, &
       'f_d_general_to_ubt', &
       [ character(len=error_message_length) :: 'Insufficient lower bandwidth in ubt', &
       'Insufficient upper bandwidth in ubt' ])

  type(routine_info), parameter :: info_c_general_to_ubt=routine_info(id_c_general_to_ubt, &
       'c_general_to_ubt', &
       [ character(len=error_message_length) :: 'n<1', &
       'Size of a and ubt not the same.' ] )

  type(routine_info), parameter :: info_f_c_general_to_ubt=routine_info(id_f_c_general_to_ubt, &
       'f_c_general_to_ubt', &
       [ character(len=error_message_length) :: 'Insufficient lower bandwidth in ubt', &
       'Insufficient upper bandwidth in ubt' ])

contains

  subroutine f_d_general_ubt(a, n, lbws, ubws, lbwmax, ubwmax, &
       numrotsu, jsu, csu, ssu, &
       numrotst, kst, cst, sst, tol, error)
    real(kind=dp), target, dimension(n,n), intent(inout) :: a
    integer(kind=int32), dimension(n,lbwmax), intent(out) :: kst
    real(kind=dp), dimension(n,lbwmax), intent(out) :: cst, sst
    integer(kind=int32), dimension(ubwmax,n), intent(out) :: jsu
    real(kind=dp), dimension(ubwmax,n), intent(out) :: csu, ssu
    real(kind=dp), intent(in) :: tol
    integer(kind=int32), dimension(n), intent(out) :: numrotst, numrotsu
    integer(kind=int32), dimension(n), intent(out) :: lbws, ubws
    type(error_info), intent(out) :: error
    integer(kind=int32), intent(in) :: n, lbwmax, ubwmax
    !
    integer(kind=int32), dimension(lbwmax,n) :: jsu0
    real(kind=dp), dimension(lbwmax,n) :: csu0, ssu0

    call f_d_general_ub(a,n,ubws,ubwmax,numrotsu,jsu,csu,ssu,tol,error)
    call ip_transpose(a)
    call f_d_general_ub(a,n,lbws,lbwmax,numrotst,jsu0,csu0,ssu0,tol,error)
    call ip_transpose(a)
    
    kst=transpose(jsu0)
    cst=transpose(csu0)
    sst=transpose(ssu0)

  end subroutine f_d_general_ubt

  ! Errors:
  ! 0: no error
  ! 1: insufficient lower bw in ubt.
  ! 2: insufficient upper bw in ubt.
  subroutine f_d_general_to_ubt(a, n, bc, lbw, ubw, lbwmax, ubwmax, &
       numrotsu, jsu, csu, ssu, &
       numrotst, kst, cst, sst, tol, error)
    real(kind=dp), target, dimension(n,n), intent(inout) :: a
    integer(kind=int32), dimension(n,lbwmax), intent(out) :: kst
    real(kind=dp), dimension(n,lbwmax), intent(out) :: cst, sst
    integer(kind=int32), dimension(ubwmax,n), intent(out) :: jsu
    real(kind=dp), dimension(ubwmax,n), intent(out) :: csu, ssu
    real(kind=dp), dimension(ubwmax+lbwmax+1,n), intent(out) :: bc
    real(kind=dp), intent(in) :: tol
    integer(kind=int32), dimension(n), intent(out) :: numrotst, numrotsu
    integer(kind=int32), intent(out) :: lbw, ubw
    type(error_info), intent(out) :: error
    integer(kind=int32), intent(in) :: n, ubwmax, lbwmax
    !
    integer(kind=int32), dimension(n) :: lbws, ubws

    call clear_error(error)

    if (n == 1) then
       numrotst=0;
       sst=0.0_dp; cst=0.0_dp; kst=0
       lbw=0
       numrotsu=0;
       ssu=0.0_dp; csu=0.0_dp; kst=0
       ubw=0
       bc(1,1)=a(1,1)
       return
    end if

    call f_d_general_ubt(a, n, lbws, ubws, lbwmax, ubwmax, &
         numrotsu, jsu, csu, ssu, &
         numrotst, kst, cst, sst, tol, error)

    lbw=maxval(lbws)
    ubw=maxval(ubws)
    if (lbw > lbwmax) then
       ! This should already have been detected in f_d_general_bt.
       call set_error(error, 1, id_f_d_general_to_ubt)
    else if (ubw > ubwmax) then
       call set_error(error, 2, id_f_d_general_to_ubt)
    else
       call d_extract_diagonals_bc(a, n, bc, lbw, ubw, lbwmax, ubwmax)
    end if
  end subroutine f_d_general_to_ubt

  ! Errors:
  ! 1: n < 1
  ! 2: n is not the same for a and ubt.
  subroutine d_general_to_ubt(a,ubt,tol,error)
    real(kind=dp), target, dimension(:,:), intent(inout) :: a
    type(d_ubt), intent(inout) :: ubt
    type(error_info), intent(out) :: error
    real(kind=dp), intent(in) :: tol
    call clear_error(error)
    if (size(a,1) < 1) then
       call set_error(error, 1, id_d_general_to_ubt); return
    end if
    if (get_n(ubt) /= size(a,1) .or. get_n(ubt) /= size(a,2)) then
       call set_error(error, 2, id_d_general_to_ubt); return
    end if
    call f_d_general_to_ubt(a,get_n(ubt),ubt%bc, ubt%lbw, ubt%ubw, get_lbwmax(ubt), &
         get_ubwmax(ubt), ubt%numrotsu, ubt%jsu, ubt%csu, ubt%ssu, &
         ubt%numrotst, ubt%kst, ubt%cst, ubt%sst, tol, error)
  end subroutine d_general_to_ubt

  subroutine f_c_general_ubt(a, n, lbws, ubws, lbwmax, ubwmax, &
       numrotsu, jsu, csu, ssu, &
       numrotst, kst, cst, sst, tol, error)
    complex(kind=dp), target, dimension(n,n), intent(inout) :: a
    integer(kind=int32), dimension(n,lbwmax), intent(out) :: kst
    complex(kind=dp), dimension(n,lbwmax), intent(out) :: cst, sst
    integer(kind=int32), dimension(ubwmax,n), intent(out) :: jsu
    complex(kind=dp), dimension(ubwmax,n), intent(out) :: csu, ssu
    real(kind=dp), intent(in) :: tol
    integer(kind=int32), dimension(n), intent(out) :: numrotst, numrotsu
    integer(kind=int32), dimension(n), intent(out) :: lbws, ubws
    type(error_info), intent(out) :: error
    integer(kind=int32), intent(in) :: n, lbwmax, ubwmax
    !
    integer(kind=int32), dimension(lbwmax,n) :: jsu0
    complex(kind=dp), dimension(lbwmax,n) :: csu0, ssu0

    call f_c_general_ub(a,n,ubws,ubwmax,numrotsu,jsu,csu,ssu,tol,error)
    call ip_transpose(a)
    call f_c_general_ub(a,n,lbws,lbwmax,numrotst,jsu0,csu0,ssu0,tol,error)
    call ip_transpose(a)
    
    kst=transpose(jsu0)
    cst=transpose(csu0)
    sst=transpose(ssu0)

  end subroutine f_c_general_ubt

  ! Errors:
  ! 0: no error
  ! 1: insufficient lower bw in ubt.
  ! 2: insufficient upper bw in ubt.
  subroutine f_c_general_to_ubt(a, n, bc, lbw, ubw, lbwmax, ubwmax, &
       numrotsu, jsu, csu, ssu, &
       numrotst, kst, cst, sst, tol, error)
    complex(kind=dp), target, dimension(n,n), intent(inout) :: a
    integer(kind=int32), dimension(n,lbwmax), intent(out) :: kst
    complex(kind=dp), dimension(n,lbwmax), intent(out) :: cst, sst
    integer(kind=int32), dimension(ubwmax,n), intent(out) :: jsu
    complex(kind=dp), dimension(ubwmax,n), intent(out) :: csu, ssu
    complex(kind=dp), dimension(ubwmax+lbwmax+1,n), intent(out) :: bc
    real(kind=dp), intent(in) :: tol
    integer(kind=int32), dimension(n), intent(out) :: numrotst, numrotsu
    integer(kind=int32), intent(out) :: lbw, ubw
    type(error_info), intent(out) :: error
    integer(kind=int32), intent(in) :: n, ubwmax, lbwmax
    !
    integer(kind=int32), dimension(n) :: lbws, ubws

    call clear_error(error)

    if (n == 1) then
       numrotst=0;
       sst=0.0_dp; cst=0.0_dp; kst=0
       lbw=0
       numrotsu=0;
       ssu=0.0_dp; csu=0.0_dp; kst=0
       ubw=0
       bc(1,1)=a(1,1)
       return
    end if

    call f_c_general_ubt(a, n, lbws, ubws, lbwmax, ubwmax, &
         numrotsu, jsu, csu, ssu, &
         numrotst, kst, cst, sst, tol, error)

    lbw=maxval(lbws)
    ubw=maxval(ubws)
    if (lbw > lbwmax) then
       ! This should already have been detected in f_c_general_bt.
       call set_error(error, 1, id_f_c_general_to_ubt)
    else if (ubw > ubwmax) then
       call set_error(error, 2, id_f_c_general_to_ubt)
    else
       call c_extract_diagonals_bc(a, n, bc, lbw, ubw, lbwmax, ubwmax)
    end if
  end subroutine f_c_general_to_ubt

  ! Errors:
  ! 1: n < 1
  ! 2: n is not the same for a and ubt.
  subroutine c_general_to_ubt(a,ubt,tol,error)
    complex(kind=dp), target, dimension(:,:), intent(inout) :: a
    type(c_ubt), intent(inout) :: ubt
    type(error_info), intent(out) :: error
    real(kind=dp), intent(in) :: tol
    call clear_error(error)
    if (size(a,1) < 1) then
       call set_error(error, 1, id_c_general_to_ubt); return
    end if
    if (get_n(ubt) /= size(a,1) .or. get_n(ubt) /= size(a,2)) then
       call set_error(error, 2, id_c_general_to_ubt); return
    end if
    call f_c_general_to_ubt(a,get_n(ubt),ubt%bc, ubt%lbw, ubt%ubw, get_lbwmax(ubt), &
         get_ubwmax(ubt), ubt%numrotsu, ubt%jsu, ubt%csu, ubt%ssu, &
         ubt%numrotst, ubt%kst, ubt%cst, ubt%sst, tol, error)
  end subroutine c_general_to_ubt


end module mod_general_ubt
