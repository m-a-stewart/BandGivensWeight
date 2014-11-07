module mod_general_bt
  use mod_prec
  use mod_error_id
  use mod_utility
  use mod_nested_types
  use mod_band_types
  use mod_general_ub
  implicit none

  private

  public :: f_general_bt, f_d_general_bt, f_c_general_bt, &
       f_lower_to_bt, f_d_lower_to_bt, f_c_lower_to_bt, &
       lower_to_bt, d_lower_to_bt, c_lower_to_bt

  public :: info_d_lower_to_bt, info_f_d_lower_to_bt, info_f_d_general_bt, &
       info_c_lower_to_bt, info_f_c_lower_to_bt

  interface f_general_bt
     module procedure f_d_general_bt, f_c_general_bt
  end interface f_general_bt

  interface f_lower_to_bt
     module procedure f_d_lower_to_bt, f_c_lower_to_bt
  end interface f_lower_to_bt

  interface lower_to_bt
     module procedure d_lower_to_bt, c_lower_to_bt
  end interface lower_to_bt

  type(routine_info), parameter :: info_d_lower_to_bt=routine_info(id_d_lower_to_bt, &
       'd_lower_to_bt', &
       [ character(len=error_message_length) :: 'n<1', 'bt%ubwmax < ubw', 'Size of a and bt not the same.' ] )

  type(routine_info), parameter :: info_f_d_lower_to_bt=routine_info(id_f_d_lower_to_bt, &
       'f_d_lower_to_bt', &
       [ character(len=error_message_length) :: 'Insufficient lower bandwidth in bt' ])

  type(routine_info), parameter :: info_f_d_general_bt=routine_info(id_f_d_general_bt, &
       'f_d_general_bt', &
       [ character(len=error_message_length) :: 'Insufficient lower bandwidth.' ])

  type(routine_info), parameter :: info_c_lower_to_bt=routine_info(id_c_lower_to_bt, &
       'c_lower_to_bt', &
       [ character(len=error_message_length) :: 'n<1', 'bt%ubwmax < ubw', 'Size of a and bt not the same.' ] )

  type(routine_info), parameter :: info_f_c_lower_to_bt=routine_info(id_f_c_lower_to_bt, &
       'f_c_lower_to_bt', &
       [ character(len=error_message_length) :: 'Insufficient lower bandwidth in bt' ])

contains

  subroutine f_d_general_bt(a, n, lbws, lbwmax, numrotst, kst, cst, sst, tol, error)
    real(kind=dp), target, dimension(n,n), intent(inout) :: a
    integer(kind=int32), dimension(n,lbwmax), intent(out) :: kst
    real(kind=dp), dimension(n,lbwmax), intent(out) :: cst, sst
    real(kind=dp), intent(in) :: tol
    integer(kind=int32), dimension(n), intent(out) :: numrotst
    integer(kind=int32), dimension(n), intent(out) :: lbws
    type(error_info), intent(out) :: error
    integer(kind=int32), intent(in) :: n, lbwmax
    !
    integer(kind=int32), dimension(lbwmax,n) :: jsu
    real(kind=dp), dimension(lbwmax,n) :: csu, ssu

    call ip_transpose(a)
    call f_d_general_ub(a,n,lbws,lbwmax,numrotst,jsu,csu,ssu,tol,error)
    call ip_transpose(a)
    
    kst=transpose(jsu)
    cst=transpose(csu)
    sst=transpose(ssu)

  end subroutine f_d_general_bt

  ! Errors:
  ! 0: no error
  ! 1: insufficient lower bw in bt%bc
  subroutine f_d_lower_to_bt(a, n, b, lbw, ubw, lbwmax, ubwmax, &
       numrotst, kst, cst, sst, tol, error)
    real(kind=dp), target, dimension(n,n), intent(inout) :: a
    integer(kind=int32), dimension(n,lbwmax), intent(out) :: kst
    real(kind=dp), dimension(n,lbwmax), intent(out) :: cst, sst
    real(kind=dp), dimension(n,ubwmax+lbwmax+1), intent(out) :: b
    real(kind=dp), intent(in) :: tol
    integer(kind=int32), dimension(n), intent(out) :: numrotst
    integer(kind=int32), intent(out) :: lbw
    type(error_info), intent(out) :: error
    integer(kind=int32), intent(in) :: n, ubw, ubwmax, lbwmax
    !
    integer(kind=int32), dimension(n) :: lbws

    call clear_error(error)

    if (n == 1) then
       numrotst=0;
       sst=0.0_dp; cst=0.0_dp; kst=0
       lbw=0
       b(1,1)=a(1,1)
       return
    end if

    call f_d_general_bt(a, n, lbws, lbwmax, numrotst, kst, cst, sst, tol, error)
    
    lbw=maxval(lbws)
    if (maxval(lbws) > lbwmax) then
       ! This should already have been detected in f_d_general_bt.
       call set_error(error, 1, id_f_d_lower_to_bt)
    else
       call d_extract_diagonals_br(a, n, b, lbw, ubw, lbwmax, ubwmax)
    end if
  end subroutine f_d_lower_to_bt

  ! Errors:
  ! 0: no error
  ! 1: n < 1
  ! 2: bt%lbwmax < lbw
  ! 3: n is not the same for a and ub or bv.
  subroutine d_lower_to_bt(a,bt,ubw,tol,error)
    real(kind=dp), target, dimension(:,:), intent(inout) :: a
    type(d_bt), intent(inout) :: bt
    type(error_info), intent(out) :: error
    real(kind=dp), intent(in) :: tol
    integer(kind=int32), intent(in) :: ubw
    call clear_error(error)
    if (size(a,1) < 1) then
       call set_error(error, 1, id_d_lower_to_bt); return
    end if
    if (get_ubwmax(bt) < ubw) then
       call set_error(error, 2, id_d_lower_to_bt); return
    end if
    if (get_n(bt) /= size(a,1) .or. get_n(bt) /= size(a,2)) then
       call set_error(error, 3, id_d_lower_to_bt); return
    end if
    call f_d_lower_to_bt(a,get_n(bt),bt%br, bt%lbw, ubw, get_lbwmax(bt), get_ubwmax(bt), &
         bt%numrotst, bt%kst, bt%cst, bt%sst, tol, error)
    bt%ubw=ubw
  end subroutine d_lower_to_bt


  subroutine f_c_general_bt(a, n, lbws, lbwmax, numrotst, kst, cst, sst, tol, error)
    complex(kind=dp), target, dimension(n,n), intent(inout) :: a
    integer(kind=int32), dimension(n,lbwmax), intent(out) :: kst
    real(kind=dp), dimension(n,lbwmax), intent(out) :: cst
    complex(kind=dp), dimension(n,lbwmax), intent(out) :: sst
    real(kind=dp), intent(in) :: tol
    integer(kind=int32), dimension(n), intent(out) :: numrotst
    integer(kind=int32), dimension(n), intent(out) :: lbws
    type(error_info), intent(out) :: error
    integer(kind=int32), intent(in) :: n, lbwmax
    !
    integer(kind=int32), dimension(lbwmax,n) :: jsu
    real(kind=dp), dimension(lbwmax,n) :: csu
    complex(kind=dp), dimension(lbwmax,n) :: ssu

    call ip_transpose(a)
    call f_c_general_ub(a,n,lbws,lbwmax,numrotst,jsu,csu,ssu,tol,error)
    call ip_transpose(a)
    
    kst=transpose(jsu)
    cst=transpose(csu)
    sst=transpose(ssu)

  end subroutine f_c_general_bt

  ! Errors:
  ! 0: no error
  ! 1: insufficient lower bw in bt%bc
  subroutine f_c_lower_to_bt(a, n, b, lbw, ubw, lbwmax, ubwmax, &
       numrotst, kst, cst, sst, tol, error)
    complex(kind=dp), target, dimension(n,n), intent(inout) :: a
    integer(kind=int32), dimension(n,lbwmax), intent(out) :: kst
    real(kind=dp), dimension(n,lbwmax), intent(out) :: cst
    complex(kind=dp), dimension(n,lbwmax), intent(out) :: sst
    complex(kind=dp), dimension(n,ubwmax+lbwmax+1), intent(out) :: b
    real(kind=dp), intent(in) :: tol
    integer(kind=int32), dimension(n), intent(out) :: numrotst
    integer(kind=int32), intent(out) :: lbw
    type(error_info), intent(out) :: error
    integer(kind=int32), intent(in) :: n, ubw, ubwmax, lbwmax
    !
    integer(kind=int32), dimension(n) :: lbws

    call clear_error(error)

    if (n == 1) then
       numrotst=0;
       sst=(0.0_dp,0.0_dp); cst=0.0_dp; kst=0
       lbw=0
       b(1,1)=a(1,1)
       return
    end if

    call f_c_general_bt(a, n, lbws, lbwmax, numrotst, kst, cst, sst, tol, error)
    
    lbw=maxval(lbws)
    if (maxval(lbws) > lbwmax) then
       ! This should already have been detected in f_c_general_bt.
       call set_error(error, 1, id_f_c_lower_to_bt)
    else
       call c_extract_diagonals_br(a, n, b, lbw, ubw, lbwmax, ubwmax)
    end if
  end subroutine f_c_lower_to_bt

  ! Errors:
  ! 0: no error
  ! 1: n < 1
  ! 2: ub%ubwmax < ubw
  ! 3: n is not the same for a and ub or bv.
  subroutine c_lower_to_bt(a,bt,ubw,tol,error)
    complex(kind=dp), target, dimension(:,:), intent(inout) :: a
    type(c_bt), intent(inout) :: bt
    type(error_info), intent(out) :: error
    real(kind=dp), intent(in) :: tol
    integer(kind=int32), intent(in) :: ubw
    call clear_error(error)
    if (size(a,1) < 1) then
       call set_error(error, 1, id_c_lower_to_bt); return
    end if
    if (get_ubwmax(bt) < ubw) then
       call set_error(error, 2, id_c_lower_to_bt); return
    end if
    if (get_n(bt) /= size(a,1) .or. get_n(bt) /= size(a,2)) then
       call set_error(error, 3, id_c_lower_to_bt); return
    end if
    call f_c_lower_to_bt(a,get_n(bt),bt%br, bt%lbw, ubw, get_lbwmax(bt), get_ubwmax(bt), &
         bt%numrotst, bt%kst, bt%cst, bt%sst, tol, error)
    bt%ubw=ubw
  end subroutine c_lower_to_bt

end module mod_general_bt
