module general_bt
  use orth
  use shift
  use rotation
  use types
  use general_ub
  use misc
  implicit none
  integer(kind=int32), private, parameter :: nullmaxits=5

  interface f_general_bt
     module procedure f_d_general_bt, f_c_general_bt
  end interface f_general_bt

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

  type(routine_info), parameter :: info_f_c_general_bt=routine_info(id_f_c_general_bt, &
       'f_c_general_bt', &
       [ character(len=error_message_length) :: 'Insufficient lower bandwidth.' ])

contains

  subroutine f_d_general_bt(a, n, lbws, lbwmax, numrotst, kst, cst, sst, tol, error)
    real(kind=dp), target, dimension(n,n), intent(inout) :: a
    integer(kind=int32), dimension(lbwmax,n), intent(out) :: kst
    real(kind=dp), dimension(lbwmax,n), intent(out) :: cst, sst
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
    sst=-transpose(ssu)

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

  subroutine f_c_general_bt(a, n, lbws, lbwmax, numrotst, kst, cst, sst, tol, error)
    complex(kind=dp), target, dimension(n,n), intent(inout) :: a
    integer(kind=int32), dimension(lbwmax,n), intent(out) :: kst
    complex(kind=dp), dimension(lbwmax,n), intent(out) :: cst, sst
    real(kind=dp), intent(in) :: tol
    integer(kind=int32), dimension(n), intent(out) :: numrotst
    integer(kind=int32), dimension(n), intent(out) :: lbws
    type(error_info), intent(out) :: error
    integer(kind=int32), intent(in) :: n, lbwmax
    !
    integer(kind=int32), dimension(lbwmax,n) :: jsu
    complex(kind=dp), dimension(lbwmax,n) :: csu, ssu

    call ip_transpose(a)
    call f_c_general_ub(a,n,lbws,lbwmax,numrotst,jsu,csu,ssu,tol,error)
    call ip_transpose(a)
    
    kst=transpose(jsu)
    cst=transpose(csu)
    sst=-conjg(transpose(ssu))

  end subroutine f_c_general_bt

  ! Errors:
  ! 0: no error
  ! 1: insufficient lower bw in bt%bc
  subroutine f_c_lower_to_bt(a, n, b, lbw, ubw, lbwmax, ubwmax, &
       numrotst, kst, cst, sst, tol, error)
    complex(kind=dp), target, dimension(n,n), intent(inout) :: a
    integer(kind=int32), dimension(n,lbwmax), intent(out) :: kst
    complex(kind=dp), dimension(n,lbwmax), intent(out) :: cst, sst
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
       sst=0.0_dp; cst=0.0_dp; kst=0
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

end module general_bt
