module mod_general_wb
  use mod_prec
  use mod_error_id
  use mod_utility
  use mod_orth_band_types
  use mod_band_types
  use mod_general_bv
  implicit none

  private
  public :: f_general_wb, f_d_general_wb, f_c_general_wb, &
       f_lower_to_wb, f_d_lower_to_wb, f_c_lower_to_wb, &
       lower_to_wb, d_lower_to_wb, c_lower_to_wb

  public :: info_d_lower_to_wb, info_f_d_lower_to_wb, info_f_d_general_wb, &
       info_c_lower_to_wb, info_f_c_lower_to_wb, info_f_c_general_wb

  interface f_general_wb
     module procedure f_d_general_wb, f_c_general_wb
  end interface f_general_wb

  interface f_lower_to_wb
     module procedure f_d_lower_to_wb, f_c_lower_to_wb
  end interface f_lower_to_wb

  interface lower_to_wb
     module procedure d_lower_to_wb, c_lower_to_wb
  end interface lower_to_wb

  type(routine_info), parameter :: info_d_lower_to_wb=routine_info(id_d_lower_to_wb, &
       'd_lower_to_wb', &
       [ character(len=error_message_length) :: 'n<1', 'bt%ubwmax < ubw', 'Size of a and bt not the same.' ] )

  type(routine_info), parameter :: info_f_d_lower_to_wb=routine_info(id_f_d_lower_to_wb, &
       'f_d_lower_to_wb', &
       [ character(len=error_message_length) :: 'Insufficient lower bandwidth in bt' ])

  type(routine_info), parameter :: info_f_d_general_wb=routine_info(id_f_d_general_wb, &
       'f_d_general_wb', &
       [ character(len=error_message_length) :: 'Insufficient lower bandwidth.' ])

  type(routine_info), parameter :: info_c_lower_to_wb=routine_info(id_c_lower_to_wb, &
       'c_lower_to_wb', &
       [ character(len=error_message_length) :: 'n<1', 'bt%ubwmax < ubw', 'Size of a and bt not the same.' ] )

  type(routine_info), parameter :: info_f_c_lower_to_wb=routine_info(id_f_c_lower_to_wb, &
       'f_c_lower_to_wb', &
       [ character(len=error_message_length) :: 'Insufficient lower bandwidth in bt' ])

  type(routine_info), parameter :: info_f_c_general_wb=routine_info(id_f_c_general_wb, &
       'f_c_general_wb', &
       [ character(len=error_message_length) :: 'Insufficient lower bandwidth.' ])

contains

  subroutine f_d_general_wb(a, n, lbws, lbwmax, numrotsw, jsw, csw, ssw, tol, error)
    real(kind=dp), target, dimension(n,n), intent(inout) :: a
    integer(kind=int32), dimension(lbwmax,n), intent(out) :: jsw
    real(kind=dp), dimension(lbwmax,n), intent(out) :: csw, ssw
    real(kind=dp), intent(in) :: tol
    integer(kind=int32), dimension(n), intent(out) :: numrotsw
    integer(kind=int32), dimension(n), intent(out) :: lbws
    type(error_info), intent(out) :: error
    integer(kind=int32), intent(in) :: n, lbwmax
    !
    integer(kind=int32), dimension(n,lbwmax) :: ksv
    real(kind=dp), dimension(n,lbwmax) :: csv, ssv

    call ip_transpose(a)
    call f_d_general_bv(a,n,lbws,lbwmax,numrotsw,ksv,csv,ssv,tol,error)
    call ip_transpose(a)
    
    jsw=transpose(ksv)
    csw=transpose(csv)
    ssw=transpose(ssv)

  end subroutine f_d_general_wb

  ! Errors:
  ! 0: no error
  ! 1: insufficient lower bw in wb%bc
  subroutine f_d_lower_to_wb(a, n, b, lbw, ubw, lbwmax, ubwmax, &
       numrotsw, jsw, csw, ssw, tol, error)
    real(kind=dp), target, dimension(n,n), intent(inout) :: a
    integer(kind=int32), dimension(lbwmax,n), intent(out) :: jsw
    real(kind=dp), dimension(lbwmax,n), intent(out) :: csw, ssw
    real(kind=dp), dimension(ubwmax+lbwmax+1,n), intent(out) :: b
    real(kind=dp), intent(in) :: tol
    integer(kind=int32), dimension(n), intent(out) :: numrotsw
    integer(kind=int32), intent(out) :: lbw
    type(error_info), intent(out) :: error
    integer(kind=int32), intent(in) :: n, ubw, ubwmax, lbwmax
    !
    integer(kind=int32), dimension(n) :: lbws

    call clear_error(error)

    if (n == 1) then
       numrotsw=0;
       ssw=0.0_dp; csw=0.0_dp; jsw=0
       lbw=0
       b(1,1)=a(1,1)
       return
    end if

    call f_d_general_wb(a, n, lbws, lbwmax, numrotsw, jsw, csw, ssw, tol, error)

    lbw=maxval(lbws)
    if (maxval(lbws) > lbwmax) then
       ! This should already have been detected in f_d_general_wb.
       call set_error(error, 1, id_f_d_lower_to_wb)
    else
       call d_extract_diagonals_bc(a, n, b, lbw, ubw, lbwmax, ubwmax)
    end if
  end subroutine f_d_lower_to_wb

  ! Errors:
  ! 0: no error
  ! 1: n < 1
  ! 2: wb%ubwmax < ubw
  ! 3: n is not the same for a and wb.
  subroutine d_lower_to_wb(a,wb,ubw,tol,error)
    real(kind=dp), target, dimension(:,:), intent(inout) :: a
    type(d_wb), intent(inout) :: wb
    type(error_info), intent(out) :: error
    real(kind=dp), intent(in) :: tol
    integer(kind=int32), intent(in) :: ubw
    call clear_error(error)
    if (size(a,1) < 1) then
       call set_error(error, 1, id_d_lower_to_wb); return
    end if
    if (get_ubwmax(wb) < ubw) then
       call set_error(error, 2, id_d_lower_to_wb); return
    end if
    if (get_n(wb) /= size(a,1) .or. get_n(wb) /= size(a,2)) then
       call set_error(error, 3, id_d_lower_to_wb); return
    end if
    call f_d_lower_to_wb(a,get_n(wb),wb%bc, wb%lbw, ubw, get_lbwmax(wb), get_ubwmax(wb), &
         wb%numrotsw, wb%jsw, wb%csw, wb%ssw, tol, error)
    wb%ubw=ubw
  end subroutine d_lower_to_wb

  subroutine f_c_general_wb(a, n, lbws, lbwmax, numrotsw, jsw, csw, ssw, tol, error)
    complex(kind=dp), target, dimension(n,n), intent(inout) :: a
    integer(kind=int32), dimension(lbwmax,n), intent(out) :: jsw
    real(kind=dp), dimension(lbwmax,n), intent(out) :: csw
    complex(kind=dp), dimension(lbwmax,n), intent(out) :: ssw
    real(kind=dp), intent(in) :: tol
    integer(kind=int32), dimension(n), intent(out) :: numrotsw
    integer(kind=int32), dimension(n), intent(out) :: lbws
    type(error_info), intent(out) :: error
    integer(kind=int32), intent(in) :: n, lbwmax
    !
    integer(kind=int32), dimension(n,lbwmax) :: ksv
    real(kind=dp), dimension(n,lbwmax) :: csv
    complex(kind=dp), dimension(n,lbwmax) :: ssv

    call ip_transpose(a)
    call f_c_general_bv(a,n,lbws,lbwmax,numrotsw,ksv,csv,ssv,tol,error)
    call ip_transpose(a)
    
    jsw=transpose(ksv)
    csw=transpose(csv)
    ssw=transpose(ssv)

  end subroutine f_c_general_wb

  ! Errors:
  ! 0: no error
  ! 1: insufficient lower bw in wb%bc
  subroutine f_c_lower_to_wb(a, n, b, lbw, ubw, lbwmax, ubwmax, &
       numrotsw, jsw, csw, ssw, tol, error)
    complex(kind=dp), target, dimension(n,n), intent(inout) :: a
    integer(kind=int32), dimension(lbwmax,n), intent(out) :: jsw
    real(kind=dp), dimension(lbwmax,n), intent(out) :: csw
    complex(kind=dp), dimension(lbwmax,n), intent(out) :: ssw
    complex(kind=dp), dimension(ubwmax+lbwmax+1,n), intent(out) :: b
    real(kind=dp), intent(in) :: tol
    integer(kind=int32), dimension(n), intent(out) :: numrotsw
    integer(kind=int32), intent(out) :: lbw
    type(error_info), intent(out) :: error
    integer(kind=int32), intent(in) :: n, ubw, ubwmax, lbwmax
    !
    integer(kind=int32), dimension(n) :: lbws

    call clear_error(error)

    if (n == 1) then
       numrotsw=0;
       ssw=(0.0_dp, 0.0_dp); csw=0.0_dp; jsw=0
       lbw=0
       b(1,1)=a(1,1)
       return
    end if

    call f_c_general_wb(a, n, lbws, lbwmax, numrotsw, jsw, csw, ssw, tol, error)

    lbw=maxval(lbws)
    if (maxval(lbws) > lbwmax) then
       ! This should already have been detected in f_c_general_wb.
       call set_error(error, 1, id_f_c_lower_to_wb)
    else
       call c_extract_diagonals_bc(a, n, b, lbw, ubw, lbwmax, ubwmax)
    end if
  end subroutine f_c_lower_to_wb

  ! Errors:
  ! 0: no error
  ! 1: n < 1
  ! 2: wb%ubwmax < ubw
  ! 3: n is not the same for a and wb.
  subroutine c_lower_to_wb(a,wb,ubw,tol,error)
    complex(kind=dp), target, dimension(:,:), intent(inout) :: a
    type(c_wb), intent(inout) :: wb
    type(error_info), intent(out) :: error
    real(kind=dp), intent(in) :: tol
    integer(kind=int32), intent(in) :: ubw
    call clear_error(error)
    if (size(a,1) < 1) then
       call set_error(error, 1, id_c_lower_to_wb); return
    end if
    if (get_ubwmax(wb) < ubw) then
       call set_error(error, 2, id_c_lower_to_wb); return
    end if
    if (get_n(wb) /= size(a,1) .or. get_n(wb) /= size(a,2)) then
       call set_error(error, 3, id_c_lower_to_wb); return
    end if
    call f_c_lower_to_wb(a,get_n(wb),wb%bc, wb%lbw, ubw, get_lbwmax(wb), get_ubwmax(wb), &
         wb%numrotsw, wb%jsw, wb%csw, wb%ssw, tol, error)
    wb%ubw=ubw
  end subroutine c_lower_to_wb


end module mod_general_wb
