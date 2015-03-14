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
       f_general_to_wb, f_d_general_to_wb, f_c_general_to_wb, &
       general_to_wb, d_general_to_wb, c_general_to_wb, &
       d_wb_of_general, c_wb_of_general, wb_of_general, wb

  interface wb_of_general
     module procedure d_wb_of_general, c_wb_of_general
  end interface wb_of_general

  interface wb
     module procedure d_wb_of_general, c_wb_of_general
  end interface wb
  
  interface f_general_wb
     module procedure f_d_general_wb, f_c_general_wb
  end interface f_general_wb

  interface f_general_to_wb
     module procedure f_d_general_to_wb, f_c_general_to_wb
  end interface f_general_to_wb

  interface general_to_wb
     module procedure d_general_to_wb, c_general_to_wb
  end interface general_to_wb

contains

  function d_wb_of_general(a, ubw, lbwmax, ubwmax, tol, error) result(wb)
    type(d_wb), allocatable :: wb
    real(kind=dp), target, dimension(:,:), intent(inout) :: a
    type(error_info), intent(inout), optional :: error
    real(kind=dp), intent(in) :: tol
    integer(kind=int32), intent(in) :: ubw, lbwmax, ubwmax
    type(routine_info), parameter :: info=info_d_wb_of_general
    integer(kind=int32) :: n

    if (failure(error)) return
    call push_id(info, error)

    n=size(a,1)
    wb=d_new_wb(n,lbwmax,ubwmax)

    call d_general_to_wb(a,wb,ubw,tol,error)

    call pop_id(error)
    
  end function d_wb_of_general

  subroutine f_d_general_wb(a, n, lbws, lbwmax, numrotsw, jsw, csw, ssw, tol, error)
    real(kind=dp), target, dimension(n,n), intent(inout) :: a
    integer(kind=int32), dimension(lbwmax,n), intent(out) :: jsw
    real(kind=dp), dimension(lbwmax,n), intent(out) :: csw, ssw
    real(kind=dp), intent(in) :: tol
    integer(kind=int32), dimension(n), intent(out) :: numrotsw
    integer(kind=int32), dimension(n), intent(out) :: lbws
    type(error_info), intent(inout), optional :: error
    integer(kind=int32), intent(in) :: n, lbwmax
    !
    integer(kind=int32), dimension(n,lbwmax) :: ksv
    real(kind=dp), dimension(n,lbwmax) :: csv, ssv
    type(routine_info), parameter :: info=info_f_d_general_wb

    if (failure(error)) return
    call push_id(info, error)

    call ip_transpose(a)
    call f_d_general_bv(a,n,lbws,lbwmax,numrotsw,ksv,csv,ssv,tol,error)

    if (success(error)) then

       call ip_transpose(a)

       jsw=transpose(ksv)
       csw=transpose(csv)
       ssw=transpose(ssv)
    end if

    call pop_id(error)

  end subroutine f_d_general_wb

  ! Errors:
  ! 0: no error
  ! 1: insufficient lower bw in wb%bc
  subroutine f_d_general_to_wb(a, n, b, lbw, ubw, lbwmax, ubwmax, &
       numrotsw, jsw, csw, ssw, tol, error)
    real(kind=dp), target, dimension(n,n), intent(inout) :: a
    integer(kind=int32), dimension(lbwmax,n), intent(out) :: jsw
    real(kind=dp), dimension(lbwmax,n), intent(out) :: csw, ssw
    real(kind=dp), dimension(ubwmax+lbwmax+1,n), intent(out) :: b
    real(kind=dp), intent(in) :: tol
    integer(kind=int32), dimension(n), intent(out) :: numrotsw
    integer(kind=int32), intent(out) :: lbw
    type(error_info), intent(inout), optional :: error
    integer(kind=int32), intent(in) :: n, ubw, ubwmax, lbwmax
    !
    integer(kind=int32), dimension(n) :: lbws
    type(routine_info), parameter :: info=info_f_d_general_to_wb
    
    if (failure(error)) return
    call push_id(info, error)

    if (n == 1) then
       numrotsw=0;
       ssw=0.0_dp; csw=0.0_dp; jsw=0
       lbw=0
       b(1,1)=a(1,1)
       return
    end if

    call f_d_general_wb(a, n, lbws, lbwmax, numrotsw, jsw, csw, ssw, tol, error)

    if (success(error)) then
       lbw=maxval(lbws)
       if (maxval(lbws) > lbwmax) then
          ! This should already have been detected in f_d_general_wb.
          call set_error(1, info, error); return
       else
          call d_extract_diagonals_bc(a, n, b, lbw, ubw, lbwmax, ubwmax)
       end if
    end if

    call pop_id(error)
  end subroutine f_d_general_to_wb

  ! Errors:
  ! 0: no error
  ! 1: n < 1
  ! 2: wb%ubwmax < ubw
  ! 3: n is not the same for a and wb.
  subroutine d_general_to_wb(a,wb,ubw,tol,error)
    real(kind=dp), target, dimension(:,:), intent(inout) :: a
    type(d_wb), intent(inout) :: wb
    type(error_info), intent(inout), optional :: error
    real(kind=dp), intent(in) :: tol
    integer(kind=int32), intent(in) :: ubw
    type(routine_info), parameter :: info=info_d_general_to_wb

    if (failure(error)) return
    call push_id(info, error)
    
    if (size(a,1) < 1) then
       call set_error(1, info, error); return
    end if
    if (get_ubwmax(wb) < ubw) then
       call set_error(2, info, error); return
    end if
    if (get_n(wb) /= size(a,1) .or. get_n(wb) /= size(a,2)) then
       call set_error(3, info, error); return
    end if
    call f_d_general_to_wb(a,get_n(wb),wb%bc, wb%lbw, ubw, get_lbwmax(wb), get_ubwmax(wb), &
         wb%numrotsw, wb%jsw, wb%csw, wb%ssw, tol, error)

    wb%ubw=ubw
    call pop_id(error)
  end subroutine d_general_to_wb

  function c_wb_of_general(a, ubw, lbwmax, ubwmax, tol, error) result(wb)
    type(c_wb), allocatable :: wb
    complex(kind=dp), target, dimension(:,:), intent(inout) :: a
    type(error_info), intent(inout), optional :: error
    real(kind=dp), intent(in) :: tol
    integer(kind=int32), intent(in) :: ubw, lbwmax, ubwmax
    type(routine_info), parameter :: info=info_c_wb_of_general
    integer(kind=int32) :: n

    if (failure(error)) return
    call push_id(info, error)

    n=size(a,1)
    wb=c_new_wb(n,lbwmax,ubwmax)

    call c_general_to_wb(a,wb,ubw,tol,error)

    call pop_id(error)
    
  end function c_wb_of_general
  

  subroutine f_c_general_wb(a, n, lbws, lbwmax, numrotsw, jsw, csw, ssw, tol, error)
    complex(kind=dp), target, dimension(n,n), intent(inout) :: a
    integer(kind=int32), dimension(lbwmax,n), intent(out) :: jsw
    real(kind=dp), dimension(lbwmax,n), intent(out) :: csw
    complex(kind=dp), dimension(lbwmax,n), intent(out) :: ssw
    real(kind=dp), intent(in) :: tol
    integer(kind=int32), dimension(n), intent(out) :: numrotsw
    integer(kind=int32), dimension(n), intent(out) :: lbws
    type(error_info), intent(inout), optional :: error
    integer(kind=int32), intent(in) :: n, lbwmax
    !
    integer(kind=int32), dimension(n,lbwmax) :: ksv
    real(kind=dp), dimension(n,lbwmax) :: csv
    complex(kind=dp), dimension(n,lbwmax) :: ssv
    type(routine_info), parameter :: info=info_f_c_general_wb

    if (failure(error)) return
    call push_id(info, error)
    
    call ip_transpose(a)
    call f_c_general_bv(a,n,lbws,lbwmax,numrotsw,ksv,csv,ssv,tol,error)

    if (success(error)) then
       
       call ip_transpose(a)

       jsw=transpose(ksv)
       csw=transpose(csv)
       ssw=transpose(ssv)
    end if
    call pop_id(error)
  end subroutine f_c_general_wb

  ! Errors:
  ! 0: no error
  ! 1: insufficient lower bw in wb%bc
  subroutine f_c_general_to_wb(a, n, b, lbw, ubw, lbwmax, ubwmax, &
       numrotsw, jsw, csw, ssw, tol, error)
    complex(kind=dp), target, dimension(n,n), intent(inout) :: a
    integer(kind=int32), dimension(lbwmax,n), intent(out) :: jsw
    real(kind=dp), dimension(lbwmax,n), intent(out) :: csw
    complex(kind=dp), dimension(lbwmax,n), intent(out) :: ssw
    complex(kind=dp), dimension(ubwmax+lbwmax+1,n), intent(out) :: b
    real(kind=dp), intent(in) :: tol
    integer(kind=int32), dimension(n), intent(out) :: numrotsw
    integer(kind=int32), intent(out) :: lbw
    type(error_info), intent(inout), optional :: error
    integer(kind=int32), intent(in) :: n, ubw, ubwmax, lbwmax
    !
    integer(kind=int32), dimension(n) :: lbws
    type(routine_info), parameter :: info=info_f_c_general_to_wb

    if (failure(error)) return
    call push_id(info, error)

    if (n == 1) then
       numrotsw=0;
       ssw=(0.0_dp, 0.0_dp); csw=0.0_dp; jsw=0
       lbw=0
       b(1,1)=a(1,1)
       return
    end if

    call f_c_general_wb(a, n, lbws, lbwmax, numrotsw, jsw, csw, ssw, tol, error)

    if (success(error)) then
       
       lbw=maxval(lbws)
       if (maxval(lbws) > lbwmax) then
          ! This should already have been detected in f_c_general_wb.
          call set_error(1, info, error); return
       else
          call c_extract_diagonals_bc(a, n, b, lbw, ubw, lbwmax, ubwmax)
       end if
    end if
    call pop_id(error)
  end subroutine f_c_general_to_wb

  ! Errors:
  ! 0: no error
  ! 1: n < 1
  ! 2: wb%ubwmax < ubw
  ! 3: n is not the same for a and wb.
  subroutine c_general_to_wb(a,wb,ubw,tol,error)
    complex(kind=dp), target, dimension(:,:), intent(inout) :: a
    type(c_wb), intent(inout) :: wb
    type(error_info), intent(inout), optional :: error
    real(kind=dp), intent(in) :: tol
    integer(kind=int32), intent(in) :: ubw
    type(routine_info), parameter :: info=info_c_general_to_wb

    if (failure(error)) return
    call push_id(info, error)
    
    if (size(a,1) < 1) then
       call set_error(1, info, error); return
    end if
    if (get_ubwmax(wb) < ubw) then
       call set_error(2, info, error); return
    end if
    if (get_n(wb) /= size(a,1) .or. get_n(wb) /= size(a,2)) then
       call set_error(3, info, error); return
    end if
    call f_c_general_to_wb(a,get_n(wb),wb%bc, wb%lbw, ubw, get_lbwmax(wb), get_ubwmax(wb), &
         wb%numrotsw, wb%jsw, wb%csw, wb%ssw, tol, error)
    wb%ubw=ubw
    call pop_id(error)
  end subroutine c_general_to_wb


end module mod_general_wb
