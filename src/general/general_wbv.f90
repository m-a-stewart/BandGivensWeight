module general_wbv
  use prec
  use error_id
  use utility
  use nested_types
  use band_types
  use general_bv
  implicit none

  private

  public :: f_general_wbv, f_d_general_wbv, f_c_general_wbv, &
       f_general_to_wbv, f_d_general_to_wbv, f_c_general_to_wbv, &
       general_to_wbv, d_general_to_wbv, c_general_to_wbv
  
  public :: info_d_general_to_wbv, info_f_d_general_to_wbv, info_c_general_to_wbv, &
       info_f_c_general_to_wbv

  interface f_general_wbv
     module procedure f_d_general_wbv, f_c_general_wbv
  end interface f_general_wbv

  interface f_general_to_wbv
     module procedure f_d_general_to_wbv, f_c_general_to_wbv
  end interface f_general_to_wbv

  interface general_to_wbv
     module procedure d_general_to_wbv, c_general_to_wbv
  end interface general_to_wbv

  type(routine_info), parameter :: info_d_general_to_wbv=routine_info(id_d_general_to_wbv, &
       'd_general_to_wbv', &
       [ character(len=error_message_length) :: 'n<1', &
       'Size of a and wbv not the same.' ] )

  type(routine_info), parameter :: info_f_d_general_to_wbv=routine_info(id_f_d_general_to_wbv, &
       'f_d_general_to_wbv', &
       [ character(len=error_message_length) :: 'Insufficient lower bandwidth in wbv', &
       'Insufficient upper bandwidth in wbv' ])

  type(routine_info), parameter :: info_c_general_to_wbv=routine_info(id_c_general_to_wbv, &
       'c_general_to_wbv', &
       [ character(len=error_message_length) :: 'n<1', &
       'Size of a and wbv not the same.' ] )

  type(routine_info), parameter :: info_f_c_general_to_wbv=routine_info(id_f_c_general_to_wbv, &
       'f_c_general_to_wbv', &
       [ character(len=error_message_length) :: 'Insufficient lower bandwidth in wbv', &
       'Insufficient upper bandwidth in wbv' ])

contains

  subroutine f_d_general_wbv(a, n, lbws, ubws, lbwmax, ubwmax, &
       numrotsw, jsw, csw, ssw, &
       numrotsv, ksv, csv, ssv, tol, error)
    real(kind=dp), target, dimension(n,n), intent(inout) :: a
    integer(kind=int32), dimension(lbwmax,n), intent(out) :: jsw
    real(kind=dp), dimension(lbwmax,n), intent(out) :: csw, ssw
    integer(kind=int32), dimension(n,ubwmax), intent(out) :: ksv
    real(kind=dp), dimension(n,ubwmax), intent(out) :: csv, ssv
    real(kind=dp), intent(in) :: tol
    integer(kind=int32), dimension(n), intent(out) :: numrotsw, numrotsv
    integer(kind=int32), dimension(n), intent(out) :: lbws, ubws
    type(error_info), intent(out) :: error
    integer(kind=int32), intent(in) :: n, lbwmax, ubwmax
    !
    integer(kind=int32), dimension(n,lbwmax) :: ksv0
    real(kind=dp), dimension(n,lbwmax) :: csv0, ssv0
    

    call f_d_general_bv(a,n,ubws,ubwmax,numrotsv,ksv,csv,ssv,tol,error)
    call ip_transpose(a)
    call f_d_general_bv(a,n,lbws,lbwmax,numrotsw,ksv0,csv0,ssv0,tol,error)
    call ip_transpose(a)
    
    jsw=transpose(ksv0)
    csw=transpose(csv0)
    ssw=transpose(ssv0)

  end subroutine f_d_general_wbv

  subroutine f_c_general_wbv(a, n, lbws, ubws, lbwmax, ubwmax, &
       numrotsw, jsw, csw, ssw, &
       numrotsv, ksv, csv, ssv, tol, error)
    complex(kind=dp), target, dimension(n,n), intent(inout) :: a
    integer(kind=int32), dimension(lbwmax,n), intent(out) :: jsw
    complex(kind=dp), dimension(lbwmax,n), intent(out) :: csw, ssw
    integer(kind=int32), dimension(n,ubwmax), intent(out) :: ksv
    complex(kind=dp), dimension(n,ubwmax), intent(out) :: csv, ssv
    real(kind=dp), intent(in) :: tol
    integer(kind=int32), dimension(n), intent(out) :: numrotsw, numrotsv
    integer(kind=int32), dimension(n), intent(out) :: lbws, ubws
    type(error_info), intent(out) :: error
    integer(kind=int32), intent(in) :: n, lbwmax, ubwmax
    !
    integer(kind=int32), dimension(n,lbwmax) :: ksv0
    complex(kind=dp), dimension(n,lbwmax) :: csv0, ssv0

    call f_c_general_bv(a,n,ubws,ubwmax,numrotsv,ksv,csv,ssv,tol,error)
    call ip_transpose(a)
    call f_c_general_bv(a,n,lbws,lbwmax,numrotsw,ksv0,csv0,ssv0,tol,error)
    call ip_transpose(a)
    
    jsw=transpose(ksv0)
    csw=transpose(csv0)
    ssw=transpose(ssv0)

  end subroutine f_c_general_wbv

  ! Errors:
  ! 0: no error
  ! 1: insufficient lower bw in wbv.
  ! 2: insufficient upper bw in wbv.
  subroutine f_d_general_to_wbv(a, n, br, lbw, ubw, lbwmax, ubwmax, &
       numrotsw, jsw, csw, ssw, &
       numrotsv, ksv, csv, ssv, tol, error)
    real(kind=dp), target, dimension(n,n), intent(inout) :: a
    real(kind=dp), dimension(n,ubwmax+lbwmax+1), intent(out) :: br
    integer(kind=int32), dimension(lbwmax,n), intent(out) :: jsw
    real(kind=dp), dimension(lbwmax,n), intent(out) :: csw, ssw
    integer(kind=int32), dimension(n,ubwmax), intent(out) :: ksv
    real(kind=dp), dimension(n,ubwmax), intent(out) :: csv, ssv
    real(kind=dp), intent(in) :: tol
    integer(kind=int32), dimension(n), intent(out) :: numrotsw, numrotsv
    integer(kind=int32), intent(out) :: lbw, ubw
    type(error_info), intent(out) :: error
    integer(kind=int32), intent(in) :: n, ubwmax, lbwmax
    !
    integer(kind=int32), dimension(n) :: lbws, ubws

    call clear_error(error)

    if (n == 1) then
       numrotsw=0;
       ssw=0.0_dp; csw=0.0_dp; jsw=0
       lbw=0
       numrotsv=0;
       ssv=0.0_dp; csv=0.0_dp; ksv=0
       ubw=0
       br(1,1)=a(1,1)
       return
    end if

    call f_d_general_wbv(a, n, lbws, ubws, lbwmax, ubwmax, &
         numrotsw, jsw, csw, ssw, &
         numrotsv, ksv, csv, ssv, tol, error)

    lbw=maxval(lbws)
    ubw=maxval(ubws)
    if (lbw > lbwmax) then
       ! This should already have been detected in f_d_general_wbv.
       call set_error(error, 1, id_f_d_general_to_wbv)
    else if (ubw > ubwmax) then
       call set_error(error, 2, id_f_d_general_to_wbv)
    else
       call d_extract_diagonals_br(a, n, br, lbw, ubw, lbwmax, ubwmax)
    end if
  end subroutine f_d_general_to_wbv

  ! Errors:
  ! 0: no error
  ! 1: insufficient lower bw in wbv.
  ! 2: insufficient upper bw in wbv.
  subroutine f_c_general_to_wbv(a, n, br, lbw, ubw, lbwmax, ubwmax, &
       numrotsw, jsw, csw, ssw, &
       numrotsv, ksv, csv, ssv, tol, error)
    complex(kind=dp), target, dimension(n,n), intent(inout) :: a
    complex(kind=dp), dimension(n,ubwmax+lbwmax+1), intent(out) :: br
    integer(kind=int32), dimension(lbwmax,n), intent(out) :: jsw
    complex(kind=dp), dimension(lbwmax,n), intent(out) :: csw, ssw
    integer(kind=int32), dimension(n,ubwmax), intent(out) :: ksv
    complex(kind=dp), dimension(n,ubwmax), intent(out) :: csv, ssv
    real(kind=dp), intent(in) :: tol
    integer(kind=int32), dimension(n), intent(out) :: numrotsw, numrotsv
    integer(kind=int32), intent(out) :: lbw, ubw
    type(error_info), intent(out) :: error
    integer(kind=int32), intent(in) :: n, ubwmax, lbwmax
    !
    integer(kind=int32), dimension(n) :: lbws, ubws

    call clear_error(error)

    if (n == 1) then
       numrotsw=0;
       ssw=0.0_dp; csw=0.0_dp; jsw=0
       lbw=0
       numrotsv=0;
       ssv=0.0_dp; csv=0.0_dp; ksv=0
       ubw=0
       br(1,1)=a(1,1)
       return
    end if

    call f_c_general_wbv(a, n, lbws, ubws, lbwmax, ubwmax, &
         numrotsw, jsw, csw, ssw, &
         numrotsv, ksv, csv, ssv, tol, error)

    lbw=maxval(lbws)
    ubw=maxval(ubws)
    if (lbw > lbwmax) then
       ! This should already have been detected in f_c_general_wbv.
       call set_error(error, 1, id_f_c_general_to_wbv)
    else if (ubw > ubwmax) then
       call set_error(error, 2, id_f_c_general_to_wbv)
    else
       call c_extract_diagonals_br(a, n, br, lbw, ubw, lbwmax, ubwmax)
    end if
  end subroutine f_c_general_to_wbv

  ! Errors:
  ! 1: n < 1
  ! 2: n is not the same for a and wbv.
  subroutine d_general_to_wbv(a,wbv,tol,error)
    real(kind=dp), target, dimension(:,:), intent(inout) :: a
    type(d_wbv), intent(inout) :: wbv
    type(error_info), intent(out) :: error
    real(kind=dp), intent(in) :: tol
    call clear_error(error)
    if (size(a,1) < 1) then
       call set_error(error, 1, id_d_general_to_wbv); return
    end if
    if (get_n(wbv) /= size(a,1) .or. get_n(wbv) /= size(a,2)) then
       call set_error(error, 2, id_d_general_to_wbv); return
    end if
    call f_d_general_to_wbv(a,get_n(wbv),wbv%br, wbv%lbw, wbv%ubw, get_lbwmax(wbv), &
         get_ubwmax(wbv), wbv%numrotsw, wbv%jsw, wbv%csw, wbv%ssw, &
         wbv%numrotsv, wbv%ksv, wbv%csv, wbv%ssv, tol, error)
  end subroutine d_general_to_wbv

  ! Errors:
  ! 1: n < 1
  ! 2: n is not the same for a and wbv.
  subroutine c_general_to_wbv(a,wbv,tol,error)
    complex(kind=dp), target, dimension(:,:), intent(inout) :: a
    type(c_wbv), intent(inout) :: wbv
    type(error_info), intent(out) :: error
    real(kind=dp), intent(in) :: tol
    call clear_error(error)
    if (size(a,1) < 1) then
       call set_error(error, 1, id_c_general_to_wbv); return
    end if
    if (get_n(wbv) /= size(a,1) .or. get_n(wbv) /= size(a,2)) then
       call set_error(error, 2, id_c_general_to_wbv); return
    end if
    call f_c_general_to_wbv(a,get_n(wbv),wbv%br, wbv%lbw, wbv%ubw, get_lbwmax(wbv), &
         get_ubwmax(wbv), wbv%numrotsw, wbv%jsw, wbv%csw, wbv%ssw, &
         wbv%numrotsv, wbv%ksv, wbv%csv, wbv%ssv, tol, error)
  end subroutine c_general_to_wbv

end module general_wbv
