module mod_general_bt
  use mod_prec
  use mod_error_id
  use mod_utility
  use mod_orth_band_types
  use mod_band_types
  use mod_general_ub
  implicit none

  private

  public :: f_general_bt, f_d_general_bt, f_z_general_bt, &
       f_general_to_bt, f_d_general_to_bt, f_z_general_to_bt, &
       general_to_bt, d_general_to_bt, z_general_to_bt, &
       d_bt_of_general, z_bt_of_general, bt_of_general, bt

  interface bt_of_general
     module procedure d_bt_of_general, z_bt_of_general
  end interface bt_of_general

  interface bt
     module procedure d_bt_of_general, z_bt_of_general
  end interface bt
  
  interface f_general_bt
     module procedure f_d_general_bt, f_z_general_bt
  end interface f_general_bt

  interface f_general_to_bt
     module procedure f_d_general_to_bt, f_z_general_to_bt
  end interface f_general_to_bt

  interface general_to_bt
     module procedure d_general_to_bt, z_general_to_bt
  end interface general_to_bt

contains

  function d_bt_of_general(a, ubw, lbwmax, ubwmax, tol, error) result(bt)
    type(d_bt), allocatable :: bt
    real(kind=dp), target, dimension(:,:), intent(inout) :: a
    type(error_info), intent(inout), optional :: error
    real(kind=dp), intent(in) :: tol
    integer(kind=int32), intent(in) :: ubw, lbwmax, ubwmax
    type(routine_info), parameter :: info=info_d_bt_of_general
    integer(kind=int32) :: n

    if (failure(error)) return
    call push_id(info, error)

    n=size(a,1)
    bt=d_new_bt(n,lbwmax,ubwmax)

    call d_general_to_bt(a,bt,ubw,tol,error)

    call pop_id(error)
    
  end function d_bt_of_general
  
  subroutine f_d_general_bt(a, n, lbws, lbwmax, numrotst, kst, cst, sst, tol, error)
    real(kind=dp), target, dimension(n,n), intent(inout) :: a
    integer(kind=int32), dimension(n,lbwmax), intent(out) :: kst
    real(kind=dp), dimension(n,lbwmax), intent(out) :: cst, sst
    real(kind=dp), intent(in) :: tol
    integer(kind=int32), dimension(n), intent(out) :: numrotst
    integer(kind=int32), dimension(n), intent(out) :: lbws
    type(error_info), intent(inout), optional :: error
    integer(kind=int32), intent(in) :: n, lbwmax
    !
    integer(kind=int32), dimension(lbwmax,n) :: jsu
    real(kind=dp), dimension(lbwmax,n) :: csu, ssu
    type(routine_info), parameter :: info=info_f_d_general_bt

    if (failure(error)) return
    call push_id(info, error)
    
    call ip_transpose(a)
    call f_d_general_ub(a,n,lbws,lbwmax,numrotst,jsu,csu,ssu,tol,error)
    if (success(error)) then
       call ip_transpose(a)
       kst=transpose(jsu)
       cst=transpose(csu)
       sst=transpose(ssu)
    endif
    call pop_id(error)

  end subroutine f_d_general_bt

  ! Errors:
  ! 0: no error
  ! 1: insufficient lower bw in bt%bc
  subroutine f_d_general_to_bt(a, n, b, lbw, ubw, lbwmax, ubwmax, &
       numrotst, kst, cst, sst, tol, error)
    real(kind=dp), target, dimension(n,n), intent(inout) :: a
    integer(kind=int32), dimension(n,lbwmax), intent(out) :: kst
    real(kind=dp), dimension(n,lbwmax), intent(out) :: cst, sst
    real(kind=dp), dimension(n,ubwmax+lbwmax+1), intent(out) :: b
    real(kind=dp), intent(in) :: tol
    integer(kind=int32), dimension(n), intent(out) :: numrotst
    integer(kind=int32), intent(out) :: lbw
    type(error_info), intent(inout), optional :: error
    integer(kind=int32), intent(in) :: n, ubw, ubwmax, lbwmax
    !
    integer(kind=int32), dimension(n) :: lbws
    type(routine_info), parameter :: info=info_f_d_general_to_bt

    if (failure(error)) return
    call push_id(info, error)

    if (n == 1) then
       numrotst=0;
       sst=0.0_dp; cst=0.0_dp; kst=0
       lbw=0
       b(1,1)=a(1,1)
       return
    end if

    call f_d_general_bt(a, n, lbws, lbwmax, numrotst, kst, cst, sst, tol, error)

    if (success(error)) then
       lbw=maxval(lbws)
       if (maxval(lbws) > lbwmax) then
          ! This should already have been detected in f_d_general_bt.
          call set_error(1, info, error); return
       else
          call d_extract_diagonals_br(a, n, b, lbw, ubw, lbwmax, ubwmax)
       end if
    end if

    call pop_id(error)
  end subroutine f_d_general_to_bt

  ! Errors:
  ! 0: no error
  ! 1: n < 1
  ! 2: bt%lbwmax < lbw
  ! 3: n is not the same for a and ub or bv.
  subroutine d_general_to_bt(a,bt,ubw,tol,error)
    real(kind=dp), target, dimension(:,:), intent(inout) :: a
    type(d_bt), intent(inout) :: bt
    type(error_info), intent(inout), optional :: error
    real(kind=dp), intent(in) :: tol
    integer(kind=int32), intent(in) :: ubw
    type(routine_info), parameter :: info=info_d_general_to_bt

    if (failure(error)) return
    call push_id(info, error)

    if (size(a,1) < 1) then
       call set_error(1, info, error); return
    end if
    if (get_ubwmax(bt) < ubw) then
       call set_error(2, info, error); return
    end if
    if (get_n(bt) /= size(a,1) .or. get_n(bt) /= size(a,2)) then
       call set_error(3, info, error); return
    end if
    call f_d_general_to_bt(a,get_n(bt),bt%br, bt%lbw, ubw, get_lbwmax(bt), get_ubwmax(bt), &
         bt%numrotst, bt%kst, bt%cst, bt%sst, tol, error)
    bt%ubw=ubw
    call pop_id(error)
  end subroutine d_general_to_bt


  function z_bt_of_general(a, ubw, lbwmax, ubwmax, tol, error) result(bt)
    type(z_bt), allocatable :: bt
    complex(kind=dp), target, dimension(:,:), intent(inout) :: a
    type(error_info), intent(inout), optional :: error
    real(kind=dp), intent(in) :: tol
    integer(kind=int32), intent(in) :: ubw, lbwmax, ubwmax
    type(routine_info), parameter :: info=info_z_bt_of_general
    integer(kind=int32) :: n

    if (failure(error)) return
    call push_id(info, error)

    n=size(a,1)
    bt=z_new_bt(n,lbwmax,ubwmax)

    call z_general_to_bt(a,bt,ubw,tol,error)

    call pop_id(error)
    
  end function z_bt_of_general

  
  subroutine f_z_general_bt(a, n, lbws, lbwmax, numrotst, kst, cst, sst, tol, error)
    complex(kind=dp), target, dimension(n,n), intent(inout) :: a
    integer(kind=int32), dimension(n,lbwmax), intent(out) :: kst
    real(kind=dp), dimension(n,lbwmax), intent(out) :: cst
    complex(kind=dp), dimension(n,lbwmax), intent(out) :: sst
    real(kind=dp), intent(in) :: tol
    integer(kind=int32), dimension(n), intent(out) :: numrotst
    integer(kind=int32), dimension(n), intent(out) :: lbws
    type(error_info), intent(inout), optional :: error
    integer(kind=int32), intent(in) :: n, lbwmax
    !
    integer(kind=int32), dimension(lbwmax,n) :: jsu
    real(kind=dp), dimension(lbwmax,n) :: csu
    complex(kind=dp), dimension(lbwmax,n) :: ssu
    type(routine_info), parameter :: info=info_f_z_general_bt

    if (failure(error)) return
    call push_id(info, error)

    call ip_transpose(a)
    call f_z_general_ub(a,n,lbws,lbwmax,numrotst,jsu,csu,ssu,tol,error)

    if (success(error)) then
       
       call ip_transpose(a)

       kst=transpose(jsu)
       cst=transpose(csu)
       sst=transpose(ssu)
    end if

    call pop_id(error)
  end subroutine f_z_general_bt

  ! Errors:
  ! 0: no error
  ! 1: insufficient lower bw in bt%bc
  subroutine f_z_general_to_bt(a, n, b, lbw, ubw, lbwmax, ubwmax, &
       numrotst, kst, cst, sst, tol, error)
    complex(kind=dp), target, dimension(n,n), intent(inout) :: a
    integer(kind=int32), dimension(n,lbwmax), intent(out) :: kst
    real(kind=dp), dimension(n,lbwmax), intent(out) :: cst
    complex(kind=dp), dimension(n,lbwmax), intent(out) :: sst
    complex(kind=dp), dimension(n,ubwmax+lbwmax+1), intent(out) :: b
    real(kind=dp), intent(in) :: tol
    integer(kind=int32), dimension(n), intent(out) :: numrotst
    integer(kind=int32), intent(out) :: lbw
    type(error_info), intent(inout), optional :: error
    integer(kind=int32), intent(in) :: n, ubw, ubwmax, lbwmax
    !
    integer(kind=int32), dimension(n) :: lbws
    type(routine_info), parameter :: info=info_f_z_general_to_bt

    if (failure(error)) return
    call push_id(info, error)

    if (n == 1) then
       numrotst=0;
       sst=(0.0_dp,0.0_dp); cst=0.0_dp; kst=0
       lbw=0
       b(1,1)=a(1,1)
       return
    end if

    call f_z_general_bt(a, n, lbws, lbwmax, numrotst, kst, cst, sst, tol, error)

    if (success(error)) then
       lbw=maxval(lbws)
       if (maxval(lbws) > lbwmax) then
          ! This should already have been detected in f_z_general_bt.
          call set_error(1, info, error); return
       else
          call z_extract_diagonals_br(a, n, b, lbw, ubw, lbwmax, ubwmax)
       end if
    end if
    call pop_id(error)
  end subroutine f_z_general_to_bt

  ! Errors:
  ! 0: no error
  ! 1: n < 1
  ! 2: ub%ubwmax < ubw
  ! 3: n is not the same for a and ub or bv.
  subroutine z_general_to_bt(a,bt,ubw,tol,error)
    complex(kind=dp), target, dimension(:,:), intent(inout) :: a
    type(z_bt), intent(inout) :: bt
    type(error_info), intent(inout), optional :: error
    real(kind=dp), intent(in) :: tol
    integer(kind=int32), intent(in) :: ubw
    type(routine_info), parameter :: info=info_z_general_to_bt

    if (failure(error)) return
    call push_id(info, error)

    if (size(a,1) < 1) then
       call set_error(1, info, error); return
    end if
    if (get_ubwmax(bt) < ubw) then
       call set_error(2, info, error); return
    end if
    if (get_n(bt) /= size(a,1) .or. get_n(bt) /= size(a,2)) then
       call set_error(3, info, error); return
    end if
    
    call f_z_general_to_bt(a,get_n(bt),bt%br, bt%lbw, ubw, get_lbwmax(bt), get_ubwmax(bt), &
         bt%numrotst, bt%kst, bt%cst, bt%sst, tol, error)
    bt%ubw=ubw
    call pop_id(error)
    
  end subroutine z_general_to_bt

end module mod_general_bt
