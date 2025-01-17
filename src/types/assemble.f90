module mod_assemble
  use mod_prec
  use mod_error_id
  use mod_band_types
  use mod_orth_band_types
  use mod_rotation
  implicit none
  ! This module contains routines for assembling a Givens weight
  ! parameterized matrix into an unstructured representation in
  ! general $n\times n$ array.

  private

  public :: ub_to_general, d_ub_to_general, z_ub_to_general, &
       d_general_of_ub, z_general_of_ub, &
       f_ub_to_general, f_d_ub_to_general, f_z_ub_to_general, &
       bt_to_general, d_bt_to_general, z_bt_to_general, &
       d_general_of_bt, z_general_of_bt, &
       f_bt_to_general, f_d_bt_to_general, f_z_bt_to_general, &
       ubt_to_general, d_ubt_to_general, z_ubt_to_general, &
       d_general_of_ubt, z_general_of_ubt, &
       f_ubt_to_general, f_d_ubt_to_general, f_z_ubt_to_general, &
       bv_to_general, d_bv_to_general, z_bv_to_general, &
       d_general_of_bv, z_general_of_bv, &
       f_bv_to_general, f_d_bv_to_general, f_z_bv_to_general, &
       wb_to_general, d_wb_to_general, z_wb_to_general, &
       d_general_of_wb, z_general_of_wb, &       
       f_wb_to_general, f_d_wb_to_general, f_z_wb_to_general, &
       wbv_to_general, d_wbv_to_general, z_wbv_to_general, &
       d_general_of_wbv, z_general_of_wbv, &
       f_wbv_to_general, f_d_wbv_to_general, f_z_wbv_to_general, &
       general, to_general

  interface to_general
     module procedure d_ub_to_general, z_ub_to_general, d_bt_to_general, z_bt_to_general, &
          d_ubt_to_general, z_ubt_to_general, d_bv_to_general, z_bv_to_general, &
          d_wb_to_general, z_wb_to_general, d_wbv_to_general, z_wbv_to_general
  end interface to_general

  interface general
     module procedure d_general_of_ub, z_general_of_ub, d_general_of_bt, z_general_of_bt, &
          d_general_of_ubt, z_general_of_ubt, d_general_of_bv, z_general_of_bv, &
          d_general_of_wb, z_general_of_wb, d_general_of_wbv, z_general_of_wbv
  end interface general

  interface ub_to_general
     module procedure d_ub_to_general, z_ub_to_general
  end interface ub_to_general

  interface f_ub_to_general
     module procedure f_d_ub_to_general, f_z_ub_to_general
  end interface f_ub_to_general

  interface bt_to_general
     module procedure d_bt_to_general, z_bt_to_general
  end interface bt_to_general

  interface f_bt_to_general
     module procedure f_d_bt_to_general, f_z_bt_to_general
  end interface f_bt_to_general

  interface ubt_to_general
     module procedure d_ubt_to_general, z_ubt_to_general
  end interface ubt_to_general

  interface f_ubt_to_general
     module procedure f_d_ubt_to_general, f_z_ubt_to_general
  end interface f_ubt_to_general

  interface bv_to_general
     module procedure d_bv_to_general, z_bv_to_general
  end interface bv_to_general

  interface f_bv_to_general
     module procedure f_d_bv_to_general, f_z_bv_to_general
  end interface f_bv_to_general

  interface wb_to_general
     module procedure d_wb_to_general, z_wb_to_general
  end interface wb_to_general

  interface f_wb_to_general
     module procedure f_d_wb_to_general, f_z_wb_to_general
  end interface f_wb_to_general

  interface wbv_to_general
     module procedure d_wbv_to_general, z_wbv_to_general
  end interface wbv_to_general

  interface f_wbv_to_general
     module procedure f_d_wbv_to_general, f_z_wbv_to_general
  end interface f_wbv_to_general

contains

  function d_general_of_ub(ub, error) result(a)
    real(kind=dp), dimension(:,:), allocatable :: a
    type(d_ub), intent(in) :: ub
    type(error_info), intent(inout), optional :: error
    type(routine_info), parameter :: info=info_d_general_of_ub
    integer(kind=int32) :: n

    if (failure(error)) return
    call push_id(info, error)
    
    n=get_n(ub)
    allocate(a(n,n))
    call d_ub_to_general(ub,a,error)
    call pop_id(error)
  end function d_general_of_ub

  ! Errors
  ! 0: no error
  ! 1: ub%n /= n
  subroutine d_ub_to_general(ub,a,error)
    type(d_ub), intent(in) :: ub
    real(kind=dp), dimension(:,:), intent(out) :: a
    type(error_info), intent(inout), optional :: error
    type(routine_info), parameter :: info=info_d_ub_to_general

    if (failure(error)) return
    call push_id(info, error)

    if (get_n(ub) /= size(a,1) .or. get_n(ub) /= size(a,2)) then
       call set_error(1, info, error); return
    end if
    call f_d_ub_to_general(ub%bc, get_n(ub), ub%lbw, ub%ubw, get_lbwmax(ub), &
         get_ubwmax(ub), ub%numrotsu, ub%jsu, ub%csu, ub%ssu, a)
    call pop_id(error)
  end subroutine d_ub_to_general

  subroutine f_d_ub_to_general(bc, n, lbw, ubw, lbwmax, ubwmax, numrotsu, jsu, csu, ssu, a)
    real(kind=dp), target, dimension(n,n), intent(out) :: a
    integer(kind=int32), dimension(ubwmax,n), intent(in) :: jsu
    real(kind=dp), dimension(ubwmax,n), intent(in) :: csu, ssu
    real(kind=dp), dimension(lbwmax+ubwmax+1,n), intent(in) :: bc
    integer(kind=int32), dimension(n), intent(in) :: numrotsu
    integer(kind=int32), intent(in) :: ubw, lbw, n, lbwmax, ubwmax
    !
    integer(kind=int32) :: j,k
    call bc_to_general(bc,lbw,ubw,a)
    if (n==1) then
       return
    end if
    do k=n-1,2,-1
       do j=1,numrotsu(k)
          call f_rotation_times_general(csu(j,k),ssu(j,k),a(:,k+1:n),jsu(j,k),jsu(j,k)+1)
       end do
    end do
  end subroutine f_d_ub_to_general

  function z_general_of_ub(ub, error) result(a)
    complex(kind=dp), dimension(:,:), allocatable :: a
    type(z_ub), intent(in) :: ub
    type(error_info), intent(inout), optional :: error
    type(routine_info), parameter :: info=info_z_general_of_ub
    integer(kind=int32) :: n

    if (failure(error)) return
    call push_id(info, error)
    
    n=get_n(ub)
    allocate(a(n,n))
    call z_ub_to_general(ub,a,error)
    call pop_id(error)
  end function z_general_of_ub

  subroutine z_ub_to_general(ub,a,error)
    type(z_ub), intent(in) :: ub
    complex(kind=dp), dimension(:,:), intent(out) :: a
    type(error_info), intent(inout), optional :: error
    type(routine_info), parameter :: info=info_z_ub_to_general

    if (failure(error)) return
    call push_id(info, error)
    
    if (get_n(ub) /= size(a,1) .or. get_n(ub) /= size(a,2)) then
       call set_error(1, info, error); return
    end if
    call f_z_ub_to_general(ub%bc, get_n(ub), ub%lbw, ub%ubw, get_lbwmax(ub), get_ubwmax(ub), &
         ub%numrotsu, ub%jsu, ub%csu, ub%ssu, a)
    call pop_id(error)
  end subroutine z_ub_to_general

  subroutine f_z_ub_to_general(bc, n, lbw, ubw, lbwmax, ubwmax, numrotsu, jsu, csu, ssu, a)
    complex(kind=dp), target, dimension(n,n), intent(out) :: a
    integer(kind=int32), dimension(ubwmax,n), intent(in) :: jsu
    real(kind=dp), dimension(ubwmax,n), intent(in) :: csu
    complex(kind=dp), dimension(ubwmax,n), intent(in) :: ssu
    complex(kind=dp), dimension(lbwmax+ubwmax+1,n), intent(in) :: bc
    integer(kind=int32), dimension(n), intent(in) :: numrotsu
    integer(kind=int32), intent(in) :: ubw, lbw, n, lbwmax, ubwmax
    !
    integer(kind=int32) :: j,k

    call bc_to_general(bc,lbw,ubw,a)
    if (n==1) then
       return
    end if
    do k=n-1,2,-1
       do j=1,numrotsu(k)
          call f_rotation_times_general(csu(j,k),ssu(j,k),a(:,k+1:n),jsu(j,k),jsu(j,k)+1)
       end do
    end do
  end subroutine f_z_ub_to_general

  function d_general_of_bt(bt, error) result(a)
    real(kind=dp), dimension(:,:), allocatable :: a
    type(d_bt), intent(in) :: bt
    type(error_info), intent(inout), optional :: error
    type(routine_info), parameter :: info=info_d_general_of_bt
    integer(kind=int32) :: n

    if (failure(error)) return
    call push_id(info, error)
    
    n=get_n(bt)
    allocate(a(n,n))
    call d_bt_to_general(bt,a,error)
    call pop_id(error)
  end function d_general_of_bt

  ! Errors
  ! 0: no error
  ! 1: bt%n /= n
  subroutine d_bt_to_general(bt,a,error)
    type(d_bt), intent(in) :: bt
    real(kind=dp), dimension(:,:), intent(out) :: a
    type(error_info), intent(inout), optional :: error
    type(routine_info), parameter :: info=info_d_bt_to_general

    if (failure(error)) return
    call push_id(info, error)
    
    if (get_n(bt) /= size(a,1) .or. get_n(bt) /= size(a,2)) then
       call set_error(1, info, error); return
    end if
    call f_d_bt_to_general(bt%br, get_n(bt), bt%lbw, bt%ubw, get_lbwmax(bt), get_ubwmax(bt), &
         bt%numrotst, bt%kst, bt%cst, bt%sst, a)
    call pop_id(error)
  end subroutine d_bt_to_general

  subroutine f_d_bt_to_general(br, n, lbw, ubw, lbwmax, ubwmax, numrotst, kst, cst, sst, a)
    real(kind=dp), target, dimension(n,n), intent(out) :: a
    integer(kind=int32), dimension(n,lbwmax), intent(in) :: kst
    real(kind=dp), dimension(n,lbwmax), intent(in) :: cst, sst
    real(kind=dp), dimension(n,lbwmax+ubwmax+1), intent(in) :: br
    integer(kind=int32), dimension(n), intent(in) :: numrotst
    integer(kind=int32), intent(in) :: ubw, lbw, n, lbwmax, ubwmax
    !
    integer(kind=int32) :: j,k
    call br_to_general(br,lbw,ubw,a)
    if (n==1) then
       return
    end if
    do k=n-1,2,-1
       do j=1,numrotst(k)
          call f_general_times_rotation(a(k+1:n,:), cst(k,j), -sst(k,j), kst(k,j),kst(k,j)+1)
       end do
    end do
  end subroutine f_d_bt_to_general

  function z_general_of_bt(bt, error) result(a)
    complex(kind=dp), dimension(:,:), allocatable :: a
    type(z_bt), intent(in) :: bt
    type(error_info), intent(inout), optional :: error
    type(routine_info), parameter :: info=info_z_general_of_bt
    integer(kind=int32) :: n

    if (failure(error)) return
    call push_id(info, error)
    
    n=get_n(bt)
    allocate(a(n,n))
    call z_bt_to_general(bt,a,error)
    call pop_id(error)
  end function z_general_of_bt

  ! Errors
  ! 0: no error
  ! 1: bt%n /= n
  subroutine z_bt_to_general(bt,a,error)
    type(z_bt), intent(in) :: bt
    complex(kind=dp), dimension(:,:), intent(out) :: a
    type(error_info), intent(inout), optional :: error
    type(routine_info), parameter :: info=info_z_bt_to_general

    if (failure(error)) return
    call push_id(info, error)

    if (get_n(bt) /= size(a,1) .or. get_n(bt) /= size(a,2)) then
       call set_error(1, info, error); return
    end if
    call f_z_bt_to_general(bt%br, get_n(bt), bt%lbw, bt%ubw, get_lbwmax(bt), get_ubwmax(bt), &
         bt%numrotst, bt%kst, bt%cst, bt%sst, a)
    call pop_id(error)
  end subroutine z_bt_to_general

  subroutine f_z_bt_to_general(br, n, lbw, ubw, lbwmax, ubwmax, numrotst, kst, cst, sst, a)
    complex(kind=dp), target, dimension(n,n), intent(out) :: a
    integer(kind=int32), dimension(n,lbwmax), intent(in) :: kst
    real(kind=dp), dimension(n,lbwmax), intent(in) :: cst
    complex(kind=dp), dimension(n,lbwmax), intent(in) :: sst
    complex(kind=dp), dimension(n,lbwmax+ubwmax+1), intent(in) :: br
    integer(kind=int32), dimension(n), intent(in) :: numrotst
    integer(kind=int32), intent(in) :: ubw, lbw, n, lbwmax, ubwmax
    !
    integer(kind=int32) :: j,k
    call br_to_general(br,lbw,ubw,a)
    if (n==1) then
       return
    end if
    do k=n-1,2,-1
       do j=1,numrotst(k)
          call f_general_times_rotation(a(k+1:n,:), cst(k,j),-sst(k,j), kst(k,j),kst(k,j)+1)          
       end do
    end do
  end subroutine f_z_bt_to_general

  function d_general_of_ubt(ubt, error) result(a)
    real(kind=dp), dimension(:,:), allocatable :: a
    type(d_ubt), intent(in) :: ubt
    type(error_info), intent(inout), optional :: error
    type(routine_info), parameter :: info=info_d_general_of_ubt
    integer(kind=int32) :: n

    if (failure(error)) return

    call push_id(info, error)
    
    n=get_n(ubt)
    allocate(a(n,n))
    call d_ubt_to_general(ubt,a,error)
    call pop_id(error)
  end function d_general_of_ubt


  ! Errors
  ! 0: no error
  ! 1: ubt%n /= n
  subroutine d_ubt_to_general(ubt,a,error)
    type(d_ubt), intent(in) :: ubt
    real(kind=dp), dimension(:,:), intent(out) :: a
    type(error_info), intent(inout), optional :: error
    type(routine_info), parameter :: info=info_d_ubt_to_general

    if (failure(error)) return
    call push_id(info, error)

    if (get_n(ubt) /= size(a,1) .or. get_n(ubt) /= size(a,2)) then
       call set_error(1, info, error); return
    end if
    call f_d_ubt_to_general(ubt%bc, get_n(ubt), ubt%lbw, ubt%ubw, get_lbwmax(ubt), &
         get_ubwmax(ubt), ubt%numrotsu, ubt%jsu, ubt%csu, ubt%ssu, &
         ubt%numrotst, ubt%kst, ubt%cst, ubt%sst, a)
    call pop_id(error)
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

    call bc_to_general(bc,lbw,ubw,a)
    if (n==1) then
       return
    end if
    do k=n-1,2,-1
       do j=1,numrotsu(k)
          call f_rotation_times_general(csu(j,k),ssu(j,k),a(:,k+1:n),jsu(j,k),jsu(j,k)+1)
       end do
    end do

    do k=n-1,2,-1
       do j=1,numrotst(k)
          call f_general_times_rotation(a(k+1:n,:), cst(k,j),-sst(k,j), kst(k,j),kst(k,j)+1)
       end do
    end do
  end subroutine f_d_ubt_to_general

  function z_general_of_ubt(ubt, error) result(a)
    complex(kind=dp), dimension(:,:), allocatable :: a
    type(z_ubt), intent(in) :: ubt
    type(error_info), intent(inout), optional :: error
    type(routine_info), parameter :: info=info_z_general_of_ubt
    integer(kind=int32) :: n

    if (failure(error)) return
    call push_id(info, error)
    
    n=get_n(ubt)
    allocate(a(n,n))
    call z_ubt_to_general(ubt,a,error)
    call pop_id(error)
  end function z_general_of_ubt

  ! Errors
  ! 0: no error
  ! 1: ubt%n /= n
  subroutine z_ubt_to_general(ubt,a,error)
    type(z_ubt), intent(in) :: ubt
    complex(kind=dp), dimension(:,:), intent(out) :: a
    type(error_info), intent(inout), optional :: error
    type(routine_info), parameter :: info=info_z_ubt_to_general

    if (failure(error)) return
    call push_id(info, error)
    
    if (get_n(ubt) /= size(a,1) .or. get_n(ubt) /= size(a,2)) then
       call set_error(1, info, error); return
    end if
    call f_z_ubt_to_general(ubt%bc, get_n(ubt), ubt%lbw, ubt%ubw, get_lbwmax(ubt), &
         get_ubwmax(ubt), ubt%numrotsu, ubt%jsu, ubt%csu, ubt%ssu, &
         ubt%numrotst, ubt%kst, ubt%cst, ubt%sst, a)
    call pop_id(error)
  end subroutine z_ubt_to_general

  subroutine f_z_ubt_to_general(bc, n, lbw, ubw, lbwmax, ubwmax, numrotsu, jsu, csu, ssu, &
       numrotst, kst, cst, sst, a)
    complex(kind=dp), target, dimension(n,n), intent(out) :: a
    integer(kind=int32), dimension(ubwmax,n), intent(in) :: jsu
    real(kind=dp), dimension(ubwmax,n), intent(in) :: csu
    complex(kind=dp), dimension(ubwmax,n), intent(in) :: ssu
    integer(kind=int32), dimension(n,lbwmax), intent(in) :: kst
    real(kind=dp), dimension(n,lbwmax), intent(in) :: cst
    complex(kind=dp), dimension(n,lbwmax), intent(in) :: sst
    complex(kind=dp), dimension(lbwmax+ubwmax+1,n), intent(in) :: bc
    integer(kind=int32), dimension(n), intent(in) :: numrotsu, numrotst
    integer(kind=int32), intent(in) :: ubw, lbw, n, lbwmax, ubwmax
    !
    integer(kind=int32) :: j,k
    call bc_to_general(bc,lbw,ubw,a)
    if (n==1) then
       return
    end if
    do k=n-1,2,-1
       do j=1,numrotsu(k)
          call f_rotation_times_general(csu(j,k),ssu(j,k),a(:,k+1:n),jsu(j,k),jsu(j,k)+1)
       end do
    end do

    do k=n-1,2,-1
       do j=1,numrotst(k)
          call f_general_times_rotation(a(k+1:n,:), cst(k,j),-sst(k,j), kst(k,j),kst(k,j)+1)          
       end do
    end do
  end subroutine f_z_ubt_to_general

  function d_general_of_bv(bv, error) result(a)
    real(kind=dp), dimension(:,:), allocatable :: a
    type(d_bv), intent(in) :: bv
    type(error_info), intent(inout), optional :: error
    type(routine_info), parameter :: info=info_d_general_of_bv
    integer(kind=int32) :: n

    if (failure(error)) return
    call push_id(info, error)
    
    n=get_n(bv)
    allocate(a(n,n))
    call d_bv_to_general(bv,a,error)
    call pop_id(error)
  end function d_general_of_bv


  subroutine d_bv_to_general(bv,a,error)
    type(d_bv), intent(in) :: bv
    real(kind=dp), dimension(:,:), intent(out) :: a
    type(error_info), intent(inout), optional :: error
    type(routine_info), parameter :: info=info_d_bv_to_general
    
    if (failure(error)) return
    call push_id(info, error)
    if (get_n(bv) /= size(a,1) .or. get_n(bv) /= size(a,2)) then
       call set_error(1, info, error); return
    end if
    call f_d_bv_to_general(bv%br, get_n(bv), bv%lbw, bv%ubw, get_lbwmax(bv), get_ubwmax(bv), &
         bv%numrotsv, bv%ksv, bv%csv, bv%ssv, a)
    call pop_id(error)
  end subroutine d_bv_to_general

  subroutine f_d_bv_to_general(br, n, lbw, ubw, lbwmax, ubwmax, numrotsv, ksv, csv, ssv, a)
    real(kind=dp), target, dimension(n,n), intent(out) :: a
    integer(kind=int32), dimension(n,ubwmax), intent(in) :: ksv
    real(kind=dp), dimension(n, ubwmax), intent(in) :: csv, ssv
    real(kind=dp), dimension(n,lbwmax+ubwmax+1), intent(in) :: br
    integer(kind=int32), dimension(n), intent(in) :: numrotsv
    integer(kind=int32), intent(in) :: ubw, lbw, n, lbwmax, ubwmax
    !
    integer(kind=int32) :: j,k
    call br_to_general(br,lbw,ubw,a)
    if (n==1) then
       return
    end if
    do j=1,n-2
       do k=1,numrotsv(j)
          call f_general_times_rotation(a(1:j,:),csv(j,k),-ssv(j,k),ksv(j,k), ksv(j,k)+1)
       end do
    end do
  end subroutine f_d_bv_to_general

  function z_general_of_bv(bv, error) result(a)
    complex(kind=dp), dimension(:,:), allocatable :: a
    type(z_bv), intent(in) :: bv
    type(error_info), intent(inout), optional :: error
    type(routine_info), parameter :: info=info_z_general_of_bv
    integer(kind=int32) :: n

    if (failure(error)) return
    call push_id(info, error)
    
    n=get_n(bv)
    allocate(a(n,n))
    call z_bv_to_general(bv,a,error)
    call pop_id(error)
  end function z_general_of_bv

  subroutine z_bv_to_general(bv,a,error)
    type(z_bv), intent(in) :: bv
    complex(kind=dp), dimension(:,:), intent(out) :: a
    type(error_info), intent(inout), optional :: error
    type(routine_info), parameter :: info=info_z_bv_to_general

    if (failure(error)) return
    call push_id(info, error)
    if (get_n(bv) /= size(a,1) .or. get_n(bv) /= size(a,2)) then
       call set_error(1, info, error); return
    end if
    call f_z_bv_to_general(bv%br, get_n(bv), bv%lbw, bv%ubw, get_lbwmax(bv), get_ubwmax(bv), &
         bv%numrotsv, bv%ksv, bv%csv, bv%ssv, a)
    call pop_id(error)
  end subroutine z_bv_to_general

  subroutine f_z_bv_to_general(br, n, lbw, ubw, lbwmax, ubwmax, numrotsv, ksv, csv, ssv, a)
    complex(kind=dp), target, dimension(n,n), intent(out) :: a
    integer(kind=int32), dimension(n,ubwmax), intent(in) :: ksv
    real(kind=dp), dimension(n, ubwmax), intent(in) :: csv
    complex(kind=dp), dimension(n, ubwmax), intent(in) :: ssv
    complex(kind=dp), dimension(n,lbwmax+ubwmax+1), intent(in) :: br
    integer(kind=int32), dimension(n), intent(in) :: numrotsv
    integer(kind=int32), intent(in) :: ubw, lbw, n, lbwmax, ubwmax
    !
    integer(kind=int32) :: j,k
    call br_to_general(br,lbw,ubw,a)
    if (n==1) then
       return
    end if
    do j=1,n-2
       do k=1,numrotsv(j)
          call f_general_times_rotation(a(1:j,:),csv(j,k),-ssv(j,k),ksv(j,k), ksv(j,k)+1)
       end do
    end do
  end subroutine f_z_bv_to_general

  function d_general_of_wb(wb, error) result(a)
    real(kind=dp), dimension(:,:), allocatable :: a
    type(d_wb), intent(in) :: wb
    type(error_info), intent(inout), optional :: error
    type(routine_info), parameter :: info=info_d_general_of_wb
    integer(kind=int32) :: n

    if (failure(error)) return
    call push_id(info, error)
    
    n=get_n(wb)
    allocate(a(n,n))
    call d_wb_to_general(wb,a,error)
    call pop_id(error)
  end function d_general_of_wb

  subroutine d_wb_to_general(wb,a,error)
    type(d_wb), intent(in) :: wb
    real(kind=dp), dimension(:,:), intent(out) :: a
    type(error_info), intent(inout), optional :: error
    type(routine_info), parameter :: info=info_d_wb_to_general

    if (failure(error)) return
    call push_id(info, error)
    if (get_n(wb) /= size(a,1) .or. get_n(wb) /= size(a,2)) then
       call set_error(1, info, error); return
    end if
    call f_d_wb_to_general(wb%bc, get_n(wb), wb%lbw, wb%ubw, get_lbwmax(wb), get_ubwmax(wb), &
         wb%numrotsw, wb%jsw, wb%csw, wb%ssw, a)
    call pop_id(error)
  end subroutine d_wb_to_general

  subroutine f_d_wb_to_general(bc, n, lbw, ubw, lbwmax, ubwmax, numrotsw, jsw, csw, ssw, a)
    real(kind=dp), target, dimension(n,n), intent(out) :: a
    integer(kind=int32), dimension(lbwmax,n), intent(in) :: jsw
    real(kind=dp), dimension(lbwmax,n), intent(in) :: csw, ssw
    real(kind=dp), dimension(lbwmax+ubwmax+1,n), intent(in) :: bc
    integer(kind=int32), dimension(n), intent(in) :: numrotsw
    integer(kind=int32), intent(in) :: ubw, lbw, n, lbwmax, ubwmax
    !
    integer(kind=int32) :: j,k
    call bc_to_general(bc,lbw,ubw,a)
    if (n==1) then
       return
    end if

    do k=1,n-2
       do j=1,numrotsw(k)
          call f_rotation_times_general(csw(j,k),ssw(j,k),a(:,1:k),jsw(j,k),jsw(j,k)+1)
       end do
    end do
  end subroutine f_d_wb_to_general

  function z_general_of_wb(wb, error) result(a)
    complex(kind=dp), dimension(:,:), allocatable :: a
    type(z_wb), intent(in) :: wb
    type(error_info), intent(inout), optional :: error
    type(routine_info), parameter :: info=info_z_general_of_wb
    integer(kind=int32) :: n
    if (failure(error)) return
    call push_id(info, error)
    
    n=get_n(wb)
    allocate(a(n,n))
    call z_wb_to_general(wb,a,error)
    call pop_id(error)
  end function z_general_of_wb

  subroutine z_wb_to_general(wb,a,error)
    type(z_wb), intent(in) :: wb
    complex(kind=dp), dimension(:,:), intent(out) :: a
    type(error_info), intent(inout), optional :: error
    type(routine_info), parameter :: info=info_z_wb_to_general
    
    if (failure(error)) return
    call push_id(info, error)
    if (get_n(wb) /= size(a,1) .or. get_n(wb) /= size(a,2)) then
       call set_error(1, info, error); return
    end if
    call f_z_wb_to_general(wb%bc, get_n(wb), wb%lbw, wb%ubw, get_lbwmax(wb), get_ubwmax(wb), &
         wb%numrotsw, wb%jsw, wb%csw, wb%ssw, a)
    call pop_id(error)
  end subroutine z_wb_to_general

  subroutine f_z_wb_to_general(bc, n, lbw, ubw, lbwmax, ubwmax, numrotsw, jsw, csw, ssw, a)
    complex(kind=dp), target, dimension(n,n), intent(out) :: a
    integer(kind=int32), dimension(lbwmax,n), intent(in) :: jsw
    real(kind=dp), dimension(lbwmax,n), intent(in) :: csw
    complex(kind=dp), dimension(lbwmax,n), intent(in) :: ssw
    complex(kind=dp), dimension(lbwmax+ubwmax+1,n), intent(in) :: bc
    integer(kind=int32), dimension(n), intent(in) :: numrotsw
    integer(kind=int32), intent(in) :: ubw, lbw, n, lbwmax, ubwmax
    !
    integer(kind=int32) :: j,k
    call bc_to_general(bc,lbw,ubw,a)
    if (n==1) then
       return
    end if

    do k=1,n-2
       do j=1,numrotsw(k)
          call f_rotation_times_general(csw(j,k),ssw(j,k),a(:,1:k),jsw(j,k),jsw(j,k)+1)
       end do
    end do
  end subroutine f_z_wb_to_general

  function d_general_of_wbv(wbv, error) result(a)
    real(kind=dp), dimension(:,:), allocatable :: a
    type(d_wbv), intent(in) :: wbv
    type(error_info), intent(inout), optional :: error
    type(routine_info), parameter :: info=info_d_general_of_wbv
    integer(kind=int32) :: n
    if (failure(error)) return
    call push_id(info, error)
    
    n=get_n(wbv)
    allocate(a(n,n))
    call d_wbv_to_general(wbv,a,error)
    call pop_id(error)
  end function d_general_of_wbv

  subroutine d_wbv_to_general(wbv,a,error)
    type(d_wbv), intent(in) :: wbv
    real(kind=dp), dimension(:,:), intent(out) :: a
    type(error_info), intent(inout), optional :: error
    type(routine_info), parameter :: info=info_d_wbv_to_general
    
    if (failure(error)) return
    call push_id(info, error)
    if (get_n(wbv) /= size(a,1) .or. get_n(wbv) /= size(a,2)) then
       call set_error(1, info, error); return
    end if
    call f_d_wbv_to_general(wbv%br, get_n(wbv), wbv%lbw, wbv%ubw, get_lbwmax(wbv), &
         get_ubwmax(wbv), wbv%numrotsw, wbv%jsw, wbv%csw, wbv%ssw, &
         wbv%numrotsv, wbv%ksv, wbv%csv, wbv%ssv, a)
    call pop_id(error)
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
    call br_to_general(br,lbw,ubw,a)
    if (n==1) then
       return
    end if

    do j=1,n-2
       do k=1,numrotsv(j)
          call f_general_times_rotation(a(1:j,:),csv(j,k),-ssv(j,k),ksv(j,k), ksv(j,k)+1)
       end do
    end do

    do k=1,n-2
       do j=1,numrotsw(k)
          call f_rotation_times_general(csw(j,k),ssw(j,k),a(:,1:k),jsw(j,k),jsw(j,k)+1)
       end do
    end do
  end subroutine f_d_wbv_to_general

  function z_general_of_wbv(wbv, error) result(a)
    complex(kind=dp), dimension(:,:), allocatable :: a
    type(z_wbv), intent(in) :: wbv
    type(error_info), intent(inout), optional :: error
    type(routine_info), parameter :: info=info_z_general_of_wbv
    integer(kind=int32) :: n
    if (failure(error)) return
    call push_id(info, error)
    
    n=get_n(wbv)
    allocate(a(n,n))
    call z_wbv_to_general(wbv,a,error)
    call pop_id(error)
  end function z_general_of_wbv

  subroutine z_wbv_to_general(wbv,a,error)
    type(z_wbv), intent(in) :: wbv
    complex(kind=dp), dimension(:,:), intent(out) :: a
    type(error_info), intent(inout), optional :: error
    type(routine_info), parameter :: info=info_z_wbv_to_general
    
    if (failure(error)) return
    call push_id(info, error)
    
    if (get_n(wbv) /= size(a,1) .or. get_n(wbv) /= size(a,2)) then
       call set_error(1, info, error); return
    end if
    call f_z_wbv_to_general(wbv%br, get_n(wbv), wbv%lbw, wbv%ubw, get_lbwmax(wbv), &
         get_ubwmax(wbv), wbv%numrotsw, wbv%jsw, wbv%csw, wbv%ssw, &
         wbv%numrotsv, wbv%ksv, wbv%csv, wbv%ssv, a)
    call pop_id(error)
  end subroutine z_wbv_to_general

  subroutine f_z_wbv_to_general(br, n, lbw, ubw, lbwmax, ubwmax, numrotsw, jsw, csw, ssw, &
       numrotsv, ksv, csv, ssv, a)
    complex(kind=dp), target, dimension(n,n), intent(out) :: a
    integer(kind=int32), dimension(lbwmax,n), intent(in) :: jsw
    real(kind=dp), dimension(lbwmax,n), intent(in) :: csw
    complex(kind=dp), dimension(lbwmax,n), intent(in) :: ssw
    integer(kind=int32), dimension(n,ubwmax), intent(in) :: ksv
    real(kind=dp), dimension(n,ubwmax), intent(in) :: csv
    complex(kind=dp), dimension(n,ubwmax), intent(in) :: ssv
    complex(kind=dp), dimension(n,ubwmax+lbwmax+1), intent(in) :: br
    integer(kind=int32), dimension(n), intent(in) :: numrotsw, numrotsv
    integer(kind=int32), intent(in) :: ubw, lbw, n, lbwmax, ubwmax
    !
    integer(kind=int32) :: j,k
    call br_to_general(br,lbw,ubw,a)
    if (n==1) then
       return
    end if

    do j=1,n-2
       do k=1,numrotsv(j)
          call f_general_times_rotation(a(1:j,:),csv(j,k),-ssv(j,k),ksv(j,k), ksv(j,k)+1)
       end do
    end do

    do k=1,n-2
       do j=1,numrotsw(k)
          call f_rotation_times_general(csw(j,k),ssw(j,k),a(:,1:k),jsw(j,k),jsw(j,k)+1)
       end do
    end do
  end subroutine f_z_wbv_to_general

end module mod_assemble
