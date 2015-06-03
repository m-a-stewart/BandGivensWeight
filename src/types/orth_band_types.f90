module mod_orth_band_types
  use mod_prec
  use mod_utility
  use mod_band_types
  implicit none
  ! This module contains derived types for UB, BV, WB, BT, UBT, and
  ! WBV decompositions.  The types include components to store a band
  ! matrix, information on the bandwidth and amount of available
  ! storage, sines and cosines for rotations, and information on the
  ! rows or columns on which the rotations act.  There are special
  ! allocation and deallocation routines for each type.
  private

  public :: d_ub, d_bt, d_ubt, z_ub, z_bt, z_ubt, &
       d_bv, d_wb, d_wbv, z_bv, z_wb, z_wbv

  public :: d_new_ub, d_new_bt, d_new_ubt, z_new_ub, z_new_bt, z_new_ubt, &
       d_new_bv, d_new_wb, d_new_wbv, z_new_bv, z_new_wb, z_new_wbv

  public :: deallocate_ub, d_deallocate_ub, z_deallocate_ub, &
       deallocate_bt, d_deallocate_bt, z_deallocate_bt, &
       deallocate_ubt, d_deallocate_ubt, z_deallocate_ubt, &
       deallocate_bv, d_deallocate_bv, z_deallocate_bv, &
       deallocate_wb, d_deallocate_wb, z_deallocate_wb, &
       deallocate_wbv, d_deallocate_wbv, z_deallocate_wbv

  public :: copy, d_copy_ub, z_copy_ub, &
       d_copy_bt, z_copy_bt, &
       d_copy_ubt, z_copy_ubt, &
       d_copy_bv, z_copy_bv, &
       d_copy_wb, z_copy_wb, &
       d_copy_wbv, z_copy_wbv

  public :: d_truncate_profile_ub, z_truncate_profile_ub, &
       d_truncate_profile_ubt, z_truncate_profile_ubt, &
       d_truncate_profile_bv, z_truncate_profile_bv, &
       d_truncate_profile_wbv, z_truncate_profile_wbv, &
       d_truncate_profile_wb, z_truncate_profile_wb, &
       d_truncate_profile_bt, z_truncate_profile_bt, &
       truncate_profile

  public :: get_n, d_ub_get_n, d_bv_get_n, z_ub_get_n, z_bv_get_n, &
       d_ubt_get_n, d_wbv_get_n, z_ubt_get_n, z_wbv_get_n, &
       d_bt_get_n, d_wb_get_n, z_bt_get_n, z_wb_get_n

  public :: get_lbwmax, d_ub_get_lbwmax, d_bv_get_lbwmax, z_ub_get_lbwmax, &
       z_bv_get_lbwmax, &
       d_ubt_get_lbwmax, d_wbv_get_lbwmax, z_ubt_get_lbwmax, z_wbv_get_lbwmax, &
       d_bt_get_lbwmax, d_wb_get_lbwmax, z_bt_get_lbwmax, z_wb_get_lbwmax

  public :: get_ubwmax, d_ub_get_ubwmax, d_bv_get_ubwmax, z_ub_get_ubwmax, &
       z_bv_get_ubwmax, &
       d_ubt_get_ubwmax, d_wbv_get_ubwmax, z_ubt_get_ubwmax, z_wbv_get_ubwmax, &
       d_bt_get_ubwmax, d_wb_get_ubwmax, z_bt_get_ubwmax, z_wb_get_ubwmax

  ! Stored by columns. (i.e. columns of A correspond to columns of B)
  type d_ub
     integer(kind=int32) :: lbw, ubw
     integer(kind=int32), private :: n, lbwmax, ubwmax
     real(kind=dp), dimension(:,:), allocatable :: bc
     integer(kind=int32), dimension(:), allocatable :: numrotsu
     integer(kind=int32), dimension(:,:), allocatable :: jsu
     real(kind=dp), dimension(:,:), allocatable :: csu, ssu
  end type d_ub

  ! Stored by rows.
  type d_bt
     integer(kind=int32) :: lbw, ubw
     integer(kind=int32), private :: n, lbwmax, ubwmax
     real(kind=dp), dimension(:,:), allocatable :: br
     integer(kind=int32), dimension(:), allocatable :: numrotst
     integer(kind=int32), dimension(:,:), allocatable :: kst
     real(kind=dp), dimension(:,:), allocatable :: cst, sst
  end type d_bt

  ! Stored by columns.
  type d_ubt
     integer(kind=int32) :: lbw, ubw
     integer(kind=int32), private :: n, lbwmax, ubwmax
     real(kind=dp), dimension(:,:), allocatable :: bc
     integer(kind=int32), dimension(:), allocatable :: numrotsu, numrotst
     integer(kind=int32), dimension(:,:), allocatable :: jsu, kst
     real(kind=dp), dimension(:,:), allocatable :: csu, ssu, cst, sst
  end type d_ubt

  ! Stored by columns.
  type z_ub
     integer(kind=int32) :: lbw, ubw
     integer(kind=int32), private :: n, lbwmax, ubwmax
     complex(kind=dp), dimension(:,:), allocatable :: bc
     integer(kind=int32), dimension(:), allocatable :: numrotsu
     integer(kind=int32), dimension(:,:), allocatable :: jsu
     complex(kind=dp), dimension(:,:), allocatable :: ssu
     real(kind=dp), dimension(:,:), allocatable :: csu
  end type z_ub

  ! Stored by rows.
  type z_bt
     integer(kind=int32) :: lbw, ubw
     integer(kind=int32), private :: n, lbwmax, ubwmax
     complex(kind=dp), dimension(:,:), allocatable :: br
     integer(kind=int32), dimension(:), allocatable :: numrotst
     integer(kind=int32), dimension(:,:), allocatable :: kst
     complex(kind=dp), dimension(:,:), allocatable :: sst
     real(kind=dp), dimension(:,:), allocatable :: cst
  end type z_bt

  ! Stored by columns.
  type z_ubt
     integer(kind=int32) :: lbw, ubw
     integer(kind=int32), private :: n, lbwmax, ubwmax
     complex(kind=dp), dimension(:,:), allocatable :: bc
     integer(kind=int32), dimension(:), allocatable :: numrotsu, numrotst
     integer(kind=int32), dimension(:,:), allocatable :: jsu, kst
     complex(kind=dp), dimension(:,:), allocatable :: ssu, sst
     real(kind=dp), dimension(:,:), allocatable :: csu, cst
  end type z_ubt

  ! Stored by rows.
  type d_bv
     integer(kind=int32) :: lbw, ubw
     integer(kind=int32), private :: n, lbwmax, ubwmax
     real(kind=dp), dimension(:,:), allocatable :: br
     integer(kind=int32), dimension(:), allocatable :: numrotsv
     integer(kind=int32), dimension(:,:), allocatable :: ksv
     real(kind=dp), dimension(:,:), allocatable :: csv, ssv
  end type d_bv

  ! Stored by columns.
  type d_wb
     integer(kind=int32) :: lbw, ubw
     integer(kind=int32), private :: n, lbwmax, ubwmax
     real(kind=dp), dimension(:,:), allocatable :: bc
     integer(kind=int32), dimension(:), allocatable :: numrotsw
     integer(kind=int32), dimension(:,:), allocatable :: jsw
     real(kind=dp), dimension(:,:), allocatable :: csw, ssw
  end type d_wb

  ! Stored by rows.
  type d_wbv
     integer(kind=int32) :: lbw, ubw
     integer(kind=int32), private :: n, lbwmax, ubwmax
     real(kind=dp), dimension(:,:), allocatable :: br
     integer(kind=int32), dimension(:), allocatable :: numrotsw, numrotsv
     integer(kind=int32), dimension(:,:), allocatable :: ksv, jsw
     real(kind=dp), dimension(:,:), allocatable :: csv, ssv, csw, ssw
  end type d_wbv

  ! Stored by rows.
  type z_bv
     integer(kind=int32) :: lbw, ubw
     integer(kind=int32), private :: n, lbwmax, ubwmax
     complex(kind=dp), dimension(:,:), allocatable :: br
     integer(kind=int32), dimension(:), allocatable :: numrotsv
     integer(kind=int32), dimension(:,:), allocatable :: ksv
     complex(kind=dp), dimension(:,:), allocatable :: ssv
     real(kind=dp), dimension(:,:), allocatable :: csv
  end type z_bv

  ! Stored by columns.
  type z_wb
     integer(kind=int32) :: lbw, ubw
     integer(kind=int32), private :: n, lbwmax, ubwmax
     complex(kind=dp), dimension(:,:), allocatable :: bc
     integer(kind=int32), dimension(:), allocatable :: numrotsw
     integer(kind=int32), dimension(:,:), allocatable :: jsw
     complex(kind=dp), dimension(:,:), allocatable :: ssw
     real(kind=dp), dimension(:,:), allocatable :: csw
  end type z_wb

  ! Stored by rows.
  type z_wbv
     integer(kind=int32) :: lbw, ubw
     integer(kind=int32), private :: n, lbwmax, ubwmax
     complex(kind=dp), dimension(:,:), allocatable :: br
     integer(kind=int32), dimension(:), allocatable :: numrotsv, numrotsw
     integer(kind=int32), dimension(:,:), allocatable :: ksv, jsw
     complex(kind=dp), dimension(:,:), allocatable :: ssv, ssw
     real(kind=dp), dimension(:,:), allocatable :: csv, csw
  end type z_wbv

  interface deallocate_ub
     module procedure d_deallocate_ub, z_deallocate_ub
  end interface deallocate_ub

  interface deallocate_bt
     module procedure d_deallocate_bt, z_deallocate_bt
  end interface deallocate_bt

  interface deallocate_ubt
     module procedure d_deallocate_ubt, z_deallocate_ubt
  end interface deallocate_ubt

  interface deallocate_bv
     module procedure d_deallocate_bv, z_deallocate_bv
  end interface deallocate_bv

  interface deallocate_wb
     module procedure d_deallocate_wb, z_deallocate_wb
  end interface deallocate_wb

  interface deallocate_wbv
     module procedure d_deallocate_wbv, z_deallocate_wbv
  end interface deallocate_wbv

  interface copy
     module procedure d_copy_ub, z_copy_ub, d_copy_bt, z_copy_bt,  &
          d_copy_ubt, z_copy_ubt, d_copy_bv, z_copy_bv, d_copy_wb, z_copy_wb, &
          d_copy_wbv, z_copy_wbv
  end interface copy
  
  interface truncate_profile
     module procedure d_truncate_profile_ub, z_truncate_profile_ub, &
          d_truncate_profile_ubt, z_truncate_profile_ubt, &
          d_truncate_profile_bv, z_truncate_profile_bv, &
          d_truncate_profile_wbv, z_truncate_profile_wbv, &
          d_truncate_profile_wb, z_truncate_profile_wb, &
          d_truncate_profile_bt, z_truncate_profile_bt
  end interface truncate_profile

  interface get_n
     module procedure d_ub_get_n, d_bv_get_n, z_ub_get_n, z_bv_get_n, &
          d_ubt_get_n, d_wbv_get_n, z_ubt_get_n, z_wbv_get_n, &
          d_bt_get_n, d_wb_get_n, z_bt_get_n, z_wb_get_n
  end interface get_n

  interface get_lbwmax
     module procedure d_ub_get_lbwmax, d_bv_get_lbwmax, z_ub_get_lbwmax, z_bv_get_lbwmax, &
          d_ubt_get_lbwmax, d_wbv_get_lbwmax, z_ubt_get_lbwmax, z_wbv_get_lbwmax, &
          d_bt_get_lbwmax, d_wb_get_lbwmax, z_bt_get_lbwmax, z_wb_get_lbwmax
  end interface get_lbwmax

  interface get_ubwmax
     module procedure d_ub_get_ubwmax, d_bv_get_ubwmax, z_ub_get_ubwmax, z_bv_get_ubwmax, &
          d_ubt_get_ubwmax, d_wbv_get_ubwmax, z_ubt_get_ubwmax, z_wbv_get_ubwmax, &
          d_bt_get_ubwmax, d_wb_get_ubwmax, z_bt_get_ubwmax, z_wb_get_ubwmax
  end interface get_ubwmax

contains

  ! UB

  type(d_ub) function d_new_ub(n,lbwmax,ubwmax) result(ub)
    integer(kind=int32), intent(in) :: n, lbwmax, ubwmax
    ub%n=n; ub%lbwmax=lbwmax; ub%ubwmax=ubwmax
    ub%ubw=0; ub%lbw=0
    allocate(ub%bc(lbwmax+ubwmax+1,n), ub%csu(ubwmax,n), ub%ssu(ubwmax,n), &
         ub%jsu(ubwmax,n), ub%numrotsu(n))
    ub%bc=0.0_dp; ub%csu=0.0_dp; ub%ssu=0.0_dp
    ub%jsu=0; ub%numrotsu=0
  end function d_new_ub

  subroutine d_deallocate_ub(ub)
    type(d_ub), intent(inout) :: ub
    call maybe_deallocate(ub%jsu)
    call maybe_deallocate(ub%numrotsu)
    call maybe_deallocate(ub%bc, ub%csu, ub%ssu)
    ub%n=0; ub%lbwmax=0; ub%ubwmax=0
  end subroutine d_deallocate_ub

  type(z_ub) function z_new_ub(n,lbwmax,ubwmax) result(ub)
    integer(kind=int32), intent(in) :: n, lbwmax, ubwmax
    ub%n=n; ub%lbwmax=lbwmax; ub%ubwmax=ubwmax
    ub%ubw=0; ub%lbw=0
    allocate(ub%bc(lbwmax+ubwmax+1,n), ub%csu(ubwmax,n), ub%ssu(ubwmax,n), &
         ub%jsu(ubwmax,n), ub%numrotsu(n))
    ub%bc=(0.0_dp,0.0_dp); ub%csu=0.0_dp; ub%ssu=(0.0_dp,0.0_dp)
    ub%jsu=0; ub%numrotsu=0
  end function z_new_ub

  subroutine z_deallocate_ub(ub)
    type(z_ub), intent(inout) :: ub
    call maybe_deallocate(ub%csu)
    call maybe_deallocate(ub%jsu)
    call maybe_deallocate(ub%numrotsu)
    call maybe_deallocate(ub%bc, ub%ssu)
    ub%n=0; ub%lbwmax=0; ub%ubwmax=0
    ub%lbw=0; ub%ubw=0
  end subroutine z_deallocate_ub

  ! BT

  type(d_bt) function d_new_bt(n,lbwmax,ubwmax) result(bt)
    integer(kind=int32), intent(in) :: n, lbwmax, ubwmax
    bt%n=n; bt%lbwmax=lbwmax; bt%ubwmax=ubwmax
    bt%ubw=0; bt%lbw=0
    allocate(bt%br(n,lbwmax+ubwmax+1), bt%cst(n,lbwmax), bt%sst(n,lbwmax), &
         bt%kst(n,lbwmax), bt%numrotst(n))
    bt%br=0.0_dp
    bt%cst=0.0_dp; bt%sst=0.0_dp
    bt%kst=0; bt%numrotst=0
  end function d_new_bt

  subroutine d_deallocate_bt(bt)
    type(d_bt), intent(inout) :: bt
    call maybe_deallocate(bt%kst)
    call maybe_deallocate(bt%numrotst)
    call maybe_deallocate(bt%br, bt%cst, bt%sst)
    bt%n=0; bt%lbwmax=0; bt%ubwmax=0;
    bt%lbw=0; bt%ubw=0
  end subroutine d_deallocate_bt

  type(z_bt) function z_new_bt(n,lbwmax,ubwmax) result(bt)
    integer(kind=int32), intent(in) :: n, lbwmax, ubwmax
    bt%n=n; bt%lbwmax=lbwmax; bt%ubwmax=ubwmax
    bt%ubw=0; bt%lbw=0
    allocate(bt%br(n,lbwmax+ubwmax+1), bt%cst(n,lbwmax), bt%sst(n,lbwmax), &
         bt%kst(n,lbwmax), bt%numrotst(n))
    bt%br=(0.0_dp,0.0_dp)
    bt%cst=0.0_dp; bt%sst=(0.0_dp,0.0_dp)
    bt%kst=0; bt%numrotst=0
  end function z_new_bt

  subroutine z_deallocate_bt(bt)
    type(z_bt), intent(inout) :: bt
    call maybe_deallocate(bt%cst)
    call maybe_deallocate(bt%kst)
    call maybe_deallocate(bt%numrotst)
    call maybe_deallocate(bt%br, bt%sst)
    bt%n=0; bt%lbwmax=0; bt%ubwmax=0;
    bt%lbw=0; bt%ubw=0
  end subroutine z_deallocate_bt


  ! UBT

  type(d_ubt) function d_new_ubt(n,lbwmax,ubwmax) result(ubt)
    integer(kind=int32), intent(in) :: n, lbwmax, ubwmax
    ubt%n=n; ubt%lbwmax=lbwmax; ubt%ubwmax=ubwmax
    ubt%ubw=0; ubt%lbw=0
    allocate(ubt%bc(lbwmax+ubwmax+1,n), ubt%csu(ubwmax,n), ubt%ssu(ubwmax,n), &
         ubt%jsu(ubwmax,n), ubt%numrotsu(n), ubt%cst(n,lbwmax), ubt%sst(n,lbwmax), &
         ubt%kst(n,lbwmax), ubt%numrotst(n))
    ubt%bc=0.0_dp
    ubt%csu=0.0_dp; ubt%ssu=0.0_dp
    ubt%jsu=0; ubt%numrotsu=0
    ubt%cst=0.0_dp; ubt%sst=0.0_dp
    ubt%kst=0; ubt%numrotst=0
  end function d_new_ubt

  subroutine d_deallocate_ubt(ubt)
    type(d_ubt), intent(inout) :: ubt
    call maybe_deallocate(ubt%jsu, ubt%kst)
    call maybe_deallocate(ubt%numrotsu, ubt%numrotst)
    call maybe_deallocate(ubt%bc, ubt%csu, ubt%ssu, ubt%cst, ubt%sst)
    ubt%n=0; ubt%lbwmax=0; ubt%ubwmax=0;
    ubt%lbw=0; ubt%ubw=0
  end subroutine d_deallocate_ubt

  type(z_ubt) function z_new_ubt(n,lbwmax,ubwmax) result(ubt)
    integer(kind=int32), intent(in) :: n, lbwmax, ubwmax
    ubt%n=n; ubt%lbwmax=lbwmax; ubt%ubwmax=ubwmax
    ubt%ubw=0; ubt%lbw=0
    allocate(ubt%bc(lbwmax+ubwmax+1,n), ubt%csu(ubwmax,n), ubt%ssu(ubwmax,n), &
         ubt%jsu(ubwmax,n), ubt%numrotsu(n), ubt%cst(n,lbwmax), ubt%sst(n,lbwmax), &
         ubt%kst(n,lbwmax), ubt%numrotst(n))
    ubt%bc=(0.0_dp,0.0_dp)
    ubt%csu=0.0_dp; ubt%ssu=(0.0_dp,0.0_dp)
    ubt%jsu=0; ubt%numrotsu=0
    ubt%cst=0.0_dp; ubt%sst=(0.0_dp,0.0_dp)
    ubt%kst=0; ubt%numrotst=0
  end function z_new_ubt

  subroutine z_deallocate_ubt(ubt)
    type(z_ubt), intent(inout) :: ubt
    call maybe_deallocate(ubt%csu, ubt%cst)
    call maybe_deallocate(ubt%jsu, ubt%kst)
    call maybe_deallocate(ubt%numrotsu, ubt%numrotst)
    call maybe_deallocate(ubt%bc, ubt%ssu, ubt%sst)
    ubt%n=0; ubt%lbwmax=0; ubt%ubwmax=0;
    ubt%lbw=0; ubt%ubw=0
  end subroutine z_deallocate_ubt

  ! BV

  type(d_bv) function d_new_bv(n,lbwmax,ubwmax) result(bv)
    integer(kind=int32), intent(in) :: n, lbwmax, ubwmax
    bv%n=n; bv%lbwmax=lbwmax; bv%ubwmax=ubwmax
    bv%ubw=0; bv%lbw=0
    allocate(bv%br(n,lbwmax+ubwmax+1), bv%csv(n,ubwmax), bv%ssv(n,ubwmax), &
         bv%ksv(n,ubwmax), bv%numrotsv(n))
    bv%br=0.0_dp; bv%csv=0.0_dp; bv%ssv=0.0_dp
    bv%ksv=0; bv%numrotsv=0
  end function d_new_bv

  subroutine d_deallocate_bv(bv)
    type(d_bv), intent(inout) :: bv
    call maybe_deallocate(bv%ksv)
    call maybe_deallocate(bv%numrotsv)
    call maybe_deallocate(bv%br, bv%csv, bv%ssv)
    bv%n=0; bv%lbwmax=0; bv%ubwmax=0
    bv%lbw=0; bv%ubw=0
  end subroutine d_deallocate_bv

  type(z_bv) function z_new_bv(n,lbwmax,ubwmax) result(bv)
    integer(kind=int32), intent(in) :: n, lbwmax, ubwmax
    bv%n=n; bv%lbwmax=lbwmax; bv%ubwmax=ubwmax
    bv%ubw=0; bv%lbw=0
    allocate(bv%br(n,lbwmax+ubwmax+1), bv%csv(n,ubwmax), bv%ssv(n,ubwmax), &
         bv%ksv(n,ubwmax), bv%numrotsv(n))
    bv%br=(0.0_dp, 0.0_dp); bv%csv=0.0_dp; bv%ssv=(0.0_dp, 0.0_dp)
    bv%ksv=0; bv%numrotsv=0
  end function z_new_bv

  subroutine z_deallocate_bv(bv)
    type(z_bv), intent(inout) :: bv
    call maybe_deallocate(bv%csv)
    call maybe_deallocate(bv%ksv)
    call maybe_deallocate(bv%numrotsv)
    call maybe_deallocate(bv%br, bv%ssv)
    bv%n=0; bv%lbwmax=0; bv%ubwmax=0
    bv%lbw=0; bv%ubw=0
  end subroutine z_deallocate_bv

  ! WB

  type(d_wb) function d_new_wb(n,lbwmax,ubwmax) result(wb)
    integer(kind=int32), intent(in) :: n, lbwmax, ubwmax
    wb%n=n; wb%lbwmax=lbwmax; wb%ubwmax=ubwmax
    wb%ubw=0; wb%lbw=0
    allocate(wb%bc(lbwmax+ubwmax+1,n), wb%csw(lbwmax,n), wb%ssw(lbwmax,n), wb%jsw(lbwmax,n), &
         wb%numrotsw(n))
    wb%bc=0.0_dp; wb%csw=0.0_dp; wb%ssw=0.0_dp
    wb%jsw=0; wb%numrotsw=0
  end function d_new_wb

  subroutine d_deallocate_wb(wb)
    type(d_wb), intent(inout) :: wb
    call maybe_deallocate(wb%jsw)
    call maybe_deallocate(wb%numrotsw)
    call maybe_deallocate(wb%bc, wb%csw, wb%ssw)
    wb%n=0; wb%lbwmax=0; wb%ubwmax=0
    wb%lbw=0; wb%ubw=0
  end subroutine d_deallocate_wb

  type(z_wb) function z_new_wb(n,lbwmax,ubwmax) result(wb)
    integer(kind=int32), intent(in) :: n, lbwmax, ubwmax
    wb%n=n; wb%lbwmax=lbwmax; wb%ubwmax=ubwmax
    wb%ubw=0; wb%lbw=0
    allocate(wb%bc(lbwmax+ubwmax+1,n), wb%csw(lbwmax,n), &
         wb%ssw(lbwmax,n), wb%jsw(lbwmax,n), wb%numrotsw(n))
    wb%bc=(0.0_dp,0.0_dp); wb%csw=0.0_dp; wb%ssw=(0.0_dp,0.0_dp)
    wb%jsw=0; wb%numrotsw=0
  end function z_new_wb

  subroutine z_deallocate_wb(wb)
    type(z_wb), intent(inout) :: wb
    call maybe_deallocate(wb%csw)
    call maybe_deallocate(wb%jsw)
    call maybe_deallocate(wb%numrotsw)
    call maybe_deallocate(wb%bc, wb%ssw)
    wb%n=0; wb%lbwmax=0; wb%ubwmax=0
    wb%lbw=0; wb%ubw=0
  end subroutine z_deallocate_wb

  ! WBV

  type(d_wbv) function d_new_wbv(n,lbwmax,ubwmax) result(wbv)
    integer(kind=int32), intent(in) :: n, lbwmax, ubwmax
    wbv%n=n; wbv%lbwmax=lbwmax; wbv%ubwmax=ubwmax
    wbv%ubw=0; wbv%lbw=0
    allocate(wbv%br(n,lbwmax+ubwmax+1), wbv%csv(n,ubwmax), wbv%ssv(n,ubwmax), &
         wbv%ksv(n,ubwmax), wbv%numrotsv(n), &
         wbv%csw(lbwmax,n), wbv%ssw(lbwmax,n), wbv%jsw(lbwmax,n), wbv%numrotsw(n))
    wbv%br=0.0_dp; wbv%csv=0.0_dp; wbv%ssv=0.0_dp
    wbv%ksv=0; wbv%numrotsv=0
    wbv%csw=0.0_dp; wbv%ssw=0.0_dp
    wbv%jsw=0; wbv%numrotsw=0
  end function d_new_wbv

  subroutine d_deallocate_wbv(wbv)
    type(d_wbv), intent(inout) :: wbv
    call maybe_deallocate(wbv%ksv, wbv%jsw)
    call maybe_deallocate(wbv%numrotsv, wbv%numrotsw)
    call maybe_deallocate(wbv%br, wbv%csv, wbv%ssv, wbv%csw, wbv%ssw)
    wbv%n=0; wbv%lbwmax=0; wbv%ubwmax=0
    wbv%lbw=0; wbv%ubw=0
  end subroutine d_deallocate_wbv

  type(z_wbv) function z_new_wbv(n,lbwmax,ubwmax) result(wbv)
    integer(kind=int32), intent(in) :: n, lbwmax, ubwmax
    wbv%n=n; wbv%lbwmax=lbwmax; wbv%ubwmax=ubwmax
    wbv%ubw=0; wbv%lbw=0
    allocate(wbv%br(n,lbwmax+ubwmax+1), wbv%csv(n,ubwmax), wbv%ssv(n,ubwmax), &
         wbv%ksv(n,ubwmax), wbv%numrotsv(n), &
         wbv%csw(lbwmax,n), wbv%ssw(lbwmax,n), wbv%jsw(lbwmax,n), wbv%numrotsw(n))
    wbv%br=(0.0_dp,0.0_dp); wbv%csv=0.0_dp; wbv%ssv=(0.0_dp,0.0_dp)
    wbv%ksv=0; wbv%numrotsv=0
    wbv%csw=0.0_dp; wbv%ssw=(0.0_dp,0.0_dp)
    wbv%jsw=0; wbv%numrotsw=0
  end function z_new_wbv

  subroutine z_deallocate_wbv(wbv)
    type(z_wbv), intent(inout) :: wbv
    call maybe_deallocate(wbv%csv, wbv%csw)
    call maybe_deallocate(wbv%ksv, wbv%jsw)
    call maybe_deallocate(wbv%numrotsv, wbv%numrotsw)
    call maybe_deallocate(wbv%br, wbv%ssv, wbv%ssw)
    wbv%n=0; wbv%lbwmax=0; wbv%ubwmax=0
    wbv%lbw=0; wbv%ubw=0
  end subroutine z_deallocate_wbv

  ! copy

  subroutine d_copy_ub(ub2,ub1)
    type(d_ub), intent(in) :: ub1
    type(d_ub), intent(inout) :: ub2
    ub2%n=ub1%n
    ub2%lbw=ub1%lbw
    ub2%ubw=ub1%ubw
    if ( ub2%ubwmax < ub1%ubw .or. ub2%lbwmax < ub1%lbw .or. &
         size(ub2%bc,2) /= ub1%n .or. not_all_allocated(ub2%bc, ub2%csu, ub2%ssu) .or. &
         not_all_allocated(ub2%jsu) .or. not_all_allocated(ub2%numrotsu)) then
       call d_deallocate_ub(ub2)
       allocate(ub2%bc(ub1%lbwmax+ub1%ubwmax+1,ub1%n), ub2%csu(ub1%ubwmax,ub1%n), &
            ub2%ssu(ub1%ubwmax,ub1%n), ub2%jsu(ub1%ubwmax,ub1%n), ub2%numrotsu(ub1%n))
       ub2%lbwmax=ub1%lbwmax
       ub2%ubwmax=ub1%ubwmax
    end if
    ub2%bc(1:ub1%lbw+ub1%ubw+1,:)=ub1%bc(1:ub1%lbw+ub1%ubw+1,:)
    ub2%numrotsu=ub1%numrotsu
    ub2%csu(1:ub1%ubw,:)=ub1%csu(1:ub1%ubw,:)
    ub2%ssu(1:ub1%ubw,:)=ub1%ssu(1:ub1%ubw,:)
    ub2%jsu(1:ub1%ubw,:)=ub1%jsu(1:ub1%ubw,:)
  end subroutine d_copy_ub

  subroutine z_copy_ub(ub2,ub1)
    type(z_ub), intent(in) :: ub1
    type(z_ub), intent(inout) :: ub2
    ub2%n=ub1%n
    ub2%lbw=ub1%lbw
    ub2%ubw=ub1%ubw
    if ( ub2%ubwmax < ub1%ubw .or. ub2%lbwmax < ub1%lbw .or. &
         size(ub2%bc,2) /= ub1%n .or. not_all_allocated(ub2%bc, ub2%ssu) .or. &
         not_all_allocated(ub2%csu) .or. not_all_allocated(ub2%jsu) .or. &
         not_all_allocated(ub2%numrotsu)) then
       call z_deallocate_ub(ub2)
       allocate(ub2%bc(ub1%lbwmax+ub1%ubwmax+1,ub1%n), ub2%csu(ub1%ubwmax,ub1%n), &
            ub2%ssu(ub1%ubwmax,ub1%n), ub2%jsu(ub1%ubwmax,ub1%n), ub2%numrotsu(ub1%n))
       ub2%lbwmax=ub1%lbwmax
       ub2%ubwmax=ub1%ubwmax
    end if
    ub2%bc(1:ub1%lbw+ub1%ubw+1,:)=ub1%bc(1:ub1%lbw+ub1%ubw+1,:)
    ub2%numrotsu=ub1%numrotsu
    ub2%csu(1:ub1%ubw,:)=ub1%csu(1:ub1%ubw,:)
    ub2%ssu(1:ub1%ubw,:)=ub1%ssu(1:ub1%ubw,:)
    ub2%jsu(1:ub1%ubw,:)=ub1%jsu(1:ub1%ubw,:)
  end subroutine z_copy_ub

  subroutine d_copy_bt(bt2,bt1)
    type(d_bt), intent(in) :: bt1
    type(d_bt), intent(inout) :: bt2
    bt2%n=bt1%n
    bt2%lbw=bt1%lbw
    bt2%ubw=bt1%ubw
    if ( bt2%ubwmax < bt1%ubw .or. bt2%lbwmax < bt1%lbw .or. &
         size(bt2%br,1) /= bt1%n .or. not_all_allocated(bt2%br, bt2%cst, bt2%sst) .or. &
         not_all_allocated(bt2%kst) .or. not_all_allocated(bt2%numrotst) ) then
       call d_deallocate_bt(bt2)
       allocate(bt2%br(bt1%n,bt1%lbwmax+bt1%ubwmax+1), &
            bt2%cst(bt1%n, bt1%lbwmax), bt2%sst(bt1%n,bt1%lbwmax), &
            bt2%kst(bt1%n, bt1%lbwmax), bt2%numrotst(bt1%n))
       bt2%lbwmax=bt1%lbwmax
       bt2%ubwmax=bt1%ubwmax
    end if
    
    bt2%br(:,1:bt1%lbw+bt1%ubw+1)=bt1%br(:,1:bt1%lbw+bt1%ubw+1)
    bt2%numrotst=bt1%numrotst
    bt2%cst(:,1:bt1%lbw)=bt1%cst(:,1:bt1%lbw)
    bt2%sst(:,1:bt1%lbw)=bt1%sst(:,1:bt1%lbw)
    bt2%kst(:,1:bt1%lbw)=bt1%kst(:,1:bt1%lbw)
  end subroutine d_copy_bt

  subroutine z_copy_bt(bt2,bt1)
    type(z_bt), intent(in) :: bt1
    type(z_bt), intent(inout) :: bt2
    bt2%n=bt1%n
    bt2%lbw=bt1%lbw
    bt2%ubw=bt1%ubw
    if ( bt2%ubwmax < bt1%ubw .or. bt2%lbwmax < bt1%lbw .or. &
         size(bt2%br,1) /= bt1%n .or. not_all_allocated(bt2%br, bt2%sst) .or. &
         not_all_allocated(bt2%cst) .or. not_all_allocated(bt2%kst) .or. &
         not_all_allocated(bt2%numrotst) ) then
       call z_deallocate_bt(bt2)
       allocate(bt2%br(bt1%n,bt1%lbwmax+bt1%ubwmax+1), &
            bt2%cst(bt1%n, bt1%lbwmax), bt2%sst(bt1%n,bt1%lbwmax), &
            bt2%kst(bt1%n, bt1%lbwmax), bt2%numrotst(bt1%n))
       bt2%lbwmax=bt1%lbwmax
       bt2%ubwmax=bt1%ubwmax
    end if
    bt2%br(:,1:bt1%lbw+bt1%ubw+1)=bt1%br(:,1:bt1%lbw+bt1%ubw+1)
    bt2%numrotst=bt1%numrotst
    bt2%cst(:,1:bt1%lbw)=bt1%cst(:,1:bt1%lbw)
    bt2%sst(:,1:bt1%lbw)=bt1%sst(:,1:bt1%lbw)
    bt2%kst(:,1:bt1%lbw)=bt1%kst(:,1:bt1%lbw)
  end subroutine z_copy_bt

  subroutine d_copy_ubt(ubt2,ubt1)
    type(d_ubt), intent(in) :: ubt1
    type(d_ubt), intent(inout) :: ubt2
    ubt2%n=ubt1%n
    ubt2%lbw=ubt1%lbw
    ubt2%ubw=ubt1%ubw
    if ( ubt2%ubwmax < ubt1%ubw .or. ubt2%lbwmax < ubt1%lbw .or. &
         size(ubt2%bc,2) /= ubt1%n .or. &
         not_all_allocated(ubt2%bc, ubt2%csu, ubt2%ssu, ubt2%cst, ubt2%sst) .or. &
         not_all_allocated(ubt2%jsu, ubt2%kst) .or. &
         not_all_allocated(ubt2%numrotsu, ubt2%numrotst) ) then
       call d_deallocate_ubt(ubt2)
       allocate(ubt2%bc(ubt1%lbwmax+ubt1%ubwmax+1,ubt1%n), ubt2%csu(ubt1%ubwmax,ubt1%n), &
            ubt2%ssu(ubt1%ubwmax,ubt1%n), ubt2%jsu(ubt1%ubwmax,ubt1%n), ubt2%numrotsu(ubt1%n), &
            ubt2%cst(ubt1%n, ubt1%lbwmax), ubt2%sst(ubt1%n, ubt1%lbwmax), &
            ubt2%kst(ubt1%n, ubt1%lbwmax), ubt2%numrotst(ubt1%n))
       ubt2%lbwmax=ubt1%lbwmax
       ubt2%ubwmax=ubt1%ubwmax
    end if
    
    ubt2%bc(1:ubt1%lbw+ubt1%ubw+1,:)=ubt1%bc(1:ubt1%lbw+ubt1%ubw+1,:)
    ubt2%numrotsu=ubt1%numrotsu
    ubt2%csu(1:ubt1%ubw,:)=ubt1%csu(1:ubt1%ubw,:)
    ubt2%ssu(1:ubt1%ubw,:)=ubt1%ssu(1:ubt1%ubw,:)
    ubt2%jsu(1:ubt1%ubw,:)=ubt1%jsu(1:ubt1%ubw,:)
    
    ubt2%numrotst=ubt1%numrotst
    ubt2%cst(:,1:ubt1%lbw)=ubt1%cst(:,1:ubt1%lbw)
    ubt2%sst(:,1:ubt1%lbw)=ubt1%sst(:,1:ubt1%lbw)
    ubt2%kst(:,1:ubt1%lbw)=ubt1%kst(:,1:ubt1%lbw)
  end subroutine d_copy_ubt

  subroutine z_copy_ubt(ubt2,ubt1)
    type(z_ubt), intent(in) :: ubt1
    type(z_ubt), intent(inout) :: ubt2
    ubt2%n=ubt1%n
    ubt2%lbw=ubt1%lbw
    ubt2%ubw=ubt1%ubw
    if ( ubt2%ubwmax < ubt1%ubw .or. ubt2%lbwmax < ubt1%lbw .or. &
         size(ubt2%bc,2) /= ubt1%n .or. &
         not_all_allocated(ubt2%bc, ubt2%ssu, ubt2%sst) .or. &
         not_all_allocated(ubt2%csu,ubt2%cst) .or. &
         not_all_allocated(ubt2%jsu, ubt2%kst) .or. &
         not_all_allocated(ubt2%numrotsu, ubt2%numrotst)) then
       call z_deallocate_ubt(ubt2)
       allocate(ubt2%bc(ubt1%lbwmax+ubt1%ubwmax+1,ubt1%n), ubt2%csu(ubt1%ubwmax,ubt1%n), &
            ubt2%ssu(ubt1%ubwmax,ubt1%n), ubt2%jsu(ubt1%ubwmax,ubt1%n), ubt2%numrotsu(ubt1%n), &
            ubt2%cst(ubt1%n, ubt1%lbwmax), ubt2%sst(ubt1%n, ubt1%lbwmax), &
            ubt2%kst(ubt1%n, ubt1%lbwmax), ubt2%numrotst(ubt1%n))
       ubt2%lbwmax=ubt1%lbwmax
       ubt2%ubwmax=ubt1%ubwmax
    end if
    
    ubt2%bc(1:ubt1%lbw+ubt1%ubw+1,:)=ubt1%bc(1:ubt1%lbw+ubt1%ubw+1,:)
    ubt2%numrotsu=ubt1%numrotsu
    ubt2%csu(1:ubt1%ubw,:)=ubt1%csu(1:ubt1%ubw,:)
    ubt2%ssu(1:ubt1%ubw,:)=ubt1%ssu(1:ubt1%ubw,:)
    ubt2%jsu(1:ubt1%ubw,:)=ubt1%jsu(1:ubt1%ubw,:)
    
    ubt2%numrotst=ubt1%numrotst
    ubt2%cst(:,1:ubt1%lbw)=ubt1%cst(:,1:ubt1%lbw)
    ubt2%sst(:,1:ubt1%lbw)=ubt1%sst(:,1:ubt1%lbw)
    ubt2%kst(:,1:ubt1%lbw)=ubt1%kst(:,1:ubt1%lbw)
  end subroutine z_copy_ubt

  subroutine d_copy_bv(bv2,bv1)
    type(d_bv), intent(in) :: bv1
    type(d_bv), intent(inout) :: bv2
    bv2%n=bv1%n
    bv2%lbw=bv1%lbw
    bv2%ubw=bv1%ubw
    if ( bv2%ubwmax < bv1%ubw .or. bv2%lbwmax < bv1%lbw .or. &
         size(bv2%br,1) /= bv1%n .or. not_all_allocated(bv2%br, bv2%csv, bv2%ssv) .or. &
         not_all_allocated(bv2%ksv) .or. not_all_allocated(bv2%numrotsv) ) then
       call d_deallocate_bv(bv2)
       allocate(bv2%br(bv1%n,bv1%lbwmax+bv1%ubwmax+1), bv2%csv(bv1%n,bv1%ubwmax), &
            bv2%ssv(bv1%n,bv1%ubwmax), bv2%ksv(bv1%n,bv1%ubwmax), bv2%numrotsv(bv1%n))
       bv2%lbwmax=bv1%lbwmax
       bv2%ubwmax=bv1%ubwmax
    end if
    bv2%br(:,1:bv1%lbw+bv1%ubw+1)=bv1%br(:,1:bv1%lbw+bv1%ubw+1)
    bv2%numrotsv=bv1%numrotsv
    bv2%csv(:,1:bv1%ubw)=bv1%csv(:,1:bv1%ubw)
    bv2%ssv(:,1:bv1%ubw)=bv1%ssv(:,1:bv1%ubw)
    bv2%ksv(:,1:bv1%ubw)=bv1%ksv(:,1:bv1%ubw)
  end subroutine d_copy_bv

  subroutine d_copy_wb(wb2,wb1)
    type(d_wb), intent(in) :: wb1
    type(d_wb), intent(inout) :: wb2
    wb2%n=wb1%n
    wb2%lbw=wb1%lbw
    wb2%ubw=wb1%ubw
    if ( wb2%ubwmax < wb1%ubw .or. wb2%lbwmax < wb1%lbw .or. &
         size(wb2%bc,2) /= wb1%n .or. not_all_allocated(wb2%bc, wb2%csw, wb2%ssw) .or. &
         not_all_allocated(wb2%jsw) .or. not_all_allocated(wb2%numrotsw) ) then
       call d_deallocate_wb(wb2)
       allocate(wb2%bc(wb1%lbwmax+wb1%ubwmax+1,wb1%n), &
            wb2%csw(wb1%lbwmax,wb1%n), wb2%ssw(wb1%lbwmax,wb1%n), &
            wb2%jsw(wb1%lbwmax,wb1%n), wb2%numrotsw(wb1%n))
       wb2%lbwmax=wb1%lbwmax
       wb2%ubwmax=wb1%ubwmax
    end if
    wb2%bc(1:wb1%lbw+wb1%ubw+1,:)=wb1%bc(1:wb1%lbw+wb1%ubw+1,:)
    wb2%numrotsw=wb1%numrotsw
    wb2%csw(1:wb1%lbw,:)=wb1%csw(1:wb1%lbw,:)
    wb2%ssw(1:wb1%lbw,:)=wb1%ssw(1:wb1%lbw,:)
    wb2%jsw(1:wb1%lbw,:)=wb1%jsw(1:wb1%lbw,:)
  end subroutine d_copy_wb

  subroutine d_copy_wbv(wbv2,wbv1)
    type(d_wbv), intent(in) :: wbv1
    type(d_wbv), intent(inout) :: wbv2
    wbv2%n=wbv1%n
    wbv2%lbw=wbv1%lbw
    wbv2%ubw=wbv1%ubw
    if ( wbv2%ubwmax < wbv1%ubw .or. wbv2%lbwmax < wbv1%lbw .or. &
         size(wbv2%br,1) /= wbv1%n .or. &
         not_all_allocated(wbv2%br, wbv2%csw, wbv2%ssw, wbv2%csv, wbv2%ssv) .or. &
         not_all_allocated(wbv2%jsw, wbv2%ksv) .or. &
         not_all_allocated(wbv2%numrotsw, wbv2%numrotsv) ) then
       call d_deallocate_wbv(wbv2)
       allocate(wbv2%br(wbv1%n,wbv1%lbwmax+wbv1%ubwmax+1), wbv2%csv(wbv1%n,wbv1%ubwmax), &
            wbv2%ssv(wbv1%n,wbv1%ubwmax), wbv2%ksv(wbv1%n,wbv1%ubwmax), wbv2%numrotsv(wbv1%n), &
            wbv2%csw(wbv1%lbwmax,wbv1%n), wbv2%ssw(wbv1%lbwmax,wbv1%n), &
            wbv2%jsw(wbv1%lbwmax,wbv1%n), wbv2%numrotsw(wbv1%n))
       wbv2%lbwmax=wbv1%lbwmax
       wbv2%ubwmax=wbv1%ubwmax
    end if
    wbv2%br(:,1:wbv1%lbw+wbv1%ubw+1)=wbv1%br(:,1:wbv1%lbw+wbv1%ubw+1)

    wbv2%numrotsv=wbv1%numrotsv
    wbv2%csv(:,1:wbv1%ubw)=wbv1%csv(:,1:wbv1%ubw)
    wbv2%ssv(:,1:wbv1%ubw)=wbv1%ssv(:,1:wbv1%ubw)
    wbv2%ksv(:,1:wbv1%ubw)=wbv1%ksv(:,1:wbv1%ubw)

    wbv2%numrotsw=wbv1%numrotsw
    wbv2%csw(1:wbv1%lbw,:)=wbv1%csw(1:wbv1%lbw,:)
    wbv2%ssw(1:wbv1%lbw,:)=wbv1%ssw(1:wbv1%lbw,:)
    wbv2%jsw(1:wbv1%lbw,:)=wbv1%jsw(1:wbv1%lbw,:)
  end subroutine d_copy_wbv

  subroutine z_copy_bv(bv2,bv1)
    type(z_bv), intent(in) :: bv1
    type(z_bv), intent(inout) :: bv2
    bv2%n=bv1%n
    bv2%lbw=bv1%lbw
    bv2%ubw=bv1%ubw

    if ( bv2%ubwmax < bv1%ubw .or. bv2%lbwmax < bv1%lbw .or. &
         size(bv2%br,1) /= bv1%n .or. not_all_allocated(bv2%br, bv2%ssv) .or. &
         not_all_allocated(bv2%csv) .or. &
         not_all_allocated(bv2%ksv) .or. not_all_allocated(bv2%numrotsv) ) then
       call z_deallocate_bv(bv2)
       allocate(bv2%br(bv1%n,bv1%lbwmax+bv1%ubwmax+1), bv2%csv(bv1%n,bv1%ubwmax), &
            bv2%ssv(bv1%n,bv1%ubwmax), bv2%ksv(bv1%n,bv1%ubwmax), bv2%numrotsv(bv1%n))
       bv2%lbwmax=bv1%lbwmax
       bv2%ubwmax=bv1%ubwmax
    end if
    bv2%br(:,1:bv1%lbw+bv1%ubw+1)=bv1%br(:,1:bv1%lbw+bv1%ubw+1)
    bv2%numrotsv=bv1%numrotsv
    bv2%csv(:,1:bv1%ubw)=bv1%csv(:,1:bv1%ubw)
    bv2%ssv(:,1:bv1%ubw)=bv1%ssv(:,1:bv1%ubw)
    bv2%ksv(:,1:bv1%ubw)=bv1%ksv(:,1:bv1%ubw)
  end subroutine z_copy_bv

  subroutine z_copy_wb(wb2,wb1)
    type(z_wb), intent(in) :: wb1
    type(z_wb), intent(inout) :: wb2
    wb2%n=wb1%n
    wb2%lbw=wb1%lbw
    wb2%ubw=wb1%ubw
    if ( wb2%ubwmax < wb1%ubw .or. wb2%lbwmax < wb1%lbw .or. &
         size(wb2%bc,2) /= wb1%n .or. not_all_allocated(wb2%bc, wb2%ssw) .or. &
         not_all_allocated(wb2%csw) .or. &
         not_all_allocated(wb2%jsw) .or. not_all_allocated(wb2%numrotsw) ) then
       call z_deallocate_wb(wb2)
       allocate(wb2%bc(wb1%lbwmax+wb1%ubwmax+1,wb1%n), &
            wb2%csw(wb1%lbwmax,wb1%n), wb2%ssw(wb1%lbwmax,wb1%n), &
            wb2%jsw(wb1%lbwmax,wb1%n), wb2%numrotsw(wb1%n))
       wb2%lbwmax=wb1%lbwmax
       wb2%ubwmax=wb1%ubwmax
    end if
    wb2%bc(1:wb1%lbw+wb1%ubw+1,:)=wb1%bc(1:wb1%lbw+wb1%ubw+1,:)
    wb2%numrotsw=wb1%numrotsw
    wb2%csw(1:wb1%lbw,:)=wb1%csw(1:wb1%lbw,:)
    wb2%ssw(1:wb1%lbw,:)=wb1%ssw(1:wb1%lbw,:)
    wb2%jsw(1:wb1%lbw,:)=wb1%jsw(1:wb1%lbw,:)
  end subroutine z_copy_wb

  subroutine z_copy_wbv(wbv2,wbv1)
    type(z_wbv), intent(in) :: wbv1
    type(z_wbv), intent(inout) :: wbv2
    wbv2%n=wbv1%n
    wbv2%lbw=wbv1%lbw
    wbv2%ubw=wbv1%ubw
    if ( wbv2%ubwmax < wbv1%ubw .or. wbv2%lbwmax < wbv1%lbw .or. &
         size(wbv2%br,1) /= wbv1%n .or. &
         not_all_allocated(wbv2%br, wbv2%ssw, wbv2%ssv) .or. &
         not_all_allocated(wbv2%jsw, wbv2%ksv) .or. &
         not_all_allocated(wbv2%csw, wbv2%csv) .or. &
         not_all_allocated(wbv2%numrotsw, wbv2%numrotsv) ) then
       call z_deallocate_wbv(wbv2)
       allocate(wbv2%br(wbv1%n,wbv1%lbwmax+wbv1%ubwmax+1), wbv2%csv(wbv1%n,wbv1%ubwmax), &
            wbv2%ssv(wbv1%n,wbv1%ubwmax), wbv2%ksv(wbv1%n,wbv1%ubwmax), wbv2%numrotsv(wbv1%n), &
            wbv2%csw(wbv1%lbwmax,wbv1%n), wbv2%ssw(wbv1%lbwmax,wbv1%n), &
            wbv2%jsw(wbv1%lbwmax,wbv1%n), wbv2%numrotsw(wbv1%n))
       wbv2%lbwmax=wbv1%lbwmax
       wbv2%ubwmax=wbv1%ubwmax
    end if
    wbv2%br(:,1:wbv1%lbw+wbv1%ubw+1)=wbv1%br(:,1:wbv1%lbw+wbv1%ubw+1)

    wbv2%numrotsv=wbv1%numrotsv
    wbv2%csv(:,1:wbv1%ubw)=wbv1%csv(:,1:wbv1%ubw)
    wbv2%ssv(:,1:wbv1%ubw)=wbv1%ssv(:,1:wbv1%ubw)
    wbv2%ksv(:,1:wbv1%ubw)=wbv1%ksv(:,1:wbv1%ubw)

    wbv2%numrotsw=wbv1%numrotsw
    wbv2%csw(1:wbv1%lbw,:)=wbv1%csw(1:wbv1%lbw,:)
    wbv2%ssw(1:wbv1%lbw,:)=wbv1%ssw(1:wbv1%lbw,:)
    wbv2%jsw(1:wbv1%lbw,:)=wbv1%jsw(1:wbv1%lbw,:)
  end subroutine z_copy_wbv

  subroutine d_truncate_profile_ub(ub,lower,upper)
    type(d_ub), intent(inout) :: ub
    integer(kind=int32), dimension(:), intent(in), optional :: lower, upper

    call d_truncate_profile_bc(ub%bc,get_n(ub), ub%ubw, get_lbwmax(ub), &
         get_ubwmax(ub), lower, upper)
  end subroutine d_truncate_profile_ub

  subroutine z_truncate_profile_ub(ub,lower,upper)
    type(z_ub), intent(inout) :: ub
    integer(kind=int32), dimension(:), intent(in), optional :: lower, upper

    call z_truncate_profile_bc(ub%bc,get_n(ub), ub%ubw, get_lbwmax(ub), &
         get_ubwmax(ub), lower, upper)
  end subroutine z_truncate_profile_ub

  subroutine d_truncate_profile_bt(bt,lower,upper)
    type(d_bt), intent(inout) :: bt
    integer(kind=int32), dimension(:), intent(in), optional :: lower, upper

    call d_truncate_profile_br(bt%br,get_n(bt), bt%lbw, get_lbwmax(bt), &
         get_ubwmax(bt), lower, upper)
  end subroutine d_truncate_profile_bt

  subroutine z_truncate_profile_bt(bt,lower,upper)
    type(z_bt), intent(inout) :: bt
    integer(kind=int32), dimension(:), intent(in), optional :: lower, upper

    call z_truncate_profile_br(bt%br,get_n(bt), bt%lbw, get_lbwmax(bt), &
         get_ubwmax(bt), lower, upper)
  end subroutine z_truncate_profile_bt

  subroutine d_truncate_profile_bv(bv,lower,upper)
    type(d_bv), intent(inout) :: bv
    integer(kind=int32), dimension(:), intent(in), optional :: lower, upper

    call d_truncate_profile_br(bv%br,get_n(bv), bv%lbw, get_lbwmax(bv), &
         get_ubwmax(bv), lower, upper)
  end subroutine d_truncate_profile_bv

  subroutine z_truncate_profile_bv(bv,lower,upper)
    type(z_bv), intent(inout) :: bv
    integer(kind=int32), dimension(:), intent(in), optional :: lower, upper

    call z_truncate_profile_br(bv%br,get_n(bv), bv%lbw, get_lbwmax(bv), &
         get_ubwmax(bv), lower, upper)
  end subroutine z_truncate_profile_bv

  subroutine d_truncate_profile_ubt(ubt,lower,upper)
    type(d_ubt), intent(inout) :: ubt
    integer(kind=int32), dimension(:), intent(in), optional :: lower, upper

    call d_truncate_profile_bc(ubt%bc,get_n(ubt), ubt%ubw, get_lbwmax(ubt), &
         get_ubwmax(ubt), lower, upper)
  end subroutine d_truncate_profile_ubt

  subroutine z_truncate_profile_ubt(ubt,lower,upper)
    type(z_ubt), intent(inout) :: ubt
    integer(kind=int32), dimension(:), intent(in), optional :: lower, upper

    call z_truncate_profile_bc(ubt%bc,get_n(ubt), ubt%ubw, get_lbwmax(ubt), &
         get_ubwmax(ubt), lower, upper)
  end subroutine z_truncate_profile_ubt
  
  subroutine d_truncate_profile_wbv(wbv,lower,upper)
    type(d_wbv), intent(inout) :: wbv
    integer(kind=int32), dimension(:), intent(in), optional :: lower, upper

    call d_truncate_profile_br(wbv%br,get_n(wbv), wbv%lbw, get_lbwmax(wbv), &
         get_ubwmax(wbv), lower, upper)
  end subroutine d_truncate_profile_wbv

  subroutine z_truncate_profile_wbv(wbv,lower,upper)
    type(z_wbv), intent(inout) :: wbv
    integer(kind=int32), dimension(:), intent(in), optional :: lower, upper

    call z_truncate_profile_br(wbv%br,get_n(wbv), wbv%lbw, get_lbwmax(wbv), &
         get_ubwmax(wbv), lower, upper)
  end subroutine z_truncate_profile_wbv

  subroutine d_truncate_profile_wb(wb,lower,upper)
    type(d_wb), intent(inout) :: wb
    integer(kind=int32), dimension(:), intent(in), optional :: lower, upper

    call d_truncate_profile_bc(wb%bc,get_n(wb), wb%ubw, get_lbwmax(wb), &
         get_ubwmax(wb), lower, upper)
  end subroutine d_truncate_profile_wb

  subroutine z_truncate_profile_wb(wb,lower,upper)
    type(z_wb), intent(inout) :: wb
    integer(kind=int32), dimension(:), intent(in), optional :: lower, upper

    call z_truncate_profile_bc(wb%bc,get_n(wb), wb%ubw, get_lbwmax(wb), &
         get_ubwmax(wb), lower, upper)
  end subroutine z_truncate_profile_wb
  
  integer(kind=int32) function d_ub_get_n(ub) result(n)
    type(d_ub) :: ub
    n=ub%n
  end function d_ub_get_n

  integer(kind=int32) function d_bt_get_n(bt) result(n)
    type(d_bt) :: bt
    n=bt%n
  end function d_bt_get_n

  integer(kind=int32) function d_ubt_get_n(ubt) result(n)
    type(d_ubt) :: ubt
    n=ubt%n
  end function d_ubt_get_n

  integer(kind=int32) function d_bv_get_n(bv) result(n)
    type(d_bv) :: bv
    n=bv%n
  end function d_bv_get_n

  integer(kind=int32) function d_wb_get_n(wb) result(n)
    type(d_wb) :: wb
    n=wb%n
  end function d_wb_get_n

  integer(kind=int32) function d_wbv_get_n(wbv) result(n)
    type(d_wbv) :: wbv
    n=wbv%n
  end function d_wbv_get_n

  integer(kind=int32) function z_ub_get_n(ub) result(n)
    type(z_ub) :: ub
    n=ub%n
  end function z_ub_get_n

  integer(kind=int32) function z_bv_get_n(bv) result(n)
    type(z_bv) :: bv
    n=bv%n
  end function z_bv_get_n

  integer(kind=int32) function z_bt_get_n(bt) result(n)
    type(z_bt) :: bt
    n=bt%n
  end function z_bt_get_n

  integer(kind=int32) function z_wb_get_n(wb) result(n)
    type(z_wb) :: wb
    n=wb%n
  end function z_wb_get_n

  integer(kind=int32) function z_ubt_get_n(ubt) result(n)
    type(z_ubt) :: ubt
    n=ubt%n
  end function z_ubt_get_n

  integer(kind=int32) function z_wbv_get_n(wbv) result(n)
    type(z_wbv) :: wbv
    n=wbv%n
  end function z_wbv_get_n

  ! lbwmax

  integer(kind=int32) function d_ub_get_lbwmax(ub) result(n)
    type(d_ub) :: ub
    n=ub%lbwmax
  end function d_ub_get_lbwmax

  integer(kind=int32) function d_bv_get_lbwmax(bv) result(n)
    type(d_bv) :: bv
    n=bv%lbwmax
  end function d_bv_get_lbwmax

  integer(kind=int32) function d_bt_get_lbwmax(bt) result(n)
    type(d_bt) :: bt
    n=bt%lbwmax
  end function d_bt_get_lbwmax

  integer(kind=int32) function d_wb_get_lbwmax(wb) result(n)
    type(d_wb) :: wb
    n=wb%lbwmax
  end function d_wb_get_lbwmax

  integer(kind=int32) function d_ubt_get_lbwmax(ubt) result(n)
    type(d_ubt) :: ubt
    n=ubt%lbwmax
  end function d_ubt_get_lbwmax

  integer(kind=int32) function d_wbv_get_lbwmax(wbv) result(n)
    type(d_wbv) :: wbv
    n=wbv%lbwmax
  end function d_wbv_get_lbwmax

  integer(kind=int32) function z_ub_get_lbwmax(ub) result(n)
    type(z_ub) :: ub
    n=ub%lbwmax
  end function z_ub_get_lbwmax

  integer(kind=int32) function z_bv_get_lbwmax(bv) result(n)
    type(z_bv) :: bv
    n=bv%lbwmax
  end function z_bv_get_lbwmax

  integer(kind=int32) function z_bt_get_lbwmax(bt) result(n)
    type(z_bt) :: bt
    n=bt%lbwmax
  end function z_bt_get_lbwmax

  integer(kind=int32) function z_wb_get_lbwmax(wb) result(n)
    type(z_wb) :: wb
    n=wb%lbwmax
  end function z_wb_get_lbwmax

  integer(kind=int32) function z_ubt_get_lbwmax(ubt) result(n)
    type(z_ubt) :: ubt
    n=ubt%lbwmax
  end function z_ubt_get_lbwmax

  integer(kind=int32) function z_wbv_get_lbwmax(wbv) result(n)
    type(z_wbv) :: wbv
    n=wbv%lbwmax
  end function z_wbv_get_lbwmax

  ! ubwmax

  integer(kind=int32) function d_ub_get_ubwmax(ub) result(n)
    type(d_ub) :: ub
    n=ub%ubwmax
  end function d_ub_get_ubwmax

  integer(kind=int32) function d_bv_get_ubwmax(bv) result(n)
    type(d_bv) :: bv
    n=bv%ubwmax
  end function d_bv_get_ubwmax

  integer(kind=int32) function d_bt_get_ubwmax(bt) result(n)
    type(d_bt) :: bt
    n=bt%ubwmax
  end function d_bt_get_ubwmax

  integer(kind=int32) function d_wb_get_ubwmax(wb) result(n)
    type(d_wb) :: wb
    n=wb%ubwmax
  end function d_wb_get_ubwmax

  integer(kind=int32) function d_ubt_get_ubwmax(ubt) result(n)
    type(d_ubt) :: ubt
    n=ubt%ubwmax
  end function d_ubt_get_ubwmax

  integer(kind=int32) function d_wbv_get_ubwmax(wbv) result(n)
    type(d_wbv) :: wbv
    n=wbv%ubwmax
  end function d_wbv_get_ubwmax

  integer(kind=int32) function z_ub_get_ubwmax(ub) result(n)
    type(z_ub) :: ub
    n=ub%ubwmax
  end function z_ub_get_ubwmax

  integer(kind=int32) function z_bv_get_ubwmax(bv) result(n)
    type(z_bv) :: bv
    n=bv%ubwmax
  end function z_bv_get_ubwmax

  integer(kind=int32) function z_bt_get_ubwmax(bt) result(n)
    type(z_bt) :: bt
    n=bt%ubwmax
  end function z_bt_get_ubwmax

  integer(kind=int32) function z_wb_get_ubwmax(wb) result(n)
    type(z_wb) :: wb
    n=wb%ubwmax
  end function z_wb_get_ubwmax

  integer(kind=int32) function z_ubt_get_ubwmax(ubt) result(n)
    type(z_ubt) :: ubt
    n=ubt%ubwmax
  end function z_ubt_get_ubwmax

  integer(kind=int32) function z_wbv_get_ubwmax(wbv) result(n)
    type(z_wbv) :: wbv
    n=wbv%ubwmax
  end function z_wbv_get_ubwmax

end module mod_orth_band_types
