module nested_types
  use prec
  implicit none

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
  type c_ub
     integer(kind=int32) :: lbw, ubw
     integer(kind=int32), private :: n, lbwmax, ubwmax
     complex(kind=dp), dimension(:,:), allocatable :: bc
     integer(kind=int32), dimension(:), allocatable :: numrotsu
     integer(kind=int32), dimension(:,:), allocatable :: jsu
     complex(kind=dp), dimension(:,:), allocatable :: csu, ssu
  end type c_ub

  ! Stored by rows.
  type c_bt
     integer(kind=int32) :: lbw, ubw
     integer(kind=int32), private :: n, lbwmax, ubwmax
     complex(kind=dp), dimension(:,:), allocatable :: br
     integer(kind=int32), dimension(:), allocatable :: numrotst
     integer(kind=int32), dimension(:,:), allocatable :: kst
     complex(kind=dp), dimension(:,:), allocatable :: cst, sst
  end type c_bt

  ! Stored by columns.
  type c_ubt
     integer(kind=int32) :: lbw, ubw
     integer(kind=int32), private :: n, lbwmax, ubwmax
     complex(kind=dp), dimension(:,:), allocatable :: bc
     integer(kind=int32), dimension(:), allocatable :: numrotsu, numrotst
     integer(kind=int32), dimension(:,:), allocatable :: jsu, kst
     complex(kind=dp), dimension(:,:), allocatable :: csu, ssu, cst, sst
  end type c_ubt

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
  type c_bv
     integer(kind=int32) :: lbw, ubw
     integer(kind=int32), private :: n, lbwmax, ubwmax
     complex(kind=dp), dimension(:,:), allocatable :: br
     integer(kind=int32), dimension(:), allocatable :: numrotsv
     integer(kind=int32), dimension(:,:), allocatable :: ksv
     complex(kind=dp), dimension(:,:), allocatable :: csv, ssv
  end type c_bv

  ! Stored by columns.
  type c_wb
     integer(kind=int32) :: lbw, ubw
     integer(kind=int32), private :: n, lbwmax, ubwmax
     complex(kind=dp), dimension(:,:), allocatable :: bc
     integer(kind=int32), dimension(:), allocatable :: numrotsw
     integer(kind=int32), dimension(:,:), allocatable :: jsw
     complex(kind=dp), dimension(:,:), allocatable :: csw,ssw
  end type c_wb

  ! Stored by rows.
  type c_wbv
     integer(kind=int32) :: lbw, ubw
     integer(kind=int32), private :: n, lbwmax, ubwmax
     complex(kind=dp), dimension(:,:), allocatable :: br
     integer(kind=int32), dimension(:), allocatable :: numrotsv, numrotsw
     integer(kind=int32), dimension(:,:), allocatable :: ksv, jsw
     complex(kind=dp), dimension(:,:), allocatable :: csv, ssv, csw,ssw
  end type c_wbv

  interface deallocate_ub
     module procedure d_deallocate_ub, c_deallocate_ub
  end interface deallocate_ub

  interface deallocate_bt
     module procedure d_deallocate_bt, c_deallocate_bt
  end interface deallocate_bt

  interface deallocate_ubt
     module procedure d_deallocate_ubt, c_deallocate_ubt
  end interface deallocate_ubt

  interface deallocate_bv
     module procedure d_deallocate_bv, c_deallocate_bv
  end interface deallocate_bv

  interface deallocate_wb
     module procedure d_deallocate_wb, c_deallocate_wb
  end interface deallocate_wb

  interface deallocate_wbv
     module procedure d_deallocate_wbv, c_deallocate_wbv
  end interface deallocate_wbv

  interface copy_ub
     module procedure d_copy_ub, c_copy_ub
  end interface copy_ub

  interface copy_bt
     module procedure d_copy_bt, c_copy_bt
  end interface copy_bt

  interface copy_ubt
     module procedure d_copy_ubt, c_copy_ubt
  end interface copy_ubt

  interface copy_bv
     module procedure d_copy_bv, c_copy_bv
  end interface copy_bv

  interface copy_wb
     module procedure d_copy_wb, c_copy_wb
  end interface copy_wb

  interface copy_wbv
     module procedure d_copy_wbv, c_copy_wbv
  end interface copy_wbv

  interface get_n
     module procedure d_ub_get_n, d_bv_get_n, c_ub_get_n, c_bv_get_n, &
          d_ubt_get_n, d_wbv_get_n, c_ubt_get_n, c_wbv_get_n, &
          d_bt_get_n, d_wb_get_n, c_bt_get_n, c_wb_get_n
  end interface get_n

  interface get_lbwmax
     module procedure d_ub_get_lbwmax, d_bv_get_lbwmax, c_ub_get_lbwmax, c_bv_get_lbwmax, &
          d_ubt_get_lbwmax, d_wbv_get_lbwmax, c_ubt_get_lbwmax, c_wbv_get_lbwmax, &
          d_bt_get_lbwmax, d_wb_get_lbwmax, c_bt_get_lbwmax, c_wb_get_lbwmax
  end interface get_lbwmax

  interface get_ubwmax
     module procedure d_ub_get_ubwmax, d_bv_get_ubwmax, c_ub_get_ubwmax, c_bv_get_ubwmax, &
          d_ubt_get_ubwmax, d_wbv_get_ubwmax, c_ubt_get_ubwmax, c_wbv_get_ubwmax, &
          d_bt_get_ubwmax, d_wb_get_ubwmax, c_bt_get_ubwmax, c_wb_get_ubwmax
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
    deallocate(ub%bc, ub%csu, ub%ssu, ub%jsu, ub%numrotsu)
    ub%n=0; ub%lbwmax=0; ub%ubwmax=0
  end subroutine d_deallocate_ub

  type(c_ub) function c_new_ub(n,lbwmax,ubwmax) result(ub)
    integer(kind=int32), intent(in) :: n, lbwmax, ubwmax
    ub%n=n; ub%lbwmax=lbwmax; ub%ubwmax=ubwmax
    ub%ubw=0; ub%lbw=0
    allocate(ub%bc(lbwmax+ubwmax+1,n), ub%csu(ubwmax,n), ub%ssu(ubwmax,n), &
         ub%jsu(ubwmax,n), ub%numrotsu(n))
    ub%bc=(0.0_dp,0.0_dp); ub%csu=(0.0_dp,0.0_dp); ub%ssu=(0.0_dp,0.0_dp)
    ub%jsu=0; ub%numrotsu=0
  end function c_new_ub

  subroutine c_deallocate_ub(ub)
    type(c_ub), intent(inout) :: ub
    deallocate(ub%bc, ub%csu, ub%ssu, ub%jsu, ub%numrotsu)
    ub%n=0; ub%lbwmax=0; ub%ubwmax=0
    ub%lbw=0; ub%ubw=0
  end subroutine c_deallocate_ub

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
    deallocate(bt%br, bt%cst, bt%sst, bt%kst, bt%numrotst)
    bt%n=0; bt%lbwmax=0; bt%ubwmax=0;
    bt%lbw=0; bt%ubw=0
  end subroutine d_deallocate_bt

  type(c_bt) function c_new_bt(n,lbwmax,ubwmax) result(bt)
    integer(kind=int32), intent(in) :: n, lbwmax, ubwmax
    bt%n=n; bt%lbwmax=lbwmax; bt%ubwmax=ubwmax
    bt%ubw=0; bt%lbw=0
    allocate(bt%br(n,lbwmax+ubwmax+1), bt%cst(n,lbwmax), bt%sst(n,lbwmax), &
         bt%kst(n,lbwmax), bt%numrotst(n))
    bt%br=(0.0_dp,0.0_dp)
    bt%cst=(0.0_dp,0.0_dp); bt%sst=(0.0_dp,0.0_dp)
    bt%kst=0; bt%numrotst=0
  end function c_new_bt

  subroutine c_deallocate_bt(bt)
    type(c_bt), intent(inout) :: bt
    deallocate(bt%br, bt%cst, bt%sst, bt%kst, bt%numrotst)
    bt%n=0; bt%lbwmax=0; bt%ubwmax=0;
    bt%lbw=0; bt%ubw=0
  end subroutine c_deallocate_bt


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
    deallocate(ubt%bc, ubt%csu, ubt%ssu, ubt%jsu, ubt%numrotsu, &
         ubt%cst, ubt%sst, ubt%kst, ubt%numrotst)
    ubt%n=0; ubt%lbwmax=0; ubt%ubwmax=0;
    ubt%lbw=0; ubt%ubw=0
  end subroutine d_deallocate_ubt

  type(c_ubt) function c_new_ubt(n,lbwmax,ubwmax) result(ubt)
    integer(kind=int32), intent(in) :: n, lbwmax, ubwmax
    ubt%n=n; ubt%lbwmax=lbwmax; ubt%ubwmax=ubwmax
    ubt%ubw=0; ubt%lbw=0
    allocate(ubt%bc(lbwmax+ubwmax+1,n), ubt%csu(ubwmax,n), ubt%ssu(ubwmax,n), &
         ubt%jsu(ubwmax,n), ubt%numrotsu(n), ubt%cst(n,lbwmax), ubt%sst(n,lbwmax), &
         ubt%kst(n,lbwmax), ubt%numrotst(n))
    ubt%bc=(0.0_dp,0.0_dp)
    ubt%csu=(0.0_dp,0.0_dp); ubt%ssu=(0.0_dp,0.0_dp)
    ubt%jsu=0; ubt%numrotsu=0
    ubt%cst=(0.0_dp,0.0_dp); ubt%sst=(0.0_dp,0.0_dp)
    ubt%kst=0; ubt%numrotst=0
  end function c_new_ubt

  subroutine c_deallocate_ubt(ubt)
    type(c_ubt), intent(inout) :: ubt
    deallocate(ubt%bc, ubt%csu, ubt%ssu, ubt%jsu, ubt%numrotsu, &
         ubt%cst, ubt%sst, ubt%kst, ubt%numrotst)
    ubt%n=0; ubt%lbwmax=0; ubt%ubwmax=0;
    ubt%lbw=0; ubt%ubw=0
  end subroutine c_deallocate_ubt

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
    deallocate(bv%br, bv%csv, bv%ssv, bv%ksv, bv%numrotsv)
    bv%n=0; bv%lbwmax=0; bv%ubwmax=0
    bv%lbw=0; bv%ubw=0
  end subroutine d_deallocate_bv

  type(c_bv) function c_new_bv(n,lbwmax,ubwmax) result(bv)
    integer(kind=int32), intent(in) :: n, lbwmax, ubwmax
    bv%n=n; bv%lbwmax=lbwmax; bv%ubwmax=ubwmax
    bv%ubw=0; bv%lbw=0
    allocate(bv%br(n,lbwmax+ubwmax+1), bv%csv(n,ubwmax), bv%ssv(n,ubwmax), &
         bv%ksv(n,ubwmax), bv%numrotsv(n))
    bv%br=(0.0_dp, 0.0_dp); bv%csv=(0.0_dp, 0.0_dp); bv%ssv=(0.0_dp, 0.0_dp)
    bv%ksv=0; bv%numrotsv=0
  end function c_new_bv

  subroutine c_deallocate_bv(bv)
    type(c_bv), intent(inout) :: bv
    deallocate(bv%br, bv%csv, bv%ssv, bv%ksv, bv%numrotsv)
    bv%n=0; bv%lbwmax=0; bv%ubwmax=0
    bv%lbw=0; bv%ubw=0
  end subroutine c_deallocate_bv

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
    deallocate(wb%bc, wb%csw, wb%ssw, wb%jsw, wb%numrotsw)
    wb%n=0; wb%lbwmax=0; wb%ubwmax=0
    wb%lbw=0; wb%ubw=0
  end subroutine d_deallocate_wb

  type(c_wb) function c_new_wb(n,lbwmax,ubwmax) result(wb)
    integer(kind=int32), intent(in) :: n, lbwmax, ubwmax
    wb%n=n; wb%lbwmax=lbwmax; wb%ubwmax=ubwmax
    wb%ubw=0; wb%lbw=0
    allocate(wb%bc(lbwmax+ubwmax+1,n), wb%csw(lbwmax,n), &
         wb%ssw(lbwmax,n), wb%jsw(lbwmax,n), wb%numrotsw(n))
    wb%bc=(0.0_dp,0.0_dp); wb%csw=(0.0_dp,0.0_dp); wb%ssw=(0.0_dp,0.0_dp)
    wb%jsw=0; wb%numrotsw=0
  end function c_new_wb

  subroutine c_deallocate_wb(wb)
    type(c_wb), intent(inout) :: wb
    deallocate(wb%bc, wb%csw, wb%ssw, wb%jsw, wb%numrotsw)
    wb%n=0; wb%lbwmax=0; wb%ubwmax=0
    wb%lbw=0; wb%ubw=0
  end subroutine c_deallocate_wb

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
    deallocate(wbv%br, wbv%csv, wbv%ssv, wbv%ksv, wbv%numrotsv, &
         wbv%csw, wbv%ssw, wbv%jsw, wbv%numrotsw)
    wbv%n=0; wbv%lbwmax=0; wbv%ubwmax=0
    wbv%lbw=0; wbv%ubw=0
  end subroutine d_deallocate_wbv

  type(c_wbv) function c_new_wbv(n,lbwmax,ubwmax) result(wbv)
    integer(kind=int32), intent(in) :: n, lbwmax, ubwmax
    wbv%n=n; wbv%lbwmax=lbwmax; wbv%ubwmax=ubwmax
    wbv%ubw=0; wbv%lbw=0
    allocate(wbv%br(n,lbwmax+ubwmax+1), wbv%csv(n,ubwmax), wbv%ssv(n,ubwmax), &
         wbv%ksv(n,ubwmax), wbv%numrotsv(n), &
         wbv%csw(lbwmax,n), wbv%ssw(lbwmax,n), wbv%jsw(lbwmax,n), wbv%numrotsw(n))
    wbv%br=(0.0_dp,0.0_dp); wbv%csv=(0.0_dp,0.0_dp); wbv%ssv=(0.0_dp,0.0_dp)
    wbv%ksv=0; wbv%numrotsv=0
    wbv%csw=(0.0_dp,0.0_dp); wbv%ssw=(0.0_dp,0.0_dp)
    wbv%jsw=0; wbv%numrotsw=0
  end function c_new_wbv

  subroutine c_deallocate_wbv(wbv)
    type(c_wbv), intent(inout) :: wbv
    deallocate(wbv%br, wbv%csv, wbv%ssv, wbv%ksv, wbv%numrotsv, &
         wbv%csw, wbv%ssw, wbv%jsw, wbv%numrotsw)
    wbv%n=0; wbv%lbwmax=0; wbv%ubwmax=0
    wbv%lbw=0; wbv%ubw=0
  end subroutine c_deallocate_wbv

  ! copy

  subroutine d_copy_ub(ub1,ub2)
    type(d_ub), intent(in) :: ub1
    type(d_ub) :: ub2
    ub2%n=ub1%n
    ub2%lbw=ub1%lbw
    ub2%ubw=ub1%ubw
    if ( ub2%ubwmax < ub1%ubwmax .or. ub2%ubwmax < ub1%ubwmax .or. &
         size(ub2%bc,2) /= ub1%n ) then
       deallocate(ub2%bc, ub2%jsu, ub2%numrotsu, ub2%csu,ub2%ssu)
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

  subroutine c_copy_ub(ub1,ub2)
    type(c_ub), intent(in) :: ub1
    type(c_ub) :: ub2
    ub2%n=ub1%n
    ub2%lbw=ub1%lbw
    ub2%ubw=ub1%ubw
    if ( ub2%ubwmax < ub1%ubwmax .or. ub2%ubwmax < ub1%ubwmax .or. &
         size(ub2%bc,2) /= ub1%n ) then
       deallocate(ub2%bc, ub2%jsu, ub2%numrotsu, ub2%csu, ub2%ssu)
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
  end subroutine c_copy_ub

  subroutine d_copy_bt(bt1,bt2)
    type(d_bt), intent(in) :: bt1
    type(d_bt) :: bt2
    bt2%n=bt1%n
    bt2%lbw=bt1%lbw
    bt2%ubw=bt1%ubw
    if ( bt2%ubwmax < bt1%ubwmax .or. bt2%ubwmax < bt1%ubwmax .or. &
         size(bt2%br,1) /= bt1%n ) then
       deallocate(bt2%br, bt2%kst, bt2%numrotst, bt2%cst, bt2%sst)
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

  subroutine c_copy_bt(bt1,bt2)
    type(c_bt), intent(in) :: bt1
    type(c_bt) :: bt2
    bt2%n=bt1%n
    bt2%lbw=bt1%lbw
    bt2%ubw=bt1%ubw
    if ( bt2%ubwmax < bt1%ubwmax .or. bt2%ubwmax < bt1%ubwmax .or. &
         size(bt2%br,1) /= bt1%n ) then
       deallocate(bt2%br, bt2%kst, bt2%numrotst, bt2%cst, bt2%sst)
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
  end subroutine c_copy_bt

  subroutine d_copy_ubt(ubt1,ubt2)
    type(d_ubt), intent(in) :: ubt1
    type(d_ubt) :: ubt2
    ubt2%n=ubt1%n
    ubt2%lbw=ubt1%lbw
    ubt2%ubw=ubt1%ubw
    if ( ubt2%ubwmax < ubt1%ubwmax .or. ubt2%ubwmax < ubt1%ubwmax .or. &
         size(ubt2%bc,2) /= ubt1%n ) then
       deallocate(ubt2%bc, ubt2%jsu, ubt2%numrotsu, ubt2%csu,ubt2%ssu, &
            ubt2%kst, ubt2%numrotst, ubt2%cst, ubt2%sst)
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

  subroutine c_copy_ubt(ubt1,ubt2)
    type(c_ubt), intent(in) :: ubt1
    type(c_ubt) :: ubt2
    ubt2%n=ubt1%n
    ubt2%lbw=ubt1%lbw
    ubt2%ubw=ubt1%ubw
    if ( ubt2%ubwmax < ubt1%ubwmax .or. ubt2%ubwmax < ubt1%ubwmax .or. &
         size(ubt2%bc,2) /= ubt1%n ) then
       deallocate(ubt2%bc, ubt2%jsu, ubt2%numrotsu, ubt2%csu,ubt2%ssu, &
            ubt2%kst, ubt2%numrotst, ubt2%cst, ubt2%sst)
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
  end subroutine c_copy_ubt

  subroutine d_copy_bv(bv1,bv2)
    type(d_bv), intent(in) :: bv1
    type(d_bv) :: bv2
    bv2%n=bv1%n
    bv2%lbw=bv1%lbw
    bv2%ubw=bv1%ubw
    if ( bv2%ubwmax < bv1%ubwmax .or. bv2%ubwmax < bv1%ubwmax .or. &
         size(bv2%br,1) /= bv1%n ) then
       deallocate(bv2%br, bv2%ksv, bv2%numrotsv, bv2%csv, bv2%ssv)
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

  subroutine d_copy_wb(wb1,wb2)
    type(d_wb), intent(in) :: wb1
    type(d_wb) :: wb2
    wb2%n=wb1%n
    wb2%lbw=wb1%lbw
    wb2%ubw=wb1%ubw
    if ( wb2%ubwmax < wb1%ubwmax .or. wb2%ubwmax < wb1%ubwmax .or. &
         size(wb2%bc,2) /= wb1%n ) then
       deallocate(wb2%bc, wb2%jsw, wb2%numrotsw, wb2%csw, wb2%ssw)
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
    wb2%jsw(1:wb1%lbw,:)=wb1%jsw(1:wb1%ubw,:)
  end subroutine d_copy_wb

  subroutine d_copy_wbv(wbv1,wbv2)
    type(d_wbv), intent(in) :: wbv1
    type(d_wbv) :: wbv2
    wbv2%n=wbv1%n
    wbv2%lbw=wbv1%lbw
    wbv2%ubw=wbv1%ubw
    if ( wbv2%ubwmax < wbv1%ubwmax .or. wbv2%ubwmax < wbv1%ubwmax .or. &
         size(wbv2%br,1) /= wbv1%n ) then
       deallocate(wbv2%br, wbv2%ksv, wbv2%numrotsv, wbv2%csv, wbv2%ssv, &
            wbv2%jsw, wbv2%numrotsw, wbv2%csw, wbv2%ssw)
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
    wbv2%jsw(1:wbv1%lbw,:)=wbv1%jsw(1:wbv1%ubw,:)
  end subroutine d_copy_wbv

  subroutine c_copy_bv(bv1,bv2)
    type(c_bv), intent(in) :: bv1
    type(c_bv) :: bv2
    bv2%n=bv1%n
    bv2%lbw=bv1%lbw
    bv2%ubw=bv1%ubw

    if ( bv2%ubwmax < bv1%ubwmax .or. bv2%ubwmax < bv1%ubwmax .or. &
         size(bv2%br,1) /= bv1%n ) then
       deallocate(bv2%br, bv2%ksv, bv2%numrotsv, bv2%csv, bv2%ssv)
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
  end subroutine c_copy_bv

  subroutine c_copy_wb(wb1,wb2)
    type(c_wb), intent(in) :: wb1
    type(c_wb) :: wb2
    wb2%n=wb1%n
    wb2%lbw=wb1%lbw
    wb2%ubw=wb1%ubw
    if ( wb2%ubwmax < wb1%ubwmax .or. wb2%ubwmax < wb1%ubwmax .or. &
         size(wb2%bc,2) /= wb1%n ) then
       deallocate(wb2%bc, wb2%jsw, wb2%numrotsw, wb2%csw, wb2%ssw)
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
    wb2%jsw(1:wb1%lbw,:)=wb1%jsw(1:wb1%ubw,:)
  end subroutine c_copy_wb

  subroutine c_copy_wbv(wbv1,wbv2)
    type(c_wbv), intent(in) :: wbv1
    type(c_wbv) :: wbv2
    wbv2%n=wbv1%n
    wbv2%lbw=wbv1%lbw
    wbv2%ubw=wbv1%ubw
    if ( wbv2%ubwmax < wbv1%ubwmax .or. wbv2%ubwmax < wbv1%ubwmax .or. &
         size(wbv2%br,1) /= wbv1%n ) then
       deallocate(wbv2%br, wbv2%ksv, wbv2%numrotsv, wbv2%csv, wbv2%ssv, &
            wbv2%jsw, wbv2%numrotsw, wbv2%csw, wbv2%ssw)
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
    wbv2%jsw(1:wbv1%lbw,:)=wbv1%jsw(1:wbv1%ubw,:)
  end subroutine c_copy_wbv

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

  integer(kind=int32) function c_ub_get_n(ub) result(n)
    type(c_ub) :: ub
    n=ub%n
  end function c_ub_get_n

  integer(kind=int32) function c_bv_get_n(bv) result(n)
    type(c_bv) :: bv
    n=bv%n
  end function c_bv_get_n

  integer(kind=int32) function c_bt_get_n(bt) result(n)
    type(c_bt) :: bt
    n=bt%n
  end function c_bt_get_n

  integer(kind=int32) function c_wb_get_n(wb) result(n)
    type(c_wb) :: wb
    n=wb%n
  end function c_wb_get_n

  integer(kind=int32) function c_ubt_get_n(ubt) result(n)
    type(c_ubt) :: ubt
    n=ubt%n
  end function c_ubt_get_n

  integer(kind=int32) function c_wbv_get_n(wbv) result(n)
    type(c_wbv) :: wbv
    n=wbv%n
  end function c_wbv_get_n

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

  integer(kind=int32) function c_ub_get_lbwmax(ub) result(n)
    type(c_ub) :: ub
    n=ub%lbwmax
  end function c_ub_get_lbwmax

  integer(kind=int32) function c_bv_get_lbwmax(bv) result(n)
    type(c_bv) :: bv
    n=bv%lbwmax
  end function c_bv_get_lbwmax

  integer(kind=int32) function c_bt_get_lbwmax(bt) result(n)
    type(c_bt) :: bt
    n=bt%lbwmax
  end function c_bt_get_lbwmax

  integer(kind=int32) function c_wb_get_lbwmax(wb) result(n)
    type(c_wb) :: wb
    n=wb%lbwmax
  end function c_wb_get_lbwmax

  integer(kind=int32) function c_ubt_get_lbwmax(ubt) result(n)
    type(c_ubt) :: ubt
    n=ubt%lbwmax
  end function c_ubt_get_lbwmax

  integer(kind=int32) function c_wbv_get_lbwmax(wbv) result(n)
    type(c_wbv) :: wbv
    n=wbv%lbwmax
  end function c_wbv_get_lbwmax

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

  integer(kind=int32) function c_ub_get_ubwmax(ub) result(n)
    type(c_ub) :: ub
    n=ub%ubwmax
  end function c_ub_get_ubwmax

  integer(kind=int32) function c_bv_get_ubwmax(bv) result(n)
    type(c_bv) :: bv
    n=bv%ubwmax
  end function c_bv_get_ubwmax

  integer(kind=int32) function c_bt_get_ubwmax(bt) result(n)
    type(c_bt) :: bt
    n=bt%ubwmax
  end function c_bt_get_ubwmax

  integer(kind=int32) function c_wb_get_ubwmax(wb) result(n)
    type(c_wb) :: wb
    n=wb%ubwmax
  end function c_wb_get_ubwmax

  integer(kind=int32) function c_ubt_get_ubwmax(ubt) result(n)
    type(c_ubt) :: ubt
    n=ubt%ubwmax
  end function c_ubt_get_ubwmax

  integer(kind=int32) function c_wbv_get_ubwmax(wbv) result(n)
    type(c_wbv) :: wbv
    n=wbv%ubwmax
  end function c_wbv_get_ubwmax

end module nested_types
