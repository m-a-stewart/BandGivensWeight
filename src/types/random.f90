module mod_random
  use mod_prec
  use mod_utility
  use mod_band_types
  use mod_orth_band_types
  use mod_error_id
  implicit none

  private

  public :: d_random_bc, c_random_bc, f_d_random_bc, f_c_random_bc, f_random_bc, &
       d_random_br, c_random_br, f_d_random_br, f_c_random_br, f_random_br


  public :: d_random_ub, c_random_ub, f_d_random_rotations_ub, f_c_random_rotations_ub, &
       d_random_bv, f_d_random_rotations_bv, c_random_bv, f_c_random_rotations_bv, &
       d_random_bt, f_d_random_rotations_bt, c_random_bt, f_c_random_rotations_bt, &
       d_random_wb, f_d_random_rotations_wb, c_random_wb, f_c_random_rotations_wb, &
       d_random_ubt, c_random_ubt, d_random_wbv, c_random_wbv, &
       f_random_rotations_ub, f_random_rotations_bv, f_random_rotations_bt, &
       f_random_rotations_wb

  interface f_random_rotations_ub
     module procedure f_d_random_rotations_ub, f_c_random_rotations_ub
  end interface f_random_rotations_ub

  interface f_random_rotations_bv
     module procedure f_d_random_rotations_bv, f_c_random_rotations_bv
  end interface f_random_rotations_bv

  interface f_random_rotations_bt
     module procedure f_d_random_rotations_bt, f_c_random_rotations_bt
  end interface f_random_rotations_bt

  interface f_random_rotations_wb
     module procedure f_d_random_rotations_wb, f_c_random_rotations_wb
  end interface f_random_rotations_wb
  
  interface f_random_bc
     module procedure f_d_random_bc, f_c_random_bc
  end interface f_random_bc

  interface f_random_br
     module procedure f_d_random_br, f_c_random_br
  end interface f_random_br

contains

  ! Random band types

  function d_random_bc(n,lbw,ubw,lbwmax0,ubwmax0,error) result(b)
    real(kind=dp), dimension(:,:), allocatable :: b
    integer(kind=int32), intent(in) :: n, lbw, ubw
    integer(kind=int32), intent(in), optional :: lbwmax0, ubwmax0
    type(error_info), optional :: error
    type(routine_info), parameter :: info=info_d_random_bc
    integer(kind=int32) :: lbwmax, ubwmax

    call clear_error(error)
    call push_id(info,error)

    lbwmax=equals_option(lbw,lbwmax0)
    ubwmax=equals_option(ubw,ubwmax0)

    if(n<1) then
       call set_error(1,info,error); return
    end if
    if(lbw > lbwmax) then
       call set_error(2,info,error); return
    end if
    if(ubw > ubwmax) then
       call set_error(3,info,error); return
    end if

    allocate(b(lbwmax+ubwmax+1,n))
    call f_d_random_bc(b,n,lbw,ubw,lbwmax,ubwmax)

    call pop_id(error)
  end function d_random_bc

  subroutine f_d_random_bc(b,n,lbw,ubw,lbwmax,ubwmax)
    real(kind=dp), dimension(lbwmax+ubwmax+1,n), intent(out) :: b
    integer(kind=int32), intent(in) :: n, lbw, ubw
    integer(kind=int32), intent(in) :: lbwmax, ubwmax
    integer(kind=int32) :: k, j0, j1, bw

    b=0.0_dp
    bw=lbw+ubw+1
    do k=1,n
       j0=max(1,ubw-k+2)
       j1=min(bw,ubw+n-k+1)
       call random_matrix(b(j0:j1,k))
    end do
  end subroutine f_d_random_bc

  function c_random_bc(n,lbw,ubw,lbwmax0,ubwmax0,error) result(b)
    complex(kind=dp), dimension(:,:), allocatable :: b
    integer(kind=int32), intent(in) :: n, lbw, ubw
    integer(kind=int32), intent(in), optional :: lbwmax0, ubwmax0
    type(error_info), optional :: error
    type(routine_info), parameter :: info=info_c_random_bc
    integer(kind=int32) :: lbwmax, ubwmax

    call clear_error(error)
    call push_id(info,error)
    lbwmax=equals_option(lbw,lbwmax0)
    ubwmax=equals_option(ubw,ubwmax0)

    if(n<1) then
       call set_error(1,info,error); return
    end if
    if(lbw > lbwmax) then
       call set_error(2,info,error); return
    end if
    if(ubw > ubwmax) then
       call set_error(3,info,error); return
    end if

    allocate(b(lbwmax+ubwmax+1,n))
    call f_c_random_bc(b,n,lbw,ubw,lbwmax,ubwmax)

    call pop_id(error)
  end function c_random_bc

  subroutine f_c_random_bc(b,n,lbw,ubw,lbwmax,ubwmax)
    complex(kind=dp), dimension(lbwmax+ubwmax+1,n), intent(out) :: b
    integer(kind=int32), intent(in) :: n, lbw, ubw
    integer(kind=int32), intent(in) :: lbwmax, ubwmax
    integer(kind=int32) :: k, j0, j1, bw

    b=(0.0_dp,0.0_dp)
    bw=lbw+ubw+1
    do k=1,n
       j0=max(1,ubw-k+2)
       j1=min(bw,ubw+n-k+1)
       call random_matrix(b(j0:j1,k))
    end do
  end subroutine f_c_random_bc

  function d_random_br(n,lbw,ubw,lbwmax0,ubwmax0,error) result(b)
    real(kind=dp), dimension(:,:), allocatable :: b
    integer(kind=int32), intent(in) :: n, lbw, ubw
    integer(kind=int32), intent(in), optional :: lbwmax0, ubwmax0
    type(error_info), optional :: error
    type(routine_info), parameter :: info=info_d_random_br
    integer(kind=int32) :: lbwmax, ubwmax

    call clear_error(error)
    call push_id(info,error)

    lbwmax=equals_option(lbw,lbwmax0)
    ubwmax=equals_option(ubw,ubwmax0)

    if(n<1) then
       call set_error(1,info,error); return
    end if
    if(lbw > lbwmax) then
       call set_error(2,info,error); return
    end if
    if(ubw > ubwmax) then
       call set_error(3,info,error); return
    end if

    allocate(b(n,lbwmax+ubwmax+1))
    call f_d_random_br(b,n,lbw,ubw,lbwmax,ubwmax)
    call pop_id(error)

  end function d_random_br

  subroutine f_d_random_br(b,n,lbw,ubw,lbwmax,ubwmax)
    real(kind=dp), dimension(n,lbwmax+ubwmax+1), intent(out) :: b
    integer(kind=int32), intent(in) :: n, lbw, ubw
    integer(kind=int32), intent(in) :: lbwmax, ubwmax
    integer(kind=int32) :: j, k0, k1, bw

    b=0.0_dp
    bw=lbw+ubw+1
    do j=1,n
       k0=max(1,lbw-j+2)
       k1=min(bw,lbw+n-j+1)
       call random_matrix(b(j,k0:k1))
    end do
  end subroutine f_d_random_br

  function c_random_br(n,lbw,ubw,lbwmax0,ubwmax0,error) result(b)
    complex(kind=dp), dimension(:,:), allocatable :: b
    integer(kind=int32), intent(in) :: n, lbw, ubw
    integer(kind=int32), intent(in), optional :: lbwmax0, ubwmax0
    type(error_info), optional :: error
    type(routine_info), parameter :: info=info_c_random_br
    integer(kind=int32) :: lbwmax, ubwmax

    call clear_error(error)
    call push_id(info,error)

    lbwmax=equals_option(lbw,lbwmax0)
    ubwmax=equals_option(ubw,ubwmax0)

    if(n<1) then
       call set_error(1,info,error); return
    end if
    if(lbw > lbwmax) then
       call set_error(2,info,error); return
    end if
    if(ubw > ubwmax) then
       call set_error(3,info,error); return
    end if

    allocate(b(n,lbwmax+ubwmax+1))
    call f_c_random_br(b,n,lbw,ubw,lbwmax,ubwmax)
    call pop_id(error)

  end function c_random_br

  subroutine f_c_random_br(b,n,lbw,ubw,lbwmax,ubwmax)
    complex(kind=dp), dimension(n,lbwmax+ubwmax+1), intent(out) :: b
    integer(kind=int32), intent(in) :: n, lbw, ubw
    integer(kind=int32), intent(in) :: lbwmax, ubwmax
    integer(kind=int32) :: j, k0, k1, bw

    b=(0.0_dp,0.0_dp)
    bw=lbw+ubw+1
    do j=1,n
       k0=max(1,lbw-j+2)
       k1=min(bw,lbw+n-j+1)
       call random_matrix(b(j,k0:k1))
    end do
  end subroutine f_c_random_br
  
  ! Random ub

  type(d_ub) function d_random_ub(n,lbw,ubw,lbwmax0,ubwmax0,error) result(ub)
    integer(kind=int32), intent(in) :: n, lbw, ubw
    integer(kind=int32), intent(in), optional :: lbwmax0, ubwmax0
    type(error_info), optional :: error
    type(routine_info), parameter :: info=info_d_random_ub
    integer(kind=int32) :: lbwmax, ubwmax

    call clear_error(error)
    call push_id(info,error)
    lbwmax=equals_option(lbw,lbwmax0)
    ubwmax=equals_option(ubw,ubwmax0)

    if(n<1) then
       call set_error(1,info,error); return
    end if
    if(lbw > lbwmax) then
       call set_error(2,info,error); return
    end if
    if(ubw > ubwmax) then
       call set_error(3,info,error); return
    end if
    
    ub=d_new_ub(n,lbwmax,ubwmax)
    ub%ubw=ubw; ub%lbw=lbw
    call f_d_random_bc(ub%bc,n,lbw,ubw,lbwmax,ubwmax)
    call f_d_random_rotations_ub(n,lbw,ubw,lbwmax,ubwmax,ub%numrotsu, ub%jsu, ub%csu, ub%ssu)
    call pop_id(error)
  end function d_random_ub

  subroutine f_d_random_rotations_ub(n,lbw,ubw,lbwmax,ubwmax,numrotsu,jsu,csu,ssu)
    integer(kind=int32), intent(in) :: n, lbw, ubw
    integer(kind=int32), intent(in) :: lbwmax, ubwmax
    integer(kind=int32), dimension(n), intent(out) :: numrotsu
    integer(kind=int32), dimension(ubwmax,n), intent(out) :: jsu
    real(kind=dp), dimension(ubwmax,n), intent(out) :: csu, ssu
    
    integer(kind=int32) :: j,k,nr, bw
    real(kind=dp) :: nrm

    bw=lbw+ubw+1
    csu=0.0_dp; ssu=0.0_dp
    jsu=0; numrotsu=0
    do k=2,n-1
       nr=min(k-1,ubw)
       numrotsu(k)=nr
       call random_matrix(csu(1:nr,k))
       call random_matrix(ssu(1:nr,k))
       jsu(1:nr,k) = (/ (j, j=k-nr, k-1) /)
       do j=1,nr
          nrm = sqrt(csu(j,k)**2 + ssu(j,k)**2)
          csu(j,k)=csu(j,k)/nrm
          ssu(j,k)=ssu(j,k)/nrm
       end do
    end do
  end subroutine f_d_random_rotations_ub

  type(c_ub) function c_random_ub(n,lbw,ubw,lbwmax0,ubwmax0,error) result(ub)
    integer(kind=int32), intent(in) :: n, lbw, ubw
    integer(kind=int32), intent(in), optional :: lbwmax0, ubwmax0
    type(error_info), optional :: error
    type(routine_info), parameter :: info=info_c_random_ub
    integer(kind=int32) :: lbwmax, ubwmax

    call clear_error(error)
    call push_id(info,error)

    lbwmax=equals_option(lbw,lbwmax0)
    ubwmax=equals_option(ubw,ubwmax0)

    if(n<1) then
       call set_error(1,info,error); return
    end if
    if(lbw > lbwmax) then
       call set_error(2,info,error); return
    end if
    if(ubw > ubwmax) then
       call set_error(3,info,error); return
    end if

    ub=c_new_ub(n,lbwmax,ubwmax)
    ub%ubw=ubw; ub%lbw=lbw
    call f_c_random_bc(ub%bc,n,lbw,ubw,lbwmax,ubwmax)
    call f_c_random_rotations_ub(n,lbw,ubw,lbwmax,ubwmax,ub%numrotsu, ub%jsu, ub%csu, ub%ssu)
    call pop_id(error)
  end function c_random_ub

  subroutine f_c_random_rotations_ub(n,lbw,ubw,lbwmax,ubwmax,numrotsu,jsu,csu,ssu)
    integer(kind=int32), intent(in) :: n, lbw, ubw
    integer(kind=int32), intent(in) :: lbwmax, ubwmax
    integer(kind=int32), dimension(n), intent(out) :: numrotsu
    integer(kind=int32), dimension(ubwmax,n), intent(out) :: jsu
    real(kind=dp), dimension(ubwmax,n), intent(out) :: csu
    complex(kind=dp), dimension(ubwmax,n), intent(out) :: ssu
    integer(kind=int32) :: j,k,nr, bw
    real(kind=dp) :: nrm

    bw=lbw+ubw+1
    csu=0.0_dp; ssu=(0.0_dp,0.0_dp)
    jsu=0; numrotsu=0
    do k=2,n-1
       nr=min(k-1,ubw)
       numrotsu(k)=nr
       call random_matrix(csu(1:nr,k))
       call random_matrix(ssu(1:nr,k))
       jsu(1:nr,k) = (/ (j, j=k-nr, k-1) /)
       do j=1,nr
          nrm = sqrt(csu(j,k)**2 + abs(ssu(j,k))**2)
          csu(j,k)=csu(j,k)/nrm
          ssu(j,k)=ssu(j,k)/nrm
       end do
    end do
  end subroutine f_c_random_rotations_ub

  ! Random bv

  type(d_bv) function d_random_bv(n,lbw,ubw,lbwmax0,ubwmax0,error) result(bv)
    integer(kind=int32), intent(in) :: n, lbw, ubw
    integer(kind=int32), intent(in), optional :: lbwmax0, ubwmax0
    type(error_info), optional :: error
    type(routine_info), parameter :: info=info_d_random_bv
    integer(kind=int32) :: lbwmax, ubwmax
    call clear_error(error)
    call push_id(info,error)
    lbwmax=equals_option(lbw,lbwmax0)
    ubwmax=equals_option(ubw,ubwmax0)

    if(n<1) then
       call set_error(1,info,error); return
    end if
    if(lbw > lbwmax) then
       call set_error(2,info,error); return
    end if
    if(ubw > ubwmax) then
       call set_error(3,info,error); return
    end if

    bv=d_new_bv(n,lbwmax,ubwmax)
    bv%ubw=ubw; bv%lbw=lbw
    call f_d_random_br(bv%br,n,lbw,ubw,lbwmax,ubwmax)
    call f_d_random_rotations_bv(n,lbw,ubw,lbwmax,ubwmax,bv%numrotsv, bv%ksv, bv%csv, bv%ssv)
    call pop_id(error)

  end function d_random_bv

  subroutine f_d_random_rotations_bv(n,lbw,ubw,lbwmax,ubwmax,numrotsv,ksv,csv,ssv)
    integer(kind=int32), intent(in) :: n, lbw, ubw
    integer(kind=int32), intent(in) :: lbwmax, ubwmax
    integer(kind=int32), dimension(n), intent(out) :: numrotsv
    integer(kind=int32), dimension(n,ubwmax), intent(out) :: ksv
    real(kind=dp), dimension(n,ubwmax), intent(out) :: csv, ssv
    
    integer(kind=int32) :: j,k,nr, bw
    real(kind=dp) :: nrm

    bw=lbw+ubw+1
    csv=0.0_dp; ssv=0.0_dp
    ksv=0; numrotsv=0
    do j=1,n-2
       nr=min(n-j-1,ubw)
       numrotsv(j)=nr
       call random_matrix(csv(j,1:nr))
       call random_matrix(ssv(j,1:nr))
       ksv(j,1:nr) = (/ (k, k=j+nr,j+1,-1) /)
       do k=1,nr
          nrm = sqrt(csv(j,k)**2 + ssv(j,k)**2)
          csv(j,k)=csv(j,k)/nrm
          ssv(j,k)=ssv(j,k)/nrm
       end do
    end do
  end subroutine f_d_random_rotations_bv

  type(c_bv) function c_random_bv(n,lbw,ubw,lbwmax0,ubwmax0,error) result(bv)
    integer(kind=int32), intent(in) :: n, lbw, ubw
    integer(kind=int32), intent(in), optional :: lbwmax0, ubwmax0
    type(error_info), optional :: error
    type(routine_info), parameter :: info=info_c_random_bv
    integer(kind=int32) :: lbwmax, ubwmax

    call clear_error(error)
    call push_id(info,error)
    lbwmax=equals_option(lbw,lbwmax0)
    ubwmax=equals_option(ubw,ubwmax0)

    if(n<1) then
       call set_error(1,info,error); return
    end if
    if(lbw > lbwmax) then
       call set_error(2,info,error); return
    end if
    if(ubw > ubwmax) then
       call set_error(3,info,error); return
    end if

    bv=c_new_bv(n,lbwmax,ubwmax)
    bv%ubw=ubw; bv%lbw=lbw
    call f_c_random_br(bv%br,n,lbw,ubw,lbwmax,ubwmax)
    call f_c_random_rotations_bv(n,lbw,ubw,lbwmax,ubwmax,bv%numrotsv, bv%ksv, bv%csv, bv%ssv)
    call pop_id(error)
  end function c_random_bv

  subroutine f_c_random_rotations_bv(n,lbw,ubw,lbwmax,ubwmax,numrotsv,ksv,csv,ssv)
    integer(kind=int32), intent(in) :: n, lbw, ubw
    integer(kind=int32), intent(in) :: lbwmax, ubwmax
    integer(kind=int32), dimension(n), intent(out) :: numrotsv
    integer(kind=int32), dimension(n,ubwmax), intent(out) :: ksv
    real(kind=dp), dimension(n,ubwmax), intent(out) :: csv
    complex(kind=dp), dimension(n,ubwmax), intent(out) :: ssv
    
    integer(kind=int32) :: j,k,nr, bw
    real(kind=dp) :: nrm

    bw=lbw+ubw+1
    csv=0.0_dp; ssv=(0.0_dp,0.0_dp)
    ksv=0; numrotsv=0
    do j=1,n-2
       nr=min(n-j-1,ubw)
       numrotsv(j)=nr
       call random_matrix(csv(j,1:nr))
       call random_matrix(ssv(j,1:nr))
       ksv(j,1:nr) = (/ (k, k=j+nr,j+1,-1) /)
       do k=1,nr
          nrm = sqrt(csv(j,k)**2 + abs(ssv(j,k))**2)
          csv(j,k)=csv(j,k)/nrm
          ssv(j,k)=ssv(j,k)/nrm
       end do
    end do
  end subroutine f_c_random_rotations_bv

  ! Random bt

  type(d_bt) function d_random_bt(n,lbw,ubw,lbwmax0,ubwmax0,error) result(bt)
    integer(kind=int32), intent(in) :: n, lbw, ubw
    integer(kind=int32), intent(in), optional :: lbwmax0, ubwmax0
    type(error_info), optional :: error
    type(routine_info), parameter :: info=info_d_random_bt
    integer(kind=int32) :: lbwmax, ubwmax

    call clear_error(error)
    call push_id(info,error)
    lbwmax=equals_option(lbw,lbwmax0)
    ubwmax=equals_option(ubw,ubwmax0)
    if(n<1) then
       call set_error(1,info,error); return
    end if
    if(lbw > lbwmax) then
       call set_error(2,info,error); return
    end if
    if(ubw > ubwmax) then
       call set_error(3,info,error); return
    end if
    bt=d_new_bt(n,lbwmax,ubwmax)
    bt%ubw=ubw; bt%lbw=lbw
    call f_d_random_br(bt%br,n,lbw,ubw,lbwmax,ubwmax)
    call f_d_random_rotations_bt(n,lbw,ubw,lbwmax,ubwmax,bt%numrotst, bt%kst, bt%cst, bt%sst)
    call pop_id(error)
  end function d_random_bt

  subroutine f_d_random_rotations_bt(n,lbw,ubw,lbwmax,ubwmax,numrotst,kst,cst,sst)
    integer(kind=int32), intent(in) :: n, lbw, ubw
    integer(kind=int32), intent(in) :: lbwmax, ubwmax
    integer(kind=int32), dimension(n), intent(out) :: numrotst
    integer(kind=int32), dimension(n,lbwmax), intent(out) :: kst
    real(kind=dp), dimension(n,lbwmax), intent(out) :: cst, sst
    
    integer(kind=int32) :: j,k,nr, bw
    real(kind=dp) :: nrm

    bw=lbw+ubw+1
    cst=0.0_dp; sst=0.0_dp
    kst=0; numrotst=0
    do j=2,n-1
       nr=min(j-1,lbw)
       numrotst(j)=nr
       call random_matrix(cst(j,1:nr))
       call random_matrix(sst(j,1:nr))
       kst(j,1:nr) = (/ (k, k=j-nr,j-1) /)
       do k=1,nr
          nrm = sqrt(cst(j,k)**2 + sst(j,k)**2)
          cst(j,k)=cst(j,k)/nrm
          sst(j,k)=sst(j,k)/nrm
       end do
    end do
  end subroutine f_d_random_rotations_bt

  type(c_bt) function c_random_bt(n,lbw,ubw,lbwmax0,ubwmax0,error) result(bt)
    integer(kind=int32), intent(in) :: n, lbw, ubw
    integer(kind=int32), intent(in), optional :: lbwmax0, ubwmax0
    type(error_info), optional :: error
    type(routine_info), parameter :: info=info_c_random_bt
    integer(kind=int32) :: lbwmax, ubwmax

    call clear_error(error)
    call push_id(info,error)
    lbwmax=equals_option(lbw,lbwmax0)
    ubwmax=equals_option(ubw,ubwmax0)
    if(n<1) then
       call set_error(1,info,error); return
    end if
    if(lbw > lbwmax) then
       call set_error(2,info,error); return
    end if
    if(ubw > ubwmax) then
       call set_error(3,info,error); return
    end if
    bt=c_new_bt(n,lbwmax,ubwmax)
    bt%ubw=ubw; bt%lbw=lbw
    call f_c_random_br(bt%br,n,lbw,ubw,lbwmax,ubwmax)
    call f_c_random_rotations_bt(n,lbw,ubw,lbwmax,ubwmax,bt%numrotst, bt%kst, bt%cst, bt%sst)
    call pop_id(error)
  end function c_random_bt

  subroutine f_c_random_rotations_bt(n,lbw,ubw,lbwmax,ubwmax,numrotst,kst,cst,sst)
    integer(kind=int32), intent(in) :: n, lbw, ubw
    integer(kind=int32), intent(in) :: lbwmax, ubwmax
    integer(kind=int32), dimension(n), intent(out) :: numrotst
    integer(kind=int32), dimension(n,lbwmax), intent(out) :: kst
    real(kind=dp), dimension(n,lbwmax), intent(out) :: cst
    complex(kind=dp), dimension(n,lbwmax), intent(out) :: sst
    
    integer(kind=int32) :: j,k,nr, bw
    real(kind=dp) :: nrm

    bw=lbw+ubw+1
    cst=0.0_dp; sst=(0.0_dp,0.0_dp)
    kst=0; numrotst=0
    do j=2,n-1
       nr=min(j-1,lbw)
       numrotst(j)=nr
       call random_matrix(cst(j,1:nr))
       call random_matrix(sst(j,1:nr))
       kst(j,1:nr) = (/ (k, k=j-nr,j-1) /)
       do k=1,nr
          nrm = sqrt(cst(j,k)**2 + abs(sst(j,k))**2)
          cst(j,k)=cst(j,k)/nrm
          sst(j,k)=sst(j,k)/nrm
       end do
    end do
  end subroutine f_c_random_rotations_bt

  ! Random wb

  type(d_wb) function d_random_wb(n,lbw,ubw,lbwmax0,ubwmax0,error) result(wb)
    integer(kind=int32), intent(in) :: n, lbw, ubw
    integer(kind=int32), intent(in), optional :: lbwmax0, ubwmax0
    type(error_info), optional :: error
    type(routine_info), parameter :: info=info_d_random_wb
    integer(kind=int32) :: lbwmax, ubwmax

    call clear_error(error)
    call push_id(info,error)
    lbwmax=equals_option(lbw,lbwmax0)
    ubwmax=equals_option(ubw,ubwmax0)
    if(n<1) then
       call set_error(1,info,error); return
    end if
    if(lbw > lbwmax) then
       call set_error(2,info,error); return
    end if
    if(ubw > ubwmax) then
       call set_error(3,info,error); return
    end if
    wb=d_new_wb(n,lbwmax,ubwmax)
    wb%ubw=ubw; wb%lbw=lbw
    call f_d_random_bc(wb%bc,n,lbw,ubw,lbwmax,ubwmax)
    call f_d_random_rotations_wb(n,lbw,ubw,lbwmax,ubwmax,wb%numrotsw, wb%jsw, wb%csw, wb%ssw)
    call pop_id(error)
  end function d_random_wb

  subroutine f_d_random_rotations_wb(n,lbw,ubw,lbwmax,ubwmax,numrotsw,jsw,csw,ssw)
    integer(kind=int32), intent(in) :: n, lbw, ubw
    integer(kind=int32), intent(in) :: lbwmax, ubwmax
    integer(kind=int32), dimension(n), intent(out) :: numrotsw
    integer(kind=int32), dimension(lbwmax,n), intent(out) :: jsw
    real(kind=dp), dimension(lbwmax,n), intent(out) :: csw, ssw
    
    integer(kind=int32) :: j,k,nr, bw
    real(kind=dp) :: nrm

    bw=lbw+ubw+1
    csw=0.0_dp; ssw=0.0_dp
    jsw=0; numrotsw=0
    do k=1,n-2
       nr=min(n-k-2,lbw)
       numrotsw(k)=nr
       call random_matrix(csw(1:nr,k))
       call random_matrix(ssw(1:nr,k))
       jsw(1:nr,k) = (/ (j, j=k+nr,k+1,-1) /)
       do j=1,nr
          nrm = sqrt(csw(j,k)**2 + ssw(j,k)**2)
          csw(j,k)=csw(j,k)/nrm
          ssw(j,k)=ssw(j,k)/nrm
       end do
    end do
  end subroutine f_d_random_rotations_wb

  type(c_wb) function c_random_wb(n,lbw,ubw,lbwmax0,ubwmax0,error) result(wb)
    integer(kind=int32), intent(in) :: n, lbw, ubw
    integer(kind=int32), intent(in), optional :: lbwmax0, ubwmax0
    type(error_info), optional :: error
    type(routine_info), parameter :: info=info_c_random_wb
    integer(kind=int32) :: lbwmax, ubwmax

    call clear_error(error)
    call push_id(info,error)
    lbwmax=equals_option(lbw,lbwmax0)
    ubwmax=equals_option(ubw,ubwmax0)
    if(n<1) then
       call set_error(1,info,error); return
    end if
    if(lbw > lbwmax) then
       call set_error(2,info,error); return
    end if
    if(ubw > ubwmax) then
       call set_error(3,info,error); return
    end if
    wb=c_new_wb(n,lbwmax,ubwmax)
    wb%ubw=ubw; wb%lbw=lbw
    call f_c_random_bc(wb%bc,n,lbw,ubw,lbwmax,ubwmax)
    call f_c_random_rotations_wb(n,lbw,ubw,lbwmax,ubwmax,wb%numrotsw, wb%jsw, wb%csw, wb%ssw)
    call pop_id(error)
  end function c_random_wb

  subroutine f_c_random_rotations_wb(n,lbw,ubw,lbwmax,ubwmax,numrotsw,jsw,csw,ssw)
    integer(kind=int32), intent(in) :: n, lbw, ubw
    integer(kind=int32), intent(in) :: lbwmax, ubwmax
    integer(kind=int32), dimension(n), intent(out) :: numrotsw
    integer(kind=int32), dimension(lbwmax,n), intent(out) :: jsw
    real(kind=dp), dimension(lbwmax,n), intent(out) :: csw
    complex(kind=dp), dimension(lbwmax,n), intent(out) :: ssw
    
    integer(kind=int32) :: j,k,nr, bw
    real(kind=dp) :: nrm

    bw=lbw+ubw+1
    csw=0.0_dp; ssw=(0.0_dp,0.0_dp)
    jsw=0; numrotsw=0
    do k=1,n-2
       nr=min(n-k-2,lbw)
       numrotsw(k)=nr
       call random_matrix(csw(1:nr,k))
       call random_matrix(ssw(1:nr,k))
       jsw(1:nr,k) = (/ (j, j=k+nr,k+1,-1) /)
       do j=1,nr
          nrm = sqrt(csw(j,k)**2 + abs(ssw(j,k))**2)
          csw(j,k)=csw(j,k)/nrm
          ssw(j,k)=ssw(j,k)/nrm
       end do
    end do
  end subroutine f_c_random_rotations_wb

  type(d_ubt) function d_random_ubt(n,lbw,ubw,lbwmax0,ubwmax0,error) result(ubt)
    integer(kind=int32), intent(in) :: n, lbw, ubw
    integer(kind=int32), intent(in), optional :: lbwmax0, ubwmax0
    type(error_info), optional :: error
    type(routine_info), parameter :: info=info_d_random_ubt
    integer(kind=int32) :: lbwmax, ubwmax

    call clear_error(error)
    call push_id(info,error)
    lbwmax=equals_option(lbw,lbwmax0)
    ubwmax=equals_option(ubw,ubwmax0)
    if(n<1) then
       call set_error(1,info,error); return
    end if
    if(lbw > lbwmax) then
       call set_error(2,info,error); return
    end if
    if(ubw > ubwmax) then
       call set_error(3,info,error); return
    end if
    ubt=d_new_ubt(n,lbwmax,ubwmax)
    ubt%ubw=ubw; ubt%lbw=lbw
    call f_d_random_bc(ubt%bc,n,lbw,ubw,lbwmax,ubwmax)
    call f_d_random_rotations_ub(n,lbw,ubw,lbwmax,ubwmax,ubt%numrotsu, ubt%jsu, ubt%csu, ubt%ssu)
    call f_d_random_rotations_bt(n,lbw,ubw,lbwmax,ubwmax,ubt%numrotst, ubt%kst, ubt%cst, ubt%sst)
    call pop_id(error)
  end function d_random_ubt

  type(c_ubt) function c_random_ubt(n,lbw,ubw,lbwmax0,ubwmax0,error) result(ubt)
    integer(kind=int32), intent(in) :: n, lbw, ubw
    integer(kind=int32), intent(in), optional :: lbwmax0, ubwmax0
    type(error_info), optional :: error
    type(routine_info), parameter :: info=info_c_random_ubt
    integer(kind=int32) :: lbwmax, ubwmax

    call clear_error(error)
    call push_id(info,error)
    lbwmax=equals_option(lbw,lbwmax0)
    ubwmax=equals_option(ubw,ubwmax0)
    if(n<1) then
       call set_error(1,info,error); return
    end if
    if(lbw > lbwmax) then
       call set_error(2,info,error); return
    end if
    if(ubw > ubwmax) then
       call set_error(3,info,error); return
    end if
    ubt=c_new_ubt(n,lbwmax,ubwmax)
    ubt%ubw=ubw; ubt%lbw=lbw
    call f_c_random_bc(ubt%bc,n,lbw,ubw,lbwmax,ubwmax)
    call f_c_random_rotations_ub(n,lbw,ubw,lbwmax,ubwmax,ubt%numrotsu, ubt%jsu, ubt%csu, ubt%ssu)
    call f_c_random_rotations_bt(n,lbw,ubw,lbwmax,ubwmax,ubt%numrotst, ubt%kst, ubt%cst, ubt%sst)
    call pop_id(error)
  end function c_random_ubt

  type(d_wbv) function d_random_wbv(n,lbw,ubw,lbwmax0,ubwmax0,error) result(wbv)
    integer(kind=int32), intent(in) :: n, lbw, ubw
    integer(kind=int32), intent(in), optional :: lbwmax0, ubwmax0
    type(error_info), optional :: error
    type(routine_info), parameter :: info=info_d_random_wbv
    integer(kind=int32) :: lbwmax, ubwmax

    call clear_error(error)
    call push_id(info,error)
    lbwmax=equals_option(lbw,lbwmax0)
    ubwmax=equals_option(ubw,ubwmax0)
    if(n<1) then
       call set_error(1,info,error); return
    end if
    if(lbw > lbwmax) then
       call set_error(2,info,error); return
    end if
    if(ubw > ubwmax) then
       call set_error(3,info,error); return
    end if
    wbv=d_new_wbv(n,lbwmax,ubwmax)
    wbv%ubw=ubw; wbv%lbw=lbw
    call f_d_random_br(wbv%br,n,lbw,ubw,lbwmax,ubwmax)
    call f_d_random_rotations_wb(n,lbw,ubw,lbwmax,ubwmax,wbv%numrotsw, wbv%jsw, wbv%csw, wbv%ssw)
    call f_d_random_rotations_bv(n,lbw,ubw,lbwmax,ubwmax,wbv%numrotsv, wbv%ksv, wbv%csv, wbv%ssv)
    call pop_id(error)
  end function d_random_wbv

  type(c_wbv) function c_random_wbv(n,lbw,ubw,lbwmax0,ubwmax0,error) result(wbv)
    integer(kind=int32), intent(in) :: n, lbw, ubw
    integer(kind=int32), intent(in), optional :: lbwmax0, ubwmax0
    type(error_info), optional :: error
    type(routine_info), parameter :: info=info_c_random_wbv
    integer(kind=int32) :: lbwmax, ubwmax

    call clear_error(error)
    call push_id(info,error)
    lbwmax=equals_option(lbw,lbwmax0)
    ubwmax=equals_option(ubw,ubwmax0)
    if(n<1) then
       call set_error(1,info,error); return
    end if
    if(lbw > lbwmax) then
       call set_error(2,info,error); return
    end if
    if(ubw > ubwmax) then
       call set_error(3,info,error); return
    end if
    wbv=c_new_wbv(n,lbwmax,ubwmax)
    wbv%ubw=ubw; wbv%lbw=lbw
    call f_c_random_br(wbv%br,n,lbw,ubw,lbwmax,ubwmax)
    call f_c_random_rotations_wb(n,lbw,ubw,lbwmax,ubwmax,wbv%numrotsw, wbv%jsw, wbv%csw, wbv%ssw)
    call f_c_random_rotations_bv(n,lbw,ubw,lbwmax,ubwmax,wbv%numrotsv, wbv%ksv, wbv%csv, wbv%ssv)
    call pop_id(error)
  end function c_random_wbv

  
end module mod_random
