module mod_products
  use mod_prec
  use mod_error_id
  use mod_band_types
  use mod_orth_band_types
  use mod_rotation
  use mod_utility
  implicit none
  ! This module contains routines for forming fast products of Givens
  ! weight representations with general matrices.

  private

  public :: ub_times_general, product_of_ub_and_general, &
       d_ub_times_general, d_product_of_ub_and_general, &
       z_ub_times_general, z_product_of_ub_and_general, &
       d_bv_times_general, d_product_of_bv_and_general, &
       z_bv_times_general, z_product_of_bv_and_general, &
       d_wb_times_general, d_product_of_wb_and_general, &
       z_wb_times_general, z_product_of_wb_and_general, &
       d_bt_times_general, d_product_of_bt_and_general, &
       z_bt_times_general, z_product_of_bt_and_general, &
       f_d_ub_times_general_plus, f_z_ub_times_general_plus, &
       f_d_bv_times_general_plus, f_z_bv_times_general_plus, &
       f_d_wb_times_general_plus, f_z_wb_times_general_plus, &
       f_d_bt_times_general_plus, f_z_bt_times_general_plus, &       
       operator(*)
       
  interface ub_times_general
     module procedure d_ub_times_general, z_ub_times_general
  end interface ub_times_general

  interface f_ub_times_general
     module procedure f_d_ub_times_general_plus, f_z_ub_times_general_plus
  end interface f_ub_times_general
  
  interface product_of_ub_and_general
     module procedure d_product_of_ub_and_general, z_product_of_ub_and_general
  end interface product_of_ub_and_general

  interface bv_times_general
     module procedure d_bv_times_general, z_bv_times_general
  end interface bv_times_general

  interface f_bv_times_general
     module procedure f_d_bv_times_general_plus, f_z_bv_times_general_plus
  end interface f_bv_times_general
  
  interface product_of_bv_and_general
     module procedure d_product_of_bv_and_general, z_product_of_bv_and_general
  end interface product_of_bv_and_general

  interface wb_times_general
     module procedure d_wb_times_general, z_wb_times_general
  end interface wb_times_general

  interface f_wb_times_general
     module procedure f_d_wb_times_general_plus, f_z_wb_times_general_plus
  end interface f_wb_times_general
  
  interface product_of_wb_and_general
     module procedure d_product_of_wb_and_general, z_product_of_wb_and_general
  end interface product_of_wb_and_general

  interface bt_times_general
     module procedure d_bt_times_general, z_bt_times_general
  end interface bt_times_general

  interface f_bt_times_general
     module procedure f_d_bt_times_general_plus, f_z_bt_times_general_plus
  end interface f_bt_times_general
  
  interface product_of_bt_and_general
     module procedure d_product_of_bt_and_general, z_product_of_bt_and_general
  end interface product_of_bt_and_general

  
  interface operator (*)
     module procedure d_product_of_ub_and_general0,  z_product_of_ub_and_general0, &
          d_product_of_bv_and_general0, z_product_of_bv_and_general0, &
          d_product_of_wb_and_general0, z_product_of_wb_and_general0, &
          d_product_of_bt_and_general0, z_product_of_bt_and_general0          
  end interface operator (*)

contains

  ! Real ub times general
  function d_product_of_ub_and_general0(ub, a) result(c)
    real(kind=dp), dimension(:,:), intent(in) :: a
    real(kind=dp), dimension(:,:), allocatable :: c
    type(d_ub), intent(in) :: ub
    
    allocate(c(size(a,1),size(a,2)))
    call d_ub_times_general(ub,a,c)
  end function d_product_of_ub_and_general0

  function d_product_of_ub_and_general(ub, a, error) result(c)
    real(kind=dp), dimension(:,:), intent(in) :: a
    real(kind=dp), dimension(:,:), allocatable :: c
    type(d_ub), intent(in) :: ub
    type(error_info), intent(inout), optional :: error
    type(routine_info), parameter :: info=info_d_product_of_ub_and_general

    if (failure(error)) return
    call push_id(info, error)
    
    allocate(c(size(a,1),size(a,2)))
    call d_ub_times_general(ub,a,c,error)
    call pop_id(error)
  end function d_product_of_ub_and_general

  ! Errors
  ! 0: no error
  ! 1 ub%n < 1
  ! 2: ub%n /= size(a,1)
  ! 3: size(a) /= size(c)
  subroutine d_ub_times_general(ub,a,c,error)
    type(d_ub), intent(in) :: ub
    real(kind=dp), target, dimension(:,:), intent(in) :: a
    real(kind=dp), target, dimension(:,:), intent(out) :: c    
    type(error_info), intent(inout), optional :: error
    type(routine_info), parameter :: info=info_d_ub_times_general

    if (failure(error)) return
    call push_id(info, error)
    if (get_n(ub) < 1) then
       call set_error(1, info, error); return
    end if
    if (get_n(ub) /= size(a,1)) then
       call set_error(2, info, error); return
    end if
    if (size(a,1) /= size(c,1) .or. size(a,2) /= size(c,2)) then
       call set_error(3, info, error); return
    end if
    c=0.0_dp
    call f_d_ub_times_general_plus(ub%bc, get_n(ub), ub%lbw, ub%ubw, get_lbwmax(ub), &
         get_ubwmax(ub), ub%numrotsu, ub%jsu, ub%csu, ub%ssu, a, size(a,2), c)
    call pop_id(error)
  end subroutine d_ub_times_general

  subroutine f_d_ub_times_general_plus(bc, n, lbw, ubw, lbwmax, ubwmax, numrotsu, &
       jsu, csu, ssu, a, na, c)
    real(kind=dp), target, dimension(n,na), intent(in) :: a
    real(kind=dp), target, dimension(n,na), intent(out) :: c    
    integer(kind=int32), dimension(ubwmax,n), intent(in) :: jsu
    real(kind=dp), dimension(ubwmax,n), intent(in) :: csu, ssu
    real(kind=dp), dimension(lbwmax+ubwmax+1,n), intent(in) :: bc
    integer(kind=int32), dimension(n), intent(in) :: numrotsu
    integer(kind=int32), intent(in) :: ubw, lbw, n, lbwmax, ubwmax, na
    !
    integer(kind=int32) :: j,k,l,lbwk,ubwk,d0
    real(kind=dp), dimension(n) :: x
    type(d_rotation) :: rot
    d0=ubw+1
    do j=1,na
       x=0.0_dp
       do k=n,3,-1
          ubwk=min(ubw,k-1) ! number of superdiagonals in column k
          lbwk=min(lbw,n-k) ! number of subdiagonals in column k
          x(k-ubwk:k+lbwk)=x(k-ubwk:k+lbwk)+bc(d0-ubwk:d0+lbwk,k)*a(k,j)
          ! apply u_{k-1} to x
          do l=1,numrotsu(k-1)
             rot%cosine=csu(l,k-1); rot%sine=ssu(l,k-1)
             call rotation_times_general(rot,x,jsu(l,k-1),jsu(l,k-1)+1)
          end do
       end do
       do k=2,1,-1
          ubwk=min(ubw,k-1) ! number of superdiagonals in column k
          lbwk=min(lbw,n-k) ! number of subdiagonals in column k
          x(k-ubwk:k+lbwk)=x(k-ubwk:k+lbwk)+bc(d0-ubwk:d0+lbwk,k)*a(k,j)
       end do
       c(:,j)=c(:,j)+x
    end do
  end subroutine f_d_ub_times_general_plus

  ! Complex ub times general

  function z_product_of_ub_and_general0(ub, a) result(c)
    complex(kind=dp), dimension(:,:), intent(in) :: a
    complex(kind=dp), dimension(:,:), allocatable :: c
    type(z_ub), intent(in) :: ub
    
    allocate(c(size(a,1),size(a,2)))
    call z_ub_times_general(ub,a,c)
  end function z_product_of_ub_and_general0

  function z_product_of_ub_and_general(ub, a, error) result(c)
    complex(kind=dp), dimension(:,:), intent(in) :: a
    complex(kind=dp), dimension(:,:), allocatable :: c
    type(z_ub), intent(in) :: ub
    type(error_info), intent(inout), optional :: error
    type(routine_info), parameter :: info=info_z_product_of_ub_and_general

    if (failure(error)) return
    call push_id(info, error)
    
    allocate(c(size(a,1),size(a,2)))
    call z_ub_times_general(ub,a,c,error)
    call pop_id(error)
  end function z_product_of_ub_and_general

  ! Errors
  ! 0: no error
  ! 1 ub%n < 1
  ! 2: ub%n /= size(a,1)
  ! 3: size(a) /= size(c)
  subroutine z_ub_times_general(ub,a,c,error)
    type(z_ub), intent(in) :: ub
    complex(kind=dp), target, dimension(:,:), intent(in) :: a
    complex(kind=dp), target, dimension(:,:), intent(out) :: c    
    type(error_info), intent(inout), optional :: error
    type(routine_info), parameter :: info=info_z_ub_times_general

    if (failure(error)) return
    call push_id(info, error)
    if (get_n(ub) < 1) then
       call set_error(1, info, error); return
    end if
    if (get_n(ub) /= size(a,1)) then
       call set_error(2, info, error); return
    end if
    if (size(a,1) /= size(c,1) .or. size(a,2) /= size(c,2)) then
       call set_error(3, info, error); return
    end if
    c=(0.0_dp,0.0_dp)
    call f_z_ub_times_general_plus(ub%bc, get_n(ub), ub%lbw, ub%ubw, get_lbwmax(ub), &
         get_ubwmax(ub), ub%numrotsu, ub%jsu, ub%csu, ub%ssu, a, size(a,2), c)
    call pop_id(error)
  end subroutine z_ub_times_general

  subroutine f_z_ub_times_general_plus(bc, n, lbw, ubw, lbwmax, ubwmax, numrotsu, &
       jsu, csu, ssu, a, na, c)
    complex(kind=dp), target, dimension(n,na), intent(in) :: a
    complex(kind=dp), target, dimension(n,na), intent(out) :: c    
    integer(kind=int32), dimension(ubwmax,n), intent(in) :: jsu
    real(kind=dp), dimension(ubwmax,n), intent(in) :: csu
    complex(kind=dp), dimension(ubwmax,n), intent(in) :: ssu    
    complex(kind=dp), dimension(lbwmax+ubwmax+1,n), intent(in) :: bc
    integer(kind=int32), dimension(n), intent(in) :: numrotsu
    integer(kind=int32), intent(in) :: ubw, lbw, n, lbwmax, ubwmax, na
    !
    integer(kind=int32) :: j,k,l,lbwk,ubwk,d0
    complex(kind=dp), dimension(n) :: x
    type(z_rotation) :: rot
    d0=ubw+1
    do j=1,na
       x=0.0_dp
       do k=n,3,-1
          ubwk=min(ubw,k-1) ! number of superdiagonals in column k
          lbwk=min(lbw,n-k) ! number of subdiagonals in column k
          x(k-ubwk:k+lbwk)=x(k-ubwk:k+lbwk)+bc(d0-ubwk:d0+lbwk,k)*a(k,j)
          ! apply u_{k-1} to x
          do l=1,numrotsu(k-1)
             rot%cosine=csu(l,k-1); rot%sine=ssu(l,k-1)
             call rotation_times_general(rot,x,jsu(l,k-1),jsu(l,k-1)+1)
          end do
       end do
       do k=2,1,-1
          ubwk=min(ubw,k-1) ! number of superdiagonals in column k
          lbwk=min(lbw,n-k) ! number of subdiagonals in column k
          x(k-ubwk:k+lbwk)=x(k-ubwk:k+lbwk)+bc(d0-ubwk:d0+lbwk,k)*a(k,j)
       end do
       c(:,j)=c(:,j)+x
    end do
  end subroutine f_z_ub_times_general_plus
  
  ! real BV times general

  function d_product_of_bv_and_general0(bv, a) result(c)
    real(kind=dp), dimension(:,:), intent(in) :: a
    real(kind=dp), dimension(:,:), allocatable :: c
    type(d_bv), intent(in) :: bv
    
    allocate(c(size(a,1),size(a,2)))
    call d_bv_times_general(bv,a,c)
  end function d_product_of_bv_and_general0

  function d_product_of_bv_and_general(bv, a, error) result(c)
    real(kind=dp), dimension(:,:), intent(in) :: a
    real(kind=dp), dimension(:,:), allocatable :: c
    type(d_bv), intent(in) :: bv
    type(error_info), intent(inout), optional :: error
    type(routine_info), parameter :: info=info_d_product_of_bv_and_general

    if (failure(error)) return
    call push_id(info, error)
    
    allocate(c(size(a,1),size(a,2)))
    call d_bv_times_general(bv,a,c,error)
    call pop_id(error)
  end function d_product_of_bv_and_general

  ! Errors
  ! 0: no error
  ! 1 bv%n < 1
  ! 2: bv%n /= size(a,1)
  ! 3: size(a) /= size(c)
  subroutine d_bv_times_general(bv,a,c,error)
    type(d_bv), intent(in) :: bv
    real(kind=dp), target, dimension(:,:), intent(in) :: a
    real(kind=dp), target, dimension(:,:), intent(out) :: c    
    type(error_info), intent(inout), optional :: error
    type(routine_info), parameter :: info=info_d_bv_times_general

    if (failure(error)) return
    call push_id(info, error)
    if (get_n(bv) < 1) then
       call set_error(1, info, error); return
    end if
    if (get_n(bv) /= size(a,1)) then
       call set_error(2, info, error); return
    end if
    if (size(a,1) /= size(c,1) .or. size(a,2) /= size(c,2)) then
       call set_error(3, info, error); return
    end if
    c=0.0_dp
    call f_d_bv_times_general_plus(bv%br, get_n(bv), bv%lbw, bv%ubw, get_lbwmax(bv), &
         get_ubwmax(bv), bv%numrotsv, bv%ksv, bv%csv, bv%ssv, a, size(a,2), c)
    call pop_id(error)
  end subroutine d_bv_times_general

  subroutine f_d_bv_times_general_plus(br, n, lbw, ubw, lbwmax, ubwmax, numrotsv, &
       ksv, csv, ssv, a, na, c)
    real(kind=dp), target, dimension(n,na), intent(in) :: a
    real(kind=dp), target, dimension(n,na), intent(out) :: c    
    integer(kind=int32), dimension(n,ubwmax), intent(in) :: ksv
    real(kind=dp), dimension(n,ubwmax), intent(in) :: csv, ssv
    real(kind=dp), dimension(n,lbwmax+ubwmax+1), intent(in) :: br
    integer(kind=int32), dimension(n), intent(in) :: numrotsv
    integer(kind=int32), intent(in) :: ubw, lbw, n, lbwmax, ubwmax, na
    !
    integer(kind=int32) :: j,k,l,lbwj,ubwj,d0
    real(kind=dp), dimension(n) :: x
    type(d_rotation) :: rot
    d0=lbw+1
    do k=1,na
       x=a(:,k)
       lbwj=min(lbw,n-1) ! number of subdiagonals in row n
       c(n,k)=c(n,k)+dot_product(br(n,d0-lbwj:d0),x(n-lbwj:n))
       lbwj=min(lbw,n-2) ! number of subdiagonals in row n-1
       ubwj=min(ubw,1)
       c(n-1,k)=c(n-1,k)+dot_product(br(n-1,d0-lbwj:d0+ubwj),x(n-1-lbwj:n-1+ubwj))
       do j=n-2,1,-1
          ! apply V_j^T to x
          do l=numrotsv(j),1,-1
             rot%cosine=csv(j,l); rot%sine=ssv(j,l)
             call rotation_times_general(trp_rot(rot),x,ksv(j,l),ksv(j,l)+1)
          end do
          lbwj=min(lbw,j-1) ! number of subdiagonals in row j
          ubwj=min(ubw,n-j)
          c(j,k)=c(j,k)+dot_product(br(j,d0-lbwj:d0+ubwj),x(j-lbwj:j+ubwj))
       end do
    end do
  end subroutine f_d_bv_times_general_plus


  ! Complex BV times general

  function z_product_of_bv_and_general0(bv, a) result(c)
    complex(kind=dp), dimension(:,:), intent(in) :: a
    complex(kind=dp), dimension(:,:), allocatable :: c
    type(z_bv), intent(in) :: bv
    
    allocate(c(size(a,1),size(a,2)))
    call z_bv_times_general(bv,a,c)
  end function z_product_of_bv_and_general0

  function z_product_of_bv_and_general(bv, a, error) result(c)
    complex(kind=dp), dimension(:,:), intent(in) :: a
    complex(kind=dp), dimension(:,:), allocatable :: c
    type(z_bv), intent(in) :: bv
    type(error_info), intent(inout), optional :: error
    type(routine_info), parameter :: info=info_z_product_of_bv_and_general

    if (failure(error)) return
    call push_id(info, error)
    
    allocate(c(size(a,1),size(a,2)))
    call z_bv_times_general(bv,a,c,error)
    call pop_id(error)
  end function z_product_of_bv_and_general

  ! Errors
  ! 0: no error
  ! 1 bv%n < 1
  ! 2: bv%n /= size(a,1)
  ! 3: size(a) /= size(c)
  subroutine z_bv_times_general(bv,a,c,error)
    type(z_bv), intent(in) :: bv
    complex(kind=dp), target, dimension(:,:), intent(in) :: a
    complex(kind=dp), target, dimension(:,:), intent(out) :: c    
    type(error_info), intent(inout), optional :: error
    type(routine_info), parameter :: info=info_z_bv_times_general

    if (failure(error)) return
    call push_id(info, error)
    if (get_n(bv) < 1) then
       call set_error(1, info, error); return
    end if
    if (get_n(bv) /= size(a,1)) then
       call set_error(2, info, error); return
    end if
    if (size(a,1) /= size(c,1) .or. size(a,2) /= size(c,2)) then
       call set_error(3, info, error); return
    end if
    c=(0.0_dp,0.0_dp)
    call f_z_bv_times_general_plus(bv%br, get_n(bv), bv%lbw, bv%ubw, get_lbwmax(bv), &
         get_ubwmax(bv), bv%numrotsv, bv%ksv, bv%csv, bv%ssv, a, size(a,2), c)
    call pop_id(error)
  end subroutine z_bv_times_general

  subroutine f_z_bv_times_general_plus(br, n, lbw, ubw, lbwmax, ubwmax, numrotsv, &
       ksv, csv, ssv, a, na, c)
    complex(kind=dp), target, dimension(n,na), intent(in) :: a
    complex(kind=dp), target, dimension(n,na), intent(out) :: c    
    integer(kind=int32), dimension(n,ubwmax), intent(in) :: ksv
    real(kind=dp), dimension(n,ubwmax), intent(in) :: csv
    complex(kind=dp), dimension(n,ubwmax), intent(in) :: ssv    
    complex(kind=dp), dimension(n,lbwmax+ubwmax+1), intent(in) :: br
    integer(kind=int32), dimension(n), intent(in) :: numrotsv
    integer(kind=int32), intent(in) :: ubw, lbw, n, lbwmax, ubwmax, na
    !
    integer(kind=int32) :: j,k,l,lbwj,ubwj,d0
    complex(kind=dp), dimension(n) :: x
    type(z_rotation) :: rot
    d0=lbw+1
    do k=1,na
       x=a(:,k)
       lbwj=min(lbw,n-1) ! number of subdiagonals in row n
       c(n,k)=c(n,k)+dot_no_conjg(br(n,d0-lbwj:d0),x(n-lbwj:n))
       lbwj=min(lbw,n-2) ! number of subdiagonals in row n-1
       ubwj=min(ubw,1)
       c(n-1,k)=c(n-1,k)+dot_no_conjg(br(n-1,d0-lbwj:d0+ubwj),x(n-1-lbwj:n-1+ubwj))
       do j=n-2,1,-1
          ! apply V_j^T to x
          do l=numrotsv(j),1,-1
             rot%cosine=csv(j,l); rot%sine=ssv(j,l)
             call rotation_times_general(trp_rot(rot),x,ksv(j,l),ksv(j,l)+1)
          end do
          lbwj=min(lbw,j-1) ! number of subdiagonals in row j
          ubwj=min(ubw,n-j)
          c(j,k)=c(j,k)+dot_no_conjg(br(j,d0-lbwj:d0+ubwj),x(j-lbwj:j+ubwj))
       end do
    end do
  end subroutine f_z_bv_times_general_plus

  ! WB times general

  function d_product_of_wb_and_general0(wb, a) result(c)
    real(kind=dp), dimension(:,:), intent(in) :: a
    real(kind=dp), dimension(:,:), allocatable :: c
    type(d_wb), intent(in) :: wb
    
    allocate(c(size(a,1),size(a,2)))
    call d_wb_times_general(wb,a,c)
  end function d_product_of_wb_and_general0

  function d_product_of_wb_and_general(wb, a, error) result(c)
    real(kind=dp), dimension(:,:), intent(in) :: a
    real(kind=dp), dimension(:,:), allocatable :: c
    type(d_wb), intent(in) :: wb
    type(error_info), intent(inout), optional :: error
    type(routine_info), parameter :: info=info_d_product_of_wb_and_general

    if (failure(error)) return
    call push_id(info, error)
    
    allocate(c(size(a,1),size(a,2)))
    call d_wb_times_general(wb,a,c,error)
    call pop_id(error)
  end function d_product_of_wb_and_general

  ! Errors
  ! 0: no error
  ! 1 wb%n < 1
  ! 2: wb%n /= size(a,1)
  ! 3: size(a) /= size(c)
  subroutine d_wb_times_general(wb,a,c,error)
    type(d_wb), intent(in) :: wb
    real(kind=dp), target, dimension(:,:), intent(in) :: a
    real(kind=dp), target, dimension(:,:), intent(out) :: c    
    type(error_info), intent(inout), optional :: error
    type(routine_info), parameter :: info=info_d_wb_times_general

    if (failure(error)) return
    call push_id(info, error)
    if (get_n(wb) < 1) then
       call set_error(1, info, error); return
    end if
    if (get_n(wb) /= size(a,1)) then
       call set_error(2, info, error); return
    end if
    if (size(a,1) /= size(c,1) .or. size(a,2) /= size(c,2)) then
       call set_error(3, info, error); return
    end if
    c=0.0_dp
    call f_d_wb_times_general_plus(wb%bc, get_n(wb), wb%lbw, wb%ubw, get_lbwmax(wb), &
         get_ubwmax(wb), wb%numrotsw, wb%jsw, wb%csw, wb%ssw, a, size(a,2), c)
    call pop_id(error)
  end subroutine d_wb_times_general
  

  subroutine f_d_wb_times_general_plus(bc, n, lbw, ubw, lbwmax, ubwmax, numrotsw, &
       jsw, csw, ssw, a, na, c)
    real(kind=dp), target, dimension(n,na), intent(in) :: a
    real(kind=dp), target, dimension(n,na), intent(out) :: c    
    integer(kind=int32), dimension(lbwmax,n), intent(in) :: jsw
    real(kind=dp), dimension(lbwmax,n), intent(in) :: csw, ssw
    real(kind=dp), dimension(lbwmax+ubwmax+1,n), intent(in) :: bc
    integer(kind=int32), dimension(n), intent(in) :: numrotsw
    integer(kind=int32), intent(in) :: ubw, lbw, n, lbwmax, ubwmax, na
    !
    integer(kind=int32) :: j,k,l,lbwk,ubwk,d0
    real(kind=dp), dimension(n) :: x
    type(d_rotation) :: rot
    d0=ubw+1
    do j=1,na
       x=0.0_dp
       do k=1,n-2
          ubwk=min(ubw,k-1) ! number of superdiagonals in column k
          lbwk=min(lbw,n-k) ! number of subdiagonals in column k
          x(k-ubwk:k+lbwk)=x(k-ubwk:k+lbwk)+bc(d0-ubwk:d0+lbwk,k)*a(k,j)
          ! apply w_{k} to x
          do l=1,numrotsw(k)
             rot%cosine=csw(l,k); rot%sine=ssw(l,k)
             call rotation_times_general(rot,x,jsw(l,k),jsw(l,k)+1)
          end do
       end do
       do k=n-1,n
          ubwk=min(ubw,k-1) ! number of superdiagonals in column k
          lbwk=min(lbw,n-k) ! number of subdiagonals in column k
          x(k-ubwk:k+lbwk)=x(k-ubwk:k+lbwk)+bc(d0-ubwk:d0+lbwk,k)*a(k,j)
       end do
       c(:,j)=c(:,j)+x
    end do
  end subroutine f_d_wb_times_general_plus

  function z_product_of_wb_and_general0(wb, a) result(c)
    complex(kind=dp), dimension(:,:), intent(in) :: a
    complex(kind=dp), dimension(:,:), allocatable :: c
    type(z_wb), intent(in) :: wb
    
    allocate(c(size(a,1),size(a,2)))
    call z_wb_times_general(wb,a,c)
  end function z_product_of_wb_and_general0

  function z_product_of_wb_and_general(wb, a, error) result(c)
    complex(kind=dp), dimension(:,:), intent(in) :: a
    complex(kind=dp), dimension(:,:), allocatable :: c
    type(z_wb), intent(in) :: wb
    type(error_info), intent(inout), optional :: error
    type(routine_info), parameter :: info=info_z_product_of_wb_and_general

    if (failure(error)) return
    call push_id(info, error)
    
    allocate(c(size(a,1),size(a,2)))
    call z_wb_times_general(wb,a,c,error)
    call pop_id(error)
  end function z_product_of_wb_and_general

  ! Errors
  ! 0: no error
  ! 1 wb%n < 1
  ! 2: wb%n /= size(a,1)
  ! 3: size(a) /= size(c)
  subroutine z_wb_times_general(wb,a,c,error)
    type(z_wb), intent(in) :: wb
    complex(kind=dp), target, dimension(:,:), intent(in) :: a
    complex(kind=dp), target, dimension(:,:), intent(out) :: c    
    type(error_info), intent(inout), optional :: error
    type(routine_info), parameter :: info=info_z_wb_times_general

    if (failure(error)) return
    call push_id(info, error)
    if (get_n(wb) < 1) then
       call set_error(1, info, error); return
    end if
    if (get_n(wb) /= size(a,1)) then
       call set_error(2, info, error); return
    end if
    if (size(a,1) /= size(c,1) .or. size(a,2) /= size(c,2)) then
       call set_error(3, info, error); return
    end if
    c=0.0_dp
    call f_z_wb_times_general_plus(wb%bc, get_n(wb), wb%lbw, wb%ubw, get_lbwmax(wb), &
         get_ubwmax(wb), wb%numrotsw, wb%jsw, wb%csw, wb%ssw, a, size(a,2), c)
    call pop_id(error)
  end subroutine z_wb_times_general
  

  subroutine f_z_wb_times_general_plus(bc, n, lbw, ubw, lbwmax, ubwmax, numrotsw, &
       jsw, csw, ssw, a, na, c)
    complex(kind=dp), target, dimension(n,na), intent(in) :: a
    complex(kind=dp), target, dimension(n,na), intent(out) :: c    
    integer(kind=int32), dimension(lbwmax,n), intent(in) :: jsw
    real(kind=dp), dimension(lbwmax,n), intent(in) :: csw
    complex(kind=dp), dimension(lbwmax,n), intent(in) :: ssw
    complex(kind=dp), dimension(lbwmax+ubwmax+1,n), intent(in) :: bc
    integer(kind=int32), dimension(n), intent(in) :: numrotsw
    integer(kind=int32), intent(in) :: ubw, lbw, n, lbwmax, ubwmax, na
    !
    integer(kind=int32) :: j,k,l,lbwk,ubwk,d0
    complex(kind=dp), dimension(n) :: x
    type(z_rotation) :: rot
    d0=ubw+1
    do j=1,na
       x=(0.0_dp,0.0_dp)
       do k=1,n-2
          ubwk=min(ubw,k-1) ! number of superdiagonals in column k
          lbwk=min(lbw,n-k) ! number of subdiagonals in column k
          x(k-ubwk:k+lbwk)=x(k-ubwk:k+lbwk)+bc(d0-ubwk:d0+lbwk,k)*a(k,j)
          ! apply w_{k} to x
          do l=1,numrotsw(k)
             rot%cosine=csw(l,k); rot%sine=ssw(l,k)
             call rotation_times_general(rot,x,jsw(l,k),jsw(l,k)+1)
          end do
       end do
       do k=n-1,n
          ubwk=min(ubw,k-1) ! number of superdiagonals in column k
          lbwk=min(lbw,n-k) ! number of subdiagonals in column k
          x(k-ubwk:k+lbwk)=x(k-ubwk:k+lbwk)+bc(d0-ubwk:d0+lbwk,k)*a(k,j)
       end do
       c(:,j)=c(:,j)+x
    end do
  end subroutine f_z_wb_times_general_plus
  
  ! BT times general

  function d_product_of_bt_and_general0(bt, a) result(c)
    real(kind=dp), dimension(:,:), intent(in) :: a
    real(kind=dp), dimension(:,:), allocatable :: c
    type(d_bt), intent(in) :: bt
    
    allocate(c(size(a,1),size(a,2)))
    call d_bt_times_general(bt,a,c)
  end function d_product_of_bt_and_general0

  function d_product_of_bt_and_general(bt, a, error) result(c)
    real(kind=dp), dimension(:,:), intent(in) :: a
    real(kind=dp), dimension(:,:), allocatable :: c
    type(d_bt), intent(in) :: bt
    type(error_info), intent(inout), optional :: error
    type(routine_info), parameter :: info=info_d_product_of_bt_and_general

    if (failure(error)) return
    call push_id(info, error)
    
    allocate(c(size(a,1),size(a,2)))
    call d_bt_times_general(bt,a,c,error)
    call pop_id(error)
  end function d_product_of_bt_and_general

  ! Errors
  ! 0: no error
  ! 1 bt%n < 1
  ! 2: bt%n /= size(a,1)
  ! 3: size(a) /= size(c)
  subroutine d_bt_times_general(bt,a,c,error)
    type(d_bt), intent(in) :: bt
    real(kind=dp), target, dimension(:,:), intent(in) :: a
    real(kind=dp), target, dimension(:,:), intent(out) :: c    
    type(error_info), intent(inout), optional :: error
    type(routine_info), parameter :: info=info_d_bt_times_general

    if (failure(error)) return
    call push_id(info, error)
    if (get_n(bt) < 1) then
       call set_error(1, info, error); return
    end if
    if (get_n(bt) /= size(a,1)) then
       call set_error(2, info, error); return
    end if
    if (size(a,1) /= size(c,1) .or. size(a,2) /= size(c,2)) then
       call set_error(3, info, error); return
    end if
    c=0.0_dp
    call f_d_bt_times_general_plus(bt%br, get_n(bt), bt%lbw, bt%ubw, get_lbwmax(bt), &
         get_ubwmax(bt), bt%numrotst, bt%kst, bt%cst, bt%sst, a, size(a,2), c)
    call pop_id(error)
  end subroutine d_bt_times_general
  

  subroutine f_d_bt_times_general_plus(br, n, lbw, ubw, lbwmax, ubwmax, numrotst, &
       kst, cst, sst, a, na, c)
    real(kind=dp), target, dimension(n,na), intent(in) :: a
    real(kind=dp), target, dimension(n,na), intent(out) :: c    
    integer(kind=int32), dimension(n,lbwmax), intent(in) :: kst
    real(kind=dp), dimension(n,lbwmax), intent(in) :: cst, sst
    real(kind=dp), dimension(n,lbwmax+ubwmax+1), intent(in) :: br
    integer(kind=int32), dimension(n), intent(in) :: numrotst
    integer(kind=int32), intent(in) :: ubw, lbw, n, lbwmax, ubwmax, na
    !
    integer(kind=int32) :: j,k,l,lbwj,ubwj,d0
    real(kind=dp), dimension(n) :: x
    type(d_rotation) :: rot
    d0=lbw+1
    do k=1,na
       x=a(:,k)
       ubwj=min(ubw,n-1) ! number of superdiagonals in row 1
       c(1,k)=c(1,k)+dot_product(br(1,d0:d0+ubwj),x(1:1+ubwj))
       lbwj=min(lbw,1) ! number of subdiagonals in row 2
       ubwj=min(ubw,n-2) ! number of superdiagonals in row 2
       c(2,k)=c(2,k)+dot_product(br(2,d0-lbwj:d0+ubwj),x(2-lbwj:2+ubwj))
       do j=3,n
          ! apply T_{j-1}^T to x
          do l=numrotst(j-1),1,-1
             rot%cosine=cst(j-1,l); rot%sine=sst(j-1,l)
             call rotation_times_general(trp_rot(rot),x,kst(j-1,l),kst(j-1,l)+1)
          end do
          lbwj=min(lbw,j-1) ! number of subdiagonals in row j
          ubwj=min(ubw,n-j)
          c(j,k)=c(j,k)+dot_product(br(j,d0-lbwj:d0+ubwj),x(j-lbwj:j+ubwj))
       end do
    end do
  end subroutine f_d_bt_times_general_plus

  function z_product_of_bt_and_general0(bt, a) result(c)
    complex(kind=dp), dimension(:,:), intent(in) :: a
    complex(kind=dp), dimension(:,:), allocatable :: c
    type(z_bt), intent(in) :: bt
    
    allocate(c(size(a,1),size(a,2)))
    call z_bt_times_general(bt,a,c)
  end function z_product_of_bt_and_general0

  function z_product_of_bt_and_general(bt, a, error) result(c)
    complex(kind=dp), dimension(:,:), intent(in) :: a
    complex(kind=dp), dimension(:,:), allocatable :: c
    type(z_bt), intent(in) :: bt
    type(error_info), intent(inout), optional :: error
    type(routine_info), parameter :: info=info_z_product_of_bt_and_general

    if (failure(error)) return
    call push_id(info, error)
    
    allocate(c(size(a,1),size(a,2)))
    call z_bt_times_general(bt,a,c,error)
    call pop_id(error)
  end function z_product_of_bt_and_general

  ! Errors
  ! 0: no error
  ! 1 bt%n < 1
  ! 2: bt%n /= size(a,1)
  ! 3: size(a) /= size(c)
  subroutine z_bt_times_general(bt,a,c,error)
    type(z_bt), intent(in) :: bt
    complex(kind=dp), target, dimension(:,:), intent(in) :: a
    complex(kind=dp), target, dimension(:,:), intent(out) :: c    
    type(error_info), intent(inout), optional :: error
    type(routine_info), parameter :: info=info_z_bt_times_general

    if (failure(error)) return
    call push_id(info, error)
    if (get_n(bt) < 1) then
       call set_error(1, info, error); return
    end if
    if (get_n(bt) /= size(a,1)) then
       call set_error(2, info, error); return
    end if
    if (size(a,1) /= size(c,1) .or. size(a,2) /= size(c,2)) then
       call set_error(3, info, error); return
    end if
    c=(0.0_dp,0.0_dp)
    call f_z_bt_times_general_plus(bt%br, get_n(bt), bt%lbw, bt%ubw, get_lbwmax(bt), &
         get_ubwmax(bt), bt%numrotst, bt%kst, bt%cst, bt%sst, a, size(a,2), c)
    call pop_id(error)
  end subroutine z_bt_times_general
  

  subroutine f_z_bt_times_general_plus(br, n, lbw, ubw, lbwmax, ubwmax, numrotst, &
       kst, cst, sst, a, na, c)
    complex(kind=dp), target, dimension(n,na), intent(in) :: a
    complex(kind=dp), target, dimension(n,na), intent(out) :: c    
    integer(kind=int32), dimension(n,lbwmax), intent(in) :: kst
    real(kind=dp), dimension(n,lbwmax), intent(in) :: cst
    complex(kind=dp), dimension(n,lbwmax), intent(in) :: sst
    complex(kind=dp), dimension(n,lbwmax+ubwmax+1), intent(in) :: br
    integer(kind=int32), dimension(n), intent(in) :: numrotst
    integer(kind=int32), intent(in) :: ubw, lbw, n, lbwmax, ubwmax, na
    !
    integer(kind=int32) :: j,k,l,lbwj,ubwj,d0
    complex(kind=dp), dimension(n) :: x
    type(z_rotation) :: rot
    d0=lbw+1
    do k=1,na
       x=a(:,k)
       ubwj=min(ubw,n-1) ! number of superdiagonals in row 1
       c(1,k)=c(1,k)+dot_no_conjg(br(1,d0:d0+ubwj),x(1:1+ubwj))
       lbwj=min(lbw,1) ! number of subdiagonals in row 2
       ubwj=min(ubw,n-2) ! number of superdiagonals in row 2
       c(2,k)=c(2,k)+dot_no_conjg(br(2,d0-lbwj:d0+ubwj),x(2-lbwj:2+ubwj))
       do j=3,n
          ! apply T_{j-1}^T to x
          do l=numrotst(j-1),1,-1
             rot%cosine=cst(j-1,l); rot%sine=sst(j-1,l)
             call rotation_times_general(trp_rot(rot),x,kst(j-1,l),kst(j-1,l)+1)
          end do
          lbwj=min(lbw,j-1) ! number of subdiagonals in row j
          ubwj=min(ubw,n-j)
          c(j,k)=c(j,k)+dot_no_conjg(br(j,d0-lbwj:d0+ubwj),x(j-lbwj:j+ubwj))
       end do
    end do
  end subroutine f_z_bt_times_general_plus
  
end module mod_products
