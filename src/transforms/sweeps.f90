module mod_sweeps
  use mod_prec
  use mod_band_types
  use mod_rotation
  use mod_utility
  implicit none

  ! Derived types for sequences of sweeps of plane rotations that act
  ! on adjacent rows.  The format is flexible with separate arrays for
  ! the cosines, sines, and the rows on which the rotations act.  See
  ! the module for more detail.  The primary use is in representing
  ! sequences of order r leading or trailing transformations used in
  ! the computation of a row compression or QR factorization.  There
  ! are also routines for applying such transformations to
  ! unstructured matrices (with or without transposition) and routines
  ! for generating random sweeps (used in test code).

  private

  public :: d_sweeps, z_sweeps

  public :: d_new_sweeps, z_new_sweeps

  public :: d_random_full_sweeps, z_random_full_sweeps

  public :: get_maxind, d_get_maxind, z_get_maxind, &
       get_minind, d_get_minind, z_get_minind

  public :: d_sweeps_general_product_of, z_sweeps_general_product_of, &
       d_general_sweeps_product_of, z_general_sweeps_product_of, &
       d_v_sweeps_general_product_of, z_v_sweeps_general_product_of, &
       d_v_general_sweeps_product_of, z_v_general_sweeps_product_of, operator(*)

  public :: trp, d_trp_sweeps, z_trp_sweeps

  public :: get_maxord, d_get_maxord, z_get_maxord

  public :: get_n, d_get_n_sweeps, z_get_n_sweeps

  public :: deallocate_sweeps, d_deallocate_sweeps, z_deallocate_sweeps

  public :: sweeps_times_general, d_sweeps_times_general, z_sweeps_times_general, &
       d_v_sweeps_times_general, z_v_sweeps_times_general

  public :: trp_sweeps_times_general, d_trp_sweeps_times_general, &
       z_trp_sweeps_times_general, &
       d_v_trp_sweeps_times_general, z_v_trp_sweeps_times_general

  public :: general_times_sweeps, d_general_times_sweeps, z_general_times_sweeps, &
       d_v_general_times_sweeps, z_v_general_times_sweeps

  public :: general_times_trp_sweeps, d_general_times_trp_sweeps, &
       z_general_times_trp_sweeps, &
       d_v_general_times_trp_sweeps, z_v_general_times_trp_sweeps

  ! The types d_sweeps and z_sweeps represent a sequence of sweeps of
  ! plane rotations.  In particular
  !
  ! Q = Q_{left} Q_{left+inc}... Q_{right}
  !
  ! where
  !
  ! Q_k = G_{numrots(k),k} ... G_{1,k}
  !
  ! is formed from rotations G_{j,k} with cosine and sine cs(j,k) and ss(j,k)
  ! acting on adjacent rows js(j,k) and js(j,k)+1.

  ! The range left..right can be increasing or decreasing, depending on inc.
  ! cs, ss, and js have bounds (maxord,minind:maxind).  So both left and
  ! right need to be in the interval [minind,maxind].

  ! The fields ord and type are used to specifically describe
  ! structure of leading and trailing sweeps.  type takes on values;
  ! -1=nothing, 0=leading, 1=trailing.  ord=-1 indicates the order is not specified.

  ! It is expected that: cs(maxord,minind:maxind)

  type d_sweeps
     integer(kind=int32), private :: minind, maxind, n, maxord
     integer(kind=int32) :: left, right, inc, ord, type
     real(kind=dp), allocatable, dimension(:,:) :: cs, ss
     integer(kind=int32), allocatable, dimension(:) :: numrots
     integer(kind=int32), allocatable, dimension(:,:) :: js
     logical :: trp
  end type d_sweeps

  type z_sweeps
     integer(kind=int32), private :: minind, maxind, n, maxord
     integer(kind=int32) :: left, right, inc, ord, type
     real(kind=dp), allocatable, dimension(:,:) :: cs
     complex(kind=dp), allocatable, dimension(:,:) :: ss
     integer(kind=int32), allocatable, dimension(:) :: numrots
     integer(kind=int32), allocatable, dimension(:,:) :: js
     logical :: trp
  end type z_sweeps

  interface get_maxind
     module procedure d_get_maxind, z_get_maxind
  end interface get_maxind

  interface get_minind
     module procedure d_get_minind, z_get_minind
  end interface get_minind

  interface get_maxord
     module procedure d_get_maxord, z_get_maxord
  end interface get_maxord

  interface get_n
     module procedure d_get_n_sweeps, z_get_n_sweeps
  end interface get_n

  interface deallocate_sweeps
     module procedure d_deallocate_sweeps, z_deallocate_sweeps
  end interface deallocate_sweeps

  interface trp
     module procedure d_trp_sweeps, z_trp_sweeps
  end interface trp

  interface operator (*)
     module procedure d_sweeps_general_product_of, z_sweeps_general_product_of, &
          d_general_sweeps_product_of, z_general_sweeps_product_of, &
          d_v_sweeps_general_product_of, z_v_sweeps_general_product_of, &
          d_v_general_sweeps_product_of, z_v_general_sweeps_product_of
  end interface operator (*)

  interface sweeps_times_general
     module procedure d_sweeps_times_general, z_sweeps_times_general, &
          d_v_sweeps_times_general, z_v_sweeps_times_general
  end interface sweeps_times_general

  interface trp_sweeps_times_general
     module procedure d_trp_sweeps_times_general, z_trp_sweeps_times_general, &
          d_v_trp_sweeps_times_general, z_v_trp_sweeps_times_general
  end interface trp_sweeps_times_general

  interface general_times_sweeps
     module procedure d_general_times_sweeps, z_general_times_sweeps, &
          d_v_general_times_sweeps, z_v_general_times_sweeps
  end interface general_times_sweeps

  interface general_times_trp_sweeps
     module procedure d_general_times_trp_sweeps, z_general_times_trp_sweeps, &
          d_v_general_times_trp_sweeps, z_v_general_times_trp_sweeps
  end interface general_times_trp_sweeps

contains

  integer(kind=int32) function d_get_maxind(sw) result(maxind)
    type(d_sweeps) :: sw
    maxind=sw%maxind
  end function d_get_maxind

  integer(kind=int32) function d_get_minind(sw) result(minind)
    type(d_sweeps) :: sw
    minind=sw%minind
  end function d_get_minind

  integer(kind=int32) function d_get_maxord(sw) result(maxord)
    type(d_sweeps) :: sw
    maxord=sw%maxord
  end function d_get_maxord

  integer(kind=int32) function z_get_maxord(sw) result(maxord)
    type(z_sweeps) :: sw
    maxord=sw%maxord
  end function z_get_maxord

  integer(kind=int32) function z_get_maxind(sw) result(maxind)
    type(z_sweeps) :: sw
    maxind=sw%maxind
  end function z_get_maxind

  integer(kind=int32) function z_get_minind(sw) result(minind)
    type(z_sweeps) :: sw
    minind=sw%minind
  end function z_get_minind

  integer(kind=int32) function d_get_n_sweeps(sw) result(n)
    type(d_sweeps) :: sw
    n=sw%n
  end function d_get_n_sweeps

  integer(kind=int32) function z_get_n_sweeps(sw) result(n)
    type(z_sweeps) :: sw
    n=sw%n
  end function z_get_n_sweeps

  type(d_sweeps) function d_new_sweeps(n,minind, maxind, maxord) result(sw)
    integer(kind=int32), intent(in) :: n, minind, maxind, maxord
    sw%n=n; sw%maxind=maxind; sw%minind=minind; sw%maxord=maxord
    allocate(sw%cs(maxord,minind:maxind), sw%ss(maxord,minind:maxind), &
         sw%js(maxord,minind:maxind), sw%numrots(minind:maxind))
    sw%left=0; sw%right=-1; sw%inc=1; sw%ord=-1
    sw%type=-1
    sw%cs=0.0_dp; sw%ss=0.0_dp
    sw%trp=.false.
  end function d_new_sweeps

  type(z_sweeps) function z_new_sweeps(n,minind, maxind, maxord) result(sw)
    integer(kind=int32), intent(in) :: n, minind, maxind, maxord
    sw%n=n; sw%maxind=maxind; sw%minind=minind; sw%maxord=maxord
    allocate(sw%cs(maxord,minind:maxind), sw%ss(maxord,minind:maxind), &
         sw%js(maxord,minind:maxind), sw%numrots(minind:maxind))
    sw%left=0; sw%right=-1; sw%inc=1; sw%ord=-1
    sw%type=-1
    sw%cs=0.0_dp; sw%ss=(0.0_dp,0.0_dp)
    sw%trp=.false.
  end function z_new_sweeps

  subroutine d_deallocate_sweeps(sw)
    type(d_sweeps), intent(inout) :: sw
    deallocate(sw%cs, sw%ss, sw%js, sw%numrots)
  end subroutine d_deallocate_sweeps

  subroutine z_deallocate_sweeps(sw)
    type(z_sweeps), intent(inout) :: sw
    deallocate(sw%cs, sw%ss, sw%js, sw%numrots)
  end subroutine z_deallocate_sweeps

  type(d_sweeps) function d_random_full_sweeps(n,l) result(sw)
    integer(kind=int32) :: n, l, j, k
    real(kind=dp) :: nrm
    sw=d_new_sweeps(n, 1, l, n-1)
    sw%left=1
    sw%right=l
    sw%numrots=n-1
    do j=1,n-1
       sw%js(j,:)=j
    end do
    call random_matrix_to(sw%cs)
    call random_matrix_to(sw%ss)
    do j=1,n-1
       do k=1,l
          nrm=sqrt(sw%cs(j,k)**2+sw%ss(j,k)**2)
          sw%cs(j,k)=sw%cs(j,k)/nrm
          sw%ss(j,k)=sw%ss(j,k)/nrm
       end do
    end do
  end function d_random_full_sweeps

  type(z_sweeps) function z_random_full_sweeps(n,l) result(sw)
    integer(kind=int32) :: n, l, j, k
    real(kind=dp) :: nrm
    sw=z_new_sweeps(n, 1, l, n-1)
    sw%left=1
    sw%right=l
    sw%numrots=n-1
    do j=1,n-1
       sw%js(j,:)=j
    end do
    call random_matrix_to(sw%cs)
    call random_matrix_to(sw%ss)
    do j=1,n-1
       do k=1,l
          nrm=sqrt(sw%cs(j,k)**2+abs(sw%ss(j,k))**2)
          sw%cs(j,k)=sw%cs(j,k)/nrm
          sw%ss(j,k)=sw%ss(j,k)/nrm
       end do
    end do
  end function z_random_full_sweeps

  ! Multiply routines

  function d_sweeps_general_product_of(sw,a) result(b)
    real(kind=dp), allocatable, dimension(:,:) :: b
    real(kind=dp), dimension(:,:), intent(in) :: a
    type(d_sweeps), intent(in) :: sw

    b=a
    if (sw%trp) then
       call d_trp_sweeps_times_general(sw,b)
    else
       call d_sweeps_times_general(sw,b)
    end if
  end function d_sweeps_general_product_of

  function z_sweeps_general_product_of(sw,a) result(b)
    complex(kind=dp), allocatable, dimension(:,:) :: b
    complex(kind=dp), dimension(:,:), intent(in) :: a
    type(z_sweeps), intent(in) :: sw

    b=a
    if (sw%trp) then
       call z_trp_sweeps_times_general(sw,b)
    else
       call z_sweeps_times_general(sw,b)
    end if
  end function z_sweeps_general_product_of

  function d_general_sweeps_product_of(a,sw) result(b)
    real(kind=dp), allocatable, dimension(:,:) :: b
    real(kind=dp), dimension(:,:), intent(in) :: a
    type(d_sweeps), intent(in) :: sw

    b=a
    if (sw%trp) then
       call d_general_times_trp_sweeps(b,sw)
    else
       call d_general_times_sweeps(b,sw)
    end if
  end function d_general_sweeps_product_of

  function z_general_sweeps_product_of(a,sw) result(b)
    complex(kind=dp), allocatable, dimension(:,:) :: b
    complex(kind=dp), dimension(:,:), intent(in) :: a
    type(z_sweeps), intent(in) :: sw

    b=a
    if (sw%trp) then
       call z_general_times_trp_sweeps(b,sw)
    else
       call z_general_times_sweeps(b,sw)
    end if
  end function z_general_sweeps_product_of

  ! vectors

  function d_v_sweeps_general_product_of(sw,a) result(b)
    real(kind=dp), allocatable, dimension(:) :: b
    real(kind=dp), dimension(:), intent(in) :: a
    type(d_sweeps), intent(in) :: sw

    b=a
    if (sw%trp) then
       call d_v_trp_sweeps_times_general(sw,b)
    else
       call d_v_sweeps_times_general(sw,b)
    end if
  end function d_v_sweeps_general_product_of

  function z_v_sweeps_general_product_of(sw,a) result(b)
    complex(kind=dp), allocatable, dimension(:) :: b
    complex(kind=dp), dimension(:), intent(in) :: a
    type(z_sweeps), intent(in) :: sw

    b=a
    if (sw%trp) then
       call z_v_trp_sweeps_times_general(sw,b)
    else
       call z_v_sweeps_times_general(sw,b)
    end if
  end function z_v_sweeps_general_product_of

  function d_v_general_sweeps_product_of(a,sw) result(b)
    real(kind=dp), allocatable, dimension(:) :: b
    real(kind=dp), dimension(:), intent(in) :: a
    type(d_sweeps), intent(in) :: sw

    b=a
    if (sw%trp) then
       call d_v_general_times_trp_sweeps(b,sw)
    else
       call d_v_general_times_sweeps(b,sw)
    end if
  end function d_v_general_sweeps_product_of

  function z_v_general_sweeps_product_of(a,sw) result(b)
    complex(kind=dp), allocatable, dimension(:) :: b
    complex(kind=dp), dimension(:), intent(in) :: a
    type(z_sweeps), intent(in) :: sw

    b=a
    if (sw%trp) then
       call z_v_general_times_trp_sweeps(b,sw)
    else
       call z_v_general_times_sweeps(b,sw)
    end if
  end function z_v_general_sweeps_product_of

  function d_trp_sweeps(swa) result(swb)
    type(d_sweeps), allocatable :: swb
    type(d_sweeps), intent(in) :: swa

    swb=swa
    swb%trp=.not. swb%trp
  end function d_trp_sweeps

  function z_trp_sweeps(swa) result(swb)
    type(z_sweeps), allocatable :: swb
    type(z_sweeps), intent(in) :: swa

    swb=swa
    swb%trp=.not. swb%trp
  end function z_trp_sweeps
  
  subroutine d_sweeps_times_general(sw,a)
    type(d_sweeps) :: sw
    real(kind=dp), dimension(:,:), intent(inout) :: a
    integer(kind=int32) :: j, m, l,jj
    m=size(a,1)
    do l=sw%right,sw%left,-sw%inc
       do j=1,sw%numrots(l)
          jj=sw%js(j,l)
          call f_d_rotation_times_general(sw%cs(j,l),sw%ss(j,l),a,jj,jj+1)
       end do
    end do
  end subroutine d_sweeps_times_general

  subroutine z_sweeps_times_general(sw,a)
    type(z_sweeps) :: sw
    complex(kind=dp), dimension(:,:), intent(inout) :: a
    integer(kind=int32) :: j, m, l,jj
    m=size(a,1)
    do l=sw%right,sw%left,-sw%inc
       do j=1,sw%numrots(l)
          jj=sw%js(j,l)
          call f_z_rotation_times_general(sw%cs(j,l),sw%ss(j,l),a,jj,jj+1)
       end do
    end do
  end subroutine z_sweeps_times_general

  subroutine d_v_sweeps_times_general(sw,a)
    type(d_sweeps) :: sw
    real(kind=dp), dimension(:), intent(inout) :: a
    integer(kind=int32) :: j, l, jj
    real(kind=dp) :: tmp, c, s 
    do l=sw%right,sw%left,-sw%inc
       do j=1,sw%numrots(l)
          jj=sw%js(j,l)
          c=sw%cs(j,l); s=sw%ss(j,l)
          tmp=a(jj)
          a(jj)=c*tmp-s*a(jj+1)
          a(jj+1)=s*tmp+c*a(jj+1)
       end do
    end do
  end subroutine d_v_sweeps_times_general

  subroutine z_v_sweeps_times_general(sw,a)
    type(z_sweeps) :: sw
    complex(kind=dp), dimension(:), intent(inout) :: a
    integer(kind=int32) :: j, l, jj
    real(kind=dp) :: c
    complex(kind=dp) :: tmp, s 

    do l=sw%right,sw%left,-sw%inc
       do j=1,sw%numrots(l)
          jj=sw%js(j,l)
          c=sw%cs(j,l); s=sw%ss(j,l)
          tmp=a(jj)
          a(jj)=c*tmp-conjg(s)*a(jj+1)
          a(jj+1)=s*tmp+c*a(jj+1)
       end do
    end do
  end subroutine z_v_sweeps_times_general

  subroutine d_trp_sweeps_times_general(sw,a)
    type(d_sweeps) :: sw
    real(kind=dp), dimension(:,:), intent(inout) :: a
    integer(kind=int32) :: j, l, jj
    do l=sw%left,sw%right,sw%inc
       do j=sw%numrots(l),1,-1
          jj=sw%js(j,l)
          call f_d_rotation_times_general(sw%cs(j,l),-sw%ss(j,l),a,jj,jj+1)
       end do
    end do
  end subroutine d_trp_sweeps_times_general

  subroutine z_trp_sweeps_times_general(sw,a)
    type(z_sweeps) :: sw
    complex(kind=dp), dimension(:,:), intent(inout) :: a
    integer(kind=int32) :: j, l, jj
    do l=sw%left,sw%right,sw%inc
       do j=sw%numrots(l),1,-1
          jj=sw%js(j,l)
          call f_z_rotation_times_general(sw%cs(j,l),-sw%ss(j,l),a,jj,jj+1)          
       end do
    end do
  end subroutine z_trp_sweeps_times_general

  subroutine d_v_trp_sweeps_times_general(sw,a)
    type(d_sweeps) :: sw
    real(kind=dp), dimension(:), intent(inout) :: a
    integer(kind=int32) :: j, l, jj
    real(kind=dp) :: tmp, c, s 
    do l=sw%left,sw%right,sw%inc
       do j=sw%numrots(l),1,-1
          jj=sw%js(j,l)
          c=sw%cs(j,l); s=sw%ss(j,l)
          tmp=a(jj)
          a(jj)=c*tmp+s*a(jj+1)
          a(jj+1)=-s*tmp+c*a(jj+1)
       end do
    end do
  end subroutine d_v_trp_sweeps_times_general

  subroutine z_v_trp_sweeps_times_general(sw,a)
    type(z_sweeps) :: sw
    complex(kind=dp), dimension(:), intent(inout) :: a
    integer(kind=int32) :: j, l, jj
    real(kind=dp) :: c
    complex(kind=dp) :: tmp, s
    do l=sw%left,sw%right,sw%inc
       do j=sw%numrots(l),1,-1
          jj=sw%js(j,l)
          c=sw%cs(j,l); s=sw%ss(j,l)
          tmp=a(jj)
          a(jj)=c*tmp+conjg(s)*a(jj+1)
          a(jj+1)=-s*tmp+c*a(jj+1)
       end do
    end do
  end subroutine z_v_trp_sweeps_times_general

  subroutine d_general_times_sweeps(a,sw)
    type(d_sweeps) :: sw
    real(kind=dp), dimension(:,:), intent(inout) :: a
    integer(kind=int32) :: j, l, kk
    do l=sw%left,sw%right,sw%inc
       do j=sw%numrots(l),1,-1
          kk=sw%js(j,l)
          call f_d_general_times_rotation(a,sw%cs(j,l),sw%ss(j,l),kk,kk+1)
       end do
    end do
  end subroutine d_general_times_sweeps

  subroutine z_general_times_sweeps(a,sw)
    type(z_sweeps) :: sw
    complex(kind=dp), dimension(:,:), intent(inout) :: a
    integer(kind=int32) :: j, l, kk
    do l=sw%left,sw%right,sw%inc
       do j=sw%numrots(l),1,-1
          kk=sw%js(j,l)
          call f_z_general_times_rotation(a,sw%cs(j,l),sw%ss(j,l),kk,kk+1)
       end do
    end do
  end subroutine z_general_times_sweeps

  subroutine d_v_general_times_sweeps(a,sw)
    type(d_sweeps) :: sw
    real(kind=dp), dimension(:), intent(inout) :: a
    integer(kind=int32) :: j, l, kk
    real(kind=dp) :: c,s,tmp
    do l=sw%left,sw%right,sw%inc
       do j=sw%numrots(l),1,-1
          kk=sw%js(j,l)
          c=sw%cs(j,l); s=sw%ss(j,l)
          tmp=a(kk)
          a(kk)=c*tmp+s*a(kk+1)
          a(kk+1)=-s*tmp+c*a(kk+1)
       end do
    end do
  end subroutine d_v_general_times_sweeps

  subroutine z_v_general_times_sweeps(a,sw)
    type(z_sweeps) :: sw
    complex(kind=dp), dimension(:), intent(inout) :: a
    integer(kind=int32) :: j, l, kk
    complex(kind=dp) :: c
    complex(kind=dp) :: s,tmp
    do l=sw%left,sw%right,sw%inc
       do j=sw%numrots(l),1,-1
          kk=sw%js(j,l)
          c=sw%cs(j,l); s=sw%ss(j,l)
          tmp=a(kk)
          a(kk)=c*tmp+s*a(kk+1)
          a(kk+1)=-conjg(s)*tmp+c*a(kk+1)
       end do
    end do
  end subroutine z_v_general_times_sweeps

  subroutine d_general_times_trp_sweeps(a,sw)
    type(d_sweeps) :: sw
    real(kind=dp), dimension(:,:), intent(inout) :: a
    integer(kind=int32) :: j, l, kk
    do l=sw%right,sw%left,-sw%inc
       do j=1,sw%numrots(l)
          kk=sw%js(j,l)
          call f_d_general_times_rotation(a,sw%cs(j,l),-sw%ss(j,l),kk,kk+1)
       end do
    end do
  end subroutine d_general_times_trp_sweeps

  subroutine z_general_times_trp_sweeps(a,sw)
    type(z_sweeps) :: sw
    complex(kind=dp), dimension(:,:), intent(inout) :: a
    integer(kind=int32) :: j, l, kk
    do l=sw%right,sw%left,-sw%inc
       do j=1,sw%numrots(l)
          kk=sw%js(j,l)
          call f_z_general_times_rotation(a, sw%cs(j,l),-sw%ss(j,l),kk,kk+1)
       end do
    end do
  end subroutine z_general_times_trp_sweeps

  subroutine d_v_general_times_trp_sweeps(a,sw)
    type(d_sweeps) :: sw
    real(kind=dp), dimension(:), intent(inout) :: a
    integer(kind=int32) :: j, l, kk
    real(kind=dp) :: c,s,tmp
    do l=sw%right,sw%left,-sw%inc
       do j=1,sw%numrots(l)
          c=sw%cs(j,l); s=sw%ss(j,l)
          kk=sw%js(j,l)
          tmp=a(kk)
          a(kk)=c*tmp-s*a(kk+1)
          a(kk+1)=s*tmp+c*a(kk+1)
       end do
    end do
  end subroutine d_v_general_times_trp_sweeps

  subroutine z_v_general_times_trp_sweeps(a,sw)
    type(z_sweeps) :: sw
    complex(kind=dp), dimension(:), intent(inout) :: a
    integer(kind=int32) :: j, l, kk
    real(kind=dp) :: c
    complex(kind=dp) :: s,tmp
    do l=sw%right,sw%left,-sw%inc
       do j=1,sw%numrots(l)
          c=sw%cs(j,l); s=sw%ss(j,l)
          kk=sw%js(j,l)
          tmp=a(kk)
          a(kk)=c*tmp-s*a(kk+1)
          a(kk+1)=conjg(s)*tmp+c*a(kk+1)
       end do
    end do
  end subroutine z_v_general_times_trp_sweeps

end module mod_sweeps
