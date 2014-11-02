module mod_sweeps
  use mod_prec
  use mod_band_types
  use mod_rotation
  implicit none

  private

  public :: d_sweeps, c_sweeps

  public :: d_new_sweeps, c_new_sweeps

  public :: get_maxind, d_get_maxind, c_get_maxind, &
       get_minind, d_get_minind, c_get_minind

  public :: get_maxord, d_get_maxord, c_get_maxord

  public :: get_n, d_get_n_sweeps, c_get_n_sweeps

  public :: deallocate_sweeps, d_deallocate_sweeps, c_deallocate_sweeps

  public :: sweeps_times_general, d_sweeps_times_general, c_sweeps_times_general, &
       d_v_sweeps_times_general, c_v_sweeps_times_general

  public :: trp_sweeps_times_general, d_trp_sweeps_times_general, c_trp_sweeps_times_general, &
       d_v_trp_sweeps_times_general, c_v_trp_sweeps_times_general

  public :: general_times_sweeps, d_general_times_sweeps, c_general_times_sweeps, &
       d_v_general_times_sweeps, c_v_general_times_sweeps

  public :: general_times_trp_sweeps, d_general_times_trp_sweeps, c_general_times_trp_sweeps, &
       d_v_general_times_trp_sweeps, c_v_general_times_trp_sweeps

  ! These types represent a sequence of sweeps of plane rotations.
  ! In particular
  !
  ! Q = Q_{left} Q_{left+inc}... Q_{right}
  !
  ! where
  !
  ! Q_k = G_{numrots(k),k} ... G_{1,k}
  !
  ! is formed from rotations G_{j,k} with cosine and sine cs(j,k) and ss(j,k)
  ! acting on adjacent rows js(j,k) and js(j,k)+1.

  ! The fields ord and type are used to specifically describe
  ! structure of leading and trailing sweeps.  type takes on values
  ! 'L', 'T', or 'N'.

  type d_sweeps
     integer(kind=int32), private :: minind, maxind, n, maxord
     integer(kind=int32) :: left, right, inc, ord
     character :: type
     real(kind=dp), allocatable, dimension(:,:) :: cs, ss
     integer(kind=int32), allocatable, dimension(:) :: numrots
     integer(kind=int32), allocatable, dimension(:,:) :: js
  end type d_sweeps

  type c_sweeps
     integer(kind=int32), private :: minind, maxind, n, maxord
     integer(kind=int32) :: left, right, inc, ord
     character :: type
     complex(kind=dp), allocatable, dimension(:,:) :: cs, ss
     integer(kind=int32), allocatable, dimension(:) :: numrots
     integer(kind=int32), allocatable, dimension(:,:) :: js
  end type c_sweeps

  interface get_maxind
     module procedure d_get_maxind, c_get_maxind
  end interface get_maxind

  interface get_minind
     module procedure d_get_minind, c_get_minind
  end interface get_minind

  interface get_maxord
     module procedure d_get_maxord, c_get_maxord
  end interface get_maxord

  interface get_n
     module procedure d_get_n_sweeps, c_get_n_sweeps
  end interface get_n

  interface deallocate_sweeps
     module procedure d_deallocate_sweeps, c_deallocate_sweeps
  end interface deallocate_sweeps

  interface sweeps_times_general
     module procedure d_sweeps_times_general, c_sweeps_times_general, &
          d_v_sweeps_times_general, c_v_sweeps_times_general
  end interface sweeps_times_general

  interface trp_sweeps_times_general
     module procedure d_trp_sweeps_times_general, c_trp_sweeps_times_general, &
          d_v_trp_sweeps_times_general, c_v_trp_sweeps_times_general
  end interface trp_sweeps_times_general

  interface general_times_sweeps
     module procedure d_general_times_sweeps, c_general_times_sweeps, &
          d_v_general_times_sweeps, c_v_general_times_sweeps
  end interface general_times_sweeps

  interface general_times_trp_sweeps
     module procedure d_general_times_trp_sweeps, c_general_times_trp_sweeps, &
          d_v_general_times_trp_sweeps, c_v_general_times_trp_sweeps
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

  integer(kind=int32) function c_get_maxord(sw) result(maxord)
    type(c_sweeps) :: sw
    maxord=sw%maxord
  end function c_get_maxord

  integer(kind=int32) function c_get_maxind(sw) result(maxind)
    type(c_sweeps) :: sw
    maxind=sw%maxind
  end function c_get_maxind

  integer(kind=int32) function c_get_minind(sw) result(minind)
    type(c_sweeps) :: sw
    minind=sw%minind
  end function c_get_minind

  integer(kind=int32) function d_get_n_sweeps(sw) result(n)
    type(d_sweeps) :: sw
    n=sw%n
  end function d_get_n_sweeps

  integer(kind=int32) function c_get_n_sweeps(sw) result(n)
    type(c_sweeps) :: sw
    n=sw%n
  end function c_get_n_sweeps

  type(d_sweeps) function d_new_sweeps(n,minind, maxind, maxord) result(sw)
    integer(kind=int32), intent(in) :: n, minind, maxind, maxord
    sw%n=n; sw%maxind=maxind; sw%minind=minind; sw%maxord=maxord
    allocate(sw%cs(maxord,minind:maxind), sw%ss(maxord,minind:maxind), &
         sw%js(maxord,minind:maxind), sw%numrots(minind:maxind))
    sw%left=0; sw%right=-1; sw%inc=1; sw%ord=-1
    sw%type='N'
    sw%cs=0.0_dp; sw%ss=0.0_dp
  end function d_new_sweeps

  type(c_sweeps) function c_new_sweeps(n,minind, maxind, maxord) result(sw)
    integer(kind=int32), intent(in) :: n, minind, maxind, maxord
    sw%n=n; sw%maxind=maxind; sw%minind=minind; sw%maxord=maxord
    allocate(sw%cs(maxord,minind:maxind), sw%ss(maxord,minind:maxind), &
         sw%js(maxord,minind:maxind), sw%numrots(minind:maxind))
    sw%left=0; sw%right=-1; sw%inc=1; sw%ord=-1
    sw%type='N'
    sw%cs=(0.0_dp,0.0_dp); sw%ss=(0.0_dp,0.0_dp)
  end function c_new_sweeps

  subroutine d_deallocate_sweeps(sw)
    type(d_sweeps), intent(inout) :: sw
    deallocate(sw%cs, sw%ss, sw%js, sw%numrots)
  end subroutine d_deallocate_sweeps

  subroutine c_deallocate_sweeps(sw)
    type(c_sweeps), intent(inout) :: sw
    deallocate(sw%cs, sw%ss, sw%js, sw%numrots)
  end subroutine c_deallocate_sweeps

  subroutine d_sweeps_times_general(sw,a)
    type(d_sweeps) :: sw
    real(kind=dp), dimension(:,:), intent(inout) :: a
    type(d_rotation) :: rot
    integer(kind=int32) :: j, m, l,jj
    m=size(a,1)
    do l=sw%right,sw%left,-sw%inc
       do j=1,sw%numrots(l)
          jj=sw%js(j,l)
          rot%cosine=sw%cs(j,l); rot%sine=sw%ss(j,l)
          call d_rotation_times_general(rot,a,jj,jj+1)
       end do
    end do
  end subroutine d_sweeps_times_general

  subroutine c_sweeps_times_general(sw,a)
    type(c_sweeps) :: sw
    complex(kind=dp), dimension(:,:), intent(inout) :: a
    type(c_rotation) :: rot
    integer(kind=int32) :: j, m, l,jj
    m=size(a,1)
    do l=sw%right,sw%left,-sw%inc
       do j=1,sw%numrots(l)
          jj=sw%js(j,l)
          rot%cosine=sw%cs(j,l); rot%sine=sw%ss(j,l)
          call c_rotation_times_general(rot,a,jj,jj+1)
       end do
    end do
  end subroutine c_sweeps_times_general

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

  subroutine c_v_sweeps_times_general(sw,a)
    type(c_sweeps) :: sw
    complex(kind=dp), dimension(:), intent(inout) :: a
    integer(kind=int32) :: j, l, jj
    complex(kind=dp) :: tmp, c, s 
    do l=sw%right,sw%left,-sw%inc
       do j=1,sw%numrots(l)
          jj=sw%js(j,l)
          c=sw%cs(j,l); s=sw%ss(j,l)
          tmp=a(jj)
          a(jj)=c*tmp-conjg(s)*a(jj+1)
          a(jj+1)=s*tmp+c*a(jj+1)
       end do
    end do
  end subroutine c_v_sweeps_times_general

  subroutine d_trp_sweeps_times_general(sw,a)
    type(d_sweeps) :: sw
    real(kind=dp), dimension(:,:), intent(inout) :: a
    type(d_rotation) :: rot
    integer(kind=int32) :: j, l, jj
    do l=sw%left,sw%right,sw%inc
       do j=sw%numrots(l),1,-1
          jj=sw%js(j,l)
          rot%cosine=sw%cs(j,l); rot%sine=sw%ss(j,l)
          call d_rotation_times_general(trp_rot(rot),a,jj,jj+1)
       end do
    end do
  end subroutine d_trp_sweeps_times_general

  subroutine c_trp_sweeps_times_general(sw,a)
    type(c_sweeps) :: sw
    complex(kind=dp), dimension(:,:), intent(inout) :: a
    type(c_rotation) :: rot
    integer(kind=int32) :: j, l, jj
    do l=sw%left,sw%right,sw%inc
       do j=sw%numrots(l),1,-1
          jj=sw%js(j,l)
          rot%cosine=sw%cs(j,l); rot%sine=sw%ss(j,l)
          call c_rotation_times_general(trp_rot(rot),a,jj,jj+1)
       end do
    end do
  end subroutine c_trp_sweeps_times_general

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

  subroutine c_v_trp_sweeps_times_general(sw,a)
    type(c_sweeps) :: sw
    complex(kind=dp), dimension(:), intent(inout) :: a
    integer(kind=int32) :: j, l, jj
    complex(kind=dp) :: tmp, c, s 
    do l=sw%left,sw%right,sw%inc
       do j=sw%numrots(l),1,-1
          jj=sw%js(j,l)
          c=sw%cs(j,l); s=sw%ss(j,l)
          tmp=a(jj)
          a(jj)=c*tmp+conjg(s)*a(jj+1)
          a(jj+1)=-s*tmp+c*a(jj+1)
       end do
    end do
  end subroutine c_v_trp_sweeps_times_general

  subroutine d_general_times_sweeps(a,sw)
    type(d_sweeps) :: sw
    real(kind=dp), dimension(:,:), intent(inout) :: a
    type(d_rotation) :: rot
    integer(kind=int32) :: j, l, kk
    do l=sw%left,sw%right,sw%inc
       do j=sw%numrots(l),1,-1
          kk=sw%js(j,l)
          rot%cosine=sw%cs(j,l); rot%sine=sw%ss(j,l)
          call d_general_times_rotation(a,rot,kk,kk+1)
       end do
    end do
  end subroutine d_general_times_sweeps

  subroutine c_general_times_sweeps(a,sw)
    type(c_sweeps) :: sw
    complex(kind=dp), dimension(:,:), intent(inout) :: a
    type(c_rotation) :: rot
    integer(kind=int32) :: j, l, kk
    do l=sw%left,sw%right,sw%inc
       do j=sw%numrots(l),1,-1
          kk=sw%js(j,l)
          rot%cosine=sw%cs(j,l); rot%sine=sw%ss(j,l)
          call c_general_times_rotation(a,rot,kk,kk+1)
       end do
    end do
  end subroutine c_general_times_sweeps

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

  subroutine c_v_general_times_sweeps(a,sw)
    type(c_sweeps) :: sw
    complex(kind=dp), dimension(:), intent(inout) :: a
    integer(kind=int32) :: j, l, kk
    complex(kind=dp) :: c,s,tmp
    do l=sw%left,sw%right,sw%inc
       do j=sw%numrots(l),1,-1
          kk=sw%js(j,l)
          c=sw%cs(j,l); s=sw%ss(j,l)
          tmp=a(kk)
          a(kk)=c*tmp+s*a(kk+1)
          a(kk+1)=-conjg(s)*tmp+c*a(kk+1)
       end do
    end do
  end subroutine c_v_general_times_sweeps

  subroutine d_general_times_trp_sweeps(a,sw)
    type(d_sweeps) :: sw
    real(kind=dp), dimension(:,:), intent(inout) :: a
    type(d_rotation) :: rot
    integer(kind=int32) :: j, l, kk
    do l=sw%right,sw%left,-sw%inc
       do j=1,sw%numrots(l)
          rot%cosine=sw%cs(j,l); rot%sine=sw%ss(j,l)
          kk=sw%js(j,l)
          call d_general_times_rotation(a, trp_rot(rot),kk,kk+1)
       end do
    end do
  end subroutine d_general_times_trp_sweeps

  subroutine c_general_times_trp_sweeps(a,sw)
    type(c_sweeps) :: sw
    complex(kind=dp), dimension(:,:), intent(inout) :: a
    type(c_rotation) :: rot
    integer(kind=int32) :: j, l, kk
    do l=sw%right,sw%left,-sw%inc
       do j=1,sw%numrots(l)
          rot%cosine=sw%cs(j,l); rot%sine=sw%ss(j,l)
          kk=sw%js(j,l)
          call c_general_times_rotation(a, trp_rot(rot),kk,kk+1)
       end do
    end do
  end subroutine c_general_times_trp_sweeps

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

  subroutine c_v_general_times_trp_sweeps(a,sw)
    type(c_sweeps) :: sw
    complex(kind=dp), dimension(:), intent(inout) :: a
    integer(kind=int32) :: j, l, kk
    complex(kind=dp) :: c,s,tmp
    do l=sw%right,sw%left,-sw%inc
       do j=1,sw%numrots(l)
          c=sw%cs(j,l); s=sw%ss(j,l)
          kk=sw%js(j,l)
          tmp=a(kk)
          a(kk)=c*tmp-s*a(kk+1)
          a(kk+1)=conjg(s)*tmp+c*a(kk+1)
       end do
    end do
  end subroutine c_v_general_times_trp_sweeps

end module mod_sweeps
