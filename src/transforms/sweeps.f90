module sweeps
  use misc
  use types
  use shift
  use conversion

  implicit none

  ! These types represent a linear transformation.
  ! Thus
  ! Q = Q_{num} Q_{num-1} ... Q_{1}
  ! where
  ! Q_k = G_{k,numrots(k)} ... G_{k,1}
  ! is formed from rotations acting on adjacent rows.

  type d_sweeps
     integer(kind=int32), private :: maxnum, n, maxord
     integer(kind=int32) :: num
     real(kind=dp), allocatable, dimension(:,:) :: cs, ss
     integer(kind=int32), allocatable, dimension(:) :: numrots
     integer(kind=int32), allocatable, dimension(:,:) :: js
  end type d_sweeps

  type c_sweeps
     integer(kind=int32), private :: maxnum, n, maxord
     integer(kind=int32) :: num
     complex(kind=dp), allocatable, dimension(:,:) :: cs, ss
     integer(kind=int32), allocatable, dimension(:) :: numrots
     integer(kind=int32), allocatable, dimension(:,:) :: js
  end type c_sweeps

  interface get_maxnum
     module procedure d_get_maxnum, c_get_maxnum
  end interface get_maxnum

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

  integer(kind=int32) function d_get_maxnum(sw) result(maxnum)
    type(d_sweeps) :: sw
    maxnum=sw%maxnum
  end function d_get_maxnum

  integer(kind=int32) function d_get_maxord(sw) result(maxord)
    type(d_sweeps) :: sw
    maxord=sw%maxord
  end function d_get_maxord

  integer(kind=int32) function c_get_maxord(sw) result(maxord)
    type(c_sweeps) :: sw
    maxord=sw%maxord
  end function c_get_maxord

  integer(kind=int32) function c_get_maxnum(sw) result(maxnum)
    type(c_sweeps) :: sw
    maxnum=sw%maxnum
  end function c_get_maxnum

  integer(kind=int32) function d_get_n_sweeps(sw) result(n)
    type(d_sweeps) :: sw
    n=sw%n
  end function d_get_n_sweeps

  integer(kind=int32) function c_get_n_sweeps(sw) result(n)
    type(c_sweeps) :: sw
    n=sw%n
  end function c_get_n_sweeps

  type(d_sweeps) function d_new_sweeps(n,maxnum,maxord) result(sw)
    integer(kind=int32), intent(in) :: n, maxnum, maxord
    sw%n=n; sw%maxnum=maxnum; sw%maxord=maxord
    allocate(sw%cs(maxord,maxnum), sw%ss(maxord,maxnum), sw%js(maxord,maxnum), sw%numrots(maxnum))
    sw%num=0
    sw%cs=1.0_dp; sw%ss=0.0_dp
  end function d_new_sweeps

  type(c_sweeps) function c_new_sweeps(n,maxnum,maxord) result(sw)
    integer(kind=int32), intent(in) :: n, maxnum, maxord
    sw%n=n; sw%maxnum=maxnum; sw%maxord=maxord
    allocate(sw%cs(maxord,maxnum), sw%ss(maxord,maxnum), sw%js(maxord,maxnum), sw%numrots(maxnum))
    sw%num=0
    sw%cs=1.0_dp; sw%ss=0.0_dp
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
!    do l=sw%num,1,-1
    do l=1,sw%num
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
!    do l=sw%num,1,-1
    do l=1,sw%num
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
    ! do l=sw%num,1,-1
    do l=1,sw%num
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
    ! do l=sw%num,1,-1
    do l=1,sw%num
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
    !    do l=1,sw%num
    do l=sw%num,1,-1
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
    ! do l=1,sw%num
    do l=sw%num,1,-1
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
    ! do l=1,sw%num
    do l=sw%num,1,-1
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
    ! do l=1,sw%num
    do l=sw%num,1,-1
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
    ! do l=1,sw%num
    do l=sw%num,1,-1
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
    ! do l=1,sw%num
    do l=sw%num,1,-1
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
    ! do l=1,sw%num
    do l=sw%num,1,-1
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
    ! do l=1,sw%num
    do l=sw%num,1,-1
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
    ! do l=sw%num,1,-1
    do l=1,sw%num
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
    ! do l=sw%num,1,-1
    do l=1,sw%num
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
    ! do l=sw%num,1,-1
    do l=1,sw%num
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
    ! do l=sw%num,1,-1
    do l=1,sw%num
       do j=1,sw%numrots(l)
          c=sw%cs(j,l); s=sw%ss(j,l)
          kk=sw%js(j,l)
          tmp=a(kk)
          a(kk)=c*tmp-s*a(kk+1)
          a(kk+1)=conjg(s)*tmp+c*a(kk+1)
       end do
    end do
  end subroutine c_v_general_times_trp_sweeps

end module sweeps
