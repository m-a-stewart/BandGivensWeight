module sweeps1
  use misc
  use types
  use shift
  use conversion

  implicit none

  !
  ! These types represent a linear transformation
  ! Q = Q_1 Q_2 ... Q_{numsweeps1}
  ! where
  ! Q_k = G_{k,1} ... G_{k,n-1}.
  ! G_{k,j} acts in rows j and j+1.  Thus Q is a product of upper Hessenberg matrices.


  type d_sweeps1
     integer(kind=int32), private :: maxsweeps1, n
     integer(kind=int32) :: numsweeps1
     real(kind=dp), allocatable, dimension(:,:) :: cs, ss
  end type d_sweeps1

  type c_sweeps1
     integer(kind=int32), private :: maxsweeps1, n
     integer(kind=int32) :: numsweeps1
     complex(kind=dp), allocatable, dimension(:,:) :: cs, ss
  end type c_sweeps1

  interface deallocate_sweeps1
     module procedure d_deallocate_sweeps1, c_deallocate_sweeps1
  end interface deallocate_sweeps1

  interface get_maxsweeps1
     module procedure d_get_maxsweeps1, c_get_maxsweeps1
  end interface get_maxsweeps1

  interface get_n
     module procedure d_get_n_sweeps1, c_get_n_sweeps1
  end interface get_n

  interface sweeps1_times_general
     module procedure d_sweeps1_times_general, c_sweeps1_times_general, &
          d_v_sweeps1_times_general, c_v_sweeps1_times_general
  end interface sweeps1_times_general

  interface trp_sweeps1_times_general
     module procedure d_trp_sweeps1_times_general, c_trp_sweeps1_times_general, &
          d_v_trp_sweeps1_times_general, c_v_trp_sweeps1_times_general
  end interface trp_sweeps1_times_general

  interface general_times_sweeps1
     module procedure d_general_times_sweeps1, c_general_times_sweeps1, &
          d_v_general_times_sweeps1, c_v_general_times_sweeps1
  end interface general_times_sweeps1

  interface general_times_trp_sweeps1
     module procedure d_general_times_trp_sweeps1, c_general_times_trp_sweeps1, &
          d_v_general_times_trp_sweeps1, c_v_general_times_trp_sweeps1
  end interface general_times_trp_sweeps1

  interface sweeps1_times_ub
     module procedure d_sweeps1_times_ub, c_sweeps1_times_ub
  end interface sweeps1_times_ub

  interface f_sweeps1_times_ub
     module procedure f_d_sweeps1_times_ub, f_c_sweeps1_times_ub
  end interface f_sweeps1_times_ub

  interface bv_times_sweeps1
     module procedure d_bv_times_sweeps1, c_bv_times_sweeps1
  end interface bv_times_sweeps1

  interface f_bv_times_sweeps1
     module procedure f_d_bv_times_sweeps1, f_c_bv_times_sweeps1
  end interface f_bv_times_sweeps1

  interface trp_sweeps1_times_bv
     module procedure d_trp_sweeps1_times_bv, c_trp_sweeps1_times_bv
  end interface trp_sweeps1_times_bv

  interface f_trp_sweeps1_times_bv
     module procedure f_d_trp_sweeps1_times_bv, f_c_trp_sweeps1_times_bv
  end interface f_trp_sweeps1_times_bv

  interface ub_times_trp_sweeps1
     module procedure d_ub_times_trp_sweeps1, c_ub_times_trp_sweeps1 
  end interface ub_times_trp_sweeps1

  interface f_ub_times_trp_sweeps1
     module procedure f_d_ub_times_trp_sweeps1, f_c_ub_times_trp_sweeps1 
  end interface f_ub_times_trp_sweeps1

  type(routine_info), parameter :: info_d_sweeps1_times_ub= &
       routine_info(id_d_sweeps1_times_ub, &
       'd_sweeps1_times_ub', &
       [ character(len=error_message_length) :: 'Sizes of bv, ub, and sw do not match.', &
       'Insufficienty lbw in bv.', 'Insufficient ubw in bv.', 'Insufficient lbw in ub', &
       'Insufficient ubw in ub', 'n<1' ])

  type(routine_info), parameter :: info_c_sweeps1_times_ub= &
       routine_info(id_c_sweeps1_times_ub, &
       'c_sweeps1_times_ub', &
       [ character(len=error_message_length) :: 'Sizes of bv, ub, and sw do not match.', &
       'Insufficienty lbw in bv.', 'Insufficient ubw in bv.', 'Insufficient lbw in ub', &
       'Insufficient ubw in ub', 'n<1' ])

  type(routine_info), parameter :: info_f_d_sweeps1_times_ub= &
       routine_info(id_f_d_sweeps1_times_ub, &
       'f_d_sweeps1_times_ub', &
       [ character(len=error_message_length) :: '' ])

  type(routine_info), parameter :: info_f_c_sweeps1_times_ub= &
       routine_info(id_f_c_sweeps1_times_ub, &
       'f_c_sweeps1_times_ub', &
       [ character(len=error_message_length) :: '' ])

  type(routine_info), parameter :: info_d_bv_times_sweeps1= &
       routine_info(id_d_bv_times_sweeps1, &
       'd_bv_times_sweeps1', &
       [ character(len=error_message_length) :: 'Sizes of bv, ub, and sw do not match.', &
       'Insufficient lbw in ub.', 'Insufficient ubw in ub.', 'Insufficient lbw in bv', &
       'Insufficient ubw in bv', 'n<1' ])

  type(routine_info), parameter :: info_f_d_bv_times_sweeps1= &
       routine_info(id_f_d_bv_times_sweeps1, &
       'f_d_bv_times_sweeps1', &
       [ character(len=error_message_length) :: '' ])

  type(routine_info), parameter :: info_d_trp_sweeps1_times_bv= &
       routine_info(id_d_trp_sweeps1_times_bv, &
       'd_trp_sweeps1_times_bv', &
       [ character(len=error_message_length) :: 'Sizes of bv, ub, and sw do not match.', &
       'Insufficient ubw in bv.', 'Insufficient ubw in ub.', 'n<1' ])

  type(routine_info), parameter :: info_f_d_trp_sweeps1_times_bv= &
       routine_info(id_f_d_trp_sweeps1_times_bv, &
       'f_d_trp_sweeps1_times_bv', &
       [ character(len=error_message_length) :: '' ])

  type(routine_info), parameter :: info_d_ub_times_trp_sweeps1= &
       routine_info(id_d_ub_times_trp_sweeps1, &
       'd_ub_times_trp_sweeps1', &
       [ character(len=error_message_length) :: 'Sizes of bv, ub, and sw do not match.', &
       'Insufficient ubw in bv.', 'Insufficient ubw in ub.', 'n<1' ])

  type(routine_info), parameter :: info_f_d_ub_times_trp_sweeps1= &
       routine_info(id_f_d_ub_times_trp_sweeps1, &
       'f_d_ub_times_trp_sweeps1', &
       [ character(len=error_message_length) :: '' ])

contains

  integer(kind=int32) function d_get_maxsweeps1(sw) result(maxsweeps1)
    type(d_sweeps1) :: sw
    maxsweeps1=sw%maxsweeps1
  end function d_get_maxsweeps1

  integer(kind=int32) function c_get_maxsweeps1(sw) result(maxsweeps1)
    type(c_sweeps1) :: sw
    maxsweeps1=sw%maxsweeps1
  end function c_get_maxsweeps1

  integer(kind=int32) function d_get_n_sweeps1(sw) result(n)
    type(d_sweeps1) :: sw
    n=sw%n
  end function d_get_n_sweeps1

  integer(kind=int32) function c_get_n_sweeps1(sw) result(n)
    type(c_sweeps1) :: sw
    n=sw%n
  end function c_get_n_sweeps1

  type(d_sweeps1) function d_new_sweeps1(n,maxsweeps1) result(sw)
    integer(kind=int32), intent(in) :: n, maxsweeps1
    sw%n=n; sw%maxsweeps1=maxsweeps1
    allocate(sw%cs(n-1,maxsweeps1), sw%ss(n-1,maxsweeps1))
    sw%numsweeps1=0
    sw%cs=1.0_dp; sw%ss=0.0_dp
  end function d_new_sweeps1

  type(c_sweeps1) function c_new_sweeps1(n,maxsweeps1) result(sw)
    integer(kind=int32), intent(in) :: n, maxsweeps1
    sw%n=n; sw%maxsweeps1=maxsweeps1
    allocate(sw%cs(n-1,maxsweeps1), sw%ss(n-1,maxsweeps1))
    sw%numsweeps1=0
    sw%cs=(1.0_dp,0.0_dp); sw%ss=(0.0_dp,0.0_dp)
  end function c_new_sweeps1

  subroutine d_deallocate_sweeps1(sw)
    type(d_sweeps1), intent(inout) :: sw
    deallocate(sw%cs, sw%ss)
  end subroutine d_deallocate_sweeps1

  subroutine c_deallocate_sweeps1(sw)
    type(c_sweeps1), intent(inout) :: sw
    deallocate(sw%cs, sw%ss)
  end subroutine c_deallocate_sweeps1

  subroutine d_sweeps1_times_general(sw,a)
    type(d_sweeps1) :: sw
    real(kind=dp), dimension(:,:), intent(inout) :: a
    type(d_rotation) :: rot
    integer(kind=int32) :: j, m, l
    do l=sw%numsweeps1,1,-1
       m=size(a,1)
       do j=m-1,1,-1
          rot%cosine=sw%cs(j,l); rot%sine=sw%ss(j,l)
          call d_rotation_times_general(rot,a,j,j+1)
       end do
    end do
  end subroutine d_sweeps1_times_general

  subroutine d_trp_sweeps1_times_general(sw,a)
    type(d_sweeps1) :: sw
    real(kind=dp), dimension(:,:), intent(inout) :: a
    type(d_rotation) :: rot
    integer(kind=int32) :: j, m, l
    do l=1,sw%numsweeps1
       m=size(a,1)
       do j=1,m-1
          rot%cosine=sw%cs(j,l); rot%sine=sw%ss(j,l)
          call d_rotation_times_general(trp_rot(rot),a,j,j+1)
       end do
    end do
  end subroutine d_trp_sweeps1_times_general

  subroutine d_v_sweeps1_times_general(sw,a)
    type(d_sweeps1) :: sw
    real(kind=dp), dimension(:), intent(inout) :: a
    integer(kind=int32) :: j, m, l
    real(kind=dp) :: tmp, c, s 
    do l=sw%numsweeps1,1,-1
       m=size(a)
       do j=m-1,1,-1
          tmp=a(j)
          c=sw%cs(j,l); s=sw%ss(j,l)
          a(j)=c*tmp-s*a(j+1)
          a(j+1)=s*tmp+c*a(j+1)
       end do
    end do
  end subroutine d_v_sweeps1_times_general

  subroutine d_v_trp_sweeps1_times_general(sw,a)
    type(d_sweeps1) :: sw
    real(kind=dp), dimension(:), intent(inout) :: a
    integer(kind=int32) :: j, m, l
    real(kind=dp) :: tmp, c, s 
    do l=1,sw%numsweeps1
       m=size(a)
       do j=1,m-1
          tmp=a(j)
          c=sw%cs(j,l); s=sw%ss(j,l)
          a(j)=c*tmp+s*a(j+1)
          a(j+1)=-s*tmp+c*a(j+1)
       end do
    end do
  end subroutine d_v_trp_sweeps1_times_general

  subroutine c_sweeps1_times_general(sw,a)
    type(c_sweeps1) :: sw
    complex(kind=dp), dimension(:,:), intent(inout) :: a
    type(c_rotation) :: rot
    integer(kind=int32) :: j, m, l
    do l=sw%numsweeps1,1,-1
       m=size(a,1)
       do j=m-1,1,-1
          rot%cosine=sw%cs(j,l); rot%sine=sw%ss(j,l)
          call c_rotation_times_general(rot,a,j,j+1)
       end do
    end do
  end subroutine c_sweeps1_times_general

  subroutine c_trp_sweeps1_times_general(sw,a)
    type(c_sweeps1) :: sw
    complex(kind=dp), dimension(:,:), intent(inout) :: a
    type(c_rotation) :: rot
    integer(kind=int32) :: j, m, l
    do l=1,sw%numsweeps1
       m=size(a,1)
       do j=1,m-1
          rot%cosine=sw%cs(j,l); rot%sine=sw%ss(j,l)
          call c_rotation_times_general(trp_rot(rot),a,j,j+1)
       end do
    end do
  end subroutine c_trp_sweeps1_times_general

  subroutine c_v_sweeps1_times_general(sw,a)
    type(c_sweeps1) :: sw
    complex(kind=dp), dimension(:), intent(inout) :: a
    integer(kind=int32) :: j, m, l
    complex(kind=dp) :: tmp, c, s 
    do l=sw%numsweeps1,1,-1
       m=size(a)
       do j=m-1,1,-1
          tmp=a(j)
          c=sw%cs(j,l); s=sw%ss(j,l)
          a(j)=c*tmp-conjg(s)*a(j+1)
          a(j+1)=s*tmp+c*a(j+1)
       end do
    end do
  end subroutine c_v_sweeps1_times_general

  subroutine c_v_trp_sweeps1_times_general(sw,a)
    type(c_sweeps1) :: sw
    complex(kind=dp), dimension(:), intent(inout) :: a
    integer(kind=int32) :: j, m, l
    complex(kind=dp) :: tmp, c, s 
    do l=1,sw%numsweeps1
       m=size(a)
       do j=1,m-1
          tmp=a(j)
          c=sw%cs(j,l); s=sw%ss(j,l)
          a(j)=c*tmp+conjg(s)*a(j+1)
          a(j+1)=-s*tmp+c*a(j+1)
       end do
    end do
  end subroutine c_v_trp_sweeps1_times_general

  subroutine d_general_times_sweeps1(a,sw)
    type(d_sweeps1) :: sw
    real(kind=dp), dimension(:,:), intent(inout) :: a
    type(d_rotation) :: rot
    integer(kind=int32) :: j, n, l
    do l=1,sw%numsweeps1
       n=size(a,2)
       do j=1,n-1
          rot%cosine=sw%cs(j,l); rot%sine=sw%ss(j,l)
          call d_general_times_rotation(a,rot,j,j+1)
       end do
    end do
  end subroutine d_general_times_sweeps1

  subroutine d_general_times_trp_sweeps1(a,sw)
    type(d_sweeps1) :: sw
    real(kind=dp), dimension(:,:), intent(inout) :: a
    type(d_rotation) :: rot
    integer(kind=int32) :: j, n, l
    do l=sw%numsweeps1,1,-1
       n=size(a,2)
       do j=n-1,1,-1
          rot%cosine=sw%cs(j,l); rot%sine=sw%ss(j,l)
          call d_general_times_rotation(a, trp_rot(rot),j,j+1)
       end do
    end do
  end subroutine d_general_times_trp_sweeps1

  subroutine d_v_general_times_sweeps1(a,sw)
    type(d_sweeps1) :: sw
    real(kind=dp), dimension(:), intent(inout) :: a
    integer(kind=int32) :: j, n, l
    real(kind=dp) :: c,s,tmp
    do l=1,sw%numsweeps1
       n=size(a)
       do j=1,n-1
          tmp=a(j)
          c=sw%cs(j,l); s=sw%ss(j,l)
          a(j)=c*tmp+s*a(j+1)
          a(j+1)=-s*tmp+c*a(j+1)
       end do
    end do
  end subroutine d_v_general_times_sweeps1

  subroutine d_v_general_times_trp_sweeps1(a,sw)
    type(d_sweeps1) :: sw
    real(kind=dp), dimension(:), intent(inout) :: a
    integer(kind=int32) :: j, n, l
    real(kind=dp) :: c,s,tmp
    do l=sw%numsweeps1,1,-1
       n=size(a)
       do j=n-1,1,-1
          tmp=a(j)
          c=sw%cs(j,l); s=sw%ss(j,l)
          a(j)=c*tmp-s*a(j+1)
          a(j+1)=s*tmp+c*a(j+1)
       end do
    end do
  end subroutine d_v_general_times_trp_sweeps1

  subroutine c_general_times_sweeps1(a,sw)
    type(c_sweeps1) :: sw
    complex(kind=dp), dimension(:,:), intent(inout) :: a
    type(c_rotation) :: rot
    integer(kind=int32) :: j, n, l
    do l=1,sw%numsweeps1
       n=size(a,2)
       do j=1,n-1
          rot%cosine=sw%cs(j,l); rot%sine=sw%ss(j,l)
          call c_general_times_rotation(a,rot,j,j+1)
       end do
    end do
  end subroutine c_general_times_sweeps1

  subroutine c_general_times_trp_sweeps1(a,sw)
    type(c_sweeps1) :: sw
    complex(kind=dp), dimension(:,:), intent(inout) :: a
    type(c_rotation) :: rot
    integer(kind=int32) :: j, n, l
    do l=sw%numsweeps1,1,-1
       n=size(a,2)
       do j=n-1,1,-1
          rot%cosine=sw%cs(j,l); rot%sine=sw%ss(j,l)
          call c_general_times_rotation(a, trp_rot(rot),j,j+1)
       end do
    end do
  end subroutine c_general_times_trp_sweeps1

  subroutine c_v_general_times_sweeps1(a,sw)
    type(c_sweeps1) :: sw
    complex(kind=dp), dimension(:), intent(inout) :: a
    integer(kind=int32) :: j, n, l
    complex(kind=dp) :: c,s,tmp
    do l=1,sw%numsweeps1
       n=size(a)
       do j=1,n-1
          tmp=a(j)
          c=sw%cs(j,l); s=sw%ss(j,l)
          a(j)=c*tmp+s*a(j+1)
          a(j+1)=-conjg(s)*tmp+c*a(j+1)
       end do
    end do
  end subroutine c_v_general_times_sweeps1

  subroutine c_v_general_times_trp_sweeps1(a,sw)
    type(c_sweeps1) :: sw
    complex(kind=dp), dimension(:), intent(inout) :: a
    integer(kind=int32) :: j, n, l
    complex(kind=dp) :: c,s,tmp
    do l=sw%numsweeps1,1,-1
       n=size(a)
       do j=n-1,1,-1
          tmp=a(j)
          c=sw%cs(j,l); s=sw%ss(j,l)
          a(j)=c*tmp-s*a(j+1)
          a(j+1)=conjg(s)*tmp+c*a(j+1)
       end do
    end do
  end subroutine c_v_general_times_trp_sweeps1


  ! sweeps1 times ub
  ! Errors:
  ! 1: Sizes of bv, ub, and sw do not match.
  ! 2: Insufficient lbw in bv.
  ! 3: Insufficient ubw in bv.
  ! 4: Insufficient lbw in ub
  ! 5: Insufficient ubw in ub.
  ! 6: n < 1.

  subroutine d_sweeps1_times_ub(sw,ub,bv,error)
    type(d_sweeps1) :: sw
    type(d_ub) :: ub
    type(d_bv) :: bv
    type(error_info) :: error
    call clear_error(error)
    if (get_n(bv) /= get_n(ub) .or. get_n(bv) /= get_n(sw)) then
       call set_error(error,1,id_d_sweeps1_times_ub); return
    end if
    if (get_lbwmax(bv) < sw%numsweeps1 + ub%lbw .and. get_lbwmax(bv) < get_n(ub)-1) then
       call set_error(error,2,id_d_sweeps1_times_ub); return
    end if
    if (get_ubwmax(bv) < sw%numsweeps1 + ub%ubw+1 .and. get_ubwmax(bv) < get_n(ub)-1) then
       call set_error(error,3,id_d_sweeps1_times_ub); return
    end if
    if (get_lbwmax(ub) < sw%numsweeps1 + ub%lbw .and. get_lbwmax(ub) < get_n(ub)-1) then
       call set_error(error,4,id_d_sweeps1_times_ub); return
    end if
    if (get_ubwmax(ub) < sw%numsweeps1 + ub%ubw+1 .and. get_ubwmax(ub) < get_n(ub)-1) then
       call set_error(error,5,id_d_sweeps1_times_ub); return
    end if
    if (get_n(bv) < 1) then
       call set_error(error,6,id_d_sweeps1_times_ub); return
    end if
    call f_d_sweeps1_times_ub(sw%cs, sw%ss, sw%numsweeps1, get_maxsweeps1(sw), get_n(sw), &
         ub%bc, ub%lbw, ub%ubw, get_lbwmax(ub), get_ubwmax(ub), &
         ub%numrotsu, ub%jsu, ub%csu, ub%ssu, &
         bv%br, bv%lbw, bv%ubw, get_lbwmax(bv), get_ubwmax(bv), &
         bv%numrotsv, bv%ksv, bv%csv, bv%ssv, error)
  end subroutine d_sweeps1_times_ub

  subroutine f_d_sweep_times_ub(cs_sw, ss_sw, n, &
       b_ub, lbw_ub, ubw_ub, lbwmax_ub, ubwmax_ub, numrotsu, jsu, csu, ssu, &
       b_bv, lbw_bv, ubw_bv, lbwmax_bv, ubwmax_bv, numrotsv, ksv, csv, ssv, error)

    real(kind=dp), dimension(n-1), intent(in) :: cs_sw, ss_sw
    integer(kind=int32), intent(in) :: n

    real(kind=dp), dimension(lbwmax_ub+ubwmax_ub+1,n), intent(inout) :: b_ub
    integer(kind=int32), dimension(n), intent(inout) :: numrotsu
    integer(kind=int32), dimension(ubwmax_ub,n), intent(inout) :: jsu
    real(kind=dp), dimension(ubwmax_ub,n), intent(inout) :: csu, ssu

    real(kind=dp), dimension(n, lbwmax_bv+ubwmax_bv+1), intent(out) :: b_bv
    integer(kind=int32), dimension(n), intent(out) :: numrotsv
    integer(kind=int32), dimension(n,ubwmax_bv), intent(out) :: ksv
    real(kind=dp), dimension(n,ubwmax_bv), intent(out) :: csv, ssv

    integer(kind=int32), intent(in) :: lbwmax_ub, ubwmax_ub
    integer(kind=int32), intent(inout) :: lbw_ub, ubw_ub
    integer(kind=int32), intent(in) :: lbwmax_bv, ubwmax_bv
    integer(kind=int32), intent(out) :: lbw_bv, ubw_bv

    type(error_info), intent(out) :: error
    integer(kind=int32) :: j, k, k0, k1, dlbw, dlbw_tmp, dubw, dubw_tmp
    type(d_rotation) :: rot
    call clear_error(error)

    if (n==1) then
       b_bv(1,1)=b_ub(1,1);
       lbw_bv=0; ubw_bv=0; numrotsv=0; return
    end if
    ! Initial expansion for first sweep: one extra subdiagonal and
    ! one extra superdiagonal to fill-in.
    ! Expand as need for later sweeps1.
    
    dlbw=1; dlbw_tmp=0
    dubw=1; dubw_tmp=1
    call f_bw_expand_bc(b_ub, n , lbw_ub, ubw_ub, lbwmax_ub, ubwmax_ub, &
         dlbw, dlbw_tmp, dubw, dubw_tmp, lbw_bv, ubw_bv)

    b_bv(:,1:lbw_bv+ubw_bv+1)=0.0_dp
    numrotsv=0
    ssv(:,1:ubw_bv)=0.0_dp
    csv(:,1:ubw_bv)=0.0_dp
    ksv(:,1:ubw_bv)=0

    do k=1,n-1
       do j=1,numrotsu(n-k)
          rot%cosine=csu(j,n-k); rot%sine=ssu(j,n-k)
          call rotation_times_tbc(rot,b_ub,n,lbw_ub,ubw_ub,n-k,0,jsu(j,n-k))
       end do
       ! Apply q_{n-k}
       rot%cosine=cs_sw(n-k); rot%sine=ss_sw(n-k)
       call rotation_times_tbc(rot,b_ub,n,lbw_ub,ubw_ub,0,0,n-k)
       ! zeros have appeared in the second superdiagonal in
       ! rows n-k-(ubw_ub-2), ... , n-k or
       ! columns n-k+2, ..., n-k+ubw_ub
       k0=max(n-k+2,ubw_ub+1)
       k1=min(n-k+ubw_ub,n)
       numrotsv(n-k)=max(k1-k0+1,0)
       do j=k0,k1
          rot=rgivens(get_el_bc(b_ub,ubw_ub,j-ubw_ub,j-1), get_el_bc(b_ub,ubw_ub,j-ubw_ub,j))
          ksv(n-k,k1-j+1)=j-1
          csv(n-k,k1-j+1)=rot%cosine; ssv(n-k,k1-j+1)=rot%sine
          call tbc_times_rotation(b_ub,n,lbw_ub,ubw_ub,0,k,rot,j-1)
       end do
    end do
    call f_bw_contract_bc(b_ub, n, lbw_ub, ubw_ub, lbwmax_ub, ubwmax_ub, &
         dlbw, dlbw_tmp, dubw, dubw_tmp)
    call bc_to_br(b_ub,b_bv,lbw_ub,ubw_ub)
  end subroutine f_d_sweep_times_ub

  subroutine f_d_sweeps1_times_ub(cs_sw, ss_sw, numsweeps1_sw, maxsweeps1_sw, n, &
       b_ub, lbw_ub, ubw_ub, lbwmax_ub, ubwmax_ub, numrotsu, jsu, csu, ssu, &
       b_bv, lbw_bv, ubw_bv, lbwmax_bv, ubwmax_bv, numrotsv, ksv, csv, ssv, error)

    real(kind=dp), dimension(n-1,maxsweeps1_sw), intent(in) :: cs_sw, ss_sw
    integer(kind=int32), intent(in) :: numsweeps1_sw, maxsweeps1_sw, n

    real(kind=dp), dimension(lbwmax_ub+ubwmax_ub+1,n), intent(inout) :: b_ub
    integer(kind=int32), dimension(n), intent(inout) :: numrotsu
    integer(kind=int32), dimension(ubwmax_ub,n), intent(inout) :: jsu
    real(kind=dp), dimension(ubwmax_ub,n), intent(inout) :: csu, ssu

    real(kind=dp), dimension(n, lbwmax_bv+ubwmax_bv+1), intent(out) :: b_bv
    integer(kind=int32), dimension(n), intent(out) :: numrotsv
    integer(kind=int32), dimension(n,ubwmax_bv), intent(out) :: ksv
    real(kind=dp), dimension(n,ubwmax_bv), intent(out) :: csv, ssv

    integer(kind=int32), intent(in) :: lbwmax_ub, ubwmax_ub
    integer(kind=int32), intent(inout) :: lbw_ub, ubw_ub
    integer(kind=int32), intent(in) :: lbwmax_bv, ubwmax_bv
    integer(kind=int32), intent(out) :: lbw_bv, ubw_bv

    type(error_info), intent(out) :: error
    integer(kind=int32) :: l

    call clear_error(error)

    if (n==1) then
       b_bv(1,1)=b_ub(1,1);
       lbw_bv=0; ubw_bv=0; numrotsv=0; return
    end if

    do l=numsweeps1_sw,1,-1
       if (l < numsweeps1_sw) then
          call f_d_convert_bv_to_ub(b_bv, n, lbw_bv, ubw_bv, lbwmax_bv, ubwmax_bv, &
               numrotsv, ksv, csv, ssv, b_ub, lbw_ub, ubw_ub, lbwmax_ub, ubwmax_ub, & 
               numrotsu, jsu, csu, ssu, error)
       end if
       
       call f_d_sweep_times_ub(cs_sw(:,l), ss_sw(:,l), n, &
            b_ub, lbw_ub, ubw_ub, lbwmax_ub, ubwmax_ub, numrotsu, jsu, csu, ssu, &
            b_bv, lbw_bv, ubw_bv, lbwmax_bv, ubwmax_bv, numrotsv, ksv, csv, ssv, error)
    end do
  end subroutine f_d_sweeps1_times_ub

  subroutine c_sweeps1_times_ub(sw,ub,bv,error)
    type(c_sweeps1) :: sw
    type(c_ub) :: ub
    type(c_bv) :: bv
    type(error_info) :: error
    call clear_error(error)
    if (get_n(bv) /= get_n(ub) .or. get_n(bv) /= get_n(sw)) then
       call set_error(error,1,id_c_sweeps1_times_ub); return
    end if
    if (get_lbwmax(bv) < sw%numsweeps1 + ub%lbw .and. get_lbwmax(bv) < get_n(ub)-1) then
       call set_error(error,2,id_c_sweeps1_times_ub); return
    end if
    if (get_ubwmax(bv) < sw%numsweeps1 + ub%ubw+1 .and. get_ubwmax(bv) < get_n(ub)-1) then
       call set_error(error,3,id_c_sweeps1_times_ub); return
    end if
    if (get_lbwmax(ub) < sw%numsweeps1 + ub%lbw .and. get_lbwmax(ub) < get_n(ub)-1) then
       call set_error(error,4,id_c_sweeps1_times_ub); return
    end if
    if (get_ubwmax(ub) < sw%numsweeps1 + ub%ubw+1 .and. get_ubwmax(ub) < get_n(ub)-1) then
       call set_error(error,5,id_c_sweeps1_times_ub); return
    end if
    if (get_n(bv) < 1) then
       call set_error(error,6,id_c_sweeps1_times_ub); return
    end if
    call f_c_sweeps1_times_ub(sw%cs, sw%ss, sw%numsweeps1, get_maxsweeps1(sw), get_n(sw), &
         ub%bc, ub%lbw, ub%ubw, get_lbwmax(ub), get_ubwmax(ub), &
         ub%numrotsu, ub%jsu, ub%csu, ub%ssu, &
         bv%br, bv%lbw, bv%ubw, get_lbwmax(bv), get_ubwmax(bv), &
         bv%numrotsv, bv%ksv, bv%csv, bv%ssv, error)
  end subroutine c_sweeps1_times_ub

  subroutine f_c_sweep_times_ub(cs_sw, ss_sw, n, &
       b_ub, lbw_ub, ubw_ub, lbwmax_ub, ubwmax_ub, numrotsu, jsu, csu, ssu, &
       b_bv, lbw_bv, ubw_bv, lbwmax_bv, ubwmax_bv, numrotsv, ksv, csv, ssv, error)

    complex(kind=dp), dimension(n-1), intent(in) :: cs_sw, ss_sw
    integer(kind=int32), intent(in) :: n

    complex(kind=dp), dimension(lbwmax_ub+ubwmax_ub+1,n), intent(inout) :: b_ub
    integer(kind=int32), dimension(n), intent(inout) :: numrotsu
    integer(kind=int32), dimension(ubwmax_ub,n), intent(inout) :: jsu
    complex(kind=dp), dimension(ubwmax_ub,n), intent(inout) :: csu, ssu

    complex(kind=dp), dimension(n, lbwmax_bv+ubwmax_bv+1), intent(out) :: b_bv
    integer(kind=int32), dimension(n), intent(out) :: numrotsv
    integer(kind=int32), dimension(n,ubwmax_bv), intent(out) :: ksv
    complex(kind=dp), dimension(n,ubwmax_bv), intent(out) :: csv, ssv

    integer(kind=int32), intent(in) :: lbwmax_ub, ubwmax_ub
    integer(kind=int32), intent(inout) :: lbw_ub, ubw_ub
    integer(kind=int32), intent(in) :: lbwmax_bv, ubwmax_bv
    integer(kind=int32), intent(out) :: lbw_bv, ubw_bv

    type(error_info), intent(out) :: error
    integer(kind=int32) :: j, k, k0, k1, dlbw, dlbw_tmp, dubw, dubw_tmp
    type(c_rotation) :: rot
    call clear_error(error)

    if (n==1) then
       b_bv(1,1)=b_ub(1,1);
       lbw_bv=0; ubw_bv=0; numrotsv=0; return
    end if
    ! Initial expansion for first sweep: one extra subdiagonal and
    ! one extra superdiagonal to fill-in.
    ! Expand as need for later sweeps1.
    
    dlbw=1; dlbw_tmp=0
    dubw=1; dubw_tmp=1
    call f_bw_expand_bc(b_ub, n , lbw_ub, ubw_ub, lbwmax_ub, ubwmax_ub, &
         dlbw, dlbw_tmp, dubw, dubw_tmp, lbw_bv, ubw_bv)

    b_bv(:,1:lbw_bv+ubw_bv+1)=(0.0_dp,0.0_dp)
    numrotsv=0
    ssv(:,1:ubw_bv)=(0.0_dp,0.0_dp)
    csv(:,1:ubw_bv)=(0.0_dp,0.0_dp)
    ksv(:,1:ubw_bv)=0

    do k=1,n-1
       do j=1,numrotsu(n-k)
          rot%cosine=csu(j,n-k); rot%sine=ssu(j,n-k)
          call rotation_times_tbc(rot,b_ub,n,lbw_ub,ubw_ub,n-k,0,jsu(j,n-k))
       end do
       ! Apply q_{n-k}
       rot%cosine=cs_sw(n-k); rot%sine=ss_sw(n-k)
       call rotation_times_tbc(rot,b_ub,n,lbw_ub,ubw_ub,0,0,n-k)
       ! zeros have appeared in the second superdiagonal in
       ! rows n-k-(ubw_ub-2), ... , n-k or
       ! columns n-k+2, ..., n-k+ubw_ub
       k0=max(n-k+2,ubw_ub+1)
       k1=min(n-k+ubw_ub,n)
       numrotsv(n-k)=max(k1-k0+1,0)
       do j=k0,k1
          rot=rgivens(get_el_bc(b_ub,ubw_ub,j-ubw_ub,j-1), get_el_bc(b_ub,ubw_ub,j-ubw_ub,j))
          ksv(n-k,k1-j+1)=j-1
          csv(n-k,k1-j+1)=rot%cosine; ssv(n-k,k1-j+1)=rot%sine
          call tbc_times_rotation(b_ub,n,lbw_ub,ubw_ub,0,k,rot,j-1)
       end do
    end do
    call f_bw_contract_bc(b_ub, n, lbw_ub, ubw_ub, lbwmax_ub, ubwmax_ub, &
         dlbw, dlbw_tmp, dubw, dubw_tmp)
    call bc_to_br(b_ub,b_bv,lbw_ub,ubw_ub)
  end subroutine f_c_sweep_times_ub

  subroutine f_c_sweeps1_times_ub(cs_sw, ss_sw, numsweeps1_sw, maxsweeps1_sw, n, &
       b_ub, lbw_ub, ubw_ub, lbwmax_ub, ubwmax_ub, numrotsu, jsu, csu, ssu, &
       b_bv, lbw_bv, ubw_bv, lbwmax_bv, ubwmax_bv, numrotsv, ksv, csv, ssv, error)

    complex(kind=dp), dimension(n-1,maxsweeps1_sw), intent(in) :: cs_sw, ss_sw
    integer(kind=int32), intent(in) :: numsweeps1_sw, maxsweeps1_sw, n

    complex(kind=dp), dimension(lbwmax_ub+ubwmax_ub+1,n), intent(inout) :: b_ub
    integer(kind=int32), dimension(n), intent(inout) :: numrotsu
    integer(kind=int32), dimension(ubwmax_ub,n), intent(inout) :: jsu
    complex(kind=dp), dimension(ubwmax_ub,n), intent(inout) :: csu, ssu

    complex(kind=dp), dimension(n, lbwmax_bv+ubwmax_bv+1), intent(out) :: b_bv
    integer(kind=int32), dimension(n), intent(out) :: numrotsv
    integer(kind=int32), dimension(n,ubwmax_bv), intent(out) :: ksv
    complex(kind=dp), dimension(n,ubwmax_bv), intent(out) :: csv, ssv

    integer(kind=int32), intent(in) :: lbwmax_ub, ubwmax_ub
    integer(kind=int32), intent(inout) :: lbw_ub, ubw_ub
    integer(kind=int32), intent(in) :: lbwmax_bv, ubwmax_bv
    integer(kind=int32), intent(out) :: lbw_bv, ubw_bv

    type(error_info), intent(out) :: error
    integer(kind=int32) :: l

    call clear_error(error)

    if (n==1) then
       b_bv(1,1)=b_ub(1,1);
       lbw_bv=0; ubw_bv=0; numrotsv=0; return
    end if

    do l=numsweeps1_sw,1,-1
       if (l < numsweeps1_sw) then
          call f_c_convert_bv_to_ub(b_bv, n, lbw_bv, ubw_bv, lbwmax_bv, ubwmax_bv, &
               numrotsv, ksv, csv, ssv, b_ub, lbw_ub, ubw_ub, lbwmax_ub, ubwmax_ub, & 
               numrotsu, jsu, csu, ssu, error)
       end if
       
       call f_c_sweep_times_ub(cs_sw(:,l), ss_sw(:,l), n, &
            b_ub, lbw_ub, ubw_ub, lbwmax_ub, ubwmax_ub, numrotsu, jsu, csu, ssu, &
            b_bv, lbw_bv, ubw_bv, lbwmax_bv, ubwmax_bv, numrotsv, ksv, csv, ssv, error)
    end do
  end subroutine f_c_sweeps1_times_ub

  ! bv times sweeps1
  ! Errors:
  ! 1: Sizes of bv, ub, and sw do not match.
  ! 2: Insufficient lbw in ub.
  ! 3: Insufficient ubw in ub.
  ! 4: Insufficient lbw in bv
  ! 5: Insufficient ubw in bv.
  ! 6: n < 1.

  subroutine d_bv_times_sweeps1(bv,sw,ub,error)
    type(d_sweeps1) :: sw
    type(d_ub) :: ub
    type(d_bv) :: bv
    type(error_info) :: error
    call clear_error(error)
    if (get_n(bv) /= get_n(ub) .or. get_n(bv) /= get_n(sw)) then
       call set_error(error,1,id_d_bv_times_sweeps1); return
    end if
    if (get_lbwmax(ub) < sw%numsweeps1 + bv%lbw .and. get_lbwmax(ub) < get_n(bv)-1) then
       call set_error(error,2,id_d_bv_times_sweeps1); return
    end if
    if (get_ubwmax(ub) < sw%numsweeps1 + bv%ubw+1 .and. get_ubwmax(ub) < get_n(bv)-1) then
       call set_error(error,3,id_d_bv_times_sweeps1); return
    end if
    if (get_lbwmax(bv) < sw%numsweeps1 + bv%lbw .and. get_lbwmax(bv) < get_n(bv)-1) then
       call set_error(error,4,id_d_bv_times_sweeps1); return
    end if
    if (get_ubwmax(bv) < sw%numsweeps1 + bv%ubw+1 .and. get_ubwmax(bv) < get_n(bv)-1) then
       call set_error(error,5,id_d_bv_times_sweeps1); return
    end if
    if (get_n(bv) < 1) then
       call set_error(error,6,id_d_bv_times_sweeps1); return
    end if
    call f_d_bv_times_sweeps1(bv%br, get_n(bv), bv%lbw, bv%ubw, get_lbwmax(bv), get_ubwmax(bv), &
         bv%numrotsv, bv%ksv, bv%csv, bv%ssv, &
         sw%cs, sw%ss, sw%numsweeps1, get_maxsweeps1(sw), &
         ub%bc, ub%lbw, ub%ubw, get_lbwmax(ub), get_ubwmax(ub), &
         ub%numrotsu, ub%jsu, ub%csu, ub%ssu, error)
  end subroutine d_bv_times_sweeps1

  subroutine f_d_bv_times_sweep(b_bv, n, lbw_bv, ubw_bv, lbwmax_bv, ubwmax_bv, numrotsv, ksv, & 
       csv, ssv, &
       cs_sw, ss_sw, numsweeps1_sw, maxsweeps1_sw, &
       b_ub, lbw_ub, ubw_ub, lbwmax_ub, ubwmax_ub, numrotsu, jsu, csu, ssu, error)

    real(kind=dp), dimension(n-1), intent(in) :: cs_sw, ss_sw
    integer(kind=int32), intent(in) :: numsweeps1_sw, maxsweeps1_sw, n

    real(kind=dp), dimension(n, lbwmax_bv+ubwmax_bv+1), intent(inout) :: b_bv
    integer(kind=int32), dimension(n), intent(inout) :: numrotsv
    integer(kind=int32), dimension(n,ubwmax_bv), intent(inout) :: ksv
    real(kind=dp), dimension(n,ubwmax_bv), intent(inout) :: csv, ssv

    real(kind=dp), dimension(lbwmax_ub+ubwmax_ub+1,n), intent(out) :: b_ub
    integer(kind=int32), dimension(n), intent(out) :: numrotsu
    integer(kind=int32), dimension(ubwmax_ub,n), intent(out) :: jsu
    real(kind=dp), dimension(ubwmax_ub,n), intent(out) :: csu, ssu

    integer(kind=int32), intent(in) :: lbwmax_ub, ubwmax_ub
    integer(kind=int32), intent(inout) :: lbw_ub, ubw_ub
    integer(kind=int32), intent(in) :: lbwmax_bv, ubwmax_bv
    integer(kind=int32), intent(out) :: lbw_bv, ubw_bv

    type(error_info), intent(out) :: error
    integer(kind=int32) :: j, k, k0, k1, dlbw, dlbw_tmp, dubw, dubw_tmp
    type(d_rotation) :: rot

    call clear_error(error)

    if (n==1) then
       b_ub(1,1)=b_bv(1,1);
       lbw_ub=0; ubw_ub=0; numrotsu=0; return
    end if

    dlbw=1; dlbw_tmp=0
    dubw=1; dubw_tmp=1
    call f_bw_expand_br(b_bv, n, lbw_bv, ubw_bv, lbwmax_bv, ubwmax_bv, &
         dlbw, dlbw_tmp, dubw, dubw_tmp, lbw_ub, ubw_ub)

    b_ub(1:lbw_ub+ubw_ub+1,:)=0.0_dp
    numrotsu=0
    ssu(1:ubw_ub,:)=0.0_dp; csu(1:ubw_ub,:)=0.0_dp
    jsu(1:ubw_ub,:)=0
    
    do k=1,n-1
       do j=1,numrotsv(k)
          rot%cosine=csv(k,j); rot%sine=ssv(k,j)
          call tbr_times_rotation(b_bv,n,lbw_bv,ubw_bv,0,n-k,trp_rot(rot),ksv(k,j))
       end do
       ! Apply q_{k}
       rot%cosine=cs_sw(k); rot%sine=ss_sw(k)
       call tbr_times_rotation(b_bv,n,lbw_bv,ubw_bv,0,0,rot,k)
       ! zeros have appeared in the second superdiagonal
       ! in columns k+1, ..., k+ubw_bv-1
       ! or rows k+1-ubw_bv, ..., k-1
       k0=max(k+1,ubw_bv+1)
       k1=min(k+ubw_bv-1,n)
       numrotsu(k)=max(k1-k0+1,0)
       do j=k1,k0,-1
          rot=lgivens2(get_el_br(b_bv,lbw_bv,j-ubw_bv,j), get_el_br(b_bv,lbw_bv,j-ubw_bv+1,j))
          jsu(j-k0+1,k)=j-ubw_bv
          csu(j-k0+1,k)=rot%cosine; ssu(j-k0+1,k)=rot%sine
          call rotation_times_tbr(trp_rot(rot), b_bv,n,lbw_bv,ubw_bv,k,0,j-ubw_bv)
       end do
    end do
    
    call f_bw_contract_br(b_bv, n, lbw_bv, ubw_bv, lbwmax_bv, ubwmax_bv, &
         dlbw, dlbw_tmp, dubw, dubw_tmp)
    call br_to_bc(b_bv,b_ub,lbw_ub,ubw_ub)
  end subroutine f_d_bv_times_sweep

  subroutine f_d_bv_times_sweeps1(b_bv, n, lbw_bv, ubw_bv, lbwmax_bv, ubwmax_bv, numrotsv, ksv, & 
       csv, ssv, &
       cs_sw, ss_sw, numsweeps1_sw, maxsweeps1_sw, &
       b_ub, lbw_ub, ubw_ub, lbwmax_ub, ubwmax_ub, numrotsu, jsu, csu, ssu, error)

    real(kind=dp), dimension(n-1,maxsweeps1_sw), intent(in) :: cs_sw, ss_sw
    integer(kind=int32), intent(in) :: numsweeps1_sw, maxsweeps1_sw, n

    real(kind=dp), dimension(n, lbwmax_bv+ubwmax_bv+1), intent(inout) :: b_bv
    integer(kind=int32), dimension(n), intent(inout) :: numrotsv
    integer(kind=int32), dimension(n,ubwmax_bv), intent(inout) :: ksv
    real(kind=dp), dimension(n,ubwmax_bv), intent(inout) :: csv, ssv

    real(kind=dp), dimension(lbwmax_ub+ubwmax_ub+1,n), intent(out) :: b_ub
    integer(kind=int32), dimension(n), intent(out) :: numrotsu
    integer(kind=int32), dimension(ubwmax_ub,n), intent(out) :: jsu
    real(kind=dp), dimension(ubwmax_ub,n), intent(out) :: csu, ssu

    integer(kind=int32), intent(in) :: lbwmax_ub, ubwmax_ub
    integer(kind=int32), intent(inout) :: lbw_ub, ubw_ub
    integer(kind=int32), intent(in) :: lbwmax_bv, ubwmax_bv
    integer(kind=int32), intent(out) :: lbw_bv, ubw_bv

    type(error_info), intent(out) :: error
    integer(kind=int32) :: l

    call clear_error(error)

    if (n==1) then
       b_ub(1,1)=b_bv(1,1);
       lbw_ub=0; ubw_ub=0; numrotsu=0; return
    end if
    do l=1, numsweeps1_sw
       if (l > 1) then
          call f_d_convert_ub_to_bv(b_ub, n, lbw_ub, ubw_ub, lbwmax_ub, ubwmax_ub, & 
               numrotsu, jsu, csu, ssu, &
               b_bv, lbw_bv, ubw_bv, lbwmax_bv, ubwmax_bv, numrotsv, ksv, csv, ssv, error)
       end if
       call f_d_bv_times_sweep(b_bv, n, lbw_bv, ubw_bv, lbwmax_bv, ubwmax_bv, numrotsv, ksv, & 
            csv, ssv, &
            cs_sw(:,l), ss_sw(:,l), numsweeps1_sw, maxsweeps1_sw, &
            b_ub, lbw_ub, ubw_ub, lbwmax_ub, ubwmax_ub, numrotsu, jsu, csu, ssu, error)
    end do

  end subroutine f_d_bv_times_sweeps1

  subroutine c_bv_times_sweeps1(bv,sw,ub,error)
    type(c_sweeps1) :: sw
    type(c_ub) :: ub
    type(c_bv) :: bv
    type(error_info) :: error
    call clear_error(error)
    if (get_n(bv) /= get_n(ub) .or. get_n(bv) /= get_n(sw)) then
       call set_error(error,1,id_c_bv_times_sweeps1); return
    end if
    if (get_lbwmax(ub) < sw%numsweeps1 + bv%lbw .and. get_lbwmax(ub) < get_n(bv)-1) then
       call set_error(error,2,id_c_bv_times_sweeps1); return
    end if
    if (get_ubwmax(ub) < sw%numsweeps1 + bv%ubw+1 .and. get_ubwmax(ub) < get_n(bv)-1) then
       call set_error(error,3,id_c_bv_times_sweeps1); return
    end if
    if (get_lbwmax(bv) < sw%numsweeps1 + bv%lbw .and. get_lbwmax(bv) < get_n(bv)-1) then
       call set_error(error,4,id_c_bv_times_sweeps1); return
    end if
    if (get_ubwmax(bv) < sw%numsweeps1 + bv%ubw+1 .and. get_ubwmax(bv) < get_n(bv)-1) then
       call set_error(error,5,id_c_bv_times_sweeps1); return
    end if
    if (get_n(bv) < 1) then
       call set_error(error,6,id_c_bv_times_sweeps1); return
    end if
    call f_c_bv_times_sweeps1(bv%br, get_n(bv), bv%lbw, bv%ubw, get_lbwmax(bv), get_ubwmax(bv), &
         bv%numrotsv, bv%ksv, bv%csv, bv%ssv, &
         sw%cs, sw%ss, sw%numsweeps1, get_maxsweeps1(sw), &
         ub%bc, ub%lbw, ub%ubw, get_lbwmax(ub), get_ubwmax(ub), &
         ub%numrotsu, ub%jsu, ub%csu, ub%ssu, error)
  end subroutine c_bv_times_sweeps1

  subroutine f_c_bv_times_sweep(b_bv, n, lbw_bv, ubw_bv, lbwmax_bv, ubwmax_bv, numrotsv, ksv, & 
       csv, ssv, &
       cs_sw, ss_sw, numsweeps1_sw, maxsweeps1_sw, &
       b_ub, lbw_ub, ubw_ub, lbwmax_ub, ubwmax_ub, numrotsu, jsu, csu, ssu, error)

    complex(kind=dp), dimension(n-1), intent(in) :: cs_sw, ss_sw
    integer(kind=int32), intent(in) :: numsweeps1_sw, maxsweeps1_sw, n

    complex(kind=dp), dimension(n, lbwmax_bv+ubwmax_bv+1), intent(inout) :: b_bv
    integer(kind=int32), dimension(n), intent(inout) :: numrotsv
    integer(kind=int32), dimension(n,ubwmax_bv), intent(inout) :: ksv
    complex(kind=dp), dimension(n,ubwmax_bv), intent(inout) :: csv, ssv

    complex(kind=dp), dimension(lbwmax_ub+ubwmax_ub+1,n), intent(out) :: b_ub
    integer(kind=int32), dimension(n), intent(out) :: numrotsu
    integer(kind=int32), dimension(ubwmax_ub,n), intent(out) :: jsu
    complex(kind=dp), dimension(ubwmax_ub,n), intent(out) :: csu, ssu

    integer(kind=int32), intent(in) :: lbwmax_ub, ubwmax_ub
    integer(kind=int32), intent(inout) :: lbw_ub, ubw_ub
    integer(kind=int32), intent(in) :: lbwmax_bv, ubwmax_bv
    integer(kind=int32), intent(out) :: lbw_bv, ubw_bv

    type(error_info), intent(out) :: error
    integer(kind=int32) :: j, k, k0, k1, dlbw, dlbw_tmp, dubw, dubw_tmp
    type(c_rotation) :: rot

    call clear_error(error)

    if (n==1) then
       b_ub(1,1)=b_bv(1,1);
       lbw_ub=0; ubw_ub=0; numrotsu=0; return
    end if

    dubw=1; dubw_tmp=1
    dlbw=1; dlbw_tmp=0
    call f_bw_expand_br(b_bv, n, lbw_bv, ubw_bv, lbwmax_bv, ubwmax_bv, &
         dlbw, dlbw_tmp, dubw, dubw_tmp, lbw_ub, ubw_ub)

    b_ub(1:lbw_ub+ubw_ub+1,:)=(0.0_dp,0.0_dp)
    numrotsu=0
    ssu(1:ubw_ub,:)=(0.0_dp,0.0_dp); csu(1:ubw_ub,:)=(0.0_dp,0.0_dp)
    jsu(1:ubw_ub,:)=0
    
    do k=1,n-1
       do j=1,numrotsv(k)
          rot%cosine=csv(k,j); rot%sine=ssv(k,j)
          call tbr_times_rotation(b_bv,n,lbw_bv,ubw_bv,0,n-k,trp_rot(rot),ksv(k,j))
       end do
       ! Apply q_{k}
       rot%cosine=cs_sw(k); rot%sine=ss_sw(k)
       call tbr_times_rotation(b_bv,n,lbw_bv,ubw_bv,0,0,rot,k)
       ! zeros have appeared in the second superdiagonal
       ! in columns k+1, ..., k+ubw_bv-1
       ! or rows k+1-ubw_bv, ..., k-1
       k0=max(k+1,ubw_bv+1)
       k1=min(k+ubw_bv-1,n)
       numrotsu(k)=max(k1-k0+1,0)
       do j=k1,k0,-1
          rot=lgivens2(get_el_br(b_bv,lbw_bv,j-ubw_bv,j), get_el_br(b_bv,lbw_bv,j-ubw_bv+1,j))
          jsu(j-k0+1,k)=j-ubw_bv
          csu(j-k0+1,k)=rot%cosine; ssu(j-k0+1,k)=rot%sine
          call rotation_times_tbr(trp_rot(rot), b_bv,n,lbw_bv,ubw_bv,k,0,j-ubw_bv)
       end do
    end do
    call f_bw_contract_br(b_bv, n, lbw_bv, ubw_bv, lbwmax_bv, ubwmax_bv, &
         dlbw, dlbw_tmp, dubw, dubw_tmp)
    call br_to_bc(b_bv,b_ub,lbw_ub,ubw_ub)
  end subroutine f_c_bv_times_sweep

  subroutine f_c_bv_times_sweeps1(b_bv, n, lbw_bv, ubw_bv, lbwmax_bv, ubwmax_bv, numrotsv, ksv, & 
       csv, ssv, &
       cs_sw, ss_sw, numsweeps1_sw, maxsweeps1_sw, &
       b_ub, lbw_ub, ubw_ub, lbwmax_ub, ubwmax_ub, numrotsu, jsu, csu, ssu, error)

    complex(kind=dp), dimension(n-1,maxsweeps1_sw), intent(in) :: cs_sw, ss_sw
    integer(kind=int32), intent(in) :: numsweeps1_sw, maxsweeps1_sw, n

    complex(kind=dp), dimension(n, lbwmax_bv+ubwmax_bv+1), intent(inout) :: b_bv
    integer(kind=int32), dimension(n), intent(inout) :: numrotsv
    integer(kind=int32), dimension(n,ubwmax_bv), intent(inout) :: ksv
    complex(kind=dp), dimension(n,ubwmax_bv), intent(inout) :: csv, ssv

    complex(kind=dp), dimension(lbwmax_ub+ubwmax_ub+1,n), intent(out) :: b_ub
    integer(kind=int32), dimension(n), intent(out) :: numrotsu
    integer(kind=int32), dimension(ubwmax_ub,n), intent(out) :: jsu
    complex(kind=dp), dimension(ubwmax_ub,n), intent(out) :: csu, ssu

    integer(kind=int32), intent(in) :: lbwmax_ub, ubwmax_ub
    integer(kind=int32), intent(inout) :: lbw_ub, ubw_ub
    integer(kind=int32), intent(in) :: lbwmax_bv, ubwmax_bv
    integer(kind=int32), intent(out) :: lbw_bv, ubw_bv

    type(error_info), intent(out) :: error
    integer(kind=int32) :: l

    call clear_error(error)

    if (n==1) then
       b_ub(1,1)=b_bv(1,1);
       lbw_ub=0; ubw_ub=0; numrotsu=0; return
    end if
    do l=1, numsweeps1_sw
       if (l > 1) then
          call f_c_convert_ub_to_bv(b_ub, n, lbw_ub, ubw_ub, lbwmax_ub, ubwmax_ub, & 
               numrotsu, jsu, csu, ssu, &
               b_bv, lbw_bv, ubw_bv, lbwmax_bv, ubwmax_bv, numrotsv, ksv, csv, ssv, error)
       end if
       call f_c_bv_times_sweep(b_bv, n, lbw_bv, ubw_bv, lbwmax_bv, ubwmax_bv, numrotsv, ksv, & 
            csv, ssv, &
            cs_sw(:,l), ss_sw(:,l), numsweeps1_sw, maxsweeps1_sw, &
            b_ub, lbw_ub, ubw_ub, lbwmax_ub, ubwmax_ub, numrotsu, jsu, csu, ssu, error)
    end do

  end subroutine f_c_bv_times_sweeps1

  ! Transposed sweeps1: The following assume bounded fill-in in the
  ! lower triangular part.  In particular, it is assumed that the lbw
  ! will not increase.  If there is partial fill-in, additional zero
  ! subdiagonals should be included to accomodate it.
  !
  ! errors:
  ! 1: Sizes of bv, ub, and sw do not match.
  ! 2: Insufficient ubw in bv.
  ! 3: Insufficient ubw in ub.
  ! 4: n < 1.

  subroutine d_trp_sweeps1_times_bv(sw,bv,ub,error)
    type(d_sweeps1) :: sw
    type(d_ub) :: ub
    type(d_bv) :: bv
    type(error_info) :: error
    call clear_error(error)
    if (get_n(bv) /= get_n(ub) .or. get_n(bv) /= get_n(sw)) then
       call set_error(error,1,id_d_trp_sweeps1_times_bv); return
    end if
    if (get_ubwmax(bv) < sw%numsweeps1 + ub%ubw+1 .and. get_ubwmax(bv) < get_n(ub)-1) then
       call set_error(error,2,id_d_trp_sweeps1_times_bv); return
    end if
    if (get_ubwmax(ub) < sw%numsweeps1 + ub%ubw+1 .and. get_ubwmax(ub) < get_n(ub)-1) then
       call set_error(error,3,id_d_trp_sweeps1_times_bv); return
    end if
    if (get_n(bv) < 1) then
       call set_error(error,4,id_d_trp_sweeps1_times_bv); return
    end if
    call f_d_trp_sweeps1_times_bv(sw%cs, sw%ss, get_n(sw), sw%numsweeps1, get_maxsweeps1(sw), &
         bv%br, bv%lbw, bv%ubw, get_lbwmax(bv), get_ubwmax(bv), &
         bv%numrotsv, bv%ksv, bv%csv, bv%ssv, &
         ub%bc, ub%lbw, ub%ubw, get_lbwmax(ub), get_ubwmax(ub), &
         ub%numrotsu, ub%jsu, ub%csu, ub%ssu, error)
  end subroutine d_trp_sweeps1_times_bv

  subroutine f_d_trp_sweep_times_bv(cs_sw, ss_sw, n, numsweeps1_sw, maxsweeps1_sw, &
       b_bv, lbw_bv, ubw_bv, lbwmax_bv, ubwmax_bv, numrotsv, ksv, csv, ssv, &
       b_ub, lbw_ub, ubw_ub, lbwmax_ub, ubwmax_ub, numrotsu, jsu, csu, ssu, error)

    real(kind=dp), dimension(n-1), intent(in) :: cs_sw, ss_sw
    integer(kind=int32), intent(in) :: numsweeps1_sw, maxsweeps1_sw, n

    real(kind=dp), dimension(n, lbwmax_bv+ubwmax_bv+1), intent(inout) :: b_bv
    integer(kind=int32), dimension(n), intent(inout) :: numrotsv
    integer(kind=int32), dimension(n,ubwmax_bv), intent(inout) :: ksv
    real(kind=dp), dimension(n,ubwmax_bv), intent(inout) :: csv, ssv

    real(kind=dp), dimension(lbwmax_ub+ubwmax_ub+1,n), intent(out) :: b_ub
    integer(kind=int32), dimension(n), intent(out) :: numrotsu
    integer(kind=int32), dimension(ubwmax_ub,n), intent(out) :: jsu
    real(kind=dp), dimension(ubwmax_ub,n), intent(out) :: csu, ssu

    integer(kind=int32), intent(in) :: lbwmax_ub, ubwmax_ub
    integer(kind=int32), intent(inout) :: lbw_ub, ubw_ub
    integer(kind=int32), intent(in) :: lbwmax_bv, ubwmax_bv
    integer(kind=int32), intent(out) :: lbw_bv, ubw_bv

    type(error_info), intent(out) :: error
    integer(kind=int32) :: j, k, k0, k1, dlbw, dlbw_tmp, dubw, dubw_tmp
    type(d_rotation) :: rot

    call clear_error(error)

    if (n==1) then
       b_ub(1,1)=b_bv(1,1);
       lbw_ub=0; ubw_ub=0; numrotsu=0; return
    end if

    dlbw=0; dlbw_tmp=0
    dubw=1; dubw_tmp=1
    
    call f_bw_expand_br(b_bv, n, lbw_bv, ubw_bv, lbwmax_bv, ubwmax_bv, &
         dlbw, dlbw_tmp, dubw, dubw_tmp, lbw_ub, ubw_ub)

    b_ub(1:lbw_ub+ubw_ub+1,:)=0.0_dp
    numrotsu=0
    ssu(1:ubw_ub,:)=0.0_dp; csu(1:ubw_ub,:)=0.0_dp
    jsu(1:ubw_ub,:)=0

    do k=1,n-1
       ! apply v_{n-k}
       do j=1,numrotsv(k)
          rot%cosine=csv(k,j); rot%sine=ssv(k,j)
          call tbr_times_rotation(b_bv,n,lbw_bv,ubw_bv,0,n-k,trp_rot(rot),ksv(k,j))
       end do
       ! Apply q_{k}
       rot%cosine=cs_sw(k); rot%sine=ss_sw(k)
       call rotation_times_tbr(trp_rot(rot),b_bv,n,lbw_bv,ubw_bv,0,0,k)
       ! zeros have appeared in the second superdiagonal
       ! in columns k+1, ..., k+ubw_bv-1
       ! or rows k+1-ubw_bv, ..., k-1
       k0=max(k+1,ubw_bv+1)
       k1=min(k+ubw_bv-1,n)
       numrotsu(k)=max(k1-k0+1,0)
       do j=k1,k0,-1
          rot=lgivens2(get_el_br(b_bv,lbw_bv,j-ubw_bv,j), get_el_br(b_bv,lbw_bv,j-ubw_bv+1,j))
          jsu(j-k0+1,k)=j-ubw_bv
          csu(j-k0+1,k)=rot%cosine; ssu(j-k0+1,k)=rot%sine
          call rotation_times_tbr(trp_rot(rot), b_bv,n,lbw_bv,ubw_bv,k,0,j-ubw_bv)
       end do
    end do

    call f_bw_contract_br(b_bv, n, lbw_bv, ubw_bv, lbwmax_bv, ubwmax_bv, &
         dlbw, dlbw_tmp, dubw, dubw_tmp)
    call br_to_bc(b_bv,b_ub,lbw_ub,ubw_ub)

  end subroutine f_d_trp_sweep_times_bv

  subroutine f_d_trp_sweeps1_times_bv(cs_sw, ss_sw, n, numsweeps1_sw, maxsweeps1_sw, &
       b_bv, lbw_bv, ubw_bv, lbwmax_bv, ubwmax_bv, numrotsv, ksv, csv, ssv, &
       b_ub, lbw_ub, ubw_ub, lbwmax_ub, ubwmax_ub, numrotsu, jsu, csu, ssu, error)

    real(kind=dp), dimension(n-1,maxsweeps1_sw), intent(in) :: cs_sw, ss_sw
    integer(kind=int32), intent(in) :: numsweeps1_sw, maxsweeps1_sw, n

    real(kind=dp), dimension(n, lbwmax_bv+ubwmax_bv+1), intent(inout) :: b_bv
    integer(kind=int32), dimension(n), intent(inout) :: numrotsv
    integer(kind=int32), dimension(n,ubwmax_bv), intent(inout) :: ksv
    real(kind=dp), dimension(n,ubwmax_bv), intent(inout) :: csv, ssv

    real(kind=dp), dimension(lbwmax_ub+ubwmax_ub+1,n), intent(out) :: b_ub
    integer(kind=int32), dimension(n), intent(out) :: numrotsu
    integer(kind=int32), dimension(ubwmax_ub,n), intent(out) :: jsu
    real(kind=dp), dimension(ubwmax_ub,n), intent(out) :: csu, ssu

    integer(kind=int32), intent(in) :: lbwmax_ub, ubwmax_ub
    integer(kind=int32), intent(inout) :: lbw_ub, ubw_ub
    integer(kind=int32), intent(in) :: lbwmax_bv, ubwmax_bv
    integer(kind=int32), intent(out) :: lbw_bv, ubw_bv

    type(error_info), intent(out) :: error
    integer(kind=int32) :: l

    call clear_error(error)

    if (n==1) then
       b_ub(1,1)=b_bv(1,1);
       lbw_ub=0; ubw_ub=0; numrotsu=0; return
    end if

    do l=1, numsweeps1_sw
       if (l > 1) then
          call f_d_convert_ub_to_bv(b_ub, n, lbw_ub, ubw_ub, lbwmax_ub, ubwmax_ub, & 
               numrotsu, jsu, csu, ssu, &
               b_bv, lbw_bv, ubw_bv, lbwmax_bv, ubwmax_bv, numrotsv, ksv, csv, ssv, error)
       end if
       call f_d_trp_sweep_times_bv(cs_sw(:,l), ss_sw(:,l), n, numsweeps1_sw, maxsweeps1_sw, &
            b_bv, lbw_bv, ubw_bv, lbwmax_bv, ubwmax_bv, numrotsv, ksv, csv, ssv, &
            b_ub, lbw_ub, ubw_ub, lbwmax_ub, ubwmax_ub, numrotsu, jsu, csu, ssu, error)
    end do
  end subroutine f_d_trp_sweeps1_times_bv

  subroutine c_trp_sweeps1_times_bv(sw,bv,ub,error)
    type(c_sweeps1) :: sw
    type(c_ub) :: ub
    type(c_bv) :: bv
    type(error_info) :: error
    call clear_error(error)
    if (get_n(bv) /= get_n(ub) .or. get_n(bv) /= get_n(sw)) then
       call set_error(error,1,id_c_trp_sweeps1_times_bv); return
    end if
    if (get_ubwmax(bv) < sw%numsweeps1 + ub%ubw+1 .and. get_ubwmax(bv) < get_n(ub)-1) then
       call set_error(error,2,id_c_trp_sweeps1_times_bv); return
    end if
    if (get_ubwmax(ub) < sw%numsweeps1 + ub%ubw+1 .and. get_ubwmax(ub) < get_n(ub)-1) then
       call set_error(error,3,id_c_trp_sweeps1_times_bv); return
    end if
    if (get_n(bv) < 1) then
       call set_error(error,4,id_c_trp_sweeps1_times_bv); return
    end if
    call f_c_trp_sweeps1_times_bv(sw%cs, sw%ss, get_n(sw), sw%numsweeps1, get_maxsweeps1(sw), &
         bv%br, bv%lbw, bv%ubw, get_lbwmax(bv), get_ubwmax(bv), &
         bv%numrotsv, bv%ksv, bv%csv, bv%ssv, &
         ub%bc, ub%lbw, ub%ubw, get_lbwmax(ub), get_ubwmax(ub), &
         ub%numrotsu, ub%jsu, ub%csu, ub%ssu, error)
  end subroutine c_trp_sweeps1_times_bv

  subroutine f_c_trp_sweep_times_bv(cs_sw, ss_sw, n, numsweeps1_sw, maxsweeps1_sw, &
       b_bv, lbw_bv, ubw_bv, lbwmax_bv, ubwmax_bv, numrotsv, ksv, csv, ssv, &
       b_ub, lbw_ub, ubw_ub, lbwmax_ub, ubwmax_ub, numrotsu, jsu, csu, ssu, error)

    complex(kind=dp), dimension(n-1), intent(in) :: cs_sw, ss_sw
    integer(kind=int32), intent(in) :: numsweeps1_sw, maxsweeps1_sw, n

    complex(kind=dp), dimension(n, lbwmax_bv+ubwmax_bv+1), intent(inout) :: b_bv
    integer(kind=int32), dimension(n), intent(inout) :: numrotsv
    integer(kind=int32), dimension(n,ubwmax_bv), intent(inout) :: ksv
    complex(kind=dp), dimension(n,ubwmax_bv), intent(inout) :: csv, ssv

    complex(kind=dp), dimension(lbwmax_ub+ubwmax_ub+1,n), intent(out) :: b_ub
    integer(kind=int32), dimension(n), intent(out) :: numrotsu
    integer(kind=int32), dimension(ubwmax_ub,n), intent(out) :: jsu
    complex(kind=dp), dimension(ubwmax_ub,n), intent(out) :: csu, ssu

    integer(kind=int32), intent(in) :: lbwmax_ub, ubwmax_ub
    integer(kind=int32), intent(inout) :: lbw_ub, ubw_ub
    integer(kind=int32), intent(in) :: lbwmax_bv, ubwmax_bv
    integer(kind=int32), intent(out) :: lbw_bv, ubw_bv

    type(error_info), intent(out) :: error
    integer(kind=int32) :: j, k, k0, k1, dlbw, dlbw_tmp, dubw, dubw_tmp
    type(c_rotation) :: rot

    call clear_error(error)

    if (n==1) then
       b_ub(1,1)=b_bv(1,1);
       lbw_ub=0; ubw_ub=0; numrotsu=0; return
    end if

    dlbw=0; dlbw_tmp=0
    dubw=1; dubw_tmp=1
    
    call f_bw_expand_br(b_bv, n, lbw_bv, ubw_bv, lbwmax_bv, ubwmax_bv, &
         dlbw, dlbw_tmp, dubw, dubw_tmp, lbw_ub, ubw_ub)

    b_ub(1:lbw_ub+ubw_ub+1,:)=(0.0_dp,0.0_dp)
    numrotsu=0
    ssu(1:ubw_ub,:)=(0.0_dp,0.0_dp); csu(1:ubw_ub,:)=(0.0_dp,0.0_dp)
    jsu(1:ubw_ub,:)=0

    do k=1,n-1
       ! apply v_{n-k}
       do j=1,numrotsv(k)
          rot%cosine=csv(k,j); rot%sine=ssv(k,j)
          call tbr_times_rotation(b_bv,n,lbw_bv,ubw_bv,0,n-k,trp_rot(rot),ksv(k,j))
       end do
       ! Apply q_{k}
       rot%cosine=cs_sw(k); rot%sine=ss_sw(k)
       call rotation_times_tbr(trp_rot(rot),b_bv,n,lbw_bv,ubw_bv,0,0,k)
       ! zeros have appeared in the second superdiagonal
       ! in columns k+1, ..., k+ubw_bv-1
       ! or rows k+1-ubw_bv, ..., k-1
       k0=max(k+1,ubw_bv+1)
       k1=min(k+ubw_bv-1,n)
       numrotsu(k)=max(k1-k0+1,0)
       do j=k1,k0,-1
          rot=lgivens2(get_el_br(b_bv,lbw_bv,j-ubw_bv,j), get_el_br(b_bv,lbw_bv,j-ubw_bv+1,j))
          jsu(j-k0+1,k)=j-ubw_bv
          csu(j-k0+1,k)=rot%cosine; ssu(j-k0+1,k)=rot%sine
          call rotation_times_tbr(trp_rot(rot), b_bv,n,lbw_bv,ubw_bv,k,0,j-ubw_bv)
       end do
    end do

    call f_bw_contract_br(b_bv, n, lbw_bv, ubw_bv, lbwmax_bv, ubwmax_bv, &
         dlbw, dlbw_tmp, dubw, dubw_tmp)
    call br_to_bc(b_bv,b_ub,lbw_ub,ubw_ub)

  end subroutine f_c_trp_sweep_times_bv

  subroutine f_c_trp_sweeps1_times_bv(cs_sw, ss_sw, n, numsweeps1_sw, maxsweeps1_sw, &
       b_bv, lbw_bv, ubw_bv, lbwmax_bv, ubwmax_bv, numrotsv, ksv, csv, ssv, &
       b_ub, lbw_ub, ubw_ub, lbwmax_ub, ubwmax_ub, numrotsu, jsu, csu, ssu, error)

    complex(kind=dp), dimension(n-1,maxsweeps1_sw), intent(in) :: cs_sw, ss_sw
    integer(kind=int32), intent(in) :: numsweeps1_sw, maxsweeps1_sw, n

    complex(kind=dp), dimension(n, lbwmax_bv+ubwmax_bv+1), intent(inout) :: b_bv
    integer(kind=int32), dimension(n), intent(inout) :: numrotsv
    integer(kind=int32), dimension(n,ubwmax_bv), intent(inout) :: ksv
    complex(kind=dp), dimension(n,ubwmax_bv), intent(inout) :: csv, ssv

    complex(kind=dp), dimension(lbwmax_ub+ubwmax_ub+1,n), intent(out) :: b_ub
    integer(kind=int32), dimension(n), intent(out) :: numrotsu
    integer(kind=int32), dimension(ubwmax_ub,n), intent(out) :: jsu
    complex(kind=dp), dimension(ubwmax_ub,n), intent(out) :: csu, ssu

    integer(kind=int32), intent(in) :: lbwmax_ub, ubwmax_ub
    integer(kind=int32), intent(inout) :: lbw_ub, ubw_ub
    integer(kind=int32), intent(in) :: lbwmax_bv, ubwmax_bv
    integer(kind=int32), intent(out) :: lbw_bv, ubw_bv

    type(error_info), intent(out) :: error
    integer(kind=int32) :: l

    call clear_error(error)

    if (n==1) then
       b_ub(1,1)=b_bv(1,1);
       lbw_ub=0; ubw_ub=0; numrotsu=0; return
    end if

    do l=1, numsweeps1_sw
       if (l > 1) then
          call f_c_convert_ub_to_bv(b_ub, n, lbw_ub, ubw_ub, lbwmax_ub, ubwmax_ub, & 
               numrotsu, jsu, csu, ssu, &
               b_bv, lbw_bv, ubw_bv, lbwmax_bv, ubwmax_bv, numrotsv, ksv, csv, ssv, error)
       end if
       call f_c_trp_sweep_times_bv(cs_sw(:,l), ss_sw(:,l), n, numsweeps1_sw, maxsweeps1_sw, &
            b_bv, lbw_bv, ubw_bv, lbwmax_bv, ubwmax_bv, numrotsv, ksv, csv, ssv, &
            b_ub, lbw_ub, ubw_ub, lbwmax_ub, ubwmax_ub, numrotsu, jsu, csu, ssu, error)
    end do
  end subroutine f_c_trp_sweeps1_times_bv

  !
  ! errors:
  ! 1: Sizes of bv, ub, and sw do not match.
  ! 2: Insufficient ubw in bv.
  ! 3: Insufficient ubw in ub.
  ! 4: n < 1.


  subroutine d_ub_times_trp_sweeps1(ub,sw,bv,error)
    type(d_sweeps1) :: sw
    type(d_ub) :: ub
    type(d_bv) :: bv
    type(error_info) :: error
    call clear_error(error)
    if (get_n(bv) /= get_n(ub) .or. get_n(bv) /= get_n(sw)) then
       call set_error(error,1,id_d_ub_times_trp_sweeps1); return
    end if
    if (get_ubwmax(bv) < sw%numsweeps1 + ub%ubw+1 .and. get_ubwmax(bv) < get_n(ub)-1) then
       call set_error(error,2,id_d_ub_times_trp_sweeps1); return
    end if
    if (get_ubwmax(ub) < sw%numsweeps1 + ub%ubw+1 .and. get_ubwmax(ub) < get_n(ub)-1) then
       call set_error(error,3,id_d_ub_times_trp_sweeps1); return
    end if
    if (get_n(bv) < 1) then
       call set_error(error,4,id_d_ub_times_trp_sweeps1); return
    end if
    call f_d_ub_times_trp_sweeps1(ub%bc, get_n(ub), ub%lbw, ub%ubw, get_lbwmax(ub), get_ubwmax(ub), &
         ub%numrotsu, ub%jsu, ub%csu, ub%ssu, &
         sw%cs, sw%ss, sw%numsweeps1, get_maxsweeps1(sw), &
         bv%br, bv%lbw, bv%ubw, get_lbwmax(bv), get_ubwmax(bv), &
         bv%numrotsv, bv%ksv, bv%csv, bv%ssv, error)
  end subroutine d_ub_times_trp_sweeps1

  subroutine f_d_ub_times_trp_sweep(b_ub, n, lbw_ub, ubw_ub, lbwmax_ub, ubwmax_ub, &
       numrotsu, jsu, csu, ssu, &
       cs_sw, ss_sw, numsweeps1_sw, maxsweeps1_sw, &
       b_bv, lbw_bv, ubw_bv, lbwmax_bv, ubwmax_bv, numrotsv, ksv, csv, ssv, &
       error)
    real(kind=dp), dimension(n-1), intent(in) :: cs_sw, ss_sw
    integer(kind=int32), intent(in) :: numsweeps1_sw, maxsweeps1_sw, n

    real(kind=dp), dimension(n, lbwmax_bv+ubwmax_bv+1), intent(inout) :: b_bv
    integer(kind=int32), dimension(n), intent(inout) :: numrotsv
    integer(kind=int32), dimension(n,ubwmax_bv), intent(inout) :: ksv
    real(kind=dp), dimension(n,ubwmax_bv), intent(inout) :: csv, ssv

    real(kind=dp), dimension(lbwmax_ub+ubwmax_ub+1,n), intent(out) :: b_ub
    integer(kind=int32), dimension(n), intent(out) :: numrotsu
    integer(kind=int32), dimension(ubwmax_ub,n), intent(out) :: jsu
    real(kind=dp), dimension(ubwmax_ub,n), intent(out) :: csu, ssu

    integer(kind=int32), intent(in) :: lbwmax_ub, ubwmax_ub
    integer(kind=int32), intent(inout) :: lbw_ub, ubw_ub
    integer(kind=int32), intent(in) :: lbwmax_bv, ubwmax_bv
    integer(kind=int32), intent(out) :: lbw_bv, ubw_bv

    type(error_info), intent(out) :: error
    integer(kind=int32) :: j, k, k0, k1, dlbw, dlbw_tmp, dubw, dubw_tmp
    type(d_rotation) :: rot

    call clear_error(error)

    if (n==1) then
       b_bv(1,1)=b_ub(1,1);
       lbw_bv=0; ubw_bv=0; numrotsv=0; return
    end if
    ! Initial expansion for first sweep:
    ! one extra superdiagonal to fill-in.
    ! Expand as need for later sweeps1.

    dlbw=0; dlbw_tmp=0
    dubw=1; dubw_tmp=1
    call f_bw_expand_bc(b_ub, n , lbw_ub, ubw_ub, lbwmax_ub, ubwmax_ub, &
         dlbw, dlbw_tmp, dubw, dubw_tmp, lbw_bv, ubw_bv)

    b_bv(:,1:lbw_bv+ubw_bv+1)=0.0_dp
    numrotsv=0
    ssv(:,1:ubw_bv)=0.0_dp
    csv(:,1:ubw_bv)=0.0_dp
    ksv(:,1:ubw_bv)=0

    do k=1,n-1
       ! Apply u_{n-k}
       do j=1,numrotsu(n-k)
          rot%cosine=csu(j,n-k); rot%sine=ssu(j,n-k)
          call rotation_times_tbc(rot,b_ub,n,lbw_ub,ubw_ub,n-k,0,jsu(j,n-k))
       end do
       rot%cosine=cs_sw(n-k); rot%sine=ss_sw(n-k)
       call tbc_times_rotation(b_ub,n,lbw_ub,ubw_ub,0,0,trp_rot(rot),n-k)
       k0=max(n-k+2,ubw_ub+1)
       k1=min(n-k+ubw_ub-1,n)
       numrotsv(n-k)=max(k1-k0+1,0)
       do j=k0,k1
          rot=rgivens(get_el_bc(b_ub,ubw_ub,j-ubw_ub,j-1), get_el_bc(b_ub,ubw_ub,j-ubw_ub,j))
          ksv(n-k,k1-j+1)=j-1
          csv(n-k,k1-j+1)=rot%cosine; ssv(n-k,k1-j+1)=rot%sine
          call tbc_times_rotation(b_ub,n,lbw_ub,ubw_ub,0,k,rot,j-1)
       end do
    end do

    call f_bw_contract_bc(b_ub, n, lbw_ub, ubw_ub, lbwmax_ub, ubwmax_ub, &
         dlbw, dlbw_tmp, dubw, dubw_tmp)

    call bc_to_br(b_ub,b_bv,lbw_ub,ubw_ub)

  end subroutine f_d_ub_times_trp_sweep

  subroutine f_d_ub_times_trp_sweeps1(b_ub, n, lbw_ub, ubw_ub, lbwmax_ub, ubwmax_ub, &
       numrotsu, jsu, csu, ssu, &
       cs_sw, ss_sw, numsweeps1_sw, maxsweeps1_sw, &
       b_bv, lbw_bv, ubw_bv, lbwmax_bv, ubwmax_bv, numrotsv, ksv, csv, ssv, &
       error)
    real(kind=dp), dimension(n-1,maxsweeps1_sw), intent(in) :: cs_sw, ss_sw
    integer(kind=int32), intent(in) :: numsweeps1_sw, maxsweeps1_sw, n

    real(kind=dp), dimension(n, lbwmax_bv+ubwmax_bv+1), intent(inout) :: b_bv
    integer(kind=int32), dimension(n), intent(inout) :: numrotsv
    integer(kind=int32), dimension(n,ubwmax_bv), intent(inout) :: ksv
    real(kind=dp), dimension(n,ubwmax_bv), intent(inout) :: csv, ssv

    real(kind=dp), dimension(lbwmax_ub+ubwmax_ub+1,n), intent(out) :: b_ub
    integer(kind=int32), dimension(n), intent(out) :: numrotsu
    integer(kind=int32), dimension(ubwmax_ub,n), intent(out) :: jsu
    real(kind=dp), dimension(ubwmax_ub,n), intent(out) :: csu, ssu

    integer(kind=int32), intent(in) :: lbwmax_ub, ubwmax_ub
    integer(kind=int32), intent(inout) :: lbw_ub, ubw_ub
    integer(kind=int32), intent(in) :: lbwmax_bv, ubwmax_bv
    integer(kind=int32), intent(out) :: lbw_bv, ubw_bv

    type(error_info), intent(out) :: error
    integer(kind=int32) :: l

    call clear_error(error)

    if (n==1) then
       b_bv(1,1)=b_ub(1,1);
       lbw_bv=0; ubw_bv=0; numrotsv=0; return
    end if

    do l=numsweeps1_sw,1,-1
       if (l < numsweeps1_sw) then
          call f_d_convert_bv_to_ub(b_bv, n, lbw_bv, ubw_bv, lbwmax_bv, ubwmax_bv, numrotsv, &
               ksv, csv, ssv, b_ub, lbw_ub, ubw_ub, lbwmax_ub, ubwmax_ub, & 
               numrotsu, jsu, csu, ssu, error)
       end if
       call f_d_ub_times_trp_sweep(b_ub, n, lbw_ub, ubw_ub, lbwmax_ub, ubwmax_ub, &
            numrotsu, jsu, csu, ssu, &
            cs_sw(:,l), ss_sw(:,l), numsweeps1_sw, maxsweeps1_sw, &
            b_bv, lbw_bv, ubw_bv, lbwmax_bv, ubwmax_bv, numrotsv, ksv, csv, ssv, &
            error)
    end do
  end subroutine f_d_ub_times_trp_sweeps1

  subroutine c_ub_times_trp_sweeps1(ub,sw,bv,error)
    type(c_sweeps1) :: sw
    type(c_ub) :: ub
    type(c_bv) :: bv
    type(error_info) :: error
    call clear_error(error)
    if (get_n(bv) /= get_n(ub) .or. get_n(bv) /= get_n(sw)) then
       call set_error(error,1,id_c_ub_times_trp_sweeps1); return
    end if
    if (get_ubwmax(bv) < sw%numsweeps1 + ub%ubw+1 .and. get_ubwmax(bv) < get_n(ub)-1) then
       call set_error(error,2,id_c_ub_times_trp_sweeps1); return
    end if
    if (get_ubwmax(ub) < sw%numsweeps1 + ub%ubw+1 .and. get_ubwmax(ub) < get_n(ub)-1) then
       call set_error(error,3,id_c_ub_times_trp_sweeps1); return
    end if
    if (get_n(bv) < 1) then
       call set_error(error,4,id_c_ub_times_trp_sweeps1); return
    end if
    call f_c_ub_times_trp_sweeps1(ub%bc, get_n(ub), ub%lbw, ub%ubw, get_lbwmax(ub), get_ubwmax(ub), &
         ub%numrotsu, ub%jsu, ub%csu, ub%ssu, &
         sw%cs, sw%ss, sw%numsweeps1, get_maxsweeps1(sw), &
         bv%br, bv%lbw, bv%ubw, get_lbwmax(bv), get_ubwmax(bv), &
         bv%numrotsv, bv%ksv, bv%csv, bv%ssv, error)
  end subroutine c_ub_times_trp_sweeps1

  subroutine f_c_ub_times_trp_sweep(b_ub, n, lbw_ub, ubw_ub, lbwmax_ub, ubwmax_ub, &
       numrotsu, jsu, csu, ssu, &
       cs_sw, ss_sw, numsweeps1_sw, maxsweeps1_sw, &
       b_bv, lbw_bv, ubw_bv, lbwmax_bv, ubwmax_bv, numrotsv, ksv, csv, ssv, &
       error)
    complex(kind=dp), dimension(n-1), intent(in) :: cs_sw, ss_sw
    integer(kind=int32), intent(in) :: numsweeps1_sw, maxsweeps1_sw, n

    complex(kind=dp), dimension(n, lbwmax_bv+ubwmax_bv+1), intent(inout) :: b_bv
    integer(kind=int32), dimension(n), intent(inout) :: numrotsv
    integer(kind=int32), dimension(n,ubwmax_bv), intent(inout) :: ksv
    complex(kind=dp), dimension(n,ubwmax_bv), intent(inout) :: csv, ssv

    complex(kind=dp), dimension(lbwmax_ub+ubwmax_ub+1,n), intent(out) :: b_ub
    integer(kind=int32), dimension(n), intent(out) :: numrotsu
    integer(kind=int32), dimension(ubwmax_ub,n), intent(out) :: jsu
    complex(kind=dp), dimension(ubwmax_ub,n), intent(out) :: csu, ssu

    integer(kind=int32), intent(in) :: lbwmax_ub, ubwmax_ub
    integer(kind=int32), intent(inout) :: lbw_ub, ubw_ub
    integer(kind=int32), intent(in) :: lbwmax_bv, ubwmax_bv
    integer(kind=int32), intent(out) :: lbw_bv, ubw_bv

    type(error_info), intent(out) :: error
    integer(kind=int32) :: j, k, k0, k1, dlbw, dlbw_tmp, dubw, dubw_tmp
    type(c_rotation) :: rot

    call clear_error(error)

    if (n==1) then
       b_bv(1,1)=b_ub(1,1);
       lbw_bv=0; ubw_bv=0; numrotsv=0; return
    end if
    ! Initial expansion for first sweep:
    ! one extra superdiagonal to fill-in.
    ! Expand as need for later sweeps1.

    dlbw=0; dlbw_tmp=0
    dubw=1; dubw_tmp=1
    call f_bw_expand_bc(b_ub, n , lbw_ub, ubw_ub, lbwmax_ub, ubwmax_ub, &
         dlbw, dlbw_tmp, dubw, dubw_tmp, lbw_bv, ubw_bv)

    b_bv(:,1:lbw_bv+ubw_bv+1)=(0.0_dp,0.0_dp)
    numrotsv=0
    ssv(:,1:ubw_bv)=(0.0_dp,0.0_dp)
    csv(:,1:ubw_bv)=(0.0_dp,0.0_dp)
    ksv(:,1:ubw_bv)=0

    do k=1,n-1
       ! Apply u_{n-k}
       do j=1,numrotsu(n-k)
          rot%cosine=csu(j,n-k); rot%sine=ssu(j,n-k)
          call rotation_times_tbc(rot,b_ub,n,lbw_ub,ubw_ub,n-k,0,jsu(j,n-k))
       end do
       rot%cosine=cs_sw(n-k); rot%sine=ss_sw(n-k)
       call tbc_times_rotation(b_ub,n,lbw_ub,ubw_ub,0,0,trp_rot(rot),n-k)
       k0=max(n-k+2,ubw_ub+1)
       k1=min(n-k+ubw_ub-1,n)
       numrotsv(n-k)=max(k1-k0+1,0)
       do j=k0,k1
          rot=rgivens(get_el_bc(b_ub,ubw_ub,j-ubw_ub,j-1), get_el_bc(b_ub,ubw_ub,j-ubw_ub,j))
          ksv(n-k,k1-j+1)=j-1
          csv(n-k,k1-j+1)=rot%cosine; ssv(n-k,k1-j+1)=rot%sine
          call tbc_times_rotation(b_ub,n,lbw_ub,ubw_ub,0,k,rot,j-1)
       end do
    end do

    call f_bw_contract_bc(b_ub, n, lbw_ub, ubw_ub, lbwmax_ub, ubwmax_ub, &
         dlbw, dlbw_tmp, dubw, dubw_tmp)

    call bc_to_br(b_ub,b_bv,lbw_ub,ubw_ub)

  end subroutine f_c_ub_times_trp_sweep

  subroutine f_c_ub_times_trp_sweeps1(b_ub, n, lbw_ub, ubw_ub, lbwmax_ub, ubwmax_ub, &
       numrotsu, jsu, csu, ssu, &
       cs_sw, ss_sw, numsweeps1_sw, maxsweeps1_sw, &
       b_bv, lbw_bv, ubw_bv, lbwmax_bv, ubwmax_bv, numrotsv, ksv, csv, ssv, &
       error)
    complex(kind=dp), dimension(n-1,maxsweeps1_sw), intent(in) :: cs_sw, ss_sw
    integer(kind=int32), intent(in) :: numsweeps1_sw, maxsweeps1_sw, n

    complex(kind=dp), dimension(n, lbwmax_bv+ubwmax_bv+1), intent(inout) :: b_bv
    integer(kind=int32), dimension(n), intent(inout) :: numrotsv
    integer(kind=int32), dimension(n,ubwmax_bv), intent(inout) :: ksv
    complex(kind=dp), dimension(n,ubwmax_bv), intent(inout) :: csv, ssv

    complex(kind=dp), dimension(lbwmax_ub+ubwmax_ub+1,n), intent(out) :: b_ub
    integer(kind=int32), dimension(n), intent(out) :: numrotsu
    integer(kind=int32), dimension(ubwmax_ub,n), intent(out) :: jsu
    complex(kind=dp), dimension(ubwmax_ub,n), intent(out) :: csu, ssu

    integer(kind=int32), intent(in) :: lbwmax_ub, ubwmax_ub
    integer(kind=int32), intent(inout) :: lbw_ub, ubw_ub
    integer(kind=int32), intent(in) :: lbwmax_bv, ubwmax_bv
    integer(kind=int32), intent(out) :: lbw_bv, ubw_bv

    type(error_info), intent(out) :: error
    integer(kind=int32) :: l

    call clear_error(error)

    if (n==1) then
       b_bv(1,1)=b_ub(1,1);
       lbw_bv=0; ubw_bv=0; numrotsv=0; return
    end if

    do l=numsweeps1_sw,1,-1
       if (l < numsweeps1_sw) then
          call f_c_convert_bv_to_ub(b_bv, n, lbw_bv, ubw_bv, lbwmax_bv, ubwmax_bv, numrotsv, &
               ksv, csv, ssv, b_ub, lbw_ub, ubw_ub, lbwmax_ub, ubwmax_ub, & 
               numrotsu, jsu, csu, ssu, error)
       end if
       call f_c_ub_times_trp_sweep(b_ub, n, lbw_ub, ubw_ub, lbwmax_ub, ubwmax_ub, &
            numrotsu, jsu, csu, ssu, &
            cs_sw(:,l), ss_sw(:,l), numsweeps1_sw, maxsweeps1_sw, &
            b_bv, lbw_bv, ubw_bv, lbwmax_bv, ubwmax_bv, numrotsv, ksv, csv, ssv, &
            error)
    end do
  end subroutine f_c_ub_times_trp_sweeps1

  type(d_sweeps1) function d_random_sweeps1(n,l) result(sw)
    integer(kind=int32) :: n, l, j, k
    real(kind=dp) :: nrm
    sw=d_new_sweeps1(n,l)
    sw%numsweeps1=l
    call random_matrix(sw%cs)
    call random_matrix(sw%ss)
    do j=1,n-1
       do k=1,l
          nrm=sqrt(sw%cs(j,k)**2+sw%ss(j,k)**2)
          sw%cs(j,k)=sw%cs(j,k)/nrm
          sw%ss(j,k)=sw%ss(j,k)/nrm
       end do
    end do
  end function d_random_sweeps1

  type(c_sweeps1) function c_random_sweeps1(n,l) result(sw)
    integer(kind=int32) :: n, l, j, k
    real(kind=dp) :: nrm, cr
    sw=c_new_sweeps1(n,l)
    sw%numsweeps1=l
    call random_matrix(sw%cs)
    call random_matrix(sw%ss)
    do j=1,n-1
       do k=1,l
          cr=real(sw%cs(j,k),kind=dp)
          sw%cs(j,k)=cmplx(cr,0.0_dp)
          nrm=sqrt(abs(sw%cs(j,k))**2+abs(sw%ss(j,k))**2)
          sw%cs(j,k)=sw%cs(j,k)/nrm
          sw%ss(j,k)=sw%ss(j,k)/nrm
       end do
    end do
  end function c_random_sweeps1

end module sweeps1
