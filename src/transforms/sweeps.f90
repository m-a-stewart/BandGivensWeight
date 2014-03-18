module sweeps
  use misc
  use types
  use shift
  use conversions

  implicit none

  !
  ! These types represent a linear transformation
  ! Q = Q_1 Q_2 ... Q_{numsweeps}
  ! where
  ! Q_k = G_{k,1} ... G_{k,n-1}
  ! Thus Q is a product of upper Hessenberg matrices.


  type d_sweeps
     integer(kind=int32), private :: maxsweeps, n
     logical :: transposed
     integer(kind=int32) :: numsweeps
     real(kind=dp), allocatable, dimension(:,:) :: cs, ss
  end type d_sweeps

  type c_sweeps
     integer(kind=int32), private :: maxsweeps, n
     logical :: transposed
     integer(kind=int32) :: numsweeps
     complex(kind=dp), allocatable, dimension(:,:) :: cs, ss
  end type c_sweeps

  interface deallocate_sweeps
     module procedure d_deallocate_sweeps, c_deallocate_sweeps
  end interface deallocate_sweeps

  interface get_maxsweeps
     module procedure d_get_maxsweeps, c_get_maxsweeps
  end interface get_maxsweeps

  interface get_n
     module procedure d_get_n_sweeps, c_get_n_sweeps
  end interface get_n

  interface trp_sweeps
     module procedure d_trp_sweeps, c_trp_sweeps
  end interface trp_sweeps

  interface sweeps_times_general
     module procedure d_sweeps_times_general, c_sweeps_times_general, &
          d_v_sweeps_times_general, c_v_sweeps_times_general
  end interface sweeps_times_general

  interface general_times_sweeps
     module procedure d_general_times_sweeps, c_general_times_sweeps, &
          d_v_general_times_sweeps, c_v_general_times_sweeps
  end interface general_times_sweeps

  interface sweeps_times_ub
     module procedure d_sweeps_times_ub, c_sweeps_times_ub
  end interface sweeps_times_ub

  interface f_sweeps_times_ub
     module procedure f_d_sweeps_times_ub, f_c_sweeps_times_ub
  end interface f_sweeps_times_ub

  type(routine_info), parameter :: info_d_sweeps_times_ub= &
       routine_info(id_d_sweeps_times_ub, &
       'd_sweeps_times_ub', &
       [ character(len=error_message_length) :: 'Sizes of bv, ub, and sw do not match.', &
       'Insufficienty lbw in bv.', 'Insufficient ubw in bv.', 'Insufficient lbw in ub', &
       'Insufficient ubw in ub', 'n<1', 'sw is transposed' ])

  type(routine_info), parameter :: info_c_sweeps_times_ub= &
       routine_info(id_c_sweeps_times_ub, &
       'c_sweeps_times_ub', &
       [ character(len=error_message_length) :: 'Sizes of bv, ub, and sw do not match.', &
       'Insufficienty lbw in bv.', 'Insufficient ubw in bv.', 'Insufficient lbw in ub', &
       'Insufficient ubw in ub', 'n<1', 'sw is transposed' ])

  type(routine_info), parameter :: info_f_d_sweeps_times_ub= &
       routine_info(id_f_d_sweeps_times_ub, &
       'f_d_sweeps_times_ub', &
       [ character(len=error_message_length) :: '' ])

  type(routine_info), parameter :: info_f_c_sweeps_times_ub= &
       routine_info(id_f_c_sweeps_times_ub, &
       'f_c_sweeps_times_ub', &
       [ character(len=error_message_length) :: '' ])


contains

  integer(kind=int32) function d_get_maxsweeps(sw) result(maxsweeps)
    type(d_sweeps) :: sw
    maxsweeps=sw%maxsweeps
  end function d_get_maxsweeps

  integer(kind=int32) function c_get_maxsweeps(sw) result(maxsweeps)
    type(c_sweeps) :: sw
    maxsweeps=sw%maxsweeps
  end function c_get_maxsweeps

  integer(kind=int32) function d_get_n_sweeps(sw) result(n)
    type(d_sweeps) :: sw
    n=sw%n
  end function d_get_n_sweeps

  integer(kind=int32) function c_get_n_sweeps(sw) result(n)
    type(c_sweeps) :: sw
    n=sw%n
  end function c_get_n_sweeps

  type(d_sweeps) function d_new_sweeps(n,maxsweeps) result(sw)
    integer(kind=int32), intent(in) :: n, maxsweeps
    sw%n=n; sw%maxsweeps=maxsweeps
    allocate(sw%cs(n,maxsweeps), sw%ss(n,maxsweeps))
    sw%transposed=.false.
    sw%numsweeps=0
    sw%cs=1.0_dp; sw%ss=0.0_dp
  end function d_new_sweeps

  type(c_sweeps) function c_new_sweeps(n,maxsweeps) result(sw)
    integer(kind=int32), intent(in) :: n, maxsweeps
    sw%n=n; sw%maxsweeps=maxsweeps
    allocate(sw%cs(n,maxsweeps), sw%ss(n,maxsweeps))
    sw%transposed=.false.
    sw%numsweeps=0
    sw%cs=(1.0_dp,0.0_dp); sw%ss=(0.0_dp,0.0_dp)
  end function c_new_sweeps

  subroutine d_deallocate_sweeps(sw)
    type(d_sweeps), intent(inout) :: sw
    deallocate(sw%cs, sw%ss)
  end subroutine d_deallocate_sweeps

  subroutine c_deallocate_sweeps(sw)
    type(c_sweeps), intent(inout) :: sw
    deallocate(sw%cs, sw%ss)
  end subroutine c_deallocate_sweeps

  subroutine d_trp_sweeps(sw)
    type(d_sweeps), intent(inout) :: sw
    if (sw%transposed) then
       sw%transposed=.false.
    else
       sw%transposed=.true.
    end if
  end subroutine d_trp_sweeps

  subroutine c_trp_sweeps(sw)
    type(c_sweeps), intent(inout) :: sw
    if (sw%transposed) then
       sw%transposed=.false.
    else
       sw%transposed=.true.
    end if
  end subroutine c_trp_sweeps

  subroutine d_sweeps_times_general(sw,a)
    type(d_sweeps) :: sw
    real(kind=dp), dimension(:,:), intent(inout) :: a
    type(d_rotation) :: rot
    integer(kind=int32) :: j, m, l
    if (sw%transposed) then
       do l=1,sw%numsweeps
          m=size(a,1)
          do j=1,m-1
             rot%cosine=sw%cs(j,l); rot%sine=sw%ss(j,l)
             call d_rotation_times_general(trp_rot(rot),a,j,j+1)
          end do
       end do
    else
       do l=sw%numsweeps,1,-1
          m=size(a,1)
          do j=m-1,1,-1
             rot%cosine=sw%cs(j,l); rot%sine=sw%ss(j,l)
             call d_rotation_times_general(rot,a,j,j+1)
          end do
       end do
    end if
  end subroutine d_sweeps_times_general

  subroutine d_v_sweeps_times_general(sw,a)
    type(d_sweeps) :: sw
    real(kind=dp), dimension(:), intent(inout) :: a
    integer(kind=int32) :: j, m, l
    real(kind=dp) :: tmp, c, s 
    if (sw%transposed) then
       do l=1,sw%numsweeps
          m=size(a)
          do j=1,m-1
             tmp=a(j)
             c=sw%cs(j,l); s=sw%ss(j,l)
             a(j)=c*tmp+s*a(j+1)
             a(j+1)=-s*tmp+c*a(j+1)
          end do
       end do
    else
       do l=sw%numsweeps,1,-1
          m=size(a)
          do j=m-1,1,-1
             tmp=a(j)
             c=sw%cs(j,l); s=sw%ss(j,l)
             a(j)=c*tmp-s*a(j+1)
             a(j+1)=s*tmp+c*a(j+1)
          end do
       end do
    end if
  end subroutine d_v_sweeps_times_general


  subroutine c_sweeps_times_general(sw,a)
    type(c_sweeps) :: sw
    complex(kind=dp), dimension(:,:), intent(inout) :: a
    type(c_rotation) :: rot
    integer(kind=int32) :: j, m, l
    if (sw%transposed) then
       do l=1,sw%numsweeps
          m=size(a,1)
          do j=1,m-1
             rot%cosine=sw%cs(j,l); rot%sine=sw%ss(j,l)
             call c_rotation_times_general(trp_rot(rot),a,j,j+1)
          end do
       end do
    else
       do l=sw%numsweeps,1,-1
          m=size(a,1)
          do j=m-1,1,-1
             rot%cosine=sw%cs(j,l); rot%sine=sw%ss(j,l)
             call c_rotation_times_general(rot,a,j,j+1)
          end do
       end do
    end if
  end subroutine c_sweeps_times_general

  subroutine c_v_sweeps_times_general(sw,a)
    type(c_sweeps) :: sw
    complex(kind=dp), dimension(:), intent(inout) :: a
    integer(kind=int32) :: j, m, l
    complex(kind=dp) :: tmp, c, s 
    if (sw%transposed) then
       do l=1,sw%numsweeps
          m=size(a)
          do j=1,m-1
             tmp=a(j)
             c=sw%cs(j,l); s=sw%ss(j,l)
             a(j)=c*tmp+conjg(s)*a(j+1)
             a(j+1)=-s*tmp+c*a(j+1)
          end do
       end do
    else
       do l=sw%numsweeps,1,-1
          m=size(a)
          do j=m-1,1,-1
             tmp=a(j)
             c=sw%cs(j,l); s=sw%ss(j,l)
             a(j)=c*tmp-conjg(s)*a(j+1)
             a(j+1)=s*tmp+c*a(j+1)
          end do
       end do
    end if
  end subroutine c_v_sweeps_times_general

  subroutine d_general_times_sweeps(a,sw)
    type(d_sweeps) :: sw
    real(kind=dp), dimension(:,:), intent(inout) :: a
    type(d_rotation) :: rot
    integer(kind=int32) :: j, n, l
    if (sw%transposed) then
       do l=sw%numsweeps,1,-1
          n=size(a,2)
          do j=n-1,1,-1
             rot%cosine=sw%cs(j,l); rot%sine=sw%ss(j,l)
             call d_general_times_rotation(a, trp_rot(rot),j,j+1)
          end do
       end do
    else
       do l=1,sw%numsweeps
          n=size(a,2)
          do j=1,n-1
             rot%cosine=sw%cs(j,l); rot%sine=sw%ss(j,l)
             call d_general_times_rotation(a,rot,j,j+1)
          end do
       end do
    end if
  end subroutine d_general_times_sweeps

  subroutine d_v_general_times_sweeps(a,sw)
    type(d_sweeps) :: sw
    real(kind=dp), dimension(:), intent(inout) :: a
    integer(kind=int32) :: j, n, l
    real(kind=dp) :: c,s,tmp
    if (sw%transposed) then
       do l=sw%numsweeps,1,-1
          n=size(a)
          do j=n-1,1,-1
             tmp=a(j)
             c=sw%cs(j,l); s=sw%ss(j,l)
             a(j)=c*tmp-s*a(j+1)
             a(j+1)=s*tmp+c*a(j+1)
          end do
       end do
    else
       do l=1,sw%numsweeps
          n=size(a)
          do j=1,n-1
             tmp=a(j)
             c=sw%cs(j,l); s=sw%ss(j,l)
             a(j)=c*tmp+s*a(j+1)
             a(j+1)=-s*tmp+c*a(j+1)
          end do
       end do
    end if
  end subroutine d_v_general_times_sweeps


  subroutine c_general_times_sweeps(a,sw)
    type(c_sweeps) :: sw
    complex(kind=dp), dimension(:,:), intent(inout) :: a
    type(c_rotation) :: rot
    integer(kind=int32) :: j, n, l
    if (sw%transposed) then
       do l=sw%numsweeps,1,-1
          n=size(a,2)
          do j=n-1,1,-1
             rot%cosine=sw%cs(j,l); rot%sine=sw%ss(j,l)
             call c_general_times_rotation(a, trp_rot(rot),j,j+1)
          end do
       end do
    else
       do l=1,sw%numsweeps
          n=size(a,2)
          do j=1,n-1
             rot%cosine=sw%cs(j,l); rot%sine=sw%ss(j,l)
             call c_general_times_rotation(a,rot,j,j+1)
          end do
       end do
    end if
  end subroutine c_general_times_sweeps

  subroutine c_v_general_times_sweeps(a,sw)
    type(c_sweeps) :: sw
    complex(kind=dp), dimension(:), intent(inout) :: a
    integer(kind=int32) :: j, n, l
    complex(kind=dp) :: c,s,tmp
    if (sw%transposed) then
       do l=sw%numsweeps,1,-1
          n=size(a)
          do j=n-1,1,-1
             tmp=a(j)
             c=sw%cs(j,l); s=sw%ss(j,l)
             a(j)=c*tmp-s*a(j+1)
             a(j+1)=conjg(s)*tmp+c*a(j+1)
          end do
       end do
    else
       do l=1,sw%numsweeps
          n=size(a)
          do j=1,n-1
             tmp=a(j)
             c=sw%cs(j,l); s=sw%ss(j,l)
             a(j)=c*tmp+s*a(j+1)
             a(j+1)=-conjg(s)*tmp+c*a(j+1)
          end do
       end do
    end if
  end subroutine c_v_general_times_sweeps

  ! sweeps times ub
  ! Errors:
  ! 1: Sizes of bv, ub, and sw do not match.
  ! 2: Insufficient lbw in bv.
  ! 3: Insufficient ubw in bv.
  ! 4: Insufficient lbw in ub
  ! 5: Insufficient ubw in ub.
  ! 6: n < 1.
  ! 7: sw is transposed.

  subroutine d_sweeps_times_ub(sw,ub,bv,error)
    type(d_sweeps) :: sw
    type(d_ub) :: ub
    type(d_bv) :: bv
    type(error_info) :: error
    call clear_error(error)
    if (get_n(bv) /= get_n(ub) .or. get_n(bv) /= get_n(sw)) then
       call set_error(error,1,id_d_sweeps_times_ub); return
    end if
    if (get_lbwmax(bv) < sw%numsweeps + ub%lbw .and. get_lbwmax(bv) < get_n(ub)-1) then
       call set_error(error,2,id_d_sweeps_times_ub); return
    end if
    if (get_ubwmax(bv) < sw%numsweeps + ub%ubw+1 .and. get_ubwmax(bv) < get_n(ub)-1) then
       call set_error(error,3,id_d_sweeps_times_ub); return
    end if
    if (get_lbwmax(ub) < sw%numsweeps + ub%lbw .and. get_lbwmax(ub) < get_n(ub)-1) then
       call set_error(error,4,id_d_sweeps_times_ub); return
    end if
    if (get_ubwmax(ub) < sw%numsweeps + ub%ubw+2 .and. get_ubwmax(ub) < get_n(ub)+1) then
       call set_error(error,5,id_d_sweeps_times_ub); return
    end if
    if (get_n(bv) < 1) then
       call set_error(error,6,id_d_sweeps_times_ub); return
    end if
    if (sw%transposed) then
       call set_error(error,7,id_d_sweeps_times_ub); return
    end if
    call f_d_sweeps_times_ub(sw%cs, sw%ss, sw%numsweeps, get_maxsweeps(sw), get_n(sw), &
         ub%b, ub%lbw, ub%ubw, get_lbwmax(ub), get_ubwmax(ub), &
         ub%numrotsu, ub%jsu, ub%csu, ub%ssu, &
         bv%b, bv%lbw, bv%ubw, get_lbwmax(bv), get_ubwmax(bv), &
         bv%numrotsv, bv%ksv, bv%csv, bv%ssv, error)
  end subroutine d_sweeps_times_ub

  subroutine f_d_sweeps_times_ub(cs_sw, ss_sw, numsweeps_sw, maxsweeps_sw, n, &
       b_ub, lbw_ub, ubw_ub, lbwmax_ub, ubwmax_ub, numrots_ub, js_ub, cs_ub, ss_ub, &
       b_bv, lbw_bv, ubw_bv, lbwmax_bv, ubwmax_bv, numrots_bv, ks_bv, cs_bv, ss_bv, error)

    real(kind=dp), dimension(n,maxsweeps_sw), intent(in) :: cs_sw, ss_sw
    integer(kind=int32), intent(in) :: numsweeps_sw, maxsweeps_sw, n

    real(kind=dp), dimension(lbwmax_ub+ubwmax_ub+1,n), intent(inout) :: b_ub
    integer(kind=int32), dimension(n), intent(inout) :: numrots_ub
    integer(kind=int32), dimension(ubwmax_ub,n), intent(inout) :: js_ub
    real(kind=dp), dimension(ubwmax_ub,n), intent(inout) :: cs_ub, ss_ub

    real(kind=dp), dimension(n, lbwmax_bv+ubwmax_bv+1), intent(out) :: b_bv
    integer(kind=int32), dimension(n), intent(out) :: numrots_bv
    integer(kind=int32), dimension(n,ubwmax_bv), intent(out) :: ks_bv
    real(kind=dp), dimension(n,ubwmax_bv), intent(out) :: cs_bv, ss_bv

    integer(kind=int32), intent(in) :: lbwmax_ub, ubwmax_ub
    integer(kind=int32), intent(inout) :: lbw_ub, ubw_ub
    integer(kind=int32), intent(in) :: lbwmax_bv, ubwmax_bv
    integer(kind=int32), intent(out) :: lbw_bv, ubw_bv

    type(error_info), intent(out) :: error
    integer(kind=int32) :: j, k, l, ubw, lbw, k0, k1, d
    type(d_rotation) :: rot
    logical :: full_ubw, full_lbw

    call clear_error(error)

    if (n==1) then
       b_bv(1,1)=b_ub(1,1); return
    end if
    ! Initial expansion for first sweep: one extra subdiagonal and
    ! one extra superdiagonal to fill-in.
    ! Expand as need for later sweeps.
    lbw=min(lbw_ub+numsweeps_sw,n-1)
    ubw=min(ubw_ub+numsweeps_sw,n-1)
    b_bv(:,1:lbw+ubw+1)=0.0_dp
    numrots_bv=0
    ss_bv(:,1:ubw)=0.0_dp
    cs_bv(:,1:ubw)=0.0_dp
    ks_bv(:,1:ubw)=0
    ubw=ubw_ub
    lbw=lbw_ub

    full_lbw=(lbw==n-1)
    full_ubw=(ubw==n-1)

    do l=numsweeps_sw,1,-1
       if (l < numsweeps_sw) then
          if (.not. full_ubw) then
             ubw=ubw-1
          end if
          call f_d_convert_bv_to_ub(b_bv, n, lbw, ubw, lbwmax_bv, ubwmax_bv, &
               numrots_bv, ks_bv, cs_bv, ss_bv, b_ub, lbw_ub, ubw_ub, lbwmax_ub, ubwmax_ub, & 
               numrots_ub, js_ub, cs_ub, ss_ub, error)
       end if
       numrots_bv=0
       ! Unless all band structure is filled, sweep will fill in one sub and one super diagonal.
       if (.not. full_ubw) then
          ! ubw==n-3 results in ubw=n-1, but this is only due to the temporary diagonal.
          ! So this case does not set full_ubw to .true.  It will be set on the next sweep.
          if (ubw <= n-3) then  
             ubw=ubw+2
             call down_shift(b_ub)
             call down_shift(b_ub)
          else if (ubw==n-2) then
             ubw=n-1
             call down_shift(b_ub)
             full_ubw=.true.
          end if
       end if
       if (.not. full_lbw) then
          lbw=lbw+1;
          full_lbw=(lbw==n-1)
          b_ub(lbw+ubw+1,:)=0.0_dp
       end if
       do k=1,n-1
          do j=1,numrots_ub(n-k)
             rot%cosine=cs_ub(j,n-k); rot%sine=ss_ub(j,n-k)
             call rotation_times_tbc(rot,b_ub,n,lbw,ubw,n-k,0,js_ub(j,n-k))
          end do
          ! Apply q_{n-k}
          rot%cosine=cs_sw(n-k,l); rot%sine=ss_sw(n-k,l)
          call rotation_times_tbc(rot,b_ub,n,lbw,ubw,0,0,n-k)
          ! zeros have appeared in the second superdiagonal in
          ! rows n-k-(ubw-2), ... , n-k or
          ! columns n-k+2, ..., n-k+ubw
          k0=max(n-k+2,ubw+1)
          k1=min(n-k+ubw,n)
          numrots_bv(k)=max(k1-k0+1,0)
          do j=k0,k1
             rot=rgivens(get_el_bc(b_ub,ubw,j-ubw,j-1), get_el_bc(b_ub,ubw,j-ubw,j))
             ks_bv(k,k1-j+1)=j-1
             cs_bv(k,k1-j+1)=rot%cosine; ss_bv(k,k1-j+1)=rot%sine
             call tbc_times_rotation(b_ub,n,lbw,ubw,0,k,rot,j-1)
          end do
       end do
       call bc_to_br(b_ub,b_bv,lbw,ubw)
    end do
    lbw_bv=lbw
    if (.not. full_ubw) then
       ubw_bv=ubw-1
    else
       ubw_bv=ubw
    end if
  end subroutine f_d_sweeps_times_ub

  subroutine c_sweeps_times_ub(sw,ub,bv,error)
    type(c_sweeps) :: sw
    type(c_ub) :: ub
    type(c_bv) :: bv
    type(error_info) :: error
    call clear_error(error)
    if (get_n(bv) /= get_n(ub) .or. get_n(bv) /= get_n(sw)) then
       call set_error(error,1,id_c_sweeps_times_ub); return
    end if
    if (get_lbwmax(bv) < sw%numsweeps + ub%lbw .and. get_lbwmax(bv) < get_n(ub)-1) then
       call set_error(error,2,id_c_sweeps_times_ub); return
    end if
    if (get_ubwmax(bv) < sw%numsweeps + ub%ubw+1 .and. get_ubwmax(bv) < get_n(ub)-1) then
       call set_error(error,3,id_c_sweeps_times_ub); return
    end if
    if (get_lbwmax(ub) < sw%numsweeps + ub%lbw .and. get_lbwmax(ub) < get_n(ub)-1) then
       call set_error(error,4,id_c_sweeps_times_ub); return
    end if
    if (get_ubwmax(ub) < sw%numsweeps + ub%ubw+2 .and. get_ubwmax(ub) < get_n(ub)+1) then
       call set_error(error,5,id_c_sweeps_times_ub); return
    end if
    if (get_n(bv) < 1) then
       call set_error(error,6,id_c_sweeps_times_ub); return
    end if
    if (sw%transposed) then
       call set_error(error,7,id_c_sweeps_times_ub); return
    end if
    call f_c_sweeps_times_ub(sw%cs, sw%ss, sw%numsweeps, get_maxsweeps(sw), get_n(sw), &
         ub%b, ub%lbw, ub%ubw, get_lbwmax(ub), get_ubwmax(ub), &
         ub%numrotsu, ub%jsu, ub%csu, ub%ssu, &
         bv%b, bv%lbw, bv%ubw, get_lbwmax(bv), get_ubwmax(bv), &
         bv%numrotsv, bv%ksv, bv%csv, bv%ssv, error)
  end subroutine c_sweeps_times_ub

  subroutine f_c_sweeps_times_ub(cs_sw, ss_sw, numsweeps_sw, maxsweeps_sw, n, &
       b_ub, lbw_ub, ubw_ub, lbwmax_ub, ubwmax_ub, numrots_ub, js_ub, cs_ub, ss_ub, &
       b_bv, lbw_bv, ubw_bv, lbwmax_bv, ubwmax_bv, numrots_bv, ks_bv, cs_bv, ss_bv, error)

    complex(kind=dp), dimension(n,maxsweeps_sw), intent(in) :: cs_sw, ss_sw
    integer(kind=int32), intent(in) :: numsweeps_sw, maxsweeps_sw, n

    complex(kind=dp), dimension(lbwmax_ub+ubwmax_ub+1,n), intent(inout) :: b_ub
    integer(kind=int32), dimension(n), intent(inout) :: numrots_ub
    integer(kind=int32), dimension(ubwmax_ub,n), intent(inout) :: js_ub
    complex(kind=dp), dimension(ubwmax_ub,n), intent(inout) :: cs_ub, ss_ub

    complex(kind=dp), dimension(n, lbwmax_bv+ubwmax_bv+1), intent(out) :: b_bv
    integer(kind=int32), dimension(n), intent(out) :: numrots_bv
    integer(kind=int32), dimension(n,ubwmax_bv), intent(out) :: ks_bv
    complex(kind=dp), dimension(n,ubwmax_bv), intent(out) :: cs_bv, ss_bv

    integer(kind=int32), intent(in) :: lbwmax_ub, ubwmax_ub
    integer(kind=int32), intent(inout) :: lbw_ub, ubw_ub
    integer(kind=int32), intent(in) :: lbwmax_bv, ubwmax_bv
    integer(kind=int32), intent(out) :: lbw_bv, ubw_bv

    type(error_info), intent(out) :: error
    integer(kind=int32) :: j, k, l, ubw, lbw, k0, k1, d
    type(c_rotation) :: rot
    logical :: full_ubw, full_lbw

    call clear_error(error)

    if (n==1) then
       b_bv(1,1)=b_ub(1,1); return
    end if
    ! Initial expansion for first sweep: one extra subdiagonal and
    ! one extra superdiagonal to fill-in.
    ! Expand as need for later sweeps.
    lbw=min(lbw_ub+numsweeps_sw,n-1)
    ubw=min(ubw_ub+numsweeps_sw,n-1)
    b_bv(:,1:lbw+ubw+1)=(0.0_dp,0.0_dp)
    numrots_bv=0
    ss_bv(:,1:ubw)=(0.0_dp,0.0_dp)
    cs_bv(:,1:ubw)=(0.0_dp,0.0_dp)
    ks_bv(:,1:ubw)=0

    ubw=ubw_ub
    lbw=lbw_ub

    full_lbw=(lbw==n-1)
    full_ubw=(ubw==n-1)

    do l=numsweeps_sw,1,-1
       if (l < numsweeps_sw) then
          if (.not. full_ubw) then
             ubw=ubw-1
          end if
          call f_c_convert_bv_to_ub(b_bv, n, lbw, ubw, lbwmax_bv, ubwmax_bv, &
               numrots_bv, ks_bv, cs_bv, ss_bv, b_ub, lbw_ub, ubw_ub, lbwmax_ub, ubwmax_ub, & 
               numrots_ub, js_ub, cs_ub, ss_ub, error)
       end if
       numrots_bv=0
       ! Unless all band structure is filled, sweep will fill in one sub and one super diagonal.
       if (.not. full_ubw) then
          ! ubw==n-3 results in ubw=n-1, but this is only due to the temporary diagonal.
          ! So this case does not set full_ubw to .true.  It will be set on the next sweep.
          if (ubw <= n-3) then  
             ubw=ubw+2
             call down_shift(b_ub)
             call down_shift(b_ub)
          else if (ubw==n-2) then
             ubw=n-1
             call down_shift(b_ub)
             full_ubw=.true.
          end if
       end if
       if (.not. full_lbw) then
          lbw=lbw+1;
          full_lbw=(lbw==n-1)
          b_ub(lbw+ubw+1,:)=(0.0_dp, 0.0_dp)
       end if
       do k=1,n-1
          do j=1,numrots_ub(n-k)
             rot%cosine=cs_ub(j,n-k); rot%sine=ss_ub(j,n-k)
             call rotation_times_tbc(rot,b_ub,n,lbw,ubw,n-k,0,js_ub(j,n-k))
          end do
          ! Apply q_{n-k}
          rot%cosine=cs_sw(n-k,l); rot%sine=ss_sw(n-k,l)
          call rotation_times_tbc(rot,b_ub,n,lbw,ubw,0,0,n-k)
          ! zeros have appeared in the second superdiagonal in
          ! rows n-k-(ubw-2), ... , n-k or
          ! columns n-k+2, ..., n-k+ubw
          k0=max(n-k+2,ubw+1)
          k1=min(n-k+ubw,n)
          numrots_bv(k)=max(k1-k0+1,0)
          do j=k0,k1
             rot=rgivens(get_el_bc(b_ub,ubw,j-ubw,j-1), get_el_bc(b_ub,ubw,j-ubw,j))
             ks_bv(k,k1-j+1)=j-1
             cs_bv(k,k1-j+1)=rot%cosine; ss_bv(k,k1-j+1)=rot%sine
             call tbc_times_rotation(b_ub,n,lbw,ubw,0,k,rot,j-1)
          end do
       end do
       call bc_to_br(b_ub,b_bv,lbw,ubw)
    end do
    lbw_bv=lbw
    if (.not. full_ubw) then
       ubw_bv=ubw-1
    else
       ubw_bv=ubw
    end if
  end subroutine f_c_sweeps_times_ub


  ! bv times sweeps
  ! Errors:
  ! 1: Sizes of bv, ub, and sw do not match.
  ! 2: Insufficient lbw in ub.
  ! 3: Insufficient ubw in ub.
  ! 4: Insufficient lbw in bv
  ! 5: Insufficient ubw in bv.
  ! 6: n < 1.
  ! 7: sw is transposed.

  subroutine d_bv_times_sweeps(bv,sw,ub,error)
    type(d_sweeps) :: sw
    type(d_ub) :: ub
    type(d_bv) :: bv
    type(error_info) :: error
    call clear_error(error)
    if (get_n(bv) /= get_n(ub) .or. get_n(bv) /= get_n(sw)) then
       call set_error(error,1,id_d_bv_times_sweeps); return
    end if
    if (get_lbwmax(ub) < sw%numsweeps + bv%lbw .and. get_lbwmax(ub) < get_n(bv)-1) then
       call set_error(error,2,id_d_bv_times_sweeps); return
    end if
    if (get_ubwmax(ub) < sw%numsweeps + bv%ubw+1 .and. get_ubwmax(ub) < get_n(bv)-1) then
       call set_error(error,3,id_d_bv_times_sweeps); return
    end if
    if (get_lbwmax(bv) < sw%numsweeps + bv%lbw .and. get_lbwmax(bv) < get_n(bv)-1) then
       call set_error(error,4,id_d_bv_times_sweeps); return
    end if
    if (get_ubwmax(bv) < sw%numsweeps + bv%ubw+2 .and. get_ubwmax(bv) < get_n(bv)+1) then
       call set_error(error,5,id_d_bv_times_sweeps); return
    end if
    if (get_n(bv) < 1) then
       call set_error(error,6,id_d_bv_times_sweeps); return
    end if
    if (sw%transposed) then
       call set_error(error,7,id_d_bv_times_sweeps); return
    end if
    call f_d_bv_times_sweeps(bv%b, get_n(bv), bv%lbw, bv%ubw, get_lbwmax(bv), get_ubwmax(bv), &
         bv%numrotsv, bv%ksv, bv%csv, bv%ssv, &
         sw%cs, sw%ss, sw%numsweeps, get_maxsweeps(sw), &
         ub%b, ub%lbw, ub%ubw, get_lbwmax(ub), get_ubwmax(ub), &
         ub%numrotsu, ub%jsu, ub%csu, ub%ssu, error)
  end subroutine d_bv_times_sweeps


  subroutine f_d_bv_times_sweeps(b_bv, n, lbw_bv, ubw_bv, lbwmax_bv, ubwmax_bv, numrots_bv, ks_bv, & 
       cs_bv, ss_bv, &
       cs_sw, ss_sw, numsweeps_sw, maxsweeps_sw, &
       b_ub, lbw_ub, ubw_ub, lbwmax_ub, ubwmax_ub, numrots_ub, js_ub, cs_ub, ss_ub, error)

    real(kind=dp), dimension(n,maxsweeps_sw), intent(in) :: cs_sw, ss_sw
    integer(kind=int32), intent(in) :: numsweeps_sw, maxsweeps_sw, n

    real(kind=dp), dimension(n, lbwmax_bv+ubwmax_bv+1), intent(inout) :: b_bv
    integer(kind=int32), dimension(n), intent(inout) :: numrots_bv
    integer(kind=int32), dimension(n,ubwmax_bv), intent(inout) :: ks_bv
    real(kind=dp), dimension(n,ubwmax_bv), intent(inout) :: cs_bv, ss_bv

    real(kind=dp), dimension(lbwmax_ub+ubwmax_ub+1,n), intent(out) :: b_ub
    integer(kind=int32), dimension(n), intent(out) :: numrots_ub
    integer(kind=int32), dimension(ubwmax_ub,n), intent(out) :: js_ub
    real(kind=dp), dimension(ubwmax_ub,n), intent(out) :: cs_ub, ss_ub

    integer(kind=int32), intent(in) :: lbwmax_ub, ubwmax_ub
    integer(kind=int32), intent(inout) :: lbw_ub, ubw_ub
    integer(kind=int32), intent(in) :: lbwmax_bv, ubwmax_bv
    integer(kind=int32), intent(out) :: lbw_bv, ubw_bv

    type(error_info), intent(out) :: error
    integer(kind=int32) :: j, k, l, ubw, lbw, k0, k1, d
    type(d_rotation) :: rot
    logical :: full_ubw, full_lbw

    call clear_error(error)

    if (n==1) then
       b_ub(1,1)=b_bv(1,1); return
    end if
    ! Initial expansion for first sweep: one extra subdiagonal and
    ! one extra superdiagonal to fill-in.
    ! Expand as need for later sweeps.
    lbw=min(lbw_bv+numsweeps_sw,n-1)
    ubw=min(ubw_bv+numsweeps_sw,n-1)
    b_ub(1:lbw+ubw+1,:)=0.0_dp
    numrots_ub=0
    ss_ub(1:ubw,:)=0.0_dp; cs_ub(1:ubw,:)=0.0_dp
    js_ub(1:ubw,:)=0

    ubw=ubw_bv
    lbw=lbw_bv

    full_lbw=(lbw==n-1)
    full_ubw=(ubw==n-1)

    do l=1, numsweeps_sw
       if (l > 1) then
          if (.not. full_ubw) then
             ubw=ubw-1
             call up_shift(b_ub)
          end if
          lbw_ub=lbw; ubw_ub=ubw
          call f_d_convert_ub_to_bv(b_ub, n, lbw_ub, ubw_ub, lbwmax_ub, ubwmax_ub, & 
               numrots_ub, js_ub, cs_ub, ss_ub, &
               b_bv, lbw, ubw, lbwmax_bv, ubwmax_bv, numrots_bv, ks_bv, cs_bv, ss_bv, error)
       end if
       numrots_ub=0
       ! Unless all band structure is filled, sweep will fill in one sub and one super diagonal.
       if (.not. full_ubw) then
          ! ubw==n-3 results in ubw=n-1, but this is only due to the temporary diagonal.
          ! So this case does not set full_ubw to .true.  It will be set on the next sweep.
          if (ubw <= n-3) then  
             ubw=ubw+2
             b_bv(:,lbw+ubw+1)=0.0_dp
             b_bv(:,lbw+ubw)=0.0_dp
          else if (ubw==n-2) then
             ubw=n-1
             b_bv(:,lbw+ubw+1)=0.0_dp
             full_ubw=.true.
          end if
       end if
       if (.not. full_lbw) then
          lbw=lbw+1;
          full_lbw=(lbw==n-1)
          call right_shift(b_bv)
       end if
       do k=1,n-1
          do j=1,numrots_bv(n-k)
             rot%cosine=cs_bv(n-k,j); rot%sine=ss_bv(n-k,j)
             call tbr_times_rotation(b_bv,n,lbw,ubw,0,n-k,trp_rot(rot),ks_bv(n-k,j))
          end do
          ! Apply q_{k}
          rot%cosine=cs_sw(k,l); rot%sine=ss_sw(k,l)
          call tbr_times_rotation(b_bv,n,lbw,ubw,0,0,rot,k)
          ! zeros have appeared in the second superdiagonal
          ! in columns k+1, ..., k+ubw-1
          ! or rows k+1-ubw, ..., k-1
          k0=max(k+1,ubw+1)
          k1=min(k+ubw-1,n)
          numrots_ub(k)=max(k1-k0+1,0)
          do j=k1,k0,-1
             rot=lgivens2(get_el_br(b_bv,lbw,j-ubw,j), get_el_br(b_bv,lbw,j-ubw+1,j))
             js_ub(j-k0+1,k)=j-ubw
             cs_ub(j-k0+1,k)=rot%cosine; ss_ub(j-k0+1,k)=rot%sine
             call rotation_times_tbr(trp_rot(rot), b_bv,n,lbw,ubw,k,0,j-ubw)
          end do
       end do
       call br_to_bc(b_bv,b_ub,lbw,ubw)
    end do
    lbw_ub=lbw
    if (.not. full_ubw) then
       ubw_ub=ubw-1
       call up_shift(b_ub)
    else
       ubw_ub=ubw
    end if
  end subroutine f_d_bv_times_sweeps

  subroutine c_bv_times_sweeps(bv,sw,ub,error)
    type(c_sweeps) :: sw
    type(c_ub) :: ub
    type(c_bv) :: bv
    type(error_info) :: error
    call clear_error(error)
    if (get_n(bv) /= get_n(ub) .or. get_n(bv) /= get_n(sw)) then
       call set_error(error,1,id_c_bv_times_sweeps); return
    end if
    if (get_lbwmax(ub) < sw%numsweeps + bv%lbw .and. get_lbwmax(ub) < get_n(bv)-1) then
       call set_error(error,2,id_c_bv_times_sweeps); return
    end if
    if (get_ubwmax(ub) < sw%numsweeps + bv%ubw+1 .and. get_ubwmax(ub) < get_n(bv)-1) then
       call set_error(error,3,id_c_bv_times_sweeps); return
    end if
    if (get_lbwmax(bv) < sw%numsweeps + bv%lbw .and. get_lbwmax(bv) < get_n(bv)-1) then
       call set_error(error,4,id_c_bv_times_sweeps); return
    end if
    if (get_ubwmax(bv) < sw%numsweeps + bv%ubw+2 .and. get_ubwmax(bv) < get_n(bv)+1) then
       call set_error(error,5,id_c_bv_times_sweeps); return
    end if
    if (get_n(bv) < 1) then
       call set_error(error,6,id_c_bv_times_sweeps); return
    end if
    if (sw%transposed) then
       call set_error(error,7,id_c_bv_times_sweeps); return
    end if
    call f_c_bv_times_sweeps(bv%b, get_n(bv), bv%lbw, bv%ubw, get_lbwmax(bv), get_ubwmax(bv), &
         bv%numrotsv, bv%ksv, bv%csv, bv%ssv, &
         sw%cs, sw%ss, sw%numsweeps, get_maxsweeps(sw), &
         ub%b, ub%lbw, ub%ubw, get_lbwmax(ub), get_ubwmax(ub), &
         ub%numrotsu, ub%jsu, ub%csu, ub%ssu, error)
  end subroutine c_bv_times_sweeps


  subroutine f_c_bv_times_sweeps(b_bv, n, lbw_bv, ubw_bv, lbwmax_bv, ubwmax_bv, numrots_bv, ks_bv, & 
       cs_bv, ss_bv, &
       cs_sw, ss_sw, numsweeps_sw, maxsweeps_sw, &
       b_ub, lbw_ub, ubw_ub, lbwmax_ub, ubwmax_ub, numrots_ub, js_ub, cs_ub, ss_ub, error)

    complex(kind=dp), dimension(n,maxsweeps_sw), intent(in) :: cs_sw, ss_sw
    integer(kind=int32), intent(in) :: numsweeps_sw, maxsweeps_sw, n

    complex(kind=dp), dimension(n, lbwmax_bv+ubwmax_bv+1), intent(inout) :: b_bv
    integer(kind=int32), dimension(n), intent(inout) :: numrots_bv
    integer(kind=int32), dimension(n,ubwmax_bv), intent(inout) :: ks_bv
    complex(kind=dp), dimension(n,ubwmax_bv), intent(inout) :: cs_bv, ss_bv

    complex(kind=dp), dimension(lbwmax_ub+ubwmax_ub+1,n), intent(out) :: b_ub
    integer(kind=int32), dimension(n), intent(out) :: numrots_ub
    integer(kind=int32), dimension(ubwmax_ub,n), intent(out) :: js_ub
    complex(kind=dp), dimension(ubwmax_ub,n), intent(out) :: cs_ub, ss_ub

    integer(kind=int32), intent(in) :: lbwmax_ub, ubwmax_ub
    integer(kind=int32), intent(inout) :: lbw_ub, ubw_ub
    integer(kind=int32), intent(in) :: lbwmax_bv, ubwmax_bv
    integer(kind=int32), intent(out) :: lbw_bv, ubw_bv

    type(error_info), intent(out) :: error
    integer(kind=int32) :: j, k, l, ubw, lbw, k0, k1, d
    type(c_rotation) :: rot
    logical :: full_ubw, full_lbw

    call clear_error(error)

    if (n==1) then
       b_ub(1,1)=b_bv(1,1); return
    end if
    ! Initial expansion for first sweep: one extra subdiagonal and
    ! one extra superdiagonal to fill-in.
    ! Expand as need for later sweeps.
    lbw=min(lbw_bv+numsweeps_sw,n-1)
    ubw=min(ubw_bv+numsweeps_sw,n-1)
    b_ub(1:lbw+ubw+1,:)=(0.0_dp,0.0_dp)
    numrots_ub=0
    ss_ub(1:ubw,:)=(0.0_dp,0.0_dp)
    cs_ub(1:ubw,:)=(0.0_dp,0.0_dp)
    js_ub(1:ubw,:)=0

    ubw=ubw_bv
    lbw=lbw_bv

    full_lbw=(lbw==n-1)
    full_ubw=(ubw==n-1)

    do l=1, numsweeps_sw
       if (l > 1) then
          if (.not. full_ubw) then
             ubw=ubw-1
             call up_shift(b_ub)
          end if
          lbw_ub=lbw; ubw_ub=ubw
          call f_c_convert_ub_to_bv(b_ub, n, lbw_ub, ubw_ub, lbwmax_ub, ubwmax_ub, & 
               numrots_ub, js_ub, cs_ub, ss_ub, &
               b_bv, lbw, ubw, lbwmax_bv, ubwmax_bv, numrots_bv, ks_bv, cs_bv, ss_bv, error)
       end if
       numrots_ub=0
       ! Unless all band structure is filled, sweep will fill in one sub and one super diagonal.
       if (.not. full_ubw) then
          ! ubw==n-3 results in ubw=n-1, but this is only due to the temporary diagonal.
          ! So this case does not set full_ubw to .true.  It will be set on the next sweep.
          if (ubw <= n-3) then  
             ubw=ubw+2
             b_bv(:,lbw+ubw+1)=(0.0_dp, 0.0_dp)
             b_bv(:,lbw+ubw)=(0.0_dp, 0.0_dp)
          else if (ubw==n-2) then
             ubw=n-1
             b_bv(:,lbw+ubw+1)=(0.0_dp, 0.0_dp)
             full_ubw=.true.
          end if
       end if
       if (.not. full_lbw) then
          lbw=lbw+1;
          full_lbw=(lbw==n-1)
          call right_shift(b_bv)
       end if
       do k=1,n-1
          do j=1,numrots_bv(n-k)
             rot%cosine=cs_bv(n-k,j); rot%sine=ss_bv(n-k,j)
             call tbr_times_rotation(b_bv,n,lbw,ubw,0,n-k,trp_rot(rot),ks_bv(n-k,j))
          end do
          ! Apply q_{k}
          rot%cosine=cs_sw(k,l); rot%sine=ss_sw(k,l)
          call tbr_times_rotation(b_bv,n,lbw,ubw,0,0,rot,k)
          ! zeros have appeared in the second superdiagonal
          ! in columns k+1, ..., k+ubw-1
          ! or rows k+1-ubw, ..., k-1
          k0=max(k+1,ubw+1)
          k1=min(k+ubw-1,n)
          numrots_ub(k)=max(k1-k0+1,0)
          do j=k1,k0,-1
             rot=lgivens2(get_el_br(b_bv,lbw,j-ubw,j), get_el_br(b_bv,lbw,j-ubw+1,j))
             js_ub(j-k0+1,k)=j-ubw
             cs_ub(j-k0+1,k)=rot%cosine; ss_ub(j-k0+1,k)=rot%sine
             call rotation_times_tbr(trp_rot(rot), b_bv,n,lbw,ubw,k,0,j-ubw)
          end do
       end do
       call br_to_bc(b_bv,b_ub,lbw,ubw)
    end do
    lbw_ub=lbw
    if (.not. full_ubw) then
       ubw_ub=ubw-1
       call up_shift(b_ub)
    else
       ubw_ub=ubw
    end if
  end subroutine f_c_bv_times_sweeps


end module sweeps
