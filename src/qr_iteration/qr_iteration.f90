module qr_iteration
  use compressions
  use conversions

  implicit none

  type c_shift
     logical :: flag
     complex(kind=dp) :: shift
  end type c_shift

  real(kind=dp), private, parameter :: tol=1.0e-15

  type(routine_info), parameter :: info_ss_qr=routine_info(id_ss_qr, &
       'ss_qr', &
       [ character(len=error_message_length) :: 'size(q,2) /= bv%n', 'bv%lbw /= 1'  ])

  type(routine_info), parameter :: info_ss_qr_iteration=routine_info(id_ss_qr_iteration, &
       'ss_qr_iteration', &
       [ character(len=error_message_length) :: 'n<1', 'Insufficient storage in bv.', &
       'ub and bv are not the same size.', 'bv is not Hessenberg.', 'Insufficient storage in ub.'])

contains

  ! Errors:
  ! 1: size(q,2) /= bv%n
  ! 2: bv%lbw /= 1
  !         
  subroutine ss_qr(bv,q,told,tol,error)
    type(c_bv), intent(inout) :: bv
    complex(kind=dp), dimension(:,:), intent(out) :: q
    type(error_info), intent(out) :: error
    real(kind=dp), intent(in) :: told, tol

    type(c_shift), dimension(:), allocatable :: shifts
    integer(kind=int32) :: k, numshifts, n, lbw, ubw
    type(c_bv) :: bv_tmp
    type(c_ub) :: ub
    complex(kind=dp) :: subd, sigma

    n=get_n(bv); lbw=bv%lbw; ubw=bv%ubw
    allocate(shifts(n))

    call clear_error(error)
    if (size(q,2) /= n) then
       call set_error(error, 1, id_ss_qr); return
    end if
    if (lbw /= 1) then
       call set_error(error, 2, id_ss_qr); return
    end if
    bv_tmp=c_new_bv(n,get_lbwmax(bv), get_ubwmax(bv))

    ub=c_new_ub(n,get_lbwmax(bv), get_ubwmax(bv))

    qr_loop: do
       ! zero any small subdiagonals
       do k=1,n-1
          subd=get_el_br(bv%b,lbw,k+1,k)
          if ( abs(subd) <= tol ) then
             call set_el_br(bv%b,lbw,k+1,k,(0.0_dp,0.0_dp))
          end if
       end do
       ! set up shifts.
       numshifts=0
       call copy_bv(bv,bv_tmp)
       call c_bv_reveal_superdiag(bv_tmp)

       ! compute shifts
       shifts%flag=.false.
       shifts%shift=(0.0_dp,0.0_dp)
       k=n-1
       numshifts=0
       shift_loop: do while (k >= 1)
          subd=get_el_br(bv%b,lbw,k+1,k)
          if (subd==(0.0_dp, 0.0_dp)) then
             k=k-1;  cycle shift_loop
          else
             ! Found a nonzero subdiagonal
             sigma=c_wilkinson_shift(get_el_br(bv%b,lbw,k,k),get_el_br(bv_tmp%b,lbw,k,k+1), &
                  subd,get_el_br(bv%b,lbw,k+1,k+1))
             numshifts=numshifts+1
             ! find the next zero subdiagonal
             find_unreduced: do
                if (k==1) then
                   ! shift starts at the top.
                   shifts(1)%flag=.true.
                   shifts(1)%shift=sigma
                   exit shift_loop
                else
                   subd=get_el_br(bv%b,lbw,k,k-1)
                   if (subd==(0.0_dp,0.0_dp)) then
                      shifts(k)%flag=.true.
                      shifts(k)%shift=sigma
                      k=k-1
                      exit find_unreduced
                   else
                      k=k-1
                      cycle find_unreduced
                   end if
                end if
             end do find_unreduced
          end if
       end do shift_loop

       if (numshifts==0) then
          exit qr_loop
       end if
       call ss_qr_iteration(bv,ub,shifts, q, error)
       call compress_ub_to_bv(ub,bv,told, tol,0,error)
    end do qr_loop
    deallocate(shifts)
    call deallocate_ub(ub)
  end subroutine ss_qr

  subroutine f_ss_qr(b_bv,n, lbw_bv, ubw_bv, lbwmax_bv, ubwmax_bv, numrots_bv, ks_bv, &
       cs_bv, ss_bv, q, p, told, tol, error)
    integer(kind=int32), intent(in) :: n, lbw_bv, lbwmax_bv, ubwmax_bv, p
    integer(kind=int32), intent(inout) :: ubw_bv
    complex(kind=dp), dimension(lbwmax_bv+ubwmax_bv+1,n), intent(inout) :: b_bv
    integer(kind=int32), dimension(n), intent(out) :: numrots_bv
    integer(kind=int32), dimension(ubwmax_bv,n), intent(out) :: ks_bv
    complex(kind=dp), dimension(ubwmax_bv,n), intent(out) :: cs_bv, ss_bv
    complex(kind=dp), dimension(p,n), intent(out) :: q
    type(error_info), intent(out) :: error
    real(kind=dp), intent(in) :: told, tol

    type (c_bv) :: bv

    call clear_error(error)
    bv=c_new_bv(n,lbwmax_bv,ubwmax_bv)

    bv%b=b_bv
    bv%lbw=lbw_bv; bv%ubw=ubw_bv
    bv%numrotsv=numrots_bv
    bv%ksv=ks_bv
    bv%csv=cs_bv; bv%ssv=ss_bv
    call ss_qr(bv,q,told, tol, error)
    b_bv=bv%b
    ubw_bv=bv%ubw
    numrots_bv=bv%numrotsv
    ks_bv=bv%ksv
    cs_bv=bv%csv; ss_bv=bv%ssv

    call deallocate_bv(bv)

  end subroutine f_ss_qr

  ! Errors:
  ! 1: n<1
  ! 2: insufficient storage in bv.
  ! 3: ub%n /= bv%n
  ! 4: Not Hessenberg: bv%lbw /= 1
  ! 5: insufficient storage in ub

  subroutine ss_qr_iteration(bv,ub,shifts, q, error)
    type(c_bv), intent(inout) :: bv
    type(c_ub) :: ub
    type(c_shift), dimension(:), intent(in) :: shifts
    type(error_info), intent(out) :: error
    complex(kind=dp), dimension(:,:), intent(out) :: q

    integer :: n, j
    integer(kind=int32), dimension(size(shifts)) :: shifts_i
    complex(kind=dp), dimension(size(shifts)) :: shifts_c

    n=get_n(ub)
    call clear_error(error)
    if (n < 1) then
       call set_error(error, 1, id_ss_qr_iteration); return
    end if
    if (get_lbwmax(bv)+get_ubwmax(bv)<bv%ubw + 4) then
       call set_error(error, 2, id_ss_qr_iteration); return
    end if
    if (get_n(bv) /= n .or. size(shifts) /= n .or. size(q,2) /= n) then
       call set_error(error, 3, id_ss_qr_iteration); return
    end if
    if (bv%lbw /= 1) then
       call set_error(error, 4, id_ss_qr_iteration); return
    end if
    if (get_lbwmax(ub) < 1 .or. get_ubwmax(ub) < bv%ubw) then
       call set_error(error, 5, id_ss_qr_iteration); return
    end if

    shifts_i=0
    do j=1,n
       shifts_c(j)=shifts(j)%shift
       if (shifts(j)%flag .eqv. .true.) then
          shifts_i(j)=1
       end if
    end do
    call f_ss_qr_iteration(bv%b, get_n(bv), bv%ubw, get_lbwmax(bv), get_ubwmax(bv), &
         bv%numrotsv, bv%ksv, bv%csv, bv%ssv, & 
         ub%b, get_lbwmax(ub), get_ubwmax(ub), ub%numrotsu, ub%jsu, &
         ub%csu, ub%ssu, shifts_i, shifts_c, q, size(q,2), error)
    ub%ubw=bv%ubw+1
    ub%lbw=1
  end subroutine ss_qr_iteration

  subroutine f_ss_qr_iteration(b_bv,n, ubw_bv, lbwmax_bv, ubwmax_bv, numrots_bv, ks_bv, &
       cs_bv, ss_bv, b_ub, lbwmax_ub, ubwmax_ub, numrots_ub, js_ub, cs_ub, &
       ss_ub, shifts_i, shifts_c, q, p, error)
    integer(kind=int32), intent(in) :: n, ubw_bv, lbwmax_bv, ubwmax_bv, lbwmax_ub, ubwmax_ub, p
    complex(kind=dp), dimension(n,lbwmax_bv+ubwmax_bv+1), intent(inout) :: b_bv
    integer(kind=int32), dimension(n), intent(in) :: numrots_bv
    integer(kind=int32), dimension(n,ubwmax_bv), intent(in) :: ks_bv
    complex(kind=dp), dimension(n,ubwmax_bv), intent(in) :: cs_bv, ss_bv

    complex(kind=dp), dimension(lbwmax_ub+ubwmax_ub+1,n), intent(out) :: b_ub
    integer(kind=int32), dimension(n), intent(out) :: numrots_ub
    integer(kind=int32), dimension(ubwmax_ub,n), intent(out) :: js_ub
    complex(kind=dp), dimension(ubwmax_ub,n), intent(out) :: cs_ub, ss_ub
    complex(kind=dp), dimension(p,n), intent(out) :: q

    integer(kind=int32), dimension(n) :: shifts_i
    complex(kind=dp), dimension(n) :: shifts_c

    type(error_info), intent(out) :: error

    integer(kind=int32) :: j, k, ubw, lbw, d, ubw_ub, k0,k1
    type(c_rotation) :: rot

    call clear_error(error)
    numrots_ub=0
    ss_ub=(0.0_dp, 0.0_dp); cs_ub=(0.0_dp, 0.0_dp)
    js_ub=0
    b_ub=(0.0_dp, 0.0_dp)
    ubw=ubw_bv+2 ! ubw increases by 1, with extra diagonal for temporary storage.
    lbw=2 ! 1 extra subdiagonal for the bulge.

    if (n == 1) then
       b_ub(1,1)=b_bv(1,1)
       return
    end if

    ! create room for the extra subdiagonal
    call right_shift(b_bv)
    do k=1,n-1
       ! apply v_{n-k}
       do j=1,numrots_bv(n-k)
          rot%cosine=cs_bv(n-k,j); rot%sine=ss_bv(n-k,j)
          call tbr_times_rotation(b_bv,n,lbw,ubw,0,n-k,trp_rot(rot),ks_bv(n-k,j))
       end do
       if (shifts_i(k)==1) then
          ! introduce a new shift.
          rot=lgivens(get_el_br(b_bv,lbw,k,k)-shifts_c(k),get_el_br(b_bv,lbw,k+1,k))
          call rotation_times_tbr(trp_rot(rot),b_bv,n,lbw,ubw,0,0,k)
          call tbr_times_rotation(b_bv,n,lbw,ubw,0,0,rot,k)
          call general_times_rotation(q,rot,k,k+1)
       else if (get_el_br(b_bv,lbw,k+1,k-1) /= (0.0_dp, 0.0_dp)) then
          ! chase the bulge.
          rot=lgivens(get_el_br(b_bv,lbw,k,k-1),get_el_br(b_bv,lbw,k+1,k-1))
          call rotation_times_tbr(trp_rot(rot),b_bv,n,lbw,ubw,0,0,k)
          call tbr_times_rotation(b_bv,n,lbw,ubw,0,0,rot,k)
          call general_times_rotation(q,rot,k,k+1)
       end if
       ! columns that have a nonzero in superdiagonal ubw
       ! Three cases:
       ! 1. For k <= ubw-1, columns ubw+1, ..., k+ubw-1 might
       !    have a nonzero in superdiagonal ubw.
       ! 2. For ubw <= k <= n-ubw+1, columns k+1, .... k+ubw-1
       !    might have a nonzero in superdiagonal ubw.
       ! 3. For n-ubw+2 <= k <= n-1, columns k+1, ..., n
       !    might have a nonzero in superdiagonal ubw.
       k0=max(k+1,ubw+1)
       k1=min(k+ubw-1,n)
       numrots_ub(k)=k1-k0+1
       do j=k1,k0,-1
          rot=lgivens2(get_el_br(b_bv,lbw,j-ubw,j), get_el_br(b_bv,lbw,j-ubw+1,j))
          js_ub(j-k0+1,k)=j-ubw
          cs_ub(j-k0+1,k)=rot%cosine; ss_ub(j-k0+1,k)=rot%sine
          call rotation_times_tbr(trp_rot(rot),b_bv,n,lbw,ubw,k,0,j-ubw)
       end do
    end do
    ubw_ub=ubw_bv+1
    call left_shift(b_bv)
    ! put diagonals in b
    do d=1,ubw_ub+1
       do k=ubw_ub-d+2,n
          b_ub(d,k) = get_el_br(b_bv,1,d+k-ubw_ub-1,k)
       end do
    end do
    do d=ubw_ub+2, ubw_ub+1+1
       do k=1,n-d+ubw_ub+1
          b_ub(d,k) = get_el_br(b_bv,1,d+k-ubw_ub-1,k)
       end do
    end do
  end subroutine f_ss_qr_iteration

  ! These function are destructive.  bv is modified so that it
  ! is no longer really a nested bv decomposition.  The
  ! superdiagonal of A is revealed in bv%b.
  subroutine c_bv_reveal_superdiag(bv)
    type(c_bv) :: bv
    call f_c_bv_reveal_superdiag(bv%b,get_n(bv),bv%lbw,bv%ubw,get_lbwmax(bv),get_ubwmax(bv), &
         bv%numrotsv, bv%ksv, bv%csv, bv%ssv)
  end subroutine c_bv_reveal_superdiag

  ! Apply v_{n-j} to row j to reveal the first superdiagonal element in each row.
  ! requires on extra superdiagonal.
  subroutine f_c_bv_reveal_superdiag(b_bv, n, lbw_bv, ubw_bv, lbwmax_bv, ubwmax_bv, numrots_bv, ks_bv, &
       cs_bv, ss_bv)
    integer(kind=int32), intent(in) :: n, lbw_bv, ubw_bv, lbwmax_bv, ubwmax_bv
    complex(kind=dp), dimension(n,lbwmax_bv+ubwmax_bv+1), intent(inout) :: b_bv
    integer(kind=int32), dimension(n), intent(in) :: numrots_bv
    integer(kind=int32), dimension(n,ubwmax_bv), intent(in) :: ks_bv
    complex(kind=dp), dimension(n,ubwmax_bv), intent(in) :: cs_bv, ss_bv

    integer(kind=int32) :: j, k
    type(c_rotation) :: rot

    do j=1,n-1
       do k=1,numrots_bv(n-j)
          rot%cosine=cs_bv(n-j,k); rot%sine=ss_bv(n-j,k)
          call tbr_times_rotation(b_bv,n,lbw_bv,ubw_bv,j-1,n-j,trp_rot(rot),ks_bv(n-j,k))
       end do
    end do
  end subroutine f_c_bv_reveal_superdiag

  ! Compute the eigenvalue of
  ! [a,b;
  !  c,d]
  ! that is closest to d.
  complex(kind=dp) function c_wilkinson_shift(a,b,c,d) result(lambda)
    complex(kind=dp), intent(in) :: a,b,c,d
    complex(kind=dp) :: tr, w, x, y
    tr=a+d
    w=sqrt(tr**2-4.0*(a*d-b*c))
    x=(tr+w)/2.0
    y=(tr-w)/2.0
    if (abs(d-y) < abs(d-x)) then
       lambda=y
    else
       lambda=x
    end if

  end function c_wilkinson_shift

end module qr_iteration
