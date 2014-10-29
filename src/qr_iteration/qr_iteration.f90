module mod_qr_iteration
  use mod_prec
  use mod_error_id
  use mod_rotation
  use mod_nested_types
  use mod_band_types
  use mod_recompress_ub_to_bv
  use mod_convert_ub_to_bv
  use mod_convert_bv_to_ub
  use mod_qr_factorization
  use mod_update
  use mod_solve
  use mod_sweeps1

  implicit none

  private
  public :: c_shift, ss_r1_qr, f_ss_r1_qr, f_r1_reorth, ss_qr_iteration, f_ss_qr_iteration, &
       c_wilkinson_shift

  public :: info_ss_r1_qr, info_ss_qr_iteration

  type c_shift
     logical :: flag
     complex(kind=dp) :: shift
  end type c_shift

  real(kind=dp), private, parameter :: tol=1.0e-16

  type(routine_info), parameter :: info_ss_r1_qr=routine_info(id_ss_r1_qr, &
       'ss_r1_qr', &
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
  subroutine ss_r1_qr(bv,u,v,q,error)
    type(c_bv), intent(inout) :: bv
    complex(kind=dp), dimension(:,:), intent(inout) :: q
    complex(kind=dp), dimension(:), intent(inout) :: u, v
    type(error_info), intent(out) :: error

    complex(kind=dp), dimension(:), allocatable :: shifts_c
    complex(kind=dp), dimension(:,:), allocatable :: cs_sw, ss_sw, b_bv_tmp
    integer(kind=int32), dimension(:), allocatable :: shifts_i
    integer(kind=int32) :: n
    type(c_ub) :: ub
    
    call clear_error(error)
    n=get_n(bv)
    allocate(shifts_c(n), shifts_i(n), cs_sw(n,4), ss_sw(n,4), &
         b_bv_tmp(n, get_lbwmax(bv) + get_ubwmax(bv)+1))

    if (size(q,2) /= n) then
       call set_error(error, 1, id_ss_r1_qr); return
    end if
    if (bv%lbw /= 1) then
       call set_error(error, 2, id_ss_r1_qr); return
    end if

    ub=c_new_ub(n,get_lbwmax(bv), get_ubwmax(bv))
    
    call f_ss_r1_qr(bv%br, n, bv%lbw, bv%ubw, get_lbwmax(bv), get_ubwmax(bv), bv%numrotsv, bv%ksv, &
         bv%csv, bv%ssv, u,v,q,size(q,1), &
         ub%bc, ub%lbw, ub%ubw, get_lbwmax(ub), get_ubwmax(ub), ub%numrotsu, ub%jsu, ub%csu, ub%ssu, &
         cs_sw, ss_sw, shifts_c, shifts_i, b_bv_tmp, error)
    deallocate(shifts_c, shifts_i, cs_sw, ss_sw, b_bv_tmp)
    call deallocate_ub(ub)

  end subroutine ss_r1_qr

  subroutine f_ss_r1_qr(b_bv, n, lbw_bv, ubw_bv, lbwmax_bv, ubwmax_bv, numrotsv, ksv, &
       csv, ssv, u, v, q, p, &
       b_ub, lbw_ub, ubw_ub, lbwmax_ub, ubwmax_ub, numrotsu, jsu, csu, ssu, &
       cs_sw, ss_sw, shifts_c, shifts_i, b_bv_tmp, error)

    integer(kind=int32), intent(in) :: n, lbwmax_bv, ubwmax_bv, lbwmax_ub, ubwmax_ub, p
    integer(kind=int32), intent(inout) :: lbw_bv, ubw_bv, lbw_ub, ubw_ub

    complex(kind=dp), dimension(n, lbwmax_bv+ubwmax_bv+1), intent(inout) :: b_bv
    integer(kind=int32), dimension(n), intent(inout) :: numrotsv
    integer(kind=int32), dimension(n,ubwmax_bv), intent(inout) :: ksv
    complex(kind=dp), dimension(n,ubwmax_bv), intent(inout) :: csv, ssv

    complex(kind=dp), dimension(lbwmax_ub+ubwmax_ub+1,n), intent(out) :: b_ub
    integer(kind=int32), dimension(n), intent(out) :: numrotsu
    integer(kind=int32), dimension(ubwmax_ub,n), intent(out) :: jsu
    complex(kind=dp), dimension(ubwmax_ub,n), intent(out) :: csu, ssu

    complex(kind=dp), dimension(n,4) :: cs_sw, ss_sw
    complex(kind=dp), dimension(n) :: shifts_c
    integer(kind=int32), dimension(n) :: shifts_i

    complex(kind=dp), dimension(p,n), intent(inout) :: q
    complex(kind=dp), dimension(n), intent(inout) :: u, v
    complex(kind=dp), dimension(n,lbwmax_bv+ubwmax_bv+1) :: b_bv_tmp
    type(error_info), intent(out) :: error

    integer(kind=int32) :: k, numshifts
    complex(kind=dp) :: subd, sigma
    type(c_rotation) :: rot
    complex(kind=dp), dimension(n) :: u_tmp, v_tmp

    call clear_error(error)

   qr_loop: do

       do k=1,n-1
          subd=get_el_br(b_bv,lbw_bv,k+1,k)
          if ( abs(subd) <= tol ) then
             call set_el_br(b_bv,lbw_bv,k+1,k,(0.0_dp,0.0_dp))
          end if
       end do
       ! set up shifts.
       numshifts=0
       b_bv_tmp(:,1:lbw_bv+ubw_bv+1)=b_bv(:,1:lbw_bv+ubw_bv+1)
       call f_c_bv_reveal_superdiag(b_bv_tmp, n, lbw_bv, ubw_bv, lbwmax_bv, ubwmax_bv, numrotsv, ksv, &
            csv, ssv)
       shifts_i=0
       shifts_c=(0.0_dp,0.0_dp)
       k=n-1
       numshifts=0
       shift_loop: do while (k >= 1)
          subd=get_el_br(b_bv,lbw_bv,k+1,k)
          if (subd==(0.0_dp, 0.0_dp)) then
             k=k-1;  cycle shift_loop
          else
             ! Found a nonzero subdiagonal
             sigma=c_wilkinson_shift(get_el_br(b_bv,lbw_bv,k,k),get_el_br(b_bv_tmp,lbw_bv,k,k+1), &
                  subd,get_el_br(b_bv,lbw_bv,k+1,k+1))
             numshifts=numshifts+1
             ! find the next zero subdiagonal
             ! find the next zero subdiagonal
             find_unreduced: do
                if (k==1) then
                   ! shift starts at the top.
                   shifts_i(1)=1
                   shifts_c(1)=sigma
                   exit shift_loop
                else
                   subd=get_el_br(b_bv,lbw_bv,k,k-1)
                   if (subd==(0.0_dp,0.0_dp)) then
                      shifts_i(k)=1
                      shifts_c(k)=sigma
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

       call f_ss_qr_iteration(b_bv,n, lbw_bv, ubw_bv, lbwmax_bv, ubwmax_bv, numrotsv, ksv, &
            csv, ssv, b_ub, lbw_ub, ubw_ub, lbwmax_ub, ubwmax_ub, numrotsu, jsu, csu, &
            ssu, shifts_i, shifts_c, cs_sw(:,1), ss_sw(:,1), error)

       do k=1,n-1
          rot%cosine=cs_sw(k,1); rot%sine=ss_sw(k,1)
          call general_times_rotation(q,rot,k,k+1)
          call rotation_times_general(trp_rot(rot), u, k,k+1)
          call rotation_times_general(trp_rot(rot), v, k,k+1)
       end do

      if (ubw_ub > 3) then

          call f_r1_reorth(b_ub, n, lbw_ub, ubw_ub, lbwmax_ub, ubwmax_ub, &
               numrotsu, jsu, csu, ssu, u, v, u_tmp, v_tmp, &
               b_bv, lbw_bv, ubw_bv, lbwmax_bv, ubwmax_bv, numrotsv, ksv, csv, ssv, &
               cs_sw, ss_sw, shifts_i, error)

       else
          call f_c_convert_ub_to_bv(b_ub, n, lbw_ub, ubw_ub, lbwmax_ub, ubwmax_ub, numrotsu, &
               jsu, csu, ssu, &
               b_bv, lbw_bv, ubw_bv, lbwmax_bv, ubwmax_bv, numrotsv, ksv, csv, ssv, error)
       end if

    end do qr_loop

  end subroutine f_ss_r1_qr


  ! reorthogonalization

  subroutine f_r1_reorth(b_ub, n, lbw_ub, ubw_ub, lbwmax_ub, ubwmax_ub, &
       numrotsu, jsu, csu, ssu, u, v, u_tmp, v_tmp, &
       b_bv, lbw_bv, ubw_bv, lbwmax_bv, ubwmax_bv, numrotsv, ksv, csv, ssv, &
       cs_sw, ss_sw, zeros_i, error)
    integer(kind=int32), intent(in) :: n, lbwmax_bv, ubwmax_bv, lbwmax_ub, ubwmax_ub
    integer(kind=int32), intent(inout) :: lbw_bv, ubw_bv, lbw_ub, ubw_ub

    complex(kind=dp), dimension(n, lbwmax_bv+ubwmax_bv+1), intent(out) :: b_bv
    integer(kind=int32), dimension(n), intent(out) :: numrotsv
    integer(kind=int32), dimension(n,ubwmax_bv), intent(out) :: ksv
    complex(kind=dp), dimension(n,ubwmax_bv), intent(out) :: csv, ssv

    complex(kind=dp), dimension(lbwmax_ub+ubwmax_ub+1,n), intent(inout) :: b_ub
    integer(kind=int32), dimension(n), intent(inout) :: numrotsu
    integer(kind=int32), dimension(ubwmax_ub, n), intent(inout) :: jsu
    complex(kind=dp), dimension(ubwmax_ub, n), intent(inout) :: csu, ssu

    complex(kind=dp), dimension(n,4) :: cs_sw, ss_sw
    complex(kind=dp), dimension(n), intent(inout) :: u, v, u_tmp, v_tmp
    integer(kind=int32), dimension(n), intent(out) :: zeros_i

    type(error_info), intent(out) :: error

    integer(kind=int32) :: j, numsweeps1
    complex(kind=dp), dimension(n) :: d
    complex(kind=dp) :: x

    zeros_i=0
    do j=2,n
       if ( abs(b_bv(j,1)) <= tol ) then
          zeros_i(j)=1
       end if
    end do
    
    numsweeps1=0
    u_tmp=u; v_tmp=v

    call f_c_r1_update_ub_to_bv(b_ub, n, lbw_ub, ubw_ub, lbwmax_ub, ubwmax_ub, numrotsu, &
         jsu, csu, ssu, u_tmp, v_tmp, &
         b_bv, lbw_bv, ubw_bv, lbwmax_bv, ubwmax_bv, numrotsv, ksv, csv, ssv, &
         cs_sw, ss_sw, numsweeps1, 4, error)
    
    b_bv(1,lbw_bv+1)=b_bv(1,lbw_bv+1) - u_tmp(1)*conjg(v_tmp(1))
    b_bv(1,lbw_bv+2)=b_bv(1,lbw_bv+2) - u_tmp(1)*conjg(v_tmp(2))

    call f_c_qr_bv_to_ub(b_bv, n, lbw_bv, ubw_bv, lbwmax_bv, ubwmax_bv, numrotsv, &
       ksv, csv, ssv, &
       b_ub, lbw_ub, ubw_ub, lbwmax_ub, ubwmax_ub, numrotsu, jsu, csu, ssu, &
       cs_sw(:,2:3), ss_sw(:,2:3), error)

    call f_c_convert_ub_to_bv(b_ub, n, lbw_ub, ubw_ub, lbwmax_ub, ubwmax_ub, numrotsu, &
         jsu, csu, ssu, &
         b_bv, lbw_bv, ubw_bv, lbwmax_bv, ubwmax_bv, numrotsv, ksv, csv, ssv, error)

    do j=1,n
       x=b_bv(j,lbw_bv+1)
       d(j)=x/abs(x)
       b_bv(j,1:lbw_bv+ubw_bv+1)=b_bv(j,1:lbw_bv+ubw_bv+1)*conjg(d(j))
    end do
    
    v_tmp=conjg(v)
    call f_c_forward_solve_bv(v, b_bv, n, lbw_bv, ubw_bv, &
         lbwmax_bv, ubwmax_bv, numrotsv, ksv, csv, ssv, &
         v_tmp, 1, error)
    v=conjg(v)

    lbw_ub=0; ubw_ub=0; numrotsu=0
    do j=1,n
       b_ub(1,j)=d(j)
    end do

    call f_c_sweeps1_times_ub(cs_sw(:,2:3), ss_sw(:,2:3), 2, 3, n, &
       b_ub, lbw_ub, ubw_ub, lbwmax_ub, ubwmax_ub, numrotsu, jsu, csu, ssu, &
       b_bv, lbw_bv, ubw_bv, lbwmax_bv, ubwmax_bv, numrotsv, ksv, csv, ssv, error)

    call f_c_convert_bv_to_ub(b_bv, n, lbw_bv, ubw_bv, lbwmax_bv, ubwmax_bv, numrotsv, ksv, &
         csv, ssv, &
         b_ub, lbw_ub, ubw_ub, lbwmax_ub, ubwmax_ub, numrotsu, &
         jsu, csu, ssu, error)
         
    v_tmp=v

    call f_c_e1v_update_ub_to_bv(b_ub, n, lbw_ub, ubw_ub, lbwmax_ub, ubwmax_ub, numrotsu, &
         jsu, csu, ssu, v_tmp, &
         b_bv, lbw_bv, ubw_bv, lbwmax_bv, ubwmax_bv, numrotsv, ksv, csv, ssv, error)

    b_bv(1,lbw_bv+1)=b_bv(1,lbw_bv+1) + u_tmp(1)*conjg(v_tmp(1))
    b_bv(1,lbw_bv+2)=b_bv(1,lbw_bv+2) + u_tmp(1)*conjg(v_tmp(2))

    call f_c_trp_sweeps1_times_bv(cs_sw(:,1), ss_sw(:,1), n, 1, 1, &
       b_bv, lbw_bv, ubw_bv, lbwmax_bv, ubwmax_bv, numrotsv, ksv, csv, ssv, &
       b_ub, lbw_ub, ubw_ub, lbwmax_ub, ubwmax_ub, numrotsu, jsu, csu, ssu, error)

    lbw_ub=1
    call f_c_convert_ub_to_bv(b_ub, n, lbw_ub, ubw_ub, lbwmax_ub, ubwmax_ub, numrotsu, &
         jsu, csu, ssu, &
         b_bv, lbw_bv, ubw_bv, lbwmax_bv, ubwmax_bv, numrotsv, ksv, csv, ssv, error)

    ! if (ubw_ub > 3) then
    !    call f_c_recompress_ub_to_bv(b_ub, n, lbw_ub, ubw_ub, lbwmax_ub, ubwmax_ub, numrotsu, &
    !         jsu, csu, ssu, &
    !         b_bv, lbw_bv, ubw_bv, lbwmax_bv, ubwmax_bv, numrotsv, ksv, csv, ssv, &
    !         0.0_dp, 0.0_dp, ubw_ub-3, error)
    ! else
    !    call f_c_convert_ub_to_bv(b_ub, n, lbw_ub, ubw_ub, lbwmax_ub, ubwmax_ub, numrotsu, &
    !         jsu, csu, ssu, &
    !         b_bv, lbw_bv, ubw_bv, lbwmax_bv, ubwmax_bv, numrotsv, ksv, csv, ssv, error)
    ! end if

    do j=2,n
       if (zeros_i(j)==1) then
          b_bv(j,1)=(0.0_dp, 0.0_dp)
       end if
    end do

  end subroutine f_r1_reorth
  

  ! Errors:
  ! 1: n<1
  ! 2: insufficient storage in bv.
  ! 3: ub%n /= bv%n
  ! 4: Not Hessenberg: bv%lbw /= 1
  ! 5: insufficient storage in ub

  subroutine ss_qr_iteration(bv,ub,shifts, sw, error)
    type(c_bv), intent(inout) :: bv
    type(c_ub) :: ub
    type(c_shift), dimension(:), intent(in) :: shifts
    type(error_info), intent(out) :: error
    type(c_sweeps1), intent(inout) :: sw

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
    if (get_n(bv) /= n .or. size(shifts) /= n .or. get_n(sw) /= n) then
       call set_error(error, 3, id_ss_qr_iteration); return
    end if
    if (bv%lbw /= 1) then
       call set_error(error, 4, id_ss_qr_iteration); return
    end if
    if (get_lbwmax(ub) < 1 .or. get_ubwmax(ub) < bv%ubw) then
       call set_error(error, 5, id_ss_qr_iteration); return
    end if

    sw%numsweeps1=1
    shifts_i=0
    do j=1,n
       shifts_c(j)=shifts(j)%shift
       if (shifts(j)%flag .eqv. .true.) then
          shifts_i(j)=1
       end if
    end do
    call f_ss_qr_iteration(bv%br, get_n(bv), bv%lbw, bv%ubw, get_lbwmax(bv), get_ubwmax(bv), &
         bv%numrotsv, bv%ksv, bv%csv, bv%ssv, & 
         ub%bc, ub%lbw, ub%ubw, get_lbwmax(ub), get_ubwmax(ub), ub%numrotsu, ub%jsu, &
         ub%csu, ub%ssu, shifts_i, shifts_c, sw%cs(:,1), sw%ss(:,1), error)
  end subroutine ss_qr_iteration

  ! A single QR iteration: Start with a B.V decompositions and apply shifts in shifts_c starting
  ! in positions given in shifts_i.   Result is a U.B decomposition and the cosines and sines for Q.
  subroutine f_ss_qr_iteration(b_bv,n, lbw_bv, ubw_bv, lbwmax_bv, ubwmax_bv, numrotsv, ksv, &
       csv, ssv, b_ub, lbw_ub, ubw_ub, lbwmax_ub, ubwmax_ub, numrotsu, jsu, csu, &
       ssu, shifts_i, shifts_c, cs_q, ss_q, error)
    integer(kind=int32), intent(in) :: n, lbwmax_bv, ubwmax_bv, lbwmax_ub, ubwmax_ub
    complex(kind=dp), dimension(n,lbwmax_bv+ubwmax_bv+1), intent(inout) :: b_bv
    integer(kind=int32), dimension(n), intent(in) :: numrotsv
    integer(kind=int32), dimension(n,ubwmax_bv), intent(in) :: ksv
    complex(kind=dp), dimension(n,ubwmax_bv), intent(in) :: csv, ssv
    integer(kind=int32), intent(inout) :: lbw_bv, ubw_bv

    complex(kind=dp), dimension(lbwmax_ub+ubwmax_ub+1,n), intent(out) :: b_ub
    integer(kind=int32), intent(out) :: lbw_ub, ubw_ub
    integer(kind=int32), dimension(n), intent(out) :: numrotsu
    integer(kind=int32), dimension(ubwmax_ub,n), intent(out) :: jsu
    complex(kind=dp), dimension(ubwmax_ub,n), intent(out) :: csu, ssu
    complex(kind=dp), dimension(n), intent(out) :: cs_q, ss_q

    integer(kind=int32), dimension(n) :: shifts_i
    complex(kind=dp), dimension(n) :: shifts_c

    type(error_info), intent(out) :: error

    integer(kind=int32) :: j, k, k0, k1, dubw, dubw_tmp, dlbw, dlbw_tmp
    type(c_rotation) :: rot

    call clear_error(error)

    if (n==1) then
       b_ub(1,1)=b_bv(1,1);
       lbw_ub=0; ubw_ub=0; numrotsu=0; return
    end if

    dubw=1; dubw_tmp=1
    dlbw=0; dlbw_tmp=1

    call f_bw_expand_br(b_bv, n, lbw_bv, ubw_bv, lbwmax_bv, ubwmax_bv, &
         dlbw, dlbw_tmp, dubw, dubw_tmp, lbw_ub, ubw_ub)

    numrotsu=0
    ssu(1:ubw_ub,:)=(0.0_dp, 0.0_dp)
    csu(1:ubw_ub,:)=(0.0_dp, 0.0_dp)
    jsu(1:ubw_ub,:)=0
    b_ub(1:lbw_ub+ubw_ub+1,:)=(0.0_dp, 0.0_dp)

    cs_q=(1.0_dp,0.0_dp)
    ss_q=(0.0_dp,0.0_dp)

    do k=1,n-1
       ! apply v_{k}
       do j=1,numrotsv(k)
          rot%cosine=csv(k,j); rot%sine=ssv(k,j)
          call tbr_times_rotation(b_bv,n,lbw_bv,ubw_bv,0,n-k,trp_rot(rot),ksv(k,j))
       end do
       if (shifts_i(k)==1) then
          ! introduce a new shift.
          rot=lgivens(get_el_br(b_bv,lbw_bv,k,k)-shifts_c(k),get_el_br(b_bv,lbw_bv,k+1,k))
          call rotation_times_tbr(trp_rot(rot),b_bv,n,lbw_bv,ubw_bv,0,0,k)
          call tbr_times_rotation(b_bv,n,lbw_bv,ubw_bv,0,0,rot,k)
          cs_q(k)=rot%cosine; ss_q(k)=rot%sine
       else if (get_el_br(b_bv,lbw_bv,k+1,k-1) /= (0.0_dp, 0.0_dp)) then
          ! chase the bulge.
          rot=lgivens(get_el_br(b_bv,lbw_bv,k,k-1),get_el_br(b_bv,lbw_bv,k+1,k-1))
          call rotation_times_tbr(trp_rot(rot),b_bv,n,lbw_bv,ubw_bv,0,0,k)
          call tbr_times_rotation(b_bv,n,lbw_bv,ubw_bv,0,0,rot,k)
          cs_q(k)=rot%cosine; ss_q(k)=rot%sine
       end if
       ! columns that have a nonzero in superdiagonal ubw
       ! Three cases:
       ! 1. For k <= ubw-1, columns ubw+1, ..., k+ubw-1 might
       !    have a nonzero in superdiagonal ubw.
       ! 2. For ubw <= k <= n-ubw+1, columns k+1, .... k+ubw-1
       !    might have a nonzero in superdiagonal ubw.
       ! 3. For n-ubw+2 <= k <= n-1, columns k+1, ..., n
       !    might have a nonzero in superdiagonal ubw.
       k0=max(k+1,ubw_bv+1)
       k1=min(k+ubw_bv-1,n)
       numrotsu(k)=max(k1-k0+1,0)
       do j=k1,k0,-1
          rot=lgivens2(get_el_br(b_bv,lbw_bv,j-ubw_bv,j), get_el_br(b_bv,lbw_bv,j-ubw_bv+1,j))
          jsu(j-k0+1,k)=j-ubw_bv
          csu(j-k0+1,k)=rot%cosine; ssu(j-k0+1,k)=rot%sine
          call rotation_times_tbr(trp_rot(rot),b_bv,n,lbw_bv,ubw_bv,k,0,j-ubw_bv)
       end do
    end do

    call f_bw_contract_br(b_bv, n, lbw_bv, ubw_bv, lbwmax_bv, ubwmax_bv, &
         dlbw, dlbw_tmp, dubw, dubw_tmp)
    call br_to_bc(b_bv,b_ub,lbw_ub,ubw_ub)
  end subroutine f_ss_qr_iteration

  ! These function are destructive.  bv is modified so that it
  ! is no longer really a nested bv decomposition.  The
  ! superdiagonal of A is revealed in bv%br.
  subroutine c_bv_reveal_superdiag(bv)
    type(c_bv) :: bv
    call f_c_bv_reveal_superdiag(bv%br,get_n(bv),bv%lbw,bv%ubw,get_lbwmax(bv),get_ubwmax(bv), &
         bv%numrotsv, bv%ksv, bv%csv, bv%ssv)
  end subroutine c_bv_reveal_superdiag

  ! Apply v_{n-j} to row j to reveal the first superdiagonal element in each row.
  ! requires on extra superdiagonal.
  subroutine f_c_bv_reveal_superdiag(b_bv, n, lbw_bv, ubw_bv, lbwmax_bv, ubwmax_bv, numrotsv, ksv, &
       csv, ssv)
    integer(kind=int32), intent(in) :: n, lbw_bv, ubw_bv, lbwmax_bv, ubwmax_bv
    complex(kind=dp), dimension(n,lbwmax_bv+ubwmax_bv+1), intent(inout) :: b_bv
    integer(kind=int32), dimension(n), intent(in) :: numrotsv
    integer(kind=int32), dimension(n,ubwmax_bv), intent(in) :: ksv
    complex(kind=dp), dimension(n,ubwmax_bv), intent(in) :: csv, ssv

    integer(kind=int32) :: j, k
    type(c_rotation) :: rot

    do j=1,n-1
       do k=1,numrotsv(j)
          rot%cosine=csv(j,k); rot%sine=ssv(j,k)
          call tbr_times_rotation(b_bv,n,lbw_bv,ubw_bv,j-1,n-j,trp_rot(rot),ksv(j,k))
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
    w=sqrt(tr**2-4.0_dp*(a*d-b*c))
    x=(tr+w)/2.0_dp
    y=(tr-w)/2.0_dp
    if (abs(d-y) < abs(d-x)) then
       lambda=y
    else
       lambda=x
    end if

  end function c_wilkinson_shift

end module mod_qr_iteration
