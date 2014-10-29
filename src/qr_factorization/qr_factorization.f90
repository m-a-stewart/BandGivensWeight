module mod_qr_factorization
  use mod_prec
  use mod_error_id
  use mod_convert_ub_to_bv
  use mod_convert_bv_to_ub
  use mod_sweeps1
  use mod_rotation
  use mod_band_types
  use mod_nested_types
  implicit none

  private

  public :: reduce_lbw_bv_to_ub, d_reduce_lbw_bv_to_ub, c_reduce_lbw_bv_to_ub, &
       f_reduce_lbw_bv_to_ub, f_d_reduce_lbw_bv_to_ub, f_c_reduce_lbw_bv_to_ub

  public qr_bv_to_ub, d_qr_bv_to_ub, c_qr_bv_to_ub, &
       f_qr_bv_to_ub, f_d_qr_bv_to_ub, f_c_qr_bv_to_ub

  public :: info_d_reduce_lbw_bv_to_ub, info_c_reduce_lbw_bv_to_ub, &
       info_d_qr_bv_to_ub, info_c_qr_bv_to_ub

  interface reduce_lbw_bv_to_ub
     module procedure d_reduce_lbw_bv_to_ub, c_reduce_lbw_bv_to_ub
  end interface reduce_lbw_bv_to_ub

  interface f_reduce_lbw_bv_to_ub
     module procedure f_d_reduce_lbw_bv_to_ub, f_c_reduce_lbw_bv_to_ub
  end interface f_reduce_lbw_bv_to_ub

  interface qr_bv_to_ub
     module procedure d_qr_bv_to_ub, c_qr_bv_to_ub
  end interface qr_bv_to_ub

  interface f_qr_bv_to_ub
     module procedure f_d_qr_bv_to_ub, f_c_qr_bv_to_ub
  end interface f_qr_bv_to_ub

  type(routine_info), parameter :: info_d_reduce_lbw_bv_to_ub=routine_info(id_d_reduce_lbw_bv_to_ub, &
       'd_reduce_lbw_bv_to_ub', &
       [ character(len=error_message_length) :: 'ub%n /= bv%n', 'bv%lbw <= 0', 'n < 1', &
       'not enough temp. storage in bv', 'not enough storage in ub', 'dim. of cs or ss /= n' ])

  type(routine_info), parameter :: info_c_reduce_lbw_bv_to_ub=routine_info(id_c_reduce_lbw_bv_to_ub, &
       'c_reduce_lbw_bv_to_ub', &
       [ character(len=error_message_length) :: 'ub%n /= bv%n', 'bv%lbw <= 0', 'n < 1', &
       'not enough temp. storage in bv', 'not enough storage in ub', 'dim. of cs or ss /= n' ])

  type(routine_info), parameter :: info_d_qr_bv_to_ub=routine_info(id_d_qr_bv_to_ub, &
       'd_qr_bv_to_ub', &
       [ character(len=error_message_length) :: 'ub%n /= bv%n or sw%n /= ub%n', &
       'Not enough stroage for sweeps1' ])

  type(routine_info), parameter :: info_c_qr_bv_to_ub=routine_info(id_c_qr_bv_to_ub, &
       'c_qr_bv_to_ub', &
       [ character(len=error_message_length) :: 'ub%n /= bv%n', 'bv%lbw <= 0', 'dim. of cs or ss /= n' ])

contains

  !
  ! Errors
  ! 0: no error
  ! 1: ub%n /= bv%n
  ! 2: bv%lbw <= 0
  ! 3: n < 1
  ! 4: not enough temp in bv
  ! 5: not enough storage in ub
  ! 6: ss or cs are the wrong size.
  !
  subroutine d_reduce_lbw_bv_to_ub(bv,ub,cs,ss,error)
    type(d_ub) :: ub
    type(d_bv) :: bv
    type(error_info), intent(out) :: error
    real(kind=dp), dimension(:) :: cs, ss
    call clear_error(error)
    if (get_n(ub) /= get_n(bv)) then
       call set_error(error, 1, id_d_reduce_lbw_bv_to_ub); return
    end if
    if (bv%lbw <= 0) then
       call set_error(error, 2, id_d_reduce_lbw_bv_to_ub); return
    end if
    if (get_n(ub) < 1) then
       call set_error(error, 3, id_d_reduce_lbw_bv_to_ub); return
    end if
    if (get_ubwmax(bv) < bv%ubw+2 .and. get_ubwmax(bv) < get_n(bv)-1) then
       call set_error(error, 4, id_d_reduce_lbw_bv_to_ub); return
    end if
    if (get_lbwmax(ub) < bv%lbw-1 .or. get_ubwmax(ub) < bv%ubw+1) then
       call set_error(error, 5, id_d_reduce_lbw_bv_to_ub); return
    end if
    if (size(cs) /= get_n(ub) - 1 .or. size(ss) /= get_n(ub) - 1 ) then
       call set_error(error, 6, id_d_reduce_lbw_bv_to_ub); return
    end if

    call f_d_reduce_lbw_bv_to_ub(bv%br, get_n(bv), bv%lbw, bv%ubw, get_lbwmax(bv), &
         get_ubwmax(bv), bv%numrotsv, bv%ksv, bv%csv, bv%ssv, & 
         ub%bc, ub%lbw, ub%ubw, get_lbwmax(ub), get_ubwmax(ub), ub%numrotsu, ub%jsu, ub%csu, ub%ssu, &
         cs, ss, error)
  end subroutine d_reduce_lbw_bv_to_ub

  subroutine f_d_reduce_lbw_bv_to_ub(b_bv, n, lbw_bv, ubw_bv, lbwmax_bv, ubwmax_bv, numrotsv, &
       ksv, csv, ssv, &
       b_ub, lbw_ub, ubw_ub, lbwmax_ub, ubwmax_ub, numrotsu, jsu, csu, ssu, cs, ss, error)
    integer(kind=int32), intent(in) :: n, lbwmax_ub, ubwmax_ub, lbwmax_bv, ubwmax_bv
    real(kind=dp), dimension(n,lbwmax_bv+ubwmax_bv+1), intent(inout) :: b_bv
    integer(kind=int32), dimension(n), intent(in) :: numrotsv
    integer(kind=int32), dimension(n,ubwmax_bv), intent(in) :: ksv
    real(kind=dp), dimension(n,ubwmax_bv), intent(in) :: csv, ssv
    integer(kind=int32), intent(inout) :: lbw_bv, ubw_bv

    real(kind=dp), dimension(lbwmax_ub+ubwmax_ub+1,n), intent(out) :: b_ub
    integer(kind=int32), dimension(n), intent(out) :: numrotsu
    integer(kind=int32), dimension(ubwmax_ub,n), intent(out) :: jsu
    real(kind=dp), dimension(ubwmax_ub,n), intent(out) :: csu, ssu

    real(kind=dp), dimension(n-1) :: cs, ss
    integer(kind=int32), intent(out) :: lbw_ub, ubw_ub
    type(error_info), intent(out) :: error

    integer(kind=int32) :: j, k, k0,k1, dubw, dubw_tmp, dlbw, dlbw_tmp
    type(d_rotation) :: rot

    call clear_error(error)

    dubw=1; dubw_tmp=1
    dlbw=0; dlbw_tmp=0

    call f_bw_expand_br(b_bv, n, lbw_bv, ubw_bv, lbwmax_bv, ubwmax_bv, &
         dlbw, dlbw_tmp, dubw, dubw_tmp, lbw_ub, ubw_ub)

    b_ub(1:lbw_ub+ubw_ub+1,:)=0.0_dp
    numrotsu=0
    ssu(1:ubw_ub,:)=0.0_dp; csu(1:ubw_ub,:)=0.0_dp
    jsu(1:ubw_ub,:)=0

    ss=0.0_dp; cs=1.0_dp

    if (n==1) then
       b_ub(1,1)=b_bv(1,1);
       lbw_ub=0; ubw_ub=0; numrotsu=0; return
    end if

    do k=1,n-1
       ! apply v_{n-k}
       do j=1,numrotsv(k)
          rot%cosine=csv(k,j); rot%sine=ssv(k,j)
          call tbr_times_rotation(b_bv,n,lbw_bv,ubw_bv,0,n-k,trp_rot(rot),ksv(k,j))
       end do
       ! eliminate subdiagonal in row k+1 and column k-lbw_bv+1
       if (k-lbw_bv+1 >= 1 .and. k+1 <= n) then
          rot=lgivens(get_el_br(b_bv,lbw_bv,k,k-lbw_bv+1),get_el_br(b_bv,lbw_bv,k+1,k-lbw_bv+1))
          call rotation_times_tbr(trp_rot(rot),b_bv,n,lbw_bv,ubw_bv,0,0,k)
          cs(k)=rot%cosine; ss(k)=rot%sine
          call set_el_br(b_bv,lbw_bv,k+1,k-lbw_bv+1,0.0_dp)
       end if
       ! Eliminate superdiagonal ubw_bv.
       ! columns that have a nonzero in superdiagonal ubw_bv
       ! Three cases:
       ! 1. For k <= ubw_bv-1, columns ubw_bv+1, ..., k+ubw_bv-1 might
       !    have a nonzero in superdiagonal ubw_bv.
       ! 2. For ubw_bv <= k <= n-ubw_bv+1, columns k+1, .... k+ubw_bv-1
       !    might have a nonzero in superdiagonal ubw_bv.
       ! 3. For n-ubw_bv+2 <= k <= n-1, columns k+1, ..., n
       !    might have a nonzero in superdiagonal ubw_bv.
       k0=max(k+1,ubw_bv+1)
       k1=min(k+ubw_bv-1,n)
       numrotsu(k)=max(k1-k0+1,0)
       do j=k1,k0,-1
          rot=lgivens2(get_el_br(b_bv,lbw_bv,j-ubw_bv,j), get_el_br(b_bv,lbw_bv,j-ubw_bv+1,j))
          jsu(j-k0+1,k)=j-ubw_bv
          csu(j-k0+1,k)=rot%cosine; ssu(j-k0+1,k)=rot%sine
          call rotation_times_tbr(trp_rot(rot),b_bv,n,lbw_bv,ubw_bv,k,0,j-ubw_bv)
          call set_el_br(b_bv,lbw_bv,j-ubw_bv,j, 0.0_dp)
       end do
    end do
    call f_bw_contract_br(b_bv, n, lbw_bv, ubw_bv, lbwmax_bv, ubwmax_bv, &
         dlbw, dlbw_tmp, dubw, dubw_tmp)
    call br_to_bc(b_bv,b_ub,lbw_ub, ubw_ub)
    lbw_ub=lbw_ub-1
  end subroutine f_d_reduce_lbw_bv_to_ub

  subroutine c_reduce_lbw_bv_to_ub(bv,ub,cs,ss,error)
    type(c_ub) :: ub
    type(c_bv) :: bv
    type(error_info), intent(out) :: error
    complex(kind=dp), dimension(:) :: cs, ss
    call clear_error(error)
    if (get_n(ub) /= get_n(bv)) then
       call set_error(error, 1, id_d_reduce_lbw_bv_to_ub); return
    end if
    if (bv%lbw <= 0) then
       call set_error(error, 2, id_d_reduce_lbw_bv_to_ub); return
    end if
    if (get_n(ub) < 1) then
       call set_error(error, 3, id_f_d_reduce_lbw_bv_to_ub); return
    end if
    if (get_ubwmax(bv) < bv%ubw+2 .and. get_ubwmax(bv) < get_n(bv)-1) then
       call set_error(error, 4, id_d_reduce_lbw_bv_to_ub); return
    end if
    if (get_lbwmax(ub) < bv%lbw-1 .or. get_ubwmax(ub) < bv%ubw+1) then
       call set_error(error, 5, id_d_reduce_lbw_bv_to_ub); return
    end if
    if (size(cs) /= get_n(ub) - 1 .or. size(ss) /= get_n(ub) - 1) then
       call set_error(error, 6, id_d_reduce_lbw_bv_to_ub); return
    end if
    call f_c_reduce_lbw_bv_to_ub(bv%br, get_n(bv), bv%lbw, bv%ubw, get_lbwmax(bv), &
         get_ubwmax(bv), bv%numrotsv, bv%ksv, bv%csv, bv%ssv, & 
         ub%bc, ub%lbw, ub%ubw, get_lbwmax(ub), get_ubwmax(ub), ub%numrotsu, ub%jsu, ub%csu, ub%ssu, &
         cs, ss, error)
  end subroutine c_reduce_lbw_bv_to_ub

  ! 1: n < 1
  ! 2: Insufficient temp storage in b_bv
  ! 3: Insufficient storage in b_ub

  subroutine f_c_reduce_lbw_bv_to_ub(b_bv, n, lbw_bv, ubw_bv, lbwmax_bv, ubwmax_bv, numrotsv, &
       ksv, csv, ssv, &
       b_ub, lbw_ub, ubw_ub, lbwmax_ub, ubwmax_ub, numrotsu, jsu, csu, ssu, cs, ss, error)
    integer(kind=int32), intent(in) :: n, lbwmax_ub, ubwmax_ub, lbwmax_bv, ubwmax_bv
    complex(kind=dp), dimension(n,lbwmax_bv+ubwmax_bv+1), intent(inout) :: b_bv
    integer(kind=int32), dimension(n), intent(in) :: numrotsv
    integer(kind=int32), dimension(n,ubwmax_bv), intent(in) :: ksv
    complex(kind=dp), dimension(n,ubwmax_bv), intent(in) :: csv, ssv
    integer(kind=int32), intent(inout) :: lbw_bv, ubw_bv

    complex(kind=dp), dimension(lbwmax_ub+ubwmax_ub+1,n), intent(out) :: b_ub
    integer(kind=int32), dimension(n), intent(out) :: numrotsu
    integer(kind=int32), dimension(ubwmax_ub,n), intent(out) :: jsu
    complex(kind=dp), dimension(ubwmax_ub,n), intent(out) :: csu, ssu

    complex(kind=dp), dimension(n-1) :: cs, ss
    integer(kind=int32), intent(out) :: lbw_ub, ubw_ub
    type(error_info), intent(out) :: error

    integer(kind=int32) :: j, k, k0,k1, dubw, dubw_tmp, dlbw, dlbw_tmp
    type(c_rotation) :: rot

    call clear_error(error)

    dubw=1; dubw_tmp=1
    dlbw=0; dlbw_tmp=0

    call f_bw_expand_br(b_bv, n, lbw_bv, ubw_bv, lbwmax_bv, ubwmax_bv, &
         dlbw, dlbw_tmp, dubw, dubw_tmp, lbw_ub, ubw_ub)

    b_ub(1:lbw_ub+ubw_ub+1,:)=(0.0_dp,0.0_dp)
    numrotsu=0
    ssu(1:ubw_ub,:)=(0.0_dp,0.0_dp); csu(1:ubw_ub,:)=(0.0_dp,0.0_dp)
    jsu(1:ubw_ub,:)=0

    ss=(0.0_dp,0.0_dp); cs=(1.0_dp,0.0_dp)

    if (n==1) then
       b_ub(1,1)=b_bv(1,1);
       lbw_ub=0; ubw_ub=0; numrotsu=0; return
    end if

    do k=1,n-1
       ! apply v_{n-k}
       do j=1,numrotsv(k)
          rot%cosine=csv(k,j); rot%sine=ssv(k,j)
          call tbr_times_rotation(b_bv,n,lbw_bv,ubw_bv,0,n-k,trp_rot(rot),ksv(k,j))
       end do
       ! eliminate subdiagonal in row k+1 and column k-lbw_bv+1
       if (k-lbw_bv+1 >= 1 .and. k+1 <= n) then
          rot=lgivens(get_el_br(b_bv,lbw_bv,k,k-lbw_bv+1),get_el_br(b_bv,lbw_bv,k+1,k-lbw_bv+1))
          call rotation_times_tbr(trp_rot(rot),b_bv,n,lbw_bv,ubw_bv,0,0,k)
          cs(k)=rot%cosine; ss(k)=rot%sine
          call set_el_br(b_bv,lbw_bv,k+1,k-lbw_bv+1,(0.0_dp,0.0_dp))
       end if
       ! Eliminate superdiagonal ubw_bv.
       ! columns that have a nonzero in superdiagonal ubw_bv
       ! Three cases:
       ! 1. For k <= ubw_bv-1, columns ubw_bv+1, ..., k+ubw_bv-1 might
       !    have a nonzero in superdiagonal ubw_bv.
       ! 2. For ubw_bv <= k <= n-ubw_bv+1, columns k+1, .... k+ubw_bv-1
       !    might have a nonzero in superdiagonal ubw_bv.
       ! 3. For n-ubw_bv+2 <= k <= n-1, columns k+1, ..., n
       !    might have a nonzero in superdiagonal ubw_bv.
       k0=max(k+1,ubw_bv+1)
       k1=min(k+ubw_bv-1,n)
       numrotsu(k)=max(k1-k0+1,0)
       do j=k1,k0,-1
          rot=lgivens2(get_el_br(b_bv,lbw_bv,j-ubw_bv,j), get_el_br(b_bv,lbw_bv,j-ubw_bv+1,j))
          jsu(j-k0+1,k)=j-ubw_bv
          csu(j-k0+1,k)=rot%cosine; ssu(j-k0+1,k)=rot%sine
          call rotation_times_tbr(trp_rot(rot),b_bv,n,lbw_bv,ubw_bv,k,0,j-ubw_bv)
          call set_el_br(b_bv,lbw_bv,j-ubw_bv,j, (0.0_dp,0.0_dp))
       end do
    end do
    call f_bw_contract_br(b_bv, n, lbw_bv, ubw_bv, lbwmax_bv, ubwmax_bv, &
         dlbw, dlbw_tmp, dubw, dubw_tmp)
    call br_to_bc(b_bv,b_ub,lbw_ub, ubw_ub)
    lbw_ub=lbw_ub-1
  end subroutine f_c_reduce_lbw_bv_to_ub

  !
  ! Errors
  ! 0: no error
  ! 1: ub%n /= bv%n .or. sw%n /= ub%n
  ! 2: Not enough storage for sweeps1.
  !
  subroutine d_qr_bv_to_ub(bv,ub,sw,error)
    type(d_ub) :: ub
    type(d_bv) :: bv
    type(d_sweeps1) :: sw
    type(error_info), intent(out) :: error

    integer(kind=int32) :: lbw
    lbw=bv%lbw
    sw%numsweeps1=lbw
    call clear_error(error)
    if (get_n(ub) /= get_n(bv) .or. get_n(ub) /= get_n(sw)) then
       call set_error(error, 1, id_d_qr_bv_to_ub); return
    end if
    if (lbw <= 0) then
       call convert_bv_to_ub(bv,ub,error); return
    end if
    if (get_maxsweeps1(sw) < lbw) then
       call set_error(error, 2, id_d_qr_bv_to_ub); return
    end if
    call f_d_qr_bv_to_ub(bv%br, get_n(bv), bv%lbw, bv%ubw, get_lbwmax(bv), &
         get_ubwmax(bv), bv%numrotsv, bv%ksv, bv%csv, bv%ssv, & 
         ub%bc, ub%lbw, ub%ubw, get_lbwmax(ub), get_ubwmax(ub), ub%numrotsu, ub%jsu, ub%csu, ub%ssu, &
         sw%cs(:,1:lbw), sw%ss(:,1:lbw), error)
  end subroutine d_qr_bv_to_ub

  subroutine f_d_qr_bv_to_ub(b_bv, n, lbw_bv, ubw_bv, lbwmax_bv, ubwmax_bv, numrotsv, &
       ksv, csv, ssv, &
       b_ub, lbw_ub, ubw_ub, lbwmax_ub, ubwmax_ub, numrotsu, jsu, csu, ssu, cs, ss, error)
    integer(kind=int32), intent(in) :: n, lbw_bv, ubw_bv, lbwmax_ub, ubwmax_ub, lbwmax_bv, ubwmax_bv
    real(kind=dp), dimension(n,lbwmax_bv+ubwmax_bv+1), intent(inout) :: b_bv
    integer(kind=int32), dimension(n), intent(inout) :: numrotsv
    integer(kind=int32), dimension(n,ubwmax_bv), intent(inout) :: ksv
    real(kind=dp), dimension(n,ubwmax_bv), intent(inout) :: csv, ssv

    real(kind=dp), dimension(lbwmax_ub+ubwmax_ub+1,n), intent(out) :: b_ub
    integer(kind=int32), dimension(n), intent(out) :: numrotsu
    integer(kind=int32), dimension(ubwmax_ub,n), intent(out) :: jsu
    real(kind=dp), dimension(ubwmax_ub,n), intent(out) :: csu, ssu

    real(kind=dp), dimension(n-1,lbw_bv) :: cs, ss
    integer(kind=int32), intent(out) :: lbw_ub, ubw_ub
    type(error_info), intent(out) :: error
    integer(kind=int32) :: j, lbw, ubw

    call clear_error(error)
    lbw=lbw_bv; ubw=ubw_bv
    if (n == 1) then
       b_ub(1,1)=b_bv(1,1);
       lbw_ub=0; ubw_ub=0; numrotsu=0; return
    end if
    if (lbw <= 0) then
       call f_convert_bv_to_ub(b_bv, n, lbw, ubw, lbwmax_bv, ubwmax_bv, &
            numrotsv, ksv, csv, ssv, b_ub, lbw_ub, ubw_ub, &
            lbwmax_ub, ubwmax_ub, numrotsu, jsu, &
            csu, ssu, error)
       return
    end if
    do j=1,lbw-1
       call f_d_reduce_lbw_bv_to_ub(b_bv, n, lbw, ubw, lbwmax_bv, ubwmax_bv, numrotsv, &
            ksv, csv, ssv, &
            b_ub, lbw_ub, ubw_ub, lbwmax_ub, ubwmax_ub, numrotsu, jsu, &
            csu, ssu, cs(:,j), ss(:,j), error)
       call f_d_convert_ub_to_bv(b_ub, n, lbw_ub, ubw_ub, lbwmax_ub, ubwmax_ub, numrotsu, &
            jsu, csu, ssu, b_bv, lbw, ubw, lbwmax_bv, ubwmax_bv, numrotsv, ksv, &
            csv, ssv, error)
    end do
    call f_d_reduce_lbw_bv_to_ub(b_bv, n, lbw, ubw, lbwmax_bv, ubwmax_bv, numrotsv, &
         ksv, csv, ssv, &
         b_ub, lbw_ub, ubw_ub, lbwmax_ub, ubwmax_ub, numrotsu, jsu, &
         csu, ssu, cs(:,j), ss(:,j), error)
  end subroutine f_d_qr_bv_to_ub


  subroutine c_qr_bv_to_ub(bv,ub,sw,error)
    type(c_ub) :: ub
    type(c_bv) :: bv
    type(c_sweeps1) :: sw
    type(error_info), intent(out) :: error

    integer(kind=int32) :: lbw
    lbw=bv%lbw
    sw%numsweeps1=lbw
    call clear_error(error)
    if (get_n(ub) /= get_n(bv) .or. get_n(ub) /= get_n(sw)) then
       call set_error(error, 1, id_d_qr_bv_to_ub); return
    end if
    if (get_n(ub) < 1) then
       call set_error(error, 2, id_f_c_qr_bv_to_ub); return
    end if
    if (lbw <= 0) then
       call convert_bv_to_ub(bv,ub,error); return
    end if
    if (get_maxsweeps1(sw) < lbw) then
       call set_error(error, 2, id_d_qr_bv_to_ub); return
    end if
    call f_c_qr_bv_to_ub(bv%br, get_n(bv), bv%lbw, bv%ubw, get_lbwmax(bv), &
         get_ubwmax(bv), bv%numrotsv, bv%ksv, bv%csv, bv%ssv, & 
         ub%bc, ub%lbw, ub%ubw, get_lbwmax(ub), get_ubwmax(ub), ub%numrotsu, ub%jsu, ub%csu, ub%ssu, &
         sw%cs(:,1:lbw), sw%ss(:,1:lbw), error)
  end subroutine c_qr_bv_to_ub

  subroutine f_c_qr_bv_to_ub(b_bv, n, lbw_bv, ubw_bv, lbwmax_bv, ubwmax_bv, numrotsv, &
       ksv, csv, ssv, &
       b_ub, lbw_ub, ubw_ub, lbwmax_ub, ubwmax_ub, numrotsu, jsu, csu, ssu, cs, ss, error)
    integer(kind=int32), intent(in) :: n, lbw_bv, ubw_bv, lbwmax_ub, ubwmax_ub, lbwmax_bv, ubwmax_bv
    complex(kind=dp), dimension(n,lbwmax_bv+ubwmax_bv+1), intent(inout) :: b_bv
    integer(kind=int32), dimension(n), intent(inout) :: numrotsv
    integer(kind=int32), dimension(n,ubwmax_bv), intent(inout) :: ksv
    complex(kind=dp), dimension(n,ubwmax_bv), intent(inout) :: csv, ssv

    complex(kind=dp), dimension(lbwmax_ub+ubwmax_ub+1,n), intent(out) :: b_ub
    integer(kind=int32), dimension(n), intent(out) :: numrotsu
    integer(kind=int32), dimension(ubwmax_ub,n), intent(out) :: jsu
    complex(kind=dp), dimension(ubwmax_ub,n), intent(out) :: csu, ssu

    complex(kind=dp), dimension(n-1,lbw_bv) :: cs, ss
    integer(kind=int32), intent(out) :: lbw_ub, ubw_ub
    type(error_info), intent(out) :: error
    integer(kind=int32) :: j, lbw, ubw

    call clear_error(error)
    lbw=lbw_bv; ubw=ubw_bv
    if (n == 1) then
       b_ub(1,1)=b_bv(1,1);
       lbw_ub=0; ubw_ub=0; numrotsu=0; return
    end if
    if (lbw <= 0) then
       call f_convert_bv_to_ub(b_bv, n, lbw, ubw, lbwmax_bv, ubwmax_bv, &
            numrotsv, ksv, csv, ssv, b_ub, lbw_ub, ubw_ub, lbwmax_ub, ubwmax_ub, numrotsu, jsu, &
            csu, ssu, error)
       return
    end if
    do j=1,lbw-1
       call f_c_reduce_lbw_bv_to_ub(b_bv, n, lbw, ubw, lbwmax_bv, ubwmax_bv, numrotsv, &
            ksv, csv, ssv, &
            b_ub, lbw_ub, ubw_ub, lbwmax_ub, ubwmax_ub, numrotsu, jsu, csu, ssu, cs(:,j), ss(:,j), error)
       call f_c_convert_ub_to_bv(b_ub, n, lbw_ub, ubw_ub, lbwmax_ub, ubwmax_ub, numrotsu, &
            jsu, csu, ssu, b_bv, lbw, ubw, lbwmax_bv, ubwmax_bv, numrotsv, ksv, &
            csv, ssv, error)
    end do
    call f_c_reduce_lbw_bv_to_ub(b_bv, n, lbw, ubw, lbwmax_bv, ubwmax_bv, numrotsv, &
         ksv, csv, ssv, &
         b_ub, lbw_ub, ubw_ub, lbwmax_ub, ubwmax_ub, numrotsu, jsu, csu, ssu, cs(:,j), ss(:,j), error)
  end subroutine f_c_qr_bv_to_ub

end module mod_qr_factorization

