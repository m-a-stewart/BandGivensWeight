module mod_qr_factorization
  use mod_prec
  use mod_error_id
  use mod_convert_bv_to_ub
  use mod_sweeps
  use mod_shift
  use mod_rotation
  use mod_band_types
  use mod_orth_band_types
  implicit none

  private

  public d_qr, c_qr
  
  public qr_bv_to_ub, d_qr_bv_to_ub, c_qr_bv_to_ub, &
       f_qr_bv_to_ub, f_d_qr_bv_to_ub, f_c_qr_bv_to_ub, &
       d_qr_of, c_qr_of, qr_of

  interface qr_of
     module procedure d_qr_of, c_qr_of
  end interface qr_of

  interface qr_bv_to_ub
     module procedure d_qr_bv_to_ub, c_qr_bv_to_ub
  end interface qr_bv_to_ub

  interface f_qr_bv_to_ub
     module procedure f_d_qr_bv_to_ub, f_c_qr_bv_to_ub
  end interface f_qr_bv_to_ub

  type d_qr
     type(d_sweeps), allocatable :: sw
     type(d_ub), allocatable :: ub
  end type d_qr

  type c_qr
     type(c_sweeps), allocatable :: sw
     type(c_ub), allocatable :: ub
  end type c_qr

contains

  function d_qr_of(bv,error) result(swub)
    type(d_qr) :: swub
    type(d_bv), intent(in) :: bv
    type(error_info), intent(inout), optional :: error
    type(routine_info), parameter :: info=info_d_qr_of
    integer(kind=int32) :: n, lbwmax, ubwmax, lbw, ubw
    type(d_bv), allocatable :: bv1

    call clear_error(error)
    call push_id(info,error)
    n=get_n(bv)
    lbwmax=get_lbwmax(bv);
    ubwmax=get_ubwmax(bv);
    lbw=bv%lbw
    ubw=bv%ubw
    swub%ub=d_new_ub(n,0,min(lbw+ubw,n-1))
    bv1=d_new_bv(n,lbw,min(lbw+ubw+1,n-1))
    call copy_bv(bv,bv1)
    swub%sw=d_new_sweeps(n, lbw+1, n+lbw-1, lbw)
    call d_qr_bv_to_ub(bv1, swub%ub, swub%sw,error)
    deallocate(bv1)
    call pop_id(error)
  end function d_qr_of
  
  ! Errors
  ! 0: no error
  ! 1: ub%n /= bv%n .or. sw%n /= ub%n
  ! 2: Not enough storage for each sweep.
  ! 3: Not enough storage for the number of sweeps.
  ! 4: Not enough storage in bv.
  ! 5: Not enough storage in ub.

  subroutine d_qr_bv_to_ub(bv,ub,sw,error)
    type(d_ub) :: ub
    type(d_bv) :: bv
    type(d_sweeps) :: sw
    type(error_info), intent(inout), optional :: error

    integer(kind=int32) :: lbw, n
    type(routine_info), parameter :: info=info_d_qr_bv_to_ub

    call clear_error(error)
    call push_id(info, error)
    
    lbw=bv%lbw
    n=get_n(bv)
    if (n /= get_n(sw) .or. n /= get_n(ub)) then
       call set_error(1, info, error); return
    end if
    if (get_maxord(sw) < bv%lbw) then
       call set_error(2, info, error); return
    end if
    if (get_maxind(sw) < n+lbw-1 .or. get_minind(sw) > lbw+1) then
       call set_error(3, info, error); return
    end if
    if (get_ubwmax(bv) < min(bv%lbw+bv%ubw+1,n-1) .or. &
         get_lbwmax(bv) < min(bv%lbw,n-1)) then
       call set_error(4, info, error); return
    end if
    if (get_ubwmax(ub) < min(bv%lbw+bv%ubw,n-1)) then
       call set_error(5, info, error); return
    end if
    call f_d_qr_bv_to_ub(bv%br, get_n(bv), bv%lbw, bv%ubw, get_lbwmax(bv), &
         get_ubwmax(bv), bv%numrotsv, bv%ksv, bv%csv, bv%ssv, & 
         ub%bc, ub%lbw, ub%ubw, get_lbwmax(ub), get_ubwmax(ub), ub%numrotsu, ub%jsu, &
         ub%csu, ub%ssu, &
         sw%left, sw%right, sw%inc, get_minind(sw), get_maxind(sw), get_maxord(sw), &
         sw%numrots, sw%js, &
         sw%cs, sw%ss)
    call pop_id(error)
  end subroutine d_qr_bv_to_ub

  subroutine f_d_qr_bv_to_ub(b_bv, n, lbw_bv, ubw_bv, lbwmax_bv, ubwmax_bv, numrotsv, &
       ksv, csv, ssv, &
       b_ub, lbw_ub, ubw_ub, lbwmax_ub, ubwmax_ub, numrotsu, jsu, csu, ssu, &
       leftq, rightq, incq, minind, maxind, maxord, numrotsq, jsq, csq, ssq)
    integer(kind=int32), intent(in) :: n, lbw_bv, ubw_bv, lbwmax_ub, ubwmax_ub, &
         lbwmax_bv, ubwmax_bv
    real(kind=dp), dimension(n,lbwmax_bv+ubwmax_bv+1), intent(inout) :: b_bv
    integer(kind=int32), dimension(n), intent(inout) :: numrotsv
    integer(kind=int32), dimension(n,ubwmax_bv), intent(inout) :: ksv
    real(kind=dp), dimension(n,ubwmax_bv), intent(inout) :: csv, ssv

    real(kind=dp), dimension(lbwmax_ub+ubwmax_ub+1,n), intent(out) :: b_ub
    integer(kind=int32), dimension(n), intent(out) :: numrotsu
    integer(kind=int32), dimension(ubwmax_ub,n), intent(out) :: jsu
    real(kind=dp), dimension(ubwmax_ub,n), intent(out) :: csu, ssu
    integer(kind=int32), intent(out) :: lbw_ub, ubw_ub

    integer(kind=int32), intent(in) :: minind, maxind, maxord
    real(kind=dp), dimension(maxord, minind:maxind), intent(out) :: csq, ssq
    integer(kind=int32), dimension(minind:maxind), intent(out) :: numrotsq
    integer(kind=int32), dimension(maxord,minind:maxind), intent(out) :: jsq
    integer(kind=int32), intent(out) :: leftq, rightq, incq

    integer(kind=int32) :: j, k, lbw, ubw, lbw1, ubw1, k0, k1
    type(d_rotation) :: rot

    if (n == 1) then
       b_ub(1,1)=b_bv(1,1);
       lbw_ub=0; ubw_ub=0; numrotsu=0; 
       numrotsq=0; leftq=0; rightq=-1; incq=1
       return
    end if

    if (lbw_bv == 0) then
       call f_d_convert_bv_to_ub(b_bv, n, lbw_bv, ubw_bv, lbwmax_bv, ubwmax_bv, &
            numrotsv, ksv, csv, ssv, b_ub, lbw_ub, ubw_ub, lbwmax_ub, ubwmax_ub, &
            numrotsu, jsu, csu, ssu)
       numrotsq=0; leftq=0; rightq=-1; incq=1
       return
    end if

    ubw=min(n-1,ubw_bv+lbw_bv)
    lbw=lbw_bv
    lbw1=lbw_bv
    b_bv(:,lbw_bv+ubw_bv+2:lbw+ubw+1)=0.0_dp
    if (ubw < n-1) then
       ubw1=ubw+1
       b_bv(:,lbw1+ubw1+1)=0.0_dp
    else
       ubw1=ubw
    end if

    numrotsu=0
    csu=0.0_dp; ssu=0.0_dp; jsu=0
    numrotsq=0
    leftq=lbw_bv+1; rightq=n+lbw_bv-1; incq=1
    csq=0.0_dp; ssq=0.0_dp; jsq=0
    
    ! Apply V_1, ..., V_{lbw_bv-1}
    do k=1,lbw_bv-1
       do j=1,numrotsv(k)
          rot%cosine=csv(k,j); rot%sine=ssv(k,j)
          call tbr_times_rotation(b_bv,n,lbw1,ubw1,0,n-k,trp_rot(rot),ksv(k,j))
       end do
    end do

    do k=lbw_bv, n-2
       ! Apply V_k^H
       do j=1,numrotsv(k)
          rot%cosine=csv(k,j); rot%sine=ssv(k,j)
          call tbr_times_rotation(b_bv,n,lbw1,ubw1,0,n-k,trp_rot(rot),ksv(k,j))
       end do
       ! Zero the subdiagonal elements in column k+1-lbw_bv using Q_{k+1}
       ! (rotation stored in csq(:,k+1-lbw_bv), etc.)
       k0=k+1-lbw_bv
       numrotsq(k+1)=lbw_bv
       do j=lbw_bv,1,-1
          rot=lgivens(get_el_br(b_bv, lbw1, k0+j-1, k0), &
               get_el_br(b_bv, lbw1, k0+j, k0))
          jsq(j,k+1)=k0+j-1
          csq(j,k+1)=rot%cosine; ssq(j,k+1)=rot%sine
          call rotation_times_tbr(trp_rot(rot), b_bv, n, lbw1, ubw1, 0, 0, k0+j-1)
       end do
       ! Now eliminate the superdiagonal elements introduced by V_k^H
       ! using U_{k+1-lbw_bv}
       k0=max(k+2,lbw_bv+ubw_bv+2)
       k1=min(n,k+ubw_bv+1)
       numrotsu(k+1-lbw_bv)=max(k1-k0+1,0)
       do j=k1,k0,-1
          rot=lgivens2(get_el_br(b_bv,lbw1,j-ubw1,j), get_el_br(b_bv,lbw1,j-ubw1+1,j))
          jsu(j-k0+1,k+1-lbw_bv)=j-ubw1
          csu(j-k0+1,k+1-lbw_bv)=rot%cosine; ssu(j-k0+1,k+1-lbw_bv)=rot%sine
          call rotation_times_tbr(trp_rot(rot),b_bv,n,lbw1,ubw1,k+1-lbw_bv,0,j-ubw1)
       end do
    end do

    ! Only need to apply Q_{k+1} for the last steps
    do k=n-1, n+lbw_bv-2
       k0=k+1-lbw_bv
       numrotsq(k+1)=n-k0
       do j=n-k0,1,-1
          rot=lgivens(get_el_br(b_bv, lbw1, k0+j-1, k0), &
               get_el_br(b_bv, lbw1, k0+j, k0))
          jsq(j,k+1)=k0+j-1
          csq(j,k+1)=rot%cosine; ssq(j,k+1)=rot%sine
          call rotation_times_tbr(trp_rot(rot), b_bv, n, lbw1, ubw1, 0, 0, k0+j-1)
       end do
    end do

    lbw_ub=0; ubw_ub=ubw
    call shift(b_bv,0,-lbw)
    call br_to_bc(b_bv,b_ub,0,ubw)
  end subroutine f_d_qr_bv_to_ub

  function c_qr_of(bv,error) result(swub)
    type(c_qr) :: swub
    type(c_bv), intent(in) :: bv
    type(error_info), intent(inout), optional :: error
    type(routine_info), parameter :: info=info_c_qr_of
    integer(kind=int32) :: n, lbwmax, ubwmax, lbw, ubw
    type(c_bv), allocatable :: bv1

    call clear_error(error)
    call push_id(info,error)
    n=get_n(bv)
    lbwmax=get_lbwmax(bv);
    ubwmax=get_ubwmax(bv);
    lbw=bv%lbw
    ubw=bv%ubw
    swub%ub=c_new_ub(n,0,min(lbw+ubw,n-1))
    bv1=c_new_bv(n,lbw,min(lbw+ubw+1,n-1))
    call copy_bv(bv,bv1)
    swub%sw=c_new_sweeps(n, lbw+1, n+lbw-1, lbw)
    call c_qr_bv_to_ub(bv1, swub%ub, swub%sw,error)
    deallocate(bv1)
    call pop_id(error)
  end function c_qr_of
  

  ! Errors
  ! 0: no error
  ! 1: ub%n /= bv%n .or. sw%n /= ub%n
  ! 2: Not enough storage for each sweep.
  ! 3: Not enough storage for the number of sweeps.
  ! 4: Not enough storage in bv.
  ! 5: Not enough storage in ub.

  subroutine c_qr_bv_to_ub(bv,ub,sw,error)
    type(c_ub) :: ub
    type(c_bv) :: bv
    type(c_sweeps) :: sw
    type(error_info), intent(inout), optional :: error
    type(routine_info), parameter :: info=info_c_qr_bv_to_ub

    integer(kind=int32) :: lbw, n

    call clear_error(error)
    call push_id(info, error)

    lbw=bv%lbw
    n=get_n(bv)
    if (n /= get_n(sw) .or. n /= get_n(ub)) then
       call set_error(1, info, error); return
    end if
    if (get_maxord(sw) < bv%lbw) then
       call set_error(2, info, error); return
    end if
    if (get_maxind(sw)  < n+lbw-1 .or. get_minind(sw) > lbw+1) then
       call set_error(3, info, error); return
    end if
    if (get_ubwmax(bv) < min(bv%lbw+bv%ubw+1,n-1) .or. &
         get_lbwmax(bv) < min(bv%lbw,n-1)) then
       call set_error(4, info, error); return
    end if
    if (get_ubwmax(ub) < min(bv%lbw+bv%ubw,n-1)) then
       call set_error(5, info, error); return
    end if
    call f_c_qr_bv_to_ub(bv%br, get_n(bv), bv%lbw, bv%ubw, get_lbwmax(bv), &
         get_ubwmax(bv), bv%numrotsv, bv%ksv, bv%csv, bv%ssv, & 
         ub%bc, ub%lbw, ub%ubw, get_lbwmax(ub), get_ubwmax(ub), ub%numrotsu, &
         ub%jsu, ub%csu, ub%ssu, sw%left, sw%right, sw%inc, get_minind(sw), &
         get_maxind(sw), get_maxord(sw), sw%numrots, sw%js, &
         sw%cs, sw%ss)
    call pop_id(error)
  end subroutine c_qr_bv_to_ub

  subroutine f_c_qr_bv_to_ub(b_bv, n, lbw_bv, ubw_bv, lbwmax_bv, ubwmax_bv, numrotsv, &
       ksv, csv, ssv, &
       b_ub, lbw_ub, ubw_ub, lbwmax_ub, ubwmax_ub, numrotsu, jsu, csu, ssu, &
       leftq, rightq, incq, minind, maxind, maxord, numrotsq, jsq, csq, ssq)
    integer(kind=int32), intent(in) :: n, lbw_bv, ubw_bv, lbwmax_ub, ubwmax_ub, &
         lbwmax_bv, ubwmax_bv
    complex(kind=dp), dimension(n,lbwmax_bv+ubwmax_bv+1), intent(inout) :: b_bv
    integer(kind=int32), dimension(n), intent(inout) :: numrotsv
    integer(kind=int32), dimension(n,ubwmax_bv), intent(inout) :: ksv
    real(kind=dp), dimension(n,ubwmax_bv), intent(inout) :: csv
    complex(kind=dp), dimension(n,ubwmax_bv), intent(inout) :: ssv

    complex(kind=dp), dimension(lbwmax_ub+ubwmax_ub+1,n), intent(out) :: b_ub
    integer(kind=int32), dimension(n), intent(out) :: numrotsu
    integer(kind=int32), dimension(ubwmax_ub,n), intent(out) :: jsu
    real(kind=dp), dimension(ubwmax_ub,n), intent(out) :: csu
    complex(kind=dp), dimension(ubwmax_ub,n), intent(out) :: ssu
    integer(kind=int32), intent(out) :: lbw_ub, ubw_ub

    integer(kind=int32), intent(in) :: minind, maxind, maxord
    real(kind=dp), dimension(maxord, minind:maxind), intent(out) :: csq
    complex(kind=dp), dimension(maxord, minind:maxind), intent(out) :: ssq
    integer(kind=int32), dimension(minind:maxind), intent(out) :: numrotsq
    integer(kind=int32), dimension(maxord,minind:maxind), intent(out) :: jsq
    integer(kind=int32), intent(out) :: leftq, rightq, incq

    integer(kind=int32) :: j, k, lbw, ubw, lbw1, ubw1, k0, k1
    type(c_rotation) :: rot

    if (n == 1) then
       b_ub(1,1)=b_bv(1,1);
       lbw_ub=0; ubw_ub=0; numrotsu=0; 
       numrotsq=0; leftq=0; rightq=-1; incq=1
       return
    end if

    if (lbw_bv == 0) then
       call f_c_convert_bv_to_ub(b_bv, n, lbw_bv, ubw_bv, lbwmax_bv, ubwmax_bv, &
            numrotsv, ksv, csv, ssv, b_ub, lbw_ub, ubw_ub, lbwmax_ub, ubwmax_ub, &
            numrotsu, jsu, csu, ssu)
       numrotsq=0; leftq=0; rightq=-1; incq=1
       return
    end if

    ubw=min(n-1,ubw_bv+lbw_bv)
    lbw=lbw_bv
    lbw1=lbw_bv
    b_bv(:,lbw_bv+ubw_bv+2:lbw+ubw+1)=(0.0_dp, 0.0_dp)
    if (ubw < n-1) then
       ubw1=ubw+1
       b_bv(:,lbw1+ubw1+1)=(0.0_dp, 0.0_dp)
    else
       ubw1=ubw
    end if

    numrotsu=0
    csu=0.0_dp; ssu=(0.0_dp, 0.0_dp); jsu=0
    numrotsq=0
    leftq=lbw_bv+1; rightq=n+lbw_bv-1; incq=1
    csq=0.0_dp; ssq=(0.0_dp, 0.0_dp); jsq=0
    
    ! Apply V_1, ..., V_{lbw_bv-1}
    do k=1,lbw_bv-1
       do j=1,numrotsv(k)
          rot%cosine=csv(k,j); rot%sine=ssv(k,j)
          call tbr_times_rotation(b_bv,n,lbw1,ubw1,0,n-k,trp_rot(rot),ksv(k,j))
       end do
    end do

    do k=lbw_bv, n-2
       ! Apply V_k^H
       do j=1,numrotsv(k)
          rot%cosine=csv(k,j); rot%sine=ssv(k,j)
          call tbr_times_rotation(b_bv,n,lbw1,ubw1,0,n-k,trp_rot(rot),ksv(k,j))
       end do
       ! Zero the subdiagonal elements in column k+1-lbw_bv using Q_{k+1}
       ! (rotation stored in csq(:,k+1-lbw_bv), etc.)
       k0=k+1-lbw_bv
       numrotsq(k+1)=lbw_bv
       do j=lbw_bv,1,-1
          rot=lgivens(get_el_br(b_bv, lbw1, k0+j-1, k0), &
               get_el_br(b_bv, lbw1, k0+j, k0))
          jsq(j,k+1)=k0+j-1
          csq(j,k+1)=rot%cosine; ssq(j,k+1)=rot%sine
          call rotation_times_tbr(trp_rot(rot), b_bv, n, lbw1, ubw1, 0, 0, k0+j-1)
       end do
       ! Now eliminate the superdiagonal elements introduced by V_k^H
       ! using U_{k+1-lbw_bv}
       k0=max(k+2,lbw_bv+ubw_bv+2)
       k1=min(n,k+ubw_bv+1)
       numrotsu(k+1-lbw_bv)=max(k1-k0+1,0)
       do j=k1,k0,-1
          rot=lgivens2(get_el_br(b_bv,lbw1,j-ubw1,j), get_el_br(b_bv,lbw1,j-ubw1+1,j))
          jsu(j-k0+1,k+1-lbw_bv)=j-ubw1
          csu(j-k0+1,k+1-lbw_bv)=rot%cosine; ssu(j-k0+1,k+1-lbw_bv)=rot%sine
          call rotation_times_tbr(trp_rot(rot),b_bv,n,lbw1,ubw1,k+1-lbw_bv,0,j-ubw1)
       end do
    end do

    ! Only need to apply Q_{k+1} for the last steps
    do k=n-1, n+lbw_bv-2
       k0=k+1-lbw_bv
       numrotsq(k+1)=n-k0
       do j=n-k0,1,-1
          rot=lgivens(get_el_br(b_bv, lbw1, k0+j-1, k0), &
               get_el_br(b_bv, lbw1, k0+j, k0))
          jsq(j,k+1)=k0+j-1
          csq(j,k+1)=rot%cosine; ssq(j,k+1)=rot%sine
          call rotation_times_tbr(trp_rot(rot), b_bv, n, lbw1, ubw1, 0, 0, k0+j-1)
       end do
    end do

    lbw_ub=0; ubw_ub=ubw
    call shift(b_bv,0,-lbw)
    call br_to_bc(b_bv,b_ub,0,ubw)
  end subroutine f_c_qr_bv_to_ub

end module mod_qr_factorization

