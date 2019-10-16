module mod_row_compress
  use mod_prec
  use mod_error_id
  use mod_band_types
  use mod_orth_band_types
  use mod_sweeps
  use mod_rotation
  implicit none

  ! Routines that when given a rank structured matrix A parameterized
  ! by a UBT decomposition, compute the BV decomposition of a row
  ! compression A=QC where C has band structure in its lower
  ! triangular part and rank structure in its upper triangular part.

  interface row_compress
     module procedure d_row_compress, z_row_compress
  end interface row_compress

  interface f_row_compress
     module procedure f_d_row_compress, f_z_row_compress
  end interface f_row_compress

  interface rc_of
     module procedure d_rc_of, z_rc_of
  end interface rc_of

  type d_rc
     type(d_sweeps), allocatable :: sw
     type(d_bv), allocatable :: bv
  end type d_rc

  type z_rc
     type(z_sweeps), allocatable :: sw
     type(z_bv), allocatable :: bv
  end type z_rc
  
contains

  function d_rc_of(ubt,error) result(swbv)
    type(d_rc) :: swbv
    type(d_ubt), intent(in) :: ubt
    type(error_info), intent(inout), optional :: error
    type(routine_info), parameter :: info=info_d_rc_of
    integer(kind=int32) :: n, lbwmax, ubwmax, lbw, ubw
    type(d_ubt), allocatable :: ubt1

    if (failure(error)) return
    call push_id(info,error)

    n=get_n(ubt)
    lbwmax=get_lbwmax(ubt);
    ubwmax=get_ubwmax(ubt);
    lbw=ubt%lbw
    ubw=ubt%ubw
    swbv%bv=d_new_bv(n,lbw,min(lbw+ubw,n-1))
    ubt1=d_new_ubt(n,min(lbw+1,n-1),min(lbw+ubw+1,n-1))
    call copy(ubt1,ubt)
    swbv%sw=d_new_sweeps(n, 1, n-2, lbw)
    call d_row_compress(ubt1, swbv%bv, swbv%sw, error)
    deallocate(ubt1)
    call pop_id(error)
  end function d_rc_of

  ! Errors
  ! 0: no error
  ! 1: ub%n /= bv%n .or. sw%n /= ub%n
  ! 2: Not enough storage for each sweep.
  ! 3: Not enough storage for the number of sweeps.
  ! 4: Not enough storage in ubt.
  ! 5: Not enough storage in bv.
  subroutine d_row_compress(ubt,bv,sw,error)
    type(d_ubt) :: ubt
    type(d_bv) :: bv
    type(d_sweeps) :: sw
    type(error_info), intent(inout), optional :: error
    type(routine_info), parameter :: info=info_d_row_compress
    integer(kind=int32) :: n, lbw
    lbw=ubt%lbw
    n=get_n(ubt)
    if (failure(error)) return
    call push_id(info, error)
    if (n /= get_n(bv) .or. n /= get_n(sw)) then
       call set_error(1, info, error); return
    end if
    if (get_maxord(sw) < ubt%lbw) then
       call set_error(2, info, error); return
    end if
    if (get_maxind(sw) < n-2 .or. get_minind(sw) > 1) then
       call set_error(3, info, error); return
    end if
    if (get_ubwmax(ubt) < min(ubt%lbw+ubt%ubw+1,n-1) .or. &
         get_lbwmax(ubt) < min(ubt%lbw+1,n-1)) then
       call set_error(4, info, error); return
    end if
    if (get_lbwmax(bv) < ubt%lbw .or. &
         get_ubwmax(bv) < min(ubt%lbw+ubt%ubw,n-1)) then
       call set_error(5, info, error); return
    end if
    call f_d_row_compress(ubt%bc, get_n(ubt), ubt%lbw, ubt%ubw, get_lbwmax(ubt), &
         get_ubwmax(ubt), ubt%numrotsu, ubt%jsu, ubt%csu, ubt%ssu, & 
         ubt%numrotst, ubt%kst, ubt%cst, ubt%sst, & 
         bv%br, bv%lbw, bv%ubw, get_lbwmax(bv), get_ubwmax(bv), bv%numrotsv, bv%ksv, &
         bv%csv, bv%ssv, sw%left, sw%right, sw%inc, get_minind(sw), get_maxind(sw), &
         get_maxord(sw), sw%numrots, sw%js, sw%cs, sw%ss)
    call pop_id(error)
  end subroutine d_row_compress

  subroutine f_d_row_compress(b_ubt, n, lbw_ubt, ubw_ubt, lbwmax_ubt, ubwmax_ubt, numrotsu, &
       jsu, csu, ssu, numrotst, kst, cst, sst, &
       b_bv, lbw_bv, ubw_bv, lbwmax_bv, ubwmax_bv, numrotsv, ksv, csv, ssv, &
       leftq, rightq, incq, minind, maxind, maxord, numrotsq, jsq, csq, ssq)
    integer(kind=int32), intent(in) :: n, lbw_ubt, ubw_ubt, lbwmax_ubt, ubwmax_ubt, &
         lbwmax_bv, ubwmax_bv
    real(kind=dp), dimension(lbwmax_ubt+ubwmax_ubt+1,n), intent(inout) :: b_ubt
    integer(kind=int32), dimension(n), intent(in) :: numrotsu, numrotst
    integer(kind=int32), dimension(ubwmax_ubt,n), intent(in) :: jsu
    real(kind=dp), dimension(ubwmax_ubt,n), intent(in) :: csu, ssu
    integer(kind=int32), dimension(n,lbwmax_ubt), intent(in) :: kst
    real(kind=dp), dimension(n,lbwmax_ubt), intent(in) :: cst, sst

    real(kind=dp), dimension(n,lbwmax_bv+ubwmax_bv+1), intent(out) :: b_bv
    integer(kind=int32), dimension(n), intent(out) :: numrotsv
    integer(kind=int32), dimension(n,ubwmax_bv), intent(out) :: ksv
    real(kind=dp), dimension(n,ubwmax_bv), intent(out) :: csv, ssv
    integer(kind=int32), intent(out) :: lbw_bv, ubw_bv
    
    integer(kind=int32), intent(in) :: minind, maxind, maxord
    real(kind=dp), dimension(maxord,minind:maxind) :: csq, ssq
    integer(kind=int32), dimension(minind:maxind), intent(out) :: numrotsq
    integer(kind=int32), dimension(maxord,minind:maxind), intent(out) :: jsq
    integer(kind=int32), intent(out) :: leftq, rightq, incq
    integer(kind=int32) :: j, k, j0, j1, k0, k1, lbw, ubw, lbw1, ubw1
    logical :: full_ubw, full_lbw
    type(d_rotation) :: rot

    if (n == 1) then
       b_bv(1,1)=b_ubt(1,1);
       lbw_bv=0; ubw_bv=0; numrotsv=0; 
       numrotsq=0; leftq=0; rightq=-1; incq=1
       return
    end if

    ubw=min(n-1,ubw_ubt+lbw_ubt)
    b_ubt=eoshift(b_ubt,ubw_ubt-ubw,dim=1_int64)
    if (ubw < n-1) then
       ubw1=ubw+1
       full_ubw=.false.
       b_ubt=eoshift(b_ubt,-1,dim=1_int64)
    else
       ubw1=ubw
       full_ubw=.true.
    end if

    lbw=lbw_ubt
    if (lbw < n-1) then
       lbw1=lbw+1
       full_lbw=.false.
       b_ubt(lbw1+ubw1+1,:) = 0.0_dp
    else
       lbw1=lbw
       full_lbw=.true.
    end if

    numrotsv=0
    csv=0.0_dp; ssv=0.0_dp; ksv=0
    numrotsq=0
    leftq=n-2; rightq=1; incq=-1

    csq=0.0_dp; ssq=0.0_dp; jsq=0
    
    do k=n-1,2,-1
       ! Apply U_k
       do j=1,numrotsu(k)
          call f_rotation_times_tbc(csu(j,k),ssu(j,k),b_ubt,n,lbw1,ubw1,k,0,jsu(j,k))
       end do
       ! Apply T_k
       do j=1,numrotst(k)
          call f_tbc_times_rotation(b_ubt,n,lbw1,ubw1,k,0,cst(k,j),-sst(k,j),kst(k,j))
       end do
       ! rows in which nonzeros have been introduced into the extra subdiagonal
       j0=max(k+1,lbw_ubt+2)
       j1=min(k+lbw_ubt,n)
       ! Apply Q_{k-1}
       numrotsq(k-1)=max(j1-j0+1,0)
       do j=j0,j1
          rot=lgivens(get_el_bc(b_ubt,ubw1,j-1,j-lbw1),get_el_bc(b_ubt,ubw1,j,j-lbw1))
          jsq(j1-j+1,k-1)=j-1
          csq(j1-j+1,k-1)=rot%cosine; ssq(j1-j+1,k-1)=rot%sine
          call rotation_times_tbc(trp_rot(rot),b_ubt,n,lbw1,ubw1,0,0,j-1)
       end do
       ! rows in which nonzeros have been introduced into the extra superdiagonal
       k0=max(k+lbw_ubt+1,lbw_ubt+ubw_ubt+2)
       k1=min(k+2*lbw_ubt+ubw_ubt,n)
       ! Apply V_{k+lbw_ubt-1}
       if (k1 >= k0) then
          numrotsv(k+lbw_ubt-1)=max(k1-k0+1,0)
          do j=k0,k1
            rot=rgivens(get_el_bc(b_ubt,ubw1,j-ubw1,j-1), get_el_bc(b_ubt,ubw1,j-ubw1,j))
            ksv(k+lbw_ubt-1,k1-j+1)=j-1
            csv(k+lbw_ubt-1,k1-j+1)=rot%cosine; ssv(k+lbw_ubt-1,k1-j+1)=rot%sine
            call tbc_times_rotation(b_ubt,n,lbw1,ubw1,0,n-(k+lbw_ubt-1),rot,j-1)
          end do
       end if
    end do
    if (.not. full_ubw) then
       b_ubt=eoshift(b_ubt,1,dim=1_int64)
    end if
    call bc_to_br(b_ubt,b_bv,lbw,ubw)
    lbw_bv=lbw; ubw_bv=ubw
  end subroutine f_d_row_compress

  function z_rc_of(ubt,error) result(swbv)
    type(z_rc) :: swbv
    type(z_ubt), intent(in) :: ubt
    type(error_info), intent(inout), optional :: error
    type(routine_info), parameter :: info=info_z_rc_of
    integer(kind=int32) :: n, lbwmax, ubwmax, lbw, ubw
    type(z_ubt), allocatable :: ubt1

    if (failure(error)) return
    call push_id(info,error)

    n=get_n(ubt)
    lbwmax=get_lbwmax(ubt);
    ubwmax=get_ubwmax(ubt);
    lbw=ubt%lbw
    ubw=ubt%ubw
    swbv%bv=z_new_bv(n,lbw,min(lbw+ubw,n-1))
    ubt1=z_new_ubt(n,min(lbw+1,n-1),min(lbw+ubw+1,n-1))
    call copy(ubt1,ubt)
    swbv%sw=z_new_sweeps(n, 1, n-2, lbw)
    call z_row_compress(ubt1, swbv%bv, swbv%sw, error)
    deallocate(ubt1)
    call pop_id(error)
  end function z_rc_of

  ! Errors
  ! 0: no error
  ! 1: ub%n /= bv%n .or. sw%n /= ub%n
  ! 2: Not enough storage for each sweep.
  ! 3: Not enough storage for the number of sweeps.
  ! 4: Not enough storage in ubt.
  ! 5: Not enough storage in bv.
  subroutine z_row_compress(ubt,bv,sw,error)
    type(z_ubt) :: ubt
    type(z_bv) :: bv
    type(z_sweeps) :: sw
    type(error_info), intent(inout), optional :: error
    type(routine_info), parameter :: info=info_z_row_compress
    integer(kind=int32) :: n, lbw
    lbw=ubt%lbw
    n=get_n(ubt)
    if (failure(error)) return
    call push_id(info, error)
    
    if (n /= get_n(bv) .or. n /= get_n(sw)) then
       call set_error(1, info, error); return
    end if
    if (get_maxord(sw) < ubt%lbw) then
       call set_error(2, info, error); return
    end if
    if (get_maxind(sw) < n-2 .or. get_minind(sw) > 1) then
       call set_error(3, info, error); return
    end if
    if (get_ubwmax(ubt) < min(ubt%lbw+ubt%ubw+1,n-1) .or. &
         get_lbwmax(ubt) < min(ubt%lbw+1,n-1)) then
       call set_error(4, info, error); return
    end if
    if (get_lbwmax(bv) < ubt%lbw .or. &
         get_ubwmax(bv) < min(ubt%lbw+ubt%ubw,n-1)) then
       call set_error(5, info, error); return
    end if
    call f_z_row_compress(ubt%bc, get_n(ubt), ubt%lbw, ubt%ubw, get_lbwmax(ubt), &
         get_ubwmax(ubt), ubt%numrotsu, ubt%jsu, ubt%csu, ubt%ssu, & 
         ubt%numrotst, ubt%kst, ubt%cst, ubt%sst, & 
         bv%br, bv%lbw, bv%ubw, get_lbwmax(bv), get_ubwmax(bv), bv%numrotsv, bv%ksv, &
         bv%csv, bv%ssv, sw%left, sw%right, sw%inc, get_minind(sw), get_maxind(sw), &
         get_maxord(sw), sw%numrots, sw%js, sw%cs, sw%ss)
    call pop_id(error)
  end subroutine z_row_compress

  subroutine f_z_row_compress(b_ubt, n, lbw_ubt, ubw_ubt, lbwmax_ubt, ubwmax_ubt, numrotsu, &
       jsu, csu, ssu, numrotst, kst, cst, sst, &
       b_bv, lbw_bv, ubw_bv, lbwmax_bv, ubwmax_bv, numrotsv, ksv, csv, ssv, &
       leftq, rightq, incq, minind, maxind, maxord, numrotsq, jsq, csq, ssq)
    integer(kind=int32), intent(in) :: n, lbw_ubt, ubw_ubt, lbwmax_ubt, ubwmax_ubt, &
         lbwmax_bv, ubwmax_bv
    complex(kind=dp), dimension(lbwmax_ubt+ubwmax_ubt+1,n), intent(inout) :: b_ubt
    integer(kind=int32), dimension(n), intent(in) :: numrotsu, numrotst
    integer(kind=int32), dimension(ubwmax_ubt,n), intent(in) :: jsu
    real(kind=dp), dimension(ubwmax_ubt,n), intent(in) :: csu
    complex(kind=dp), dimension(ubwmax_ubt,n), intent(in) :: ssu
    integer(kind=int32), dimension(n,lbwmax_ubt), intent(in) :: kst
    real(kind=dp), dimension(n,lbwmax_ubt), intent(in) :: cst
    complex(kind=dp), dimension(n,lbwmax_ubt), intent(in) :: sst

    complex(kind=dp), dimension(n,lbwmax_bv+ubwmax_bv+1), intent(out) :: b_bv
    integer(kind=int32), dimension(n), intent(out) :: numrotsv
    integer(kind=int32), dimension(n,ubwmax_bv), intent(out) :: ksv
    real(kind=dp), dimension(n,ubwmax_bv), intent(out) :: csv
    complex(kind=dp), dimension(n,ubwmax_bv), intent(out) :: ssv
    integer(kind=int32), intent(out) :: lbw_bv, ubw_bv
    
    integer(kind=int32), intent(in) :: minind, maxind, maxord
    real(kind=dp), dimension(maxord,minind:maxind) :: csq
    complex(kind=dp), dimension(maxord,minind:maxind) :: ssq
    integer(kind=int32), dimension(minind:maxind), intent(out) :: numrotsq
    integer(kind=int32), dimension(maxord,minind:maxind), intent(out) :: jsq
    integer(kind=int32), intent(out) :: leftq, rightq, incq
    integer(kind=int32) :: j, k, j0, j1, k0, k1, lbw, ubw, lbw1, ubw1
    logical :: full_ubw, full_lbw
    type(z_rotation) :: rot

    if (n == 1) then
       b_bv(1,1)=b_ubt(1,1);
       lbw_bv=0; ubw_bv=0; numrotsv=0; 
       numrotsq=0; leftq=0; rightq=-1; incq=1;
       return
    end if

    ubw=min(n-1,ubw_ubt+lbw_ubt)
    b_ubt=eoshift(b_ubt,ubw_ubt-ubw,dim=1_int64)
    if (ubw < n-1) then
       ubw1=ubw+1
       full_ubw=.false.
       b_ubt=eoshift(b_ubt,-1,dim=1_int64)
    else
       ubw1=ubw
       full_ubw=.true.
    end if

    lbw=lbw_ubt
    if (lbw < n-1) then
       lbw1=lbw+1
       full_lbw=.false.
       b_ubt(lbw1+ubw1+1,:) = (0.0_dp, 0.0_dp)
    else
       lbw1=lbw
       full_lbw=.true.
    end if

    numrotsv=0
    csv=0.0_dp; ssv=(0.0_dp, 0.0_dp); ksv=0
    numrotsq=0
    leftq=n-2; rightq=1; incq=-1
    csq=0.0_dp; ssq=(0.0_dp, 0.0_dp); jsq=0
    
    do k=n-1,2,-1
       ! Apply U_k
       do j=1,numrotsu(k)
          call f_rotation_times_tbc(csu(j,k),ssu(j,k),b_ubt,n,lbw1,ubw1,k,0,jsu(j,k))
       end do
       ! Apply T_k
       do j=1,numrotst(k)
          call f_tbc_times_rotation(b_ubt,n,lbw1,ubw1,k,0,cst(k,j),-sst(k,j),kst(k,j))
       end do
       ! rows in which nonzeros have been introduced into the extra subdiagonal
       j0=max(k+1,lbw_ubt+2)
       j1=min(k+lbw_ubt,n)
       ! Apply Q_{k-1}
       numrotsq(k-1)=max(j1-j0+1,0)
       do j=j0,j1
          rot=lgivens(get_el_bc(b_ubt,ubw1,j-1,j-lbw1),get_el_bc(b_ubt,ubw1,j,j-lbw1))
          jsq(j1-j+1,k-1)=j-1
          csq(j1-j+1,k-1)=rot%cosine; ssq(j1-j+1,k-1)=rot%sine
          call rotation_times_tbc(trp_rot(rot),b_ubt,n,lbw1,ubw1,0,0,j-1)
       end do
       ! rows in which nonzeros have been introduced into the extra superdiagonal
       k0=max(k+lbw_ubt+1,lbw_ubt+ubw_ubt+2)
       k1=min(k+2*lbw_ubt+ubw_ubt,n)
       ! Apply V_{k+lbw_ubt-1}
       if (k1 >= k0) then
          numrotsv(k+lbw_ubt-1)=max(k1-k0+1,0)
          do j=k0,k1
            rot=rgivens(get_el_bc(b_ubt,ubw1,j-ubw1,j-1), get_el_bc(b_ubt,ubw1,j-ubw1,j))
            ksv(k+lbw_ubt-1,k1-j+1)=j-1
            csv(k+lbw_ubt-1,k1-j+1)=rot%cosine; ssv(k+lbw_ubt-1,k1-j+1)=rot%sine
            call tbc_times_rotation(b_ubt,n,lbw1,ubw1,0,n-(k+lbw_ubt-1),rot,j-1)
          end do
       end if
    end do

    if (.not. full_ubw) then
       b_ubt=eoshift(b_ubt,1,dim=1_int64)
    end if
    call bc_to_br(b_ubt,b_bv,lbw,ubw)
    lbw_bv=lbw; ubw_bv=ubw
  end subroutine f_z_row_compress

end module mod_row_compress
