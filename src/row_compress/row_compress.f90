module mod_row_compress
  use mod_misc
  use mod_sweeps
  use mod_rotation
  use mod_types
  use mod_shift
  implicit none

  interface row_compress
     module procedure d_row_compress, c_row_compress
  end interface row_compress

  interface f_row_compress
     module procedure f_d_row_compress, f_c_row_compress
  end interface f_row_compress

contains

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
    type(error_info), intent(out) :: error

    integer(kind=int32) :: n, lbw
    lbw=ubt%lbw
    n=get_n(ubt)
    call clear_error(error)
    if (n /= get_n(bv) .or. n /= get_n(sw)) then
       call set_error(error, 1, id_d_row_compress); return
    end if
    if (get_maxord(sw) < ubt%lbw) then
       call set_error(error, 2, id_d_row_compress); return
    end if
    if (get_maxind(sw) < n-2 .or. get_minind(sw) > 1) then
       call set_error(error, 3, id_d_row_compress); return
    end if
    if (get_ubwmax(ubt) < min(ubt%lbw+ubt%ubw+1,n-1) .or. &
         get_lbwmax(ubt) < min(ubt%lbw+1,n-1)) then
       call set_error(error, 4, id_d_row_compress); return
    end if
    if (get_lbwmax(bv) < ubt%lbw .or. &
         get_ubwmax(bv) < min(ubt%lbw+ubt%ubw,n-1)) then
       call set_error(error, 5, id_d_row_compress); return
    end if
    call f_d_row_compress(ubt%bc, get_n(ubt), ubt%lbw, ubt%ubw, get_lbwmax(ubt), &
         get_ubwmax(ubt), ubt%numrotsu, ubt%jsu, ubt%csu, ubt%ssu, & 
         ubt%numrotst, ubt%kst, ubt%cst, ubt%sst, & 
         bv%br, bv%lbw, bv%ubw, get_lbwmax(bv), get_ubwmax(bv), bv%numrotsv, bv%ksv, &
         bv%csv, bv%ssv, sw%left, sw%right, sw%inc, get_minind(sw), get_maxind(sw), &
         get_maxord(sw), sw%numrots, sw%js, sw%cs, sw%ss)
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
    call shift2(b_ubt,ubw-ubw_ubt,0)
    if (ubw < n-1) then
       ubw1=ubw+1
       full_ubw=.false.
       call shift2(b_ubt,1,0)
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
          rot%cosine=csu(j,k); rot%sine=ssu(j,k)
          call rotation_times_tbc(rot,b_ubt,n,lbw1,ubw1,k,0,jsu(j,k))
       end do
       ! Apply T_k
       do j=1,numrotst(k)
          rot%cosine=cst(k,j); rot%sine=sst(k,j)
          call tbc_times_rotation(b_ubt,n,lbw1,ubw1,k,0,trp_rot(rot),kst(k,j))
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
       call shift2(b_ubt,-1,0)
    end if
    call bc_to_br(b_ubt,b_bv,lbw,ubw)
    lbw_bv=lbw; ubw_bv=ubw
  end subroutine f_d_row_compress

  ! Errors
  ! 0: no error
  ! 1: ub%n /= bv%n .or. sw%n /= ub%n
  ! 2: Not enough storage for each sweep.
  ! 3: Not enough storage for the number of sweeps.
  ! 4: Not enough storage in ubt.
  ! 5: Not enough storage in bv.
  subroutine c_row_compress(ubt,bv,sw,error)
    type(c_ubt) :: ubt
    type(c_bv) :: bv
    type(c_sweeps) :: sw
    type(error_info), intent(out) :: error

    integer(kind=int32) :: n, lbw
    lbw=ubt%lbw
    n=get_n(ubt)
    call clear_error(error)
    if (n /= get_n(bv) .or. n /= get_n(sw)) then
       call set_error(error, 1, id_c_row_compress); return
    end if
    if (get_maxord(sw) < ubt%lbw) then
       call set_error(error, 2, id_c_row_compress); return
    end if
    if (get_maxind(sw) < n-2 .or. get_minind(sw) > 1) then
       call set_error(error, 3, id_c_row_compress); return
    end if
    if (get_ubwmax(ubt) < min(ubt%lbw+ubt%ubw+1,n-1) .or. &
         get_lbwmax(ubt) < min(ubt%lbw+1,n-1)) then
       call set_error(error, 4, id_c_row_compress); return
    end if
    if (get_lbwmax(bv) < ubt%lbw .or. &
         get_ubwmax(bv) < min(ubt%lbw+ubt%ubw,n-1)) then
       call set_error(error, 5, id_d_row_compress); return
    end if
    call f_c_row_compress(ubt%bc, get_n(ubt), ubt%lbw, ubt%ubw, get_lbwmax(ubt), &
         get_ubwmax(ubt), ubt%numrotsu, ubt%jsu, ubt%csu, ubt%ssu, & 
         ubt%numrotst, ubt%kst, ubt%cst, ubt%sst, & 
         bv%br, bv%lbw, bv%ubw, get_lbwmax(bv), get_ubwmax(bv), bv%numrotsv, bv%ksv, &
         bv%csv, bv%ssv, sw%left, sw%right, sw%inc, get_minind(sw), get_maxind(sw), &
         get_maxord(sw), sw%numrots, sw%js, sw%cs, sw%ss)
  end subroutine c_row_compress

  subroutine f_c_row_compress(b_ubt, n, lbw_ubt, ubw_ubt, lbwmax_ubt, ubwmax_ubt, numrotsu, &
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
    type(c_rotation) :: rot

    if (n == 1) then
       b_bv(1,1)=b_ubt(1,1);
       lbw_bv=0; ubw_bv=0; numrotsv=0; 
       numrotsq=0; leftq=0; rightq=-1; incq=1;
       return
    end if

    ubw=min(n-1,ubw_ubt+lbw_ubt)
    call shift2(b_ubt,ubw-ubw_ubt,0)
    if (ubw < n-1) then
       ubw1=ubw+1
       full_ubw=.false.
       call shift2(b_ubt,1,0)
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
          rot%cosine=csu(j,k); rot%sine=ssu(j,k)
          call rotation_times_tbc(rot,b_ubt,n,lbw1,ubw1,k,0,jsu(j,k))
       end do
       ! Apply T_k
       do j=1,numrotst(k)
          rot%cosine=cst(k,j); rot%sine=sst(k,j)
          call tbc_times_rotation(b_ubt,n,lbw1,ubw1,k,0,trp_rot(rot),kst(k,j))
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
       call shift2(b_ubt,-1,0)
    end if
    call bc_to_br(b_ubt,b_bv,lbw,ubw)
    lbw_bv=lbw; ubw_bv=ubw
  end subroutine f_c_row_compress

end module mod_row_compress

