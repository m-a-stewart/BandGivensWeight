module mod_update
  use mod_prec
  use mod_error_id
  use mod_nested_types
  use mod_band_types
  use mod_shift
  use mod_rotation
  use mod_sweeps1
  use mod_convert_bv_to_ub
  use mod_convert_ub_to_bv
  implicit none

  private

  public :: r1_update_ub_to_bv, d_r1_update_ub_to_bv, c_r1_update_ub_to_bv, &
       f_r1_update_ub_to_bv, f_d_r1_update_ub_to_bv, f_c_r1_update_ub_to_bv

  public :: e1v_update_ub_to_bv, d_e1v_update_ub_to_bv, c_e1v_update_ub_to_bv, &
       f_e1v_update_ub_to_bv, f_d_e1v_update_ub_to_bv, f_c_e1v_update_ub_to_bv       

  public :: info_d_r1_update_ub_to_bv, info_c_r1_update_ub_to_bv, info_d_e1v_update_ub_to_bv, &
       info_c_e1v_update_ub_to_bv

  interface r1_update_ub_to_bv
     module procedure d_r1_update_ub_to_bv, c_r1_update_ub_to_bv
  end interface r1_update_ub_to_bv

  interface f_r1_update_ub_to_bv
     module procedure f_d_r1_update_ub_to_bv, f_c_r1_update_ub_to_bv
  end interface f_r1_update_ub_to_bv

  interface e1v_update_ub_to_bv
     module procedure d_e1v_update_ub_to_bv, c_e1v_update_ub_to_bv
  end interface e1v_update_ub_to_bv

  interface f_e1v_update_ub_to_bv
     module procedure f_d_e1v_update_ub_to_bv, f_c_e1v_update_ub_to_bv
  end interface f_e1v_update_ub_to_bv

  type(routine_info), parameter :: info_d_r1_update_ub_to_bv=routine_info(id_d_r1_update_ub_to_bv, &
       'd_r1_update_ub_to_bv', &
       [ character(len=error_message_length) :: 'n<1', 'Size mismatch.', 'Insufficient storage for a sweep', &
       'Insufficient lbw.', 'Insufficient ubw' ] )

  type(routine_info), parameter :: info_c_r1_update_ub_to_bv=routine_info(id_c_r1_update_ub_to_bv, &
       'c_r1_update_ub_to_bv', &
       [ character(len=error_message_length) :: 'n<1', 'Size mismatch.', 'Insufficient storage for a sweep', &
       'Insufficient lbw.', 'Insufficient ubw' ] )

  type(routine_info), parameter :: info_d_e1v_update_ub_to_bv=routine_info(id_d_e1v_update_ub_to_bv, &
       'd_e1v_update_ub_to_bv', &
       [ character(len=error_message_length) :: 'n<1', 'Size mismatch.', &
       'Insufficient lbw.', 'Insufficient ubw' ] )

  type(routine_info), parameter :: info_c_e1v_update_ub_to_bv=routine_info(id_c_e1v_update_ub_to_bv, &
       'c_e1v_update_ub_to_bv', &
       [ character(len=error_message_length) :: 'n<1', 'Size mismatch.', &
       'Insufficient lbw.', 'Insufficient ubw' ] )

contains

  !
  ! Compute a factorization Q ( U . B_0 + u v^T) = (B_1 + b e_1 (c
  ! e_1^T + d e_2^T) . V Once this transformation is done, a rank one
  ! update is completed simply by adding bc to B_1(1,1) and bd to B_1(1,2).
  !
  ! Errors:
  ! 1: n < 1
  ! 2: size mismatch.
  ! 3: insufficient storage for a sweep.
  ! 4: insufficient lbw in bv or ub.
  ! 5: insufficient ubw in bv or ub.


  subroutine d_r1_update_ub_to_bv(ub,u,v,sw,bv,error)
    type(d_bv) :: bv
    type(d_ub) :: ub
    type(d_sweeps1) :: sw
    real(kind=dp), dimension(:) :: u, v
    type(error_info) :: error

    call clear_error(error)
    
    if(get_n(bv) < 1) then
       call set_error(error,1,id_d_r1_update_ub_to_bv)
    end if
    if(get_n(bv) /= get_n(ub) .or. get_n(bv) /= size(u) &
         .or. get_n(bv) /= size(v)) then
       call set_error(error,2,id_d_r1_update_ub_to_bv)
    end if
    if (get_maxsweeps1(sw) - sw%numsweeps1 < 1) then
       call set_error(error,3,id_d_r1_update_ub_to_bv)
    end if
    if (get_lbwmax(ub) < ub%lbw+1 .or. &
         get_lbwmax(bv) < ub%lbw+1) then
       call set_error(error,4,id_d_r1_update_ub_to_bv)
    end if
    if (get_ubwmax(ub) < ub%ubw+3 .or. &
         get_ubwmax(bv) < ub%ubw+3) then
       call set_error(error,5,id_d_r1_update_ub_to_bv)
    end if
    call f_d_r1_update_ub_to_bv(ub%bc, get_n(ub), ub%lbw, ub%ubw, get_lbwmax(ub), get_ubwmax(ub), &
         ub%numrotsu, ub%jsu, ub%csu, ub%ssu, u, v, &
         bv%br, bv%lbw, bv%ubw, get_lbwmax(bv), get_ubwmax(bv), &
         bv%numrotsv, bv%ksv, bv%csv, bv%ssv, &
         sw%cs, sw%ss, sw%numsweeps1, get_maxsweeps1(sw), &
         error)
  end subroutine d_r1_update_ub_to_bv

  subroutine f_d_r1_update_ub_to_bv(b_ub, n, lbw_ub, ubw_ub, lbwmax_ub, ubwmax_ub, numrotsu, &
       jsu, csu, ssu, u, v, &
       b_bv, lbw_bv, ubw_bv, lbwmax_bv, ubwmax_bv, numrotsv, ksv, csv, ssv, &
       cs_sw, ss_sw, numsweeps1_sw, maxsweeps1_sw, &
       error)

    real(kind=dp), dimension(lbwmax_ub+ubwmax_ub+1,n), intent(inout) :: b_ub
    integer(kind=int32), dimension(n), intent(out) :: numrotsu
    integer(kind=int32), dimension(ubwmax_ub,n), intent(inout) :: jsu
    real(kind=dp), dimension(ubwmax_ub,n), intent(inout) :: csu, ssu

    real(kind=dp), dimension(n,maxsweeps1_sw), intent(out) :: cs_sw, ss_sw
    integer(kind=int32), intent(inout) :: numsweeps1_sw
    integer(kind=int32), intent(in) :: maxsweeps1_sw, n

    real(kind=dp), dimension(n, lbwmax_bv+ubwmax_bv+1), intent(out) :: b_bv
    integer(kind=int32), dimension(n), intent(out) :: numrotsv
    integer(kind=int32), dimension(n,ubwmax_bv), intent(out) :: ksv
    real(kind=dp), dimension(n,ubwmax_bv), intent(out) :: csv, ssv

    integer(kind=int32), intent(in) :: lbwmax_ub, ubwmax_ub
    integer(kind=int32), intent(inout) :: lbw_ub, ubw_ub
    integer(kind=int32), intent(in) :: lbwmax_bv, ubwmax_bv
    integer(kind=int32), intent(out) :: lbw_bv, ubw_bv

    real(kind=dp), dimension(:), intent(inout) :: u, v 

    type(error_info), intent(out) :: error
    integer(kind=int32) :: j
    type(d_rotation) :: rot

    call clear_error(error)
    if (n==1) then
       b_bv(1,1)=b_ub(1,1);
       lbw_bv=0; ubw_bv=0; numrotsv=0; return
    end if

    b_bv(:,1:lbw_ub+ubw_ub+4)=0.0_dp; numrotsv=0
    ssv(:,1:ubw_ub+3)=0.0_dp; csv(:,1:ubw_ub+3)=0.0_dp
    ksv(:,1:ubw_ub+3)=0

    ! zero elements in u
    if (numsweeps1_sw >= 1) then
       call shift2(cs_sw(:,1:numsweeps1_sw+1),0,1)
       call shift2(ss_sw(:,1:numsweeps1_sw+1),0,1)
    end if
    numsweeps1_sw=numsweeps1_sw+1
    do j=n,2,-1
       rot = lgivens(u(j-1), u(j))
       call rotation_times_general(trp_rot(rot), u, j-1, j)
       u(j)=0.0_dp
       cs_sw(j-1,1)=rot%cosine
       ss_sw(j-1,1)=-rot%sine
    end do

    ! apply to ub, filling in one subdiagonal and increasing ubw by one.
    call f_d_sweeps1_times_ub(cs_sw, ss_sw, 1, maxsweeps1_sw, n, &
         b_ub, lbw_ub, ubw_ub, lbwmax_ub, ubwmax_ub, numrotsu, jsu, csu, ssu, &
         b_bv, lbw_bv, ubw_bv, lbwmax_bv, ubwmax_bv, numrotsv, ksv, csv, ssv, error)

    call f_d_convert_bv_to_ub(b_bv, n, lbw_bv, ubw_bv, lbwmax_bv, ubwmax_bv, &
         numrotsv, ksv, csv, ssv, b_ub, lbw_ub, ubw_ub, lbwmax_ub, ubwmax_ub, numrotsu, jsu, &
         csu, ssu, error)

    call f_d_e1v_update_ub_to_bv(b_ub, n, lbw_ub, ubw_ub, lbwmax_ub, ubwmax_ub, numrotsu, &
         jsu, csu, ssu, v, &
         b_bv, lbw_bv, ubw_bv, lbwmax_bv, ubwmax_bv, numrotsv, ksv, csv, ssv, &
         error)
  end subroutine f_d_r1_update_ub_to_bv

  subroutine d_e1v_update_ub_to_bv(ub,v,bv,error)
    type(d_bv) :: bv
    type(d_ub) :: ub
    real(kind=dp), dimension(:) :: v
    type(error_info) :: error

    call clear_error(error)
    
    if (get_n(bv) < 1) then
       call set_error(error,1,id_d_e1v_update_ub_to_bv)
    end if
    if (get_n(bv) /= get_n(ub) .or. get_n(bv) /= size(v)) then
       call set_error(error,2,id_d_e1v_update_ub_to_bv)
    end if
    if (get_lbwmax(bv) < ub%lbw) then
       call set_error(error,3,id_d_e1v_update_ub_to_bv)
    end if
    if (get_ubwmax(ub) < ub%ubw+2 .or. &
         get_ubwmax(bv) < ub%ubw+2) then
       call set_error(error,4,id_d_e1v_update_ub_to_bv)
    end if
    call f_d_e1v_update_ub_to_bv(ub%bc, get_n(ub), ub%lbw, ub%ubw, get_lbwmax(ub), get_ubwmax(ub), &
         ub%numrotsu, ub%jsu, ub%csu, ub%ssu, v, &
         bv%br, bv%lbw, bv%ubw, get_lbwmax(bv), get_ubwmax(bv), &
         bv%numrotsv, bv%ksv, bv%csv, bv%ssv, &
         error)
  end subroutine d_e1v_update_ub_to_bv

  subroutine f_d_e1v_update_ub_to_bv(b_ub, n, lbw_ub, ubw_ub, lbwmax_ub, ubwmax_ub, numrotsu, &
       jsu, csu, ssu, v, &
       b_bv, lbw_bv, ubw_bv, lbwmax_bv, ubwmax_bv, numrotsv, ksv, csv, ssv, &
       error)

    real(kind=dp), dimension(lbwmax_ub+ubwmax_ub+1,n), intent(inout) :: b_ub
    integer(kind=int32), dimension(n), intent(out) :: numrotsu
    integer(kind=int32), dimension(ubwmax_ub,n), intent(inout) :: jsu
    real(kind=dp), dimension(ubwmax_ub,n), intent(inout) :: csu, ssu

    integer(kind=int32), intent(in) :: n

    real(kind=dp), dimension(n, lbwmax_bv+ubwmax_bv+1), intent(out) :: b_bv
    integer(kind=int32), dimension(n), intent(out) :: numrotsv
    integer(kind=int32), dimension(n,ubwmax_bv), intent(out) :: ksv
    real(kind=dp), dimension(n,ubwmax_bv), intent(out) :: csv, ssv

    integer(kind=int32), intent(in) :: lbwmax_ub, ubwmax_ub
    integer(kind=int32), intent(inout) :: lbw_ub, ubw_ub
    integer(kind=int32), intent(in) :: lbwmax_bv, ubwmax_bv
    integer(kind=int32), intent(out) :: lbw_bv, ubw_bv

    real(kind=dp), dimension(:), intent(inout) :: v

    type(error_info), intent(out) :: error
    integer(kind=int32) :: j, k, k0, k1, dubw, dubw_tmp, dlbw, dlbw_tmp
    type(d_rotation) :: rot

    call clear_error(error)
    if (n==1) then
       b_bv(1,1)=b_ub(1,1);
       lbw_bv=0; ubw_bv=0; numrotsv=0; return
    end if
    
    dubw=1; dubw_tmp=1
    dlbw=0; dlbw_tmp=0

    call f_bw_expand_bc(b_ub, n, lbw_ub, ubw_ub, lbwmax_ub, ubwmax_ub, &
         dlbw, dlbw_tmp, dubw, dubw_tmp, lbw_bv, ubw_bv)

    b_bv(:,1:lbw_bv+ubw_bv+1)=0.0_dp
    numrotsv=0
    ssv(:,1:ubw_bv)=0.0_dp; csv(:,1:ubw_bv)=0.0_dp
    ksv(:,1:ubw_bv)=0

    ! Require d=min(ubw_ub+2,n-1) - ubw_ub extra diagonals
    ! Now compute a bv decomposition in which transformations are applied to v^T.
    do k=1,n-2
       ! Apply u_{n-k}
       do j=1,numrotsu(n-k)
          rot%cosine=csu(j,n-k); rot%sine=ssu(j,n-k)
          call rotation_times_tbc(rot,b_ub,n,lbw_ub,ubw_ub,n-k,0,jsu(j,n-k))
       end do
       k0=max(n-k+2,ubw_ub+1)
       k1=min(n-k+ubw_ub+1,n)
       numrotsv(n-(k+1))=k1-k0+2
       rot=rgivens(v(n-k),v(n-k+1))
       call general_times_rotation(v,rot,n-k,n-k+1)
       call tbc_times_rotation(b_ub,n,lbw_ub,ubw_ub,0,k+1,rot,n-k)
       ksv(n-(k+1),k1-k0+2)=n-k
       csv(n-(k+1),k1-k0+2)=rot%cosine
       ssv(n-(k+1),k1-k0+2)=rot%sine
       do j=k0,k1 
          rot=rgivens(get_el_bc(b_ub,ubw_ub,j-ubw_ub,j-1), get_el_bc(b_ub,ubw_ub,j-ubw_ub,j))
          ksv(n-(k+1),k1-j+1)=j-1
          csv(n-(k+1),k1-j+1)=rot%cosine; ssv(n-(k+1),k1-j+1)=rot%sine
          call tbc_times_rotation(b_ub,n,lbw_ub,ubw_ub,0,k+1,rot,j-1)
       end do
    end do
    call f_bw_contract_bc(b_ub, n, lbw_ub, ubw_ub, lbwmax_ub, ubwmax_ub, &
         dlbw, dlbw_tmp, dubw, dubw_tmp)
    call bc_to_br(b_ub,b_bv,lbw_bv,ubw_bv)
  end subroutine f_d_e1v_update_ub_to_bv



  subroutine c_r1_update_ub_to_bv(ub,u,v,sw,bv,error)
    type(c_bv) :: bv
    type(c_ub) :: ub
    type(c_sweeps1) :: sw
    complex(kind=dp), dimension(:) :: u, v
    type(error_info) :: error

    call clear_error(error)
    
    if(get_n(bv) < 1) then
       call set_error(error,1,id_c_r1_update_ub_to_bv)
    end if
    if(get_n(bv) /= get_n(ub) .or. get_n(bv) /= size(u) &
         .or. get_n(bv) /= size(v)) then
       call set_error(error,2,id_c_r1_update_ub_to_bv)
    end if
    if (get_maxsweeps1(sw) - sw%numsweeps1 < 1) then
       call set_error(error,3,id_c_r1_update_ub_to_bv)
    end if
    if (get_lbwmax(ub) < ub%lbw+1 .or. &
         get_lbwmax(bv) < ub%lbw+1) then
       call set_error(error,4,id_c_r1_update_ub_to_bv)
    end if
    if (get_ubwmax(ub) < ub%ubw+3 .or. &
         get_ubwmax(bv) < ub%ubw+3) then
       call set_error(error,5,id_c_r1_update_ub_to_bv)
    end if
    call f_c_r1_update_ub_to_bv(ub%bc, get_n(ub), ub%lbw, ub%ubw, get_lbwmax(ub), get_ubwmax(ub), &
         ub%numrotsu, ub%jsu, ub%csu, ub%ssu, u, v, &
         bv%br, bv%lbw, bv%ubw, get_lbwmax(bv), get_ubwmax(bv), &
         bv%numrotsv, bv%ksv, bv%csv, bv%ssv, &
         sw%cs, sw%ss, sw%numsweeps1, get_maxsweeps1(sw), &
         error)
  end subroutine c_r1_update_ub_to_bv

  subroutine f_c_r1_update_ub_to_bv(b_ub, n, lbw_ub, ubw_ub, lbwmax_ub, ubwmax_ub, numrotsu, &
       jsu, csu, ssu, u, v, &
       b_bv, lbw_bv, ubw_bv, lbwmax_bv, ubwmax_bv, numrotsv, ksv, csv, ssv, &
       cs_sw, ss_sw, numsweeps1_sw, maxsweeps1_sw, &
       error)

    complex(kind=dp), dimension(lbwmax_ub+ubwmax_ub+1,n), intent(inout) :: b_ub
    integer(kind=int32), dimension(n), intent(out) :: numrotsu
    integer(kind=int32), dimension(ubwmax_ub,n), intent(inout) :: jsu
    real(kind=dp), dimension(ubwmax_ub,n), intent(inout) :: csu
    complex(kind=dp), dimension(ubwmax_ub,n), intent(inout) :: ssu

    real(kind=dp), dimension(n,maxsweeps1_sw), intent(out) :: cs_sw
    complex(kind=dp), dimension(n,maxsweeps1_sw), intent(out) :: ss_sw
    integer(kind=int32), intent(inout) :: numsweeps1_sw
    integer(kind=int32), intent(in) :: maxsweeps1_sw, n

    complex(kind=dp), dimension(n, lbwmax_bv+ubwmax_bv+1), intent(out) :: b_bv
    integer(kind=int32), dimension(n), intent(out) :: numrotsv
    integer(kind=int32), dimension(n,ubwmax_bv), intent(out) :: ksv
    real(kind=dp), dimension(n,ubwmax_bv), intent(out) :: csv
    complex(kind=dp), dimension(n,ubwmax_bv), intent(out) :: ssv

    integer(kind=int32), intent(in) :: lbwmax_ub, ubwmax_ub
    integer(kind=int32), intent(inout) :: lbw_ub, ubw_ub
    integer(kind=int32), intent(in) :: lbwmax_bv, ubwmax_bv
    integer(kind=int32), intent(out) :: lbw_bv, ubw_bv

    complex(kind=dp), dimension(:), intent(inout) :: u, v 

    type(error_info), intent(out) :: error
    integer(kind=int32) :: j
    type(c_rotation) :: rot

    call clear_error(error)
    if (n==1) then
       b_bv(1,1)=b_ub(1,1);
       lbw_bv=0; ubw_bv=0; numrotsv=0; return
    end if

    b_bv(:,1:lbw_ub+ubw_ub+4)=(0.0_dp,0.0_dp); numrotsv=0
    ssv(:,1:ubw_ub+3)=(0.0_dp,0.0_dp); csv(:,1:ubw_ub+3)=0.0_dp
    ksv(:,1:ubw_ub+3)=0

    ! zero elements in u
    if (numsweeps1_sw >= 1) then
       call shift2(cs_sw(:,1:numsweeps1_sw+1),0,1)
       call shift2(ss_sw(:,1:numsweeps1_sw+1),0,1)
    end if
    numsweeps1_sw=numsweeps1_sw+1
    do j=n,2,-1
       rot = lgivens(u(j-1), u(j))
       call rotation_times_general(trp_rot(rot), u, j-1, j)
       u(j)=(0.0_dp,0.0_dp)
       cs_sw(j-1,1)=rot%cosine
       ss_sw(j-1,1)=-rot%sine
    end do

    ! apply to ub, filling in one subdiagonal and increasing ubw by one.
    call f_c_sweeps1_times_ub(cs_sw, ss_sw, 1, maxsweeps1_sw, n, &
         b_ub, lbw_ub, ubw_ub, lbwmax_ub, ubwmax_ub, numrotsu, jsu, csu, ssu, &
         b_bv, lbw_bv, ubw_bv, lbwmax_bv, ubwmax_bv, numrotsv, ksv, csv, ssv, error)

    call f_c_convert_bv_to_ub(b_bv, n, lbw_bv, ubw_bv, lbwmax_bv, ubwmax_bv, &
         numrotsv, ksv, csv, ssv, b_ub, lbw_ub, ubw_ub, lbwmax_ub, ubwmax_ub, numrotsu, jsu, &
         csu, ssu, error)

    call f_c_e1v_update_ub_to_bv(b_ub, n, lbw_ub, ubw_ub, lbwmax_ub, ubwmax_ub, numrotsu, &
         jsu, csu, ssu, v, &
         b_bv, lbw_bv, ubw_bv, lbwmax_bv, ubwmax_bv, numrotsv, ksv, csv, ssv, &
         error)
  end subroutine f_c_r1_update_ub_to_bv

  subroutine c_e1v_update_ub_to_bv(ub,v,bv,error)
    type(c_bv) :: bv
    type(c_ub) :: ub
    complex(kind=dp), dimension(:) :: v
    type(error_info) :: error

    call clear_error(error)
    
    if (get_n(bv) < 1) then
       call set_error(error,1,id_c_e1v_update_ub_to_bv)
    end if
    if (get_n(bv) /= get_n(ub) .or. get_n(bv) /= size(v)) then
       call set_error(error,2,id_c_e1v_update_ub_to_bv)
    end if
    if (get_lbwmax(bv) < ub%lbw) then
       call set_error(error,3,id_c_e1v_update_ub_to_bv)
    end if
    if (get_ubwmax(ub) < ub%ubw+2 .or. &
         get_ubwmax(bv) < ub%ubw+2) then
       call set_error(error,4,id_c_e1v_update_ub_to_bv)
    end if
    call f_c_e1v_update_ub_to_bv(ub%bc, get_n(ub), ub%lbw, ub%ubw, get_lbwmax(ub), get_ubwmax(ub), &
         ub%numrotsu, ub%jsu, ub%csu, ub%ssu, v, &
         bv%br, bv%lbw, bv%ubw, get_lbwmax(bv), get_ubwmax(bv), &
         bv%numrotsv, bv%ksv, bv%csv, bv%ssv, &
         error)
  end subroutine c_e1v_update_ub_to_bv

  subroutine f_c_e1v_update_ub_to_bv(b_ub, n, lbw_ub, ubw_ub, lbwmax_ub, ubwmax_ub, numrotsu, &
       jsu, csu, ssu, v, &
       b_bv, lbw_bv, ubw_bv, lbwmax_bv, ubwmax_bv, numrotsv, ksv, csv, ssv, &
       error)

    complex(kind=dp), dimension(lbwmax_ub+ubwmax_ub+1,n), intent(inout) :: b_ub
    integer(kind=int32), dimension(n), intent(out) :: numrotsu
    integer(kind=int32), dimension(ubwmax_ub,n), intent(inout) :: jsu
    real(kind=dp), dimension(ubwmax_ub,n), intent(inout) :: csu
    complex(kind=dp), dimension(ubwmax_ub,n), intent(inout) :: ssu

    integer(kind=int32), intent(in) :: n

    complex(kind=dp), dimension(n, lbwmax_bv+ubwmax_bv+1), intent(out) :: b_bv
    integer(kind=int32), dimension(n), intent(out) :: numrotsv
    integer(kind=int32), dimension(n,ubwmax_bv), intent(out) :: ksv
    real(kind=dp), dimension(n,ubwmax_bv), intent(out) :: csv
    complex(kind=dp), dimension(n,ubwmax_bv), intent(out) :: ssv

    integer(kind=int32), intent(in) :: lbwmax_ub, ubwmax_ub
    integer(kind=int32), intent(inout) :: lbw_ub, ubw_ub
    integer(kind=int32), intent(in) :: lbwmax_bv, ubwmax_bv
    integer(kind=int32), intent(out) :: lbw_bv, ubw_bv

    complex(kind=dp), dimension(:), intent(inout) :: v

    type(error_info), intent(out) :: error
    integer(kind=int32) :: j, k, k0, k1, dubw, dubw_tmp, dlbw, dlbw_tmp
    type(c_rotation) :: rot

    call clear_error(error)
    if (n==1) then
       b_bv(1,1)=b_ub(1,1);
       lbw_bv=0; ubw_bv=0; numrotsv=0; return
    end if
    
    dubw=1; dubw_tmp=1
    dlbw=0; dlbw_tmp=0

    call f_bw_expand_bc(b_ub, n, lbw_ub, ubw_ub, lbwmax_ub, ubwmax_ub, &
         dlbw, dlbw_tmp, dubw, dubw_tmp, lbw_bv, ubw_bv)

    b_bv(:,1:lbw_bv+ubw_bv+1)=(0.0_dp, 0.0_dp)
    numrotsv=0
    ssv(:,1:ubw_bv)=(0.0_dp,0.0_dp); csv(:,1:ubw_bv)=0.0_dp
    ksv(:,1:ubw_bv)=0

    ! Now compute a bv decomposition in which transformations are applied to v^T.
    v=conjg(v)
    do k=1,n-2
       ! Apply u_{n-k}
       do j=1,numrotsu(n-k)
          rot%cosine=csu(j,n-k); rot%sine=ssu(j,n-k)
          call rotation_times_tbc(rot,b_ub,n,lbw_ub,ubw_ub,n-k,0,jsu(j,n-k))
       end do
       k0=max(n-k+2,ubw_ub+1)
       k1=min(n-k+ubw_ub+1,n)
       numrotsv(n-(k+1))=k1-k0+2
       rot=rgivens(v(n-k),v(n-k+1))
       call general_times_rotation(v,rot,n-k,n-k+1)
       call tbc_times_rotation(b_ub,n,lbw_ub,ubw_ub,0,k+1,rot,n-k)
       ksv(n-(k+1),k1-k0+2)=n-k
       csv(n-(k+1),k1-k0+2)=rot%cosine
       ssv(n-(k+1),k1-k0+2)=rot%sine
       do j=k0,k1 
          rot=rgivens(get_el_bc(b_ub,ubw_ub,j-ubw_ub,j-1), get_el_bc(b_ub,ubw_ub,j-ubw_ub,j))
          ksv(n-(k+1),k1-j+1)=j-1
          csv(n-(k+1),k1-j+1)=rot%cosine; ssv(n-(k+1),k1-j+1)=rot%sine
          call tbc_times_rotation(b_ub,n,lbw_ub,ubw_ub,0,k+1,rot,j-1)
       end do
    end do
    call f_bw_contract_bc(b_ub, n, lbw_ub, ubw_ub, lbwmax_ub, ubwmax_ub, &
         dlbw, dlbw_tmp, dubw, dubw_tmp)
    call bc_to_br(b_ub,b_bv,lbw_bv,ubw_bv)
    v=conjg(v)
  end subroutine f_c_e1v_update_ub_to_bv

End module mod_update
