module conversions_ub_to_bv
use prec
use shift
use rotation
use band_types
use nested_types
implicit none

interface convert_ub_to_bv
   module procedure d_convert_ub_to_bv, c_convert_ub_to_bv
end interface convert_ub_to_bv

interface f_convert_ub_to_bv
   module procedure f_d_convert_ub_to_bv, f_c_convert_ub_to_bv
end interface f_convert_ub_to_bv

contains

  ! Errors:
  ! 0: no error
  ! 1: n<1
  ! 3: ub%n /= bv%n

  subroutine d_convert_ub_to_bv(ub, bv, error)
    type(d_ub) :: ub
    type(d_bv) :: bv
    integer(kind=int32), intent(out) :: error
    if (get_n(ub) /= get_n(bv)) then
       error = 3; return
    end if
    call f_d_convert_ub_to_bv(ub%b, get_n(ub), ub%lbw, ub%ubw, get_lbwmax(ub), &
         get_ubwmax(ub), ub%numrotsu, ub%jsu, ub%csu, ub%ssu, bv%b, get_lbwmax(bv), &
         get_ubwmax(bv), bv%numrotsv, bv%ksv, bv%csv, bv%ssv, error)
    bv%lbw=ub%lbw; bv%ubw=ub%ubw
  end subroutine d_convert_ub_to_bv

  subroutine f_d_convert_ub_to_bv(b_ub, n, lbw, ubw, lbwmax_ub, ubwmax_ub, numrots_ub, &
       js_ub, cs_ub, ss_ub, b_bv, lbwmax_bv, ubwmax_bv, numrots_bv, ks_bv, cs_bv, ss_bv, error)
    real(kind=dp), dimension(lbwmax_ub+ubwmax_ub+1,n), intent(inout) :: b_ub
    integer(kind=int32), intent(in) :: n, lbw, ubw, lbwmax_ub, ubwmax_ub, lbwmax_bv, ubwmax_bv
    integer(kind=int32), dimension(n), intent(in) :: numrots_ub
    integer(kind=int32), dimension(ubwmax_ub,n), intent(in) :: js_ub
    real(kind=dp), dimension(ubwmax_ub,n), intent(in) :: cs_ub, ss_ub

    real(kind=dp), dimension(n,lbwmax_bv+ubwmax_bv+1), intent(out) :: b_bv
    integer(kind=int32), dimension(n), intent(out) :: numrots_bv
    integer(kind=int32), dimension(n,ubwmax_bv), intent(out) :: ks_bv
    real(kind=dp), dimension(n,ubwmax_bv), intent(out) :: cs_bv, ss_bv
    integer(kind=int32), intent(out) :: error

    integer(kind=int32) :: j, k, k0, k1, ubw1, lbw1, d
    type(d_rotation) :: rot

    error = 0
    b_bv=0.0_dp; numrots_bv=0
    ss_bv=0.0_dp; cs_bv=0.0_dp
    ks_bv=0
    ubw1=ubw+1; lbw1=lbw
    
    if (n < 1) then
       error = 1
       return
    end if
    if (n == 1) then
       b_bv(1,1)=b_ub(1,1)
       return
    end if
    ! must allow for temporary fill-in
    if (lbwmax_ub+ubwmax_ub+1<ubw1+lbw1+1) then
       error = 2
       return
    end if
    call down_shift(b_ub)
    do k=1,n-1
       do j=1,numrots_ub(n-k)
          rot%cosine=cs_ub(j,n-k); rot%sine=ss_ub(j,n-k)
          call rotation_times_tbc(rot,b_ub,n,lbw1,ubw1,n-k,0,js_ub(j,n-k))
       end do
       k0=max(n-k+1,ubw+2)
       k1=min(n-k+ubw,n)
       numrots_bv(k+1)=k1-k0+1
       do j=k0,k1
          rot=rgivens(get_el_bc(b_ub,ubw1,j-ubw1,j-1), get_el_bc(b_ub,ubw1,j-ubw1,j))
          ks_bv(k+1,k1-j+1)=j-1
          cs_bv(k+1,k1-j+1)=rot%cosine; ss_bv(k+1,k1-j+1)=rot%sine
          call tbc_times_rotation(b_ub,n,lbw1,ubw1,0,k+1,rot,j-1)
       end do
    end do
    ! Store the results in b_bv
    call up_shift(b_ub)
    do d=1,lbw+1
       do j=1,n-lbw+d-1
          b_bv(lbw-d+j+1,d)=b_ub(ubw+lbw+2-d,j)
       end do
    end do
    do d=lbw+2,lbw+ubw+1
       do j=1,n-d+lbw+1
          b_bv(j,d)=b_ub(ubw+lbw+2-d,j+d-lbw-1)
       end do
    end do
  end subroutine f_d_convert_ub_to_bv

  subroutine c_convert_ub_to_bv(ub, bv, error)
    type(c_ub) :: ub
    type(c_bv) :: bv
    integer(kind=int32), intent(out) :: error
    if (get_n(ub) /= get_n(bv)) then
       error = 3; return
    end if
    call f_c_convert_ub_to_bv(ub%b, get_n(ub), ub%lbw, ub%ubw, get_lbwmax(ub), &
         get_ubwmax(ub), ub%numrotsu, ub%jsu, ub%csu, ub%ssu, bv%b, get_lbwmax(bv), &
         get_ubwmax(bv), bv%numrotsv, bv%ksv, bv%csv, bv%ssv, error)
    bv%lbw=ub%lbw; bv%ubw=ub%ubw
  end subroutine c_convert_ub_to_bv

  subroutine f_c_convert_ub_to_bv(b_ub, n, lbw, ubw, lbwmax_ub, ubwmax_ub, numrots_ub, &
       js_ub, cs_ub, ss_ub, b_bv, lbwmax_bv, ubwmax_bv, numrots_bv, ks_bv, cs_bv, &
       ss_bv, error)
    complex(kind=dp), dimension(lbwmax_ub+ubwmax_bv+1,n), intent(inout) :: b_ub
    integer(kind=int32), intent(in) :: n, lbw, ubw, lbwmax_ub, ubwmax_ub, lbwmax_bv, ubwmax_bv
    integer(kind=int32), dimension(n), intent(in) :: numrots_ub
    integer(kind=int32), dimension(ubwmax_ub,n), intent(in) :: js_ub
    complex(kind=dp), dimension(ubwmax_ub,n), intent(in) :: cs_ub, ss_ub

    complex(kind=dp), dimension(n,lbwmax_bv+ubwmax_bv+1), intent(out) :: b_bv
    integer(kind=int32), dimension(n), intent(out) :: numrots_bv
    integer(kind=int32), dimension(n,ubwmax_bv), intent(out) :: ks_bv
    complex(kind=dp), dimension(n,ubwmax_bv), intent(out) :: cs_bv, ss_bv
    integer(kind=int32), intent(out) :: error

    integer(kind=int32) :: j, k, k0, k1, ubw1, lbw1, d
    type(c_rotation) :: rot

    error = 0
    b_bv=(0.0_dp,0.0_dp); numrots_bv=0
    ss_bv=(0.0_dp, 0.0_dp); cs_bv=(0.0_dp, 0.0_dp)
    ks_bv=0
    ubw1=ubw+1; lbw1=lbw
    
    if (n < 1) then
       error = 1
       return
    end if
    if (n == 1) then
       b_bv(1,1)=b_ub(1,1)
       return
    end if
    ! must allow for temporary fill-in
    if (lbwmax_ub+ubwmax_ub+1<ubw1+lbw1+1) then
       error = 2
       return
    end if
    call down_shift(b_ub)
    do k=1,n-1
       do j=1,numrots_ub(n-k)
          rot%cosine=cs_ub(j,n-k); rot%sine=ss_ub(j,n-k)
          call rotation_times_tbc(rot,b_ub,n,lbw1,ubw1,n-k,0,js_ub(j,n-k))
       end do
       k0=max(n-k+1,ubw+2)
       k1=min(n-k+ubw,n)
       numrots_bv(k+1)=k1-k0+1
       do j=k0,k1
          rot=rgivens(get_el_bc(b_ub,ubw1,j-ubw1,j-1), get_el_bc(b_ub,ubw1,j-ubw1,j))
          ks_bv(k+1,k1-j+1)=j-1
          cs_bv(k+1,k1-j+1)=rot%cosine; ss_bv(k+1,k1-j+1)=rot%sine
          call tbc_times_rotation(b_ub,n,lbw1,ubw1,0,k+1,rot,j-1)
       end do
    end do
    ! Store the results in b_bv
    call up_shift(b_ub)
    do d=1,lbw+1
       do j=1,n-lbw+d-1
          b_bv(lbw-d+j+1,d)=b_ub(ubw+lbw+2-d,j)
       end do
    end do
    do d=lbw+2,lbw+ubw+1
       do j=1,n-d+lbw+1
          b_bv(j,d)=b_ub(ubw+lbw+2-d,j+d-lbw-1)
       end do
    end do
  end subroutine f_c_convert_ub_to_bv


end module conversions_ub_to_bv
