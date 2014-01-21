module conversions_bv_to_ub
use prec
use shift
use rotation
use band_types
use nested_types
implicit none

interface convert_bv_to_ub
   module procedure d_convert_bv_to_ub, c_convert_bv_to_ub
end interface convert_bv_to_ub

interface f_convert_bv_to_ub
   module procedure f_d_convert_bv_to_ub, f_c_convert_bv_to_ub
end interface f_convert_bv_to_ub

contains

  ! Errors:
  ! 0: no error
  ! 1: n<1
  ! 3: ub%n /= bv%n
  subroutine d_convert_bv_to_ub(bv, ub, error)
    type(d_bv) :: bv
    type(d_ub) :: ub
    integer(kind=int32), intent(out) :: error
    if (ub%n /= bv%n) then
       error = 3; return
    end if
    call f_d_convert_bv_to_ub(bv%b, bv%n, bv%lbw, bv%ubw, bv%lbwmax, bv%ubwmax, bv%numrotsv, bv%k1sv, bv%k2sv, &
         bv%csv, bv%ssv, ub%b,  ub%lbwmax, ub%ubwmax, ub%numrotsu, ub%j1su, ub%j2su, &
         ub%csu, ub%ssu, error)
    ub%lbw=bv%lbw; ub%ubw=bv%ubw; ub%n=bv%n
  end subroutine d_convert_bv_to_ub

  subroutine f_d_convert_bv_to_ub(b_bv, n, lbw, ubw, lbwmax_bv, ubwmax_bv, numrots_bv, k1s_bv, k2s_bv, cs_bv, ss_bv, &
       b_ub, lbwmax_ub, ubwmax_ub, numrots_ub, j1s_ub, j2s_ub, cs_ub, ss_ub, error)
    real(kind=dp), dimension(n,lbwmax_bv+ubwmax_bv+1), intent(inout) :: b_bv
    integer(kind=int32), intent(in) :: n, lbw, ubw, lbwmax_ub, ubwmax_ub, lbwmax_bv, ubwmax_bv
    integer(kind=int32), dimension(n), intent(in) :: numrots_bv
    integer(kind=int32), dimension(n,ubwmax_bv), intent(in) :: k1s_bv, k2s_bv
    real(kind=dp), dimension(n,ubwmax_bv), intent(in) :: cs_bv, ss_bv

    real(kind=dp), dimension(lbwmax_ub+ubwmax_ub+1,n), intent(out) :: b_ub
    integer(kind=int32), dimension(n), intent(out) :: numrots_ub
    integer(kind=int32), dimension(ubwmax_ub,n), intent(out) :: j1s_ub, j2s_ub
    real(kind=dp), dimension(ubwmax_ub,n), intent(out) :: cs_ub, ss_ub
    integer(kind=int32), intent(out) :: error

    integer(kind=int32) :: j, k, l, roffs, coffs, ubw1, lbw1, d
    type(d_rotation) :: rot

    error = 0
    b_ub=0.0_dp; numrots_ub=0
    ss_ub=0.0_dp; cs_ub=0.0_dp
    j1s_ub=0; j2s_ub=0
    ubw1=ubw+1; lbw1=lbw
    l = 1 ! The last row of b_bv is already revealed in b_ub
    
    if (n < 1) then
       error = 1
       return
    end if
    if (n == 1) then
       b_ub(1,1)=b_bv(1,1)
       return
    end if
    ! must allow for temporary fill-in
    if (lbwmax_bv+ubwmax_bv+1<ubw1+lbw1+1) then
       error = 2
       return
    end if
    b_bv(:,ubw1+lbw1+1)=0.0_dp
    ! apply v_{n-1}
    do j=1, numrots_bv(n-1)
       rot%cosine=cs_bv(n-1,j); rot%sine=ss_bv(n-1,j)
       call tbr_times_rotation(b_bv,n,lbw1,ubw1,n-1,trp_rot(rot),k1s_bv(n-1,j))
    end do
    ! Generate a (ubw+1) x ubw lower triangular matrix.
    do k=2,ubw1-1
       coffs=ubw1
       numrots_ub(k)=k-1
       do j=k-1,1,-1
          rot=lgivens2(get_el_br(b_bv,lbw1,j,coffs+j), get_el_br(b_bv,lbw1,j+1,coffs+j))
          cs_ub(j,k)=rot%cosine; ss_ub(j,k)=rot%sine
          j1s_ub(j,k)=j; j2s_ub(j,k)=j+1
          call rotation_times_tbr(trp_rot(rot),b_bv,n,lbw1,ubw1,k,j)
          call set_el_br(b_bv,lbw1,j,coffs+j, 0.0_dp)
       end do
       do j=1,numrots_bv(n-k)
          rot%cosine=cs_bv(n-k,j); rot%sine=ss_bv(n-k,j)
          call tbr_times_rotation(b_bv,n,lbw1,ubw1,n-k,trp_rot(rot),k1s_bv(n-k,j))
       end do
    end do
    do k=ubw1,n-ubw1
       roffs=k-ubw1; coffs=k
       numrots_ub(k)=ubw1-1
       do j=ubw1-1,1,-1
          rot=lgivens2(get_el_br(b_bv,lbw1,roffs+j,coffs+j), get_el_br(b_bv,lbw1,roffs+j+1,coffs+j))
          cs_ub(j,k)=rot%cosine; ss_ub(j,k)=rot%sine
          j1s_ub(j,k)=roffs+j; j2s_ub(j,k)=roffs+j+1
          call rotation_times_tbr(trp_rot(rot), b_bv,n,lbw1,ubw1, k, roffs+j)
          call set_el_br(b_bv,lbw1,roffs+j,coffs+j, 0.0_dp)
       end do
       do j=1,numrots_bv(n-k)
          rot%cosine=cs_bv(n-k,j); rot%sine=ss_bv(n-k,j)
          call tbr_times_rotation(b_bv,n,lbw1,ubw1,n-k,trp_rot(rot),k1s_bv(n-k,j))
       end do
    end do
    if (n-ubw>ubw) then ! if n-ubw=ubw, there is a square upper triangular matrix already.
       do k=n-ubw,n-1
          numrots_ub(k)=n-k
          coffs=k
          roffs=k-ubw1
          do j=n-k,1,-1
             rot=lgivens2(get_el_br(b_bv,lbw1,roffs+j,coffs+j),get_el_br(b_bv,lbw1,roffs+j+1,coffs+j))
             cs_ub(j,k)=rot%cosine; ss_ub(j,k)=rot%sine
             j1s_ub(j,k)=roffs+j; j2s_ub(j,k)=roffs+j+1
             call rotation_times_tbr(trp_rot(rot), b_bv,n,lbw1,ubw1, k, roffs+j)
             call set_el_br(b_bv,lbw1,roffs+j,coffs+j, 0.0_dp)
          end do
          do j=1,numrots_bv(n-k)
             rot%cosine=cs_bv(n-k,j); rot%sine=ss_bv(n-k,j)
             call tbr_times_rotation(b_bv,n,lbw1,ubw1,n-k,trp_rot(rot),k1s_bv(n-k,j))
          end do
       end do
    end if
    ! store the results in b_ub
    do d=1,lbw+1
       do j=1,n-lbw+d-1
          b_ub(ubw+lbw+2-d,j)=b_bv(lbw-d+j+1,d)
       end do
    end do
    do d=lbw+2,lbw+ubw+1
       do j=1,n-d+lbw+1
          b_ub(ubw+lbw+2-d,j+d-lbw-1)=b_bv(j,d)
       end do
    end do
  end subroutine f_d_convert_bv_to_ub

  subroutine c_convert_bv_to_ub(bv, ub, error)
    type(c_bv) :: bv
    type(c_ub) :: ub
    integer(kind=int32), intent(out) :: error
    if (ub%n /= bv%n) then
       error = 3; return
    end if
    call f_c_convert_bv_to_ub(bv%b, bv%n, bv%lbw, bv%ubw, bv%lbwmax, bv%ubwmax, bv%numrotsv, bv%k1sv, bv%k2sv, &
         bv%csv, bv%ssv, ub%b,  ub%lbwmax, ub%ubwmax, ub%numrotsu, ub%j1su, ub%j2su, &
         ub%csu, ub%ssu, error)
    ub%lbw=bv%lbw; ub%ubw=bv%ubw; ub%n=bv%n
  end subroutine c_convert_bv_to_ub

  subroutine f_c_convert_bv_to_ub(b_bv, n, lbw, ubw, lbwmax_bv, ubwmax_bv, numrots_bv, k1s_bv, k2s_bv, cs_bv, ss_bv, &
       b_ub, lbwmax_ub, ubwmax_ub, numrots_ub, j1s_ub, j2s_ub, cs_ub, ss_ub, error)
    complex(kind=dp), dimension(n,lbwmax_bv+ubwmax_bv+1), intent(inout) :: b_bv
    integer(kind=int32), intent(in) :: n, lbw, ubw, lbwmax_bv, ubwmax_bv, lbwmax_ub, ubwmax_ub
    integer(kind=int32), dimension(n), intent(in) :: numrots_bv
    integer(kind=int32), dimension(n,ubwmax_bv), intent(in) :: k1s_bv, k2s_bv
    complex(kind=dp), dimension(n,ubwmax_bv), intent(in) :: cs_bv, ss_bv

    complex(kind=dp), dimension(lbwmax_ub+ubwmax_ub+1,n), intent(out) :: b_ub
    integer(kind=int32), dimension(n), intent(out) :: numrots_ub
    integer(kind=int32), dimension(ubwmax_ub,n), intent(out) :: j1s_ub, j2s_ub
    complex(kind=dp), dimension(ubwmax_ub,n), intent(out) :: cs_ub, ss_ub
    integer(kind=int32), intent(out) :: error

    integer(kind=int32) :: j, k, l, roffs, coffs, ubw1, lbw1, d
    type(c_rotation) :: rot

    error = 0
    b_ub=(0.0_dp, 0.0_dp); numrots_ub=0
    ss_ub=(0.0_dp, 0.0_dp); cs_ub=(0.0_dp, 0.0_dp)
    j1s_ub=0; j2s_ub=0
    ubw1=ubw+1; lbw1=lbw
    l = 1 ! The last row of b_bv is already revealed in b_ub
    
    if (n < 1) then
       error = 1
       return
    end if
    if (n == 1) then
       b_ub(1,1)=b_bv(1,1)
       return
    end if
    ! must allow for temporary fill-in
    if (lbwmax_bv+ubwmax_bv+1<ubw1+lbw1+1) then
       error = 2
       return
    end if
    b_bv(:,ubw1+lbw1+1)=(0.0_dp, 0.0_dp)
    ! apply v_{n-1}
    do j=1, numrots_bv(n-1)
       rot%cosine=cs_bv(n-1,j); rot%sine=ss_bv(n-1,j)
       call tbr_times_rotation(b_bv,n,lbw1,ubw1,n-1,trp_rot(rot),k1s_bv(n-1,j))
    end do
    ! Generate a (ubw+1) x ubw lower triangular matrix.
    do k=2,ubw1-1
       coffs=ubw1
       numrots_ub(k)=k-1
       do j=k-1,1,-1
          rot=lgivens2(get_el_br(b_bv,lbw1,j,coffs+j), get_el_br(b_bv,lbw1,j+1,coffs+j))
          cs_ub(j,k)=rot%cosine; ss_ub(j,k)=rot%sine
          j1s_ub(j,k)=j; j2s_ub(j,k)=j+1
          call rotation_times_tbr(trp_rot(rot),b_bv,n,lbw1,ubw1,k,j)
          call set_el_br(b_bv,lbw1,j,coffs+j, (0.0_dp, 0.0_dp))
       end do
       do j=1,numrots_bv(n-k)
          rot%cosine=cs_bv(n-k,j); rot%sine=ss_bv(n-k,j)
          call tbr_times_rotation(b_bv,n,lbw1,ubw1,n-k,trp_rot(rot),k1s_bv(n-k,j))
       end do
    end do
    do k=ubw1,n-ubw1
       roffs=k-ubw1; coffs=k
       numrots_ub(k)=ubw1-1
       do j=ubw1-1,1,-1
          rot=lgivens2(get_el_br(b_bv,lbw1,roffs+j,coffs+j), get_el_br(b_bv,lbw1,roffs+j+1,coffs+j))
          cs_ub(j,k)=rot%cosine; ss_ub(j,k)=rot%sine
          j1s_ub(j,k)=roffs+j; j2s_ub(j,k)=roffs+j+1
          call rotation_times_tbr(trp_rot(rot), b_bv,n,lbw1,ubw1, k, roffs+j)
          call set_el_br(b_bv,lbw1,roffs+j,coffs+j, (0.0_dp, 0.0_dp))
       end do
       do j=1,numrots_bv(n-k)
          rot%cosine=cs_bv(n-k,j); rot%sine=ss_bv(n-k,j)
          call tbr_times_rotation(b_bv,n,lbw1,ubw1,n-k,trp_rot(rot),k1s_bv(n-k,j))
       end do
    end do
    if (n-ubw>ubw) then ! if n-ubw=ubw, there is a square upper triangular matrix already.
       do k=n-ubw,n-1
          numrots_ub(k)=n-k
          coffs=k
          roffs=k-ubw1
          do j=n-k,1,-1
             rot=lgivens2(get_el_br(b_bv,lbw1,roffs+j,coffs+j),get_el_br(b_bv,lbw1,roffs+j+1,coffs+j))
             cs_ub(j,k)=rot%cosine; ss_ub(j,k)=rot%sine
             j1s_ub(j,k)=roffs+j; j2s_ub(j,k)=roffs+j+1
             call rotation_times_tbr(trp_rot(rot), b_bv,n,lbw1,ubw1, k, roffs+j)
             call set_el_br(b_bv,lbw1,roffs+j,coffs+j, (0.0_dp, 0.0_dp))
          end do
          do j=1,numrots_bv(n-k)
             rot%cosine=cs_bv(n-k,j); rot%sine=ss_bv(n-k,j)
             call tbr_times_rotation(b_bv,n,lbw1,ubw1,n-k,trp_rot(rot),k1s_bv(n-k,j))
          end do
       end do
    end if
    ! store the results in b_ub
    do d=1,lbw+1
       do j=1,n-lbw+d-1
          b_ub(ubw+lbw+2-d,j)=b_bv(lbw-d+j+1,d)
       end do
    end do
    do d=lbw+2,lbw+ubw+1
       do j=1,n-d+lbw+1
          b_ub(ubw+lbw+2-d,j+d-lbw-1)=b_bv(j,d)
       end do
    end do
  end subroutine f_c_convert_bv_to_ub

end module conversions_bv_to_ub
