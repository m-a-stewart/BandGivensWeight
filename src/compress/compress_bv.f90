module compress_bv
use prec
use shift
use rotation
use band_types
use nested_types
use nullvec
implicit none

integer(kind=int32), private, parameter :: nullmaxits=5

interface compress_bv_to_ub_1
   module procedure d_compress_bv_to_ub_1, c_compress_bv_to_ub_1
end interface compress_bv_to_ub_1

interface f_compress_bv_to_ub_1
   module procedure f_d_compress_bv_to_ub_1, f_c_compress_bv_to_ub_1
end interface f_compress_bv_to_ub_1

contains

  ! Errors:
  ! 0: no error
  ! 1: n<1
  ! 3: ub%n /= bv%n
  ! 4: compression error, failed to get a null vector at some point.
  ! 
  ! told governs whether a diagonal of L is considered small enough to move to the
  ! lower right corner of L.  If this tolerance is not met, the algorithm looks
  ! for a null vector.  If both are zero the result is a forced compression.
  ! The algorithm will use the maximum number of iterations to get the best
  ! null vector possible and truncate whatever appears in the lower right corner
  ! of L, whether it is small or not.
  !
  subroutine d_compress_bv_to_ub_1(bv, ub, told, tol, error)
    type(d_ub) :: ub
    type(d_bv) :: bv
    integer(kind=int32), intent(out) :: error
    real(kind=dp), intent(in) :: tol, told
    if (ub%n /= bv%n) then
       error = 3; return
    end if
    call f_d_compress_bv_to_ub_1(bv%b, bv%n, bv%lbw, bv%ubw, bv%lbwmax, bv%ubwmax, bv%numrotsv, bv%k1sv, bv%k2sv, &
         bv%csv, bv%ssv, ub%b, ub%lbwmax, ub%ubwmax, ub%numrotsu, ub%j1su, ub%j2su, ub%csu, ub%ssu, told, tol, &
         error)
    if (error == 0) then
       ub%lbw=bv%lbw; ub%ubw=bv%ubw-1; ub%n=bv%n
    end if
  end subroutine d_compress_bv_to_ub_1

  subroutine f_d_compress_bv_to_ub_1(b_bv, n, lbw, ubw, lbwmax_bv, ubwmax_bv, numrots_bv, &
       k1s_bv, k2s_bv, cs_bv, ss_bv, &
       b_ub, lbwmax_ub, ubwmax_ub, numrots_ub, j1s_ub, j2s_ub, cs_ub, ss_ub, told, tol, error)
    integer(kind=int32), intent(in) :: n, lbw, ubw, lbwmax_ub, ubwmax_ub, lbwmax_bv, ubwmax_bv
    real(kind=dp), dimension(n,lbwmax_bv+ubwmax_bv+1), intent(inout) :: b_bv
    integer(kind=int32), dimension(n), intent(in) :: numrots_bv
    integer(kind=int32), dimension(n,ubwmax_bv), intent(in) :: k1s_bv, k2s_bv
    real(kind=dp), dimension(n,ubwmax_bv), intent(in) :: cs_bv, ss_bv

    real(kind=dp), dimension(lbwmax_ub+ubwmax_ub+1,n), intent(out) :: b_ub
    integer(kind=int32), dimension(n), intent(out) :: numrots_ub
    integer(kind=int32), dimension(ubwmax_ub,n), intent(out) :: j1s_ub, j2s_ub
    real(kind=dp), dimension(ubwmax_ub,n), intent(out) :: cs_ub, ss_ub

    real(kind=dp), intent(in) :: tol, told
    integer(kind=int32), intent(out) :: error

    integer(kind=int32) :: j, k, jj, roffs, coffs, ubw1, lbw1, nullerr, d, minindex
    type(d_rotation) :: rot
    real(kind=dp), target, dimension(ubw+1,ubw+1) :: q
    real(kind=dp), dimension(ubw+1,ubw+1) :: l
    real(kind=dp), dimension(ubw+1) :: x
    real(kind=dp) :: nrma, tmp, mindiag
    real(kind=dp), pointer, dimension(:,:) :: pq

    error = 0
    numrots_ub=0
    ss_ub=0.0_dp; cs_ub=0.0_dp
    j1s_ub=0; j2s_ub=0
    ubw1=ubw+1; lbw1=lbw
    nrma = maxabs(b_bv)*sqrt(real(n))
    
    ! store a copy in case compression fails.
    b_ub(1:ubw+lbw+1,1:n)=transpose(b_bv(1:n,1:ubw+lbw+1))
    if (n < 1) then
       error = 1
       return
    end if
    if (n == 1) then
        b_ub(1,1)=b_bv(1,1)
        return
    end if
    ! must allow for temporary fill-in of one extra superdiagonal in b_bv.
    if (lbwmax_bv+ubwmax_bv+1<ubw1+lbw1+1) then
       error = 2
       return
    end if
    q=0.0_dp
    do j=1,ubw1
       q(j,j)=1.0_dp
    end do
    ! provide working space in a superdiagonal on the right of b_bv
    b_bv(:,ubw1+lbw1+1)=0.0_dp
    ! apply v_{n-1}
    do k=1,numrots_bv(n-1)
       rot%cosine=cs_bv(n-1,k); rot%sine=ss_bv(n-1,k)
       call tbr_times_rotation(b_bv,n,lbw1,ubw1,n-1,trp_rot(rot),k1s_bv(n-1,k))
    end do
    ! Apply v_{n-2} to v_{n-ubw+1}
    ! Generate a ubw x ubw lower Hessenberg matrix
    do j=2,ubw-1
       numrots_ub(j) = j-1
       do k=j,2,-1
          rot=lgivens2(get_el_br(b_bv,lbw1,k-1,ubw1+k-1),get_el_br(b_bv,lbw1,k,ubw1+k-1))
          cs_ub(k-1,j)=rot%cosine; ss_ub(k-1,j)=rot%sine
          j1s_ub(k-1,j)=k-1; j2s_ub(k-1,j)=k
          call rotation_times_tbr(trp_rot(rot),b_bv,n,lbw1,ubw1,j,k-1)
          call set_el_br(b_bv,lbw1,k-1,ubw1+k-1,0.0_dp)
       end do
       do k=1,numrots_bv(n-j)
          rot%cosine=cs_bv(n-j,k); rot%sine=ss_bv(n-j,k)
          call tbr_times_rotation(b_bv,n,lbw1,ubw1,n-j,trp_rot(rot),k1s_bv(n-j,k))
       end do
    end do
    ! At this point there is a ubw x ubw lower Hessenberg matrix in
    ! A(1:ubw,ubw+1:2*ubw)
    ! Compute an LQ decomposition.
    coffs=ubw
    do k=1,ubw-1
       rot=rgivens(get_el_br(b_bv,lbw1,k,coffs+k),get_el_br(b_bv,lbw1,k,coffs+k+1))
       call tbr_times_rotation(b_bv,n,lbw1,ubw1,n-ubw,rot,coffs+k)
       call set_el_br(b_bv,lbw1,k,coffs+k+1,0.0_dp)
       call rotation_times_general(trp_rot(rot),q,k,k+1)
    end do
    ! main loop
    do k=ubw,n-ubw-1
       ! There should be a lower triangular matrix in
       ! A(k-ubw+1:k,k+1:k+ubw) representing the L factor in an LQ decomposition.
       roffs=k-ubw; coffs=k
       mindiag=abs(get_el_br(b_bv,lbw1,roffs+ubw,coffs+ubw))
       minindex=ubw
       do j=ubw-1,1,-1
          tmp=abs(get_el_br(b_bv,lbw1,roffs+j,coffs+j))
          if (tmp <= mindiag) then
             minindex=j
             mindiag=tmp
          end if
       end do
       if (mindiag < told*nrma) then
          numrots_ub(k) = minindex-1
          do j=minindex,2,-1
             rot=lgivens2(get_el_br(b_bv,lbw1,roffs+j-1,coffs+j-1),get_el_br(b_bv,lbw1,roffs+j,coffs+j-1))
             call rotation_times_tbr(trp_rot(rot), b_bv,n,lbw1,ubw1,coffs,roffs+j-1)
             cs_ub(j-1,k)=rot%cosine; ss_ub(j-1,k)=rot%sine
             j1s_ub(j-1,k)=roffs+j-1; j2s_ub(j-1,k)=roffs+j
             call set_el_br(b_bv,lbw1,roffs+j-1,coffs+j-1,0.0_dp)
             ! swap columns in L
             do jj=j-1,ubw
                tmp=get_el_br(b_bv,lbw1,roffs+jj,coffs+j)
                call set_el_br(b_bv,lbw1,roffs+jj,coffs+j, &
                     get_el_br(b_bv,lbw1,roffs+jj,coffs+j-1))
                call set_el_br(b_bv,lbw1,roffs+jj,coffs+j-1,tmp)
             end do
             do jj=1,ubw1
                tmp=q(j,jj)
                q(j,jj)=q(j-1,jj)
                q(j-1,jj)=tmp
             end do
          end do
          call set_el_br(b_bv,lbw1,roffs+1,coffs+1,0.0_dp)
       else
          ! Need a more general null vector of L.
          call submatrix_br(b_bv,lbw1,ubw1, roffs+1,roffs+ubw,coffs+1,coffs+ubw,l(1:ubw,1:ubw))
          call f_d_lower_left_nullvec(x(1:ubw),l(1:ubw,1:ubw),tol*nrma,nullmaxits, nullerr)
          if (nullerr >= 0) then ! nullvector found
             ! Introduce a zero while preserving the LQ factorization.
             numrots_ub(k)=ubw-1
             do j=ubw-1,1,-1
                rot=rgivens(x(j),x(j+1))
                call general_times_rotation(x,rot,j,j+1)
                call rotation_times_tbr(trp_rot(rot),b_bv,n,lbw1,ubw1,coffs,roffs+j)
                cs_ub(j,k)=rot%cosine; ss_ub(j,k)=rot%sine
                j1s_ub(j,k)=roffs+j; j2s_ub(j,k)=roffs+j+1
                rot=rgivens(get_el_br(b_bv,lbw1,roffs+j,coffs+j),get_el_br(b_bv,lbw1,roffs+j,coffs+j+1))
                call tbr_times_rotation(b_bv,n,lbw1,ubw1,n-k,rot,coffs+j)
                call rotation_times_general(trp_rot(rot),q,j,j+1)
                call set_el_br(b_bv,lbw1,roffs+j,coffs+j+1,0.0_dp)
             end do
             call set_el_br(b_bv,lbw1,roffs+1,coffs+1,0.0_dp)
          else
             ! Compression has failed.  restore b_ub and return.
             error=4
             b_bv(1:n,1:ubw+lbw+1)=transpose(b_ub(1:ubw+lbw+1,1:n))
             return
          end if
       end if
       ! Apply v_{n-k} to Q.
       do j=1,numrots_bv(n-k)
          rot%cosine=cs_bv(n-k,j); rot%sine=ss_bv(n-k,j)
          call general_times_rotation(q,trp_rot(rot),k1s_bv(n-k,j)-coffs,k2s_bv(n-k,j)-coffs)
       end do
       ! include row k+1 in the LQ factorization
       b_bv(k+1,lbw1+1:lbw1+ubw1)=matmul(b_bv(k+1,lbw1+1:lbw1+ubw1),transpose(q))
       ! downdate
       do j=ubw,1,-1
          rot=lgivens(q(j,1),q(j+1,1))
          call rotation_times_general(trp_rot(rot),q,j,j+1)
          q(j+1,1)=0.0_dp
          call tbr_times_rotation(b_bv,n,lbw1,ubw1,n-k-1,rot,coffs+j)
       end do
       do j=1,ubw
          tmp=get_el_br(b_bv,lbw1,roffs+j+1,coffs+1)
          call set_el_br(b_bv,lbw1,roffs+j+1,coffs+1,tmp*q(1,1))
       end do
       call up_left_shift(q)
       q(ubw1,ubw1)=1.0_dp
       coffs=coffs+1
       roffs=roffs+1
       do j=1,ubw-1
          rot=rgivens(get_el_br(b_bv,lbw1,roffs+j,coffs+j),get_el_br(b_bv,lbw1,roffs+j,coffs+j+1))
          call tbr_times_rotation(b_bv,n,lbw1,ubw1,n-k-1,rot,coffs+j)
          call rotation_times_general(trp_rot(rot),q,j,j+1)
          call set_el_br(b_bv,lbw1,roffs+j,coffs+j+1,0.0_dp)
       end do
    end do
    !
    ! Intermediate step: L is still square, as is q.  We need to find a null vector.
    ! corresponds to j=n-ubw, roffs=0, coffs=ubw.  At the end of this section
    ! L will be trapezoidal.
    !
    roffs=n-2*ubw; coffs=n-ubw
    pq => q(1:ubw,1:ubw)
    call submatrix_br(b_bv,lbw1,ubw1,roffs+1,roffs+ubw,coffs+1,coffs+ubw,l)
    call f_d_lower_left_nullvec(x(1:ubw),l(1:ubw,1:ubw),tol*nrma,nullmaxits, nullerr)
    if (nullerr >= 0) then
       numrots_ub(n-ubw)=ubw-1
       do j=ubw,2,-1
          rot=rgivens(x(j-1),x(j))
          call general_times_rotation(x,rot,j-1,j)
          call rotation_times_tbr(trp_rot(rot),b_bv,n,lbw1,ubw1,n-ubw,roffs+j-1)
          cs_ub(j-1,n-ubw)=rot%cosine; ss_ub(j-1,n-ubw)=rot%sine
          j1s_ub(j-1,n-ubw)=roffs+j-1; j2s_ub(j-1,n-ubw)=roffs+j
          rot=rgivens(get_el_br(b_bv,lbw1,roffs+j-1,coffs+j-1), &
               get_el_br(b_bv,lbw1,roffs+j-1,coffs+j))
          call tbr_times_rotation(b_bv,n,lbw1,ubw1,ubw,rot,coffs+j-1)
          call rotation_times_general(trp_rot(rot), pq, j-1,j)
          call set_el_br(b_bv,lbw1,roffs+j-1,coffs+j,0.0_dp)
       end do
       call set_el_br(b_bv,lbw1,roffs+1,coffs+1,0.0_dp)
       ! apply v_{ubw} to q
       do j=1,numrots_bv(ubw)
          rot%cosine=cs_bv(ubw,j); rot%sine=ss_bv(ubw,j)
          call general_times_rotation(pq,trp_rot(rot),k1s_bv(ubw,j)-coffs,k2s_bv(ubw,j)-coffs)
       end do
       ! Include row n-ubw+1
       b_bv(n-ubw+1,lbw1+1:lbw1+ubw)=matmul(b_bv(n-ubw+1,lbw1+1:lbw1+ubw),transpose(pq))
       ! restore triangularity
       roffs=roffs+1
       do j=1,ubw-1
          rot=rgivens(get_el_br(b_bv,lbw1,roffs+j,coffs+j),get_el_br(b_bv,lbw1,roffs+j,coffs+j+1))
          call tbr_times_rotation(b_bv,n,lbw1,ubw1,ubw-1,rot,coffs+j)
          call rotation_times_general(trp_rot(rot),pq,j,j+1)
       end do
       ! downdate
       do j=ubw-1,1,-1
          rot=lgivens(pq(j,1),pq(j+1,1))
          call rotation_times_general(trp_rot(rot),pq,j,j+1)
          pq(j+1,1)=0.0_dp
          call tbr_times_rotation(b_bv,n,lbw1,ubw1,ubw-1,rot,coffs+j)
       end do
       do j=1,ubw
          tmp=get_el_br(b_bv,lbw1,roffs+j,coffs+1)
          call set_el_br(b_bv,lbw1,roffs+j,coffs+1,tmp*pq(1,1))
       end do
       call up_left_shift(pq)
    else
       ! Compression has failed.  restore b_ub and return.
       error=4
       b_bv(1:n,1:ubw+lbw+1)=transpose(b_ub(1:ubw+lbw+1,1:n))
       return
    end if
    ! Final stage:
    do k=n-ubw+1,n-1
       coffs=k; roffs=k-ubw
       pq=>q(1:n-k,1:n-k)
       ! apply u_k
       numrots_ub(k)=n-k
       do j=n-k+1,2,-1
          rot=lgivens2(get_el_br(b_bv,lbw1,roffs+j-1,coffs+j-1),get_el_br(b_bv,lbw1,roffs+j,coffs+j-1))
          call rotation_times_tbr(trp_rot(rot),b_bv,n,lbw1,ubw1,k,roffs+j-1)
          cs_ub(j-1,k)=rot%cosine; ss_ub(j-1,k)=rot%sine
          j1s_ub(j-1,k)=roffs+j-1; j2s_ub(j-1,k)=roffs+j
       end do
       ! Apply v_{n-k} to Q
       do j=1,numrots_bv(n-k)
          rot%cosine=cs_bv(n-k,j); rot%sine=ss_bv(n-k,j)
          call general_times_rotation(pq,trp_rot(rot),k1s_bv(n-k,j)-coffs,k2s_bv(n-k,j)-coffs)
       end do
       ! include row k+1 in the decomposition
       b_bv(k+1,lbw1+1:lbw1+n-k)=matmul(b_bv(k+1,lbw1+1:lbw1+n-k),transpose(pq))
       ! downdate
       do j=n-k-1,1,-1
          rot=lgivens(pq(j,1),pq(j+1,1))
          call rotation_times_general(trp_rot(rot),pq,j,j+1)
          call tbr_times_rotation(b_bv,n,lbw1,ubw1,n-k-1,rot,coffs+j)
       end do
       do j=2,ubw1
          tmp=get_el_br(b_bv,lbw1,roffs+j,coffs+1)
          call set_el_br(b_bv,lbw1,roffs+j,coffs+1,pq(1,1)*tmp)
       end do
       call up_left_shift(pq)
    end do
    do d=1,lbw+1
       do j=1,n-lbw+d-1
          b_ub(ubw-1+lbw+2-d,j)=b_bv(lbw-d+j+1,d)
       end do
    end do
    do d=lbw+2,lbw+ubw
       do j=1,n-d+lbw+1
          b_ub(ubw-1+lbw+2-d,j+d-lbw-1)=b_bv(j,d)
       end do
    end do
  end subroutine f_d_compress_bv_to_ub_1


! Complex

  subroutine c_compress_bv_to_ub_1(bv, ub, told, tol, error)
    type(c_ub) :: ub
    type(c_bv) :: bv
    integer(kind=int32), intent(out) :: error
    real(kind=dp), intent(in) :: tol, told
    if (ub%n /= bv%n) then
       error = 3; return
    end if
    call f_c_compress_bv_to_ub_1(bv%b, bv%n, bv%lbw, bv%ubw, bv%lbwmax, bv%ubwmax, bv%numrotsv, bv%k1sv, bv%k2sv, &
         bv%csv, bv%ssv, ub%b, ub%lbwmax, ub%ubwmax, ub%numrotsu, ub%j1su, ub%j2su, ub%csu, ub%ssu, told, tol, &
         error)
    if (error == 0) then
       ub%lbw=bv%lbw; ub%ubw=bv%ubw-1; ub%n=bv%n
    end if
  end subroutine c_compress_bv_to_ub_1

  subroutine f_c_compress_bv_to_ub_1(b_bv, n, lbw, ubw, lbwmax_bv, ubwmax_bv, numrots_bv, &
       k1s_bv, k2s_bv, cs_bv, ss_bv, &
       b_ub, lbwmax_ub, ubwmax_ub, numrots_ub, j1s_ub, j2s_ub, cs_ub, ss_ub, told, tol, error)
    integer(kind=int32), intent(in) :: n, lbw, ubw, lbwmax_ub, ubwmax_ub, lbwmax_bv, ubwmax_bv
    complex(kind=dp), dimension(n,lbwmax_bv+ubwmax_bv+1), intent(inout) :: b_bv
    integer(kind=int32), dimension(n), intent(in) :: numrots_bv
    integer(kind=int32), dimension(n,ubwmax_bv), intent(in) :: k1s_bv, k2s_bv
    complex(kind=dp), dimension(n,ubwmax_bv), intent(in) :: cs_bv, ss_bv

    complex(kind=dp), dimension(lbwmax_ub+ubwmax_ub+1,n), intent(out) :: b_ub
    integer(kind=int32), dimension(n), intent(out) :: numrots_ub
    integer(kind=int32), dimension(ubwmax_ub,n), intent(out) :: j1s_ub, j2s_ub
    complex(kind=dp), dimension(ubwmax_ub,n), intent(out) :: cs_ub, ss_ub

    real(kind=dp), intent(in) :: tol, told
    integer(kind=int32), intent(out) :: error

    integer(kind=int32) :: j, k, jj, roffs, coffs, ubw1, lbw1, nullerr, d, minindex
    type(c_rotation) :: rot
    complex(kind=dp), target, dimension(ubw+1,ubw+1) :: q
    complex(kind=dp), dimension(ubw+1,ubw+1) :: l
    complex(kind=dp), dimension(ubw+1) :: x
    complex(kind=dp) :: tmp
    real(kind=dp) :: nrma, tmpr, mindiag
    complex(kind=dp), pointer, dimension(:,:) :: pq

    error = 0
    numrots_ub=0
    ss_ub=(0.0_dp,0.0_dp); cs_ub=(0.0_dp, 0.0_dp)
    j1s_ub=0; j2s_ub=0
    ubw1=ubw+1; lbw1=lbw
    nrma = maxabs(b_bv)*sqrt(real(n))
    
    ! store a copy in case compression fails.
    b_ub(1:ubw+lbw+1,1:n)=transpose(b_bv(1:n,1:ubw+lbw+1))
    if (n < 1) then
       error = 1
       return
    end if
    if (n == 1) then
        b_ub(1,1)=b_bv(1,1)
        return
    end if
    ! must allow for temporary fill-in of one extra superdiagonal in b_bv.
    if (lbwmax_bv+ubwmax_bv+1<ubw1+lbw1+1) then
       error = 2
       return
    end if
    q=(0.0_dp,0.0_dp)
    do j=1,ubw1
       q(j,j)=(1.0_dp,0.0_dp)
    end do
    ! provide working space in a superdiagonal on the right of b_bv
    b_bv(:,ubw1+lbw1+1)=(0.0_dp,0.0_dp)
    ! apply v_{n-1}
    do k=1,numrots_bv(n-1)
       rot%cosine=cs_bv(n-1,k); rot%sine=ss_bv(n-1,k)
       call tbr_times_rotation(b_bv,n,lbw1,ubw1,n-1,trp_rot(rot),k1s_bv(n-1,k))
    end do
    ! Apply v_{n-2} to v_{n-ubw+1}
    ! Generate a ubw x ubw lower Hessenberg matrix
    do j=2,ubw-1
       numrots_ub(j) = j-1
       do k=j,2,-1
          rot=lgivens2(get_el_br(b_bv,lbw1,k-1,ubw1+k-1),get_el_br(b_bv,lbw1,k,ubw1+k-1))
          cs_ub(k-1,j)=rot%cosine; ss_ub(k-1,j)=rot%sine
          j1s_ub(k-1,j)=k-1; j2s_ub(k-1,j)=k
          call rotation_times_tbr(trp_rot(rot),b_bv,n,lbw1,ubw1,j,k-1)
          call set_el_br(b_bv,lbw1,k-1,ubw1+k-1,(0.0_dp,0.0_dp))
       end do
       do k=1,numrots_bv(n-j)
          rot%cosine=cs_bv(n-j,k); rot%sine=ss_bv(n-j,k)
          call tbr_times_rotation(b_bv,n,lbw1,ubw1,n-j,trp_rot(rot),k1s_bv(n-j,k))
       end do
    end do
    ! At this point there is a ubw x ubw lower Hessenberg matrix in
    ! A(1:ubw,ubw+1:2*ubw)
    ! Compute an LQ decomposition.
    coffs=ubw
    do k=1,ubw-1
       rot=rgivens(get_el_br(b_bv,lbw1,k,coffs+k),get_el_br(b_bv,lbw1,k,coffs+k+1))
       call tbr_times_rotation(b_bv,n,lbw1,ubw1,n-ubw,rot,coffs+k)
       call set_el_br(b_bv,lbw1,k,coffs+k+1,(0.0_dp,0.0_dp))
       call rotation_times_general(trp_rot(rot),q,k,k+1)
    end do
    ! main loop
    do k=ubw,n-ubw-1
       ! There should be a lower triangular matrix in
       ! A(k-ubw+1:k,k+1:k+ubw) representing the L factor in an LQ decomposition.
       roffs=k-ubw; coffs=k
       mindiag=abs(get_el_br(b_bv,lbw1,roffs+ubw,coffs+ubw))
       minindex=ubw
       do j=ubw-1,1,-1
          tmpr=abs(get_el_br(b_bv,lbw1,roffs+j,coffs+j))
          if (tmpr <= mindiag) then
             minindex=j
             mindiag=tmpr
          end if
       end do
       if (mindiag < told*nrma) then
          numrots_ub(k) = minindex-1
          do j=minindex,2,-1
             rot=lgivens2(get_el_br(b_bv,lbw1,roffs+j-1,coffs+j-1),get_el_br(b_bv,lbw1,roffs+j,coffs+j-1))
             call rotation_times_tbr(trp_rot(rot), b_bv,n,lbw1,ubw1,coffs,roffs+j-1)
             cs_ub(j-1,k)=rot%cosine; ss_ub(j-1,k)=rot%sine
             j1s_ub(j-1,k)=roffs+j-1; j2s_ub(j-1,k)=roffs+j
             call set_el_br(b_bv,lbw1,roffs+j-1,coffs+j-1,(0.0_dp,0.0_dp))
             ! swap columns in L
             do jj=j-1,ubw
                tmp=get_el_br(b_bv,lbw1,roffs+jj,coffs+j)
                call set_el_br(b_bv,lbw1,roffs+jj,coffs+j, &
                     get_el_br(b_bv,lbw1,roffs+jj,coffs+j-1))
                call set_el_br(b_bv,lbw1,roffs+jj,coffs+j-1,tmp)
             end do
             do jj=1,ubw1
                tmp=q(j,jj)
                q(j,jj)=q(j-1,jj)
                q(j-1,jj)=tmp
             end do
          end do
          call set_el_br(b_bv,lbw1,roffs+1,coffs+1,(0.0_dp,0.0_dp))
       else
          ! Need a more general null vector of L.
          call submatrix_br(b_bv,lbw1,ubw1, roffs+1,roffs+ubw,coffs+1,coffs+ubw,l(1:ubw,1:ubw))
          call f_c_lower_left_nullvec(x(1:ubw),l(1:ubw,1:ubw),tol*nrma,nullmaxits, nullerr)
          if (nullerr >= 0) then ! nullvector found
             ! Introduce a zero while preserving the LQ factorization.
             numrots_ub(k)=ubw-1
             do j=ubw-1,1,-1
                rot=rgivens(x(j),x(j+1))
                call general_times_rotation(x,rot,j,j+1)
                call rotation_times_tbr(trp_rot(rot),b_bv,n,lbw1,ubw1,coffs,roffs+j)
                cs_ub(j,k)=rot%cosine; ss_ub(j,k)=rot%sine
                j1s_ub(j,k)=roffs+j; j2s_ub(j,k)=roffs+j+1
                rot=rgivens(get_el_br(b_bv,lbw1,roffs+j,coffs+j),get_el_br(b_bv,lbw1,roffs+j,coffs+j+1))
                call tbr_times_rotation(b_bv,n,lbw1,ubw1,n-k,rot,coffs+j)
                call rotation_times_general(trp_rot(rot),q,j,j+1)
                call set_el_br(b_bv,lbw1,roffs+j,coffs+j+1,(0.0_dp,0.0_dp))
             end do
             call set_el_br(b_bv,lbw1,roffs+1,coffs+1,(0.0_dp,0.0_dp))
          else
             ! Compression has failed.  restore b_ub and return.
             error=4
             b_bv(1:n,1:ubw+lbw+1)=transpose(b_ub(1:ubw+lbw+1,1:n))
             return
          end if
       end if
       ! Apply v_{n-k} to Q.
       do j=1,numrots_bv(n-k)
          rot%cosine=cs_bv(n-k,j); rot%sine=ss_bv(n-k,j)
          call general_times_rotation(q,trp_rot(rot),k1s_bv(n-k,j)-coffs,k2s_bv(n-k,j)-coffs)
       end do
       ! include row k+1 in the LQ factorization
       b_bv(k+1,lbw1+1:lbw1+ubw1)=matmul(b_bv(k+1,lbw1+1:lbw1+ubw1),conjg(transpose(q)))
       ! downdate
       do j=ubw,1,-1
          rot=lgivens(q(j,1),q(j+1,1))
          call rotation_times_general(trp_rot(rot),q,j,j+1)
          q(j+1,1)=(0.0_dp,0.0_dp)
          call tbr_times_rotation(b_bv,n,lbw1,ubw1,n-k-1,rot,coffs+j)
       end do
       do j=1,ubw
          tmp=get_el_br(b_bv,lbw1,roffs+j+1,coffs+1)
          call set_el_br(b_bv,lbw1,roffs+j+1,coffs+1,tmp*q(1,1))
       end do
       call up_left_shift(q)
       q(ubw1,ubw1)=(1.0_dp,0.0_dp)
       coffs=coffs+1
       roffs=roffs+1
       do j=1,ubw-1
          rot=rgivens(get_el_br(b_bv,lbw1,roffs+j,coffs+j),get_el_br(b_bv,lbw1,roffs+j,coffs+j+1))
          call tbr_times_rotation(b_bv,n,lbw1,ubw1,n-k-1,rot,coffs+j)
          call rotation_times_general(trp_rot(rot),q,j,j+1)
          call set_el_br(b_bv,lbw1,roffs+j,coffs+j+1,(0.0_dp,0.0_dp))
       end do
    end do
    !
    ! Intermediate step: L is still square, as is q.  We need to find a null vector.
    ! corresponds to j=n-ubw, roffs=0, coffs=ubw.  At the end of this section
    ! L will be trapezoidal.
    !
    roffs=n-2*ubw; coffs=n-ubw
    pq => q(1:ubw,1:ubw)
    call submatrix_br(b_bv,lbw1,ubw1,roffs+1,roffs+ubw,coffs+1,coffs+ubw,l)
    call f_c_lower_left_nullvec(x(1:ubw),l(1:ubw,1:ubw),tol*nrma,nullmaxits, nullerr)
    if (nullerr >= 0) then
       numrots_ub(n-ubw)=ubw-1
       do j=ubw,2,-1
          rot=rgivens(x(j-1),x(j))
          call general_times_rotation(x,rot,j-1,j)
          call rotation_times_tbr(trp_rot(rot),b_bv,n,lbw1,ubw1,n-ubw,roffs+j-1)
          cs_ub(j-1,n-ubw)=rot%cosine; ss_ub(j-1,n-ubw)=rot%sine
          j1s_ub(j-1,n-ubw)=roffs+j-1; j2s_ub(j-1,n-ubw)=roffs+j
          rot=rgivens(get_el_br(b_bv,lbw1,roffs+j-1,coffs+j-1), &
               get_el_br(b_bv,lbw1,roffs+j-1,coffs+j))
          call tbr_times_rotation(b_bv,n,lbw1,ubw1,ubw,rot,coffs+j-1)
          call rotation_times_general(trp_rot(rot), pq, j-1,j)
          call set_el_br(b_bv,lbw1,roffs+j-1,coffs+j,(0.0_dp,0.0_dp))
       end do
       call set_el_br(b_bv,lbw1,roffs+1,coffs+1,(0.0_dp,0.0_dp))
       ! apply v_{ubw} to q
       do j=1,numrots_bv(ubw)
          rot%cosine=cs_bv(ubw,j); rot%sine=ss_bv(ubw,j)
          call general_times_rotation(pq,trp_rot(rot),k1s_bv(ubw,j)-coffs,k2s_bv(ubw,j)-coffs)
       end do
       ! Include row n-ubw+1
       b_bv(n-ubw+1,lbw1+1:lbw1+ubw)=matmul(b_bv(n-ubw+1,lbw1+1:lbw1+ubw),conjg(transpose(pq)))
       ! restore triangularity
       roffs=roffs+1
       do j=1,ubw-1
          rot=rgivens(get_el_br(b_bv,lbw1,roffs+j,coffs+j),get_el_br(b_bv,lbw1,roffs+j,coffs+j+1))
          call tbr_times_rotation(b_bv,n,lbw1,ubw1,ubw-1,rot,coffs+j)
          call rotation_times_general(trp_rot(rot),pq,j,j+1)
       end do
       ! downdate
       do j=ubw-1,1,-1
          rot=lgivens(pq(j,1),pq(j+1,1))
          call rotation_times_general(trp_rot(rot),pq,j,j+1)
          pq(j+1,1)=(0.0_dp,0.0_dp)
          call tbr_times_rotation(b_bv,n,lbw1,ubw1,ubw-1,rot,coffs+j)
       end do
       do j=1,ubw
          tmp=get_el_br(b_bv,lbw1,roffs+j,coffs+1)
          call set_el_br(b_bv,lbw1,roffs+j,coffs+1,tmp*pq(1,1))
       end do
       call up_left_shift(pq)
    else
       ! Compression has failed.  restore b_ub and return.
       error=4
       b_bv(1:n,1:ubw+lbw+1)=transpose(b_ub(1:ubw+lbw+1,1:n))
       return
    end if
    ! Final stage:
    do k=n-ubw+1,n-1
       coffs=k; roffs=k-ubw
       pq=>q(1:n-k,1:n-k)
       ! apply u_k
       numrots_ub(k)=n-k
       do j=n-k+1,2,-1
          rot=lgivens2(get_el_br(b_bv,lbw1,roffs+j-1,coffs+j-1),get_el_br(b_bv,lbw1,roffs+j,coffs+j-1))
          call rotation_times_tbr(trp_rot(rot),b_bv,n,lbw1,ubw1,k,roffs+j-1)
          cs_ub(j-1,k)=rot%cosine; ss_ub(j-1,k)=rot%sine
          j1s_ub(j-1,k)=roffs+j-1; j2s_ub(j-1,k)=roffs+j
       end do
       ! Apply v_{n-k} to Q
       do j=1,numrots_bv(n-k)
          rot%cosine=cs_bv(n-k,j); rot%sine=ss_bv(n-k,j)
          call general_times_rotation(pq,trp_rot(rot),k1s_bv(n-k,j)-coffs,k2s_bv(n-k,j)-coffs)
       end do
       ! include row k+1 in the decomposition
       b_bv(k+1,lbw1+1:lbw1+n-k)=matmul(b_bv(k+1,lbw1+1:lbw1+n-k),conjg(transpose(pq)))
       ! downdate
       do j=n-k-1,1,-1
          rot=lgivens(pq(j,1),pq(j+1,1))
          call rotation_times_general(trp_rot(rot),pq,j,j+1)
          call tbr_times_rotation(b_bv,n,lbw1,ubw1,n-k-1,rot,coffs+j)
       end do
       do j=2,ubw1
          tmp=get_el_br(b_bv,lbw1,roffs+j,coffs+1)
          call set_el_br(b_bv,lbw1,roffs+j,coffs+1,pq(1,1)*tmp)
       end do
       call up_left_shift(pq)
    end do
    do d=1,lbw+1
       do j=1,n-lbw+d-1
          b_ub(ubw-1+lbw+2-d,j)=b_bv(lbw-d+j+1,d)
       end do
    end do
    do d=lbw+2,lbw+ubw
       do j=1,n-d+lbw+1
          b_ub(ubw-1+lbw+2-d,j+d-lbw-1)=b_bv(j,d)
       end do
    end do
  end subroutine f_c_compress_bv_to_ub_1

end module compress_bv
