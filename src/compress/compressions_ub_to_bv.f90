module compressions_ub_to_bv
use prec
use shift
use rotation
use band_types
use nested_types
use nullvec
implicit none

integer(kind=int32), private, parameter :: nullmaxits=5

interface compress_ub_to_bv_1
   module procedure d_compress_ub_to_bv_1, c_compress_ub_to_bv_1
end interface compress_ub_to_bv_1

interface f_compress_ub_to_bv_1
   module procedure f_d_compress_ub_to_bv_1, f_c_compress_ub_to_bv_1
end interface f_compress_ub_to_bv_1

contains

  ! Rank one compression.

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
  subroutine d_compress_ub_to_bv_1(ub, bv, told, tol, error)
    type(d_ub) :: ub
    type(d_bv) :: bv
    integer(kind=int32), intent(out) :: error
    real(kind=dp), intent(in) :: tol, told
    if (ub%n /= bv%n) then
       error = 3; return
    end if
    call f_d_compress_ub_to_bv_1(ub%b, ub%n, ub%lbw, ub%ubw, ub%lbwmax, ub%ubwmax, ub%numrotsu, ub%j1su, ub%j2su, &
         ub%csu, ub%ssu, bv%b, bv%lbwmax, bv%ubwmax, bv%numrotsv, bv%k1sv, bv%k2sv, bv%csv, bv%ssv, told, tol, &
         error)
    if (error == 0) then
       bv%lbw=ub%lbw; bv%ubw=ub%ubw-1; bv%n=ub%n
    end if
  end subroutine d_compress_ub_to_bv_1

  subroutine f_d_compress_ub_to_bv_1(b_ub, n, lbw, ubw, lbwmax_ub, ubwmax_ub, numrots_ub, &
       j1s_ub, j2s_ub, cs_ub, ss_ub, &
       b_bv, lbwmax_bv, ubwmax_bv, numrots_bv, k1s_bv, k2s_bv, cs_bv, ss_bv, told, tol, error)
    real(kind=dp), dimension(lbwmax_ub+ubwmax_ub+1,n), intent(inout) :: b_ub
    integer(kind=int32), intent(in) :: n, lbw, ubw, lbwmax_ub, ubwmax_ub, lbwmax_bv, ubwmax_bv
    integer(kind=int32), dimension(n), intent(in) :: numrots_ub
    integer(kind=int32), dimension(ubwmax_ub,n), intent(in) :: j1s_ub, j2s_ub
    real(kind=dp), dimension(ubwmax_ub,n), intent(in) :: cs_ub, ss_ub
    real(kind=dp), intent(in) :: tol, told

    real(kind=dp), dimension(n,lbwmax_bv+ubwmax_bv+1), intent(out) :: b_bv
    integer(kind=int32), dimension(n), intent(out) :: numrots_bv
    integer(kind=int32), dimension(n,ubwmax_bv), intent(out) :: k1s_bv, k2s_bv
    real(kind=dp), dimension(n,ubwmax_bv), intent(out) :: cs_bv, ss_bv
    integer(kind=int32), intent(out) :: error

    integer(kind=int32) :: j, k, kk, roffs, coffs, ubw1, lbw1, nullerr, d, minindex
    type(d_rotation) :: rot
    real(kind=dp), target, dimension(ubw+1,ubw+1) :: q
    real(kind=dp), dimension(ubw+1,ubw+1) :: l
    real(kind=dp), dimension(ubw+1) :: x
    real(kind=dp) :: nrma, tmp, mindiag
    real(kind=dp), pointer, dimension(:,:) :: pq


    error = 0
    numrots_bv=0
    ss_bv=0.0_dp; cs_bv=0.0_dp
    k1s_bv=0; k2s_bv=0
    ubw1=ubw+1; lbw1=lbw
    nrma = maxabs(b_ub)*sqrt(real(n))
    
    ! store a copy in case compression fails.
    b_bv(1:n,1:ubw+lbw+1)=transpose(b_ub(1:ubw+lbw+1,1:n))

    if (n < 1) then
       error = 1
       return
    end if
    if (n == 1) then
        b_bv(1,1)=b_ub(1,1)
        return
    end if
    ! must allow for temporary fill-in of one extra superdiagonal.
    if (lbwmax_ub+ubwmax_ub+1<ubw1+lbw1+1) then
       error = 2
       return
    end if
    q=0.0_dp
    do j=1,ubw1
       q(j,j)=1.0_dp
    end do
    ! provide working space in a superdiagonal
    call down_shift(b_ub)

    ! apply u_{n-1}
    do j=1,numrots_ub(n-1)
       rot%cosine=cs_ub(j,n-1); rot%sine=ss_ub(j,n-1)
       call rotation_times_tbc(rot,b_ub,n,lbw1,ubw1,n-1,j1s_ub(j,n-1))
    end do
    ! Apply u_{n-2} to u_{n-ubw+1}
    ! Generate a ubw x ubw lower Hessenberg matrix
    do k=2,ubw-1
       roffs=n-k-ubw1; coffs=n-k
       numrots_bv(k) = k-1
       do j=2,k
          rot=rgivens(get_el_bc(b_ub,ubw1,roffs+j,coffs+j-1),get_el_bc(b_ub,ubw1,roffs+j,coffs+j))
          cs_bv(k,k+1-j)=rot%cosine; ss_bv(k,k+1-j)=rot%sine
          k1s_bv(k,k+1-j)=coffs+j-1; k2s_bv(k,k+1-j)=coffs+j
          call tbc_times_rotation(b_ub,n,lbw1,ubw1,k,rot,coffs+j-1)
          call set_el_bc(b_ub, ubw1, roffs+j,coffs+j,0.0_dp)
       end do
       do j=1,numrots_ub(n-k)
          rot%cosine=cs_ub(j,n-k); rot%sine=ss_ub(j,n-k)
          call rotation_times_tbc(rot,b_ub,n,lbw1,ubw1,n-k,j1s_ub(j,n-k))
       end do
    end do
    ! At this point there is a ubw x ubw lower Hessenberg matrix in
    ! A(n-2*ubw+1:n-ubw, n-ubw+1:n).
    ! Triangularize using Q:
    roffs=n-2*ubw; coffs=n-ubw
    do j=ubw,2,-1
       rot=lgivens2(get_el_bc(b_ub,ubw1,roffs+j-1,coffs+j), get_el_bc(b_ub,ubw1,roffs+j,coffs+j))
       call rotation_times_tbc(trp_rot(rot),b_ub,n,lbw1,ubw1,coffs,roffs+j-1)
       call set_el_bc(b_ub,ubw1,roffs+j-1,coffs+j,0.0_dp)
       call general_times_rotation(q,rot,j,j+1)
    end do
    ! main loop
    do j=ubw,n-ubw-1
       ! There should be a triangular matrix in
       ! B(n-k-ubw+1:n-k,n-k+1:n-k+ubw) representing the L factor
       ! in a QL decomposition.
       !
       ! Try to find a small diagonal of L and move it to the lower right.
       roffs=n-j-ubw; coffs=n-j
       mindiag=abs(get_el_bc(b_ub,ubw1,roffs+1,coffs+1))
       minindex=1
       do k=2,ubw
          !          tmp=abs(get_el_bc(b_ub,ubw1,roffs+k,coffs+k))
          tmp=abs(b_ub(roffs-coffs+ubw1+1,coffs+k))
          if (tmp <= mindiag) then
             minindex=k
             mindiag=tmp
          end if
       end do
       if (mindiag < told*nrma) then
          numrots_bv(j) = ubw - minindex
          do k=minindex,ubw-1
             rot=rgivens(get_el_bc(b_ub,ubw1,roffs+k+1,coffs+k),get_el_bc(b_ub,ubw1,roffs+k+1,coffs+k+1))
             call tbc_times_rotation(b_ub,n,lbw1,ubw1,j,rot,coffs+k)
             cs_bv(j,ubw-k)=rot%cosine; ss_bv(j,ubw-k)=rot%sine
             k1s_bv(j,ubw-k)=coffs+k; k2s_bv(j,ubw-k)=coffs+k+1
             call set_el_bc(b_ub,ubw1,roffs+k+1,coffs+k+1,0.0_dp)
             ! swap rows in l
             do kk=1,k+1
                tmp=get_el_bc(b_ub,ubw1,roffs+k+1,coffs+kk)
                call set_el_bc(b_ub,ubw1,roffs+k+1,coffs+kk, &
                     get_el_bc(b_ub,ubw1,roffs+k,coffs+kk))
                call set_el_bc(b_ub,ubw1,roffs+k,coffs+kk,tmp)
             end do
             do kk=1,ubw1
                tmp=q(kk,k+1)
                q(kk,k+1)=q(kk,k+2)
                q(kk,k+2)=tmp
             end do
          end do
       else
          ! Need a more general null vector of L.
          call submatrix_bc(b_ub,lbw1,ubw1, roffs+1,roffs+ubw,coffs+1,coffs+ubw,l(1:ubw,1:ubw))
          call f_d_lower_right_nullvec(x(1:ubw),l(1:ubw,1:ubw),tol*nrma,nullmaxits, nullerr)
          if (nullerr >= 0) then ! nullvector found
             ! Introduce a zero while preserving the QL factorization.
             numrots_bv(j)=ubw-1
             do k=1,ubw-1
                rot=lgivens2(x(k),x(k+1))
                call rotation_times_general(trp_rot(rot),x,k,k+1)
                call tbc_times_rotation(b_ub,n,lbw1,ubw1,j,rot,coffs+k)
                cs_bv(j,ubw-k)=rot%cosine; ss_bv(j,ubw-k)=rot%sine
                k1s_bv(j,ubw-k)=coffs+k; k2s_bv(j,ubw-k)=coffs+k+1
                rot=lgivens2(get_el_bc(b_ub,ubw1,roffs+k,coffs+k+1),get_el_bc(b_ub,ubw1,roffs+k+1,coffs+k+1))
                call rotation_times_tbc(trp_rot(rot),b_ub,n,lbw1,ubw1,coffs,roffs+k)
                call general_times_rotation(q,rot,k+1,k+2)
                call set_el_bc(b_ub,ubw1,roffs+k,coffs+k+1,0.0_dp)
             end do
             call set_el_bc(b_ub,ubw1,roffs+ubw,coffs+ubw,0.0_dp)
          else
             ! Compression has failed.  restore b_ub and return.
             error=4
             b_ub(1:ubw+lbw+1,1:n)=transpose(b_bv(1:n,1:ubw+lbw+1))
             return
          end if
       end if
       ! include an extra row in L in anticipation of including column n-j
       roffs=roffs-1
       ! Apply u_{n-j}: This is the biggest cost.
       do k=1,numrots_ub(n-j)
          rot%cosine=cs_ub(k,n-j); rot%sine=ss_ub(k,n-j)
          call rotation_times_general(rot,q,j1s_ub(k,n-j)-roffs,j2s_ub(k,n-j)-roffs)
       end do
       ! include column n-j in the LQ factorization
       b_ub(2:ubw+2,n-j)=matmul(transpose(q),b_ub(2:ubw+2,n-j))
       coffs=coffs-1
       ! downdate
       do k=1,ubw
          rot=rgivens2(q(ubw1,k),q(ubw1,k+1))
          call general_times_rotation(q,rot,k,k+1)
          q(ubw1,k)=0.0_dp
          call rotation_times_tbc(trp_rot(rot),b_ub,n,lbw1,ubw1,coffs,roffs+k)
       end do
       do k=1,ubw
          ! tmp=get_el_bc(b_ub,ubw1,n-j,n-j+k-1)
          ! call set_el_bc(b_ub,ubw1,n-j,n-j+k-1,tmp*q(ubw1,ubw1))
          b_ub(n-j-(n-j+k-1)+ubw1+1,n-j+k-1) = q(ubw1,ubw1)*b_ub(n-j-(n-j+k-1)+ubw1+1,n-j+k-1)
       end do
       call down_right_shift(q)
       q(1,1)=1.0_dp
       roffs=roffs-1
       ! triangularize a Hessenberg matrix.
       do k=ubw,2,-1
          rot=lgivens2(get_el_bc(b_ub,ubw1,roffs+k,coffs+k),get_el_bc(b_ub,ubw1,roffs+k+1,coffs+k))
          call rotation_times_tbc(trp_rot(rot),b_ub,n,lbw1,ubw1,coffs,roffs+k)
          call general_times_rotation(q,rot,k,k+1)
          call set_el_bc(b_ub,ubw1,roffs+k,coffs+k,0.0_dp)
       end do
    end do
    !
    ! Intermediate step: L is still square, as is q.  We need to find a null vector.
    ! corresponds to j=n-ubw, roffs=0, coffs=ubw.  At the end of this section
    ! L will be trapezoidal.
    !
    call up_left_shift(q)
    pq=>q(1:ubw,1:ubw)
    call submatrix_bc(b_ub,lbw1,ubw1, 1,ubw,ubw+1,2*ubw,l(1:ubw,1:ubw))
    call f_d_lower_right_nullvec(x(1:ubw),l(1:ubw,1:ubw),tol*nrma,nullmaxits, nullerr)
    if (nullerr >= 0) then ! nullvector found
       ! Introduce a zero while preserving the QL factorization.
       numrots_bv(n-ubw)=ubw-1
       do k=1,ubw-1
          rot=lgivens2(x(k),x(k+1))
          call rotation_times_general(trp_rot(rot),x,k,k+1)
          call tbc_times_rotation(b_ub,n,lbw1,ubw1,n-ubw,rot,ubw+k)
          cs_bv(n-ubw,ubw-k)=rot%cosine; ss_bv(n-ubw,ubw-k)=rot%sine
          k1s_bv(n-ubw,ubw-k)=ubw+k; k2s_bv(n-ubw,ubw-k)=ubw+k+1
          rot=lgivens2(get_el_bc(b_ub,ubw1,k,ubw+k+1),get_el_bc(b_ub,ubw1,k+1,ubw+k+1))
          call rotation_times_tbc(trp_rot(rot),b_ub,n,lbw1,ubw1,ubw,k)
          call general_times_rotation(pq,rot,k,k+1)
          call set_el_bc(b_ub,ubw1,k,ubw+k+1,0.0_dp)
       end do
       call set_el_bc(b_ub,ubw1,ubw,2*ubw,0.0_dp)
       roffs=roffs-1
       ! Apply u_{ubw}
       do k=1,numrots_ub(ubw)
          rot%cosine=cs_ub(k,ubw); rot%sine=ss_ub(k,ubw)
          call rotation_times_general(rot,pq,j1s_ub(k,ubw),j2s_ub(k,n-j))
       end do
       ! include column n-j in the LQ factorization
       b_ub(3:ubw+2,ubw)=matmul(transpose(pq),b_ub(3:ubw+2,ubw))
       !restore triangularity
       do k=ubw,2,-1
          rot=lgivens2(get_el_bc(b_ub,ubw1,k-1,k+ubw-1), get_el_bc(b_ub,ubw1,k,k+ubw-1))
          call rotation_times_tbc(trp_rot(rot),b_ub,n,lbw1,ubw1,ubw-1,k-1)
          call general_times_rotation(pq,rot,k-1,k)
       end do
       ! downdate
       do k=1,ubw-1
          rot=rgivens2(pq(ubw,k),q(ubw,k+1))
          call general_times_rotation(pq,rot,k,k+1)
          pq(ubw,k)=0.0_dp
          call rotation_times_tbc(trp_rot(rot),b_ub,n,lbw1,ubw1,ubw-1,k)
       end do
       do k=1,ubw
          tmp=get_el_bc(b_ub,ubw1,ubw,ubw+k-1)
          call set_el_bc(b_ub,ubw1,ubw,ubw+k-1,tmp*pq(ubw,ubw))
       end do
    else
       ! Compression has failed.  restore b_ub and return.
       error=4
       b_ub(1:ubw+lbw+1,1:n)=transpose(b_bv(1:n,1:ubw+lbw+1))
       return
    end if
    ! Final Stage:
    ! Now L is a lower trapezoidal matrix of size (ubw-1) x ubw
    do j=n-ubw+1,n-1
       coffs=n-j
       pq => q(1:coffs,1:coffs)
       ! apply v_j
       numrots_bv(j)=coffs
       do k=1,coffs
          rot=rgivens(get_el_bc(b_ub,ubw1,k,k+ubw-1),get_el_bc(b_ub,ubw1,k,k+ubw))
          call tbc_times_rotation(b_ub,n,lbw1,ubw1,j,rot,ubw+k-1)
          cs_bv(j,coffs-k+1)=rot%cosine; ss_bv(j,coffs-k+1)=rot%sine
          k1s_bv(j,coffs-k+1)=ubw+k-1; k2s_bv(j,coffs-k+1)=ubw+k
       end do
       ! apply u_{n-j}
       do k=1,numrots_ub(n-j)
          rot%cosine=cs_ub(k,n-j); rot%sine=ss_ub(k,n-j)
          call rotation_times_general(rot,pq,j1s_ub(k,n-j),j2s_ub(k,n-j))
       end do
       ! Include column n-j in the QL decomposition
       b_ub(ubw+2-(n-j)+1:ubw+2-(n-j)+coffs,n-j)=matmul(transpose(pq),b_ub(ubw+2-(n-j)+1:ubw+2-(n-j)+coffs,n-j))
       ! downdate
       do k=1,coffs-1
          rot=rgivens2(pq(coffs,k),pq(coffs,k+1))
          call general_times_rotation(pq,rot,k,k+1)
          call rotation_times_tbc(trp_rot(rot),b_ub,n,lbw1,ubw1,coffs-1,k)
       end do
       ! TODO: double check this.
       do k=1,ubw
          tmp=get_el_bc(b_ub,ubw1,coffs,coffs+k-1)
          call set_el_bc(b_ub,ubw1,coffs,coffs+k-1,pq(coffs,coffs)*tmp)
       end do
    end do
    call up_shift(b_ub)
    call up_shift(b_ub)
    do d=1,lbw+1
       do j=1,n-lbw+d-1
          b_bv(lbw-d+j+1,d)=b_ub(ubw-1+lbw+2-d,j)
       end do
    end do
    do d=lbw+2,lbw+ubw
       do j=1,n-d+lbw+1
          b_bv(j,d)=b_ub(ubw+lbw+1-d,j+d-lbw-1)
       end do
    end do
  end subroutine f_d_compress_ub_to_bv_1

  subroutine c_compress_ub_to_bv_1(ub, bv, told, tol, error)
    type(c_ub) :: ub
    type(c_bv) :: bv
    integer(kind=int32), intent(out) :: error
    real(kind=dp), intent(in) :: tol, told
    if (ub%n /= bv%n) then
       error = 3; return
    end if
    call f_c_compress_ub_to_bv_1(ub%b, ub%n, ub%lbw, ub%ubw, ub%lbwmax, ub%ubwmax, ub%numrotsu, ub%j1su, ub%j2su, &
         ub%csu, ub%ssu, bv%b, bv%lbwmax, bv%ubwmax, bv%numrotsv, bv%k1sv, bv%k2sv, bv%csv, bv%ssv, told, tol, &
         error)
    if (error == 0) then
       bv%lbw=ub%lbw; bv%ubw=ub%ubw-1; bv%n=ub%n
    end if
  end subroutine c_compress_ub_to_bv_1

  subroutine f_c_compress_ub_to_bv_1(b_ub, n, lbw, ubw, lbwmax_ub, ubwmax_ub, numrots_ub, &
       j1s_ub, j2s_ub, cs_ub, ss_ub, &
       b_bv, lbwmax_bv, ubwmax_bv, numrots_bv, k1s_bv, k2s_bv, cs_bv, ss_bv, told, tol, error)
    complex(kind=dp), dimension(lbwmax_ub+ubwmax_ub+1,n), intent(inout) :: b_ub
    integer(kind=int32), intent(in) :: n, lbw, ubw, lbwmax_ub, ubwmax_ub, lbwmax_bv, ubwmax_bv
    integer(kind=int32), dimension(n), intent(in) :: numrots_ub
    integer(kind=int32), dimension(ubwmax_ub,n), intent(in) :: j1s_ub, j2s_ub
    complex(kind=dp), dimension(ubwmax_ub,n), intent(in) :: cs_ub, ss_ub
    real(kind=dp), intent(in) :: tol, told

    complex(kind=dp), dimension(n,lbwmax_bv+ubwmax_bv+1), intent(out) :: b_bv
    integer(kind=int32), dimension(n), intent(out) :: numrots_bv
    integer(kind=int32), dimension(n,ubwmax_bv), intent(out) :: k1s_bv, k2s_bv
    complex(kind=dp), dimension(n,ubwmax_bv), intent(out) :: cs_bv, ss_bv
    integer(kind=int32), intent(out) :: error

    integer(kind=int32) :: j, k, kk, roffs, coffs, ubw1, lbw1, nullerr, d, minindex
    type(c_rotation) :: rot
    complex(kind=dp), target, dimension(ubw+1,ubw+1) :: q
    complex(kind=dp), dimension(ubw+1,ubw+1) :: l
    complex(kind=dp), dimension(ubw+1) :: x
    real(kind=dp) :: nrma, tmp, mindiag
    complex(kind=dp) :: tmpc
    complex(kind=dp), pointer, dimension(:,:) :: pq


    error = 0
    numrots_bv=0
    ss_bv=(0.0_dp,0.0_dp); cs_bv=(0.0_dp,0.0_dp)
    k1s_bv=0; k2s_bv=0
    ubw1=ubw+1; lbw1=lbw
    nrma = maxabs(b_ub)*sqrt(real(n))
    
    ! store a copy in case compression fails.
    b_bv(1:n,1:ubw+lbw+1)=transpose(b_ub(1:ubw+lbw+1,1:n))

    if (n < 1) then
       error = 1
       return
    end if
    if (n == 1) then
        b_bv(1,1)=b_ub(1,1)
        return
    end if
    ! must allow for temporary fill-in of one extra superdiagonal.
    if (lbwmax_ub+ubwmax_ub+1<ubw1+lbw1+1) then
       error = 2
       return
    end if
    q=(0.0_dp,0.0_dp)
    do j=1,ubw1
       q(j,j)=(1.0_dp,0.0_dp)
    end do
    ! provide working space in a superdiagonal
    call down_shift(b_ub)

    ! apply u_{n-1}
    do j=1,numrots_ub(n-1)
       rot%cosine=cs_ub(j,n-1); rot%sine=ss_ub(j,n-1)
       call rotation_times_tbc(rot,b_ub,n,lbw1,ubw1,n-1,j1s_ub(j,n-1))
    end do
    ! Apply u_{n-2} to u_{n-ubw+1}
    ! Generate a ubw x ubw lower Hessenberg matrix
    do k=2,ubw-1
       roffs=n-k-ubw1; coffs=n-k
       numrots_bv(k) = k-1
       do j=2,k
          rot=rgivens(get_el_bc(b_ub,ubw1,roffs+j,coffs+j-1),get_el_bc(b_ub,ubw1,roffs+j,coffs+j))
          cs_bv(k,k+1-j)=rot%cosine; ss_bv(k,k+1-j)=rot%sine
          k1s_bv(k,k+1-j)=coffs+j-1; k2s_bv(k,k+1-j)=coffs+j
          call tbc_times_rotation(b_ub,n,lbw1,ubw1,k,rot,coffs+j-1)
          call set_el_bc(b_ub, ubw1, roffs+j,coffs+j,(0.0_dp,0.0_dp))
       end do
       do j=1,numrots_ub(n-k)
          rot%cosine=cs_ub(j,n-k); rot%sine=ss_ub(j,n-k)
          call rotation_times_tbc(rot,b_ub,n,lbw1,ubw1,n-k,j1s_ub(j,n-k))
       end do
    end do
    ! At this point there is a ubw x ubw lower Hessenberg matrix in
    ! A(n-2*ubw+1:n-ubw, n-ubw+1:n).
    ! Triangularize using Q:
    roffs=n-2*ubw; coffs=n-ubw
    do j=ubw,2,-1
       rot=lgivens2(get_el_bc(b_ub,ubw1,roffs+j-1,coffs+j), get_el_bc(b_ub,ubw1,roffs+j,coffs+j))
       call rotation_times_tbc(trp_rot(rot),b_ub,n,lbw1,ubw1,coffs,roffs+j-1)
       call set_el_bc(b_ub,ubw1,roffs+j-1,coffs+j,(0.0_dp,0.0_dp))
       call general_times_rotation(q,rot,j,j+1)
    end do
    ! main loop
    do j=ubw,n-ubw-1
       ! There should be a triangular matrix in
       ! B(n-k-ubw+1:n-k,n-k+1:n-k+ubw) representing the L factor
       ! in a QL decomposition.
       !
       ! Try to find a small diagonal of L.
       roffs=n-j-ubw; coffs=n-j
       mindiag=abs(get_el_bc(b_ub,ubw1,roffs+1,coffs+1))
       minindex=1
       do k=2,ubw
          !          tmp=abs(get_el_bc(b_ub,ubw1,roffs+k,coffs+k))
          tmp=abs(b_ub(roffs-coffs+ubw1+1,coffs+k))
          if (tmp <= mindiag) then
             minindex=k
             mindiag=tmp
          end if
       end do
       if (mindiag < told*nrma) then
          numrots_bv(j) = ubw - minindex
          do k=minindex,ubw-1
             rot=rgivens(get_el_bc(b_ub,ubw1,roffs+k+1,coffs+k),get_el_bc(b_ub,ubw1,roffs+k+1,coffs+k+1))
             call tbc_times_rotation(b_ub,n,lbw1,ubw1,j,rot,coffs+k)
             cs_bv(j,ubw-k)=rot%cosine; ss_bv(j,ubw-k)=rot%sine
             k1s_bv(j,ubw-k)=coffs+k; k2s_bv(j,ubw-k)=coffs+k+1
             call set_el_bc(b_ub,ubw1,roffs+k+1,coffs+k+1,(0.0_dp,0.0_dp))
             ! swap rows in l
             do kk=1,k+1
                tmpc=get_el_bc(b_ub,ubw1,roffs+k+1,coffs+kk)
                call set_el_bc(b_ub,ubw1,roffs+k+1,coffs+kk, &
                     get_el_bc(b_ub,ubw1,roffs+k,coffs+kk))
                call set_el_bc(b_ub,ubw1,roffs+k,coffs+kk,tmpc)
             end do
             do kk=1,ubw1
                tmpc=q(kk,k+1)
                q(kk,k+1)=q(kk,k+2)
                q(kk,k+2)=tmpc
             end do
          end do
       else
          ! Need a more general null vector of L.
          call submatrix_bc(b_ub,lbw1,ubw1, roffs+1,roffs+ubw,coffs+1,coffs+ubw,l(1:ubw,1:ubw))
          call f_c_lower_right_nullvec(x(1:ubw),l(1:ubw,1:ubw),tol*nrma,nullmaxits, nullerr)
          if (nullerr >= 0) then ! nullvector found
             ! Introduce a zero while preserving the QL factorization.
             numrots_bv(j)=ubw-1
             do k=1,ubw-1
                rot=lgivens2(x(k),x(k+1))
                call rotation_times_general(trp_rot(rot),x,k,k+1)
                call tbc_times_rotation(b_ub,n,lbw1,ubw1,j,rot,coffs+k)
                cs_bv(j,ubw-k)=rot%cosine; ss_bv(j,ubw-k)=rot%sine
                k1s_bv(j,ubw-k)=coffs+k; k2s_bv(j,ubw-k)=coffs+k+1
                rot=lgivens2(get_el_bc(b_ub,ubw1,roffs+k,coffs+k+1),get_el_bc(b_ub,ubw1,roffs+k+1,coffs+k+1))
                call rotation_times_tbc(trp_rot(rot),b_ub,n,lbw1,ubw1,coffs,roffs+k)
                call general_times_rotation(q,rot,k+1,k+2)
                call set_el_bc(b_ub,ubw1,roffs+k,coffs+k+1,(0.0_dp,0.0_dp))
             end do
             call set_el_bc(b_ub,ubw1,roffs+ubw,coffs+ubw,(0.0_dp,0.0_dp))
          else
             ! Compression has failed.  restore b_ub and return.
             error=4
             b_ub(1:ubw+lbw+1,1:n)=transpose(b_bv(1:n,1:ubw+lbw+1))
             return
          end if
       end if
       ! include an extra row in L in anticipation of including column n-j
       roffs=roffs-1
       ! Apply u_{n-j}
       do k=1,numrots_ub(n-j)
          rot%cosine=cs_ub(k,n-j); rot%sine=ss_ub(k,n-j)
          call rotation_times_general(rot,q,j1s_ub(k,n-j)-roffs,j2s_ub(k,n-j)-roffs)
       end do
       ! include column n-j in the LQ factorization
       b_ub(2:ubw+2,n-j)=matmul(conjg(transpose(q)),b_ub(2:ubw+2,n-j))
       coffs=coffs-1
       ! downdate
       do k=1,ubw
          rot=rgivens2(q(ubw1,k),q(ubw1,k+1))
          call general_times_rotation(q,rot,k,k+1)
          q(ubw1,k)=(0.0_dp,0.0_dp)
          call rotation_times_tbc(trp_rot(rot),b_ub,n,lbw1,ubw1,coffs,roffs+k)
       end do
       do k=1,ubw
          ! tmp=get_el_bc(b_ub,ubw1,n-j,n-j+k-1)
          ! call set_el_bc(b_ub,ubw1,n-j,n-j+k-1,tmp*q(ubw1,ubw1))
          b_ub(n-j-(n-j+k-1)+ubw1+1,n-j+k-1) = q(ubw1,ubw1)*b_ub(n-j-(n-j+k-1)+ubw1+1,n-j+k-1)
       end do
       call down_right_shift(q)
       q(1,1)=(1.0_dp,0.0_dp)
       roffs=roffs-1
       ! triangularize a Hessenberg matrix.
       do k=ubw,2,-1
          rot=lgivens2(get_el_bc(b_ub,ubw1,roffs+k,coffs+k),get_el_bc(b_ub,ubw1,roffs+k+1,coffs+k))
          call rotation_times_tbc(trp_rot(rot),b_ub,n,lbw1,ubw1,coffs,roffs+k)
          call general_times_rotation(q,rot,k,k+1)
          call set_el_bc(b_ub,ubw1,roffs+k,coffs+k,(0.0_dp,0.0_dp))
       end do
    end do
    !
    ! Intermediate step: L is still square, as is q.  We need to find a null vector.
    ! corresponds to j=n-ubw, roffs=0, coffs=ubw.  At the end of this section
    ! L will be trapezoidal.
    !
    call up_left_shift(q)
    pq=>q(1:ubw,1:ubw)
    call submatrix_bc(b_ub,lbw1,ubw1, 1,ubw,ubw+1,2*ubw,l(1:ubw,1:ubw))
    call f_c_lower_right_nullvec(x(1:ubw),l(1:ubw,1:ubw),tol*nrma,nullmaxits, nullerr)
    if (nullerr >= 0) then ! nullvector found
       ! Introduce a zero while preserving the QL factorization.
       numrots_bv(j)=ubw-1
       do k=1,ubw-1
          rot=lgivens2(x(k),x(k+1))
          call rotation_times_general(trp_rot(rot),x,k,k+1)
          call tbc_times_rotation(b_ub,n,lbw1,ubw1,n-ubw,rot,ubw+k)
          cs_bv(n-ubw,ubw-k)=rot%cosine; ss_bv(n-ubw,ubw-k)=rot%sine
          k1s_bv(n-ubw,ubw-k)=ubw+k; k2s_bv(n-ubw,ubw-k)=ubw+k+1
          rot=lgivens2(get_el_bc(b_ub,ubw1,k,ubw+k+1),get_el_bc(b_ub,ubw1,k+1,ubw+k+1))
          call rotation_times_tbc(trp_rot(rot),b_ub,n,lbw1,ubw1,ubw,k)
          call general_times_rotation(pq,rot,k,k+1)
          call set_el_bc(b_ub,ubw1,k,ubw+k+1,(0.0_dp,0.0_dp))
       end do
       call set_el_bc(b_ub,ubw1,ubw,2*ubw,(0.0_dp,0.0_dp))
       roffs=roffs-1
       ! Apply u_{ubw}
       do k=1,numrots_ub(ubw)
          rot%cosine=cs_ub(k,ubw); rot%sine=ss_ub(k,ubw)
          call rotation_times_general(rot,pq,j1s_ub(k,ubw),j2s_ub(k,n-j))
       end do
       ! include column n-j in the LQ factorization
       b_ub(3:ubw+2,ubw)=matmul(conjg(transpose(pq)),b_ub(3:ubw+2,ubw))
       !restore triangularity
       do k=ubw,2,-1
          rot=lgivens2(get_el_bc(b_ub,ubw1,k-1,k+ubw-1), get_el_bc(b_ub,ubw1,k,k+ubw-1))
          call rotation_times_tbc(trp_rot(rot),b_ub,n,lbw1,ubw1,ubw-1,k-1)
          call general_times_rotation(pq,rot,k-1,k)
       end do
       ! downdate
       do k=1,ubw-1
          rot=rgivens2(pq(ubw,k),q(ubw,k+1))
          call general_times_rotation(pq,rot,k,k+1)
          pq(ubw,k)=(0.0_dp,0.0_dp)
          call rotation_times_tbc(trp_rot(rot),b_ub,n,lbw1,ubw1,ubw-1,k)
       end do
       do k=1,ubw
          tmpc=get_el_bc(b_ub,ubw1,ubw,ubw+k-1)
          call set_el_bc(b_ub,ubw1,ubw,ubw+k-1,tmpc*pq(ubw,ubw))
       end do
    else
       ! Compression has failed.  restore b_ub and return.
       error=4
       b_ub(1:ubw+lbw+1,1:n)=transpose(b_bv(1:n,1:ubw+lbw+1))
       return
    end if
    ! Final Stage:
    ! Now L is a lower trapezoidal matrix of size (ubw-1) x ubw
    do j=n-ubw+1,n-1
       coffs=n-j
       pq => q(1:coffs,1:coffs)
       ! apply v_j
       numrots_bv(j)=coffs
       do k=1,coffs
          rot=rgivens(get_el_bc(b_ub,ubw1,k,k+ubw-1),get_el_bc(b_ub,ubw1,k,k+ubw))
          call tbc_times_rotation(b_ub,n,lbw1,ubw1,j,rot,ubw+k-1)
          cs_bv(j,coffs-k+1)=rot%cosine; ss_bv(j,coffs-k+1)=rot%sine
          k1s_bv(j,coffs-k+1)=ubw+k-1; k2s_bv(j,coffs-k+1)=ubw+k
       end do
       ! apply u_{n-j}
       do k=1,numrots_ub(n-j)
          rot%cosine=cs_ub(k,n-j); rot%sine=ss_ub(k,n-j)
          call rotation_times_general(rot,pq,j1s_ub(k,n-j),j2s_ub(k,n-j))
       end do
       ! Include column n-j in the QL decomposition
       b_ub(ubw+2-(n-j)+1:ubw+2-(n-j)+coffs,n-j)=matmul(conjg(transpose(pq)), &
            b_ub(ubw+2-(n-j)+1:ubw+2-(n-j)+coffs,n-j))
       ! downdate
       do k=1,coffs-1
          rot=rgivens2(pq(coffs,k),pq(coffs,k+1))
          call general_times_rotation(pq,rot,k,k+1)
          call rotation_times_tbc(trp_rot(rot),b_ub,n,lbw1,ubw1,coffs-1,k)
       end do
       do k=1,ubw
          tmpc=get_el_bc(b_ub,ubw1,coffs,coffs+k-1)
          call set_el_bc(b_ub,ubw1,coffs,coffs+k-1,pq(coffs,coffs)*tmpc)
       end do
    end do
    call up_shift(b_ub)
    call up_shift(b_ub)
    do d=1,lbw+1
       do j=1,n-lbw+d-1
          b_bv(lbw-d+j+1,d)=b_ub(ubw-1+lbw+2-d,j)
       end do
    end do
    do d=lbw+2,lbw+ubw
       do j=1,n-d+lbw+1
          b_bv(j,d)=b_ub(ubw+lbw+1-d,j+d-lbw-1)
       end do
    end do
  end subroutine f_c_compress_ub_to_bv_1


end module compressions_ub_to_bv
