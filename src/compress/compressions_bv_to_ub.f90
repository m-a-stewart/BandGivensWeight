module compressions_bv_to_ub
use prec
use shift
use rotation
use band_types
use nested_types
use nullvec
implicit none

integer(kind=int32), private, parameter :: nullmaxits=5

interface compress_bv_to_ub
   module procedure d_compress_bv_to_ub, c_compress_bv_to_ub
end interface compress_bv_to_ub

interface f_compress_bv_to_ub
   module procedure f_d_compress_bv_to_ub, f_c_compress_bv_to_ub
end interface f_compress_bv_to_ub

contains

!
! Errors:
! 0: no error
! 1: n<1
! 3: ub%n /= bv%n
! 4: compression error, no compression achieved
!    or less compression within the specified tolderance than requested.
! 
! told governs whether a diagonal of L is considered small enough to move to the
! lower right corner of L.  If this tolerance is not met, the algorithm looks
! for a better null vector.  If both are zero the result is a forced compression
! that drops the upper bandwidth by dr.  Forced compression always computes
! a null vector.

  subroutine d_compress_bv_to_ub(bv, ub, told, tol, dr, error)
    type(d_ub) :: ub
    type(d_bv) :: bv
    integer(kind=int32), intent(out) :: error
    integer(kind=int32), intent(in)  :: dr
    real(kind=dp), intent(in) :: tol, told
    if (ub%n /= bv%n) then
       error = 3; return
    end if
    call f_d_compress_bv_to_ub(bv%b, bv%n, bv%lbw, bv%ubw, bv%lbwmax, bv%ubwmax, bv%numrotsv, bv%k1sv, bv%k2sv, &
         bv%csv, bv%ssv, ub%b, ub%ubw, ub%lbwmax, ub%ubwmax, ub%numrotsu, ub%j1su, ub%j2su, ub%csu, ub%ssu, &
         told, tol, dr, error)
    if (error == 0) then
       ub%lbw=bv%lbw; ub%n=bv%n
    end if
  end subroutine d_compress_bv_to_ub

  subroutine f_d_compress_bv_to_ub(b_bv, n, lbw, ubw, lbwmax_bv, ubwmax_bv, numrots_bv, &
       k1s_bv, k2s_bv, cs_bv, ss_bv, &
       b_ub, ubw_ub, lbwmax_ub, ubwmax_ub, numrots_ub, j1s_ub, j2s_ub, cs_ub, ss_ub, told, tol, dr, error)
    integer(kind=int32), intent(in) :: n, lbw, ubw, lbwmax_ub, ubwmax_ub, lbwmax_bv, ubwmax_bv, dr
    real(kind=dp), dimension(n,lbwmax_bv+ubwmax_bv+1), intent(inout) :: b_bv
    integer(kind=int32), dimension(n), intent(in) :: numrots_bv
    integer(kind=int32), dimension(n,ubwmax_bv), intent(in) :: k1s_bv, k2s_bv
    real(kind=dp), dimension(n,ubwmax_bv), intent(in) :: cs_bv, ss_bv

    real(kind=dp), dimension(lbwmax_ub+ubwmax_ub+1,n), intent(out) :: b_ub
    integer(kind=int32), dimension(n), intent(out) :: numrots_ub
    integer(kind=int32), dimension(ubwmax_ub,n), intent(out) :: j1s_ub, j2s_ub
    real(kind=dp), dimension(ubwmax_ub,n), intent(out) :: cs_ub, ss_ub

    real(kind=dp), intent(in) :: tol, told
    integer(kind=int32), intent(out) :: error, ubw_ub

    integer(kind=int32) :: j, k, jj, roffs, coffs, ubw2, nullerr, d, minindex, ml, nl, dnl, nq
    type(d_rotation) :: rot
    real(kind=dp), target, dimension(ubw+1,ubw+1) :: q
    real(kind=dp), dimension(ubw+1,ubw+1) :: l
    real(kind=dp), dimension(ubw+1) :: x
    real(kind=dp) :: nrma, tmp, mindiag
    integer(kind=int32), dimension(n) :: ubws

    error = 0
    numrots_ub=0
    ss_ub=0.0_dp; cs_ub=0.0_dp
    j1s_ub=0; j2s_ub=0
    b_ub=0.0_dp
    ubw2=ubw+2
    nrma = maxabs(b_bv)*sqrt(real(n))
    ubws=0

    if (n < 1) then
       error = 1
       return
    end if
    if (n == 1) then
        b_ub(1,1)=b_bv(1,1)
        return
    end if
    ! must allow for temporary fill-in of one extra superdiagonal in b_bv.
    if (lbwmax_bv+ubwmax_bv+1<ubw2+lbw+1) then
       error = 2
       return
    end if
    nl=1
    ml=1
    q=0.0_dp
    do j=1,ubw+1
       q(j,j)=1.0_dp
    end do
    ! provide working space in a superdiagonal on the right of b_bv
    b_bv(:,ubw+lbw+2:ubw+lbw+3)=0.0_dp
    ! Initial LQ factorization.
    roffs=0; coffs=1
    do j=ubw-1,1,-1
       rot=rgivens(get_el_br(b_bv,lbw,1,coffs+j),get_el_br(b_bv,lbw,1,coffs+j+1))
       call tbr_times_rotation(b_bv,n,lbw,ubw2,n-1,rot,coffs+j)
       call rotation_times_general(trp_rot(rot), q, j,j+1)
       call set_el_br(b_bv,lbw,1,coffs+j+1,0.0_dp)
    end do
    ! Apply v_{n-1} to q.
    do j=1,numrots_bv(n-1)
       rot%cosine=cs_bv(n-1,j); rot%sine=ss_bv(n-1,j)
       call general_times_rotation(q,trp_rot(rot),k1s_bv(n-1,j)-1,k2s_bv(n-1,j)-1)
    end do
    ! main loop
    kloop: do k=1,n
       ! Current, possibly singular, L should be contained in
       ! b_bv(k-nl+1:k,k+1:k+nl)
       nq=min(ubw+1,n-k)
       roffs=k-nl; coffs=k
       mindiag=abs(get_el_br(b_bv,lbw,roffs+nl,coffs+nl))
       minindex=nl
       do j=nl-1,1,-1
          tmp=abs(get_el_br(b_bv,lbw,roffs+j,coffs+j))
          if (tmp <= mindiag) then
             minindex=j
             mindiag=tmp
          end if
       end do
       if (mindiag <= told*nrma) then
          dnl=0  ! don't increase nl
          ubws(k)=nl-1
          numrots_ub(k) = minindex-1
          do j=minindex,2,-1
             rot=lgivens2(get_el_br(b_bv,lbw,roffs+j-1,coffs+j-1),get_el_br(b_bv,lbw,roffs+j,coffs+j-1))
             call rotation_times_tbr(trp_rot(rot), b_bv,n,lbw,ubw2,coffs,roffs+j-1)
             cs_ub(j-1,k)=rot%cosine; ss_ub(j-1,k)=rot%sine
             j1s_ub(j-1,k)=roffs+j-1; j2s_ub(j-1,k)=roffs+j
             call set_el_br(b_bv,lbw,roffs+j-1,coffs+j-1,0.0_dp)
             ! swap columns in L
             do jj=j-1,nl
                tmp=get_el_br(b_bv,lbw,roffs+jj,coffs+j)
                call set_el_br(b_bv,lbw,roffs+jj,coffs+j, &
                     get_el_br(b_bv,lbw,roffs+jj,coffs+j-1))
                call set_el_br(b_bv,lbw,roffs+jj,coffs+j-1,tmp)
             end do
             do jj=1,ubw+1
                tmp=q(j,jj)
                q(j,jj)=q(j-1,jj)
                q(j-1,jj)=tmp
             end do
          end do
          call set_el_br(b_bv,lbw,roffs+1,coffs+1,0.0_dp)
       else ! find a null vector
          call submatrix_br(b_bv,lbw,ubw2, roffs+1,roffs+nl,coffs+1,coffs+nl,l(1:nl,1:nl))
          call f_d_lower_left_nullvec(x(1:nl),l(1:nl,1:nl),tol*nrma,nullmaxits, nullerr)
          if ((nullerr >= 0 .and. tol > 0.0_dp) .or. &
               (tol==0.0_dp .and. told == 0.0_dp .and. nl > ubw-dr)) then ! null vector found
             dnl=0
             ubws(k)=nl-1
             ! Introduce a zero while preserving the triangularity of L
             numrots_ub(k)=nl-1
             do j=nl-1,1,-1
                rot=rgivens(x(j),x(j+1))
                call general_times_rotation(x,rot,j,j+1)
                call rotation_times_tbr(trp_rot(rot),b_bv,n,lbw,ubw2,coffs,roffs+j)
                cs_ub(j,k)=rot%cosine; ss_ub(j,k)=rot%sine
                j1s_ub(j,k)=roffs+j; j2s_ub(j,k)=roffs+j+1
                rot=rgivens(get_el_br(b_bv,lbw,roffs+j,coffs+j),get_el_br(b_bv,lbw,roffs+j,coffs+j+1))
                call tbr_times_rotation(b_bv,n,lbw,ubw2,n-k,rot,coffs+j)
                call rotation_times_general(trp_rot(rot),q,j,j+1)
                call set_el_br(b_bv,lbw,roffs+j,coffs+j+1,0.0_dp)
             end do
             call set_el_br(b_bv,lbw,roffs+1,coffs+1,0.0_dp)
          else
             ! Compression has failed; increase nl.
             ubws(k)=nl
             dnl=1
          end if
       end if
       !
       ! Downdate
       ! 
       do j=nq-1,1,-1
          rot=lgivens(q(j,1),q(j+1,1))
          call rotation_times_general(trp_rot(rot),q,j,j+1)
          q(j+1,1)=0.0_dp
          call tbr_times_rotation(b_bv,n,lbw,ubw2,n-k,rot,coffs+j)
       end do
       do j=1,nl
          tmp=get_el_br(b_bv,lbw,roffs+j,coffs+1)
          call set_el_br(b_bv,lbw,roffs+j,coffs+1,tmp*q(1,1))
       end do
       call up_left_shift(q)
       q(nq,nq)=1.0_dp
       !
       ! Apply Q to row k+1
       !
       b_bv(k+1,lbw+2:lbw+nq)=matmul(b_bv(k+1,lbw+2:lbw+nq),transpose(q(1:nq-1,1:nq-1)))
       ! remove a superdiagonal from L that formed during downdating.
       if (dnl==0) then
          do j=2,nl
             rot=rgivens(get_el_br(b_bv,lbw,roffs+j,coffs+j),get_el_br(b_bv,lbw,roffs+j,coffs+j+1))
             call tbr_times_rotation(b_bv,n,lbw,ubw2,n-k-1,rot,coffs+j)
             call rotation_times_general(trp_rot(rot),q,j-1,j)
             call set_el_br(b_bv,lbw,roffs+j,coffs+j+1,0.0_dp)
          end do
       end if
       ! compress row k+1
       do j=nq-1,nl+1+dnl,-1
          rot=rgivens(get_el_br(b_bv,lbw,k+1,coffs+j), get_el_br(b_bv,lbw,k+1,coffs+j+1))
          call tbr_times_rotation(b_bv,n,lbw,ubw2,n-k-1,rot,coffs+j)
          call rotation_times_general(trp_rot(rot),q,j-1,j)
          call set_el_br(b_bv,lbw,k+1,coffs+j+1,0.0_dp)
       end do
       !
       ! Apply v_{n-k-1} to Q
       !
       if (n-k-1>0) then
          do j=1,numrots_bv(n-k-1)
             rot%cosine=cs_bv(n-k-1,j); rot%sine=ss_bv(n-k-1,j)
             call general_times_rotation(q,trp_rot(rot),k1s_bv(n-k-1,j)-coffs-1,k2s_bv(n-k-1,j)-coffs-1)
          end do
       end if
       ! Termination cases:
       !
       ! 1. If dnl==1 and k+nl==n-1 then b(k-nl+1:k+1,k+2:n) is lower trapezoidal.
       ! 2. If dnl==0 then k+nl==n and b(k-nl+2:k+1,k+2:n) is lower trapezoidal.
       ! 3. If dnl==1 and k+nl==n then b(k-nl+1:k+1,k+2:n) is square lower trapezoidal.
       !
       coffs = k+1
       if (dnl==1 .and. k+nl==n-1) then
          roffs=k-nl
          ml=nl+1
          exit kloop
       else if (k+nl==n .and. dnl==0) then
          roffs=k-nl+1
          ml=nl; nl=nl-1
          exit kloop
       else if (k+nl==n .and. dnl==1) then
          roffs=k-nl
          ml=nl+1; nl=nl-1
          exit kloop
       end if
       nl=nl+dnl
       ml=nl
    end do kloop
    ! Repeatedly eliminate the diagonal of L to finish the decomposition.
    do k=n-nl,n-1
       ubws(k)=ml-1
       nl=n-k
       coffs=k
       ! apply u_k
       numrots_ub(k) = nl
       do j=nl,1,-1
          rot=lgivens2(get_el_br(b_bv,lbw,roffs+j,coffs+j),get_el_br(b_bv,lbw,roffs+j+1,coffs+j))
          call rotation_times_tbr(trp_rot(rot), b_bv,n,lbw,ubw2,coffs,roffs+j)
          cs_ub(j,k)=rot%cosine; ss_ub(j,k)=rot%sine
          j1s_ub(j,k)=roffs+j; j2s_ub(j,k)=roffs+j+1
          call set_el_br(b_bv,lbw,roffs+j,coffs+j,0.0_dp)
       end do
       ! downdate
       do j=nl-1,1,-1
          rot=lgivens(q(j,1),q(j+1,1))
          call rotation_times_general(trp_rot(rot),q,j,j+1)
          q(j+1,1)=0.0_dp
          call tbr_times_rotation(b_bv,n,lbw,ubw2,n-k,rot,coffs+j)
       end do
       do j=1,ml
          tmp=get_el_br(b_bv,lbw,roffs+j,coffs+1)
          call set_el_br(b_bv,lbw,roffs+j,coffs+1,tmp*q(1,1))
       end do
       call up_left_shift(q)
       if (nl==1) then 
          exit
       end if
       ! add an extra row on to L
       b_bv(k+1,lbw+2:lbw+nl)=matmul(b_bv(k+1,lbw+2:lbw+nl),transpose(q(1:nl-1,1:nl-1)))
       ! apply v_{n-k-1} to Q
       do j=1,numrots_bv(n-k-1)
          rot%cosine=cs_bv(n-k-1,j); rot%sine=ss_bv(n-k-1,j)
          call general_times_rotation(q,trp_rot(rot),k1s_bv(n-k-1,j)-coffs-1,k2s_bv(n-k-1,j)-coffs-1)
       end do
       roffs=roffs+1
    end do
    !
    ! Extract diagonals
    !
    ubw_ub=maxval(ubws)
    ! put diagonals in b
    do d=1,ubw_ub+1
       do k=ubw_ub-d+2,n
          b_ub(d,k) = get_el_br(b_bv,lbw,d+k-ubw_ub-1,k)
       end do
    end do
    do d=ubw_ub+2, ubw_ub+lbw+1
       do k=1,n-d+ubw_ub+1
          b_ub(d,k) = get_el_br(b_bv,lbw,d+k-ubw_ub-1,k)
       end do
    end do
  end subroutine f_d_compress_bv_to_ub

  subroutine c_compress_bv_to_ub(bv, ub, told, tol, dr, error)
    type(c_ub) :: ub
    type(c_bv) :: bv
    integer(kind=int32), intent(out) :: error
    integer(kind=int32), intent(in)  :: dr
    real(kind=dp), intent(in) :: tol, told
    if (ub%n /= bv%n) then
       error = 3; return
    end if
    call f_c_compress_bv_to_ub(bv%b, bv%n, bv%lbw, bv%ubw, bv%lbwmax, bv%ubwmax, bv%numrotsv, bv%k1sv, bv%k2sv, &
         bv%csv, bv%ssv, ub%b, ub%ubw, ub%lbwmax, ub%ubwmax, ub%numrotsu, ub%j1su, ub%j2su, ub%csu, ub%ssu, &
         told, tol, dr, error)
    if (error == 0) then
       ub%lbw=bv%lbw; ub%n=bv%n
    end if
  end subroutine c_compress_bv_to_ub

  subroutine f_c_compress_bv_to_ub(b_bv, n, lbw, ubw, lbwmax_bv, ubwmax_bv, numrots_bv, &
       k1s_bv, k2s_bv, cs_bv, ss_bv, &
       b_ub, ubw_ub, lbwmax_ub, ubwmax_ub, numrots_ub, j1s_ub, j2s_ub, cs_ub, ss_ub, told, tol, dr, error)
    integer(kind=int32), intent(in) :: n, lbw, ubw, lbwmax_ub, ubwmax_ub, lbwmax_bv, ubwmax_bv, dr
    complex(kind=dp), dimension(n,lbwmax_bv+ubwmax_bv+1), intent(inout) :: b_bv
    integer(kind=int32), dimension(n), intent(in) :: numrots_bv
    integer(kind=int32), dimension(n,ubwmax_bv), intent(in) :: k1s_bv, k2s_bv
    complex(kind=dp), dimension(n,ubwmax_bv), intent(in) :: cs_bv, ss_bv

    complex(kind=dp), dimension(lbwmax_ub+ubwmax_ub+1,n), intent(out) :: b_ub
    integer(kind=int32), dimension(n), intent(out) :: numrots_ub
    integer(kind=int32), dimension(ubwmax_ub,n), intent(out) :: j1s_ub, j2s_ub
    complex(kind=dp), dimension(ubwmax_ub,n), intent(out) :: cs_ub, ss_ub

    real(kind=dp), intent(in) :: tol, told
    integer(kind=int32), intent(out) :: error, ubw_ub

    integer(kind=int32) :: j, k, jj, roffs, coffs, ubw2, nullerr, d, minindex, ml, nl, dnl, nq
    type(c_rotation) :: rot
    complex(kind=dp), target, dimension(ubw+1,ubw+1) :: q
    complex(kind=dp), dimension(ubw+1,ubw+1) :: l
    complex(kind=dp), dimension(ubw+1) :: x
    real(kind=dp) :: nrma, tmpr, mindiag
    complex(kind=dp) :: tmp
    integer(kind=int32), dimension(n) :: ubws

    error = 0
    numrots_ub=0
    ss_ub=(0.0_dp,0.0_dp); cs_ub=(0.0_dp,0.0_dp)
    j1s_ub=0; j2s_ub=0
    b_ub=(0.0_dp,0.0_dp)
    ubw2=ubw+2
    nrma = maxabs(b_bv)*sqrt(real(n))
    ubws=0

    if (n < 1) then
       error = 1
       return
    end if
    if (n == 1) then
        b_ub(1,1)=b_bv(1,1)
        return
    end if
    ! must allow for temporary fill-in of one extra superdiagonal in b_bv.
    if (lbwmax_bv+ubwmax_bv+1<ubw2+lbw+1) then
       error = 2
       return
    end if
    nl=1
    ml=1
    q=(0.0_dp,0.0_dp)
    do j=1,ubw+1
       q(j,j)=(1.0_dp,0.0_dp)
    end do
    ! provide working space in a superdiagonal on the right of b_bv
    b_bv(:,ubw+lbw+2:ubw+lbw+3)=(0.0_dp,0.0_dp)
    ! Initial LQ factorization.
    roffs=0; coffs=1
    do j=ubw-1,1,-1
       rot=rgivens(get_el_br(b_bv,lbw,1,coffs+j),get_el_br(b_bv,lbw,1,coffs+j+1))
       call tbr_times_rotation(b_bv,n,lbw,ubw2,n-1,rot,coffs+j)
       call rotation_times_general(trp_rot(rot), q, j,j+1)
       call set_el_br(b_bv,lbw,1,coffs+j+1,(0.0_dp,0.0_dp))
    end do
    ! Apply v_{n-1} to q.
    do j=1,numrots_bv(n-1)
       rot%cosine=cs_bv(n-1,j); rot%sine=ss_bv(n-1,j)
       call general_times_rotation(q,trp_rot(rot),k1s_bv(n-1,j)-1,k2s_bv(n-1,j)-1)
    end do
    ! main loop
    kloop: do k=1,n
       ! Current, possibly singular, L should be contained in
       ! b_bv(k-nl+1:k,k+1:k+nl)
       nq=min(ubw+1,n-k)
       roffs=k-nl; coffs=k
       mindiag=abs(get_el_br(b_bv,lbw,roffs+nl,coffs+nl))
       minindex=nl
       do j=nl-1,1,-1
          tmpr=abs(get_el_br(b_bv,lbw,roffs+j,coffs+j))
          if (tmpr <= mindiag) then
             minindex=j
             mindiag=tmpr
          end if
       end do
       if (mindiag <= told*nrma) then
          dnl=0  ! don't increase nl
          ubws(k)=nl-1
          numrots_ub(k) = minindex-1
          do j=minindex,2,-1
             rot=lgivens2(get_el_br(b_bv,lbw,roffs+j-1,coffs+j-1),get_el_br(b_bv,lbw,roffs+j,coffs+j-1))
             call rotation_times_tbr(trp_rot(rot), b_bv,n,lbw,ubw2,coffs,roffs+j-1)
             cs_ub(j-1,k)=rot%cosine; ss_ub(j-1,k)=rot%sine
             j1s_ub(j-1,k)=roffs+j-1; j2s_ub(j-1,k)=roffs+j
             call set_el_br(b_bv,lbw,roffs+j-1,coffs+j-1,(0.0_dp,0.0_dp))
             ! swap columns in L
             do jj=j-1,nl
                tmp=get_el_br(b_bv,lbw,roffs+jj,coffs+j)
                call set_el_br(b_bv,lbw,roffs+jj,coffs+j, &
                     get_el_br(b_bv,lbw,roffs+jj,coffs+j-1))
                call set_el_br(b_bv,lbw,roffs+jj,coffs+j-1,tmp)
             end do
             do jj=1,ubw+1
                tmp=q(j,jj)
                q(j,jj)=q(j-1,jj)
                q(j-1,jj)=tmp
             end do
          end do
          call set_el_br(b_bv,lbw,roffs+1,coffs+1,(0.0_dp,0.0_dp))
       else ! find a null vector
          call submatrix_br(b_bv,lbw,ubw2, roffs+1,roffs+nl,coffs+1,coffs+nl,l(1:nl,1:nl))
          call f_c_lower_left_nullvec(x(1:nl),l(1:nl,1:nl),tol*nrma,nullmaxits, nullerr)
          if ((nullerr >= 0 .and. tol > 0.0_dp) .or. &
               (tol==0.0_dp .and. told == 0.0_dp .and. nl > ubw-dr)) then ! null vector found
             dnl=0
             ubws(k)=nl-1
             ! Introduce a zero while preserving the triangularity of L
             numrots_ub(k)=nl-1
             do j=nl-1,1,-1
                rot=rgivens(x(j),x(j+1))
                call general_times_rotation(x,rot,j,j+1)
                call rotation_times_tbr(trp_rot(rot),b_bv,n,lbw,ubw2,coffs,roffs+j)
                cs_ub(j,k)=rot%cosine; ss_ub(j,k)=rot%sine
                j1s_ub(j,k)=roffs+j; j2s_ub(j,k)=roffs+j+1
                rot=rgivens(get_el_br(b_bv,lbw,roffs+j,coffs+j),get_el_br(b_bv,lbw,roffs+j,coffs+j+1))
                call tbr_times_rotation(b_bv,n,lbw,ubw2,n-k,rot,coffs+j)
                call rotation_times_general(trp_rot(rot),q,j,j+1)
                call set_el_br(b_bv,lbw,roffs+j,coffs+j+1,(0.0_dp,0.0_dp))
             end do
             call set_el_br(b_bv,lbw,roffs+1,coffs+1,(0.0_dp,0.0_dp))
          else
             ! Compression has failed; increase nl.
             ubws(k)=nl
             dnl=1
          end if
       end if
       !
       ! Downdate
       ! 
       do j=nq-1,1,-1
          rot=lgivens(q(j,1),q(j+1,1))
          call rotation_times_general(trp_rot(rot),q,j,j+1)
          q(j+1,1)=(0.0_dp,0.0_dp)
          call tbr_times_rotation(b_bv,n,lbw,ubw2,n-k,rot,coffs+j)
       end do
       do j=1,nl
          tmp=get_el_br(b_bv,lbw,roffs+j,coffs+1)
          call set_el_br(b_bv,lbw,roffs+j,coffs+1,tmp*q(1,1))
       end do
       call up_left_shift(q)
       q(nq,nq)=(1.0_dp,0.0_dp)
       !
       ! Apply Q to row k+1
       !
       b_bv(k+1,lbw+2:lbw+nq)=matmul(b_bv(k+1,lbw+2:lbw+nq),transpose(conjg(q(1:nq-1,1:nq-1))))
       ! remove a superdiagonal from L that formed during downdating.
       if (dnl==0) then
          do j=2,nl
             rot=rgivens(get_el_br(b_bv,lbw,roffs+j,coffs+j),get_el_br(b_bv,lbw,roffs+j,coffs+j+1))
             call tbr_times_rotation(b_bv,n,lbw,ubw2,n-k-1,rot,coffs+j)
             call rotation_times_general(trp_rot(rot),q,j-1,j)
             call set_el_br(b_bv,lbw,roffs+j,coffs+j+1,(0.0_dp,0.0_dp))
          end do
       end if
       ! compress row k+1
       do j=nq-1,nl+1+dnl,-1
          rot=rgivens(get_el_br(b_bv,lbw,k+1,coffs+j), get_el_br(b_bv,lbw,k+1,coffs+j+1))
          call tbr_times_rotation(b_bv,n,lbw,ubw2,n-k-1,rot,coffs+j)
          call rotation_times_general(trp_rot(rot),q,j-1,j)
          call set_el_br(b_bv,lbw,k+1,coffs+j+1,(0.0_dp,0.0_dp))
       end do
       !
       ! Apply v_{n-k-1} to Q
       !
       if (n-k-1>0) then
          do j=1,numrots_bv(n-k-1)
             rot%cosine=cs_bv(n-k-1,j); rot%sine=ss_bv(n-k-1,j)
             call general_times_rotation(q,trp_rot(rot),k1s_bv(n-k-1,j)-coffs-1,k2s_bv(n-k-1,j)-coffs-1)
          end do
       end if
       ! Termination cases:
       !
       ! 1. If dnl==1 and k+nl==n-1 then b(k-nl+1:k+1,k+2:n) is lower trapezoidal.
       ! 2. If dnl==0 then k+nl==n and b(k-nl+2:k+1,k+2:n) is lower trapezoidal.
       ! 3. If dnl==1 and k+nl==n then b(k-nl+1:k+1,k+2:n) is square lower trapezoidal.
       !
       coffs = k+1
       if (dnl==1 .and. k+nl==n-1) then
          roffs=k-nl
          ml=nl+1
          exit kloop
       else if (k+nl==n .and. dnl==0) then
          roffs=k-nl+1
          ml=nl; nl=nl-1
          exit kloop
       else if (k+nl==n .and. dnl==1) then
          roffs=k-nl
          ml=nl+1; nl=nl-1
          exit kloop
       end if
       nl=nl+dnl
       ml=nl
    end do kloop
    ! Repeatedly eliminate the diagonal of L to finish the decomposition.
    do k=n-nl,n-1
       ubws(k)=ml-1
       nl=n-k
       coffs=k
       ! apply u_k
       numrots_ub(k) = nl
       do j=nl,1,-1
          rot=lgivens2(get_el_br(b_bv,lbw,roffs+j,coffs+j),get_el_br(b_bv,lbw,roffs+j+1,coffs+j))
          call rotation_times_tbr(trp_rot(rot), b_bv,n,lbw,ubw2,coffs,roffs+j)
          cs_ub(j,k)=rot%cosine; ss_ub(j,k)=rot%sine
          j1s_ub(j,k)=roffs+j; j2s_ub(j,k)=roffs+j+1
          call set_el_br(b_bv,lbw,roffs+j,coffs+j,(0.0_dp,0.0_dp))
       end do
       ! downdate
       do j=nl-1,1,-1
          rot=lgivens(q(j,1),q(j+1,1))
          call rotation_times_general(trp_rot(rot),q,j,j+1)
          q(j+1,1)=(0.0_dp,0.0_dp)
          call tbr_times_rotation(b_bv,n,lbw,ubw2,n-k,rot,coffs+j)
       end do
       do j=1,ml
          tmp=get_el_br(b_bv,lbw,roffs+j,coffs+1)
          call set_el_br(b_bv,lbw,roffs+j,coffs+1,tmp*q(1,1))
       end do
       call up_left_shift(q)
       if (nl==1) then 
          exit
       end if
       ! add an extra row on to L
       b_bv(k+1,lbw+2:lbw+nl)=matmul(b_bv(k+1,lbw+2:lbw+nl),transpose(conjg(q(1:nl-1,1:nl-1))))
       ! apply v_{n-k-1} to Q
       do j=1,numrots_bv(n-k-1)
          rot%cosine=cs_bv(n-k-1,j); rot%sine=ss_bv(n-k-1,j)
          call general_times_rotation(q,trp_rot(rot),k1s_bv(n-k-1,j)-coffs-1,k2s_bv(n-k-1,j)-coffs-1)
       end do
       roffs=roffs+1
    end do
    !
    ! Extract diagonals
    !
    ubw_ub=maxval(ubws)
    ! put diagonals in b
    do d=1,ubw_ub+1
       do k=ubw_ub-d+2,n
          b_ub(d,k) = get_el_br(b_bv,lbw,d+k-ubw_ub-1,k)
       end do
    end do
    do d=ubw_ub+2, ubw_ub+lbw+1
       do k=1,n-d+ubw_ub+1
          b_ub(d,k) = get_el_br(b_bv,lbw,d+k-ubw_ub-1,k)
       end do
    end do
  end subroutine f_c_compress_bv_to_ub


end module compressions_bv_to_ub
