module compressions_ub_to_bv
use prec
use shift
use rotation
use band_types
use nested_types
use nullvec
implicit none

integer(kind=int32), private, parameter :: nullmaxits=5

interface compress_ub_to_bv
   module procedure d_compress_ub_to_bv, c_compress_ub_to_bv
end interface compress_ub_to_bv

interface f_compress_ub_to_bv
   module procedure f_d_compress_ub_to_bv, f_c_compress_ub_to_bv
end interface f_compress_ub_to_bv

contains

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

  subroutine d_compress_ub_to_bv(ub, bv, told, tol, dr, error)
    type(d_ub) :: ub
    type(d_bv) :: bv
    integer(kind=int32), intent(out) :: error
    integer(kind=int32), intent(in) :: dr
    real(kind=dp), intent(in) :: tol, told
    if (ub%n /= bv%n) then
       error = 3; return
    end if
    call f_d_compress_ub_to_bv(ub%b, ub%n, ub%lbw, ub%ubw, ub%lbwmax, ub%ubwmax, ub%numrotsu, ub%j1su, ub%j2su, &
         ub%csu, ub%ssu, bv%b, bv%ubw, bv%lbwmax, bv%ubwmax, bv%numrotsv, bv%k1sv, bv%k2sv, bv%csv, bv%ssv, &
         told, tol, dr, error)
    if (error == 0) then
       bv%lbw=ub%lbw; bv%n=ub%n
    end if
  end subroutine d_compress_ub_to_bv


  subroutine f_d_compress_ub_to_bv(b_ub, n, lbw, ubw, lbwmax_ub, ubwmax_ub, numrots_ub, &
       j1s_ub, j2s_ub, cs_ub, ss_ub, &
       b_bv, ubw_bv, lbwmax_bv, ubwmax_bv, numrots_bv, k1s_bv, k2s_bv, cs_bv, ss_bv, told, tol, dr, error)
    real(kind=dp), dimension(lbwmax_ub+ubwmax_ub+1,n), intent(inout) :: b_ub
    integer(kind=int32), intent(in) :: n, lbw, ubw, lbwmax_ub, ubwmax_ub, lbwmax_bv, ubwmax_bv, dr
    integer(kind=int32), dimension(n), intent(in) :: numrots_ub
    integer(kind=int32), dimension(ubwmax_ub,n), intent(in) :: j1s_ub, j2s_ub
    real(kind=dp), dimension(ubwmax_ub,n), intent(in) :: cs_ub, ss_ub
    real(kind=dp), intent(in) :: tol, told

    real(kind=dp), dimension(n,lbwmax_bv+ubwmax_bv+1), intent(out) :: b_bv
    integer(kind=int32), dimension(n), intent(out) :: numrots_bv
    integer(kind=int32), dimension(n,ubwmax_bv), intent(out) :: k1s_bv, k2s_bv
    real(kind=dp), dimension(n,ubwmax_bv), intent(out) :: cs_bv, ss_bv
    integer(kind=int32), intent(out) :: error, ubw_bv

    integer(kind=int32) :: j, k, kk, roffs, coffs, qoffs, ubw2, nullerr, d, minindex, ml, nl, dnl, mq
    type(d_rotation) :: rot
    real(kind=dp), target, dimension(ubw+1,ubw+1) :: q
    real(kind=dp), pointer, dimension(:,:) :: pq
    real(kind=dp), dimension(ubw+1,ubw+1) :: l
    real(kind=dp), dimension(ubw+1) :: x
    real(kind=dp) :: nrma, tmp, mindiag
    integer(kind=int32), dimension(n) :: ubws

    error = 0
    numrots_bv=0
    ss_bv=0.0_dp; cs_bv=0.0_dp
    k1s_bv=0; k2s_bv=0
    b_bv=0.0_dp
    ubw2=ubw+2
    nrma = maxabs(b_ub)*sqrt(real(n))
    ubws=0

    if (n < 1) then
       error = 1
       return
    end if
    if (n == 1) then
        b_bv(1,1)=b_ub(1,1)
        ubw_bv=0
        return
    end if

    ! must allow for temporary fill-in of two extra superdiagonals.
    if (lbwmax_ub+ubwmax_ub+1<ubw+lbw+3) then
       error = 2
       return
    end if
    ! provide working space in an extra superdiagonal on top of b_ub
    call down_shift(b_ub)
    call down_shift(b_ub)
    nl=1
    ml=1
    ! initialize q
    q=0.0_dp
    do j=1,ubw+1
       q(j,j)=1.0_dp
    end do
    roffs=n-ubw-1; coffs=n-1
    qoffs=n-ubw-2
    mq=min(ubw+1,n-1)
    ! Initial QL Factorization
    do k=1,ubw-1
       rot=lgivens2(get_el_bc(b_ub,ubw2,roffs+k,coffs+1),get_el_bc(b_ub,ubw2,roffs+k+1,coffs+1))
       call rotation_times_tbc(trp_rot(rot),b_ub,n,lbw,ubw2,n-1,roffs+k)
       call general_times_rotation(q,rot,roffs-qoffs+k,roffs-qoffs+k+1)
       call set_el_bc(b_ub,ubw2,roffs+k,coffs+1,0.0_dp)
    end do
    ! Apply u_{n-1} to q
    do k=1,numrots_ub(n-1)
       rot%cosine=cs_ub(k,n-1); rot%sine=ss_ub(k,n-1)
       call rotation_times_general(rot,q,j1s_ub(k,n-1)-qoffs,j2s_ub(k,n-1)-qoffs)
    end do
    ! main loop
    jloop: do j=1,n
       ! Current, possibly singular, square L should be in
       ! b_ub(n-j-nl+1:n-j,n-j+1:n-j+nl)
       mq=min(ubw+1,n-j)
       pq => q(ubw+2-mq:ubw+1,ubw+2-mq:ubw+1)
       qoffs=n-j-mq ! qoffs+1 is the first row acted on by q
       roffs=n-j-nl; coffs=n-j ! (roffs+1,coffs+1) is the upper left corner of L
       mindiag=abs(get_el_bc(b_ub,ubw2,roffs+1,coffs+1))
       minindex=1
       do k=2,nl
          tmp=abs(get_el_bc(b_ub,ubw2,roffs+k,coffs+k))
          if (tmp <= mindiag) then
             minindex=k
             mindiag=tmp
          end if
       end do
       if (mindiag <= told * nrma) then
          dnl=0
          ubws(j)=nl-1
          numrots_bv(j)=nl-minindex
          do k=minindex,nl-1
             rot=rgivens(get_el_bc(b_ub,ubw2,roffs+k+1,coffs+k),get_el_bc(b_ub,ubw2,roffs+k+1,coffs+k+1))
             call tbc_times_rotation(b_ub,n,lbw,ubw2,j,rot,coffs+k)
             cs_bv(j,nl-k)=rot%cosine; ss_bv(j,nl-k)=rot%sine
             k1s_bv(j,nl-k)=coffs+k; k2s_bv(j,nl-k)=coffs+k+1
             call set_el_bc(b_ub,ubw2,roffs+k+1,coffs+k+1,0.0_dp)
             ! swap rows in L
             do kk=1,k+1
                tmp=get_el_bc(b_ub,ubw2,roffs+k,coffs+kk)
                call set_el_bc(b_ub,ubw2,roffs+k,coffs+kk, &
                     get_el_bc(b_ub,ubw2,roffs+k+1,coffs+kk))
                call set_el_bc(b_ub,ubw2,roffs+k+1,coffs+kk,tmp)
             end do
             do kk=1,mq
                tmp=pq(kk,roffs-qoffs+k)
                pq(kk,roffs-qoffs+k)=pq(kk,roffs-qoffs+k+1)
                pq(kk,roffs-qoffs+k+1)=tmp
             end do
          end do
          call set_el_bc(b_ub,ubw2,roffs+nl,coffs+nl,0.0_dp)
       else ! find a null vector
          call submatrix_bc(b_ub,lbw,ubw2,roffs+1,roffs+nl,coffs+1,coffs+nl,l(1:nl,1:nl))
          call f_d_lower_right_nullvec(x(1:nl),l(1:nl,1:nl),tol*nrma,nullmaxits,nullerr)
          if ((nullerr>=0 .and. tol > 0.0_dp) .or. &
               (tol==0.0_dp .and. told == 0.0_dp .and. nl > ubw-dr)) then ! null vector found
             dnl = 0
             ubws(j)=nl-1
             ! introduce a zero while preserving triangularity of L
             numrots_bv(j)=nl-1
             do k=1,nl-1
                rot=lgivens2(x(k),x(k+1))
                call rotation_times_general(trp_rot(rot),x,k,k+1)
                call tbc_times_rotation(b_ub,n,lbw,ubw2,j,rot,coffs+k)
                cs_bv(j,nl-k)=rot%cosine; ss_bv(j,nl-k)=rot%sine
                k1s_bv(j,nl-k)=coffs+k; k2s_bv(j,nl-k)=coffs+k+1
                rot=lgivens2(get_el_bc(b_ub,ubw2,roffs+k,coffs+k+1),get_el_bc(b_ub,ubw2,roffs+k+1,coffs+k+1))
                call rotation_times_tbc(trp_rot(rot),b_ub,n,lbw,ubw2,n-j,roffs+k)
                call general_times_rotation(pq,rot,roffs-qoffs+k,roffs-qoffs+k+1)
                call set_el_bc(b_ub,ubw2,roffs+k,coffs+k+1,0.0_dp)
             end do
             call set_el_bc(b_ub,ubw2,roffs+nl,coffs+nl,0.0_dp)
          else
             ! Compression has failed; increase nl.
             ubws(j)=nl
             dnl=1
          end if
       end if
       !
       ! Downdate
       ! 
       do k=1,mq-1
          rot=rgivens2(pq(mq,k),pq(mq,k+1))
          call general_times_rotation(pq,rot,k,k+1)
          pq(mq,k)=0.0_dp
          call rotation_times_tbc(trp_rot(rot),b_ub,n,lbw,ubw2,n-j,qoffs+k)
       end do
       do k=1,nl
          tmp=get_el_bc(b_ub,ubw2,roffs+nl,coffs+k)
          call set_el_bc(b_ub,ubw2,roffs+nl,coffs+k,tmp*pq(mq,mq))
       end do
       call down_right_shift(pq)
       pq(1,1)=1.0_dp
       qoffs=qoffs-1
       ! Eliminate a diagonal.
       if (dnl==0) then
          do k=nl-1,1,-1
             rot=lgivens2(get_el_bc(b_ub,ubw2,roffs+k-1,coffs+k),get_el_bc(b_ub,ubw2,roffs+k,coffs+k))
             call rotation_times_tbc(trp_rot(rot),b_ub,n,lbw,ubw2,n-j,roffs+k-1)
             call general_times_rotation(pq,rot,roffs-qoffs+k-1,roffs-qoffs+k)
             call set_el_bc(b_ub,ubw2,roffs+k-1,coffs+k,0.0_dp)
          end do
       end if
       ! apply q to column n-j
       b_ub(ubw2-mq+1:ubw2,n-j)=matmul(transpose(pq),b_ub(ubw2-mq+1:ubw2,n-j))
       ! compress column n-j
       do k=-(mq-nl)+2,-dnl
          rot=lgivens2(get_el_bc(b_ub,ubw2,roffs+k-1,n-j), get_el_bc(b_ub,ubw2,roffs+k,n-j))
          call rotation_times_tbc(trp_rot(rot),b_ub,n,lbw,ubw2,n-j-1,roffs+k-1)
          call general_times_rotation(pq,rot,roffs-qoffs+k-1,roffs-qoffs+k)
          call set_el_bc(b_ub,ubw2,roffs+k-1,n-j,0.0_dp)
       end do
       ! Apply u_{n-j-1} to Q
       if (n-j-1 >0) then
          do k=1,numrots_ub(n-j-1)
             rot%cosine=cs_ub(k,n-j-1); rot%sine=ss_ub(k,n-j-1)
             call rotation_times_general(rot,pq,j1s_ub(k,n-j-1)-qoffs,j2s_ub(k,n-j-1)-qoffs)
          end do
       end if
       ! Termination cases:
       ! 1. If dnl==1 and j+nl==n-1 then b(1:nl,nl+1:2*nl+1) is lower trapezoidal
       ! 2. If dnl==0 and j+nl==n then b(1:nl-1,nl:2*nl-1) is lower trapezoidal
       ! 3. If dnl==1 and j+nl==n then b(1:nl-1,nl:2*nl) is lower trapezoidal
       roffs=0
       if (dnl==1 .and. j+nl==n-1) then
          nl=ml+1
          exit jloop
       else if (j+nl==n .and. dnl==0) then
          ml=nl-1
          exit jloop
       else if (j+nl==n .and. dnl==1) then
          ml=nl-1; nl=nl+1
          exit jloop
       end if
       nl=nl+dnl
       ml=nl
    end do jloop
    ! Repeatedly eliminate the super diagonal to finish the decomposition.
    qoffs=0; roffs=0
    do j=n-ml,n-1
       ubws(j)=nl-1
       ml=n-j
       ! Note that mq=ml
       pq => q(ubw+2-ml:ubw+1,ubw+2-ml:ubw+1)
       coffs=n-j
       ! apply v_j
       numrots_bv(j)=ml
       do k=1,ml
          rot=rgivens(get_el_bc(b_ub,ubw2,k,coffs+nl-ml+k-1),get_el_bc(b_ub,ubw2,k,coffs+nl-ml+k))
          call tbc_times_rotation(b_ub,n,lbw,ubw2,j,rot,coffs+nl-ml+k-1)
          cs_bv(j,ml-k+1)=rot%cosine; ss_bv(j,ml-k+1)=rot%sine
          k1s_bv(j,ml-k+1)=coffs+nl-ml+k-1; k2s_bv(j,ml-k+1)=coffs+nl-ml+k
          call set_el_bc(b_ub,ubw2,k,coffs+nl-ml+k,0.0_dp)
       end do
       ! downdate
       do k=1,ml-1
          rot=rgivens2(pq(ml,k),pq(ml,k+1))
          call general_times_rotation(pq,rot,k,k+1)
          pq(ml,k)=0.0_dp
          call rotation_times_tbc(trp_rot(rot),b_ub,n,lbw,ubw2,ml,k)
       end do
       do k=1,nl
          tmp=get_el_bc(b_ub,ubw2,ml,ml+k)
          call set_el_bc(b_ub,ubw2,ml,ml+k,tmp*pq(ml,ml))
       end do
       if (ml==1) then
          exit
       end if
       ! add an extra column to L
       b_ub((ubw2-(n-j)+1)+1:(ubw2-(n-j)+1)+ml-1,n-j)=matmul(transpose(pq(1:ml-1,1:ml-1)), &
            b_ub((ubw2-(n-j)+1)+1:(ubw2-(n-j)+1)+ml-1,n-j))
       ! apply u_{n-j-1} to Q
       do k=1,numrots_ub(n-j-1)
          rot%cosine=cs_ub(k,n-j-1); rot%sine=ss_ub(k,n-j-1)
          call rotation_times_general(rot,pq,j1s_ub(k,n-j-1),j2s_ub(k,n-j-1))
       end do
       call down_right_shift(pq)
    end do
    ubw_bv=maxval(ubws)
    ! Extract diagonals
    do d=1,lbw+1
       do j=lbw-d+2, n
          b_bv(j,d)=get_el_bc(b_ub,ubw2,j,d+j-lbw-1)
       end do
    end do
    do d=lbw+2,ubw_bv+lbw+1
       do j=1, n-d+lbw+1
          b_bv(j,d)=get_el_bc(b_ub,ubw2,j,d+j-lbw-1)
       end do
    end do
  end subroutine f_d_compress_ub_to_bv



  subroutine c_compress_ub_to_bv(ub, bv, told, tol, dr, error)
    type(c_ub) :: ub
    type(c_bv) :: bv
    integer(kind=int32), intent(out) :: error
    integer(kind=int32), intent(in) :: dr
    real(kind=dp), intent(in) :: tol, told
    if (ub%n /= bv%n) then
       error = 3; return
    end if
    call f_c_compress_ub_to_bv(ub%b, ub%n, ub%lbw, ub%ubw, ub%lbwmax, ub%ubwmax, ub%numrotsu, ub%j1su, ub%j2su, &
         ub%csu, ub%ssu, bv%b, bv%ubw, bv%lbwmax, bv%ubwmax, bv%numrotsv, bv%k1sv, bv%k2sv, bv%csv, bv%ssv, &
         told, tol, dr, error)
    if (error == 0) then
       bv%lbw=ub%lbw; bv%n=ub%n
    end if
  end subroutine c_compress_ub_to_bv

  subroutine f_c_compress_ub_to_bv(b_ub, n, lbw, ubw, lbwmax_ub, ubwmax_ub, numrots_ub, &
       j1s_ub, j2s_ub, cs_ub, ss_ub, &
       b_bv, ubw_bv, lbwmax_bv, ubwmax_bv, numrots_bv, k1s_bv, k2s_bv, cs_bv, ss_bv, told, tol, dr, error)
    complex(kind=dp), dimension(lbwmax_ub+ubwmax_ub+1,n), intent(inout) :: b_ub
    integer(kind=int32), intent(in) :: n, lbw, ubw, lbwmax_ub, ubwmax_ub, lbwmax_bv, ubwmax_bv, dr
    integer(kind=int32), dimension(n), intent(in) :: numrots_ub
    integer(kind=int32), dimension(ubwmax_ub,n), intent(in) :: j1s_ub, j2s_ub
    complex(kind=dp), dimension(ubwmax_ub,n), intent(in) :: cs_ub, ss_ub
    real(kind=dp), intent(in) :: tol, told

    complex(kind=dp), dimension(n,lbwmax_bv+ubwmax_bv+1), intent(out) :: b_bv
    integer(kind=int32), dimension(n), intent(out) :: numrots_bv
    integer(kind=int32), dimension(n,ubwmax_bv), intent(out) :: k1s_bv, k2s_bv
    complex(kind=dp), dimension(n,ubwmax_bv), intent(out) :: cs_bv, ss_bv
    integer(kind=int32), intent(out) :: error, ubw_bv

    integer(kind=int32) :: j, k, kk, roffs, coffs, qoffs, ubw2, nullerr, d, minindex, ml, nl, dnl, mq
    type(c_rotation) :: rot
    complex(kind=dp), target, dimension(ubw+1,ubw+1) :: q
    complex(kind=dp), pointer, dimension(:,:) :: pq
    complex(kind=dp), dimension(ubw+1,ubw+1) :: l
    complex(kind=dp), dimension(ubw+1) :: x
    real(kind=dp) :: nrma, tmpr, mindiag
    complex(kind=dp) :: tmp
    integer(kind=int32), dimension(n) :: ubws

    error = 0
    numrots_bv=0
    ss_bv=(0.0_dp, 0.0_dp); cs_bv=(0.0_dp, 0.0_dp)
    k1s_bv=0; k2s_bv=0
    b_bv=(0.0_dp, 0.0_dp)
    ubw2=ubw+2
    nrma = maxabs(b_ub)*sqrt(real(n))
    ubws=0

    if (n < 1) then
       error = 1
       return
    end if
    if (n == 1) then
        b_bv(1,1)=b_ub(1,1)
        return
    end if

    ! must allow for temporary fill-in of two extra superdiagonals.
    if (lbwmax_ub+ubwmax_ub+1<ubw+lbw+3) then
       error = 2
       return
    end if
    ! provide working space in an extra superdiagonal on top of b_ub
    call down_shift(b_ub)
    call down_shift(b_ub)
    nl=1
    ml=1
    ! initialize q
    q=(0.0_dp, 0.0_dp)
    do j=1,ubw+1
       q(j,j)=(1.0_dp,0.0_dp)
    end do
    roffs=n-ubw-1; coffs=n-1
    qoffs=n-ubw-2
    mq=min(ubw+1,n-1)
    ! Initial QL Factorization
    do k=1,ubw-1
       rot=lgivens2(get_el_bc(b_ub,ubw2,roffs+k,coffs+1),get_el_bc(b_ub,ubw2,roffs+k+1,coffs+1))
       call rotation_times_tbc(trp_rot(rot),b_ub,n,lbw,ubw2,n-1,roffs+k)
       call general_times_rotation(q,rot,roffs-qoffs+k,roffs-qoffs+k+1)
       call set_el_bc(b_ub,ubw2,roffs+k,coffs+1,(0.0_dp, 0.0_dp))
    end do
    ! Apply u_{n-1} to q
    do k=1,numrots_ub(n-1)
       rot%cosine=cs_ub(k,n-1); rot%sine=ss_ub(k,n-1)
       call rotation_times_general(rot,q,j1s_ub(k,n-1)-qoffs,j2s_ub(k,n-1)-qoffs)
    end do
    ! main loop
    jloop: do j=1,n
       ! Current, possibly singular, square L should be in
       ! b_ub(n-j-nl+1:n-j,n-j+1:n-j+nl)
       mq=min(ubw+1,n-j)
       pq => q(ubw+2-mq:ubw+1,ubw+2-mq:ubw+1)
       qoffs=n-j-mq ! qoffs+1 is the first row acted on by q
       roffs=n-j-nl; coffs=n-j ! (roffs+1,coffs+1) is the upper left corner of L
       mindiag=abs(get_el_bc(b_ub,ubw2,roffs+1,coffs+1))
       minindex=1
       do k=2,nl
          tmpr=abs(get_el_bc(b_ub,ubw2,roffs+k,coffs+k))
          if (tmpr <= mindiag) then
             minindex=k
             mindiag=tmpr
          end if
       end do
       if (mindiag <= told * nrma) then
          dnl=0
          ubws(j)=nl-1
          numrots_bv(j)=nl-minindex
          do k=minindex,nl-1
             rot=rgivens(get_el_bc(b_ub,ubw2,roffs+k+1,coffs+k),get_el_bc(b_ub,ubw2,roffs+k+1,coffs+k+1))
             call tbc_times_rotation(b_ub,n,lbw,ubw2,j,rot,coffs+k)
             cs_bv(j,nl-k)=rot%cosine; ss_bv(j,nl-k)=rot%sine
             k1s_bv(j,nl-k)=coffs+k; k2s_bv(j,nl-k)=coffs+k+1
             call set_el_bc(b_ub,ubw2,roffs+k+1,coffs+k+1,(0.0_dp, 0.0_dp))
             ! swap rows in L
             do kk=1,k+1
                tmp=get_el_bc(b_ub,ubw2,roffs+k,coffs+kk)
                call set_el_bc(b_ub,ubw2,roffs+k,coffs+kk, &
                     get_el_bc(b_ub,ubw2,roffs+k+1,coffs+kk))
                call set_el_bc(b_ub,ubw2,roffs+k+1,coffs+kk,tmp)
             end do
             do kk=1,mq
                tmp=pq(kk,roffs-qoffs+k)
                pq(kk,roffs-qoffs+k)=pq(kk,roffs-qoffs+k+1)
                pq(kk,roffs-qoffs+k+1)=tmp
             end do
          end do
          call set_el_bc(b_ub,ubw2,roffs+nl,coffs+nl,(0.0_dp, 0.0_dp))
       else ! find a null vector
          call submatrix_bc(b_ub,lbw,ubw2,roffs+1,roffs+nl,coffs+1,coffs+nl,l(1:nl,1:nl))
          call f_c_lower_right_nullvec(x(1:nl),l(1:nl,1:nl),tol*nrma,nullmaxits,nullerr)
          if ((nullerr>=0 .and. tol > 0.0_dp) .or. &
               (tol==0.0_dp .and. told == 0.0_dp .and. nl > ubw-dr)) then ! null vector found
             dnl = 0
             ubws(j)=nl-1
             ! introduce a zero while preserving triangularity of L
             numrots_bv(j)=nl-1
             do k=1,nl-1
                rot=lgivens2(x(k),x(k+1))
                call rotation_times_general(trp_rot(rot),x,k,k+1)
                call tbc_times_rotation(b_ub,n,lbw,ubw2,j,rot,coffs+k)
                cs_bv(j,nl-k)=rot%cosine; ss_bv(j,nl-k)=rot%sine
                k1s_bv(j,nl-k)=coffs+k; k2s_bv(j,nl-k)=coffs+k+1
                rot=lgivens2(get_el_bc(b_ub,ubw2,roffs+k,coffs+k+1),get_el_bc(b_ub,ubw2,roffs+k+1,coffs+k+1))
                call rotation_times_tbc(trp_rot(rot),b_ub,n,lbw,ubw2,n-j,roffs+k)
                call general_times_rotation(pq,rot,roffs-qoffs+k,roffs-qoffs+k+1)
                call set_el_bc(b_ub,ubw2,roffs+k,coffs+k+1,(0.0_dp, 0.0_dp))
             end do
             call set_el_bc(b_ub,ubw2,roffs+nl,coffs+nl,(0.0_dp, 0.0_dp))
          else
             ! Compression has failed; increase nl.
             ubws(j)=nl
             dnl=1
          end if
       end if
       !
       ! Downdate
       ! 
       do k=1,mq-1
          rot=rgivens2(pq(mq,k),pq(mq,k+1))
          call general_times_rotation(pq,rot,k,k+1)
          pq(mq,k)=(0.0_dp, 0.0_dp)
          call rotation_times_tbc(trp_rot(rot),b_ub,n,lbw,ubw2,n-j,qoffs+k)
       end do
       do k=1,nl
          tmp=get_el_bc(b_ub,ubw2,roffs+nl,coffs+k)
          call set_el_bc(b_ub,ubw2,roffs+nl,coffs+k,tmp*pq(mq,mq))
       end do
       call down_right_shift(pq)
       pq(1,1)=(1.0_dp,0.0_dp)
       qoffs=qoffs-1
       ! Eliminate a diagonal.
       if (dnl==0) then
          do k=nl-1,1,-1
             rot=lgivens2(get_el_bc(b_ub,ubw2,roffs+k-1,coffs+k),get_el_bc(b_ub,ubw2,roffs+k,coffs+k))
             call rotation_times_tbc(trp_rot(rot),b_ub,n,lbw,ubw2,n-j,roffs+k-1)
             call general_times_rotation(pq,rot,roffs-qoffs+k-1,roffs-qoffs+k)
             call set_el_bc(b_ub,ubw2,roffs+k-1,coffs+k,(0.0_dp, 0.0_dp))
          end do
       end if
       ! apply q to column n-j
       b_ub(ubw2-mq+1:ubw2,n-j)=matmul(transpose(conjg(pq)),b_ub(ubw2-mq+1:ubw2,n-j))
       ! compress column n-j
       do k=-(mq-nl)+2,-dnl
          rot=lgivens2(get_el_bc(b_ub,ubw2,roffs+k-1,n-j), get_el_bc(b_ub,ubw2,roffs+k,n-j))
          call rotation_times_tbc(trp_rot(rot),b_ub,n,lbw,ubw2,n-j-1,roffs+k-1)
          call general_times_rotation(pq,rot,roffs-qoffs+k-1,roffs-qoffs+k)
          call set_el_bc(b_ub,ubw2,roffs+k-1,n-j,(0.0_dp, 0.0_dp))
       end do
       ! Apply u_{n-j-1} to Q
       if (n-j-1 >0) then
          do k=1,numrots_ub(n-j-1)
             rot%cosine=cs_ub(k,n-j-1); rot%sine=ss_ub(k,n-j-1)
             call rotation_times_general(rot,pq,j1s_ub(k,n-j-1)-qoffs,j2s_ub(k,n-j-1)-qoffs)
          end do
       end if
       ! Termination cases:
       ! 1. If dnl==1 and j+nl==n-1 then b(1:nl,nl+1:2*nl+1) is lower trapezoidal
       ! 2. If dnl==0 and j+nl==n then b(1:nl-1,nl:2*nl-1) is lower trapezoidal
       ! 3. If dnl==1 and j+nl==n then b(1:nl-1,nl:2*nl) is lower trapezoidal
       roffs=0
       if (dnl==1 .and. j+nl==n-1) then
          nl=ml+1
          exit jloop
       else if (j+nl==n .and. dnl==0) then
          ml=nl-1
          exit jloop
       else if (j+nl==n .and. dnl==1) then
          ml=nl-1; nl=nl+1
          exit jloop
       end if
       nl=nl+dnl
       ml=nl
    end do jloop
    ! Repeatedly eliminate the super diagonal to finish the decomposition.
    qoffs=0; roffs=0
    do j=n-ml,n-1
       ubws(j)=nl-1
       ml=n-j
       ! Note that mq=ml
       pq => q(ubw+2-ml:ubw+1,ubw+2-ml:ubw+1)
       coffs=n-j
       ! apply v_j
       numrots_bv(j)=ml
       do k=1,ml
          rot=rgivens(get_el_bc(b_ub,ubw2,k,coffs+nl-ml+k-1),get_el_bc(b_ub,ubw2,k,coffs+nl-ml+k))
          call tbc_times_rotation(b_ub,n,lbw,ubw2,j,rot,coffs+nl-ml+k-1)
          cs_bv(j,ml-k+1)=rot%cosine; ss_bv(j,ml-k+1)=rot%sine
          k1s_bv(j,ml-k+1)=coffs+nl-ml+k-1; k2s_bv(j,ml-k+1)=coffs+nl-ml+k
          call set_el_bc(b_ub,ubw2,k,coffs+nl-ml+k,(0.0_dp, 0.0_dp))
       end do
       ! downdate
       do k=1,ml-1
          rot=rgivens2(pq(ml,k),pq(ml,k+1))
          call general_times_rotation(pq,rot,k,k+1)
          pq(ml,k)=(0.0_dp, 0.0_dp)
          call rotation_times_tbc(trp_rot(rot),b_ub,n,lbw,ubw2,ml,k)
       end do
       do k=1,nl
          tmp=get_el_bc(b_ub,ubw2,ml,ml+k)
          call set_el_bc(b_ub,ubw2,ml,ml+k,tmp*pq(ml,ml))
       end do
       if (ml==1) then
          exit
       end if
       ! add an extra column to L
       b_ub((ubw2-(n-j)+1)+1:(ubw2-(n-j)+1)+ml-1,n-j)=matmul(transpose(conjg(pq(1:ml-1,1:ml-1))), &
            b_ub((ubw2-(n-j)+1)+1:(ubw2-(n-j)+1)+ml-1,n-j))
       ! apply u_{n-j-1} to Q
       do k=1,numrots_ub(n-j-1)
          rot%cosine=cs_ub(k,n-j-1); rot%sine=ss_ub(k,n-j-1)
          call rotation_times_general(rot,pq,j1s_ub(k,n-j-1),j2s_ub(k,n-j-1))
       end do
       call down_right_shift(pq)
    end do
    ubw_bv=maxval(ubws)
    ! Extract diagonals
    do d=1,lbw+1
       do j=lbw-d+2, n
          b_bv(j,d)=get_el_bc(b_ub,ubw2,j,d+j-lbw-1)
       end do
    end do
    do d=lbw+2,ubw_bv+lbw+1
       do j=1, n-d+lbw+1
          b_bv(j,d)=get_el_bc(b_ub,ubw2,j,d+j-lbw-1)
       end do
    end do
  end subroutine f_c_compress_ub_to_bv



end module compressions_ub_to_bv
