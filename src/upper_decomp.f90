module upper_decomp
  use prec
  use nullvec
  use utility
  use rotation
  use shift
  implicit none
  real(kind=dp), parameter :: eta=1.414_dp
  integer(kind=int32), parameter :: nullmaxits=5
  integer(kind=int32), parameter :: orthmaxits=5
contains
  ! Updating procedure to compute a UB factorization from a general matrix.
  subroutine d_upper_general_to_upper_ub(a, n, ubwmax, tol, b, lbw, ubw, &
       j1s, j2s, cs, ss, numrots, error)
    real(kind=dp), target, dimension(n,n), intent(inout) :: a
    integer(kind=int32), dimension(ubwmax,n), intent(out) :: j1s, j2s
    real(kind=dp), dimension(ubwmax,n), intent(out) :: cs, ss
    real(kind=dp), dimension(ubwmax+lbw+1,n), intent(out) :: b
    real(kind=dp), intent(in) :: tol
    integer(kind=int32), dimension(n), intent(out) :: numrots
    integer(kind=int32), intent(out) :: ubw, error
    integer(kind=int32), intent(in) :: n, ubwmax, lbw
    !
    real(kind=dp), target, dimension(ubwmax+1,n) :: q
    real(kind=dp), dimension(ubwmax+2) :: x
    real(kind=dp) :: nrma, nrmq12
    real(kind=dp), pointer, dimension(:,:) :: pl, pq
    integer(kind=int32) :: i, j, k, nullerr, ortherror, roffs, nl, klast
    integer(kind=int32), dimension(n) :: ubws
    type(d_rotation) :: rot
    !
    q=0.0_dp; numrots=0;
    ss=0.0_dp; cs=0.0_dp; j2s=0; j1s=0
    ubw=0; ubws=0
    error = 0
    nrma = maxabs(a)*sqrt(real(n))
    !
    if (n < 1) then
       error = 1
       return
    end if
    if (n == 1) then
       b(1,1)=a(1,1)
       return
    end if
    ! Compute an initial trivial LQ factorization
    q(1,1:n-1) = a(1,2:n)
    a(1,2:n)=0.0_dp
    a(1,2)=norm2(q(1,1:n-1))
    if (a(1,2) == 0.0_dp) then
       q(1,1:n-1) = 1.0_dp
       q(1,1:n-1)=q(1,1:n-1)/norm2(q(1,1:n-1))
    else
       q(1,1:n-1)=q(1,1:n-1)/a(1,2)
    end if
    nl=1
    klast=n
    kloop: do k=1,n
       ! Current, possibly singular, L should be contained in
       ! a(k-nl+1:k,k+1:k+nl)
       ! At the start of the loop nl is a possible overestimate that
       ! can be decreased by finding a null vector.
       roffs=k-nl
       pl => a(roffs+1:k, k+1:k+nl)
       pq => q(1:nl, 1:n-k)
       if (nl==1) then
          if (abs(pl(1,1)) < tol*nrma) then ! null vector
             ubws(k)=0;  pl(1,1)=0.0_dp
             if (k==n-1) then ! k+nl=n
                klast=k+1
                exit kloop
             else if (k==n-2) then ! k+nl=n-1
                q(1,1)=1.0_dp
                klast=k+1
                exit kloop
             else
                ! set up a 1x1 LQ factorization for computing the next U_k
                ! and start with the next k.
                q(1,1:n-k-1) = a(k+1, k+2:n)
                a(k+1,k+2:n)=0.0_dp
                a(k+1,k+2)=norm2(q(1,1:n-k-1))
                if (a(k+1,k+2) == 0.0_dp) then
                   q(1,1:n-k-1) = 1.0_dp
                   q(1,1:n-k-1)=q(1,1:n-k-1)/norm2(q(1,1:n-k-1))
                else
                   q(1,1:n-k-1)=q(1,1:n-k-1)/a(k+1,k+2)
                end if
                cycle kloop ! next L might be rank deficient; no need to extend it to 2x2.
             end if
          else ! no null vector
             ubws(k)=1
             if (k==n-1) then
                pl(1,1)=pl(1,1)*pq(1,1)
                klast=k+1
                exit kloop
             else
                pl => a(k:k, k+1:k+2) ! extend pl to the right
                pl(1,2) = pl(1,1) ! shift
                pl(1,1)=pl(1,1)*pq(1,1) ! reveal column k+1
                nrmq12 = norm2(pq(1,2:n-k))
                if (nrmq12 == 0.0_dp) then
                   pl(1,2)=0.0_dp;   pq(1,2:n-k)=1.0_dp
                   nrmq12 = norm2(pq(1,2:n-k))
                   pq(1,2:n-k)=pq(1,2:n-k)/nrmq12
                else
                   pq(1,2:n-k)=pq(1,2:n-k)/nrmq12
                   pl(1,2) = pl(1,2)*nrmq12
                end if
                call left_shift(pq)
                if (k < n-2) then
                   pq => q(1:2,1:n-k-1)
                   pq(2,:)=a(k+1,k+2:n)
                   pl => a(k:k+1,k+2:k+3)
                   call d_extend_gs_lq(pq(1:1,:), pl(2,1:1), pl(2,2), pq(2,:), ortherror)
                   if (ortherror == 1) then
                      error = 1; return
                   end if
                   if (k+4 <= n) then
                      a(k+1,k+4:n)=0.0_dp
                   end if
                   nl=nl+1
                else ! k==n-2; only needed to extend L down before exiting.
                   a(k+1,k+2)=a(k+1,k+2)/pq(1,1)
                   klast=k+1
                   exit kloop                   
                end if
             end if
          end if
       else ! nl > 1
          call lower_left_nullvec(x(1:nl),pl,tol*nrma,nullmaxits,nullerr)
          if (nullerr >= 0) then ! if there is a left null vector then introduce a zero row.
             ubws(k)=nl-1
             if (nullerr == 1) then
                pl(1,1)=0.0_dp
                numrots(k)=0
             else if (nullerr > 1) then
                numrots(k)=nullerr-1
                pl(nullerr,nullerr)=0.0_dp
                do j=nullerr-1,1,-1
                   rot=lgivens2(pl(j,j),pl(j+1,j))
                   call rotation_times_general(trp_rot(rot), pl(:,1:j), j,j+1)
                   pl(j,j)=0.0_dp
                   cs(j,k)=rot%cosine; ss(j,k)=rot%sine
                   j1s(j,k)=roffs+j; j2s(j,k)=roffs+j+1
                end do
             else ! nullerr=0
                numrots(k)=nl-1;
                do j=nl,2,-1 ! apply u_k while preserving the triangular structure of L
                   rot=rgivens(x(j-1),x(j))
                   call general_times_rotation(x,rot,j-1,j)
                   call rotation_times_general(trp_rot(rot), pl,j-1,j)
                   cs(j-1,k)=rot%cosine; ss(j-1,k)=rot%sine
                   j1s(j-1,k)=roffs+j-1; j2s(j-1,k)=roffs+j
                   rot=rgivens(pl(j-1,j-1),pl(j-1,j))
                   call general_times_rotation(pl,rot,j-1,j)
                   call rotation_times_general(trp_rot(rot),pq,j-1,j)
                   pl(j-1,j)=0.0_dp
                end do
                pl(1,1)=0.0_dp
             end if
             do j=2,nl ! compress
                rot=rgivens2(pl(j,1),pl(j,j))
                call general_times_rotation(pl(j:nl,:), rot, 1,j)
                call rotation_times_general(trp_rot(rot), pq, 1,j)
                pl(j,1)=0.0_dp
             end do
             if (k+nl==n) then ! square case
                ! reveal column k+1
                do j=nl,2,-1
                   rot=lgivens(pq(1,1),pq(j,1))
                   call rotation_times_general(trp_rot(rot), pq, 1,j)
                   call general_times_rotation(pl(j:nl,:),rot,1,j)
                end do
                pl(:,1)=pl(:,1)*pq(1,1)
                call upper_left_shift(pq)
                ! extend one row
                x(1:nl-1)=0.0_dp
                do j=1,nl-1
                   do i=1,nl-1
                      x(j)=x(j)+a(k+1,k+1+i)*pq(j,i)
                   end do
                end do
                a(k+1,k+2:n)=x(1:nl-1)
                klast=k+1
                exit kloop ! terminate
             else
                pq(1,:)=0.0_dp;     pq(1,1)=1.0_dp
                call d_extend_gs_lq(pq(2:nl,:), x(1:nl-1), x(nl), pq(1,:), ortherror) ! orthogonalize
                if (ortherror == 1) then
                   error = 1; return
                end if
                do j=nl,2,-1
                   rot=lgivens(pq(1,1),pq(j,1))
                   call rotation_times_general(trp_rot(rot), pq, 1,j)
                   call general_times_rotation(pl(j:nl,:),rot,1,j)
                end do
                pl(:,1)=pl(:,1)*pq(1,1)
                call upper_left_shift(pq)
                ! extend the LQ factorization
                pq => q(1:nl,1:n-k-1)
                pq(nl,:)=a(k+1,k+2:n)
                pl => a(roffs+2:k+1,k+2:k+nl+1) ! nl by nl
                call d_extend_gs_lq(pq(1:nl-1,:), pl(nl,1:nl-1), pl(nl,nl), pq(nl,:), ortherror)
                if (ortherror == 1) then
                   error = 1; return
                end if
                if (k+nl+2 <= n) then
                   a(k+1,k+nl+2:n)=0.0_dp
                end if
             end if
          else
             ! no null vector found.  Simply reveal column k+1 if there is room.  Otherwise terminate
             !  with square L.
             ubws(k)=nl
             if (k+nl==n) then
                klast=k ! terminate with square L.
                exit kloop
             else
                pl => a(roffs+1:k, k+1:k+nl+1) ! extend pl to the right. (note this requires k+nl < n)
                call right_shift(pl)
                pq => q(1:nl+1, 1:n-k)
                call down_shift(pq)
                pq(1,:)=0.0_dp
                pq(1,1)=1.0_dp
                ! orthogonalizing.  Note x is used as workspace for coefficients.
                call d_extend_gs_lq(pq(2:nl+1,:), x(1:nl), x(nl+1), pq(1,:), ortherror)
                if (ortherror == 1) then
                   error = 1; return
                end if
                do j=nl+1,2,-1
                   rot=lgivens(pq(1,1),pq(j,1))
                   call rotation_times_general(trp_rot(rot), pq, 1,j)
                   call general_times_rotation(pl(j-1:nl,:),rot,1,j)
                end do
                pl(:,1)=pl(:,1)*pq(1,1)
                call upper_left_shift(pq)
                if (k+nl==n-1) then ! q is now nl by nl.  Extend the LQ factorization down one row before stopping.
                   pq => q(1:nl,1:n-k-1)
                   x(1:nl)=0.0_dp
                   do j=1,nl
                      do i=1,nl
                         x(j)=x(j)+a(k+1,k+1+i)*pq(j,i)
                      end do
                   end do
                   a(k+1,k+2:n)=x(1:nl)
                   klast=k+1
                   exit kloop
                else ! q is not square.  Make L (nl+1)x(nl+1)
                   pq => q(1:nl+1,1:n-k-1)
                   pq(nl+1,:)=a(k+1,k+2:n)
                   pl => a(roffs+1:k+1,k+2:k+nl+2)
                   call d_extend_gs_lq(pq(1:nl,:), pl(nl+1,1:nl), pl(nl+1,nl+1), pq(nl+1,:), ortherror)
                   if (ortherror == 1) then
                      error = 1; return
                   end if
                   if (k+nl+3 <= n) then
                      a(k+1,k+nl+3:n)=0.0_dp
                   end if
                   nl=nl+1
                end if
             end if
          end if ! null vector check
       end if
    end do kloop
    ! If klast = k+1, we need to terminate on an nl+1 by nl matrix L
    ! contained in a(klast-(n-klast):klast,klast+1:n).  If klast=k, then L is square and
    ! contained in a(k-(n-k)+1:k,k+1:n).  The former happens when the last step
    ! of kloop found a null vector.  The latter happens when it didn't.
    if (klast == k) then ! square termination
       ubws(k:n-1)=n-k
       nl=n-k
       if (n-k > 1) then
          do i=1,n-k-1 ! nl-1
             nl=n-k-i+1
             roffs=k-nl
             pl => a(roffs+1:k,k+i:n)
             pq => q(i:n-k,i:n-k)
             ! reveal column k+i
             do j=nl,2,-1
                rot=lgivens(pq(j-1,1),pq(j,1))
                call rotation_times_general(trp_rot(rot), pq, j-1,j)
                call general_times_rotation(pl, rot, j-1,j)
             end do
             ! triangularize L
             pl(:,1)=pl(:,1)*pq(1,1)
             numrots(k+i)=nl-1
             do j=nl-1,1,-1
                rot=lgivens2(pl(j,j+1),pl(j+1,j+1))
                call rotation_times_general(trp_rot(rot),pl(:,2:nl),j,j+1)
                cs(j,k+i)=rot%cosine; ss(j,k+i)=rot%sine
                j1s(j,k+i)=roffs+j; j2s(j,k+i)=roffs+j+1
                pl(j,j+1)=0.0_dp
             end do
          end do
          pl(2,2)=pl(2,2)*pq(2,2)
       end if
    else ! rectangular termination
       k=k+1
       ubws(k:n-1)=n-k
       nl=n-k
       if (nl > 0) then
          do i=1,n-k
             nl=n-k-i+1
             roffs=k-nl-1
             pl=>a(roffs+1:k,k+i:n)
             pq=>q(i:n-k,i:n-k)
             numrots(k+i-1)=nl
             ! Triangularize L
             do j=nl,1,-1
                rot=lgivens2(pl(j,j),pl(j+1,j))
                cs(j,k+i-1)=rot%cosine; ss(j,k+i-1)=rot%sine
                j1s(j,k+i-1)=roffs+j; j2s(j,k+i-1)=roffs+j+1
                call rotation_times_general(trp_rot(rot), pl,j,j+1)
                pl(j,j)=0.0_dp
             end do
             ! reveal column k + i
             if (nl > 1) then
                do j=nl,2,-1
                   rot=lgivens(pq(j-1,1),pq(j,1))
                   call rotation_times_general(trp_rot(rot), pq, j-1,j)
                   call general_times_rotation(pl, rot, j-1,j)
                end do
             end if
             pl(:,1)=pl(:,1)*pq(1,1)
          end do
       end if
    end if
    call d_extract_diagonals_upper_ub(a,n,ubwmax,lbw,b,ubws,ubw)
  end subroutine d_upper_general_to_upper_ub

  subroutine d_extract_diagonals_upper_ub(a, n, ubwmax, lbw, b, ubws, ubw)
    real(kind=dp), target, dimension(n,n), intent(in) :: a
    real(kind=dp), dimension(ubwmax+lbw+1,n), intent(out) :: b
    integer(kind=int32), intent(out) :: ubw
    integer(kind=int32), intent(in) :: n, ubwmax, lbw
    integer(kind=int32), dimension(n), intent(in) :: ubws
    !
    integer(kind=int32) :: k
    ubw=maxval(ubws)
    ! put diagonals in b
    b=0.0_dp
    do k = 1,ubw
       b(ubw+2-k:ubw+lbw+1,k)=a(1:lbw+k,k)
    end do
    do k = 1,lbw
       b(1:ubw+k,n-k+1)=a(n-ubw-k+1:n,n-k+1)
    end do
    do k=ubw+1,n-lbw
       b(1:ubw+lbw+1,k)=a(k-ubw:k+lbw,k)
    end do
  end subroutine d_extract_diagonals_upper_ub


  ! Extend a Gram-Schmidt LQ decomposition
  subroutine d_extend_gs_lq(q1, l, rho, q2, error)
    ! orthogonalize q2 against q1.  Return error = 1 if
    ! too many iterations are required.
    real(kind=dp), dimension (:,:), intent(in) :: q1
    real(kind=dp), intent(out) :: rho
    real(kind=dp), dimension(:), intent(out) :: l
    real(kind=dp), dimension(:), intent(inout) :: q2
    integer(kind=int32), intent(out) :: error
    !
    real(kind=dp), dimension(size(l)) :: x
    real(kind=dp) :: nrm1, nrm2
    integer(kind=int32) :: k
    !
    nrm1=norm2(q2)
    l = matmul(q1,q2)
    q2=q2-matmul(l,q1)
    nrm2=norm2(q2)
    rho=nrm2
    ! reorthogonalize as needed
    error=0
    k=0
    do while (nrm1 > eta * nrm2 .and. nrm2 > 0 .and. k < orthmaxits)
       k=k+1
       x=matmul(q1,q2)
       l=l+x
       nrm1=norm2(q2)
       q2=q2-matmul(x,q1)
       nrm2=norm2(q2)
       rho=nrm2
    end do
    do while (nrm2 == 0 .and. k < orthmaxits)
       ! got the zero vector.  Orthogonalize something random.
       k=0
       call random_number(q2)
       nrm2 = 1.0_dp
       nrm1=2.0_dp*eta
       do while (nrm1 > eta*nrm2 .and. nrm2 > 0.0_dp .and. k < orthmaxits)
          k=k+1
          nrm1=norm2(q2)
          x=matmul(q1,q2)
          q2=q2-matmul(x,q1)
          nrm2=norm2(q2)
       end do
    end do
    q2=q2/nrm2
    if (k >= orthmaxits) then
       error = 1
    end if
  end subroutine d_extend_gs_lq

  subroutine d_form_upper_ub(a, mb, n, b, lbw, ubw, j1s, j2s, cs, ss, numrots)
    real(kind=dp), target, dimension(n,n), intent(out) :: a
    integer(kind=int32), dimension(mb-lbw-1,n), intent(in) :: j1s, j2s
    real(kind=dp), dimension(mb-lbw-1,n), intent(in) :: cs, ss
    real(kind=dp), dimension(mb,n), intent(in) :: b
    integer(kind=int32), dimension(n), intent(in) :: numrots
    integer(kind=int32), intent(in) :: ubw, lbw, n, mb
    !
    integer(kind=int32) :: j,k
    type(d_rotation) :: rot
    a=0.0_dp
    do k=n-1,n-lbw,-1
       a(k+1-ubw:n,k+1)=b(1:n-k+ubw,k+1)
       do j=1,numrots(k)
          rot%cosine=cs(j,k); rot%sine=ss(j,k)
          call rotation_times_general(rot,a(:,k+1:n),j1s(j,k),j2s(j,k))
       end do
    end do
    do k=n-lbw-1, ubw,-1
       a(k+1-ubw:k+1+lbw,k+1)=b(1:lbw+ubw+1,k+1)
       do j=1,numrots(k)
          rot%cosine=cs(j,k); rot%sine=ss(j,k)
          call rotation_times_general(rot,a(:,k+1:n),j1s(j,k),j2s(j,k))
       end do
    end do
    do k=ubw-1,1,-1
       a(1:lbw+k+1,k+1)=b(ubw-k+1:ubw+lbw+1,k+1)
       do j=1,numrots(k)
          rot%cosine=cs(j,k); rot%sine=ss(j,k)
          call rotation_times_general(rot,a(:,k+1:n),j1s(j,k),j2s(j,k))
       end do
    end do
    a(1:lbw+1,1)=b(ubw+1:ubw+lbw+1,1)
  end subroutine d_form_upper_ub

  !
  ! Complex.
  !
  ! Updating procedure to compute a UB factorization from a general matrix.
  !
  subroutine c_upper_general_to_upper_ub(a, n, ubwmax, tol, b, lbw, ubw, &
       j1s, j2s, cs, ss, numrots, error)
    complex(kind=dp), target, dimension(n,n), intent(inout) :: a
    integer(kind=int32), dimension(ubwmax,n), intent(out) :: j1s, j2s
    complex(kind=dp), dimension(ubwmax,n), intent(out) :: cs, ss
    complex(kind=dp), dimension(ubwmax+lbw+1,n), intent(out) :: b
    real(kind=dp), intent(in) :: tol
    integer(kind=int32), dimension(n), intent(out) :: numrots
    integer(kind=int32), intent(out) :: ubw, error
    integer(kind=int32), intent(in) :: n, ubwmax, lbw
    !
    complex(kind=dp), target, dimension(ubwmax+1,n) :: q
    complex(kind=dp), dimension(ubwmax+2) :: x
    real(kind=dp) :: nrma, nrmq12
    complex(kind=dp), pointer, dimension(:,:) :: pl, pq
    integer(kind=int32) :: i, j, k, nullerr, ortherror, roffs, nl, klast
    integer(kind=int32), dimension(n) :: ubws
    type(c_rotation) :: rot
    !
    q=(0.0_dp, 0.0_dp); numrots=0;
    ss=(0.0_dp, 0.0_dp); cs=(0.0_dp, 0.0_dp); j2s=0; j1s=0
    ubw=0; ubws=0
    error = 0
    nrma = maxabs(a)*sqrt(real(n))
    !
    if (n < 1) then
       error = 1
       return
    end if
    if (n == 1) then
       b(1,1)=a(1,1)
       return
    end if
    ! Compute an initial trivial LQ factorization
    q(1,1:n-1) = a(1,2:n)
    a(1,2:n)=(0.0_dp, 0.0_dp)
    a(1,2)=norm2(q(1,1:n-1))
    if (a(1,2) == (0.0_dp, 0.0_dp)) then
       q(1,1:n-1) = (1.0_dp,0.0_dp)
       q(1,1:n-1)=q(1,1:n-1)/norm2(q(1,1:n-1))
    else
       q(1,1:n-1)=q(1,1:n-1)/a(1,2)
    end if
    nl=1
    klast=n
    kloop: do k=1,n
       ! Current, possibly singular, L should be contained in
       ! a(k-nl+1:k,k+1:k+nl)
       ! At the start of the loop nl is a possible overestimate that
       ! can be decreased by finding a null vector.
       roffs=k-nl
       pl => a(roffs+1:k, k+1:k+nl)
       pq => q(1:nl, 1:n-k)
       if (nl==1) then
          if (abs(pl(1,1)) < tol*nrma) then ! null vector
             ubws(k)=0;  pl(1,1)=(0.0_dp, 0.0_dp)
             if (k==n-1) then ! k+nl=n
                klast=k+1
                exit kloop
             else if (k==n-2) then ! k+nl=n-1
                q(1,1)=(1.0_dp,0.0_dp)
                klast=k+1
                exit kloop
             else
                ! set up a 1x1 LQ factorization for computing the next U_k
                ! and start with the next k.
                q(1,1:n-k-1) = a(k+1, k+2:n)
                a(k+1,k+2:n)=(0.0_dp, 0.0_dp)
                a(k+1,k+2)=norm2(q(1,1:n-k-1))
                if (a(k+1,k+2) == (0.0_dp, 0.0_dp)) then
                   q(1,1:n-k-1) = (1.0_dp,0.0_dp)
                   q(1,1:n-k-1)=q(1,1:n-k-1)/norm2(q(1,1:n-k-1))
                else
                   q(1,1:n-k-1)=q(1,1:n-k-1)/a(k+1,k+2)
                end if
                cycle kloop ! next L might be rank deficient; no need to extend it to 2x2.
             end if
          else ! no null vector
             ubws(k)=1
             if (k==n-1) then
                pl(1,1)=pl(1,1)*pq(1,1)
                klast=k+1
                exit kloop
             else
                pl => a(k:k, k+1:k+2) ! extend pl to the right
                pl(1,2) = pl(1,1) ! shift
                pl(1,1)=pl(1,1)*pq(1,1) ! reveal column k+1
                nrmq12 = norm2(pq(1,2:n-k))
                if (nrmq12 == 0.0_dp) then
                   pl(1,2)=(0.0_dp, 0.0_dp);   pq(1,2:n-k)=(1.0_dp,0.0_dp)
                   nrmq12 = norm2(pq(1,2:n-k))
                   pq(1,2:n-k)=pq(1,2:n-k)/nrmq12
                else
                   pq(1,2:n-k)=pq(1,2:n-k)/nrmq12
                   pl(1,2) = pl(1,2)*nrmq12
                end if
                call left_shift(pq)
                if (k < n-2) then
                   pq => q(1:2,1:n-k-1)
                   pq(2,:)=a(k+1,k+2:n)
                   pl => a(k:k+1,k+2:k+3)
                   call c_extend_gs_lq(pq(1:1,:), pl(2,1:1), pl(2,2), pq(2,:), ortherror)
                   if (ortherror == 1) then
                      error = 1; return
                   end if
                   if (k+4 <= n) then
                      a(k+1,k+4:n)=(0.0_dp, 0.0_dp)
                   end if
                   nl=nl+1
                else ! k==n-2; only needed to extend L down before exiting.
                   a(k+1,k+2)=a(k+1,k+2)/pq(1,1)
                   klast=k+1
                   exit kloop                   
                end if
             end if
          end if
       else ! nl > 1
          call lower_left_nullvec(x(1:nl),pl,tol*nrma,nullmaxits,nullerr)
          if (nullerr >= 0) then ! if there is a left null vector then introduce a zero row.
             ubws(k)=nl-1
             if (nullerr == 1) then
                pl(1,1)=(0.0_dp, 0.0_dp)
                numrots(k)=0
             else if (nullerr > 1) then
                numrots(k)=nullerr-1
                pl(nullerr,nullerr)=(0.0_dp, 0.0_dp)
                do j=nullerr-1,1,-1
                   rot=lgivens2(pl(j,j),pl(j+1,j))
                   call rotation_times_general(trp_rot(rot), pl(:,1:j), j,j+1)
                   pl(j,j)=(0.0_dp, 0.0_dp)
                   cs(j,k)=rot%cosine; ss(j,k)=rot%sine
                   j1s(j,k)=roffs+j; j2s(j,k)=roffs+j+1
                end do
             else ! nullerr=0
                numrots(k)=nl-1;
                do j=nl,2,-1 ! apply u_k while preserving the triangular structure of L
                   rot=rgivens(x(j-1),x(j))
                   call general_times_rotation(x,rot,j-1,j)
                   call rotation_times_general(trp_rot(rot), pl,j-1,j)
                   cs(j-1,k)=rot%cosine; ss(j-1,k)=rot%sine
                   j1s(j-1,k)=roffs+j-1; j2s(j-1,k)=roffs+j
                   rot=rgivens(pl(j-1,j-1),pl(j-1,j))
                   call general_times_rotation(pl,rot,j-1,j)
                   call rotation_times_general(trp_rot(rot),pq,j-1,j)
                   pl(j-1,j)=(0.0_dp, 0.0_dp)
                end do
                pl(1,1)=(0.0_dp, 0.0_dp)
             end if
             do j=2,nl ! compress
                rot=rgivens2(pl(j,1),pl(j,j))
                call general_times_rotation(pl(j:nl,:), rot, 1,j)
                call rotation_times_general(trp_rot(rot), pq, 1,j)
                pl(j,1)=(0.0_dp, 0.0_dp)
             end do
             if (k+nl==n) then ! square case
                ! reveal column k+1
                do j=nl,2,-1
                   rot=lgivens(pq(1,1),pq(j,1))
                   call rotation_times_general(trp_rot(rot), pq, 1,j)
                   call general_times_rotation(pl(j:nl,:),rot,1,j)
                end do
                pl(:,1)=pl(:,1)*pq(1,1)
                call upper_left_shift(pq)
                ! extend one row
                x(1:nl-1)=(0.0_dp, 0.0_dp)
                do j=1,nl-1
                   do i=1,nl-1
                      x(j)=x(j)+a(k+1,k+1+i)*conjg(pq(j,i))
                   end do
                end do
                a(k+1,k+2:n)=x(1:nl-1)
                klast=k+1
                exit kloop ! terminate
             else
                pq(1,:)=(0.0_dp, 0.0_dp);     pq(1,1)=(1.0_dp,0.0_dp)
                call c_extend_gs_lq(pq(2:nl,:), x(1:nl-1), x(nl), pq(1,:), ortherror) ! orthogonalize
                if (ortherror == 1) then
                   error = 1; return
                end if
                do j=nl,2,-1
                   rot=lgivens(pq(1,1),pq(j,1))
                   call rotation_times_general(trp_rot(rot), pq, 1,j)
                   call general_times_rotation(pl(j:nl,:),rot,1,j)
                end do
                pl(:,1)=pl(:,1)*pq(1,1)
                call upper_left_shift(pq)
                ! extend the LQ factorization
                pq => q(1:nl,1:n-k-1)
                pq(nl,:)=a(k+1,k+2:n)
                pl => a(roffs+2:k+1,k+2:k+nl+1) ! nl by nl
                call c_extend_gs_lq(pq(1:nl-1,:), pl(nl,1:nl-1), pl(nl,nl), pq(nl,:), ortherror)
                if (ortherror == 1) then
                   error = 1; return
                end if
                if (k+nl+2 <= n) then
                   a(k+1,k+nl+2:n)=(0.0_dp, 0.0_dp)
                end if
             end if
          else
             ! no null vector found.  Simply reveal column k+1 if there is room.  Otherwise terminate
             !  with square L.
             ubws(k)=nl
             if (k+nl==n) then
                klast=k ! terminate with square L.
                exit kloop
             else
                pl => a(roffs+1:k, k+1:k+nl+1) ! extend pl to the right. (note this requires k+nl < n)
                call right_shift(pl)
                pq => q(1:nl+1, 1:n-k)
                call down_shift(pq)
                pq(1,:)=(0.0_dp, 0.0_dp)
                pq(1,1)=(1.0_dp, 0.0_dp)
                ! orthogonalizing.  Note x is used as workspace for coefficients.
                call c_extend_gs_lq(pq(2:nl+1,:), x(1:nl), x(nl+1), pq(1,:), ortherror)
                if (ortherror == 1) then
                   error = 1; return
                end if
                do j=nl+1,2,-1
                   rot=lgivens(pq(1,1),pq(j,1))
                   call rotation_times_general(trp_rot(rot), pq, 1,j)
                   call general_times_rotation(pl(j-1:nl,:),rot,1,j)
                end do
                pl(:,1)=pl(:,1)*pq(1,1)
                call upper_left_shift(pq)
                if (k+nl==n-1) then ! q is now nl by nl.  Extend the LQ factorization down one row before stopping.
                   pq => q(1:nl,1:n-k-1)
                   x(1:nl)=(0.0_dp, 0.0_dp)
                   do j=1,nl
                      do i=1,nl
                         x(j)=x(j)+a(k+1,k+1+i)*conjg(pq(j,i))
                      end do
                   end do
                   a(k+1,k+2:n)=x(1:nl)
                   klast=k+1
                   exit kloop
                else ! q is not square.  Make L (nl+1)x(nl+1)
                   pq => q(1:nl+1,1:n-k-1)
                   pq(nl+1,:)=a(k+1,k+2:n)
                   pl => a(roffs+1:k+1,k+2:k+nl+2)
                   call c_extend_gs_lq(pq(1:nl,:), pl(nl+1,1:nl), pl(nl+1,nl+1), pq(nl+1,:), ortherror)
                   if (ortherror == 1) then
                      error = 1; return
                   end if
                   if (k+nl+3 <= n) then
                      a(k+1,k+nl+3:n)=(0.0_dp, 0.0_dp)
                   end if
                   nl=nl+1
                end if
             end if
          end if ! null vector check
       end if
    end do kloop
    ! If klast = k+1, we need to terminate on an nl+1 by nl matrix L
    ! contained in a(klast-(n-klast):klast,klast+1:n).  If klast=k, then L is square and
    ! contained in a(k-(n-k)+1:k,k+1:n).  The former happens when the last step
    ! of kloop found a null vector.  The latter happens when it didn't.
    if (klast == k) then ! square termination
       ubws(k:n-1)=n-k
       nl=n-k
       if (n-k > 1) then
          do i=1,n-k-1 ! nl-1
             nl=n-k-i+1
             roffs=k-nl
             pl => a(roffs+1:k,k+i:n)
             pq => q(i:n-k,i:n-k)
             ! reveal column k+i
             do j=nl,2,-1
                rot=lgivens(pq(j-1,1),pq(j,1))
                call rotation_times_general(trp_rot(rot), pq, j-1,j)
                call general_times_rotation(pl, rot, j-1,j)
             end do
             ! triangularize L
             pl(:,1)=pl(:,1)*pq(1,1)
             numrots(k+i)=nl-1
             do j=nl-1,1,-1
                rot=lgivens2(pl(j,j+1),pl(j+1,j+1))
                call rotation_times_general(trp_rot(rot),pl(:,2:nl),j,j+1)
                cs(j,k+i)=rot%cosine; ss(j,k+i)=rot%sine
                j1s(j,k+i)=roffs+j; j2s(j,k+i)=roffs+j+1
                pl(j,j+1)=(0.0_dp, 0.0_dp)
             end do
          end do
          pl(2,2)=pl(2,2)*pq(2,2)
       end if
    else ! rectangular termination
       k=k+1
       ubws(k:n-1)=n-k
       nl=n-k
       if (nl > 0) then
          do i=1,n-k
             nl=n-k-i+1
             roffs=k-nl-1
             pl=>a(roffs+1:k,k+i:n)
             pq=>q(i:n-k,i:n-k)
             numrots(k+i-1)=nl
             ! Triangularize L
             do j=nl,1,-1
                rot=lgivens2(pl(j,j),pl(j+1,j))
                cs(j,k+i-1)=rot%cosine; ss(j,k+i-1)=rot%sine
                j1s(j,k+i-1)=roffs+j; j2s(j,k+i-1)=roffs+j+1
                call rotation_times_general(trp_rot(rot), pl,j,j+1)
                pl(j,j)=(0.0_dp, 0.0_dp)
             end do
             ! reveal column k + i
             if (nl > 1) then
                do j=nl,2,-1
                   rot=lgivens(pq(j-1,1),pq(j,1))
                   call rotation_times_general(trp_rot(rot), pq, j-1,j)
                   call general_times_rotation(pl, rot, j-1,j)
                end do
             end if
             pl(:,1)=pl(:,1)*pq(1,1)
          end do
       end if
    end if
    call c_extract_diagonals_upper_ub(a,n,ubwmax,lbw,b,ubws,ubw)
  end subroutine c_upper_general_to_upper_ub

  subroutine c_extract_diagonals_upper_ub(a, n, ubwmax, lbw, b, ubws, ubw)
    complex(kind=dp), target, dimension(n,n), intent(in) :: a
    complex(kind=dp), dimension(ubwmax+lbw+1,n), intent(out) :: b
    integer(kind=int32), intent(out) :: ubw
    integer(kind=int32), intent(in) :: n, ubwmax, lbw
    integer(kind=int32), dimension(n), intent(in) :: ubws
    !
    integer(kind=int32) :: k
    ubw=maxval(ubws)
    ! put diagonals in b
    b=(0.0_dp, 0.0_dp)
    do k = 1,ubw
       b(ubw+2-k:ubw+lbw+1,k)=a(1:lbw+k,k)
    end do
    do k = 1,lbw
       b(1:ubw+k,n-k+1)=a(n-ubw-k+1:n,n-k+1)
    end do
    do k=ubw+1,n-lbw
       b(1:ubw+lbw+1,k)=a(k-ubw:k+lbw,k)
    end do
  end subroutine c_extract_diagonals_upper_ub


  ! Extend a Gram-Schmidt LQ decomposition
  subroutine c_extend_gs_lq(q1, l, rho, q2, error)
    ! orthogonalize q2 against q1.  Return error = 1 if
    ! too many iterations are required.
    complex(kind=dp), dimension (:,:), intent(in) :: q1
    complex(kind=dp), intent(out) :: rho
    complex(kind=dp), dimension(:), intent(out) :: l
    complex(kind=dp), dimension(:), intent(inout) :: q2
    integer(kind=int32), intent(out) :: error
    !
    complex(kind=dp), dimension(size(l)) :: x
    real(kind=dp) :: nrm1, nrm2
    integer(kind=int32) :: k
    !
    nrm1=norm2(q2)
    l = matmul(q1,conjg(q2))
    l=conjg(l)
    q2=q2-matmul(l,q1)
    nrm2=norm2(q2)
    rho=nrm2
    ! reorthogonalize as needed
    error=0
    k=0
    do while (nrm1 > eta * nrm2 .and. nrm2 > 0 .and. k < orthmaxits)
       k=k+1
       x=matmul(q1,conjg(q2))
       x=conjg(x)
       l=l+x
       nrm1=norm2(q2)
       q2=q2-matmul(x,q1)
       nrm2=norm2(q2)
       rho=nrm2
    end do
    do while (nrm2 == 0 .and. k < orthmaxits)
       ! got the zero vector.  Orthogonalize something random.
       k=0
       call random_complex(q2)
       nrm2 = 1.0_dp
       nrm1=2.0_dp*eta
       do while (nrm1 > eta*nrm2 .and. nrm2 > 0.0_dp .and. k < orthmaxits)
          k=k+1
          nrm1=norm2(q2)
          x=matmul(q1,conjg(q2))
          x=conjg(x)
          q2=q2-matmul(x,q1)
          nrm2=norm2(q2)
       end do
    end do
    q2=q2/nrm2
    if (k >= orthmaxits) then
       error = 1
    end if
  end subroutine c_extend_gs_lq

  subroutine c_form_upper_ub(a, mb, n, b, lbw, ubw, j1s, j2s, cs, ss, numrots)
    complex(kind=dp), target, dimension(n,n), intent(out) :: a
    integer(kind=int32), dimension(mb-lbw-1,n), intent(in) :: j1s, j2s
    complex(kind=dp), dimension(mb-lbw-1,n), intent(in) :: cs, ss
    complex(kind=dp), dimension(mb,n), intent(in) :: b
    integer(kind=int32), dimension(n), intent(in) :: numrots
    integer(kind=int32), intent(in) :: ubw, lbw, n, mb
    !
    integer(kind=int32) :: j,k
    type(c_rotation) :: rot
    a=(0.0_dp, 0.0_dp)
    do k=n-1,n-lbw,-1
       a(k+1-ubw:n,k+1)=b(1:n-k+ubw,k+1)
       do j=1,numrots(k)
          rot%cosine=cs(j,k); rot%sine=ss(j,k)
          call rotation_times_general(rot,a(:,k+1:n),j1s(j,k),j2s(j,k))
       end do
    end do
    do k=n-lbw-1, ubw,-1
       a(k+1-ubw:k+1+lbw,k+1)=b(1:lbw+ubw+1,k+1)
       do j=1,numrots(k)
          rot%cosine=cs(j,k); rot%sine=ss(j,k)
          call rotation_times_general(rot,a(:,k+1:n),j1s(j,k),j2s(j,k))
       end do
    end do
    do k=ubw-1,1,-1
       a(1:lbw+k+1,k+1)=b(ubw-k+1:ubw+lbw+1,k+1)
       do j=1,numrots(k)
          rot%cosine=cs(j,k); rot%sine=ss(j,k)
          call rotation_times_general(rot,a(:,k+1:n),j1s(j,k),j2s(j,k))
       end do
    end do
    a(1:lbw+1,1)=b(ubw+1:ubw+lbw+1,1)
  end subroutine c_form_upper_ub

end module upper_decomp
