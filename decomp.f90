module decomp
  use prec
  use nullvec
  use utility
  use rotation
  implicit none
  real(kind=dp), parameter :: eta=1.414_dp
  integer(kind=int32), parameter :: nullmaxits=5
  integer(kind=int32), parameter :: orthmaxits=5
contains
  ! Updating procedure to compute a UB factorization from a general matrix.
  subroutine d_upper_general_to_upper_ub(a, n, ubwmax, lbw, tol, j1s, j2s, &
       cs, ss, numrots, b, ubw, error)
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
    integer(kind=int32) :: i, j, k, nullerr, ortherror
    integer(kind=int32), dimension(n) :: ubws
    type(d_rotation) :: rot
    !
    q=0.0_dp; numrots=0;
    ss=0.0_dp; cs=0.0_dp; j2s=0; j1s=0
    ubw=0; ubws=0
    error = 0
    nrma = normf(a)
    !
    if (n < 1) then
       error = 1
       return
    end if
    if (n == 1) then
       b(1,1)=a(1,1)
       return
    end if
    ! Store the first column unmodified.
    b(ubwmax+1:ubwmax+lbw+1,1) = a(1:lbw+1,1)
    ! Compute an initial trivial LQ factorization
    q(1,1:n-1) = a(1,2:n)
    a(1,2:n)=0.0_dp
    a(1,2)=norm2(q(1,2:n))
    if (a(1,2) == 0.0_dp) then
       q(1,1:n-1) = 1.0_dp
       q(1,1:n-1)=q(1,1:n-1)/norm2(q(1,1:n-1))
    else
       q(1,1:n-1)=q(1,1:n-1)/a(1,2)
    end if
    ubw=1
    kloop: do k=1,n
       ! Current, possibly singular, L should be contained in
       ! a(k-ubw+1:k,k+1:k+ubw)
       ! At the start of the loop ubw is a possible overestimate that
       ! can be decreased by finding a null vector.
       ! At the start of the loop we require k+ubw < n and ubw > 0
       pl => a(k-ubw+1:k, k+1:k+ubw)
       pq => q(1:ubw, 1:n-k)
       if (k+ubw == n) then
          exit kloop
       end if
       if (ubw==1) then
          if (abs(pl(1,1)) < tol*nrma) then ! null vector
             ubws(k)=0;  ubw=0;  pl(1,1)=0.0_dp
             if (k < n-1) then ! At least one more iteration required (otherwise k+1+ubw == n.)
                ! set up a 1x1 LQ factorization for computing the next U_k
                ! and start with the next k.
                q(1,1:n-k-1) = a(k+1, k+2:n)
                a(k+1,k+2:n)=0.0_dp
                a(k+1,k+2)=norm2(q(1,1:n-k-1))
                if (a(k+1,k+2) == 0) then
                   q(1,1:n-k-1) = 1.0_dp
                   q(1,1:n-k-1)=q(1,1:n-k-1)/norm2(q(1,1:n-k-1))
                else
                   q(1,1:n-k-1)=q(1,1:n-k-1)/a(k+1,k+2)
                end if
                ubw=1
             end if
             cycle kloop ! next L might be rank deficient; no need to extend it to 2x2.
          else ! no null vector
             ubws(k)=1
             if (k < n-1) then ! reveal column k+1 and shift pl to the right
                pl => a(k:k, k+1:k+2) ! extend pl to the right
                pl(1,2) = pl(1,1) ! shift
                pl(1,1)=0.0_dp
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
                pl => a(k:k, k+2:k+2)
                pq => q(1:1,1:n-k-1) ! If k+1+ubw < n (k<n-2), these will be extended below to make L 2x2.
             else ! Simply compute a(n-1,n)
                ubw=0 ! L no longer exists after this iteration.
                a(n-1,n)=pl(1,1)*pq(1,1)
                cycle kloop ! The loop terminates after this since k+1+ubw==n
             end if
          end if
       else ! (ubw > 1 and k+ubw < n (k <= n-2))
          call left_nullvec(x(1:ubw),pl,tol*nrma,nullmaxits,nullerr)
          if (nullerr == 0) then ! if there is a left null vector then introduce a zero row.
             numrots(k)=ubw-1;   ubws(k)=ubw-1
             do j=ubw,2,-1 ! apply u_k while preserving the triangular structure of L
                rot=rgivens(x(j-1),x(j))
                call general_times_rotation(x,rot,j-1,j)
                call rotation_times_general(trp_rot(rot), pl,j-1,j)
                cs(j-1,k)=rot%cosine; ss(j-1,k)=rot%sine
                j1s(j-1,k)=j-1+k-ubw; j2s(j-1,k)=j+k-ubw
                rot=rgivens(pl(j-1,j-1),pl(j-1,j))
                call general_times_rotation(pl,rot,j-1,j)
                call rotation_times_general(trp_rot(rot),pq,j-1,j)
                pl(j-1,j)=0.0_dp
             end do
             pl(1,1)=0.0_dp
             do j=2,ubw ! compress
                rot=rgivens2(pl(j,1),pl(j,j))
                call general_times_rotation(pl(j:ubw,:), rot, 1,j)
                call rotation_times_general(trp_rot(rot), pq, 1,j)
                pl(j,1)=0.0_dp
             end do
             ! reveal column k+1
             pq(1,:)=0.0_dp;     pq(1,1)=1.0_dp
             call d_extend_gs(pq(2:ubw,:), x(1:ubw), x(ubw+1), pq(1,:), ortherror) ! orthogonalize
             if (ortherror == 1) then
                error = 1; return
             end if
             do j=ubw,2,-1
                rot=lgivens(pq(1,1),pq(j,1))
                call rotation_times_general(trp_rot(rot), pq, 1,j)
                call general_times_rotation(pl(j:ubw,:),rot,1,j)
             end do
             pl(:,1)=pl(:,1)*pq(1,1)
             call upper_left_shift(pq)
             ubw=ubw-1
             pl => a(k-ubw+1:k,k+2:k+ubw+1)
             pq => q(1:ubw-1,1:n-k-1)
          else ! no null vector found.  Simply reveal column k+1
             ubws(k)=ubw
             pl => a(k-ubw+1:k, k+1:k+ubw+1) ! move pl to the right in (note this requires k+ubw+1 <= n)
             call right_shift(pl)
             pq => q(1:ubw+1, 1:n-k)
             call down_shift(pq)
             pq(1,:)=0.0_dp
             pq(1,1)=1.0_dp
             ! orthogonalizing.  Note x is used as workspace for coefficients.
             call d_extend_gs(pq(2:ubw+1,:), x(1:ubw+1), x(ubw+2), pq(1,:), ortherror)
             if (ortherror == 1) then
                error = 1; return
             end if
             do j=ubw+1,2,-1
                rot=lgivens(pq(1,1),pq(j,1))
                call rotation_times_general(trp_rot(rot), pq, 1,j)
                call general_times_rotation(pl(j:ubw,:),rot,1,j)
             end do
             pl(:,1)=pl(:,1)*pq(1,1)
             call upper_left_shift(pq)
             pl => a(k-ubw+1:k, k+2:k+ubw+1)
             pq => q(1:ubw,1:n-k-1)
          end if ! null vector check
       end if
       ! Currently pl => a(k-ubw+1:k, k+2:k+ubw+1).
       ! At this point ubw >= 1.  If a null vector was found, ubw was decreased by 1.
       ! If no null vector was found, ubw stayed the same.  Either way, L needs to be expanded
       ! unless it has reached the last column of a (k+ubw+1==n).
       if (k+1+ubw < n) then ! if k+1+ubw==n, termination at the top of the next iteration.
          pq => q(1:ubw+1,1:n-k-1) ! There is room to expand l.
          pq(ubw+1,:)=a(k+1,k+2:n)
          pl => a(k-ubw+1:k+1,k+2:k+ubw+2)
          call d_extend_gs(pq(1:ubw,:), pl(ubw+1,1:ubw), pl(ubw+1,ubw+1), pq(ubw+1,:), ortherror)
          if (ortherror == 1) then
             error = 1; return
          end if
          if (k+ubw+3 <= n) then
             a(k+1,k+ubw+3:n)=0.0_dp
          end if
          ubw=ubw+1
       else if (ubw > 0) then
          ! k+1+ubw = n implies pl and pq are square and ubw x ubw.  The
          ! iteration will terminate.  Add an extra row to pl to make
          ! pl ubw+1 x ubw.
          x(1:ubw)=0.0_dp
          do j=1,ubw
             do i=1,ubw
                x(j)=x(j)+a(k+1,k+1+i)*pq(i,j)
             end do
          end do
          a(k+1,k+2:n)=x(1:ubw)
       end if
    end do kloop
    ! kloop has terminated with k+ubw==n.  The LQ factorization
    ! has pl => a(k-ubw:k, k+1:n) and
    !     pq => q(1:ubw, 1:ubw)
    ! pl is ubw+1 x ubw and pq is ubw x ubw
    ! (Althought the pointers have not been updated to have these
    ! dimensions yet).
    if (ubw > 0) then
       ubws(k:n-1)=ubw
       do i=1,ubw
          pl=>a(k-ubw+i-1:k,k+i:n)
          pq=>q(i:ubw,i:ubw)
          numrots(k+i-1)=ubw
          ! Triangularize L
          do j=ubw,1,-1
             rot=lgivens2(pl(j,j),pl(j+1,j))
             cs(j,k+i-1)=rot%cosine; ss(j,k+i-1)=rot%sine
             j1s(j,k+i-1)=j+k-ubw-1; j2s(j,k+i-1)=j+k-ubw
             call rotation_times_general(trp_rot(rot), pl,j,j+1)
             pl(j,j)=0.0_dp
          end do
          ! reveal column k + i
          if (ubw - i + 1 > 1) then
             do j=ubw-i+1,2,-1
                rot=lgivens(pq(j-1,1),pq(j,1))
                call rotation_times_general(trp_rot(rot), pq, j-1,j)
                call general_times_rotation(pl, rot, j-1,j)
             end do
          end if
          pl(:,1)=pl(:,1)*pq(1,1)
       end do
    end if
    call d_regularize_upper_ub(a,n,ubwmax,lbw,b,ubws,ubw)
  end subroutine d_upper_general_to_upper_ub

  subroutine d_regularize_upper_ub(a, n, ubwmax, lbw, b, ubws, ubw)
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
  end subroutine d_regularize_upper_ub
    

  ! Extend a Gram-Schmidt LQ decomposition
  subroutine d_extend_gs(q1, l, rho, q2, error)
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
    q2=q2-matmul(q1,l)
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
       q2=q2-matmul(q1,x)
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
          q2=q2-matmul(q1,x)
          nrm2=norm2(q2)
       end do
    end do
    q2=q2/nrm2
    if (k >= orthmaxits) then
       error = 1
    end if
  end subroutine d_extend_gs

end module decomp
