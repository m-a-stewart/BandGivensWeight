module general_bv
  use orth
  use transforms
  use types
  implicit none
  integer(kind=int32), private, parameter :: nullmaxits=5

  interface upper_to_bv
     module procedure d_upper_to_bv, c_upper_to_bv
  end interface upper_to_bv

  interface f_upper_to_bv
     module procedure f_d_upper_to_bv, f_c_upper_to_bv
  end interface f_upper_to_bv

type(routine_info), parameter :: info_d_upper_to_bv=routine_info(id_d_upper_to_bv, &
     'd_upper_to_bv', &
     [ character(len=error_message_length) :: '', '', 'ub%n /= bv%n' ] )

type(routine_info), parameter :: info_f_d_upper_to_bv=routine_info(id_f_d_upper_to_bv, &
     'f_d_upper_to_bv', &
     [ character(len=error_message_length) :: 'n<1', 'Insufficient Upper Bandwidth in bv', &
     'Insufficient Lower Bandwidth in bv' ] )

type(routine_info), parameter :: info_c_upper_to_bv=routine_info(id_c_upper_to_bv, &
     'c_upper_to_bv', &
     [ character(len=error_message_length) :: '', '', 'ub%n /= bv%n' ] )

type(routine_info), parameter :: info_f_c_upper_to_bv=routine_info(id_f_c_upper_to_bv, &
     'f_c_upper_to_bv', &
     [ character(len=error_message_length) :: 'n<1', 'Insufficient Upper Bandwidth in bv', &
     'Insufficient Lower Bandwidth in bv' ] )

contains
  
  ! Errors:
  ! 0: no error
  ! 3: n is not the same for a and ub or bv.
  subroutine d_upper_to_bv(a,bv,lbw,tol,error)
    real(kind=dp), target, dimension(:,:), intent(inout) :: a
    type(d_bv), intent(inout) :: bv
    type(error_info), intent(out) :: error
    real(kind=dp), intent(in) :: tol
    integer(kind=int32), intent(in) :: lbw
    if (get_n(bv) /= size(a,1) .or. get_n(bv) /= size(a,2)) then
       call set_error(error, 3, id_d_upper_to_bv)
    end if
    call f_d_upper_to_bv(a,get_n(bv),bv%b, lbw, bv%ubw, get_lbwmax(bv), get_ubwmax(bv), &
         bv%numrotsv, bv%ksv, bv%csv, bv%ssv, tol, error)
    bv%lbw=lbw
  end subroutine d_upper_to_bv

  ! Errors:
  ! 0: no error
  ! 1: n<1
  ! 2: insufficient upper bw in bv%b
  ! 3: insufficient lower bw in bv%b
  subroutine f_d_upper_to_bv(a, n, b, lbw, ubw, lbwmax, ubwmax, &
       numrots, ks, cs, ss, tol, error)
    real(kind=dp), target, dimension(n,n), intent(inout) :: a
    integer(kind=int32), dimension(n,ubwmax), intent(out) :: ks
    real(kind=dp), dimension(n,ubwmax), intent(out) :: cs, ss
    real(kind=dp), dimension(n,lbwmax+ubwmax+1), intent(out) :: b
    real(kind=dp), intent(in) :: tol
    integer(kind=int32), dimension(n), intent(out) :: numrots
    integer(kind=int32), intent(out) :: ubw
    type(error_info) :: error
    integer(kind=int32), intent(in) :: n, lbw, lbwmax, ubwmax
    !
    real(kind=dp), target, dimension(n,ubwmax+1) :: q
    real(kind=dp), dimension(ubwmax+1) :: x
    real(kind=dp) :: nrma, nrmq1
    real(kind=dp), pointer, dimension(:,:) :: pl, pq
    integer(kind=int32) :: i, j, k, roffs, coffs, nl, klast
    integer(kind=int32), dimension(n) :: ubws
    type(d_rotation) :: rot
    !
    q=0.0_dp; numrots=0;
    ss=0.0_dp; cs=0.0_dp; ks=0
    ubw=0; ubws=0
    call clear_error(error)
    nrma = maxabs(a)*sqrt(real(n))
    !
    if (n < 1) then
       call set_error(error, 1, id_f_d_upper_to_bv); return
    end if
    if (n == 1) then
       b(1,1)=a(1,1)
       return
    end if
    if (lbwmax < lbw) then
       call set_error(error, 3, id_f_d_upper_to_bv); return
    end if

    ! Compute an initial trivial QL factorization
    q(1:n-1,1) = a(1:n-1,n)
    a(1:n-1,n)=0.0_dp
    a(n-1,n)=norm2(q(1:n-1,1))
    if (a(n-1,n) == 0.0_dp) then
       q(1:n-1,1) = 1.0_dp
       q(1:n-1,1)=q(1:n-1,1)/norm2(q(1:n-1,1))
    else
       q(1:n-1,1)=q(1:n-1,1)/a(n-1,n)
    end if
    nl=1
    klast=n
    kloop: do k=1,n
       ! Current, possibly singular, L should be contained in
       ! a(n-k-nl+1:n-k,n-k+1:n-k+nl)
       ! At the start of the loop nl is a possible overestimate that
       ! can be decreased by finding a null vector.
       roffs=n-k-nl
       coffs=n-k
       pl => a(roffs+1:roffs+nl, coffs+1:coffs+nl)
       pq => q(1:n-k,1:nl)
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
                ! set up a 1x1 LQ factorization for computing the next V_k
                q(1:n-k-1,1) = a(1:n-k-1,coffs)
                a(1:n-k-1,coffs)=0.0_dp
                a(n-k-1,coffs)=norm2(q(1:n-k-1,1))
                if (a(n-k-1,coffs) == 0.0_dp) then
                   q(1:n-k-1,1)=1.0_dp
                   q(1:n-k-1,1)=q(1:n-k-1,1)/norm2(q(1:n-k-1,1))
                else
                   q(1:n-k-1,1)=q(1:n-k-1,1)/a(n-k-1,coffs)
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
                pl => a(roffs:roffs+1,coffs+1:coffs+1) ! extend pl up
                pl(1,1)=pl(2,1) ! shift
                pl(2,1)=pl(1,1)*pq(n-k,1) ! reveal row n-k
                nrmq1 = norm2(pq(1:n-k-1,1))
                if (nrmq1 == 0.0_dp) then
                   pl(1,1)=0.0_dp;   pq(1:n-k-1,1)=1.0_dp
                   nrmq1 = norm2(pq(1:n-k-1,1))
                   pq(1:n-k-1,1)=pq(1:n-k-1,1)/nrmq1
                else
                   pq(1:n-k-1,1)=pq(1:n-k-1,1)/nrmq1
                   pl(1,1) = pl(1,1)*nrmq1
                end if
                pq=>q(1:n-k-1,1:2)
                if (k < n-2) then
                   call right_shift(pq)
                   pq(:,1)=a(1:n-k-1,coffs)
                   pl => a(roffs-1:roffs,coffs:coffs+1)
                   call extend_gs_columns(pq(:,2:2),pl(2:2,1), pl(1,1), pq(:,1), error)
                   if (error%code > 0) then
                      call add_id(error,id_f_d_upper_to_bv)
                      return
                   end if
                   if (n-k-3 >= 1) then
                      a(1:n-k-3,coffs)=0.0_dp
                   end if
                   nl=nl+1
                else ! k==n-2; only needed to extend L down before exiting.
                   a(1,2)=a(1,2)/pq(1,1)
                   klast=k+1
                   exit kloop                   
                end if
             end if
          end if
       else ! nl > 1
          call lower_right_nullvec(x(1:nl),pl,tol*nrma,nullmaxits,error)
          if (error%code <= 0) then ! if there is a left null vector then introduce a zero row.
             ubws(k)=nl-1
             if (error%code == -nl) then
                pl(nl,nl)=0.0_dp
                numrots(k)=0
             call clear_error(error)
             else if (error%code <= -1) then
                numrots(k)=nl+error%code
                pl(-error%code,-error%code)=0.0_dp
                do j=-error%code,nl-1
                   rot=rgivens(pl(j+1,j),pl(j+1,j+1))
                   call general_times_rotation(pl(j+1:nl,:),rot,j,j+1)
                   pl(j+1,j+1)=0.0_dp
                   cs(k,nl-j)=rot%cosine; ss(k,nl-j)=rot%sine
                   ks(k,nl-j)=coffs+j
                end do
             call clear_error(error)
             else ! error=0
                numrots(k)=nl-1;
                do j=2,nl ! apply v_k while preserving the triangular structure of L
                   rot=lgivens2(x(j-1),x(j))
                   call rotation_times_general(trp_rot(rot),x, j-1,j)
                   call general_times_rotation(pl,rot,j-1,j)
                   cs(k,nl-j+1)=rot%cosine; ss(k,nl-j+1)=rot%sine
                   ks(k,nl-j+1)=coffs+j-1
                   rot=lgivens2(pl(j-1,j),pl(j,j))
                   call rotation_times_general(trp_rot(rot), pl(:,1:j), j-1,j)
                   call general_times_rotation(pq,rot,j-1,j)
                   pl(j-1,j)=0.0_dp
                end do
                pl(nl,nl)=0.0_dp
             end if
             do j=nl-1,1,-1 ! compress (TODO: this is not necessary if zero is in pl(1,1))
                rot=lgivens(pl(j,j),pl(nl,j))
                call rotation_times_general(trp_rot(rot),pl(:,1:j),j,nl)
                call general_times_rotation(pq,rot,j,nl)
                pl(nl,j)=0.0_dp
             end do
             if (k+nl==n) then ! q is square and nl by nl
                ! reveal row nl
                do j=1,nl-1
                   rot=rgivens2(pq(nl,j),pq(nl,nl))
                   call general_times_rotation(pq,rot,j,nl)
                   call rotation_times_general(trp_rot(rot),pl(:,1:j),j,nl)
                end do
                pl(nl,:)=pq(nl,nl)*pl(nl,:)
                ! extend one row
                x(1:nl-1)=0.0_dp
                do j=1,nl-1
                   do i=1,nl-1
                      x(j)=x(j)+pq(i,j)*a(i,coffs)
                   end do
                end do
                a(1:nl-1,coffs)=x(1:nl-1)
                klast=k+1
                exit kloop ! terminate
             else
                pq(:,nl)=0.0_dp; pq(n-k,nl)=1.0_dp
                call extend_gs_columns(pq(:,1:nl-1), x(1:nl-1), x(nl), pq(:,nl), error)
                if (error%code > 0) then
                   call add_id(error,id_f_d_upper_to_bv)
                   return
                end if
                do j=1,nl-1
                   rot=rgivens2(pq(n-k,j),pq(n-k,nl))
                   call general_times_rotation(pq,rot,j,nl)
                   call rotation_times_general(trp_rot(rot),pl(:,1:j),j,nl)
                end do
                pl(nl,:)=pq(n-k,nl)*pl(nl,:)
                ! extend the LQ factorization with a(1:n-k-1,coffs)
                call right_shift(pq)
                pq => q(1:n-k-1,1:nl)
                pq(:,1)=a(1:n-k-1,coffs)
                pl => a(roffs:roffs+nl-1, coffs:coffs+nl-1)
                call extend_gs_columns(pq(:,2:nl),pl(2:nl,1),pl(1,1),pq(:,1),error)
                if (error%code > 0) then
                   call add_id(error,id_f_d_upper_to_bv)
                   return
                end if
                if (roffs>1) then
                   a(1:roffs-1,coffs)=0.0_dp
                end if
             end if
          else
             ! no null vector found.  Simply reveal row nl if there is room.
             ! Otherwise terminate with square L
             ubws(k)=nl
             if (k+nl == n) then
                klast=k ! terminate with square L
                exit kloop
             else
                pl => a(roffs:roffs+nl, coffs+1:coffs+nl) ! extend pl up
                call up_shift(pl)
                pq => q(1:n-k,1:nl+1)
                pq(:,nl+1)=0.0_dp; pq(n-k,nl+1)=1.0_dp
                call extend_gs_columns(pq(:,1:nl),x(1:nl), x(nl+1),pq(:,nl+1),error)
                if (error%code > 0) then
                   call add_id(error,id_f_d_upper_to_bv)
                   return
                end if
                do j=1,nl
                   rot=rgivens2(pq(n-k,j),pq(n-k,nl+1))
                   call general_times_rotation(pq,rot,j,nl+1)
                   call rotation_times_general(trp_rot(rot),pl(:,1:j),j,nl+1)
                end do
                pl(nl+1,:)=pl(nl+1,:)*pq(n-k,nl+1)
                if (nl==n-k-1) then ! q is now nl x nl
                   pq => q(1:nl,1:nl)
                   x(1:nl)=0.0_dp
                   do j=1,nl
                      do i=1,nl
                         x(j)=x(j)+pq(i,j)*a(i,coffs)
                      end do
                   end do
                   a(1:nl,coffs)=x(1:nl)
                   klast=k+1
                   exit kloop
                else ! q is not square.  Make L (nl+1)x(nl+1)
                   pq => q(1:n-k-1,1:nl+1)
                   call right_shift(pq)
                   pq(:,1)=a(1:n-k-1,coffs)
                   pl => a(roffs-1:roffs-1+nl,coffs:coffs+nl)
                   call extend_gs_columns(pq(:,2:nl+1), pl(2:nl+1,1), pl(1,1), pq(:,1), error)
                   if (error%code>0) then
                      call add_id(error,id_f_d_upper_to_bv)
                      return
                   end if
                   if (roffs-2 >= 1) then
                      a(1:roffs-2,coffs)=0.0_dp
                   end if
                   nl=nl+1
                end if
             end if
          end if ! null vector check
       end if
    end do kloop
    ! If klast = k+1, we need to terminate on an nl by nl+1 matrix L
    ! contained in a(1:nl,n-klast+1:n-klast+nl+1)
    ! If klast=k, then L is square and
    ! contained in a(1:n-k,n-k+1:n-k+(n-k)).  The former happens when the last step
    ! of kloop found a null vector.  The latter happens when it didn't.
    if (klast == k) then ! square termination
       ubws(k:n-1)=n-k
       nl = n-k
       coffs=n-k
       if (n-k > 1) then
          do i=1,n-k-1 ! 1,nl-1
             nl = n-k-i+1
             pl => a(1:nl,coffs+1:coffs+nl)
             pq => q(1:nl,1:nl)
             ! reveal row nl
             do j=1,nl-1
                rot=rgivens2(pq(nl,j),pq(nl,j+1))
                call general_times_rotation(pq,rot,j,j+1)
                call rotation_times_general(trp_rot(rot), pl(:,1:j+1), j, j+1)
             end do
             pl(nl,:)=pq(nl,nl)*pl(nl,:)
             numrots(k+i)=nl-1
             do j=1,nl-1
                rot=rgivens(pl(j,j),pl(j,j+1))
                call general_times_rotation(pl(j:nl-1,:),rot,j,j+1)
                pl(j,j+1)=0.0_dp
                cs(k+i,nl-j)=rot%cosine; ss(k+i,nl-j)=rot%sine
                ks(k+i,nl-j)=coffs+j
             end do
          end do
          pl(1,1)=pq(1,1)*pl(1,1)
       end if
    else ! rectangular termination
       k=k+1
       ubws(k:n-1)=n-k
       nl=n-k
       coffs=n-k
       do i=1,n-k
          nl=n-k-i+1
          pl=>a(1:nl,coffs+1:coffs+nl+1)
          pq=>q(1:nl,1:nl)
          numrots(k+i-1)=nl
          do j=1,nl
             rot=rgivens(pl(j,j),pl(j,j+1))
             cs(k+i-1,nl-j+1)=rot%cosine; ss(k+i-1,nl-j+1)=rot%sine
             ks(k+i-1,nl-j+1)=coffs+j
             call general_times_rotation(pl(j:nl,:),rot,j,j+1)
             pl(j,j+1)=0.0_dp
          end do
          ! reveal row nl
          do j=1,nl-1
             rot=rgivens2(pq(nl,j),pq(nl,j+1))
             call general_times_rotation(pq,rot,j,j+1)
             call rotation_times_general(trp_rot(rot),pl(:,1:j+1),j,j+1)
          end do
          pl(nl,:)=pq(nl,nl)*pl(nl,:)
       end do
    end if
    if (maxval(ubws) > ubwmax) then
       call set_error(error, 2, id_f_d_upper_to_bv)
    else
       call d_extract_diagonals_bv(a, n, b, lbw, ubw, lbwmax, ubwmax, ubws)
    end if

  end subroutine f_d_upper_to_bv

  subroutine d_extract_diagonals_bv(a, n, b, lbw, ubw, lbwmax, ubwmax, ubws)
    real(kind=dp), target, dimension(n,n), intent(in) :: a
    real(kind=dp), dimension(n,lbwmax+ubwmax+1), intent(out) :: b
    integer(kind=int32), intent(out) :: ubw
    integer(kind=int32), intent(in) :: n, lbw, lbwmax, ubwmax
    integer(kind=int32), dimension(n), intent(in) :: ubws
    !
    integer(kind=int32) :: j, d
    ubw=maxval(ubws)
    ! put diagonals in b
    b=0.0_dp
    do d=1,lbw+1
       do j=lbw-d+2, n
          b(j,d)=a(j,d+j-lbw-1)
       end do
    end do
    do d=lbw+2,ubw+lbw+1
       do j=1, n-d+lbw+1
          b(j,d)=a(j,d+j-lbw-1)
       end do
    end do
  end subroutine d_extract_diagonals_bv

! complex BV

  subroutine c_upper_to_bv(a,bv,lbw,tol,error)
    complex(kind=dp), target, dimension(:,:), intent(inout) :: a
    type(c_bv), intent(inout) :: bv
    type(error_info), intent(out) :: error
    real(kind=dp), intent(in) :: tol
    integer(kind=int32), intent(in) :: lbw
    if (get_n(bv) /= size(a,1) .or. get_n(bv) /= size(a,2)) then
       call set_error(error, 3, id_c_upper_to_bv); return
    end if
    call f_c_upper_to_bv(a,get_n(bv),bv%b, lbw, bv%ubw, get_lbwmax(bv), get_ubwmax(bv), &
         bv%numrotsv, bv%ksv, bv%csv, bv%ssv, tol, error)
    bv%lbw=lbw
  end subroutine c_upper_to_bv


  subroutine f_c_upper_to_bv(a, n, b, lbw, ubw, lbwmax, ubwmax, &
       numrots, ks, cs, ss, tol, error)
    complex(kind=dp), target, dimension(n,n), intent(inout) :: a
    integer(kind=int32), dimension(n,ubwmax), intent(out) :: ks
    complex(kind=dp), dimension(n,ubwmax), intent(out) :: cs, ss
    complex(kind=dp), dimension(n,lbwmax+ubwmax+1), intent(out) :: b
    real(kind=dp), intent(in) :: tol
    integer(kind=int32), dimension(n), intent(out) :: numrots
    integer(kind=int32), intent(out) :: ubw
    type(error_info), intent(out) :: error
    integer(kind=int32), intent(in) :: n, lbw, lbwmax, ubwmax
    !
    complex(kind=dp), target, dimension(n,ubwmax+1) :: q
    complex(kind=dp), dimension(ubwmax+1) :: x
    real(kind=dp) :: nrma, nrmq1
    complex(kind=dp), pointer, dimension(:,:) :: pl, pq
    integer(kind=int32) :: i, j, k, roffs, coffs, nl, klast
    integer(kind=int32), dimension(n) :: ubws
    type(c_rotation) :: rot
    !
    q=(0.0_dp, 0.0_dp); numrots=0;
    ss=(0.0_dp, 0.0_dp); cs=(0.0_dp,0.0_dp); ks=0
    ubw=0; ubws=0
    call clear_error(error)
    nrma = maxabs(a)*sqrt(real(n))
    !
    if (n < 1) then
       call set_error(error, 1, id_f_c_upper_to_bv); return
    end if
    if (n == 1) then
       b(1,1)=a(1,1)
       return
    end if
    if (lbwmax < lbw) then
       call set_error(error, 3, id_f_c_upper_to_bv); return
    end if
    ! Compute an initial trivial QL factorization
    q(1:n-1,1) = a(1:n-1,n)
    a(1:n-1,n)=(0.0_dp, 0.0_dp)
    a(n-1,n)=norm2(q(1:n-1,1))
    if (a(n-1,n) == (0.0_dp, 0.0_dp)) then
       q(1:n-1,1) = (1.0_dp, 0.0_dp)
       q(1:n-1,1)=q(1:n-1,1)/norm2(q(1:n-1,1))
    else
       q(1:n-1,1)=q(1:n-1,1)/a(n-1,n)
    end if
    nl=1
    klast=n
    kloop: do k=1,n
       ! Current, possibly singular, L should be contained in
       ! a(n-k-nl+1:n-k,n-k+1:n-k+nl)
       ! At the start of the loop nl is a possible overestimate that
       ! can be decreased by finding a null vector.
       roffs=n-k-nl
       coffs=n-k
       pl => a(roffs+1:roffs+nl, coffs+1:coffs+nl)
       pq => q(1:n-k,1:nl)
       if (nl==1) then
          if (abs(pl(1,1)) < tol*nrma) then ! null vector
             ubws(k)=0;  pl(1,1)=(0.0_dp, 0.0_dp)
             if (k==n-1) then ! k+nl=n
                klast=k+1
                exit kloop
             else if (k==n-2) then ! k+nl=n-1
                q(1,1)=(1.0_dp, 0.0_dp)
                klast=k+1
                exit kloop
             else
                ! set up a 1x1 LQ factorization for computing the next V_k
                q(1:n-k-1,1) = a(1:n-k-1,coffs)
                a(1:n-k-1,coffs)=(0.0_dp, 0.0_dp)
                a(n-k-1,coffs)=norm2(q(1:n-k-1,1))
                if (a(n-k-1,coffs) == (0.0_dp, 0.0_dp)) then
                   q(1:n-k-1,1)=(1.0_dp, 0.0_dp)
                   q(1:n-k-1,1)=q(1:n-k-1,1)/norm2(q(1:n-k-1,1))
                else
                   q(1:n-k-1,1)=q(1:n-k-1,1)/a(n-k-1,coffs)
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
                pl => a(roffs:roffs+1,coffs+1:coffs+1) ! extend pl up
                pl(1,1)=pl(2,1) ! shift
                pl(2,1)=pl(1,1)*pq(n-k,1) ! reveal row n-k
                nrmq1 = norm2(pq(1:n-k-1,1))
                if (nrmq1 == 0.0_dp) then
                   pl(1,1)=(0.0_dp, 0.0_dp);   pq(1:n-k-1,1)=(1.0_dp, 0.0_dp)
                   nrmq1 = norm2(pq(1:n-k-1,1))
                   pq(1:n-k-1,1)=pq(1:n-k-1,1)/nrmq1
                else
                   pq(1:n-k-1,1)=pq(1:n-k-1,1)/nrmq1
                   pl(1,1) = pl(1,1)*nrmq1
                end if
                pq=>q(1:n-k-1,1:2)
                if (k < n-2) then
                   call right_shift(pq)
                   pq(:,1)=a(1:n-k-1,coffs)
                   pl => a(roffs-1:roffs,coffs:coffs+1)
                   call extend_gs_columns(pq(:,2:2),pl(2:2,1), pl(1,1), pq(:,1), error)
                   if (error%code > 0 ) then
                      call add_id(error,id_f_c_upper_to_bv)
                      return
                   end if
                   if (n-k-3 >= 1) then
                      a(1:n-k-3,coffs)=(0.0_dp, 0.0_dp)
                   end if
                   nl=nl+1
                else ! k==n-2; only needed to extend L down before exiting.
                   a(1,2)=a(1,2)/pq(1,1)
                   klast=k+1
                   exit kloop                   
                end if
             end if
          end if
       else ! nl > 1
          call lower_right_nullvec(x(1:nl),pl,tol*nrma,nullmaxits,error)
          if (error%code <= 0) then ! if there is a left null vector then introduce a zero row.
             ubws(k)=nl-1
             if (error%code == -nl) then
                pl(nl,nl)=(0.0_dp, 0.0_dp)
                numrots(k)=0
                call clear_error(error)
             else if (error%code <= -1) then
                numrots(k)=nl+error%code
                pl(-error%code,-error%code)=(0.0_dp, 0.0_dp)
                do j=-error%code,nl-1
                   rot=rgivens(pl(j+1,j),pl(j+1,j+1))
                   call general_times_rotation(pl(j+1:nl,:),rot,j,j+1)
                   pl(j+1,j+1)=(0.0_dp, 0.0_dp)
                   cs(k,nl-j)=rot%cosine; ss(k,nl-j)=rot%sine
                   ks(k,nl-j)=coffs+j
                end do
                call clear_error(error)
             else ! nullerr=0
                numrots(k)=nl-1;
                do j=2,nl ! apply v_k while preserving the triangular structure of L
                   rot=lgivens2(x(j-1),x(j))
                   call rotation_times_general(trp_rot(rot),x, j-1,j)
                   call general_times_rotation(pl,rot,j-1,j)
                   cs(k,nl-j+1)=rot%cosine; ss(k,nl-j+1)=rot%sine
                   ks(k,nl-j+1)=coffs+j-1
                   rot=lgivens2(pl(j-1,j),pl(j,j))
                   call rotation_times_general(trp_rot(rot), pl(:,1:j), j-1,j)
                   call general_times_rotation(pq,rot,j-1,j)
                   pl(j-1,j)=(0.0_dp, 0.0_dp)
                end do
                pl(nl,nl)=(0.0_dp, 0.0_dp)
             end if
             do j=nl-1,1,-1 ! compress (TODO: this is no necessary if zero is in pl(1,1))
                rot=lgivens(pl(j,j),pl(nl,j))
                call rotation_times_general(trp_rot(rot),pl(:,1:j),j,nl)
                call general_times_rotation(pq,rot,j,nl)
                pl(nl,j)=(0.0_dp, 0.0_dp)
             end do
             if (k+nl==n) then ! q is square and nl by nl
                ! reveal row nl
                do j=1,nl-1
                   rot=rgivens2(pq(nl,j),pq(nl,nl))
                   call general_times_rotation(pq,rot,j,nl)
                   call rotation_times_general(trp_rot(rot),pl(:,1:j),j,nl)
                end do
                pl(nl,:)=pq(nl,nl)*pl(nl,:)
                ! extend one row
                x(1:nl-1)=(0.0_dp, 0.0_dp)
                do j=1,nl-1
                   do i=1,nl-1
                      x(j)=x(j)+conjg(pq(i,j))*a(i,coffs)
                   end do
                end do
                a(1:nl-1,coffs)=x(1:nl-1)
                klast=k+1
                exit kloop ! terminate
             else
                pq(:,nl)=(0.0_dp, 0.0_dp); pq(n-k,nl)=(1.0_dp, 0.0_dp)
                call extend_gs_columns(pq(:,1:nl-1), x(1:nl-1), x(nl), pq(:,nl), error)
                if (error%code > 0) then
                   call add_id(error,id_f_c_upper_to_bv)
                   return
                end if
                do j=1,nl-1
                   rot=rgivens2(pq(n-k,j),pq(n-k,nl))
                   call general_times_rotation(pq,rot,j,nl)
                   call rotation_times_general(trp_rot(rot),pl(:,1:j),j,nl)
                end do
                pl(nl,:)=pq(n-k,nl)*pl(nl,:)
                ! extend the LQ factorization with a(1:n-k-1,coffs)
                call right_shift(pq)
                pq => q(1:n-k-1,1:nl)
                pq(:,1)=a(1:n-k-1,coffs)
                pl => a(roffs:roffs+nl-1, coffs:coffs+nl-1)
                call extend_gs_columns(pq(:,2:nl),pl(2:nl,1),pl(1,1),pq(:,1),error)
                if (error%code > 0) then
                   call add_id(error,id_f_c_upper_to_bv)
                   return
                end if
                if (roffs>1) then
                   a(1:roffs-1,coffs)=(0.0_dp, 0.0_dp)
                end if
             end if
          else
             ! no null vector found.  Simply reveal row nl if there is room.
             ! Otherwise terminate with square L
             ubws(k)=nl
             if (k+nl == n) then
                klast=k ! terminate with square L
                exit kloop
             else
                pl => a(roffs:roffs+nl, coffs+1:coffs+nl) ! extend pl up
                call up_shift(pl)
                pq => q(1:n-k,1:nl+1)
                pq(:,nl+1)=(0.0_dp, 0.0_dp); pq(n-k,nl+1)=(1.0_dp, 0.0_dp)
                call extend_gs_columns(pq(:,1:nl),x(1:nl), x(nl+1),pq(:,nl+1),error)
                if (error%code > 0) then
                   call add_id(error,id_f_c_upper_to_bv)
                   return
                end if
                do j=1,nl
                   rot=rgivens2(pq(n-k,j),pq(n-k,nl+1))
                   call general_times_rotation(pq,rot,j,nl+1)
                   call rotation_times_general(trp_rot(rot),pl(:,1:j),j,nl+1)
                end do
                pl(nl+1,:)=pl(nl+1,:)*pq(n-k,nl+1)
                if (nl==n-k-1) then ! q is now nl x nl
                   pq => q(1:nl,1:nl)
                   x(1:nl)=(0.0_dp, 0.0_dp)
                   do j=1,nl
                      do i=1,nl
                         x(j)=x(j)+conjg(pq(i,j))*a(i,coffs)
                      end do
                   end do
                   a(1:nl,coffs)=x(1:nl)
                   klast=k+1
                   exit kloop
                else ! q is not square.  Make L (nl+1)x(nl+1)
                   pq => q(1:n-k-1,1:nl+1)
                   call right_shift(pq)
                   pq(:,1)=a(1:n-k-1,coffs)
                   pl => a(roffs-1:roffs-1+nl,coffs:coffs+nl)
                   call extend_gs_columns(pq(:,2:nl+1), pl(2:nl+1,1), pl(1,1), pq(:,1), error)
                   if (error%code > 0) then
                      call add_id(error,id_f_c_upper_to_bv)
                      return
                   end if
                   if (roffs-2 >= 1) then
                      a(1:roffs-2,coffs)=(0.0_dp, 0.0_dp)
                   end if
                   nl=nl+1
                end if
             end if
          end if ! null vector check
       end if
    end do kloop
    ! If klast = k+1, we need to terminate on an nl by nl+1 matrix L
    ! contained in a(1:nl,n-klast+1:n-klast+nl+1)
    ! If klast=k, then L is square and
    ! contained in a(1:n-k,n-k+1:n-k+(n-k)).  The former happens when the last step
    ! of kloop found a null vector.  The latter happens when it didn't.
    if (klast == k) then ! square termination
       ubws(k:n-1)=n-k
       nl = n-k
       coffs=n-k
       if (n-k > 1) then
          do i=1,n-k-1 ! 1,nl-1
             nl = n-k-i+1
             pl => a(1:nl,coffs+1:coffs+nl)
             pq => q(1:nl,1:nl)
             ! reveal row nl
             do j=1,nl-1
                rot=rgivens2(pq(nl,j),pq(nl,j+1))
                call general_times_rotation(pq,rot,j,j+1)
                call rotation_times_general(trp_rot(rot), pl(:,1:j+1), j, j+1)
             end do
             pl(nl,:)=pq(nl,nl)*pl(nl,:)
             numrots(k+i)=nl-1
             do j=1,nl-1
                rot=rgivens(pl(j,j),pl(j,j+1))
                call general_times_rotation(pl(j:nl-1,:),rot,j,j+1)
                pl(j,j+1)=(0.0_dp, 0.0_dp)
                cs(k+i,nl-j)=rot%cosine; ss(k+i,nl-j)=rot%sine
                ks(k+i,nl-j)=coffs+j
             end do
          end do
          pl(1,1)=pq(1,1)*pl(1,1)
       end if
    else ! rectangular termination
       k=k+1
       ubws(k:n-1)=n-k
       nl=n-k
       coffs=n-k
       do i=1,n-k
          nl=n-k-i+1
          pl=>a(1:nl,coffs+1:coffs+nl+1)
          pq=>q(1:nl,1:nl)
          numrots(k+i-1)=nl
          do j=1,nl
             rot=rgivens(pl(j,j),pl(j,j+1))
             cs(k+i-1,nl-j+1)=rot%cosine; ss(k+i-1,nl-j+1)=rot%sine
             ks(k+i-1,nl-j+1)=coffs+j
             call general_times_rotation(pl(j:nl,:),rot,j,j+1)
             pl(j,j+1)=(0.0_dp, 0.0_dp)
          end do
          ! reveal row nl
          do j=1,nl-1
             rot=rgivens2(pq(nl,j),pq(nl,j+1))
             call general_times_rotation(pq,rot,j,j+1)
             call rotation_times_general(trp_rot(rot),pl(:,1:j+1),j,j+1)
          end do
          pl(nl,:)=pq(nl,nl)*pl(nl,:)
       end do
    end if
    if (maxval(ubws) > ubwmax) then
       call set_error(error, 2, id_f_c_upper_to_bv)
    else
       call c_extract_diagonals_bv(a,n,b,lbw,ubw,lbwmax, ubwmax,ubws)
    end if
  end subroutine f_c_upper_to_bv

  subroutine c_extract_diagonals_bv(a, n, b, lbw, ubw, lbwmax, ubwmax, ubws)
    complex(kind=dp), target, dimension(n,n), intent(in) :: a
    complex(kind=dp), dimension(n,lbwmax+ubwmax+1), intent(out) :: b
    integer(kind=int32), intent(out) :: ubw
    integer(kind=int32), intent(in) :: n, lbw, lbwmax, ubwmax
    integer(kind=int32), dimension(n), intent(in) :: ubws
    !
    integer(kind=int32) :: j, d
    ubw=maxval(ubws)
    ! put diagonals in b
    b=(0.0_dp, 0.0_dp)
    do d=1,lbw+1
       do j=lbw-d+2, n
          b(j,d)=a(j,d+j-lbw-1)
       end do
    end do
    do d=lbw+2,ubw+lbw+1
       do j=1, n-d+lbw+1
          b(j,d)=a(j,d+j-lbw-1)
       end do
    end do
  end subroutine c_extract_diagonals_bv

end module general_bv
