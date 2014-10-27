module general_ub
  use orth
  use shift
  use rotation
  use types
  implicit none
  integer(kind=int32), private, parameter :: nullmaxits=5

  interface upper_to_ub
     module procedure d_upper_to_ub, c_upper_to_ub
  end interface upper_to_ub

  interface f_upper_to_ub
     module procedure f_d_upper_to_ub, f_c_upper_to_ub
  end interface f_upper_to_ub

  interface f_general_ub
     module procedure f_d_general_ub, f_c_general_ub
  end interface f_general_ub

  type(routine_info), parameter :: info_d_upper_to_ub=routine_info(id_d_upper_to_ub, &
       'd_upper_to_ub', &
       [ character(len=error_message_length) :: 'n<1', 'ub%lbwmax < lbw', 'Size of a and ub not the same.' ] )

  type(routine_info), parameter :: info_f_d_upper_to_ub=routine_info(id_f_d_upper_to_ub, &
       'f_d_upper_to_ub', &
       [ character(len=error_message_length) :: 'Insufficient Upper Bandwidth in ub' ])

  type(routine_info), parameter :: info_f_d_general_ub=routine_info(id_f_d_general_ub, &
       'f_d_general_ub', &
       [ character(len=error_message_length) :: 'Insufficient Upper Bandwidth.' ])

  type(routine_info), parameter :: info_c_upper_to_ub=routine_info(id_c_upper_to_ub, &
       'c_upper_to_ub', &
       [ character(len=error_message_length) :: 'n<1', 'ub%lbwmax < lbw', 'Size of a and ub not the same.' ] )

  type(routine_info), parameter :: info_f_c_upper_to_ub=routine_info(id_f_c_upper_to_ub, &
       'f_c_upper_to_ub', &
       [ character(len=error_message_length) :: 'Insufficient Upper Bandwidth in ub' ])

  type(routine_info), parameter :: info_f_c_general_ub=routine_info(id_f_c_general_ub, &
       'f_c_general_ub', &
       [ character(len=error_message_length) :: 'Insufficient Upper Bandwidth.' ])

contains

  ! Errors:
  ! 0: no error
  ! 1: n < 1
  ! 2: ub%lbwmax < lbw
  ! 3: n is not the same for a and ub or bv.
  subroutine d_upper_to_ub(a,ub,lbw,tol,error)
    real(kind=dp), target, dimension(:,:), intent(inout) :: a
    type(d_ub), intent(inout) :: ub
    type(error_info), intent(out) :: error
    real(kind=dp), intent(in) :: tol
    integer(kind=int32), intent(in) :: lbw
    call clear_error(error)
    if (size(a,1) < 1) then
       call set_error(error, 1, id_d_upper_to_ub); return
    end if
    if (get_lbwmax(ub) < lbw) then
       call set_error(error, 2, id_d_upper_to_ub); return
    end if
    if (get_n(ub) /= size(a,1) .or. get_n(ub) /= size(a,2)) then
       call set_error(error, 3, id_d_upper_to_ub); return
    end if
    call f_d_upper_to_ub(a,get_n(ub),ub%bc, lbw, ub%ubw, get_lbwmax(ub), get_ubwmax(ub), &
         ub%numrotsu, ub%jsu, ub%csu, ub%ssu, tol, error)
    ub%lbw=lbw
  end subroutine d_upper_to_ub


  ! Errors:
  ! 0: no error
  ! 1: insufficient upper bw in ub%bc
  subroutine f_d_upper_to_ub(a, n, b, lbw, ubw, lbwmax, ubwmax, &
       numrotsu, jsu, csu, ssu, tol, error)
    real(kind=dp), target, dimension(n,n), intent(inout) :: a
    integer(kind=int32), dimension(ubwmax,n), intent(out) :: jsu
    real(kind=dp), dimension(ubwmax,n), intent(out) :: csu, ssu
    real(kind=dp), dimension(ubwmax+lbwmax+1,n), intent(out) :: b
    real(kind=dp), intent(in) :: tol
    integer(kind=int32), dimension(n), intent(out) :: numrotsu
    integer(kind=int32), intent(out) :: ubw
    type(error_info), intent(out) :: error
    integer(kind=int32), intent(in) :: n, lbw, ubwmax, lbwmax
    !
    integer(kind=int32), dimension(n) :: ubws

    call clear_error(error)

    if (n == 1) then
       numrotsu=0;
       ssu=0.0_dp; csu=0.0_dp; jsu=0
       ubw=0
       b(1,1)=a(1,1)
       return
    end if

    call f_d_general_ub(a, n, ubws, ubwmax, numrotsu, jsu, csu, ssu, tol, error)
    
    ubw=maxval(ubws)
    if (ubw > ubwmax) then
       ! This should already have been detected in f_d_general_ub.
       call set_error(error, 1, id_f_d_upper_to_ub)
    else
       call d_extract_diagonals_bc(a, n, b, lbw, ubw, lbwmax, ubwmax)
    end if
  end subroutine f_d_upper_to_ub

  subroutine f_d_general_ub(a, n, ubws, ubwmax, numrotsu, jsu, csu, ssu, tol, error)
    real(kind=dp), target, dimension(n,n), intent(inout) :: a
    integer(kind=int32), dimension(ubwmax,n), intent(out) :: jsu
    real(kind=dp), dimension(ubwmax,n), intent(out) :: csu, ssu
    real(kind=dp), intent(in) :: tol
    integer(kind=int32), dimension(n), intent(out) :: numrotsu
    integer(kind=int32), dimension(n), intent(out) :: ubws
    type(error_info), intent(out) :: error
    integer(kind=int32), intent(in) :: n, ubwmax
    !
    real(kind=dp), target, dimension(ubwmax+1,n) :: q
    real(kind=dp), dimension(ubwmax+1) :: x
    real(kind=dp) :: nrma, nrmq12
    real(kind=dp), pointer, dimension(:,:) :: pl, pq
    integer(kind=int32) :: i, j, k, roffs, nl, klast
    type(d_rotation) :: rot
    !
    call clear_error(error)
    q=0.0_dp; numrotsu=0;
    ssu=0.0_dp; csu=0.0_dp; jsu=0
    ubws=0
    nrma = maxabs(a)*sqrt(real(n))
    !
    if (n == 1) then
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
             if (ubwmax==0) then
                call set_error(error, 1, id_f_d_general_ub); return
             end if
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
                call shift2(pq,0,-1)
                if (k < n-2) then
                   pq => q(1:2,1:n-k-1)
                   pq(2,:)=a(k+1,k+2:n)
                   pl => a(k:k+1,k+2:k+3)
                   call extend_gs_rows(pq(1:1,:), pl(2,1:1), pl(2,2), pq(2,:), error)
                   if (error%code > 0) then
                      call add_id(error, id_f_d_general_ub); return
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
          call lower_left_nullvec(x(1:nl),pl,tol*nrma,nullmaxits,error)
          if (error%code <= 0) then ! if there is a left null vector then introduce a zero row.
             ubws(k)=nl-1
             if (error%code == -1) then
                pl(1,1)=0.0_dp
                numrotsu(k)=0
                call clear_error(error)
             else if (error%code < -1) then
                numrotsu(k)=-error%code-1
                pl(-error%code,-error%code)=0.0_dp
                do j=-error%code-1,1,-1
                   rot=lgivens2(pl(j,j),pl(j+1,j))
                   call rotation_times_general(trp_rot(rot), pl(:,1:j), j,j+1)
                   pl(j,j)=0.0_dp
                   csu(j,k)=rot%cosine; ssu(j,k)=rot%sine
                   jsu(j,k)=roffs+j
                end do
                call clear_error(error)
             else ! error=0
                numrotsu(k)=nl-1;
                do j=nl,2,-1 ! apply u_k while preserving the triangular structure of L
                   rot=rgivens(x(j-1),x(j))
                   call general_times_rotation(x,rot,j-1,j)
                   call rotation_times_general(trp_rot(rot), pl,j-1,j)
                   csu(j-1,k)=rot%cosine; ssu(j-1,k)=rot%sine
                   jsu(j-1,k)=roffs+j-1
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
                call shift2(pq,-1,-1)
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
                call extend_gs_rows(pq(2:nl,:), x(1:nl-1), x(nl), pq(1,:), error) ! orthogonalize
                if (error%code > 0) then
                   call add_id(error, id_f_d_general_ub); return
                end if
                do j=nl,2,-1
                   rot=lgivens(pq(1,1),pq(j,1))
                   call rotation_times_general(trp_rot(rot), pq, 1,j)
                   call general_times_rotation(pl(j:nl,:),rot,1,j)
                end do
                pl(:,1)=pl(:,1)*pq(1,1)
                call shift2(pq,-1,-1)
                ! extend the LQ factorization
                pq => q(1:nl,1:n-k-1)
                pq(nl,:)=a(k+1,k+2:n)
                pl => a(roffs+2:k+1,k+2:k+nl+1) ! nl by nl
                call extend_gs_rows(pq(1:nl-1,:), pl(nl,1:nl-1), pl(nl,nl), pq(nl,:), error)
                if (error%code > 0) then
                   call add_id(error, id_f_d_general_ub); return
                end if
                if (k+nl+2 <= n) then
                   a(k+1,k+nl+2:n)=0.0_dp
                end if
             end if
          else
             call clear_error(error)
             ! no null vector found.  Simply reveal column k+1 if there is room.  Otherwise terminate
             !  with square L.
             if (ubwmax<nl) then
                call set_error(error, 1, id_f_d_general_ub); return
             end if
             ubws(k)=nl
             if (k+nl==n) then
                klast=k ! terminate with square L.
                exit kloop
             else
                pl => a(roffs+1:k, k+1:k+nl+1) ! extend pl to the right. (note this requires k+nl < n)
                call shift2(pl,0,1)
                pq => q(1:nl+1, 1:n-k)
                call shift2(pq,1,0)
                pq(1,:)=0.0_dp
                pq(1,1)=1.0_dp
                ! orthogonalizing.  Note x is used as workspace for coefficients.
                call extend_gs_rows(pq(2:nl+1,:), x(1:nl), x(nl+1), pq(1,:), error)
                if (error%code > 0) then
                   call add_id(error, id_f_d_general_ub); return
                end if
                do j=nl+1,2,-1
                   rot=lgivens(pq(1,1),pq(j,1))
                   call rotation_times_general(trp_rot(rot), pq, 1,j)
                   call general_times_rotation(pl(j-1:nl,:),rot,1,j)
                end do
                pl(:,1)=pl(:,1)*pq(1,1)
                call shift2(pq,-1,-1)
                if (nl==n-k-1) then ! q is now nl by nl.  Extend the LQ factorization down one row before stopping.
                   pq => q(1:nl,1:nl)
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
                   call extend_gs_rows(pq(1:nl,:), pl(nl+1,1:nl), pl(nl+1,nl+1), pq(nl+1,:), error)
                   if (error%code > 0) then
                      call add_id(error, id_f_d_general_ub); return
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
       if (ubwmax<n-k) then
          call set_error(error, 1, id_f_d_general_ub); return
       end if
       ubws(k:n-1)=n-k
       nl=n-k
       if (nl > 1) then
          do i=1,nl-1 ! nl-1
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
             numrotsu(k+i)=nl-1
             do j=nl-1,1,-1
                rot=lgivens2(pl(j,j+1),pl(j+1,j+1))
                call rotation_times_general(trp_rot(rot),pl(:,2:nl),j,j+1)
                csu(j,k+i)=rot%cosine; ssu(j,k+i)=rot%sine
                jsu(j,k+i)=roffs+j
                pl(j,j+1)=0.0_dp
             end do
          end do
          pl(2,2)=pl(2,2)*pq(2,2)
       end if
    else ! rectangular termination
       k=k+1
       if (ubwmax < n-k) then
          call set_error(error, 1, id_f_d_general_ub); return
       end if
       ubws(k:n-1)=n-k
       nl=n-k
       do i=1,n-k
          nl=n-k-i+1
          roffs=k-nl-1
          pl=>a(roffs+1:k,k+i:n)
          pq=>q(i:n-k,i:n-k)
          numrotsu(k+i-1)=nl
          ! Triangularize L
          do j=nl,1,-1
             rot=lgivens2(pl(j,j),pl(j+1,j))
             csu(j,k+i-1)=rot%cosine; ssu(j,k+i-1)=rot%sine
             jsu(j,k+i-1)=roffs+j
             call rotation_times_general(trp_rot(rot), pl,j,j+1)
             pl(j,j)=0.0_dp
          end do
          ! reveal column k + i
          do j=nl,2,-1
             rot=lgivens(pq(j-1,1),pq(j,1))
             call rotation_times_general(trp_rot(rot), pq, j-1,j)
             call general_times_rotation(pl, rot, j-1,j)
          end do
          pl(:,1)=pl(:,1)*pq(1,1)
       end do
    end if
  end subroutine f_d_general_ub

  !
  ! Complex.
  !
  ! Updating procedure to compute a UB factorization from a general matrix.
  !

  subroutine c_upper_to_ub(a,ub,lbw,tol,error)
    complex(kind=dp), target, dimension(:,:), intent(inout) :: a
    type(c_ub), intent(inout) :: ub
    type(error_info), intent(out) :: error
    real(kind=dp), intent(in) :: tol
    integer(kind=int32), intent(in) :: lbw
    call clear_error(error)
    if (size(a,1) < 1) then
       call set_error(error, 1, id_c_upper_to_ub); return
    end if
    if (get_lbwmax(ub) < lbw) then
       call set_error(error, 2, id_c_upper_to_ub); return
    end if
    if (get_n(ub) /= size(a,1) .or. get_n(ub) /= size(a,2)) then
       call set_error(error, 3, id_c_upper_to_ub); return
    end if
    call f_c_upper_to_ub(a,get_n(ub),ub%bc, lbw, ub%ubw, get_lbwmax(ub), get_ubwmax(ub), &
         ub%numrotsu, ub%jsu, ub%csu, ub%ssu, tol, error)
    ub%lbw=lbw
  end subroutine c_upper_to_ub


  subroutine f_c_upper_to_ub(a, n, b, lbw, ubw, lbwmax, ubwmax, &
       numrotsu, jsu, csu, ssu, tol, error)
    complex(kind=dp), target, dimension(n,n), intent(inout) :: a
    integer(kind=int32), dimension(ubwmax,n), intent(out) :: jsu
    complex(kind=dp), dimension(ubwmax,n), intent(out) :: csu, ssu
    complex(kind=dp), dimension(lbwmax+ubwmax+1,n), intent(out) :: b
    real(kind=dp), intent(in) :: tol
    integer(kind=int32), dimension(n), intent(out) :: numrotsu
    integer(kind=int32), intent(out) :: ubw
    type(error_info), intent(out) :: error
    integer(kind=int32), intent(in) :: n, lbw, lbwmax, ubwmax
    !
    integer(kind=int32), dimension(n) :: ubws
    !
    call clear_error(error)
    !
    if (n == 1) then
       numrotsu=0;
       ssu=(0.0_dp, 0.0_dp); csu=(0.0_dp, 0.0_dp); jsu=0
       ubw=0
       b(1,1)=a(1,1)
       return
    end if

    call f_c_general_ub(a, n, ubws, ubwmax, numrotsu, jsu, csu, ssu, tol, error)
    ubw=maxval(ubws)
    if (ubw > ubwmax) then
       ! This should already have been detected in f_d_general_ub.
       call set_error(error, 1, id_f_c_upper_to_ub)
    else
       call c_extract_diagonals_bc(a,n,b,lbw,ubw,lbwmax, ubwmax)
    end if
  end subroutine f_c_upper_to_ub

  subroutine f_c_general_ub(a, n, ubws, ubwmax, numrotsu, jsu, csu, ssu, tol, error)
    complex(kind=dp), target, dimension(n,n), intent(inout) :: a
    integer(kind=int32), dimension(ubwmax,n), intent(out) :: jsu
    complex(kind=dp), dimension(ubwmax,n), intent(out) :: csu, ssu
    real(kind=dp), intent(in) :: tol
    integer(kind=int32), dimension(n), intent(out) :: numrotsu
    integer(kind=int32), dimension(n), intent(out) :: ubws
    type(error_info), intent(out) :: error
    integer(kind=int32), intent(in) :: n, ubwmax
    !
    complex(kind=dp), target, dimension(ubwmax+1,n) :: q
    complex(kind=dp), dimension(ubwmax+1) :: x
    real(kind=dp) :: nrma, nrmq12
    complex(kind=dp), pointer, dimension(:,:) :: pl, pq
    integer(kind=int32) :: i, j, k, roffs, nl, klast
    type(c_rotation) :: rot

    call clear_error(error)
    q=(0.0_dp, 0.0_dp); numrotsu=0;
    ssu=(0.0_dp, 0.0_dp); csu=(0.0_dp, 0.0_dp); jsu=0
    ubws=0
    nrma = maxabs(a)*sqrt(real(n))

    if (n == 1) then
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
             if (ubwmax==0) then
                call set_error(error, 1, id_f_c_general_ub); return
             end if
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
                call shift2(pq,0,-1)
                if (k < n-2) then
                   pq => q(1:2,1:n-k-1)
                   pq(2,:)=a(k+1,k+2:n)
                   pl => a(k:k+1,k+2:k+3)
                   call extend_gs_rows(pq(1:1,:), pl(2,1:1), pl(2,2), pq(2,:), error)
                   if (error%code > 0) then
                      call add_id(error, id_f_c_general_ub); return
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
          call lower_left_nullvec(x(1:nl),pl,tol*nrma,nullmaxits,error)
          if (error%code <= 0) then ! if there is a left null vector then introduce a zero row.
             ubws(k)=nl-1
             if (error%code == -1) then
                pl(1,1)=(0.0_dp, 0.0_dp)
                numrotsu(k)=0
                call clear_error(error)
             else if (error%code < -1) then
                numrotsu(k)=-error%code-1
                pl(-error%code,-error%code)=(0.0_dp, 0.0_dp)
                do j=-error%code-1,1,-1
                   rot=lgivens2(pl(j,j),pl(j+1,j))
                   call rotation_times_general(trp_rot(rot), pl(:,1:j), j,j+1)
                   pl(j,j)=(0.0_dp, 0.0_dp)
                   csu(j,k)=rot%cosine; ssu(j,k)=rot%sine
                   jsu(j,k)=roffs+j
                end do
                call clear_error(error)
             else ! error=0
                numrotsu(k)=nl-1;
                do j=nl,2,-1 ! apply u_k while preserving the triangular structure of L
                   rot=rgivens(x(j-1),x(j))
                   call general_times_rotation(x,rot,j-1,j)
                   call rotation_times_general(trp_rot(rot), pl,j-1,j)
                   csu(j-1,k)=rot%cosine; ssu(j-1,k)=rot%sine
                   jsu(j-1,k)=roffs+j-1
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
                call shift2(pq,-1,-1)
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
                call extend_gs_rows(pq(2:nl,:), x(1:nl-1), x(nl), pq(1,:), error) ! orthogonalize
                if (error%code > 0) then
                   call add_id(error, id_f_c_general_ub); return
                end if
                do j=nl,2,-1
                   rot=lgivens(pq(1,1),pq(j,1))
                   call rotation_times_general(trp_rot(rot), pq, 1,j)
                   call general_times_rotation(pl(j:nl,:),rot,1,j)
                end do
                pl(:,1)=pl(:,1)*pq(1,1)
                call shift2(pq,-1,-1)
                ! extend the LQ factorization
                pq => q(1:nl,1:n-k-1)
                pq(nl,:)=a(k+1,k+2:n)
                pl => a(roffs+2:k+1,k+2:k+nl+1) ! nl by nl
                call extend_gs_rows(pq(1:nl-1,:), pl(nl,1:nl-1), pl(nl,nl), pq(nl,:), error)
                if (error%code > 0) then
                   call add_id(error, id_f_c_general_ub); return
                end if
                if (k+nl+2 <= n) then
                   a(k+1,k+nl+2:n)=(0.0_dp, 0.0_dp)
                end if
             end if
          else
             ! no null vector found.  Simply reveal column k+1 if there is room.  Otherwise terminate
             !  with square L.
             call clear_error(error)
             if (ubwmax<nl) then
                call set_error(error, 1, id_f_c_general_ub); return
             end if
             ubws(k)=nl
             if (k+nl==n) then
                klast=k ! terminate with square L.
                exit kloop
             else
                pl => a(roffs+1:k, k+1:k+nl+1) ! extend pl to the right. (note this requires k+nl < n)
                call shift2(pl,0,1)
                pq => q(1:nl+1, 1:n-k)
                call shift2(pq,1,0)
                pq(1,:)=(0.0_dp, 0.0_dp)
                pq(1,1)=(1.0_dp, 0.0_dp)
                ! orthogonalizing.  Note x is used as workspace for coefficients.
                call extend_gs_rows(pq(2:nl+1,:), x(1:nl), x(nl+1), pq(1,:), error)
                if (error%code > 0) then
                   call add_id(error, id_f_c_general_ub); return
                end if
                do j=nl+1,2,-1
                   rot=lgivens(pq(1,1),pq(j,1))
                   call rotation_times_general(trp_rot(rot), pq, 1,j)
                   call general_times_rotation(pl(j-1:nl,:),rot,1,j)
                end do
                pl(:,1)=pl(:,1)*pq(1,1)
                call shift2(pq,-1,-1)
                if (nl==n-k-1) then ! q is now nl by nl.  Extend the LQ factorization down one row before stopping.
                   pq => q(1:nl,1:nl)
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
                   call extend_gs_rows(pq(1:nl,:), pl(nl+1,1:nl), pl(nl+1,nl+1), pq(nl+1,:), error)
                   if (error%code > 0) then
                      call add_id(error, id_f_c_general_ub); return
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
       if (ubwmax<n-k) then
          call set_error(error, 1, id_f_d_general_ub); return
       end if
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
             numrotsu(k+i)=nl-1
             do j=nl-1,1,-1
                rot=lgivens2(pl(j,j+1),pl(j+1,j+1))
                call rotation_times_general(trp_rot(rot),pl(:,2:nl),j,j+1)
                csu(j,k+i)=rot%cosine; ssu(j,k+i)=rot%sine
                jsu(j,k+i)=roffs+j
                pl(j,j+1)=(0.0_dp, 0.0_dp)
             end do
          end do
          pl(2,2)=pl(2,2)*pq(2,2)
       end if
    else ! rectangular termination
       k=k+1
       if (ubwmax < n-k) then
          call set_error(error, 1, id_f_d_general_ub); return
       end if
       ubws(k:n-1)=n-k
       nl=n-k
       do i=1,n-k
          nl=n-k-i+1
          roffs=k-nl-1
          pl=>a(roffs+1:k,k+i:n)
          pq=>q(i:n-k,i:n-k)
          numrotsu(k+i-1)=nl
          ! Triangularize L
          do j=nl,1,-1
             rot=lgivens2(pl(j,j),pl(j+1,j))
             csu(j,k+i-1)=rot%cosine; ssu(j,k+i-1)=rot%sine
             jsu(j,k+i-1)=roffs+j
             call rotation_times_general(trp_rot(rot), pl,j,j+1)
             pl(j,j)=(0.0_dp, 0.0_dp)
          end do
          ! reveal column k + i
          do j=nl,2,-1
             rot=lgivens(pq(j-1,1),pq(j,1))
             call rotation_times_general(trp_rot(rot), pq, j-1,j)
             call general_times_rotation(pl, rot, j-1,j)
          end do
          pl(:,1)=pl(:,1)*pq(1,1)
       end do
    end if
  end subroutine f_c_general_ub

end module general_ub
