module mod_general_bv
  use mod_prec
  use mod_error_id
  use mod_utility
  use mod_gs
  use mod_nullvec
  use mod_shift
  use mod_rotation
  use mod_orth_band_types
  use mod_band_types
  implicit none
  integer(kind=int32), parameter :: nullmaxits=5

  private

  public :: general_to_bv, d_general_to_bv, c_general_to_bv, &
       f_general_to_bv, f_d_general_to_bv, f_c_general_to_bv, &
       f_general_bv, f_d_general_bv, f_c_general_bv, &
       d_bv_of_general, c_bv_of_general, bv_of_general, bv

  interface bv_of_general
     module procedure d_bv_of_general, c_bv_of_general
  end interface bv_of_general

  interface bv
     module procedure d_bv_of_general, c_bv_of_general
  end interface bv

  interface general_to_bv
     module procedure d_general_to_bv, c_general_to_bv
  end interface general_to_bv

  interface f_general_to_bv
     module procedure f_d_general_to_bv, f_c_general_to_bv
  end interface f_general_to_bv

  interface f_general_bv
     module procedure f_d_general_bv, f_c_general_bv
  end interface f_general_bv

contains

  function d_bv_of_general(a, lbw, lbwmax, ubwmax, tol, error) result(bv)
    type(d_bv), allocatable :: bv
    real(kind=dp), target, dimension(:,:), intent(inout) :: a
    type(error_info), intent(inout), optional :: error
    real(kind=dp), intent(in) :: tol
    integer(kind=int32), intent(in) :: lbw, lbwmax, ubwmax
    type(routine_info), parameter :: info=info_d_bv_of_general
    integer(kind=int32) :: n

    if (failure(error)) return
    call push_id(info, error)

    n=size(a,1)
    bv=d_new_bv(n,lbwmax,ubwmax)

    call d_general_to_bv(a,bv,lbw,tol,error)

    call pop_id(error)
    
  end function d_bv_of_general


  ! Errors:
  ! 0: no error
  ! 1: n<1
  ! 2: bv%lbwmax < lbw
  ! 3: n is not the same for a and bv.
  subroutine d_general_to_bv(a,bv,lbw,tol,error)
    real(kind=dp), target, dimension(:,:), intent(inout) :: a
    type(d_bv), intent(inout) :: bv
    type(error_info), intent(inout), optional :: error
    real(kind=dp), intent(in) :: tol
    integer(kind=int32), intent(in) :: lbw
    type(routine_info), parameter :: info=info_d_general_to_bv

    if (failure(error)) return
    call push_id(info, error)
    
    if (size(a,1) < 1) then
       call set_error(1, info, error); return
    end if
    if (get_lbwmax(bv) < lbw) then
       call set_error(2, info, error); return
    end if
    if (get_n(bv) /= size(a,1) .or. get_n(bv) /= size(a,2)) then
       call set_error(3, info, error); return
    end if
    call f_d_general_to_bv(a,get_n(bv),bv%br, lbw, bv%ubw, get_lbwmax(bv), get_ubwmax(bv), &
         bv%numrotsv, bv%ksv, bv%csv, bv%ssv, tol, error)
    bv%lbw=lbw

    call pop_id(error)
    
  end subroutine d_general_to_bv

  ! Errors:
  ! 0: no error
  ! 1: insufficient upper bw in bv%br
  subroutine f_d_general_to_bv(a, n, b, lbw, ubw, lbwmax, ubwmax, &
       numrotsv, ksv, csv, ssv, tol, error)
    real(kind=dp), target, dimension(n,n), intent(inout) :: a
    integer(kind=int32), dimension(n,ubwmax), intent(out) :: ksv
    real(kind=dp), dimension(n,ubwmax), intent(out) :: csv, ssv
    real(kind=dp), dimension(n,lbwmax+ubwmax+1), intent(out) :: b
    real(kind=dp), intent(in) :: tol
    integer(kind=int32), dimension(n), intent(out) :: numrotsv
    integer(kind=int32), intent(out) :: ubw
    type(error_info), optional :: error
    integer(kind=int32), intent(in) :: n, lbw, lbwmax, ubwmax
    !
    integer(kind=int32), dimension(n) :: ubws
    type(routine_info), parameter :: info=info_f_d_general_to_bv
    !
    if (failure(error)) return
    call push_id(info, error)

    if (n == 1) then
       numrotsv=0
       ssv=0.0_dp; csv=0.0_dp; ksv=0
       ubw=0
       b(1,1)=a(1,1)
       return
    end if

    call f_d_general_bv(a, n, ubws, ubwmax, numrotsv, ksv, csv, ssv, tol, error)

    if (success(error)) then
       ubw=maxval(ubws)
       if (ubw > ubwmax) then
          call set_error(1, info, error); return
       else
          call d_extract_diagonals_br(a, n, b, lbw, ubw, lbwmax, ubwmax)
       end if
    end if

    call pop_id(error)
  end subroutine f_d_general_to_bv

  ! Errors:
  ! 0: no error
  ! 1: insufficient upper bw in bv%br
  subroutine f_d_general_bv(a, n, ubws, ubwmax, &
       numrotsv, ksv, csv, ssv, tol, error)
    real(kind=dp), target, dimension(n,n), intent(inout) :: a
    integer(kind=int32), dimension(n,ubwmax), intent(out) :: ksv
    real(kind=dp), dimension(n,ubwmax), intent(out) :: csv, ssv
    real(kind=dp), intent(in) :: tol
    integer(kind=int32), dimension(n), intent(out) :: numrotsv
    integer(kind=int32), dimension(n), intent(out) :: ubws
    type(error_info), optional :: error
    integer(kind=int32), intent(in) :: n, ubwmax
    !
    real(kind=dp), target, dimension(n,ubwmax+1) :: q
    real(kind=dp), dimension(ubwmax+1) :: x
    real(kind=dp) :: nrma
    real(kind=dp), pointer, dimension(:,:) :: pl, pq
    integer(kind=int32) :: i, j, k, roffs, coffs, nl, p
    type(d_rotation) :: rot
    type(routine_info), parameter :: info=info_f_d_general_bv
    type(error_info) :: errornv
    logical :: null
    !
    if (failure(error)) return
    call push_id(info, error)
    
    q=0.0_dp; numrotsv=0;
    ssv=0.0_dp; csv=0.0_dp; ksv=0
    ubws=0
    nrma = maxabs(a)*sqrt(real(n))
    !
    if (n == 1) then
       return
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
    ! k is the trailing principal submatrix size.
    kloop: do k=1,n-1
       ! Current, possibly singular, L should be contained in
       ! a(n-k-nl+1:n-k,n-k+1:n-k+nl)
       ! At the start of the loop nl is a possible overestimate that
       ! can be decreased by finding a null vector.
       roffs=n-k-nl
       coffs=n-k
       pl => a(roffs+1:roffs+nl, coffs+1:coffs+nl)
       pq => q(1:n-k,1:nl)
       call clear_error(errornv)
       call clear_routines(errornv)
       errornv%halt=.false.
       call lower_right_nullvec(x(1:nl),pl,tol*nrma,nullmaxits,p,errornv)
       if (success(errornv)) then
          null=.true.
          ! if there is a left null vector then introduce a zero row.
          ubws(k)=nl-1
          if (p >= 1) then
             numrotsv(n-k)=nl-p
             pl(p,p)=0.0_dp
             do j=p,nl-1
                rot=rgivens(pl(j+1,j),pl(j+1,j+1))
                call general_times_rotation(pl(j+1:nl,:),rot,j,j+1)
                pl(j+1,j+1)=0.0_dp
                csv(n-k,nl-j)=rot%cosine; ssv(n-k,nl-j)=rot%sine
                ksv(n-k,nl-j)=coffs+j
             end do
          else ! p==0
             numrotsv(n-k)=nl-1;
             do j=2,nl ! apply v_k while preserving the triangular structure of L
                rot=lgivens2(x(j-1),x(j))
                call rotation_times_general(trp_rot(rot),x, j-1,j)
                call general_times_rotation(pl,rot,j-1,j)
                csv(n-k,nl-j+1)=rot%cosine; ssv(n-k,nl-j+1)=rot%sine
                ksv(n-k,nl-j+1)=coffs+j-1
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
             exit kloop ! terminate
          else
             pq(:,nl)=0.0_dp; pq(n-k,nl)=1.0_dp
             call extend_gs_columns(pq(:,1:nl-1), x(1:nl-1), x(nl), pq(:,nl), error)
             if (failure(error)) return
             do j=1,nl-1
                rot=rgivens2(pq(n-k,j),pq(n-k,nl))
                call general_times_rotation(pq,rot,j,nl)
                call rotation_times_general(trp_rot(rot),pl(:,1:j),j,nl)
             end do
             pl(nl,:)=pq(n-k,nl)*pl(nl,:)
             ! extend the LQ factorization with a(1:n-k-1,coffs)
             call shift(pq,0,1)
             pq => q(1:n-k-1,1:nl)
             pq(:,1)=a(1:n-k-1,coffs)
             pl => a(roffs:roffs+nl-1, coffs:coffs+nl-1)
             call extend_gs_columns(pq(:,2:nl),pl(2:nl,1),pl(1,1),pq(:,1),error)
             if (failure(error)) return
             a(1:roffs-1,coffs)=0.0_dp
          end if
       else
          ! no null vector found.  Simply reveal row nl if there is room.
          ! Otherwise terminate with square L
          null=.false.
          if (ubwmax < nl) then
             call set_error(1, info, error); return
          end if
          ubws(k)=nl
          if (k+nl == n) then
             exit kloop
          else
             pl => a(roffs:roffs+nl, coffs+1:coffs+nl) ! extend pl up
             call shift(pl,-1,0)
             pq => q(1:n-k,1:nl+1)
             pq(:,nl+1)=0.0_dp; pq(n-k,nl+1)=1.0_dp
             call extend_gs_columns(pq(:,1:nl),x(1:nl), x(nl+1),pq(:,nl+1),error)
             if (failure(error)) return
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
                null=.true.
                exit kloop
             else ! q is not square.  Make L (nl+1)x(nl+1)
                pq => q(1:n-k-1,1:nl+1)
                call shift(pq,0,1)
                pq(:,1)=a(1:n-k-1,coffs)
                pl => a(roffs-1:roffs-1+nl,coffs:coffs+nl)
                call extend_gs_columns(pq(:,2:nl+1), pl(2:nl+1,1), &
                     pl(1,1), pq(:,1), error)
                if (failure(error)) return
                a(1:roffs-2,coffs)=0.0_dp
                nl=nl+1
             end if
          end if
       end if ! null vector check
    end do kloop
    ! If null=.true., we need to terminate on an nl by nl+1 matrix L
    ! contained in a(1:nl,n-k:n-k+nl)
    ! If null=.false., then L is square and
    ! contained in a(1:n-k,n-k+1:n-k+(n-k)).  The former happens when the last step
    ! of kloop found a null vector.  The latter happens when it didn't.
    if (.not. null) then ! square termination
       if (ubwmax<n-k) then
          call set_error(1, info, error); return
       end if
       ubws(k:n-1)=n-k       
       nl = n-k
       coffs=n-k
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
          numrotsv(n-(k+i))=nl-1
          do j=1,nl-1
             rot=rgivens(pl(j,j),pl(j,j+1))
             call general_times_rotation(pl(j:nl-1,:),rot,j,j+1)
             pl(j,j+1)=0.0_dp
             csv(n-(k+i),nl-j)=rot%cosine; ssv(n-(k+i),nl-j)=rot%sine
             ksv(n-(k+i),nl-j)=coffs+j
          end do
       end do
       pl(1,1)=pq(1,1)*pl(1,1)
    else ! rectangular termination
       k=k+1
       if (ubwmax<n-k) then
          call set_error(1, info, error); return
       end if
       ubws(k:n-1)=n-k
       nl=n-k
       coffs=n-k
       do i=1,n-k
          nl=n-k-i+1
          pl=>a(1:nl,coffs+1:coffs+nl+1)
          pq=>q(1:nl,1:nl)
          numrotsv(n-(k+i-1))=nl
          do j=1,nl
             rot=rgivens(pl(j,j),pl(j,j+1))
             csv(n-(k+i-1),nl-j+1)=rot%cosine; ssv(n-(k+i-1),nl-j+1)=rot%sine
             ksv(n-(k+i-1),nl-j+1)=coffs+j
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
    call pop_id(error)
  end subroutine f_d_general_bv

  ! complex BV

  function c_bv_of_general(a, lbw, lbwmax, ubwmax, tol, error) result(bv)
    type(c_bv), allocatable :: bv
    complex(kind=dp), target, dimension(:,:), intent(inout) :: a
    type(error_info), intent(inout), optional :: error
    real(kind=dp), intent(in) :: tol
    integer(kind=int32), intent(in) :: lbw, lbwmax, ubwmax
    type(routine_info), parameter :: info=info_c_bv_of_general
    integer(kind=int32) :: n

    if (failure(error)) return
    call push_id(info, error)

    n=size(a,1)
    bv=c_new_bv(n,lbwmax,ubwmax)

    call c_general_to_bv(a,bv,lbw,tol,error)

    call pop_id(error)
    
  end function c_bv_of_general

  subroutine c_general_to_bv(a,bv,lbw,tol,error)
    complex(kind=dp), target, dimension(:,:), intent(inout) :: a
    type(c_bv), intent(inout) :: bv
    type(error_info), intent(inout), optional :: error
    real(kind=dp), intent(in) :: tol
    integer(kind=int32), intent(in) :: lbw
    type(routine_info), parameter :: info=info_c_general_to_bv
    
    if (failure(error)) return
    call push_id(info, error)
    !
    if (size(a,1) < 1) then
       call set_error(1, info, error); return
    end if
    if (get_lbwmax(bv) < lbw) then
       call set_error(2, info, error); return
    end if
    if (get_n(bv) /= size(a,1) .or. get_n(bv) /= size(a,2)) then
       call set_error(3, info, error); return
    end if
    call f_c_general_to_bv(a,get_n(bv),bv%br, lbw, bv%ubw, &
         get_lbwmax(bv), get_ubwmax(bv), &
         bv%numrotsv, bv%ksv, bv%csv, bv%ssv, tol, error)

    bv%lbw=lbw
    call pop_id(error)
    
  end subroutine c_general_to_bv

  subroutine f_c_general_to_bv(a, n, b, lbw, ubw, lbwmax, ubwmax, &
       numrotsv, ksv, csv, ssv, tol, error)
    complex(kind=dp), target, dimension(n,n), intent(inout) :: a
    integer(kind=int32), dimension(n,ubwmax), intent(out) :: ksv
    real(kind=dp), dimension(n,ubwmax), intent(out) :: csv
    complex(kind=dp), dimension(n,ubwmax), intent(out) :: ssv
    complex(kind=dp), dimension(n,lbwmax+ubwmax+1), intent(out) :: b
    real(kind=dp), intent(in) :: tol
    integer(kind=int32), dimension(n), intent(out) :: numrotsv
    integer(kind=int32), intent(out) :: ubw
    type(error_info), intent(inout), optional :: error
    integer(kind=int32), intent(in) :: n, lbw, lbwmax, ubwmax
    !
    integer(kind=int32), dimension(n) :: ubws
    type(routine_info), parameter :: info=info_f_c_general_to_bv
    !

    if (failure(error)) return
    call push_id(info, error)

    if (n == 1) then
       numrotsv=0
       ssv=(0.0_dp,0.0_dp); csv=0.0_dp; ksv=0
       ubw=0
       b(1,1)=a(1,1)
       return
    end if

    call f_c_general_bv(a, n, ubws, ubwmax, numrotsv, ksv, csv, ssv, tol, error)

    if (success(error)) then
       ubw=maxval(ubws)
       if (ubw > ubwmax) then
          call set_error(2, info, error); return
       else
          call c_extract_diagonals_br(a,n,b,lbw,ubw,lbwmax, ubwmax)
       end if
    end if
    call pop_id(error)
  end subroutine f_c_general_to_bv

  subroutine f_c_general_bv(a, n, ubws, ubwmax, &
       numrotsv, ksv, csv, ssv, tol, error)
    complex(kind=dp), target, dimension(n,n), intent(inout) :: a
    integer(kind=int32), dimension(n,ubwmax), intent(out) :: ksv
    real(kind=dp), dimension(n,ubwmax), intent(out) :: csv
    complex(kind=dp), dimension(n,ubwmax), intent(out) :: ssv
    real(kind=dp), intent(in) :: tol
    integer(kind=int32), dimension(n), intent(out) :: ubws
    integer(kind=int32), dimension(n), intent(out) :: numrotsv
    type(error_info), intent(inout), optional :: error
    integer(kind=int32), intent(in) :: n, ubwmax
    !
    complex(kind=dp), target, dimension(n,ubwmax+1) :: q
    complex(kind=dp), dimension(ubwmax+1) :: x
    real(kind=dp) :: nrma
    complex(kind=dp), pointer, dimension(:,:) :: pl, pq
    integer(kind=int32) :: i, j, k, roffs, coffs, nl, p
    type(c_rotation) :: rot
    type(routine_info), parameter :: info=info_f_c_general_bv
    type(error_info) :: errornv
    logical :: null
    !
    q=(0.0_dp, 0.0_dp); numrotsv=0;
    ssv=(0.0_dp, 0.0_dp); csv=0.0_dp; ksv=0
    ubws=0
    if (failure(error)) return
    call push_id(info, error)

    nrma = maxabs(a)*sqrt(real(n))

    if (n == 1) then
       return
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
    kloop: do k=1,n-1
       ! Current, possibly singular, L should be contained in
       ! a(n-k-nl+1:n-k,n-k+1:n-k+nl)
       ! At the start of the loop nl is a possible overestimate that
       ! can be decreased by finding a null vector.
       roffs=n-k-nl
       coffs=n-k
       pl => a(roffs+1:roffs+nl, coffs+1:coffs+nl)
       pq => q(1:n-k,1:nl)
       call clear_error(errornv)
       call clear_routines(errornv)
       errornv%halt=.false.
       call lower_right_nullvec(x(1:nl),pl,tol*nrma,nullmaxits,p,errornv)
       if (success(errornv)) then
          ! if there is a left null vector then introduce a zero row.
          null=.true.
          ubws(k)=nl-1
          if (p >= 1) then
             numrotsv(n-k)=nl-p
             pl(p,p)=(0.0_dp, 0.0_dp)
             do j=p,nl-1
                rot=rgivens(pl(j+1,j),pl(j+1,j+1))
                call general_times_rotation(pl(j+1:nl,:),rot,j,j+1)
                pl(j+1,j+1)=(0.0_dp, 0.0_dp)
                csv(n-k,nl-j)=rot%cosine; ssv(n-k,nl-j)=rot%sine
                ksv(n-k,nl-j)=coffs+j
             end do
          else ! p==0
             numrotsv(n-k)=nl-1;
             do j=2,nl
                ! apply v_k while preserving the triangular structure of L
                rot=lgivens2(x(j-1),x(j))
                call rotation_times_general(trp_rot(rot),x, j-1,j)
                call general_times_rotation(pl,rot,j-1,j)
                csv(n-k,nl-j+1)=rot%cosine; ssv(n-k,nl-j+1)=rot%sine
                ksv(n-k,nl-j+1)=coffs+j-1
                rot=lgivens2(pl(j-1,j),pl(j,j))
                call rotation_times_general(trp_rot(rot), pl(:,1:j), j-1,j)
                call general_times_rotation(pq,rot,j-1,j)
                pl(j-1,j)=(0.0_dp, 0.0_dp)
             end do
             pl(nl,nl)=(0.0_dp, 0.0_dp)
          end if
          do j=nl-1,1,-1
             ! compress (TODO: this is not necessary if zero is in pl(1,1))
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
             exit kloop ! terminate
          else
             pq(:,nl)=(0.0_dp, 0.0_dp); pq(n-k,nl)=(1.0_dp, 0.0_dp)
             call extend_gs_columns(pq(:,1:nl-1), x(1:nl-1), x(nl), pq(:,nl), error)
             if (failure(error)) return
             do j=1,nl-1
                rot=rgivens2(pq(n-k,j),pq(n-k,nl))
                call general_times_rotation(pq,rot,j,nl)
                call rotation_times_general(trp_rot(rot),pl(:,1:j),j,nl)
             end do
             pl(nl,:)=pq(n-k,nl)*pl(nl,:)
             ! extend the LQ factorization with a(1:n-k-1,coffs)
             call shift(pq,0,1)
             pq => q(1:n-k-1,1:nl)
             pq(:,1)=a(1:n-k-1,coffs)
             pl => a(roffs:roffs+nl-1, coffs:coffs+nl-1)
             call extend_gs_columns(pq(:,2:nl),pl(2:nl,1),pl(1,1),pq(:,1),error)
             if (failure(error)) return
             a(1:roffs-1,coffs)=(0.0_dp, 0.0_dp)
          end if
       else
          ! no null vector found.  Simply reveal row nl if there is room.
          ! Otherwise terminate with square L
          null=.false.
          if (ubwmax < nl) then
             call set_error(1, info, error); return
          end if
          ubws(k)=nl
          if (k+nl == n) then
             exit kloop
          else
             pl => a(roffs:roffs+nl, coffs+1:coffs+nl) ! extend pl up
             call shift(pl,-1,0)
             pq => q(1:n-k,1:nl+1)
             pq(:,nl+1)=(0.0_dp, 0.0_dp); pq(n-k,nl+1)=(1.0_dp, 0.0_dp)
             call extend_gs_columns(pq(:,1:nl),x(1:nl), x(nl+1),pq(:,nl+1),error)
             if (failure(error)) return
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
                null=.true.
                exit kloop
             else ! q is not square.  Make L (nl+1)x(nl+1)
                pq => q(1:n-k-1,1:nl+1)
                call shift(pq,0,1)
                pq(:,1)=a(1:n-k-1,coffs)
                pl => a(roffs-1:roffs-1+nl,coffs:coffs+nl)
                call extend_gs_columns(pq(:,2:nl+1), pl(2:nl+1,1), &
                     pl(1,1), pq(:,1), error)
                if (failure(error)) return
                a(1:roffs-2,coffs)=(0.0_dp, 0.0_dp)
                nl=nl+1
             end if
          end if
       end if ! null vector check
    end do kloop
    ! If null=.true., we need to terminate on an nl by nl+1 matrix L
    ! contained in a(1:nl,n-k:n-k+nl)
    ! If null=.false., then L is square and
    ! contained in a(1:n-k,n-k+1:n-k+(n-k)).  The former happens when the last step
    ! of kloop found a null vector.  The latter happens when it didn't.
    if (.not. null) then ! square termination
       if (ubwmax < n-k) then
          call set_error(1, info, error); return
       end if
       ubws(k:n-1)=n-k
       nl = n-k
       coffs=n-k
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
          numrotsv(n-(k+i))=nl-1
          do j=1,nl-1
             rot=rgivens(pl(j,j),pl(j,j+1))
             call general_times_rotation(pl(j:nl-1,:),rot,j,j+1)
             pl(j,j+1)=(0.0_dp, 0.0_dp)
             csv(n-(k+i),nl-j)=rot%cosine; ssv(n-(k+i),nl-j)=rot%sine
             ksv(n-(k+i),nl-j)=coffs+j
          end do
       end do
       pl(1,1)=pq(1,1)*pl(1,1)
    else ! rectangular termination
       k=k+1
       if (ubwmax < n-k) then
          call set_error(1, info, error); return
       end if
       ubws(k:n-1)=n-k
       nl=n-k
       coffs=n-k
       do i=1,n-k
          nl=n-k-i+1
          pl=>a(1:nl,coffs+1:coffs+nl+1)
          pq=>q(1:nl,1:nl)
          numrotsv(n-(k+i-1))=nl
          do j=1,nl
             rot=rgivens(pl(j,j),pl(j,j+1))
             csv(n-(k+i-1),nl-j+1)=rot%cosine; ssv(n-(k+i-1),nl-j+1)=rot%sine
             ksv(n-(k+i-1),nl-j+1)=coffs+j
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
    call pop_id(error)
  end subroutine f_c_general_bv

end module mod_general_bv
