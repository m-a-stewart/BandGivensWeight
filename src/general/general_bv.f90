module mod_general_bv
  use mod_prec
  use mod_error_id
  use mod_utility
  use mod_gs
  use mod_cond_triangular
  use mod_shift
  use mod_rotation
  use mod_orth_band_types
  use mod_band_types
  implicit none
  integer(kind=int32), parameter :: nullmaxits=5

  private

  public :: general_to_bv, d_general_to_bv, z_general_to_bv, &
       f_general_to_bv, f_d_general_to_bv, f_z_general_to_bv, &
       f_general_bv, f_d_general_bv, f_z_general_bv, &
       d_bv_of_general, z_bv_of_general, bv_of_general, bv

  interface bv_of_general
     module procedure d_bv_of_general, z_bv_of_general
  end interface bv_of_general

  interface bv
     module procedure d_bv_of_general, z_bv_of_general
  end interface bv

  interface general_to_bv
     module procedure d_general_to_bv, z_general_to_bv
  end interface general_to_bv

  interface f_general_to_bv
     module procedure f_d_general_to_bv, f_z_general_to_bv
  end interface f_general_to_bv

  interface f_general_bv
     module procedure f_d_general_bv, f_z_general_bv
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
    real(kind=dp), pointer, dimension(:,:) :: pl, pq
    integer(kind=int32) :: i, j, k, roffs, coffs, nl, p
    type(d_rotation) :: rot
    type(routine_info), parameter :: info=info_f_d_general_bv
    type(error_info) :: errornv
    integer(kind=int32) :: coffs1, k1, ubw1, ml, ml1
    !
    if (failure(error)) return
    call push_id(info, error)
    
    q=0.0_dp; numrotsv=0;
    ssv=0.0_dp; csv=0.0_dp; ksv=0
    ubws=0
    !
    if (n == 1) return
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
    ! k is the leading principal submatrix size.
    leading_loop: do k=n-1,1,-1
       ! Current L should be contained in
       ! a(k-nl+1:k,k+1:k+nl) or a(roffs+1:roffs+nl,coffs+1:coffs+nl)
       roffs=k-nl
       coffs=k
       pl => a(roffs+1:roffs+nl, coffs+1:coffs+nl)
       pq => q(1:k,1:nl)
       call clear_error(errornv)
       call clear_routines(errornv)
       errornv%halt=.false.
       call lower_right_nullvec(x(1:nl),pl,tol,nullmaxits,p,errornv)
       if (success(errornv)) then
          ! if there is a right null vector then introduce a zero column.
          ubws(k)=nl-1
          ! null vec in x(p:nl)
          numrotsv(k)=nl-p;
          do j=p+1,nl ! apply v_k while preserving the triangular structure of L
             rot=lgivens2(x(j-1),x(j))
             call rotation_times_general(trp_rot(rot),x, j-1,j)
             call general_times_rotation(pl,rot,j-1,j)
             csv(k,nl-j+1)=rot%cosine; ssv(k,nl-j+1)=rot%sine
             ksv(k,nl-j+1)=coffs+j-1
             rot=lgivens2(pl(j-1,j),pl(j,j))
             call rotation_times_general(trp_rot(rot), pl(:,1:j), j-1,j)
             call general_times_rotation(pq,rot,j-1,j)
             pl(j-1,j)=0.0_dp
          end do
          pl(nl,nl)=0.0_dp

          do j=nl-1,1,-1 ! compress
             rot=lgivens(pl(j,j),pl(nl,j))
             call rotation_times_general(trp_rot(rot),pl(:,1:j),j,nl)
             call general_times_rotation(pq,rot,j,nl)
             pl(nl,j)=0.0_dp
          end do

          ! now L is (nl-1) x (nl-1) with a zero row stored
          ! directly below it in a and q is k x (nl-1).
          ! reveal row k of a of into the zero row below L.
          pq(:,nl)=0.0_dp; pq(k,nl)=1.0_dp
          call extend_gs_columns(pq(:,1:nl-1), x(1:nl-1), x(nl), pq(:,nl), error)
          if (failure(error)) return
          do j=1,nl-1
             rot=rgivens2(pq(k,j),pq(k,nl))
             call general_times_rotation(pq,rot,j,nl)
             call rotation_times_general(trp_rot(rot),pl(:,1:j),j,nl)
          end do
          pl(nl,:)=pq(k,nl)*pl(nl,:)
          ! now L is (nl-1) x (nl-1) stored in
          ! A(k-nl+1:k-1,k+1:k+nl-1) and Q is (k-1) x (nl-1).
          if (k==nl) then          
             ! q is square and (nl-1) x (nl-1)
             ! extend the QL decomposition by one column and exit
             ! with (nl-1) x nl L contained in
             ! a(1:nl-1,nl:2*nl-1) or
             x(1:nl-1)=matmul(transpose(pq(1:nl-1,1:nl-1)),a(1:nl-1,coffs))
             a(1:nl-1,coffs)=x(1:nl-1)
             k1 = k-1 ! v_{k-1} is next to be applied
             coffs1 = nl-1
             ubw1=nl-1
             ml1=nl-1
             exit leading_loop ! terminate
          else
             ! Q is rectangular and (k-1) x (nl-1).
             ! Extend the QL decomposition with an additional column
             ! A(1:k-1,coffs) so that L is nl x nl and in
             ! A(k-nl:k-1,k:k+nl-1) and q is (k-1) x nl
             call shift(pq,0,1)
             pq => q(1:k-1,1:nl)
             pq(:,1)=a(1:k-1,coffs)
             pl => a(roffs:roffs+nl-1, coffs:coffs+nl-1)
             call extend_gs_columns(pq(:,2:nl),pl(2:nl,1),pl(1,1),pq(:,1),error)
             if (failure(error)) return
             a(1:roffs-1,coffs)=0.0_dp
          end if
       else
          ! no null vector found.
          ! L is nl x nl in a(k-nl+1:k,k+1:k+nl) and
          ! q is k x nl.
          if (ubwmax < nl) then
             call set_error(1, info, error); return
          end if
          ubws(k)=nl
          if (k==nl) then          
             ! q is square and nl x nl.  L is nl x nl
             ! and in A(1:nl,nl+1:2*nl).
             ! Downdate row nl of A so that
             ! L is nl-1 x nl
             do j=2,nl
                rot=rgivens2(pq(nl,j-1),pq(nl,j))
                call general_times_rotation(pq,rot,j-1,j)
                call rotation_times_general(trp_rot(rot),pl,j-1,j)
             end do
             pl(nl,:)=pq(nl,nl)*pl(nl,:)
             k1=k-1 ! v_{k-1} is the next transformation.
             coffs1=nl
             ubw1=nl
             ml1=nl-1
             exit leading_loop
          end if

          ! Downdate a row of A so that L is nl x nl and
          ! stored in A(k-nl:k-1)
          pl => a(roffs:roffs+nl, coffs+1:coffs+nl) ! extend pl up
          call shift(pl,-1,0)
          pq => q(1:k,1:nl+1)
          pq(:,nl+1)=0.0_dp; pq(k,nl+1)=1.0_dp
          call extend_gs_columns(pq(:,1:nl),x(1:nl), x(nl+1),pq(:,nl+1),error)
          if (failure(error)) return
          do j=1,nl
             rot=rgivens2(pq(k,j),pq(k,nl+1))
             call general_times_rotation(pq,rot,j,nl+1)
             call rotation_times_general(trp_rot(rot),pl(:,1:j),j,nl+1)
          end do
          pl(nl+1,:)=pl(nl+1,:)*pq(k,nl+1)
          if (nl==k-1) then
             ! q is now nl x nl.  Extend the QL factorization left
             ! one column before stopping so that l is nl x (nl+1)
             ! and stored in 
             pq => q(1:nl,1:nl)
             x(1:nl)=matmul(transpose(pq),a(1:nl,coffs))
             a(1:nl,coffs)=x(1:nl)
             k1=k-1 ! v_{k-1} is the next transformation
             coffs1=nl
             ubw1=nl
             ml1=nl
             exit leading_loop
          else ! q is not square.  Make L (nl+1)x(nl+1)
             pq => q(1:k-1,1:nl+1)
             call shift(pq,0,1)
             pq(:,1)=a(1:k-1,coffs)
             pl => a(roffs-1:roffs-1+nl,coffs:coffs+nl)
             call extend_gs_columns(pq(:,2:nl+1), pl(2:nl+1,1), &
                  pl(1,1), pq(:,1), error)
             if (failure(error)) return
             a(1:roffs-2,coffs)=0.0_dp
             nl=nl+1
          end if
       end if ! null vector check
    end do leading_loop

    ! There is an ml1 x ml1+1 lower trapezoidal matrix in
    ! A(1:ml1,coffs1+1:coffs1+ml1+1).  The upper bandwidth
    ! of all later blocks is ubw1.  The next transformation
    ! should be v_{k1}.
    if (ubwmax<ubw1) then
       call set_error(1, info, error); return
    end if
    ubws(1:k)=ubw1
    coffs=coffs1
    do i=1,ml1
       ml=ml1-i+1
       !       k=k1+i-1
       k=k1-i+1       
       pl=>a(1:ml,coffs+1:coffs+ml+1)
       pq=>q(1:ml,1:ml)
       numrotsv(k)=ml
       do j=1,ml
          rot=rgivens(pl(j,j),pl(j,j+1))
          csv(k,ml-j+1)=rot%cosine; ssv(k,ml-j+1)=rot%sine
          ksv(k,ml-j+1)=coffs+j
          call general_times_rotation(pl(j:ml,:),rot,j,j+1)
          pl(j,j+1)=0.0_dp
       end do
       ! reveal row ml
       do j=1,ml-1
          rot=rgivens2(pq(ml,j),pq(ml,j+1))
          call general_times_rotation(pq,rot,j,j+1)
          call rotation_times_general(trp_rot(rot),pl(:,1:j+1),j,j+1)
       end do
       pl(ml,:)=pq(ml,ml)*pl(ml,:)
    end do
    call pop_id(error)
  end subroutine f_d_general_bv

  ! complex BV

  function z_bv_of_general(a, lbw, lbwmax, ubwmax, tol, error) result(bv)
    type(z_bv), allocatable :: bv
    complex(kind=dp), target, dimension(:,:), intent(inout) :: a
    type(error_info), intent(inout), optional :: error
    real(kind=dp), intent(in) :: tol
    integer(kind=int32), intent(in) :: lbw, lbwmax, ubwmax
    type(routine_info), parameter :: info=info_z_bv_of_general
    integer(kind=int32) :: n

    if (failure(error)) return
    call push_id(info, error)

    n=size(a,1)
    bv=z_new_bv(n,lbwmax,ubwmax)

    call z_general_to_bv(a,bv,lbw,tol,error)

    call pop_id(error)
    
  end function z_bv_of_general

  subroutine z_general_to_bv(a,bv,lbw,tol,error)
    complex(kind=dp), target, dimension(:,:), intent(inout) :: a
    type(z_bv), intent(inout) :: bv
    type(error_info), intent(inout), optional :: error
    real(kind=dp), intent(in) :: tol
    integer(kind=int32), intent(in) :: lbw
    type(routine_info), parameter :: info=info_z_general_to_bv
    
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
    call f_z_general_to_bv(a,get_n(bv),bv%br, lbw, bv%ubw, &
         get_lbwmax(bv), get_ubwmax(bv), &
         bv%numrotsv, bv%ksv, bv%csv, bv%ssv, tol, error)

    bv%lbw=lbw
    call pop_id(error)
    
  end subroutine z_general_to_bv

  subroutine f_z_general_to_bv(a, n, b, lbw, ubw, lbwmax, ubwmax, &
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
    type(routine_info), parameter :: info=info_f_z_general_to_bv
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

    call f_z_general_bv(a, n, ubws, ubwmax, numrotsv, ksv, csv, ssv, tol, error)

    if (success(error)) then
       ubw=maxval(ubws)
       if (ubw > ubwmax) then
          call set_error(2, info, error); return
       else
          call z_extract_diagonals_br(a,n,b,lbw,ubw,lbwmax, ubwmax)
       end if
    end if
    call pop_id(error)
  end subroutine f_z_general_to_bv

  ! Errors:
  ! 0: no error
  ! 1: insufficient upper bw in bv%br
  subroutine f_z_general_bv(a, n, ubws, ubwmax, &
       numrotsv, ksv, csv, ssv, tol, error)
    complex(kind=dp), target, dimension(n,n), intent(inout) :: a
    integer(kind=int32), dimension(n,ubwmax), intent(out) :: ksv
    real(kind=dp), dimension(n,ubwmax), intent(out) :: csv
    complex(kind=dp), dimension(n,ubwmax), intent(out) :: ssv
    real(kind=dp), intent(in) :: tol
    integer(kind=int32), dimension(n), intent(out) :: numrotsv
    integer(kind=int32), dimension(n), intent(out) :: ubws
    type(error_info), optional :: error
    integer(kind=int32), intent(in) :: n, ubwmax
    !
    complex(kind=dp), target, dimension(n,ubwmax+1) :: q
    complex(kind=dp), dimension(ubwmax+1) :: x
    complex(kind=dp), pointer, dimension(:,:) :: pl, pq
    integer(kind=int32) :: i, j, k, roffs, coffs, nl, p
    type(z_rotation) :: rot
    type(routine_info), parameter :: info=info_f_z_general_bv
    type(error_info) :: errornv
    integer(kind=int32) :: coffs1, k1, ubw1, ml, ml1
    !
    if (failure(error)) return
    call push_id(info, error)
    
    q=(0.0_dp,0.0_dp); numrotsv=0;
    ssv=(0.0_dp,0.0_dp); csv=0.0_dp; ksv=0
    ubws=0
    !
    if (n == 1) return
    ! Compute an initial trivial QL factorization
    q(1:n-1,1) = a(1:n-1,n)
    a(1:n-1,n)=(0.0_dp,0.0_dp)
    a(n-1,n)=norm2(q(1:n-1,1))
    if (a(n-1,n) == (0.0_dp,0.0_dp)) then
       q(1:n-1,1) = (1.0_dp,0.0_dp)
       q(1:n-1,1)=q(1:n-1,1)/norm2(q(1:n-1,1))
    else
       q(1:n-1,1)=q(1:n-1,1)/a(n-1,n)
    end if
    nl=1
    ! k is the leading principal submatrix size.
    leading_loop: do k=n-1,1,-1
       ! Current L should be contained in
       ! a(k-nl+1:k,k+1:k+nl) or a(roffs+1:roffs+nl,coffs+1:coffs+nl)
       roffs=k-nl
       coffs=k
       pl => a(roffs+1:roffs+nl, coffs+1:coffs+nl)
       pq => q(1:k,1:nl)
       call clear_error(errornv)
       call clear_routines(errornv)
       errornv%halt=.false.
       call lower_right_nullvec(x(1:nl),pl,tol,nullmaxits,p,errornv)
       if (success(errornv)) then
          ! if there is a right null vector then introduce a zero column.
          ubws(k)=nl-1
          ! null vec in x(p:nl)
          numrotsv(k)=nl-p;
          do j=p+1,nl ! apply v_k while preserving the triangular structure of L
             rot=lgivens2(x(j-1),x(j))
             call rotation_times_general(trp_rot(rot),x, j-1,j)
             call general_times_rotation(pl,rot,j-1,j)
             csv(k,nl-j+1)=rot%cosine; ssv(k,nl-j+1)=rot%sine
             ksv(k,nl-j+1)=coffs+j-1
             rot=lgivens2(pl(j-1,j),pl(j,j))
             call rotation_times_general(trp_rot(rot), pl(:,1:j), j-1,j)
             call general_times_rotation(pq,rot,j-1,j)
             pl(j-1,j)=(0.0_dp,0.0_dp)
          end do
          pl(nl,nl)=(0.0_dp,0.0_dp)

          do j=nl-1,1,-1 ! compress
             rot=lgivens(pl(j,j),pl(nl,j))
             call rotation_times_general(trp_rot(rot),pl(:,1:j),j,nl)
             call general_times_rotation(pq,rot,j,nl)
             pl(nl,j)=(0.0_dp,0.0_dp)
          end do

          ! now L is (nl-1) x (nl-1) with a zero row stored
          ! directly below it in a and q is k x (nl-1).
          ! reveal row k of a of into the zero row below L.
          pq(:,nl)=(0.0_dp,0.0_dp); pq(k,nl)=(1.0_dp,0.0_dp)
          call extend_gs_columns(pq(:,1:nl-1), x(1:nl-1), x(nl), pq(:,nl), error)
          if (failure(error)) return
          do j=1,nl-1
             rot=rgivens2(pq(k,j),pq(k,nl))
             call general_times_rotation(pq,rot,j,nl)
             call rotation_times_general(trp_rot(rot),pl(:,1:j),j,nl)
          end do
          pl(nl,:)=pq(k,nl)*pl(nl,:)
          ! now L is (nl-1) x (nl-1) stored in
          ! A(k-nl+1:k-1,k+1:k+nl-1) and Q is (k-1) x (nl-1).
          if (k==nl) then          
             ! q is square and (nl-1) x (nl-1)
             ! extend the QL decomposition by one column and exit
             ! with (nl-1) x nl L contained in
             ! a(1:nl-1,nl:2*nl-1) or
             x(1:nl-1)=matmul(transpose(conjg(pq(1:nl-1,1:nl-1))),a(1:nl-1,coffs))
             a(1:nl-1,coffs)=x(1:nl-1)
             k1 = k-1 ! v_{k-1} is next to be applied
             coffs1 = nl-1
             ubw1=nl-1
             ml1=nl-1
             exit leading_loop ! terminate
          else
             ! Q is rectangular and (k-1) x (nl-1).
             ! Extend the QL decomposition with an additional column
             ! A(1:k-1,coffs) so that L is nl x nl and in
             ! A(k-nl:k-1,k:k+nl-1) and q is (k-1) x nl
             call shift(pq,0,1)
             pq => q(1:k-1,1:nl)
             pq(:,1)=a(1:k-1,coffs)
             pl => a(roffs:roffs+nl-1, coffs:coffs+nl-1)
             call extend_gs_columns(pq(:,2:nl),pl(2:nl,1),pl(1,1),pq(:,1),error)
             if (failure(error)) return
             a(1:roffs-1,coffs)=(0.0_dp,0.0_dp)
          end if
       else
          ! no null vector found.
          ! L is nl x nl in a(k-nl+1:k,k+1:k+nl) and
          ! q is k x nl.
          if (ubwmax < nl) then
             call set_error(1, info, error); return
          end if
          ubws(k)=nl
          if (k==nl) then          
             ! q is square and nl x nl.  L is nl x nl
             ! and in A(1:nl,nl+1:2*nl).
             ! Downdate row nl of A so that
             ! L is nl-1 x nl
             do j=2,nl
                rot=rgivens2(pq(nl,j-1),pq(nl,j))
                call general_times_rotation(pq,rot,j-1,j)
                call rotation_times_general(trp_rot(rot),pl,j-1,j)
             end do
             pl(nl,:)=pq(nl,nl)*pl(nl,:)
             k1=k-1 ! v_{k-1} is the next transformation.
             coffs1=nl
             ubw1=nl
             ml1=nl-1
             exit leading_loop
          end if

          ! Downdate a row of A so that L is nl x nl and
          ! stored in A(k-nl:k-1)
          pl => a(roffs:roffs+nl, coffs+1:coffs+nl) ! extend pl up
          call shift(pl,-1,0)
          pq => q(1:k,1:nl+1)
          pq(:,nl+1)=(0.0_dp,0.0_dp); pq(k,nl+1)=(1.0_dp,0.0_dp)
          call extend_gs_columns(pq(:,1:nl),x(1:nl), x(nl+1),pq(:,nl+1),error)
          if (failure(error)) return
          do j=1,nl
             rot=rgivens2(pq(k,j),pq(k,nl+1))
             call general_times_rotation(pq,rot,j,nl+1)
             call rotation_times_general(trp_rot(rot),pl(:,1:j),j,nl+1)
          end do
          pl(nl+1,:)=pl(nl+1,:)*pq(k,nl+1)
          if (nl==k-1) then
             ! q is now nl x nl.  Extend the QL factorization left
             ! one column before stopping so that l is nl x (nl+1)
             ! and stored in 
             pq => q(1:nl,1:nl)
             x(1:nl)=matmul(transpose(conjg(pq)),a(1:nl,coffs))
             a(1:nl,coffs)=x(1:nl)
             k1=k-1 ! v_{k-1} is the next transformation
             coffs1=nl
             ubw1=nl
             ml1=nl
             exit leading_loop
          else ! q is not square.  Make L (nl+1)x(nl+1)
             pq => q(1:k-1,1:nl+1)
             call shift(pq,0,1)
             pq(:,1)=a(1:k-1,coffs)
             pl => a(roffs-1:roffs-1+nl,coffs:coffs+nl)
             call extend_gs_columns(pq(:,2:nl+1), pl(2:nl+1,1), &
                  pl(1,1), pq(:,1), error)
             if (failure(error)) return
             a(1:roffs-2,coffs)=(0.0_dp,0.0_dp)
             nl=nl+1
          end if
       end if ! null vector check
    end do leading_loop

    ! There is an ml1 x ml1+1 lower trapezoidal matrix in
    ! A(1:ml1,coffs1+1:coffs1+ml1+1).  The upper bandwidth
    ! of all later blocks is ubw1.  The next transformation
    ! should be v_{k1}.
    if (ubwmax<ubw1) then
       call set_error(1, info, error); return
    end if
    ubws(1:k)=ubw1
    coffs=coffs1
    do i=1,ml1
       ml=ml1-i+1
       !       k=k1+i-1
       k=k1-i+1       
       pl=>a(1:ml,coffs+1:coffs+ml+1)
       pq=>q(1:ml,1:ml)
       numrotsv(k)=ml
       do j=1,ml
          rot=rgivens(pl(j,j),pl(j,j+1))
          csv(k,ml-j+1)=rot%cosine; ssv(k,ml-j+1)=rot%sine
          ksv(k,ml-j+1)=coffs+j
          call general_times_rotation(pl(j:ml,:),rot,j,j+1)
          pl(j,j+1)=(0.0_dp,0.0_dp)
       end do
       ! reveal row ml
       do j=1,ml-1
          rot=rgivens2(pq(ml,j),pq(ml,j+1))
          call general_times_rotation(pq,rot,j,j+1)
          call rotation_times_general(trp_rot(rot),pl(:,1:j+1),j,j+1)
       end do
       pl(ml,:)=pq(ml,ml)*pl(ml,:)
    end do
    call pop_id(error)
  end subroutine f_z_general_bv

end module mod_general_bv
