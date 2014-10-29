module mod_solve
  use mod_prec
  use mod_error_id
  use mod_rotation
  use mod_nested_types
  use mod_band_types
  implicit none

  private

  public :: back_solve_ub, d_back_solve_ub, c_back_solve_ub, &
          d_v_back_solve_ub, c_v_back_solve_ub, &
          f_back_solve_ub, f_d_back_solve_ub, f_c_back_solve_ub, &
          f_d_v_back_solve_ub, f_c_v_back_solve_ub

  public :: forward_solve_bv, d_forward_solve_bv, c_forward_solve_bv, &
          d_v_forward_solve_bv, c_v_forward_solve_bv, &
          f_forward_solve_bv, f_d_forward_solve_bv, f_c_forward_solve_bv, &
          f_d_v_forward_solve_bv, f_c_v_forward_solve_bv

  public :: info_d_back_solve_ub, info_c_back_solve_ub, &
       info_d_v_back_solve_ub, info_c_v_back_solve_ub, &
       info_d_forward_solve_bv, info_c_forward_solve_bv, &
       info_d_v_forward_solve_bv, info_c_v_forward_solve_bv

  interface back_solve_ub
     module procedure d_back_solve_ub, c_back_solve_ub, &
          d_v_back_solve_ub, c_v_back_solve_ub
  end interface back_solve_ub

  interface f_back_solve_ub
     module procedure f_d_back_solve_ub, f_c_back_solve_ub, &
          f_d_v_back_solve_ub, f_c_v_back_solve_ub
  end interface f_back_solve_ub

  interface forward_solve_bv
     module procedure d_forward_solve_bv, c_forward_solve_bv, &
          d_v_forward_solve_bv, c_v_forward_solve_bv
  end interface forward_solve_bv

  interface f_forward_solve_bv
     module procedure f_d_forward_solve_bv, f_c_forward_solve_bv, &
          f_d_v_forward_solve_bv, f_c_v_forward_solve_bv
  end interface f_forward_solve_bv


  type(routine_info), parameter :: info_d_back_solve_ub= &
       routine_info(id_d_back_solve_ub, &
       'd_back_solve_ub', &
       [ character(len=error_message_length) :: 'ub%n /= size(c,1)', 'ub%lbw /= 0', 'n < 1', &
       'x is not the same size as c.' ])

  type(routine_info), parameter :: info_c_back_solve_ub= &
       routine_info(id_c_back_solve_ub, &
       'd_back_solve_ub', &
       [ character(len=error_message_length) :: 'ub%n /= size(c,1)', 'ub%lbw /= 0', 'n < 1', &
       'x is not the same size as c.' ])

  type(routine_info), parameter :: info_d_v_back_solve_ub= &
       routine_info(id_d_v_back_solve_ub, &
       'd_v_back_solve_ub', &
       [ character(len=error_message_length) :: 'ub%n /= size(c)', 'ub%lbw /= 0', 'n < 1', &
       'x is not the same size as c.' ])

  type(routine_info), parameter :: info_c_v_back_solve_ub= &
       routine_info(id_c_v_back_solve_ub, &
       'd_back_solve_ub', &
       [ character(len=error_message_length) :: 'ub%n /= size(c)', 'ub%lbw /= 0', 'n < 1', &
       'x is not the same size as c.' ])


  type(routine_info), parameter :: info_d_forward_solve_bv= &
       routine_info(id_d_forward_solve_bv, &
       'd_forward_solve_bv', &
       [ character(len=error_message_length) :: 'bv%n /= size(c,2)', 'bv%lbw /= 0', 'n < 1', &
       'x is not the same size as c.' ])

  type(routine_info), parameter :: info_c_forward_solve_bv= &
       routine_info(id_c_forward_solve_bv, &
       'd_forward_solve_bv', &
       [ character(len=error_message_length) :: 'bv%n /= size(c,2)', 'bv%lbw /= 0', 'n < 1', &
       'x is not the same size as c.' ])

  type(routine_info), parameter :: info_d_v_forward_solve_bv= &
       routine_info(id_d_v_forward_solve_bv, &
       'd_v_forward_solve_bv', &
       [ character(len=error_message_length) :: 'bv%n /= size(c,2)', 'bv%lbw /= 0', 'n < 1', &
       'x is not the same size as c.' ])

  type(routine_info), parameter :: info_c_v_forward_solve_bv= &
       routine_info(id_c_v_forward_solve_bv, &
       'd_forward_solve_bv', &
       [ character(len=error_message_length) :: 'bv%n /= size(c,1)', 'bv%lbw /= 0', 'n < 1', &
       'x is not the same size as c.' ])

contains

  ! Back solve routines.
  ! ub should represent an upper triangular matrix.  Solve ub*x=c.
  ! c is overwritten.

  ! Errors
  ! 0: no error
  ! 1: ub%n /= size(c,1)
  ! 2: bv%lbw /= 0
  ! 3: n < 1
  ! 4: size(x,1) /= n or size(x,2) /= size(c,2)
  subroutine d_back_solve_ub(ub,x,c,error)
    type(d_ub) :: ub
    real(kind=dp), dimension(:,:), intent(inout) :: c
    real(kind=dp), dimension(:,:), intent(out) :: x
    type(error_info), intent(out) :: error

    integer(kind=int32) :: n

    call clear_error(error)
    n=size(c,1)
    if (get_n(ub) /= n) then
       call set_error(error,1,id_d_back_solve_ub); return
    end if
    if (ub%lbw /= 0) then
       call set_error(error,2,id_d_back_solve_ub); return
    end if
    if (n < 1) then
       call set_error(error, 3, id_d_back_solve_ub); return
    end if
    if (size(x,1)/=n .or. size(x,2) /= size(c,2)) then
       call set_error(error, 4, id_d_back_solve_ub); return
    end if
    call f_d_back_solve_ub(ub%bc, n, ub%lbw, ub%ubw, get_lbwmax(ub), get_ubwmax(ub), &
         ub%numrotsu, ub%jsu, ub%csu, ub%ssu, x, c, size(c,2), error)
  end subroutine d_back_solve_ub

  subroutine d_v_back_solve_ub(ub,x,c,error)
    type(d_ub) :: ub
    real(kind=dp), dimension(:), intent(inout) :: c
    real(kind=dp), dimension(:), intent(out) :: x
    type(error_info), intent(out) :: error

    integer(kind=int32) :: n

    call clear_error(error)
    n=size(c)
    if (get_n(ub) /= n) then
       call set_error(error,1,id_d_v_back_solve_ub); return
    end if
    if (ub%lbw /= 0) then
       call set_error(error,2,id_d_v_back_solve_ub); return
    end if
    if (n < 1) then
       call set_error(error, 3, id_d_v_back_solve_ub); return
    end if
    if (size(x)/=n) then
       call set_error(error, 4, id_d_v_back_solve_ub); return
    end if
    call f_d_v_back_solve_ub(ub%bc, n, ub%lbw, ub%ubw, get_lbwmax(ub), get_ubwmax(ub), &
         ub%numrotsu, ub%jsu, ub%csu, ub%ssu, x, c, error)
  end subroutine d_v_back_solve_ub


  subroutine f_d_back_solve_ub(b_ub, n, lbw_ub, ubw_ub, lbwmax_ub, ubwmax_ub, numrotsu, jsu, &
       csu, ssu, x, c, nc, error)
    real(kind=dp), dimension(lbwmax_ub+ubwmax_ub+1,n), intent(in) :: b_ub
    integer(kind=int32), dimension(n), intent(in) :: numrotsu
    integer(kind=int32), dimension(ubwmax_ub,n), intent(in) :: jsu
    real(kind=dp), dimension(ubwmax_ub,n), intent(in) :: csu, ssu
    integer(kind=int32), intent(in) :: n, lbwmax_ub, ubwmax_ub, lbw_ub, ubw_ub, nc
    real(kind=dp), dimension(n,nc), intent(inout) :: c
    real(kind=dp), dimension(n,nc), intent(out) :: x
    type(error_info), intent(out) :: error

    integer(kind=int32) :: j, k, l
    type(d_rotation) :: rot

    call clear_error(error)
    if (n==1) then
       x=c/b_ub(1,1); return
    end if
    x=0.0_dp
    ! Apply all the u_k^T to c.
    do k=1,n-1
       do j=numrotsu(k),1,-1
          rot%cosine=csu(j,k); rot%sine=ssu(j,k)
          call rotation_times_general(trp_rot(rot), c, jsu(j,k), jsu(j,k)+1)
       end do
    end do
    ! diagonals of b are in row ubw_ub+1 of b_ub.
    x(n,:)=c(n,:)/b_ub(ubw_ub+1,n)
    do k=n-1,1,-1
       l=min(ubw_ub,k)       ! number of superdiagonals in column k+1
       do j=1,nc
          c(k+1-l:k+1,j)=c(k+1-l:k+1,j) - x(k+1,j) * b_ub(ubw_ub+1-l:ubw_ub+1 ,k+1)
       end do
       ! Apply u_k to c
       do j=1,numrotsu(k)
          rot%cosine=csu(j,k); rot%sine=ssu(j,k)
          call rotation_times_general(rot, c, jsu(j,k), jsu(j,k)+1)
       end do
       x(k,:)=c(k,:)/b_ub(ubw_ub+1,k)
    end do
  end subroutine f_d_back_solve_ub

  subroutine f_d_v_back_solve_ub(b_ub, n, lbw_ub, ubw_ub, lbwmax_ub, ubwmax_ub, numrotsu, jsu, &
       csu, ssu, x, c, error)
    real(kind=dp), dimension(lbwmax_ub+ubwmax_ub+1,n), intent(in) :: b_ub
    integer(kind=int32), dimension(n), intent(in) :: numrotsu
    integer(kind=int32), dimension(ubwmax_ub,n), intent(in) :: jsu
    real(kind=dp), dimension(ubwmax_ub,n), intent(in) :: csu, ssu
    integer(kind=int32), intent(in) :: n, lbwmax_ub, ubwmax_ub, lbw_ub, ubw_ub
    real(kind=dp), dimension(n), intent(inout) :: c
    real(kind=dp), dimension(n), intent(out) :: x
    type(error_info), intent(out) :: error

    integer(kind=int32) :: j, k, l
    type(d_rotation) :: rot

    call clear_error(error)
    if (n==1) then
       x=c/b_ub(1,1); return
    end if
    x=0.0_dp
    ! Apply all the u_k^T to c.
    do k=1,n-1
       do j=numrotsu(k),1,-1
          rot%cosine=csu(j,k); rot%sine=ssu(j,k)
          call rotation_times_general(trp_rot(rot), c, jsu(j,k), jsu(j,k)+1)
       end do
    end do
    ! diagonals of b are in row ubw_ub+1 of b_ub.
    x(n)=c(n)/b_ub(ubw_ub+1,n)
    do k=n-1,1,-1
       l=min(ubw_ub,k)       ! number of superdiagonals in column k+1
       c(k+1-l:k+1)=c(k+1-l:k+1) - x(k+1) * b_ub(ubw_ub+1-l:ubw_ub+1 ,k+1)
       ! Apply u_k to c
       do j=1,numrotsu(k)
          rot%cosine=csu(j,k); rot%sine=ssu(j,k)
          call rotation_times_general(rot, c, jsu(j,k), jsu(j,k)+1)
       end do
       x(k)=c(k)/b_ub(ubw_ub+1,k)
    end do
  end subroutine f_d_v_back_solve_ub



  subroutine c_back_solve_ub(ub,x,c,error)
    type(c_ub) :: ub
    complex(kind=dp), dimension(:,:), intent(inout) :: c
    complex(kind=dp), dimension(:,:), intent(out) :: x
    type(error_info), intent(out) :: error

    integer(kind=int32) :: n

    call clear_error(error)
    n=size(c,1)
    if (get_n(ub) /= n) then
       call set_error(error,1,id_d_back_solve_ub); return
    end if
    if (ub%lbw /= 0) then
       call set_error(error,2,id_d_back_solve_ub); return
    end if
    if (n < 1) then
       call set_error(error, 3, id_d_back_solve_ub); return
    end if
    if (size(x,1)/=n .or. size(x,2) /= size(c,2)) then
       call set_error(error, 4, id_d_back_solve_ub); return
    end if
    call f_c_back_solve_ub(ub%bc, n, ub%lbw, ub%ubw, get_lbwmax(ub), get_ubwmax(ub), &
         ub%numrotsu, ub%jsu, ub%csu, ub%ssu, x, c, size(c,2), error)
  end subroutine c_back_solve_ub

  subroutine c_v_back_solve_ub(ub,x,c,error)
    type(c_ub) :: ub
    complex(kind=dp), dimension(:), intent(inout) :: c
    complex(kind=dp), dimension(:), intent(out) :: x
    type(error_info), intent(out) :: error

    integer(kind=int32) :: n

    call clear_error(error)
    n=size(c)
    if (get_n(ub) /= n) then
       call set_error(error,1,id_d_v_back_solve_ub); return
    end if
    if (ub%lbw /= 0) then
       call set_error(error,2,id_d_v_back_solve_ub); return
    end if
    if (n < 1) then
       call set_error(error, 3, id_d_v_back_solve_ub); return
    end if
    if (size(x)/=n) then
       call set_error(error, 4, id_d_v_back_solve_ub); return
    end if
    call f_c_v_back_solve_ub(ub%bc, n, ub%lbw, ub%ubw, get_lbwmax(ub), get_ubwmax(ub), &
         ub%numrotsu, ub%jsu, ub%csu, ub%ssu, x, c, error)
  end subroutine c_v_back_solve_ub


  subroutine f_c_back_solve_ub(b_ub, n, lbw_ub, ubw_ub, lbwmax_ub, ubwmax_ub, numrotsu, jsu, &
       csu, ssu, x, c, nc, error)
    complex(kind=dp), dimension(lbwmax_ub+ubwmax_ub+1,n), intent(in) :: b_ub
    integer(kind=int32), dimension(n), intent(in) :: numrotsu
    integer(kind=int32), dimension(ubwmax_ub,n), intent(in) :: jsu
    complex(kind=dp), dimension(ubwmax_ub,n), intent(in) :: csu, ssu
    integer(kind=int32), intent(in) :: n, lbwmax_ub, ubwmax_ub, lbw_ub, ubw_ub, nc
    complex(kind=dp), dimension(n,nc), intent(inout) :: c
    complex(kind=dp), dimension(n,nc), intent(out) :: x
    type(error_info), intent(out) :: error

    integer(kind=int32) :: j, k, l
    type(c_rotation) :: rot

    call clear_error(error)
    if (n==1) then
       x=c/b_ub(1,1); return
    end if
    x=(0.0_dp,0.0_dp)
    ! Apply all the u_k^T to c.
    do k=1,n-1
       do j=numrotsu(k),1,-1
          rot%cosine=csu(j,k); rot%sine=ssu(j,k)
          call rotation_times_general(trp_rot(rot), c, jsu(j,k), jsu(j,k)+1)
       end do
    end do
    ! diagonals of b are in row ubw_ub+1 of b_ub.
    x(n,:)=c(n,:)/b_ub(ubw_ub+1,n)
    do k=n-1,1,-1
       l=min(ubw_ub,k)       ! number of superdiagonals in column k+1
       do j=1,nc
          c(k+1-l:k+1,j)=c(k+1-l:k+1,j) - x(k+1,j) * b_ub(ubw_ub+1-l:ubw_ub+1 ,k+1)
       end do
       ! Apply u_k to c
       do j=1,numrotsu(k)
          rot%cosine=csu(j,k); rot%sine=ssu(j,k)
          call rotation_times_general(rot, c, jsu(j,k), jsu(j,k)+1)
       end do
       x(k,:)=c(k,:)/b_ub(ubw_ub+1,k)
    end do
  end subroutine f_c_back_solve_ub

  subroutine f_c_v_back_solve_ub(b_ub, n, lbw_ub, ubw_ub, lbwmax_ub, ubwmax_ub, numrotsu, jsu, &
       csu, ssu, x, c, error)
    complex(kind=dp), dimension(lbwmax_ub+ubwmax_ub+1,n), intent(in) :: b_ub
    integer(kind=int32), dimension(n), intent(in) :: numrotsu
    integer(kind=int32), dimension(ubwmax_ub,n), intent(in) :: jsu
    complex(kind=dp), dimension(ubwmax_ub,n), intent(in) :: csu, ssu
    integer(kind=int32), intent(in) :: n, lbwmax_ub, ubwmax_ub, lbw_ub, ubw_ub
    complex(kind=dp), dimension(n), intent(inout) :: c
    complex(kind=dp), dimension(n), intent(out) :: x
    type(error_info), intent(out) :: error

    integer(kind=int32) :: j, k, l
    type(c_rotation) :: rot

    call clear_error(error)
    if (n==1) then
       x=c/b_ub(1,1); return
    end if
    x=0.0_dp
    ! Apply all the u_k^T to c.
    do k=1,n-1
       do j=numrotsu(k),1,-1
          rot%cosine=csu(j,k); rot%sine=ssu(j,k)
          call rotation_times_general(trp_rot(rot), c, jsu(j,k), jsu(j,k)+1)
       end do
    end do
    ! diagonals of b are in row ubw_ub+1 of b_ub.
    x(n)=c(n)/b_ub(ubw_ub+1,n)
    do k=n-1,1,-1
       l=min(ubw_ub,k)       ! number of superdiagonals in column k+1
       c(k+1-l:k+1)=c(k+1-l:k+1) - x(k+1) * b_ub(ubw_ub+1-l:ubw_ub+1 ,k+1)
       ! Apply u_k to c
       do j=1,numrotsu(k)
          rot%cosine=csu(j,k); rot%sine=ssu(j,k)
          call rotation_times_general(rot, c, jsu(j,k), jsu(j,k)+1)
       end do
       x(k)=c(k)/b_ub(ubw_ub+1,k)
    end do
  end subroutine f_c_v_back_solve_ub


  ! Forward solve:
  ! bv should represent an upper triangular matrix.  Solve x*ub=c.
  ! c is overwritten.

  ! Errors
  ! 0: no error
  ! 1: ub%n /= size(c,1)
  ! 2: bv%lbw /= 0
  ! 3: n < 1
  ! 4: size(x,1) /= n or size(x,2) /= size(c,2)
  subroutine d_forward_solve_bv(x,bv,c,error)
    type(d_bv) :: bv
    real(kind=dp), dimension(:,:), intent(inout) :: c
    real(kind=dp), dimension(:,:), intent(out) :: x
    type(error_info), intent(out) :: error

    integer(kind=int32) :: n

    call clear_error(error)
    n=size(c,2)
    if (get_n(bv) /= n) then
       call set_error(error,1,id_d_forward_solve_bv); return
    end if
    if (bv%lbw /= 0) then
       call set_error(error,2,id_d_forward_solve_bv); return
    end if
    if (n < 1) then
       call set_error(error, 3,id_d_forward_solve_bv); return
    end if
    if (size(x,2)/=n .or. size(x,1) /= size(c,1)) then
       call set_error(error, 4,id_d_forward_solve_bv); return
    end if
    call f_d_forward_solve_bv(x, bv%br, n, bv%lbw, bv%ubw, get_lbwmax(bv), get_ubwmax(bv), &
         bv%numrotsv, bv%ksv, bv%csv, bv%ssv, c, size(c,1), error)

  end subroutine d_forward_solve_bv

  subroutine f_d_forward_solve_bv(x, b_bv, n, lbw_bv, ubw_bv, lbwmax_bv, ubwmax_bv, numrotsv, ksv, &
       csv, ssv, c, mc, error)
    real(kind=dp), dimension(n, lbwmax_bv+ubwmax_bv+1), intent(in) :: b_bv
    integer(kind=int32), dimension(n), intent(in) :: numrotsv
    integer(kind=int32), dimension(n,ubwmax_bv), intent(in) :: ksv
    real(kind=dp), dimension(n,ubwmax_bv), intent(in) :: csv, ssv
    integer(kind=int32), intent(in) :: n, lbwmax_bv, ubwmax_bv, lbw_bv, ubw_bv, mc
    real(kind=dp), dimension(mc,n), intent(inout) :: c
    real(kind=dp), dimension(mc,n), intent(out) :: x
    type(error_info), intent(out) :: error

    integer(kind=int32) :: j, k, l
    type(d_rotation) :: rot

    call clear_error(error)
    if (n==1) then
       x=c/b_bv(1,1); return
    end if
    x=0.0_dp
    ! Apply all the v_k to c
    do j=1,n-1
       do k=numrotsv(n-j),1,-1
          rot%cosine=csv(n-j,k); rot%sine=ssv(n-j,k)
          call general_times_rotation(c,rot,ksv(n-j,k),ksv(n-j,k)+1)
       end do
    end do
    ! diagonals of b are in column lbw_bv+1 of b_bv
    x(:,1)=c(:,1)/b_bv(1,lbw_bv+1)
    do j=2,n
       l =  min(ubw_bv,n-j+1) ! number of superdiagonals in row j-1
       do k=1,mc
          c(k,j-1:j-1+l)=c(k,j-1:j-1+l) - b_bv(j-1,lbw_bv+1:lbw_bv+l+1) * x(k,j-1)
       end do
       ! Apply v_{n-j+1} to c
       do k=1,numrotsv(j-1)
          rot%cosine=csv(j-1,k); rot%sine=ssv(j-1,k)
          call general_times_rotation(c,trp_rot(rot),ksv(j-1,k), ksv(j-1,k)+1)
       end do
       x(:,j)=c(:,j)/b_bv(j,lbw_bv+1)
    end do
  end subroutine f_d_forward_solve_bv

  subroutine d_v_forward_solve_bv(x,bv,c,error)
    type(d_bv) :: bv
    real(kind=dp), dimension(:), intent(inout) :: c
    real(kind=dp), dimension(:), intent(out) :: x
    type(error_info), intent(out) :: error

    integer(kind=int32) :: n

    call clear_error(error)
    n=size(c)
    if (get_n(bv) /= n) then
       call set_error(error,1,id_d_forward_solve_bv); return
    end if
    if (bv%lbw /= 0) then
       call set_error(error,2,id_d_forward_solve_bv); return
    end if
    if (n < 1) then
       call set_error(error, 3,id_d_forward_solve_bv); return
    end if
    if (size(x)/=n) then
       call set_error(error, 4,id_d_forward_solve_bv); return
    end if
    call f_d_v_forward_solve_bv(x, bv%br, n, bv%lbw, bv%ubw, get_lbwmax(bv), get_ubwmax(bv), &
         bv%numrotsv, bv%ksv, bv%csv, bv%ssv, c, error)
  end subroutine d_v_forward_solve_bv

  subroutine f_d_v_forward_solve_bv(x, b_bv, n, lbw_bv, ubw_bv, lbwmax_bv, ubwmax_bv, numrotsv, ksv, &
       csv, ssv, c, error)
    real(kind=dp), dimension(n, lbwmax_bv+ubwmax_bv+1), intent(in) :: b_bv
    integer(kind=int32), dimension(n), intent(in) :: numrotsv
    integer(kind=int32), dimension(n,ubwmax_bv), intent(in) :: ksv
    real(kind=dp), dimension(n,ubwmax_bv), intent(in) :: csv, ssv
    integer(kind=int32), intent(in) :: n, lbwmax_bv, ubwmax_bv, lbw_bv, ubw_bv
    real(kind=dp), dimension(n), intent(inout) :: c
    real(kind=dp), dimension(n), intent(out) :: x
    type(error_info), intent(out) :: error

    integer(kind=int32) :: j, k, l
    type(d_rotation) :: rot

    call clear_error(error)
    if (n==1) then
       x=c/b_bv(1,1); return
    end if
    x=0.0_dp
    ! Apply all the v_k to c
    do j=1,n-1
       do k=numrotsv(n-j),1,-1
          rot%cosine=csv(n-j,k); rot%sine=ssv(n-j,k)
          call general_times_rotation(c,rot,ksv(n-j,k),ksv(n-j,k)+1)
       end do
    end do
    ! diagonals of b are in column lbw_bv+1 of b_bv
    x(1)=c(1)/b_bv(1,lbw_bv+1)
    do j=2,n
       l =  min(ubw_bv,n-j+1) ! number of superdiagonals in row j-1
       c(j-1:j-1+l)=c(j-1:j-1+l) - b_bv(j-1,lbw_bv+1:lbw_bv+l+1) * x(j-1)
       ! Apply v_{n-j+1} to c
       do k=1,numrotsv(j-1)
          rot%cosine=csv(j-1,k); rot%sine=ssv(j-1,k)
          call general_times_rotation(c,trp_rot(rot),ksv(j-1,k), ksv(j-1,k)+1)
       end do
       x(j)=c(j)/b_bv(j,lbw_bv+1)
    end do
  end subroutine f_d_v_forward_solve_bv

  ! Complex forward

  subroutine c_forward_solve_bv(x,bv,c,error)
    type(c_bv) :: bv
    complex(kind=dp), dimension(:,:), intent(inout) :: c
    complex(kind=dp), dimension(:,:), intent(out) :: x
    type(error_info), intent(out) :: error

    integer(kind=int32) :: n
    n=size(c,2)
    call clear_error(error)
    if (get_n(bv) /= n) then
       call set_error(error,1,id_c_forward_solve_bv); return
    end if
    if (bv%lbw /= 0) then
       call set_error(error,2,id_c_forward_solve_bv); return
    end if
    if (n < 1) then
       call set_error(error, 3,id_c_forward_solve_bv); return
    end if
    if (size(x,2)/=n .or. size(x,1) /= size(c,1)) then
       call set_error(error, 4,id_c_forward_solve_bv); return
    end if
    call f_c_forward_solve_bv(x, bv%br, n, bv%lbw, bv%ubw, get_lbwmax(bv), get_ubwmax(bv), &
         bv%numrotsv, bv%ksv, bv%csv, bv%ssv, c, size(c,1), error)
  end subroutine c_forward_solve_bv

  subroutine f_c_forward_solve_bv(x, b_bv, n, lbw_bv, ubw_bv, lbwmax_bv, ubwmax_bv, numrotsv, ksv, &
       csv, ssv, c, mc, error)
    complex(kind=dp), dimension(n, lbwmax_bv+ubwmax_bv+1), intent(in) :: b_bv
    integer(kind=int32), dimension(n), intent(in) :: numrotsv
    integer(kind=int32), dimension(n,ubwmax_bv), intent(in) :: ksv
    complex(kind=dp), dimension(n,ubwmax_bv), intent(in) :: csv, ssv
    integer(kind=int32), intent(in) :: n, lbwmax_bv, ubwmax_bv, lbw_bv, ubw_bv, mc
    complex(kind=dp), dimension(mc,n), intent(inout) :: c
    complex(kind=dp), dimension(mc,n), intent(out) :: x
    type(error_info), intent(out) :: error

    integer(kind=int32) :: j, k, l
    type(c_rotation) :: rot

    call clear_error(error)
    if (n==1) then
       x=c/b_bv(1,1); return
    end if
    x=(0.0_dp,0.0_dp)
    ! Apply all the v_k to c
    do j=1,n-1
       do k=numrotsv(n-j),1,-1
          rot%cosine=csv(n-j,k); rot%sine=ssv(n-j,k)
          call general_times_rotation(c,rot,ksv(n-j,k),ksv(n-j,k)+1)
       end do
    end do
    ! diagonals of b are in column lbw_bv+1 of b_bv
    x(:,1)=c(:,1)/b_bv(1,lbw_bv+1)
    do j=2,n
       l =  min(ubw_bv,n-j+1) ! number of superdiagonals in row j-1
       do k=1,mc
          c(k,j-1:j-1+l)=c(k,j-1:j-1+l) - b_bv(j-1,lbw_bv+1:lbw_bv+l+1) * x(k,j-1)
       end do
       ! Apply v_{n-j+1} to c
       do k=1,numrotsv(j-1)
          rot%cosine=csv(j-1,k); rot%sine=ssv(j-1,k)
          call general_times_rotation(c,trp_rot(rot),ksv(j-1,k), ksv(j-1,k)+1)
       end do
       x(:,j)=c(:,j)/b_bv(j,lbw_bv+1)
    end do
  end subroutine f_c_forward_solve_bv

  subroutine c_v_forward_solve_bv(x,bv,c,error)
    type(c_bv) :: bv
    complex(kind=dp), dimension(:), intent(inout) :: c
    complex(kind=dp), dimension(:), intent(out) :: x
    type(error_info), intent(out) :: error

    integer(kind=int32) :: n
    n=size(c)
    call clear_error(error)
    if (get_n(bv) /= n) then
       call set_error(error,1,id_c_forward_solve_bv); return
    end if
    if (bv%lbw /= 0) then
       call set_error(error,2,id_c_forward_solve_bv); return
    end if
    if (n < 1) then
       call set_error(error, 3,id_c_forward_solve_bv); return
    end if
    if (size(x)/=n) then
       call set_error(error, 4,id_c_forward_solve_bv); return
    end if
    call f_c_v_forward_solve_bv(x, bv%br, n, bv%lbw, bv%ubw, get_lbwmax(bv), get_ubwmax(bv), &
         bv%numrotsv, bv%ksv, bv%csv, bv%ssv, c, error)
  end subroutine c_v_forward_solve_bv

  subroutine f_c_v_forward_solve_bv(x, b_bv, n, lbw_bv, ubw_bv, lbwmax_bv, ubwmax_bv, numrotsv, ksv, &
       csv, ssv, c, error)
    complex(kind=dp), dimension(n, lbwmax_bv+ubwmax_bv+1), intent(in) :: b_bv
    integer(kind=int32), dimension(n), intent(in) :: numrotsv
    integer(kind=int32), dimension(n,ubwmax_bv), intent(in) :: ksv
    complex(kind=dp), dimension(n,ubwmax_bv), intent(in) :: csv, ssv
    integer(kind=int32), intent(in) :: n, lbwmax_bv, ubwmax_bv, lbw_bv, ubw_bv
    complex(kind=dp), dimension(n), intent(inout) :: c
    complex(kind=dp), dimension(n), intent(out) :: x
    type(error_info), intent(out) :: error

    integer(kind=int32) :: j, k, l
    type(c_rotation) :: rot

    call clear_error(error)
    if (n==1) then
       x=c/b_bv(1,1); return
    end if
    x=(0.0_dp,0.0_dp)
    ! Apply all the v_k to c
    do j=1,n-1
       do k=numrotsv(n-j),1,-1
          rot%cosine=csv(n-j,k); rot%sine=ssv(n-j,k)
          call general_times_rotation(c,rot,ksv(n-j,k),ksv(n-j,k)+1)
       end do
    end do
    ! diagonals of b are in column lbw_bv+1 of b_bv
    x(1)=c(1)/b_bv(1,lbw_bv+1)
    do j=2,n
       l =  min(ubw_bv,n-j+1) ! number of superdiagonals in row j-1
       c(j-1:j-1+l)=c(j-1:j-1+l) - b_bv(j-1,lbw_bv+1:lbw_bv+l+1) * x(j-1)
       ! Apply v_{n-j+1} to c
       do k=1,numrotsv(j-1)
          rot%cosine=csv(j-1,k); rot%sine=ssv(j-1,k)
          call general_times_rotation(c,trp_rot(rot),ksv(j-1,k), ksv(j-1,k)+1)
       end do
       x(j)=c(j)/b_bv(j,lbw_bv+1)
    end do
  end subroutine f_c_v_forward_solve_bv


end module mod_solve
