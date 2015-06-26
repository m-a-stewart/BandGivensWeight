module mod_solve
  use mod_prec
  use mod_error_id
  use mod_rotation
  use mod_orth_band_types
  use mod_band_types
  implicit none
  ! Various routines for fast backward and forward substitution for
  ! upper triangular matrices represented by BV or UB decompositions.
  ! The backward substitution solves $Rx=b$ while the forward
  ! substitution solve $x^T R = b^T$.
  private

  public :: back_solve_ub, d_back_solve_ub, z_back_solve_ub, &
          d_v_back_solve_ub, z_v_back_solve_ub, &
          f_back_solve_ub, f_d_back_solve_ub, f_z_back_solve_ub, &
          f_d_v_back_solve_ub, f_z_v_back_solve_ub, &
          d_solve_ub, d_v_solve_ub, z_solve_ub, z_v_solve_ub

  public :: forward_solve_bv, d_forward_solve_bv, z_forward_solve_bv, &
          d_v_forward_solve_bv, z_v_forward_solve_bv, &
          f_forward_solve_bv, f_d_forward_solve_bv, f_z_forward_solve_bv, &
          f_d_v_forward_solve_bv, f_z_v_forward_solve_bv, &
          d_solve_bv, d_v_solve_bv, z_solve_bv, z_v_solve_bv, solve

  interface solve
     module procedure d_solve_ub, d_v_solve_ub, z_solve_ub, z_v_solve_ub, &
          d_solve_bv, d_v_solve_bv, z_solve_bv, z_v_solve_bv
  end interface solve

  interface back_solve_ub
     module procedure d_back_solve_ub, z_back_solve_ub, &
          d_v_back_solve_ub, z_v_back_solve_ub
  end interface back_solve_ub

  interface f_back_solve_ub
     module procedure f_d_back_solve_ub, f_z_back_solve_ub, &
          f_d_v_back_solve_ub, f_z_v_back_solve_ub
  end interface f_back_solve_ub

  interface forward_solve_bv
     module procedure d_forward_solve_bv, z_forward_solve_bv, &
          d_v_forward_solve_bv, z_v_forward_solve_bv
  end interface forward_solve_bv

  interface f_forward_solve_bv
     module procedure f_d_forward_solve_bv, f_z_forward_solve_bv, &
          f_d_v_forward_solve_bv, f_z_v_forward_solve_bv
  end interface f_forward_solve_bv

contains

  function d_solve_ub(ub,c,error) result(x)
    real(kind=dp), dimension(:,:), allocatable :: x
    real(kind=dp), dimension(:,:), intent(in) :: c
    type(d_ub), intent(in) :: ub
    type(error_info), intent(inout), optional :: error

    real(kind=dp), dimension(:,:), allocatable :: c1
    type(routine_info), parameter :: info=info_d_solve_ub

    if (failure(error)) return
    call push_id(info,error)

    c1=c
    allocate(x(size(c,1),size(c,2)))
    call d_back_solve_ub(ub,x,c1,error)
    deallocate(c1)
    
    call pop_id(error)

  end function d_solve_ub

  
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
    type(d_ub), intent(in) :: ub
    real(kind=dp), dimension(:,:), intent(inout) :: c
    real(kind=dp), dimension(:,:), intent(out) :: x
    type(error_info), intent(inout), optional :: error
    type(routine_info), parameter :: info=info_d_back_solve_ub

    integer(kind=int32) :: n

    if (failure(error)) return
    call push_id(info, error)

    n=size(c,1)
    if (get_n(ub) /= n) then
       call set_error(1, info, error); return
    end if
    if (ub%lbw /= 0) then
       call set_error(2, info, error); return
    end if
    if (n < 1) then
       call set_error(3, info, error); return
    end if
    if (size(x,1)/=n .or. size(x,2) /= size(c,2)) then
       call set_error(4, info, error); return
    end if
    call f_d_back_solve_ub(ub%bc, n, ub%ubw, get_lbwmax(ub), get_ubwmax(ub), &
         ub%numrotsu, ub%jsu, ub%csu, ub%ssu, x, c, size(c,2))
    call pop_id(error)
  end subroutine d_back_solve_ub

  function d_v_solve_ub(ub,c,error) result(x)
    real(kind=dp), dimension(:), allocatable :: x
    real(kind=dp), dimension(:), intent(in) :: c
    type(d_ub), intent(in) :: ub
    type(error_info), intent(inout), optional :: error

    real(kind=dp), dimension(:), allocatable :: c1
    type(routine_info), parameter :: info=info_d_v_solve_ub

    if (failure(error)) return
    call push_id(info,error)

    c1=c
    allocate(x(size(c)))
    call d_v_back_solve_ub(ub,x,c1,error)
    deallocate(c1)
    
    call pop_id(error)

  end function d_v_solve_ub

  subroutine d_v_back_solve_ub(ub,x,c,error)
    type(d_ub), intent(in) :: ub
    real(kind=dp), dimension(:), intent(inout) :: c
    real(kind=dp), dimension(:), intent(out) :: x
    type(error_info), intent(inout), optional :: error
    type(routine_info), parameter :: info=info_d_v_back_solve_ub

    integer(kind=int32) :: n

    if (failure(error)) return
    call push_id(info, error)

    n=size(c)
    if (get_n(ub) /= n) then
       call set_error(1, info, error); return
    end if
    if (ub%lbw /= 0) then
       call set_error(2, info, error); return
    end if
    if (n < 1) then
       call set_error(3, info, error); return
    end if
    if (size(x)/=n) then
       call set_error(4, info, error); return
    end if
    call f_d_v_back_solve_ub(ub%bc, n, ub%ubw, get_lbwmax(ub), get_ubwmax(ub), &
         ub%numrotsu, ub%jsu, ub%csu, ub%ssu, x, c)
    call pop_id(error)
  end subroutine d_v_back_solve_ub


  subroutine f_d_back_solve_ub(b_ub, n, ubw_ub, lbwmax_ub, ubwmax_ub, numrotsu, jsu, &
       csu, ssu, x, c, nc)
    real(kind=dp), dimension(lbwmax_ub+ubwmax_ub+1,n), intent(in) :: b_ub
    integer(kind=int32), dimension(n), intent(in) :: numrotsu
    integer(kind=int32), dimension(ubwmax_ub,n), intent(in) :: jsu
    real(kind=dp), dimension(ubwmax_ub,n), intent(in) :: csu, ssu
    integer(kind=int32), intent(in) :: n, lbwmax_ub, ubwmax_ub, ubw_ub, nc
    real(kind=dp), dimension(n,nc), intent(inout) :: c
    real(kind=dp), dimension(n,nc), intent(out) :: x

    integer(kind=int32) :: j, k, l
    type(d_rotation) :: rot

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

  subroutine f_d_v_back_solve_ub(b_ub, n, ubw_ub, lbwmax_ub, ubwmax_ub, numrotsu, jsu, &
       csu, ssu, x, c)
    real(kind=dp), dimension(lbwmax_ub+ubwmax_ub+1,n), intent(in) :: b_ub
    integer(kind=int32), dimension(n), intent(in) :: numrotsu
    integer(kind=int32), dimension(ubwmax_ub,n), intent(in) :: jsu
    real(kind=dp), dimension(ubwmax_ub,n), intent(in) :: csu, ssu
    integer(kind=int32), intent(in) :: n, lbwmax_ub, ubwmax_ub, ubw_ub
    real(kind=dp), dimension(n), intent(inout) :: c
    real(kind=dp), dimension(n), intent(out) :: x

    integer(kind=int32) :: j, k, l
    type(d_rotation) :: rot

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

  function z_solve_ub(ub,c,error) result(x)
    complex(kind=dp), dimension(:,:), allocatable :: x
    complex(kind=dp), dimension(:,:), intent(in) :: c
    type(z_ub), intent(in) :: ub
    type(error_info), intent(inout), optional :: error

    complex(kind=dp), dimension(:,:), allocatable :: c1
    type(routine_info), parameter :: info=info_z_solve_ub

    if (failure(error)) return
    call push_id(info,error)

    c1=c
    allocate(x(size(c,1),size(c,2)))
    call z_back_solve_ub(ub,x,c1,error)
    deallocate(c1)
    
    call pop_id(error)

  end function z_solve_ub

  subroutine z_back_solve_ub(ub,x,c,error)
    type(z_ub), intent(in) :: ub
    complex(kind=dp), dimension(:,:), intent(inout) :: c
    complex(kind=dp), dimension(:,:), intent(out) :: x
    type(error_info), intent(inout), optional :: error
    type(routine_info), parameter :: info=info_z_back_solve_ub
    integer(kind=int32) :: n

    if (failure(error)) return
    call push_id(info, error)
    n=size(c,1)
    if (get_n(ub) /= n) then
       call set_error(1, info, error); return
    end if
    if (ub%lbw /= 0) then
       call set_error(2, info, error); return
    end if
    if (n < 1) then
       call set_error(3, info, error); return
    end if
    if (size(x,1)/=n .or. size(x,2) /= size(c,2)) then
       call set_error(4, info, error); return
    end if
    call f_z_back_solve_ub(ub%bc, n, ub%ubw, get_lbwmax(ub), get_ubwmax(ub), &
         ub%numrotsu, ub%jsu, ub%csu, ub%ssu, x, c, size(c,2))
    call pop_id(error)
  end subroutine z_back_solve_ub

  function z_v_solve_ub(ub,c,error) result(x)
    complex(kind=dp), dimension(:), allocatable :: x
    complex(kind=dp), dimension(:), intent(in) :: c
    type(z_ub), intent(in) :: ub
    type(error_info), intent(inout), optional :: error

    complex(kind=dp), dimension(:), allocatable :: c1
    type(routine_info), parameter :: info=info_z_v_solve_ub

    if (failure(error)) return
    call push_id(info,error)

    c1=c
    allocate(x(size(c)))
    call z_v_back_solve_ub(ub,x,c1,error)
    deallocate(c1)
    
    call pop_id(error)

  end function z_v_solve_ub

  subroutine z_v_back_solve_ub(ub,x,c,error)
    type(z_ub), intent(in) :: ub
    complex(kind=dp), dimension(:), intent(inout) :: c
    complex(kind=dp), dimension(:), intent(out) :: x
    type(error_info), intent(inout), optional :: error
    type(routine_info), parameter :: info=info_z_v_back_solve_ub

    integer(kind=int32) :: n

    if (failure(error)) return
    call push_id(info, error)
    n=size(c)
    if (get_n(ub) /= n) then
       call set_error(1, info, error); return
    end if
    if (ub%lbw /= 0) then
       call set_error(2, info, error); return
    end if
    if (n < 1) then
       call set_error(3, info, error); return
    end if
    if (size(x)/=n) then
       call set_error(4, info, error); return
    end if
    call f_z_v_back_solve_ub(ub%bc, n, ub%ubw, get_lbwmax(ub), get_ubwmax(ub), &
         ub%numrotsu, ub%jsu, ub%csu, ub%ssu, x, c)
    call pop_id(error)
  end subroutine z_v_back_solve_ub


  subroutine f_z_back_solve_ub(b_ub, n, ubw_ub, lbwmax_ub, ubwmax_ub, numrotsu, jsu, &
       csu, ssu, x, c, nc)
    complex(kind=dp), dimension(lbwmax_ub+ubwmax_ub+1,n), intent(in) :: b_ub
    integer(kind=int32), dimension(n), intent(in) :: numrotsu
    integer(kind=int32), dimension(ubwmax_ub,n), intent(in) :: jsu
    real(kind=dp), dimension(ubwmax_ub,n), intent(in) :: csu
    complex(kind=dp), dimension(ubwmax_ub,n), intent(in) :: ssu
    integer(kind=int32), intent(in) :: n, lbwmax_ub, ubwmax_ub, ubw_ub, nc
    complex(kind=dp), dimension(n,nc), intent(inout) :: c
    complex(kind=dp), dimension(n,nc), intent(out) :: x

    integer(kind=int32) :: j, k, l
    type(z_rotation) :: rot

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
  end subroutine f_z_back_solve_ub

  subroutine f_z_v_back_solve_ub(b_ub, n, ubw_ub, lbwmax_ub, ubwmax_ub, numrotsu, jsu, &
       csu, ssu, x, c)
    complex(kind=dp), dimension(lbwmax_ub+ubwmax_ub+1,n), intent(in) :: b_ub
    integer(kind=int32), dimension(n), intent(in) :: numrotsu
    integer(kind=int32), dimension(ubwmax_ub,n), intent(in) :: jsu
    real(kind=dp), dimension(ubwmax_ub,n), intent(in) :: csu
    complex(kind=dp), dimension(ubwmax_ub,n), intent(in) :: ssu
    integer(kind=int32), intent(in) :: n, lbwmax_ub, ubwmax_ub, ubw_ub
    complex(kind=dp), dimension(n), intent(inout) :: c
    complex(kind=dp), dimension(n), intent(out) :: x

    integer(kind=int32) :: j, k, l
    type(z_rotation) :: rot

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
  end subroutine f_z_v_back_solve_ub

  ! Forward solve:
  ! bv should represent an upper triangular matrix.  Solve x^T*bv=c^T.
  ! c is overwritten.

  function d_solve_bv(bv,c,error) result(x)
    real(kind=dp), dimension(:,:), allocatable :: x
    real(kind=dp), dimension(:,:), intent(in) :: c
    type(d_bv), intent(in) :: bv
    type(error_info), intent(inout), optional :: error

    real(kind=dp), dimension(:,:), allocatable :: c1
    type(routine_info), parameter :: info=info_d_solve_bv

    if (failure(error)) return
    call push_id(info,error)

    c1=c
    allocate(x(size(c,1),size(c,2)))
    call d_forward_solve_bv(x,bv,c1,error)
    deallocate(c1)
    
    call pop_id(error)

  end function d_solve_bv

  ! Errors
  ! 0: no error
  ! 1: ub%n /= size(c,1)
  ! 2: bv%lbw /= 0
  ! 3: n < 1
  ! 4: size(x,1) /= n or size(x,2) /= size(c,2)
  subroutine d_forward_solve_bv(x,bv,c,error)
    type(d_bv), intent(in) :: bv
    real(kind=dp), dimension(:,:), intent(inout) :: c
    real(kind=dp), dimension(:,:), intent(out) :: x
    type(error_info), intent(inout), optional :: error
    type(routine_info), parameter :: info=info_d_forward_solve_bv

    integer(kind=int32) :: n

    if (failure(error)) return
    call push_id(info, error)

    n=size(c,2)
    if (get_n(bv) /= n) then
       call set_error(1, info, error); return
    end if
    if (bv%lbw /= 0) then
       call set_error(2, info, error); return
    end if
    if (n < 1) then
       call set_error(3, info, error); return
    end if
    if (size(x,2)/=n .or. size(x,1) /= size(c,1)) then
       call set_error(4, info, error); return
    end if
    call f_d_forward_solve_bv(x, bv%br, n, bv%ubw, get_lbwmax(bv), get_ubwmax(bv), &
         bv%numrotsv, bv%ksv, bv%csv, bv%ssv, c, size(c,1))
    call pop_id(error)
  end subroutine d_forward_solve_bv

  subroutine f_d_forward_solve_bv(x, b_bv, n, ubw_bv, lbwmax_bv, ubwmax_bv, numrotsv, ksv, &
       csv, ssv, c, mc)
    real(kind=dp), dimension(n, lbwmax_bv+ubwmax_bv+1), intent(in) :: b_bv
    integer(kind=int32), dimension(n), intent(in) :: numrotsv
    integer(kind=int32), dimension(n,ubwmax_bv), intent(in) :: ksv
    real(kind=dp), dimension(n,ubwmax_bv), intent(in) :: csv, ssv
    integer(kind=int32), intent(in) :: n, lbwmax_bv, ubwmax_bv, ubw_bv, mc
    real(kind=dp), dimension(mc,n), intent(inout) :: c
    real(kind=dp), dimension(mc,n), intent(out) :: x

    integer(kind=int32) :: j, k, l
    type(d_rotation) :: rot

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
    ! diagonals of b are in column 1 of b_bv
    x(:,1)=c(:,1)/b_bv(1,1)
    do j=2,n
       l =  min(ubw_bv,n-j+1) ! number of superdiagonals in row j-1
       do k=1,mc
          c(k,j-1:j-1+l)=c(k,j-1:j-1+l) - b_bv(j-1,1:l+1) * x(k,j-1)
       end do
       ! Apply v_{n-j+1} to c
       do k=1,numrotsv(j-1)
          rot%cosine=csv(j-1,k); rot%sine=ssv(j-1,k)
          call general_times_rotation(c,trp_rot(rot),ksv(j-1,k), ksv(j-1,k)+1)
       end do
       x(:,j)=c(:,j)/b_bv(j,1)
    end do
  end subroutine f_d_forward_solve_bv

  function d_v_solve_bv(bv,c,error) result(x)
    real(kind=dp), dimension(:), allocatable :: x
    real(kind=dp), dimension(:), intent(in) :: c
    type(d_bv), intent(in) :: bv
    type(error_info), intent(inout), optional :: error

    real(kind=dp), dimension(:), allocatable :: c1
    type(routine_info), parameter :: info=info_d_v_solve_bv

    if (failure(error)) return
    call push_id(info,error)

    c1=c
    allocate(x(size(c)))
    call d_v_forward_solve_bv(x,bv,c1,error)
    deallocate(c1)
    
    call pop_id(error)
  end function d_v_solve_bv

  subroutine d_v_forward_solve_bv(x,bv,c,error)
    type(d_bv), intent(in) :: bv
    real(kind=dp), dimension(:), intent(inout) :: c
    real(kind=dp), dimension(:), intent(out) :: x
    type(error_info), intent(inout), optional :: error
    type(routine_info), parameter :: info=info_d_v_forward_solve_bv
    integer(kind=int32) :: n

    if (failure(error)) return
    call push_id(info, error)

    n=size(c)
    if (get_n(bv) /= n) then
       call set_error(1, info, error); return
    end if
    if (bv%lbw /= 0) then
       call set_error(2, info, error); return
    end if
    if (n < 1) then
       call set_error(3, info, error); return
    end if
    if (size(x)/=n) then
       call set_error(4, info, error); return
    end if
    call f_d_v_forward_solve_bv(x, bv%br, n, bv%ubw, get_lbwmax(bv), get_ubwmax(bv), &
         bv%numrotsv, bv%ksv, bv%csv, bv%ssv, c)
    call pop_id(error)
  end subroutine d_v_forward_solve_bv

  subroutine f_d_v_forward_solve_bv(x, b_bv, n, ubw_bv, lbwmax_bv, ubwmax_bv, numrotsv, ksv, &
       csv, ssv, c)
    real(kind=dp), dimension(n, lbwmax_bv+ubwmax_bv+1), intent(in) :: b_bv
    integer(kind=int32), dimension(n), intent(in) :: numrotsv
    integer(kind=int32), dimension(n,ubwmax_bv), intent(in) :: ksv
    real(kind=dp), dimension(n,ubwmax_bv), intent(in) :: csv, ssv
    integer(kind=int32), intent(in) :: n, lbwmax_bv, ubwmax_bv, ubw_bv
    real(kind=dp), dimension(n), intent(inout) :: c
    real(kind=dp), dimension(n), intent(out) :: x

    integer(kind=int32) :: j, k, l
    type(d_rotation) :: rot

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
    ! diagonals of b are in column 1 of b_bv
    x(1)=c(1)/b_bv(1,1)
    do j=2,n
       l =  min(ubw_bv,n-j+1) ! number of superdiagonals in row j-1
       c(j-1:j-1+l)=c(j-1:j-1+l) - b_bv(j-1,1:l+1) * x(j-1)
       ! Apply v_{n-j+1} to c
       do k=1,numrotsv(j-1)
          rot%cosine=csv(j-1,k); rot%sine=ssv(j-1,k)
          call general_times_rotation(c,trp_rot(rot),ksv(j-1,k), ksv(j-1,k)+1)
       end do
       x(j)=c(j)/b_bv(j,1)
    end do
  end subroutine f_d_v_forward_solve_bv

  ! Complex forward

  function z_solve_bv(bv,c,error) result(x)
    complex(kind=dp), dimension(:,:), allocatable :: x
    complex(kind=dp), dimension(:,:), intent(in) :: c
    type(z_bv), intent(in) :: bv
    type(error_info), intent(inout), optional :: error

    complex(kind=dp), dimension(:,:), allocatable :: c1
    type(routine_info), parameter :: info=info_z_solve_bv

    if (failure(error)) return
    call push_id(info,error)

    c1=c
    allocate(x(size(c,1),size(c,2)))
    call z_forward_solve_bv(x,bv,c1,error)
    deallocate(c1)
    
    call pop_id(error)

  end function z_solve_bv


  subroutine z_forward_solve_bv(x,bv,c,error)
    type(z_bv), intent(in) :: bv
    complex(kind=dp), dimension(:,:), intent(inout) :: c
    complex(kind=dp), dimension(:,:), intent(out) :: x
    type(error_info), intent(inout), optional :: error
    type(routine_info), parameter :: info=info_z_forward_solve_bv

    integer(kind=int32) :: n
    n=size(c,2)
    if (failure(error)) return
    call push_id(info, error)
    if (get_n(bv) /= n) then
       call set_error(1, info, error); return
    end if
    if (bv%lbw /= 0) then
       call set_error(2, info, error); return
    end if
    if (n < 1) then
       call set_error(3, info, error); return
    end if
    if (size(x,2)/=n .or. size(x,1) /= size(c,1)) then
       call set_error(4, info, error); return
    end if
    call f_z_forward_solve_bv(x, bv%br, n, bv%ubw, get_lbwmax(bv), get_ubwmax(bv), &
         bv%numrotsv, bv%ksv, bv%csv, bv%ssv, c, size(c,1))
    call pop_id(error)
  end subroutine z_forward_solve_bv

  subroutine f_z_forward_solve_bv(x, b_bv, n, ubw_bv, lbwmax_bv, ubwmax_bv, numrotsv, ksv, &
       csv, ssv, c, mc)
    complex(kind=dp), dimension(n, lbwmax_bv+ubwmax_bv+1), intent(in) :: b_bv
    integer(kind=int32), dimension(n), intent(in) :: numrotsv
    integer(kind=int32), dimension(n,ubwmax_bv), intent(in) :: ksv
    real(kind=dp), dimension(n,ubwmax_bv), intent(in) :: csv
    complex(kind=dp), dimension(n,ubwmax_bv), intent(in) :: ssv
    integer(kind=int32), intent(in) :: n, lbwmax_bv, ubwmax_bv, ubw_bv, mc
    complex(kind=dp), dimension(mc,n), intent(inout) :: c
    complex(kind=dp), dimension(mc,n), intent(out) :: x

    integer(kind=int32) :: j, k, l
    type(z_rotation) :: rot

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
    ! diagonals of b are in column 1 of b_bv
    x(:,1)=c(:,1)/b_bv(1,1)
    do j=2,n
       l =  min(ubw_bv,n-j+1) ! number of superdiagonals in row j-1
       do k=1,mc
          c(k,j-1:j-1+l)=c(k,j-1:j-1+l) - b_bv(j-1,1:l+1) * x(k,j-1)
       end do
       ! Apply v_{n-j+1} to c
       do k=1,numrotsv(j-1)
          rot%cosine=csv(j-1,k); rot%sine=ssv(j-1,k)
          call general_times_rotation(c,trp_rot(rot),ksv(j-1,k), ksv(j-1,k)+1)
       end do
       x(:,j)=c(:,j)/b_bv(j,1)
    end do
  end subroutine f_z_forward_solve_bv

  function z_v_solve_bv(bv,c,error) result(x)
    complex(kind=dp), dimension(:), allocatable :: x
    complex(kind=dp), dimension(:), intent(in) :: c
    type(z_bv), intent(in) :: bv
    type(error_info), intent(inout), optional :: error

    complex(kind=dp), dimension(:), allocatable :: c1
    type(routine_info), parameter :: info=info_z_v_solve_bv

    if (failure(error)) return
    call push_id(info,error)

    c1=c
    allocate(x(size(c)))
    call z_v_forward_solve_bv(x,bv,c1,error)
    deallocate(c1)
    
    call pop_id(error)
  end function z_v_solve_bv

  subroutine z_v_forward_solve_bv(x,bv,c,error)
    type(z_bv), intent(in) :: bv
    complex(kind=dp), dimension(:), intent(inout) :: c
    complex(kind=dp), dimension(:), intent(out) :: x
    type(error_info), intent(inout), optional :: error
    type(routine_info), parameter :: info=info_z_v_forward_solve_bv

    integer(kind=int32) :: n
    n=size(c)
    if (failure(error)) return
    call push_id(info, error)
    if (get_n(bv) /= n) then
       call set_error(1, info, error); return
    end if
    if (bv%lbw /= 0) then
       call set_error(2, info, error); return
    end if
    if (n < 1) then
       call set_error(3, info, error); return
    end if
    if (size(x)/=n) then
       call set_error(4, info, error); return
    end if
    call f_z_v_forward_solve_bv(x, bv%br, n, bv%ubw, get_lbwmax(bv), get_ubwmax(bv), &
         bv%numrotsv, bv%ksv, bv%csv, bv%ssv, c)
    call pop_id(error)
  end subroutine z_v_forward_solve_bv

  subroutine f_z_v_forward_solve_bv(x, b_bv, n, ubw_bv, lbwmax_bv, ubwmax_bv, numrotsv, ksv, &
       csv, ssv, c)
    complex(kind=dp), dimension(n, lbwmax_bv+ubwmax_bv+1), intent(in) :: b_bv
    integer(kind=int32), dimension(n), intent(in) :: numrotsv
    integer(kind=int32), dimension(n,ubwmax_bv), intent(in) :: ksv
    real(kind=dp), dimension(n,ubwmax_bv), intent(in) :: csv
    complex(kind=dp), dimension(n,ubwmax_bv), intent(in) :: ssv
    integer(kind=int32), intent(in) :: n, lbwmax_bv, ubwmax_bv, ubw_bv
    complex(kind=dp), dimension(n), intent(inout) :: c
    complex(kind=dp), dimension(n), intent(out) :: x

    integer(kind=int32) :: j, k, l
    type(z_rotation) :: rot

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
    ! diagonals of b are in column 1 of b_bv
    x(1)=c(1)/b_bv(1,1)
    do j=2,n
       l =  min(ubw_bv,n-j+1) ! number of superdiagonals in row j-1
       c(j-1:j-1+l)=c(j-1:j-1+l) - b_bv(j-1,1:l+1) * x(j-1)
       ! Apply v_{n-j+1} to c
       do k=1,numrotsv(j-1)
          rot%cosine=csv(j-1,k); rot%sine=ssv(j-1,k)
          call general_times_rotation(c,trp_rot(rot),ksv(j-1,k), ksv(j-1,k)+1)
       end do
       x(j)=c(j)/b_bv(j,1)
    end do
  end subroutine f_z_v_forward_solve_bv


end module mod_solve
