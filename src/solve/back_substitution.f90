module back_substitution
  use transforms
  use misc
  use types
  implicit none

  interface back_substitution_ub
     module procedure d_back_substitution_ub, c_back_substitution_ub, &
          d_v_back_substitution_ub, c_v_back_substitution_ub
  end interface back_substitution_ub

  interface f_back_substitution_ub
     module procedure f_d_back_substitution_ub, f_c_back_substitution_ub, &
          f_d_v_back_substitution_ub, f_c_v_back_substitution_ub
  end interface f_back_substitution_ub

  type(routine_info), parameter :: info_d_back_substitution_ub= &
       routine_info(id_d_back_substitution_ub, &
       'd_back_substitution_ub', &
       [ character(len=error_message_length) :: 'ub%n /= size(c,1)', 'bv%lbw /= 0', 'n < 1', &
       'x is not the same size as c.' ])

  type(routine_info), parameter :: info_c_back_substitution_ub= &
       routine_info(id_c_back_substitution_ub, &
       'd_back_substitution_ub', &
       [ character(len=error_message_length) :: 'ub%n /= size(c,1)', 'bv%lbw /= 0', 'n < 1', &
       'x is not the same size as c.' ])

  type(routine_info), parameter :: info_d_v_back_substitution_ub= &
       routine_info(id_d_v_back_substitution_ub, &
       'd_v_back_substitution_ub', &
       [ character(len=error_message_length) :: 'ub%n /= size(c,1)', 'bv%lbw /= 0', 'n < 1', &
       'x is not the same size as c.' ])

  type(routine_info), parameter :: info_c_v_back_substitution_ub= &
       routine_info(id_c_v_back_substitution_ub, &
       'd_back_substitution_ub', &
       [ character(len=error_message_length) :: 'ub%n /= size(c,1)', 'bv%lbw /= 0', 'n < 1', &
       'x is not the same size as c.' ])

contains

  !
  ! ub should represent an upper triangular matrix.  Solve ub*x=c.
  ! c is overwritten.

  ! Errors
  ! 0: no error
  ! 1: ub%n /= size(c,1)
  ! 2: bv%lbw /= 0
  ! 3: n < 1
  ! 4: size(x,1) /= n or size(x,2) /= size(c,2)
  subroutine d_back_substitution_ub(ub,x,c,error)
    type(d_ub) :: ub
    real(kind=dp), dimension(:,:), intent(inout) :: c
    real(kind=dp), dimension(:,:), intent(out) :: x
    type(error_info), intent(out) :: error

    integer(kind=int32) :: n
    n=size(c,1)
    if (get_n(ub) /= n) then
       call set_error(error,1,id_d_back_substitution_ub); return
    end if
    if (ub%lbw /= 0) then
       call set_error(error,2,id_d_back_substitution_ub); return
    end if
    if (n < 1) then
       call set_error(error, 3, id_d_back_substitution_ub); return
    end if
    if (size(x,1)/=n .or. size(x,2) /= size(c,2)) then
       call set_error(error, 4, id_d_back_substitution_ub); return
    end if
    call f_d_back_substitution_ub(ub%b, n, ub%lbw, ub%ubw, get_lbwmax(ub), get_ubwmax(ub), &
         ub%numrotsu, ub%jsu, ub%csu, ub%ssu, x, c, size(c,2), error)
  end subroutine d_back_substitution_ub

  subroutine d_v_back_substitution_ub(ub,x,c,error)
    type(d_ub) :: ub
    real(kind=dp), dimension(:), intent(inout) :: c
    real(kind=dp), dimension(:), intent(out) :: x
    type(error_info), intent(out) :: error

    integer(kind=int32) :: n
    n=size(c)
    if (get_n(ub) /= n) then
       call set_error(error,1,id_d_v_back_substitution_ub); return
    end if
    if (ub%lbw /= 0) then
       call set_error(error,2,id_d_v_back_substitution_ub); return
    end if
    if (n < 1) then
       call set_error(error, 3, id_d_v_back_substitution_ub); return
    end if
    if (size(x)/=n) then
       call set_error(error, 4, id_d_v_back_substitution_ub); return
    end if
    call f_d_v_back_substitution_ub(ub%b, n, ub%lbw, ub%ubw, get_lbwmax(ub), get_ubwmax(ub), &
         ub%numrotsu, ub%jsu, ub%csu, ub%ssu, x, c, error)
  end subroutine d_v_back_substitution_ub


  subroutine f_d_back_substitution_ub(b_ub, n, lbw_ub, ubw_ub, lbwmax_ub, ubwmax_ub, numrots_ub, js_ub, &
       cs_ub, ss_ub, x, c, nc, error)
    real(kind=dp), dimension(lbwmax_ub+ubwmax_ub+1,n), intent(in) :: b_ub
    integer(kind=int32), dimension(n), intent(in) :: numrots_ub
    integer(kind=int32), dimension(ubwmax_ub,n), intent(in) :: js_ub
    real(kind=dp), dimension(ubwmax_ub,n), intent(in) :: cs_ub, ss_ub
    integer(kind=int32), intent(in) :: n, lbwmax_ub, ubwmax_ub, lbw_ub, ubw_ub, nc
    real(kind=dp), dimension(n,nc), intent(inout) :: c
    real(kind=dp), dimension(n,nc), intent(out) :: x
    type(error_info), intent(out) :: error

    integer(kind=int32) :: j, k, l
    type(d_rotation) :: rot

    if (n==1) then
       x=c/b_ub(1,1); return
    end if
    x=0.0_dp
    ! Apply all the u_k^T to c.
    do k=1,n-1
       do j=numrots_ub(k),1,-1
          rot%cosine=cs_ub(j,k); rot%sine=ss_ub(j,k)
          call rotation_times_general(trp_rot(rot), c, js_ub(j,k), js_ub(j,k)+1)
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
       do j=1,numrots_ub(k)
          rot%cosine=cs_ub(j,k); rot%sine=ss_ub(j,k)
          call rotation_times_general(rot, c, js_ub(j,k), js_ub(j,k)+1)
       end do
       x(k,:)=c(k,:)/b_ub(ubw_ub+1,k)
    end do
  end subroutine f_d_back_substitution_ub

  subroutine f_d_v_back_substitution_ub(b_ub, n, lbw_ub, ubw_ub, lbwmax_ub, ubwmax_ub, numrots_ub, js_ub, &
       cs_ub, ss_ub, x, c, error)
    real(kind=dp), dimension(lbwmax_ub+ubwmax_ub+1,n), intent(in) :: b_ub
    integer(kind=int32), dimension(n), intent(in) :: numrots_ub
    integer(kind=int32), dimension(ubwmax_ub,n), intent(in) :: js_ub
    real(kind=dp), dimension(ubwmax_ub,n), intent(in) :: cs_ub, ss_ub
    integer(kind=int32), intent(in) :: n, lbwmax_ub, ubwmax_ub, lbw_ub, ubw_ub
    real(kind=dp), dimension(n), intent(inout) :: c
    real(kind=dp), dimension(n), intent(out) :: x
    type(error_info), intent(out) :: error

    integer(kind=int32) :: j, k, l
    type(d_rotation) :: rot

    if (n==1) then
       x=c/b_ub(1,1); return
    end if
    x=0.0_dp
    ! Apply all the u_k^T to c.
    do k=1,n-1
       do j=numrots_ub(k),1,-1
          rot%cosine=cs_ub(j,k); rot%sine=ss_ub(j,k)
          call rotation_times_general(trp_rot(rot), c, js_ub(j,k), js_ub(j,k)+1)
       end do
    end do
    ! diagonals of b are in row ubw_ub+1 of b_ub.
    x(n)=c(n)/b_ub(ubw_ub+1,n)
    do k=n-1,1,-1
       l=min(ubw_ub,k)       ! number of superdiagonals in column k+1
       c(k+1-l:k+1)=c(k+1-l:k+1) - x(k+1) * b_ub(ubw_ub+1-l:ubw_ub+1 ,k+1)
       ! Apply u_k to c
       do j=1,numrots_ub(k)
          rot%cosine=cs_ub(j,k); rot%sine=ss_ub(j,k)
          call rotation_times_general(rot, c, js_ub(j,k), js_ub(j,k)+1)
       end do
       x(k)=c(k)/b_ub(ubw_ub+1,k)
    end do
  end subroutine f_d_v_back_substitution_ub



  subroutine c_back_substitution_ub(ub,x,c,error)
    type(c_ub) :: ub
    complex(kind=dp), dimension(:,:), intent(inout) :: c
    complex(kind=dp), dimension(:,:), intent(out) :: x
    type(error_info), intent(out) :: error

    integer(kind=int32) :: n
    n=size(c,1)
    if (get_n(ub) /= n) then
       call set_error(error,1,id_d_back_substitution_ub); return
    end if
    if (ub%lbw /= 0) then
       call set_error(error,2,id_d_back_substitution_ub); return
    end if
    if (n < 1) then
       call set_error(error, 3, id_d_back_substitution_ub); return
    end if
    if (size(x,1)/=n .or. size(x,2) /= size(c,2)) then
       call set_error(error, 4, id_d_back_substitution_ub); return
    end if
    call f_c_back_substitution_ub(ub%b, n, ub%lbw, ub%ubw, get_lbwmax(ub), get_ubwmax(ub), &
         ub%numrotsu, ub%jsu, ub%csu, ub%ssu, x, c, size(c,2), error)
  end subroutine c_back_substitution_ub

  subroutine c_v_back_substitution_ub(ub,x,c,error)
    type(c_ub) :: ub
    complex(kind=dp), dimension(:), intent(inout) :: c
    complex(kind=dp), dimension(:), intent(out) :: x
    type(error_info), intent(out) :: error

    integer(kind=int32) :: n
    n=size(c)
    if (get_n(ub) /= n) then
       call set_error(error,1,id_d_v_back_substitution_ub); return
    end if
    if (ub%lbw /= 0) then
       call set_error(error,2,id_d_v_back_substitution_ub); return
    end if
    if (n < 1) then
       call set_error(error, 3, id_d_v_back_substitution_ub); return
    end if
    if (size(x)/=n) then
       call set_error(error, 4, id_d_v_back_substitution_ub); return
    end if
    call f_c_v_back_substitution_ub(ub%b, n, ub%lbw, ub%ubw, get_lbwmax(ub), get_ubwmax(ub), &
         ub%numrotsu, ub%jsu, ub%csu, ub%ssu, x, c, error)
  end subroutine c_v_back_substitution_ub


  subroutine f_c_back_substitution_ub(b_ub, n, lbw_ub, ubw_ub, lbwmax_ub, ubwmax_ub, numrots_ub, js_ub, &
       cs_ub, ss_ub, x, c, nc, error)
    complex(kind=dp), dimension(lbwmax_ub+ubwmax_ub+1,n), intent(in) :: b_ub
    integer(kind=int32), dimension(n), intent(in) :: numrots_ub
    integer(kind=int32), dimension(ubwmax_ub,n), intent(in) :: js_ub
    complex(kind=dp), dimension(ubwmax_ub,n), intent(in) :: cs_ub, ss_ub
    integer(kind=int32), intent(in) :: n, lbwmax_ub, ubwmax_ub, lbw_ub, ubw_ub, nc
    complex(kind=dp), dimension(n,nc), intent(inout) :: c
    complex(kind=dp), dimension(n,nc), intent(out) :: x
    type(error_info), intent(out) :: error

    integer(kind=int32) :: j, k, l
    type(c_rotation) :: rot

    if (n==1) then
       x=c/b_ub(1,1); return
    end if
    x=(0.0_dp,0.0_dp)
    ! Apply all the u_k^T to c.
    do k=1,n-1
       do j=numrots_ub(k),1,-1
          rot%cosine=cs_ub(j,k); rot%sine=ss_ub(j,k)
          call rotation_times_general(trp_rot(rot), c, js_ub(j,k), js_ub(j,k)+1)
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
       do j=1,numrots_ub(k)
          rot%cosine=cs_ub(j,k); rot%sine=ss_ub(j,k)
          call rotation_times_general(rot, c, js_ub(j,k), js_ub(j,k)+1)
       end do
       x(k,:)=c(k,:)/b_ub(ubw_ub+1,k)
    end do
  end subroutine f_c_back_substitution_ub

  subroutine f_c_v_back_substitution_ub(b_ub, n, lbw_ub, ubw_ub, lbwmax_ub, ubwmax_ub, numrots_ub, js_ub, &
       cs_ub, ss_ub, x, c, error)
    complex(kind=dp), dimension(lbwmax_ub+ubwmax_ub+1,n), intent(in) :: b_ub
    integer(kind=int32), dimension(n), intent(in) :: numrots_ub
    integer(kind=int32), dimension(ubwmax_ub,n), intent(in) :: js_ub
    complex(kind=dp), dimension(ubwmax_ub,n), intent(in) :: cs_ub, ss_ub
    integer(kind=int32), intent(in) :: n, lbwmax_ub, ubwmax_ub, lbw_ub, ubw_ub
    complex(kind=dp), dimension(n), intent(inout) :: c
    complex(kind=dp), dimension(n), intent(out) :: x
    type(error_info), intent(out) :: error

    integer(kind=int32) :: j, k, l
    type(c_rotation) :: rot

    if (n==1) then
       x=c/b_ub(1,1); return
    end if
    x=0.0_dp
    ! Apply all the u_k^T to c.
    do k=1,n-1
       do j=numrots_ub(k),1,-1
          rot%cosine=cs_ub(j,k); rot%sine=ss_ub(j,k)
          call rotation_times_general(trp_rot(rot), c, js_ub(j,k), js_ub(j,k)+1)
       end do
    end do
    ! diagonals of b are in row ubw_ub+1 of b_ub.
    x(n)=c(n)/b_ub(ubw_ub+1,n)
    do k=n-1,1,-1
       l=min(ubw_ub,k)       ! number of superdiagonals in column k+1
       c(k+1-l:k+1)=c(k+1-l:k+1) - x(k+1) * b_ub(ubw_ub+1-l:ubw_ub+1 ,k+1)
       ! Apply u_k to c
       do j=1,numrots_ub(k)
          rot%cosine=cs_ub(j,k); rot%sine=ss_ub(j,k)
          call rotation_times_general(rot, c, js_ub(j,k), js_ub(j,k)+1)
       end do
       x(k)=c(k)/b_ub(ubw_ub+1,k)
    end do
  end subroutine f_c_v_back_substitution_ub


end module back_substitution

