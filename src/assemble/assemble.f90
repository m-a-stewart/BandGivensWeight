module assemble
use prec
use error_id
use rotation
use band_types
use nested_types
implicit none

interface ub_to_upper
   module procedure d_ub_to_upper, c_ub_to_upper
end interface ub_to_upper

interface f_ub_to_upper
   module procedure f_d_ub_to_upper, f_c_ub_to_upper
end interface f_ub_to_upper

interface bv_to_upper
   module procedure d_bv_to_upper, c_bv_to_upper
end interface bv_to_upper

interface f_bv_to_upper
   module procedure f_d_bv_to_upper, f_c_bv_to_upper
end interface f_bv_to_upper

type(routine_info), parameter :: info_d_ub_to_upper=routine_info(id_d_ub_to_upper, 'd_ub_to_upper', &
     [ character(len=error_message_length) :: 'Size error in A.' ] )
type(routine_info), parameter :: info_f_d_ub_to_upper=routine_info(id_f_d_ub_to_upper, 'f_d_ub_to_upper', &
     [ character(len=error_message_length) :: '' ] )
type(routine_info), parameter :: info_c_ub_to_upper=routine_info(id_c_ub_to_upper, 'c_ub_to_upper', &
     [ character(len=error_message_length) :: 'Size error in A.' ] )
type(routine_info), parameter :: info_f_c_ub_to_upper=routine_info(id_f_c_ub_to_upper, 'f_c_ub_to_upper', &
     [ character(len=error_message_length) :: '' ] )

type(routine_info), parameter :: info_d_bv_to_upper=routine_info(id_d_bv_to_upper, 'd_bv_to_upper', &
     [ character(len=error_message_length) :: 'Size error in A.' ] )
type(routine_info), parameter :: info_f_d_bv_to_upper=routine_info(id_f_d_bv_to_upper, 'f_d_bv_to_upper', &
     [ character(len=error_message_length) :: '' ] )
type(routine_info), parameter :: info_c_bv_to_upper=routine_info(id_c_bv_to_upper, 'c_bv_to_upper', &
     [ character(len=error_message_length) :: 'Size error in A.' ] )
type(routine_info), parameter :: info_f_c_bv_to_upper=routine_info(id_f_c_bv_to_upper, 'f_c_bv_to_upper', &
     [ character(len=error_message_length) :: '' ] )

contains

  ! Errors
  ! 0: no error
  ! 1: ub%n /= n
  subroutine d_ub_to_upper(ub,a,error)
    type(d_ub), intent(in) :: ub
    real(kind=dp), dimension(:,:), intent(out) :: a
    type(error_info), intent(out) :: error
    call clear_error(error)
    if (get_n(ub) /= size(a,1) .or. get_n(ub) /= size(a,2)) then
       call set_error(error, 1, id_d_ub_to_upper); return
    end if
    call f_d_ub_to_upper(ub%b, get_n(ub), ub%lbw, ub%ubw, get_lbwmax(ub), get_ubwmax(ub), &
         ub%numrotsu, ub%jsu, ub%csu, ub%ssu, a)
  end subroutine d_ub_to_upper

  subroutine f_d_ub_to_upper(b, n, lbw, ubw, lbwmax, ubwmax, numrots, js, cs, ss, a)
    real(kind=dp), target, dimension(n,n), intent(out) :: a
    integer(kind=int32), dimension(ubwmax,n), intent(in) :: js
    real(kind=dp), dimension(ubwmax,n), intent(in) :: cs, ss
    real(kind=dp), dimension(lbwmax+ubwmax+1,n), intent(in) :: b
    integer(kind=int32), dimension(n), intent(in) :: numrots
    integer(kind=int32), intent(in) :: ubw, lbw, n, lbwmax, ubwmax
    !
    integer(kind=int32) :: j,k
    type(d_rotation) :: rot
    call bc_to_general(b,lbw,ubw,a)
    if (n==1) then
       return
    end if
    do k=n-1,1,-1
       do j=1,numrots(k)
          rot%cosine=cs(j,k); rot%sine=ss(j,k)
          call rotation_times_general(rot,a(:,k+1:n),js(j,k),js(j,k)+1)
       end do
    end do
  end subroutine f_d_ub_to_upper

  subroutine c_ub_to_upper(ub,a,error)
    type(c_ub), intent(in) :: ub
    complex(kind=dp), dimension(:,:), intent(out) :: a
    type(error_info), intent(out) :: error
    call clear_error(error)
    if (get_n(ub) /= size(a,1) .or. get_n(ub) /= size(a,2)) then
       call set_error(error, 1, id_c_ub_to_upper); return
    end if
    call f_c_ub_to_upper(ub%b, get_n(ub), ub%lbw, ub%ubw, get_lbwmax(ub), get_ubwmax(ub), &
         ub%numrotsu, ub%jsu, ub%csu, ub%ssu, a)
  end subroutine c_ub_to_upper

  subroutine f_c_ub_to_upper(b, n, lbw, ubw, lbwmax, ubwmax, numrots, js, cs, ss, a)
    complex(kind=dp), target, dimension(n,n), intent(out) :: a
    integer(kind=int32), dimension(ubwmax,n), intent(in) :: js
    complex(kind=dp), dimension(ubwmax,n), intent(in) :: cs, ss
    complex(kind=dp), dimension(lbwmax+ubwmax+1,n), intent(in) :: b
    integer(kind=int32), dimension(n), intent(in) :: numrots
    integer(kind=int32), intent(in) :: ubw, lbw, n, lbwmax, ubwmax
    !
    integer(kind=int32) :: j,k
    type(c_rotation) :: rot

    call bc_to_general(b,lbw,ubw,a)
    if (n==1) then
       return
    end if
    do k=n-1,1,-1
       do j=1,numrots(k)
          rot%cosine=cs(j,k); rot%sine=ss(j,k)
          call rotation_times_general(rot,a(:,k+1:n),js(j,k),js(j,k)+1)
       end do
    end do
  end subroutine f_c_ub_to_upper

  subroutine d_bv_to_upper(bv,a,error)
    type(d_bv), intent(in) :: bv
    real(kind=dp), dimension(:,:), intent(out) :: a
    type(error_info), intent(out) :: error
    call clear_error(error)
    if (get_n(bv) /= size(a,1) .or. get_n(bv) /= size(a,2)) then
       call set_error(error, 1, id_d_bv_to_upper); return
    end if
    call f_d_bv_to_upper(bv%b, get_n(bv), bv%lbw, bv%ubw, get_lbwmax(bv), get_ubwmax(bv), &
         bv%numrotsv, bv%ksv, bv%csv, bv%ssv, a)
  end subroutine d_bv_to_upper

  subroutine f_d_bv_to_upper(b, n, lbw, ubw, lbwmax, ubwmax, numrots, ks, cs, ss, a)
    real(kind=dp), target, dimension(n,n), intent(out) :: a
    integer(kind=int32), dimension(n,ubwmax), intent(in) :: ks
    real(kind=dp), dimension(n, ubwmax), intent(in) :: cs, ss
    real(kind=dp), dimension(n,lbwmax+ubwmax+1), intent(in) :: b
    integer(kind=int32), dimension(n), intent(in) :: numrots
    integer(kind=int32), intent(in) :: ubw, lbw, n, lbwmax, ubwmax
    !
    integer(kind=int32) :: j,k
    type(d_rotation) :: rot
    call br_to_general(b,lbw,ubw,a)
    if (n==1) then
       return
    end if
    do j=1,n-2
       do k=1,numrots(n-j)
          rot%cosine=cs(n-j,k); rot%sine=ss(n-j,k)
          call general_times_rotation(a(1:j,:),trp_rot(rot),ks(n-j,k), ks(n-j,k)+1)
       end do
    end do
  end subroutine f_d_bv_to_upper

  subroutine c_bv_to_upper(bv,a,error)
    type(c_bv), intent(in) :: bv
    complex(kind=dp), dimension(:,:), intent(out) :: a
    type(error_info), intent(out) :: error
    call clear_error(error)
    if (get_n(bv) /= size(a,1) .or. get_n(bv) /= size(a,2)) then
       call set_error(error, 1, id_c_bv_to_upper); return
    end if
    call f_c_bv_to_upper(bv%b, get_n(bv), bv%lbw, bv%ubw, get_lbwmax(bv), get_ubwmax(bv), &
         bv%numrotsv, bv%ksv, bv%csv, bv%ssv, a)
  end subroutine c_bv_to_upper

  subroutine f_c_bv_to_upper(b, n, lbw, ubw, lbwmax, ubwmax, numrots, ks, cs, ss, a)
    complex(kind=dp), target, dimension(n,n), intent(out) :: a
    integer(kind=int32), dimension(n,ubwmax), intent(in) :: ks
    complex(kind=dp), dimension(n, ubwmax), intent(in) :: cs, ss
    complex(kind=dp), dimension(n,lbwmax+ubwmax+1), intent(in) :: b
    integer(kind=int32), dimension(n), intent(in) :: numrots
    integer(kind=int32), intent(in) :: ubw, lbw, n, lbwmax, ubwmax
    !
    integer(kind=int32) :: j,k
    type(c_rotation) :: rot
    call br_to_general(b,lbw,ubw,a)
    if (n==1) then
       return
    end if
    do j=1,n-2
       do k=1,numrots(n-j)
          rot%cosine=cs(n-j,k); rot%sine=ss(n-j,k)
          call general_times_rotation(a(1:j,:),trp_rot(rot),ks(n-j,k), ks(n-j,k)+1)
       end do
    end do
  end subroutine f_c_bv_to_upper
    
end module assemble
