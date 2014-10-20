module assemble
use misc
use band_types
use nested_types
implicit none

interface ub_to_upper
   module procedure d_ub_to_upper, c_ub_to_upper
end interface ub_to_upper

interface f_ub_to_upper
   module procedure f_d_ub_to_upper, f_c_ub_to_upper
end interface f_ub_to_upper

interface bt_to_lower
   module procedure d_bt_to_lower!, c_bt_to_lower
end interface bt_to_lower

interface f_bt_to_lower
   module procedure f_d_bt_to_lower!, f_c_bt_to_lower
end interface f_bt_to_lower

interface bv_to_upper
   module procedure d_bv_to_upper, c_bv_to_upper
end interface bv_to_upper

interface f_bv_to_upper
   module procedure f_d_bv_to_upper, f_c_bv_to_upper
end interface f_bv_to_upper

type(routine_info), parameter :: info_d_ub_to_upper=routine_info(id_d_ub_to_upper, 'd_ub_to_upper', &
     [ character(len=error_message_length) :: 'Size error in A.' ] )
type(routine_info), parameter :: info_c_ub_to_upper=routine_info(id_c_ub_to_upper, 'c_ub_to_upper', &
     [ character(len=error_message_length) :: 'Size error in A.' ] )

type(routine_info), parameter :: info_d_bt_to_lower=routine_info(id_d_bt_to_lower, 'd_bt_to_lower', &
     [ character(len=error_message_length) :: 'Size error in A.' ] )
type(routine_info), parameter :: info_c_bt_to_lower=routine_info(id_c_bt_to_lower, 'c_bt_to_lower', &
     [ character(len=error_message_length) :: 'Size error in A.' ] )

type(routine_info), parameter :: info_d_bv_to_upper=routine_info(id_d_bv_to_upper, 'd_bv_to_upper', &
     [ character(len=error_message_length) :: 'Size error in A.' ] )
type(routine_info), parameter :: info_c_bv_to_upper=routine_info(id_c_bv_to_upper, 'c_bv_to_upper', &
     [ character(len=error_message_length) :: 'Size error in A.' ] )

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
    call f_d_ub_to_upper(ub%bc, get_n(ub), ub%lbw, ub%ubw, get_lbwmax(ub), get_ubwmax(ub), &
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

  ! Errors
  ! 0: no error
  ! 1: bt%n /= n
  subroutine d_bt_to_lower(bt,a,error)
    type(d_bt), intent(in) :: bt
    real(kind=dp), dimension(:,:), intent(out) :: a
    type(error_info), intent(out) :: error
    call clear_error(error)
    if (get_n(bt) /= size(a,1) .or. get_n(bt) /= size(a,2)) then
       call set_error(error, 1, id_d_bt_to_lower); return
    end if
    call f_d_bt_to_lower(bt%br, get_n(bt), bt%lbw, bt%ubw, get_lbwmax(bt), get_ubwmax(bt), &
         bt%numrotst, bt%kst, bt%cst, bt%sst, a)
  end subroutine d_bt_to_lower

  subroutine f_d_bt_to_lower(br, n, lbw, ubw, lbwmax, ubwmax, numrotst, kst, cst, sst, a)
    real(kind=dp), target, dimension(n,n), intent(out) :: a
    integer(kind=int32), dimension(n,lbwmax), intent(in) :: kst
    real(kind=dp), dimension(n,lbwmax), intent(in) :: cst, sst
    real(kind=dp), dimension(n,lbwmax+ubwmax+1), intent(in) :: br
    integer(kind=int32), dimension(n), intent(in) :: numrotst
    integer(kind=int32), intent(in) :: ubw, lbw, n, lbwmax, ubwmax
    !
    integer(kind=int32) :: j,k
    type(d_rotation) :: rot
    call br_to_general(br,lbw,ubw,a)
    if (n==1) then
       return
    end if
    do k=n-1,1,-1
       do j=1,numrotst(k)
          rot%cosine=cst(j,k); rot%sine=sst(j,k)
          call rotation_times_general(trp_rot(rot),a(k+1:n,:),kst(j,k),kst(j,k)+1)
       end do
    end do
  end subroutine f_d_bt_to_lower

  subroutine c_ub_to_upper(ub,a,error)
    type(c_ub), intent(in) :: ub
    complex(kind=dp), dimension(:,:), intent(out) :: a
    type(error_info), intent(out) :: error
    call clear_error(error)
    if (get_n(ub) /= size(a,1) .or. get_n(ub) /= size(a,2)) then
       call set_error(error, 1, id_c_ub_to_upper); return
    end if
    call f_c_ub_to_upper(ub%bc, get_n(ub), ub%lbw, ub%ubw, get_lbwmax(ub), get_ubwmax(ub), &
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
    call f_d_bv_to_upper(bv%br, get_n(bv), bv%lbw, bv%ubw, get_lbwmax(bv), get_ubwmax(bv), &
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
    call f_c_bv_to_upper(bv%br, get_n(bv), bv%lbw, bv%ubw, get_lbwmax(bv), get_ubwmax(bv), &
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
