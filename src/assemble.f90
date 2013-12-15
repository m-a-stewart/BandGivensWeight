module assemble
use prec
use rotation
use band
use decomp_types
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

contains

  ! Errors
  ! 0: no error
  ! 3: ub%n /= n
  subroutine d_ub_to_upper(ub,a,error)
    type(d_ub), intent(in) :: ub
    real(kind=dp), dimension(:,:), intent(out) :: a
    integer(kind=int32), intent(out) :: error
    error=0
    if (ub%n /= size(a,1) .or. ub%n /= size(a,2)) then
       error=3; return
    end if
    call f_d_ub_to_upper(ub%b, ub%n, ub%lbw, ub%ubw, ub%lbwmax, ub%ubwmax, ub%numrotsu, &
         ub%j1su, ub%j2su, ub%csu, ub%ssu, a)
  end subroutine d_ub_to_upper

  subroutine f_d_ub_to_upper(b, n, lbw, ubw, lbwmax, ubwmax, numrots, j1s, j2s, cs, ss, a)
    real(kind=dp), target, dimension(n,n), intent(out) :: a
    integer(kind=int32), dimension(ubwmax,n), intent(in) :: j1s, j2s
    real(kind=dp), dimension(ubwmax,n), intent(in) :: cs, ss
    real(kind=dp), dimension(lbwmax+ubwmax+1,n), intent(in) :: b
    integer(kind=int32), dimension(n), intent(in) :: numrots
    integer(kind=int32), intent(in) :: ubw, lbw, n, lbwmax, ubwmax
    !
    integer(kind=int32) :: j,k
    type(d_rotation) :: rot
    call bc_to_general(b,lbw,ubw,a)
    do k=n-1,n-lbw,-1
       do j=1,numrots(k)
          rot%cosine=cs(j,k); rot%sine=ss(j,k)
          call rotation_times_general(rot,a(:,k+1:n),j1s(j,k),j2s(j,k))
       end do
    end do
    do k=n-lbw-1, ubw,-1
       do j=1,numrots(k)
          rot%cosine=cs(j,k); rot%sine=ss(j,k)
          call rotation_times_general(rot,a(:,k+1:n),j1s(j,k),j2s(j,k))
       end do
    end do
    do k=ubw-1,1,-1
       do j=1,numrots(k)
          rot%cosine=cs(j,k); rot%sine=ss(j,k)
          call rotation_times_general(rot,a(:,k+1:n),j1s(j,k),j2s(j,k))
       end do
    end do
  end subroutine f_d_ub_to_upper

  subroutine c_ub_to_upper(ub,a,error)
    type(c_ub), intent(in) :: ub
    complex(kind=dp), dimension(:,:), intent(out) :: a
    integer(kind=int32), intent(out) :: error
    error=0
    if (ub%n /= size(a,1) .or. ub%n /= size(a,2)) then
       error=3; return
    end if
    call f_c_ub_to_upper(ub%b, ub%n, ub%lbw, ub%ubw, ub%lbwmax, ub%ubwmax, ub%numrotsu, &
         ub%j1su, ub%j2su, ub%csu, ub%ssu, a)
  end subroutine c_ub_to_upper

  subroutine f_c_ub_to_upper(b, n, lbw, ubw, lbwmax, ubwmax, numrots, j1s, j2s, cs, ss, a)
    complex(kind=dp), target, dimension(n,n), intent(out) :: a
    integer(kind=int32), dimension(ubwmax,n), intent(in) :: j1s, j2s
    complex(kind=dp), dimension(ubwmax,n), intent(in) :: cs, ss
    complex(kind=dp), dimension(lbwmax+ubwmax+1,n), intent(in) :: b
    integer(kind=int32), dimension(n), intent(in) :: numrots
    integer(kind=int32), intent(in) :: ubw, lbw, n, lbwmax, ubwmax
    !
    integer(kind=int32) :: j,k
    type(c_rotation) :: rot

    call bc_to_general(b,lbw,ubw,a)
    do k=n-1,n-lbw,-1
       do j=1,numrots(k)
          rot%cosine=cs(j,k); rot%sine=ss(j,k)
          call rotation_times_general(rot,a(:,k+1:n),j1s(j,k),j2s(j,k))
       end do
    end do
    do k=n-lbw-1, ubw,-1
       do j=1,numrots(k)
          rot%cosine=cs(j,k); rot%sine=ss(j,k)
          call rotation_times_general(rot,a(:,k+1:n),j1s(j,k),j2s(j,k))
       end do
    end do
    do k=ubw-1,1,-1
       do j=1,numrots(k)
          rot%cosine=cs(j,k); rot%sine=ss(j,k)
          call rotation_times_general(rot,a(:,k+1:n),j1s(j,k),j2s(j,k))
       end do
    end do
  end subroutine f_c_ub_to_upper

  subroutine d_bv_to_upper(bv,a,error)
    type(d_bv), intent(in) :: bv
    real(kind=dp), dimension(:,:), intent(out) :: a
    integer(kind=int32), intent(out) :: error
    error=0
    if (bv%n /= size(a,1) .or. bv%n /= size(a,2)) then
       error=3; return
    end if
    call f_d_bv_to_upper(bv%b, bv%n, bv%lbw, bv%ubw, bv%lbwmax, bv%ubwmax, bv%numrotsv, &
         bv%k1sv, bv%k2sv, bv%csv, bv%ssv, a)
  end subroutine d_bv_to_upper

  subroutine f_d_bv_to_upper(b, n, lbw, ubw, lbwmax, ubwmax, numrots, k1s, k2s, cs, ss, a)
    real(kind=dp), target, dimension(n,n), intent(out) :: a
    integer(kind=int32), dimension(n,ubwmax), intent(in) :: k1s, k2s
    real(kind=dp), dimension(n, ubwmax), intent(in) :: cs, ss
    real(kind=dp), dimension(n,lbwmax+ubwmax+1), intent(in) :: b
    integer(kind=int32), dimension(n), intent(in) :: numrots
    integer(kind=int32), intent(in) :: ubw, lbw, n, lbwmax, ubwmax
    !
    integer(kind=int32) :: j,k
    type(d_rotation) :: rot
    call br_to_general(b,lbw,ubw,a)
    do j=1,lbw
       do k=1,numrots(n-j)
          rot%cosine=cs(n-j,k); rot%sine=ss(n-j,k)
          call general_times_rotation(a(1:j,:),trp_rot(rot),k1s(n-j,k), k2s(n-j,k))
       end do
    end do
    do j=lbw+1,n-ubw
       do k=1,numrots(n-j)
          rot%cosine=cs(n-j,k); rot%sine=ss(n-j,k)
          call general_times_rotation(a(1:j,:),trp_rot(rot),k1s(n-j,k), k2s(n-j,k))
       end do
    end do
    do j=n-ubw+1,n-1
       do k=1,numrots(n-j)
          rot%cosine=cs(n-j,k); rot%sine=ss(n-j,k)
          call general_times_rotation(a(1:j,:),trp_rot(rot),k1s(n-j,k), k2s(n-j,k))
       end do
    end do
  end subroutine f_d_bv_to_upper

  subroutine c_bv_to_upper(bv,a,error)
    type(c_bv), intent(in) :: bv
    complex(kind=dp), dimension(:,:), intent(out) :: a
    integer(kind=int32), intent(out) :: error
    error=0
    if (bv%n /= size(a,1) .or. bv%n /= size(a,2)) then
       error=3; return
    end if
    call f_c_bv_to_upper(bv%b, bv%n, bv%lbw, bv%ubw, bv%lbwmax, bv%ubwmax, bv%numrotsv, &
         bv%k1sv, bv%k2sv, bv%csv, bv%ssv, a)
  end subroutine c_bv_to_upper

  subroutine f_c_bv_to_upper(b, n, lbw, ubw, lbwmax, ubwmax, numrots, k1s, k2s, cs, ss, a)
    complex(kind=dp), target, dimension(n,n), intent(out) :: a
    integer(kind=int32), dimension(n,ubwmax), intent(in) :: k1s, k2s
    complex(kind=dp), dimension(n, ubwmax), intent(in) :: cs, ss
    complex(kind=dp), dimension(n,lbwmax+ubwmax+1), intent(in) :: b
    integer(kind=int32), dimension(n), intent(in) :: numrots
    integer(kind=int32), intent(in) :: ubw, lbw, n, lbwmax, ubwmax
    !
    integer(kind=int32) :: j,k
    type(c_rotation) :: rot
    call br_to_general(b,lbw,ubw,a)
    do j=1,lbw
       do k=1,numrots(n-j)
          rot%cosine=cs(n-j,k); rot%sine=ss(n-j,k)
          call general_times_rotation(a(1:j,:),trp_rot(rot),k1s(n-j,k), k2s(n-j,k))
       end do
    end do
    do j=lbw+1,n-ubw
       do k=1,numrots(n-j)
          rot%cosine=cs(n-j,k); rot%sine=ss(n-j,k)
          call general_times_rotation(a(1:j,:),trp_rot(rot),k1s(n-j,k), k2s(n-j,k))
       end do
    end do
    do j=n-ubw+1,n-1
       do k=1,numrots(n-j)
          rot%cosine=cs(n-j,k); rot%sine=ss(n-j,k)
          call general_times_rotation(a(1:j,:),trp_rot(rot),k1s(n-j,k), k2s(n-j,k))
       end do
    end do
  end subroutine f_c_bv_to_upper
    
end module assemble
