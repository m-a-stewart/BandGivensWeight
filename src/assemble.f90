module assemble
use prec
use rotation
implicit none

interface form_upper_ub
   module procedure d_form_upper_ub, c_form_upper_ub
end interface form_upper_ub

interface form_upper_bv
   module procedure d_form_upper_bv, c_form_upper_bv
end interface form_upper_bv

contains

  subroutine d_form_upper_ub(a, n, b, mb, lbw, ubw, numrots, j1s, j2s, cs, ss)
    real(kind=dp), target, dimension(n,n), intent(out) :: a
    integer(kind=int32), dimension(mb-lbw-1,n), intent(in) :: j1s, j2s
    real(kind=dp), dimension(mb-lbw-1,n), intent(in) :: cs, ss
    real(kind=dp), dimension(mb,n), intent(in) :: b
    integer(kind=int32), dimension(n), intent(in) :: numrots
    integer(kind=int32), intent(in) :: ubw, lbw, n, mb
    !
    integer(kind=int32) :: j,k
    type(d_rotation) :: rot
    a=0.0_dp
    do k=n-1,n-lbw,-1
       a(k+1-ubw:n,k+1)=b(1:n-k+ubw,k+1)
       do j=1,numrots(k)
          rot%cosine=cs(j,k); rot%sine=ss(j,k)
          call rotation_times_general(rot,a(:,k+1:n),j1s(j,k),j2s(j,k))
       end do
    end do
    do k=n-lbw-1, ubw,-1
       a(k+1-ubw:k+1+lbw,k+1)=b(1:lbw+ubw+1,k+1)
       do j=1,numrots(k)
          rot%cosine=cs(j,k); rot%sine=ss(j,k)
          call rotation_times_general(rot,a(:,k+1:n),j1s(j,k),j2s(j,k))
       end do
    end do
    do k=ubw-1,1,-1
       a(1:lbw+k+1,k+1)=b(ubw-k+1:ubw+lbw+1,k+1)
       do j=1,numrots(k)
          rot%cosine=cs(j,k); rot%sine=ss(j,k)
          call rotation_times_general(rot,a(:,k+1:n),j1s(j,k),j2s(j,k))
       end do
    end do
    a(1:lbw+1,1)=b(ubw+1:ubw+lbw+1,1)
  end subroutine d_form_upper_ub

  subroutine c_form_upper_ub(a, n, b, mb, lbw, ubw, numrots, j1s, j2s, cs, ss)
    complex(kind=dp), target, dimension(n,n), intent(out) :: a
    integer(kind=int32), dimension(mb-lbw-1,n), intent(in) :: j1s, j2s
    complex(kind=dp), dimension(mb-lbw-1,n), intent(in) :: cs, ss
    complex(kind=dp), dimension(mb,n), intent(in) :: b
    integer(kind=int32), dimension(n), intent(in) :: numrots
    integer(kind=int32), intent(in) :: ubw, lbw, n, mb
    !
    integer(kind=int32) :: j,k
    type(c_rotation) :: rot
    a=(0.0_dp, 0.0_dp)
    do k=n-1,n-lbw,-1
       a(k+1-ubw:n,k+1)=b(1:n-k+ubw,k+1)
       do j=1,numrots(k)
          rot%cosine=cs(j,k); rot%sine=ss(j,k)
          call rotation_times_general(rot,a(:,k+1:n),j1s(j,k),j2s(j,k))
       end do
    end do
    do k=n-lbw-1, ubw,-1
       a(k+1-ubw:k+1+lbw,k+1)=b(1:lbw+ubw+1,k+1)
       do j=1,numrots(k)
          rot%cosine=cs(j,k); rot%sine=ss(j,k)
          call rotation_times_general(rot,a(:,k+1:n),j1s(j,k),j2s(j,k))
       end do
    end do
    do k=ubw-1,1,-1
       a(1:lbw+k+1,k+1)=b(ubw-k+1:ubw+lbw+1,k+1)
       do j=1,numrots(k)
          rot%cosine=cs(j,k); rot%sine=ss(j,k)
          call rotation_times_general(rot,a(:,k+1:n),j1s(j,k),j2s(j,k))
       end do
    end do
    a(1:lbw+1,1)=b(ubw+1:ubw+lbw+1,1)
  end subroutine c_form_upper_ub

  subroutine d_form_upper_bv(a, n, b, nb, lbw, ubw, numrots, k1s, k2s, cs, ss)
    real(kind=dp), target, dimension(n,n), intent(out) :: a
    integer(kind=int32), dimension(n,nb-lbw-1), intent(in) :: k1s, k2s
    real(kind=dp), dimension(n, nb-lbw-1), intent(in) :: cs, ss
    real(kind=dp), dimension(n,nb), intent(in) :: b
    integer(kind=int32), dimension(n), intent(in) :: numrots
    integer(kind=int32), intent(in) :: ubw, lbw, n, nb
    !
    integer(kind=int32) :: j,k
    type(d_rotation) :: rot
    a=0.0_dp
    do j=1,lbw
       a(j,1:ubw+j)=b(j,lbw-j+2:lbw+ubw+1)
       do k=1,numrots(n-j)
          rot%cosine=cs(n-j,k); rot%sine=ss(n-j,k)
          call general_times_rotation(a(1:j,:),trp_rot(rot),k1s(n-j,k), k2s(n-j,k))
       end do
    end do
    do j=lbw+1,n-ubw
       a(j,j-lbw:j+ubw)=b(j,1:lbw+ubw+1)
       do k=1,numrots(n-j)
          rot%cosine=cs(n-j,k); rot%sine=ss(n-j,k)
          call general_times_rotation(a(1:j,:),trp_rot(rot),k1s(n-j,k), k2s(n-j,k))
       end do
    end do
    do j=n-ubw+1,n-1
       a(j,j-lbw:n)=b(j,1:n-j+lbw+1)
       do k=1,numrots(n-j)
          rot%cosine=cs(n-j,k); rot%sine=ss(n-j,k)
          call general_times_rotation(a(1:j,:),trp_rot(rot),k1s(n-j,k), k2s(n-j,k))
       end do
    end do
    a(n,n-lbw:n)=b(n,1:lbw+1)
  end subroutine d_form_upper_bv

  subroutine c_form_upper_bv(a, n, b, nb, lbw, ubw, numrots, k1s, k2s, cs, ss)
    complex(kind=dp), target, dimension(n,n), intent(out) :: a
    integer(kind=int32), dimension(n,nb-lbw-1), intent(in) :: k1s, k2s
    complex(kind=dp), dimension(n, nb-lbw-1), intent(in) :: cs, ss
    complex(kind=dp), dimension(n,nb), intent(in) :: b
    integer(kind=int32), dimension(n), intent(in) :: numrots
    integer(kind=int32), intent(in) :: ubw, lbw, n, nb
    !
    integer(kind=int32) :: j,k
    type(c_rotation) :: rot
    a=(0.0_dp, 0.0_dp)
    do j=1,lbw
       a(j,1:ubw+j)=b(j,lbw-j+2:lbw+ubw+1)
       do k=1,numrots(n-j)
          rot%cosine=cs(n-j,k); rot%sine=ss(n-j,k)
          call general_times_rotation(a(1:j,:),trp_rot(rot),k1s(n-j,k), k2s(n-j,k))
       end do
    end do
    do j=lbw+1,n-ubw
       a(j,j-lbw:j+ubw)=b(j,1:lbw+ubw+1)
       do k=1,numrots(n-j)
          rot%cosine=cs(n-j,k); rot%sine=ss(n-j,k)
          call general_times_rotation(a(1:j,:),trp_rot(rot),k1s(n-j,k), k2s(n-j,k))
       end do
    end do
    do j=n-ubw+1,n-1
       a(j,j-lbw:n)=b(j,1:n-j+lbw+1)
       do k=1,numrots(n-j)
          rot%cosine=cs(n-j,k); rot%sine=ss(n-j,k)
          call general_times_rotation(a(1:j,:),trp_rot(rot),k1s(n-j,k), k2s(n-j,k))
       end do
    end do
    a(n,n-lbw:n)=b(n,1:lbw+1)
  end subroutine c_form_upper_bv
    
end module assemble
