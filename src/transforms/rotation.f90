module mod_rotation
  use mod_prec
  implicit none
 ! A simple derived type for rotations, routines for computing
 ! rotations that introduce zeros, and routines for applying the
 ! rotations to matrices and vectors.
  private

  public :: d_rotation, z_rotation

  public :: lgivens, d_lgivens, z_lgivens, lgivens2, d_lgivens2, z_lgivens2, &
       rgivens, d_rgivens, z_rgivens, rgivens2, d_rgivens2, z_rgivens2

  public :: trp_rot, d_trp_rot, z_trp_rot

  public :: rotation_times_general, d_rotation_times_general, z_rotation_times_general, &
       d_rotation_times_general_v, z_rotation_times_general_v

  public :: general_times_rotation, d_general_times_rotation, z_general_times_rotation, &
       d_general_times_rotation_v, z_general_times_rotation_v

  type d_rotation
     real(kind=dp) :: cosine
     real(kind=dp) :: sine
  end type d_rotation

  type z_rotation
     real(kind=dp) :: cosine
     complex(kind=dp) :: sine
  end type z_rotation

  interface lgivens
     module procedure d_lgivens, z_lgivens
  end interface lgivens

  interface lgivens2
     module procedure d_lgivens2, z_lgivens2
  end interface lgivens2

  interface rgivens
     module procedure d_rgivens, z_rgivens
  end interface rgivens

  interface rgivens2
     module procedure d_rgivens2, z_rgivens2
  end interface rgivens2


  interface trp_rot
     module procedure d_trp_rot, z_trp_rot
  end interface trp_rot

  interface rotation_times_general
     module procedure d_rotation_times_general, z_rotation_times_general, &
          d_rotation_times_general_v, z_rotation_times_general_v
  end interface rotation_times_general

  interface general_times_rotation
     module procedure d_general_times_rotation, z_general_times_rotation, &
          d_general_times_rotation_v, z_general_times_rotation_v
  end interface general_times_rotation

contains

  ! compute a rotation r=[c,-s;s,c] such that r^T * [x;y]=[z;0]
  type(d_rotation) function d_lgivens(x,y) result(r)
    real(kind=dp), intent(in) :: x, y
    real(kind=dp) :: w, z
    z=max(abs(x), abs(y))
    if (z == 0.0_dp) then
       r%cosine=1.0_dp
       r%sine=0.0_dp
    else
       w=z*sqrt((x/z)**2+(y/z)**2)
       r%cosine=x/w
       r%sine=y/w
    end if
  end function d_lgivens

  ! rotation such that r=[c,-conj(s); s, c] and r^H * [x;y] = [z;0]
  type(z_rotation) function z_lgivens(x,y) result(r)
    complex(kind=dp), intent(in) :: x, y
    real(kind=dp) :: w, z, xabs
    xabs=abs(x)
    z=max(xabs, abs(y))
    if (z == 0.0_dp) then
       r%cosine=1.0_dp
       r%sine=(0.0_dp, 0.0_dp)
    else if (xabs == 0.0_dp) then
       r%cosine=0.0_dp
       r%sine=(1.0_dp,0.0_dp)
    else
       w=z*sqrt(abs(x/z)**2+abs(y/z)**2)
       r%cosine=xabs/w
       r%sine=y * conjg(x)/xabs/w
    end if
  end function z_lgivens

  ! compute a rotation r=[c,-s;s,c] such that r^T * [x;y]=[0;z]
  type(d_rotation) function d_lgivens2(x,y) result(r)
    real(kind=dp), intent(in) :: x, y
    real(kind=dp) :: w, z
    z=max(abs(x), abs(y))
    if (z == 0.0_dp) then
       r%cosine=1.0_dp
       r%sine=0.0_dp
    else
       w=z*sqrt((x/z)**2+(y/z)**2)
       r%cosine=y/w
       r%sine=-x/w
    end if
  end function d_lgivens2

  ! rotation such that r=[c,-conj(s); s, c] and r^H * [x;y] = [0;z]
  type(z_rotation) function z_lgivens2(x,y) result(r)
    complex(kind=dp), intent(in) :: x, y
    real(kind=dp) :: w, z, yabs
    yabs=abs(y)
    z=max(yabs, abs(x))
    if (z == 0.0_dp) then
       r%cosine=1.0_dp
       r%sine=(0.0_dp,0.0_dp)
    else if (yabs == 0.0_dp) then
       r%cosine=0.0_dp
       r%sine=(1.0_dp,0.0_dp)
    else
       w=z*sqrt(abs(x/z)**2+abs(y/z)**2)
       r%cosine=yabs/w
       r%sine=-conjg(x) * y/yabs/w
    end if
  end function z_lgivens2


  ! rotations r=[c,-s; s, c] such that [x,y] * r = [z, 0]
  type(d_rotation) function d_rgivens(x,y) result(r)
    real(kind=dp), intent(in) :: x, y
    real(kind=dp) :: w, z
    z=max(abs(x), abs(y))
    if (z == 0.0_dp) then
       r%cosine=1.0_dp
       r%sine=0.0_dp
    else
       w=z*sqrt((x/z)**2+(y/z)**2)
       r%cosine=x/w
       r%sine=y/w
    end if
  end function d_rgivens

  ! rotations r=[c,-conj(s); s, c] such that [x,y] * r = [z, 0]
  type(z_rotation) function z_rgivens(x,y) result(r)
    complex(kind=dp), intent(in) :: x, y
    real(kind=dp) :: w, z, xabs
    xabs=abs(x)
    z=max(xabs, abs(y))
    if (z == 0.0_dp) then
       r%cosine=1.0_dp
       r%sine=(0.0_dp,0.0_dp)
    else if (xabs == 0.0_dp) then
       r%cosine=0.0_dp
       r%sine=(1.0_dp,0.0_dp)
    else
       w=z*sqrt(abs(x/z)**2+abs(y/z)**2)
       r%cosine=xabs/w
       r%sine=conjg(y) * x/xabs/w
    end if
  end function z_rgivens

  ! rotations r=[c,-s; s, c] such that [x,y] * r = [0, z]
  type(d_rotation) function d_rgivens2(x,y) result(r)
    real(kind=dp), intent(in) :: x, y
    real(kind=dp) :: w, z
    z=max(abs(x), abs(y))
    if (z == 0.0_dp) then
       r%cosine=1.0_dp
       r%sine=0.0_dp
    else
       w=z*sqrt((x/z)**2+(y/z)**2)
       r%cosine=y/w
       r%sine=-x/w
    end if
  end function d_rgivens2

  ! rotations r=[c,-conj(s); s, c] such that [x,y] * r = [0, z]
  type(z_rotation) function z_rgivens2(x,y) result(r)
    complex(kind=dp), intent(in) :: x, y
    real(kind=dp) :: w, z, yabs
    yabs=abs(y)
    z=max(yabs, abs(x))
    if (z == 0.0_dp) then
       r%cosine=1.0_dp
       r%sine=(0.0_dp,0.0_dp)
    else if (yabs == 0.0_dp) then
       r%cosine=0.0_dp
       r%sine=(1.0_dp,0.0_dp)
    else
       w=z*sqrt(abs(x/z)**2+abs(y/z)**2)
       r%cosine=yabs/w
       r%sine=-x * conjg(y)/yabs/w
    end if
  end function z_rgivens2

  ! r*a with r=[c,-s;s,c]
  subroutine d_rotation_times_general(r,a,j1,j2)
    real(kind=dp), dimension(:,:), intent(inout) :: a
    type(d_rotation), intent(in) :: r
    integer(kind=int32), intent(in) :: j1, j2
    !
    integer(kind=int32) :: n, k
    real(kind=dp) :: c, s, tmp
    c=r%cosine
    s=r%sine
    n=size(a,2)
    do k=1,n
       tmp=a(j1, k)
       a(j1,k)=c*tmp-s*a(j2,k)
       a(j2,k)=s*tmp+c*a(j2,k)
    end do
  end subroutine d_rotation_times_general

  subroutine d_rotation_times_general_v(r,a,j1,j2)
    real(kind=dp), dimension(:), intent(inout) :: a
    type(d_rotation), intent(in) :: r
    integer(kind=int32), intent(in) :: j1, j2
    !
    real(kind=dp) :: c, s, tmp
    c=r%cosine
    s=r%sine
    tmp=a(j1)
    a(j1)=c*tmp-s*a(j2)
    a(j2)=s*tmp+c*a(j2)
  end subroutine d_rotation_times_general_v

  ! r*a with r=[c,-conj(s);s,c]
  subroutine z_rotation_times_general(r,a,j1,j2)
    complex(kind=dp), dimension(:,:), intent(inout) :: a
    type(z_rotation), intent(in) :: r
    integer(kind=int32), intent(in) :: j1, j2
    !
    integer(kind=int32) :: n, k
    complex(kind=dp) :: s, tmp
    real(kind=dp) :: c
    c=r%cosine
    s=r%sine
    n=size(a,2)
    do k=1,n
       tmp=a(j1, k)
       a(j1,k)=c*tmp-conjg(s)*a(j2,k)
       a(j2,k)=s*tmp+c*a(j2,k)
    end do
  end subroutine z_rotation_times_general

  subroutine z_rotation_times_general_v(r,a,j1,j2)
    complex(kind=dp), dimension(:), intent(inout) :: a
    type(z_rotation), intent(in) :: r
    integer(kind=int32), intent(in) :: j1, j2
    !
    complex(kind=dp) :: s, tmp
    real(kind=dp) :: c
    c=r%cosine
    s=r%sine
    tmp=a(j1)
    a(j1)=c*tmp-conjg(s)*a(j2)
    a(j2)=s*tmp+c*a(j2)
  end subroutine z_rotation_times_general_v

  ! a*r with r=[c,-s;s,c]
  subroutine d_general_times_rotation(a,r,k1,k2)
    real(kind=dp), dimension(:,:), intent(inout) :: a
    type(d_rotation), intent(in) :: r
    integer(kind=int32), intent(in) :: k1, k2
    !
    integer(kind=int32) :: m, j
    real(kind=dp) :: c, s, tmp
    c=r%cosine
    s=r%sine
    m=size(a,1)
    do j=1,m
       tmp=a(j,k1)
       a(j,k1)=c*tmp+s*a(j,k2)
       a(j,k2)=-s*tmp+c*a(j,k2)
    end do
  end subroutine d_general_times_rotation

  subroutine d_general_times_rotation_v(a,r,k1,k2)
    real(kind=dp), dimension(:), intent(inout) :: a
    type(d_rotation), intent(in) :: r
    integer(kind=int32), intent(in) :: k1, k2
    !
    real(kind=dp) :: c, s, tmp
    c=r%cosine
    s=r%sine
    tmp=a(k1)
    a(k1)=c*tmp+s*a(k2)
    a(k2)=-s*tmp+c*a(k2)
  end subroutine d_general_times_rotation_v


  ! a*r with r=[c,-conj(s);s,c]
  subroutine z_general_times_rotation(a,r,k1,k2)
    complex(kind=dp), dimension(:,:), intent(inout) :: a
    type(z_rotation), intent(in) :: r
    integer(kind=int32), intent(in) :: k1, k2
    !
    integer(kind=int32) :: m, j
    complex(kind=dp) :: s, tmp
    real(kind=dp) :: c
    c=r%cosine
    s=r%sine
    m=size(a,1)
    do j=1,m
       tmp=a(j,k1)
       a(j,k1)=c*tmp+s*a(j,k2)
       a(j,k2)=-conjg(s)*tmp+c*a(j,k2)
    end do
  end subroutine z_general_times_rotation

  subroutine z_general_times_rotation_v(a,r,k1,k2)
    complex(kind=dp), dimension(:), intent(inout) :: a
    type(z_rotation), intent(in) :: r
    integer(kind=int32), intent(in) :: k1, k2
    !
    complex(kind=dp) :: s, tmp
    real(kind=dp) :: c
    c=r%cosine
    s=r%sine
    tmp=a(k1)
    a(k1)=c*tmp+s*a(k2)
    a(k2)=-conjg(s)*tmp+c*a(k2)
  end subroutine z_general_times_rotation_v

  type(d_rotation) function d_trp_rot(r) result(tr)
    type(d_rotation), intent(in) :: r
    tr%cosine=r%cosine
    tr%sine=-r%sine
  end function d_trp_rot

  type(z_rotation) function z_trp_rot(r) result(tr)
    type(z_rotation), intent(in) :: r
    tr%cosine=r%cosine
    tr%sine=-r%sine
  end function z_trp_rot

end module mod_rotation
