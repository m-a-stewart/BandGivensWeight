program exp
  use mod_orb
  implicit none

  type(error_info) :: error
  integer(kind=int32) :: na, lbwa, ubwa
  real(kind=dp), dimension(:,:), allocatable :: a, r
  real(kind=dp), dimension(:), allocatable :: x, rhs, res, rhs0, ur, vr, svres
  type(d_qr), allocatable :: swub
  type(d_rc), allocatable :: swbv
  type(d_ubt), allocatable :: ubta
  real(kind=dp) :: t0, t1, scale, sigma1, sigman

  na=2000; lbwa=5; ubwa=5
  allocate(ur(na),vr(na),svres(na))
  rhs=d_random_vector(na)
  rhs0=rhs
  ubta=d_random_ubt(na,lbwa,ubwa,error=error)
  a = general(ubta,error)
  swbv=rc(ubta,error)
  rhs=trp(swbv%sw) * rhs
  swub=qr(swbv%bv,error)
  r=general(swub%ub,error)
  sigma1=d_upper_max_sv(r,ur,vr,svres,1e-8_dp,1000,error)
  print *, sigma1, norm2(svres)
  sigman=d_upper_min_sv(r,ur,vr,svres,1e-16_dp,1000,error)
  print *, sigman, norm2(svres)
  rhs=trp(swub%sw) * rhs
  x=solve(swub%ub, rhs, error)
  res=matmul(a,x)-rhs0
  print *, "Relative residual: ", norm2(res)/(sigma1*norm2(x)+norm2(rhs0))
  
end program exp
