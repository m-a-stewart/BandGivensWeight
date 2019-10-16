program exp
  use mod_orrb
  implicit none

  type(error_info) :: error
  integer(kind=int32) :: na, lbwa, ubwa, j
  real(kind=dp), dimension(:), allocatable :: x, rhs, res, rhs0, ur, vr, svres
  type(d_qr), allocatable :: swub
  type(d_rc), allocatable :: swbv
  type(d_ubt), allocatable :: ubta
  real(kind=dp) :: t0, t1, sigma1, sigman

  integer(kind=int32), parameter :: runs=50
  real(kind=dp) :: cond_sum, rc_time, qr_time, q_rc_time, q_qr_time, res_sum, &
       max_res, solve_time

  call initialize_errors
  ! na=100000; lbwa=5; ubwa=5
  ! na=1000000; lbwa=5; ubwa=5
  na=1000000; lbwa=5; ubwa=6
  cond_sum=0.0_dp; rc_time=0.0_dp; qr_time=0.0_dp;
  q_rc_time=0.0_dp; q_qr_time=0.0_dp; res_sum=0.0_dp
  solve_time=0.0_dp
  max_res=0.0_dp
  do j=1,runs
     print *, "j: ", j
     allocate(ur(na),vr(na),svres(na))
     rhs=d_random_vector(na)
     rhs0=rhs
     ubta=d_random_ubt(na,lbwa,ubwa,error=error)
     call cpu_time(t0)
     swbv=rc_of(ubta,error)
     call cpu_time(t1)
     rc_time=rc_time+t1-t0
     print *, "Row compression time: ", t1-t0
     call cpu_time(t0)
     rhs=trp(swbv%sw) * rhs
     call cpu_time(t1)
     q_rc_time=q_rc_time+t1-t0
     print *, "RHS times Q_rc: ", t1-t0
     call cpu_time(t0)
     swub=qr_of(swbv%bv,error)
     call cpu_time(t1)
     qr_time=qr_time+t1-t0
     print *, "QR time: ", t1-t0
     call cpu_time(t0)
     rhs=trp(swub%sw) * rhs
     call cpu_time(t1)
     q_qr_time=q_qr_time+t1-t0
     print *, "RHS times Q_qr: ", t1-t0
     call cpu_time(t0)
     x=solve(swub%ub, rhs, error)
     call cpu_time(t1)
     solve_time=solve_time+t1-t0
     print *, "solve: ", t1-t0
     sigma1=ub_max_sv(swub%ub,ur,vr,svres,1e-2_dp,1000,error)
     print *, "sigma1: ", sigma1, norm2(svres)
     sigman=ub_min_sv(swub%ub,ur,vr,svres,1e-16_dp,1000,error)
     print *, "sigman: ", sigman, norm2(svres)
     print *, norm2(swub%ub*reshape(vr,[na,1])-sigman*reshape(ur,[na,1])), norm2(swub%ub*reshape(vr,[na,1]))
     write (*, "(' Cond2: ', ES8.2)") sigma1/sigman
     cond_sum=cond_sum + sigma1/sigman
     res=reshape(ubta*reshape(x,[na,1])-reshape(rhs0,[na,1]),[na])
     print *, "Relative residual: ", norm2(res)/(sigma1*norm2(x)+norm2(rhs0))
     res_sum=res_sum+norm2(res)/(sigma1*norm2(x)+norm2(rhs0))
     max_res=max(max_res,norm2(res)/(sigma1*norm2(x)+norm2(rhs0)))
     deallocate(ur,vr,svres)
     print *
  end do
  print *, "Average residual: ", res_sum/runs
  print *, "Max residual: ", max_res
  write (*, "(' Cond2: ', ES8.2)") cond_sum/runs
  print *, "Average row compress time: ", rc_time/runs
  print *, "Averge Q_rc multiply time: ", q_rc_time/runs
  print *, "Average qr compress time: ", qr_time/runs
  print *, "Average Q_qr multiply time: ", q_qr_time/runs
  print *, "Average solve time: ", solve_time/runs
  

end program exp
