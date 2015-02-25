program test_solve
  use mod_orth_rank
  use mod_test_data
  implicit none
  real(kind=dp) :: t0, t1, scale
  integer(kind=int32) :: na, lbwa, ubwa, nc
  type(error_info) :: error
  real(kind=dp), parameter :: tol=1e-15, c=5.0
  !

  real(kind=dp), dimension(:,:), allocatable :: a_d, x_d, rhs_d, rhs0_d
  complex(kind=dp), dimension(:,:), allocatable :: a_c, x_c, rhs_c, rhs0_c

  real(kind=dp), dimension(:), allocatable :: rhs_v_d, rhs0_v_d, x_v_d
  complex(kind=dp), dimension(:), allocatable :: rhs_v_c, rhs0_v_c, x_v_c

  type(d_ubt), allocatable :: ubt_d
  type(c_ubt), allocatable :: ubt_c
  type(d_rc), allocatable :: swbv_d
  type(c_rc), allocatable :: swbv_c
  type(d_qr), allocatable :: swub_d
  type(c_qr), allocatable :: swub_c
  type(d_bv), allocatable :: bv_d
  type(c_bv), allocatable :: bv_c
  
  call initialize_errors
  print *
  print *, "--------------------------------"
  print *
  print *, "Back Solver Tests (Timings for back subs.)"
  print *

  na=40; lbwa=5; ubwa=7; nc=3
  rhs_d=d_random_matrix(na,nc)
  rhs0_d=rhs_d
  ubt_d=d_random_ubt(na,lbwa,ubwa,error=error)
  a_d = general(ubt_d,error)
  swbv_d=rc(ubt_d,error)
  rhs_d=trp(swbv_d%sw) * rhs_d
  swub_d=qr(swbv_d%bv,error)
  rhs_d=trp(swub_d%sw) * rhs_d
  call cpu_time(t0)
  x_d=solve(swub_d%ub, rhs_d, error)
  call cpu_time(t1)
  rhs_d=matmul(a_d,x_d)
  scale=maxabs(a_d)*maxabs(x_d)
  test_name = "Real Back Subs. (n=40)"
  call d_output_result_upper(test_name,rhs0_d/scale,rhs_d/scale,ubwa,ubwa,t0,t1,c*tol,error)

  na=40; lbwa=5; ubwa=7
  rhs_v_d=d_random_vector(na)
  rhs0_v_d=rhs_v_d
  ubt_d=d_random_ubt(na,lbwa,ubwa,error=error)
  a_d = general(ubt_d,error)
  swbv_d=rc(ubt_d,error)
  rhs_v_d=trp(swbv_d%sw) * rhs_v_d
  swub_d=qr(swbv_d%bv,error)
  rhs_v_d=trp(swub_d%sw) * rhs_v_d
  call cpu_time(t0)
  x_v_d=solve(swub_d%ub, rhs_v_d, error)
  call cpu_time(t1)
  rhs_v_d=matmul(a_d,x_v_d)
  scale=maxabs(a_d)*maxabs(x_v_d)
  test_name = "Real Vector Back Subs. (n=40)"
  call d_output_result_upper(test_name,reshape(rhs0_v_d/scale,[na,1]),&
       reshape(rhs_v_d/scale,[na,1]),ubwa,ubwa,t0,t1,c*tol,error)

  na=40; lbwa=5; ubwa=7; nc=3
  rhs_c=c_random_matrix(na,nc)
  rhs0_c=rhs_c
  ubt_c=c_random_ubt(na,lbwa,ubwa,error=error)
  a_c = general(ubt_c,error)
  swbv_c=rc(ubt_c,error)
  rhs_c=trp(swbv_c%sw) * rhs_c
  swub_c=qr(swbv_c%bv,error)
  rhs_c=trp(swub_c%sw) * rhs_c
  call cpu_time(t0)
  x_c=solve(swub_c%ub, rhs_c, error)
  call cpu_time(t1)
  rhs_c=matmul(a_c,x_c)
  scale=maxabs(a_c)*maxabs(x_c)
  test_name = "Complex Back Subs. (n=40)"
  call c_output_result_upper(test_name,rhs0_c/scale,rhs_c/scale,ubwa,ubwa,t0,t1,c*tol,error)

  na=40; lbwa=5; ubwa=7
  rhs_v_c=c_random_vector(na)
  rhs0_v_c=rhs_v_c
  ubt_c=c_random_ubt(na,lbwa,ubwa,error=error)
  a_c = general(ubt_c,error)
  swbv_c=rc(ubt_c,error)
  rhs_v_c=trp(swbv_c%sw) * rhs_v_c
  swub_c=qr(swbv_c%bv,error)
  rhs_v_c=trp(swub_c%sw) * rhs_v_c
  call cpu_time(t0)
  x_v_c=solve(swub_c%ub, rhs_v_c, error)
  call cpu_time(t1)
  rhs_v_c=matmul(a_c,x_v_c)
  scale=maxabs(a_c)*maxabs(x_v_c)
  test_name = "Complex Vector Back Subs. (n=40)"
  call c_output_result_upper(test_name,reshape(rhs0_v_c/scale,[na,1]),&
       reshape(rhs_v_c/scale,[na,1]),ubwa,ubwa,t0,t1,c*tol,error)

  print *
  print *, "--------------------------------"
  print *
  print *, "Complex Forward Solver Tests (Timings for forward subs.)"
  print *

  na=40; lbwa=5; ubwa=7; nc=3
  rhs_d=d_random_matrix(nc,na)
  rhs0_d=rhs_d
  ubt_d=d_random_ubt(na,lbwa,ubwa,error=error)
  a_d = general(ubt_d,error)
  swbv_d=rc(ubt_d,error)
  swub_d=qr(swbv_d%bv,error)
  bv_d=bv(swub_d%ub,error)
  call cpu_time(t0)
  x_d=solve(bv_d, rhs_d, error)
  call cpu_time(t1)
  x_d = x_d * trp(swub_d%sw)
  x_d = x_d * trp(swbv_d%sw)
  rhs_d=matmul(x_d,a_d)
  scale=maxabs(a_d)*maxabs(x_d)
  test_name = "Real Forward Subs. (n=40)"
  call d_output_result_upper(test_name,rhs0_d/scale,rhs_d/scale,ubwa,ubwa,t0,t1,c*tol,error)

  na=40; lbwa=5; ubwa=7
  rhs_v_d=d_random_vector(na)
  rhs0_v_d=rhs_v_d
  ubt_d=d_random_ubt(na,lbwa,ubwa,error=error)
  a_d = general(ubt_d,error)
  swbv_d=rc(ubt_d,error)
  swub_d=qr(swbv_d%bv,error)
  bv_d=bv(swub_d%ub,error)
  call cpu_time(t0)
  x_v_d=solve(bv_d, rhs_v_d, error)
  call cpu_time(t1)
  x_v_d = x_v_d * trp(swub_d%sw)
  x_v_d = x_v_d * trp(swbv_d%sw)
  rhs_v_d=matmul(x_v_d,a_d)
  scale=maxabs(a_d)*maxabs(x_v_d)
  test_name = "Real Vector Forward Subs. (n=40)"
  call d_output_result_upper(test_name,reshape(rhs0_v_d/scale,[1,na]),&
       reshape(rhs_v_d/scale,[1,na]),ubwa,ubwa,t0,t1,c*tol,error)

  na=40; lbwa=5; ubwa=7; nc=3
  rhs_c=c_random_matrix(nc,na)
  rhs0_c=rhs_c
  ubt_c=c_random_ubt(na,lbwa,ubwa,error=error)
  a_c = general(ubt_c,error)
  swbv_c=rc(ubt_c,error)
  swub_c=qr(swbv_c%bv,error)
  bv_c=bv(swub_c%ub,error)
  call cpu_time(t0)
  x_c=solve(bv_c, rhs_c, error)
  call cpu_time(t1)
  x_c = x_c * trp(swub_c%sw)
  x_c = x_c * trp(swbv_c%sw)
  rhs_c=matmul(x_c,a_c)
  scale=maxabs(a_c)*maxabs(x_c)
  test_name = "Complex Forward Subs. (n=40)"
  call c_output_result_upper(test_name,rhs0_c/scale,rhs_c/scale,ubwa,ubwa,t0,t1,c*tol,error)

  na=40; lbwa=5; ubwa=7
  rhs_v_c=c_random_vector(na)
  rhs0_v_c=rhs_v_c
  ubt_c=c_random_ubt(na,lbwa,ubwa,error=error)
  a_c = general(ubt_c,error)
  swbv_c=rc(ubt_c,error)
  swub_c=qr(swbv_c%bv,error)
  bv_c=bv(swub_c%ub,error)
  call cpu_time(t0)
  x_v_c=solve(bv_c, rhs_v_c, error)
  call cpu_time(t1)
  x_v_c = x_v_c * trp(swub_c%sw)
  x_v_c = x_v_c * trp(swbv_c%sw)
  rhs_v_c=matmul(x_v_c,a_c)
  scale=maxabs(a_c)*maxabs(x_v_c)
  test_name = "Complex Vector Forward Subs. (n=40)"
  call c_output_result_upper(test_name,reshape(rhs0_v_c/scale,[1,na]),&
       reshape(rhs_v_c/scale,[1,na]),ubwa,ubwa,t0,t1,c*tol,error)
  
end program test_solve
