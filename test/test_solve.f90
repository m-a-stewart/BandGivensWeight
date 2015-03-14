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
  complex(kind=dp), dimension(:,:), allocatable :: a_z, x_z, rhs_z, rhs0_z

  real(kind=dp), dimension(:), allocatable :: rhs_v_d, rhs0_v_d, x_v_d
  complex(kind=dp), dimension(:), allocatable :: rhs_v_z, rhs0_v_z, x_v_z

  type(d_ubt), allocatable :: ubt_d
  type(z_ubt), allocatable :: ubt_z
  type(d_rc), allocatable :: swbv_d
  type(z_rc), allocatable :: swbv_z
  type(d_qr), allocatable :: swub_d
  type(z_qr), allocatable :: swub_z
  type(d_bv), allocatable :: bv_d
  type(z_bv), allocatable :: bv_z
  
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
  rhs_z=z_random_matrix(na,nc)
  rhs0_z=rhs_z
  ubt_z=z_random_ubt(na,lbwa,ubwa,error=error)
  a_z = general(ubt_z,error)
  swbv_z=rc(ubt_z,error)
  rhs_z=trp(swbv_z%sw) * rhs_z
  swub_z=qr(swbv_z%bv,error)
  rhs_z=trp(swub_z%sw) * rhs_z
  call cpu_time(t0)
  x_z=solve(swub_z%ub, rhs_z, error)
  call cpu_time(t1)
  rhs_z=matmul(a_z,x_z)
  scale=maxabs(a_z)*maxabs(x_z)
  test_name = "Complex Back Subs. (n=40)"
  call z_output_result_upper(test_name,rhs0_z/scale,rhs_z/scale,ubwa,ubwa,t0,t1,c*tol,error)

  na=40; lbwa=5; ubwa=7
  rhs_v_z=z_random_vector(na)
  rhs0_v_z=rhs_v_z
  ubt_z=z_random_ubt(na,lbwa,ubwa,error=error)
  a_z = general(ubt_z,error)
  swbv_z=rc(ubt_z,error)
  rhs_v_z=trp(swbv_z%sw) * rhs_v_z
  swub_z=qr(swbv_z%bv,error)
  rhs_v_z=trp(swub_z%sw) * rhs_v_z
  call cpu_time(t0)
  x_v_z=solve(swub_z%ub, rhs_v_z, error)
  call cpu_time(t1)
  rhs_v_z=matmul(a_z,x_v_z)
  scale=maxabs(a_z)*maxabs(x_v_z)
  test_name = "Complex Vector Back Subs. (n=40)"
  call z_output_result_upper(test_name,reshape(rhs0_v_z/scale,[na,1]),&
       reshape(rhs_v_z/scale,[na,1]),ubwa,ubwa,t0,t1,c*tol,error)

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
  rhs_z=z_random_matrix(nc,na)
  rhs0_z=rhs_z
  ubt_z=z_random_ubt(na,lbwa,ubwa,error=error)
  a_z = general(ubt_z,error)
  swbv_z=rc(ubt_z,error)
  swub_z=qr(swbv_z%bv,error)
  bv_z=bv(swub_z%ub,error)
  call cpu_time(t0)
  x_z=solve(bv_z, rhs_z, error)
  call cpu_time(t1)
  x_z = x_z * trp(swub_z%sw)
  x_z = x_z * trp(swbv_z%sw)
  rhs_z=matmul(x_z,a_z)
  scale=maxabs(a_z)*maxabs(x_z)
  test_name = "Complex Forward Subs. (n=40)"
  call z_output_result_upper(test_name,rhs0_z/scale,rhs_z/scale,ubwa,ubwa,t0,t1,c*tol,error)

  na=40; lbwa=5; ubwa=7
  rhs_v_z=z_random_vector(na)
  rhs0_v_z=rhs_v_z
  ubt_z=z_random_ubt(na,lbwa,ubwa,error=error)
  a_z = general(ubt_z,error)
  swbv_z=rc(ubt_z,error)
  swub_z=qr(swbv_z%bv,error)
  bv_z=bv(swub_z%ub,error)
  call cpu_time(t0)
  x_v_z=solve(bv_z, rhs_v_z, error)
  call cpu_time(t1)
  x_v_z = x_v_z * trp(swub_z%sw)
  x_v_z = x_v_z * trp(swbv_z%sw)
  rhs_v_z=matmul(x_v_z,a_z)
  scale=maxabs(a_z)*maxabs(x_v_z)
  test_name = "Complex Vector Forward Subs. (n=40)"
  call z_output_result_upper(test_name,reshape(rhs0_v_z/scale,[1,na]),&
       reshape(rhs_v_z/scale,[1,na]),ubwa,ubwa,t0,t1,c*tol,error)
  
end program test_solve
