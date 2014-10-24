program test_general
  use nested
  use test_data
  implicit none
  real(kind=dp) :: t0, t1
  type(error_info) :: error
  integer, parameter :: n=50, rmax=13, ubwmax=rmax+5, lbw=2, lbwmax=10
  real(kind=dp), parameter :: tol=1e-14, tol1=1e-14, tol2=1e-10
  !
  real(kind=dp), dimension(n,n) :: a, a0, a1
  real(kind=dp), dimension(n,rmax) :: u, u0
  real(kind=dp), dimension(rmax,n) :: v, v0
  real(kind=dp), dimension(n) :: d
  complex(kind=dp), dimension(n,n) :: a_c, a0_c, a1_c
  complex(kind=dp), dimension(n,rmax) :: u_c, u0_c
  complex(kind=dp), dimension(rmax,n) :: v_c, v0_c
  complex(kind=dp), dimension(n) :: dc
  type(d_ub) :: ub_d
  type(c_ub) :: ub_c
  type(d_bv) :: bv_d
  type(c_bv) :: bv_c
  ub_d=d_new_ub(n,lbwmax,ubwmax)
  ub_c=c_new_ub(n,lbwmax,ubwmax)
  bv_d=d_new_bv(n,lbwmax,ubwmax)
  bv_c=c_new_bv(n,lbwmax,ubwmax)

  call random_seed
  call random_matrix(u)
  call random_matrix(v)
  call random_matrix(d)
  u0=u; v0=v
  call random_matrix(u_c)
  call random_matrix(v_c)
  call random_matrix(dc)
  u0_c=u_c; v0_c=v_c

  ! test one
  print *
  print *, "--------------------------------"
  print *
  print *, "Real UB Decomposition Tests"
  print *
  call d_assemble_upper(a,u,v,d,lbw)
  a0=a
  call cpu_time(t0)
  call upper_to_ub(a,ub_d,lbw,tol,error)
  call cpu_time(t1)
  call ub_to_upper(ub_d,a1,error)
  test_name="Real UB;"
  call d_output_result_upper(test_name,a0,a1,rmax,ub_d%ubw,t0,t0,tol2,error)
  ! square termination.
  u=u0; v=v0;
  u(1:n-rmax-1,rmax)=0.0_dp
  call d_assemble_upper(a,u,v,d,lbw)
  a0=a
  call cpu_time(t0)
  call upper_to_ub(a,ub_d,lbw,tol,error)
  call cpu_time(t1)
  call ub_to_upper(ub_d,a1,error)
  test_name="Real Sq. term. UB;"
  call d_output_result_upper(test_name,a0,a1,rmax,ub_d%ubw,t0,t0,tol2,error)
  print *
  print *, "--------------------------------"
  print *
  print *, "Real BV Decomposition Tests"
  print *
  u=u0; v=v0
  call d_assemble_upper(a,u,v,d,lbw)
  a0=a
  call cpu_time(t0)
  call upper_to_bv(a,bv_d,lbw,tol,error)
  call cpu_time(t1)
  call bv_to_upper(bv_d,a1,error)
  test_name="Real BV;"
  call d_output_result_upper(test_name,a0,a1,rmax,bv_d%ubw,t0,t0,tol2,error)
  ! square termination BV
  u=u0; v=v0
  v(rmax,rmax+2:n)=0.0_dp
  call d_assemble_upper(a,u,v,d,lbw)
  a0=a
  call cpu_time(t0)
  call upper_to_bv(a,bv_d,lbw,tol,error)
  call cpu_time(t1)
  call bv_to_upper(bv_d,a1,error)
  test_name="Real Sq. Term BV;"
  call d_output_result_upper(test_name,a0,a1,rmax,bv_d%ubw,t0,t0,tol2,error)
  !
  ! Complex UB test
  !
  print *
  print *, "--------------------------------"
  print *
  print *, "Complex UB Decomposition Tests"
  print *
  u_c=u0_c; v_c=v0_c
  call c_assemble_upper(a_c,u_c,v_c,dc,lbw)
  a0_c=a_c
  call cpu_time(t0)
  call upper_to_ub(a_c,ub_c,lbw,tol,error)
  call cpu_time(t1)
  call ub_to_upper(ub_c,a1_c,error)
  test_name="Complex UB;"
  call c_output_result_upper(test_name,a0_c,a1_c,rmax,ub_c%ubw,t0,t0,tol2,error)
  ! Complex UB square termination test
  u_c=u0_c; v_c=v0_c
  u_c(1:n-rmax-1,rmax)=0.0_dp
  call c_assemble_upper(a_c,u_c,v_c,dc,lbw)
  a0_c=a_c
  call cpu_time(t0)
  call upper_to_ub(a_c,ub_c,lbw,tol,error)
  call cpu_time(t1)
  call ub_to_upper(ub_c,a1_c,error)
  test_name="Complex square term. UB;"
  call c_output_result_upper(test_name,a0_c,a1_c,rmax,ub_c%ubw,t0,t0,tol2,error)
  print *
  print *, "--------------------------------"
  print *
  print *, "Complex BV Decomposition Tests"
  print *
  u_c=u0_c; v_c=v0_c
  call c_assemble_upper(a_c,u_c,v_c,dc,lbw)
  a0_c=a_c
  call cpu_time(t0)
  call upper_to_bv(a_c,bv_c,lbw,tol,error)
  call cpu_time(t1)
  call bv_to_upper(bv_c,a1_c,error)
  test_name="Complex BV;"
  call c_output_result_upper(test_name,a0_c,a1_c,rmax,bv_c%ubw,t0,t0,tol2,error)
  !
  u_c=u0_c; v_c=v0_c
  v_c(rmax,rmax+2:n)=0.0_dp
  call c_assemble_upper(a_c,u_c,v_c,dc,lbw)
  a0_c=a_c
  call cpu_time(t0)
  call upper_to_bv(a_c,bv_c,lbw,tol,error)
  call cpu_time(t1)
  call bv_to_upper(bv_c,a1_c,error)
  test_name="Complex Sq. Term. BV;"
  call c_output_result_upper(test_name,a0_c,a1_c,rmax,bv_c%ubw,t0,t0,tol2,error)
  print *

end program test_general
