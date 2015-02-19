program test_general
  use mod_orth_rank
  use mod_test_data
  implicit none
  real(kind=dp) :: t0, t1
  type(error_info) :: error
  integer(kind=int32), parameter :: n=50, rmax=13, ubwmax=rmax+5, lbw=2, lbwmax=10
  integer(kind=int32) :: na, lbwa, ubwa
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
  type(d_ub), allocatable :: ub_d
  type(c_ub), allocatable :: ub_c
  type(d_bv), allocatable :: bv_d
  type(c_bv), allocatable :: bv_c
  real(kind=dp), dimension(:,:), allocatable :: a_a
  ub_d=d_new_ub(n,lbwmax,ubwmax)
  ub_c=c_new_ub(n,lbwmax,ubwmax)
  bv_d=d_new_bv(n,lbwmax,ubwmax)
  bv_c=c_new_bv(n,lbwmax,ubwmax)

  call initialize_errors

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
  call general_to_ub(a,ub_d,lbw,tol,error)
  call cpu_time(t1)
  call ub_to_general(ub_d,a1,error)
  test_name="Real UB;"
  call d_output_result_upper(test_name,a0,a1,rmax,ub_d%ubw,t0,t0,tol2,error)
  ! square termination.
  u=u0; v=v0;
  u(1:n-rmax-1,rmax)=0.0_dp
  call d_assemble_upper(a,u,v,d,lbw)
  a0=a
  call cpu_time(t0)
  call general_to_ub(a,ub_d,lbw,tol,error)
  call cpu_time(t1)
  call ub_to_general(ub_d,a1,error)
  test_name="Real Sq. term. UB;"
  call d_output_result_upper(test_name,a0,a1,rmax,ub_d%ubw,t0,t0,tol2,error)

  na=40
  lbwa=3; ubwa=5
  ub_d=d_random_ub(na,lbwa,ubwa)
  a(1:na,1:na)=general_of(ub_d)
  a1(1:na,1:na)=a(1:na,1:na)
  ub_d=ub_of_general(a(1:na,1:na),lbwa,lbwmax,ubwmax,tol)
  a0(1:na,1:na) = general_of(ub_d)
  test_name="Random Real UB;"  
  call d_output_result_upper(test_name,a0(1:na,1:na),a1(1:na,1:na),ubwa,ub_d%ubw,t0,t0,tol2,error)
  deallocate(ub_d)

  print *
  print *, "--------------------------------"
  print *
  print *, "Real BV Decomposition Tests"
  print *
  u=u0; v=v0
  call d_assemble_upper(a,u,v,d,lbw)
  a0=a
  call cpu_time(t0)
  call general_to_bv(a,bv_d,lbw,tol,error)
  call cpu_time(t1)
  call bv_to_general(bv_d,a1,error)
  test_name="Real BV;"
  call d_output_result_upper(test_name,a0,a1,rmax,bv_d%ubw,t0,t0,tol2,error)
  ! square termination BV
  u=u0; v=v0
  v(rmax,rmax+2:n)=0.0_dp
  call d_assemble_upper(a,u,v,d,lbw)
  a0=a
  call cpu_time(t0)
  call general_to_bv(a,bv_d,lbw,tol,error)
  call cpu_time(t1)
  call bv_to_general(bv_d,a1,error)
  test_name="Real Sq. Term BV;"
  call d_output_result_upper(test_name,a0,a1,rmax,bv_d%ubw,t0,t0,tol2,error)
  !
  na=40
  lbwa=3; ubwa=5
  bv_d=d_random_bv(na,lbwa,ubwa)
  call bv_to_general(bv_d,a(1:na,1:na))
  a1(1:na,1:na)=a(1:na,1:na)
  call general_to_bv(a(1:na,1:na),bv_d,lbwa,tol)
  call bv_to_general(bv_d,a0(1:na,1:na))
  test_name="Random Real BV;"  
  call d_output_result_upper(test_name,a0(1:na,1:na),a1(1:na,1:na),ubwa,bv_d%ubw,t0,t0,tol2,error)
  deallocate(bv_d)

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
  call general_to_ub(a_c,ub_c,lbw,tol,error)
  call cpu_time(t1)
  
  call ub_to_general(ub_c,a1_c,error)
  test_name="Complex UB;"
  call c_output_result_upper(test_name,a0_c,a1_c,rmax,ub_c%ubw,t0,t0,tol2,error)
  ! Complex UB square termination test
  u_c=u0_c; v_c=v0_c
  u_c(1:n-rmax-1,rmax)=0.0_dp
  call c_assemble_upper(a_c,u_c,v_c,dc,lbw)
  a0_c=a_c
  call cpu_time(t0)
  call general_to_ub(a_c,ub_c,lbw,tol,error)
  call cpu_time(t1)
  call ub_to_general(ub_c,a1_c,error)
  test_name="Complex square term. UB;"
  call c_output_result_upper(test_name,a0_c,a1_c,rmax,ub_c%ubw,t0,t0,tol2,error)

  na=40
  lbwa=3; ubwa=5
  ub_c=c_random_ub(na,lbwa,ubwa)
  a_c(1:na,1:na)=general_of(ub_c)
  a1_c(1:na,1:na)=a_c(1:na,1:na)
  ub_c=ub_of_general(a_c(1:na,1:na), lbwa, lbwmax, ubwmax, tol)
  a0_c(1:na,1:na)=general_of(ub_c)
  test_name="Random Complex UB;"  
  call c_output_result_upper(test_name,a0_c(1:na,1:na),a1_c(1:na,1:na),ubwa,ub_c%ubw,t0,t0,tol2,error)
  deallocate(ub_c)
  print *
  
  print *
  print *, "--------------------------------"
  print *
  print *, "Complex BV Decomposition Tests"
  print *
  u_c=u0_c; v_c=v0_c
  call c_assemble_upper(a_c,u_c,v_c,dc,lbw)
  a0_c=a_c
  call cpu_time(t0)
  call general_to_bv(a_c,bv_c,lbw,tol,error)
  call cpu_time(t1)
  call bv_to_general(bv_c,a1_c,error)
  test_name="Complex BV;"
  call c_output_result_upper(test_name,a0_c,a1_c,rmax,bv_c%ubw,t0,t0,tol2,error)
  !
  u_c=u0_c; v_c=v0_c
  v_c(rmax,rmax+2:n)=0.0_dp
  call c_assemble_upper(a_c,u_c,v_c,dc,lbw)
  a0_c=a_c
  call cpu_time(t0)
  call general_to_bv(a_c,bv_c,lbw,tol,error)
  call cpu_time(t1)
  call bv_to_general(bv_c,a1_c,error)
  test_name="Complex Sq. Term. BV;"
  call c_output_result_upper(test_name,a0_c,a1_c,rmax,bv_c%ubw,t0,t0,tol2,error)

  na=40
  lbwa=3; ubwa=5
  bv_c=c_random_bv(na,lbwa,ubwa)
  call bv_to_general(bv_c,a_c(1:na,1:na))
  a1_c(1:na,1:na)=a_c(1:na,1:na)
  call general_to_bv(a_c(1:na,1:na),bv_c,lbwa,tol)
  call bv_to_general(bv_c,a0_c(1:na,1:na))
  test_name="Random Complex BV;"  
  call c_output_result_upper(test_name,a0_c(1:na,1:na),a1_c(1:na,1:na),ubwa,bv_c%ubw,t0,t0,tol2,error)
  deallocate(bv_c)
  print *

end program test_general
