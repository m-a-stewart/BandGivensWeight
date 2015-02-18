program test_general_wb
  use mod_orth_rank
  use mod_test_data
  implicit none
  real(kind=dp) :: t0, t1
  type(error_info) :: error
  integer(kind=int32), parameter :: n=50, rmax=13, lbwmax=rmax+5, ubw=2, ubwmax=10
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
  complex(kind=dp), dimension(n) :: d_c
  type(d_wb), allocatable :: wb_d
  type(c_wb), allocatable :: wb_c

  call initialize_errors
  
  wb_d=d_new_wb(n,lbwmax,ubwmax)
  wb_c=c_new_wb(n,lbwmax,ubwmax)

  call random_seed
  call random_matrix(u)
  call random_matrix(v)
  call random_matrix(d)
  u0=u; v0=v
  call random_matrix(u_c)
  call random_matrix(v_c)
  call random_matrix(d_c)
  u0_c=u_c; v0_c=v_c

  ! test one
  print *
  print *, "--------------------------------"
  print *
  print *, "Real WB Decomposition Tests"
  print *
  call d_assemble_lower(a,u,v,d,ubw)
  a0=a
  call cpu_time(t0)
  call general_to_wb(a,wb_d,ubw,tol,error)
  call cpu_time(t1)
  call wb_to_general(wb_d,a1,error)
  test_name="Real WB;"
  call d_output_result_lower(test_name,a0,a1,rmax,wb_d%lbw,t0,t0,tol2,error)
  deallocate(wb_d)

  na=40
  lbwa=3; ubwa=5
  wb_d=d_random_wb(na,lbwa,ubwa)
  call wb_to_general(wb_d,a(1:na,1:na))
  a1(1:na,1:na)=a(1:na,1:na)
  call general_to_wb(a(1:na,1:na),wb_d,ubwa,tol)
  call wb_to_general(wb_d,a0(1:na,1:na))
  test_name="Random Real WB;"  
  call d_output_result_upper(test_name,a0(1:na,1:na),a1(1:na,1:na),lbwa,wb_d%lbw,t0,t0,tol2,error)
  deallocate(wb_d)

  print *
  print *, "--------------------------------"
  print *
  print *, "Complex WB Decomposition Tests"
  print *
  call c_assemble_lower(a_c,u_c,v_c,d_c,ubw)
  a0_c=a_c
  call cpu_time(t0)
  call general_to_wb(a_c,wb_c,ubw,tol,error)
  call cpu_time(t1)
  call wb_to_general(wb_c,a1_c,error)
  test_name="Complex WB;"
  call c_output_result_lower(test_name,a0_c,a1_c,rmax,wb_c%lbw,t0,t0,tol2,error)
  deallocate(wb_c)

  na=40
  lbwa=3; ubwa=5
  wb_c=c_random_wb(na,lbwa,ubwa)
  call wb_to_general(wb_c,a_c(1:na,1:na))
  a1_c(1:na,1:na)=a_c(1:na,1:na)
  call general_to_wb(a_c(1:na,1:na),wb_c,ubwa,tol)
  call wb_to_general(wb_c,a0_c(1:na,1:na))
  test_name="Random Complex WB;"  
  call c_output_result_upper(test_name,a0_c(1:na,1:na),a1_c(1:na,1:na),&
       lbwa,wb_c%lbw,t0,t0,tol2,error)
  deallocate(wb_c)

end program test_general_wb
