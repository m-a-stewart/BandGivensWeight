program test_general_wb
  use mod_orth_rank
  use mod_test_data
  implicit none
  real(kind=dp) :: t0, t1
  type(error_info) :: error
  integer, parameter :: n=50, rmax=13, lbwmax=rmax+5, ubw=2, ubwmax=10
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
  type(d_wb) :: wb_d
  type(c_wb) :: wb_c
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
  call lower_to_wb(a,wb_d,ubw,tol,error)
  call cpu_time(t1)
  call wb_to_lower(wb_d,a1,error)
  test_name="Real WB;"
  call d_output_result_lower(test_name,a0,a1,rmax,wb_d%lbw,t0,t0,tol2,error)


  print *
  print *, "--------------------------------"
  print *
  print *, "Complex WB Decomposition Tests"
  print *
  call c_assemble_lower(a_c,u_c,v_c,d_c,ubw)
  a0_c=a_c
  call cpu_time(t0)
  call lower_to_wb(a_c,wb_c,ubw,tol,error)
  call cpu_time(t1)
  call wb_to_lower(wb_c,a1_c,error)
  test_name="Complex WB;"
  call c_output_result_lower(test_name,a0_c,a1_c,rmax,wb_c%lbw,t0,t0,tol2,error)

end program test_general_wb
