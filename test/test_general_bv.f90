program test_general_bv
  use mod_orth_rank
  use mod_test_data
  implicit none
  real(kind=dp) :: t0, t1
  type(error_info) :: error
  integer(kind=int32) :: na, lbwa, ubwa, j
  real(kind=dp), parameter :: tol=1e-14, c=1.5
  !
  real(kind=dp), dimension(:,:), allocatable :: a_d, a0_d, a1_d
  complex(kind=dp), dimension(:,:), allocatable :: a_c, a0_c, a1_c
  type(d_bv), allocatable :: bv_d
  type(c_bv), allocatable :: bv_c

  call initialize_errors

  print *
  print *, "--------------------------------"
  print *
  print *, "Real BV Decomposition Tests"
  print *
  
  na=40
  lbwa=3; ubwa=5
  bv_d=d_random_bv(na,lbwa,ubwa,error=error)
  a_d=general(bv_d,error)
  a1_d=a_d
  call cpu_time(t0)
  bv_d=bv(a_d,lbwa,lbwa,ubwa+1,tol,error)
  call cpu_time(t1)
  a0_d=general(bv_d,error)
  test_name="Random Real BV;"  
  call d_output_result_upper(test_name,a0_d,a1_d,ubwa,bv_d%ubw,t0,t1,c*tol,error)

  na=40
  lbwa=3; ubwa=0
  bv_d=d_random_bv(na,lbwa,ubwa,error=error)
  a_d=general(bv_d,error)
  a1_d=a_d
  call cpu_time(t0)
  bv_d=bv(a_d,lbwa,lbwa,ubwa+1,tol,error)
  call cpu_time(t1)
  a0_d=general(bv_d,error)
  test_name="Random Real BV, ubwa=0;"  
  call d_output_result_upper(test_name,a0_d,a1_d,ubwa,bv_d%ubw,t0,t1,c*tol,error)

  na=40
  lbwa=3; ubwa=1
  bv_d=d_random_bv(na,lbwa,ubwa,error=error)
  a_d=general(bv_d,error)
  a1_d=a_d
  call cpu_time(t0)
  bv_d=bv(a_d,lbwa,lbwa,ubwa+1,tol,error)
  call cpu_time(t1)
  a0_d=general(bv_d,error)
  test_name="Random Real BV, ubwa=1;"
  call d_output_result_upper(test_name,a0_d,a1_d,ubwa,bv_d%ubw,t0,t1,c*tol,error)

  na=1
  lbwa=0; ubwa=0
  bv_d=d_random_bv(na,lbwa,ubwa,error=error)
  a_d=general(bv_d,error)
  a1_d=a_d
  call cpu_time(t0)
  bv_d=bv(a_d,lbwa,lbwa,ubwa+1,tol,error)
  call cpu_time(t1)
  a0_d=general(bv_d,error)
  test_name="Random Real BV, na=1;"  
  call d_output_result_upper(test_name,a0_d,a1_d,ubwa,bv_d%ubw,t0,t1,c*tol,error)

  na=2
  lbwa=1; ubwa=1
  bv_d=d_random_bv(na,lbwa,ubwa,error=error)
  a_d=general(bv_d,error)
  a1_d=a_d
  call cpu_time(t0)
  bv_d=bv(a_d,lbwa,lbwa,ubwa+1,tol,error)
  call cpu_time(t1)
  a0_d=general(bv_d,error)
  test_name="Random Real BV, na=2;"
  call d_output_result_upper(test_name,a0_d,a1_d,ubwa,bv_d%ubw,t0,t1,c*tol,error)

  na=3
  lbwa=1; ubwa=1
  bv_d=d_random_bv(na,lbwa,ubwa,error=error)
  a_d=general(bv_d,error)
  a1_d=a_d
  call cpu_time(t0)
  bv_d=bv(a_d,lbwa,lbwa,ubwa+1,tol,error)
  call cpu_time(t1)
  a0_d=general(bv_d,error)
  test_name="Random Real BV, na=3;"
  call d_output_result_upper(test_name,a0_d,a1_d,ubwa,bv_d%ubw,t0,t1,c*tol,error)

  na=4
  lbwa=2; ubwa=2
  bv_d=d_random_bv(na,lbwa,ubwa,error=error)
  a_d=general(bv_d,error)
  a1_d=a_d
  call cpu_time(t0)
  bv_d=bv(a_d,lbwa,lbwa,ubwa+1,tol,error)
  call cpu_time(t1)
  a0_d=general(bv_d,error)
  test_name="Random Real BV, na=4;"
  call d_output_result_upper(test_name,a0_d,a1_d,ubwa,bv_d%ubw,t0,t1,c*tol,error)
  
  na=50
  lbwa=3; ubwa=13
  bv_d=d_random_bv(na,[ (lbwa, j=1,na) ], &
       [ (ubwa, j=1,ubwa-1), (ubwa-1, j=ubwa,na) ], error=error )
  a_d=general(bv_d, error)
  a1_d=a_d
  call cpu_time(t0)
  bv_d=bv(a_d,lbwa,lbwa,ubwa+1,tol, error)
  call cpu_time(t1)
  a0_d=general(bv_d, error)
  test_name="Random Real Square Termination BV;"  
  call d_output_result_upper(test_name,a0_d,a1_d,ubwa,bv_d%ubw,t0,t1,c*tol,error)

  print *
  print *, "--------------------------------"
  print *
  print *, "Complex BV Decomposition Tests"
  print *
  na=40
  lbwa=3; ubwa=5
  bv_c=c_random_bv(na,lbwa,ubwa,error=error)
  a_c=general(bv_c,error)
  a1_c=a_c
  call cpu_time(t0)
  bv_c=bv(a_c,lbwa,lbwa,ubwa+1,tol,error)
  call cpu_time(t1)
  a0_c=general(bv_c)
  test_name="Random Complex BV;"  
  call c_output_result_upper(test_name,a0_c,a1_c,ubwa,bv_c%ubw,t0,t1,c*tol,error)

  na=40
  lbwa=3; ubwa=0
  bv_c=c_random_bv(na,lbwa,ubwa,error=error)
  a_c=general(bv_c,error)
  a1_c=a_c
  call cpu_time(t0)
  bv_c=bv(a_c,lbwa,lbwa,ubwa+1,tol,error)
  call cpu_time(t1)
  a0_c=general(bv_c)
  test_name="Random Complex BV, ubwa=0;"  
  call c_output_result_upper(test_name,a0_c,a1_c,ubwa,bv_c%ubw,t0,t1,c*tol,error)

  na=40
  lbwa=3; ubwa=1
  bv_c=c_random_bv(na,lbwa,ubwa,error=error)
  a_c=general(bv_c,error)
  a1_c=a_c
  call cpu_time(t0)
  bv_c=bv(a_c,lbwa,lbwa,ubwa+1,tol,error)
  call cpu_time(t1)
  a0_c=general(bv_c)
  test_name="Random Complex BV, ubwa=1;"
  call c_output_result_upper(test_name,a0_c,a1_c,ubwa,bv_c%ubw,t0,t1,c*tol,error)

  na=1
  lbwa=0; ubwa=0
  bv_c=c_random_bv(na,lbwa,ubwa,error=error)
  a_c=general(bv_c,error)
  a1_c=a_c
  call cpu_time(t0)
  bv_c=bv(a_c,lbwa,lbwa,ubwa+1,tol,error)
  call cpu_time(t1)
  a0_c=general(bv_c,error)
  test_name="Random Complex BV, na=1;"  
  call c_output_result_upper(test_name,a0_c,a1_c,ubwa,bv_c%ubw,t0,t1,c*tol,error)

  na=2
  lbwa=1; ubwa=1
  bv_c=c_random_bv(na,lbwa,ubwa,error=error)
  a_c=general(bv_c,error)
  a1_c=a_c
  call cpu_time(t0)
  bv_c=bv(a_c,lbwa,lbwa,ubwa+1,tol,error)
  call cpu_time(t1)
  a0_c=general(bv_c,error)
  test_name="Random Complex BV, na=2;"
  call c_output_result_upper(test_name,a0_c,a1_c,ubwa,bv_c%ubw,t0,t1,c*tol,error)

  na=3
  lbwa=1; ubwa=1
  bv_c=c_random_bv(na,lbwa,ubwa,error=error)
  a_c=general(bv_c,error)
  a1_c=a_c
  call cpu_time(t0)
  bv_c=bv(a_c,lbwa,lbwa,ubwa+1,tol,error)
  call cpu_time(t1)
  a0_c=general(bv_c,error)
  test_name="Random Complex BV, na=3;"
  call c_output_result_upper(test_name,a0_c,a1_c,ubwa,bv_c%ubw,t0,t1,c*tol,error)

  na=4
  lbwa=2; ubwa=2
  bv_c=c_random_bv(na,lbwa,ubwa,error=error)
  a_c=general(bv_c,error)
  a1_c=a_c
  call cpu_time(t0)
  bv_c=bv(a_c,lbwa,lbwa,ubwa+1,tol,error)
  call cpu_time(t1)
  a0_c=general(bv_c,error)
  test_name="Random Complex BV, na=4;"
  call c_output_result_upper(test_name,a0_c,a1_c,ubwa,bv_c%ubw,t0,t1,c*tol,error)
  
  na=50
  lbwa=3; ubwa=13
  bv_c=c_random_bv(na,[ (lbwa, j=1,na) ], &
       [ (ubwa, j=1,ubwa-1), (ubwa-1, j=ubwa,na) ], error=error )
  a_c=general(bv_c, error)
  a1_c=a_c
  call cpu_time(t0)
  bv_c=bv(a_c,lbwa,lbwa,ubwa+1,tol,error)
  call cpu_time(t1)
  a0_c=general(bv_c,error)
  test_name="Random Complex Square Termination BV;"  
  call c_output_result_upper(test_name,a0_c,a1_c,ubwa,bv_c%ubw,t0,t1,c*tol,error)

end program test_general_bv
