program test_convert_ub_and_bv
  use mod_orth_rank
  use mod_test_data
  implicit none

  !
  ! Note: This does not cover large rank cases (e.g. ubw=n-1).
  !

  real(kind=dp) :: t0, t1
  integer(kind=int32) :: na, lbwa, ubwa
  type(error_info) :: error
  real(kind=dp), parameter :: tol=1e-15, c=2.0
  !
  real(kind=dp), dimension(:,:), allocatable :: a0_d, a1_d
  complex(kind=dp), dimension(:,:), allocatable :: a0_c, a1_c

  type(d_ub), allocatable :: ub_d
  type(c_ub), allocatable :: ub_c
  type(d_bv), allocatable :: bv_d
  type(c_bv), allocatable :: bv_c

  call initialize_errors

  print *
  print *, "--------------------------------"
  print *
  print *, "Real UB and BV Conversion Tests:"
  print *
  na=40
  lbwa=3; ubwa=5
  ub_d=d_random_ub(na,lbwa,ubwa,error=error)
  a0_d=general_of(ub_d,error)
  call cpu_time(t0)
  bv_d=bv(ub_d,error)
  call cpu_time(t1)
  a1_d=general_of(bv_d,error)
  test_name = "Real UB to BV, n=40;"
  call d_output_result_upper(test_name,a0_d,a1_d,ubwa,bv_d%ubw,t0,t1,c*tol,error)

  na=1
  lbwa=0; ubwa=0
  ub_d=d_random_ub(na,lbwa,ubwa,error=error)
  a0_d=general_of(ub_d,error)
  call cpu_time(t0)
  bv_d=bv(ub_d,error)
  call cpu_time(t1)
  a1_d=general_of(bv_d,error)
  test_name = "Real UB to BV, n=1;"
  call d_output_result_upper(test_name,a0_d,a1_d,ubwa,bv_d%ubw,t0,t1,c*tol,error)

  na=2
  lbwa=1; ubwa=1
  ub_d=d_random_ub(na,lbwa,ubwa,error=error)
  a0_d=general_of(ub_d,error)
  call cpu_time(t0)
  bv_d=bv(ub_d,error)
  call cpu_time(t1)
  a1_d=general_of(bv_d,error)
  test_name = "Real UB to BV, n=2;"
  call d_output_result_upper(test_name,a0_d,a1_d,ubwa,bv_d%ubw,t0,t1,c*tol,error)

  na=3
  lbwa=1; ubwa=1
  ub_d=d_random_ub(na,lbwa,ubwa,error=error)
  a0_d=general_of(ub_d,error)
  call cpu_time(t0)
  bv_d=bv(ub_d,error)
  call cpu_time(t1)
  a1_d=general_of(bv_d,error)
  test_name = "Real UB to BV, n=3;"
  call d_output_result_upper(test_name,a0_d,a1_d,ubwa,bv_d%ubw,t0,t1,c*tol,error)

  na=4
  lbwa=2; ubwa=2
  ub_d=d_random_ub(na,lbwa,ubwa,error=error)
  a0_d=general_of(ub_d,error)
  call cpu_time(t0)
  bv_d=bv(ub_d,error)
  call cpu_time(t1)
  a1_d=general_of(bv_d,error)
  test_name = "Real UB to BV, n=4;"
  call d_output_result_upper(test_name,a0_d,a1_d,ubwa,bv_d%ubw,t0,t1,c*tol,error)
  print *

  na=40
  lbwa=3; ubwa=5
  bv_d=d_random_bv(na,lbwa,ubwa,error=error)
  a0_d=general_of(bv_d,error)
  call cpu_time(t0)
  ub_d=ub(bv_d,error)
  call cpu_time(t1)
  a1_d=general_of(ub_d,error)
  test_name = "Real BV to UB, n=40;"
  call d_output_result_upper(test_name,a0_d,a1_d,ubwa,ub_d%ubw,t0,t1,c*tol,error)

  na=1
  lbwa=0; ubwa=0
  bv_d=d_random_bv(na,lbwa,ubwa,error=error)
  a0_d=general_of(bv_d,error)
  call cpu_time(t0)
  ub_d=ub(bv_d,error)
  call cpu_time(t1)
  a1_d=general_of(ub_d,error)
  test_name = "Real BV to UB, n=1;"
  call d_output_result_upper(test_name,a0_d,a1_d,ubwa,ub_d%ubw,t0,t1,c*tol,error)

  na=2
  lbwa=1; ubwa=1
  bv_d=d_random_bv(na,lbwa,ubwa,error=error)
  a0_d=general_of(bv_d,error)
  call cpu_time(t0)
  ub_d=ub(bv_d,error)
  call cpu_time(t1)
  a1_d=general_of(ub_d,error)
  test_name = "Real BV to UB, n=2;"
  call d_output_result_upper(test_name,a0_d,a1_d,ubwa,ub_d%ubw,t0,t1,c*tol,error)

  na=3
  lbwa=1; ubwa=1
  bv_d=d_random_bv(na,lbwa,ubwa,error=error)
  a0_d=general_of(bv_d,error)
  call cpu_time(t0)
  ub_d=ub(bv_d,error)
  call cpu_time(t1)
  a1_d=general_of(ub_d,error)
  test_name = "Real BV to UB, n=3;"
  call d_output_result_upper(test_name,a0_d,a1_d,ubwa,ub_d%ubw,t0,t1,c*tol,error)

  na=4
  lbwa=2; ubwa=2
  bv_d=d_random_bv(na,lbwa,ubwa,error=error)
  a0_d=general_of(bv_d,error)
  call cpu_time(t0)
  ub_d=ub(bv_d,error)
  call cpu_time(t1)
  a1_d=general_of(ub_d,error)
  test_name = "Real BV to UB, n=4;"
  call d_output_result_upper(test_name,a0_d,a1_d,ubwa,ub_d%ubw,t0,t1,c*tol,error)
  

  
  print *
  print *, "--------------------------------"
  print *
  print *, "Complex UB and BV Conversion Tests:"
  print *

  na=40
  lbwa=3; ubwa=5
  ub_c=c_random_ub(na,lbwa,ubwa,error=error)
  a0_c=general_of(ub_c,error)
  call cpu_time(t0)
  bv_c=bv(ub_c,error)
  call cpu_time(t1)
  a1_c=general_of(bv_c,error)
  test_name = "Complex UB to BV, n=40;"
  call c_output_result_upper(test_name,a0_c,a1_c,ubwa,bv_c%ubw,t0,t1,c*tol,error)

  na=1
  lbwa=0; ubwa=0
  ub_c=c_random_ub(na,lbwa,ubwa,error=error)
  a0_c=general_of(ub_c,error)
  call cpu_time(t0)
  bv_c=bv(ub_c,error)
  call cpu_time(t1)
  a1_c=general_of(bv_c,error)
  test_name = "Complex UB to BV, n=1;"
  call c_output_result_upper(test_name,a0_c,a1_c,ubwa,bv_c%ubw,t0,t1,c*tol,error)

  na=2
  lbwa=1; ubwa=1
  ub_c=c_random_ub(na,lbwa,ubwa,error=error)
  a0_c=general_of(ub_c,error)
  call cpu_time(t0)
  bv_c=bv(ub_c,error)
  call cpu_time(t1)
  a1_c=general_of(bv_c,error)
  test_name = "Complex UB to BV, n=2;"
  call c_output_result_upper(test_name,a0_c,a1_c,ubwa,bv_c%ubw,t0,t1,c*tol,error)

  na=3
  lbwa=1; ubwa=1
  ub_c=c_random_ub(na,lbwa,ubwa,error=error)
  a0_c=general_of(ub_c,error)
  call cpu_time(t0)
  bv_c=bv(ub_c,error)
  call cpu_time(t1)
  a1_c=general_of(bv_c,error)
  test_name = "Complex UB to BV, n=3;"
  call c_output_result_upper(test_name,a0_c,a1_c,ubwa,bv_c%ubw,t0,t1,c*tol,error)

  na=4
  lbwa=2; ubwa=2
  ub_c=c_random_ub(na,lbwa,ubwa,error=error)
  a0_c=general_of(ub_c,error)
  call cpu_time(t0)
  bv_c=bv(ub_c,error)
  call cpu_time(t1)
  a1_c=general_of(bv_c,error)
  test_name = "Complex UB to BV, n=4;"
  call c_output_result_upper(test_name,a0_c,a1_c,ubwa,bv_c%ubw,t0,t1,c*tol,error)
  print *

  na=40
  lbwa=3; ubwa=5
  bv_c=c_random_bv(na,lbwa,ubwa,error=error)
  a0_c=general_of(bv_c,error)
  call cpu_time(t0)
  ub_c=ub(bv_c,error)
  call cpu_time(t1)
  a1_c=general_of(ub_c,error)
  test_name = "Complex BV to UB, n=40;"
  call c_output_result_upper(test_name,a0_c,a1_c,ubwa,ub_c%ubw,t0,t1,c*tol,error)

  na=1
  lbwa=0; ubwa=0
  bv_c=c_random_bv(na,lbwa,ubwa,error=error)
  a0_c=general_of(bv_c,error)
  call cpu_time(t0)
  ub_c=ub(bv_c,error)
  call cpu_time(t1)
  a1_c=general_of(ub_c,error)
  test_name = "Complex BV to UB, n=1;"
  call c_output_result_upper(test_name,a0_c,a1_c,ubwa,ub_c%ubw,t0,t1,c*tol,error)

  na=2
  lbwa=1; ubwa=1
  bv_c=c_random_bv(na,lbwa,ubwa,error=error)
  a0_c=general_of(bv_c,error)
  call cpu_time(t0)
  ub_c=ub(bv_c,error)
  call cpu_time(t1)
  a1_c=general_of(ub_c,error)
  test_name = "Complex BV to UB, n=2;"
  call c_output_result_upper(test_name,a0_c,a1_c,ubwa,ub_c%ubw,t0,t1,c*tol,error)

  na=3
  lbwa=1; ubwa=1
  bv_c=c_random_bv(na,lbwa,ubwa,error=error)
  a0_c=general_of(bv_c,error)
  call cpu_time(t0)
  ub_c=ub(bv_c,error)
  call cpu_time(t1)
  a1_c=general_of(ub_c,error)
  test_name = "Complex BV to UB, n=3;"
  call c_output_result_upper(test_name,a0_c,a1_c,ubwa,ub_c%ubw,t0,t1,c*tol,error)

  na=4
  lbwa=2; ubwa=2
  bv_c=c_random_bv(na,lbwa,ubwa,error=error)
  a0_c=general_of(bv_c,error)
  call cpu_time(t0)
  ub_c=ub(bv_c,error)
  call cpu_time(t1)
  a1_c=general_of(ub_c,error)
  test_name = "Complex BV to UB, n=4;"
  call c_output_result_upper(test_name,a0_c,a1_c,ubwa,ub_c%ubw,t0,t1,c*tol,error)

end program test_convert_ub_and_bv
