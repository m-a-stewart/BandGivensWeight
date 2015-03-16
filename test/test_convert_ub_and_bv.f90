program test_convert_ub_and_bv
  use mod_orb
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
  complex(kind=dp), dimension(:,:), allocatable :: a0_z, a1_z

  type(d_ub), allocatable :: ub_d
  type(z_ub), allocatable :: ub_z
  type(d_bv), allocatable :: bv_d
  type(z_bv), allocatable :: bv_z

  call initialize_errors

  print *
  print *, "--------------------------------"
  print *
  print *, "Real UB and BV Conversion Tests:"
  print *
  na=40
  lbwa=3; ubwa=5
  ub_d=d_random_ub(na,lbwa,ubwa,error=error)
  a0_d=general(ub_d,error)
  call cpu_time(t0)
  bv_d=bv(ub_d,error)
  call cpu_time(t1)
  a1_d=general(bv_d,error)
  test_name = "Real UB to BV, n=40;"
  call d_output_result_upper(test_name,a0_d,a1_d,ubwa,bv_d%ubw,t0,t1,c*tol,error)

  na=1
  lbwa=0; ubwa=0
  ub_d=d_random_ub(na,lbwa,ubwa,error=error)
  a0_d=general(ub_d,error)
  call cpu_time(t0)
  bv_d=bv(ub_d,error)
  call cpu_time(t1)
  a1_d=general(bv_d,error)
  test_name = "Real UB to BV, n=1;"
  call d_output_result_upper(test_name,a0_d,a1_d,ubwa,bv_d%ubw,t0,t1,c*tol,error)

  na=2
  lbwa=1; ubwa=1
  ub_d=d_random_ub(na,lbwa,ubwa,error=error)
  a0_d=general(ub_d,error)
  call cpu_time(t0)
  bv_d=bv(ub_d,error)
  call cpu_time(t1)
  a1_d=general(bv_d,error)
  test_name = "Real UB to BV, n=2;"
  call d_output_result_upper(test_name,a0_d,a1_d,ubwa,bv_d%ubw,t0,t1,c*tol,error)

  na=3
  lbwa=1; ubwa=1
  ub_d=d_random_ub(na,lbwa,ubwa,error=error)
  a0_d=general(ub_d,error)
  call cpu_time(t0)
  bv_d=bv(ub_d,error)
  call cpu_time(t1)
  a1_d=general(bv_d,error)
  test_name = "Real UB to BV, n=3;"
  call d_output_result_upper(test_name,a0_d,a1_d,ubwa,bv_d%ubw,t0,t1,c*tol,error)

  na=4
  lbwa=2; ubwa=2
  ub_d=d_random_ub(na,lbwa,ubwa,error=error)
  a0_d=general(ub_d,error)
  call cpu_time(t0)
  bv_d=bv(ub_d,error)
  call cpu_time(t1)
  a1_d=general(bv_d,error)
  test_name = "Real UB to BV, n=4;"
  call d_output_result_upper(test_name,a0_d,a1_d,ubwa,bv_d%ubw,t0,t1,c*tol,error)
  print *

  na=40
  lbwa=3; ubwa=5
  bv_d=d_random_bv(na,lbwa,ubwa,error=error)
  a0_d=general(bv_d,error)
  call cpu_time(t0)
  ub_d=ub(bv_d,error)
  call cpu_time(t1)
  a1_d=general(ub_d,error)
  test_name = "Real BV to UB, n=40;"
  call d_output_result_upper(test_name,a0_d,a1_d,ubwa,ub_d%ubw,t0,t1,c*tol,error)

  na=1
  lbwa=0; ubwa=0
  bv_d=d_random_bv(na,lbwa,ubwa,error=error)
  a0_d=general(bv_d,error)
  call cpu_time(t0)
  ub_d=ub(bv_d,error)
  call cpu_time(t1)
  a1_d=general(ub_d,error)
  test_name = "Real BV to UB, n=1;"
  call d_output_result_upper(test_name,a0_d,a1_d,ubwa,ub_d%ubw,t0,t1,c*tol,error)

  na=2
  lbwa=1; ubwa=1
  bv_d=d_random_bv(na,lbwa,ubwa,error=error)
  a0_d=general(bv_d,error)
  call cpu_time(t0)
  ub_d=ub(bv_d,error)
  call cpu_time(t1)
  a1_d=general(ub_d,error)
  test_name = "Real BV to UB, n=2;"
  call d_output_result_upper(test_name,a0_d,a1_d,ubwa,ub_d%ubw,t0,t1,c*tol,error)

  na=3
  lbwa=1; ubwa=1
  bv_d=d_random_bv(na,lbwa,ubwa,error=error)
  a0_d=general(bv_d,error)
  call cpu_time(t0)
  ub_d=ub(bv_d,error)
  call cpu_time(t1)
  a1_d=general(ub_d,error)
  test_name = "Real BV to UB, n=3;"
  call d_output_result_upper(test_name,a0_d,a1_d,ubwa,ub_d%ubw,t0,t1,c*tol,error)

  na=4
  lbwa=2; ubwa=2
  bv_d=d_random_bv(na,lbwa,ubwa,error=error)
  a0_d=general(bv_d,error)
  call cpu_time(t0)
  ub_d=ub(bv_d,error)
  call cpu_time(t1)
  a1_d=general(ub_d,error)
  test_name = "Real BV to UB, n=4;"
  call d_output_result_upper(test_name,a0_d,a1_d,ubwa,ub_d%ubw,t0,t1,c*tol,error)
  

  
  print *
  print *, "--------------------------------"
  print *
  print *, "Complex UB and BV Conversion Tests:"
  print *

  na=40
  lbwa=3; ubwa=5
  ub_z=z_random_ub(na,lbwa,ubwa,error=error)
  a0_z=general(ub_z,error)
  call cpu_time(t0)
  bv_z=bv(ub_z,error)
  call cpu_time(t1)
  a1_z=general(bv_z,error)
  test_name = "Complex UB to BV, n=40;"
  call z_output_result_upper(test_name,a0_z,a1_z,ubwa,bv_z%ubw,t0,t1,c*tol,error)

  na=1
  lbwa=0; ubwa=0
  ub_z=z_random_ub(na,lbwa,ubwa,error=error)
  a0_z=general(ub_z,error)
  call cpu_time(t0)
  bv_z=bv(ub_z,error)
  call cpu_time(t1)
  a1_z=general(bv_z,error)
  test_name = "Complex UB to BV, n=1;"
  call z_output_result_upper(test_name,a0_z,a1_z,ubwa,bv_z%ubw,t0,t1,c*tol,error)

  na=2
  lbwa=1; ubwa=1
  ub_z=z_random_ub(na,lbwa,ubwa,error=error)
  a0_z=general(ub_z,error)
  call cpu_time(t0)
  bv_z=bv(ub_z,error)
  call cpu_time(t1)
  a1_z=general(bv_z,error)
  test_name = "Complex UB to BV, n=2;"
  call z_output_result_upper(test_name,a0_z,a1_z,ubwa,bv_z%ubw,t0,t1,c*tol,error)

  na=3
  lbwa=1; ubwa=1
  ub_z=z_random_ub(na,lbwa,ubwa,error=error)
  a0_z=general(ub_z,error)
  call cpu_time(t0)
  bv_z=bv(ub_z,error)
  call cpu_time(t1)
  a1_z=general(bv_z,error)
  test_name = "Complex UB to BV, n=3;"
  call z_output_result_upper(test_name,a0_z,a1_z,ubwa,bv_z%ubw,t0,t1,c*tol,error)

  na=4
  lbwa=2; ubwa=2
  ub_z=z_random_ub(na,lbwa,ubwa,error=error)
  a0_z=general(ub_z,error)
  call cpu_time(t0)
  bv_z=bv(ub_z,error)
  call cpu_time(t1)
  a1_z=general(bv_z,error)
  test_name = "Complex UB to BV, n=4;"
  call z_output_result_upper(test_name,a0_z,a1_z,ubwa,bv_z%ubw,t0,t1,c*tol,error)
  print *

  na=40
  lbwa=3; ubwa=5
  bv_z=z_random_bv(na,lbwa,ubwa,error=error)
  a0_z=general(bv_z,error)
  call cpu_time(t0)
  ub_z=ub(bv_z,error)
  call cpu_time(t1)
  a1_z=general(ub_z,error)
  test_name = "Complex BV to UB, n=40;"
  call z_output_result_upper(test_name,a0_z,a1_z,ubwa,ub_z%ubw,t0,t1,c*tol,error)

  na=1
  lbwa=0; ubwa=0
  bv_z=z_random_bv(na,lbwa,ubwa,error=error)
  a0_z=general(bv_z,error)
  call cpu_time(t0)
  ub_z=ub(bv_z,error)
  call cpu_time(t1)
  a1_z=general(ub_z,error)
  test_name = "Complex BV to UB, n=1;"
  call z_output_result_upper(test_name,a0_z,a1_z,ubwa,ub_z%ubw,t0,t1,c*tol,error)

  na=2
  lbwa=1; ubwa=1
  bv_z=z_random_bv(na,lbwa,ubwa,error=error)
  a0_z=general(bv_z,error)
  call cpu_time(t0)
  ub_z=ub(bv_z,error)
  call cpu_time(t1)
  a1_z=general(ub_z,error)
  test_name = "Complex BV to UB, n=2;"
  call z_output_result_upper(test_name,a0_z,a1_z,ubwa,ub_z%ubw,t0,t1,c*tol,error)

  na=3
  lbwa=1; ubwa=1
  bv_z=z_random_bv(na,lbwa,ubwa,error=error)
  a0_z=general(bv_z,error)
  call cpu_time(t0)
  ub_z=ub(bv_z,error)
  call cpu_time(t1)
  a1_z=general(ub_z,error)
  test_name = "Complex BV to UB, n=3;"
  call z_output_result_upper(test_name,a0_z,a1_z,ubwa,ub_z%ubw,t0,t1,c*tol,error)

  na=4
  lbwa=2; ubwa=2
  bv_z=z_random_bv(na,lbwa,ubwa,error=error)
  a0_z=general(bv_z,error)
  call cpu_time(t0)
  ub_z=ub(bv_z,error)
  call cpu_time(t1)
  a1_z=general(ub_z,error)
  test_name = "Complex BV to UB, n=4;"
  call z_output_result_upper(test_name,a0_z,a1_z,ubwa,ub_z%ubw,t0,t1,c*tol,error)

end program test_convert_ub_and_bv
