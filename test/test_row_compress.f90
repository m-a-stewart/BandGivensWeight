program test_row_compress
  use mod_orth_rank
  use mod_test_data
  implicit none

  real(kind=dp) :: t0, t1
  integer(kind=int32) :: na, lbwa, ubwa
  type(error_info) :: error
  real(kind=dp), parameter :: tol=1e-15, c=4.0
  real(kind=dp), dimension(:,:), allocatable :: a_d, a0_d
  complex(kind=dp), dimension(:,:), allocatable :: a_c, a0_c
  !
  type(d_ubt), allocatable :: ubt_d
  type(c_ubt), allocatable :: ubt_c
  type(d_rc), allocatable :: swbv_d
  type(c_rc), allocatable :: swbv_c
  type(d_sweeps), allocatable :: sw_d
  type(c_sweeps), allocatable :: sw_c

  call initialize_errors

  print *
  print *, "--------------------------------"
  print *
  print *, "Real Row Compression Tests:"
  print *

  na=40; lbwa=7; ubwa=5
  ubt_d=d_random_ubt(na,lbwa,ubwa,error=error)
  a_d = general(ubt_d,error)
  a0_d = a_d
  call cpu_time(t0)
  swbv_d=rc(ubt_d,error)
  call cpu_time(t1)  
  a_d = general(swbv_d%bv,error)
  a_d = swbv_d%sw * a_d
  test_name = "Random Real Row Compression, n=40"
  call d_output_result_lower_upper(test_name,a0_d,a_d,lbwa, &
       swbv_d%bv%lbw, min(ubwa+lbwa,na-1),swbv_d%bv%ubw,t0,t1,c*tol,error)

  na=1; lbwa=0; ubwa=0
  ubt_d=d_random_ubt(na,lbwa,ubwa,error=error)
  a_d = general(ubt_d,error)
  a0_d = a_d
  call cpu_time(t0)
  swbv_d=rc(ubt_d,error)
  call cpu_time(t1)  
  a_d = general(swbv_d%bv,error)
  a_d = swbv_d%sw * a_d
  test_name = "Random Real Row Compression, n=1"
  call d_output_result_lower_upper(test_name,a0_d,a_d,lbwa, &
       swbv_d%bv%lbw, min(ubwa+lbwa,na-1),swbv_d%bv%ubw,t0,t1,c*tol,error)

  na=2; lbwa=1; ubwa=1
  ubt_d=d_random_ubt(na,lbwa,ubwa,error=error)
  a_d = general(ubt_d,error)
  a0_d = a_d
  call cpu_time(t0)
  swbv_d=rc(ubt_d,error)
  call cpu_time(t1)  
  a_d = general(swbv_d%bv,error)
  a_d = swbv_d%sw * a_d
  test_name = "Random Real Row Compression, n=2"
  call d_output_result_lower_upper(test_name,a0_d,a_d,lbwa, &
       swbv_d%bv%lbw, min(ubwa+lbwa,na-1),swbv_d%bv%ubw,t0,t1,c*tol,error)

  na=3; lbwa=1; ubwa=1
  ubt_d=d_random_ubt(na,lbwa,ubwa,error=error)
  a_d = general(ubt_d,error)
  a0_d = a_d
  call cpu_time(t0)
  swbv_d=rc(ubt_d,error)
  call cpu_time(t1)  
  a_d = general(swbv_d%bv,error)
  a_d = swbv_d%sw * a_d
  test_name = "Random Real Row Compression, n=3"
  call d_output_result_lower_upper(test_name,a0_d,a_d,lbwa, &
       swbv_d%bv%lbw, min(ubwa+lbwa,na-1),swbv_d%bv%ubw,t0,t1,c*tol,error)

  na=4; lbwa=2; ubwa=2
  ubt_d=d_random_ubt(na,lbwa,ubwa,error=error)
  a_d = general(ubt_d,error)
  a0_d = a_d
  call cpu_time(t0)
  swbv_d=rc(ubt_d,error)
  call cpu_time(t1)  
  a_d = general(swbv_d%bv,error)
  a_d = swbv_d%sw * a_d
  test_name = "Random Real Row Compression, n=4"
  call d_output_result_lower_upper(test_name,a0_d,a_d,lbwa, &
       swbv_d%bv%lbw, min(ubwa+lbwa,na-1),swbv_d%bv%ubw,t0,t1,c*tol,error)
  
  print *
  print *, "--------------------------------"
  print *
  print *, "Complex Row Compression Tests:"
  print *

  na=40; lbwa=7; ubwa=5
  ubt_c=c_random_ubt(na,lbwa,ubwa,error=error)
  a_c = general(ubt_c,error)
  a0_c = a_c
  call cpu_time(t0)
  swbv_c=rc(ubt_c,error)
  call cpu_time(t1)  
  a_c = general(swbv_c%bv,error)
  a_c = swbv_c%sw * a_c
  test_name = "Random Complex Row Compression, n=40"
  call c_output_result_lower_upper(test_name,a0_c,a_c,lbwa, &
       swbv_c%bv%lbw, min(ubwa+lbwa,na-1),swbv_c%bv%ubw,t0,t1,c*tol,error)

  na=1; lbwa=0; ubwa=0
  ubt_c=c_random_ubt(na,lbwa,ubwa,error=error)
  a_c = general(ubt_c,error)
  a0_c = a_c
  call cpu_time(t0)
  swbv_c=rc(ubt_c,error)
  call cpu_time(t1)  
  a_c = general(swbv_c%bv,error)
  a_c = swbv_c%sw * a_c
  test_name = "Random Complex Row Compression, n=1"
  call c_output_result_lower_upper(test_name,a0_c,a_c,lbwa, &
       swbv_c%bv%lbw, min(ubwa+lbwa,na-1),swbv_c%bv%ubw,t0,t1,c*tol,error)

  na=2; lbwa=1; ubwa=1
  ubt_c=c_random_ubt(na,lbwa,ubwa,error=error)
  a_c = general(ubt_c,error)
  a0_c = a_c
  call cpu_time(t0)
  swbv_c=rc(ubt_c,error)
  call cpu_time(t1)  
  a_c = general(swbv_c%bv,error)
  a_c = swbv_c%sw * a_c
  test_name = "Random Complex Row Compression, n=2"
  call c_output_result_lower_upper(test_name,a0_c,a_c,lbwa, &
       swbv_c%bv%lbw, min(ubwa+lbwa,na-1),swbv_c%bv%ubw,t0,t1,c*tol,error)

  na=3; lbwa=1; ubwa=1
  ubt_c=c_random_ubt(na,lbwa,ubwa,error=error)
  a_c = general(ubt_c,error)
  a0_c = a_c
  call cpu_time(t0)
  swbv_c=rc(ubt_c,error)
  call cpu_time(t1)  
  a_c = general(swbv_c%bv,error)
  a_c = swbv_c%sw * a_c
  test_name = "Random Complex Row Compression, n=3"
  call c_output_result_lower_upper(test_name,a0_c,a_c,lbwa, &
       swbv_c%bv%lbw, min(ubwa+lbwa,na-1),swbv_c%bv%ubw,t0,t1,c*tol,error)

  na=4; lbwa=2; ubwa=2
  ubt_c=c_random_ubt(na,lbwa,ubwa,error=error)
  a_c = general(ubt_c,error)
  a0_c = a_c
  call cpu_time(t0)
  swbv_c=rc(ubt_c,error)
  call cpu_time(t1)  
  a_c = general(swbv_c%bv,error)
  a_c = swbv_c%sw * a_c
  test_name = "Random Complex Row Compression, n=4"
  call c_output_result_lower_upper(test_name,a0_c,a_c,lbwa, &
       swbv_c%bv%lbw, min(ubwa+lbwa,na-1),swbv_c%bv%ubw,t0,t1,c*tol,error)
  
end program test_row_compress
