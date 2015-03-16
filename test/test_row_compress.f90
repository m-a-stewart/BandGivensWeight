program test_row_compress
  use mod_orb
  use mod_test_data
  implicit none

  real(kind=dp) :: t0, t1
  integer(kind=int32) :: na, lbwa, ubwa
  type(error_info) :: error
  real(kind=dp), parameter :: tol=1e-15, c=4.0
  real(kind=dp), dimension(:,:), allocatable :: a_d, a0_d
  complex(kind=dp), dimension(:,:), allocatable :: a_z, a0_z
  !
  type(d_ubt), allocatable :: ubt_d
  type(z_ubt), allocatable :: ubt_z
  type(d_rc), allocatable :: swbv_d
  type(z_rc), allocatable :: swbv_z

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
  ubt_z=z_random_ubt(na,lbwa,ubwa,error=error)
  a_z = general(ubt_z,error)
  a0_z = a_z
  call cpu_time(t0)
  swbv_z=rc(ubt_z,error)
  call cpu_time(t1)  
  a_z = general(swbv_z%bv,error)
  a_z = swbv_z%sw * a_z
  test_name = "Random Complex Row Compression, n=40"
  call z_output_result_lower_upper(test_name,a0_z,a_z,lbwa, &
       swbv_z%bv%lbw, min(ubwa+lbwa,na-1),swbv_z%bv%ubw,t0,t1,c*tol,error)

  na=1; lbwa=0; ubwa=0
  ubt_z=z_random_ubt(na,lbwa,ubwa,error=error)
  a_z = general(ubt_z,error)
  a0_z = a_z
  call cpu_time(t0)
  swbv_z=rc(ubt_z,error)
  call cpu_time(t1)  
  a_z = general(swbv_z%bv,error)
  a_z = swbv_z%sw * a_z
  test_name = "Random Complex Row Compression, n=1"
  call z_output_result_lower_upper(test_name,a0_z,a_z,lbwa, &
       swbv_z%bv%lbw, min(ubwa+lbwa,na-1),swbv_z%bv%ubw,t0,t1,c*tol,error)

  na=2; lbwa=1; ubwa=1
  ubt_z=z_random_ubt(na,lbwa,ubwa,error=error)
  a_z = general(ubt_z,error)
  a0_z = a_z
  call cpu_time(t0)
  swbv_z=rc(ubt_z,error)
  call cpu_time(t1)  
  a_z = general(swbv_z%bv,error)
  a_z = swbv_z%sw * a_z
  test_name = "Random Complex Row Compression, n=2"
  call z_output_result_lower_upper(test_name,a0_z,a_z,lbwa, &
       swbv_z%bv%lbw, min(ubwa+lbwa,na-1),swbv_z%bv%ubw,t0,t1,c*tol,error)

  na=3; lbwa=1; ubwa=1
  ubt_z=z_random_ubt(na,lbwa,ubwa,error=error)
  a_z = general(ubt_z,error)
  a0_z = a_z
  call cpu_time(t0)
  swbv_z=rc(ubt_z,error)
  call cpu_time(t1)  
  a_z = general(swbv_z%bv,error)
  a_z = swbv_z%sw * a_z
  test_name = "Random Complex Row Compression, n=3"
  call z_output_result_lower_upper(test_name,a0_z,a_z,lbwa, &
       swbv_z%bv%lbw, min(ubwa+lbwa,na-1),swbv_z%bv%ubw,t0,t1,c*tol,error)

  na=4; lbwa=2; ubwa=2
  ubt_z=z_random_ubt(na,lbwa,ubwa,error=error)
  a_z = general(ubt_z,error)
  a0_z = a_z
  call cpu_time(t0)
  swbv_z=rc(ubt_z,error)
  call cpu_time(t1)  
  a_z = general(swbv_z%bv,error)
  a_z = swbv_z%sw * a_z
  test_name = "Random Complex Row Compression, n=4"
  call z_output_result_lower_upper(test_name,a0_z,a_z,lbwa, &
       swbv_z%bv%lbw, min(ubwa+lbwa,na-1),swbv_z%bv%ubw,t0,t1,c*tol,error)
  
end program test_row_compress
