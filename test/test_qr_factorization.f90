program test_qr_factorization
  use mod_orb
  use mod_test_data
  implicit none
  real(kind=dp) :: t0, t1
  integer(kind=int32) :: na, lbwa, ubwa
  type(error_info) :: error
  real(kind=dp), parameter :: tol=1e-15, c=2.5
  !

  real(kind=dp), dimension(:,:), allocatable :: a_d, a0_d
  complex(kind=dp), dimension(:,:), allocatable :: a_z, a0_z

  type(d_bv), allocatable :: bv_d
  type(z_bv), allocatable :: bv_z
  type(d_qr), allocatable :: swub_d
  type(z_qr), allocatable :: swub_z

  call initialize_errors
  print *
  print *, "--------------------------------"
  print *
  print *, "Real QR Factorization Tests"
  print *
  !
  ! full qr factorization
  !
  na=40; lbwa=3; ubwa=5
  bv_d=d_random_bv(na,lbwa,ubwa,error=error)
  a_d = general(bv_d,error)
  a0_d = a_d
  call cpu_time(t0)
  swub_d=qr(bv_d,error)
  call cpu_time(t1)  
  a_d = general(swub_d%ub,error)
  a_d = swub_d%sw * a_d
  test_name = "Random Real QR Factorization"
  call d_output_result_lower_upper(test_name,a0_d,a_d,0, &
       swub_d%ub%lbw, min(ubwa+lbwa,na-1),swub_d%ub%ubw,t0,t1,c*tol,error)

  na=1; lbwa=0; ubwa=0
  bv_d=d_random_bv(na,lbwa,ubwa,error=error)
  a_d = general(bv_d,error)
  a0_d = a_d
  call cpu_time(t0)
  swub_d=qr(bv_d,error)
  call cpu_time(t1)  
  a_d = general(swub_d%ub,error)
  a_d = swub_d%sw * a_d
  test_name = "Random Real QR Factorization, n=1"
  call d_output_result_lower_upper(test_name,a0_d,a_d,0, &
       swub_d%ub%lbw, min(ubwa+lbwa,na-1), swub_d%ub%ubw,t0,t1,c*tol,error)

  na=2; lbwa=1; ubwa=1
  bv_d=d_random_bv(na,lbwa,ubwa,error=error)
  a_d = general(bv_d,error)
  a0_d = a_d
  call cpu_time(t0)
  swub_d=qr(bv_d,error)
  call cpu_time(t1)  
  a_d = general(swub_d%ub,error)
  a_d = swub_d%sw * a_d
  test_name = "Random Real QR Factorization, n=2"
  call d_output_result_lower_upper(test_name,a0_d,a_d,0, &
       swub_d%ub%lbw, min(ubwa+lbwa,na-1), swub_d%ub%ubw,t0,t1,c*tol,error)

  na=3; lbwa=1; ubwa=1
  bv_d=d_random_bv(na,lbwa,ubwa,error=error)
  a_d = general(bv_d,error)
  a0_d = a_d
  call cpu_time(t0)
  swub_d=qr(bv_d,error)
  call cpu_time(t1)  
  a_d = general(swub_d%ub,error)
  a_d = swub_d%sw * a_d
  test_name = "Random Real QR Factorization, n=3"
  call d_output_result_lower_upper(test_name,a0_d,a_d,0, &
       swub_d%ub%lbw, min(ubwa+lbwa,na-1), swub_d%ub%ubw,t0,t1,c*tol,error)

  na=4; lbwa=2; ubwa=2
  bv_d=d_random_bv(na,lbwa,ubwa,error=error)
  a_d = general(bv_d,error)
  a0_d = a_d
  call cpu_time(t0)
  swub_d=qr(bv_d,error)
  call cpu_time(t1)  
  a_d = general(swub_d%ub,error)
  a_d = swub_d%sw * a_d
  test_name = "Random Real QR Factorization, n=4"
  call d_output_result_lower_upper(test_name,a0_d,a_d,0, &
       swub_d%ub%lbw, min(ubwa+lbwa,na-1), swub_d%ub%ubw,t0,t1,c*tol,error)
  
  !
  print *
  print *, "--------------------------------"
  print *
  print *, "Complex QR Factorization Tests"
  print *

  na=40; lbwa=3; ubwa=5
  bv_z=z_random_bv(na,lbwa,ubwa,error=error)
  a_z = general(bv_z,error)
  a0_z = a_z
  call cpu_time(t0)
  swub_z=qr(bv_z,error)
  call cpu_time(t1)  
  a_z = general(swub_z%ub,error)
  a_z = swub_z%sw * a_z
  test_name = "Random Complex QR Factorization"
  call z_output_result_lower_upper(test_name,a0_z,a_z,0, &
       swub_z%ub%lbw, min(ubwa+lbwa,na-1),swub_z%ub%ubw,t0,t1,c*tol,error)

  na=1; lbwa=0; ubwa=0
  bv_z=z_random_bv(na,lbwa,ubwa,error=error)
  a_z = general(bv_z,error)
  a0_z = a_z
  call cpu_time(t0)
  swub_z=qr(bv_z,error)
  call cpu_time(t1)  
  a_z = general(swub_z%ub,error)
  a_z = swub_z%sw * a_z
  test_name = "Random Complex QR Factorization, n=1"
  call z_output_result_lower_upper(test_name,a0_z,a_z,0, &
       swub_z%ub%lbw, min(ubwa+lbwa,na-1), swub_z%ub%ubw,t0,t1,c*tol,error)

  na=2; lbwa=1; ubwa=1
  bv_z=z_random_bv(na,lbwa,ubwa,error=error)
  a_z = general(bv_z,error)
  a0_z = a_z
  call cpu_time(t0)
  swub_z=qr(bv_z,error)
  call cpu_time(t1)  
  a_z = general(swub_z%ub,error)
  a_z = swub_z%sw * a_z
  test_name = "Random Complex QR Factorization, n=2"
  call z_output_result_lower_upper(test_name,a0_z,a_z,0, &
       swub_z%ub%lbw, min(ubwa+lbwa,na-1), swub_z%ub%ubw,t0,t1,c*tol,error)

  na=3; lbwa=1; ubwa=1
  bv_z=z_random_bv(na,lbwa,ubwa,error=error)
  a_z = general(bv_z,error)
  a0_z = a_z
  call cpu_time(t0)
  swub_z=qr(bv_z,error)
  call cpu_time(t1)  
  a_z = general(swub_z%ub,error)
  a_z = swub_z%sw * a_z
  test_name = "Random Complex QR Factorization, n=3"
  call z_output_result_lower_upper(test_name,a0_z,a_z,0, &
       swub_z%ub%lbw, min(ubwa+lbwa,na-1), swub_z%ub%ubw,t0,t1,c*tol,error)

  na=4; lbwa=2; ubwa=2
  bv_z=z_random_bv(na,lbwa,ubwa,error=error)
  a_z = general(bv_z,error)
  a0_z = a_z
  call cpu_time(t0)
  swub_z=qr(bv_z,error)
  call cpu_time(t1)  
  a_z = general(swub_z%ub,error)
  a_z = swub_z%sw * a_z
  test_name = "Random Complex QR Factorization, n=4"
  call z_output_result_lower_upper(test_name,a0_z,a_z,0, &
       swub_z%ub%lbw, min(ubwa+lbwa,na-1), swub_z%ub%ubw,t0,t1,c*tol,error)

end program test_qr_factorization
