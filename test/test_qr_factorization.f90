program test_qr_factorization
  use mod_orth_rank
  use mod_test_data
  implicit none
  real(kind=dp) :: t0, t1
  integer(kind=int32) :: na, lbwa, ubwa
  type(error_info) :: error
  real(kind=dp), parameter :: tol=1e-15, c=2.5
  !

  real(kind=dp), dimension(:,:), allocatable :: a_d, a0_d
  complex(kind=dp), dimension(:,:), allocatable :: a_c, a0_c

  type(d_bv), allocatable :: bv_d
  type(c_bv), allocatable :: bv_c
  type(d_qr), allocatable :: swub_d
  type(c_qr), allocatable :: swub_c

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
  bv_c=c_random_bv(na,lbwa,ubwa,error=error)
  a_c = general(bv_c,error)
  a0_c = a_c
  call cpu_time(t0)
  swub_c=qr(bv_c,error)
  call cpu_time(t1)  
  a_c = general(swub_c%ub,error)
  a_c = swub_c%sw * a_c
  test_name = "Random Complex QR Factorization"
  call c_output_result_lower_upper(test_name,a0_c,a_c,0, &
       swub_c%ub%lbw, min(ubwa+lbwa,na-1),swub_c%ub%ubw,t0,t1,c*tol,error)

  na=1; lbwa=0; ubwa=0
  bv_c=c_random_bv(na,lbwa,ubwa,error=error)
  a_c = general(bv_c,error)
  a0_c = a_c
  call cpu_time(t0)
  swub_c=qr(bv_c,error)
  call cpu_time(t1)  
  a_c = general(swub_c%ub,error)
  a_c = swub_c%sw * a_c
  test_name = "Random Complex QR Factorization, n=1"
  call c_output_result_lower_upper(test_name,a0_c,a_c,0, &
       swub_c%ub%lbw, min(ubwa+lbwa,na-1), swub_c%ub%ubw,t0,t1,c*tol,error)

  na=2; lbwa=1; ubwa=1
  bv_c=c_random_bv(na,lbwa,ubwa,error=error)
  a_c = general(bv_c,error)
  a0_c = a_c
  call cpu_time(t0)
  swub_c=qr(bv_c,error)
  call cpu_time(t1)  
  a_c = general(swub_c%ub,error)
  a_c = swub_c%sw * a_c
  test_name = "Random Complex QR Factorization, n=2"
  call c_output_result_lower_upper(test_name,a0_c,a_c,0, &
       swub_c%ub%lbw, min(ubwa+lbwa,na-1), swub_c%ub%ubw,t0,t1,c*tol,error)

  na=3; lbwa=1; ubwa=1
  bv_c=c_random_bv(na,lbwa,ubwa,error=error)
  a_c = general(bv_c,error)
  a0_c = a_c
  call cpu_time(t0)
  swub_c=qr(bv_c,error)
  call cpu_time(t1)  
  a_c = general(swub_c%ub,error)
  a_c = swub_c%sw * a_c
  test_name = "Random Complex QR Factorization, n=3"
  call c_output_result_lower_upper(test_name,a0_c,a_c,0, &
       swub_c%ub%lbw, min(ubwa+lbwa,na-1), swub_c%ub%ubw,t0,t1,c*tol,error)

  na=4; lbwa=2; ubwa=2
  bv_c=c_random_bv(na,lbwa,ubwa,error=error)
  a_c = general(bv_c,error)
  a0_c = a_c
  call cpu_time(t0)
  swub_c=qr(bv_c,error)
  call cpu_time(t1)  
  a_c = general(swub_c%ub,error)
  a_c = swub_c%sw * a_c
  test_name = "Random Complex QR Factorization, n=4"
  call c_output_result_lower_upper(test_name,a0_c,a_c,0, &
       swub_c%ub%lbw, min(ubwa+lbwa,na-1), swub_c%ub%ubw,t0,t1,c*tol,error)

end program test_qr_factorization
