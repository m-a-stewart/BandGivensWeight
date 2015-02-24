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
  allocate(a_d(na,na), a0_d(na,na))
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
       swub_d%ub%lbw, ubwa+lbwa,swub_d%ub%ubw,t0,t1,c*tol,error)

  !
  print *
  print *, "--------------------------------"
  print *
  print *, "Complex QR Factorization Tests"
  print *

  na=40; lbwa=3; ubwa=5
  allocate(a_c(na,na), a0_c(na,na))
  bv_c=c_random_bv(na,lbwa,ubwa)
  a_c = general(bv_c)
  a0_c = a_c
  call cpu_time(t0)
  swub_c=qr(bv_c)
  call cpu_time(t1)
  a_c(1:na,1:na) = general(swub_c%ub)
  a_c=swub_c%sw * a_c
  test_name = "Random Complex QR Factorization"
  call c_output_result_lower_upper(test_name,a0_c,a_c,0, &
       swub_c%ub%lbw, ubwa+lbwa,swub_c%ub%ubw,t0,t1,c*tol,error)

end program test_qr_factorization
