program test_convert_ubt_and_wbv
  use mod_orth_rank
  use mod_test_data
  implicit none
  !
  ! Note: This does not cover large rank cases (e.g. ubw=n-1).
  !

  real(kind=dp) :: t0, t1
  integer(kind=int32) :: na, ubwa, lbwa
  type(error_info) :: error
  real(kind=dp), parameter :: tol=1e-15, c=2.0
  !
  real(kind=dp), dimension(:,:), allocatable :: a0_d, a1_d
  complex(kind=dp), dimension(:,:), allocatable :: a0_c, a1_c

  type(d_wbv), allocatable :: wbv_d
  type(c_wbv), allocatable :: wbv_c
  type(d_ubt), allocatable :: ubt_d
  type(c_ubt), allocatable :: ubt_c

  call initialize_errors

  print *
  print *, "--------------------------------"
  print *
  print *, "Real WBV and UBT Conversion Tests:"
  print *
  na=40
  lbwa=5; ubwa=3
  wbv_d=d_random_wbv(na,lbwa,ubwa,error=error)
  a0_d=general(wbv_d,error)
  call cpu_time(t0)
  ubt_d=ubt(wbv_d,error)
  call cpu_time(t1)
  a1_d=general(ubt_d,error)
  test_name = "Real WBV to UBT, n=40;"
  call d_output_result_lower_upper(test_name,a0_d,a1_d,lbwa,ubt_d%lbw, &
       ubwa, ubt_d%ubw, t0,t1,c*tol,error)

  na=1
  lbwa=0; ubwa=0
  wbv_d=d_random_wbv(na,lbwa,ubwa,error=error)
  a0_d=general(wbv_d,error)
  call cpu_time(t0)
  ubt_d=ubt(wbv_d,error)
  call cpu_time(t1)
  a1_d=general(ubt_d,error)
  test_name = "Real WBV to UBT, n=1;"
  call d_output_result_lower_upper(test_name,a0_d,a1_d,lbwa,ubt_d%lbw, &
       ubwa, ubt_d%ubw, t0,t1,c*tol,error)


  na=2
  lbwa=1; ubwa=1
  wbv_d=d_random_wbv(na,lbwa,ubwa,error=error)
  a0_d=general(wbv_d,error)
  call cpu_time(t0)
  ubt_d=ubt(wbv_d,error)
  call cpu_time(t1)
  a1_d=general(ubt_d,error)
  test_name = "Real WBV to UBT, n=2;"
  call d_output_result_lower_upper(test_name,a0_d,a1_d,lbwa,ubt_d%lbw, &
       ubwa, ubt_d%ubw, t0,t1,c*tol,error)

  
  na=3
  lbwa=1; ubwa=1
  wbv_d=d_random_wbv(na,lbwa,ubwa,error=error)
  a0_d=general(wbv_d,error)
  call cpu_time(t0)
  ubt_d=ubt(wbv_d,error)
  call cpu_time(t1)
  a1_d=general(ubt_d,error)
  test_name = "Real WBV to UBT, n=3;"
  call d_output_result_lower_upper(test_name,a0_d,a1_d,lbwa,ubt_d%lbw, &
       ubwa, ubt_d%ubw, t0,t1,c*tol,error)


  na=4
  lbwa=2; ubwa=2
  wbv_d=d_random_wbv(na,lbwa,ubwa,error=error)
  a0_d=general(wbv_d,error)
  call cpu_time(t0)
  ubt_d=ubt(wbv_d,error)
  call cpu_time(t1)
  a1_d=general(ubt_d,error)
  test_name = "Real WBV to UBT, n=4;"
  call d_output_result_lower_upper(test_name,a0_d,a1_d,lbwa,ubt_d%lbw, &
       ubwa, ubt_d%ubw, t0,t1,c*tol,error)
  print *
  
  na=40
  lbwa=5; ubwa=3
  ubt_d=d_random_ubt(na,lbwa,ubwa,error=error)
  a0_d=general(ubt_d,error)
  call cpu_time(t0)
  wbv_d=wbv(ubt_d,error)
  call cpu_time(t1)
  a1_d=general(wbv_d,error)
  test_name = "Real UBT to WBV, n=40;"
  call d_output_result_lower_upper(test_name,a0_d,a1_d,lbwa,wbv_d%lbw, &
       ubwa, wbv_d%ubw, t0,t1,c*tol,error)

  na=1
  lbwa=0; ubwa=0
  ubt_d=d_random_ubt(na,lbwa,ubwa,error=error)
  a0_d=general(ubt_d,error)
  call cpu_time(t0)
  wbv_d=wbv(ubt_d,error)
  call cpu_time(t1)
  a1_d=general(wbv_d,error)
  test_name = "Real UBT to WBV, n=1;"
  call d_output_result_lower_upper(test_name,a0_d,a1_d,lbwa,wbv_d%lbw, &
       ubwa, wbv_d%ubw, t0,t1,c*tol,error)


  na=2
  lbwa=1; ubwa=1
  ubt_d=d_random_ubt(na,lbwa,ubwa,error=error)
  a0_d=general(ubt_d,error)
  call cpu_time(t0)
  wbv_d=wbv(ubt_d,error)
  call cpu_time(t1)
  a1_d=general(wbv_d,error)
  test_name = "Real UBT to WBV, n=2;"
  call d_output_result_lower_upper(test_name,a0_d,a1_d,lbwa,wbv_d%lbw, &
       ubwa, wbv_d%ubw, t0,t1,c*tol,error)


  na=3
  lbwa=1; ubwa=1
  ubt_d=d_random_ubt(na,lbwa,ubwa,error=error)
  a0_d=general(ubt_d,error)
  call cpu_time(t0)
  wbv_d=wbv(ubt_d,error)
  call cpu_time(t1)
  a1_d=general(wbv_d,error)
  test_name = "Real UBT to WBV, n=3;"
  call d_output_result_lower_upper(test_name,a0_d,a1_d,lbwa,wbv_d%lbw, &
       ubwa, wbv_d%ubw, t0,t1,c*tol,error)

  na=4
  lbwa=2; ubwa=2
  ubt_d=d_random_ubt(na,lbwa,ubwa,error=error)
  a0_d=general(ubt_d,error)
  call cpu_time(t0)
  wbv_d=wbv(ubt_d,error)
  call cpu_time(t1)
  a1_d=general(wbv_d,error)
  test_name = "Real UBT to WBV, n=4;"
  call d_output_result_lower_upper(test_name,a0_d,a1_d,lbwa,wbv_d%lbw, &
       ubwa, wbv_d%ubw, t0,t1,c*tol,error)

  print *
  print *, "--------------------------------"
  print *
  print *, "Complex WBV and UBT Conversion Tests:"
  print *

  na=40
  lbwa=5; ubwa=3
  wbv_c=c_random_wbv(na,lbwa,ubwa,error=error)
  a0_c=general(wbv_c,error)
  call cpu_time(t0)
  ubt_c=ubt(wbv_c,error)
  call cpu_time(t1)
  a1_c=general(ubt_c,error)
  test_name = "Complex WBV to UBT, n=40;"
  call c_output_result_lower_upper(test_name,a0_c,a1_c,lbwa,ubt_c%lbw, &
       ubwa, ubt_c%ubw, t0,t1,c*tol,error)

  na=1
  lbwa=0; ubwa=0
  wbv_c=c_random_wbv(na,lbwa,ubwa,error=error)
  a0_c=general(wbv_c,error)
  call cpu_time(t0)
  ubt_c=ubt(wbv_c,error)
  call cpu_time(t1)
  a1_c=general(ubt_c,error)
  test_name = "Complex WBV to UBT, n=1;"
  call c_output_result_lower_upper(test_name,a0_c,a1_c,lbwa,ubt_c%lbw, &
       ubwa, ubt_c%ubw, t0,t1,c*tol,error)

  na=2
  lbwa=1; ubwa=1
  wbv_c=c_random_wbv(na,lbwa,ubwa,error=error)
  a0_c=general(wbv_c,error)
  call cpu_time(t0)
  ubt_c=ubt(wbv_c,error)
  call cpu_time(t1)
  a1_c=general(ubt_c,error)
  test_name = "Complex WBV to UBT, n=2;"
  call c_output_result_lower_upper(test_name,a0_c,a1_c,lbwa,ubt_c%lbw, &
       ubwa, ubt_c%ubw, t0,t1,c*tol,error)
  
  na=3
  lbwa=1; ubwa=1
  wbv_c=c_random_wbv(na,lbwa,ubwa,error=error)
  a0_c=general(wbv_c,error)
  call cpu_time(t0)
  ubt_c=ubt(wbv_c,error)
  call cpu_time(t1)
  a1_c=general(ubt_c,error)
  test_name = "Complex WBV to UBT, n=3;"
  call c_output_result_lower_upper(test_name,a0_c,a1_c,lbwa,ubt_c%lbw, &
       ubwa, ubt_c%ubw, t0,t1,c*tol,error)

  na=4
  lbwa=2; ubwa=2
  wbv_c=c_random_wbv(na,lbwa,ubwa,error=error)
  a0_c=general(wbv_c,error)
  call cpu_time(t0)
  ubt_c=ubt(wbv_c,error)
  call cpu_time(t1)
  a1_c=general(ubt_c,error)
  test_name = "Complex WBV to UBT, n=4;"
  call c_output_result_lower_upper(test_name,a0_c,a1_c,lbwa,ubt_c%lbw, &
       ubwa, ubt_c%ubw, t0,t1,c*tol,error)
  print *
  
  na=40
  lbwa=5; ubwa=3
  ubt_c=c_random_ubt(na,lbwa,ubwa,error=error)
  a0_c=general(ubt_c,error)
  call cpu_time(t0)
  wbv_c=wbv(ubt_c,error)
  call cpu_time(t1)
  a1_c=general(wbv_c,error)
  test_name = "Complex UBT to WBV, n=40;"
  call c_output_result_lower_upper(test_name,a0_c,a1_c,lbwa,wbv_c%lbw, &
       ubwa,wbv_c%ubw,t0,t1,c*tol,error)

  na=1
  lbwa=0; ubwa=0
  ubt_c=c_random_ubt(na,lbwa,ubwa,error=error)
  a0_c=general(ubt_c,error)
  call cpu_time(t0)
  wbv_c=wbv(ubt_c,error)
  call cpu_time(t1)
  a1_c=general(wbv_c,error)
  test_name = "Complex UBT to WBV, n=1;"
  call c_output_result_lower_upper(test_name,a0_c,a1_c,lbwa,wbv_c%lbw, &
       ubwa,wbv_c%ubw,t0,t1,c*tol,error)

  na=2
  lbwa=1; ubwa=1
  ubt_c=c_random_ubt(na,lbwa,ubwa,error=error)
  a0_c=general(ubt_c,error)
  call cpu_time(t0)
  wbv_c=wbv(ubt_c,error)
  call cpu_time(t1)
  a1_c=general(wbv_c,error)
  test_name = "Complex UBT to WBV, n=2;"
  call c_output_result_lower_upper(test_name,a0_c,a1_c,lbwa,wbv_c%lbw, &
       ubwa,wbv_c%ubw,t0,t1,c*tol,error)

  na=3
  lbwa=1; ubwa=1
  ubt_c=c_random_ubt(na,lbwa,ubwa,error=error)
  a0_c=general(ubt_c,error)
  call cpu_time(t0)
  wbv_c=wbv(ubt_c,error)
  call cpu_time(t1)
  a1_c=general(wbv_c,error)
  test_name = "Complex UBT to WBV, n=3;"
  call c_output_result_lower_upper(test_name,a0_c,a1_c,lbwa,wbv_c%lbw, &
       ubwa,wbv_c%ubw,t0,t1,c*tol,error)

  na=4
  lbwa=2; ubwa=2
  ubt_c=c_random_ubt(na,lbwa,ubwa,error=error)
  a0_c=general(ubt_c,error)
  call cpu_time(t0)
  wbv_c=wbv(ubt_c,error)
  call cpu_time(t1)
  a1_c=general(wbv_c,error)
  test_name = "Complex UBT to WBV, n=4;"
  call c_output_result_lower_upper(test_name,a0_c,a1_c,lbwa,wbv_c%lbw, &
       ubwa,wbv_c%ubw,t0,t1,c*tol,error)


end program test_convert_ubt_and_wbv
