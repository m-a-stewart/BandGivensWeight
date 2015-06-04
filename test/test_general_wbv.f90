program test_general_wbv
  use mod_orrb
  use mod_test_data
  implicit none
  real(kind=dp) :: t0, t1
  type(error_info) :: error
  integer(kind=int32) :: na, lbwa, ubwa, j
  real(kind=dp), parameter :: tol=1e-14, c=1.5
  !
  real(kind=dp), dimension(:,:), allocatable :: a_d, a0_d, a1_d
  complex(kind=dp), dimension(:,:), allocatable :: a_z, a0_z, a1_z
  type(d_wbv), allocatable :: wbv_d
  type(z_wbv), allocatable :: wbv_z

  call initialize_errors

  ! test one
  print *
  print *, "--------------------------------"
  print *
  print *, "Real WBV Decomposition Tests"
  print *

  na=40
  lbwa=5; ubwa=3
  wbv_d=d_random_wbv(na,lbwa,ubwa,error=error)
  a_d=general(wbv_d,error)
  a1_d=a_d
  call cpu_time(t0)
  wbv_d=wbv(a_d,lbwa+1,ubwa+1,tol,error)
  call cpu_time(t1)
  a0_d = general(wbv_d,error)
  test_name="Random Real WBV;"  
  call d_output_result_lower_upper(test_name,a0_d,a1_d,lbwa,wbv_d%lbw,&
       ubwa,wbv_d%ubw,t0,t1,c*tol,error)

  na=50
  lbwa=13; ubwa=3
  wbv_d=d_random_wbv(na,[ (lbwa-1, j=1,na-lbwa-1), (lbwa, j=na-lbwa,na) ], &
       [ (ubwa, j=1,na) ], error=error )
  a_d=general(wbv_d,error)
  a1_d=a_d
  wbv_d=wbv(a_d,lbwa+1,ubwa+1,tol,error)
  call cpu_time(t0)
  a0_d = general(wbv_d,error)
  call cpu_time(t1)
  test_name="Random Real Square Termination WBV;"  
  call d_output_result_lower_upper(test_name,a0_d,a1_d,lbwa,wbv_d%lbw,&
       ubwa,wbv_d%ubw,t0,t1,c*tol,error)

  !
  ! Complex WBV test
  !
  print *
  print *, "--------------------------------"
  print *
  print *, "Complex WBV Decomposition Tests"
  print *

  na=40
  lbwa=5; ubwa=3
  wbv_z=z_random_wbv(na,lbwa,ubwa,error=error)
  a_z=general(wbv_z,error)
  a1_z=a_z
  call cpu_time(t0)
  wbv_z=wbv(a_z, lbwa+1, ubwa+1, tol,error)
  call cpu_time(t1)
  a0_z=general(wbv_z,error)
  test_name="Random Complex WBV;"  
  call z_output_result_lower_upper(test_name,a0_z,a1_z,lbwa,wbv_z%lbw, &
       ubwa,wbv_z%ubw,t0,t1,c*tol,error)

  na=50
  lbwa=13; ubwa=3
  wbv_z=z_random_wbv(na,[ (lbwa-1, j=1,na-lbwa-1), (lbwa, j=na-lbwa,na) ], &
       [ (ubwa, j=1,na) ], error=error )
  a_z=general(wbv_z,error)
  a1_z=a_z
  call cpu_time(t0)
  wbv_z=wbv(a_z,lbwa+1,ubwa+1,tol,error)
  call cpu_time(t1)
  a0_z = general(wbv_z,error)
  test_name="Random Complex Square Termination WBV;"  
  call z_output_result_lower_upper(test_name,a0_z,a1_z,lbwa,wbv_z%lbw, &
       ubwa,wbv_z%ubw,t0,t1,c*tol,error)
  print *

end program test_general_wbv
