program test_general_ubt
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
  type(d_ubt), allocatable :: ubt_d
  type(c_ubt), allocatable :: ubt_c

  call initialize_errors

  ! test one
  print *
  print *, "--------------------------------"
  print *
  print *, "Real UBT Decomposition Tests"
  print *

  na=40
  lbwa=5; ubwa=3
  ubt_d=d_random_ubt(na,lbwa,ubwa,error=error)
  a_d=general(ubt_d,error)
  a1_d=a_d
  call cpu_time(t0)
  ubt_d=ubt(a_d,lbwa+1,ubwa+1,tol,error)
  call cpu_time(t1)
  a0_d = general(ubt_d,error)
  test_name="Random Real UBT;"  
  call d_output_result_lower_upper(test_name,a0_d,a1_d,lbwa,ubt_d%lbw,&
       ubwa,ubt_d%ubw,t0,t1,c*tol,error)

  na=50
  lbwa=13; ubwa=3
  ubt_d=d_random_ubt(na,[(lbwa-1, j=1,na-lbwa-1), (lbwa, j=na-lbwa,na) ], &
       [ (ubwa, j=1,na) ], error=error )
  a_d=general(ubt_d,error)
  a1_d=a_d
  ubt_d=ubt(a_d,lbwa+1,ubwa+1,tol,error)
  call cpu_time(t0)
  a0_d = general(ubt_d,error)
  call cpu_time(t1)
  test_name="Random Real Square Termination UBT;"  
  call d_output_result_lower_upper(test_name,a0_d,a1_d,lbwa,ubt_d%lbw,&
       ubwa,ubt_d%ubw,t0,t1,c*tol,error)

  !
  ! Complex UBT test
  !
  print *
  print *, "--------------------------------"
  print *
  print *, "Complex UBT Decomposition Tests"
  print *

  na=40
  lbwa=5; ubwa=3
  ubt_c=c_random_ubt(na,lbwa,ubwa,error=error)
  a_c=general(ubt_c,error)
  a1_c=a_c
  call cpu_time(t0)
  ubt_c=ubt(a_c, lbwa+1, ubwa+1, tol,error)
  call cpu_time(t1)
  a0_c=general(ubt_c,error)
  test_name="Random Complex UBT;"  
  call c_output_result_lower_upper(test_name,a0_c,a1_c,lbwa,ubt_c%lbw, &
       ubwa,ubt_d%ubw,t0,t1,c*tol,error)

  na=50
  lbwa=13; ubwa=3
  ubt_c=c_random_ubt(na,[ (lbwa-1, j=1,na-lbwa-1), (lbwa, j=na-lbwa,na) ], &
       [ (ubwa, j=1,na) ], error=error )
  a_c=general(ubt_c,error)
  a1_c=a_c
  call cpu_time(t0)
  ubt_c=ubt(a_c,lbwa+1,ubwa+1,tol,error)
  call cpu_time(t1)
  a0_c = general(ubt_c,error)
  test_name="Random Complex Square Termination UBT;"  
  call c_output_result_lower_upper(test_name,a0_c,a1_c,lbwa,ubt_c%lbw, &
       ubwa,ubt_d%ubw,t0,t1,c*tol,error)
  print *

end program test_general_ubt
