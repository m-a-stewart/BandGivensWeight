program test_general_bt

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
  type(d_bt), allocatable :: bt_d
  type(c_bt), allocatable :: bt_c

  call initialize_errors

  ! test one
  print *
  print *, "--------------------------------"
  print *
  print *, "Real BT Decomposition Tests"
  print *

  na=40
  lbwa=5; ubwa=3
  bt_d=d_random_bt(na,lbwa,ubwa,error=error)
  a_d=general_of(bt_d,error)
  a1_d=a_d
  call cpu_time(t0)
  bt_d=bt_of_general(a_d,ubwa,lbwa+1,ubwa,tol,error)
  call cpu_time(t1)
  a0_d = general_of(bt_d,error)
  test_name="Random Real BT;"  
  call d_output_result_lower(test_name,a0_d,a1_d,lbwa,bt_d%lbw,t0,t1,c*tol,error)

  na=50
  lbwa=13; ubwa=3
  bt_d=d_random_bt(na,(/ (lbwa-1, j=1,na-lbwa), (lbwa, j=na-lbwa+1,na) /), &
       (/ (ubwa, j=1,na) /), error=error )
  a_d=general_of(bt_d,error)
  a1_d=a_d
  bt_d=bt_of_general(a_d,ubwa,lbwa+1,ubwa,tol,error)
  call cpu_time(t0)
  a0_d = general_of(bt_d,error)
  call cpu_time(t1)
  test_name="Random Real Square Termination BT;"  
  call d_output_result_lower(test_name,a0_d,a1_d,lbwa,bt_d%lbw,t0,t1,c*tol,error)
  deallocate(bt_d)

  !
  ! Complex BT test
  !
  print *
  print *, "--------------------------------"
  print *
  print *, "Complex BT Decomposition Tests"
  print *

  na=40
  lbwa=5; ubwa=3
  bt_c=c_random_bt(na,lbwa,ubwa,error=error)
  a_c=general_of(bt_c,error)
  a1_c=a_c
  call cpu_time(t0)
  bt_c=bt_of_general(a_c, ubwa, lbwa+1, ubwa, tol,error)
  call cpu_time(t1)
  a0_c=general_of(bt_c,error)
  test_name="Random Complex BT;"  
  call c_output_result_lower(test_name,a0_c,a1_c,lbwa,bt_c%lbw,t0,t1,c*tol,error)

  na=50
  lbwa=13; ubwa=3
  bt_c=c_random_bt(na,(/ (lbwa-1, j=1,na-lbwa), (lbwa, j=na-lbwa+1,na) /), &
       (/ (ubwa, j=1,na) /), error=error )
  a_c=general_of(bt_c,error)
  a1_c=a_c
  call cpu_time(t0)
  bt_c=bt_of_general(a_c,ubwa,lbwa+1,ubwa,tol,error)
  call cpu_time(t1)
  a0_c = general_of(bt_c,error)
  test_name="Random Complex Square Termination BT;"  
  call c_output_result_lower(test_name,a0_c,a1_c,lbwa,bt_c%lbw,t0,t1,c*tol,error)
  print *

  
end program test_general_bt
