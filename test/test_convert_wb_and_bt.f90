program test_convert_wb_and_bt
  use mod_orrb
  implicit none

  !
  ! Note: This does not cover large rank cases (e.g. ubw=n-1).
  !

  character(len=40) :: test_name
  character(len=*), parameter :: fmt="(A40, 'Time: ',ES8.2,', lbw: ',I3,', error: ',ES8.2, ', ', A10)"

  real(kind=dp) :: t0, t1
  integer(kind=int32) :: na, ubwa, lbwa
  type(error_info) :: error
  real(kind=dp), parameter :: tol=1e-15, c=2.0
  !
  real(kind=dp), dimension(:,:), allocatable :: a0_d, a1_d
  complex(kind=dp), dimension(:,:), allocatable :: a0_z, a1_z

  type(d_wb), allocatable :: wb_d
  type(z_wb), allocatable :: wb_z
  type(d_bt), allocatable :: bt_d
  type(z_bt), allocatable :: bt_z

  call initialize_errors

  print *
  print *, "--------------------------------"
  print *
  print *, "Real WB and BT Conversion Tests:"
  print *
  na=40
  lbwa=5; ubwa=3
  wb_d=d_random_wb(na,lbwa,ubwa,error=error)
  a0_d=general(wb_d,error)
  call cpu_time(t0)
  bt_d=bt(wb_d,error)
  call cpu_time(t1)
  a1_d=general(bt_d,error)
  test_name = "Real WB to BT, n=40;"
  call d_output_result(test_name,a0_d,a1_d,lbwa,bt_d%lbw,t0,t1,c*tol,error)

  na=1
  lbwa=0; ubwa=0
  wb_d=d_random_wb(na,lbwa,ubwa,error=error)
  a0_d=general(wb_d,error)
  call cpu_time(t0)
  bt_d=bt(wb_d,error)
  call cpu_time(t1)
  a1_d=general(bt_d,error)
  test_name = "Real WB to BT, n=1;"
  call d_output_result(test_name,a0_d,a1_d,lbwa,bt_d%lbw,t0,t1,c*tol,error)

  na=2
  lbwa=1; ubwa=1
  wb_d=d_random_wb(na,lbwa,ubwa,error=error)
  a0_d=general(wb_d,error)
  call cpu_time(t0)
  bt_d=bt(wb_d,error)
  call cpu_time(t1)
  a1_d=general(bt_d,error)
  test_name = "Real WB to BT, n=2;"
  call d_output_result(test_name,a0_d,a1_d,lbwa,bt_d%lbw,t0,t1,c*tol,error)
  
  na=3
  lbwa=1; ubwa=1
  wb_d=d_random_wb(na,lbwa,ubwa,error=error)
  a0_d=general(wb_d,error)
  call cpu_time(t0)
  bt_d=bt(wb_d,error)
  call cpu_time(t1)
  a1_d=general(bt_d,error)
  test_name = "Real WB to BT, n=3;"
  call d_output_result(test_name,a0_d,a1_d,lbwa,bt_d%lbw,t0,t1,c*tol,error)

  na=4
  lbwa=2; ubwa=2
  wb_d=d_random_wb(na,lbwa,ubwa,error=error)
  a0_d=general(wb_d,error)
  call cpu_time(t0)
  bt_d=bt(wb_d,error)
  call cpu_time(t1)
  a1_d=general(bt_d,error)
  test_name = "Real WB to BT, n=4;"
  call d_output_result(test_name,a0_d,a1_d,lbwa,bt_d%lbw,t0,t1,c*tol,error)
  print *
  
  na=40
  lbwa=5; ubwa=3
  bt_d=d_random_bt(na,lbwa,ubwa,error=error)
  a0_d=general(bt_d,error)
  call cpu_time(t0)
  wb_d=wb(bt_d,error)
  call cpu_time(t1)
  a1_d=general(wb_d,error)
  test_name = "Real BT to WB, n=40;"
  call d_output_result(test_name,a0_d,a1_d,lbwa,wb_d%lbw,t0,t1,c*tol,error)

  na=1
  lbwa=0; ubwa=0
  bt_d=d_random_bt(na,lbwa,ubwa,error=error)
  a0_d=general(bt_d,error)
  call cpu_time(t0)
  wb_d=wb(bt_d,error)
  call cpu_time(t1)
  a1_d=general(wb_d,error)
  test_name = "Real BT to WB, n=1;"
  call d_output_result(test_name,a0_d,a1_d,lbwa,wb_d%lbw,t0,t1,c*tol,error)

  na=2
  lbwa=1; ubwa=1
  bt_d=d_random_bt(na,lbwa,ubwa,error=error)
  a0_d=general(bt_d,error)
  call cpu_time(t0)
  wb_d=wb(bt_d,error)
  call cpu_time(t1)
  a1_d=general(wb_d,error)
  test_name = "Real BT to WB, n=2;"
  call d_output_result(test_name,a0_d,a1_d,lbwa,wb_d%lbw,t0,t1,c*tol,error)

  na=3
  lbwa=1; ubwa=1
  bt_d=d_random_bt(na,lbwa,ubwa,error=error)
  a0_d=general(bt_d,error)
  call cpu_time(t0)
  wb_d=wb(bt_d,error)
  call cpu_time(t1)
  a1_d=general(wb_d,error)
  test_name = "Real BT to WB, n=3;"
  call d_output_result(test_name,a0_d,a1_d,lbwa,wb_d%lbw,t0,t1,c*tol,error)

  na=4
  lbwa=2; ubwa=2
  bt_d=d_random_bt(na,lbwa,ubwa,error=error)
  a0_d=general(bt_d,error)
  call cpu_time(t0)
  wb_d=wb(bt_d,error)
  call cpu_time(t1)
  a1_d=general(wb_d,error)
  test_name = "Real BT to WB, n=4;"
  call d_output_result(test_name,a0_d,a1_d,lbwa,wb_d%lbw,t0,t1,c*tol,error)
  

  print *
  print *, "--------------------------------"
  print *
  print *, "Complex WB and BT Conversion Tests:"
  print *

  na=40
  lbwa=5; ubwa=3
  wb_z=z_random_wb(na,lbwa,ubwa,error=error)
  a0_z=general(wb_z,error)
  call cpu_time(t0)
  bt_z=bt(wb_z,error)
  call cpu_time(t1)
  a1_z=general(bt_z,error)
  test_name = "Complex WB to BT, n=40;"
  call z_output_result(test_name,a0_z,a1_z,lbwa,bt_z%lbw,t0,t1,c*tol,error)

  na=1
  lbwa=0; ubwa=0
  wb_z=z_random_wb(na,lbwa,ubwa,error=error)
  a0_z=general(wb_z,error)
  call cpu_time(t0)
  bt_z=bt(wb_z,error)
  call cpu_time(t1)
  a1_z=general(bt_z,error)
  test_name = "Complex WB to BT, n=1;"
  call z_output_result(test_name,a0_z,a1_z,lbwa,bt_z%lbw,t0,t1,c*tol,error)

  na=2
  lbwa=1; ubwa=1
  wb_z=z_random_wb(na,lbwa,ubwa,error=error)
  a0_z=general(wb_z,error)
  call cpu_time(t0)
  bt_z=bt(wb_z,error)
  call cpu_time(t1)
  a1_z=general(bt_z,error)
  test_name = "Complex WB to BT, n=2;"
  call z_output_result(test_name,a0_z,a1_z,lbwa,bt_z%lbw,t0,t1,c*tol,error)
  
  na=3
  lbwa=1; ubwa=1
  wb_z=z_random_wb(na,lbwa,ubwa,error=error)
  a0_z=general(wb_z,error)
  call cpu_time(t0)
  bt_z=bt(wb_z,error)
  call cpu_time(t1)
  a1_z=general(bt_z,error)
  test_name = "Complex WB to BT, n=3;"
  call z_output_result(test_name,a0_z,a1_z,lbwa,bt_z%lbw,t0,t1,c*tol,error)

  na=4
  lbwa=2; ubwa=2
  wb_z=z_random_wb(na,lbwa,ubwa,error=error)
  a0_z=general(wb_z,error)
  call cpu_time(t0)
  bt_z=bt(wb_z,error)
  call cpu_time(t1)
  a1_z=general(bt_z,error)
  test_name = "Complex WB to BT, n=4;"
  call z_output_result(test_name,a0_z,a1_z,lbwa,bt_z%lbw,t0,t1,c*tol,error)
  print *
  
  na=40
  lbwa=5; ubwa=3
  bt_z=z_random_bt(na,lbwa,ubwa,error=error)
  a0_z=general(bt_z,error)
  call cpu_time(t0)
  wb_z=wb(bt_z,error)
  call cpu_time(t1)
  a1_z=general(wb_z,error)
  test_name = "Complex BT to WB, n=40;"
  call z_output_result(test_name,a0_z,a1_z,lbwa,wb_z%lbw,t0,t1,c*tol,error)

  na=1
  lbwa=0; ubwa=0
  bt_z=z_random_bt(na,lbwa,ubwa,error=error)
  a0_z=general(bt_z,error)
  call cpu_time(t0)
  wb_z=wb(bt_z,error)
  call cpu_time(t1)
  a1_z=general(wb_z,error)
  test_name = "Complex BT to WB, n=1;"
  call z_output_result(test_name,a0_z,a1_z,lbwa,wb_z%lbw,t0,t1,c*tol,error)

  na=2
  lbwa=1; ubwa=1
  bt_z=z_random_bt(na,lbwa,ubwa,error=error)
  a0_z=general(bt_z,error)
  call cpu_time(t0)
  wb_z=wb(bt_z,error)
  call cpu_time(t1)
  a1_z=general(wb_z,error)
  test_name = "Complex BT to WB, n=2;"
  call z_output_result(test_name,a0_z,a1_z,lbwa,wb_z%lbw,t0,t1,c*tol,error)

  na=3
  lbwa=1; ubwa=1
  bt_z=z_random_bt(na,lbwa,ubwa,error=error)
  a0_z=general(bt_z,error)
  call cpu_time(t0)
  wb_z=wb(bt_z,error)
  call cpu_time(t1)
  a1_z=general(wb_z,error)
  test_name = "Complex BT to WB, n=3;"
  call z_output_result(test_name,a0_z,a1_z,lbwa,wb_z%lbw,t0,t1,c*tol,error)

  na=4
  lbwa=2; ubwa=2
  bt_z=z_random_bt(na,lbwa,ubwa,error=error)
  a0_z=general(bt_z,error)
  call cpu_time(t0)
  wb_z=wb(bt_z,error)
  call cpu_time(t1)
  a1_z=general(wb_z,error)
  test_name = "Complex BT to WB, n=4;"
  call z_output_result(test_name,a0_z,a1_z,lbwa,wb_z%lbw,t0,t1,c*tol,error)

contains

    subroutine d_output_result(name,a0,a1,lbw0,lbw1,t0,t1,bnd,error)
    character(len=*) :: name
    real(kind=dp), dimension(:,:) :: a0, a1     
    real(kind=dp) :: bnd, t0, t1
    integer(kind=int32) :: lbw0, lbw1
    type(error_info) :: error

    real(kind=dp) :: berr
    character(len=10) :: test_result

    if (error%code > 0) then
       print *, "Calling error in test: ", name
    else
       berr = maxabs(a1-a0)
       if (lbw0==lbw1 .and. berr < bnd) then
          test_result="PASSED"
       else
          test_result="    FAILED"
       end if
       write (*,fmt) name, t1-t0, lbw1, berr, test_result
    end if
  end subroutine d_output_result

  subroutine z_output_result(name,a0,a1,lbw0,lbw1,t0,t1,bnd,error)
    character(len=*) :: name
    complex(kind=dp), dimension(:,:) :: a0, a1     
    real(kind=dp) :: bnd, t0, t1
    integer(kind=int32) :: lbw0, lbw1
    type(error_info) :: error

    real(kind=dp) :: berr
    character(len=10) :: test_result

    if (error%code > 0) then
       print *, "Calling error in test: ", name
    else
       berr = maxabs(a1-a0)
       if (lbw0==lbw1 .and. berr < bnd) then
          test_result="PASSED"
       else
          test_result="    FAILED"
       end if
       write (*,fmt) name, t1-t0, lbw1, berr, test_result
    end if
  end subroutine z_output_result

end program test_convert_wb_and_bt
