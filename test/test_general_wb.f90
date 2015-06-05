program test_general_wb
  use mod_orrb
  implicit none

  character(len=40) :: test_name
  character(len=*), parameter :: fmt="(A40, 'Time: ',ES8.2,', lbw: ',I3,', error: ',ES8.2, ', ', A10)"
  
  real(kind=dp) :: t0, t1
  type(error_info) :: error
  integer(kind=int32) :: na, lbwa, ubwa, j
  real(kind=dp), parameter :: tol=1e-14, c=1.5
  !
  real(kind=dp), dimension(:,:), allocatable :: a_d, a0_d, a1_d
  complex(kind=dp), dimension(:,:), allocatable :: a_z, a0_z, a1_z
  type(d_wb), allocatable :: wb_d
  type(z_wb), allocatable :: wb_z

  call initialize_errors

  ! test one
  print *
  print *, "--------------------------------"
  print *
  print *, "Real WB Decomposition Tests"
  print *

  na=40
  lbwa=5; ubwa=3
  wb_d=d_random_wb(na,lbwa,ubwa,error=error)
  a_d=general(wb_d,error)
  a1_d=a_d
  call cpu_time(t0)
  wb_d=wb(a_d,ubwa,lbwa+1,ubwa,tol,error)
  call cpu_time(t1)
  a0_d = general(wb_d,error)
  test_name="Random Real WB;"  
  call d_output_result(test_name,a0_d,a1_d,lbwa,wb_d%lbw,t0,t1,c*tol,error)

  na=50
  lbwa=13; ubwa=3
  wb_d=d_random_wb(na,[ (lbwa, j=1,lbwa-1), (lbwa-1, j=lbwa,na) ], &
      [ (ubwa, j=1,na) ], error=error )  
  a_d=general(wb_d,error)
  a1_d=a_d
  wb_d=wb(a_d,ubwa,lbwa+1,ubwa,tol,error)
  call cpu_time(t0)
  a0_d = general(wb_d,error)
  call cpu_time(t1)
  test_name="Random Real Square Termination WB;"  
  call d_output_result(test_name,a0_d,a1_d,lbwa,wb_d%lbw,t0,t1,c*tol,error)
  deallocate(wb_d)

  !
  ! Complex WB test
  !
  print *
  print *, "--------------------------------"
  print *
  print *, "Complex WB Decomposition Tests"
  print *

  na=40
  lbwa=5; ubwa=3
  wb_z=z_random_wb(na,lbwa,ubwa,error=error)
  a_z=general(wb_z,error)
  a1_z=a_z
  call cpu_time(t0)
  wb_z=wb(a_z, ubwa, lbwa+1, ubwa, tol,error)
  call cpu_time(t1)
  a0_z=general(wb_z,error)
  test_name="Random Complex WB;"  
  call z_output_result(test_name,a0_z,a1_z,lbwa,wb_z%lbw,t0,t1,c*tol,error)

  na=50
  lbwa=13; ubwa=3
  wb_z=z_random_wb(na,[ (lbwa, j=1,lbwa-1), (lbwa-1, j=lbwa,na) ], &
      [ (ubwa, j=1,na) ], error=error )  
  a_z=general(wb_z,error)
  a1_z=a_z
  call cpu_time(t0)
  wb_z=wb(a_z,ubwa,lbwa+1,ubwa,tol,error)
  call cpu_time(t1)
  a0_z = general(wb_z,error)
  test_name="Random Complex Square Termination WB;"  
  call z_output_result(test_name,a0_z,a1_z,lbwa,wb_z%lbw,t0,t1,c*tol,error)
  print *

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
  
end program test_general_wb
