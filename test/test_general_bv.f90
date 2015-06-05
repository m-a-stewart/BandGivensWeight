program test_general_bv
  use mod_orrb
  implicit none

  character(len=40) :: test_name
  character(len=*), parameter :: fmt="(A40, 'Time: ',ES8.2,', ubw: ',I3,', error: ',ES8.2, ', ', A10)"

  real(kind=dp) :: t0, t1
  type(error_info) :: error
  integer(kind=int32) :: na, lbwa, ubwa, j
  real(kind=dp), parameter :: tol=1e-14, c=1.5
  !
  real(kind=dp), dimension(:,:), allocatable :: a_d, a0_d, a1_d
  complex(kind=dp), dimension(:,:), allocatable :: a_z, a0_z, a1_z
  type(d_bv), allocatable :: bv_d
  type(z_bv), allocatable :: bv_z

  call initialize_errors

  print *
  print *, "--------------------------------"
  print *
  print *, "Real BV Decomposition Tests"
  print *

  na=40
  lbwa=3; ubwa=5
  bv_d=d_random_bv(na,lbwa,ubwa,error=error)
  a_d=general(bv_d,error)
  a1_d=a_d
  call cpu_time(t0)
  bv_d=bv(a_d,lbwa,lbwa,ubwa+1,tol,error)
  call cpu_time(t1)
  a0_d=general(bv_d,error)
  test_name="Random Real BV;"  
  call d_output_result(test_name,a0_d,a1_d,ubwa,bv_d%ubw,t0,t1,c*tol,error)

  na=40
  lbwa=3; ubwa=0
  bv_d=d_random_bv(na,lbwa,ubwa,error=error)
  a_d=general(bv_d,error)
  a1_d=a_d
  call cpu_time(t0)
  bv_d=bv(a_d,lbwa,lbwa,ubwa+1,tol,error)
  call cpu_time(t1)
  a0_d=general(bv_d,error)
  test_name="Random Real BV, ubwa=0;"  
  call d_output_result(test_name,a0_d,a1_d,ubwa,bv_d%ubw,t0,t1,c*tol,error)

  na=40
  lbwa=3; ubwa=1
  bv_d=d_random_bv(na,lbwa,ubwa,error=error)
  a_d=general(bv_d,error)
  a1_d=a_d
  call cpu_time(t0)
  bv_d=bv(a_d,lbwa,lbwa,ubwa+1,tol,error)
  call cpu_time(t1)
  a0_d=general(bv_d,error)
  test_name="Random Real BV, ubwa=1;"
  call d_output_result(test_name,a0_d,a1_d,ubwa,bv_d%ubw,t0,t1,c*tol,error)

  na=1
  lbwa=0; ubwa=0
  bv_d=d_random_bv(na,lbwa,ubwa,error=error)
  a_d=general(bv_d,error)
  a1_d=a_d
  call cpu_time(t0)
  bv_d=bv(a_d,lbwa,lbwa,ubwa+1,tol,error)
  call cpu_time(t1)
  a0_d=general(bv_d,error)
  test_name="Random Real BV, na=1;"  
  call d_output_result(test_name,a0_d,a1_d,ubwa,bv_d%ubw,t0,t1,c*tol,error)

  na=2
  lbwa=1; ubwa=1
  bv_d=d_random_bv(na,lbwa,ubwa,error=error)
  a_d=general(bv_d,error)
  a1_d=a_d
  call cpu_time(t0)
  bv_d=bv(a_d,lbwa,lbwa,ubwa+1,tol,error)
  call cpu_time(t1)
  a0_d=general(bv_d,error)
  test_name="Random Real BV, na=2;"
  call d_output_result(test_name,a0_d,a1_d,ubwa,bv_d%ubw,t0,t1,c*tol,error)

  na=3
  lbwa=1; ubwa=1
  bv_d=d_random_bv(na,lbwa,ubwa,error=error)
  a_d=general(bv_d,error)
  a1_d=a_d
  call cpu_time(t0)
  bv_d=bv(a_d,lbwa,lbwa,ubwa+1,tol,error)
  call cpu_time(t1)
  a0_d=general(bv_d,error)
  test_name="Random Real BV, na=3;"
  call d_output_result(test_name,a0_d,a1_d,ubwa,bv_d%ubw,t0,t1,c*tol,error)

  na=4
  lbwa=2; ubwa=2
  bv_d=d_random_bv(na,lbwa,ubwa,error=error)
  a_d=general(bv_d,error)
  a1_d=a_d
  call cpu_time(t0)
  bv_d=bv(a_d,lbwa,lbwa,ubwa+1,tol,error)
  call cpu_time(t1)
  a0_d=general(bv_d,error)
  test_name="Random Real BV, na=4;"
  call d_output_result(test_name,a0_d,a1_d,ubwa,bv_d%ubw,t0,t1,c*tol,error)
  
  na=50
  lbwa=3; ubwa=13
  bv_d=d_random_bv(na,[ (lbwa, j=1,na) ], &
       [ (ubwa, j=1,ubwa-1), (ubwa-1, j=ubwa,na) ], error=error )
  a_d=general(bv_d, error)
  a1_d=a_d
  call cpu_time(t0)
  bv_d=bv(a_d,lbwa,lbwa,ubwa+1,tol, error)
  call cpu_time(t1)
  a0_d=general(bv_d, error)
  test_name="Random Real Square Termination BV;"  
  call d_output_result(test_name,a0_d,a1_d,ubwa,bv_d%ubw,t0,t1,c*tol,error)

  print *
  print *, "--------------------------------"
  print *
  print *, "Complex BV Decomposition Tests"
  print *
  na=40
  lbwa=3; ubwa=5
  bv_z=z_random_bv(na,lbwa,ubwa,error=error)
  a_z=general(bv_z,error)
  a1_z=a_z
  call cpu_time(t0)
  bv_z=bv(a_z,lbwa,lbwa,ubwa+1,tol,error)
  call cpu_time(t1)
  a0_z=general(bv_z)
  test_name="Random Complex BV;"  
  call z_output_result(test_name,a0_z,a1_z,ubwa,bv_z%ubw,t0,t1,c*tol,error)

  na=40
  lbwa=3; ubwa=0
  bv_z=z_random_bv(na,lbwa,ubwa,error=error)
  a_z=general(bv_z,error)
  a1_z=a_z
  call cpu_time(t0)
  bv_z=bv(a_z,lbwa,lbwa,ubwa+1,tol,error)
  call cpu_time(t1)
  a0_z=general(bv_z)
  test_name="Random Complex BV, ubwa=0;"  
  call z_output_result(test_name,a0_z,a1_z,ubwa,bv_z%ubw,t0,t1,c*tol,error)

  na=40
  lbwa=3; ubwa=1
  bv_z=z_random_bv(na,lbwa,ubwa,error=error)
  a_z=general(bv_z,error)
  a1_z=a_z
  call cpu_time(t0)
  bv_z=bv(a_z,lbwa,lbwa,ubwa+1,tol,error)
  call cpu_time(t1)
  a0_z=general(bv_z)
  test_name="Random Complex BV, ubwa=1;"
  call z_output_result(test_name,a0_z,a1_z,ubwa,bv_z%ubw,t0,t1,c*tol,error)

  na=1
  lbwa=0; ubwa=0
  bv_z=z_random_bv(na,lbwa,ubwa,error=error)
  a_z=general(bv_z,error)
  a1_z=a_z
  call cpu_time(t0)
  bv_z=bv(a_z,lbwa,lbwa,ubwa+1,tol,error)
  call cpu_time(t1)
  a0_z=general(bv_z,error)
  test_name="Random Complex BV, na=1;"  
  call z_output_result(test_name,a0_z,a1_z,ubwa,bv_z%ubw,t0,t1,c*tol,error)

  na=2
  lbwa=1; ubwa=1
  bv_z=z_random_bv(na,lbwa,ubwa,error=error)
  a_z=general(bv_z,error)
  a1_z=a_z
  call cpu_time(t0)
  bv_z=bv(a_z,lbwa,lbwa,ubwa+1,tol,error)
  call cpu_time(t1)
  a0_z=general(bv_z,error)
  test_name="Random Complex BV, na=2;"
  call z_output_result(test_name,a0_z,a1_z,ubwa,bv_z%ubw,t0,t1,c*tol,error)

  na=3
  lbwa=1; ubwa=1
  bv_z=z_random_bv(na,lbwa,ubwa,error=error)
  a_z=general(bv_z,error)
  a1_z=a_z
  call cpu_time(t0)
  bv_z=bv(a_z,lbwa,lbwa,ubwa+1,tol,error)
  call cpu_time(t1)
  a0_z=general(bv_z,error)
  test_name="Random Complex BV, na=3;"
  call z_output_result(test_name,a0_z,a1_z,ubwa,bv_z%ubw,t0,t1,c*tol,error)

  na=4
  lbwa=2; ubwa=2
  bv_z=z_random_bv(na,lbwa,ubwa,error=error)
  a_z=general(bv_z,error)
  a1_z=a_z
  call cpu_time(t0)
  bv_z=bv(a_z,lbwa,lbwa,ubwa+1,tol,error)
  call cpu_time(t1)
  a0_z=general(bv_z,error)
  test_name="Random Complex BV, na=4;"
  call z_output_result(test_name,a0_z,a1_z,ubwa,bv_z%ubw,t0,t1,c*tol,error)
  
  na=50
  lbwa=3; ubwa=13
  bv_z=z_random_bv(na,[ (lbwa, j=1,na) ], &
       [ (ubwa, j=1,ubwa-1), (ubwa-1, j=ubwa,na) ], error=error )
  a_z=general(bv_z, error)
  a1_z=a_z
  call cpu_time(t0)
  bv_z=bv(a_z,lbwa,lbwa,ubwa+1,tol,error)
  call cpu_time(t1)
  a0_z=general(bv_z,error)
  test_name="Random Complex Square Termination BV;"  
  call z_output_result(test_name,a0_z,a1_z,ubwa,bv_z%ubw,t0,t1,c*tol,error)

contains

  subroutine d_output_result(name,a0,a1,ubw0,ubw1,t0,t1,bnd,error)
    character(len=*) :: name
    real(kind=dp), dimension(:,:) :: a0, a1     
    real(kind=dp) :: bnd, t0, t1
    integer(kind=int32) :: ubw0, ubw1
    type(error_info) :: error

    real(kind=dp) :: berr
    character(len=10) :: test_result

    if (error%code > 0) then
       print *, "Calling error in test: ", name
    else
       berr = maxabs(a1-a0)
       if (ubw0==ubw1 .and. berr < bnd) then
          test_result="PASSED"
       else
          test_result="    FAILED"
       end if
       write (*,fmt) name, t1-t0, ubw1, berr, test_result
    end if
  end subroutine d_output_result

  subroutine z_output_result(name,a0,a1,ubw0,ubw1,t0,t1,bnd,error)
    character(len=*) :: name
    complex(kind=dp), dimension(:,:) :: a0, a1     
    real(kind=dp) :: bnd, t0, t1
    integer(kind=int32) :: ubw0, ubw1
    type(error_info) :: error

    real(kind=dp) :: berr
    character(len=10) :: test_result

    if (error%code > 0) then
       print *, "Calling error in test: ", name
    else
       berr = maxabs(a1-a0)
       if (ubw0==ubw1 .and. berr < bnd) then
          test_result="PASSED"
       else
          test_result="    FAILED"
       end if
       write (*,fmt) name, t1-t0, ubw1, berr, test_result
    end if
  end subroutine z_output_result
  
end program test_general_bv
