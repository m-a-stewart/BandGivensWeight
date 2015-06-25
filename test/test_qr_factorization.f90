program test_qr_factorization
  use mod_orrb
  implicit none

  character(len=40) :: test_name
  character(len=*), parameter :: fmt="(A40, 'Time: ',ES8.2,', lbw: ',I3,', ubw: ', I3,', error: ',ES8.2, ', ', A10)"

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
  swub_d=qr_of(bv_d,error)
  call cpu_time(t1)  
  a_d = general(swub_d%ub,error)
  a_d = swub_d%sw * a_d
  test_name = "Random Real QR Factorization"
  call d_output_result(test_name,a0_d,a_d,0, &
       swub_d%ub%lbw, min(ubwa+lbwa,na-1),swub_d%ub%ubw,t0,t1,c*tol,error)

  na=1; lbwa=0; ubwa=0
  bv_d=d_random_bv(na,lbwa,ubwa,error=error)
  a_d = general(bv_d,error)
  a0_d = a_d
  call cpu_time(t0)
  swub_d=qr_of(bv_d,error)
  call cpu_time(t1)  
  a_d = general(swub_d%ub,error)
  a_d = swub_d%sw * a_d
  test_name = "Random Real QR Factorization, n=1"
  call d_output_result(test_name,a0_d,a_d,0, &
       swub_d%ub%lbw, min(ubwa+lbwa,na-1), swub_d%ub%ubw,t0,t1,c*tol,error)

  na=2; lbwa=1; ubwa=1
  bv_d=d_random_bv(na,lbwa,ubwa,error=error)
  a_d = general(bv_d,error)
  a0_d = a_d
  call cpu_time(t0)
  swub_d=qr_of(bv_d,error)
  call cpu_time(t1)  
  a_d = general(swub_d%ub,error)
  a_d = swub_d%sw * a_d
  test_name = "Random Real QR Factorization, n=2"
  call d_output_result(test_name,a0_d,a_d,0, &
       swub_d%ub%lbw, min(ubwa+lbwa,na-1), swub_d%ub%ubw,t0,t1,c*tol,error)

  na=3; lbwa=1; ubwa=1
  bv_d=d_random_bv(na,lbwa,ubwa,error=error)
  a_d = general(bv_d,error)
  a0_d = a_d
  call cpu_time(t0)
  swub_d=qr_of(bv_d,error)
  call cpu_time(t1)  
  a_d = general(swub_d%ub,error)
  a_d = swub_d%sw * a_d
  test_name = "Random Real QR Factorization, n=3"
  call d_output_result(test_name,a0_d,a_d,0, &
       swub_d%ub%lbw, min(ubwa+lbwa,na-1), swub_d%ub%ubw,t0,t1,c*tol,error)

  na=4; lbwa=2; ubwa=2
  bv_d=d_random_bv(na,lbwa,ubwa,error=error)
  a_d = general(bv_d,error)
  a0_d = a_d
  call cpu_time(t0)
  swub_d=qr_of(bv_d,error)
  call cpu_time(t1)  
  a_d = general(swub_d%ub,error)
  a_d = swub_d%sw * a_d
  test_name = "Random Real QR Factorization, n=4"
  call d_output_result(test_name,a0_d,a_d,0, &
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
  swub_z=qr_of(bv_z,error)
  call cpu_time(t1)  
  a_z = general(swub_z%ub,error)
  a_z = swub_z%sw * a_z
  test_name = "Random Complex QR Factorization"
  call z_output_result(test_name,a0_z,a_z,0, &
       swub_z%ub%lbw, min(ubwa+lbwa,na-1),swub_z%ub%ubw,t0,t1,c*tol,error)

  na=1; lbwa=0; ubwa=0
  bv_z=z_random_bv(na,lbwa,ubwa,error=error)
  a_z = general(bv_z,error)
  a0_z = a_z
  call cpu_time(t0)
  swub_z=qr_of(bv_z,error)
  call cpu_time(t1)  
  a_z = general(swub_z%ub,error)
  a_z = swub_z%sw * a_z
  test_name = "Random Complex QR Factorization, n=1"
  call z_output_result(test_name,a0_z,a_z,0, &
       swub_z%ub%lbw, min(ubwa+lbwa,na-1), swub_z%ub%ubw,t0,t1,c*tol,error)

  na=2; lbwa=1; ubwa=1
  bv_z=z_random_bv(na,lbwa,ubwa,error=error)
  a_z = general(bv_z,error)
  a0_z = a_z
  call cpu_time(t0)
  swub_z=qr_of(bv_z,error)
  call cpu_time(t1)  
  a_z = general(swub_z%ub,error)
  a_z = swub_z%sw * a_z
  test_name = "Random Complex QR Factorization, n=2"
  call z_output_result(test_name,a0_z,a_z,0, &
       swub_z%ub%lbw, min(ubwa+lbwa,na-1), swub_z%ub%ubw,t0,t1,c*tol,error)

  na=3; lbwa=1; ubwa=1
  bv_z=z_random_bv(na,lbwa,ubwa,error=error)
  a_z = general(bv_z,error)
  a0_z = a_z
  call cpu_time(t0)
  swub_z=qr_of(bv_z,error)
  call cpu_time(t1)  
  a_z = general(swub_z%ub,error)
  a_z = swub_z%sw * a_z
  test_name = "Random Complex QR Factorization, n=3"
  call z_output_result(test_name,a0_z,a_z,0, &
       swub_z%ub%lbw, min(ubwa+lbwa,na-1), swub_z%ub%ubw,t0,t1,c*tol,error)

  na=4; lbwa=2; ubwa=2
  bv_z=z_random_bv(na,lbwa,ubwa,error=error)
  a_z = general(bv_z,error)
  a0_z = a_z
  call cpu_time(t0)
  swub_z=qr_of(bv_z,error)
  call cpu_time(t1)  
  a_z = general(swub_z%ub,error)
  a_z = swub_z%sw * a_z
  test_name = "Random Complex QR Factorization, n=4"
  call z_output_result(test_name,a0_z,a_z,0, &
       swub_z%ub%lbw, min(ubwa+lbwa,na-1), swub_z%ub%ubw,t0,t1,c*tol,error)

contains

  subroutine d_output_result(name,a0,a1,lbw0,lbw1,ubw0,ubw1,t0,t1,bnd,error)
    character(len=*) :: name
    real(kind=dp), dimension(:,:) :: a0, a1     
    real(kind=dp) :: bnd, t0, t1
    integer(kind=int32) :: lbw0, lbw1, ubw0, ubw1
    type(error_info) :: error

    real(kind=dp) :: berr
    character(len=10) :: test_result

    if (error%code > 0) then
       print *, "Calling error in test: ", name
    else
       berr = maxabs(a1-a0)
       if (lbw0==lbw1 .and. ubw0==ubw1 .and. berr < bnd) then
          test_result="PASSED"
       else
          test_result="    FAILED"
       end if
       write (*,fmt) name, t1-t0, lbw1, ubw1, berr, test_result
    end if
  end subroutine d_output_result

  subroutine z_output_result(name,a0,a1,lbw0,lbw1,ubw0,ubw1,t0,t1,bnd,error)
    character(len=*) :: name
    complex(kind=dp), dimension(:,:) :: a0, a1     
    real(kind=dp) :: bnd, t0, t1
    integer(kind=int32) :: lbw0, lbw1, ubw0, ubw1
    type(error_info) :: error

    real(kind=dp) :: berr
    character(len=10) :: test_result

    if (error%code > 0) then
       print *, "Calling error in test: ", name
    else
       berr = maxabs(a1-a0)
       if (lbw0==lbw1 .and. ubw0==ubw1 .and. berr < bnd) then
          test_result="PASSED"
       else
          test_result="    FAILED"
       end if
       write (*,fmt) name, t1-t0, lbw1, ubw1, berr, test_result
    end if
  end subroutine z_output_result
  
end program test_qr_factorization
