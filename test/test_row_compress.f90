program test_row_compress
  use mod_orrb
  implicit none

  character(len=40) :: test_name
  character(len=*), parameter :: fmt="(A40, 'Time: ',ES8.2,', lbw: ',I3,', ubw: ', I3,', error: ',ES8.2, ', ', A10)"
  
  real(kind=dp) :: t0, t1
  integer(kind=int32) :: na, lbwa, ubwa
  type(error_info) :: error
  real(kind=dp), parameter :: tol=1e-15, c=4.0
  real(kind=dp), dimension(:,:), allocatable :: a_d, a0_d
  complex(kind=dp), dimension(:,:), allocatable :: a_z, a0_z
  !
  type(d_ubt), allocatable :: ubt_d
  type(z_ubt), allocatable :: ubt_z
  type(d_rc), allocatable :: swbv_d
  type(z_rc), allocatable :: swbv_z

  call initialize_errors

  print *
  print *, "--------------------------------"
  print *
  print *, "Real Row Compression Tests:"
  print *

  na=40; lbwa=7; ubwa=5
  ubt_d=d_random_ubt(na,lbwa,ubwa,error=error)
  a_d = general(ubt_d,error)
  a0_d = a_d
  call cpu_time(t0)
  swbv_d=rc(ubt_d,error)
  call cpu_time(t1)  
  a_d = general(swbv_d%bv,error)
  a_d = swbv_d%sw * a_d
  test_name = "Random Real Row Compression, n=40"
  call d_output_result(test_name,a0_d,a_d,lbwa, &
       swbv_d%bv%lbw, min(ubwa+lbwa,na-1),swbv_d%bv%ubw,t0,t1,c*tol,error)

  na=1; lbwa=0; ubwa=0
  ubt_d=d_random_ubt(na,lbwa,ubwa,error=error)
  a_d = general(ubt_d,error)
  a0_d = a_d
  call cpu_time(t0)
  swbv_d=rc(ubt_d,error)
  call cpu_time(t1)  
  a_d = general(swbv_d%bv,error)
  a_d = swbv_d%sw * a_d
  test_name = "Random Real Row Compression, n=1"
  call d_output_result(test_name,a0_d,a_d,lbwa, &
       swbv_d%bv%lbw, min(ubwa+lbwa,na-1),swbv_d%bv%ubw,t0,t1,c*tol,error)

  na=2; lbwa=1; ubwa=1
  ubt_d=d_random_ubt(na,lbwa,ubwa,error=error)
  a_d = general(ubt_d,error)
  a0_d = a_d
  call cpu_time(t0)
  swbv_d=rc(ubt_d,error)
  call cpu_time(t1)  
  a_d = general(swbv_d%bv,error)
  a_d = swbv_d%sw * a_d
  test_name = "Random Real Row Compression, n=2"
  call d_output_result(test_name,a0_d,a_d,lbwa, &
       swbv_d%bv%lbw, min(ubwa+lbwa,na-1),swbv_d%bv%ubw,t0,t1,c*tol,error)

  na=3; lbwa=1; ubwa=1
  ubt_d=d_random_ubt(na,lbwa,ubwa,error=error)
  a_d = general(ubt_d,error)
  a0_d = a_d
  call cpu_time(t0)
  swbv_d=rc(ubt_d,error)
  call cpu_time(t1)  
  a_d = general(swbv_d%bv,error)
  a_d = swbv_d%sw * a_d
  test_name = "Random Real Row Compression, n=3"
  call d_output_result(test_name,a0_d,a_d,lbwa, &
       swbv_d%bv%lbw, min(ubwa+lbwa,na-1),swbv_d%bv%ubw,t0,t1,c*tol,error)

  na=4; lbwa=2; ubwa=2
  ubt_d=d_random_ubt(na,lbwa,ubwa,error=error)
  a_d = general(ubt_d,error)
  a0_d = a_d
  call cpu_time(t0)
  swbv_d=rc(ubt_d,error)
  call cpu_time(t1)  
  a_d = general(swbv_d%bv,error)
  a_d = swbv_d%sw * a_d
  test_name = "Random Real Row Compression, n=4"
  call d_output_result(test_name,a0_d,a_d,lbwa, &
       swbv_d%bv%lbw, min(ubwa+lbwa,na-1),swbv_d%bv%ubw,t0,t1,c*tol,error)
  
  print *
  print *, "--------------------------------"
  print *
  print *, "Complex Row Compression Tests:"
  print *

  na=40; lbwa=7; ubwa=5
  ubt_z=z_random_ubt(na,lbwa,ubwa,error=error)
  a_z = general(ubt_z,error)
  a0_z = a_z
  call cpu_time(t0)
  swbv_z=rc(ubt_z,error)
  call cpu_time(t1)  
  a_z = general(swbv_z%bv,error)
  a_z = swbv_z%sw * a_z
  test_name = "Random Complex Row Compression, n=40"
  call z_output_result(test_name,a0_z,a_z,lbwa, &
       swbv_z%bv%lbw, min(ubwa+lbwa,na-1),swbv_z%bv%ubw,t0,t1,c*tol,error)

  na=1; lbwa=0; ubwa=0
  ubt_z=z_random_ubt(na,lbwa,ubwa,error=error)
  a_z = general(ubt_z,error)
  a0_z = a_z
  call cpu_time(t0)
  swbv_z=rc(ubt_z,error)
  call cpu_time(t1)  
  a_z = general(swbv_z%bv,error)
  a_z = swbv_z%sw * a_z
  test_name = "Random Complex Row Compression, n=1"
  call z_output_result(test_name,a0_z,a_z,lbwa, &
       swbv_z%bv%lbw, min(ubwa+lbwa,na-1),swbv_z%bv%ubw,t0,t1,c*tol,error)

  na=2; lbwa=1; ubwa=1
  ubt_z=z_random_ubt(na,lbwa,ubwa,error=error)
  a_z = general(ubt_z,error)
  a0_z = a_z
  call cpu_time(t0)
  swbv_z=rc(ubt_z,error)
  call cpu_time(t1)  
  a_z = general(swbv_z%bv,error)
  a_z = swbv_z%sw * a_z
  test_name = "Random Complex Row Compression, n=2"
  call z_output_result(test_name,a0_z,a_z,lbwa, &
       swbv_z%bv%lbw, min(ubwa+lbwa,na-1),swbv_z%bv%ubw,t0,t1,c*tol,error)

  na=3; lbwa=1; ubwa=1
  ubt_z=z_random_ubt(na,lbwa,ubwa,error=error)
  a_z = general(ubt_z,error)
  a0_z = a_z
  call cpu_time(t0)
  swbv_z=rc(ubt_z,error)
  call cpu_time(t1)  
  a_z = general(swbv_z%bv,error)
  a_z = swbv_z%sw * a_z
  test_name = "Random Complex Row Compression, n=3"
  call z_output_result(test_name,a0_z,a_z,lbwa, &
       swbv_z%bv%lbw, min(ubwa+lbwa,na-1),swbv_z%bv%ubw,t0,t1,c*tol,error)

  na=4; lbwa=2; ubwa=2
  ubt_z=z_random_ubt(na,lbwa,ubwa,error=error)
  a_z = general(ubt_z,error)
  a0_z = a_z
  call cpu_time(t0)
  swbv_z=rc(ubt_z,error)
  call cpu_time(t1)  
  a_z = general(swbv_z%bv,error)
  a_z = swbv_z%sw * a_z
  test_name = "Random Complex Row Compression, n=4"
  call z_output_result(test_name,a0_z,a_z,lbwa, &
       swbv_z%bv%lbw, min(ubwa+lbwa,na-1),swbv_z%bv%ubw,t0,t1,c*tol,error)

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

end program test_row_compress
