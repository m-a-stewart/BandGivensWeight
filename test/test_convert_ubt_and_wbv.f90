program test_convert_ubt_and_wbv
  use mod_orrb
  implicit none
  !
  ! Note: This does not cover large rank cases (e.g. ubw=n-1).
  !
  character(len=40) :: test_name
  character(len=*), parameter :: fmt="(A40, 'Time: ',ES8.2,', lbw: ',I3,', ubw: ', I3,', error: ',ES8.2, ', ', A10)"

  real(kind=dp) :: t0, t1
  integer(kind=int32) :: na, ubwa, lbwa
  type(error_info) :: error
  real(kind=dp), parameter :: tol=3e-15, c=2.0
  !
  real(kind=dp), dimension(:,:), allocatable :: a0_d, a1_d
  complex(kind=dp), dimension(:,:), allocatable :: a0_z, a1_z

  type(d_wbv), allocatable :: wbv_d
  type(z_wbv), allocatable :: wbv_z
  type(d_ubt), allocatable :: ubt_d
  type(z_ubt), allocatable :: ubt_z

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
  ubt_d=ubt_of(wbv_d,error)
  call cpu_time(t1)
  a1_d=general(ubt_d,error)
  test_name = "Real WBV to UBT, n=40;"
  call d_output_result(test_name,a0_d,a1_d,lbwa,ubt_d%lbw, &
       ubwa, ubt_d%ubw, t0,t1,c*tol,error)

  na=1
  lbwa=0; ubwa=0
  wbv_d=d_random_wbv(na,lbwa,ubwa,error=error)
  a0_d=general(wbv_d,error)
  call cpu_time(t0)
  ubt_d=ubt_of(wbv_d,error)
  call cpu_time(t1)
  a1_d=general(ubt_d,error)
  test_name = "Real WBV to UBT, n=1;"
  call d_output_result(test_name,a0_d,a1_d,lbwa,ubt_d%lbw, &
       ubwa, ubt_d%ubw, t0,t1,c*tol,error)


  na=2
  lbwa=1; ubwa=1
  wbv_d=d_random_wbv(na,lbwa,ubwa,error=error)
  a0_d=general(wbv_d,error)
  call cpu_time(t0)
  ubt_d=ubt_of(wbv_d,error)
  call cpu_time(t1)
  a1_d=general(ubt_d,error)
  test_name = "Real WBV to UBT, n=2;"
  call d_output_result(test_name,a0_d,a1_d,lbwa,ubt_d%lbw, &
       ubwa, ubt_d%ubw, t0,t1,c*tol,error)

  
  na=3
  lbwa=1; ubwa=1
  wbv_d=d_random_wbv(na,lbwa,ubwa,error=error)
  a0_d=general(wbv_d,error)
  call cpu_time(t0)
  ubt_d=ubt_of(wbv_d,error)
  call cpu_time(t1)
  a1_d=general(ubt_d,error)
  test_name = "Real WBV to UBT, n=3;"
  call d_output_result(test_name,a0_d,a1_d,lbwa,ubt_d%lbw, &
       ubwa, ubt_d%ubw, t0,t1,c*tol,error)


  na=4
  lbwa=2; ubwa=2
  wbv_d=d_random_wbv(na,lbwa,ubwa,error=error)
  a0_d=general(wbv_d,error)
  call cpu_time(t0)
  ubt_d=ubt_of(wbv_d,error)
  call cpu_time(t1)
  a1_d=general(ubt_d,error)
  test_name = "Real WBV to UBT, n=4;"
  call d_output_result(test_name,a0_d,a1_d,lbwa,ubt_d%lbw, &
       ubwa, ubt_d%ubw, t0,t1,c*tol,error)
  print *
  
  na=40
  lbwa=5; ubwa=3
  ubt_d=d_random_ubt(na,lbwa,ubwa,error=error)
  a0_d=general(ubt_d,error)
  call cpu_time(t0)
  wbv_d=wbv_of(ubt_d,error)
  call cpu_time(t1)
  a1_d=general(wbv_d,error)
  test_name = "Real UBT to WBV, n=40;"
  call d_output_result(test_name,a0_d,a1_d,lbwa,wbv_d%lbw, &
       ubwa, wbv_d%ubw, t0,t1,c*tol,error)

  na=1
  lbwa=0; ubwa=0
  ubt_d=d_random_ubt(na,lbwa,ubwa,error=error)
  a0_d=general(ubt_d,error)
  call cpu_time(t0)
  wbv_d=wbv_of(ubt_d,error)
  call cpu_time(t1)
  a1_d=general(wbv_d,error)
  test_name = "Real UBT to WBV, n=1;"
  call d_output_result(test_name,a0_d,a1_d,lbwa,wbv_d%lbw, &
       ubwa, wbv_d%ubw, t0,t1,c*tol,error)


  na=2
  lbwa=1; ubwa=1
  ubt_d=d_random_ubt(na,lbwa,ubwa,error=error)
  a0_d=general(ubt_d,error)
  call cpu_time(t0)
  wbv_d=wbv_of(ubt_d,error)
  call cpu_time(t1)
  a1_d=general(wbv_d,error)
  test_name = "Real UBT to WBV, n=2;"
  call d_output_result(test_name,a0_d,a1_d,lbwa,wbv_d%lbw, &
       ubwa, wbv_d%ubw, t0,t1,c*tol,error)


  na=3
  lbwa=1; ubwa=1
  ubt_d=d_random_ubt(na,lbwa,ubwa,error=error)
  a0_d=general(ubt_d,error)
  call cpu_time(t0)
  wbv_d=wbv_of(ubt_d,error)
  call cpu_time(t1)
  a1_d=general(wbv_d,error)
  test_name = "Real UBT to WBV, n=3;"
  call d_output_result(test_name,a0_d,a1_d,lbwa,wbv_d%lbw, &
       ubwa, wbv_d%ubw, t0,t1,c*tol,error)

  na=4
  lbwa=2; ubwa=2
  ubt_d=d_random_ubt(na,lbwa,ubwa,error=error)
  a0_d=general(ubt_d,error)
  call cpu_time(t0)
  wbv_d=wbv_of(ubt_d,error)
  call cpu_time(t1)
  a1_d=general(wbv_d,error)
  test_name = "Real UBT to WBV, n=4;"
  call d_output_result(test_name,a0_d,a1_d,lbwa,wbv_d%lbw, &
       ubwa, wbv_d%ubw, t0,t1,c*tol,error)

  print *
  print *, "--------------------------------"
  print *
  print *, "Complex WBV and UBT Conversion Tests:"
  print *

  na=40
  lbwa=5; ubwa=3
  wbv_z=z_random_wbv(na,lbwa,ubwa,error=error)
  a0_z=general(wbv_z,error)
  call cpu_time(t0)
  ubt_z=ubt_of(wbv_z,error)
  call cpu_time(t1)
  a1_z=general(ubt_z,error)
  test_name = "Complex WBV to UBT, n=40;"
  call z_output_result(test_name,a0_z,a1_z,lbwa,ubt_z%lbw, &
       ubwa, ubt_z%ubw, t0,t1,c*tol,error)

  na=1
  lbwa=0; ubwa=0
  wbv_z=z_random_wbv(na,lbwa,ubwa,error=error)
  a0_z=general(wbv_z,error)
  call cpu_time(t0)
  ubt_z=ubt_of(wbv_z,error)
  call cpu_time(t1)
  a1_z=general(ubt_z,error)
  test_name = "Complex WBV to UBT, n=1;"
  call z_output_result(test_name,a0_z,a1_z,lbwa,ubt_z%lbw, &
       ubwa, ubt_z%ubw, t0,t1,c*tol,error)

  na=2
  lbwa=1; ubwa=1
  wbv_z=z_random_wbv(na,lbwa,ubwa,error=error)
  a0_z=general(wbv_z,error)
  call cpu_time(t0)
  ubt_z=ubt_of(wbv_z,error)
  call cpu_time(t1)
  a1_z=general(ubt_z,error)
  test_name = "Complex WBV to UBT, n=2;"
  call z_output_result(test_name,a0_z,a1_z,lbwa,ubt_z%lbw, &
       ubwa, ubt_z%ubw, t0,t1,c*tol,error)
  
  na=3
  lbwa=1; ubwa=1
  wbv_z=z_random_wbv(na,lbwa,ubwa,error=error)
  a0_z=general(wbv_z,error)
  call cpu_time(t0)
  ubt_z=ubt_of(wbv_z,error)
  call cpu_time(t1)
  a1_z=general(ubt_z,error)
  test_name = "Complex WBV to UBT, n=3;"
  call z_output_result(test_name,a0_z,a1_z,lbwa,ubt_z%lbw, &
       ubwa, ubt_z%ubw, t0,t1,c*tol,error)

  na=4
  lbwa=2; ubwa=2
  wbv_z=z_random_wbv(na,lbwa,ubwa,error=error)
  a0_z=general(wbv_z,error)
  call cpu_time(t0)
  ubt_z=ubt_of(wbv_z,error)
  call cpu_time(t1)
  a1_z=general(ubt_z,error)
  test_name = "Complex WBV to UBT, n=4;"
  call z_output_result(test_name,a0_z,a1_z,lbwa,ubt_z%lbw, &
       ubwa, ubt_z%ubw, t0,t1,c*tol,error)
  print *
  
  na=40
  lbwa=5; ubwa=3
  ubt_z=z_random_ubt(na,lbwa,ubwa,error=error)
  a0_z=general(ubt_z,error)
  call cpu_time(t0)
  wbv_z=wbv_of(ubt_z,error)
  call cpu_time(t1)
  a1_z=general(wbv_z,error)
  test_name = "Complex UBT to WBV, n=40;"
  call z_output_result(test_name,a0_z,a1_z,lbwa,wbv_z%lbw, &
       ubwa,wbv_z%ubw,t0,t1,c*tol,error)

  na=1
  lbwa=0; ubwa=0
  ubt_z=z_random_ubt(na,lbwa,ubwa,error=error)
  a0_z=general(ubt_z,error)
  call cpu_time(t0)
  wbv_z=wbv_of(ubt_z,error)
  call cpu_time(t1)
  a1_z=general(wbv_z,error)
  test_name = "Complex UBT to WBV, n=1;"
  call z_output_result(test_name,a0_z,a1_z,lbwa,wbv_z%lbw, &
       ubwa,wbv_z%ubw,t0,t1,c*tol,error)

  na=2
  lbwa=1; ubwa=1
  ubt_z=z_random_ubt(na,lbwa,ubwa,error=error)
  a0_z=general(ubt_z,error)
  call cpu_time(t0)
  wbv_z=wbv_of(ubt_z,error)
  call cpu_time(t1)
  a1_z=general(wbv_z,error)
  test_name = "Complex UBT to WBV, n=2;"
  call z_output_result(test_name,a0_z,a1_z,lbwa,wbv_z%lbw, &
       ubwa,wbv_z%ubw,t0,t1,c*tol,error)

  na=3
  lbwa=1; ubwa=1
  ubt_z=z_random_ubt(na,lbwa,ubwa,error=error)
  a0_z=general(ubt_z,error)
  call cpu_time(t0)
  wbv_z=wbv_of(ubt_z,error)
  call cpu_time(t1)
  a1_z=general(wbv_z,error)
  test_name = "Complex UBT to WBV, n=3;"
  call z_output_result(test_name,a0_z,a1_z,lbwa,wbv_z%lbw, &
       ubwa,wbv_z%ubw,t0,t1,c*tol,error)

  na=4
  lbwa=2; ubwa=2
  ubt_z=z_random_ubt(na,lbwa,ubwa,error=error)
  a0_z=general(ubt_z,error)
  call cpu_time(t0)
  wbv_z=wbv_of(ubt_z,error)
  call cpu_time(t1)
  a1_z=general(wbv_z,error)
  test_name = "Complex UBT to WBV, n=4;"
  call z_output_result(test_name,a0_z,a1_z,lbwa,wbv_z%lbw, &
       ubwa,wbv_z%ubw,t0,t1,c*tol,error)

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


end program test_convert_ubt_and_wbv
