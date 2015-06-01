module mod_test_data
  use mod_error_id
  use mod_utility
  implicit none

  character(len=40) :: test_name
  character(len=*), parameter :: fmt_upper="(A40, 'Time: ',ES8.2,', ubw: ',I3,', error: ',ES8.2, ', ', A10)"
  character(len=*), parameter :: fmt_lower="(A40, 'Time: ',ES8.2,', lbw: ',I3,', error: ',ES8.2, ', ', A10)"
  character(len=*), parameter :: fmt_lower_upper="(A40, 'Time: ',ES8.2,', lbw: ',I3,', ubw: ', I3,', error: ',ES8.2, ', ', A10)"
  character(len=*), parameter :: fmt_qr="(A25, 'Time: ',ES8.2,', ubw: ',I3,', errors: ',ES8.2, ', ', ES8.2, ', ', ES8.2, ', ', A10)"

contains

  subroutine d_output_result_upper(name,a0,a1,ubw0,ubw1,t0,t1,bnd,error)
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
       write (*,fmt_upper) name, t1-t0, ubw1, berr, test_result
    end if
  end subroutine d_output_result_upper

  subroutine z_output_result_upper(name,a0,a1,ubw0,ubw1,t0,t1,bnd,error)
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
       write (*,fmt_upper) name, t1-t0, ubw1, berr, test_result
    end if
  end subroutine z_output_result_upper

  subroutine d_output_result_lower(name,a0,a1,lbw0,lbw1,t0,t1,bnd,error)
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
       write (*,fmt_lower) name, t1-t0, lbw1, berr, test_result
    end if
  end subroutine d_output_result_lower

  subroutine z_output_result_lower(name,a0,a1,lbw0,lbw1,t0,t1,bnd,error)
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
       write (*,fmt_lower) name, t1-t0, lbw1, berr, test_result
    end if
  end subroutine z_output_result_lower

  subroutine d_output_result_lower_upper(name,a0,a1,lbw0,lbw1,ubw0,ubw1,t0,t1,bnd,error)
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
       write (*,fmt_lower_upper) name, t1-t0, lbw1, ubw1, berr, test_result
    end if
  end subroutine d_output_result_lower_upper

  subroutine z_output_result_lower_upper(name,a0,a1,lbw0,lbw1,ubw0,ubw1,t0,t1,bnd,error)
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
       write (*,fmt_lower_upper) name, t1-t0, lbw1, ubw1, berr, test_result
    end if
  end subroutine z_output_result_lower_upper

end module mod_test_data
