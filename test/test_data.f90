module test_data
use misc
implicit none

integer, parameter :: n=50, rmax=13, ubwmax=rmax+1, lbw=2, lbwmax=10
real(kind=dp), parameter :: tol=1e-14, tol1=1e-14, tol2=1e-10
character(len=40) :: test_name
character(len=*), parameter :: fmt="(A40, 'Time: ',ES8.2,', ubw: ',I3,', error: ',ES8.2, ', ', A10)"
character(len=*), parameter :: fmt_qr="(A25, 'Time: ',ES8.2,', ubw: ',I3,', errors: ',ES8.2, ', ', ES8.2, ', ', ES8.2, ', ', A10)"
contains
  
  subroutine d_assemble_a(a,u,v,d,lbw)
    real(kind=dp), dimension(:,:), intent(out) :: a
    real(kind=dp), dimension(:,:), intent(in) :: u, v
    real(kind=dp), dimension(:), intent(in) :: d
    integer(kind=int32), intent(in) :: lbw

    integer(kind=int32) :: j,k,n
    n=size(a,1)
    a=matmul(u,v)
    do k=1,n-lbw-1
       do j=k+lbw+1,n
          a(j,k)=0.0_dp
       end do
    end do
    do j=1,n
       a(j,j)=d(j)
    end do
  end subroutine d_assemble_a

  subroutine c_assemble_a(a,u,v,d,lbw)
    complex(kind=dp), dimension(:,:), intent(out) :: a
    complex(kind=dp), dimension(:,:), intent(in) :: u, v
    complex(kind=dp), dimension(:), intent(in) :: d
    integer(kind=int32), intent(in) :: lbw

    integer(kind=int32) :: j,k,n
    n=size(a,1)
    a=matmul(u,v)
    do k=1,n-lbw-1
       do j=k+lbw+1,n
          a(j,k)=(0.0_dp, 0.0_dp)
       end do
    end do
    do j=1,n
       a(j,j)=d(j)
    end do
  end subroutine c_assemble_a

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

  subroutine c_output_result(name,a0,a1,ubw0,ubw1,t0,t1,bnd,error)
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
  end subroutine c_output_result

  subroutine c_output_result_qr(name,a0,a1,q,ubw0,ubw1,t0,t1,bnd,error)
    character(len=*) :: name
    complex(kind=dp), dimension(:,:) :: a0, a1, q
    real(kind=dp) :: bnd, t0, t1
    integer(kind=int32) :: ubw0, ubw1
    type(error_info) :: error

    real(kind=dp) :: berr, qerr, subd_err
    character(len=10) :: test_result
    complex(kind=dp), dimension(size(a0,1),size(a0,2)) :: f
    integer(kind=int32) :: j
    if (error%code > 0) then
       print *, "Calling error in test: ", name
    else
       berr = maxabs(matmul(matmul(transpose(conjg(q)),a0),q) - a1)
       f = matmul(transpose(conjg(q)),q)
       do j=1,size(q,2)
          f(j,j)=f(j,j)-(1.0_dp,0.0_dp)
       end do
       qerr=maxabs(f)
       subd_err=0.0_dp
       do j=1,size(a0,1)-1
          if (abs(a1(j+1,j)) > subd_err) then
             subd_err=abs(a1(j+1,j))
          end if
       end do
       if (ubw0<=ubw1 .and. berr < bnd .and. qerr < bnd .and. subd_err==0.0_dp) then
          test_result="PASSED"
       else
          test_result="    FAILED"
       end if
       write (*,fmt_qr) name, t1-t0, ubw1, berr, qerr, subd_err, test_result
    end if
  end subroutine c_output_result_qr

end module test_data
