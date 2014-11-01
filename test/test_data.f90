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

  subroutine d_assemble_upper(a,u,v,d,lbw)
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
  end subroutine d_assemble_upper

  subroutine c_assemble_upper(a,u,v,d,lbw)
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
  end subroutine c_assemble_upper

  subroutine d_assemble_lower(a,u,v,d,ubw)
    real(kind=dp), dimension(:,:), intent(out) :: a
    real(kind=dp), dimension(:,:), intent(in) :: u, v
    real(kind=dp), dimension(:), intent(in) :: d
    integer(kind=int32), intent(in) :: ubw

    integer(kind=int32) :: j,k,n
    n=size(a,1)
    a=matmul(u,v)
    do j=1,n-ubw-1
       do k=j+ubw+1,n
          a(j,k)=0.0_dp
       end do
    end do
    do j=1,n
       a(j,j)=d(j)
    end do
  end subroutine d_assemble_lower

  subroutine c_assemble_lower(a,u,v,d,ubw)
    complex(kind=dp), dimension(:,:), intent(out) :: a
    complex(kind=dp), dimension(:,:), intent(in) :: u, v
    complex(kind=dp), dimension(:), intent(in) :: d
    integer(kind=int32), intent(in) :: ubw

    integer(kind=int32) :: j,k,n
    n=size(a,1)
    a=matmul(u,v)
    do j=1,n-ubw-1
       do k=j+ubw+1,n
          a(j,k)=(0.0_dp,0.0_dp)
       end do
    end do
    do j=1,n
       a(j,j)=d(j)
    end do
  end subroutine c_assemble_lower

  subroutine d_assemble_general(a,ul,vl,uu,vu,d)
    real(kind=dp), dimension(:,:), intent(out) :: a
    real(kind=dp), dimension(:,:), intent(in) :: ul, vl, uu, vu
    real(kind=dp), dimension(:), intent(in) :: d

    integer(kind=int32) :: j,k,n
    n=size(a,1)
    a=matmul(ul,vl)
    do k=2,n
       do j=1,k-1
          a(j,k)=dot_product(uu(j,:),vu(:,k))
       end do
    end do
    do j=1,n
       a(j,j)=d(j)
    end do
  end subroutine d_assemble_general

  subroutine c_assemble_general(a,ul,vl,uu,vu,d)
    complex(kind=dp), dimension(:,:), intent(out) :: a
    complex(kind=dp), dimension(:,:), intent(in) :: ul, vl, uu, vu
    complex(kind=dp), dimension(:), intent(in) :: d

    integer(kind=int32) :: j,k,n
    n=size(a,1)
    a=matmul(ul,vl)
    do k=2,n
       do j=1,k-1
          a(j,k)=dot_product(uu(j,:),vu(:,k))
       end do
    end do
    do j=1,n
       a(j,j)=d(j)
    end do
  end subroutine c_assemble_general

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

  subroutine c_output_result_upper(name,a0,a1,ubw0,ubw1,t0,t1,bnd,error)
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
  end subroutine c_output_result_upper

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

  subroutine c_output_result_lower(name,a0,a1,lbw0,lbw1,t0,t1,bnd,error)
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
  end subroutine c_output_result_lower

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

  subroutine c_output_result_lower_upper(name,a0,a1,lbw0,lbw1,ubw0,ubw1,t0,t1,bnd,error)
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
  end subroutine c_output_result_lower_upper



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
       write (*,fmt_qr) name, t1-t0, ubw0, berr, qerr, subd_err, test_result
    end if
  end subroutine c_output_result_qr

end module mod_test_data
