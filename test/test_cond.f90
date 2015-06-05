program test_cond
  use mod_orrb
  implicit none
  integer(kind=int32) :: n, j, k
  type(error_info) :: error
  real(kind=dp), parameter :: tol=1.0e-15_dp, tol1=10*tol
  character(len=*), parameter :: fmt_sv="(A15, '; SV: ',ES8.2, ', Res1: ', ES8.2, ', Res2: ', ES8.2, ', Res3: ', ES8.2, ', ', A10)"
  character(len=40) :: test_name
  
  real(kind=dp), dimension(:,:), allocatable :: l_d, r_d
  real(kind=dp), dimension(:), allocatable :: u_d, v_d, res_d, res_d1
  real(kind=dp) :: sigma_d, tmp

  complex(kind=dp), dimension(:,:), allocatable :: l_z, r_z
  complex(kind=dp), dimension(:), allocatable :: u_z, v_z, res_z, res_z1
  real(kind=dp) :: sigma_z
  
  print *
  print *, "--------------------------------"
  print *
  print *, "Real Condition Estimation Tests"
  print *

  n=30
  allocate(l_d(n,n),u_d(n),v_d(n),res_d(n),res_d1(n))
  call random_matrix_to(l_d)
  do j=1,n
     do k=j+1,n
        l_d(j,k)=0.0_dp
     end do
  end do
  sigma_d=lower_min_sv(l_d,u_d,v_d,res_d,tol,1000,error)
  tmp=norm2(res_d)
  call lower_left_multiply(l_d,v_d,res_d)
  res_d=res_d-sigma_d*u_d
  call lower_tr_left_multiply(l_d,u_d,res_d1)
  res_d1=res_d1-sigma_d*v_d
  test_name='Lower Min SV'
  call output_result(test_name,sigma_d,tmp,norm2(res_d),norm2(res_d1),tol1,error)
  
  sigma_d=lower_max_sv(l_d,u_d,v_d,res_d,tol,1000,error)
  tmp=norm2(res_d)
  call lower_left_multiply(l_d,v_d,res_d)
  res_d=res_d-sigma_d*u_d
  call lower_tr_left_multiply(l_d,u_d,res_d1)
  res_d1=res_d1-sigma_d*v_d
  test_name='Lower Max SV'
  call output_result(test_name,sigma_d,tmp,norm2(res_d),norm2(res_d1),tol1,error)
  
  l_d(5,5)=2.3e-40_dp
  l_d(10,10)=2.3e-40_dp
  sigma_d=lower_min_sv(l_d,u_d,v_d,res_d,tol,1000,error)
  tmp=norm2(res_d)
  call lower_left_multiply(l_d,v_d,res_d)
  res_d=res_d-sigma_d*u_d
  call lower_tr_left_multiply(l_d,u_d,res_d1)
  res_d1=res_d1-sigma_d*v_d
  test_name='Lower Min SV'
  call output_result(test_name,sigma_d,tmp,norm2(res_d),norm2(res_d1),tol1,error)
  deallocate(l_d,u_d,v_d,res_d,res_d1)

  !upper
  n=20
  allocate(r_d(n,n),u_d(n),v_d(n),res_d(n),res_d1(n))
  call random_matrix_to(r_d)
  do j=1,n
     do k=1,j-1
        r_d(j,k)=0.0_dp
     end do
  end do
  sigma_d=upper_min_sv(r_d,u_d,v_d,res_d,tol,1000,error)
  tmp=norm2(res_d)
  call upper_left_multiply(r_d,v_d,res_d)
  res_d=res_d-sigma_d*u_d
  call upper_tr_left_multiply(r_d,u_d,res_d1)
  res_d1=res_d1-sigma_d*v_d
  test_name='Upper Min SV'
  call output_result(test_name,sigma_d,tmp,norm2(res_d),norm2(res_d1),tol1,error)

  sigma_d=upper_max_sv(r_d,u_d,v_d,res_d,tol,1000,error)
  tmp=norm2(res_d)
  call upper_left_multiply(r_d,v_d,res_d)
  res_d=res_d-sigma_d*u_d
  call upper_tr_left_multiply(r_d,u_d,res_d1)
  res_d1=res_d1-sigma_d*v_d

  test_name='Upper Max SV'
  call output_result(test_name,sigma_d,tmp,norm2(res_d),norm2(res_d1),tol1,error)
  
  r_d(5,5)=2.3e-40_dp
  r_d(10,10)=2.3e-40_dp
  sigma_d=upper_min_sv(r_d,u_d,v_d,res_d,tol,1000,error)
  tmp=norm2(res_d)
  call upper_left_multiply(r_d,v_d,res_d)
  res_d=res_d-sigma_d*u_d
  call upper_tr_left_multiply(r_d,u_d,res_d1)
  res_d1=res_d1-sigma_d*v_d

  test_name='Upper Min SV'
  call output_result(test_name,sigma_d,tmp,norm2(res_d),norm2(res_d1),tol1,error)
  deallocate(r_d,u_d,v_d,res_d,res_d1)
  
  print *
  print *, "--------------------------------"
  print *
  print *, "Complex Condition Estimation Tests"
  print *
  n=30
  allocate(l_z(n,n),u_z(n),v_z(n),res_z(n),res_z1(n))
  call random_matrix_to(l_z)
  do j=1,n
     do k=j+1,n
        l_z(j,k)=(0.0_dp,0.0_dp)
     end do
  end do

  sigma_z=lower_min_sv(l_z,u_z,v_z,res_z,1.0e-14_dp,1000,error)

  tmp=norm2(res_z)
  call lower_left_multiply(l_z,v_z,res_z)
  res_z=res_z-sigma_z*u_z
  call lower_tr_left_multiply(l_z,u_z,res_z1)
  res_z1=res_z1-sigma_z*v_z
  test_name='Lower Min SV'
  call output_result(test_name,sigma_z,tmp,norm2(res_z),norm2(res_z1),tol1,error)

  sigma_z=lower_max_sv(l_z,u_z,v_z,res_z,tol,1000,error)
  tmp=norm2(res_z)
  call lower_left_multiply(l_z,v_z,res_z)
  res_z=res_z-sigma_z*u_z
  call lower_tr_left_multiply(l_z,u_z,res_z1)
  res_z1=res_z1-sigma_z*v_z
  test_name='Lower Max SV'
  call output_result(test_name,sigma_z,tmp,norm2(res_z),norm2(res_z1),tol1,error)

  l_z(5,5)=(2.3e-40_dp,5e-40_dp)
  l_z(10,10)=(2.3e-40_dp,5e-40_dp)
  sigma_z=lower_min_sv(l_z,u_z,v_z,res_z,1.0e-14_dp,1000,error)
  tmp=norm2(res_z)
  call lower_left_multiply(l_z,v_z,res_z)
  res_z=res_z-sigma_z*u_z
  call lower_tr_left_multiply(l_z,u_z,res_z1)
  res_z1=res_z1-sigma_z*v_z
  test_name='Lower Min SV'
  call output_result(test_name,sigma_z,tmp,norm2(res_z),norm2(res_z1),tol1,error)
  deallocate(l_z,u_z,v_z,res_z,res_z1)
  
  !upper
  n=20
  allocate(r_z(n,n),u_z(n),v_z(n),res_z(n),res_z1(n))
  call random_matrix_to(r_z)
  do j=1,n
     do k=1,j-1
        r_z(j,k)=0.0_dp
     end do
  end do
  sigma_z=upper_min_sv(r_z,u_z,v_z,res_z,tol,1000,error)
  tmp=norm2(res_z)
  call upper_left_multiply(r_z,v_z,res_z)
  res_z=res_z-sigma_z*u_z
  call upper_tr_left_multiply(r_z,u_z,res_z1)
  res_z1=res_z1-sigma_z*v_z
  test_name='Upper Min SV'
  call output_result(test_name,sigma_z,tmp,norm2(res_z),norm2(res_z1),tol1,error)

  sigma_z=upper_max_sv(r_z,u_z,v_z,res_z,tol,1000,error)
  tmp=norm2(res_z)
  call upper_left_multiply(r_z,v_z,res_z)
  res_z=res_z-sigma_z*u_z
  call upper_tr_left_multiply(r_z,u_z,res_z1)
  res_z1=res_z1-sigma_z*v_z
  test_name='Upper Max SV'
  call output_result(test_name,sigma_z,tmp,norm2(res_z),norm2(res_z1),tol1,error)

  r_z(5,5)=(2.3e-40_dp,5e-40_dp)
  r_z(10,10)=(2.3e-40_dp,5e-40_dp)
  sigma_z=upper_min_sv(r_z,u_z,v_z,res_z,tol,1000,error)
  tmp=norm2(res_z)
  call upper_left_multiply(r_z,v_z,res_z)
  res_z=res_z-sigma_z*u_z
  call upper_tr_left_multiply(r_z,u_z,res_z1)
  res_z1=res_z1-sigma_z*v_z
  test_name='Upper Min SV'
  call output_result(test_name,sigma_z,tmp,norm2(res_z),norm2(res_z1),tol1,error)
  deallocate(r_z,u_z,v_z,res_z,res_z1)
  print *

contains
  
  subroutine output_result(name,sv,res1,res2,res3,bnd,error)
    character(len=*) :: name
    real(kind=dp) :: bnd, sv, res1, res2, res3
    type(error_info) :: error
    character(len=10) :: test_result

    if (error%code > 0) then
       print *, "Calling error in test: ", name
    else
       if (res1 < bnd .and. res2 < bnd .and. res3 < bnd) then
          test_result="PASSED"
       else
          test_result="    FAILED"
       end if
       write (*,fmt_sv) name, sv, res1, res2, res3, test_result
    end if
  end subroutine output_result

end program test_cond
