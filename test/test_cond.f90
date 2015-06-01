program test_cond
  use mod_orb
  implicit none

  integer(kind=int32) :: n, j, k
  type(error_info) :: error
  real(kind=dp), parameter :: tol=1.0e-15_dp, tol1=10*tol
  character(len=*), parameter :: fmt_sv="(A25, '; SV: ',ES8.2, ', Res1: ', ES8.2, ', Res2: ', ES8.2, ', Res3: ', ES8.2, ', ', A10)"

  real(kind=dp), dimension(:,:), allocatable :: l_d, r_d
  real(kind=dp), dimension(:), allocatable :: u_d, v_d, res_d, res_d1
  real(kind=dp) :: sigma_d, tmp

  complex(kind=dp), dimension(:,:), allocatable :: l_c, r_c
  complex(kind=dp), dimension(:), allocatable :: u_c, v_c, res_c, res_c1
  real(kind=dp) :: sigma_c
  
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
  if (tmp < tol1 .and. norm2(res_d) < tol1 .and. norm2(res_d1) < tol1) then
     write (*, fmt_sv) "Lower Min SV", sigma_d, tmp, norm2(res_d), norm2(res_d1), "PASSED"
  else
     write (*, fmt_sv) "Lower Min SV", sigma_d, tmp, norm2(res_d), norm2(res_d1), "FAILED"
  end if

  sigma_d=lower_max_sv(l_d,u_d,v_d,res_d,tol,1000,error)
  tmp=norm2(res_d)
  call lower_left_multiply(l_d,v_d,res_d)
  res_d=res_d-sigma_d*u_d
  call lower_tr_left_multiply(l_d,u_d,res_d1)
  res_d1=res_d1-sigma_d*v_d
  if (tmp < tol1 .and. norm2(res_d) < tol1 .and. norm2(res_d1) < tol1) then
     write (*, fmt_sv) "Lower Max SV", sigma_d, tmp, norm2(res_d), norm2(res_d1), "PASSED"
  else
     write (*, fmt_sv) "Lower Max SV", sigma_d, tmp, norm2(res_d), norm2(res_d1), "FAILED"
  end if
  
  l_d(5,5)=2.3e-40_dp
  l_d(10,10)=2.3e-40_dp
  sigma_d=lower_min_sv(l_d,u_d,v_d,res_d,tol,1000,error)
  tmp=norm2(res_d)
  call lower_left_multiply(l_d,v_d,res_d)
  res_d=res_d-sigma_d*u_d
  call lower_tr_left_multiply(l_d,u_d,res_d1)
  res_d1=res_d1-sigma_d*v_d
  if (tmp < tol1 .and. norm2(res_d) < tol1 .and. norm2(res_d1) < tol1) then
     write (*, fmt_sv) "Lower Min SV", sigma_d, tmp, norm2(res_d), norm2(res_d1), "PASSED"
  else
     write (*, fmt_sv) "Lower Min SV", sigma_d, tmp, norm2(res_d), norm2(res_d1), "FAILED"
  end if
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
  if (tmp < tol1 .and. norm2(res_d) < tol1 .and. norm2(res_d1) < tol1) then
     write (*, fmt_sv) "Upper Min SV", sigma_d, tmp, norm2(res_d), norm2(res_d1), "PASSED"
  else
     write (*, fmt_sv) "Upper Min SV", sigma_d, tmp, norm2(res_d), norm2(res_d1), "FAILED"
  end if

  sigma_d=upper_max_sv(r_d,u_d,v_d,res_d,tol,1000,error)
  tmp=norm2(res_d)
  call upper_left_multiply(r_d,v_d,res_d)
  res_d=res_d-sigma_d*u_d
  call upper_tr_left_multiply(r_d,u_d,res_d1)
  res_d1=res_d1-sigma_d*v_d
  if (tmp < tol1 .and. norm2(res_d) < tol1 .and. norm2(res_d1) < tol1) then
     write (*, fmt_sv) "Upper Max SV", sigma_d, tmp, norm2(res_d), norm2(res_d1), "PASSED"
  else
     write (*, fmt_sv) "Upper Max SV", sigma_d, tmp, norm2(res_d), norm2(res_d1), "FAILED"
  end if
  
  r_d(5,5)=2.3e-40_dp
  r_d(10,10)=2.3e-40_dp
  sigma_d=upper_min_sv(r_d,u_d,v_d,res_d,tol,1000,error)
  tmp=norm2(res_d)
  call upper_left_multiply(r_d,v_d,res_d)
  res_d=res_d-sigma_d*u_d
  call upper_tr_left_multiply(r_d,u_d,res_d1)
  res_d1=res_d1-sigma_d*v_d
  if (tmp < tol1 .and. norm2(res_d) < tol1 .and. norm2(res_d1) < tol1) then
     write (*, fmt_sv) "Upper Min SV", sigma_d, tmp, norm2(res_d), norm2(res_d1), "PASSED"
  else
     write (*, fmt_sv) "Upper Min SV", sigma_d, tmp, norm2(res_d), norm2(res_d1), "FAILED"
  end if
  deallocate(r_d,u_d,v_d,res_d,res_d1)
  
  print *
  print *, "--------------------------------"
  print *
  print *, "Complex Condition Estimation Tests"
  print *
  n=30
  allocate(l_c(n,n),u_c(n),v_c(n),res_c(n),res_c1(n))
  call random_matrix_to(l_c)
  do j=1,n
     do k=j+1,n
        l_c(j,k)=(0.0_dp,0.0_dp)
     end do
  end do

  sigma_c=lower_min_sv(l_c,u_c,v_c,res_c,1.0e-14_dp,1000,error)

  tmp=norm2(res_c)
  call lower_left_multiply(l_c,v_c,res_c)
  res_c=res_c-sigma_c*u_c
  call lower_tr_left_multiply(l_c,u_c,res_c1)
  res_c1=res_c1-sigma_c*v_c
  if (tmp < tol1 .and. norm2(res_c) < tol1 .and. norm2(res_c1) < tol1) then
     write (*, fmt_sv) "Lower Min SV", sigma_c, tmp, norm2(res_c), norm2(res_c1), "PASSED"
  else
     write (*, fmt_sv) "Lower Min SV", sigma_c, tmp, norm2(res_c), norm2(res_c1), "FAILED"
  end if

  sigma_c=lower_max_sv(l_c,u_c,v_c,res_c,tol,1000,error)
  tmp=norm2(res_c)
  call lower_left_multiply(l_c,v_c,res_c)
  res_c=res_c-sigma_c*u_c
  call lower_tr_left_multiply(l_c,u_c,res_c1)
  res_c1=res_c1-sigma_c*v_c
  if (tmp < tol1 .and. norm2(res_c) < tol1 .and. norm2(res_c1) < tol1) then
     write (*, fmt_sv) "Lower Max SV", sigma_c, tmp, norm2(res_c), norm2(res_c1), "PASSED"
  else
     write (*, fmt_sv) "Lower Max SV", sigma_c, tmp, norm2(res_c), norm2(res_c1), "FAILED"
  end if

  l_c(5,5)=(2.3e-40_dp,5e-40_dp)
  l_c(10,10)=(2.3e-40_dp,5e-40_dp)
  sigma_c=lower_min_sv(l_c,u_c,v_c,res_c,1.0e-14_dp,1000,error)

  tmp=norm2(res_c)
  call lower_left_multiply(l_c,v_c,res_c)
  res_c=res_c-sigma_c*u_c
  call lower_tr_left_multiply(l_c,u_c,res_c1)
  res_c1=res_c1-sigma_c*v_c
  if (tmp < tol1 .and. norm2(res_c) < tol1 .and. norm2(res_c1) < tol1) then
     write (*, fmt_sv) "Lower Min SV", sigma_c, tmp, norm2(res_c), norm2(res_c1), "PASSED"
  else
     write (*, fmt_sv) "Lower Min SV", sigma_c, tmp, norm2(res_c), norm2(res_c1), "FAILED"
  end if
  deallocate(l_c,u_c,v_c,res_c,res_c1)
  
  !upper
  n=20
  allocate(r_c(n,n),u_c(n),v_c(n),res_c(n),res_c1(n))
  call random_matrix_to(r_c)
  do j=1,n
     do k=1,j-1
        r_c(j,k)=0.0_dp
     end do
  end do
  sigma_c=upper_min_sv(r_c,u_c,v_c,res_c,tol,1000,error)
  tmp=norm2(res_c)
  call upper_left_multiply(r_c,v_c,res_c)
  res_c=res_c-sigma_c*u_c
  call upper_tr_left_multiply(r_c,u_c,res_c1)
  res_c1=res_c1-sigma_c*v_c
  if (tmp < tol1 .and. norm2(res_c) < tol1 .and. norm2(res_c1) < tol1) then
     write (*, fmt_sv) "Upper Min SV", sigma_c, tmp, norm2(res_c), norm2(res_c1), "PASSED"
  else
     write (*, fmt_sv) "Upper Min SV", sigma_c, tmp, norm2(res_c), norm2(res_c1), "FAILED"
  end if

  sigma_c=upper_max_sv(r_c,u_c,v_c,res_c,tol,1000,error)
  tmp=norm2(res_c)
  call upper_left_multiply(r_c,v_c,res_c)
  res_c=res_c-sigma_c*u_c
  call upper_tr_left_multiply(r_c,u_c,res_c1)
  res_c1=res_c1-sigma_c*v_c
  if (tmp < tol1 .and. norm2(res_c) < tol1 .and. norm2(res_c1) < tol1) then
     write (*, fmt_sv) "Upper Max SV", sigma_c, tmp, norm2(res_c), norm2(res_c1), "PASSED"
  else
     write (*, fmt_sv) "Upper Max SV", sigma_c, tmp, norm2(res_c), norm2(res_c1), "FAILED"
  end if

  r_c(5,5)=(2.3e-40_dp,5e-40_dp)
  r_c(10,10)=(2.3e-40_dp,5e-40_dp)
  sigma_c=upper_min_sv(r_c,u_c,v_c,res_c,tol,1000,error)
  tmp=norm2(res_c)
  call upper_left_multiply(r_c,v_c,res_c)
  res_c=res_c-sigma_c*u_c
  call upper_tr_left_multiply(r_c,u_c,res_c1)
  res_c1=res_c1-sigma_c*v_c
  if (tmp < tol1 .and. norm2(res_c) < tol1 .and. norm2(res_c1) < tol1) then
     write (*, fmt_sv) "Upper Min SV", sigma_c, tmp, norm2(res_c), norm2(res_c1), "PASSED"
  else
     write (*, fmt_sv) "Upper Min SV", sigma_c, tmp, norm2(res_c), norm2(res_c1), "FAILED"
  end if
  deallocate(r_c,u_c,v_c,res_c,res_c1)
  print *

end program test_cond
