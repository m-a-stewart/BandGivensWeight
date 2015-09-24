program test_cond_orth_band
  use mod_orrb
  implicit none
  integer(kind=int32) :: na, lbw, ubw
  type(error_info) :: error
  real(kind=dp), parameter :: tol=1.0e-15_dp, tol1=10*tol
  character(len=*), parameter :: fmt_sv="(A15, '; SV: ',ES8.2, ', Res1: ', ES8.2, ', Res2: ', ES8.2, ', Res3: ', ES8.2, ', ', A10)"
  character(len=40) :: test_name
  
  real(kind=dp), dimension(:), allocatable :: u_d, v_d, res_d, res_d1
  real(kind=dp) :: sigma_d, tmp
  real(kind=dp), allocatable, dimension(:,:) :: a_d
  complex(kind=dp), allocatable, dimension(:,:) :: a_z

  complex(kind=dp), dimension(:), allocatable :: u_z, v_z, res_z, res_z1
  real(kind=dp) :: sigma_z

  type(d_ub), allocatable :: ub_d
  type(z_ub), allocatable :: ub_z

  call initialize_errors
  print *
  print *, "--------------------------------"
  print *
  print *, "Real UB Condition Estimation Tests"
  print *

  na=50; lbw=0; ubw=5
  allocate(u_d(na),v_d(na),res_d(na),res_d1(na))
  ub_d=d_random_ub(na,lbw,ubw,error=error)
  a_d=general(ub_d, error)
  sigma_d=d_ub_min_sv(ub_d,u_d,v_d,res_d,tol,1000,error)
  tmp=norm2(res_d)
  res_d=matmul(a_d,v_d)-sigma_d*u_d
  res_d1=matmul(u_d,a_d)-sigma_d*v_d
  test_name='UB Min SV'
  call output_result(test_name,sigma_d,tmp,norm2(res_d),norm2(res_d1),tol1,error)
  deallocate(u_d,v_d,res_d,res_d1)

  na=10; lbw=3; ubw=5
  allocate(u_d(na),v_d(na),res_d(na),res_d1(na))
  ub_d=d_random_ub(na,lbw,ubw,error=error)
  a_d=general(ub_d, error)
  sigma_d=d_ub_max_sv(ub_d,u_d,v_d,res_d,tol,1000,error)
  tmp=norm2(res_d)
  res_d=matmul(a_d,v_d)-sigma_d*u_d
  res_d1=matmul(u_d,a_d)-sigma_d*v_d
  test_name='UB Max SV'
  call output_result(test_name,sigma_d,tmp,norm2(res_d),norm2(res_d1),tol1,error)
  deallocate(u_d,v_d,res_d,res_d1)
  
  print *
  print *, "--------------------------------"
  print *
  print *, "Complex UB Condition Estimation Tests"
  print *

  na=50; lbw=0; ubw=5
  allocate(u_z(na),v_z(na),res_z(na),res_z1(na))
  ub_z=z_random_ub(na,lbw,ubw,error=error)
  a_z=general(ub_z, error)
  sigma_z=z_ub_min_sv(ub_z,u_z,v_z,res_z,tol,100,error)
  tmp=norm2(res_z)
  res_z=matmul(a_z,v_z)-sigma_z*u_z
  res_z1=matmul(conjg(u_z),a_z)-sigma_z*conjg(v_z)
  test_name='UB Min SV'
  call output_result(test_name,sigma_z,tmp,norm2(res_z),norm2(res_z1),tol1,error)
  deallocate(u_z,v_z,res_z,res_z1)

  na=10; lbw=3; ubw=5
  allocate(u_z(na),v_z(na),res_z(na),res_z1(na))
  ub_z=z_random_ub(na,lbw,ubw,error=error)
  a_z=general(ub_z, error)
  sigma_z=z_ub_max_sv(ub_z,u_z,v_z,res_z,tol,1000,error)
  tmp=norm2(res_z)
  res_z=matmul(a_z,v_z)-sigma_z*u_z
  res_z1=matmul(conjg(u_z),a_z)-sigma_z*conjg(v_z)
  test_name='UB Max SV'
  call output_result(test_name,sigma_z,tmp,norm2(res_z),norm2(res_z1),tol1,error)
  deallocate(u_z,v_z,res_z,res_z1)
  
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

end program test_cond_orth_band
