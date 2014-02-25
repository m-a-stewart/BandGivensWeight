program test_qr_iteration
  use general_ub
  use general_bv
  use utility
  use assemble
  use conversions_ub_to_bv
  use conversions_bv_to_ub
  use compressions_ub_to_bv
  use compressions_bv_to_ub
  use band_types
  use nested_types
  use qr_iteration
  use test_data
  implicit none
  real(kind=dp) :: t0, t1
  integer(kind=int32) :: j,k, na
  type(error_info) :: error
  integer(kind=int32) :: lbwmaxa=2, ubwmaxa=5
  !
  complex(kind=dp), dimension(:,:), allocatable :: a_c, a0_c, a1_c, q

  type(c_ub) :: ub_c
  type(c_bv) :: bv_c

  
  !! Single Shift Complex QR Tests
  print *
  print *, "--------------------------------"
  print *
  print *, "Single Shift Complex QR Tests"
  print *
  na=3
  allocate(a_c(na,na), q(na,na), a0_c(na,na), a1_c(na,na))
  ub_c=c_new_ub(na,lbwmaxa,ubwmaxa)
  bv_c=c_new_bv(na,lbwmaxa,ubwmaxa)
  a_c=transpose(reshape([ complex(kind=dp) :: &
       0, 0, -24, &
       1, 0 ,50, &
       0, 1, -35 ], shape(a_c)))
  q=reshape([ complex(kind=dp) :: ((c_delta(j,k), j=1,na), k=1,na) ], shape(q))
  a0_c=a_c
  call upper_to_bv(a_c,bv_c,1,tol2,error)
  call cpu_time(t0)
  call ss_qr(bv_c,q,tol2, tol2, error)
  call cpu_time(t1)
  call bv_to_upper(bv_c,a_c,error)
  test_name = "Companion SSQR, n=3;"
  a1_c=matmul(matmul(transpose(conjg(q)),a0_c),q)
  call c_output_result_qr(test_name, a0_c, a_c,q,bv_c%ubw,3,t0,t1,tol2,error)
  deallocate(a_c, q, a0_c, a1_c)
  call deallocate_ub(ub_c); call deallocate_bv(bv_c)
  !
  !
  na=4
  allocate(a_c(na,na), q(na,na), a0_c(na,na), a1_c(na,na))
  ub_c=c_new_ub(na,lbwmaxa,ubwmaxa)
  bv_c=c_new_bv(na,lbwmaxa,ubwmaxa)
  a_c=transpose(reshape([ complex(kind=dp) :: &
       0, 0, 0, -24, &
       1, 0 ,0, 50, &
       0, 1, 0, -35, &
       0, 0, 1, 10 ], shape(a_c)))
  q=reshape([ complex(kind=dp) :: ((c_delta(j,k), j=1,na), k=1,na) ], shape(q))
  a0_c=a_c
  call upper_to_bv(a_c,bv_c,1,tol2,error)
  call cpu_time(t0)
  call ss_qr(bv_c,q,tol2, tol2, error)
  call cpu_time(t1)
  call bv_to_upper(bv_c,a_c,error)
  test_name = "Companion SSQR, n=4;"
  a1_c=matmul(matmul(transpose(conjg(q)),a0_c),q)
  call c_output_result_qr(test_name, a0_c, a_c,q,bv_c%ubw,3,t0,t1,tol2,error)
  deallocate(a_c, q, a0_c, a1_c)
  call deallocate_ub(ub_c); call deallocate_bv(bv_c)
  !
  na=6
  allocate(a_c(na,na), q(na,na), a0_c(na,na), a1_c(na,na))
  ub_c=c_new_ub(na,lbwmaxa,ubwmaxa)
  bv_c=c_new_bv(na,lbwmaxa,ubwmaxa)
  a_c=transpose(reshape([ complex(kind=dp) :: &
       0, 0, 0, 0, 0, -24, &
       1, 0 ,0, 0, 0, 50, &
       0, 1, 0, 0, 0, -35, &
       0, 0, 1, 0, 0, 10, &
       0, 0 , 0, 1, 0, 5, &
       0 , 0 , 0, 0, 1, 1], shape(a_c)))
  q=reshape([ complex(kind=dp) :: ((c_delta(j,k), j=1,na), k=1,na) ], shape(q))
  a0_c=a_c
  call upper_to_bv(a_c,bv_c,1,tol2,error)
  call cpu_time(t0)
  call ss_qr(bv_c,q,tol2, tol2, error)
  call cpu_time(t1)
  call bv_to_upper(bv_c,a_c,error)
  test_name = "Companion SSQR, n=6;"
  a1_c=matmul(matmul(transpose(conjg(q)),a0_c),q)
  call c_output_result_qr(test_name, a0_c, a_c,q,bv_c%ubw,3,t0,t1,tol2,error)
  deallocate(a_c, q, a0_c, a1_c)
  call deallocate_ub(ub_c); call deallocate_bv(bv_c)
  !
  na=50
  allocate(a_c(na,na), q(na,na), a0_c(na,na), a1_c(na,na))
  ub_c=c_new_ub(na,lbwmaxa,ubwmaxa)
  bv_c=c_new_bv(na,lbwmaxa,ubwmaxa)
  a_c=reshape([ complex(kind=dp) :: ((c_delta(j-1,k), j=1,na), k=1,na) ], shape(a_c))
  call random_complex(a_c(:,na))
  q=reshape([ complex(kind=dp) :: ((c_delta(j,k), j=1,na), k=1,na) ], shape(q))
  a0_c=a_c
  call upper_to_bv(a_c,bv_c,1,tol2,error)
  call cpu_time(t0)
  call ss_qr(bv_c,q,tol2, tol2, error)
  call cpu_time(t1)
  call bv_to_upper(bv_c,a_c,error)
  test_name = "Companion SSQR, n=50;"
  a1_c=matmul(matmul(transpose(conjg(q)),a0_c),q)
  call c_output_result_qr(test_name, a0_c, a_c,q,bv_c%ubw,3,t0,t1,na*tol2,error)
  deallocate(a_c, q, a0_c, a1_c)
  call deallocate_ub(ub_c); call deallocate_bv(bv_c)
  !
  na=50
  allocate(a_c(na,na), q(na,na), a0_c(na,na), a1_c(na,na))
  ub_c=c_new_ub(na,lbwmaxa,ubwmaxa)
  bv_c=c_new_bv(na,lbwmaxa,ubwmaxa)
  a_c=reshape([ complex(kind=dp) :: ((c_delta(j-1,k), j=1,na), k=1,na) ], shape(a_c))
  call random_complex(a_c(:,na))
  a_c(26,25)=(0.0_dp,0.0_dp)
  q=reshape([ complex(kind=dp) :: ((c_delta(j,k), j=1,na), k=1,na) ], shape(q))
  a0_c=a_c
  call upper_to_bv(a_c,bv_c,1,tol2,error)
  call cpu_time(t0)
  call ss_qr(bv_c,q,tol2, tol2, error)
  call cpu_time(t1)
  call bv_to_upper(bv_c,a_c,error)
  test_name = "Reduced Comp. SSQR, n=50;"
  a1_c=matmul(matmul(transpose(conjg(q)),a0_c),q)
  call c_output_result_qr(test_name, a0_c, a_c,q,bv_c%ubw,3,t0,t1,na*tol2,error)
  deallocate(a_c, q, a0_c, a1_c)
  call deallocate_ub(ub_c); call deallocate_bv(bv_c)

end program test_qr_iteration
