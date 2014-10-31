program test_qr_iteration
  use mod_nested
  use mod_test_data
  implicit none
  real(kind=dp) :: t0, t1
  integer(kind=int32) :: j,k, na
  type(error_info) :: error
  integer(kind=int32) :: lbwmaxa=3, ubwmaxa=30
  real(kind=dp), parameter :: tol=1e-14, tol1=1e-14, tol2=1e-12
  !
  complex(kind=dp), dimension(:,:), allocatable :: a_c, a0_c, a1_c, q
  complex(kind=dp), dimension(:), allocatable :: u_c, v_c

  type(c_ub) :: ub_c
  type(c_bv) :: bv_c

  type(c_shift), dimension(4) :: shifts
  type(c_sweeps1) :: sw_c

  shifts%flag=.false.
  shifts(1)%flag=.true.
  shifts%shift=(0.0_dp, 0.0_dp)

  
  !! Single Shift Complex QR Tests
  print *
  print *, "--------------------------------"
  print *
  print *, "Single Shift Complex QR Tests"
  print *
  na=3
  allocate(a_c(na,na), q(na,na), a0_c(na,na), a1_c(na,na),u_c(na), v_c(na))
  ub_c=c_new_ub(na,lbwmaxa,ubwmaxa)
  bv_c=c_new_bv(na,lbwmaxa,ubwmaxa)
  a_c=transpose(reshape([ complex(kind=dp) :: &
       0, 0, -24, &
       1, 0 ,50, &
       0, 1, -35 ], shape(a_c)))
  u_c=a_c(:,na); u_c(1)=u_c(1)-1
  v_c=(0.0_dp,0.0_dp); v_c(na)=1
  q=reshape([ complex(kind=dp) :: ((c_delta(j,k), j=1,na), k=1,na) ], shape(q))
  call upper_to_bv(a_c,bv_c,1,tol2,error)
  call bv_to_upper(bv_c, a_c, error)
  a0_c=a_c
  call cpu_time(t0)
  call ss_r1_qr(bv_c,u_c,v_c,q,error)
  call cpu_time(t1)
  call bv_to_upper(bv_c,a_c,error)
  test_name = "Companion SSQR, n=3;"
  a1_c=matmul(matmul(transpose(conjg(q)),a0_c),q)
  call c_output_result_qr(test_name, a0_c, a_c,q,bv_c%ubw,bv_c%ubw,t0,t1,tol2,error)
  deallocate(a_c, q, a0_c, a1_c, u_c, v_c)
  call deallocate_ub(ub_c); call deallocate_bv(bv_c)
  
  na=4
  allocate(a_c(na,na), q(na,na), a0_c(na,na), a1_c(na,na), u_c(na), v_c(na))
  ub_c=c_new_ub(na,lbwmaxa,ubwmaxa)
  bv_c=c_new_bv(na,lbwmaxa,ubwmaxa)
  sw_c=c_new_sweeps1(4,4)
  a_c=transpose(reshape([ complex(kind=dp) :: &
       0, 0, 0, -24, &
       1, 0 ,0, 50, &
       0, 1, 0, -35, &
       0, 0, 1, 10 ], shape(a_c)))
  u_c=a_c(:,na); u_c(1)=u_c(1)-1
  v_c=(0.0_dp,0.0_dp); v_c(na)=1
  q=reshape([ complex(kind=dp) :: ((c_delta(j,k), j=1,na), k=1,na) ], shape(q))
  call upper_to_bv(a_c,bv_c,1,tol2,error)
  call bv_to_upper(bv_c, a_c, error)
  a0_c=a_c
  call cpu_time(t0)
  call ss_r1_qr(bv_c,u_c,v_c,q,error)
  call cpu_time(t1)
  call bv_to_upper(bv_c,a_c,error)
  test_name = "Companion SSQR, n=4;"
  a1_c=matmul(matmul(transpose(conjg(q)),a0_c),q)
  call c_output_result_qr(test_name, a0_c, a_c,q,bv_c%ubw,bv_c%ubw,t0,t1,tol2,error)
  deallocate(a_c, q, a0_c, a1_c, u_c, v_c)
  call deallocate_ub(ub_c); call deallocate_bv(bv_c)
  
  na=6
  allocate(a_c(na,na), q(na,na), a0_c(na,na), a1_c(na,na), u_c(na), v_c(na))
  ub_c=c_new_ub(na,lbwmaxa,ubwmaxa)
  bv_c=c_new_bv(na,lbwmaxa,ubwmaxa)
  a_c=transpose(reshape([ complex(kind=dp) :: &
       0, 0, 0, 0, 0, -24, &
       1, 0 ,0, 0, 0, 50, &
       0, 1, 0, 0, 0, -35, &
       0, 0, 1, 0, 0, 10, &
       0, 0 , 0, 1, 0, 5, &
       0 , 0 , 0, 0, 1, 1], shape(a_c)))
  u_c=a_c(:,na); u_c(1)=u_c(1)-1
  v_c=(0.0_dp,0.0_dp); v_c(na)=1
  q=reshape([ complex(kind=dp) :: ((c_delta(j,k), j=1,na), k=1,na) ], shape(q))
  call upper_to_bv(a_c,bv_c,1,tol2,error)
  call bv_to_upper(bv_c, a_c, error)
  a0_c=a_c
  call cpu_time(t0)
  call ss_r1_qr(bv_c,u_c,v_c,q,error)
  call cpu_time(t1)
  call bv_to_upper(bv_c,a_c,error)
  test_name = "Companion SSQR, n=6;"
  a1_c=matmul(matmul(transpose(conjg(q)),a0_c),q)
  call c_output_result_qr(test_name, a0_c, a_c,q,bv_c%ubw,bv_c%ubw,t0,t1,10*na*tol2,error)
  deallocate(a_c, q, a0_c, a1_c, u_c, v_c)
  call deallocate_ub(ub_c); call deallocate_bv(bv_c)
  
  ! na=200
  ! allocate(a_c(na,na), q(na,na), a0_c(na,na), a1_c(na,na), u_c(na), v_c(na))
  ! ub_c=c_new_ub(na,lbwmaxa,ubwmaxa)
  ! bv_c=c_new_bv(na,lbwmaxa,ubwmaxa)
  ! a_c=reshape([ complex(kind=dp) :: ((c_delta(j-1,k), j=1,na), k=1,na) ], shape(a_c))
  ! call random_matrix(a_c(:,na))
  ! u_c=a_c(:,na); u_c(1)=u_c(1)-1
  ! v_c=(0.0_dp,0.0_dp); v_c(na)=1
  ! q=reshape([ complex(kind=dp) :: ((c_delta(j,k), j=1,na), k=1,na) ], shape(q))
  ! call upper_to_bv(a_c,bv_c,1,tol2,error)
  ! call bv_to_upper(bv_c, a_c, error)
  ! a0_c=a_c
  ! call cpu_time(t0)
  ! call ss_r1_qr(bv_c,u_c,v_c,q,error)
  ! call cpu_time(t1)
  ! call bv_to_upper(bv_c,a_c,error)
  ! test_name = "Companion SSQR, n=200;"
  ! a1_c=matmul(matmul(transpose(conjg(q)),a0_c),q)
  ! call c_output_result_qr(test_name, a0_c, a_c,q,bv_c%ubw,bv_c%ubw,t0,t1,na*tol2,error)
  ! deallocate(a_c, q, a0_c, a1_c, u_c, v_c)
  ! call deallocate_ub(ub_c); call deallocate_bv(bv_c)

  !! Modified shift
  na=50
  allocate(a_c(na,na), q(na,na), a0_c(na,na), a1_c(na,na), u_c(na), v_c(na))
  ub_c=c_new_ub(na,lbwmaxa,ubwmaxa)
  bv_c=c_new_bv(na,lbwmaxa,ubwmaxa)
  a_c=reshape([ complex(kind=dp) :: ((c_delta(j-1,k), j=1,na), k=1,na) ], shape(a_c))
  !  call random_matrix(a_c(:,na))
  a_c(1,na)=1
  a_c(26,25)=(0.0_dp,0.0_dp)
  u_c=(0.0_dp,0.0_dp); u_c(26)=(-1.0_dp,0.0_dp)
  v_c=(0.0_dp,0.0_dp); v_c(25)=1
  q=reshape([ complex(kind=dp) :: ((c_delta(j,k), j=1,na), k=1,na) ], shape(q))
  call upper_to_bv(a_c,bv_c,1,tol2,error)
  call bv_to_upper(bv_c, a_c, error)
  a0_c=a_c
  call cpu_time(t0)
  call ss_r1_qr(bv_c,u_c,v_c,q,error)
  call cpu_time(t1)
  call bv_to_upper(bv_c,a_c,error)
  test_name = "Reduced Shift SSQR, n=50;"
  a1_c=matmul(matmul(transpose(conjg(q)),a0_c),q)
  call c_output_result_qr(test_name, a0_c, a_c,q,bv_c%ubw,bv_c%ubw,t0,t1,na*tol2,error)
  deallocate(a_c, q, a0_c, a1_c, u_c, v_c)
  call deallocate_ub(ub_c); call deallocate_bv(bv_c)

  !! Reduced matrix
  ! na=50
  ! allocate(a_c(na,na), q(na,na), a0_c(na,na), a1_c(na,na), u_c(na), v_c(na))
  ! ub_c=c_new_ub(na,lbwmaxa,ubwmaxa)
  ! bv_c=c_new_bv(na,lbwmaxa,ubwmaxa)
  ! ! rank one modification of the identity.
  ! a_c=reshape([ complex(kind=dp) :: ((c_delta(j,k), j=1,na), k=1,na) ], shape(a_c))
  ! call random_matrix(u_c)
  ! a_c(:,na)=a_c(:,na)+u_c
  ! sw_c=c_random_sweeps1(na/2,1)
  ! call sweeps1_times_general(sw_c, a_c(1:na/2,:))
  ! call sweeps1_times_general(sw_c, u_c(1:na/2))
  ! call deallocate_sweeps1(sw_c)
  ! sw_c=c_random_sweeps1(na/2,1)
  ! call sweeps1_times_general(sw_c, a_c(na/2+1:na,:))
  ! call sweeps1_times_general(sw_c, u_c(na/2+1:na))
  ! v_c=(0.0_dp,0.0_dp); v_c(na)=1
  ! q=reshape([ complex(kind=dp) :: ((c_delta(j,k), j=1,na), k=1,na) ], shape(q))
  ! call upper_to_bv(a_c,bv_c,1,tol2,error)
  ! call bv_to_upper(bv_c, a_c, error)
  ! a0_c=a_c
  ! call cpu_time(t0)
  ! call ss_r1_qr(bv_c,u_c,v_c,q,error)
  ! call cpu_time(t1)
  ! call bv_to_upper(bv_c,a_c,error)
  ! test_name = "Reduced Hessenberg SSQR, n=50;"
  ! a1_c=matmul(matmul(transpose(conjg(q)),a0_c),q)
  ! call c_output_result_qr(test_name, a0_c, a_c,q,bv_c%ubw,bv_c%ubw,t0,t1,na*tol2,error)
  ! deallocate(a_c, q, a0_c, a1_c, u_c, v_c)
  ! call deallocate_ub(ub_c); call deallocate_bv(bv_c)


end program test_qr_iteration
