program test_recompress
  use nested
  use transforms
  use test_data
  implicit none
  real(kind=dp) :: t0, t1
  integer(kind=int32) :: na, lbwa, numsweeps
  type(error_info) :: error
  integer, parameter :: n=50, rmax=13, ubwmax=rmax+5, lbw=2, lbwmax=10
  real(kind=dp), parameter :: tol=1e-14, tol1=1e-14, tol2=1e-10
  !
  real(kind=dp), dimension(n,n) :: a, a0, a1
  real(kind=dp), dimension(n,rmax) :: u, u0
  real(kind=dp), dimension(rmax,n) :: v, v0
  real(kind=dp), dimension(n) :: d, d0
  complex(kind=dp), dimension(n,n) :: a_c, a0_c, a1_c
  complex(kind=dp), dimension(n,rmax) :: u_c, u0_c
  complex(kind=dp), dimension(rmax,n) :: v_c, v0_c
  complex(kind=dp), dimension(n) :: d_c, d0_c

  type(d_ub) :: ub_d, ub_na_d
  type(c_ub) :: ub_c, ub_na_c
  type(d_bv) :: bv_d, bv_na_d
  type(c_bv) :: bv_c, bv_na_c

  type(d_sweeps) :: sw_d
  type(c_sweeps) :: sw_c

  ub_d=d_new_ub(n,lbwmax,ubwmax)
  ub_c=c_new_ub(n,lbwmax,ubwmax)
  bv_d=d_new_bv(n,lbwmax,ubwmax)
  bv_c=c_new_bv(n,lbwmax,ubwmax)

  ! real
  call random_seed
  call random_matrix(u0)
  call random_matrix(v0)
  call random_matrix(d0)


  !! UB to BV tests
  print *
  print *, "--------------------------------"
  print *
  print *, "Real UB to BV Recompression Tests:"
  print *
  ! Test 3; case 1 termination test.
  u=u0; v=v0; d=d0
  v(rmax-2:rmax,:)=0.1*tol2*v(rmax-2:rmax,:)
  v(rmax-3,rmax:n)=0.1_dp*tol2*v(rmax-3,rmax:n)
  call d_assemble_upper(a,u,v,d,lbw)
  a0=a
  call upper_to_ub(a,ub_d,lbw,tol1,error)
  call cpu_time(t0)
  call recompress_ub_to_bv(ub_d, bv_d,tol2,tol2,0,error)
  call cpu_time(t1)
  call bv_to_upper(bv_d,a1,error)
  test_name = "Real UB to BV (term. 1);"
  call d_output_result_upper(test_name,a0,a1,rmax-3,bv_d%ubw,t0,t1,tol2,error)
  ! Test 4
  u=u0; v=v0; d=d0
  u(:,rmax-3:rmax)=0.1*tol2*u(:,rmax-3:rmax)
  call d_assemble_upper(a,u,v,d,lbw)
  a0=a
  call upper_to_ub(a,ub_d,lbw,tol1,error)
  call cpu_time(t0)
  call recompress_ub_to_bv(ub_d, bv_d,tol2,tol2,0,error)
  call cpu_time(t1)
  call bv_to_upper(bv_d,a1,error)
  test_name = "Real UB to BV (term. 2);"
  call d_output_result_upper(test_name,a0,a1,rmax-4,bv_d%ubw,t0,t1,tol2,error)
  ! Test 5
  ! case 3 termination test
  u=u0; v=v0; d=d0
  v(rmax-2:rmax,:)=0.1*tol2*v(rmax-2:rmax,:)
  v(rmax-3,rmax+1:n)=0.0_dp
  call d_assemble_upper(a,u,v,d,lbw)
  a0=a
  call upper_to_ub(a,ub_d,lbw, tol1,error)
  call cpu_time(t0)
  call recompress_ub_to_bv(ub_d, bv_d,tol2,tol2,0,error)
  call cpu_time(t1)
  call bv_to_upper(bv_d,a1,error)
  test_name = "Real UB to BV (term. 3);"
  call d_output_result_upper(test_name,a0,a1,rmax-3,bv_d%ubw,t0,t1,tol2,error)
  ! test 6: recompress to one test.
  u=u0; v=v0; d=d0
  u(:,2:rmax)=0.1*tol2*u(:,2:rmax)
  call d_assemble_upper(a,u,v,d,lbw)
  a0=a
  call upper_to_ub(a,ub_d,lbw,tol1,error)
  call cpu_time(t0)
  call recompress_ub_to_bv(ub_d, bv_d,tol2,tol2,0,error)
  call cpu_time(t1)
  call bv_to_upper(bv_d,a1,error)
  test_name = "Real UB to BV (recompress to 1);"
  call d_output_result_upper(test_name,a0,a1,1,bv_d%ubw,t0,t1,tol2,error)
  ! test 7: recompress to zero test.
  u=u0; v=v0; d=d0
  u=0.1*tol2*u
  call d_assemble_upper(a,u,v,d,lbw)
  a0=a
  call upper_to_ub(a,ub_d,lbw,tol1,error)
  call cpu_time(t0)
  call recompress_ub_to_bv(ub_d, bv_d,tol2,tol2,0,error)
  call cpu_time(t1)
  call bv_to_upper(bv_d,a1,error)
  test_name = "Real UB to BV (recompress to 0);"
  call d_output_result_upper(test_name,a0,a1,0,bv_d%ubw,t0,t1,tol2,error)
  ! test 8: no recompression
  u=u0; v=v0; d=d0
  call d_assemble_upper(a,u,v,d,lbw)
  a0=a
  call upper_to_ub(a,ub_d,lbw,tol1,error)
  call cpu_time(t0)
  call recompress_ub_to_bv(ub_d, bv_d,tol2,tol2,0,error)
  call cpu_time(t1)
  call bv_to_upper(bv_d,a1,error)
  test_name = "Real UB to BV (no recompression);"
  call d_output_result_upper(test_name,a0,a1,rmax,bv_d%ubw,t0,t1,tol2,error)
  ! small matrices
  na=1
  lbwa=min(na-1,lbw)
  ub_na_d=d_new_ub(na,lbwmax,ubwmax)
  bv_na_d=d_new_bv(na,lbwmax,ubwmax)
  u=u0; v=v0; d=d0
  call d_assemble_upper(a(1:na,1:na),u(1:na,:),v(:,1:na),d(1:na),lbwa)
  a0(1:na,1:na)=a(1:na,1:na)
  call upper_to_ub(a(1:na,1:na),ub_na_d,lbwa, tol1,error)
  call cpu_time(t0)
  call recompress_ub_to_bv(ub_na_d, bv_na_d,tol2,tol2,0,error)
  call cpu_time(t1)
  call bv_to_upper(bv_na_d,a1(1:na,1:na),error)
  test_name = "Real UB to BV (n=1);"
  call d_output_result_upper(test_name,a0(1:na,1:na),a1(1:na,1:na),0,bv_na_d%ubw,t0,t1,tol2,error)
  call d_deallocate_ub(ub_na_d)
  call d_deallocate_bv(bv_na_d)
  na=2
  lbwa=min(na-1,lbw)
  ub_na_d=d_new_ub(na,lbwmax,ubwmax)
  bv_na_d=d_new_bv(na,lbwmax,ubwmax)
  u=u0; v=v0; d=d0
  call d_assemble_upper(a(1:na,1:na),u(1:na,:),v(:,1:na),d(1:na),lbwa)
  a0(1:na,1:na)=a(1:na,1:na)
  call upper_to_ub(a(1:na,1:na),ub_na_d,lbwa, tol1,error)
  call cpu_time(t0)
  call recompress_ub_to_bv(ub_na_d, bv_na_d,tol2,tol2,0,error)
  call cpu_time(t1)
  call bv_to_upper(bv_na_d,a1(1:na,1:na),error)
  test_name = "Real UB to BV (n=2);"
  call d_output_result_upper(test_name,a0(1:na,1:na),a1(1:na,1:na),1,bv_na_d%ubw,t0,t1,tol2,error)
  call d_deallocate_ub(ub_na_d)
  call d_deallocate_bv(bv_na_d)
  !
  na=3
  lbwa=min(na-1,lbw)
  ub_na_d=d_new_ub(na,lbwmax,ubwmax)
  bv_na_d=d_new_bv(na,lbwmax,ubwmax)
  u=u0; v=v0; d=d0
  call d_assemble_upper(a(1:na,1:na),u(1:na,:),v(:,1:na),d(1:na),lbwa)
  a0(1:na,1:na)=a(1:na,1:na)
  call upper_to_ub(a(1:na,1:na),ub_na_d,lbwa, tol1,error)
  call cpu_time(t0)
  call recompress_ub_to_bv(ub_na_d, bv_na_d,tol2,tol2,0,error)
  call cpu_time(t1)
  call bv_to_upper(bv_na_d,a1(1:na,1:na),error)
  test_name = "Real UB to BV (n=3);"
  call d_output_result_upper(test_name,a0(1:na,1:na),a1(1:na,1:na),1,bv_na_d%ubw,t0,t1,tol2,error)
  call d_deallocate_ub(ub_na_d)
  call d_deallocate_bv(bv_na_d)
  !
  na=4
  lbwa=min(na-1,lbw)
  ub_na_d=d_new_ub(na,lbwmax,ubwmax)
  bv_na_d=d_new_bv(na,lbwmax,ubwmax)
  u=u0; v=v0; d=d0
  call d_assemble_upper(a(1:na,1:na),u(1:na,:),v(:,1:na),d(1:na),lbwa)
  a0(1:na,1:na)=a(1:na,1:na)
  call upper_to_ub(a(1:na,1:na),ub_na_d,lbwa, tol1,error)
  call cpu_time(t0)
  call recompress_ub_to_bv(ub_na_d, bv_na_d,tol2,tol2,0,error)
  call cpu_time(t1)
  call bv_to_upper(bv_na_d,a1(1:na,1:na),error)
  test_name = "Real UB to BV (n=4);"
  call d_output_result_upper(test_name,a0(1:na,1:na),a1(1:na,1:na),2,bv_na_d%ubw,t0,t1,tol2,error)
  call d_deallocate_ub(ub_na_d)
  call d_deallocate_bv(bv_na_d)
  ! maximum upper and lower bandwidth
  na=10; lbwa=6; numsweeps=7
  sw_d=d_random_sweeps(na,numsweeps)
  ub_na_d=d_new_ub(na,lbwmax,ubwmax)
  bv_na_d=d_new_bv(na,lbwmax,ubwmax)
  u=u0; v=v0; d=d0
  call d_assemble_upper(a(1:na,1:na),u(1:na,:),v(:,1:na),d(1:na),lbwa)
  call upper_to_bv(a(1:na,1:na),bv_na_d,lbwa, tol1,error)
  call bv_times_sweeps(bv_na_d,sw_d,ub_na_d,error)
  call ub_to_upper(ub_na_d,a(1:na,1:na),error)
  a0(1:na,1:na)=a(1:na,1:na)
  call cpu_time(t0)
  call recompress_ub_to_bv(ub_na_d, bv_na_d,tol2,tol2,0,error)
  call cpu_time(t1)
  call bv_to_upper(bv_na_d,a1(1:na,1:na),error)
  test_name = "Real UB to BV (n=10; full ubw);"
  call d_output_result_upper(test_name,a0(1:na,1:na),a1(1:na,1:na),na/2,bv_na_d%ubw,t0,t1,tol2,error)
  call d_deallocate_ub(ub_na_d)
  call d_deallocate_bv(bv_na_d)
  call d_deallocate_sweeps(sw_d)

  !! BV to UB tests
  print *
  print *, "--------------------------------"
  print *
  print *, "Real BV to UB Recompression Tests:"
  print *
  ! Test 9; case 1 termination test.
  u=u0; v=v0; d=d0
  u(:,rmax-2:rmax)=0.1*tol2*u(:,rmax-2:rmax)
  u(1:n-rmax+1,rmax-3)=0.1_dp*tol2*u(1:n-rmax+1,rmax-3)
  call d_assemble_upper(a,u,v,d,lbw)
  a0=a
  call upper_to_bv(a,bv_d,lbw,tol1,error)
  call cpu_time(t0)
  call recompress_bv_to_ub(bv_d, ub_d,tol2,tol2,0,error)
  call cpu_time(t1)
  call ub_to_upper(ub_d,a1,error)
  test_name = "Real BV to UB (term. 1);"
  call d_output_result_upper(test_name,a0,a1,rmax-3,ub_d%ubw,t0,t1,tol2,error)
  ! Test 10; case 2 termination test.
  u=u0; v=v0; d=d0
  u(:,rmax-2:rmax)=0.1*tol2*u(:,rmax-2:rmax)
  call d_assemble_upper(a,u,v,d,lbw)
  a0=a
  call upper_to_bv(a,bv_d,lbw,tol1,error)
  call cpu_time(t0)
  call recompress_bv_to_ub(bv_d, ub_d,tol2,tol2,0,error)
  call cpu_time(t1)
  call ub_to_upper(ub_d,a1,error)
  test_name = "Real BV to UB (term. 2);"
  call d_output_result_upper(test_name,a0,a1,rmax-3,ub_d%ubw,t0,t1,tol2,error)
  ! Test 11;
  u=u0; v=v0; d=d0
  u(:,rmax-2:rmax)=0.1*tol2*u(:,rmax-2:rmax)
  u(1:n-rmax+2,rmax-3)=0.1_dp*tol2*u(1:n-rmax+2,rmax-3)
  call d_assemble_upper(a,u,v,d,lbw)
  a0=a
  call upper_to_bv(a,bv_d,lbw,tol1,error)
  call cpu_time(t0)
  call recompress_bv_to_ub(bv_d, ub_d,tol2,tol2,0,error)
  call cpu_time(t1)
  call ub_to_upper(ub_d,a1,error)
  test_name = "Real BV to UB (term. 3);"
  call d_output_result_upper(test_name,a0,a1,rmax-3,ub_d%ubw,t0,t1,tol2,error)
  ! Test 12; recompress to one test.  
  u=u0; v=v0; d=d0
  u(:,2:rmax)=0.1*tol2*u(:,2:rmax)
  call d_assemble_upper(a,u,v,d,lbw)
  a0=a
  call upper_to_bv(a,bv_d,lbw,tol1,error)
  call cpu_time(t0)
  call recompress_bv_to_ub(bv_d, ub_d,tol2,tol2,0,error)
  call cpu_time(t1)
  call ub_to_upper(ub_d,a1,error)
  test_name = "Real BV to UB (recompress to 1);"
  call d_output_result_upper(test_name,a0,a1,1,ub_d%ubw,t0,t1,tol2,error)
  ! Test 13; recompress to zero test.  
  u=u0; v=v0; d=d0
  u=0.1*tol2*u
  call d_assemble_upper(a,u,v,d,lbw)
  a0=a
  call upper_to_bv(a,bv_d,lbw,tol1,error)
  call cpu_time(t0)
  call recompress_bv_to_ub(bv_d, ub_d,tol2,tol2,0,error)
  call cpu_time(t1)
  call ub_to_upper(ub_d,a1,error)
  test_name = "Real BV to UB (recompress to 0);"
  call d_output_result_upper(test_name,a0,a1,0,ub_d%ubw,t0,t1,tol2,error)
  ! Test 14; no recompression
  u=u0; v=v0; d=d0
  call d_assemble_upper(a,u,v,d,lbw)
  a0=a
  call upper_to_bv(a,bv_d,lbw,tol1,error)
  call cpu_time(t0)
  call recompress_bv_to_ub(bv_d, ub_d,tol2,tol2,0,error)
  call cpu_time(t1)
  call ub_to_upper(ub_d,a1,error)
  test_name = "Real BV to UB (no recompression);"
  call d_output_result_upper(test_name,a0,a1,rmax,ub_d%ubw,t0,t1,tol2,error)
  ! small matrices
  na=1
  lbwa=min(na-1,lbw)
  ub_na_d=d_new_ub(na,lbwmax,ubwmax)
  bv_na_d=d_new_bv(na,lbwmax,ubwmax)
  u=u0; v=v0; d=d0
  call d_assemble_upper(a(1:na,1:na),u(1:na,:),v(:,1:na),d(1:na),lbwa)
  a0(1:na,1:na)=a(1:na,1:na)
  call upper_to_bv(a(1:na,1:na),bv_na_d,lbwa, tol1,error)
  call cpu_time(t0)
  call recompress_bv_to_ub(bv_na_d, ub_na_d,tol2,tol2,0,error)
  call cpu_time(t1)
  call ub_to_upper(ub_na_d,a1(1:na,1:na),error)
  test_name = "Real BV to UB (n=1);"
  call d_output_result_upper(test_name,a0(1:na,1:na),a1(1:na,1:na),0,ub_na_d%ubw,t0,t1,tol2,error)
  call d_deallocate_ub(ub_na_d)
  call d_deallocate_bv(bv_na_d)
  na=2
  lbwa=min(na-1,lbw)
  ub_na_d=d_new_ub(na,lbwmax,ubwmax)
  bv_na_d=d_new_bv(na,lbwmax,ubwmax)
  u=u0; v=v0; d=d0
  call d_assemble_upper(a(1:na,1:na),u(1:na,:),v(:,1:na),d(1:na),lbwa)
  a0(1:na,1:na)=a(1:na,1:na)
  call upper_to_bv(a(1:na,1:na),bv_na_d,lbwa, tol1,error)
  call cpu_time(t0)
  call recompress_bv_to_ub(bv_na_d, ub_na_d,tol2,tol2,0,error)
  call cpu_time(t1)
  call ub_to_upper(ub_na_d,a1(1:na,1:na),error)
  test_name = "Real BV to UB (n=2);"
  call d_output_result_upper(test_name,a0(1:na,1:na),a1(1:na,1:na),1,ub_na_d%ubw,t0,t1,tol2,error)
  call d_deallocate_ub(ub_na_d)
  call d_deallocate_bv(bv_na_d)
  na=3
  lbwa=min(na-1,lbw)
  ub_na_d=d_new_ub(na,lbwmax,ubwmax)
  bv_na_d=d_new_bv(na,lbwmax,ubwmax)
  u=u0; v=v0; d=d0
  call d_assemble_upper(a(1:na,1:na),u(1:na,:),v(:,1:na),d(1:na),lbwa)
  a0(1:na,1:na)=a(1:na,1:na)
  call upper_to_bv(a(1:na,1:na),bv_na_d,lbwa, tol1,error)
  call cpu_time(t0)
  call recompress_bv_to_ub(bv_na_d, ub_na_d,tol2,tol2,0,error)
  call cpu_time(t1)
  call ub_to_upper(ub_na_d,a1(1:na,1:na),error)
  test_name = "Real BV to UB (n=3);"
  call d_output_result_upper(test_name,a0(1:na,1:na),a1(1:na,1:na),1,ub_na_d%ubw,t0,t1,tol2,error)
  call d_deallocate_ub(ub_na_d)
  call d_deallocate_bv(bv_na_d)

  na=4
  lbwa=min(na-1,lbw)
  ub_na_d=d_new_ub(na,lbwmax,ubwmax)
  bv_na_d=d_new_bv(na,lbwmax,ubwmax)
  u=u0; v=v0; d=d0
  call d_assemble_upper(a(1:na,1:na),u(1:na,:),v(:,1:na),d(1:na),lbwa)
  a0(1:na,1:na)=a(1:na,1:na)
  call upper_to_bv(a(1:na,1:na),bv_na_d,lbwa, tol1,error)

  call cpu_time(t0)
  call recompress_bv_to_ub(bv_na_d, ub_na_d,tol2,tol2,0,error)
  call cpu_time(t1)
  call ub_to_upper(ub_na_d,a1(1:na,1:na),error)
  test_name = "Real BV to UB (n=4);"
  call d_output_result_upper(test_name,a0(1:na,1:na),a1(1:na,1:na),2,ub_na_d%ubw,t0,t1,tol2,error)
  call d_deallocate_ub(ub_na_d)
  call d_deallocate_bv(bv_na_d)
  !
  na=10; lbwa=6; numsweeps=7
  ub_na_d=d_new_ub(na,lbwmax,ubwmax)
  bv_na_d=d_new_bv(na,lbwmax,ubwmax)
  sw_d=d_random_sweeps(na,numsweeps)
  u=u0; v=v0; d=d0
  call d_assemble_upper(a(1:na,1:na),u(1:na,:),v(:,1:na),d(1:na),lbwa)
  call upper_to_ub(a(1:na,1:na),ub_na_d,lbwa, tol1,error)
  call sweeps_times_ub(sw_d,ub_na_d,bv_na_d,error)
  call bv_to_upper(bv_na_d,a(1:na,1:na),error)
  a0(1:na,1:na)=a(1:na,1:na)
  call cpu_time(t0)
  call recompress_bv_to_ub(bv_na_d, ub_na_d,tol2,tol2,0,error)
  call cpu_time(t1)
  call ub_to_upper(ub_na_d,a1(1:na,1:na),error)
  test_name = "Real BV to UB (n=10; full ubw);"
  call d_output_result_upper(test_name,a0(1:na,1:na),a1(1:na,1:na),na/2,ub_na_d%ubw,t0,t1,tol2,error)
  call d_deallocate_ub(ub_na_d)
  call d_deallocate_bv(bv_na_d)
  call d_deallocate_sweeps(sw_d)

  !
  ! complex
  !
  call random_matrix(u0_c)
  call random_matrix(v0_c)
  call random_matrix(d0_c)
  u_c=u0_c; v_c=v0_c; d_c=d0_c

  ! !! UB to BV tests
  print *
  print *, "--------------------------------"
  print *
  print *, "Complex UB to BV Recompression Tests"
  print *
  ! case 1 termination test.
  u_c=u0_c; v_c=v0_c; d_c=d0_c
  v_c(rmax-2:rmax,:)=0.1*tol2*v_c(rmax-2:rmax,:)
  v_c(rmax-3,rmax:n)=0.1_dp*tol2*v_c(rmax-3,rmax:n)
  call c_assemble_upper(a_c,u_c,v_c,d_c,lbw)
  a0_c=a_c
  call upper_to_ub(a_c,ub_c,lbw,tol1,error)
  call cpu_time(t0)
  call recompress_ub_to_bv(ub_c, bv_c,tol2,tol2,0,error)
  call cpu_time(t1)
  call bv_to_upper(bv_c,a1_c,error)
  test_name = "Complex UB to BV (term. 1);"
  call c_output_result_upper(test_name,a0_c,a1_c,rmax-3,bv_c%ubw,t0,t1,tol2,error)
  ! 
  u_c=u0_c; v_c=v0_c; d_c=d0_c
  u_c(:,rmax-3:rmax)=0.1*tol2*u_c(:,rmax-3:rmax)
  call c_assemble_upper(a_c,u_c,v_c,d_c,lbw)
  a0_c=a_c
  call upper_to_ub(a_c,ub_c,lbw,tol1,error)
  call cpu_time(t0)
  call recompress_ub_to_bv(ub_c, bv_c,tol2,tol2,0,error)
  call cpu_time(t1)
  call bv_to_upper(bv_c,a1_c,error)
  test_name = "Complex UB to BV (term. 2);"
  call c_output_result_upper(test_name,a0_c,a1_c,rmax-4,bv_c%ubw,t0,t1,tol2,error)
  ! case 3 termination test
  u_c=u0_c; v_c=v0_c; d_c=d0_c
  v_c(rmax-2:rmax,:)=0.1*tol2*v_c(rmax-2:rmax,:)
  v_c(rmax-3,rmax-1:n)=0.1*tol2*v_c(rmax-3,rmax-1:n)
  call c_assemble_upper(a_c,u_c,v_c,d_c,lbw)
  a0_c=a_c
  call upper_to_ub(a_c,ub_c,lbw,tol1,error)
  call cpu_time(t0)
  call recompress_ub_to_bv(ub_c, bv_c,tol2,tol2,0,error)
  call cpu_time(t1)
  call bv_to_upper(bv_c,a1_c,error)
  test_name = "Complex UB to BV (term. 3);"
  call c_output_result_upper(test_name,a0_c,a1_c,rmax-3,bv_c%ubw,t0,t1,tol2,error)
  ! recompress to one test.
  u_c=u0_c; v_c=v0_c; d_c=d0_c
  u_c(:,2:rmax)=0.1*tol2*u_c(:,2:rmax)
  call c_assemble_upper(a_c,u_c,v_c,d_c,lbw)
  a0_c=a_c
  call upper_to_ub(a_c,ub_c,lbw,tol1,error)
  call cpu_time(t0)
  call recompress_ub_to_bv(ub_c, bv_c,tol2,tol2,0,error)
  call cpu_time(t1)
  call bv_to_upper(bv_c,a1_c,error)
  test_name = "Complex UB to BV (recompress to 1);"
  call c_output_result_upper(test_name,a0_c,a1_c,1,bv_c%ubw,t0,t1,tol2,error)
  ! recompress to zero test.
  u_c=u0_c; v_c=v0_c; d_c=d0_c
  u_c=0.05*tol2*u_c
  call c_assemble_upper(a_c,u_c,v_c,d_c,lbw)
  a0_c=a_c
  call upper_to_ub(a_c,ub_c,lbw,tol1,error)
  call cpu_time(t0)
  call recompress_ub_to_bv(ub_c, bv_c,tol2,tol2,0,error)
  call cpu_time(t1)
  call bv_to_upper(bv_c,a1_c,error)
  test_name = "Complex UB to BV (recompress to 0);"
  call c_output_result_upper(test_name,a0_c,a1_c,0,bv_c%ubw,t0,t1,tol2,error)
  ! test 8: no recompression
  u_c=u0_c; v_c=v0_c; d_c=d0_c
  call c_assemble_upper(a_c,u_c,v_c,d_c,lbw)
  a0_c=a_c
  call upper_to_ub(a_c,ub_c,lbw,tol1,error)
  call cpu_time(t0)
  call recompress_ub_to_bv(ub_c, bv_c,tol2,tol2,0,error)
  call cpu_time(t1)
  call bv_to_upper(bv_c,a1_c,error)
  test_name = "Complex UB to BV (no recompression);"
  call c_output_result_upper(test_name,a0_c,a1_c,rmax,bv_c%ubw,t0,t1,tol2,error)
  ! small matrices
  na=1
  lbwa=min(na-1,lbw)
  ub_na_c=c_new_ub(na,lbwmax,ubwmax)
  bv_na_c=c_new_bv(na,lbwmax,ubwmax)
  u_c=u0_c; v_c=v0_c; d_c=d0_c
  call c_assemble_upper(a_c(1:na,1:na),u_c(1:na,:),v_c(:,1:na),d_c(1:na),lbwa)
  a0_c(1:na,1:na)=a_c(1:na,1:na)
  call upper_to_ub(a_c(1:na,1:na),ub_na_c,lbwa, tol1,error)
  call cpu_time(t0)
  call recompress_ub_to_bv(ub_na_c, bv_na_c,tol2,tol2,0,error)
  call cpu_time(t1)
  call bv_to_upper(bv_na_c,a1_c(1:na,1:na),error)
  test_name = "Complex UB to BV (n=1);"
  call c_output_result_upper(test_name,a0_c(1:na,1:na),a1_c(1:na,1:na),0,bv_na_c%ubw,t0,t1,tol2,error)
  call c_deallocate_ub(ub_na_c)
  call c_deallocate_bv(bv_na_c)
  na=2
  lbwa=min(na-1,lbw)
  ub_na_c=c_new_ub(na,lbwmax,ubwmax)
  bv_na_c=c_new_bv(na,lbwmax,ubwmax)
  u_c=u0_c; v_c=v0_c; d_c=d0_c
  call c_assemble_upper(a_c(1:na,1:na),u_c(1:na,:),v_c(:,1:na),d_c(1:na),lbwa)
  a0_c(1:na,1:na)=a_c(1:na,1:na)
  call upper_to_ub(a_c(1:na,1:na),ub_na_c,lbwa, tol1,error)
  call cpu_time(t0)
  call recompress_ub_to_bv(ub_na_c, bv_na_c,tol2,tol2,0,error)
  call cpu_time(t1)
  call bv_to_upper(bv_na_c,a1_c(1:na,1:na),error)
  test_name = "Complex UB to BV (n=2);"
  call c_output_result_upper(test_name,a0_c(1:na,1:na),a1_c(1:na,1:na),1,bv_na_c%ubw,t0,t1,tol2,error)
  call c_deallocate_ub(ub_na_c)
  call c_deallocate_bv(bv_na_c)
  na=3
  lbwa=min(na-1,lbw)
  ub_na_c=c_new_ub(na,lbwmax,ubwmax)
  bv_na_c=c_new_bv(na,lbwmax,ubwmax)
  u_c=u0_c; v_c=v0_c; d_c=d0_c
  call c_assemble_upper(a_c(1:na,1:na),u_c(1:na,:),v_c(:,1:na),d_c(1:na),lbwa)
  a0_c(1:na,1:na)=a_c(1:na,1:na)
  call upper_to_ub(a_c(1:na,1:na),ub_na_c,lbwa, tol1,error)
  call cpu_time(t0)
  call recompress_ub_to_bv(ub_na_c, bv_na_c,tol2,tol2,0,error)
  call cpu_time(t1)
  call bv_to_upper(bv_na_c,a1_c(1:na,1:na),error)
  test_name = "Complex UB to BV (n=3);"
  call c_output_result_upper(test_name,a0_c(1:na,1:na),a1_c(1:na,1:na),1,bv_na_c%ubw,t0,t1,tol2,error)
  call c_deallocate_ub(ub_na_c)
  call c_deallocate_bv(bv_na_c)
  na=4
  lbwa=min(na-1,lbw)
  ub_na_c=c_new_ub(na,lbwmax,ubwmax)
  bv_na_c=c_new_bv(na,lbwmax,ubwmax)
  u_c=u0_c; v_c=v0_c; d_c=d0_c
  call c_assemble_upper(a_c(1:na,1:na),u_c(1:na,:),v_c(:,1:na),d_c(1:na),lbwa)
  a0_c(1:na,1:na)=a_c(1:na,1:na)
  call upper_to_ub(a_c(1:na,1:na),ub_na_c,lbwa, tol1,error)
  call cpu_time(t0)
  call recompress_ub_to_bv(ub_na_c, bv_na_c,tol2,tol2,0,error)
  call cpu_time(t1)
  call bv_to_upper(bv_na_c,a1_c(1:na,1:na),error)
  test_name = "Complex UB to BV (n=4);"
  call c_output_result_upper(test_name,a0_c(1:na,1:na),a1_c(1:na,1:na),2,bv_na_c%ubw,t0,t1,tol2,error)
  call c_deallocate_ub(ub_na_c)
  call c_deallocate_bv(bv_na_c)
  !

  na=10; lbwa=6; numsweeps=7
  sw_c=c_random_sweeps(na,numsweeps)
  ub_na_c=c_new_ub(na,lbwmax,ubwmax)
  bv_na_c=c_new_bv(na,lbwmax,ubwmax)
  u_c=u0_c; v_c=v0_c; d_c=d0_c
  call c_assemble_upper(a_c(1:na,1:na),u_c(1:na,:),v_c(:,1:na),d_c(1:na),lbwa)
  call upper_to_bv(a_c(1:na,1:na),bv_na_c,lbwa, tol1,error)
  call bv_times_sweeps(bv_na_c,sw_c,ub_na_c,error)
  call ub_to_upper(ub_na_c,a_c(1:na,1:na),error)
  a0_c(1:na,1:na)=a_c(1:na,1:na)
  call cpu_time(t0)
  call recompress_ub_to_bv(ub_na_c, bv_na_c,tol2,tol2,0,error)
  call cpu_time(t1)
  call bv_to_upper(bv_na_c,a1_c(1:na,1:na),error)
  test_name = "Complex UB to BV (n=10; full ubw);"
  call c_output_result_upper(test_name,a0_c(1:na,1:na),a1_c(1:na,1:na),na/2,bv_na_c%ubw,t0,t1,tol2,error)
  call c_deallocate_ub(ub_na_c)
  call c_deallocate_bv(bv_na_c)
  call c_deallocate_sweeps(sw_c)


  !! BV to UB tests


  print *
  print *, "--------------------------------"
  print *
  print *, "Complex BV to UB Recompression Tests"
  print *
  ! case 1 termination test.
  u_c=u0_c; v_c=v0_c; d_c=d0_c
  u_c(:,rmax-2:rmax)=0.1*tol2*u_c(:,rmax-2:rmax)
  u_c(1:n-rmax+1,rmax-3)=0.1_dp*tol2*u_c(1:n-rmax+1,rmax-3)
  call c_assemble_upper(a_c,u_c,v_c,d_c,lbw)
  a0_c=a_c
  call upper_to_bv(a_c,bv_c,lbw,tol1,error)
  call cpu_time(t0)
  call recompress_bv_to_ub(bv_c, ub_c,tol2,tol2,0,error)
  call cpu_time(t1)
  call ub_to_upper(ub_c,a1_c,error)
  test_name = "Complex BV to UB (term. 1);"
  call c_output_result_upper(test_name,a0_c,a1_c,rmax-3,ub_c%ubw,t0,t1,tol2,error)
  ! case 2 termination test.
  u_c=u0_c; v_c=v0_c; d_c=d0_c
  u_c(:,rmax-2:rmax)=0.1*tol2*u_c(:,rmax-2:rmax)
  call c_assemble_upper(a_c,u_c,v_c,d_c,lbw)
  a0_c=a_c
  call upper_to_bv(a_c,bv_c,lbw,tol1,error)
  call cpu_time(t0)
  call recompress_bv_to_ub(bv_c, ub_c,tol2,tol2,0,error)
  call cpu_time(t1)
  call ub_to_upper(ub_c,a1_c,error)
  test_name = "Complex BV to UB (term. 2);"
  call c_output_result_upper(test_name,a0_c,a1_c,rmax-3,ub_c%ubw,t0,t1,tol2,error)
  ! case 3 termination
  u_c=u0_c; v_c=v0_c; d_c=d0_c
  u_c(:,rmax-2:rmax)=0.1*tol2*u_c(:,rmax-2:rmax)
  u_c(1:n-rmax+2,rmax-3)=0.1_dp*tol2*u_c(1:n-rmax+2,rmax-3)
  call c_assemble_upper(a_c,u_c,v_c,d_c,lbw)
  a0_c=a_c
  call upper_to_bv(a_c,bv_c,lbw,tol1,error)
  call cpu_time(t0)
  call recompress_bv_to_ub(bv_c, ub_c,tol2,tol2,0,error)
  call cpu_time(t1)
  call ub_to_upper(ub_c,a1_c,error)
  test_name = "Complex BV to UB (term. 3);"
  call c_output_result_upper(test_name,a0_c,a1_c,rmax-3,ub_c%ubw,t0,t1,tol2,error)
  ! recompress to one test.  
  u_c=u0_c; v_c=v0_c; d_c=d0_c
  u_c(:,2:rmax)=0.1*tol2*u_c(:,2:rmax)
  call c_assemble_upper(a_c,u_c,v_c,d_c,lbw)
  a0_c=a_c
  call upper_to_bv(a_c,bv_c,lbw,tol1,error)
  call cpu_time(t0)
  call recompress_bv_to_ub(bv_c, ub_c,tol2,tol2,0,error)
  call cpu_time(t1)
  call ub_to_upper(ub_c,a1_c,error)
  test_name = "Complex BV to UB (recompress to 1);"
  call c_output_result_upper(test_name,a0_c,a1_c,1,ub_c%ubw,t0,t1,tol2,error)


  ! recompress to zero test.  
  u_c=u0_c; v_c=v0_c; d_c=d0_c
  u_c=0.05*tol2*u_c
  call c_assemble_upper(a_c,u_c,v_c,d_c,lbw)
  a0_c=a_c
  call upper_to_bv(a_c,bv_c,lbw,tol1,error)
  call cpu_time(t0)
  call recompress_bv_to_ub(bv_c, ub_c,tol2,tol2,0,error)
  call cpu_time(t1)
  call ub_to_upper(ub_c,a1_c,error)
  test_name = "Complex BV to UB (recompress to 0);"
  call c_output_result_upper(test_name,a0_c,a1_c,0,ub_c%ubw,t0,t1,tol2,error)
  ! no recompression
  u_c=u0_c; v_c=v0_c; d_c=d0_c
  call c_assemble_upper(a_c,u_c,v_c,d_c,lbw)
  a0_c=a_c
  call upper_to_bv(a_c,bv_c,lbw,tol1,error)
  call cpu_time(t0)
  call recompress_bv_to_ub(bv_c, ub_c,tol2,tol2,0,error)
  call cpu_time(t1)
  call ub_to_upper(ub_c,a1_c,error)
  test_name = "Complex BV to UB (no recompression);"
  call c_output_result_upper(test_name,a0_c,a1_c,rmax,ub_c%ubw,t0,t1,tol2,error)
  ! small matrices

  na=1
  lbwa=min(na-1,lbw)
  ub_na_c=c_new_ub(na,lbwmax,ubwmax)
  bv_na_c=c_new_bv(na,lbwmax,ubwmax)
  u_c=u0_c; v_c=v0_c; d_c=d0_c
  call c_assemble_upper(a_c(1:na,1:na),u_c(1:na,:),v_c(:,1:na),d_c(1:na),lbwa)
  a0_c(1:na,1:na)=a_c(1:na,1:na)
  call upper_to_bv(a_c(1:na,1:na),bv_na_c,lbwa, tol1,error)
  call cpu_time(t0)
  call recompress_bv_to_ub(bv_na_c, ub_na_c,tol2,tol2,0,error)
  call cpu_time(t1)
  call ub_to_upper(ub_na_c,a1_c(1:na,1:na),error)
  test_name = "Complex BV to UB (n=1);"
  call c_output_result_upper(test_name,a0_c(1:na,1:na),a1_c(1:na,1:na),0,ub_na_c%ubw,t0,t1,tol2,error)
  call c_deallocate_ub(ub_na_c)
  call c_deallocate_bv(bv_na_c)
  na=2
  lbwa=min(na-1,lbw)
  ub_na_c=c_new_ub(na,lbwmax,ubwmax)
  bv_na_c=c_new_bv(na,lbwmax,ubwmax)
  u_c=u0_c; v_c=v0_c; d_c=d0_c
  call c_assemble_upper(a_c(1:na,1:na),u_c(1:na,:),v_c(:,1:na),d_c(1:na),lbwa)
  a0_c(1:na,1:na)=a_c(1:na,1:na)
  call upper_to_bv(a_c(1:na,1:na),bv_na_c,lbwa, tol1,error)
  call cpu_time(t0)
  call recompress_bv_to_ub(bv_na_c, ub_na_c,tol2,tol2,0,error)
  call cpu_time(t1)
  call ub_to_upper(ub_na_c,a1_c(1:na,1:na),error)
  test_name = "Complex BV to UB (n=2);"
  call c_output_result_upper(test_name,a0_c(1:na,1:na),a1_c(1:na,1:na),1,ub_na_c%ubw,t0,t1,tol2,error)
  call c_deallocate_ub(ub_na_c)
  call c_deallocate_bv(bv_na_c)
  na=3
  lbwa=min(na-1,lbw)
  ub_na_c=c_new_ub(na,lbwmax,ubwmax)
  bv_na_c=c_new_bv(na,lbwmax,ubwmax)
  u_c=u0_c; v_c=v0_c; d_c=d0_c
  call c_assemble_upper(a_c(1:na,1:na),u_c(1:na,:),v_c(:,1:na),d_c(1:na),lbwa)
  a0_c(1:na,1:na)=a_c(1:na,1:na)
  call upper_to_bv(a_c(1:na,1:na),bv_na_c,lbwa, tol1,error)
  call cpu_time(t0)
  call recompress_bv_to_ub(bv_na_c, ub_na_c,tol2,tol2,0,error)
  call cpu_time(t1)
  call ub_to_upper(ub_na_c,a1_c(1:na,1:na),error)
  test_name = "Complex BV to UB (n=3);"
  call c_output_result_upper(test_name,a0_c(1:na,1:na),a1_c(1:na,1:na),1,ub_na_c%ubw,t0,t1,tol2,error)
  call c_deallocate_ub(ub_na_c)
  call c_deallocate_bv(bv_na_c)


  na=4
  lbwa=min(na-1,lbw)
  ub_na_c=c_new_ub(na,lbwmax,ubwmax)
  bv_na_c=c_new_bv(na,lbwmax,ubwmax)
  u_c=u0_c; v_c=v0_c; d_c=d0_c
  call c_assemble_upper(a_c(1:na,1:na),u_c(1:na,:),v_c(:,1:na),d_c(1:na),lbwa)
  a0_c(1:na,1:na)=a_c(1:na,1:na)
  call upper_to_bv(a_c(1:na,1:na),bv_na_c,lbwa, tol1,error)
  call cpu_time(t0)
  call recompress_bv_to_ub(bv_na_c, ub_na_c,tol2,tol2,0,error)
  call cpu_time(t1)
  call ub_to_upper(ub_na_c,a1_c(1:na,1:na),error)
  test_name = "Complex BV to UB (n=4);"
  call c_output_result_upper(test_name,a0_c(1:na,1:na),a1_c(1:na,1:na),2,ub_na_c%ubw,t0,t1,tol2,error)
  call c_deallocate_ub(ub_na_c)
  call c_deallocate_bv(bv_na_c)
  !
  na=10; lbwa=6; numsweeps=7
  ub_na_c=c_new_ub(na,lbwmax,ubwmax)
  bv_na_c=c_new_bv(na,lbwmax,ubwmax)
  sw_c=c_random_sweeps(na,numsweeps)
  u_c=u0_c; v_c=v0_c; d_c=d0_c
  call c_assemble_upper(a_c(1:na,1:na),u_c(1:na,:),v_c(:,1:na),d_c(1:na),lbwa)
  call upper_to_ub(a_c(1:na,1:na),ub_na_c,lbwa, tol1,error)
  call sweeps_times_ub(sw_c,ub_na_c,bv_na_c,error)
  call bv_to_upper(bv_na_c,a_c(1:na,1:na),error)
  a0_c(1:na,1:na)=a_c(1:na,1:na)
  call cpu_time(t0)
  call recompress_bv_to_ub(bv_na_c, ub_na_c,tol2,tol2,0,error)
  call cpu_time(t1)
  call ub_to_upper(ub_na_c,a1_c(1:na,1:na),error)
  test_name = "Complex BV to UB (n=10; full ubw);"
  call c_output_result_upper(test_name,a0_c(1:na,1:na),a1_c(1:na,1:na),na/2,ub_na_c%ubw,t0,t1,tol2,error)
  call c_deallocate_ub(ub_na_c)
  call c_deallocate_bv(bv_na_c)
  call c_deallocate_sweeps(sw_c)
  print *

end program test_recompress
