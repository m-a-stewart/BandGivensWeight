program test_convert
  use nested
  use test_data
  implicit none

  !
  ! Note: This does not cover large rank cases (e.g. ubw=n-1).  But the conversion
  ! routines are tested for these cases by test_sweeps.f90, since sweeps.f90 uses
  ! the conversion routines.
  !

  real(kind=dp) :: t0, t1
  integer(kind=int32) :: na, lbwa
  type(error_info) :: error
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
  ub_d=d_new_ub(n,lbwmax,ubwmax)
  ub_c=c_new_ub(n,lbwmax,ubwmax)
  bv_d=d_new_bv(n,lbwmax,ubwmax)
  bv_c=c_new_bv(n,lbwmax,ubwmax)

  ! real
  call random_seed
  call random_matrix(u)
  call random_matrix(v)
  call random_matrix(d)
  u0=u; v0=v; d0=d
  print *
  print *, "--------------------------------"
  print *
  print *, "Real UB and BV Conversion Tests:"
  print *
  ! ub to bv
  u=u0; v=v0; d=d0
  call d_assemble_a(a,u,v,d,lbw)
  a0=a
  call upper_to_ub(a,ub_d,lbw, tol,error)
  call cpu_time(t0)
  call convert_ub_to_bv(ub_d, bv_d,error)
  call cpu_time(t1)
  call bv_to_upper(bv_d,a1,error)
  test_name = "Real UB to BV;"
  call d_output_result(test_name,a0,a1,rmax,bv_d%ubw,t0,t1,tol2,error)
  !
  na=1
  lbwa=min(na-1,lbw)
  ub_na_d=d_new_ub(na,lbwmax,ubwmax)
  bv_na_d=d_new_bv(na,lbwmax,ubwmax)
  u=u0; v=v0; d=d0
  call d_assemble_a(a(1:na,1:na),u(1:na,:),v(:,1:na),d(1:na),lbwa)
  a0(1:na,1:na)=a(1:na,1:na)
  call upper_to_ub(a(1:na,1:na),ub_na_d,lbwa, tol1,error)
  call cpu_time(t0)
  call convert_ub_to_bv(ub_na_d, bv_na_d,error)
  call cpu_time(t1)
  call bv_to_upper(bv_na_d,a1(1:na,1:na),error)
  test_name = "Real UB to BV (n=1);"
  call d_output_result(test_name,a0(1:na,1:na),a1(1:na,1:na),0,bv_na_d%ubw,t0,t1,tol2,error)
  call d_deallocate_ub(ub_na_d)
  call d_deallocate_bv(bv_na_d)
  !
  na=2
  lbwa=min(na-1,lbw)
  ub_na_d=d_new_ub(na,lbwmax,ubwmax)
  bv_na_d=d_new_bv(na,lbwmax,ubwmax)
  u=u0; v=v0; d=d0
  call d_assemble_a(a(1:na,1:na),u(1:na,:),v(:,1:na),d(1:na),lbwa)
  a0(1:na,1:na)=a(1:na,1:na)
  call upper_to_ub(a(1:na,1:na),ub_na_d,lbwa, tol1,error)
  call cpu_time(t0)
  call convert_ub_to_bv(ub_na_d, bv_na_d,error)
  call cpu_time(t1)
  call bv_to_upper(bv_na_d,a1(1:na,1:na),error)
  test_name = "Real UB to BV (n=2);"
  call d_output_result(test_name,a0(1:na,1:na),a1(1:na,1:na),1,bv_na_d%ubw,t0,t1,tol2,error)
  call d_deallocate_ub(ub_na_d)
  call d_deallocate_bv(bv_na_d)
  !
  na=3
  lbwa=min(na-1,lbw)
  ub_na_d=d_new_ub(na,lbwmax,ubwmax)
  bv_na_d=d_new_bv(na,lbwmax,ubwmax)
  u=u0; v=v0; d=d0
  call d_assemble_a(a(1:na,1:na),u(1:na,:),v(:,1:na),d(1:na),lbwa)
  a0(1:na,1:na)=a(1:na,1:na)
  call upper_to_ub(a(1:na,1:na),ub_na_d,lbwa, tol1,error)
  call cpu_time(t0)
  call convert_ub_to_bv(ub_na_d, bv_na_d,error)
  call cpu_time(t1)
  call bv_to_upper(bv_na_d,a1(1:na,1:na),error)
  test_name = "Real UB to BV (n=3);"
  call d_output_result(test_name,a0(1:na,1:na),a1(1:na,1:na),1,bv_na_d%ubw,t0,t1,tol2,error)
  call d_deallocate_ub(ub_na_d)
  call d_deallocate_bv(bv_na_d)
  !
  na=4
  lbwa=min(na-1,lbw)
  ub_na_d=d_new_ub(na,lbwmax,ubwmax)
  bv_na_d=d_new_bv(na,lbwmax,ubwmax)
  u=u0; v=v0; d=d0
  call d_assemble_a(a(1:na,1:na),u(1:na,:),v(:,1:na),d(1:na),lbwa)
  a0(1:na,1:na)=a(1:na,1:na)
  call upper_to_ub(a(1:na,1:na),ub_na_d,lbwa, tol1,error)
  call cpu_time(t0)
  call convert_ub_to_bv(ub_na_d, bv_na_d,error)
  call cpu_time(t1)
  call bv_to_upper(bv_na_d,a1(1:na,1:na),error)
  test_name = "Real UB to BV (n=4);"
  call d_output_result(test_name,a0(1:na,1:na),a1(1:na,1:na),2,bv_na_d%ubw,t0,t1,tol2,error)
  call d_deallocate_ub(ub_na_d)
  call d_deallocate_bv(bv_na_d)
  ! bv to ub
  u=u0; v=v0; d=d0
  call d_assemble_a(a,u,v,d,lbw)
  a0=a
  call upper_to_bv(a,bv_d,lbw,tol,error)
  call cpu_time(t0)
  call convert_bv_to_ub(bv_d, ub_d, error)
  call cpu_time(t1)
  call ub_to_upper(ub_d,a1,error)
  test_name="Real BV to UB;"
  call d_output_result(test_name,a0,a1,rmax,ub_d%ubw,t0,t1,tol2,error)
  !
  na=1
  lbwa=min(na-1,lbw)
  ub_na_d=d_new_ub(na,lbwmax,ubwmax)
  bv_na_d=d_new_bv(na,lbwmax,ubwmax)
  u=u0; v=v0; d=d0
  call d_assemble_a(a(1:na,1:na),u(1:na,:),v(:,1:na),d(1:na),lbwa)
  a0(1:na,1:na)=a(1:na,1:na)
  call upper_to_bv(a(1:na,1:na),bv_na_d, lbwa, tol1, error)
  call cpu_time(t0)
  call convert_bv_to_ub(bv_na_d, ub_na_d,error)
  call cpu_time(t1)
  call ub_to_upper(ub_na_d,a1(1:na,1:na),error)
  test_name = "Real BV to UB (n=1);"
  call d_output_result(test_name,a0(1:na,1:na),a1(1:na,1:na),0,ub_na_d%ubw,t0,t1,tol2,error)
  call d_deallocate_ub(ub_na_d)
  call d_deallocate_bv(bv_na_d)
  !
  na=2
  lbwa=min(na-1,lbw)
  ub_na_d=d_new_ub(na,lbwmax,ubwmax)
  bv_na_d=d_new_bv(na,lbwmax,ubwmax)
  u=u0; v=v0; d=d0
  call d_assemble_a(a(1:na,1:na),u(1:na,:),v(:,1:na),d(1:na),lbwa)
  a0(1:na,1:na)=a(1:na,1:na)
  call upper_to_bv(a(1:na,1:na),bv_na_d, lbwa, tol1, error)
  call cpu_time(t0)
  call convert_bv_to_ub(bv_na_d, ub_na_d,error)
  call cpu_time(t1)
  call ub_to_upper(ub_na_d,a1(1:na,1:na),error)
  test_name = "Real BV to UB (n=2);"
  call d_output_result(test_name,a0(1:na,1:na),a1(1:na,1:na),1,ub_na_d%ubw,t0,t1,tol2,error)
  call d_deallocate_ub(ub_na_d)
  call d_deallocate_bv(bv_na_d)
  !
  na=3
  lbwa=min(na-1,lbw)
  ub_na_d=d_new_ub(na,lbwmax,ubwmax)
  bv_na_d=d_new_bv(na,lbwmax,ubwmax)
  u=u0; v=v0; d=d0
  call d_assemble_a(a(1:na,1:na),u(1:na,:),v(:,1:na),d(1:na),lbwa)
  a0(1:na,1:na)=a(1:na,1:na)
  call upper_to_bv(a(1:na,1:na),bv_na_d, lbwa, tol1, error)
  call cpu_time(t0)
  call convert_bv_to_ub(bv_na_d, ub_na_d,error)
  call cpu_time(t1)
  call ub_to_upper(ub_na_d,a1(1:na,1:na),error)
  test_name = "Real BV to UB (n=3);"
  call d_output_result(test_name,a0(1:na,1:na),a1(1:na,1:na),1,ub_na_d%ubw,t0,t1,tol2,error)
  call d_deallocate_ub(ub_na_d)
  call d_deallocate_bv(bv_na_d)
  !
  na=4
  lbwa=min(na-1,lbw)
  ub_na_d=d_new_ub(na,lbwmax,ubwmax)
  bv_na_d=d_new_bv(na,lbwmax,ubwmax)
  u=u0; v=v0; d=d0
  call d_assemble_a(a(1:na,1:na),u(1:na,:),v(:,1:na),d(1:na),lbwa)
  a0(1:na,1:na)=a(1:na,1:na)
  call upper_to_bv(a(1:na,1:na),bv_na_d, lbwa, tol1, error)
  call cpu_time(t0)
  call convert_bv_to_ub(bv_na_d, ub_na_d,error)
  call cpu_time(t1)
  call ub_to_upper(ub_na_d,a1(1:na,1:na),error)
  test_name = "Real BV to UB (n=4);"
  call d_output_result(test_name,a0(1:na,1:na),a1(1:na,1:na),2,ub_na_d%ubw,t0,t1,tol2,error)
  call d_deallocate_ub(ub_na_d)
  call d_deallocate_bv(bv_na_d)

  print *
  print *, "--------------------------------"
  print *
  print *, "Complex UB and BV Conversion Tests:"
  print *
  call random_matrix(u_c)
  call random_matrix(v_c)
  call random_matrix(d_c)
  u0_c=u_c; v0_c=v_c; d0_c=d_c
  call c_assemble_a(a_c, u_c, v_c, d_c, lbw)
  a0_c=a_c
  call upper_to_ub(a_c, ub_c, lbw, tol, error)
  call cpu_time(t0)
  call convert_ub_to_bv(ub_c, bv_c, error)
  call cpu_time(t1)
  call bv_to_upper(bv_c, a1_c, error)
  test_name="Complex UB to BV;"
  call c_output_result(test_name,a0_c,a1_c,rmax,bv_c%ubw,t0,t1,tol2,error)
  !
  na=1
  lbwa=min(na-1,lbw)
  ub_na_c=c_new_ub(na,lbwmax,ubwmax)
  bv_na_c=c_new_bv(na,lbwmax,ubwmax)
  u_c=u0_c; v_c=v0_c; d_c=d0_c
  call c_assemble_a(a_c(1:na,1:na),u_c(1:na,:),v_c(:,1:na),d_c(1:na),lbwa)
  a0_c(1:na,1:na)=a_c(1:na,1:na)
  call upper_to_ub(a_c(1:na,1:na),ub_na_c,lbwa, tol1,error)
  call cpu_time(t0)
  call convert_ub_to_bv(ub_na_c, bv_na_c,error)
  call cpu_time(t1)
  call bv_to_upper(bv_na_c,a1_c(1:na,1:na),error)
  test_name = "Complex UB to BV (n=1);"
  call c_output_result(test_name,a0_c(1:na,1:na),a1_c(1:na,1:na),0,bv_na_c%ubw,t0,t1,tol2,error)
  call c_deallocate_ub(ub_na_c)
  call c_deallocate_bv(bv_na_c)
  !
  na=2
  lbwa=min(na-1,lbw)
  ub_na_c=c_new_ub(na,lbwmax,ubwmax)
  bv_na_c=c_new_bv(na,lbwmax,ubwmax)
  u_c=u0_c; v_c=v0_c; d_c=d0_c
  call c_assemble_a(a_c(1:na,1:na),u_c(1:na,:),v_c(:,1:na),d_c(1:na),lbwa)
  a0_c(1:na,1:na)=a_c(1:na,1:na)
  call upper_to_ub(a_c(1:na,1:na),ub_na_c,lbwa, tol1,error)
  call cpu_time(t0)
  call convert_ub_to_bv(ub_na_c, bv_na_c,error)
  call cpu_time(t1)
  call bv_to_upper(bv_na_c,a1_c(1:na,1:na),error)
  test_name = "Complex UB to BV (n=2);"
  call c_output_result(test_name,a0_c(1:na,1:na),a1_c(1:na,1:na),1,bv_na_c%ubw,t0,t1,tol2,error)
  call c_deallocate_ub(ub_na_c)
  call c_deallocate_bv(bv_na_c)
  !
  na=3
  lbwa=min(na-1,lbw)
  ub_na_c=c_new_ub(na,lbwmax,ubwmax)
  bv_na_c=c_new_bv(na,lbwmax,ubwmax)
  u_c=u0_c; v_c=v0_c; d_c=d0_c
  call c_assemble_a(a_c(1:na,1:na),u_c(1:na,:),v_c(:,1:na),d_c(1:na),lbwa)
  a0_c(1:na,1:na)=a_c(1:na,1:na)
  call upper_to_ub(a_c(1:na,1:na),ub_na_c,lbwa, tol1,error)
  call cpu_time(t0)
  call convert_ub_to_bv(ub_na_c, bv_na_c,error)
  call cpu_time(t1)
  call bv_to_upper(bv_na_c,a1_c(1:na,1:na),error)
  test_name = "Complex UB to BV (n=3);"
  call c_output_result(test_name,a0_c(1:na,1:na),a1_c(1:na,1:na),1,bv_na_c%ubw,t0,t1,tol2,error)
  call c_deallocate_ub(ub_na_c)
  call c_deallocate_bv(bv_na_c)
  !
  na=4
  lbwa=min(na-1,lbw)
  ub_na_c=c_new_ub(na,lbwmax,ubwmax)
  bv_na_c=c_new_bv(na,lbwmax,ubwmax)
  u_c=u0_c; v_c=v0_c; d_c=d0_c
  call c_assemble_a(a_c(1:na,1:na),u_c(1:na,:),v_c(:,1:na),d_c(1:na),lbwa)
  a0_c(1:na,1:na)=a_c(1:na,1:na)
  call upper_to_ub(a_c(1:na,1:na),ub_na_c,lbwa, tol1,error)
  call cpu_time(t0)
  call convert_ub_to_bv(ub_na_c, bv_na_c,error)
  call cpu_time(t1)
  call bv_to_upper(bv_na_c,a1_c(1:na,1:na),error)
  test_name = "Complex UB to BV (n=4);"
  call c_output_result(test_name,a0_c(1:na,1:na),a1_c(1:na,1:na),2,bv_na_c%ubw,t0,t1,tol2,error)
  call c_deallocate_ub(ub_na_c)
  call c_deallocate_bv(bv_na_c)
  ! bv to ub
  u_c=u0_c; v_c=v0_c; d_c=d0_c
  call c_assemble_a(a_c, u_c, v_c, d_c, lbw)
  a0_c=a_c
  call upper_to_bv(a_c, bv_c, lbw, tol, error)
  call cpu_time(t0)
  call convert_bv_to_ub(bv_c, ub_c, error)
  call cpu_time(t1)
  call ub_to_upper(ub_c, a1_c, error)
  test_name="Complex BV to UB;"
  call c_output_result(test_name,a0_c,a1_c,rmax,ub_c%ubw,t0,t1,tol2,error)
  !
  na=1
  lbwa=min(na-1,lbw)
  ub_na_c=c_new_ub(na,lbwmax,ubwmax)
  bv_na_c=c_new_bv(na,lbwmax,ubwmax)
  u_c=u0_c; v_c=v0_c; d_c=d0_c
  call c_assemble_a(a_c(1:na,1:na),u_c(1:na,:),v_c(:,1:na),d_c(1:na),lbwa)
  a0_c(1:na,1:na)=a_c(1:na,1:na)
  call upper_to_bv(a_c(1:na,1:na),bv_na_c,lbwa, tol1,error)
  call cpu_time(t0)
  call convert_bv_to_ub(bv_na_c, ub_na_c,error)
  call cpu_time(t1)
  call ub_to_upper(ub_na_c,a1_c(1:na,1:na),error)
  test_name = "Complex BV to UB (n=1);"
  call c_output_result(test_name,a0_c(1:na,1:na),a1_c(1:na,1:na),0,ub_na_c%ubw,t0,t1,tol2,error)
  call c_deallocate_ub(ub_na_c)
  call c_deallocate_bv(bv_na_c)
  !
  na=2
  lbwa=min(na-1,lbw)
  ub_na_c=c_new_ub(na,lbwmax,ubwmax)
  bv_na_c=c_new_bv(na,lbwmax,ubwmax)
  u_c=u0_c; v_c=v0_c; d_c=d0_c
  call c_assemble_a(a_c(1:na,1:na),u_c(1:na,:),v_c(:,1:na),d_c(1:na),lbwa)
  a0_c(1:na,1:na)=a_c(1:na,1:na)
  call upper_to_bv(a_c(1:na,1:na),bv_na_c,lbwa, tol1,error)
  call cpu_time(t0)
  call convert_bv_to_ub(bv_na_c, ub_na_c,error)
  call cpu_time(t1)
  call ub_to_upper(ub_na_c,a1_c(1:na,1:na),error)
  test_name = "Complex BV to UB (n=2);"
  call c_output_result(test_name,a0_c(1:na,1:na),a1_c(1:na,1:na),1,ub_na_c%ubw,t0,t1,tol2,error)
  call c_deallocate_ub(ub_na_c)
  call c_deallocate_bv(bv_na_c)
  !
  na=3
  lbwa=min(na-1,lbw)
  ub_na_c=c_new_ub(na,lbwmax,ubwmax)
  bv_na_c=c_new_bv(na,lbwmax,ubwmax)
  u_c=u0_c; v_c=v0_c; d_c=d0_c
  call c_assemble_a(a_c(1:na,1:na),u_c(1:na,:),v_c(:,1:na),d_c(1:na),lbwa)
  a0_c(1:na,1:na)=a_c(1:na,1:na)
  call upper_to_bv(a_c(1:na,1:na),bv_na_c,lbwa, tol1,error)
  call cpu_time(t0)
  call convert_bv_to_ub(bv_na_c, ub_na_c,error)
  call cpu_time(t1)
  call ub_to_upper(ub_na_c,a1_c(1:na,1:na),error)
  test_name = "Complex BV to UB (n=3);"
  call c_output_result(test_name,a0_c(1:na,1:na),a1_c(1:na,1:na),1,ub_na_c%ubw,t0,t1,tol2,error)
  call c_deallocate_ub(ub_na_c)
  call c_deallocate_bv(bv_na_c)
  !
  na=4
  lbwa=min(na-1,lbw)
  ub_na_c=c_new_ub(na,lbwmax,ubwmax)
  bv_na_c=c_new_bv(na,lbwmax,ubwmax)
  u_c=u0_c; v_c=v0_c; d_c=d0_c
  call c_assemble_a(a_c(1:na,1:na),u_c(1:na,:),v_c(:,1:na),d_c(1:na),lbwa)
  a0_c(1:na,1:na)=a_c(1:na,1:na)
  call upper_to_bv(a_c(1:na,1:na),bv_na_c,lbwa, tol1,error)
  call cpu_time(t0)
  call convert_bv_to_ub(bv_na_c, ub_na_c,error)
  call cpu_time(t1)
  call ub_to_upper(ub_na_c,a1_c(1:na,1:na),error)
  test_name = "Complex BV to UB (n=4);"
  call c_output_result(test_name,a0_c(1:na,1:na),a1_c(1:na,1:na),2,ub_na_c%ubw,t0,t1,tol2,error)
  call c_deallocate_ub(ub_na_c)
  call c_deallocate_bv(bv_na_c)
  print *

end program test_convert
