program test_convert_wb_and_bt
  use mod_orth_rank
  use mod_test_data
  implicit none

  !
  ! Note: This does not cover large rank cases (e.g. ubw=n-1).  But the conversion
  ! routines are tested for these cases by test_sweeps1.f90, since sweeps1.f90 uses
  ! the conversion routines.
  !

  real(kind=dp) :: t0, t1
  integer(kind=int32) :: na, ubwa
  type(error_info) :: error
  integer, parameter :: n=50, rmax=13, lbwmax=rmax+5, ubw=2, ubwmax=10
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

  type(d_wb) :: wb_d, wb_na_d
  type(c_wb) :: wb_c, wb_na_c
  type(d_bt) :: bt_d, bt_na_d
  type(c_bt) :: bt_c, bt_na_c
  wb_d=d_new_wb(n,lbwmax,ubwmax)
  wb_c=c_new_wb(n,lbwmax,ubwmax)
  bt_d=d_new_bt(n,lbwmax,ubwmax)
  bt_c=c_new_bt(n,lbwmax,ubwmax)

  ! real
  call random_seed
  call random_matrix(u)
  call random_matrix(v)
  call random_matrix(d)
  u0=u; v0=v; d0=d
  print *
  print *, "--------------------------------"
  print *
  print *, "Real WB and BT Conversion Tests:"
  print *
  ! wb to bt
  u=u0; v=v0; d=d0
  call d_assemble_lower(a,u,v,d,ubw)
  a0=a
  call lower_to_wb(a,wb_d, ubw, tol,error)
  call cpu_time(t0)
  call convert_wb_to_bt(wb_d, bt_d,error)
  call cpu_time(t1)
  call bt_to_lower(bt_d,a1,error)
  test_name = "Real WB to BT;"
  call d_output_result_lower(test_name,a0,a1,rmax,bt_d%lbw,t0,t1,tol2,error)
  !
  na=1
  ubwa=min(na-1,ubw)
  wb_na_d=d_new_wb(na,lbwmax,ubwmax)
  bt_na_d=d_new_bt(na,lbwmax,ubwmax)
  u=u0; v=v0; d=d0
  call d_assemble_lower(a(1:na,1:na),u(1:na,:),v(:,1:na),d(1:na),ubwa)
  a0(1:na,1:na)=a(1:na,1:na)
  call lower_to_wb(a(1:na,1:na),wb_na_d,ubwa, tol1,error)
  call cpu_time(t0)
  call convert_wb_to_bt(wb_na_d, bt_na_d,error)
  call cpu_time(t1)
  call bt_to_lower(bt_na_d,a1(1:na,1:na),error)
  test_name = "Real WB to BT (n=1);"
  call d_output_result_lower(test_name,a0(1:na,1:na),a1(1:na,1:na),0,bt_na_d%lbw,t0,t1,tol2,error)
  call d_deallocate_wb(wb_na_d)
  call d_deallocate_bt(bt_na_d)
  !
  na=2
  ubwa=min(na-1,ubw)
  wb_na_d=d_new_wb(na,lbwmax,ubwmax)
  bt_na_d=d_new_bt(na,lbwmax,ubwmax)
  u=u0; v=v0; d=d0
  call d_assemble_lower(a(1:na,1:na),u(1:na,:),v(:,1:na),d(1:na),ubwa)
  a0(1:na,1:na)=a(1:na,1:na)
  call lower_to_wb(a(1:na,1:na),wb_na_d,ubwa, tol1,error)
  call cpu_time(t0)
  call convert_wb_to_bt(wb_na_d, bt_na_d,error)
  call cpu_time(t1)
  call bt_to_lower(bt_na_d,a1(1:na,1:na),error)
  test_name = "Real WB to BT (n=2);"
  call d_output_result_lower(test_name,a0(1:na,1:na),a1(1:na,1:na),1,bt_na_d%lbw,t0,t1,tol2,error)
  call d_deallocate_wb(wb_na_d)
  call d_deallocate_bt(bt_na_d)
  !
  na=3
  ubwa=min(na-1,ubw)
  wb_na_d=d_new_wb(na,lbwmax,ubwmax)
  bt_na_d=d_new_bt(na,lbwmax,ubwmax)
  u=u0; v=v0; d=d0
  call d_assemble_lower(a(1:na,1:na),u(1:na,:),v(:,1:na),d(1:na),ubwa)
  a0(1:na,1:na)=a(1:na,1:na)
  call lower_to_wb(a(1:na,1:na),wb_na_d,ubwa, tol1,error)
  call cpu_time(t0)
  call convert_wb_to_bt(wb_na_d, bt_na_d,error)
  call cpu_time(t1)
  call bt_to_lower(bt_na_d,a1(1:na,1:na),error)
  test_name = "Real WB to BT (n=3);"
  call d_output_result_lower(test_name,a0(1:na,1:na),a1(1:na,1:na),1,bt_na_d%lbw,t0,t1,tol2,error)
  call d_deallocate_wb(wb_na_d)
  call d_deallocate_bt(bt_na_d)
  !
  na=4
  ubwa=min(na-1,ubw)
  wb_na_d=d_new_wb(na,lbwmax,ubwmax)
  bt_na_d=d_new_bt(na,lbwmax,ubwmax)
  u=u0; v=v0; d=d0
  call d_assemble_lower(a(1:na,1:na),u(1:na,:),v(:,1:na),d(1:na),ubwa)
  a0(1:na,1:na)=a(1:na,1:na)
  call lower_to_wb(a(1:na,1:na),wb_na_d,ubwa, tol1,error)
  call cpu_time(t0)
  call convert_wb_to_bt(wb_na_d, bt_na_d,error)
  call cpu_time(t1)
  call bt_to_lower(bt_na_d,a1(1:na,1:na),error)
  test_name = "Real WB to BT (n=4);"
  call d_output_result_lower(test_name,a0(1:na,1:na),a1(1:na,1:na),2,bt_na_d%lbw,t0,t1,tol2,error)
  call d_deallocate_wb(wb_na_d)
  call d_deallocate_bt(bt_na_d)
  
  ! bt to wb
  u=u0; v=v0; d=d0
  call d_assemble_lower(a,u,v,d,ubw)
  a0=a
  call lower_to_bt(a,bt_d, ubw, tol,error)
  call cpu_time(t0)
  call convert_bt_to_wb(bt_d, wb_d,error)
  call cpu_time(t1)
  call wb_to_lower(wb_d,a1,error)
  test_name = "Real BT to WB;"
  call d_output_result_lower(test_name,a0,a1,rmax,wb_d%lbw,t0,t1,tol2,error)

  na=1
  ubwa=min(na-1,ubw)
  wb_na_d=d_new_wb(na,lbwmax,ubwmax)
  bt_na_d=d_new_bt(na,lbwmax,ubwmax)
  u=u0; v=v0; d=d0
  call d_assemble_lower(a(1:na,1:na),u(1:na,:),v(:,1:na),d(1:na),ubwa)
  a0(1:na,1:na)=a(1:na,1:na)
  call lower_to_bt(a(1:na,1:na),bt_na_d,ubwa, tol1, error)
  call cpu_time(t0)
  call convert_bt_to_wb(bt_na_d, wb_na_d,error)
  call cpu_time(t1)
  call wb_to_lower(wb_na_d,a1(1:na,1:na),error)
  test_name = "Real BT to WB (n=1);"
  call d_output_result_lower(test_name,a0(1:na,1:na),a1(1:na,1:na),0,wb_na_d%lbw,t0,t1,tol2,error)
  call d_deallocate_wb(wb_na_d)
  call d_deallocate_bt(bt_na_d)

  na=2
  ubwa=min(na-1,ubw)
  wb_na_d=d_new_wb(na,lbwmax,ubwmax)
  bt_na_d=d_new_bt(na,lbwmax,ubwmax)
  u=u0; v=v0; d=d0
  call d_assemble_lower(a(1:na,1:na),u(1:na,:),v(:,1:na),d(1:na),ubwa)
  a0(1:na,1:na)=a(1:na,1:na)
  call lower_to_bt(a(1:na,1:na),bt_na_d,ubwa, tol1, error)
  call cpu_time(t0)
  call convert_bt_to_wb(bt_na_d, wb_na_d,error)
  call cpu_time(t1)
  call wb_to_lower(wb_na_d,a1(1:na,1:na),error)
  test_name = "Real BT to WB (n=2);"
  call d_output_result_lower(test_name,a0(1:na,1:na),a1(1:na,1:na),1,wb_na_d%lbw,t0,t1,tol2,error)
  call d_deallocate_wb(wb_na_d)
  call d_deallocate_bt(bt_na_d)

  na=3
  ubwa=min(na-1,ubw)
  wb_na_d=d_new_wb(na,lbwmax,ubwmax)
  bt_na_d=d_new_bt(na,lbwmax,ubwmax)
  u=u0; v=v0; d=d0
  call d_assemble_lower(a(1:na,1:na),u(1:na,:),v(:,1:na),d(1:na),ubwa)
  a0(1:na,1:na)=a(1:na,1:na)
  call lower_to_bt(a(1:na,1:na),bt_na_d,ubwa, tol1, error)
  call cpu_time(t0)
  call convert_bt_to_wb(bt_na_d, wb_na_d,error)
  call cpu_time(t1)
  call wb_to_lower(wb_na_d,a1(1:na,1:na),error)
  test_name = "Real BT to WB (n=3);"
  call d_output_result_lower(test_name,a0(1:na,1:na),a1(1:na,1:na),1,wb_na_d%lbw,t0,t1,tol2,error)
  call d_deallocate_wb(wb_na_d)
  call d_deallocate_bt(bt_na_d)

  na=4
  ubwa=min(na-1,ubw)
  wb_na_d=d_new_wb(na,lbwmax,ubwmax)
  bt_na_d=d_new_bt(na,lbwmax,ubwmax)
  u=u0; v=v0; d=d0
  call d_assemble_lower(a(1:na,1:na),u(1:na,:),v(:,1:na),d(1:na),ubwa)
  a0(1:na,1:na)=a(1:na,1:na)
  call lower_to_bt(a(1:na,1:na),bt_na_d,ubwa, tol1, error)
  call cpu_time(t0)
  call convert_bt_to_wb(bt_na_d, wb_na_d,error)
  call cpu_time(t1)
  call wb_to_lower(wb_na_d,a1(1:na,1:na),error)
  test_name = "Real BT to WB (n=4);"
  call d_output_result_lower(test_name,a0(1:na,1:na),a1(1:na,1:na),2,wb_na_d%lbw,t0,t1,tol2,error)
  call d_deallocate_wb(wb_na_d)
  call d_deallocate_bt(bt_na_d)

  print *
  print *, "--------------------------------"
  print *
  print *, "Complex WB and BT Conversion Tests:"
  print *
  call random_matrix(u_c)
  call random_matrix(v_c)
  call random_matrix(d_c)
  u0_c=u_c; v0_c=v_c; d0_c=d_c

  call c_assemble_lower(a_c,u_c,v_c,d_c,ubw)
  a0_c=a_c
  call lower_to_wb(a_c,wb_c, ubw, tol,error)
  call cpu_time(t0)
  call convert_wb_to_bt(wb_c, bt_c,error)
  call cpu_time(t1)
  call bt_to_lower(bt_c,a1_c,error)
  test_name = "Complex WB to BT;"
  call c_output_result_lower(test_name,a0_c,a1_c,rmax,bt_c%lbw,t0,t1,tol2,error)
  !
  na=1
  ubwa=min(na-1,ubw)
  wb_na_c=c_new_wb(na,lbwmax,ubwmax)
  bt_na_c=c_new_bt(na,lbwmax,ubwmax)
  u_c=u0_c; v_c=v0_c; d_c=d0_c
  call c_assemble_lower(a_c(1:na,1:na),u_c(1:na,:),v_c(:,1:na),d_c(1:na),ubwa)
  a0_c(1:na,1:na)=a_c(1:na,1:na)
  call lower_to_wb(a_c(1:na,1:na),wb_na_c,ubwa, tol1,error)
  call cpu_time(t0)
  call convert_wb_to_bt(wb_na_c, bt_na_c,error)
  call cpu_time(t1)
  call bt_to_lower(bt_na_c,a1_c(1:na,1:na),error)
  test_name = "Complex WB to BT (n=1);"
  call c_output_result_lower(test_name,a0_c(1:na,1:na),a1_c(1:na,1:na),0,bt_na_c%lbw,t0,t1,tol2,error)
  call c_deallocate_wb(wb_na_c)
  call c_deallocate_bt(bt_na_c)

  na=2
  ubwa=min(na-1,ubw)
  wb_na_c=c_new_wb(na,lbwmax,ubwmax)
  bt_na_c=c_new_bt(na,lbwmax,ubwmax)
  u_c=u0_c; v_c=v0_c; d_c=d0_c
  call c_assemble_lower(a_c(1:na,1:na),u_c(1:na,:),v_c(:,1:na),d_c(1:na),ubwa)
  a0_c(1:na,1:na)=a_c(1:na,1:na)
  call lower_to_wb(a_c(1:na,1:na),wb_na_c,ubwa, tol1,error)
  call cpu_time(t0)
  call convert_wb_to_bt(wb_na_c, bt_na_c,error)
  call cpu_time(t1)
  call bt_to_lower(bt_na_c,a1_c(1:na,1:na),error)
  test_name = "Complex WB to BT (n=2);"
  call c_output_result_lower(test_name,a0_c(1:na,1:na),a1_c(1:na,1:na),1,bt_na_c%lbw,t0,t1,tol2,error)
  call c_deallocate_wb(wb_na_c)
  call c_deallocate_bt(bt_na_c)

  na=3
  ubwa=min(na-1,ubw)
  wb_na_c=c_new_wb(na,lbwmax,ubwmax)
  bt_na_c=c_new_bt(na,lbwmax,ubwmax)
  u_c=u0_c; v_c=v0_c; d_c=d0_c
  call c_assemble_lower(a_c(1:na,1:na),u_c(1:na,:),v_c(:,1:na),d_c(1:na),ubwa)
  a0_c(1:na,1:na)=a_c(1:na,1:na)
  call lower_to_wb(a_c(1:na,1:na),wb_na_c,ubwa, tol1,error)
  call cpu_time(t0)
  call convert_wb_to_bt(wb_na_c, bt_na_c,error)
  call cpu_time(t1)
  call bt_to_lower(bt_na_c,a1_c(1:na,1:na),error)
  test_name = "Complex WB to BT (n=3);"
  call c_output_result_lower(test_name,a0_c(1:na,1:na),a1_c(1:na,1:na),1,bt_na_c%lbw,t0,t1,tol2,error)
  call c_deallocate_wb(wb_na_c)
  call c_deallocate_bt(bt_na_c)

  na=4
  ubwa=min(na-1,ubw)
  wb_na_c=c_new_wb(na,lbwmax,ubwmax)
  bt_na_c=c_new_bt(na,lbwmax,ubwmax)
  u_c=u0_c; v_c=v0_c; d_c=d0_c
  call c_assemble_lower(a_c(1:na,1:na),u_c(1:na,:),v_c(:,1:na),d_c(1:na),ubwa)
  a0_c(1:na,1:na)=a_c(1:na,1:na)
  call lower_to_wb(a_c(1:na,1:na),wb_na_c,ubwa, tol1,error)
  call cpu_time(t0)
  call convert_wb_to_bt(wb_na_c, bt_na_c,error)
  call cpu_time(t1)
  call bt_to_lower(bt_na_c,a1_c(1:na,1:na),error)
  test_name = "Complex WB to BT (n=4);"
  call c_output_result_lower(test_name,a0_c(1:na,1:na),a1_c(1:na,1:na),2,bt_na_c%lbw,t0,t1,tol2,error)
  call c_deallocate_wb(wb_na_c)
  call c_deallocate_bt(bt_na_c)

  ! BT to WB

  ! bt to wb
  u_c=u0_c; v_c=v0_c; d_c=d0_c
  call c_assemble_lower(a_c,u_c,v_c,d_c,ubw)
  a0_c=a_c
  call lower_to_bt(a_c,bt_c, ubw, tol, error)
  call cpu_time(t0)
  call convert_bt_to_wb(bt_c, wb_c, error)
  call cpu_time(t1)
  call wb_to_lower(wb_c,a1_c,error)
  test_name = "Complex BT to WB;"
  call c_output_result_lower(test_name,a0_c,a1_c,rmax,wb_c%lbw,t0,t1,tol2,error)

  na=1
  ubwa=min(na-1,ubw)
  wb_na_c=c_new_wb(na,lbwmax,ubwmax)
  bt_na_c=c_new_bt(na,lbwmax,ubwmax)
  u_c=u0_c; v_c=v0_c; d_c=d0_c
  call c_assemble_lower(a_c(1:na,1:na),u_c(1:na,:),v_c(:,1:na),d_c(1:na),ubwa)
  a0_c(1:na,1:na)=a_c(1:na,1:na)
  call lower_to_bt(a_c(1:na,1:na),bt_na_c,ubwa, tol1, error)
  call cpu_time(t0)
  call convert_bt_to_wb(bt_na_c, wb_na_c,error)
  call cpu_time(t1)
  call wb_to_lower(wb_na_c,a1_c(1:na,1:na),error)
  test_name = "Complex BT to WB (n=1);"
  call c_output_result_lower(test_name,a0_c(1:na,1:na),a1_c(1:na,1:na),0,wb_na_c%lbw,t0,t1,tol2,error)
  call c_deallocate_wb(wb_na_c)
  call c_deallocate_bt(bt_na_c)

  na=2
  ubwa=min(na-1,ubw)
  wb_na_c=c_new_wb(na,lbwmax,ubwmax)
  bt_na_c=c_new_bt(na,lbwmax,ubwmax)
  u_c=u0_c; v_c=v0_c; d_c=d0_c
  call c_assemble_lower(a_c(1:na,1:na),u_c(1:na,:),v_c(:,1:na),d_c(1:na),ubwa)
  a0_c(1:na,1:na)=a_c(1:na,1:na)
  call lower_to_bt(a_c(1:na,1:na),bt_na_c,ubwa, tol1, error)
  call cpu_time(t0)
  call convert_bt_to_wb(bt_na_c, wb_na_c,error)
  call cpu_time(t1)
  call wb_to_lower(wb_na_c,a1_c(1:na,1:na),error)
  test_name = "Complex BT to WB (n=2);"
  call c_output_result_lower(test_name,a0_c(1:na,1:na),a1_c(1:na,1:na),1,wb_na_c%lbw,t0,t1,tol2,error)
  call c_deallocate_wb(wb_na_c)
  call c_deallocate_bt(bt_na_c)

  na=3
  ubwa=min(na-1,ubw)
  wb_na_c=c_new_wb(na,lbwmax,ubwmax)
  bt_na_c=c_new_bt(na,lbwmax,ubwmax)
  u_c=u0_c; v_c=v0_c; d_c=d0_c
  call c_assemble_lower(a_c(1:na,1:na),u_c(1:na,:),v_c(:,1:na),d_c(1:na),ubwa)
  a0_c(1:na,1:na)=a_c(1:na,1:na)
  call lower_to_bt(a_c(1:na,1:na),bt_na_c,ubwa, tol1, error)
  call cpu_time(t0)
  call convert_bt_to_wb(bt_na_c, wb_na_c,error)
  call cpu_time(t1)
  call wb_to_lower(wb_na_c,a1_c(1:na,1:na),error)
  test_name = "Complex BT to WB (n=3);"
  call c_output_result_lower(test_name,a0_c(1:na,1:na),a1_c(1:na,1:na),1,wb_na_c%lbw,t0,t1,tol2,error)
  call c_deallocate_wb(wb_na_c)
  call c_deallocate_bt(bt_na_c)

  na=4
  ubwa=min(na-1,ubw)
  wb_na_c=c_new_wb(na,lbwmax,ubwmax)
  bt_na_c=c_new_bt(na,lbwmax,ubwmax)
  u_c=u0_c; v_c=v0_c; d_c=d0_c
  call c_assemble_lower(a_c(1:na,1:na),u_c(1:na,:),v_c(:,1:na),d_c(1:na),ubwa)
  a0_c(1:na,1:na)=a_c(1:na,1:na)
  call lower_to_bt(a_c(1:na,1:na),bt_na_c,ubwa, tol1, error)
  call cpu_time(t0)
  call convert_bt_to_wb(bt_na_c, wb_na_c,error)
  call cpu_time(t1)
  call wb_to_lower(wb_na_c,a1_c(1:na,1:na),error)
  test_name = "Complex BT to WB (n=4);"
  call c_output_result_lower(test_name,a0_c(1:na,1:na),a1_c(1:na,1:na),2,wb_na_c%lbw,t0,t1,tol2,error)
  call c_deallocate_wb(wb_na_c)
  call c_deallocate_bt(bt_na_c)

end program test_convert_wb_and_bt
