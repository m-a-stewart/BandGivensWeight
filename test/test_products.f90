program test_products
  use mod_orrb
  implicit none

  character(len=*), parameter :: fmt="(A40, 'Time: ',ES8.2,', error: ',ES8.2, ', ', A10)"
  character(len=40) :: test_name
  
  real(kind=dp) :: t0, t1
  integer(kind=int32) :: ma, na, lbwb, ubwb
  type(error_info) :: error
  real(kind=dp), parameter :: tol=1e-15, c=7.0
  !
  real(kind=dp), dimension(:,:), allocatable :: b_d, a_d, c0_d, c1_d
  complex(kind=dp), dimension(:,:), allocatable :: b_z, a_z, c0_z, c1_z

  type(d_ub), allocatable :: ub_d
  type(d_bv), allocatable :: bv_d
  type(d_wb), allocatable :: wb_d
  type(d_bt), allocatable :: bt_d
  type(d_ubt), allocatable :: ubt_d
  type(d_wbv), allocatable :: wbv_d
  
  type(z_ub), allocatable :: ub_z
  type(z_bv), allocatable :: bv_z
  type(z_wb), allocatable :: wb_z
  type(z_bt), allocatable :: bt_z
  type(z_ubt), allocatable :: ubt_z
  type(z_wbv), allocatable :: wbv_z
  
  call initialize_errors
  print *
  print *, "--------------------------------"
  print *
  print *, "Tests of Real Products"
  print *

  ma=40; na=10; lbwb=5; ubwb=7
  a_d=d_random_matrix(ma,na)
  ub_d=d_random_ub(ma,lbwb,ubwb,error=error)
  b_d = general(ub_d,error)
  call cpu_time(t0)
  c0_d=d_product_of_ub_and_general(ub_d,a_d,error)
  call cpu_time(t1)
  c1_d=matmul(b_d,a_d)
  test_name = "Real UB Times Gen. (40x10)"
  call d_output_result(test_name,c0_d,c1_d,t0,t1,c*tol,error)

  ma=40; na=10; lbwb=5; ubwb=7
  a_d=d_random_matrix(ma,na)
  bv_d=d_random_bv(ma,lbwb,ubwb,error=error)
  b_d = general(bv_d,error)
  call cpu_time(t0)
  c0_d=d_product_of_bv_and_general(bv_d,a_d,error)
  call cpu_time(t1)
  c1_d=matmul(b_d,a_d)
  test_name = "Real BV Times Gen. (40x10)"
  call d_output_result(test_name,c0_d,c1_d,t0,t1,c*tol,error)

  ma=40; na=10; lbwb=5; ubwb=7
  a_d=d_random_matrix(ma,na)
  wb_d=d_random_wb(ma,lbwb,ubwb,error=error)
  b_d = general(wb_d,error)
  call cpu_time(t0)
  c0_d=d_product_of_wb_and_general(wb_d,a_d,error)
  call cpu_time(t1)
  c1_d=matmul(b_d,a_d)
  test_name = "Real WB Times Gen. (40x10)"
  call d_output_result(test_name,c0_d,c1_d,t0,t1,c*tol,error)

  ma=40; na=10; lbwb=5; ubwb=7
  a_d=d_random_matrix(ma,na)
  bt_d=d_random_bt(ma,lbwb,ubwb,error=error)
  b_d = general(bt_d,error)
  call cpu_time(t0)
  c0_d=d_product_of_bt_and_general(bt_d,a_d,error)
  call cpu_time(t1)
  c1_d=matmul(b_d,a_d)
  test_name = "Real BT Times Gen. (40x10)"
  call d_output_result(test_name,c0_d,c1_d,t0,t1,c*tol,error)

  ma=40; na=10; lbwb=5; ubwb=7
  a_d=d_random_matrix(ma,na)
  ubt_d=d_random_ubt(ma,lbwb,ubwb,error=error)
  b_d = general(ubt_d,error)
  call cpu_time(t0)
  c0_d=d_product_of_ubt_and_general(ubt_d,a_d,error)
  call cpu_time(t1)
  c1_d=matmul(b_d,a_d)
  test_name = "Real UBT Times Gen. (40x10)"
  call d_output_result(test_name,c0_d,c1_d,t0,t1,c*tol,error)

  ma=40; na=10; lbwb=5; ubwb=7
  a_d=d_random_matrix(ma,na)
  wbv_d=d_random_wbv(ma,lbwb,ubwb,error=error)
  b_d = general(wbv_d,error)
  call cpu_time(t0)
  c0_d=d_product_of_wbv_and_general(wbv_d,a_d,error)
  call cpu_time(t1)
  c1_d=matmul(b_d,a_d)
  test_name = "Real WBV Times Gen. (40x10)"
  call d_output_result(test_name,c0_d,c1_d,t0,t1,c*tol,error)

  ma=10; na=40; lbwb=5; ubwb=7
  a_d=d_random_matrix(ma,na)
  ub_d=d_random_ub(na,lbwb,ubwb,error=error)
  b_d = general(ub_d,error)
  call cpu_time(t0)
  c0_d=d_product_of_general_and_ub(a_d, ub_d,error)
  call cpu_time(t1)
  c1_d=matmul(a_d,b_d)
  test_name = "Real Gen. Times UB. (10x40)"
  call d_output_result(test_name,c0_d,c1_d,t0,t1,c*tol,error)
  
  
  print *
  print *, "--------------------------------"
  print *
  print *, "Tests of Complex Products"
  print *

  ma=40; na=10; lbwb=5; ubwb=7
  a_z=z_random_matrix(ma,na)
  ub_z=z_random_ub(ma,lbwb,ubwb,error=error)
  b_z = general(ub_z,error)
  call cpu_time(t0)
  c0_z=z_product_of_ub_and_general(ub_z,a_z,error)
  call cpu_time(t1)
  c1_z=matmul(b_z,a_z)
  test_name = "Complex UB Times Gen. (40x10)"
  call z_output_result(test_name,c0_z,c1_z,t0,t1,c*tol,error)
  
  ma=40; na=10; lbwb=5; ubwb=7
  a_z=z_random_matrix(ma,na)
  bv_z=z_random_bv(ma,lbwb,ubwb,error=error)
  b_z = general(bv_z,error)
  call cpu_time(t0)
  c0_z=z_product_of_bv_and_general(bv_z,a_z,error)
  call cpu_time(t1)
  c1_z=matmul(b_z,a_z)
  test_name = "Complex BV Times Gen. (40x10)"
  call z_output_result(test_name,c0_z,c1_z,t0,t1,c*tol,error)

  ma=40; na=10; lbwb=5; ubwb=7
  a_z=z_random_matrix(ma,na)
  wb_z=z_random_wb(ma,lbwb,ubwb,error=error)
  b_z = general(wb_z,error)
  call cpu_time(t0)
  c0_z=z_product_of_wb_and_general(wb_z,a_z,error)
  call cpu_time(t1)
  c1_z=matmul(b_z,a_z)
  test_name = "Complex WB Times Gen. (40x10)"
  call z_output_result(test_name,c0_z,c1_z,t0,t1,c*tol,error)
  
  ma=40; na=10; lbwb=5; ubwb=7
  a_z=z_random_matrix(ma,na)
  bt_z=z_random_bt(ma,lbwb,ubwb,error=error)
  b_z = general(bt_z,error)
  call cpu_time(t0)
  c0_z=z_product_of_bt_and_general(bt_z,a_z,error)
  call cpu_time(t1)
  c1_z=matmul(b_z,a_z)
  test_name = "Complex BT Times Gen. (40x10)"
  call z_output_result(test_name,c0_z,c1_z,t0,t1,c*tol,error)

  ma=40; na=10; lbwb=5; ubwb=7
  a_z=z_random_matrix(ma,na)
  ubt_z=z_random_ubt(ma,lbwb,ubwb,error=error)
  b_z = general(ubt_z,error)
  call cpu_time(t0)
  c0_z=z_product_of_ubt_and_general(ubt_z,a_z,error)
  call cpu_time(t1)
  c1_z=matmul(b_z,a_z)
  test_name = "Complex UBT Times Gen. (40x10)"
  call z_output_result(test_name,c0_z,c1_z,t0,t1,c*tol,error)

  ma=40; na=10; lbwb=5; ubwb=7
  a_z=z_random_matrix(ma,na)
  wbv_z=z_random_wbv(ma,lbwb,ubwb,error=error)
  b_z = general(wbv_z,error)
  call cpu_time(t0)
  c0_z=z_product_of_wbv_and_general(wbv_z,a_z,error)
  call cpu_time(t1)
  c1_z=matmul(b_z,a_z)
  test_name = "Complex WBV Times Gen. (40x10)"
  call z_output_result(test_name,c0_z,c1_z,t0,t1,c*tol,error)
  
  ma=10; na=40; lbwb=5; ubwb=7
  a_z=z_random_matrix(ma,na)
  ub_z=z_random_ub(na,lbwb,ubwb,error=error)
  b_z = general(ub_z,error)
  call cpu_time(t0)
  c0_z=z_product_of_general_and_ub(a_z, ub_z,error)
  call cpu_time(t1)
  c1_z=matmul(a_z,b_z)
  test_name = "Complex Gen. Times UB. (10x40)"
  call z_output_result(test_name,c0_z,c1_z,t0,t1,c*tol,error)

contains

  subroutine d_output_result(name,a0,a1,t0,t1,bnd,error)
    character(len=*) :: name
    real(kind=dp), dimension(:,:) :: a0, a1     
    real(kind=dp) :: bnd, t0, t1
    type(error_info) :: error

    real(kind=dp) :: err
    character(len=10) :: test_result

    if (error%code > 0) then
       print *, "Calling error in test: ", name
    else
       err = maxabs(a1-a0)
       if (err < bnd) then
          test_result="PASSED"
       else
          test_result="    FAILED"
       end if
       write (*,fmt) name, t1-t0, err, test_result
    end if
  end subroutine d_output_result

  subroutine z_output_result(name,a0,a1,t0,t1,bnd,error)
    character(len=*) :: name
    complex(kind=dp), dimension(:,:) :: a0, a1     
    real(kind=dp) :: bnd, t0, t1
    type(error_info) :: error

    real(kind=dp) :: err
    character(len=10) :: test_result

    if (error%code > 0) then
       print *, "Calling error in test: ", name
    else
       err = maxabs(a1-a0)
       if (err < bnd) then
          test_result="PASSED"
       else
          test_result="    FAILED"
       end if
       write (*,fmt) name, t1-t0, err, test_result
    end if
  end subroutine z_output_result
  
end program test_products
