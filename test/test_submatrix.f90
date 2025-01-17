program test_submatrix
  use mod_orrb
  implicit none

  character(len=*), parameter :: fmt="(A45, 'error: ',ES8.2, ', ', A10)"
  character(len=45) :: test_name
  
  integer(kind=int32) :: na, lbw, ubw, la, k0, k1, j0, j1
  type(error_info) :: error
  real(kind=dp), parameter :: tol=1e-15, c=7.0
  !
  real(kind=dp), dimension(:,:), allocatable :: a_d, a0_d
  complex(kind=dp), dimension(:,:), allocatable :: a_z, a0_z

  type(d_ub), allocatable :: ub_d, ub0_d
  type(d_bv), allocatable :: bv_d, bv0_d
  type(d_bt), allocatable :: bt_d, bt0_d
  type(d_wb), allocatable :: wb_d, wb0_d
  type(d_ubt), allocatable :: ubt_d, ubt0_d
  type(d_wbv), allocatable :: wbv_d, wbv0_d
  
  type(z_ub), allocatable :: ub_z, ub0_z
  type(z_bv), allocatable :: bv_z, bv0_z
  type(z_bt), allocatable :: bt_z, bt0_z
  type(z_wb), allocatable :: wb_z, wb0_z
  type(z_ubt), allocatable :: ubt_z, ubt0_z
  type(z_wbv), allocatable :: wbv_z, wbv0_z
  
  call initialize_errors
  print *
  print *, "--------------------------------"
  print *
  print *, "Tests of Real Structured Leading and Trailing Submatrices"
  print *

  na=50; la=10; lbw=3; ubw=4
  ub_d=d_random_ub(na,lbw,ubw,error=error)
  a_d = general(ub_d,error)
  ub0_d  = leading(ub_d,la,error=error)
  a0_d=general(ub0_d,error)
  test_name = "Real UB Submatrix na=50, la=10:"
  call d_output_result(test_name,a_d(1:la,1:la),a0_d,c*tol,error)

  na=50; la=50; lbw=3; ubw=4
  ub_d=d_random_ub(na,lbw,ubw,error=error)
  a_d = general(ub_d,error)
  ub0_d  = leading(ub_d,la,error=error)
  a0_d=general(ub0_d,error)
  test_name = "Real UB Submatrix na=50, la=50:"
  call d_output_result(test_name,a_d(1:la,1:la),a0_d,c*tol,error)

  na=50; la=1; lbw=3; ubw=4
  ub_d=d_random_ub(na,lbw,ubw,error=error)
  a_d = general(ub_d,error)
  ub0_d  = leading(ub_d,la,error=error)
  a0_d=general(ub0_d,error)
  test_name = "Real UB Submatrix na=50, la=1:"
  call d_output_result(test_name,a_d(1:la,1:la),a0_d,c*tol,error)

  na=50; la=2; lbw=3; ubw=4
  ub_d=d_random_ub(na,lbw,ubw,error=error)
  a_d = general(ub_d,error)
  ub0_d  = leading(ub_d,la,error=error)
  a0_d=general(ub0_d,error)
  test_name = "Real UB Submatrix na=50, la=2:"
  call d_output_result(test_name,a_d(1:la,1:la),a0_d,c*tol,error)
  
  na=50; la=10; lbw=3; ubw=4
  bv_d=d_random_bv(na,lbw,ubw,error=error)
  a_d = general(bv_d,error)
  bv0_d  = trailing(bv_d,la,error=error)
  a0_d=general(bv0_d,error)
  test_name = "Real BV Submatrix na=50, la=10:"
  call d_output_result(test_name,a_d(na-la+1:na,na-la+1:na),a0_d,c*tol,error)

  na=50; la=50; lbw=3; ubw=4
  bv_d=d_random_bv(na,lbw,ubw,error=error)
  a_d = general(bv_d,error)
  bv0_d  = trailing(bv_d,la,error=error)
  a0_d=general(bv0_d,error)
  test_name = "Real BV Submatrix na=50, la=50:"
  call d_output_result(test_name,a_d(na-la+1:na,na-la+1:na),a0_d,c*tol,error)

  na=50; la=1; lbw=3; ubw=4
  bv_d=d_random_bv(na,lbw,ubw,error=error)
  a_d = general(bv_d,error)
  bv0_d  = trailing(bv_d,la,error=error)
  a0_d=general(bv0_d,error)
  test_name = "Real BV Submatrix na=50, la=1:"
  call d_output_result(test_name,a_d(na-la+1:na,na-la+1:na),a0_d,c*tol,error)
  
  na=50; la=2; lbw=3; ubw=4
  bv_d=d_random_bv(na,lbw,ubw,error=error)
  a_d = general(bv_d,error)
  bv0_d  = trailing(bv_d,la,error=error)
  a0_d=general(bv0_d,error)
  test_name = "Real BV Submatrix na=50, la=2:"
  call d_output_result(test_name,a_d(na-la+1:na,na-la+1:na),a0_d,c*tol,error)
  
  na=50; la=10; lbw=3; ubw=4
  bt_d=d_random_bt(na,lbw,ubw,error=error)
  a_d = general(bt_d,error)
  bt0_d  = leading(bt_d,la,error=error)
  a0_d=general(bt0_d,error)
  test_name = "Real BT Submatrix na=50, la=10:"
  call d_output_result(test_name,a_d(1:la,1:la),a0_d,c*tol,error)

  na=50; la=50; lbw=3; ubw=4
  bt_d=d_random_bt(na,lbw,ubw,error=error)
  a_d = general(bt_d,error)
  bt0_d  = leading(bt_d,la,error=error)
  a0_d=general(bt0_d,error)
  test_name = "Real BT Submatrix na=50, la=50:"
  call d_output_result(test_name,a_d(1:la,1:la),a0_d,c*tol,error)

  na=50; la=1; lbw=3; ubw=4
  bt_d=d_random_bt(na,lbw,ubw,error=error)
  a_d = general(bt_d,error)
  bt0_d  = leading(bt_d,la,error=error)
  a0_d=general(bt0_d,error)
  test_name = "Real BT Submatrix na=50, la=1:"
  call d_output_result(test_name,a_d(1:la,1:la),a0_d,c*tol,error)

  na=50; la=2; lbw=3; ubw=4
  bt_d=d_random_bt(na,lbw,ubw,error=error)
  a_d = general(bt_d,error)
  bt0_d  = leading(bt_d,la,error=error)
  a0_d=general(bt0_d,error)
  test_name = "Real BT Submatrix na=50, la=2:"
  call d_output_result(test_name,a_d(1:la,1:la),a0_d,c*tol,error)
  
  na=50; la=10; lbw=3; ubw=4
  wb_d=d_random_wb(na,lbw,ubw,error=error)
  a_d = general(wb_d,error)
  wb0_d  = trailing(wb_d,la,error=error)
  a0_d=general(wb0_d,error)
  test_name = "Real WB Submatrix na=50, la=10:"
  call d_output_result(test_name,a_d(na-la+1:na,na-la+1:na),a0_d,c*tol,error)

  na=50; la=50; lbw=3; ubw=4
  wb_d=d_random_wb(na,lbw,ubw,error=error)
  a_d = general(wb_d,error)
  wb0_d  = trailing(wb_d,la,error=error)
  a0_d=general(wb0_d,error)
  test_name = "Real WB Submatrix na=50, la=50:"
  call d_output_result(test_name,a_d(na-la+1:na,na-la+1:na),a0_d,c*tol,error)

  na=50; la=1; lbw=3; ubw=4
  wb_d=d_random_wb(na,lbw,ubw,error=error)
  a_d = general(wb_d,error)
  wb0_d  = trailing(wb_d,la,error=error)
  a0_d=general(wb0_d,error)
  test_name = "Real WB Submatrix na=50, la=1:"
  call d_output_result(test_name,a_d(na-la+1:na,na-la+1:na),a0_d,c*tol,error)
  
  na=50; la=2; lbw=3; ubw=4
  wb_d=d_random_wb(na,lbw,ubw,error=error)
  a_d = general(wb_d,error)
  wb0_d  = trailing(wb_d,la,error=error)
  a0_d=general(wb0_d,error)
  test_name = "Real WB Submatrix na=50, la=2:"
  call d_output_result(test_name,a_d(na-la+1:na,na-la+1:na),a0_d,c*tol,error)
  
  na=50; la=10; lbw=3; ubw=4
  ubt_d=d_random_ubt(na,lbw,ubw,error=error)
  a_d = general(ubt_d,error)
  ubt0_d  = leading(ubt_d,la,error=error)
  a0_d=general(ubt0_d,error)
  test_name = "Real UBT Submatrix na=50, la=10:"
  call d_output_result(test_name,a_d(1:la,1:la),a0_d,c*tol,error)

  na=50; la=50; lbw=3; ubw=4
  ubt_d=d_random_ubt(na,lbw,ubw,error=error)
  a_d = general(ubt_d,error)
  ubt0_d  = leading(ubt_d,la,error=error)
  a0_d=general(ubt0_d,error)
  test_name = "Real UBT Submatrix na=50, la=50:"
  call d_output_result(test_name,a_d(1:la,1:la),a0_d,c*tol,error)

  na=50; la=1; lbw=3; ubw=4
  ubt_d=d_random_ubt(na,lbw,ubw,error=error)
  a_d = general(ubt_d,error)
  ubt0_d  = leading(ubt_d,la,error=error)
  a0_d=general(ubt0_d,error)
  test_name = "Real UBT Submatrix na=50, la=1:"
  call d_output_result(test_name,a_d(1:la,1:la),a0_d,c*tol,error)

  na=50; la=2; lbw=3; ubw=4
  ubt_d=d_random_ubt(na,lbw,ubw,error=error)
  a_d = general(ubt_d,error)
  ubt0_d  = leading(ubt_d,la,error=error)
  a0_d=general(ubt0_d,error)
  test_name = "Real UBT Submatrix na=50, la=2:"
  call d_output_result(test_name,a_d(1:la,1:la),a0_d,c*tol,error)
  
  na=50; la=10; lbw=3; ubw=4
  wbv_d=d_random_wbv(na,lbw,ubw,error=error)
  a_d = general(wbv_d,error)
  wbv0_d  = trailing(wbv_d,la,error=error)
  a0_d=general(wbv0_d,error)
  test_name = "Real WBV Submatrix na=50, la=10:"
  call d_output_result(test_name,a_d(na-la+1:na,na-la+1:na),a0_d,c*tol,error)

  na=50; la=50; lbw=3; ubw=4
  wbv_d=d_random_wbv(na,lbw,ubw,error=error)
  a_d = general(wbv_d,error)
  wbv0_d  = trailing(wbv_d,la,error=error)
  a0_d=general(wbv0_d,error)
  test_name = "Real WBV Submatrix na=50, la=50:"
  call d_output_result(test_name,a_d(na-la+1:na,na-la+1:na),a0_d,c*tol,error)

  na=50; la=1; lbw=3; ubw=4
  wbv_d=d_random_wbv(na,lbw,ubw,error=error)
  a_d = general(wbv_d,error)
  wbv0_d  = trailing(wbv_d,la,error=error)
  a0_d=general(wbv0_d,error)
  test_name = "Real WBV Submatrix na=50, la=1:"
  call d_output_result(test_name,a_d(na-la+1:na,na-la+1:na),a0_d,c*tol,error)

  na=50; la=2; lbw=3; ubw=4
  wbv_d=d_random_wbv(na,lbw,ubw,error=error)
  a_d = general(wbv_d,error)
  wbv0_d  = trailing(wbv_d,la,error=error)
  a0_d=general(wbv0_d,error)
  test_name = "Real WBV Submatrix na=50, la=2:"
  call d_output_result(test_name,a_d(na-la+1:na,na-la+1:na),a0_d,c*tol,error)
  
  
  print *
  print *, "--------------------------------"
  print *
  print *, "Tests of Complex Structured Leading and Trailing Submatrices"
  print *

  na=50; la=10; lbw=3; ubw=4
  ub_z=z_random_ub(na,lbw,ubw,error=error)
  a_z = general(ub_z,error)
  ub0_z  = leading(ub_z,la,error=error)
  a0_z=general(ub0_z,error)
  test_name = "Complex UB Submatrix na=50, la=10:"
  call z_output_result(test_name,a_z(1:la,1:la),a0_z,c*tol,error)
  
  na=50; la=50; lbw=3; ubw=4
  ub_z=z_random_ub(na,lbw,ubw,error=error)
  a_z = general(ub_z,error)
  ub0_z  = leading(ub_z,la,error=error)
  a0_z=general(ub0_z,error)
  test_name = "Complex UB Submatrix na=50, la=50:"
  call z_output_result(test_name,a_z(1:la,1:la),a0_z,c*tol,error)

  na=50; la=1; lbw=3; ubw=4
  ub_z=z_random_ub(na,lbw,ubw,error=error)
  a_z = general(ub_z,error)
  ub0_z  = leading(ub_z,la,error=error)
  a0_z=general(ub0_z,error)
  test_name = "Complex UB Submatrix na=50, la=1:"
  call z_output_result(test_name,a_z(1:la,1:la),a0_z,c*tol,error)

  na=50; la=2; lbw=3; ubw=4
  ub_z=z_random_ub(na,lbw,ubw,error=error)
  a_z = general(ub_z,error)
  ub0_z  = leading(ub_z,la,error=error)
  a0_z=general(ub0_z,error)
  test_name = "Complex UB Submatrix na=50, la=2:"
  call z_output_result(test_name,a_z(1:la,1:la),a0_z,c*tol,error)

  na=50; la=10; lbw=3; ubw=4
  bv_z=z_random_bv(na,lbw,ubw,error=error)
  a_z = general(bv_z,error)
  bv0_z  = trailing(bv_z,la,error=error)
  a0_z=general(bv0_z,error)
  test_name = "Complex BV Submatrix na=50, la=10:"
  call z_output_result(test_name,a_z(na-la+1:na,na-la+1:na),a0_z,c*tol,error)

  na=50; la=50; lbw=3; ubw=4
  bv_z=z_random_bv(na,lbw,ubw,error=error)
  a_z = general(bv_z,error)
  bv0_z  = trailing(bv_z,la,error=error)
  a0_z=general(bv0_z,error)
  test_name = "Complex BV Submatrix na=50, la=50:"
  call z_output_result(test_name,a_z(na-la+1:na,na-la+1:na),a0_z,c*tol,error)
  
  na=50; la=1; lbw=3; ubw=4
  bv_z=z_random_bv(na,lbw,ubw,error=error)
  a_z = general(bv_z,error)
  bv0_z  = trailing(bv_z,la,error=error)
  a0_z=general(bv0_z,error)
  test_name = "Complex BV Submatrix na=50, la=1:"
  call z_output_result(test_name,a_z(na-la+1:na,na-la+1:na),a0_z,c*tol,error)

  na=50; la=2; lbw=3; ubw=4
  bv_z=z_random_bv(na,lbw,ubw,error=error)
  a_z = general(bv_z,error)
  bv0_z  = trailing(bv_z,la,error=error)
  a0_z=general(bv0_z,error)
  test_name = "Complex BV Submatrix na=50, la=2:"
  call z_output_result(test_name,a_z(na-la+1:na,na-la+1:na),a0_z,c*tol,error)

  na=50; la=10; lbw=3; ubw=4
  bt_z=z_random_bt(na,lbw,ubw,error=error)
  a_z = general(bt_z,error)
  bt0_z  = leading(bt_z,la,error=error)
  a0_z=general(bt0_z,error)
  test_name = "Complex BT Submatrix na=50, la=10:"
  call z_output_result(test_name,a_z(1:la,1:la),a0_z,c*tol,error)

  na=50; la=50; lbw=3; ubw=4
  bt_z=z_random_bt(na,lbw,ubw,error=error)
  a_z = general(bt_z,error)
  bt0_z  = leading(bt_z,la,error=error)
  a0_z=general(bt0_z,error)
  test_name = "Complex BT Submatrix na=50, la=50:"
  call z_output_result(test_name,a_z(1:la,1:la),a0_z,c*tol,error)

  na=50; la=1; lbw=3; ubw=4
  bt_z=z_random_bt(na,lbw,ubw,error=error)
  a_z = general(bt_z,error)
  bt0_z  = leading(bt_z,la,error=error)
  a0_z=general(bt0_z,error)
  test_name = "Complex BT Submatrix na=50, la=1:"
  call z_output_result(test_name,a_z(1:la,1:la),a0_z,c*tol,error)

  na=50; la=2; lbw=3; ubw=4
  bt_z=z_random_bt(na,lbw,ubw,error=error)
  a_z = general(bt_z,error)
  bt0_z  = leading(bt_z,la,error=error)
  a0_z=general(bt0_z,error)
  test_name = "Complex BT Submatrix na=50, la=2:"
  call z_output_result(test_name,a_z(1:la,1:la),a0_z,c*tol,error)

  na=50; la=10; lbw=3; ubw=4
  wb_z=z_random_wb(na,lbw,ubw,error=error)
  a_z = general(wb_z,error)
  wb0_z  = trailing(wb_z,la,error=error)
  a0_z=general(wb0_z,error)
  test_name = "Complex WB Submatrix na=50, la=10:"
  call z_output_result(test_name,a_z(na-la+1:na,na-la+1:na),a0_z,c*tol,error)

  na=50; la=50; lbw=3; ubw=4
  wb_z=z_random_wb(na,lbw,ubw,error=error)
  a_z = general(wb_z,error)
  wb0_z  = trailing(wb_z,la,error=error)
  a0_z=general(wb0_z,error)
  test_name = "Complex WB Submatrix na=50, la=50:"
  call z_output_result(test_name,a_z(na-la+1:na,na-la+1:na),a0_z,c*tol,error)

  na=50; la=1; lbw=3; ubw=4
  wb_z=z_random_wb(na,lbw,ubw,error=error)
  a_z = general(wb_z,error)
  wb0_z  = trailing(wb_z,la,error=error)
  a0_z=general(wb0_z,error)
  test_name = "Complex WB Submatrix na=50, la=1:"
  call z_output_result(test_name,a_z(na-la+1:na,na-la+1:na),a0_z,c*tol,error)
  
  na=50; la=2; lbw=3; ubw=4
  wb_z=z_random_wb(na,lbw,ubw,error=error)
  a_z = general(wb_z,error)
  wb0_z  = trailing(wb_z,la,error=error)
  a0_z=general(wb0_z,error)
  test_name = "Complex WB Submatrix na=50, la=2:"
  call z_output_result(test_name,a_z(na-la+1:na,na-la+1:na),a0_z,c*tol,error)

  na=50; la=10; lbw=3; ubw=4
  ubt_z=z_random_ubt(na,lbw,ubw,error=error)
  a_z = general(ubt_z,error)
  ubt0_z  = leading(ubt_z,la,error=error)
  a0_z=general(ubt0_z,error)
  test_name = "Complex UBT Submatrix na=50, la=10:"
  call z_output_result(test_name,a_z(1:la,1:la),a0_z,c*tol,error)

  na=50; la=50; lbw=3; ubw=4
  ubt_z=z_random_ubt(na,lbw,ubw,error=error)
  a_z = general(ubt_z,error)
  ubt0_z  = leading(ubt_z,la,error=error)
  a0_z=general(ubt0_z,error)
  test_name = "Complex UBT Submatrix na=50, la=50:"
  call z_output_result(test_name,a_z(1:la,1:la),a0_z,c*tol,error)

  na=50; la=1; lbw=3; ubw=4
  ubt_z=z_random_ubt(na,lbw,ubw,error=error)
  a_z = general(ubt_z,error)
  ubt0_z  = leading(ubt_z,la,error=error)
  a0_z=general(ubt0_z,error)
  test_name = "Complex UBT Submatrix na=50, la=1:"
  call z_output_result(test_name,a_z(1:la,1:la),a0_z,c*tol,error)

  na=50; la=2; lbw=3; ubw=4
  ubt_z=z_random_ubt(na,lbw,ubw,error=error)
  a_z = general(ubt_z,error)
  ubt0_z  = leading(ubt_z,la,error=error)
  a0_z=general(ubt0_z,error)
  test_name = "Complex UBT Submatrix na=50, la=2:"
  call z_output_result(test_name,a_z(1:la,1:la),a0_z,c*tol,error)
  
  na=50; la=10; lbw=3; ubw=4
  wbv_z=z_random_wbv(na,lbw,ubw,error=error)
  a_z = general(wbv_z,error)
  wbv0_z  = trailing(wbv_z,la,error=error)
  a0_z=general(wbv0_z,error)
  test_name = "Complex WBV Submatrix na=50, la=10:"
  call z_output_result(test_name,a_z(na-la+1:na,na-la+1:na),a0_z,c*tol,error)

  na=50; la=50; lbw=3; ubw=4
  wbv_z=z_random_wbv(na,lbw,ubw,error=error)
  a_z = general(wbv_z,error)
  wbv0_z  = trailing(wbv_z,la,error=error)
  a0_z=general(wbv0_z,error)
  test_name = "Complex WBV Submatrix na=50, la=50:"
  call z_output_result(test_name,a_z(na-la+1:na,na-la+1:na),a0_z,c*tol,error)

  na=50; la=1; lbw=3; ubw=4
  wbv_z=z_random_wbv(na,lbw,ubw,error=error)
  a_z = general(wbv_z,error)
  wbv0_z  = trailing(wbv_z,la,error=error)
  a0_z=general(wbv0_z,error)
  test_name = "Complex WBV Submatrix na=50, la=1:"
  call z_output_result(test_name,a_z(na-la+1:na,na-la+1:na),a0_z,c*tol,error)

  na=50; la=2; lbw=3; ubw=4
  wbv_z=z_random_wbv(na,lbw,ubw,error=error)
  a_z = general(wbv_z,error)
  wbv0_z  = trailing(wbv_z,la,error=error)
  a0_z=general(wbv0_z,error)
  test_name = "Complex WBV Submatrix na=50, la=2:"
  call z_output_result(test_name,a_z(na-la+1:na,na-la+1:na),a0_z,c*tol,error)

  print *
  print *, "--------------------------------"
  print *
  print *, "Tests of Real Unstructured Columns and Rows"
  print *

  na=50; lbw=3; ubw=4; k0=25; k1=30
  ub_d=d_random_ub(na,lbw,ubw,error=error)
  a_d = general(ub_d,error)
  a0_d=columns_of_ub(k0,k1,ub_d,error)
  test_name = "Real columns of UB na=50, k0=25, k1=30:"
  call d_output_result(test_name,a_d(:,k0:k1),a0_d,c*tol,error)

  na=50; lbw=3; ubw=4; k0=1; k1=50
  ub_d=d_random_ub(na,lbw,ubw,error=error)
  a_d = general(ub_d,error)
  a0_d=columns_of_ub(k0,k1,ub_d,error)
  test_name = "Real columns of UB na=50, k0=1, k1=50:"
  call d_output_result(test_name,a_d(:,k0:k1),a0_d,c*tol,error)

  na=50; lbw=3; ubw=4; k0=5; k1=5
  ub_d=d_random_ub(na,lbw,ubw,error=error)
  a_d = general(ub_d,error)
  a0_d=columns_of_ub(k0,k1,ub_d,error)
  test_name = "Real columns of UB na=50, k0=5, k1=5:"
  call d_output_result(test_name,a_d(:,k0:k1),a0_d,c*tol,error)

  na=50; lbw=3; ubw=4; k0=1; k1=1
  ub_d=d_random_ub(na,lbw,ubw,error=error)
  a_d = general(ub_d,error)
  a0_d=columns_of_ub(k0,k1,ub_d,error)
  test_name = "Real columns of UB na=50, k0=1, k1=1:"
  call d_output_result(test_name,a_d(:,k0:k1),a0_d,c*tol,error)

  na=50; lbw=3; ubw=4; k0=50; k1=50
  ub_d=d_random_ub(na,lbw,ubw,error=error)
  a_d = general(ub_d,error)
  a0_d=columns_of_ub(k0,k1,ub_d,error)
  test_name = "Real columns of UB na=50, k0=50, k1=50:"
  call d_output_result(test_name,a_d(:,k0:k1),a0_d,c*tol,error)

  na=50; lbw=3; ubw=4; j0=25; j1=30
  bv_d=d_random_bv(na,lbw,ubw,error=error)
  a_d = general(bv_d,error)
  a0_d=rows_of_bv(j0,j1,bv_d,error)
  test_name = "Real rows of BV na=50, j0=25, j1=30:"
  call d_output_result(test_name,a_d(j0:j1,:),a0_d,c*tol,error)

  na=50; lbw=3; ubw=4; j0=1; j1=50
  bv_d=d_random_bv(na,lbw,ubw,error=error)
  a_d = general(bv_d,error)
  a0_d=rows_of_bv(j0,j1,bv_d,error)
  test_name = "Real rows of BV na=50, j0=1, j1=50:"
  call d_output_result(test_name,a_d(j0:j1,:),a0_d,c*tol,error)

  na=50; lbw=3; ubw=4; j0=5; j1=5
  bv_d=d_random_bv(na,lbw,ubw,error=error)
  a_d = general(bv_d,error)
  a0_d=rows_of_bv(j0,j1,bv_d,error)
  test_name = "Real rows of BV na=50, j0=5, j1=5:"
  call d_output_result(test_name,a_d(j0:j1,:),a0_d,c*tol,error)

  na=50; lbw=3; ubw=4; j0=1; j1=1
  bv_d=d_random_bv(na,lbw,ubw,error=error)
  a_d = general(bv_d,error)
  a0_d=rows_of_bv(j0,j1,bv_d,error)
  test_name = "Real rows of BV na=50, j0=1, j1=1:"
  call d_output_result(test_name,a_d(j0:j1,:),a0_d,c*tol,error)

  na=50; lbw=3; ubw=4; j0=50; j1=50
  bv_d=d_random_bv(na,lbw,ubw,error=error)
  a_d = general(bv_d,error)
  a0_d=rows_of_bv(j0,j1,bv_d,error)
  test_name = "Real rows of BV na=50, j0=50, j1=50:"
  call d_output_result(test_name,a_d(j0:j1,:),a0_d,c*tol,error)
  
  na=50; lbw=3; ubw=4; k0=25; k1=30
  wb_d=d_random_wb(na,lbw,ubw,error=error)
  a_d = general(wb_d,error)
  a0_d=columns_of_wb(k0,k1,wb_d,error)
  test_name = "Real columns of WB na=50, k0=25, k1=30:"
  call d_output_result(test_name,a_d(:,k0:k1),a0_d,c*tol,error)

  na=50; lbw=3; ubw=4; k0=1; k1=50
  wb_d=d_random_wb(na,lbw,ubw,error=error)
  a_d = general(wb_d,error)
  a0_d=columns_of_wb(k0,k1,wb_d,error)
  test_name = "Real columns of WB na=50, k0=1, k1=50:"
  call d_output_result(test_name,a_d(:,k0:k1),a0_d,c*tol,error)

  na=50; lbw=3; ubw=4; k0=5; k1=5
  wb_d=d_random_wb(na,lbw,ubw,error=error)
  a_d = general(wb_d,error)
  a0_d=columns_of_wb(k0,k1,wb_d,error)
  test_name = "Real columns of WB na=50, k0=5, k1=5:"
  call d_output_result(test_name,a_d(:,k0:k1),a0_d,c*tol,error)
  
  na=50; lbw=3; ubw=4; k0=1; k1=1
  wb_d=d_random_wb(na,lbw,ubw,error=error)
  a_d = general(wb_d,error)
  a0_d=columns_of_wb(k0,k1,wb_d,error)
  test_name = "Real columns of WB na=50, k0=1, k1=1:"
  call d_output_result(test_name,a_d(:,k0:k1),a0_d,c*tol,error)

  na=50; lbw=3; ubw=4; k0=50; k1=50
  wb_d=d_random_wb(na,lbw,ubw,error=error)
  a_d = general(wb_d,error)
  a0_d=columns_of_wb(k0,k1,wb_d,error)
  test_name = "Real columns of WB na=50, k0=50, k1=50:"
  call d_output_result(test_name,a_d(:,k0:k1),a0_d,c*tol,error)

  na=50; lbw=3; ubw=4; j0=25; j1=30
  bt_d=d_random_bt(na,lbw,ubw,error=error)
  a_d = general(bt_d,error)
  a0_d=rows_of_bt(j0,j1,bt_d,error)
  test_name = "Real rows of BT na=50, j0=25, j1=30:"
  call d_output_result(test_name,a_d(j0:j1,:),a0_d,c*tol,error)

  na=50; lbw=3; ubw=4; j0=1; j1=50
  bt_d=d_random_bt(na,lbw,ubw,error=error)
  a_d = general(bt_d,error)
  a0_d=rows_of_bt(j0,j1,bt_d,error)
  test_name = "Real rows of BT na=50, j0=1, j1=50:"
  call d_output_result(test_name,a_d(j0:j1,:),a0_d,c*tol,error)

  na=50; lbw=3; ubw=4; j0=5; j1=5
  bt_d=d_random_bt(na,lbw,ubw,error=error)
  a_d = general(bt_d,error)
  a0_d=rows_of_bt(j0,j1,bt_d,error)
  test_name = "Real rows of BT na=50, j0=5, j1=5:"
  call d_output_result(test_name,a_d(j0:j1,:),a0_d,c*tol,error)
  
  na=50; lbw=3; ubw=4; j0=1; j1=1
  bt_d=d_random_bt(na,lbw,ubw,error=error)
  a_d = general(bt_d,error)
  a0_d=rows_of_bt(j0,j1,bt_d,error)
  test_name = "Real rows of BT na=50, j0=1, j1=1:"
  call d_output_result(test_name,a_d(j0:j1,:),a0_d,c*tol,error)

  na=50; lbw=3; ubw=4; j0=50; j1=50
  bt_d=d_random_bt(na,lbw,ubw,error=error)
  a_d = general(bt_d,error)
  a0_d=rows_of_bt(j0,j1,bt_d,error)
  test_name = "Real rows of BT na=50, j0=50, j1=50:"
  call d_output_result(test_name,a_d(j0:j1,:),a0_d,c*tol,error)
  
  
  print *
  print *, "--------------------------------"
  print *
  print *, "Tests of Complex Unstructured Columns and Rows"
  print *

  na=50; lbw=3; ubw=4; k0=25; k1=30
  ub_z=z_random_ub(na,lbw,ubw,error=error)
  a_z = general(ub_z,error)
  a0_z=columns_of_ub(k0,k1,ub_z,error)
  test_name = "Complex columns of UB na=50, k0=25, k1=30:"
  call z_output_result(test_name,a_z(:,k0:k1),a0_z,c*tol,error)

  na=50; lbw=3; ubw=4; k0=1; k1=50
  ub_z=z_random_ub(na,lbw,ubw,error=error)
  a_z = general(ub_z,error)
  a0_z=columns_of_ub(k0,k1,ub_z,error)
  test_name = "Complex columns of UB na=50, k0=1, k1=50:"
  call z_output_result(test_name,a_z(:,k0:k1),a0_z,c*tol,error)

  na=50; lbw=3; ubw=4; k0=5; k1=5
  ub_z=z_random_ub(na,lbw,ubw,error=error)
  a_z = general(ub_z,error)
  a0_z=columns_of_ub(k0,k1,ub_z,error)
  test_name = "Complex columns of UB na=50, k0=5, k1=5:"
  call z_output_result(test_name,a_z(:,k0:k1),a0_z,c*tol,error)

  na=50; lbw=3; ubw=4; k0=1; k1=1
  ub_z=z_random_ub(na,lbw,ubw,error=error)
  a_z = general(ub_z,error)
  a0_z=columns_of_ub(k0,k1,ub_z,error)
  test_name = "Complex columns of UB na=50, k0=1, k1=1:"
  call z_output_result(test_name,a_z(:,k0:k1),a0_z,c*tol,error)

  na=50; lbw=3; ubw=4; k0=50; k1=50
  ub_z=z_random_ub(na,lbw,ubw,error=error)
  a_z = general(ub_z,error)
  a0_z=columns_of_ub(k0,k1,ub_z,error)
  test_name = "Complex columns of UB na=50, k0=50, k1=50:"
  call z_output_result(test_name,a_z(:,k0:k1),a0_z,c*tol,error)
  
  na=50; lbw=3; ubw=4; j0=25; j1=30
  bv_z=z_random_bv(na,lbw,ubw,error=error)
  a_z = general(bv_z,error)
  a0_z=rows_of_bv(j0,j1,bv_z,error)
  test_name = "Complex rows of BV na=50, j0=25, j1=30:"
  call z_output_result(test_name,a_z(j0:j1,:),a0_z,c*tol,error)

  na=50; lbw=3; ubw=4; j0=1; j1=50
  bv_z=z_random_bv(na,lbw,ubw,error=error)
  a_z = general(bv_z,error)
  a0_z=rows_of_bv(j0,j1,bv_z,error)
  test_name = "Complex rows of BV na=50, j0=1, j1=50:"
  call z_output_result(test_name,a_z(j0:j1,:),a0_z,c*tol,error)

  na=50; lbw=3; ubw=4; j0=5; j1=5
  bv_z=z_random_bv(na,lbw,ubw,error=error)
  a_z = general(bv_z,error)
  a0_z=rows_of_bv(j0,j1,bv_z,error)
  test_name = "Complex rows of BV na=50, j0=5, j1=5:"
  call z_output_result(test_name,a_z(j0:j1,:),a0_z,c*tol,error)

  na=50; lbw=3; ubw=4; j0=1; j1=1
  bv_z=z_random_bv(na,lbw,ubw,error=error)
  a_z = general(bv_z,error)
  a0_z=rows_of_bv(j0,j1,bv_z,error)
  test_name = "Complex rows of BV na=50, j0=1, j1=1:"
  call z_output_result(test_name,a_z(j0:j1,:),a0_z,c*tol,error)

  na=50; lbw=3; ubw=4; j0=50; j1=50
  bv_z=z_random_bv(na,lbw,ubw,error=error)
  a_z = general(bv_z,error)
  a0_z=rows_of_bv(j0,j1,bv_z,error)
  test_name = "Complex rows of BV na=50, j0=50, j1=50:"
  call z_output_result(test_name,a_z(j0:j1,:),a0_z,c*tol,error)

  na=50; lbw=3; ubw=4; k0=25; k1=30
  wb_z=z_random_wb(na,lbw,ubw,error=error)
  a_z = general(wb_z,error)
  a0_z=columns_of_wb(k0,k1,wb_z,error)
  test_name = "Complex columns of WB na=50, k0=25, k1=30:"
  call z_output_result(test_name,a_z(:,k0:k1),a0_z,c*tol,error)

  na=50; lbw=3; ubw=4; k0=1; k1=50
  wb_z=z_random_wb(na,lbw,ubw,error=error)
  a_z = general(wb_z,error)
  a0_z=columns_of_wb(k0,k1,wb_z,error)
  test_name = "Complex columns of WB na=50, k0=1, k1=50:"
  call z_output_result(test_name,a_z(:,k0:k1),a0_z,c*tol,error)

  na=50; lbw=3; ubw=4; k0=5; k1=5
  wb_z=z_random_wb(na,lbw,ubw,error=error)
  a_z = general(wb_z,error)
  a0_z=columns_of_wb(k0,k1,wb_z,error)
  test_name = "Complex columns of WB na=50, k0=5, k1=5:"
  call z_output_result(test_name,a_z(:,k0:k1),a0_z,c*tol,error)
  
  na=50; lbw=3; ubw=4; k0=1; k1=1
  wb_z=z_random_wb(na,lbw,ubw,error=error)
  a_z = general(wb_z,error)
  a0_z=columns_of_wb(k0,k1,wb_z,error)
  test_name = "Complex columns of WB na=50, k0=1, k1=1:"
  call z_output_result(test_name,a_z(:,k0:k1),a0_z,c*tol,error)

  na=50; lbw=3; ubw=4; k0=50; k1=50
  wb_z=z_random_wb(na,lbw,ubw,error=error)
  a_z = general(wb_z,error)
  a0_z=columns_of_wb(k0,k1,wb_z,error)
  test_name = "Complex columns of WB na=50, k0=50, k1=50:"
  call z_output_result(test_name,a_z(:,k0:k1),a0_z,c*tol,error)
  
  na=50; lbw=3; ubw=4; j0=1; j1=50
  bt_z=z_random_bt(na,lbw,ubw,error=error)
  a_z = general(bt_z,error)
  a0_z=rows_of_bt(j0,j1,bt_z,error)
  test_name = "Complex rows of BT na=50, j0=1, j1=50:"
  call z_output_result(test_name,a_z(j0:j1,:),a0_z,c*tol,error)

  na=50; lbw=3; ubw=4; j0=5; j1=5
  bt_z=z_random_bt(na,lbw,ubw,error=error)
  a_z = general(bt_z,error)
  a0_z=rows_of_bt(j0,j1,bt_z,error)
  test_name = "Complex rows of BT na=50, j0=5, j1=5:"
  call z_output_result(test_name,a_z(j0:j1,:),a0_z,c*tol,error)
  
  na=50; lbw=3; ubw=4; j0=1; j1=1
  bt_z=z_random_bt(na,lbw,ubw,error=error)
  a_z = general(bt_z,error)
  a0_z=rows_of_bt(j0,j1,bt_z,error)
  test_name = "Complex rows of BT na=50, j0=1, j1=1:"
  call z_output_result(test_name,a_z(j0:j1,:),a0_z,c*tol,error)

  na=50; lbw=3; ubw=4; j0=50; j1=50
  bt_z=z_random_bt(na,lbw,ubw,error=error)
  a_z = general(bt_z,error)
  a0_z=rows_of_bt(j0,j1,bt_z,error)
  test_name = "Complex rows of BT na=50, j0=50, j1=50:"
  call z_output_result(test_name,a_z(j0:j1,:),a0_z,c*tol,error)

contains

  subroutine d_output_result(name,a0,a1,bnd,error)
    character(len=*) :: name
    real(kind=dp), dimension(:,:) :: a0, a1     
    real(kind=dp) :: bnd
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
       write (*,fmt) name, err, test_result
    end if
  end subroutine d_output_result

  subroutine z_output_result(name,a0,a1,bnd,error)
    character(len=*) :: name
    complex(kind=dp), dimension(:,:) :: a0, a1     
    real(kind=dp) :: bnd
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
       write (*,fmt) name, err, test_result
    end if
  end subroutine z_output_result
  
end program test_submatrix
