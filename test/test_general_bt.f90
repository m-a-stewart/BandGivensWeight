program test_general_bt
  use mod_orth_rank
  use mod_test_data
  implicit none
  real(kind=dp) :: t0, t1
  type(error_info) :: error
  integer, parameter :: n=50, rmax=13, lbwmax=rmax+5, ubw=2, ubwmax=10
  real(kind=dp), parameter :: tol=1e-14, tol1=1e-14, tol2=1e-10
  !
  real(kind=dp), dimension(n,n) :: a, a0, a1
  real(kind=dp), dimension(n,rmax) :: u, u0
  real(kind=dp), dimension(rmax,n) :: v, v0
  real(kind=dp), dimension(n) :: d
  complex(kind=dp), dimension(n,n) :: a_c, a0_c, a1_c
  complex(kind=dp), dimension(n,rmax) :: u_c, u0_c
  complex(kind=dp), dimension(rmax,n) :: v_c, v0_c
  complex(kind=dp), dimension(n) :: d_c
  type(d_bt), allocatable :: bt_d
  type(c_bt), allocatable :: bt_c
  bt_d=d_new_bt(n,lbwmax,ubwmax)
  bt_c=c_new_bt(n,lbwmax,ubwmax)

  call random_seed
  call random_matrix(u)
  call random_matrix(v)
  call random_matrix(d)
  u0=u; v0=v
  call random_matrix(u_c)
  call random_matrix(v_c)
  call random_matrix(d_c)
  u0_c=u_c; v0_c=v_c

  ! test one
  print *
  print *, "--------------------------------"
  print *
  print *, "Real BT Decomposition Tests"
  print *
  call d_assemble_lower(a,u,v,d,ubw)
  a0=a
  call cpu_time(t0)
  call lower_to_bt(a,bt_d,ubw,tol,error)
  call cpu_time(t1)
  call bt_to_lower(bt_d,a1,error)
  test_name="Real BT;"
  call d_output_result_lower(test_name,a0,a1,rmax,bt_d%lbw,t0,t0,tol2,error)


  print *
  print *, "--------------------------------"
  print *
  print *, "Complex BT Decomposition Tests"
  print *
  call c_assemble_lower(a_c,u_c,v_c,d_c,ubw)
  a0_c=a_c
  call cpu_time(t0)
  call lower_to_bt(a_c,bt_c,ubw,tol,error)
  call cpu_time(t1)
  call bt_to_lower(bt_c,a1_c,error)
  test_name="Complex BT;"
  call c_output_result_lower(test_name,a0_c,a1_c,rmax,bt_c%lbw,t0,t0,tol2,error)

end program test_general_bt
