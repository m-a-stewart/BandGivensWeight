program test_convert
  use general_ub
  use general_bv
  use utility
  use assemble
  use conversions_ub_to_bv
  use conversions_bv_to_ub
  use band_types
  use nested_types
  use test_data
  implicit none

  real(kind=dp) :: t0, t1
  integer(kind=int32) :: error
  character(len=*), parameter :: fmt="(A40, 'Time: ',ES8.2,', ubw: ',I3,', error: ',ES8.2)"
  !
  real(kind=dp), dimension(n,n) :: a, a0, a1
  real(kind=dp), dimension(n,rmax) :: u, u0
  real(kind=dp), dimension(rmax,n) :: v, v0
  real(kind=dp), dimension(n) :: d
  complex(kind=dp), dimension(n,n) :: a_c, a0_c, a1_c
  complex(kind=dp), dimension(n,rmax) :: u_c, u0_c
  complex(kind=dp), dimension(rmax,n) :: v_c, v0_c
  complex(kind=dp), dimension(n) :: d_c

  type(d_ub) :: ub_d
  type(c_ub) :: ub_c
  type(d_bv) :: bv_d
  type(c_bv) :: bv_c
  ub_d=d_new_ub(n,lbwmax,ubwmax)
  ub_c=c_new_ub(n,lbwmax,ubwmax)
  bv_d=d_new_bv(n,lbwmax,ubwmax)
  bv_c=c_new_bv(n,lbwmax,ubwmax)

  ! real
  call random_seed
  call random_number(u)
  call random_number(v)
  call random_number(d)
  u0=u; v0=v
  print *
  print *, "--------------------------------"
  print *
  print *, "Real UB and BV Conversion Tests:"
  print *
  ! ub to bv
  call d_assemble_a(a,u,v,d,lbw)
  a0=a
  call upper_to_ub(a,ub_d,lbw, tol,error)
  call cpu_time(t0)
  call convert_ub_to_bv(ub_d, bv_d,error)
  call cpu_time(t1)
  call bv_to_upper(bv_d,a1,error)
  test_name = "Real UB to BV;"
  call d_output_result(test_name,a0,a1,rmax,bv_d%ubw,t0,t1,tol2,error)
  ! bv to ub
  a=a0
  call upper_to_bv(a,bv_d,lbw,tol,error)
  call cpu_time(t0)
  call convert_bv_to_ub(bv_d, ub_d, error)
  call cpu_time(t1)
  call ub_to_upper(ub_d,a1,error)
  test_name="Real BV to UB;"
  call d_output_result(test_name,a0,a1,rmax,ub_d%ubw,t0,t1,tol2,error)
  print *
  print *, "--------------------------------"
  print *
  print *, "Real UB and BV Conversion Tests:"
  print *
  call random_complex(u_c)
  call random_complex(v_c)
  call random_complex(d_c)
  u0_c=u_c; v0_c=v_c
  call c_assemble_a(a_c, u_c, v_c, d_c, lbw)
  a0_c=a_c
  call upper_to_ub(a_c, ub_c, lbw, tol, error)
  call cpu_time(t0)
  call convert_ub_to_bv(ub_c, bv_c, error)
  call cpu_time(t1)
  call bv_to_upper(bv_c, a1_c, error)
  test_name="Complex UB to BV;"
  call c_output_result(test_name,a0_c,a1_c,rmax,bv_c%ubw,t0,t1,tol2,error)
  ! bv to ub
  a_c=a0_c
  call upper_to_bv(a_c, bv_c, lbw, tol, error)
  call cpu_time(t0)
  call convert_bv_to_ub(bv_c, ub_c, error)
  call cpu_time(t1)
  call ub_to_upper(ub_c, a1_c, error)
  test_name="Complex BV to UB;"
  call c_output_result(test_name,a0_c,a1_c,rmax,ub_c%ubw,t0,t1,tol2,error)
  print *

end program test_convert
