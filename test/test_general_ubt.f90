program test_general_ubt
  use mod_orth_rank
  use mod_test_data
  implicit none
  real(kind=dp) :: t0, t1
  type(error_info) :: error
  integer, parameter :: n=50, rmaxl=13, lbwmax=rmaxl+5, rmaxu=11, ubwmax=rmaxu+5
  real(kind=dp), parameter :: tol=1e-14, tol1=1e-14, tol2=1e-10
  !
  real(kind=dp), dimension(n,n) :: a, a0, a1
  real(kind=dp), dimension(n,rmaxl) :: ul, ul0
  real(kind=dp), dimension(n,rmaxu) :: uu, uu0
  real(kind=dp), dimension(rmaxl,n) :: vl, vl0
  real(kind=dp), dimension(rmaxu,n) :: vu, vu0 
  real(kind=dp), dimension(n) :: d
  complex(kind=dp), dimension(n,n) :: a_c, a0_c, a1_c
  complex(kind=dp), dimension(n,rmaxl) :: ul_c, ul0_c
  complex(kind=dp), dimension(n,rmaxu) :: uu_c, uu0_c
  complex(kind=dp), dimension(rmaxl,n) :: vl_c, vl0_c
  complex(kind=dp), dimension(rmaxu,n) :: vu_c, vu0_c
  complex(kind=dp), dimension(n) :: d_c
  type(d_ubt) :: ubt_d
  type(c_ubt) :: ubt_c
  ubt_d=d_new_ubt(n,lbwmax,ubwmax)
  ubt_c=c_new_ubt(n,lbwmax,ubwmax)

  call random_seed
  call random_matrix(ul)
  call random_matrix(vl)
  call random_matrix(uu)
  call random_matrix(vu)
  call random_matrix(d)
  ul0=ul; vl0=vl
  uu0=uu; vu0=vu
  call random_matrix(ul_c)
  call random_matrix(vl_c)
  call random_matrix(uu_c)
  call random_matrix(vu_c)
  call random_matrix(d_c)
  ul0_c=ul_c; vl0_c=vl_c
  uu0_c=uu_c; vu0_c=vu_c

  ! test one
  print *
  print *, "--------------------------------"
  print *
  print *, "Real UBT Decomposition Tests"
  print *
  call d_assemble_general(a,ul,vl,uu,vu,d)
  a0=a
  call cpu_time(t0)
  call general_to_ubt(a,ubt_d,tol,error)
  call cpu_time(t1)
  call ubt_to_general(ubt_d,a1,error)
  test_name="Real UBT;"
  call d_output_result_lower_upper(test_name,a0,a1,rmaxl,ubt_d%lbw,rmaxu,ubt_d%ubw,t0,t0,tol2,error)


  print *
  print *, "--------------------------------"
  print *
  print *, "Complex UBT Decomposition Tests"
  print *
  call c_assemble_general(a_c,ul_c,vl_c,uu_c,vu_c,d_c)
  a0_c=a_c
  call cpu_time(t0)
  call general_to_ubt(a_c,ubt_c,tol,error)
  call cpu_time(t1)
  call ubt_to_general(ubt_c,a1_c,error)
  test_name="Complex UBT;"
  call c_output_result_lower_upper(test_name,a0_c,a1_c,rmaxl,ubt_c%lbw,rmaxu,ubt_c%ubw,t0,t0,tol2,error)

end program test_general_ubt
