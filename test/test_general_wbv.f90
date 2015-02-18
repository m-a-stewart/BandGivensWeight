program test_general_wbv
  use mod_orth_rank
  use mod_test_data
  implicit none
  real(kind=dp) :: t0, t1
  type(error_info) :: error
  integer(kind=int32), parameter :: n=50, rmaxl=13, lbwmax=rmaxl+5, rmaxu=11, ubwmax=rmaxu+5
  integer(kind=int32) :: na, lbwa, ubwa  
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
  type(d_wbv), allocatable :: wbv_d
  type(c_wbv), allocatable :: wbv_c

  call initialize_errors
  
  wbv_d=d_new_wbv(n,lbwmax,ubwmax)
  wbv_c=c_new_wbv(n,lbwmax,ubwmax)

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
  print *, "Real WBV Decomposition Tests"
  print *
  call d_assemble_general(a,ul,vl,uu,vu,d)
  a0=a
  call cpu_time(t0)
  call general_to_wbv(a,wbv_d,tol,error)
  call cpu_time(t1)
  call wbv_to_general(wbv_d,a1,error)
  test_name="Real WBV;"
  call d_output_result_lower_upper(test_name,a0,a1,rmaxl,wbv_d%lbw,rmaxu,wbv_d%ubw,t0,t0,tol2,error)
  deallocate(wbv_d)
  
  na=40
  lbwa=3; ubwa=5
  wbv_d=d_random_wbv(na,lbwa,ubwa)
  call wbv_to_general(wbv_d,a(1:na,1:na))
  a1(1:na,1:na)=a(1:na,1:na)
  call general_to_wbv(a(1:na,1:na),wbv_d,tol)
  call wbv_to_general(wbv_d,a0(1:na,1:na))
  test_name="Random Real WBV;"  
  call d_output_result_lower_upper(test_name,a0(1:na,1:na),a1(1:na,1:na), &
       lbwa,wbv_d%lbw,ubwa,wbv_d%ubw,t0,t0,tol2,error)
  deallocate(wbv_d)

  print *
  print *, "--------------------------------"
  print *
  print *, "Complex WBV Decomposition Tests"
  print *
  call c_assemble_general(a_c,ul_c,vl_c,uu_c,vu_c,d_c)
  a0_c=a_c
  call cpu_time(t0)
  call general_to_wbv(a_c,wbv_c,tol,error)
  call cpu_time(t1)
  call wbv_to_general(wbv_c,a1_c,error)
  test_name="Complex WBV;"
  call c_output_result_lower_upper(test_name,a0_c,a1_c,rmaxl,wbv_c%lbw,rmaxu, &
       wbv_c%ubw,t0,t0,tol2,error)
  deallocate(wbv_c)
  
  na=40
  lbwa=3; ubwa=5
  wbv_c=c_random_wbv(na,lbwa,ubwa)
  call wbv_to_general(wbv_c,a_c(1:na,1:na))
  a1_c(1:na,1:na)=a_c(1:na,1:na)
  call general_to_wbv(a_c(1:na,1:na),wbv_c,tol)
  call wbv_to_general(wbv_c,a0_c(1:na,1:na))
  test_name="Random Complex WBV;"  
  call c_output_result_lower_upper(test_name,a0_c(1:na,1:na),a1_c(1:na,1:na), &
       lbwa,wbv_c%lbw,ubwa,wbv_c%ubw,t0,t0,tol2,error)
  deallocate(wbv_c)

end program test_general_wbv
