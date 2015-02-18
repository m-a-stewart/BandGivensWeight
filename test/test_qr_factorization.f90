program test_qr_factorization
  use mod_orth_rank
  use mod_test_data
  implicit none
  real(kind=dp) :: t0, t1
  integer(kind=int32) :: na, lbwa, ubwa, j, k
  ! integer(kind=int32) :: lbwmaxa, ubwmaxa, numsweeps
  integer(kind=int32), parameter :: nmax=1000
  type(error_info) :: error
  integer, parameter :: n=50, rmax=13, ubwmax=rmax+5, lbw=2, lbwmax=10
  real(kind=dp), parameter :: tol=1e-14, tol1=1e-14, tol2=1e-10
  !
  real(kind=dp), dimension(nmax,rmax) :: u_d, u0_d
  real(kind=dp), dimension(rmax,nmax) :: v_d, v0_d
  real(kind=dp), dimension(nmax) :: d_d, d0_d
  complex(kind=dp), dimension(nmax,rmax) :: u_c, u0_c
  complex(kind=dp), dimension(rmax,nmax) :: v_c, v0_c
  complex(kind=dp), dimension(nmax) :: d_c, d0_c

  real(kind=dp), dimension(:,:), allocatable :: a_d, a0_d, a1_d, q_d
  complex(kind=dp), dimension(:,:), allocatable :: a_c, a0_c, a1_c, q_c

  type(d_ub), allocatable :: ub_d
  type(c_ub), allocatable :: ub_c
  type(d_bv), allocatable :: bv_d
  type(c_bv), allocatable :: bv_c
  type(d_sweeps), allocatable :: sw_d ! , rsw_d
  type(c_sweeps), allocatable :: sw_c ! , rsw_c

  call initialize_errors

  call random_seed
  call random_matrix(u_d)
  call random_matrix(v_d)
  call random_matrix(d_d)
  u0_d=u_d; v0_d=v_d; d0_d=d_d

  print *
  print *, "--------------------------------"
  print *
  print *, "Real QR Factorization Tests"
  print *
  !
  ! full qr factorization
  !
  na=10; lbwa=3; ubwa=2
  u_d=u0_d; v_d=v0_d; d_d=d0_d
  allocate(a_d(na,na), a0_d(na,na), a1_d(na,na), q_d(na,na))
  ub_d=d_new_ub(na,lbwmax,ubwmax)
  bv_d=d_new_bv(na,lbwmax,ubwmax)
  sw_d=d_new_sweeps(na,lbwa+1, na+lbwa-1,lbwmax)
  call d_assemble_upper(a_d,u_d(1:na,1:ubwa),v_d(1:ubwa,1:na),d_d(1:na),lbwa)
  a0_d=a_d
  call upper_to_bv(a_d,bv_d,lbwa, tol,error)
  call cpu_time(t0)
  call d_qr_bv_to_ub(bv_d,ub_d,sw_d,error)
  call cpu_time(t1)
  call ub_to_upper(ub_d,a1_d,error)
  q_d=reshape([ real(kind=dp) :: ((d_delta(j,k), j=1,na), k=1,na) ], shape(q_d))
  call sweeps_times_general(sw_d,a1_d)
  test_name = "Real QR Factorization"
  call d_output_result_lower_upper(test_name,a0_d,a1_d,0,ub_d%lbw,ubwa+lbwa,ub_d%ubw,t0,t1,tol2,error)
  deallocate(a_d, a0_d, a1_d, q_d, ub_d, bv_d, sw_d)
  
  !
  print *
  print *, "--------------------------------"
  print *
  print *, "Complex QR Factorization Tests"
  print *
  call random_matrix(u_c)
  call random_matrix(v_c)
  call random_matrix(d_c)
  u0_c=u_c; v0_c=v_c; d0_c=d_c
  !
  ! full qr factorization
  !
  na=100; lbwa=7; ubwa=5
  u_c=u0_c; v_c=v0_c; d_c=d0_c
  allocate(a_c(na,na), a0_c(na,na), a1_c(na,na), q_c(na,na))
  ub_c=c_new_ub(na,lbwmax,ubwmax)
  bv_c=c_new_bv(na,lbwmax,ubwmax)
  sw_c=c_new_sweeps(na,lbwa+1,na+lbwa-1,lbwmax)
  call c_assemble_upper(a_c,u_c(1:na,1:ubwa),v_c(1:ubwa,1:na),d_c(1:na),lbwa)
  a0_c=a_c
  call upper_to_bv(a_c,bv_c,lbwa, tol,error)
  call cpu_time(t0)
  call c_qr_bv_to_ub(bv_c,ub_c,sw_c,error)
  call cpu_time(t1)
  call ub_to_upper(ub_c,a1_c,error)
  q_c=reshape([ complex(kind=dp) :: ((c_delta(j,k), j=1,na), k=1,na) ], shape(q_c))
  call sweeps_times_general(sw_c,a1_c)
  test_name = "Complex QR Factorization"
  call c_output_result_lower_upper(test_name,a0_c,a1_c,0,ub_c%lbw,ubwa+lbwa,ub_c%ubw,t0,t1,tol2,error)
  deallocate(a_c, a0_c, a1_c, q_c, ub_c, bv_c, sw_c)

end program test_qr_factorization
