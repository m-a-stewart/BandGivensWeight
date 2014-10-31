program test_qr_factorization
  use mod_nested
  use mod_test_data
  implicit none
  real(kind=dp) :: t0, t1
  integer(kind=int32) :: na, lbwa, ubwa, j, k, numsweeps, lbwmaxa, ubwmaxa
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
  real(kind=dp), dimension(:), allocatable :: cs_d, ss_d
  complex(kind=dp), dimension(:,:), allocatable :: a_c, a0_c, a1_c, q_c
  complex(kind=dp), dimension(:), allocatable :: cs_c, ss_c

  type(d_ub) :: ub_d
  type(c_ub) :: ub_c
  type(d_bv) :: bv_d
  type(c_bv) :: bv_c
  type(d_sweeps) :: sw_d, rsw_d
  type(c_sweeps) :: sw_c, rsw_c

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
  sw_d=d_new_sweeps(na,na,lbwmax)
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
  deallocate(a_d, a0_d, a1_d, q_d)
  call deallocate_ub(ub_d); call deallocate_bv(bv_d)
  call deallocate_sweeps(sw_d)
  
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
  sw_c=c_new_sweeps(na,na,lbwmax)
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
  deallocate(a_c, a0_c, a1_c, q_c)
  call deallocate_ub(ub_c); call deallocate_bv(bv_c)
  call deallocate_sweeps(sw_c)
  ! !
  ! ! full qr factorization with full lbw and ubw.
  ! !
  ! na=10; lbwa=5; ubwa=5; numsweeps=6
  ! u_c=u0_c; v_c=v0_c; d_c=d0_c
  ! allocate(a_c(na,na), a0_c(na,na), a1_c(na,na), q_c(na,na))
  ! ub_c=c_new_ub(na,lbwmax,ubwmax)
  ! bv_c=c_new_bv(na,lbwmax,ubwmax)
  ! sw_c=c_new_sweeps(na,lbwmax)
  ! rsw_c=c_random_sweeps(na,numsweeps)
  ! call c_assemble_upper(a_c,u_c(1:na,1:ubwa),v_c(1:ubwa,1:na),d_c(1:na),lbwa)
  ! call upper_to_ub(a_c,ub_c,lbwa,tol,error)
  ! call sweeps_times_ub(rsw_c,ub_c,bv_c,error)
  ! call bv_to_upper(bv_c,a_c,error)
  ! a0_c=a_c
  ! call cpu_time(t0)
  ! call c_qr_bv_to_ub(bv_c,ub_c,sw_c,error)
  ! call cpu_time(t1)
  ! call ub_to_upper(ub_c,a1_c(1:na,1:na),error)
  ! q_c=reshape([ complex(kind=dp) :: ((c_delta(j,k), j=1,na), k=1,na) ], shape(q_c))
  ! call sweeps_times_general(sw_c,a1_c)
  ! test_name = "Complex QR Factorization (full lbw, ubw);"
  ! call c_output_result_upper(test_name,a0_c,a1_c,na-1,ub_c%ubw,t0,t1,tol2,error)
  ! deallocate(a_c, a0_c, a1_c, q_c)
  ! call deallocate_ub(ub_c); call deallocate_bv(bv_c)
  ! call deallocate_sweeps(sw_c)
  ! call deallocate_sweeps(rsw_c)

  ! na=5
  ! lbwa=2; ubwa=2
  ! numsweeps=lbwa
  ! ubwmaxa=na+1; lbwmaxa=na+1
  ! u_c=u0_c; v_c=v0_c; d_c=d0_c
  ! allocate(a_c(na,na), a0_c(na,na), a1_c(na,na), q_c(na,na))
  ! ub_c=c_new_ub(na,lbwmaxa,ubwmaxa)
  ! bv_c=c_new_bv(na,lbwmaxa,ubwmaxa)
  ! sw_c=c_new_sweeps(na,numsweeps)
  ! call c_assemble_upper(a_c,u_c(1:na,1:ubwa),v_c(1:ubwa,1:na),d_c(1:na),lbwa)
  ! a0_c=a_c
  ! call upper_to_bv(a_c,bv_c,lbwa, tol,error)
  ! call bv_to_upper(bv_c, a1_c,error)
  ! call c_qr_bv_to_ub(bv_c,ub_c,sw_c,error)
  ! call ub_to_upper(ub_c,a1_c, error)
  ! q_c=reshape([ complex(kind=dp) :: ((c_delta(j,k), j=1,na), k=1,na) ], shape(q_c))
  ! call sweeps_times_general(sw_c,q_c)
  ! a_c=matmul(q_c,a1_c)
  ! test_name = "Complex QR Fact. (n=5, lbw=2, ubw=2)"
  ! call c_output_result_upper(test_name,a_c,a0_c,ubwa+lbwa,min(ub_c%ubw,ubwa+numsweeps), &
  !      t0,t1,tol2,error)
  ! deallocate(a_c, a0_c, a1_c, q_c)
  ! call deallocate_ub(ub_c); call deallocate_bv(bv_c)
  ! call deallocate_sweeps(sw_c)

end program test_qr_factorization
