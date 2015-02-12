program test_sweeps1
  use mod_orth_rank
  use mod_test_data
  implicit none
  real(kind=dp) :: t0, t1
  integer(kind=int32) :: na, lbwa, ubwa, lbwmaxa, ubwmaxa, numsweeps1
  integer(kind=int32), parameter :: nmax=1000, rmaxa=50
  type(error_info) :: error
  real(kind=dp), parameter :: tol=1e-14, tol1=1e-14, tol2=1e-10
  !
  real(kind=dp), dimension(nmax,rmaxa) :: u_d, u0_d
  real(kind=dp), dimension(rmaxa,nmax) :: v_d, v0_d
  real(kind=dp), dimension(nmax) :: d_d, d0_d
  complex(kind=dp), dimension(nmax,rmaxa) :: u_c, u0_c
  complex(kind=dp), dimension(rmaxa,nmax) :: v_c, v0_c
  complex(kind=dp), dimension(nmax) :: d_c, d0_c

  real(kind=dp), dimension(:,:), allocatable :: a_d, a0_d, a1_d, q_d
  complex(kind=dp), dimension(:,:), allocatable :: a_c, a0_c, a1_c, q_c

  type(d_ub) :: ub_d
  type(c_ub) :: ub_c
  type(d_bv) :: bv_d
  type(c_bv) :: bv_c
  type(d_sweeps1) :: sw_d
  type(c_sweeps1) :: sw_c

  call random_seed
  call random_matrix(u_d)
  call random_matrix(v_d)
  call random_matrix(d_d)
  u0_d=u_d; v0_d=v_d; d0_d=d_d

  call random_matrix(u_c)
  call random_matrix(v_c)
  call random_matrix(d_c)
  u0_c=u_c; v0_c=v_c; d0_c=d_c



  print *
  print *, "--------------------------------"
  print *
  print *, "Real Sweeps1 Tests."
  print *

  na=20; numsweeps1=2
  lbwa=5; ubwa=5
  ubwmaxa=min(na-1,numsweeps1+ubwa+1)
  lbwmaxa=min(na-1,numsweeps1+lbwa)
  u_d=u0_d; v_d=v0_d; d_d=d0_d
  allocate(a_d(na,na), a0_d(na,na), a1_d(na,na), q_d(na,na))
  ub_d=d_new_ub(na,lbwmaxa,ubwmaxa)
  bv_d=d_new_bv(na,lbwmaxa,ubwmaxa)
  sw_d=d_random_sweeps1(na,numsweeps1)
  call d_assemble_upper(a_d,u_d(1:na,1:ubwa),v_d(1:ubwa,1:na),d_d(1:na),lbwa)
  a0_d=a_d
  call upper_to_ub(a_d,ub_d,lbwa, tol,error)
  a1_d=a0_d
  call sweeps1_times_general(sw_d,a1_d)
  call cpu_time(t0)
  call d_sweeps1_times_ub(sw_d,ub_d,bv_d,error)
  call cpu_time(t1)
  call d_bv_to_upper(bv_d,a_d,error)
  test_name = "Real Sweeps1 Times UB"
  call d_output_result_upper(test_name,a_d,a1_d,min(bv_d%ubw,ubwa+numsweeps1),min(bv_d%ubw,ubwa+numsweeps1), &
       t0,t1,tol2,error)
  deallocate(a_d, a0_d, a1_d, q_d)
  call deallocate_ub(ub_d); call deallocate_bv(bv_d)
  call deallocate_sweeps1(sw_d)

  na=20; numsweeps1=20
  lbwa=10; ubwa=10
  ubwmaxa=min(na-1,numsweeps1+ubwa+1)
  lbwmaxa=min(na-1,numsweeps1+lbwa)
  u_d=u0_d; v_d=v0_d; d_d=d0_d
  allocate(a_d(na,na), a0_d(na,na), a1_d(na,na), q_d(na,na))
  ub_d=d_new_ub(na,lbwmaxa,ubwmaxa)
  bv_d=d_new_bv(na,lbwmaxa,ubwmaxa)
  sw_d=d_random_sweeps1(na,numsweeps1)
  call d_assemble_upper(a_d,u_d(1:na,1:ubwa),v_d(1:ubwa,1:na),d_d(1:na),lbwa)
  a0_d=a_d
  call upper_to_ub(a_d,ub_d,lbwa, tol,error)
  a1_d=a0_d
  call sweeps1_times_general(sw_d,a1_d)
  call cpu_time(t0)
  call d_sweeps1_times_ub(sw_d,ub_d,bv_d,error)
  call cpu_time(t1)
  call d_bv_to_upper(bv_d,a_d,error)
  test_name = "Real Sweeps1 Times UB (complete fill)"
  a_d = a_d/maxabs(a_d)
  a1_d = a1_d/maxabs(a1_d)
  call d_output_result_upper(test_name,a_d,a1_d,min(bv_d%ubw,ubwa+numsweeps1),min(bv_d%ubw,ubwa+numsweeps1), &
       t0,t1,tol2,error)
  deallocate(a_d, a0_d, a1_d, q_d)
  call deallocate_ub(ub_d); call deallocate_bv(bv_d)
  call deallocate_sweeps1(sw_d)

  ! Real BV
  na=20; numsweeps1=7
  lbwa=5; ubwa=5
  ubwmaxa=min(na-1,numsweeps1+ubwa+1)
  lbwmaxa=min(na-1,numsweeps1+lbwa)
  u_d=u0_d; v_d=v0_d; d_d=d0_d
  allocate(a_d(na,na), a0_d(na,na), a1_d(na,na), q_d(na,na))
  ub_d=d_new_ub(na,lbwmaxa,ubwmaxa)
  bv_d=d_new_bv(na,lbwmaxa,ubwmaxa)
  sw_d=d_random_sweeps1(na,numsweeps1)
  call d_assemble_upper(a_d,u_d(1:na,1:ubwa),v_d(1:ubwa,1:na),d_d(1:na),lbwa)
  a0_d=a_d
  call upper_to_bv(a_d,bv_d,lbwa, tol,error)
  a1_d=a0_d
  call general_times_sweeps1(a1_d,sw_d)
  call cpu_time(t0)
  call d_bv_times_sweeps1(bv_d, sw_d, ub_d, error)
  call cpu_time(t1)
  call d_ub_to_upper(ub_d,a_d,error)
  test_name = "Real BV Times Sweeps1"
  call d_output_result_upper(test_name,a_d,a1_d,min(ub_d%ubw,ubwa+numsweeps1),min(ub_d%ubw,ubwa+numsweeps1), &
       t0,t1,tol2,error)
  deallocate(a_d, a0_d, a1_d, q_d)
  call deallocate_ub(ub_d); call deallocate_bv(bv_d)
  call deallocate_sweeps1(sw_d)


  na=20; numsweeps1=20
  lbwa=5; ubwa=5
  ubwmaxa=min(na-1,numsweeps1+ubwa+1)
  lbwmaxa=min(na-1,numsweeps1+lbwa)
  u_d=u0_d; v_d=v0_d; d_d=d0_d
  allocate(a_d(na,na), a0_d(na,na), a1_d(na,na), q_d(na,na))
  ub_d=d_new_ub(na,lbwmaxa,ubwmaxa)
  bv_d=d_new_bv(na,lbwmaxa,ubwmaxa)
  sw_d=d_random_sweeps1(na,numsweeps1)
  call d_assemble_upper(a_d,u_d(1:na,1:ubwa),v_d(1:ubwa,1:na),d_d(1:na),lbwa)
  a0_d=a_d
  call upper_to_bv(a_d,bv_d,lbwa, tol,error)
  a1_d=a0_d
  call general_times_sweeps1(a1_d,sw_d)
  call cpu_time(t0)
  call d_bv_times_sweeps1(bv_d, sw_d, ub_d, error)
  call cpu_time(t1)
  call d_ub_to_upper(ub_d,a_d,error)
  test_name = "Real BV Times Sweeps1 (complete fill)"
  call d_output_result_upper(test_name,a_d,a1_d,min(ub_d%ubw,ubwa+numsweeps1),min(ub_d%ubw,ubwa+numsweeps1), &
       t0,t1,tol2,error)
  deallocate(a_d, a0_d, a1_d, q_d)
  call deallocate_ub(ub_d); call deallocate_bv(bv_d)
  call deallocate_sweeps1(sw_d)

  ! Complex

  print *
  print *, "--------------------------------"
  print *
  print *, "Complex Sweeps1 Tests."
  print *

  na=20; numsweeps1=2
  lbwa=5; ubwa=5
  ubwmaxa=min(na-1,numsweeps1+ubwa+1)
  lbwmaxa=min(na-1,numsweeps1+lbwa)
  u_c=u0_c; v_c=v0_c; d_c=d0_c
  allocate(a_c(na,na), a0_c(na,na), a1_c(na,na), q_c(na,na))
  ub_c=c_new_ub(na,lbwmaxa,ubwmaxa)
  bv_c=c_new_bv(na,lbwmaxa,ubwmaxa)
  sw_c=c_random_sweeps1(na,numsweeps1)
  call c_assemble_upper(a_c,u_c(1:na,1:ubwa),v_c(1:ubwa,1:na),d_c(1:na),lbwa)
  a0_c=a_c
  call upper_to_ub(a_c,ub_c,lbwa, tol,error)
  a1_c=a0_c
  call sweeps1_times_general(sw_c,a1_c)
  call cpu_time(t0)
  call c_sweeps1_times_ub(sw_c,ub_c,bv_c,error)
  call cpu_time(t1)
  call c_bv_to_upper(bv_c,a_c,error)
  test_name = "Complex Sweeps1 Times UB"
  call c_output_result_upper(test_name,a_c,a1_c,min(bv_c%ubw,ubwa+numsweeps1),min(bv_c%ubw,ubwa+numsweeps1), &
       t0,t1,tol2,error)
  deallocate(a_c, a0_c, a1_c, q_c)
  call deallocate_ub(ub_c); call deallocate_bv(bv_c)
  call deallocate_sweeps1(sw_c)

  na=20; numsweeps1=20
  lbwa=10; ubwa=10
  ubwmaxa=min(na-1,numsweeps1+ubwa+1)
  lbwmaxa=min(na-1,numsweeps1+lbwa)
  u_c=u0_c; v_c=v0_c; d_c=d0_c
  allocate(a_c(na,na), a0_c(na,na), a1_c(na,na), q_c(na,na))
  ub_c=c_new_ub(na,lbwmaxa,ubwmaxa)
  bv_c=c_new_bv(na,lbwmaxa,ubwmaxa)
  sw_c=c_random_sweeps1(na,numsweeps1)
  call c_assemble_upper(a_c,u_c(1:na,1:ubwa),v_c(1:ubwa,1:na),d_c(1:na),lbwa)
  a0_c=a_c
  call upper_to_ub(a_c,ub_c,lbwa, tol,error)
  a1_c=a0_c
  call sweeps1_times_general(sw_c,a1_c)
  call cpu_time(t0)
  call c_sweeps1_times_ub(sw_c,ub_c,bv_c,error)
  call cpu_time(t1)
  call c_bv_to_upper(bv_c,a_c,error)
  test_name = "Complex Sweeps1 Times UB"
  a_c = a_c/maxabs(a_c)
  a1_c = a1_c/maxabs(a1_c)
  call c_output_result_upper(test_name,a_c,a1_c,min(bv_c%ubw,ubwa+numsweeps1),min(bv_c%ubw,ubwa+numsweeps1), &
       t0,t1,tol2,error)
  deallocate(a_c, a0_c, a1_c, q_c)
  call deallocate_ub(ub_c); call deallocate_bv(bv_c)
  call deallocate_sweeps1(sw_c)

  ! BV

  na=20; numsweeps1=7
  lbwa=5; ubwa=5
  ubwmaxa=min(na-1,numsweeps1+ubwa+1)
  lbwmaxa=min(na-1,numsweeps1+lbwa)
  u_c=u0_c; v_c=v0_c; d_c=d0_c
  allocate(a_c(na,na), a0_c(na,na), a1_c(na,na), q_c(na,na))
  ub_c=c_new_ub(na,lbwmaxa,ubwmaxa)
  bv_c=c_new_bv(na,lbwmaxa,ubwmaxa)
  sw_c=c_random_sweeps1(na,numsweeps1)
  call c_assemble_upper(a_c,u_c(1:na,1:ubwa),v_c(1:ubwa,1:na),d_c(1:na),lbwa)
  a0_c=a_c
  call upper_to_bv(a_c,bv_c,lbwa, tol,error)
  a1_c=a0_c
  call general_times_sweeps1(a1_c,sw_c)
  call cpu_time(t0)
  call c_bv_times_sweeps1(bv_c, sw_c, ub_c, error)
  call cpu_time(t1)
  call c_ub_to_upper(ub_c,a_c,error)
  test_name = "Complex BV Times Sweeps1"
  call c_output_result_upper(test_name,a_c,a1_c,min(ub_c%ubw,ubwa+numsweeps1),min(ub_c%ubw,ubwa+numsweeps1), &
       t0,t1,tol2,error)
  deallocate(a_c, a0_c, a1_c, q_c)
  call deallocate_ub(ub_c); call deallocate_bv(bv_c)
  call deallocate_sweeps1(sw_c)

  na=20; numsweeps1=20
  lbwa=5; ubwa=5
  ubwmaxa=min(na-1,numsweeps1+ubwa+1)
  lbwmaxa=min(na-1,numsweeps1+lbwa)
  u_c=u0_c; v_c=v0_c; d_c=d0_c
  allocate(a_c(na,na), a0_c(na,na), a1_c(na,na), q_c(na,na))
  ub_c=c_new_ub(na,lbwmaxa,ubwmaxa)
  bv_c=c_new_bv(na,lbwmaxa,ubwmaxa)
  sw_c=c_random_sweeps1(na,numsweeps1)
  call c_assemble_upper(a_c,u_c(1:na,1:ubwa),v_c(1:ubwa,1:na),d_c(1:na),lbwa)
  a0_c=a_c
  call upper_to_bv(a_c,bv_c,lbwa, tol,error)
  a1_c=a0_c
  call general_times_sweeps1(a1_c,sw_c)
  call cpu_time(t0)
  call c_bv_times_sweeps1(bv_c, sw_c, ub_c, error)
  call cpu_time(t1)
  call c_ub_to_upper(ub_c,a_c,error)
  test_name = "Complex BV Times Sweeps1 (Complete fill)"
  call c_output_result_upper(test_name,a_c,a1_c,min(ub_c%ubw,ubwa+numsweeps1),min(ub_c%ubw,ubwa+numsweeps1), &
       t0,t1,tol2,error)
  deallocate(a_c, a0_c, a1_c, q_c)
  call deallocate_ub(ub_c); call deallocate_bv(bv_c)
  call deallocate_sweeps1(sw_c)

  !
  ! Transposed variants
  !
  print *
  print *, "--------------------------------"
  print *
  print *, "Real Transposed Sweeps1 Tests."
  print *

  ! Real BV
  na=10;
  lbwa=3; ubwa=3
  numsweeps1=lbwa
  ubwmaxa=min(na-1,numsweeps1+ubwa+1)
  lbwmaxa=min(na-1,numsweeps1+lbwa)
  u_d=u0_d; v_d=v0_d; d_d=d0_d
  allocate(a_d(na,na), a0_d(na,na), a1_d(na,na))
  ub_d=d_new_ub(na,lbwmaxa,ubwmaxa)
  bv_d=d_new_bv(na,lbwmaxa,ubwmaxa)
  sw_d=d_random_sweeps1(na,numsweeps1)
  call d_assemble_upper(a_d,u_d(1:na,1:ubwa),v_d(1:ubwa,1:na),d_d(1:na),lbwa)
  a0_d=a_d
  call upper_to_bv(a_d,bv_d,lbwa, tol,error)
  call d_qr_bv_to_ub(bv_d,ub_d,sw_d,error)
  call ub_to_upper(ub_d,a1_d, error)
  a_d=a0_d
  call upper_to_bv(a_d,bv_d,lbwa, tol,error)
  call cpu_time(t0)
  call trp_sweeps1_times_bv(sw_d, bv_d, ub_d, error)
  call cpu_time(t1)
  call ub_to_upper(ub_d,a_d,error)
  test_name = "Trp Sw, Real BV (na=10,lbwa=3,ubwa=3)"
  call d_output_result_upper(test_name,a_d,a1_d,min(ubwa+lbwa,na-1),ub_d%ubw, &
       t0,t1,tol2,error)
  deallocate(a_d, a0_d, a1_d)
  call deallocate_ub(ub_d); call deallocate_bv(bv_d)
  call deallocate_sweeps1(sw_d)

  na=5
  lbwa=2; ubwa=2
  numsweeps1=lbwa
  ubwmaxa=min(na-1,numsweeps1+ubwa+1)
  lbwmaxa=min(na-1,numsweeps1+lbwa)
  u_d=u0_d; v_d=v0_d; d_d=d0_d
  allocate(a_d(na,na), a0_d(na,na), a1_d(na,na))
  ub_d=d_new_ub(na,lbwmaxa,ubwmaxa)
  bv_d=d_new_bv(na,lbwmaxa,ubwmaxa)
  sw_d=d_new_sweeps1(na,numsweeps1)
  call d_assemble_upper(a_d,u_d(1:na,1:ubwa),v_d(1:ubwa,1:na),d_d(1:na),lbwa)
  a0_d=a_d
  call upper_to_bv(a_d,bv_d,lbwa, tol,error)
  call d_qr_bv_to_ub(bv_d,ub_d,sw_d,error)
  call ub_to_upper(ub_d,a1_d, error)
  a_d=a0_d
  call upper_to_bv(a_d,bv_d,lbwa, tol,error)
  call cpu_time(t0)
  call trp_sweeps1_times_bv(sw_d, bv_d, ub_d, error)
  call cpu_time(t1)
  call ub_to_upper(ub_d,a_d,error)
  test_name = "Trp Sw, Real BV (na=5,lbwa=2,ubwa=2)"
  call d_output_result_upper(test_name,a_d,a1_d,min(ubwa+lbwa,na-1),ub_d%ubw, &
       t0,t1,tol2,error)
  deallocate(a_d, a0_d, a1_d)
  call deallocate_ub(ub_d); call deallocate_bv(bv_d)
  call deallocate_sweeps1(sw_d)

  na=10
  lbwa=9; ubwa=9
  numsweeps1=lbwa
  ubwmaxa=min(na-1,numsweeps1+ubwa+1)
  lbwmaxa=min(na-1,numsweeps1+lbwa)
  u_d=u0_d; v_d=v0_d; d_d=d0_d
  allocate(a_d(na,na), a0_d(na,na), a1_d(na,na))
  ub_d=d_new_ub(na,lbwmaxa,ubwmaxa)
  bv_d=d_new_bv(na,lbwmaxa,ubwmaxa)
  sw_d=d_new_sweeps1(na,numsweeps1)
  call d_assemble_upper(a_d,u_d(1:na,1:ubwa),v_d(1:ubwa,1:na),d_d(1:na),lbwa)
  a0_d=a_d
  call upper_to_bv(a_d,bv_d,lbwa, tol,error)
  call d_qr_bv_to_ub(bv_d,ub_d,sw_d,error)
  call ub_to_upper(ub_d,a1_d, error)
  a_d=a0_d
  call upper_to_bv(a_d,bv_d,lbwa, tol,error)
  call cpu_time(t0)
  call trp_sweeps1_times_bv(sw_d, bv_d, ub_d, error)
  call cpu_time(t1)
  call ub_to_upper(ub_d,a_d,error)
  test_name = "Trp Sw, Real BV (na=10,lbwa=9,ubwa=9)"
  call d_output_result_upper(test_name,a_d,a1_d,min(ubwa+lbwa,na-1),ub_d%ubw, &
       t0,t1,tol2,error)
  deallocate(a_d, a0_d, a1_d)
  call deallocate_ub(ub_d); call deallocate_bv(bv_d)
  call deallocate_sweeps1(sw_d)

  na=100
  lbwa=3; ubwa=3
  numsweeps1=lbwa
  ubwmaxa=min(na-1,numsweeps1+ubwa+1)
  lbwmaxa=min(na-1,numsweeps1+lbwa)
  u_d=u0_d; v_d=v0_d; d_d=d0_d
  allocate(a_d(na,na), a0_d(na,na), a1_d(na,na))
  ub_d=d_new_ub(na,lbwmaxa,ubwmaxa)
  bv_d=d_new_bv(na,lbwmaxa,ubwmaxa)
  sw_d=d_new_sweeps1(na,numsweeps1)
  call d_assemble_upper(a_d,u_d(1:na,1:ubwa),v_d(1:ubwa,1:na),d_d(1:na),lbwa)
  a0_d=a_d
  call upper_to_bv(a_d,bv_d,lbwa, tol,error)
  call d_qr_bv_to_ub(bv_d,ub_d,sw_d,error)
  call ub_to_upper(ub_d,a1_d, error)
  a_d=a0_d
  call upper_to_bv(a_d,bv_d,lbwa, tol,error)
  call cpu_time(t0)
  call trp_sweeps1_times_bv(sw_d, bv_d, ub_d, error)
  call cpu_time(t1)
  call ub_to_upper(ub_d,a_d,error)
  test_name = "Trp Sw, Real BV (na=100,lbwa=3,ubwa=3)"
  call d_output_result_upper(test_name,a_d,a1_d,min(ubwa+lbwa,na-1),ub_d%ubw, &
       t0,t1,tol2,error)
  deallocate(a_d, a0_d, a1_d)
  call deallocate_ub(ub_d); call deallocate_bv(bv_d)
  call deallocate_sweeps1(sw_d)

  !UB
  na=10;
  lbwa=3; ubwa=3
  numsweeps1=2
  ubwmaxa=min(na-1,numsweeps1+ubwa+1)
  lbwmaxa=min(na-1,numsweeps1+lbwa)
  u_d=u0_d; v_d=v0_d; d_d=d0_d
  allocate(a_d(na,na), a0_d(na,na), a1_d(na,na))
  ub_d=d_new_ub(na,lbwmaxa,ubwmaxa)
  bv_d=d_new_bv(na,lbwmaxa,ubwmaxa)
  sw_d=d_random_sweeps1(na,numsweeps1)
  call d_assemble_upper(a_d,u_d(1:na,1:ubwa-numsweeps1),v_d(1:ubwa-numsweeps1,1:na), &
       d_d(1:na),lbwa-numsweeps1)
  a0_d=a_d
  call general_times_sweeps1(a_d, sw_d)
  call upper_to_ub(a_d,ub_d,lbwa, tol,error)
  call d_ub_times_trp_sweeps1(ub_d,sw_d,bv_d,error)
  call bv_to_upper(bv_d,a1_d,error)
  test_name = "Trp Sw, Real UB (na=10,lbwa=3,ubwa=3)"
  call d_output_result_upper(test_name,a0_d,a1_d,min(ubwa+numsweeps1,na-1),bv_d%ubw, &
       t0,t1,tol2,error)
  deallocate(a_d, a0_d, a1_d)
  call deallocate_ub(ub_d); call deallocate_bv(bv_d)
  call deallocate_sweeps1(sw_d)

  na=50;
  lbwa=10; ubwa=10
  numsweeps1=7
  ubwmaxa=min(na-1,numsweeps1+ubwa+1)
  lbwmaxa=min(na-1,numsweeps1+lbwa)
  u_d=u0_d; v_d=v0_d; d_d=d0_d
  allocate(a_d(na,na), a0_d(na,na), a1_d(na,na))
  ub_d=d_new_ub(na,lbwmaxa,ubwmaxa)
  bv_d=d_new_bv(na,lbwmaxa,ubwmaxa)
  sw_d=d_random_sweeps1(na,numsweeps1)
  call d_assemble_upper(a_d,u_d(1:na,1:ubwa-numsweeps1),v_d(1:ubwa-numsweeps1,1:na), &
       d_d(1:na),lbwa-numsweeps1)
  a0_d=a_d
  call general_times_sweeps1(a_d, sw_d)
  call upper_to_ub(a_d,ub_d,lbwa, tol,error)
  call d_ub_times_trp_sweeps1(ub_d,sw_d,bv_d,error)
  call bv_to_upper(bv_d,a1_d,error)
  test_name = "Trp Sw, Real UB (na=50,lbwa=10,ubwa=10)"
  call d_output_result_upper(test_name,a0_d,a1_d,min(ubwa+numsweeps1,na-1),bv_d%ubw, &
       t0,t1,tol2,error)
  deallocate(a_d, a0_d, a1_d)
  call deallocate_ub(ub_d); call deallocate_bv(bv_d)
  call deallocate_sweeps1(sw_d)

  na=20;
  lbwa=19; ubwa=19
  numsweeps1=18
  ubwmaxa=min(na-1,numsweeps1+ubwa+1)
  lbwmaxa=min(na-1,numsweeps1+lbwa)
  u_d=u0_d; v_d=v0_d; d_d=d0_d
  allocate(a_d(na,na), a0_d(na,na), a1_d(na,na))
  ub_d=d_new_ub(na,lbwmaxa,ubwmaxa)
  bv_d=d_new_bv(na,lbwmaxa,ubwmaxa)
  sw_d=d_random_sweeps1(na,numsweeps1)
  call d_assemble_upper(a_d,u_d(1:na,1:ubwa-numsweeps1),v_d(1:ubwa-numsweeps1,1:na), &
       d_d(1:na),lbwa-numsweeps1)
  a0_d=a_d
  call general_times_sweeps1(a_d, sw_d)
  call upper_to_ub(a_d,ub_d,lbwa, tol,error)
  call d_ub_times_trp_sweeps1(ub_d,sw_d,bv_d,error)
  call bv_to_upper(bv_d,a1_d,error)
  test_name = "Trp Sw, Real UB (na=20,lbwa=19,ubwa=19)"
  call d_output_result_upper(test_name,a0_d,a1_d,min(ubwa+numsweeps1,na-1),bv_d%ubw, &
       t0,t1,tol2,error)
  deallocate(a_d, a0_d, a1_d)
  call deallocate_ub(ub_d); call deallocate_bv(bv_d)
  call deallocate_sweeps1(sw_d)

  print *
  print *, "--------------------------------"
  print *
  print *, "Complex Transposed Sweeps1 Tests."
  print *

  na=10;
  lbwa=3; ubwa=3
  numsweeps1=lbwa
  ubwmaxa=min(na-1,numsweeps1+ubwa+1)
  lbwmaxa=min(na-1,numsweeps1+lbwa)
  u_c=u0_c; v_c=v0_c; d_c=d0_c
  allocate(a_c(na,na), a0_c(na,na), a1_c(na,na))
  ub_c=c_new_ub(na,lbwmaxa,ubwmaxa)
  bv_c=c_new_bv(na,lbwmaxa,ubwmaxa)
  sw_c=c_random_sweeps1(na,numsweeps1)
  call c_assemble_upper(a_c,u_c(1:na,1:ubwa),v_c(1:ubwa,1:na),d_c(1:na),lbwa)
  a0_c=a_c
  call upper_to_bv(a_c,bv_c,lbwa, tol,error)
  call c_qr_bv_to_ub(bv_c,ub_c,sw_c,error)
  call ub_to_upper(ub_c,a1_c, error)
  a_c=a0_c
  call upper_to_bv(a_c,bv_c,lbwa, tol,error)
  call cpu_time(t0)
  call trp_sweeps1_times_bv(sw_c, bv_c, ub_c, error)
  call cpu_time(t1)
  call ub_to_upper(ub_c,a_c,error)
  test_name = "Trp Sw, Cpx BV (na=10, lbwa=3, ubwa=3)"
  call c_output_result_upper(test_name,a_c,a1_c,min(ubwa+lbwa,na-1),ub_c%ubw, &
       t0,t1,tol2,error)
  deallocate(a_c, a0_c, a1_c)
  call deallocate_ub(ub_c); call deallocate_bv(bv_c)
  call deallocate_sweeps1(sw_c)

  na=5
  lbwa=2; ubwa=2
  numsweeps1=lbwa
  ubwmaxa=min(na-1,numsweeps1+ubwa+1)
  lbwmaxa=min(na-1,numsweeps1+lbwa)
  u_c=u0_c; v_c=v0_c; d_c=d0_c
  allocate(a_c(na,na), a0_c(na,na), a1_c(na,na))
  ub_c=c_new_ub(na,lbwmaxa,ubwmaxa)
  bv_c=c_new_bv(na,lbwmaxa,ubwmaxa)
  sw_c=c_new_sweeps1(na,numsweeps1)
  call c_assemble_upper(a_c,u_c(1:na,1:ubwa),v_c(1:ubwa,1:na),d_c(1:na),lbwa)
  a0_c=a_c
  call upper_to_bv(a_c,bv_c,lbwa, tol,error)
  call c_qr_bv_to_ub(bv_c,ub_c,sw_c,error)
  call ub_to_upper(ub_c,a1_c, error)
  a_c=a0_c
  call upper_to_bv(a_c,bv_c,lbwa, tol,error)
  call cpu_time(t0)
  call trp_sweeps1_times_bv(sw_c, bv_c, ub_c, error)
  call cpu_time(t1)
  call ub_to_upper(ub_c,a_c,error)
  test_name = "Trp Sw, Cpx BV (na=5, lbwa=2, ubwa=2)"
  call c_output_result_upper(test_name,a_c,a1_c,min(ubwa+lbwa,na-1),ub_c%ubw, &
       t0,t1,tol2,error)
  deallocate(a_c, a0_c, a1_c)
  call deallocate_ub(ub_c); call deallocate_bv(bv_c)
  call deallocate_sweeps1(sw_c)

  na=10
  lbwa=9; ubwa=9
  numsweeps1=lbwa
  ubwmaxa=min(na-1,numsweeps1+ubwa+1)
  lbwmaxa=min(na-1,numsweeps1+lbwa)
  u_c=u0_c; v_c=v0_c; d_c=d0_c
  allocate(a_c(na,na), a0_c(na,na), a1_c(na,na))
  ub_c=c_new_ub(na,lbwmaxa,ubwmaxa)
  bv_c=c_new_bv(na,lbwmaxa,ubwmaxa)
  sw_c=c_new_sweeps1(na,numsweeps1)
  call c_assemble_upper(a_c,u_c(1:na,1:ubwa),v_c(1:ubwa,1:na),d_c(1:na),lbwa)
  a0_c=a_c
  call upper_to_bv(a_c,bv_c,lbwa, tol,error)
  call c_qr_bv_to_ub(bv_c,ub_c,sw_c,error)
  call ub_to_upper(ub_c,a1_c, error)
  a_c=a0_c
  call upper_to_bv(a_c,bv_c,lbwa, tol,error)
  call cpu_time(t0)
  call trp_sweeps1_times_bv(sw_c, bv_c, ub_c, error)
  call cpu_time(t1)
  call ub_to_upper(ub_c,a_c,error)
  test_name = "Trp Sw, Cpx BV (na=10, lbwa=9, ubwa=9)"
  call c_output_result_upper(test_name,a_c,a1_c,min(ubwa+lbwa,na-1),ub_c%ubw, &
       t0,t1,tol2,error)
  deallocate(a_c, a0_c, a1_c)
  call deallocate_ub(ub_c); call deallocate_bv(bv_c)
  call deallocate_sweeps1(sw_c)

  na=100
  lbwa=3; ubwa=3
  numsweeps1=lbwa
  ubwmaxa=min(na-1,numsweeps1+ubwa+1)
  lbwmaxa=min(na-1,numsweeps1+lbwa)
  u_c=u0_c; v_c=v0_c; d_c=d0_c
  allocate(a_c(na,na), a0_c(na,na), a1_c(na,na))
  ub_c=c_new_ub(na,lbwmaxa,ubwmaxa)
  bv_c=c_new_bv(na,lbwmaxa,ubwmaxa)
  sw_c=c_new_sweeps1(na,numsweeps1)
  call c_assemble_upper(a_c,u_c(1:na,1:ubwa),v_c(1:ubwa,1:na),d_c(1:na),lbwa)
  a0_c=a_c
  call upper_to_bv(a_c,bv_c,lbwa, tol,error)
  call c_qr_bv_to_ub(bv_c,ub_c,sw_c,error)
  call ub_to_upper(ub_c,a1_c, error)
  a_c=a0_c
  call upper_to_bv(a_c,bv_c,lbwa, tol,error)
  call cpu_time(t0)
  call trp_sweeps1_times_bv(sw_c, bv_c, ub_c, error)
  call cpu_time(t1)
  call ub_to_upper(ub_c,a_c,error)
  test_name = "Trp Sw, Cpx BV (na=100, lbwa=3, ubwa=3)"
  call c_output_result_upper(test_name,a_c,a1_c,min(ubwa+lbwa,na-1),ub_c%ubw, &
       t0,t1,tol2,error)
  deallocate(a_c, a0_c, a1_c)
  call deallocate_ub(ub_c); call deallocate_bv(bv_c)
  call deallocate_sweeps1(sw_c)

  !UB
  na=10
  lbwa=3; ubwa=3
  numsweeps1=2
  ubwmaxa=min(na-1,numsweeps1+ubwa+1)
  lbwmaxa=min(na-1,numsweeps1+lbwa)
  u_c=u0_c; v_c=v0_c; d_c=d0_c
  allocate(a_c(na,na), a0_c(na,na), a1_c(na,na))
  ub_c=c_new_ub(na,lbwmaxa,ubwmaxa)
  bv_c=c_new_bv(na,lbwmaxa,ubwmaxa)
  sw_c=c_random_sweeps1(na,numsweeps1)
  call c_assemble_upper(a_c,u_c(1:na,1:ubwa-numsweeps1),v_c(1:ubwa-numsweeps1,1:na), &
       d_c(1:na),lbwa-numsweeps1)
  a0_c=a_c
  call general_times_sweeps1(a_c, sw_c)
  call upper_to_ub(a_c,ub_c,lbwa, tol,error)
  call c_ub_times_trp_sweeps1(ub_c,sw_c,bv_c,error)
  call bv_to_upper(bv_c,a1_c,error)
  test_name = "Trp Sw, Cpx UB (na=10,lbwa=3,ubwa=3)"
  call c_output_result_upper(test_name,a0_c,a1_c,min(ubwa+numsweeps1,na-1),bv_c%ubw, &
       t0,t1,tol2,error)
  deallocate(a_c, a0_c, a1_c)
  call deallocate_ub(ub_c); call deallocate_bv(bv_c)
  call deallocate_sweeps1(sw_c)

  na=50;
  lbwa=10; ubwa=10
  numsweeps1=7
  ubwmaxa=min(na-1,numsweeps1+ubwa+1)
  lbwmaxa=min(na-1,numsweeps1+lbwa)
  u_c=u0_c; v_c=v0_c; d_c=d0_c
  allocate(a_c(na,na), a0_c(na,na), a1_c(na,na))
  ub_c=c_new_ub(na,lbwmaxa,ubwmaxa)
  bv_c=c_new_bv(na,lbwmaxa,ubwmaxa)
  sw_c=c_random_sweeps1(na,numsweeps1)
  call c_assemble_upper(a_c,u_c(1:na,1:ubwa-numsweeps1),v_c(1:ubwa-numsweeps1,1:na), &
       d_c(1:na),lbwa-numsweeps1)
  a0_c=a_c
  call general_times_sweeps1(a_c, sw_c)
  call upper_to_ub(a_c,ub_c,lbwa, tol,error)
  call c_ub_times_trp_sweeps1(ub_c,sw_c,bv_c,error)
  call bv_to_upper(bv_c,a1_c,error)
  test_name = "Trp Sw, Cpx UB (na=50,lbwa=10,ubwa=10)"
  call c_output_result_upper(test_name,a0_c,a1_c,min(ubwa+numsweeps1,na-1),bv_c%ubw, &
       t0,t1,tol2,error)
  deallocate(a_c, a0_c, a1_c)
  call deallocate_ub(ub_c); call deallocate_bv(bv_c)
  call deallocate_sweeps1(sw_c)

  na=20;
  lbwa=19; ubwa=19
  numsweeps1=18
  ubwmaxa=min(na-1,numsweeps1+ubwa+1)
  lbwmaxa=min(na-1,numsweeps1+lbwa)
  u_c=u0_c; v_c=v0_c; d_c=d0_c
  allocate(a_c(na,na), a0_c(na,na), a1_c(na,na))
  ub_c=c_new_ub(na,lbwmaxa,ubwmaxa)
  bv_c=c_new_bv(na,lbwmaxa,ubwmaxa)
  sw_c=c_random_sweeps1(na,numsweeps1)
  call c_assemble_upper(a_c,u_c(1:na,1:ubwa-numsweeps1),v_c(1:ubwa-numsweeps1,1:na), &
       d_c(1:na),lbwa-numsweeps1)
  a0_c=a_c
  call general_times_sweeps1(a_c, sw_c)
  call upper_to_ub(a_c,ub_c,lbwa, tol,error)
  call c_ub_times_trp_sweeps1(ub_c,sw_c,bv_c,error)
  call bv_to_upper(bv_c,a1_c,error)
  test_name = "Trp Sw, Cpx UB (na=20,lbwa=19,ubwa=19)"
  call c_output_result_upper(test_name,a0_c,a1_c,min(ubwa+numsweeps1,na-1),bv_c%ubw, &
       t0,t1,tol2,error)
  deallocate(a_c, a0_c, a1_c)
  call deallocate_ub(ub_c); call deallocate_bv(bv_c)
  call deallocate_sweeps1(sw_c)

end program test_sweeps1
