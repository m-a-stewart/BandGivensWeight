program test_solve
  use nested
  use test_data
  implicit none
  real(kind=dp) :: t0, t1, scale
  integer(kind=int32) :: na, lbwa, ubwa, j, k, nc
  integer(kind=int32), parameter :: nmax=1000
  type(error_info) :: error
  !
  real(kind=dp), dimension(nmax,rmax) :: u_d, u0_d
  real(kind=dp), dimension(rmax,nmax) :: v_d, v0_d
  real(kind=dp), dimension(nmax) :: d_d, d0_d
  complex(kind=dp), dimension(nmax,rmax) :: u_c, u0_c
  complex(kind=dp), dimension(rmax,nmax) :: v_c, v0_c
  complex(kind=dp), dimension(nmax) :: d_c, d0_c

  real(kind=dp), dimension(:,:), allocatable :: a_d, a0_d, a1_d, x_d, rhs_d, rhs0_d
  complex(kind=dp), dimension(:,:), allocatable :: a_c, a0_c, a1_c, x_c, rhs_c, rhs0_c

  real(kind=dp), dimension(:), allocatable :: rhs_v_d, rhs0_v_d, x_v_d
  complex(kind=dp), dimension(:), allocatable :: rhs_v_c, rhs0_v_c, x_v_c

  type(d_ub) :: ub_d
  type(c_ub) :: ub_c
  type(d_bv) :: bv_d
  type(c_bv) :: bv_c
  type(d_sweeps) :: sw_d
  type(c_sweeps) :: sw_c

  call random_seed
  call random_number(u_d)
  call random_number(v_d)
  call random_number(d_d)
  u0_d=u_d; v0_d=v_d; d0_d=d_d

  print *
  print *, "--------------------------------"
  print *
  print *, "Real Back Solver Tests (Timings for back substitution)"
  print *

  na=100; lbwa=5; ubwa=7; nc=3
  u_d=u0_d; v_d=v0_d; d_d=d0_d
  allocate(a_d(na,na), a0_d(na,na), a1_d(na,na),  &
       rhs_d(na,nc), rhs0_d(na,nc), x_d(na,nc))
  call random_number(rhs_d)
  rhs0_d=rhs_d
  ub_d=d_new_ub(na,lbwmax,ubwmax)
  bv_d=d_new_bv(na,lbwmax,ubwmax)
  sw_d=d_new_sweeps(na,lbwmax)
  call d_assemble_a(a_d,u_d(1:na,1:ubwa),v_d(1:ubwa,1:na),d_d(1:na),lbwa)
  a0_d=a_d
  call upper_to_bv(a_d,bv_d,lbwa, tol,error)
  call qr_bv_to_ub(bv_d,ub_d,sw_d,error)
  call trp_sweeps(sw_d)
  call sweeps_times_general(sw_d,rhs_d)
  call cpu_time(t0)
  call back_substitution_ub(ub_d,x_d,rhs_d,error)
  call cpu_time(t1)
  rhs_d=matmul(a0_d,x_d)
  test_name = "Real Back Solver (n=100)"
  scale=maxabs(a0_d)*maxabs(x_d)
  call d_output_result(test_name,rhs0_d/scale,rhs_d/scale,ub_d%ubw,ub_d%ubw,t0,t1,tol2,error)
  deallocate(a_d, a0_d, a1_d, rhs_d, rhs0_d, x_d)
  call deallocate_ub(ub_d); call deallocate_bv(bv_d)
  call deallocate_sweeps(sw_d)

  na=1000; lbwa=3; ubwa=5; nc=3
  u_d=u0_d; v_d=v0_d; d_d=d0_d
  allocate(a_d(na,na), a0_d(na,na), a1_d(na,na), &
       rhs_d(na,nc), rhs0_d(na,nc), x_d(na,nc))
  call random_number(rhs_d)
  rhs0_d=rhs_d
  ub_d=d_new_ub(na,lbwmax,ubwmax)
  bv_d=d_new_bv(na,lbwmax,ubwmax)
  sw_d=d_new_sweeps(na,lbwmax)
  call d_assemble_a(a_d,u_d(1:na,1:ubwa),v_d(1:ubwa,1:na),d_d(1:na),lbwa)
  a0_d=a_d
  call upper_to_bv(a_d,bv_d,lbwa, tol,error)
  call qr_bv_to_ub(bv_d,ub_d,sw_d,error)
  call trp_sweeps(sw_d)
  call sweeps_times_general(sw_d,rhs_d)
  call cpu_time(t0)
  call back_substitution_ub(ub_d,x_d,rhs_d,error)
  call cpu_time(t1)
  rhs_d=matmul(a0_d,x_d)
  test_name = "Real Back Solver (n=1000)"
  scale=maxabs(a0_d)*maxabs(x_d)
  call d_output_result(test_name,rhs0_d/scale,rhs_d/scale,ub_d%ubw,ub_d%ubw,t0,t1,tol2,error)
  deallocate(a_d, a0_d, a1_d, rhs_d, rhs0_d, x_d)
  call deallocate_ub(ub_d); call deallocate_bv(bv_d)
  call deallocate_sweeps(sw_d)

  na=1000; lbwa=3; ubwa=5; nc=3
  u_d=u0_d; v_d=v0_d; d_d=d0_d
  allocate(a_d(na,na), a0_d(na,na), a1_d(na,na), &
       rhs_v_d(na), rhs0_v_d(na), x_v_d(na))
  call random_number(rhs_v_d)
  rhs0_v_d=rhs_v_d
  ub_d=d_new_ub(na,lbwmax,ubwmax)
  bv_d=d_new_bv(na,lbwmax,ubwmax)
  sw_d=d_new_sweeps(na,lbwmax)
  call d_assemble_a(a_d,u_d(1:na,1:ubwa),v_d(1:ubwa,1:na),d_d(1:na),lbwa)
  a0_d=a_d
  call upper_to_bv(a_d,bv_d,lbwa, tol,error)
  call qr_bv_to_ub(bv_d,ub_d,sw_d,error)
  call trp_sweeps(sw_d)
  call sweeps_times_general(sw_d,rhs_v_d)
  call cpu_time(t0)
  call back_substitution_ub(ub_d,x_v_d,rhs_v_d,error)
  call cpu_time(t1)
  rhs_v_d=matmul(a0_d,x_v_d)
  test_name = "Real Back Solver, vector rhs (n=1000)"
  scale=maxabs(a0_d)*maxabs(x_v_d)
  call d_output_result(test_name,reshape(rhs0_v_d/scale,[na,1]),reshape(rhs_v_d/scale,[na,1]), &
       ub_d%ubw,ub_d%ubw,t0,t1,tol2,error)
  deallocate(a_d, a0_d, a1_d, rhs_v_d, rhs0_v_d, x_v_d)
  call deallocate_ub(ub_d); call deallocate_bv(bv_d)
  call deallocate_sweeps(sw_d)

  print *
  print *, "--------------------------------"
  print *
  print *, "Real Forward Solver Tests (Timings for forward substitution)"
  print *

  na=100; lbwa=0; ubwa=7; nc=3
  u_d=u0_d; v_d=v0_d; d_d=d0_d
  allocate(a_d(na,na), a0_d(na,na), a1_d(na,na),  &
       rhs_d(nc,na), rhs0_d(nc,na), x_d(nc,na))
  call random_number(rhs_d)
  rhs0_d=rhs_d
  ub_d=d_new_ub(na,lbwmax,ubwmax)
  bv_d=d_new_bv(na,lbwmax,ubwmax)
  sw_d=d_new_sweeps(na,lbwmax)
  call d_assemble_a(a_d,u_d(1:na,1:ubwa),v_d(1:ubwa,1:na),d_d(1:na),lbwa)
  a0_d=a_d
  call upper_to_bv(a_d,bv_d,lbwa, tol,error)
  call cpu_time(t0)
  call d_forward_substitution_bv(x_d,bv_d,rhs_d,error)
  call cpu_time(t1)
  rhs_d=matmul(x_d,a0_d)
  test_name = "Real Forward Solver (n=100)"
  scale=maxabs(a0_d)*maxabs(x_d)
  call d_output_result(test_name,rhs0_d/scale,rhs_d/scale,bv_d%ubw,bv_d%ubw,t0,t1,tol2,error)
  deallocate(a_d, a0_d, a1_d, rhs_d, rhs0_d, x_d)
  call deallocate_ub(ub_d); call deallocate_bv(bv_d)
  call deallocate_sweeps(sw_d)

  na=1000; lbwa=0; ubwa=7; nc=3
  u_d=u0_d; v_d=v0_d; d_d=d0_d
  allocate(a_d(na,na), a0_d(na,na), a1_d(na,na),  &
       rhs_d(nc,na), rhs0_d(nc,na), x_d(nc,na))
  call random_number(rhs_d)
  rhs0_d=rhs_d
  ub_d=d_new_ub(na,lbwmax,ubwmax)
  bv_d=d_new_bv(na,lbwmax,ubwmax)
  sw_d=d_new_sweeps(na,lbwmax)
  call d_assemble_a(a_d,u_d(1:na,1:ubwa),v_d(1:ubwa,1:na),d_d(1:na),lbwa)
  a0_d=a_d
  call upper_to_bv(a_d,bv_d,lbwa, tol,error)
  call cpu_time(t0)
  call d_forward_substitution_bv(x_d,bv_d,rhs_d,error)
  call cpu_time(t1)
  rhs_d=matmul(x_d,a0_d)
  test_name = "Real Forward Solver (n=1000)"
  scale=maxabs(a0_d)*maxabs(x_d)
  call d_output_result(test_name,rhs0_d/scale,rhs_d/scale,bv_d%ubw,bv_d%ubw,t0,t1,tol2,error)
  deallocate(a_d, a0_d, a1_d, rhs_d, rhs0_d, x_d)
  call deallocate_ub(ub_d); call deallocate_bv(bv_d)
  call deallocate_sweeps(sw_d)

  na=1000; lbwa=0; ubwa=7
  u_d=u0_d; v_d=v0_d; d_d=d0_d
  allocate(a_d(na,na), a0_d(na,na), a1_d(na,na),  &
       rhs_v_d(na), rhs0_v_d(na), x_v_d(na))
  call random_number(rhs_v_d)
  rhs0_v_d=rhs_v_d
  ub_d=d_new_ub(na,lbwmax,ubwmax)
  bv_d=d_new_bv(na,lbwmax,ubwmax)
  sw_d=d_new_sweeps(na,lbwmax)
  call d_assemble_a(a_d,u_d(1:na,1:ubwa),v_d(1:ubwa,1:na),d_d(1:na),lbwa)
  a0_d=a_d
  call upper_to_bv(a_d,bv_d,lbwa, tol,error)
  call cpu_time(t0)
  call forward_substitution_bv(x_v_d,bv_d,rhs_v_d,error)
  call cpu_time(t1)
  rhs_v_d=matmul(x_v_d,a0_d)
  test_name = "Real Forward Solver, Vector (n=1000)"
  scale=maxabs(a0_d)*maxabs(x_v_d)
  call d_output_result(test_name,reshape(rhs0_v_d/scale,[1,na]),&
       reshape(rhs_v_d/scale,[1,na]),bv_d%ubw,bv_d%ubw,t0,t1,tol2,error)
  deallocate(a_d, a0_d, a1_d, rhs_v_d, rhs0_v_d, x_v_d)
  call deallocate_ub(ub_d); call deallocate_bv(bv_d)
  call deallocate_sweeps(sw_d)


  print *
  print *, "--------------------------------"
  print *
  print *, "Complex Solver Tests (Timings for back substitution)"
  print *
  call random_complex(u_c)
  call random_complex(v_c)
  call random_complex(d_c)
  u0_c=u_c; v0_c=v_c; d0_c=d_c

  na=100; lbwa=5; ubwa=7; nc=3
  u_c=u0_c; v_c=v0_c; d_c=d0_c
  allocate(a_c(na,na), a0_c(na,na), a1_c(na,na), &
       rhs_c(na,nc), rhs0_c(na,nc), x_c(na,nc))
  call random_complex(rhs_c)
  rhs0_c=rhs_c
  ub_c=c_new_ub(na,lbwmax,ubwmax)
  bv_c=c_new_bv(na,lbwmax,ubwmax)
  sw_c=c_new_sweeps(na,lbwmax)
  call c_assemble_a(a_c,u_c(1:na,1:ubwa),v_c(1:ubwa,1:na),d_c(1:na),lbwa)
  a0_c=a_c
  call upper_to_bv(a_c,bv_c,lbwa, tol,error)
  call qr_bv_to_ub(bv_c,ub_c,sw_c,error)
  call trp_sweeps(sw_c)
  call sweeps_times_general(sw_c,rhs_c)
  call cpu_time(t0)
  call back_substitution_ub(ub_c,x_c,rhs_c,error)
  call cpu_time(t1)
  rhs_c=matmul(a0_c,x_c)
  test_name = "Complex Back Solver (n=100)"
  scale=maxabs(a0_c)*maxabs(x_c)
  call c_output_result(test_name,rhs0_c/scale,rhs_c/scale,ub_c%ubw,ub_c%ubw,t0,t1,tol2,error)
  deallocate(a_c, a0_c, a1_c, rhs_c, rhs0_c, x_c)
  call deallocate_ub(ub_c); call deallocate_bv(bv_c)
  call deallocate_sweeps(sw_c)

  na=1000; lbwa=5; ubwa=5; nc=3
  u_c=u0_c; v_c=v0_c; d_c=d0_c
  allocate(a_c(na,na), a0_c(na,na), a1_c(na,na), &
       rhs_c(na,nc), rhs0_c(na,nc), x_c(na,nc))
  call random_complex(rhs_c)
  rhs0_c=rhs_c
  ub_c=c_new_ub(na,lbwmax,ubwmax)
  bv_c=c_new_bv(na,lbwmax,ubwmax)
  sw_c=c_new_sweeps(na,lbwmax)
  call c_assemble_a(a_c,u_c(1:na,1:ubwa),v_c(1:ubwa,1:na),d_c(1:na),lbwa)
  a0_c=a_c
  call upper_to_bv(a_c,bv_c,lbwa, tol,error)
  call qr_bv_to_ub(bv_c,ub_c,sw_c,error)
  call cpu_time(t1)
  call trp_sweeps(sw_c)
  call sweeps_times_general(sw_c,rhs_c)
  call cpu_time(t0)
  call back_substitution_ub(ub_c,x_c,rhs_c,error)
  call cpu_time(t1)
  rhs_c=matmul(a0_c,x_c)
  test_name = "Complex Back Solver (n=1000)"
  scale=maxabs(a0_c)*maxabs(x_c)
  call c_output_result(test_name,rhs0_c/scale,rhs_c/scale,ub_c%ubw,ub_c%ubw,t0,t1,tol2,error)
  deallocate(a_c, a0_c, a1_c, rhs_c, rhs0_c, x_c)
  call deallocate_ub(ub_c); call deallocate_bv(bv_c)
  call deallocate_sweeps(sw_c)

  na=1000; lbwa=3; ubwa=5; nc=3
  u_c=u0_c; v_c=v0_c; d_c=d0_c
  allocate(a_c(na,na), a0_c(na,na), a1_c(na,na), &
       rhs_v_c(na), rhs0_v_c(na), x_v_c(na))
  call random_complex(rhs_v_c)
  rhs0_v_c=rhs_v_c
  ub_c=c_new_ub(na,lbwmax,ubwmax)
  bv_c=c_new_bv(na,lbwmax,ubwmax)
  sw_c=c_new_sweeps(na,lbwmax)
  call c_assemble_a(a_c,u_c(1:na,1:ubwa),v_c(1:ubwa,1:na),d_c(1:na),lbwa)
  a0_c=a_c
  call upper_to_bv(a_c,bv_c,lbwa, tol,error)
  call qr_bv_to_ub(bv_c,ub_c,sw_c,error)
  call trp_sweeps(sw_c)
  call sweeps_times_general(sw_c,rhs_v_c)
  call cpu_time(t0)
  call back_substitution_ub(ub_c,x_v_c,rhs_v_c,error)
  call cpu_time(t1)
  rhs_v_c=matmul(a0_c,x_v_c)
  test_name = "Complex Back Solver, vector (n=1000)"
  scale=maxabs(a0_c)*maxabs(x_v_c)
  call c_output_result(test_name,reshape(rhs0_v_c/scale,[na,1]),reshape(rhs_v_c/scale,[na,1]), &
       ub_c%ubw,ub_c%ubw,t0,t1,tol2,error)
  deallocate(a_c, a0_c, a1_c, rhs_v_c, rhs0_v_c, x_v_c)
  call deallocate_ub(ub_c); call deallocate_bv(bv_c)
  call deallocate_sweeps(sw_c)

  print *
  print *, "--------------------------------"
  print *
  print *, "Complex Forward Solver Tests (Timings for forward substitution)"
  print *

  na=100; lbwa=0; ubwa=7; nc=3
  u_c=u0_c; v_c=v0_c; d_c=d0_c
  allocate(a_c(na,na), a0_c(na,na), a1_c(na,na),  &
       rhs_c(nc,na), rhs0_c(nc,na), x_c(nc,na))
  call random_complex(rhs_c)
  rhs0_c=rhs_c
  ub_c=c_new_ub(na,lbwmax,ubwmax)
  bv_c=c_new_bv(na,lbwmax,ubwmax)
  sw_c=c_new_sweeps(na,lbwmax)
  call c_assemble_a(a_c,u_c(1:na,1:ubwa),v_c(1:ubwa,1:na),d_c(1:na),lbwa)
  a0_c=a_c
  call upper_to_bv(a_c,bv_c,lbwa, tol,error)
  call cpu_time(t0)
  call c_forward_substitution_bv(x_c,bv_c,rhs_c,error)
  call cpu_time(t1)
  rhs_c=matmul(x_c,a0_c)
  test_name = "Complex Forward Solver (n=100)"
  scale=maxabs(a0_c)*maxabs(x_c)
  call c_output_result(test_name,rhs0_c/scale,rhs_c/scale,bv_c%ubw,bv_c%ubw,t0,t1,tol2,error)
  deallocate(a_c, a0_c, a1_c, rhs_c, rhs0_c, x_c)
  call deallocate_ub(ub_c); call deallocate_bv(bv_c)
  call deallocate_sweeps(sw_c)

  na=1000; lbwa=0; ubwa=7; nc=3
  u_c=u0_c; v_c=v0_c; d_c=d0_c
  allocate(a_c(na,na), a0_c(na,na), a1_c(na,na),  &
       rhs_c(nc,na), rhs0_c(nc,na), x_c(nc,na))
  call random_complex(rhs_c)
  rhs0_c=rhs_c
  ub_c=c_new_ub(na,lbwmax,ubwmax)
  bv_c=c_new_bv(na,lbwmax,ubwmax)
  sw_c=c_new_sweeps(na,lbwmax)
  call c_assemble_a(a_c,u_c(1:na,1:ubwa),v_c(1:ubwa,1:na),d_c(1:na),lbwa)
  a0_c=a_c
  call upper_to_bv(a_c,bv_c,lbwa, tol,error)
  call cpu_time(t0)
  call c_forward_substitution_bv(x_c,bv_c,rhs_c,error)
  call cpu_time(t1)
  rhs_c=matmul(x_c,a0_c)
  test_name = "Complex Forward Solver (n=1000)"
  scale=maxabs(a0_c)*maxabs(x_c)
  call c_output_result(test_name,rhs0_c/scale,rhs_c/scale,bv_c%ubw,bv_c%ubw,t0,t1,tol2,error)
  deallocate(a_c, a0_c, a1_c, rhs_c, rhs0_c, x_c)
  call deallocate_ub(ub_c); call deallocate_bv(bv_c)
  call deallocate_sweeps(sw_c)

  na=1000; lbwa=0; ubwa=7
  u_c=u0_c; v_c=v0_c; d_c=d0_c
  allocate(a_c(na,na), a0_c(na,na), a1_c(na,na),  &
       rhs_v_c(na), rhs0_v_c(na), x_v_c(na))
  call random_complex(rhs_v_c)
  rhs0_v_c=rhs_v_c
  ub_c=c_new_ub(na,lbwmax,ubwmax)
  bv_c=c_new_bv(na,lbwmax,ubwmax)
  sw_c=c_new_sweeps(na,lbwmax)
  call c_assemble_a(a_c,u_c(1:na,1:ubwa),v_c(1:ubwa,1:na),d_c(1:na),lbwa)
  a0_c=a_c
  call upper_to_bv(a_c,bv_c,lbwa, tol,error)
  call cpu_time(t0)
  call forward_substitution_bv(x_v_c,bv_c,rhs_v_c,error)
  call cpu_time(t1)
  rhs_v_c=matmul(x_v_c,a0_c)
  test_name = "Complex Forward Solver, Vector (n=1000)"
  scale=maxabs(a0_c)*maxabs(x_v_c)
  call c_output_result(test_name,reshape(rhs0_v_c/scale,[1,na]),&
       reshape(rhs_v_c/scale,[1,na]),bv_c%ubw,bv_c%ubw,t0,t1,tol2,error)
  deallocate(a_c, a0_c, a1_c, rhs_v_c, rhs0_v_c, x_v_c)
  call deallocate_ub(ub_c); call deallocate_bv(bv_c)
  call deallocate_sweeps(sw_c)

end program test_solve
