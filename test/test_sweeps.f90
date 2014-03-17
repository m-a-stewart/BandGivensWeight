program test_sweeps
  use nested
  use test_data
  implicit none
  real(kind=dp) :: t0, t1
  integer(kind=int32) :: na, lbwa, ubwa, j, k, lbwmaxa, ubwmaxa, numsweeps
  integer(kind=int32), parameter :: nmax=1000, rmaxa=50
  type(error_info) :: error
  !
  real(kind=dp), dimension(nmax,rmaxa) :: u_d, u0_d
  real(kind=dp), dimension(rmaxa,nmax) :: v_d, v0_d
  real(kind=dp), dimension(nmax) :: d_d, d0_d
  complex(kind=dp), dimension(nmax,rmaxa) :: u_c, u0_c
  complex(kind=dp), dimension(rmaxa,nmax) :: v_c, v0_c
  complex(kind=dp), dimension(nmax) :: d_c, d0_c

  real(kind=dp), dimension(:,:), allocatable :: a_d, a0_d, a1_d, q_d
  real(kind=dp), dimension(:), allocatable :: cs_d, ss_d
  complex(kind=dp), dimension(:,:), allocatable :: a_c, a0_c, a1_c, q_c
  complex(kind=dp), dimension(:), allocatable :: cs_c, ss_c

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

  call random_complex(u_c)
  call random_complex(v_c)
  call random_complex(d_c)
  u0_c=u_c; v0_c=v_c; d0_c=d_c



  print *
  print *, "--------------------------------"
  print *
  print *, "Real Sweeps Tests."
  print *

  na=20; numsweeps=2
  lbwa=5; ubwa=5
  ubwmaxa=na+1; lbwmaxa=na+1
  u_d=u0_d; v_d=v0_d; d_d=d0_d
  allocate(a_d(na,na), a0_d(na,na), a1_d(na,na), q_d(na,na))
  ub_d=d_new_ub(na,lbwmaxa,ubwmaxa)
  bv_d=d_new_bv(na,lbwmaxa,ubwmaxa)
  sw_d=d_random_sweeps(na,numsweeps)
  call d_assemble_a(a_d,u_d(1:na,1:ubwa),v_d(1:ubwa,1:na),d_d(1:na),lbwa)
  a0_d=a_d
  call upper_to_ub(a_d,ub_d,lbwa, tol,error)
  a1_d=a0_d
  call sweeps_times_general(sw_d,a1_d)
  call cpu_time(t0)
  call d_sweeps_times_ub(sw_d,ub_d,bv_d,error)
  call cpu_time(t1)
  call d_bv_to_upper(bv_d,a_d,error)
  test_name = "Real Sweeps Times UB"
  call d_output_result(test_name,a_d,a1_d,min(bv_d%ubw,ubwa+numsweeps),min(bv_d%ubw,ubwa+numsweeps), &
       t0,t1,tol2,error)
  deallocate(a_d, a0_d, a1_d, q_d)
  call deallocate_ub(ub_d); call deallocate_bv(bv_d)
  call deallocate_sweeps(sw_d)

  na=20; numsweeps=20
  lbwa=10; ubwa=10
  ubwmaxa=na+1; lbwmaxa=na+1
  u_d=u0_d; v_d=v0_d; d_d=d0_d
  allocate(a_d(na,na), a0_d(na,na), a1_d(na,na), q_d(na,na))
  ub_d=d_new_ub(na,lbwmaxa,ubwmaxa)
  bv_d=d_new_bv(na,lbwmaxa,ubwmaxa)
  sw_d=d_random_sweeps(na,numsweeps)
  call d_assemble_a(a_d,u_d(1:na,1:ubwa),v_d(1:ubwa,1:na),d_d(1:na),lbwa)
  a0_d=a_d
  call upper_to_ub(a_d,ub_d,lbwa, tol,error)
  a1_d=a0_d
  call sweeps_times_general(sw_d,a1_d)
  call cpu_time(t0)
  call d_sweeps_times_ub(sw_d,ub_d,bv_d,error)
  call cpu_time(t1)
  call d_bv_to_upper(bv_d,a_d,error)
  test_name = "Real Sweeps Times UB (complete fill)"
  a_d = a_d/maxabs(a_d)
  a1_d = a1_d/maxabs(a1_d)
  call d_output_result(test_name,a_d,a1_d,min(bv_d%ubw,ubwa+numsweeps),min(bv_d%ubw,ubwa+numsweeps), &
       t0,t1,tol2,error)
  deallocate(a_d, a0_d, a1_d, q_d)
  call deallocate_ub(ub_d); call deallocate_bv(bv_d)
  call deallocate_sweeps(sw_d)

  ! Real BV
  na=20; numsweeps=7
  lbwa=5; ubwa=5
  ubwmaxa=na+1; lbwmaxa=na+1
  u_d=u0_d; v_d=v0_d; d_d=d0_d
  allocate(a_d(na,na), a0_d(na,na), a1_d(na,na), q_d(na,na))
  ub_d=d_new_ub(na,lbwmaxa,ubwmaxa)
  bv_d=d_new_bv(na,lbwmaxa,ubwmaxa)
  sw_d=d_random_sweeps(na,numsweeps)
  call d_assemble_a(a_d,u_d(1:na,1:ubwa),v_d(1:ubwa,1:na),d_d(1:na),lbwa)
  a0_d=a_d
  call upper_to_bv(a_d,bv_d,lbwa, tol,error)
  a1_d=a0_d
  call general_times_sweeps(a1_d,sw_d)
  call cpu_time(t0)
  call d_bv_times_sweeps(bv_d, sw_d, ub_d, error)
  call cpu_time(t1)
  call d_ub_to_upper(ub_d,a_d,error)
  test_name = "Real BV Times Sweeps"
  call d_output_result(test_name,a_d,a1_d,min(ub_d%ubw,ubwa+numsweeps),min(ub_d%ubw,ubwa+numsweeps), &
       t0,t1,tol2,error)
  deallocate(a_d, a0_d, a1_d, q_d)
  call deallocate_ub(ub_d); call deallocate_bv(bv_d)
  call deallocate_sweeps(sw_d)


  na=20; numsweeps=20
  lbwa=5; ubwa=5
  ubwmaxa=na+1; lbwmaxa=na+1
  u_d=u0_d; v_d=v0_d; d_d=d0_d
  allocate(a_d(na,na), a0_d(na,na), a1_d(na,na), q_d(na,na))
  ub_d=d_new_ub(na,lbwmaxa,ubwmaxa)
  bv_d=d_new_bv(na,lbwmaxa,ubwmaxa)
  sw_d=d_random_sweeps(na,numsweeps)
  call d_assemble_a(a_d,u_d(1:na,1:ubwa),v_d(1:ubwa,1:na),d_d(1:na),lbwa)
  a0_d=a_d
  call upper_to_bv(a_d,bv_d,lbwa, tol,error)
  a1_d=a0_d
  call general_times_sweeps(a1_d,sw_d)
  call cpu_time(t0)
  call d_bv_times_sweeps(bv_d, sw_d, ub_d, error)
  call cpu_time(t1)
  call d_ub_to_upper(ub_d,a_d,error)
  test_name = "Real BV Times Sweeps (complete fill)"
  call d_output_result(test_name,a_d,a1_d,min(ub_d%ubw,ubwa+numsweeps),min(ub_d%ubw,ubwa+numsweeps), &
       t0,t1,tol2,error)
  deallocate(a_d, a0_d, a1_d, q_d)
  call deallocate_ub(ub_d); call deallocate_bv(bv_d)
  call deallocate_sweeps(sw_d)

  ! Complex

  print *
  print *, "--------------------------------"
  print *
  print *, "Complex Sweeps Tests."
  print *

  na=20; numsweeps=2
  lbwa=5; ubwa=5
  ubwmaxa=na+1; lbwmaxa=na+1
  u_c=u0_c; v_c=v0_c; d_c=d0_c
  allocate(a_c(na,na), a0_c(na,na), a1_c(na,na), q_c(na,na))
  ub_c=c_new_ub(na,lbwmaxa,ubwmaxa)
  bv_c=c_new_bv(na,lbwmaxa,ubwmaxa)
  sw_c=c_random_sweeps(na,numsweeps)
  call c_assemble_a(a_c,u_c(1:na,1:ubwa),v_c(1:ubwa,1:na),d_c(1:na),lbwa)
  a0_c=a_c
  call upper_to_ub(a_c,ub_c,lbwa, tol,error)
  a1_c=a0_c
  call sweeps_times_general(sw_c,a1_c)
  call cpu_time(t0)
  call c_sweeps_times_ub(sw_c,ub_c,bv_c,error)
  call cpu_time(t1)
  call c_bv_to_upper(bv_c,a_c,error)
  test_name = "Complex Sweeps Times UB"
  call c_output_result(test_name,a_c,a1_c,min(bv_c%ubw,ubwa+numsweeps),min(bv_c%ubw,ubwa+numsweeps), &
       t0,t1,tol2,error)
  deallocate(a_c, a0_c, a1_c, q_c)
  call deallocate_ub(ub_c); call deallocate_bv(bv_c)
  call deallocate_sweeps(sw_c)

  na=20; numsweeps=20
  lbwa=10; ubwa=10
  ubwmaxa=na+1; lbwmaxa=na+1
  u_c=u0_c; v_c=v0_c; d_c=d0_c
  allocate(a_c(na,na), a0_c(na,na), a1_c(na,na), q_c(na,na))
  ub_c=c_new_ub(na,lbwmaxa,ubwmaxa)
  bv_c=c_new_bv(na,lbwmaxa,ubwmaxa)
  sw_c=c_random_sweeps(na,numsweeps)
  call c_assemble_a(a_c,u_c(1:na,1:ubwa),v_c(1:ubwa,1:na),d_c(1:na),lbwa)
  a0_c=a_c
  call upper_to_ub(a_c,ub_c,lbwa, tol,error)
  a1_c=a0_c
  call sweeps_times_general(sw_c,a1_c)
  call cpu_time(t0)
  call c_sweeps_times_ub(sw_c,ub_c,bv_c,error)
  call cpu_time(t1)
  call c_bv_to_upper(bv_c,a_c,error)
  test_name = "Complex Sweeps Times UB"
  a_c = a_c/maxabs(a_c)
  a1_c = a1_c/maxabs(a1_c)
  call c_output_result(test_name,a_c,a1_c,min(bv_c%ubw,ubwa+numsweeps),min(bv_c%ubw,ubwa+numsweeps), &
       t0,t1,tol2,error)
  deallocate(a_c, a0_c, a1_c, q_c)
  call deallocate_ub(ub_c); call deallocate_bv(bv_c)
  call deallocate_sweeps(sw_c)

  ! BV

  na=20; numsweeps=7
  lbwa=5; ubwa=5
  ubwmaxa=na+1; lbwmaxa=na+1
  u_c=u0_c; v_c=v0_c; d_c=d0_c
  allocate(a_c(na,na), a0_c(na,na), a1_c(na,na), q_c(na,na))
  ub_c=c_new_ub(na,lbwmaxa,ubwmaxa)
  bv_c=c_new_bv(na,lbwmaxa,ubwmaxa)
  sw_c=c_random_sweeps(na,numsweeps)
  call c_assemble_a(a_c,u_c(1:na,1:ubwa),v_c(1:ubwa,1:na),d_c(1:na),lbwa)
  a0_c=a_c
  call upper_to_bv(a_c,bv_c,lbwa, tol,error)
  a1_c=a0_c
  call general_times_sweeps(a1_c,sw_c)
  call cpu_time(t0)
  call c_bv_times_sweeps(bv_c, sw_c, ub_c, error)
  call cpu_time(t1)
  call c_ub_to_upper(ub_c,a_c,error)
  test_name = "Complex BV Times Sweeps"
  call c_output_result(test_name,a_c,a1_c,min(ub_c%ubw,ubwa+numsweeps),min(ub_c%ubw,ubwa+numsweeps), &
       t0,t1,tol2,error)
  deallocate(a_c, a0_c, a1_c, q_c)
  call deallocate_ub(ub_c); call deallocate_bv(bv_c)
  call deallocate_sweeps(sw_c)

  na=20; numsweeps=20
  lbwa=5; ubwa=5
  ubwmaxa=na+1; lbwmaxa=na+1
  u_c=u0_c; v_c=v0_c; d_c=d0_c
  allocate(a_c(na,na), a0_c(na,na), a1_c(na,na), q_c(na,na))
  ub_c=c_new_ub(na,lbwmaxa,ubwmaxa)
  bv_c=c_new_bv(na,lbwmaxa,ubwmaxa)
  sw_c=c_random_sweeps(na,numsweeps)
  call c_assemble_a(a_c,u_c(1:na,1:ubwa),v_c(1:ubwa,1:na),d_c(1:na),lbwa)
  a0_c=a_c
  call upper_to_bv(a_c,bv_c,lbwa, tol,error)
  a1_c=a0_c
  call general_times_sweeps(a1_c,sw_c)
  call cpu_time(t0)
  call c_bv_times_sweeps(bv_c, sw_c, ub_c, error)
  call cpu_time(t1)
  call c_ub_to_upper(ub_c,a_c,error)
  test_name = "Complex BV Times Sweeps (Complete fill)"
  call c_output_result(test_name,a_c,a1_c,min(ub_c%ubw,ubwa+numsweeps),min(ub_c%ubw,ubwa+numsweeps), &
       t0,t1,tol2,error)
  deallocate(a_c, a0_c, a1_c, q_c)
  call deallocate_ub(ub_c); call deallocate_bv(bv_c)
  call deallocate_sweeps(sw_c)


contains

type(d_sweeps) function d_random_sweeps(n,l) result(sw)
  integer(kind=int32) :: n, l, j, k
  real(kind=dp) :: nrm
  sw=d_new_sweeps(n,l)
  sw%numsweeps=l
  call random_number(sw%cs)
  call random_number(sw%ss)
  do j=1,n
     do k=1,l
        nrm=sqrt(sw%cs(j,k)**2+sw%ss(j,k)**2)
        sw%cs(j,k)=sw%cs(j,k)/nrm
        sw%ss(j,k)=sw%ss(j,k)/nrm
     end do
  end do
end function d_random_sweeps

type(c_sweeps) function c_random_sweeps(n,l) result(sw)
  integer(kind=int32) :: n, l, j, k
  real(kind=dp) :: nrm
  sw=c_new_sweeps(n,l)
  sw%numsweeps=l
  call random_complex(sw%cs)
  call random_complex(sw%ss)
  do j=1,n
     do k=1,l
        nrm=sqrt(abs(sw%cs(j,k))**2+abs(sw%ss(j,k))**2)
        sw%cs(j,k)=sw%cs(j,k)/nrm
        sw%ss(j,k)=sw%ss(j,k)/nrm
     end do
  end do
end function c_random_sweeps


end program test_sweeps
