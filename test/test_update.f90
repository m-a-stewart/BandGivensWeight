program test_update
  use nested
  use test_data
  implicit none
  real(kind=dp) :: t0, t1
  integer(kind=int32) :: na, lbwa, ubwa, j, k, lbwmaxa, ubwmaxa, numsweeps
  integer(kind=int32), parameter :: nmax=1000, rmaxa=50
  type(error_info) :: error
  real(kind=dp), parameter :: tol=1e-14, tol1=1e-14, tol2=1e-10
  !
  real(kind=dp), dimension(nmax,rmaxa) :: u_d, u0_d
  real(kind=dp), dimension(rmaxa,nmax) :: v_d, v0_d
  real(kind=dp), dimension(nmax) :: d_d, d0_d
  real(kind=dp), dimension(nmax), target :: w_d, x_d, w0_d, x0_d
  real(kind=dp), dimension(:), pointer :: wp_d, xp_d
  complex(kind=dp), dimension(nmax,rmaxa) :: u_c, u0_c
  complex(kind=dp), dimension(rmaxa,nmax) :: v_c, v0_c
  complex(kind=dp), dimension(nmax) :: d_c, d0_c
  complex(kind=dp), dimension(nmax), target :: w_c, x_c, w0_c, x0_c
  complex(kind=dp), dimension(:), pointer :: wp_c, xp_c

  real(kind=dp), dimension(:,:), allocatable :: a_d, a0_d, a1_d
  complex(kind=dp), dimension(:,:), allocatable :: a_c, a0_c, a1_c

  real(kind=dp) :: b11_d, b12_d
  complex(kind=dp) :: b11_c, b12_c

  type(d_ub) :: ub_d
  type(c_ub) :: ub_c
  type(d_bv) :: bv_d



  type(c_bv) :: bv_c
  type(d_sweeps) :: sw_d
  type(c_sweeps) :: sw_c

  call random_seed
  call random_matrix(u_d)
  call random_matrix(v_d)
  call random_matrix(d_d)
  call random_matrix(w_d)
  call random_matrix(x_d)
  u0_d=u_d; v0_d=v_d; d0_d=d_d
  w0_d=w_d; x0_d=x_d

  call random_matrix(u_c)
  call random_matrix(v_c)
  call random_matrix(d_c)
  call random_matrix(w_c)
  call random_matrix(x_c)
  u0_c=u_c; v0_c=v_c; d0_c=d_c
  w0_c=w_c; x0_c=x_c

  print *
  print *, "--------------------------------"
  print *
  print *, "Real Rank One Update Tests."
  print *

  na=10; numsweeps=1
  lbwa=2; ubwa=5
  ubwmaxa=ubwa+3
  lbwmaxa=lbwa+1
  u_d=u0_d; v_d=v0_d; d_d=d0_d
  w_d=w0_d; x_d=x0_d
  allocate(a_d(na,na), a0_d(na,na), a1_d(na,na))
  ub_d=d_new_ub(na,lbwmaxa,ubwmaxa)
  bv_d=d_new_bv(na,lbwmaxa,ubwmaxa)
  sw_d=d_new_sweeps(na,numsweeps)
  call d_assemble_upper(a_d,u_d(1:na,1:ubwa),v_d(1:ubwa,1:na),d_d(1:na),lbwa)
  a0_d=a_d
  wp_d => w_d(1:na); xp_d => x_d(1:na)
  call upper_to_ub(a_d,ub_d,lbwa, tol,error)
  call cpu_time(t0)
  call d_r1_update_ub_to_bv(ub_d, wp_d, xp_d, sw_d, bv_d, error)
  call cpu_time(t1)
  b11_d=get_el_br(bv_d%br, bv_d%lbw,1,1)
  call set_el_br(bv_d%br,bv_d%lbw,1,1,b11_d+wp_d(1)*xp_d(1))
  b12_d=get_el_br(bv_d%br, bv_d%lbw,1,2)
  call set_el_br(bv_d%br,bv_d%lbw,1,2,b12_d+wp_d(1)*xp_d(2))
  do j=1,na
     do k=1,na
        a1_d(j,k)=a0_d(j,k)+w0_d(j)* x0_d(k)
     end do
  end do
  call sweeps_times_general(sw_d,a1_d)
  call d_bv_to_upper(bv_d,a_d,error)  
  test_name = "Real Rank 1 Update"
  call d_output_result_upper(test_name,a_d,a1_d,bv_d%ubw,bv_d%ubw,t0,t1,tol2,error)
  deallocate(a_d, a0_d, a1_d)
  call deallocate_ub(ub_d); call deallocate_bv(bv_d)
  call deallocate_sweeps(sw_d)
  !
  ! full lbw and ubw
  !
  na=10; numsweeps=10
  lbwa=na-1; ubwa=na-1
  ubwmaxa=ubwa+3
  lbwmaxa=lbwa+1
  u_d=u0_d; v_d=v0_d; d_d=d0_d
  w_d=w0_d; x_d=x0_d
  allocate(a_d(na,na), a0_d(na,na), a1_d(na,na))
  ub_d=d_new_ub(na,lbwmaxa,ubwmaxa)
  bv_d=d_new_bv(na,lbwmaxa,ubwmaxa)
  sw_d=d_random_sweeps(na,numsweeps)
  call d_assemble_upper(a_d,u_d(1:na,1:ubwa),v_d(1:ubwa,1:na),d_d(1:na),lbwa)
  wp_d => w_d(1:na); xp_d => x_d(1:na)
  call upper_to_ub(a_d,ub_d,lbwa, tol,error)
  call sweeps_times_ub(sw_d, ub_d, bv_d, error)
  call convert_bv_to_ub(bv_d, ub_d, error)
  call ub_to_upper(ub_d,a_d,error)
  a0_d=a_d
  sw_d%numsweeps=0
  call cpu_time(t0)
  call d_r1_update_ub_to_bv(ub_d, wp_d, xp_d, sw_d, bv_d, error)
  call cpu_time(t1)
  b11_d=get_el_br(bv_d%br, bv_d%lbw,1,1)
  call set_el_br(bv_d%br,bv_d%lbw,1,1,b11_d+wp_d(1)*xp_d(1))
  b12_d=get_el_br(bv_d%br, bv_d%lbw,1,2)
  call set_el_br(bv_d%br,bv_d%lbw,1,2,b12_d+wp_d(1)*xp_d(2))
  do j=1,na
     do k=1,na
        a1_d(j,k)=a0_d(j,k)+w0_d(j)* x0_d(k)
     end do
  end do
  call sweeps_times_general(sw_d,a1_d)
  call d_bv_to_upper(bv_d,a_d,error)  
  test_name = "Real Rank 1 Update (full lbw and ubw)"
  call d_output_result_upper(test_name,a_d,a1_d,bv_d%ubw,bv_d%ubw,t0,t1,tol2,error)
  deallocate(a_d, a0_d, a1_d)
  call deallocate_ub(ub_d); call deallocate_bv(bv_d)
  call deallocate_sweeps(sw_d)


  print *
  print *, "--------------------------------"
  print *
  print *, "Complex Rank One Update Tests."
  print *

  na=10; numsweeps=1
  lbwa=2; ubwa=5
  ubwmaxa=ubwa+3
  lbwmaxa=lbwa+1
  u_c=u0_c; v_c=v0_c; d_c=d0_c
  w_c=w0_c; x_c=x0_c
  allocate(a_c(na,na), a0_c(na,na), a1_c(na,na))
  ub_c=c_new_ub(na,lbwmaxa,ubwmaxa)
  bv_c=c_new_bv(na,lbwmaxa,ubwmaxa)
  sw_c=c_new_sweeps(na,numsweeps)
  call c_assemble_upper(a_c,u_c(1:na,1:ubwa),v_c(1:ubwa,1:na),d_c(1:na),lbwa)
  a0_c=a_c
  wp_c => w_c(1:na); xp_c => x_c(1:na)
  call upper_to_ub(a_c,ub_c,lbwa, tol,error)
  call cpu_time(t0)
  call c_r1_update_ub_to_bv(ub_c, wp_c, xp_c, sw_c, bv_c, error)
  call cpu_time(t1)
  b11_c=get_el_br(bv_c%br, bv_c%lbw,1,1)
  call set_el_br(bv_c%br,bv_c%lbw,1,1,b11_c+wp_c(1)*conjg(xp_c(1)))
  b12_c=get_el_br(bv_c%br, bv_c%lbw,1,2)
  call set_el_br(bv_c%br,bv_c%lbw,1,2,b12_c+wp_c(1)*conjg(xp_c(2)))
  do j=1,na
     do k=1,na
        a1_c(j,k)=a0_c(j,k)+w0_c(j)* conjg(x0_c(k))
     end do
  end do
  call sweeps_times_general(sw_c,a1_c)
  call c_bv_to_upper(bv_c,a_c,error)  
  test_name = "Real Rank 1 Update"
  call c_output_result_upper(test_name,a_c,a1_c,bv_c%ubw,bv_c%ubw,t0,t1,tol2,error)
  deallocate(a_c, a0_c, a1_c)
  call deallocate_ub(ub_c); call deallocate_bv(bv_c)
  call deallocate_sweeps(sw_c)
  !
  ! full lbw and ubw
  !
  na=10; numsweeps=10
  lbwa=na-1; ubwa=na-1
  ubwmaxa=ubwa+3
  lbwmaxa=lbwa+1
  u_c=u0_c; v_c=v0_c; d_c=d0_c
  w_c=w0_c; x_c=x0_c
  allocate(a_c(na,na), a0_c(na,na), a1_c(na,na))
  ub_c=c_new_ub(na,lbwmaxa,ubwmaxa)
  bv_c=c_new_bv(na,lbwmaxa,ubwmaxa)
  sw_c=c_random_sweeps(na,numsweeps)
  call c_assemble_upper(a_c,u_c(1:na,1:ubwa),v_c(1:ubwa,1:na),d_c(1:na),lbwa)
  wp_c => w_c(1:na); xp_c => x_c(1:na)
  call upper_to_ub(a_c,ub_c,lbwa, tol,error)
  call sweeps_times_ub(sw_c, ub_c, bv_c, error)
  call convert_bv_to_ub(bv_c, ub_c, error)
  call ub_to_upper(ub_c,a_c,error)
  a0_c=a_c
  sw_c%numsweeps=0
  call cpu_time(t0)
  call c_r1_update_ub_to_bv(ub_c, wp_c, xp_c, sw_c, bv_c, error)
  call cpu_time(t1)
  b11_c=get_el_br(bv_c%br, bv_c%lbw,1,1)
  call set_el_br(bv_c%br,bv_c%lbw,1,1,b11_c+wp_c(1)*conjg(xp_c(1)))
  b12_c=get_el_br(bv_c%br, bv_c%lbw,1,2)
  call set_el_br(bv_c%br,bv_c%lbw,1,2,b12_c+wp_c(1)*conjg(xp_c(2)))
  do j=1,na
     do k=1,na
        a1_c(j,k)=a0_c(j,k)+w0_c(j)* conjg(x0_c(k))
     end do
  end do
  call sweeps_times_general(sw_c,a1_c)
  call c_bv_to_upper(bv_c,a_c,error)  
  test_name = "Real Rank 1 Update (full lbw and ubw)"
  call c_output_result_upper(test_name,a_c,a1_c,bv_c%ubw,bv_c%ubw,t0,t1,tol2,error)
  deallocate(a_c, a0_c, a1_c)
  call deallocate_ub(ub_c); call deallocate_bv(bv_c)
  call deallocate_sweeps(sw_c)

end program test_update
