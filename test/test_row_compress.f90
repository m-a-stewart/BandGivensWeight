program test_row_compress
  use mod_orth_rank
  use mod_test_data
  implicit none

  ! Note: This does not cover large rank cases (e.g. ubw=n-1).  But the conversion
  ! routines are tested for these cases by test_sweeps1.f90, since sweeps1.f90 uses
  ! the conversion routines.

  real(kind=dp) :: t0, t1
  integer(kind=int32) :: na
  type(error_info) :: error
  integer, parameter :: n=50, rmaxl=13, lbwmax=rmaxl+1, rmaxu=11, ubwmax=rmaxu+rmaxl+1
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

  type(d_ubt) :: ubt_d, ubt_na_d
  type(c_ubt) :: ubt_c, ubt_na_c
  type(d_bv) :: bv_d, bv_na_d
  type(c_bv) :: bv_c, bv_na_c
  type(d_sweeps) :: sw_d, sw_na_d
  type(c_sweeps) :: sw_c, sw_na_c
  ubt_d=d_new_ubt(n,lbwmax,ubwmax)
  ubt_c=c_new_ubt(n,lbwmax,ubwmax)
  bv_d=d_new_bv(n,lbwmax,lbwmax+ubwmax)
  bv_c=c_new_bv(n,lbwmax,lbwmax+ubwmax)
  sw_d=d_new_sweeps(n,1,n-2,lbwmax)
  sw_c=c_new_sweeps(n,1,n-2,lbwmax)

  ! real
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
  print *
  print *, "--------------------------------"
  print *
  print *, "Real Row Compression Tests:"
  print *

  call d_assemble_general(a,ul,vl,uu,vu,d)
  a0=a
  call general_to_ubt(a,ubt_d,tol,error)
  call cpu_time(t0)
  call row_compress(ubt_d, bv_d, sw_d, error)
  call cpu_time(t1)
  call bv_to_upper(bv_d,a1,error)
  call sweeps_times_general(sw_d,a1)
  test_name = "Real Row Compression;"
  call d_output_result_lower_upper(test_name,a0,a1,rmaxl,bv_d%lbw, &
       rmaxu+rmaxl,bv_d%ubw,t0,t1,tol2,error)

  na=1
  ubt_na_d=d_new_ubt(na,lbwmax,ubwmax)
  bv_na_d=d_new_bv(na,lbwmax,ubwmax)
  sw_na_d=d_new_sweeps(na,1,na-2,lbwmax)
  ul=ul0; vl=vl0
  uu=uu0; vu=vu0
  call d_assemble_general(a(1:na,1:na),ul(1:na,:),vl(:,1:na),uu(1:na,:),vu(:,1:na),&
       d(1:na))
  a0(1:na,1:na)=a(1:na,1:na)
  call general_to_ubt(a(1:na,1:na),ubt_na_d,tol1,error)
  call cpu_time(t0)
  call row_compress(ubt_na_d, bv_na_d, sw_na_d,error)
  call cpu_time(t1)
  call bv_to_upper(bv_na_d,a1(1:na,1:na),error)
  call sweeps_times_general(sw_na_d,a1)
  test_name = "Complex Row Compression (n=1);"
  call d_output_result_lower_upper(test_name,a0(1:na,1:na),a1(1:na,1:na),0,bv_na_d%lbw,0, &
       bv_na_d%ubw,t0,t1,tol2,error)
  call d_deallocate_ubt(ubt_na_d)
  call d_deallocate_bv(bv_na_d)
  call d_deallocate_sweeps(sw_na_d)

  na=2
  ubt_na_d=d_new_ubt(na,lbwmax,ubwmax)
  bv_na_d=d_new_bv(na,lbwmax,ubwmax)
  sw_na_d=d_new_sweeps(na,1,na-2,lbwmax)
  ul=ul0; vl=vl0
  uu=uu0; vu=vu0
  call d_assemble_general(a(1:na,1:na),ul(1:na,:),vl(:,1:na),uu(1:na,:),vu(:,1:na),&
       d(1:na))
  a0(1:na,1:na)=a(1:na,1:na)
  call general_to_ubt(a(1:na,1:na),ubt_na_d,tol1,error)
  call cpu_time(t0)
  call row_compress(ubt_na_d, bv_na_d, sw_na_d,error)
  call cpu_time(t1)
  call bv_to_upper(bv_na_d,a1(1:na,1:na),error)
  call sweeps_times_general(sw_na_d,a1)
  test_name = "Complex Row Compression (n=2);"
  call d_output_result_lower_upper(test_name,a0(1:na,1:na),a1(1:na,1:na),1,bv_na_d%lbw,1, &
       bv_na_d%ubw,t0,t1,tol2,error)
  call d_deallocate_ubt(ubt_na_d)
  call d_deallocate_bv(bv_na_d)
  call d_deallocate_sweeps(sw_na_d)

  na=3
  ubt_na_d=d_new_ubt(na,lbwmax,ubwmax)
  bv_na_d=d_new_bv(na,lbwmax,ubwmax)
  sw_na_d=d_new_sweeps(na,1,na-2,lbwmax)
  ul=ul0; vl=vl0
  uu=uu0; vu=vu0
  call d_assemble_general(a(1:na,1:na),ul(1:na,:),vl(:,1:na),uu(1:na,:),vu(:,1:na),&
       d(1:na))
  a0(1:na,1:na)=a(1:na,1:na)
  call general_to_ubt(a(1:na,1:na),ubt_na_d,tol1,error)
  call cpu_time(t0)
  call row_compress(ubt_na_d, bv_na_d, sw_na_d,error)
  call cpu_time(t1)
  call bv_to_upper(bv_na_d,a1(1:na,1:na),error)
  call sweeps_times_general(sw_na_d,a1(1:na,1:na))
  test_name = "Real Row Compression (n=3);"
  call d_output_result_lower_upper(test_name,a0(1:na,1:na),a1(1:na,1:na),1,bv_na_d%lbw,2, &
       bv_na_d%ubw,t0,t1,tol2,error)

  call d_deallocate_ubt(ubt_na_d)
  call d_deallocate_bv(bv_na_d)
  call d_deallocate_sweeps(sw_na_d)

  na=4
  ubt_na_d=d_new_ubt(na,lbwmax,ubwmax)
  bv_na_d=d_new_bv(na,lbwmax,ubwmax)
  sw_na_d=d_new_sweeps(na,1,na-2,lbwmax)
  ul=ul0; vl=vl0
  uu=uu0; vu=vu0
  call d_assemble_general(a(1:na,1:na),ul(1:na,:),vl(:,1:na),uu(1:na,:),vu(:,1:na),&
       d(1:na))
  a0(1:na,1:na)=a(1:na,1:na)
  call general_to_ubt(a(1:na,1:na),ubt_na_d,tol1,error)
  call cpu_time(t0)
  call row_compress(ubt_na_d, bv_na_d, sw_na_d,error)
  call cpu_time(t1)
  call bv_to_upper(bv_na_d,a1(1:na,1:na),error)
  call sweeps_times_general(sw_na_d,a1(1:na,1:na))
  test_name = "Real Row Compression (n=4);"
  call d_output_result_lower_upper(test_name,a0(1:na,1:na),a1(1:na,1:na),2,bv_na_d%lbw,3, &
       bv_na_d%ubw,t0,t1,tol2,error)

  call d_deallocate_ubt(ubt_na_d)
  call d_deallocate_bv(bv_na_d)
  call d_deallocate_sweeps(sw_na_d)

  

  print *
  print *, "--------------------------------"
  print *
  print *, "Complex Row Compression Tests:"
  print *

  call c_assemble_general(a_c,ul_c,vl_c,uu_c,vu_c,d_c)
  a0_c=a_c
  call general_to_ubt(a_c,ubt_c,tol,error)
  call cpu_time(t0)
  call row_compress(ubt_c, bv_c, sw_c, error)
  call cpu_time(t1)
  call bv_to_upper(bv_c,a1_c,error)
  call sweeps_times_general(sw_c,a1_c)
  test_name = "Complex Row Compression;"
  call c_output_result_lower_upper(test_name,a0_c,a1_c,rmaxl,bv_c%lbw, &
       rmaxu+rmaxl,bv_c%ubw,t0,t1,tol2,error)
  
  na=1
  ubt_na_c=c_new_ubt(na,lbwmax,ubwmax)
  bv_na_c=c_new_bv(na,lbwmax,ubwmax)
  sw_na_c=c_new_sweeps(na,1,na-2,lbwmax)
  ul_c=ul0_c; vl_c=vl0_c
  uu_c=uu0_c; vu_c=vu0_c
  call c_assemble_general(a_c(1:na,1:na),ul_c(1:na,:),vl_c(:,1:na),uu_c(1:na,:),vu_c(:,1:na),&
       d_c(1:na))
  a0_c(1:na,1:na)=a_c(1:na,1:na)
  call general_to_ubt(a_c(1:na,1:na),ubt_na_c,tol1,error)
  call cpu_time(t0)
  call row_compress(ubt_na_c, bv_na_c, sw_na_c,error)
  call cpu_time(t1)
  call bv_to_upper(bv_na_c,a1_c(1:na,1:na),error)
  call sweeps_times_general(sw_na_c,a1_c(1:na,1:na))
  test_name = "Complex Row Compression (n=1);"
  call c_output_result_lower_upper(test_name,a0_c(1:na,1:na),a1_c(1:na,1:na),0,bv_na_c%lbw,0, &
       bv_na_c%ubw,t0,t1,tol2,error)
  call c_deallocate_ubt(ubt_na_c)
  call c_deallocate_bv(bv_na_c)
  call c_deallocate_sweeps(sw_na_c)

  na=2
  ubt_na_c=c_new_ubt(na,lbwmax,ubwmax)
  bv_na_c=c_new_bv(na,lbwmax,ubwmax)
  sw_na_c=c_new_sweeps(na,1,na-2,lbwmax)
  ul_c=ul0_c; vl_c=vl0_c
  uu_c=uu0_c; vu_c=vu0_c
  call c_assemble_general(a_c(1:na,1:na),ul_c(1:na,:),vl_c(:,1:na),uu_c(1:na,:),vu_c(:,1:na),&
       d_c(1:na))
  a0_c(1:na,1:na)=a_c(1:na,1:na)
  call general_to_ubt(a_c(1:na,1:na),ubt_na_c,tol1,error)
  call cpu_time(t0)
  call row_compress(ubt_na_c, bv_na_c, sw_na_c,error)
  call cpu_time(t1)
  call bv_to_upper(bv_na_c,a1_c(1:na,1:na),error)
  call sweeps_times_general(sw_na_c,a1_c(1:na,1:na))
  test_name = "Complex Row Compression (n=2);"
  call c_output_result_lower_upper(test_name,a0_c(1:na,1:na),a1_c(1:na,1:na),1,bv_na_c%lbw,1, &
       bv_na_c%ubw,t0,t1,tol2,error)
  call c_deallocate_ubt(ubt_na_c)
  call c_deallocate_bv(bv_na_c)
  call c_deallocate_sweeps(sw_na_c)

  na=3
  ubt_na_c=c_new_ubt(na,lbwmax,ubwmax)
  bv_na_c=c_new_bv(na,lbwmax,ubwmax)
  sw_na_c=c_new_sweeps(na,1,na-2,lbwmax)
  ul_c=ul0_c; vl_c=vl0_c
  uu_c=uu0_c; vu_c=vu0_c
  call c_assemble_general(a_c(1:na,1:na),ul_c(1:na,:),vl_c(:,1:na),uu_c(1:na,:),vu_c(:,1:na),&
       d_c(1:na))
  a0_c(1:na,1:na)=a_c(1:na,1:na)
  call general_to_ubt(a_c(1:na,1:na),ubt_na_c,tol1,error)
  call cpu_time(t0)
  call row_compress(ubt_na_c, bv_na_c, sw_na_c,error)
  call cpu_time(t1)
  call bv_to_upper(bv_na_c,a1_c(1:na,1:na),error)
  call sweeps_times_general(sw_na_c,a1_c(1:na,1:na))
  test_name = "Complex Row Compression (n=2);"
  call c_output_result_lower_upper(test_name,a0_c(1:na,1:na),a1_c(1:na,1:na),1,bv_na_c%lbw,2, &
       bv_na_c%ubw,t0,t1,tol2,error)
  call c_deallocate_ubt(ubt_na_c)
  call c_deallocate_bv(bv_na_c)
  call c_deallocate_sweeps(sw_na_c)

  na=4
  ubt_na_c=c_new_ubt(na,lbwmax,ubwmax)
  bv_na_c=c_new_bv(na,lbwmax,ubwmax)
  sw_na_c=c_new_sweeps(na,1,na-2,lbwmax)
  ul_c=ul0_c; vl_c=vl0_c
  uu_c=uu0_c; vu_c=vu0_c
  call c_assemble_general(a_c(1:na,1:na),ul_c(1:na,:),vl_c(:,1:na),uu_c(1:na,:),vu_c(:,1:na),&
       d_c(1:na))
  a0_c(1:na,1:na)=a_c(1:na,1:na)
  call general_to_ubt(a_c(1:na,1:na),ubt_na_c,tol1,error)
  call cpu_time(t0)
  call row_compress(ubt_na_c, bv_na_c, sw_na_c,error)
  call cpu_time(t1)
  call bv_to_upper(bv_na_c,a1_c(1:na,1:na),error)
  call sweeps_times_general(sw_na_c,a1_c(1:na,1:na))
  test_name = "Complex Row Compression (n=4);"
  call c_output_result_lower_upper(test_name,a0_c(1:na,1:na),a1_c(1:na,1:na),2,bv_na_c%lbw,3, &
       bv_na_c%ubw,t0,t1,tol2,error)
  call c_deallocate_ubt(ubt_na_c)
  call c_deallocate_bv(bv_na_c)
  call c_deallocate_sweeps(sw_na_c)



end program test_row_compress
