program test_convert_ubt_and_wbv
  use mod_orth_rank
  use mod_test_data
  implicit none

  !
  ! Note: This does not cover large rank cases (e.g. ubw=n-1).  But the conversion
  ! routines are tested for these cases by test_sweeps1.f90, since sweeps1.f90 uses
  ! the conversion routines.
  !

  real(kind=dp) :: t0, t1
  integer(kind=int32) :: na
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

  type(d_ubt) :: ubt_d, ubt_na_d
  type(c_ubt) :: ubt_c, ubt_na_c
  type(d_wbv) :: wbv_d, wbv_na_d
  type(c_wbv) :: wbv_c, wbv_na_c
  ubt_d=d_new_ubt(n,lbwmax,ubwmax)
  ubt_c=c_new_ubt(n,lbwmax,ubwmax)
  wbv_d=d_new_wbv(n,lbwmax,ubwmax)
  wbv_c=c_new_wbv(n,lbwmax,ubwmax)

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
  print *, "Real UBT and WBV Conversion Tests:"
  print *
  ! ubt to wbv
  call d_assemble_general(a,ul,vl,uu,vu,d)
  a0=a
  call general_to_ubt(a,ubt_d,tol,error)
  call cpu_time(t0)
  call convert_ubt_to_wbv(ubt_d, wbv_d, error)
  call cpu_time(t1)
  call wbv_to_general(wbv_d,a1,error)
  test_name = "Real UBT to WBV;"
  call d_output_result_lower_upper(test_name,a0,a1,rmaxl,wbv_d%lbw, &
       rmaxu,wbv_d%ubw,t0,t1,tol2,error)
  
  na=1
  ubt_na_d=d_new_ubt(na,lbwmax,ubwmax)
  wbv_na_d=d_new_wbv(na,lbwmax,ubwmax)
  ul=ul0; vl=vl0
  uu=uu0; vu=vu0
  call d_assemble_general(a(1:na,1:na),ul(1:na,:),vl(:,1:na),uu(1:na,:),vu(:,1:na),d(1:na))
  a0(1:na,1:na)=a(1:na,1:na)
  call general_to_ubt(a(1:na,1:na),ubt_na_d,tol1,error)
  call cpu_time(t0)
  call convert_ubt_to_wbv(ubt_na_d, wbv_na_d,error)
  call cpu_time(t1)
  call wbv_to_general(wbv_na_d,a1(1:na,1:na),error)
  test_name = "Real UBT to WBV (n=1);"
  call d_output_result_lower_upper(test_name,a0(1:na,1:na),a1(1:na,1:na),0,wbv_na_d%lbw,0, &
       wbv_na_d%ubw,t0,t1,tol2,error)
  call d_deallocate_ubt(ubt_na_d)
  call d_deallocate_wbv(wbv_na_d)

  na=2
  ubt_na_d=d_new_ubt(na,lbwmax,ubwmax)
  wbv_na_d=d_new_wbv(na,lbwmax,ubwmax)
  ul=ul0; vl=vl0
  uu=uu0; vu=vu0
  call d_assemble_general(a(1:na,1:na),ul(1:na,:),vl(:,1:na),uu(1:na,:),vu(:,1:na),d(1:na))
  a0(1:na,1:na)=a(1:na,1:na)
  call general_to_ubt(a(1:na,1:na),ubt_na_d,tol1,error)
  call cpu_time(t0)
  call convert_ubt_to_wbv(ubt_na_d, wbv_na_d,error)
  call cpu_time(t1)
  call wbv_to_general(wbv_na_d,a1(1:na,1:na),error)
  test_name = "Real UBT to WBV (n=2);"
  call d_output_result_lower_upper(test_name,a0(1:na,1:na),a1(1:na,1:na),1,wbv_na_d%lbw,1, &
       wbv_na_d%ubw,t0,t1,tol2,error)
  call d_deallocate_ubt(ubt_na_d)
  call d_deallocate_wbv(wbv_na_d)

  na=3
  ubt_na_d=d_new_ubt(na,lbwmax,ubwmax)
  wbv_na_d=d_new_wbv(na,lbwmax,ubwmax)
  ul=ul0; vl=vl0
  uu=uu0; vu=vu0
  call d_assemble_general(a(1:na,1:na),ul(1:na,:),vl(:,1:na),uu(1:na,:),vu(:,1:na),d(1:na))
  a0(1:na,1:na)=a(1:na,1:na)
  call general_to_ubt(a(1:na,1:na),ubt_na_d,tol1,error)
  call cpu_time(t0)
  call convert_ubt_to_wbv(ubt_na_d, wbv_na_d,error)
  call cpu_time(t1)
  call wbv_to_general(wbv_na_d,a1(1:na,1:na),error)
  test_name = "Real UBT to WBV (n=3);"
  call d_output_result_lower_upper(test_name,a0(1:na,1:na),a1(1:na,1:na),1,wbv_na_d%lbw,1, &
       wbv_na_d%ubw,t0,t1,tol2,error)
  call d_deallocate_ubt(ubt_na_d)
  call d_deallocate_wbv(wbv_na_d)

  na=4
  ubt_na_d=d_new_ubt(na,lbwmax,ubwmax)
  wbv_na_d=d_new_wbv(na,lbwmax,ubwmax)
  ul=ul0; vl=vl0
  uu=uu0; vu=vu0
  call d_assemble_general(a(1:na,1:na),ul(1:na,:),vl(:,1:na),uu(1:na,:),vu(:,1:na),d(1:na))
  a0(1:na,1:na)=a(1:na,1:na)
  call general_to_ubt(a(1:na,1:na),ubt_na_d,tol1,error)
  call cpu_time(t0)
  call convert_ubt_to_wbv(ubt_na_d, wbv_na_d,error)
  call cpu_time(t1)
  call wbv_to_general(wbv_na_d,a1(1:na,1:na),error)
  test_name = "Real UBT to WBV (n=2);"
  call d_output_result_lower_upper(test_name,a0(1:na,1:na),a1(1:na,1:na),2,wbv_na_d%lbw,2, &
       wbv_na_d%ubw,t0,t1,tol2,error)
  call d_deallocate_ubt(ubt_na_d)
  call d_deallocate_wbv(wbv_na_d)

  ! wbv to ubt
  ul=ul0; vl=vl0
  uu=uu0; vu=vu0
  call d_assemble_general(a,ul,vl,uu,vu,d)
  a0=a
  call general_to_wbv(a,wbv_d,tol,error)
  call cpu_time(t0)
  call convert_wbv_to_ubt(wbv_d, ubt_d, error)
  call cpu_time(t1)
  call ubt_to_general(ubt_d,a1,error)
  test_name = "Real WBV to UBT;"
  call d_output_result_lower_upper(test_name,a0,a1,rmaxl,ubt_d%lbw, &
       rmaxu,ubt_d%ubw,t0,t1,tol2,error)

  na=1
  ubt_na_d=d_new_ubt(na,lbwmax,ubwmax)
  wbv_na_d=d_new_wbv(na,lbwmax,ubwmax)
  ul=ul0; vl=vl0
  uu=uu0; vu=vu0
  call d_assemble_general(a(1:na,1:na),ul(1:na,:),vl(:,1:na),uu(1:na,:),vu(:,1:na),d(1:na))
  a0(1:na,1:na)=a(1:na,1:na)
  call general_to_wbv(a(1:na,1:na),wbv_na_d,tol1,error)
  call cpu_time(t0)
  call convert_wbv_to_ubt(wbv_na_d, ubt_na_d,error)
  call cpu_time(t1)
  call ubt_to_general(ubt_na_d,a1(1:na,1:na),error)
  test_name = "Real WBV to UBT (n=1);"
  call d_output_result_lower_upper(test_name,a0(1:na,1:na),a1(1:na,1:na),0,ubt_na_d%lbw,0, &
       ubt_na_d%ubw,t0,t1,tol2,error)
  call d_deallocate_ubt(ubt_na_d)
  call d_deallocate_wbv(wbv_na_d)

  na=2
  ubt_na_d=d_new_ubt(na,lbwmax,ubwmax)
  wbv_na_d=d_new_wbv(na,lbwmax,ubwmax)
  ul=ul0; vl=vl0
  uu=uu0; vu=vu0
  call d_assemble_general(a(1:na,1:na),ul(1:na,:),vl(:,1:na),uu(1:na,:),vu(:,1:na),d(1:na))
  a0(1:na,1:na)=a(1:na,1:na)
  call general_to_wbv(a(1:na,1:na),wbv_na_d,tol1,error)
  call cpu_time(t0)
  call convert_wbv_to_ubt(wbv_na_d, ubt_na_d,error)
  call cpu_time(t1)
  call ubt_to_general(ubt_na_d,a1(1:na,1:na),error)
  test_name = "Real WBV to UBT (n=2);"
  call d_output_result_lower_upper(test_name,a0(1:na,1:na),a1(1:na,1:na),1,ubt_na_d%lbw,1, &
       ubt_na_d%ubw,t0,t1,tol2,error)
  call d_deallocate_ubt(ubt_na_d)
  call d_deallocate_wbv(wbv_na_d)

  na=3
  ubt_na_d=d_new_ubt(na,lbwmax,ubwmax)
  wbv_na_d=d_new_wbv(na,lbwmax,ubwmax)
  ul=ul0; vl=vl0
  uu=uu0; vu=vu0
  call d_assemble_general(a(1:na,1:na),ul(1:na,:),vl(:,1:na),uu(1:na,:),vu(:,1:na),d(1:na))
  a0(1:na,1:na)=a(1:na,1:na)
  call general_to_wbv(a(1:na,1:na),wbv_na_d,tol1,error)
  call cpu_time(t0)
  call convert_wbv_to_ubt(wbv_na_d, ubt_na_d,error)
  call cpu_time(t1)
  call ubt_to_general(ubt_na_d,a1(1:na,1:na),error)
  test_name = "Real WBV to UBT (n=3);"
  call d_output_result_lower_upper(test_name,a0(1:na,1:na),a1(1:na,1:na),1,ubt_na_d%lbw,1, &
       ubt_na_d%ubw,t0,t1,tol2,error)
  call d_deallocate_ubt(ubt_na_d)
  call d_deallocate_wbv(wbv_na_d)

  na=4
  ubt_na_d=d_new_ubt(na,lbwmax,ubwmax)
  wbv_na_d=d_new_wbv(na,lbwmax,ubwmax)
  ul=ul0; vl=vl0
  uu=uu0; vu=vu0
  call d_assemble_general(a(1:na,1:na),ul(1:na,:),vl(:,1:na),uu(1:na,:),vu(:,1:na),d(1:na))
  a0(1:na,1:na)=a(1:na,1:na)
  call general_to_wbv(a(1:na,1:na),wbv_na_d,tol1,error)
  call cpu_time(t0)
  call convert_wbv_to_ubt(wbv_na_d, ubt_na_d,error)
  call cpu_time(t1)
  call ubt_to_general(ubt_na_d,a1(1:na,1:na),error)
  test_name = "Real WBV to UBT (n=4);"
  call d_output_result_lower_upper(test_name,a0(1:na,1:na),a1(1:na,1:na),2,ubt_na_d%lbw,2, &
       ubt_na_d%ubw,t0,t1,tol2,error)
  call d_deallocate_ubt(ubt_na_d)
  call d_deallocate_wbv(wbv_na_d)

  print *
  print *, "--------------------------------"
  print *
  print *, "Complex UBT and WBV Conversion Tests:"
  print *

  call c_assemble_general(a_c,ul_c,vl_c,uu_c,vu_c,d_c)
  a0_c=a_c
  call general_to_ubt(a_c,ubt_c,tol,error)
  call cpu_time(t0)
  call convert_ubt_to_wbv(ubt_c, wbv_c, error)
  call cpu_time(t1)
  call wbv_to_general(wbv_c,a1_c,error)
  test_name = "Complex UBT to WBV;"
  call c_output_result_lower_upper(test_name,a0_c,a1_c,rmaxl,wbv_c%lbw, &
       rmaxu,wbv_c%ubw,t0,t1,tol2,error)
  
  na=1
  ubt_na_c=c_new_ubt(na,lbwmax,ubwmax)
  wbv_na_c=c_new_wbv(na,lbwmax,ubwmax)
  ul_c=ul0_c; vl_c=vl0_c
  uu_c=uu0_c; vu_c=vu0_c
  call c_assemble_general(a_c(1:na,1:na),ul_c(1:na,:),vl_c(:,1:na),uu_c(1:na,:),vu_c(:,1:na),&
       d_c(1:na))
  a0_c(1:na,1:na)=a_c(1:na,1:na)
  call general_to_ubt(a_c(1:na,1:na),ubt_na_c,tol1,error)
  call cpu_time(t0)
  call convert_ubt_to_wbv(ubt_na_c, wbv_na_c,error)
  call cpu_time(t1)
  call wbv_to_general(wbv_na_c,a1_c(1:na,1:na),error)
  test_name = "Complex UBT to WBV (n=1);"
  call c_output_result_lower_upper(test_name,a0_c(1:na,1:na),a1_c(1:na,1:na),0,wbv_na_c%lbw,0, &
       wbv_na_c%ubw,t0,t1,tol2,error)
  call c_deallocate_ubt(ubt_na_c)
  call c_deallocate_wbv(wbv_na_c)

  na=2
  ubt_na_c=c_new_ubt(na,lbwmax,ubwmax)
  wbv_na_c=c_new_wbv(na,lbwmax,ubwmax)
  ul_c=ul0_c; vl_c=vl0_c
  uu_c=uu0_c; vu_c=vu0_c
  call c_assemble_general(a_c(1:na,1:na),ul_c(1:na,:),vl_c(:,1:na),uu_c(1:na,:),vu_c(:,1:na),&
       d_c(1:na))
  a0_c(1:na,1:na)=a_c(1:na,1:na)
  call general_to_ubt(a_c(1:na,1:na),ubt_na_c,tol1,error)
  call cpu_time(t0)
  call convert_ubt_to_wbv(ubt_na_c, wbv_na_c,error)
  call cpu_time(t1)
  call wbv_to_general(wbv_na_c,a1_c(1:na,1:na),error)
  test_name = "Complex UBT to WBV (n=2);"
  call c_output_result_lower_upper(test_name,a0_c(1:na,1:na),a1_c(1:na,1:na),1,wbv_na_c%lbw,1, &
       wbv_na_c%ubw,t0,t1,tol2,error)
  call c_deallocate_ubt(ubt_na_c)
  call c_deallocate_wbv(wbv_na_c)

  na=3
  ubt_na_c=c_new_ubt(na,lbwmax,ubwmax)
  wbv_na_c=c_new_wbv(na,lbwmax,ubwmax)
  ul_c=ul0_c; vl_c=vl0_c
  uu_c=uu0_c; vu_c=vu0_c
  call c_assemble_general(a_c(1:na,1:na),ul_c(1:na,:),vl_c(:,1:na),uu_c(1:na,:),vu_c(:,1:na),&
       d_c(1:na))
  a0_c(1:na,1:na)=a_c(1:na,1:na)
  call general_to_ubt(a_c(1:na,1:na),ubt_na_c,tol1,error)
  call cpu_time(t0)
  call convert_ubt_to_wbv(ubt_na_c, wbv_na_c,error)
  call cpu_time(t1)
  call wbv_to_general(wbv_na_c,a1_c(1:na,1:na),error)
  test_name = "Complex UBT to WBV (n=3);"
  call c_output_result_lower_upper(test_name,a0_c(1:na,1:na),a1_c(1:na,1:na),1,wbv_na_c%lbw,1, &
       wbv_na_c%ubw,t0,t1,tol2,error)
  call c_deallocate_ubt(ubt_na_c)
  call c_deallocate_wbv(wbv_na_c)

  na=4
  ubt_na_c=c_new_ubt(na,lbwmax,ubwmax)
  wbv_na_c=c_new_wbv(na,lbwmax,ubwmax)
  ul_c=ul0_c; vl_c=vl0_c
  uu_c=uu0_c; vu_c=vu0_c
  call c_assemble_general(a_c(1:na,1:na),ul_c(1:na,:),vl_c(:,1:na),uu_c(1:na,:),vu_c(:,1:na),&
       d_c(1:na))
  a0_c(1:na,1:na)=a_c(1:na,1:na)
  call general_to_ubt(a_c(1:na,1:na),ubt_na_c,tol1,error)
  call cpu_time(t0)
  call convert_ubt_to_wbv(ubt_na_c, wbv_na_c,error)
  call cpu_time(t1)
  call wbv_to_general(wbv_na_c,a1_c(1:na,1:na),error)
  test_name = "Complex UBT to WBV (n=4);"
  call c_output_result_lower_upper(test_name,a0_c(1:na,1:na),a1_c(1:na,1:na),2,wbv_na_c%lbw,2, &
       wbv_na_c%ubw,t0,t1,tol2,error)
  call c_deallocate_ubt(ubt_na_c)
  call c_deallocate_wbv(wbv_na_c)

  ! wbv to ubt
  ul_c=ul0_c; vl_c=vl0_c
  uu_c=uu0_c; vu_c=vu0_c
  call c_assemble_general(a_c,ul_c,vl_c,uu_c,vu_c,d_c)
  a0_c=a_c
  call general_to_wbv(a_c,wbv_c,tol,error)
  call cpu_time(t0)
  call convert_wbv_to_ubt(wbv_c, ubt_c, error)
  call cpu_time(t1)
  call ubt_to_general(ubt_c,a1_c,error)
  test_name = "Complex WBV to UBT;"
  call c_output_result_lower_upper(test_name,a0_c,a1_c,rmaxl,ubt_c%lbw, &
       rmaxu,ubt_c%ubw,t0,t1,tol2,error)

  na=1
  ubt_na_c=c_new_ubt(na,lbwmax,ubwmax)
  wbv_na_c=c_new_wbv(na,lbwmax,ubwmax)
  ul_c=ul0_c; vl_c=vl0_c
  uu_c=uu0_c; vu_c=vu0_c
  call c_assemble_general(a_c(1:na,1:na),ul_c(1:na,:),vl_c(:,1:na),uu_c(1:na,:),vu_c(:,1:na),d_c(1:na))
  a0_c(1:na,1:na)=a_c(1:na,1:na)
  call general_to_wbv(a_c(1:na,1:na),wbv_na_c,tol1,error)
  call cpu_time(t0)
  call convert_wbv_to_ubt(wbv_na_c, ubt_na_c,error)
  call cpu_time(t1)
  call ubt_to_general(ubt_na_c,a1_c(1:na,1:na),error)
  test_name = "Complex WBV to UBT (n=1);"
  call c_output_result_lower_upper(test_name,a0_c(1:na,1:na),a1_c(1:na,1:na),0,ubt_na_c%lbw,0, &
       ubt_na_c%ubw,t0,t1,tol2,error)
  call c_deallocate_ubt(ubt_na_c)
  call c_deallocate_wbv(wbv_na_c)

  na=2
  ubt_na_c=c_new_ubt(na,lbwmax,ubwmax)
  wbv_na_c=c_new_wbv(na,lbwmax,ubwmax)
  ul_c=ul0_c; vl_c=vl0_c
  uu_c=uu0_c; vu_c=vu0_c
  call c_assemble_general(a_c(1:na,1:na),ul_c(1:na,:),vl_c(:,1:na),uu_c(1:na,:),vu_c(:,1:na),d_c(1:na))
  a0_c(1:na,1:na)=a_c(1:na,1:na)
  call general_to_wbv(a_c(1:na,1:na),wbv_na_c,tol1,error)
  call cpu_time(t0)
  call convert_wbv_to_ubt(wbv_na_c, ubt_na_c,error)
  call cpu_time(t1)
  call ubt_to_general(ubt_na_c,a1_c(1:na,1:na),error)
  test_name = "Complex WBV to UBT (n=2);"
  call c_output_result_lower_upper(test_name,a0_c(1:na,1:na),a1_c(1:na,1:na),1,ubt_na_c%lbw,1, &
       ubt_na_c%ubw,t0,t1,tol2,error)
  call c_deallocate_ubt(ubt_na_c)
  call c_deallocate_wbv(wbv_na_c)

  na=3
  ubt_na_c=c_new_ubt(na,lbwmax,ubwmax)
  wbv_na_c=c_new_wbv(na,lbwmax,ubwmax)
  ul_c=ul0_c; vl_c=vl0_c
  uu_c=uu0_c; vu_c=vu0_c
  call c_assemble_general(a_c(1:na,1:na),ul_c(1:na,:),vl_c(:,1:na),uu_c(1:na,:),vu_c(:,1:na),d_c(1:na))
  a0_c(1:na,1:na)=a_c(1:na,1:na)
  call general_to_wbv(a_c(1:na,1:na),wbv_na_c,tol1,error)
  call cpu_time(t0)
  call convert_wbv_to_ubt(wbv_na_c, ubt_na_c,error)
  call cpu_time(t1)
  call ubt_to_general(ubt_na_c,a1_c(1:na,1:na),error)
  test_name = "Complex WBV to UBT (n=3);"
  call c_output_result_lower_upper(test_name,a0_c(1:na,1:na),a1_c(1:na,1:na),1,ubt_na_c%lbw,1, &
       ubt_na_c%ubw,t0,t1,tol2,error)
  call c_deallocate_ubt(ubt_na_c)
  call c_deallocate_wbv(wbv_na_c)

  na=4
  ubt_na_c=c_new_ubt(na,lbwmax,ubwmax)
  wbv_na_c=c_new_wbv(na,lbwmax,ubwmax)
  ul_c=ul0_c; vl_c=vl0_c
  uu_c=uu0_c; vu_c=vu0_c
  call c_assemble_general(a_c(1:na,1:na),ul_c(1:na,:),vl_c(:,1:na),uu_c(1:na,:),vu_c(:,1:na),d_c(1:na))
  a0_c(1:na,1:na)=a_c(1:na,1:na)
  call general_to_wbv(a_c(1:na,1:na),wbv_na_c,tol1,error)
  call cpu_time(t0)
  call convert_wbv_to_ubt(wbv_na_c, ubt_na_c,error)
  call cpu_time(t1)
  call ubt_to_general(ubt_na_c,a1_c(1:na,1:na),error)
  test_name = "Complex WBV to UBT (n=4);"
  call c_output_result_lower_upper(test_name,a0_c(1:na,1:na),a1_c(1:na,1:na),2,ubt_na_c%lbw,2, &
       ubt_na_c%ubw,t0,t1,tol2,error)
  call c_deallocate_ubt(ubt_na_c)
  call c_deallocate_wbv(wbv_na_c)

end program test_convert_ubt_and_wbv
