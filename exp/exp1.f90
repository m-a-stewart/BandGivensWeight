program exp1
  use mod_orb
  implicit none

  real(kind=dp) :: t0, t1
  integer(kind=int32) :: n, j, k
  type(error_info) :: error

  real(kind=dp), dimension(:,:), allocatable :: l_d
  real(kind=dp), dimension(:), allocatable :: u_d, v_d, res_d, res_d1
  real(kind=dp) :: sigma_d

  complex(kind=dp), dimension(:,:), allocatable :: l_c
  complex(kind=dp), dimension(:), allocatable :: u_c, v_c, res_c, res_c1
  real(kind=dp) :: sigma_c
  
  print *, "Running real numerical experiment."
  n=20
  allocate(l_d(n,n),u_d(n),v_d(n),res_d(n),res_d1(n))
  call random_matrix_to(l_d)
  do j=1,n
     do k=j+1,n
        l_d(j,k)=0.0_dp
     end do
  end do

  l_d(5,5)=2.3e-40_dp
  l_d(10,10)=2.3e-40_dp

  sigma_d=lower_min_sv(l_d,u_d,v_d,res_d,1.0e-14_dp,1000,error)

  print *, sigma_d
  print *, norm2(res_d)
  call lower_left_multiply(l_d,v_d,res_d)
  res_d=res_d-sigma_d*u_d
  call lower_tr_left_multiply(l_d,u_d,res_d1)
  res_d1=res_d1-sigma_d*v_d
  print *, norm2(res_d), norm2(res_d1)
  print *
  deallocate(l_d,u_d,v_d,res_d,res_d1)

  print *, "Running complex numerical experiment."
  n=20
  allocate(l_c(n,n),u_c(n),v_c(n),res_c(n),res_c1(n))
  call random_matrix_to(l_c)
  do j=1,n
     do k=j+1,n
        l_c(j,k)=(0.0_dp,0.0_dp)
     end do
  end do

  ! l_c(5,5)=(2.3e-40_dp,5e-40_dp)
  ! l_c(10,10)=(2.3e-40_dp,5e-40_dp)

  sigma_c=lower_min_sv(l_c,u_c,v_c,res_c,1.0e-14_dp,1000,error)

  print *, sigma_c
  print *, norm2(res_c)
  call lower_left_multiply(l_c,v_c,res_c)
  res_c=res_c-sigma_c*u_c
  call lower_tr_left_multiply(l_c,u_c,res_c1)
  res_c1=res_c1-sigma_c*v_c
  print *, norm2(res_c), norm2(res_c1)
  print *
  deallocate(l_c,u_c,v_c,res_c,res_c1)
  
  
end program exp1
