program test_decomp
  use general_ub
  use general_bv
  use utility
  use assemble
  use nested_types
  use test_data
  implicit none
  real(kind=dp) :: t1, t2
  integer(kind=int32) :: error
  character(len=*), parameter :: fmt="(A30, 'Time: ',ES8.2,', ubw: ',I3,', error: ',ES8.2)"
  character(len=30) :: test_name
  !
  real(kind=dp), dimension(n,n) :: a, a0, a1
  real(kind=dp), dimension(n,rmax) :: u, u0
  real(kind=dp), dimension(rmax,n) :: v, v0
  real(kind=dp), dimension(n) :: d
  complex(kind=dp), dimension(n,n) :: ac, a0c, a1c
  complex(kind=dp), dimension(n,rmax) :: uc, u0c
  complex(kind=dp), dimension(rmax,n) :: vc, v0c
  complex(kind=dp), dimension(n) :: dc
  type(d_ub) :: ub_d
  type(c_ub) :: ub_c
  type(d_bv) :: bv_d
  type(c_bv) :: bv_c
  ub_d=d_new_ub(n,lbwmax,ubwmax)
  ub_c=c_new_ub(n,lbwmax,ubwmax)
  bv_d=d_new_bv(n,lbwmax,ubwmax)
  bv_c=c_new_bv(n,lbwmax,ubwmax)
  ub_d%lbw=lbw;   ub_c%lbw=lbw
  bv_d%lbw=lbw;   bv_c%lbw=lbw

  call random_seed
  call random_number(u)
  call random_number(v)
  call random_number(d)
  u0=u; v0=v
  call random_complex(uc)
  call random_complex(vc)
  call random_complex(dc)
  u0c=uc; v0c=vc

  ! test one

  call d_assemble_a(a,u,v,d,lbw)
  a0=a
  call cpu_time(t1)
  call upper_to_ub(a,ub_d,tol,error)
  call cpu_time(t2)
  if (error == 1) then
     print *, "orthogonalization error in real UB"
  else
     call ub_to_upper(ub_d,a1,error)
     test_name="Real UB;"
     write(*,fmt) test_name, t2-t1, ub_d%ubw, maxabs(a1-a0)
  end if

  ! test 2: square termination.

  u=u0; v=v0;
  u(1:n-rmax-1,rmax)=0.0_dp
  call d_assemble_a(a,u,v,d,lbw)
  a0=a
  call cpu_time(t1)
  call upper_to_ub(a,ub_d,tol,error)
  call cpu_time(t2)
  if (error == 1) then
     print *, "orthogonalization error in real square termination UB"
  else
     call ub_to_upper(ub_d,a1,error)
     test_name="Real Sq. term. UB;"
     write (*,fmt) test_name, t2-t1, ub_d%ubw, maxabs(a1-a0)
  end if
  ! BV
  u=u0; v=v0
  call d_assemble_a(a,u,v,d,lbw)
  a0=a
  call cpu_time(t1)
  call upper_to_bv(a,bv_d,tol,error)
  call cpu_time(t2)
  if (error == 1) then
     print *, "orthogonalization error in real BV"
  else
     call bv_to_upper(bv_d,a1,error)
     test_name="Real BV;"
     write(*,fmt) test_name, t2-t1, bv_d%ubw, maxabs(a1-a0)
  end if
  ! square termination BV
  u=u0; v=v0
  v(rmax,rmax+2:n)=0.0_dp
  call d_assemble_a(a,u,v,d,lbw)
  a0=a
  call cpu_time(t1)
  call upper_to_bv(a,bv_d,tol,error)
  call cpu_time(t2)
  if (error == 1) then
     print *, "orthogonalization error in real square termination BV"
  else
     call bv_to_upper(bv_d,a1,error)
     test_name="Real Sq. Term BV;"
     write (*,fmt) test_name, t2-t1, bv_d%ubw, maxabs(a1-a0)
  end if
  !
  ! Complex UB test
  !
  uc=u0c; vc=v0c
  call c_assemble_a(ac,uc,vc,dc,lbw)
  a0c=ac
  call cpu_time(t1)
  call upper_to_ub(ac,ub_c,tol,error)
  call cpu_time(t2)
  if (error == 1) then
     print *, "orthogonalization error in complex UB"
  else
     call ub_to_upper(ub_c,a1c,error)
     test_name="Complex UB;"
     write (*,fmt) test_name, t2-t1, ub_c%ubw, maxabs(a1c-a0c)
  end if
  !
  ! Complex UB square termination test
  !
  uc=u0c; vc=v0c
  uc(1:n-rmax-1,rmax)=(0.0_dp, 0.0_dp)
  call c_assemble_a(ac,uc,vc,dc,lbw)
  a0c=ac
  call cpu_time(t1)
  call upper_to_ub(ac,ub_c,tol,error)
  call cpu_time(t2)
  if (error == 1) then
     print *, "orthogonalization error in complex UB"
  else
     call ub_to_upper(ub_c,a1c,error)
     test_name="Complex square term. UB;"
     write (*,fmt) test_name, t2-t1, ub_c%ubw, maxabs(a1c-a0c)
  end if
  !
  ! Complex BV
  !
  uc=u0c; vc=v0c
  uc(1:n-rmax-1,rmax)=(0.0_dp, 0.0_dp)
  call c_assemble_a(ac,uc,vc,dc,lbw)
  a0c=ac
  call cpu_time(t1)
  call upper_to_bv(ac,bv_c,tol,error)
  call cpu_time(t2)
  if (error == 1) then
     print *, "orthogonalization error in complex BV"
  else
     call bv_to_upper(bv_c,a1c,error)
     test_name="Complex BV;"
     write (*,fmt) test_name, t2-t1, bv_c%ubw, maxabs(a1c-a0c)
  end if
  print *
contains
  
  subroutine d_assemble_a(a,u,v,d,lbw)
    real(kind=dp), dimension(:,:), intent(out) :: a
    real(kind=dp), dimension(:,:), intent(in) :: u, v
    real(kind=dp), dimension(:), intent(in) :: d
    integer(kind=int32), intent(in) :: lbw

    integer(kind=int32) :: j,k,n
    n=size(a,1)
    a=matmul(u,v)
    do k=1,n-lbw-1
       do j=k+lbw+1,n
          a(j,k)=0.0_dp
       end do
    end do
    do j=1,n
       a(j,j)=d(j)
    end do
  end subroutine d_assemble_a

  subroutine c_assemble_a(a,u,v,d,lbw)
    complex(kind=dp), dimension(:,:), intent(out) :: a
    complex(kind=dp), dimension(:,:), intent(in) :: u, v
    complex(kind=dp), dimension(:), intent(in) :: d
    integer(kind=int32), intent(in) :: lbw

    integer(kind=int32) :: j,k,n
    n=size(a,1)
    a=matmul(u,v)
    do k=1,n-lbw-1
       do j=k+lbw+1,n
          a(j,k)=(0.0_dp, 0.0_dp)
       end do
    end do
    do j=1,n
       a(j,j)=d(j)
    end do
  end subroutine c_assemble_a



end program test_decomp
