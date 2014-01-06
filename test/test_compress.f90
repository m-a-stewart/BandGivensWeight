program test_compress
  use general_ub
  use general_bv
  use utility
  use assemble
  use convert_ub
  use convert_bv
  use compress_ub
  use band_types
  use nested_types
  implicit none
  integer(kind=int32), parameter :: n=10000, rmax=4, ubwmax=rmax+1, lbw=3, lbwmax=3
  real(kind=dp), parameter :: tol=1e-14, tol1=1e-14, tol2=1e-10
  real(kind=dp) :: t1, t2
  integer(kind=int32) :: error
  character(len=*), parameter :: fmt="(A30, 'Time: ',ES8.2,', ubw: ',I3,', error: ',ES8.2)"
  character(len=30) :: test_name
  !
  real(kind=dp), dimension(n,n) :: a, a0, a1
  real(kind=dp), dimension(n,rmax) :: u, u0
  real(kind=dp), dimension(rmax,n) :: v, v0
  real(kind=dp), dimension(n) :: d
  complex(kind=dp), dimension(n,n) :: a_c, a0_c, a1_c
  complex(kind=dp), dimension(n,rmax) :: u_c, u0_c
  complex(kind=dp), dimension(rmax,n) :: v_c, v0_c
  complex(kind=dp), dimension(n) :: d_c

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

  ! real
  call random_seed
  call random_number(u)
  call random_number(v)
  call random_number(d)
  u(:,rmax)=0.1*tol2*u(:,rmax)
  u0=u; v0=v
  ! test one
  call d_assemble_a(a,u,v,d,lbw)
  a0=a
  call upper_to_ub(a,ub_d,tol1,error)
  print *, ub_d%ubw
  if (error == 1) then
     print *, "orthogonalization error in real UB"
  else
     call cpu_time(t1)
     call compress_ub_to_bv(ub_d, bv_d,tol2,tol2,error)
     call cpu_time(t2)
     if (error > 0) then
        print *, "Error in compress_ub_to_bv:", error
     else
        call bv_to_upper(bv_d,a1,error)
        test_name = "Real UB to BV;"
        write (*,fmt) test_name, t2-t1, bv_d%ubw, maxabs(a1-a0)
     end if
  end if
  !
  ! complex
  call random_complex(u_c)
  call random_complex(v_c)
  call random_complex(d_c)
  u_c(:,rmax)=0.1*tol2*u_c(:,rmax)
  u0_c=u_c; v0_c=v_c
  call c_assemble_a(a_c, u_c, v_c, d_c, lbw)
  a0_c=a_c
  call upper_to_ub(a_c, ub_c, tol, error)
  if (error == 1) then
     print *, "orthogonalization error in complex UB"
  else
     call cpu_time(t1)
     call compress_ub_to_bv(ub_c, bv_c, tol2, tol2, error)
     call cpu_time(t2)
     if (error > 0) then
        print *, "Error in convert_ub_to_bv:", error
     else
        call bv_to_upper(bv_c, a1_c, error)
        test_name="Complex UB to BV;"
        write (*,fmt) test_name, t2-t1, bv_c%ubw, maxabs(a1_c-a0_c)/maxabs(a0_c)
     end if
  end if

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

end program test_compress
