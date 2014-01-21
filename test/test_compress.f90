program test_compress
  use general_ub
  use general_bv
  use utility
  use assemble
  use convert_ub
  use convert_bv
  use compress_ub
  use compress_bv
  use band_types
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
  if (error == 1) then
     print *, "orthogonalization error in real UB"
  else
     call cpu_time(t1)
     call compress_ub_to_bv_1(ub_d, bv_d,tol2,tol2,error)
     call cpu_time(t2)
     if (error > 0) then
        print *, "Error in compress_ub_to_bv:", error
     else
        call bv_to_upper(bv_d,a1,error)
        test_name = "Real UB to BV;"
        write (*,fmt) test_name, t2-t1, bv_d%ubw, maxabs(a1-a0)
     end if
  end if
  a=a0
  call upper_to_bv(a,bv_d,tol1,error)
  if (error == 1) then
     print *, "orthogonalization error in real BV"
  else
     call cpu_time(t1)
     call compress_bv_to_ub_1(bv_d, ub_d,tol2,tol2,error)
     call cpu_time(t2)
     if (error > 0) then
        print *, "Error in compress_bv_to_ub:", error
     else
        call ub_to_upper(ub_d,a1,error)
        test_name = "Real BV to UB;"
        write (*,fmt) test_name, t2-t1, ub_d%ubw, maxabs(a1-a0)
     end if
  end if
  !
  ! complex
  !
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
     call compress_ub_to_bv_1(ub_c, bv_c, tol2, tol2, error)
     call cpu_time(t2)
     if (error > 0) then
        print *, "Error in complex convert_ub_to_bv:", error
     else
        call bv_to_upper(bv_c, a1_c, error)
        test_name="Complex UB to BV;"
        write (*,fmt) test_name, t2-t1, bv_c%ubw, maxabs(a1_c-a0_c)/maxabs(a0_c)
     end if
  end if
  ! BV to UB
  a_c=a0_c
  call upper_to_bv(a_c, bv_c, tol, error)
  if (error == 1) then
     print *, "orthogonalization error in complex UB"
  else
     call cpu_time(t1)
     call compress_bv_to_ub_1(bv_c, ub_c, tol2, tol2, error)
     call cpu_time(t2)
     if (error > 0) then
        print *, "Error in complex convert_bv_to_ub:", error
     else
        call ub_to_upper(ub_c, a1_c, error)
        test_name="Complex BV to UB;"
        write (*,fmt) test_name, t2-t1, ub_c%ubw, maxabs(a1_c-a0_c)/maxabs(a0_c)
     end if
  end if

end program test_compress
