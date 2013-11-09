program test_convert
  use upper_decomp
  use utility
  use assemble
  use upper_convert
  use band
  implicit none
  integer(kind=int32), parameter :: n=500, rmax=4, ubwmax=rmax, lbw=3, mb=ubwmax+lbw+5, nb=mb
  real(kind=dp), parameter :: tol=1e-15
  real(kind=dp) :: t1, t2
  integer(kind=int32) :: j, k, error
  !
  real(kind=dp), dimension(n,n) :: a, a0, a1
  real(kind=dp), dimension(mb,n) :: b_ub
  real(kind=dp), dimension(n,mb) :: b_bv
  real(kind=dp), dimension(n,rmax) :: u, u0
  real(kind=dp), dimension(rmax,n) :: v, v0
  real(kind=dp), dimension(n) :: d
  real(kind=dp), dimension(mb-lbw-1, n) :: cs_ub, ss_ub
  real(kind=dp), dimension(n,mb-lbw-1) :: cs_bv, ss_bv
  integer(kind=int32) :: ubw
  integer(kind=int32), dimension(mb-lbw-1, n) :: j1s, j2s
  integer(kind=int32), dimension(n,mb-lbw-1) :: k1s, k2s
  integer(kind=int32), dimension(n) :: numrots_ub, numrots_bv
  !
  complex(kind=dp), dimension(n,n) :: a_c, a0_c, a1_c
  complex(kind=dp), dimension(mb,n) :: b_ub_c
  complex(kind=dp), dimension(n,mb) :: b_bv_c
  complex(kind=dp), dimension(n,rmax) :: u_c, u0_c
  complex(kind=dp), dimension(rmax,n) :: v_c, v0_c
  complex(kind=dp), dimension(n) :: d_c
  complex(kind=dp), dimension(mb-lbw-1, n) :: cs_ub_c, ss_ub_c
  complex(kind=dp), dimension(n,mb-lbw-1) :: cs_bv_c, ss_bv_c
  integer(kind=int32) :: ubw_c
  integer(kind=int32), dimension(mb-lbw-1, n) :: j1s_c, j2s_c
  integer(kind=int32), dimension(n,mb-lbw-1) :: k1s_c, k2s_c
  integer(kind=int32), dimension(n) :: numrots_ub_c, numrots_bv_c

  ! real
  call random_seed
  call random_number(u)
  call random_number(v)
  call random_number(d)
  u0=u; v0=v
  ! test one
  a=matmul(u,v)
  do k=1,n-lbw-1
     do j=k+lbw+1,n
        a(j,k)=0.0_dp
     end do
  end do
  do j=1,n
     a(j,j)=d(j)
  end do
  a0=a
  call upper_general_to_upper_ub(a,n, b_ub, mb, lbw, ubw, numrots_ub, j1s, j2s,  &
       cs_ub, ss_ub, tol, error)
  if (error == 1) then
     print *, "orthogonalization error in real UB"
  else
     call cpu_time(t1)
     call d_upper_ub_to_bv(b_ub, mb, n, lbw, ubw, numrots_ub, j1s, j2s, cs_ub, ss_ub, &
          b_bv, mb, numrots_bv, k1s, k2s, cs_bv, ss_bv, error)
     call cpu_time(t2)
     if (error > 0) then
        print *, "Error in upper_ub_to_bv:", error
     else
        call form_upper_bv(a1, n, b_bv, mb, lbw, ubw, numrots_bv, k1s, k2s, cs_bv, ss_bv)
        write (*, "('Real UB to BV;                     ', 'Time: ',ES8.2,', ubw: ',I3,', error: ',ES8.2)") &
             t2-t1, ubw, maxabs(a1-a0)
     end if
  end if



  ! bv to ub
  a=a0
  call upper_general_to_upper_bv(a,n, b_bv, mb, lbw, ubw, numrots_bv, k1s, k2s,  &
       cs_bv, ss_bv, tol, error)
  if (error == 1) then
     print *, "orthogonalization error in real BV"
  else
     call cpu_time(t1)
     call d_upper_bv_to_ub(b_bv, mb, n, lbw, ubw, numrots_bv, k1s, k2s, cs_bv, ss_bv, &
          b_ub, mb, numrots_ub, j1s, j2s, cs_ub, ss_ub, error)
     call cpu_time(t2)
     if (error > 0) then
        print *, "Error in upper_bv_to_ub:", error
     else
        call form_upper_ub(a1, n, b_ub, mb, lbw, ubw, numrots_ub, j1s, j2s, cs_ub, ss_ub)
        write (*, "('Real BV to UB;                     ', 'Time: ',ES8.2,', ubw: ',I3,', error: ',ES8.2)") &
             t2-t1, ubw, maxabs(a1-a0)/maxabs(a0)
     end if
  end if

! complex

  call random_complex(u_c)
  call random_complex(v_c)
  call random_complex(d_c)
  u0_c=u_c; v0_c=v_c
  a_c=matmul(u_c,v_c)
  do k=1,n-lbw-1
     do j=k+lbw+1,n
        a_c(j,k)=(0.0_dp, 0.0_dp)
     end do
  end do
  do j=1,n
     a_c(j,j)=d_c(j)
  end do
  ! ub to bv
  a0_c=a_c
  call upper_general_to_upper_ub(a_c,n, b_ub_c, mb, lbw, ubw_c, numrots_ub_c, j1s_c, j2s_c,  &
       cs_ub_c, ss_ub_c, tol, error)
  if (error == 1) then
     print *, "orthogonalization error in complex UB"
  else
     call cpu_time(t1)
     call c_upper_ub_to_bv(b_ub_c, mb, n, lbw, ubw_c, numrots_ub_c, j1s_c, j2s_c, cs_ub_c, ss_ub_c, &
          b_bv_c, mb, numrots_bv_c, k1s_c, k2s_c, cs_bv_c, ss_bv_c, error)
     call cpu_time(t2)
     if (error > 0) then
        print *, "Error in upper_ub_to_bv:", error
     else
        call form_upper_bv(a1_c, n, b_bv_c, mb, lbw, ubw_c, numrots_bv_c, k1s_c, k2s_c, cs_bv_c, ss_bv_c)
        write (*, "('Complex UB to BV;                     ', 'Time: ',ES8.2,', ubw: ',I3,', error: ',ES8.2)") &
             t2-t1, ubw_c, maxabs(a1_c-a0_c)/maxabs(a0_c)
     end if
  end if
  ! bv to ub
  a_c=a0_c
  call upper_general_to_upper_bv(a_c,n, b_bv_c, mb, lbw, ubw_c, numrots_bv_c, k1s_c, k2s_c,  &
       cs_bv_c, ss_bv_c, tol, error)
  if (error == 1) then
     print *, "orthogonalization error in complex BV"
  else
     call cpu_time(t1)
     call c_upper_bv_to_ub(b_bv_c, mb, n, lbw, ubw_c, numrots_bv_c, k1s_c, k2s_c, cs_bv_c, ss_bv_c, &
          b_ub_c, mb, numrots_ub_c, j1s_c, j2s_c, cs_ub_c, ss_ub_c, error)
     call cpu_time(t2)
     if (error > 0) then
        print *, "Error in upper_bv_to_ub:", error
     else
        call form_upper_ub(a1_c, n, b_ub_c, mb, lbw, ubw_c, numrots_ub_c, j1s_c, j2s_c, cs_ub_c, ss_ub_c)
        write (*, "('Complex BV to UB;                     ', 'Time: ',ES8.2,', ubw: ',I3,', error: ',ES8.2)") &
             t2-t1, ubw_c, maxabs(a1_c-a0_c)/maxabs(a0_c)
     end if
  end if

end program test_convert
