program test_decomp
  use upper_decomp
  use utility
  use assemble
  implicit none
  integer(kind=int32), parameter :: n=1000, rmax=5, ubwmax=rmax, lbw=1, lbwmax=lbw, mb=ubwmax+lbwmax+1, nb=mb
  real(kind=dp), parameter :: tol=1e-14
  real(kind=dp) :: t1, t2
  integer(kind=int32) :: j, k, error
  !
  real(kind=dp), dimension(n,n) :: a, a0, a1
  real(kind=dp), dimension(mb,n) :: b_ub
  real(kind=dp), dimension(n,mb) :: b_bv
  real(kind=dp), dimension(n,rmax) :: u, u0
  real(kind=dp), dimension(rmax,n) :: v, v0
  real(kind=dp), dimension(n) :: d
  integer(kind=int32) :: ubw
  integer(kind=int32), dimension(ubwmax, n) :: j1s, j2s
  integer(kind=int32), dimension(n,ubwmax) :: k1s, k2s
  integer(kind=int32), dimension(n) :: numrots
  real(kind=dp), dimension(ubwmax, n) :: cs_ub, ss_ub
  real(kind=dp), dimension(n,ubwmax) :: cs_bv, ss_bv
  !
  complex(kind=dp), dimension(n,n) :: ac, a0c, a1c
  complex(kind=dp), dimension(mb,n) :: bc_ub
  complex(kind=dp), dimension(n,mb) :: bc_bv
  complex(kind=dp), dimension(n,rmax) :: uc, u0c
  complex(kind=dp), dimension(rmax,n) :: vc, v0c
  complex(kind=dp), dimension(n) :: dc
  integer(kind=int32) :: ubwc
  integer(kind=int32), dimension(ubwmax, n) :: j1sc, j2sc
  integer(kind=int32), dimension(n,ubwmax) :: k1sc, k2sc
  integer(kind=int32), dimension(n) :: numrotsc
  complex(kind=dp), dimension(ubwmax, n) :: csc_ub, ssc_ub
  complex(kind=dp), dimension(n, ubwmax) :: csc_bv, ssc_bv

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
  call cpu_time(t1)
  call f_upper_to_ub(a,n, b_ub, lbw, ubw, lbwmax, ubwmax, numrots, j1s, j2s, cs_ub, ss_ub, tol, error)
  call cpu_time(t2)
  if (error == 1) then
     print *, "orthogonalization error in real UB"
  else
     call f_ub_to_upper(b_ub, n, lbw, ubw, lbwmax, ubwmax, numrots, j1s, j2s, cs_ub, ss_ub, a1)
     write (*,"('Real UB;                     ', 'Time: ',ES8.2,', ubw: ',I3,', error: ',ES8.2)") &
          t2-t1, ubw, maxabs(a1-a0)
  end if
  !
  ! test 2: square termination.
  !
  u(1:n-rmax-1,rmax)=0.0_dp
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
  call cpu_time(t1)
  call f_upper_to_ub(a, n, b_ub, lbw, ubw, lbwmax, ubwmax, numrots, j1s, j2s, cs_ub, ss_ub, tol, error)
  call cpu_time(t2)
  if (error == 1) then
     print *, "orthogonalization error in real square termination UB"
  else
     call f_ub_to_upper(b_ub, n, lbw, ubw, lbwmax, ubwmax, numrots, j1s, j2s, cs_ub, ss_ub, a1)
     write (*,"('Real Sq. term. UB;           ', 'Time: ',ES8.2,', ubw: ',I3,', error: ',ES8.2)") &
          t2-t1, ubw, maxabs(a1-a0)
  end if
  ! BV
  u=u0; v=v0
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
  call cpu_time(t1)
  call f_upper_to_bv(a, n, b_bv, lbw, ubw, lbwmax, ubwmax, numrots, k1s, k2s, cs_bv, ss_bv, tol, error)
  call cpu_time(t2)
  if (error == 1) then
     print *, "orthogonalization error in real BV"
  else
     call f_bv_to_upper(b_bv, n, lbw, ubw, lbwmax, ubwmax, numrots, k1s, k2s, cs_bv, ss_bv, a1)
     write (*,"('Real BV;                     ','Time: ',ES8.2,', ubw: ',I3,', error: ',ES8.2)") &
          t2-t1, ubw, maxabs(a1-a0)
  end if
  ! square termination BV
  u=u0; v=v0
  v(rmax,rmax+2:n)=0.0_dp
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
  call cpu_time(t1)
  call f_upper_to_bv(a,n, b_bv, lbw, ubw, lbwmax, ubwmax, numrots, k1s, k2s, cs_bv, ss_bv, tol, error)
  call cpu_time(t2)
  if (error == 1) then
     print *, "orthogonalization error in real square termination BV"
  else
     call f_bv_to_upper(b_bv, n, lbw, ubw, lbwmax, ubwmax, numrots, k1s, k2s, cs_bv, ss_bv, a1)
     write (*,"('Real Sq. Term BV;            ','Time: ',ES8.2,', ubw: ',I3,', error: ',ES8.2)") &
          t2-t1, ubw, maxabs(a1-a0)
  end if
  !
  ! Complex UB test
  !
  ac=matmul(uc,vc)
  do k=1,n-lbw-1
     do j=k+lbw+1,n
        ac(j,k)=(0.0_dp, 0.0_dp)
     end do
  end do
  do j=1,n
     ac(j,j)=dc(j)
  end do
  a0c=ac
  call cpu_time(t1)
  call f_upper_to_ub(ac,n, bc_ub, lbw, ubwc, lbwmax, ubwmax, numrotsc, j1sc, j2sc, &
       csc_ub, ssc_ub, tol, error)
  call cpu_time(t2)
  if (error == 1) then
     print *, "orthogonalization error in complex UB"
  else
     call f_ub_to_upper(bc_ub, n, lbw, ubwc, lbwmax, ubwmax, numrotsc, j1sc, j2sc, csc_ub, ssc_ub, a1c)
     write (*,"('Complex UB;                  ','Time: ',ES8.2,', ubw: ',I3,', error: ',ES8.2)") &
          t2-t1, ubwc, maxabs(a1c-a0c)
  end if
  !
  ! Complex UB square termination test
  !
  uc=u0c; vc=v0c
  uc(1:n-rmax-1,rmax)=(0.0_dp, 0.0_dp)
  ac=matmul(uc,vc)
  do k=1,n-lbw-1
     do j=k+lbw+1,n
        ac(j,k)=(0.0_dp, 0.0_dp)
     end do
  end do
  do j=1,n
     ac(j,j)=dc(j)
  end do
  a0c=ac
  call cpu_time(t1)
  call f_upper_to_ub(ac,n, bc_ub, lbw, ubwc, lbwmax, ubwmax, numrotsc, j1sc, j2sc, &
       csc_ub, ssc_ub, tol, error)
  call cpu_time(t2)
  if (error == 1) then
     print *, "orthogonalization error in complex UB"
  else
     call f_ub_to_upper(bc_ub, n, lbw, ubwc, lbwmax, ubwmax, numrotsc, j1sc, j2sc, csc_ub, ssc_ub, a1c)
     write (*,"('Complex square term. UB;     ','Time: ',ES8.2,', ubw: ',I3,', error: ',ES8.2)") &
          t2-t1, ubwc, maxabs(a1c-a0c)
  end if
  !
  ! Complex BV
  !
  uc=u0c; vc=v0c
  uc(1:n-rmax-1,rmax)=(0.0_dp, 0.0_dp)
  ac=matmul(uc,vc)
  do k=1,n-lbw-1
     do j=k+lbw+1,n
        ac(j,k)=(0.0_dp, 0.0_dp)
     end do
  end do
  do j=1,n
     ac(j,j)=dc(j)
  end do
  a0c=ac
  call cpu_time(t1)
  call f_upper_to_bv(ac,n, bc_bv, lbw, ubwc, lbwmax, ubwmax, numrots, k1sc, k2sc, csc_bv, ssc_bv, &
       tol, error)
  call cpu_time(t2)
  if (error == 1) then
     print *, "orthogonalization error in complex BV"
  else
     call f_bv_to_upper(bc_bv, n, lbw, ubwc, lbwmax, ubwmax, numrots, k1sc, k2sc, csc_bv, ssc_bv, a1c)
     write (*,"('Real BV;                     ','Time: ',ES8.2,', ubw: ',I3,', error: ',ES8.2)") &
          t2-t1, ubw, maxabs(a1c-a0c)
  end if
end program test_decomp
