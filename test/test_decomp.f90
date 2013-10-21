program test_decomp
  use decomp
  use utility
  implicit none
  integer(kind=int32), parameter :: n=1000, rmax=2, ubwmax=100, lbw=3, mb=ubwmax+lbw+1
  real(kind=dp), parameter :: tol=1e-14
  real(kind=dp) :: t1, t2
  integer(kind=int32) :: j, k, error
  !
  real(kind=dp), dimension(n,n) :: a, a0, a1
  real(kind=dp), dimension(mb,n) :: b
  real(kind=dp), dimension(n,rmax) :: u
  real(kind=dp), dimension(rmax,n) :: v
  real(kind=dp), dimension(n) :: d
  integer(kind=int32) :: ubw
  integer(kind=int32), dimension(ubwmax, n) :: j1s, j2s
  integer(kind=int32), dimension(n) :: numrots
  real(kind=dp), dimension(ubwmax, n) :: cs, ss
  !
  complex(kind=dp), dimension(n,n) :: ac, a0c, a1c
  complex(kind=dp), dimension(mb,n) :: bc
  complex(kind=dp), dimension(n,rmax) :: uc
  complex(kind=dp), dimension(rmax,n) :: vc
  complex(kind=dp), dimension(n) :: dc
  integer(kind=int32) :: ubwc
  integer(kind=int32), dimension(ubwmax, n) :: j1sc, j2sc
  integer(kind=int32), dimension(n) :: numrotsc
  complex(kind=dp), dimension(ubwmax, n) :: csc, ssc

  call random_seed
  call random_number(u)
  call random_number(v)
  call random_number(d)
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
  print *, "extracting"
  call cpu_time(t1)
  call d_upper_general_to_upper_ub(a,n,ubwmax, tol, b, lbw, ubw, j1s, j2s, cs, ss, numrots, error)
  call cpu_time(t2)
  print *, t2-t1, "seconds"
  print *, "reconstructing"
  call d_form_upper_ub(a1, mb, n, b, lbw, ubw, j1s, j2s, cs, ss, numrots)
  print *, ubw  
  print *, "error: ", maxabs(a1-a0)
  ! Complex test
  call random_complex(uc)
  call random_complex(vc)
  call random_complex(dc)
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
  print *, "Complex extracting"
  call cpu_time(t1)
  call c_upper_general_to_upper_ub(ac,n,ubwmax, tol, bc, lbw, ubwc, j1sc, j2sc, csc, ssc, numrotsc, error)
  call cpu_time(t2)
  print *, t2-t1, "seconds"
  print *, "reconstructing"
  call c_form_upper_ub(a1c, mb, n, bc, lbw, ubwc, j1sc, j2sc, csc, ssc, numrotsc)
  print *, ubwc  
  print *, "error: ", maxabs(a1c-a0c)

end program test_decomp
