program test_decomp
  use decomp
  use utility
  implicit none
  integer(kind=int32), parameter :: n=100, rmax=2, ubwmax=100, lbw=1, mb=ubwmax+lbw+1
  real(kind=dp), parameter :: tol=1e-18
  real(kind=dp), dimension(n,n) :: a, a0, a1
  real(kind=dp), dimension(mb,n) :: b
  real(kind=dp), dimension(n,rmax) :: u
  real(kind=dp), dimension(rmax,n) :: v
  real(kind=dp), dimension(n) :: d
  integer(kind=int32) :: j, k, error, ubw
  integer(kind=int32), dimension(ubwmax, n) :: j1s, j2s
  integer(kind=int32), dimension(n) :: numrots
  real(kind=dp), dimension(ubwmax, n) :: cs, ss
  real(kind=dp) :: t1, t2

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
!  print *
  ! call print_matrix(a0-a1)
  ! print *
  ! call print_matrix(a1)
  ! call print_matrix(a)
  ! print *, "numrots"
  ! print *, numrots
  ! print *, "j1s"
  ! call print_matrix(j1s)
  ! print *, "j2s"
  ! call print_matrix(j2s)
  ! print *, "cs"
  ! call print_matrix(cs)
  ! print *, "ss"
  ! call print_matrix(ss)
  ! print *, "b"
  ! call print_matrix(b)
  ! print *
  ! call print_matrix(a)
end program test_decomp
