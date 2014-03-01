module gs
  use misc
  implicit none
  real(kind=dp), parameter :: eta=1.414_dp
  integer(kind=int32), parameter :: orthmaxits=5

  interface extend_gs_rows
     module procedure d_extend_gs_rows, c_extend_gs_rows
  end interface extend_gs_rows

  interface extend_gs_columns
     module procedure d_extend_gs_columns, c_extend_gs_columns
  end interface extend_gs_columns

  type(routine_info), parameter :: info_d_extend_gs_rows=routine_info(id_d_extend_gs_rows, &
       'd_extend_gs_rows', &
       [ character(len=error_message_length) :: 'MGS Orthogonalization Error' ])

  type(routine_info), parameter :: info_c_extend_gs_rows=routine_info(id_c_extend_gs_rows, &
       'c_extend_gs_rows', &
       [ character(len=error_message_length) :: 'MGS Orthogonalization Error' ])

  type(routine_info), parameter :: info_d_extend_gs_columns=routine_info(id_d_extend_gs_columns, &
       'd_extend_gs_columns', &
       [ character(len=error_message_length) :: 'MGS Orthogonalization Error' ])

  type(routine_info), parameter :: info_c_extend_gs_columns=routine_info(id_c_extend_gs_columns, &
       'c_extend_gs_columns', &
       [ character(len=error_message_length) :: 'MGS Orthogonalization Error' ])


contains

  ! Error 1: failure to orthogonalize
  ! Extend a Gram-Schmidt LQ decomposition

  subroutine d_extend_gs_rows(q1, l, rho, q2, error)
    ! orthogonalize q2 against q1.  Return error = 1 if
    ! too many iterations are required.
    real(kind=dp), dimension (:,:), intent(in) :: q1
    real(kind=dp), intent(out) :: rho
    real(kind=dp), dimension(:), intent(out) :: l
    real(kind=dp), dimension(:), intent(inout) :: q2
    type(error_info), intent(out) :: error
    !
    real(kind=dp), dimension(size(l)) :: x
    real(kind=dp) :: nrm1, nrm2
    integer(kind=int32) :: k
    !
    nrm1=norm2(q2)
    l = matmul(q1,q2)
    q2=q2-matmul(l,q1)
    nrm2=norm2(q2)
    rho=nrm2
    ! reorthogonalize as needed
    call clear_error(error)
    k=0
    do while (nrm1 > eta * nrm2 .and. nrm2 > 0 .and. k < orthmaxits)
       k=k+1
       x=matmul(q1,q2)
       l=l+x
       nrm1=norm2(q2)
       q2=q2-matmul(x,q1)
       nrm2=norm2(q2)
       rho=nrm2
    end do
    do while (nrm2 == 0 .and. k < orthmaxits)
       ! got the zero vector.  Orthogonalize something random.
       k=0
       call random_number(q2)
       nrm2 = 1.0_dp
       nrm1=2.0_dp*eta
       do while (nrm1 > eta*nrm2 .and. nrm2 > 0.0_dp .and. k < orthmaxits)
          k=k+1
          nrm1=norm2(q2)
          x=matmul(q1,q2)
          q2=q2-matmul(x,q1)
          nrm2=norm2(q2)
       end do
    end do
    q2=q2/nrm2
    if (k >= orthmaxits) then
       call set_error(error, 1, id_d_extend_gs_rows)
    end if
  end subroutine d_extend_gs_rows

  ! Extend a Gram-Schmidt LQ decomposition
  subroutine c_extend_gs_rows(q1, l, rho, q2, error)
    ! orthogonalize q2 against q1.  Return error = 1 if
    ! too many iterations are required.
    complex(kind=dp), dimension (:,:), intent(in) :: q1
    complex(kind=dp), intent(out) :: rho
    complex(kind=dp), dimension(:), intent(out) :: l
    complex(kind=dp), dimension(:), intent(inout) :: q2
    type(error_info), intent(out) :: error
    !
    complex(kind=dp), dimension(size(l)) :: x
    real(kind=dp) :: nrm1, nrm2
    integer(kind=int32) :: k
    !
    nrm1=norm2(q2)
    l = matmul(q1,conjg(q2))
    l=conjg(l)
    q2=q2-matmul(l,q1)
    nrm2=norm2(q2)
    rho=nrm2
    ! reorthogonalize as needed
    call clear_error(error)
    k=0
    do while (nrm1 > eta * nrm2 .and. nrm2 > 0 .and. k < orthmaxits)
       k=k+1
       x=matmul(q1,conjg(q2))
       x=conjg(x)
       l=l+x
       nrm1=norm2(q2)
       q2=q2-matmul(x,q1)
       nrm2=norm2(q2)
       rho=nrm2
    end do
    do while (nrm2 == 0 .and. k < orthmaxits)
       ! got the zero vector.  Orthogonalize something random.
       k=0
       call random_complex(q2)
       nrm2 = 1.0_dp
       nrm1=2.0_dp*eta
       do while (nrm1 > eta*nrm2 .and. nrm2 > 0.0_dp .and. k < orthmaxits)
          k=k+1
          nrm1=norm2(q2)
          x=matmul(q1,conjg(q2))
          x=conjg(x)
          q2=q2-matmul(x,q1)
          nrm2=norm2(q2)
       end do
    end do
    q2=q2/nrm2
    if (k >= orthmaxits) then
       call set_error(error, 1, id_c_extend_gs_rows)
    end if
  end subroutine c_extend_gs_rows

  ! Extend a Gram-Schmidt QL decomposition
  subroutine d_extend_gs_columns(q1, l, rho, q2, error)
    ! orthogonalize q2 against q1.  Return error = 1 if
    ! too many iterations are required.
    real(kind=dp), dimension (:,:), intent(in) :: q1
    real(kind=dp), intent(out) :: rho
    real(kind=dp), dimension(:), intent(out) :: l
    real(kind=dp), dimension(:), intent(inout) :: q2
    type(error_info), intent(out) :: error
    !
    real(kind=dp), dimension(size(l)) :: x
    real(kind=dp) :: nrm1, nrm2
    integer(kind=int32) :: k
    !
    nrm1=norm2(q2)
    l = matmul(q2,q1)
    q2=q2-matmul(q1,l)
    nrm2=norm2(q2)
    rho=nrm2
    ! reorthogonalize as needed
    call clear_error(error)
    k=0
    do while (nrm1 > eta * nrm2 .and. nrm2 > 0 .and. k < orthmaxits)
       k=k+1
       x=matmul(q2,q1)
       l=l+x
       nrm1=norm2(q2)
       q2=q2-matmul(q1,x)
       nrm2=norm2(q2)
       rho=nrm2
    end do
    do while (nrm2 == 0 .and. k < orthmaxits)
       ! got the zero vector.  Orthogonalize something random.
       k=0
       call random_number(q2)
       nrm2 = 1.0_dp
       nrm1=2.0_dp*eta
       do while (nrm1 > eta*nrm2 .and. nrm2 > 0.0_dp .and. k < orthmaxits)
          k=k+1
          nrm1=norm2(q2)
          x=matmul(q2,q1)
          q2=q2-matmul(q1,x)
          nrm2=norm2(q2)
       end do
    end do
    q2=q2/nrm2
    if (k >= orthmaxits) then
       call set_error(error, 1, id_d_extend_gs_columns)
    end if
  end subroutine d_extend_gs_columns

  ! Extend a Gram-Schmidt QL decomposition
  subroutine c_extend_gs_columns(q1, l, rho, q2, error)
    ! orthogonalize q2 against q1.  Return error = 1 if
    ! too many iterations are required.
    complex(kind=dp), dimension (:,:), intent(in) :: q1
    complex(kind=dp), intent(out) :: rho
    complex(kind=dp), dimension(:), intent(out) :: l
    complex(kind=dp), dimension(:), intent(inout) :: q2
    type(error_info), intent(out) :: error
    !
    complex(kind=dp), dimension(size(l)) :: x
    real(kind=dp) :: nrm1, nrm2
    integer(kind=int32) :: k
    !
    nrm1=norm2(q2)
    l = matmul(conjg(q2),q1)
    l=conjg(l)
    q2=q2-matmul(q1,l)
    nrm2=norm2(q2)
    rho=nrm2
    ! reorthogonalize as needed
    call clear_error(error)
    k=0
    do while (nrm1 > eta * nrm2 .and. nrm2 > 0 .and. k < orthmaxits)
       k=k+1
       x=matmul(conjg(q2),q1)
       x=conjg(x)
       l=l+x
       nrm1=norm2(q2)
       q2=q2-matmul(q1,x)
       nrm2=norm2(q2)
       rho=nrm2
    end do
    do while (nrm2 == 0 .and. k < orthmaxits)
       ! got the zero vector.  Orthogonalize something random.
       k=0
       call random_complex(q2)
       nrm2 = 1.0_dp
       nrm1=2.0_dp*eta
       do while (nrm1 > eta*nrm2 .and. nrm2 > 0.0_dp .and. k < orthmaxits)
          k=k+1
          nrm1=norm2(q2)
          x=matmul(conjg(q2),q1)
          x=conjg(x)
          q2=q2-matmul(q1,x)
          nrm2=norm2(q2)
       end do
    end do
    q2=q2/nrm2
    if (k >= orthmaxits) then
       call set_error(error, 1, id_c_extend_gs_columns)
    end if
  end subroutine c_extend_gs_columns


end module gs
