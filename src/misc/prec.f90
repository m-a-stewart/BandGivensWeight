module mod_prec

  integer, parameter :: dp = kind(0.0d0)
  integer, parameter :: int32 = selected_int_kind(8)
  real(kind=dp), parameter :: eps=epsilon(1.0_dp)/2
  real(kind=dp), parameter :: inf=1.0_dp/0.0_dp
  
  public

  interface is_number
     module procedure d_s_is_number, d_v_is_number, d_is_number, &
          z_s_is_number, z_v_is_number, z_is_number
  end interface is_number
contains

  type(logical) function d_s_is_number(a) result (numflag)
    real(kind=dp) :: a
    numflag=.true.
    if (a/=a .or. a==inf .or. a==-inf) numflag=.false.
  end function d_s_is_number

  type(logical) function d_v_is_number(a) result (numflag)
    real(kind=dp), dimension(:) :: a
    integer(kind=int32) :: j
    numflag=.true.
    do j=1,size(a)
       if (a(j)/=a(j) .or. a(j)==inf .or. a(j)==-inf) then
          numflag=.false.
          return
       end if
    end do
  end function d_v_is_number

  type(logical) function d_is_number(a) result (numflag)
    real(kind=dp), dimension(:,:) :: a
    integer(kind=int32) :: j,k
    numflag=.true.
    do j=1,size(a,1)
       do k=1,size(a,2)
          if (a(j,k)/=a(j,k) .or. a(j,k)==inf .or. a(j,k)==-inf) then
             numflag=.false.
             return
          end if
       end do
    end do
  end function d_is_number

  type(logical) function z_s_is_number(a) result (numflag)
    complex(kind=dp) :: a
    numflag=.true.
    if (a/=a .or. real(a)==inf .or. aimag(a)==inf .or. &
         real(a)==-inf .or. aimag(a)==-inf) numflag=.false.
  end function z_s_is_number

  type(logical) function z_v_is_number(a) result (numflag)
    complex(kind=dp), dimension(:) :: a
    integer(kind=int32) :: j
    numflag=.true.
    do j=1,size(a)
       if (a(j)/=a(j) .or. real(a(j))==inf .or. aimag(a(j))==inf .or. &
            real(a(j))==-inf .or. aimag(a(j))==-inf) then
          numflag=.false.
          return
       end if
    end do
  end function z_v_is_number

  type(logical) function z_is_number(a) result (numflag)
    complex(kind=dp), dimension(:,:) :: a
    integer(kind=int32) :: j,k
    numflag=.true.
    do j=1,size(a,1)
       do k=1,size(a,2)
          if (a(j,k)/=a(j,k) .or. real(a(j,k))==inf .or. aimag(a(j,k))==inf .or. &
               real(a(j,k))==-inf .or. aimag(a(j,k))==-inf) then
             numflag=.false.
             return
          end if
       end do
    end do
  end function z_is_number
  
  type(logical) function d_includes_nan(a) result(nanflag)
    real(kind=dp), dimension(:,:) :: a
    integer(kind=int32) :: m,n,j,k
    nanflag=.false.
    m=size(a,1)
    n=size(a,2)
    do j=1,m
       do k=1,n
          if (a(j,k) /= a(j,k)) then
             nanflag=.true.
          end if
       end do
    end do
  end function d_includes_nan

  type(logical) function d_v_includes_nan(a) result(nanflag)
    real(kind=dp), dimension(:) :: a
    integer(kind=int32) :: n,j
    nanflag=.false.
    n=size(a)
    do j=1,n
       if (a(j) /= a(j)) then
          nanflag=.true.
       end if
    end do
  end function d_v_includes_nan

  type(logical) function c_includes_nan(a) result(nanflag)
    complex(kind=dp), dimension(:,:) :: a
    integer(kind=int32) :: m,n,j,k
    nanflag=.false.
    m=size(a,1)
    n=size(a,2)
    do j=1,m
       do k=1,n
          if (a(j,k) /= a(j,k)) then
             nanflag=.true.
          end if
       end do
    end do
  end function c_includes_nan

  type(logical) function c_v_includes_nan(a) result(nanflag)
    complex(kind=dp), dimension(:) :: a
    integer(kind=int32) :: n,j
    nanflag=.false.
    n=size(a)
    do j=1,n
       if (a(j) /= a(j)) then
          nanflag=.true.
       end if
    end do
  end function c_v_includes_nan

end module mod_prec
