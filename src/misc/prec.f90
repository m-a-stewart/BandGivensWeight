module mod_prec

  integer, parameter :: dp = kind(0.0d0)
  integer, parameter :: int32 = selected_int_kind(8)
  real(kind=dp), parameter :: eps=epsilon(1.0_dp)/2

  public

end module mod_prec
