module prec

integer, parameter :: dp = selected_real_kind(12)
integer, parameter :: int32 = selected_int_kind(8)
real(kind=dp), parameter :: eps=epsilon(1.0_dp)/2

end module prec
