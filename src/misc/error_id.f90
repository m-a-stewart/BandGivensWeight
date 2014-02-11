module error_id
use prec
implicit none

integer(kind=int32), parameter :: routine_name_length=30
integer(kind=int32), parameter :: max_errors=15
integer(kind=int32), parameter :: error_message_length=50
integer(kind=int32), parameter :: max_routines=10

type routine_info
   integer(kind=int32) :: routine_id
   character(len=routine_name_length) :: routine_name
   character(len=error_message_length), dimension(max_errors) :: error_messages=''
end type routine_info

type error_info
   integer(kind=int32) :: code=0, rptr=1
   integer(kind=int32), dimension(max_routines) :: routines=-1
end type error_info

! src/assemble/assemble.f90, 000s
integer(int32), parameter :: id_d_ub_to_upper=000
integer(int32), parameter :: id_f_d_ub_to_upper=001
integer(int32), parameter :: id_c_ub_to_upper=002
integer(int32), parameter :: id_f_c_ub_to_upper=003
integer(int32), parameter :: id_d_bv_to_upper=004
integer(int32), parameter :: id_f_d_bv_to_upper=005
integer(int32), parameter :: id_c_bv_to_upper=006
integer(int32), parameter :: id_f_c_bv_to_upper=007

! src/compress/compressions_bv_to_ub, 010s

integer(int32), parameter :: id_d_compress_bv_to_ub=010
integer(int32), parameter :: id_f_d_compress_bv_to_ub=011
integer(int32), parameter :: id_c_compress_bv_to_ub=012
integer(int32), parameter :: id_f_c_compress_bv_to_ub=013

! src/compress/compressions_ub_to_bv, 020s

integer(int32), parameter :: id_d_compress_ub_to_bv=020
integer(int32), parameter :: id_f_d_compress_ub_to_bv=021
integer(int32), parameter :: id_c_compress_ub_to_bv=022
integer(int32), parameter :: id_f_c_compress_ub_to_bv=023

! src/convert/conversions_bv_to_ub, 030s

integer(int32), parameter :: id_d_convert_bv_to_ub=030
integer(int32), parameter :: id_f_d_convert_bv_to_ub=031
integer(int32), parameter :: id_c_convert_bv_to_ub=032
integer(int32), parameter :: id_f_c_convert_bv_to_ub=033

! src/convert/conversions_ub_to_bv, 040s

integer(int32), parameter :: id_d_convert_ub_to_bv=040
integer(int32), parameter :: id_f_d_convert_ub_to_bv=041
integer(int32), parameter :: id_c_convert_ub_to_bv=042
integer(int32), parameter :: id_f_c_convert_ub_to_bv=043

! src/decomp/general_bv, 050s
integer(int32), parameter :: id_d_upper_to_bv=050
integer(int32), parameter :: id_f_d_upper_to_bv=051
integer(int32), parameter :: id_c_upper_to_bv=052
integer(int32), parameter :: id_f_c_upper_to_bv=053

! src/decomp/general_ub, 060s
integer(int32), parameter :: id_d_upper_to_ub=060
integer(int32), parameter :: id_f_d_upper_to_ub=061
integer(int32), parameter :: id_c_upper_to_ub=062
integer(int32), parameter :: id_f_c_upper_to_ub=063

! src/orth/gs 070s
integer(int32), parameter :: id_d_extend_gs_rows=070
integer(int32), parameter :: id_c_extend_gs_rows=071
integer(int32), parameter :: id_d_extend_gs_columns=072
integer(int32), parameter :: id_c_extend_gs_columns=073

! src/orth/nullvec 080s
integer(int32), parameter :: id_d_lower_left_nullvec=080
integer(int32), parameter :: id_c_lower_left_nullvec=081
integer(int32), parameter :: id_d_lower_right_nullvec=082
integer(int32), parameter :: id_c_lower_right_nullvec=083

! src/qr_iteration/qr_iteration 090s
integer(int32), parameter :: id_ss_qr=090
integer(int32), parameter :: id_f_ss_qr=091
integer(int32), parameter :: id_ss_qr_iteration=092
integer(int32), parameter :: id_f_ss_qr_iteration=093

contains

subroutine add_id(err,id)
  type(error_info), intent(inout) :: err
  integer(kind=int32), intent(in) :: id
  if (err%rptr <= max_routines) then
     err%routines(err%rptr)=id
     err%rptr=err%rptr+1
  end if
end subroutine add_id

subroutine set_error(err,code, id)
  type(error_info), intent(inout) :: err
  integer(kind=int32), intent(in) :: code, id
  err%rptr=2
  err%code=code
  err%routines(1)=id
end subroutine set_error

subroutine clear_error(err)
  type(error_info), intent(inout) :: err
  err%rptr=1
  err%code=0
end subroutine clear_error

end module error_id
