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

  ! src/assemble/assemble.f90, 000s, 010s, and 020s
  integer(int32), parameter :: id_d_ub_to_upper=000
  integer(int32), parameter :: id_f_d_ub_to_upper=001
  integer(int32), parameter :: id_c_ub_to_upper=002
  integer(int32), parameter :: id_f_c_ub_to_upper=003
  integer(int32), parameter :: id_d_bt_to_lower=004
  integer(int32), parameter :: id_f_d_bt_to_lower=005
  integer(int32), parameter :: id_c_bt_to_lower=006
  integer(int32), parameter :: id_f_c_bt_to_lower=007
  integer(int32), parameter :: id_d_ubt_to_general=008
  integer(int32), parameter :: id_f_d_ubt_to_general=009
  integer(int32), parameter :: id_c_ubt_to_general=010
  integer(int32), parameter :: id_f_c_ubt_to_general=011
  integer(int32), parameter :: id_d_bv_to_upper=012
  integer(int32), parameter :: id_f_d_bv_to_upper=013
  integer(int32), parameter :: id_c_bv_to_upper=014
  integer(int32), parameter :: id_f_c_bv_to_upper=015
  integer(int32), parameter :: id_d_wb_to_lower=016
  integer(int32), parameter :: id_f_d_wb_to_lower=017
  integer(int32), parameter :: id_c_wb_to_lower=018
  integer(int32), parameter :: id_f_c_wb_to_lower=019
  integer(int32), parameter :: id_d_wbv_to_general=020
  integer(int32), parameter :: id_f_d_wbv_to_general=021
  integer(int32), parameter :: id_c_wbv_to_general=022
  integer(int32), parameter :: id_f_c_wbv_to_general=023

  ! src/recompress/recompression_bv_to_ub, 030s

  integer(int32), parameter :: id_d_recompress_bv_to_ub=030
  integer(int32), parameter :: id_f_d_recompress_bv_to_ub=031
  integer(int32), parameter :: id_c_recompress_bv_to_ub=032
  integer(int32), parameter :: id_f_c_recompress_bv_to_ub=033

  ! src/recompress/recompression_ub_to_bv, 040s

  integer(int32), parameter :: id_d_recompress_ub_to_bv=040
  integer(int32), parameter :: id_f_d_recompress_ub_to_bv=041
  integer(int32), parameter :: id_c_recompress_ub_to_bv=042
  integer(int32), parameter :: id_f_c_recompress_ub_to_bv=043

  ! src/convert/conversion_bv_to_ub, 050s

  integer(int32), parameter :: id_d_convert_bv_to_ub=050
  integer(int32), parameter :: id_f_d_convert_bv_to_ub=051
  integer(int32), parameter :: id_c_convert_bv_to_ub=052
  integer(int32), parameter :: id_f_c_convert_bv_to_ub=053

  ! src/convert/conversion_ub_to_bv, 060s

  integer(int32), parameter :: id_d_convert_ub_to_bv=060
  integer(int32), parameter :: id_f_d_convert_ub_to_bv=061
  integer(int32), parameter :: id_c_convert_ub_to_bv=062
  integer(int32), parameter :: id_f_c_convert_ub_to_bv=063

  ! src/general/general_bv, 070s
  integer(int32), parameter :: id_d_upper_to_bv=070
  integer(int32), parameter :: id_f_d_upper_to_bv=071
  integer(int32), parameter :: id_f_d_general_bv=072
  integer(int32), parameter :: id_c_upper_to_bv=073
  integer(int32), parameter :: id_f_c_upper_to_bv=074
  integer(int32), parameter :: id_f_c_general_bv=075

  ! src/general/general_ub, 080s
  integer(int32), parameter :: id_d_upper_to_ub=080
  integer(int32), parameter :: id_f_d_upper_to_ub=081
  integer(int32), parameter :: id_f_d_general_ub=082
  integer(int32), parameter :: id_c_upper_to_ub=083
  integer(int32), parameter :: id_f_c_upper_to_ub=084
  integer(int32), parameter :: id_f_c_general_ub=085

  ! /src/general/general_bt 090s
  integer(int32), parameter :: id_d_lower_to_bt=090
  integer(int32), parameter :: id_f_d_lower_to_bt=091
  integer(int32), parameter :: id_f_d_general_bt=092
  integer(int32), parameter :: id_c_lower_to_bt=093
  integer(int32), parameter :: id_f_c_lower_to_bt=094
  integer(int32), parameter :: id_f_c_general_bt=095

  ! /src/general/general_wb 100s
  integer(int32), parameter :: id_d_lower_to_wb=100
  integer(int32), parameter :: id_f_d_lower_to_wb=101
  integer(int32), parameter :: id_f_d_general_wb=102
  integer(int32), parameter :: id_c_lower_to_wb=103
  integer(int32), parameter :: id_f_c_lower_to_wb=104
  integer(int32), parameter :: id_f_c_general_wb=105

  ! src/general/general_ubt, 110s
  integer(int32), parameter :: id_d_general_to_ubt=110
  integer(int32), parameter :: id_f_d_general_to_ubt=111
  integer(int32), parameter :: id_f_d_general_ubt=112
  integer(int32), parameter :: id_c_general_to_ubt=113
  integer(int32), parameter :: id_f_c_general_to_ubt=114
  integer(int32), parameter :: id_f_c_general_ubt=115

  ! src/general/general_wbv, 120s
  integer(int32), parameter :: id_d_general_to_wbv=120
  integer(int32), parameter :: id_f_d_general_to_wbv=121
  integer(int32), parameter :: id_f_d_general_wbv=122
  integer(int32), parameter :: id_c_general_to_wbv=123
  integer(int32), parameter :: id_f_c_general_to_wbv=124
  integer(int32), parameter :: id_f_c_general_wbv=125

  ! src/orth/gs 130s
  integer(int32), parameter :: id_d_extend_gs_rows=130
  integer(int32), parameter :: id_c_extend_gs_rows=131
  integer(int32), parameter :: id_d_extend_gs_columns=132
  integer(int32), parameter :: id_c_extend_gs_columns=133

  ! src/orth/nullvec 140s
  integer(int32), parameter :: id_d_lower_left_nullvec=140
  integer(int32), parameter :: id_c_lower_left_nullvec=141
  integer(int32), parameter :: id_d_lower_right_nullvec=142
  integer(int32), parameter :: id_c_lower_right_nullvec=143

  ! src/qr_iteration/qr_iteration 150s
  integer(int32), parameter :: id_ss_r1_qr=150
  integer(int32), parameter :: id_f_ss_r1_qr=151
  integer(int32), parameter :: id_ss_qr_iteration=152
  integer(int32), parameter :: id_f_ss_qr_iteration=153

  ! src/qr_factorization/qr_factorization 160s
  integer(int32), parameter :: id_d_reduce_lbw_bv_to_ub=160
  integer(int32), parameter :: id_f_d_reduce_lbw_bv_to_ub=161
  integer(int32), parameter :: id_c_reduce_lbw_bv_to_ub=162
  integer(int32), parameter :: id_f_c_reduce_lbw_bv_to_ub=163
  integer(int32), parameter :: id_d_qr_bv_to_ub=164
  integer(int32), parameter :: id_f_d_qr_bv_to_ub=165
  integer(int32), parameter :: id_c_qr_bv_to_ub=166
  integer(int32), parameter :: id_f_c_qr_bv_to_ub=167

  ! src/solve/back_substitution 170s and 180s
  integer(int32), parameter :: id_d_back_substitution_ub=170
  integer(int32), parameter :: id_f_d_back_substitution_ub=171
  integer(int32), parameter :: id_c_back_substitution_ub=172
  integer(int32), parameter :: id_f_c_back_substitution_ub=173
  integer(int32), parameter :: id_d_v_back_substitution_ub=174
  integer(int32), parameter :: id_f_d_v_back_substitution_ub=175
  integer(int32), parameter :: id_c_v_back_substitution_ub=176
  integer(int32), parameter :: id_f_c_v_back_substitution_ub=177

  integer(int32), parameter :: id_d_forward_substitution_bv=178
  integer(int32), parameter :: id_f_d_forward_substitution_bv=179
  integer(int32), parameter :: id_d_v_forward_substitution_bv=180
  integer(int32), parameter :: id_f_d_v_forward_substitution_bv=181
  integer(int32), parameter :: id_c_forward_substitution_bv=118
  integer(int32), parameter :: id_f_c_forward_substitution_bv=119
  integer(int32), parameter :: id_c_v_forward_substitution_bv=180
  integer(int32), parameter :: id_f_c_v_forward_substitution_bv=181

  ! src/sweeps 190s and 200s
  integer(int32), parameter :: id_d_sweeps_times_ub=190
  integer(int32), parameter :: id_f_d_sweeps_times_ub=191
  integer(int32), parameter :: id_c_sweeps_times_ub=192
  integer(int32), parameter :: id_f_c_sweeps_times_ub=193
  integer(int32), parameter :: id_d_bv_times_sweeps=194
  integer(int32), parameter :: id_f_d_bv_times_sweeps=195
  integer(int32), parameter :: id_c_bv_times_sweeps=196
  integer(int32), parameter :: id_f_c_bv_times_sweeps=197
  integer(int32), parameter :: id_d_trp_sweeps_times_bv=198
  integer(int32), parameter :: id_f_d_trp_sweeps_times_bv=199
  integer(int32), parameter :: id_c_trp_sweeps_times_bv=200
  integer(int32), parameter :: id_f_c_trp_sweeps_times_bv=201
  integer(int32), parameter :: id_d_ub_times_trp_sweeps=202
  integer(int32), parameter :: id_f_d_ub_times_trp_sweeps=203
  integer(int32), parameter :: id_c_ub_times_trp_sweeps=204
  integer(int32), parameter :: id_f_c_ub_times_trp_sweeps=205

  ! src/update 210s
  integer(int32), parameter :: id_d_r1_update_ub_to_bv=210
  integer(int32), parameter :: id_f_d_r1_update_ub_to_bv=211
  integer(int32), parameter :: id_c_r1_update_ub_to_bv=212
  integer(int32), parameter :: id_f_c_r1_update_ub_to_bv=213
  integer(int32), parameter :: id_d_e1v_update_ub_to_bv=214
  integer(int32), parameter :: id_f_d_e1v_update_ub_to_bv=215
  integer(int32), parameter :: id_c_e1v_update_ub_to_bv=216
  integer(int32), parameter :: id_f_c_e1v_update_ub_to_bv=217

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
    err%routines(1:err%rptr)=-1
    err%rptr=1
    err%code=0

  end subroutine clear_error

end module error_id
