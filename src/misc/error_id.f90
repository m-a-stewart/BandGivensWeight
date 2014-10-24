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

  ! src/convert/conversion_wb_to_bt, 070s
  integer(int32), parameter :: id_d_convert_wb_to_bt=070
  integer(int32), parameter :: id_f_d_convert_wb_to_bt=071
  integer(int32), parameter :: id_c_convert_wb_to_bt=072
  integer(int32), parameter :: id_f_c_convert_wb_to_bt=073

  ! src/convert/conversion_bt_to_wb, 080s
  integer(int32), parameter :: id_d_convert_bt_to_wb=080
  integer(int32), parameter :: id_f_d_convert_bt_to_wb=081
  integer(int32), parameter :: id_c_convert_bt_to_wb=082
  integer(int32), parameter :: id_f_c_convert_bt_to_wb=083

  ! src/convert/conversion_wbv_to_ubt, 090s
  integer(int32), parameter :: id_d_convert_wbv_to_ubt=090
  integer(int32), parameter :: id_f_d_convert_wbv_to_ubt=091
  integer(int32), parameter :: id_c_convert_wbv_to_ubt=092
  integer(int32), parameter :: id_f_c_convert_wbv_to_ubt=093

  ! src/convert/conversion_ubt_to_wbv, 090s
  integer(int32), parameter :: id_d_convert_ubt_to_wbv=100
  integer(int32), parameter :: id_f_d_convert_ubt_to_wbv=101
  integer(int32), parameter :: id_c_convert_ubt_to_wbv=102
  integer(int32), parameter :: id_f_c_convert_ubt_to_wbv=103

  ! src/general/general_bv, 110s
  integer(int32), parameter :: id_d_upper_to_bv=110
  integer(int32), parameter :: id_f_d_upper_to_bv=111
  integer(int32), parameter :: id_f_d_general_bv=112
  integer(int32), parameter :: id_c_upper_to_bv=113
  integer(int32), parameter :: id_f_c_upper_to_bv=114
  integer(int32), parameter :: id_f_c_general_bv=115

  ! src/general/general_ub, 120s
  integer(int32), parameter :: id_d_upper_to_ub=120
  integer(int32), parameter :: id_f_d_upper_to_ub=121
  integer(int32), parameter :: id_f_d_general_ub=122
  integer(int32), parameter :: id_c_upper_to_ub=123
  integer(int32), parameter :: id_f_c_upper_to_ub=124
  integer(int32), parameter :: id_f_c_general_ub=125

  ! /src/general/general_bt 130s
  integer(int32), parameter :: id_d_lower_to_bt=130
  integer(int32), parameter :: id_f_d_lower_to_bt=131
  integer(int32), parameter :: id_f_d_general_bt=132
  integer(int32), parameter :: id_c_lower_to_bt=133
  integer(int32), parameter :: id_f_c_lower_to_bt=134
  integer(int32), parameter :: id_f_c_general_bt=135

  ! /src/general/general_wb 140s
  integer(int32), parameter :: id_d_lower_to_wb=140
  integer(int32), parameter :: id_f_d_lower_to_wb=141
  integer(int32), parameter :: id_f_d_general_wb=142
  integer(int32), parameter :: id_c_lower_to_wb=143
  integer(int32), parameter :: id_f_c_lower_to_wb=144
  integer(int32), parameter :: id_f_c_general_wb=145

  ! src/general/general_ubt, 150s
  integer(int32), parameter :: id_d_general_to_ubt=150
  integer(int32), parameter :: id_f_d_general_to_ubt=151
  integer(int32), parameter :: id_f_d_general_ubt=152
  integer(int32), parameter :: id_c_general_to_ubt=153
  integer(int32), parameter :: id_f_c_general_to_ubt=154
  integer(int32), parameter :: id_f_c_general_ubt=155

  ! src/general/general_wbv, 160s
  integer(int32), parameter :: id_d_general_to_wbv=160
  integer(int32), parameter :: id_f_d_general_to_wbv=161
  integer(int32), parameter :: id_f_d_general_wbv=162
  integer(int32), parameter :: id_c_general_to_wbv=163
  integer(int32), parameter :: id_f_c_general_to_wbv=164
  integer(int32), parameter :: id_f_c_general_wbv=165

  ! src/orth/gs 170s
  integer(int32), parameter :: id_d_extend_gs_rows=170
  integer(int32), parameter :: id_c_extend_gs_rows=171
  integer(int32), parameter :: id_d_extend_gs_columns=172
  integer(int32), parameter :: id_c_extend_gs_columns=173

  ! src/orth/nullvec 180s
  integer(int32), parameter :: id_d_lower_left_nullvec=180
  integer(int32), parameter :: id_c_lower_left_nullvec=181
  integer(int32), parameter :: id_d_lower_right_nullvec=182
  integer(int32), parameter :: id_c_lower_right_nullvec=183

  ! src/qr_iteration/qr_iteration 190s
  integer(int32), parameter :: id_ss_r1_qr=190
  integer(int32), parameter :: id_f_ss_r1_qr=191
  integer(int32), parameter :: id_ss_qr_iteration=192
  integer(int32), parameter :: id_f_ss_qr_iteration=193

  ! src/qr_factorization/qr_factorization 200s
  integer(int32), parameter :: id_d_reduce_lbw_bv_to_ub=200
  integer(int32), parameter :: id_f_d_reduce_lbw_bv_to_ub=201
  integer(int32), parameter :: id_c_reduce_lbw_bv_to_ub=202
  integer(int32), parameter :: id_f_c_reduce_lbw_bv_to_ub=203
  integer(int32), parameter :: id_d_qr_bv_to_ub=204
  integer(int32), parameter :: id_f_d_qr_bv_to_ub=205
  integer(int32), parameter :: id_c_qr_bv_to_ub=206
  integer(int32), parameter :: id_f_c_qr_bv_to_ub=207

  ! src/solve/back_substitution 210s and 220s
  integer(int32), parameter :: id_d_back_substitution_ub=210
  integer(int32), parameter :: id_f_d_back_substitution_ub=211
  integer(int32), parameter :: id_c_back_substitution_ub=212
  integer(int32), parameter :: id_f_c_back_substitution_ub=213
  integer(int32), parameter :: id_d_v_back_substitution_ub=214
  integer(int32), parameter :: id_f_d_v_back_substitution_ub=215
  integer(int32), parameter :: id_c_v_back_substitution_ub=216
  integer(int32), parameter :: id_f_c_v_back_substitution_ub=217

  integer(int32), parameter :: id_d_forward_substitution_bv=218
  integer(int32), parameter :: id_f_d_forward_substitution_bv=219
  integer(int32), parameter :: id_d_v_forward_substitution_bv=220
  integer(int32), parameter :: id_f_d_v_forward_substitution_bv=221
  integer(int32), parameter :: id_c_forward_substitution_bv=222
  integer(int32), parameter :: id_f_c_forward_substitution_bv=223
  integer(int32), parameter :: id_c_v_forward_substitution_bv=224
  integer(int32), parameter :: id_f_c_v_forward_substitution_bv=225

  ! src/sweeps 230s and 240s
  integer(int32), parameter :: id_d_sweeps_times_ub=230
  integer(int32), parameter :: id_f_d_sweeps_times_ub=231
  integer(int32), parameter :: id_c_sweeps_times_ub=232
  integer(int32), parameter :: id_f_c_sweeps_times_ub=233
  integer(int32), parameter :: id_d_bv_times_sweeps=234
  integer(int32), parameter :: id_f_d_bv_times_sweeps=235
  integer(int32), parameter :: id_c_bv_times_sweeps=236
  integer(int32), parameter :: id_f_c_bv_times_sweeps=237
  integer(int32), parameter :: id_d_trp_sweeps_times_bv=238
  integer(int32), parameter :: id_f_d_trp_sweeps_times_bv=239
  integer(int32), parameter :: id_c_trp_sweeps_times_bv=240
  integer(int32), parameter :: id_f_c_trp_sweeps_times_bv=241
  integer(int32), parameter :: id_d_ub_times_trp_sweeps=242
  integer(int32), parameter :: id_f_d_ub_times_trp_sweeps=243
  integer(int32), parameter :: id_c_ub_times_trp_sweeps=244
  integer(int32), parameter :: id_f_c_ub_times_trp_sweeps=245

  ! src/update 250s
  integer(int32), parameter :: id_d_r1_update_ub_to_bv=250
  integer(int32), parameter :: id_f_d_r1_update_ub_to_bv=251
  integer(int32), parameter :: id_c_r1_update_ub_to_bv=252
  integer(int32), parameter :: id_f_c_r1_update_ub_to_bv=253
  integer(int32), parameter :: id_d_e1v_update_ub_to_bv=254
  integer(int32), parameter :: id_f_d_e1v_update_ub_to_bv=255
  integer(int32), parameter :: id_c_e1v_update_ub_to_bv=256
  integer(int32), parameter :: id_f_c_e1v_update_ub_to_bv=257

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
