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

  ! src/recompress/recompression_bv_to_ub, 010s

  integer(int32), parameter :: id_d_recompress_bv_to_ub=010
  integer(int32), parameter :: id_f_d_recompress_bv_to_ub=011
  integer(int32), parameter :: id_c_recompress_bv_to_ub=012
  integer(int32), parameter :: id_f_c_recompress_bv_to_ub=013

  ! src/recompress/recompression_ub_to_bv, 020s

  integer(int32), parameter :: id_d_recompress_ub_to_bv=020
  integer(int32), parameter :: id_f_d_recompress_ub_to_bv=021
  integer(int32), parameter :: id_c_recompress_ub_to_bv=022
  integer(int32), parameter :: id_f_c_recompress_ub_to_bv=023

  ! src/convert/conversion_bv_to_ub, 030s

  integer(int32), parameter :: id_d_convert_bv_to_ub=030
  integer(int32), parameter :: id_f_d_convert_bv_to_ub=031
  integer(int32), parameter :: id_c_convert_bv_to_ub=032
  integer(int32), parameter :: id_f_c_convert_bv_to_ub=033

  ! src/convert/conversion_ub_to_bv, 040s

  integer(int32), parameter :: id_d_convert_ub_to_bv=040
  integer(int32), parameter :: id_f_d_convert_ub_to_bv=041
  integer(int32), parameter :: id_c_convert_ub_to_bv=042
  integer(int32), parameter :: id_f_c_convert_ub_to_bv=043

  ! src/general/general_bv, 050s
  integer(int32), parameter :: id_d_upper_to_bv=050
  integer(int32), parameter :: id_f_d_upper_to_bv=051
  integer(int32), parameter :: id_f_d_general_bv=052
  integer(int32), parameter :: id_c_upper_to_bv=053
  integer(int32), parameter :: id_f_c_upper_to_bv=054
  integer(int32), parameter :: id_f_c_general_bv=055

  ! src/general/general_ub, 060s
  integer(int32), parameter :: id_d_upper_to_ub=060
  integer(int32), parameter :: id_f_d_upper_to_ub=061
  integer(int32), parameter :: id_f_d_general_ub=062
  integer(int32), parameter :: id_c_upper_to_ub=063
  integer(int32), parameter :: id_f_c_upper_to_ub=064
  integer(int32), parameter :: id_f_c_general_ub=065

  ! /src/general/general_bt 070s
  integer(int32), parameter :: id_d_lower_to_bt=070
  integer(int32), parameter :: id_f_d_lower_to_bt=071
  integer(int32), parameter :: id_f_d_general_bt=072
  integer(int32), parameter :: id_c_lower_to_bt=073
  integer(int32), parameter :: id_f_c_lower_to_bt=074
  integer(int32), parameter :: id_f_c_general_bt=075

  ! /src/general/general_wb 080s
  integer(int32), parameter :: id_d_lower_to_wb=080
  integer(int32), parameter :: id_f_d_lower_to_wb=081
  integer(int32), parameter :: id_f_d_general_wb=082
  integer(int32), parameter :: id_c_lower_to_wb=083
  integer(int32), parameter :: id_f_c_lower_to_wb=084
  integer(int32), parameter :: id_f_c_general_wb=085

  ! src/general/general_ubt, 090s
  integer(int32), parameter :: id_d_general_to_ubt=090
  integer(int32), parameter :: id_f_d_general_to_ubt=091
  integer(int32), parameter :: id_f_d_general_ubt=092
  integer(int32), parameter :: id_c_general_to_ubt=093
  integer(int32), parameter :: id_f_c_general_to_ubt=094
  integer(int32), parameter :: id_f_c_general_ubt=095

  ! src/general/general_wbv, 100s
  integer(int32), parameter :: id_d_general_to_wbv=100
  integer(int32), parameter :: id_f_d_general_to_wbv=101
  integer(int32), parameter :: id_f_d_general_wbv=102
  integer(int32), parameter :: id_c_general_to_wbv=103
  integer(int32), parameter :: id_f_c_general_to_wbv=104
  integer(int32), parameter :: id_f_c_general_wbv=105

  ! src/orth/gs 110s
  integer(int32), parameter :: id_d_extend_gs_rows=110
  integer(int32), parameter :: id_c_extend_gs_rows=111
  integer(int32), parameter :: id_d_extend_gs_columns=112
  integer(int32), parameter :: id_c_extend_gs_columns=113

  ! src/orth/nullvec 120s
  integer(int32), parameter :: id_d_lower_left_nullvec=120
  integer(int32), parameter :: id_c_lower_left_nullvec=121
  integer(int32), parameter :: id_d_lower_right_nullvec=122
  integer(int32), parameter :: id_c_lower_right_nullvec=123

  ! src/qr_iteration/qr_iteration 130s
  integer(int32), parameter :: id_ss_r1_qr=130
  integer(int32), parameter :: id_f_ss_r1_qr=131
  integer(int32), parameter :: id_ss_qr_iteration=132
  integer(int32), parameter :: id_f_ss_qr_iteration=133

  ! src/qr_factorization/qr_factorization 140s
  integer(int32), parameter :: id_d_reduce_lbw_bv_to_ub=140
  integer(int32), parameter :: id_f_d_reduce_lbw_bv_to_ub=141
  integer(int32), parameter :: id_c_reduce_lbw_bv_to_ub=142
  integer(int32), parameter :: id_f_c_reduce_lbw_bv_to_ub=143
  integer(int32), parameter :: id_d_qr_bv_to_ub=144
  integer(int32), parameter :: id_f_d_qr_bv_to_ub=145
  integer(int32), parameter :: id_c_qr_bv_to_ub=146
  integer(int32), parameter :: id_f_c_qr_bv_to_ub=147

  ! src/solve/back_substitution 150s and 160s
  integer(int32), parameter :: id_d_back_substitution_ub=150
  integer(int32), parameter :: id_f_d_back_substitution_ub=151
  integer(int32), parameter :: id_c_back_substitution_ub=152
  integer(int32), parameter :: id_f_c_back_substitution_ub=153
  integer(int32), parameter :: id_d_v_back_substitution_ub=154
  integer(int32), parameter :: id_f_d_v_back_substitution_ub=155
  integer(int32), parameter :: id_c_v_back_substitution_ub=156
  integer(int32), parameter :: id_f_c_v_back_substitution_ub=157

  integer(int32), parameter :: id_d_forward_substitution_bv=158
  integer(int32), parameter :: id_f_d_forward_substitution_bv=159
  integer(int32), parameter :: id_d_v_forward_substitution_bv=160
  integer(int32), parameter :: id_f_d_v_forward_substitution_bv=161
  integer(int32), parameter :: id_c_forward_substitution_bv=118
  integer(int32), parameter :: id_f_c_forward_substitution_bv=119
  integer(int32), parameter :: id_c_v_forward_substitution_bv=160
  integer(int32), parameter :: id_f_c_v_forward_substitution_bv=161

  ! src/sweeps 170s and 180s
  integer(int32), parameter :: id_d_sweeps_times_ub=170
  integer(int32), parameter :: id_f_d_sweeps_times_ub=171
  integer(int32), parameter :: id_c_sweeps_times_ub=172
  integer(int32), parameter :: id_f_c_sweeps_times_ub=173
  integer(int32), parameter :: id_d_bv_times_sweeps=174
  integer(int32), parameter :: id_f_d_bv_times_sweeps=175
  integer(int32), parameter :: id_c_bv_times_sweeps=176
  integer(int32), parameter :: id_f_c_bv_times_sweeps=177
  integer(int32), parameter :: id_d_trp_sweeps_times_bv=178
  integer(int32), parameter :: id_f_d_trp_sweeps_times_bv=179
  integer(int32), parameter :: id_c_trp_sweeps_times_bv=180
  integer(int32), parameter :: id_f_c_trp_sweeps_times_bv=181
  integer(int32), parameter :: id_d_ub_times_trp_sweeps=182
  integer(int32), parameter :: id_f_d_ub_times_trp_sweeps=183
  integer(int32), parameter :: id_c_ub_times_trp_sweeps=184
  integer(int32), parameter :: id_f_c_ub_times_trp_sweeps=185

  ! src/update 190s
  integer(int32), parameter :: id_d_r1_update_ub_to_bv=190
  integer(int32), parameter :: id_f_d_r1_update_ub_to_bv=191
  integer(int32), parameter :: id_c_r1_update_ub_to_bv=192
  integer(int32), parameter :: id_f_c_r1_update_ub_to_bv=193
  integer(int32), parameter :: id_d_e1v_update_ub_to_bv=194
  integer(int32), parameter :: id_f_d_e1v_update_ub_to_bv=195
  integer(int32), parameter :: id_c_e1v_update_ub_to_bv=196
  integer(int32), parameter :: id_f_c_e1v_update_ub_to_bv=197

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
