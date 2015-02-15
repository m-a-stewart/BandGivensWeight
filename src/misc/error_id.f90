module mod_error_id
  use mod_prec
  implicit none

  integer(kind=int32), parameter :: routine_name_length=30
  integer(kind=int32), parameter :: max_errors=15
  integer(kind=int32), parameter :: error_message_length=50
  integer(kind=int32), parameter :: max_routines=10

  integer(kind=int32), parameter :: max_error_code=300

  public

  type routine_info
     integer(kind=int32) :: routine_id
     character(len=routine_name_length) :: routine_name
     character(len=error_message_length), dimension(max_errors) :: error_messages=''
  end type routine_info

  type error_info
     integer(kind=int32) :: code=0, rptr=1
     integer(kind=int32), dimension(max_routines) :: routines=-1
  end type error_info

  type(routine_info), parameter :: info_empty=routine_info(0, &
       '', [ character(len=error_message_length) :: '' ] )

  type(routine_info), dimension(:), allocatable :: info_index

  ! src/assemble/assemble.f90, 000s, 010s, and 020s
  integer(int32), parameter :: id_d_ub_to_upper=001
  integer(int32), parameter :: id_f_d_ub_to_upper=002
  integer(int32), parameter :: id_c_ub_to_upper=003
  integer(int32), parameter :: id_f_c_ub_to_upper=004

  integer(int32), parameter :: id_d_bt_to_lower=005
  integer(int32), parameter :: id_f_d_bt_to_lower=006
  integer(int32), parameter :: id_c_bt_to_lower=007
  integer(int32), parameter :: id_f_c_bt_to_lower=008

  integer(int32), parameter :: id_d_ubt_to_general=009
  integer(int32), parameter :: id_f_d_ubt_to_general=010
  integer(int32), parameter :: id_c_ubt_to_general=011
  integer(int32), parameter :: id_f_c_ubt_to_general=012

  integer(int32), parameter :: id_d_bv_to_upper=013
  integer(int32), parameter :: id_f_d_bv_to_upper=014
  integer(int32), parameter :: id_c_bv_to_upper=015
  integer(int32), parameter :: id_f_c_bv_to_upper=016

  integer(int32), parameter :: id_d_wb_to_lower=017
  integer(int32), parameter :: id_f_d_wb_to_lower=018
  integer(int32), parameter :: id_c_wb_to_lower=019
  integer(int32), parameter :: id_f_c_wb_to_lower=020

  integer(int32), parameter :: id_d_wbv_to_general=021
  integer(int32), parameter :: id_f_d_wbv_to_general=022
  integer(int32), parameter :: id_c_wbv_to_general=023
  integer(int32), parameter :: id_f_c_wbv_to_general=024

  type(routine_info), parameter :: info_d_ub_to_upper=routine_info(id_d_ub_to_upper, &
       'd_ub_to_upper', [ character(len=error_message_length) :: 'Size error in A.' ] )

  type(routine_info), parameter :: info_c_ub_to_upper=routine_info(id_c_ub_to_upper, &
       'c_ub_to_upper', [ character(len=error_message_length) :: 'Size error in A.' ] )

  type(routine_info), parameter :: info_d_bt_to_lower=routine_info(id_d_bt_to_lower, &
       'd_bt_to_lower', [ character(len=error_message_length) :: 'Size error in A.' ] )

  type(routine_info), parameter :: info_c_bt_to_lower=routine_info(id_c_bt_to_lower, &
       'c_bt_to_lower', [ character(len=error_message_length) :: 'Size error in A.' ] )

  type(routine_info), parameter :: info_d_ubt_to_general=routine_info(id_d_ubt_to_general, &
       'd_ubt_to_general', [ character(len=error_message_length) :: 'Size error in A.' ] )

  type(routine_info), parameter :: info_c_ubt_to_general=routine_info(id_c_ubt_to_general, &
       'c_ubt_to_general', [ character(len=error_message_length) :: 'Size error in A.' ] )

  type(routine_info), parameter :: info_d_bv_to_upper=routine_info(id_d_bv_to_upper, &
       'd_bv_to_upper', [ character(len=error_message_length) :: 'Size error in A.' ] )
  
  type(routine_info), parameter :: info_c_bv_to_upper=routine_info(id_c_bv_to_upper, &
       'c_bv_to_upper', [ character(len=error_message_length) :: 'Size error in A.' ] )

  type(routine_info), parameter :: info_d_wb_to_lower=routine_info(id_d_wb_to_lower, &
       'd_wb_to_lower', [ character(len=error_message_length) :: 'Size error in A.' ] )
  
  type(routine_info), parameter :: info_c_wb_to_lower=routine_info(id_c_wb_to_lower, &
       'c_wb_to_lower', [ character(len=error_message_length) :: 'Size error in A.' ] )

  type(routine_info), parameter :: info_d_wbv_to_general=routine_info(id_d_wbv_to_general, &
       'd_wbv_to_general', [ character(len=error_message_length) :: 'Size error in A.' ] )
  
  type(routine_info), parameter :: info_c_wbv_to_general=routine_info(id_c_wbv_to_general, &
       'c_wbv_to_general', [ character(len=error_message_length) :: 'Size error in A.' ] )

  ! src/convert/convert_bv_to_ub, 050s
  integer(int32), parameter :: id_d_convert_bv_to_ub=050
  integer(int32), parameter :: id_f_d_convert_bv_to_ub=051
  integer(int32), parameter :: id_c_convert_bv_to_ub=052
  integer(int32), parameter :: id_f_c_convert_bv_to_ub=053

  type(routine_info), parameter :: info_d_convert_bv_to_ub=routine_info(id_d_convert_bv_to_ub, &
       'd_convert_bv_to_ub', &
       [ character(len=error_message_length) :: 'n<1', 'Insufficient storage in bv.', &
       'Insufficient Storage in ub.', 'ub%n /= bv%n' ] )

  type(routine_info), parameter :: info_c_convert_bv_to_ub=routine_info(id_c_convert_bv_to_ub, &
       'c_convert_bv_to_ub', &
       [ character(len=error_message_length) :: 'n<1', 'Insufficient storage in bv.', &
       'Insufficient Storage in ub.', 'ub%n /= bv%n' ] )

  ! src/convert/convert_ub_to_bv, 060s
  integer(int32), parameter :: id_d_convert_ub_to_bv=060
  integer(int32), parameter :: id_f_d_convert_ub_to_bv=061
  integer(int32), parameter :: id_c_convert_ub_to_bv=062
  integer(int32), parameter :: id_f_c_convert_ub_to_bv=063

  type(routine_info), parameter :: info_d_convert_ub_to_bv=routine_info(id_d_convert_ub_to_bv, &
       'd_convert_ub_to_bv', &
       [ character(len=error_message_length) :: 'n<1', 'Insufficient storage in ub', &
       'Insufficient storage in bv.', 'ub%n /= bv%n' ] )

  type(routine_info), parameter :: info_c_convert_ub_to_bv=routine_info(id_c_convert_ub_to_bv, &
       'c_convert_ub_to_bv', &
       [ character(len=error_message_length) :: 'Insufficient storage in ub', &
       'Insufficient storage in bv.', 'ub%n /= bv%n' ] )

  ! src/convert/convert_wb_to_bt, 070s
  integer(int32), parameter :: id_d_convert_wb_to_bt=070
  integer(int32), parameter :: id_f_d_convert_wb_to_bt=071
  integer(int32), parameter :: id_c_convert_wb_to_bt=072
  integer(int32), parameter :: id_f_c_convert_wb_to_bt=073

  type(routine_info), parameter :: info_d_convert_wb_to_bt=routine_info(id_d_convert_wb_to_bt, &
       'd_convert_wb_to_bt', &
       [ character(len=error_message_length) :: 'n<1', 'Insufficient storage in wb.', &
       'Insufficient Storage in bt.', 'bt%n /= wb%n' ] )

  type(routine_info), parameter :: info_c_convert_wb_to_bt=routine_info(id_c_convert_wb_to_bt, &
       'c_convert_wb_to_bt', &
       [ character(len=error_message_length) :: 'n<1', 'Insufficient storage in wb.', &
       'Insufficient Storage in bt.', 'bt%n /= wb%n' ] )

  ! src/convert/convert_bt_to_wb, 080s
  integer(int32), parameter :: id_d_convert_bt_to_wb=080
  integer(int32), parameter :: id_f_d_convert_bt_to_wb=081
  integer(int32), parameter :: id_c_convert_bt_to_wb=082
  integer(int32), parameter :: id_f_c_convert_bt_to_wb=083

  type(routine_info), parameter :: info_d_convert_bt_to_wb=routine_info(id_d_convert_bt_to_wb, &
       'd_convert_bt_to_wb', &
       [ character(len=error_message_length) :: 'n<1', 'Insufficient storage in bt.', &
       'Insufficient Storage in wb.', 'wb%n /= bt%n' ] )

  type(routine_info), parameter :: info_c_convert_bt_to_wb=routine_info(id_c_convert_bt_to_wb, &
       'c_convert_bt_to_wb', &
       [ character(len=error_message_length) :: 'n<1', 'Insufficient storage in bt.', &
       'Insufficient Storage in wb.', 'wb%n /= bt%n' ] )

  ! src/convert/convert_wbv_to_ubt, 090s
  integer(int32), parameter :: id_d_convert_wbv_to_ubt=090
  integer(int32), parameter :: id_f_d_convert_wbv_to_ubt=091
  integer(int32), parameter :: id_c_convert_wbv_to_ubt=092
  integer(int32), parameter :: id_f_c_convert_wbv_to_ubt=093

  type(routine_info), parameter :: info_d_convert_wbv_to_ubt=routine_info(id_d_convert_wbv_to_ubt, &
       'd_convert_wbv_to_ubt', &
       [ character(len=error_message_length) :: 'n<1', 'Insufficient storage in wbv', &
       'Insufficient storage in ubt.', 'wbv%n /= ubt%n' ] )

  type(routine_info), parameter :: info_c_convert_wbv_to_ubt=routine_info(id_c_convert_wbv_to_ubt, &
       'c_convert_wbv_to_ubt', &
       [ character(len=error_message_length) :: 'Insufficient storage in wbv', &
       'Insufficient storage in ubt.', 'wbv%n /= ubt%n' ] )

  ! src/convert/convert_ubt_to_wbv, 090s
  integer(int32), parameter :: id_d_convert_ubt_to_wbv=100
  integer(int32), parameter :: id_f_d_convert_ubt_to_wbv=101
  integer(int32), parameter :: id_c_convert_ubt_to_wbv=102
  integer(int32), parameter :: id_f_c_convert_ubt_to_wbv=103

  type(routine_info), parameter :: info_d_convert_ubt_to_wbv=routine_info(id_d_convert_ubt_to_wbv, &
       'd_convert_ubt_to_wbv', &
       [ character(len=error_message_length) :: 'n<1', 'Insufficient storage in ubt', &
       'Insufficient storage in wbv.', 'ubt%n /= wbv%n' ] )

  type(routine_info), parameter :: info_c_convert_ubt_to_wbv=routine_info(id_c_convert_ubt_to_wbv, &
       'c_convert_ubt_to_wbv', &
       [ character(len=error_message_length) :: 'Insufficient storage in ubt', &
       'Insufficient storage in wbv.', 'ubt%n /= wbv%n' ] )

  ! src/general/general_bv, 110s
  integer(int32), parameter :: id_d_upper_to_bv=110
  integer(int32), parameter :: id_f_d_upper_to_bv=111
  integer(int32), parameter :: id_f_d_general_bv=112
  integer(int32), parameter :: id_c_upper_to_bv=113
  integer(int32), parameter :: id_f_c_upper_to_bv=114
  integer(int32), parameter :: id_f_c_general_bv=115

  type(routine_info), parameter :: info_d_upper_to_bv=routine_info(id_d_upper_to_bv, &
       'd_upper_to_bv', &
       [ character(len=error_message_length) :: 'n<1', 'bv%lbwmax < lbw', &
       'Size of a and bv not the same.' ] )

  type(routine_info), parameter :: info_f_d_upper_to_bv=routine_info(id_f_d_upper_to_bv, &
       'f_d_upper_to_bv', &
       [ character(len=error_message_length) :: 'Insufficient Upper Bandwidth in bv' ])

  type(routine_info), parameter :: info_f_d_general_bv=routine_info(id_f_d_general_bv, &
       'f_d_general_bv', &
       [ character(len=error_message_length) :: 'Insufficient Upper Bandwidth.' ])

  type(routine_info), parameter :: info_c_upper_to_bv=routine_info(id_c_upper_to_bv, &
       'c_upper_to_bv', &
       [ character(len=error_message_length) :: 'n<1', 'bv%lbwmax < lbw', &
       'Size of a and bv not the same.' ] )

  type(routine_info), parameter :: info_f_c_upper_to_bv=routine_info(id_f_c_upper_to_bv, &
       'f_c_upper_to_bv', &
       [ character(len=error_message_length) :: 'n<1', 'Insufficient Upper Bandwidth in bv' ])

  type(routine_info), parameter :: info_f_c_general_bv=routine_info(id_f_c_general_bv, &
       'f_c_general_bv', &
       [ character(len=error_message_length) :: 'Insufficient Upper Bandwidth.' ])
  
  !
  ! src/general/general_ub, 120s
  !
  integer(int32), parameter :: id_d_upper_to_ub=120
  integer(int32), parameter :: id_f_d_upper_to_ub=121
  integer(int32), parameter :: id_f_d_general_ub=122
  integer(int32), parameter :: id_c_upper_to_ub=123
  integer(int32), parameter :: id_f_c_upper_to_ub=124
  integer(int32), parameter :: id_f_c_general_ub=125

  type(routine_info), parameter :: info_d_upper_to_ub=routine_info(id_d_upper_to_ub, &
       'd_upper_to_ub', &
       [ character(len=error_message_length) :: 'n<1', 'ub%lbwmax < lbw', &
       'Size of a and ub not the same.' ] )

  type(routine_info), parameter :: info_f_d_upper_to_ub=routine_info(id_f_d_upper_to_ub, &
       'f_d_upper_to_ub', &
       [ character(len=error_message_length) :: 'Insufficient Upper Bandwidth in ub' ])

  type(routine_info), parameter :: info_f_d_general_ub=routine_info(id_f_d_general_ub, &
       'f_d_general_ub', &
       [ character(len=error_message_length) :: 'Insufficient Upper Bandwidth.' ])

  type(routine_info), parameter :: info_c_upper_to_ub=routine_info(id_c_upper_to_ub, &
       'c_upper_to_ub', &
       [ character(len=error_message_length) :: 'n<1', 'ub%lbwmax < lbw', &
       'Size of a and ub not the same.' ] )

  type(routine_info), parameter :: info_f_c_upper_to_ub=routine_info(id_f_c_upper_to_ub, &
       'f_c_upper_to_ub', &
       [ character(len=error_message_length) :: 'Insufficient Upper Bandwidth in ub' ])

  type(routine_info), parameter :: info_f_c_general_ub=routine_info(id_f_c_general_ub, &
       'f_c_general_ub', &
       [ character(len=error_message_length) :: 'Insufficient Upper Bandwidth.' ])

  ! /src/general/general_bt 130s
  integer(int32), parameter :: id_d_lower_to_bt=130
  integer(int32), parameter :: id_f_d_lower_to_bt=131
  integer(int32), parameter :: id_f_d_general_bt=132
  integer(int32), parameter :: id_c_lower_to_bt=133
  integer(int32), parameter :: id_f_c_lower_to_bt=134
  integer(int32), parameter :: id_f_c_general_bt=135


  type(routine_info), parameter :: info_d_lower_to_bt=routine_info(id_d_lower_to_bt, &
       'd_lower_to_bt', &
       [ character(len=error_message_length) :: 'n<1', 'bt%ubwmax < ubw', &
       'Size of a and bt not the same.' ] )

  type(routine_info), parameter :: info_f_d_lower_to_bt=routine_info(id_f_d_lower_to_bt, &
       'f_d_lower_to_bt', &
       [ character(len=error_message_length) :: 'Insufficient lower bandwidth in bt' ])

  type(routine_info), parameter :: info_f_d_general_bt=routine_info(id_f_d_general_bt, &
       'f_d_general_bt', &
       [ character(len=error_message_length) :: 'Insufficient lower bandwidth.' ])

  type(routine_info), parameter :: info_c_lower_to_bt=routine_info(id_c_lower_to_bt, &
       'c_lower_to_bt', &
       [ character(len=error_message_length) :: 'n<1', &
       'bt%ubwmax < ubw', 'Size of a and bt not the same.' ] )

  type(routine_info), parameter :: info_f_c_lower_to_bt=routine_info(id_f_c_lower_to_bt, &
       'f_c_lower_to_bt', &
       [ character(len=error_message_length) :: 'Insufficient lower bandwidth in bt' ])

  type(routine_info), parameter :: info_f_c_general_bt=routine_info(id_f_c_general_bt, &
       'f_c_general_bt', &
       [ character(len=error_message_length) :: 'Insufficient lower bandwidth.' ])


  ! /src/general/general_wb 140s
  integer(int32), parameter :: id_d_lower_to_wb=140
  integer(int32), parameter :: id_f_d_lower_to_wb=141
  integer(int32), parameter :: id_f_d_general_wb=142
  integer(int32), parameter :: id_c_lower_to_wb=143
  integer(int32), parameter :: id_f_c_lower_to_wb=144
  integer(int32), parameter :: id_f_c_general_wb=145
  
  type(routine_info), parameter :: info_d_lower_to_wb=routine_info(id_d_lower_to_wb, &
       'd_lower_to_wb', &
       [ character(len=error_message_length) :: 'n<1', 'bt%ubwmax < ubw', &
       'Size of a and bt not the same.' ] )

  type(routine_info), parameter :: info_f_d_lower_to_wb=routine_info(id_f_d_lower_to_wb, &
       'f_d_lower_to_wb', &
       [ character(len=error_message_length) :: 'Insufficient lower bandwidth in bt' ])

  type(routine_info), parameter :: info_f_d_general_wb=routine_info(id_f_d_general_wb, &
       'f_d_general_wb', &
       [ character(len=error_message_length) :: 'Insufficient lower bandwidth.' ])

  type(routine_info), parameter :: info_c_lower_to_wb=routine_info(id_c_lower_to_wb, &
       'c_lower_to_wb', &
       [ character(len=error_message_length) :: 'n<1', 'bt%ubwmax < ubw', &
       'Size of a and bt not the same.' ] )

  type(routine_info), parameter :: info_f_c_lower_to_wb=routine_info(id_f_c_lower_to_wb, &
       'f_c_lower_to_wb', &
       [ character(len=error_message_length) :: 'Insufficient lower bandwidth in bt' ])

  type(routine_info), parameter :: info_f_c_general_wb=routine_info(id_f_c_general_wb, &
       'f_c_general_wb', &
       [ character(len=error_message_length) :: 'Insufficient lower bandwidth.' ])

  ! src/general/general_ubt, 150s
  integer(int32), parameter :: id_d_general_to_ubt=150
  integer(int32), parameter :: id_f_d_general_to_ubt=151
  integer(int32), parameter :: id_f_d_general_ubt=152
  integer(int32), parameter :: id_c_general_to_ubt=153
  integer(int32), parameter :: id_f_c_general_to_ubt=154
  integer(int32), parameter :: id_f_c_general_ubt=155

  type(routine_info), parameter :: info_d_general_to_ubt=routine_info(id_d_general_to_ubt, &
       'd_general_to_ubt', &
       [ character(len=error_message_length) :: 'n<1', &
       'Size of a and ubt not the same.' ] )

  type(routine_info), parameter :: info_f_d_general_to_ubt=routine_info(id_f_d_general_to_ubt, &
       'f_d_general_to_ubt', &
       [ character(len=error_message_length) :: 'Insufficient lower bandwidth in ubt', &
       'Insufficient upper bandwidth in ubt' ])

  type(routine_info), parameter :: info_c_general_to_ubt=routine_info(id_c_general_to_ubt, &
       'c_general_to_ubt', &
       [ character(len=error_message_length) :: 'n<1', &
       'Size of a and ubt not the same.' ] )

  type(routine_info), parameter :: info_f_c_general_to_ubt=routine_info(id_f_c_general_to_ubt, &
       'f_c_general_to_ubt', &
       [ character(len=error_message_length) :: 'Insufficient lower bandwidth in ubt', &
       'Insufficient upper bandwidth in ubt' ])

  ! src/general/general_wbv, 160s
  integer(int32), parameter :: id_d_general_to_wbv=160
  integer(int32), parameter :: id_f_d_general_to_wbv=161
  integer(int32), parameter :: id_f_d_general_wbv=162
  integer(int32), parameter :: id_c_general_to_wbv=163
  integer(int32), parameter :: id_f_c_general_to_wbv=164
  integer(int32), parameter :: id_f_c_general_wbv=165

  type(routine_info), parameter :: info_d_general_to_wbv=routine_info(id_d_general_to_wbv, &
       'd_general_to_wbv', &
       [ character(len=error_message_length) :: 'n<1', &
       'Size of a and wbv not the same.' ] )

  type(routine_info), parameter :: info_f_d_general_to_wbv=routine_info(id_f_d_general_to_wbv, &
       'f_d_general_to_wbv', &
       [ character(len=error_message_length) :: 'Insufficient lower bandwidth in wbv', &
       'Insufficient upper bandwidth in wbv' ])

  type(routine_info), parameter :: info_c_general_to_wbv=routine_info(id_c_general_to_wbv, &
       'c_general_to_wbv', &
       [ character(len=error_message_length) :: 'n<1', &
       'Size of a and wbv not the same.' ] )

  type(routine_info), parameter :: info_f_c_general_to_wbv=routine_info(id_f_c_general_to_wbv, &
       'f_c_general_to_wbv', &
       [ character(len=error_message_length) :: 'Insufficient lower bandwidth in wbv', &
       'Insufficient upper bandwidth in wbv' ])

  ! src/orth/gs 170s
  integer(int32), parameter :: id_d_extend_gs_rows=170
  integer(int32), parameter :: id_c_extend_gs_rows=171
  integer(int32), parameter :: id_d_extend_gs_columns=172
  integer(int32), parameter :: id_c_extend_gs_columns=173

  type(routine_info), parameter :: info_d_extend_gs_rows=routine_info(id_d_extend_gs_rows, &
       'd_extend_gs_rows', &
       [ character(len=error_message_length) :: 'GS Orthogonalization Error' ])

  type(routine_info), parameter :: info_c_extend_gs_rows=routine_info(id_c_extend_gs_rows, &
       'c_extend_gs_rows', &
       [ character(len=error_message_length) :: 'GS Orthogonalization Error' ])

  type(routine_info), parameter :: info_d_extend_gs_columns=routine_info(id_d_extend_gs_columns, &
       'd_extend_gs_columns', &
       [ character(len=error_message_length) :: 'GS Orthogonalization Error' ])

  type(routine_info), parameter :: info_c_extend_gs_columns=routine_info(id_c_extend_gs_columns, &
       'c_extend_gs_columns', &
       [ character(len=error_message_length) :: 'GS Orthogonalization Error' ])

  ! src/orth/nullvec 180s
  integer(int32), parameter :: id_f_d_lower_left_nullvec=180
  integer(int32), parameter :: id_f_c_lower_left_nullvec=181
  integer(int32), parameter :: id_f_d_lower_right_nullvec=182
  integer(int32), parameter :: id_f_c_lower_right_nullvec=183

  type(routine_info), parameter :: info_f_d_lower_left_nullvec= &
       routine_info(id_f_d_lower_left_nullvec, 'f_d_lower_left_nullvec', &
       [ character(len=error_message_length) :: 'Failure to find a null vector.' ])
  type(routine_info), parameter :: info_f_d_lower_right_nullvec= &
       routine_info(id_f_d_lower_right_nullvec, 'f_d_lower_right_nullvec', &
       [ character(len=error_message_length) :: 'Failure to find a null vector.' ])
  type(routine_info), parameter :: info_f_c_lower_left_nullvec= &
       routine_info(id_f_c_lower_left_nullvec, 'f_c_lower_left_nullvec', &
       [ character(len=error_message_length) :: 'Failure to find a null vector.' ])
  type(routine_info), parameter :: info_f_c_lower_right_nullvec= &
       routine_info(id_f_c_lower_right_nullvec, 'f_c_lower_right_nullvec', &
       [ character(len=error_message_length) :: 'Failure to find a null vector.' ])

  ! src/qr_factorization/qr_factorization 200s
  integer(int32), parameter :: id_d_qr_bv_to_ub=200
  integer(int32), parameter :: id_f_d_qr_bv_to_ub=201
  integer(int32), parameter :: id_c_qr_bv_to_ub=202
  integer(int32), parameter :: id_f_c_qr_bv_to_ub=203

  type(routine_info), parameter :: info_d_qr_bv_to_ub=routine_info(id_d_qr_bv_to_ub, &
       'd_qr_bv_to_ub', &
       [ character(len=error_message_length) :: 'ub%n /= bv%n or sw%n /= ub%n', &
       'Not enough stroage for sweeps' ])

  type(routine_info), parameter :: info_c_qr_bv_to_ub=routine_info(id_c_qr_bv_to_ub, &
       'c_qr_bv_to_ub', &
       [ character(len=error_message_length) :: 'ub%n /= bv%n', 'bv%lbw <= 0', &
       'dim. of cs or ss /= n' ])

  ! src/solve/back_solve 210s and 220s
  integer(int32), parameter :: id_d_back_solve_ub=210
  integer(int32), parameter :: id_f_d_back_solve_ub=211
  integer(int32), parameter :: id_c_back_solve_ub=212
  integer(int32), parameter :: id_f_c_back_solve_ub=213
  integer(int32), parameter :: id_d_v_back_solve_ub=214
  integer(int32), parameter :: id_f_d_v_back_solve_ub=215
  integer(int32), parameter :: id_c_v_back_solve_ub=216
  integer(int32), parameter :: id_f_c_v_back_solve_ub=217

  type(routine_info), parameter :: info_d_back_solve_ub= &
       routine_info(id_d_back_solve_ub, &
       'd_back_solve_ub', &
       [ character(len=error_message_length) :: 'ub%n /= size(c,1)', 'ub%lbw /= 0', 'n < 1', &
       'x is not the same size as c.' ])

  type(routine_info), parameter :: info_c_back_solve_ub= &
       routine_info(id_c_back_solve_ub, &
       'd_back_solve_ub', &
       [ character(len=error_message_length) :: 'ub%n /= size(c,1)', 'ub%lbw /= 0', 'n < 1', &
       'x is not the same size as c.' ])

  type(routine_info), parameter :: info_d_v_back_solve_ub= &
       routine_info(id_d_v_back_solve_ub, &
       'd_v_back_solve_ub', &
       [ character(len=error_message_length) :: 'ub%n /= size(c)', 'ub%lbw /= 0', 'n < 1', &
       'x is not the same size as c.' ])

  type(routine_info), parameter :: info_c_v_back_solve_ub= &
       routine_info(id_c_v_back_solve_ub, &
       'd_back_solve_ub', &
       [ character(len=error_message_length) :: 'ub%n /= size(c)', 'ub%lbw /= 0', 'n < 1', &
       'x is not the same size as c.' ])

  integer(int32), parameter :: id_d_forward_solve_bv=218
  integer(int32), parameter :: id_f_d_forward_solve_bv=219
  integer(int32), parameter :: id_d_v_forward_solve_bv=220
  integer(int32), parameter :: id_f_d_v_forward_solve_bv=221
  integer(int32), parameter :: id_c_forward_solve_bv=222
  integer(int32), parameter :: id_f_c_forward_solve_bv=223
  integer(int32), parameter :: id_c_v_forward_solve_bv=224
  integer(int32), parameter :: id_f_c_v_forward_solve_bv=225

  type(routine_info), parameter :: info_d_forward_solve_bv= &
       routine_info(id_d_forward_solve_bv, &
       'd_forward_solve_bv', &
       [ character(len=error_message_length) :: 'bv%n /= size(c,2)', 'bv%lbw /= 0', 'n < 1', &
       'x is not the same size as c.' ])

  type(routine_info), parameter :: info_c_forward_solve_bv= &
       routine_info(id_c_forward_solve_bv, &
       'd_forward_solve_bv', &
       [ character(len=error_message_length) :: 'bv%n /= size(c,2)', 'bv%lbw /= 0', 'n < 1', &
       'x is not the same size as c.' ])

  type(routine_info), parameter :: info_d_v_forward_solve_bv= &
       routine_info(id_d_v_forward_solve_bv, &
       'd_v_forward_solve_bv', &
       [ character(len=error_message_length) :: 'bv%n /= size(c,2)', 'bv%lbw /= 0', 'n < 1', &
       'x is not the same size as c.' ])

  type(routine_info), parameter :: info_c_v_forward_solve_bv= &
       routine_info(id_c_v_forward_solve_bv, &
       'd_forward_solve_bv', &
       [ character(len=error_message_length) :: 'bv%n /= size(c,1)', 'bv%lbw /= 0', 'n < 1', &
       'x is not the same size as c.' ])

  ! src/row_compress 260s
  integer(int32), parameter :: id_d_row_compress=260
  integer(int32), parameter :: id_f_d_row_compress=261
  integer(int32), parameter :: id_c_row_compress=262
  integer(int32), parameter :: id_f_c_row_compress=263

  type(routine_info), parameter :: info_d_row_compress=routine_info(id_d_row_compress, &
       'd_row_compress', &
       [ character(len=error_message_length) :: 'ub%n /= bv%n or sw%n /= ub%n', &
       'Not enough stroage for each sweeps', 'Not enough storage for number of sweeps', &
       'Insufficient storage in ubt.', 'Insufficient storage in bv.'])

  type(routine_info), parameter :: info_c_row_compress=routine_info(id_c_row_compress, &
       'c_row_compress', &
       [ character(len=error_message_length) :: 'ub%n /= bv%n or sw%n /= ub%n', &
       'Not enough stroage for each sweeps', 'Not enough storage for number of sweeps', &
       'Insufficient storage in ubt.', 'Insufficient storage in bv.'])

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

  subroutine initialize_errors
    if (.not. allocated(info_index)) then
       allocate(info_index(max_error_code))
       info_index = info_empty
       ! assemble.f90
       info_index(info_d_ub_to_upper%routine_id)=info_d_ub_to_upper
       info_index(info_c_ub_to_upper%routine_id)=info_c_ub_to_upper

       info_index(info_d_bt_to_lower%routine_id)=info_d_bt_to_lower
       info_index(info_c_bt_to_lower%routine_id)=info_c_bt_to_lower

       info_index(info_d_ubt_to_general%routine_id)=info_d_ubt_to_general
       info_index(info_c_ubt_to_general%routine_id)=info_c_ubt_to_general

       info_index(info_d_bv_to_upper%routine_id)=info_d_bv_to_upper
       info_index(info_c_bv_to_upper%routine_id)=info_c_bv_to_upper

       info_index(info_d_wb_to_lower%routine_id)=info_d_wb_to_lower
       info_index(info_c_wb_to_lower%routine_id)=info_c_wb_to_lower

       info_index(info_d_wbv_to_general%routine_id)=info_d_wbv_to_general
       info_index(info_c_wbv_to_general%routine_id)=info_c_wbv_to_general

       ! convert_bv_to_ub
       info_index(info_d_convert_bv_to_ub%routine_id)=info_d_convert_bv_to_ub
       info_index(info_c_convert_bv_to_ub%routine_id)=info_c_convert_bv_to_ub       

       ! convert_ub_to_bv
       info_index(info_d_convert_ub_to_bv%routine_id)=info_d_convert_ub_to_bv
       info_index(info_c_convert_ub_to_bv%routine_id)=info_c_convert_ub_to_bv
       
       ! convert_wb_to_bt
       info_index(info_d_convert_wb_to_bt%routine_id)=info_d_convert_wb_to_bt
       info_index(info_c_convert_wb_to_bt%routine_id)=info_c_convert_wb_to_bt       

       ! convert_bt_to_wb
       info_index(info_d_convert_bt_to_wb%routine_id)=info_d_convert_bt_to_wb
       info_index(info_c_convert_bt_to_wb%routine_id)=info_c_convert_bt_to_wb

       ! convert_wbv_to_ubt
       info_index(info_d_convert_wbv_to_ubt%routine_id)=info_d_convert_wbv_to_ubt
       info_index(info_c_convert_wbv_to_ubt%routine_id)=info_c_convert_wbv_to_ubt

       ! convert_ubt_to_wbv
       info_index(info_d_convert_ubt_to_wbv%routine_id)=info_d_convert_ubt_to_wbv
       info_index(info_c_convert_ubt_to_wbv%routine_id)=info_c_convert_ubt_to_wbv
       
       ! general_bv

       info_index(info_d_upper_to_bv%routine_id)=info_d_upper_to_bv
       info_index(info_f_d_upper_to_bv%routine_id)=info_f_d_upper_to_bv
       info_index(info_f_d_general_bv%routine_id)=info_f_d_general_bv
       info_index(info_c_upper_to_bv%routine_id)=info_c_upper_to_bv
       info_index(info_f_c_upper_to_bv%routine_id)=info_f_c_upper_to_bv
       info_index(info_f_c_general_bv%routine_id)=info_f_c_general_bv

       ! general_ub

       info_index(info_d_upper_to_ub%routine_id)=info_d_upper_to_ub
       info_index(info_f_d_upper_to_ub%routine_id)=info_f_d_upper_to_ub
       info_index(info_f_d_general_ub%routine_id)=info_f_d_general_ub
       info_index(info_c_upper_to_ub%routine_id)=info_c_upper_to_ub
       info_index(info_f_c_upper_to_ub%routine_id)=info_f_c_upper_to_ub
       info_index(info_f_c_general_ub%routine_id)=info_f_c_general_ub

       ! general_bt

       info_index(info_d_lower_to_bt%routine_id)=info_d_lower_to_bt
       info_index(info_f_d_lower_to_bt%routine_id)=info_f_d_lower_to_bt
       info_index(info_f_d_general_bt%routine_id)=info_f_d_general_bt
       info_index(info_c_lower_to_bt%routine_id)=info_c_lower_to_bt
       info_index(info_f_c_lower_to_bt%routine_id)=info_f_c_lower_to_bt
       info_index(info_f_c_general_bt%routine_id)=info_f_c_general_bt

       ! general_wb

       info_index(info_d_lower_to_wb%routine_id)=info_d_lower_to_wb
       info_index(info_f_d_lower_to_wb%routine_id)=info_f_d_lower_to_wb
       info_index(info_f_d_general_wb%routine_id)=info_f_d_general_wb
       info_index(info_c_lower_to_wb%routine_id)=info_c_lower_to_wb
       info_index(info_f_c_lower_to_wb%routine_id)=info_f_c_lower_to_wb
       info_index(info_f_c_general_wb%routine_id)=info_f_c_general_wb

       ! general_ubt

       info_index(info_d_general_to_ubt%routine_id)=info_d_general_to_ubt
       info_index(info_f_d_general_to_ubt%routine_id)=info_f_d_general_to_ubt
       info_index(info_c_general_to_ubt%routine_id)=info_c_general_to_ubt
       info_index(info_f_c_general_to_ubt%routine_id)=info_f_c_general_to_ubt

       ! general_wbv

       info_index(info_d_general_to_wbv%routine_id)=info_d_general_to_wbv
       info_index(info_f_d_general_to_wbv%routine_id)=info_f_d_general_to_wbv
       info_index(info_c_general_to_wbv%routine_id)=info_c_general_to_wbv
       info_index(info_f_c_general_to_wbv%routine_id)=info_f_c_general_to_wbv

       ! gs
       info_index(info_d_extend_gs_rows%routine_id)=info_d_extend_gs_rows
       info_index(info_c_extend_gs_rows%routine_id)=info_c_extend_gs_rows
       info_index(info_d_extend_gs_columns%routine_id)=info_d_extend_gs_columns
       info_index(info_c_extend_gs_columns%routine_id)=info_c_extend_gs_columns
       
       ! nullvec
       info_index(info_f_d_lower_left_nullvec%routine_id)=info_f_d_lower_left_nullvec
       info_index(info_f_c_lower_left_nullvec%routine_id)=info_f_c_lower_left_nullvec
       info_index(info_f_d_lower_right_nullvec%routine_id)=info_f_d_lower_right_nullvec
       info_index(info_f_c_lower_right_nullvec%routine_id)=info_f_c_lower_right_nullvec

       ! qr_factorization

       info_index(info_d_qr_bv_to_ub%routine_id)=info_d_qr_bv_to_ub
       info_index(info_c_qr_bv_to_ub%routine_id)=info_c_qr_bv_to_ub

       ! back_solve
       info_index(info_d_back_solve_ub%routine_id)=info_d_back_solve_ub
       info_index(info_c_back_solve_ub%routine_id)=info_c_back_solve_ub
       info_index(info_d_v_back_solve_ub%routine_id)=info_d_v_back_solve_ub
       info_index(info_c_v_back_solve_ub%routine_id)=info_c_v_back_solve_ub

       info_index(info_d_forward_solve_bv%routine_id)=info_d_forward_solve_bv
       info_index(info_c_forward_solve_bv%routine_id)=info_c_forward_solve_bv
       info_index(info_d_v_forward_solve_bv%routine_id)=info_d_v_forward_solve_bv
       info_index(info_c_v_forward_solve_bv%routine_id)=info_c_v_forward_solve_bv

       ! row_compress

       info_index(info_d_row_compress%routine_id)=info_d_row_compress
       info_index(info_c_row_compress%routine_id)=info_c_row_compress

    end if
  end subroutine initialize_errors

end module mod_error_id
