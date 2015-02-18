module mod_error_id
  use mod_prec
  use, intrinsic :: iso_fortran_env, only : error_unit
  implicit none
  integer(kind=int32), parameter :: routine_name_length=30
  integer(kind=int32), parameter :: max_errors=15
  integer(kind=int32), parameter :: error_message_length=50
  integer(kind=int32), parameter :: max_routines=10

  integer(kind=int32), parameter :: max_error_code=5000

  public

  type routine_info
     integer(kind=int32) :: routine_id
     character(len=routine_name_length) :: routine_name
     character(len=error_message_length), dimension(max_errors) :: error_messages=''
  end type routine_info

  type error_info
     integer(kind=int32) :: code=0, rix=1
     integer(kind=int32), dimension(max_routines) :: routines=0
  end type error_info

  type(routine_info), parameter :: info_empty=routine_info(0, &
       '', [ character(len=error_message_length) :: '' ] )

  type(routine_info), dimension(:), allocatable :: info_index

  ! src/assemble/assemble.f90, 0000s
  integer(int32), parameter :: mod_id_assemble=0
  integer(int32), parameter :: id_d_ub_to_general=mod_id_assemble + 1
  integer(int32), parameter :: id_f_d_ub_to_general=mod_id_assemble + 0002
  integer(int32), parameter :: id_d_general_of_ub=mod_id_assemble + 0003
  integer(int32), parameter :: id_c_ub_to_general=mod_id_assemble + 0004
  integer(int32), parameter :: id_f_c_ub_to_general=mod_id_assemble + 0005
  integer(int32), parameter :: id_c_general_of_ub=mod_id_assemble + 0006

  integer(int32), parameter :: id_d_bt_to_general=mod_id_assemble + 0007
  integer(int32), parameter :: id_f_d_bt_to_general=mod_id_assemble + 0008
  integer(int32), parameter :: id_d_general_of_bt=mod_id_assemble + 0009
  integer(int32), parameter :: id_c_bt_to_general=mod_id_assemble + 0010
  integer(int32), parameter :: id_f_c_bt_to_general=mod_id_assemble + 0011
  integer(int32), parameter :: id_c_general_of_bt=mod_id_assemble + 0012
  
  integer(int32), parameter :: id_d_ubt_to_general=mod_id_assemble + 0013
  integer(int32), parameter :: id_f_d_ubt_to_general=mod_id_assemble + 0014
  integer(int32), parameter :: id_d_general_of_ubt=mod_id_assemble + 0015
  integer(int32), parameter :: id_c_ubt_to_general=mod_id_assemble + 0016
  integer(int32), parameter :: id_f_c_ubt_to_general=mod_id_assemble + 0017
  integer(int32), parameter :: id_c_general_of_ubt=mod_id_assemble + 0018

  integer(int32), parameter :: id_d_bv_to_general=mod_id_assemble + 0019
  integer(int32), parameter :: id_f_d_bv_to_general=mod_id_assemble + 0020
  integer(int32), parameter :: id_d_general_of_bv=mod_id_assemble + 0021
  integer(int32), parameter :: id_c_bv_to_general=mod_id_assemble + 0022
  integer(int32), parameter :: id_f_c_bv_to_general=mod_id_assemble + 0023
  integer(int32), parameter :: id_c_general_of_bv=mod_id_assemble + 0024

  integer(int32), parameter :: id_d_wb_to_general=mod_id_assemble + 0025
  integer(int32), parameter :: id_f_d_wb_to_general=mod_id_assemble + 0026
  integer(int32), parameter :: id_d_general_of_wb=mod_id_assemble + 0027
  integer(int32), parameter :: id_c_wb_to_general=mod_id_assemble + 0028
  integer(int32), parameter :: id_f_c_wb_to_general=mod_id_assemble + 0029
  integer(int32), parameter :: id_c_general_of_wb=mod_id_assemble + 0030
  
  integer(int32), parameter :: id_d_wbv_to_general=mod_id_assemble + 0031
  integer(int32), parameter :: id_f_d_wbv_to_general=mod_id_assemble + 0032
  integer(int32), parameter :: id_d_general_of_wbv=mod_id_assemble + 0033
  integer(int32), parameter :: id_c_wbv_to_general=mod_id_assemble + 0034
  integer(int32), parameter :: id_f_c_wbv_to_general=mod_id_assemble + 0035
  integer(int32), parameter :: id_c_general_of_wbv=mod_id_assemble + 0036

  type(routine_info), parameter :: info_d_ub_to_general=routine_info(id_d_ub_to_general, &
       'd_ub_to_general', [ character(len=error_message_length) :: 'Size error in A.' ] )

  type(routine_info), parameter :: info_c_ub_to_general=routine_info(id_c_ub_to_general, &
       'c_ub_to_general', [ character(len=error_message_length) :: 'Size error in A.' ] )

  type(routine_info), parameter :: info_d_bt_to_general=routine_info(id_d_bt_to_general, &
       'd_bt_to_general', [ character(len=error_message_length) :: 'Size error in A.' ] )

  type(routine_info), parameter :: info_c_bt_to_general=routine_info(id_c_bt_to_general, &
       'c_bt_to_general', [ character(len=error_message_length) :: 'Size error in A.' ] )

  type(routine_info), parameter :: info_d_ubt_to_general=routine_info(id_d_ubt_to_general, &
       'd_ubt_to_general', [ character(len=error_message_length) :: 'Size error in A.' ] )

  type(routine_info), parameter :: info_c_ubt_to_general=routine_info(id_c_ubt_to_general, &
       'c_ubt_to_general', [ character(len=error_message_length) :: 'Size error in A.' ] )

  type(routine_info), parameter :: info_d_bv_to_general=routine_info(id_d_bv_to_general, &
       'd_bv_to_general', [ character(len=error_message_length) :: 'Size error in A.' ] )
  
  type(routine_info), parameter :: info_c_bv_to_general=routine_info(id_c_bv_to_general, &
       'c_bv_to_general', [ character(len=error_message_length) :: 'Size error in A.' ] )

  type(routine_info), parameter :: info_d_wb_to_general=routine_info(id_d_wb_to_general, &
       'd_wb_to_general', [ character(len=error_message_length) :: 'Size error in A.' ] )
  
  type(routine_info), parameter :: info_c_wb_to_general=routine_info(id_c_wb_to_general, &
       'c_wb_to_general', [ character(len=error_message_length) :: 'Size error in A.' ] )

  type(routine_info), parameter :: info_d_wbv_to_general=routine_info(id_d_wbv_to_general, &
       'd_wbv_to_general', [ character(len=error_message_length) :: 'Size error in A.' ] )
  
  type(routine_info), parameter :: info_c_wbv_to_general=routine_info(id_c_wbv_to_general, &
       'c_wbv_to_general', [ character(len=error_message_length) :: 'Size error in A.' ] )

  ! src/convert/convert_bv_to_ub, 0100s
  integer(int32), parameter :: mod_id_convert_bv_to_ub=100
  integer(int32), parameter :: id_d_convert_bv_to_ub=mod_id_convert_bv_to_ub+0
  integer(int32), parameter :: id_f_d_convert_bv_to_ub=mod_id_convert_bv_to_ub+2
  integer(int32), parameter :: id_c_convert_bv_to_ub=mod_id_convert_bv_to_ub+3
  integer(int32), parameter :: id_f_c_convert_bv_to_ub=mod_id_convert_bv_to_ub+4

  type(routine_info), parameter :: info_d_convert_bv_to_ub=routine_info(id_d_convert_bv_to_ub, &
       'd_convert_bv_to_ub', &
       [ character(len=error_message_length) :: 'n<1', 'Insufficient storage in bv.', &
       'Insufficient Storage in ub.', 'ub%n /= bv%n' ] )

  type(routine_info), parameter :: info_c_convert_bv_to_ub=routine_info(id_c_convert_bv_to_ub, &
       'c_convert_bv_to_ub', &
       [ character(len=error_message_length) :: 'n<1', 'Insufficient storage in bv.', &
       'Insufficient Storage in ub.', 'ub%n /= bv%n' ] )

  ! src/convert/convert_ub_to_bv, 200s
  integer(int32), parameter :: mod_id_convert_ub_to_bv=200
  integer(int32), parameter :: id_d_convert_ub_to_bv=mod_id_convert_ub_to_bv + 0
  integer(int32), parameter :: id_f_d_convert_ub_to_bv=mod_id_convert_ub_to_bv + 1
  integer(int32), parameter :: id_c_convert_ub_to_bv=mod_id_convert_ub_to_bv + 2
  integer(int32), parameter :: id_f_c_convert_ub_to_bv=mod_id_convert_ub_to_bv + 3

  type(routine_info), parameter :: info_d_convert_ub_to_bv=routine_info(id_d_convert_ub_to_bv, &
       'd_convert_ub_to_bv', &
       [ character(len=error_message_length) :: 'n<1', 'Insufficient storage in ub', &
       'Insufficient storage in bv.', 'ub%n /= bv%n' ] )

  type(routine_info), parameter :: info_c_convert_ub_to_bv=routine_info(id_c_convert_ub_to_bv, &
       'c_convert_ub_to_bv', &
       [ character(len=error_message_length) :: 'Insufficient storage in ub', &
       'Insufficient storage in bv.', 'ub%n /= bv%n' ] )

  ! src/convert/convert_wb_to_bt, 300
  integer(int32), parameter :: mod_id_convert_wb_to_bt=300
  integer(int32), parameter :: id_d_convert_wb_to_bt=mod_id_convert_wb_to_bt + 0
  integer(int32), parameter :: id_f_d_convert_wb_to_bt=mod_id_convert_wb_to_bt + 1
  integer(int32), parameter :: id_c_convert_wb_to_bt=mod_id_convert_wb_to_bt + 2
  integer(int32), parameter :: id_f_c_convert_wb_to_bt=mod_id_convert_wb_to_bt + 3

  type(routine_info), parameter :: info_d_convert_wb_to_bt=routine_info(id_d_convert_wb_to_bt, &
       'd_convert_wb_to_bt', &
       [ character(len=error_message_length) :: 'n<1', 'Insufficient storage in wb.', &
       'Insufficient Storage in bt.', 'bt%n /= wb%n' ] )

  type(routine_info), parameter :: info_c_convert_wb_to_bt=routine_info(id_c_convert_wb_to_bt, &
       'c_convert_wb_to_bt', &
       [ character(len=error_message_length) :: 'n<1', 'Insufficient storage in wb.', &
       'Insufficient Storage in bt.', 'bt%n /= wb%n' ] )

  ! src/convert/convert_bt_to_wb, 400
  integer(int32), parameter :: mod_id_convert_bt_to_wb=400
  integer(int32), parameter :: id_d_convert_bt_to_wb=mod_id_convert_bt_to_wb + 0
  integer(int32), parameter :: id_f_d_convert_bt_to_wb=mod_id_convert_bt_to_wb + 1
  integer(int32), parameter :: id_c_convert_bt_to_wb=mod_id_convert_bt_to_wb + 2
  integer(int32), parameter :: id_f_c_convert_bt_to_wb=mod_id_convert_bt_to_wb + 3

  type(routine_info), parameter :: info_d_convert_bt_to_wb=routine_info(id_d_convert_bt_to_wb, &
       'd_convert_bt_to_wb', &
       [ character(len=error_message_length) :: 'n<1', 'Insufficient storage in bt.', &
       'Insufficient Storage in wb.', 'wb%n /= bt%n' ] )

  type(routine_info), parameter :: info_c_convert_bt_to_wb=routine_info(id_c_convert_bt_to_wb, &
       'c_convert_bt_to_wb', &
       [ character(len=error_message_length) :: 'n<1', 'Insufficient storage in bt.', &
       'Insufficient Storage in wb.', 'wb%n /= bt%n' ] )

  ! src/convert/convert_wbv_to_ubt, 500
  integer(int32), parameter :: mod_id_convert_wbv_to_ubt=500
  integer(int32), parameter :: id_d_convert_wbv_to_ubt=mod_id_convert_wbv_to_ubt + 0
  integer(int32), parameter :: id_f_d_convert_wbv_to_ubt=mod_id_convert_wbv_to_ubt + 1
  integer(int32), parameter :: id_c_convert_wbv_to_ubt=mod_id_convert_wbv_to_ubt + 2
  integer(int32), parameter :: id_f_c_convert_wbv_to_ubt=mod_id_convert_wbv_to_ubt + 3

  type(routine_info), parameter :: info_d_convert_wbv_to_ubt=routine_info(id_d_convert_wbv_to_ubt, &
       'd_convert_wbv_to_ubt', &
       [ character(len=error_message_length) :: 'n<1', 'Insufficient storage in wbv', &
       'Insufficient storage in ubt.', 'wbv%n /= ubt%n' ] )

  type(routine_info), parameter :: info_c_convert_wbv_to_ubt=routine_info(id_c_convert_wbv_to_ubt, &
       'c_convert_wbv_to_ubt', &
       [ character(len=error_message_length) :: 'Insufficient storage in wbv', &
       'Insufficient storage in ubt.', 'wbv%n /= ubt%n' ] )

  ! src/convert/convert_ubt_to_wbv, 600
  integer(int32), parameter :: mod_id_convert_ubt_to_wbv=600
  integer(int32), parameter :: id_d_convert_ubt_to_wbv=mod_id_convert_ubt_to_wbv + 0
  integer(int32), parameter :: id_f_d_convert_ubt_to_wbv=mod_id_convert_ubt_to_wbv + 1
  integer(int32), parameter :: id_c_convert_ubt_to_wbv=mod_id_convert_ubt_to_wbv + 2
  integer(int32), parameter :: id_f_c_convert_ubt_to_wbv=mod_id_convert_ubt_to_wbv + 3

  type(routine_info), parameter :: info_d_convert_ubt_to_wbv=routine_info(id_d_convert_ubt_to_wbv, &
       'd_convert_ubt_to_wbv', &
       [ character(len=error_message_length) :: 'n<1', 'Insufficient storage in ubt', &
       'Insufficient storage in wbv.', 'ubt%n /= wbv%n' ] )

  type(routine_info), parameter :: info_c_convert_ubt_to_wbv=routine_info(id_c_convert_ubt_to_wbv, &
       'c_convert_ubt_to_wbv', &
       [ character(len=error_message_length) :: 'Insufficient storage in ubt', &
       'Insufficient storage in wbv.', 'ubt%n /= wbv%n' ] )

  ! src/general/general_bv, 700
  integer(int32), parameter :: mod_id_general_bv=700
  integer(int32), parameter :: id_d_general_to_bv=mod_id_general_bv + 0
  integer(int32), parameter :: id_f_d_general_to_bv=mod_id_general_bv + 1
  integer(int32), parameter :: id_f_d_general_bv=mod_id_general_bv + 2
  integer(int32), parameter :: id_c_general_to_bv=mod_id_general_bv + 3
  integer(int32), parameter :: id_f_c_general_to_bv=mod_id_general_bv + 4
  integer(int32), parameter :: id_f_c_general_bv=mod_id_general_bv + 5

  type(routine_info), parameter :: info_d_general_to_bv=routine_info(id_d_general_to_bv, &
       'd_general_to_bv', &
       [ character(len=error_message_length) :: 'n<1', 'bv%lbwmax < lbw', &
       'Size of a and bv not the same.' ] )

  type(routine_info), parameter :: info_f_d_general_to_bv=routine_info(id_f_d_general_to_bv, &
       'f_d_general_to_bv', &
       [ character(len=error_message_length) :: 'Insufficient General Bandwidth in bv' ])

  type(routine_info), parameter :: info_f_d_general_bv=routine_info(id_f_d_general_bv, &
       'f_d_general_bv', &
       [ character(len=error_message_length) :: 'Insufficient General Bandwidth.' ])

  type(routine_info), parameter :: info_c_general_to_bv=routine_info(id_c_general_to_bv, &
       'c_general_to_bv', &
       [ character(len=error_message_length) :: 'n<1', 'bv%lbwmax < lbw', &
       'Size of a and bv not the same.' ] )

  type(routine_info), parameter :: info_f_c_general_to_bv=routine_info(id_f_c_general_to_bv, &
       'f_c_general_to_bv', &
       [ character(len=error_message_length) :: 'n<1', 'Insufficient General Bandwidth in bv' ])

  type(routine_info), parameter :: info_f_c_general_bv=routine_info(id_f_c_general_bv, &
       'f_c_general_bv', &
       [ character(len=error_message_length) :: 'Insufficient General Bandwidth.' ])
  
  !
  ! src/general/general_ub, 800
  !
  integer(int32), parameter :: mod_id_general_ub=800
  integer(int32), parameter :: id_d_general_to_ub=mod_id_general_ub + 0
  integer(int32), parameter :: id_f_d_general_to_ub=mod_id_general_ub + 1
  integer(int32), parameter :: id_f_d_general_ub=mod_id_general_ub + 2
  integer(int32), parameter :: id_c_general_to_ub=mod_id_general_ub + 3
  integer(int32), parameter :: id_f_c_general_to_ub=mod_id_general_ub + 4
  integer(int32), parameter :: id_f_c_general_ub=mod_id_general_ub + 5

  type(routine_info), parameter :: info_d_general_to_ub=routine_info(id_d_general_to_ub, &
       'd_general_to_ub', &
       [ character(len=error_message_length) :: 'n<1', 'ub%lbwmax < lbw', &
       'Size of a and ub not the same.' ] )

  type(routine_info), parameter :: info_f_d_general_to_ub=routine_info(id_f_d_general_to_ub, &
       'f_d_general_to_ub', &
       [ character(len=error_message_length) :: 'Insufficient General Bandwidth in ub' ])

  type(routine_info), parameter :: info_f_d_general_ub=routine_info(id_f_d_general_ub, &
       'f_d_general_ub', &
       [ character(len=error_message_length) :: 'Insufficient General Bandwidth.' ])

  type(routine_info), parameter :: info_c_general_to_ub=routine_info(id_c_general_to_ub, &
       'c_general_to_ub', &
       [ character(len=error_message_length) :: 'n<1', 'ub%lbwmax < lbw', &
       'Size of a and ub not the same.' ] )

  type(routine_info), parameter :: info_f_c_general_to_ub=routine_info(id_f_c_general_to_ub, &
       'f_c_general_to_ub', &
       [ character(len=error_message_length) :: 'Insufficient General Bandwidth in ub' ])

  type(routine_info), parameter :: info_f_c_general_ub=routine_info(id_f_c_general_ub, &
       'f_c_general_ub', &
       [ character(len=error_message_length) :: 'Insufficient General Bandwidth.' ])

  ! /src/general/general_bt 900
  integer(int32), parameter :: mod_id_general_bt=900
  integer(int32), parameter :: id_d_general_to_bt=mod_id_general_bt + 0
  integer(int32), parameter :: id_f_d_general_to_bt=mod_id_general_bt + 1
  integer(int32), parameter :: id_f_d_general_bt=mod_id_general_bt + 2
  integer(int32), parameter :: id_c_general_to_bt=mod_id_general_bt + 3
  integer(int32), parameter :: id_f_c_general_to_bt=mod_id_general_bt + 4
  integer(int32), parameter :: id_f_c_general_bt=mod_id_general_bt + 5


  type(routine_info), parameter :: info_d_general_to_bt=routine_info(id_d_general_to_bt, &
       'd_general_to_bt', &
       [ character(len=error_message_length) :: 'n<1', 'bt%ubwmax < ubw', &
       'Size of a and bt not the same.' ] )

  type(routine_info), parameter :: info_f_d_general_to_bt=routine_info(id_f_d_general_to_bt, &
       'f_d_general_to_bt', &
       [ character(len=error_message_length) :: 'Insufficient general bandwidth in bt' ])

  type(routine_info), parameter :: info_f_d_general_bt=routine_info(id_f_d_general_bt, &
       'f_d_general_bt', &
       [ character(len=error_message_length) :: 'Insufficient general bandwidth.' ])

  type(routine_info), parameter :: info_c_general_to_bt=routine_info(id_c_general_to_bt, &
       'c_general_to_bt', &
       [ character(len=error_message_length) :: 'n<1', &
       'bt%ubwmax < ubw', 'Size of a and bt not the same.' ] )

  type(routine_info), parameter :: info_f_c_general_to_bt=routine_info(id_f_c_general_to_bt, &
       'f_c_general_to_bt', &
       [ character(len=error_message_length) :: 'Insufficient general bandwidth in bt' ])

  type(routine_info), parameter :: info_f_c_general_bt=routine_info(id_f_c_general_bt, &
       'f_c_general_bt', &
       [ character(len=error_message_length) :: 'Insufficient general bandwidth.' ])


  ! /src/general/general_wb 140s
  integer(int32), parameter :: mod_id_general_wb=1000
  integer(int32), parameter :: id_d_general_to_wb=mod_id_general_wb + 0
  integer(int32), parameter :: id_f_d_general_to_wb=mod_id_general_wb + 1
  integer(int32), parameter :: id_f_d_general_wb=mod_id_general_wb + 2
  integer(int32), parameter :: id_c_general_to_wb=mod_id_general_wb + 3
  integer(int32), parameter :: id_f_c_general_to_wb=mod_id_general_wb + 4
  integer(int32), parameter :: id_f_c_general_wb=mod_id_general_wb + 5
  
  type(routine_info), parameter :: info_d_general_to_wb=routine_info(id_d_general_to_wb, &
       'd_general_to_wb', &
       [ character(len=error_message_length) :: 'n<1', 'bt%ubwmax < ubw', &
       'Size of a and bt not the same.' ] )

  type(routine_info), parameter :: info_f_d_general_to_wb=routine_info(id_f_d_general_to_wb, &
       'f_d_general_to_wb', &
       [ character(len=error_message_length) :: 'Insufficient general bandwidth in bt' ])

  type(routine_info), parameter :: info_f_d_general_wb=routine_info(id_f_d_general_wb, &
       'f_d_general_wb', &
       [ character(len=error_message_length) :: 'Insufficient general bandwidth.' ])

  type(routine_info), parameter :: info_c_general_to_wb=routine_info(id_c_general_to_wb, &
       'c_general_to_wb', &
       [ character(len=error_message_length) :: 'n<1', 'bt%ubwmax < ubw', &
       'Size of a and bt not the same.' ] )

  type(routine_info), parameter :: info_f_c_general_to_wb=routine_info(id_f_c_general_to_wb, &
       'f_c_general_to_wb', &
       [ character(len=error_message_length) :: 'Insufficient general bandwidth in bt' ])

  type(routine_info), parameter :: info_f_c_general_wb=routine_info(id_f_c_general_wb, &
       'f_c_general_wb', &
       [ character(len=error_message_length) :: 'Insufficient general bandwidth.' ])

  ! src/general/general_ubt, 1100
  integer(int32), parameter :: mod_id_general_ubt=1100
  integer(int32), parameter :: id_d_general_to_ubt=mod_id_general_ubt + 0
  integer(int32), parameter :: id_f_d_general_to_ubt=mod_id_general_ubt + 1
  integer(int32), parameter :: id_f_d_general_ubt=mod_id_general_ubt + 2
  integer(int32), parameter :: id_c_general_to_ubt=mod_id_general_ubt + 3
  integer(int32), parameter :: id_f_c_general_to_ubt=mod_id_general_ubt + 4
  integer(int32), parameter :: id_f_c_general_ubt=mod_id_general_ubt + 5

  type(routine_info), parameter :: info_d_general_to_ubt=routine_info(id_d_general_to_ubt, &
       'd_general_to_ubt', &
       [ character(len=error_message_length) :: 'n<1', &
       'Size of a and ubt not the same.' ] )

  type(routine_info), parameter :: info_f_d_general_to_ubt=routine_info(id_f_d_general_to_ubt, &
       'f_d_general_to_ubt', &
       [ character(len=error_message_length) :: 'Insufficient general bandwidth in ubt', &
       'Insufficient general bandwidth in ubt' ])

  type(routine_info), parameter :: info_f_d_general_ubt=routine_info(id_f_d_general_ubt, &
       'f_d_general_ubt', &
       [ character(len=error_message_length) :: '' ])

  type(routine_info), parameter :: info_c_general_to_ubt=routine_info(id_c_general_to_ubt, &
       'c_general_to_ubt', &
       [ character(len=error_message_length) :: 'n<1', &
       'Size of a and ubt not the same.' ] )

  type(routine_info), parameter :: info_f_c_general_to_ubt=routine_info(id_f_c_general_to_ubt, &
       'f_c_general_to_ubt', &
       [ character(len=error_message_length) :: 'Insufficient general bandwidth in ubt', &
       'Insufficient general bandwidth in ubt' ])

  type(routine_info), parameter :: info_f_c_general_ubt=routine_info(id_f_c_general_ubt, &
       'f_c_general_ubt', &
       [ character(len=error_message_length) :: '' ])

  ! src/general/general_wbv, 1200
  integer(int32), parameter :: mod_id_general_wbv=1200
  integer(int32), parameter :: id_d_general_to_wbv=mod_id_general_wbv + 0
  integer(int32), parameter :: id_f_d_general_to_wbv=mod_id_general_wbv + 1
  integer(int32), parameter :: id_f_d_general_wbv=mod_id_general_wbv + 2
  integer(int32), parameter :: id_c_general_to_wbv=mod_id_general_wbv + 3
  integer(int32), parameter :: id_f_c_general_to_wbv=mod_id_general_wbv + 4
  integer(int32), parameter :: id_f_c_general_wbv=mod_id_general_wbv + 5

  type(routine_info), parameter :: info_d_general_to_wbv=routine_info(id_d_general_to_wbv, &
       'd_general_to_wbv', &
       [ character(len=error_message_length) :: 'n<1', &
       'Size of a and wbv not the same.' ] )

  type(routine_info), parameter :: info_f_d_general_to_wbv=routine_info(id_f_d_general_to_wbv, &
       'f_d_general_to_wbv', &
       [ character(len=error_message_length) :: 'Insufficient general bandwidth in wbv', &
       'Insufficient general bandwidth in wbv' ])

  type(routine_info), parameter :: info_f_d_general_wbv=routine_info(id_f_d_general_wbv, &
       'f_d_general_wbv', &
       [ character(len=error_message_length) :: '' ])
  
  type(routine_info), parameter :: info_c_general_to_wbv=routine_info(id_c_general_to_wbv, &
       'c_general_to_wbv', &
       [ character(len=error_message_length) :: 'n<1', &
       'Size of a and wbv not the same.' ] )

  type(routine_info), parameter :: info_f_c_general_to_wbv=routine_info(id_f_c_general_to_wbv, &
       'f_c_general_to_wbv', &
       [ character(len=error_message_length) :: 'Insufficient general bandwidth in wbv', &
       'Insufficient general bandwidth in wbv' ])

  type(routine_info), parameter :: info_f_c_general_wbv=routine_info(id_f_c_general_wbv, &
       'f_c_general_wbv', &
       [ character(len=error_message_length) :: '' ])

  ! src/orth/gs 170s
  integer(int32), parameter :: mod_id_gs=1300
  integer(int32), parameter :: id_d_extend_gs_rows=mod_id_gs + 0
  integer(int32), parameter :: id_c_extend_gs_rows=mod_id_gs + 1
  integer(int32), parameter :: id_d_extend_gs_columns=mod_id_gs + 2
  integer(int32), parameter :: id_c_extend_gs_columns=mod_id_gs + 3

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

  ! src/orth/nullvec 1400
  integer(int32), parameter :: mod_id_nullvec=1400
  integer(int32), parameter :: id_f_d_lower_left_nullvec=mod_id_nullvec + 0
  integer(int32), parameter :: id_f_c_lower_left_nullvec=mod_id_nullvec + 1
  integer(int32), parameter :: id_f_d_lower_right_nullvec=mod_id_nullvec + 2
  integer(int32), parameter :: id_f_c_lower_right_nullvec=mod_id_nullvec + 3

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

  ! src/qr_factorization/qr_factorization 1500
  integer(int32), parameter :: mod_id_qr_factorization=1500
  integer(int32), parameter :: id_d_qr_bv_to_ub=mod_id_qr_factorization + 0
  integer(int32), parameter :: id_f_d_qr_bv_to_ub=mod_id_qr_factorization + 1
  integer(int32), parameter :: id_c_qr_bv_to_ub=mod_id_qr_factorization + 2
  integer(int32), parameter :: id_f_c_qr_bv_to_ub=mod_id_qr_factorization + 3

  type(routine_info), parameter :: info_d_qr_bv_to_ub=routine_info(id_d_qr_bv_to_ub, &
       'd_qr_bv_to_ub', &
       [ character(len=error_message_length) :: 'ub%n /= bv%n or sw%n /= ub%n', &
       'Not enough stroage for sweeps' ])

  type(routine_info), parameter :: info_c_qr_bv_to_ub=routine_info(id_c_qr_bv_to_ub, &
       'c_qr_bv_to_ub', &
       [ character(len=error_message_length) :: 'ub%n /= bv%n', 'bv%lbw <= 0', &
       'dim. of cs or ss /= n' ])

  ! src/solve/back_solve 210s and 220s
  integer(int32), parameter :: mod_id_back_solve=1600
  integer(int32), parameter :: id_d_back_solve_ub=mod_id_back_solve + 0
  integer(int32), parameter :: id_f_d_back_solve_ub=mod_id_back_solve + 1
  integer(int32), parameter :: id_c_back_solve_ub=mod_id_back_solve + 2
  integer(int32), parameter :: id_f_c_back_solve_ub=mod_id_back_solve + 3
  integer(int32), parameter :: id_d_v_back_solve_ub=mod_id_back_solve + 4
  integer(int32), parameter :: id_f_d_v_back_solve_ub=mod_id_back_solve + 5
  integer(int32), parameter :: id_c_v_back_solve_ub=mod_id_back_solve + 6
  integer(int32), parameter :: id_f_c_v_back_solve_ub=mod_id_back_solve + 7

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

  integer(int32), parameter :: id_d_forward_solve_bv=mod_id_back_solve + 8
  integer(int32), parameter :: id_f_d_forward_solve_bv=mod_id_back_solve + 9
  integer(int32), parameter :: id_d_v_forward_solve_bv=mod_id_back_solve + 10
  integer(int32), parameter :: id_f_d_v_forward_solve_bv=mod_id_back_solve + 11
  integer(int32), parameter :: id_c_forward_solve_bv=mod_id_back_solve + 12
  integer(int32), parameter :: id_f_c_forward_solve_bv=mod_id_back_solve + 13
  integer(int32), parameter :: id_c_v_forward_solve_bv=mod_id_back_solve + 14
  integer(int32), parameter :: id_f_c_v_forward_solve_bv=mod_id_back_solve + 15

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

  ! src/row_compress 1700
  integer(int32), parameter :: mod_id_row_compress=1700
  integer(int32), parameter :: id_d_row_compress=mod_id_row_compress + 0
  integer(int32), parameter :: id_f_d_row_compress=mod_id_row_compress + 1
  integer(int32), parameter :: id_c_row_compress=mod_id_row_compress + 2
  integer(int32), parameter :: id_f_c_row_compress=mod_id_row_compress + 3

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

  ! src/types/random 1800
  integer(int32), parameter :: mod_id_random=1800
  integer(int32), parameter :: id_d_random_bc=mod_id_random + 0
  integer(int32), parameter :: id_d_random_br=mod_id_random + 1
  integer(int32), parameter :: id_d_random_ub=mod_id_random + 2
  integer(int32), parameter :: id_d_random_bv=mod_id_random + 3
  integer(int32), parameter :: id_d_random_bt=mod_id_random + 4
  integer(int32), parameter :: id_d_random_wb=mod_id_random + 5
  integer(int32), parameter :: id_d_random_ubt=mod_id_random + 6
  integer(int32), parameter :: id_d_random_wbv=mod_id_random + 7

  integer(int32), parameter :: id_c_random_bc=mod_id_random + 8
  integer(int32), parameter :: id_c_random_br=mod_id_random + 9
  integer(int32), parameter :: id_c_random_ub=mod_id_random + 10
  integer(int32), parameter :: id_c_random_bv=mod_id_random + 11
  integer(int32), parameter :: id_c_random_bt=mod_id_random + 12
  integer(int32), parameter :: id_c_random_wb=mod_id_random + 13
  integer(int32), parameter :: id_c_random_ubt=mod_id_random + 14
  integer(int32), parameter :: id_c_random_wbv=mod_id_random + 15

  type(routine_info), parameter :: info_d_random_bc=routine_info(id_d_random_bc, &
       'd_random_bc', &
       [ character(len=error_message_length) :: 'n < 1', &
       'lbw > lbwmax', 'ubw > ubwmax'])

  type(routine_info), parameter :: info_d_random_br=routine_info(id_d_random_br, &
       'd_random_br', &
       [ character(len=error_message_length) :: 'n < 1', &
       'lbw > lbwmax', 'ubw > ubwmax'])

  type(routine_info), parameter :: info_d_random_ub=routine_info(id_d_random_ub, &
       'd_random_ub', &
       [ character(len=error_message_length) :: 'n < 1', &
       'lbw > lbwmax', 'ubw > ubwmax'])

  type(routine_info), parameter :: info_d_random_bv=routine_info(id_d_random_bv, &
       'd_random_bv', &
       [ character(len=error_message_length) :: 'n < 1', &
       'lbw > lbwmax', 'ubw > ubwmax'])

  type(routine_info), parameter :: info_d_random_bt=routine_info(id_d_random_bt, &
       'd_random_bt', &
       [ character(len=error_message_length) :: 'n < 1', &
       'lbw > lbwmax', 'ubw > ubwmax'])

  type(routine_info), parameter :: info_d_random_wb=routine_info(id_d_random_wb, &
       'd_random_wb', &
       [ character(len=error_message_length) :: 'n < 1', &
       'lbw > lbwmax', 'ubw > ubwmax'])

  type(routine_info), parameter :: info_d_random_ubt=routine_info(id_d_random_ubt, &
       'd_random_ubt', &
       [ character(len=error_message_length) :: 'n < 1', &
       'lbw > lbwmax', 'ubw > ubwmax'])
  
  type(routine_info), parameter :: info_d_random_wbv=routine_info(id_d_random_wbv, &
       'd_random_wbv', &
       [ character(len=error_message_length) :: 'n < 1', &
       'lbw > lbwmax', 'ubw > ubwmax'])

  type(routine_info), parameter :: info_c_random_bc=routine_info(id_c_random_bc, &
       'c_random_bc', &
       [ character(len=error_message_length) :: 'n < 1', &
       'lbw > lbwmax', 'ubw > ubwmax'])

  type(routine_info), parameter :: info_c_random_br=routine_info(id_c_random_br, &
       'c_random_br', &
       [ character(len=error_message_length) :: 'n < 1', &
       'lbw > lbwmax', 'ubw > ubwmax'])

  type(routine_info), parameter :: info_c_random_ub=routine_info(id_c_random_ub, &
       'c_random_ub', &
       [ character(len=error_message_length) :: 'n < 1', &
       'lbw > lbwmax', 'ubw > ubwmax'])

  type(routine_info), parameter :: info_c_random_bv=routine_info(id_c_random_bv, &
       'c_random_bv', &
       [ character(len=error_message_length) :: 'n < 1', &
       'lbw > lbwmax', 'ubw > ubwmax'])

  type(routine_info), parameter :: info_c_random_bt=routine_info(id_c_random_bt, &
       'c_random_bt', &
       [ character(len=error_message_length) :: 'n < 1', &
       'lbw > lbwmax', 'ubw > ubwmax'])

  type(routine_info), parameter :: info_c_random_wb=routine_info(id_c_random_wb, &
       'c_random_wb', &
       [ character(len=error_message_length) :: 'n < 1', &
       'lbw > lbwmax', 'ubw > ubwmax'])

  type(routine_info), parameter :: info_c_random_ubt=routine_info(id_c_random_ubt, &
       'c_random_ubt', &
       [ character(len=error_message_length) :: 'n < 1', &
       'lbw > lbwmax', 'ubw > ubwmax'])
  
  type(routine_info), parameter :: info_c_random_wbv=routine_info(id_c_random_wbv, &
       'c_random_wbv', &
       [ character(len=error_message_length) :: 'n < 1', &
       'lbw > lbwmax', 'ubw > ubwmax'])
  
contains

  subroutine push_id(info,err)
    type(error_info), intent(inout), optional :: err
    type(routine_info), intent(in) :: info    
    if (present(err)) then
       if (err%rix <= max_routines) then
          err%routines(err%rix)=info%routine_id
          err%rix=err%rix+1
       end if
    end if
  end subroutine push_id

  subroutine set_error(code, info, err)
    type(error_info), intent(inout), optional :: err
    type(routine_info), intent(in) :: info
    integer(kind=int32), intent(in) :: code
    if (present(err)) then
       err%code=code
    else
       write(error_unit,*) info%error_messages(code)
       stop
    end if
  end subroutine set_error

  subroutine clear_error(err)
    type(error_info), intent(inout), optional :: err
    if (present(err)) then
       err%code=0
    end if
  end subroutine clear_error

  subroutine clear_routines(err)
    type(error_info), intent(inout), optional :: err
    if (present(err)) then
       err%routines(1:err%rix)=0
       err%rix=1
    end if
  end subroutine clear_routines

  function failure(err) result(f)
    logical :: f
    type(error_info), intent(inout), optional :: err
    f=.false.
    if (present(err)) then
       if (err%code > 0) then
          f=.true.
       end if
    end if
  end function failure

  function success(err) result(s)
    logical :: s
    type(error_info), intent(inout), optional :: err
    s=.true.
    if (present(err)) then
       if (err%code > 0) then
          s=.false.
       end if
    end if
  end function success

  subroutine pop_id(err)
    type(error_info), intent(inout), optional :: err
    if (present(err)) then
       if (err%code <= 0 .and. err%rix > 1) then
          err%rix=err%rix-1
          err%routines(err%rix)=0
       end if
    end if
  end subroutine pop_id

  subroutine initialize_errors
    if (.not. allocated(info_index)) then
       allocate(info_index(max_error_code))
       info_index = info_empty
       ! assemble.f90
       info_index(info_d_ub_to_general%routine_id)=info_d_ub_to_general
       info_index(info_c_ub_to_general%routine_id)=info_c_ub_to_general

       info_index(info_d_bt_to_general%routine_id)=info_d_bt_to_general
       info_index(info_c_bt_to_general%routine_id)=info_c_bt_to_general

       info_index(info_d_ubt_to_general%routine_id)=info_d_ubt_to_general
       info_index(info_c_ubt_to_general%routine_id)=info_c_ubt_to_general

       info_index(info_d_bv_to_general%routine_id)=info_d_bv_to_general
       info_index(info_c_bv_to_general%routine_id)=info_c_bv_to_general

       info_index(info_d_wb_to_general%routine_id)=info_d_wb_to_general
       info_index(info_c_wb_to_general%routine_id)=info_c_wb_to_general

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

       info_index(info_d_general_to_bv%routine_id)=info_d_general_to_bv
       info_index(info_f_d_general_to_bv%routine_id)=info_f_d_general_to_bv
       info_index(info_f_d_general_bv%routine_id)=info_f_d_general_bv
       info_index(info_c_general_to_bv%routine_id)=info_c_general_to_bv
       info_index(info_f_c_general_to_bv%routine_id)=info_f_c_general_to_bv
       info_index(info_f_c_general_bv%routine_id)=info_f_c_general_bv

       ! general_ub

       info_index(info_d_general_to_ub%routine_id)=info_d_general_to_ub
       info_index(info_f_d_general_to_ub%routine_id)=info_f_d_general_to_ub
       info_index(info_f_d_general_ub%routine_id)=info_f_d_general_ub
       info_index(info_c_general_to_ub%routine_id)=info_c_general_to_ub
       info_index(info_f_c_general_to_ub%routine_id)=info_f_c_general_to_ub
       info_index(info_f_c_general_ub%routine_id)=info_f_c_general_ub

       ! general_bt

       info_index(info_d_general_to_bt%routine_id)=info_d_general_to_bt
       info_index(info_f_d_general_to_bt%routine_id)=info_f_d_general_to_bt
       info_index(info_f_d_general_bt%routine_id)=info_f_d_general_bt
       info_index(info_c_general_to_bt%routine_id)=info_c_general_to_bt
       info_index(info_f_c_general_to_bt%routine_id)=info_f_c_general_to_bt
       info_index(info_f_c_general_bt%routine_id)=info_f_c_general_bt

       ! general_wb

       info_index(info_d_general_to_wb%routine_id)=info_d_general_to_wb
       info_index(info_f_d_general_to_wb%routine_id)=info_f_d_general_to_wb
       info_index(info_f_d_general_wb%routine_id)=info_f_d_general_wb
       info_index(info_c_general_to_wb%routine_id)=info_c_general_to_wb
       info_index(info_f_c_general_to_wb%routine_id)=info_f_c_general_to_wb
       info_index(info_f_c_general_wb%routine_id)=info_f_c_general_wb

       ! general_ubt

       info_index(info_d_general_to_ubt%routine_id)=info_d_general_to_ubt
       info_index(info_f_d_general_to_ubt%routine_id)=info_f_d_general_to_ubt
       info_index(info_f_d_general_ubt%routine_id)=info_f_d_general_ubt
       info_index(info_c_general_to_ubt%routine_id)=info_c_general_to_ubt
       info_index(info_f_c_general_to_ubt%routine_id)=info_f_c_general_to_ubt
       info_index(info_f_c_general_ubt%routine_id)=info_f_c_general_ubt

       ! general_wbv

       info_index(info_d_general_to_wbv%routine_id)=info_d_general_to_wbv
       info_index(info_f_d_general_to_wbv%routine_id)=info_f_d_general_to_wbv
       info_index(info_f_d_general_wbv%routine_id)=info_f_d_general_wbv
       info_index(info_c_general_to_wbv%routine_id)=info_c_general_to_wbv
       info_index(info_f_c_general_to_wbv%routine_id)=info_f_c_general_to_wbv
       info_index(info_f_c_general_wbv%routine_id)=info_f_c_general_wbv

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

       ! random

       info_index(info_d_random_bc%routine_id)=info_d_random_bc
       info_index(info_d_random_br%routine_id)=info_d_random_br
       info_index(info_d_random_ub%routine_id)=info_d_random_ub
       info_index(info_d_random_bv%routine_id)=info_d_random_bv
       info_index(info_d_random_bt%routine_id)=info_d_random_bt
       info_index(info_d_random_wb%routine_id)=info_d_random_wb
       info_index(info_d_random_ubt%routine_id)=info_d_random_ubt
       info_index(info_d_random_wbv%routine_id)=info_d_random_wbv

       info_index(info_c_random_bc%routine_id)=info_c_random_bc
       info_index(info_c_random_br%routine_id)=info_c_random_br
       info_index(info_c_random_ub%routine_id)=info_c_random_ub
       info_index(info_c_random_bv%routine_id)=info_c_random_bv
       info_index(info_c_random_bt%routine_id)=info_c_random_bt
       info_index(info_c_random_wb%routine_id)=info_c_random_wb
       info_index(info_c_random_ubt%routine_id)=info_c_random_ubt
       info_index(info_c_random_wbv%routine_id)=info_c_random_wbv
    end if
  end subroutine initialize_errors

end module mod_error_id
