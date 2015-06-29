module mod_error_id
  use mod_prec
  use mod_shift
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

  ! Error type.  error%routines and the index error%rix act as a
  ! stack for keeping track of the calls that led up to the error.
  type error_info
     integer(kind=int32) :: code=0, rix=1
     integer(kind=int32), dimension(max_routines) :: routines=0
     logical :: halt=.true. ! If .true., halt at the first error.
  end type error_info

  type(routine_info), parameter :: info_empty=routine_info(0, &
       '', [ character(len=error_message_length) :: '' ] )

  type(routine_info), dimension(:), allocatable :: info_index

  ! src/assemble/assemble.f90, 0000
  integer(int32), parameter :: mod_id_assemble=0
  integer(int32), parameter :: id_d_ub_to_general=mod_id_assemble + 1
  integer(int32), parameter :: id_f_d_ub_to_general=mod_id_assemble + 0002
  integer(int32), parameter :: id_d_general_of_ub=mod_id_assemble + 0003
  integer(int32), parameter :: id_z_ub_to_general=mod_id_assemble + 0004
  integer(int32), parameter :: id_f_z_ub_to_general=mod_id_assemble + 0005
  integer(int32), parameter :: id_z_general_of_ub=mod_id_assemble + 0006

  integer(int32), parameter :: id_d_bt_to_general=mod_id_assemble + 0007
  integer(int32), parameter :: id_f_d_bt_to_general=mod_id_assemble + 0008
  integer(int32), parameter :: id_d_general_of_bt=mod_id_assemble + 0009
  integer(int32), parameter :: id_z_bt_to_general=mod_id_assemble + 0010
  integer(int32), parameter :: id_f_z_bt_to_general=mod_id_assemble + 0011
  integer(int32), parameter :: id_z_general_of_bt=mod_id_assemble + 0012
  
  integer(int32), parameter :: id_d_ubt_to_general=mod_id_assemble + 0013
  integer(int32), parameter :: id_f_d_ubt_to_general=mod_id_assemble + 0014
  integer(int32), parameter :: id_d_general_of_ubt=mod_id_assemble + 0015
  integer(int32), parameter :: id_z_ubt_to_general=mod_id_assemble + 0016
  integer(int32), parameter :: id_f_z_ubt_to_general=mod_id_assemble + 0017
  integer(int32), parameter :: id_z_general_of_ubt=mod_id_assemble + 0018

  integer(int32), parameter :: id_d_bv_to_general=mod_id_assemble + 0019
  integer(int32), parameter :: id_f_d_bv_to_general=mod_id_assemble + 0020
  integer(int32), parameter :: id_d_general_of_bv=mod_id_assemble + 0021
  integer(int32), parameter :: id_z_bv_to_general=mod_id_assemble + 0022
  integer(int32), parameter :: id_f_z_bv_to_general=mod_id_assemble + 0023
  integer(int32), parameter :: id_z_general_of_bv=mod_id_assemble + 0024

  integer(int32), parameter :: id_d_wb_to_general=mod_id_assemble + 0025
  integer(int32), parameter :: id_f_d_wb_to_general=mod_id_assemble + 0026
  integer(int32), parameter :: id_d_general_of_wb=mod_id_assemble + 0027
  integer(int32), parameter :: id_z_wb_to_general=mod_id_assemble + 0028
  integer(int32), parameter :: id_f_z_wb_to_general=mod_id_assemble + 0029
  integer(int32), parameter :: id_z_general_of_wb=mod_id_assemble + 0030
  
  integer(int32), parameter :: id_d_wbv_to_general=mod_id_assemble + 0031
  integer(int32), parameter :: id_f_d_wbv_to_general=mod_id_assemble + 0032
  integer(int32), parameter :: id_d_general_of_wbv=mod_id_assemble + 0033
  integer(int32), parameter :: id_z_wbv_to_general=mod_id_assemble + 0034
  integer(int32), parameter :: id_f_z_wbv_to_general=mod_id_assemble + 0035
  integer(int32), parameter :: id_z_general_of_wbv=mod_id_assemble + 0036

  type(routine_info), parameter :: info_d_ub_to_general=routine_info(id_d_ub_to_general, &
       'd_ub_to_general', [ character(len=error_message_length) :: 'Size error in A.' ] )

  type(routine_info), parameter :: info_z_ub_to_general=routine_info(id_z_ub_to_general, &
       'z_ub_to_general', [ character(len=error_message_length) :: 'Size error in A.' ] )

  type(routine_info), parameter :: info_d_general_of_ub=routine_info(id_d_general_of_ub, &
       'd_general_of_ub', [ character(len=error_message_length) :: '' ] )

  type(routine_info), parameter :: info_z_general_of_ub=routine_info(id_z_general_of_ub, &
       'z_general_of_ub', [ character(len=error_message_length) :: '' ] )
  
  type(routine_info), parameter :: info_d_bt_to_general=routine_info(id_d_bt_to_general, &
       'd_bt_to_general', [ character(len=error_message_length) :: 'Size error in A.' ] )

  type(routine_info), parameter :: info_z_bt_to_general=routine_info(id_z_bt_to_general, &
       'z_bt_to_general', [ character(len=error_message_length) :: 'Size error in A.' ] )

  type(routine_info), parameter :: info_d_general_of_bt=routine_info(id_d_general_of_bt, &
       'd_general_of_bt', [ character(len=error_message_length) :: '' ] )

  type(routine_info), parameter :: info_z_general_of_bt=routine_info(id_z_general_of_bt, &
       'z_general_of_bt', [ character(len=error_message_length) :: '' ] )

  type(routine_info), parameter :: info_d_ubt_to_general=routine_info(id_d_ubt_to_general, &
       'd_ubt_to_general', [ character(len=error_message_length) :: 'Size error in A.' ] )

  type(routine_info), parameter :: info_z_ubt_to_general=routine_info(id_z_ubt_to_general, &
       'z_ubt_to_general', [ character(len=error_message_length) :: 'Size error in A.' ] )

  type(routine_info), parameter :: info_d_general_of_ubt=routine_info(id_d_general_of_ubt, &
       'd_general_of_ubt', [ character(len=error_message_length) :: '' ] )

  type(routine_info), parameter :: info_z_general_of_ubt=routine_info(id_z_general_of_ubt, &
       'z_general_of_ubt', [ character(len=error_message_length) :: '' ] )

  type(routine_info), parameter :: info_d_bv_to_general=routine_info(id_d_bv_to_general, &
       'd_bv_to_general', [ character(len=error_message_length) :: 'Size error in A.' ] )
  
  type(routine_info), parameter :: info_z_bv_to_general=routine_info(id_z_bv_to_general, &
       'z_bv_to_general', [ character(len=error_message_length) :: 'Size error in A.' ] )

  type(routine_info), parameter :: info_d_general_of_bv=routine_info(id_d_general_of_bv, &
       'd_general_of_bv', [ character(len=error_message_length) :: '' ] )

  type(routine_info), parameter :: info_z_general_of_bv=routine_info(id_z_general_of_bv, &
       'z_general_of_bv', [ character(len=error_message_length) :: '' ] )

  type(routine_info), parameter :: info_d_wb_to_general=routine_info(id_d_wb_to_general, &
       'd_wb_to_general', [ character(len=error_message_length) :: 'Size error in A.' ] )
  
  type(routine_info), parameter :: info_z_wb_to_general=routine_info(id_z_wb_to_general, &
       'z_wb_to_general', [ character(len=error_message_length) :: 'Size error in A.' ] )

  type(routine_info), parameter :: info_d_general_of_wb=routine_info(id_d_general_of_wb, &
       'd_general_of_wb', [ character(len=error_message_length) :: '' ] )

  type(routine_info), parameter :: info_z_general_of_wb=routine_info(id_z_general_of_wb, &
       'z_general_of_wb', [ character(len=error_message_length) :: '' ] )
  
  type(routine_info), parameter :: info_d_wbv_to_general=routine_info(id_d_wbv_to_general, &
       'd_wbv_to_general', [ character(len=error_message_length) :: 'Size error in A.' ] )
  
  type(routine_info), parameter :: info_z_wbv_to_general=routine_info(id_z_wbv_to_general, &
       'z_wbv_to_general', [ character(len=error_message_length) :: 'Size error in A.' ] )

  type(routine_info), parameter :: info_d_general_of_wbv=routine_info(id_d_general_of_wbv, &
       'd_general_of_wbv', [ character(len=error_message_length) :: '' ] )

  type(routine_info), parameter :: info_z_general_of_wbv=routine_info(id_z_general_of_wbv, &
       'z_general_of_wbv', [ character(len=error_message_length) :: '' ] )

  ! src/convert/convert_bv_to_ub, 0100
  integer(int32), parameter :: mod_id_convert_bv_to_ub=100
  integer(int32), parameter :: id_d_convert_bv_to_ub=mod_id_convert_bv_to_ub+0
  integer(int32), parameter :: id_f_d_convert_bv_to_ub=mod_id_convert_bv_to_ub+1
  integer(int32), parameter :: id_d_ub_of_bv=mod_id_convert_bv_to_ub+2
  integer(int32), parameter :: id_z_convert_bv_to_ub=mod_id_convert_bv_to_ub+3
  integer(int32), parameter :: id_f_z_convert_bv_to_ub=mod_id_convert_bv_to_ub+4
  integer(int32), parameter :: id_z_ub_of_bv=mod_id_convert_bv_to_ub+5

  type(routine_info), parameter :: info_d_convert_bv_to_ub=routine_info(id_d_convert_bv_to_ub, &
       'd_convert_bv_to_ub', &
       [ character(len=error_message_length) :: 'n<1', 'Insufficient storage in bv.', &
       'Insufficient Storage in ub.', 'ub%n /= bv%n' ] )

  type(routine_info), parameter :: info_d_ub_of_bv=routine_info(id_d_ub_of_bv, &
       'd_ub_of_bv', &
       [ character(len=error_message_length) :: '' ] )

  type(routine_info), parameter :: info_z_convert_bv_to_ub=routine_info(id_z_convert_bv_to_ub, &
       'z_convert_bv_to_ub', &
       [ character(len=error_message_length) :: 'n<1', 'Insufficient storage in bv.', &
       'Insufficient Storage in ub.', 'ub%n /= bv%n' ] )

  type(routine_info), parameter :: info_z_ub_of_bv=routine_info(id_z_ub_of_bv, &
       'z_ub_of_bv', &
       [ character(len=error_message_length) :: '' ] )

  ! src/convert/convert_ub_to_bv, 200
  integer(int32), parameter :: mod_id_convert_ub_to_bv=200
  integer(int32), parameter :: id_d_convert_ub_to_bv=mod_id_convert_ub_to_bv + 0
  integer(int32), parameter :: id_f_d_convert_ub_to_bv=mod_id_convert_ub_to_bv + 1
  integer(int32), parameter :: id_d_bv_of_ub=mod_id_convert_ub_to_bv + 2
  integer(int32), parameter :: id_z_convert_ub_to_bv=mod_id_convert_ub_to_bv + 3
  integer(int32), parameter :: id_f_z_convert_ub_to_bv=mod_id_convert_ub_to_bv + 4
  integer(int32), parameter :: id_z_bv_of_ub=mod_id_convert_ub_to_bv + 5

  type(routine_info), parameter :: info_d_convert_ub_to_bv=routine_info(id_d_convert_ub_to_bv, &
       'd_convert_ub_to_bv', &
       [ character(len=error_message_length) :: 'n<1', 'Insufficient storage in ub', &
       'Insufficient storage in bv.', 'ub%n /= bv%n' ] )

  type(routine_info), parameter :: info_d_bv_of_ub=routine_info(id_d_bv_of_ub, &
       'd_bv_of_ub', [ character(len=error_message_length) :: '' ] )
  
  type(routine_info), parameter :: info_z_convert_ub_to_bv=routine_info(id_z_convert_ub_to_bv, &
       'z_convert_ub_to_bv', &
       [ character(len=error_message_length) :: 'Insufficient storage in ub', &
       'Insufficient storage in bv.', 'ub%n /= bv%n' ] )

  type(routine_info), parameter :: info_z_bv_of_ub=routine_info(id_z_bv_of_ub, &
       'z_bv_of_ub', [ character(len=error_message_length) :: '' ] )

  ! src/convert/convert_wb_to_bt, 300
  integer(int32), parameter :: mod_id_convert_wb_to_bt=300
  integer(int32), parameter :: id_d_convert_wb_to_bt=mod_id_convert_wb_to_bt + 0
  integer(int32), parameter :: id_f_d_convert_wb_to_bt=mod_id_convert_wb_to_bt + 1
  integer(int32), parameter :: id_d_bt_of_wb=mod_id_convert_wb_to_bt + 2
  integer(int32), parameter :: id_z_convert_wb_to_bt=mod_id_convert_wb_to_bt + 3
  integer(int32), parameter :: id_f_z_convert_wb_to_bt=mod_id_convert_wb_to_bt + 4
  integer(int32), parameter :: id_z_bt_of_wb=mod_id_convert_wb_to_bt + 5

  type(routine_info), parameter :: info_d_convert_wb_to_bt=routine_info(id_d_convert_wb_to_bt, &
       'd_convert_wb_to_bt', &
       [ character(len=error_message_length) :: 'n<1', 'Insufficient storage in wb.', &
       'Insufficient Storage in bt.', 'bt%n /= wb%n' ] )

  type(routine_info), parameter :: info_d_bt_of_wb=routine_info(id_d_bt_of_wb, &
       'd_bt_of_wb', &
       [ character(len=error_message_length) :: '' ] )

  type(routine_info), parameter :: info_z_convert_wb_to_bt=routine_info(id_z_convert_wb_to_bt, &
       'z_convert_wb_to_bt', &
       [ character(len=error_message_length) :: 'n<1', 'Insufficient storage in wb.', &
       'Insufficient Storage in bt.', 'bt%n /= wb%n' ] )

  type(routine_info), parameter :: info_z_bt_of_wb=routine_info(id_z_bt_of_wb, &
       'z_bt_of_wb', &
       [ character(len=error_message_length) :: '' ] )

  ! src/convert/convert_bt_to_wb, 400
  integer(int32), parameter :: mod_id_convert_bt_to_wb=400
  integer(int32), parameter :: id_d_convert_bt_to_wb=mod_id_convert_bt_to_wb + 0
  integer(int32), parameter :: id_f_d_convert_bt_to_wb=mod_id_convert_bt_to_wb + 1
  integer(int32), parameter :: id_d_wb_of_bt=mod_id_convert_bt_to_wb + 2
  integer(int32), parameter :: id_z_convert_bt_to_wb=mod_id_convert_bt_to_wb + 3
  integer(int32), parameter :: id_f_z_convert_bt_to_wb=mod_id_convert_bt_to_wb + 4
  integer(int32), parameter :: id_z_wb_of_bt=mod_id_convert_bt_to_wb + 5
  
  type(routine_info), parameter :: info_d_convert_bt_to_wb=routine_info(id_d_convert_bt_to_wb, &
       'd_convert_bt_to_wb', &
       [ character(len=error_message_length) :: 'n<1', 'Insufficient storage in bt.', &
       'Insufficient Storage in wb.', 'wb%n /= bt%n' ] )

  type(routine_info), parameter :: info_d_wb_of_bt=routine_info(id_d_wb_of_bt, &
       'd_wb_of_bt', &
       [ character(len=error_message_length) :: '' ] )
  
  type(routine_info), parameter :: info_z_convert_bt_to_wb=routine_info(id_z_convert_bt_to_wb, &
       'z_convert_bt_to_wb', &
       [ character(len=error_message_length) :: 'n<1', 'Insufficient storage in bt.', &
       'Insufficient Storage in wb.', 'wb%n /= bt%n' ] )

  type(routine_info), parameter :: info_z_wb_of_bt=routine_info(id_z_wb_of_bt, &
       'z_wb_of_bt', &
       [ character(len=error_message_length) :: '' ] )

  ! src/convert/convert_wbv_to_ubt, 500
  integer(int32), parameter :: mod_id_convert_wbv_to_ubt=500
  integer(int32), parameter :: id_d_convert_wbv_to_ubt=mod_id_convert_wbv_to_ubt + 0
  integer(int32), parameter :: id_f_d_convert_wbv_to_ubt=mod_id_convert_wbv_to_ubt + 1
  integer(int32), parameter :: id_d_ubt_of_wbv=mod_id_convert_wbv_to_ubt + 2
  integer(int32), parameter :: id_z_convert_wbv_to_ubt=mod_id_convert_wbv_to_ubt + 3
  integer(int32), parameter :: id_f_z_convert_wbv_to_ubt=mod_id_convert_wbv_to_ubt + 4
  integer(int32), parameter :: id_z_ubt_of_wbv=mod_id_convert_wbv_to_ubt + 5

  type(routine_info), parameter :: info_d_convert_wbv_to_ubt=routine_info(id_d_convert_wbv_to_ubt, &
       'd_convert_wbv_to_ubt', &
       [ character(len=error_message_length) :: 'n<1', 'Insufficient storage in wbv', &
       'Insufficient storage in ubt.', 'wbv%n /= ubt%n' ] )

  type(routine_info), parameter :: info_d_ubt_of_wbv=routine_info(id_d_ubt_of_wbv, &
       'd_ubt_of_wbv', &
       [ character(len=error_message_length) :: '' ] )
  
  type(routine_info), parameter :: info_z_convert_wbv_to_ubt=routine_info(id_z_convert_wbv_to_ubt, &
       'z_convert_wbv_to_ubt', &
       [ character(len=error_message_length) :: 'Insufficient storage in wbv', &
       'Insufficient storage in ubt.', 'wbv%n /= ubt%n' ] )

  type(routine_info), parameter :: info_z_ubt_of_wbv=routine_info(id_z_ubt_of_wbv, &
       'z_ubt_of_wbv', &
       [ character(len=error_message_length) :: '' ] )

  ! src/convert/convert_ubt_to_wbv, 600
  integer(int32), parameter :: mod_id_convert_ubt_to_wbv=600
  integer(int32), parameter :: id_d_convert_ubt_to_wbv=mod_id_convert_ubt_to_wbv + 0
  integer(int32), parameter :: id_f_d_convert_ubt_to_wbv=mod_id_convert_ubt_to_wbv + 1
  integer(int32), parameter :: id_d_wbv_of_ubt=mod_id_convert_ubt_to_wbv + 2
  integer(int32), parameter :: id_z_convert_ubt_to_wbv=mod_id_convert_ubt_to_wbv + 3
  integer(int32), parameter :: id_f_z_convert_ubt_to_wbv=mod_id_convert_ubt_to_wbv + 4
  integer(int32), parameter :: id_z_wbv_of_ubt=mod_id_convert_ubt_to_wbv + 5
  
  type(routine_info), parameter :: info_d_convert_ubt_to_wbv=routine_info(id_d_convert_ubt_to_wbv, &
       'd_convert_ubt_to_wbv', &
       [ character(len=error_message_length) :: 'n<1', 'Insufficient storage in ubt', &
       'Insufficient storage in wbv.', 'ubt%n /= wbv%n' ] )

  type(routine_info), parameter :: info_d_wbv_of_ubt=routine_info(id_d_wbv_of_ubt, &
       'd_wbv_of_ubt', &
       [ character(len=error_message_length) :: '' ] )
  
  type(routine_info), parameter :: info_z_convert_ubt_to_wbv=routine_info(id_z_convert_ubt_to_wbv, &
       'z_convert_ubt_to_wbv', &
       [ character(len=error_message_length) :: 'Insufficient storage in ubt', &
       'Insufficient storage in wbv.', 'ubt%n /= wbv%n' ] )

  type(routine_info), parameter :: info_z_wbv_of_ubt=routine_info(id_z_wbv_of_ubt, &
       'z_wbv_of_ubt', &
       [ character(len=error_message_length) :: '' ] )

  ! src/general/general_bv, 700
  integer(int32), parameter :: mod_id_general_bv=700
  integer(int32), parameter :: id_d_bv_of_general=mod_id_general_bv + 0
  integer(int32), parameter :: id_d_general_to_bv=mod_id_general_bv + 1
  integer(int32), parameter :: id_f_d_general_to_bv=mod_id_general_bv + 2
  integer(int32), parameter :: id_f_d_general_bv=mod_id_general_bv + 3
  integer(int32), parameter :: id_z_bv_of_general=mod_id_general_bv + 4
  integer(int32), parameter :: id_z_general_to_bv=mod_id_general_bv + 5
  integer(int32), parameter :: id_f_z_general_to_bv=mod_id_general_bv + 6
  integer(int32), parameter :: id_f_z_general_bv=mod_id_general_bv + 7

  type(routine_info), parameter :: info_d_bv_of_general=routine_info(id_d_bv_of_general, &
       'd_bv_of_general', &
       [ character(len=error_message_length) :: '' ] )

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

  type(routine_info), parameter :: info_z_bv_of_general=routine_info(id_z_bv_of_general, &
       'z_bv_of_general', &
       [ character(len=error_message_length) :: '' ] )
  
  type(routine_info), parameter :: info_z_general_to_bv=routine_info(id_z_general_to_bv, &
       'z_general_to_bv', &
       [ character(len=error_message_length) :: 'n<1', 'bv%lbwmax < lbw', &
       'Size of a and bv not the same.' ] )

  type(routine_info), parameter :: info_f_z_general_to_bv=routine_info(id_f_z_general_to_bv, &
       'f_z_general_to_bv', &
       [ character(len=error_message_length) :: 'n<1', 'Insufficient General Bandwidth in bv' ])

  type(routine_info), parameter :: info_f_z_general_bv=routine_info(id_f_z_general_bv, &
       'f_z_general_bv', &
       [ character(len=error_message_length) :: 'Insufficient General Bandwidth.' ])
  
  !
  ! src/general/general_ub, 800
  !
  integer(int32), parameter :: mod_id_general_ub=800
  integer(int32), parameter :: id_d_ub_of_general=mod_id_general_ub + 0
  integer(int32), parameter :: id_d_general_to_ub=mod_id_general_ub + 1
  integer(int32), parameter :: id_f_d_general_to_ub=mod_id_general_ub + 2
  integer(int32), parameter :: id_f_d_general_ub=mod_id_general_ub + 3
  integer(int32), parameter :: id_z_ub_of_general=mod_id_general_ub + 4
  integer(int32), parameter :: id_z_general_to_ub=mod_id_general_ub + 5
  integer(int32), parameter :: id_f_z_general_to_ub=mod_id_general_ub + 6
  integer(int32), parameter :: id_f_z_general_ub=mod_id_general_ub + 7

  type(routine_info), parameter :: info_d_ub_of_general=routine_info(id_d_ub_of_general, &
       'd_ub_of_general', &
       [ character(len=error_message_length) :: '' ] )

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

  type(routine_info), parameter :: info_z_ub_of_general=routine_info(id_z_ub_of_general, &
       'z_ub_of_general', &
       [ character(len=error_message_length) :: '' ] )

  type(routine_info), parameter :: info_z_general_to_ub=routine_info(id_z_general_to_ub, &
       'z_general_to_ub', &
       [ character(len=error_message_length) :: 'n<1', 'ub%lbwmax < lbw', &
       'Size of a and ub not the same.' ] )

  type(routine_info), parameter :: info_f_z_general_to_ub=routine_info(id_f_z_general_to_ub, &
       'f_z_general_to_ub', &
       [ character(len=error_message_length) :: 'Insufficient General Bandwidth in ub' ])

  type(routine_info), parameter :: info_f_z_general_ub=routine_info(id_f_z_general_ub, &
       'f_z_general_ub', &
       [ character(len=error_message_length) :: 'Insufficient General Bandwidth.' ])

  ! /src/general/general_bt 900
  integer(int32), parameter :: mod_id_general_bt=900
  integer(int32), parameter :: id_d_bt_of_general=mod_id_general_bt + 0
  integer(int32), parameter :: id_d_general_to_bt=mod_id_general_bt + 1
  integer(int32), parameter :: id_f_d_general_to_bt=mod_id_general_bt + 2
  integer(int32), parameter :: id_f_d_general_bt=mod_id_general_bt + 3
  integer(int32), parameter :: id_z_bt_of_general=mod_id_general_bt + 4
  integer(int32), parameter :: id_z_general_to_bt=mod_id_general_bt + 5
  integer(int32), parameter :: id_f_z_general_to_bt=mod_id_general_bt + 6
  integer(int32), parameter :: id_f_z_general_bt=mod_id_general_bt + 7

  type(routine_info), parameter :: info_d_bt_of_general=routine_info(id_d_bt_of_general, &
       'd_bt_of_general', &
       [ character(len=error_message_length) :: '' ] )

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

  type(routine_info), parameter :: info_z_bt_of_general=routine_info(id_z_bt_of_general, &
       'z_bt_of_general', &
       [ character(len=error_message_length) :: '' ] )
  
  type(routine_info), parameter :: info_z_general_to_bt=routine_info(id_z_general_to_bt, &
       'z_general_to_bt', &
       [ character(len=error_message_length) :: 'n<1', &
       'bt%ubwmax < ubw', 'Size of a and bt not the same.' ] )

  type(routine_info), parameter :: info_f_z_general_to_bt=routine_info(id_f_z_general_to_bt, &
       'f_z_general_to_bt', &
       [ character(len=error_message_length) :: 'Insufficient general bandwidth in bt' ])

  type(routine_info), parameter :: info_f_z_general_bt=routine_info(id_f_z_general_bt, &
       'f_z_general_bt', &
       [ character(len=error_message_length) :: 'Insufficient general bandwidth.' ])


  ! /src/general/general_wb 1000
  integer(int32), parameter :: mod_id_general_wb=1000
  integer(int32), parameter :: id_d_wb_of_general=mod_id_general_wb + 0
  integer(int32), parameter :: id_d_general_to_wb=mod_id_general_wb + 1
  integer(int32), parameter :: id_f_d_general_to_wb=mod_id_general_wb + 2
  integer(int32), parameter :: id_f_d_general_wb=mod_id_general_wb + 3
  integer(int32), parameter :: id_z_wb_of_general=mod_id_general_wb + 4
  integer(int32), parameter :: id_z_general_to_wb=mod_id_general_wb + 5
  integer(int32), parameter :: id_f_z_general_to_wb=mod_id_general_wb + 6
  integer(int32), parameter :: id_f_z_general_wb=mod_id_general_wb + 7
  
  type(routine_info), parameter :: info_d_wb_of_general=routine_info(id_d_wb_of_general, &
       'd_wb_of_general', &
       [ character(len=error_message_length) :: '' ] )

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

  type(routine_info), parameter :: info_z_wb_of_general=routine_info(id_z_wb_of_general, &
       'z_wb_of_general', &
       [ character(len=error_message_length) :: '' ] )

  type(routine_info), parameter :: info_z_general_to_wb=routine_info(id_z_general_to_wb, &
       'z_general_to_wb', &
       [ character(len=error_message_length) :: 'n<1', 'bt%ubwmax < ubw', &
       'Size of a and bt not the same.' ] )

  type(routine_info), parameter :: info_f_z_general_to_wb=routine_info(id_f_z_general_to_wb, &
       'f_z_general_to_wb', &
       [ character(len=error_message_length) :: 'Insufficient general bandwidth in bt' ])

  type(routine_info), parameter :: info_f_z_general_wb=routine_info(id_f_z_general_wb, &
       'f_z_general_wb', &
       [ character(len=error_message_length) :: 'Insufficient general bandwidth.' ])

  ! src/general/general_ubt, 1100
  integer(int32), parameter :: mod_id_general_ubt=1100
  integer(int32), parameter :: id_d_ubt_of_general=mod_id_general_ubt + 0
  integer(int32), parameter :: id_d_general_to_ubt=mod_id_general_ubt + 1
  integer(int32), parameter :: id_f_d_general_to_ubt=mod_id_general_ubt + 2
  integer(int32), parameter :: id_f_d_general_ubt=mod_id_general_ubt + 3
  integer(int32), parameter :: id_z_ubt_of_general=mod_id_general_ubt + 4
  integer(int32), parameter :: id_z_general_to_ubt=mod_id_general_ubt + 5
  integer(int32), parameter :: id_f_z_general_to_ubt=mod_id_general_ubt + 6
  integer(int32), parameter :: id_f_z_general_ubt=mod_id_general_ubt + 7

  type(routine_info), parameter :: info_d_ubt_of_general=routine_info(id_d_ubt_of_general, &
       'd_ubt_of_general', &
       [ character(len=error_message_length) :: '' ] )
  
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

  type(routine_info), parameter :: info_z_ubt_of_general=routine_info(id_z_ubt_of_general, &
       'z_ubt_of_general', &
       [ character(len=error_message_length) :: '' ] )

  type(routine_info), parameter :: info_z_general_to_ubt=routine_info(id_z_general_to_ubt, &
       'z_general_to_ubt', &
       [ character(len=error_message_length) :: 'n<1', &
       'Size of a and ubt not the same.' ] )

  type(routine_info), parameter :: info_f_z_general_to_ubt=routine_info(id_f_z_general_to_ubt, &
       'f_z_general_to_ubt', &
       [ character(len=error_message_length) :: 'Insufficient general bandwidth in ubt', &
       'Insufficient general bandwidth in ubt' ])

  type(routine_info), parameter :: info_f_z_general_ubt=routine_info(id_f_z_general_ubt, &
       'f_z_general_ubt', &
       [ character(len=error_message_length) :: '' ])

  ! src/general/general_wbv, 1200
  integer(int32), parameter :: mod_id_general_wbv=1200
  integer(int32), parameter :: id_d_wbv_of_general=mod_id_general_wbv + 0
  integer(int32), parameter :: id_d_general_to_wbv=mod_id_general_wbv + 1
  integer(int32), parameter :: id_f_d_general_to_wbv=mod_id_general_wbv + 2
  integer(int32), parameter :: id_f_d_general_wbv=mod_id_general_wbv + 3
  integer(int32), parameter :: id_z_wbv_of_general=mod_id_general_wbv + 4
  integer(int32), parameter :: id_z_general_to_wbv=mod_id_general_wbv + 5
  integer(int32), parameter :: id_f_z_general_to_wbv=mod_id_general_wbv + 6
  integer(int32), parameter :: id_f_z_general_wbv=mod_id_general_wbv + 7

  type(routine_info), parameter :: info_d_wbv_of_general=routine_info(id_d_wbv_of_general, &
       'd_wbv_of_general', &
       [ character(len=error_message_length) :: '' ] )

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
  
  type(routine_info), parameter :: info_z_wbv_of_general=routine_info(id_z_wbv_of_general, &
       'z_wbv_of_general', &
       [ character(len=error_message_length) :: '' ] )

  type(routine_info), parameter :: info_z_general_to_wbv=routine_info(id_z_general_to_wbv, &
       'z_general_to_wbv', &
       [ character(len=error_message_length) :: 'n<1', &
       'Size of a and wbv not the same.' ] )

  type(routine_info), parameter :: info_f_z_general_to_wbv=routine_info(id_f_z_general_to_wbv, &
       'f_z_general_to_wbv', &
       [ character(len=error_message_length) :: 'Insufficient general bandwidth in wbv', &
       'Insufficient general bandwidth in wbv' ])

  type(routine_info), parameter :: info_f_z_general_wbv=routine_info(id_f_z_general_wbv, &
       'f_z_general_wbv', &
       [ character(len=error_message_length) :: '' ])

  ! src/orth/gs 1300
  integer(int32), parameter :: mod_id_gs=1300
  integer(int32), parameter :: id_d_extend_gs_rows=mod_id_gs + 0
  integer(int32), parameter :: id_z_extend_gs_rows=mod_id_gs + 1
  integer(int32), parameter :: id_d_extend_gs_columns=mod_id_gs + 2
  integer(int32), parameter :: id_z_extend_gs_columns=mod_id_gs + 3

  type(routine_info), parameter :: info_d_extend_gs_rows=routine_info(id_d_extend_gs_rows, &
       'd_extend_gs_rows', &
       [ character(len=error_message_length) :: 'GS Orthogonalization Error' ])

  type(routine_info), parameter :: info_z_extend_gs_rows=routine_info(id_z_extend_gs_rows, &
       'z_extend_gs_rows', &
       [ character(len=error_message_length) :: 'GS Orthogonalization Error' ])

  type(routine_info), parameter :: info_d_extend_gs_columns=routine_info(id_d_extend_gs_columns, &
       'd_extend_gs_columns', &
       [ character(len=error_message_length) :: 'GS Orthogonalization Error' ])

  type(routine_info), parameter :: info_z_extend_gs_columns=routine_info(id_z_extend_gs_columns, &
       'z_extend_gs_columns', &
       [ character(len=error_message_length) :: 'GS Orthogonalization Error' ])

  ! src/orth/cond_triangular 1400
  integer(int32), parameter :: mod_id_cond_triangular=1400
  integer(int32), parameter :: id_d_lower_left_nullvec=mod_id_cond_triangular + 0
  integer(int32), parameter :: id_z_lower_left_nullvec=mod_id_cond_triangular + 1
  integer(int32), parameter :: id_d_lower_right_nullvec=mod_id_cond_triangular + 2
  integer(int32), parameter :: id_z_lower_right_nullvec=mod_id_cond_triangular + 3
  integer(int32), parameter :: id_d_lower_min_sv=mod_id_cond_triangular + 4
  integer(int32), parameter :: id_z_lower_min_sv=mod_id_cond_triangular + 5
  integer(int32), parameter :: id_d_upper_min_sv=mod_id_cond_triangular + 6
  integer(int32), parameter :: id_z_upper_min_sv=mod_id_cond_triangular + 7
  integer(int32), parameter :: id_d_lower_max_sv=mod_id_cond_triangular + 8
  integer(int32), parameter :: id_z_lower_max_sv=mod_id_cond_triangular + 9
  integer(int32), parameter :: id_d_upper_max_sv=mod_id_cond_triangular + 10
  integer(int32), parameter :: id_z_upper_max_sv=mod_id_cond_triangular + 11

  type(routine_info), parameter :: info_d_lower_left_nullvec= &
       routine_info(id_d_lower_left_nullvec, 'd_lower_left_nullvec', &
       [ character(len=error_message_length) :: 'Failure to find a null vector.' ])
  type(routine_info), parameter :: info_d_lower_right_nullvec= &
       routine_info(id_d_lower_right_nullvec, 'd_lower_right_nullvec', &
       [ character(len=error_message_length) :: 'Failure to find a null vector.' ])
  type(routine_info), parameter :: info_z_lower_left_nullvec= &
       routine_info(id_z_lower_left_nullvec, 'z_lower_left_nullvec', &
       [ character(len=error_message_length) :: 'Failure to find a null vector.' ])
  type(routine_info), parameter :: info_z_lower_right_nullvec= &
       routine_info(id_z_lower_right_nullvec, 'z_lower_right_nullvec', &
       [ character(len=error_message_length) :: 'Failure to find a null vector.' ])

  type(routine_info), parameter :: info_d_lower_min_sv= &
       routine_info(id_d_lower_min_sv, 'd_lower_min_sv', &
       [ character(len=error_message_length) :: 'n<1', 'A is not square.', &
       'Singular vector size error', 'Failure to converge.'])
  type(routine_info), parameter :: info_z_lower_min_sv= &
       routine_info(id_z_lower_min_sv, 'z_lower_min_sv', &
       [ character(len=error_message_length) :: 'n<1', 'A is not square.', &
       'Singular vector size error', 'Failure to converge.'])
  type(routine_info), parameter :: info_d_upper_min_sv= &
       routine_info(id_d_upper_min_sv, 'd_upper_min_sv', &
       [ character(len=error_message_length) :: 'n<1', 'A is not square.', &
       'Singular vector size error', 'Failure to converge.'])
  type(routine_info), parameter :: info_z_upper_min_sv= &
       routine_info(id_z_upper_min_sv, 'z_upper_min_sv', &
       [ character(len=error_message_length) :: 'n<1', 'A is not square.', &
       'Singular vector size error', 'Failure to converge.'])

  type(routine_info), parameter :: info_d_lower_max_sv= &
       routine_info(id_d_lower_max_sv, 'd_lower_max_sv', &
       [ character(len=error_message_length) :: 'n<1', 'A is not square.', &
       'Singular vector size error', 'Failure to converge.'])
  type(routine_info), parameter :: info_z_lower_max_sv= &
       routine_info(id_z_lower_max_sv, 'z_lower_max_sv', &
       [ character(len=error_message_length) :: 'n<1', 'A is not square.', &
       'Singular vector size error', 'Failure to converge.'])
  type(routine_info), parameter :: info_d_upper_max_sv= &
       routine_info(id_d_upper_max_sv, 'd_upper_max_sv', &
       [ character(len=error_message_length) :: 'n<1', 'A is not square.', &
       'Singular vector size error', 'Failure to converge.'])
  type(routine_info), parameter :: info_z_upper_max_sv= &
       routine_info(id_z_upper_max_sv, 'z_upper_max_sv', &
       [ character(len=error_message_length) :: 'n<1', 'A is not square.', &
       'Singular vector size error', 'Failure to converge.'])

  ! src/qr_factorization/qr_factorization 1500
  integer(int32), parameter :: mod_id_qr_factorization=1500
  integer(int32), parameter :: id_d_qr_of=mod_id_qr_factorization + 0
  integer(int32), parameter :: id_d_qr_bv_to_ub=mod_id_qr_factorization + 1
  integer(int32), parameter :: id_f_d_qr_bv_to_ub=mod_id_qr_factorization + 2
  integer(int32), parameter :: id_z_qr_of=mod_id_qr_factorization + 3
  integer(int32), parameter :: id_z_qr_bv_to_ub=mod_id_qr_factorization + 4
  integer(int32), parameter :: id_f_z_qr_bv_to_ub=mod_id_qr_factorization + 5

  type(routine_info), parameter :: info_d_qr_of=routine_info(id_d_qr_of, &
       'd_qr_of', [ character(len=error_message_length) :: '' ])

  type(routine_info), parameter :: info_d_qr_bv_to_ub=routine_info(id_d_qr_bv_to_ub, &
       'd_qr_bv_to_ub', &
       [ character(len=error_message_length) :: 'ub%n /= bv%n or sw%n /= ub%n', &
       'Not enough stroage for sweeps', 'Maxind or Minind out of bounds for sweeps.', &
       'Insufficient storage in bv', 'Insufficient storage in ub'])

  type(routine_info), parameter :: info_z_qr_of=routine_info(id_z_qr_of, &
       'z_qr_of', [ character(len=error_message_length) :: '' ])

  type(routine_info), parameter :: info_z_qr_bv_to_ub=routine_info(id_z_qr_bv_to_ub, &
       'z_qr_bv_to_ub', &
       [ character(len=error_message_length) :: 'ub%n /= bv%n or sw%n /= ub%n', &
       'Not enough stroage for sweeps', 'Maxind or Minind out of bounds for sweeps.', &
       'Insufficient storage in bv', 'Insufficient storage in ub' ])

  ! src/solve/solve 1600
  integer(int32), parameter :: mod_id_solve=1600
  integer(int32), parameter :: id_d_solve_ub=mod_id_solve + 0
  integer(int32), parameter :: id_d_back_solve_ub=mod_id_solve + 1
  integer(int32), parameter :: id_f_d_back_solve_ub=mod_id_solve + 2
  integer(int32), parameter :: id_z_solve_ub=mod_id_solve + 3
  integer(int32), parameter :: id_z_back_solve_ub=mod_id_solve + 4
  integer(int32), parameter :: id_f_z_back_solve_ub=mod_id_solve + 5
  integer(int32), parameter :: id_d_v_solve_ub=mod_id_solve + 6
  integer(int32), parameter :: id_d_v_back_solve_ub=mod_id_solve + 7
  integer(int32), parameter :: id_f_d_v_back_solve_ub=mod_id_solve + 8
  integer(int32), parameter :: id_z_v_solve_ub=mod_id_solve + 9
  integer(int32), parameter :: id_z_v_back_solve_ub=mod_id_solve + 10
  integer(int32), parameter :: id_f_z_v_back_solve_ub=mod_id_solve + 11

  type(routine_info), parameter :: info_d_solve_ub= &
       routine_info(id_d_solve_ub, &
       'd_solve_ub', &
       [ character(len=error_message_length) :: '' ])

  type(routine_info), parameter :: info_d_back_solve_ub= &
       routine_info(id_d_back_solve_ub, &
       'd_back_solve_ub', &
       [ character(len=error_message_length) :: 'ub%n /= size(c,1)', 'ub%lbw /= 0', 'n < 1', &
       'x is not the same size as c.' ])

  type(routine_info), parameter :: info_z_solve_ub= &
       routine_info(id_z_solve_ub, &
       'd_solve_ub', &
       [ character(len=error_message_length) :: '' ])

  type(routine_info), parameter :: info_z_back_solve_ub= &
       routine_info(id_z_back_solve_ub, &
       'd_back_solve_ub', &
       [ character(len=error_message_length) :: 'ub%n /= size(c,1)', 'ub%lbw /= 0', 'n < 1', &
       'x is not the same size as c.' ])

  type(routine_info), parameter :: info_d_v_solve_ub= &
       routine_info(id_d_v_solve_ub, &
       'd_v_solve_ub', &
       [ character(len=error_message_length) :: '' ])

  type(routine_info), parameter :: info_d_v_back_solve_ub= &
       routine_info(id_d_v_back_solve_ub, &
       'd_v_back_solve_ub', &
       [ character(len=error_message_length) :: 'ub%n /= size(c)', 'ub%lbw /= 0', 'n < 1', &
       'x is not the same size as c.' ])

  type(routine_info), parameter :: info_z_v_solve_ub= &
       routine_info(id_z_v_solve_ub, &
       'd_solve_ub', &
       [ character(len=error_message_length) :: '' ])

  type(routine_info), parameter :: info_z_v_back_solve_ub= &
       routine_info(id_z_v_back_solve_ub, &
       'd_back_solve_ub', &
       [ character(len=error_message_length) :: 'ub%n /= size(c)', 'ub%lbw /= 0', 'n < 1', &
       'x is not the same size as c.' ])

  integer(int32), parameter :: id_d_solve_bv=mod_id_solve + 12
  integer(int32), parameter :: id_d_forward_solve_bv=mod_id_solve + 13
  integer(int32), parameter :: id_f_d_forward_solve_bv=mod_id_solve + 14
  integer(int32), parameter :: id_d_v_solve_bv=mod_id_solve + 15
  integer(int32), parameter :: id_d_v_forward_solve_bv=mod_id_solve + 16
  integer(int32), parameter :: id_f_d_v_forward_solve_bv=mod_id_solve + 17
  integer(int32), parameter :: id_z_solve_bv=mod_id_solve + 18
  integer(int32), parameter :: id_z_forward_solve_bv=mod_id_solve + 19
  integer(int32), parameter :: id_f_z_forward_solve_bv=mod_id_solve + 20
  integer(int32), parameter :: id_z_v_solve_bv=mod_id_solve + 21
  integer(int32), parameter :: id_z_v_forward_solve_bv=mod_id_solve + 22
  integer(int32), parameter :: id_f_z_v_forward_solve_bv=mod_id_solve + 23

  type(routine_info), parameter :: info_d_solve_bv= &
       routine_info(id_d_solve_bv, &
       'd_solve_bv', &
       [ character(len=error_message_length) :: '' ])

  type(routine_info), parameter :: info_d_forward_solve_bv= &
       routine_info(id_d_forward_solve_bv, &
       'd_forward_solve_bv', &
       [ character(len=error_message_length) :: 'bv%n /= size(c,2)', 'bv%lbw /= 0', 'n < 1', &
       'x is not the same size as c.' ])

  type(routine_info), parameter :: info_z_solve_bv= &
       routine_info(id_z_solve_bv, &
       'd_solve_bv', &
       [ character(len=error_message_length) :: '' ])
  
  type(routine_info), parameter :: info_z_forward_solve_bv= &
       routine_info(id_z_forward_solve_bv, &
       'd_forward_solve_bv', &
       [ character(len=error_message_length) :: 'bv%n /= size(c,2)', 'bv%lbw /= 0', 'n < 1', &
       'x is not the same size as c.' ])

  type(routine_info), parameter :: info_d_v_solve_bv= &
       routine_info(id_d_v_solve_bv, &
       'd_v_solve_bv', &
       [ character(len=error_message_length) :: '' ])

  type(routine_info), parameter :: info_d_v_forward_solve_bv= &
       routine_info(id_d_v_forward_solve_bv, &
       'd_v_forward_solve_bv', &
       [ character(len=error_message_length) :: 'bv%n /= size(c,2)', 'bv%lbw /= 0', 'n < 1', &
       'x is not the same size as c.' ])

  type(routine_info), parameter :: info_z_v_solve_bv= &
       routine_info(id_z_v_solve_bv, &
       'd_solve_bv', &
       [ character(len=error_message_length) :: '' ])

  type(routine_info), parameter :: info_z_v_forward_solve_bv= &
       routine_info(id_z_v_forward_solve_bv, &
       'd_forward_solve_bv', &
       [ character(len=error_message_length) :: 'bv%n /= size(c,1)', 'bv%lbw /= 0', 'n < 1', &
       'x is not the same size as c.' ])

  ! src/row_compress 1700
  integer(int32), parameter :: mod_id_row_compress=1700
  integer(int32), parameter :: id_d_rc_of=mod_id_row_compress + 0  
  integer(int32), parameter :: id_d_row_compress=mod_id_row_compress + 1
  integer(int32), parameter :: id_f_d_row_compress=mod_id_row_compress + 2
  integer(int32), parameter :: id_z_rc_of=mod_id_row_compress + 3
  integer(int32), parameter :: id_z_row_compress=mod_id_row_compress + 4
  integer(int32), parameter :: id_f_z_row_compress=mod_id_row_compress + 5

  type(routine_info), parameter :: info_d_rc_of=routine_info(id_d_rc_of, &
       'd_rc_of', [ character(len=error_message_length) :: '' ])

  type(routine_info), parameter :: info_d_row_compress=routine_info(id_d_row_compress, &
       'd_row_compress', &
       [ character(len=error_message_length) :: 'ub%n /= bv%n or sw%n /= ub%n', &
       'Not enough stroage for each sweeps', 'Not enough storage for number of sweeps', &
       'Insufficient storage in ubt.', 'Insufficient storage in bv.'])

  type(routine_info), parameter :: info_z_rc_of=routine_info(id_z_rc_of, &
       'z_rc_of', [ character(len=error_message_length) :: '' ])

  type(routine_info), parameter :: info_z_row_compress=routine_info(id_z_row_compress, &
       'z_row_compress', &
       [ character(len=error_message_length) :: 'ub%n /= bv%n or sw%n /= ub%n', &
       'Not enough stroage for each sweeps', 'Not enough storage for number of sweeps', &
       'Insufficient storage in ubt.', 'Insufficient storage in bv.'])

  ! src/types/random 1800
  integer(int32), parameter :: mod_id_random=1800
  integer(int32), parameter :: id_d_random_bc=mod_id_random + 0
  integer(int32), parameter :: id_d_random_br=mod_id_random + 1
  integer(int32), parameter :: id_d_random_ub0=mod_id_random + 2
  integer(int32), parameter :: id_d_random_bv0=mod_id_random + 3
  integer(int32), parameter :: id_d_random_bt0=mod_id_random + 4
  integer(int32), parameter :: id_d_random_wb0=mod_id_random + 5
  integer(int32), parameter :: id_d_random_ubt0=mod_id_random + 6
  integer(int32), parameter :: id_d_random_wbv0=mod_id_random + 7

  integer(int32), parameter :: id_d_random_ub1=mod_id_random + 8
  integer(int32), parameter :: id_d_random_bv1=mod_id_random + 9
  integer(int32), parameter :: id_d_random_bt1=mod_id_random + 10
  integer(int32), parameter :: id_d_random_wb1=mod_id_random + 11
  integer(int32), parameter :: id_d_random_ubt1=mod_id_random + 12
  integer(int32), parameter :: id_d_random_wbv1=mod_id_random + 13
  

  integer(int32), parameter :: id_z_random_bc=mod_id_random + 14
  integer(int32), parameter :: id_z_random_br=mod_id_random + 15
  integer(int32), parameter :: id_z_random_ub0=mod_id_random + 16
  integer(int32), parameter :: id_z_random_bv0=mod_id_random + 17
  integer(int32), parameter :: id_z_random_bt0=mod_id_random + 18
  integer(int32), parameter :: id_z_random_wb0=mod_id_random + 19
  integer(int32), parameter :: id_z_random_ubt0=mod_id_random + 20
  integer(int32), parameter :: id_z_random_wbv0=mod_id_random + 21

  integer(int32), parameter :: id_z_random_ub1=mod_id_random + 22
  integer(int32), parameter :: id_z_random_bv1=mod_id_random + 23
  integer(int32), parameter :: id_z_random_bt1=mod_id_random + 24
  integer(int32), parameter :: id_z_random_wb1=mod_id_random + 25
  integer(int32), parameter :: id_z_random_ubt1=mod_id_random + 26
  integer(int32), parameter :: id_z_random_wbv1=mod_id_random + 27
  
  type(routine_info), parameter :: info_d_random_bc=routine_info(id_d_random_bc, &
       'd_random_bc', &
       [ character(len=error_message_length) :: 'n < 1', &
       'lbw > lbwmax', 'ubw > ubwmax'])

  type(routine_info), parameter :: info_d_random_br=routine_info(id_d_random_br, &
       'd_random_br', &
       [ character(len=error_message_length) :: 'n < 1', &
       'lbw > lbwmax', 'ubw > ubwmax'])

  type(routine_info), parameter :: info_d_random_ub0=routine_info(id_d_random_ub0, &
       'd_random_ub0', &
       [ character(len=error_message_length) :: 'n < 1', &
       'lbw > lbwmax', 'ubw > ubwmax'])

  type(routine_info), parameter :: info_d_random_bv0=routine_info(id_d_random_bv0, &
       'd_random_bv0', &
       [ character(len=error_message_length) :: 'n < 1', &
       'lbw > lbwmax', 'ubw > ubwmax'])

  type(routine_info), parameter :: info_d_random_bt0=routine_info(id_d_random_bt0, &
       'd_random_bt0', &
       [ character(len=error_message_length) :: 'n < 1', &
       'lbw > lbwmax', 'ubw > ubwmax'])

  type(routine_info), parameter :: info_d_random_wb0=routine_info(id_d_random_wb0, &
       'd_random_wb0', &
       [ character(len=error_message_length) :: 'n < 1', &
       'lbw > lbwmax', 'ubw > ubwmax'])

  type(routine_info), parameter :: info_d_random_ubt0=routine_info(id_d_random_ubt0, &
       'd_random_ubt0', &
       [ character(len=error_message_length) :: 'n < 1', &
       'lbw > lbwmax', 'ubw > ubwmax'])
  
  type(routine_info), parameter :: info_d_random_wbv0=routine_info(id_d_random_wbv0, &
       'd_random_wbv0', &
       [ character(len=error_message_length) :: 'n < 1', &
       'lbw > lbwmax', 'ubw > ubwmax'])

  ! array lbw ubw
  type(routine_info), parameter :: info_d_random_ub1=routine_info(id_d_random_ub1, &
       'd_random_ub1', &
       [ character(len=error_message_length) :: 'n < 1', &
       'lbw > lbwmax', 'ubw > ubwmax'])

  type(routine_info), parameter :: info_d_random_bv1=routine_info(id_d_random_bv1, &
       'd_random_bv1', &
       [ character(len=error_message_length) :: 'n < 1', &
       'lbw > lbwmax', 'ubw > ubwmax'])

  type(routine_info), parameter :: info_d_random_bt1=routine_info(id_d_random_bt1, &
       'd_random_bt1', &
       [ character(len=error_message_length) :: 'n < 1', &
       'lbw > lbwmax', 'ubw > ubwmax'])

  type(routine_info), parameter :: info_d_random_wb1=routine_info(id_d_random_wb1, &
       'd_random_wb1', &
       [ character(len=error_message_length) :: 'n < 1', &
       'lbw > lbwmax', 'ubw > ubwmax'])

  type(routine_info), parameter :: info_d_random_ubt1=routine_info(id_d_random_ubt1, &
       'd_random_ubt1', &
       [ character(len=error_message_length) :: 'n < 1', &
       'lbw > lbwmax', 'ubw > ubwmax'])
  
  type(routine_info), parameter :: info_d_random_wbv1=routine_info(id_d_random_wbv1, &
       'd_random_wbv1', &
       [ character(len=error_message_length) :: 'n < 1', &
       'lbw > lbwmax', 'ubw > ubwmax'])

  ! complex
  
  type(routine_info), parameter :: info_z_random_bc=routine_info(id_z_random_bc, &
       'z_random_bc', &
       [ character(len=error_message_length) :: 'n < 1', &
       'lbw > lbwmax', 'ubw > ubwmax'])

  type(routine_info), parameter :: info_z_random_br=routine_info(id_z_random_br, &
       'z_random_br', &
       [ character(len=error_message_length) :: 'n < 1', &
       'lbw > lbwmax', 'ubw > ubwmax'])

  type(routine_info), parameter :: info_z_random_ub0=routine_info(id_z_random_ub0, &
       'z_random_ub0', &
       [ character(len=error_message_length) :: 'n < 1', &
       'lbw > lbwmax', 'ubw > ubwmax'])

  type(routine_info), parameter :: info_z_random_bv0=routine_info(id_z_random_bv0, &
       'z_random_bv0', &
       [ character(len=error_message_length) :: 'n < 1', &
       'lbw > lbwmax', 'ubw > ubwmax'])

  type(routine_info), parameter :: info_z_random_bt0=routine_info(id_z_random_bt0, &
       'z_random_bt0', &
       [ character(len=error_message_length) :: 'n < 1', &
       'lbw > lbwmax', 'ubw > ubwmax'])

  type(routine_info), parameter :: info_z_random_wb0=routine_info(id_z_random_wb0, &
       'z_random_wb0', &
       [ character(len=error_message_length) :: 'n < 1', &
       'lbw > lbwmax', 'ubw > ubwmax'])

  type(routine_info), parameter :: info_z_random_ubt0=routine_info(id_z_random_ubt0, &
       'z_random_ubt0', &
       [ character(len=error_message_length) :: 'n < 1', &
       'lbw > lbwmax', 'ubw > ubwmax'])
  
  type(routine_info), parameter :: info_z_random_wbv0=routine_info(id_z_random_wbv0, &
       'z_random_wbv0', &
       [ character(len=error_message_length) :: 'n < 1', &
       'lbw > lbwmax', 'ubw > ubwmax'])

  ! array versions
  
  type(routine_info), parameter :: info_z_random_ub1=routine_info(id_z_random_ub1, &
       'z_random_ub1', &
       [ character(len=error_message_length) :: 'n < 1', &
       'lbw > lbwmax', 'ubw > ubwmax'])

  type(routine_info), parameter :: info_z_random_bv1=routine_info(id_z_random_bv1, &
       'z_random_bv1', &
       [ character(len=error_message_length) :: 'n < 1', &
       'lbw > lbwmax', 'ubw > ubwmax'])

  type(routine_info), parameter :: info_z_random_bt1=routine_info(id_z_random_bt1, &
       'z_random_bt1', &
       [ character(len=error_message_length) :: 'n < 1', &
       'lbw > lbwmax', 'ubw > ubwmax'])

  type(routine_info), parameter :: info_z_random_wb1=routine_info(id_z_random_wb1, &
       'z_random_wb1', &
       [ character(len=error_message_length) :: 'n < 1', &
       'lbw > lbwmax', 'ubw > ubwmax'])

  type(routine_info), parameter :: info_z_random_ubt1=routine_info(id_z_random_ubt1, &
       'z_random_ubt1', &
       [ character(len=error_message_length) :: 'n < 1', &
       'lbw > lbwmax', 'ubw > ubwmax'])
  
  type(routine_info), parameter :: info_z_random_wbv1=routine_info(id_z_random_wbv1, &
       'z_random_wbv1', &
       [ character(len=error_message_length) :: 'n < 1', &
       'lbw > lbwmax', 'ubw > ubwmax'])

  ! src/types/products 1900
  integer(int32), parameter :: mod_id_products=1900
  integer(int32), parameter :: id_d_ub_times_general=mod_id_products + 0
  integer(int32), parameter :: id_d_product_of_ub_and_general=mod_id_products + 1
  integer(int32), parameter :: id_z_ub_times_general=mod_id_products + 2
  integer(int32), parameter :: id_z_product_of_ub_and_general=mod_id_products + 3

  integer(int32), parameter :: id_d_bv_times_general=mod_id_products + 4
  integer(int32), parameter :: id_d_product_of_bv_and_general=mod_id_products + 5
  integer(int32), parameter :: id_z_bv_times_general=mod_id_products + 6
  integer(int32), parameter :: id_z_product_of_bv_and_general=mod_id_products + 7

  integer(int32), parameter :: id_d_wb_times_general=mod_id_products + 8
  integer(int32), parameter :: id_d_product_of_wb_and_general=mod_id_products + 9
  integer(int32), parameter :: id_z_wb_times_general=mod_id_products + 10
  integer(int32), parameter :: id_z_product_of_wb_and_general=mod_id_products + 11

  integer(int32), parameter :: id_d_bt_times_general=mod_id_products + 12
  integer(int32), parameter :: id_d_product_of_bt_and_general=mod_id_products + 13
  integer(int32), parameter :: id_z_bt_times_general=mod_id_products + 14
  integer(int32), parameter :: id_z_product_of_bt_and_general=mod_id_products + 15

  integer(int32), parameter :: id_d_ubt_times_general=mod_id_products + 16
  integer(int32), parameter :: id_d_product_of_ubt_and_general=mod_id_products + 17
  integer(int32), parameter :: id_z_ubt_times_general=mod_id_products + 18
  integer(int32), parameter :: id_z_product_of_ubt_and_general=mod_id_products + 19

  integer(int32), parameter :: id_d_wbv_times_general=mod_id_products + 20
  integer(int32), parameter :: id_d_product_of_wbv_and_general=mod_id_products + 21
  integer(int32), parameter :: id_z_wbv_times_general=mod_id_products + 22
  integer(int32), parameter :: id_z_product_of_wbv_and_general=mod_id_products + 23

  integer(int32), parameter :: id_d_general_times_ub=mod_id_products + 24
  integer(int32), parameter :: id_d_product_of_general_and_ub=mod_id_products + 25
  integer(int32), parameter :: id_z_general_times_ub=mod_id_products + 26
  integer(int32), parameter :: id_z_product_of_general_and_ub=mod_id_products + 27
  
  type(routine_info), parameter :: info_d_ub_times_general=routine_info(id_d_ub_times_general, &
       'd_ub_times_general', &
       [ character(len=error_message_length) :: 'ub%n<1', &
       'ub%n /= size(a,1)', 'size(a) /= size(c)'])

  type(routine_info), parameter :: info_d_product_of_ub_and_general= &
       routine_info(id_d_product_of_ub_and_general, &
       'd_product_of_ub_and_general', &
       [ character(len=error_message_length) :: 'ub%n<1', &
       'ub%n /= size(a,1)', 'size(a) /= size(c)'])

  type(routine_info), parameter :: info_z_ub_times_general=routine_info(id_z_ub_times_general, &
       'z_ub_times_general', &
       [ character(len=error_message_length) :: 'ub%n<1', &
       'ub%n /= size(a,1)', 'size(a) /= size(c)'])

  type(routine_info), parameter :: info_z_product_of_ub_and_general= &
       routine_info(id_z_product_of_ub_and_general, &
       'z_product_of_ub_and_general', &
       [ character(len=error_message_length) :: 'ub%n<1', &
       'ub%n /= size(a,1)', 'size(a) /= size(c)'])

  type(routine_info), parameter :: info_d_bv_times_general=routine_info(id_d_bv_times_general, &
       'd_bv_times_general', &
       [ character(len=error_message_length) :: 'bv%n<1', &
       'bv%n /= size(a,1)', 'size(a) /= size(c)'])

  type(routine_info), parameter :: info_d_product_of_bv_and_general= &
       routine_info(id_d_product_of_bv_and_general, &
       'd_product_of_bv_and_general', &
       [ character(len=error_message_length) :: 'bv%n<1', &
       'bv%n /= size(a,1)', 'size(a) /= size(c)'])

  type(routine_info), parameter :: info_z_bv_times_general=routine_info(id_z_bv_times_general, &
       'z_bv_times_general', &
       [ character(len=error_message_length) :: 'bv%n<1', &
       'bv%n /= size(a,1)', 'size(a) /= size(c)'])

  type(routine_info), parameter :: info_z_product_of_bv_and_general= &
       routine_info(id_z_product_of_bv_and_general, &
       'z_product_of_bv_and_general', &
       [ character(len=error_message_length) :: 'bv%n<1', &
       'bv%n /= size(a,1)', 'size(a) /= size(c)'])

  type(routine_info), parameter :: info_d_wb_times_general=routine_info(id_d_wb_times_general, &
       'd_wb_times_general', &
       [ character(len=error_message_length) :: 'wb%n<1', &
       'wb%n /= size(a,1)', 'size(a) /= size(c)'])

  type(routine_info), parameter :: info_d_product_of_wb_and_general= &
       routine_info(id_d_product_of_wb_and_general, &
       'd_product_of_wb_and_general', &
       [ character(len=error_message_length) :: 'wb%n<1', &
       'wb%n /= size(a,1)', 'size(a) /= size(c)'])

  type(routine_info), parameter :: info_z_wb_times_general=routine_info(id_z_wb_times_general, &
       'z_wb_times_general', &
       [ character(len=error_message_length) :: 'wb%n<1', &
       'wb%n /= size(a,1)', 'size(a) /= size(c)'])

  type(routine_info), parameter :: info_z_product_of_wb_and_general= &
       routine_info(id_z_product_of_wb_and_general, &
       'z_product_of_wb_and_general', &
       [ character(len=error_message_length) :: 'wb%n<1', &
       'wb%n /= size(a,1)', 'size(a) /= size(c)'])

  type(routine_info), parameter :: info_d_bt_times_general=routine_info(id_d_bt_times_general, &
       'd_bt_times_general', &
       [ character(len=error_message_length) :: 'bt%n<1', &
       'bt%n /= size(a,1)', 'size(a) /= size(c)'])

  type(routine_info), parameter :: info_d_product_of_bt_and_general= &
       routine_info(id_d_product_of_bt_and_general, &
       'd_product_of_bt_and_general', &
       [ character(len=error_message_length) :: 'bt%n<1', &
       'bt%n /= size(a,1)', 'size(a) /= size(c)'])

  type(routine_info), parameter :: info_z_bt_times_general=routine_info(id_z_bt_times_general, &
       'z_bt_times_general', &
       [ character(len=error_message_length) :: 'bt%n<1', &
       'bt%n /= size(a,1)', 'size(a) /= size(c)'])

  type(routine_info), parameter :: info_z_product_of_bt_and_general= &
       routine_info(id_z_product_of_bt_and_general, &
       'z_product_of_bt_and_general', &
       [ character(len=error_message_length) :: 'bt%n<1', &
       'bt%n /= size(a,1)', 'size(a) /= size(c)'])
  
  type(routine_info), parameter :: info_d_ubt_times_general=routine_info(id_d_ubt_times_general, &
       'd_ubt_times_general', &
       [ character(len=error_message_length) :: 'ubt%n<1', &
       'ubt%n /= size(a,1)', 'size(a) /= size(c)'])

  type(routine_info), parameter :: info_d_product_of_ubt_and_general= &
       routine_info(id_d_product_of_ubt_and_general, &
       'd_product_of_ubt_and_general', &
       [ character(len=error_message_length) :: 'ubt%n<1', &
       'ubt%n /= size(a,1)', 'size(a) /= size(c)'])

  type(routine_info), parameter :: info_z_ubt_times_general=routine_info(id_z_ubt_times_general, &
       'z_ubt_times_general', &
       [ character(len=error_message_length) :: 'ubt%n<1', &
       'ubt%n /= size(a,1)', 'size(a) /= size(c)'])

  type(routine_info), parameter :: info_z_product_of_ubt_and_general= &
       routine_info(id_z_product_of_ubt_and_general, &
       'z_product_of_ubt_and_general', &
       [ character(len=error_message_length) :: 'ubt%n<1', &
       'ubt%n /= size(a,1)', 'size(a) /= size(c)'])

  type(routine_info), parameter :: info_d_wbv_times_general=routine_info(id_d_wbv_times_general, &
       'd_wbv_times_general', &
       [ character(len=error_message_length) :: 'wbv%n<1', &
       'wbv%n /= size(a,1)', 'size(a) /= size(c)'])

  type(routine_info), parameter :: info_d_product_of_wbv_and_general= &
       routine_info(id_d_product_of_wbv_and_general, &
       'd_product_of_wbv_and_general', &
       [ character(len=error_message_length) :: 'wbv%n<1', &
       'wbv%n /= size(a,1)', 'size(a) /= size(c)'])

  type(routine_info), parameter :: info_z_wbv_times_general=routine_info(id_z_wbv_times_general, &
       'z_wbv_times_general', &
       [ character(len=error_message_length) :: 'wbv%n<1', &
       'wbv%n /= size(a,1)', 'size(a) /= size(c)'])

  type(routine_info), parameter :: info_z_product_of_wbv_and_general= &
       routine_info(id_z_product_of_wbv_and_general, &
       'z_product_of_wbv_and_general', &
       [ character(len=error_message_length) :: 'wbv%n<1', &
       'wbv%n /= size(a,1)', 'size(a) /= size(c)'])

  type(routine_info), parameter :: info_d_general_times_ub=routine_info(id_d_general_times_ub, &
       'd_general_times_ub', &
       [ character(len=error_message_length) :: 'ub%n<1', &
       'ub%n /= size(a,1)', 'size(a) /= size(c)'])

  type(routine_info), parameter :: info_d_product_of_general_and_ub= &
       routine_info(id_d_product_of_general_and_ub, &
       'd_product_of_general_and_ub', &
       [ character(len=error_message_length) :: 'ub%n<1', &
       'ub%n /= size(a,1)', 'size(a) /= size(c)'])

  type(routine_info), parameter :: info_z_general_times_ub=routine_info(id_z_general_times_ub, &
       'z_general_times_ub', &
       [ character(len=error_message_length) :: 'ub%n<1', &
       'ub%n /= size(a,1)', 'size(a) /= size(c)'])

  type(routine_info), parameter :: info_z_product_of_general_and_ub= &
       routine_info(id_z_product_of_general_and_ub, &
       'z_product_of_general_and_ub', &
       [ character(len=error_message_length) :: 'ub%n<1', &
       'ub%n /= size(a,1)', 'size(a) /= size(c)'])
  

  ! src/types/cond_orth_band 2000
  integer(int32), parameter :: mod_id_cond_orth_band=2000
  integer(int32), parameter :: id_d_ub_min_sv=mod_id_cond_orth_band + 0
  integer(int32), parameter :: id_d_ub_max_sv=mod_id_cond_orth_band + 1
  integer(int32), parameter :: id_z_ub_min_sv=mod_id_cond_orth_band + 2
  integer(int32), parameter :: id_z_ub_max_sv=mod_id_cond_orth_band + 3

  type(routine_info), parameter :: info_d_ub_min_sv= &
       routine_info(id_d_ub_min_sv, 'd_ub_min_sv', &
       [ character(len=error_message_length) :: 'n<1', &
       'Singular vector size error', 'Failure to converge.'])

  type(routine_info), parameter :: info_d_ub_max_sv= &
       routine_info(id_d_ub_max_sv, 'd_ub_max_sv', &
       [ character(len=error_message_length) :: 'n<1', &
       'Singular vector size error', 'Failure to converge.'])

  type(routine_info), parameter :: info_z_ub_min_sv= &
       routine_info(id_z_ub_min_sv, 'z_ub_min_sv', &
       [ character(len=error_message_length) :: 'n<1', &
       'Singular vector size error', 'Failure to converge.'])

  type(routine_info), parameter :: info_z_ub_max_sv= &
       routine_info(id_z_ub_max_sv, 'z_ub_max_sv', &
       [ character(len=error_message_length) :: 'n<1', &
       'Singular vector size error', 'Failure to converge.'])
  
  ! src/types/submatrix 2100
  integer(int32), parameter :: mod_id_submatrix=1900
  integer(int32), parameter :: id_d_ub_leading=mod_id_submatrix + 0
  integer(int32), parameter :: id_z_ub_leading=mod_id_submatrix + 1
  integer(int32), parameter :: id_d_bv_trailing=mod_id_submatrix + 2
  integer(int32), parameter :: id_z_bv_trailing=mod_id_submatrix + 3
  integer(int32), parameter :: id_d_bt_leading=mod_id_submatrix + 4
  integer(int32), parameter :: id_z_bt_leading=mod_id_submatrix + 5
  integer(int32), parameter :: id_d_wb_trailing=mod_id_submatrix + 6
  integer(int32), parameter :: id_z_wb_trailing=mod_id_submatrix + 7
  integer(int32), parameter :: id_d_ubt_leading=mod_id_submatrix + 8
  integer(int32), parameter :: id_z_ubt_leading=mod_id_submatrix + 9
  integer(int32), parameter :: id_d_wbv_trailing=mod_id_submatrix + 10
  integer(int32), parameter :: id_z_wbv_trailing=mod_id_submatrix + 11

  integer(int32), parameter :: id_d_ub_to_columns=mod_id_submatrix + 12
  integer(int32), parameter :: id_f_d_ub_to_columns=mod_id_submatrix + 13
  integer(int32), parameter :: id_d_columns_of_ub=mod_id_submatrix + 14
  integer(int32), parameter :: id_z_ub_to_columns=mod_id_submatrix + 15
  integer(int32), parameter :: id_f_z_ub_to_columns=mod_id_submatrix + 16
  integer(int32), parameter :: id_z_columns_of_ub=mod_id_submatrix + 17

  integer(int32), parameter :: id_d_bv_to_rows=mod_id_submatrix + 18
  integer(int32), parameter :: id_f_d_bv_to_rows=mod_id_submatrix + 19
  integer(int32), parameter :: id_d_rows_of_bv=mod_id_submatrix + 20
  integer(int32), parameter :: id_z_bv_to_rows=mod_id_submatrix + 21
  integer(int32), parameter :: id_f_z_bv_to_rows=mod_id_submatrix + 22
  integer(int32), parameter :: id_z_rows_of_bv=mod_id_submatrix + 23

  integer(int32), parameter :: id_d_wb_to_columns=mod_id_submatrix + 24
  integer(int32), parameter :: id_f_d_wb_to_columns=mod_id_submatrix + 25
  integer(int32), parameter :: id_d_columns_of_wb=mod_id_submatrix + 26
  integer(int32), parameter :: id_z_wb_to_columns=mod_id_submatrix + 27
  integer(int32), parameter :: id_f_z_wb_to_columns=mod_id_submatrix + 28
  integer(int32), parameter :: id_z_columns_of_wb=mod_id_submatrix + 29

  integer(int32), parameter :: id_d_bt_to_rows=mod_id_submatrix + 30
  integer(int32), parameter :: id_f_d_bt_to_rows=mod_id_submatrix + 31
  integer(int32), parameter :: id_d_rows_of_bt=mod_id_submatrix + 32
  integer(int32), parameter :: id_z_bt_to_rows=mod_id_submatrix + 33
  integer(int32), parameter :: id_f_z_bt_to_rows=mod_id_submatrix + 34
  integer(int32), parameter :: id_z_rows_of_bt=mod_id_submatrix + 35

  type(routine_info), parameter :: info_d_ub_leading=routine_info(id_d_ub_leading, &
       'd_ub_leading', &
       [ character(len=error_message_length) :: 'ub%n<1', &
       'lbwmaxl0 < lbw or ubwmaxu0 < ubw'])

  type(routine_info), parameter :: info_z_ub_leading=routine_info(id_z_ub_leading, &
       'z_ub_leading', &
       [ character(len=error_message_length) :: 'ub%n<1', &
       'lbwmaxl0 < lbw or ubwmaxu0 < ubw'])

  type(routine_info), parameter :: info_d_bv_trailing=routine_info(id_d_bv_trailing, &
       'd_bv_trailing', &
       [ character(len=error_message_length) :: 'ub%n<1', &
       'lbwmaxl0 < lbw or ubwmaxu0 < ubw'])

  type(routine_info), parameter :: info_z_bv_trailing=routine_info(id_z_bv_trailing, &
       'z_bv_trailing', &
       [ character(len=error_message_length) :: 'ub%n<1', &
       'lbwmaxl0 < lbw or ubwmaxu0 < ubw'])

  type(routine_info), parameter :: info_d_bt_leading=routine_info(id_d_bt_leading, &
       'd_bt_leading', &
       [ character(len=error_message_length) :: 'bt%n<1', &
       'lbwmaxl0 < lbw or ubwmaxu0 < ubw'])

  type(routine_info), parameter :: info_z_bt_leading=routine_info(id_z_bt_leading, &
       'z_bt_leading', &
       [ character(len=error_message_length) :: 'bt%n<1', &
       'lbwmaxl0 < lbw or ubwmaxu0 < ubw'])

  type(routine_info), parameter :: info_d_wb_trailing=routine_info(id_d_wb_trailing, &
       'd_wb_trailing', &
       [ character(len=error_message_length) :: 'ub%n<1', &
       'lbwmaxl0 < lbw or ubwmaxu0 < ubw'])

  type(routine_info), parameter :: info_z_wb_trailing=routine_info(id_z_wb_trailing, &
       'z_wb_trailing', &
       [ character(len=error_message_length) :: 'ub%n<1', &
       'lbwmaxl0 < lbw or ubwmaxu0 < ubw'])

  type(routine_info), parameter :: info_d_ubt_leading=routine_info(id_d_ubt_leading, &
       'd_ubt_leading', &
       [ character(len=error_message_length) :: 'ubt%n<1', &
       'lbwmaxl0 < lbw or ubwmaxu0 < ubw'])

  type(routine_info), parameter :: info_z_ubt_leading=routine_info(id_z_ubt_leading, &
       'z_ubt_leading', &
       [ character(len=error_message_length) :: 'ubt%n<1', &
       'lbwmaxl0 < lbw or ubwmaxu0 < ubw'])

  type(routine_info), parameter :: info_d_wbv_trailing=routine_info(id_d_wbv_trailing, &
       'd_wbv_trailing', &
       [ character(len=error_message_length) :: 'wbv%n<1', &
       'lbwmaxl0 < lbw or ubwmaxu0 < ubw'])

  type(routine_info), parameter :: info_z_wbv_trailing=routine_info(id_z_wbv_trailing, &
       'z_wbv_trailing', &
       [ character(len=error_message_length) :: 'wbv%n<1', &
       'lbwmaxl0 < lbw or ubwmaxu0 < ubw'])
  
  type(routine_info), parameter :: info_d_ub_to_columns=routine_info(id_d_ub_to_columns, &
       'd_ub_to_columns', [ character(len=error_message_length) :: 'Size error in A.' ] )

  type(routine_info), parameter :: info_z_ub_to_columns=routine_info(id_z_ub_to_columns, &
       'z_ub_to_columns', [ character(len=error_message_length) :: 'Size error in A.' ] )

  type(routine_info), parameter :: info_d_columns_of_ub=routine_info(id_d_columns_of_ub, &
       'd_columns_of_ub', [ character(len=error_message_length) :: '' ] )

  type(routine_info), parameter :: info_z_columns_of_ub=routine_info(id_z_columns_of_ub, &
       'z_columns_of_ub', [ character(len=error_message_length) :: '' ] )

  type(routine_info), parameter :: info_d_wb_to_columns=routine_info(id_d_wb_to_columns, &
       'd_wb_to_columns', [ character(len=error_message_length) :: 'Size error in A.' ] )

  type(routine_info), parameter :: info_z_wb_to_columns=routine_info(id_z_wb_to_columns, &
       'z_wb_to_columns', [ character(len=error_message_length) :: 'Size error in A.' ] )

  type(routine_info), parameter :: info_d_columns_of_wb=routine_info(id_d_columns_of_wb, &
       'd_columns_of_wb', [ character(len=error_message_length) :: '' ] )

  type(routine_info), parameter :: info_z_columns_of_wb=routine_info(id_z_columns_of_wb, &
       'z_columns_of_wb', [ character(len=error_message_length) :: '' ] )
  
  type(routine_info), parameter :: info_d_bv_to_rows=routine_info(id_d_bv_to_rows, &
       'd_bv_to_rows', [ character(len=error_message_length) :: 'Size error in A.' ] )
  
  type(routine_info), parameter :: info_z_bv_to_rows=routine_info(id_z_bv_to_rows, &
       'z_bv_to_rows', [ character(len=error_message_length) :: 'Size error in A.' ] )

  type(routine_info), parameter :: info_d_rows_of_bv=routine_info(id_d_rows_of_bv, &
       'd_rows_of_bv', [ character(len=error_message_length) :: '' ] )

  type(routine_info), parameter :: info_z_rows_of_bv=routine_info(id_z_rows_of_bv, &
       'z_rows_of_bv', [ character(len=error_message_length) :: '' ] )

  type(routine_info), parameter :: info_d_bt_to_rows=routine_info(id_d_bt_to_rows, &
       'd_bt_to_rows', [ character(len=error_message_length) :: 'Size error in A.' ] )
  
  type(routine_info), parameter :: info_z_bt_to_rows=routine_info(id_z_bt_to_rows, &
       'z_bt_to_rows', [ character(len=error_message_length) :: 'Size error in A.' ] )

  type(routine_info), parameter :: info_d_rows_of_bt=routine_info(id_d_rows_of_bt, &
       'd_rows_of_bt', [ character(len=error_message_length) :: '' ] )

  type(routine_info), parameter :: info_z_rows_of_bt=routine_info(id_z_rows_of_bt, &
       'z_rows_of_bt', [ character(len=error_message_length) :: '' ] )
  
contains

  subroutine push_id(info,err)
    type(error_info), intent(inout), optional :: err
    type(routine_info), intent(in) :: info    
    if (present(err)) then
       if (err%rix <= max_routines) then
          err%routines(err%rix)=info%routine_id
          err%rix=err%rix+1
       else
          call shift(err%routines,-1)
          err%routines(max_routines)=info%routine_id
          err%rix=max_routines+1
       end if
    end if
  end subroutine push_id

  subroutine set_error(code, info, err)
    type(error_info), intent(inout), optional :: err
    type(routine_info), intent(in) :: info
    integer(kind=int32), intent(in) :: code
    integer(kind=int32) :: j
    if (present(err)) then
       if (err%halt) then
          write(error_unit,*) info%error_messages(code)
          write(error_unit,*) "Routines called: "
          do j=1,err%rix-1
             write(error_unit,*) info_index(err%routines(j))%routine_name
          end do
          stop
       else
          err%code=code
       end if
    else
       write(error_unit,*) info%error_messages(code)
       write(error_unit,*) "Routine called: "
       write(error_unit,*) info%routine_name
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
       info_index(info_z_ub_to_general%routine_id)=info_z_ub_to_general
       info_index(info_d_general_of_ub%routine_id)=info_d_general_of_ub
       info_index(info_z_general_of_ub%routine_id)=info_z_general_of_ub

       info_index(info_d_bt_to_general%routine_id)=info_d_bt_to_general
       info_index(info_z_bt_to_general%routine_id)=info_z_bt_to_general
       info_index(info_d_general_of_bt%routine_id)=info_d_general_of_bt
       info_index(info_z_general_of_bt%routine_id)=info_z_general_of_bt

       info_index(info_d_ubt_to_general%routine_id)=info_d_ubt_to_general
       info_index(info_z_ubt_to_general%routine_id)=info_z_ubt_to_general
       info_index(info_d_general_of_ubt%routine_id)=info_d_general_of_ubt
       info_index(info_z_general_of_ubt%routine_id)=info_z_general_of_ubt

       info_index(info_d_bv_to_general%routine_id)=info_d_bv_to_general
       info_index(info_z_bv_to_general%routine_id)=info_z_bv_to_general
       info_index(info_d_general_of_bv%routine_id)=info_d_general_of_bv
       info_index(info_z_general_of_bv%routine_id)=info_z_general_of_bv

       info_index(info_d_wb_to_general%routine_id)=info_d_wb_to_general
       info_index(info_z_wb_to_general%routine_id)=info_z_wb_to_general
       info_index(info_d_general_of_wb%routine_id)=info_d_general_of_wb
       info_index(info_z_general_of_wb%routine_id)=info_z_general_of_wb

       info_index(info_d_wbv_to_general%routine_id)=info_d_wbv_to_general
       info_index(info_z_wbv_to_general%routine_id)=info_z_wbv_to_general
       info_index(info_d_general_of_wbv%routine_id)=info_d_general_of_wbv
       info_index(info_z_general_of_wbv%routine_id)=info_z_general_of_wbv

       ! convert_bv_to_ub
       info_index(info_d_convert_bv_to_ub%routine_id)=info_d_convert_bv_to_ub
       info_index(info_d_ub_of_bv%routine_id)=info_d_ub_of_bv
       info_index(info_z_convert_bv_to_ub%routine_id)=info_z_convert_bv_to_ub
       info_index(info_z_ub_of_bv%routine_id)=info_z_ub_of_bv

       ! convert_ub_to_bv
       info_index(info_d_convert_ub_to_bv%routine_id)=info_d_convert_ub_to_bv
       info_index(info_d_bv_of_ub%routine_id)=info_d_bv_of_ub
       info_index(info_z_convert_ub_to_bv%routine_id)=info_z_convert_ub_to_bv
       info_index(info_z_bv_of_ub%routine_id)=info_z_bv_of_ub
       
       ! convert_wb_to_bt
       info_index(info_d_convert_wb_to_bt%routine_id)=info_d_convert_wb_to_bt
       info_index(info_d_bt_of_wb%routine_id)=info_d_bt_of_wb
       info_index(info_z_convert_wb_to_bt%routine_id)=info_z_convert_wb_to_bt
       info_index(info_z_bt_of_wb%routine_id)=info_z_bt_of_wb

       ! convert_bt_to_wb
       info_index(info_d_convert_bt_to_wb%routine_id)=info_d_convert_bt_to_wb
       info_index(info_d_wb_of_bt%routine_id)=info_d_wb_of_bt
       info_index(info_z_convert_bt_to_wb%routine_id)=info_z_convert_bt_to_wb
       info_index(info_z_wb_of_bt%routine_id)=info_z_wb_of_bt

       ! convert_wbv_to_ubt
       info_index(info_d_convert_wbv_to_ubt%routine_id)=info_d_convert_wbv_to_ubt
       info_index(info_d_ubt_of_wbv%routine_id)=info_d_ubt_of_wbv
       info_index(info_z_convert_wbv_to_ubt%routine_id)=info_z_convert_wbv_to_ubt
       info_index(info_z_ubt_of_wbv%routine_id)=info_z_ubt_of_wbv
       
       ! convert_ubt_to_wbv
       info_index(info_d_convert_ubt_to_wbv%routine_id)=info_d_convert_ubt_to_wbv
       info_index(info_d_wbv_of_ubt%routine_id)=info_d_wbv_of_ubt
       info_index(info_z_convert_ubt_to_wbv%routine_id)=info_z_convert_ubt_to_wbv
       info_index(info_z_wbv_of_ubt%routine_id)=info_z_wbv_of_ubt
       
       ! general_bv

       info_index(info_d_general_to_bv%routine_id)=info_d_general_to_bv
       info_index(info_f_d_general_to_bv%routine_id)=info_f_d_general_to_bv
       info_index(info_f_d_general_bv%routine_id)=info_f_d_general_bv
       info_index(info_z_general_to_bv%routine_id)=info_z_general_to_bv
       info_index(info_f_z_general_to_bv%routine_id)=info_f_z_general_to_bv
       info_index(info_f_z_general_bv%routine_id)=info_f_z_general_bv

       ! general_ub

       info_index(info_d_general_to_ub%routine_id)=info_d_general_to_ub
       info_index(info_f_d_general_to_ub%routine_id)=info_f_d_general_to_ub
       info_index(info_f_d_general_ub%routine_id)=info_f_d_general_ub
       info_index(info_z_general_to_ub%routine_id)=info_z_general_to_ub
       info_index(info_f_z_general_to_ub%routine_id)=info_f_z_general_to_ub
       info_index(info_f_z_general_ub%routine_id)=info_f_z_general_ub

       ! general_bt

       info_index(info_d_general_to_bt%routine_id)=info_d_general_to_bt
       info_index(info_f_d_general_to_bt%routine_id)=info_f_d_general_to_bt
       info_index(info_f_d_general_bt%routine_id)=info_f_d_general_bt
       info_index(info_z_general_to_bt%routine_id)=info_z_general_to_bt
       info_index(info_f_z_general_to_bt%routine_id)=info_f_z_general_to_bt
       info_index(info_f_z_general_bt%routine_id)=info_f_z_general_bt

       ! general_wb

       info_index(info_d_general_to_wb%routine_id)=info_d_general_to_wb
       info_index(info_f_d_general_to_wb%routine_id)=info_f_d_general_to_wb
       info_index(info_f_d_general_wb%routine_id)=info_f_d_general_wb
       info_index(info_z_general_to_wb%routine_id)=info_z_general_to_wb
       info_index(info_f_z_general_to_wb%routine_id)=info_f_z_general_to_wb
       info_index(info_f_z_general_wb%routine_id)=info_f_z_general_wb

       ! general_ubt

       info_index(info_d_general_to_ubt%routine_id)=info_d_general_to_ubt
       info_index(info_f_d_general_to_ubt%routine_id)=info_f_d_general_to_ubt
       info_index(info_f_d_general_ubt%routine_id)=info_f_d_general_ubt
       info_index(info_z_general_to_ubt%routine_id)=info_z_general_to_ubt
       info_index(info_f_z_general_to_ubt%routine_id)=info_f_z_general_to_ubt
       info_index(info_f_z_general_ubt%routine_id)=info_f_z_general_ubt

       ! general_wbv

       info_index(info_d_general_to_wbv%routine_id)=info_d_general_to_wbv
       info_index(info_f_d_general_to_wbv%routine_id)=info_f_d_general_to_wbv
       info_index(info_f_d_general_wbv%routine_id)=info_f_d_general_wbv
       info_index(info_z_general_to_wbv%routine_id)=info_z_general_to_wbv
       info_index(info_f_z_general_to_wbv%routine_id)=info_f_z_general_to_wbv
       info_index(info_f_z_general_wbv%routine_id)=info_f_z_general_wbv

       ! gs
       info_index(info_d_extend_gs_rows%routine_id)=info_d_extend_gs_rows
       info_index(info_z_extend_gs_rows%routine_id)=info_z_extend_gs_rows
       info_index(info_d_extend_gs_columns%routine_id)=info_d_extend_gs_columns
       info_index(info_z_extend_gs_columns%routine_id)=info_z_extend_gs_columns
       
       ! cond
       info_index(info_d_lower_left_nullvec%routine_id)=info_d_lower_left_nullvec
       info_index(info_z_lower_left_nullvec%routine_id)=info_z_lower_left_nullvec
       info_index(info_d_lower_right_nullvec%routine_id)=info_d_lower_right_nullvec
       info_index(info_z_lower_right_nullvec%routine_id)=info_z_lower_right_nullvec

       info_index(info_d_lower_min_sv%routine_id)=info_d_lower_min_sv
       info_index(info_z_lower_min_sv%routine_id)=info_z_lower_min_sv
       info_index(info_d_upper_min_sv%routine_id)=info_d_upper_min_sv
       info_index(info_z_upper_min_sv%routine_id)=info_z_upper_min_sv

       info_index(info_d_lower_max_sv%routine_id)=info_d_lower_max_sv
       info_index(info_z_lower_max_sv%routine_id)=info_z_lower_max_sv
       info_index(info_d_upper_max_sv%routine_id)=info_d_upper_max_sv
       info_index(info_z_upper_max_sv%routine_id)=info_z_upper_max_sv
       
       ! qr_factorization

       info_index(info_d_qr_bv_to_ub%routine_id)=info_d_qr_bv_to_ub
       info_index(info_d_qr_of%routine_id)=info_d_qr_of
       info_index(info_z_qr_bv_to_ub%routine_id)=info_z_qr_bv_to_ub
       info_index(info_z_qr_of%routine_id)=info_z_qr_of
       
       ! solve
       info_index(info_d_solve_ub%routine_id)=info_d_solve_ub
       info_index(info_d_back_solve_ub%routine_id)=info_d_back_solve_ub
       info_index(info_z_solve_ub%routine_id)=info_z_solve_ub       
       info_index(info_z_back_solve_ub%routine_id)=info_z_back_solve_ub
       info_index(info_d_v_solve_ub%routine_id)=info_d_v_solve_ub
       info_index(info_d_v_back_solve_ub%routine_id)=info_d_v_back_solve_ub
       info_index(info_z_v_solve_ub%routine_id)=info_z_v_solve_ub
       info_index(info_z_v_back_solve_ub%routine_id)=info_z_v_back_solve_ub

       info_index(info_d_solve_bv%routine_id)=info_d_solve_bv
       info_index(info_d_forward_solve_bv%routine_id)=info_d_forward_solve_bv
       info_index(info_z_solve_bv%routine_id)=info_z_solve_bv       
       info_index(info_z_forward_solve_bv%routine_id)=info_z_forward_solve_bv
       info_index(info_d_v_solve_bv%routine_id)=info_d_v_solve_bv
       info_index(info_d_v_forward_solve_bv%routine_id)=info_d_v_forward_solve_bv
       info_index(info_z_v_solve_bv%routine_id)=info_z_v_solve_bv
       info_index(info_z_v_forward_solve_bv%routine_id)=info_z_v_forward_solve_bv

       ! row_compress
       info_index(info_d_rc_of%routine_id)=info_d_rc_of
       info_index(info_d_row_compress%routine_id)=info_d_row_compress
       info_index(info_z_rc_of%routine_id)=info_z_rc_of
       info_index(info_z_row_compress%routine_id)=info_z_row_compress

       ! random

       info_index(info_d_random_bc%routine_id)=info_d_random_bc
       info_index(info_d_random_br%routine_id)=info_d_random_br
       info_index(info_d_random_ub0%routine_id)=info_d_random_ub0
       info_index(info_d_random_bv0%routine_id)=info_d_random_bv0
       info_index(info_d_random_bt0%routine_id)=info_d_random_bt0
       info_index(info_d_random_wb0%routine_id)=info_d_random_wb0
       info_index(info_d_random_ubt0%routine_id)=info_d_random_ubt0
       info_index(info_d_random_wbv0%routine_id)=info_d_random_wbv0

       info_index(info_d_random_ub1%routine_id)=info_d_random_ub1
       info_index(info_d_random_bv1%routine_id)=info_d_random_bv1
       info_index(info_d_random_bt1%routine_id)=info_d_random_bt1
       info_index(info_d_random_wb1%routine_id)=info_d_random_wb1
       info_index(info_d_random_ubt1%routine_id)=info_d_random_ubt1
       info_index(info_d_random_wbv1%routine_id)=info_d_random_wbv1
       
       info_index(info_z_random_bc%routine_id)=info_z_random_bc
       info_index(info_z_random_br%routine_id)=info_z_random_br
       info_index(info_z_random_ub0%routine_id)=info_z_random_ub0
       info_index(info_z_random_bv0%routine_id)=info_z_random_bv0
       info_index(info_z_random_bt0%routine_id)=info_z_random_bt0
       info_index(info_z_random_wb0%routine_id)=info_z_random_wb0
       info_index(info_z_random_ubt0%routine_id)=info_z_random_ubt0
       info_index(info_z_random_wbv0%routine_id)=info_z_random_wbv0

       info_index(info_z_random_ub1%routine_id)=info_z_random_ub1
       info_index(info_z_random_bv1%routine_id)=info_z_random_bv1
       info_index(info_z_random_bt1%routine_id)=info_z_random_bt1
       info_index(info_z_random_wb1%routine_id)=info_z_random_wb1
       info_index(info_z_random_ubt1%routine_id)=info_z_random_ubt1
       info_index(info_z_random_wbv1%routine_id)=info_z_random_wbv1

       ! products.f90

       info_index(info_d_ub_times_general%routine_id)=info_d_ub_times_general
       info_index(info_d_product_of_ub_and_general%routine_id)=info_d_product_of_ub_and_general
       info_index(info_z_ub_times_general%routine_id)=info_z_ub_times_general
       info_index(info_z_product_of_ub_and_general%routine_id)=info_z_product_of_ub_and_general

       info_index(info_d_bv_times_general%routine_id)=info_d_bv_times_general
       info_index(info_d_product_of_bv_and_general%routine_id)=info_d_product_of_bv_and_general
       info_index(info_z_bv_times_general%routine_id)=info_z_bv_times_general
       info_index(info_z_product_of_bv_and_general%routine_id)=info_z_product_of_bv_and_general

       info_index(info_d_wb_times_general%routine_id)=info_d_wb_times_general
       info_index(info_d_product_of_wb_and_general%routine_id)=info_d_product_of_wb_and_general
       info_index(info_z_wb_times_general%routine_id)=info_z_wb_times_general
       info_index(info_z_product_of_wb_and_general%routine_id)=info_z_product_of_wb_and_general

       info_index(info_d_bt_times_general%routine_id)=info_d_bt_times_general
       info_index(info_d_product_of_bt_and_general%routine_id)=info_d_product_of_bt_and_general
       info_index(info_z_bt_times_general%routine_id)=info_z_bt_times_general
       info_index(info_z_product_of_bt_and_general%routine_id)=info_z_product_of_bt_and_general

       ! condition estimation for orth. band types.
       
       info_index(info_d_ub_min_sv%routine_id)=info_d_ub_min_sv
       info_index(info_z_ub_min_sv%routine_id)=info_z_ub_min_sv
       info_index(info_d_ub_max_sv%routine_id)=info_d_ub_max_sv
       info_index(info_z_ub_max_sv%routine_id)=info_z_ub_max_sv

       ! submatrix

       info_index(info_d_ub_leading%routine_id)=info_d_ub_leading
       info_index(info_z_ub_leading%routine_id)=info_z_ub_leading       
       info_index(info_d_bv_trailing%routine_id)=info_d_bv_trailing
       info_index(info_z_bv_trailing%routine_id)=info_z_bv_trailing
       info_index(info_d_bt_leading%routine_id)=info_d_bt_leading
       info_index(info_z_bt_leading%routine_id)=info_z_bt_leading       
       info_index(info_d_wb_trailing%routine_id)=info_d_wb_trailing
       info_index(info_z_wb_trailing%routine_id)=info_z_wb_trailing
       info_index(info_d_ubt_leading%routine_id)=info_d_ubt_leading
       info_index(info_z_ubt_leading%routine_id)=info_z_ubt_leading       
       info_index(info_d_wbv_trailing%routine_id)=info_d_wbv_trailing
       info_index(info_z_wbv_trailing%routine_id)=info_z_wbv_trailing

       info_index(info_d_ub_to_columns%routine_id)=info_d_ub_to_columns
       info_index(info_z_ub_to_columns%routine_id)=info_z_ub_to_columns
       info_index(info_d_columns_of_ub%routine_id)=info_d_columns_of_ub
       info_index(info_z_columns_of_ub%routine_id)=info_z_columns_of_ub

       info_index(info_d_bt_to_rows%routine_id)=info_d_bt_to_rows
       info_index(info_z_bt_to_rows%routine_id)=info_z_bt_to_rows
       info_index(info_d_rows_of_bt%routine_id)=info_d_rows_of_bt
       info_index(info_z_rows_of_bt%routine_id)=info_z_rows_of_bt

       info_index(info_d_bv_to_rows%routine_id)=info_d_bv_to_rows
       info_index(info_z_bv_to_rows%routine_id)=info_z_bv_to_rows
       info_index(info_d_rows_of_bv%routine_id)=info_d_rows_of_bv
       info_index(info_z_rows_of_bv%routine_id)=info_z_rows_of_bv

       info_index(info_d_wb_to_columns%routine_id)=info_d_wb_to_columns
       info_index(info_z_wb_to_columns%routine_id)=info_z_wb_to_columns
       info_index(info_d_columns_of_wb%routine_id)=info_d_columns_of_wb
       info_index(info_z_columns_of_wb%routine_id)=info_z_columns_of_wb

    end if
  end subroutine initialize_errors

end module mod_error_id
