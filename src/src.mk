# module dependencies

_ORRB = general assemble products qr_factorization solve \
	row_compress transforms convert types misc orth

_MISC = error_id utility

_ERROR_ID = prec shift

_UTILITY = prec

_TYPES = band_types orth_band_types random
_BAND_TYPES = rotation utility shift prec
_ORTH_BAND_TYPES = prec utility band_types
_RANDOM = prec utility band_types orth_band_types
_ASSEMBLE = prec error_id band_types orth_band_types rotation
_PRODUCTS = prec error_id band_types orth_band_types rotation

_TRANSFORMS = sweeps rotation shift
_ROTATION = prec
_SWEEPS = band_types prec utility
_SHIFT = prec

_ORTH = triangular gs cond
_TRIANGULAR = prec
_COND = prec error_id utility triangular
_GS = utility prec error_id


_GENERAL = general_bv general_ub general_bt \
	general_wb general_ubt general_wbv
_GENERAL_BV = prec utility error_id cond gs shift rotation \
	orth_band_types band_types
_GENERAL_UB = prec utility error_id cond gs shift rotation \
	orth_band_types band_types
_GENERAL_BT = prec error_id utility \
	orth_band_types band_types general_ub
_GENERAL_WB = prec error_id utility \
	orth_band_types band_types general_bv
_GENERAL_UBT = prec error_id utility \
	orth_band_types band_types general_ub
_GENERAL_WBV = prec error_id utility \
	orth_band_types band_types general_bv


_CONVERT = convert_bv_to_ub convert_ub_to_bv \
	convert_wb_to_bt convert_bt_to_wb \
	convert_wbv_to_ubt convert_ubt_to_wbv 
_CONVERT_BV_TO_UB = prec error_id \
	rotation orth_band_types band_types
_CONVERT_UB_TO_BV = prec error_id shift \
	rotation orth_band_types band_types
_CONVERT_WB_TO_BT = prec error_id \
	rotation orth_band_types band_types
_CONVERT_BT_TO_WB = prec error_id shift \
	rotation orth_band_types band_types
_CONVERT_WBV_TO_UBT = prec error_id shift \
	rotation orth_band_types band_types
_CONVERT_UBT_TO_WBV = prec error_id shift \
	rotation orth_band_types band_types

_ROW_COMPRESS = prec error_id band_types orth_band_types \
	sweeps rotation shift

_QR_FACTORIZATION = prec error_id \
	convert_bv_to_ub convert_ub_to_bv \
	sweeps shift rotation \
	band_types orth_band_types

_SOLVE = prec error_id orth_band_types \
	band_types rotation

ORRB=$(patsubst %,$(OBJDIR)/%.o,$(_ORRB))
MISC=$(patsubst %,$(OBJDIR)/%.o,$(_MISC))
ERROR_ID=$(patsubst %,$(OBJDIR)/%.o,$(_ERROR_ID))
UTILITY=$(patsubst %,$(OBJDIR)/%.o,$(_UTILITY))
TYPES=$(patsubst %,$(OBJDIR)/%.o,$(_TYPES))
BAND_TYPES=$(patsubst %,$(OBJDIR)/%.o,$(_BAND_TYPES))
ORTH_BAND_TYPES=$(patsubst %,$(OBJDIR)/%.o,$(_ORTH_BAND_TYPES))
RANDOM=$(patsubst %,$(OBJDIR)/%.o,$(_RANDOM))
TRANSFORMS=$(patsubst %,$(OBJDIR)/%.o,$(_TRANSFORMS))
ROTATION=$(patsubst %,$(OBJDIR)/%.o,$(_ROTATION))
SWEEPS=$(patsubst %,$(OBJDIR)/%.o,$(_SWEEPS))
SHIFT=$(patsubst %,$(OBJDIR)/%.o,$(_SHIFT))
ORTH=$(patsubst %,$(OBJDIR)/%.o,$(_ORTH))
TRIANGULAR=$(patsubst %,$(OBJDIR)/%.o,$(_TRIANGULAR))
COND=$(patsubst %,$(OBJDIR)/%.o,$(_COND))
GS=$(patsubst %,$(OBJDIR)/%.o,$(_GS))
ASSEMBLE=$(patsubst %,$(OBJDIR)/%.o,$(_ASSEMBLE))
PRODUCTS=$(patsubst %,$(OBJDIR)/%.o,$(_PRODUCTS))
GENERAL=$(patsubst %,$(OBJDIR)/%.o,$(_GENERAL))
GENERAL_BV=$(patsubst %,$(OBJDIR)/%.o,$(_GENERAL_BV))
GENERAL_UB=$(patsubst %,$(OBJDIR)/%.o,$(_GENERAL_UB))
GENERAL_BT=$(patsubst %,$(OBJDIR)/%.o,$(_GENERAL_BT))
GENERAL_WB=$(patsubst %,$(OBJDIR)/%.o,$(_GENERAL_WB))
GENERAL_UBT=$(patsubst %,$(OBJDIR)/%.o,$(_GENERAL_UBT))
GENERAL_WBV=$(patsubst %,$(OBJDIR)/%.o,$(_GENERAL_WBV))
CONVERT=$(patsubst %,$(OBJDIR)/%.o,$(_CONVERT))
CONVERT_BV_TO_UB=$(patsubst %,$(OBJDIR)/%.o,$(_CONVERT_BV_TO_UB))
CONVERT_UB_TO_BV=$(patsubst %,$(OBJDIR)/%.o,$(_CONVERT_UB_TO_BV))
CONVERT_WB_TO_BT=$(patsubst %,$(OBJDIR)/%.o,$(_CONVERT_WB_TO_BT))
CONVERT_BT_TO_WB=$(patsubst %,$(OBJDIR)/%.o,$(_CONVERT_BT_TO_WB))
CONVERT_WBV_TO_UBT=$(patsubst %,$(OBJDIR)/%.o,$(_CONVERT_WBV_TO_UBT))
CONVERT_UBT_TO_WBV=$(patsubst %,$(OBJDIR)/%.o,$(_CONVERT_UBT_TO_WBV))
ROW_COMPRESS=$(patsubst %,$(OBJDIR)/%.o,$(_ROW_COMPRESS))
QR_FACTORIZATION=$(patsubst %,$(OBJDIR)/%.o,$(_QR_FACTORIZATION))
SOLVE=$(patsubst %,$(OBJDIR)/%.o,$(_SOLVE))

$(OBJDIR)/orrb.o : $(ORRB)

$(OBJDIR)/misc.o : $(MISC)

$(OBJDIR)/error_id.o : $(ERROR_ID)
$(OBJDIR)/utility.o : $(UTILITY)

$(OBJDIR)/types.o : $(TYPES)
$(OBJDIR)/band_types.o : $(BAND_TYPES)
$(OBJDIR)/orth_band_types.o : $(ORTH_BAND_TYPES)
$(OBJDIR)/random.o : $(RANDOM)

$(OBJDIR)/transforms.o :$(TRANSFORMS)
$(OBJDIR)/rotation.o :$(ROTATION)
$(OBJDIR)/sweeps.o : $(SWEEPS)
$(OBJDIR)/shift.o : $(SHIFT)


$(OBJDIR)/orth.o : $(ORTH)
$(OBJDIR)/triangular.o : $(TRIANGULAR)
$(OBJDIR)/cond.o : $(COND)
$(OBJDIR)/gs.o : $(GS)


$(OBJDIR)/assemble.o : $(ASSEMBLE)
$(OBJDIR)/products.o : $(PRODUCTS)

$(OBJDIR)/general.o : $(GENERAL)
$(OBJDIR)/general_bv.o : $(GENERAL_BV)
$(OBJDIR)/general_ub.o : $(GENERAL_UB)
$(OBJDIR)/general_bt.o : $(GENERAL_BT)
$(OBJDIR)/general_wb.o : $(GENERAL_WB)
$(OBJDIR)/general_ubt.o : $(GENERAL_UBT)
$(OBJDIR)/general_wbv.o : $(GENERAL_WBV)

$(OBJDIR)/convert.o : $(CONVERT)
$(OBJDIR)/convert_bv_to_ub.o : $(CONVERT_BV_TO_UB)
$(OBJDIR)/convert_ub_to_bv.o : $(CONVERT_UB_TO_BV)
$(OBJDIR)/convert_wb_to_bt.o : $(CONVERT_WB_TO_BT)
$(OBJDIR)/convert_bt_to_wb.o : $(CONVERT_BT_TO_WB)
$(OBJDIR)/convert_wbv_to_ubt.o : $(CONVERT_WBV_TO_UBT)
$(OBJDIR)/convert_ubt_to_wbv.o : $(CONVERT_UBT_TO_WBV)

$(OBJDIR)/row_compress.o : $(ROW_COMPRESS)

$(OBJDIR)/qr_factorization.o : $(QR_FACTORIZATION)

$(OBJDIR)/solve.o : $(SOLVE)

# List of object files for linking

OBJS=$(patsubst %,$(OBJDIR)/%.o,$(_OBJS))

_MISCOBJS = misc prec utility error_id

_TYPESOBJS = types band_types orth_band_types \
	random assemble products

_TRANSFORMSOBJS = transforms rotation shift \
	sweeps

_ORTHOBJS = triangular cond gs

_GENERALOBJS = general general_bv general_ub \
	general_bt general_wb general_ubt \
	general_wbv

_CONVERTOBJS = convert convert_ub_to_bv \
	convert_bv_to_ub convert_wb_to_bt \
	convert_bt_to_wb convert_wbv_to_ubt \
	convert_ubt_to_wbv 

_OBJS =  $(_MISCOBJS) $(_TYPESOBJS) $(_TRANSFORMSOBJS) $(_ORTHOBJS) $(_GENERALOBJS) \
	$(_CONVERTOBJS) qr_factorization solve orrb row_compress
