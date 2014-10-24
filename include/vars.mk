FC = gfortran
OBJDIR = ../mod
PROFLFLAGS = -pg
PROFCFLAGS = -pg
LFLAGS = 
CFLAGS = -fbounds-check -mcmodel=medium -Waliasing -Wampersand -Wconversion -Wsurprising \
	-Wintrinsics-std -Wno-tabs -Wintrinsic-shadow -Wline-truncation -Wreal-q-constant -Wunused -O2
#CFLAGS = -Wall -O2 -mcmodel=medium

MISCOBJS = $(OBJDIR)/misc.o $(OBJDIR)/prec.o $(OBJDIR)/utility.o $(OBJDIR)/error_id.o

TYPESOBJS = $(OBJDIR)/types.o $(OBJDIR)/band_types.o $(OBJDIR)/nested_types.o

TRANSFORMSOBJS = $(OBJDIR)/transforms.o $(OBJDIR)/rotation.o $(OBJDIR)/shift.o $(OBJDIR)/sweeps1.o

ORTHOBJS = $(OBJDIR)/triangular.o $(OBJDIR)/nullvec.o $(OBJDIR)/gs.o

GENERALOBJS = $(OBJDIR)/general.o $(OBJDIR)/general_bv.o $(OBJDIR)/general_ub.o \
	$(OBJDIR)/general_bt.o $(OBJDIR)/general_wb.o $(OBJDIR)/general_ubt.o \
	$(OBJDIR)/general_wbv.o

CONVERSIONOBJS = $(OBJDIR)/conversion.o $(OBJDIR)/conversion_ub_to_bv.o \
	$(OBJDIR)/conversion_bv_to_ub.o $(OBJDIR)/conversion_wb_to_bt.o \
	$(OBJDIR)/conversion_bt_to_wb.o $(OBJDIR)/conversion_wbv_to_ubt.o \
	$(OBJDIR)/conversion_ubt_to_wbv.o 

RECOMPRESSIONOBJS = $(OBJDIR)/recompression.o $(OBJDIR)/recompression_ub_to_bv.o \
	$(OBJDIR)/recompression_bv_to_ub.o


OBJS =  $(MISCOBJS) $(TYPESOBJS) $(TRANSFORMSOBJS) $(ORTHOBJS) $(GENERALOBJS) \
	$(CONVERSIONOBJS) $(OBJDIR)/assemble.o $(RECOMPRESSIONOBJS) \
	$(OBJDIR)/qr_iteration.o $(OBJDIR)/qr_factorization.o $(OBJDIR)/substitution.o \
	$(OBJDIR)/nested.o $(OBJDIR)/update.o $(OBJDIR)/row_compression.o

TEST_EXECS =  $(OBJDIR)/test_general $(OBJDIR)/test_general_bt $(OBJDIR)/test_general_wb \
	$(OBJDIR)/test_general_ubt $(OBJDIR)/test_general_wbv \
	$(OBJDIR)/test_convert_ub_and_bv $(OBJDIR)/test_convert_wb_and_bt \
	$(OBJDIR)/test_convert_ubt_and_wbv $(OBJDIR)/test_recompress \
	$(OBJDIR)/test_qr_iteration $(OBJDIR)/test_qr_factorization $(OBJDIR)/test_solve \
	$(OBJDIR)/test_sweeps1 $(OBJDIR)/test_update
