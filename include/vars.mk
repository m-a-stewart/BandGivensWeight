FC = gfortran
OBJDIR = ../mod
PROFLFLAGS = -pg
PROFCFLAGS = -pg
LFLAGS = 
CFLAGS = -fbounds-check -mcmodel=medium -Waliasing -Wampersand -Wconversion -Wsurprising \
	-Wintrinsics-std -Wno-tabs -Wintrinsic-shadow -Wline-truncation -Wreal-q-constant -Wunused
#CFLAGS = -Wall -O2 -mcmodel=medium

MISCOBJS = $(OBJDIR)/misc.o $(OBJDIR)/prec.o $(OBJDIR)/utility.o $(OBJDIR)/error_id.o

TYPESOBJS = $(OBJDIR)/types.o $(OBJDIR)/band_types.o $(OBJDIR)/nested_types.o

TRANSFORMSOBJS = $(OBJDIR)/transforms.o $(OBJDIR)/rotation.o $(OBJDIR)/shift.o $(OBJDIR)/sweeps.o

ORTHOBJS = $(OBJDIR)/triangular.o $(OBJDIR)/nullvec.o $(OBJDIR)/gs.o

DECOMPOBJS = $(OBJDIR)/decomp.o $(OBJDIR)/general_bv.o $(OBJDIR)/general_ub.o

CONVERSIONOBJS = $(OBJDIR)/conversions.o $(OBJDIR)/conversions_ub_to_bv.o \
	$(OBJDIR)/conversions_bv_to_ub.o

COMPRESSIONOBJS = $(OBJDIR)/compressions.o $(OBJDIR)/compressions_ub_to_bv.o \
	$(OBJDIR)/compressions_bv_to_ub.o


OBJS =  $(MISCOBJS) $(TYPESOBJS) $(TRANSFORMSOBJS) $(ORTHOBJS) $(DECOMPOBJS) \
	$(CONVERSIONOBJS) $(OBJDIR)/assemble.o $(COMPRESSIONOBJS) \
	$(OBJDIR)/qr_iteration.o $(OBJDIR)/qr_factorization.o $(OBJDIR)/substitution.o \
	$(OBJDIR)/nested.o $(OBJDIR)/update.o

TEST_EXECS =  $(OBJDIR)/test_decomp $(OBJDIR)/test_convert $(OBJDIR)/test_compress \
	$(OBJDIR)/test_qr_iteration $(OBJDIR)/test_qr_factorization $(OBJDIR)/test_solve \
	$(OBJDIR)/test_sweeps $(OBJDIR)/test_update
