FC = gfortran
OBJDIR = ../mod
PROFLFLAGS = -pg
PROFCFLAGS = -pg
LFLAGS = 
CFLAGS = -Wall -fbounds-check -mcmodel=medium
#CFLAGS = -Wall -O2 -mcmodel=medium


OBJS =  $(OBJDIR)/utility.o $(OBJDIR)/shift.o $(OBJDIR)/nullvec.o $(OBJDIR)/prec.o $(OBJDIR)/error_id.o \
	$(OBJDIR)/rotation.o $(OBJDIR)/general_ub.o $(OBJDIR)/general_bv.o $(OBJDIR)/triangular.o \
	$(OBJDIR)/gs.o $(OBJDIR)/assemble.o $(OBJDIR)/band_types.o $(OBJDIR)/conversions_ub_to_bv.o \
	$(OBJDIR)/conversions_bv_to_ub.o $(OBJDIR)/nested_types.o \
	$(OBJDIR)/compressions_ub_to_bv.o $(OBJDIR)/compressions_bv_to_ub.o $(OBJDIR)/qr_iteration.o \
	$(OBJDIR)/qr_factorization.o $(OBJDIR)/sweeps.o $(OBJDIR)/back_substitution.o

TEST_EXECS =  $(OBJDIR)/test_decomp $(OBJDIR)/test_convert $(OBJDIR)/test_compress \
	$(OBJDIR)/test_qr_iteration $(OBJDIR)/test_qr_factorization $(OBJDIR)/test_solve
