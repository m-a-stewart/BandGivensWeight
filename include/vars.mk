FC = gfortran
OBJDIR = ../mod
PROFLFLAGS = -pg
PROFCFLAGS = -pg
LFLAGS = 
CFLAGS = -Wall -fbounds-check
#CFLAGS = -Wall -O2


OBJS =  $(OBJDIR)/utility.o $(OBJDIR)/shift.o $(OBJDIR)/nullvec.o $(OBJDIR)/prec.o \
	$(OBJDIR)/rotation.o $(OBJDIR)/general_ub.o $(OBJDIR)/general_bv.o $(OBJDIR)/triangular.o \
	$(OBJDIR)/gs.o $(OBJDIR)/assemble.o $(OBJDIR)/band_types.o $(OBJDIR)/convert_ub.o \
	$(OBJDIR)/convert_bv.o $(OBJDIR)/nested_types.o
