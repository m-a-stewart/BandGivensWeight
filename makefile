OBJDIR = mod
BINDIR = bin
SRCDIR = src
TESTDIR = test
SRCDIRS = src src/convert src/general src/misc src/orth \
	src/solve src/transforms src/types
FC = gfortran
LFLAGS =
#CFLAGS = -fbounds-check -mcmodel=medium -Waliasing -Wampersand -Wconversion -Wsurprising \
#	-Wintrinsics-std -Wno-tabs -Wintrinsic-shadow -Wline-truncation -Wreal-q-constant -Wunused
CFLAGS = -Wall -Wno-maybe-uninitialized -O2 -mcmodel=medium
#PROFLFLAGS = -pg
#PROFCFLAGS = -pg
PROFLFLAGS =
PROFCFLAGS =

vpath %.f90 $(SRCDIRS) test

.PHONY : all
all : $(OBJDIR)/orb.o

include $(SRCDIR)/src.mk

include $(TESTDIR)/test.mk

$(OBJDIR)/%.o : %.f90
	$(FC) $(CFLAGS) $(PROFCFLAGS) -I$(OBJDIR) -J$(OBJDIR) -c $< -o $@

$(BINDIR)/% : $(OBJDIR)/%.o
	$(FC) $(PROFLFLAGS) $(LFLAGS) $< $(OBJS) $(OBJDIR)/test_data.o -o $@

.PHONY : clean
clean :
	rm $(OBJDIR)/* $(BINDIR)/*

# prof : test_general
# 	gprof -c -q --no-time='__general_ub_MOD_f_d_upper_to_ub' \
# 		-no-time='__assemble_MOD_d_bv_to_upper' test_general | less

