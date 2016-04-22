OBJDIR = mod
BINDIR = bin
SRCDIR = src
TESTDIR = test
EXPDIR = exp
SRCDIRS = src src/convert src/general src/misc src/orth \
	src/solve src/transforms src/types
FC = gfortran
LFLAGS =
CFLAGS = -Wall -Wno-maybe-uninitialized -O2 -mcmodel=medium -fno-range-check
#CFLAGS = -Wall -Wno-maybe-uninitialized -mcmodel=medium -fbounds-check -fno-range-check
#PROFLFLAGS = -pg
#PROFCFLAGS = -pg
VIEWER = evince

vpath %.f90 $(SRCDIRS) test exp

.PHONY : all
all : $(OBJDIR)/orrb.o

notes.pdf : notes.md
	pandoc -o notes.pdf notes.md

.PHONY : view_notes
view_notes : notes.pdf
	$(VIEWER) notes.pdf

include $(SRCDIR)/src.mk
include $(TESTDIR)/test.mk
include $(EXPDIR)/exp.mk

$(OBJDIR)/%.o : %.f90
	$(FC) $(CFLAGS) $(PROFCFLAGS) -I$(OBJDIR) -J$(OBJDIR) -c $< -o $@

$(BINDIR)/% : $(OBJDIR)/%.o
	$(FC) $(PROFLFLAGS) $< $(OBJS) $(LFLAGS) -o $@

.PHONY : clean
clean :
	rm $(OBJDIR)/* $(BINDIR)/*
