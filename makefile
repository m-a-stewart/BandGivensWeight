include ./include/vars.mk
vpath %.f90 $(SRCDIRS) test

.PHONY : all
all : $(OBJDIR)/orth_rank.o

include $(SRCDIR)/src.mk

include $(TESTDIR)/test.mk

$(OBJDIR)/%.o : %.f90
	$(FC) $(CFLAGS) $(PROFCFLAGS) -I$(OBJDIR) -J$(OBJDIR) -c $< -o $@

$(BINDIR)/% : $(OBJDIR)/%.o
	$(FC) $(LFLAGS) $< $(OBJS) $(OBJDIR)/test_data.o -o $@

.PHONY : clean
clean :
	rm $(OBJDIR)/* $(BINDIR)/*

# prof : test_general
# 	gprof -c -q --no-time='__general_ub_MOD_f_d_upper_to_ub' \
# 		-no-time='__assemble_MOD_d_bv_to_upper' test_general | less

