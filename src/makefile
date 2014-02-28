include ../include/vars.mk
all : $(OBJS)

$(OBJDIR)/misc.o : $(OBJDIR)/error_id.o $(OBJDIR)/utility.o
$(OBJDIR)/error_id.o : $(OBJDIR)/prec.o
$(OBJDIR)/utility.o : $(OBJDIR)/prec.o

$(OBJDIR)/types.o : $(OBJDIR)/band_types.o $(OBJDIR)/nested_types.o
$(OBJDIR)/band_types.o : $(OBJDIR)/rotation.o $(OBJDIR)/utility.o
$(OBJDIR)/nested_types.o : $(OBJDIR)/prec.o

$(OBJDIR)/transforms.o : $(OBJDIR)/sweeps.o $(OBJDIR)/rotation.o $(OBJDIR)/shift.o
$(OBJDIR)/rotation.o : $(OBJDIR)/prec.o
$(OBJDIR)/sweeps.o : $(OBJDIR)/rotation.o
$(OBJDIR)/shift.o : $(OBJDIR)/prec.o

$(OBJDIR)/orth.o : $(OBJDIR)/triangular.o $(OBJDIR)/gs.o $(OBJDIR)/nullvec.o
$(OBJDIR)/triangular.o : $(OBJDIR)/prec.o
$(OBJDIR)/nullvec.o : $(OBJDIR)/misc.o
$(OBJDIR)/gs.o : $(OBJDIR)/misc.o

$(OBJDIR)/assemble.o : $(OBJDIR)/band_types.o $(OBJDIR)/nested_types.o $(OBJDIR)/misc.o

$(OBJDIR)/decomp.o : $(OBJDIR)/general_bv.o $(OBJDIR)/general_ub.o
$(OBJDIR)/general_bv.o : $(OBJDIR)/orth.o $(OBJDIR)/transforms.o $(OBJDIR)/types.o
$(OBJDIR)/general_ub.o : $(OBJDIR)/orth.o $(OBJDIR)/transforms.o $(OBJDIR)/types.o

$(OBJDIR)/conversions.o : $(OBJDIR)/conversions_bv_to_ub.o $(OBJDIR)/conversions_ub_to_bv.o
$(OBJDIR)/conversions_bv_to_ub.o : $(OBJDIR)/misc.o $(OBJDIR)/transforms.o $(OBJDIR)/types.o
$(OBJDIR)/conversions_ub_to_bv.o : $(OBJDIR)/misc.o $(OBJDIR)/transforms.o $(OBJDIR)/types.o

$(OBJDIR)/compressions.o : $(OBJDIR)/compressions_bv_to_ub.o $(OBJDIR)/compressions_ub_to_bv.o
$(OBJDIR)/compressions_bv_to_ub.o : $(OBJDIR)/misc.o $(OBJDIR)/transforms.o $(OBJDIR)/types.o \
	$(OBJDIR)/orth.o
$(OBJDIR)/compressions_ub_to_bv.o : $(OBJDIR)/misc.o $(OBJDIR)/transforms.o $(OBJDIR)/types.o \
	$(OBJDIR)/orth.o

$(OBJDIR)/qr_iteration.o : $(OBJDIR)/compressions.o $(OBJDIR)/conversions.o

$(OBJDIR)/qr_factorization.o : $(OBJDIR)/conversions.o

$(OBJDIR)/back_substitution.o : $(OBJDIR)/transforms.o $(OBJDIR)/misc.o $(OBJDIR)/types.o

$(OBJDIR)/nested.o : $(OBJDIR)/decomp.o $(OBJDIR)/assemble.o $(OBJDIR)/qr_iteration.o \
	$(OBJDIR)/qr_factorization.o $(OBJDIR)/back_substitution.o

$(OBJDIR)/%.o : %.f90
	$(FC) $(CFLAGS) $(PROFCFLAGS) -I$(OBJDIR) -J$(OBJDIR) -c $< -o $@

$(OBJDIR)/%.o : transforms/%.f90
	$(FC) $(CFLAGS) $(PROFCFLAGS) -I$(OBJDIR) -J$(OBJDIR) -c $< -o $@

$(OBJDIR)/%.o : orth/%.f90
	$(FC) $(CFLAGS) $(PROFCFLAGS) -I$(OBJDIR) -J$(OBJDIR) -c $< -o $@

$(OBJDIR)/%.o : misc/%.f90
	$(FC) $(CFLAGS) $(PROFCFLAGS) -I$(OBJDIR) -J$(OBJDIR) -c $< -o $@

$(OBJDIR)/%.o : decomp/%.f90
	$(FC) $(CFLAGS) $(PROFCFLAGS) -I$(OBJDIR) -J$(OBJDIR) -c $< -o $@

$(OBJDIR)/%.o : types/%.f90
	$(FC) $(CFLAGS) $(PROFCFLAGS) -I$(OBJDIR) -J$(OBJDIR) -c $< -o $@

$(OBJDIR)/%.o : convert/%.f90
	$(FC) $(CFLAGS) $(PROFCFLAGS) -I$(OBJDIR) -J$(OBJDIR) -c $< -o $@

$(OBJDIR)/%.o : compress/%.f90
	$(FC) $(CFLAGS) $(PROFCFLAGS) -I$(OBJDIR) -J$(OBJDIR) -c $< -o $@

$(OBJDIR)/%.o : qr_iteration/%.f90
	$(FC) $(CFLAGS) $(PROFCFLAGS) -I$(OBJDIR) -J$(OBJDIR) -c $< -o $@

$(OBJDIR)/%.o : qr_factorization/%.f90
	$(FC) $(CFLAGS) $(PROFCFLAGS) -I$(OBJDIR) -J$(OBJDIR) -c $< -o $@

$(OBJDIR)/%.o : solve/%.f90
	$(FC) $(CFLAGS) $(PROFCFLAGS) -I$(OBJDIR) -J$(OBJDIR) -c $< -o $@

$(OBJDIR)/%.o : assemble/%.f90
	$(FC) $(CFLAGS) $(PROFCFLAGS) -I$(OBJDIR) -J$(OBJDIR) -c $< -o $@

$(OBJDIR)/%.o : nested/%.f90
	$(FC) $(CFLAGS) $(PROFCFLAGS) -I$(OBJDIR) -J$(OBJDIR) -c $< -o $@

clean :
	rm $(OBJDIR)/*.o $(OBJDIR)/*.mod $(TEST_EXECS)

kill :
	kill -9 `ps -e | grep test_nest | sed -e "s/ pts.*//"`

ps :
	ps -e | grep test_nest

