include ../include/vars.mk
all : $(OBJS)

$(OBJDIR)/general_ub.o : $(OBJDIR)/error_id.o $(OBJDIR)/nullvec.o $(OBJDIR)/rotation.o $(OBJDIR)/shift.o \
		           $(OBJDIR)/gs.o $(OBJDIR)/nested_types.o
$(OBJDIR)/general_bv.o : $(OBJDIR)/error_id.o $(OBJDIR)/nullvec.o $(OBJDIR)/rotation.o $(OBJDIR)/shift.o \
		           $(OBJDIR)/gs.o $(OBJDIR)/nested_types.o
$(OBJDIR)/nullvec.o : $(OBJDIR)/utility.o $(OBJDIR)/triangular.o $(OBJDIR)/error_id.o
$(OBJDIR)/utility.o : $(OBJDIR)/prec.o
$(OBJDIR)/triangular.o : $(OBJDIR)/prec.o
$(OBJDIR)/shift.o : $(OBJDIR)/prec.o
$(OBJDIR)/error_id.o : $(OBJDIR)/prec.o
$(OBJDIR)/rotation.o : $(OBJDIR)/prec.o
$(OBJDIR)/nested_types.o : $(OBJDIR)/prec.o
$(OBJDIR)/gs.o : $(OBJDIR)/prec.o $(OBJDIR)/utility.o
$(OBJDIR)/assemble.o : $(OBJDIR)/error_id.o $(OBJDIR)/rotation.o $(OBJDIR)/band_types.o $(OBJDIR)/nested_types.o
$(OBJDIR)/band_types.o : $(OBJDIR)/prec.o $(OBJDIR)/rotation.o
$(OBJDIR)/conversions_ub_to_bv.o : $(OBJDIR)/prec.o $(OBJDIR)/band_types.o $(OBJDIR)/shift.o $(OBJDIR)/nested_types.o \
				$(OBJDIR)/rotation.o
$(OBJDIR)/conversions_bv_to_ub.o : $(OBJDIR)/prec.o $(OBJDIR)/band_types.o $(OBJDIR)/shift.o $(OBJDIR)/nested_types.o \
				$(OBJDIR)/rotation.o
$(OBJDIR)/compressions_ub_to_bv.o : $(OBJDIR)/prec.o $(OBJDIR)/band_types.o $(OBJDIR)/shift.o $(OBJDIR)/nested_types.o \
			$(OBJDIR)/rotation.o
$(OBJDIR)/compressions_bv_to_ub.o : $(OBJDIR)/prec.o $(OBJDIR)/band_types.o $(OBJDIR)/shift.o $(OBJDIR)/nested_types.o \
			$(OBJDIR)/rotation.o
$(OBJDIR)/qr_iteration.o : $(OBJDIR)/prec.o $(OBJDIR)/band_types.o $(OBJDIR)/shift.o $(OBJDIR)/nested_types.o \
			$(OBJDIR)/rotation.o



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



$(OBJDIR)/%.o : assemble/%.f90
	$(FC) $(CFLAGS) $(PROFCFLAGS) -I$(OBJDIR) -J$(OBJDIR) -c $< -o $@

clean :
	rm $(OBJDIR)/*.o $(OBJDIR)/*.mod

kill :
	kill -9 `ps -e | grep test_nest | sed -e "s/ pts.*//"`

ps :
	ps -e | grep test_nest

