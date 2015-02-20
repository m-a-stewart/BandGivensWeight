include ../include/vars.mk
all : $(OBJS)

$(OBJDIR)/orth_rank.o : $(OBJDIR)/general.o $(OBJDIR)/assemble.o \
	$(OBJDIR)/qr_factorization.o $(OBJDIR)/solve.o \
	$(OBJDIR)/row_compress.o $(OBJDIR)/transforms.o \
	$(OBJDIR)/convert.o $(OBJDIR)/types.o $(OBJDIR)/misc.o


$(OBJDIR)/misc.o : $(OBJDIR)/error_id.o $(OBJDIR)/utility.o
$(OBJDIR)/error_id.o : $(OBJDIR)/prec.o $(OBJDIR)/shift.o
$(OBJDIR)/utility.o : $(OBJDIR)/prec.o


$(OBJDIR)/types.o : $(OBJDIR)/band_types.o $(OBJDIR)/orth_band_types.o \
	$(OBJDIR)/random.o
$(OBJDIR)/band_types.o : $(OBJDIR)/rotation.o $(OBJDIR)/utility.o $(OBJDIR)/shift.o \
	$(OBJDIR)/prec.o
$(OBJDIR)/orth_band_types.o : $(OBJDIR)/prec.o $(OBJDIR)/utility.o
$(OBJDIR)/random.o : $(OBJDIR)/prec.o $(OBJDIR)/utility.o \
	$(OBJDIR)/band_types.o $(OBJDIR)/orth_band_types.o

$(OBJDIR)/transforms.o : $(OBJDIR)/sweeps.o $(OBJDIR)/rotation.o $(OBJDIR)/shift.o
$(OBJDIR)/rotation.o : $(OBJDIR)/prec.o
$(OBJDIR)/sweeps.o : $(OBJDIR)/band_types.o $(OBJDIR)/prec.o $(OBJDIR)/utility.o
$(OBJDIR)/shift.o : $(OBJDIR)/prec.o


$(OBJDIR)/orth.o : $(OBJDIR)/triangular.o $(OBJDIR)/gs.o $(OBJDIR)/nullvec.o
$(OBJDIR)/triangular.o : $(OBJDIR)/prec.o
$(OBJDIR)/nullvec.o : $(OBJDIR)/prec.o $(OBJDIR)/error_id.o $(OBJDIR)/utility.o \
	$(OBJDIR)/triangular.o
$(OBJDIR)/gs.o : $(OBJDIR)/utility.o $(OBJDIR)/prec.o $(OBJDIR)/error_id.o


$(OBJDIR)/assemble.o : $(OBJDIR)/prec.o $(OBJDIR)/error_id.o $(OBJDIR)/band_types.o \
	$(OBJDIR)/orth_band_types.o $(OBJDIR)/rotation.o

$(OBJDIR)/general.o : $(OBJDIR)/general_bv.o $(OBJDIR)/general_ub.o $(OBJDIR)/general_bt.o \
	$(OBJDIR)/general_wb.o $(OBJDIR)/general_ubt.o $(OBJDIR)/general_wbv.o
$(OBJDIR)/general_bv.o : $(OBJDIR)/orth.o $(OBJDIR)/shift.o $(OBJDIR)/rotation.o \
	$(OBJDIR)/prec.o $(OBJDIR)/error_id.o $(OBJDIR)/utility.o \
	$(OBJDIR)/orth_band_types.o $(OBJDIR)/band_types.o
$(OBJDIR)/general_ub.o : $(OBJDIR)/orth.o $(OBJDIR)/shift.o $(OBJDIR)/rotation.o \
	$(OBJDIR)/prec.o $(OBJDIR)/error_id.o $(OBJDIR)/utility.o \
	$(OBJDIR)/orth_band_types.o $(OBJDIR)/band_types.o
$(OBJDIR)/general_bt.o : $(OBJDIR)/prec.o $(OBJDIR)/error_id.o $(OBJDIR)/utility.o \
	$(OBJDIR)/orth_band_types.o $(OBJDIR)/band_types.o $(OBJDIR)/general_ub.o
$(OBJDIR)/general_wb.o : $(OBJDIR)/prec.o $(OBJDIR)/error_id.o $(OBJDIR)/utility.o \
	$(OBJDIR)/orth_band_types.o $(OBJDIR)/band_types.o $(OBJDIR)/general_bv.o
$(OBJDIR)/general_ubt.o : $(OBJDIR)/prec.o $(OBJDIR)/error_id.o $(OBJDIR)/utility.o \
	$(OBJDIR)/orth_band_types.o $(OBJDIR)/band_types.o $(OBJDIR)/general_ub.o
$(OBJDIR)/general_wbv.o : $(OBJDIR)/prec.o $(OBJDIR)/error_id.o $(OBJDIR)/utility.o \
	$(OBJDIR)/orth_band_types.o $(OBJDIR)/band_types.o $(OBJDIR)/general_bv.o


$(OBJDIR)/convert.o : $(OBJDIR)/convert_bv_to_ub.o $(OBJDIR)/convert_ub_to_bv.o \
	$(OBJDIR)/convert_wb_to_bt.o $(OBJDIR)/convert_bt_to_wb.o \
	$(OBJDIR)/convert_wbv_to_ubt.o $(OBJDIR)/convert_ubt_to_wbv.o 
$(OBJDIR)/convert_bv_to_ub.o : $(OBJDIR)/prec.o $(OBJDIR)/error_id.o \
	$(OBJDIR)/rotation.o $(OBJDIR)/orth_band_types.o $(OBJDIR)/band_types.o
$(OBJDIR)/convert_ub_to_bv.o : $(OBJDIR)/prec.o $(OBJDIR)/error_id.o $(OBJDIR)/shift.o \
	$(OBJDIR)/rotation.o $(OBJDIR)/orth_band_types.o $(OBJDIR)/band_types.o
$(OBJDIR)/convert_wb_to_bt.o : $(OBJDIR)/prec.o $(OBJDIR)/error_id.o \
	$(OBJDIR)/rotation.o $(OBJDIR)/orth_band_types.o $(OBJDIR)/band_types.o
$(OBJDIR)/convert_bt_to_wb.o : $(OBJDIR)/prec.o $(OBJDIR)/error_id.o $(OBJDIR)/shift.o \
	$(OBJDIR)/rotation.o $(OBJDIR)/orth_band_types.o $(OBJDIR)/band_types.o
$(OBJDIR)/convert_wbv_to_ubt.o : $(OBJDIR)/prec.o $(OBJDIR)/error_id.o $(OBJDIR)/shift.o \
	$(OBJDIR)/rotation.o $(OBJDIR)/orth_band_types.o $(OBJDIR)/band_types.o
$(OBJDIR)/convert_ubt_to_wbv.o : $(OBJDIR)/prec.o $(OBJDIR)/error_id.o $(OBJDIR)/shift.o \
	$(OBJDIR)/rotation.o $(OBJDIR)/orth_band_types.o $(OBJDIR)/band_types.o


$(OBJDIR)/row_compress.o : $(OBJDIR)/misc.o $(OBJDIR)/sweeps.o 

$(OBJDIR)/qr_factorization.o : $(OBJDIR)/prec.o $(OBJDIR)/error_id.o \
	$(OBJDIR)/convert_bv_to_ub.o $(OBJDIR)/convert_ub_to_bv.o \
	$(OBJDIR)/sweeps.o $(OBJDIR)/shift.o $(OBJDIR)/rotation.o \
	$(OBJDIR)/band_types.o $(OBJDIR)/orth_band_types.o

$(OBJDIR)/solve.o : $(OBJDIR)/prec.o $(OBJDIR)/error_id.o $(OBJDIR)/orth_band_types.o \
	$(OBJDIR)/band_types.o $(OBJDIR)/rotation.o

$(OBJDIR)/%.o : %.f90
	$(FC) $(CFLAGS) $(PROFCFLAGS) -I$(OBJDIR) -J$(OBJDIR) -c $< -o $@

$(OBJDIR)/%.o : transforms/%.f90
	$(FC) $(CFLAGS) $(PROFCFLAGS) -I$(OBJDIR) -J$(OBJDIR) -c $< -o $@

$(OBJDIR)/%.o : orth/%.f90
	$(FC) $(CFLAGS) $(PROFCFLAGS) -I$(OBJDIR) -J$(OBJDIR) -c $< -o $@

$(OBJDIR)/%.o : misc/%.f90
	$(FC) $(CFLAGS) $(PROFCFLAGS) -I$(OBJDIR) -J$(OBJDIR) -c $< -o $@

$(OBJDIR)/%.o : general/%.f90
	$(FC) $(CFLAGS) $(PROFCFLAGS) -I$(OBJDIR) -J$(OBJDIR) -c $< -o $@

$(OBJDIR)/%.o : types/%.f90
	$(FC) $(CFLAGS) $(PROFCFLAGS) -I$(OBJDIR) -J$(OBJDIR) -c $< -o $@

$(OBJDIR)/%.o : convert/%.f90
	$(FC) $(CFLAGS) $(PROFCFLAGS) -I$(OBJDIR) -J$(OBJDIR) -c $< -o $@

$(OBJDIR)/%.o : qr_factorization/%.f90
	$(FC) $(CFLAGS) $(PROFCFLAGS) -I$(OBJDIR) -J$(OBJDIR) -c $< -o $@

$(OBJDIR)/%.o : solve/%.f90
	$(FC) $(CFLAGS) $(PROFCFLAGS) -I$(OBJDIR) -J$(OBJDIR) -c $< -o $@

$(OBJDIR)/%.o : assemble/%.f90
	$(FC) $(CFLAGS) $(PROFCFLAGS) -I$(OBJDIR) -J$(OBJDIR) -c $< -o $@

$(OBJDIR)/%.o : orth_rank/%.f90
	$(FC) $(CFLAGS) $(PROFCFLAGS) -I$(OBJDIR) -J$(OBJDIR) -c $< -o $@

$(OBJDIR)/%.o : row_compress/%.f90
	$(FC) $(CFLAGS) $(PROFCFLAGS) -I$(OBJDIR) -J$(OBJDIR) -c $< -o $@


clean :
	rm $(OBJDIR)/*.o $(OBJDIR)/*.mod $(TEST_EXECS)
