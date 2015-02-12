include ../include/vars.mk

run_all : all
	$(BINDIR)/test_general
	$(BINDIR)/test_general_bt
	$(BINDIR)/test_general_wb
	$(BINDIR)/test_general_ubt
	$(BINDIR)/test_general_wbv
	$(BINDIR)/test_convert_ub_and_bv
	$(BINDIR)/test_convert_wb_and_bt
	$(BINDIR)/test_convert_ubt_and_wbv
	$(BINDIR)/test_row_compress
	$(BINDIR)/test_recompress
	$(BINDIR)/test_qr_factorization
	$(BINDIR)/test_qr_iteration
	$(BINDIR)/test_solve
	$(BINDIR)/test_update

all : ORTH_RANK $(BINDIR)/test_recompress $(BINDIR)/test_general $(BINDIR)/test_general_bt \
		$(BINDIR)/test_general_wb $(BINDIR)/test_general_ubt $(BINDIR)/test_general_wbv \
		$(BINDIR)/test_convert_ub_and_bv $(BINDIR)/test_convert_wb_and_bt  \
		$(BINDIR)/test_convert_ubt_and_wbv $(BINDIR)/test_row_compress \
		$(BINDIR)/test_qr_iteration $(BINDIR)/test_qr_factorization $(BINDIR)/test_solve \
		$(BINDIR)/test_update

ORTH_RANK :
	cd ../src; make

run_row_compress : ORTH_RANK $(BINDIR)/test_row_compress
	$(BINDIR)/test_row_compress

run_qr_factorization : ORTH_RANK $(BINDIR)/test_qr_factorization
	$(BINDIR)/test_qr_factorization

$(BINDIR)/% : $(OBJDIR)/%.o $(OBJDIR)/test_data.o
	$(FC) $(LFLAGS) $(PROFLFLAGS) $(OBJDIR)/$< $(OBJS) \
	$(OBJDIR)/test_data.o -o $@

$(OBJDIR)/test_data.o : test_data.f90

$(OBJDIR)/%.o : %.f90 $(OBJDIR)/test_data.o
	$(FC) $(CFLAGS) $(PROFCFLAGS) -I$(OBJDIR) -J$(OBJDIR) -c $< -o $@

clean :
	rm $(OBJDIR)/*.o $(OBJDIR)/*.mod $(BINDIR)/*

prof : test_recompress
	gprof -c -q --no-time='__general_ub_MOD_f_d_upper_to_ub' \
		-no-time='__assemble_MOD_d_bv_to_upper' test_recompress | less
