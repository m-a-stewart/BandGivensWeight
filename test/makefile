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
	$(BINDIR)/test_qr_factorization
	$(BINDIR)/test_solve

all : $(BINDIR)/test_general $(BINDIR)/test_general_bt \
	$(BINDIR)/test_general_wb $(BINDIR)/test_general_ubt $(BINDIR)/test_general_wbv \
	$(BINDIR)/test_convert_ub_and_bv $(BINDIR)/test_convert_wb_and_bt  \
	$(BINDIR)/test_convert_ubt_and_wbv $(BINDIR)/test_row_compress \
	$(BINDIR)/test_qr_factorization $(BINDIR)/test_solve

orth_rank : force_look
	cd ../src; $(MAKE)

force_look :
	true

run_row_compress : $(BINDIR)/test_row_compress
	$(BINDIR)/test_row_compress

run_qr_factorization : $(BINDIR)/test_qr_factorization
	$(BINDIR)/test_qr_factorization

$(BINDIR)/% : $(OBJDIR)/%.o
	$(FC) $(LFLAGS) $(PROFLFLAGS) $(OBJDIR)/$< $(OBJS) \
	$(OBJDIR)/test_data.o -o $@

$(OBJDIR)/test_data.o : orth_rank
$(OBJDIR)/test_general.o : $(OBJDIR)/test_data.o
$(OBJDIR)/test_general_bt.o : $(OBJDIR)/test_data.o
$(OBJDIR)/test_general_wb.o : $(OBJDIR)/test_data.o
$(OBJDIR)/test_general_ubt.o : $(OBJDIR)/test_data.o
$(OBJDIR)/test_general_wbv.o : $(OBJDIR)/test_data.o
$(OBJDIR)/test_convert_ub_and_bv.o : $(OBJDIR)/test_data.o
$(OBJDIR)/test_convert_wb_and_bt.o : $(OBJDIR)/test_data.o
$(OBJDIR)/test_convert_ubt_and_wbv.o : $(OBJDIR)/test_data.o
$(OBJDIR)/test_row_compress.o : $(OBJDIR)/test_data.o
$(OBJDIR)/test_qr_factorization.o : $(OBJDIR)/test_data.o
$(OBJDIR)/test_solve.o : $(OBJDIR)/test_data.o

$(OBJDIR)/%.o : %.f90
	$(FC) $(CFLAGS) $(PROFCFLAGS) -I$(OBJDIR) -J$(OBJDIR) -c $< -o $@

clean :
	rm $(OBJDIR)/*.o $(OBJDIR)/*.mod $(BINDIR)/*

prof : test_general
	gprof -c -q --no-time='__general_ub_MOD_f_d_upper_to_ub' \
		-no-time='__assemble_MOD_d_bv_to_upper' test_general | less
