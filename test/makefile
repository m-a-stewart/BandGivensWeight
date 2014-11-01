include ../include/vars.mk

run_all : $(OBJDIR)/test_recompress $(OBJDIR)/test_general $(OBJDIR)/test_general_bt \
		$(OBJDIR)/test_general_wb $(OBJDIR)/test_general_ubt $(OBJDIR)/test_general_wbv \
		$(OBJDIR)/test_convert_ub_and_bv $(OBJDIR)/test_convert_wb_and_bt  \
		$(OBJDIR)/test_convert_ubt_and_wbv $(OBJDIR)/test_row_compress \
		$(OBJDIR)/test_qr_iteration $(OBJDIR)/test_qr_factorization $(OBJDIR)/test_solve \
		$(OBJDIR)/test_update
	$(OBJDIR)/test_general
	$(OBJDIR)/test_general_bt
	$(OBJDIR)/test_general_wb
	$(OBJDIR)/test_general_ubt
	$(OBJDIR)/test_general_wbv
	$(OBJDIR)/test_convert_ub_and_bv
	$(OBJDIR)/test_convert_wb_and_bt
	$(OBJDIR)/test_convert_ubt_and_wbv
	$(OBJDIR)/test_row_compress
	$(OBJDIR)/test_recompress
	$(OBJDIR)/test_qr_factorization
	$(OBJDIR)/test_qr_iteration
	$(OBJDIR)/test_solve
	$(OBJDIR)/test_update

run_row_compress : $(OBJDIR)/test_row_compress
	$(OBJDIR)/test_row_compress

run_qr_factorization : NESTED $(BINDIR)/test_qr_factorization
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

kill :
	kill -9 `ps -e | grep test_nest | sed -e "s/ pts.*//"`

ps :
	ps -e | grep test_nest

