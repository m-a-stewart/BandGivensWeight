include ../include/vars.mk

run_all : all
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

all : NESTED $(OBJDIR)/test_recompress $(OBJDIR)/test_general $(OBJDIR)/test_general_bt \
		$(OBJDIR)/test_general_wb $(OBJDIR)/test_general_ubt $(OBJDIR)/test_general_wbv \
		$(OBJDIR)/test_convert_ub_and_bv $(OBJDIR)/test_convert_wb_and_bt  \
		$(OBJDIR)/test_convert_ubt_and_wbv $(OBJDIR)/test_row_compress \
		$(OBJDIR)/test_qr_iteration $(OBJDIR)/test_qr_factorization $(OBJDIR)/test_solve \
		$(OBJDIR)/test_update

NESTED :
	cd ../src; make

run_row_compress : $(OBJDIR)/test_row_compress
	$(OBJDIR)/test_row_compress

run_qr_factorization : $(OBJDIR)/test_qr_factorization
	$(OBJDIR)/test_qr_factorization

$(OBJDIR)/test_general : $(OBJDIR)/test_general.o $(OBJS) $(OBJDIR)/test_data.o
	$(FC) $(LFLAGS) $(PROFLFLAGS) $(OBJDIR)/test_general.o $(OBJS) \
	$(OBJDIR)/test_data.o -o $(OBJDIR)/test_general

$(OBJDIR)/test_general_bt : $(OBJDIR)/test_general_bt.o $(OBJS) $(OBJDIR)/test_data.o
	$(FC) $(LFLAGS) $(PROFLFLAGS) $(OBJDIR)/test_general_bt.o $(OBJS) \
	$(OBJDIR)/test_data.o -o $(OBJDIR)/test_general_bt

$(OBJDIR)/test_general_wb : $(OBJDIR)/test_general_wb.o $(OBJS) $(OBJDIR)/test_data.o
	$(FC) $(LFLAGS) $(PROFLFLAGS) $(OBJDIR)/test_general_wb.o $(OBJS) \
	$(OBJDIR)/test_data.o -o $(OBJDIR)/test_general_wb

$(OBJDIR)/test_general_ubt : $(OBJDIR)/test_general_ubt.o $(OBJS) $(OBJDIR)/test_data.o
	$(FC) $(LFLAGS) $(PROFLFLAGS) $(OBJDIR)/test_general_ubt.o $(OBJS) \
	$(OBJDIR)/test_data.o -o $(OBJDIR)/test_general_ubt

$(OBJDIR)/test_general_wbv : $(OBJDIR)/test_general_wbv.o $(OBJS) $(OBJDIR)/test_data.o
	$(FC) $(LFLAGS) $(PROFLFLAGS) $(OBJDIR)/test_general_wbv.o $(OBJS) \
	$(OBJDIR)/test_data.o -o $(OBJDIR)/test_general_wbv

$(OBJDIR)/test_convert_ub_and_bv : $(OBJDIR)/test_convert_ub_and_bv.o $(OBJS) $(OBJDIR)/test_data.o
	$(FC) $(LFLAGS) $(PROFLFLAGS) $(OBJDIR)/test_convert_ub_and_bv.o $(OBJS) \
	$(OBJDIR)/test_data.o -o $(OBJDIR)/test_convert_ub_and_bv

$(OBJDIR)/test_convert_wb_and_bt : $(OBJDIR)/test_convert_wb_and_bt.o $(OBJS) $(OBJDIR)/test_data.o
	$(FC) $(LFLAGS) $(PROFLFLAGS) $(OBJDIR)/test_convert_wb_and_bt.o $(OBJS) \
	$(OBJDIR)/test_data.o -o $(OBJDIR)/test_convert_wb_and_bt

$(OBJDIR)/test_convert_ubt_and_wbv : $(OBJDIR)/test_convert_ubt_and_wbv.o $(OBJS) $(OBJDIR)/test_data.o
	$(FC) $(LFLAGS) $(PROFLFLAGS) $(OBJDIR)/test_convert_ubt_and_wbv.o $(OBJS) \
	$(OBJDIR)/test_data.o -o $(OBJDIR)/test_convert_ubt_and_wbv

$(OBJDIR)/test_row_compress : $(OBJDIR)/test_row_compress.o $(OBJS) $(OBJDIR)/test_data.o
	$(FC) $(LFLAGS) $(PROFLFLAGS) $(OBJDIR)/test_row_compress.o $(OBJS) \
	$(OBJDIR)/test_data.o -o $(OBJDIR)/test_row_compress

$(OBJDIR)/test_recompress : $(OBJDIR)/test_recompress.o $(OBJS) $(OBJDIR)/test_data.o
	$(FC) $(LFLAGS) $(PROFLFLAGS) $(OBJDIR)/test_recompress.o $(OBJS) \
	$(OBJDIR)/test_data.o -o $(OBJDIR)/test_recompress

$(OBJDIR)/test_qr_iteration : $(OBJDIR)/test_qr_iteration.o $(OBJS) $(OBJDIR)/test_data.o
	$(FC) $(LFLAGS) $(PROFLFLAGS) $(OBJDIR)/test_qr_iteration.o $(OBJS) \
	$(OBJDIR)/test_data.o -o $(OBJDIR)/test_qr_iteration

$(OBJDIR)/test_qr_factorization : $(OBJDIR)/test_qr_factorization.o $(OBJS) $(OBJDIR)/test_data.o
	$(FC) $(LFLAGS) $(PROFLFLAGS) $(OBJDIR)/test_qr_factorization.o $(OBJS) \
	$(OBJDIR)/test_data.o -o $(OBJDIR)/test_qr_factorization

$(OBJDIR)/test_solve : $(OBJDIR)/test_solve.o $(OBJS) $(OBJDIR)/test_data.o
	$(FC) $(LFLAGS) $(PROFLFLAGS) $(OBJDIR)/test_solve.o $(OBJS) \
	$(OBJDIR)/test_data.o -o $(OBJDIR)/test_solve

$(OBJDIR)/test_sweeps1 : $(OBJDIR)/test_sweeps1.o $(OBJS) $(OBJDIR)/test_data.o
	$(FC) $(LFLAGS) $(PROFLFLAGS) $(OBJDIR)/test_sweeps1.o $(OBJS) \
	$(OBJDIR)/test_data.o -o $(OBJDIR)/test_sweeps1

$(OBJDIR)/test_update : $(OBJDIR)/test_update.o $(OBJS) $(OBJDIR)/test_data.o
	$(FC) $(LFLAGS) $(PROFLFLAGS) $(OBJDIR)/test_update.o $(OBJS) \
	$(OBJDIR)/test_data.o -o $(OBJDIR)/test_update

$(OBJDIR)/test_data.o : test_data.f90

$(OBJDIR)/%.o : %.f90 $(OBJDIR)/test_data.o
	$(FC) $(CFLAGS) $(PROFCFLAGS) -I$(OBJDIR) -J$(OBJDIR) -c $< -o $@

clean :
	rm $(OBJDIR)/*.o $(OBJDIR)/*.mod $(OBJDIR)/$(TEST_EXECS)

prof : test_recompress
	gprof -c -q --no-time='__general_ub_MOD_f_d_upper_to_ub' \
		-no-time='__assemble_MOD_d_bv_to_upper' test_recompress | less

kill :
	kill -9 `ps -e | grep test_nest | sed -e "s/ pts.*//"`

ps :
	ps -e | grep test_nest

