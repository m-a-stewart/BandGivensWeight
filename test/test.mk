_TESTS = test_general test_general_bt test_general_wb \
	test_general_ubt test_general_wbv \
	test_convert_ub_and_bv test_convert_wb_and_bt \
	test_convert_ubt_and_wbv test_row_compress \
	test_qr_factorization test_solve
TESTS = $(patsubst %,$(BINDIR)/%,$(_TESTS))

.PHONY : all_tests
all_tests : $(TESTS)

.PHONY : run_all_tests
run_all_tests : $(TESTS)
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

.PHONY : run_test_general
run_test_general : $(BINDIR)/test_general
	$(BINDIR)/test_general

.PHONY : run_test_general_bt
run_test_general_bt : $(BINDIR)/test_general_bt
	$(BINDIR)/test_general_bt

.PHONY : run_test_general_wb
run_test_general_wb : $(BINDIR)/test_general_wb
	$(BINDIR)/test_general_wb

.PHONY : run_test_general_ubt
run_test_general_ubt : $(BINDIR)/test_general_ubt
	$(BINDIR)/test_general_ubt

.PHONY : run_test_general_wbv
run_test_general_wbv : $(BINDIR)/test_general_wbv
	$(BINDIR)/test_general_wbv

.PHONY : run_test_convert_ub_and_bv
run_test_convert_ub_and_bv : $(BINDIR)/test_convert_ub_and_bv
	$(BINDIR)/test_convert_ub_and_bv

.PHONY : run_test_convert_wb_and_bt
run_test_convert_wb_and_bt : $(BINDIR)/test_convert_wb_and_bt
	$(BINDIR)/test_convert_wb_and_bt

.PHONY : run_test_convert_ubt_and_wbv
run_test_convert_ubt_and_wbv : $(BINDIR)/test_convert_ubt_and_wbv
	$(BINDIR)/test_convert_ubt_and_wbv

.PHONY : run_row_compress
run_row_compress : $(BINDIR)/test_row_compress
	$(BINDIR)/test_row_compress

.PHONY : run_qr_factorization
run_qr_factorization : $(BINDIR)/test_qr_factorization
	$(BINDIR)/test_qr_factorization

.PHONY : run_solve
run_solve : $(BINDIR)/test_solve
	$(BINDIR)/test_solve

$(OBJDIR)/test_data.o : $(OBJDIR)/orth_rank.o
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