_TESTS = test_general_ub test_general_bv test_general_bt test_general_wb \
	test_general_ubt test_general_wbv \
	test_convert_ub_and_bv test_convert_wb_and_bt \
	test_convert_ubt_and_wbv test_row_compress \
	test_qr_factorization test_solve test_cond_triangular test_products \
	test_submatrix test_cond_orth_band
TESTS = $(patsubst %,$(BINDIR)/%,$(_TESTS))

.PHONY : all_tests
all_tests : $(TESTS)

.PHONY : run_all_tests
run_all_tests : $(TESTS)
	$(BINDIR)/test_cond_triangular
	$(BINDIR)/test_cond_orth_band
	$(BINDIR)/test_products
	$(BINDIR)/test_submatrix
	$(BINDIR)/test_general_ub
	$(BINDIR)/test_general_bv
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

.PHONY : run_cond_triangular
run_cond_triangular : $(BINDIR)/test_cond_triangular
	$(BINDIR)/test_cond_triangular

.PHONY : run_cond_orth_band
run_cond_orth_band : $(BINDIR)/test_cond_orth_band
	$(BINDIR)/test_cond_orth_band

.PHONY : run_products
run_products : $(BINDIR)/test_products
	$(BINDIR)/test_products

.PHONY : run_submatrix
run_submatrix : $(BINDIR)/test_submatrix
	$(BINDIR)/test_submatrix

.PHONY : run_general
run_general : $(BINDIR)/test_general_ub $(BINDIR)/test_general_bv $(BINDIR)/test_general_bt \
	$(BINDIR)/test_general_wb $(BINDIR)/test_general_wbv $(BINDIR)/test_general_ubt
	$(BINDIR)/test_general_ub
	$(BINDIR)/test_general_bv
	$(BINDIR)/test_general_bt
	$(BINDIR)/test_general_wb
	$(BINDIR)/test_general_ubt
	$(BINDIR)/test_general_wbv

.PHONY : run_general_ub
run_general_ub : $(BINDIR)/test_general_ub
	$(BINDIR)/test_general_ub

.PHONY : run_general_bv
run_general_bv : $(BINDIR)/test_general_bv
	$(BINDIR)/test_general_bv

.PHONY : run_general_bt
run_general_bt : $(BINDIR)/test_general_bt
	$(BINDIR)/test_general_bt

.PHONY : run_general_wb
run_general_wb : $(BINDIR)/test_general_wb
	$(BINDIR)/test_general_wb

.PHONY : run_general_ubt
run_general_ubt : $(BINDIR)/test_general_ubt
	$(BINDIR)/test_general_ubt

.PHONY : run_general_wbv
run_general_wbv : $(BINDIR)/test_general_wbv
	$(BINDIR)/test_general_wbv

.PHONY : run_convert_ub_and_bv
run_convert_ub_and_bv : $(BINDIR)/test_convert_ub_and_bv
	$(BINDIR)/test_convert_ub_and_bv

.PHONY : run_convert_wb_and_bt
run_convert_wb_and_bt : $(BINDIR)/test_convert_wb_and_bt
	$(BINDIR)/test_convert_wb_and_bt

.PHONY : run_convert_ubt_and_wbv
run_convert_ubt_and_wbv : $(BINDIR)/test_convert_ubt_and_wbv
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

$(OBJDIR)/test_cond_triangular.o : $(OBJDIR)/orrb.o
$(OBJDIR)/test_cond_orth_band.o : $(OBJDIR)/orrb.o
$(OBJDIR)/test_products.o : $(OBJDIR)/orrb.o
$(OBJDIR)/test_submatrix.o : $(OBJDIR)/orrb.o
$(OBJDIR)/test_general_ub.o : $(OBJDIR)/orrb.o
$(OBJDIR)/test_general_bv.o : $(OBJDIR)/orrb.o
$(OBJDIR)/test_general_bt.o : $(OBJDIR)/orrb.o
$(OBJDIR)/test_general_wb.o : $(OBJDIR)/orrb.o
$(OBJDIR)/test_general_ubt.o : $(OBJDIR)/orrb.o
$(OBJDIR)/test_general_wbv.o : $(OBJDIR)/orrb.o
$(OBJDIR)/test_convert_ub_and_bv.o : $(OBJDIR)/orrb.o
$(OBJDIR)/test_convert_wb_and_bt.o : $(OBJDIR)/orrb.o
$(OBJDIR)/test_convert_ubt_and_wbv.o : $(OBJDIR)/orrb.o
$(OBJDIR)/test_row_compress.o : $(OBJDIR)/orrb.o
$(OBJDIR)/test_qr_factorization.o : $(OBJDIR)/orrb.o
$(OBJDIR)/test_solve.o : $(OBJDIR)/orrb.o
