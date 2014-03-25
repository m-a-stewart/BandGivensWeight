include ../include/vars.mk

run_test : $(OBJDIR)/test_compress $(OBJDIR)/test_decomp $(OBJDIR)/test_convert \
		$(OBJDIR)/test_qr_iteration $(OBJDIR)/test_qr_factorization $(OBJDIR)/test_solve \
		$(OBJDIR)/test_sweeps $(OBJDIR)/test_update
	$(OBJDIR)/test_decomp
	$(OBJDIR)/test_convert
	$(OBJDIR)/test_compress
	$(OBJDIR)/test_qr_iteration
	$(OBJDIR)/test_qr_factorization
	$(OBJDIR)/test_solve
	$(OBJDIR)/test_sweeps
	$(OBJDIR)/test_update


$(OBJDIR)/test_decomp : $(OBJDIR)/test_decomp.o $(OBJS) $(OBJDIR)/test_data.o
	$(FC) $(LFLAGS) $(PROFLFLAGS) $(OBJDIR)/test_decomp.o $(OBJS) \
	$(OBJDIR)/test_data.o -o $(OBJDIR)/test_decomp

$(OBJDIR)/test_convert : $(OBJDIR)/test_convert.o $(OBJS) $(OBJDIR)/test_data.o
	$(FC) $(LFLAGS) $(PROFLFLAGS) $(OBJDIR)/test_convert.o $(OBJS) \
	$(OBJDIR)/test_data.o -o $(OBJDIR)/test_convert

$(OBJDIR)/test_compress : $(OBJDIR)/test_compress.o $(OBJS) $(OBJDIR)/test_data.o
	$(FC) $(LFLAGS) $(PROFLFLAGS) $(OBJDIR)/test_compress.o $(OBJS) \
	$(OBJDIR)/test_data.o -o $(OBJDIR)/test_compress

$(OBJDIR)/test_qr_iteration : $(OBJDIR)/test_qr_iteration.o $(OBJS) $(OBJDIR)/test_data.o
	$(FC) $(LFLAGS) $(PROFLFLAGS) $(OBJDIR)/test_qr_iteration.o $(OBJS) \
	$(OBJDIR)/test_data.o -o $(OBJDIR)/test_qr_iteration

$(OBJDIR)/test_qr_factorization : $(OBJDIR)/test_qr_factorization.o $(OBJS) $(OBJDIR)/test_data.o
	$(FC) $(LFLAGS) $(PROFLFLAGS) $(OBJDIR)/test_qr_factorization.o $(OBJS) \
	$(OBJDIR)/test_data.o -o $(OBJDIR)/test_qr_factorization

$(OBJDIR)/test_solve : $(OBJDIR)/test_solve.o $(OBJS) $(OBJDIR)/test_data.o
	$(FC) $(LFLAGS) $(PROFLFLAGS) $(OBJDIR)/test_solve.o $(OBJS) \
	$(OBJDIR)/test_data.o -o $(OBJDIR)/test_solve

$(OBJDIR)/test_sweeps : $(OBJDIR)/test_sweeps.o $(OBJS) $(OBJDIR)/test_data.o
	$(FC) $(LFLAGS) $(PROFLFLAGS) $(OBJDIR)/test_sweeps.o $(OBJS) \
	$(OBJDIR)/test_data.o -o $(OBJDIR)/test_sweeps

$(OBJDIR)/test_update : $(OBJDIR)/test_update.o $(OBJS) $(OBJDIR)/test_data.o
	$(FC) $(LFLAGS) $(PROFLFLAGS) $(OBJDIR)/test_update.o $(OBJS) \
	$(OBJDIR)/test_data.o -o $(OBJDIR)/test_update


$(OBJDIR)/test_data.o : test_data.f90 $(OBJDIR)/prec.o $(OBJS)
$(OBJDIR)/test_decomp.o : test_decomp.f90 $(OBJDIR)/prec.o $(OBJDIR)/test_data.o $(OBJS)
$(OBJDIR)/test_compress.o : test_compress.f90 $(OBJDIR)/prec.o $(OBJDIR)/test_data.o $(OBJS)
$(OBJDIR)/test_convert.o : test_convert.f90 $(OBJDIR)/prec.o $(OBJDIR)/test_data.o $(OBJS)
$(OBJDIR)/test_qr_factorization.o : test_qr_factorization.f90 $(OBJDIR)/prec.o $(OBJDIR)/test_data.o $(OBJS)
$(OBJDIR)/test_qr_iteration.o : test_qr_iteration.f90 $(OBJDIR)/prec.o $(OBJDIR)/test_data.o $(OBJS)
$(OBJDIR)/test_solve.o : test_solve.f90 $(OBJDIR)/prec.o $(OBJDIR)/test_data.o $(OBJS)
$(OBJDIR)/test_sweeps.o : test_sweeps.f90 $(OBJDIR)/prec.o $(OBJDIR)/test_data.o $(OBJS)
$(OBJDIR)/test_update.o : test_update.f90 $(OBJDIR)/prec.o $(OBJDIR)/test_data.o $(OBJS)

$(OBJDIR)/prec.o : force_look
	cd ../src; $(MAKE)

$(OBJDIR)/%.o : %.f90
	$(FC) $(CFLAGS) $(PROFCFLAGS) -I$(OBJDIR) -J$(OBJDIR) -c $< -o $@


force_look :
	true

clean :
	rm $(OBJDIR)/*.o $(OBJDIR)/*.mod $(OBJDIR)/$(TEST_EXECS)

prof : test_compress
	gprof -c -q --no-time='__general_ub_MOD_f_d_upper_to_ub' \
		-no-time='__assemble_MOD_d_bv_to_upper' test_compress | less

kill :
	kill -9 `ps -e | grep test_nest | sed -e "s/ pts.*//"`

ps :
	ps -e | grep test_nest

