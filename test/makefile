include ../include/vars.mk

run_test : test_compress #test_decomp test_convert 
#	./test_decomp
#	./test_convert
	./test_compress

test_decomp : $(OBJDIR)/test_decomp.o $(OBJS)
	$(FC) $(LFLAGS) $(PROFLFLAGS) $(OBJDIR)/test_decomp.o $(OBJS) -o test_decomp

test_convert : $(OBJDIR)/test_convert.o $(OBJS)
	$(FC) $(LFLAGS) $(PROFLFLAGS) $(OBJDIR)/test_convert.o $(OBJS) -o test_convert

test_compress : $(OBJDIR)/test_compress.o $(OBJS)
	$(FC) $(LFLAGS) $(PROFLFLAGS) $(OBJDIR)/test_compress.o $(OBJS) -o test_compress


$(OBJDIR)/test_decomp.o : $(OBJDIR)/prec.o
$(OBJDIR)/test_compress.o : $(OBJDIR)/prec.o
$(OBJDIR)/test_convert.o : $(OBJDIR)/prec.o

$(OBJDIR)/prec.o : force_look
	cd ../src; $(MAKE)

$(OBJDIR)/%.o : %.f90
	$(FC) $(CFLAGS) $(PROFCFLAGS) -I$(OBJDIR) -J$(OBJDIR) -c $< -o $@


force_look :
	true

clean :
	rm $(OBJDIR)/*.o $(OBJDIR)/*.mod test_decomp test_convert test_compress

prof : run_test
	gprof -c -q test_decomp gmon.out | less

kill :
	kill -9 `ps -e | grep test_nest | sed -e "s/ pts.*//"`

ps :
	ps -e | grep test_nest

