FC = gfortran
OBJDIR = ../mod
#LFLAGS = -warn all -O -pg
#CFLAGS = -warn all -O -pg
LFLAGS = -pg
#CFLAGS = -Wall -O -pg
CFLAGS = -Wall -O2
OBJS =  $(OBJDIR)/utility.o $(OBJDIR)/shift.o $(OBJDIR)/nullvec.o $(OBJDIR)/prec.o $(OBJDIR)/rotation.o $(OBJDIR)/upper_decomp.o $(OBJDIR)/triangular.o $(OBJDIR)/gs.o $(OBJDIR)/assemble.o $(OBJDIR)/band.o $(OBJDIR)/upper_convert.o

run_test : test_decomp test_convert
	./test_decomp
	./test_convert

test_decomp : $(OBJDIR)/test_decomp.o $(OBJS)
	$(FC) $(LFLAGS) $(OBJDIR)/test_decomp.o $(OBJS) -o test_decomp

test_convert : $(OBJDIR)/test_convert.o $(OBJS)
	$(FC) $(LFLAGS) $(OBJDIR)/test_convert.o $(OBJS) -o test_convert

$(OBJDIR)/test_decomp.o : $(OBJDIR)/decomp.o

$(OBJDIR)/decomp.o : force_look
	cd ../src; $(MAKE)

$(OBJDIR)/%.o : %.f90
	$(FC) $(CFLAGS) -I$(OBJDIR) -J$(OBJDIR) -c $< -o $@


force_look :
	true

clean :
	rm $(OBJDIR)/*.o $(OBJDIR)/*.mod test_decomp

prof : run_test
	gprof -c -q test_decomp gmon.out | less

kill :
	kill -9 `ps -e | grep test_nest | sed -e "s/ pts.*//"`

ps :
	ps -e | grep test_nest

