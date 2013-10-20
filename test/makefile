FC = gfortran
OBJDIR = ../mod
#LFLAGS = -warn all -O -pg
#CFLAGS = -warn all -O -pg
LFLAGS =
#CFLAGS = -Wall
CFLAGS = 
OBJS =  $(OBJDIR)/utility.o $(OBJDIR)/shift.o $(OBJDIR)/nullvec.o $(OBJDIR)/prec.o $(OBJDIR)/rotation.o $(OBJDIR)/decomp.o

run_test : test_decomp
	./test_decomp

test_decomp : $(OBJDIR)/test_decomp.o $(OBJS)
	$(FC) $(LFLAGS) $(OBJDIR)/test_decomp.o $(OBJS) -o test_decomp

$(OBJDIR)/test_decomp.o : $(OBJDIR)/decomp.o

$(OBJDIR)/decomp.o : force_look
	cd ../src; $(MAKE)

$(OBJDIR)/%.o : %.f90
	$(FC) $(CFLAGS) -I$(OBJDIR) -J$(OBJDIR) -c $< -o $@


force_look :
	true

clean :
	rm $(OBJDIR)/*.o $(OBJDIR)/*.mod test_decomp

prof :
	gprof -c test_nest gmon.out | less

kill :
	kill -9 `ps -e | grep test_nest | sed -e "s/ pts.*//"`

ps :
	ps -e | grep test_nest

