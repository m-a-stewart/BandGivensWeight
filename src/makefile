FC = gfortran
OBJDIR = ../mod
#LFLAGS = -warn all -O -pg
#CFLAGS = -warn all -O -pg
LFLAGS = 
#CFLAGS = -Wall -fbounds-check
#CFLAGS = -Wall -O2

$(OBJDIR)/%.o : %.f90
	$(FC) $(CFLAGS) -I$(OBJDIR) -J$(OBJDIR) -c $< -o $@

$(OBJDIR)/decomp.o : $(OBJDIR)/prec.o $(OBJDIR)/nullvec.o $(OBJDIR)/rotation.o $(OBJDIR)/shift.o
$(OBJDIR)/nullvec.o : $(OBJDIR)/utility.o $(OBJDIR)/prec.o
$(OBJDIR)/utility.o : $(OBJDIR)/prec.o
$(OBJDIR)/shift.o : $(OBJDIR)/prec.o
$(OBJDIR)/rotation.o : $(OBJDIR)/prec.o



clean :
	rm $(OBJDIR)/*.o $(OBJDIR)/*.mod

prof :
	gprof -c test_nest gmon.out | less

kill :
	kill -9 `ps -e | grep test_nest | sed -e "s/ pts.*//"`

ps :
	ps -e | grep test_nest

