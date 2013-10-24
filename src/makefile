FC = gfortran
OBJDIR = ../mod
#LFLAGS = -warn all -O -pg
#CFLAGS = -warn all -O -pg
LFLAGS = 
#CFLAGS = -Wall -fbounds-check
CFLAGS = -Wall -O2
#CFLAGS = -Wall -O -pg

$(OBJDIR)/%.o : %.f90
	$(FC) $(CFLAGS) -I$(OBJDIR) -J$(OBJDIR) -c $< -o $@

$(OBJDIR)/upper_decomp.o : $(OBJDIR)/prec.o $(OBJDIR)/nullvec.o $(OBJDIR)/rotation.o $(OBJDIR)/shift.o
$(OBJDIR)/nullvec.o : $(OBJDIR)/utility.o $(OBJDIR)/triangular.o $(OBJDIR)/prec.o
$(OBJDIR)/utility.o : $(OBJDIR)/prec.o
$(OBJDIR)/triangular.o : $(OBJDIR)/prec.o
$(OBJDIR)/shift.o : $(OBJDIR)/prec.o
$(OBJDIR)/rotation.o : $(OBJDIR)/prec.o

clean :
	rm $(OBJDIR)/*.o $(OBJDIR)/*.mod

kill :
	kill -9 `ps -e | grep test_nest | sed -e "s/ pts.*//"`

ps :
	ps -e | grep test_nest

