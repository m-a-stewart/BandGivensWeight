FC = gfortran
OBJS = decomp.o
#LFLAGS = -warn all -O -pg
#CFLAGS = -warn all -O -pg
LFLAGS =
CFLAGS = -Wall

#test_nest : decomp.o test_nest.o
#	$(FC) $(LFLAGS) -o test_nest $(OBJS)

all: decomp.o nullvec.o utility.o rotation.o shift.o

decomp.o : prec.o nullvec.o rotation.o shift.o
nullvec.o : utility.o prec.o
utility.o : prec.o
shift.o: prec.o

%.o : %.f90
	$(FC) $(CFLAGS) -c $<

clean :
	rm *.o *.mod

prof :
	gprof -c test_nest gmon.out | less

kill :
	kill -9 `ps -e | grep test_nest | sed -e "s/ pts.*//"`

ps :
	ps -e | grep test_nest

