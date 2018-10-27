FC = gfortran
FCOPTS = -O3
LN = $(FC)
LAPACK =  -llapack -lblas

FORTRANPATH=/home/litman/Yair/gitlab/lapackinterface

# Object files
OBJS = routines.o

all: $(OBJS)
	f2py -c routines.f90 $(LAPACK) tete.f90 -m tete

routines.o: routines.f90
	gfortran -c $(LAPACK) $(FLAGS)$<


clean:
	rm -f *.o *.so *.mod


