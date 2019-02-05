PROG = m2dpt

SRCS = $(wildcard *.f90) 
OBJS = $(SRCS:.f90=.o)    
    
F90 = mpif90
OPT = -O3
F90FLAGS = $(OPT) 
LDFLAGS = $(OPT) 

all: $(PROG)

m2dpt: $(OBJS)
	$(F90) $(LDFLAGS) -o $@ $^ $(LIBS)

clean:
	rm -rf $(PROG) *.o *.mod

.SUFFIXES: $(SUFFIXES) .f90

.f90.o:
	$(F90) $(F90FLAGS) -c $<

.mod.o:
	$(F90) $(F90FLAGS) -c $*.f90

m2dpt.o: physics.o numerics.o
