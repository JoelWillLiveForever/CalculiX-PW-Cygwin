
CFLAGS = -Wall -O3  -I ../../spooles.2.2 -DARCH="Linux" -DSPOOLES -DARPACK -DMATRIXSTORAGE -DNETWORKOUT
FFLAGS = -Wall -O3 -fallow-argument-mismatch

CC=cc
FC=gfortran

.c.o :
	$(CC) $(CFLAGS) -c $<
.f.o :
	$(FC) $(FFLAGS) -c $<

include Makefile.inc

SCCXMAIN = ccx_2.20.c

OCCXF = $(SCCXF:.f=.o)
OCCXC = $(SCCXC:.c=.o)
OCCXMAIN = $(SCCXMAIN:.c=.o)

DIR=../../spooles.2.2

LIBS = \
       $(DIR)/spooles.a \
	../../ARPACK/libarpack_linux.a \
       -lpthread -lm -lc

ccx_2.20.exe: $(OCCXMAIN) ccx_2.20.a  $(LIBS)
	./date.pl; $(CC) $(CFLAGS) -c ccx_2.20.c; $(FC)  -Wall -O3 -fallow-argument-mismatch -o $@ $(OCCXMAIN) ccx_2.20.a $(LIBS) -fopenmp

ccx_2.20.a: $(OCCXF) $(OCCXC)
	ar vr $@ $?

clean:
	rm -f *.a *.o

