
GSLDIR = /opt/local/

CC_COMMAND  = gcc
CCFLAGS     = -g
INCFLAGS    = -I$(GSLDIR)/include
LIBPATHS    = -L$(GSLDIR)/lib
CC          = $(CC_COMMAND) $(CCFLAGS) $(INCFLAGS)
CC_LOAD     = $(CC_COMMAND) $(CCFLAGS) $(LIBPATHS)

LIBS  = -lgsl -lgslcblas -lm


#OBJ = driver.o CoordTransf.o Equations.o FuncAndJacobian.o Newton.o ParamCheck.o TP_utilities.o TwoPunctures.o
OBJ = driver.o Newton.o TP_utilities.o CoordTransf.o Equations.o FuncAndJacobian.o TwoPunctures.o approx.o

.c.o:
	$(CC) -c $*.c

tpid: $(OBJ)
	$(CC_LOAD) $^ $(LIBS) -o $@

clean:
	/bin/rm *.o tpid *_x.dat *_y.dat
