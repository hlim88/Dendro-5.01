
GSLDIR = $(HOME)/apps/gsl

CXX_COMMAND  = g++-6
CXXFLAGS     = -g
INCFLAGS    = -I$(GSLDIR)/include
LIBPATHS    = -L$(GSLDIR)/lib
CXX         = $(CXX_COMMAND) $(CXXFLAGS) $(INCFLAGS)
CXX_LOAD    = $(CXX_COMMAND) $(CXXFLAGS) $(LIBPATHS)

LIBS  = -lgsl -lgslcblas -lm


#OBJ = driver.o CoordTransf.o Equations.o FuncAndJacobian.o Newton.o ParamCheck.o TP_utilities.o TwoPunctures.o
OBJ = driver.o params.o Newton.o TP_utilities.o CoordTransf.o Equations.o FuncAndJacobian.o TwoPunctures.o approx.o

.cpp.o:
	$(CXX) -c $*.cpp

tpid: $(OBJ)
	$(CXX_LOAD) $^ $(LIBS) -o $@

clean:
	/bin/rm *.o tpid *_x.dat *_y.dat
