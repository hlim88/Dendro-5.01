.IGNORE:
.SUFFIXES : .F90 .f90 .F


F90_      = $(MPIF90)  $(F90_FREEFORM_FLAG) -g
F90_LOAD  = $(MPIF90)  $(LOCALIB) -O3 -g


.f90.o:
	$(F90_) -c $*.f90
.f.o:
	$(F90_) -c $*.f

EXECUTABLES = harmonic3 processc2p2 processclm uniformize unifandextrap

all: $(EXECUTABLES)

LIBS  =   -lbbhutil -lvutil $(LMPI) $(LSV)
LIBS2 =   -lbbhutil -lutilio -lvutil $(LMPI) $(LSV)
#LIBS =   -L/usr/local/intel81/lib -lrnpl -lxvs -L/usr/local/intel81/lib -lvutil -lodepack -llinpack -llapack -lblas $(LMPI) -lbbhutil $(LSV)

############################################################
# Note explicit statement of dependence of .o files 
# on various .inc files as well as corresponding .f files.
############################################################
#harmonic.o:    harmonic.f
#processc2p2.o: processc2p2.f90

harmonic: harmonic.o 
	$(F90_LOAD) harmonic.o $(LIBS) -o harmonic

harmonic2: harmonic2.o 
	$(F90_LOAD) harmonic2.o $(LIBS2) -o harmonic2

harmonic3: harmonic3.o
	$(F90_LOAD) harmonic3.o $(LIBS2) -o harmonic3

computeh: computeh.o 
	$(F90_LOAD) computeh.o $(LIBS2) -o computeh

processc2p2: processc2p2.o
	$(F90_LOAD) processc2p2.o $(LIBS2) -o processc2p2

processclm: processclm.o
	$(F90_LOAD) processclm.o $(LIBS2) -o processclm

uniformize: uniformize.o
	$(F90_LOAD) uniformize.o $(LIBS2) -o uniformize 

unifandextrap: unifandextrap.o
	$(F90_LOAD) unifandextrap.o $(LIBS2) -o unifandextrap

clean:
	rm *.o
	rm *.dat
	rm $(EXECUTABLES)
