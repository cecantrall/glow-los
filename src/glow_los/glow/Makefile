#
FC = gfortran

 FFLAGS = -O3 -ffree-line-length-none  -std=legacy -g -O0 -fbacktrace -I $(INC_NETCDF) $(LIBS) -I .
#FFLAGS = -g $(DBGFLAGS)

DBGFLAGS = -debug full -traceback
DBGFLAGS += -check bounds -check format -check output_conversion -check pointers -check uninit
DBGFLAGS += -fpe-all=0 # this traps all floating point exceptions

.SUFFIXES: .o .F .F90 .f .f90

GLOWdir = .
%.o: %.f90
	$(FC) $(FFLAGS) -c  -o $@ $<
%.o: %.f
	$(FC) $(FFLAGS) -c  -o $@ $<

# Sources (in order of dependency):
#
SOURCES = $(GLOWdir)/cglow.f90 $(GLOWdir)/glowdriver.f90 $(GLOWdir)/glow.f90 $(GLOWdir)/bands.f90 $(GLOWdir)/egrid.f90 $(GLOWdir)/ephoto.f90 $(GLOWdir)/etrans.f90 $(GLOWdir)/exsect.f $(GLOWdir)/fieldm.f $(GLOWdir)/gchem.f90 $(GLOWdir)/maxt.f90 $(GLOWdir)/qback.f90 $(GLOWdir)/rcolum.f90 $(GLOWdir)/solzen.f90 $(GLOWdir)/ssflux.f90 

OBJS := $(addsuffix .o, $(basename $(SOURCES)))
EXEC = testLBHcol.exe

$(EXEC): $(OBJS)
	$(FC) -o $@ $(OBJS) $(LIBS) $(LDFLAGS)

#LIB_NETCDF = /opt/homebrew/Cellar/netcdf-fortran/4.6.0/lib
#INC_NETCDF = /opt/homebrew/Cellar/netcdf-fortran/4.6.0/include

LIB_NETCDF = /opt/homebrew/Cellar/netcdf-fortran/4.6.2/lib
INC_NETCDF = /opt/homebrew/Cellar/netcdf-fortran/4.6.2/include

LIB_NETCDF = /opt/homebrew/lib
INC_NETCDF = /opt/homebrew/include

LIBS       = -L$(LIB_NETCDF) -lnetcdff

clean:
	rm -f $(GLOWdir)/*.o $(GLOWdir)/*.mod $(EXEC)
	rm -f *.o *.mod $(EXEC)

