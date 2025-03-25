## Created on 2024-07-21 at 13:42:32 CEST by David Gaspard <gaspard.dave@gmail.com> 
## Makefile of the Fortran program Ebsolve (v2) aiming at solving the Eilenberger equation.
##########################
## FILE PARAMETERS
##########################
PROGNAME = ebsolve
SRCDIR   = src
BINDIR   = bin

## Sources in dependency 'use' order (from low to high level):
SOURCES   = constants.f90 \
			base_utils.f90 \
			gauss_quad.f90 \
			modal_mesh.f90 \
			integrator.f90 \
			eilenberger_system.f90 \
			eilenberger_solver.f90 

MAINFILE = main.f90

SRCLIST = $(SOURCES:%.f90=src/%.f90)
BINLIST = $(SOURCES:%.f90=bin/%.o)

##########################
## COMPILATION OPTIONS
##########################
FORT    = gfortran
STD     = -std=f2008
OMP     = -fopenmp
DEBUG   = -g -O0 -fbacktrace -fcheck=all
## -fsanitize=address
FFLAGS  = -J$(BINDIR) -Wall -Wno-tabs $(OMP)
LIBS    = -llapack

all: $(PROGNAME)

$(PROGNAME): $(BINLIST) $(MAINFILE:%.f90=bin/%.o)
	$(FORT) $(FFLAGS) $^ $(LIBS) -o $@

$(BINDIR)/%.o: $(SRCDIR)/%.f90
	$(FORT) $(FFLAGS) -c $< -o $@

##########################
## TEST BUILD FUNCTIONS
##########################
TESTSRCLIST = $(shell find $(SRCDIR) -name "*.test.f90")
TESTEXELIST = $(TESTSRCLIST:$(SRCDIR)/%.f90=%)

test: $(BINLIST) $(TESTEXELIST)

%.test: $(BINLIST) bin/testing.o $(BINDIR)/%.test.o
	$(FORT) $(FFLAGS) $^ $(LIBS) -o $@

##################################
## CLEAN ALL BUILDS AND TESTS
##################################
clean:
	rm -rfv $(PROGNAME) $(BINDIR)/*.o $(BINDIR)/*.mod $(TESTEXELIST)

###### END OF FILE ######
