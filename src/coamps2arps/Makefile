#
#     ##################################################################
#     ##################################################################
#     ######                                                      ######
#     ######      Advanced Regional Prediction System (ARPS)      ######
#     ######                   Version 5.2                        ######
#     ######                                                      ######
#     ######                     Developed by                     ######
#     ######     Center for Analysis and Prediction of Storms     ######
#     ######                University of Oklahoma                ######
#     ######                                                      ######
#     ##################################################################
#     ##################################################################
#
#=======================================================================
#
#  PURPOSE: This makefile generates the coamps2arps executable
#
#  AUTHOR:  Yunheng Wang
#           05/06/2011
#
#  Modification history:
#
#
#  OTHER INFORMATION:
#       See the makearps command.
#
#=======================================================================

INCFILE = /dev/null

#include $(INCFILE)

#-----------------------------------------------------------------------
#
# Default names of the loader and tar, can be replaced on the command
# line for other system
#
#-----------------------------------------------------------------------

FTN = ifort
FTNMP = $(FTN)
LDR = $(FTN)
CC  = gcc
CPP = cpp

ARPS_LD = $(LDR)

TAR = tar
AWK = awk
RM  = rm
LN  = ln

TOPDIR     = ../..
BINDIR     = $(TOPDIR)/bin
INCLDIR    = $(TOPDIR)/include
ARPSLIBDIR = $(TOPDIR)/lib
MODDIR     = $(TOPDIR)/modules

ENSDIR = $(TOPDIR)/src/arpsens
A2CDIR = $(TOPDIR)/src/arps2coamps
COMDIR = $(TOPDIR)/src/common

#-----------------------------------------------------------------------
#
# Compiler Flag of Options. The default is for Intel Fortran on Linux.
#
#-----------------------------------------------------------------------

FFLAGS = -I$(INCLDIR)
CFLAGS =
LDFLAGS =

FIXFLAGS  =
FREEFLAGS =
CPPFLAGS  = -C -P

#-----------------------------------------------------------------------
#
# Dependencies
#
#-----------------------------------------------------------------------

.SUFFIXES: $(SUFFIXES) .f90

.f.o:
	$(FTN) $(FFLAGS) $(FIXFLAGS) -c $<
.f90.o:
	$(FTN) $(FFLAGS) $(FREEFLAGS) -c $<
.c.o:
	$(CC)  $(CFLAGS) -c $<
.F.o:
	$(CPP) $(CPPFLAGS) $< $(addsuffix .f90, $(basename $<));
	$(FTNMP) $(FFLAGS) $(FREEFLAGS) -c $(addsuffix .f90, $(basename $<))
.mod.o:        # .mod --> .o dependency dropped

#-----------------------------------------------------------------------
#
# Executables to be generated by this make file:
#
# C2AEXE   = coamps2arps          COAMPS2ARPS executable
# C2AMPEXE = coamps2arps_mpi      COAMPS2ARPS mpi executable
#
#-----------------------------------------------------------------------

C2AEXE   = coamps2arps
C2AMPEXE = coamps2arps_mpi

#-----------------------------------------------------------------------
#
# ARPS library and WRF library to be linked
#
# LIBARPS = libarps      ARPS shared library
#
#-----------------------------------------------------------------------

LIBARPS = libarps
LIBADAS = libadas

#-----------------------------------------------------------------------
#
# List of machine-dependent object codes
#
#-----------------------------------------------------------------------

MPIOBJ_ON  = coamps_parallel_module.o
MPIOBJ_OFF = coamps_parallel_module.o

MPIOBJS    = $(MPIOBJ_OFF)

#-----------------------------------------------------------------------
#
# List of code objects
#
#-----------------------------------------------------------------------

C2AMAINOBJ = coamps2arps.o
C2AMODULES = module_coamps2arps_namelist.o module_coampsgrid.o   \
             module_interpolation.o
C2AOBJS1   = process_file.o interplib.o

EXTMODULES = module_arpsgrid_constants.o module_dump_arpsgrid.o \
             module_coamps_constants.o

C2AIMPORT  = coamps_mapprojs_module.o

GEMOBJ_ON  = gemio.o
GEMOBJ_OFF = nogemio.o
GEMOBJ     = $(GEMOBJ_OFF)

CRTM_ON   = crtmpost.o
CRTM_OFF  = crtm_no.o
CRTMOBJ   = $(CRTM_OFF)

CITM_ON   =
CITM_OFF  = $(ENSDIR)/citm_no.o
CITMOBJ   = $(CITM_OFF)

POSTOBJS   = $(ENSDIR)/postcore.o $(ENSDIR)/arpsenslib.o      \
             $(ENSDIR)/arpspostlib.o $(ENSDIR)/extrefleclib.o \
             $(ENSDIR)/$(GEMOBJ) $(ENSDIR)/$(CRTMOBJ) $(CITMOBJ)

C2AOBJS = $(C2AIMPORT) $(EXTMODULES) $(C2AMODULES) \
          $(MPIOBJS)   $(C2AOBJS1)   $(C2AMAINOBJ)  $(POSTOBJS)

#-----------------------------------------------------------------------
#
# Set Default
#
#-----------------------------------------------------------------------

default: $(C2AEXE)

#-----------------------------------------------------------------------
#
# Compile and link executable
#
#-----------------------------------------------------------------------

$(C2AEXE):      $(BINDIR)/$(C2AEXE)
	@$(RM) coamps_parallel_module.f90 coamps_parallel_module.o
	ls -l $(BINDIR)/$(C2AEXE)

$(C2AMPEXE):    $(BINDIR)/$(C2AMPEXE)
	@$(RM) coamps_parallel_module.f90 coamps_parallel_module.o
	ls -l $(BINDIR)/$(C2AMPEXE)

$(BINDIR)/$(C2AMPEXE): $(C2AOBJS) $(ARPSLIBDIR)/$(LIBARPS).a
	$(ARPS_LD) $(LDFLAGS) -o $@ $(C2AOBJS)            \
                $(ARPSLIBDIR)/$(LIBADAS).a $(ARPSLIBDIR)/$(LIBARPS).a       \
        	$(LIBS)

$(BINDIR)/$(C2AEXE): $(C2AOBJS) $(ARPSLIBDIR)/$(LIBARPS).a
	$(ARPS_LD) $(LDFLAGS) -o $@ $(C2AOBJS)            \
                $(ARPSLIBDIR)/$(LIBADAS).a $(ARPSLIBDIR)/$(LIBARPS).a       \
        	$(LIBS)

#-----------------------------------------------------------------------
#
# Remove the object code for individual programs
#
#-----------------------------------------------------------------------

clean.coamps2arps:
	-$(RM) -f $(BINDIR)/$(C2AEXE) $(BINDIR)/$(C2AMPEXE)  \
                  $(C2AOBJS) $(MPIOBJ_ON) $(MPIOBJ_OFF)

#-----------------------------------------------------------------------
#
# Object code dependency list:
#
#-----------------------------------------------------------------------

coamps2arps.o : coamps2arps.f90 module_coamps_constants.o module_coamps2arps_namelist.o

module_coamps2arps_namelist.o  : module_coamps2arps_namelist.f90

module_coampsgrid.o : module_coampsgrid.f90 module_coamps_constants.o          \
                      coamps_parallel_module.o coamps_mapprojs_module.o

module_coamps_constants.o    : $(A2CDIR)/module_coamps_constants.f90
	$(FTN) $(FFLAGS) $(FREEFLAGS) -c $<
coamps_mapprojs_module.o : $(A2CDIR)/coamps_mapprojs_module.f90
	$(FTN) $(FFLAGS) $(FREEFLAGS) -c $<
coamps_parallel_module.o : $(A2CDIR)/coamps_parallel_module.F
	$(CPP) $(CPPFLAGS) $< coamps_parallel_module.f90
	$(FTNMP) $(FFLAGS) $(FREEFLAGS) -c coamps_parallel_module.f90

module_arpsgrid_constants.o : $(COMDIR)/module_arpsgrid_constants.f90
	$(FTN) $(FFLAGS) $(FREEFLAGS) -c $<

module_dump_arpsgrid.o : $(COMDIR)/module_dump_arpsgrid.f90
	$(FTN) $(FFLAGS) $(FREEFLAGS) -c $<

$(ENSDIR)/postcore.o    : $(ENSDIR)/postcore.f90
	$(FTN) $(FFLAGS) $(FREEFLAGS) -o $@ -c $<
$(ENSDIR)/arpsenslib.o  : $(ENSDIR)/arpsenslib.f90
	$(FTN) $(FFLAGS) $(FREEFLAGS) -o $@ -c $<
$(ENSDIR)/arpspostlib.o : $(ENSDIR)/arpspostlib.f90
	$(FTN) $(FFLAGS) $(FREEFLAGS) -o $@ -c $<
$(ENSDIR)/extrefleclib.o : $(ENSDIR)/extrefleclib.f90
	$(FTN) $(FFLAGS) $(FREEFLAGS) -o $@ -c $<

$(ENSDIR)/gemio.o       : $(ENSDIR)/gemio.f90
	@$(RM) -f $(ENSDIR)/GEMPRM.PRM
	@if [ "$(NAWIPS)" != "" ]; then \
        sed -e 's/^C/!/;/MPAORT/s/= 3,/= 3, \&/;/MPALAM/s/+//;s/INCLUDE/!INCLUDE/' $(NAWIPS)/gempak/include/GEMPRM.PRM > $(ENSDIR)/GEMPRM.PRM; \
        fi
	$(FTN) $(FFLAGS) $(FREEFLAGS) -o $@ -c $<

$(ENSDIR)/nogemio.o     : $(ENSDIR)/nogemio.f90
	$(FTN) $(FFLAGS) $(FREEFLAGS) -o $@ -c $<
$(ENSDIR)/crtmpost.o : $(ENSDIR)/crtmpost.f90
	$(FTN) $(FFLAGS) $(FREEFLAGS) -o $@ -c $<
$(ENSDIR)/crtm_no.o  : $(ENSDIR)/crtm_no.f90
	$(FTN) $(FFLAGS) $(FREEFLAGS) -o $@ -c $<
$(ENSDIR)/citm_no.o  : $(ENSDIR)/citm_no.f90
	$(FTN) $(FFLAGS) $(FREEFLAGS) -o $@ -c $<

