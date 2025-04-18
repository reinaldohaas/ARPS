#
#     ##################################################################
#     ##################################################################
#     ######                                                      ######
#     ######      Advanced Regional Prediction System (ARPS)      ######
#     ######                   Version 4.3                        ######
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
#  PURPOSE: This makefile generates the ARPSVTOBS executable
#
#  AUTHOR:  Yuhe Liu
#           7/22/1997
#
#  Modification history:
#
#
#  OTHER INFORMATION:
#       See the makearps command.
#
#=======================================================================

#-----------------------------------------------------------------------
#
# Default shell
#
#-----------------------------------------------------------------------

# SHELL=/bin/csh

#-----------------------------------------------------------------------
#
# Default names of the loader and tar, can be replaced on the command
# line for other system
#
#-----------------------------------------------------------------------

FTN = f90
LDR = f90
CC  = cc

ARPS_LD = $(LDR)

TAR = tar
AWK = awk
RM  = rm
LN  = ln

TOPDIR =
BINDIR = $(TOPDIR)
INCLDIR = $(TOPDIR)/include
LIBDIR = $(TOPDIR)/lib

INCFILE = /dev/null

include $(INCFILE)

#-----------------------------------------------------------------------
#
# Compiler Flag of Options. The default is for AIX Fortran xlf.
#
#-----------------------------------------------------------------------

FFLAGS =
CFLAGS =
LDFLAGS =

FIXFLAGS  =
FREEFLAGS =

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

#-----------------------------------------------------------------------
#
# Executable to be generated by this make file:
#
# ARPSCVTOBSEXE = arpscvtobs     ARPSCVTOBS executable
#
#-----------------------------------------------------------------------

ARPSCVTOBSEXE = arpscvtobs

#-----------------------------------------------------------------------
#
# Shared library to be linked
#
# LIBARPS = libarps      ARPS shared library
#
#-----------------------------------------------------------------------

LIBARPS = libarps

#-----------------------------------------------------------------------
#
# List of machine-dependent object codes
#
#-----------------------------------------------------------------------

ARPSCVTOBSOBJ = arpscvtobs.o module_verifhdf.o wrtlsohdf.o

#-----------------------------------------------------------------------
#
# Set Default
#
#-----------------------------------------------------------------------

default: $(ARPSCVTOBSEXE)

#-----------------------------------------------------------------------
#
# Compile and link ARPS model executable, arps
#
#-----------------------------------------------------------------------

$(ARPSCVTOBSEXE): $(ARPSCVTOBSOBJ)
	$(ARPS_LD) $(LDFLAGS) -o $(BINDIR)/$@ $(ARPSCVTOBSOBJ) \
			       $(LIBS)

#-----------------------------------------------------------------------
#
# Remove the object code for individual programs
#
#-----------------------------------------------------------------------

clean.arpscvtobs:
	-$(RM) -f $(BINDIR)/$(ARPSCVTOBSEXE) $(ARPSCVTOBSOBJ)

#-----------------------------------------------------------------------
#
# Object code dependency list:
#
#-----------------------------------------------------------------------

arpscvtobs.o : arpscvtobs.f90     $(INCLDIR)/globcst.inc \
                              $(INCLDIR)/phycst.inc \
			      $(INCLDIR)/grid.inc \
                              $(INCLDIR)/adas.inc

#module_verifhdf.o : module_verifhdf.f90		$(INCLDIR)/hdf.f90

wrtlsohdf.o:	wrtlsohdf.f90
