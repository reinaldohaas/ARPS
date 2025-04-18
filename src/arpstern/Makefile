#
#     ##################################################################
#     ##################################################################
#     ######                                                      ######
#     ######      Advanced Regional Prediction System (ARPS)      ######
#     ######                   Version 5.0                        ######
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
#  PURPOSE: This makefile generates the ARPSTER executables.
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
TRNDIR = $(TOPDIR)/src/arpstrn

#-----------------------------------------------------------------------
#
# Compiler Flag of Options
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
# Executable commands to be generated by this make file:
#
# ARPSTRNEXE = arpstern      ARPS data converter.
#
# DIR1DEGEXE  = dir1deg    Convertor of terrain data to direct access
#                          file for 1 degree resolution
# DIR5MINEXE  = dir5min    Convertor for 5 minutes resolution
# DIR30SECEXE = dir30sec   Convertor for 30 seconds resolution
#
#-----------------------------------------------------------------------

ARPSTRNEXE  = arpstern

DIR1DEGEXE  = dir1deg
DIR5MINEXE  = dir5min
DIR30SECEXE = dir30sec

#-----------------------------------------------------------------------
#
# ARPS shared library
#
#-----------------------------------------------------------------------

LIBARPS = libarps

#-----------------------------------------------------------------------
#
# List of Makefiles
#
#-----------------------------------------------------------------------

MAKEFILE = Makefile

#-----------------------------------------------------------------------
#
# List of machine-dependent object codes
#
#-----------------------------------------------------------------------

ARPSTRNOBJ = arpstern.o 

DIR1DEGOBJ  = dir1deg.o
DIR5MINOBJ  = dir5min.o
DIR30SECOBJ = dir30sec.o

#-----------------------------------------------------------------------
#
# Set Default
#
#-----------------------------------------------------------------------

default: $(ARPSTRNEXE)

#-----------------------------------------------------------------------
#
# Compile and link the executable, arpstern
#
#-----------------------------------------------------------------------

$(ARPSTRNEXE):  $(BINDIR)/$(ARPSTRNEXE)
	ls -l $(BINDIR)/$(ARPSTRNEXE)

$(BINDIR)/$(ARPSTRNEXE): $(ARPSTRNOBJ) $(LIBDIR)/$(LIBARPS).a 
	$(ARPS_LD) $(LDFLAGS) -o $@ $(ARPSTRNOBJ)                     \
                               $(LIBDIR)/$(LIBARPS).a $(LIBS)

$(DIR1DEGEXE): $(DIR1DEGOBJ)
	$(ARPS_LD) $(LDFLAGS) -o $(BINDIR)/$@ $(DIR1DEGOBJ)

$(DIR5MINEXE): $(DIR5MINOBJ)
	$(ARPS_LD) $(LDFLAGS) -o $(BINDIR)/$@ $(DIR5MINOBJ)

$(DIR30SECEXE): $(DIR30SECOBJ)
	$(ARPS_LD) $(LDFLAGS) -o $(BINDIR)/$@ $(DIR30SECOBJ)

#-----------------------------------------------------------------------
#
# Remove the object code for individual programs
#
#-----------------------------------------------------------------------

clean: clean.arpstern

clean.arpstern:
	-$(RM) -f $(BINDIR)/$(ARPSTRNEXE)  $(ARPSTRNOBJ)

clean.dirtern:
	-$(RM) -f $(BINDIR)/$(DIR1DEGEXE)  $(DIR1DEGOBJ)  \
                  $(BINDIR)/$(DIR5MINEXE)  $(DIR5MINOBJ)  \
                  $(BINDIR)/$(DIR30SECEXE) $(DIR30SECOBJ)

#-----------------------------------------------------------------------
#
# Object code dependency list:
#
#-----------------------------------------------------------------------

arpstern.o  : arpstern.f90

dir1deg.o   : dir1deg.f90
dir5min.o   : dir5min.f90
dir30sec.o  : dir30sec.f90

