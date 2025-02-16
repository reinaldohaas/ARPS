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
#  PURPOSE: This makefile generates the Zmosaic2arps executable.
#
#
#  AUTHOR:  Yunheng Wang
#           02/01/2006
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

FTN   = f90
FTNMP = $(FTN)
LDR   = $(FTN)
CC    = cc
CPP   = cpp

ARPS_LD = $(LDR)

TAR = tar
AWK = awk
RM  = rm
LN  = ln
AR  = ar

TOPDIR  =
BINDIR  = $(TOPDIR)
INCLDIR = $(TOPDIR)/include
LIBDIR  = $(TOPDIR)/lib

#-----------------------------------------------------------------------
#
# Compiler Flag of Options. The default is for AIX Fortran xlf.
#
#-----------------------------------------------------------------------

FFLAGS =
FFLAGS_main =
CFLAGS =
LDFLAGS =
ICEFLAG =
RADFLAG =
THMFLAG =

CPPFLAGS  = -C -P
FIXFLAGS  =
FREEFLAGS =

ARFLAG =


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
.mod.o:     # .mod --> .o dependency dropped


#-----------------------------------------------------------------------
#
# ARPS library to be linked
#
# LIBARPS = libarps      ARPS shared library
#
#-----------------------------------------------------------------------

LIBARPS = libarps

#----- object files -----

SRCSF90  = Zmosaic2arps.f90 Zmosaiclib.f90

OBJS   	 = ${SRCSF90:.f90=.o}

CMD	= Zmosaic2arps

#----- build commands -----

$(CMD):     $(BINDIR)/$(CMD)
	ls -l $(BINDIR)/$(CMD)

$(BINDIR)/$(CMD):	$(OBJS) $(LIBDIR)/$(LIBARPS).a
	$(ARPS_LD) $(LDFLAGS) -o $@ $(OBJS) $(LIBDIR)/$(LIBARPS).a $(LIBS)

clean:
	 -$(RM) -f $(OBJS) $(BINDIR)/$(CMD)
