#
#     ##################################################################
#     ##################################################################
#     ######                                                      ######
#     ######      Advanced Regional Prediction System (ARPS)      ######
#     ######            (ZXPLOT library Version 3.0               ######
#     ######                                                      ######
#     ######                     Developed by                     ######
#     ######              Ming Xue and Zuojun Zhang               ######
#     ######                                                      ######
#     ##################################################################
#     ##################################################################
#
#=======================================================================
#
#  PURPOSE: This makefile generates the ZXPLOT libraries for ARPS 
#           executables, such as arpspltpost, arpspltncar, arpstrn etc.
#
#  AUTHOR:  Yunheng Wang
#           11/17/2003
#
#  Modification history:
#
#  Option "-O0" was hard-coded.
#
#=======================================================================

#-----------------------------------------------------------------------
#
# Default names of the loader and tar, can be replaced on the command
# line for other system
#
#-----------------------------------------------------------------------

FTN   = f90
FTNMP = f90
LDR   = f90
CC    = cc

ARFLAG =

ARPS_LD = $(LDR)

TAR   = tar
AWK   = awk
RM    = rm
LN    = ln
AR    = ar

TOPDIR  = 
BINDIR  = $(TOPDIR)
INCLDIR = $(TOPDIR)/include
LIBDIR  = $(TOPDIR)/lib

#-----------------------------------------------------------------------
#
# Compiler Flag of Options. The default is for Intel Fortran Compiler
# on Linux platform
#
#-----------------------------------------------------------------------

FFLAGS =
FFLAGS_main =
CFLAGS =
LDFLAGS =
ICEFLAG =
RADFLAG =
THMFLAG =

FIXFLAGS  =
FREEFLAGS =

#-----------------------------------------------------------------------
#
# Dependencies
#
#-----------------------------------------------------------------------

.SUFFIXES: $(SUFFIXES) .f90

.f.o:
	$(FTN) $(FFLAGS) $(FIXFLAGS)  -c $<
.f90.o:
	$(FTN) $(FFLAGS) $(FREEFLAGS) -c $<
.c.o:
	$(CC)  $(CFLAGS) -c $<

#-----------------------------------------------------------------------
#
# Object library to be generated for ARPS executables:
#
# LIBZXPOST    = libzxpost       ZXPLOT library for PS file format
# LIBZXNCAR    = libzxncar       ZXPLOT library for Ncarg format
#
#-----------------------------------------------------------------------

LIBZXPOST    = libzxpost
LIBZXNCAR    = libzxncar

#-----------------------------------------------------------------------
#
# List of machine-dependent object codes
#
#-----------------------------------------------------------------------

LIBZXPOSTOBJS = xpost3.o zxplot3.o zxpslib3.o

LIBZXNCAROBJS = xncar3.o zxnglib3.o zxplot3.o

#-----------------------------------------------------------------------
#
# Set Default
#
#-----------------------------------------------------------------------

default: $(LIBZXPOST)

#-----------------------------------------------------------------------
#
# Make ZXPLOT libraries.
#
#-----------------------------------------------------------------------

$(LIBZXPOST): $(LIBDIR)/$(LIBZXPOST)
	ls -l $(LIBDIR)/$(LIBZXPOST).a

$(LIBDIR)/$(LIBZXPOST): $(LIBZXPOSTOBJS)
	$(RM) -f $@.a
	$(AR) $(ARFLAG) rc $@.a $(LIBZXPOSTOBJS)

$(LIBZXNCAR): $(LIBDIR)/$(LIBZXNCAR)
	ls -l $(LIBDIR)/$(LIBZXNCAR).a

$(LIBDIR)/$(LIBZXNCAR): $(LIBZXNCAROBJS)
	$(RM) -f $@.a
	$(AR) $(ARFLAG) rc $@.a $(LIBZXNCAROBJS)

#-----------------------------------------------------------------------
#
# Remove the object code for individual programs
#
#-----------------------------------------------------------------------

clean: clean.libzxpost clean.libzxncar

clean.libzxpost:
	-$(RM) -f $(LIBZXPOSTOBJS) $(LIBDIR)/$(LIBZXPOST).a      \
                  *.mod *.d work.pc work.pcl

clean.libzxncar:
	-$(RM) -f $(LIBZXNCAROBJS) $(LIBDIR)/$(LIBZXNCAR).a      \
                  *.mod *.d work.pc work.pcl

#-----------------------------------------------------------------------
#
# Object code dependency list:
#
#	$(FTN) $(FIXFLAGS) -O0 -c $<
#
#-----------------------------------------------------------------------

xpost3.o    : xpost3.f
zxplot3.o   : zxplot3.f
zxpslib3.o  : zxpslib3.f
xncar3.o    : xncar3.f
zxnglib3.o  : zxnglib3.f
