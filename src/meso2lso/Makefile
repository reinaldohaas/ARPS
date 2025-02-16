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


#-----------------------------------------------------------------------
#
# Default names of the loader and tar, can be replaced on the command
# line for other system
#
#-----------------------------------------------------------------------

FTN   = ifort
LDR   = $(FTN)
CC    = cc

ARPS_LD = $(LDR)

RM  = rm

TOPDIR  = ./
BINDIR  = $(TOPDIR)/bin
INCLDIR = $(TOPDIR)/include
LIBDIR  = $(TOPDIR)/lib
MODDIR  = $(TOPDIR)/modules

#-----------------------------------------------------------------------
#
# Compiler Flag of Options. 
#
#-----------------------------------------------------------------------

FFLAGS =
CFLAGS = -DUNDERSCORE
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
# List of Objects
#
#-----------------------------------------------------------------------

OBJS = readmeso.o readmdf.o readwtx.o    \
       mkitime.o getime.o cvtitoi.o citime.o strlnth.o \
       upcase.o

#-----------------------------------------------------------------------
#
# Set Default
#
#-----------------------------------------------------------------------

PROG = meso2lso

default: $(PROG)

#-----------------------------------------------------------------------
#
# Compile and link executable, meso2lso
#
#-----------------------------------------------------------------------

$(PROG):   $(BINDIR)/$(PROG)
	ls -l $(BINDIR)/$(PROG)

$(BINDIR)/$(PROG): $(OBJS) $(PROG).o 
	$(ARPS_LD) $(LDFLAGS) -o $@ $(PROG).o $(OBJS) $(LIBS)
clean:
	$(RM) -f $(BINDIR)/$(PROG) $(PROG).o $(OBJS)
