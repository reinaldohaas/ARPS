#makefile to build io_int that does binary i/o

OBJSL   = io_int.o module_internal_header_util.o pack_utils.o
OBJS    = $(OBJSL)
OPTS    =
FREEFLAGS = 
FFLAGS  =  $(OPTS)
CFLAGS  = 
CC      =
FC      =
LIBS    = 
CPP     = cpp
CPPFLAGS = -C -P
M4      = m4 -Uinclude -Uindex -Ulen
AR      = ar
#RANLIB	= ranlib
RANLIB  = echo
RM      = rm -f

#-----------------------------------------------------------
#
# Added by Yunheng to include into ARPS
#
#   Build internal static library to ARPS lib directory
#
#-----------------------------------------------------------
CFLAGS     =
ARFLAG     =
 
TOPDIR =
LIBDIR = $(TOPDIR)/lib
 
LIBWRFINT = libwrfio_int
 
#------ END of WYH -----------------------------------------

.SUFFIXES:      .F90 .F .f90 .o .code .c

all : $(LIBWRFINT)

$(LIBWRFINT): $(LIBDIR)/$(LIBWRFINT).a
	ls -l $(LIBDIR)/$(LIBWRFINT).a

$(LIBDIR)/$(LIBWRFINT).a:	$(OBJS)
			/bin/rm -f $(LIBDIR)/$(LIBWRFINT).a
			$(AR) $(ARFLAG) cr $(LIBDIR)/$(LIBWRFINT).a $(OBJSL)
			$(RANLIB) $(LIBDIR)/$(LIBWRFINT).a

io_int.o:               io_int.F module_internal_header_util.o
			$(CPP) $(CPPFLAGS) -I../inc io_int.F | $(M4) - > io_int.f90
			$(FC) $(FFLAGS) $(FREEFLAGS) -I../inc -c io_int.f90

pack_utils.o:           pack_utils.c
	$(RM) $@
	$(CC) -c $(CFLAGS) $*.c

# this is just so we can get the .mod file, it will be recompiled in frame again as part of WRF frmwk

module_internal_header_util.o : module_internal_header_util.F
			$(CPP) $(CPPFLAGS) -I../inc module_internal_header_util.F > module_internal_header_util.f90
			$(FC) $(FFLAGS) $(FREEFLAGS) -I../inc -c module_internal_header_util.f90
			
superclean:
			/bin/rm -f *.f90 *.o *.mod $(LIBDIR)/$(LIBWRFINT).a
