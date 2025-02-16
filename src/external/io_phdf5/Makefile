#makefile to build a wrf-phdf5 with netCDF

OBJSL   = wrf-phdf5.o wrf-phdf5attr.o wrf-phdf5support.o wrf-phdf5bdy.o
OBJS    = $(OBJSL)
OPTS    = -bmaxdata:0x80000000
FFLAGS  =  $(OPTS) -I$(PHDF5PATH)/lib 
FORTRANLIB=-I$(PHDF5PATH)/lib $(PHDF5PATH)/lib/libhdf5_fortran.a
LIBSHDF = $(FORTRANLIB) $(PHDF5PATH)/lib/libhdf5.a
LIB     = -lm -lz
CPP     = /lib/cpp 
CPPFLAGS = -C -P
M4      = m4 -Uinclude -Uindex -Ulen
AR      = ar
RANLIB  = ranlib

#-----------------------------------------------------------
#
# Added by Yunheng to include into ARPS
#
#  o Build PHDF5 static library to ARPS lib directory
#  o Add a new object wrf-phdf5bdy.o
#
#-----------------------------------------------------------

CFLAGS     =
ARFLAG     =
 
TOPDIR =
LIBDIR = $(TOPDIR)/lib

LIBPHDF5 = libwrfio_phdf5
 
#------ END of WYH -----------------------------------------

.SUFFIXES:      .F90 .f .o .code

all : $(LIBPHDF5)

$(LIBPHDF5): $(LIBDIR)/$(LIBPHDF5).a
	ls -l $(LIBDIR)/$(LIBPHDF5).a

$(LIBDIR)/$(LIBPHDF5).a:	$(OBJS)
			/bin/rm -f $(LIBDIR)/$(LIBPHDF5).a
			$(AR) $(ARFLAG) cr $(LIBDIR)/$(LIBPHDF5).a $(OBJSL)
			$(RANLIB) $(LIBDIR)/$(LIBPHDF5).a

wrf-phdf5support.o:      wrf-phdf5support.F90 
			$(CPP) $(CPPFLAGS) -I../inc wrf-phdf5support.F90 | $(M4) - > wrf-phdf5support.f
			$(FC) $(FFLAGS) -I../inc -c wrf-phdf5support.f
wrf-phdf5attr.o:         wrf-phdf5attr.F90 wrf-phdf5support.o
			$(CPP) $(CPPFLAGS)  wrf-phdf5attr.F90 | $(M4) - > wrf-phdf5attr.f
			$(FC) $(FFLAGS) -I../inc -c wrf-phdf5attr.f

wrf-phdf5.o:             wrf-phdf5.F90 wrf-phdf5attr.o wrf-phdf5support.o
			$(CPP) $(CPPFLAGS) wrf-phdf5.F90 | $(M4) - > wrf-phdf5.f
			$(FC) $(FFLAGS) -I../inc -c wrf-phdf5.f

wrf-phdf5bdy.o:          wrf-phdf5bdy.F90 wrf-phdf5attr.o wrf-phdf5support.o
			$(CPP) $(CPPFLAGS) wrf-phdf5bdy.F90 | $(M4) - > wrf-phdf5bdy.f
			$(FC) $(FFLAGS) -I../inc -c wrf-phdf5bdy.f

superclean:
			/bin/rm -f *.f *.o *.mod $(LIBDIR)/$(LIBPHDF5).a
