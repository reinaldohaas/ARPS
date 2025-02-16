#makefile to build a wrf_io with netCDF

OBJSL   = wrf_io.o field_routines.o module_wrfsi_static.o
OBJS    = $(OBJSL)
CODE    = ext_ncd_get_dom_ti.code ext_ncd_get_var_td.code ext_ncd_get_var_ti.code ext_ncd_put_dom_ti.code ext_ncd_put_var_td.code ext_ncd_put_var_ti.code transpose.code
FFLAGS  = $(FFLAGS) $(FREEFLAGS)
LIBS    = 
CPPFLAGS   = -C -P
M4      = m4 -Uinclude -Uindex -Ulen
AR      = ar
ARFLAG  =

TOPDIR  =
LIBDIR  = $(TOPDIR)/lib
LIBWRFNET = libwrfio_net

.SUFFIXES:      .F90 .f .o .code

all : $(LIBWRFNET)

$(LIBWRFNET): $(LIBDIR)/$(LIBWRFNET).a
	ls -l $(LIBDIR)/$(LIBWRFNET).a

$(LIBDIR)/$(LIBWRFNET).a: $(OBJS) bitwise_operators.o $(CODE)
	/bin/rm -f $@
	$(AR) $(ARFLAG) cr $(LIBDIR)/$(LIBWRFNET).a $(OBJSL) bitwise_operators.o

wrf_io.o:   wrf_io.F bitwise_operators.o $(CODE)
	$(CPP) $(CPPFLAGS) -I../inc -DWRFIO_NCD_LARGE_FILE_SUPPORT wrf_io.F | $(M4) - > wrf_io.f90
	$(FC) $(FFLAGS) -I../inc -c wrf_io.f90

module_wrfsi_static.o:  module_wrfsi_static.F
	$(CPP) $(CPPFLAGS) -I../inc module_wrfsi_static.F > module_wrfsi_static.f90
	$(FC) $(FFLAGS) -c module_wrfsi_static.f90

bitwise_operators.o : bitwise_operators.c
	$(CC) -c $(CFLAGS) bitwise_operators.c

field_routines.o: field_routines.F
	$(CPP) $(CPPFLAGS) -I../inc field_routines.F > field_routines.f90
	$(FC) $(FFLAGS) -I../inc -c field_routines.f90

superclean:
	/bin/rm -f *.f90 *.f *.o *.mod $(LIBDIR)/$(LIBWRFNET).a
