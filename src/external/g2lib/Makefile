SHELL=/bin/sh

#  Make sure one of the following options appears in your CFLAGS
#  variable to indicate which system you are on.  Used to call
#  "C" routines from Fortran.
#  -DLINUX, -DLINUXG95, -DSGI, -DHP, -DCRAY90, -DAIX, -DLINUXF90, -DVPP5000

#  If you want to enable support for PNG or JPEG2000 encoding/decoding,
#  you must specify -DUSE_PNG and/or -DUSE_JPEG2000 in the FDEFS variable
#  for the Fortran pre-processor
#  -DUSE_PNG requires libpng.a and libz.a
#                     ( and png.h pngconf.h zconf.h zlib.h include files).
#  -DUSE_JPEG2000 requires libjasper.a
#                     ( and all the jasper/*.h include files).
#
#  In addition, INCDIR must include all directories where the above
#  mentioned include files can be found.
CPPFLAGS = -C -P
#FDEFS=-DUSE_PNG -DUSE_JPEG2000 $(CPPFLAGS)

TOPDIR = 
LIBDIR = $(TOPDIR)/lib
INCDIR = -I$(TOPDIR)/include

LIBG2  = libg2
MODDIR = $(TOPDIR)/lib

#--------------------------------------
# The following was used for XLF on AIX
#DEFS=-DAIX -DHAVE_SYS_TYPES_H=1
#FC=xlf
#CC=xlc
#CPP=/usr/ccs/lib/cpp -P
#MODDIR=../g2mod
#FFLAGS=-O3 -g -qnosave -qarch=auto -qmoddir=$(MODDIR) -I $(MODDIR)
#CFLAGS=-O3 -q64 -g -qcpluscmt -qarch=auto $(DEFS) $(INCDIR)
#ARFLAGS=-X64
#--------------------------------------
# The following was used for G95 on LINUX
#
#DEFS=-DLINUXG95
FC=g95
CC=cc
CPP=cpp
FFLAGS=-O3 -I $(MODDIR)
CFLAGS=-O3 $(INCDIR)
ARFLAGS=
FIXFLAGS =
FREEFLAGS =

#--------------------------------------

.SUFFIXES: .a .f .F .c

.F.f:
	$(CPP) $(CPPFLAGS) $*.F $*.f
.f.o:
	$(FC) -c $(FFLAGS) $(FIXFLAGS) $<
.c.o:
	$(CC) -c $(CFLAGS) $<

OBJS =  gridtemplates.o pdstemplates.o drstemplates.o \
        gribmod.o realloc.o addfield.o addgrid.o addlocal.o \
        getfield.o gb_info.o gf_getfld.o gf_free.o \
	gf_unpack1.o gf_unpack2.o gf_unpack3.o gf_unpack4.o \
	gf_unpack5.o gf_unpack6.o gf_unpack7.o \
	gettemplates.o getlocal.o getdim.o getpoly.o \
	gribcreate.o gribend.o gribinfo.o mkieee.o \
	rdieee.o simpack.o simunpack.o cmplxpack.o \
	compack.o misspack.o pack_gp.o reduce.o \
	comunpack.o specpack.o specunpack.o jpcpack.o \
	jpcunpack.o enc_jpeg2000.o dec_jpeg2000.o \
	pngpack.o pngunpack.o enc_png.o dec_png.o \
	gbytesc.o skgb.o ixgb2.o getidx.o getg2i.o \
	getg2ir.o getgb2s.o getgb2r.o getgb2l.o \
	getgb2.o getgb2p.o getgb2rp.o putgb2.o g2grids.o \
	gdt2gds.o params.o params_ecmwf.o mova2i.o baciof.o bacio.v1.3.o

$(LIBG2) : $(LIBDIR)/$(LIBG2).a
	ls -l $(LIBDIR)/$(LIBG2).a
$(LIBDIR)/$(LIBG2).a : $(OBJS)
	rm -f $(LIBDIR)/$(LIBG2).a
	ar $(ARFLAG) cr $@ $(OBJS)
clean:
	rm -f $(LIBDIR)/$(LIBG2).a $(LIBDIR)/*.mod $(OBJS) *.mod
	rm -f addfield.f90 gf_unpack7.f90

gf_unpack7.o : gf_unpack7.F
	$(CPP) $(CPPFLAGS) gf_unpack7.F gf_unpack7.f90
	$(FC) -c $(FFLAGS) $(FIXFLAGS) gf_unpack7.f90
addfield.o : addfield.F
	$(CPP) $(CPPFLAGS) addfield.F addfield.f90
	$(FC) -c $(FFLAGS) $(FIXFLAGS) addfield.f90
gribmod.o : gribmod.f
	$(FC) -c $(FFLAGS) $(FIXFLAGS) $<
	rm -f $(TOPDIR)/src/ext2arps/grib_mod.mod
	ln -s `pwd`/grib_mod.mod $(TOPDIR)/src/ext2arps
