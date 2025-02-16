#
#	A2IO Makefile --- build liba2io.a
#
#	"Makefile_sun %W%	%G%"
#

#-----------------------------------------------------------
#
# Added by Yunheng to include into ARPS
#
#-----------------------------------------------------------

BZLIB_DIR  = /usr/freeware/include
CFLAGS     =
ARFLAG     = 

TOPDIR = 
LIBDIR = $(TOPDIR)/lib

.c.o:
	$(CC) $(CFLAGS) -c $< 

#------ END of WYH -----------------------------------------

CC=		cc 

AR=		ar

OBJ=		a2io.o \
		compress.o \
		decode_data.o endian.o find_key.o funct_init.o \
		get_data.o get_data_a2io.o get_data_tdwr.o \
		get_record.o \
		radar.o radar_a2io.o radar_tdwr.o \
		radarinfo.o read_radial.o read_record.o \
		read_record_tape.o \
		read_record_raw.o \
		read_record_raw_disk.o \
		read_record_raw_disk_bzip2.o \
		read_record_raw_disk_compress.o \
		read_record_raw_disk_dir.o \
		read_record_raw_disk_file.o \
		read_record_raw_disk_gzip.o \
		read_record_raw_disk_ldm.o \
		read_record_raw_disk_mmap.o \
		read_record_raw_mem.o \
		a2io_fortran.o \
		get_data_fortran.o \
		get_record_fortran.o \
		radar_fortran.o \
		radarinfo_fortran.o \
		read_radial_fortran.o \
		read_record_fortran.o \
		user.o

CFILES=		a2io.c compress.c \
		decode_data.c endian.c find_key.c funct_init.c \
		get_data.c get_data_a2io.c get_data_tdwr.c \
		get_record.c \
		radar.c radar_a2io.c radar_tdwr.c \
		radarinfo.c read_radial.c read_record.c \
		read_record_tape.c \
		read_record_raw.c \
		read_record_raw_disk.c \
		read_record_raw_disk_bzip2.c \
		read_record_raw_disk_compress.c \
		read_record_raw_disk_dir.c \
		read_record_raw_disk_file.c \
		read_record_raw_disk_gzip.c \
		read_record_raw_disk_ldm.c \
		read_record_raw_disk_mmap.c \
		read_record_raw_mem.c \
		a2io_fortran.c \
		get_data_fortran.c \
		get_record_fortran.c \
		radar_fortran.c \
		radarinfo_fortran.c \
		read_radial_fortran.c \
		read_record_fortran.c \
		user.c

#HFILES=	a2defs.h a2io.h const.h defs.h extern.h ncdc.h \
#		read_record.h vcp.h

HFILES=		*.h

MFILES=		Makefile

DFILES=		a2io.doc

DIST=		a2io.tar

LIBA2IO=	liba2io

ALL=		$(LIBA2IO)


all:		$(ALL)

$(LIBA2IO): $(LIBDIR)/$(LIBA2IO).a
	ls -l $(LIBDIR)/$(LIBA2IO).a

$(LIBDIR)/$(LIBA2IO).a:	$(OBJ)
		$(AR) $(ARFLAG) rc $(LIBDIR)/$(LIBA2IO).a $(OBJ)
#		ranlib $(LIBA2IO)

lint:
		lint -ach -DLINT $(CPPFLAGS) $(DEFINES) $(CFILES)

clean:
		/bin/rm -rf $(OBJ) $(LIBDIR)/$(LIBA2IO).a $(DIST) core

dist:
		tar cf $(DIST) $(CFILES) $(HFILES) $(MFILES) $(DFILES)

depend:
		makedepend $(CPPFLAGS) $(DEFINES) $(CFILES)

# DO NOT DELETE THIS LINE -- make depend depends on it.

a2io.o: a2io.h a2defs.h radardefs.h extern.h proto.h
compress.o: config.h proto.h
decode_data.o: a2io.h a2defs.h radardefs.h a2io_message_31.h const.h extern.h
decode_data.o: vcp.h
endian.o: config.h a2io.h a2defs.h radardefs.h a2io_message_31.h extern.h
endian.o: proto.h tdwrio.h tdwrdefs.h tdwr_extern.h
funct_init.o: config.h funct_init.h proto.h
get_data.o: a2io.h a2defs.h radardefs.h extern.h extern_funct.h config.h
get_data.o: proto.h
get_data_a2io.o: a2io.h a2defs.h radardefs.h a2io_message_31.h extern.h
get_data_a2io.o: extern_funct.h config.h proto.h
get_data_tdwr.o: config.h a2defs.h radardefs.h extern.h a2io.h tdwrio.h
get_data_tdwr.o: tdwrdefs.h tdwr_extern.h proto.h
get_record.o: a2io.h a2defs.h radardefs.h extern.h proto.h
index.o: a2defs.h radardefs.h
radar.o: a2defs.h radardefs.h const.h extern.h extern_funct.h config.h
radar.o: proto.h
radar_a2io.o: config.h a2io.h a2defs.h radardefs.h a2io_message_31.h const.h
radar_a2io.o: extern.h proto.h
radar_tdwr.o: config.h a2defs.h radardefs.h const.h extern.h proto.h tdwrio.h
radar_tdwr.o: tdwrdefs.h tdwr_extern.h tdwr_nyq.h
radarinfo.o: a2defs.h radardefs.h proto.h
read_radial.o: config.h a2defs.h radardefs.h tdwrdefs.h proto.h
read_record.o: config.h a2io.h a2defs.h radardefs.h a2io_message_31.h const.h
read_record.o: extern.h ncdc.h read_record.h tdwrio.h tdwrdefs.h
read_record.o: tdwr_struct.h vcp.h proto.h
read_record_raw.o: a2io.h a2defs.h radardefs.h extern.h proto.h
read_record_raw_disk_bzip2.o: config.h proto.h a2io.h a2defs.h radardefs.h
read_record_raw_disk_bzip2.o: tdwrio.h tdwrdefs.h tdwr_extern.h
read_record_raw_disk_compress.o: config.h proto.h a2io.h a2defs.h radardefs.h
read_record_raw_disk_compress.o: tdwrio.h tdwrdefs.h tdwr_extern.h
read_record_raw_disk_dir.o: a2io.h a2defs.h radardefs.h config.h extern.h
read_record_raw_disk_file.o: a2io.h a2defs.h radardefs.h a2io_message_31.h
read_record_raw_disk_file.o: extern.h ncdc.h config.h tdwrio.h tdwrdefs.h
read_record_raw_disk_file.o: tdwr_extern.h proto.h
read_record_raw_disk_gzip.o: config.h proto.h
read_record_raw_disk_ldm.o: a2io.h a2defs.h radardefs.h config.h proto.h
read_record_raw_disk_mmap.o: config.h proto.h
read_record_raw_mem.o: config.h proto.h mem.h
read_record_tape.o: config.h a2io.h a2defs.h radardefs.h a2io_message_31.h
read_record_tape.o: const.h extern.h ncdc.h tdwrio.h tdwrdefs.h tdwr_extern.h
read_record_tape.o: vcp.h proto.h
user.o: extern.h
