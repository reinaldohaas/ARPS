
#-----------------------------------------------------------
#
# For including into the ARPS package
#
#-----------------------------------------------------------

CC     = gcc
CFLAGS = -g -I.
AR     = ar
ARFLAG =

TOPDIR =
LIBDIR = $(TOPDIR)/lib

LIBSOLOR = libsoloread
LIBSOLOD = libsolodump

#-----------------------------------------------------------------------
#
# Default rules
#
#-----------------------------------------------------------------------

.c.o:
	$(CC) $(CFLAGS) -c $<

#-----------------------------------------------------------------------
#
# Objects
#
#-----------------------------------------------------------------------

INC = alloc.h  dorade.h  errMsg.h  geoLines.h  geography.h  \
      getVals.h  hash.h  radar.h  sizes.h

OBJ = dorade_lib.o errMsg.o hash.o geography.o geoLines.o  \
      getVals.o geoTime.o radar.o mydoradelib.o

A2OBJ = a2ioapi.o a2ioapi_f.o dumpsweeps.o dumpsweeps_f.o $(OBJ)

DPOBJ = dumpsweeps.o dumpsweeps_f.o $(OBJ)

#-----------------------------------------------------------------------
#
# Rules
#
#-----------------------------------------------------------------------

$(LIBSOLOR): $(LIBDIR)/$(LIBSOLOR).a
	ls -l $(LIBDIR)/$(LIBSOLOR).a

$(LIBDIR)/$(LIBSOLOR).a: $(A2OBJ)
	$(AR) $(ARFLAG) rc $(LIBDIR)/$(LIBSOLOR).a $(A2OBJ)
#       ranlib $(LIBA2IO)

$(LIBSOLOD): $(LIBDIR)/$(LIBSOLOD).a
	ls -l $(LIBDIR)/$(LIBSOLOD).a

$(LIBDIR)/$(LIBSOLOD).a: $(DPOBJ)
	$(AR) $(ARFLAG) rc $(LIBDIR)/$(LIBSOLOD).a $(DPOBJ)
#       ranlib $(LIBA2IO)

clean :
	rm -f $(LIBDIR)/$(LIBSOLOR).a $(LIBDIR)/$(LIBSOLOD).a core *.o

#-----------------------------------------------------------------------
#
# Dependency
#
#-----------------------------------------------------------------------

dorade_lib.o : dorade_lib.c dorade.h
	${CC} ${CFLAGS} -c dorade_lib.c

