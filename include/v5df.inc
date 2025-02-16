! v5df.inc

! Include file for using v5d functions from FORTRAN programs


! Function prototypes.  See the README file for details.  These are
! the functions you'll want to use for writing v5d file converters.

  INTEGER :: v5dcreate

  INTEGER :: v5dcreatesimple

  INTEGER :: v5dwrite

  INTEGER :: v5dmcfile

  INTEGER :: v5dclose


! 5-D grid limits, must match those in v5d.h!!!
  INTEGER :: maxvars, maxtimes, maxrows, maxcolumns, maxlevels

  PARAMETER (maxvars=200)
  PARAMETER (maxtimes=400)
  PARAMETER (maxrows=500)
  PARAMETER (maxcolumns=500)
  PARAMETER (maxlevels=400)

! Missing values
  REAL    :: missing
  INTEGER :: imissing

  PARAMETER (missing=1.0E35)
  PARAMETER (imissing=-987654)

