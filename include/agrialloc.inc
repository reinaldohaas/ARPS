!
!  lstore is the size of the main storage array.  It can be
!  changed iat any time if needed.  statistics on its use
!  are output on restart and exiting of the interface.
!  ihighwater is the maximum amount of a(*) used by the interface,
!  and ntemp denotes the boundary between grid storage and
!  free (temp) space that is available for use by any routine
!
  INTEGER :: lstore,ntemp,ihighwater

  PARAMETER( lstore = 99900000 ) ! for 1 67x67x35 grids

  REAL :: a(lstore)
  COMMON /cmalloc/ ntemp,a,ihighwater
