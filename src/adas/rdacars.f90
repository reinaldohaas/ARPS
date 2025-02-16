!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINE RDACARS                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE rdacars(nvar,mxsng,acarsfile,                                &
           stnsng,latsng,lonsng,xsng,ysng,hgtsng,                       &
           obsng,qualsng,isrcsng,                                       &
           rmiss,nprev,ntotal,istatus)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Read ACARS data.
!  This is a dummy subroutine provided as an example of
!  reading non-surface single-level data.  The data is appended
!  to surface data already collected (nprev).
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster, CAPS
!  October, 1997
!
!  MODIFICATION HISTORY:
!
!
!-----------------------------------------------------------------------
!
!  Variable Declarations:
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
  INTEGER :: nvar,mxsng

  CHARACTER (LEN=*) :: acarsfile
  CHARACTER (LEN=5) :: stnsng(mxsng)
  REAL :: latsng(mxsng)
  REAL :: lonsng(mxsng)
  REAL :: xsng(mxsng)
  REAL :: ysng(mxsng)
  REAL :: hgtsng(mxsng)
  REAL :: obsng(nvar,mxsng)
  INTEGER :: qualsng(nvar,mxsng)
  INTEGER :: isrcsng(mxsng)
  REAL :: rmiss
  INTEGER :: nprev
  INTEGER :: ntotal
  INTEGER :: istatus
!
!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  ntotal=nprev
  istatus=0
  RETURN
END SUBROUTINE rdacars
