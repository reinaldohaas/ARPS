
!########################################################################
!########################################################################
!######                                                            ######
!######                  SUBROUTINE QPFMCORN                       ######
!######                                                            ######
!######                      Developed by                          ######
!######         Center for Analysis and Prediction of Storms       ######
!######                   University of Oklahoma                   ######
!######                                                            ######
!########################################################################
!########################################################################

SUBROUTINE QPFMCORN(intmapproj,intscale,inttrulat1,inttrulat2, &
                    inttrulon,intdx,intdy,nx,ny,lat2d,lon2d, &
                    imin,imax,jmin,jmax, &
                    mapswlat,mapswlon,mapnwlat,mapnwlon, &
                    mapselat,mapselon,mapnelat,mapnelon)

!-----------------------------------------------------------------------
!
! PURPOSE:  
!
! Calculates latitude/longitude coordinates of four corners of scalar
! model grid for NCL map plotting.  ARPS map projection software
! is invoked in this subroutine.  Designed for use with INTQPF.
!
! AUTHOR:
!
! Eric Kemp, March 2000
!
!----------------------------------------------------------------------- 
!
! Variable declarations
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: nx,ny,imin,imax,jmin,jmax
  INTEGER, INTENT(IN) :: intmapproj
  REAL, INTENT(IN) :: intscale,inttrulat1,inttrulat2,inttrulon
  REAL, INTENT(IN) :: intdx,intdy
  REAL, INTENT(IN) :: lat2d(nx,ny),lon2d(nx,ny)
  REAL, INTENT(OUT) :: mapswlat,mapswlon,mapnwlat,mapnwlon, &
                       mapselat,mapselon,mapnelat,mapnelon

!----------------------------------------------------------------------- 
!
! Miscellaneous Variables
!
!-----------------------------------------------------------------------

  REAL :: inttrulat(2)
  REAL :: x_ext(nx),y_ext(ny)
  INTEGER :: i,j
  REAL :: fx0,fy0

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  inttrulat(1) = inttrulat1
  inttrulat(2) = inttrulat2
  
  CALL setmapr(intmapproj,intscale,inttrulat,inttrulon)
  CALL lltoxy(1,1,lat2d(1,1),lon2d(1,1),fx0,fy0)
  CALL setorig(1,fx0,fy0)

  DO i = 1,nx
    x_ext(i) = (i-1)*intdx
  END DO

  DO j = 1,ny
    y_ext(j) = (j-1)*intdy
  END DO

  CALL mcorner(nx,ny,x_ext,y_ext,imin,imax,jmin,jmax,&
               intdx,intdy,mapswlat,mapswlon, &
               mapnwlat,mapnwlon,mapselat,mapselon,mapnelat, &
               mapnelon)

END SUBROUTINE QPFMCORN

!########################################################################
!########################################################################
!######                                                            ######
!######                   SUBROUTINE MCORNER                       ######
!######                                                            ######
!######                      Developed by                          ######
!######         Center for Analysis and Prediction of Storms       ######
!######                   University of Oklahoma                   ######
!######                                                            ######
!########################################################################
!########################################################################

SUBROUTINE MCORNER(nx,ny,x_ext,y_ext,imin,imax,jmin,jmax, &
                   dx,dy,mapswlat,mapswlon, &
                   mapnwlat,mapnwlon,mapselat,mapselon, &
                   mapnelat,mapnelon)

!-----------------------------------------------------------------------
!
! PURPOSE:  
!
! Calculates latitude/longitude coordinates of four corners of scalar
! model grid for NCL map plotting.  Note that the map projection must
! have already been set before calling this subroutine.
!
! AUTHOR:
!
! Eric Kemp, March 2000
!
!----------------------------------------------------------------------- 
!
! Variable declarations
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: nx,ny,imin,imax,jmin,jmax
  REAL, INTENT(IN) :: x_ext(nx),y_ext(ny)
  REAL, INTENT(IN) :: dx,dy

  REAL, INTENT(OUT) :: mapswlat,mapswlon,mapnwlat,mapnwlon, &
                       mapselat,mapselon,mapnelat,mapnelon

!----------------------------------------------------------------------- 
!
! Miscellaneous variables
!
!-----------------------------------------------------------------------

  REAL :: x0,y0

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  x0 = x_ext(imin) - (dx/2.)
  y0 = y_ext(jmin) - (dy/2.)
  CALL xytoll (1,1,x0,y0,mapswlat,mapswlon)

  x0 = x_ext(imin) - (dx/2.)
  y0 = y_ext(jmax) + (dy/2.)
  CALL xytoll (1,1,x0,y0,mapnwlat,mapnwlon)

  x0 = x_ext(imax) + (dx/2.)
  y0 = y_ext(jmin) - (dy/2.)
  CALL xytoll (1,1,x0,y0,mapselat,mapselon)

  x0 = x_ext(imax) + (dx/2.)
  y0 = y_ext(jmax) + (dy/2.)
  CALL xytoll (1,1,x0,y0,mapnelat,mapnelon)

END SUBROUTINE MCORNER
