!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE PLOT                       ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE plot(ntx,nty,nptx,iproj,alatpro,alonpro,cint,                &
           latul,lonul,latlr,lonlr, glab,h)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Plot an array given a projection type (iproj), the lat/lon of the
!  corner points and the graphics label (glab) using NCAR graphics
!  subroutines.
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!
!  AUTHOR:  Dan Weber
!  1/12/94
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    ntx      Dimension in the e-w direction for the plot array
!    nty      Dimension in the n-s direction for the plot array
!    nptx     Limiting dimension in the e-w direction for the plot array
!    iproj    Projection type for use in supmap
!    alatpro  Projection latitudes (true at)
!    alonpro  Projection longitude (true at)
!    cint     Contour interval to be passed to conrec.
!
!    latul    Latitude of the upper left corner of the area
!    lonul    Longitude of the upper left corner of the area
!    latlr    Latitude of the lower right corner of the area
!    lonlr    Longitude of the lower right corner of the area
!    glab     Label for the difference field plot
!
!    h        Variable to be contoured
!
!  COMMON block variables:
!
!    alatpro  Projection latitudes (true at)
!    alonpro  Projection longitude (true at)
!
!  SUPMAP variables:
!
!    rota     Rotation from the standard longitude
!    jproj    Projection type in supmap
!    plat     Standard Latitude of the projection
!    plon     Standard Longitude of the projection
!    jlts     Plot corners to be specified
!    jgrd     Control parameter for plotting lat/lon lines
!    iout     Control parameter for drawing US state outlines
!    idot     Control parameter for line type drawn
!    ierr
!
!  OUTPUT:
!
!             Plot of h using supmap and conrec
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!-----------------------------------------------------------------------
!

  INTEGER :: ntx               ! Number of points to be plotted in the
                               ! e-w direction
  INTEGER :: nty               ! Number of points to be plotted in the
                               ! n-s direction
  INTEGER :: nptx              ! i dimension of the h array
  INTEGER :: iproj             ! Projection type defined by user
  REAL :: alatpro(2)           ! Projection latitudes (true at)
  REAL :: alonpro              ! Projection longitude (true at)
  REAL :: cint                 ! Contour interval to be passed to conrec

  REAL :: latul                ! Latitude of the upper left corner
  REAL :: lonul                ! Longitude of the upper left corner
  REAL :: latlr                ! Latitude of the lower right corner
  REAL :: lonlr                ! Longitude of the lower right corner
  CHARACTER (LEN=65) :: glab         ! Label for the difference field plot

  REAL :: h(nptx,nty)          ! Variable to be contoured
!
!-----------------------------------------------------------------------
!
!  OUTPUT:
!
!    Plot of h using supmap and conrec
!
!-----------------------------------------------------------------------
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  WRITE (6, '(a)') 'No NCAR Graphics plotting.'

  RETURN
END SUBROUTINE plot
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE PLOTINT                    ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE plotint
  RETURN
END SUBROUTINE plotint
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE opngks                     ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE opngks
  RETURN
END SUBROUTINE opngks
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE ezxy                       ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE ezxy
  RETURN
END SUBROUTINE ezxy
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE agseti                     ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE agseti
  RETURN
END SUBROUTINE agseti
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE agsetc                     ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE agsetc
  RETURN
END SUBROUTINE agsetc
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE clsgks                     ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE clsgks
  RETURN
END SUBROUTINE clsgks
