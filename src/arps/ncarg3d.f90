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
                latul,lonul,latlr,lonlr,glab,var)
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
!  MODIFICATIONS:
!
!  K. Brewster 1/30/95
!  Calculates appropriate jgrd.
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
!    var      Variable to be contoured
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
!
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

  REAL :: var(nptx,nty)        ! Variable to be contoured
!
!-----------------------------------------------------------------------
!
!  SUPMAP variables:
!
!-----------------------------------------------------------------------
!
  REAL :: rota                 ! rotation from the standard longitude
  INTEGER :: jproj             ! projection type in supmap
  INTEGER :: jlts,jgrd,iout,idot
  REAL :: plat,plon            ! standard latitudes of the projection
  PARAMETER (jlts = 2,      & ! corners are specified
         iout = 4,      & ! draw all outlines
         idot = 1 )     ! continuous lines -- idot=1 dotted
  REAL :: plm1,plm2,plm3,plm4
!
!-----------------------------------------------------------------------
!
!  Local miscellanenous variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: ierr
  REAL :: degr
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  plm1=latul
  plm2=lonul
  plm3=latlr
  plm4=lonlr
  ierr=100

  degr=AMAX1(ABS(plm2-plm4),ABS(plm1-plm3))
  IF(degr < 10.) THEN
    jgrd=1
  ELSE IF(degr < 20.) THEN
    jgrd=2
  ELSE IF(degr < 50.) THEN
    jgrd=5
  ELSE
    jgrd=10
  END IF

  PRINT *,'entering plot'
  IF(iproj == 1) THEN       ! stereographic projection...
    jproj=iproj
    plat=90.
    plon=0.
    rota=alonpro
  ELSE IF(iproj == 8) THEN  ! cylindrical equidistant projection...
    jproj=iproj
    plat=0.
    plon=0.
    rota=0.
  ELSE IF(iproj == 2) THEN  ! lambert conformal projection...
    jproj=3
    plat=alatpro(1)
    plon=alonpro
    rota=alatpro(2)
  ELSE IF(iproj == 3) THEN  ! mercator conformal projection...
    jproj=9
    plat=0.
    plon=0.
    rota=0.
  ELSE                      ! Projection type not supported by program...
    PRINT *,'projection type not supported by the program'
    PRINT *,'program will print iproject and stop'
    PRINT *,'iproject=',iproj
    CALL arpsstop('arpsstop called from PLOT improper projection ',1)

  END IF

  PRINT *,'corners of the plot grid are',plm1,plm2,plm3,plm4
  CALL set(0.,1.,0.,1.,0.1,0.9,0.1,0.9,1)
  CALL wtstr(.5,.88,glab,15,0,0)
  CALL supmap(jproj,plat,plon,rota,plm1,plm2,                           &
              plm3,plm4,jlts,jgrd,iout,idot,ierr)

  PRINT *,'entering conrec'
  CALL conrec(var,nptx,ntx,nty,0.0,0.0,cint,-1,0,0)

  CALL frame

  PRINT *,'leaving plot'

  RETURN
END SUBROUTINE plot
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE PLOTINT                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE plotint(ntx,nty,nptx,iproj,alatpro,alonpro,                  &
           latul,lonul,latlr,lonlr,glab,ivar)
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
!  MODIFICATIONS:
!
!  K. Brewster 1/30/95
!  Calculates appropriate jgrd.
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
!
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

  REAL :: latul                ! Latitude of the upper left corner
  REAL :: lonul                ! Longitude of the upper left corner
  REAL :: latlr                ! Latitude of the lower right corner
  REAL :: lonlr                ! Longitude of the lower right corner
  CHARACTER (LEN=65) :: glab         ! Label for the difference field plot

  INTEGER :: ivar(nptx,nty)    ! Integer variable to be contoured
!
!-----------------------------------------------------------------------
!
!  SUPMAP variables:
!
!-----------------------------------------------------------------------
!
  REAL :: rota                 ! rotation from the standard longitude
  INTEGER :: jproj             ! projection type in supmap
  INTEGER :: jlts,jgrd,iout,idot
  REAL :: plat,plon            ! standard latitudes of the projection
  PARAMETER (jlts = 2,      & ! corners are specified
         iout = 4,      & ! draw all outlines
         idot = 1 )     ! continuous lines -- idot=1 dotted
  REAL :: plm1,plm2,plm3,plm4
!
!-----------------------------------------------------------------------
!
!  Local miscellanenous variables
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=6) :: FMT
  CHARACTER (LEN=8) :: ch
  INTEGER :: i,j
  REAL :: xx1,xx2,yy1,yy2,zz1,zz2,ww1,ww2
  INTEGER :: ll
  INTEGER :: ierr
  INTEGER :: ifmt

  REAL :: degr

  REAL :: rgb(3,14)
  INTEGER :: ioc(14)
  DATA ioc / 6,2,5,12,10,11,1,3,4,8,9,7,13,14 /
  DATA rgb / 0.00 , 0.00 , 0.00 ,                                       &
             0.75 , 0.50 , 1.00 ,                                       &
             0.50 , 0.00 , 1.00 ,                                       &
             0.00 , 0.00 , 1.00 ,                                       &
             0.00 , 1.00 , 0.00 ,                                       &
             0.25 , 0.25 , 0.25 ,                                       &
             0.00 , 1.00 , 0.60 ,                                       &
             0.00 , 0.50 , 1.00 ,                                       &
             0.70 , 1.00 , 0.00 ,                                       &
             0.00 , 0.00 , 1.00 ,                                       &
             1.00 , 0.75 , 0.00 ,                                       &
             1.00 , 0.00 , 0.00 ,                                       &
             1.00 , 0.00 , 0.75 ,                                       &
             0.00 , 0.80 , 0.10 /
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  plm1=latul
  plm2=lonul
  plm3=latlr
  plm4=lonlr
  ierr=100

  degr=AMAX1(ABS(plm2-plm4),ABS(plm1-plm3))
  IF(degr < 10.) THEN
    jgrd=1
  ELSE IF(degr < 20.) THEN
    jgrd=2
  ELSE IF(degr < 50.) THEN
    jgrd=5
  ELSE
    jgrd=10
  END IF

  PRINT *,'entering plot'
  IF(iproj == 1) THEN       ! stereographic projection...
    jproj=iproj
    plat=90.
    plon=0.
    rota=alonpro
  ELSE IF(iproj == 8) THEN  ! cylindrical equidistant projection...
    jproj=iproj
    plat=0.
    plon=0.
    rota=0.
  ELSE IF(iproj == 2) THEN  ! lambert conformal projection...
    jproj=3
    plat=alatpro(1)
    plon=alonpro
    rota=alatpro(2)
  ELSE IF(iproj == 3) THEN  ! mercator conformal projection...
    jproj=9
    plat=0.
    plon=0.
    rota=0.
  ELSE                      ! Projection type not supported by program...
    PRINT *,'projection type not supported by the program'
    PRINT *,'program will print iproject and stop'
    PRINT *,'iproject=',iproj
    CALL arpsstop('arpsstop called from PLOTINT improper projection ',1)

  END IF

  DO j=1,14
    i=ioc(j)
    CALL gscr(1,j,rgb(1,i),rgb(2,i),rgb(3,i))
  END DO

  CALL set(0.,1.,0.,1.,0.1,0.9,0.1,0.9,1)
  CALL wtstr(.5,.88,glab,15,0,0)
  CALL supmap(jproj,plat,plon,rota,plm1,plm2,                           &
              plm3,plm4,jlts,jgrd,iout,idot,ierr)

  PRINT *,'corners of the plot grid are',plm1,plm2,plm3,plm4

  CALL getset(xx1,xx2,yy1,yy2,zz1,zz2,ww1,ww2,ll)
  CALL set(xx1,xx2,yy1,yy2,1.0,FLOAT(nptx+1),1.0,FLOAT(nty),1)

  DO j=1,nty
    DO i=1,nptx
      ifmt = nint( LOG10(FLOAT(ivar(i,j))) ) + 1
      WRITE(FMT,'(a,i1,a)') '(i',ifmt,')'
      WRITE(ch,FMT) ivar(i,j)
      CALL gsplci(ivar(i,j))
      CALL pwritx(FLOAT(i),FLOAT(j),ch(1:ifmt),2,8,0,0)
    END DO
  END DO

  CALL frame

  PRINT *,'leaving plot'

  RETURN
END SUBROUTINE plotint
