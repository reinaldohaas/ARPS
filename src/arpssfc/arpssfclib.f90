!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE GET_GED                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
!

SUBROUTINE get_ged(dtype,filenm,nx,ny,xdims,ydims,resl,glon,glat,  &
                   iout,rout)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Read in the surface data set and determine the region just
!  covering the model domain.
!
!  For soil data, the resolution is   1 x 1   degrees of lat x lon
!  For NDVI data, the resolution is 1/6 x 1/6 degrees of lat x lon
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yuhe Liu
!
!  1/28/94
!
!  MODIFICATIONS:
!
!  2001/07/26 Gene Bassett
!  Moved calls to getstyp1, getvtyp1, getndvi1 inside this subroutine.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!  dtype    Data flag: 1 for soil type data
!                      2 for vegetation type data
!  filenm   File name of usrface data set
!  nx       Number of grid points in the x-direction
!  ny       Number of grid points in the y-direction
!  xdims    X-dimsion size of the data set
!  ydims    Y-dimsion size of the data set
!
!  glon     Longitude values of model grid points
!  glat     Latitude  values of model grid points
!
!  OUTPUT:
!
!  dtlat    Latitude  values of data grid points
!  dtlon    Longitude values of data grid points
!
!  dtarr    Integer data array
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: dtype        ! Data flag: 1 for soil; 2 for ndvi
  CHARACTER (LEN=*  ) :: filenm ! File name of usrface data set
  INTEGER :: nx           ! Number of grid points in the x-direction
  INTEGER :: ny           ! Number of grid points in the y-direction
  INTEGER :: xdims        ! Columns of/ the data set
  INTEGER :: ydims        ! Rows of the data set

  REAL :: resl         ! Resolution of data
  INTEGER :: dtproj       ! Projection of data

  REAL :: glon(nx,ny)  ! Longitude values of model grid points
  REAL :: glat(nx,ny)  ! Latitude  values of model grid points

  INTEGER :: iout(1)   ! Points to vegtyp with dtype=2
  REAL    :: rout(1)   ! Points to soiltyp when dtype=1, ndvi when dtype=3
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  REAL,    ALLOCATABLE :: dtlat(:) ! Latitude  values of data grid points
  REAL,    ALLOCATABLE :: dtlon(:) ! Longitude values of data grid points
  INTEGER, ALLOCATABLE :: dtarr(:,:)  ! Data array

  INTEGER :: lenfl        ! Length of filenm
  INTEGER :: flunit       ! IO unit of filenm
  CHARACTER (LEN=20) :: title     ! Title of data
  CHARACTER (LEN=20) :: temtit    ! Temporary string of data title
  INTEGER :: colmn        ! Columns of data
  INTEGER :: row          ! Rows of data
  REAL    :: resl0        ! Resolution of data

  INTEGER :: i,j
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'arpssfc.inc'
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  ALLOCATE (dtlat(ydims))
  ALLOCATE (dtlon(xdims))
  ALLOCATE (dtarr(xdims,ydims))

  lenfl = LEN( filenm )
  CALL strlnth( filenm, lenfl )
  CALL getunit( flunit )

  WRITE(6,'(a,a)') '    Opening ',filenm(1:lenfl)
  OPEN (UNIT = flunit, FILE = filenm(1:lenfl), STATUS = 'old')
!    :      form = 'unformatted', access = 'sequential')

  READ (flunit,'(a20,i8,i8,e20.10,i8)')                                 &
       title, colmn, row, resl, dtproj

  IF ( dtype == 1 ) THEN
    temtit = 'W & H-S Soil Type   '
    resl0 = 1.0
  ELSE IF ( dtype == 2 )THEN
    temtit = 'World Ecosystems    '
    resl0 = 1.0/6.0
  ELSE IF ( dtype == 3 )THEN
    temtit = 'NDVI data           '
    resl0 = 1.0/6.0
  ELSE
    WRITE (6, '(a/a)')                                                  &
        'Data type is not correct.',                                    &
        'Program stops here in subroutine GET_GED.'
    STOP
  END IF

  IF ( title /= temtit ) THEN
    WRITE (6, '(a/a/a/a/a)')                                            &
        'Data title is not correct.',                                   &
        'The file opened perhaps was not what you wanted.',             &
        'title required = '//temtit,                                    &
        'title in the file = '//title,                                  &
        'Program stops here in subroutine GET_GED.'
    STOP
  ELSE IF ( colmn /= xdims .OR. row /= ydims ) THEN
    WRITE (6, '(a/a/a,i5,a,i5/a,i5,a,i5/a)')                            &
        'Data dimension size is not correct.',                          &
        'The file opened perhaps was not what you wanted.',             &
        'column required = ',xdims, ', column in the file = ', colmn,   &
        'row required = ',ydims, ', row in the file = ', row,           &
        'Program stops here in subroutine GET_GED.'
    STOP
  ELSE IF ( ABS( ( resl - resl0 ) / resl ) > 0.01 ) THEN
    WRITE (6, '(a/a/a,i5/a,i5/a)')                                      &
        'Data resolution is not correct.',                              &
        'The file opened perhaps was not what you wanted.',             &
        'resolution required = ',resl0,                                 &
        'resolution in the file = ', resl,                              &
        'Program stops here in subroutine GET_GED.'
    STOP
  ELSE IF ( dtproj /= 1 ) THEN
    WRITE (6, '(a/a/a,i5/a,i5/a)')                                      &
        'Data projection is not correct.',                              &
        'The file opened perhaps was not what you wanted.',             &
        'projection required = ', 1,                                    &
        'projection in the file = ', dtproj,                            &
        'Program stops here in subroutine GET_GED.'
    STOP
  END IF

  WRITE(6,'(a)') '    Reading data set for GED '//title
  READ (flunit,'(20i4)') dtarr

  DO i = 1, xdims
    dtlon(i) = resl*(FLOAT(i)-.5)
  END DO

  DO j = 1, ydims
    dtlat(j) = resl*(FLOAT(j)-.5) - 90.0
  END DO

  CLOSE ( flunit )
  CALL retunit ( flunit )

  IF (dtype == 1) THEN

!-----------------------------------------------------------------------
!  Translate the soil classes data (stypdat) into the USDA soil type
!  category and fill them into the model grid array soiltyp (in iout).
!-----------------------------------------------------------------------

    CALL getstyp2(xdims,ydims,nx,ny,dtlon,dtlat,resl,glon,glat,dtarr,iout)

  ELSE IF (dtype == 2) THEN

!-----------------------------------------------------------------------
!  Transfer the World Ecosystem Classes data into the 12 vegetation
!  type categories and fill into model grid array, vegtyp(nx,ny) (in iout)
!-----------------------------------------------------------------------

    CALL getvtyp2(xdims,ydims,nx,ny,dtlon,dtlat,resl,glon,glat,dtarr,iout)

  ELSE ! dtype = 3

!-----------------------------------------------------------------------
!  Calculate the real-NDVI data from byte-NDVI and fill into model
!  grid array, ndvi(nx,ny) (in rout)
!-----------------------------------------------------------------------

    CALL getndvi2(xdims,ydims,nx,ny,dtlon,dtlat,resl,glon,glat,dtarr,rout)

  ENDIF

  DEALLOCATE (dtlat)
  DEALLOCATE (dtlon)
  DEALLOCATE (dtarr)

  RETURN
END SUBROUTINE get_ged
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE GET_STATSGO                ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
!
SUBROUTINE get_statsgo(flag,filenm1,filenm2,nx,ny,nstyps,x1,y1,        &
           col,row,glat,glon,resl,  &
           iout,rout)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Read in the surface property data set & determine the region just
!  covering the model domain.
!
!  the resolution is   1 * 1 km  of grid
!
!-----------------------------------------------------------------------
!
!  AUTHOR:
!
!  Leilei Wang ,Vince Wong
!  3/18/97
!
!  MODIFICATIONS:
!
!  2001/07/26 Gene Bassett
!  Moved calls to getstyp1, getvtyp1, getndvi1 inside this subroutine.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!  filenm   File name of usrface data set
!
!  OUTPUT:
!
!  iout     Integer data array (soiltyp when flag=1, vegtyp when flag=2)
!  rout     Real data array (stypfrct when flag=1, ndvi when flag=3)
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
  INTEGER :: flag         ! Flag for what kind of data
  CHARACTER (LEN=* ) :: filenm1 ! File name of usrface data set
  CHARACTER (LEN=* ) :: filenm2 ! File name of GED     data set
  INTEGER :: nx
  INTEGER :: ny
  INTEGER :: nstyps

  INTEGER :: dcol,drow
  INTEGER :: col,row      ! Columns of the data set

  INTEGER :: xdims        ! Actual colmns in data array
  INTEGER :: ydims        ! Actual rows in data array
!
  REAL :: resl              ! Resolution of data
  REAL :: resl_ged          ! Resolution of GED soil data

  REAL :: dctrlat           ! LAT. of projection center of data
  REAL :: dctrlon           ! LON. of projection center of data

  REAL :: x1(nx,ny)         ! X_coordinate of arps model grid
  REAL :: y1(nx,ny)         ! Y_coordinate of arps model grid
  REAL :: glat(nx,ny)       ! latitude   of arps model grid
  REAL :: glon(nx,ny)       ! longitude   of arps model grid
  REAL, ALLOCATABLE :: xdat(:)  ! X_coordinate of data grid covering model
                                ! domain
  REAL, ALLOCATABLE :: ydat(:)  ! Y_coordinate of data grid covering model
                                ! domain

  INTEGER, ALLOCATABLE :: idtarr (:,:)       ! Output data array

  INTEGER :: iout(1)   ! Points to soiltyp when flag=1, vegtyp when flag=2
  REAL :: rout(1)      ! Points to stypfrct when flag=1, ndvi when flag=3
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: lenfl
  INTEGER :: flunit
  INTEGER :: xmin,ymin
  INTEGER :: xmax,ymax
  REAL :: dnwx,dnwy,dsex,dsey
  REAL :: xxmin,yymax
  REAL :: xorgmin,xorgmax
  REAL :: yorgmin,yorgmax
  INTEGER :: i,j,m,n

  INTEGER, ALLOCATABLE :: item1(:,:) ! add by wyh
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'arpssfc.inc'
  INCLUDE 'globcst.inc'
  INCLUDE 'grid.inc'
  INCLUDE 'phycst.inc'
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

  lenfl = 70
  CALL strlnth( filenm1, lenfl )
  CALL getunit( flunit )

  resl = 1000.
  dctrlon = -100.00

  IF (flag == 1) THEN
    dcol = 4587
    drow = 2889
    dctrlat =    45.0
    dnwx =  - 2050500
    dnwy =    752500
    dsex =    2536500
    dsey =  - 2136500
  ELSE
    dcol = 9223
    drow = 8996
    dctrlat =    50.0
    dnwx =  - 4486550
    dnwy =    4479868
    dsex =    4735450
    dsey =  - 4515132
  END IF

  DO i = 1,nx
    DO j = 1,ny
      CALL lae_lltoxy(1,1,dctrlat,dctrlon,glat(i,j),glon(i,j),          &
                      x1(i,j),y1(i,j))
    END DO
  END DO
!
  xorgmin = x1(1, 1)
  xorgmax = x1(nx,1)
  yorgmin = y1(1, 1)
  yorgmax = y1(1,ny)

  DO i = 1, nx
    DO j = 1, ny
      IF(x1(i,j) < xorgmin) THEN
        xorgmin = x1(i,j)
      END IF
      IF(x1(i,j) > xorgmax) THEN
        xorgmax = x1(i,j)
      END IF
      IF(y1(i,j) < yorgmin) THEN
        yorgmin = y1(i,j)
      END IF
      IF(y1(i,j) > yorgmax) THEN
        yorgmax = y1(i,j)
      END IF
    END DO
  END DO
!
  IF(dx < resl) THEN
    xmin = INT( (xorgmin-dnwx-2*resl)/resl )
    xmax = INT( (xorgmax-dnwx+2*resl)/resl )+1
  ELSE
    xmin = INT( (xorgmin-dnwx-2*dx)/resl )
    xmax = INT( (xorgmax-dnwx+2*dx)/resl )+1
  END IF

  IF(dy < resl) THEN
    ymin = INT( (dnwy-yorgmax-2*resl)/resl )
    ymax = INT( (dnwy-yorgmin+2*resl)/resl )+1
  ELSE
    ymin = INT( (dnwy-yorgmax-2*dy)/resl )
    ymax = INT( (dnwy-yorgmin+2*dy)/resl )+1
  END IF

  CALL get_ged(flag,filenm2,nx,ny,col,row,resl_ged,glon,glat,iout,rout)

  IF ( xmin < 1   ) xmin = 1
  IF ( xmax > dcol) xmax = dcol
  IF ( ymin < 1   ) ymin = 1
  IF ( ymax > drow) ymax = drow

  xdims = xmax - xmin + 1
  ydims = ymax - ymin + 1

  ALLOCATE(xdat(xdims))
  ALLOCATE(ydat(ydims))
  ALLOCATE(idtarr(xdims,ydims))
  ALLOCATE(item1(xdims, ydims))

  xxmin = dnwx + (xmin-1) * resl
  yymax = dnwy - (ymax-1) * resl

  item1 = 0

  CALL readbsq (filenm1, flunit, 1, drow,dcol, 1,                       &
      1, ymin, xmin, 1, ydims , xdims,item1,xdims,ydims)

  DO j=1,ydims
    DO i=1,xdims
      idtarr(i,j) = item1(i,ydims-j+1)
    END DO
  END DO

  DEALLOCATE(item1)

  DO i = 1 ,xdims
    xdat( i ) = xxmin + (i-1)*resl
  END DO

  DO j = 1 ,ydims
    ydat( j ) = yymax + (j-1)*resl
  END DO

  CLOSE ( flunit )
  CALL retunit ( flunit )

  IF (flag == 1) THEN

!-----------------------------------------------------------------------
!  Translate the soil classes data storged in stypdat into the USDA soil
!  type categories. nstyp soil types with the highest percentages are
!  then brought to the model grid points, with each type being associated
!  with a value of soil type fraction. (iout is soiltyp, rout is stypfrct)
!-----------------------------------------------------------------------

    CALL getstyp1(xdims,ydims, nx,ny,nstyps,glon,glat,              &
                  x1,y1,idtarr,xdat,ydat,iout,rout)

  ELSE IF (flag == 2) THEN

!-----------------------------------------------------------------------
!  Calculate vegetation type array using statistical method to get
!  vegtyp(nx,ny) (in iout)
!-----------------------------------------------------------------------

    CALL getvtyp1(xdims,ydims, nx,ny,glon,glat,                     &
                  x1,y1,idtarr,xdat,ydat,iout)

  ELSE  ! flag = 3

!-----------------------------------------------------------------------
!  Calculate real NDVI type data using statistical method to get
!  ndvi(nx,ny) (in rout)
!-----------------------------------------------------------------------

    CALL getndvi1(xdims,ydims,nx,ny,                              &
                  x1,y1,idtarr,xdat,ydat,rout)

  ENDIF

  DEALLOCATE(xdat)
  DEALLOCATE(ydat)
  DEALLOCATE(idtarr)

  RETURN
END SUBROUTINE get_statsgo
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE GETSTYP1                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE getstyp1(xdims,ydims,nx,ny,nstyps,glon,glat,                 &
           x1,y1,stypdat,xdat,ydat,soiltyp,stypfrct)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Transfer the soil classes data (soiltyp) into the USDA soil type
!  category and fill them into the model grid.
!
!-----------------------------------------------------------------------
!
!  AUTHOR:
!  Leilei Wang ,Vince Wong
!  3/20/97
!
!-----------------------------------------------------------------------
!  MODIFICATIONS:
!
!  Yunheng Wang
!  8/04/01
!  Rewrite some code in Fortran 90.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!  nx       Number of grid points in the x-direction
!  ny       Number of grid points in the y-direction
!  xdims    X-dimsion size of the data set
!  ydims    Y-dimsion size of the data set
!
!  glon     Longitude values of model grid points
!  glat     Latitude  values of model grid points
!
!  dtlon    Longitude values of data grid points
!  dtlat    Latitude  values of data grid points
!
!  stypdat  Soil data array
!  sfrctdat Soil type fraction data array
!
!  OUTPUT:
!
!  soiltyp  Soil type array
!  stypfrct Soil type fraction array
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'arpssfc.inc'
  INCLUDE 'globcst.inc'
  INCLUDE 'grid.inc'
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  INTEGER :: nx           ! Number of grid points in the x-direction
  INTEGER :: ny           ! Number of grid points in the y-direction
  INTEGER :: xdims        ! X-dimension size of the data set
  INTEGER :: ydims        ! Y-dimension size of the data set
  REAL :: x1  (nx,ny)  ! X_coordinate of model grid
  REAL :: y1  (nx,ny)  ! X_coordinate of model grid

  INTEGER :: nstyps       ! Number of soil types in each grid point
!
  REAL :: xdat(xdims)  ! X_coor. of soil data
  REAL :: ydat(ydims)  ! Y_coor. of soil data
  REAL :: glon(nx,ny)  ! Lon. of  model grid
  REAL :: glat(nx,ny)  ! Lat. of  model grid

  INTEGER :: stypdat (xdims,ydims)        ! Soil type data array
  INTEGER :: soiltyp (nx,ny,nstyps)       ! Soil type array
  REAL :: stypfrct(nx,ny,nstyps)       ! Soil type fraction array
  REAL :: dresl                        ! data resolution
  INTEGER :: item1  (nsoiltyp)            ! Temporary array
  INTEGER :: item2  (nsoiltyp)            ! Temporary array
  REAL :: tem3  (nsoiltyp)             ! Temporary array
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,im, iis,is,ii,jj
  INTEGER :: imax
  INTEGER :: count
  INTEGER :: xleft,xright,yleft,yright
  REAL :: tot,temp1
  INTEGER :: itemp2
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!-----------------------------------------------------------------------
!
!  Fill the soil types into the model grid using statistics results
!  by chosing soil types of highest weights and their corresponding
!  fractions for each model grid point
!
!-----------------------------------------------------------------------
!
  dresl =1000.

  DO j = 1,ydims
    DO i = 1,xdims

      SELECT CASE (stypdat(i,j))
        CASE (5:12)
          stypdat(i,j) = stypdat(i,j) -1
        CASE (13)
          stypdat(i,j) = 4
        CASE (14)
          stypdat(i,j) = 13
        CASE (15, 16)
          stypdat(i,j) = 11
      END SELECT

    END DO
  END DO

  DO i = 1,nx
    DO j = 1,ny

      item1(:) = 0
      item2(:) = 0
      tem3(:)  = 0.0

      count = 0
      tot   = 0.0

      IF(dx < dresl) THEN
        xleft = INT( (x1(i,j) - xdat(1)-dresl/2)/dresl )+1
        xright= INT( (x1(i,j) - xdat(1)+dresl/2)/dresl )
      ELSE
        xleft = INT( (x1(i,j) - xdat(1)-dx/2)/dresl )+1
        xright= INT( (x1(i,j) - xdat(1)+dx/2)/dresl )
      END IF

      IF(dy < 1000.0) THEN
        yleft = INT( (y1(i,j) - ydat(1)-dresl/2)/dresl )+1
        yright= INT( (y1(i,j) - ydat(1)+dresl/2)/dresl )
      ELSE
        yleft = INT( (y1(i,j) - ydat(1)-dy/2)/dresl )+1
        yright= INT( (y1(i,j) - ydat(1)+dy/2)/dresl )
      END IF

      DO ii = xleft, xright
        DO jj = yleft, yright
          IF(xleft < 1.OR.yleft < 1) GO TO 800
          IF(xright > xdims.OR.yright > ydims) GO TO 800
          IF(stypdat(ii,jj) /= 0) count=count + 1
          DO is=1,nsoiltyp
            IF(stypdat(ii,jj) == is) THEN
              item1(is) = item1(is) + 1
              item2(is) = is
            END IF
          END DO
        800     CONTINUE
        END DO
      END DO

      IF(count == 0) GOTO 100

      tem3(:) = FLOAT(item1(:))/count

      DO iis=1,nsoiltyp
        imax= iis
        DO im=iis, nsoiltyp
          IF (tem3(im) > tem3(imax)) imax = im
        END DO
        temp1=tem3(iis)
        tem3(iis)=tem3(imax)
        tem3(imax)=temp1
        itemp2=item2(iis)
        item2(iis)=item2(imax)
        item2(imax)=itemp2
      END DO

      tot = SUM(tem3(1:nstyp))

      100   CONTINUE

      DO is=1,nstyp

        IF(tot /= 0.0)THEN
          soiltyp(i,j,is) = item2(is)
          stypfrct(i,j,is)= tem3(is)/tot
        ELSE
          stypfrct(i,j,is) = 0.0
          stypfrct(i,j,1) = 1.0
        END IF

      END DO

    END DO
  END DO

  DO is = 2, nstyp
     WHERE(soiltyp(:,:,is) == 0)
       soiltyp(:,:,is) = soiltyp (:,:,is-1)
       stypfrct(:,:,is) = 0.0
     END WHERE

     WHERE (soiltyp(:,:,is) == soiltyp(:,:,is-1))
       stypfrct(:,:,is) = 0.0
     END WHERE
  END DO

!
!-----------------------------------------------------------------------
!
!  The following DO loop is just for the Vortex experiment. On the
!  south boundary Vortex domain the vegetation type changes sharply,
!  that will pretty much influence the boundary conditions. Therefore
!  the following DO loop will use the values of inner grids to those
!  boundary grids.
!
!-----------------------------------------------------------------------
!
!  DO 300 j = 1, 7
!  DO 300 i = 1, nx
!    soiltyp(i,j) = soiltyp(i,8)
!300  CONTINUE

  RETURN
END SUBROUTINE getstyp1
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE GETSTYP2                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
!

SUBROUTINE getstyp2(xdims,ydims, nx,ny,                                 &
           dtlon,dtlat,resl, glon,glat,                                 &
           stypdat, soiltyp)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Translate the soil classes data (stypdat) into the USDA soil type
!  category and fill them into the model grid array soiltyp.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yuhe Liu
!
!  1/29/94
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!  nx       Number of grid points in the x-direction
!  ny       Number of grid points in the y-direction
!  xdims    X-dimsion size of the data set
!  ydims    Y-dimsion size of the data set
!
!  glon     Longitude values of model grid points
!  glat     Latitude  values of model grid points
!
!  dtlon    Longitude values of data grid points
!  dtlat    Latitude  values of data grid points
!
!  stypdat  Soil data array
!
!  OUTPUT:
!
!  soiltyp  Soil type array
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx           ! Number of grid points in the x-direction
  INTEGER :: ny           ! Number of grid points in the y-direction
  INTEGER :: xdims        ! X-dimsion size of the data set
  INTEGER :: ydims        ! Y-dimsion size of the data set

  REAL :: glon(nx,ny)  ! Longitude values of model grid points
  REAL :: glat(nx,ny)  ! Latitude  values of model grid points

  REAL :: dtlon(xdims) ! Longitude values of data grid points
  REAL :: dtlat(ydims) ! Latitude  values of data grid points
  REAL :: resl         ! Data resolution

  INTEGER :: stypdat(xdims,ydims)  ! Soil data array

  INTEGER :: soiltyp(nx,ny)        ! Soil type array

!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j, ii, jj, ix, jy
  INTEGER :: ixmin,ixmax,jymin,jymax
  REAL    :: r,r1

  LOGICAL :: find
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'arpssfc.inc'
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

  DO j = 1, ydims
    DO i = 1, xdims

      SELECT CASE (stypdat(i,j))

        CASE ( 11, 17, 23 )

          stypdat(i,j) = 2

        CASE ( 14, 20, 26, 27 )

          stypdat(i,j) = 3

        CASE ( 12, 18, 24 )

          stypdat(i,j) = 5

        CASE ( 15, 21, 28 )

          stypdat(i,j) = 6

        CASE ( 13 )

          stypdat(i,j) = 8

        CASE ( 19, 25 )

          stypdat(i,j) = 9

        CASE ( 16, 22 )

          stypdat(i,j) = 10

        CASE ( 29, 30, 31 )

          stypdat(i,j) = 11

        CASE ( 34 )

          stypdat(i,j) = 12

        CASE ( 00 )

          stypdat(i,j) = 13

        CASE DEFAULT

          WRITE (6, '(a/a,i4,a,i4,a,i2/a)')                           &
            'The value for soil type is invalid.',                    &
            'stypdat(', i, ',', j, ') = ', stypdat(i,j),              &
            'Program stops in subroutine GETSTYP2.'

          STOP

        END SELECT

    END DO
  END DO

!
!-----------------------------------------------------------------------
!
!  Fill the soil type into the model grid by chosing the values
!  at the nearest data grid points.
!
!-----------------------------------------------------------------------
!

!J.Case, ENSCO Inc -- (8/10/2004) -- Bug fix from version 4.5.2. Incorrect
!       soil/veg/NDVI types get assigned on hemispheric/large domains without
!       this fix.  Originally made by J. Manobianco.
!
  ixmin = xdims-1
  ixmax = 1
  DO j = 1, ny
    DO i = 1, nx
      ixmin = MIN( ixmin, INT( ( glon(i,j)-dtlon(1) ) / resl ) + 1)
      ixmax = MAX( ixmax, INT( ( glon(i,j)-dtlon(1) ) / resl ) + 1)
    END DO
  END DO
  ixmin = MAX( 1,     ixmin - 5 )
  ixmax = MIN( xdims, ixmax + 5 )

  jymin = ydims-1
  jymax = 1
  DO j = 1 , ny
    DO i = 1 , nx
      jymin = MIN( jymin, INT( ( glat(i,j)-dtlat(1) ) / resl ) + 1 )
      jymax = MAX( jymax, INT( ( glat(i,j)-dtlat(1) ) / resl ) + 1 )
    END DO
  END DO
  jymin = MAX( 1,     jymin - 5 )
  jymax = MIN( ydims, jymax + 5 )


  DO j = 1, ny
    DO i = 1, nx

      r1 = 1.0E20
      find = .FALSE.
      ii = 0
      jj = 0

      IF (ixmin .GT. ixmax) THEN       ! Greenwich Meridian in domain
        DO jy = jymin, jymax
          DO ix = ixmin, xdims
             r = ( dtlon(ix) - glon(i,j) ) ** 2                         &
               + ( dtlat(jy) - glat(i,j) ) ** 2
             IF(.NOT. find .AND. r < r1) find = .TRUE.
             IF(.NOT. find) write(6,*) 'Warnning in GETSTYP2 ......'
             IF ( r < r1 ) THEN
               r1 = r
               ii = ix
               jj = jy
             END IF
          END DO
          DO ix = 1,ixmax
             r = ( dtlon(ix) - glon(i,j) ) ** 2                         &
               + ( dtlat(jy) - glat(i,j) ) ** 2
             IF(.NOT. find .AND. r < r1) find = .TRUE.
             IF(.NOT. find) write(6,*) 'Warnning in GETSTYP2 ......'
             IF ( r < r1 ) THEN
               r1 = r
               ii = ix
               jj = jy
             END IF
          END DO
        END DO

      ELSE                             ! Greenwich Meridian not in domain

        DO jy = jymin, jymax
          DO ix = ixmin, ixmax
            r = ( dtlon(ix) - glon(i,j) ) ** 2                         &
              + ( dtlat(jy) - glat(i,j) ) ** 2
            IF(.NOT. find .AND. r < r1) find = .TRUE.
            IF(.NOT. find) write(6,*) 'Warnning in GETSTYP2 ......'
            IF ( r < r1 ) THEN
              r1 = r
              ii = ix
              jj = jy
            END IF
          END DO
        END DO
      END IF

      IF(find) THEN
        soiltyp(i,j) = stypdat(ii,jj)
      ELSE
        WRITE(6,*) 'GETSTYP2: Error in getstyp2 for some machine',    &
                   ' Please try again: makearps -opt 1 arpssfc'
        STOP
      END IF

    END DO
  END DO

!
!-----------------------------------------------------------------------
!
!  The following DO loop is just for the Vortex experiment. On the
!  south boundary Vortex domain the vegetation type changes sharply,
!  that will pretty much influence the boundary conditions. Therefore
!  the following DO loop will use the values of inner grids to those
!  boundary grids.
!
!-----------------------------------------------------------------------
!
!  DO 300 j = 1, 7
!  DO 300 i = 1, nx
!    soiltyp(i,j) = soiltyp(i,8)
!300  CONTINUE

  RETURN
END SUBROUTINE getstyp2
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE GETVTYP1                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE getvtyp1(xdims,ydims,nx,ny,glon,glat,                        &
           x1,y1,vtypdat,xdat,ydat,vegtyp)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Converse the OWE classes of vegetation type (vegtyp)into ARPS
!  type category and fill them into the model grid.
!
!-----------------------------------------------------------------------
!
!  AUTHOR:
!  Leilei Wang ,Vince Wong
!  8/11/97
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!  nx       Number of grid points in the x-direction
!  ny       Number of grid points in the y-direction
!  xdims    X-dimsion size of the data set
!  ydims    Y-dimsion size of the data set
!
!  glon     Longitude values of model grid points
!  glat     Latitude  values of model grid points
!
!  dtlon    Longitude values of data grid points
!  dtlat    Latitude  values of data grid points
!
!  nvegtyp  Number of soil types in STATESGO class
!  vtypdat  Vegetation data array
!
!  OUTPUT:
!
!  vegtyp   Vegetation type array
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'arpssfc.inc'
  INCLUDE 'globcst.inc'
  INCLUDE 'grid.inc'
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  INTEGER :: nx           ! Number of grid points in the x-direction
  INTEGER :: ny           ! Number of grid points in the y-direction
  INTEGER :: xdims        ! X-dimension size of the data set
  INTEGER :: ydims        ! Y-dimension size of the data set
  REAL :: x1  (nx,ny)  ! X_coordinate of model grid
  REAL :: y1  (nx,ny)  ! X_coordinate of model grid

  REAL :: xdat(xdims)  ! X_coor. of soil data
  REAL :: ydat(ydims)  ! Y_coor. of soil data
  REAL :: glon(nx,ny)  ! Lon. of  model grid
  REAL :: glat(nx,ny)  ! Lat. of  model grid

  REAL :: dresl                        ! data resolution
  INTEGER :: vtypdat (xdims,ydims)        ! Veg. type data array
  INTEGER :: vegtyp  (nx,ny)              ! Veg. type array

  INTEGER :: item2(nvegtyp)               ! Temporary array
  INTEGER :: item3(nvegtyp)               ! Temporary array
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j, iis,is,ii,jj
  INTEGER :: imax,count
  INTEGER :: xleft,xright,yleft,yright
  INTEGER, ALLOCATABLE :: item1(:,:)     ! add by wyh
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!-----------------------------------------------------------------------
!
!  Fill the veg. types into the model grid using statistics results
!  by chosing veg. types of highest weights and their corresponding
!  fractions for each model grid point
!
!-----------------------------------------------------------------------
!
  ALLOCATE(item1(xdims, ydims))

  dresl =1000.

  DO j = 1,ydims
    DO i = 1,xdims

      SELECT CASE (vtypdat(i,j))

        CASE ( 01, 08, 50, 69, 71 )

          item1(i,j) = 1

        CASE ( 42, 09, 53 )

          item1(i,j) = 2

        CASE ( 16, 30, 37, 40, 10, 76, 93, 52 )

          item1(i,j) = 3

        CASE ( 02, 41, 43, 49, 07, 82, 83, 85, 94, 87 )

          item1(i,j) = 4

        CASE ( 18, 91, 58 )

          item1(i,j) = 5

        CASE ( 24:27, 29, 56, 61, 03:05, 19, 78, 88, 89 )

          item1(i,j) = 6

        CASE ( 06, 20:23, 46:48, 57, 60, 62, 34, 77 )

          item1(i,j) = 7

        CASE ( 32:33, 54, 79, 86, 90 )

          item1(i,j) = 8

        CASE ( 17, 70 )

          item1(i,j) = 9

        CASE ( 28, 31, 36, 38, 39, 55, 35, 92 )

          item1(i,j) = 10

        CASE ( 44, 45, 13, 74, 75, 80 )

          item1(i,j) = 11

        CASE ( 59, 63, 12, 64 )

          item1(i,j) = 12

        CASE ( 51, 81, 84, 11 )

          item1(i,j) = 13

        CASE ( 00, 65:68, 72, 73, 14, 15 )

          item1(i,j) = 14

        CASE DEFAULT

          WRITE (6, '(a/a,i4,a,i4,a,i2/a)')                               &
            'The value for vegetation type is invalid.',                &
            'vtypdat(', i, ',', j, ') = ', vtypdat(i,j),                &
            'Program stops in subroutine GETVTYP1.'

          STOP

      END SELECT

    END DO
  END DO

  DO j = 1,ny
    DO i = 1,nx

      item2=0
      item3=0

      count = 0
      imax = 0

      IF(dx < dresl) THEN
        xleft = INT( (x1(i,j) - xdat(1)-dresl/2)/dresl )+1
        xright= INT( (x1(i,j) - xdat(1)+dresl/2)/dresl )
      ELSE
        xleft = INT( (x1(i,j) - xdat(1)-dx/2)/dresl )+1
        xright= INT( (x1(i,j) - xdat(1)+dx/2)/dresl )
      END IF

      IF(dy < 1000.0) THEN
        yleft = INT( (y1(i,j) - ydat(1)-dresl/2)/dresl )+1
        yright= INT( (y1(i,j) - ydat(1)+dresl/2)/dresl )
      ELSE
        yleft = INT( (y1(i,j) - ydat(1)-dy/2)/dresl )+1
        yright= INT( (y1(i,j) - ydat(1)+dy/2)/dresl )
      END IF

      DO ii = xleft, xright
        DO jj = yleft, yright

          IF(xleft < 1.OR.yleft < 1) GO TO 800
          IF(xright > xdims.OR.yright > ydims) GO TO 800

          IF(item1(ii,jj) /= 0) count=count + 1
          DO is=1,nvegtyp
            IF(item1(ii,jj) == is) THEN
              item2  (is) = item2(is) + 1
              item3  (is) = is
            END IF
          END DO
800     CONTINUE
        END DO
      END DO

      IF(count /= 0)  THEN
        imax = 1
        DO iis = 1, nvegtyp
          IF(item2(iis) >= item2(imax))  imax=iis
        END DO
      END IF

      IF (imax /= 0) vegtyp(i,j) = item3(imax)

    END DO
  END DO
!
!-----------------------------------------------------------------------
!
!  The following DO loop is just for the Vortex experiment. On the
!  south boundary Vortex domain the vegetation type changes sharply,
!  that will pretty much influence the boundary conditions. Therefore
!  the following DO loop will use the values of inner grids to those
!  boundary grids.
!
!-----------------------------------------------------------------------
!
!  DO 300 j = 1, 7
!  DO 300 i = 1, nx
!    vegtyp(i,j) = vegtyp(i,8)
!300  CONTINUE

  DEALLOCATE(item1)

  RETURN
END SUBROUTINE getvtyp1
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE GETVTYP2                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
!

SUBROUTINE getvtyp2(xdims,ydims, nx,ny,                                 &
           dtlon,dtlat,resl, glon,glat,                                 &
           vtypdat, vegtyp)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Transit the Olson World Ecosystem Classes of vegetation data into
!  the simple vegetation type categories and fill them into the model
!  domain by chosing the values at the nearest data grid points.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yuhe Liu
!
!  2/20/94
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!  nx       Number of grid points in the x-direction
!  ny       Number of grid points in the y-direction
!  xdims    X-dimsion size of the data set
!  ydims    Y-dimsion size of the data set
!
!  glon     Longitude values of model grid points
!  glat     Latitude  values of model grid points
!
!  dtlon    Longitude values of data grid points
!  dtlat    Latitude  values of data grid points
!
!  vtyp     Vegetation data array
!
!  OUTPUT:
!
!  vegtyp   Vegetation type array
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
  INTEGER :: nx           ! Number of grid points in the x-direction
  INTEGER :: ny           ! Number of grid points in the y-direction
  INTEGER :: xdims        ! X-dimsion size of the data set
  INTEGER :: ydims        ! Y-dimsion size of the data set
!
  REAL :: glon(nx,ny)  ! Longitude values of model grid points
  REAL :: glat(nx,ny)  ! Latitude  values of model grid points
!
  REAL :: dtlon(xdims) ! Longitude values of data grid points
  REAL :: dtlat(ydims) ! Latitude  values of data grid points
  REAL :: resl         ! Data resolution
!
  INTEGER :: vtypdat(xdims,ydims)  ! Soil data array
!
  INTEGER :: vegtyp(nx,ny)      ! Soil type array

!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j, ii, jj, ix, jy
  INTEGER :: ixmin,ixmax,jymin,jymax
  REAL :: r,r1

!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'arpssfc.inc'
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

  DO j = 1, ydims
    DO i = 1, xdims

      IF ( vtypdat(i,j) == 01 .OR.  vtypdat(i,j) == 08 .OR.             &
             vtypdat(i,j) == 50 .OR.                                    &
             vtypdat(i,j) == 69 .OR.                                    &
             vtypdat(i,j) == 71 ) THEN

        vtypdat(i,j) = 1

      ELSE IF ( vtypdat(i,j) == 42 .OR.                                 &
                 vtypdat(i,j) == 53 ) THEN

        vtypdat(i,j) = 2

      ELSE IF ( vtypdat(i,j) == 16 .OR.                                 &
                 vtypdat(i,j) == 30 .OR.                                &
                 vtypdat(i,j) == 37 .OR.                                &
                 vtypdat(i,j) == 40 .OR.                                &
                 vtypdat(i,j) == 52 ) THEN

        vtypdat(i,j) = 3

      ELSE IF ( vtypdat(i,j) == 02 .OR.                                 &
                 vtypdat(i,j) == 41 .OR.                                &
                 vtypdat(i,j) == 43 .OR.                                &
                 vtypdat(i,j) == 49 ) THEN

        vtypdat(i,j) = 4

      ELSE IF ( vtypdat(i,j) == 58 .OR.                                 &
                 vtypdat(i,j) == 47 ) THEN

        vtypdat(i,j) = 5

      ELSE IF ( vtypdat(i,j) == 24 .OR.                                 &
                 vtypdat(i,j) == 25 .OR.                                &
                 vtypdat(i,j) == 26 .OR.                                &
                 vtypdat(i,j) == 27 .OR.                                &
                 vtypdat(i,j) == 29 .OR.                                &
                 vtypdat(i,j) == 56 .OR.                                &
                 vtypdat(i,j) == 61 ) THEN

        vtypdat(i,j) = 6

      ELSE IF ( vtypdat(i,j) == 06 .OR.                                 &
                 vtypdat(i,j) == 20 .OR.                                &
                 vtypdat(i,j) == 21 .OR.                                &
                 vtypdat(i,j) == 22 .OR.                                &
                 vtypdat(i,j) == 23 .OR.                                &
                 vtypdat(i,j) == 46 .OR.                                &
                 vtypdat(i,j) == 48 .OR.                                &
                 vtypdat(i,j) == 57 .OR.                                &
                 vtypdat(i,j) == 60 .OR.                                &
                 vtypdat(i,j) == 62 ) THEN

        vtypdat(i,j) = 7

      ELSE IF ( vtypdat(i,j) == 32 .OR.                                 &
                 vtypdat(i,j) == 33 .OR.                                &
                 vtypdat(i,j) == 54 ) THEN

        vtypdat(i,j) = 8

      ELSE IF ( vtypdat(i,j) == 17 .OR.                                 &
                 vtypdat(i,j) == 70 ) THEN

        vtypdat(i,j) = 9

      ELSE IF ( vtypdat(i,j) == 28 .OR.                                 &
                 vtypdat(i,j) == 31 .OR.                                &
                 vtypdat(i,j) == 36 .OR.                                &
                 vtypdat(i,j) == 38 .OR.                                &
                 vtypdat(i,j) == 39 .OR.                                &
                 vtypdat(i,j) == 55 ) THEN

        vtypdat(i,j) = 10

      ELSE IF ( vtypdat(i,j) == 44 .OR.                                 &
                 vtypdat(i,j) == 45 ) THEN

        vtypdat(i,j) = 11

      ELSE IF ( vtypdat(i,j) == 59 .OR.                                 &
                 vtypdat(i,j) == 63 .OR.                                &
                 vtypdat(i,j) == 64 ) THEN

        vtypdat(i,j) = 12

      ELSE IF ( vtypdat(i,j) == 51 ) THEN

        vtypdat(i,j) = 13

      ELSE IF ( vtypdat(i,j) == 00 .OR.                                 &
                 vtypdat(i,j) == 65 .OR.                                &
                 vtypdat(i,j) == 66 .OR.                                &
                 vtypdat(i,j) == 67 .OR.                                &
                 vtypdat(i,j) == 68 .OR.                                &
                 vtypdat(i,j) == 72 .OR.                                &
                 vtypdat(i,j) == 73 ) THEN

        vtypdat(i,j) = 14
      ELSE

        WRITE (6, '(a/a,i4,a,i4,a,i2/a)')                               &
            'The value for vegetation type is invalid.',                &
            'vtypdat(', i, ',', j, ') = ', vtypdat(i,j),                &
            'Program stops in subroutine GETVTYP.'

        STOP

      END IF

    END DO
  END DO


!
!-----------------------------------------------------------------------
!
!  Fill the vegetation type into the model grid by chosing the values
!  at the nearest data grid points.
!
!-----------------------------------------------------------------------
!

!J.Case, ENSCO Inc. -- (8/10/2004) -- Bug fix from version 4.5.2. Incorrect soil/veg/NDVI types
!  get assigned on hemispheric/large domains without this fix.  Originally made by J. Manobianco.
!
  ixmin = xdims-1
  ixmax = 1
  DO j = 1, ny
    DO i = 1, nx
      ixmin = MIN( ixmin, INT( ( glon(i,j)-dtlon(1) ) / resl ) + 1 )
      ixmax = MAX( ixmax, INT( ( glon(i,j)-dtlon(1) ) / resl ) + 1 )
    END DO
  END DO
  ixmin = MAX( 1,     ixmin - 5 )
  ixmax = MIN( xdims, ixmax + 5 )

  jymin = ydims-1
  jymax = 1
  DO j = 1, ny
    DO i = 1, nx
      jymin = MIN( jymin, INT( ( glat(i,j)-dtlat(1) ) / resl ) + 1 )
      jymax = MAX( jymax, INT( ( glat(i,j)-dtlat(1) ) / resl ) + 1 )
    END DO
  END DO
  jymin = MAX( 1,     jymin - 5 )
  jymax = MIN( ydims, jymax + 5 )


  DO j = 1, ny
    DO i = 1, nx

      r1 = 1.0E20

      IF (ixmin .GT. ixmax) THEN       ! Greenwich Meridian in domain
        DO jy = jymin, jymax
          DO ix = ixmin, xdims
             r = ( dtlon(ix) - glon(i,j) ) ** 2                         &
               + ( dtlat(jy) - glat(i,j) ) ** 2
             IF ( r < r1 ) THEN
               r1 = r
               ii = ix
               jj = jy
             END IF
          END DO
          DO ix = 1,ixmax
             r = ( dtlon(ix) - glon(i,j) ) ** 2                         &
               + ( dtlat(jy) - glat(i,j) ) ** 2
             IF ( r < r1 ) THEN
               r1 = r
               ii = ix
               jj = jy
             END IF
          END DO
        END DO

      ELSE                             ! Greenwich Meridian not in domain
        DO jy = jymin, jymax
          DO ix = ixmin, ixmax
            r = ( dtlon(ix) - glon(i,j) ) ** 2                            &
              + ( dtlat(jy) - glat(i,j) ) ** 2
            IF ( r < r1 ) THEN
              r1 = r
              ii = ix
              jj = jy
            END IF
          END DO
        END DO
      END IF

      vegtyp(i,j) = vtypdat(ii,jj)

    END DO
  END DO

!
!-----------------------------------------------------------------------
!
!  The following DO loop is just for the Vortex experiment. On the
!  south boundary Vortex domain the vegetation type changes sharply,
!  that will pretty much influence the boundary conditions. Therefore
!  the following DO loop will use the values of inner grids to those
!  boundary grids.
!
!-----------------------------------------------------------------------
!
!  DO 300 j = 1, 7
!  DO 300 i = 1, nx
!    vegtyp(i,j) = vegtyp(i,8)
!300  CONTINUE

  RETURN
END SUBROUTINE getvtyp2
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE GETLAI                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
!

SUBROUTINE getlai( nx,ny, vtypdat, ndvi, lai )
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate the Leaf Area Index (lai) from NDVI data and interpolate
!  them into the model domain.
!
!  The calculation of LAI depends on vegetation type: grass or tree.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yuhe Liu
!
!  3/22/94
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!  nx       Number of grid points in the x-direction
!  ny       Number of grid points in the y-direction
!
!  vtyp     Vegetation type;
!  ndvi     NDVI data array
!
!  OUTPUT:
!
!  lai      LAI array
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
  INTEGER :: nx           ! Number of grid points in the x-direction
  INTEGER :: ny           ! Number of grid points in the y-direction
!
  INTEGER :: vtypdat(nx,ny)  ! Vegetation type: 1 -- grass; 2 -- tree
  REAL :: ndvi(nx,ny)  ! Soil data array
!
  REAL :: lai (nx,ny)  ! Soil type array
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'arpssfc.inc'
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  DO j = 1, ny
    DO i = 1, nx
      IF ( vtypdat(i,j) == 09 .OR.             & ! Water & ice
           vtypdat(i,j) == 14  ) THEN
        lai(i,j) = 0.

      ELSE IF ( vtypdat(i,j) == 01 .OR.                                 &
                vtypdat(i,j) == 02 .OR.                                 &
                vtypdat(i,j) == 03 .OR.                                 &
                vtypdat(i,j) == 04 .OR.                                 &
                vtypdat(i,j) == 05 .OR.                                 &
                vtypdat(i,j) == 09 .OR.                                 &
                vtypdat(i,j) == 10 .OR.                                 &
                vtypdat(i,j) == 11 .OR.                                 &
                vtypdat(i,j) == 12 .OR.                                 &
                vtypdat(i,j) == 13 .OR.                                 &
                vtypdat(i,j) == 14 ) THEN

        lai(i,j) = - LOG( ( 1. - ndvi(i,j)/.915 ) / .83 ) / 0.96

!    ELSEIF ( vtypdat(i,j) .eq. 04 .or.  ! uncomment this block for 4 & 5
!    :           vtypdat(i,j) .eq. 05 ) THEN
!
!      lai(i,j) = 0.5
!    :             * ( - log( ( 1. - ndvi(i,j)/.915 ) / .83 ) / 0.96
!    :               + 1.625 * exp( ndvi(i,j) / 0.34 ) )

    ELSE

      lai(i,j) = 1.625 * EXP( ndvi(i,j) / 0.34 )

    END IF

    lai(i,j) = MAX( lai(i,j), 0. )

  END DO
  END DO

  RETURN
END SUBROUTINE getlai
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE GETNDVI1                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE getndvi1(xdims,ydims,nx,ny,                                  &
           x1,y1,ndvibyt,xdat,ydat,ndvi)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate the NDVI data (ndvi) from byte-NDVI data and
!  fill them into the model domain by chosing the average
!  value of most frequently appeared in a grid.
!
!-----------------------------------------------------------------------
!
!  AUTHOR:
!  Leilei Wang ,Vince Wong
!  8/11/97
!
!  01/20/2008 (Keith Brewster)
!  Corrected out-of-bounds array condition on tem2 and
!  cleaned up logic nearby to eliminate GO-TO's.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!  nx       Number of grid points in the x-direction
!  ny       Number of grid points in the y-direction
!  xdims    X-dimsion size of the data set
!  ydims    Y-dimsion size of the data set
!
!  glon     Longitude values of model grid points
!  glat     Latitude  values of model grid points
!
!  dtlon    Longitude values of data grid points
!  dtlat    Latitude  values of data grid points
!
!  ndvibyt  Ndvi data array
!
!  OUTPUT:
!
!  ndvi     Ndvi type array
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'arpssfc.inc'
  INCLUDE 'globcst.inc'
  INCLUDE 'grid.inc'
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  INTEGER :: nx           ! Number of grid points in the x-direction
  INTEGER :: ny           ! Number of grid points in the y-direction
  INTEGER :: xdims        ! X-dimension size of the data set
  INTEGER :: ydims        ! Y-dimension size of the data set
  REAL :: x1  (nx,ny)  ! X_coordinate of model grid
  REAL :: y1  (nx,ny)  ! X_coordinate of model grid

  REAL :: xdat(xdims)  ! X_coor. of soil data
  REAL :: ydat(ydims)  ! Y_coor. of soil data

  INTEGER :: ndvibyt (xdims,ydims)        ! byte_NDVI data array
  REAL :: ndvi    (nx,ny)              ! NDVI  array
  REAL :: dresl                        ! data resolution

  REAL :: tem2    (20)                 ! Temporary array
  REAL :: tem3    (20)                 ! Temporary array
  INTEGER :: item4(20)
  REAL, ALLOCATABLE :: tem1(:,:)                         ! add by wyh
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j, iis,is,ii,jj
  INTEGER :: imax,count
  INTEGER :: xleft,xright,yleft,yright
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

  ALLOCATE(tem1(xdims, ydims))                        ! add by wyh

  DO j = 1, ydims
    DO i = 1, xdims

      tem1(i,j) = MAX( 0.0, ( ndvibyt(i,j) - 100. ) / 100. )

    END DO
  END DO

!-----------------------------------------------------------------------
!
!  Fill the NDVI into the model grid using statistics results
!  by chosing NDVI  of highest weights  for each model grid point
!
!-----------------------------------------------------------------------
!
  dresl =1000.

  DO j = 1,ny
    DO i = 1,nx

      DO is = 1,20
        tem2(is) = -1.0+(is-1)*0.1
        tem3(is) = 0.0
        item4(is) = 0
      END DO

      count = 0
      imax = 0

      IF(dx < dresl) THEN
        xleft = INT( (x1(i,j) - xdat(1)-dresl/2)/dresl )+1
        xright= INT( (x1(i,j) - xdat(1)+dresl/2)/dresl )
      ELSE
        xleft = INT( (x1(i,j) - xdat(1)-dx/2)/dresl )+1
        xright= INT( (x1(i,j) - xdat(1)+dx/2)/dresl )
      END IF

      IF(dy < 1000.0) THEN
        yleft = INT( (y1(i,j) - ydat(1)-dresl/2)/dresl )+1
        yright= INT( (y1(i,j) - ydat(1)+dresl/2)/dresl )
      ELSE
        yleft = INT( (y1(i,j) - ydat(1)-dy/2)/dresl )+1
        yright= INT( (y1(i,j) - ydat(1)+dy/2)/dresl )
      END IF

      DO ii = xleft, xright
        DO jj = yleft, yright

          IF(xleft < 1.OR.yleft < 1) CYCLE
          IF(xright > xdims.OR.yright > ydims) CYCLE

          IF(tem1(ii,jj) >= -1.0.AND.tem1(ii,jj) <= 1.0) count=count + 1
          DO is=1,19
            IF(tem1(ii,jj) >= tem2(is) .AND. tem1(ii,jj) < tem2(is+1)) THEN
              tem3(is) = tem3(is) + tem1(ii,jj)
              item4(is) = item4(is) + 1
            END IF
          END DO
          IF(tem1(ii,jj) >= tem2(20)) THEN
            tem3(20) = tem3(20) + tem1(ii,jj)
            item4(20) = item4(20) + 1
          END IF
        END DO
      END DO

      IF(count /= 0) THEN
        imax = 1
        DO iis=1,20
          IF(item4(iis) >= item4(imax)) THEN
            imax=iis
          END IF
        END DO
      END IF

      IF(imax /= 0) ndvi(i,j) = tem3(imax)/item4(imax)

    END DO
  END DO

  DEALLOCATE(tem1)                    ! add by wyh

  RETURN
END SUBROUTINE getndvi1
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE GETNDVI2                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
!

SUBROUTINE getndvi2(xdims,ydims, nx,ny,                                 &
           dtlon,dtlat,resl, glon,glat,                                 &
           ndvibyt, ndvi)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate the NDVI data (ndvi) from byte-NDVI data and
!  fill them into the model domain by chosing the values at the
!  nearest data grid points.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yuhe Liu
!
!  2/20/94
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!  nx       Number of grid points in the x-direction
!  ny       Number of grid points in the y-direction
!  xdims    X-dimsion size of the data set
!  ydims    Y-dimsion size of the data set
!
!  glon     Longitude values of model grid points
!  glat     Latitude  values of model grid points
!
!  dtlon    Longitude values of data grid points
!  dtlat    Latitude  values of data grid points
!
!  ndvibyt  Byte-NDVI data
!
!  OUTPUT:
!
!  ndvi     Vegetation type array
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
  INTEGER :: nx           ! Number of grid points in the x-direction
  INTEGER :: ny           ! Number of grid points in the y-direction
  INTEGER :: xdims        ! X-dimsion size of the data set
  INTEGER :: ydims        ! Y-dimsion size of the data set
!
  REAL :: glon(nx,ny)  ! Longitude values of model grid points
  REAL :: glat(nx,ny)  ! Latitude  values of model grid points
!
  REAL :: dtlon(xdims) ! Longitude values of data grid points
  REAL :: dtlat(ydims) ! Latitude  values of data grid points
  REAL :: resl         ! Data resolution
!
  INTEGER :: ndvibyt(xdims,ydims)  ! Byte-NDVI data array
!
  REAL :: ndvi(nx,ny)  ! NDVI data array

!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j, ii, jj, ix,jy
  INTEGER :: ixmin,ixmax,jymin,jymax
  REAL :: r,r1

  REAL, ALLOCATABLE :: tem1(:,:)                ! add by wyh
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'arpssfc.inc'
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

  ALLOCATE(tem1(xdims, ydims))                     ! add by wyh

  WHERE (ndvibyt > 100)
    tem1 = (ndvibyt - 100. ) / 100.
  ELSEWHERE
    tem1 = 0.0
  END WHERE

!
!-----------------------------------------------------------------------
!
!  Fill the vegetation type into the model grid by chosing the values
!  at the nearest data grid points.
!
!-----------------------------------------------------------------------
!

!J.Case, ENSCO Inc. -- (8/10/2004) -- Bug fix from version 4.5.2. Incorrect soil/veg/NDVI types
!  get assigned on hemispheric/large domains without this fix.  Originally made by J. Manobianco.
!
  ixmin = xdims-1
  ixmax = 1
  DO j = 1, ny
    DO i = 1, nx
      ixmin = MIN( ixmin, INT( ( glon(i,j)-dtlon(1) ) / resl ) + 1)
      ixmax = MAX( ixmax, INT( ( glon(i,j)-dtlon(1) ) / resl ) + 1)
    END DO
  END DO
  ixmin = MAX( 1,     ixmin - 5 )
  ixmax = MIN( xdims, ixmax + 5 )

  jymin = ydims-1
  jymax = 1
  DO j = 1, ny
    DO i = 1, nx
      jymin = MIN( jymin, INT( ( glat(i,j)-dtlat(1) ) / resl ) + 1 )
      jymax = MAX( jymax, INT( ( glat(i,j)-dtlat(1) ) / resl ) + 1)
    END DO
  END DO
  jymin = MAX( 1,     jymin - 5 )
  jymax = MIN( ydims, jymax + 5 )

  DO j = 1, ny
    DO i = 1, nx

      r1 = 1.0E20

      IF (ixmin > ixmax) THEN       ! Greenwich Meridian in domain
        DO jy = jymin, jymax
          DO ix = ixmin, xdims
             r = ( dtlon(ix) - glon(i,j) ) ** 2                         &
               + ( dtlat(jy) - glat(i,j) ) ** 2
             IF ( r < r1 ) THEN
               r1 = r
               ii = ix
               jj = jy
             END IF
          END DO
          DO ix = 1,ixmax
             r = ( dtlon(ix) - glon(i,j) ) ** 2                         &
               + ( dtlat(jy) - glat(i,j) ) ** 2
             IF ( r < r1 ) THEN
               r1 = r
               ii = ix
               jj = jy
             END IF
          END DO
        END DO

      ELSE                             ! Greenwich Meridian not in domain

        DO jy = jymin, jymax
          DO ix = ixmin, ixmax
            r = ( dtlon(ix) - glon(i,j) ) ** 2                            &
              + ( dtlat(jy) - glat(i,j) ) ** 2
            IF ( r < r1 ) THEN
              r1 = r
              ii = ix
              jj = jy
            END IF
          END DO
        END DO

      END IF

      ndvi(i,j) = tem1(ii,jj)

    END DO
  END DO

  DEALLOCATE(tem1)
  RETURN
END SUBROUTINE getndvi2
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE GETRFNS                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
!

SUBROUTINE getrfns( nx,ny, vtypdat, roufns )
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Generate the surface roughness (roufns) for ARPS model from the
!  table vs vegetation type.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yuhe Liu
!
!  4/19/94
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!  nx       Number of grid points in the x-direction
!  ny       Number of grid points in the y-direction
!
!  vtyp     Vegetation type;
!
!  OUTPUT:
!
!  roufns   LAI array
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
  INTEGER :: nx           ! Number of grid points in the x-direction
  INTEGER :: ny           ! Number of grid points in the y-direction
!
  INTEGER :: vtypdat(nx,ny)  ! Vegetation type: 1 -- grass; 2 -- tree
!
  REAL :: roufns(nx,ny)  ! Surface roughness

  REAL :: rfns(14)       ! Table data for surface roughness vs vtyp

  DATA rfns/ 0.011, 0.076, 0.075, 0.238, 0.563, 0.826, 1.089,           &
             2.653, 0.011, 0.075, 0.100, 0.856, 0.065, 0.002/

!  DATA rfns/ 0.002, 0.030, 0.030, 0.100, 0.200, 0.500, 0.750,           &
!             1.000, 0.005, 0.150, 0.100, 0.100, 0.050, 0.002/

!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'arpssfc.inc'
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  DO j = 1, ny
    DO i = 1, nx

      roufns(i,j) = rfns(vtypdat(i,j))

    END DO
  END DO

  RETURN
END SUBROUTINE getrfns
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE GETVEG                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
!

SUBROUTINE getveg( nx,ny, vtypdat, veg )
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Generate the vegetation fraction (roufns) for ARPS model from the
!  table vs vegetation type.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yuhe Liu
!
!  02/07/95
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!  nx       Number of grid points in the x-direction
!  ny       Number of grid points in the y-direction
!
!  vtypdat  Vegetation type;
!
!  OUTPUT:
!
!  veg      veg array
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
  INTEGER :: nx           ! Number of grid points in the x-direction
  INTEGER :: ny           ! Number of grid points in the y-direction
!
  INTEGER :: vtypdat(nx,ny)  ! Vegetation type: 1 -- grass; 2 -- tree
!
  REAL :: veg(nx,ny)      ! Vegetation fraction

  REAL :: vegtab(14)      ! Table data for surface roughness vs vtyp
  DATA vegtab/ 0.10, 0.10, 0.60, 0.40, 0.40, 0.90, 0.99,                &
               0.99, 0.01, 0.30, 0.99, 0.40, 0.20, 0.00/
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'arpssfc.inc'
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  DO j = 1, ny
    DO i = 1, nx

      veg(i,j) = vegtab(vtypdat(i,j))

    END DO
  END DO

  RETURN
END SUBROUTINE getveg
!
!##################################################################
!##################################################################
!######                                                      ######
!######               SUBROUTINE GTLKUPTBL                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
!

SUBROUTINE gtlkuptbl( lkupfl, vegtbl, laitbl, z0tbl )
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Read in vegtbl, a table of veg, roufns, & lai as a funtion of
!  vegitation type.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Gene Bassett
!
!  02/20/97
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!  lkupfl   File name of the vegitation table.
!
!  OUTPUT:
!
!  vegtbl   Mapping of vegitation fraction as a fuction of
!           vegitation type.
!
!  laitbl   Mapping of leaf area index as a fuction of
!           vegitation type.
!
!  vegtbl   Mapping of vegetation fraction as a fuction of
!           vegetation type.
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
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'arpssfc.inc'
!
  REAL :: vegtbl(nvegtyp) ! veg as a fuction of vtype.
  REAL :: laitbl(nvegtyp) ! lai as a fuction of vtype.
  REAL :: z0tbl(nvegtyp)  ! z0 (roufns) as a fuction of vtype.

  CHARACTER (LEN=*) :: lkupfl ! File name of vegetation table.

  INTEGER :: iv, vin
  REAL    :: vegin, laiin, z0in
  CHARACTER (LEN=80) :: instr
  INTEGER :: n, nskip
  INTEGER :: lenfl, flunit
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

!
!-----------------------------------------------------------------------
!
!  Initialize variable.
!
!-----------------------------------------------------------------------
!
  DO iv = 1,nvegtyp

    vegtbl(iv) = 0.0
    laitbl(iv) = 0.0
    z0tbl(iv) = 0.000001

  END DO

  nskip = 0
!
!-----------------------------------------------------------------------
!
!  Read in the vegetation table.
!
!-----------------------------------------------------------------------
!
!  lenfl = 70
!  CALL strlnth( lkupfl, lenfl )
  lenfl = LEN_TRIM(lkupfl)

  CALL getunit( flunit )
  WRITE(6,'(a)') '    Opening vegetation lookup table file '            &
                 //lkupfl(1:lenfl)
  OPEN (UNIT = flunit, FILE = lkupfl(1:lenfl), STATUS = 'old')

  100   CONTINUE

  READ (flunit,'(a80)') instr

  IF ((instr(1:1) /= "c") .AND. (instr(1:1) /= "c") .AND.               &
      (instr(1:1) /= "&") .AND. (instr(1:1) /= "!") .AND.               &
      (instr(1:1) /= "#")) GO TO 700

  nskip = nskip + 1
  GO TO 100

  700   CONTINUE

  CLOSE ( flunit )
  OPEN (UNIT = flunit, FILE = lkupfl(1:lenfl), STATUS = 'old')

  DO n = 1,nskip

    READ (flunit,'(a80)') instr

  END DO

  500   READ (flunit,*,END=200) vin, vegin, laiin, z0in

  IF ( vin > nvegtyp .OR. vin < 1 ) THEN
    WRITE (*,*) "ERROR, entry skipped.  Veg type out of bounds:",vin
  END IF

  vegtbl(vin) = vegin
  laitbl(vin) = laiin
  z0tbl(vin) = z0in

  GO TO 500

  200   CONTINUE

  CLOSE ( flunit )
  CALL retunit ( flunit )

  WRITE (*,*) "Vegataion table used:"
  WRITE (*,*) "vtype   veg     lai     z0"

  DO n = 1,nvegtyp

    WRITE (*,1000) n,vegtbl(n),laitbl(n),z0tbl(n)

  END DO

  1000  FORMAT(i2,7X,3F8.5)

  RETURN
END SUBROUTINE gtlkuptbl
!
!##################################################################
!##################################################################
!######                                                      ######
!######               SUBROUTINE GTLAITBL                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
!

SUBROUTINE gtlaitbl( nx,ny, vtyp, laitbl, lai )
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Assign the Leaf Area Index (lai) as a function of vegetation type.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Gene Bassett & Richard Carpenter
!
!  02/20/97
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!  nx       Number of grid points in the x-direction
!  ny       Number of grid points in the y-direction
!
!  vtyp     Vegetation type;
!  laitbl
!
!  OUTPUT:
!
!  lai      LAI array
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
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'arpssfc.inc'

  INTEGER :: nx           ! Number of grid points in the x-direction
  INTEGER :: ny           ! Number of grid points in the y-direction

  INTEGER :: vtyp(nx,ny)  ! Vegetation type
  REAL :: laitbl(nvegtyp)

  REAL :: lai (nx,ny)  ! Soil type array
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

  DO j = 1, ny
    DO i = 1, nx

      lai(i,j) = laitbl(vtyp(i,j))

    END DO
  END DO

  RETURN
END SUBROUTINE gtlaitbl
!
!##################################################################
!##################################################################
!######                                                      ######
!######               SUBROUTINE GTVEGTBL                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
!

SUBROUTINE gtvegtbl( nx,ny, vtyp, vegtbl, veg )
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Assign the vegetation fraction (veg) as a function of vegetation type.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Gene Bassett & Richard Carpenter
!
!  02/20/97
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!  nx       Number of grid points in the x-direction
!  ny       Number of grid points in the y-direction
!
!  vtyp     Vegetation type;
!  vegtbl
!
!  OUTPUT:
!
!  veg      vegetation fraction
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
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'arpssfc.inc'

  INTEGER :: nx           ! Number of grid points in the x-direction
  INTEGER :: ny           ! Number of grid points in the y-direction

  INTEGER :: vtyp(nx,ny)  ! Vegetation type
  REAL :: vegtbl(nvegtyp)

  REAL :: veg (nx,ny)  ! Soil type array
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

  DO j = 1, ny
    DO i = 1, nx

      veg(i,j) = vegtbl(vtyp(i,j))

    END DO
  END DO

  RETURN
END SUBROUTINE gtvegtbl
!
!##################################################################
!##################################################################
!######                                                      ######
!######               SUBROUTINE GTRFNSTBL                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
!

SUBROUTINE gtrfnstbl( nx,ny, vtyp, z0tbl, roufns )
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Assign the roughness length (roufns) as a function of vegetation type.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Gene Bassett & Richard Carpenter
!
!  02/20/97
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!  nx       Number of grid points in the x-direction
!  ny       Number of grid points in the y-direction
!
!  vtyp     Vegetation type;
!  z0tbl
!
!  OUTPUT:
!
!  roufns      roughness length
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
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'arpssfc.inc'

  INTEGER :: nx           ! Number of grid points in the x-direction
  INTEGER :: ny           ! Number of grid points in the y-direction

  INTEGER :: vtyp(nx,ny)  ! Vegetation type
  REAL :: z0tbl(nvegtyp)

  REAL :: roufns (nx,ny)  ! Soil type array
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

  DO j = 1, ny
    DO i = 1, nx

      roufns(i,j) = z0tbl(vtyp(i,j))

    END DO
  END DO

  RETURN
END SUBROUTINE gtrfnstbl
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE READBSQ                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE readbsq (filename, lun, inband, inrow, incol, inbyte,        &
           iband, irow, icol, nband, nrow, ncol, outarray,colomn,roww)
!
!-----------------------------------------------------------------------
!
!  Purpose:
!  This subroutine extracts a specified window (rectangular subset)
!  from a band-sequential (BSQ) binary data array and places it into
!  an INTEGER array specified by the calling program.
!
!-----------------------------------------------------------------------
!
!  01/96   Initial version, R. A. White, PSU/ESSC,
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!  INTEGER     inrow     Number of rows in input data array.
!  INTEGER     incol     Number of columns in input data array.
!  INTEGER     inbyte    Number of bytes per input data value.
!  INTEGER     iband     Initial band for which data are desired.
!  INTEGER     irow      Initial row for which data are desired:
!  INTEGER     icol      Initial column for which data are
!  INTEGER     nband     Number of consecutive bands to be
!                          extracted from input data array
!  INTEGER     nrow      Number of consecutive rows to be
!                          extracted from each band of data array.
!  INTEGER     ncol      Number of consecutive columns to be
!                          extracted from each row of data array.
!  OUTPUT:
!
!  INTEGER     outarray(ncol,nrow,nband)
!                            Array to receive the data from input
!                            file.  If nband = 1, the array may be
!                            declared as outarray(ncol,nrow).
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
!  PARAMETER (maxinbyte =16384)              !  comment out by wyh

  CHARACTER (LEN=*) :: filename
  INTEGER :: lun, inband, inrow, incol, inbyte, iband, irow, icol
  INTEGER :: nband, nrow, ncol
  INTEGER :: colomn,roww
  INTEGER :: outarray(colomn,roww)

  INTEGER :: band, row, col
  CHARACTER (LEN=1) :: inbuf(incol*inbyte)    !  add by wyh

  INTEGER :: lenfl
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  lenfl = LEN( filename )
  CALL strlnth( filename,lenfl )

  WRITE(6,'(a,a)') '    Opening ',filename(1:lenfl)
  OPEN (UNIT=lun, FILE=filename(1:lenfl), ACCESS='DIRECT',              &
                  FORM='UNFORMATTED',     STATUS='OLD',                 &
                  RECL=incol*inbyte)

  DO band=iband,iband+nband-1
    DO row=irow,irow+nrow-1
      READ (lun, REC=row+(band-1)*inrow) (inbuf(i), i=1,incol*inbyte)
      DO col=icol,icol+ncol-1
        outarray(col-icol+1, row-irow+1) = ICHAR(inbuf(col))
      END DO
    END DO
  END DO
  CLOSE (lun)

  RETURN
END SUBROUTINE readbsq
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE LAE_LLTOXY                 ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE lae_lltoxy(m,n,ctrlat,ctrlon,lat,lon,x,y)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Determine x,y from latitude, longitude coordinates in
!  Lambert Azimuthal equal area map projection
!
!-----------------------------------------------------------------------
!  AUTHOR: Leilei Wang , Vince Wong and Ming Xue
!
!   3/25/97.
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  REAL :: lat,lon
  REAL :: x,y
  INTEGER :: m,n
  REAL :: d2rad,eradius
  PARAMETER (d2rad=3.141592654/180.,                                    &
             eradius = 6371000. )      ! mean earth radius in m

  REAL :: ctrlat,ctrlon                   ! Lat. & Lon. of projection
                                          ! center
  REAL :: const1,k

  const1=SIN(ctrlat*d2rad)*SIN(lat*d2rad)+COS(ctrlat*d2rad)*            &
         COS(lat*d2rad)*COS((lon-ctrlon)*d2rad)
  k = SQRT(2.0/(1+const1))
  x = eradius*k*COS(lat*d2rad)*SIN((lon-ctrlon)*d2rad)
  y = eradius*k*( COS(ctrlat*d2rad)*SIN(lat*d2rad)-                     &
      SIN(ctrlat*d2rad)*COS(lat*d2rad)*COS((lon-ctrlon)*d2rad) )

  RETURN
END SUBROUTINE lae_lltoxy
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE LAE_XYTOLL                 ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE lae_xytoll(m,n,x,y,lat,lon)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Determine latitude, longitude coordinates from  the given x,y in
!  Lambert Azimuthal equal area map projection
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Leilei Wang , Vince Wong and Ming Xue
!
!   3/25/97.
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  REAL    :: lat,lon
  REAL    :: x,y
  INTEGER :: m,n
  REAL    :: r, c
  REAL, PARAMETER :: d2rad   = 3.141592654/180.,                        &
                     eradius = 6371000.        ! mean earth radius in m
  REAL, PARAMETER :: ctrlat = 45.0,                                     &
                     ctrlon = -100.0 ! center Lat. & Lon. of projection

  r = SQRT(x**2+y**2)
  c = 2* ASIN( r/(2.0*eradius) )

  lat=ASIN( COS(c)*SIN(ctrlat*d2rad)+y*SIN(c)*COS(ctrlat*d2rad)/r)
  lon=ctrlon*d2rad + ATAN (x*SIN(c)/(r*COS(ctrlat*d2rad)*COS(c)         &
       -y*SIN(ctrlat*d2rad)*SIN(c)) )
  lat = lat/d2rad
  lon = lon/d2rad

  RETURN
END SUBROUTINE lae_xytoll
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                 PROGRAM GVEGFRAC                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE gvegfrac(nx,ny,vfrcdr,glat,glon,vegin1,vegin2,veg)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Prepare the green vegetation fraction array veg for use in the
!  ARPS soil model.  Calling this subroutine will overwrite and
!  previously defined values of veg that are within the data rich
!  regions.  The pre-existing data is generated using the veg
!  type/table conversion.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Dan Weber
!  04/15/98
!
!  MODIFICATIONS:
!
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!  month    month of the forecast initialization time.
!
!  day      day of the initialization time.
!
!  OUTPUT:
!
!  veg      vegetation fraction.
!
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT   NONE

  INCLUDE   'arpssfc.inc'      ! Includes dimension size for data

  INTEGER :: nx                ! Number of model grid points in the
                               ! x direction.
  INTEGER :: ny                ! Number of model grid points in the
                               ! y direction.

  REAL :: glat (nx,ny)       ! Lat. of model grid
  REAL :: glon (nx,ny)       ! Lon. of model grid

  REAL :: veg    (nx,ny)     ! Vegetation fraction in model domain
  REAL :: vegin1 (gvnx,gvny) ! Green Vegetation fraction data.
                             ! for use in reading in the NESDIS
                             ! green vegetation data set. month#1
  REAL :: vegin2 (gvnx,gvny) ! Green Vegetation fraction data.
                             ! for use in reading in the NESDIS
                             ! green vegetation data set. month#2
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j
  REAL    :: tema, temb, temc, temd, t1, t2, temx, temy
  INTEGER :: itema,itemb,m1,m2,ltdir
  CHARACTER (LEN=256) :: nfile1,nfile2 ! for reading in the green vegetation
                                       ! fraction data set.
  CHARACTER (LEN=*)   :: vfrcdr        ! vegetation fraction directory name.
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'      ! Includes month and day information
!
!-----------------------------------------------------------------------
!
!
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

!
!  determine the monthly data sets to use....
!

!  ltdir=LEN(vfrcdr)
!  CALL strlnth( vfrcdr , ltdir)
  ltdir=LEN_TRIM(vfrcdr)
  IF( ltdir  == 0 ) THEN
    vfrcdr = '.'
    ltdir  =1
  END IF

  IF( vfrcdr(ltdir:ltdir) /= '/') THEN
    ltdir=ltdir+1
    vfrcdr(ltdir:ltdir)='/'
  END IF


  IF(month == 2)THEN  ! february...

    IF(day <= 14)THEN
      m1 = month-1
      m2 = month
    ELSE
      m1 = month
      m2 = month+1
    END IF

  ELSE                ! all other months...

    IF(day < 15)THEN
      m1 = month-1
      m2 = month
    ELSE
      m1 = month
      m2 = month+1
    END IF

  END IF               ! end of month definition if block

  IF (m1 == 0) m1 = 12
  IF (m2 == 13) m2 = 1

!
!  approximate the green vegetation data set for a particular day
!  from monthly averaged data sets.  (not all possible combinations
!  are included in the logic, just the majority)....
!

  IF(month == 2)THEN   !  28 day months...
    IF(day == 29)THEN  !  leap year...
      t1 = 0.5
      t2 = 1.0-t1
    ELSE IF(day > 14)THEN
      t1 = 1.0-FLOAT(day-14)/29.
      t2 = 1.0-t1
    ELSE IF(day <= 14)THEN
      t1 = FLOAT(day+16)/30.
      t2 = 1.0-t1
    END IF
  ELSE IF(month == 4.OR.month == 6.OR.month == 9.OR.                    &
           month == 11)THEN       ! 30 day months...
    IF(day > 15)THEN
      t1 = 1.0-FLOAT(day-15)/30.
      t2 = 1.0-t1
    ELSE IF(day <= 15)THEN
      t2 = FLOAT(day+16)/31.
      t1 = 1.0-t2
    END IF

  ELSE                 !  31 day months....
    IF(day > 15)THEN
      t1 = 1.0-FLOAT(day-15)/31.
      t2 = 1.0-t1
    ELSE IF(day <= 15)THEN
      t2 = FLOAT(day+15)/30.
      t1 = 1.0-t2
    END IF
  END IF    ! end of month if block.....
  PRINT *,' t1 and t1 in gvegfract are ',t1,t2
!
!  read in the green vegetation fraction (range 0.0 to 1.0)
!

!  finish writing the data set filename...

  nfile1(1:ltdir) = vfrcdr(1:ltdir)
  WRITE(nfile1(ltdir+1:ltdir+2),'(i2.2)') m1
  WRITE(nfile1(ltdir+3:ltdir+12),'(a10)') '.gvegfract'

  nfile2(1:ltdir) = vfrcdr(1:ltdir)
  WRITE(nfile2(ltdir+1:ltdir+2),'(i2.2)') m2
  WRITE(nfile2(ltdir+3:ltdir+12),'(a10)') '.gvegfract'

  WRITE (*,*) "Opening green veg file 1: ",nfile1(1:ltdir+12)," ..."
  OPEN (12,FILE=nfile1(1:ltdir+12),STATUS='old',FORM='formatted')
  DO j=1,gvny
    DO i=1,gvnx
      READ (12,'(i3)') itema
      vegin1(i,j) = FLOAT(itema)/100.0
    END DO
  END DO
  CLOSE(12)

  WRITE (*,*) "Opening green veg file 2: ",nfile2(1:ltdir+12)," ..."
  OPEN (13,FILE=nfile2(1:ltdir+12),STATUS='old',FORM='formatted')
  DO j=1,gvny
    DO i=1,gvnx
      READ (13,'(i3)') itema
      vegin2(i,j) = FLOAT(itema)/100.0
    END DO
  END DO
  CLOSE(13)

!
!  Compute the veg field by using a linear interpolation (x,y and t).
!  Note the lower left data point is -55.0N and 0W.
!  The j-0.5 is correct (j is an integer).
!  The first box is at lat=75.016, which is considered the center.
!

  temx = 75.016-0.5*0.144
  temy = -55.00
  DO j=1,ny
    DO i=1,nx
      IF (glat(i,j) <= temx.AND.glat(i,j) >= temy)THEN
                                            ! use green veg
! data base
        itema  = glon(i,j)/0.144
        temb   = glon(i,j)-itema*0.144
        tema   = 1.0-temb

        itemb  =(glat(i,j)+55.00)/0.144
        temd   = glat(i,j)-(-55.00+itemb*0.144)
        temc   = 1.0-temd

        veg(i,j) = t1*(temc*(tema*vegin1(itema,itemb)+                  &
                             temb*vegin1(itema+1,itemb))+               &
                       temd*(tema*vegin1(itema,itemb+1)+                &
                             temb*vegin1(itema+1,itemb+1)))             &
                 + t2*(temc*(tema*vegin2(itema,itemb)+                  &
                             temb*vegin2(itema+1,itemb))+               &
                       temd*(tema*vegin2(itema,itemb+1)+                &
                             temb*vegin2(itema+1,itemb+1)))

      END IF
    END DO
  END DO
  PRINT *,'interpolated veg fraction at ',2,2,' is ',veg(2,2)

  RETURN
END SUBROUTINE gvegfrac

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!TINA - the next two subroutines were added to read in the 30second USGS
!       landuse datasets - modified from the code of CC Lam (Queenie)
!       in Hong Kong
!Andreas - included 30s option by Queenie
!
!Procedures: GET_30s
!            MAPTY30s
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE GET_30S                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE GET_30S(dtype,filedir,latbgn,latend,lonbgn,lonend,           &
                   xdims,ydims,resl,dtlon,dtlat,dtarr)
!
!-----------------------------------------------------------------------
!
! PURPOSE:
!
!   Read in 30s data.
!
!-----------------------------------------------------------------------
!
! AUTHOR:
!   Original copy from Tina Chow.
!
! MODIFICATION:
!
! 07/06/2005 (Y. Wang)
! Code was improved to handle general grid.
!
!-----------------------------------------------------------------------
!
! INPUT:
!
!     dtype    Data flag: 1 for soil type data
!                         2 for vegetation type data
!     filenm   File name of usrface data set
!     xdims    X-dimsion size of the data set
!     ydims    Y-dimsion size of the data set
!
! OUTPUT:
!
!     dtlat    Latitude  values of data grid points
!     dtlon    Longitude values of data grid points
!
!     dtarr    Integer data array
!
!-----------------------------------------------------------------------
!
! Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER,   INTENT(IN) :: dtype        ! Data flag: 1 for soil; 2 for veg
  CHARACTER(LEN=*), INTENT(IN) :: filedir   ! Directory name which contains surface data set
  INTEGER,   INTENT(IN) :: latbgn,latend ! in degree
  INTEGER,   INTENT(IN) :: lonbgn,lonend ! in degree

  INTEGER,   INTENT(IN) :: xdims        ! Columns of the data set
  INTEGER,   INTENT(IN) :: ydims        ! Rows of the data set

  REAL,      INTENT(OUT) :: dtlat(ydims) ! Latitude  values of data grid points
  REAL,      INTENT(OUT) :: dtlon(xdims) ! Longitude values of data grid points
  REAL,      INTENT(OUT) :: resl         ! Resolution of data

  INTEGER,   INTENT(OUT) :: dtarr(xdims,ydims)  ! Data array

!
!-----------------------------------------------------------------------
!
!     Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER          :: lenfl        ! Length of filenm
  INTEGER          :: flunit       ! IO unit of filenm

  INTEGER          :: i,j,r,c
  INTEGER          :: ro,co,ii,jj
  INTEGER          :: clon
  CHARACTER(LEN=8)   :: tname       ! tile name
  CHARACTER(LEN=256) :: flname
  INTEGER*1          :: temp(1200,1200)

  CHARACTER(LEN=1) :: prefix(2) = (/'O','V'/)
  CHARACTER(LEN=1) :: latfix, lonfix
!
!-----------------------------------------------------------------------
!
! Include files:
!
!-----------------------------------------------------------------------
!
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!-----------------------------------------------------------------------
!
! Do some checking
!
!-----------------------------------------------------------------------
  IF (dtype /= 1 .AND. dtype /= 2) THEN
    WRITE(6,'(2a,I2)') 'ERROR in get_30s. dtype must be either ',       &
                       '1 for soil or 2 for veg. dtype = ',dtype
    WRITE(6,'(a)') 'Program stopping ... '
    STOP
  END IF

  IF (MOD(latbgn,10) /= 0 .OR. MOD(latend,10) /= 0) THEN
    WRITE(6,'(2a,I3,a,I3)') 'ERROR in get_30s. latbgn & latend must ',  &
                    'be in 10 degree increasement. latbgn = ',latbgn,   &
                    ' latend = ', latend
    WRITE(6,'(a)') 'Program stopping ... '
    STOP
  END IF

  IF (MOD(lonbgn,10) /= 0 .OR. MOD(lonend,10) /= 0) THEN
    WRITE(6,'(2a,I3,a,I3)') 'ERROR in get_30s. lonbgn & lonend must ',  &
                    'be in 10 degree increasement. lonbgn = ',lonbgn,   &
                    ' lonend = ', lonend
    WRITE(6,'(a)') 'Program stopping ... '
    STOP
  END IF

!-----------------------------------------------------------------------
!
! Read file and fill array
!
!-----------------------------------------------------------------------

  lenfl = LEN( filedir )
  CALL strlnth( filedir, lenfl )

  CALL getunit( flunit )
  DO r = latbgn,latend, 10
    DO c =  lonbgn, lonend, 10
      clon = c
      IF (c >= 180) clon = c - 360  ! map 180-360 to -180 -- -10

      latfix = 'N'
      lonfix = 'E'
      IF (r < 0)    latfix = 'S'
      IF (clon < 0) lonfix = 'W'

      WRITE(tname,'(A1,I2.2,A1,I3.3,A1)')                              &
                          prefix(dtype),abs(r),latfix,abs(clon),lonfix

      WRITE(flname,'(3a)') filedir(1:lenfl),'/',tname

      WRITE(6,'(2A)') 'ARPSSFC reading ',TRIM(flname)

      OPEN (UNIT=flunit, FILE = TRIM(flname), ACTION='READ',            &
                         ACCESS = 'DIRECT', RECL=1200*1200)

      READ (flunit,REC=1) temp

      ro = (r-latbgn)/10
      co = (c-lonbgn)/10
      DO i = 1, 1200
         DO j = 1, 1200
            ii = i+co*1200
            jj = (1200-j+1)+ro*1200
            dtarr(ii,jj) = temp(i,j)
         END DO
      END DO

      CLOSE ( flunit )

    END DO
  END DO
  CALL retunit ( flunit )

  resl = 1./120.
  DO i = 1, xdims
    dtlon(i) = resl*(float(i-1)) + lonbgn
  END DO
  dtlon(:) = MOD(dtlon(:)+360.,360.)

  DO j = 1, ydims
    dtlat(j) = resl*(float(j-1)) + latbgn
  END DO

  RETURN
END SUBROUTINE get_30s
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE MAPTY30S                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE mapty30s(dtype,xdims,ydims,nx,ny,resl,dtlon,dtlat,glon,glat, &
                    srcdat,tartyp,item1)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Transit the 30s USGS landuse data into the simple vegetation type
!  categories and fill them into the model domain by chosing the values
!  at the nearest data grid points.
!
!-----------------------------------------------------------------------
!
!  AUTHOR:
!    Original copy was obtained from Tina Chow.
!
!  MODIFICATIONS:
!
!  07/06/2005 (Y. Wang)
!  Interpolation scheme was improved both for efficiency and correctness.
!
!-----------------------------------------------------------------------
!
! INPUT:
!
!     nx       Number of grid points in the x-direction
!     ny       Number of grid points in the y-direction
!     xdims    X-dimsion size of the data set
!     ydims    Y-dimsion size of the data set
!
!     glon     Longitude values of model grid points
!     glat     Latitude  values of model grid points
!
!     dtlon    Longitude values of data grid points
!     dtlat    Latitude  values of data grid points
!
!     vtyp     Vegetation data array
!     dtype    Data type. 1:soil 2:veg
!
! OUTPUT:
!
!     tartyp   Vegetation type array
!
!-----------------------------------------------------------------------
!
!     Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER   dtype        ! data type

  INTEGER   nx           ! Number of grid points in the x-direction
  INTEGER   ny           ! Number of grid points in the y-direction
  INTEGER   xdims        ! X-dimsion size of the data set
  INTEGER   ydims        ! Y-dimsion size of the data set

  REAL      resl         ! Data resolution

  REAL      glon(nx,ny)  ! Longitude values of model grid points
  REAL      glat(nx,ny)  ! Latitude  values of model grid points

  REAL      dtlon(xdims) ! Longitude values of data grid points
  REAL      dtlat(ydims) ! Latitude  values of data grid points

  INTEGER   srcdat(xdims,xdims)  ! veg data array, input

  INTEGER   tartyp(nx,ny)        ! veg type array, output

  INTEGER   item1(xdims,ydims)   ! Temporary array, work array
!
!-----------------------------------------------------------------------
!
!     Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j, ii, jj, ix, jy
  INTEGER :: ixmin,ixmax,jymin,jymax
  INTEGER :: latmin, latmax, lonmin, lonmax
  REAL    :: r,r1
  LOGICAL :: lonneg
!
!-----------------------------------------------------------------------
!
!     Include files:
!
!-----------------------------------------------------------------------
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!     Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

!-----------------------------------------------------------------------
!
! rematch soil
!
!                     Soil type match table
!
!                   Source index     ARPS index
!                  ==============    ===========
!                         1               1
!                         2               2
!                         3               3
!                         4               4
!                         5               4
!                         6               5
!                         7               6
!                         8               7
!                         9               8
!                        10               9
!                        11              10
!                        12              11
!                        13               5
!                        14              13
!                        15              12
!                        16              12
!
!-----------------------------------------------------------------------

  IF (dtype == 1) THEN

    DO j = 1, ydims
      DO i = 1, xdims

        SELECT CASE (srcdat(i,j))
        CASE (1:4)
          item1(i,j) = srcdat(i,j)
        CASE (5:12, 14)     ! srcdat(i,j) = 4 & 5 to the same item1(i,j)
          item1(i,j) = srcdat(i,j)-1
        CASE (13)
          item1(i,j) = 5
        CASE (15, 16)
          item1(i,j) = 12
        CASE DEFAULT

          WRITE (6, '(a/a,i4,a,i4,a,i2/a)')                 &
             'The value for vegetation type is invalid.',   &
             'srcdat(', i, ',', j, ') = ', srcdat(i,j),     &
             'Program stops in subroutine mapty30s.'

          STOP

        END SELECT

      END DO
    END DO

  END IF

!-----------------------------------------------------------------------
!
! rematch vegetation
!
!                     Vegetation match table
!
!                      USGS
!                   Source index     ARPS index
!                  ==============    ===========
!                         1              13
!                         2              10
!                         3              10
!                         4              10
!                         5               3
!                         6               5
!                         7               3
!                         8              12
!                         9               4
!                        10               4
!                        11               6
!                        12               6
!                        13               7
!                        14               7
!                        15               6
!                        16              14
!                        17              11
!                        18              11
!                        19              13
!                        20               2
!                        21               2
!                        22               2
!                        23               1
!                        24               9
!                        25
!
!-----------------------------------------------------------------------

  IF (dtype == 2) THEN

    DO j = 1, ydims
      DO i = 1, xdims

        SELECT CASE (srcdat(i,j))
        CASE (23)
          item1(i,j) = 1
        CASE (20:22)
          item1(i,j) = 2
        CASE (5,7)
          item1(i,j) = 3
        CASE (9,10)
          item1(i,j) = 4
        CASE (6)
          item1(i,j) = 5
        CASE (11,12,15)
          item1(i,j) = 6
        CASE (13,14)
          item1(i,j) = 7
        CASE (24)
          item1(i,j) = 9
        CASE (2:4)
          item1(i,j) = 10
        CASE (17,18)
          item1(i,j) = 11
        CASE (8)
          item1(i,j) = 12
        CASE (1, 19)
          item1(i,j) = 13
        CASE (16)
          item1(i,j) = 14
        CASE DEFAULT

          WRITE (6, '(a/a,i4,a,i4,a,i2/a)')                  &
              'The value for vegetation type is invalid.',   &
              'srcdat(', i, ',', j, ') = ', srcdat(i,j),     &
              'Program stops in subroutine mapty30s.'
          STOP

        END SELECT

      END DO
    END DO

  END IF

!
!-----------------------------------------------------------------------
!
!  Fill the vegetation type into the model grid by chosing the values
!  at the nearest data grid points.
!
!-----------------------------------------------------------------------
!
!  ixmin = min( xdims-1,                                                 &
!                  int( ( glon(1,1) -dtlon(1) ) / resl ) + 1,            &
!                  int( ( glon(1,ny)-dtlon(1) ) / resl ) + 1 )
!
!  ixmin = max( 1, ixmin - 5 )
!
!  jymin = min( ydims-1,                                                 &
!                  int( ( glat(1,1) -dtlat(1) ) / resl ) + 1,            &
!                  int( ( glat(nx,1)-dtlat(1) ) / resl ) + 1 )
!
!  jymin = max( 1, jymin - 5 )
!
!  ixmax = min( xdims,                                                   &
!                  nint( ( glon(nx,1) -dtlon(1) ) / resl ) + 1,          &
!                  nint( ( glon(nx,ny)-dtlon(1) ) / resl ) + 1 )
!
!  ixmax = min( xdims, ixmax + 5 )
!
!  jymax = min( ydims,                                                   &
!                  nint( ( glat(1,ny) -dtlat(1) ) / resl ) + 1,          &
!                  nint( ( glat(nx,ny)-dtlat(1) ) / resl ) + 1 )
!
!  jymax = min( ydims, jymax + 5 )

  lonneg = .FALSE.   ! longitude changes from 0 to 360
  IF (dtlon(xdims) < dtlon(1) .AND. dtlon(1) >= 180) THEN ! Meridian in between
    lonneg = .TRUE.
    DO i = 1,xdims
      IF (dtlon(i) >= 180) dtlon(i) = dtlon(i) - 360.
    END DO

    DO j = 1,ny
      DO i = 1,nx
        IF (glon(i,j) >= 180) glon(i,j) = glon(i,j) - 360.
      END DO
    END DO
  END IF

  DO j = 1, ny
    DO i = 1, nx

      ! John Manobianco, Ph.D. Vice President/Operations of MESO, Inc.
      ! suggested a more efficient algorithm.

      ii = NINT((abs(dtlon(1)-glon(i,j))/resl)+1.)
      jj = NINT((abs(dtlat(1)-glat(i,j))/resl)+1.)

      !!latmin = FLOOR  (glat(i,j))-1
      !!latmax = CEILING(glat(i,j))+1
      !!IF (latmax > dtlat(ydims)) latmax = dtlat(ydims)
      !!IF (latmin < dtlat(1)    ) latmin = dtlat(1)
      !!
      !!lonmin = FLOOR  (glon(i,j))-1
      !!lonmax = CEILING(glon(i,j))+1
      !!IF (lonmax > dtlon(xdims)) lonmax = dtlon(xdims)
      !!IF (lonmin < dtlon(1)    ) lonmin = dtlon(1)
      !!
      !!jymin = nint( (latmin-dtlat(1))/ resl ) + 1
      !!jymax = nint( (latmax-dtlat(1))/ resl )
      !!ixmin = nint( (lonmin-dtlon(1))/ resl ) + 1
      !!ixmax = nint( (lonmax-dtlon(1))/ resl )
      !!
      !!r1 = 1.0e20
      !!ii = 1
      !!jj = 1
      !!DO jy = jymin, jymax          ! nearest neighbor interpolation
      !!  DO ix = ixmin, ixmax
      !!    r = ( dtlon(ix) - glon(i,j) ) ** 2            &
      !!       +( dtlat(jy) - glat(i,j) ) ** 2
      !!    IF ( r < r1 ) THEN
      !!      r1 = r
      !!      ii = ix
      !!      jj = jy
      !!    ENDIF
      !!  END DO
      !!END DO

      tartyp(i,j) = item1(ii,jj)

    END DO
  END DO

  IF (lonneg) THEN  ! recover back to 0 -- 360
    dtlon(:)  = MOD(dtlon(:)+360,   360.)
    glon(:,:) = MOD(glon(:,:)+360., 360.)
  END IF

  RETURN
END SUBROUTINE mapty30s
