
!########################################################################
!########################################################################
!#########                                                      #########
!#########               SUBROUTINE setrastergrid               #########
!#########                                                      #########
!#########                     Developed by                     #########
!#########     Center for Analysis and Prediction of Storms     #########
!#########                University of Oklahoma                #########
!#########                                                      #########
!########################################################################
!########################################################################

SUBROUTINE setrastergrid(rasterdim,radar_lat,radar_lon,rasterdx,        &
                         rasterlat,rasterlon,rasterx,rastery,           &
                         rasterdata,missing)

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! Calculate latitudes/longitudes and Cartesian coordinates (relative to
! radar) of NIDS raster grid points surrounding a radar.  Also, resets
! raster data to missing if raster point is beyond radar range.
!
!-----------------------------------------------------------------------
!
! AUTHOR:  Eric Kemp, 5 April 2001
!
! MODIFICATION HISTORY:
!
! Eric Kemp, 1 August 2001
! Changed rasterx and rastery to automatic arrays, and added USE
! statement to module n2aconst.
!
! Eric Kemp, 9 August 2001
! Reversed earlier changes -- rasterx and rastery are now arguments
! with INTENT(OUT).  Also, added rasterdata array and missing flag
! as arguments.
!
! Keith Brewster, 4 April 2010
! Modification to accomodate conversion of n2aconst to nidscst.inc
!
!-----------------------------------------------------------------------
!
! Force explicit declarations
!
!----------------------------------------------------------------------- 

  IMPLICIT NONE

!-----------------------------------------------------------------------
! 
! Declare arguments
! 
!-----------------------------------------------------------------------
 
  INTEGER, INTENT(IN) :: rasterdim ! Dimension of raster grid
  REAL,INTENT(IN) :: radar_lat     ! Radar latitude (positive deg north)
  REAL,INTENT(IN) :: radar_lon     ! Radar longitude (positive deg east)
  REAL,INTENT(IN) :: rasterdx      ! Raster grid spacing (m)

  REAL,INTENT(OUT) :: rasterlat(rasterdim,rasterdim) ! Latitude of raster
                                                     ! grid points
                                                     ! (positive deg
                                                     !  north)
  REAL,INTENT(OUT) :: rasterlon(rasterdim,rasterdim) ! Longitude of raster
                                                     ! grid points
                                                     ! (positive deg east)
  REAL,INTENT(OUT) :: rasterx(rasterdim,rasterdim)   ! x-coordinate of
                                                     ! grid points
  REAL,INTENT(OUT) :: rastery(rasterdim,rasterdim)   ! y-coordinate of
                                                     ! grid points
  REAL,INTENT(INOUT) :: rasterdata(rasterdim,rasterdim) ! Raster data
                                                        ! array.
  REAL,INTENT(IN) :: missing                            ! Missing value
                                                        ! flag.
                  
!-----------------------------------------------------------------------
!
! Internal variables
!
!-----------------------------------------------------------------------

  REAL :: dist,head

  INTEGER :: istatus
  INTEGER :: i,j

!-----------------------------------------------------------------------
!
! Include file
!
!-----------------------------------------------------------------------

  INCLUDE 'nidscst.inc'

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!-----------------------------------------------------------------------
!
! Construct raster Cartesian grid centered on and relative to the radar.
!
!-----------------------------------------------------------------------
  
  DO j = 1,rasterdim ! Column
    DO i = 1,rasterdim ! Row
      rasterx(i,j) = rasterdx*(REAL(i) - REAL(rasterdim*0.5) - 0.5) ! m  
      rastery(i,j) = rasterdx*(REAL(j) - REAL(rasterdim*0.5) - 0.5) ! m
    END DO ! Row
  END DO ! Column

!-----------------------------------------------------------------------
!
! Calculate latitudes and longitudes of each raster point by following
! the great circle path from the radar (lat and lon) to the raster
! point.  Also, set rasterdata to missing if it is beyond 230 km from 
! the radar.
!
!-----------------------------------------------------------------------

  DO j = 1,rasterdim ! Column
    DO i = 1,rasterdim ! Row
      dist = SQRT( (rasterx(i,j)*rasterx(i,j)) +                        &
                   (rastery(i,j)*rastery(i,j)) )

      IF (dist > 230000.) THEN
        rasterdata(i,j) = missing
      END IF

      IF (rasterx(i,j) == 0 .AND. rastery(i,j) == 0) THEN
        head = 0.
      ELSE
        head = ATAN2(rasterx(i,j),rastery(i,j))*rad2deg
      END IF
      CALL gcircle(radar_lat,radar_lon,head,dist,                       &
                   rasterlat(i,j),rasterlon(i,j))
    END DO ! Row
  END DO ! Column
      
  RETURN
END SUBROUTINE setrastergrid

!########################################################################
!########################################################################
!#########                                                      #########
!#########                SUBROUTINE remapraster                #########
!#########                                                      #########
!#########                     Developed by                     #########
!#########     Center for Analysis and Prediction of Storms     #########
!#########                University of Oklahoma                #########
!#########                                                      #########
!########################################################################
!########################################################################
    
SUBROUTINE remapraster(nx,ny,xs2d,ys2d,dx,rasterchek,datamiss,xrad,yrad,&
                       remapopt,radius,rnx,rny,rasterdata,              &
                       radx,rady,remdata,istatus)
  
!-----------------------------------------------------------------------
!
! PURPOSE:
!
! Remap 2d raster data to ARPS grid.  Several remapping methods are 
! available: 
!
! 1.  Bi-quadratic interpolation.  IMPORTANT:  Raster and ARPS
!     Cartesian coordinates must be in radar projection space.
!
! 2.  Bi-quadratic regression (requires radius of influence).
!
!     a.  If at least eight radar data points are matched to an ARPS
!         scalar point, a bi-quadratic regression equation of the
!         form V = A + B*x + C*y + D*x*x + E*y*y is solved for
!         to determine the data value at the scalar point.  A check
!         is made to make sure that the solution is within the
!         maximum and minimum radar point values; if this criteria
!         is not met, the minimum radar point value is substituted
!         for the regression solution.
!
!     b.  If between one and eight radar data points are matched to
!         an ARPS scalar point, the radar point values are simply
!         averaged.
!
! 3.  Using maximum raster value (requires radius of influence).
!     
!-----------------------------------------------------------------------
!         
! AUTHOR:  Eric Kemp, 27 April 2001
!          Based on subroutine 'remap' by Keith Brewster and
!          interpolation procedure used in program arpsextsnd.
!
! MODIFICATION HISTORY:
!
! Eric Kemp, 1 August 2001
! Reorganized code and updated documentation.
!          
! Eric Kemp, 9 August 2001
! Changed ARPS scalar grid point arguments to 2-dimensions.
!
!-----------------------------------------------------------------------
!
! Force explicit declarations
!
!-----------------------------------------------------------------------

  IMPLICIT NONE
      
!-----------------------------------------------------------------------
!             
! Declare arguments
!     
!-----------------------------------------------------------------------
          
  INTEGER, INTENT(IN) :: nx,ny     ! Dimensions of ARPS grid
  INTEGER, INTENT(IN) :: rnx, rny  ! Dimensions of raster data
  INTEGER, INTENT(IN) :: remapopt  ! Remapping option.
                                   ! 1 = Bi-quadratic interpolation
                                   ! 2 = Bi-quadratic regreesion/mean
                                   ! 3 = Maximum raster value
  REAL, INTENT(IN) :: radius       ! Radius of influence used with
                                   ! remapopt = 2 or 3.
  REAL, INTENT(IN) :: dx           ! ARPS horizontal grid resolution (m)
  REAL, INTENT(IN) :: rasterchek   ! Flag for missing raster data.
  REAL, INTENT(IN) :: datamiss     ! Flag for missing remapped data.
  REAL, INTENT(IN) :: xs2d(nx,ny),ys2d(nx,ny) ! Scalar coordinates of
                                              ! ARPS grid points (m).
  REAL, INTENT(IN) :: xrad,yrad    ! Scalar coordinates of radar (m).
  REAL, INTENT(IN) :: rasterdata(rnx,rny) ! Raster data (VIL, echo top, 
                                          ! etc.)
  
  REAL, INTENT(IN) :: radx(rnx,rny)   ! x coordinate of raster point.
  REAL, INTENT(IN) :: rady(rnx,rny)   ! y coordinate of raster point.
  REAL, INTENT(OUT) :: remdata(nx,ny) ! Remapped radar data.
  INTEGER, INTENT(OUT) :: istatus     ! Status variable
                                   
!-----------------------------------------------------------------------
!
!  Declare functions
!
!-----------------------------------------------------------------------
  
   REAL :: pntint2d2d
  
!-----------------------------------------------------------------------
!
!  Interpolation arrays
!
!-----------------------------------------------------------------------
  
  INTEGER, ALLOCATABLE :: ipt(:,:)
  INTEGER, ALLOCATABLE :: jpt(:,:) 

  REAL, ALLOCATABLE :: dxfld(:,:)
  REAL, ALLOCATABLE :: dyfld(:,:)
  REAL, ALLOCATABLE :: rdxfld(:,:)
  REAL, ALLOCATABLE :: rdyfld(:,:)
  REAL, ALLOCATABLE :: slopey(:,:)
  REAL, ALLOCATABLE :: alphay(:,:)
  REAL, ALLOCATABLE :: betay(:,:)

!-----------------------------------------------------------------------
!  
! Declare internal variables for least-squares curve fitting.
!
!-----------------------------------------------------------------------
  
  INTEGER, PARAMETER :: n = 5   ! Quadratic

  REAL, ALLOCATABLE :: araster(:,:)
  REAL, ALLOCATABLE :: rhsraster(:)
  REAL, ALLOCATABLE :: sol(:)
  REAL, ALLOCATABLE :: work(:,:)
  REAL, ALLOCATABLE :: work1d(:)

!-----------------------------------------------------------------------
!
! Other internal variables
!
!-----------------------------------------------------------------------
   
  REAL :: ddx,ddy
  REAL :: dxthr
  LOGICAL :: rastergood
  REAL :: rastermax, rastermin, rastermean
  INTEGER :: i,j,ii,jj
  REAL :: delx, dely
  REAL :: slrange
  INTEGER :: knt
  
  REAL, PARAMETER :: eps = 1.0E-25
  
  INTEGER, PARAMETER :: iorder = 2 ! Bi-quadratic interpolation

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! 
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  IF (remapopt < 1 .OR. remapopt > 3) THEN
    WRITE(6,*)'remapraster:  ERROR -- Invalid remap option.'
    WRITE(6,*)'remapraster:  remapopt = ',remapopt
    WRITE(6,*)'remapraster:  Aborting...'
    istatus = 1
    RETURN
  END IF  

  remdata(:,:) = datamiss
  
!-----------------------------------------------------------------------
!
! Allocate necessary arrays and coefficients.
!
!-----------------------------------------------------------------------
  
  IF (remapopt == 1) THEN ! Bi-quadratic interpolation

    ALLOCATE(ipt(nx,ny),STAT=istatus)
    ALLOCATE(jpt(nx,ny),STAT=istatus)
    ALLOCATE(dxfld(rnx,rny),STAT=istatus)
    ALLOCATE(dyfld(rny,rny),STAT=istatus)
    ALLOCATE(rdxfld(rnx,rny),STAT=istatus)
    ALLOCATE(rdyfld(rnx,rny),STAT=istatus)
    ALLOCATE(slopey(rnx,rny),STAT=istatus)
    ALLOCATE(alphay(rnx,rny),STAT=istatus)
    ALLOCATE(betay(rnx,rny),STAT=istatus)
    
    CALL setijloc2(rnx,rny,nx,ny,ipt,jpt,radx,rady,xs2d,ys2d)
  
    CALL setdxdy2d(rnx,rny,                                              &
               1,rnx,1,rny,                                              &
               radx,rady,dxfld,dyfld,rdxfld,rdyfld)
    
    CALL setdrvy2d(rnx,rny,1,                                            &
               1,rnx,1,rny,1,1,                                          &
               dyfld,rdyfld,rasterdata,                                  &
               datamiss,slopey,alphay,betay)

  ELSE

    ALLOCATE(araster(n,n),STAT=istatus)
    ALLOCATE(rhsraster(n),STAT=istatus)
    ALLOCATE(sol(n),STAT=istatus)
    ALLOCATE(work(n,n+1),STAT=istatus)
    ALLOCATE(work1d(n+1),STAT=istatus)

  END IF ! IF (remapopt == 1)
        
  DO j=1,ny
    DO i=1,nx

!-----------------------------------------------------------------------
!
!     Proceed with the bi-quadratic interpolation.
!
!-----------------------------------------------------------------------
               
      IF (remapopt == 1) THEN
        IF (ipt(i,j) <= rnx-1 .AND. ipt(i,j) >= 1 .AND.                 &
            jpt(i,j) <= rny-1 .AND. jpt(i,j) >= 1)  THEN
  
          remdata(i,j)=pntint2d2d(rnx,rny,                              &
               1,rnx-1,1,rny-1,                                         &
               iorder,radx,rady,xs2d(i,j),ys2d(i,j),                    &
               ipt(i,j),jpt(i,j),rasterdata,                            &
               dxfld,dyfld,rdxfld,rdyfld,                               &
               datamiss,slopey,alphay,betay)

          IF (remdata(i,j) /= datamiss) THEN
            remdata(i,j) = MAX(0.,remdata(i,j))
          END IF
        END IF
      ELSE ! remapopt /= 1
 
!-----------------------------------------------------------------------
!
!       Using an option that requires the radius of influence.  First,
!       determine the distance between the current ARPS grid point and 
!       the radar.              
!      
!----------------------------------------------------------------------- 

        delx=xs2d(i,j)-xrad
        dely=ys2d(i,j)-yrad
        slrange=SQRT(delx*delx + dely*dely)
    
        IF (slrange > 0.) THEN
               
          rastermax = -9999.
          rastermin = 9999.
        
          araster(:,:) = 0.
          rhsraster(:) = 0.

!-----------------------------------------------------------------------
!       
!         Loop through each raster data point, searching for valid
!         data.
!
!-----------------------------------------------------------------------
  
          DO jj = 1,rny
            DO ii = 1,rnx
              rastergood = (rasterdata(ii,jj) > rasterchek)
        
              IF (rastergood) THEN

!-----------------------------------------------------------------------
!       
!               For a "valid" raster data point, find the distance
!               between the data point and the current ARPS scalar
!               grid point.
!
!-----------------------------------------------------------------------

                ddx=radx(ii,jj)-xs2d(i,j)
                ddy=rady(ii,jj)-ys2d(i,j)

!-----------------------------------------------------------------------
!
!               Compare radius of influence with actual horizontal 
!               ddx and ddy) distances between the raster point and the
!               ARPS scalar point.  If actual distance is within the
!               threshold, we'll process the raster point.
!
!-----------------------------------------------------------------------

                dxthr = radius

                IF ( SQRT((ddx*ddx) + (ddy*ddy)) < dxthr ) THEN
                  rastermax=MAX(rastermax,rasterdata(ii,jj))
                  rastermin=MIN(rastermin,rasterdata(ii,jj))

                  IF (remapopt == 2) THEN ! Regression/average option

!-----------------------------------------------------------------------
!
!                   Save raster data and 2-D raster point to grid point
!                   data.  These data will be used to determine the
!                   grid point value using least squares curve fitting
!                   (Data = A + B*ddx + C*ddy + D*ddx*ddx + E*ddy*ddy).
!               
!-----------------------------------------------------------------------

                    rhsraster(1)=rhsraster(1)+rasterdata(ii,jj) 
                    rhsraster(2)=rhsraster(2)+rasterdata(ii,jj)*ddx
                    rhsraster(3)=rhsraster(3)+rasterdata(ii,jj)*ddy
                    rhsraster(4)=rhsraster(4)+rasterdata(ii,jj)*ddx*ddx
                    rhsraster(5)=rhsraster(5)+rasterdata(ii,jj)*ddy*ddy

                    araster(1,1)=araster(1,1)+1.  
                    araster(1,2)=araster(1,2)+ddx
                    araster(1,3)=araster(1,3)+ddy 
                    araster(1,4)=araster(1,4)+(ddx*ddx)
                    araster(1,5)=araster(1,5)+(ddy*ddy)
                    
                    araster(2,1)=araster(2,1)+ddx 
                    araster(2,2)=araster(2,2)+ddx*ddx 
                    araster(2,3)=araster(2,3)+ddx*ddy 
                    araster(2,4)=araster(2,4)+ddx*ddx*ddx
                    araster(2,5)=araster(2,5)+ddx*ddy*ddy

                    araster(3,1)=araster(3,1)+ddy
                    araster(3,2)=araster(3,2)+ddy*ddx
                    araster(3,3)=araster(3,3)+ddy*ddy
                    araster(3,4)=araster(3,4)+ddy*ddx*ddx
                    araster(3,5)=araster(3,5)+ddy*ddy*ddy

                    araster(4,1)=araster(4,1)+ddx*ddx
                    araster(4,2)=araster(4,2)+ddx*ddx*ddx
                    araster(4,3)=araster(4,3)+ddx*ddx*ddy
                    araster(4,4)=araster(4,4)+ddx*ddx*ddx*ddx
                    araster(4,5)=araster(4,5)+ddx*ddx*ddy*ddy
                    
                    araster(5,1)=araster(5,1)+ddy*ddy
                    araster(5,2)=araster(5,2)+ddy*ddy*ddx
                    araster(5,3)=araster(5,3)+ddy*ddy*ddy
                    araster(5,4)=araster(5,4)+ddy*ddy*ddx*ddx
                    araster(5,5)=araster(5,5)+ddy*ddy*ddy*ddy

                  END IF ! remapopt
                END IF ! (ABS(ddx) < dxthr .AND. ABS(ddy) < dxthr)
              END IF ! Is rastergood true?
            END DO ! ii loop
          END DO ! jj loop

!-----------------------------------------------------------------------
!
!         At this point we've gone through all the radar pixels to see
!         which one(s) can be matched to an ARPS grid point.  Now we'll
!         process the data we've saved.
!                   
!-----------------------------------------------------------------------
                    
          IF (remapopt == 3) THEN ! Use maximum value
            remdata(i,j) = rastermax

          ELSE ! remapopt == 2
            knt = NINT(araster(1,1)) ! Number of raster points matched
                                     ! to grid point.

            IF (knt >= 8) THEN
!            IF (knt >= 5) THEN

!-----------------------------------------------------------------------
!         
!             If at least eight raster points have been matched to the
!             grid point, use Gauss-Jordan elimination to calculate the 
!             constants in the least square curve equation.  The
!             resulting grid point value is then checked with the
!             maximum and minimum raster point values.
!
!-----------------------------------------------------------------------
                  
              CALL gjelim(n,araster,rhsraster,sol,work,work1d,eps,      &
                          istatus)
              remdata(i,j) = MIN(rastermax,MAX(rastermin,sol(1)))

!            ELSE IF (knt >= 4) THEN
            ELSE IF (knt >= 1) THEN

!-----------------------------------------------------------------------
!
!             If at least one (but less than eight) raster points have  
!             been matched to the grid point, simply calculate the
!             average raster value and assign it to the grid point.
!             
!-----------------------------------------------------------------------

              rastermean=rhsraster(1)/araster(1,1)
              remdata(i,j)=rastermean

            END IF ! knt
          END IF ! remapopt = 2 vs. 3
        END IF ! slrange > 0.
      END IF ! if remapopt = 1
    END DO ! i loop
  END DO ! j loop

!-----------------------------------------------------------------------
!             
! Deallocate arrays and return.
!
!-----------------------------------------------------------------------
              
  IF (remapopt == 1) THEN
    DEALLOCATE(ipt,jpt,dxfld,dyfld,rdxfld,rdyfld,slopey,alphay,betay,   &
               STAT=istatus)
  ELSE
    DEALLOCATE(araster,rhsraster,sol,work,work1d,STAT=istatus)
  END IF

  istatus = 0
  RETURN
END SUBROUTINE remapraster

!########################################################################
!########################################################################
!#########                                                      #########
!#########                 SUBROUTINE setdxdy2d                 #########
!#########                                                      #########
!#########                     Developed by                     #########
!#########     Center for Analysis and Prediction of Storms     #########
!#########                University of Oklahoma                #########
!#########                                                      #########
!########################################################################
!########################################################################

SUBROUTINE setdxdy2d(nx,ny,                                             &
           ibeg,iend,jbeg,jend,                                         &
           x2d,y2d,dxfld,dyfld,rdxfld,rdyfld)

!-----------------------------------------------------------------------
!
!  PURPOSE:
!    Calculate the local delta-x, delta-y and their inverses.
!    Precalculating these variables speeds up later calculations.  
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster, CAPS, November, 1996
!
!  MODIFICATION HISTORY:
!
!  Eric Kemp, 27 April 2001
!  Modified to allow for two-dimensioned arrays.
!
!  Eric Kemp, 9 August 2001
!  Added INTENT statements.
!
!----------------------------------------------------------------------- 
!
! Force explicit declarations
!
!----------------------------------------------------------------------- 

  IMPLICIT NONE

!----------------------------------------------------------------------- 
!
!  INPUT:
!    nx       Number of model grid points in the x-direction (east/west)
!    ny       Number of model grid points in the y-direction (north/south)
!  
!    ibeg,iend   Range of x index to do interpolation
!    jbeg,jend   Range of y index to do interpolation
!
!    x2d     Array of x-coordinate grid locations (m)
!    y2d     Array of y-coordinate grid locations (m)
!  
!  OUTPUT:
!    dxfld    Vector of delta-x (m) of field to be interpolated
!    dyfld    Vector of delta-y (m) of field to be interpolated
!    rdxfld   Vector of 1./delta-x (1/m) of field to be interpolated
!    rdyfld   Vector of 1./delta-y (1/m) of field to be interpolated
!
!----------------------------------------------------------------------- 

  INTEGER, INTENT(IN) :: nx,ny
  INTEGER, INTENT(IN) :: ibeg,iend
  INTEGER, INTENT(IN) :: jbeg,jend
  REAL, INTENT(IN) :: x2d(nx,ny)
  REAL, INTENT(IN) :: y2d(nx,ny)
  REAL, INTENT(OUT) :: dxfld(nx,ny)
  REAL, INTENT(OUT) :: dyfld(nx,ny)
  REAL, INTENT(OUT) :: rdxfld(nx,ny)
  REAL, INTENT(OUT) :: rdyfld(nx,ny)

!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------
  
  INTEGER :: i,j,istop,jstop

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  
  istop=MIN((iend-1),(nx-1))
  jstop=MIN((jend-1),(ny-1))
 
  DO j=jbeg,jstop
    DO i=ibeg,istop
      dxfld(i,j)=(x2d(i+1,j)-x2d(i,j))
      rdxfld(i,j)=1./(x2d(i+1,j)-x2d(i,j))
      dyfld(i,j)=(y2d(i,j+1)-y2d(i,j))
      rdyfld(i,j)=1./(y2d(i,j+1)-y2d(i,j))
    END DO ! i loop
  END DO ! j loop

  RETURN
END SUBROUTINE setdxdy2d

!########################################################################
!########################################################################
!#########                                                      #########
!#########                 SUBROUTINE setdrvy2d                 #########
!#########                                                      #########
!#########                     Developed by                     #########
!#########     Center for Analysis and Prediction of Storms     #########
!#########                University of Oklahoma                #########
!#########                                                      #########
!########################################################################
!########################################################################

SUBROUTINE setdrvy2d(nx,ny,nz,                                          &
           ibeg,iend,jbeg,jend,kbeg,kend,                               &
           dyfld,rdyfld,var,                                            &
           missing,slopey,alphay,betay)
    
!-----------------------------------------------------------------------
!
! PURPOSE:
!    Calculate the coefficients of interpolating polynomials
!    in the y-direction.
!
!-----------------------------------------------------------------------
!
! AUTHOR: Keith Brewster, CAPS, November, 1996
!
! MODIFICATION HISTORY:
!
! Eric Kemp, 27 April 2001
! Modified to allow for 2D arrays
!
! Eric Kemp, 7 August 2001
! Modified to check for missing data, and added INTENT statements.
!
!----------------------------------------------------------------------- 
!
! Force explicit declarations.
!
!----------------------------------------------------------------------- 

  IMPLICIT NONE

!----------------------------------------------------------------------- 
!
! INPUT:
!    nx       Number of model grid points in the x-direction (east/west)
!    ny       Number of model grid points in the y-direction (north/south)
!    nz       Number of model grid points in the vertical
!    
!    ibeg,iend   Range of x index to do interpolation
!    jbeg,jend   Range of y index to do interpolation
!    kbeg,kend   Range of z index to do interpolation
!
!    dyfld    Vector of delta-y (m) of field to be interpolated
!    rdyfld   Vector of 1./delta-y (1/m) of field to be interpolated
!  
!    var      variable to be interpolated
!    missing  Value of missing data
!  
!    slopey   Piecewise linear df/dy
!    alphay   Coefficient of y-squared term in y quadratic interpolator
!    betay    Coefficient of y term in y quadratic interpolator
!
!-----------------------------------------------------------------------
   
  INTEGER, INTENT(IN) :: nx,ny,nz
  INTEGER, INTENT(IN) :: ibeg,iend,jbeg,jend,kbeg,kend
  REAL, INTENT(IN) :: dyfld(nx,ny)
  REAL, INTENT(IN) :: rdyfld(nx,ny)
  REAL, INTENT(IN) :: var(nx,ny,nz)
  REAL, INTENT(IN) :: missing
  REAL, INTENT(OUT) :: slopey(nx,ny,nz)
  REAL, INTENT(OUT) :: alphay(nx,ny,nz)
  REAL, INTENT(OUT) :: betay(nx,ny,nz)  
    
!-----------------------------------------------------------------------
!    
! Misc. local variables
!
!-----------------------------------------------------------------------

  INTEGER :: i,j,k   
  INTEGER :: jstart,jstop
  REAL :: rtwody

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  jstart=MAX(jbeg,2)
  jstop=MIN((jend-1),(ny-2))
  DO k=kbeg,kend
    DO j=jstart,jstop
      DO i=ibeg,iend 
        IF (var(i,j+1,k) == missing .OR.                                &
            var(i,j,k) == missing) THEN
          slopey(i,j,k) = missing
        ELSE
          slopey(i,j,k)=(var(i,j+1,k)-var(i,j,k))*rdyfld(i,j)
        END IF

        rtwody=1./(dyfld(i,j-1)+dyfld(i,j))

        IF (var(i,j+1,k) == missing .OR.                                &
            var(i,j,k) == missing .OR.                                  &
            var(i,j-1,k) == missing) THEN
          alphay(i,j,k) = missing
        ELSE
          alphay(i,j,k)=((var(i,j+1,k)-var(i,j,k))*rdyfld(i,j) +        &
                 (var(i,j-1,k)-var(i,j,k))*rdyfld(i,j-1))*rtwody
        END IF

        IF (var(i,j+1,k) == missing .OR.                                &
            var(i,j,k) == missing  .OR.                                 &
            alphay(i,j,k) == missing ) THEN
          betay(i,j,k)=missing
        ELSE
          betay(i,j,k)=(var(i,j+1,k)-var(i,j,k))*rdyfld(i,j) -          &
                   dyfld(i,j)*alphay(i,j,k)
        END IF

      END DO
    END DO
  END DO
  RETURN
END SUBROUTINE setdrvy2d

!########################################################################
!########################################################################
!#########                                                      #########
!#########                 FUNCTION pntint2d2d                  #########
!#########                                                      #########
!#########                     Developed by                     #########
!#########     Center for Analysis and Prediction of Storms     #########
!#########                University of Oklahoma                #########
!#########                                                      #########
!########################################################################
!########################################################################

REAL FUNCTION pntint2d2d(vnx,vny,ivbeg,ivend,jvbeg,jvend,iorder,vx,vy,  &
                         xpnt,ypnt,iloc,jloc,var,dxfld,dyfld,rdxfld,    &
                         rdyfld,missing,slopey,alphay,betay)

!-----------------------------------------------------------------------
!
! PURPOSE:
!
!    Interpolate a 2-d field for a single point on that plane.
!
!-----------------------------------------------------------------------
!
! AUTHOR: Keith Brewster, CAPS, November, 1996
!
! MODIFICATION HISTORY:
!   
! Eric Kemp, 27 April 2001
! Modified to allow for 2D arrays
!
! Eric Kemp, 9 August 2001
! Modified to check for missing data, and added INTENT statements.
!
!----------------------------------------------------------------------- 
!
! Force explicit declarations
!
!----------------------------------------------------------------------- 

  IMPLICIT NONE

!----------------------------------------------------------------------- 
!
! INPUT:
!    vnx       Number of model grid points in the x-direction (east/west)
!    vny       Number of model grid points in the y-direction
!              (north/south)
!
!    ivbeg,ivend   Range of x index to use in verification array
!    jvbeg,jvend   Range of y index to use in verification array
!
!    iorder   Interpolation parameter.
!             iorder specifies the order of interpolation
!             1 = bi-linear
!             2 = bi-quadratic
!   
!    vx       x coordinate of verif scalar grid points in physical space
!             (m)
!    vy       y coordinate of verif scalar grid points in physical space
!             (m)
!
!    xpnt     x coordinate (m) of interpolation point
!    ypnt     y coordinate (m) of interpolation point
!
!    iloc     I-index of interpolation point in field to be interpolated
!    jloc     J-index of interpolation point in field to be interpolated
!    dxfld    Vector of delta-x (m) of field to be interpolated 
!    dyfld    Vector of delta-y (m) of field to be interpolated 
!    rdxfld   Vector of 1./delta-x (1/m) of field to be interpolated
!    rdyfld   Vector of 1./delta-y (1/m) of field to be interpolated
!    missing  Value of missing data
!             
!    slopey   Piecewise linear df/dy
!    alphay   Coefficient of y-squared term in y quadratic interpolator
!    betay    Coefficient of y term in y quadratic interpolator
!    
!-----------------------------------------------------------------------

  INTEGER, INTENT(IN) :: vnx,vny
  INTEGER, INTENT(IN) :: ivbeg,ivend,jvbeg,jvend
  INTEGER, INTENT(IN) :: iorder
  REAL, INTENT(IN) :: vx(vnx,vny)
  REAL, INTENT(IN) :: vy(vnx,vny)
  REAL, INTENT(IN) :: xpnt
  REAL, INTENT(IN) :: ypnt
  INTEGER, INTENT(IN) :: iloc
  INTEGER, INTENT(IN) :: jloc
  REAL, INTENT(IN) :: var(vnx,vny)
  REAL, INTENT(IN) :: dxfld(vnx,vny)
  REAL, INTENT(IN) :: dyfld(vnx,vny)
  REAL, INTENT(IN) :: rdxfld(vnx,vny)
  REAL, INTENT(IN) :: rdyfld(vnx,vny)
  REAL, INTENT(IN) :: slopey(vnx,vny)
  REAL, INTENT(IN) :: alphay(vnx,vny)
  REAL, INTENT(IN) :: betay(vnx,vny)
  REAL, INTENT(IN) :: missing

!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------

  INTEGER :: ii,jj   
  REAL :: delx,dely  
  REAL :: alpha,beta,rtwodx
  REAL :: varm1,var00,varp1
  REAL :: varint

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!-----------------------------------------------------------------------
!
! Compute bilinear interpolated value
!
!-----------------------------------------------------------------------

  IF(iorder == 1) THEN
    ii=MIN(MAX(iloc,ivbeg),(ivend-1))
    jj=MIN(MAX(jloc,jvbeg),(jvend-1))
  
    delx=(xpnt-vx(ii,jj))
    dely=(ypnt-vy(ii,jj))

    IF (var(ii,jj) == missing .OR.                                      &
        slopey(ii,jj) == missing .OR.                                   &
        var(ii+1,jj) == missing .OR.                                    &
        slopey(ii+1,jj) == missing) THEN
      varint = missing
    ELSE
      varint=(1.-delx*rdxfld(ii,jj))*                                   &
               (var(ii  ,jj)+slopey(ii  ,jj)*dely)+                     &
               (delx*rdxfld(ii,jj))*                                    &
               (var(ii+1,jj)+slopey(ii+1,jj)*dely)
    END IF

!-----------------------------------------------------------------------
!
! Compute biquadratic
! 
!-----------------------------------------------------------------------

  ELSE

    ii=MIN(MAX(iloc,(ivbeg+1)),(ivend-1))
    jj=MIN(MAX(jloc,(jvbeg+1)),(jvend-1))
    
    delx=(xpnt-vx(ii,jj))
    dely=(ypnt-vy(ii,jj))

!-----------------------------------------------------------------------
!
!   Stencil is ii-1 to ii+1 and jj-1 to jj + 1
!
!   Interpolate in y.
! 
!-----------------------------------------------------------------------

    IF (alphay(ii-1,jj) == missing .OR. betay(ii-1,jj) == missing .OR.  &
        var(ii-1,jj) == missing) THEN
      varm1 = missing
    ELSE
      varm1=(alphay(ii-1,jj)*dely+betay(ii-1,jj))*dely+var(ii-1,jj)
    END IF

    IF (alphay(ii,jj) == missing .OR. betay(ii,jj) == missing .OR.      &
        var(ii,jj) == missing) THEN
      varm1 = missing
    ELSE
      var00=(alphay(ii  ,jj)*dely+betay(ii  ,jj))*dely+var(ii  ,jj)
    END IF

    IF (alphay(ii+1,jj) == missing .OR. betay(ii+1,jj) == missing .OR.  &
        var(ii+1,jj) == missing) THEN
      varm1 = missing
    ELSE
      varp1=(alphay(ii+1,jj)*dely+betay(ii+1,jj))*dely+var(ii+1,jj)
    END IF

!-----------------------------------------------------------------------
!    
!   Interpolate intermediate results in x.
!
!-----------------------------------------------------------------------

    rtwodx=1./(dxfld(ii-1,jj)+dxfld(ii,jj))
    IF (varp1 == missing .OR. var00 == missing .OR.                     &
        varm1 == missing) THEN
      varint = missing
    ELSE
      alpha=((varp1-var00)*rdxfld(ii  ,jj) +                            &
           (varm1-var00)*rdxfld(ii-1,jj))*rtwodx
      beta=(varp1-var00)*rdxfld(ii,jj) -                                &
             dxfld(ii,jj)*alpha
      varint=(alpha*delx+beta)*delx+var00
    END IF       
  END IF
  pntint2d2d=varint
  RETURN
END FUNCTION pntint2d2d

!########################################################################
!########################################################################
!#########                                                      #########
!#########                 SUBROUTINE setijloc2                 #########
!#########                                                      #########
!#########                     Developed by                     #########
!#########     Center for Analysis and Prediction of Storms     #########
!#########                University of Oklahoma                #########
!#########                                                      #########
!########################################################################
!########################################################################

SUBROUTINE setijloc2(fnx,fny,anx,any,iloc,jloc,fx2d,fy2d,ax2d,ay2d)

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! Find i,j indicies in forecast grid of each analysis grid point.
! No indicies are returned if the analysis grid point is outside
! of the forecast grid.
!
! NOTE:  The forecast and analysis grid points are assumed to have 
! the same origin and map projection.
!
!-----------------------------------------------------------------------
!
! AUTHOR:  Eric Kemp, November 1999.
!
! Based on subroutine setijloc
!
! MODIFICATION HISTORY:
! 
! Eric Kemp, 9 August 2001
! Added INTENT statements.
!
!-----------------------------------------------------------------------
!
! Force explicit declarations.
! 
!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------
!
! Arguments
! 
!-----------------------------------------------------------------------

  INTEGER,INTENT(IN) :: fnx,fny,anx,any
  REAL,INTENT(IN) :: fx2d(fnx,fny),fy2d(fnx,fny)
  REAL,INTENT(IN) :: ax2d(anx,any),ay2d(anx,any)
  INTEGER,INTENT(INOUT) :: iloc(anx,any),jloc(anx,any)

!-----------------------------------------------------------------------
! 
! Miscellaneous variables
! 
!-----------------------------------------------------------------------

  INTEGER :: i,j,m,n
  INTEGER :: fimid,fjmid
  REAL :: fxmid,fymid

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!-----------------------------------------------------------------------
!
! Initialize iloc and jloc to -9999 (i.e., outside the forecast grid)   
! 
!-----------------------------------------------------------------------
  
  DO n = 1,any
    DO m = 1,anx
      iloc(m,n) = -9999
      jloc(m,n) = -9999 
    END DO
  END DO

!-----------------------------------------------------------------------
! 
! Determine the scalar coordinates of the middle of the forecast field.
!
!-----------------------------------------------------------------------

  fimid = fnx/2
  fjmid = fny/2
  
  fxmid = fx2d(fimid,fjmid)
  fymid = fy2d(fimid,fjmid)
  
!-----------------------------------------------------------------------
!
! Find the forecast i,j of each analysis grid point.
!
!-----------------------------------------------------------------------

  DO n = 1,any-1
    analysism: DO m = 1,anx-1

      IF (ax2d(m,n) < fxmid) THEN
        IF (ay2d(m,n) < fymid) THEN
          DO j = fjmid,2,-1
            DO i = fimid,2,-1  
              IF (fx2d(i,j) <= ax2d(m,n) .AND.                          &
                  fx2d(i+1,j) > ax2d(m,n) .AND.                         &
                  fy2d(i,j) <= ay2d(m,n) .AND.                          &
                  fy2d(i,j+1) > ay2d(m,n)) THEN
                iloc(m,n) = i
                jloc(m,n) = j
                CYCLE analysism
              END IF
            END DO
          END DO
        ELSE
          DO j = fjmid,fny-1
            DO i = fimid,2,-1
              IF (fx2d(i,j) <= ax2d(m,n) .AND.                          &
                  fx2d(i+1,j) > ax2d(m,n) .AND.                         &
                  fy2d(i,j-1) < ay2d(m,n) .AND.                         &
                  fy2d(i,j) >= ay2d(m,n)) THEN
                iloc(m,n) = i
                jloc(m,n) = j-1
                CYCLE analysism
              END IF
            END DO
          END DO
        END IF
      ELSE
        IF (ay2d(m,n) < fymid) THEN
          DO j = fjmid,2,-1
            DO i = fimid,fnx-1
              IF (fx2d(i-1,j) < ax2d(m,n) .AND.                         &
                  fx2d(i,j) >= ax2d(m,n) .AND.                          &
                  fy2d(i,j) <= ay2d(m,n) .AND.                          &
                  fy2d(i,j+1) > ay2d(m,n)) THEN
                iloc(m,n) = i-1
                jloc(m,n) = j
                CYCLE analysism
              END IF
            END DO
          END DO
        ELSE
          DO j = fjmid,fny-1
            DO i = fimid,fnx-1
              IF (fx2d(i-1,j) < ax2d(m,n) .AND.                         &
                  fx2d(i,j) >= ax2d(m,n) .AND.                          &
                  fy2d(i,j-1) < ay2d(m,n) .AND.                         &
                  fy2d(i,j) >= ay2d(m,n)) THEN
                iloc(m,n) = i-1
                jloc(m,n) = j-1
                CYCLE analysism
              END IF
            END DO
          END DO
        END IF
      END IF
    END DO analysism
  END DO              
END SUBROUTINE setijloc2

!########################################################################
!########################################################################
!#########                                                      #########
!#########                 SUBROUTINE dpalatlon                 #########
!#########                                                      #########
!#########                     Developed by                     #########
!#########     Center for Analysis and Prediction of Storms     #########
!#########                University of Oklahoma                #########
!#########                                                      #########
!########################################################################
!########################################################################

SUBROUTINE dpalatlon(nrow,ncol,radarlat,radarlon,dpalat,dpalon)

!------------------------------------------------------------------------
!
! PURPOSE:
!
! Determines latitude and longitude of each NIDS Digital Precipitation
! Array (DPA) point for a given radar latitude and longitude.  The DPA 
! is on a "1/40 LFM grid", i.e., a polar stereographic grid nested within 
! the old Limited Fine Mesh model grid.  Equations for calculating
! latitudes and longitudes on the DPA grid are taken from Appendix C 
! of Fulton (1998).
!
! IMPORTANT:  The NEXRAD equations assume that the DPA grid is such that 
! i increases to the right and j increases *downward*.  An easy 
! correction is made so that the output arrays follow the ARPS grid, 
! i.e., i increases to the right and j increases *upward*.
!
!------------------------------------------------------------------------
!
! REFERENCE:
!
! Fulton, R. A., 1998:  WSR-88D Polar-to-HRAP Mapping.  Technical
!   Memorandum, Hydrologic Research Laboratory, 33 pp.  Available at:
!   http://hsp.nws.noaa.gov/oh/hrl/papers/wsr88d/hrapmap.pdf
!
!------------------------------------------------------------------------
!
! AUTHOR:
! 
! Eric Kemp, 3 August 2001.
!
! MODIFICATION HISTORY:
!
! Eric Kemp, 9 August 2001.
! Now uses module N2ACONST.
!
! Keith Brewster, 4 April 2010
! Modification to accomodate conversion of n2aconst to nidscst.inc
!
!------------------------------------------------------------------------

  IMPLICIT NONE

!------------------------------------------------------------------------
!
! Arguments
!
!------------------------------------------------------------------------

  INTEGER, INTENT(IN) :: nrow,ncol       ! Dimensions of DPA product
  REAL, INTENT(IN) :: radarlat,radarlon  ! Lat/Lon of radar.
  REAL, INTENT(OUT) :: dpalat(nrow,ncol) ! Lat of DPA points.
  REAL, INTENT(OUT) :: dpalon(nrow,ncol) ! Lon of DPA points.

!------------------------------------------------------------------------
!
! Grid scale factor.
!
!------------------------------------------------------------------------

  REAL, PARAMETER :: Kc = 249.6348607

!------------------------------------------------------------------------
!
! Global grid coordinates of the North Pole in units of 1/4 LFM
! boxes.
!
!------------------------------------------------------------------------

  INTEGER, PARAMETER :: Ip = 433
  INTEGER, PARAMETER :: Jp = 433

!------------------------------------------------------------------------
!
! Global grid coordinates of radar site in units of 1/4 LFM boxes.
!
!------------------------------------------------------------------------

  REAL :: GIs, GJs

!------------------------------------------------------------------------
!
! Global grid box number for 0,0 of 1/40 LFM grid
!
!------------------------------------------------------------------------

  INTEGER :: Is3, Js3

!------------------------------------------------------------------------
!
! Coefficients used to calculate coordinates of box centers in units
! of 1/4 LFM boxes.
!
!------------------------------------------------------------------------

  REAL :: AI3, AJ3
  REAL, PARAMETER :: B3 = 0.1

!------------------------------------------------------------------------
!
! Coordinates of box centers in units of 1/4 LFM boxes.
!
!------------------------------------------------------------------------

  REAL :: CIi, CJj
 
!------------------------------------------------------------------------
!
! Misc. variables.
!
!------------------------------------------------------------------------

  REAL :: sinls, cosls, coslambdaplus, sinlambdaplus
  INTEGER :: i,j

  REAL :: tmplat, tmplon

!-----------------------------------------------------------------------
!
! Include file
!
!-----------------------------------------------------------------------

  INCLUDE 'nidscst.inc'

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! 
! Beginning of executable code...
!                       
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!------------------------------------------------------------------------
!
! Calculate GIs, GJs.  See Appendix C, Sections 30.6.2 and 30.6.3 of
! Fulton (1998).
!
!------------------------------------------------------------------------

  sinls = SIN(radarlat*deg2rad)
  cosls = COS(radarlat*deg2rad)
  coslambdaplus = COS((radarlon + 105.)*deg2rad)
  sinlambdaplus = SIN((radarlon + 105.)*deg2rad)

  GIs = 1. + sinls
  GIs = Kc*cosls*sinlambdaplus/GIs
  GIs = GIs + REAL(Ip)

  GJs = 1 + sinls
  GJs = Kc*cosls*coslambdaplus/GJs
  GJs = GJs + REAL(Jp)

!------------------------------------------------------------------------
!
! Calculate Is3 and Js3.  See Appendix C, Section 30.6.3 of Fulton 
! (1998).  Note that GIs and GJs are multiplied by 10 to make Is3 and 
! Js3 valid on the 1/40 LFM grid.  This corrects an error in the
! NEXRAD documentation.
!
!------------------------------------------------------------------------

  Is3 = INT(10.*GIs) - 66
  Js3 = INT(10.*GJs) - 66
 
!------------------------------------------------------------------------
!
! Calculate AI3 and AJ3 (B3 is already set as a parameter).  See
! Appendix C, Section 30.6.4 of Fulton (1998).  
!
!------------------------------------------------------------------------

  AI3 = -10.*REAL(Ip)
  AI3 = (REAL(Is3) + AI3 + 0.5)/10.

  AJ3 = -10.*REAL(Jp)
  AJ3 =  (REAL(Js3) + AJ3 + 0.5)/10.

!------------------------------------------------------------------------
!
! Begin looping through the DPA points.
!
!------------------------------------------------------------------------

  DO j = 1, ncol
    DO i = 1, nrow

!------------------------------------------------------------------------
!
!     Calculate CIi and CJj for a given DPA point.  See Appendix C, 
!     Section 30.6.4 of Fulton (1998).
!
!------------------------------------------------------------------------

      CIi = AI3 + (B3*REAL(i))
      CJj = AJ3 + (B3*REAL(j))

!------------------------------------------------------------------------
!
!     Now calculate latitude and longitude of a given DPA point.  See
!     Section 30.6.4 of Fulton (1998).  Note that the output array will
!     follow the ARPS grid, i.e., i increases to the right and j 
!     increases *upward*.
!
!------------------------------------------------------------------------

      tmplat = (CIi*CIi) + (CJj*CJj)
      tmplat = SQRT(tmplat)/Kc
      tmplat = 2.*ATAN(tmplat)*rad2deg

      tmplat = 90. - tmplat
      dpalat(i,ncol-j+1) = tmplat

      tmplon = -105. + (ATAN2(CIi,CJj)*rad2deg)
      dpalon(i,ncol-j+1) = tmplon

    END DO ! i loop
  END DO ! j loop

!------------------------------------------------------------------------
!
! The end.
!
!------------------------------------------------------------------------

  RETURN
END SUBROUTINE dpalatlon

!########################################################################
!########################################################################
!#########                                                      #########
!#########               SUBROUTINE getarpsrasterxy             #########
!#########                                                      #########
!#########                     Developed by                     #########
!#########     Center for Analysis and Prediction of Storms     #########
!#########                University of Oklahoma                #########
!#########                                                      #########
!########################################################################
!########################################################################

SUBROUTINE getarpsrasterxy(nx,ny,arpslat,arpslon,radarlat,radarlon,      &
                           arpsrasterx,arpsrastery)

!------------------------------------------------------------------------
!
! PURPOSE:
!
! Calculates Cartesian coordinates of ARPS scalar points relative to
! a radar using spherical geometry.
!
!------------------------------------------------------------------------
!
! AUTHOR:
!
! Eric Kemp, 8 August 2001.
!    
! Keith Brewster, 4 April 2010
! Modification to accomodate conversion of n2aconst to nidscst.inc
!
!------------------------------------------------------------------------

  IMPLICIT NONE

!------------------------------------------------------------------------
!
! Declare arguments
!  
!------------------------------------------------------------------------

  INTEGER, INTENT(IN) :: nx,ny              ! Dimensions of ARPS grid  
  REAL, INTENT(IN) :: arpslat(nx,ny)        ! Latitude of ARPS point
  REAL, INTENT(IN) :: arpslon(nx,ny)        ! Longitude of ARPS point
  REAL, INTENT(IN) :: radarlat,radarlon     ! Lat/Lon of radar.
  REAL, INTENT(OUT) :: arpsrasterx(nx,ny)   ! x-coordinate of ARPS point.
  REAL, INTENT(OUT) :: arpsrastery(nx,ny)   ! y-coordinate of ARPS point.

!------------------------------------------------------------------------
!
! Internal variables.
!
!------------------------------------------------------------------------

  INTEGER :: i,j
  REAL :: headng,dist

!-----------------------------------------------------------------------
!
! Include file
!
!-----------------------------------------------------------------------

  INCLUDE 'nidscst.inc'

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  DO j = 1,ny-1
    DO i = 1,nx-1
      CALL disthead(radarlat,radarlon,arpslat(i,j),arpslon(i,j),         &
                    headng,dist)
      arpsrasterx(i,j) = dist*SIN(deg2rad*headng)
      arpsrastery(i,j) = dist*COS(deg2rad*headng)
    END DO ! i loop
  END DO ! j loop

  RETURN
END SUBROUTINE getarpsrasterxy

!########################################################################
!########################################################################
!#########                                                      #########
!#########                SUBROUTINE getarpsdpaxy               #########
!#########                                                      #########
!#########                     Developed by                     #########
!#########     Center for Analysis and Prediction of Storms     #########
!#########                University of Oklahoma                #########
!#########                                                      #########
!########################################################################
!########################################################################

SUBROUTINE getarpsdpaxy(nx,ny,arpslat,arpslon,radarlat,radarlon,         &
                           arpsdpax,arpsdpay)

!------------------------------------------------------------------------
!
! PURPOSE:
!
! Calculates Cartesian coordinates of ARPS scalar points relative to
! a radar in HRAP (1/40 LFM) polar stereographic grid projection (used
! by NIDS DPA product).
!
!------------------------------------------------------------------------
!
! AUTHOR:
!
! Eric Kemp, 8 August 2001.
!    
!------------------------------------------------------------------------
!
! REFERENCE:
!
! Fulton, R. A., 1998:  WSR-88D Polar-to-HRAP Mapping.  Technical
!   Memorandum, Hydrologic Research Laboratory, 33 pp.  Available at:
!   http://hsp.nws.noaa.gov/oh/hrl/papers/wsr88d/hrapmap.pdf
!
! Keith Brewster, 4 April 2010
! Modification to accomodate conversion of n2aconst to nidscst.inc
!
!------------------------------------------------------------------------

  IMPLICIT NONE

!------------------------------------------------------------------------
!
! Declare arguments.
! 
!------------------------------------------------------------------------

  INTEGER, INTENT(IN) :: nx,ny              ! Dimensions of ARPS grid  
  REAL, INTENT(IN) :: arpslat(nx,ny)        ! Latitude of ARPS point
  REAL, INTENT(IN) :: arpslon(nx,ny)        ! Longitude of ARPS point
  REAL, INTENT(IN) :: radarlat,radarlon     ! Lat/Lon of radar.
  REAL, INTENT(OUT) :: arpsdpax(nx,ny)   ! x-coordinate of ARPS point.
  REAL, INTENT(OUT) :: arpsdpay(nx,ny)   ! y-coordinate of ARPS point.

!------------------------------------------------------------------------
!
! Internal variables
! 
!------------------------------------------------------------------------
  
  INTEGER :: i,j
!  INTEGER :: Is3, Js3
  REAL :: GIs, GJs
  INTEGER, PARAMETER :: Ip = 433
  INTEGER, PARAMETER :: Jp = 433
  REAL, PARAMETER :: Kc = 249.6348607
  REAL :: GI, GJ
  REAL :: sinls, cosls, coslambdaplus, sinlambdaplus
  REAL :: sindeltalambda,cosdeltalambda,deltalambda, sinl,cosl
  REAL :: Rl

!-----------------------------------------------------------------------
!
! Include file
!
!-----------------------------------------------------------------------

  INCLUDE 'nidscst.inc'

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  

!------------------------------------------------------------------------
!
! Calculate GIs, GJs.  See Appendix C, Sections 30.6.2 and 30.6.3 of
! Fulton (1998).
!
!------------------------------------------------------------------------

  sinls = SIN(radarlat*deg2rad)
  cosls = COS(radarlat*deg2rad)
  coslambdaplus = COS((radarlon + 105.)*deg2rad)
  sinlambdaplus = SIN((radarlon + 105.)*deg2rad)

  GIs = 1. + sinls
  GIs = Kc*cosls*sinlambdaplus/GIs
  GIs = GIs + REAL(Ip)

  GJs = 1 + sinls
  GJs = Kc*cosls*coslambdaplus/GJs
  GJs = GJs + REAL(Jp)

!------------------------------------------------------------------------
!
! Calculate Is3 and Js3.  See Appendix C, Section 30.6.3 of Fulton 
! (1998).  Note that GIs and GJs are multiplied by 10 to make Is3 and 
! Js3 valid on the 1/40 LFM grid.  This corrects an error in the
! NEXRAD documentation.
!
!------------------------------------------------------------------------

  DO j = 1,ny-1
    DO i = 1,nx-1

      deltalambda = arpslon(i,j) - radarlon
      sinl = SIN(arpslat(i,j)*deg2rad)
      cosl = COS(arpslat(i,j)*deg2rad)
    
      Rl = Kc*cosl/(1. + sinl)

      sindeltalambda = SIN(deltalambda*deg2rad)
      cosdeltalambda = COS(deltalambda*deg2rad)
             
      GI = Rl*sindeltalambda*coslambdaplus
      GI = GI + (Rl*cosdeltalambda*sinlambdaplus) + REAL(Ip)

      GJ = Rl*cosdeltalambda*coslambdaplus
      GJ = GJ + (Rl*sindeltalambda*sinlambdaplus) + REAL(Jp)

      arpsdpax(i,j) = DPAdx*(GI - GIs)*10.
      arpsdpay(i,j) = DPAdx*(GJ - GJs)*10.

    END DO ! i loop
  END DO ! j loop

  RETURN
END SUBROUTINE getarpsdpaxy
