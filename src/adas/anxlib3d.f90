!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINE FINDLC                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE findlc(nx,ny,xs,ys,xpt,ypt,ipt,jpt,ireturn)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Searches in x and y to find i,j location of xpt, ypt.
!
!  X and Y do not have to be on a regular grid, however it is
!  assumed that x and y are monotonically increasing as i and j
!  indices increase.
!
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster
!  April 1992.
!
!  MODIFICATION HISTORY:
!
!  February, 1993 (K. Brewster)
!  Additional documentation for ARPS 3.1 release
!
!  October, 1994 (K. Brewster)
!  Changed to reference scalar points.
!
!  July, 1995 (K. Brewster)
!  Changed to return error if extrapolation is required.
!
!  07/02/2001 (K. Brewster)
!  Cleaned up code to be more consistent with Fortran-90 conventions.
!
!  12/05/2005 (K. W. Thomas)
!  Prevent overlap in MPI mode.
!
!  01/10/2008 (Y. Wang)
!  Added handling of overlaps in MPI mode.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    xs       x coordinate of scalar points in physical/comp. space (m)
!    ys       y coordinate of scalar points in physical/comp. space (m)
!
!    xpt      location to find in x coordinate (m)
!    ypt      location to find in y coordinate (m)
!
!  OUTPUT:
!
!    ipt      i index to the west of desired location
!    jpt      j index to the south of desired location
!    ireturn  status indicator, 0/1/2/3 = good
!                              -3/-1    = extrapolation in x
!                              -2       = extrapolation in y
!
!  NOTE: 
!     1. ireturn based on the location of (xpt, ypt)
!                                            + -30
!      ny  -----------------------------     +----
!         | 21|          20          |32|    + 20
!    ny-1 |---|----------------------|--|    +----        
!         |   |                      |  |    +    
!         |   |                      |  |    +
!         |  1|           0          |12|    +  0
!         |   |                      |  |    +
!         |   |                      |  |    +
!       2 |--------------------------|--|    +----
!         | 11|          10          |22|    + 10
!       1 |___|______________________|__|    +____
!                                            +
!         1   2                    nx-1 nx   + -30
!                                            +
!    =========================================
!    -30  | 1 |           0          |12 | -30
!
!     2. Low and left boundary are included, but up and right boundary i
!        are exlusive. So only computations with index i+1, j+1 are 
!        expected. Computation for index i-1, j-1, however, is not supported.
!   
!-----------------------------------------------------------------------
!
!  Arguments
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx,ny          ! Dimensions of ARPS grids
  REAL :: xs(nx)            ! x coordinate of scalar grid points in
                            ! physical/comp. space (m)
  REAL :: ys(ny)            ! y coordinate of grid points in
                            ! physical/comp. space (m)

  REAL :: xpt               ! location to find in x coordinate
  REAL :: ypt               ! location to find in y coordinate
  INTEGER :: ipt            ! i index to the west of desired
                            ! location
  INTEGER :: jpt            ! j index to the south of desired
                            ! location
  INTEGER :: ireturn        ! status
!
!-----------------------------------------------------------------------
!
!  Misc. local variables
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
  ireturn=0

  IF(xpt >= xs(nx-1)) THEN
    ipt     = nx
    ireturn = -30
  ELSE IF (xpt < xs(1)) THEN
    ipt     = 1
    ireturn = -30
!  ELSE IF(xpt >= xs(nx-1)) THEN
!    ipt     = nx-1
!    ireturn = 12 
  ELSE IF (xpt < xs(2)) THEN
    ipt     = 1
    ireturn = 1
  ELSE
    DO i=3,nx-2
      IF(xpt < xs(i)) EXIT
    END DO
    ipt=i-1
  END IF

  IF(ypt >= ys(ny-1)) THEN
    jpt     = ny
    ireturn = -30 + ireturn
  ELSE IF (ypt < ys(1)) THEN
    jpt     =  1
    ireturn = -30 + ireturn
!  ELSE IF(ypt >= ys(ny-1)) THEN
!    jpt     = ny-1
!    ireturn =  20 + ireturn
  ELSE IF (ypt < ys(2)) THEN
    jpt = 1
    ireturn =  10 + ireturn
  ELSE
    DO j=3,ny-2
      IF(ypt < ys(j)) EXIT
    END DO
    jpt=j-1
  END IF

  RETURN
END SUBROUTINE findlc
!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINE COLINTA                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE colinta(nx,ny,nz,nvar,                                       &
           xs,ys,zp,xpt,ypt,ipt,jpt,anx,                                &
           su,sv,stheta,spres,shght,sqv,selev,                          &
           nlevs)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Interpolates ARPS history data in the horizontal to create
!  a column of data located at point xpt, ypt.
!
!  Bilinear interpolation is used.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster
!  April 1992.
!
!  MODIFICATION HISTORY:
!
!  October, 1992 (K. Brewster)
!  Conversion to ARPS 3.0.
!
!  October, 1994 (K. Brewster)
!  Conversion to ARPS 4.0.
!
!  Dec, 1998 (K. Brewster)
!  Changed RHstar to qv.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx,ny,nz Dimensions of ARPS grids.
!
!    xs       x coordinate of scalar points in physical/comp. space (m)
!    ys       y coordinate of scalar points in physical/comp. space (m)
!    zp       z coordinate of scalar grid points in physical space (m)
!
!    xpt      x coordinate of desired sounding (m)
!    ypt      y coordinate of desired sounding (m)
!
!    ipt      i index of grid point just west of xpt,ypt
!    jpt      j index of grid point just south of xpt,ypt
!
!    anx      Background field
!
!  OUTPUT:
!
!    su       Interpolated u wind component.  (m/s)
!    sv       Interpolated v wind component.  (m/s)
!    stheta   Interpolated potential temperature (K).
!    spres    Interpolated pressure. (Pascals)
!    shght    Interpolated height (meters)
!    sqv      Interpolated specific humidity (kg/kg)
!    selev    Interpolated surface elevation (m)
!    nlevs    Number of above-ground sounding levels.
!
!-----------------------------------------------------------------------
!
!  Variable Declarations:
!
!-----------------------------------------------------------------------
!
!  Arguments -- location data
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx,ny,nz,nvar     ! Dimensions of ARPS grids.
  REAL :: xs(nx)               ! x coordinate of grid points in
                               ! physical/comp. space (m)
  REAL :: ys(ny)               ! y coordinate of grid points in
                               ! physical/comp. space (m)
  REAL :: zp(nx,ny,nz)         ! z coordinate of grid points in
                               ! physical space (m)
  REAL :: xpt                  ! location to find in x coordinate (m)
  REAL :: ypt                  ! location to find in y coordinate (m)
  INTEGER :: ipt               ! i index to the west of desired
                               ! location
  INTEGER :: jpt               ! j index to the south of desired
                               ! location
!
!-----------------------------------------------------------------------
!
!  Arguments -- background field
!
!-----------------------------------------------------------------------
!
  REAL :: anx(nx,ny,nz,nvar)
!
!-----------------------------------------------------------------------
!
!  Arguments -- Extracted sounding variables
!
!-----------------------------------------------------------------------
!
  REAL :: su(nz),sv(nz),stheta(nz),sqv(nz)
  REAL :: spres(nz),shght(nz)
  REAL :: selev
  INTEGER :: nlevs
!
!-----------------------------------------------------------------------
!
!  Functions called
!
!-----------------------------------------------------------------------
!
  REAL :: aint2d
!
!-----------------------------------------------------------------------
!
!  Include files
!
!-----------------------------------------------------------------------
!
  INCLUDE 'phycst.inc'
!
!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: k,in,jn
  REAL :: delx,ddx,dely,ddy,w1,w2,w3,w4
  REAL :: t2,t3,hmid,tmid
!
!-----------------------------------------------------------------------
!
!  Find corner weights
!
!-----------------------------------------------------------------------
!
  in=ipt+1
  delx=xs(in)-xs(ipt)
  IF(ABS(delx) > 0.) THEN
    ddx=(xpt-xs(ipt))/delx
  ELSE
    ddx=0.
  END IF

  jn=jpt+1
  dely=ys(jn)-ys(jpt)
  IF(ABS(dely) > 0.) THEN
    ddy=(ypt-ys(jpt))/dely
  ELSE
    ddy=0.
  END IF

  w1=(1.-ddx)*(1.-ddy)
  w2=ddx*(1.-ddy)
  w3=ddx*ddy
  w4=(1.-ddx)*ddy
!
!-----------------------------------------------------------------------
!
!  Interpolate all variables at all levels.
!
!-----------------------------------------------------------------------
!
  nlevs=nz-1
  DO k=1,nz-1
    shght(k)=                                                           &
        aint2d(nx,ny,nz,    zp,ipt,jpt,k,in,jn,w1,w2,w3,w4)
    su(k)=                                                              &
        aint2d(nx,ny,nz, anx(1,1,1,1),ipt,jpt,k,in,jn,w1,w2,w3,w4)
    sv(k)=                                                              &
        aint2d(nx,ny,nz, anx(1,1,1,2),ipt,jpt,k,in,jn,w1,w2,w3,w4)
    spres(k)=                                                           &
        aint2d(nx,ny,nz, anx(1,1,1,3),ipt,jpt,k,in,jn,w1,w2,w3,w4)
    stheta(k)=                                                          &
        aint2d(nx,ny,nz, anx(1,1,1,4),ipt,jpt,k,in,jn,w1,w2,w3,w4)
    sqv(k)=                                                             &
        aint2d(nx,ny,nz, anx(1,1,1,5),ipt,jpt,k,in,jn,w1,w2,w3,w4)
  END DO
!
!-----------------------------------------------------------------------
!
!  Get height at scalar points, since zp was defined at w points.
!
!-----------------------------------------------------------------------
!
  selev=shght(2)
  DO k=1,nz-1
    shght(k)=0.5*(shght(k+1)+shght(k))
  END DO
!
!-----------------------------------------------------------------------
!
!  Get a value at the surface, by linearly interpolating
!  between the 1st and second levels.
!
!-----------------------------------------------------------------------
!
  w2=(selev-shght(1))/(shght(2)-shght(1))
  w1=1.-w2
  su(1)=w1*    su(1) + w2*    su(2)
  sv(1)=w1*    sv(1) + w2*    sv(2)
  stheta(1)=w1*stheta(1) + w2*stheta(2)
  sqv(1)=w1*   sqv(1) + w2*   sqv(2)
  shght(1)=selev
!
!-----------------------------------------------------------------------
!
!  Integrate downward to get the pressure at level 1.
!
!-----------------------------------------------------------------------
!
  t3=stheta(3)*(spres(3)/100000.)**rddcp
  t2=stheta(2)*(spres(2)/100000.)**rddcp
  hmid=0.5*(shght(2)+shght(1))
  tmid=t3+((shght(3)-hmid)/(shght(3)-shght(2)))*(t2-t3)
  spres(1)=spres(2)*EXP(g*(shght(2)-shght(1))/(rd*tmid))
  RETURN
END SUBROUTINE colinta
!
!##################################################################
!##################################################################
!######                                                      ######
!######                 FUNCTION AINT2D                      ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

  FUNCTION aint2d(nx,ny,nz, a,im,jm,k,in,jn,w1,w2,w3,w4)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Applies bilinear interpolation to array "a" to find value of "a"
!  at a point between x(im) and x(in) and y(im) and y(in).
!
!  Weights are determined outside this function so that the weights
!  can be determined once, then this function used for all variables.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster
!  April 1992.
!
!  MODIFICATION HISTORY:
!
!  October, 1992 (K. Brewster)
!  Conversion to ARPS 3.0.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!  nx,ny,nz    Dimensions of data arrays.
!  a           Array to be interpolated.
!  im          First array index -- left of desired point.
!  in          First array index -- right of desired point.
!  jm          Second array index -- below desired point.
!  jn          Second array index -- above deaired point.
!  k           Third array index.
!  w1,w2,w3,w4    Weights of surrounding pts.
!
!  Location of surrounding points
!
!         im,jn       in,jn
!
!             xpt,ypt
!               +
!
!         im,jm       in,jm
!
!  Weights for surrounding points
!
!           w4          w3
!
!             xpt,ypt
!               +
!
!           w1          w2
!
!-----------------------------------------------------------------------
!
!  Variable declarations
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  REAL :: aint2d
  INTEGER :: nx,ny,nz          ! Dimensions of data arrays

  REAL :: a(nx,ny,nz)          ! Array to be interpolated
  INTEGER :: im                ! First array index -- left of desired
                               ! point
  INTEGER :: jm                ! Second array index -- below desired
                               ! point
  INTEGER :: k                 ! Third array index
  INTEGER :: in                ! First array index -- right of desired
                               ! point
  INTEGER :: jn                ! Second array index -- above deaired
                               ! point
  REAL :: w1,w2,w3,w4          ! Weights of surrounding pts
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  aint2d=w1*a(im,jm,k) +                                                &
         w2*a(in,jm,k) +                                                &
         w3*a(in,jn,k) +                                                &
         w4*a(im,jn,k)

  RETURN
  END FUNCTION aint2d
!
!##################################################################
!##################################################################
!######                                                      ######
!######                 FUNCTION ALTTOSTPR                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

  FUNCTION alttostpr(altim,elev)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Use altimeter equation (standard atmosphere) to
!  change altimeter setting to pressure at the
!  station elevation (elev).
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster
!  July, 1992
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    altim     altimeter setting in arbitrary pressure units
!    elev      station elevation (m MSL)
!
!  OUTPUT:
!
!    alttostpr   station pressure in the same units as the
!                input altimeter setting
!
!-----------------------------------------------------------------------
!
!  Variable Declarations:
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
!  Arguments
!
!-----------------------------------------------------------------------
!
  REAL :: alttostpr   ! station pressure
  REAL :: altim       ! altimeter setting
  REAL :: elev        ! elevation in meters above sea level
!
!-----------------------------------------------------------------------
!
!  Physical constants
!
!-----------------------------------------------------------------------
!
  REAL :: TO,gamusd,rdgas,g,c1
  PARAMETER ( TO = 288., & ! degrees K sea-level temp in US std atmos
      gamusd = 0.0065, & ! K /m std atmos lapse rate
       rdgas = 287.04, & ! J/(K*kg) gas constant for dry air
           g = 9.80616,    & ! m/(s*s)
          c1 = (g/(gamusd*rdgas)) )

!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  alttostpr=altim*(1.-(elev*gamusd/TO)) ** c1

  RETURN
  END FUNCTION alttostpr
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                 FUNCTION MSLTOSTPR                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

  FUNCTION msltostpr(msl,tk,elev)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Use altimeter equation (standard atmosphere) to
!  change altimeter setting to pressure at the
!  station elevation (elev).
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster
!  July, 1992
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    msl       mean-sea-level pressure in arbitrary pressure units
!    tk        surface temperature (Kelvin)
!    elev      station elevation (m MSL)
!
!  OUTPUT:
!
!    msltostpr   station pressure in the same units as the
!                mean sea-level pressure
!
!-----------------------------------------------------------------------
!
!  Variable Declarations:
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
!  Arguments
!
!-----------------------------------------------------------------------
!
  REAL :: msltostpr   ! station pressure
  REAL :: msl         ! mean sea-level pressure
  REAL :: tk          ! surface temperature (K)
  REAL :: elev        ! elevation in meters above sea level
!
!-----------------------------------------------------------------------
!
!  Physical constants
!
!-----------------------------------------------------------------------
!
  REAL :: gamusd,rdgas,g,c1
  PARAMETER ( gamusd = 0.0065, & ! K /m std atmos lapse rate
       rdgas = 287.04, & ! J/(K*kg) gas constant for dry air
           g = 9.80616,    & ! m/(s*s)
          c1 = (g/(gamusd*rdgas)) )
!
  msltostpr=msl*((tk/(tk+gamusd*elev)) ** c1)
!
  RETURN
  END FUNCTION msltostpr
!
!##################################################################
!##################################################################
!######                                                      ######
!######                 FUNCTION STPRTOHPR                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

  FUNCTION stprtohpr(stpr,elev,hgt)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Use altimeter equation (standard atmosphere) to
!  change station pressure (stpr) to pressure at another
!  height (hgt).   Heights in meters.  stpr  and output
!  are in the same units.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster
!  July, 1992
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    stpr     Station pressure (arbitrary pressure units)
!    elev     Station elevation (m MSL)
!    hgt      Standard atmosphere height requested (m MSL)
!
!  OUTPUT
!
!    stprtohpr  pressure (same units as input stpr)
!
!-----------------------------------------------------------------------
!
!  Variable Declarations:
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
!  Arguments
!
!-----------------------------------------------------------------------
!
  REAL :: stprtohpr,stpr,elev,hgt
!
!-----------------------------------------------------------------------
!
!  Physical constants
!
!-----------------------------------------------------------------------
!
  REAL :: TO,gamusd,rdgas,g,c1
  PARAMETER  ( TO = 288., & ! degrees K sea-level temp in US std atmos
           gamusd = 0.0065, & ! K /m std atmos lapse rate
           rdgas = 287.04, & ! J/(K*kg) gas constant for dry air
               g = 9.80616,    & ! m/(s*s)
              c1 = (g/(gamusd*rdgas)) )
  REAL :: telev
!
  telev=TO - (gamusd*elev)
  stprtohpr=stpr*(1.-((hgt-elev)*gamusd/telev)) ** c1
  RETURN
  END FUNCTION stprtohpr
!
!##################################################################
!##################################################################
!######                                                      ######
!######                 FUNCTION STPRTOPMSL                  ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

  FUNCTION stprtopmsl(stpr,tk,elev)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Use altimeter equation (standard atmosphere) to
!  change station pressure a given station elevation (elev).
!  to mean sea-level pressure.
!
!-----------------------------------------------------------------------
!
!  INPUT
!    stpr    station pressure (arbitrary pressure units)
!    elev    elevation (meters MSL)
!
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
!  Arguments
!
!-----------------------------------------------------------------------
!
  REAL :: stprtopmsl,stpr,tk,elev
!
!-----------------------------------------------------------------------
!
!  Physical constants
!
!-----------------------------------------------------------------------
!
  REAL :: TO,gamusd,rdgas,g,c1
  PARAMETER ( TO = 288., & ! degrees K sea-level temp in US std atmos
      gamusd = 0.0065, & ! K /m std atmos lapse rate
       rdgas = 287.04, & ! J/(K*kg) gas constant for dry air
           g = 9.80616,    & ! m/(s*s)
          c1 = (g/(gamusd*rdgas)) )
!
  stprtopmsl=stpr*(1.+(elev*gamusd/tk)) ** c1
  RETURN
  END FUNCTION stprtopmsl
!
!##################################################################
!##################################################################
!######                                                      ######
!######                   FUNCTION ZTOPSA                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

  FUNCTION ztopsa(z)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  This routine converts a height in meters into a pressure
!  in a standard atmosphere in millibars.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Tracy Smith, Stan Benjamin     NOAA/FSL
!
!  MODIFICATION HISTORY:
!
!    1999/11/19 (D. Hou)
!      Incorporated into ADAS for use by MDCRS data (in prepsng.f).
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    z         altimeter height (m MSL)
!
!  OUTPUT:
!
!    ztopsa      station pressure in millibars
!
!-----------------------------------------------------------------------
!
!  Variable Declarations:
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
!  Arguments
!
!-----------------------------------------------------------------------
!
  REAL :: ztopsa      ! station pressure
  REAL :: z           ! altimeter height (m MSL)
!
!-----------------------------------------------------------------------
!
!  Physical constants
!
!-----------------------------------------------------------------------
!
  REAL :: t0,gamma,p0,p11,z11,c1,c2,flag,flg

  DATA flag,flg/99999.,99998./
  DATA t0,gamma,p0/288.,.0065,1013.2/
  DATA c1,c2/5.256,14600./
  DATA z11,p11/11000.,226.0971/
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  IF (z > flg) THEN
    ztopsa = flag
  ELSE IF (z < z11) THEN
    ztopsa = p0*((t0-gamma*z)/t0)**c1
  ELSE
    ztopsa = p11*10.**((z11-z)/c2)
  END IF

  RETURN
  END FUNCTION ztopsa
