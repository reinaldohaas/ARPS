!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE V3DDUMP                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Coastal Meteorology Research Project  (CMRP)     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE v5ddump(nx,ny,nz,nzsoil,nstyps, filnam, istgr,               &
           u,v,w,ptprt,pprt,qv,qscalar,tke,kmh,kmv,                     &
           ubar,vbar,wbar,ptbar,pbar,rhobar,qvbar,                      &
           x,y,z,zp,zpsoil,                                             &
           soiltyp,stypfrct,vegtyp,lai,roufns,veg,                      &
           tsoil,qsoil,wetcanp,snowdpth,                                &
           raing,rainc,prcrate,                                         &
           radfrc,radsw,rnflx,radswnet,radlwin,                         &
           usflx,vsflx,ptsflx,qvsflx,                                   &
           tem1,tem2)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Write history data into the file "filnam" in the Vis5D format.
!
!  All output data are located at the grid cell centers.
!
!       X_center(i) = ( X_edge(i) + X_edge(i+1) )/2,  i = 1, n-1, 1
!
!  The last edge value were kept unchanged so that it can be used to
!  restore the stagger grid point values.
!
!       X_edge(n) = X_center(n)
!       X_edge(i) = 2*X_center(i) - X_edge(i+1),      i = n-1, 1, -1
!
!  Since total and perturbation variables will be dumped, the base
!  state variables can be obtained by Xbar = X - Xprt.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Fanyou kong
!  4/28/97.
!
!  1/31/2002 (Ming Xue)
!  Fixed a problem that occurs when nstyp>1.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!    nzsoil   Number of grid points in the soil
!
!    nstyps   number of soil types
!
!    u        x component of velocity at a given time level (m/s)
!    v        y component of velocity at a given time level (m/s)
!    w        Vertical component of Cartesian velocity at a given
!             time level (m/s)
!    ptprt    Perturbation potential temperature at a given time
!             level (K)
!    pprt     Perturbation pressure at a given time level (Pascal)
!    qv       Water vapor specific humidity at a given time level
!             (kg/kg)
!    qc       Cloud water mixing ratio at a given time level (kg/kg)
!    qr       Rainwater mixing ratio at a given time level (kg/kg)
!    qi       Cloud ice mixing ratio at a given time level (kg/kg)
!    qs       Snow mixing ratio at a given time level (kg/kg)
!    qh       Hail mixing ratio at a given time level (kg/kg)
!    tke      Turbulent Kinetic Energy ((m/s)**2)
!
!    kmh      Horizontal turb. mixing coef. for momentum ( m**2/s )
!    kmv      Vertical turb. mixing coef. for momentum ( m**2/s )
!
!    ubar     Base state zonal velocity component (m/s)
!    vbar     Base state meridional velocity component (m/s)
!    ptbar    Base state potential temperature (K)
!    pbar     Base state pressure (Pascal)
!    rhobar   Base state density (kg/m**3)
!    qvbar    Base state water vapor specific humidity (kg/kg)
!
!    x        x coordinate of grid points in physical/comp. space (m)
!    y        y coordinate of grid points in physical/comp. space (m)
!    z        z coordinate of grid points in computational space (m)
!    zp       Vertical coordinate of grid points in physical space(m)
!    zpsoil   Vertical coordinate of grid points in the soil (m)
!
!    soiltyp  Soil type
!    vegtyp   Vegetation type
!    lai      Leaf Area Index
!    roufns   Surface roughness
!    veg      Vegetation fraction
!
!    tsoil    Soil temperature (K)
!    qsoil    Soil moisture (m**3/m**3)
!    wetcanp  Canopy water amount
!
!    raing    Grid supersaturation rain
!    rainc    Cumulus convective rain
!
!  OUTPUT:
!
!    None.
!
!  WORK ARRAY:
!
!    tem1     Temporary work array. Used to store data for Vis5D
!             write function, which has different order of nx,ny
!    tem2     Temporary work array.
!
!
!-----------------------------------------------------------------------
!
!  The following parameters are passed into this subroutine through
!  a common block in globcst.inc, and they determine which
!  variables are output.
!
!  varout =0 or 1. If varout=0, model perturbation variables are not
!                  dumped.
!  mstout =0 or 1. If mstout=0, water variables are not dumped.
!  rainout=0 or 1. If rainout=0, rain variables are not dumped.
!  prcout =0 or 1. If prcout=0, precipitation rates are not dumped.
!  iceout =0 or 1. If iceout=0, qi, qs and qh are not dumped.
!  tkeout =0 or 1. If tkeout=0, tke is not dumped.
!  trbout =0 or 1. If trbout=0, kmh and kmv are not dumped
!  sfcout =0 or 1. If sfcout=0, surface variables are not dumped.
!  landout=0 or 1. If landout=0, surface property arrays are not dumped.
!  radout =0 or 1. If radout=0, radiation arrays are not dumped.
!  flxout =0 or 1. If flxout=0, precipitation rates are not dumped.
!
!-----------------------------------------------------------------------

  IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'grid.inc'
  INCLUDE 'v5df.inc'

!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  INTEGER :: nx,ny,nz          ! Number of grid points in 3 directions
  INTEGER :: nzsoil            ! Number of grid points in the soil

  INTEGER :: istgr             ! Flag for dumping stager point data

  REAL :: u     (nx,ny,nz)     ! Total u-velocity (m/s)
  REAL :: v     (nx,ny,nz)     ! Total v-velocity (m/s)
  REAL :: w     (nx,ny,nz)     ! Total w-velocity (m/s)
  REAL :: ptprt (nx,ny,nz)     ! Perturbation potential temperature (K)
  REAL :: pprt  (nx,ny,nz)     ! Perturbation pressure (Pascal)

  REAL :: qv    (nx,ny,nz)     ! Water vapor specific humidity (kg/kg)

  REAL :: qscalar (nx,ny,nz,nscalar)

  REAL :: tke   (nx,ny,nz)     ! Turbulent Kinetic Energy ((m/s)**2)
  REAL :: kmh   (nx,ny,nz)     ! Horizontal turb. mixing coef. for
                               ! momentum. ( m**2/s )
  REAL :: kmv   (nx,ny,nz)     ! Vertical turb. mixing coef. for
                               ! momentum. ( m**2/s )

  REAL :: ubar  (nx,ny,nz)     ! Base state u-velocity (m/s)
  REAL :: vbar  (nx,ny,nz)     ! Base state v-velocity (m/s)
  REAL :: wbar  (nx,ny,nz)     ! Base state w-velocity (m/s)
  REAL :: ptbar (nx,ny,nz)     ! Base state potential temperature (K)
  REAL :: pbar  (nx,ny,nz)     ! Base state pressure (Pascal)
  REAL :: rhobar(nx,ny,nz)     ! Base state air density (kg/m**3)
  REAL :: qvbar (nx,ny,nz)     ! Base state water vapor specific
                               ! humidity(kg/kg)

  REAL :: x     (nx)           ! x-coord. of the physical and
                               ! computational grid. Defined at u-point.
  REAL :: y     (ny)           ! y-coord. of the physical and
                               ! computational grid. Defined at v-point.
  REAL :: z     (nz)           ! z-coord. of the computational grid.
                               ! Defined at w-point on the staggered grid.
  REAL :: zp    (nx,ny,nz)     ! Physical height coordinate defined at
                               ! w-point of the staggered grid.
  REAL :: zpsoil(nx,ny,nzsoil) ! Physical height coordinate defined at
                               ! w-point of the soil.

  INTEGER :: nstyps
  INTEGER :: soiltyp(nx,ny,nstyps) ! Soil type
  REAL    :: stypfrct(nx,ny,nstyps)! Fraction of soil types
  INTEGER :: vegtyp (nx,ny)        ! Vegetation type
  REAL :: lai    (nx,ny)           ! Leaf Area Index
  REAL :: roufns (nx,ny)           ! Surface roughness
  REAL :: veg    (nx,ny)           ! Vegetation fraction

  REAL :: tsoil (nx,ny,nzsoil,0:nstyps) ! Soil temperature (K)
  REAL :: qsoil (nx,ny,nzsoil,0:nstyps) ! Soil moisture
  REAL :: wetcanp(nx,ny,0:nstyps)       ! Canopy water amount
  REAL :: snowdpth(nx,ny)               ! Snow depth (m)

  REAL :: raing(nx,ny)                ! Grid supersaturation rain
  REAL :: rainc(nx,ny)                ! Cumulus convective rain
  REAL :: prcrate(nx,ny,4)     ! precipitation rates (kg/(m**2*s))
                               ! prcrate(1,1,1) = total precip. rate
                               ! prcrate(1,1,2) = grid scale precip. rate
                               ! prcrate(1,1,3) = cumulus precip. rate
                               ! prcrate(1,1,4) = microphysics precip. rate

  REAL :: radfrc(nx,ny,nz)     ! Radiation forcing (K/s)
  REAL :: radsw (nx,ny)        ! Solar radiation reaching the surface
  REAL :: rnflx (nx,ny)        ! Net radiation flux absorbed by surface
  REAL :: radswnet(nx,ny)      ! Net shortwave radiation
  REAL :: radlwin(nx,ny)       ! Incoming longwave radiation

  REAL :: usflx (nx,ny)        ! Surface flux of u-momentum (kg/(m*s**2))
  REAL :: vsflx (nx,ny)        ! Surface flux of v-momentum (kg/(m*s**2))
  REAL :: ptsflx(nx,ny)        ! Surface heat flux (K*kg/(m*s**2))
  REAL :: qvsflx(nx,ny)        ! Surface moisture flux (kg/(m**2*s))

  REAL :: tem1  (ny,nx,nz)            ! Temporary work array (note the
                                      ! order of indexes)
  REAL :: tem2  (nx,ny,nz)            ! Temporary work array

!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k
  INTEGER :: nq
  INTEGER :: m
  INTEGER :: ierr,is, it,iv, nl(maxvars),nr,nc, maxnl
  INTEGER :: varnum,numtimes
  INTEGER :: dates(maxtimes),times(maxtimes)

  INTEGER :: compressmode

  CHARACTER (LEN=*)  :: filnam
  CHARACTER (LEN=10) :: varnam(maxvars)
  CHARACTER (LEN=20) :: chrstr

  REAL :: proj_args(100)
  REAL :: vert_args(maxlevels)

  INTEGER :: projection
  INTEGER :: vertical

  REAL :: latnot(2), x0,y0
  REAL :: xbgn,ybgn,zbgn,xinc,yinc,zinc,xmin,xmax,ymin,ymax
  INTEGER :: yy,date0,time0,idatime,hh,mm,ss
!  INTEGER :: yy,mmon,dd,hh,mm,ss
  INTEGER :: abstsec,jdy

  CHARACTER(LEN=40) :: qnames_upper(20)
  CHARACTER(LEN=40) :: upcase
  REAL :: factor

  INTEGER :: mndys(12)                 ! days for each months
  DATA mndys/0,31,59,90,120,151,181,212,243,273,304,334/


!
!-----------------------------------------------------------------------

  DATA numtimes/1/
  DATA varnam(1) /'ZP'/
  DATA compressmode/1/

  DO nq = 1,nscalar
    qnames_upper(nq) = upcase(qnames(nq))
  END DO
!
!-----------------------------------------------------------------------
!
!  Calculate time and grid parameters
!
!-----------------------------------------------------------------------
!
  varnum = 1

  yy = year - 100*(year/100)
  date0 = 1000*yy + mndys(month) +day
  time0 = 10000*hour + 100*minute + second
  dates(1) = date0

  idatime = INT(curtim)
  hh = idatime/3600
  idatime = idatime - 3600*hh
  mm = idatime/60
  ss = idatime - 60*mm

  hh = hh + hour
  IF(hh > 24) THEN
    hh = hh - 24
    dates(1) = dates(1) + 1
  END IF

  times(1) = 10000*hh+100*mm+ss

!  yy=year
!  mmon=month
!  dd=day
!  hh=hour
!  mm=minute
!  ss=second
!  CALL ctim2abss( yy,mmon,dd,hh,mm,ss, abstsec )
!  abstsec=abstsec+curtim
!  CALL abss2ctim( abstsec,yy,mmon,dd,hh,mm,ss)

!  yy = yy - 100*(yy/100)
!  CALL julday(yy, mmon, dd, jdy)
!  dates(1)=1000*yy+jdy
!  times(1)=10000*hh + 100*mm + ss

  nr=ny
  nc=nx

  DO i=1, varnum
    nl(i) = nz
  END DO

  xbgn =  x(2)
  ybgn =  y(2)
  zbgn =  z(2)

  xinc = (x(2) - x(1))
  yinc = (y(2) - y(1))

  xmin = xbgn
  xmax = x(nx) + xinc
  ymin = ybgn
  ymax = y(ny) + yinc

!-----------------------------------------------------------------------
!
!   Determine projection parameters used in the program
!
!-----------------------------------------------------------------------

  IF(mapproj /= 0) THEN
    latnot(1) = trulat1
    latnot(2) = trulat2


    CALL getmapr(mapproj, sclfct, latnot,trulon,x0,y0)
    WRITE(6,*) 'mapproj, sclfct,trulat1,trulat2,trulon:'
    WRITE(6,*) mapproj, sclfct,trulat1,trulat2,trulon
    WRITE(6,*) 'ctrlat,ctrlon:',ctrlat,ctrlon
    WRITE(6,*) 'x0,y0:', x0,y0

    IF( trulon > 180.0) trulon = trulon - 360.0
    IF(-trulon > 180.0) trulon = trulon + 360.0
    trulon = -trulon

  END IF

  IF(mapproj > 2) THEN
    WRITE(6,*) 'mapproj=',mapproj,                                      &
              ' is not a supported map projection by Vis5D, STOP'
    CALL EXIT(1)
  END IF

  IF(mapproj == 0) THEN    ! Linear, rectangular, generic units

    projection = 0

    proj_args(1) = ymax
    proj_args(2) = xmin
    proj_args(3) = yinc
!      proj_args(4) = -xinc
    proj_args(4) = xinc

  ELSE IF (mapproj == 1) THEN ! Stereographic projection
                              ! In Vis5D only azimuthal stereographic
                              ! is supported, so the use of this
                              ! projection in ARPS data is only an
                              ! approximation which needs trulat1
                              ! not less than 75 degree

    IF(ABS(trulat1) < 75.0) THEN
      WRITE(6,*) 'trulat1=',trulat1,' is less than 75deg',              &
                 ' which may cause large error'
      WRITE(6,*) 'in polar stereographic projection mode, STOP'
      CALL EXIT(1)
    END IF

    projection = 3

    proj_args(1) = 90.0
    proj_args(2) = trulon
    proj_args(3) = y0/yinc + FLOAT(ny-2)
    proj_args(4) = -x0/xinc + 2.0
    proj_args(5) = xinc/1000.0

  ELSE IF (mapproj == 2) THEN ! Lambert conformal projection

    projection = 2

    IF(trulat1 >= trulat2) THEN
      proj_args(1) = trulat1
      proj_args(2) = trulat2
    ELSE
      proj_args(1) = trulat2
      proj_args(2) = trulat1
    END IF
    proj_args(3) = y0/yinc + FLOAT(ny-2)
    proj_args(4) = -x0/xinc + 1.0
    proj_args(5) = trulon
    proj_args(6) = xinc/1000.0

  END IF

  WRITE(6,*) (proj_args(i),i=1,6)

!-----------------------------------------------------------------------
!
!   Determine vertical coordinate parameters used in the program
!
!-----------------------------------------------------------------------

  vertical = 2                ! unequally spaced (km)

  zinc = (zp(1,1,2)-zp(1,1,1))/2.0
  DO k = 1,nz
    vert_args(k) = (zp(1,1,k)+zinc)/1000.0
  END DO

!     write(6,*) 'vertical:',vertical
!     write(6,*) 'vert_args(k):'
!     write(6,'(i10,f12.4/)') (k,vert_args(k),k=1,nz)

!-----------------------------------------------------------------------
!
!   Determine the output fields according to input switches
!
!-----------------------------------------------------------------------

  IF ( varout == 1 ) THEN

    varnum = varnum + 1
    varnam(varnum) = 'U'
    nl(varnum) = nz

    varnum = varnum + 1
    varnam(varnum) = 'V'
    nl(varnum) = nz

    varnum = varnum + 1
    varnam(varnum) = 'W'
    nl(varnum) = nz

    varnum = varnum + 1
    varnam(varnum) = 'UPRT'
    nl(varnum) = nz

    varnum = varnum + 1
    varnam(varnum) = 'VPRT'
    nl(varnum) = nz

    varnum = varnum + 1
    varnam(varnum) = 'PT'
    nl(varnum) = nz

    varnum = varnum + 1
    varnam(varnum) = 'PTPRT'
    nl(varnum) = nz

    varnum = varnum + 1
    varnam(varnum) = 'P'
    nl(varnum) = nz

    varnum = varnum + 1
    varnam(varnum) = 'PPRT'
    nl(varnum) = nz

    varnum = varnum + 1
    varnam(varnum) = 'VORT'
    nl(varnum) = nz

    varnum = varnum + 1
    varnam(varnum) = 'DIV'
    nl(varnum) = nz

  END IF

  IF ( mstout == 1 ) THEN

    varnum = varnum + 1
    varnam(varnum) = 'QV'
    nl(varnum) = nz

    varnum = varnum + 1
    varnam(varnum) = 'QVPRT'
    nl(varnum) = nz

    print*,'Inside v5ddump: nscalar = ',nscalar

    DO nq = 1,nscalar
      varnum = varnum + 1
      varnam(varnum) = TRIM(qnames_upper(nq))
      nl(varnum) = nz
    END DO

    IF ( rainout == 1 ) THEN

      varnum = varnum + 1
      varnam(varnum) = 'RAING'
      nl(varnum) = 1

      varnum = varnum + 1
      varnam(varnum) = 'RAINC'
      nl(varnum) = 1

    END IF

    IF ( prcout == 1 ) THEN

      varnum = varnum + 1
      varnam(varnum) = 'PRCRATE1'
      nl(varnum) = 1

      varnum = varnum + 1
      varnam(varnum) = 'PRCRATE2'
      nl(varnum) = 1

      varnum = varnum + 1
      varnam(varnum) = 'PRCRATE3'
      nl(varnum) = 1

      varnum = varnum + 1
      varnam(varnum) = 'PRCRATE4'
      nl(varnum) = 1

    END IF
  END IF

  IF ( tkeout == 1 ) THEN

    varnum = varnum + 1
    varnam(varnum) = 'TKE'
    nl(varnum) = nz

  END IF

  IF ( trbout == 1 ) THEN

    varnum = varnum + 1
    varnam(varnum) = 'KMH'
    nl(varnum) = nz

    varnum = varnum + 1
    varnam(varnum) = 'KMV'
    nl(varnum) = nz

  END IF

  IF ( sfcout == 1 ) THEN

    IF ( nstyp <= 1 ) THEN

      varnum = varnum + 1
      varnam(varnum) = 'TSOIL'
      nl(varnum) = nzsoil

      varnum = varnum + 1
      varnam(varnum) = 'QSOIL'
      nl(varnum) = nzsoil

      varnum = varnum + 1
      varnam(varnum) = 'WR'
      nl(varnum) = 1

    ELSE

      DO is=0,nstyp
        WRITE (chrstr,'(i1)') is

        varnum = varnum + 1
        varnam(varnum) = 'TSOIL'//chrstr(1:1)
        nl(varnum) = nzsoil

        varnum = varnum + 1
        varnam(varnum) = 'QSOIL'//chrstr(1:1)
        nl(varnum) = nzsoil

        varnum = varnum + 1
        varnam(varnum) = 'WR'//chrstr(1:1)
        nl(varnum) = 1

      END DO

    END IF

    varnum = varnum + 1
    varnam(varnum) = 'SC'
    nl(varnum) = 1

  END IF

  IF(landout == 1) THEN

    IF (nstyp <= 1) THEN

      varnum = varnum + 1
      varnam(varnum) = 'SOILTY'
      nl(varnum) = 1

    ELSE
      DO is=1,nstyp
        WRITE (chrstr,'(i1)') is

        varnum = varnum + 1
        varnam(varnum) = 'SOILTY'//chrstr(1:1)
        nl(varnum) = 1

        varnum = varnum + 1
        varnam(varnum) = 'SOILFR'//chrstr(1:1)
        nl(varnum) = 1

      END DO
    END IF

    varnum = varnum + 1
    varnam(varnum) = 'VEGTYP'
    nl(varnum) = 1

    varnum = varnum + 1
    varnam(varnum) = 'LAI'
    nl(varnum) = 1

    varnum = varnum + 1
    varnam(varnum) = 'ROUFNS'
    nl(varnum) = 1

    varnum = varnum + 1
    varnam(varnum) = 'VEG'
    nl(varnum) = 1

  END IF

  IF ( radout == 1 ) THEN

    varnum = varnum + 1
    varnam(varnum) = 'RADFRC'
    nl(varnum) = nz

    varnum = varnum + 1
    varnam(varnum) = 'RADSW'
    nl(varnum) = 1

    varnum = varnum + 1
    varnam(varnum) = 'RNFLX'
    nl(varnum) = 1

  END IF

  IF ( flxout == 1 ) THEN

    varnum = varnum + 1
    varnam(varnum) = 'USFLX'
    nl(varnum) = 1

    varnum = varnum + 1
    varnam(varnum) = 'VSFLX'
    nl(varnum) = 1

    varnum = varnum + 1
    varnam(varnum) = 'PTSFLX'
    nl(varnum) = 1

    varnum = varnum + 1
    varnam(varnum) = 'QVSFLX'
    nl(varnum) = 1

  END IF

  WRITE(6,*) 'Number of variables to be written: ',varnum

!----------------------------------------------------------------------
!
!  Open Vis5D file & write header to it
!
!-----------------------------------------------------------------------

  m = v5dcreate (filnam, numtimes, varnum, nr, nc, nl,                  &
                 varnam, times, dates, compressmode,                    &
                 projection, proj_args, vertical, vert_args)


  IF (m == 0) CALL EXIT(1)

  maxnl = nl(1)
  DO i = 1, varnum
    IF(nl(i) > maxnl) THEN
      maxnl = nl(i)
    END IF
  END DO

!----------------------------------------------------------------------
!
!  Output data to Vis5D file
!
!-----------------------------------------------------------------------

  DO it = 1, numtimes
    DO iv = 1, varnum

      DO k=1,nz
        DO i=1,nx
          DO j=1,ny
            tem1(j,i,k) = 0.0
            tem2(i,j,k) = 0.0
          END DO
        END DO
      END DO

      IF(varnam(iv) == 'ZP') THEN
        CALL v5dsetunits(iv, "m\0")
        CALL edgfill(zp  ,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz)

        IF ( istgr == 0 ) THEN
          DO j=1,ny
            DO i=1,nx
              DO k=1,nz-1
                tem2(i,j,k) = .5 * ( zp(i,j,k) + zp(i,j,k+1) )
              END DO
              tem2(i,j,nz) = zp(i,j,nz)
            END DO
          END DO
        ELSE
          DO j=1,ny
            DO i=1,nx
              DO k=1,nz
                tem2(i,j,k) = zp(i,j,k)
              END DO
            END DO
          END DO
        END IF

        DO k=1,nl(iv)
          DO i=1,nc
            DO j=1,nr
              tem1(j,i,k) = tem2(i,nr-j+1,k)
            END DO
          END DO
        END DO

      ELSE IF(varnam(iv) == 'U') THEN
        CALL v5dsetunits(iv, "m/s\0")
        CALL edgfill(u,   1,nx,1,nx, 1,ny,1,ny-1, 1,nz,1,nz-1)

        IF ( istgr == 0 ) THEN
          DO k=1,nz
            DO j=1,ny
              DO i=1,nx-1
                tem2(i,j,k) = .5 * ( u(i,j,k) + u(i+1,j,k) )
              END DO
              tem2(nx,j,k) = u(nx,j,k)
            END DO
          END DO
        ELSE
          DO k=1,nz
            DO j=1,ny
              DO i=1,nx
                tem2(i,j,k) = u(i,j,k)
              END DO
            END DO
          END DO
        END IF

        DO k=1,nl(iv)
          DO i=1,nc
            DO j=1,nr
              tem1(j,i,k) = tem2(i,nr-j+1,k)
            END DO
          END DO
        END DO

      ELSE IF(varnam(iv) == 'V') THEN
        CALL v5dsetunits(iv, "m/s\0")
        CALL edgfill(v,   1,nx,1,nx-1, 1,ny,1,ny, 1,nz,1,nz-1)

        IF ( istgr == 0 ) THEN
          DO k=1,nz
            DO i=1,nx
              DO j=1,ny-1
                tem2(i,j,k) = .5 * ( v(i,j,k) + v(i,j+1,k) )
              END DO
              tem2(i,ny,k) = v(i,ny,k)
            END DO
          END DO
        ELSE
          DO k=1,nz
            DO j=1,ny
              DO i=1,nx
                tem2(i,j,k) = v(i,j,k)
              END DO
            END DO
          END DO
        END IF

        DO k=1,nl(iv)
          DO i=1,nc
            DO j=1,nr
              tem1(j,i,k) = tem2(i,nr-j+1,k)
            END DO
          END DO
        END DO

      ELSE IF(varnam(iv) == 'UPRT') THEN
        CALL v5dsetunits(iv, "m/s\0")
        CALL edgfill(ubar,1,nx,1,nx, 1,ny,1,ny-1, 1,nz,1,nz-1)

        IF ( istgr == 0 ) THEN
          DO k=1,nz
            DO j=1,ny
              DO i=1,nx-1
                tem2(i,j,k) = .5 * ( u(i,j,k) + u(i+1,j,k) )            &
                            - .5 * ( ubar(i,j,k) + ubar(i+1,j,k) )
              END DO
              tem2(nx,j,k) = u(nx,j,k) - ubar(nx,j,k)
            END DO
          END DO
        ELSE
          DO k=1,nz
            DO j=1,ny
              DO i=1,nx
                tem2(i,j,k) = u(i,j,k) - ubar(i,j,k)
              END DO
            END DO
          END DO
        END IF

        DO k=1,nl(iv)
          DO i=1,nc
            DO j=1,nr
              tem1(j,i,k) = tem2(i,nr-j+1,k)
            END DO
          END DO
        END DO

      ELSE IF(varnam(iv) == 'VPRT') THEN
        CALL v5dsetunits(iv, "m/s\0")
        CALL edgfill(vbar,1,nx,1,nx-1, 1,ny,1,ny, 1,nz,1,nz-1)

        IF ( istgr == 0 ) THEN
          DO k=1,nz
            DO i=1,nx
              DO j=1,ny-1
                tem2(i,j,k) = .5 * ( v(i,j,k) + v(i,j+1,k) )            &
                            - .5 * ( vbar(i,j,k) + vbar(i,j+1,k) )
              END DO
              tem2(i,ny,k) = v(i,ny,k) - vbar(i,ny,k)
            END DO
          END DO
        ELSE
          DO k=1,nz
            DO j=1,ny
              DO i=1,nx
                tem2(i,j,k) = v(i,j,k) - vbar(i,j,k)
              END DO
            END DO
          END DO
        END IF

        DO k=1,nl(iv)
          DO i=1,nc
            DO j=1,nr
              tem1(j,i,k) = tem2(i,nr-j+1,k)
            END DO
          END DO
        END DO

      ELSE IF(varnam(iv) == 'W') THEN
        CALL v5dsetunits(iv, "m/s\0")
        CALL edgfill(w   ,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz)

        IF ( istgr == 0 ) THEN
          DO j=1,ny
            DO i=1,nx
              DO k=1,nz-1
                tem2(i,j,k) = .5 * ( w(i,j,k) + w(i,j,k+1) )
              END DO
              tem2(i,j,nz) = w(i,j,nz)
            END DO
          END DO
        ELSE
          DO j=1,ny
            DO i=1,nx
              DO k=1,nz
                tem2(i,j,k) = w(i,j,k)
              END DO
            END DO
          END DO
        END IF

        DO k=1,nl(iv)
          DO i=1,nc
            DO j=1,nr
              tem1(j,i,k) = tem2(i,nr-j+1,k)
            END DO
          END DO
        END DO

      ELSE IF(varnam(iv) == 'PT') THEN
        CALL v5dsetunits(iv, "K\0")
        CALL edgfill(ptprt,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)
        CALL edgfill(ptbar,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)
        DO k=1,nl(iv)
          DO i=1,nc
            DO j=1,nr
              tem1(j,i,k) = ptprt(i,nr-j+1,k)+ptbar(i,nr-j+1,k)
            END DO
          END DO
        END DO

      ELSE IF(varnam(iv) == 'PTPRT') THEN
        CALL v5dsetunits(iv, "K\0")
        DO k=1,nl(iv)
          DO i=1,nc
            DO j=1,nr
              tem1(j,i,k) = ptprt(i,nr-j+1,k)
            END DO
          END DO
        END DO

      ELSE IF(varnam(iv) == 'P') THEN
        CALL v5dsetunits(iv, "mb\0")
        CALL edgfill(pprt,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)
        CALL edgfill(pbar,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)
        DO k=1,nl(iv)
          DO i=1,nc
            DO j=1,nr
              tem1(j,i,k) = 0.01*(pprt(i,nr-j+1,k)+pbar(i,nr-j+1,k))
            END DO
          END DO
        END DO

      ELSE IF(varnam(iv) == 'PPRT') THEN
        CALL v5dsetunits(iv, "mb\0")
        DO k=1,nl(iv)
          DO i=1,nc
            DO j=1,nr
              tem1(j,i,k) = 0.01*pprt(i,nr-j+1,k)
            END DO
          END DO
        END DO

      ELSE IF(varnam(iv) == 'VORT') THEN
        CALL v5dsetunits(iv, "1/s\0")
        DO k=2,nz-2
          DO j=2,ny-2
            DO i=2,nx-2
              tem2(i,j,k)=                                              &
                  (v(i+1,j,k)-v(i-1,j,k)+v(i+1,j+1,k)-v(i-1,j+1,k))/    &
                  (4*(x(i+1)-x(i)))-                                    &
                  (u(i,j+1,k)-u(i,j-1,k)+u(i+1,j+1,k)-u(i+1,j-1,k))/    &
                  (4*(y(j+1)-y(j)))
            END DO
          END DO
        END DO

        DO j=2,ny-2
          DO i=2,nx-2
            tem2(i,j,   1)=tem2(i,j,   2)
            tem2(i,j,nz-1)=tem2(i,j,nz-2)
          END DO
        END DO

        DO k=1,nz-1
          DO j=2,ny-2
            tem2(   1,j,k)=tem2(   2,j,k)
            tem2(nx-1,j,k)=tem2(nx-2,j,k)
          END DO
        END DO

        DO k=1,nz-1
          DO i=1,nx-1
            tem2(i,   1,k)=tem2(i,   2,k)
            tem2(i,ny-1,k)=tem2(i,ny-2,k)
          END DO
        END DO

        CALL edgfill(tem2,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)

        DO k=1,nl(iv)
          DO i=1,nc
            DO j=1,nr
              tem1(j,i,k) = tem2(i,nr-j+1,k)
            END DO
          END DO
        END DO

      ELSE IF(varnam(iv) == 'DIV') THEN
        CALL v5dsetunits(iv, "1/s\0")
        DO k=1,nz-1
          DO j=1,ny-1
            DO i=1,nx-1
              tem2(i,j,k) = (u(i+1,j,k)-u(i,j,k))/(x(i+1)-x(i))         &
                          + (v(i,j+1,k)-v(i,j,k))/(y(j+1)-y(j))
            END DO
          END DO
        END DO

        CALL edgfill(tem2,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)
        DO k=1,nl(iv)
          DO i=1,nc
            DO j=1,nr
              tem1(j,i,k) = tem2(i,nr-j+1,k)
            END DO
          END DO
        END DO

      ELSE IF(varnam(iv) == 'QV') THEN
        CALL v5dsetunits(iv, "g/kg\0")
        CALL edgfill(qv,   1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)
        DO k=1,nl(iv)
          DO i=1,nc
            DO j=1,nr
              tem1(j,i,k) = 1000.*qv(i,nr-j+1,k)
            END DO
          END DO
        END DO

      ELSE IF(varnam(iv) == 'QVPRT') THEN
        CALL v5dsetunits(iv, "g/kg\0")
        CALL edgfill(qvbar,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)
        DO k=1,nl(iv)
          DO i=1,nc
            DO j=1,nr
              tem1(j,i,k) = 1000.*(qv(i,nr-j+1,k) - qvbar(i,j,k))
            END DO
          END DO
        END DO

      END IF

      DO nq=1,nscalar
        IF(varnam(iv) == TRIM(qnames_upper(nq))) THEN
          IF(varnam(iv)(1:1) == 'Q') THEN
            CALL v5dsetunits(iv, "g/kg\0")
            factor = 1000.
          ELSE IF(varnam(iv)(1:1) == 'N') THEN
            CALL v5dsetunits(iv, "#/m^3\0")
            factor = 1.
          ELSE IF(varnam(iv)(1:1) == 'Z') THEN
            CALL v5dsetunits(iv, "m^3\0")
            factor = 1.
          END IF
          CALL edgfill(qscalar(:,:,:,nq),1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)
          DO k=1,nl(iv)
            DO i=1,nc
              DO j=1,nr
                tem1(j,i,k) = factor*qscalar(i,nr-j+1,k,nq)
              END DO
            END DO
          END DO
        END IF
      END DO

!      ELSE IF(varnam(iv) == 'QC') THEN
!        CALL v5dsetunits(iv, "g/kg\0")
!        CALL edgfill(qscalar(:,:,:,P_QC),1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)
!        DO k=1,nl(iv)
!          DO i=1,nc
!            DO j=1,nr
!              tem1(j,i,k) = 1000.*qscalar(i,nr-j+1,k,P_QC)
!            END DO
!          END DO
!        END DO

!      ELSE IF(varnam(iv) == 'QR') THEN
!        CALL v5dsetunits(iv, "g/kg\0")
!        CALL edgfill(qscalar(:,:,:,P_QR),1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)
!        DO k=1,nl(iv)
!          DO i=1,nc
!            DO j=1,nr
!              tem1(j,i,k) = 1000.*qscalar(i,nr-j+1,k,P_QR)
!            END DO
!          END DO
!        END DO

!      ELSE IF(varnam(iv) == 'QI') THEN
!        CALL v5dsetunits(iv, "g/kg\0")
!        CALL edgfill(qscalar(:,:,:,P_QI),1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)
!        DO k=1,nl(iv)
!          DO i=1,nc
!            DO j=1,nr
!              tem1(j,i,k) = 1000.*qscalar(i,nr-j+1,k,P_QI)
!            END DO
!          END DO
!        END DO

!      ELSE IF(varnam(iv) == 'QS') THEN
!        CALL v5dsetunits(iv, "g/kg\0")
!        CALL edgfill(qscalar(:,:,:,P_QS),1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)
!        DO k=1,nl(iv)
!          DO i=1,nc
!            DO j=1,nr
!              tem1(j,i,k) = 1000.*qscalar(i,nr-j+1,k,P_QS)
!            END DO
!          END DO
!        END DO

!      ELSE IF(varnam(iv) == 'QH') THEN
!        CALL v5dsetunits(iv, "g/kg\0")
!        CALL edgfill(qscalar(:,:,:,P_QH),1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)
!        DO k=1,nl(iv)
!          DO i=1,nc
!            DO j=1,nr
!              tem1(j,i,k) = 1000.*qscalar(i,nr-j+1,k,P_QH)
!            END DO
!          END DO
!        END DO

      IF(varnam(iv) == 'RAING') THEN
        CALL v5dsetunits(iv, "mm\0")
        CALL edgfill(raing,1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)
        DO k=1,nl(iv)
          DO i=1,nc
            DO j=1,nr
              tem1(j,i,k) = raing(i,nr-j+1)
            END DO
          END DO
        END DO

      ELSE IF(varnam(iv) == 'RAINC') THEN
        CALL v5dsetunits(iv, "mm\0")
        CALL edgfill(rainc,1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)
        DO k=1,nl(iv)
          DO i=1,nc
            DO j=1,nr
              tem1(j,i,k) = rainc(i,nr-j+1)
            END DO
          END DO
        END DO

      ELSE IF(varnam(iv) == 'PRCRATE1') THEN
        CALL v5dsetunits(iv, "kg/m**2/s\0")
        CALL edgfill(prcrate(1,1,1),1,nx,1,nx-1,1,ny,1,ny-1,1,1,1,1)
        DO k=1,nl(iv)
          DO i=1,nc
            DO j=1,nr
              tem1(j,i,k) = prcrate(i,nr-j+1,1)
            END DO
          END DO
        END DO

      ELSE IF(varnam(iv) == 'PRCRATE2') THEN
        CALL v5dsetunits(iv, "kg/m**2/s\0")
        CALL edgfill(prcrate(1,1,2),1,nx,1,nx-1,1,ny,1,ny-1,1,1,1,1)
        DO k=1,nl(iv)
          DO i=1,nc
            DO j=1,nr
              tem1(j,i,k) = prcrate(i,nr-j+1,2)
            END DO
          END DO
        END DO

      ELSE IF(varnam(iv) == 'PRCRATE3') THEN
        CALL edgfill(prcrate(1,1,3),1,nx,1,nx-1,1,ny,1,ny-1,1,1,1,1)
        CALL v5dsetunits(iv, "kg/m**2/s\0")
        DO k=1,nl(iv)
          DO i=1,nc
            DO j=1,nr
              tem1(j,i,k) = prcrate(i,nr-j+1,3)
            END DO
          END DO
        END DO

      ELSE IF(varnam(iv) == 'PRCRATE4') THEN
        CALL v5dsetunits(iv, "kg/m**2/s\0")
        CALL edgfill(prcrate(1,1,4),1,nx,1,nx-1,1,ny,1,ny-1,1,1,1,1)
        DO k=1,nl(iv)
          DO i=1,nc
            DO j=1,nr
              tem1(j,i,k) = prcrate(i,nr-j+1,4)
            END DO
          END DO
        END DO

      ELSE IF(varnam(iv) == 'TKE') THEN
        CALL v5dsetunits(iv, "m**2/s\0")
        CALL edgfill(tke,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)
        DO k=1,nl(iv)
          DO i=1,nc
            DO j=1,nr
              tem1(j,i,k) = tke(i,nr-j+1,k)
            END DO
          END DO
        END DO

      ELSE IF(varnam(iv) == 'KMH') THEN
        CALL v5dsetunits(iv, "m**2/s\0")
        CALL edgfill(kmh,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)
        DO k=1,nl(iv)
          DO i=1,nc
            DO j=1,nr
              tem1(j,i,k) = kmh(i,nr-j+1,k)
            END DO
          END DO
        END DO

      ELSE IF(varnam(iv) == 'KMV') THEN
        CALL v5dsetunits(iv, "m**2/s\0")
        CALL edgfill(kmv,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)
        DO k=1,nl(iv)
          DO i=1,nc
            DO j=1,nr
              tem1(j,i,k) = kmv(i,nr-j+1,k)
            END DO
          END DO
        END DO

      END IF

      IF ( sfcout == 1 ) THEN

        IF (nstyp <= 1) THEN

          IF(varnam(iv) == 'TSOIL') THEN
            CALL v5dsetunits(iv, "K\0")
            CALL edgfill(tsoil(1,1,1,0),   1,nx,1,nx-1, 1,ny,1,ny-1,       &
                         1,nzsoil,1,nzsoil)
            DO k=1,nl(iv)
              DO i=1,nc
                DO j=1,nr
!                  tem1(j,i,k) = tsoil(i,nr-j+1,0)       ! Origin, ?-- WYH
                  tem1(j,i,k) = tsoil(i,nr-j+1,k,0)
                END DO
              END DO
            END DO

          ELSE IF(varnam(iv) == 'QSOIL') THEN
            CALL v5dsetunits(iv, "m**3/m**3\0")
            CALL edgfill(qsoil(1,1,1,0),   1,nx,1,nx-1, 1,ny,1,ny-1,     &
                         1,nzsoil,1,nzsoil)
            DO k=1,nl(iv)
              DO i=1,nc
                DO j=1,nr
!                  tem1(j,i,k) = qsoil(i,nr-j+1,0)       ! Origin, ?-- WYH
                  tem1(j,i,k) = qsoil(i,nr-j+1,k,0)
                END DO
              END DO
            END DO

          ELSE IF(varnam(iv) == 'WR') THEN
            CALL v5dsetunits(iv, "m\0")
            CALL edgfill(wetcanp(1,1,0),   1,nx,1,nx-1, 1,ny,1,ny-1,    &
                         1,1,1,1)
            DO k=1,nl(iv)
              DO i=1,nc
                DO j=1,nr
                  tem1(j,i,k) = wetcanp(i,nr-j+1,0)
                END DO
              END DO
            END DO

          END IF

        ELSE

          DO is=0,nstyp

            IF(varnam(iv)(1:5) == 'TSOIL') THEN
!              READ(varnam(iv)(6:6), '(i1)') l
              CALL v5dsetunits(iv, "K\0")
              CALL edgfill(tsoil(1,1,1,is),   1,nx,1,nx-1, 1,ny,1,ny-1,       &
                           1,nzsoil,1,nzsoil)
              DO k=1,nl(iv)
                DO i=1,nc
                  DO j=1,nr
                    tem1(j,i,k) = tsoil(i,nr-j+1,k,is)
                  END DO
                END DO
              END DO

            ELSE IF(varnam(iv)(1:5) == 'QSOIL') THEN
!              READ(varnam(iv)(6:6), '(i1)') l
              CALL v5dsetunits(iv, "m**3/m**3\0")
              CALL edgfill(qsoil(1,1,1,is),   1,nx,1,nx-1, 1,ny,1,ny-1,     &
                           1,nzsoil,1,nzsoil)
              DO k=1,nl(iv)
                DO i=1,nc
                  DO j=1,nr
                    tem1(j,i,k) = qsoil(i,nr-j+1,k,is)
                  END DO
                END DO
              END DO

            ELSE IF(varnam(iv)(1:2) == 'WR') THEN
!              READ(varnam(iv)(3:3), '(i1)') l
              CALL v5dsetunits(iv, "m\0")
              CALL edgfill(wetcanp(1,1,is),   1,nx,1,nx-1, 1,ny,1,ny-1,    &
                           1,1,1,1)
              DO k=1,nl(iv)
                DO i=1,nc
                  DO j=1,nr
                    tem1(j,i,k) = wetcanp(i,nr-j+1,is)
                  END DO
                END DO
              END DO

            END IF

          END DO   ! is

        END IF

        IF(snowout == 1) THEN
          IF (varnam(iv)(1:2) == 'SD') THEN
            CALL v5dsetunits(iv, "index\0")
            CALL edgfill(snowdpth,1,nx,1,nx-1, 1,ny,1,ny-1,1,1,1,1)
            DO k=1,nl(iv)
              DO i=1,nc
                DO j=1,nr
                  tem1(j,i,k) = snowdpth(i,nr-j+1)
                END DO
              END DO
            END DO
          END IF
        END IF

      END IF       ! End sfcout

      IF(landout == 1) THEN

        IF (nstyp <= 1) THEN

          IF(varnam(iv) == 'SOILTY') THEN
            CALL v5dsetunits(iv, "index\0")
            CALL iedgfill(soiltyp(1,1,1),1,nx,1,nx-1, 1,ny,1,ny-1,      &
                         1,1,1,1)
            DO k=1,nl(iv)
              DO i=1,nc
                DO j=1,nr
                  tem1(j,i,k) = soiltyp(i,nr-j+1,1)
                END DO
              END DO
            END DO
          END IF

        ELSE

          DO is = 1, nstyp

            IF(varnam(iv)(1:6) == 'SOILTY') THEN
              !READ(varnam(iv)(7:7), '(i1)') is
              CALL v5dsetunits(iv, "index\0")
              CALL iedgfill(soiltyp(1,1,is), 1,nx,1,nx-1, 1,ny,1,ny-1,    &
                           1,1,1,1)
              DO k=1,nl(iv)
                DO i=1,nc
                  DO j=1,nr
                    tem1(j,i,k) = soiltyp(i,nr-j+1,is)
                  END DO
                END DO
              END DO

            ELSE IF(varnam(iv)(1:6) == 'SOILFR') THEN
              !READ(varnam(iv)(7:7), '(i1)') is
              CALL v5dsetunits(iv, "fraction\0")
              CALL edgfill(stypfrct(1,1,is),   1,nx,1,nx-1, 1,ny,1,ny-1,   &
                           1,1,1,1)
              DO k=1,nl(iv)
                DO i=1,nc
                  DO j=1,nr
                    tem1(j,i,k) = stypfrct(i,nr-j+1,is)
                  END DO
                END DO
              END DO

            END IF

          END DO

        END IF

        IF(varnam(iv) == 'VEGTYP') THEN
          CALL v5dsetunits(iv, "index\0")
          CALL iedgfill(vegtyp ,1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)
          DO k=1,nl(iv)
            DO i=1,nc
              DO j=1,nr
                tem1(j,i,k) = vegtyp(i,nr-j+1)
              END DO
            END DO
          END DO

        ELSE IF(varnam(iv) == 'LAI') THEN
          CALL edgfill(lai ,1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)
          DO k=1,nl(iv)
            DO i=1,nc
              DO j=1,nr
                tem1(j,i,k) = lai(i,nr-j+1)
              END DO
            END DO
          END DO

        ELSE IF(varnam(iv) == 'ROUFNS') THEN
          CALL v5dsetunits(iv, "m\0")
          CALL edgfill(roufns ,1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)
          DO k=1,nl(iv)
            DO i=1,nc
              DO j=1,nr
                tem1(j,i,k) = roufns(i,nr-j+1)
              END DO
            END DO
          END DO

        ELSE IF(varnam(iv) == 'VEG') THEN
          CALL v5dsetunits(iv, "fraction\0")
          CALL edgfill(veg ,1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)
          DO k=1,nl(iv)
            DO i=1,nc
              DO j=1,nr
                tem1(j,i,k) = veg(i,nr-j+1)
              END DO
            END DO
          END DO
        END IF

      END IF     ! End landout

      IF ( varnam(iv) == 'RADFRC' ) THEN
        CALL v5dsetunits(iv, "K/s\0")
        CALL edgfill(radfrc ,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)
        DO k=1,nl(iv)
          DO i=1,nc
            DO j=1,nr
              tem1(j,i,k) = radfrc(i,nr-j+1,k)
            END DO
          END DO
        END DO

      ELSE IF ( varnam(iv) == 'RADSW' ) THEN
        CALL v5dsetunits(iv, "W/m**2\0")
        CALL edgfill(radsw ,1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)
        DO k=1,nl(iv)
          DO i=1,nc
            DO j=1,nr
              tem1(j,i,k) = radsw(i,nr-j+1)
            END DO
          END DO
        END DO
      ELSE IF ( varnam(iv) == 'RNFLX' ) THEN
        CALL v5dsetunits(iv, "W/m**2\0")
        CALL edgfill(rnflx ,1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)
        DO k=1,nl(iv)
          DO i=1,nc
            DO j=1,nr
              tem1(j,i,k) = rnflx(i,nr-j+1)
            END DO
          END DO
        END DO

      ELSE IF ( varnam(iv) == 'USFLX' ) THEN
        CALL v5dsetunits(iv, "kg/(m*s**2)\0")
        CALL edgfill(usflx ,1,nx,1,nx, 1,ny,1,ny-1, 1,1,1,1)
        DO k=1,nl(iv)
          DO i=1,nc
            DO j=1,nr
              tem1(j,i,k) = usflx(i,nr-j+1)
            END DO
          END DO
        END DO

      ELSE IF ( varnam(iv) == 'VSFLX' ) THEN
        CALL v5dsetunits(iv, "kg/(m*s**2)\0")
        CALL edgfill(vsflx ,1,nx,1,nx-1, 1,ny,1,ny, 1,1,1,1)
        DO k=1,nl(iv)
          DO i=1,nc
            DO j=1,nr
              tem1(j,i,k) = vsflx(i,nr-j+1)
            END DO
          END DO
        END DO

      ELSE IF ( varnam(iv) == 'PTSFLX' ) THEN
        CALL v5dsetunits(iv, "K*kg/(m**2*s)\0")
        CALL edgfill(ptsflx ,1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)
        DO k=1,nl(iv)
          DO i=1,nc
            DO j=1,nr
              tem1(j,i,k) = ptsflx(i,nr-j+1)
            END DO
          END DO
        END DO

      ELSE IF ( varnam(iv) == 'QVSFLX' ) THEN
        CALL v5dsetunits(iv, "kg/(m**2*s)\0")
        CALL edgfill(qvsflx ,1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)
        DO k=1,nl(iv)
          DO i=1,nc
            DO j=1,nr
              tem1(j,i,k) = qvsflx(i,nr-j+1)
            END DO
          END DO
        END DO

!      ELSE

!        WRITE(6,'(3a,/)') 'Unknown variable name - ',varnam(iv),' - inside v5dump.'
!        CALL EXIT(1)

      END IF

      WRITE(6,'(2a,/)') '  Writing varaible: ',varnam(iv)

      m = v5dwrite(it,iv,tem1)
      IF(m == 0) CALL EXIT(1)

    END DO   ! enddo of iv
  END DO   ! enddo of it
!
!----------------------------------------------------------------------
!
!  End of putput, close the file & return
!
!-----------------------------------------------------------------------

  m = v5dclose()
  IF(m == 0) CALL EXIT(1)

  RETURN
END SUBROUTINE v5ddump
