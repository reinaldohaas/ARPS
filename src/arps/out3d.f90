!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE INITOUT                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE initout(mptr,nx,ny,nz,nzsoil,nstyps,exbcbufsz,               &
           u,v,w,wcont,ptprt,pprt,qv,qscalar,tke,                       &
           ubar,vbar,ptbar,pbar,rhostr,qvbar,kmh,kmv,                   &
           x,y,z,zp,zpsoil,hterain, mapfct, j1,j2,j3,j3inv,             &
           soiltyp,stypfrct,vegtyp,lai,roufns,veg,                      &
           tsoil,qsoil,wetcanp,snowdpth,                                &
           raing,rainc,prcrate,exbcbuf,                                 &
           radfrc,radsw,rnflx,radswnet,radlwin,                         &
           usflx,vsflx,ptsflx,qvsflx,                                   &
           rhobar,tem1,tem2,tem3,tem4,tem5, tem6,                       &
           tem7,tem8,tem9,tem10)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Handles the model data output at the initial time.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  11/13/91.
!
!  MODIFICATION HISTORY:
!
!  6/06/92 (M. Xue)
!  Added full documentation.
!
!  7/13/92 (K. Brewster)
!  Added comment line variables to arg list to be passed to the
!  data dumping routines.
!
!  7/20/92 (M. Xue)
!  Added call to energy and base state dump routines.
!
!  1/24/94 (Y. Liu)
!  Added surface variables to the data dumping routines.
!
!  12/09/1998 (Donghai Wang)
!  Added the snow cover.
!
!  05/14/2002 (J. Brotzge)
!  Added variables, modified call statements to allow for multiple
!  soil schemes
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    mptr     Grid identifier.
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!    nzsoil   Number of grid points in the soil
!
!    u        x component of velocity (m/s)
!    v        y component of velocity (m/s)
!    w        Vertical component of Cartesian velocity (m/s)
!    wcont    Contravariant vertical velocity (m/s)
!    ptprt    Perturbation potential temperature (K)
!    pprt     Perturbation pressure (Pascal)
!    qv       Water vapor specific humidity (kg/kg)
!    qc       Cloud water mixing ratio (kg/kg)
!    qr       Rainwater mixing ratio (kg/kg)
!    qi       Cloud ice mixing ratio (kg/kg)
!    qs       Snow mixing ratio (kg/kg)
!    qh       Hail mixing ratio (kg/kg)
!    tke      Turbulent Kinetic Energy ((m/s)**2)
!
!    ubar     Base state zonal velocity component (m/s)
!    vbar     Base state meridional velocity component (m/s)
!    ptbar    Base state potential temperature (K)
!    pbar     Base state pressure (Pascal)
!    rhostr   Base state density * j3 (kg/m**3)
!    qvbar    Base state water vapor specific humidity (kg/kg)
!
!    kmh      Horizontal turb. mixing coef. for momentum ( m**2/s )
!    kmv      Vertical turb. mixing coef. for momentum ( m**2/s )
!
!    x        x coordinate of grid points in physical/comp. space (m)
!    y        y coordinate of grid points in physical/comp. space (m)
!    z        z coordinate of grid points in computational space (m)
!    zp       Vertical coordinate of grid points in physical space (m)
!    zpsoil   Vertical coordinate of grid points in the soil (m)
!    hterain  Terrain height (m)
!    mapfct   Map factors at scalar, u and v points
!
!    j1       Coordinate transformation Jacobian -d(zp)/dx
!    j2       Coordinate transformation Jacobian -d(zp)/dy
!    j3       Coordinate transformation Jacobian  d(zp)/dz
!
!    soiltyp  Soil type
!    stypfrct  Soil type fraction
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
!    prcrate  Precipitation rates
!
!    radfrc   Radiation forcing (K/s)
!    radsw    Solar radiation reaching the surface
!    rnflx    Net radiation flux absorbed by surface
!    radswnet Net solar radiation
!    radlwin  Incoming longwave radiation
!
!    usflx    Surface flux of u-momentum (kg/(m*s**2))
!    vsflx    Surface flux of v-momentum (kg/(m*s**2))
!    ptsflx   Surface heat flux (K*kg/(m**2 * s ))
!    qvsflx   Surface moisture flux of (kg/(m**2 * s))
!
!  OUTPUT:
!
!    None.
!
!  WORK ARRAYS:
!
!    rhobar   Base state density (kg/m**3)
!    tem1     Temporary work array.
!    tem2     Temporary work array.
!    tem3     Temporary work array.
!    tem4     Temporary work array.
!    tem5     Temporary work array.
!    tem6     Temporary work array.
!    tem7     Temporary work array.
!    tem8     Temporary work array.
!    tem9     Temporary work array.
!    tem10    Temporary work array.
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'phycst.inc'
  INCLUDE 'globcst.inc'
  INCLUDE 'grid.inc'          ! Grid & map parameters.
  INCLUDE 'bndry.inc'
  INCLUDE 'exbc.inc'
  INCLUDE 'mp.inc'            ! Message passing parameters.
  INCLUDE 'timelvls.inc'

!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  INTEGER :: mptr              ! Grid identifier.


  INTEGER :: nx,ny,nz          ! Number of grid points in 3 directions
  INTEGER :: nzsoil            ! Number of grid points in the soil


  REAL :: u     (nx,ny,nz,nt)  ! Total u-velocity (m/s)
  REAL :: v     (nx,ny,nz,nt)  ! Total v-velocity (m/s)
  REAL :: w     (nx,ny,nz,nt)  ! Total w-velocity (m/s)
  REAL :: wcont (nx,ny,nz)     ! Contravariant vertical velocity (m/s)
  REAL :: ptprt (nx,ny,nz,nt)  ! Perturbation potential temperature (K)
  REAL :: pprt  (nx,ny,nz,nt)  ! Perturbation pressure (Pascal)

  REAL :: qv    (nx,ny,nz,nt)  ! Water vapor specific humidity (kg/kg)

  REAL :: qscalar(nx,ny,nz,nt,nscalar)

  REAL :: tke   (nx,ny,nz,nt)  ! Turbulent Kinetic Energy ((m/s)**2)

  REAL :: ubar  (nx,ny,nz)     ! Base state u-velocity (m/s)
  REAL :: vbar  (nx,ny,nz)     ! Base state v-velocity (m/s)
  REAL :: ptbar (nx,ny,nz)     ! Base state potential temperature (K)
  REAL :: pbar  (nx,ny,nz)     ! Base state pressure (Pascal)
  REAL :: rhostr(nx,ny,nz)     ! Base state air density * j3 (kg/m**3)
  REAL :: qvbar (nx,ny,nz)     ! Base state water vapor specific humidity
                               ! (kg/kg)

  REAL :: kmh   (nx,ny,nz)     ! Horizontal turb. mixing coef. for
                               ! momentum. ( m**2/s )
  REAL :: kmv   (nx,ny,nz)     ! Vertical turb. mixing coef. for
                               ! momentum. ( m**2/s )

  REAL :: x     (nx)           ! The x-coord. of the physical and
                               ! computational grid. Defined at u-point.
  REAL :: y     (ny)           ! The y-coord. of the physical and
                               ! computational grid. Defined at v-point.
  REAL :: z     (nz)           ! The z-coord. of the computational grid.
                               ! Defined at w-point on the staggered grid.
  REAL :: zp    (nx,ny,nz)     ! The physical height coordinate defined at
                               ! w-point of the staggered grid.
  REAL :: zpsoil(nx,ny,nzsoil) ! The physical height coordinate defined at
                               ! w-point of the soil

  REAL :: hterain(nx,ny)       ! Terrain height.
  REAL :: mapfct(nx,ny,8)      ! Map factors at scalar, u and v points

  REAL :: j1    (nx,ny,nz)     ! Coordinate transformation Jacobian defined
                               ! as - d( zp )/d( x )
  REAL :: j2    (nx,ny,nz)     ! Coordinate transformation Jacobian defined
                               ! as - d( zp )/d( y )
  REAL :: j3    (nx,ny,nz)     ! Coordinate transformation Jacobian defined
                               ! as d( zp )/d( z )
  REAL :: j3inv (nx,ny,nz)     ! Coordinate transformation Jacobian defined
                               ! as d( zp )/d( z )

  REAL :: rhobar(nx,ny,nz)     ! Base state air density (kg/m**3)

  INTEGER :: nstyps                 ! Number of soil types
  INTEGER :: soiltyp (nx,ny,nstyps) ! Soil type
  REAL    :: stypfrct(nx,ny,nstyps) ! Soil type fraction
  INTEGER :: vegtyp (nx,ny)         ! Vegetation type
  REAL    :: lai    (nx,ny)            ! Leaf Area Index
  REAL    :: roufns (nx,ny)            ! Surface roughness
  REAL    :: veg    (nx,ny)            ! Vegetation fraction

  REAL :: tsoil (nx,ny,nzsoil,0:nstyps) ! Soil temperature (K)
  REAL :: qsoil (nx,ny,nzsoil,0:nstyps) ! Soil moisture (m**3/m**3)
  REAL :: wetcanp(nx,ny,0:nstyps)      ! Canopy water amount
  REAL :: snowdpth(nx,ny)              ! Snow depth (m)

  REAL :: raing(nx,ny)         ! Grid supersaturation rain
  REAL :: rainc(nx,ny)         ! Cumulus convective rain
  REAL :: prcrate(nx,ny,4)     ! precipitation rates (kg/(m**2*s))
                               ! prcrate(1,1,1) = total precip. rate
                               ! prcrate(1,1,2) = grid scale precip. rate
                               ! prcrate(1,1,3) = cumulus precip. rate
                               ! prcrate(1,1,4) = microphysics precip. rate

  REAL :: radfrc(nx,ny,nz)     ! Radiation forcing (K/s)
  REAL :: radsw (nx,ny)        ! Solar radiation reaching the surface
  REAL :: rnflx (nx,ny)        ! Net radiation flux absorbed by surface
  REAL :: radswnet(nx,ny)      ! Net solar radiation
  REAL :: radlwin(nx,ny)       ! Incoming longwave radiation

  REAL :: usflx (nx,ny)        ! Surface flux of u-momentum (kg/(m*s**2))
  REAL :: vsflx (nx,ny)        ! Surface flux of v-momentum (kg/(m*s**2))
  REAL :: ptsflx(nx,ny)        ! Surface heat flux (K*kg/(m*s**2))
  REAL :: qvsflx(nx,ny)        ! Surface moisture flux (kg/(m**2*s))

  INTEGER :: exbcbufsz         ! Size of EXBC buffer
  REAL :: exbcbuf(exbcbufsz)   ! External boundary arrays

  REAL :: tem1  (nx,ny,nz)     ! Temporary work array
  REAL :: tem2  (nx,ny,nz)     ! Temporary work array
  REAL :: tem3  (nx,ny,nz)     ! Temporary work array
  REAL :: tem4  (nx,ny,nz)     ! Temporary work array
  REAL :: tem5  (nx,ny,nz)     ! Temporary work array
  REAL :: tem6  (nx,ny,nz)     ! Temporary work array
  REAL :: tem7  (nx,ny,nz)     ! Temporary work array
  REAL :: tem8  (nx,ny,nz)     ! Temporary work array
  REAL :: tem9  (nx,ny,nz)     ! Temporary work array
  REAL :: tem10 (nx,ny,nz)     ! Temporary work array
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: tlevel ,tim       ! Time level at which data are printed.
  CHARACTER (LEN=256) :: basdmpfn
  CHARACTER (LEN=256) :: exbcdmpfn
  INTEGER :: lexdmpf
  INTEGER :: lbasdmpf
  INTEGER :: grdbas
  REAL    :: amin, amax
  INTEGER :: i,j,k, nq
  INTEGER :: hdmpfmt1
  CHARACTER (LEN=256) :: ternfn,sfcoutfl,soiloutfl,temchar
  INTEGER :: lternfn,lfn
  INTEGER :: istatus

  CHARACTER (LEN=256) :: savename, outdirname
  CHARACTER (LEN=80 ) :: timsnd
  INTEGER :: tmstrln, oldirnam
  INTEGER :: nunit

!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  DO k=1,nz-1
    DO j=1,ny-1
      DO i=1,nx-1
        rhobar(i,j,k) = rhostr(i,j,k)*j3inv(i,j,k)
      END DO
    END DO
  END DO

  mgrid = mptr
!
!-----------------------------------------------------------------------
!
!  For the initial time step, print out the max/min of all varaibles
!  to check if their values are correct.
!
!-----------------------------------------------------------------------
!
  IF(myproc == 0) &
    WRITE(6,'(/1x,a,i2/)') &
      'Min. and max. of data arrays at initial time on grid ',mgrid

  CALL a3dmax0(x,1,nx,1,nx,1,1,1,1, 1,1,1,1, amax,amin)
  IF (myproc == 0) &
    WRITE(6,'(/1x,2(a,e13.6))') 'xmin    = ', amin,',  xmax    =',amax

  CALL a3dmax0(y,1,ny,1,ny,1,1,1,1, 1,1,1,1, amax,amin)
  IF (myproc == 0)  &
    WRITE(6,'(1x,2(a,e13.6))') 'ymin    = ', amin,',  ymax    =',amax

  CALL a3dmax0(z,1,nz,1,nz,1,1,1,1, 1,1,1,1, amax,amin)
  IF (myproc == 0)  &
    WRITE(6,'(1x,2(a,e13.6))') 'zmin    = ', amin,',  zmax    =',amax

  CALL a3dmax0(zp,1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz, amax,amin)
  IF (myproc == 0)  &
    WRITE(6,'(1x,2(a,e13.6))') 'zpmin   = ', amin,',  zpmax   =',amax

  CALL a3dmax0(hterain,1,nx,1,nx-1,1,ny,1,ny-1,1,1,1,1, amax,amin)
  IF (myproc == 0)  &
    WRITE(6,'(1x,2(a,e13.6))') 'hmin    = ', amin,',  hmax    =',amax

  CALL a3dmax0(ubar,1,nx,1,nx,1,ny,1,ny-1,1,nz,1,nz-1, amax,amin)
  IF (myproc == 0)  &
    WRITE(6,'(1x,2(a,e13.6))') 'ubarmin = ', amin,',  ubarmax =',amax

  CALL a3dmax0(vbar,1,nx,1,nx-1,1,ny,1,ny,1,nz,1,nz-1, amax,amin)
  IF (myproc == 0)  &
    WRITE(6,'(1x,2(a,e13.6))') 'vbarmin = ', amin,',  vbarmax =',amax

  CALL a3dmax0(ptbar,1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz-1,amax,amin)
  IF (myproc == 0)  &
    WRITE(6,'(1x,2(a,e13.6))') 'ptbarmin= ', amin,',  ptbarmax=',amax

  CALL a3dmax0(pbar,1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz-1, amax,amin)
  IF (myproc == 0)  &
    WRITE(6,'(1x,2(a,e13.6))') 'pbarmin = ', amin,',  pbarmax =',amax

  CALL a3dmax0(qvbar,1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz-1,amax,amin)
  IF (myproc == 0)  &
    WRITE(6,'(1x,2(a,e13.6))') 'qvbarmin= ', amin,',  qvbarmax=',amax


  IF (myproc == 0)  &
    WRITE(6,'(/1x,a/)') 'Min/max of fields at tpresent:'
  tim = tpresent

  CALL a3dmax0(u(1,1,1,tim),1,nx,1,nx,1,ny,1,ny-1,1,nz,1,nz-1,          &
               amax,amin)
  IF (myproc == 0)  &
    WRITE(6,'(1x,2(a,e13.6))') 'umin    = ', amin,',  umax    =',amax

  CALL a3dmax0(v(1,1,1,tim),1,nx,1,nx-1,1,ny,1,ny,1,nz,1,nz-1,          &
               amax,amin)
  IF (myproc == 0)  &
    WRITE(6,'(1x,2(a,e13.6))') 'vmin    = ', amin,',  vmax    =',amax

  CALL a3dmax0(w(1,1,1,tim),1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz,          &
               amax,amin)
  IF (myproc == 0)  &
    WRITE(6,'(1x,2(a,e13.6))') 'wmin    = ', amin,',  wmax    =',amax

  CALL a3dmax0(ptprt(1,1,1,tim),1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz-1,    &
               amax,amin)
  IF (myproc == 0)  &
    WRITE(6,'(1x,2(a,e13.6))') 'ptprtmin= ', amin,',  ptprtmax=',amax

  CALL a3dmax0(pprt(1,1,1,tim),1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz-1,     &
               amax,amin)
  IF (myproc == 0)  &
    WRITE(6,'(1x,2(a,e13.6))') 'pprtmin = ', amin,',  pprtmax =',amax

  CALL a3dmax0(qv(1,1,1,tim),1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz-1,       &
               amax,amin)
  IF (myproc == 0)  &
    WRITE(6,'(1x,2(a,e13.6))') 'qvmin   = ', amin,',  qvmax   =',amax

  DO nq = 1,nscalar

    CALL a3dmax0(qscalar(:,:,:,tim,nq),1,nx,1,nx-1,1,ny,1,ny-1,         &
                 1,nz,1,nz-1,amax,amin)

    IF (myproc == 0)  &
      WRITE(6,'(1x,2(a,e13.6))') TRIM(qnames(nq))//'min   = ', amin,    &
                          ',  '//TRIM(qnames(nq))//'max   =',  amax

  END DO

  IF (myproc == 0)  &
    WRITE(6,'(/1x,a/)') 'Min/max of fields at tpast:'

  tim = tpast

  CALL a3dmax0(u(1,1,1,tim),1,nx,1,nx,1,ny,1,ny-1,1,nz,1,nz-1,          &
               amax,amin)
  IF (myproc == 0)  &
    WRITE(6,'(1x,2(a,e13.6))') 'umin    = ', amin,',  umax    =',amax

  CALL a3dmax0(v(1,1,1,tim),1,nx,1,nx-1,1,ny,1,ny,1,nz,1,nz-1,          &
               amax,amin)
  IF (myproc == 0)  &
    WRITE(6,'(1x,2(a,e13.6))') 'vmin    = ', amin,',  vmax    =',amax

  CALL a3dmax0(w(1,1,1,tim),1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz,          &
               amax,amin)
  IF (myproc == 0)  &
    WRITE(6,'(1x,2(a,e13.6))') 'wmin    = ', amin,',  wmax    =',amax

  CALL a3dmax0(ptprt(1,1,1,tim),1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz-1,    &
               amax,amin)
  IF (myproc == 0)  &
    WRITE(6,'(1x,2(a,e13.6))') 'ptprtmin= ', amin,',  ptprtmax=',amax

  CALL a3dmax0(pprt(1,1,1,tim),1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz-1,     &
               amax,amin)
  IF (myproc == 0)  &
    WRITE(6,'(1x,2(a,e13.6))') 'pprtmin = ', amin,',  pprtmax =',amax

  CALL a3dmax0(qv(1,1,1,tim),1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz-1,       &
               amax,amin)
  IF (myproc == 0)  &
    WRITE(6,'(1x,2(a,e13.6))') 'qvmin   = ', amin,',  qvmax   =',amax

  DO nq = 1,nscalar

    CALL a3dmax0(qscalar(:,:,:,tim,nq),1,nx,1,nx-1,1,ny,1,ny-1,         &
                 1,nz,1,nz-1,amax,amin)

    IF (myproc == 0)  &
      WRITE(6,'(1x,2(a,e13.6))') TRIM(qnames(nq))//'min   = ', amin,    &
                          ',  '//TRIM(qnames(nq))//'max   =',  amax

  END DO

  tlevel = tpresent
!
!-----------------------------------------------------------------------
!
!  Calculation of the max./min. statistics and printing of initial fields
!
!-----------------------------------------------------------------------
!
  IF( nmaxmin > 0 ) THEN

    CALL maxmin(mptr,nx,ny,nz,nzsoil,tlevel,rhobar,                     &
                u,v,w,wcont,ptprt,pprt,qv,qscalar,tke,kmh,kmv,          &
                x,y,z,zp,zpsoil,mapfct,                                 &
                tsoil(1,1,1,0),qsoil(1,1,1,0),                          &
                wetcanp(1,1,0),                                         &
                tem1,tem2,tem3,tem4)

  END IF

!
!-----------------------------------------------------------------------
!
!  Calculation of the energy statistics and printing at the initial time
!
!-----------------------------------------------------------------------
!
  IF( nenergy > 0 ) THEN

    CALL energy(nx,ny,nz,                                               &
                u,v,w,ptprt,pprt,qv,qscalar,rhobar,                     &
                x,y,z,zp,hterain, j3,                                   &
                tem1,tem2,tem3,tem4,tem5)

  END IF

  IF( nfmtprt > 0 ) THEN
!
!-----------------------------------------------------------------------
!
!  Formatted printing of base state variables
!
!-----------------------------------------------------------------------
!
    CALL basprt(nx,ny,nz,ubar,vbar,ptbar,pbar,rhobar,qvbar,             &
                zp,hterain)
!
!-----------------------------------------------------------------------
!
!  Formatted printing of time dependent fields
!
!-----------------------------------------------------------------------
!
    IF (myproc == 0) THEN
      WRITE(6,'(/'' THE INITIAL FIELDS:''/)')

      CALL fmtprt(nx,ny,nz, tlevel,                                     &
                  u,v,w,ptprt,pprt,qv,qscalar,tke,kmh,kmv,              &
                  x,y,z,zp,hterain, j1,j2,j3, tem1)
    END IF
  END IF

!
!-----------------------------------------------------------------------
!
!  Check to see if any history data dump is to be produced.
!
!-----------------------------------------------------------------------
!

  IF( hdmpfmt == 0 .OR. nhisdmp <= 0 ) THEN

    IF (myproc == 0) WRITE(6,'(1x,a)')  &
         'History data dump option=0, no data is dumped.'
    GO TO 800

  END IF

!
!-----------------------------------------------------------------------
!
!  Data dump of the model grid and base state arrays:
!
!-----------------------------------------------------------------------
!
  IF( hdmpfmt == 5 )  GOTO 700
  IF( hdmpfmt == 9 )  GOTO 700
  IF( hdmpfmt == 11 ) GOTO 700

  IF (grdbasfout == 0) GOTO 700  ! Skip grid & base file output
!
!-----------------------------------------------------------------------
!
!  Find a unique name basdmpfn(1:lbasdmpf) for the grid and
!  base state array dump file
!
!-----------------------------------------------------------------------
!
  CALL get_output_dirname(1,dirname,-1.0E-12,1,outdirname,istatus)
  oldirnam = LEN_TRIM(outdirname)

  CALL gtbasfn(runname(1:lfnkey),outdirname,oldirnam,hdmpfmt,           &
               mgrid,nestgrd,basdmpfn,lbasdmpf)

  IF (myproc == 0) WRITE(6,'(1x,a,a)')                                  &
       'Data dump of grid and base state arrays into file ',            &
        basdmpfn(1:lbasdmpf)

  grdbas = 1                ! Dump out grid and base state arrays only

  tim = tpresent

  DO k=1,nz
    DO j=1,ny
      DO i=1,nx
        tem1(i,j,k)=0.0
      END DO
    END DO
  END DO

!    blocking inserted for ordering i/o for message passing
  DO i=0,nprocs-1,dumpstride
    IF(myproc >= i .AND. myproc <= i+dumpstride-1)THEN

      CALL dtadump(nx,ny,nz,nzsoil,nstyps,                              &
                 hdmpfmt,nchdmp,basdmpfn(1:lbasdmpf),                   &
                 grdbas,filcmprs,                                       &
                 u(1,1,1,tim),v(1,1,1,tim),                             &
                 w(1,1,1,tim),ptprt(1,1,1,tim),                         &
                 pprt(1,1,1,tim),qv(1,1,1,tim),                         &
                 qscalar(:,:,:,tim,:),tke(1,1,1,tim),kmh,kmv,           &
                 ubar,vbar,tem1,ptbar,pbar,rhobar,qvbar,                &
                 x,y,z,zp,zpsoil,                                       &
                 soiltyp,stypfrct,vegtyp,lai,roufns,veg,                &
                 tsoil,qsoil,wetcanp,snowdpth,                          &
                 raing,rainc,prcrate,                                   &
                 radfrc,radsw,rnflx,radswnet,radlwin,                   &
                 usflx,vsflx,ptsflx,qvsflx,                             &
                 tem2,tem3,tem6)


    END IF
    IF (mp_opt > 0) CALL mpbarrier
  END DO

  700   CONTINUE
!
!-----------------------------------------------------------------------
!
!  Data dump of the model time dependent arrays at initial time:
!
!-----------------------------------------------------------------------
!
  CALL get_output_dirname(1,dirname,curtim,1,outdirname,istatus)
  oldirnam = LEN_TRIM(outdirname)
!
!-----------------------------------------------------------------------
!
!  Find a unique name hdmpfn(1:ldmpf) for history dump data set
!  at time 'curtim'.
!
!-----------------------------------------------------------------------
!
!
  CALL gtdmpfn(runname(1:lfnkey),outdirname,                            &
               oldirnam,curtim,hdmpfmt,                                 &
               mgrid,nestgrd, hdmpfn, ldmpf)
  IF (myproc == 0)  &
    WRITE(6,'(1x,a,a)') 'History data dump in file ',hdmpfn(1:ldmpf)

  grdbas = 0                ! No base state or grid array is dumped.

  tim = tpresent

  DO k=1,nz
    DO j=1,ny
      DO i=1,nx
        tem1(i,j,k)=0.0
      END DO
    END DO
  END DO


!  blocking inserted for ordering i/o for message passing
  DO i=0,nprocs-1,dumpstride
    IF(myproc >= i.AND.myproc <= i+dumpstride-1)THEN

      CALL dtadump(nx,ny,nz,nzsoil,nstyps,                              &
               hdmpfmt,nchdmp,hdmpfn(1:ldmpf),                          &
               grdbas,filcmprs,                                         &
               u(1,1,1,tim),v(1,1,1,tim),                               &
               w(1,1,1,tim),ptprt(1,1,1,tim),                           &
               pprt(1,1,1,tim),qv(1,1,1,tim),                           &
               qscalar(:,:,:,tim,:),tke(1,1,1,tim),kmh,kmv,             &
               ubar,vbar,tem1,ptbar,pbar,rhobar,qvbar,                  &
               x,y,z,zp,zpsoil,                                         &
               soiltyp,stypfrct,vegtyp,lai,roufns,veg,                  &
               tsoil,qsoil,wetcanp,snowdpth,                            &
               raing,rainc,prcrate,                                     &
               radfrc,radsw,rnflx,radswnet,radlwin,                     &
               usflx,vsflx,ptsflx,qvsflx,                               &
               tem2,tem3,tem6)

    END IF
    IF (mp_opt > 0) CALL mpbarrier
  END DO

!
!-----------------------------------------------------------------------
!
!  Call EXBCDUMP to dump the external boundary fields.
!
!-----------------------------------------------------------------------
!
  IF ( extdadmp == 1 .AND. lbcopt == 2 ) THEN
    IF ( hdmpfmt == 5 .OR. hdmpfmt == 9 .OR. hdmpfmt == 11 ) THEN
      hdmpfmt1 = 1
    ELSE
      hdmpfmt1 = hdmpfmt
    END IF

    CALL gtdmpfn(runname(1:lfnkey),outdirname,oldirnam,curtim,1,        &
                 mgrid,nestgrd, exbcdmpfn, lexdmpf)

    exbcdmpfn(lexdmpf+1:lexdmpf+5) = '.exbc'

!  blocking inserted for ordering i/o for message passing
    DO i=0,nprocs-1,dumpstride
      IF(myproc >= i.AND.myproc <= i+dumpstride-1)THEN

        CALL exbcdump(nx,ny,nz,nzsoil,nstyps,                           &
                  hdmpfmt1,exbcdmpfn,grdbas,filcmprs,                   &
                  tem2,tem3,tem4,tem5,tem6,tem7,                        &
                  qscalar(:,:,:,tim,:),tke(1,1,1,tim),kmh,kmv,          &
                  ubar,vbar,tem1,ptbar,pbar,                            &
                  rhobar,qvbar, x,y,z,zp,zpsoil,                        &
                  soiltyp,stypfrct,vegtyp,lai,roufns,veg,               &
                  tsoil,qsoil,wetcanp,snowdpth,                         &
                  raing,rainc,prcrate,                                  &
                  radfrc,radsw,rnflx,radswnet,radlwin,                  &
                  usflx,vsflx,ptsflx,qvsflx,                            &
                  exbcbuf(nu0exb), exbcbuf(nv0exb),                     &
                  exbcbuf(nw0exb), exbcbuf(npt0exb),                    &
                  exbcbuf(npr0exb),exbcbuf(nqv0exb),                    &
                  exbcbuf(nqscalar0exb(1)),                             &
                  exbcbuf(nudtexb), exbcbuf(nvdtexb),                   &
                  exbcbuf(nwdtexb), exbcbuf(nptdtexb),                  &
                  exbcbuf(nprdtexb),exbcbuf(nqvdtexb),                  &
                  exbcbuf(nqscalardtexb(1)),                            &
                  tem8,tem9,tem10)


      END IF
      IF (mp_opt > 0) CALL mpbarrier
    END DO
  END IF

!-----------------------------------------------------------------------
!
!  Write out soil model variable file
!
!-----------------------------------------------------------------------

  IF (soildmp > 0) THEN

    CALL cvttsnd( curtim, timsnd, tmstrln )

    soiloutfl = runname(1:lfnkey)//".soilvar."//timsnd(1:tmstrln)

    temchar = soiloutfl
    soiloutfl = TRIM(outdirname)//'/'//temchar

    IF(myproc == 0) WRITE (6,*) 'Write soil initial data to ',TRIM(soiloutfl)

!    blocking inserted for ordering i/o for message passing
    DO i=0,nprocs-1,dumpstride
      IF(myproc >= i.AND.myproc <= i+dumpstride-1)THEN
        IF(mp_opt > 0 .AND. joindmp(FINDX_S) > 0) THEN
          CALL wrtjoinsoil(nx,ny,nzsoil,nstyps, soiloutfl,              &
                   dx,dy,zpsoil,                                        &
                   mapproj,trulat1,trulat2,trulon,sclfct,ctrlat,ctrlon, &
                   1,1,1,1,1,                                           &
                   tsoil,qsoil,wetcanp,snowdpth,soiltyp)
        ELSE
          CALL wrtsoil(nx,ny,nzsoil,nstyps, soiloutfl,                  &
                   dx,dy,zpsoil,                                        &
                   mapproj,trulat1,trulat2,trulon,sclfct,ctrlat,ctrlon, &
                   1,1,1,1,1,                                           &
                   tsoil,qsoil,wetcanp,snowdpth,soiltyp)
        END IF
      END IF
      IF (mp_opt > 0) CALL mpbarrier
    END DO

  END IF

!-----------------------------------------------------------------------
!
!  Write out terrain data
!
!-----------------------------------------------------------------------

  IF (terndmp > 0) THEN

    CALL getunit( nunit )

    ternfn = runname(1:lfnkey)//".trndata"
    lternfn = lfnkey + 8

    temchar = ternfn
    ternfn = TRIM(outdirname)//'/'//temchar
    lternfn  = LEN_TRIM(ternfn)


    IF (mp_opt > 0 .AND. joindmp(FINDX_T) == 0) THEN
      savename(1:256) = ternfn(1:256)
      CALL gtsplitfn(savename,1,1,loc_x,loc_y,1,1,0,0,0,lvldbg,ternfn,istatus)
      lternfn = LEN_TRIM(ternfn)
    END IF

    CALL fnversn(ternfn, lternfn )

    IF(myproc == 0) WRITE (6,*) 'Write terrain data to ',ternfn(1:lternfn)

!    blocking inserted for ordering i/o for message passing
    DO i=0,nprocs-1,dumpstride
      IF(myproc >= i.AND.myproc <= i+dumpstride-1)THEN

      IF(mp_opt >0 .AND. joindmp(FINDX_T) > 0) THEN
        CALL writjointrn(nx,ny,ternfn(1:lternfn), dx,dy,                &
            mapproj,trulat1,trulat2,trulon,sclfct,                      &
            ctrlat,ctrlon,hterain)
      ELSE
        CALL writtrn(nx,ny,ternfn(1:lternfn), dx,dy,                    &
            mapproj,trulat1,trulat2,trulon,sclfct,                      &
            ctrlat,ctrlon,hterain)
      END IF

      END IF
      IF (mp_opt > 0) CALL mpbarrier
    END DO
  END IF

!-----------------------------------------------------------------------
!
!  Write out surface property data file: sfcoutfl .
!
!-----------------------------------------------------------------------

  IF (sfcdmp > 0) THEN

    sfcoutfl = runname(1:lfnkey)//".sfcdata"
    lfn = lfnkey + 8

    temchar = sfcoutfl
    sfcoutfl = TRIM(outdirname)//'/'//temchar
    lfn  = LEN_TRIM(sfcoutfl)

    IF (mp_opt > 0 .AND. joindmp(FINDX_T) == 0) THEN
      savename(1:256) = sfcoutfl(1:256)
      CALL gtsplitfn(savename,1,1,loc_x,loc_y,1,1,0,0,0,lvldbg,sfcoutfl,istatus)
      lfn = LEN_TRIM(sfcoutfl)
    END IF

    CALL fnversn(sfcoutfl, lfn)

    IF(myproc == 0) WRITE (6,*) 'Write surface property data in ',sfcoutfl(1:lfn)

!    blocking inserted for ordering i/o for message passing
    DO i=0,nprocs-1,dumpstride
      IF(myproc >= i.AND.myproc <= i+dumpstride-1)THEN

        IF(mp_opt > 0 .AND. joindmp(FINDX_T) > 0) THEN
        CALL wrtjoinsfcdt(nx,ny,nstyps,sfcoutfl(1:lfn), dx,dy,          &
              mapproj,trulat1,trulat2,trulon,sclfct,ctrlat,ctrlon,      &
              1,1,1,1,1,0,                                              &
              soiltyp,stypfrct,vegtyp,lai,roufns,veg,veg)
        ELSE
        CALL wrtsfcdt(nx,ny,nstyps,sfcoutfl(1:lfn), dx,dy,              &
              mapproj,trulat1,trulat2,trulon,sclfct,ctrlat,ctrlon,      &
              1,1,1,1,1,0,                                              &
              soiltyp,stypfrct,vegtyp,lai,roufns,veg,veg)
        END IF

      END IF
      IF (mp_opt > 0) CALL mpbarrier
    END DO

  END IF       ! landin.eq.1

  800   CONTINUE                  ! Entry if hdmpfmt=0

  RETURN

END SUBROUTINE initout

!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE OUTPUT                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE output(mptr,nx,ny,nz,nzsoil,nstyps,exbcbufsz,                &
           u,v,w,wcont,ptprt,pprt,qv,qscalar,tke,                       &
           udteb, udtwb, vdtnb, vdtsb,                                  &
           pdteb ,pdtwb ,pdtnb ,pdtsb,                                  &
           ubar,vbar,ptbar,pbar,rhostr,qvbar,kmh,kmv,                   &
           x,y,z,zp,zpsoil,hterain, mapfct,j1,j2,j3,j3soil,j3inv,       &
           j3soilinv,soiltyp,stypfrct,vegtyp,lai,roufns,veg,            &
           tsoil,qsoil,wetcanp,snowdpth,qvsfc,                          &
           ptcumsrc,qcumsrc,w0avg,nca,kfraincv,                         &
           cldefi,xland,bmjraincv,                                      &
           raing,rainc,prcrate,exbcbuf,                                 &
           radfrc,radsw,rnflx,radswnet,radlwin,                         &
           usflx,vsflx,ptsflx,qvsflx,                                   &
           tem1,tem2,tem3,tem4,tem5,tem6,                               &
           tem7,tem8,tem9,tem10,tem11)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Coordinate the output of model data at a particular time.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  11/13/91.
!
!  MODIFICATION HISTORY:
!
!  6/06/92 (M. Xue)
!  Added full documentation.
!
!  7/13/92 (K. Brewster)
!  Added comment line variables to arg list to be passed to the
!  data dumping routines.
!
!  1/24/94 (Y. Liu)
!  Added surface variables to arg list to be passed to the
!  data dumping routines.
!
!  02/07/1995 (Yuhe Liu)
!  Added a new 2-D permanent array, veg(nx,ny), to the argument list
!
!  12/6/95 (J. Zong and M. Xue)
!  Added qs and qh to the argument list of REFLEC when it is called
!  to include the contributions of qs and qh to reflectivity.
!
!  2/2/96 (Donghai Wang & Yuhe Liu)
!  Added a 3-D array, mapfct, for map projection factor.
!
!  08/01/97 (Zonghui Huo)
!  Added Kain-fritsch cumulus parameterization scheme.
!
!  11/06/97 (D. Weber)
!  Added three additional levels to the mapfct array.  The three
!  levels (4,5,6) represent the inverse of the first three in order.
!  The inverse map factors are computed to improve efficiency.
!
!  4/15/1998 (Donghai Wang)
!  Added the source terms to the right hand terms of the qc,qr,qi,qs
!  equations due to the K-F cumulus parameterization.
!
!  4/15/1998 (Donghai Wang)
!  Added the running average vertical velocity (array w0avg)
!  for the K-F cumulus parameterization scheme.
!
!  12/09/1998 (Donghai Wang)
!  Added the snow cover.
!
!  03/13/2002 (Eric Kemp)
!  Added arrays for WRF BMJ cumulus scheme.
!
!  05/13/2002 (Jerry Brotzge)
!  Added arrays to allow for multiple soil schemes
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    mptr     Grid identifier.
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!    nzsoil   Number of grid points in the soil
!
!    u        x component of velocity (m/s)
!    v        y component of velocity (m/s)
!    w        Vertical component of Cartesian velocity (m/s)
!    wcont    Contravariant vertical velocity (m/s)
!    ptprt    Perturbation potential temperature (K)
!    pprt     Perturbation pressure (Pascal)
!    qv       Water vapor specific humidity (kg/kg)
!    qscalar
!    tke      Turbulent Kinetic Energy ((m/s)**2)
!
!    udteb    Time tendency of u field at east boundary (m/s**2)
!    udtwb    Time tendency of u field at west boundary (m/s**2)
!
!    vdtnb    Time tendency of v field at north boundary (m/s**2)
!    vdtsb    Time tendency of v field at south boundary (m/s**2)
!
!    pdteb    Time tendency of pprt field at east boundary (PASCAL/s)
!    pdtwb    Time tendency of pprt field at west boundary (PASCAL/s)
!    pdtnb    Time tendency of pprt field at north boundary (PASCAL/s)
!    pdtsb    Time tendency of pprt field at south boundary (PASCAL/s)
!
!    ubar     Base state zonal velocity component (m/s)
!    vbar     Base state meridional velocity component (m/s)
!    ptbar    Base state potential temperature (K)
!    pbar     Base state pressure (Pascal)
!    rhostr   Base state density * j3 (kg/m**3)
!    qvbar    Base state water vapor specific humidity (kg/kg)
!
!    kmh      Horizontal turb. mixing coef. for momentum ( m**2/s )
!    kmv      Vertical turb. mixing coef. for momentum ( m**2/s )
!
!    x        x coordinate of grid points in physical/comp. space (m)
!    y        y coordinate of grid points in physical/comp. space (m)
!    z        z coordinate of grid points in computational space (m)
!    zp       Vertical coordinate of grid points in physical space (m)
!    zpsoil   Vertical coordinate of grid points in the soil (m)
!    hterain  Terrain height (m)
!
!    mapfct   Map factors at scalar, u and v points
!
!    j1       Coordinate transformation Jacobian -d(zp)/dx
!    j2       Coordinate transformation Jacobian -d(zp)/dy
!    j3       Coordinate transformation Jacobian  d(zp)/dz
!    j3soil   Coordinate transformation Jacobian  d(zpsoil)/dz
!
!    soiltyp  Soil type at the horizontal grid points
!    stypfrct  Soil type fraction
!    vegtyp   Vegetation type at the horizontal grid points
!    lai      Leaf Area Index
!    roufns   Surface roughness
!    veg      Vegetation fraction
!
!    tsoil    Soil temperature (K)
!    qsoil    Soil moisture (m**3/m**3)
!    wetcanp  Canopy water amount
!    qvsfc    Effective qv at sfc.
!
!    ptcumsrc Source term in pt-equation due to cumulus parameterization
!    qcumsrc Source term in water equation due to cumulus parameterization
!
!    kfraincv K-F convective rainfall (cm)
!    nca      K-F counter for CAPE release
!    cldefi   BMJ cloud efficiency
!    xland    BMJ land/sea mask
!    bmjraincv BMJ convective rainfall (cm)
!
!    raing    Grid supersaturation rain
!    rainc    Cumulus convective rain
!    prcrate  Precipitation rates
!
!    radfrc   Radiation forcing (K/s)
!    radsw    Solar radiation reaching the surface
!    rnflx    Net radiation flux absorbed by surface
!    radswnet Net solar radiation
!    radlwin  Incoming longwave radiation
!
!    usflx    Surface flux of u-momentum (kg/(m*s**2))
!    vsflx    Surface flux of v-momentum (kg/(m*s**2))
!    ptsflx   Surface heat flux (K*kg/(m**2 * s ))
!    qvsflx   Surface moisture flux of (kg/(m**2 * s))
!
!  OUTPUT:
!
!    None.
!
!  WORK ARRAYS:
!
!    rhobar   Base state density (kg/m**3)
!    tem1     Temporary work array.
!    tem2     Temporary work array.
!    tem3     Temporary work array.
!    tem4     Temporary work array.
!    tem5     Temporary work array.
!    tem6     Temporary work array.
!    tem7     Temporary work array.
!    tem8     Temporary work array.
!    tem9     Temporary work array.
!    tem10    Temporary work array.
!    tem11    Temporary work array.
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
  INCLUDE 'phycst.inc'
  INCLUDE 'globcst.inc'
  INCLUDE 'grid.inc'          ! Grid & map parameters.
  INCLUDE 'bndry.inc'
  INCLUDE 'exbc.inc'
  INCLUDE 'mp.inc'            ! Message passing parameters.
  INCLUDE 'timelvls.inc'

!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------

  INTEGER :: mptr              ! Grid identifier.

  INTEGER :: nx,ny,nz          ! Number of grid points in 3 directions
  INTEGER :: nzsoil            ! Number of grid points in the soil

  REAL :: u     (nx,ny,nz,nt)  ! Total u-velocity (m/s)
  REAL :: v     (nx,ny,nz,nt)  ! Total v-velocity (m/s)
  REAL :: w     (nx,ny,nz,nt)  ! Total w-velocity (m/s)
  REAL :: wcont (nx,ny,nz)     ! Contravariant vertical velocity (m/s)
  REAL :: ptprt (nx,ny,nz,nt)  ! Perturbation potential temperature (K)
  REAL :: pprt  (nx,ny,nz,nt)  ! Perturbation pressure (Pascal)

  REAL :: qv    (nx,ny,nz,nt)  ! Water vapor specific humidity (kg/kg)

  REAL :: qscalar(nx,ny,nz,nt,nscalar)

  REAL :: tke   (nx,ny,nz,nt)  ! Turbulent Kinetic Energy ((m/s)**2)

  REAL :: udteb (ny,nz)        ! T-tendency of u at e-boundary (m/s**2)
  REAL :: udtwb (ny,nz)        ! T-tendency of u at w-boundary (m/s**2)

  REAL :: vdtnb (nx,nz)        ! T-tendency of v at n-boundary (m/s**2)
  REAL :: vdtsb (nx,nz)        ! T-tendency of v at s-boundary (m/s**2)

  REAL :: pdteb (ny,nz)        ! T-tendency of pprt at e-boundary (PASCAL/s)
  REAL :: pdtwb (ny,nz)        ! T-tendency of pprt at w-boundary (PASCAL/s)
  REAL :: pdtnb (nx,nz)        ! T-tendency of pprt at n-boundary (PASCAL/s)
  REAL :: pdtsb (nx,nz)        ! T-tendency of pprt at s-boundary (PASCAL/s)

  REAL :: ubar  (nx,ny,nz)     ! Base state u-velocity (m/s)
  REAL :: vbar  (nx,ny,nz)     ! Base state v-velocity (m/s)
  REAL :: ptbar (nx,ny,nz)     ! Base state potential temperature (K)
  REAL :: pbar  (nx,ny,nz)     ! Base state pressure (Pascal)
  REAL :: rhostr(nx,ny,nz)     ! Base state air density * j3 (kg/m**3)
  REAL :: qvbar (nx,ny,nz)     ! Base state water vapor specific humidity
                               ! (kg/kg)

  REAL :: kmh   (nx,ny,nz)     ! Horizontal turb. mixing coef. for
                               ! momentum. ( m**2/s )
  REAL :: kmv   (nx,ny,nz)     ! Vertical turb. mixing coef. for
                               ! momentum. ( m**2/s )

  REAL :: x     (nx)           ! The x-coord. of the physical and
                               ! computational grid. Defined at u-point.
  REAL :: y     (ny)           ! The y-coord. of the physical and
                               ! computational grid. Defined at v-point.
  REAL :: z     (nz)           ! The z-coord. of the computational grid.
                               ! Defined at w-point on the staggered grid.
  REAL :: zp    (nx,ny,nz)     ! The physical height coordinate defined at
                               ! w-point of the staggered grid.
  REAL :: zpsoil(nx,ny,nzsoil)       ! The physical height coordinate defined at
                               ! w-point of the soil.

  REAL :: hterain(nx,ny)       ! Terrain height.

  REAL :: mapfct(nx,ny,8)      ! Map factors at scalar, u and v points

  REAL :: j1    (nx,ny,nz)     ! Coordinate transformation Jacobian defined
                               ! as - d( zp )/d( x )
  REAL :: j2    (nx,ny,nz)     ! Coordinate transformation Jacobian defined
                               ! as - d( zp )/d( y )
  REAL :: j3    (nx,ny,nz)     ! Coordinate transformation Jacobian defined
                               ! as d( zp )/d( z )

  REAL :: j3soil (nx,ny,nzsoil)   ! Coordinate transformation Jacobian defined
                               ! as d( zpsoil )/d( z )

  REAL :: j3inv (nx,ny,nz)     ! Coordinate transformation Jacobian defined
                               ! as d( zp )/d( z )

  REAL :: j3soilinv (nx,ny,nzsoil) ! Coordinate transformation Jacobian defined
                                   ! as d( zpsoil )/d( z )

  REAL :: rhobar(nx,ny,nz)     ! Base state air density (kg/m**3)

  INTEGER :: nstyps                    ! Number of soil types
  INTEGER :: soiltyp (nx,ny,nstyps)    ! Soil type at each point
  REAL :: stypfrct(nx,ny,nstyps)    ! Soil type fraction
  INTEGER :: vegtyp (nx,ny)            ! Vegetation type at each point
  REAL :: lai    (nx,ny)            ! Leaf Area Index
  REAL :: roufns (nx,ny)            ! Surface roughness
  REAL :: veg    (nx,ny)            ! Vegetation fraction

  REAL :: qvsfc  (nx,ny,0:nstyps)      ! Effective qv at sfc.
  REAL :: tsoil  (nx,ny,nzsoil,0:nstyps) ! Deep soil temperature (K)
  REAL :: qsoil  (nx,ny,nzsoil,0:nstyps) ! Deep soil moisture (m**3/m**3)
  REAL :: wetcanp(nx,ny,0:nstyps)      ! Canopy water amount
  REAL :: snowdpth(nx,ny)            ! Snow depth (m)

  REAL :: ptcumsrc(nx,ny,nz)   ! Source term in pt-equation due
                               ! to cumulus parameterization
  REAL :: qcumsrc(nx,ny,nz,5)  ! Source term in water equations due
                               ! to cumulus parameterization:
                               ! qcumsrc(1,1,1,1) for qv equation
                               ! qcumsrc(1,1,1,2) for qc equation
                               ! qcumsrc(1,1,1,3) for qr equation
                               ! qcumsrc(1,1,1,4) for qi equation
                               ! qcumsrc(1,1,1,5) for qs equation
  REAL :: w0avg(nx,ny,nz)      ! a closing running average vertical
                               ! velocity in 10min for K-F scheme
  REAL :: kfraincv(nx,ny)      ! K-F convective rainfall (cm)
  INTEGER :: nca(nx,ny)        ! K-F counter for CAPE release

!EMK BMJ
  REAL,INTENT(IN) :: cldefi(nx,ny) ! BMJ cloud efficiency
  REAL,INTENT(IN) :: xland(nx,ny)  ! BMJ land mask
                                   ! (1.0 = land, 2.0 = sea)
  REAL,INTENT(IN) :: bmjraincv(nx,ny) ! BMJ convective rainfall (cm)
!EMK END

  REAL :: raing(nx,ny)         ! Grid supersaturation rain
  REAL :: rainc(nx,ny)         ! Cumulus convective rain
  REAL :: prcrate(nx,ny,4)     ! precipitation rate (kg/(m**2*s))
                               ! prcrate(1,1,1) = total precipitation rate
                               ! prcrate(1,1,2) = grid scale precip. rate
                               ! prcrate(1,1,3) = cumulus precip. rate
                               ! prcrate(1,1,4) = microphysics precip. rate

  REAL :: radfrc(nx,ny,nz)     ! Radiation forcing (K/s)
  REAL :: radsw (nx,ny)        ! Solar radiation reacing the surface
  REAL :: rnflx (nx,ny)        ! Net absorbed radiation by the surface
  REAL :: radswnet (nx,ny)     ! Net solar radiation, SWin - SWout
  REAL :: radlwin  (nx,ny)     ! Incoming longwave radiation

  REAL :: usflx (nx,ny)        ! Surface flux of u-momentum (kg/(m*s**2))
  REAL :: vsflx (nx,ny)        ! Surface flux of v-momentum (kg/(m*s**2))
  REAL :: ptsflx(nx,ny)        ! Surface heat flux (K*kg/(m*s**2))
  REAL :: qvsflx(nx,ny)        ! Surface moisture flux (kg/(m**2*s))

  INTEGER :: exbcbufsz         ! Size of external buffer array
  REAL :: exbcbuf ( exbcbufsz) ! External buffer array

  REAL :: tem1  (nx,ny,nz)     ! Temporary work array
  REAL :: tem2  (nx,ny,nz)     ! Temporary work array
  REAL :: tem3  (nx,ny,nz)     ! Temporary work array
  REAL :: tem4  (nx,ny,nz)     ! Temporary work array
  REAL :: tem5  (nx,ny,nz)     ! Temporary work array
  REAL :: tem6  (nx,ny,nz)     ! Temporary work array
  REAL :: tem7  (nx,ny,nz)     ! Temporary work array
  REAL :: tem8  (nx,ny,nz)     ! Temporary work array
  REAL :: tem9  (nx,ny,nz)     ! Temporary work array
  REAL :: tem10 (nx,ny,nz)     ! Temporary work array
  REAL :: tem11 (nx,ny,nz)     ! Temporary work array

!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: tlevel, tim       ! Time level at which data are printed.
  CHARACTER (LEN=256) :: exbcdmpfn
  INTEGER :: lexdmpf
  INTEGER :: grdbas
  REAL    :: timage
  INTEGER :: i,j,k,istatwrt

  LOGICAL :: dumphist

  CHARACTER (LEN=256) :: soiloutfl,temchar
  CHARACTER (LEN=80 ) :: timsnd
  INTEGER :: lfn,tmstrln

  CHARACTER (LEN=256) :: savename, outdirname
  INTEGER             :: oldirnam
  INTEGER :: istatus
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  DO k=1,nz-1
    DO j=1,ny-1
      DO i=1,nx-1
        rhobar(i,j,k) = rhostr(i,j,k)*j3inv(i,j,k)
      END DO
    END DO
  END DO

  mgrid = mptr
!
!-----------------------------------------------------------------------
!
!  Printing, diagnostics and history dump of data at time tpresent
!
!-----------------------------------------------------------------------
!
  tlevel = tpresent

!
!-----------------------------------------------------------------------
!
!  Dump restart data files every nrstout number of time steps.
!
!-----------------------------------------------------------------------
!
  IF( nrstout > 0 .AND. MOD(nstep,nrstout) == 0) THEN

!    blocking inserted for ordering i/o for message passing
    DO i=0,nprocs-1,dumpstride
      IF(myproc >= i .AND. myproc <= i+dumpstride-1)THEN

        IF (mp_opt > 0 .AND. joindmp(FINDX_R) > 0) THEN

        CALL rstjoinout(nx,ny,nz,nzsoil,nstyps, exbcbufsz,              &
                u,v,w,ptprt,pprt,qv,qscalar,tke,                        &
                udteb, udtwb, vdtnb, vdtsb,                             &
                pdteb ,pdtwb ,pdtnb ,pdtsb,                             &
                ubar,vbar,ptbar,pbar,rhostr,qvbar,                      &
                x,y,z,zp,zpsoil,hterain,mapfct,                         &
                soiltyp,stypfrct,vegtyp,lai,roufns,veg,                 &
                tsoil,qsoil,wetcanp,snowdpth, qvsfc,                    &
                ptcumsrc,qcumsrc,w0avg,nca,kfraincv,                    &
                cldefi,xland,bmjraincv,                                 &
                radfrc,radsw,rnflx,radswnet,radlwin,                    &
                raing,rainc,prcrate,exbcbuf, tem1)
        ELSE

        CALL rstout(nx,ny,nz,nzsoil,nstyps, exbcbufsz,                  &
                u,v,w,ptprt,pprt,qv,qscalar,tke,                        &
                udteb, udtwb, vdtnb, vdtsb,                             &
                pdteb ,pdtwb ,pdtnb ,pdtsb,                             &
                ubar,vbar,ptbar,pbar,rhostr,qvbar,                      &
                x,y,z,zp,zpsoil,hterain,mapfct,                         &
                soiltyp,stypfrct,vegtyp,lai,roufns,veg,                 &
                tsoil,qsoil,wetcanp,snowdpth, qvsfc,                    &
                ptcumsrc,qcumsrc,w0avg,nca,kfraincv,                    &
                cldefi,xland,bmjraincv,                                 &
                radfrc,radsw,rnflx,radswnet,radlwin,                    &
                raing,rainc,prcrate,exbcbuf, tem1)

        END IF  ! mp_opt >0 AND joindmp >0

      END IF
      IF (mp_opt > 0) CALL mpbarrier
    END DO

  END IF
!
!-----------------------------------------------------------------------
!
!  Calculate the max and min values of variables and write them
!  into a separate file.
!
!-----------------------------------------------------------------------
!
  IF( nmaxmin > 0) THEN
    IF(MOD(nstep,nmaxmin) == 0) THEN

      CALL maxmin(mptr,nx,ny,nz,nzsoil,tlevel,rhobar,                   &
                  u,v,w,wcont,ptprt,pprt,qv,qscalar,tke,kmh,kmv,        &
                  x,y,z,zp,zpsoil,mapfct,                               &
                  tsoil(1,1,1,0),qsoil(1,1,1,0),wetcanp(1,1,0),         &
                  tem1,tem2,tem3,tem4)

    END IF
  END IF
!
!-----------------------------------------------------------------------
!
!  Calculation of energy statistics and printing
!
!-----------------------------------------------------------------------
!
  IF( nenergy > 0) THEN
    IF(MOD(nstep,nenergy) == 0) THEN

      CALL energy(nx,ny,nz,                                             &
                  u,v,w,ptprt,pprt,qv,qscalar,rhobar,                   &
                  x,y,z,zp,hterain, j3,                                 &
                  tem1,tem2,tem3,tem4,tem5)

    END IF
  END IF
!
!-----------------------------------------------------------------------
!
!  Print out formatted 2-d slices of 3-d arrays
!
!-----------------------------------------------------------------------
!
  IF( nfmtprt > 0) THEN
    IF(MOD(nstep,nfmtprt) == 0) THEN

      WRITE(6,'(1x,a,f13.3,a,i10)')                                     &
           'Print fields at time=',curtim,'(s) with nstep=',nstep

      CALL fmtprt(nx,ny,nz,tlevel,                                      &
                  u,v,w,ptprt,pprt,qv,qscalar,tke,kmh,kmv,              &
                  x,y,z,zp,hterain, j1,j2,j3, tem1)

    END IF
  END IF

!
!-----------------------------------------------------------------------
!
!  Produce history data dump every nhisdmp number of time steps.
!
!-----------------------------------------------------------------------
!
  dumphist = .false.

  IF ( nhisdmp > 0 .AND. hdmpfmt /= 0) THEN
    IF ( hdmpopt == 1 .AND. nstep >= nstrtdmp ) THEN
      dumphist = MOD(nstep-nstrtdmp,nhisdmp) == 0
    ELSE IF ( hdmpopt == 2 .AND. nhisdmp <= numhdmp ) THEN
      DO WHILE ( nstep > hdmpstp(nhisdmp) .AND. nhisdmp < numhdmp )
         nhisdmp = nhisdmp + 1
      END DO
      dumphist = ( nstep == hdmpstp(nhisdmp) )
      IF ( dumphist ) nhisdmp = nhisdmp + 1
    END IF
  END IF

  IF ( dumphist ) THEN

    CALL get_output_dirname(1,dirname,curtim,1,outdirname,istatus)
    oldirnam = LEN_TRIM(outdirname)
!
!-----------------------------------------------------------------------
!
!  Find a unique name hdmpfn(1:ldmpf) for the history dump data
!  set at time 'curtim'.
!
!  For the savi3D data dump case, the file name is specified only
!  once in INITOUT.
!
!-----------------------------------------------------------------------
!
    IF( hdmpfmt /= 5 .AND. hdmpfmt /= 9 ) THEN

      CALL gtdmpfn(runname(1:lfnkey),outdirname,                        &
                   oldirnam,curtim,hdmpfmt,                             &
                   mgrid,nestgrd, hdmpfn, ldmpf)

    END IF
    IF (myproc == 0) WRITE(6,'(1x,a,a)')  &
         'History data dump in file ',hdmpfn(1:ldmpf)

    grdbas = 0      ! No base state or grid array is dumped.

    tim = tpresent

!    blocking inserted for ordering i/o for message passing
    DO i=0,nprocs-1,dumpstride
      IF(myproc >= i .AND. myproc <= i+dumpstride-1)THEN

        CALL dtadump(nx,ny,nz,nzsoil,nstyps,                            &
                 hdmpfmt,nchdmp,hdmpfn(1:ldmpf),                        &
                 grdbas,filcmprs,                                       &
                 u(1,1,1,tim),v(1,1,1,tim),                             &
                 w(1,1,1,tim),                                          &
                 ptprt(1,1,1,tim),                                      &
                 pprt(1,1,1,tim),qv(1,1,1,tim),                         &
                 qscalar(:,:,:,tim,:),tke(1,1,1,tim),kmh,kmv,           &
                 ubar,vbar,tem1,ptbar,pbar,rhobar,qvbar,                &
                 x,y,z,zp,zpsoil,                                       &
                 soiltyp,stypfrct,vegtyp,lai,roufns,veg,                &
                 tsoil,qsoil,wetcanp,snowdpth,                          &
                 raing,rainc,prcrate,                                   &
                 radfrc,radsw,rnflx,radswnet,radlwin,                   &
                 usflx,vsflx,ptsflx,qvsflx,                             &
                 tem2,tem3,tem6)

      END IF
      IF (mp_opt > 0) CALL mpbarrier
    END DO

!
!-----------------------------------------------------------------------
!
!  Dump the fields of external boundary conditions for diagnostic
!  purpose.
!
!-----------------------------------------------------------------------
!
    IF ( extdadmp == 1 .AND. lbcopt == 2 ) THEN

!
!-----------------------------------------------------------------------
!
!  Call EXBCDUMP to dump the external boundary fields.
!
!-----------------------------------------------------------------------
!
      CALL gtdmpfn(runname(1:lfnkey),outdirname,oldirnam,curtim,1,      &
                   mgrid,nestgrd, exbcdmpfn, lexdmpf)

      exbcdmpfn(lexdmpf+1:lexdmpf+5) = '.exbc'

!    blocking inserted for ordering i/o for message passing
      DO i=0,nprocs-1,dumpstride
        IF(myproc >= i.AND.myproc <= i+dumpstride-1)THEN

          CALL exbcdump( nx,ny,nz,nzsoil,nstyps,                        &
                    1,exbcdmpfn,grdbas,filcmprs,                        &
                    tem2,tem3,tem4,tem5,tem6,tem7,                      &
                    qscalar(:,:,:,tim,:),tke(:,:,:,tim),kmh,kmv,        &
                    ubar,vbar,tem1,ptbar,pbar,                          &
                    rhobar,qvbar,                                       &
                    x,y,z,zp,zpsoil,                                    &
                    soiltyp,stypfrct,vegtyp,lai,roufns,veg,             &
                    tsoil,qsoil,wetcanp,snowdpth,                       &
                    raing,rainc,prcrate,                                &
                    radfrc,radsw,rnflx,radswnet,radlwin,                &
                    usflx,vsflx,ptsflx,qvsflx,                          &
                    exbcbuf(nu0exb), exbcbuf(nv0exb),                   &
                    exbcbuf(nw0exb), exbcbuf(npt0exb),                  &
                    exbcbuf(npr0exb),exbcbuf(nqv0exb),                  &
                    exbcbuf(nqscalar0exb(1)),                           &
                    exbcbuf(nudtexb), exbcbuf(nvdtexb),                 &
                    exbcbuf(nwdtexb), exbcbuf(nptdtexb),                &
                    exbcbuf(nprdtexb),exbcbuf(nqvdtexb),                &
                    exbcbuf(nqscalardtexb(1)),                          &
                    tem8,tem9,tem10)

        END IF
        IF (mp_opt > 0) CALL mpbarrier
      END DO
    END IF


!-----------------------------------------------------------------------
!
!  Write out soil model variable file
!
!-----------------------------------------------------------------------

    IF (soildmp > 0) THEN

      CALL cvttsnd( curtim, timsnd, tmstrln )

      soiloutfl = runname(1:lfnkey)//".soilvar."//timsnd(1:tmstrln)
      lfn = lfnkey + 9 + tmstrln

      temchar = soiloutfl
      soiloutfl = TRIM(outdirname)//'/'//temchar
      lfn  = LEN_TRIM(soiloutfl)

      !IF (mp_opt > 0 .AND. joindmp == 0) THEN
      !  savename(1:256) = soiloutfl(1:256)
      !  CALL gtsplitfn(savename,1,1,loc_x,loc_y,1,1,0,0,0,lvldbg,soiloutfl,istatwrt)
      !  lfn = LEN_TRIM(soiloutfl)
      !END IF

      !CALL fnversn(soiloutfl, lfn)

      !IF(myproc == 0) WRITE (6,*) 'Write soil initial data to ',soiloutfl(1:lfn)

      ! blocking inserted for ordering i/o for message passing
      DO i=0,nprocs-1,dumpstride
        IF(myproc >= i.AND.myproc <= i+dumpstride-1)THEN

          IF(mp_opt > 0 .AND. joindmp(FINDX_S) > 0) THEN
            CALL wrtjoinsoil(nx,ny,nzsoil,nstyps, soiloutfl(1:lfn),     &
                     dx,dy,zpsoil,                                      &
                     mapproj,trulat1,trulat2,trulon,sclfct,ctrlat,ctrlon, &
                     1,1,1,1,1,                                         &
                     tsoil,qsoil,wetcanp,snowdpth,soiltyp)
          ELSE
            CALL wrtsoil(nx,ny,nzsoil,nstyps, soiloutfl(1:lfn),         &
                     dx,dy,zpsoil,                                      &
                     mapproj,trulat1,trulat2,trulon,sclfct,ctrlat,ctrlon, &
                     1,1,1,1,1,                                         &
                     tsoil,qsoil,wetcanp,snowdpth,soiltyp)
          END IF

        END IF
        IF (mp_opt > 0) CALL mpbarrier
      END DO

    END IF

  END IF ! (dumphist)

!-----------------------------------------------------------------------
!
!  Produces HDF images of Radar reflectivity factor at z=.25 and 4 km.
!
!-----------------------------------------------------------------------

  IF( nimgdmp /= 0 ) THEN
    IF( imgopt==1 .AND. (MOD(nstep,nimgdmp)==0 .OR. nstep==1)) THEN

      tim = tpresent
      timage = curtim
      IF(nstep == 1) THEN
        timage=0.
        tim=tpast
      END IF

      IF (P_QR > 0) THEN
        tem2(:,:,:) = qscalar(:,:,:,tim,P_QR)
      ELSE
        tem2(:,:,:) = 0.0
      END IF

      IF (P_QS > 0) THEN
        tem3(:,:,:) = qscalar(:,:,:,tim,P_QS)
      ELSE
        tem3(:,:,:) = 0.0
      END IF

      IF (P_QH > 0) THEN
        tem4(:,:,:) = qscalar(:,:,:,tim,P_QH)
      ELSE
        tem4(:,:,:) = 0.0
      END IF

      CALL reflec(nx,ny,nz, rhobar, tem2, tem3, tem4, tem1)

      CALL wrtvar2(nx,ny,1,tem1(1,1,7),                                &
              'k7refl','Reflectivity k=7','dBZ',                       &
              timage,runname,outdirname,hdmpfmt,hdfcompr,joindmp(FINDX_A),istatwrt)
      tem1=0.
      CALL avgx(u(1,1,1,tim), 0, nx,ny,nz,1,nx-1,1,ny-1,1,nz-1,tem1)

      CALL wrtvar2(nx,ny,1,tem1(1,1,7),                                &
              'k7uwnd','U wind comp k=7','m/s',                        &
              timage,runname,outdirname,hdmpfmt,hdfcompr,joindmp(FINDX_A),istatwrt)
      tem2=0.
      CALL avgy(v(1,1,1,tim), 0, nx,ny,nz,1,nx-1,1,ny-1,1,nz-1,tem2)

      CALL wrtvar2(nx,ny,1,tem2(1,1,7),                                &
              'k7vwnd','V wind component','m/s',                       &
              timage,runname,outdirname,hdmpfmt,hdfcompr,joindmp(FINDX_A),istatwrt)
!
!-----------------------------------------------------------------------
!
!  Compute the vertical vorticity by doing 2*dx derivatives across
!  a scalar point.
!
!-----------------------------------------------------------------------
!
      DO k=2,nz-1
        DO j=2,ny-2
          DO i=2,nx-2
            tem3(i,j,k)=0.5E05 *                                      &
               (((tem2(i+1,j,k)-tem2(i-1,j,k))*mapfct(i,j,1)*dxinv) - &
                ((tem1(i,j+1,k)-tem1(i,j-1,k))*mapfct(i,j,1)*dyinv))
          END DO
        END DO
      END DO
!
!-----------------------------------------------------------------------
!
!  Fill edges and corners
!
!-----------------------------------------------------------------------
!
      DO k=2,nz-1
        DO j=2,ny-2
          tem3( 1,j,k)=tem3(   2,j,k)
          tem3(nx-1,j,k)=tem3(nx-2,j,k)
          tem3(nx,j,k)=tem3(nx-2,j,k)
        END DO
      END DO

      DO k=2,nz-1
        DO i=1,nx
          tem3(i, 1,k)=tem3(i,   2,k)
          tem3(i,ny-1,k)=tem3(i,ny-2,k)
          tem3(i,ny,k)=tem3(i,ny-2,k)
        END DO
      END DO

      DO j=1,ny
        DO i=1,nx
          tem3(i,j, 1)=tem3(i,j,   2)
          tem3(i,j,nz)=tem3(i,j,nz-1)
        END DO
      END DO

      CALL wrtvar2(nx,ny,1,tem3(1,1,7),                                &
           'k7vort','Vertical Vort k=7','10^5 1/s',                    &
           timage,runname,outdirname,hdmpfmt,hdfcompr,joindmp(FINDX_A),istatwrt)

!    CALL img3d0(tem1,nx,1,nx-1,-2, ny,1,ny-1,-2, nz,1,nz-1,-2,
!    :              75.0 , 0.0 , 'rf' ,runname(1:lfnkey),curtim, 1, 2 ,
!    :              1 ,1 ,
!    :              mptr, nestgrd, tem1)

!    CALL avgz(tem1 , 1 ,
!    :       nx,ny,nz, 1,nx-1, 1,ny-1, 2,nz-1, tem2)

!      CALL img3d0(tem1,nx,1,nx-1,-2, ny,1,ny-1,-2, nz,1,nz-1,-2,        &
!                  75.0 , 0.0 , 'rf' ,runname(1:lfnkey),curtim, 1, 24 ,  &
!                  1 ,1 , mptr, nestgrd, tem1)

!      CALL img3d0(ptprt(1,1,1,tim)                                      &
!                  ,nx,1,nx-1,-2, ny,1,ny-1,-2, nz,1,nz-1,-2,            &
!                  0.0 , -12.0 , 'pt' ,runname(1:lfnkey),curtim, 1, 2 ,  &
!                  1 ,1 , mptr, nestgrd, tem1)

!      CALL img3d0(w(1,1,1,tim)                                          &
!                  ,nx,1,nx-1,-2, ny,1,ny-1,-2, nz,1,nz-1,-2,            &
!                  10.0 , -10.0 , 'w' ,runname(1:lfnkey),curtim, 1, 5 ,  &
!                  1 ,1 , mptr, nestgrd, tem1)

!      CALL img3d0(w(1,1,1,tim)                                          &
!                  ,nx,1,nx-1,-2, ny,1,ny-1,-2, nz,1,nz-1,-2,            &
!                  60.0 , -30.0 , 'w' ,runname(1:lfnkey),curtim, 1,24 ,  &
!                  1 ,1 , mptr, nestgrd, tem1)


    END IF
  END IF

  IF( nplots /= 0) THEN

    IF( MOD(nstep,nplots) == 0 .AND. pltopt == 1) THEN

      tim = tpresent

      CALL wrtxyslic(nx,ny,nz, u(1,1,1,tim),v(1,1,1,tim),w(1,1,1,tim),  &
                     ptprt(1,1,1,tim),pprt(1,1,1,tim),qv(1,1,1,tim),    &
                     qscalar(:,:,:,tim,:),ubar,vbar,ptbar,pbar,qvbar,   &
                     x,y,zp,   50.0,runname(1:lfnkey),curtim,           &
                     tem1,tem2)

      CALL wrtxyslic(nx,ny,nz, u(1,1,1,tim),v(1,1,1,tim),w(1,1,1,tim),  &
                     ptprt(1,1,1,tim),pprt(1,1,1,tim),qv(1,1,1,tim),    &
                     qscalar(:,:,:,tim,:),ubar,vbar,ptbar,pbar,qvbar,   &
                     x,y,zp, 1000.0,runname(1:lfnkey),curtim,           &
                     tem1,tem2)

      CALL wrtxyslic(nx,ny,nz, u(1,1,1,tim),v(1,1,1,tim),w(1,1,1,tim),  &
                     ptprt(1,1,1,tim),pprt(1,1,1,tim),qv(1,1,1,tim),    &
                     qscalar(:,:,:,tim,:),ubar,vbar,ptbar,pbar,qvbar,   &
                     x,y,zp, 2000.0,runname(1:lfnkey),curtim,           &
                     tem1,tem2)

      CALL wrtxyslic(nx,ny,nz, u(1,1,1,tim),v(1,1,1,tim),w(1,1,1,tim),  &
                     ptprt(1,1,1,tim),pprt(1,1,1,tim),qv(1,1,1,tim),    &
                     qscalar(:,:,:,tim,:),ubar,vbar,ptbar,pbar,qvbar,   &
                     x,y,zp, 4000.0,runname(1:lfnkey),curtim,           &
                     tem1,tem2)

      CALL wrtxyslic(nx,ny,nz, u(1,1,1,tim),v(1,1,1,tim),w(1,1,1,tim),  &
                     ptprt(1,1,1,tim),pprt(1,1,1,tim),qv(1,1,1,tim),    &
                     qscalar(:,:,:,tim,:),ubar,vbar,ptbar,pbar,qvbar,   &
                     x,y,zp, 8000.0,runname(1:lfnkey),curtim,           &
                     tem1,tem2)

    END IF
  END IF

  RETURN
END SUBROUTINE output

!######################################################################

SUBROUTINE dfilter(mptr,nx,ny,nz,nzsoil,nstyps,                         &
                   u,v,w,wcont,ptprt,pprt,qv,qscalar,tke,               &
                   tsoil,qsoil,wetcanp,snowdpth,                        &
                   numdfvar_atm3d,numdfvar_soil,                        &
                   numdfvar_soil2d,numdfvar_soil2ds,                    &
                   idfstp,ndfstp,rwght, addwght,                        &
                   df_store_atm3d,df_store_soil,                        &
                   df_store_soil2d,df_store_soil2ds)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!    Digital filtering atmospheric fields
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ting Lei
!  03/23/2011.
!  Modified from output for doing digtal filtering
!
!  MODIFICATION HISTORY:
!
!  3/23/11 (Y. Wang)
!  Added full documentation and modified based on ARPS coding conventions.
!  Adopted for arps5.3.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!  OUTPUT:
!
!    None.
!
!  WORK ARRAYS:
!
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'phycst.inc'
  INCLUDE 'globcst.inc'
  INCLUDE 'grid.inc'          ! Grid & map parameters.
  INCLUDE 'bndry.inc'
  INCLUDE 'mp.inc'            ! Message passing parameters.

  INCLUDE 'timelvls.inc'

  INTEGER :: mptr              ! Grid identifier.
  INTEGER :: nx,ny,nz          ! Number of grid points in 3 directions
  INTEGER :: nzsoil            ! Number of grid points in the soil

  REAL :: u     (nx,ny,nz,nt)  ! Total u-velocity (m/s)
  REAL :: v     (nx,ny,nz,nt)  ! Total v-velocity (m/s)
  REAL :: w     (nx,ny,nz,nt)  ! Total w-velocity (m/s)
  REAL :: wcont (nx,ny,nz)     ! Contravariant vertical velocity (m/s)
  REAL :: ptprt (nx,ny,nz,nt)  ! Perturbation potential temperature (K)
  REAL :: pprt  (nx,ny,nz,nt)  ! Perturbation pressure (Pascal)

  REAL :: qv    (nx,ny,nz,nt)  ! Water vapor specific humidity (kg/kg)
  REAL :: qscalar(nx,ny,nz,nt,nscalar)
!  REAL :: qc    (nx,ny,nz,nt)  ! Cloud water mixing ratio (kg/kg)
!  REAL :: qr    (nx,ny,nz,nt)  ! Rain water mixing ratio (kg/kg)
!  REAL :: qi    (nx,ny,nz,nt)  ! Cloud ice mixing ratio (kg/kg)
!  REAL :: qs    (nx,ny,nz,nt)  ! Snow mixing ratio (kg/kg)
!  REAL :: qh    (nx,ny,nz,nt)  ! Hail mixing ratio (kg/kg)
  REAL :: tke   (nx,ny,nz,nt)  ! Turbulent Kinetic Energy ((m/s)**2)



  INTEGER :: nstyps                    ! Number of soil types
  REAL :: veg    (nx,ny)            ! Vegetation fraction

  REAL :: tsoil  (nx,ny,nzsoil,0:nstyps) ! Deep soil temperature (K)
  REAL :: qsoil  (nx,ny,nzsoil,0:nstyps) ! Deep soil moisture (m**3/m**3)
  REAL :: wetcanp(nx,ny,0:nstyps)      ! Canopy water amount
  REAL :: snowdpth(nx,ny)            ! Snow depth (m)
  INTEGER :: numdfvar_atm3d,numdfvar_soil,numdfvar_soil2d,numdfvar_soil2ds
  REAL :: df_store_atm3d(nx,ny,nz,numdfvar_atm3d)
  REAL :: df_store_soil(nx,ny,nzsoil,0:nstyps,numdfvar_soil)
  REAL :: df_store_soil2d(nx,ny,0:nstyps,numdfvar_soil2d)
  REAL :: df_store_soil2ds(nx,ny,numdfvar_soil2ds)
  INTEGER :: idfstp,ndfstp
  REAL :: addwght
  REAL :: rwght

!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: nq

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  addwght=addwght+rwght
  WRITE(*,*)'rwght,addwght= ',rwght,addwght
  CALL store_atmvar(nx,ny,nz,u(:,:,:,3),df_store_atm3d(:,:,:,1),rwght)
  CALL store_atmvar(nx,ny,nz,v(:,:,:,3),df_store_atm3d(:,:,:,2),rwght)
  CALL store_atmvar(nx,ny,nz,w(:,:,:,3),df_store_atm3d(:,:,:,3),rwght)
  CALL store_atmvar(nx,ny,nz,wcont(:,:,:),df_store_atm3d(:,:,:,4),rwght)
  CALL store_atmvar(nx,ny,nz,ptprt(:,:,:,3),df_store_atm3d(:,:,:,5),rwght)
  CALL store_atmvar(nx,ny,nz,pprt(:,:,:,3),df_store_atm3d(:,:,:,6),rwght)
  CALL store_atmvar(nx,ny,nz,tke(:,:,:,3), df_store_atm3d(:,:,:,7),rwght)
  CALL store_atmvar(nx,ny,nz,qv(:,:,:,3),  df_store_atm3d(:,:,:,8),rwght)
  DO nq = 1, nscalarq
    CALL store_atmvar(nx,ny,nz,qscalar(:,:,:,3,nq),df_store_atm3d(:,:,:,8+nq),rwght)
  END DO

  CALL store_soilvar(nx,ny,nzsoil,nstyps,tsoil(:,:,:,:),df_store_soil(:,:,:,:,1),rwght)
  CALL store_soilvar(nx,ny,nzsoil,nstyps,qsoil(:,:,:,:),df_store_soil(:,:,:,:,2),rwght)
  CALL store_soil2dvar(nx,ny,nstyps,wetcanp(:,:,:),df_store_soil2d(:,:,:,1),rwght)
  CALL store_soil2dvars(nx,ny,snowdpth(:,:),df_store_soil2ds(:,:,1),rwght)

  IF (idfstp == ndfstp)  THEN ! the last step of df
    IF(ABS(addwght-1) > 1.0E-4) THEN
      WRITE(*,*)'wrong, in dflter, the addwght doesn"t eq. 1 ,stop'
      CALL arpsstop('ERROR: addwght not equal 1.',1)
    ENDIF
    CALL re_store_atmvar(nx,ny,nz,nt,u,    df_store_atm3d(:,:,:,1))
    CALL re_store_atmvar(nx,ny,nz,nt,v,    df_store_atm3d(:,:,:,2))
    CALL re_store_atmvar(nx,ny,nz,nt,w,    df_store_atm3d(:,:,:,3))
    CALL re_store_atmvar(nx,ny,nz,nt,wcont,df_store_atm3d(:,:,:,4))
    CALL re_store_atmvar(nx,ny,nz,nt,ptprt,df_store_atm3d(:,:,:,5))
    CALL re_store_atmvar(nx,ny,nz,nt,pprt, df_store_atm3d(:,:,:,6))
    CALL re_store_atmvar(nx,ny,nz,nt,tke,  df_store_atm3d(:,:,:,7))
    CALL re_store_atmvar(nx,ny,nz,nt,qv,   df_store_atm3d(:,:,:,8))
    DO nq = 1, nscalarq
      CALL re_store_atmvar(nx,ny,nz,nt,qscalar(:,:,:,1,nq),df_store_atm3d(:,:,:,8+nq))
    END DO
    CALL re_store_soilvar   (nx,ny,nzsoil,nstyps,tsoil,df_store_soil(:,:,:,:,1))
    CALL re_store_soilvar   (nx,ny,nzsoil,nstyps,qsoil,df_store_soil(:,:,:,:,2))
    CALL re_store_soil2dvar (nx,ny,nstyps,wetcanp,     df_store_soil2d(:,:,:,1))
    CALL re_store_soil2dvars(nx,ny,snowdpth,           df_store_soil2ds(:,:,1))
  ENDIF

  RETURN
END SUBROUTINE dfilter

subroutine store_atmvar(nx,ny,nz,ain,astore,rwght)
  implicit none
  integer nx,ny,nz
  integer i,j,k
  real rwght
  real ain(nx,ny,nz),astore(nx,ny,nz)
!clt
  do k=1,nz
   do j=1,ny
    do i=1,nx
    astore(i,j,k)=astore(i,j,k)+ain(i,j,k)*rwght
    enddo
    enddo
   enddo
!
end subroutine store_atmvar

subroutine store_soilvar(nx,ny,nzsoil,nstyps,ain,astore,rwght)
  implicit none
  integer nx,ny,nzsoil,nstyps
  integer i,j,k
  real rwght
  real ain(nx,ny,nzsoil,nstyps),astore(nx,ny,nzsoil,nstyps)
  integer istyp
!clt
  do istyp=0,nstyps
  do k=1,nzsoil
   do j=1,ny
    do i=1,nx
    astore(i,j,k,istyp)=astore(i,j,k,istyp)+ain(i,j,k,istyp)*rwght
    enddo
    enddo
   enddo
   enddo
end subroutine store_soilvar

subroutine store_soil2dvar(nx,ny,nstyps,ain,astore,rwght)
  implicit none
  integer nx,ny,nstyps
  integer i,j,k
  real rwght
  real ain(nx,ny,nstyps),astore(nx,ny,nstyps)
  integer istyp
!clt
  do istyp=0,nstyps
   do j=1,ny
    do i=1,nx
    astore(i,j,istyp)=astore(i,j,istyp)+ain(i,j,istyp)*rwght
    enddo
   enddo
   enddo
end subroutine store_soil2dvar

subroutine store_soil2dvars(nx,ny,ain,astore,rwght)
  implicit none
  integer nx,ny
  integer i,j
  real rwght
  real ain(nx,ny),astore(nx,ny)
  integer istyp
!clt
   do j=1,ny
    do i=1,nx
    astore(i,j)=astore(i,j)+ain(i,j)*rwght
    enddo
   enddo
end subroutine store_soil2dvars

subroutine re_store_atmvar(nx,ny,nz,nt,aout,astore)
  implicit none
  integer nx,ny,nz,nt
  integer i,j,k
  real rwght
  integer int
  real aout(nx,ny,nz,nt),astore(nx,ny,nz)
!clt
  do int=1,nt
  do k=1,nz
   do j=1,ny
    do i=1,nx
    aout(i,j,k,int)=astore(i,j,k)
    enddo
    enddo
   enddo
   enddo
!
end subroutine re_store_atmvar

subroutine re_store_soilvar(nx,ny,nzsoil,nstyps,aout,astore)
  implicit none
  integer nx,ny,nzsoil,nstyps
  integer i,j,k
  real rwght
  real aout(nx,ny,nzsoil,nstyps),astore(nx,ny,nzsoil,nstyps)
  integer istyp
!clt

  do istyp=0,nstyps
  do k=1,nzsoil
   do j=1,ny
    do i=1,nx
    aout(i,j,k,istyp) =astore(i,j,k,istyp)
    enddo
    enddo
   enddo
   enddo
!
end subroutine re_store_soilvar

subroutine re_store_soil2dvar(nx,ny,nstyps,aout,astore)
  implicit none
  integer nx,ny,nstyps
  integer i,j,k
  real rwght
  real aout(nx,ny,nstyps),astore(nx,ny,nstyps)
  integer istyp
!clt

  do istyp=0,nstyps
   do j=1,ny
    do i=1,nx
    aout(i,j,istyp) =astore(i,j,istyp)
    enddo
   enddo
   enddo
!
end subroutine re_store_soil2dvar

subroutine re_store_soil2dvars(nx,ny,aout,astore)
  implicit none
  integer nx,ny
  integer i,j
  real rwght
  real aout(nx,ny),astore(nx,ny)
!clt

   do j=1,ny
    do i=1,nx
    aout(i,j) =astore(i,j)
    enddo
   enddo
!
end subroutine re_store_soil2dvars
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE CHKSTAB                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE chkstab(mptr,nx,ny,nz,nzsoil,nstyps,                         &
           u,v,w,wcont,ptprt,pprt,qv,qscalar,tke,                       &
           ubar,vbar,ptbar,pbar,rhobar,qvbar,kmh,kmv,                   &
           x,y,z,zp,zpsoil,hterain,mapfct,j1,j2,j3,j3soil,              &
           soiltyp,stypfrct,vegtyp,lai,roufns,veg,                      &
           tsoil,qsoil,wetcanp,snowdpth,                                &
           raing,rainc,prcrate,                                         &
           radfrc,radsw,rnflx,radswnet,radlwin,                         &
           usflx,vsflx,ptsflx,qvsflx,                                   &
           tem1,tem2,tem3,tem4)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Check the stability of the time integration. If unstable, dump
!  out the model data and stop the model run.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  11/13/91.
!
!  MODIFICATION HISTORY:
!
!  6/06/92 (M. Xue)
!  Added full documentation.
!
!  7/13/92 (K. Brewster)
!  Added comment line variables to arg list to be passed to the
!  data dumping routines.
!
!  12/09/1998 (Donghai Wang)
!  Added the snow cover.
!
!  05/13/2002 (J. Brotzge)
!  Added additional soil scheme arrays.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    mptr     Grid identifier.
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!    nzsoil   Number of grid points in the soil
!
!    u        x component of velocity (m/s)
!    v        y component of velocity (m/s)
!    w        Vertical component of Cartesian velocity (m/s)
!    wcont    Contravariant vertical velocity (m/s)
!    ptprt    Perturbation potential temperature (K)
!    pprt     Perturbation pressure (Pascal)
!    qv       Water vapor specific humidity (kg/kg)
!    qc       Cloud water mixing ratio (kg/kg)
!    qr       Rainwater mixing ratio (kg/kg)
!    qi       Cloud ice mixing ratio (kg/kg)
!    qs       Snow mixing ratio (kg/kg)
!    qh       Hail mixing ratio (kg/kg)
!    tke      Turbulent Kinetic Energy ((m/s)**2)
!
!    ubar     Base state zonal velocity component (m/s)
!    vbar     Base state meridional velocity component (m/s)
!    ptbar    Base state potential temperature (K)
!    pbar     Base state pressure (Pascal)
!    rhobar   Base state density (kg/m**3)
!    qvbar    Base state water vapor specific humidity (kg/kg)
!
!    kmh      Horizontal turb. mixing coef. for momentum ( m**2/s )
!    kmv      Vertical turb. mixing coef. for momentum ( m**2/s )
!
!    x        x coordinate of grid points in physical/comp. space (m)
!    y        y coordinate of grid points in physical/comp. space (m)
!    z        z coordinate of grid points in computational space (m)
!    zp       Vertical coordinate of grid points in physical space (m)
!    zpsoil   Vertical coordinate of grid points in the soil (m)
!    hterain  Terrain height (m)
!    mapfct   Map factors at scalar, u and v points
!
!    j1       Coordinate transformation Jacobian -d(zp)/dx
!    j2       Coordinate transformation Jacobian -d(zp)/dy
!    j3       Coordinate transformation Jacobian  d(zp)/dz
!    j3soil   Coordinate transformation Jacobian -d(zpsoil)/dz
!
!    soiltyp  Soil type
!    stypfrct  Soil type fraction
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
!    prcrate  Precipitation rates
!
!    radfrc   Radiation forcing (K/s)
!    radsw    Solar radiation reaching the surface
!    rnflx    Net radiation flux absorbed by surface
!    radswnet Net solar radiation
!    radlwin  Incoming longwave radiation
!
!    usflx    Surface flux of u-momentum (kg/(m*s**2))
!    vsflx    Surface flux of v-momentum (kg/(m*s**2))
!    ptsflx   Surface heat flux (K*kg/(m**2 * s ))
!    qvsflx   Surface moisture flux of (kg/(m**2 * s))
!
!  OUTPUT:
!
!    None.
!
!  WORK ARRAYS:
!
!    tem1     Temporary work array.
!    tem2     Temporary work array.
!    tem3     Temporary work array.
!    tem4     Temporary work array.
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
  INCLUDE 'phycst.inc'
  INCLUDE 'globcst.inc'

  INCLUDE 'timelvls.inc'
  INCLUDE 'mp.inc'

!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  INTEGER :: mptr              ! Grid identifier.

  INTEGER :: nx,ny,nz          ! Number of grid points in 3 directions
  INTEGER :: nzsoil            ! Number of grid points in the soil

  REAL :: u     (nx,ny,nz,nt)  ! Total u-velocity (m/s)
  REAL :: v     (nx,ny,nz,nt)  ! Total v-velocity (m/s)
  REAL :: w     (nx,ny,nz,nt)  ! Total w-velocity (m/s)
  REAL :: wcont (nx,ny,nz)     ! Contravariant vertical velocity (m/s)
  REAL :: ptprt (nx,ny,nz,nt)  ! Perturbation potential temperature (K)
  REAL :: pprt  (nx,ny,nz,nt)  ! Perturbation pressure (Pascal)

  REAL :: qv    (nx,ny,nz,nt)  ! Water vapor specific humidity (kg/kg)

  REAL :: qscalar(nx,ny,nz,nt,nscalar)

  REAL :: tke   (nx,ny,nz,nt)  ! turbulent kinetic energy

  REAL :: ubar  (nx,ny,nz)     ! Base state u-velocity (m/s)
  REAL :: vbar  (nx,ny,nz)     ! Base state v-velocity (m/s)
  REAL :: ptbar (nx,ny,nz)     ! Base state potential temperature (K)
  REAL :: pbar  (nx,ny,nz)     ! Base state pressure (Pascal)
  REAL :: rhobar(nx,ny,nz)     ! Base state air density (kg/m**3)
  REAL :: qvbar (nx,ny,nz)     ! Base state water vapor specific humidity
                               ! (kg/kg)

  REAL :: kmh   (nx,ny,nz)     ! Horizontal turb. mixing coef. for
                               ! momentum. ( m**2/s )
  REAL :: kmv   (nx,ny,nz)     ! Vertical turb. mixing coef. for
                               ! momentum. ( m**2/s )

  REAL :: x     (nx)           ! The x-coord. of the physical and
                               ! computational grid. Defined at u-point.
  REAL :: y     (ny)           ! The y-coord. of the physical and
                               ! computational grid. Defined at v-point.
  REAL :: z     (nz)           ! The z-coord. of the computational grid.
                               ! Defined at w-point on the staggered grid.
  REAL :: zp    (nx,ny,nz)     ! The physical height coordinate defined at
                               ! w-point of the staggered grid.
  REAL :: zpsoil (nx,ny,nzsoil)      ! The physical height coordinate defined at
                               ! w-point of the soil.

  REAL :: hterain(nx,ny)       ! Terrain height.
  REAL :: mapfct(nx,ny,8)      ! Map factors at scalar, u and v points

  REAL :: j1    (nx,ny,nz)     ! Coordinate transformation Jacobian defined as
                               ! - d( zp )/d( x )
  REAL :: j2    (nx,ny,nz)     ! Coordinate transformation Jacobian defined as
                               ! - d( zp )/d( y )
  REAL :: j3    (nx,ny,nz)     ! Coordinate transformation Jacobian defined as
                               ! d( zp )/d( z )
  REAL :: j3soil (nx,ny,nzsoil)! Coordinate transformation Jacobian defined as
                               ! d( zpsoil )/d( z )

  INTEGER :: nstyps                   ! Number of soil types
  INTEGER :: soiltyp (nx,ny,nstyps)   ! Soil type
  REAL    :: stypfrct(nx,ny,nstyps)   ! Soil type fraction
  INTEGER :: vegtyp (nx,ny)           ! Vegetation type
  REAL    :: lai    (nx,ny)           ! Leaf Area Index
  REAL    :: roufns (nx,ny)           ! Surface roughness
  REAL    :: veg    (nx,ny)           ! Vegetation fraction

  REAL :: tsoil  (nx,ny,nzsoil,0:nstyps)     ! Deep soil temperature (K)
  REAL :: qsoil  (nx,ny,nzsoil,0:nstyps)     ! Deep soil moisture (m**3/m**3)
  REAL :: wetcanp(nx,ny,0:nstyps)     ! Canopy water amount
  REAL :: snowdpth(nx,ny)             ! Snow depth (m)

  REAL :: raing(nx,ny)         ! Grid supersaturation rain
  REAL :: rainc(nx,ny)         ! Cumulus convective rain
  REAL :: prcrate(nx,ny,4)     ! precipitation rate (kg/(m**2*s))
                               ! prcrate(1,1,1) = total precip. rate
                               ! prcrate(1,1,2) = grid scale precip. rate
                               ! prcrate(1,1,3) = cumulus precip. rate
                               ! prcrate(1,1,4) = microphysics precip. rate

  REAL :: radfrc(nx,ny,nz)     ! Radiation forcing (K/s)
  REAL :: radsw (nx,ny)        ! Solar radiation reaching the surface
  REAL :: rnflx (nx,ny)        ! Net radiation flux absorbed by surface
  REAL :: radswnet (nx,ny)     ! Net solar radiation, SWin - SWout
  REAL :: radlwin  (nx,ny)     ! Incoming longwave radiation


  REAL :: usflx (nx,ny)        ! Surface flux of u-momentum (kg/(m*s**2))
  REAL :: vsflx (nx,ny)        ! Surface flux of v-momentum (kg/(m*s**2))
  REAL :: ptsflx(nx,ny)        ! Surface heat flux (K*kg/(m*s**2))
  REAL :: qvsflx(nx,ny)        ! Surface moisture flux (kg/(m**2*s))

  REAL :: tem1  (nx,ny,nz)     ! Temporary work array
  REAL :: tem2  (nx,ny,nz)     ! Temporary work array
  REAL :: tem3  (nx,ny,nz)     ! Temporary work array
  REAL :: tem4  (nx,ny,nz)     ! Temporary work array
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k
  INTEGER :: tlevel
  REAL    :: absmax,vellmt
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  mgrid = mptr

  tlevel = tpresent

!
!-----------------------------------------------------------------------
!
!  Check for instabilty of model integration.
!
!  When instability is detected, stop the model run and dump out the
!  model data.
!
!-----------------------------------------------------------------------
!
  absmax = MAX(ABS(u(1,1,1,tlevel)),                                    &
           ABS(v(1,1,1,tlevel)),ABS(w(1,1,1,tlevel)))
!
  DO k=1,nz-1
    DO j=1,ny-1
      DO i=1,nx
        absmax = MAX(absmax,ABS(u(i,j,k,tlevel)))
      END DO
    END DO
  END DO

  DO k=1,nz-1
    DO j=1,ny
      DO i=1,nx-1
        absmax = MAX(absmax,ABS(v(i,j,k,tlevel)))
      END DO
    END DO
  END DO

  DO k=1,nz
    DO j=1,ny-1
      DO i=1,nx-1
        absmax = MAX(absmax,ABS(w(i,j,k,tlevel)))
      END DO
    END DO
  END DO
!
!-----------------------------------------------------------------------
!
!  Set an upper limit on maximum velocity component, above which
!  the integeration is regarded as unstable.
!
!-----------------------------------------------------------------------
!
  vellmt = 150.0

  CALL mpmaxr(absmax)

  IF( absmax > vellmt ) THEN

    IF (myproc == 0) WRITE(6,'(/1x,a,f5.1,a,/1x,a,i6,/1x,a)')           &
        'Maximum velocity component exceeded ',vellmt,' m/s,',          &
        'Time integration is unstable. Job aborted at the end of step', &
        nstep,' Fields at this time dumped out.'

!    IF (mp_opt > 0) THEN
!      CALL arpsstop("arpsstop called from CHKSTAB for MP version",1)
!    END IF

    CALL set_acct(output_acct)
    CALL abortdmp(mptr,nx,ny,nz,nzsoil,nstyps,                          &
                  u,v,w,wcont,ptprt,pprt,qv,qscalar,tke,                &
                  ubar,vbar,ptbar,pbar,rhobar,qvbar,kmh,kmv,            &
                  x,y,z,zp,hterain,zpsoil,mapfct, j1,j2,j3,j3soil,      &
                  soiltyp,stypfrct,vegtyp,lai,roufns,veg,               &
                  tsoil,qsoil,wetcanp,snowdpth,                         &
                  raing,rainc,prcrate,                                  &
                  radfrc,radsw,rnflx,radswnet,radlwin,                  &
                  usflx,vsflx,ptsflx,qvsflx,                            &
                  tem1,tem2,tem3, tem4)

    CALL arpsstop("arpsstop called from CHKSTAB after abortdmp",1)

  END IF

  RETURN
END SUBROUTINE chkstab
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE MAXMIN                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE maxmin(mptr,nx,ny,nz,nzsoil,tlevel,rhobar,                   &
           u,v,w,wcont,ptprt,pprt,qv,qscalar,tke,kmh,kmv,               &
           x,y,z,zp,zpsoil,mapfct,                                      &
           tsoil,qsoil,wetcanp,                                         &
           tem1,tem2,tem3,tem4)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate the maximum and minimum of the time dependent fields and
!  write them into a file named runname(1:lfnkey).maxmin.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!
!
!  MODIFICATION HISTORY:
!
!  6/06/92 (M. Xue)
!  Added full documentation.
!
!  4/10/93 (K. Droegemeier and MX)
!  Added max.&min. calcualtions of reflectivity, vorticity and max.
!  low-level winds.
!
!  12/6/95 (J. Zong amd M. Xue)
!  Call subroutine REFLEC to calculate reflectivity (in dBz) rather
!  than use an explicit expression to ease modifications to reflectivity
!  formula (e.g., including contributions of snow and graupel/hail to
!  reflectivity). While doing so, the reflectivity in dBz must be
!  converted to Z in calculating VIL.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    mptr     Grid identifier.
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!    nzsoil   Number of grid points in the soil
!
!    tlevel   The time level at which the data are printed.
!
!    rhobar   Base state density (kg/m**3)
!    u        x component of velocity (m/s)
!    v        y component of velocity (m/s)
!    w        Vertical component of Cartesian velocity (m/s)
!    wcont    Contravariant vertical velocity (m/s)
!    ptprt    Perturbation potential temperature (K)
!    pprt     Perturbation pressure (Pascal)
!    qv       Water vapor specific humidity (kg/kg)
!    qc       Cloud water mixing ratio (kg/kg)
!    qr       Rainwater mixing ratio (kg/kg)
!    qi       Cloud ice mixing ratio (kg/kg)
!    qs       Snow mixing ratio (kg/kg)
!    qh       Hail mixing ratio (kg/kg)
!    tke      Turbulent Kinetic Energy ((m/s)**2)
!
!    x        x coordinate of grid points in physical/comp. space (m)
!    y        y coordinate of grid points in physical/comp. space (m)
!    z        z coordinate of grid points in computational space (m)
!    zp       Vertical coordinate of grid points in physical space (m)
!    zpsoil   Vertical coordinate of grid points in the soil (m)
!    mapfct   Map factors at scalar, u and v points
!
!    tsoil    Soil temperature (K)
!    qsoil    Soil moisture (m**3/m**3)
!    wetcanp  Canopy water amount
!
!  OUTPUT:
!
!    None.
!
!  WORK ARRAYS:
!
!    tem1     Radar reflectivity factor (dBz)
!    tem2     Vertical vorticity (per second)
!    tem3     Vertically-Integrated Liquid (kg/m**2)
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'grid.inc'          ! Grid & map parameters.
  INCLUDE 'mp.inc'            ! Message passing parameters.
  INCLUDE 'timelvls.inc'
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  INTEGER :: mptr              ! Grid identifier.
  INTEGER :: nx,ny,nz          ! Number of grid points in 3 directions.
  INTEGER :: nzsoil            ! Number of grid points in the soil.
  INTEGER :: tlevel            ! Time level at which data are processed.

  REAL :: rhobar(nx,ny,nz)     ! Base state density (kg/m**3)
  REAL :: u     (nx,ny,nz,nt)  ! Total u-velocity (m/s)
  REAL :: v     (nx,ny,nz,nt)  ! Total v-velocity (m/s)
  REAL :: w     (nx,ny,nz,nt)  ! Total w-velocity (m/s)
  REAL :: wcont (nx,ny,nz)     ! Contravariant vertical velocity (m/s)
  REAL :: ptprt (nx,ny,nz,nt)  ! Perturbation potential temperature (K)
  REAL :: pprt  (nx,ny,nz,nt)  ! Perturbation pressure (Pascal)

  REAL :: qv    (nx,ny,nz,nt)  ! Water vapor specific humidity (kg/kg)

  REAL :: qscalar(nx,ny,nz,nt,nscalar)

  REAL :: tke   (nx,ny,nz,nt)  ! Turbulent Kinetic Energy ((m/s)**2)
  REAL :: kmh   (nx,ny,nz)     ! Horizontal turb. mixing coef. for
                               ! momentum. ( m**2/s )
  REAL :: kmv   (nx,ny,nz)     ! Vertical turb. mixing coef. for
                               ! momentum. ( m**2/s )

  REAL :: x     (nx)           ! The x-coord. of the physical and
                               ! computational grid. Defined at u-point.
  REAL :: y     (ny)           ! The y-coord. of the physical and
                               ! computational grid. Defined at v-point.
  REAL :: z     (nz)           ! The z-coord. of the computational grid.
                               ! Defined at w-point on the staggered grid.
  REAL :: zp    (nx,ny,nz)     ! The physical height coordinate defined at
                               ! w-point of the staggered grid.
  REAL :: zpsoil(nx,ny,nzsoil)       ! The physical height coordinate defined at
                               ! w-point of the soil.
  REAL :: mapfct(nx,ny,8)      ! Map factors at scalar, u and v points

  REAL :: tsoil (nx,ny,nzsoil) ! Soil temperature (K)
  REAL :: qsoil (nx,ny,nzsoil) ! Soil moisture (m**3/m**3)
  REAL :: wetcanp(nx,ny)       ! Canopy water amount

  REAL :: tem1  (nx,ny,nz)     ! Radar reflectivity factor (dBz)
  REAL :: tem2  (nx,ny,nz)     ! Vertical vorticity (per second)
  REAL :: tem3  (nx,ny,nz)     ! Vertically-integrated liquid (kg/m**2)
  REAL :: tem4  (nx,ny,nz)     ! temporary work array
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=256      ) :: maxfn, outdirname  ! Name of max./min. file
  CHARACTER (LEN=40       ) :: fmtstr
  INTEGER :: lmaxfn, tnq
  INTEGER :: i,j,k
  INTEGER :: nq
  INTEGER :: istat
  INTEGER :: imax,jmax,kmax,imin,jmin,kmin, klevel
  REAL    :: umax,vmax,wmax,ptmax,pmax,qvmax,qscalarmax(nscalar)
  REAL    :: tkemin,tkemax
  REAL    :: amin, amax
  REAL    :: umin,vmin,wmin,ptmin,pmin,qvmin,qscalarmin(nscalar)
  REAL    :: zrefmax,refmax,refmin
  REAL    :: udtdxmax,udtdxmin,vdtdymax,vdtdymin,wdtdzmax,wdtdzmin
  REAL    :: uvtem,utem,vtem,ucnmax,vcnmax,wcnmax,uvcnmax,uvcnmin
  INTEGER :: iumax, jumax, kumax
  INTEGER :: ivmax, jvmax, kvmax
  INTEGER :: iwmax, jwmax, kwmax
  INTEGER :: iuvmax, juvmax, kuvmax

  REAL :: uvmax, uvmaxdr, uvmin, eps, pi
  REAL :: weight,wreflk,frsvth
  REAL :: zlevel, pastim, denwater,vtr
  REAL :: acumrain             ! Accumulated surface rain fall
                               ! in entire domain (mm)

  REAL    :: qcmax,qrmax,qimax,qsmax,qhmax  ! introduced for generating
                               ! similar maxmin file as previous version

  REAL :: vorlx,vorln,vorux,vorun
                               ! Lower and upper level vorticity max/min.
  REAL :: vil                  ! Vertical integrated max reflectivity
                               ! averaged over 5 grid points
  REAL :: gumove,gvmove

  INTEGER :: ncalls(mgrdmax), nchmax1(mgrdmax)
  SAVE acumrain,pastim
  SAVE ncalls,nchmax1
  DATA pastim,acumrain /0.0,0.0/
  DATA ncalls /mgrdmax*0/

!  INTEGER :: iproc, jproc, iloc, jloc, isource

  CHARACTER(LEN=4) :: upcase
  INTEGER :: istatus
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

!-----------------------------------------------------------------------
!
!  Find the first scalar level whose height is equal to or higher than
!  a given value zlevel. We assume that all grid points on a given
!  level have the same height.
!
!-----------------------------------------------------------------------

  zlevel = 2000.0

  DO k=2,nz-1
    IF( 0.5*( zp(2,2,k)+zp(2,2,k+1) ) >= zlevel ) EXIT
  END DO

  klevel = MIN(nz-1,MAX(2, k-1))
  CALL mpupdatei(klevel,1)

!
!-----------------------------------------------------------------------
!  Compute Courant numbers for monitoring time integration stability
!-----------------------------------------------------------------------
!
  IF (myproc == 0 ) WRITE(6,*) ' '

!-----------------------------------------------------------------------
!   Max Courant number in x direction
!-----------------------------------------------------------------------

  DO k=2,nz-2
    DO j=2,ny-2
      DO i=2,nx-1
        tem3(i,j,k) = dtbig*ABS(u(i,j,k,2))/dx * mapfct(i,j,2)
      END DO
    END DO
  END DO

  CALL a3dmax(tem3,1,nx,2,nx-1,1,ny,2,ny-2,1,nz,2,nz-2,                 &
              udtdxmax,udtdxmin, iumax,jumax,kumax, imin,jmin,kmin)

  ucnmax = udtdxmax
  IF (myproc == 0) WRITE(6,'((1x,a,f10.4,3(a,i4),a,f10.3))')            &
       'ucourant =',udtdxmax,' at i=',iumax,', j=',jumax,', k=',kumax
!       ', u=',u(iumax,jumax,kumax,2)  ! comment out for PGF90

!-----------------------------------------------------------------------
!   Max Courant number in y direction
!-----------------------------------------------------------------------

  DO k=2,nz-2
    DO j=2,ny-1
      DO i=2,nx-2
        tem3(i,j,k) = dtbig*ABS(v(i,j,k,2))/dy * mapfct(i,j,3)
      END DO
    END DO
  END DO
  CALL a3dmax(tem3,1,nx,2,nx-2,1,ny,2,ny-1,1,nz,2,nz-2  ,               &
              vdtdymax,vdtdymin, ivmax,jvmax,kvmax, imin,jmin,kmin)

  vcnmax = vdtdymax
  IF (myproc == 0) WRITE(6,'((1x,a,f10.4,3(a,i4),a,f10.3))')            &
       'vcourant =',vdtdymax,' at i=',ivmax,', j=',jvmax,', k=',kvmax
!       ', v=',v(ivmax,jvmax,kvmax,2)  ! comment our for PGF90

!-----------------------------------------------------------------------
!   Max Courant number in z direction for each layer
!-----------------------------------------------------------------------

  DO k=2,nz-1
    DO j=2,ny-2
      DO i=2,nx-2
        tem3(i,j,k) = dtbig*ABS(wcont(i,j,k))/dz
      END DO
    END DO
  END DO
  CALL a3dmax(tem3,1,nx,2,nx-2,1,ny,2,ny-2,1,nz,2,nz-1,                 &
              wdtdzmax,wdtdzmin, iwmax,jwmax,kwmax, imin,jmin,kmin)

  wcnmax = wdtdzmax
  IF (myproc == 0) WRITE(6,'((1x,a,f10.4,3(a,i4),2(a,f10.3)))')         &
       'wcourant =',wdtdzmax,' at i=',iwmax,', j=',jwmax,', k=',kwmax
!       ', wcont=',wcont(iwmax,jwmax,kwmax),', w=',w(iwmax,jwmax,kwmax,2)
! Commented out for PGF90 out of bound error

!-----------------------------------------------------------------------
!   Max Courant number in horizontal direction for each layer
!-----------------------------------------------------------------------

  DO k=2,nz-2
    DO j=2,ny-2
      DO i=2,nx-2
        utem = 0.5*(u(i,j,k,2)+u(i+1,j,k,2))
        vtem = 0.5*(v(i,j,k,2)+v(i,j+1,k,2))
        uvtem = sqrt(utem*utem+vtem*vtem)
        tem3(i,j,k) = dtbig*utem*vtem/(dx*max(1.0,uvtem))*mapfct(i,j,1)
      END DO
    END DO
  END DO
  CALL a3dmax(tem3,1,nx,2,nx-2,1,ny,2,ny-2,1,nz,2,nz-2,                 &
              uvcnmax,uvcnmin, iuvmax,juvmax,kuvmax,imin,jmin,kmin)

  IF (myproc == 0) WRITE(6,'((1x,a,f10.4,3(a,i4),a,f10.4))')            &
       'uvcourant=',uvcnmax,' at i=',iuvmax,', j=',juvmax,', k=',kuvmax
!       ', mapfct=',mapfct(iuvmax,juvmax,1)

!-----------------------------------------------------------------------
! Domain maximum Courant number in all three directions
!-----------------------------------------------------------------------
!
! IF (myproc == 0)  &
! WRITE(6,'(4(a,f10.4))') 'ucnmax =', ucnmax, ', vcnmax =', vcnmax, ',   &
!   wcnmax =', wcnmax,', uvcnmax=',uvcnmax
!
!-----------------------------------------------------------------------
!
!  Compute the radar reflectivity factor following Kessler (1969).
!  Here, arg=Z (mm**6/m**3), and dBz = 10log10 (arg).
!
!-----------------------------------------------------------------------
!
  IF (P_QR > 0) THEN
    tem2(:,:,:) = qscalar(:,:,:,tlevel,P_QR)
  ELSE
    tem2(:,:,:) = 0.0
  END IF

  IF (P_QS > 0) THEN
    tem3(:,:,:) = qscalar(:,:,:,tlevel,P_QS)
  ELSE
    tem3(:,:,:) = 0.0
  END IF

  IF (P_QH > 0) THEN
    tem4(:,:,:) = qscalar(:,:,:,tlevel,P_QH)
  ELSE
    tem4(:,:,:) = 0.0
  END IF

  CALL reflec(nx,ny,nz, rhobar, tem2, tem3, tem4, tem1)

!
!-----------------------------------------------------------------------
!
!  Compute the verticallly-integrated liquid water using the
!  formula given by Stewart (1991) - NOAA Tech Memo NWR SR-136,
!  National Weather Service, 20 pp.
!
!  VIL is obtained by averaging over a 3 x 3 km square centered on
!  the domain-maximum reflectivity.  It is computed only where the
!  reflectivity is non-zero.
!
!-----------------------------------------------------------------------
!
!  First, find the domain-maximum reflectivity factor
!
!-----------------------------------------------------------------------

  CALL a3dmax(tem1,1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz-1,                 &
              refmax,refmin, imax,jmax,kmax, imin,jmin,kmin)


  IF (mp_opt == 0) THEN
    zrefmax=MAX(0.0,(zp(imax,jmax,kmax)+zp(imax,jmax,kmax+1))*0.5)

  END IF

  IF (myproc == 0) WRITE(6,'((1x,a,g25.12,3(a,i4)))')                   &
        'refmax=',refmax,' at i=',imax,', j=',jmax,', k=',kmax
!
!-----------------------------------------------------------------------
!
!  Compute a horizontally-averaged reflectivity centered on the max,
!  and sum vertically. wreflk = weighted reflectivity at level k.
!
!  First, we recompute the "arg" function in mm**6/m**3.  VIL does
!  not use reflectivity in dBz, but rather "arg", or Z.  At this
!  point, tem1 is overwritten and no longer contains 10log10 (Z).
!
!-----------------------------------------------------------------------
!

  DO k = 1,nz-1
    DO j = 1,ny-1
      DO i = 1,nx-1
        tem1(i,j,k) = 10.0**( 0.1*tem1(i,j,k) )
      END DO
    END DO
  END DO

  frsvth = 4./7.
  weight = 1./9.
  vil = 0.0

  imax = MIN( nx-2,MAX(2,imax) )
  jmax = MIN( ny-2,MAX(2,jmax) )
  kmax = MIN( nz-2,MAX(2,kmax) )

  DO k = 2,nz-2
    wreflk  = weight*(tem1(imax,jmax,k)                                 &
              +tem1(imax+1,jmax,k)+tem1(imax-1,jmax,k)                  &
              +tem1(imax,jmax+1,k)+tem1(imax,jmax-1,k)                  &
              +tem1(imax+1,jmax+1,k)+tem1(imax-1,jmax+1,k)              &
              +tem1(imax+1,jmax-1,k)+tem1(imax-1,jmax-1,k))
    vil = vil+dz*3.44E-6* wreflk**frsvth
  END DO

!-----------------------------------------------------------------------
!
!  Compute the surface accumulated precipitation.  The depth of water
!  accumulated per timestep is given by:
!
!  depth (m) = terminal velocity * qr * dt * (rhobar/denwater)
!
!  This equation is derived by temporally integrating the vertical
!  flux of rain through the lowest model level.  The total is given
!  for the entire computational domain.
!
!-----------------------------------------------------------------------

  denwater = 1000.0         ! Density of liquid water.

  IF (P_QR > 0) THEN
    DO j=1,ny-1
      DO i=1,nx-1

        vtr = 36.34 * (0.001*rhobar(i,j,2)*                             &
                       MAX(0.0,qscalar(i,j,2,tlevel,P_QR)))**0.1364     &
            * SQRT(rho0/rhobar(i,j,2))

        acumrain = acumrain + vtr*qscalar(i,j,2,tlevel,P_QR)*           &
                   (curtim-pastim)*rhobar(i,j,2)/denwater * 1000.0
      END DO
    END DO
  END IF
!
!-----------------------------------------------------------------------
!
!  Compute the vertical vorticity, which is defined at the corner
!  point of grid box.
!
!-----------------------------------------------------------------------

  DO k=2,nz-1
    DO j=2,ny-1
      DO i=2,nx-1
        tem2(i,j,k)=                                                    &
            (v(i,j,k,tlevel)-v(i-1,j,k,tlevel))/dx-                     &
            (u(i,j,k,tlevel)-u(i,j-1,k,tlevel))/dy
      END DO
    END DO
  END DO


  DO k=2,nz-1
    DO j=2,ny-1
      tem2( 1,j,k)=tem2(   2,j,k)
      tem2(nx,j,k)=tem2(nx-1,j,k)
    END DO
  END DO

  DO k=2,nz-1
    DO i=1,nx
      tem2(i, 1,k)=tem2(i,   2,k)
      tem2(i,ny,k)=tem2(i,ny-1,k)
    END DO
  END DO

  DO j=1,ny
    DO i=1,nx
      tem2(i,j, 1)=tem2(i,j,   2)
      tem2(i,j,nz)=tem2(i,j,nz-1)
    END DO
  END DO


!-----------------------------------------------------------------------
!
!  Pull out the domain-maximum vertical vorticity below z = zlevel.
!  vrbx = max (x) vorticity at or below z = zlevel.
!
!-----------------------------------------------------------------------
!
  CALL a3dmax(tem2,1,nx,2,nx-1,1,ny,2,ny-1,1,nz,2,klevel,               &
              vorlx,vorln,imax,jmax,kmax, imin,jmin,kmin)

  IF (myproc == 0) WRITE(6,'((1x,a,g25.12,3(a,i4)))')                   &
       'vorlx=',vorlx,' at i=',imax,', j=',jmax,', k=',kmax

!-----------------------------------------------------------------------
!
!  Pull out the domain-maximum vertical vorticity above z = zlevel.
!  vrbx = max (x) vorticity at or above z = zlevel.
!
!-----------------------------------------------------------------------

  CALL a3dmax(tem2,1,nx,2,nx-1,1,ny,2,ny-1,1,nz,klevel+1,nz-2,          &
              vorux,vorun,imax,jmax,kmax, imin,jmin,kmin)

  IF (myproc == 0) WRITE(6,'((1x,a,g25.12,3(a,i4)))')                   &
       'vorux=',vorux,' at i=',imax,', j=',jmax,', k=',kmax
!
!-----------------------------------------------------------------------
!
!  The maximum ground-relative surface winds below z = zlevel.
!
!-----------------------------------------------------------------------
!
  IF ( grdtrns == 0 ) THEN
    gumove = 0.0
    gvmove = 0.0
  ELSE
    gumove = umove
    gvmove = vmove
  END IF

  eps = 1.0E-30
  pi = ATAN( 1.0 )*4
  DO k = 1,nz-1
    DO j = 1,ny-1
      DO i = 1,nx-1
        tem1(i,j,k)= SQRT(                                              &
             ((u(i,j,k,tlevel)+u(i+1,j,k,tlevel))*0.5+gumove)**2+       &
             ((v(i,j,k,tlevel)+v(i,j+1,k,tlevel))*0.5+gvmove)**2)
      END DO
    END DO
  END DO

  CALL a3dmax(tem1,1,nx,1,nx-1,1,ny,1,ny-1,1,nz,2,klevel,               &
              uvmax,uvmin, imax,jmax,kmax, imin,jmin,kmin)
!
!-----------------------------------------------------------------------
!
!  The direction of the maximum ground-relative surface wind
!
!-----------------------------------------------------------------------
!
  IF (mp_opt == 0) THEN

    uvmaxdr = 180.0/pi * ASIN( MAX(-1.0, MIN(1.0,                       &
          ((u(imax,jmax,kmax,tlevel)+u(imax+1,jmax,kmax,tlevel))        &
          *0.5+gumove)/(tem1(imax,jmax,kmax)+eps))) )

    IF( ((v(imax,jmax,kmax,tlevel)+v(imax,jmax+1,kmax,tlevel))          &
        *0.5+gvmove) < 0.0 ) uvmaxdr = uvmaxdr + 360.0

  ELSE

! It is too expensive for just output, if you really desired same output
! turn the following block on

!    CALL mpupdatei(imax,1)
!    CALL mpupdatei(jmax,1)
!
!    iproc = (imax-2)/(nx-3)+1
!    jproc = (jmax-2)/(ny-3)+1
!    IF (iproc < 1) iproc = 1
!    IF (jproc < 1) jproc = 1
!    isource = (jproc-1)*nproc_x+iproc-1
!
!    iloc = MOD((imax-2),(nx-3)) + 2
!    jloc = MOD((jmax-2),(ny-3)) + 2
!
!    IF (loc_x == iproc .AND. loc_y == jproc) THEN
!      uvmaxdr = 180.0/pi * ASIN( MAX(-1.0, MIN(1.0,                     &
!          ((u(iloc,jloc,kmax,tlevel)+u(iloc+1,jloc,kmax,tlevel))        &
!          *0.5+gumove)/(tem1(iloc,jloc,kmax)+eps))) )
!
!      IF( ((v(iloc,jloc,kmax,tlevel)+v(iloc,jloc+1,kmax,tlevel))        &
!          *0.5+gvmove) < 0.0 ) uvmaxdr = uvmaxdr + 360.0
!      CALL mpsendr(uvmaxdr,1,0,300,istat)
!    END IF
!
!    IF (myproc == 0 .AND. isource /= 0) CALL mprecvr(uvmaxdr,1,isource,300,istat)

    uvmaxdr = 90.0

  END IF

!
!-----------------------------------------------------------------------
!
!  Make the wind direction conform to the weather report standard
!  i.e. zero degree for the northerly wind.
!
!-----------------------------------------------------------------------
!
  uvmaxdr = 90.0 - uvmaxdr
  IF( uvmaxdr < 0.0) uvmaxdr = 360.0 + uvmaxdr


  IF (myproc == 0) THEN
    WRITE(6,'((1x,a,g25.12,3(a,i4)))')                                  &
       'uvmax=',uvmax,' at i=',imax,', j=',jmax,', k=',kmax

    WRITE(6,'((1x,a,g25.12,3(a,i4)))')                                  &
       'uvmaxdr=',uvmaxdr,' at i=',imax,', j=',jmax,', k=',kmax
  END IF


  CALL a3dmax(u(1,1,1,tlevel),1,nx,1,nx,                                &
              1,ny,1,ny-1,1,nz,1,nz-1,                                  &
              umax,umin, imax,jmax,kmax, imin,jmin,kmin)

  IF (myproc == 0) WRITE(6,'(2(1x,a,g25.12,3(a,i4)))')                  &
       'umin =',umin,' at i=',imin,', j=',jmin,', k=',kmin,             &
       'umax =',umax,' at i=',imax,', j=',jmax,', k=',kmax

  CALL a3dmax(v(1,1,1,tlevel),1,nx,1,nx-1,                              &
              1,ny,1,ny,1,nz,1,nz-1,                                    &
              vmax,vmin, imax,jmax,kmax, imin,jmin,kmin)

  IF (myproc == 0) WRITE(6,'(2(1x,a,g25.12,3(a,i4)))')                  &
       'vmin =',vmin,' at i=',imin,', j=',jmin,', k=',kmin,             &
       'vmax =',vmax,' at i=',imax,', j=',jmax,', k=',kmax

  CALL a3dmax(w(1,1,1,tlevel),1,nx,1,nx-1,                              &
              1,ny,1,ny-1,1,nz,1,nz,                                    &
              wmax,wmin, imax,jmax,kmax, imin,jmin,kmin)

  IF (myproc == 0) WRITE(6,'(2(1x,a,g25.12,3(a,i4)))')                  &
       'wmin =',wmin,' at i=',imin,', j=',jmin,', k=',kmin,             &
       'wmax =',wmax,' at i=',imax,', j=',jmax,', k=',kmax

  CALL a3dmax(ptprt(1,1,1,tlevel),1,nx,1,nx-1,                          &
              1,ny,1,ny-1,1,nz,1,nz-1,                                  &
              ptmax,ptmin, imax,jmax,kmax, imin,jmin,kmin)

  IF (myproc == 0) WRITE(6,'(2(1x,a,g25.12,3(a,i4)))')                  &
       'ptprtmin=',ptmin,' at i=',imin,', j=',jmin,', k=',kmin,         &
       'ptprtmax=',ptmax,' at i=',imax,', j=',jmax,', k=',kmax

  CALL a3dmax(pprt(1,1,1,tlevel),1,nx,1,nx-1,                           &
              1,ny,1,ny-1,1,nz,1,nz-1,                                  &
              pmax,pmin, imax,jmax,kmax, imin,jmin,kmin)


  IF (myproc == 0) WRITE(6,'(2(1x,a,g25.10,3(a,i4)))')                  &
       'pprtmin =',pmin,' at i=',imin,', j=',jmin,', k=',kmin,          &
       'pprtmax =',pmax,' at i=',imax,', j=',jmax,', k=',kmax

  CALL a3dmax(qv(1,1,1,tlevel),1,nx,1,nx-1,                             &
              1,ny,1,ny-1,1,nz,1,nz-1,                                  &
              qvmax,qvmin, imax,jmax,kmax, imin,jmin,kmin)

  IF (myproc == 0) WRITE(6,'(2(1x,a,g25.12,3(a,i4)))')                  &
       'qvmin=',qvmin,' at i=',imin,', j=',jmin,', k=',kmin,            &
       'qvmax=',qvmax,' at i=',imax,', j=',jmax,', k=',kmax

  DO nq = 1,nscalar
    CALL a3dmax(qscalar(1,1,1,tlevel,nq),1,nx,1,nx-1,                   &
                1,ny,1,ny-1,1,nz,1,nz-1,                                &
                qscalarmax(nq),qscalarmin(nq),imax,jmax,kmax, imin,jmin,kmin)

    IF (myproc == 0) WRITE(6,'(2(1x,a,g25.12,3(a,i4)))')                &
  TRIM(qnames(nq))//'min=',qscalarmin(nq),' at i=',imin,', j=',jmin,', k=',kmin, &
  TRIM(qnames(nq))//'max=',qscalarmax(nq),' at i=',imax,', j=',jmax,', k=',kmax

  END DO

  IF( tkeout == 1 ) THEN
    CALL a3dmax(tke(1,1,1,tlevel),1,nx,1,nx-1,                          &
                1,ny,1,ny-1,1,nz,1,nz-1,                                &
                tkemax,tkemin, imax,jmax,kmax, imin,jmin,kmin)
    IF (myproc == 0) WRITE(6,'(2(1x,a,g25.12,3(a,i4)))')                &
         'tkemin=',tkemin,' at i=',imin,', j=',jmin,', k=',kmin,        &
         'tkemax=',tkemax,' at i=',imax,', j=',jmax,', k=',kmax
  END IF

  IF( trbout == 1 ) THEN
    CALL a3dmax(kmh,1,nx,1,nx-1,                                        &
                1,ny,1,ny-1,1,nz,1,nz-1,                                &
                amax,amin, imax,jmax,kmax, imin,jmin,kmin)
    IF (myproc == 0) WRITE(6,'(2(1x,a,g25.12,3(a,i4)))')                &
         'kmhmin=',amin,' at i=',imin,', j=',jmin,', k=',kmin,          &
         'kmhmax=',amax,' at i=',imax,', j=',jmax,', k=',kmax

    CALL a3dmax(kmv,1,nx,1,nx-1,                                        &
                1,ny,1,ny-1,1,nz,1,nz-1,                                &
                amax,amin, imax,jmax,kmax, imin,jmin,kmin)
    IF (myproc == 0) WRITE(6,'(2(1x,a,g25.12,3(a,i4)))')                &
         'kmvmin=',amin,' at i=',imin,', j=',jmin,', k=',kmin,          &
         'kmvmax=',amax,' at i=',imax,', j=',jmax,', k=',kmax
  END IF

  IF( sfcphy /= 0 ) THEN

    CALL a3dmax(tsoil,1,nx,1,nx-1,1,ny,1,ny-1,1,nzsoil,1,nzsoil,        &
                amax,amin, imax,jmax,kmax, imin,jmin,kmin)
    IF (myproc == 0) WRITE(6,'(2(1x,a,g25.12,3(a,i4)))')                &
         'tsoilmin=',amin,' at i=',imin,', j=',jmin,', k=',kmin,        &
         'tsoilmax=',amax,' at i=',imax,', j=',jmax,', k=',kmax

    CALL a3dmax(qsoil,1,nx,1,nx-1,1,ny,1,ny-1,1,nzsoil,1,nzsoil,        &
                amax,amin, imax,jmax,kmax, imin,jmin,kmin)
    IF (myproc == 0) WRITE(6,'(2(1x,a,g25.12,3(a,i4)))')                &
         'qsoilmin=',amin,' at i=',imin,', j=',jmin,', k=',kmin,        &
         'qsoilmax=',amax,' at i=',imax,', j=',jmax,', k=',kmax

    CALL a3dmax(wetcanp,1,nx,1,nx-1,1,ny,1,ny-1,1,1,1,1,                &
                amax,amin, imax,jmax,kmax, imin,jmin,kmin)
    IF (myproc == 0) WRITE(6,'(2(1x,a,g25.12,2(a,i4)))')                &
         'wetcmin=',amin,' at i=',imin,', j=',jmin,                     &
         'wetcmax=',amax,' at i=',imax,', j=',jmax

  END IF

  IF (myproc == 0) THEN

    WRITE(6,*) ' '

    IF ( ncalls(mptr) == 0 ) THEN

      CALL get_output_dirname(0,dirname,curtim,1,outdirname,istatus)

      maxfn  = TRIM(outdirname)//runname(1:lfnkey)//'.maxmin'
      lmaxfn = LEN_TRIM(maxfn)

      IF(nestgrd == 1) THEN
        WRITE(maxfn((lmaxfn+1):(lmaxfn+4)),'(a,i2.2)')'.g',mptr
        lmaxfn =  lmaxfn + 4
      END IF

      WRITE(6,'(1x,a,a,a/,1x,a)')                                       &
          'Check to see if file ',maxfn(1:lmaxfn),' already exists.',   &
          'If so, append a version number to the filename.'

      CALL fnversn( maxfn, lmaxfn )

      CALL getunit ( nchmax1(mptr) )

      OPEN(nchmax1(mptr),FORM='formatted',STATUS='new',                 &
           FILE=maxfn(1:lmaxfn),IOSTAT=istat)

      IF( istat /= 0) THEN

        WRITE(6,'(/a,i2,/a/)')                                          &
            ' Error occured when opening file '//maxfn(1:lmaxfn)//      &
            ' using FORTRAN unit ',nchmax1(mptr),                       &
            ' Program stopped in MAXMIN.'
        CALL arpsstop("arpsstop called from MAXMIN with istat problem",1)

      END IF

      WRITE(nchmax1(mptr),'(a)') ''''//runname//''''

      WRITE(nchmax1(mptr),'(t2,a,t15,a,t25,a,t35,a,t45,                 &
      &   a,t55,a,t65,a,t75,a,                                          &
      &   t85,a,t95,a,t105,a,t115,a,t125,a,t135,                        &
      &   a,t145,a,t155,a,t165,                                         &
      &   a,t175,a,t185,a,t195,a,t204,a)')                              &
          '''TIME','UMIN','UMAX','VMIN','VMAX',                         &
          'WMIN','WMAX','PTMIN',                                        &
          'PTMAX','PMIN','PMAX','QVMAX','QCMAX',                        &
          'QRMAX','QIMAX',                                              &
          'QSMAX','QHMAX','UVMAX','VORLX','VORUX',                      &
          'REFMAX'//''''

      ncalls(mptr) = 1

    END IF

!    WRITE(nchmax1(mptr),'(f9.2,7f10.5,f10.5,2f10.2,f10.5)',ADVANCE='NO')&
!                        curtim,umin,umax,vmin,vmax,wmin,wmax,ptmin,     &
!                        ptmax,pmin,pmax,qvmax*1000
!    DO nq = 1,nscalar
!      WRITE(nchmax1(mptr),'(f10.5)',ADVANCE='NO') qscalarmax(nq)*1000.0
!    END DO
!    WRITE(nchmax1(mptr),'(4f10.5)') uvmax,vorlx,vorux,refmax

    qcmax = 0.0; qrmax = 0.0; qimax = 0.0; qsmax = 0.0; qhmax = 0.0

    IF (P_QC > 0) qcmax = qscalarmax(P_QC)
    IF (P_QR > 0) qrmax = qscalarmax(P_QR)
    IF (P_QI > 0) qimax = qscalarmax(P_QI)
    IF (P_QS > 0) qsmax = qscalarmax(P_QS)
    IF (P_QH > 0) THEN
       qhmax = qscalarmax(P_QH)
    ELSE IF (P_QG > 0) THEN
       qhmax = qscalarmax(P_QG)
    END IF

    WRITE(nchmax1(mptr),'(f9.2,7f10.5, f10.5,2f10.2,6f10.5,4f10.5)')    &
        curtim,umin,umax,vmin,vmax,wmin,wmax,ptmin,ptmax,               &
        pmin,pmax,qvmax*1000,qcmax*1000.0,qrmax*1000.0,                 &
        qimax*1000,qsmax*1000.0,qhmax*1000.0,                           &
        uvmax,vorlx,vorux,refmax

  END IF

  RETURN
END SUBROUTINE maxmin
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE BASPRT                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE basprt(nx,ny,nz,                                             &
           ubar,vbar,ptbar,pbar,rhobar,qvbar, zp,hterain)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Print out base state fields in FORTRAN unit nch=6.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  11/13/91.
!
!  MODIFICATION HISTORY:
!
!  6/06/92 (M. Xue)
!  Added full documentation.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!    ubar     Base state zonal velocity component (m/s)
!    vbar     Base state meridional velocity component (m/s)
!    ptbar    Base state potential temperature (K)
!    pbar     Base state pressure (Pascal)
!    rhobar   Base state density (kg/m**3)
!    qvbar    Base state water vapor specific humidity (kg/kg)
!    zp       Vertical coordinate of grid points in physical space (m)
!    hterain  Terrain height (m)
!
!  OUTPUT:
!
!    None.
!
!-----------------------------------------------------------------------
!

!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx,ny,nz          ! Number of grid points in 3 directions.

  REAL :: ubar  (nx,ny,nz)     ! Base state u-velocity (m/s)
  REAL :: vbar  (nx,ny,nz)     ! Base state v-velocity (m/s)
  REAL :: ptbar (nx,ny,nz)     ! Base state potential temperature (K)
  REAL :: pbar  (nx,ny,nz)     ! Base state pressure (Pascal)
  REAL :: rhobar(nx,ny,nz)     ! Base state air density (kg/m**3)
  REAL :: qvbar (nx,ny,nz)     ! Base state water vapor specific humidity
                               ! (kg/kg)
  REAL :: zp    (nx,ny,nz)     ! The physical height coordinate defined at
                               ! w-point of staggered grid.
  REAL :: hterain(nx,ny)       ! The height of the terrain.
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k,mode
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

  WRITE(6,'(/1x,a/)') 'Print out of base state fields.'
!
!-----------------------------------------------------------------------
!
!  Formatted printing of x-y slice at the model mid-level
!
!-----------------------------------------------------------------------
!

  mode = 1

  k = (nz-1)/2+1

  CALL wrigar(ubar,1,nx,1,ny,1,nz,1,nx,1,ny-1,k,k,                      &
              'ubar xy',0.0, mode )

  CALL wrigar(vbar,1,nx,1,ny,1,nz,1,nx-1,1,ny,k,k,                      &
              'vbar xy',0.0, mode )

  CALL wrigar(ptbar,1,nx,1,ny,1,nz,1,nx-1,1,ny-1,k,k,                   &
              'ptbar xy',0.00, mode )

  CALL wrigar(pbar,1,nx,1,ny,1,nz,1,nx-1,1,ny-1,k,k,                    &
              'pbar xy',0.0, mode )

  CALL wrigar(rhobar,1,nx,1,ny,1,nz,1,nx-1,1,ny-1,k,k,                  &
              'rhobar xy',0.00, mode )

  CALL wrigar(qvbar,1,nx,1,ny,1,nz,1,nx-1,1,ny-1,k,k,                   &
              'qvbar xy',0.0, mode )

!
!-----------------------------------------------------------------------
!
!  Formatted printing of x-z slice through the center of the model
!  domain.
!
!-----------------------------------------------------------------------
!

  j = (ny-1)/2+1

  CALL wrigar(zp,1,nx,1,ny,1,nz,1,nx-1,j,j,1,nz,'zp x-z',0.0, 2 )

  CALL wrigar(ubar,1,nx,1,ny,1,nz,1,nx,j,j,1,nz-1                       &
              ,'ubar x-z',0.0, 2 )

  CALL wrigar(vbar,1,nx,1,ny,1,nz,1,nx-1,j,j,1,nz-1                     &
              ,'vbar x-z',0.0, 2 )

  CALL wrigar(ptbar,1,nx,1,ny,1,nz,1,nx-1,j,j,1,nz-1                    &
              ,'ptbar x-z',0.00, 2 )

  CALL wrigar(pbar,1,nx,1,ny,1,nz,1,nx-1,j,j,1,nz-1                     &
              ,'pbar x-z',0.0, 2 )

  CALL wrigar(rhobar,1,nx,1,ny,1,nz,1,nx-1,j,j,1,nz-1                   &
              ,'rhobar x-z',0.0, 2 )

  CALL wrigar(qvbar,1,nx,1,ny,1,nz,1,nx-1,j,j,1,nz-1                    &
              ,'qvbar x-z',0.0, 2 )

!
!-----------------------------------------------------------------------
!
!  Formatted printing of y-z slice through the model domain center
!
!-----------------------------------------------------------------------
!

  i = (nx-1)/2+1

  CALL wrigar(ubar,1,nx,1,ny,1,nz,i,i,1,ny-1,1,nz-1                     &
              ,'ubar y-z',0.0, 3 )

  CALL wrigar(vbar,1,nx,1,ny,1,nz,i,i,1,ny  ,1,nz-1                     &
              ,'vbar y-z',0.0, 3 )


  CALL wrigar(ptbar,1,nx,1,ny,1,nz,i,i,1,ny-1,1,nz-1                    &
              ,'ptbar y-z',0.00, 3 )

  CALL wrigar(pbar,1,nx,1,ny,1,nz,i,i,1,ny-1,1,nz-1                     &
              ,'pbar y-z',0.0, 3 )

  CALL wrigar(rhobar,1,nx,1,ny,1,nz,i,i,1,ny-1,1,nz-1                   &
              ,'rhobar y-z',0.0, 3 )

  CALL wrigar(qvbar,1,nx,1,ny,1,nz,i,i,1,ny-1,1,nz-1                    &
              ,'qvbar y-z',0.0, 3 )

  RETURN
END SUBROUTINE basprt
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE FMTPRT                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE fmtprt(nx,ny,nz, tlevel,                                     &
                  u, v, w, ptprt, pprt, qv,qscalar,tke,kmh,kmv,         &
                  x,y,z,zp,hterain, j1,j2,j3, tem1)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Produce formatted print out of array tables in FORTRAN unit nch=6.
!  By default, the x-y, x-z and y-z slices through the domain center
!  are printed.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  5/29/91.
!
!  MODIFICATION HISTORY:
!
!  6/06/92 (M. Xue)
!  Added full documentation.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!    tlevel   The time level at which the data are printed.
!
!    u        x component of velocity (m/s)
!    v        y component of velocity (m/s)
!    w        Vertical component of Cartesian velocity (m/s)
!    ptprt    Perturbation potential temperature (K)
!    pprt     Perturbation pressure (Pascal)
!    qv       Water vapor specific humidity (kg/kg)
!    qc       Cloud water mixing ratio (kg/kg)
!    qr       Rainwater mixing ratio (kg/kg)
!    qi       Cloud ice mixing ratio (kg/kg)
!    qs       Snow mixing ratio (kg/kg)
!    qh       Hail mixing ratio (kg/kg)
!    tke      Turbulent Kinetic Energy ((m/s)**2)
!
!    kmh      Horizontal turb. mixing coef. for momentum ( m**2/s )
!    kmv      Vertical turb. mixing coef. for momentum ( m**2/s )
!
!    x        x coordinate of grid points in physical/comp. space (m)
!    y        y coordinate of grid points in physical/comp. space (m)
!    z        z coordinate of grid points in computational space (m)
!    zp       Vertical coordinate of grid points in physical space (m)
!    hterain  Terrain height (m)
!
!    j1       Coordinate transformation Jacobian -d(zp)/dx
!    j2       Coordinate transformation Jacobian -d(zp)/dy
!    j3       Coordinate transformation Jacobian  d(zp)/dz
!
!
!  OUTPUT:
!
!    None.
!
!  WORK ARRAYS:
!
!    tem1     Temporary work array.
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
  INCLUDE 'globcst.inc'
  INCLUDE 'timelvls.inc'
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  INTEGER :: nx,ny,nz          ! Number of grid points in 3 directions

  INTEGER :: tlevel            ! Time level at which data are printed.

  REAL :: u     (nx,ny,nz,nt)  ! Total u-velocity (m/s)
  REAL :: v     (nx,ny,nz,nt)  ! Total v-velocity (m/s)
  REAL :: w     (nx,ny,nz,nt)  ! Total w-velocity (m/s)
  REAL :: ptprt (nx,ny,nz,nt)  ! Perturbation potential temperature (K)
  REAL :: pprt  (nx,ny,nz,nt)  ! Perturbation pressure (Pascal)

  REAL :: qv    (nx,ny,nz,nt)  ! Water vapor specific humidity (kg/kg)

  REAL :: qscalar(nx,ny,nz,nt,nscalar)

  REAL :: tke   (nx,ny,nz,nt)  ! Turbulent Kinetic Energy ((m/s)**2)

  REAL :: kmh   (nx,ny,nz)     ! Horizontal turb. mixing coef. for
                               ! momentum. ( m**2/s )
  REAL :: kmv   (nx,ny,nz)     ! Vertical turb. mixing coef. for
                               ! momentum. ( m**2/s )

  REAL :: x     (nx)           ! The x-coord. of the physical and
                               ! computational grid. Defined at u-point.
  REAL :: y     (ny)           ! The y-coord. of the physical and
                               ! computational grid. Defined at v-point.
  REAL :: z     (nz)           ! The z-coord. of the computational grid.
                               ! Defined at w-point on the staggered grid.
  REAL :: zp    (nx,ny,nz)     ! The physical height coordinate defined at
                               ! w-point of the staggered grid.

  REAL :: hterain(nx,ny)       ! Terrain height.

  REAL :: j1    (nx,ny,nz)     ! Coordinate transformation Jacobian defined
                               ! as - d( zp )/d( x )
  REAL :: j2    (nx,ny,nz)     ! Coordinate transformation Jacobian defined
                               ! as - d( zp )/d( y )
  REAL :: j3    (nx,ny,nz)     ! Coordinate transformation Jacobian defined
                               ! as d( zp )/d( z )

  REAL :: tem1  (nx,ny,nz)     ! Temporary work array.
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i, j, k, nq
  INTEGER :: mode      ! Control the printing of an individual slice
                       ! mode = 1: print x-y slice
                       !      = 2: print x-z slice
                       !      = 3: print y-z slice
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
!  The following parameters are passed into this subroutine through
!  a common block in globcst.inc. They control the output of
!  water and ice variables.
!
!  mstout =0 or 1. If mstout=0, water variables are not printed.
!  iceout =0 or 1. If iceout=0, printing of qi, qs and qh is skipped.
!  basout =0 or 1. If basout=0, base state variables are not dumped.
!
!-----------------------------------------------------------------------
!

!
!-----------------------------------------------------------------------
!
!  Formatted printing of x-y slice at the model mid-levels.
!
!-----------------------------------------------------------------------
!
!  return

  mode = 1
  k = (nz-1)/2+1

  CALL wrigar(u(1,1,1,tlevel),1,nx,1,ny,1,nz,1,nx,1,ny-1,k,k,           &
              'u x-y',0.0, mode )

  CALL wrigar(v(1,1,1,tlevel),1,nx,1,ny,1,nz,1,nx-1,1,ny,k,k,           &
              'v x-y',0.0, mode )

  CALL wrigar(w(1,1,1,tlevel),1,nx,1,ny,1,nz,1,nx-1,1,ny-1,k+1,k+1,     &
              'w x-y',0.0, mode )

  CALL wrigar(ptprt(1,1,1,tlevel),1,nx,1,ny,1,nz,1,nx-1,1,ny-1,k,k,     &
              'ptprt xy',0.00, mode )

  CALL wrigar(pprt(1,1,1,tlevel),1,nx,1,ny,1,nz,1,nx-1,1,ny-1,k,k,      &
              'pprt xy',0.0, mode )

  IF( mstout == 1) THEN
  CALL wrigar(qv(1,1,1,tlevel),1,nx,1,ny,1,nz,1,nx-1,1,ny-1,k,k,      &
              'qv x-y',0.0, mode )

  DO nq = 1,nscalar

    CALL wrigar(qscalar(:,:,:,tlevel,nq),1,nx,1,ny,1,nz,1,nx-1,1,ny-1,k,k,      &
                TRIM(qnames(nq))//' x-y',0.0, mode )

  END DO
  END IF

  CALL wrigar(kmh,1,nx,1,ny,1,nz,1,nx-1,1,ny-1,k,k,'kmhx-y',.0,mode)
  CALL wrigar(kmv,1,nx,1,ny,1,nz,1,nx-1,1,ny-1,k,k,'kmvx-y',.0,mode)
!
!-----------------------------------------------------------------------
!
!  Formatted printing of x-z slice through the model domain center.
!
!-----------------------------------------------------------------------
!
  mode = 2         ! Print x-z slices

  j = (ny-1)/2+1

  CALL wrigar(u(1,1,1,tlevel),1,nx,1,ny,1,nz,1,nx,j,j,1,nz-1            &
               ,'u x-z',0.0, mode )

  CALL wrigar(v(1,1,1,tlevel),1,nx,1,ny,1,nz,1,nx-1,j,j,1,nz-1          &
               ,'v x-z',0.0, mode )

  CALL wrigar(w(1,1,1,tlevel),1,nx,1,ny,1,nz,1,nx-1,j,j,1,nz            &
               ,'w x-z',0.0, mode )

  CALL wrigar(ptprt(1,1,1,tlevel),1,nx,1,ny,1,nz,1,nx-1,j,j,1,nz-1      &
               ,'ptprt x-z',0.00, mode )

  CALL wrigar(pprt(1,1,1,tlevel),1,nx,1,ny,1,nz,1,nx-1,j,j,1,nz-1       &
               ,'pprt x-z',0.0, mode )

  IF( mstout == 1) THEN
  CALL wrigar(qv(1,1,1,tlevel),1,nx,1,ny,1,nz,1,nx-1,j,j,1,nz-1         &
               ,'qv x-z',0.0, mode )

  DO nq = 1,nscalar

    CALL wrigar(qscalar(:,:,:,tlevel,nq),1,nx,1,ny,1,nz,1,nx-1,j,j,     &
                1,nz-1,TRIM(qnames(nq))//' x-z',0.0, mode )

  END DO
  END IF

  CALL wrigar(kmh ,1,nx,1,ny,1,nz,1,nx-1,j,j,1,nz-1                     &
               ,'kmhx-z',0.0, mode )
  CALL wrigar(kmv ,1,nx,1,ny,1,nz,1,nx-1,j,j,1,nz-1                     &
               ,'kmvx-z',0.0, mode )
!
!
!-----------------------------------------------------------------------
!
!  Formatted printing of y-z slice through the model domain center.
!
!-----------------------------------------------------------------------
!
  mode = 3         ! Print y-z slices

  i = (nx-1)/2+1

  CALL wrigar(u(1,1,1,tlevel),1,nx,1,ny,1,nz,i,i,1,ny-1,1,nz-1          &
              ,'u y-z',0.0, mode )

  CALL wrigar(v(1,1,1,tlevel),1,nx,1,ny,1,nz,i,i,1,ny  ,1,nz-1          &
              ,'v y-z',0.0, mode )

  CALL wrigar(w(1,1,1,tlevel),1,nx,1,ny,1,nz,i,i,1,ny-1,1,nz            &
              ,'w y-z',0.0, mode )

  CALL wrigar(ptprt(1,1,1,tlevel),1,nx,1,ny,1,nz,i,i,1,ny-1,1,nz-1      &
              ,'ptprt y-z',0.00, mode )

  CALL wrigar(pprt(1,1,1,tlevel),1,nx,1,ny,1,nz,i,i,1,ny-1,1,nz-1       &
              ,'pprt y-z',0.0, mode )

  IF( mstout == 1) THEN
  CALL wrigar(qv(1,1,1,tlevel),1,nx,1,ny,1,nz,i,i,1,ny-1,1,nz-1         &
              ,'qv y-z',0.0, mode )

  DO nq = 1,nscalar

    CALL wrigar(qscalar(:,:,:,tlevel,nq),1,nx,1,ny,1,nz,i,i,1,ny-1,     &
                1,nz-1,TRIM(qnames(nq))//' y-z',0.0, mode )

  END DO
  END IF

  CALL wrigar(kmh ,1,nx,1,ny,1,nz,i,i,1,ny-1,1,nz-1                     &
              ,'kmhy-z',0.0, mode )
  CALL wrigar(kmv ,1,nx,1,ny,1,nz,i,i,1,ny-1,1,nz-1                     &
              ,'kmvy-z',0.0, mode )


  RETURN
END SUBROUTINE fmtprt

!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE ABORTDMP                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE abortdmp(mptr,nx,ny,nz,nzsoil,nstyps,                        &
           u,v,w,wcont,ptprt,pprt,qv,qscalar,tke,                       &
           ubar,vbar,ptbar,pbar,rhobar,qvbar,kmh,kmv,                   &
           x,y,z,zp,hterain,zpsoil,mapfct,j1,j2,j3,j3soil,              &
           soiltyp,stypfrct,vegtyp,lai,roufns,veg,                      &
           tsoil,qsoil,wetcanp,snowdpth,                                &
           raing,rainc,prcrate,                                         &
           radfrc,radsw,rnflx,radswnet,radlwin,                         &
           usflx,vsflx,ptsflx,qvsflx,                                   &
           tem1,tem2,tem3, tem4)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Write out a number of files including the history and
!  restart dumps when the model aborts.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  8/8/95.
!
!  MODIFICATION HISTORY:
!
!  Based on CHKSTAB.
!
!  05/13/2002  (J. Brotzge)
!  Added additional soil scheme arrays, soil layers
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    mptr     Grid identifier.
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!    nzsoil   Number of grid points in the soil
!
!    u        x component of velocity (m/s)
!    v        y component of velocity (m/s)
!    w        Vertical component of Cartesian velocity (m/s)
!    wcont    Contravariant vertical velocity (m/s)
!    ptprt    Perturbation potential temperature (K)
!    pprt     Perturbation pressure (Pascal)
!    qv       Water vapor specific humidity (kg/kg)
!    qc       Cloud water mixing ratio (kg/kg)
!    qr       Rainwater mixing ratio (kg/kg)
!    qi       Cloud ice mixing ratio (kg/kg)
!    qs       Snow mixing ratio (kg/kg)
!    qh       Hail mixing ratio (kg/kg)
!    tke      Turbulent Kinetic Energy ((m/s)**2)
!
!    ubar     Base state zonal velocity component (m/s)
!    vbar     Base state meridional velocity component (m/s)
!    ptbar    Base state potential temperature (K)
!    pbar     Base state pressure (Pascal)
!    rhobar   Base state density (kg/m**3)
!    qvbar    Base state water vapor specific humidity (kg/kg)
!
!    kmh      Horizontal turb. mixing coef. for momentum ( m**2/s )
!    kmv      Vertical turb. mixing coef. for momentum ( m**2/s )
!
!    x        x coordinate of grid points in physical/comp. space (m)
!    y        y coordinate of grid points in physical/comp. space (m)
!    z        z coordinate of grid points in computational space (m)
!    zp       Vertical coordinate of grid points in physical space (m)
!    zpsoil   Vertical coordinate of grid points in the soil (m)
!    hterain  Terrain height (m)
!    mapfct   Map factors at scalar, u and v points
!    j1       Coordinate transformation Jacobian -d(zp)/dx
!    j2       Coordinate transformation Jacobian -d(zp)/dy
!    j3       Coordinate transformation Jacobian  d(zp)/dz
!    j3soil   Coordinate transformation Jacobian  d(zpsoil)/dz
!
!    soiltyp  Soil type
!    stypfrct  Soil type fraction
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
!    prcrate  Precipitation rates
!
!    radfrc   Radiation forcing (K/s)
!    radsw    Solar radiation reaching the surface
!    rnflx    Net radiation flux absorbed by surface
!    radswnet Net solar radiation
!    radlwin  Incoming longwave radiation
!
!    usflx    Surface flux of u-momentum (kg/(m*s**2))
!    vsflx    Surface flux of v-momentum (kg/(m*s**2))
!    ptsflx   Surface heat flux (K*kg/(m**2 * s ))
!    qvsflx   Surface moisture flux of (kg/(m**2 * s))
!
!  OUTPUT:
!
!    None.
!
!  WORK ARRAYS:
!
!    tem1     Temporary work array.
!    tem2     Temporary work array.
!    tem3     Temporary work array.
!    tem4     Temporary work array.
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
  INCLUDE 'phycst.inc'
  INCLUDE 'globcst.inc'
  INCLUDE 'mp.inc'            ! Message passing parameters.

  INCLUDE 'timelvls.inc'
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  INTEGER :: mptr              ! Grid identifier.

  INTEGER :: nx,ny,nz          ! Number of grid points in 3 directions
  INTEGER :: nzsoil            ! Number of grid points in the soil

  REAL :: u     (nx,ny,nz,nt)  ! Total u-velocity (m/s)
  REAL :: v     (nx,ny,nz,nt)  ! Total v-velocity (m/s)
  REAL :: w     (nx,ny,nz,nt)  ! Total w-velocity (m/s)
  REAL :: wcont (nx,ny,nz)     ! Contravariant vertical velocity (m/s)
  REAL :: ptprt (nx,ny,nz,nt)  ! Perturbation potential temperature (K)
  REAL :: pprt  (nx,ny,nz,nt)  ! Perturbation pressure (Pascal)

  REAL :: qv    (nx,ny,nz,nt)  ! Water vapor specific humidity (kg/kg)
  REAL :: qscalar(nx,ny,nz,nt,nscalar)
  REAL :: tke   (nx,ny,nz,nt)  ! turbulent kinetic energy

  REAL :: ubar  (nx,ny,nz)     ! Base state u-velocity (m/s)
  REAL :: vbar  (nx,ny,nz)     ! Base state v-velocity (m/s)
  REAL :: ptbar (nx,ny,nz)     ! Base state potential temperature (K)
  REAL :: pbar  (nx,ny,nz)     ! Base state pressure (Pascal)
  REAL :: rhobar(nx,ny,nz)     ! Base state air density (kg/m**3)
  REAL :: qvbar (nx,ny,nz)     ! Base state water vapor specific humidity
                               ! (kg/kg)

  REAL :: kmh   (nx,ny,nz)     ! Horizontal turb. mixing coef. for
                               ! momentum. ( m**2/s )
  REAL :: kmv   (nx,ny,nz)     ! Vertical turb. mixing coef. for
                               ! momentum. ( m**2/s )

  REAL :: x     (nx)           ! The x-coord. of the physical and
                               ! computational grid. Defined at u-point.
  REAL :: y     (ny)           ! The y-coord. of the physical and
                               ! computational grid. Defined at v-point.
  REAL :: z     (nz)           ! The z-coord. of the computational grid.
                               ! Defined at w-point on the staggered grid.
  REAL :: zp    (nx,ny,nz)     ! The physical height coordinate defined at
                               ! w-point of the staggered grid.
  REAL :: zpsoil (nx,ny,nzsoil)      ! The physical height coordinate defined at
                               ! w-point of the soil layers.

  REAL :: hterain(nx,ny)       ! Terrain height.
  REAL :: mapfct(nx,ny,8)      ! Map factors at scalar, u and v points

  REAL :: j1    (nx,ny,nz)     ! Coordinate transformation Jacobian defined as
                               ! - d( zp )/d( x )
  REAL :: j2    (nx,ny,nz)     ! Coordinate transformation Jacobian defined as
                               ! - d( zp )/d( y )
  REAL :: j3    (nx,ny,nz)     ! Coordinate transformation Jacobian defined as
                               ! d( zp )/d( z )
  REAL :: j3soil (nx,ny,nzsoil)      ! Coordinate transformation Jacobian defined as
                               ! d( zpsoil)/d( z )

  INTEGER :: nstyps                    ! Number of soil types
  INTEGER :: soiltyp (nx,ny,nstyps)    ! Soil type
  REAL    :: stypfrct(nx,ny,nstyps)    ! Soil type fraction
  INTEGER :: vegtyp (nx,ny)            ! Vegetation type
  REAL    :: lai    (nx,ny)            ! Leaf Area Index
  REAL    :: roufns (nx,ny)            ! Surface roughness
  REAL    :: veg    (nx,ny)            ! Vegetation fraction

  REAL :: tsoil  (nx,ny,nzsoil,0:nstyps) ! Soil temperature (K)
  REAL :: qsoil  (nx,ny,nzsoil,0:nstyps) ! Soil moisture (m**3/m**3)
  REAL :: wetcanp(nx,ny,0:nstyps)      ! Canopy water amount
  REAL :: snowdpth(nx,ny)              ! Snow depth (m)

  REAL :: raing(nx,ny)         ! Grid supersaturation rain
  REAL :: rainc(nx,ny)         ! Cumulus convective rain
  REAL :: prcrate(nx,ny,4)     ! precipitation rate (kg/(m**2*s))
                               ! prcrate(1,1,1) = total precip. rate
                               ! prcrate(1,1,2) = grid scale precip. rate
                               ! prcrate(1,1,3) = cumulus precip. rate
                               ! prcrate(1,1,4) = microphysics precip. rate

  REAL :: radfrc(nx,ny,nz)     ! Radiation forcing (K/s)
  REAL :: radsw (nx,ny)        ! Solar radiation reaching the surface
  REAL :: rnflx (nx,ny)        ! Net radiation flux absorbed by surface
  REAL :: radswnet (nx,ny)     ! Net solar radiation, SWin - SWout
  REAL :: radlwin  (nx,ny)     ! Incoming longwave radiation


  REAL :: usflx (nx,ny)        ! Surface flux of u-momentum (kg/(m*s**2))
  REAL :: vsflx (nx,ny)        ! Surface flux of v-momentum (kg/(m*s**2))
  REAL :: ptsflx(nx,ny)        ! Surface heat flux (K*kg/(m*s**2))
  REAL :: qvsflx(nx,ny)        ! Surface moisture flux (kg/(m**2*s))

  REAL :: tem1  (nx,ny,nz)     ! Temporary work array
  REAL :: tem2  (nx,ny,nz)     ! Temporary work array
  REAL :: tem3  (nx,ny,nz)     ! Temporary work array
  REAL :: tem4  (nx,ny,nz)     ! Temporary work array
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: ierr
  INTEGER :: tlevel, tim
  INTEGER :: grdbas,i

  CHARACTER(LEN=256) :: outdirname
  INTEGER            :: oldirnam
  INTEGER :: istatus
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  mgrid = mptr
  tlevel = tpresent
!
!-----------------------------------------------------------------------
!
!  Max./min. statistics calculation and printing of initial fields
!
!-----------------------------------------------------------------------
!
  CALL maxmin(mptr,nx,ny,nz,nzsoil,tlevel,rhobar,                       &
              u,v,w,wcont,ptprt,pprt,qv,qscalar,tke,kmh,kmv,            &
              x,y,z,zp,zpsoil,mapfct,                                   &
              tsoil(1,1,1,0),qsoil(1,1,1,0),wetcanp(1,1,0),             &
              tem1,tem2,tem3,tem4)
!
!-----------------------------------------------------------------------
!
!  Formatted printing of time dependent variables
!
!-----------------------------------------------------------------------
!
  CALL fmtprt(nx,ny,nz, tlevel,                                         &
              u,v,w,ptprt,pprt,qv,qscalar,tke,kmh,kmv,                  &
              x,y,z,zp,hterain, j1,j2,j3, tem1)
!
!-----------------------------------------------------------------------
!
!  Check to see if any history data dump is to be produced.
!
!-----------------------------------------------------------------------
!
  IF( hdmpfmt == 0 ) THEN

    IF (myproc == 0) WRITE(6,'(1x,a)')  &
         'History data dump option was 0, no data is dumped.'
    GO TO 800

  END IF

!
!-----------------------------------------------------------------------
!
!  Find a unique name hdmpfn(1:ldmpf) for the history dump data set
!  at time 'curtim'.
!
!  For the savi3D data dump case, the file name is specified only
!  once in INITOUT.
!
!-----------------------------------------------------------------------
!
  CALL get_output_dirname(1,dirname,curtim,1,outdirname,istatus)
  oldirnam = LEN_TRIM(dirname)

  IF( hdmpfmt /= 5 .AND. hdmpfmt /= 9 ) THEN

    CALL gtdmpfn(runname(1:lfnkey),outdirname,                          &
                 oldirnam,curtim,hdmpfmt,                               &
                 mgrid,nestgrd, hdmpfn, ldmpf)

  END IF

  IF (myproc == 0) WRITE(6,'(1x,a,a)')  &
       'History data dump in file ',hdmpfn(1:ldmpf)
  grdbas = 0      ! No base state or grid array is dumped.

  tim = tpresent

!   blocking inserted for ordering i/o for message passing
  DO i=0,nprocs-1,dumpstride
    IF(myproc >= i.AND.myproc <= i+dumpstride-1)THEN

      CALL dtadump(nx,ny,nz,nzsoil,nstyps,                              &
               hdmpfmt,nchdmp,hdmpfn(1:ldmpf),                          &
               grdbas,filcmprs,                                         &
               u(1,1,1,tim),v(1,1,1,tim),w(1,1,1,tim),                  &
               ptprt(1,1,1,tim),pprt(1,1,1,tim),qv(1,1,1,tim),          &
               qscalar(:,:,:,tim,:),tke(1,1,1,tim),kmh,kmv,             &
               ubar,vbar,tem1,ptbar,pbar,rhobar,qvbar,                  &
               x,y,z,zp,zpsoil,                                         &
               soiltyp,stypfrct,vegtyp,lai,roufns,veg,                  &
               tsoil,qsoil,wetcanp,snowdpth,                            &
               raing,rainc,prcrate,                                     &
               radfrc,radsw,rnflx,radswnet,radlwin,                     &
               usflx,vsflx,ptsflx,qvsflx,                               &
               tem2,tem3,tem4)

    END IF
    IF (mp_opt > 0) CALL mpbarrier
  END DO

!
!-----------------------------------------------------------------------
!
!  If ARPS is dumping Savi3D GRAF format history files,
!  call grafclose to close the GRAF file before program stop.
!
!-----------------------------------------------------------------------
!

  IF (hdmpfmt == 5) THEN
    CALL mclosescheme (gridid, ierr)
    CALL mclosedataset (dsindex, ierr)
  END IF

  IF (hdmpfmt == 9) CLOSE (nchdmp)

  800   CONTINUE

  RETURN
END SUBROUTINE abortdmp

!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE WRTXYSLIC                  ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE wrtxyslic(nx,ny,nz, u,v,w,ptprt,pprt,qv,qscalar,             &
           ubar,vbar,ptbar,pbar,qvbar, x,y,zp, zslice,                  &
           fnkey,time,zs,tem1)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Write out 2-D x-y slice of a 3-D array for given k or height z.
!
!  We assume the grid levels are flat, and base state variables
!  are constant on each level.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!    u        x component of velocity at a given time level (m/s)
!    v        y component of velocity at a given time level (m/s)
!    w        Vertical component of Cartesian velocity at a given
!             time level (m/s)
!    ptprt    Perturbation potential temperature at a given time
!             level (K)
!    pprt     Perturbation pressure at  a given time level (Pascal)
!    qv       Water vapor specific humidity (kg/kg)
!    qc       Cloud water mixing ratio at a given time level (kg/kg)
!    qr       Rainwater mixing ratio at a given time level (kg/kg)
!
!    ubar     Base state x-velocity component (m/s)
!    vbar     Base state y-velocity component (m/s)
!    ptbar    Base state potential temperature (K)
!    pbar     Base state pressure (Pascal)
!    qvbar    Base state water vapor specific humidity (kg/kg)
!
!    x        x coordinate of grid points in physical/comp. space (m)
!    y        y coordinate of grid points in physical/comp. space (m)
!    zp       z coordinate of grid points in physcial space (m)
!
!    zslice   The height of the x-y slice to be written out.
!
!    fnkey    File name key
!    time     time of data in seconds
!
!    mgrid    The grid number
!    nestgrd  Flag for nested grid run.
!
!  WORK ARRAYS:
!
!    zs       Temporary work array, zp averaged to scalar point.
!    tem1     Temporary work array.
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------

  INCLUDE 'globcst.inc'
  INCLUDE 'mp.inc'            ! Message passing parameters.

!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------

  INTEGER :: nx,ny,nz

  REAL :: u     (nx,ny,nz)     ! Total u-velocity (m/s)
  REAL :: v     (nx,ny,nz)     ! Total v-velocity (m/s)
  REAL :: w     (nx,ny,nz)     ! Total w-velocity (m/s)
  REAL :: ptprt (nx,ny,nz)     ! Perturbation potential temperature (K)
  REAL :: pprt  (nx,ny,nz)     ! Perturbation pressure (Pascal)

  REAL :: qv    (nx,ny,nz)     ! Water vapor specific humidity (kg/kg)
  REAL :: qscalar(nx,ny,nz,nscalar)

  REAL :: ubar  (nx,ny,nz)     ! Base state x-velocity (m/s)
  REAL :: vbar  (nx,ny,nz)     ! Base state y-velocity (m/s)
  REAL :: ptbar (nx,ny,nz)     ! Base state potential temperature (K)
  REAL :: pbar  (nx,ny,nz)     ! Base state pressure (Pascal)
  REAL :: qvbar (nx,ny,nz)     ! Base state Water vapor specific
                               ! humidity (kg/kg)

  REAL :: x     (nx)           ! The x-coord. of the physical and
                               ! computational grid. Defined at u-point.
  REAL :: y     (ny)           ! The y-coord. of the physical and
                               ! computational grid. Defined at v-point.
  REAL :: zp    (nx,ny,nz)     ! The physical height coordinate defined at
                               ! w-point of the staggered grid.

  REAL :: zslice               ! The height of the x-y slice to be written out.

  CHARACTER (LEN=80 ) :: timsnd
  CHARACTER (LEN=*  ) :: fnkey ! The height of the x-y slice to be written out
  INTEGER :: tmstrln           ! length of timsnd
  REAL    :: time              ! time of data in seconds

  INTEGER :: kslice

  REAL :: zs  (nx,ny,nz)       ! Temporary work array, zp averaged
                               ! to scalar point
  REAL :: tem1(nx,ny)

!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------

  INTEGER :: i,j,k,length,nunit,istat,ierr
  CHARACTER (LEN=256) :: xysfn
  CHARACTER (LEN=35 ) :: gridnum
  DATA gridnum /'123456789abcdefghijklmnopqrstuvwxyz'/

  CHARACTER (LEN=256) :: savename
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  CALL cvttsnd( time, timsnd, tmstrln )

  xysfn = fnkey
  length = LEN( fnkey )
  WRITE(xysfn(length+1:length+7),'(a,i5.5)') '.z',nint(zslice)
  xysfn(length+8:length+8+tmstrln) = 't'//timsnd(1:tmstrln)
  length = length + 8 + tmstrln

  IF( nestgrd == 1 ) THEN
    WRITE(xysfn(length+1:length+3),'(''.G'',A)')                        &
          gridnum(mgrid:mgrid)
    length = length + 3
  END IF

  CALL getunit( nunit )

  IF (mp_opt > 0) THEN
    savename(1:256) = xysfn(1:256)
    CALL gtsplitfn(savename,1,1,loc_x,loc_y,1,1,0,0,0,2,xysfn,istat)
    length = LEN_TRIM(xysfn)
  END IF

  CALL fnversn(xysfn, length )

  CALL asnctl ('NEWLOCAL', 1, ierr)
  CALL asnfile(xysfn(1:length), '-F f77 -N ieee', ierr)

  OPEN(UNIT=nunit,FILE=trim(xysfn(1:length)),STATUS='unknown',          &
       FORM='unformatted',IOSTAT= istat )

  IF (mp_opt > 0) THEN
    xysfn(1:256) = savename(1:256)
    length = length - 5
  END IF

  CALL avgz(zp, 0 ,                                                     &
            nx,ny,nz, 1,nx, 1,ny, 1,nz-1, zs)

  kslice = nz-2
  DO k=2,nz-1
    IF( zs(1,1,k) > zslice ) THEN
      kslice = k-1
      EXIT
    END IF
  END DO

  k = kslice

  WRITE(nunit) fnkey
  WRITE(nunit) nx,ny
  WRITE(nunit) time
  WRITE(nunit) zslice

  WRITE(nunit) 'x    '
  WRITE(nunit) x

  WRITE(nunit) 'y    '
  WRITE(nunit) y

  DO j=1,ny
    DO i=1,nx
      tem1(i,j)=u(i,j,k)+ (u(i,j,k+1)-u(i,j,k))*                        &
          (zslice-zs(i,j,k))/(zs(i,j,k+1)-zs(i,j,k))
    END DO
  END DO

  WRITE(nunit) 'u    '
  WRITE(nunit) ubar(1,1,k)+ (ubar(1,1,k+1)-ubar(1,1,k))*                &
               (zslice-zs(1,1,k))/(zs(1,1,k+1)-zs(1,1,k))
  WRITE(nunit) tem1

  DO j=1,ny
    DO i=1,nx
      tem1(i,j)=v(i,j,k)+ (v(i,j,k+1)-v(i,j,k))*                        &
          (zslice-zs(i,j,k))/(zs(i,j,k+1)-zs(i,j,k))
    END DO
  END DO

  WRITE(nunit) 'v    '
  WRITE(nunit) vbar(1,1,k)+ (vbar(1,1,k+1)-vbar(1,1,k))*                &
               (zslice-zs(1,1,k))/(zs(1,1,k+1)-zs(1,1,k))
  WRITE(nunit) tem1

  DO j=1,ny
    DO i=1,nx
      tem1(i,j)=0.5*(w(i,j,k)+w(i,j,k+1))+                              &
                0.5*(w(i,j,k+2)-w(i,j,k))*                              &
                (zslice-zs(i,j,k))/(zs(i,j,k+1)-zs(i,j,k))
    END DO
  END DO

  WRITE(nunit) 'w    '
  WRITE(nunit) 0.0
  WRITE(nunit) tem1

  DO j=1,ny
    DO i=1,nx
      tem1(i,j)=ptprt(i,j,k)+ (ptprt(i,j,k+1)-ptprt(i,j,k))*            &
          (zslice-zs(i,j,k))/(zs(i,j,k+1)-zs(i,j,k))
    END DO
  END DO

  WRITE(nunit) 'ptprt'
  WRITE(nunit) ptbar(1,1,k)+ (ptbar(1,1,k+1)-ptbar(1,1,k))*             &
               (zslice-zs(1,1,k))/(zs(1,1,k+1)-zs(1,1,k))
  WRITE(nunit) tem1

  DO j=1,ny
    DO i=1,nx
      tem1(i,j)=pprt(i,j,k)+ (pprt(i,j,k+1)-pprt(i,j,k))*               &
          (zslice-zs(i,j,k))/(zs(i,j,k+1)-zs(i,j,k))
    END DO
  END DO

  WRITE(nunit) 'pprt '
  WRITE(nunit) pbar(1,1,k)+ (pbar(1,1,k+1)-pbar(1,1,k))*                &
               (zslice-zs(1,1,k))/(zs(1,1,k+1)-zs(1,1,k))
  WRITE(nunit) tem1

  DO j=1,ny
    DO i=1,nx
      tem1(i,j)= qv(i,j,k)+ (qv(i,j,k+1)-qv(i,j,k))*                    &
          (zslice-zs(i,j,k))/(zs(i,j,k+1)-zs(i,j,k))
    END DO
  END DO

  WRITE(nunit) 'qv   '
  WRITE(nunit) qvbar(1,1,k)+ (qvbar(1,1,k+1)-qvbar(1,1,k))*             &
               (zslice-zs(1,1,k))/(zs(1,1,k+1)-zs(1,1,k))
  WRITE(nunit) tem1

  IF (P_QC > 0) THEN
    DO j=1,ny
      DO i=1,nx
        tem1(i,j)= qscalar(i,j,k,P_QC)                                  &
                 + (qscalar(i,j,k+1,P_QC)-qscalar(i,j,k,P_QC))*         &
                   (zslice-zs(i,j,k))/(zs(i,j,k+1)-zs(i,j,k))
      END DO
    END DO

    WRITE(nunit) 'qc   '
    WRITE(nunit) 0.0
    WRITE(nunit) tem1
  END IF

  IF (P_QR > 0) THEN
    DO j=1,ny
      DO i=1,nx
        tem1(i,j)= qscalar(i,j,k,P_QR)                                  &
                 + (qscalar(i,j,k+1,P_QR)-qscalar(i,j,k,P_QR))*         &
                   (zslice-zs(i,j,k))/(zs(i,j,k+1)-zs(i,j,k))
      END DO
    END DO

    WRITE(nunit) 'qr   '
    WRITE(nunit) 0.0
    WRITE(nunit) tem1
  END IF

  CLOSE(UNIT = nunit )
  CALL retunit( nunit )

  RETURN
END SUBROUTINE wrtxyslic
