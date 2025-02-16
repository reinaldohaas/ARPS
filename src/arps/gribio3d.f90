!
!##################################################################
!##################################################################
!######                                                      ######
!######               SUBROUTINE GRIBDUMP                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE gribdump(nx,ny,nz,nzsoil,nstyps, nchanl,filnam, grdbas,      &
           u,v,w,ptprt,pprt,qv,qscalar,tke,kmh,kmv,                     &
           ubar,vbar,ptbar,pbar,rhobar,qvbar,                           &
           x,y,z,zp,zpsoil,                                             &
           soiltyp,stypfrct,vegtyp,lai,roufns,veg,                      &
           tsoil,qsoil,wetcanp,snowdpth,                                &
           raing,rainc,prcrate,                                         &
           radfrc,radsw,rnflx,radswnet,radlwin,                         &
           usflx,vsflx,ptsflx,qvsflx,                                   &
           tem1,temxy1,temxy2)
!           tem1,item1,temxy1,temxy2,itemxy1,itemxy2)    ! original
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Write history data into channel nchanl as GRIB data.
!
!  All data read in are located at the original staggered grid points.
!
!  Note: coordinate fields are dumped as 3 dimensional fields which
!  have been converted from meters to kilometers.  This is for the
!  convenience of the plotting applications.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yuhe Liu
!  11/01/1995
!
!  MODIFICATION HISTORY:
!
!  12/09/1998 (Donghai Wang)
!  Added the snow cover.
!
!  1/07/2000 (Eric Kemp)
!  Snow cover is now written as GRIB variable 66 (snow depth (m))
!
!  05/23/2002 (J. Brotzge)
!  Added soil variables to allow for multiple soil schemes.  
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
!    nchanl   FORTRAN I/O channel number for history data output.
!    filnam   File name to store the GRIB messages
!    grdbas   Flag indicating if this is a call for the data dump
!             of grid and base state arrays only. If so, grdbas=1.
!
!    u        x component of velocity at a given time level (m/s)
!    v        y component of velocity at a given time level (m/s)
!    w        Vertical component of Cartesian velocity at a given
!             time level (m/s)
!    ptprt    Perturbation potential temperature at a given time
!             level (K)
!    pprt     Perturbation pressure at a given time level (Pascal)
!    qv       Water vapor specific humidity at a given time level (kg/kg)
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
!    zp       Vertical coordinate of grid points in physical space (m)
!    zpsoil   Vertical coordinate of grid points in the soil (m)
!
!    soiltyp  Soil type
!    stypfrct  Soil type fraction
!    vegtyp   Vegetation type
!    lai      Leaf Area Index
!    roufns   Surface roughness
!    veg      Vegetation fraction
!
!    tsoil    Soil temperature (K) 
!    qsoil    Soil moisture
!    wetcanp  Canopy water amount
!
!    raing    Grid supersaturation rain
!    rainc    Cumulus convective rain
!    prcrate  Precipitation rates
!
!    radfrc   Radiation forcing (K/s)
!    radsw    Solar radiation reaching the surface
!    rnflx    Net radiation flux absorbed by surface
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
!  WORK ARRAY:
!
!    tem1     Temporary work array.
!    item1    Integer temporary work array.
!    temxy1   2-D temporary work array.
!    temxy2   2-D temporary work array.
!    itemxy1  2-D integer temporary work array.
!    itemxy2  2-D integer temporary work array.
!
!
!-----------------------------------------------------------------------
!
!  The following parameters are passed into this subroutine through
!  a common block in globcst.inc, and they determine which
!  variables are output.
!
!  grdout =0 or 1. If grdout =0, grid variables are not dumped.
!  basout =0 or 1. If basout =0, base state variables are not dumped.
!  varout =0 or 1. If varout =0, model variables are not dumped.
!  mstout =0 or 1. If mstout =0, water variables are not dumped.
!  rainout=0 or 1. If rainout=0, rain variables are not dumped.
!  prcout =0 or 1. If prcout =0, precipitation rates are not dumped.
!  iceout =0 or 1. If iceout =0, qi, qs and qh are not dumped.
!  tkeout =0 or 1. If tkeout =0, tke is not dumped.
!  trbout =0 or 1. If trbout =0, kmh and kmv are not dumped.
!  sfcout =0 or 1. If sfcout =0, surface variables are not dumped.
!  landout=0 or 1. If landout=0, surface propertty arrays are not dumped.
!  radout =0 or 1. If radout =0, precipitation rates are not dumped.
!  flxout =0 or 1. If flxout =0, precipitation rates are not dumped.
!
!  These following parameters are also passed in through common
!  blocks in globcst.inc.
!
!  runname,curtim,umove,vmove,xgrdorg,ygrdorg
!
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'mp.inc'          ! Message passing parameters.

!  include 'gribcst.inc'
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  INTEGER :: nx,ny,nz        ! Number of grid points in 3 directions
  INTEGER :: nzsoil          ! Number of grid points in the soil.  

  INTEGER :: nchanl            ! FORTRAN I/O channel number for output
  CHARACTER (LEN=*     ) :: filnam ! File name
  INTEGER :: grdbas            ! If this is a grid/base state array dump

  REAL :: u     (nx,ny,nz)     ! Total u-velocity (m/s)
  REAL :: v     (nx,ny,nz)     ! Total v-velocity (m/s)
  REAL :: w     (nx,ny,nz)     ! Total w-velocity (m/s)
  REAL :: ptprt (nx,ny,nz)     ! Perturbation potential temperature (K)
  REAL :: pprt  (nx,ny,nz)     ! Perturbation pressure (Pascal)

  REAL :: qv    (nx,ny,nz)     ! Water vapor specific humidity (kg/kg)

  REAL :: qscalar    (nx,ny,nz,nscalar)

  REAL :: tke   (nx,ny,nz)     ! Turbulent Kinetic Energy ((m/s)**2)

  REAL :: kmh   (nx,ny,nz)     ! Horizontal turb. mixing coef. for
                               ! momentum. ( m**2/s )
  REAL :: kmv   (nx,ny,nz)     ! Vertical turb. mixing coef. for
                               ! momentum. ( m**2/s )

  REAL :: ubar  (nx,ny,nz)     ! Base state u-velocity (m/s)
  REAL :: vbar  (nx,ny,nz)     ! Base state v-velocity (m/s)
  REAL :: ptbar (nx,ny,nz)     ! Base state potential temperature (K)
  REAL :: pbar  (nx,ny,nz)     ! Base state pressure (Pascal)
  REAL :: rhobar(nx,ny,nz)     ! Base state air density (kg/m**3)
  REAL :: qvbar (nx,ny,nz)     ! Base state water vapor specific humidity
                               ! (kg/kg)

  REAL :: x     (nx)           ! The x-coord. of the physical and
                               ! computational grid. Defined at u-point.
  REAL :: y     (ny)           ! The y-coord. of the physical and
                               ! computational grid. Defined at v-point.
  REAL :: z     (nz)           ! The z-coord. of the computational grid.
                               ! Defined at w-point on the staggered grid.
  REAL :: zp    (nx,ny,nz)     ! The physical height coordinate defined at
                               ! w-point of the staggered grid.
  REAL :: zpsoil (nx,ny,nzsoil) ! The physical height coordinate defined at
                                ! w-point of the soil.  

  INTEGER :: nstyps                  ! Number of soil types
  INTEGER :: soiltyp(nx,ny,nstyps)   ! Soil type
  REAL :: stypfrct(nx,ny,nstyps)     ! Soil type fraction
  INTEGER :: vegtyp (nx,ny)          ! Vegetation type
  REAL :: lai    (nx,ny)          ! Leaf Area Index
  REAL :: roufns (nx,ny)          ! Surface roughness
  REAL :: veg    (nx,ny)          ! Vegetation fraction

  REAL :: tsoil (nx,ny,nzsoil,0:nstyps) ! Soil temperature (K)
  REAL :: qsoil (nx,ny,nzsoil,0:nstyps) ! Soil moisture
  REAL :: wetcanp(nx,ny,0:nstyps)    ! Canopy water amount
  REAL :: snowdpth(nx,ny)            ! Snow depth (m)

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
  REAL :: radswnet(nx,ny)      ! Net shortwave radiation
  REAL :: radlwin(nx,ny)       ! Incoming longwave radiation

  REAL :: usflx (nx,ny)        ! Surface flux of u-momentum (kg/(m*s**2))
  REAL :: vsflx (nx,ny)        ! Surface flux of v-momentum (kg/(m*s**2))
  REAL :: ptsflx(nx,ny)        ! Surface heat flux (K*kg/(m*s**2))
  REAL :: qvsflx(nx,ny)        ! Surface moisture flux (kg/(m**2*s))

  REAL :: tem1(nx,ny,nz)       ! Temporary work array
  REAL :: temxy1(nx,ny)        ! temporary work array
  REAL :: temxy2(nx,ny)        ! temporary work array

!
!-----------------------------------------------------------------------
!
!  Parameters for GRIB packing
!
!-----------------------------------------------------------------------
!
  INTEGER :: nbufsz
  PARAMETER   ( nbufsz = 800000 )  ! Size of GRIB buffer

  INTEGER :: ipds(25)          ! PDS integer array
  INTEGER :: igdss(25)         ! GDS integer array for s-grid
  INTEGER :: igdsu(25)         ! GDS integer array for u-grid
  INTEGER :: igdsv(25)         ! GDS integer array for v-grid

  INTEGER :: ibdshd(4)         ! BDS header

  CHARACTER (LEN=1) :: pds(28)       ! PDS ( GRIB Section 1)
  CHARACTER (LEN=1) :: gds(42)       ! GDS ( GRIB Section 2)
  CHARACTER (LEN=1) :: bds(nbufsz)   ! BDS ( GRIB Section 4)
  INTEGER :: ibds(nbufsz/4)    ! identical to BDS

  EQUIVALENCE ( bds,ibds )

  CHARACTER (LEN=1) :: mgrib(nbufsz) ! Buffer to carry GRIB messages

  INTEGER :: nbytes            ! Byte position to start read or write
  INTEGER :: wdlen             ! Length of machine word
!
!-----------------------------------------------------------------------
!
!  Parameters describing this routine
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=40) :: fmtver
  PARAMETER (fmtver='005.30 ARPS Data Format 10')
!
!-----------------------------------------------------------------------
!
!  Variables to determine the length of machine word, wdlen
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=8) :: ctema,ctemb

  INTEGER :: itema,itemb

  EQUIVALENCE ( itema,ctema )
  EQUIVALENCE ( itemb,ctemb )
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k,is
  INTEGER :: istatus
  INTEGER :: itype     ! Array type: 0 - floating, 1 - integer
  INTEGER :: dflag     ! Dimension flag: 0 - 3-D, 1 - 2-D

  INTEGER :: nq
  INTEGER :: nqipds(10)
  INTEGER, ALLOCATABLE :: item1  (:,:,:)  ! Integer temporary work array
  INTEGER, ALLOCATABLE :: itemxy1(:,:)    ! Integer temporary work array
  INTEGER, ALLOCATABLE :: itemxy2(:,:)    ! Integer temporary work array
!
!-----------------------------------------------------------------------
!
!  Save ipds & igds for future use
!
!-----------------------------------------------------------------------
!
  LOGICAL :: firstcall(0:mgrdmax)

  SAVE firstcall, wdlen
  DATA firstcall/ mgrdmax*.true.,.true. /
  DATA ctema / '12345678' /
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

  WRITE(6,'(/2a/,2(a/,2a//),2(a/))') 'WARNING: Using of GRIB format',   &
     ' is discouraged for ARPS history files. this is because',         &
     '         If grbpkbit = 16, the data loses precision badly.',      &
     '         If you prefer to GRIB format, grbpkbit = 32 is ',        &
     'recommended.', &
     '         ARPS supported data format are binary, HDF4 and NetCDF', &
     '         GrADS and GRIB format should be used with GrADS',        &
     ' for display purpose only.',                                      &
     'To use GRIB format, please go to file "src/arps/gribio3d.f90" ',  &
     'and comment out the first "arpsstop" call in subroutine gribdump.'
  CALL arpsstop('GRIB format not support. Please check output file for detail.',1)

  ALLOCATE(item1(nx,ny,nz), STAT = istatus)
  ALLOCATE(itemxy1(nx,ny),  STAT = istatus)
  ALLOCATE(itemxy2(nx,ny),  STAT = istatus)

  item1(:,:,:) = 0
  itemxy1(:,:) = 0
  itemxy2(:,:) = 0

!
!-----------------------------------------------------------------------
!
!   Initialize the integer array for BDS header
!
!-----------------------------------------------------------------------
!
  ibdshd(1) = 0
  ibdshd(2) = 0
  ibdshd(3) = 0
  ibdshd(4) = 0
!
!-----------------------------------------------------------------------
!
!  Initialize the integer array, ipds, for PDS.
!
!-----------------------------------------------------------------------
!
  CALL mkipds( nx,ny,nz, ipds )
!
!-----------------------------------------------------------------------
!
!  Initialize the GDS integer arrays for different grid type, igdsu,
!  igdsv, and igdss
!
!-----------------------------------------------------------------------
!
  CALL mkigds( nx,ny,nz, 0, igdss )
  CALL mkigds( nx,ny,nz, 1, igdsu )
  CALL mkigds( nx,ny,nz, 2, igdsv )

  IF (mgrid < 0 .or. mgrid > mgrdmax) THEN
    mgrid = 0
  ENDIF

  IF ( firstcall(mgrid) ) THEN
!
!-----------------------------------------------------------------------
!
!  Calculate the length of machine word, wdlen, in bytes, which will
!  be used to pack the data. Assume the length has only two possible
!  value: 4 and 8
!
!-----------------------------------------------------------------------
!
    itemb = itema
    IF ( ctema /= ctemb ) THEN
      wdlen = 4
    ELSE
      wdlen = 8
    END IF
    IF (myproc == 0) WRITE (6,'(a,i2,a)')  &
        'The length of machine word is ',wdlen,' bytes'
!
!-----------------------------------------------------------------------
!
!  Check if number of bits per datum fits size of computer word.
!
!-----------------------------------------------------------------------
!
    IF ( wdlen*8 < grbpkbit ) THEN
      IF (myproc == 0) WRITE (6,'(a,i3,a,i3/a,a/a)')  &
          'Number of bits per datum, ',grbpkbit,                        &
          ', exceeds the machine word length, ',wdlen*8,                &
          'Please change grbpkbit in input namelist &output ', &
          'to fit the machine word length.',                            &
          'Program stoped in GRIBDUMP'
      CALL arpsstop(" ",1)
    END IF
!
!-----------------------------------------------------------------------
!
!  Write the sample of GrADS control file used to display the ARPS
!  GRIB files. This control file is to display the variables on the
!  scalar points. The GRIB grid type was given by ipds(5)
!
!-----------------------------------------------------------------------
!
    CALL gribcntl( nx,ny,nz,nzsoil,x,y,z,zp,zpsoil,ipds(5), temxy1,temxy2 )

    firstcall(mgrid) = .false.

  END IF

  IF (myproc == 0) WRITE(6,'(1x,a,f13.3/)') 'Writing history data at time=', curtim
!
!-----------------------------------------------------------------------
!
!  Write the header of ARPS GRIB file
!
!-----------------------------------------------------------------------
!
  CALL wrthishd(nx,ny,nz,nzsoil,nchanl,grdbas,fmtver, x,y,z,             &
                wdlen,nbufsz,mgrib,nbytes)
!
!-----------------------------------------------------------------------
!
!  If grdout=1 or grdbas=1, write out grid variables
!
!-----------------------------------------------------------------------
!
  IF(grdout == 1 .OR. grdbas == 1 ) THEN

    itype = 0       ! floating array
    ibdshd(3) = itype

    dflag = 0       ! 3-D array
    ipds(8)  = 8    ! geometric height (m)
    ipds(9)  = 103  ! z-coorinates
    ipds(18) = 0    ! forecast time = 0 for grid and base variables

    IF (myproc == 0) WRITE(6,'(a)') 'Writing zp'
    CALL edgfill(zp,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz)
    CALL wrtgrib(nx,ny,nz,1,z,zp,item1,                                 &
                  nchanl,nbytes,wdlen,itype,dflag,                      &
                  ipds,igdss,ibdshd,pds,gds,                            &
                  nbufsz,bds,ibds,                                      &
                  mgrib,temxy1,temxy2,itemxy1,itemxy2)

    IF (myproc == 0) WRITE(6,'(a)') 'Writing zpsoil'
    CALL edgfill(zpsoil,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nzsoil,1,nzsoil)
    DO k = 1,nzsoil
      CALL wrtgrib(nx,ny,1,1,z,zpsoil(1,1,k),item1,                     &
                  nchanl,nbytes,wdlen,itype,dflag,                      &
                  ipds,igdss,ibdshd,pds,gds,                            &
                  nbufsz,bds,ibds,                                      &
                  mgrib,temxy1,temxy2,itemxy1,itemxy2)
    END DO
  END IF    ! grdout
!
!-----------------------------------------------------------------------
!
!  If basout=1, write out base state variables.
!
!-----------------------------------------------------------------------
!
  IF(basout == 1 .OR. grdbas == 1 ) THEN

    itype = 0       ! floating array
    ibdshd(3) = itype

    dflag = 0       ! 3-D array
    ipds(9)  = 103  ! z-coorinates
    ipds(18) = 0    ! forecast time = 0 for grid and base variables

    ipds(8)  = 33   ! u wind (m/s)
    IF (myproc == 0) WRITE(6,'(a,f10.5)') 'Writing ubar ',ubar(1,1,1)
    CALL edgfill(ubar,1,nx,1,nx, 1,ny,1,ny-1, 1,nz,1,nz-1)
    CALL wrtgrib(nx,ny,nz,0,z,ubar,item1,                               &
                  nchanl,nbytes,wdlen,itype,dflag,                      &
                  ipds,igdsu,ibdshd,pds,gds,                            &
                  nbufsz,bds,ibds,                                      &
                  mgrib,temxy1,temxy2,itemxy1,itemxy2)

    ipds(8)  = 34   ! v wind (m/s)
    IF (myproc == 0) WRITE(6,'(a,f10.5)') 'Writing vbar ',vbar(1,1,1)
    CALL edgfill(vbar,1,nx,1,nx-1, 1,ny,1,ny, 1,nz,1,nz-1)
    CALL wrtgrib(nx,ny,nz,0,z,vbar,item1,                               &
                  nchanl,nbytes,wdlen,itype,dflag,                      &
                  ipds,igdsv,ibdshd,pds,gds,                            &
                  nbufsz,bds,ibds,                                      &
                  mgrib,temxy1,temxy2,itemxy1,itemxy2)


    ipds(8)  = 13   ! potential temperature (K)
    IF (myproc == 0) WRITE(6,'(a,f10.2)') 'Writing ptbar ',ptbar(1,1,1)
    CALL edgfill(ptbar,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)
    CALL wrtgrib(nx,ny,nz,0,z,ptbar,item1,                              &
                  nchanl,nbytes,wdlen,itype,dflag,                      &
                  ipds,igdss,ibdshd,pds,gds,                            &
                  nbufsz,bds,ibds,                                      &
                  mgrib,temxy1,temxy2,itemxy1,itemxy2)

    ipds(8)  = 1    ! pressure (Pascal)
    IF (myproc == 0) WRITE(6,'(a,f10.2)') 'Writing pbar ', pbar(1,1,1)
    CALL edgfill(pbar,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)
    CALL wrtgrib(nx,ny,nz,0,z,pbar,item1,                               &
                  nchanl,nbytes,wdlen,itype,dflag,                      &
                  ipds,igdss,ibdshd,pds,gds,                            &
                  nbufsz,bds,ibds,                                      &
                  mgrib,temxy1,temxy2,itemxy1,itemxy2)

    IF(mstout == 1) THEN

      ipds(8)  = 51   ! specific humidity (kg/kg)
      IF (myproc == 0) WRITE(6,'(a)') 'Writing qvbar'
      CALL edgfill(qvbar,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)
      CALL wrtgrib(nx,ny,nz,0,z,qvbar,item1,                            &
                    nchanl,nbytes,wdlen,itype,dflag,                    &
                    ipds,igdss,ibdshd,pds,gds,                          &
                    nbufsz,bds,ibds,                                    &
                    mgrib,temxy1,temxy2,itemxy1,itemxy2)

    END IF

    IF(landout == 1) THEN

      dflag = 1         ! 2-D array
      ipds(9)  = 111    ! depth below land surface

      IF( nstyp <= 1 ) THEN

        itype = 1       ! integer array
        ibdshd(3) = itype

        ipds(11) = 0    ! 0 centimeter for surface

        ipds(8)  = 224  ! soil type
        IF (myproc == 0) WRITE(6,'(a)') 'Writing soiltyp '

        CALL iedgfill(soiltyp(1,1,1),1,nx,1,nx-1, 1,ny,1,ny-1,          &
                     1,1,1,1)
        CALL wrtgrib(nx,ny,1,0,z,tem1(1,1,1),soiltyp(1,1,1),            &
                     nchanl,nbytes,wdlen,itype,dflag,                   &
                     ipds,igdss,ibdshd,pds,gds,                         &
                     nbufsz,bds,ibds,                                   &
                     mgrib,temxy1,temxy2,itemxy1,itemxy2)

      ELSE

        itype = 1       ! integer array
        ibdshd(3) = itype
        ipds(8)  = 224  ! soil type

        IF (myproc == 0) WRITE(6,'(a)') 'Writing soiltyp '

        DO is=1,nstyp
          ipds(11) = is ! 0 centimeter for surface

          CALL iedgfill(soiltyp(1,1,is),1,nx,1,nx-1, 1,ny,1,ny-1,       &
                        1,1,1,1)
          CALL wrtgrib(nx,ny,1,0,z,tem1(1,1,1),soiltyp(1,1,is),         &
                       nchanl,nbytes,wdlen,itype,dflag,                 &
                       ipds,igdss,ibdshd,pds,gds,                       &
                       nbufsz,bds,ibds,                                 &
                       mgrib,temxy1,temxy2,itemxy1,itemxy2)
        END DO

        itype = 0       ! floating array
        ibdshd(3) = itype
        ipds(8)  = 226  ! redefine 226 for fraction of soil types

        IF (myproc == 0) WRITE(6,'(a)') 'Writing stypfrct'

        DO is=1,nstyp
          ipds(11) = is ! 0 centimeter for surface

          CALL edgfill(stypfrct(1,1,is),1,nx,1,nx-1, 1,ny,1,ny-1,       &
                       1,1,1,1)
          CALL wrtgrib(nx,ny,1,0,z,stypfrct(1,1,is),item1(1,1,1),       &
                       nchanl,nbytes,wdlen,itype,dflag,                 &
                       ipds,igdss,ibdshd,pds,gds,                       &
                       nbufsz,bds,ibds,                                 &
                       mgrib,temxy1,temxy2,itemxy1,itemxy2)
        END DO

      END IF

      itype = 1       ! integer array
      ibdshd(3) = itype

      ipds(11) = 0
      ipds(8)  = 225  ! vegtyp
      IF (myproc == 0) WRITE(6,'(a)') 'Writing vegtyp'
      CALL iedgfill(vegtyp ,1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)
      CALL wrtgrib(nx,ny,1,0,z,tem1(1,1,1),vegtyp,                      &
                    nchanl,nbytes,wdlen,itype,dflag,                    &
                    ipds,igdss,ibdshd,pds,gds,                          &
                    nbufsz,bds,ibds,                                    &
                    mgrib,temxy1,temxy2,itemxy1,itemxy2)

      itype = 0       ! floating array
      ibdshd(3) = itype

      ipds(8)  = 227  ! redefine 227 for LAI
      IF (myproc == 0) WRITE(6,'(a)') 'Writing lai'
      CALL edgfill(lai    ,1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)
      CALL wrtgrib(nx,ny,1,0,z,lai,item1(1,1,1),                        &
                    nchanl,nbytes,wdlen,itype,dflag,                    &
                    ipds,igdss,ibdshd,pds,gds,                          &
                    nbufsz,bds,ibds,                                    &
                    mgrib,temxy1,temxy2,itemxy1,itemxy2)

      ipds(8)  = 83   ! surface roughness (m)
      IF (myproc == 0) WRITE(6,'(a)') 'Writing roufns'
      CALL edgfill(roufns ,1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)
      CALL wrtgrib(nx,ny,1,0,z,roufns,item1(1,1,1),                     &
                    nchanl,nbytes,wdlen,itype,dflag,                    &
                    ipds,igdss,ibdshd,pds,gds,                          &
                    nbufsz,bds,ibds,                                    &
                    mgrib,temxy1,temxy2,itemxy1,itemxy2)

      ipds(8)  = 87   ! vegetation (%)
      IF (myproc == 0) WRITE(6,'(a)') 'Writing veg'
      CALL edgfill(veg    ,1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)
      DO j = 1,ny
        DO i = 1,nx
          tem1(i,j,1) = veg(i,j) * 100.
        END DO
      END DO
      CALL wrtgrib(nx,ny,1,0,z,tem1(1,1,1),item1(1,1,1),                &
                    nchanl,nbytes,wdlen,itype,dflag,                    &
                    ipds,igdss,ibdshd,pds,gds,                          &
                    nbufsz,bds,ibds,                                    &
                    mgrib,temxy1,temxy2,itemxy1,itemxy2)

    END IF

  END IF

  IF ( grdbas == 1 ) RETURN

  ipds(18) = nint(curtim/60) ! forecast time in minutes
!
!-----------------------------------------------------------------------
!
!  If varout = 1, Write out u, v, w, pt, and p.
!
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!
!  Write out u,v and w.
!
!-----------------------------------------------------------------------
!
  IF(varout == 1) THEN

    itype = 0                  ! floating array
    ibdshd(3) = itype

    dflag = 0                  ! 3-D array
    ipds(9)  = 103             ! z-coorinates

    ipds(8)  = 33   ! u wind (m/s)
    IF (myproc == 0) WRITE(6,'(a)') 'Writing u'
    CALL edgfill(u,1,nx,1,nx, 1,ny,1,ny-1, 1,nz,1,nz-1)
    CALL wrtgrib(nx,ny,nz,0,z,u,item1,                                  &
                  nchanl,nbytes,wdlen,itype,dflag,                      &
                  ipds,igdsu,ibdshd,pds,gds,                            &
                  nbufsz,bds,ibds,                                      &
                  mgrib,temxy1,temxy2,itemxy1,itemxy2)

    ipds(8)  = 34   ! v wind (m/s)
    IF (myproc == 0) WRITE(6,'(a)') 'Writing v'
    CALL edgfill(v,1,nx,1,nx-1, 1,ny,1,ny, 1,nz,1,nz-1)
    CALL wrtgrib(nx,ny,nz,0,z,v,item1,                                  &
                  nchanl,nbytes,wdlen,itype,dflag,                      &
                  ipds,igdsv,ibdshd,pds,gds,                            &
                  nbufsz,bds,ibds,                                      &
                  mgrib,temxy1,temxy2,itemxy1,itemxy2)

    ipds(8)  = 40   ! w wind (m/s)
    IF (myproc == 0) WRITE(6,'(a)') 'Writing w'
    CALL edgfill(w,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz)
    CALL wrtgrib(nx,ny,nz,1,z,w,item1,                                  &
                  nchanl,nbytes,wdlen,itype,dflag,                      &
                  ipds,igdss,ibdshd,pds,gds,                            &
                  nbufsz,bds,ibds,                                      &
                  mgrib,temxy1,temxy2,itemxy1,itemxy2)
!
!-----------------------------------------------------------------------
!
!  Write out scalars
!
!-----------------------------------------------------------------------
!
    ipds(8)  = 13   ! potential temperature (K)
    DO k = 1,nz-1
      DO j = 1,ny-1
        DO i = 1,nx-1
          tem1(i,j,k) = ptbar(i,j,k) + ptprt(i,j,k)
        END DO
      END DO
    END DO
    IF (myproc == 0) WRITE(6,'(a,f10.2)') 'Writing pt ', tem1(1,1,1)
    CALL edgfill(tem1,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)
    CALL wrtgrib(nx,ny,nz,0,z,tem1,item1,                               &
                  nchanl,nbytes,wdlen,itype,dflag,                      &
                  ipds,igdss,ibdshd,pds,gds,                            &
                  nbufsz,bds,ibds,                                      &
                  mgrib,temxy1,temxy2,itemxy1,itemxy2)

    ipds(8)  = 1    ! pressure (Pascal)
    DO k = 1,nz-1
      DO j = 1,ny-1
        DO i = 1,nx-1
          tem1(i,j,k) = pbar(i,j,k) + pprt(i,j,k)
        END DO
      END DO
    END DO
    CALL edgfill(tem1,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)
    IF (myproc == 0) WRITE(6,'(a,f10.2)') 'Writing p ', tem1(1,1,1)
    CALL wrtgrib(nx,ny,nz,0,z,tem1,item1,                               &
                  nchanl,nbytes,wdlen,itype,dflag,                      &
                  ipds,igdss,ibdshd,pds,gds,                            &
                  nbufsz,bds,ibds,                                      &
                  mgrib,temxy1,temxy2,itemxy1,itemxy2)

  END IF     ! varout
!
!-----------------------------------------------------------------------
!
!  If mstout = 1, write out moisture scalars.
!
!-----------------------------------------------------------------------
!
  IF(mstout == 1) THEN

    itype = 0                  ! floating array
    ibdshd(3) = itype

    dflag = 0                  ! 3-D array
    ipds(9)  = 103             ! z-coorinates

    ipds(8)  = 51   ! specific humidity (kg/kg)
    IF (myproc == 0) WRITE(6,'(a)') 'Writing qv'
    CALL edgfill(qv,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)
    CALL wrtgrib(nx,ny,nz,0,z,qv,item1,                                 &
                  nchanl,nbytes,wdlen,itype,dflag,                      &
                  ipds,igdss,ibdshd,pds,gds,                            &
                  nbufsz,bds,ibds,                                      &
                  mgrib,temxy1,temxy2,itemxy1,itemxy2)

    nqipds(:) = 255
    IF (P_QC > 0) nqipds(P_QC) = 153   ! cloud water mixing ratio (kg/kg)
    IF (P_QR > 0) nqipds(P_QR) = 170   ! rain water mixing ratio (kg/kg)
    IF (P_QI > 0) nqipds(P_QI) = 178   ! ice mixing ratio (kg/kg)
    IF (P_QS > 0) nqipds(P_QS) = 171   ! snow mixing ratio (kg/kg)
    IF (P_QH > 0) nqipds(P_QH) = 179   ! graupel hail mixing ratio (kg/kg)

    DO nq = 1,nscalar
      ipds(8)  = nqipds(nq)  
      IF (myproc == 0) WRITE(6,'(2a)') 'Writing ',TRIM(qnames(nq))
      CALL edgfill(qscalar(:,:,:,nq),1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)
      CALL wrtgrib(nx,ny,nz,0,z,qscalar(:,:,:,nq),item1,                &
                   nchanl,nbytes,wdlen,itype,dflag,                     &
                   ipds,igdss,ibdshd,pds,gds,                           &
                   nbufsz,bds,ibds,                                     &
                   mgrib,temxy1,temxy2,itemxy1,itemxy2)
    END DO

    IF(rainout == 1) THEN

      dflag = 1       ! 2-D
      ipds(9)  = 1    ! at surface

      ipds(8)  = 62   ! large scale precipitation for raing (kg/m**2)
      ipds(11) = 0    ! at surface
      IF (myproc == 0) WRITE(6,'(a)') 'Writing raing'
      CALL edgfill(raing,   1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)
      CALL wrtgrib(nx,ny,1,0,z,raing,item1,                             &
                    nchanl,nbytes,wdlen,itype,dflag,                    &
                    ipds,igdss,ibdshd,pds,gds,                          &
                    nbufsz,bds,ibds,                                    &
                    mgrib,temxy1,temxy2,itemxy1,itemxy2)

      ipds(8)  = 63   ! convective precipitation for rainc (kg/m**2)
      ipds(11) = 0    ! at surface
      IF (myproc == 0) WRITE(6,'(a)') 'Writing rainc'
      CALL edgfill(rainc,1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)
      CALL wrtgrib(nx,ny,1,0,z,rainc,item1,                             &
                    nchanl,nbytes,wdlen,itype,dflag,                    &
                    ipds,igdss,ibdshd,pds,gds,                          &
                    nbufsz,bds,ibds,                                    &
                    mgrib,temxy1,temxy2,itemxy1,itemxy2)

    END IF   !rainout

    IF ( prcout == 1 ) THEN

      CALL edgfill(prcrate, 1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)

      dflag = 1       ! 2-D
      ipds(8)  = 59   ! precipitation rates (kg/m**2/s)
      ipds(9)  = 111  ! use 111 to save more than one precip. rate

      ipds(11) = 1    ! total precipitation rate
      IF (myproc == 0) WRITE(6,'(a)') 'Writing prcrate1'
      CALL wrtgrib(nx,ny,1,0,z,prcrate(1,1,1),item1,                    &
                    nchanl,nbytes,wdlen,itype,dflag,                    &
                    ipds,igdss,ibdshd,pds,gds,                          &
                    nbufsz,bds,ibds,                                    &
                    mgrib,temxy1,temxy2,itemxy1,itemxy2)

      ipds(11)  = 2   ! grid scale precipitation rates
      IF (myproc == 0) WRITE(6,'(a)') 'Writing prcrate2'
      CALL wrtgrib(nx,ny,1,0,z,prcrate(1,1,2),item1,                    &
                    nchanl,nbytes,wdlen,itype,dflag,                    &
                    ipds,igdss,ibdshd,pds,gds,                          &
                    nbufsz,bds,ibds,                                    &
                    mgrib,temxy1,temxy2,itemxy1,itemxy2)

      ipds(11) = 3   ! cumulus precipitation rates
      IF (myproc == 0) WRITE(6,'(a)') 'Writing prcrate3'
      CALL wrtgrib(nx,ny,1,0,z,prcrate(1,1,3),item1,                    &
                    nchanl,nbytes,wdlen,itype,dflag,                    &
                    ipds,igdss,ibdshd,pds,gds,                          &
                    nbufsz,bds,ibds,                                    &
                    mgrib,temxy1,temxy2,itemxy1,itemxy2)

      ipds(11) = 4   ! microphysics precipitation rates
      IF (myproc == 0) WRITE(6,'(a)') 'Writing prcrate4'
      CALL wrtgrib(nx,ny,1,0,z,prcrate(1,1,4),item1,                    &
                    nchanl,nbytes,wdlen,itype,dflag,                    &
                    ipds,igdss,ibdshd,pds,gds,                          &
                    nbufsz,bds,ibds,                                    &
                    mgrib,temxy1,temxy2,itemxy1,itemxy2)

    END IF   ! prcout

  END IF   !mstout

  IF( tkeout == 1 ) THEN
    itype = 0                  ! floating array
    ibdshd(3) = itype

    dflag = 0       ! 3-D
    ipds(9)  = 103  ! z-coorinates
    ipds(8)  = 158  ! tke (J kg**-1)
    CALL edgfill(tke,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)
    CALL wrtgrib(nx,ny,nz,0,z,tke,item1,                                &
                  nchanl,nbytes,wdlen,itype,dflag,                      &
                  ipds,igdss,ibdshd,pds,gds,                            &
                  nbufsz,bds,ibds,                                      &
                  mgrib,temxy1,temxy2,itemxy1,itemxy2)

  END IF
!
!-----------------------------------------------------------------------
!
!  If trbout = 1, write out the turbulence parameter, km.
!
!-----------------------------------------------------------------------
!
  IF( trbout == 1 ) THEN
    itype = 0                  ! floating array
    ibdshd(3) = itype

    dflag = 0       ! 3-D
    ipds(9)  = 103  ! z-coorinates
    ipds(8)  = 185  ! redefine 185 for kmh
    IF (myproc == 0) WRITE(6,'(a)') 'Writing kmh'
    CALL edgfill(kmh,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)
    CALL wrtgrib(nx,ny,nz,0,z,kmh,item1,                                &
                  nchanl,nbytes,wdlen,itype,dflag,                      &
                  ipds,igdss,ibdshd,pds,gds,                            &
                  nbufsz,bds,ibds,                                      &
                  mgrib,temxy1,temxy2,itemxy1,itemxy2)

    dflag = 0       ! 3-D
    ipds(9)  = 103  ! z-coorinates
    ipds(8)  = 186  ! redefine 186 for kmv
    CALL edgfill(kmv,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)
    CALL wrtgrib(nx,ny,nz,0,z,kmv,item1,                                &
                  nchanl,nbytes,wdlen,itype,dflag,                      &
                  ipds,igdss,ibdshd,pds,gds,                            &
                  nbufsz,bds,ibds,                                      &
                  mgrib,temxy1,temxy2,itemxy1,itemxy2)

  END IF   ! trbout

!
!-----------------------------------------------------------------------
!
!  If sfcout = 1, write out the surface variables,
!  tsfc, tsoil, wetsfc, wetdp, and wetcanp.
!
!-----------------------------------------------------------------------
!
  IF ( sfcout == 1 ) THEN
    itype = 0                  ! floating array
    ibdshd(3) = itype

    dflag = 1         ! 2-D
    ipds(9)  = 111    ! below land surface

    IF( nstyp <= 1 ) THEN

      ipds(8)  = 85   ! soil temperature (K)
!      ipds(11) = 1    ! 1 cm below the surface
      ipds(11) = 0    ! At the surface
      IF (myproc == 0) WRITE(6,'(a)') 'Writing tsoil'
      CALL edgfill(tsoil(1,1,1,0),1,nx,1,nx-1,1,ny,1,ny-1,1,nzsoil,     &
                   1,nzsoil)
!      CALL wrtgrib(nx,ny,nzsoil,0,zpsoil,tsoil(1,1,1,0),item1,          &
!                   nchanl,nbytes,wdlen,itype,dflag,                     &
!                   ipds,igdss,ibdshd,pds,gds,                           &
!                   nbufsz,bds,ibds,                                     &
!                   mgrib,temxy1,temxy2,itemxy1,itemxy2)
      CALL wrtgrib(nx,ny,nzsoil,1,zpsoil,tsoil(1,1,1,0),item1,          &
                   nchanl,nbytes,wdlen,itype,dflag,                     &
                   ipds,igdss,ibdshd,pds,gds,                           &
                   nbufsz,bds,ibds,                                     &
                   mgrib,temxy1,temxy2,itemxy1,itemxy2)

      ipds(8)  = 144   ! Soil moisture (m**3/m**3)
!      ipds(11) = 100  ! 100 cm below the surface
      ipds(11) = 0    ! At the surface
      IF (myproc == 0) WRITE(6,'(a)') 'Writing qsoil'
      CALL edgfill(qsoil(1,1,1,0),1,nx,1,nx-1,1,ny,1,ny-1,1,nzsoil,1,   &
                     nzsoil)        
!      CALL wrtgrib(nx,ny,nzsoil,0,zpsoil,qsoil(1,1,1,0),item1,          &
!                   nchanl,nbytes,wdlen,itype,dflag,                     &
!                   ipds,igdss,ibdshd,pds,gds,                           &
!                   nbufsz,bds,ibds,                                     &
!                   mgrib,temxy1,temxy2,itemxy1,itemxy2)
      CALL wrtgrib(nx,ny,nzsoil,1,zpsoil,qsoil(1,1,1,0),item1,          &
                   nchanl,nbytes,wdlen,itype,dflag,                     &
                   ipds,igdss,ibdshd,pds,gds,                           &
                   nbufsz,bds,ibds,                                     &
                   mgrib,temxy1,temxy2,itemxy1,itemxy2)


      ipds(8)  = 223  ! amount of liquid water retained on canopy (m)
      ipds(11) = 0    ! 0 cm below the surface
      IF (myproc == 0) WRITE(6,'(a)') 'Writing wetcnpy'
      CALL edgfill(wetcanp(1,1,0),1,nx,1,nx-1,1,ny,1,ny-1,1,1,1,1)
      CALL wrtgrib(nx,ny,1,0,z,wetcanp(1,1,0),item1,                    &
                   nchanl,nbytes,wdlen,itype,dflag,                     &
                   ipds,igdss,ibdshd,pds,gds,                           &
                   nbufsz,bds,ibds,                                     &
                   mgrib,temxy1,temxy2,itemxy1,itemxy2)

    ELSE

      DO is=0,nstyp
        ipds(8)  = 85   ! soil temperature (K)
!        ipds(11) = 1+is ! 1 cm below the surface
        ipds(11) = is    ! 'is' cm below the surface
        IF (myproc == 0) WRITE(6,'(a)') 'Writing tsoil'
        CALL edgfill(tsoil(1,1,1,is),1,nx,1,nx-1,1,ny,1,ny-1,           &
                     1,nzsoil,1,nzsoil)
!        CALL wrtgrib(nx,ny,nzsoil,0,zpsoil,tsoil(1,1,1,is),item1,       &
!                     nchanl,nbytes,wdlen,itype,dflag,                   &
!                     ipds,igdss,ibdshd,pds,gds,                         &
!                     nbufsz,bds,ibds,                                   &
!                     mgrib,temxy1,temxy2,itemxy1,itemxy2)
        CALL wrtgrib(nx,ny,nzsoil,1,zpsoil,tsoil(1,1,1,is),item1,       &
                     nchanl,nbytes,wdlen,itype,dflag,                   &
                     ipds,igdss,ibdshd,pds,gds,                         &
                     nbufsz,bds,ibds,                                   &
                     mgrib,temxy1,temxy2,itemxy1,itemxy2)

        ipds(8)  = 144  ! Soil moisture in (m**3/m**3)
!        ipds(11) = 1+is ! 1 cm below the surface
        ipds(11) = is    ! 'is' cm below the surface
        IF (myproc == 0) WRITE(6,'(a)') 'Writing qsoil'
        CALL edgfill(qsoil(1,1,1,is), 1,nx,1,nx-1, 1,ny,1,ny-1,          &
                     1,nzsoil,1,nzsoil)
!        CALL wrtgrib(nx,ny,nzsoil,0,zpsoil,qsoil(1,1,1,is),item1,       &
!                     nchanl,nbytes,wdlen,itype,dflag,                   &
!                     ipds,igdss,ibdshd,pds,gds,                         &
!                     nbufsz,bds,ibds,                                   &
!                     mgrib,temxy1,temxy2,itemxy1,itemxy2)
        CALL wrtgrib(nx,ny,nzsoil,1,zpsoil,qsoil(1,1,1,is),item1,       &
                     nchanl,nbytes,wdlen,itype,dflag,                   &
                     ipds,igdss,ibdshd,pds,gds,                         &
                     nbufsz,bds,ibds,                                   &
                     mgrib,temxy1,temxy2,itemxy1,itemxy2)

        ipds(8)  = 223  ! amount of liquid water retained on canopy (m)
        ipds(11) = is   ! 0 cm below the surface
        IF (myproc == 0) WRITE(6,'(a)') 'Writing wetcnpy'
        CALL edgfill(wetcanp(1,1,is),1,nx,1,nx-1,1,ny,1,ny-1,           &
                     1,1,1,1)
        CALL wrtgrib(nx,ny,1,0,z,wetcanp(1,1,is),item1,                 &
                     nchanl,nbytes,wdlen,itype,dflag,                   &
                     ipds,igdss,ibdshd,pds,gds,                         &
                     nbufsz,bds,ibds,                                   &
                     mgrib,temxy1,temxy2,itemxy1,itemxy2)
      END DO
    END IF

    IF (snowout == 1) THEN
      itype = 0                  ! floating array
      ibdshd(3) = itype
      ipds(8)  = 66   ! Snow depth (m)
      ipds(9)  = 1    ! at surface
      ipds(11) = 0    ! at surface
      IF (myproc == 0) WRITE(6,'(a)') 'Writing snowdpth'
      CALL edgfill(snowdpth,1,nx,1,nx-1,1,ny,1,ny-1, 1,1,1,1)
      CALL wrtgrib(nx,ny,1,0,z,snowdpth,item1,                          &
                   nchanl,nbytes,wdlen,itype,dflag,                     &
                   ipds,igdss,ibdshd,pds,gds,                           &
                   nbufsz,bds,ibds,                                     &
                   mgrib,temxy1,temxy2,itemxy1,itemxy2)
    END IF

  END IF
!
!-----------------------------------------------------------------------
!
!  If radout = 1, write out radiation arrays
!
!-----------------------------------------------------------------------
!
  IF( radout == 1 ) THEN
    itype = 0       ! integer array
    ibdshd(3) = itype

    dflag = 0       ! 3-D
    ipds(9)  = 103  ! z-coorinates
    ipds(8)  = 216  ! radiation forcing, radfrc
    CALL edgfill(radfrc,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)
    CALL wrtgrib(nx,ny,nz,0,z,radfrc,item1,                             &
                  nchanl,nbytes,wdlen,itype,dflag,                      &
                  ipds,igdss,ibdshd,pds,gds,                            &
                  nbufsz,bds,ibds,                                      &
                  mgrib,temxy1,temxy2,itemxy1,itemxy2)

    dflag = 1       ! 2-D
    ipds(9)  = 1    ! at surface
    ipds(8)  = 204  ! downward short wave radiation flux
    ipds(11) = 0    ! 0 cm below the surface
    CALL edgfill(radsw,1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)
    CALL wrtgrib(nx,ny,1,0,z,radsw,item1,                               &
                  nchanl,nbytes,wdlen,itype,dflag,                      &
                  ipds,igdss,ibdshd,pds,gds,                            &
                  nbufsz,bds,ibds,                                      &
                  mgrib,temxy1,temxy2,itemxy1,itemxy2)

    dflag = 1       ! 2-D
    ipds(9)  = 1    ! at surface
    ipds(8)  = 232  ! downward total (net) radiation flux
    CALL edgfill(rnflx,1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)
    CALL wrtgrib(nx,ny,1,0,z,rnflx,item1,                               &
                  nchanl,nbytes,wdlen,itype,dflag,                      &
                  ipds,igdss,ibdshd,pds,gds,                            &
                  nbufsz,bds,ibds,                                      &
                  mgrib,temxy1,temxy2,itemxy1,itemxy2)

    dflag = 1       ! 2-D
    ipds(9)  = 1    ! at surface
    ipds(8)  = 111  ! Net shortwave radiation
    CALL edgfill(radswnet,1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)
    CALL wrtgrib(nx,ny,1,0,z,radswnet,item1,                            &
                  nchanl,nbytes,wdlen,itype,dflag,                      &
                  ipds,igdss,ibdshd,pds,gds,                            &
                  nbufsz,bds,ibds,                                      &
                  mgrib,temxy1,temxy2,itemxy1,itemxy2)

    dflag = 1       ! 2-D
    ipds(9)  = 1    ! at surface
    ipds(8)  = 205  ! Incoming longwave radiation
    CALL edgfill(radlwin,1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)
    CALL wrtgrib(nx,ny,1,0,z,radlwin,item1,                             &
                  nchanl,nbytes,wdlen,itype,dflag,                      &
                  ipds,igdss,ibdshd,pds,gds,                            &
                  nbufsz,bds,ibds,                                      &
                  mgrib,temxy1,temxy2,itemxy1,itemxy2)

  END IF
!
!-----------------------------------------------------------------------
!
!  If flxout = 1, write out surface fluxes
!
!-----------------------------------------------------------------------
!
  IF( flxout == 1 ) THEN
    itype = 0       ! integer array
    ibdshd(3) = itype

    dflag = 1       ! 2-D
    ipds(9)  = 1    ! at surface
    ipds(8)  = 124  ! u flux
    ipds(11) = 0    ! 0 cm below the surface
    CALL edgfill(usflx,1,nx,1,nx, 1,ny,1,ny-1, 1,1,1,1)
    CALL wrtgrib(nx,ny,1,0,z,usflx,item1,                               &
                  nchanl,nbytes,wdlen,itype,dflag,                      &
                  ipds,igdss,ibdshd,pds,gds,                            &
                  nbufsz,bds,ibds,                                      &
                  mgrib,temxy1,temxy2,itemxy1,itemxy2)

    ipds(8)  = 125  ! v flux
    CALL edgfill(vsflx,1,nx,1,nx-1, 1,ny,1,ny, 1,1,1,1)
    CALL wrtgrib(nx,ny,1,0,z,vsflx,item1,                               &
                  nchanl,nbytes,wdlen,itype,dflag,                      &
                  ipds,igdss,ibdshd,pds,gds,                            &
                  nbufsz,bds,ibds,                                      &
                  mgrib,temxy1,temxy2,itemxy1,itemxy2)

    ipds(8)  = 122  ! sensible heat flux for ptsflx
    CALL edgfill(ptsflx,1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)
    CALL wrtgrib(nx,ny,1,0,z,ptsflx,item1,                              &
                  nchanl,nbytes,wdlen,itype,dflag,                      &
                  ipds,igdss,ibdshd,pds,gds,                            &
                  nbufsz,bds,ibds,                                      &
                  mgrib,temxy1,temxy2,itemxy1,itemxy2)

    ipds(8)  = 121  ! latent heat flux for qvsflx
    CALL edgfill(qvsflx,1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)
    CALL wrtgrib(nx,ny,1,0,z,qvsflx,item1,                              &
                  nchanl,nbytes,wdlen,itype,dflag,                      &
                  ipds,igdss,ibdshd,pds,gds,                            &
                  nbufsz,bds,ibds,                                      &
                  mgrib,temxy1,temxy2,itemxy1,itemxy2)

  END IF

  DEALLOCATE(item1)
  DEALLOCATE(itemxy1,itemxy2)

  RETURN
END SUBROUTINE gribdump
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE GRIBREAD                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE gribread(nx,ny,nz,nzsoil,nstyps, grdbas,inch,time,x,y,z,     &
           zp, zpsoil,uprt, vprt, wprt, ptprt, pprt,                    &
           qvprt, qscalar, tke,kmh,kmv,                                 &
           ubar, vbar, wbar, ptbar, pbar, rhobar, qvbar,                &
           soiltyp,stypfrct,vegtyp,lai,roufns,veg,                      &
           tsoil,qsoil,wetcanp,snowdpth,                                &
           raing,rainc,prcrate,                                         &
           radfrc,radsw,rnflx,radswnet,radlwin,                         &
           usflx,vsflx,ptsflx,qvsflx)

!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Read in GRIB data set created by ARPS using history dump format
!  No. 10.
!  All data read in are located at the original staggered grid points
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yuhe Liu
!  11/01/1995
!
!  MODIFICATION HISTORY:
!
!  12/09/1998 (Donghai Wang)
!  Added the snow cover.
!
!  1/07/2000 (Eric Kemp)
!  Subroutine now reads in GRIB variable 66 (snow depth (m))
!
!  05/23/2002 (J. Brotzge)
!  Added soil variables to allow for multiple soil schemes  
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!    nzsoil   Number of grid points in the soil  
!
!    grdbas   Data read flag.
!             =1, only grid and base state arrays will be read
!             =0, all arrays will be read based on data
!                 parameter setting.
!    inch     Channel number for binary reading.
!             This channel must be opened for unformatted reading
!             by the calling routine.
!
!  OUTPUT:
!
!    time     Time in seconds of data read from "filename"
!
!    x        x coordinate of grid points in physical/comp. space (m)
!    y        y coordinate of grid points in physical/comp. space (m)
!    z        z coordinate of grid points in computational space (m)
!    zp       z coordinate of grid points in physical space (m)
!    zpsoil   z coordinate of grid points in the soil (m)
!
!    uprt     x component of perturbation velocity (m/s)
!    vprt     y component of perturbation velocity (m/s)
!    wprt     Vertical component of perturbation velocity in
!             Cartesian coordinates (m/s).
!    ptprt    Perturbation potential temperature (K)
!    pprt     Perturbation pressure (Pascal)
!
!    qvprt    Perturbation water vapor mixing ratio (kg/kg)
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
!    ubar     Base state x velocity component (m/s)
!    vbar     Base state y velocity component (m/s)
!    wbar     Base state z velocity component (m/s)
!    ptbar    Base state potential temperature (K)
!    pbar     Base state pressure (Pascal)
!    rhobar   Base state air density (kg/m**3)
!    qvbar    Base state water vapor mixing ratio (kg/kg)
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
!    radswnet Net shortwave radiation, SWin - SWout
!    radlwin  Incoming longwave radiation
!
!    usflx    Surface flux of u-momentum (kg/(m*s**2))
!    vsflx    Surface flux of v-momentum (kg/(m*s**2))
!    ptsflx   Surface heat flux (K*kg/(m**2 * s ))
!    qvsflx   Surface moisture flux of (kg/(m**2 * s))
!
!    ireturn  Return status indicator
!             =0, successful read of all data
!             =1, error reading data
!             =2, end-of-file reached during read attempt
!
!  TEMPORARY:
!
!    tem1     Temporary work array.
!    item1    Integer temporary work array.
!    temxy1   2-D temporary work array.
!    itemxy1  2-D integer temporary work array.
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
  INCLUDE 'indtflg.inc'
  INCLUDE 'globcst.inc'
  INCLUDE 'grid.inc'          ! Grid & map parameters.
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  INTEGER :: nx,ny,nz          ! Number of grid points in 3 directions
  INTEGER :: nzsoil            ! Number of grid points in the soil  

  REAL :: x(nx)                ! x-coord. of the physical and compu
                               ! -tational grid. Defined at u-point(m).
  REAL :: y(ny)                ! y-coord. of the physical and compu
                               ! -tational grid. Defined at v-point(m).
  REAL :: z(nz)                ! z-coord. of the computational grid.
                               ! Defined at w-point on the staggered
                               ! grid(m).
  REAL :: zp(nx,ny,nz)         ! Physical height coordinate defined at
                               ! w-point of the staggered grid(m).
  REAL :: zpsoil(nx,ny,nzsoil) ! Physical height coordinate defined at
                               ! w-point of the soil (m).

  INTEGER :: grdbas            ! Data read flag.
  INTEGER :: inch              ! Channel number for binary reading
  REAL :: time                 ! Time in seconds of data read
                               ! from "filename"

  REAL :: uprt  (nx,ny,nz)     ! Perturbation u-velocity (m/s)
  REAL :: vprt  (nx,ny,nz)     ! Perturbation v-velocity (m/s)
  REAL :: wprt  (nx,ny,nz)     ! Perturbation w-velocity (m/s)
  REAL :: ptprt (nx,ny,nz)     ! Perturbation potential temperature (K)
  REAL :: pprt  (nx,ny,nz)     ! Perturbation pressure (Pascal)

  REAL :: qvprt (nx,ny,nz)     ! Perturbation water vapor mixing
                               ! ratio (kg/kg)

  REAL :: qscalar    (nx,ny,nz,nscalar)

  REAL :: kmh   (nx,ny,nz)     ! Horizontal turb. mixing coef. for
  REAL :: tke   (nx,ny,nz)     ! Turbulent Kinetic Energy ((m/s)**2)
                               ! momentum. ( m**2/s )
  REAL :: kmv   (nx,ny,nz)     ! Vertical turb. mixing coef. for
                               ! momentum. ( m**2/s )

  REAL :: ubar  (nx,ny,nz)     ! Base state u-velocity (m/s)
  REAL :: vbar  (nx,ny,nz)     ! Base state v-velocity (m/s)
  REAL :: wbar  (nx,ny,nz)     ! Base state w-velocity (m/s)
  REAL :: ptbar (nx,ny,nz)     ! Base state potential temperature (K)
  REAL :: pbar  (nx,ny,nz)     ! Base state pressure (Pascal)
  REAL :: rhobar(nx,ny,nz)     ! Base state air density (kg/m**3)
  REAL :: qvbar (nx,ny,nz)     ! Base state water vapor mixing ratio

  INTEGER :: nstyps                  ! Number of soil types
  INTEGER :: soiltyp(nx,ny,nstyps)   ! Soil type
  REAL :: stypfrct(nx,ny,nstyps)  ! Soil type fraction
  INTEGER :: vegtyp (nx,ny)          ! Vegetation type
  REAL :: lai    (nx,ny)          ! Leaf Area Index
  REAL :: roufns (nx,ny)          ! Surface roughness
  REAL :: veg    (nx,ny)          ! Vegetation fraction

  REAL :: tsoil (nx,ny,nzsoil,0:nstyps) ! Soil temperature (K)
  REAL :: qsoil (nx,ny,nzsoil,0:nstyps) ! Soil moisture (m**3/m**3)
  REAL :: wetcanp(nx,ny,0:nstyps)    ! Canopy water amount
  REAL :: snowdpth(nx,ny)            ! Snow depth (m)

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
  REAL :: radswnet(nx,ny)      ! Net shortwave radiation
  REAL :: radlwin(nx,ny)       ! Incoming longwave radiation

  REAL :: usflx (nx,ny)        ! Surface flux of u-momentum (kg/(m*s**2))
  REAL :: vsflx (nx,ny)        ! Surface flux of v-momentum (kg/(m*s**2))
  REAL :: ptsflx(nx,ny)        ! Surface heat flux (K*kg/(m*s**2))
  REAL :: qvsflx(nx,ny)        ! Surface moisture flux (kg/(m**2*s))

  REAL, allocatable :: tem1(:,:,:)     ! Temporary work array
  REAL, allocatable :: temxy1(:,:)     ! temporary work array

  INTEGER, allocatable :: item1(:,:,:) ! Integer temporary work array
  INTEGER, allocatable :: itemxy1(:,:) ! Integer temporary work array

  INTEGER :: ireturn
!
!-----------------------------------------------------------------------
!
!  Parameters describing routine that wrote the gridded data
!
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!
!  Parameters for GRIB packing
!
!-----------------------------------------------------------------------
!
  INTEGER :: nbufsz
  PARAMETER ( nbufsz = 800000 )  ! Size of GRIB buffer

  INTEGER :: ipds(25)          ! PDS integer array
  INTEGER :: igdss(25)         ! GDS integer array for s-grid
  INTEGER :: igdsu(25)         ! GDS integer array for u-grid
  INTEGER :: igdsv(25)         ! GDS integer array for v-grid
  INTEGER :: ibdshd(4)         ! Integer array for BDS header

  CHARACTER (LEN=1) :: bds(nbufsz)   ! BDS ( GRIB Section 4)
  INTEGER :: ibds(nbufsz/4)    ! identical S ( GRIB Section 4)

  EQUIVALENCE ( bds,ibds )

  CHARACTER (LEN=1) :: mgrib(nbufsz) ! Buffer to carry GRIB messages

  INTEGER :: itype             ! Data type indicator
  INTEGER :: dflag             ! Dimension flag
!
!-----------------------------------------------------------------------
!
!  Variables to determine the length of machine word, wdlen
!
!-----------------------------------------------------------------------
!
  INTEGER :: wdlen             ! Length of machine word

  CHARACTER (LEN=8) :: ctema,ctemb

  INTEGER :: itema,itemb

  EQUIVALENCE ( itema,ctema )
  EQUIVALENCE ( itemb,ctemb )

  LOGICAL :: firstcall

  SAVE firstcall, wdlen
  DATA firstcall/ .true. /
  DATA ctema / '12345678' /
!
!-----------------------------------------------------------------------
!
!  Parameters describing this routine
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=40) :: fmtverin
  CHARACTER (LEN=40) :: fmtver410,fmtver530
  INTEGER            :: inverin
  INTEGER            :: inver410, inver530
  PARAMETER (fmtver410='004.10 ARPS Data Format 10',inver410=410)
  PARAMETER (fmtver530='005.30 ARPS Data Format 10',inver530=530)
!
!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: lchanl
  PARAMETER (lchanl=6)      ! Channel number for formatted printing.

  INTEGER :: i,j,k,is
  INTEGER :: nq
  INTEGER :: nqipds(10)

  REAL :: tema1, tema2, temb1, temb2
  REAL :: alatnot(2)

  INTEGER :: istat, nstyp1, nstyp2
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  ALLOCATE (tem1(nx,ny,nz),stat=istat)
  IF (istat /= 0) THEN
    WRITE (6,*) "GRIBREAD: ERROR allocating tem1, returning"
    RETURN
  END IF
  ALLOCATE (item1(nx,ny,nz),stat=istat)
  IF (istat /= 0) THEN
    WRITE (6,*) "GRIBREAD: ERROR allocating item1, returning"
    RETURN
  END IF
  ALLOCATE (temxy1(nx,ny),stat=istat)
  IF (istat /= 0) THEN
    WRITE (6,*) "GRIBREAD: ERROR allocating temxy1, returning"
    RETURN
  END IF
  ALLOCATE (itemxy1(nx,ny),stat=istat)
  IF (istat /= 0) THEN
    WRITE (6,*) "GRIBREAD: ERROR allocating itemxy1, returning"
    RETURN
  END IF

  IF ( firstcall ) THEN

    itemb = itema
    IF ( ctema /= ctemb ) THEN
      wdlen = 4
    ELSE
      wdlen = 8
    END IF

    WRITE (6,'(a,i2,a)')                                                &
        'The machine word length is ',wdlen,' bytes'

    firstcall = .false.

  END IF
!
!-----------------------------------------------------------------------
!
!  Check header info
!
!-----------------------------------------------------------------------
!
  nstyp2 = nstyp

  CALL rdhishd( nx,ny,nz,nzsoil,inch,grdbas,fmtverin,time,x,y,z,        &
                nbufsz,mgrib )

  IF ( fmtverin == fmtver410 ) THEN
    inverin = inver410
  ELSE IF (fmtverin == fmtver530) THEN
    inverin = inver530
  ELSE
    WRITE (6,'(/1x,a,/1x,2a,/1x,a)')                                   &
        'Data format incompatible with the data reader.',               &
        'Format of data is ',fmtverin, '. Job stopped.'
    CALL arpsstop('arpsstop called from gribread with fmtver',1)
  END IF

  nstyp1 = nstyp ! value read in from data
  nstyp = nstyp2 ! original value
  IF ( nstyp1 < 1 ) THEN
    nstyp1 = 1
  END IF
!
!-----------------------------------------------------------------------
!
!  Calculate the length of machine word, wdlen, in bytes, which will
!  be used to pack the data. Assume the length has only two possible
!  value: 4 and 8
!
!-----------------------------------------------------------------------
!
  dx = x(2) - x(1)
  dy = y(2) - y(1)

  alatnot(1) = trulat1
  alatnot(2) = trulat2

  CALL setmapr(mapproj, 1.0, alatnot, trulon)
  CALL lltoxy( 1,1, ctrlat,ctrlon, tema1, tema2 )

  temb1 = tema1 - (FLOAT(nx-3)/2.) * dx
  temb2 = tema2 - (FLOAT(ny-3)/2.) * dy

  CALL setorig( 1, temb1, temb2)

  CALL setcornerll( nx,ny, x,y )               ! set corner lat/lon
!
!-----------------------------------------------------------------------
!
!   Initialize the integer array for BDS header
!
!-----------------------------------------------------------------------
!
  ibdshd(1) = 0
  ibdshd(2) = 0
  ibdshd(3) = 0
  ibdshd(4) = 0
!
!-----------------------------------------------------------------------
!
!  Get the integer array, ipds, from the header.
!
!-----------------------------------------------------------------------
!
  CALL mkipds( nx,ny,nz, ipds )
!
!-----------------------------------------------------------------------
!
!  Get the GDS integer arrays for different grid type, igdsu, igdsv,
!  and igdss, from the header
!
!-----------------------------------------------------------------------
!
  CALL mkigds( nx,ny,nz, 0, igdss )
  CALL mkigds( nx,ny,nz, 1, igdsu )
  CALL mkigds( nx,ny,nz, 2, igdsv )
!
!-----------------------------------------------------------------------
!
!  Get zp
!
!----------------------------------------------------------------------
!
  IF( grdin == 1 .OR. grdbas == 1 ) THEN

    itype = 0       ! floating array
    ibdshd(3) = itype

    dflag = 0       ! 3-D array
    ipds(8)  = 8    ! geometric height (m)
    ipds(9)  = 103  ! z-coorinates
    ipds(18) = 0    ! forecast time = 0 for grid and base variables

    CALL rdgrib(nx,ny,nz,1,z,dflag,inch, wdlen,                         &
                ipds,igdss,ibdshd,                                      &
                nbufsz,bds,ibds,mgrib,                                  &
                zp,item1, temxy1,itemxy1)

    WRITE(lchanl,910) ipds(8),' zp.'

    DO k = 1,nzsoil
      CALL rdgrib(nx,ny,1,1,z,dflag,inch, wdlen,                        &
                ipds,igdss,ibdshd,                                      &
                nbufsz,bds,ibds,mgrib,                                  &
                zpsoil(1,1,k),item1, temxy1,itemxy1)
    END DO
    WRITE(lchanl,910) ipds(8),' zpsoil.'

  END IF  ! grdin
!
!-----------------------------------------------------------------------
!
!  Read in base state fields
!
!----------------------------------------------------------------------
!
  IF( basin == 1 .OR. grdbas == 1 ) THEN

    itype = 0       ! floating array
    ibdshd(3) = itype

    dflag = 0       ! 3-D array
    ipds(9)  = 103  ! z-coorinates
    ipds(18) = 0    ! forecast time = 0 for grid and base variables

    ipds(8)  = 33   ! u wind (m/s)
    CALL rdgrib(nx,ny,nz,0,z,dflag,inch, wdlen,                         &
                ipds,igdsu,ibdshd,                                      &
                nbufsz,bds,ibds,mgrib,                                  &
                ubar,item1, temxy1,itemxy1)

    WRITE(lchanl,910) ipds(8),' ubar.'

    ipds(8)  = 34   ! v wind (m/s)
    CALL rdgrib(nx,ny,nz,0,z,dflag,inch, wdlen,                         &
                ipds,igdsv,ibdshd,                                      &
                nbufsz,bds,ibds,mgrib,                                  &
                vbar,item1, temxy1,itemxy1)

    WRITE(lchanl,910) ipds(8),' vbar.'

    ipds(8)  = 13   ! potential temperature (K)
    CALL rdgrib(nx,ny,nz,0,z,dflag,inch, wdlen,                         &
                ipds,igdss,ibdshd,                                      &
                nbufsz,bds,ibds,mgrib,                                  &
                ptbar,item1, temxy1,itemxy1)

    WRITE(lchanl,910) ipds(8),' ptbar.'

    ipds(8)  = 1    ! pressure (Pascal)
    CALL rdgrib(nx,ny,nz,0,z,dflag,inch, wdlen,                         &
                ipds,igdss,ibdshd,                                      &
                nbufsz,bds,ibds,mgrib,                                  &
                pbar,item1, temxy1,itemxy1)

    WRITE(lchanl,910) ipds(8),' pbar.'

    IF( mstin == 1) THEN
      ipds(8)  = 51   ! specific humidity (kg/kg)
      CALL rdgrib(nx,ny,nz,0,z,dflag,inch, wdlen,                       &
                  ipds,igdss,ibdshd,                                    &
                  nbufsz,bds,ibds,mgrib,                                &
                  qvbar,item1, temxy1,itemxy1)

      WRITE(lchanl,910) ipds(8),' qvbar.'
    END IF

    IF(landin == 1) THEN

      dflag = 1       ! 2-D array
      ipds(9)  = 111  ! depth below land surface

      IF( nstyp1 <= 1 ) THEN
        itype = 1       ! integer array
        ibdshd(3) = itype

        ipds(8)  = 224  ! soil type
        ipds(11) = 0    ! 0 centimeter for surface

        CALL rdgrib(nx,ny,1,0,z,dflag,inch, wdlen,                      &
                    ipds,igdss,ibdshd,                                  &
                    nbufsz,bds,ibds,mgrib,                              &
                    tem1,soiltyp(1,1,1), temxy1,itemxy1)
        WRITE(lchanl,910) ipds(8),' soiltyp.'

      ELSE

        itype = 1       ! integer array
        ibdshd(3) = itype
        ipds(8)  = 224  ! soil type

        DO is=1,nstyp1

          ipds(11) = is ! number to identify soil types

          CALL rdgrib(nx,ny,1,0,z,dflag,inch, wdlen,                    &
                      ipds,igdss,ibdshd,                                &
                      nbufsz,bds,ibds,mgrib,                            &
                      tem1,soiltyp(1,1,is), temxy1,itemxy1)
          WRITE(lchanl,910) ipds(8),' soiltyp.'

        END DO

        itype = 0       ! floating array
        ibdshd(3) = itype
        ipds(8)  = 226  ! redefine 226 for fraction of soil types

        DO is=1,nstyp1

          ipds(11) = is ! fraction of different soil types

          CALL rdgrib(nx,ny,1,0,z,dflag,inch, wdlen,                    &
                      ipds,igdss,ibdshd,                                &
                      nbufsz,bds,ibds,mgrib,                            &
                      stypfrct(1,1,is),item1, temxy1,itemxy1)
          WRITE(lchanl,910) ipds(8),' stypfrct.'

        END DO

      END IF

      CALL fix_stypfrct_nstyp(nx,ny,nstyp1,nstyp,stypfrct)

      itype = 1       ! integer array
      ibdshd(3) = itype

      ipds(8)  = 225  ! vegtyp
      ipds(11) = 0
      CALL rdgrib(nx,ny,1,0,z,dflag,inch, wdlen,                        &
                  ipds,igdss,ibdshd,                                    &
                  nbufsz,bds,ibds,mgrib,                                &
                  tem1,vegtyp, temxy1,itemxy1)

      WRITE(lchanl,910) ipds(8),' vegtyp.'

      itype = 0       ! floating array
      ibdshd(3) = itype

      ipds(8)  = 227  ! redefine 227 for LAI
      CALL rdgrib(nx,ny,1,0,z,dflag,inch, wdlen,                        &
                  ipds,igdss,ibdshd,                                    &
                  nbufsz,bds,ibds,mgrib,                                &
                  lai,item1, temxy1,itemxy1)

      WRITE(lchanl,910) ipds(8),' lai.'

      ipds(8)  = 83   ! surface roughness (m)
      CALL rdgrib(nx,ny,1,0,z,dflag,inch, wdlen,                        &
                  ipds,igdss,ibdshd,                                    &
                  nbufsz,bds,ibds,mgrib,                                &
                  roufns,item1, temxy1,itemxy1)

      WRITE(lchanl,910) ipds(8),' roufns.'

      ipds(8)  = 87   ! vegetation (%)
      CALL rdgrib(nx,ny,1,0,z,dflag,inch, wdlen,                        &
                  ipds,igdss,ibdshd,                                    &
                  nbufsz,bds,ibds,mgrib,                                &
                  tem1,item1, temxy1,itemxy1)

      DO j = 1,ny
        DO i = 1,nx
          veg(i,j) = tem1(i,j,1) / 100.
        END DO
      END DO
      WRITE(lchanl,910) ipds(8),' veg.'

    END IF

  END IF

  IF( grdbas == 1 ) GO TO 930

  ipds(18) = nint(time/60) ! forecast time in minutes

  IF( varin == 1 ) THEN
!
!-----------------------------------------------------------------------
!
!
!----------------------------------------------------------------------
!
    itype = 0                  ! floating array
    ibdshd(3) = itype

    dflag = 0                  ! 3-D array
    ipds(9)  = 103             ! z-coorinates

    ipds(8)  = 33   ! u wind (m/s)
    CALL rdgrib(nx,ny,nz,0,z,dflag,inch, wdlen,                         &
                ipds,igdsu,ibdshd,                                      &
                nbufsz,bds,ibds,mgrib,                                  &
                tem1,item1, temxy1,itemxy1)

    DO k = 1,nz
      DO j = 1,ny
        DO i = 1,nx
          uprt(i,j,k) = tem1(i,j,k) - ubar(i,j,k)
        END DO
      END DO
    END DO
    WRITE(lchanl,910) ipds(8),' uprt.'

    ipds(8)  = 34   ! v wind (m/s)
    CALL rdgrib(nx,ny,nz,0,z,dflag,inch, wdlen,                         &
                ipds,igdsv,ibdshd,                                      &
                nbufsz,bds,ibds,mgrib,                                  &
                tem1,item1, temxy1,itemxy1)

    DO k = 1,nz
      DO j = 1,ny
        DO i = 1,nx
          vprt(i,j,k) = tem1(i,j,k) - vbar(i,j,k)
        END DO
      END DO
    END DO
    WRITE(lchanl,910) ipds(8),' vprt.'

    ipds(8)  = 40   ! w wind (m/s)
    CALL rdgrib(nx,ny,nz,1,z,dflag,inch, wdlen,                         &
                ipds,igdss,ibdshd,                                      &
                nbufsz,bds,ibds,mgrib,                                  &
                wprt,item1, temxy1,itemxy1)

    WRITE(lchanl,910) ipds(8),' wprt.'
!
!-----------------------------------------------------------------------
!
!  Read in scalars
!
!----------------------------------------------------------------------
!
    ipds(8)  = 13   ! potential temperature (K)
    CALL rdgrib(nx,ny,nz,0,z,dflag,inch, wdlen,                         &
                ipds,igdss,ibdshd,                                      &
                nbufsz,bds,ibds,mgrib,                                  &
                tem1,item1, temxy1,itemxy1)

    DO k = 1,nz
      DO j = 1,ny
        DO i = 1,nx
          ptprt(i,j,k) = tem1(i,j,k) - ptbar(i,j,k)
        END DO
      END DO
    END DO
    WRITE(lchanl,910) ipds(8),' ptprt.'

    ipds(8)  = 1    ! pressure (Pascal)
    CALL rdgrib(nx,ny,nz,0,z,dflag,inch, wdlen,                         &
                ipds,igdss,ibdshd,                                      &
                nbufsz,bds,ibds,mgrib,                                  &
                tem1,item1, temxy1,itemxy1)

    DO k = 1,nz
      DO j = 1,ny
        DO i = 1,nx
          pprt(i,j,k) = tem1(i,j,k) - pbar(i,j,k)
        END DO
      END DO
    END DO
    WRITE(lchanl,910) ipds(8),' pprt.'

  END IF
!
!-----------------------------------------------------------------------
!
!  Read in moisture variables
!
!----------------------------------------------------------------------
!
  IF( mstin == 1 ) THEN

    itype = 0                  ! floating array
    ibdshd(3) = itype

    dflag = 0                  ! 3-D array
    ipds(9)  = 103             ! z-coorinates

    ipds(8)  = 51   ! specific humidity (kg/kg)
    CALL rdgrib(nx,ny,nz,0,z,dflag,inch, wdlen,                         &
                ipds,igdss,ibdshd,                                      &
                nbufsz,bds,ibds,mgrib,                                  &
                tem1,item1, temxy1,itemxy1)

    DO k = 1,nz
      DO j = 1,ny
        DO i = 1,nx
          qvprt(i,j,k) = tem1(i,j,k) - qvbar(i,j,k)
        END DO
      END DO
    END DO
    WRITE(lchanl,910) ipds(8),' qvprt.'

    IF (inverin >= inver530) THEN
      nqipds(:) = 255
      IF (p_qcin > 0) nqipds(p_qc) = 153   ! cloud water mixing ratio (kg/kg)
      IF (p_qrin > 0) nqipds(p_qr) = 170   ! rain water mixing ratio (kg/kg)
      IF (p_qiin > 0) nqipds(p_qi) = 178   ! ice mixing ratio (kg/kg)
      IF (p_qsin > 0) nqipds(p_qs) = 171   ! snow mixing ratio (kg/kg)
      IF (p_qhin > 0) nqipds(p_qh) = 179   ! graupel hail mixing ratio (kg/kg)

      DO nq = 1,nscalarin
        ipds(8)  = nqipds(nq)
        CALL rdgrib(nx,ny,nz,0,z,dflag,inch, wdlen,                     &
                    ipds,igdss,ibdshd,                                  &
                    nbufsz,bds,ibds,mgrib,                              &
                    tem1,item1, temxy1,itemxy1)

        IF (p_qcin == nq .AND. P_QC > 0 ) THEN
          qscalar(:,:,:,P_QC) = tem1(:,:,:)
        ELSE IF (p_qrin == nq .AND. P_QR > 0) THEN
          qscalar(:,:,:,P_QR) = tem1(:,:,:)
        ELSE IF (p_qiin == nq .AND. P_QI > 0) THEN
          qscalar(:,:,:,P_QI) = tem1(:,:,:)
        ELSE IF (p_qsin == nq .AND. P_QS > 0) THEN
          qscalar(:,:,:,P_QS) = tem1(:,:,:)
        ELSE IF (p_qhin == nq .AND. P_QH > 0) THEN
          qscalar(:,:,:,P_QH) = tem1(:,:,:)
        END IF
            
      END DO

    ELSE
      ipds(8)  = 153  ! cloud water mixing ratio (kg/kg)
      CALL rdgrib(nx,ny,nz,0,z,dflag,inch, wdlen,                       &
                  ipds,igdss,ibdshd,                                    &
                  nbufsz,bds,ibds,mgrib,                                &
                  tem1,item1, temxy1,itemxy1)

      IF (P_QC > 0) qscalar(:,:,:,P_QC) = tem1(:,:,:)
  
      WRITE(lchanl,910) ipds(8),' qc.'
  
      ipds(8)  = 170  ! rain water mixing ratio (kg/kg)
      CALL rdgrib(nx,ny,nz,0,z,dflag,inch, wdlen,                       &
                  ipds,igdss,ibdshd,                                    &
                  nbufsz,bds,ibds,mgrib,                                &
                  tem1,item1, temxy1,itemxy1)

      IF (P_QR > 0) qscalar(:,:,:,P_QR) = tem1(:,:,:)
  
      WRITE(lchanl,910) ipds(8),' qr.'
  
      IF( icein == 1 ) THEN
  
        dflag = 0       ! 3-D
        ipds(9)  = 103  ! z-coorinates
  
        ipds(8)  = 178  ! ice mixing ratio (kg/kg)
        CALL rdgrib(nx,ny,nz,0,z,dflag,inch, wdlen,                     &
                    ipds,igdss,ibdshd,                                  &
                    nbufsz,bds,ibds,mgrib,                              &
                    tem1,item1, temxy1,itemxy1)
  
        IF (P_QI > 0) qscalar(:,:,:,P_QI) = tem1(:,:,:)

        WRITE(lchanl,910) ipds(8),' qi.'
  
        ipds(8)  = 171  ! snow mixing ratio (kg/kg)
        CALL rdgrib(nx,ny,nz,0,z,dflag,inch, wdlen,                     &
                    ipds,igdss,ibdshd,                                  &
                    nbufsz,bds,ibds,mgrib,                              &
                    tem1,item1, temxy1,itemxy1)
  
        IF (P_QS > 0) qscalar(:,:,:,P_QS) = tem1(:,:,:)

        WRITE(lchanl,910) ipds(8),' qs.'
  
        ipds(8)  = 179  ! graupel hail mixing ratio (kg/kg)
        CALL rdgrib(nx,ny,nz,0,z,dflag,inch, wdlen,                     &
                    ipds,igdss,ibdshd,                                  &
                    nbufsz,bds,ibds,mgrib,                              &
                    tem1,item1, temxy1,itemxy1)
  
        IF (P_QH > 0) qscalar(:,:,:,P_QH) = tem1(:,:,:)

        WRITE(lchanl,910) ipds(8),' qh.'
  
      END IF
    END IF

    IF( rainin == 1 ) THEN
      dflag = 1       ! 2-D
      ipds(9)  = 1    ! at surface

      ipds(8)  = 62   ! large scale precipitation for raing (kg/m**2)
      ipds(11) = 0    ! at surface
      CALL rdgrib(nx,ny,1,0,z,dflag,inch, wdlen,                        &
                  ipds,igdss,ibdshd,                                    &
                  nbufsz,bds,ibds,mgrib,                                &
                  raing,item1, temxy1,itemxy1)

      WRITE(lchanl,910) ipds(8),' raing.'

      ipds(8)  = 63   ! convective precipitation for rainc (kg/m**2)
      ipds(11) = 0    ! at surface
      CALL rdgrib(nx,ny,1,0,z,dflag,inch, wdlen,                        &
                  ipds,igdss,ibdshd,                                    &
                  nbufsz,bds,ibds,mgrib,                                &
                  rainc,item1, temxy1,itemxy1)

      WRITE(lchanl,910) ipds(8),' rainc.'
    END IF

    IF( prcin == 1 ) THEN
      dflag = 1       ! 2-D
      ipds(8)  = 59   ! precipitation rates (kg/m**2/s)
      ipds(9)  = 111  ! use 111 for saving more than one precip.

      ipds(11) = 1    ! total precipitation rate
      CALL rdgrib(nx,ny,1,0,z,dflag,inch, wdlen,                        &
                  ipds,igdss,ibdshd,                                    &
                  nbufsz,bds,ibds,mgrib,                                &
                  prcrate(1,1,1),item1, temxy1,itemxy1)
      WRITE(lchanl,910) ipds(8),' prcrate1.'

      ipds(11) = 2    ! grid scale precipitation rate
      CALL rdgrib(nx,ny,1,0,z,dflag,inch, wdlen,                        &
                  ipds,igdss,ibdshd,                                    &
                  nbufsz,bds,ibds,mgrib,                                &
                  prcrate(1,1,2),item1, temxy1,itemxy1)
      WRITE(lchanl,910) ipds(8),' prcrate2.'

      ipds(11) = 3    ! cumulus precipitation rate
      CALL rdgrib(nx,ny,1,0,z,dflag,inch, wdlen,                        &
                  ipds,igdss,ibdshd,                                    &
                  nbufsz,bds,ibds,mgrib,                                &
                  prcrate(1,1,3),item1, temxy1,itemxy1)
      WRITE(lchanl,910) ipds(8),' prcrate3.'

      ipds(11) = 4    ! microphysics precipitation rate
      CALL rdgrib(nx,ny,1,0,z,dflag,inch, wdlen,                        &
                  ipds,igdss,ibdshd,                                    &
                  nbufsz,bds,ibds,mgrib,                                &
                  prcrate(1,1,4),item1, temxy1,itemxy1)
      WRITE(lchanl,910) ipds(8),' prcrate4.'

    END IF

  END IF

  IF( tkein == 1 ) THEN

    dflag = 0       ! 3-D
    ipds(9)  = 103  ! z-coorinates

    ipds(8)  = 158  ! tke (J kg**-1)
    CALL rdgrib(nx,ny,nz,0,z,dflag,inch, wdlen,                         &
                ipds,igdss,ibdshd,                                      &
                nbufsz,bds,ibds,mgrib,                                  &
                tke,item1, temxy1,itemxy1)

    WRITE(lchanl,910) ipds(8),' tke.'

  END IF

  IF( trbin == 1 ) THEN

    dflag = 0       ! 3-D
    ipds(9)  = 103  ! z-coorinates

    ipds(8)  = 185  ! redefine 185 for kmh
    CALL rdgrib(nx,ny,nz,0,z,dflag,inch, wdlen,                         &
                ipds,igdss,ibdshd,                                      &
                nbufsz,bds,ibds,mgrib,                                  &
                kmh,item1, temxy1,itemxy1)

    WRITE(lchanl,910) ipds(8),' kmh.'

    dflag = 0       ! 3-D
    ipds(9)  = 103  ! z-coorinates

    ipds(8)  = 186  ! redefine 186 for kmv
    CALL rdgrib(nx,ny,nz,0,z,dflag,inch, wdlen,                         &
                ipds,igdss,ibdshd,                                      &
                nbufsz,bds,ibds,mgrib,                                  &
                kmv,item1, temxy1,itemxy1)

    WRITE(lchanl,910) ipds(8),' kmv.'

  END IF

  IF( sfcin == 1 ) THEN

    dflag = 1         ! 2-D
    ipds(9)  = 111    ! below land surface

    IF( nstyp1 <= 1 ) THEN

      ipds(8)  = 85   ! soil temperature (K)
!      ipds(11) = 1    ! 1 cm below the surface
      ipds(11) = 0    ! 0 cm below the surface

      CALL rdgrib(nx,ny,nzsoil,0,zpsoil,dflag,inch, wdlen,              &
                  ipds,igdss,ibdshd,                                    &
                  nbufsz,bds,ibds,mgrib,                                &
                  tsoil(1,1,1,0),item1, temxy1,itemxy1)
      WRITE(lchanl,910) ipds(8),' tsoil.'

      ipds(8)  = 144  ! Soil moisture in (m**3/m**3)
!      ipds(11) = 1    ! 1 cm below the surface
      ipds(11) = 0    ! 0 cm below the surface
      CALL rdgrib(nx,ny,nzsoil,0,zpsoil,dflag,inch, wdlen,              &
                  ipds,igdss,ibdshd,                                    &
                  nbufsz,bds,ibds,mgrib,                                &
                  qsoil(1,1,1,0),item1, temxy1,itemxy1)
      WRITE(lchanl,910) ipds(8),' qsoil.'

      ipds(8)  = 223  ! amount of liquid water retained on canopy (m)
      ipds(11) = 0    ! 0 cm below the surface
      CALL rdgrib(nx,ny,1,0,z,dflag,inch, wdlen,                        &
                  ipds,igdss,ibdshd,                                    &
                  nbufsz,bds,ibds,mgrib,                                &
                  wetcanp(1,1,0),item1, temxy1,itemxy1)
      WRITE(lchanl,910) ipds(8),' wetcanp.'

    ELSE

      DO is=0,nstyp1

        ipds(8)  = 85   ! soil temperature (K)
!        ipds(11) = 1+is ! 1 cm below the surface
        ipds(11) = is ! is cm below the surface
        CALL rdgrib(nx,ny,nzsoil,0,zpsoil,dflag,inch, wdlen,            &
                    ipds,igdss,ibdshd,                                  &
                    nbufsz,bds,ibds,mgrib,                              &
                    tsoil(1,1,1,is),item1, temxy1,itemxy1)
        WRITE(lchanl,910) ipds(8),' tsoil.'

        ipds(8)  = 144  ! Soil moisture in (m**3/m**3)
!        ipds(11) = 1+is ! 1 cm below the surface
        ipds(11) = is ! is cm below the surface
        CALL rdgrib(nx,ny,nzsoil,0,zpsoil,dflag,inch, wdlen,            &
                    ipds,igdss,ibdshd,                                  &
                    nbufsz,bds,ibds,mgrib,                              &
                    qsoil(1,1,1,is),item1, temxy1,itemxy1)
        WRITE(lchanl,910) ipds(8),' qsoil.'

        ipds(8)  = 223  ! amount of liquid water retained on canopy (m)
        ipds(11) = is   ! 0 cm below the surface
        CALL rdgrib(nx,ny,1,0,z,dflag,inch, wdlen,                      &
                    ipds,igdss,ibdshd,                                  &
                    nbufsz,bds,ibds,mgrib,                              &
                    wetcanp(1,1,is),item1, temxy1,itemxy1)
        WRITE(lchanl,910) ipds(8),' wetcanp.'

      END DO

    END IF

    CALL fix_soil_nstyp(nx,ny,nzsoil,nstyp1,nstyp,tsoil,               &
          qsoil,wetcanp)

    IF (snowcin == 1) THEN  ! Snow cover is depricated
      itype = 1
      ibdshd(3) = itype
      ipds(8)  = 143  ! Snow cover
      ipds(9)  = 111
      ipds(11) = 0
      CALL rdgrib(nx,ny,1,0,z,dflag,inch, wdlen,                        &
                  ipds,igdss,ibdshd,                                    &
                  nbufsz,bds,ibds,mgrib,                                &
                  tem1,item1, temxy1,itemxy1)
      WRITE(lchanl,910) ipds(8),' snowcvr -- discarding.'
    END IF

    IF (snowin == 1) THEN
      itype = 0       ! floating array
      ibdshd(3) = itype
      ipds(8)  = 66   ! Snow depth (m)
      ipds(9)  = 1    ! at surface
      ipds(11) = 0    ! at surface
      CALL rdgrib(nx,ny,1,0,z,dflag,inch, wdlen,                        &
                  ipds,igdss,ibdshd,                                    &
                  nbufsz,bds,ibds,mgrib,                                &
                  snowdpth,item1, temxy1,itemxy1)
      WRITE(lchanl,910) ipds(8),' snowdpth.'
    END IF

  END IF

  IF( radin == 1 ) THEN
    itype = 0       ! floating array
    ibdshd(3) = itype

    dflag = 0       ! 3-D
    ipds(9)  = 103  ! z-coorinates

    ipds(8)  = 216  ! radiation forcing, radfrc
    CALL rdgrib(nx,ny,nz,0,z,dflag,inch, wdlen,                         &
                ipds,igdss,ibdshd,                                      &
                nbufsz,bds,ibds,mgrib,                                  &
                radfrc,item1, temxy1,itemxy1)
    WRITE(lchanl,910) ipds(8),' radfrc.'

    dflag = 1       ! 2-D
    ipds(9)  = 1    ! at surface
    ipds(11) = 0    ! at surface

    ipds(8)  = 204  ! downward short wave radiation flux
    CALL rdgrib(nx,ny,1,0,z,dflag,inch, wdlen,                          &
                ipds,igdss,ibdshd,                                      &
                nbufsz,bds,ibds,mgrib,                                  &
                radsw,item1, temxy1,itemxy1)
    WRITE(lchanl,910) ipds(8),' radsw.'

    ipds(8)  = 232  ! downward total (net) radiation flux
    CALL rdgrib(nx,ny,1,0,z,dflag,inch, wdlen,                          &
                ipds,igdss,ibdshd,                                      &
                nbufsz,bds,ibds,mgrib,                                  &
                rnflx,item1, temxy1,itemxy1)
    WRITE(lchanl,910) ipds(8),' rnflx.'

    ipds(8)  = 111  ! downward total (net) radiation flux
    CALL rdgrib(nx,ny,1,0,z,dflag,inch, wdlen,                          &
                ipds,igdss,ibdshd,                                      &
                nbufsz,bds,ibds,mgrib,                                  &
                radswnet,item1, temxy1,itemxy1)
    WRITE(lchanl,910) ipds(8),' radswnet.'

    ipds(8)  = 205  ! downward total (net) radiation flux
    CALL rdgrib(nx,ny,1,0,z,dflag,inch, wdlen,                          &
                ipds,igdss,ibdshd,                                      &
                nbufsz,bds,ibds,mgrib,                                  &
                radlwin,item1, temxy1,itemxy1)
    WRITE(lchanl,910) ipds(8),' radlwin.'

  END IF

  IF( flxin == 1 ) THEN
    itype = 0       ! integer array
    ibdshd(3) = itype

    ipds(11) = 0    ! at surface
    ipds(8)  = 124  ! u flux
    CALL rdgrib(nx,ny,1,0,z,dflag,inch, wdlen,                          &
                ipds,igdss,ibdshd,                                      &
                nbufsz,bds,ibds,mgrib,                                  &
                usflx,item1, temxy1,itemxy1)
    WRITE(lchanl,910) ipds(8),' usflx.'

    ipds(8)  = 125  ! v flux
    CALL rdgrib(nx,ny,1,0,z,dflag,inch, wdlen,                          &
                ipds,igdss,ibdshd,                                      &
                nbufsz,bds,ibds,mgrib,                                  &
                vsflx,item1, temxy1,itemxy1)
    WRITE(lchanl,910) ipds(8),' vsflx.'

    ipds(8)  = 122  ! sensible heat flux for ptsflx
    CALL rdgrib(nx,ny,1,0,z,dflag,inch, wdlen,                          &
                ipds,igdss,ibdshd,                                      &
                nbufsz,bds,ibds,mgrib,                                  &
                ptsflx,item1, temxy1,itemxy1)
    WRITE(lchanl,910) ipds(8),' ptsflx.'

    ipds(8)  = 121  ! latent heat flux for qvsflx
    CALL rdgrib(nx,ny,1,0,z,dflag,inch, wdlen,                          &
                ipds,igdss,ibdshd,                                      &
                nbufsz,bds,ibds,mgrib,                                  &
                qvsflx,item1, temxy1,itemxy1)
    WRITE(lchanl,910) ipds(8),' qvsflx.'

  END IF

  910   FORMAT(1X,'Field ',i4,' was read into array',a)

!
!-----------------------------------------------------------------------
!
!  Friendly exit message
!
!----------------------------------------------------------------------
!
  930   CONTINUE

  WRITE(6,'(/a,F8.1,a/)')                                               &
      ' Data at time=', time/60,' (min) were successfully read.'

  ireturn = 0

  DEALLOCATE (tem1,stat=istat)
  DEALLOCATE (item1,stat=istat)
  DEALLOCATE (temxy1,stat=istat)
  DEALLOCATE (itemxy1,stat=istat)

  RETURN

END SUBROUTINE gribread
!
!##################################################################
!##################################################################
!######                                                      ######
!######               SUBROUTINE WRTGRIB                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE wrtgrib(nx,ny,nz,wstag,z,fvar,ivar,                   &
           nchanl,nbytes,wdlen,itype,dflag,                      &
           ipds,igds,ibdshd,pds,gds,                             &
           nbufsz,bds,ibds,                                      &
           mgrib,tem1,tem2,item1,item2)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Write variable, fvar for floating array or ivar for integer array,
!  into file pointer nchanl which points to a GRIB file.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yuhe Liu
!  11/01/1995
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!    wstag    Staggering flag of w-grid
!    z        z coordinate of grid points in computational space (m)
!
!    fvar     Floating array to be written into GRIB file
!    ivar     Integer  array to be written into GRIB file
!
!    nchanl   FILE pointer of GRIB stream file which was opened by
!             a C program, GOPEN
!
!    itype    Data type: 0 for floating array and 1 for integer array
!    dflag    Dimension flag: 0 for 3-D array and 1 for 2-D array
!    nbytes   Byte position at which the GRIB message starts
!    wdlen    Length of machine word (i.e. size of integer) in bytes
!
!    ipds     Integer array to carry the parameters for generating
!             PDS (Section 1).
!    igds     Integer array to carry the parameters for generating
!             GDS (Section 3).
!
!    pds      PDS (GRIB Section 1)
!    gds      GDS (GRIB Section 3)
!
!    nbufsz   Size of buffer of a GRIB message array
!    bds      BDS (GRIB Section 4)
!    mgrib    Working buffer of GRIB message array
!
!  OUTPUT:
!
!    nbtyes   Position at which the next GRIB message starts
!
!  WORK ARRAY:
!
!    tem1     Temporary work array.
!    item1    Integer temporary work array.
!    tem2     Temporary work array.
!    item2    Integer temporary work array.
!
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx,ny,nz          ! Number of grid points in x,y,z dir.

  INTEGER :: wstag             ! Staggering flag for w-grid

  REAL :: z(nz)                ! The z-coord. of the computational grid.
                               ! Defined at w-point on the staggered grid.
  REAL :: fvar(nx,ny,nz)    ! Floating array
  INTEGER :: ivar(nx,ny,nz) ! Integer  array

  INTEGER :: nchanl            ! FILE pointer indicates the GRIB file
  INTEGER :: itype             ! Data type: 0 - floating, 1 - integer
  INTEGER :: dflag             ! Dimension flag: 0 - 3-D, 1 - 2-D
  INTEGER :: nbytes            ! Starting byte position
  INTEGER :: wdlen             ! Length of machine word

  INTEGER :: ipds(25)          ! ipds
  INTEGER :: igds(25)          ! igds
  INTEGER :: ibdshd(4)         ! BDS header

  CHARACTER (LEN=1) :: pds(28)       ! PDS
  CHARACTER (LEN=1) :: gds(42)       ! GDS

  INTEGER :: nbufsz            ! Size of GRIB message array
  INTEGER :: ibds(nbufsz/wdlen) ! Identical to BDS
  CHARACTER (LEN=1) :: bds(nbufsz)   ! BDS
  CHARACTER (LEN=1) :: mgrib(nbufsz) ! Working buffer

  REAL :: tem1 (nx,ny)         ! Temporary work array
  REAL :: tem2 (nx,ny)         ! Temporary work array

  INTEGER :: item1(nx,ny)      ! Integer temporary work array
  INTEGER :: item2(nx,ny)      ! Integer temporary work array
!
!-----------------------------------------------------------------------
!
!  GRIB parameters
!
!-----------------------------------------------------------------------
!
  INTEGER :: msglen    ! Length of GRIB message
  INTEGER :: npts
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k, wstag1
  REAL :: tema
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

  npts = nx*ny

  IF ( dflag == 1 ) THEN

    DO j = 1,ny
      DO i = 1,nx
        tem1 (i,j) = fvar(i,j,1)
        item1(i,j) = ivar(i,j,1)
      END DO
    END DO

    CALL gribenc(itype,wdlen,grbpkbit,npts,tem1,item1,                  &
                 ipds,igds,ibdshd,pds,gds,                              &
                 nbufsz,bds,ibds,msglen,mgrib,                          &
                 tem2,item2)

    WRITE (nchanl) (mgrib(i),i=1,msglen)

    RETURN
  END IF

  IF ( wstag > 1 .OR. wstag < 0 ) THEN
    WRITE(6,'(a,i1)')                                                   &
        'Do not know the w-grid type, wstag = ',wstag,                  &
        'Reset it to scalar grid.'
    wstag1 = 0
  ELSE
    wstag1 = wstag
  END IF

  DO k = 1,nz-1
    DO j = 1,ny
      DO i = 1,nx
        tem1 (i,j) = fvar(i,j,k)
        item1(i,j) = ivar(i,j,k)
      END DO
    END DO

    IF ( wstag1 == 0 ) THEN
      tema = 0.5 * (z(k) + z(k+1))
    ELSE
      tema = z(k)
    END IF
    ipds(11) = nint(tema)

    CALL gribenc(itype,wdlen,grbpkbit,npts,tem1,item1,                  &
                 ipds,igds,ibdshd,pds,gds,                              &
                 nbufsz,bds,ibds,msglen,mgrib,                          &
                 tem2,item2)

    WRITE (nchanl) (mgrib(i),i=1,msglen)

  END DO

  DO j = 1,ny
    DO i = 1,nx
      tem1 (i,j) = fvar(i,j,nz)
      item1(i,j) = ivar(i,j,nz)
    END DO
  END DO

  ipds(11) = nint(z(nz))

  CALL gribenc(itype,wdlen,grbpkbit,npts,tem1,item1,                    &
               ipds,igds,ibdshd,pds,gds,                                &
               nbufsz,bds,ibds,msglen,mgrib,                            &
               tem2,item2)

  WRITE (nchanl) (mgrib(i),i=1,msglen)

  RETURN
END SUBROUTINE wrtgrib
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE RDGRIB                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE rdgrib(nx,ny,nz,wstag,z,dflag,nchanl, wdlen,          &
           ipds,igds,ibdshd,                                            &
           nbufsz,bds,ibds,mgrib,                                       &
           fvar,ivar, tem1,item1)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Write variable, fvar for floating array or ivar for integer array,
!  into file pointer nchanl which points to a GRIB file.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yuhe Liu
!  11/01/1995
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       number of grid points in the vertical
!
!    dflag    Dimension flag, 0: 3-D, 1: 2-D
!    wstag    Staggering flag of w-grid
!
!    nchanl   FILE pointer of GRIB stream file which was opened by
!             a C program, GOPEN
!
!  OUTPUT:
!
!    fvar     Floating array to be written into GRIB file
!    ivar     Integer  array to be written into GRIB file
!
!    wdlen    Length of machine word (i.e. size of integer) in bytes
!
!    ipds     Integer array to carry the parameters for generating
!             PDS (Section 1).
!    igds     Integer array to carry the parameters for generating
!             GDS (Section 3).
!
!    nbufsz   Size of buffer of a GRIB message array
!    bds      BDS (GRIB Section 4)
!    mgrib    Working buffer of GRIB message array
!
!    nbtyes   Position at which the next GRIB message starts
!
!  WORK ARRAY:
!
!    tem1     Temporary work array.
!    item1    Integer temporary work array.
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx,ny,nz          ! Number of grid points in x,y,z dir.

  INTEGER :: wstag             ! Staggering flag for w-grid

  REAL :: z(nz)                ! Vertical coordinates

  INTEGER :: dflag             ! Dimension flag, 0: 3-D, 1: 2-D

  INTEGER :: nchanl            ! FILE pointer indicates the GRIB file

  INTEGER :: wdlen             ! Length of machine word

  INTEGER :: ipds(25)          ! ipds
  INTEGER :: igds(25)          ! igds
  INTEGER :: ibdshd(4)         ! BDS header

  INTEGER :: nbufsz            ! Size of GRIB message array
  CHARACTER (LEN=1) :: mgrib(nbufsz) ! GRIB message

  REAL :: fvar(nx,ny,nz)    ! Floating array
  INTEGER :: ivar(nx,ny,nz)    ! Integer  array

  CHARACTER (LEN=1) :: bds(nbufsz)   ! BDS
  INTEGER :: ibds(nbufsz/wdlen) ! Identical to BDS

  REAL :: tem1 (nx,ny)         ! Temporary work array
  INTEGER :: item1(nx,ny)      ! Integer temporary work array
!
!-----------------------------------------------------------------------
!
!  GRIB parameters
!
!-----------------------------------------------------------------------
!
  INTEGER :: msglen    ! Length of GRIB message
  INTEGER :: npts
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k, wstag1
  REAL :: tema
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
!
!-----------------------------------------------------------------------
!
!  External function
!
!-----------------------------------------------------------------------
!
  INTEGER :: char2i
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  npts = nx*ny

  IF ( dflag == 1 ) THEN
    READ (nchanl) (mgrib(i),i=1,8)

    msglen = char2i(mgrib(5))*65536                                      &
           + char2i(mgrib(6))*256                                        &
           + char2i(mgrib(7))

    BACKSPACE (nchanl)

    READ (nchanl) (mgrib(i),i=1,msglen)

    CALL gribdec(npts,nz,                                               &
                 wdlen,ipds,igds,ibdshd,nbufsz,mgrib,                   &
                 tem1,item1, bds,ibds)

    DO j = 1,ny
      DO i = 1,nx
        fvar(i,j,1) = tem1 (i,j)
        ivar(i,j,1) = item1(i,j)
      END DO
    END DO

    RETURN
  END IF

  IF ( wstag > 1 .OR. wstag < 0 ) THEN
    WRITE(6,'(a,i1)')                                                   &
        'Do not know the w-grid type, wstag = ',wstag,                  &
        'Reset it to scalar grid.'
    wstag1 = 0
  ELSE
    wstag1 = wstag
  END IF

  DO k = 1,nz-1

    IF ( wstag1 == 0 ) THEN
      tema = 0.5 * (z(k) + z(k+1))
    ELSE
      tema = z(k)
    END IF

    ipds(11) = nint(tema)

    READ (nchanl) (mgrib(i),i=1,8)

    msglen = char2i(mgrib(5))*65536                                      &
           + char2i(mgrib(6))*256                                        &
           + char2i(mgrib(7))

    BACKSPACE (nchanl)

    READ (nchanl) (mgrib(i),i=1,msglen)

    CALL gribdec(npts,nz,                                               &
                 wdlen,ipds,igds,ibdshd,nbufsz,mgrib,                   &
                 tem1,item1, bds,ibds)

    DO j = 1,ny
      DO i = 1,nx
        fvar(i,j,k) = tem1 (i,j)
        ivar(i,j,k) = item1(i,j)
      END DO
    END DO

  END DO

  ipds(11) = nint(z(nz))

  READ (nchanl) (mgrib(i),i=1,8)

  msglen = char2i(mgrib(5))*65536                                        &
         + char2i(mgrib(6))*256                                          &
         + char2i(mgrib(7))

  BACKSPACE (nchanl)

  READ (nchanl) (mgrib(i),i=1,msglen)

  CALL gribdec(npts,nz,                                                 &
               wdlen,ipds,igds,ibdshd,nbufsz,mgrib,                     &
               tem1,item1, bds,ibds)

  DO j = 1,ny
    DO i = 1,nx
      fvar(i,j,nz) = tem1 (i,j)
      ivar(i,j,nz) = item1(i,j)
    END DO
  END DO

  RETURN
END SUBROUTINE rdgrib
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE MKIPDS                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE mkipds(nx,ny,nz, ipds)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Make integer product data array ipds
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yuhe Liu
!  11/27/1995
!
!  MODIFICATION HISTORY:
!
!  2003-12-29  Richard Carpenter
!  Assigned values for Center and Process/Model.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!  OUTPUT:
!
!    ipds     IPDS array
!
!  WORK ARRAY:
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx,ny,nz          ! Number of grid points in 3 directions

  INTEGER :: ipds(25)
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!wdt - ipds(3,4)
  ipds(1) = 28          ! number of bytes in PDS section
  ipds(2) = 0           ! parameter table version number, unknown
  ipds(3) = 192         ! ID number of originating center. CAPS?
  ipds(4) = 192         ! ID number of model. ARPS/CAPS?
  ipds(5) = 255         ! grid ID, 255 for self-defined GDS
  ipds(6) = 1           ! GDS flag, 1 = GDS section included
  ipds(7) = 0           ! BMS flag, 0 = no BMS section included
  ipds(8) = 0           ! indicator of parameter and unit
                        ! (table 2), depends on variables
  ipds(9) = 103         ! indicator of type of level (table 3)
                        ! depends on which variable
                        ! 103 for z coordinates in meters
                        ! 111 for soil layers in centimeters
! both use two octets: ipds(10) & ipds(11)
  ipds(10) = 0          ! value 1 of level, N/A
  ipds(11) = 0          ! value 2 of level (value of level)

  IF ( year <= 2000 ) THEN
    ipds(12) = year-1900 ! year of century
    ipds(23) = 20        ! century (20, change to 21 on Jan 1, 2001)
  ELSE
    ipds(12) = year-2000 ! year of century
    ipds(23) = 21        ! century (20, change to 21 on Jan 1, 2001)
  END IF

  ipds(13) = month      ! month of year
  ipds(14) = day        ! day of month
  ipds(15) = hour       ! hour of day
  ipds(16) = minute     ! minute of hour
  ipds(17) = 0          ! forecast time unit, minutes
  ipds(18) = 0          ! forecast time P1, or current time, curtim
  ipds(19) = 0          ! forecast time P2, N/A
  ipds(20) = 10         ! time range indicator,
                        ! 10 = use two octets from ipds(18)
  ipds(21) = 0          ! number include in average, N/A
  ipds(22) = 0          ! number missing from average, N/A
  ipds(24) = 0          ! subcenter ID number
  ipds(25) = 0          ! scaling power of 10,
                        ! depends on which variable

  RETURN
END SUBROUTINE mkipds
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE MKIGDS                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE mkigds(nx,ny,nz,stagger,igds)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Make integer grid description array IGDS.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yuhe Liu
!  11/01/1995
!
!  MODIFICATION HISTORY:
!
!  2003-12-02  Richard Carpenter
!  Fixed projection center flag [GDS(27)=igds(12)] for S Hemi Lamb Conf.
!
!  2003-12-29  Richard Carpenter
!  Fixed lat/lon increments for Mercator and Lat/Lon grids.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!    stagger  Flag for staggered grid
!           = 0, s-grid, for scalar grid
!           = 1, u-grid, for u-vector grid
!           = 0, v-grid, for v-vector grid
!
!  OUTPUT:
!
!    igds     IPDS array
!
!  WORK ARRAY:
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx,ny,nz          ! Number of grid points in 3 directions

  INTEGER :: stagger

  INTEGER :: igds(25)
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i

  REAL :: swlat,swlon, nelat,nelon
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'grid.inc'          ! Grid & map parameters.
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  IF ( stagger == 0 ) THEN
    swlat = swlats
    swlon = swlons
    nelat = nelats
    nelon = nelons
  ELSE IF ( stagger == 1 ) THEN
    swlat = swlatu
    swlon = swlonu
    nelat = nelatu
    nelon = nelonu
  ELSE IF ( stagger == 2 ) THEN
    swlat = swlatv
    swlon = swlonv
    nelat = nelatv
    nelon = nelonv
  ELSE
    WRITE(6,'(a)')                                                      &
        'Do not know the grid type. Set the grid type to scalar grid'
  END IF

  IF ( ABS(mapproj) == 1 ) THEN       ! Polar stereographic
    igds(1) = nz        ! number of vertical (z) coordinates
    igds(2) = 255       ! location (in octets) of vertical
                        ! coordinate list
    igds(3) = 5         ! data representation type (table 6)
    igds(4) = nx        ! number of x-coordinates
    igds(5) = ny        ! number of y-coordinates
    igds(6) = nint(swlat*1000)  ! latitude  at southwest corner
    igds(7) = nint(swlon*1000)  ! longitude at southwest corner
    igds(8) = 8         ! u & v relative to x & y direction
    igds(9) = nint(trulon*1000) ! true longitude of projection
    igds(10) = dx
    igds(11) = dy

    IF ( trulat1 >= 0. ) THEN
      igds(12) = 0      ! northern pole on plane
    ELSE
      igds(12) = 1      ! southern pole on plane
    END IF

    igds(13) = 64       ! scanning flag, 64 for +i & +j direction
    DO i = 14,25
      igds(i) = 0       ! unused
    END DO

  ELSE IF ( ABS(mapproj) == 2 .OR. mapproj == 0 ) THEN ! Lambert conformal
                                                  ! and no projection too

    igds( 1) = nz       ! number of vertical (z) coordinates
    igds( 2) = 255      ! location (in octets) of vertical
                        ! coordinate list
    igds( 3) = 3        ! data representation type (table 6)
    igds( 4) = nx       ! number of x-coordinates
    igds( 5) = ny       ! number of y-coordinates
    igds( 6) = nint(swlat*1000)  ! latitude  at southwest corner
    igds( 7) = nint(swlon*1000)  ! longitude at southwest corner
    igds( 8) = 8        ! u & v relative to x & y direction
    igds( 9) = nint(trulon*1000) ! true longitude of projection
    igds(10) = dx
    igds(11) = dy

    IF ( trulat1 >= 0. ) THEN
      igds(12) = 0      ! northern pole on plane
    ELSE
      igds(12) = 128    ! southern pole on plane
    END IF

    igds(13) = 64       ! scanning flag, 64 for +i & +j direction
    igds(14) = 0        ! unused
    igds(15) = nint(trulat1*1000)  ! first true latitude (millidegree)
    igds(16) = nint(trulat2*1000)  ! second true latitude (millidegree)
    igds(17) = -90000   ! latitude  of south pole (millidegree)
    igds(18) = nint(trulon*1000)   ! longitude of south pole (millidegree)

    DO i = 19, 25
      igds(i) = 0       ! unused
    END DO

  ELSE IF (ABS(mapproj) == 3 .OR. ABS(mapproj) == 4) THEN  ! Mercator or Lat/Lon

    igds( 1) = nz       ! number of vertical (z) coordinates
    igds( 2) = 255      ! location (in octets) of vertical
                        ! coordinate list
    IF ( ABS(mapproj) == 4 ) THEN
      igds( 3) = 0      ! data representation type (table 6)
    ELSE
      igds( 3) = 1      ! data representation type (table 6)
    END IF
    igds( 4) = nx       ! number of x-coordinates
    igds( 5) = ny       ! number of y-coordinates
    igds( 6) = NINT(swlat*1000)    ! latitude  at southwest corner
    igds( 7) = NINT(swlon*1000)    ! longitude at southwest corner
    igds( 8) = 0        ! res-comp flag: 0: u & v relative to E/N
    igds( 9) = NINT(nelat*1000)    ! latitude  of last point
    igds(10) = NINT(nelon*1000)    ! longitude of last point
    igds(11) = ABS(NINT((nelat-swlat)*1000./(ny-1)))  ! latitude  increment
    igds(12) = ABS(NINT((nelon-swlon)*1000./(nx-1)))  ! longitude increment
    igds(13) = NINT(trulat1*1000)  ! true latitude
    igds(14) = 64       ! scanning flag, 64 for +i & +j direction

    DO i = 15,25
      igds(i) = 0       ! unused
    END DO

  ELSE
    WRITE (6,'(a,i2/a)')                                                &
        'ARPS does support the type of map projection: ',mapproj,       &
        'Model stopped in MKIGDS.'
    CALL arpsstop(' ',1)
  END IF

  RETURN
END SUBROUTINE mkigds
!
!##################################################################
!##################################################################
!######                                                      ######
!######               SUBROUTINE GRIBCNTL                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE gribcntl(nx,ny,nz,nzsoil,x,y,z,zp,zpsoil,gdsid,temxy1,temxy2)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  ARPS GRIB data description file for GrADS display
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yuhe Liu
!  11/27/1995
!
!  MODIFICATION HISTORY:
!
!  05/23/2002  (J. Brotzge) 
!  Added variables to allow for soil depth grid. 
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
!    x        The x-coord. of the physical and computational grid.
!             Defined at u-point.
!    y        The y-coord. of the physical and computational grid.
!             Defined at v-point.
!    z        The z-coord. of the computational grid.
!             Defined at w-point on the staggered grid.
!
!  WORK ARRAY:
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx,ny,nz          ! Number of grid points in 3 directions
  INTEGER :: nzsoil            ! Number of grid points in the soil


  REAL :: x     (nx)           ! The x-coord. of the physical and
                               ! computational grid. Defined at u-point.
  REAL :: y     (ny)           ! The y-coord. of the physical and
                               ! computational grid. Defined at v-point.
  REAL :: z     (nz)           ! The z-coord. of the computational grid.
                               ! Defined at w-point on the staggered grid.
  REAL :: zp    (nx,ny,nz)     ! The physical height coordinate defined at
                               ! w-point of the staggered grid.

  REAL :: zpsoil (nx,ny,nzsoil) ! The physical height coordinate defined at
                               ! w-point of the soil.

  INTEGER :: gdsid             ! ID number of GRIB grid.
                               ! (should always be 255)

  REAL :: temxy1(nx,ny)        ! Temporary array
  REAL :: temxy2(nx,ny)        ! Temporary array
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=256) :: grbctlfl
  CHARACTER (LEN=256) :: grbdatfl
  CHARACTER (LEN=15) :: chrstr, chr1
  CHARACTER (LEN=50) :: vartit(50), varnam(50)

  INTEGER :: varlev(50), varparam(50,4)

  INTEGER :: varnum
  INTEGER :: nchout0

  CHARACTER (LEN=3) :: monnam(12)            ! Name of months
  DATA monnam/'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun',                 &
              'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'/

  CHARACTER (LEN=2) :: dtunit
  INTEGER :: ntm,tinc

  REAL :: lonmin,latmin,lonmax,latmax
  REAL :: xbgn,ybgn,zbgn
  REAL :: xinc,yinc,z0
  REAL :: lat11,lat12,lon11,lon12,lat21,lat22,lon21,lon22
  REAL :: latinc,loninc

  INTEGER :: fnlen
  INTEGER :: i,k, is
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'grid.inc'          ! Grid & map parameters.
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
!  Open the GrADS data control file: runname.ctl
!
!-----------------------------------------------------------------------
!
  IF ( ABS(mapproj) /= 2 ) RETURN

  fnlen = lfnkey + 9
  grbctlfl(1:fnlen) = runname(1:lfnkey)//'.gribcntl'

  CALL fnversn( grbctlfl, fnlen )
  CALL getunit (nchout0)
  WRITE (6,'(a)') 'The GrADS data control file is '                     &
                    //grbctlfl(1:fnlen)

  OPEN (nchout0, FILE = grbctlfl(1:fnlen), STATUS = 'unknown')

  IF ( ldirnam == 0 ) THEN
    ldirnam = ldirnam + 1
    dirname(1:ldirnam) = '.'
  END IF

  fnlen = ldirnam + lfnkey + 14
  grbdatfl(1:fnlen) = dirname(1:ldirnam)//'/'                           &
                    //runname(1:lfnkey)//'.grb000000'

  xbgn = 0.5 * (x(1) + x(2))
  ybgn = 0.5 * (y(1) + y(2))
  zbgn = 0.5 * (z(1) + z(2))

  xinc = (x(2) - x(1))
  yinc = (y(2) - y(1))
  z0 = (z(2)-z(1))/2.

  CALL xytoll(nx,ny,x,y,temxy1,temxy2)

  CALL xytoll(1,1,xbgn,ybgn,lat11,lon11)

  CALL a3dmax0(temxy1,1,nx,1,nx,1,ny,1,ny-1,1,1,1,1,                    &
               latmax,latmin)
  CALL a3dmax0(temxy2,1,nx,1,nx,1,ny,1,ny-1,1,1,1,1,                    &
               lonmax,lonmin)

  latinc = (latmax-latmin)/(ny-1)
  loninc = (lonmax-lonmin)/(nx-1)

  WRITE (6,'(a,f10.4,a,f10.4,a,f10.4)')                                 &
           'latmin:latmax:latinc = ',                                   &
            latmin,':',latmax,':',latinc
  WRITE (6,'(a,f10.4,a,f10.4,a,f10.4)')                                 &
           'lonmin:lonmax:loninc = ',                                   &
           lonmin,':',lonmax,':',loninc

  WRITE (chrstr,'(i2.2,a,i2.2,a,i2.2,a3,i4.4)')                         &
      hour,':',minute,'Z',day,monnam(month),year

  IF ( nint(thisdmp) <= 0 ) THEN
    ntm = 1
  ELSE
    ntm = nint(tstop/thisdmp) + 1
  END IF

  IF ( nint(thisdmp) < 60 ) THEN
    WRITE (6, '(/a/a)')                                                 &
        'GrADS reqiures the smallest uint minute for time interval.',   &
        'Return to the caller.'
    tinc = nint(thisdmp)
    RETURN
  ELSE IF (thisdmp < 3600.) THEN
    tinc = nint(thisdmp/60.)
    dtunit = 'MN'
  ELSE IF (thisdmp < 86400.) THEN
    tinc = nint(thisdmp/3600.)
    dtunit = 'HR'
  ELSE
    tinc = nint(thisdmp/86400.)
    dtunit = 'DY'
  END IF

  varnum = 0

  IF ( varout == 1 ) THEN

    varnum = varnum + 1
    varnam(varnum) = 'u'
    vartit(varnum) = 'X-velocity total wind u (m/s)'
    varlev(varnum) = nz
    varparam(varnum,1) = 33
    varparam(varnum,2) = 103
    varparam(varnum,3) = z0

    varnum = varnum + 1
    varnam(varnum) = 'v'
    vartit(varnum) = 'Y-velocity total wind v (m/s)'
    varlev(varnum) = nz
    varparam(varnum,1) = 34
    varparam(varnum,2) = 103
    varparam(varnum,3) = z0

    varnum = varnum + 1
    varnam(varnum) = 'w'
    vartit(varnum) = 'Z-velocity total wind w (m/s)'
    varlev(varnum) = nz
    varparam(varnum,1) = 40
    varparam(varnum,2) = 103
    varparam(varnum,3) = z0

    varnum = varnum + 1
    varnam(varnum) = 'pt'
    vartit(varnum) = 'Potential Temperature (K)'
    varlev(varnum) = nz
    varparam(varnum,1) = 13
    varparam(varnum,2) = 103
    varparam(varnum,3) = z0

    varnum = varnum + 1
    varnam(varnum) = 'p'
    vartit(varnum) = 'Pressure (mb)'
    varlev(varnum) = nz
    varparam(varnum,1) = 1
    varparam(varnum,2) = 103
    varparam(varnum,3) = z0

  END IF

  IF ( mstout == 1 ) THEN

    varnum = varnum + 1
    varnam(varnum) = 'qv'
    vartit(varnum) = 'Water Vapor Mixing Ratio (g/kg)'
    varlev(varnum) = nz
    varparam(varnum,1) = 51
    varparam(varnum,2) = 103
    varparam(varnum,3) = z0

    varnum = varnum + 1
    varnam(varnum) = 'qc'
    vartit(varnum) = 'Cloud Water Mixing Ratio (g/kg)'
    varlev(varnum) = nz
    varparam(varnum,1) = 153
    varparam(varnum,2) = 103
    varparam(varnum,3) = z0

    varnum = varnum + 1
    varnam(varnum) = 'qr'
    vartit(varnum) = 'Rain Water Mixing Ratio (g/kg)'
    varlev(varnum) = nz
    varparam(varnum,1) = 170
    varparam(varnum,2) = 103
    varparam(varnum,3) = z0

    IF ( iceout == 1 ) THEN

      varnum = varnum + 1
      varnam(varnum) = 'qi'
      vartit(varnum) = 'Cloud Ice Mixing Ratio (g/kg)'
      varlev(varnum) = nz
      varparam(varnum,1) = 178
      varparam(varnum,2) = 103
      varparam(varnum,3) = z0

      varnum = varnum + 1
      varnam(varnum) = 'qs'
      vartit(varnum) = 'Snow Mixing Ratio (g/kg)'
      varlev(varnum) = nz
      varparam(varnum,1) = 171
      varparam(varnum,2) = 103
      varparam(varnum,3) = z0

      varnum = varnum + 1
      varnam(varnum) = 'qh'
      vartit(varnum) = 'Hail Mixing Ratio (g/kg)'
      varlev(varnum) = nz
      varparam(varnum,1) = 179
      varparam(varnum,2) = 103
      varparam(varnum,3) = z0

    END IF

    IF ( rainout == 1 ) THEN

      varnum = varnum + 1
      varnam(varnum) = 'raing'
      vartit(varnum) = 'Grid Supersaturation Rain '
      varlev(varnum) = 0
      varparam(varnum,1) = 62
      varparam(varnum,2) = 1
      varparam(varnum,3) = 0

      varnum = varnum + 1
      varnam(varnum) = 'rainc'
      vartit(varnum) = 'Cumulus Convection Rain '
      varlev(varnum) = 0
      varparam(varnum,1) = 63
      varparam(varnum,2) = 1
      varparam(varnum,3) = 0

    END IF

    IF ( prcout == 1 ) THEN

      varnum = varnum + 1
      varnam(varnum) = 'prcrt1'
      vartit(varnum) = 'Total precipitation rate '
      varlev(varnum) = 1
      varparam(varnum,1) = 59
      varparam(varnum,2) = 111
      varparam(varnum,3) = 1

      varnum = varnum + 1
      varnam(varnum) = 'prcrt2'
      vartit(varnum) = 'Grid scale precipitation rate '
      varlev(varnum) = 2
      varparam(varnum,1) = 59
      varparam(varnum,2) = 111
      varparam(varnum,3) = 2

      varnum = varnum + 1
      varnam(varnum) = 'prcrt3'
      vartit(varnum) = 'cumulus precipitation rate '
      varlev(varnum) = 3
      varparam(varnum,1) = 59
      varparam(varnum,2) = 111
      varparam(varnum,3) = 3

      varnum = varnum + 1
      varnam(varnum) = 'prcrt4'
      vartit(varnum) = 'Microphysics precipitation rate '
      varlev(varnum) = 4
      varparam(varnum,1) = 59
      varparam(varnum,2) = 111
      varparam(varnum,3) = 4

    END IF
  END IF

  IF ( tkeout == 1 ) THEN

    varnum = varnum + 1
    varnam(varnum) = 'tke'
    vartit(varnum) = 'Turbulent kinetic energy (m**2/s)'
    varlev(varnum) = nz
    varparam(varnum,1) = 158
    varparam(varnum,2) = 103
    varparam(varnum,3) = z0

  END IF

  IF ( trbout == 1 ) THEN

    varnum = varnum + 1
    varnam(varnum) = 'kmh'
    vartit(varnum) = 'Horiz. Turb. Mixing Coeff. for '                  &
                     //'Momentum (m**2/s)'
    varlev(varnum) = nz
    varparam(varnum,1) = 185
    varparam(varnum,2) = 103
    varparam(varnum,3) = z0

    varnum = varnum + 1
    varnam(varnum) = 'kmv'
    vartit(varnum) = 'Vertical Turb. Mixing Coeff. for '                &
                     //'Momentum (m**2/s)'
    varlev(varnum) = nz
    varparam(varnum,1) = 186
    varparam(varnum,2) = 103
    varparam(varnum,3) = z0

  END IF

  IF ( sfcout == 1 ) THEN

    IF ( nstyp <= 1 ) THEN
      varnum = varnum + 1
      varnam(varnum) = 'tsoil'
      vartit(varnum) = 'Soil Temperature (K)'
      varlev(varnum) = nzsoil  
      varparam(varnum,1) = 85
      varparam(varnum,2) = 111
      varparam(varnum,3) = 1

      varnum = varnum + 1
      varnam(varnum) = 'qsoil'
      vartit(varnum) = 'Soil moisture (m**3/m**3)'
      varlev(varnum) = nzsoil  
      varparam(varnum,1) = 144
      varparam(varnum,2) = 111
      varparam(varnum,3) = 1

      varnum = varnum + 1
      varnam(varnum) = 'wr'
      vartit(varnum) = 'Canopy Soil moisture'
      varlev(varnum) = 0
      varparam(varnum,1) = 223
      varparam(varnum,2) = 111
      varparam(varnum,3) = 0

    ELSE

      DO is=0,nstyp
        WRITE (chr1,'(i1)') is

        varnum = varnum + 1
        varnam(varnum) = 'tsoil'//chr1
        vartit(varnum) = 'Soil Temperature (K)'
        varlev(varnum) = nzsoil 
        varparam(varnum,1) = 85
        varparam(varnum,2) = 111
        varparam(varnum,3) = 1+is

        varnum = varnum + 1
        varnam(varnum) = 'qsoil'//chr1
        vartit(varnum) = 'Soil moisture (m**3/m**3)'
        varlev(varnum) = nzsoil  
        varparam(varnum,1) = 144
        varparam(varnum,2) = 111
        varparam(varnum,3) = 1+is

        varnum = varnum + 1
        varnam(varnum) = 'wr'//chr1
        vartit(varnum) = 'Canopy Soil moisture'
        varlev(varnum) = 0
        varparam(varnum,1) = 223
        varparam(varnum,2) = 111
        varparam(varnum,3) = is
      END DO

    END IF

    IF ( snowout == 1 ) THEN
      varnum = varnum + 1
      varnam(varnum) = 'snowd'
      vartit(varnum) = 'Snow depth (m)'
      varlev(varnum) = 0
      varparam(varnum,1) = 66
      varparam(varnum,2) = 1
      varparam(varnum,3) = 0
    END IF

  END IF

  IF( radout == 1 ) THEN
    varnum = varnum + 1
    varnam(varnum) = 'radfrc'
    vartit(varnum) = 'Radiation forcing (K/s)'
    varlev(varnum) = nz
    varparam(varnum,1) = 216
    varparam(varnum,2) = 103
    varparam(varnum,3) = z0

    varnum = varnum + 1
    varnam(varnum) = 'radsw'
    vartit(varnum) = 'Net short wave radiation flux (W/m**2)'
    varlev(varnum) = 0
    varparam(varnum,1) = 204
    varparam(varnum,2) = 1
    varparam(varnum,3) = 0

    varnum = varnum + 1
    varnam(varnum) = 'rnflx'
    vartit(varnum) = 'Total radiation flux (W/m**2)'
    varlev(varnum) = 0
    varparam(varnum,1) = 232
    varparam(varnum,2) = 1
    varparam(varnum,3) = 0
  END IF

  IF( flxout == 1 ) THEN
    varnum = varnum + 1
    varnam(varnum) = 'uflx'
    vartit(varnum) = 'u flux'
    varlev(varnum) = 0
    varparam(varnum,1) = 124
    varparam(varnum,2) = 1
    varparam(varnum,3) = 0

    varnum = varnum + 1
    varnam(varnum) = 'vflx'
    vartit(varnum) = 'v flux'
    varlev(varnum) = 0
    varparam(varnum,1) = 125
    varparam(varnum,2) = 1
    varparam(varnum,3) = 0

    varnum = varnum + 1
    varnam(varnum) = 'ptflx'
    vartit(varnum) = 'pt flux'
    varlev(varnum) = 0
    varparam(varnum,1) = 122
    varparam(varnum,2) = 1
    varparam(varnum,3) = 0

    varnum = varnum + 1
    varnam(varnum) = 'qvflx'
    vartit(varnum) = 'qv flux'
    varlev(varnum) = 0
    varparam(varnum,1) = 121
    varparam(varnum,2) = 1
    varparam(varnum,3) = 0
  END IF

  DO i=1,varnum
    varparam(i,4) = 0
  END DO

  WRITE (nchout0,'(a/a)')                                               &
      'TITLE   ARPS 4.0 GRIB Control for GrADS Display of '             &
      //runname(1:lfnkey)//' (a sample)','*'

  WRITE (nchout0,'(a,a)')                                               &
      'DSET    ',grbdatfl(1:fnlen)

  WRITE (nchout0,'(a,i4)')                                              &
      'DTYPE   GRIB ',gdsid

  WRITE (nchout0,'(a)')                                                 &
      'OPTIONS template'

  WRITE (nchout0,'(a)')                                                 &
      'INDEX   '//runname(1:lfnkey)//'.gribmap'

  WRITE (nchout0,'(a/a)')                                               &
      'UNDEF   -9.e+33','*'

  IF ( mapproj == 2 ) THEN

    WRITE (nchout0,'(a)')                                               &
        '* For lat-lon-lev display, umcomment the following 4 lines.'

    WRITE (nchout0,'(a,1x,i8,1x,i3,a,2f12.6,2i3,3f12.6,2f12.2)')        &
        'PDEF',nx,ny,' LCC',lat11,lon11,1,1,                            &
              trulat1,trulat2,trulon,xinc,yinc

  END IF

  WRITE (nchout0,'(a,1x,i8,a,f10.4,1x,f10.4)')                          &
      'XDEF',nx, '  LINEAR  ',lonmin,loninc

  WRITE (nchout0,'(a,1x,i8,a,f10.4,1x,f10.4)')                          &
      'YDEF',ny, '  LINEAR  ',latmin,latinc

  WRITE (nchout0,'(a,1x,i8,a)')                                         &
      'ZDEF',nz,'  LEVELS  '
  WRITE (nchout0,'(8f10.2)')                                            &
      ((z(k)+z(k+1))/2.,k=1,nz-1),z(nz)
  WRITE (nchout0,'(a/a/a)')                                             &
  '* WARNING& ! The vertical levels were set to computational ',        &
  '*          coordinates, z(k), because zp is not uniform in ',        &
  '*          horizontal when terrain option was turned on'

!  WRITE (nchout0,'(a,1x,i8,a)')                                         &
!      'ZSOILDEF',nzsoil,'  LEVELS  '
!  WRITE (nchout0,'(8f10.2)')                                            &
!      ((zpsoil(k)+zpsoil(k+1))/2.,k=1,nzsoil-1),zpsoil(nzsoil)
!  WRITE (nchout0,'(a/a/a)')                                             &
!  '* WARNING& ! The vertical levels were set to computational ',        &
!  '*          coordinates, zpsoil(k), because zpsoil is not uniform in ',  &
!  '*          horizontal when terrain option was turned on'


  WRITE (nchout0,'(a,1x,i8,a,a,3x,i2.2,a/a)')                           &
      'TDEF',ntm,'  LINEAR  ',chrstr,tinc,dtunit,'*'

  WRITE (nchout0,'(a,1x,i3)')                                           &
      'VARS',varnum

  DO i = 1, varnum
    WRITE (nchout0,'(a6,1x,i3,2x,i3.3,a,i4.4,a,i3.3,2x,a)')             &
        varnam(i),varlev(i), varparam(i,1),',',varparam(i,2),',',       &
        varparam(i,3),vartit(i)
  END DO

  WRITE (nchout0,'(a)')                                                 &
      'ENDVARS'

  CLOSE (nchout0)
  CALL retunit(nchout0)

  RETURN
END SUBROUTINE gribcntl
