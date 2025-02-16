!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE SVIDUMP                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE svidump(nx,ny,nz,nzsoil,nstyps, nchout,graffn, grdbas,       &
           u,v,w,ptprt,pprt,qv,qc,qr,qi,qs,qh,tke,kmh,kmv,              &
           ubar,vbar,wbar,ptbar,pbar,rhobar,qvbar,                      &
           x,y,z,zp,zpsoil,                                             &
           soiltyp,stypfrct,vegtyp,lai,roufns,veg,                      &
           tsoil,qsoil,wetcanp,snowdpth,                                &
           raing,rainc,prcrate,                                         &
           radfrc,radsw,rnflx,radswnet,radlwin,                         &
           usflx,vsflx,ptsflx,qvsflx,                                   &
           tem1, tem2)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Dump a data file for the visualization program Savi3D.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Jason J. Levit
!  07/20/92
!
!  MODIFICATION HISTORY:
!
!  08/02/92 (J. Levit)
!  Added full documentation and performed a clean-up.
!
!  08/04/92 (M. Xue)
!  Subroutine streamlined to conform to data dump format standards.
!
!  8/23/92 (M. Xue)
!  Modify to perform the dumping of both base and t-dependent arrays
!  and added control on grid staggering.
!
!  9/18/92 (J. Levit, M. Xue)
!  Added code to produce a configuration file for Savi3D.
!
!  11/2/92 (M. Xue)
!
!  Major overhaul. grafwritescalarpoint is called rather than
!  grafwritesclararray. This elliminates the need to define the
!  grid work array and make the data writing more flexible.
!
!  The capability to write out part of a data array implemented.
!
!  09/02/94 (J. Levit & Y. Lu)
!  Cleaned up documentation.
!
!  11/10/94 (Liping Sun & Min Zou)
!  Upgraded to version 1.2.2 of Savi3D. The file format was changed
!  from GRAF to MeRAF.
!
!  11/10/94 (Y. Liu)
!  Merged the upgraded version with the documentation cleaned up
!  version.
!
!  12/09/1998 (Donghai Wang)
!  Added the snow cover.
!
!  05/23/2002 (J. Brotzge)
!  Added additional variables for soil levels and new soil schemes.  
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
!    graffn   Name of the Savi3D MeRAF file.
!    grdbas   Flag indicating if this is a call for the data dump
!             of grid and base state arrays only. If so, grdbas=1
!             (not used in this routine).
!
!    u        x component of velocity at a given time level (m/s)
!    v        y component of velocity at a given time level (m/s)
!    w        Vertical component of Cartesian velocity at a given
!             time level (m/s)
!    ptprt    Perturbation potential temperature at a given time
!             level (K)
!    pprt     Perturbation pressure at  a given time level (Pascal)
!    qv       Water vapor specific humidity at a given time level
!             (kg/kg)
!    qc       Cloud water mixing ratio at a given time level (kg/kg)
!    qr       Rainwater mixing ratio at a given time level (kg/kg)
!    qi       Cloud ice mixing ratio at a given time level (kg/kg)
!    qs       Snow mixing ratio at a given time level (kg/kg)
!    qh       Hail mixing ratio at a given time level (kg/kg)
!
!    tke      Turbulent Kinetic Energy ((m/s)**2)
!    kmh      Horizontal turb. mixing coef. for momentum ( m**2/s )
!    kmv      Vertical turb. mixing coef. for momentum ( m**2/s )
!
!    ubar     Base state zonal velocity component (m/s)
!    vbar     Base state meridional velocity component (m/s)
!    wbar     Base state vertial velocity component (m/s)
!    ptbar    Base state potential temperature (K)
!    pbar     Base state pressure (Pascal)
!    rhobar   Base state density (kg/m**3)
!    qvbar    Base state water vapor specific humidity (kg/kg)
!
!    x        x coordinate of grid points in physical/comp. space (m)
!    y        y coordinate of grid points in physical/comp. space (m)
!    z        z coordinate of grid points in computational space (m)
!    zp       Vertical coordinate of grid points in physical space(m)
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
!    radswnet Net shortwave radiation
!    radlwin  Incoming longwave radiation
!
!    usflx    Surface flux of u-momentum (kg/(m*s**2))
!    vsflx    Surface flux of v-momentum (kg/(m*s**2))
!    ptsflx   Surface heat flux (K*kg/(m**2 * s ))
!    qvsflx   Surface moisture flux of (kg/(m**2 * s))
!
!    grafgrid Passed from a dummy variable, used to define the grid
!             in Savi3D
!
!  OUTPUT:
!
!    None.
!
!  WORK ARRAY:
!
!    tem1     Temporary work array.
!    tem2     Temporary work array.
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
!
  INTEGER :: nx,ny,nz          ! Number of grid points in 3 directions
  INTEGER :: nzsoil            ! Number of grid points in the soil  
!
  CHARACTER (LEN=*     ) :: graffn ! Name of the Savi3D MeRAF file
  INTEGER :: grdbas            ! If this is a grid/base state dump
!
  REAL :: u     (nx,ny,nz)     ! Total u-velocity (m/s)
  REAL :: v     (nx,ny,nz)     ! Total v-velocity (m/s)
  REAL :: w     (nx,ny,nz)     ! Total w-velocity (m/s)
  REAL :: ptprt (nx,ny,nz)     ! Perturbation potential temperature (K)
  REAL :: pprt  (nx,ny,nz)     ! Perturbation pressure (Pascal)
!
  REAL :: qv    (nx,ny,nz)     ! Water vapor specific humidity (kg/kg)
  REAL :: qc    (nx,ny,nz)     ! Cloud water mixing ratio (kg/kg)
  REAL :: qr    (nx,ny,nz)     ! Rain water mixing ratio (kg/kg)
  REAL :: qi    (nx,ny,nz)     ! Cloud ice mixing ratio (kg/kg)
  REAL :: qs    (nx,ny,nz)     ! Snow mixing ratio (kg/kg)
  REAL :: qh    (nx,ny,nz)     ! Hail mixing ratio (kg/kg)
  REAL :: tke   (nx,ny,nz)     ! Turbulent Kinetic Energy ((m/s)**2)
!
  REAL :: kmh   (nx,ny,nz)     ! Horizontal turb. mixing coef. for
                               ! momentum. ( m**2/s )
  REAL :: kmv   (nx,ny,nz)     ! Vertical turb. mixing coef. for
                               ! momentum. ( m**2/s )
!
  REAL :: ubar  (nx,ny,nz)     ! Base state u-velocity (m/s)
  REAL :: vbar  (nx,ny,nz)     ! Base state v-velocity (m/s)
  REAL :: wbar  (nx,ny,nz)     ! Base state w-velocity (m/s)
  REAL :: ptbar (nx,ny,nz)     ! Base state potential temperature (K)
  REAL :: pbar  (nx,ny,nz)     ! Base state pressure (Pascal)
  REAL :: rhobar(nx,ny,nz)     ! Base state air density (kg/m**3)
  REAL :: qvbar (nx,ny,nz)     ! Base state water vapor specific
                               ! humidity (kg/kg)

  REAL :: x     (nx)           ! The x-coord. of the physical and
                               ! computational grid. Defined at u-point.
  REAL :: y     (ny)           ! The y-coord. of the physical and
                               ! computational grid. Defined at v-point.
  REAL :: z     (nz)           ! The z-coord. of the computational grid.
                               ! Defined at w-point on the staggered grid.
  REAL :: zp    (nx,ny,nz)     ! The physical height coordinate defined
                               ! at w-point of the staggered grid.
  REAL :: zpsoil (nx,ny,nzsoil) ! The physical height coordinate defined
                               ! at w-point of the soil.  

  INTEGER :: nstyps
  INTEGER :: soiltyp (nx,ny,nstyps)    ! Soil type
  REAL :: stypfrct(nx,ny,nstyps)       ! Soil type
  INTEGER :: vegtyp(nx,ny)             ! Vegetation type
  REAL :: lai    (nx,ny)            ! Leaf Area Index
  REAL :: roufns (nx,ny)            ! Surface roughness
  REAL :: veg    (nx,ny)            ! Vegetation fraction

  REAL :: tsoil (nx,ny,nzsoil,0:nstyps) ! Soil temperature (K)
  REAL :: qsoil (nx,ny,nzsoil,0:nstyps) ! Soil moisture (m**3/m**3) 
  REAL :: wetcanp(nx,ny,0:nstyps)       ! Canopy water amount
  REAL :: snowdpth(nx,ny)               ! Snow depth (m)

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

  REAL :: tem1  (nx*ny*nz)     ! Temporary work array
  REAL :: tem2  (nx,ny,nz)     ! Temporary work array
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!

  INCLUDE 'meraf.inc'
  INCLUDE 'globcst.inc'
!  include 'grafibm.inc'
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: nxout,nyout,nzout ! The size of array to be written out.
  INTEGER :: nzsoilout ! The size of array to be written out.

  INTEGER :: ist ,ind ,isk ,jst ,jnd ,jsk ,kst ,knd ,ksk
  INTEGER :: ist1,ind1,isk1,jst1,jnd1,jsk1,kst1,knd1,ksk1
  INTEGER :: lst, lnd, lsk, lst1, lnd1, lsk1 

  INTEGER :: uid, vid, wid, uprtid, vprtid, wprtid        !*****
  INTEGER :: qvprtid, qcid, qrid, qwid, qiid, qsid, qhid  !*****
  INTEGER :: vortid, divid, ubarid, vbarid, wbarid        !*****
  INTEGER :: pbarid, rhobarid, qvbarid   !*****
  INTEGER :: windid, totalwindid         !*****
  INTEGER :: ptprtid, ptbarid, pprtid    !*****
  INTEGER :: value             ! *****
  INTEGER :: frame             ! *****

  INTEGER :: ierr              ! Used as an int'l error code by Savi3D.
  INTEGER :: i,j,k,l             ! Used by do loops.
  INTEGER :: ii,jj,kk,ll 
  INTEGER :: nchout            ! Unused.
  INTEGER :: ishf, jshf, kshf, lshf 
  CHARACTER (LEN=50) :: configname   ! Used to create Savi3D config file.

  CHARACTER (LEN=20) :: schemename   !
  CHARACTER (LEN=40) :: errorstring  !
  REAL*8 xbase,ybase,zbase  ! ****

  REAL :: conx, cony       ! Used to create Savi3D config file.
  INTEGER :: gbwrtn            ! See if grid and base state
                               ! parameter/arrays have been written
                               ! into the data file
  INTEGER :: ncalls
  DATA gbwrtn,ncalls /0,0/
  SAVE gbwrtn,ncalls
  CHARACTER (LEN=7) :: chtem2
  CHARACTER (LEN=7) :: chtem1
                                ! Used to create Savi3D config file.
  CHARACTER (LEN=6) :: timhms

!  integer year,month,day,hour,minute,second,node ! def. in globcst.inc
  INTEGER :: node
  REAL :: tem
  REAL :: xcord,ycord,zcord
  REAL*8 second1            ! Used for Savi3D
  INTEGER :: nchout0         ! Used to open Savi3D config file.

  INTEGER :: setdomn,setskip
  SAVE setdomn, setskip
  SAVE ist,ind,isk,jst,jnd,jsk,kst,knd,ksk,lst,lnd,lsk 
  DATA setdomn/0/, setskip /0/
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  IF( setdomn == 0 ) THEN
    ist = 1
    ind = nx-1
    jst = 1
    jnd = ny-1
    kst = 1
    knd = nz-1
    lst = 1
    lnd = nzsoil-1
  END IF

  IF( setskip == 0 ) THEN
    isk = 1
    jsk = 1
    ksk = 1
    lsk = 1 
  END IF

  WRITE(6,'(/1x,a,f8.1,a/)')                                            &
       'Writing Savi3D data set at time=',curtim,' s.'
!
!-----------------------------------------------------------------------
!
!  The dimension of the array to be written out:
!
!-----------------------------------------------------------------------
!
  nxout = (ind-ist)/isk +1
  nyout = (jnd-jst)/jsk +1
  nzout = (knd-kst)/ksk +1
  nzsoilout = (lnd-lst)/lsk +1

  IF( ncalls == 0) THEN
!
!-----------------------------------------------------------------------
!
!  Create the Savi3D configuration file.
!
!-----------------------------------------------------------------------
!
    CALL getunit( nchout0 )
    CALL gtlfnkey(runname, lfnkey)

    configname=runname(1:lfnkey)//'.sviconfig'
    OPEN (UNIT=nchout0,                                                 &
          FILE=configname(1:10+lfnkey),STATUS='unknown',                &
          FORM='formatted')

    cony=(y(ind)-y(ist))/111111.0
    conx=(x(jnd)-x(jst))/111111.0

    WRITE(chtem1,'(f7.3)') cony
    WRITE(chtem2,'(f7.3)') conx

    DO i=1,7
      IF ( chtem1(i:i) == ' ') chtem1(i:i) = '0'
      IF ( chtem2(i:i) == ' ') chtem2(i:i) = '0'
    END DO

    WRITE (nchout0,'(4a)')                                              &
        'toplat=',chtem1,',btmlat=000.000,lftlon=000.000,rhtlon=',chtem2
    WRITE (nchout0,'(a,a)') 'MeRAF_gridded=',graffn
    WRITE (nchout0,'(a)') 'initzstretch=1'

    WRITE (6,'(a,a)') 'MeRAF_gridded=',graffn
    WRITE (6,'(a)') 'initzstretch=1'

    CLOSE (nchout0)
    CALL retunit( nchout0 )

    ncalls = 1
  END IF

  ishf = 1
  jshf = 1
  kshf = 1
  lshf = 1 
!
!-----------------------------------------------------------------------
!
!  Setup for the MeRAF interface.  This is only done at time=0,
!  and calls are made to Savi3D subroutines to define the type of
!  grid which the data are being written to.
!
!  Also, the scalar variables are defined for Savi3D here.
!
!-----------------------------------------------------------------------
!
  IF ( gbwrtn == 0 ) THEN

    PRINT*,' Opening Savi3D file ',graffn

!
!-----------------------------------------------------------------------
!
!  Create the data set
!
!-----------------------------------------------------------------------
!

    CALL mcreatedataset (graffn,'ARPS 4.0',true,dsindex,ierr)
    IF ( ierr == meerr ) THEN
      CALL mgeterror (errorstring)
      PRINT*, 'Error: Unable to create data set'
      PRINT*, 'Error Message: ',errorstring
    END IF

!
!-----------------------------------------------------------------------
!
!  Create a Grid
!
!-----------------------------------------------------------------------
!

    CALL mcreategrid(dsindex,'ARPS 4.0',' ','3-D Model',                &
                     3, nxout, nyout, nzout, ' ',' ', gridid, ierr)
    IF ( ierr == meerr ) THEN
      CALL mgeterror (errorstring)
      PRINT*, 'Error: Unable to create grid scheme'
      PRINT*, 'Error Message: ',errorstring
    END IF

!
!-----------------------------------------------------------------------
!
!  Set up Grid
!
!-----------------------------------------------------------------------
!

    CALL mconfigurelocations(gridid, me_non_time_var,                   &
                             mecartesian, ' ', ierr)
    IF ( ierr == meerr ) THEN
      CALL mgeterror (errorstring)
      PRINT*, 'Error: Unable to configue locations'
      PRINT*, 'Error Message: ',errorstring
    END IF

    frame=0
    xbase=0.0D0
    ybase=0.0D0
    zbase=0.0D0

    CALL msetxyzbase (gridid, xbase, ybase, zbase, ierr)
    IF ( ierr == meerr ) THEN
      CALL mgeterror (errorstring)
      PRINT*, 'Error: Unable to set base point'
      PRINT*, 'Error Message: ',errorstring
    END IF

!
!-----------------------------------------------------------------------
!
!   Define Scalars
!
!-----------------------------------------------------------------------
!

    IF (varout == 1) THEN

      CALL mdefinescalar (gridid,'u','X-velocity total wind',           &
           'm/s', 'X-velocity total wind',                              &
           me_time_var, uid)
      IF ( uid == meerr ) THEN
        CALL mgeterror (errorstring)
        PRINT*, 'Error: Unable to define scalar', 'u'
        PRINT*, 'Error Message: ',errorstring
      END IF

      CALL mdefinescalar(gridid,'v','Y-velocity total wind',            &
            'm/s', 'Y-velocity total wind',                             &
            me_time_var, vid)
      IF ( vid == meerr ) THEN
        CALL mgeterror (errorstring)
        PRINT*, 'Error: Unable to define scalar', 'v'
        PRINT*, 'Error Message: ',errorstring
      END IF

      CALL mdefinescalar(gridid,'w','Z-velocity total wind',            &
            'm/s', 'Z-velocity total wind',                             &
            me_time_var, wid)
      IF ( wid == meerr ) THEN
        CALL mgeterror (errorstring)
        PRINT*, 'Error: Unable to define scalar', 'w'
        PRINT*, 'Error Message: ',errorstring
      END IF

      CALL mdefinescalar(gridid,'uprt','X-velocity perturbation',       &
            'm/s', 'X-velocity perturbation',                           &
            me_time_var, uprtid)
      IF ( uprtid == meerr ) THEN
        CALL mgeterror (errorstring)
        PRINT*, 'Error: Unable to define scalar', 'uprt'
        PRINT*, 'Error Message: ',errorstring
      END IF

      CALL mdefinescalar(gridid,'vprt','Y-velocity perturbation',       &
            'm/s', 'Y-velocity perturbation',                           &
            me_time_var, vprtid)
      IF ( vprtid == meerr ) THEN
        CALL mgeterror (errorstring)
        PRINT*, 'Error: Unable to define scalar', 'vprt'
        PRINT*, 'Error Message: ',errorstring
      END IF

      CALL mdefinescalar(gridid,'wprt','Z-velocity perturbation',       &
            'm/s', 'Z-velocity perturbation',                           &
            me_time_var, wprtid)
      IF ( wprtid == meerr ) THEN
        CALL mgeterror (errorstring)
        PRINT*, 'Error: Unable to define scalar', 'wprt'
        PRINT*, 'Error Message: ',errorstring
      END IF

      CALL mdefinescalar (gridid,'ptprt',                               &
           'Perturbation Potential Temperature','deg K',                &
           'Perturbation Potential Temperature',                        &
           me_time_var, ptprtid)
      IF ( ptprtid == meerr ) THEN
        CALL mgeterror (errorstring)
        PRINT*, 'Error: Unable to define scalar', 'ptprt'
        PRINT*, 'Error Message: ',errorstring
      END IF

      CALL mdefinescalar (gridid,'pprt',                                &
            'Perturbation Pressure','mb',                               &
            'Perturbation Pressure',                                    &
            me_time_var, pprtid)
      IF ( pprtid == meerr ) THEN
        CALL mgeterror (errorstring)
        PRINT*, 'Error: Unable to define scalar', 'pprt'
        PRINT*, 'Error Message: ',errorstring
      END IF

    END IF

    IF (mstout == 1) THEN

      CALL mdefinescalar (gridid,'qvprt',                               &
           'Water Vapor Mixing Ratio Perturbation','g/kg',              &
           'Water Vapor Mixing Ratio Perturbation',                     &
           me_time_var, qvprtid)
      IF ( qvprtid == meerr ) THEN
        CALL mgeterror (errorstring)
        PRINT*, 'Error: Unable to define scalar', 'qvprt'
        PRINT*, 'Error Message: ',errorstring
      END IF

      CALL mdefinescalar (gridid,'qc',                                  &
           'Cloud Water Mixing Ratio','g/kg',                           &
           'Cloud Water Mixing Ratio',                                  &
           me_time_var, qcid)
      IF ( qcid == meerr ) THEN
        CALL mgeterror (errorstring)
        PRINT*, 'Error: Unable to define scalar', 'qc'
        PRINT*, 'Error Message: ',errorstring
      END IF

      CALL mdefinescalar (gridid,'qr',                                  &
           'Rain Water Mixing Ratio','g/kg',                            &
           'Rain Water Mixing Ratio',                                   &
           me_time_var, qrid)
      IF ( qrid == meerr ) THEN
        CALL mgeterror (errorstring)
        PRINT*, 'Error: Unable to define scalar', 'qr'
        PRINT*, 'Error Message: ',errorstring
      END IF

      CALL mdefinescalar (gridid,'qw',                                  &
           'Total Water Mixing Ratio','g/kg',                           &
           'Total Water Mixing Ratio',                                  &
           me_time_var, qwid)
      IF ( qwid == meerr ) THEN
        CALL mgeterror (errorstring)
        PRINT*, 'Error: Unable to define scalar', 'qw'
        PRINT*, 'Error Message: ',errorstring
      END IF

      IF (iceout == 1) THEN

        CALL mdefinescalar (gridid,'qi',                                &
             'Cloud Ice Mixing Ratio','g/kg',                           &
             'Cloud Ice Mixing Ratio',                                  &
             me_time_var, qiid)
        IF ( qiid == meerr ) THEN
          CALL mgeterror (errorstring)
          PRINT*, 'Error: Unable to define scalar', 'qi'
          PRINT*, 'Error Message: ',errorstring
        END IF

        CALL mdefinescalar (gridid,'qs',                                &
             'Snow Mixing Ratio','g/kg',                                &
             'Snow Mixing Ratio',                                       &
             me_time_var, qsid)
        IF (qsid == meerr) THEN
          CALL mgeterror (errorstring)
          PRINT*, 'Error: Unable to define scalar','qs'
          PRINT*, 'Error Message: ',errorstring
        END IF

        CALL mdefinescalar (gridid,'qh',                                &
             'Hail Mixing Ratio','g/kg',                                &
             'Hail Mixing Ratio',                                       &
             me_time_var, qhid)
        IF (qhid == meerr) THEN
          CALL mgeterror (errorstring)
          PRINT*, 'Error: Unable to define scalar','qh'
          PRINT*, 'Error Message: ',errorstring
        END IF

      END IF

    END IF

    CALL mdefinescalar (gridid,'Vort','Vertical vorticity',             &
          '1/s', 'Vertical vorticity',                                  &
         me_time_var, vortid)
    IF (vortid == meerr) THEN
      CALL mgeterror (errorstring)
      PRINT*, 'Error: Unable to define scalar','Vort'
      PRINT*, 'Error Message: ',errorstring
    END IF

    CALL mdefinescalar (gridid,'Div','Horizontal divergence',           &
          '1/s', 'Horizontal divergence',                               &
         me_time_var, divid)
    IF (divid == meerr) THEN
      CALL mgeterror (errorstring)
      PRINT*, 'Error: Unable to define scalar','Div'
      PRINT*, 'Error Message: ',errorstring
    END IF

!
!-----------------------------------------------------------------------
!
!  The definitions for the time invariant scalars for Savi3D.
!  They are only written once to the Savi3D MeRAF file (at time=0).
!
!-----------------------------------------------------------------------
!

    IF (basout == 1) THEN

      CALL mdefinescalar (gridid,'ubar',                                &
           'Base State X Wind Velocity ','m/s',                         &
           'Base State X Wind Velocity',                                &
           me_time_var, ubarid)
      IF (ubarid == meerr) THEN
        CALL mgeterror (errorstring)
        PRINT*, 'Error: Unable to define scalar','ubar'
        PRINT*, 'Error Message: ',errorstring
      END IF

      CALL mdefinescalar (gridid,'vbar',                                &
           'Base State Y Wind Velocity ','m/s',                         &
           'Base State Y Wind Velocity',                                &
           me_time_var, vbarid)
      IF (vbarid == meerr) THEN
        CALL mgeterror (errorstring)
        PRINT*, 'Error: Unable to define scalar','vbar'
        PRINT*, 'Error Message: ',errorstring
      END IF

      CALL mdefinescalar (gridid,'wbar',                                &
           'Base State Z Wind Velocity ','m/s',                         &
           'Base State Z Wind Velocity',                                &
           me_time_var, wbarid)
      IF (wbarid == meerr) THEN
        CALL mgeterror (errorstring)
        PRINT*, 'Error: Unable to define scalar','wbar'
        PRINT*, 'Error Message: ',errorstring
      END IF

      CALL mdefinescalar (gridid,'ptbar',                               &
           'Base State Potential Temperature ','deg K',                 &
           'Base State Potential Temperature',                          &
            me_time_var, ptbarid)
      IF (ptbarid == meerr) THEN
        CALL mgeterror (errorstring)
        PRINT*, 'Error: Unable to define scalar','ptbar'
        PRINT*, 'Error Message: ',errorstring
      END IF

      CALL mdefinescalar (gridid,'pbar',                                &
           'Base State Pressure ','mb',                                 &
           'Base State Pressure',                                       &
           me_time_var, pbarid)
      IF (pbarid == meerr) THEN
        CALL mgeterror (errorstring)
        PRINT*, 'Error: Unable to define scalar','pbar'
        PRINT*, 'Error Message: ',errorstring
      END IF

      CALL mdefinescalar (gridid,'rhobar',                              &
           'Base State Air Density ','kg/m**3',                         &
           'Base State Air Density',                                    &
           me_time_var, rhobarid)
      IF (rhobarid == meerr) THEN
        CALL mgeterror (errorstring)
        PRINT*, 'Error: Unable to define scalar','rhobar'
        PRINT*, 'Error Message: ',errorstring
      END IF

      IF (mstout == 1) THEN

        CALL mdefinescalar (gridid,'qvbar',                             &
               'Base State Water Vapor Mixing Ratio ','g/kg',           &
               'Base State Water Vapor Mixing Ratio',                   &
               me_time_var, qvbarid)
        IF (qvbarid == meerr) THEN
          CALL mgeterror (errorstring)
          PRINT*, 'Error: Unable to define scalar','qvbar'
          PRINT*, 'Error Message: ',errorstring
        END IF

      END IF

    END IF

!
!-----------------------------------------------------------------------
!
!  Vector definitions for Savi3D.  Also, define the starting date,
!  time, and time step for Savi3D.
!
!-----------------------------------------------------------------------
!

!
!-----------------------------------------------------------------------
!
!  Define Vectors
!
!-----------------------------------------------------------------------
!

    IF (varout == 1) THEN

      CALL mdefinevectori(gridid,'Pertubation_Wind',                    &
          'Pertubation_Wind',                                           &
          uprtid, vprtid, wprtid, windid)
      IF (windid == meerr) THEN
        CALL mgeterror (errorstring)
        PRINT*, 'Error: Unable to define vector','Wind'
        PRINT*, 'Error Message: ',errorstring
      END IF

      CALL mdefinevectori(gridid,'Total_Wind','Total_Wind',             &
                          uid, vid, wid, totalwindid)
      IF (totalwindid == meerr) THEN
        CALL mgeterror (errorstring)
        PRINT*, 'Error: Unable to define vector','Total_Wind'
        PRINT*, 'Error Message: ',errorstring
      END IF

    END IF

  END IF

!
!-----------------------------------------------------------------------
!
!  Write buffers
!
!-----------------------------------------------------------------------
!

  CALL msetbufferwrite (gridid, ierr)
  IF (ierr == meerr) THEN
    CALL mgeterror (errorstring)
    PRINT*, 'Error: Unable to set buffer write mode'
    PRINT*, 'Error Message: ',errorstring
  END IF

!
!-----------------------------------------------------------------------
!
!   Set Grid Point Locations
!
!-----------------------------------------------------------------------
!

  node = 0
  DO k=kst,knd,ksk
    DO j=jst,jnd,jsk
      DO i=ist,ind,isk
        ii=(i-ist)/isk+1
        jj=(j-jst)/jsk+1
        kk=(k-kst)/ksk+1
        xcord =(x(i)+x(i+ishf)-x(ist)-x(ist+ishf) )
        ycord =(y(j)+y(j+jshf)-y(jst)-y(jst+jshf) )
        zcord =(zp(i,j,k)+zp (i,j,k+kshf))
        CALL msetxyzlocation( gridid, ii, jj, kk, xcord,ycord,          &
                              zcord, ierr)
        IF( ierr == meerr) THEN
          CALL mgeterror(errorstring)
          PRINT*,'Error: Unable to set Locations'
          PRINT*,'Error Message: ',errorstring
        END IF
        node=node+1
      END DO
    END DO
  END DO

  IF ( node /= nxout*nyout*nzout ) THEN
    WRITE(6,'(1x,a,a)') 'nxout*nyout*nzout value incorrect.',           &
                        ' Job stopped in SVIDUMP.'
    CALL arpsstop(' ',1)
  END IF

!-----------------------------------------------------------------------
!
!  Write the base state time invariant scalars to the Savi3D file.
!  Data is written the first time this routine is called.
!
!-----------------------------------------------------------------------
!

  CALL mstartframew (gridid, ierr)
  IF ( ierr == meerr ) THEN
    CALL mgeterror (errorstring)
    PRINT*, 'Error: Unable to start frame write'
    PRINT*, 'Error Message: ',errorstring
  END IF

  frame=frame+1
  PRINT*, 'frame=', frame

  year  = 1994
  month = 7
  day   = 1

  CALL cvttim(curtim, timhms)
  READ(timhms,'(3i2)') hour, minute,second
  second1=second
  PRINT*,'year=',year,',month=',month,',day=',day
  PRINT*,'hour=',hour,',minute=',minute,',second=',second1

!  CALL GrafTimeStart(grafh,year,month,day,hour,minute,second)

  CALL msettimestampu( gridid, year, month, day, hour, minute,          &
                       second1, 0, 0, ierr)
  IF ( ierr == meerr ) THEN
    CALL mgeterror (errorstring)
    PRINT*, 'Error: Unable to set time stamp'
    PRINT*, 'Error Message: ',errorstring
  END IF

  IF ( gbwrtn == 0 ) THEN

    IF (basout == 1 ) THEN

      node = 0
      DO k=kst,knd,ksk
        DO j=jst,jnd,jsk
          DO i=ist,ind,isk
            node = node + 1
            tem1(node) = ((ubar(i,j,k)+ubar(i+ishf,j,k))*0.5)
          END DO
        END DO
      END DO

      CALL mwritescalararrayi(gridid,ubarid,1,nxout,                    &
             1,nyout,1,nzout,tem1,ierr)
      IF (ierr == meerr) THEN
        CALL mgeterror (errorstring)
        PRINT*, 'Error: Unable to write scalar. ubarID=', ubarid
        PRINT*, 'Error Message: ',errorstring
      END IF

      node = 0
      DO k=kst,knd,ksk
        DO j=jst,jnd,jsk
          DO i=ist,ind,isk
            node = node + 1
            tem1(node) = ((vbar(i,j,k)+vbar(i+ishf,j,k))*0.5)
          END DO
        END DO
      END DO

      CALL mwritescalararrayi(gridid,vbarid,1,nxout,                    &
              1,nyout,1,nzout,tem1,ierr)
      IF (ierr == meerr) THEN
        CALL mgeterror (errorstring)
        PRINT*, 'Error: Unable to write scalar. vbarID=', vbarid
        PRINT*, 'Error Message: ',errorstring
      END IF

      node = 0
      DO k=kst,knd,ksk
        DO j=jst,jnd,jsk
          DO i=ist,ind,isk
            node = node + 1
            tem1(node) = ((wbar(i,j,k)+wbar(i+ishf,j,k))*0.5)
          END DO
        END DO
      END DO

      CALL mwritescalararrayi(gridid,wbarid,1,nxout,                    &
                              1,nyout,1,nzout,tem1,ierr)
      IF (ierr == meerr) THEN
        CALL mgeterror (errorstring)
        PRINT*, 'Error: Unable to write scalar. wbarID=', wbarid
        PRINT*, 'Error Message: ',errorstring
      END IF

      node = 0
      DO k=kst,knd,ksk
        DO j=jst,jnd,jsk
          DO i=ist,ind,isk
            node = node + 1
            tem1(node) = ptbar(i,j,k)
          END DO
        END DO
      END DO

      CALL mwritescalararrayi(gridid,ptbarid,1,nxout,                   &
                              1,nyout,1,nzout,tem1,ierr)
      IF (ierr == meerr) THEN
        CALL mgeterror (errorstring)
        PRINT*, 'Error: Unable to write scalar. ptbar=', ptbarid
        PRINT*, 'Error Message: ',errorstring
      END IF

      node = 0
      DO k=kst,knd,ksk
        DO j=jst,jnd,jsk
          DO i=ist,ind,isk
            node = node + 1
            tem1(node) = pbar(i,j,k)/100.0
          END DO
        END DO
      END DO

      CALL mwritescalararrayi(gridid,pbarid,1,nxout,                    &
                              1,nyout,1,nzout,tem1,ierr)
      IF (ierr == meerr) THEN
        CALL mgeterror (errorstring)
        PRINT*, 'Error: Unable to write scalar. pbarID=', pbarid
        PRINT*, 'Error Message: ',errorstring
      END IF

      node = 0
      DO k=kst,knd,ksk
        DO j=jst,jnd,jsk
          DO i=ist,ind,isk
            node = node + 1
            tem1(node) = rhobar(i,j,k)
          END DO
        END DO
      END DO

      CALL mwritescalararrayi(gridid,rhobarid,1,nxout,                  &
                            1,nyout,1,nzout,tem1,ierr)
      IF (ierr == meerr) THEN
        CALL mgeterror (errorstring)
        PRINT*, 'Error: Unable to write scalar. rhobarID=',rhobarid
        PRINT*, 'Error Message: ',errorstring
      END IF

      IF (mstout == 1) THEN

        node = 0
        DO k=kst,knd,ksk
          DO j=jst,jnd,jsk
            DO i=ist,ind,isk
              node = node + 1
              tem1(node) = qvbar(i,j,k)*1000.0
            END DO
          END DO
        END DO

        CALL mwritescalararrayi(gridid,qvbarid,1,nxout,                 &
                                1,nyout,1,nzout,tem1,ierr)
        IF (ierr == meerr) THEN
          CALL mgeterror (errorstring)
          PRINT*, 'Error: Unable to write scalar. qvbarID=', qvbarid
          PRINT*, 'Error Message: ',errorstring
        END IF

      END IF

    END IF

  END IF

!
!-----------------------------------------------------------------------
!
!  Write the scalar arrays to the Savi3D file, and mark the end
!  of the frame to the file.
!
!-----------------------------------------------------------------------
!

  IF (varout == 1) THEN

    node = 0
    DO k=kst,knd,ksk
      DO j=jst,jnd,jsk
        DO i=ist,ind,isk
          node = node + 1
          tem1(node) = (u(i,j,k)+u(i+ishf,j,k))*0.5
        END DO
      END DO
    END DO

    CALL mwritescalararrayi (gridid, uid,1,nxout,                       &
                             1,nyout,1,nzout,tem1,ierr)
    IF (ierr == meerr) THEN
      CALL mgeterror (errorstring)
      PRINT*, 'Error: Unable to write scalar. uID=', uid
      PRINT*, 'Error Message: ',errorstring
    END IF

    node = 0
    DO k=kst,knd,ksk
      DO j=jst,jnd,jsk
        DO i=ist,ind,isk
          node = node + 1
          tem1(node) =(v(i,j,k)+v(i,j+jshf,k))*0.5
        END DO
      END DO
    END DO

    CALL mwritescalararrayi (gridid, vid,1,nxout,                       &
                             1,nyout,1,nzout,tem1,ierr)
    IF (ierr == meerr) THEN
      CALL mgeterror (errorstring)
      PRINT*, 'Error: Unable to write scalar. vID=', vid
      PRINT*, 'Error Message: ',errorstring
    END IF

    node = 0
    DO k=kst,knd,ksk
      DO j=jst,jnd,jsk
        DO i=ist,ind,isk
          node = node + 1
          tem1(node) =(w(i,j,k)+w(i,j,k+kshf))*0.5
        END DO
      END DO
    END DO

    CALL mwritescalararrayi (gridid, wid,1,nxout,                       &
          1,nyout,1,nzout,tem1,ierr)
    IF (ierr == meerr) THEN
      CALL mgeterror (errorstring)
      PRINT*, 'Error: Unable to write scalar. wID=', wid
      PRINT*, 'Error Message: ',errorstring
    END IF

    node = 0
    DO k=kst,knd,ksk
      DO j=jst,jnd,jsk
        DO i=ist,ind,isk
          node = node + 1
          tem1(node) =(u(i,j,k)+u(i+ishf,j,k)                           &
              -ubar(i,j,k)-ubar(i+ishf,j,k))*.5
        END DO
      END DO
    END DO

    CALL mwritescalararrayi (gridid, uprtid,1,nxout,                    &
                             1,nyout,1,nzout,tem1,ierr)
    IF (ierr == meerr) THEN
      CALL mgeterror (errorstring)
      PRINT*, 'Error: Unable to write scalar. uprtID=', uprtid
      PRINT*, 'Error Message: ',errorstring
    END IF

    node = 0
    DO k=kst,knd,ksk
      DO j=jst,jnd,jsk
        DO i=ist,ind,isk
          node = node + 1
          tem1(node) =(v(i,j,k)+v(i,j+jshf,k)                           &
              -vbar(i,j,k)-vbar(i,j+jshf,k))*.5
        END DO
      END DO
    END DO

    CALL mwritescalararrayi (gridid, vprtid,1,nxout,                    &
                             1,nyout,1,nzout,tem1,ierr)
    IF (ierr == meerr) THEN
      CALL mgeterror (errorstring)
      PRINT*, 'Error: Unable to write scalar. vprtID=', vprtid
      PRINT*, 'Error Message: ',errorstring
    END IF

    node = 0
    DO k=kst,knd,ksk
      DO j=jst,jnd,jsk
        DO i=ist,ind,isk
          node = node + 1
          tem1(node) =(w(i,j,k)+w(i,j,k+kshf))*0.5
        END DO
      END DO
    END DO

    CALL mwritescalararrayi (gridid, wprtid,1,nxout,                    &
                             1,nyout,1,nzout,tem1,ierr)
    IF (ierr == meerr) THEN
      CALL mgeterror (errorstring)
      PRINT*, 'Error: Unable to write scalar. wprID=', wprtid
      PRINT*, 'Error Message: ',errorstring
    END IF

    node = 0
    DO k=kst,knd,ksk
      DO j=jst,jnd,jsk
        DO i=ist,ind,isk
          node = node + 1
          tem1(node) = ptprt(i,j,k)
        END DO
      END DO
    END DO

    CALL mwritescalararrayi (gridid,ptprtid,1,nxout,                    &
                             1,nyout,1,nzout,tem1,ierr)
    IF (ierr == meerr) THEN
      CALL mgeterror (errorstring)
      PRINT*, 'Error: Unable to write scalar. ptprtID=', ptprtid
      PRINT*, 'Error Message: ',errorstring
    END IF

    node = 0
    DO k=kst,knd,ksk
      DO j=jst,jnd,jsk
        DO i=ist,ind,isk
          node = node + 1
          tem1(node) = pprt(i,j,k)*0.01
        END DO
      END DO
    END DO

    CALL mwritescalararrayi (gridid, pprtid,1,nxout,                    &
                             1,nyout,1,nzout,tem1,ierr)
    IF (ierr == meerr) THEN
      CALL mgeterror (errorstring)
      PRINT*, 'Error: Unable to write scalar. pprtID=', pprtid
      PRINT*, 'Error Message: ',errorstring
    END IF

  END IF

  IF (mstout == 1) THEN

    node = 0
    DO k=kst,knd,ksk
      DO j=jst,jnd,jsk
        DO i=ist,ind,isk
          node = node + 1
          tem1(node) =(qv(i,j,k)-qvbar(i,j,k))*1000.0
        END DO
      END DO
    END DO

    CALL mwritescalararrayi(gridid,qvprtid,1,nxout,                     &
                            1,nyout,1,nzout,tem1,ierr)
    IF (ierr == meerr) THEN
      CALL mgeterror (errorstring)
      PRINT*, 'Error: Unable to write scalar. qvprtID=', qvprtid
      PRINT*, 'Error Message: ',errorstring
    END IF

    node = 0
    DO k=kst,knd,ksk
      DO j=jst,jnd,jsk
        DO i=ist,ind,isk
          node = node + 1
          tem1(node) = qc(i,j,k)*1000.0
        END DO
      END DO
    END DO

    CALL mwritescalararrayi (gridid, qcid,1,nxout,                      &
                             1,nyout,1,nzout,tem1,ierr)
    IF (ierr == meerr) THEN
      CALL mgeterror (errorstring)
      PRINT*, 'Error: Unable to write scalar. qcID=', qcid
      PRINT*, 'Error Message: ',errorstring
    END IF

    node = 0
    DO k=kst,knd,ksk
      DO j=jst,jnd,jsk
        DO i=ist,ind,isk
          node = node + 1
          tem1(node) = qr(i,j,k)*1000.0
        END DO
      END DO
    END DO

    CALL mwritescalararrayi (gridid, qrid,1,nxout,                      &
                             1,nyout,1,nzout,tem1,ierr)
    IF (ierr == meerr) THEN
      CALL mgeterror (errorstring)
      PRINT*, 'Error: Unable to write scalar. qrID=', qrid
      PRINT*, 'Error Message: ',errorstring
    END IF

    node = 0
    DO k=kst,knd,ksk
      DO j=jst,jnd,jsk
        DO i=ist,ind,isk
          node = node + 1
          tem1(node) = (qc(i,j,k)+qr(i,j,k))*1000.0
        END DO
      END DO
    END DO

    CALL mwritescalararrayi (gridid, qwid,1,nxout,                      &
                             1,nyout,1,nzout,tem1,ierr)
    IF (ierr == meerr) THEN
      CALL mgeterror (errorstring)
      PRINT*, 'Error: Unable to write scalar. qwID=', qwid
      PRINT*, 'Error Message: ',errorstring
    END IF

    IF (iceout == 1) THEN

      node = 0
      DO k=kst,knd,ksk
        DO j=jst,jnd,jsk
          DO i=ist,ind,isk
            node = node + 1
            tem1(node) = qi(i,j,k)*1000.0
          END DO
        END DO
      END DO

      CALL mwritescalararrayi (gridid, qiid,1,nxout,                    &
                               1,nyout,1,nzout,tem1,ierr)
      IF (ierr == meerr) THEN
        CALL mgeterror (errorstring)
        PRINT*, 'Error: Unable to write scalar. qiID=', qiid
        PRINT*, 'Error Message: ',errorstring
      END IF

      node = 0
      DO k=kst,knd,ksk
        DO j=jst,jnd,jsk
          DO i=ist,ind,isk
            node = node + 1
            tem1(node) = qs(i,j,k)*1000.0
          END DO
        END DO
      END DO

      CALL mwritescalararrayi (gridid, qsid,1,nxout,                    &
                               1,nyout,1,nzout,tem1,ierr)
      IF (ierr == meerr) THEN
        CALL mgeterror (errorstring)
        PRINT*, 'Error: Unable to write scalar. qsID=', qsid
        PRINT*, 'Error Message: ',errorstring
      END IF

      node = 0
      DO k=kst,knd,ksk
        DO j=jst,jnd,jsk
          DO i=ist,ind,isk
            node = node + 1
            tem1(node) = qh(i,j,k)*1000.0
          END DO
        END DO
      END DO

      CALL mwritescalararrayi (gridid, qhid,1,nxout,                    &
                               1,nyout,1,nzout,tem1,ierr)
      IF ( ierr == meerr ) THEN
        CALL mgeterror (errorstring)
        PRINT*, 'Error: Unable to write scalar. qhID=', qhid
        PRINT*, 'Error Message: ',errorstring
      END IF

    END IF

  END IF

!
!  Vorticity:
!
  DO k=2,nz-2
    DO j=2,ny-2
      DO i=2,nx-2
        tem2(i,j,k)=                                                    &
            (v(i+1,j,k)-v(i-1,j,k)+v(i+1,j+1,k)-v(i-1,j+1,k))/          &
            (4*(x(i+1)-x(i)))-                                          &
            (u(i,j+1,k)-u(i,j-1,k)+u(i+1,j+1,k)-u(i+1,j-1,k))/          &
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

  node = 0
  DO k=kst,knd,ksk
    DO j=jst,jnd,jsk
      DO i=ist,ind,isk
        node = node + 1
        tem1(node) = tem2(i,j,k)
      END DO
    END DO
  END DO

  CALL mwritescalararrayi (gridid, vortid,1,nxout,                      &
                           1,nyout,1,nzout,tem1,ierr)
  IF ( ierr == meerr ) THEN
    CALL mgeterror (errorstring)
    PRINT*, 'Error: Unable to write scalar. VortID=', vortid
    PRINT*, 'Error Message: ',errorstring
  END IF

!
!  Divergernce:
!
  DO k=1,nz-1
    DO j=1,ny-1
      DO i=1,nx-1
        tem2(i,j,k)=                                                    &
            (u(i+1,j,k)-u(i,j,k))/(x(i+1)-x(i))+                        &
            (v(i,j+1,k)-v(i,j,k))/(y(j+1)-y(j))
      END DO
    END DO
  END DO

  node = 0
  DO k=kst,knd,ksk
!   DO 4200 i=ist,ind,isk
    DO j=jst,jnd,jsk
      DO i=ist,ind,isk
        node = node + 1
        tem1(node) = tem2(i,j,k)
      END DO
    END DO
  END DO

  CALL mwritescalararrayi (gridid, divid,1,nxout,                       &
                           1,nyout,1,nzout,tem1,ierr)
  IF ( ierr == meerr ) THEN
    CALL mgeterror (errorstring)
    PRINT*, 'Error: Unable to write scalar. DivID=', divid
    PRINT*, 'Error Message: ',errorstring
  END IF

  CALL mendcurrentframew (gridid, ierr)

  gbwrtn = 1 ! This routine has been called once.

  RETURN

  ENTRY sdmpdomn(ist1,ind1,jst1,jnd1,kst1,knd1)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  To set the start and end indicies of the model subdomain
!  in which the data is dumped out.
!
!-----------------------------------------------------------------------
!
  ist = ist1
  jst = jst1
  kst = kst1
  ind = ind1
  jnd = jnd1
  knd = knd1

  setdomn = 1

  PRINT*,'setdomn in sdmpdomn',     setdomn

  RETURN

  ENTRY sdmpskip(isk1, jsk1, ksk1)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  To set data skip parameters for data dump.
!
!-----------------------------------------------------------------------
!

  isk = isk1
  jsk = jsk1
  ksk = ksk1

  setskip = 1

  RETURN
END SUBROUTINE svidump
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE GWRISCALARARRAY            ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
!  SUBROUTINE GWriScalarArray(gridh,variable_name,input_array,nxyz)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Handles the call to Savi3D routine GrafWriteScalarArray
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  01/28/1993
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------


!   implicit none

!   integer gridh
!   character variable_name*(*)
!   integer nxyz
!   real    input_array(nxyz)
!   integer ierr

!   CALL GrafWriteScalarArray(gridh,variable_name,input_array,ierr)

!   RETURN
!
  END
