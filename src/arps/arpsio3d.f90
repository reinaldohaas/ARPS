!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE DTAREAD                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE dtaread(nx,ny,nz,nzsoil,nstyps,                              &
           hinfmt, nchanl,grdbasfn,lengbf,datafn,lendtf, time,          &
           x,y,z,zp,zpsoil, uprt ,vprt ,wprt ,ptprt, pprt ,             &
           qvprt, qscalar, tke, kmh,kmv,                                &
           ubar, vbar, wbar, ptbar, pbar, rhobar, qvbar,                &
           soiltyp,stypfrct,vegtyp,lai,roufns,veg,                      &
           tsoil,qsoil,wetcanp,snowdpth,                                &
           raing,rainc,prcrate,                                         &
           radfrc,radsw,rnflx,radswnet,radlwin,                         &
           usflx,vsflx,ptsflx,qvsflx,                                   &
           ireturn, tem1, tem2, tem3)

!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Coordinate the reading of history data of various formats.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!    2/19/1992.
!
!  MODIFICATION HISTORY:
!
!  10/06/1992  (K. Brewster)
!  Added ENTRY SETGBRD to allow the reading of
!  more than one grid-base file in a program.
!
!  4/23/93 (M. Xue)
!  New data format.
!
!  3/30/94 (M.Xue)
!  GrADS data format and surface fields.
!
!  12/01/95 (Yuhe Liu)
!  Fixed a bug in netread call for dumping grid and base data file.
!  The file name was previously given to datafn, instead of grdbasfn.
!
!  12/01/95 (Yuhe Liu)
!  Changed the order of reading grid and base file and history files.
!  Previously the code read history file first and then checked if
!  need to read the grid and base file. In order to cooperate with
!  GRIB data format which no longer contains perturbation variables,
!  the base fields have to ready so that the perturbation fields can
!  be derived from total fields substracting the base fields.
!
!  12/07/95 (Yuhe Liu)
!  Added a new data dumping format, GRIB. The parameter definitions
!  is mostly in consistance with NMC standard, despite those
!  variables which are not in the WMO and NMC's tables (Table 2)
!
!  01/22/1996 (Yuhe Liu)
!  Changed the ARPS history dumpings to dump the total values of
!  ARPS variables, instead of the perturbated values.
!
!  12/09/1998 (Donghai Wang)
!  Added the snow cover.
!
!  05/15/2002 (J. Brotzge)
!  Added to allow for multiple soil schemes.
!
!  ***NOTE: Only modified for Bin, HDF, and GRADS, as of 05/15/2002.****
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx,ny,nz The dimension of data arrays
!    nzsoil   The dimension of data arrays
!
!    hinfmt   The format of the history data dump
!             =1, machine dependent unformatted binary dump,
!             =2, formatted ascii dump,
!             =3, NCSA HDF file format.
!             =4, machine dependent packed unformatted binary dump
!
!    nchanl   FORTRAN I/O channel number for history data output.
!
!    grdbasfn Name of the grid/base state array file
!    lengbf   Length of the grid/base state data file name string
!    datafn   Name of the other time dependent data file
!    lendtf   Length of the data file name string
!
!  DATA ARRAYS READ IN:
!
!    x        x coordinate of grid points in physical/comp. space (m)
!    y        y coordinate of grid points in physical/comp. space (m)
!    z        z coordinate of grid points in computational space (m)
!    zp       z coordinate of grid points in physical space (m)
!    zpsoil   z coordinate of grid points in the soil (m)
!
!    uprt     x component of perturbation velocity (m/s)
!    vprt     y component of perturbation velocity (m/s)
!    wprt     vertical component of perturbation velocity
!             in Cartesian coordinates (m/s).
!
!    ptprt    perturbation potential temperature (K)
!    pprt     perturbation pressure (Pascal)
!    qvprt    perturbation water vapor mixing ratio (kg/kg)
!
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
!  OUTPUT:
!
!    time     The time of the input data (s)
!    x        x coordinate of grid points in physical/comp. space (m)
!    y        y coordinate of grid points in physical/comp. space (m)
!    z        z coordinate of grid points in computational space (m)
!    zp       z coordinate of grid points in physical space (m)
!    zpsoil   z coordinate of grid points in the soil (m)
!
!    uprt     x component of perturbation velocity (m/s)
!    vprt     y component of perturbation velocity (m/s)
!    wprt     vertical component of perturbation velocity in
!             Cartesian coordinates (m/s).
!
!    ptprt    perturbation potential temperature (K)
!    pprt     perturbation pressure (Pascal)
!
!    qvprt    perturbation water vapor mixing ratio (kg/kg)
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
!
!    usflx    Surface flux of u-momentum (kg/(m*s**2))
!    vsflx    Surface flux of v-momentum (kg/(m*s**2))
!    ptsflx   Surface heat flux (K*kg/(m**2 * s ))
!    qvsflx   Surface moisture flux of (kg/(m**2 * s))
!
!    ireturn  Return status indicator
!
!  WORK ARRAYS:
!
!    tem1     Temporary work array.
!    tem2     Temporary work array.
!    tem3     Temporary work array.
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'indtflg.inc'
  INCLUDE 'globcst.inc'
  INCLUDE 'grid.inc'          ! Grid & map parameters.
  INCLUDE 'mp.inc'            ! mpi parameters.
  INCLUDE 'phycst.inc'

!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  INTEGER :: nx,ny,nz          ! Number of grid points in 3 directions
  INTEGER :: nzsoil            ! Number of grid points in the soil

  REAL    :: time              ! The time of the input data (s)
  INTEGER :: hinfmt            ! The format of the history data dump
  INTEGER :: lengbf            ! Length of the grid/base state data
                               ! file name string
  INTEGER :: lendtf            ! Length of the data file name string
  CHARACTER(LEN=*) :: grdbasfn ! Name of the grid/base state array file
  CHARACTER(LEN=*) :: datafn   ! Name of the other time dependent
                               ! data file
  REAL :: x     (nx)           ! x-coord. of the physical and compu
                               ! -tational grid. Defined at u-point(m).
  REAL :: y     (ny)           ! y-coord. of the physical and compu
                               ! -tational grid. Defined at v-point(m).
  REAL :: z     (nz)           ! z-coord. of the computational grid.
                               ! Defined at w-point on the staggered
                               ! grid(m).
  REAL :: zp    (nx,ny,nz)     ! Physical height coordinate defined at
                               ! w-point of the staggered grid(m).
  REAL :: zpsoil (nx,ny,nzsoil)      ! Physical height coordinate defined at
                               ! w-point of the soil (m).

  REAL :: uprt (nx,ny,nz)      ! x component of perturbation velocity
                               ! (m/s)
  REAL :: vprt (nx,ny,nz)      ! y component of perturbation velocity
                               ! (m/s)
  REAL :: wprt (nx,ny,nz)      ! vertical component of perturbation
                               ! velocity in Cartesian coordinates
                               ! (m/s).
  REAL :: ptprt(nx,ny,nz)      ! perturbation potential temperature (K)
  REAL :: pprt (nx,ny,nz)      ! perturbation pressure (Pascal)
  REAL :: qvprt(nx,ny,nz)      ! perturbation water vapor mixing ratio
                               ! (kg/kg)
  REAL :: qscalar(nx,ny,nz,nscalar)

  REAL :: tke  (nx,ny,nz)      ! Turbulent Kinetic Energy ((m/s)**2)

  REAL :: kmh   (nx,ny,nz)     ! Horizontal turb. mixing coef. for
                               ! momentum. ( m**2/s )
  REAL :: kmv   (nx,ny,nz)     ! Vertical turb. mixing coef. for
                               ! momentum. ( m**2/s )
  REAL :: ubar  (nx,ny,nz)     ! Base state x velocity component (m/s)
  REAL :: vbar  (nx,ny,nz)     ! Base state y velocity component (m/s)
  REAL :: wbar  (nx,ny,nz)     ! Base state z velocity component (m/s)
  REAL :: ptbar (nx,ny,nz)     ! Base state potential temperature (K)
  REAL :: rhobar(nx,ny,nz)     ! Base state air density (kg/m**3)
  REAL :: pbar  (nx,ny,nz)     ! Base state pressure (Pascal)
  REAL :: qvbar (nx,ny,nz)     ! Base state water vapor mixing ratio
                               ! (kg/kg)

  INTEGER :: nstyps                    ! Number of soil type

  INTEGER :: soiltyp (nx,ny,nstyps)    ! Soil type
  REAL    :: stypfrct(nx,ny,nstyps)    ! Soil type
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
  REAL :: radswnet(nx,ny)      ! Net shortwave radiation
  REAL :: radlwin(nx,ny)       ! Incoming longwave radiation

  REAL :: usflx (nx,ny)        ! Surface flux of u-momentum (kg/(m*s**2))
  REAL :: vsflx (nx,ny)        ! Surface flux of v-momentum (kg/(m*s**2))
  REAL :: ptsflx(nx,ny)        ! Surface heat flux (K*kg/(m*s**2))
  REAL :: qvsflx(nx,ny)        ! Surface moisture flux (kg/(m**2*s))

  REAL :: tem1(nx,ny,nz)       ! Temporary work array
  REAL :: tem2(nx,ny,nz)       ! Temporary work array
  REAL :: tem3(nx,ny,nz)       ! Temporary work array

  INTEGER :: grdread,iread
  SAVE grdread

  INTEGER :: nchanl            ! FORTRAN I/O channel number for output
  INTEGER :: istat
  INTEGER :: ireturn           ! Return status indicator
  INTEGER :: grdbas            ! Wether this is a grid/base state
                               ! array dump

!-----------------------------------------------------------------------
!
! Misc. Local variables
!
!-----------------------------------------------------------------------

  REAL    :: btime                ! The time of the base state data
  REAL    :: p0inv, amin, amax
  INTEGER :: i,j,k
  INTEGER :: ierr
  LOGICAL :: fexist
  INTEGER :: packed
  CHARACTER (LEN=80) :: rname
  INTEGER :: nq, is

  INTEGER, SAVE :: itime
  INTEGER       :: istop

  DATA grdread /0/

  LOGICAL :: FOPEN            ! Whether to open file, use with binary format
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  ireturn = 0
  p0inv=1./p0

  IF (lengbf < 1) lengbf = len_trim(grdbasfn)
  IF (lendtf < 1) lendtf = len_trim(datafn)

!
!-----------------------------------------------------------------------
!
!  Open and read grid and base state data file depending on the
!  values of parameters grdin and basin, which are read in from the
!  time dependent data set. If grdin or basin is zero, the grid and
!  base state arrays have to be read in from a separate file.
!
!-----------------------------------------------------------------------
!
  IF ( hinfmt == 9 ) GO TO 500

  FOPEN = .FALSE.         ! donot open file
  IF((mp_opt >0 .AND. readsplit(FINDX_H) <= 0) .OR. myproc == 0) THEN
    FOPEN = .TRUE.
  END IF

  IF( grdread == 0 ) THEN

    grdbas = 1
    rname = runname

    IF(FOPEN) THEN

      INQUIRE(FILE=trim(grdbasfn(1:lengbf)), EXIST = fexist )
      IF( .NOT. fexist ) THEN

        INQUIRE(FILE=trim(grdbasfn(1:lengbf))//'.Z', EXIST = fexist )
        IF( fexist ) THEN
          CALL uncmprs( trim(grdbasfn(1:lengbf))//'.Z' )
        ELSE

          INQUIRE(FILE=trim(grdbasfn(1:lengbf))//'.gz', EXIST = fexist )
          IF( fexist ) THEN
            CALL uncmprs( trim(grdbasfn(1:lengbf))//'.gz' )
          ELSE
            ireturn = -1
            WRITE(6,'(/1x,a,/1x,a/)')                                   &
                'File '//trim(grdbasfn(1:lengbf))//                     &
                ' or its compressed version not found.',                &
                'Program stopped in DTAREAD.'
          END IF

        END IF
      END IF

    END IF

    !CALL mpmini(ireturn)  ! note that may not all PEs are reading because of readstride
    IF (ireturn /= 0) THEN
      CALL arpsstop('arpsstop called from dtaread during base state read',1)
    END IF

!
!-----------------------------------------------------------------------
!
!  Read grid and base state fields.
!
!-----------------------------------------------------------------------
!
    IF( hinfmt == 1 ) THEN

      IF (FOPEN) THEN

        CALL getunit( nchanl )
!
!-----------------------------------------------------------------------
!
!  Cray routines to force binary data file to be in the IEEE format
!
!-----------------------------------------------------------------------
!

        CALL asnctl ('NEWLOCAL', 1, ierr)
        CALL asnfile(trim(grdbasfn(1:lengbf)), '-F f77 -N ieee', ierr)

        OPEN(UNIT=nchanl,FILE=TRIM(grdbasfn(1:lengbf)),                 &
             STATUS='old',FORM='unformatted',ACTION='READ',IOSTAT=istat)

      END IF
      IF (readsplit(FINDX_H) > 0) CALL mpupdatei(istat, 1)

      IF( istat /= 0 ) GO TO 999

      IF (mp_opt > 0 .AND. readsplit(FINDX_H) > 0) THEN

      CALL binreadsplit(nx,ny,nz,nzsoil,nstyps, grdbas, nchanl, btime,  &
                   x,y,z,zp,zpsoil,                                     &
                   uprt ,vprt ,wprt ,ptprt, pprt,                       &
                   qvprt, qscalar, tke,kmh,kmv,                         &
                   ubar, vbar, wbar, ptbar, pbar, rhobar, qvbar,        &
                   soiltyp,stypfrct,vegtyp,lai,roufns,veg,              &
                   tsoil,qsoil,wetcanp,snowdpth,                        &
                   raing,rainc,prcrate,                                 &
                   radfrc,radsw,rnflx,radswnet,radlwin,                 &
                   usflx,vsflx,ptsflx,qvsflx,                           &
                   ireturn)

      ELSE

      CALL binread(nx,ny,nz,nzsoil,nstyps, grdbas, nchanl, btime,       &
                   x,y,z,zp,zpsoil,                                     &
                   uprt ,vprt ,wprt ,ptprt, pprt,                       &
                   qvprt, qscalar, tke,kmh,kmv,                         &
                   ubar, vbar, wbar, ptbar, pbar, rhobar, qvbar,        &
                   soiltyp,stypfrct,vegtyp,lai,roufns,veg,              &
                   tsoil,qsoil,wetcanp,snowdpth,                        &
                   raing,rainc,prcrate,                                 &
                   radfrc,radsw,rnflx,radswnet,radlwin,                 &
                   usflx,vsflx,ptsflx,qvsflx,                           &
                   ireturn)

      END IF

      IF (FOPEN) THEN
        CLOSE(UNIT=nchanl)
        CALL retunit( nchanl )
      END IF

    ELSE IF( hinfmt == 2 ) THEN

      CALL getunit( nchanl )
      OPEN(UNIT=nchanl,FILE=trim(trim(grdbasfn(1:lengbf))),             &
           STATUS='old',FORM='formatted',IOSTAT=istat)

      IF( istat /= 0 ) GO TO 999

      CALL ascread(nx,ny,nz,nzsoil,nstyps,grdbas,nchanl,                &
                   btime,x,y,z,zp,zpsoil,                               &
                   uprt ,vprt ,wprt ,ptprt, pprt,                       &
                   qvprt, qscalar, tke,kmh,kmv,                         &
                   ubar, vbar, wbar, ptbar, pbar, rhobar, qvbar,        &
                   soiltyp,stypfrct,vegtyp,lai,roufns,veg,              &
                   tsoil,qsoil,wetcanp,snowdpth,                        &
                   raing,rainc,prcrate,                                 &
                   radfrc,radsw,rnflx,radswnet,radlwin,                 &
                   usflx,vsflx,ptsflx,qvsflx,                           &
                   ireturn)

      CLOSE(UNIT=nchanl)
      CALL retunit( nchanl )

    ELSE IF( hinfmt == 3 ) THEN

      IF(mp_opt > 0 .AND. readsplit(FINDX_H) > 0) THEN

      CALL hdfreadsplit(nx,ny,nz,nzsoil,nstyps,grdbas,                  &
                   trim(grdbasfn(1:lengbf)),                            &
                   btime,x,y,z,zp,zpsoil,                               &
                   uprt ,vprt ,wprt ,ptprt, pprt,                       &
                   qvprt, qscalar, tke,kmh,kmv,                         &
                   ubar, vbar, wbar, ptbar, pbar, rhobar, qvbar,        &
                   soiltyp,stypfrct,vegtyp,lai,roufns,veg,              &
                   tsoil,qsoil,wetcanp,snowdpth,                        &
                   raing,rainc,prcrate,                                 &
                   radfrc,radsw,rnflx,radswnet,radlwin,                 &
                   usflx,vsflx,ptsflx,qvsflx,                           &
                   ireturn, tem1)
      ELSE

      CALL hdfread(nx,ny,nz,nzsoil,nstyps,grdbas,                       &
                   trim(grdbasfn(1:lengbf)),                            &
                   btime,x,y,z,zp,zpsoil,                               &
                   uprt ,vprt ,wprt ,ptprt, pprt,                       &
                   qvprt, qscalar, tke,kmh,kmv,                         &
                   ubar, vbar, wbar, ptbar, pbar, rhobar, qvbar,        &
                   soiltyp,stypfrct,vegtyp,lai,roufns,veg,              &
                   tsoil,qsoil,wetcanp,snowdpth,                        &
                   raing,rainc,prcrate,                                 &
                   radfrc,radsw,rnflx,radswnet,radlwin,                 &
                   usflx,vsflx,ptsflx,qvsflx,                           &
                   ireturn, tem1)
      END IF

    ELSE IF( hinfmt == 4 ) THEN

      CALL getunit( nchanl )
      OPEN(UNIT=nchanl,FILE=trim(trim(grdbasfn(1:lengbf))),              &
           STATUS='old',FORM='unformatted',IOSTAT=istat)

      IF( istat /= 0 ) GO TO 999

!      CALL pakread(nx,ny,nz,nzsoil,nstyps,grdbas,nchanl,                &
!                   btime,x,y,z,zp,zpsoil,                               &
!                   uprt ,vprt ,wprt ,ptprt, pprt,                       &
!                   qvprt, qc, qr, qi, qs, qh, tke,kmh,kmv,              &
!                   ubar, vbar, wbar, ptbar, pbar, rhobar, qvbar,        &
!                   soiltyp,stypfrct,vegtyp,lai,roufns,veg,              &
!                   tsoil,qsoil,wetcanp,snowdpth,                        &
!                   raing,rainc,prcrate,                                 &
!                   radfrc,radsw,rnflx,radswnet,radlwin,                 &
!                   usflx,vsflx,ptsflx,qvsflx,                           &
!                   ireturn, tem1,tem2)
      CALL pakread

      CLOSE(UNIT=nchanl)
      CALL retunit( nchanl )

    ELSE IF( hinfmt == 6 ) THEN

      CALL getunit( nchanl )
      OPEN(UNIT=nchanl,FILE=trim(trim(grdbasfn(1:lengbf))),             &
           STATUS='old',FORM='unformatted',IOSTAT=istat)

      IF( istat /= 0 ) GO TO 999

!      CALL bn2read(nx,ny,nz,nzsoil,nstyps,grdbas,nchanl,                &
!                   btime,x,y,z,zp,zpsoil,                               &
!                   uprt ,vprt ,wprt ,ptprt, pprt,                       &
!                   qvprt, qc, qr, qi, qs, qh, tke,kmh,kmv,              &
!                   ubar, vbar, wbar, ptbar, pbar, rhobar, qvbar,        &
!                   soiltyp,stypfrct,vegtyp,lai,roufns,veg,              &
!                   tsoil,qsoil,wetcanp,snowdpth,                        &
!                   raing,rainc,prcrate,                                 &
!                   radfrc,radsw,rnflx,radswnet,radlwin,                 &
!                   usflx,vsflx,ptsflx,qvsflx,                           &
!                   ireturn)

      CLOSE(UNIT=nchanl)
      CALL retunit( nchanl )

    ELSE IF (hinfmt == 7) THEN     ! NetCDF format

      packed = 0

      IF(myproc == 0 .OR. readsplit(FINDX_H) == 0)                      &
           CALL netopen(TRIM(grdbasfn(1:lengbf)),'R',nchanl)

      IF(mp_opt > 0 .AND. readsplit(FINDX_H) > 0) THEN

      CALL netreadsplit(nchanl,packed,1,grdbas,btime,                   &
                   nx,ny,nz,nzsoil,nstyps, x, y, z, zp,zpsoil,          &
                   uprt, vprt, wprt, ptprt, pprt, qvprt,                &
                   qscalar, tke,kmh,kmv,                                &
                   ubar, vbar, wbar, ptbar, pbar, rhobar, qvbar,        &
                   soiltyp,stypfrct,vegtyp,lai,roufns,veg,              &
                   tsoil,qsoil,wetcanp,snowdpth,                        &
                   raing,rainc,prcrate,                                 &
                   radfrc,radsw,rnflx,radswnet,radlwin,                 &
                   usflx,vsflx,ptsflx,qvsflx,                           &
                   tem1,ireturn)
      ELSE

      CALL netread(nchanl,packed,1,grdbas,btime,                        &
                   nx,ny,nz,nzsoil,nstyps, x, y, z, zp,zpsoil,          &
                   uprt, vprt, wprt, ptprt, pprt, qvprt,                &
                   qscalar, tke,kmh,kmv,                                &
                   ubar, vbar, wbar, ptbar, pbar, rhobar, qvbar,        &
                   soiltyp,stypfrct,vegtyp,lai,roufns,veg,              &
                   tsoil,qsoil,wetcanp,snowdpth,                        &
                   raing,rainc,prcrate,                                 &
                   radfrc,radsw,rnflx,radswnet,radlwin,                 &
                   usflx,vsflx,ptsflx,qvsflx,                           &
                   tem1,ireturn)
      END IF

      IF(myproc == 0 .OR. readsplit(FINDX_H) == 0) CALL netclose(nchanl)

    ELSE IF (hinfmt == 8) THEN     ! NetCDF packed format

      ! Do nothing, base state will be read in from data file below

    ELSE IF( hinfmt == 10 ) THEN
!
!-----------------------------------------------------------------------
!
!  Cray routines to force binary data file to be in the IEEE format
!
!-----------------------------------------------------------------------
!
      CALL asnctl ('NEWLOCAL', 1, ierr)
      CALL asnfile(trim(grdbasfn(1:lengbf)), '-F f77 -N ieee', ierr)

      CALL getunit( nchanl )

      OPEN(UNIT=nchanl,FILE=trim(trim(grdbasfn(1:lengbf))),STATUS='old',&
           FORM='unformatted',IOSTAT= istat )

      CALL gribread(nx,ny,nz,nzsoil,nstyps,grdbas,nchanl,               &
                   btime,x,y,z,zp,zpsoil,                               &
                   uprt ,vprt ,wprt ,ptprt, pprt ,                      &
                   qvprt, qscalar, tke,kmh,kmv,                         &
                   ubar, vbar, wbar, ptbar, pbar, rhobar, qvbar,        &
                   soiltyp,stypfrct,vegtyp,lai,roufns,veg,              &
                   tsoil,qsoil,wetcanp,snowdpth,                        &
                   raing,rainc,prcrate,                                 &
                   radfrc,radsw,rnflx,radswnet,radlwin,                 &
                   usflx,vsflx,ptsflx,qvsflx)

      CLOSE(UNIT=nchanl)
      CALL retunit( nchanl )

    ELSE

      WRITE(6,'(a,i3,a)')                                               &
          ' Data format flag had an invalid value ',                    &
            hinfmt ,' program stopped.'
      CALL arpsstop('arpsstop called from dtaread during base state     &
                    & read',1)

    END IF

    grdread = 1

    runname = rname

  END IF

  500   CONTINUE

!
!-----------------------------------------------------------------------
!
!  Read data fields.
!
!-----------------------------------------------------------------------
!
  grdbas = 0

  IF(FOPEN) THEN
    INQUIRE(FILE=trim(datafn(1:lendtf)), EXIST = fexist )
    IF( fexist ) GO TO 100

    INQUIRE(FILE=trim(datafn(1:lendtf))//'.Z', EXIST = fexist )
    IF( fexist ) THEN
      CALL uncmprs( trim(datafn(1:lendtf))//'.Z' )
      GO TO 100
    END IF

    INQUIRE(FILE=trim(datafn(1:lendtf))//'.gz', EXIST = fexist )
    IF( fexist ) THEN
      CALL uncmprs( trim(datafn(1:lendtf))//'.gz' )
      GO TO 100
    END IF

    WRITE(6,'(/1x,a,/1x,a/)')                                        &
       'File '//trim(datafn(1:lendtf))                               &
       //' or its compressed version not found.',                    &
       'Program stopped in DTAREAD.'
    CALL arpsstop('arpsstop called from dtaread during base read-2',1)
  END IF

  100   CONTINUE

  IF( hinfmt == 1 ) THEN

!
!-----------------------------------------------------------------------
!
!  Cray routines to force binary data file to be in the IEEE format
!
!-----------------------------------------------------------------------
!
    istat = 0
    IF (FOPEN) THEN
      CALL asnctl ('NEWLOCAL', 1, ierr)
      CALL asnfile(trim(datafn(1:lendtf)), '-F f77 -N ieee', ierr)

      CALL getunit( nchanl )
      OPEN(UNIT=nchanl,FILE=trim(trim(datafn(1:lendtf))),               &
           STATUS='old',FORM='unformatted',ACTION='READ',IOSTAT=istat)
    END IF

    IF( istat /= 0 ) GO TO 998

    IF (mp_opt > 0 .AND. readsplit(FINDX_H) > 0) THEN

    CALL binreadsplit(nx,ny,nz,nzsoil,nstyps,grdbas, nchanl,            &
                 time, x,y,z,zp,zpsoil,                                 &
                 uprt ,vprt ,wprt ,ptprt, pprt,                         &
                 qvprt, qscalar, tke,kmh,kmv,                           &
                 ubar, vbar, wbar, ptbar, pbar, rhobar, qvbar,          &
                 soiltyp,stypfrct,vegtyp,lai,roufns,veg,                &
                 tsoil,qsoil,wetcanp,snowdpth,                          &
                 raing,rainc,prcrate,                                   &
                 radfrc,radsw,rnflx,radswnet,radlwin,                   &
                 usflx,vsflx,ptsflx,qvsflx,                             &
                 ireturn)
    ELSE

    CALL binread(nx,ny,nz,nzsoil,nstyps,grdbas, nchanl,                 &
                 time, x,y,z,zp,zpsoil,                                 &
                 uprt ,vprt ,wprt ,ptprt, pprt,                         &
                 qvprt, qscalar, tke,kmh,kmv,                           &
                 ubar, vbar, wbar, ptbar, pbar, rhobar, qvbar,          &
                 soiltyp,stypfrct,vegtyp,lai,roufns,veg,                &
                 tsoil,qsoil,wetcanp,snowdpth,                          &
                 raing,rainc,prcrate,                                   &
                 radfrc,radsw,rnflx,radswnet,radlwin,                   &
                 usflx,vsflx,ptsflx,qvsflx,                             &
                 ireturn)
    END IF

    IF (FOPEN) THEN
      CLOSE(UNIT=nchanl)
      CALL retunit( nchanl )
    END IF

  ELSE IF( hinfmt == 2 ) THEN

    CALL getunit( nchanl )
    OPEN(UNIT=nchanl,FILE=trim(trim(datafn(1:lendtf))),                 &
         STATUS='old',FORM='formatted',IOSTAT=istat)

    IF( istat /= 0 ) GO TO 998

!    CALL ascread(nx,ny,nz,nzsoil,nstyps,grdbas,nchanl,time,x,y,z,zp,    &
!                 zpsoil,uprt ,vprt ,wprt ,ptprt, pprt,                  &
!                 qvprt, qscalar, tke,kmh,kmv,                           &
!                 ubar, vbar, wbar, ptbar, pbar, rhobar, qvbar,          &
!                 soiltyp,stypfrct,vegtyp,lai,roufns,veg,                &
!                 tsoil,qsoil,wetcanp,snowdpth,                          &
!                 raing,rainc,prcrate,                                   &
!                 radfrc,radsw,rnflx,radswnet,radlwin,                   &
!                 usflx,vsflx,ptsflx,qvsflx,                             &
!                 ireturn)

    CLOSE(UNIT=nchanl)
    CALL retunit( nchanl )

  ELSE IF( hinfmt == 3 ) THEN

    IF(mp_opt > 0 .AND. readsplit(FINDX_H) > 0) THEN

      CALL hdfreadsplit(nx,ny,nz,nzsoil,nstyps,grdbas,                  &
                 trim(datafn(1:lendtf)),                                &
                 time,x,y,z,zp,zpsoil,                                  &
                 uprt ,vprt ,wprt ,ptprt, pprt,                         &
                 qvprt, qscalar, tke,kmh,kmv,                           &
                 ubar, vbar, wbar, ptbar, pbar, rhobar, qvbar,          &
                 soiltyp,stypfrct,vegtyp,lai,roufns,veg,                &
                 tsoil,qsoil,wetcanp,snowdpth,                          &
                 raing,rainc,prcrate,                                   &
                 radfrc,radsw,rnflx,radswnet,radlwin,                   &
                 usflx,vsflx,ptsflx,qvsflx,                             &
                 ireturn, tem1)

    ELSE

      CALL hdfread(nx,ny,nz,nzsoil,nstyps,grdbas,                       &
                 trim(datafn(1:lendtf)),                                &
                 time,x,y,z,zp,zpsoil,                                  &
                 uprt ,vprt ,wprt ,ptprt, pprt,                         &
                 qvprt, qscalar, tke,kmh,kmv,                           &
                 ubar, vbar, wbar, ptbar, pbar, rhobar, qvbar,          &
                 soiltyp,stypfrct,vegtyp,lai,roufns,veg,                &
                 tsoil,qsoil,wetcanp,snowdpth,                          &
                 raing,rainc,prcrate,                                   &
                 radfrc,radsw,rnflx,radswnet,radlwin,                   &
                 usflx,vsflx,ptsflx,qvsflx,                             &
                 ireturn, tem1)

    END IF

  ELSE IF( hinfmt == 4 ) THEN

    CALL getunit( nchanl )
    OPEN(UNIT=nchanl,FILE=trim(trim(datafn(1:lendtf))),                 &
         STATUS='old',FORM='unformatted',IOSTAT=istat)

    IF( istat /= 0 ) GO TO 998

!    CALL pakread(nx,ny,nz,nzsoil,nstyps,grdbas,nchanl,time,x,y,z,zp,    &
!                 zpsoil,uprt ,vprt ,wprt ,ptprt, pprt,                  &
!                 qvprt, qc, qr, qi, qs, qh, tke,kmh,kmv,                &
!                 ubar, vbar, wbar, ptbar, pbar, rhobar, qvbar,          &
!                 soiltyp,stypfrct,vegtyp,lai,roufns,veg,                &
!                 tsoil,qsoil,wetcanp,snowdpth,                          &
!                 raing,rainc,prcrate,                                   &
!                 radfrc,radsw,rnflx,radswnet,radlwin,                   &
!                 usflx,vsflx,ptsflx,qvsflx,                             &
!                 ireturn, tem1,tem2)
    CALL pakread

    CLOSE(UNIT=nchanl)
    CALL retunit( nchanl )

  ELSE IF( hinfmt == 6 ) THEN

    CALL getunit( nchanl )
    OPEN(UNIT=nchanl,FILE=trim(trim(datafn(1:lendtf))),                  &
            STATUS='old',FORM='unformatted',IOSTAT=istat)

    IF( istat /= 0 ) GO TO 998

!    CALL bn2read(nx,ny,nz,nzsoil,nstyps,grdbas, nchanl, time, x,y,z,zp, &
!                 zpsoil,uprt ,vprt ,wprt ,ptprt, pprt ,                 &
!                 qvprt, qc, qr, qi, qs, qh, tke,kmh,kmv,                &
!                 ubar, vbar, wbar, ptbar, pbar, rhobar, qvbar,          &
!                 soiltyp,stypfrct,vegtyp,lai,roufns,veg,                &
!                 tsoil,qsoil,wetcanp,snowdpth,                          &
!                 raing,rainc,prcrate,                                   &
!                 radfrc,radsw,rnflx,radswnet,radlwin,                   &
!                 usflx,vsflx,ptsflx,qvsflx,                             &
!                 ireturn)

    CLOSE(UNIT=nchanl)
    CALL retunit( nchanl )

  ELSE IF (hinfmt == 7) THEN     ! NetCDF format

    packed = 0

    IF(myproc == 0 .OR. readsplit(FINDX_H) == 0)                        &
       CALL netopen(TRIM(datafn(1:lendtf)),'R',nchanl)

    IF (mp_opt >0 .AND. readsplit(FINDX_H) > 0) THEN

      CALL netreadsplit(nchanl,packed,1,grdbas,time,                    &
                   nx,ny,nz,nzsoil,nstyps, x, y, z, zp,zpsoil,          &
                   uprt, vprt, wprt, ptprt, pprt, qvprt,                &
                   qscalar, tke,kmh,kmv,                                &
                   ubar, vbar, wbar, ptbar, pbar, rhobar, qvbar,        &
                   soiltyp,stypfrct,vegtyp,lai,roufns,veg,              &
                   tsoil,qsoil,wetcanp,snowdpth,                        &
                   raing,rainc,prcrate,                                 &
                   radfrc,radsw,rnflx,radswnet,radlwin,                 &
                   usflx,vsflx,ptsflx,qvsflx,                           &
                   tem1,ireturn)
    ELSE

    CALL netread(nchanl,packed,1,grdbas,time,                           &
                 nx,ny,nz,nzsoil,nstyps, x, y, z, zp,zpsoil,            &
                 uprt, vprt, wprt, ptprt, pprt, qvprt,                  &
                 qscalar, tke,kmh,kmv,                                  &
                 ubar, vbar, wbar, ptbar, pbar, rhobar, qvbar,          &
                 soiltyp,stypfrct,vegtyp,lai,roufns,veg,                &
                 tsoil,qsoil,wetcanp,snowdpth,                          &
                 raing,rainc,prcrate,                                   &
                 radfrc,radsw,rnflx,radswnet,radlwin,                   &
                 usflx,vsflx,ptsflx,qvsflx,                             &
                 tem1, ireturn)
    END IF

    IF(myproc == 0 .OR. readsplit(FINDX_H) == 0) CALL netclose(nchanl)

  ELSE IF (hinfmt == 8) THEN     ! NetCDF in one file

    packed = 0
    itime = itime + 1

    IF( itime == 1 .AND. (myproc == 0 .OR. readsplit(FINDX_H) == 0) )   &
       CALL netopen(TRIM(datafn(1:lendtf)),'R',nchanl)

    IF (mp_opt >0 .AND. readsplit(FINDX_H) > 0) THEN

      CALL netreadsplit(nchanl,packed,itime,grdbas,time,                &
                   nx,ny,nz,nzsoil,nstyps, x, y, z, zp,zpsoil,          &
                   uprt, vprt, wprt, ptprt, pprt, qvprt,                &
                   qscalar, tke,kmh,kmv,                                &
                   ubar, vbar, wbar, ptbar, pbar, rhobar, qvbar,        &
                   soiltyp,stypfrct,vegtyp,lai,roufns,veg,              &
                   tsoil,qsoil,wetcanp,snowdpth,                        &
                   raing,rainc,prcrate,                                 &
                   radfrc,radsw,rnflx,radswnet,radlwin,                 &
                   usflx,vsflx,ptsflx,qvsflx,                           &
                   tem1,ireturn)
    ELSE

      CALL netread(nchanl,packed,itime,grdbas,time,                     &
                   nx,ny,nz,nzsoil,nstyps, x, y, z, zp,zpsoil,          &
                   uprt, vprt, wprt, ptprt, pprt, qvprt,                &
                   qscalar, tke,kmh,kmv,                                &
                   ubar, vbar, wbar, ptbar, pbar, rhobar, qvbar,        &
                   soiltyp,stypfrct,vegtyp,lai,roufns,veg,              &
                   tsoil,qsoil,wetcanp,snowdpth,                        &
                   raing,rainc,prcrate,                                 &
                   radfrc,radsw,rnflx,radswnet,radlwin,                 &
                   usflx,vsflx,ptsflx,qvsflx,                           &
                   tem1, ireturn)
    END IF

    istop = NINT((tstop-tstart)/thisdmp) + 1
    IF(itime >= istop .AND. (myproc == 0 .OR. readsplit(FINDX_H) == 0)) &
      CALL netclose(nchanl)

  ELSE IF( hinfmt == 9 ) THEN
!
!-----------------------------------------------------------------------
!
!  hinfmt = 9, read GrADS format data file. Since the DTAREAD
!  subroutine doesn't define u, v, w, pt, p, qv, we have to use
!  those base variables such ubar, vbar, wbar, ptbar, pbar, and
!  qvbar as temporary arrays to store u, v, w, pt, p, qv.
!
!-----------------------------------------------------------------------
!
    CALL gradsread(nx,ny,nz,nzsoil,nstyps,                            &
                   nchanl, trim(datafn(1:lendtf)), time,              &
                   x,y,z,zp,zpsoil,uprt,vprt,wprt,ptprt,pprt,         &
                   qvprt,qscalar, tke,kmh,kmv,                        &
                   ubar,vbar,ptbar,pbar,rhobar,qvbar,                 &
                   soiltyp,stypfrct,vegtyp,lai,roufns,veg,            &
                   tsoil,qsoil,wetcanp,snowdpth,                      &
                   raing,rainc,prcrate,                               &
                   radfrc,radsw,rnflx,radswnet,radlwin,               &
                   usflx,vsflx,ptsflx,qvsflx,                         &
                   ireturn, tem1)

    DO k = 1, nz
      DO j = 1, ny
        DO i = 1, nx
          ubar (i,j,k) = ubar (i,j,k) - uprt (i,j,k)
          vbar (i,j,k) = vbar (i,j,k) - vprt (i,j,k)
          wbar (i,j,k) = 0.0
          pbar (i,j,k) = pbar (i,j,k) - pprt (i,j,k)
          ptbar(i,j,k) = ptbar(i,j,k) - ptprt(i,j,k)
          qvbar(i,j,k) = qvbar(i,j,k) - qvprt(i,j,k)
        END DO
      END DO
    END DO

    grdin = 1

  ELSE IF( hinfmt == 10 ) THEN
!
!-----------------------------------------------------------------------
!
!  hinfmt = 10, read GRIB format data file.
!
!-----------------------------------------------------------------------
!
    CALL asnctl ('NEWLOCAL', 1, ierr)
    CALL asnfile(trim(datafn(1:lendtf)), '-F f77 -N ieee', ierr)

    CALL getunit( nchanl )

    OPEN(UNIT=nchanl,FILE=trim(trim(datafn(1:lendtf))),STATUS='old',    &
         FORM='unformatted',IOSTAT= istat )

    CALL gribread(nx,ny,nz,nzsoil,nstyps, grdbas, nchanl,time,x,y,z,zp, &
                  zpsoil,uprt ,vprt ,wprt ,ptprt, pprt ,                &
                  qvprt, qscalar, tke,kmh,kmv,                          &
                  ubar, vbar, wbar, ptbar, pbar, rhobar, qvbar,         &
                  soiltyp,stypfrct,vegtyp,lai,roufns,veg,               &
                  tsoil,qsoil,wetcanp,snowdpth,                         &
                  raing,rainc,prcrate,                                  &
                  radfrc,radsw,rnflx,radswnet,radlwin,                  &
                  usflx,vsflx,ptsflx,qvsflx)

    CLOSE(UNIT=nchanl)
    CALL retunit( nchanl )

  ELSE

    WRITE(6,'(a,i3,a)')                                                 &
          ' Data format flag had an invalid value ',                    &
           hinfmt ,' program stopped.'
    CALL arpsstop('arpsstop called from dtaread during read',1)
  END IF

  DO k=1,nz-1
    DO j=1,ny-1
      DO i=1,nx-1
        rhobar(i,j,k)= pbar(i,j,k)/                                     &
               ( rd * ptbar(i,j,k)*(pbar(i,j,k)*p0inv)**rddcp )
      END DO
    END DO
  END DO

!
!-----------------------------------------------------------------------
!
!  Set dx, dy and dz from read-in data
!
!-----------------------------------------------------------------------
!
  dx = x(2) - x(1)
  dy = y(2) - y(1)
  dz = z(2) - z(1)
!
  IF (lvldbg > 90) THEN
    IF (myproc == 0)  &
      WRITE(6,'(/1x,a/)') 'Min. and max. of the data arrays read in:'

    CALL edgfill(x,1,nx,1,nx,1,1,1,1,1,1,1,1)
    CALL a3dmax0(x,1,nx,1,nx,1,1,1,1,1,1,1,1,amax,amin)
    IF (myproc == 0)  &
       WRITE(6,'(/1x,2(a,e13.6))') 'xmin    = ', amin,',  xmax    =',amax

    CALL edgfill(y,1,ny,1,ny,1,1,1,1,1,1,1,1)
    CALL a3dmax0(y,1,ny,1,ny,1,1,1,1,1,1,1,1,amax,amin)
    IF (myproc == 0)  &
       WRITE(6,'(1x,2(a,e13.6))') 'ymin    = ', amin,',  ymax    =',amax

    CALL edgfill(z,1,nz,1,nz,1,1,1,1,1,1,1,1)
    CALL a3dmax0(z,1,nz,1,nz,1,1,1,1,1,1,1,1,amax,amin)
    IF (myproc == 0)  &
       WRITE(6,'(1x,2(a,e13.6))') 'zmin    = ', amin,',  zmax    =',amax

    CALL edgfill(zp,1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz)
    CALL a3dmax0(zp,1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz,amax,amin)
    IF (myproc == 0)  &
       WRITE(6,'(1x,2(a,e13.6))') 'zpmin   = ', amin,',  zpmax   =',amax

    CALL edgfill(ubar,1,nx,1,nx,1,ny,1,ny-1,1,nz,1,nz-1)
    CALL a3dmax0(ubar,1,nx,1,nx,1,ny,1,ny-1,1,nz,1,nz-1,amax,amin)
    IF (myproc == 0)  &
       WRITE(6,'(1x,2(a,e13.6))') 'ubarmin = ', amin,',  ubarmax =',amax

    CALL edgfill(vbar,1,nx,1,nx-1,1,ny,1,ny,1,nz,1,nz-1)
    CALL a3dmax0(vbar,1,nx,1,nx-1,1,ny,1,ny,1,nz,1,nz-1,amax,amin)
    IF (myproc == 0)  &
       WRITE(6,'(1x,2(a,e13.6))') 'vbarmin = ', amin,',  vbarmax =',amax

    CALL edgfill(wbar,1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz)
    CALL a3dmax0(wbar,1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz,amax,amin)
    IF (myproc == 0)  &
       WRITE(6,'(1x,2(a,e13.6))') 'wbarmin = ', amin,',  wbarmax =',amax

    CALL edgfill(ptbar,1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz-1)
    CALL a3dmax0(ptbar,1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz-1,amax,amin)
    IF (myproc == 0)  &
       WRITE(6,'(1x,2(a,e13.6))') 'ptbarmin= ', amin,',  ptbarmax=',amax

    CALL edgfill(pbar,1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz-1)
    CALL a3dmax0(pbar,1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz-1,amax,amin)
    IF (myproc == 0)  &
       WRITE(6,'(1x,2(a,e13.6))') 'pbarmin = ', amin,',  pbarmax =',amax

    IF( mstin == 1) THEN

      CALL edgfill(qvbar,1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz-1)
      CALL a3dmax0(qvbar,1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz-1,amax,amin)
      IF (myproc == 0)  &
         WRITE(6,'(1x,2(a,e13.6))') 'qvbarmin= ', amin, ',  qvbarmax=',amax
    END IF

    CALL edgfill(uprt,1,nx,1,nx,1,ny,1,ny-1,1,nz,1,nz-1)
    CALL a3dmax0(uprt,1,nx,1,nx,1,ny,1,ny-1,1,nz,1,nz-1,amax,amin)
    IF (myproc == 0)  &
       WRITE(6,'(1x,2(a,e13.6))') 'uprtmin = ', amin,',  uprtmax =',amax

    CALL edgfill(vprt,1,nx,1,nx-1,1,ny,1,ny,1,nz,1,nz-1)
    CALL a3dmax0(vprt,1,nx,1,nx-1,1,ny,1,ny,1,nz,1,nz-1,amax,amin)
    IF (myproc == 0)  &
       WRITE(6,'(1x,2(a,e13.6))') 'vprtmin = ', amin,',  vprtmax =',amax

    CALL edgfill(wprt,1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz)
    CALL a3dmax0(wprt,1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz,amax,amin)
    IF (myproc == 0)  &
       WRITE(6,'(1x,2(a,e13.6))') 'wprtmin = ', amin,',  wprtmax =',amax

    CALL edgfill(ptprt,1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz-1)
    CALL a3dmax0(ptprt,1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz-1,amax,amin)
    IF (myproc == 0)  &
       WRITE(6,'(1x,2(a,e13.6))') 'ptprtmin= ', amin,',  ptprtmax=',amax

    CALL edgfill(pprt,1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz-1)
    CALL a3dmax0(pprt,1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz-1,amax,amin)
    IF (myproc == 0)  &
       WRITE(6,'(1x,2(a,e13.6))') 'pprtmin = ', amin,',  pprtmax =',amax

    IF( mstin == 1) THEN

      CALL edgfill(qvprt,1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz-1)
      CALL a3dmax0lcl(qvprt,1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz-1,amax,amin)
      IF (myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                         &
                               'qvprtmin= ', amin, ',  qvprtmax= ', amax

      DO nq = 1,nscalar
        CALL edgfill   (qscalar(:,:,:,nq),1,nx,1,nx-1,1,ny,1,ny-1,        &
                        1,nz,1,nz-1)
        CALL a3dmax0lcl(qscalar(:,:,:,nq),1,nx,1,nx-1,1,ny,1,ny-1,        &
                        1,nz,1,nz-1,amax,amin)

        IF (myproc == 0)  &
           WRITE(6,'(1x,2(a,e13.6))') TRIM(qnames(nq))//'min   = ', amin, &
                               ',  '//TRIM(qnames(nq))//'max   = ', amax
      END DO

      IF(rainin == 1) THEN

        CALL edgfill(raing,1,nx,1,nx-1,1,ny,1,ny-1,1,1,1,1)
        CALL a3dmax0(raing,1,nx,1,nx-1,1,ny,1,ny-1,1,1,1,1,amax,amin)
        IF (myproc == 0)  &
           WRITE(6,'(1x,2(a,e13.6))')                                        &
            'raingmin= ', amin,',  raingmax=',amax

        CALL edgfill(rainc,1,nx,1,nx-1,1,ny,1,ny-1,1,1,1,1)
        CALL a3dmax0(rainc,1,nx,1,nx-1,1,ny,1,ny-1,1,1,1,1,amax,amin)
        IF (myproc == 0)  &
           WRITE(6,'(1x,2(a,e13.6))')                                        &
            'raincmin= ', amin,',  raincmax=',amax

      END IF

      IF( prcin == 1 ) THEN

        CALL edgfill(prcrate(1,1,1),1,nx,1,nx-1,1,ny,1,ny-1,1,1,1,1)
        CALL a3dmax0(prcrate(1,1,1),1,nx,1,nx-1,1,ny,1,ny-1,1,1,1,1,amax,amin)
        IF (myproc == 0)  &
           WRITE(6,'(1x,2(a,e13.6))')                                        &
            'prcr1min= ', amin,',  prcr1max=',amax

        CALL edgfill(prcrate(1,1,2),1,nx,1,nx-1,1,ny,1,ny-1,1,1,1,1)
        CALL a3dmax0(prcrate(1,1,2),1,nx,1,nx-1,1,ny,1,ny-1,1,1,1,1,amax,amin)
        IF (myproc == 0)  &
           WRITE(6,'(1x,2(a,e13.6))')                                        &
            'prcr2min= ', amin,',  prcr2max=',amax

        CALL edgfill(prcrate(1,1,3),1,nx,1,nx-1,1,ny,1,ny-1,1,1,1,1)
        CALL a3dmax0(prcrate(1,1,3),1,nx,1,nx-1,1,ny,1,ny-1,1,1,1,1,amax,amin)
        IF (myproc == 0)  &
           WRITE(6,'(1x,2(a,e13.6))')                                        &
            'prcr3min= ', amin,',  prcr3max=',amax

        CALL edgfill(prcrate(1,1,4),1,nx,1,nx-1,1,ny,1,ny-1,1,1,1,1)
        CALL a3dmax0(prcrate(1,1,4),1,nx,1,nx-1,1,ny,1,ny-1,1,1,1,1,amax,amin)
        IF (myproc == 0)  &
           WRITE(6,'(1x,2(a,e13.6))')                                        &
            'prcr4min= ', amin,',  prcr4max=',amax

      END IF

    END IF

    IF (tkein == 1) THEN

      CALL edgfill(tke,1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz-1)
      CALL a3dmax0(tke,1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz-1,amax,amin)
      IF (myproc == 0)  &
         WRITE(6,'(1x,2(a,e13.4))')                                          &
          'tkemin  = ', amin,',  tkemax  =',amax

    END IF

    IF (trbin == 1) THEN

      CALL edgfill(kmh,1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz-1)
      CALL a3dmax0(kmh,1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz-1,amax,amin)
      IF (myproc == 0)  &
        WRITE(6,'(1x,2(a,e13.4))')                                          &
          'kmhmin  = ', amin,',  kmhmax  =',amax

      CALL edgfill(kmv,1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz-1)
      CALL a3dmax0(kmv,1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz-1,amax,amin)
      IF (myproc == 0)  &
        WRITE(6,'(1x,2(a,e13.4))')                                          &
          'kmvmin  = ', amin,',  kmvmax  =',amax

    END IF

    IF( sfcin == 1 ) THEN

      DO is = 0, nstyp

        IF(myproc == 0)  &
           PRINT*,'Max/min of soil model variables for type index ', is

        CALL edgfill(tsoil(1,1,1,is),1,nx,1,nx-1,1,ny,1,ny-1,1,nzsoil,1,     &
              nzsoil)
        CALL a3dmax0(tsoil(1,1,1,is),1,nx,1,nx-1,1,ny,1,ny-1,1,nzsoil,1,  &
              nzsoil,amax,amin)
        IF(myproc == 0)  &
           WRITE(6,'(1x,2(a,e13.6))')                                        &
            'tsoilmin = ', amin,',  tsoilmax =',amax

        CALL edgfill(qsoil(1,1,1,is),1,nx,1,nx-1,1,ny,1,ny-1,1,nzsoil, &
               1,nzsoil)
        CALL a3dmax0(qsoil(1,1,1,is),1,nx,1,nx-1,1,ny,1,ny-1,1,nzsoil, &
               1,nzsoil,amax,amin)
        IF(myproc == 0)  &
           WRITE(6,'(1x,2(a,e13.6))')                                        &
            'qsoilmin= ', amin,',  qsoilmax=',amax

        CALL edgfill(wetcanp(1,1,is),1,nx,1,nx-1,1,ny,1,ny-1,1,1,1,1)
        CALL a3dmax0(wetcanp(1,1,is),1,nx,1,nx-1,1,ny,1,ny-1,1,1,1,1,amax,amin)
        IF(myproc == 0)  &
           WRITE(6,'(1x,2(a,e13.6))')                                        &
            'wetcmin = ', amin,',  wetcmax =',amax

      END DO

    END IF

    IF ( radin == 1 ) THEN

      CALL edgfill(radfrc,1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz-1)
      CALL a3dmax0(radfrc,1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz-1,amax,amin)
      IF (myproc == 0)  &
         WRITE(6,'(1x,2(a,e13.4))')                                          &
          'radfmin = ', amin,',  radfmax =',amax

      CALL edgfill(radsw,1,nx,1,nx-1,1,ny,1,ny-1,1,1,1,1)
      CALL a3dmax0(radsw,1,nx,1,nx-1,1,ny,1,ny-1,1,1,1,1,amax,amin)
      IF (myproc == 0)  &
         WRITE(6,'(1x,2(a,e13.4))')                                          &
          'radswmin= ', amin,',  radswmax=',amax

      CALL edgfill(rnflx,1,nx,1,nx-1,1,ny,1,ny-1,1,1,1,1)
      CALL a3dmax0(rnflx,1,nx,1,nx-1,1,ny,1,ny-1,1,1,1,1,amax,amin)
      IF (myproc == 0)  &
         WRITE(6,'(1x,2(a,e13.4))')                                          &
          'rnflxmin= ', amin,',  rnflxmax=',amax

      CALL edgfill(radswnet,1,nx,1,nx-1,1,ny,1,ny-1,1,1,1,1)
      CALL a3dmax0(radswnet,1,nx,1,nx-1,1,ny,1,ny-1,1,1,1,1,amax,amin)
      IF (myproc == 0)  &
         WRITE(6,'(1x,2(a,e13.4))')                                          &
          'radswnetmin= ', amin,',  radswnetmax=',amax

      CALL edgfill(radlwin,1,nx,1,nx-1,1,ny,1,ny-1,1,1,1,1)
      CALL a3dmax0(radlwin,1,nx,1,nx-1,1,ny,1,ny-1,1,1,1,1,amax,amin)
      IF (myproc == 0)  &
         WRITE(6,'(1x,2(a,e13.4))')                                          &
          'radlwinmin= ', amin,',  radlwinmax=',amax


    END IF

    IF ( flxin == 1 ) THEN

      CALL edgfill(usflx,1,nx,1,nx-1,1,ny,1,ny-1,1,1,1,1)
      CALL a3dmax0(usflx,1,nx,1,nx-1,1,ny,1,ny-1,1,1,1,1,amax,amin)
      IF (myproc == 0)  &
         WRITE(6,'(1x,2(a,e13.4))')                                          &
          'usflxmin= ', amin,',  usflxmax=',amax

      CALL edgfill(vsflx,1,nx,1,nx-1,1,ny,1,ny-1,1,1,1,1)
      CALL a3dmax0(vsflx,1,nx,1,nx-1,1,ny,1,ny-1,1,1,1,1,amax,amin)
      IF (myproc == 0)  &
         WRITE(6,'(1x,2(a,e13.4))')                                          &
          'vsflxmin= ', amin,',  vsflxmax=',amax

      CALL edgfill(ptsflx,1,nx,1,nx-1,1,ny,1,ny-1,1,1,1,1)
      CALL a3dmax0(ptsflx,1,nx,1,nx-1,1,ny,1,ny-1,1,1,1,1,amax,amin)
      IF (myproc == 0)  &
         WRITE(6,'(1x,2(a,e13.4))')                                          &
          'ptflxmin= ', amin,',  ptflxmax=',amax

      CALL edgfill(qvsflx,1,nx,1,nx-1,1,ny,1,ny-1,1,1,1,1)
      CALL a3dmax0(qvsflx,1,nx,1,nx-1,1,ny,1,ny-1,1,1,1,1,amax,amin)
      IF (myproc == 0)  &
         WRITE(6,'(1x,2(a,e13.4))')                                          &
          'qvflxmin= ', amin,',  qvflxmax=',amax

    END IF

  END IF

  RETURN


!##################################################################
!##################################################################
!######                                                      ######
!######                ENTRY SETGBRD                         ######
!######                                                      ######
!##################################################################
!##################################################################

  ENTRY setgbrd (iread)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Set the gridread parameter so that the grid/base state file
!  will (iread=0) or will not (iread=1) be read. This is useful
!  if more than one history file is to be accessed during a run
!  (primarily in data assimilation experiments).
!
!  As a default, the grid/base file is read on the first call to
!  dtaread and not on subsequent calls.  dtaread resets grdread
!  to 1 after each time a grid/base file is read.
!
!-----------------------------------------------------------------------
!

  grdread=iread

  RETURN

  WRITE(6,'(1x,a,a,/1x,i3,a)')                                          &
      'Error occured when opening GRIB file ',trim(datafn(1:lendtf)),   &
      'with the FILE pointer ',nchanl,' Program stopped in DTAREAD.'
  CALL arpsstop('arpsstop called from dtaread during Grid read',1)

  998   CONTINUE
  WRITE(6,'(1x,a,a,/1x,i3,a)')                                          &
      'Error occured when opening file ',trim(datafn(1:lendtf)),        &
      'using FORTRAN unit ',nchanl,' Program stopped in DTAREAD.'
  CALL arpsstop('arpsstop called from dtaread during file read',1)

  999   CONTINUE
  WRITE(6,'(1x,a,a,/1x,i3,a)')                                          &
      'Error occured when opening file ',trim(grdbasfn(1:lengbf)),      &
      'using FORTRAN unit ',nchanl,' Program stopped in DTAREAD.'

  CALL arpsstop('arpsstop called from dtaread during file read',1)

END SUBROUTINE dtaread
!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINE DTADUMP                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE dtadump(nx,ny,nz,nzsoil,nstyps,                              &
           houtfmt,nchout,filnam,grdbas,flcmprs,                        &
           u,v,w,ptprt,pprt,qv,qscalar,tke,kmh,kmv,                     &
           ubar,vbar,wbar,ptbar,pbar,rhobar,qvbar,                      &
           x,y,z,zp,zpsoil,                                             &
           soiltyp,stypfrct,vegtyp,lai,roufns,veg,                      &
           tsoil,qsoil,wetcanp,snowdpth,                                &
           raing,rainc,prcrate,                                         &
           radfrc,radsw,rnflx,radswnet,radlwin,                         &
           usflx,vsflx,ptsflx,qvsflx,                                   &
           tem1,tem2,tem3)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Coordinate the history data dump in various data formats by
!  calling the appropriate history dump routine.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  8/23/1992
!
!  MODIFICATION HISTORY:
!
!  4/4/93  (M. Xue)
!  Modified, so that data on the original staggered grid are written
!  out. Averaging to the volume center is no longer done.
!
!  9/1/94 (Y. Lu)
!  Cleaned up documentation.
!
!  4/22/1995 (M. Xue)
!  Modified external boundary condition file naming format, so
!  that the files can be directed to directory 'dirname'.
!
!  4/24/1995 (M. Xue)
!  Changed boundary file naming convention.
!
!  12/07/95 (Yuhe Liu)
!  Added a new data dumping format, GRIB. The parameter definitions
!  is mostly in consistance with NMC standard, despite those
!  variables which are not in the WMO and NMC's tables (Table 2)
!
!  01/22/1996 (Yuhe Liu)
!  Changed the ARPS history dumpings to dump the total values of
!  ARPS variables, instead of the perturbated values.
!
!  04/30/1997 (Fanyou Kong -- CMRP)
!  Add Vis5D format conversion
!
!  12/09/1998 (Donghai Wang)
!  Added the snow cover.
!
!  05/13/2002 (J. Brotzge)
!  Added additional arrays for soil schemes
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
!    nstyps   number of soil type
!
!    houtfmt  The format flag of history data dump
!    nchout   FORTRAN I/O channel number for history data output.
!    filnam   The name of history data dump file
!    grdbas   Flag indicating if this is a call for the data dump
!             of grid and base state arrays only. If so, grdbas=1.
!    flcmprs  Option for file compression (0 or 1), not used in this
!             routine. A identical variable, filcmprs, is passed
!             through include file globcst.inc
!
!    u        x component of velocity at a given time level (m/s)
!    v        y component of velocity at a given time level (m/s)
!    w        Vertical component of Cartesian velocity at a given
!             time level (m/s)
!    ptprt    Perturbation potential temperature at a given time
!             level (K)
!    pprt     Perturbation pressure at  a given time level (Pascal)
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
!    ubar     Base state x-velocity component (m/s)
!    vbar     Base state y-velocity component (m/s)
!    wbar     Base state vertical velocity component (m/s)
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
!    radswnet Net shortwave radiation, SWin-SWout
!    radlwin  Incoming longwave radiation
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
!    tem2     Temporary work array.
!    tem3     Temporary work array.
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
  INCLUDE 'grid.inc'          ! Grid & map parameters.
  INCLUDE 'bndry.inc'
  INCLUDE 'mp.inc'            ! Message passing parameters.
  INCLUDE 'exbc.inc'

!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------

  INTEGER :: nx,ny,nz          ! Number of grid points in 3 directions
  INTEGER :: nzsoil            ! Number of grid points in the soil

  INTEGER :: houtfmt           ! Flag of history data dump format
  INTEGER :: nchout            ! FORTRAN I/O channel number for output
  CHARACTER(LEN=*) :: filnam   ! The name of history dump data file
  INTEGER :: grdbas            ! If this is a grid/base state array dump
  INTEGER :: flcmprs           ! Option for file compression (0 or 1).

  REAL :: u     (nx,ny,nz)     ! Total u-velocity (m/s)
  REAL :: v     (nx,ny,nz)     ! Total v-velocity (m/s)
  REAL :: w     (nx,ny,nz)     ! Total w-velocity (m/s)
  REAL :: ptprt (nx,ny,nz)     ! Perturbation potential temperature (K)
  REAL :: pprt  (nx,ny,nz)     ! Perturbation pressure (Pascal)

  REAL :: qv    (nx,ny,nz)     ! Water vapor specific humidity (kg/kg)

  REAL :: qscalar(nx,ny,nz,nscalar)

  REAL :: tke   (nx,ny,nz)     ! Turbulent Kinetic Energy ((m/s)**2)

  REAL :: kmh   (nx,ny,nz)     ! Horizontal turb. mixing coef. for
                               ! momentum. ( m**2/s )
  REAL :: kmv   (nx,ny,nz)     ! Vertical turb. mixing coef. for
                               ! momentum. ( m**2/s )

  REAL :: ubar  (nx,ny,nz)     ! Base state x-velocity (m/s)
  REAL :: vbar  (nx,ny,nz)     ! Base state y-velocity (m/s)
  REAL :: wbar  (nx,ny,nz)     ! Base state z-velocity (m/s)
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
  REAL :: zpsoil (nx,ny,nzsoil)! The physical height coordinate defined at
                               ! w-point of the soil.

  INTEGER :: nstyps                  ! Number of soil types
  INTEGER :: soiltyp(nx,ny,nstyps)   ! Soil type
  REAL    :: stypfrct(nx,ny,nstyps)  ! Soil type fractions
  INTEGER :: vegtyp (nx,ny)          ! Vegetation type
  REAL :: lai    (nx,ny)             ! Leaf Area Index
  REAL :: roufns (nx,ny)             ! Surface roughness
  REAL :: veg    (nx,ny)             ! Vegetation fraction

  REAL :: tsoil  (nx,ny,nzsoil,0:nstyps) ! Soil temperature (K)
  REAL :: qsoil  (nx,ny,nzsoil,0:nstyps) ! Soil moisture (m**3/m**3)
  REAL :: wetcanp(nx,ny,0:nstyps)        ! Canopy water amount
  REAL :: snowdpth(nx,ny)                ! Snow depth (m)

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
  REAL :: radlwin(nx,ny)       ! Incominging longwave radiation

  REAL :: usflx (nx,ny)        ! Surface flux of u-momentum (kg/(m*s**2))
  REAL :: vsflx (nx,ny)        ! Surface flux of v-momentum (kg/(m*s**2))
  REAL :: ptsflx(nx,ny)        ! Surface heat flux (K*kg/(m*s**2))
  REAL :: qvsflx(nx,ny)        ! Surface moisture flux (kg/(m**2*s))

  REAL :: tem1  (nx,ny,nz)     ! Temporary work array
  REAL :: tem2  (nx,ny,nz)     ! Temporary work array
  REAL :: tem3  (nx,ny,nz)     ! Temporary work array

  CHARACTER (LEN=256) :: filnamr

  INTEGER :: istat, nchout0, nchout1
  INTEGER :: ierr
  INTEGER :: packed

  INTEGER, SAVE :: itime = 0
  INTEGER       :: istop

  REAL :: zpmax, rmin
!
!-----------------------------------------------------------------------
!
!  Declaration for dumping the external boundary variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: i, j, k
  INTEGER :: itema, lenstr
  INTEGER :: iyr, imon, idy, ihr, imin, isec
  CHARACTER (LEN=15 ) :: ctime
  CHARACTER (LEN=256) :: exbcfn
  CHARACTER (LEN=256) :: temchar, outdirname
  INTEGER :: istatus
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  istager = 0 ! For VIS5D and GRADS output, write out data on non-staggered
              ! grid

  IF (sfcout == 1) snowout = 1   ! snowout only independetly zero
                                 ! for pre-snow cover versions.
!
!  IF( houtfmt.eq.0 ) RETURN   ! No data dump
!
!
  IF( houtfmt == 1 ) THEN     ! Unformatted binary data dump.
!
!-----------------------------------------------------------------------
!
!  History dump of unformatted binary data
!
!-----------------------------------------------------------------------
!

    CALL getunit( nchout )

!
!-----------------------------------------------------------------------
!
!  Cray routines to force binary data file to be in the IEEE format
!
!-----------------------------------------------------------------------
!

    CALL asnctl ('NEWLOCAL', 1, ierr)
    CALL asnfile(filnam, '-F f77 -N ieee', ierr)

    IF (mp_opt > 0 .AND. joindmp(FINDX_H) > 0) THEN

      IF (myproc == 0) THEN
        OPEN(UNIT=nchout,FILE=trim(filnam),STATUS='new',               &
             FORM='unformatted',IOSTAT= istat )

        IF( istat /= 0) GO TO 999
      END IF

      CALL binjoindump(nx,ny,nz,nzsoil,nstyps,nchout, grdbas,          &
                 u,v,w,ptprt,pprt,qv,qscalar,tke,                      &
                 kmh,kmv, ubar,vbar,ptbar,pbar,rhobar,qvbar,           &
                 x,y,z,zp,zpsoil,                                      &
                 soiltyp,stypfrct,vegtyp,lai,roufns,veg,               &
                 tsoil,qsoil,wetcanp,snowdpth,                         &
                 raing,rainc,prcrate,                                  &
                 radfrc,radsw,rnflx,radswnet,radlwin,                  &
                 usflx,vsflx,ptsflx,qvsflx,                            &
                 tem1)

      IF (myproc == 0) CLOSE(UNIT=nchout)

    ELSE

      OPEN(UNIT=nchout,FILE=trim(filnam),STATUS='new',                 &
         FORM='unformatted',IOSTAT= istat )

      IF( istat /= 0) GO TO 999

      CALL bindump(nx,ny,nz,nzsoil,nstyps,nchout, grdbas,              &
                 u,v,w,ptprt,pprt,qv,qscalar,tke,                      &
                 kmh,kmv, ubar,vbar,ptbar,pbar,rhobar,qvbar,           &
                 x,y,z,zp,zpsoil,                                      &
                 soiltyp,stypfrct,vegtyp,lai,roufns,veg,               &
                 tsoil,qsoil,wetcanp,snowdpth,                         &
                 raing,rainc,prcrate,                                  &
                 radfrc,radsw,rnflx,radswnet,radlwin,                  &
                 usflx,vsflx,ptsflx,qvsflx,                            &
                 tem1)

      CLOSE(UNIT=nchout)

    END IF

    CALL retunit( nchout )

  ELSE IF( houtfmt == 2 ) THEN     ! Formatted ASCII data dump
!
!-----------------------------------------------------------------------
!
!  History dump of formatted ASCII data
!
!-----------------------------------------------------------------------
!

    CALL getunit( nchout )
    OPEN(UNIT=nchout,FILE=trim(filnam),STATUS='new',                    &
         FORM='formatted',IOSTAT= istat )

    IF( istat /= 0) GO TO 999

    CALL ascdump(nx,ny,nz,nzsoil,nstyps, nchout, grdbas,               &
                 u,v,w,ptprt,pprt,qv,qscalar,tke,                      &
                 kmh,kmv, ubar,vbar,ptbar,pbar,rhobar,qvbar,           &
                 x,y,z,zp,zpsoil,                                      &
                 soiltyp,stypfrct,vegtyp,lai,roufns,veg,               &
                 tsoil,qsoil,wetcanp,snowdpth,                         &
                 raing,rainc,prcrate,                                  &
                 radfrc,radsw,rnflx,radswnet,radlwin,                  &
                 usflx,vsflx,ptsflx,qvsflx,                            &
                 tem1)

    CLOSE(UNIT=nchout)
    CALL retunit( nchout )

  ELSE IF( houtfmt == 3 ) THEN    ! HDF4 format dump.
!
!-----------------------------------------------------------------------
!
!  History data dump in NCSA HDF4 file format.
!
!-----------------------------------------------------------------------
!

    IF (mp_opt > 0 .AND. joindmp(FINDX_H) > 0) THEN

    CALL hdfjoindump(nx,ny,nz,nzsoil,nstyps,filnam, grdbas,            &
                 u,v,w,ptprt,pprt,qv,qscalar,tke,                      &
                 kmh,kmv,ubar,vbar,ptbar,pbar,rhobar,qvbar,            &
                 x,y,z,zp,zpsoil,                                      &
                 soiltyp,stypfrct,vegtyp,lai,roufns,veg,               &
                 tsoil,qsoil,wetcanp,snowdpth,                         &
                 raing,rainc,prcrate,                                  &
                 radfrc,radsw,rnflx,radswnet,radlwin,                  &
                 usflx,vsflx,ptsflx,qvsflx,                            &
                 tem1)

    ELSE

    CALL hdfdump(nx,ny,nz,nzsoil,nstyps,filnam, grdbas,                &
                 u,v,w,ptprt,pprt,qv,qscalar,tke,                      &
                 kmh,kmv,ubar,vbar,ptbar,pbar,rhobar,qvbar,            &
                 x,y,z,zp,zpsoil,                                      &
                 soiltyp,stypfrct,vegtyp,lai,roufns,veg,               &
                 tsoil,qsoil,wetcanp,snowdpth,                         &
                 raing,rainc,prcrate,                                  &
                 radfrc,radsw,rnflx,radswnet,radlwin,                  &
                 usflx,vsflx,ptsflx,qvsflx,                            &
                 tem1)
   END IF

  ELSE IF( houtfmt == 4 ) THEN     ! Packed binary data dump.
!
!-----------------------------------------------------------------------
!
!  History dump of packed unformatted binary data
!
!-----------------------------------------------------------------------
!

    CALL getunit( nchout )
    OPEN(UNIT=nchout,FILE=trim(filnam),STATUS='new',                    &
         FORM='unformatted',IOSTAT= istat )

    IF( istat /= 0) GO TO 999

!    CALL pakdump(nx,ny,nz,nzsoil,nstyps,nchout, grdbas,                 &
!                 u,v,w,ptprt,pprt,qv,qc,qr,qi,qs,qh,tke,                &
!                 kmh,kmv, ubar,vbar,ptbar,pbar,rhobar,qvbar,            &
!                 x,y,z,zp,zpsoil,hterain, j1,j2,j3,j3soil,              &
!                 soiltyp,stypfrct,vegtyp,lai,roufns,veg,                &
!                 tsoil,qsoil,wetcanp,snowdpth,                          &
!                 raing,rainc,prcrate,                                   &
!                 radfrc,radsw,rnflx,radswnet,radlwin,                   &
!                 usflx,vsflx,ptsflx,qvsflx,                             &
!                 tem1, tem2)
    CALL pakdump

    CLOSE(UNIT=nchout)
    CALL retunit( nchout )

  ELSE IF( houtfmt == 5 ) THEN     ! Data dump for Savi3D
!
!-----------------------------------------------------------------------
!
!  History dump for Savi3D.
!  Data at all time levels are dumped into a single file.
!
!-----------------------------------------------------------------------
!
    nchout0 = 79

    CALL svidump(nx,ny,nz,nzsoil,nstyps,nchout0,filnam, grdbas,         &
                 u,v,w,ptprt,pprt,qv,qscalar,tke,                       &
                 kmh,kmv,ubar,vbar,wbar,ptbar,pbar,rhobar,qvbar,        &
                 x,y,z,zp,zpsoil,                                       &
                 soiltyp,stypfrct,vegtyp,lai,roufns,veg,                &
                 tsoil,qsoil,wetcanp,snowdpth,                          &
                 raing,rainc,prcrate,                                   &
                 radfrc,radsw,rnflx,radswnet,radlwin,                   &
                 usflx,vsflx,ptsflx,qvsflx,                             &
                 tem1,tem2)

  ELSE IF( houtfmt == 6 ) THEN     ! Unformatted binary data dump.
!
!-----------------------------------------------------------------------
!
!  History dump of unformatted binary data
!
!-----------------------------------------------------------------------
!

    CALL getunit( nchout )
    OPEN(UNIT=nchout,FILE=trim(filnam),STATUS='new',                    &
         FORM='unformatted',IOSTAT= istat )

    IF( istat /= 0) GO TO 999

    CALL bn2dump(nx,ny,nz,nzsoil,nstyps, nchout, grdbas,                &
                 u,v,w,ptprt,pprt,qv,qscalar,tke,                       &
                 kmh,kmv, ubar,vbar,ptbar,pbar,rhobar,qvbar,            &
                 x,y,z,zp,zpsoil,                                       &
                 soiltyp,stypfrct,vegtyp,lai,roufns,veg,                &
                 tsoil,qsoil,wetcanp,snowdpth,                          &
                 raing,rainc,prcrate,                                   &
                 radfrc,radsw,rnflx,radswnet,radlwin,                   &
                 usflx,vsflx,ptsflx,qvsflx,                             &
                 tem1)

    CLOSE(UNIT=nchout)
    CALL retunit( nchout )
!
!-----------------------------------------------------------------------
!
!  NetCDF format dump.
!
!-----------------------------------------------------------------------
!
  ELSE IF (houtfmt == 7) THEN    ! NetCDF format

    packed = 0

    IF ( ( nproc_x_out > nproc_x .OR. nproc_y_out > nproc_y ) .AND. splitdmp > 0 ) THEN

      CALL netsplitdump(filnam,1,packed,nx,ny,nz,nzsoil,nstyps,grdbas,  &
                   u,v,w,ptprt,pprt,qv,qscalar,tke,                     &
                   kmh,kmv,ubar,vbar,ptbar,pbar,rhobar,qvbar,           &
                   x,y,z,zp,zpsoil,raing,rainc,prcrate,                 &
                   soiltyp,stypfrct,vegtyp,lai,roufns,veg,              &
                   tsoil,qsoil,wetcanp,snowdpth,                        &
                   radfrc,radsw,rnflx,radswnet,radlwin,                 &
                   usflx,vsflx,ptsflx,qvsflx )
    ELSE

      IF(myproc == 0 .OR. joindmp(FINDX_H) == 0) CALL netopen(TRIM(filnam),'C',nchout)

      IF (mp_opt > 0 .AND. joindmp(FINDX_H) > 0) THEN

        CALL netjoindump(nchout,1,packed,nx,ny,nz,nzsoil,nstyps,grdbas, &
                     u,v,w,ptprt,pprt,qv,qscalar,tke,                   &
                     kmh,kmv,ubar,vbar,ptbar,pbar,rhobar,qvbar,         &
                     x,y,z,zp,zpsoil,raing,rainc,prcrate,               &
                     soiltyp,stypfrct,vegtyp,lai,roufns,veg,            &
                     tsoil,qsoil,wetcanp,snowdpth,                      &
                     radfrc,radsw,rnflx,radswnet,radlwin,               &
                     usflx,vsflx,ptsflx,qvsflx,                         &
                     tem1,tem2,tem3)
      ELSE

        CALL netdump(nchout,1,packed,nx,ny,nz,nzsoil,nstyps,grdbas,     &
                     u,v,w,ptprt,pprt,qv,qscalar,tke,                   &
                     kmh,kmv,ubar,vbar,ptbar,pbar,rhobar,qvbar,         &
                     x,y,z,zp,zpsoil,raing,rainc,prcrate,               &
                     soiltyp,stypfrct,vegtyp,lai,roufns,veg,            &
                     tsoil,qsoil,wetcanp,snowdpth,                      &
                     radfrc,radsw,rnflx,radswnet,radlwin,               &
                     usflx,vsflx,ptsflx,qvsflx,                         &
                     tem1,tem2,tem3)
     END IF

     IF(myproc == 0 .OR. joindmp(FINDX_H) == 0) CALL netclose(nchout)
   END IF
!
!-----------------------------------------------------------------------
!
!  Packed NetCDF format dump.
!
!-----------------------------------------------------------------------
!
  ELSE IF (houtfmt == 8) THEN    ! Write in one NetCDF file

    packed = 0
    IF (grdbas == 1) THEN
      basout = 1
      RETURN      ! Do not write grid & base file, set basout = 1 instead.
    END IF

    itime = itime + 1;
    istop = NINT((tstop-tstart)/thisdmp) + 1

    IF ( itime == 1 .AND. (myproc == 0 .OR. joindmp(FINDX_H) == 0) ) THEN
       CALL netopen(TRIM(filnam),'C',nchout)
    END IF

    IF (mp_opt > 0 .AND. joindmp(FINDX_H) > 0) THEN

    CALL netjoindump(nchout,itime,packed,nx,ny,nz,nzsoil,nstyps,grdbas, &
                 u,v,w,ptprt,pprt,qv,qscalar,tke,                       &
                 kmh,kmv,ubar,vbar,ptbar,pbar,rhobar,qvbar,             &
                 x,y,z,zp,zpsoil,raing,rainc,prcrate,                   &
                 soiltyp,stypfrct,vegtyp,lai,roufns,veg,                &
                 tsoil,qsoil,wetcanp,snowdpth,                          &
                 radfrc,radsw,rnflx,radswnet,radlwin,                   &
                 usflx,vsflx,ptsflx,qvsflx,                             &
                 tem1,tem2,tem3)
    ELSE

    CALL netdump(nchout,itime,packed,nx,ny,nz,nzsoil,nstyps,grdbas,     &
                 u,v,w,ptprt,pprt,qv,qscalar,tke,                       &
                 kmh,kmv,ubar,vbar,ptbar,pbar,rhobar,qvbar,             &
                 x,y,z,zp,zpsoil,raing,rainc,prcrate,                   &
                 soiltyp,stypfrct,vegtyp,lai,roufns,veg,                &
                 tsoil,qsoil,wetcanp,snowdpth,                          &
                 radfrc,radsw,rnflx,radswnet,radlwin,                   &
                 usflx,vsflx,ptsflx,qvsflx,                             &
                 tem1,tem2,tem3)
   END IF

   IF ( itime >= istop .AND. (myproc == 0 .OR. joindmp(FINDX_H) == 0) ) THEN
     CALL netclose(nchout)
   END IF

  ELSE IF( houtfmt == 9 ) THEN     ! Unformatted binary data dump.

    IF (mp_opt > 0 .AND. joindmp(FINDX_H) > 0) THEN

    CALL gradsjoindump(nx,ny,nz,nzsoil,nstyps,nchout, filnam, istager, &
                   u,v,w,ptprt,pprt,qv,qscalar,tke,                    &
                   kmh,kmv,ubar,vbar,wbar,ptbar,pbar,rhobar,qvbar,     &
                   x,y,z,zp,zpsoil,                                    &
                   soiltyp,stypfrct,vegtyp,lai,roufns,veg,             &
                   tsoil,qsoil,wetcanp,snowdpth,                       &
                   raing,rainc,prcrate,                                &
                   radfrc,radsw,rnflx,radswnet,radlwin,                &
                   usflx,vsflx,ptsflx,qvsflx,                          &
                   tem1,tem2)

    ELSE

    CALL gradsdump(nx,ny,nz,nzsoil,nstyps,nchout, filnam, istager,     &
                   u,v,w,ptprt,pprt,qv,qscalar,tke,                    &
                   kmh,kmv,ubar,vbar,wbar,ptbar,pbar,rhobar,qvbar,     &
                   x,y,z,zp,zpsoil,                                    &
                   soiltyp,stypfrct,vegtyp,lai,roufns,veg,             &
                   tsoil,qsoil,wetcanp,snowdpth,                       &
                   raing,rainc,prcrate,                                &
                   radfrc,radsw,rnflx,radswnet,radlwin,                &
                   usflx,vsflx,ptsflx,qvsflx,                          &
                   tem1,tem2)
    END IF

  ELSE IF( houtfmt == 10 ) THEN     ! GRIB data dump.
!
!-----------------------------------------------------------------------
!
!  History dump of GRIB format data
!
!-----------------------------------------------------------------------
!
!    lenstr = 80
!    CALL strlnth( filnam, lenstr )
    lenstr = LEN_TRIM(filnam)

    CALL asnctl ('NEWLOCAL', 1, ierr)
    CALL asnfile(filnam, '-F f77 -N ieee', ierr)

    CALL getunit( nchout )

    OPEN(UNIT=nchout,FILE=trim(filnam),STATUS='new',                    &
         FORM='unformatted',IOSTAT= istat )

    IF( istat /= 0) GO TO 999

!    CALL gribdump(nx,ny,nz,nzsoil,nstyps,                               &
!                  nchout, filnam(1:lenstr),grdbas,                      &
!                  u,v,w,ptprt,pprt,qv,qc,qr,qi,qs,qh,tke,kmh,kmv,       &
!                  ubar,vbar,ptbar,pbar,rhobar,qvbar,                    &
!                  x,y,z,zp,zpsoil,                                      &
!                  soiltyp,stypfrct,vegtyp,lai,roufns,veg,               &
!                  tsoil,qsoil,wetcanp,snowdpth,                         &
!                  raing,rainc,prcrate,                                  &
!                  radfrc,radsw,rnflx,radswnet,radlwin,                  &
!                  usflx,vsflx,ptsflx,qvsflx,                            &
!                  tem1,tem2,                                            &
!                  tem3(1,1,1),tem3(1,1,2),                              &
!                  tem3(1,1,3),tem3(1,1,4))
!
! Originally, The last two actual parameters should be INTEGER arrays. Since
! no INTEGER 2D arrays available, so they were changed to be allocated
! inside the subroutine gribdump.   -- WYH

    CALL gribdump(nx,ny,nz,nzsoil,nstyps,                               &
                  nchout, filnam(1:lenstr),grdbas,                      &
                  u,v,w,ptprt,pprt,qv,qscalar,tke,kmh,kmv,              &
                  ubar,vbar,ptbar,pbar,rhobar,qvbar,                    &
                  x,y,z,zp,zpsoil,                                      &
                  soiltyp,stypfrct,vegtyp,lai,roufns,veg,               &
                  tsoil,qsoil,wetcanp,snowdpth,                         &
                  raing,rainc,prcrate,                                  &
                  radfrc,radsw,rnflx,radswnet,radlwin,                  &
                  usflx,vsflx,ptsflx,qvsflx,                            &
                  tem1,tem3(1,1,1),tem3(1,1,2))

    CLOSE(UNIT=nchout)

    CALL retunit( nchout )

  ELSE IF( houtfmt == 11 ) THEN     ! Vis5D data dump.

!    lenstr = 80
!    CALL strlnth( filnam, lenstr )
    lenstr = LEN_TRIM(filnam)

    CALL v5ddump(nx,ny,nz,nzsoil,nstyps,filnam(1:lenstr), istager,      &
                 u,v,w,ptprt,pprt,qv,qscalar,tke,                       &
                 kmh,kmv,ubar,vbar,wbar,ptbar,pbar,rhobar,qvbar,        &
                 x,y,z,zp,zpsoil,                                       &
                 soiltyp,stypfrct,vegtyp,lai,roufns,veg,                &
                 tsoil,qsoil,wetcanp,snowdpth,                          &
                 raing,rainc,prcrate,                                   &
                 radfrc,radsw,rnflx,radswnet,radlwin,                   &
                 usflx,vsflx,ptsflx,qvsflx,                             &
                 tem1,tem2)


  END IF   ! houtfmt switches
!
!-----------------------------------------------------------------------
!
!  Compress file except for Savi3D dump.
!
!-----------------------------------------------------------------------
!

  IF( houtfmt /= 0 .AND. houtfmt /= 5  .AND.                            &
      houtfmt /= 9 .AND. houtfmt /= 11 .AND. filcmprs == 1 ) THEN

     IF (joindmp(FINDX_H) == 0 .OR. myproc == 0) CALL cmprs( filnam )
  END IF

!
!-----------------------------------------------------------------------
!
!  Create ready file, indicating history dump writing is complete
!
!-----------------------------------------------------------------------
!
  IF( readyfl == 1 .AND. houtfmt /= 0 .AND.                             &
      (myproc ==0 .OR. joindmp(FINDX_H) == 0) ) THEN
    WRITE (filnamr,'(a)') trim(filnam) // "_ready"
    CALL getunit( nchout1 )
    OPEN (UNIT=nchout1,FILE=trim(filnamr))
    WRITE (nchout1,'(a)') trim(filnam)
    CLOSE (nchout1)
    CALL retunit ( nchout1 )
  END IF
!
!-----------------------------------------------------------------------
!
!  Write u,v,w,pt,p, and qv to an additional file for
!  external boundary conditions.
!
!-----------------------------------------------------------------------
!
  IF ( exbcdmp /= 0 .AND. grdbas /= 1 ) THEN
!
!-----------------------------------------------------------------------
!
!  Write the ARPS predicted fields into external boundary format
!  files. The purpose for doing so is to produce external boundary
!  conditions from ARPS. They can be used to test the external
!  BC code, as well as to run the experiments with the ARPS generated
!  external external boundary conditions.
!
!-----------------------------------------------------------------------
!
    CALL ctim2abss( year,month,day,hour,minute,second, itema )

    itema = itema + INT(curtim)

    CALL abss2ctim( itema,iyr,imon,idy,ihr,imin,isec )

    WRITE (ctime,'(i4.4,2i2.2,a,3i2.2)')                                &
          iyr,imon,idy,'.',ihr,imin,isec

    DO k = 1, nz-1
      DO j = 1, ny-1
        DO i = 1, nx-1
          tem1(i,j,k) = ptbar(i,j,k) + ptprt(i,j,k)
          tem2(i,j,k) = pbar (i,j,k) + pprt (i,j,k)
        END DO
      END DO
    END DO

    CALL edgfill(u,   1,nx,1,nx,   1,ny,1,ny-1, 1,nz,1,nz-1)
    CALL edgfill(v,   1,nx,1,nx-1, 1,ny,1,ny,   1,nz,1,nz-1)
    CALL edgfill(w,   1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz  )
    CALL edgfill(tem1,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)
    CALL edgfill(tem2,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)
    CALL edgfill(qv,  1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)

    exbcfn = runname(1:lfnkey)//'.'//ctime
    lenstr = lfnkey + 16

    !IF (mp_opt > 0 .AND. joindmp <= 0) THEN
    !  CALL gtsplitfn(runname(1:lfnkey)//'.'//ctime,1,1,loc_x,loc_y,1,1, &
    !                 0,0,0,lvldbg,exbcfn,istat)
    !  lenstr = LEN_TRIM(exbcfn)
    !END IF

    CALL get_output_dirname(1,dirname,curtim,1,outdirname,istatus)
    temchar = exbcfn
    exbcfn  = TRIM(outdirname)//'/'//temchar
    lenstr  = LEN_TRIM(exbcfn)

    !CALL fnversn(exbcfn,lenstr)

    IF (myproc == 0) WRITE(6,'(1x,a,a)')                                &
         'Dumping to the external boundary format file: ',              &
         exbcfn(1:lenstr)

    ! Calculate rayklow, needed for exbchdfcompr > 4 format, if not computed already
    IF ( exbchdfcompr > 4 .and. rayklow <= 0) THEN
      rayklow = nz-1
      DO k=nz-1,2,-1
        zpmax = zp(1,1,k)
        DO j=1,ny-1
          DO i=1,nx-1
            zpmax = MAX( zp(i,j,k), zpmax )
          END DO
        END DO

        ! for bit-for-bit accuracy with MP version:
        rmin = zpmax
        call mpmax0(zpmax,rmin)

        IF( zpmax < zbrdmp ) THEN
          rayklow = MAX(2, k+1)
          EXIT
        END IF
      END DO
    END IF

    IF (mp_opt > 0 .AND. joindmp(FINDX_B) > 0) THEN
    CALL writejoinexbc(nx,ny,nz,exbcfn,lenstr,ctime,                    &
                  1,1,1,1,1,1,                                          &
                  qcexout,qrexout,qiexout,qsexout,qhexout,              &
                  qgexout,nqexout,zqexout,u,v,w,tem1,tem2,qv,qscalar)
    ELSE
    CALL writexbc(nx,ny,nz,exbcfn,lenstr,ctime,                         &
                  1,1,1,1,1,1,                                          &
                  qcexout,qrexout,qiexout,qsexout,qhexout,              &
                  qgexout,nqexout,zqexout,u,v,w,tem1,tem2,qv,qscalar)
    END IF

  END IF

  RETURN

  999   CONTINUE
  WRITE(6,'(1x,a,a,a,/1x,i3)')                                          &
      'Error occured when opening file ',filnam,                        &
      'using FORTRAN unit ',nchout
  CALL arpsstop(' Program stopped in DTADUMP.',1)

END SUBROUTINE dtadump
!
!##################################################################
!##################################################################
!######                                                      ######
!######            SUBROUTINE GET_DIMS_FROM_DATA             ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE get_dims_from_data(hinfmt,hisfile, nx,ny,nz,nzsoil,nstyps, &
                              ireturn)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Read in grid dimensions from base state/grid history data.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  7/17/2000
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    hinfmt   The format of the history data dump
!    grdbasfn Name of the grid/base state array file
!
!  OUTPUT:
!
!    nx,ny,nz The dimension of data arrays
!    nzsoil   Number of soil levels
!    nstyps   The number of soil types
!    ireturn  Return status indicator
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: hinfmt            ! The format of the history data dump
  CHARACTER (LEN=*  ) :: hisfile  ! Name of the grid/base state array file
  INTEGER :: nx,ny,nz          ! Number of grid points in 3 directions
  INTEGER :: nzsoil            ! Number of soil levels
  INTEGER :: nstyps            ! Number of soil types
  INTEGER :: ireturn,nchanl,ierr,istat,packed
  LOGICAL :: fexist
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  ireturn = 0
!
!-----------------------------------------------------------------------
!
!  Open and read grid and base state data file depending on the
!  values of parameters grdin and basin, which are read in from the
!  time dependent data set. If grdin or basin is zero, the grid and
!  base state arrays have to be read in from a separate file.
!
!-----------------------------------------------------------------------
!
  IF ( hinfmt == 9 ) THEN
    WRITE(6,'(/2(1x,a/))')                                         &
    'Read in from a GrADS file is currently not supported. ',      &
    'Program stopped in get_dims_from_data.'
    CALL arpsstop('arpsstop called from get_dims_from_data',1)
  END IF

  INQUIRE(FILE=trim(hisfile), EXIST = fexist )
  IF( fexist ) GO TO 200

  INQUIRE(FILE=trim(hisfile)//'.Z', EXIST = fexist )
  IF( fexist ) THEN
    CALL uncmprs( trim(hisfile)//'.Z' )
    GO TO 200
  END IF

  INQUIRE(FILE=trim(hisfile)//'.gz', EXIST = fexist )
  IF( fexist ) THEN
    CALL uncmprs( trim(hisfile)//'.gz' )
    GO TO 200
  END IF

  WRITE(6,'(/2(1x,a/))')                                              &
    'File '//trim(hisfile)//                                          &
    ' or its compressed version not found.',                          &
    'Program stopped in get_dims_from_data.'
  CALL arpsstop('arpsstop called from get_dims_from_data',1)

  200     CONTINUE

!
!-----------------------------------------------------------------------
!
!  Read grid and base state fields.
!
!-----------------------------------------------------------------------
!
  IF( hinfmt == 1 .OR. hinfmt == 6 ) THEN

    CALL getunit( nchanl )
!
!-----------------------------------------------------------------------
!
!  Cray routines to force binary data file to be in the IEEE format
!
!-----------------------------------------------------------------------
!
    CALL asnctl ('NEWLOCAL', 1, ierr)
    CALL asnfile(trim(hisfile), '-F f77 -N ieee', ierr)

    OPEN(UNIT=nchanl,FILE=trim(hisfile),                           &
         STATUS='old',FORM='unformatted',IOSTAT=istat)

    IF( istat /= 0 ) GO TO 999

    CALL bin_getdims(nchanl, nx,ny,nz,nzsoil,nstyps, ireturn)

    CLOSE(UNIT=nchanl)
    CALL retunit( nchanl )

!  ELSE IF( hinfmt == 2 ) THEN
!
!    CALL getunit( nchanl )
!    OPEN(UNIT=nchanl,FILE=trim(grdbasfn),                           &
!         STATUS='old',FORM='formatted',IOSTAT=istat)
!
!    IF( istat /= 0 ) GO TO 999
!
!    CALL asc_getdims(nchanl, nx,ny,nz,nzsoil,nstyps, ireturn)
!
!    CLOSE(UNIT=nchanl)
!    CALL retunit( nchanl )

  ELSE IF( hinfmt == 3 ) THEN

    CALL hdf_getdims(trim(hisfile),nx,ny,nz,nzsoil,nstyps, ireturn)

!  ELSE IF( hinfmt == 4 ) THEN
!
!    CALL getunit( nchanl )
!    OPEN(UNIT=nchanl,FILE=trim(grdbasfn),                           &
!         STATUS='old',FORM='unformatted',IOSTAT=istat)
!
!    IF( istat /= 0 ) GO TO 999
!
!    CALL pak_getdims(nchanl, nx,ny,nz,nzsoil,nstyps, ireturn)
!
!    CLOSE(UNIT=nchanl)
!    CALL retunit( nchanl )

  ELSE IF (hinfmt == 7 .OR. hinfmt == 8) THEN ! NetCDF format

    packed = 0
    CALL netopen(TRIM(hisfile),'R',nchanl)

    CALL net_getdims(nchanl,nx,ny,nz,nzsoil,nstyps,ireturn)

    CALL netclose(nchanl)

!
!-----------------------------------------------------------------------
!
!  Cray routines to force binary data file to be in the IEEE format
!
!-----------------------------------------------------------------------
!
!  ELSE IF( hinfmt == 10 ) THEN
!    CALL asnctl ('NEWLOCAL', 1, ierr)
!    CALL asnfile(trim(grdbasfn), '-F f77 -N ieee', ierr)
!
!    CALL getunit( nchanl )
!
!    OPEN(UNIT=nchanl,FILE=trim(grdbasfn),STATUS='old',              &
!         FORM='unformatted',IOSTAT= istat )
!
!    CALL grib_getdims(nchanl, nx,ny,nz,nzsoil,nstyps,ireturn)
!
!    CLOSE(UNIT=nchanl)
!    CALL retunit( nchanl )

  ELSE

    WRITE(6,'(a,i3,a)')                                                 &
        ' Data format flag had an invalid value ',                      &
          hinfmt ,' program stopped.'
    CALL arpsstop('arpsstop called from get_dims_from_data wrong flag',1)

  END IF

  RETURN

  999   CONTINUE
  WRITE(6,'(1x,a,a,/1x,i3,a)')                                          &
      'Error occured when opening file ',trim(hisfile),                 &
      'using FORTRAN unit ',nchanl,' Program stopped in get_dims_from_data.'

  CALL arpsstop('arpsstop called from get_dims_from_data opening file',1)

END SUBROUTINE get_dims_from_data
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE BIN_GETDIMS                ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE bin_getdims(inch, nx,ny,nz,nzsoil,nstyps, ireturn)

!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!  Read in grid dimensions from history data.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  7/17/2000.
!
!  MODIFICATION HISTORY:
!
!  09/12/2006 (Y. Wang)
!  Changed to read ARPS history file instead of grid/base file. Beside
!  ARPS domains, it also set the microphysics scalar indices.
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    inch     Channel number for binary reading.
!
!  OUTPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!    nzsoil   Number of soil levels in the vertical
!
!    nstyps   Number of soil types
!
!    ireturn  Return status indicator
!             =0, successful read of all data
!             =1, error reading data
!             =2, end-of-file reached during read attempt
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx,ny,nz          ! Number of grid points in 3 directions
  INTEGER :: nzsoil            ! Number of soil levels
  INTEGER :: nstyps            ! Number of soil types
  INTEGER :: inch              ! Channel number for binary reading
  INTEGER :: ireturn           ! Return status indicator
!
! 06/28/2002 Zuwen He
!
! fmtver??: to label each data a version.
! intver??: an integer to allow faster comparison than fmtver??,
!           which are strings.
!
! Verion 5.00: significant change in soil variables since version 4.10.
!
  CHARACTER (LEN=40) :: fmtver320,fmtver400,fmtver410,fmtver500,fmtver530
  INTEGER  :: intver,intver320,intver400,intver410,intver500,intver530

  PARAMETER (fmtver320='003.20 Binary Data',intver320=320)
  PARAMETER (fmtver400='004.00 Binary Data',intver400=400)
  PARAMETER (fmtver410='004.10 Binary Data',intver410=410)
  PARAMETER (fmtver500='005.00 Binary Data',intver500=500)
  PARAMETER (fmtver530='005.30 Binary Data',intver530=530)

  CHARACTER (LEN=40) :: fmtverin

  CHARACTER (LEN=10) :: tmunit
  INTEGER :: i
  REAL    :: time

  INTEGER :: idummy,totin,mstin,icein, nscalarin
  REAL(4) :: rdummy

  INCLUDE   'globcst.inc'
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  READ(inch,ERR=110,END=120) fmtverin

  IF (fmtverin == fmtver320) THEN
    intver=intver320
  ELSE IF (fmtverin == fmtver400) THEN
    intver=intver400
  ELSE IF (fmtverin == fmtver410) THEN
    intver=intver410
  ELSE IF (fmtverin == fmtver500) THEN
    intver=intver500
  ELSE IF (fmtverin == fmtver530) THEN
    intver=intver530
  ELSE
    WRITE(6,'(/2x,a,a,a/)')                        &
        'Incoming data format, fmtverin=',fmtverin,             &
        ', not found. The Job stopped.'
    CALL arpsstop('arpstop called from bin_getdims. ',1)
  END IF

  WRITE(6,'(/2x,a,a/)') 'Incoming data format, fmtverin=',fmtverin

  READ(inch,ERR=110,END=120) runname
  WRITE(6,'(//''  THE NAME OF THIS RUN IS:  '',A//)') runname

  READ(inch,ERR=110,END=120) nocmnt

  IF( nocmnt > 0 ) THEN
    DO i=1,nocmnt
      READ(inch,ERR=110,END=120)
    END DO
  END IF

  READ(inch,ERR=110,END=120) rdummy,tmunit
  time = rdummy

  IF (intver <= intver410) THEN

    READ(inch,ERR=110,END=120) nx, ny, nz
    nzsoil = 2  ! for version prior to 410, it is a two-layer soil model

  ELSE IF (intver >= intver500) THEN

    READ(inch,ERR=110,END=120) nx, ny, nz,nzsoil

  END IF

  WRITE(6,'(1x,a,4i5,a)') 'nx,ny,nz,nzsoil read in from data are ',     &
                           nx,ny,nz,nzsoil,' respectively.'

  READ(inch,ERR=110,END=120) idummy,idummy,idummy, mstin,  icein,       &
                             idummy,idummy,idummy,idummy,  totin,       &
                             idummy,nscalarin,idummy,idummy,idummy,     &
                             idummy,idummy,idummy,idummy,idummy

  IF (intver >= intver530) THEN
    READ(inch,ERR=110,END=120) p_qc,  p_qr,  p_qi,  p_qs,  p_qg,  p_qh, &
                               p_nc,  p_nr,  p_ni,  p_ns,  p_ng,  p_nh, &
                               p_zr,  p_zi,  p_zs,  p_zg,  p_zh,idummy, &
                             idummy,idummy,idummy,idummy,idummy,idummy, &
                             idummy,idummy,idummy,idummy,  p_cc,idummy
  ELSE
    p_qc=0; p_qr=0; p_qi=0; p_qs=0; p_qg=0; p_qh=0
    p_nc=0; p_nr=0; p_ni=0; p_ns=0; p_nh=0; p_ng=0
            p_zr=0; p_zi=0; p_zs=0; p_zh=0; p_zg=0
    p_cc=0;

    nscalar = 0
    IF (mstin == 1) THEN
      p_qc=1
      p_qr=2
      nscalar = nscalar + 2
      IF (icein == 2) THEN
        p_qi=3
        p_qs=4
        p_qh=5
        nscalar = nscalar + 3
      END IF
    END IF

  END IF

  IF (totin /= 0) THEN
    READ(inch,ERR=110,END=120) ! block of 20 REALs ...
    READ(inch,ERR=110,END=120) nstyps,idummy,idummy,idummy,idummy,  &
                               idummy,idummy,idummy,idummy,idummy,  &
                               idummy,idummy,idummy,idummy,idummy,  &
                               idummy,idummy,idummy,idummy,idummy
    WRITE (6,*) " nstyps read in as: ",nstyps
  ELSE
    WRITE (6,*) " nstyps not defined in data, set to 4"
    nstyps = 4
  ENDIF

  qnames(:)= ' '; qdescp(:)= ' '
  nscalar = 0; nscalarq = 0
  IF(p_qc > 0) THEN
    nscalarq = nscalarq + 1
    qnames(P_QC) = 'qc'; qdescp(P_QC) = 'Cloud water mixing ratio (kg/kg)'
    nscalar =  nscalar + 1; nscalarq = nscalarq + 1
  END IF
  IF(p_qr > 0) THEN
    nscalarq = nscalarq + 1
    qnames(P_QR) = 'qr'; qdescp(P_QR) = 'Rain  water mixing ratio (kg/kg)'
    nscalar =  nscalar + 1; nscalarq = nscalarq + 1
  END IF
  IF(p_qi > 0) THEN
    nscalarq = nscalarq + 1
    qnames(P_QI) = 'qi'; qdescp(P_QI) = 'Cloud ice   mixing ratio (kg/kg)'
    nscalar =  nscalar + 1; nscalarq = nscalarq + 1
  END IF
  IF(p_qs > 0) THEN
    nscalarq = nscalarq + 1
    qnames(P_QS) = 'qs'; qdescp(P_QS) = 'Snow mixing ratio (kg/kg)'
    nscalar =  nscalar + 1; nscalarq = nscalarq + 1
  END IF
  IF(p_qg > 0) THEN
    nscalarq = nscalarq + 1
    print*,'teste',nscalarq,P_QG
    P_QG      = 5
    qnames(P_QG) = 'qg'; qdescp(P_QG) = 'Graupel mixing ratio (kg/kg)'
    nscalar =  nscalar + 1; nscalarq = nscalarq + 1
  END IF
  IF(p_qh > 0) THEN
    nscalarq = nscalarq + 1
    qnames(P_QH) = 'qh'; qdescp(P_QH) = 'Hail mixing ratio (kg/kg)'
    nscalar =  nscalar + 1; nscalarq = nscalarq + 1
  END IF
  IF(p_nc > 0) THEN
    qnames(P_NC) = 'nc'; qdescp(P_NC) = 'Cloud water concentrations (#/m3)'
    nscalar =  nscalar + 1
  END IF
  IF(p_nr > 0) THEN
    qnames(P_NR) = 'nr'; qdescp(P_NR) = 'Reain water concentrations (#/m3)'
    nscalar =  nscalar + 1
  END IF
  IF(p_ni > 0) THEN
    qnames(P_NI) = 'ni'; qdescp(P_NI) = 'Cloud ice concentrations (#/m3)'
    nscalar =  nscalar + 1
  END IF
  IF(p_ns > 0) THEN
    qnames(P_NS) = 'ns'; qdescp(P_NS) = 'Snow concentrations (#/m3)'
    nscalar =  nscalar + 1
  END IF
  IF(p_ng > 0) THEN
    qnames(P_NG) = 'ng'; qdescp(P_NG) = 'Graupel concentrations (#/m3)'
    nscalar =  nscalar + 1
  END IF
  IF(p_nh > 0) THEN
    qnames(P_NH) = 'nh'; qdescp(P_NH) = 'Hail concentrations (#/m3)'
    nscalar =  nscalar + 1
  END IF
  IF(p_zr > 0) THEN
    qnames(P_ZR) = 'zr'; qdescp(P_ZR) = 'Rain reflectivity (m6/m3)'
    nscalar =  nscalar + 1
  END IF
  IF(p_zi > 0) THEN
    qnames(P_ZI) = 'zi'; qdescp(P_ZI) = 'Ice reflectivity (m6/m3)'
    nscalar =  nscalar + 1
  END IF
  IF(p_zs > 0) THEN
    qnames(P_ZS) = 'zs'; qdescp(P_ZS) = 'Snow reflectivity (m6/m3)'
    nscalar =  nscalar + 1
  END IF
  IF(p_zg > 0) THEN
    qnames(P_ZG) = 'zg'; qdescp(P_ZG) = 'Graupel reflectivity (m6/m3)'
    nscalar =  nscalar + 1
  END IF
  IF(p_zh > 0) THEN
    qnames(P_ZH) = 'zh'; qdescp(P_ZH) = 'Hail reflectivity (m6/m3)'
    nscalar =  nscalar + 1
  END IF

  IF(p_cc > 0) THEN
    qnames(P_CC) = 'cc'; qdescp(P_CC) = 'Concentration (-)'
    nscalar = nscalar + 1
  END IF

  WRITE(6,'(1x,a,2(I3,a),/,1x,I3,2a)') 'The data file contains ',       &
         nscalarin,  ' microphysics variables, ',nscalar,               &
         ' of those will be read with the ARPS IO interface and ',      &
         nscalarq,' are moist mixing ratio.',                           &
         'Those variables are: '

  DO i = 1,nscalar
    WRITE(6,FMT='(5x,2a)',ADVANCE='NO') TRIM(qnames(i)),','
    IF (MOD(i,6) == 0) WRITE(6,*)
  END DO
  WRITE(6,'(/)')

  ireturn = 0

  RETURN

  110   CONTINUE
  WRITE(6,'(/a/)') ' Error reading data in BIN_GETDIMS'
  ireturn=1
  RETURN

  120   CONTINUE
  WRITE(6,'(/a/)') ' End of file reached in BIN_GETDIMS'
  ireturn=2
  RETURN

END SUBROUTINE bin_getdims
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE ASC_GETDIMS                ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

!SUBROUTINE asc_getdims(inch, nx,ny,nz,nstyps, ireturn)
SUBROUTINE asc_getdims(inch, nx,ny,nz,nzsoil,nstyps, ireturn)

!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  7/17/2000.
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    inch     Channel number for ASCII reading.
!
!  OUTPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!    nzsoil   Number of soil levels!
!    nstyps   Number of soil types
!
!    ireturn  Return status indicator
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
  INTEGER :: nx,ny,nz          ! Number of grid points in 3 directions
  INTEGER :: nzsoil            ! Number of grid levels
  INTEGER :: nstyps            ! Number of soil types
  INTEGER :: inch              ! Channel number for binary reading
  INTEGER :: ireturn           ! Return status indicator
  CHARACTER (LEN=40) :: fmtverin
  CHARACTER (LEN=10) :: tmunit
  INTEGER :: i,nocmnt
  REAL :: time
  CHARACTER (LEN=80) :: runname

  INTEGER :: idummy,totin
  REAL::     rdummy
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  READ(inch,'(1x,a40)',ERR=110,END=120) fmtverin
  READ(inch,'(1x,a80)',ERR=110,END=120) runname

  WRITE(6,'(//''  THE NAME OF THIS RUN IS:  '',A//)') runname

  READ(inch,'(1x,i4)',ERR=110,END=120) nocmnt
  IF( nocmnt > 0 ) THEN
    DO i=1,nocmnt
      READ(inch,'(1x,a80)',ERR=110,END=120)
    END DO
  END IF

  READ(inch,'(1x,e16.8,1x,a10)',ERR=110,END=120) time,tmunit

  READ(inch,'(1x, 4i12)',ERR=110,END=120) nx,ny,nz,nzsoil

  WRITE(6,'(a,4i5,a)') '  nx,ny,nz,nzsoil read in from data are ',&
                        nx,ny,nz,nzsoil, ' respectively.'

  READ(inch,'(1x,10i8)',ERR=110,END=120)                          &
                             idummy,idummy,idummy,idummy,idummy,  &
                             idummy,idummy,idummy,idummy,totin,   &   !should be modified
                             idummy,idummy,idummy,idummy,idummy,  &
                             idummy,idummy,idummy,idummy,idummy


  READ(inch,'(1x,8E16.10)',ERR=110,END=120)                       &
                             rdummy,rdummy,rdummy,rdummy,rdummy,  &
                             rdummy,rdummy,rdummy,rdummy,rdummy,  &
                             rdummy,rdummy,rdummy,rdummy,rdummy,  &
                             rdummy,rdummy,rdummy,rdummy,rdummy

  IF (totin /= 0) THEN
    READ(inch,'(1x,10i8)',ERR=110,END=120)                          &
                               nstyps,idummy,idummy,idummy,idummy,  &
                               idummy,idummy,idummy,idummy,idummy,  &
                               idummy,idummy,idummy,idummy,idummy,  &
                               idummy,idummy,idummy,idummy,idummy
    READ(inch,'(1x,8E16.10)',ERR=110,END=120) ! block of 20 REALs ...
    WRITE (6,*) "nstyps read in as: ",nstyps
  ELSE
    WRITE (6,*) "nstyps not defined in data, set to 4"
    nstyps = 4
  ENDIF

  WRITE(6,*)
  ireturn = 0

  RETURN

  110   CONTINUE
  WRITE(6,'(/a/)') ' Error reading data in ASC_GETDIMS.'
  ireturn=1
  RETURN

  120   CONTINUE
  WRITE(6,'(/a/)') ' End of file reached in ASC_GETDIMS.'
  ireturn=2
  RETURN
END SUBROUTINE asc_getdims

!##################################################################
!##################################################################
!######                                                      ######
!######             SUBROUTINE HDF_GETDIMS                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################

SUBROUTINE hdf_getdims(filename, nx,ny,nz,nzsoil,nstyps, ireturn)

!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!  Read in grid dimensions from base state/grid history data.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  7/17/2000.
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    filename Channel number for binary reading.
!
!  OUTPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!    nzsoil   Number of soil levels!
!    ireturn  Return status indicator
!             =0, successful read of all data
!             =1, error reading data
!             =2, end-of-file reached during read attempt
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER :: stat, sd_id

  CHARACTER (LEN=*) :: filename

  INTEGER :: nx,ny,nz          ! Number of grid points in 3 directions
  INTEGER :: nzsoil            ! Number of soil levels
  INTEGER :: nstyps
  INTEGER :: ireturn           ! Return status indicator
!
! 06/28/2002 Zuwen He
!
! fmtver??: to label each data a version.
! intver??: an integer to allow faster comparison than fmtver??,
!           which are strings.
!
! Verion 5.00: significant change in soil variables since version 4.10.
!
  CHARACTER (LEN=40) :: fmtver410,fmtver500,fmtver530
  INTEGER  :: intver,intver410,intver500,intver530

  PARAMETER (fmtver410='004.10 HDF4 Coded Data',intver410=410)
  PARAMETER (fmtver500='005.00 HDF4 Coded Data',intver500=500)
  PARAMETER (fmtver530='005.30 HDF4 Coded Data',intver530=530)

  CHARACTER (LEN=40) :: fmtverin

  INTEGER            :: mstin, icein, i, nscalarin
  INTEGER            :: istat

  INCLUDE 'globcst.inc'

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  CALL hdfopen(filename,1,sd_id)

  IF (sd_id < 0) THEN
    WRITE (6,*) "HDF_GETDIMS: ERROR opening ",                              &
                 trim(filename)," for reading."
    GO TO 110
  END IF

!-----------------------------------------------------------------------
!
!  Get dimensions of data in binary file and check against
!  the dimensions passed to HDFREAD
!
!-----------------------------------------------------------------------

  CALL hdfrdc(sd_id,40,"fmtver",fmtverin,istat)

  IF (fmtverin == fmtver410) THEN
    intver=intver410
  ELSE IF (fmtverin == fmtver500) THEN
    intver=intver500
  ELSE IF (fmtverin == fmtver530) THEN
    intver=intver530
  ELSE
    intver=intver500
    WRITE(6,'(/1x,a,a,a/)')                        &
        'Incoming data format, fmtverin=',fmtverin,                 &
        ', not found. The Job stopped.'
    CALL arpsstop('arpstop called from HDF_GETDIMS. ',1)
  END IF

  WRITE(6,'(/2x,a,a/)')                        &
      'Incoming data format, fmtverin=',fmtverin

  CALL hdfrdi(sd_id,"nx",nx,istat)
  CALL hdfrdi(sd_id,"ny",ny,istat)
  CALL hdfrdi(sd_id,"nz",nz,istat)

  IF (intver <= intver410) THEN
    nzsoil = 2   ! for versions earlier than 410, it is actually
                 ! a 2 level soil model.
  ELSE IF (intver >= intver500) THEN
    CALL hdfrdi(sd_id,"nzsoil",nzsoil,istat)
  END IF

  CALL hdfrdi(sd_id,"nstyp",nstyps,istat)

  WRITE(6,'(1x,a,4i5,a)') 'nx,ny,nz,nzsoil read in from data are ', &
       nx,ny,nz,nzsoil,' respectively.'

  WRITE(6,'(1x,a,I3)') "nstyps read in as: ",nstyps

  IF (intver >= intver530) THEN       ! New version

    CALL hdfrdi(sd_id, "nscalar", nscalarin, istat)
    CALL hdfrdi(sd_id, 'P_QC',    P_QC, istat)
    CALL hdfrdi(sd_id, 'P_QR',    P_QR, istat)
    CALL hdfrdi(sd_id, 'P_QI',    P_QI, istat)
    CALL hdfrdi(sd_id, 'P_QS',    P_QS, istat)
    CALL hdfrdi(sd_id, 'P_QG',    P_QG, istat)
    CALL hdfrdi(sd_id, 'P_QH',    P_QH, istat)

    CALL hdfrdi(sd_id, 'P_NC',    P_NC, istat)
    CALL hdfrdi(sd_id, 'P_NR',    P_NR, istat)
    CALL hdfrdi(sd_id, 'P_NI',    P_NI, istat)
    CALL hdfrdi(sd_id, 'P_NS',    P_NS, istat)
    CALL hdfrdi(sd_id, 'P_NG',    P_NG, istat)
    CALL hdfrdi(sd_id, 'P_NH',    P_NH, istat)

    CALL hdfrdi(sd_id, 'P_ZR',    P_ZR, istat)
    CALL hdfrdi(sd_id, 'P_ZI',    P_ZI, istat)
    CALL hdfrdi(sd_id, 'P_ZS',    P_ZS, istat)
    CALL hdfrdi(sd_id, 'P_ZG',    P_ZG, istat)
    CALL hdfrdi(sd_id, 'P_ZH',    P_ZH, istat)

    CALL hdfrdi(sd_id, 'P_CC',    P_CC, istat)

  ELSE                        ! maybe old version

    CALL hdfrdi(sd_id,"mstflg",mstin,istat)
    CALL hdfrdi(sd_id,"iceflg",icein,istat)

    nscalar = 0
    IF (mstin == 1) THEN
      P_QC = 1
      P_QR = 2
      nscalar = nscalar + 2
      IF (icein == 1) THEN
        P_QI = 3
        P_QS = 4
        P_QH = 5
        nscalar = nscalar + 3
      END IF
    END IF
    nscalarin = nscalar
  END IF

  qnames(:)= ' '; qdescp(:)= ' '
  nscalar = 0; nscalarq = 0
  IF(p_qc > 0) THEN
    IF (lvldbg > 0) WRITE(*,'(1x,a)') 'Found QC ...'
    qnames(P_QC) = 'qc'; qdescp(P_QC) = 'Cloud water mixing ratio (kg/kg)'
    nscalar = nscalar + 1;     nscalarq = nscalarq + 1
  END IF
  IF(p_qr > 0) THEN
    IF (lvldbg > 0) WRITE(*,'(1x,a)') 'Found QR ...'
    qnames(P_QR) = 'qr'; qdescp(P_QR) = 'Rain  water mixing ratio (kg/kg)'
    nscalar = nscalar + 1;     nscalarq = nscalarq + 1
  END IF
  IF(p_qi > 0) THEN
    IF (lvldbg > 0) WRITE(*,'(1x,a)') 'Found QI ...'
    qnames(P_QI) = 'qi'; qdescp(P_QI) = 'Cloud ice   mixing ratio (kg/kg)'
    nscalar = nscalar + 1;    nscalarq = nscalarq + 1
  END IF
  IF(p_qs > 0) THEN
    IF (lvldbg > 0) WRITE(*,'(1x,a)') 'Found QS ...'
    qnames(P_QS) = 'qs'; qdescp(P_QS) = 'Snow mixing ratio (kg/kg)'
    nscalar = nscalar + 1;    nscalarq = nscalarq + 1
  END IF
  IF(p_qg > 0) THEN
    IF (lvldbg > 0) WRITE(*,'(1x,a)') 'Found QG ...'
    qnames(P_QG) = 'qg'; qdescp(P_QG) = 'Graupel mixing ratio (kg/kg)'
    nscalar = nscalar + 1;    nscalarq = nscalarq + 1
  END IF
  IF(p_qh > 0) THEN
    IF (lvldbg > 0) WRITE(*,'(1x,a)') 'Found QH ...'
    qnames(P_QH) = 'qh'; qdescp(P_QH) = 'Hail mixing ratio (kg/kg)'
    nscalar = nscalar + 1;    nscalarq = nscalarq + 1
  END IF
  IF(p_nc > 0) THEN
    IF (lvldbg > 0) WRITE(*,'(1x,a)') 'Found NC ...'
    qnames(P_NC) = 'nc'; qdescp(P_NC) = 'Cloud water concentrations (#/m3)'
    nscalar = nscalar + 1
  END IF
  IF(p_nr > 0) THEN
    IF (lvldbg > 0) WRITE(*,'(1x,a)') 'Found NR ...'
    qnames(P_NR) = 'nr'; qdescp(P_NR) = 'Reain water concentrations (#/m3)'
    nscalar = nscalar + 1
  END IF
  IF(p_ni > 0) THEN
    IF (lvldbg > 0) WRITE(*,'(1x,a)') 'Found NI ...'
    qnames(P_NI) = 'ni'; qdescp(P_NI) = 'Cloud ice concentrations (#/m3)'
    nscalar = nscalar + 1
  END IF
  IF(p_ns > 0) THEN
    IF (lvldbg > 0) WRITE(*,'(1x,a)') 'Found NS ...'
    qnames(P_NS) = 'ns'; qdescp(P_NS) = 'Snow concentrations (#/m3)'
    nscalar = nscalar + 1
  END IF
  IF(p_ng > 0) THEN
    IF (lvldbg > 0) WRITE(*,'(1x,a)') 'Found NG ...'
    qnames(P_NG) = 'ng'; qdescp(P_NG) = 'Graupel concentrations (#/m3)'
    nscalar = nscalar + 1
  END IF
  IF(p_nh > 0) THEN
    IF (lvldbg > 0) WRITE(*,'(1x,a)') 'Found NH ...'
    qnames(P_NH) = 'nh'; qdescp(P_NH) = 'Hail concentrations (#/m3)'
    nscalar = nscalar + 1
  END IF
  IF(p_zr > 0) THEN
    IF (lvldbg > 0) WRITE(*,'(1x,a)') 'Found ZR ...'
    qnames(P_ZR) = 'zr'; qdescp(P_ZR) = 'Rain reflectivity (m6/m3)'
    nscalar = nscalar + 1
  END IF
  IF(p_zi > 0) THEN
    IF (lvldbg > 0) WRITE(*,'(1x,a)') 'Found ZI ...'
    qnames(P_ZI) = 'zi'; qdescp(P_ZI) = 'Ice reflectivity (m6/m3)'
    nscalar = nscalar + 1
  END IF
  IF(p_zs > 0) THEN
    IF (lvldbg > 0) WRITE(*,'(1x,a)') 'Found ZS ...'
    qnames(P_ZS) = 'zs'; qdescp(P_ZS) = 'Snow reflectivity (m6/m3)'
    nscalar = nscalar + 1
  END IF
  IF(p_zg > 0) THEN
    IF (lvldbg > 0) WRITE(*,'(1x,a)') 'Found ZG ...'
    qnames(P_ZG) = 'zg'; qdescp(P_ZG) = 'Graupel reflectivity (m6/m3)'
    nscalar = nscalar + 1
  END IF
  IF(p_zh > 0) THEN
    IF (lvldbg > 0) WRITE(*,'(1x,a)') 'Found ZH ...'
    qnames(P_ZH) = 'zh'; qdescp(P_ZH) = 'Hail reflectivity (m6/m3)'
    nscalar = nscalar + 1
  END IF

  IF(p_cc > 0) THEN
    IF (lvldbg > 0) WRITE(*,'(1x,a)') 'Found CC ...'
    qnames(P_CC) = 'cc'; qdescp(P_CC) = 'Concentration (-)'
    nscalar = nscalar + 1
  END IF

  WRITE(6,'(1x,a,I3,a,/,2(1x,I3,a,/),1x,a)')                            &
       'The data file contains ', nscalarin,' microphysics variables',  &
       nscalar,' of those will be read with the ARPS IO interface and', &
       nscalarq,' are moist mixing ratios. ',                           &
       'Those variables are: '

  DO i = 1,nscalar
    WRITE(6,FMT='(5x,2a)',ADVANCE='NO') TRIM(qnames(i)),','
    IF (MOD(i,6) == 0) WRITE(6,*)
  END DO
  WRITE(6,'(/)')

  ireturn = 0
  GO TO 130

!-----------------------------------------------------------------------
!
!  Error during read
!
!-----------------------------------------------------------------------

  110   CONTINUE
  WRITE(6,'(/a/)') ' Error reading data in HDF_GETDIMS.'
  ireturn=1

  GO TO 130

!-----------------------------------------------------------------------
!
!  End-of-file during read
!
!-----------------------------------------------------------------------

!  120   CONTINUE
!  WRITE(6,'(/a/)') ' End of file reached in HDF_GETDIMS.'
!  ireturn=2

  130   CONTINUE

!tmp  stat = sfendacc(sd_id)   ! is this necessary?
  CALL hdfclose(sd_id,stat)

  RETURN
END SUBROUTINE hdf_getdims
!
!##################################################################
!##################################################################
!######                                                      ######
!######             SUBROUTINE NET_GETDIMS                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE net_getdims(ncid, nx,ny,nz,nzsoil,nstyps, istat)

!
!-----------------------------------------------------------------------
!
!  PURPOSE:  Read in grid dimensions from history data.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yunheng Wang
!  9/15/2006.
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    filename Channel number for binary reading.
!
!  OUTPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!    nzsoil   Number of soil levels!
!    ireturn  Return status indicator
!             =0, successful read of all data
!             =1, error reading data
!             =2, end-of-file reached during read attempt
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: ncid

  INTEGER :: nx,ny,nz          ! Number of grid points in 3 directions
  INTEGER :: nzsoil            ! Number of soil levels
  INTEGER :: nstyps
  INTEGER :: istat           ! Return status indicator

  CHARACTER(LEN=17), PARAMETER :: fmtver500 = '005.10 NetCDF 3.0'
  CHARACTER(LEN=17), PARAMETER :: fmtver530 = '005.30 NetCDF 3.0'

  INTEGER                      :: intver,intver500,intver530

  CHARACTER (LEN=40) :: fmtverin

  INTEGER            :: mstin, icein, n, nscalarin

  INCLUDE 'globcst.inc'

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


!-----------------------------------------------------------------------
!
!  Get dimensions of data in netCDF file
!
!-----------------------------------------------------------------------

  CALL netreadattstr(ncid,"FMTVER",fmtverin,istat)

  IF (fmtverin(1:17) == fmtver500) THEN
    intver=intver500
  ELSE IF (fmtverin(1:17) == fmtver530) THEN
    intver=intver530
  ELSE
    WRITE(6,'(/1x,a,a,a/)')                                             &
              'Incoming data format, fmtverin = ',TRIM(fmtverin),       &
              ', not found. The Job stopped.'
    CALL arpsstop('arpstop called from net_getdims.',1)

  END IF

  WRITE(6,'(/1x,a,a/)')                        &
      'Incoming data format, fmtverin=',fmtverin

  CALL netreaddims(ncid,nx,ny,nz,nzsoil,nstyps,istat)

  WRITE(6,'(1x,a,4i5,a)') 'nx,ny,nz,nzsoil read in from data are ', &
       nx,ny,nz,nzsoil,' respectively.'

  WRITE(6,'(1x,a,I3)') "nstyps read in as: ",nstyps

  IF (intver >= intver530) THEN       ! New version

    CALL netreadatti(ncid,'nscalar',nscalarin,istat)
    CALL netreadatti(ncid, 'P_QC',    P_QC, istat)
    CALL netreadatti(ncid, 'P_QR',    P_QR, istat)
    CALL netreadatti(ncid, 'P_QI',    P_QI, istat)
    CALL netreadatti(ncid, 'P_QS',    P_QS, istat)
    CALL netreadatti(ncid, 'P_QG',    P_QG, istat)
    CALL netreadatti(ncid, 'P_QH',    P_QH, istat)

    CALL netreadatti(ncid, 'P_NC',    P_NC, istat)
    CALL netreadatti(ncid, 'P_NR',    P_NR, istat)
    CALL netreadatti(ncid, 'P_NI',    P_NI, istat)
    CALL netreadatti(ncid, 'P_NS',    P_NS, istat)
    CALL netreadatti(ncid, 'P_NG',    P_NG, istat)
    CALL netreadatti(ncid, 'P_NH',    P_NH, istat)

    CALL netreadatti(ncid, 'P_ZR',    P_ZR, istat)
    CALL netreadatti(ncid, 'P_ZI',    P_ZI, istat)
    CALL netreadatti(ncid, 'P_ZS',    P_ZS, istat)
    CALL netreadatti(ncid, 'P_ZG',    P_ZG, istat)
    CALL netreadatti(ncid, 'P_ZH',    P_ZH, istat)

    CALL netreadatti(ncid, 'P_CC',    P_CC, istat)

  ELSE                        ! maybe old version

    CALL netreadatti(ncid,"MSTFLG",mstin,istat)
    CALL netreadatti(ncid,"ICEFLG",icein,istat)

    nscalar = 0
    IF (mstin == 1) THEN
      P_QC = 1
      P_QR = 2
      nscalar = nscalar + 2
      IF (icein == 1) THEN
        P_QI = 3
        P_QS = 4
        P_QH = 5
        nscalar = nscalar + 3
      END IF
    END IF
  END IF

  qnames(:)= ' '; qdescp(:)= ' '
  nscalar = 0; nscalarq = nscalarq + 1
  IF(p_qc > 0) THEN
    qnames(P_QC) = 'qc'; qdescp(P_QC) = 'Cloud water mixing ratio (kg/kg)'
    nscalar =  nscalar + 1; nscalarq = nscalarq + 1
  END IF
  IF(p_qr > 0) THEN
    qnames(P_QR) = 'qr'; qdescp(P_QR) = 'Rain  water mixing ratio (kg/kg)'
    nscalar =  nscalar + 1; nscalarq = nscalarq + 1
  END IF
  IF(p_qi > 0) THEN
    qnames(P_QI) = 'qi'; qdescp(P_QI) = 'Cloud ice   mixing ratio (kg/kg)'
    nscalar =  nscalar + 1; nscalarq = nscalarq + 1
  END IF
  IF(p_qs > 0) THEN
    qnames(P_QS) = 'qs'; qdescp(P_QS) = 'Snow mixing ratio (kg/kg)'
    nscalar =  nscalar + 1; nscalarq = nscalarq + 1
  END IF
  IF(p_qg > 0) THEN
    qnames(P_QG) = 'qg'; qdescp(P_QG) = 'Graupel mixing ratio (kg/kg)'
    nscalar =  nscalar + 1; nscalarq = nscalarq + 1
  END IF
  IF(p_qh > 0) THEN
    qnames(P_QH) = 'qh'; qdescp(P_QH) = 'Hail mixing ratio (kg/kg)'
    nscalar =  nscalar + 1; nscalarq = nscalarq + 1
  END IF
  IF(p_nc > 0) THEN
    qnames(P_NC) = 'nc'; qdescp(P_NC) = 'Cloud water concentrations (#/m3)'
    nscalar =  nscalar + 1
  END IF
  IF(p_nr > 0) THEN
    qnames(P_NR) = 'nr'; qdescp(P_NR) = 'Reain water concentrations (#/m3)'
    nscalar =  nscalar + 1
  END IF
  IF(p_ni > 0) THEN
    qnames(P_NI) = 'ni'; qdescp(P_NI) = 'Cloud ice concentrations (#/m3)'
    nscalar =  nscalar + 1
  END IF
  IF(p_ns > 0) THEN
    qnames(P_NS) = 'ns'; qdescp(P_NS) = 'Snow concentrations (#/m3)'
    nscalar =  nscalar + 1
  END IF
  IF(p_ng > 0) THEN
    qnames(P_NG) = 'ng'; qdescp(P_NG) = 'Graupel concentrations (#/m3)'
    nscalar =  nscalar + 1
  END IF
  IF(p_nh > 0) THEN
    qnames(P_NH) = 'nh'; qdescp(P_NH) = 'Hail concentrations (#/m3)'
    nscalar =  nscalar + 1
  END IF
  IF(p_zr > 0) THEN
    qnames(P_ZR) = 'zr'; qdescp(P_ZR) = 'Rain reflectivity (m6/m3)'
    nscalar =  nscalar + 1
  END IF
  IF(p_zi > 0) THEN
    qnames(P_ZI) = 'zi'; qdescp(P_ZI) = 'Ice reflectivity (m6/m3)'
    nscalar =  nscalar + 1
  END IF
  IF(p_zs > 0) THEN
    qnames(P_ZS) = 'zs'; qdescp(P_ZS) = 'Snow reflectivity (m6/m3)'
    nscalar =  nscalar + 1
  END IF
  IF(p_zg > 0) THEN
    qnames(P_ZG) = 'zg'; qdescp(P_ZG) = 'Graupel reflectivity (m6/m3)'
    nscalar =  nscalar + 1
  END IF
  IF(p_zh > 0) THEN
    qnames(P_ZH) = 'zh'; qdescp(P_ZH) = 'Hail reflectivity (m6/m3)'
    nscalar =  nscalar + 1
  END IF

  IF(p_cc > 0) THEN
    IF (lvldbg > 0) WRITE(*,'(1x,a)') 'Found CC ...'
    qnames(P_CC) = 'cc'; qdescp(P_CC) = 'Concentration (-)'
    nscalar = nscalar + 1
  END IF

  WRITE(6,'(1x,a,2(I3,a),/,1x,I3,2a)') 'The data file contains ',       &
         nscalarin,  ' microphysics variables, ',nscalar,               &
         ' of those will be read with the ARPS IO interface and ',      &
         nscalarq,' are moist mixing ratio.',                           &
         'Those variables are: '

  DO n = 1,nscalar
    WRITE(6,FMT='(5x,2a)',ADVANCE='NO') TRIM(qnames(n)),','
    IF (MOD(n,6) == 0) WRITE(6,*)
  END DO
  WRITE(6,'(/)')

  RETURN
END SUBROUTINE net_getdims
!
!##################################################################
!##################################################################
!######                                                      ######
!######             SUBROUTINE PAK_GETDIMS                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE pak_getdims(inch, nx,ny,nz,nzsoil,nstyps, ireturn)

!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!  Read in grid dimensions from base state/grid history data.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  7/17/2000.
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    inch     Channel number for binary reading.
!
!  OUTPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!    nzsoil   Number of soil levels
!    ireturn  Return status indicator
!             =0, successful read of all data
!             =1, error reading data
!             =2, end-of-file reached during read attempt
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx,ny,nz          ! Number of grid points in 3 directions
  INTEGER :: nzsoil            ! Number of soil levels
  INTEGER :: nstyps            ! Number of soil types
  INTEGER :: inch              ! Channel number for binary reading
  INTEGER :: ireturn           ! Return status indicator
  CHARACTER (LEN=10) :: tmunit
  CHARACTER (LEN=40) :: fmtverin
  INTEGER :: i,nocmnt
  REAL :: time
  CHARACTER (LEN=80) :: runname
  INTEGER :: bgrdin,bbasin,bvarin,mstin,brainin,bicein,btkein,          &
           btrbin,bsfcin,landin,totin,prcin,radin,flxin,                &
           snowcin,snowin
  INTEGER :: ihdlen,ilen
  PARAMETER (ihdlen=5000)
  INTEGER :: ihead(ihdlen)
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  READ(inch,ERR=110,END=120) ihead

!  CALL decdhdr(nx,ny,nz,ihdlen,ihead,ilen,time,                         &
!               bgrdin,bbasin,bvarin,mstin,brainin,bicein,btkein,        &
!               btrbin,bsfcin,landin,totin,prcin,radin,flxin,            &
!               snowcin,snowin,fmtverin,tmunit)
   CALL decdhdr

  WRITE (6,*) "nstyps not defined in data, set to 4"
  nstyps = 4

  ireturn =0

  RETURN

  110   CONTINUE
  WRITE(6,'(/a/)') ' Error reading data in PAK_GETDIMS.'
  ireturn=1
  RETURN

  120   CONTINUE
  WRITE(6,'(/a/)') ' End of file reached in PAK_GETDIMS.'
  ireturn=2

  RETURN
END SUBROUTINE pak_getdims
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE GRIB_GETDIMS               ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE grib_getdims(inch, nx,ny,nz,nzsoil,nstyps, ireturn)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!  Read in grid dimensions from base state/grid history data.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  7/17/2000.
!
!  MODIFICATION HISTORY:
!
!  1 June 2002 Eric Kemp
!  Added ability to read nstyps.
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    inch     Channel number for binary reading.
!
!  OUTPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!    nzsoil   Number of soil levels
!    nstyps   NUmber of soil types.
!
!    ireturn  Return status indicator
!             =0, successful read of all data
!             =1, error reading data
!             =2, end-of-file reached during read attempt
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: inch
  INTEGER :: nx,ny,nz          ! Number of grid points in 3 directions
  INTEGER :: nzsoil            ! Number of soil levels
  INTEGER :: nstyps            ! Number of soil types
  INTEGER :: ireturn           ! Return status indicator

  CHARACTER (LEN=10) :: tmunit
  CHARACTER (LEN=40) :: fmtverin
  INTEGER :: i,j,nocmnt
  REAL :: time
  CHARACTER (LEN=80) :: runname
  CHARACTER (LEN=80 ) :: cmnt(50)  ! String of comments on this model run
!
!-----------------------------------------------------------------------
!
!  Parameters for GRIB packing
!
!-----------------------------------------------------------------------
!
  INTEGER :: nbufsz
  PARAMETER ( nbufsz = 800000 )  ! Size of GRIB buffer
  CHARACTER (LEN=1) :: mgrib(nbufsz) ! Buffer to carry GRIB messages
  INTEGER :: npos,hhdlen
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  READ (inch,ERR=110,END=120) (mgrib(i),i=1,3)

  npos = 0

  hhdlen = ICHAR(mgrib(npos+1))*65536                                   &
         + ICHAR(mgrib(npos+2))*256                                     &
         + ICHAR(mgrib(npos+3))

  BACKSPACE (inch)

  READ (inch,ERR=110,END=120) (mgrib(i),i=1,hhdlen)

  npos = npos + 3

  DO i=1,40
    fmtverin(i:i) = mgrib(npos+i)
  END DO

  npos = npos + 40

  DO i=1,80
    runname(i:i) = mgrib(npos+i)
  END DO

  npos = npos + 80

  nocmnt = ICHAR(mgrib(npos+1))
  npos = npos + 1

  IF( nocmnt > 0 ) THEN
    DO j=1,nocmnt
      DO  i=1,80
        cmnt(j)(i:i) = mgrib(npos+i)
      END DO
      npos = npos + 80
    END DO
  END IF

  CALL ibm2flt( mgrib(npos+1), time )
  npos = npos + 4

  nx = ICHAR(mgrib(npos+1))*256 + ICHAR(mgrib(npos+2))
  npos = npos + 2
  ny = ICHAR(mgrib(npos+1))*256 + ICHAR(mgrib(npos+2))
  npos = npos + 2
  nz = ICHAR(mgrib(npos+1))*256 + ICHAR(mgrib(npos+2))
  npos = npos + 2
  nzsoil = ICHAR(mgrib(npos+1))*256 + ICHAR(mgrib(npos+2))
  npos = npos + 2


!  WRITE (6,*) "nstyps not defined in data, set to 4"
!  nstyps = 4

! Code extraction from RDHISHD

  npos = npos + 21 ! Skip mapproj and time

  npos = npos + 4 ! Skip umove
  npos = npos + 4 ! Skip vmove
  npos = npos + 4 ! Skip xgrdorg
  npos = npos + 4 ! Skip ygrdorg
  npos = npos + 4 ! Skip trulat1
  npos = npos + 4 ! Skip trulat2
  npos = npos + 4 ! Skip trulon
  npos = npos + 4 ! Skip sclfct
  npos = npos + 4 ! Skip rdummy
  npos = npos + 4 ! Skip rdummy
  npos = npos + 4 ! Skip rdummy
  npos = npos + 4 ! Skip rdummy
  npos = npos + 4 ! Skip rdummy
  npos = npos + 4 ! Skip rdummy
  npos = npos + 4 ! Skip rdummy
  npos = npos + 4 ! Skip tstop
  npos = npos + 4 ! Skip thisdmp
  npos = npos + 4 ! Skip latitud
  npos = npos + 4 ! Skip ctrlat
  npos = npos + 4 ! Skip ctrlon

  nstyps = ICHAR(mgrib(npos+1))
  IF ( nstyps < 0 ) THEN
    nstyps = 1
  END IF

  ireturn = 0

  RETURN

  110   CONTINUE
  WRITE(6,'(/a/)') ' Error reading data in GRIB_GETDIMS'
  ireturn=1
  RETURN

  120   CONTINUE
  WRITE(6,'(/a/)') ' End of file reached in GRIB_GETDIMS'
  ireturn=2
  RETURN

  RETURN
END SUBROUTINE grib_getdims
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE EDGFILL                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE edgfill(a,nx1,nx2,ibgn,iend,ny1,ny2,jbgn,jend,               &
           nz1,nz2,kbgn,kend)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Fill in the edges of a data array from the valid interior grid
!  points so that the arrays are completely filled. This is done primarily
!  for the benefit of the compression routine which must have a valid
!  value in each array location.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster
!  7/14/92.
!
!  9/1/94 (Y. Lu)
!  Cleaned up documentation.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    a        Array to have edges filled in
!
!    nx1:nx2  Dimensioned range of the  first index of array "a"
!    ny1:ny2  Dimensioned range of the second index of array "a"
!    nz1:nz2  Dimensioned range of the  third index of array "a"
!
!    ibgn,iend  Valued range of the  first index of array "a"
!    jgbn,jend  Valued range of the second index of array "a"
!    kbgn,kend  Valued range of the  third index of array "a"
!
  IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  INTEGER :: nx1,nx2,ny1,ny2,nz1,nz2
  REAL :: a(nx1:nx2,ny1:ny2,nz1:nz2)
  INTEGER :: ibgn,iend,jbgn,jend,kbgn,kend
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k
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
!  Fill in top and bottom grid levels, as needed.
!
!-----------------------------------------------------------------------
!
  IF( kbgn > nz1 ) THEN
    DO k=nz1,(kbgn-1)
      DO j=jbgn,jend
        DO i=ibgn,iend
          a(i,j,k)=a(i,j,kbgn)
        END DO
      END DO
    END DO
  END IF

  IF( kend < nz2 ) THEN
    DO k=(kend+1),nz2
      DO j=jbgn,jend
        DO i=ibgn,iend
          a(i,j,k)=a(i,j,kend)
        END DO
      END DO
    END DO
  END IF
!
!-----------------------------------------------------------------------
!
!  Fill in north and south edges, as needed.
!
!-----------------------------------------------------------------------
!
  IF( jbgn > ny1 ) THEN
    DO k=nz1,nz2
      DO j=ny1,(jbgn-1)
        DO i=ibgn,iend
          a(i,j,k)=a(i,jbgn,k)
        END DO
      END DO
    END DO
  END IF

  IF( jend < ny2 ) THEN
    DO k=nz1,nz2
      DO j=(jend+1),ny2
        DO i=ibgn,iend
          a(i,j,k)=a(i,jend,k)
        END DO
      END DO
    END DO
  END IF
!
!-----------------------------------------------------------------------
!
!  Fill in east and west edges, as needed.
!
!-----------------------------------------------------------------------
!
  IF( ibgn > nx1 ) THEN
    DO k=nz1,nz2
      DO j=ny1,ny2
        DO i=nx1,(ibgn-1)
          a(i,j,k)=a(ibgn,j,k)
        END DO
      END DO
    END DO
  END IF

  IF( iend < nx2 ) THEN
    DO k=nz1,nz2
      DO j=ny1,ny2
        DO i=(iend+1),nx2
          a(i,j,k)=a(iend,j,k)
        END DO
      END DO
    END DO
  END IF

  RETURN
END SUBROUTINE edgfill
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE IEDGFILL                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE iedgfill(a,nx1,nx2,ibgn,iend,ny1,ny2,jbgn,jend,              &
           nz1,nz2,kbgn,kend)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Fill in the edges of a data array from the valid interior grid
!  points so that the arrays are completely filled. This is done primarily
!  for the benefit of the compression routine which must have a valid
!  value in each array location.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster
!  7/14/92.
!
!  9/1/94 (Y. Lu)
!  Cleaned up documentation.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    a        Array to have edges filled in
!
!    nx1:nx2  Dimensioned range of the  first index of array "a"
!    ny1:ny2  Dimensioned range of the second index of array "a"
!    nz1:nz2  Dimensioned range of the  third index of array "a"
!
!    ibgn,iend  Valued range of the  first index of array "a"
!    jgbn,jend  Valued range of the second index of array "a"
!    kbgn,kend  Valued range of the  third index of array "a"
!
  IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  INTEGER :: nx1,nx2,ny1,ny2,nz1,nz2
  INTEGER :: a(nx1:nx2,ny1:ny2,nz1:nz2)
  INTEGER :: ibgn,iend,jbgn,jend,kbgn,kend
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k
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
!  Fill in top and bottom grid levels, as needed.
!
!-----------------------------------------------------------------------
!
  IF( kbgn > nz1 ) THEN
    DO k=nz1,(kbgn-1)
      DO j=jbgn,jend
        DO i=ibgn,iend
          a(i,j,k)=a(i,j,kbgn)
        END DO
      END DO
    END DO
  END IF

  IF( kend < nz2 ) THEN
    DO k=(kend+1),nz2
      DO j=jbgn,jend
        DO i=ibgn,iend
          a(i,j,k)=a(i,j,kend)
        END DO
      END DO
    END DO
  END IF
!
!-----------------------------------------------------------------------
!
!  Fill in north and south edges, as needed.
!
!-----------------------------------------------------------------------
!
  IF( jbgn > ny1 ) THEN
    DO k=nz1,nz2
      DO j=ny1,(jbgn-1)
        DO i=ibgn,iend
          a(i,j,k)=a(i,jbgn,k)
        END DO
      END DO
    END DO
  END IF

  IF( jend < ny2 ) THEN
    DO k=nz1,nz2
      DO j=(jend+1),ny2
        DO i=ibgn,iend
          a(i,j,k)=a(i,jend,k)
        END DO
      END DO
    END DO
  END IF
!
!-----------------------------------------------------------------------
!
!  Fill in east and west edges, as needed.
!
!-----------------------------------------------------------------------
!
  IF( ibgn > nx1 ) THEN
    DO k=nz1,nz2
      DO j=ny1,ny2
        DO i=nx1,(ibgn-1)
          a(i,j,k)=a(ibgn,j,k)
        END DO
      END DO
    END DO
  END IF

  IF( iend < nx2 ) THEN
    DO k=nz1,nz2
      DO j=ny1,ny2
        DO i=(iend+1),nx2
          a(i,j,k)=a(iend,j,k)
        END DO
      END DO
    END DO
  END IF

  RETURN
END SUBROUTINE iedgfill
!
!##################################################################
!##################################################################
!######                                                      ######
!######             SUBROUTINE GET_GRIDXYZP                  ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################

SUBROUTINE get_gridxyzzp(nx,ny,nz,filenamein,filefmt,nprocx_in,nprocy_in, &
                         x,y,z,zp,istat)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!  Read in grid variables from base state/grid history data.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster, after a similar program developed by Ming Xue
!  02/06/2005.
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    filename HDF file name of grid/base file.
!
!  OUTPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!    x        x-coordinate data
!    y        y-coordinate data
!    z        z-coordinate data
!    zp       zp-coordinate data
!
!  SCRATCH:
!    itmp     Temporary array for hdf compression
!    hmin     hmin temporary array
!    hmax     hmax temporary array
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

  CHARACTER (LEN=256), INTENT(IN) :: filenamein
  INTEGER, INTENT(IN) :: filefmt
  INTEGER, INTENT(IN) :: nprocx_in, nprocy_in

  INTEGER, INTENT(IN) :: nx,ny,nz          ! Number of grid points in 3 directions
  REAL,    TARGET     :: x(nx)
  REAL,    TARGET     :: y(ny)
  REAL,    TARGET     :: z(nz)
  REAL,    TARGET     :: zp(nx,ny,nz)

  INTEGER, INTENT(OUT) :: istat           ! Return status indicator

!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------

  INTEGER :: nxpatch, nypatch
  INTEGER :: nxdta, nydta

  REAL, POINTER :: xin(:), yin(:), zin(:), zpin(:,:,:)
!
!-----------------------------------------------------------------------
!
!  Include files
!
!-----------------------------------------------------------------------
!
  INCLUDE 'mp.inc'

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  IF (mp_opt > 0 .AND. (nprocx_in > 1 .OR. nprocy_in > 1)) THEN

    IF ( MOD(nprocx_in, nproc_x) /= 0 .OR. MOD(nprocy_in, nproc_y) /= 0) THEN
      IF (myproc == 0) WRITE(6,'(3x,a/,2(3x,2(a,I2)/))')                &
      'nprocx_in (nprocy_in) must be a multiple of nproc_x(nproc_y)', &
      'nprocx_in = ',nprocx_in, ', nprocy_in = ',nprocy_in,             &
      'nproc_x   = ',nproc_x,   ', nproc_y   = ', nproc_y
      CALL arpsstop('ERROR: nprocx_in, nprocy_in.',1)
    END IF
    nxpatch = nprocx_in / nproc_x
    nypatch = nprocy_in / nproc_y
  ELSE
    nxpatch = nprocx_in
    nypatch = nprocy_in
  END IF

  IF (mp_opt > 0 .AND. readsplit(FINDX_H) > 0) THEN
    nxdta = (nx-3)*nproc_x+3
    nydta = (ny-3)*nproc_y+3
  ELSE
    nxdta = (nx-3)/nxpatch+3
    nydta = (ny-3)/nypatch+3
  END IF

  IF (nxdta /= nx .OR. nydta /= ny) THEN
    ALLOCATE(xin(nxdta),           STAT=istat)
    ALLOCATE(yin(nydta),           STAT=istat)
    ALLOCATE(zin(nz),              STAT=istat)
    ALLOCATE(zpin(nxdta,nydta,nz), STAT=istat)
    CALL check_alloc_status(istat, "GET_GRID:zpin")
  ELSE
    xin  => x
    yin  => y
    zin  => z
    zpin => zp
  END IF

  IF (myproc == 0) WRITE(6,'(/1x,2a)') 'GET_GRIDXYZZP: Reading file: ', trim(filenamein)

  IF (mp_opt > 0 .AND. readsplit(FINDX_H) > 0) THEN  ! mpi split
    CALL get_grid_split(filenamein,filefmt,                             &
                        nxdta,nydta,nz,xin,yin,zin,zpin,                &
                        nx,ny,x,y,z,zp,istat)
  ELSE  ! nompi nojoin; nompi join; mpi nojoin, mpi join
    CALL get_grid_join (filenamein,filefmt,nxpatch,nypatch,             &
                        nxdta,nydta,nz,xin,yin,zin,zpin,                &
                        nx,ny,x,y,z,zp,istat)
  END IF

  IF (istat /= 0) CALL arpsstop('Error in get_gridxyzzp',1)

  IF (nxdta /= nx .OR. nydta /= ny) DEALLOCATE(xin, yin, zin, zpin)

  RETURN
END SUBROUTINE get_gridxyzzp

!
!##################################################################
!##################################################################
!######                                                      ######
!######             SUBROUTINE GET_GRID_JOIN                 ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################

SUBROUTINE get_grid_join(filenamein,filefmt,nxpatch,nypatch,            &
                         nxdta,nydta,nz,xin,yin,zin,zpin,               &
                         nx,ny,x,y,z,zp,ireturn)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!  Read in grid variables from base state/grid history data.
!
!  This subroutine handles the following scenarios
!
!  No-mpi mode
!      Read one file
!      Read patch files and join together
!
!  MPI mode
!      Read same patches
!      Read smaller patches and join together
!
!  If no join is needed, then variable ending with in (for example xin) should
!  point to the smae memorry as the corresponding variable (for example x).
!
!  If join is done, then variables ending with in are just temporary working
!  arrays, the actual returning variable are in the correponding variables
!  without "in".
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster, after a similar program developed by Ming Xue
!  02/06/2005.
!
!  MODIFICATION HISTORY:
!  03/29/2010  Keith Brewster
!  Moved checks for nxin, nyin, nzin earlier in logic to avoid i
!  possible problems reading-in some arrays.
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    filename HDF file name of grid/base file.
!
!  OUTPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!    x        x-coordinate data
!    y        y-coordinate data
!    z        z-coordinate data
!    zp       zp-coordinate data
!
!  SCRATCH:
!    itmp     Temporary array for hdf compression
!    hmin     hmin temporary array
!    hmax     hmax temporary array
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

  CHARACTER (LEN=256), INTENT(IN) :: filenamein
  INTEGER, INTENT(IN) :: filefmt
  INTEGER, INTENT(IN) :: nxpatch, nypatch

  INTEGER, INTENT(IN) :: nxdta,nydta,nz    ! Number of grid points in 3 directions
  REAL,   INTENT(OUT) :: xin(nxdta)
  REAL,   INTENT(OUT) :: yin(nydta)
  REAL,   INTENT(OUT) :: zin(nz)
  REAL,   INTENT(OUT) :: zpin(nxdta,nydta,nz)

  INTEGER, INTENT(IN) :: nx,ny          ! Number of grid points in 3 directions
  REAL,   INTENT(OUT) :: x(nx)
  REAL,   INTENT(OUT) :: y(ny)
  REAL,   INTENT(OUT) :: z(nz)
  REAL,   INTENT(OUT) :: zp(nx,ny,nz)

  INTEGER, INTENT(OUT) :: ireturn           ! Return status indicator

!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------

  INTEGER (KIND=selected_int_kind(4)), ALLOCATABLE :: itmp(:,:,:) ! Temporary array
  REAL, ALLOCATABLE :: hmin(:)
  REAL, ALLOCATABLE :: hmax(:)

  INTEGER, PARAMETER :: lchanl = 6
  INTEGER :: nchanl, sd_id, istat
  REAL    :: alatpro(2)
  REAL    :: sclf,dxscl,dyscl,ctrx,ctry,swx,swy

  INTEGER :: grdbas,totin
  INTEGER :: nxin, nyin, nzin, nzsoilin

  INTEGER :: idummy
  REAL    :: rdummy

  REAL    :: time
  CHARACTER (LEN=40) :: fmtverin
  CHARACTER (LEN=12) :: label
  CHARACTER (LEN=10) :: tmunit

  CHARACTER (LEN=256) :: filename
  INTEGER :: iloc, jloc
  INTEGER :: i,j,k,ia, ja

  LOGICAL :: success
!
!-----------------------------------------------------------------------
!
!  Include files
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'grid.inc'
  INCLUDE 'mp.inc'

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  IF (myproc == 0) WRITE(6,'(/1x,2a)') 'GET_GRID_JOIN: Reading file: ', trim(filenamein)

  DO jloc = 1, nypatch
    DO iloc = 1, nxpatch

      IF ( mp_opt > 0 .OR. (nxpatch > 1 .OR. nypatch > 1) ) THEN
        CALL gtsplitfn(filenamein,nxpatch,nypatch,loc_x,loc_y,iloc,jloc, &
                       0,0,1,lvldbg,filename,istat)
      ELSE
        WRITE(filename,'(a)') TRIM(filenamein)
      END IF

      IF (filefmt == 1) THEN

        grdbas = 1

        CALL getunit( nchanl )

        OPEN(UNIT=nchanl,FILE=TRIM(filename),                           &
             STATUS='old',FORM='unformatted',ACTION='READ',IOSTAT=istat)

        READ(nchanl,ERR=110,END=120) fmtverin

        IF (myproc == 0) WRITE(6,'(1x,a,a/)')                           &
          'Incoming data format, fmtverin = ',TRIM(fmtverin)

        READ(nchanl,ERR=110,END=120) runname

        IF (myproc == 0) WRITE(6,'(1x,2a,/)')                           &
          'THE NAME OF THIS RUN IS: ', TRIM(runname)

        READ(nchanl,ERR=110,END=120) nocmnt
        IF( nocmnt > 0 ) THEN
          DO i=1,nocmnt
            READ(nchanl,ERR=110,END=120) cmnt(i)
          END DO
        END IF

        IF( nocmnt > 0 .AND. myproc == 0 ) THEN
          DO i=1,nocmnt
            WRITE(6,'(1x,a)') cmnt(i)
          END DO
        END IF

        READ(nchanl,ERR=110,END=120) time,tmunit

        READ(nchanl,ERR=110,END=120) nxin, nyin, nzin,nzsoilin

        IF ( nxin /= nxdta .OR. nyin /= nydta .OR. nzin /= nz ) THEN
          IF (myproc == 0) THEN
            WRITE(6,'(1x,a)') ' Dimensions in GET_GRID_JOIN inconsistent with data.'
            WRITE(6,'(1x,a,3I15)') ' Read were: ',  nxin,  nyin, nzin
            WRITE(6,'(1x,a,3I15)') ' Expected:  ', nxdta, nydta, nz
            WRITE(6,'(1x,a)') ' Program aborted in GET_GRID_JOIN.'
            CALL flush(6)
          END IF
          CALL arpsstop('arpsstop called from GET_GRID_JOIN due to nxin...',1)
        END IF

        READ(nchanl,ERR=110,END=120)                                    &
                  idummy,idummy,idummy,idummy,idummy,                   &
                  idummy,idummy,idummy,idummy,totin,                    &
                  idummy,idummy,idummy,mapproj, month,                  &
                  day, year, hour, minute, second

        READ(nchanl,ERR=110,END=120)                                    &
                  rdummy, rdummy, rdummy, rdummy, trulat1,              &
                  trulat2,trulon, sclfct, rdummy, rdummy,               &
                  rdummy, rdummy, rdummy, rdummy, rdummy,               &
                  rdummy, rdummy, rdummy, ctrlat, ctrlon

        IF ( totin /= 0 ) THEN
          READ(nchanl,ERR=110,END=120)                                  &
                  idummy,idummy,idummy,idummy,idummy,                   &
                  idummy,idummy,idummy,idummy,idummy,                   &
                  idummy,idummy,idummy,idummy,idummy,                   &
                  idummy,idummy,idummy,idummy,idummy

          READ(nchanl,ERR=110,END=120)                                  &
                  rdummy,rdummy,rdummy,rdummy,rdummy,                   &
                  rdummy,rdummy,rdummy,rdummy,rdummy,                   &
                  rdummy,rdummy,rdummy,rdummy,rdummy,                   &
                  rdummy,rdummy,rdummy,rdummy,rdummy
        END IF

        IF( grdbas == 1 ) THEN
          READ(nchanl,ERR=110,END=120) label
          READ(nchanl,ERR=110,END=120) xin
          IF (myproc == 0) WRITE(lchanl,910) label,' x.'

          READ(nchanl,ERR=110,END=120) label
          READ(nchanl,ERR=110,END=120) yin
          IF (myproc == 0) WRITE(lchanl,910) label,' y.'

          READ(nchanl,ERR=110,END=120) label
          READ(nchanl,ERR=110,END=120) zin
          IF (myproc == 0) WRITE(lchanl,910) label,' z.'

          READ(nchanl,ERR=110,END=120) label
          READ(nchanl,ERR=110,END=120) zpin
          IF (myproc == 0) WRITE(lchanl,910) label,' zp.'

        ELSE

          WRITE(6,'(1x,2(a/))') 'ERROR: It is not a grid base file.','Program stoping ...'
          CALL arpsstop('Expect a grid & base file inside GET_GRID_JOIN.',1)

        END IF  ! grdbas

        GO TO 130

        110   CONTINUE
        WRITE(6,'(/a/)') ' Error reading data in GET_GRID_JOIN'
        ireturn=-1
        CALL arpsstop('get_grid_join error.',1)

        120   CONTINUE
        WRITE(6,'(/a/)') ' End of file reached in GET_GRID_JOIN'
        ireturn=-2
        CALL arpsstop('get_grid_join end reached.',1)

        130 CONTINUE
        CLOSE(UNIT=nchanl)
        CALL retunit( nchanl )
        WRITE(6,*)

        910   FORMAT(1X,'Field ',a12,' was read into array',a)

      ELSE IF (filefmt == 3) THEN

        CALL hdfopen(filename,1,sd_id)

        IF (sd_id < 0) THEN
          WRITE (6,*) "get_grid_join: ERROR opening ",              &
                      trim(filename)," for reading."
          ireturn=-1
          RETURN
        ELSE
          WRITE(6,*) 'File ',TRIM(filename),' opened.'
        END IF

        success=.true.

        IF (iloc == 1 .AND. jloc == 1) THEN
          CALL hdfrdc(sd_id,80,"runname",runname,istat)
          CALL hdfrdi(sd_id,"mapproj",mapproj,istat)
          CALL hdfrdr(sd_id,"trulat1",trulat1,istat)
          CALL hdfrdr(sd_id,"trulat2",trulat2,istat)
          CALL hdfrdr(sd_id,"trulon",trulon,istat)
          CALL hdfrdr(sd_id,"sclfct",sclfct,istat)
          CALL hdfrdr(sd_id,"ctrlat",ctrlat,istat)
          CALL hdfrdr(sd_id,"ctrlon",ctrlon,istat)
          CALL hdfrdi(sd_id,"year",year,istat)
          CALL hdfrdi(sd_id,"month",month,istat)
          CALL hdfrdi(sd_id,"day",day,istat)
          CALL hdfrdi(sd_id,"hour",hour,istat)
          CALL hdfrdi(sd_id,"minute",minute,istat)
          CALL hdfrdi(sd_id,"second",second,istat)
        END IF

        CALL hdfrdi(sd_id,"nx",nxin,istat)
        CALL hdfrdi(sd_id,"ny",nyin,istat)
        CALL hdfrdi(sd_id,"nz",nzin,istat)

        IF ( nxin /= nxdta .OR. nyin /= nydta .OR. nzin /= nz ) THEN
          IF (myproc == 0) THEN
            WRITE(6,'(1x,a)') ' Dimensions in GET_GRID_JOIN inconsistent with data.'
            WRITE(6,'(1x,a,3I15)') ' Read were: ',  nxin,  nyin, nzin
            WRITE(6,'(1x,a,3I15)') ' Expected:  ', nxdta, nydta, nz
            WRITE(6,'(1x,a)') ' Program aborted in GET_GRID_JOIN.'
            CALL flush(6)
          END IF
          CALL arpsstop('arpsstop called from GET_GRID_JOIN due to nxin...',1)
        END IF

        CALL hdfrd1d(sd_id,"x",nxin,xin,istat)
        IF (istat /= 0) success=.false.

        CALL hdfrd1d(sd_id,"y",nyin,yin,istat)
        IF (istat /= 0) success=.false.

        CALL hdfrd1d(sd_id,"z",nz,zin,istat)
        IF (istat /= 0) success=.false.

        IF (.NOT. ALLOCATED(itmp)) THEN
          ALLOCATE(itmp(nxin,nyin,nz), STAT = istat)
          ALLOCATE(hmin(nz),       STAT = istat)
          ALLOCATE(hmax(nz),       STAT = istat)
        END IF

        CALL hdfrd3d(sd_id,"zp",nxin,nyin,nz,zpin,istat,itmp,hmax,hmin)
        IF (istat /= 0) success=.false.

        IF (success) THEN
          ireturn = 0
        ELSE
          WRITE(6,'(/a/)') ' Error reading data in GET_GRID_JOIN.'
          ireturn=-2
          CALL arpsstop('get_grid_join ERROR.',1)
        END IF

        CALL hdfclose(sd_id,istat)

      ELSE  ! other format
        WRITE(6,'(1x,a,I2,a,/)') 'The file formt - ',filefmt,' is not yet supported.'
        ireturn = -9
        RETURN
      END IF

      IF (nxpatch > 1 .OR. nypatch > 1) THEN
        DO i = 1, nxin
          ia = (iloc-1)*(nxin-3)+i
          x(ia) = xin(i)
        END DO

        DO j = 1, nyin
          ja = (jloc-1)*(nyin-3)+j
          y(ja) = yin(j)
        END DO

        DO k = 1, nzin
          z(k) = zin(k)
        END DO

        DO k = 1,nz
          DO j = 1, nyin
            ja = (jloc-1)*(nyin-3)+j
            DO i = 1, nxin
              ia = (iloc-1)*(nxin-3)+i
              zp(ia,ja,k) = zpin(i,j,k)
            END DO
          END DO
        END DO
      END IF

    END DO
  END DO

  CALL gtlfnkey( runname, lfnkey )

  mgrid = 1    ! Grid number 1
  nestgrd = 0  ! Not nested grid data

  IF (ALLOCATED(itmp)) DEALLOCATE(itmp, hmin, hmax)

  RETURN
END SUBROUTINE get_grid_join
!
!##################################################################
!##################################################################
!######                                                      ######
!######             SUBROUTINE GET_GRIDXYZP                  ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################

SUBROUTINE get_grid_split(filename,filefmt,                             &
                          nxlg,nylg,nz,xin,yin,zin,zpin,                &
                          nx,ny,x,y,z,zp,ireturn)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!  Read in grid variables from base state/grid history data.
!
!  This subroutine handles the scenario when reading one join file in
!  MPI mode. The reading data is split and propagate to all processors.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster, after a similar program developed by Ming Xue
!  02/06/2005.
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    filename HDF file name of grid/base file.
!
!  OUTPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!    x        x-coordinate data
!    y        y-coordinate data
!    z        z-coordinate data
!    zp       zp-coordinate data
!
!  SCRATCH:
!    itmp     Temporary array for hdf compression
!    hmin     hmin temporary array
!    hmax     hmax temporary array
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

  CHARACTER (LEN=256), INTENT(IN) :: filename
  INTEGER, INTENT(IN) :: filefmt

  INTEGER,  INTENT(IN) :: nxlg,nylg,nz          ! Number of grid points in 3 directions
  REAL,    INTENT(OUT) :: xin(nxlg)
  REAL,    INTENT(OUT) :: yin(nylg)
  REAL,    INTENT(OUT) :: zin(nz)
  REAL,    INTENT(OUT) :: zpin(nxlg,nylg,nz)

  INTEGER, INTENT(IN)  :: nx,ny          ! Number of grid points in 3 directions
  REAL,    INTENT(OUT) :: x(nx)
  REAL,    INTENT(OUT) :: y(ny)
  REAL,    INTENT(OUT) :: z(nz)
  REAL,    INTENT(OUT) :: zp(nx,ny,nz)

  INTEGER, INTENT(OUT) :: ireturn           ! Return status indicator

!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------

  INTEGER (KIND=selected_int_kind(4)), ALLOCATABLE :: itmp(:,:,:) ! Temporary array
  REAL, ALLOCATABLE :: hmin(:)
  REAL, ALLOCATABLE :: hmax(:)

  INTEGER, PARAMETER :: lchanl = 6
  INTEGER :: nchanl, sd_id, istat
  REAL    :: alatpro(2)
  REAL    :: sclf,dxscl,dyscl,ctrx,ctry,swx,swy

  INTEGER :: grdbas,totin
  INTEGER :: nxin, nyin, nzin, nzsoilin

  INTEGER :: idummy
  REAL    :: rdummy

  REAL    :: time
  CHARACTER (LEN=40) :: fmtverin
  CHARACTER (LEN=12) :: label
  CHARACTER (LEN=10) :: tmunit

  INTEGER :: i,j,k

  LOGICAL :: success
!
!-----------------------------------------------------------------------
!
!  Include files
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'grid.inc'
  INCLUDE 'mp.inc'

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  IF (myproc == 0) WRITE(6,'(/1x,2a)') 'get_grid_split: Reading file: ', trim(filenamein)

  IF (myproc == 0 ) THEN

    IF (filefmt == 1) THEN

      grdbas = 1

      CALL getunit( nchanl )

      OPEN(UNIT=nchanl,FILE=TRIM(filename),                           &
           STATUS='old',FORM='unformatted',ACTION='READ',IOSTAT=istat)

      READ(nchanl,ERR=110,END=120) fmtverin

      IF (myproc == 0) WRITE(6,'(1x,a,a/)')                           &
        'Incoming data format, fmtverin = ',TRIM(fmtverin)

      READ(nchanl,ERR=110,END=120) runname

      IF (myproc == 0) WRITE(6,'(1x,2a,/)')                           &
        'THE NAME OF THIS RUN IS: ', TRIM(runname)

      READ(nchanl,ERR=110,END=120) nocmnt
      IF( nocmnt > 0 ) THEN
        DO i=1,nocmnt
          READ(nchanl,ERR=110,END=120) cmnt(i)
        END DO
      END IF

      IF( nocmnt > 0 .AND. myproc == 0 ) THEN
        DO i=1,nocmnt
          WRITE(6,'(1x,a)') cmnt(i)
        END DO
      END IF

      READ(nchanl,ERR=110,END=120) time,tmunit

      READ(nchanl,ERR=110,END=120) nxin, nyin, nzin,nzsoilin

      READ(nchanl,ERR=110,END=120)                                    &
                idummy,idummy,idummy,idummy,idummy,                   &
                idummy,idummy,idummy,idummy,totin,                    &
                idummy,idummy,idummy,mapproj, month,                  &
                day, year, hour, minute, second

      READ(nchanl,ERR=110,END=120)                                    &
                rdummy, rdummy, rdummy, rdummy, trulat1,              &
                trulat2,trulon, sclfct, rdummy, rdummy,               &
                rdummy, rdummy, rdummy, rdummy, rdummy,               &
                rdummy, rdummy, rdummy, ctrlat, ctrlon

      IF ( totin /= 0 ) THEN
        READ(nchanl,ERR=110,END=120)                                  &
                idummy,idummy,idummy,idummy,idummy,                   &
                idummy,idummy,idummy,idummy,idummy,                   &
                idummy,idummy,idummy,idummy,idummy,                   &
                idummy,idummy,idummy,idummy,idummy

        READ(nchanl,ERR=110,END=120)                                  &
                rdummy,rdummy,rdummy,rdummy,rdummy,                   &
                rdummy,rdummy,rdummy,rdummy,rdummy,                   &
                rdummy,rdummy,rdummy,rdummy,rdummy,                   &
                rdummy,rdummy,rdummy,rdummy,rdummy
      END IF

      IF( grdbas == 1 ) THEN
        READ(nchanl,ERR=110,END=120) label
        READ(nchanl,ERR=110,END=120) xin
        IF (myproc == 0) WRITE(lchanl,910) label,' x.'

        READ(nchanl,ERR=110,END=120) label
        READ(nchanl,ERR=110,END=120) yin
        IF (myproc == 0) WRITE(lchanl,910) label,' y.'

        READ(nchanl,ERR=110,END=120) label
        READ(nchanl,ERR=110,END=120) zin
        IF (myproc == 0) WRITE(lchanl,910) label,' z.'

        READ(nchanl,ERR=110,END=120) label
        READ(nchanl,ERR=110,END=120) zpin
        IF (myproc == 0) WRITE(lchanl,910) label,' zp.'

      ELSE

        WRITE(6,'(1x,2(a/))') 'ERROR: It is not a grid base file.','Program stoping ...'
!        CALL arpsstop('Expect a grid & base file inside get_grid_split.',1)
        GO TO 999

      END IF  ! grdbas

      GO TO 130

      110   CONTINUE
      WRITE(6,'(/a/)') ' Error reading data in READUV'
      ireturn=-1
      GO TO 999
!      CALL arpsstop('get_grid error.',1)

      120   CONTINUE
      WRITE(6,'(/a/)') ' End of file reached in READUV'
      ireturn=-2
      GO TO 999
!      CALL arpsstop('get_grid end reached.',1)

      130 CONTINUE
      CLOSE(UNIT=nchanl)
      CALL retunit( nchanl )
      WRITE(6,*)

      910   FORMAT(1X,'Field ',a12,' was read into array',a)

    ELSE IF (filefmt == 3) THEN

      CALL hdfopen(filename,1,sd_id)

      IF (sd_id < 0) THEN
        WRITE (6,*) "get_grid_split: ERROR opening ",              &
                    trim(filename)," for reading."
        ireturn=-1
        GO TO 999
!        RETURN
      ELSE
        WRITE(6,*) 'File ',TRIM(filename),' opened.'
      END IF

      success=.true.

      CALL hdfrdc(sd_id,80,"runname",runname,istat)
      CALL hdfrdi(sd_id,"mapproj",mapproj,istat)
      CALL hdfrdr(sd_id,"trulat1",trulat1,istat)
      CALL hdfrdr(sd_id,"trulat2",trulat2,istat)
      CALL hdfrdr(sd_id,"trulon",trulon,istat)
      CALL hdfrdr(sd_id,"sclfct",sclfct,istat)
      CALL hdfrdr(sd_id,"ctrlat",ctrlat,istat)
      CALL hdfrdr(sd_id,"ctrlon",ctrlon,istat)
      CALL hdfrdi(sd_id,"year",year,istat)
      CALL hdfrdi(sd_id,"month",month,istat)
      CALL hdfrdi(sd_id,"day",day,istat)
      CALL hdfrdi(sd_id,"hour",hour,istat)
      CALL hdfrdi(sd_id,"minute",minute,istat)
      CALL hdfrdi(sd_id,"second",second,istat)

      CALL hdfrdi(sd_id,"nx",nxin,istat)
      CALL hdfrdi(sd_id,"ny",nyin,istat)
      CALL hdfrdi(sd_id,"nz",nzin,istat)

      CALL hdfrd1d(sd_id,"x",nxin,xin,istat)
      IF (istat /= 0) success=.false.

      CALL hdfrd1d(sd_id,"y",nyin,yin,istat)
      IF (istat /= 0) success=.false.

      CALL hdfrd1d(sd_id,"z",nz,zin,istat)
      IF (istat /= 0) success=.false.

      IF (.NOT. ALLOCATED(itmp)) THEN
        ALLOCATE(itmp(nxin,nyin,nz), STAT = istat)
        ALLOCATE(hmin(nz),       STAT = istat)
        ALLOCATE(hmax(nz),       STAT = istat)
      END IF

      CALL hdfrd3d(sd_id,"zp",nxin,nyin,nz,zpin,istat,itmp,hmax,hmin)
      IF (istat /= 0) success=.false.

      IF (success) THEN
        ireturn = 0
      ELSE
        WRITE(6,'(/a/)') ' Error reading data in GET_GRID_SPLIT.'
        ireturn=-2
        GO TO 999
!        CALL arpsstop('get_grid_split ERROR.',1)
      END IF

      CALL hdfclose(sd_id,istat)

    ELSE  ! other format
      WRITE(6,'(1x,a,I2,a,/)') 'The file formt - ',filefmt,' is not yet supported.'
      ireturn = -9
      GO TO 999
      !RETURN
    END IF

    IF ( nxin /= nxlg .OR. nyin /= nylg .OR. nzin /= nz ) THEN
      WRITE(6,'(1x,a)') ' Dimensions in get_grid_split inconsistent with data.'
      WRITE(6,'(1x,a,3I15)') ' Read were: ',  nxin,  nyin, nzin
      WRITE(6,'(1x,a,3I15)') ' Expected:  ', nxlg, nylg, nz
      WRITE(6,'(1x,a)') ' Program aborted in get_grid_split.'
      CALL flush(6)
      ireturn = -10
      !CALL arpsstop('arpsstop called from get_grid_split due to nxin...',1)
    END IF

    999 CONTINUE
    IF (ALLOCATED(itmp)) DEALLOCATE(itmp, hmin, hmax)
  END IF
  CALL mpupdatei(ireturn,1)
  IF (ireturn /= 0) CALL arpsstop('arpsstop called from get_grid_split',1)

  CALL mpupdatec(runname,80)
  CALL mpupdatei(mapproj,1)
  CALL mpupdater(trulat1,1)
  CALL mpupdater(trulat2,1)
  CALL mpupdater(trulon,1)
  CALL mpupdater(sclfct,1)
  CALL mpupdater(ctrlat,1)
  CALL mpupdater(ctrlon,1)
  CALL mpupdatei(year,1)
  CALL mpupdatei(month,1)
  CALL mpupdatei(day,1)
  CALL mpupdatei(hour,1)
  CALL mpupdatei(minute,1)
  CALL mpupdatei(second,1)

  CALL gtlfnkey( runname, lfnkey )

  mgrid = 1    ! Grid number 1
  nestgrd = 0  ! Not nested grid data

  CALL mpisplit1dx(xin,nx,x)
  CALL mpisplit1dy(yin,ny,y)
  z = zin
  CALL mpupdater(z,nz)
  CALL mpisplit3d(zpin,nx,ny,nz,zp)

  RETURN
END SUBROUTINE get_grid_split
