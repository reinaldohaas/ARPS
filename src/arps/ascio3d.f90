!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE ASCREAD                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE ascread(nx,ny,nz,nzsoil,nstyps,grdbas,inch,time,x,y,z,zp,    &
           zpsoil,uprt, vprt, wprt, ptprt, pprt,                        &
           qvprt, qc, qr, qi, qs, qh, tke,kmh,kmv,                      &
           ubar,vbar,wbar,ptbar,pbar,rhobar,qvbar,                      &
           soiltyp,stypfrct,vegtyp,lai,roufns,veg,                      &
           tsoil,qsoil,wetcanp,snowdpth,                                &
           raing,rainc,prcrate,                                         &
           radfrc,radsw,rnflx,radswnet,radlwin,                         &
           usflx,vsflx,ptsflx,qvsflx,                                   &
           ireturn)

!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Read history data from channel nchanl in ASCII format.
!  format.
!
!  All data read in are located at the original staggered grid points.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster
!  6/02/92.
!
!  MODIFICATION HISTORY:
!
!  7/14/92 (K. Brewster)
!  Added runname, comment and version number reading
!
!  8/20/92 (M. Xue)
!  Added data reading of computational z coordinate array z.
!
!  4/23/93 (M. Xue)
!  New data format.
!
!  02/06/95 (Y. Liu)
!  Added map projection parameters into the ASCII dumping
!
!  05/31/95 (Y. Liu)
!  Changed the integer dumping format from 20I3 to 10I8.
!
!  12/09/1998 (Donghai Wang)
!  Added the snow cover.
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
!                          parameter setting.
!    inch     Channel number for ASCII reading.
!             This channel must be opened for formatted reading
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
!    radswnet  Reflected shortwave radiation
!    radlwin  Incoming longwave radiation
!
!    usflx    Surface flux of u-momentum (kg/(m*s**2))
!    vsflx    Surface flux of v-momentum (kg/(m*s**2))
!    ptsflx   Surface heat flux (K*kg/(m**2 * s ))
!    qvsflx   Surface moisture flux of (kg/(m**2 * s))
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

  INTEGER :: grdbas            ! Data read flag.
  INTEGER :: inch              ! Channel number for binary reading
  REAL :: time                 ! Time in seconds of data read
                               ! from "filename"
  REAL :: x     (nx)           ! x-coord. of the physical and compu
                               ! -tational grid. Defined at u-point(m).
  REAL :: y     (ny)           ! y-coord. of the physical and compu
                               ! -tational grid. Defined at v-point(m).
  REAL :: z     (nz)           ! z-coord. of the computational grid.
                               ! Defined at w-point on the staggered
                               ! grid(m).
  REAL :: zp    (nx,ny,nz)     ! Physical height coordinate defined at
                               ! w-point of the staggered grid(m).
  REAL :: zpsoil(nx,ny,nzsoil) ! Physical height coordinate defined at
                               ! w-point in the soil (m).

  REAL :: uprt  (nx,ny,nz)     ! Perturbation u-velocity (m/s)
  REAL :: vprt  (nx,ny,nz)     ! Perturbation v-velocity (m/s)
  REAL :: wprt  (nx,ny,nz)     ! Perturbation w-velocity (m/s)
  REAL :: ptprt (nx,ny,nz)     ! Perturbation potential temperature (K)
  REAL :: pprt  (nx,ny,nz)     ! Perturbation pressure (Pascal)

  REAL :: qvprt (nx,ny,nz)     ! Water vapor mixing ratio (kg/kg)
  REAL :: qc    (nx,ny,nz)     ! Cloud water mixing ratio (kg/kg)
  REAL :: qr    (nx,ny,nz)     ! Rain water mixing ratio (kg/kg)
  REAL :: qi    (nx,ny,nz)     ! Cloud ice mixing ratio (kg/kg)
  REAL :: qs    (nx,ny,nz)     ! Snow mixing ratio (kg/kg)
  REAL :: qh    (nx,ny,nz)     ! Hail mixing ratio (kg/kg)
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
  REAL :: qvbar (nx,ny,nz)     ! Base state water vapor mixing ratio

  INTEGER :: nstyps             ! Number of soil type

  INTEGER :: soiltyp(nx,ny,nstyps)    ! Soil type
  REAL :: stypfrct(nx,ny,nstyps)    ! Soil type
  INTEGER :: vegtyp (nx,ny)            ! Vegetation type
  REAL :: lai    (nx,ny)            ! Leaf Area Index
  REAL :: roufns (nx,ny)            ! Surface roughness
  REAL :: veg    (nx,ny)            ! Vegetation fraction

  REAL :: tsoil (nx,ny,nzsoil,0:nstyps) ! Soil temperature (K)
  REAL :: qsoil (nx,ny,nzsoil,0:nstyps) ! Soil moisture (m**3/m**3)
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
  REAL :: radswnet (nx,ny)     ! Net solar radiation
  REAL :: radlwin (nx,ny)      ! Incoming longwave radiation

  REAL :: usflx (nx,ny)        ! Surface flux of u-momentum (kg/(m*s**2))
  REAL :: vsflx (nx,ny)        ! Surface flux of v-momentum (kg/(m*s**2))
  REAL :: ptsflx(nx,ny)        ! Surface heat flux (K*kg/(m*s**2))
  REAL :: qvsflx(nx,ny)        ! Surface moisture flux (kg/(m**2*s))

  INTEGER :: ireturn           ! Return status indicator
!
!-----------------------------------------------------------------------
!
!  Parameters describing routine that wrote the gridded data
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=40) :: fmtver0,fmtver1,fmtverin
  PARAMETER (fmtver0='004.10 ASCII Formatted Data')
  PARAMETER (fmtver1='004.10 ASCII Formatted Data')
  CHARACTER (LEN=10) :: tmunit
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
  CHARACTER (LEN=12) :: label
  INTEGER :: nxin,nyin,nzin
  INTEGER :: nzsoilin

  INTEGER :: bgrdin,bbasin,bvarin,bicein,btkein,btrbin,idummy
  REAL :: rdummy
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
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

!
!-----------------------------------------------------------------------
!
!  Read header info
!
!-----------------------------------------------------------------------
!
  READ(inch,'(1x,a40)',ERR=110,END=120) fmtverin

  IF( fmtverin /= fmtver0 .AND. fmtverin /= fmtver1 ) THEN
    WRITE(6,'(/1x,a/1x,2a/1x,2a/1x,2a/1x,a)')                           &
        'Data format incompatible with the data reader.',               &
        'Format of data is ',fmtverin,' Format of reader is ',fmtver1,  &
        'compitable to ',fmtver0, '. Job stopped.'
    CALL arpsstop('arpsstop called from asciiread',1)
  END IF

  READ(inch,'(1x,a80)',ERR=110,END=120) runname
  WRITE(6,'(//''  THE NAME OF THIS RUN IS:  '',A//)') runname

  READ(inch,'(1x,i4)',ERR=110,END=120) nocmnt
  IF( nocmnt > 0 ) THEN
    DO i=1,nocmnt
      READ(inch,'(1x,a80)',ERR=110,END=120) cmnt(i)
    END DO
  END IF

  IF( nocmnt > 0 ) THEN
    DO i=1,nocmnt
      WRITE(6,'(1x,a)') cmnt(i)
    END DO
  END IF
!
  READ(inch,'(1x,e16.8,1x,a10)',ERR=110,END=120) time,tmunit
!
!-----------------------------------------------------------------------
!
!  Get dimensions of data in ASCII file and check against
!     the dimensions passed to ASCREAD
!
!-----------------------------------------------------------------------
!
  READ(inch,'(1x, 4i12)',ERR=110,END=120) nxin,nyin,nzin,nzsoilin 
  IF( nxin /= nx .OR. nyin /= ny .OR. nzin /= nz .OR. nzsoilin /= nzsoil) &
          THEN
    WRITE(6,'(1x,a)')                                                   &
         ' Dimensions in ASCREAD inconsistent with data.'
    WRITE(6,'(1x,a,4I15)') ' Read were: ', nxin, nyin, nzin, nzsoilin 
    WRITE(6,'(1x,a)')                                                   &
         ' Program aborted in ASCREAD.'
    CALL arpsstop('arpsstop called from asciiread while reading nx..',1)
  END IF
!
!-----------------------------------------------------------------------
!
!  Read in flags for different data groups
!
!-----------------------------------------------------------------------
!

  IF( grdbas == 1 ) THEN ! Read grid and base state arrays

    WRITE(lchanl,'(1x,a,f8.1,a,f8.3,a/)')                               &
         'To read grid and base state data at time ', time,             &
         ' secs = ',(time/60.),' mins.'

    READ(inch,'(1x,10i8)',ERR=110,END=120)                              &
         bgrdin,bbasin,bvarin,mstin,bicein,                             &
         btrbin,idummy,idummy,landin,totin,                             &
         btkein,idummy,idummy,mapproj,month,                            &
         day,year,hour,minute,second

  ELSE ! Normal data reading

    WRITE(lchanl,'(1x,a,f8.1,a,f8.3,a/)')'To read data for time:',      &
         time,' secs = ',(time/60.),' mins.'

    READ(inch,'(1x,10i8)',ERR=110,END=120)                              &
         grdin,basin,varin,mstin,icein,                                 &
         trbin,sfcin,rainin,landin,totin,                               &
         tkein,idummy,idummy,mapproj,month,                             &
         day,year,hour,minute,second

  END IF

  READ(inch,910,ERR=110,END=120)                                        &
      umove,vmove,xgrdorg,ygrdorg,trulat1,                              &
      trulat2,trulon,sclfct,rdummy,rdummy,                              &
      rdummy,rdummy,rdummy,rdummy,rdummy,                               &
      tstop,thisdmp,latitud,ctrlat,ctrlon

  IF ( totin /= 0 ) THEN
!
!-----------------------------------------------------------------------
!
!  Read in additional parameters for ARPS history dump 4.0 or later
!  version.
!
!-----------------------------------------------------------------------
!
    READ(inch,'(1x,10i8)',ERR=110,END=120)                              &
         nstyp,  prcin, radin, flxin,snowcin,                           &
         snowin,idummy,idummy,idummy,idummy,                            &
         idummy,idummy,idummy,idummy,idummy,                            &
         idummy,idummy,idummy,idummy,idummy

    IF ( nstyp < 1 ) THEN
      nstyp = 1
    END IF

    READ(inch,910,ERR=110,END=120)                                      &
         rdummy,rdummy,rdummy,rdummy,rdummy,                            &
         rdummy,rdummy,rdummy,rdummy,rdummy,                            &
         rdummy,rdummy,rdummy,rdummy,rdummy,                            &
         rdummy,rdummy,rdummy,rdummy,rdummy
  END IF
!
!-----------------------------------------------------------------------
!
!  Read in x,y and z at grid cell centers (scalar points).
!
!----------------------------------------------------------------------
!
  IF( grdin == 1 .OR. grdbas == 1 ) THEN
    READ(inch,900,ERR=110,END=120) label
    READ(inch,910,ERR=110,END=120) x
    WRITE(lchanl,920) label,' x.'

    READ(inch,900,ERR=110,END=120) label
    READ(inch,910,ERR=110,END=120) y
    WRITE(lchanl,920) label,' y.'

    READ(inch,900,ERR=110,END=120) label
    READ(inch,910,ERR=110,END=120) z
    WRITE(lchanl,920) label,' z.'

    READ(inch,900,ERR=110,END=120) label
    READ(inch,910,ERR=110,END=120) zp
    WRITE(lchanl,920) label,' zp.'

    READ(inch,900,ERR=110,END=120) label
    READ(inch,910,ERR=110,END=120) zpsoil
    WRITE(lchanl,920) label,' zpsoil.'

  END IF
!
!-----------------------------------------------------------------------
!
!  Read in base state fields
!
!----------------------------------------------------------------------
!
  IF( basin == 1 .OR. grdbas == 1 ) THEN

    READ(inch,900,ERR=110,END=120) label
    READ(inch,910,ERR=110,END=120) ubar
    WRITE(lchanl,920) label,' ubar.'

    READ(inch,900,ERR=110,END=120) label
    READ(inch,910,ERR=110,END=120) vbar
    WRITE(lchanl,920) label,' vbar.'

    READ(inch,900,ERR=110,END=120) label
    READ(inch,910,ERR=110,END=120) wbar
    WRITE(lchanl,920) label,' wbar.'

    READ(inch,900,ERR=110,END=120) label
    READ(inch,910,ERR=110,END=120) ptbar
    WRITE(lchanl,920) label,' ptbar.'

    READ(inch,900,ERR=110,END=120) label
    READ(inch,910,ERR=110,END=120) pbar
    WRITE(lchanl,920) label,' pbar.'

    IF( mstin == 1 ) THEN
      READ(inch,900,ERR=110,END=120) label
      READ(inch,910,ERR=110,END=120) qvbar
      WRITE(lchanl,920) label,' qvbar.'
    END IF

    IF (landin == 1) THEN

      IF (nstyp <= 1) THEN

        READ(inch,900,ERR=110,END=120) label
        READ(inch,911,ERR=110,END=120)                                  &
              ((soiltyp(i,j,1),i=1,nx),j=1,ny)
        WRITE(lchanl,920) label,' soiltyp.'

      ELSE

        DO is=1,nstyp
          READ(inch,900,ERR=110,END=120) label
          READ(inch,911,ERR=110,END=120)                                &
                ((soiltyp(i,j,is),i=1,nx),j=1,ny)
          WRITE(lchanl,920) label,' soiltyp.'

          READ(inch,900,ERR=110,END=120) label
          READ(inch,910,ERR=110,END=120)                                &
                ((stypfrct(i,j,is),i=1,nx),j=1,ny)
          WRITE(lchanl,920) label,'stypfrct.'
        END DO

      END IF

      READ(inch,900,ERR=110,END=120) label
      READ(inch,911,ERR=110,END=120) vegtyp
      WRITE(lchanl,920) label,' vegtyp.'

      READ(inch,900,ERR=110,END=120) label
      READ(inch,910,ERR=110,END=120) lai
      WRITE(lchanl,920) label,' lai.'

      READ(inch,900,ERR=110,END=120) label
      READ(inch,910,ERR=110,END=120) roufns
      WRITE(lchanl,920) label,' roufns.'

      READ(inch,900,ERR=110,END=120) label
      READ(inch,910,ERR=110,END=120) veg
      WRITE(lchanl,920) label,' veg.'

    END IF

  END IF

  IF( grdbas == 1 ) GO TO 930

  IF( varin == 1 ) THEN

    IF ( totin == 0 ) THEN
!
!-----------------------------------------------------------------------
!
!  Read in uprt, vprt, and wprt
!
!----------------------------------------------------------------------
!
      READ(inch,900,ERR=110,END=120) label
      READ(inch,910,ERR=110,END=120) uprt
      WRITE(lchanl,920) label,' uprt.'

      READ(inch,900,ERR=110,END=120) label
      READ(inch,910,ERR=110,END=120) vprt
      WRITE(lchanl,920) label,' vprt.'

      READ(inch,900,ERR=110,END=120) label
      READ(inch,910,ERR=110,END=120) wprt
      WRITE(lchanl,920) label,' wprt.'
!
!-----------------------------------------------------------------------
!
!  Read in scalars
!
!----------------------------------------------------------------------
!
      READ(inch,900,ERR=110,END=120) label
      READ(inch,910,ERR=110,END=120) ptprt
      WRITE(lchanl,920) label,' ptprt.'

      READ(inch,900,ERR=110,END=120) label
      READ(inch,910,ERR=110,END=120) pprt
      WRITE(lchanl,920) label,' pprt.'

    ELSE
!
!-----------------------------------------------------------------------
!
!  Read in total u, v, and w, and then derive uprt, vprt, and wprt
!
!-----------------------------------------------------------------------
!
      READ(inch,900,ERR=110,END=120) label
      READ(inch,910,ERR=110,END=120) uprt
      WRITE(lchanl,920) label,' u.'
      DO k=1,nz-1
        DO j=1,ny-1
          DO i=1,nx
            uprt(i,j,k) = uprt(i,j,k) - ubar(i,j,k)
          END DO
        END DO
      END DO

      READ(inch,900,ERR=110,END=120) label
      READ(inch,910,ERR=110,END=120) vprt
      WRITE(lchanl,920) label,' v.'
      DO k=1,nz-1
        DO j=1,ny
          DO i=1,nx-1
            vprt(i,j,k) = vprt(i,j,k) - vbar(i,j,k)
          END DO
        END DO
      END DO

      READ(inch,900,ERR=110,END=120) label
      READ(inch,910,ERR=110,END=120) wprt
      WRITE(lchanl,920) label,' w.'

      READ(inch,900,ERR=110,END=120) label
      READ(inch,910,ERR=110,END=120) ptprt
      WRITE(lchanl,920) label,' pt.'
      DO k=1,nz-1
        DO j=1,ny-1
          DO i=1,nx-1
            ptprt(i,j,k) = ptprt(i,j,k) - ptbar(i,j,k)
          END DO
        END DO
      END DO

      READ(inch,900,ERR=110,END=120) label
      READ(inch,910,ERR=110,END=120) pprt
      WRITE(lchanl,920) label,' p.'
      DO k=1,nz-1
        DO j=1,ny-1
          DO i=1,nx-1
            pprt(i,j,k) = pprt(i,j,k) - pbar(i,j,k)
          END DO
        END DO
      END DO

    END IF

  END IF

  IF( mstin == 1 ) THEN

    IF ( totin == 0 ) THEN

      READ(inch,900,ERR=110,END=120) label
      READ(inch,910,ERR=110,END=120) qvprt
      WRITE(lchanl,920) label,' qvprt.'

    ELSE

      READ(inch,900,ERR=110,END=120) label
      READ(inch,910,ERR=110,END=120) qvprt
      WRITE(lchanl,920) label,' qv.'
      DO k=1,nz-1
        DO j=1,ny-1
          DO i=1,nx-1
            qvprt(i,j,k) = qvprt(i,j,k) - qvbar(i,j,k)
          END DO
        END DO
      END DO

    END IF

    READ(inch,900,ERR=110,END=120) label
    READ(inch,910,ERR=110,END=120) qc
    WRITE(lchanl,920) label,' qc.'

    READ(inch,900,ERR=110,END=120) label
    READ(inch,910,ERR=110,END=120) qr
    WRITE(lchanl,920) label,' qr.'

    IF( rainin == 1 ) THEN

      READ(inch,900,ERR=110,END=120) label
      READ(inch,910,ERR=110,END=120) raing
      WRITE(lchanl,920) label,' raing.'

      READ(inch,900,ERR=110,END=120) label
      READ(inch,910,ERR=110,END=120) rainc
      WRITE(lchanl,920) label,' rainc.'

    END IF

    IF( prcin == 1 ) THEN

      READ(inch,900,ERR=110,END=120) label
      READ(inch,910,ERR=110,END=120)                                    &
          ((prcrate(i,j,1),i=1,nx),j=1,ny)
      WRITE(lchanl,920) label,' prcrate1.'

      READ(inch,900,ERR=110,END=120) label
      READ(inch,910,ERR=110,END=120)                                    &
          ((prcrate(i,j,2),i=1,nx),j=1,ny)
      WRITE(lchanl,920) label,' prcrate2.'

      READ(inch,900,ERR=110,END=120) label
      READ(inch,910,ERR=110,END=120)                                    &
          ((prcrate(i,j,3),i=1,nx),j=1,ny)
      WRITE(lchanl,920) label,' prcrate3.'

      READ(inch,900,ERR=110,END=120) label
      READ(inch,910,ERR=110,END=120)                                    &
          ((prcrate(i,j,4),i=1,nx),j=1,ny)
      WRITE(lchanl,920) label,' prcrate4.'

    END IF

    IF( icein == 1 ) THEN

      READ(inch,900,ERR=110,END=120) label
      READ(inch,910,ERR=110,END=120) qi
      WRITE(lchanl,920) label,' qi.'

      READ(inch,900,ERR=110,END=120) label
      READ(inch,910,ERR=110,END=120) qs
      WRITE(lchanl,920) label,' qs.'

      READ(inch,900,ERR=110,END=120) label
      READ(inch,910,ERR=110,END=120) qh
      WRITE(lchanl,920) label,' qh.'

    END IF

  END IF

  IF( tkein == 1 ) THEN

    READ(inch,900,ERR=110,END=120) label
    READ(inch,910,ERR=110,END=120) tke
    WRITE(lchanl,920) label,' tke.'

  END IF

  IF( trbin == 1 ) THEN

    READ(inch,900,ERR=110,END=120) label
    READ(inch,910,ERR=110,END=120) kmh
    WRITE(lchanl,920) label,' kmh.'

    READ(inch,900,ERR=110,END=120) label
    READ(inch,910,ERR=110,END=120) kmv
    WRITE(lchanl,920) label,' kmv.'

  END IF

  IF( sfcin == 1) THEN

    IF (nstyp <= 1) THEN

      READ(inch,900,ERR=110,END=120) label
      READ(inch,910,ERR=110,END=120)                                    &
            (((tsoil(i,j,k,0),i=1,nx),j=1,ny),k=1,nzsoil) 
      WRITE(lchanl,920) label,' tsoil.'

      READ(inch,900,ERR=110,END=120) label
      READ(inch,910,ERR=110,END=120)                                    &
            (((qsoil(i,j,k,0),i=1,nx),j=1,ny),k=1,nzsoil) 
      WRITE(lchanl,920) label,' qsoil.'

      READ(inch,900,ERR=110,END=120) label
      READ(inch,910,ERR=110,END=120)                                    &
            ((wetcanp(i,j,0),i=1,nx),j=1,ny)
      WRITE(lchanl,920) label,' wetcanp.'

    ELSE

      DO is=0,nstyp
        READ(inch,900,ERR=110,END=120) label
        READ(inch,910,ERR=110,END=120)                                  &
              (((tsoil(i,j,k,is),i=1,nx),j=1,ny),k=1,nzsoil) 
        WRITE(lchanl,920) label,' tsoil.'

        READ(inch,900,ERR=110,END=120) label
        READ(inch,910,ERR=110,END=120)                                  &
              (((qsoil(i,j,k,is),i=1,nx),j=1,ny),k=1,nzsoil) 
        WRITE(lchanl,920) label,' qsoil.'

        READ(inch,900,ERR=110,END=120) label
        READ(inch,910,ERR=110,END=120)                                  &
              ((wetcanp(i,j,is),i=1,nx),j=1,ny)
        WRITE(lchanl,920) label,' wetcanp.'
      END DO

    END IF

    IF(snowcin == 1) THEN
      READ(inch,900,ERR=110,END=120) label
      READ(inch,910,ERR=110,END=120)
      WRITE(lchanl,920) label,' snowcvr -- discarding.'
    END IF

    IF(snowin == 1) THEN
      READ(inch,900,ERR=110,END=120) label
      READ(inch,910,ERR=110,END=120)                                    &
          ((snowdpth(i,j),i=1,nx),j=1,ny)
      WRITE(lchanl,920) label,' snowdpth.'
    END IF

  END IF

  IF( radin == 1 ) THEN

    READ(inch,900,ERR=110,END=120) label
    READ(inch,910,ERR=110,END=120) radfrc
    WRITE(lchanl,920) label,' radfrc.'

    READ(inch,900,ERR=110,END=120) label
    READ(inch,910,ERR=110,END=120) radsw
    WRITE(lchanl,920) label,' radsw.'

    READ(inch,900,ERR=110,END=120) label
    READ(inch,910,ERR=110,END=120) rnflx
    WRITE(lchanl,920) label,' rnflx.'

    READ(inch,900,ERR=110,END=120) label
    READ(inch,910,ERR=110,END=120) radswnet
    WRITE(lchanl,920) label,' radswnet.'

    READ(inch,900,ERR=110,END=120) label
    READ(inch,910,ERR=110,END=120) radlwin
    WRITE(lchanl,920) label,' radlwin.'

  END IF

  IF( flxin == 1 ) THEN

    READ(inch,900,ERR=110,END=120) label
    READ(inch,910,ERR=110,END=120) usflx
    WRITE(lchanl,920) label,' usflx.'

    READ(inch,900,ERR=110,END=120) label
    READ(inch,910,ERR=110,END=120) vsflx
    WRITE(lchanl,920) label,' vsflx.'

    READ(inch,900,ERR=110,END=120) label
    READ(inch,910,ERR=110,END=120) ptsflx
    WRITE(lchanl,920) label,' ptsflx.'

    READ(inch,900,ERR=110,END=120) label
    READ(inch,910,ERR=110,END=120) qvsflx
    WRITE(lchanl,920) label,' qvsflx.'

  END IF
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

  RETURN

  900   FORMAT(1X,a)
  910   FORMAT(1X,8E16.9)
  911   FORMAT(1X,10I8)
  920   FORMAT(1X,'Field ',a12,' was read into array ',a)
!
!-----------------------------------------------------------------------
!
!  Error during read
!
!----------------------------------------------------------------------
!

  110   CONTINUE
  WRITE(6,'(/a/)') ' Error reading data in ASCREAD'
  ireturn=1
  RETURN
!
!-----------------------------------------------------------------------
!
!  End-of-file during read
!
!----------------------------------------------------------------------
!

  120   CONTINUE
  WRITE(6,'(/a/)') ' End of file reached in ASCREAD'
  ireturn=2
  RETURN
END SUBROUTINE ascread
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE ASCDUMP                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE ascdump(nx,ny,nz,nzsoil,nstyps,nchanl, grdbas,               &
           u,v,w,ptprt,pprt,qv,qc,qr,qi,qs,qh,tke,kmh,kmv,              &
           ubar,vbar,ptbar,pbar,rhobar,qvbar,                           &
           x,y,z,zp,zpsoil,                                             &
           soiltyp,stypfrct,vegtyp,lai,roufns,veg,                      &
           tsoil,qsoil,wetcanp,snowdpth,                                &
           raing,rainc,prcrate,                                         &
           radfrc,radsw,rnflx,radswnet,radlwin,                         &
           usflx,vsflx,ptsflx,qvsflx,                                   &
           tem1)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Write history data into channel nchanl in ASCII format.
!
!  All data read in are located at the original staggered grid points.
!
!  Note: coordinate fields are dumped as 3 dimensional fields which
!  have been converted from meters to kilometers.  This is for the
!  convenience of the plotting applications.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  5/31/92.
!
!  MODIFICATION HISTORY:
!
!  7/13/92 (K. Brewster)
!  Added runname, comment and version number writing
!
!  8/23/92 (M. Xue)
!  Modify to perform the dumping of both base and t-dependent arrays
!  and added control on grid staggering.
!
!  4/4/93  (M. Xue)
!  Modified, so that data on the original staggered grid are written
!  out. Averaging to the volume center is no longer done.
!
!  9/1/94 (Y. Lu)
!  Cleaned up documentation.
!
!  02/06/95 (Y. Liu)
!  Added map projection parameters into the ASCII dumping
!
!  05/31/95 (Y. Liu)
!  Changed the integer dumping format from 20I4 to 10I8.
!
!  12/09/1998 (Donghai Wang)
!  Added the snow cover.
!
!  05/23/2002 (J. Brotzge)
!  Added radiation, soil variables to allow for multiple soil schemes.  
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
!    grdbas   Flag indicating if this is a call for the data dump
!             of grid and base state arrays only. If so, grdbas=1.
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
!    vegtyp   Vegetation type
!    lai      Leaf Area Index
!    roufns   Surface roughness
!    veg      Vegetation fraction
!
!    tsoil    Soil temperature (K) (in deep 1 m layer)
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
!  OUTPUT:
!
!    None.
!
!  WORK ARRAY:
!
!    tem1     Temporary work array.
!
!-----------------------------------------------------------------------
!
!  The following parameters are passed into this subroutine through
!  a common block in globcst.inc, and they determine which
!  variables are output.
!
!  grdout =0 or 1. If grdout=0, grid variables are not dumped.
!  basout =0 or 1. If basout=0, base state variables are not dumped.
!  varout =0 or 1. If varout=0, model perturbation variables are not dumped.
!  mstout =0 or 1. If mstout=0, water variables are not dumped.
!  rainout=0 or 1. If rainout=0, rain variables are not dumped.
!  prcout =0 or 1. If prcout=0, precipitation rates are not dumped.
!  iceout =0 or 1. If iceout=0, qi, qs and qh are not dumped.
!  trbout =0 or 1. If trbout=0, turbulence parameter km is not dumped.
!  tkeout =0 or 1. If tkeout=0, tke is not dumped.
!  radout =0 or 1. If radout=0, radiation arrays are not dumped.
!  flxout =0 or 1. If flxout=0, surface fluxes are not dumped.
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

  INTEGER :: nchanl            ! FORTRAN I/O channel number for output
  INTEGER :: grdbas            ! If this is a grid/base state array dump

  REAL :: u     (nx,ny,nz)     ! Total u-velocity (m/s)
  REAL :: v     (nx,ny,nz)     ! Total v-velocity (m/s)
  REAL :: w     (nx,ny,nz)     ! Total w-velocity (m/s)
  REAL :: ptprt (nx,ny,nz)     ! Perturbation potential temperature (K)
  REAL :: pprt  (nx,ny,nz)     ! Perturbation pressure (Pascal)

  REAL :: qv    (nx,ny,nz)     ! Water vapor specific humidity (kg/kg)
  REAL :: qc    (nx,ny,nz)     ! Cloud water mixing ratio (kg/kg)
  REAL :: qr    (nx,ny,nz)     ! Rain water mixing ratio (kg/kg)
  REAL :: qi    (nx,ny,nz)     ! Cloud ice mixing ratio (kg/kg)
  REAL :: qs    (nx,ny,nz)     ! Snow mixing ratio (kg/kg)
  REAL :: qh    (nx,ny,nz)     ! Hail mixing ratio (kg/kg)
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
  REAL :: stypfrct(nx,ny,nstyps)     ! Soil type
  INTEGER :: vegtyp (nx,ny)          ! Vegetation type
  REAL :: lai    (nx,ny)             ! Leaf Area Index
  REAL :: roufns (nx,ny)             ! Surface roughness
  REAL :: veg    (nx,ny)             ! Vegetation fraction

  REAL :: tsoil (nx,ny,nzsoil,0:nstyps) ! Soil temperature (K)
  REAL :: qsoil (nx,ny,nzsoil,0:nstyps) ! Soil moisture (m**3/m**3) 
  REAL :: wetcanp(nx,ny,0:nstyps)       ! Canopy water amount
  REAL :: snowdpth(nx,ny)               ! Snow depth (m)


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

  REAL :: tem1  (nx,ny,nz)     ! Temporary work array
!
!-----------------------------------------------------------------------
!
!  Parameters describing this routine
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=40) :: fmtver
  PARAMETER (fmtver='004.10 ASCII Formatted Data')
  CHARACTER (LEN=10) :: tmunit
  PARAMETER (tmunit='seconds   ')
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k,l, idummy,is
  REAL :: rdummy
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


  WRITE(6,'(1x,a,f13.3/)') 'Writing history data at time=', curtim

  WRITE(nchanl,'(1x,a)') fmtver
  WRITE(nchanl,901) runname

  WRITE(nchanl,'(1x,i4)') nocmnt
  IF(nocmnt > 0) THEN
    DO l=1,nocmnt
      WRITE(nchanl,901) cmnt(l)
    END DO
  END IF

  WRITE(nchanl,'(1x,e16.8,1x,a10)') curtim,tmunit

  WRITE(nchanl,'(1x,4i12 )') nx,ny,nz,nzsoil 
!
!-----------------------------------------------------------------------
!
!  Write the flags for different data groups.
!
!-----------------------------------------------------------------------
!
  idummy = 0

  IF( grdbas == 1 ) THEN

    WRITE(nchanl,'(1x,10i8)')                                           &
                       1,      1,      0, mstout,      0,               &
                       0,      0,      0, landout, totout,              &
                  idummy, idummy, idummy, mapproj, month,               &
                     day,   year,   hour, minute, second

  ELSE

    WRITE(nchanl,'(1x,10i8)')                                           &
                  grdout, basout, varout, mstout, iceout,               &
                  trbout, sfcout, rainout,landout,totout,               &
                  tkeout, idummy, idummy, mapproj, month,               &
                     day,   year,   hour, minute, second

  END IF

  rdummy = 0.0
  WRITE(nchanl,900)                                                     &
                  umove,   vmove, xgrdorg, ygrdorg, trulat1,            &
                trulat2,  trulon,  sclfct,  rdummy,  rdummy,            &
                 rdummy,  rdummy,  rdummy,  rdummy,  rdummy,            &
                  tstop, thisdmp,  latitud, ctrlat,  ctrlon

  IF ( totout == 1 ) THEN
    WRITE(nchanl,'(1x,10i8)')                                           &
               nstyp,  prcout, radout, flxout,      0,  & ! 0 for snowcvr
           snowout,idummy, idummy, idummy, idummy,                      &
               idummy, idummy, idummy, idummy, idummy,                  &
               idummy, idummy, idummy, idummy, idummy

    WRITE(nchanl,900)                                                   &
               rdummy, rdummy, rdummy, rdummy, rdummy,                  &
               rdummy, rdummy, rdummy, rdummy, rdummy,                  &
               rdummy, rdummy, rdummy, rdummy, rdummy,                  &
               rdummy, rdummy, rdummy, rdummy, rdummy

  END IF

  IF( grdout == 1  .OR. grdbas == 1 ) THEN

    WRITE(nchanl,'(1x,a)') 'x coordinate'
    WRITE(nchanl,900) x

    WRITE(nchanl,'(1x,a)') 'y coordinate'
    WRITE(nchanl,900) y

    WRITE(nchanl,'(1x,a)') 'z coordinate'
    WRITE(nchanl,900) z

    CALL edgfill(zp,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz)
    WRITE(nchanl,'(1x,a)') 'zp coordinat'
    WRITE(nchanl,900) zp

    CALL edgfill(zp,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nzsoil,1,nzsoil)
    WRITE(nchanl,'(1x,a)') 'zpsoil coordinat'
    WRITE(nchanl,900) zpsoil 

  END IF
!
!-----------------------------------------------------------------------
!
!  If basout=1, write out base state variables.
!
!-----------------------------------------------------------------------
!
  IF( basout == 1  .OR. grdbas == 1 ) THEN

    CALL edgfill(ubar,1,nx,1,nx, 1,ny,1,ny-1, 1,nz,1,nz-1)
    WRITE(nchanl,'(1x,a)') 'ubar        '
    WRITE(nchanl,900) ubar

    CALL edgfill(vbar,1,nx,1,nx-1, 1,ny,1,ny, 1,nz,1,nz-1)
    WRITE(nchanl,'(1x,a)') 'vbar        '
    WRITE(nchanl,900) vbar

    DO k=1,nz
      DO j=1,ny
        DO i=1,nx
          tem1(i,j,k) = 0.0
        END DO
      END DO
    END DO
    WRITE(nchanl,'(1x,a)') 'wbar        '
    WRITE(nchanl,900) tem1

    CALL edgfill(ptbar,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)
    WRITE(nchanl,'(1x,a)') 'ptbar       '
    WRITE(nchanl,900) ptbar

    CALL edgfill(pbar,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)
    WRITE(nchanl,'(1x,a)') 'pbar        '
    WRITE(nchanl,900) pbar

    IF(mstout == 1) THEN

      CALL edgfill(qvbar,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)
      WRITE(nchanl,'(1x,a)') 'qvbar       '
      WRITE(nchanl,900) qvbar

    END IF

    IF(landout == 1) THEN

      IF( nstyp <= 1 ) THEN
        CALL iedgfill(soiltyp(1,1,1),1,nx,1,nx-1, 1,ny,1,ny-1,          &
                      1,1,1,1)
        WRITE(nchanl,'(1x,a)') 'soiltyp     '
        WRITE(nchanl,902) ((soiltyp(i,j,1),i=1,nx),j=1,ny)
      ELSE
        DO is = 1,nstyp
          CALL iedgfill(soiltyp(1,1,is),1,nx,1,nx-1, 1,ny,1,ny-1,       &
                        1,1,1,1)
          WRITE(nchanl,'(1x,a)') 'soiltyp     '
          WRITE(nchanl,902) ((soiltyp(i,j,is),i=1,nx),j=1,ny)

          CALL edgfill(stypfrct(1,1,is),1,nx,1,nx-1, 1,ny,1,ny-1,       &
                       1,1,1,1)
          WRITE(nchanl,'(1x,a)') 'stypfrct     '
          WRITE(nchanl,900) ((stypfrct(i,j,is),i=1,nx),j=1,ny)
        END DO

      END IF

      CALL iedgfill(vegtyp ,1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)
      WRITE(nchanl,'(1x,a)') 'vegtyp      '
      WRITE(nchanl,902) vegtyp

      CALL edgfill(lai    ,1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)
      WRITE(nchanl,'(1x,a)') 'lai         '
      WRITE(nchanl,900) lai

      CALL edgfill(roufns ,1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)
      WRITE(nchanl,'(1x,a)') 'roufns      '
      WRITE(nchanl,900) roufns

      CALL edgfill(veg    ,1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)
      WRITE(nchanl,'(1x,a)') 'veg         '
      WRITE(nchanl,900) veg

    END IF

  END IF

  IF ( grdbas == 1 ) RETURN
!
!-----------------------------------------------------------------------
!
!  If varout = 1, Write out uprt, vprt, wprt, ptprt, pprt.
!
!-----------------------------------------------------------------------
!
  IF( varout == 1 ) THEN
!
!-----------------------------------------------------------------------
!
!  Write out u, v and w
!
!-----------------------------------------------------------------------
!
    IF ( totout == 0 ) THEN
!
!-----------------------------------------------------------------------
!
!  Write out perturbations to history dump
!
!-----------------------------------------------------------------------
!
      DO k=1,nz-1
        DO j=1,ny-1
          DO i=1,nx
            tem1(i,j,k)=u(i,j,k)-ubar(i,j,k)
          END DO
        END DO
      END DO
      CALL edgfill(tem1,1,nx,1,nx, 1,ny,1,ny-1, 1,nz,1,nz-1)
      WRITE(nchanl,'(1x,a)') 'uprt        '
      WRITE(nchanl,900) tem1

      DO k=1,nz-1
        DO i=1,nx-1
          DO j=1,ny
            tem1(i,j,k)=v(i,j,k)-vbar(i,j,k)
          END DO
        END DO
      END DO
      CALL edgfill(tem1,1,nx,1,nx-1, 1,ny,1,ny, 1,nz,1,nz-1)
      WRITE(nchanl,'(1x,a)') 'vprt        '
      WRITE(nchanl,900) tem1

      CALL edgfill(w,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz)
      WRITE(nchanl,'(1x,a)') 'wprt        '
      WRITE(nchanl,900) w
!
!-----------------------------------------------------------------------
!
!    Write out scalars
!
!-----------------------------------------------------------------------
!

      CALL edgfill(ptprt,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)
      WRITE(nchanl,'(1x,a)') 'ptprt       '
      WRITE(nchanl,900) ptprt

      CALL edgfill(pprt,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)
      WRITE(nchanl,'(1x,a)') 'pprt        '
      WRITE(nchanl,900) pprt

    ELSE
!
!-----------------------------------------------------------------------
!
!  Write out total values to history dump
!
!-----------------------------------------------------------------------
!
      CALL edgfill(u,1,nx,1,nx, 1,ny,1,ny-1, 1,nz,1,nz-1)
      WRITE(nchanl,'(1x,a)') 'u           '
      WRITE(nchanl,900) u

      CALL edgfill(v,1,nx,1,nx-1, 1,ny,1,ny, 1,nz,1,nz-1)
      WRITE(nchanl,'(1x,a)') 'v           '
      WRITE(nchanl,900) v

      CALL edgfill(w,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz)
      WRITE(nchanl,'(1x,a)') 'w           '
      WRITE(nchanl,900) w
!
!-----------------------------------------------------------------------
!
!  Write out scalars
!
!-----------------------------------------------------------------------
!
      DO k=1,nz-1
        DO j=1,ny-1
          DO i=1,nx-1
            tem1(i,j,k) = ptbar(i,j,k) + ptprt(i,j,k)
          END DO
        END DO
      END DO
      CALL edgfill(tem1,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)
      WRITE(nchanl,'(1x,a)') 'pt          '
      WRITE(nchanl,900) tem1

      DO k=1,nz-1
        DO j=1,ny-1
          DO i=1,nx-1
            tem1(i,j,k) = pbar(i,j,k) + pprt(i,j,k)
          END DO
        END DO
      END DO
      CALL edgfill(tem1,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)
      WRITE(nchanl,'(1x,a)') 'p           '
      WRITE(nchanl,900) tem1

    END IF

  END IF
!
!-----------------------------------------------------------------------
!
!  If mstout = 1, Write out moisture variables
!
!-----------------------------------------------------------------------
!
  IF( mstout == 1 ) THEN

    IF ( totout == 0 ) THEN
!
!-----------------------------------------------------------------------
!
!  Write out perturbations to history dump
!
!-----------------------------------------------------------------------
!
      DO k=1,nz-1
        DO j=1,ny-1
          DO i=1,nx-1
            tem1(i,j,k) = qv(i,j,k)-qvbar(i,j,k)
          END DO
        END DO
      END DO
      CALL edgfill(tem1,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)
      WRITE(nchanl,'(1x,a)') 'qvprt       '
      WRITE(nchanl,900) tem1

    ELSE
!
!-----------------------------------------------------------------------
!
!  Write out total values to history dump
!
!-----------------------------------------------------------------------
!
      CALL edgfill(qv,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)
      WRITE(nchanl,'(1x,a)') 'qv          '
      WRITE(nchanl,900) qv

    END IF

    CALL edgfill(qc,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)
    WRITE(nchanl,'(1x,a)') 'qc          '
    WRITE(nchanl,900) qc

    CALL edgfill(qr,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)
    WRITE(nchanl,'(1x,a)') 'qr          '
    WRITE(nchanl,900) qr

    IF( rainout == 1 ) THEN

      CALL edgfill(raing,1,nx,1,nx-1, 1,ny,1,ny-1,1,1,1,1)
      WRITE(nchanl,'(1x,a)') 'raing       '
      WRITE(nchanl,900) raing

      CALL edgfill(rainc,1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)
      WRITE(nchanl,'(1x,a)') 'rainc       '
      WRITE(nchanl,900) rainc

    END IF    ! rainout

    IF( prcout == 1 ) THEN

      CALL edgfill(prcrate,1,nx,1,nx-1, 1,ny,1,ny-1,1,4,1,4)

      WRITE(nchanl,'(1x,a)') 'prcrate1    '
      WRITE(nchanl,900) ((prcrate(i,j,1),i=1,nx),j=1,ny)
      WRITE(nchanl,'(1x,a)') 'prcrate2    '
      WRITE(nchanl,900) ((prcrate(i,j,2),i=1,nx),j=1,ny)
      WRITE(nchanl,'(1x,a)') 'prcrate3    '
      WRITE(nchanl,900) ((prcrate(i,j,3),i=1,nx),j=1,ny)
      WRITE(nchanl,'(1x,a)') 'prcrate4    '
      WRITE(nchanl,900) ((prcrate(i,j,4),i=1,nx),j=1,ny)

    END IF    ! prcout

    IF( iceout == 1 ) THEN

      CALL edgfill(qi,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)
      WRITE(nchanl,'(1x,a)') 'qi          '
      WRITE(nchanl,900) qi

      CALL edgfill(qs,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)
      WRITE(nchanl,'(1x,a)') 'qs          '
      WRITE(nchanl,900) qs

      CALL edgfill(qh,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)
      WRITE(nchanl,'(1x,a)') 'qh          '
      WRITE(nchanl,900) qh

    END IF    ! iceout

  END IF  ! mstout

  IF( tkeout == 1 ) THEN

    CALL edgfill(tke,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)
    WRITE(nchanl,'(1x,a)') 'tke         '
    WRITE(nchanl,900) tke

  END IF   ! tkeout


  IF( trbout == 1 ) THEN

    CALL edgfill(kmh,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)
    WRITE(nchanl,'(1x,a)') 'kmh         '
    WRITE(nchanl,900) kmh

    CALL edgfill(kmv,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)
    WRITE(nchanl,'(1x,a)') 'kmv         '
    WRITE(nchanl,900) kmv

  END IF   ! trbout

!
!-----------------------------------------------------------------------
!
!  If sfcout = 1, write out the surface variables, tsfc, tsoil,
!  wetsfc, wetdp, and wetcanp.
!
!-----------------------------------------------------------------------
!
  IF ( sfcout == 1) THEN

    IF ( nstyp <= 1 ) THEN

      CALL edgfill(tsoil(1,1,1,0),1,nx,1,nx-1,1,ny,1,ny-1,1,nzsoil,1,  &
           nzsoil)
      WRITE(nchanl,'(1x,a)') 'tsoil         '
      WRITE(nchanl,900) (((tsoil(i,j,k,0),i=1,nx),j=1,ny),k=1,nzsoil) 

      CALL edgfill(qsoil(1,1,1,0),1,nx,1,nx-1,1,ny,1,ny-1,1,nzsoil,1, &
           nzsoil)
      WRITE(nchanl,'(1x,a)') 'qsoil         '
      WRITE(nchanl,900) (((qsoil(i,j,k,0),i=1,nx),j=1,ny),k=1,nzsoil) 

      CALL edgfill(wetcanp(1,1,0),1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)
      WRITE(nchanl,'(1x,a)') 'wetcanp      '
      WRITE(nchanl,900) ((wetcanp(i,j,0),i=1,nx),j=1,ny)

    ELSE

      DO is=0,nstyp
        CALL edgfill(tsoil(1,1,1,is),1,nx,1,nx-1, 1,ny,1,ny-1,             &
                     1,nzsoil,1,nzsoil)
        WRITE(nchanl,'(1x,a,I2.2)') 'tsoil       ',is
        WRITE(nchanl,900) (((tsoil(i,j,k,is),i=1,nx),j=1,ny),k=1,nzsoil)

        CALL edgfill(qsoil(1,1,1,is),1,nx,1,nx-1, 1,ny,1,ny-1,            &
                     1,nzsoil,1,nzsoil)
        WRITE(nchanl,'(1x,a,I2.2)') 'qsoil       ',is
        WRITE(nchanl,900) (((qsoil(i,j,k,is),i=1,nx),j=1,ny),k=1,nzsoil) 

        CALL edgfill(wetcanp(1,1,is),1,nx,1,nx-1, 1,ny,1,ny-1,          &
                     1,1,1,1)
        WRITE(nchanl,'(1x,a,I2.2)') 'wetcanp    ',is
        WRITE(nchanl,900) ((wetcanp(i,j,is),i=1,nx),j=1,ny)
      END DO

    END IF

    IF (snowout == 1) THEN

      CALL edgfill(snowdpth,1,nx,1,nx-1, 1,ny,1,ny-1,                   &
                    1,1,1,1)
      WRITE(nchanl,'(1x,a)') 'snowdpth     '
      WRITE(nchanl,900) ((snowdpth(i,j),i=1,nx),j=1,ny)

    END IF

  END IF   ! sfcout done
!
!-----------------------------------------------------------------------
!
!  Write out radiation ararys to history dump
!
!-----------------------------------------------------------------------
!
  IF( radout == 1 ) THEN

    CALL edgfill(radfrc,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)
    WRITE(nchanl,'(1x,a)') 'radfrc      '
    WRITE(nchanl,900) radfrc

    CALL edgfill(radsw,1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)
    WRITE(nchanl,'(1x,a)') 'radsw       '
    WRITE(nchanl,900) radsw

    CALL edgfill(rnflx,1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)
    WRITE(nchanl,'(1x,a)') 'rnflx       '
    WRITE(nchanl,900) rnflx

    CALL edgfill(radswnet,1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)
    WRITE(nchanl,'(1x,a)') 'radswnet       '
    WRITE(nchanl,900) radswnet 

    CALL edgfill(radlwin,1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)
    WRITE(nchanl,'(1x,a)') 'radlwin       '
    WRITE(nchanl,900) radlwin


  END IF   ! radout

!
!-----------------------------------------------------------------------
!
!  Write out surface fluxes to history dump
!
!-----------------------------------------------------------------------
!
  IF( flxout == 1 ) THEN

    CALL edgfill(usflx,1,nx,1,nx, 1,ny,1,ny-1, 1,1,1,1)
    WRITE(nchanl,'(1x,a)') 'usflx       '
    WRITE(nchanl,900) usflx

    CALL edgfill(vsflx,1,nx,1,nx-1, 1,ny,1,ny, 1,1,1,1)
    WRITE(nchanl,'(1x,a)') 'vsflx       '
    WRITE(nchanl,900) vsflx

    CALL edgfill(ptsflx,1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)
    WRITE(nchanl,'(1x,a)') 'ptsflx      '
    WRITE(nchanl,900) ptsflx

    CALL edgfill(qvsflx,1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)
    WRITE(nchanl,'(1x,a)') 'qvsflx      '
    WRITE(nchanl,900) qvsflx

  END IF   ! flxout

  900   FORMAT(1X,8E16.9)
  901   FORMAT(1X,a80)
  902   FORMAT(1X,10I8)

  RETURN
END SUBROUTINE ascdump
