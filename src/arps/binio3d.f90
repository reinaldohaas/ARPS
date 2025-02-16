!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE BINREAD                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE binread(nx,ny,nz,nzsoil,nstyps,grdbas,inch,time,x,y,z,zp,    &
           zpsoil,uprt, vprt, wprt, ptprt, pprt,                        &
           qvprt, qscalar, tke,kmh,kmv,                                 &
           ubar, vbar, wbar, ptbar, pbar, rhobar, qvbar,                &
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
!  Read in binary data set created by ARPS using history dump format
!  No. 1.
!  All data read in are located at the original staggered grid points
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  2/27/92.
!
!  MODIFICATION HISTORY:
!
!  6/08/92
!  Added full documentation (K. Brewster)
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
!  Added map projection parameters into the binary dumping
!
!  03/26/96 (G. Bassett)
!  Backwards compatibility added for ARPS 3.2 and ARPS 4.0 binary
!  history dump formats.
!
!  12/09/1998 (Donghai Wang)
!  Added the snow cover.
!
!  05/15/2002 (J. Brotzge)
!  Added additional variables to allow for multiple soil schemes
!
!  1 June 2002 Eric Kemp
!  Bug fixes for new soil variables.
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
!    qscalar
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
!    radswnet Net shortwave radiation
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
  INCLUDE 'phycst.inc'
  INCLUDE 'mp.inc'            ! mpi parameters.
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  INTEGER :: nx,ny,nz          ! Number of grid points in 3 directions
  INTEGER :: nzsoil            ! Number of grid points in the soil

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
                               ! w-point of the soil (m)


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
  REAL :: qscalar(nx,ny,nz,nscalar)

  REAL :: tke  (nx,ny,nz)      ! Turbulent Kinetic Energy ((m/s)**2)
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

  INTEGER :: nstyps                    ! Number of soil type
  INTEGER :: soiltyp(nx,ny,nstyps)     ! Soil type
  REAL    :: stypfrct(nx,ny,nstyps)    ! Soil type fraction
  INTEGER :: vegtyp(nx,ny)             ! Vegetation type
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
  REAL :: ptsflx(nx,ny)        ! Surface heat flux (K*kg/(m**2*s))
  REAL :: qvsflx(nx,ny)        ! Surface moisture flux (kg/(m**2*s))

  INTEGER :: ireturn           ! Return status indicator
!
!-----------------------------------------------------------------------
!
!  Parameters describing routine that wrote the gridded data
!
!-----------------------------------------------------------------------
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
  INTEGER            :: intver320,intver400,intver410,intver500,intver530

  PARAMETER (fmtver320='003.20 Binary Data',intver320=320)
  PARAMETER (fmtver400='004.00 Binary Data',intver400=400)
  PARAMETER (fmtver410='004.10 Binary Data',intver410=410)
  PARAMETER (fmtver500='005.00 Binary Data',intver500=500)
  PARAMETER (fmtver530='005.30 Binary Data',intver530=530)

  CHARACTER (LEN=40) :: fmtverin
  INTEGER            :: intver

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
  INTEGER :: nstyp1
  CHARACTER (LEN=12) :: label
  INTEGER :: nxin,nyin,nzin,nzsoilin
  INTEGER :: bgrdin,bbasin,bvarin,bicein,btkein,btrbin
  INTEGER :: idummy

  INTEGER :: nq, nqin
  INTEGER :: nqscalarin(nscalar)

  REAL(4) :: rdummy
  REAL(4) :: umove4,   vmove4
  REAL(4) :: xgrdorg4, ygrdorg4
  REAL(4) :: trulat14, trulat24, trulon4
  REAL(4) :: sclfct4
  REAL(4) :: tstop4,   thisdmp4
  REAL(4) :: latitud4, ctrlat4,  ctrlon4
  REAL(4) :: n0rain4, n0snow4, n0hail4, rhosnow4, rhohail4
  REAL(4) :: ntcloud4, n0grpl4, rhoice4, rhogrpl4
  REAL(4) :: alpharain4,alphaice4,alphasnow4,alphagrpl4,alphahail4

  REAL(4), ALLOCATABLE :: invar1d(:)
  REAL(4), ALLOCATABLE :: invar2d(:,:)
  REAL(4), ALLOCATABLE :: invar3d(:,:,:)
  REAL(4), ALLOCATABLE :: invar3dsoil(:,:,:)
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  ALLOCATE(invar1d(MAX(nx,ny,nz)), STAT = idummy)
  CALL check_alloc_status(idummy, "BINREAD:invar1d")
  ALLOCATE(invar2d(nx,ny),         STAT = idummy)
  CALL check_alloc_status(idummy, "BINREAD:invar2d")
  ALLOCATE(invar3d(nx,ny,nz),      STAT = idummy)
  CALL check_alloc_status(idummy, "BINREAD:invar3d")
  ALLOCATE(invar3dsoil(nx,ny,nzsoil), STAT = idummy)
  CALL check_alloc_status(idummy, "BINREAD:invar3dsoil")
!
!-----------------------------------------------------------------------
!
!  Read header info
!
!-----------------------------------------------------------------------
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
    IF (myproc == 0) WRITE(6,'(/1x,a,a,a/)')                            &
                  'Incoming data format, fmtverin=',fmtverin,           &
                  ', not found. The Job stopped.'
    CALL arpsstop('arpstop called from binread. ',1)
  END IF

  IF (myproc == 0) WRITE(6,'(/1x,a,a/)')                                &
                         'Incoming data format, fmtverin=',fmtverin

  READ(inch,ERR=110,END=120) runname
  READ(inch,ERR=110,END=120) nocmnt
  IF( nocmnt > 0 ) THEN
    DO i=1,nocmnt
      READ(inch,ERR=110,END=120) cmnt(i)
    END DO
  END IF

  IF (myproc == 0)  &
     WRITE(6,'(//''  THE NAME OF THIS RUN IS:  '',A//)') runname

  IF( nocmnt > 0 ) THEN
    IF(myproc == 0)THEN
      DO i=1,nocmnt
        WRITE(6,'(1x,a)') cmnt(i)
      END DO
    END IF
  END IF

  READ(inch,ERR=110,END=120) rdummy,tmunit
  time = rdummy
!
!-----------------------------------------------------------------------
!
!  Get dimensions of data in binary file and check against
!  the dimensions passed to BINREAD
!
!-----------------------------------------------------------------------
!
  IF (intver <= intver410) THEN

    READ(inch,ERR=110,END=120) nxin, nyin, nzin
    nzsoilin = 2  ! for version prior to 410, it is a two-layer soil model

  ELSE IF (intver >= intver500) THEN

    READ(inch,ERR=110,END=120) nxin, nyin, nzin,nzsoilin

  END IF
!
! Data validation: dimensions
!
  IF( nxin /= nx .OR. nyin /= ny .OR. nzin /= nz) THEN
    IF(myproc == 0)THEN
      WRITE(6,'(1x,a)')                                                   &
         ' Dimensions in BINREAD inconsistent with data.'
      WRITE(6,'(1x,a,3I15)') ' Read were: ', nxin, nyin, nzin
      WRITE(6,'(1x,a,3I15)') ' Expected:  ', nx, ny, nz
      WRITE(6,'(1x,a)')      ' Program aborted in BINREAD.'
    END IF
    CALL arpsstop('arpstop called from binread nx-ny-nz read ',1)
  END IF

  IF(nzsoilin /= nzsoil) THEN

    IF (intver <= intver410) THEN
      IF(myproc == 0) WRITE(6,'(1x,a,a/,2(1x,a/))')           &
        ' The incoming data version is ', fmtverin,            &
        ' In the input file, nzsoil must be set to 2. ',       &
        ' Program aborted in BINREAD.'

    ELSE IF (intver >= intver500) THEN

      IF(myproc == 0)THEN
        WRITE(6,'(1x,a)') &
                ' Dimensions in BINREAD inconsistent with data.'
        WRITE(6,'(1x,a,I15)') ' Read were: ', nzsoilin
        WRITE(6,'(1x,a,I15)') ' Expected:  ', nzsoil
        WRITE(6,'(1x,a)') ' Program aborted in BINREAD.'
      END IF
    END IF
    CALL arpsstop('arpstop called from binread nzsoil read ',1)

  END IF
!
!-----------------------------------------------------------------------
!
!  Read in flags for different data groups
!
!-----------------------------------------------------------------------
!
  IF( grdbas == 1 ) THEN ! Read grid and base state arrays

    IF (myproc == 0)  &
       WRITE(lchanl,'(1x,a,f8.1,a,f8.3,a/)')                            &
             'To read grid and base state data at time ', time,         &
             ' secs = ',(time/60.),' mins.'

    READ(inch,ERR=110,END=120)                                          &
         bgrdin,bbasin,bvarin,mstin,bicein,                             &
         btrbin,idummy,idummy,landin,totin,                             &
         btkein,idummy,idummy,mapproj,month,                            &
         day, year, hour, minute, second

  ELSE                   ! time-dependent data reading

    IF (myproc == 0) WRITE(lchanl,'(1x,a,f8.1,a,f8.3,a/)')              &
           'To read data for time:',time,' secs = ',(time/60.),' mins.'

    READ(inch,ERR=110,END=120)                                          &
         grdin,basin,varin,mstin,icein,                                 &
         trbin, sfcin,rainin,landin,totin,                              &
         tkein,nscalarin,idummy,mapproj, month,                         &
         day, year, hour, minute, second

    IF (intver >= intver530) THEN   ! newer version 5.3 or later
      READ(inch,ERR=110,END=120)                                        &
         p_qcin,p_qrin,p_qiin,p_qsin,p_qgin,p_qhin,                     &
         p_ncin,p_nrin,p_niin,p_nsin,p_ngin,p_nhin,                     &
         p_zrin,p_ziin,p_zsin,p_zsin,p_zhin,idummy,                     &
         idummy,idummy,idummy,idummy,idummy,idummy,                     &
         idummy,idummy,idummy,idummy,p_ccin,idummy
    ELSE                      ! version earlier than 5.3
      p_qcin=0; p_qrin=0; p_qiin=0; p_qsin=0; p_qgin=0; p_qhin=0
      p_ncin=0; p_nrin=0; p_niin=0; p_nsin=0; p_nhin=0; p_ngin=0
                p_zrin=0; p_ziin=0; p_zsin=0; p_zhin=0; p_zgin=0
      p_ccin=0

      nscalarin = 0
      IF (mstin == 1) THEN
        p_qcin=1
        p_qrin=2
        nscalarin = nscalarin + 2
        IF (icein == 2) THEN
          p_qiin=3
          p_qsin=4
          p_qhin=5
          nscalarin = nscalarin + 3
        END IF
      END IF
    END IF

    nqscalarin(:) = 0
    IF (P_QC > 0) nqscalarin(P_QC) = p_qcin
    IF (P_QR > 0) nqscalarin(P_QR) = p_qrin
    IF (P_QI > 0) nqscalarin(P_QI) = p_qiin
    IF (P_QS > 0) nqscalarin(P_QS) = p_qsin
    IF (P_QG > 0) nqscalarin(P_QG) = p_qgin
    IF (P_QH > 0) nqscalarin(P_QH) = p_qhin

    IF (P_NC > 0) nqscalarin(P_NC) = p_ncin
    IF (P_NR > 0) nqscalarin(P_NR) = p_nrin
    IF (P_NI > 0) nqscalarin(P_NI) = p_niin
    IF (P_NS > 0) nqscalarin(P_NS) = p_nsin
    IF (P_NG > 0) nqscalarin(P_NG) = p_ngin
    IF (P_NH > 0) nqscalarin(P_NH) = p_nhin

    IF (P_ZR > 0) nqscalarin(P_ZR) = p_zrin
    IF (P_ZI > 0) nqscalarin(P_ZI) = p_ziin
    IF (P_ZS > 0) nqscalarin(P_ZS) = p_zsin
    IF (P_ZG > 0) nqscalarin(P_ZG) = p_zgin
    IF (P_ZH > 0) nqscalarin(P_ZH) = p_zhin

    IF (P_CC > 0) nqscalarin(P_CC) = p_ccin

  END IF

  IF (intver < intver530) THEN

    READ(inch,ERR=110,END=120)                                          &
                    umove4,  vmove4,  xgrdorg4,ygrdorg4,trulat14,       &
                    trulat24,trulon4, sclfct4, n0rain4, n0snow4,        &
                    n0hail4, rhosnow4,rhohail4,rdummy,  rdummy,         &
                    tstop4,  thisdmp4,latitud4,ctrlat4, ctrlon4
  ELSE
    READ(inch,ERR=110,END=120)                                          &
                    umove4,  vmove4,  xgrdorg4,ygrdorg4,trulat14,       &
                    trulat24,trulon4, sclfct4, ntcloud4, n0rain4,       &
                    n0snow4, n0grpl4, n0hail4, rhoice4, rhosnow4,       &
                    rhogrpl4,rhohail4,alpharain4,alphaice4,alphasnow4,  &
                    alphagrpl4,alphahail4,rdummy, rdummy,  rdummy,      &
                    tstop4,  thisdmp4,latitud4,ctrlat4, ctrlon4

   ntcloud = ntcloud4
   n0grpl  = n0grpl4
   rhoice  = rhoice4
   rhosnow = rhosnow4
   rhogrpl = rhogrpl4
   alpharain = alpharain4
   alphaice  = alphaice4
   alphasnow = alphasnow4
   alphagrpl = alphagrpl4
   alphahail = alphahail4
  END IF

  umove = umove4
  vmove = vmove4
  xgrdorg = xgrdorg4
  ygrdorg = ygrdorg4
  trulat1 = trulat14
  trulat2 = trulat24
  trulon  = trulon4
  sclfct  = sclfct4
  tstop   = tstop4
  thisdmp = thisdmp4
  latitud = latitud4
  ctrlat  = ctrlat4
  ctrlon  = ctrlon4
  n0rain  = n0rain4
  n0snow  = n0snow4
  n0hail  = n0hail4
  rhosnow = rhosnow4
  rhohail = rhohail4

  IF ( totin /= 0 ) THEN
!
!-----------------------------------------------------------------------
!
!  Read in additional parameters for ARPS history dump 4.0 or later
!  version.
!
!-----------------------------------------------------------------------
!
    READ(inch,ERR=110,END=120)                                          &
         nstyp1, prcin, radin, flxin,snowcin,                           &
         snowin,idummy,idummy,idummy,idummy,                            &
         idummy,idummy,idummy,idummy,idummy,                            &
         idummy,idummy,idummy,idummy,idummy

    IF ( nstyp1 < 1 ) THEN
      nstyp1 = 1
    END IF

    READ(inch,ERR=110,END=120)                                          &
         rdummy,rdummy,rdummy,rdummy,rdummy,                            &
         rdummy,rdummy,rdummy,rdummy,rdummy,                            &
         rdummy,rdummy,rdummy,rdummy,rdummy,                            &
         rdummy,rdummy,rdummy,rdummy,rdummy
  END IF
!
!-----------------------------------------------------------------------
!
!  Read in x, y, z and zp arrays.
!
!----------------------------------------------------------------------
!
  IF( grdin == 1 .OR. grdbas == 1 ) THEN
    READ(inch,ERR=110,END=120) label
    READ(inch,ERR=110,END=120) invar1d(1:nx)
    x(:) = invar1d(1:nx)
    IF (myproc == 0) WRITE(lchanl,910) label,' x.'

    READ(inch,ERR=110,END=120) label
    READ(inch,ERR=110,END=120) invar1d(1:ny)
    y(:) = invar1d(1:ny)
    IF (myproc == 0) WRITE(lchanl,910) label,' y.'

    READ(inch,ERR=110,END=120) label
    READ(inch,ERR=110,END=120) invar1d(1:nz)
    z(:) = invar1d(1:nz)
    IF (myproc == 0) WRITE(lchanl,910) label,' z.'

    READ(inch,ERR=110,END=120) label
    READ(inch,ERR=110,END=120) invar3d
    zp(:,:,:) = invar3d(:,:,:)
    IF (myproc == 0) WRITE(lchanl,910) label,' zp.'

    IF (intver <= intver410) THEN
!
! nzsoil must equal to 2, 06/28/2002, Zuwen
! assume zpsoil(,,2) is one meter below the surface.
!
      DO j=1,ny
        DO i=1,nx
          zpsoil(i,j,1)=zp(i,j,2)
          zpsoil(i,j,2)=zpsoil(i,j,1)-1.
        END DO
      END DO

      IF (myproc == 0) THEN
        WRITE(lchanl,910) ' Assign zpsoil. '
        WRITE(lchanl,910) ' Assume zpsoil(,,1) is zp(,,2). '
        WRITE(lchanl,910) ' Assume zpsoil(,,2) is zp(,,2)-1. '
      END IF

    ELSE IF (intver >= intver500) THEN

      READ(inch,ERR=110,END=120) label
      READ(inch,ERR=110,END=120) invar3dsoil
      zpsoil(:,:,:) = invar3dsoil(:,:,:)
      IF (myproc == 0)  WRITE(lchanl,910) label,' zpsoil.'

    END IF  ! intver

  END IF  ! grdin
!
!-----------------------------------------------------------------------
!
!  Read in base state fields
!
!----------------------------------------------------------------------
!
  IF( basin == 1 .OR. grdbas == 1 ) THEN

    READ(inch,ERR=110,END=120) label
    READ(inch,ERR=110,END=120) invar3d
    ubar(:,:,:) = invar3d(:,:,:)
    IF (myproc == 0) WRITE(lchanl,910) label,' ubar.'

    READ(inch,ERR=110,END=120) label
    READ(inch,ERR=110,END=120) invar3d
    vbar(:,:,:) = invar3d(:,:,:)
    IF (myproc == 0) WRITE(lchanl,910) label,' vbar.'

    READ(inch,ERR=110,END=120) label
    READ(inch,ERR=110,END=120) invar3d
    wbar(:,:,:) = invar3d(:,:,:)
    IF (myproc == 0) WRITE(lchanl,910) label,' wbar.'

    READ(inch,ERR=110,END=120) label
    READ(inch,ERR=110,END=120) invar3d
    ptbar(:,:,:) = invar3d(:,:,:)
    IF (myproc == 0)  WRITE(lchanl,910) label,' ptbar.'

    READ(inch,ERR=110,END=120) label
    READ(inch,ERR=110,END=120) invar3d
    pbar(:,:,:) = invar3d(:,:,:)
    IF (myproc == 0)  WRITE(lchanl,910) label,' pbar.'

    IF( mstin == 1) THEN
      READ(inch,ERR=110,END=120) label
      READ(inch,ERR=110,END=120) invar3d
      qvbar(:,:,:) = invar3d(:,:,:)
      IF (myproc == 0)  WRITE(lchanl,910) label,' qvbar.'
    END IF

    IF (landin == 1) THEN

      IF (nstyp1 <= 1) THEN

        READ(inch,ERR=110,END=120) label
        READ(inch,ERR=110,END=120) ((soiltyp(i,j,1),i=1,nx),j=1,ny)
        IF (myproc == 0) WRITE(lchanl,910) label,' soiltyp.'

      ELSE

        DO is=1,nstyp1
          IF (is <= nstyps) THEN
            READ(inch,ERR=110,END=120) label
            READ(inch,ERR=110,END=120) ((soiltyp(i,j,is),i=1,nx),j=1,ny)
            IF (myproc == 0)  WRITE(lchanl,910) label,' soiltyp.'

            READ(inch,ERR=110,END=120) label
            READ(inch,ERR=110,END=120) invar2d
            stypfrct(:,:,is) = invar2d(:,:)
            IF (myproc == 0) WRITE(lchanl,910) label,' stypfrct.'
          ELSE
            READ(inch,ERR=110,END=120) label
            IF (myproc == 0) WRITE(lchanl,910) label,'skipping soiltyp'
            READ(inch,ERR=110,END=120)
            READ(inch,ERR=110,END=120) label
            IF (myproc == 0)   &
              WRITE(lchanl,910) label,'skipping stypfrct.'
            READ(inch,ERR=110,END=120)
          ENDIF
        END DO

      END IF

      CALL fix_stypfrct_nstyp(nx,ny,nstyp1,nstyp,stypfrct)

      READ(inch,ERR=110,END=120) label
      READ(inch,ERR=110,END=120) vegtyp
      IF (myproc == 0) WRITE(lchanl,910) label,' vegtyp.'

      READ(inch,ERR=110,END=120) label
      READ(inch,ERR=110,END=120) invar2d
      lai(:,:) = invar2d(:,:)
      IF (myproc == 0) WRITE(lchanl,910) label,' lai.'

      READ(inch,ERR=110,END=120) label
      READ(inch,ERR=110,END=120) invar2d
      roufns(:,:) = invar2d(:,:)
      IF (myproc == 0) WRITE(lchanl,910) label,' roufns.'

      READ(inch,ERR=110,END=120) label
      READ(inch,ERR=110,END=120) invar2d
      veg(:,:) = invar2d(:,:)
      IF (myproc == 0) WRITE(lchanl,910) label,' veg.'

    END IF

  END IF

  IF( grdbas == 1 ) GO TO 930

  IF( varin == 1 ) THEN

    IF ( totin == 0 ) THEN
!
!-----------------------------------------------------------------------
!
!  Read in perturbations from history dump
!
!-----------------------------------------------------------------------
!
      READ(inch,ERR=110,END=120) label
      READ(inch,ERR=110,END=120) invar3d
      uprt(:,:,:) = invar3d(:,:,:)
      IF (myproc == 0) WRITE(lchanl,910) label,' uprt.'

      READ(inch,ERR=110,END=120) label
      READ(inch,ERR=110,END=120) invar3d
      vprt(:,:,:) = invar3d(:,:,:)
      IF (myproc == 0) WRITE(lchanl,910) label,' vprt.'

      READ(inch,ERR=110,END=120) label
      READ(inch,ERR=110,END=120) invar3d
      wprt(:,:,:) = invar3d(:,:,:)
      IF (myproc == 0) WRITE(lchanl,910) label,' wprt.'
!
!-----------------------------------------------------------------------
!
!  Read in scalars
!
!----------------------------------------------------------------------
!
      READ(inch,ERR=110,END=120) label
      READ(inch,ERR=110,END=120) invar3d
      ptprt(:,:,:) = invar3d(:,:,:)
      IF (myproc == 0) WRITE(lchanl,910) label,' ptprt.'

      READ(inch,ERR=110,END=120) label
      READ(inch,ERR=110,END=120) invar3d
      pprt(:,:,:) = invar3d(:,:,:)
      IF (myproc == 0) WRITE(lchanl,910) label,' pprt.'

    ELSE
!
!-----------------------------------------------------------------------
!
!  Read in total values of variables from history dump
!
!----------------------------------------------------------------------
!
      READ(inch,ERR=110,END=120) label
      READ(inch,ERR=110,END=120) invar3d
      IF (myproc == 0) WRITE(lchanl,910) label,' u.'
      DO k=1,nz-1
        DO j=1,ny-1
          DO i=1,nx
            uprt(i,j,k) = invar3d(i,j,k) - ubar(i,j,k)
          END DO
        END DO
      END DO

      READ(inch,ERR=110,END=120) label
      READ(inch,ERR=110,END=120) invar3d
      IF (myproc == 0) WRITE(lchanl,910) label,' v.'
      DO k=1,nz-1
        DO j=1,ny
          DO i=1,nx-1
            vprt(i,j,k) = invar3d(i,j,k) - vbar(i,j,k)
          END DO
        END DO
      END DO

      READ(inch,ERR=110,END=120) label
      READ(inch,ERR=110,END=120) invar3d
      wprt(:,:,:) = invar3d(:,:,:)
      IF (myproc == 0) WRITE(lchanl,910) label,' w.'

      READ(inch,ERR=110,END=120) label
      READ(inch,ERR=110,END=120) invar3d
      IF (myproc == 0) WRITE(lchanl,910) label,' pt.'
      DO k=1,nz-1
        DO j=1,ny-1
          DO i=1,nx-1
            ptprt(i,j,k) = invar3d(i,j,k) - ptbar(i,j,k)
          END DO
        END DO
      END DO

      READ(inch,ERR=110,END=120) label
      READ(inch,ERR=110,END=120) invar3d
      IF (myproc == 0) WRITE(lchanl,910) label,' p.'
      DO k=1,nz-1
        DO j=1,ny-1
          DO i=1,nx-1
            pprt(i,j,k) = invar3d(i,j,k) - pbar(i,j,k)
          END DO
        END DO
      END DO

    END IF

  END IF
!
!-----------------------------------------------------------------------
!
!  Read in moisture variables
!
!-----------------------------------------------------------------------
!
  IF( mstin == 1 ) THEN

    READ(inch,ERR=110,END=120) label
    READ(inch,ERR=110,END=120) invar3d

    IF ( totin == 0 ) THEN
      IF (myproc == 0) WRITE(lchanl,910) label,' qvprt.'
      qvprt(:,:,:) = invar3d(:,:,:)
    ELSE
      IF (myproc == 0) WRITE(lchanl,910) label,' qv.'
      DO k=1,nz-1
        DO j=1,ny-1
          DO i=1,nx-1
            qvprt(i,j,k) = invar3d(i,j,k) - qvbar(i,j,k)
          END DO
        END DO
      END DO
    END IF

    IF (intver >= intver530) THEN

      DO nqin = 1,nscalarin
        READ(inch,ERR=110,END=120) label
        READ(inch,ERR=110,END=120) invar3d
        IF (myproc == 0) WRITE(lchanl,910) label,' qscalar.'

        DO nq = 1,nscalar
          IF (nqin == nqscalarin(nq)) THEN
            qscalar(:,:,:,nq) = invar3d(:,:,:)
            EXIT
          END IF
        END DO
      END DO

    ELSE

      READ(inch,ERR=110,END=120) label
      READ(inch,ERR=110,END=120) invar3d
      IF (P_QC > 0) qscalar(:,:,:,P_QC) = invar3d(:,:,:)
      IF (myproc == 0) WRITE(lchanl,910) label,' qc.'

      READ(inch,ERR=110,END=120) label
      READ(inch,ERR=110,END=120) invar3d
      IF (P_QR > 0) qscalar(:,:,:,P_QR) = invar3d(:,:,:)
      IF (myproc == 0) WRITE(lchanl,910) label,' qr.'

    END IF

    IF( rainin == 1 ) THEN
      READ(inch,ERR=110,END=120) label
      READ(inch,ERR=110,END=120) invar2d
      raing(:,:) = invar2d(:,:)
      IF (myproc == 0) WRITE(lchanl,910) label,' raing.'

      READ(inch,ERR=110,END=120) label
      READ(inch,ERR=110,END=120) invar2d
      rainc(:,:) = invar2d(:,:)
      IF (myproc == 0) WRITE(lchanl,910) label,' rainc.'
    END IF

    IF( prcin == 1 ) THEN
      READ(inch,ERR=110,END=120) label
      READ(inch,ERR=110,END=120) invar2d
      prcrate(:,:,1) = invar2d(:,:)
      IF (myproc == 0) WRITE(lchanl,910) label,' prcrate1.'

      READ(inch,ERR=110,END=120) label
      READ(inch,ERR=110,END=120) invar2d
      prcrate(:,:,2) = invar2d(:,:)
      IF (myproc == 0) WRITE(lchanl,910) label,' prcrate2.'

      READ(inch,ERR=110,END=120) label
      READ(inch,ERR=110,END=120) invar2d
      prcrate(:,:,3) = invar2d(:,:)
      IF (myproc == 0) WRITE(lchanl,910) label,' prcrate3.'

      READ(inch,ERR=110,END=120) label
      READ(inch,ERR=110,END=120) invar2d
      prcrate(:,:,4) = invar2d(:,:)
      IF (myproc == 0) WRITE(lchanl,910) label,' prcrate4.'
    END IF

    IF( icein == 1 .AND. intver < intver530 ) THEN

      READ(inch,ERR=110,END=120) label
      READ(inch,ERR=110,END=120) invar3d
      IF (P_QI > 0) qscalar(:,:,:,P_QI) = invar3d(:,:,:)
      IF (myproc == 0) WRITE(lchanl,910) label,' qi.'

      READ(inch,ERR=110,END=120) label
      READ(inch,ERR=110,END=120) invar3d
      IF (P_QS > 0) qscalar(:,:,:,P_QS) = invar3d(:,:,:)
      IF (myproc == 0) WRITE(lchanl,910) label,' qs.'

      READ(inch,ERR=110,END=120) label
      READ(inch,ERR=110,END=120) invar3d
      IF (P_QH > 0) qscalar(:,:,:,P_QH) = invar3d(:,:,:)
      IF (myproc == 0) WRITE(lchanl,910) label,' qh.'

    END IF
  END IF

  IF( tkein == 1 ) THEN

    READ(inch,ERR=110,END=120) label
    READ(inch,ERR=110,END=120) invar3d
    tke(:,:,:) = invar3d(:,:,:)
    IF (myproc == 0) WRITE(lchanl,910) label,' tke.'

  END IF

  IF( trbin == 1 ) THEN

    READ(inch,ERR=110,END=120) label
    READ(inch,ERR=110,END=120) invar3d
    kmh(:,:,:) = invar3d(:,:,:)
    IF (myproc == 0) WRITE(lchanl,910) label,' kmh.'

    IF ( intver >= intver410 ) THEN

      READ(inch,ERR=110,END=120) label
      READ(inch,ERR=110,END=120) invar3d
      kmv(:,:,:) = invar3d(:,:,:)
      IF (myproc == 0) WRITE(lchanl,910) label,' kmv.'

    END IF ! intver

  END IF

  IF( sfcin == 1 ) THEN

    IF (nstyp1 <= 1) THEN

      IF (intver <= intver410) THEN

        READ(inch,ERR=110,END=120) label
        READ(inch,ERR=110,END=120) invar2d
        tsoil(:,:,1,0) = invar2d(:,:)
        IF (myproc == 0) WRITE(lchanl,910) label,' tsfc.'

        READ(inch,ERR=110,END=120) label
        READ(inch,ERR=110,END=120) invar2d
        tsoil(:,:,2,0) = invar2d(:,:)
        IF (myproc == 0) WRITE(lchanl,910) label,' tsoil.'

        READ(inch,ERR=110,END=120) label
        READ(inch,ERR=110,END=120) invar2d
        qsoil(:,:,1,0) = invar2d(:,:)
        IF (myproc == 0) WRITE(lchanl,910) label,' wetsfc.'

        READ(inch,ERR=110,END=120) label
        READ(inch,ERR=110,END=120) invar2d
        qsoil(:,:,2,0) = invar2d(:,:)
        IF (myproc == 0) WRITE(lchanl,910) label,' wetdp.'

      ELSE IF (intver >= intver500) THEN

        READ(inch,ERR=110,END=120) label
        READ(inch,ERR=110,END=120) invar3dsoil
        tsoil(:,:,:,0) = invar3dsoil(:,:,:)
        IF (myproc == 0) WRITE(lchanl,910) label,' tsoil.'

        READ(inch,ERR=110,END=120) label
        READ(inch,ERR=110,END=120) invar3dsoil
        qsoil(:,:,:,0) = invar3dsoil(:,:,:)
        IF (myproc == 0) WRITE(lchanl,910) label,' qsoil.'

      END IF  ! intver

      READ(inch,ERR=110,END=120) label
      READ(inch,ERR=110,END=120) invar2d
      wetcanp(:,:,0) = invar2d(:,:)
      IF (myproc == 0) WRITE(lchanl,910) label,' wetcanp.'

    ELSE

      DO is=0,nstyp1
        IF (is <= nstyps) THEN

          IF (intver <= intver410) THEN

            READ(inch,ERR=110,END=120) label
            READ(inch,ERR=110,END=120) invar2d
            tsoil(:,:,1,is) = invar2d(:,:)
            IF (myproc == 0) WRITE(lchanl,910) label,' tsfc.'

            READ(inch,ERR=110,END=120) label
            READ(inch,ERR=110,END=120) invar2d
            tsoil(:,:,2,is) = invar2d(:,:)
            IF (myproc == 0) WRITE(lchanl,910) label,' tsoil.'

            READ(inch,ERR=110,END=120) label
            READ(inch,ERR=110,END=120) invar2d
            qsoil(:,:,1,is) = invar2d(:,:)
            IF (myproc == 0) WRITE(lchanl,910) label,' wetsfc.'

            READ(inch,ERR=110,END=120) label
            READ(inch,ERR=110,END=120) invar2d
            tsoil(:,:,2,is) = invar2d(:,:)
            IF (myproc == 0) WRITE(lchanl,910) label,' wetdp.'

          ELSE IF (intver >= intver500) THEN

            READ(inch,ERR=110,END=120) label
            READ(inch,ERR=110,END=120) invar3dsoil
            tsoil(:,:,:,is) = invar3dsoil(:,:,:)
            IF (myproc == 0) WRITE(lchanl,910) label,' tsoil.'

            READ(inch,ERR=110,END=120) label
            READ(inch,ERR=110,END=120) invar3dsoil
            qsoil(:,:,:,is) = invar3dsoil(:,:,:)
            IF (myproc == 0) WRITE(lchanl,910) label,' qsoil.'

          END IF  ! intver

          READ(inch,ERR=110,END=120) label
          READ(inch,ERR=110,END=120) invar2d
          wetcanp(:,:,is) = invar2d(:,:)
          IF (myproc == 0) WRITE(lchanl,910) label,' wetcanp.'

        ELSE

          READ(inch,ERR=110,END=120) label
          READ(inch,ERR=110,END=120)
          IF (myproc == 0)  &
             WRITE(lchanl,910) label,'skipping tsoil.'

          READ(inch,ERR=110,END=120) label
          READ(inch,ERR=110,END=120)
          IF (myproc == 0)  &
             WRITE(lchanl,910) label,'skipping qsoil.'

          READ(inch,ERR=110,END=120) label
          READ(inch,ERR=110,END=120)
          IF (myproc == 0)  &
             WRITE(lchanl,910) label,'skipping wetcanp.'

        ENDIF

      END DO

    END IF

    CALL fix_soil_nstyp(nx,ny,nzsoil,nstyp1,nstyp,tsoil,qsoil,wetcanp)

    IF (snowcin == 1) THEN
      READ(inch,ERR=110,END=120) label
      READ(inch,ERR=110,END=120)
      IF (myproc == 0) WRITE(lchanl,910) label,' snowcvr -- discarding.'
    END IF

    IF (snowin == 1) THEN
      READ(inch,ERR=110,END=120) label
      READ(inch,ERR=110,END=120) invar2d
      snowdpth(:,:) = invar2d(:,:)
      IF (myproc == 0) WRITE(lchanl,910) label,' snowdpth.'
    END IF

  END IF

  IF( radin == 1 ) THEN

    READ(inch,ERR=110,END=120) label
    READ(inch,ERR=110,END=120) invar3d
    radfrc(:,:,:) = invar3d(:,:,:)
    IF (myproc == 0) WRITE(lchanl,910) label,' radfrc.'

    READ(inch,ERR=110,END=120) label
    READ(inch,ERR=110,END=120) invar2d
    radsw(:,:) = invar2d(:,:)
    IF (myproc == 0) WRITE(lchanl,910) label,' radsw.'

    READ(inch,ERR=110,END=120) label
    READ(inch,ERR=110,END=120) invar2d
    rnflx(:,:) = invar2d(:,:)
    IF (myproc == 0) WRITE(lchanl,910) label,' rnflx.'

    IF (intver <= intver410) THEN
!
! 06/28/2002 Zuwen He
!
! Do not know how to initialized radswnet, radlwin, yet,
! at this moment.
!
      radswnet = 0.
      radlwin = 0.

    ELSE IF (intver >= intver500) THEN

      READ(inch,ERR=110,END=120) label
      READ(inch,ERR=110,END=120) invar2d
      radswnet(:,:) = invar2d(:,:)
      WRITE(lchanl,910) label,' radswnet.'

      READ(inch,ERR=110,END=120) label
      READ(inch,ERR=110,END=120) invar2d
      radlwin(:,:) = invar2d(:,:)
      WRITE(lchanl,910) label,' radlwin.'

    END IF  ! intver

  END IF

  IF( flxin == 1 ) THEN

    READ(inch,ERR=110,END=120) label
    READ(inch,ERR=110,END=120) invar2d
    usflx(:,:) = invar2d(:,:)
    IF (myproc == 0) WRITE(lchanl,910) label,' usflx.'

    READ(inch,ERR=110,END=120) label
    READ(inch,ERR=110,END=120) invar2d
    vsflx(:,:) = invar2d(:,:)
    IF (myproc == 0) WRITE(lchanl,910) label,' vsflx.'

    READ(inch,ERR=110,END=120) label
    READ(inch,ERR=110,END=120) invar2d
    ptsflx(:,:) = invar2d(:,:)
    IF (myproc == 0) WRITE(lchanl,910) label,' ptsflx.'

    READ(inch,ERR=110,END=120) label
    READ(inch,ERR=110,END=120) invar2d
    qvsflx(:,:) = invar2d(:,:)
    IF (myproc == 0) WRITE(lchanl,910) label,' qvsflx.'

  END IF

  910   FORMAT(1X,'Field ',a12,' was read into array',a)

!
!-----------------------------------------------------------------------
!
!  Friendly exit message
!
!----------------------------------------------------------------------
!
  930   CONTINUE

  IF (myproc == 0)  WRITE(6,'(/a,F8.1,a/)')                             &
      ' Data at time=', time/60,' (min) were successfully read.'

  DEALLOCATE(invar1d, STAT = idummy)
  DEALLOCATE(invar2d, STAT = idummy)
  DEALLOCATE(invar3d, STAT = idummy)
  DEALLOCATE(invar3dsoil, STAT = idummy)

  ireturn = 0

  RETURN
!
!-----------------------------------------------------------------------
!
!  Error during read
!
!----------------------------------------------------------------------
!

  110   CONTINUE
  WRITE(6,'(/a/)') ' Error reading data in BINREAD'
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
  WRITE(6,'(/a/)') ' End of file reached in BINREAD'
  ireturn=2
  RETURN
END SUBROUTINE binread
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE BINREADSPLIT               ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE binreadsplit(nx,ny,nz,nzsoil,nstyps,grdbas,inch,time,x,y,z,  &
           zp,zpsoil,uprt, vprt, wprt, ptprt, pprt,                     &
           qvprt, qscalar, tke,kmh,kmv,                                 &
           ubar, vbar, wbar, ptbar, pbar, rhobar, qvbar,                &
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
!  Read in binary data set created by ARPS using history dump format
!  No. 1.
!  All data read in are located at the original staggered grid points
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yunheng Wang
!  9/03/02.
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
!    qscalar
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
!    radswnet Net shortwave radiation
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
  INCLUDE 'phycst.inc'
  INCLUDE 'mp.inc'            ! mpi parameters.
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  INTEGER :: nx,ny,nz          ! Number of grid points in 3 directions
  INTEGER :: nzsoil            ! Number of grid points in the soil

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
                               ! w-point of the soil (m)


  INTEGER :: grdbas            ! Data read flag.
  INTEGER :: inch              ! Channel number for binary reading
  REAL    :: time              ! Time in seconds of data read
                               ! from "filename"

  REAL :: uprt  (nx,ny,nz)     ! Perturbation u-velocity (m/s)
  REAL :: vprt  (nx,ny,nz)     ! Perturbation v-velocity (m/s)
  REAL :: wprt  (nx,ny,nz)     ! Perturbation w-velocity (m/s)
  REAL :: ptprt (nx,ny,nz)     ! Perturbation potential temperature (K)
  REAL :: pprt  (nx,ny,nz)     ! Perturbation pressure (Pascal)

  REAL :: qvprt (nx,ny,nz)     ! Perturbation water vapor mixing
                               ! ratio (kg/kg)

  REAL :: qscalar(nx,ny,nz,nscalar)

  REAL :: tke  (nx,ny,nz)      ! Turbulent Kinetic Energy ((m/s)**2)
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

  INTEGER :: nstyps                    ! Number of soil type
  INTEGER :: soiltyp(nx,ny,nstyps)     ! Soil type
  REAL    :: stypfrct(nx,ny,nstyps)    ! Soil type fraction
  INTEGER :: vegtyp(nx,ny)             ! Vegetation type
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
  REAL :: ptsflx(nx,ny)        ! Surface heat flux (K*kg/(m**2*s))
  REAL :: qvsflx(nx,ny)        ! Surface moisture flux (kg/(m**2*s))

  INTEGER :: ireturn           ! Return status indicator
!
!-----------------------------------------------------------------------
!
!  Parameters describing routine that wrote the gridded data
!
!-----------------------------------------------------------------------
!
! fmtver??: to label each data a version.
! intver??: an integer to allow faster comparison than fmtver??,
!           which are strings.
!
! Verion 5.00: significant change in soil variables since version 4.10.
!
  CHARACTER (LEN=40) :: fmtver320,fmtver400,fmtver410,fmtver500,fmtver530
  INTEGER            :: intver320,intver400,intver410,intver500,intver530

  PARAMETER (fmtver320='003.20 Binary Data',intver320=320)
  PARAMETER (fmtver400='004.00 Binary Data',intver400=400)
  PARAMETER (fmtver410='004.10 Binary Data',intver410=410)
  PARAMETER (fmtver500='005.00 Binary Data',intver500=500)
  PARAMETER (fmtver530='005.30 Binary Data',intver530=530)

  CHARACTER (LEN=40) :: fmtverin
  INTEGER            :: intver

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
  INTEGER :: nstyp1
  CHARACTER (LEN=12) :: label
  INTEGER :: nxin,nyin,nzin,nzsoilin
  INTEGER :: bgrdin,bbasin,bvarin,bicein,btkein,btrbin
  INTEGER :: idummy
  REAL    :: rdummy

  INTEGER :: nq, nqin
  INTEGER :: nqscalarin(nscalar)

  INTEGER :: nxlg, nylg, nzlg, nzsoillg, n3d
  INTEGER, ALLOCATABLE :: var3di(:,:,:)
  REAL,    ALLOCATABLE :: var2d(:,:),   var3d(:,:,:)
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  nxlg     = (nx-3)*nproc_x+3
  nylg     = (ny-3)*nproc_y+3
  nzlg     = nz
  nzsoillg = nzsoil
  n3d      = MAX(nz, nzsoil, nstyps+1, 4)

  ALLOCATE(var2d(nxlg,nylg))
  ALLOCATE(var3d(nxlg,nylg,n3d))
  ALLOCATE(var3di(nxlg,nylg,n3d))

  IF (myproc == 0) THEN

!-----------------------------------------------------------------------
!
!  Read header info
!
!-----------------------------------------------------------------------
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
      WRITE(6,'(/1x,a,a,a/)')                                           &
          'Incoming data format, fmtverin=',fmtverin,                   &
          ', not found. The Job stopped.'
      CALL arpsstop('arpstop called from BINREADSPLIT. ',1)
    END IF

    WRITE(6,'(/1x,a,a/)')                        &
        'Incoming data format, fmtverin=',fmtverin

    READ(inch,ERR=110,END=120) runname
    READ(inch,ERR=110,END=120) nocmnt
    IF( nocmnt > 0 ) THEN
      DO i=1,nocmnt
        READ(inch,ERR=110,END=120) cmnt(i)
      END DO
    END IF

    WRITE(6,'(//''  THE NAME OF THIS RUN IS:  '',A//)') runname

    IF( nocmnt > 0 ) THEN
      DO i=1,nocmnt
        WRITE(6,'(1x,a)') cmnt(i)
      END DO
    END IF

    READ(inch,ERR=110,END=120) time,tmunit
!
!-----------------------------------------------------------------------
!
!  Get dimensions of data in binary file and check against
!  the dimensions passed to BINREADSPLIT
!
!-----------------------------------------------------------------------
!
    IF (intver <= intver410) THEN

      READ(inch,ERR=110,END=120) nxin, nyin, nzin
      nzsoilin = 2  ! for version prior to 410, it is a two-layer soil model

    ELSE IF (intver >= intver500) THEN

      READ(inch,ERR=110,END=120) nxin, nyin, nzin,nzsoilin

    END IF
  !
  ! Data validation: dimensions
  !
    IF( nxin /= nxlg .OR. nyin /= nylg .OR. nzin /= nzlg) THEN
      WRITE(6,'(1x,a)')                                                   &
           ' Dimensions in BINREADSPLIT inconsistent with data.'
      WRITE(6,'(1x,a,3I15)') ' Read were: ', nxin, nyin, nzin
      WRITE(6,'(1x,a,3I15)') ' Expected:  ', nx, ny, nz
      WRITE(6,'(1x,a)')      ' Program aborted in BINREADSPLIT.'
      CALL arpsstop('arpstop called from BINREADSPLIT nx-ny-nz read ',1)
    END IF

    IF(nzsoilin /= nzsoillg) THEN

      IF (intver <= intver410) THEN
        WRITE(6,'(1x,a,a/,2(1x,a/))')                            &
          ' The incoming data version is ', fmtverin,            &
          ' In the input file, nzsoil must be set to 2. ',       &
          ' Program aborted in BINREADSPLIT.'
      ELSE IF (intver >= intver500) THEN
        WRITE(6,'(1x,a)') &
             ' Dimensions in BINREADSPLIT inconsistent with data.'
        WRITE(6,'(1x,a,I15)') ' Read were: ', nzsoilin
        WRITE(6,'(1x,a,I15)') ' Expected:  ', nzsoillg
        WRITE(6,'(1x,a)') ' Program aborted in BINREADSPLIT.'
      END IF
      CALL arpsstop('arpstop called from BINREADSPLIT nzsoil read ',1)

    END IF

    nxin = (nxin-3)/nproc_x+3
    nyin = (nyin-3)/nproc_y+3

  END IF

  CALL mpupdatei(intver, 1)
  CALL mpupdatec(runname, 40)
  CALL mpupdater(time, 1)
  CALL mpupdatec(tmunit, 10)

  CALL mpupdatei(nxin, 1)
  CALL mpupdatei(nyin, 1)
  CALL mpupdatei(nzin, 1)
  CALL mpupdatei(nzsoilin, 1)
!
!-----------------------------------------------------------------------
!
!  Read in flags for different data groups
!
!-----------------------------------------------------------------------
!
  IF (myproc == 0)  THEN
    IF( grdbas == 1 ) THEN ! Read grid and base state arrays

       WRITE(lchanl,'(1x,a,f8.1,a,f8.3,a/)')                            &
         'To read grid and base state data at time ', time,             &
         ' secs = ',(time/60.),' mins.'

      READ(inch,ERR=110,END=120)                                        &
           bgrdin,bbasin,bvarin,mstin,bicein,                           &
           btrbin,idummy,idummy,landin,totin,                           &
           btkein,idummy,idummy,mapproj,month,                          &
           day, year, hour, minute, second

    ELSE ! Normal data reading

       WRITE(lchanl,'(1x,a,f8.1,a,f8.3,a/)')'To read data for time:',   &
         time,' secs = ',(time/60.),' mins.'

      READ(inch,ERR=110,END=120)                                        &
           grdin,basin,varin,mstin,icein,                               &
           trbin, sfcin,rainin,landin,totin,                            &
           tkein,nscalarin,idummy,mapproj, month,                       &
           day, year, hour, minute, second

      IF (intver >= intver530 ) THEN    ! V5.3 or later
        READ(inch,ERR=110,END=120)                                      &
           p_qcin,p_qrin,p_qiin,p_qsin,p_qgin,p_qhin,                   &
           p_ncin,p_nrin,p_niin,p_nsin,p_ngin,p_nhin,                   &
           p_zrin,p_ziin,p_zsin,p_zsin,p_zhin,idummy,                   &
           idummy,idummy,idummy,idummy,idummy,idummy,                   &
           idummy,idummy,idummy,idummy,p_ccin,idummy
      ELSE
        p_qcin=0; p_qrin=0; p_qiin=0; p_qsin=0; p_qgin=0; p_qhin=0
        p_ncin=0; p_nrin=0; p_niin=0; p_nsin=0; p_nhin=0; p_ngin=0
                  p_zrin=0; p_ziin=0; p_zsin=0; p_zhin=0; p_zgin=0
        p_ccin=0

        nscalarin = 0
        IF (mstin == 1) THEN
          p_qcin=1
          p_qrin=2
          nscalarin = nscalarin + 2
          IF (icein == 2) THEN
            p_qiin=3
            p_qsin=4
            p_qhin=5
            nscalarin = nscalarin + 3
          END IF
        END IF
      END IF

      nqscalarin(:) = 0
      IF (P_QC > 0) nqscalarin(P_QC) = p_qcin
      IF (P_QR > 0) nqscalarin(P_QR) = p_qrin
      IF (P_QI > 0) nqscalarin(P_QI) = p_qiin
      IF (P_QS > 0) nqscalarin(P_QS) = p_qsin
      IF (P_QG > 0) nqscalarin(P_QG) = p_qgin
      IF (P_QH > 0) nqscalarin(P_QH) = p_qhin

      IF (P_NC > 0) nqscalarin(P_NC) = p_ncin
      IF (P_NR > 0) nqscalarin(P_NR) = p_nrin
      IF (P_NI > 0) nqscalarin(P_NI) = p_niin
      IF (P_NS > 0) nqscalarin(P_NS) = p_nsin
      IF (P_NG > 0) nqscalarin(P_NG) = p_ngin
      IF (P_NH > 0) nqscalarin(P_NH) = p_nhin

      IF (P_ZR > 0) nqscalarin(P_ZR) = p_zrin
      IF (P_ZI > 0) nqscalarin(P_ZI) = p_ziin
      IF (P_ZS > 0) nqscalarin(P_ZS) = p_zsin
      IF (P_ZG > 0) nqscalarin(P_ZG) = p_zgin
      IF (P_ZH > 0) nqscalarin(P_ZH) = p_zhin

      IF (P_CC > 0) nqscalarin(P_CC) = p_ccin
    END IF

    IF (intver < intver530) THEN

      READ(inch,ERR=110,END=120)                                        &
                      umove,vmove,xgrdorg,ygrdorg,trulat1,              &
                      trulat2,trulon,sclfct,n0rain,n0snow,              &
                      n0hail,rhosnow,rhohail,rdummy,rdummy,             &
                      tstop,thisdmp,latitud,ctrlat,ctrlon

    ELSE
      READ(inch,ERR=110,END=120)                                        &
                      umove,  vmove,  xgrdorg,ygrdorg,trulat1,          &
                      trulat2,trulon, sclfct, ntcloud, n0rain,          &
                      n0snow, n0grpl, n0hail, rhoice, rhosnow,          &
                      rhogrpl,rhohail,alpharain,alphaice,alphasnow,     &
                      alphagrpl,alphahail,rdummy, rdummy,  rdummy,      &
                      tstop,thisdmp,latitud,ctrlat,ctrlon
    END IF

    IF ( totin /= 0 ) THEN
!
!-----------------------------------------------------------------------
!
!  Read in additional parameters for ARPS history dump 4.0 or later
!  version.
!
!-----------------------------------------------------------------------
!
       READ(inch,ERR=110,END=120)                                       &
           nstyp1,  prcin, radin, flxin,snowcin,                        &
           snowin,idummy,idummy,idummy,idummy,                          &
           idummy,idummy,idummy,idummy,idummy,                          &
           idummy,idummy,idummy,idummy,idummy

      IF ( nstyp1 < 1 ) THEN
        nstyp1 = 1
      END IF

      READ(inch,ERR=110,END=120)                                        &
           rdummy,rdummy,rdummy,rdummy,rdummy,                          &
           rdummy,rdummy,rdummy,rdummy,rdummy,                          &
           rdummy,rdummy,rdummy,rdummy,rdummy,                          &
           rdummy,rdummy,rdummy,rdummy,rdummy
    END IF

  END IF  ! myproc == 0
  CALL mpupdatei(mstin,1)
  CALL mpupdatei(landin,1)
  CALL mpupdatei(totin,1)
  CALL mpupdatei(mapproj,1)
  CALL mpupdatei(month,1)
  CALL mpupdatei(day,1)
  CALL mpupdatei(year,1)
  CALL mpupdatei(hour,1)
  CALL mpupdatei(minute,1)
  CALL mpupdatei(second,1)
  IF(grdbas == 1) THEN
    CALL mpupdatei(bgrdin,1)
    CALL mpupdatei(bbasin,1)
    CALL mpupdatei(bvarin,1)
    CALL mpupdatei(btrbin,1)
    CALL mpupdatei(btkein,1)
  ELSE
    CALL mpupdatei(grdin,1)
    CALL mpupdatei(basin,1)
    CALL mpupdatei(varin,1)
    CALL mpupdatei(trbin,1)
    CALL mpupdatei(tkein,1)
    CALL mpupdatei(icein,1)
    CALL mpupdatei(sfcin,1)
    CALL mpupdatei(rainin,1)

    CALL mpupdatei(nscalarin,1)
    CALL mpupdatei(nqscalarin,nscalar)
  END IF
  CALL mpupdater(umove,1)
  CALL mpupdater(vmove,1)
  CALL mpupdater(xgrdorg,1)
  CALL mpupdater(ygrdorg,1)
  CALL mpupdater(trulat1,1)
  CALL mpupdater(trulat2,1)
  CALL mpupdater(trulon,1)
  CALL mpupdater(sclfct,1)
  CALL mpupdater(tstop,1)
  CALL mpupdater(thisdmp,1)
  CALL mpupdater(latitud,1)
  CALL mpupdater(ctrlat,1)
  CALL mpupdater(ctrlon,1)
  CALL mpupdater(n0rain,1)
  CALL mpupdater(n0snow,1)
  CALL mpupdater(n0hail,1)
  CALL mpupdater(rhosnow,1)
  CALL mpupdater(rhohail,1)
  IF(totin /= 0) THEN
    CALL mpupdatei(nstyp1,1)
    CALL mpupdatei(prcin,1)
    CALL mpupdatei(radin,1)
    CALL mpupdatei(flxin,1)
    CALL mpupdatei(snowcin,1)
    CALL mpupdatei(snowin,1)
  END IF
  IF (intver >= intver530) THEN
    CALL mpupdater(ntcloud,1)
    CALL mpupdater(n0grpl, 1)
    CALL mpupdater(rhoice, 1)
    CALL mpupdater(rhogrpl,1)
    CALL mpupdater(alpharain,1)
    CALL mpupdater(alphasnow,1)
    CALL mpupdater(alphaice, 1)
    CALL mpupdater(alphagrpl,1)
    CALL mpupdater(alphahail,1)
  END IF
!
!-----------------------------------------------------------------------
!
!  Read in x, y, z and zp arrays.
!
!----------------------------------------------------------------------
!
  IF( grdin == 1 .OR. grdbas == 1 ) THEN
    IF(myproc == 0) THEN
      READ(inch,ERR=110,END=120) label
      READ(inch,ERR=110,END=120) (var2d(i,1),i=1,nxlg)
      WRITE(lchanl,910) label,' x.'
    END IF
    CALL mpisplit1dx(var2d(:,1),nx,x)

    IF(myproc == 0) THEN
      READ(inch,ERR=110,END=120) label
      READ(inch,ERR=110,END=120) (var2d(1,j),j=1,nylg)
      WRITE(lchanl,910) label,' y.'
    END IF
    CALL mpisplit1dy(var2d(1,:),ny,y)

    IF(myproc == 0) THEN
      READ(inch,ERR=110,END=120) label
      READ(inch,ERR=110,END=120) z
      WRITE(lchanl,910) label,' z.'
    END IF
    CALL mpupdater(z,nzlg)

    IF(myproc == 0) THEN
      READ(inch,ERR=110,END=120) label
      READ(inch,ERR=110,END=120) (((var3d(i,j,k),i=1,nxlg),j=1,nylg),k=1,nzlg)
      WRITE(lchanl,910) label,' zp.'
    END IF
    CALL mpisplit3d(var3d,nx,ny,nz,zp)

    IF (intver <= intver410) THEN
!
! nzsoil must equal to 2, 06/28/2002, Zuwen
! assume zpsoil(,,2) is one meter below the surface.
!
      DO j=1,ny
        DO i=1,nx
          zpsoil(i,j,1)=zp(i,j,2)
          zpsoil(i,j,2)=zpsoil(i,j,1)-1.
        END DO
      END DO

      IF (myproc == 0) THEN
        WRITE(lchanl,'(1x,a)') ' Assign zpsoil. '
        WRITE(lchanl,'(1x,a)') ' Assume zpsoil(,,1) is zp(:,:,2). '
        WRITE(lchanl,'(1x,a)') ' Assume zpsoil(,,2) is zp(:,:,2)-1. '
      END IF

    ELSE IF (intver >= intver500) THEN

      IF(myproc == 0) THEN
        READ(inch,ERR=110,END=120) label
        READ(inch,ERR=110,END=120) (((var3d(i,j,k),i=1,nxlg),j=1,nylg),k=1,nzsoillg)
        WRITE(lchanl,910) label,' zpsoil.'
      END IF
      CALL mpisplit3d(var3d,nx,ny,nzsoil,zpsoil)

    END IF  ! intver

  END IF  ! grdin
!
!-----------------------------------------------------------------------
!
!  Read in base state fields
!
!----------------------------------------------------------------------
!
  IF( basin == 1 .OR. grdbas == 1 ) THEN

    IF(myproc == 0) THEN
      READ(inch,ERR=110,END=120) label
      READ(inch,ERR=110,END=120) (((var3d(i,j,k),i=1,nxlg),j=1,nylg),k=1,nzlg)
      WRITE(lchanl,910) label,' ubar.'
    END IF
    CALL mpisplit3d(var3d,nx,ny,nz,ubar)

    IF(myproc == 0) THEN
      READ(inch,ERR=110,END=120) label
      READ(inch,ERR=110,END=120) (((var3d(i,j,k),i=1,nxlg),j=1,nylg),k=1,nzlg)
      WRITE(lchanl,910) label,' vbar.'
    END IF
    CALL mpisplit3d(var3d,nx,ny,nz,vbar)

    IF(myproc == 0) THEN
      READ(inch,ERR=110,END=120) label
      READ(inch,ERR=110,END=120) (((var3d(i,j,k),i=1,nxlg),j=1,nylg),k=1,nzlg)
      WRITE(lchanl,910) label,' wbar.'
    END IF
    CALL mpisplit3d(var3d,nx,ny,nz,wbar)

    IF(myproc == 0) THEN
      READ(inch,ERR=110,END=120) label
      READ(inch,ERR=110,END=120) (((var3d(i,j,k),i=1,nxlg),j=1,nylg),k=1,nzlg)
      WRITE(lchanl,910) label,' ptbar.'
    END IF
    CALL mpisplit3d(var3d,nx,ny,nz,ptbar)

    IF(myproc == 0) THEN
      READ(inch,ERR=110,END=120) label
      READ(inch,ERR=110,END=120) (((var3d(i,j,k),i=1,nxlg),j=1,nylg),k=1,nzlg)
      WRITE(lchanl,910) label,' pbar.'
    END IF
    CALL mpisplit3d(var3d,nx,ny,nz,pbar)

    IF( mstin == 1) THEN
      IF(myproc == 0) THEN
        READ(inch,ERR=110,END=120) label
        READ(inch,ERR=110,END=120) (((var3d(i,j,k),i=1,nxlg),j=1,nylg),k=1,nzlg)
        WRITE(lchanl,910) label,' qvbar.'
      END IF
      CALL mpisplit3d(var3d,nx,ny,nz,qvbar)
    END IF

    IF (landin == 1) THEN

      IF (nstyp1 <= 1) THEN

        IF(myproc == 0) THEN
          READ(inch,ERR=110,END=120) label
          READ(inch,ERR=110,END=120) ((var3di(i,j,1),i=1,nxlg),j=1,nylg)
          WRITE(lchanl,910) label,' soiltyp.'
        END IF
        CALL mpisplit3di(var3di,nx,ny,1,soiltyp(:,:,1))

      ELSE

        IF(myproc == 0) THEN
          var3di(:,:,:) = 0
          var3d(:,:,:)  = 0.0
          DO is=1,nstyp1
            IF (is <= nstyps) THEN
              READ(inch,ERR=110,END=120) label
              READ(inch,ERR=110,END=120) ((var3di(i,j,is),i=1,nxlg),j=1,nylg)
              WRITE(lchanl,910) label,' soiltyp.'

              READ(inch,ERR=110,END=120) label
              READ(inch,ERR=110,END=120) ((var3d(i,j,is),i=1,nxlg),j=1,nylg)
              WRITE(lchanl,910) label,' stypfrct.'
            ELSE
              READ(inch,ERR=110,END=120) label
              WRITE(lchanl,910) label,'skipping soiltyp'
              READ(inch,ERR=110,END=120)
              READ(inch,ERR=110,END=120) label
              WRITE(lchanl,910) label,'skipping stypfrct.'
              READ(inch,ERR=110,END=120)
            END IF
          END DO
        ENDIF
        CALL mpisplit3di(var3di,nx,ny,nstyps,soiltyp)
        CALL mpisplit3d(var3d,nx,ny,nstyps,stypfrct)
      END IF

      CALL fix_stypfrct_nstyp(nx,ny,nstyp1,nstyp,stypfrct)

      IF(myproc == 0) THEN
        READ(inch,ERR=110,END=120) label
        READ(inch,ERR=110,END=120) ((var3di(i,j,1),i=1,nxlg),j=1,nylg)
        WRITE(lchanl,910) label,' vegtyp.'
      END IF
      CALL mpisplit2di(var3di(:,:,1),nx,ny,vegtyp)

      IF(myproc == 0) THEN
        READ(inch,ERR=110,END=120) label
        READ(inch,ERR=110,END=120) var2d
        WRITE(lchanl,910) label,' lai.'
      END IF
      CALL mpisplit2d(var2d,nx,ny,lai)

      IF(myproc == 0) THEN
        READ(inch,ERR=110,END=120) label
        READ(inch,ERR=110,END=120) var2d
        WRITE(lchanl,910) label,' roufns.'
      END IF
      CALL mpisplit2d(var2d,nx,ny,roufns)

      IF(myproc == 0) THEN
        READ(inch,ERR=110,END=120) label
        READ(inch,ERR=110,END=120) var2d
        WRITE(lchanl,910) label,' veg.'
      END IF
      CALL mpisplit2d(var2d,nx,ny,veg)

    END IF

  END IF

  IF( grdbas == 1 ) GO TO 930

  IF( varin == 1 ) THEN

    IF ( totin == 0 ) THEN
!
!-----------------------------------------------------------------------
!
!  Read in perturbations from history dump
!
!-----------------------------------------------------------------------
!
      IF(myproc == 0) THEN
        READ(inch,ERR=110,END=120) label
        READ(inch,ERR=110,END=120) (((var3d(i,j,k),i=1,nxlg),j=1,nylg),k=1,nzlg)
        WRITE(lchanl,910) label,' uprt.'
      END IF
      CALL mpisplit3d(var3d,nx,ny,nz,uprt)

      IF(myproc == 0) THEN
        READ(inch,ERR=110,END=120) label
        READ(inch,ERR=110,END=120) (((var3d(i,j,k),i=1,nxlg),j=1,nylg),k=1,nzlg)
        WRITE(lchanl,910) label,' vprt.'
      END IF
      CALL mpisplit3d(var3d,nx,ny,nz,vprt)

      IF(myproc == 0) THEN
        READ(inch,ERR=110,END=120) label
        READ(inch,ERR=110,END=120) (((var3d(i,j,k),i=1,nxlg),j=1,nylg),k=1,nzlg)
        WRITE(lchanl,910) label,' wprt.'
      END IF
      CALL mpisplit3d(var3d,nx,ny,nz,wprt)
!
!-----------------------------------------------------------------------
!
!  Read in scalars
!
!----------------------------------------------------------------------
!
      IF(myproc == 0) THEN
        READ(inch,ERR=110,END=120) label
        READ(inch,ERR=110,END=120) (((var3d(i,j,k),i=1,nxlg),j=1,nylg),k=1,nzlg)
        WRITE(lchanl,910) label,' ptprt.'
      END IF
      CALL mpisplit3d(var3d,nx,ny,nz,ptprt)

      IF(myproc == 0) THEN
        READ(inch,ERR=110,END=120) label
        READ(inch,ERR=110,END=120) (((var3d(i,j,k),i=1,nxlg),j=1,nylg),k=1,nzlg)
        WRITE(lchanl,910) label,' pprt.'
      END IF
      CALL mpisplit3d(var3d,nx,ny,nz,pprt)

    ELSE
!
!-----------------------------------------------------------------------
!
!  Read in total values of variables from history dump
!
!----------------------------------------------------------------------
!
      IF(myproc == 0) THEN
        READ(inch,ERR=110,END=120) label
        READ(inch,ERR=110,END=120) (((var3d(i,j,k),i=1,nxlg),j=1,nylg),k=1,nzlg)
        WRITE(lchanl,910) label,' u.'
      END IF
      CALL mpisplit3d(var3d,nx,ny,nz,uprt)
      DO k=1,nz-1
        DO j=1,ny-1
          DO i=1,nx
            uprt(i,j,k) = uprt(i,j,k) - ubar(i,j,k)
          END DO
        END DO
      END DO

      IF(myproc == 0) THEN
        READ(inch,ERR=110,END=120) label
        READ(inch,ERR=110,END=120) (((var3d(i,j,k),i=1,nxlg),j=1,nylg),k=1,nzlg)
        WRITE(lchanl,910) label,' v.'
      END IF
      CALL mpisplit3d(var3d,nx,ny,nz,vprt)
      DO k=1,nz-1
        DO j=1,ny
          DO i=1,nx-1
            vprt(i,j,k) = vprt(i,j,k) - vbar(i,j,k)
          END DO
        END DO
      END DO

      IF(myproc == 0) THEN
        READ(inch,ERR=110,END=120) label
        READ(inch,ERR=110,END=120) (((var3d(i,j,k),i=1,nxlg),j=1,nylg),k=1,nzlg)
        WRITE(lchanl,910) label,' w.'
      END IF
      CALL mpisplit3d(var3d,nx,ny,nz,wprt)

      IF(myproc == 0) THEN
        READ(inch,ERR=110,END=120) label
        READ(inch,ERR=110,END=120) (((var3d(i,j,k),i=1,nxlg),j=1,nylg),k=1,nzlg)
        WRITE(lchanl,910) label,' pt.'
      END IF
      CALL mpisplit3d(var3d,nx,ny,nz,ptprt)
      DO k=1,nz-1
        DO j=1,ny-1
          DO i=1,nx-1
            ptprt(i,j,k) = ptprt(i,j,k) - ptbar(i,j,k)
          END DO
        END DO
      END DO

      IF(myproc == 0) THEN
        READ(inch,ERR=110,END=120) label
        READ(inch,ERR=110,END=120) (((var3d(i,j,k),i=1,nxlg),j=1,nylg),k=1,nzlg)
        WRITE(lchanl,910) label,' p.'
      END IF
      CALL mpisplit3d(var3d,nx,ny,nz,pprt)
      DO k=1,nz-1
        DO j=1,ny-1
          DO i=1,nx-1
            pprt(i,j,k) = pprt(i,j,k) - pbar(i,j,k)
          END DO
        END DO
      END DO

    END IF

  END IF
!
!-----------------------------------------------------------------------
!
!  Read in moisture variables
!
!-----------------------------------------------------------------------
!
  IF( mstin == 1 ) THEN

    IF ( totin == 0 ) THEN

      IF(myproc == 0) THEN
        READ(inch,ERR=110,END=120) label
        READ(inch,ERR=110,END=120) (((var3d(i,j,k),i=1,nxlg),j=1,nylg),k=1,nzlg)
        WRITE(lchanl,910) label,' qvprt.'
      END IF
      CALL mpisplit3d(var3d,nx,ny,nz,qvprt)

    ELSE

      IF(myproc == 0) THEN
        READ(inch,ERR=110,END=120) label
        READ(inch,ERR=110,END=120) (((var3d(i,j,k),i=1,nxlg),j=1,nylg),k=1,nzlg)
        WRITE(lchanl,910) label,' qv.'
      END IF
      CALL mpisplit3d(var3d,nx,ny,nz,qvprt)
      DO k=1,nz-1
        DO j=1,ny-1
          DO i=1,nx-1
            qvprt(i,j,k) = qvprt(i,j,k) - qvbar(i,j,k)
          END DO
        END DO
      END DO

    END IF

    IF (intver >= intver530) THEN

      DO nqin = 1,nscalarin

        IF(myproc == 0) THEN
          READ(inch,ERR=110,END=120) label
          READ(inch,ERR=110,END=120) (((var3d(i,j,k),i=1,nxlg),j=1,nylg),k=1,nzlg)
          WRITE(lchanl,910) label,' qscalar.'
        END IF
        CALL mpupdatec(label,12)

        DO nq = 1,nscalar
          IF (nqin == nqscalarin(nq)) THEN
            CALL mpisplit3d(var3d,nx,ny,nz,qscalar(:,:,:,nq))
            EXIT
          END IF
        END DO

      END DO

    ELSE

      IF(myproc == 0) THEN
        READ(inch,ERR=110,END=120) label
        READ(inch,ERR=110,END=120) (((var3d(i,j,k),i=1,nxlg),j=1,nylg),k=1,nzlg)
        WRITE(lchanl,910) label,' qc.'
      END IF
      IF (P_QC > 0) CALL mpisplit3d(var3d,nx,ny,nz,qscalar(:,:,:,P_QC))

      IF(myproc == 0) THEN
        READ(inch,ERR=110,END=120) label
        READ(inch,ERR=110,END=120) (((var3d(i,j,k),i=1,nxlg),j=1,nylg),k=1,nzlg)
        WRITE(lchanl,910) label,' qr.'
      END IF
      IF (P_QR > 0) CALL mpisplit3d(var3d,nx,ny,nz,qscalar(:,:,:,P_QR))

    END IF

    IF( rainin == 1 ) THEN
      IF(myproc == 0) THEN
        READ(inch,ERR=110,END=120) label
        READ(inch,ERR=110,END=120) var2d
        WRITE(lchanl,910) label,' raing.'
      END IF
      CALL mpisplit2d(var2d,nx,ny,raing)

      IF(myproc == 0) THEN
        READ(inch,ERR=110,END=120) label
        READ(inch,ERR=110,END=120) var2d
        WRITE(lchanl,910) label,' rainc.'
      END IF
      CALL mpisplit2d(var2d,nx,ny,rainc)
    END IF

    IF( prcin == 1 ) THEN
      IF(myproc == 0) THEN
        READ(inch,ERR=110,END=120) label
        READ(inch,ERR=110,END=120) ((var3d(i,j,1),i=1,nxlg),j=1,nylg)
        WRITE(lchanl,910) label,' prcrate1.'

        READ(inch,ERR=110,END=120) label
        READ(inch,ERR=110,END=120) ((var3d(i,j,2),i=1,nxlg),j=1,nylg)
        WRITE(lchanl,910) label,' prcrate2.'

        READ(inch,ERR=110,END=120) label
        READ(inch,ERR=110,END=120) ((var3d(i,j,3),i=1,nxlg),j=1,nylg)
        WRITE(lchanl,910) label,' prcrate3.'

        READ(inch,ERR=110,END=120) label
        READ(inch,ERR=110,END=120) ((var3d(i,j,4),i=1,nxlg),j=1,nylg)
        WRITE(lchanl,910) label,' prcrate4.'
      END IF
      CALL mpisplit3d(var3d,nx,ny,4,prcrate)
    END IF

    IF( icein == 1 .AND. intver < intver530 ) THEN

      IF(myproc == 0) THEN
        READ(inch,ERR=110,END=120) label
        READ(inch,ERR=110,END=120) (((var3d(i,j,k),i=1,nxlg),j=1,nylg),k=1,nzlg)
        WRITE(lchanl,910) label,' qi.'
      END IF
      IF (P_QI > 0) CALL mpisplit3d(var3d,nx,ny,nz,qscalar(:,:,:,P_QI))

      IF(myproc == 0) THEN
        READ(inch,ERR=110,END=120) label
        READ(inch,ERR=110,END=120) (((var3d(i,j,k),i=1,nxlg),j=1,nylg),k=1,nzlg)
        WRITE(lchanl,910) label,' qs.'
      END IF
      IF (P_QS > 0) CALL mpisplit3d(var3d,nx,ny,nz,qscalar(:,:,:,P_QS))

      IF(myproc == 0) THEN
        READ(inch,ERR=110,END=120) label
        READ(inch,ERR=110,END=120) (((var3d(i,j,k),i=1,nxlg),j=1,nylg),k=1,nzlg)
        WRITE(lchanl,910) label,' qh.'
      END IF
      IF (P_QH > 0) CALL mpisplit3d(var3d,nx,ny,nz,qscalar(:,:,:,P_QH))

    END IF
  END IF

  IF( tkein == 1 ) THEN

    IF(myproc == 0) THEN
      READ(inch,ERR=110,END=120) label
      READ(inch,ERR=110,END=120) (((var3d(i,j,k),i=1,nxlg),j=1,nylg),k=1,nzlg)
      WRITE(lchanl,910) label,' tke.'
    END IF
    CALL mpisplit3d(var3d,nx,ny,nz,tke)

  END IF

  IF( trbin == 1 ) THEN

    IF(myproc == 0) THEN
      READ(inch,ERR=110,END=120) label
      READ(inch,ERR=110,END=120) (((var3d(i,j,k),i=1,nxlg),j=1,nylg),k=1,nzlg)
      WRITE(lchanl,910) label,' kmh.'
    END IF
    CALL mpisplit3d(var3d,nx,ny,nz,kmh)

    IF ( intver >= intver410 ) THEN

      IF(myproc == 0) THEN
        READ(inch,ERR=110,END=120) label
        READ(inch,ERR=110,END=120) (((var3d(i,j,k),i=1,nxlg),j=1,nylg),k=1,nzlg)
        WRITE(lchanl,910) label,' kmv.'
      END IF
      CALL mpisplit3d(var3d,nx,ny,nz,kmv)

    END IF ! intver

  END IF

  IF( sfcin == 1 ) THEN

    IF (nstyp1 <= 1) THEN

      IF (intver <= intver410) THEN

        IF(myproc == 0) THEN
          READ(inch,ERR=110,END=120) label
          READ(inch,ERR=110,END=120) var2d
          WRITE(lchanl,910) label,' tsfc.'
        END IF
        CALL mpisplit2d(var2d,nx,ny,tsoil(:,:,1,0))

        IF(myproc == 0) THEN
          READ(inch,ERR=110,END=120) label
          READ(inch,ERR=110,END=120) var2d
          WRITE(lchanl,910) label,' tsoil.'
        END IF
        CALL mpisplit2d(var2d,nx,ny,tsoil(:,:,2,0))

        IF(myproc == 0) THEN
          READ(inch,ERR=110,END=120) label
          READ(inch,ERR=110,END=120) var2d
          WRITE(lchanl,910) label,' wetsfc.'
        END IF
        CALL mpisplit2d(var2d,nx,ny,qsoil(:,:,1,0))

        IF(myproc == 0) THEN
          READ(inch,ERR=110,END=120) label
          READ(inch,ERR=110,END=120) var2d
          WRITE(lchanl,910) label,' wetdp.'
        END IF
        CALL mpisplit2d(var2d,nx,ny,qsoil(:,:,2,0))

      ELSE IF (intver >= intver500) THEN

        IF(myproc == 0) THEN
          READ(inch,ERR=110,END=120) label
          READ(inch,ERR=110,END=120) (((var3d(i,j,k),i=1,nxlg),j=1,nylg),k=1,nzsoillg)
          WRITE(lchanl,910) label,' tsoil.'
        END IF
        CALL mpisplit3d(var3d,nx,ny,nzsoil,tsoil(:,:,:,0))

        IF(myproc == 0) THEN
          READ(inch,ERR=110,END=120) label
          READ(inch,ERR=110,END=120) (((var3d(i,j,k),i=1,nxlg),j=1,nylg),k=1,nzsoillg)
          WRITE(lchanl,910) label,' qsoil.'
        END IF
        CALL mpisplit3d(var3d,nx,ny,nzsoil,qsoil(:,:,:,0))

      END IF  ! intver

      IF(myproc == 0) THEN
        READ(inch,ERR=110,END=120) label
        READ(inch,ERR=110,END=120) var2d
        WRITE(lchanl,910) label,' wetcanp.'
      END IF
      CALL mpisplit2d(var2d,nx,ny,wetcanp(:,:,0))

    ELSE

      DO is=0,nstyp1

        IF (is <= nstyps) THEN

          IF (intver <= intver410) THEN

            IF(myproc == 0) THEN
              READ(inch,ERR=110,END=120) label
              READ(inch,ERR=110,END=120) var2d
              WRITE(lchanl,910) label,' tsfc.'
            END IF
            CALL mpisplit2d(var2d,nx,ny,tsoil(:,:,1,is))

            IF(myproc == 0) THEN
            READ(inch,ERR=110,END=120) label
            READ(inch,ERR=110,END=120) var2d
            WRITE(lchanl,910) label,' tsoil.'
            END IF
            CALL mpisplit2d(var2d,nx,ny,tsoil(:,:,2,is))

            IF(myproc == 0) THEN
              READ(inch,ERR=110,END=120) label
              READ(inch,ERR=110,END=120) var2d
              WRITE(lchanl,910) label,' wetsfc.'
            END IF
            CALL mpisplit2d(var2d,nx,ny,qsoil(:,:,1,is))

            IF(myproc == 0) THEN
            READ(inch,ERR=110,END=120) label
            READ(inch,ERR=110,END=120) var2d
            WRITE(lchanl,910) label,' wetdp.'
            END IF
            CALL mpisplit2d(var2d,nx,ny,qsoil(:,:,2,is))

          ELSE IF (intver >= intver500) THEN

            IF(myproc == 0) THEN
              READ(inch,ERR=110,END=120) label
              READ(inch,ERR=110,END=120) (((var3d(i,j,k),i=1,nxlg),j=1,nylg),k=1,nzsoillg)
              WRITE(lchanl,910) label,' tsoil.'
            END IF
            CALL mpisplit3d(var3d,nx,ny,nzsoil,tsoil(:,:,:,is))

            IF(myproc == 0) THEN
            READ(inch,ERR=110,END=120) label
            READ(inch,ERR=110,END=120) (((var3d(i,j,k),i=1,nxlg),j=1,nylg),k=1,nzsoillg)
            WRITE(lchanl,910) label,' qsoil.'
            END IF
            CALL mpisplit3d(var3d,nx,ny,nzsoil,qsoil(:,:,:,is))

          END IF  ! intver

          IF(myproc == 0) THEN
          READ(inch,ERR=110,END=120) label
          READ(inch,ERR=110,END=120) var2d
          WRITE(lchanl,910) label,' wetcanp.'
          END IF
          CALL mpisplit2d(var2d,nx,ny,wetcanp(:,:,is))

        ELSE

          IF(myproc == 0) THEN
            READ(inch,ERR=110,END=120) label
            READ(inch,ERR=110,END=120)
            WRITE(lchanl,910) label,'skipping tsoil.'

            READ(inch,ERR=110,END=120) label
            READ(inch,ERR=110,END=120)
            WRITE(lchanl,910) label,'skipping qsoil.'

            READ(inch,ERR=110,END=120) label
            READ(inch,ERR=110,END=120)
            WRITE(lchanl,910) label,'skipping wetcanp.'
          END IF

        ENDIF

      END DO

    END IF

    CALL fix_soil_nstyp(nx,ny,nzsoil,nstyp1,nstyp,tsoil,qsoil,wetcanp)

    IF (snowcin == 1) THEN
      IF(myproc == 0) THEN
        READ(inch,ERR=110,END=120) label
        READ(inch,ERR=110,END=120)
        WRITE(lchanl,910) label,' snowcvr -- discarding.'
      END IF
    END IF

    IF (snowin == 1) THEN
      IF(myproc == 0) THEN
        READ(inch,ERR=110,END=120) label
        READ(inch,ERR=110,END=120) var2d
        WRITE(lchanl,910) label,' snowdpth.'
      END IF
      CALL mpisplit2d(var2d,nx,ny,snowdpth)
    END IF

  END IF

  IF( radin == 1 ) THEN

    IF(myproc == 0) THEN
      READ(inch,ERR=110,END=120) label
      READ(inch,ERR=110,END=120) (((var3d(i,j,k),i=1,nxlg),j=1,nylg),k=1,nzlg)
      WRITE(lchanl,910) label,' radfrc.'
    END IF
    CALL mpisplit3d(var3d,nx,ny,nz,radfrc)

    IF(myproc == 0) THEN
    READ(inch,ERR=110,END=120) label
    READ(inch,ERR=110,END=120) var2d
    WRITE(lchanl,910) label,' radsw.'
    END IF
    CALL mpisplit2d(var2d,nx,ny,radsw)

    IF(myproc == 0) THEN
      READ(inch,ERR=110,END=120) label
      READ(inch,ERR=110,END=120) var2d
      WRITE(lchanl,910) label,' rnflx.'
    END IF
    CALL mpisplit2d(var2d,nx,ny,rnflx)

    IF (intver <= intver410) THEN
!
! Do not know how to initialized radswnet, radlwin, yet,
! at this moment.
!
      radswnet = 0.
      radlwin = 0.

    ELSE IF (intver >= intver500) THEN

      IF(myproc == 0) THEN
        READ(inch,ERR=110,END=120) label
        READ(inch,ERR=110,END=120) var2d
        WRITE(lchanl,910) label,' radswnet.'
      END IF
      CALL mpisplit2d(var2d,nx,ny,radswnet)

      IF(myproc == 0) THEN
        READ(inch,ERR=110,END=120) label
        READ(inch,ERR=110,END=120) var2d
        WRITE(lchanl,910) label,' radlwin.'
      END IF
      CALL mpisplit2d(var2d,nx,ny,radlwin)

    END IF  ! intver

  END IF

  IF( flxin == 1 ) THEN

    IF(myproc == 0) THEN
      READ(inch,ERR=110,END=120) label
      READ(inch,ERR=110,END=120) var2d
      WRITE(lchanl,910) label,' usflx.'
    END IF
    CALL mpisplit2d(var2d,nx,ny,usflx)

    IF(myproc == 0) THEN
      READ(inch,ERR=110,END=120) label
      READ(inch,ERR=110,END=120) var2d
      WRITE(lchanl,910) label,' vsflx.'
    END IF
    CALL mpisplit2d(var2d,nx,ny,vsflx)

    IF(myproc == 0) THEN
      READ(inch,ERR=110,END=120) label
      READ(inch,ERR=110,END=120) var2d
      WRITE(lchanl,910) label,' ptsflx.'
    END IF
    CALL mpisplit2d(var2d,nx,ny,ptsflx)

    IF(myproc == 0) THEN
      READ(inch,ERR=110,END=120) label
      READ(inch,ERR=110,END=120) var2d
      WRITE(lchanl,910) label,' qvsflx.'
    END IF
    CALL mpisplit2d(var2d,nx,ny,qvsflx)

  END IF

  910   FORMAT(1X,'Field ',a12,' was read into array',a)
!
!-----------------------------------------------------------------------
!
!  Friendly exit message
!
!----------------------------------------------------------------------
!
  930   CONTINUE

  IF (myproc == 0)  &
     WRITE(6,'(/a,F8.1,a/)')                                          &
        ' Data at time=', time/60,' (min) were successfully read.'

  ireturn = 0

  DEALLOCATE(var2d, var3d, var3di)

  RETURN
!
!-----------------------------------------------------------------------
!
!  Error during read
!
!----------------------------------------------------------------------
!

  110   CONTINUE
  WRITE(6,'(/a/)') ' Error reading data in BINREADSPLIT'
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
  WRITE(6,'(/a/)') ' End of file reached in BINREADSPLIT'
  ireturn=2
  RETURN
END SUBROUTINE binreadsplit
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE BN2READ                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE bn2read(nx,ny,nz,nzsoil,nstyps,grdbas,inch,time,x,y,z,zp,    &
           zpsoil,uprt, vprt, wprt, ptprt, pprt,                        &
           qvprt, qc, qr, qi, qs, qh, tke,kmh,kmv,                      &
           ubar, vbar, wbar, ptbar, pbar, rhobar, qvbar,                &
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
!  Read in binary data set created by ARPS using history dump format
!  No.2.
!  All data read in are located at the original staggered grid points
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  2/27/92.
!
!  MODIFICATION HISTORY:
!
!  6/08/92  Added full documentation (K. Brewster)
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
!  Added map projection parameters into the second binary dumping
!
!  12/09/1998 (Donghai Wang)
!  Added the snow cover.
!
!  05/15/2002 (J. Brotzge)
!  Added variables to allow for multiple soil schemes.
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
!    radswnet Net shortwave radiation
!    radlwin  Incoming longwave radiation
!
!    usflx    Surface flux of u-momentum (kg/(m*s**2))
!    vsflx    Surface flux of v-momentum (kg/(m*s**2))
!    ptsflx   Surface heat flux (K*kg/(m**2 * s ))
!    qvsflx   Surface moisture flux of (kg/(m**2 * s))
!
!    ireturn  Return status indicator
!             =0   successful read of all data
!             =1   error reading data
!             =2   end-of-file reached during read attempt
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx,ny,nz          ! Number of grid points in 3 directions
  INTEGER :: nzsoil            ! Number of grid points in 3 directions

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
                               ! w-point of the soil (m)

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
  REAL :: qc    (nx,ny,nz)     ! Cloud water mixing ratio (kg/kg)
  REAL :: qr    (nx,ny,nz)     ! Rain water mixing ratio (kg/kg)
  REAL :: qi    (nx,ny,nz)     ! Cloud ice mixing ratio (kg/kg)
  REAL :: qs    (nx,ny,nz)     ! Snow mixing ratio (kg/kg)
  REAL :: qh    (nx,ny,nz)     ! Hail mixing ratio (kg/kg)
  REAL :: tke  (nx,ny,nz)      ! Turbulent Kinetic Energy ((m/s)**2)
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

  INTEGER :: nstyps                   ! Number of soil type
  INTEGER :: soiltyp(nx,ny,nstyps)    ! Soil type
  REAL :: stypfrct(nx,ny,nstyps)   ! Soil type
  INTEGER :: vegtyp(nx,ny)            ! Vegetation type
  REAL :: lai    (nx,ny)           ! Leaf Area Index
  REAL :: roufns (nx,ny)           ! Surface roughness
  REAL :: veg    (nx,ny)           ! Vegetation fraction

  REAL :: tsoil (nx,ny,nzsoil,0:nstyps) ! Soil temperature (K)
  REAL :: qsoil (nx,ny,nzsoil,0:nstyps) ! Soil moisture (m**3/m**3)
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
  REAL :: radswnet(nx,ny)      ! Net shortwave radiation
  REAL :: radlwin(nx,ny)       ! Incoming longwave radiation

  REAL :: usflx (nx,ny)        ! Surface flux of u-momentum (kg/(m*s**2))
  REAL :: vsflx (nx,ny)        ! Surface flux of v-momentum (kg/(m*s**2))
  REAL :: ptsflx(nx,ny)        ! Surface heat flux (K*kg/(m**2*s))
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
  PARAMETER (fmtver0='004.10 2nd Binary Data')
  PARAMETER (fmtver1='004.10 2nd Binary Data')
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
  INTEGER :: nstyp1
  CHARACTER (LEN=12) :: label
!  INTEGER :: nxin,nyin,nzin
  INTEGER :: nxin,nyin,nzin,nzsoilin
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
  READ(inch,ERR=110,END=120) fmtverin

  IF( fmtverin /= fmtver0 .AND. fmtverin /= fmtver1 ) THEN
    WRITE(6,'(/1x,a/1x,2a/1x,2a/1x,2a/1x,a)')                           &
        'Data format incompatible with the data reader.',               &
        'Format of data is ',fmtverin,' Format of reader is ',fmtver1,  &
        'compitable to ',fmtver0, '. Job stopped.'
    CALL arpsstop('arpstop called from bn2read header read ',1)
  END IF

  READ(inch,ERR=110,END=120) runname
  READ(inch,ERR=110,END=120) nocmnt
  IF( nocmnt > 0 ) THEN
    DO i=1,nocmnt
      READ(inch,ERR=110,END=120) cmnt(i)
    END DO
  END IF

  WRITE(6,'(//''  THE NAME OF THIS RUN IS:  '',A//)') runname

  IF( nocmnt > 0 ) THEN
    DO i=1,nocmnt
      WRITE(6,'(1x,a)') cmnt(i)
    END DO
  END IF

  READ(inch,ERR=110,END=120) time,tmunit
!
!-----------------------------------------------------------------------
!
!  Get dimensions of data in binary file and check against
!  the dimensions passed to BN2READ
!
!-----------------------------------------------------------------------
!
!  READ(inch,ERR=110,END=120) nxin, nyin, nzin
  READ(inch,ERR=110,END=120) nxin, nyin, nzin, nzsoilin

!  IF( nxin /= nx .OR. nyin /= ny .OR. nzin /= nz ) THEN
  IF( nxin /= nx .OR. nyin /= ny .OR. nzin /= nz .OR. nzsoilin /= nzsoil) THEN
    WRITE(6,'(1x,a)')                                                   &
         ' Dimensions in BN2READ inconsistent with data.'
!    WRITE(6,'(1x,a,3I15)') ' Read were: ', nxin, nyin, nzin
    WRITE(6,'(1x,a,4I15)') ' Read were: ', nxin, nyin, nzin, nzsoilin
    WRITE(6,'(1x,a)')                                                   &
         ' Program aborted in BN2READ.'
!    CALL arpsstop('arpstop called from bn2read nx-ny-nz read',1)
    CALL arpsstop('arpstop called from bn2read nx-ny-nz-nzsoil read',1)
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

    READ(inch,ERR=110,END=120)                                          &
         bgrdin,bbasin,bvarin,mstin,bicein,                             &
         btrbin,idummy,idummy,landin,totin,                             &
         btkein,idummy,idummy,mapproj,month,                            &
         day,year,hour,minute,second

  ELSE ! Normal data reading

    WRITE(lchanl,'(1x,a,f8.1,a,f8.3,a/)')'To read data for time:',      &
         time,' secs = ',(time/60.),' mins.'

    READ(inch,ERR=110,END=120)                                          &
         grdin,basin,varin,mstin,icein,                                 &
         trbin,sfcin,rainin,landin,totin,                               &
         tkein,idummy,idummy,mapproj,month,                             &
         day,year,hour,minute,second

  END IF

  READ(inch,ERR=110,END=120)                                            &
                  umove,vmove,xgrdorg,ygrdorg,trulat1,                  &
                  trulat2,trulon,sclfct,rdummy,rdummy,                  &
                  rdummy,rdummy,rdummy,rdummy,rdummy,                   &
                  tstop,thisdmp,latitud,ctrlat,ctrlon

  IF ( totin /= 0 ) THEN
!
!-----------------------------------------------------------------------
!
!  Read in additional parameters for ARPS history dump 4.1 or later
!  version.
!
!-----------------------------------------------------------------------
!
    READ(inch,ERR=110,END=120)                                          &
         nstyp1,  prcin, radin, flxin,snowcin,                           &
         snowin,idummy,idummy,idummy,idummy,                            &
         idummy,idummy,idummy,idummy,idummy,                            &
         idummy,idummy,idummy,idummy,idummy

    IF ( nstyp1 < 1 ) THEN
      nstyp1 = 1
    END IF

    READ(inch,ERR=110,END=120)                                          &
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
    READ(inch,ERR=110,END=120) label
    READ(inch,ERR=110,END=120) (x(i),i=1,nx)
    WRITE(lchanl,910) label,' x.'

    READ(inch,ERR=110,END=120) label
    READ(inch,ERR=110,END=120) (y(j),j=1,ny)
    WRITE(lchanl,910) label,' y.'

    READ(inch,ERR=110,END=120) label
    READ(inch,ERR=110,END=120) (z(k),k=1,nz)
    WRITE(lchanl,910) label,' z.'

    READ(inch,ERR=110,END=120) label
    READ(inch,ERR=110,END=120)                                          &
        (((zp(i,j,k),i=1,nx),j=1,ny),k=1,nz)
    WRITE(lchanl,910) label,' zp.'

    READ(inch,ERR=110,END=120) label
    READ(inch,ERR=110,END=120)                                          &
        (((zpsoil(i,j,k),i=1,nx),j=1,ny),k=1,nzsoil)
    WRITE(lchanl,910) label,' zpsoil.'

  END IF  ! grdin
!
!-----------------------------------------------------------------------
!
!  Read in base state fields
!
!----------------------------------------------------------------------
!
  IF( basin == 1 .OR. grdbas == 1 ) THEN

    READ(inch,ERR=110,END=120) label
    READ(inch,ERR=110,END=120)                                          &
        (((ubar(i,j,k),i=1,nx),j=1,ny),k=1,nz)
    WRITE(lchanl,910) label,' ubar.'

    READ(inch,ERR=110,END=120) label
    READ(inch,ERR=110,END=120)                                          &
        (((vbar(i,j,k),i=1,nx),j=1,ny),k=1,nz)
    WRITE(lchanl,910) label,' vbar.'

    READ(inch,ERR=110,END=120) label
    READ(inch,ERR=110,END=120)                                          &
        (((wbar(i,j,k),i=1,nx),j=1,ny),k=1,nz)
    WRITE(lchanl,910) label,' wbar.'

    READ(inch,ERR=110,END=120) label
    READ(inch,ERR=110,END=120)                                          &
        (((ptbar(i,j,k),i=1,nx),j=1,ny),k=1,nz)
    WRITE(lchanl,910) label,' ptbar.'

    READ(inch,ERR=110,END=120) label
    READ(inch,ERR=110,END=120)                                          &
        (((pbar(i,j,k),i=1,nx),j=1,ny),k=1,nz)
    WRITE(lchanl,910) label,' pbar.'


    IF( mstin == 1) THEN
      READ(inch,ERR=110,END=120) label
      READ(inch,ERR=110,END=120)                                        &
          (((qvbar(i,j,k),i=1,nx),j=1,ny),k=1,nz)
      WRITE(lchanl,910) label,' qvbar.'
    END IF

    IF(landin == 1) THEN

      IF( nstyp1 <= 1 ) THEN

        READ(inch,ERR=110,END=120) label
        READ(inch,ERR=110,END=120)                                      &
              ((soiltyp(i,j,1),i=1,nx),j=1,ny)
        WRITE(lchanl,910) label,' soiltyp.'

      ELSE

        DO is=1,nstyp1
          IF (is <= nstyps) THEN
            READ(inch,ERR=110,END=120) label
            READ(inch,ERR=110,END=120)                                  &
                  ((soiltyp(i,j,is),i=1,nx),j=1,ny)
            WRITE(lchanl,910) label,' soiltyp.'

            READ(inch,ERR=110,END=120) label
            READ(inch,ERR=110,END=120)                                  &
                  ((stypfrct(i,j,is),i=1,nx),j=1,ny)
            WRITE(lchanl,910) label,' stypfrct.'
          ELSE
            READ(inch,ERR=110,END=120) label
            READ(inch,ERR=110,END=120)
            WRITE(lchanl,910) label,'skipping soiltyp.'

            READ(inch,ERR=110,END=120) label
            READ(inch,ERR=110,END=120)
            WRITE(lchanl,910) label,'skipping stypfrct.'
          ENDIF
        END DO

      END IF

      CALL fix_stypfrct_nstyp(nx,ny,nstyp1,nstyp,stypfrct)

      READ(inch,ERR=110,END=120) label
      READ(inch,ERR=110,END=120) ((vegtyp (i,j),i=1,nx),j=1,ny)
      WRITE(lchanl,910) label,' vegtyp.'

      READ(inch,ERR=110,END=120) label
      READ(inch,ERR=110,END=120) ((lai    (i,j),i=1,nx),j=1,ny)
      WRITE(lchanl,910) label,' lai.'

      READ(inch,ERR=110,END=120) label
      READ(inch,ERR=110,END=120) ((roufns (i,j),i=1,nx),j=1,ny)
      WRITE(lchanl,910) label,' roufns.'

      READ(inch,ERR=110,END=120) label
      READ(inch,ERR=110,END=120) ((veg    (i,j),i=1,nx),j=1,ny)
      WRITE(lchanl,910) label,' veg.'

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
      READ(inch,ERR=110,END=120) label
      READ(inch,ERR=110,END=120)                                        &
          (((uprt(i,j,k),i=1,nx),j=1,ny),k=1,nz)
      WRITE(lchanl,910) label,' uprt.'

      READ(inch,ERR=110,END=120) label
      READ(inch,ERR=110,END=120)                                        &
          (((vprt(i,j,k),i=1,nx),j=1,ny),k=1,nz)
      WRITE(lchanl,910) label,' vprt.'

      READ(inch,ERR=110,END=120) label
      READ(inch,ERR=110,END=120)                                        &
          (((wprt(i,j,k),i=1,nx),j=1,ny),k=1,nz)
      WRITE(lchanl,910) label,' wprt.'
!
!-----------------------------------------------------------------------
!
!  Read in scalars
!
!----------------------------------------------------------------------
!
      READ(inch,ERR=110,END=120) label
      READ(inch,ERR=110,END=120)                                        &
          (((ptprt(i,j,k),i=1,nx),j=1,ny),k=1,nz)
      WRITE(lchanl,910) label,' ptprt.'

      READ(inch,ERR=110,END=120) label
      READ(inch,ERR=110,END=120)                                        &
          (((pprt(i,j,k),i=1,nx),j=1,ny),k=1,nz)
      WRITE(lchanl,910) label,' pprt.'

    ELSE

      READ(inch,ERR=110,END=120) label
      READ(inch,ERR=110,END=120)                                        &
          (((uprt(i,j,k),i=1,nx),j=1,ny),k=1,nz)
      WRITE(lchanl,910) label,' u.'
      DO k=1,nz-1
        DO j=1,ny-1
          DO i=1,nx
            uprt(i,j,k) = uprt(i,j,k) - ubar(i,j,k)
          END DO
        END DO
      END DO

      READ(inch,ERR=110,END=120) label
      READ(inch,ERR=110,END=120)                                        &
          (((vprt(i,j,k),i=1,nx),j=1,ny),k=1,nz)
      WRITE(lchanl,910) label,' v.'
      DO k=1,nz-1
        DO j=1,ny
          DO i=1,nx-1
            vprt(i,j,k) = vprt(i,j,k) - vbar(i,j,k)
          END DO
        END DO
      END DO

      READ(inch,ERR=110,END=120) label
      READ(inch,ERR=110,END=120)                                        &
          (((wprt(i,j,k),i=1,nx),j=1,ny),k=1,nz)
      WRITE(lchanl,910) label,' w.'
      READ(inch,ERR=110,END=120) label
      READ(inch,ERR=110,END=120)                                        &
          (((ptprt(i,j,k),i=1,nx),j=1,ny),k=1,nz)
      WRITE(lchanl,910) label,' pt.'
      DO k=1,nz-1
        DO j=1,ny-1
          DO i=1,nx-1
            ptprt(i,j,k) = ptprt(i,j,k) - ptbar(i,j,k)
          END DO
        END DO
      END DO

      READ(inch,ERR=110,END=120) label
      READ(inch,ERR=110,END=120)                                        &
          (((pprt(i,j,k),i=1,nx),j=1,ny),k=1,nz)
      WRITE(lchanl,910) label,' p.'
      DO k=1,nz-1
        DO j=1,ny-1
          DO i=1,nx-1
            pprt(i,j,k) = pprt(i,j,k) - pbar(i,j,k)
          END DO
        END DO
      END DO

    END IF

  END IF
!
!-----------------------------------------------------------------------
!
!  Read in moisture variables
!
!----------------------------------------------------------------------
!
  IF( mstin == 1 ) THEN

    IF ( totin == 0 ) THEN

      READ(inch,ERR=110,END=120) label
      READ(inch,ERR=110,END=120)                                        &
          (((qvprt(i,j,k),i=1,nx),j=1,ny),k=1,nz)
      WRITE(lchanl,910) label,' qvprt.'

    ELSE

      READ(inch,ERR=110,END=120) label
      READ(inch,ERR=110,END=120)                                        &
          (((qvprt(i,j,k),i=1,nx),j=1,ny),k=1,nz)
      WRITE(lchanl,910) label,' qv.'
      DO k=1,nz-1
        DO j=1,ny-1
          DO i=1,nx-1
            qvprt(i,j,k) = qvprt(i,j,k) - qvbar(i,j,k)
          END DO
        END DO
      END DO

    END IF

    READ(inch,ERR=110,END=120) label
    READ(inch,ERR=110,END=120)                                          &
        (((qc(i,j,k),i=1,nx),j=1,ny),k=1,nz)
    WRITE(lchanl,910) label,' qc.'

    READ(inch,ERR=110,END=120) label
    READ(inch,ERR=110,END=120)                                          &
        (((qr(i,j,k),i=1,nx),j=1,ny),k=1,nz)
    WRITE(lchanl,910) label,' qr.'

    IF( rainin == 1 ) THEN

      READ(inch,ERR=110,END=120) label
      READ(inch,ERR=110,END=120)                                        &
            ((raing(i,j),i=1,nx),j=1,ny)
      WRITE(lchanl,910) label,' raing.'

      READ(inch,ERR=110,END=120) label
      READ(inch,ERR=110,END=120)                                        &
            ((rainc(i,j),i=1,nx),j=1,ny)
      WRITE(lchanl,910) label,' rainc.'

    END IF

    IF( prcin == 1 ) THEN

      READ(inch,ERR=110,END=120) label
      READ(inch,ERR=110,END=120)                                        &
            ((prcrate(i,j,1),i=1,nx),j=1,ny)
      WRITE(lchanl,910) label,' prcrate1.'

      READ(inch,ERR=110,END=120) label
      READ(inch,ERR=110,END=120)                                        &
            ((prcrate(i,j,2),i=1,nx),j=1,ny)
      WRITE(lchanl,910) label,' prcrate2.'

      READ(inch,ERR=110,END=120) label
      READ(inch,ERR=110,END=120)                                        &
            ((prcrate(i,j,3),i=1,nx),j=1,ny)
      WRITE(lchanl,910) label,' prcrate3.'

      READ(inch,ERR=110,END=120) label
      READ(inch,ERR=110,END=120)                                        &
            ((prcrate(i,j,4),i=1,nx),j=1,ny)
      WRITE(lchanl,910) label,' prcrate4.'

    END IF

    IF( icein == 1 ) THEN

      READ(inch,ERR=110,END=120) label
      READ(inch,ERR=110,END=120)                                        &
          (((qi(i,j,k),i=1,nx),j=1,ny),k=1,nz)
      WRITE(lchanl,910) label,' qi.'

      READ(inch,ERR=110,END=120) label
      READ(inch,ERR=110,END=120)                                        &
          (((qs(i,j,k),i=1,nx),j=1,ny),k=1,nz)
      WRITE(lchanl,910) label,' qs.'

      READ(inch,ERR=110,END=120) label
      READ(inch,ERR=110,END=120)                                        &
          (((qh(i,j,k),i=1,nx),j=1,ny),k=1,nz)
      WRITE(lchanl,910) label,' qh.'

    END IF
  END IF

  IF( tkein == 1 ) THEN

    READ(inch,ERR=110,END=120) label
    READ(inch,ERR=110,END=120)                                          &
          (((tke(i,j,k),i=1,nx),j=1,ny),k=1,nz)
    WRITE(lchanl,910) label,' tke.'

  END IF

  IF( trbin == 1 ) THEN

    READ(inch,ERR=110,END=120) label
    READ(inch,ERR=110,END=120)                                          &
          (((kmh(i,j,k),i=1,nx),j=1,ny),k=1,nz)
    WRITE(lchanl,910) label,' kmh.'

    READ(inch,ERR=110,END=120) label
    READ(inch,ERR=110,END=120)                                          &
          (((kmv(i,j,k),i=1,nx),j=1,ny),k=1,nz)
    WRITE(lchanl,910) label,' kmv.'

  END IF

  IF( sfcin == 1 ) THEN

    IF (nstyp1 <= 1) THEN

      READ(inch,ERR=110,END=120) label
      READ(inch,ERR=110,END=120) (((tsoil(i,j,k,0),i=1,nx),j=1,ny),k=1,nzsoil)
      WRITE(lchanl,910) label,' tsoil.'

      READ(inch,ERR=110,END=120) label
      READ(inch,ERR=110,END=120) (((qsoil(i,j,k,0),i=1,nx),j=1,ny),k=1,nzsoil)
      WRITE(lchanl,910) label,' qsoil.'

      READ(inch,ERR=110,END=120) label
      READ(inch,ERR=110,END=120) ((wetcanp(i,j,0),i=1,nx),j=1,ny)
      WRITE(lchanl,910) label,' wetcanp.'

    ELSE

      DO is=0,nstyp1
        IF (is <= nstyps) THEN
          READ(inch,ERR=110,END=120) label
          READ(inch,ERR=110,END=120)(((tsoil(i,j,k,is),i=1,nx),j=1,ny), &
               k=1,nzsoil)
          WRITE(lchanl,910) label,' tsoil.'

          READ(inch,ERR=110,END=120) label
          READ(inch,ERR=110,END=120)(((qsoil(i,j,k,is),i=1,nx),j=1,ny), &
               k=1,nzsoil)
          WRITE(lchanl,910) label,' qsoil.'

          READ(inch,ERR=110,END=120) label
          READ(inch,ERR=110,END=120)((wetcanp(i,j,is),i=1,nx),j=1,ny)
          WRITE(lchanl,910) label,' wetcanp.'
        ELSE
          READ(inch,ERR=110,END=120) label
          READ(inch,ERR=110,END=120)
          WRITE(lchanl,910) label,'skipping tsoil.'

          READ(inch,ERR=110,END=120) label
          READ(inch,ERR=110,END=120)
          WRITE(lchanl,910) label,'skipping qsoil.'

          READ(inch,ERR=110,END=120) label
          READ(inch,ERR=110,END=120)
          WRITE(lchanl,910) label,'skipping wetcanp.'
        ENDIF
      END DO

    END IF

    CALL fix_soil_nstyp(nx,ny,nzsoil,nstyp1,nstyp,tsoil,qsoil,wetcanp)

    IF(snowcin == 1) THEN
      READ(inch,ERR=110,END=120) label
      READ(inch,ERR=110,END=120)
      WRITE(lchanl,910) label,' snowcvr -- discarding.'
    END IF

    IF(snowin == 1) THEN
      READ(inch,ERR=110,END=120) label
      READ(inch,ERR=110,END=120)((snowdpth(i,j),i=1,nx),j=1,ny)
      WRITE(lchanl,910) label,' snowdpth.'
    END IF

  END IF

  IF( radin == 1 ) THEN

    READ(inch,ERR=110,END=120) label
    READ(inch,ERR=110,END=120)                                          &
          (((radfrc(i,j,k),i=1,nx),j=1,ny),k=1,nz)
    WRITE(lchanl,910) label,' radfrc.'

    READ(inch,ERR=110,END=120) label
    READ(inch,ERR=110,END=120) ((radsw(i,j),i=1,nx),j=1,ny)
    WRITE(lchanl,910) label,' radsw.'

    READ(inch,ERR=110,END=120) label
    READ(inch,ERR=110,END=120) ((rnflx(i,j),i=1,nx),j=1,ny)
    WRITE(lchanl,910) label,' rnflx.'

    READ(inch,ERR=110,END=120) label
    READ(inch,ERR=110,END=120) ((radswnet(i,j),i=1,nx),j=1,ny)
    WRITE(lchanl,910) label,' radswnet.'

    READ(inch,ERR=110,END=120) label
    READ(inch,ERR=110,END=120) ((radlwin(i,j),i=1,nx),j=1,ny)
    WRITE(lchanl,910) label,' radlwin.'


  END IF

  IF( flxin == 1 ) THEN

    READ(inch,ERR=110,END=120) label
    READ(inch,ERR=110,END=120) ((usflx(i,j),i=1,nx),j=1,ny)
    WRITE(lchanl,910) label,' usflx.'

    READ(inch,ERR=110,END=120) label
    READ(inch,ERR=110,END=120) ((vsflx(i,j),i=1,nx),j=1,ny)
    WRITE(lchanl,910) label,' vsflx.'

    READ(inch,ERR=110,END=120) label
    READ(inch,ERR=110,END=120) ((ptsflx(i,j),i=1,nx),j=1,ny)
    WRITE(lchanl,910) label,' ptsflx.'

    READ(inch,ERR=110,END=120) label
    READ(inch,ERR=110,END=120) ((qvsflx(i,j),i=1,nx),j=1,ny)
    WRITE(lchanl,910) label,' qvsflx.'

  END IF

  910   FORMAT(1X,'Field ',a12,' was read into array',a)

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
!
!-----------------------------------------------------------------------
!
!  Error during read
!
!----------------------------------------------------------------------
!

  110   CONTINUE
  WRITE(6,'(/a/)') ' Error reading data in BN2READ'
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
  WRITE(6,'(/a/)') ' End of file reached in BN2READ'
  ireturn=2
  RETURN
END SUBROUTINE bn2read
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE BINDUMP                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE bindump(nx,ny,nz,nzsoil,nstyps, nchanl, grdbas,              &
           u,v,w,ptprt,pprt,qv,qscalar,tke,kmh,kmv,                     &
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
!  Write history data into channel nchanl as binary data.
!
!  All data read in are located at the original staggered grid points.
!
!  Note: coordinate fields are dumped as 3 dimensional fields which
!  have been converted from meters to kilometers.  This is for the
!  convenience of the plotting applications.
!
!  The last 4 characters of the 12 character label written out
!  with each 1-,2-, or 3-d array is used by the splitdump and
!  joinfiles subroutines used by the message passing version of the
!  ARPS (and also by some auxiliary ARPS I/O routines)
!  to determine the data type of the array.
!  Key to the labels:
!
!    'nnnnnnn tdds'
!
!     n  - characters containing the name of the variable.
!     t  - type of variable:  "r" for real and "i" for integer.
!     dd - number of dimensions:  "1d" "2d" or "3d".
!     s  - staggered dimension: "0" for centered,
!                               "1" for staggered in x,
!                               "2" for staggered in y,
!                               "3" for staggered in z.
!
!-----------------------------------------------------------------------
!
!
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  3/10/92.
!
!  MODIFICATION HISTORY:
!
!  6/06/92 (M. Xue)
!  Added full documentation.
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
!  Added map projection parameters into the binary dumping
!
!  03/26/96 (G. Bassett)
!  Labels were modified to include information about array type.
!  This information is used by splitdump and joinfiles subroutines.
!
!  12/09/1998 (Donghai Wang)
!  Added the snow cover.
!
!  05/15/2002 (J. Brotzge)
!  Added to allow for multiple soil schemes.
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
!    pprt     Perturbation pressure at a given time level (Pascal)
!    qv       Water vapor specific humidity at a given time level (kg/kg)
!    qscalar
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
!    ptsflx   Surface heat flux (K*kg/(m**2*s))
!    qvsflx   Surface moisture flux (kg/(m**2*s))
!
!  OUTPUT:
!
!    None.
!
!  WORK ARRAY:
!
!    tem1     Temporary work array.
!
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
!  tkeout =0 or 1. If tkeout=0, tke is not dumped.
!  trbout =0 or 1. If trbout=0, turbulence parameter km is not dumped.
!  sfcout =0 or 1. If sfcout=0, surface variables are not dumped.
!  landout=0 or 1. If landout=0, surface propertty arrays are not dumped.
!  radout =0 or 1. If radout =0, radiation arrays are not dumped.
!  flxout =0 or 1. If flxout =0, surface flux arrays are not dumped.
!
!  These following parameters are also passed in through common
!  blocks in globcst.inc.
!
!  runname,curtim,umove,vmove,xgrdorg,ygrdorg
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
  INCLUDE 'grid.inc'           ! Grid & map parameters.
  INCLUDE 'phycst.inc'
  INCLUDE 'mp.inc'             ! mpi parameters.

!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
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

  REAL :: qscalar(nx,ny,nz,nscalar)

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
  REAL :: zpsoil(nx,ny,nzsoil) ! The physical height coordinate defined at
                               ! w-point of the soil

  INTEGER :: nstyps                  ! Number of soil types
  INTEGER :: soiltyp(nx,ny,nstyps)   ! Soil type
  REAL    :: stypfrct(nx,ny,nstyps)  ! Soil type fractions
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
  REAL :: ptsflx(nx,ny)        ! Surface heat flux (K*kg/(m**2*s))
  REAL :: qvsflx(nx,ny)        ! Surface moisture flux (kg/(m**2*s))

  REAL :: tem1  (nx,ny,nz)     ! Temporary work array
!
!-----------------------------------------------------------------------
!
!  Parameters describing this routine
!
!-----------------------------------------------------------------------
!
! 06/28/2002 Zuwen He
!
! fmtver??: to label each data a version.
! intver??: an integer to allow faster comparison than fmtver??,
!           which are strings.
!
! Version 5.00: significant change in soil variables since version 4.10.
! Version 5.30: significant changes in microphysics variables
!
  CHARACTER (LEN=40) :: fmtver,fmtver410,fmtver500,fmtver530
  INTEGER            :: intver,intver410,intver500,intver530
  PARAMETER (fmtver410='004.10 Binary Data',intver410=410)
  PARAMETER (fmtver500='005.00 Binary Data',intver500=500)
  PARAMETER (fmtver530='005.30 Binary Data',intver530=530)

  CHARACTER (LEN=10) :: tmunit
  PARAMETER (tmunit='seconds   ')
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k,l,is,nq
  INTEGER :: idummy

  REAL(4) :: rdummy
  REAL(4) :: umove4,   vmove4,   xgrdorg4, ygrdorg4
  REAL(4) :: trulat14, trulat24, trulon4,  sclfct4
  REAL(4) :: tstop4,   thisdmp4, latitud4
  REAL(4) :: ctrlat4,  ctrlon4
  REAL(4) :: n0rain4, n0snow4, n0hail4, rhosnow4, rhohail4
  REAL(4) :: ntcloud4,n0grpl4,rhoice4,rhogrpl4,alpharain4
  REAL(4) :: alphaice4,alphasnow4,alphagrpl4,alphahail4

  REAL(4), ALLOCATABLE :: out1d(:)
  REAL(4), ALLOCATABLE :: out2d(:,:)
  REAL(4), ALLOCATABLE :: out3d(:,:,:)
  REAL(4), ALLOCATABLE :: out3dsoil(:,:,:)

!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!---------------------------------------------------------------
!
! Allocate 4 byte arrays
!
!---------------------------------------------------------------

  ALLOCATE(out1d    (MAX(nx,ny,nz)), STAT = idummy)
  ALLOCATE(out2d    (nx,ny),         STAT = idummy)
  ALLOCATE(out3d    (nx,ny,nz),      STAT = idummy)
  ALLOCATE(out3dsoil(nx,ny,nzsoil),  STAT = idummy)


  intver = intver530  !  for the time being, in the future, we will
                      !  allow to dump data in the different version
                      !  intver will be assigned from input file
  IF (intver == intver410) THEN
    fmtver=fmtver410
  ELSE IF (intver == intver500) THEN
    fmtver=fmtver500
  ELSE IF (intver == intver530) THEN
    fmtver=fmtver530
  ELSE
    IF (myproc == 0) WRITE(6,'(/1x,a,i10,a/)')                          &
        'Data format, intver=',intver,', not found. The Job stopped.'
    CALL arpsstop('arpstop called from bindump.',1)
  END IF

  IF (myproc == 0) WRITE(6,'(/1x,a,a,a/)')                              &
                            'Data dump format, fmtver=',fmtver,'. '

  IF (myproc == 0)  &
     WRITE(6,'(1x,a,f13.3/)') 'Writing history data at time=', curtim
!
!-----------------------------------------------------------------------
!
!  Write header info
!
!-----------------------------------------------------------------------
!
  WRITE(nchanl) fmtver
  WRITE(nchanl) runname
  WRITE(nchanl) nocmnt
  IF( nocmnt > 0 ) THEN
    DO l=1,nocmnt
      WRITE(nchanl) cmnt(l)
    END DO
  END IF

  rdummy = curtim
  WRITE(nchanl) rdummy,tmunit

  IF (intver == intver410) THEN
    WRITE(nchanl) nx,ny,nz
  ELSE IF (intver >= intver500) THEN
    WRITE(nchanl) nx,ny,nz,nzsoil
  END IF
!
!-----------------------------------------------------------------------
!
!  Write the flags for different data groups.
!
!-----------------------------------------------------------------------
!
  idummy = 0

  IF( grdbas == 1 ) THEN

    WRITE(nchanl)      1,      1,      0, mstout,      0,               &
                       0,      0,      0, landout,totout,               &
                       0, idummy, idummy, mapproj, month,               &
                     day,   year,   hour, minute, second

  ELSE

    WRITE(nchanl) grdout, basout, varout, mstout, iceout,               &
                  trbout, sfcout, rainout,landout,totout,               &
                  tkeout,nscalar, idummy, mapproj, month,               &
                     day,   year,   hour, minute, second

    WRITE(nchanl)   p_qc,  p_qr,  p_qi,  p_qs,  p_qg, p_qh,             &
                    p_nc,  p_nr,  p_ni,  p_ns,  p_ng, p_nh,             &
                    p_zr,  p_zi,  p_zs,  p_zg,  p_zh,idummy,            &
                  idummy,idummy,idummy,idummy,idummy,idummy,            &
                  idummy,idummy,idummy,idummy,  p_cc,mphyopt

  END IF

  rdummy   = 0.0
  umove4   = umove
  vmove4   = vmove
  xgrdorg4 = xgrdorg
  ygrdorg4 = ygrdorg
  trulat14 = trulat1
  trulat24 = trulat2
  trulon4  = trulon
  sclfct4  = sclfct
  tstop4   = tstop
  thisdmp4 = thisdmp
  latitud4 = latitud
  ctrlat4  = ctrlat
  ctrlon4  = ctrlon
  ntcloud4 = ntcloud
  n0rain4  = n0rain
  n0snow4  = n0snow
  n0grpl4  = n0grpl
  n0hail4  = n0hail
  rhoice4 = rhoice
  rhosnow4 = rhosnow
  rhogrpl4 = rhogrpl
  rhohail4 = rhohail
  alpharain4 = alpharain
  alphaice4 = alphaice
  alphasnow4 = alphasnow
  alphagrpl4 = alphagrpl
  alphahail4 = alphahail

  WRITE(nchanl)   umove4,     vmove4,   xgrdorg4,  ygrdorg4,   trulat14, &
                trulat24,    trulon4,    sclfct4,   ntcloud,    n0rain4, &
                 n0snow4,    n0grpl4,    n0hail4,   rhoice4,   rhosnow4, &
                rhogrpl4,   rhohail4, alpharain4, alphaice4, alphasnow4, &
              alphagrpl4, alphahail4,     rdummy,    rdummy,     rdummy, &
                  tstop4,   thisdmp4,   latitud4,   ctrlat4,    ctrlon4
!
!-----------------------------------------------------------------------
!
!  If totout=1, write additional parameters to history dump files.
!  This is for ARPS version 4.1.2 or later.
!
!-----------------------------------------------------------------------
!
  IF ( totout == 1 ) THEN
    WRITE(nchanl) nstyp,  prcout, radout, flxout,      0,  & ! 0 for snowcvr
              snowout,idummy, idummy, idummy, idummy,                   &
                  idummy, idummy, idummy, idummy, idummy,               &
                  idummy, idummy, idummy, idummy, idummy

    WRITE(nchanl) rdummy, rdummy, rdummy, rdummy, rdummy,               &
                  rdummy, rdummy, rdummy, rdummy, rdummy,               &
                  rdummy, rdummy, rdummy, rdummy, rdummy,               &
                  rdummy, rdummy, rdummy, rdummy, rdummy
  END IF
!
!-----------------------------------------------------------------------
!
!  If grdout=1 or grdbas=1, write out grid variables
!
!-----------------------------------------------------------------------
!
  IF(grdout == 1 .OR. grdbas == 1 ) THEN

    out1d(1:nx) = x
    WRITE(nchanl) 'x coord r1d1'
    WRITE(nchanl) out1d(1:nx)

    out1d(1:ny) = y
    WRITE(nchanl) 'y coord r1d2'
    WRITE(nchanl) out1d(1:ny)

    out1d(1:nz) = z
    WRITE(nchanl) 'z coord r1d3'
    WRITE(nchanl) out1d(1:nz)

    CALL edgfill(zp,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz)
    out3d(:,:,:) = zp(:,:,:)
    WRITE(nchanl) 'zp coor r3d0'
    WRITE(nchanl) out3d

    IF (intver >= intver500) THEN
      CALL edgfill(zpsoil,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nzsoil,1,nzsoil)
      out3dsoil(:,:,:) = zpsoil(:,:,:)
      WRITE(nchanl) 'zpsoil rs3d0'
      WRITE(nchanl) out3dsoil
    END IF

  END IF    ! grdout
!
!-----------------------------------------------------------------------
!
!  If basout=1, write out base state variables.
!
!-----------------------------------------------------------------------
!
  IF(basout == 1 .OR. grdbas == 1 ) THEN

    CALL edgfill(ubar,1,nx,1,nx, 1,ny,1,ny-1, 1,nz,1,nz-1)
    out3d(:,:,:) = ubar(:,:,:)
    WRITE(nchanl) 'ubar    r3d1'
    WRITE(nchanl) out3d

    CALL edgfill(vbar,1,nx,1,nx-1, 1,ny,1,ny, 1,nz,1,nz-1)
    WRITE(nchanl) 'vbar    r3d2'
    out3d(:,:,:) = vbar(:,:,:)
    WRITE(nchanl) out3d

    DO k=1,nz
      DO j=1,ny
        DO i=1,nx
          out3d(i,j,k) = 0.0
        END DO
      END DO
    END DO
    WRITE(nchanl) 'wbar    r3d3'
    WRITE(nchanl) out3d

    CALL edgfill(ptbar,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)
    out3d(:,:,:) = ptbar(:,:,:)
    WRITE(nchanl) 'ptbar   r3d0'
    WRITE(nchanl)  out3d

    CALL edgfill(pbar,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)
    out3d(:,:,:) = pbar(:,:,:)
    WRITE(nchanl) 'pbar    r3d0'
    WRITE(nchanl)  out3d

    IF(mstout == 1) THEN

      CALL edgfill(qvbar,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)
      out3d(:,:,:) = qvbar(:,:,:)
      WRITE(nchanl) 'qvbar   r3d0'
      WRITE(nchanl)  out3d

    END IF

    IF(landout == 1) THEN

      IF( nstyp <= 1 ) THEN

        CALL iedgfill(soiltyp(1,1,1),1,nx,1,nx-1, 1,ny,1,ny-1,          &
                      1,1,1,1)
        WRITE(nchanl) 'soiltyp i2d0'
        WRITE(nchanl) ((soiltyp(i,j,1),i=1,nx),j=1,ny)

      ELSE
        DO is=1,nstyp

          CALL iedgfill(soiltyp(1,1,is),1,nx,1,nx-1, 1,ny,1,ny-1,       &
                        1,1,1,1)
          WRITE(nchanl) 'soiltyp i2d0'
          WRITE(nchanl) ((soiltyp(i,j,is),i=1,nx),j=1,ny)

          CALL edgfill(stypfrct(1,1,is),1,nx,1,nx-1, 1,ny,1,ny-1,       &
                       1,1,1,1)
          out2d(:,:) = stypfrct(:,:,is)
          WRITE(nchanl) 'stypfrc r2d0'
          WRITE(nchanl)  out2d

        END DO
      END IF

      CALL iedgfill(vegtyp ,1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)
      WRITE(nchanl) 'vegtyp  i2d0'
      WRITE(nchanl)  vegtyp

      CALL edgfill(lai    ,1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)
      WRITE(nchanl) 'lai     r2d0'
      out2d(:,:) = lai(:,:)
      WRITE(nchanl)  out2d

      CALL edgfill(roufns ,1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)
      WRITE(nchanl) 'roufns  r2d0'
      out2d(:,:) = roufns(:,:)
      WRITE(nchanl)  out2d

      CALL edgfill(veg    ,1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)
      WRITE(nchanl) 'veg     r2d0'
      out2d(:,:) = veg(:,:)
      WRITE(nchanl)  out2d

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
!
!-----------------------------------------------------------------------
!
!  Write out u,v and w.
!
!-----------------------------------------------------------------------
!
  IF(varout == 1) THEN

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
      out3d(:,:,:) = tem1(:,:,:)
      WRITE(nchanl) 'uprt    r3d1'
      WRITE(nchanl) out3d

      DO k=1,nz-1
        DO i=1,nx-1
          DO j=1,ny
            tem1(i,j,k)=v(i,j,k)-vbar(i,j,k)
          END DO
        END DO
      END DO
      CALL edgfill(tem1,1,nx,1,nx-1, 1,ny,1,ny, 1,nz,1,nz-1)
      out3d(:,:,:) = tem1(:,:,:)
      WRITE(nchanl) 'vprt    r3d2'
      WRITE(nchanl) out3d

      CALL edgfill(w,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz)
      out3d(:,:,:) = w(:,:,:)
      WRITE(nchanl) 'wprt    r3d3'
      WRITE(nchanl) out3d
!
!-----------------------------------------------------------------------
!
!  Write out scalars
!
!-----------------------------------------------------------------------
!
      CALL edgfill(ptprt,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)
      out3d(:,:,:) = ptprt(:,:,:)
      WRITE(nchanl) 'ptprt   r3d0'
      WRITE(nchanl) out3d

      CALL edgfill(pprt,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)
      out3d(:,:,:) = pprt(:,:,:)
      WRITE(nchanl) 'pprt    r3d0'
      WRITE(nchanl) out3d

    ELSE
!
!-----------------------------------------------------------------------
!
!  Write out total values to history dump
!
!-----------------------------------------------------------------------
!
      CALL edgfill(u,1,nx,1,nx, 1,ny,1,ny-1, 1,nz,1,nz-1)
      out3d(:,:,:) = u(:,:,:)
      WRITE(nchanl) 'u       r3d1'
      WRITE(nchanl) out3d

      CALL edgfill(v,1,nx,1,nx-1, 1,ny,1,ny, 1,nz,1,nz-1)
      out3d(:,:,:) = v(:,:,:)
      WRITE(nchanl) 'v       r3d2'
      WRITE(nchanl) out3d

      CALL edgfill(w,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz)
      out3d(:,:,:) = w(:,:,:)
      WRITE(nchanl) 'w       r3d3'
      WRITE(nchanl) out3d
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
      out3d(:,:,:) = tem1(:,:,:)
      WRITE(nchanl) 'pt      r3d0'
      WRITE(nchanl) out3d

      DO k=1,nz-1
        DO j=1,ny-1
          DO i=1,nx-1
            tem1(i,j,k) = pbar(i,j,k) + pprt(i,j,k)
          END DO
        END DO
      END DO
      CALL edgfill(tem1,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)
      out3d(:,:,:) = tem1(:,:,:)
      WRITE(nchanl) 'p       r3d0'
      WRITE(nchanl) out3d

    END IF

  END IF     ! varout
!
!-----------------------------------------------------------------------
!
!  If mstout = 1, write out moisture scalars.
!
!-----------------------------------------------------------------------
!
  IF(mstout == 1) THEN

    IF( totout == 0 ) THEN
!
!-----------------------------------------------------------------------
!
!  Write out perturbation to history dump
!
!-----------------------------------------------------------------------
!
      DO k=1,nz-1
        DO j=1,ny-1
          DO i=1,nx-1
            tem1(i,j,k)=qv(i,j,k)-qvbar(i,j,k)
          END DO
        END DO
      END DO

      CALL edgfill(tem1,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)
      out3d(:,:,:) = tem1(:,:,:)
      WRITE(nchanl) 'qvprt   r3d0'
      WRITE(nchanl)  out3d

    ELSE
!
!-----------------------------------------------------------------------
!
!  Write out total values to history dump
!
!-----------------------------------------------------------------------
!
      CALL edgfill(qv,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)
      out3d(:,:,:) = qv(:,:,:)
      WRITE(nchanl) 'qv      r3d0'
      WRITE(nchanl) out3d

    END IF

    DO nq = 1,nscalar
      CALL edgfill(qscalar(:,:,:,nq),1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)
      out3d(:,:,:) = qscalar(:,:,:,nq)
      WRITE(nchanl) qnames(nq)//'    r3d0'
      WRITE(nchanl) out3d
    END DO

    IF(rainout == 1) THEN

      CALL edgfill(raing,   1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)
      WRITE(nchanl) 'raing   r2d0'
      out2d(:,:) = raing(:,:)
      WRITE(nchanl)  out2d

      CALL edgfill(rainc,1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)
      WRITE(nchanl) 'rainc   r2d0'
      out2d(:,:) = rainc(:,:)
      WRITE(nchanl)  out2d

    END IF   !rainout

    IF ( prcout == 1 ) THEN
      CALL edgfill(prcrate,1,nx,1,nx-1, 1,ny,1,ny-1, 1,4,1,4)

      out2d(:,:) = prcrate(:,:,1)
      WRITE(nchanl) 'prcrat1 r2d0'
      WRITE(nchanl)  out2d

      out2d(:,:) = prcrate(:,:,2)
      WRITE(nchanl) 'prcrat2 r2d0'
      WRITE(nchanl)  out2d

      out2d(:,:) = prcrate(:,:,3)
      WRITE(nchanl) 'prcrat3 r2d0'
      WRITE(nchanl)  out2d

      out2d(:,:) = prcrate(:,:,4)
      WRITE(nchanl) 'prcrat4 r2d0'
      WRITE(nchanl)  out2d
    END IF

  END IF   !mstout
!
!-----------------------------------------------------------------------
!
!  If tkeout = 1, write out tke.
!
!-----------------------------------------------------------------------
!
  IF( tkeout == 1 ) THEN

    CALL edgfill(tke,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)
    out3d(:,:,:) = tke(:,:,:)
    WRITE(nchanl) 'tke     r3d0'
    WRITE(nchanl)  out3d

  END IF
!
!-----------------------------------------------------------------------
!
!  If trbout = 1, write out the turbulence parameter, km.
!
!-----------------------------------------------------------------------
!
  IF( trbout == 1 ) THEN

    CALL edgfill(kmh,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)
    out3d(:,:,:) = kmh(:,:,:)
    WRITE(nchanl) 'kmh     r3d0'
    WRITE(nchanl)  out3d

    CALL edgfill(kmv,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)
    out3d(:,:,:) = kmv(:,:,:)
    WRITE(nchanl) 'kmv     r3d0'
    WRITE(nchanl)  out3d

  END IF   ! trbout

!
!-----------------------------------------------------------------------
!
!  If sfcout = 1, write out the surface variables,
!  tsoil, qsoil, and wetcanp.
!
!-----------------------------------------------------------------------
!
  IF( sfcout == 1) THEN

    IF( nstyp <= 1 ) THEN

      CALL edgfill(tsoil(1,1,1,0),  1,nx,1,nx-1, 1,ny,1,ny-1,             &
                   1,nzsoil,1,nzsoil)
      CALL edgfill(qsoil(1,1,1,0), 1,nx,1,nx-1, 1,ny,1,ny-1,             &
                   1,nzsoil,1,nzsoil)

      IF (intver == intver410) THEN
!
! 06/28/2002 Zuwen He
!
! In version 410, tsoil is the average of tsoil from 2 to nzsoil in later
! version, and wetdp is similar.
!
        out3d(:,:,1)=0.

        DO k=2,nzsoil
          DO j=1,ny
            DO i=1,nx
              out3d(i,j,1)=out3d(i,j,1)+tsoil(i,j,k,0)
            END DO
          END DO
        END DO

        DO j=1,ny
          DO i=1,nx
            out3d(i,j,1)=out3d(i,j,1)/float(nzsoil-1)
          END DO
        END DO

        out2d(:,:) = tsoil(:,:,1,0)
        WRITE(nchanl) 'tsfc   r2d0'
        WRITE(nchanl) out2d
        WRITE(nchanl) 'tsoil  r2d0'
        WRITE(nchanl) ((out3d(i,j,1),i=1,nx),j=1,ny)

        out3d(:,:,1) = 0.
        DO k=2,nzsoil
          DO j=1,ny
            DO i=1,nx
              out3d(i,j,1)=out3d(i,j,1)+qsoil(i,j,k,0)
            END DO
          END DO
        END DO

        DO j=1,ny
          DO i=1,nx
            out3d(i,j,1)=out3d(i,j,1)/float(nzsoil-1)
          END DO
        END DO

        out2d(:,:) = qsoil(:,:,1,0)
        WRITE(nchanl) 'wetsfc r2d0'
        WRITE(nchanl) out2d
        WRITE(nchanl) 'wetdp  r2d0'
        WRITE(nchanl) ((out3d(i,j,1),i=1,nx),j=1,ny)

      ELSE IF (intver >= intver500) THEN

        out3dsoil(:,:,:) = tsoil(:,:,:,0)
        WRITE(nchanl) 'tsoil  rs3d0'
        WRITE(nchanl)  out3dsoil
        out3dsoil(:,:,:) = qsoil(:,:,:,0)
        WRITE(nchanl) 'qsoil  rs3d0'
        WRITE(nchanl)  out3dsoil

      END IF  ! intver

      CALL edgfill(wetcanp(1,1,0),1,nx,1,nx-1, 1,ny,1,ny-1,             &
                   1,1,1,1)
      out2d(:,:) = wetcanp(:,:,0)
      WRITE(nchanl) 'wetcanp r2d0'
      WRITE(nchanl)  out2d

    ELSE

      DO is=0,nstyp

        CALL edgfill(tsoil(1,1,1,is),  1,nx,1,nx-1, 1,ny,1,ny-1,        &
                     1,nzsoil,1,nzsoil)
        CALL edgfill(qsoil(1,1,1,is), 1,nx,1,nx-1, 1,ny,1,ny-1,         &
                     1,nzsoil,1,nzsoil)

      IF (intver == intver410) THEN
!
! 06/28/2002 Zuwen He
!
! In version 410, tsoil is the average of tsoil from 2 to nzsoil in later
! version, and wetdp is similar.
!
        out3d(:,:,1)=0.

        DO k=2,nzsoil
          DO j=1,ny
            DO i=1,nx
              out3d(i,j,1)=out3d(i,j,1)+tsoil(i,j,k,is)
            END DO
          END DO
        END DO

        DO j=1,ny
          DO i=1,nx
            out3d(i,j,1)=out3d(i,j,1)/float(nzsoil-1)
          END DO
        END DO

        out2d(:,:) = tsoil(:,:,1,is)
        WRITE(nchanl) 'tsfc  r2d0'
        WRITE(nchanl)  out2d
        WRITE(nchanl) 'tsoil  r2d0'
        WRITE(nchanl) ((out3d(i,j,1),i=1,nx),j=1,ny)

        out3d(:,:,1)=0.
        DO k=2,nzsoil
          DO j=1,ny
            DO i=1,nx
              out3d(i,j,1)=out3d(i,j,1)+qsoil(i,j,k,is)
            END DO
          END DO
        END DO

        DO j=1,ny
          DO i=1,nx
            out3d(i,j,1)=out3d(i,j,1)/float(nzsoil-1)
          END DO
        END DO

        out2d(:,:) = qsoil(:,:,1,is)
        WRITE(nchanl) 'wetsfc r2d0'
        WRITE(nchanl)  out2d
        WRITE(nchanl) 'wetdp  r2d0'
        WRITE(nchanl) ((out3d(i,j,1),i=1,nx),j=1,ny)

      ELSE IF (intver >= intver500) THEN

        out3dsoil(:,:,:) = tsoil(:,:,:,is)
        WRITE(nchanl) 'tsoil  rs3d0'
        WRITE(nchanl)  out3dsoil
        out3dsoil(:,:,:) = qsoil(:,:,:,is)
        WRITE(nchanl) 'qsoil  rs3d0'
        WRITE(nchanl)  out3dsoil

      END IF  ! intver

      CALL edgfill(wetcanp(1,1,is),1,nx,1,nx-1, 1,ny,1,ny-1,          &
                     1,1,1,1)
      out2d(:,:) = wetcanp(:,:,is)
      WRITE(nchanl) 'wetcanp r2d0'
      WRITE(nchanl)  out2d

      END DO
    END IF

    IF (snowout == 1) THEN

      CALL edgfill(snowdpth,1,nx,1,nx-1, 1,ny,1,ny-1,                   &
                    1,1,1,1)
      out2d(:,:) = snowdpth(:,:)
      WRITE(nchanl) 'snowdpthr2d0'
      WRITE(nchanl)  out2d
    END IF

  END IF
!
!-----------------------------------------------------------------------
!
!  If radout = 1, write out the radiation arrays
!
!-----------------------------------------------------------------------
!
  IF( radout == 1 ) THEN

    CALL edgfill(radfrc,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)
    out3d(:,:,:) = radfrc(:,:,:)
    WRITE(nchanl) 'radfrc  r3d0'
    WRITE(nchanl) out3d

    CALL edgfill(radsw,1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)
    out2d(:,:) = radsw(:,:)
    WRITE(nchanl) 'radsw   r2d0'
    WRITE(nchanl) out2d

    CALL edgfill(rnflx,1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)
    out2d(:,:) = rnflx(:,:)
    WRITE(nchanl) 'rnflx   r2d0'
    WRITE(nchanl) out2d

    IF (intver >= intver500) THEN

      CALL edgfill(radswnet,1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)
      out2d(:,:) = radswnet(:,:)
      WRITE(nchanl) 'radswnetr2d0'
      WRITE(nchanl) out2d

      CALL edgfill(radlwin,1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)
      out2d(:,:) = radlwin(:,:)
      WRITE(nchanl) 'radlwin r2d0'
      WRITE(nchanl) out2d

    END IF  ! intver

  END IF   ! radout
!
!-----------------------------------------------------------------------
!
!  If flxout = 1, write out the surface fluxes
!
!-----------------------------------------------------------------------
!
  IF( flxout == 1 ) THEN

    CALL edgfill(usflx,1,nx,1,nx, 1,ny,1,ny-1, 1,1,1,1)
    out2d(:,:) = usflx(:,:)
    WRITE(nchanl) 'usflx   r2d0'
    WRITE(nchanl) out2d

    CALL edgfill(vsflx,1,nx,1,nx-1, 1,ny,1,ny, 1,1,1,1)
    out2d(:,:) = vsflx(:,:)
    WRITE(nchanl) 'vsflx   r2d0'
    WRITE(nchanl) out2d

    CALL edgfill(ptsflx,1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)
    out2d(:,:) = ptsflx(:,:)
    WRITE(nchanl) 'ptsflx  r2d0'
    WRITE(nchanl) out2d

    CALL edgfill(qvsflx,1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)
    out2d(:,:) = qvsflx(:,:)
    WRITE(nchanl) 'qvsflx  r2d0'
    WRITE(nchanl) out2d

  END IF   ! flxout
!---------------------------------------------------------------
!
! Deallocate 4 byte working arrays
!
!---------------------------------------------------------------

  DEALLOCATE(out1d)
  DEALLOCATE(out2d)
  DEALLOCATE(out3d)
  DEALLOCATE(out3dsoil)

  RETURN
END SUBROUTINE bindump
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE BINJOINDUMP                ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE binjoindump(nx,ny,nz,nzsoil,nstyps, nchanl, grdbas,         &
           u,v,w,ptprt,pprt,qv,qscalar,tke,kmh,kmv,                    &
           ubar,vbar,ptbar,pbar,rhobar,qvbar,                          &
           x,y,z,zp,zpsoil,                                            &
           soiltyp,stypfrct,vegtyp,lai,roufns,veg,                     &
           tsoil,qsoil,wetcanp,snowdpth,                               &
           raing,rainc,prcrate,                                        &
           radfrc,radsw,rnflx,radswnet,radlwin,                        &
           usflx,vsflx,ptsflx,qvsflx,                                  &
           tem1)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Write history data into channel nchanl as binary data.
!
!  All data read in are located at the original staggered grid points.
!
!  Note: coordinate fields are dumped as 3 dimensional fields which
!  have been converted from meters to kilometers.  This is for the
!  convenience of the plotting applications.
!
!  The last 4 characters of the 12 character label written out
!  with each 1-,2-, or 3-d array is used by the splitdump and
!  joinfiles subroutines used by the message passing version of the
!  ARPS (and also by some auxiliary ARPS I/O routines)
!  to determine the data type of the array.
!  Key to the labels:
!
!    'nnnnnnn tdds'
!
!     n  - characters containing the name of the variable.
!     t  - type of variable:  "r" for real and "i" for integer.
!     dd - number of dimensions:  "1d" "2d" or "3d".
!     s  - staggered dimension: "0" for centered,
!                               "1" for staggered in x,
!                               "2" for staggered in y,
!                               "3" for staggered in z.
!
!-----------------------------------------------------------------------
!
!
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yunheng Wang
!  8/16/02.
!  Based on subroutine bindump.
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
!    pprt     Perturbation pressure at a given time level (Pascal)
!    qv       Water vapor specific humidity at a given time level (kg/kg)
!    qscalar
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
!    ptsflx   Surface heat flux (K*kg/(m**2*s))
!    qvsflx   Surface moisture flux (kg/(m**2*s))
!
!  OUTPUT:
!
!    None.
!
!  WORK ARRAY:
!
!    tem1     Temporary work array.
!
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
!  tkeout =0 or 1. If tkeout=0, tke is not dumped.
!  trbout =0 or 1. If trbout=0, turbulence parameter km is not dumped.
!  sfcout =0 or 1. If sfcout=0, surface variables are not dumped.
!  landout=0 or 1. If landout=0, surface propertty arrays are not dumped.
!  radout =0 or 1. If radout =0, radiation arrays are not dumped.
!  flxout =0 or 1. If flxout =0, surface flux arrays are not dumped.
!
!  These following parameters are also passed in through common
!  blocks in globcst.inc.
!
!  runname,curtim,umove,vmove,xgrdorg,ygrdorg
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
  INCLUDE 'grid.inc'           ! Grid & map parameters.
  INCLUDE 'phycst.inc'
  INCLUDE 'mp.inc'             ! mpi parameters.

!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
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

  REAL :: qscalar(nx,ny,nz,nscalar)

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
  REAL :: zpsoil(nx,ny,nzsoil) ! The physical height coordinate defined at
                               ! w-point of the soil

  INTEGER :: nstyps                  ! Number of soil types
  INTEGER :: soiltyp(nx,ny,nstyps)   ! Soil type
  REAL    :: stypfrct(nx,ny,nstyps)  ! Soil type fractions
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
  REAL :: ptsflx(nx,ny)        ! Surface heat flux (K*kg/(m**2*s))
  REAL :: qvsflx(nx,ny)        ! Surface moisture flux (kg/(m**2*s))

  REAL :: tem1  (nx,ny,nz)     ! Temporary work array
!
!-----------------------------------------------------------------------
!
!  Parameters describing this routine
!
!-----------------------------------------------------------------------
!
! 06/28/2002 Zuwen He
!
! fmtver??: to label each data a version.
! intver??: an integer to allow faster comparison than fmtver??,
!           which are strings.
!
! Verion 5.00: significant change in soil variables since version 4.10.
!
  CHARACTER (LEN=40) :: fmtver530
  PARAMETER (fmtver530='005.30 Binary Data')
!  PARAMETER (fmtver410='004.10 Binary Data')
!  PARAMETER (fmtver500='005.00 Binary Data')

  CHARACTER(LEN=10), PARAMETER :: tmunit = 'seconds   '
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k,l,is,nq
  INTEGER :: idummy
  REAL    :: rdummy

  INTEGER :: nxlg, nylg
  INTEGER :: n3d, istat

  REAL,    ALLOCATABLE :: out1d(:)
  REAL,    ALLOCATABLE :: out3d(:,:,:)
  INTEGER, ALLOCATABLE :: out3di(:,:,:)
  REAL,    ALLOCATABLE :: outtsoil(:,:,:,:)
  REAL,    ALLOCATABLE :: outqsoil(:,:,:,:)

!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  nxlg = nproc_x*(nx-3)+3
  nylg = nproc_y*(ny-3)+3

  n3d = MAX(nz, nzsoil, nstyps+1, 4)   ! 3rd dimenson for out3d

  IF (myproc == 0) THEN

     WRITE(6,'(/1x,a,a,a/)') 'Data dump format, fmtver = ',fmtver530,'. '

     WRITE(6,'(1x,a,f13.3/)') 'Writing history data at time = ', curtim
  !
  !-----------------------------------------------------------------------
  !
  !  Write header info
  !
  !-----------------------------------------------------------------------
  !
    WRITE(nchanl) fmtver530
    WRITE(nchanl) runname
    WRITE(nchanl) nocmnt
    IF( nocmnt > 0 ) THEN
      DO l=1,nocmnt
        WRITE(nchanl) cmnt(l)
      END DO
    END IF

    WRITE(nchanl) curtim,tmunit

    WRITE(nchanl) nxlg,nylg,nz,nzsoil
  !
  !-----------------------------------------------------------------------
  !
  !  Write the flags for different data groups.
  !
  !-----------------------------------------------------------------------
  !
    idummy = 0

    IF( grdbas == 1 ) THEN

      WRITE(nchanl)      1,      1,      0, mstout,      0,             &
                         0,      0,      0, landout,totout,             &
                         0, idummy, idummy, mapproj, month,             &
                       day,   year,   hour, minute, second

    ELSE

      WRITE(nchanl) grdout,  basout,varout, mstout, iceout,             &
                    trbout,  sfcout,rainout,landout,totout,             &
                    tkeout, nscalar,idummy, mapproj, month,             &
                       day,    year,  hour,  minute, second

    WRITE(nchanl)   p_qc,  p_qr,  p_qi,  p_qs,  p_qg, p_qh,             &
                    p_nc,  p_nr,  p_ni,  p_ns,  p_ng, p_nh,             &
                    p_zr,  p_zi,  p_zs,  p_zg,  p_zh,idummy,            &
                  idummy,idummy,idummy,idummy,idummy,idummy,            &
                  idummy,idummy,idummy,idummy,  p_cc,mphyopt


    END IF

    rdummy = 0.0
    WRITE(nchanl)   umove,   vmove, xgrdorg, ygrdorg, trulat1,          &
                  trulat2,  trulon,  sclfct,  ntcloud, n0rain,          &
                   n0snow, n0grpl, n0hail, rhoice, rhosnow,             &
                   rhogrpl, rhohail,  alpharain, alphaice, alphasnow,   &
                   alphagrpl, alphahail, rdummy,  rdummy, rdummy,       &
                    tstop, thisdmp, latitud, ctrlat,  ctrlon
  !
  !-----------------------------------------------------------------------
  !
  !  If totout=1, write additional parameters to history dump files.
  !  This is for ARPS version 4.1.2 or later.
  !
  !-----------------------------------------------------------------------
  !
    IF ( totout == 1 ) THEN
      WRITE(nchanl)  nstyp, prcout, radout, flxout,      0,  & ! 0 for snowcvr
                   snowout, idummy, idummy, idummy, idummy,             &
                    idummy, idummy, idummy, idummy, idummy,             &
                    idummy, idummy, idummy, idummy, idummy

      WRITE(nchanl) rdummy, rdummy, rdummy, rdummy, rdummy,             &
                    rdummy, rdummy, rdummy, rdummy, rdummy,             &
                    rdummy, rdummy, rdummy, rdummy, rdummy,             &
                    rdummy, rdummy, rdummy, rdummy, rdummy
    END IF

  END IF  ! myproc == 0

  ALLOCATE (out1d( MAX(nxlg,nylg) ),stat=istat)
  CALL check_alloc_status(istat, "binjoindump:out1d")

  ALLOCATE (out3d( nxlg,nylg, n3d ),stat=istat)
  CALL check_alloc_status(istat, "binjoindump:out3d")

  ALLOCATE (out3di( nxlg,nylg, nstyps ),stat=istat)
  CALL check_alloc_status(istat, "binjoindump:out3di")

  ALLOCATE (outtsoil( nxlg,nylg, nzsoil, 0:nstyps ),stat=istat)
  CALL check_alloc_status(istat, "binjoindump:outtsoil")

  ALLOCATE (outqsoil( nxlg,nylg, nzsoil, 0:nstyps ),stat=istat)
  CALL check_alloc_status(istat, "binjoindump:outqsoil")

!-----------------------------------------------------------------------
!
!  If grdout=1 or grdbas=1, write out grid variables
!
!-----------------------------------------------------------------------
!
  IF(grdout == 1 .OR. grdbas == 1 ) THEN

    CALL mpimerge1dx(x,nx,out1d)
    IF (myproc == 0) THEN
      WRITE(nchanl) 'x coord r1d1'
      WRITE(nchanl) out1d(1:nxlg)
    END IF

    CALL mpimerge1dy(y,ny,out1d)
    IF (myproc == 0) THEN
      WRITE(nchanl) 'y coord r1d2'
      WRITE(nchanl) out1d(1:nylg)
    END IF

    IF (myproc == 0) THEN
      WRITE(nchanl) 'z coord r1d3'
      WRITE(nchanl) z
    END IF

    CALL edgfill(zp,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz)
    CALL mpimerge3d(zp,nx,ny,nz,out3d)
    IF (myproc == 0) THEN
      WRITE(nchanl) 'zp coor r3d0'
      WRITE(nchanl) out3d(1:nxlg,1:nylg,1:nz)
    END IF

    CALL edgfill(zpsoil,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nzsoil,1,nzsoil)
    CALL mpimerge3d(zpsoil,nx,ny,nzsoil,out3d)
    IF (myproc == 0) THEN
      WRITE(nchanl) 'zpsoil rs3d0'
      WRITE(nchanl) out3d(1:nxlg,1:nylg,1:nzsoil)
    END IF

  END IF    ! grdout
!
!-----------------------------------------------------------------------
!
!  If basout=1, write out base state variables.
!
!-----------------------------------------------------------------------
!
  IF(basout == 1 .OR. grdbas == 1 ) THEN

    CALL edgfill(ubar,1,nx,1,nx, 1,ny,1,ny-1, 1,nz,1,nz-1)
    CALL mpimerge3d(ubar,nx,ny,nz,out3d)
    IF (myproc == 0) THEN
      WRITE(nchanl) 'ubar    r3d1'
      WRITE(nchanl) out3d(1:nxlg,1:nylg,1:nz)
    END IF

    CALL edgfill(vbar,1,nx,1,nx-1, 1,ny,1,ny, 1,nz,1,nz-1)
    CALL mpimerge3d(vbar,nx,ny,nz,out3d)
    IF (myproc == 0) THEN
      WRITE(nchanl) 'vbar    r3d2'
      WRITE(nchanl) out3d(1:nxlg,1:nylg,1:nz)
    END IF

    IF (myproc == 0) THEN
      out3d(:,:,:) = 0.0
      WRITE(nchanl) 'wbar    r3d3'
      WRITE(nchanl) out3d(1:nxlg,1:nylg,1:nz)
    END IF

    CALL edgfill(ptbar,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)
    CALL mpimerge3d(ptbar,nx,ny,nz,out3d)
    IF (myproc == 0) THEN
      WRITE(nchanl) 'ptbar   r3d0'
      WRITE(nchanl) out3d(1:nxlg,1:nylg,1:nz)
    END IF

    CALL edgfill(pbar,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)
    CALL mpimerge3d(pbar,nx,ny,nz,out3d)
    IF (myproc == 0) THEN
      WRITE(nchanl) 'pbar    r3d0'
      WRITE(nchanl) out3d(1:nxlg,1:nylg,1:nz)
    END IF

    IF(mstout == 1) THEN

      CALL edgfill(qvbar,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)
      CALL mpimerge3d(qvbar,nx,ny,nz,out3d)
      IF (myproc == 0) THEN
        WRITE(nchanl) 'qvbar   r3d0'
        WRITE(nchanl) out3d(1:nxlg,1:nylg,1:nz)
      END IF

    END IF

    IF(landout == 1) THEN

      IF( nstyp <= 1 ) THEN

        CALL iedgfill(soiltyp(1,1,1),1,nx,1,nx-1, 1,ny,1,ny-1,          &
                      1,1,1,1)

        CALL mpimerge2di(soiltyp(:,:,1),nx,ny,out3di)
        IF (myproc == 0) THEN
          WRITE(nchanl) 'soiltyp i2d0'
          WRITE(nchanl) out3di(1:nxlg,1:nylg,1)
        END IF

      ELSE
        DO is=1,nstyp

          CALL iedgfill(soiltyp(1,1,is),1,nx,1,nx-1, 1,ny,1,ny-1,       &
                        1,1,1,1)

          CALL edgfill(stypfrct(1,1,is),1,nx,1,nx-1, 1,ny,1,ny-1,       &
                       1,1,1,1)
        END DO

        CALL mpimerge3di(soiltyp,nx,ny,nstyp,out3di)
        CALL mpimerge3d (stypfrct,nx,ny,nstyp,out3d)

        IF (myproc == 0) THEN
          DO is=1,nstyp
            WRITE(nchanl) 'soiltyp i2d0'
            WRITE(nchanl) ((out3di(i,j,is),i=1,nxlg),j=1,nylg)

            WRITE(nchanl) 'stypfrc r2d0'
            WRITE(nchanl) ((out3d(i,j,is),i=1,nxlg),j=1,nylg)

          END DO
        END IF

      END IF

      CALL iedgfill(vegtyp ,1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)
      CALL mpimerge2di(vegtyp,nx,ny,out3di)
      IF (myproc == 0) THEN
        WRITE(nchanl) 'vegtyp  i2d0'
        WRITE(nchanl) out3di(1:nxlg,1:nylg,1)
      END IF

      CALL edgfill(lai    ,1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)
      CALL mpimerge2d(lai,nx,ny,out3d)
      IF (myproc == 0) THEN
        WRITE(nchanl) 'lai     r2d0'
        WRITE(nchanl) out3d(1:nxlg,1:nylg,1)
      END IF

      CALL edgfill(roufns ,1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)
      CALL mpimerge2d(roufns,nx,ny,out3d)
      IF (myproc == 0) THEN
        WRITE(nchanl) 'roufns  r2d0'
        WRITE(nchanl) out3d(1:nxlg,1:nylg,1)
      END IF

      CALL edgfill(veg    ,1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)
      CALL mpimerge2d(veg,nx,ny,out3d)
      IF (myproc == 0) THEN
        WRITE(nchanl) 'veg     r2d0'
        WRITE(nchanl) out3d(1:nxlg,1:nylg,1)
      END IF

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
!
!-----------------------------------------------------------------------
!
!  Write out u,v and w.
!
!-----------------------------------------------------------------------
!
  IF(varout == 1) THEN

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
            tem1(i,j,k) = u(i,j,k) - ubar(i,j,k)
          END DO
        END DO
      END DO
      CALL edgfill(tem1,1,nx,1,nx, 1,ny,1,ny-1, 1,nz,1,nz-1)
      CALL mpimerge3d(tem1,nx,ny,nz,out3d)
      IF (myproc == 0) THEN
        WRITE(nchanl) 'uprt    r3d1'
        WRITE(nchanl) out3d(1:nxlg,1:nylg,1:nz)
      END IF

      DO k=1,nz-1
        DO j=1,ny
          DO i=1,nx-1
            tem1(i,j,k) = v(i,j,k) - vbar(i,j,k)
          END DO
        END DO
      END DO
      CALL edgfill(tem1,1,nx,1,nx-1, 1,ny,1,ny, 1,nz,1,nz-1)
      CALL mpimerge3d(tem1,nx,ny,nz,out3d)
      IF (myproc == 0) THEN
        WRITE(nchanl) 'vprt    r3d2'
        WRITE(nchanl) out3d(1:nxlg,1:nylg,1:nz)
      END IF

      CALL edgfill(w,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz)
      CALL mpimerge3d(w,nx,ny,nz,out3d)
      IF (myproc == 0) THEN
        WRITE(nchanl) 'wprt    r3d3'
        WRITE(nchanl) out3d(1:nxlg,1:nylg,1:nz)
      END IF
!
!-----------------------------------------------------------------------
!
!  Write out scalars
!
!-----------------------------------------------------------------------
!
      CALL edgfill(ptprt,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)
      CALL mpimerge3d(ptprt,nx,ny,nz,out3d)
      IF (myproc == 0) THEN
        WRITE(nchanl) 'ptprt   r3d0'
        WRITE(nchanl) out3d(1:nxlg,1:nylg,1:nz)
      END IF

      CALL edgfill(pprt,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)
      CALL mpimerge3d(pprt,nx,ny,nz,out3d)
      IF (myproc == 0) THEN
        WRITE(nchanl) 'pprt    r3d0'
        WRITE(nchanl) out3d(1:nxlg,1:nylg,1:nz)
      END IF

    ELSE
!
!-----------------------------------------------------------------------
!
!  Write out total values to history dump
!
!-----------------------------------------------------------------------
!
      CALL edgfill(u,1,nx,1,nx, 1,ny,1,ny-1, 1,nz,1,nz-1)
      CALL mpimerge3d(u,nx,ny,nz,out3d)
      IF (myproc == 0) THEN
        WRITE(nchanl) 'u       r3d1'
        WRITE(nchanl) out3d(1:nxlg,1:nylg,1:nz)
      END IF

      CALL edgfill(v,1,nx,1,nx-1, 1,ny,1,ny, 1,nz,1,nz-1)
      CALL mpimerge3d(v,nx,ny,nz,out3d)
      IF (myproc == 0) THEN
        WRITE(nchanl) 'v       r3d2'
        WRITE(nchanl) out3d(1:nxlg,1:nylg,1:nz)
      END IF

      CALL edgfill(w,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz)
      CALL mpimerge3d(w,nx,ny,nz,out3d)
      IF (myproc == 0) THEN
        WRITE(nchanl) 'w       r3d3'
        WRITE(nchanl) out3d(1:nxlg,1:nylg,1:nz)
      END IF
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
      CALL mpimerge3d(tem1,nx,ny,nz,out3d)
      IF (myproc == 0) THEN
        WRITE(nchanl) 'pt      r3d0'
        WRITE(nchanl) out3d(1:nxlg,1:nylg,1:nz)
      END IF

      DO k=1,nz-1
        DO j=1,ny-1
          DO i=1,nx-1
            tem1(i,j,k) = pbar(i,j,k) + pprt(i,j,k)
          END DO
        END DO
      END DO
      CALL edgfill(tem1,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)
      CALL mpimerge3d(tem1,nx,ny,nz,out3d)
      IF (myproc == 0) THEN
        WRITE(nchanl) 'p       r3d0'
        WRITE(nchanl) out3d(1:nxlg,1:nylg,1:nz)
      END IF

    END IF

  END IF     ! varout
!
!-----------------------------------------------------------------------
!
!  If mstout = 1, write out moisture scalars.
!
!-----------------------------------------------------------------------
!
  IF(mstout == 1) THEN

    IF( totout == 0 ) THEN
!
!-----------------------------------------------------------------------
!
!  Write out perturbation to history dump
!
!-----------------------------------------------------------------------
!
      DO k=1,nz-1
        DO j=1,ny-1
          DO i=1,nx-1
            tem1(i,j,k)=qv(i,j,k)-qvbar(i,j,k)
          END DO
        END DO
      END DO

      CALL edgfill(tem1,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)
      CALL mpimerge3d(tem1,nx,ny,nz,out3d)
      IF (myproc == 0) THEN
        WRITE(nchanl) 'qvprt   r3d0'
        WRITE(nchanl) out3d(1:nxlg,1:nylg,1:nz)
      END IF

    ELSE
!
!-----------------------------------------------------------------------
!
!  Write out total values to history dump
!
!-----------------------------------------------------------------------
!
      CALL edgfill(qv,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)
      CALL mpimerge3d(qv,nx,ny,nz,out3d)
      IF (myproc == 0) THEN
        WRITE(nchanl) 'qv      r3d0'
        WRITE(nchanl) out3d(1:nxlg,1:nylg,1:nz)
      END IF

    END IF

    DO nq = 1,nscalar
      CALL edgfill(qscalar(:,:,:,nq),1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)
      CALL mpimerge3d(qscalar(:,:,:,nq),nx,ny,nz,out3d)
      IF (myproc == 0) THEN
        WRITE(nchanl) qnames(nq)//'    r3d0'
        WRITE(nchanl) out3d(1:nxlg,1:nylg,1:nz)
      END IF
    END DO

    IF(rainout == 1) THEN

      CALL edgfill(raing,   1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)
      CALL mpimerge2d(raing,nx,ny,out3d)
      IF (myproc == 0) THEN
        WRITE(nchanl) 'raing   r2d0'
        WRITE(nchanl) out3d(1:nxlg,1:nylg,1)
      END IF

      CALL edgfill(rainc,1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)
      CALL mpimerge2d(rainc,nx,ny,out3d)
      IF (myproc == 0) THEN
        WRITE(nchanl) 'rainc   r2d0'
        WRITE(nchanl) out3d(1:nxlg,1:nylg,1)
      END IF

    END IF   !rainout

    IF ( prcout == 1 ) THEN
      CALL edgfill(prcrate,1,nx,1,nx-1, 1,ny,1,ny-1, 1,4,1,4)
      CALL mpimerge3d(prcrate,nx,ny,4,out3d)
      IF (myproc == 0) THEN
        WRITE(nchanl) 'prcrat1 r2d0'
        WRITE(nchanl) out3d(1:nxlg,1:nylg,1)
        WRITE(nchanl) 'prcrat2 r2d0'
        WRITE(nchanl) out3d(1:nxlg,1:nylg,2)
        WRITE(nchanl) 'prcrat3 r2d0'
        WRITE(nchanl) out3d(1:nxlg,1:nylg,3)
        WRITE(nchanl) 'prcrat4 r2d0'
        WRITE(nchanl) out3d(1:nxlg,1:nylg,4)
      END IF
    END IF

  END IF   !mstout
!
!-----------------------------------------------------------------------
!
!  If tkeout = 1, write out tke.
!
!-----------------------------------------------------------------------
!
  IF( tkeout == 1 ) THEN

    CALL edgfill(tke,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)
    CALL mpimerge3d(tke,nx,ny,nz,out3d)
    IF (myproc == 0) THEN
      WRITE(nchanl) 'tke     r3d0'
      WRITE(nchanl) out3d(1:nxlg,1:nylg,1:nz)
    END IF

  END IF
!
!-----------------------------------------------------------------------
!
!  If trbout = 1, write out the turbulence parameter, km.
!
!-----------------------------------------------------------------------
!
  IF( trbout == 1 ) THEN

    CALL edgfill(kmh,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)
    CALL mpimerge3d(kmh,nx,ny,nz,out3d)
    IF (myproc == 0) THEN
      WRITE(nchanl) 'kmh     r3d0'
      WRITE(nchanl) out3d(1:nxlg,1:nylg,1:nz)
    END IF

    CALL edgfill(kmv,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)
    CALL mpimerge3d(kmv,nx,ny,nz,out3d)
    IF (myproc == 0) THEN
      WRITE(nchanl) 'kmv     r3d0'
      WRITE(nchanl) out3d(1:nxlg,1:nylg,1:nz)
    END IF

  END IF   ! trbout

!
!-----------------------------------------------------------------------
!
!  If sfcout = 1, write out the surface variables,
!  tsoil, qsoil, and wetcanp.
!
!-----------------------------------------------------------------------
!
  IF( sfcout == 1) THEN

    IF( nstyp <= 1 ) THEN

      CALL edgfill(tsoil(1,1,1,0),  1,nx,1,nx-1, 1,ny,1,ny-1,             &
                   1,nzsoil,1,nzsoil)
      CALL edgfill(qsoil(1,1,1,0), 1,nx,1,nx-1, 1,ny,1,ny-1,              &
                   1,nzsoil,1,nzsoil)

      CALL mpimerge3d(tsoil(:,:,:,0),nx,ny,nzsoil,out3d)
      IF (myproc == 0) THEN
        WRITE(nchanl) 'tsoil  rs3d0'
        WRITE(nchanl) out3d(1:nxlg,1:nylg,1:nzsoil)
      END IF

      CALL mpimerge3d(qsoil(:,:,:,0),nx,ny,nzsoil,out3d)
      IF (myproc == 0) THEN
        WRITE(nchanl) 'qsoil  rs3d0'
        WRITE(nchanl) out3d(1:nxlg,1:nylg,1:nzsoil)
      END IF


      CALL edgfill(wetcanp(1,1,0),1,nx,1,nx-1, 1,ny,1,ny-1,             &
                   1,1,1,1)
      CALL mpimerge2d(wetcanp(:,:,0),nx,ny,out3d)
      IF (myproc == 0) THEN
        WRITE(nchanl) 'wetcanp r2d0'
        WRITE(nchanl) out3d(1:nxlg,1:nylg,1)
      END IF

    ELSE
      DO is=0,nstyp

         CALL edgfill(tsoil(1,1,1,is),  1,nx,1,nx-1, 1,ny,1,ny-1,        &
                      1,nzsoil,1,nzsoil)
         CALL edgfill(qsoil(1,1,1,is),  1,nx,1,nx-1, 1,ny,1,ny-1,        &
                      1,nzsoil,1,nzsoil)
         CALL edgfill(wetcanp(1,1,is),  1,nx,1,nx-1, 1,ny,1,ny-1,        &
                      1,1,1,1)

      END DO

      CALL mpimerge4d(tsoil,nx,ny,nzsoil,nstyp+1,outtsoil)
      CALL mpimerge4d(qsoil,nx,ny,nzsoil,nstyp+1,outqsoil)
      CALL mpimerge3d(wetcanp,nx,ny,nstyp+1,out3d)

      IF (myproc == 0) THEN

        DO is=0,nstyp

          WRITE(nchanl) 'tsoil  rs3d0'
          WRITE(nchanl) outtsoil(1:nxlg,1:nylg,1:nzsoil,is)
          WRITE(nchanl) 'qsoil  rs3d0'
          WRITE(nchanl) outqsoil(1:nxlg,1:nylg,1:nzsoil,is)

          WRITE(nchanl) 'wetcanp r2d0'
          WRITE(nchanl) out3d(1:nxlg,1:nylg,is+1)
        END DO

      END IF

    END IF

    IF (snowout == 1) THEN

      CALL edgfill(snowdpth,1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)
      CALL mpimerge2d(snowdpth,nx,ny,out3d)
      IF (myproc == 0) THEN
        WRITE(nchanl) 'snowdpthr2d0'
        WRITE(nchanl) out3d(1:nxlg,1:nylg,1)
      END IF

    END IF

  END IF
!
!-----------------------------------------------------------------------
!
!  If radout = 1, write out the radiation arrays
!
!-----------------------------------------------------------------------
!
  IF( radout == 1 ) THEN

    CALL edgfill(radfrc,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)
    CALL mpimerge3d(radfrc,nx,ny,nz,out3d)
    IF (myproc == 0) THEN
      WRITE(nchanl) 'radfrc  r3d0'
      WRITE(nchanl) out3d(1:nxlg,1:nylg,1:nz)
    END IF

    CALL edgfill(radsw,1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)
    CALL mpimerge2d(radsw,nx,ny,out3d)
    IF (myproc == 0) THEN
      WRITE(nchanl) 'radsw   r2d0'
      WRITE(nchanl) out3d(1:nxlg,1:nylg,1)
    END IF

    CALL edgfill(rnflx,1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)
    CALL mpimerge2d(rnflx,nx,ny,out3d)
    IF (myproc == 0) THEN
      WRITE(nchanl) 'rnflx   r2d0'
      WRITE(nchanl) out3d(1:nxlg,1:nylg,1)
    END IF

    CALL edgfill(radswnet,1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)
    CALL mpimerge2d(radswnet,nx,ny,out3d)
    IF (myproc == 0) THEN
      WRITE(nchanl) 'radswnetr2d0'
      WRITE(nchanl) out3d(1:nxlg,1:nylg,1)
    END IF

    CALL edgfill(radlwin,1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)
    CALL mpimerge2d(radlwin,nx,ny,out3d)
    IF (myproc == 0) THEN
      WRITE(nchanl) 'radlwin r2d0'
      WRITE(nchanl) out3d(1:nxlg,1:nylg,1)
    END IF

  END IF   ! radout
!
!-----------------------------------------------------------------------
!
!  If flxout = 1, write out the surface fluxes
!
!-----------------------------------------------------------------------
!
  IF( flxout == 1 ) THEN

    CALL edgfill(usflx,1,nx,1,nx, 1,ny,1,ny-1, 1,1,1,1)
    CALL mpimerge2d(usflx,nx,ny,out3d)
    IF (myproc == 0) THEN
      WRITE(nchanl) 'usflx   r2d0'
      WRITE(nchanl) out3d(1:nxlg,1:nylg,1)
    END IF

    CALL edgfill(vsflx,1,nx,1,nx-1, 1,ny,1,ny, 1,1,1,1)
    CALL mpimerge2d(vsflx,nx,ny,out3d)
    IF (myproc == 0) THEN
      WRITE(nchanl) 'vsflx   r2d0'
      WRITE(nchanl) out3d(1:nxlg,1:nylg,1)
    END IF

    CALL edgfill(ptsflx,1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)
    CALL mpimerge2d(ptsflx,nx,ny,out3d)
    IF (myproc == 0) THEN
      WRITE(nchanl) 'ptsflx  r2d0'
      WRITE(nchanl) out3d(1:nxlg,1:nylg,1)
    END IF

    CALL edgfill(qvsflx,1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)
    CALL mpimerge2d(qvsflx,nx,ny,out3d)
    IF (myproc == 0) THEN
      WRITE(nchanl) 'qvsflx  r2d0'
      WRITE(nchanl) out3d(1:nxlg,1:nylg,1)
    END IF

  END IF   ! flxout

  DEALLOCATE(out1d)
  DEALLOCATE(out3d, out3di)
  DEALLOCATE(outtsoil, outqsoil)

  RETURN
END SUBROUTINE binjoindump
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE BN2DUMP                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE bn2dump(nx,ny,nz,nzsoil,nstyps, nchanl, grdbas,              &
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
!  Write history data into channel nchanl as binary data.
!
!  This routine can dump the data arrays in a model subdomain and
!  at selected data points.
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
!  3/10/92.
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
!  02/06/95 (Y. Liu)
!  Added map projection parameters into the second binary dumping
!
!  12/09/1998 (Donghai Wang)
!  Added the snow cover.
!
!  05/15/2002 (J. Brotzge)
!  Added to allow for multiple soil schemes.
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
!    zpsoil   Vertical coordinate of grid points in soil (m)
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
!
!-----------------------------------------------------------------------
!
!  The following parameters are passed into this subroutine through
!  a common block in globcst.inc.  These parameters determine which
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
  REAL :: zpsoil(nx,ny,nzsoil) ! The physical height coordinate defined at
                               ! w-point of the soil

  INTEGER :: nstyps                  ! Number of soil types
  INTEGER :: soiltyp(nx,ny,nstyps)   ! Soil type
  REAL :: stypfrct(nx,ny,nstyps)     ! Soil type
  INTEGER :: vegtyp(nx,ny)           ! Vegetation type
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
  REAL :: ptsflx(nx,ny)        ! Surface heat flux (K*kg/(m**2*s))
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
  PARAMETER (fmtver='004.10 2nd Binary Data')
  CHARACTER (LEN=10) :: tmunit
  PARAMETER (tmunit='seconds   ')
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k,l,is, idummy
  REAL :: rdummy

  INTEGER :: nxout,nyout,nzout ! The size of array to be written out.
  INTEGER :: nzsoilout         ! The size of array to be written out.


  INTEGER :: ist ,ind ,isk ,jst ,jnd ,jsk ,kst ,knd ,ksk
  INTEGER :: ist1,ind1,isk1,jst1,jnd1,jsk1,kst1,knd1,ksk1

  INTEGER :: lst, lnd, lsk

  INTEGER :: setdomn,setskip
  SAVE setdomn, setskip
  SAVE ist,ind,isk,jst,jnd,jsk,kst,knd,ksk,lst,lnd,lsk

  DATA setdomn, setskip /0,0/
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

  IF( setdomn == 0) THEN  ! If these parameters are nevers set ...

    ist = 1
    jst = 1
    kst = 1
    lst = 1
    ind = nx
    jnd = ny
    knd = nz
    lnd = nzsoil

  END IF

  IF( setskip == 0) THEN  ! If these parameters are nevers set ...

    isk = 1
    jsk = 1
    ksk = 1
    lsk = 1

  END IF

  WRITE(6,'(1x,a,f13.3/)') 'Writing history data at time=', curtim
!
!-----------------------------------------------------------------------
!
!  Write header info
!
!-----------------------------------------------------------------------
!
  WRITE(nchanl) fmtver
  WRITE(nchanl) runname
  WRITE(nchanl) nocmnt
  IF( nocmnt > 0 ) THEN
    DO l=1,nocmnt
      WRITE(nchanl) cmnt(l)
    END DO
  END IF
  WRITE(nchanl) curtim,tmunit

  nxout = ist+(ind-ist)/isk
  nyout = jst+(jnd-jst)/jsk
  nzout = kst+(knd-kst)/ksk
  nzsoilout = lst+(lnd-lst)/lsk

  WRITE(nchanl) nxout, nyout, nzout, nzsoilout
  PRINT*,'nxout= ',nxout,'nyout= ',nyout,'nzout= ',nzout,            &
         'nzsoilout= ',nzsoilout
!
!-----------------------------------------------------------------------
!
!  Write the flags for different data groups.
!
!-----------------------------------------------------------------------
!
  idummy = 0

  IF( grdbas == 1 ) THEN

    WRITE(nchanl)      1,      1,      0, mstout,      0,               &
                       0, idummy, idummy, landout, totout,              &
                  idummy, idummy, idummy, mapproj, month,               &
                     day,   year,   hour, minute, second

  ELSE

    WRITE(nchanl) grdout, basout, varout, mstout, iceout,               &
                  trbout, sfcout, rainout,landout,totout,               &
                  tkeout, idummy, idummy, mapproj, month,               &
                     day,   year,   hour, minute, second

  END IF

  rdummy = 0.0
  WRITE(nchanl)   umove,   vmove, xgrdorg, ygrdorg, trulat1,            &
                trulat2,  trulon,  sclfct,  rdummy,  rdummy,            &
                 rdummy,  rdummy,  rdummy,  rdummy,  rdummy,            &
                  tstop, thisdmp,  latitud, ctrlat,  ctrlon

  IF ( totout /= 0 ) THEN
!
!-----------------------------------------------------------------------
!
!  Add new parameters for new version of history data dump
!
!-----------------------------------------------------------------------
!
    WRITE(nchanl) nstyp,  prcout, radout, flxout,      0,  & ! 0 for snowcvr
              snowout,idummy, idummy, idummy, idummy,                   &
                  idummy, idummy, idummy, idummy, idummy,               &
                  idummy, idummy, idummy, idummy, idummy

    WRITE(nchanl) rdummy, rdummy, rdummy, rdummy, rdummy,               &
                  rdummy, rdummy, rdummy, rdummy, rdummy,               &
                  rdummy, rdummy, rdummy, rdummy, rdummy,               &
                  rdummy, rdummy, rdummy, rdummy, rdummy

  END IF
!
!-----------------------------------------------------------------------
!
!  If grdout=1 or grdbas=1, write out grid variables
!
!-----------------------------------------------------------------------
!
  IF(grdout == 1 .OR. grdbas == 1 ) THEN

    WRITE(nchanl) 'x coordinate'
    WRITE(nchanl) ( x(i), i=ist,ind,isk)

    WRITE(nchanl) 'y coordinate'
    WRITE(nchanl) ( y(j), j=jst,jnd,jsk)

    WRITE(nchanl) 'z coordinate'
    WRITE(nchanl) ( z(k), k=kst,knd,ksk)

    WRITE(nchanl) 'zp coord    '
    WRITE(nchanl) ((( zp(i,j,k),                                        &
                  i=ist,ind,isk),j=jst,jnd,jsk),k=kst,knd,ksk)

    WRITE(nchanl) 'zpsoil coord    '
    WRITE(nchanl) (((zpsoil(i,j,l),i=ist,ind,isk),j=jst,jnd,jsk),       &
                  l=lst,lnd,lsk)

  END IF    ! grdout
!
!-----------------------------------------------------------------------
!
!  If basout=1, write out base state variables.
!
!-----------------------------------------------------------------------
!
  IF(basout == 1 .OR. grdbas == 1 ) THEN

    WRITE(nchanl) 'ubar        '
    WRITE(nchanl) ((( ubar(i,j,k),                                      &
                  i=ist,ind,isk),j=jst,jnd,jsk),k=kst,knd,ksk)

    WRITE(nchanl) 'vbar        '
    WRITE(nchanl) ((( vbar(i,j,k),                                      &
                  i=ist,ind,isk),j=jst,jnd,jsk),k=kst,knd,ksk)

    DO k=kst,knd,ksk
      DO j=jst,jnd,jsk
        DO i=ist,ind,isk
          tem1(i,j,k) = 0.0
        END DO
      END DO
    END DO
    WRITE(nchanl) 'wbar        '
    WRITE(nchanl) (((tem1(i,j,k),                                       &
                  i=ist,ind,isk),j=jst,jnd,jsk),k=kst,knd,ksk)

    WRITE(nchanl) 'ptbar       '
    WRITE(nchanl) (((ptbar(i,j,k),                                      &
                  i=ist,ind,isk),j=jst,jnd,jsk),k=kst,knd,ksk)

    WRITE(nchanl) 'pbar        '
    WRITE(nchanl) (((pbar(i,j,k),                                       &
                  i=ist,ind,isk),j=jst,jnd,jsk),k=kst,knd,ksk)

    IF(mstout == 1) THEN

      WRITE(nchanl) 'qvbar       '
      WRITE(nchanl) (((qvbar(i,j,k),                                    &
                  i=ist,ind,isk),j=jst,jnd,jsk),k=kst,knd,ksk)

    END IF

    IF(landout == 1) THEN

      IF( nstyp <= 1 ) THEN
        WRITE(nchanl) 'soiltyp     '
        WRITE(nchanl) ((soiltyp(i,j,1),i=ist,ind,isk),                  &
                        j=jst,jnd,jsk)
      ELSE
        DO is=1,nstyp
          WRITE(nchanl) 'soiltyp     '
          WRITE(nchanl) ((soiltyp(i,j,is),i=ist,ind,isk),               &
                        j=jst,jnd,jsk)

          WRITE(nchanl) 'stypfrct    '
          WRITE(nchanl) ((stypfrct(i,j,is),i=ist,ind,isk),              &
                        j=jst,jnd,jsk)
        END DO
      END IF

      WRITE(nchanl) 'vegtyp      '
      WRITE(nchanl) ((vegtyp (i,j),i=ist,ind,isk),j=jst,jnd,jsk)

      WRITE(nchanl) 'lai         '
      WRITE(nchanl) ((lai    (i,j),i=ist,ind,isk),j=jst,jnd,jsk)

      WRITE(nchanl) 'roufns      '
      WRITE(nchanl) ((roufns (i,j),i=ist,ind,isk),j=jst,jnd,jsk)

      WRITE(nchanl) 'veg         '
      WRITE(nchanl) ((veg    (i,j),i=ist,ind,isk),j=jst,jnd,jsk)

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
!
!-----------------------------------------------------------------------
!
!  Write out u, v and w
!
!-----------------------------------------------------------------------
!
  IF(varout == 1) THEN

    IF ( totout == 0 ) THEN
!
!-----------------------------------------------------------------------
!
!    Write out perturbatios to history dump
!
!-----------------------------------------------------------------------
!
      WRITE(nchanl) 'uprt        '
      WRITE(nchanl) ((( u(i,j,k)-ubar(i,j,k),                           &
                    i=ist,ind,isk),j=jst,jnd,jsk),k=kst,knd,ksk)

      WRITE(nchanl) 'vprt        '
      WRITE(nchanl) ((( v(i,j,k)-vbar(i,j,k),                           &
                    i=ist,ind,isk),j=jst,jnd,jsk),k=kst,knd,ksk)

      WRITE(nchanl) 'wprt        '
      WRITE(nchanl) ((( w(i,j,k),                                       &
                    i=ist,ind,isk),j=jst,jnd,jsk),k=kst,knd,ksk)
!
!-----------------------------------------------------------------------
!
!  Write out scalars
!
!-----------------------------------------------------------------------
!
      WRITE(nchanl) 'ptprt       '
      WRITE(nchanl) ((( ptprt(i,j,k),                                   &
                    i=ist,ind,isk),j=jst,jnd,jsk),k=kst,knd,ksk)

      WRITE(nchanl) 'pprt        '
      WRITE(nchanl) ((( pprt(i,j,k),                                    &
                    i=ist,ind,isk),j=jst,jnd,jsk),k=kst,knd,ksk)

    ELSE
!
!-----------------------------------------------------------------------
!
!    Write out total values to history dump
!
!-----------------------------------------------------------------------
!
      WRITE(nchanl) 'u           '
      WRITE(nchanl) ((( u(i,j,k),                                       &
                    i=ist,ind,isk),j=jst,jnd,jsk),k=kst,knd,ksk)

      WRITE(nchanl) 'v           '
      WRITE(nchanl) ((( v(i,j,k),                                       &
                    i=ist,ind,isk),j=jst,jnd,jsk),k=kst,knd,ksk)

      WRITE(nchanl) 'w           '
      WRITE(nchanl) ((( w(i,j,k),                                       &
                    i=ist,ind,isk),j=jst,jnd,jsk),k=kst,knd,ksk)
!
!-----------------------------------------------------------------------
!
!  Write out scalars
!
!-----------------------------------------------------------------------
!
      WRITE(nchanl) 'pt          '
      WRITE(nchanl) ((( ptprt(i,j,k)+ptbar(i,j,k),                      &
                    i=ist,ind,isk),j=jst,jnd,jsk),k=kst,knd,ksk)

      WRITE(nchanl) 'p           '
      WRITE(nchanl) ((( pprt(i,j,k)+pbar(i,j,k),                        &
                    i=ist,ind,isk),j=jst,jnd,jsk),k=kst,knd,ksk)

    END IF

  END IF     ! varout
!
!-----------------------------------------------------------------------
!
!  If mstout = 1, write out moisture scalars.
!
!-----------------------------------------------------------------------
!
  IF(mstout == 1) THEN

    IF ( totout == 0 ) THEN
!
!-----------------------------------------------------------------------
!
!    Write out perturbations to history dump
!
!-----------------------------------------------------------------------
!
      WRITE(nchanl) 'qvprt       '
      WRITE(nchanl) ((( qv(i,j,k)-qvbar(i,j,k),                         &
                    i=ist,ind,isk),j=jst,jnd,jsk),k=kst,knd,ksk)

    ELSE
!
!-----------------------------------------------------------------------
!
!    Write out total values to history dump
!
!-----------------------------------------------------------------------
!
      WRITE(nchanl) 'qv          '
      WRITE(nchanl) ((( qv(i,j,k),                                      &
                    i=ist,ind,isk),j=jst,jnd,jsk),k=kst,knd,ksk)

    END IF

    WRITE(nchanl) 'qc          '
    WRITE(nchanl) ((( qc(i,j,k),                                        &
                  i=ist,ind,isk),j=jst,jnd,jsk),k=kst,knd,ksk)

    WRITE(nchanl) 'qr          '
    WRITE(nchanl) ((( qr(i,j,k),                                        &
                  i=ist,ind,isk),j=jst,jnd,jsk),k=kst,knd,ksk)

    IF(rainout == 1) THEN

      WRITE(nchanl) 'raing       '
      WRITE(nchanl) ((raing(i,j),i=ist,ind,isk),j=jst,jnd,jsk)

      WRITE(nchanl) 'rainc       '
      WRITE(nchanl) ((rainc(i,j),i=ist,ind,isk),j=jst,jnd,jsk)

    END IF   !rainout

    IF ( prcout == 1 ) THEN

      WRITE(nchanl) 'prcrate1    '
      WRITE(nchanl) ((prcrate(i,j,1),i=ist,ind,isk),j=jst,jnd,jsk)
      WRITE(nchanl) 'prcrate2    '
      WRITE(nchanl) ((prcrate(i,j,2),i=ist,ind,isk),j=jst,jnd,jsk)
      WRITE(nchanl) 'prcrate3    '
      WRITE(nchanl) ((prcrate(i,j,3),i=ist,ind,isk),j=jst,jnd,jsk)
      WRITE(nchanl) 'prcrate4    '
      WRITE(nchanl) ((prcrate(i,j,4),i=ist,ind,isk),j=jst,jnd,jsk)

    END IF   ! prcout

    IF(iceout == 1) THEN

      WRITE(nchanl) 'qi          '
      WRITE(nchanl) ((( qi(i,j,k),                                      &
                  i=ist,ind,isk),j=jst,jnd,jsk),k=kst,knd,ksk)

      WRITE(nchanl) 'qs          '
      WRITE(nchanl) ((( qs(i,j,k),                                      &
                  i=ist,ind,isk),j=jst,jnd,jsk),k=kst,knd,ksk)

      WRITE(nchanl) 'qh          '
      WRITE(nchanl) ((( qh(i,j,k),                                      &
                  i=ist,ind,isk),j=jst,jnd,jsk),k=kst,knd,ksk)

    END IF   !iceout

  END IF   !mstout
!
!-----------------------------------------------------------------------
!
!  If tkeout = 1, write out turbulence parameter, km.
!
!-----------------------------------------------------------------------
!
  IF( tkeout == 1 ) THEN

    WRITE(nchanl) 'tke          '
    WRITE(nchanl) ((( tke(i,j,k),                                       &
                  i=ist,ind,isk),j=jst,jnd,jsk),k=kst,knd,ksk)

  END IF   ! tkeout

!
!-----------------------------------------------------------------------
!
!  If trbout = 1, write out turbulence parameter, km.
!
!-----------------------------------------------------------------------
!
  IF( trbout == 1 ) THEN

    WRITE(nchanl) 'kmh          '
    WRITE(nchanl) ((( kmh(i,j,k),                                       &
                  i=ist,ind,isk),j=jst,jnd,jsk),k=kst,knd,ksk)

    WRITE(nchanl) 'kmv          '
    WRITE(nchanl) ((( kmv(i,j,k),                                       &
                  i=ist,ind,isk),j=jst,jnd,jsk),k=kst,knd,ksk)

  END IF   ! trbout
!
!-----------------------------------------------------------------------
!
!  If sfcout = 1, write out the surface variables, tsoil,
!  qsoil, and wetcanp.
!
!-----------------------------------------------------------------------
!
  IF( sfcout == 1) THEN

    IF( nstyp <= 1 ) THEN

      WRITE(nchanl) 'tsoil         '
      WRITE(nchanl) (((tsoil(i,j,l,0),i=ist,ind,isk),j=jst,jnd,jsk),  &
                     l=lst,lnd,lsk)

      WRITE(nchanl) 'qsoil       '
      WRITE(nchanl) (((qsoil(i,j,l,0),i=ist,ind,isk),j=jst,jnd,jsk),  &
                     l=lst,lnd,lsk)

      WRITE(nchanl) 'wetcanp      '
      WRITE(nchanl) ((wetcanp(i,j,0),i=ist,ind,isk),j=jst,jnd,jsk)

    ELSE

      DO is=0,nstyp
        WRITE(nchanl) 'tsoil         '
        WRITE(nchanl) (((tsoil(i,j,l,is),i=ist,ind,isk),                   &
                      j=jst,jnd,jsk),l=lst,lnd,lsk)

        WRITE(nchanl) 'qsoil       '
        WRITE(nchanl) (((qsoil(i,j,l,is),i=ist,ind,isk),                  &
                      j=jst,jnd,jsk),l=lst,lnd,lsk)

        WRITE(nchanl) 'wetcanp      '
        WRITE(nchanl) ((wetcanp(i,j,is),i=ist,ind,isk),                 &
                      j=jst,jnd,jsk)
      END DO

    END IF

    IF(snowout == 1) THEN

      WRITE(nchanl) 'snowdpth     '
      WRITE(nchanl) ((snowdpth(i,j),i=ist,ind,isk),                     &
                    j=jst,jnd,jsk)
    END IF

  END IF   ! sfcout done
!
!-----------------------------------------------------------------------
!
!  If radout = 1, write out radiation arrays, radfrc, radsw and
!  rnflx.
!
!-----------------------------------------------------------------------
!
  IF( radout == 1 ) THEN

    WRITE(nchanl) 'radfrc       '
    WRITE(nchanl) ((( radfrc(i,j,k),                                    &
                  i=ist,ind,isk),j=jst,jnd,jsk),k=kst,knd,ksk)

    WRITE(nchanl) 'radsw        '
    WRITE(nchanl) (( radsw(i,j),                                        &
                  i=ist,ind,isk),j=jst,jnd,jsk)

    WRITE(nchanl) 'rnflx'
    WRITE(nchanl) (( rnflx(i,j),                                        &
                  i=ist,ind,isk),j=jst,jnd,jsk)

    WRITE(nchanl) 'radswnet        '
    WRITE(nchanl) (( radswnet(i,j),                                     &
                  i=ist,ind,isk),j=jst,jnd,jsk)

    WRITE(nchanl) 'radlwin        '
    WRITE(nchanl) (( radlwin(i,j),                                      &
                  i=ist,ind,isk),j=jst,jnd,jsk)


  END IF   ! radout
!
!-----------------------------------------------------------------------
!
!  If flxout = 1, write out surface fluxes
!
!-----------------------------------------------------------------------
!
  IF( flxout == 1 ) THEN

    WRITE(nchanl) 'usflx        '
    WRITE(nchanl) (( usflx(i,j),i=ist,ind,isk),j=jst,jnd,jsk)

    WRITE(nchanl) 'vsflx        '
    WRITE(nchanl) (( vsflx(i,j),i=ist,ind,isk),j=jst,jnd,jsk)

    WRITE(nchanl) 'ptsflx       '
    WRITE(nchanl) (( ptsflx(i,j),i=ist,ind,isk),j=jst,jnd,jsk)

    WRITE(nchanl) 'qvsflx       '
    WRITE(nchanl) (( qvsflx(i,j),i=ist,ind,isk),j=jst,jnd,jsk)

  END IF   ! flxout

  RETURN

  ENTRY bdmpdomn(ist1,ind1,jst1,jnd1,kst1,knd1)
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

  RETURN

  ENTRY bdmpskip(isk1, jsk1, ksk1)
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
END SUBROUTINE bn2dump
