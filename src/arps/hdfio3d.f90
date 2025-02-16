!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE HDFREAD                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################

SUBROUTINE hdfread(nx,ny,nz,nzsoil,nstyps, grdbas, filename, time,      &
           x,y,z,zp,zpsoil,                                             &
           uprt, vprt, wprt, ptprt, pprt,                               &
           qvprt, qscalar, tke, kmh,kmv,                                &
           ubar, vbar, wbar, ptbar, pbar, rhobar, qvbar,                &
           soiltyp,stypfrct,vegtyp,lai,roufns,veg,                      &
           tsoil,qsoil,wetcanp,snowdpth,                                &
           raing,rainc,prcrate,                                         &
           radfrc,radsw,rnflx,radswnet,radlwin,                         &
           usflx,vsflx,ptsflx,qvsflx,                                   &
           ireturn, tem1)

!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Read in history data in the NCSA HDF4 format.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Gene Bassett
!  2000/04/15
!
!  MODIFICATION HISTORY:
!
!  05/15/2002  (J. Brotzge)
!  Added to allow for multiple soil schemes
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
!    filename  Character variable nhming the input HDF file
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
!    wprt     Vertical component of perturbation velocity in Cartesian
!             coordinates (m/s).
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
!
!    WORK ARRAY
!
!    tem1
!    tem2     work array.
!
!-----------------------------------------------------------------------

  USE arps_precision

  IMPLICIT NONE

!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------

  INCLUDE 'globcst.inc'
  INCLUDE 'grid.inc'        ! Grid parameters
  INCLUDE 'indtflg.inc'
  INCLUDE 'alloc.inc'       ! allocation parameters & declarations
  INCLUDE 'phycst.inc'
  INCLUDE 'mp.inc'

!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------

  INTEGER :: nx,ny,nz
  INTEGER :: nzsoil

  INTEGER :: grdbas
  CHARACTER (LEN=*) :: filename
  REAL :: time

  REAL :: x     (nx)           ! x coord.
  REAL :: y     (ny)           ! y coord.
  REAL :: z     (nz)           ! z coord.
  REAL :: zp    (nx,ny,nz)     ! physical x coord.
  REAL :: zpsoil(nx,ny,nzsoil) ! physical x coord. for soil (m)

  REAL :: uprt  (nx,ny,nz)     ! Perturbation u-velocity (m/s)
  REAL :: vprt  (nx,ny,nz)     ! Perturbation v-velocity (m/s)
  REAL :: wprt  (nx,ny,nz)     ! Perturbation w-velocity (m/s)
  REAL :: ptprt (nx,ny,nz)     ! Perturbation potential temperature (K)
  REAL :: pprt  (nx,ny,nz)     ! Perturbation pressure (Pascal)

  REAL :: qvprt (nx,ny,nz)     ! Perturbation water vapor mixing ratio (kg/kg)

  REAL :: qscalar    (nx,ny,nz,nscalar)

  REAL :: tke   (nx,ny,nz)     ! Turbulent Kinetic Energy ((m/s)**2)
  REAL :: kmh   (nx,ny,nz)     ! Horizontal turb. mixing coef. for
                               ! momentum. ( m**2/s )
  REAL :: kmv   (nx,ny,nz)     ! Vertical turb. mixing coef. for
                               ! momentum. ( m**2/s )

  REAL :: ubar  (nx,ny,nz)     ! Base state u-velocity (m/s)
  REAL :: vbar  (nx,ny,nz)     ! Base state v-velocity (m/s)
  REAL :: wbar  (nx,ny,nz)     ! Base state w-velocity (m/s)
  REAL :: ptbar (nx,ny,nz)     ! Base state potential temperature (K)
  REAL :: pbar  (nx,ny,nz)     ! Base state pressure (Pascal).
  REAL :: rhobar(nx,ny,nz)     ! Base state air density (kg/m**3)
  REAL :: qvbar (nx,ny,nz)     ! Base state water vapor mixing ratio

  INTEGER :: nstyps                    ! Number of soil type
  INTEGER :: soiltyp (nx,ny,nstyps)    ! Soil type
  REAL :: stypfrct(nx,ny,nstyps)       ! Fraction of soil types
  INTEGER :: vegtyp(nx,ny)             ! Vegetation type
  REAL :: lai    (nx,ny)        ! Leaf Area Index
  REAL :: roufns (nx,ny)        ! Surface roughness
  REAL :: veg    (nx,ny)        ! Vegetation fraction

  REAL :: tsoil (nx,ny,nzsoil,0:nstyps)  ! Soil temperature (K)
  REAL :: qsoil (nx,ny,nzsoil,0:nstyps)  ! Soil moisture (m**3/m**3)
  REAL :: wetcanp(nx,ny,0:nstyps)   ! Canopy water amount
  REAL :: snowdpth(nx,ny)           ! Snow depth (m)

  REAL :: raing(nx,ny)         ! Grid supersaturation rain
  REAL :: rainc(nx,ny)         ! Cumulus convective rain
  REAL :: prcrate(nx,ny,4)     ! precipitation rates (kg/(m**2*s))
                               ! prcrate(1,1,1) = total precip. rate
                               ! prcrate(1,1,2) = grid scale precip. rate
                               ! prcrate(1,1,3) = cumulative precip. rate
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

  INTEGER (KIND=INT16), ALLOCATABLE :: itmp(:,:,:)
  INTEGER (KIND=INT16), ALLOCATABLE :: itmpsoil(:,:,:,:)
                                                ! Temporary array
  REAL(SP), ALLOCATABLE :: hmax(:), hmin(:)         ! Temporary array
  REAL(SP), ALLOCATABLE :: hmaxsoil(:), hminsoil(:) ! Temporary array
!
! 06/28/2002  Zuwen He
!
! Create a tem2 which will be used when dump data in previous version.
! It may be better to pass the array from the argument list
! in the future.
!
  REAL, ALLOCATABLE :: tem2(:,:,:) ! Temporary array

  INTEGER :: ireturn

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
!-----------------------------------------------------------------------

  CHARACTER (LEN=40) :: fmtver410,fmtver500,fmtver530
  INTEGER  :: intver,intver410,intver500,intver530

  PARAMETER (fmtver410='004.10 HDF4 Coded Data',intver410=410)
  PARAMETER (fmtver500='005.00 HDF4 Coded Data',intver500=500)
  PARAMETER (fmtver530='005.30 HDF4 Coded Data',intver530=530)

  CHARACTER (LEN=40) :: fmtverin

  CHARACTER (LEN=10) :: tmunit

!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------

  INTEGER :: lchanl
  PARAMETER (lchanl=6)      ! Channel number for formatted printing.

  INTEGER :: i,j,k,is
  INTEGER :: nxin,nyin,nzin,nzsoilin

  INTEGER :: bgrdin,bbasin,bvarin,bicein,btrbin,btkein

  INTEGER :: istat, sd_id

  INTEGER :: nstyp1,nstypin

  INTEGER :: nztmp
  INTEGER :: nq, nqscalarin(nscalar)

  CHARACTER(LEN=256) :: cmntstr

!-----------------------------------------------------------------------
!
!  Functions
!
!-----------------------------------------------------------------------

  CHARACTER(LEN=40) :: upcase

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  IF(myproc == 0) &
  WRITE(*,*) 'HDFREAD: Reading HDF file: ', trim(filename)

!
! If "nz" is a small number, the arrays may not be big enough for the
! soil family of variables.
!
  nztmp = max(nz,nstyps+1)

  ALLOCATE (itmp(nx,ny,nztmp),stat=istat)
  CALL check_alloc_status(istat, "HDFREAD:itmp")
  ALLOCATE (hmax(nztmp),stat=istat)
  CALL check_alloc_status(istat, "HDFREAD:hmax")
  ALLOCATE (hmin(nztmp),stat=istat)
  CALL check_alloc_status(istat, "HDFREAD:hmin")
  ALLOCATE (itmpsoil(nx,ny,nzsoil,0:nstyps),stat=istat)
  CALL check_alloc_status(istat, "HDFREAD:itmpsoil")
  ALLOCATE (hmaxsoil(nzsoil),stat=istat)
  CALL check_alloc_status(istat, "HDFREAD:hmaxsoil")
  ALLOCATE (hminsoil(nzsoil),stat=istat)
  CALL check_alloc_status(istat, "HDFREAD:hminsoil")

!-----------------------------------------------------------------------
!
!  Read header info
!
!-----------------------------------------------------------------------

  CALL hdfopen(filename,1,sd_id)
  IF (sd_id < 0) THEN
    IF(myproc == 0) &
    WRITE (6,*) "HDFREAD: ERROR opening ",                              &
                 trim(filename)," for reading."
    GO TO 110
  END IF

  CALL hdfrdc(sd_id,40,"fmtver",fmtverin,istat)

  IF (TRIM(fmtverin) == fmtver410) THEN
    intver=intver410
  ELSE IF (TRIM(fmtverin) == fmtver500) THEN
    intver=intver500
  ELSE IF (TRIM(fmtverin) == fmtver530) THEN
    intver=intver530
  ELSE
    intver=intver500
!    IF (myproc == 0) WRITE(6,'(/1x,a,a,a/)')                        &
!        'Incoming data format, fmtverin=',fmtverin,                 &
!        ', not found. The Job stopped.'
!    CALL arpsstop('arpstop called from HDFREAD. ',1)
  END IF

  IF (myproc == 0) WRITE(6,'(/1x,a,a/)')                                &
                              'Incoming data format, fmtverin=',fmtverin

  CALL hdfrdc(sd_id,80,"runname",runname,istat)
  CALL hdfrdi(sd_id,"nocmnt",nocmnt,istat)
  IF( nocmnt > 0 ) THEN
    IF (80*nocmnt > 256) THEN
      WRITE(*,'(1x,a,I0,a)') "ERROR: comment array is ",80*nocmnt,      &
                           ", but comment str is only 256 characters."
      CALL arpsstop('ERROR:Comment string is too short.',1)
    END IF
    CALL hdfrdc(sd_id,80*nocmnt,"cmnt",cmntstr,istat)
    DO i = 1, nocmnt
      cmnt(i) = cmntstr( ((i-1)*80+1):(i*80) )
    END DO
  END IF

  IF(myproc == 0) THEN
    WRITE(6,'(//''  THE NAME OF THIS RUN IS:  '',A//)') trim(runname)

    WRITE (6,*) "Comments:"

    IF( nocmnt > 0 ) THEN
      DO i=1,nocmnt
        WRITE(6,'(1x,a)') cmnt(i)
      END DO
    END IF
    WRITE (6,*) " "

  END IF

  CALL hdfrdc(sd_id,10,"tmunit",tmunit,istat)
  CALL hdfrdr(sd_id,"time",time,istat)

!-----------------------------------------------------------------------
!
!  Get dimensions of data in binary file and check against
!  the dimensions passed to HDFREAD
!
!-----------------------------------------------------------------------

  CALL hdfrdi(sd_id,"nx",nxin,istat)
  CALL hdfrdi(sd_id,"ny",nyin,istat)
  CALL hdfrdi(sd_id,"nz",nzin,istat)

  IF ( nxin /= nx .OR. nyin /= ny .OR. nzin /= nz ) THEN
    WRITE(6,'(1x,a)') ' Dimensions in HDFREAD inconsistent with data.'
    WRITE(6,'(1x,a,3I15)') ' Read were: ', nxin, nyin, nzin
    WRITE(6,'(1x,a,3I15)') ' Expected:  ', nx, ny, nz
    WRITE(6,'(1x,a)') ' Program aborted in HDFREAD.'
    CALL arpsstop('arpsstop called from HDFREAD due to nxin...',1)
  END IF

  IF (intver <= intver410) THEN
    nzsoilin = 2   ! for versions earlier than 410, it is actually
                   ! a 2 level soil model.
  ELSE IF (intver >= intver500) THEN
    CALL hdfrdi(sd_id,"nzsoil",nzsoilin,istat)
  END IF

  IF (nzsoilin /= nzsoil) THEN

    IF (intver <= intver410) THEN

      IF(myproc == 0) WRITE(6,'(1x,a,a/,2(1x,a/))')           &
        ' The incoming data version is ', fmtverin,            &
        ' In the input file, nzsoil must be set to 2. ',       &
        ' Program aborted in HDFREAD.'

    ELSE IF (intver >= intver500) THEN

      WRITE(6,'(1x,a)') &
                ' Dimensions in HDFREAD inconsistent with data.'
      WRITE(6,'(1x,a,I15)') ' Read were: ', nzsoilin
      WRITE(6,'(1x,a,I15)') ' Expected:  ', nzsoil
      WRITE(6,'(1x,a)') ' Program aborted in HDFREAD.'
      WRITE(6,*)' myproc = ',myproc

    END IF
    CALL arpsstop('arpsstop called from HDFREAD due to nzsoilin...',1)
  END IF

!-----------------------------------------------------------------------
!
!  Read in flags for different data groups
!
!-----------------------------------------------------------------------

  IF ( grdbas == 1 ) THEN   ! Read grid and base state arrays

    IF(myproc == 0) &
    WRITE(lchanl,'(1x,a,f8.1,a,f8.3,a/)')                               &
         'To read grid and base state data at time ', time,             &
         ' secs = ',(time/60.),' mins.'

    CALL hdfrdi(sd_id,"grdflg",bgrdin,istat)
    CALL hdfrdi(sd_id,"basflg",bbasin,istat)
    CALL hdfrdi(sd_id,"varflg",bvarin,istat)
    CALL hdfrdi(sd_id,"mstflg",mstin,istat)
    CALL hdfrdi(sd_id,"iceflg",bicein,istat)
    CALL hdfrdi(sd_id,"trbflg",btrbin,istat)
    CALL hdfrdi(sd_id,"landflg",landin,istat)
    CALL hdfrdi(sd_id,"totflg",totin,istat)
    CALL hdfrdi(sd_id,"tkeflg",btkein,istat)

  ELSE ! Normal data reading

    IF(myproc == 0) &
    WRITE(lchanl,'(1x,a,f8.1,a,f8.3,a/)')'To read data for time:',      &
         time,' secs = ',(time/60.),' mins.'

    CALL hdfrdi(sd_id,"grdflg",grdin,istat)
    CALL hdfrdi(sd_id,"basflg",basin,istat)
    CALL hdfrdi(sd_id,"varflg",varin,istat)
    CALL hdfrdi(sd_id,"mstflg",mstin,istat)
    CALL hdfrdi(sd_id,"iceflg",icein,istat)
    CALL hdfrdi(sd_id,"trbflg",trbin,istat)
    CALL hdfrdi(sd_id,"sfcflg",sfcin,istat)
    CALL hdfrdi(sd_id,"rainflg",rainin,istat)
    !Don't read landin (land data only written to base file,
    !landin flag set there).
    !CALL hdfrdi(sd_id,"landflg",landin,istat)
    CALL hdfrdi(sd_id,"totflg",totin,istat)
    CALL hdfrdi(sd_id,"tkeflg",tkein,istat)

    IF (intver >= intver530) THEN       ! New version

      CALL hdfrdi(sd_id,"nscalar",nscalarin,istat)
      DO nq = 1,nscalar
        CALL hdfrdi(sd_id,'P_'//upcase(qnames(nq)),nqscalarin(nq),istat)
      END DO

    ELSE                        ! maybe old version
      nscalarin = 0
      nqscalarin(:) = 0

      IF (mstin == 1) THEN

        IF (P_QC > 0) THEN
           nqscalarin(P_QC) = 1
           nscalarin = nscalarin + 1
        END IF
        IF (P_QR > 0) THEN
           nqscalarin(P_QR) = 2
           nscalarin = nscalarin + 1
        END IF

        IF (icein == 1) THEN
          IF (P_QI > 0) THEN
             nqscalarin(P_QI) = 3
             nscalarin = nscalarin + 1
          END IF
          IF (P_QS > 0) THEN
             nqscalarin(P_QS) = 4
             nscalarin = nscalarin + 1
          END IF
          IF (P_QH > 0) THEN
             nqscalarin(P_QH) = 5
             nscalarin = nscalarin + 1
          END IF
        END IF

      END IF

    END IF                   ! new version or old version

    IF(myproc == 0) WRITE(6,*) 'nscalarin = ',nscalarin
    DO nq = 1,nscalar
      IF(myproc == 0) THEN
        WRITE(6,FMT='(1x,a,I3,a)',ADVANCE='NO')           &
                  TRIM(qnames(nq))//'in = ',nqscalarin(nq), ','
        IF (MOD(nq,6) == 0) WRITE(6,*)
      END IF
    END DO
    IF(myproc == 0) WRITE(6,'(/)')

  END IF

  CALL hdfrdi(sd_id,"nstyp",nstyp1,istat)

  IF ( nstyp1 < 1 ) THEN
    nstyp1 = 1
  END IF
  IF (nstyp1 > nstyp) THEN
    IF(myproc == 0) WRITE (6,'(1x,a,I2,a,/,17x,2(a,I2,a))')             &
       'HDFREAD: INFO - nstyp in data file (',nstyp1,') is',            &
       'greater than that specified in namelist (',nstyp,') ',          &
       'using only nstyp = ',nstyp,'.'
    nstypin = nstyp
  ELSE
    nstypin = nstyp1
  ENDIF

  CALL hdfrdi(sd_id,"prcflg",prcin,istat)
  CALL hdfrdi(sd_id,"radflg",radin,istat)
  CALL hdfrdi(sd_id,"flxflg",flxin,istat)
  CALL hdfrdi(sd_id,"snowflg",snowin,istat)

  CALL hdfrdi(sd_id,"month",month,istat)
  CALL hdfrdi(sd_id,"day",day,istat)
  CALL hdfrdi(sd_id,"year",year,istat)
  CALL hdfrdi(sd_id,"hour",hour,istat)
  CALL hdfrdi(sd_id,"minute",minute,istat)
  CALL hdfrdi(sd_id,"second",second,istat)

  CALL hdfrdr(sd_id,"umove",umove,istat)
  CALL hdfrdr(sd_id,"vmove",vmove,istat)
  CALL hdfrdr(sd_id,"xgrdorg",xgrdorg,istat)
  CALL hdfrdr(sd_id,"ygrdorg",ygrdorg,istat)

  CALL hdfrdi(sd_id,"mapproj",mapproj,istat)
  CALL hdfrdr(sd_id,"trulat1",trulat1,istat)
  CALL hdfrdr(sd_id,"trulat2",trulat2,istat)
  CALL hdfrdr(sd_id,"trulon",trulon,istat)
  CALL hdfrdr(sd_id,"sclfct",sclfct,istat)
  CALL hdfrdr(sd_id,"tstop",tstop,istat)
  CALL hdfrdr(sd_id,"thisdmp",thisdmp,istat)
  CALL hdfrdr(sd_id,"latitud",latitud,istat)
  CALL hdfrdr(sd_id,"ctrlat",ctrlat,istat)
  CALL hdfrdr(sd_id,"ctrlon",ctrlon,istat)

  CALL hdfrdr(sd_id,"n0rain", n0rain,istat)
  CALL hdfrdr(sd_id,"n0snow", n0snow,istat)
  CALL hdfrdr(sd_id,"n0hail", n0hail,istat)
  CALL hdfrdr(sd_id,"rhosnow",rhosnow,istat)
  CALL hdfrdr(sd_id,"rhohail",rhohail,istat)

  IF (intver >= intver530) THEN
    CALL hdfrdr(sd_id,'ntcloud',ntcloud,istat)
    CALL hdfrdr(sd_id,'n0grpl', n0grpl,istat)
    CALL hdfrdr(sd_id,'rhoice', rhoice,istat)
    CALL hdfrdr(sd_id,'rhogrpl',rhogrpl,istat)
    CALL hdfrdr(sd_id,'alpharain',alpharain,istat)
    CALL hdfrdr(sd_id,'alphaice', alphaice, istat)
    CALL hdfrdr(sd_id,'alphasnow',alphasnow,istat)
    CALL hdfrdr(sd_id,'alphagrpl',alphagrpl,istat)
    CALL hdfrdr(sd_id,'alphahail',alphahail,istat)
  END IF
!-----------------------------------------------------------------------
!
!  Read in x,y and z at grid cell centers (scalar points).
!
!-----------------------------------------------------------------------

  IF( grdin == 1 .OR. grdbas == 1 ) THEN

    CALL hdfrd1d(sd_id,"x",nx,x,istat)
    IF (istat /= 0) GO TO 110

    CALL hdfrd1d(sd_id,"y",ny,y,istat)
    IF (istat /= 0) GO TO 110

    CALL hdfrd1d(sd_id,"z",nz,z,istat)
    IF (istat /= 0) GO TO 110

    CALL hdfrd3d(sd_id,"zp",nx,ny,nz,zp,istat,itmp,hmax,hmin)
    IF (istat /= 0) GO TO 110

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
        WRITE(6,*) ' Assign zpsoil. '
        WRITE(6,*) ' Assume zpsoil(,,1) is zp(,,2). '
        WRITE(6,*) ' Assume zpsoil(,,2) is zp(,,2)-1. '
      END IF

    ELSE IF (intver >= intver500) THEN
      CALL hdfrd3d(sd_id,"zpsoil",nx,ny,nzsoil,zpsoil,istat, &
                 itmpsoil,hmaxsoil,hminsoil)
    END IF

    IF (istat /= 0) GO TO 110

  END IF  ! grdin

!-----------------------------------------------------------------------
!
!  Read in base state fields
!
!-----------------------------------------------------------------------

  IF( basin == 1 .OR. grdbas == 1 ) THEN

    CALL hdfrd3d(sd_id,"ubar",nx,ny,nz,ubar,istat,itmp,hmax,hmin)
    IF (istat /= 0) GO TO 110

    CALL hdfrd3d(sd_id,"vbar",nx,ny,nz,vbar,istat,itmp,hmax,hmin)
    IF (istat /= 0) GO TO 110

    CALL hdfrd3d(sd_id,"wbar",nx,ny,nz,wbar,istat,itmp,hmax,hmin)
    IF (istat /= 0) GO TO 110

    CALL hdfrd3d(sd_id,"ptbar",nx,ny,nz,ptbar,istat,itmp,hmax,hmin)
    IF (istat /= 0) GO TO 110

    CALL hdfrd3d(sd_id,"pbar",nx,ny,nz,pbar,istat,itmp,hmax,hmin)
    IF (istat /= 0) GO TO 110

    IF( mstin == 1) THEN
      CALL hdfrd3d(sd_id,"qvbar",nx,ny,nz,qvbar,istat,itmp,hmax,hmin)
      IF (istat /= 0) GO TO 110
    END IF

    IF (landin == 1) THEN

      CALL hdfrd3di(sd_id,"soiltyp",nx,ny,nstypin,soiltyp(1,1,1),istat)
      IF (istat /= 0) GO TO 110

      CALL hdfrd3d(sd_id,"stypfrct",nx,ny,nstypin,                      &
          stypfrct(1,1,1),istat,itmp,hmax,hmin)
      IF (istat /= 0) GO TO 110

      CALL fix_stypfrct_nstyp(nx,ny,nstyp1,nstyp,stypfrct)

      CALL hdfrd2di(sd_id,"vegtyp",nx,ny,vegtyp,istat)
      IF (istat /= 0) GO TO 110

      CALL hdfrd2d(sd_id,"lai",nx,ny,lai,istat,itmp)
      IF (istat /= 0) GO TO 110

      CALL hdfrd2d(sd_id,"roufns",nx,ny,roufns,istat,itmp)
      IF (istat /= 0) GO TO 110

      CALL hdfrd2d(sd_id,"veg",nx,ny,veg,istat,itmp)
      IF (istat /= 0) GO TO 110

    END IF

  END IF

  IF( grdbas == 1 ) GO TO 930

  IF( varin == 1 ) THEN

    IF ( totin == 0 ) THEN

!-----------------------------------------------------------------------
!
!  Read in perturbations from history dump
!
!-----------------------------------------------------------------------

      CALL hdfrd3d(sd_id,"uprt",nx,ny,nz,uprt,istat,itmp,hmax,hmin)
      IF (istat /= 0) GO TO 110

      CALL hdfrd3d(sd_id,"vprt",nx,ny,nz,vprt,istat,itmp,hmax,hmin)
      IF (istat /= 0) GO TO 110

      CALL hdfrd3d(sd_id,"wprt",nx,ny,nz,wprt,istat,itmp,hmax,hmin)
      IF (istat /= 0) GO TO 110

!-----------------------------------------------------------------------
!
!  Read in scalars
!
!-----------------------------------------------------------------------

      CALL hdfrd3d(sd_id,"ptprt",nx,ny,nz,ptprt,istat,itmp,hmax,hmin)
      IF (istat /= 0) GO TO 110

      CALL hdfrd3d(sd_id,"pprt",nx,ny,nz,pprt,istat,itmp,hmax,hmin)
      IF (istat /= 0) GO TO 110

    ELSE

!-----------------------------------------------------------------------
!
!  Read in total values of variables from history dump
!
!-----------------------------------------------------------------------

      CALL hdfrd3d(sd_id,"u",nx,ny,nz,uprt,istat,itmp,hmax,hmin)
      IF (istat /= 0) GO TO 110
      DO k=1,nz-1
        DO j=1,ny-1
          DO i=1,nx
            uprt(i,j,k) = uprt(i,j,k) - ubar(i,j,k)
          END DO
        END DO
      END DO

      CALL hdfrd3d(sd_id,"v",nx,ny,nz,vprt,istat,itmp,hmax,hmin)
      IF (istat /= 0) GO TO 110
      DO k=1,nz-1
        DO j=1,ny
          DO i=1,nx-1
            vprt(i,j,k) = vprt(i,j,k) - vbar(i,j,k)
          END DO
        END DO
      END DO

      CALL hdfrd3d(sd_id,"w",nx,ny,nz,wprt,istat,itmp,hmax,hmin)
      IF (istat /= 0) GO TO 110

      CALL hdfrd3d(sd_id,"pt",nx,ny,nz,ptprt,istat,itmp,hmax,hmin)
      IF (istat /= 0) GO TO 110
      DO k=1,nz-1
        DO j=1,ny-1
          DO i=1,nx-1
            ptprt(i,j,k) = ptprt(i,j,k) - ptbar(i,j,k)
          END DO
        END DO
      END DO

      CALL hdfrd3d(sd_id,"p",nx,ny,nz,pprt,istat,itmp,hmax,hmin)
      IF (istat /= 0) GO TO 110
      DO k=1,nz-1
        DO j=1,ny-1
          DO i=1,nx-1
            pprt(i,j,k) = pprt(i,j,k) - pbar(i,j,k)
          END DO
        END DO
      END DO

    END IF

  END IF

!-----------------------------------------------------------------------
!
!  Read in moisture variables
!
!-----------------------------------------------------------------------

  IF( mstin == 1 ) THEN

    IF ( totin == 0 ) THEN

      CALL hdfrd3d(sd_id,"qvprt",nx,ny,nz,qvprt,istat,itmp,hmax,hmin)
      IF (istat /= 0) GO TO 110

    ELSE

      CALL hdfrd3d(sd_id,"qv",nx,ny,nz,qvprt,istat,itmp,hmax,hmin)
      IF (istat /= 0) GO TO 110
      DO k=1,nz-1
        DO j=1,ny-1
          DO i=1,nx-1
            qvprt(i,j,k) = qvprt(i,j,k) - qvbar(i,j,k)
          END DO
        END DO
      END DO

    END IF

    IF (intver < intver530) THEN

      IF (P_QC > 0 )   &
        CALL hdfrd3d(sd_id,"qc",nx,ny,nz,qscalar(:,:,:,P_QC),istat,itmp,hmax,hmin)
      IF (istat /= 0) GO TO 110

      IF (P_QR > 0 )   &
        CALL hdfrd3d(sd_id,"qr",nx,ny,nz,qscalar(:,:,:,P_QR),istat,itmp,hmax,hmin)
      IF (istat /= 0) GO TO 110

    ELSE

      DO nq = 1,nscalar
        IF (nqscalarin(nq) > 0 )  THEN
          CALL hdfrd3d(sd_id,qnames(nq),nx,ny,nz,qscalar(:,:,:,nq),istat,itmp,hmax,hmin)
          IF (istat /= 0) GO TO 110
        END IF
      END DO

    END IF

    IF( rainin == 1 ) THEN
      CALL hdfrd2d(sd_id,"raing",nx,ny,raing,istat,itmp)
      IF (istat /= 0) GO TO 110

      CALL hdfrd2d(sd_id,"rainc",nx,ny,rainc,istat,itmp)
      IF (istat /= 0) GO TO 110
    END IF

    IF( prcin == 1 ) THEN
      CALL hdfrd2d(sd_id,"prcrate1",nx,ny,prcrate(1,1,1),istat,itmp)
      IF (istat /= 0) GO TO 110

      CALL hdfrd2d(sd_id,"prcrate2",nx,ny,prcrate(1,1,2),istat,itmp)
      IF (istat /= 0) GO TO 110

      CALL hdfrd2d(sd_id,"prcrate3",nx,ny,prcrate(1,1,3),istat,itmp)
      IF (istat /= 0) GO TO 110

      CALL hdfrd2d(sd_id,"prcrate4",nx,ny,prcrate(1,1,4),istat,itmp)
      IF (istat /= 0) GO TO 110
    END IF

    IF( icein == 1 .AND. intver < intver530) THEN

      IF (P_QI > 0 )   &
        CALL hdfrd3d(sd_id,"qi",nx,ny,nz,qscalar(:,:,:,P_QI),istat,itmp,hmax,hmin)
      IF (istat /= 0) GO TO 110

      IF (P_QS > 0 )   &
        CALL hdfrd3d(sd_id,"qs",nx,ny,nz,qscalar(:,:,:,P_QS),istat,itmp,hmax,hmin)
      IF (istat /= 0) GO TO 110

      IF (P_QH > 0 )   &
        CALL hdfrd3d(sd_id,"qh",nx,ny,nz,qscalar(:,:,:,P_QH),istat,itmp,hmax,hmin)
      IF (istat /= 0) GO TO 110

    END IF
  END IF

  IF( tkein == 1 ) THEN

    CALL hdfrd3d(sd_id,"tke",nx,ny,nz,tke,istat,itmp,hmax,hmin)
    IF (istat /= 0) GO TO 110

  END IF

  IF( trbin == 1 ) THEN

    CALL hdfrd3d(sd_id,"kmh",nx,ny,nz,kmh,istat,itmp,hmax,hmin)
    IF (istat /= 0) GO TO 110

    CALL hdfrd3d(sd_id,"kmv",nx,ny,nz,kmv,istat,itmp,hmax,hmin)
    IF (istat /= 0) GO TO 110

  END IF

  IF( sfcin == 1 ) THEN

    IF (intver <= intver410) THEN

      ALLOCATE (tem2(nx,ny,0:nstyps),stat=istat)

      tem2=0.
      CALL hdfrd3d(sd_id,"tsfc",nx,ny,nstypin+1,tem2,istat,itmp,hmax,hmin)
      IF (istat /= 0) GO TO 110

      DO is=0,nstypin
        DO j=1,ny
          DO i=1,nx
            tsoil(i,j,1,is)=tem2(i,j,is)
          END DO
        END DO
      END DO

      tem2=0.
      CALL hdfrd3d(sd_id,"tsoil",nx,ny,nstypin+1,tem2,istat,itmp,hmax,hmin)
      IF (istat /= 0) GO TO 110

      DO is=0,nstypin
        DO j=1,ny
          DO i=1,nx
            tsoil(i,j,2,is)=tem2(i,j,is)
          END DO
        END DO
      END DO

      tem2=0.
      CALL hdfrd3d(sd_id,"wetsfc",nx,ny,nstypin+1,tem2,istat,itmp,hmax,hmin)
      IF (istat /= 0) GO TO 110

      DO is=0,nstypin
        DO j=1,ny
          DO i=1,nx
            qsoil(i,j,1,is)=tem2(i,j,is)
          END DO
        END DO
      END DO

      tem2=0.
      CALL hdfrd3d(sd_id,"wetdp",nx,ny,nstypin+1,tem2,istat,itmp,hmax,hmin)
      IF (istat /= 0) GO TO 110

      DO is=0,nstypin
        DO j=1,ny
          DO i=1,nx
            qsoil(i,j,2,is)=tem2(i,j,is)
          END DO
        END DO
      END DO

      DEALLOCATE (tem2,stat=istat)

    ELSE IF (intver >= intver500) THEN

      CALL hdfrd4d(sd_id,"tsoil",nx,ny,nzsoil,nstypin+1,tsoil,istat, &
                       itmpsoil,hmaxsoil,hminsoil)
      IF (istat /= 0) GO TO 110

      CALL hdfrd4d(sd_id,"qsoil",nx,ny,nzsoil,nstypin+1,qsoil,istat, &
                       itmpsoil,hmaxsoil,hminsoil)
      IF (istat /= 0) GO TO 110

    END IF

    CALL hdfrd3d(sd_id,"wetcanp",nx,ny,nstypin+1,wetcanp,istat,itmp,hmax,hmin)
    IF (istat /= 0) GO TO 110

    CALL fix_soil_nstyp(nx,ny,nzsoil,nstyp1,nstyp,tsoil,qsoil,wetcanp)

    IF (snowin == 1) THEN
      CALL hdfrd2d(sd_id,"snowdpth",nx,ny,snowdpth,istat,itmp)
      IF (istat /= 0) GO TO 110
    END IF

  END IF

  IF( radin == 1 ) THEN

    CALL hdfrd3d(sd_id,"radfrc",nx,ny,nz,radfrc,istat,                   &
                   itmp,hmax,hmin)
    IF (istat /= 0) GO TO 110

    CALL hdfrd2d(sd_id,"radsw",nx,ny,radsw,istat,itmp)
    IF (istat /= 0) GO TO 110

    CALL hdfrd2d(sd_id,"rnflx",nx,ny,rnflx,istat,itmp)
    IF (istat /= 0) GO TO 110

    IF (intver <= intver410) THEN

      radswnet=0.
      radlwin=0.

    ELSE IF (intver >= intver500) THEN

      CALL hdfrd2d(sd_id,"radswnet",nx,ny,radswnet,istat,itmp)
      IF (istat /= 0) GO TO 110

      CALL hdfrd2d(sd_id,"radlwin",nx,ny,radlwin,istat,itmp)
      IF (istat /= 0) GO TO 110

    END IF

  END IF

  IF( flxin == 1 ) THEN

    CALL hdfrd2d(sd_id,"usflx",nx,ny,usflx,istat,itmp)
    IF (istat /= 0) GO TO 110

    CALL hdfrd2d(sd_id,"vsflx",nx,ny,vsflx,istat,itmp)
    IF (istat /= 0) GO TO 110

    CALL hdfrd2d(sd_id,"ptsflx",nx,ny,ptsflx,istat,itmp)
    IF (istat /= 0) GO TO 110

    CALL hdfrd2d(sd_id,"qvsflx",nx,ny,qvsflx,istat,itmp)
    IF (istat /= 0) GO TO 110

  END IF

!-----------------------------------------------------------------------
!
!  Friendly exit message
!
!-----------------------------------------------------------------------

  930   CONTINUE

  IF(myproc == 0) &
  WRITE(6,'(/a,F8.1,a/)')                                               &
      ' Data at time=', time/60,' (min) were successfully read.'

  ireturn = 0

  GO TO 130

!-----------------------------------------------------------------------
!
!  Error during read
!
!-----------------------------------------------------------------------

  110   CONTINUE
  IF(myproc == 0) &
  WRITE(6,'(/a/)') ' Error reading data in HDFREAD'
  ireturn=1

  GO TO 130

!-----------------------------------------------------------------------
!
!  End-of-file during read
!
!-----------------------------------------------------------------------

!  120   CONTINUE
!  IF(myproc == 0) &
!  WRITE(6,'(/a/)') ' End of file reached in HDFREAD'
!  ireturn=2

  130   CONTINUE

!tmp      istat = sfendacc(sd_id)   ! is this necessary?
  CALL hdfclose(sd_id,istat)
  IF (ireturn == 0) THEN
    IF (istat == 0) THEN
      IF(myproc == 0) &
      WRITE(*,*) "HDFREAD: Successfully read ", trim(filename)
    ELSE
      IF(myproc == 0) &
      WRITE(*,*) "HDFREAD: ERROR (status=", istat, ") closing ", trim(filename)
    END IF
  END IF

  DEALLOCATE (itmp,stat=istat)
  DEALLOCATE (hmax,stat=istat)
  DEALLOCATE (hmin,stat=istat)
  DEALLOCATE (itmpsoil,stat=istat)
  DEALLOCATE (hmaxsoil,stat=istat)
  DEALLOCATE (hminsoil,stat=istat)

  RETURN
END SUBROUTINE hdfread

!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE HDFREADSPLIT               ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################

SUBROUTINE hdfreadsplit(nx,ny,nz,nzsoil,nstyps, grdbas, filename, time, &
           x,y,z,zp,zpsoil,                                             &
           uprt, vprt, wprt, ptprt, pprt,                               &
           qvprt, qscalar, tke, kmh,kmv,                                &
           ubar, vbar, wbar, ptbar, pbar, rhobar, qvbar,                &
           soiltyp,stypfrct,vegtyp,lai,roufns,veg,                      &
           tsoil,qsoil,wetcanp,snowdpth,                                &
           raing,rainc,prcrate,                                         &
           radfrc,radsw,rnflx,radswnet,radlwin,                         &
           usflx,vsflx,ptsflx,qvsflx,                                   &
           ireturn, tem1)

!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Read in history data in the NCSA HDF4 format, split and scatter to
!  each processors using message passing.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yunheng Wang
!  2002/08/26
!  Based on subroutine hdfread.
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
!    filename  Character variable nhming the input HDF file
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
!    wprt     Vertical component of perturbation velocity in Cartesian
!             coordinates (m/s).
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
!
!    WORK ARRAY
!
!    tem1
!    tem2     work array.
!
!-----------------------------------------------------------------------

  USE arps_precision

  IMPLICIT NONE

!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------

  INCLUDE 'globcst.inc'
  INCLUDE 'grid.inc'        ! Grid parameters
  INCLUDE 'indtflg.inc'
  INCLUDE 'alloc.inc'       ! allocation parameters & declarations
  INCLUDE 'phycst.inc'
  INCLUDE 'mp.inc'

!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------

  INTEGER :: nx,ny,nz
  INTEGER :: nzsoil

  INTEGER :: grdbas
  CHARACTER (LEN=*) :: filename

  REAL :: x     (nx)           ! x coord.
  REAL :: y     (ny)           ! y coord.
  REAL :: z     (nz)           ! z coord.
  REAL :: zp    (nx,ny,nz)     ! physical x coord.
  REAL :: zpsoil(nx,ny,nzsoil) ! physical x coord. for soil (m)

  REAL :: uprt  (nx,ny,nz)     ! Perturbation u-velocity (m/s)
  REAL :: vprt  (nx,ny,nz)     ! Perturbation v-velocity (m/s)
  REAL :: wprt  (nx,ny,nz)     ! Perturbation w-velocity (m/s)
  REAL :: ptprt (nx,ny,nz)     ! Perturbation potential temperature (K)
  REAL :: pprt  (nx,ny,nz)     ! Perturbation pressure (Pascal)

  REAL :: qvprt (nx,ny,nz)     ! Perturbation water vapor mixing ratio (kg/kg)

  REAL :: qscalar    (nx,ny,nz,nscalar)

  REAL :: tke   (nx,ny,nz)     ! Turbulent Kinetic Energy ((m/s)**2)
  REAL :: kmh   (nx,ny,nz)     ! Horizontal turb. mixing coef. for
                               ! momentum. ( m**2/s )
  REAL :: kmv   (nx,ny,nz)     ! Vertical turb. mixing coef. for
                               ! momentum. ( m**2/s )

  REAL :: ubar  (nx,ny,nz)     ! Base state u-velocity (m/s)
  REAL :: vbar  (nx,ny,nz)     ! Base state v-velocity (m/s)
  REAL :: wbar  (nx,ny,nz)     ! Base state w-velocity (m/s)
  REAL :: ptbar (nx,ny,nz)     ! Base state potential temperature (K)
  REAL :: pbar  (nx,ny,nz)     ! Base state pressure (Pascal).
  REAL :: rhobar(nx,ny,nz)     ! Base state air density (kg/m**3)
  REAL :: qvbar (nx,ny,nz)     ! Base state water vapor mixing ratio

  INTEGER :: nstyps                    ! Number of soil type
  INTEGER :: soiltyp (nx,ny,nstyps)    ! Soil type
  REAL :: stypfrct(nx,ny,nstyps)       ! Fraction of soil types
  INTEGER :: vegtyp(nx,ny)             ! Vegetation type
  REAL :: lai    (nx,ny)        ! Leaf Area Index
  REAL :: roufns (nx,ny)        ! Surface roughness
  REAL :: veg    (nx,ny)        ! Vegetation fraction

  REAL :: tsoil (nx,ny,nzsoil,0:nstyps)  ! Soil temperature (K)
  REAL :: qsoil (nx,ny,nzsoil,0:nstyps)  ! Soil moisture (m**3/m**3)
  REAL :: wetcanp(nx,ny,0:nstyps)   ! Canopy water amount
  REAL :: snowdpth(nx,ny)           ! Snow depth (m)

  REAL :: raing(nx,ny)         ! Grid supersaturation rain
  REAL :: rainc(nx,ny)         ! Cumulus convective rain
  REAL :: prcrate(nx,ny,4)     ! precipitation rates (kg/(m**2*s))
                               ! prcrate(1,1,1) = total precip. rate
                               ! prcrate(1,1,2) = grid scale precip. rate
                               ! prcrate(1,1,3) = cumulative precip. rate
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

  INTEGER (KIND=INT16),  ALLOCATABLE :: itmp(:,:,:)              ! Temporary array
  REAL(SP),              ALLOCATABLE :: hmax(:), hmin(:)         ! Temporary array
  REAL(SP),              ALLOCATABLE :: hmaxsoil(:), hminsoil(:) ! Temporary array
  INTEGER (KIND=INT16),  ALLOCATABLE :: itmpsoil(:,:,:,:)
!
!
! tem2 will be used when read data in previous version.
! It may be better to pass the array from the argument list
! in the future.
!
  REAL, allocatable :: tem2(:,:,:)              ! Temporary array

  INTEGER :: ireturn
  REAL    :: time

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
!-----------------------------------------------------------------------

  CHARACTER (LEN=40) :: fmtver410,fmtver500,fmtver530
  INTEGER  :: intver,intver410,intver500,intver530

  PARAMETER (fmtver410='004.10 HDF4 Coded Data',intver410=410)
  PARAMETER (fmtver500='005.00 HDF4 Coded Data',intver500=500)
  PARAMETER (fmtver530='005.30 HDF4 Coded Data',intver530=530)

  CHARACTER (LEN=40) :: fmtverin

  CHARACTER (LEN=10) :: tmunit

!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------

  INTEGER :: lchanl
  PARAMETER (lchanl=6)      ! Channel number for formatted printing.

  INTEGER :: i,j,k,is
  INTEGER :: nxin,nyin,nzin,nzsoilin

  INTEGER :: bgrdin,bbasin,bvarin,bicein,btrbin,btkein

  INTEGER :: istat, sd_id

  INTEGER :: nstyp1,nstypin

  INTEGER :: nq, nqscalarin(nscalar)

  INTEGER :: nxlg, nylg, nzlg, nzsoillg
  REAL,    ALLOCATABLE :: var1d(:), var3d(:,:,:), var4d(:,:,:,:)
  INTEGER, ALLOCATABLE :: var3di(:,:,:)
                       ! temporary variable only needed for processor 0

  INTEGER :: nztmp

  CHARACTER(LEN=256) :: cmntstr

!-----------------------------------------------------------------------
!
!  Functions
!
!-----------------------------------------------------------------------

  CHARACTER(LEN=40) :: upcase

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  nxlg = nproc_x*(nx-3)+3
  nylg = nproc_y*(ny-3)+3
  nzlg = nz
  nzsoillg = nzsoil

!
! If "nz" is a small number, the arrays may not be big enough for the
! soil family of variables.
!
  nztmp = max(nzlg,nstyps+1)
  ALLOCATE (itmp(nxlg,nylg,nzlg),stat=istat)
  CALL check_alloc_status(istat, "hdfreadsplit:itmp")

  ALLOCATE (hmax(nzlg),stat=istat)
  CALL check_alloc_status(istat, "hdfreadsplit:hmax")

  ALLOCATE (hmin(nzlg),stat=istat)
  CALL check_alloc_status(istat, "hdfreadsplit:hmin")

  ALLOCATE (itmpsoil(nxlg,nylg,nzsoillg,0:nstyps),stat=istat)
  CALL check_alloc_status(istat, "hdfreadsplit:itmpsoil")

  ALLOCATE (hmaxsoil(nzsoillg),stat=istat)
  CALL check_alloc_status(istat, "hdfreadsplit:hmaxsoil")

  ALLOCATE (hminsoil(nzsoillg),stat=istat)
  CALL check_alloc_status(istat, "hdfreadsplit:hminsoil")

  IF(myproc == 0) THEN

    WRITE(*,*) 'HDFREADSPLIT: Reading HDF file: ', trim(filename)

!-----------------------------------------------------------------------
!
!  Read header info
!
!-----------------------------------------------------------------------

    CALL hdfopen(filename,1,sd_id)
    IF (sd_id < 0) THEN
      WRITE (6,*) "HDFREADSPLIT: ERROR opening ",                     &
                   trim(filename)," for reading."
      GO TO 110
    END IF

    CALL hdfrdc(sd_id,40,"fmtver",fmtverin,istat)

    IF (TRIM(fmtverin) == fmtver410) THEN
      intver=intver410
    ELSE IF (TRIM(fmtverin) == fmtver500) THEN
      intver=intver500
    ELSE IF (TRIM(fmtverin) == fmtver530) THEN
      intver=intver530
    ELSE
      WRITE(6,'(/1x,a,a,a/)')                                           &
          'Incoming data format, fmtverin=',fmtverin,                   &
          ', not found. The Job stopped.'
      CALL arpsstop('arpstop called from HDFREADSPLIT. ',1)
    END IF

    WRITE(6,'(/1x,a,a/)')                                               &
        'Incoming data format, fmtverin=',fmtverin

    CALL hdfrdc(sd_id,80,"runname",runname,istat)
    CALL hdfrdi(sd_id,"nocmnt",nocmnt,istat)
    IF( nocmnt > 0 ) THEN
      IF (80*nocmnt > 256) THEN
        WRITE(*,'(1x,a,I0,a)') "ERROR: comment array is ",80*nocmnt,    &
                             ", but comment str is only 256 characters."
        CALL arpsstop('ERROR:Comment string is too short.',1)
      END IF
      CALL hdfrdc(sd_id,80*nocmnt,"cmnt",cmntstr,istat)
      DO i = 1, nocmnt
        cmnt(i) = cmntstr((i-1)*80+1:i*80)
      END DO
    END IF

    WRITE(6,'(//''  THE NAME OF THIS RUN IS:  '',A//)') trim(runname)

    WRITE (6,*) "Comments:"
    IF( nocmnt > 0 ) THEN
      DO i=1,nocmnt
        WRITE(6,'(1x,a)') cmnt(i)
      END DO
    END IF
    WRITE (6,*) " "

    CALL hdfrdc(sd_id,10,"tmunit",tmunit,istat)
    CALL hdfrdr(sd_id,"time",time,istat)

  END IF ! myproc == 0

  CALL mpupdatei(intver,  1)
  CALL mpupdatec(runname, 80)
  CALL mpupdatec(tmunit,  10)
  CALL mpupdater(time,    1)

!-----------------------------------------------------------------------
!
!  Get dimensions of data in binary file and check against
!  the dimensions passed to HDFREADSPLIT
!
!-----------------------------------------------------------------------

  IF (myproc == 0) THEN
    CALL hdfrdi(sd_id,"nx",nxin,istat)
    CALL hdfrdi(sd_id,"ny",nyin,istat)
    CALL hdfrdi(sd_id,"nz",nzin,istat)

    IF ( nxin /= nxlg .OR. nyin /= nylg .OR. nzin /= nzlg ) THEN
      WRITE(6,'(1x,a)') ' Dimensions in HDFREADSPLIT inconsistent with data.'
      WRITE(6,'(1x,a,3I15)') ' Read were: ', nxin, nyin, nzin
      WRITE(6,'(1x,a,3I15)') ' Expected:  ', nxlg, nylg, nzlg
      WRITE(6,'(1x,a)') ' Program aborted in HDFREADSPLIT.'
      CALL arpsstop('arpsstop called from HDFREADSPLIT due to nxin...',1)
    END IF

    IF (intver <= intver410) THEN
      nzsoilin = 2   ! for versions earlier than 410, it is actually
                     ! a 2 level soil model.
    ELSE IF (intver >= intver500) THEN
      CALL hdfrdi(sd_id,"nzsoil",nzsoilin,istat)
    END IF

    IF (nzsoilin /= nzsoillg) THEN

      IF (intver <= intver410) THEN

        WRITE(6,'(1x,a,a/,2(1x,a/))')                            &
          ' The incoming data version is ', fmtverin,            &
          ' In the input file, nzsoil must be set to 2. ',       &
          ' Program aborted in HDFREADSPLIT.'

      ELSE IF (intver >= intver500) THEN

        WRITE(6,'(1x,a)') &
                  ' Dimensions in HDFREAD inconsistent with data.'
        WRITE(6,'(1x,a,I15)') ' Read were: ', nzsoilin
        WRITE(6,'(1x,a,I15)') ' Expected:  ', nzsoillg
        WRITE(6,'(1x,a)') ' Program aborted in HDFREADSPLIT.'

      END IF
      CALL arpsstop('arpsstop called from HDFREADSPLIT due to nzsoilin...',1)
    END IF

  END IF  ! myproc == 0

!-----------------------------------------------------------------------
!
!  Read in flags for different data groups
!
!-----------------------------------------------------------------------

  IF ( grdbas == 1 ) THEN   ! Read grid and base state arrays

    IF(myproc == 0) THEN
      WRITE(lchanl,'(1x,a,f8.1,a,f8.3,a/)')                             &
         'To read grid and base state data at time ', time,             &
         ' secs = ',(time/60.),' mins.'

      CALL hdfrdi(sd_id,"grdflg",bgrdin,istat)
      CALL hdfrdi(sd_id,"basflg",bbasin,istat)
      CALL hdfrdi(sd_id,"varflg",bvarin,istat)
      CALL hdfrdi(sd_id,"mstflg",mstin,istat)
      CALL hdfrdi(sd_id,"iceflg",bicein,istat)
      CALL hdfrdi(sd_id,"trbflg",btrbin,istat)
      CALL hdfrdi(sd_id,"landflg",landin,istat)
      CALL hdfrdi(sd_id,"totflg",totin,istat)
      CALL hdfrdi(sd_id,"tkeflg",btkein,istat)
    END IF

    CALL mpupdatei(bgrdin, 1)
    CALL mpupdatei(bbasin, 1)
    CALL mpupdatei(bvarin, 1)
    CALL mpupdatei(mstin, 1)
    CALL mpupdatei(bicein, 1)
    CALL mpupdatei(btrbin, 1)
    CALL mpupdatei(landin, 1)
    CALL mpupdatei(totin, 1)
    CALL mpupdatei(btkein, 1)

  ELSE ! Normal data reading

    IF(myproc == 0) THEN
      WRITE(lchanl,'(1x,a,f8.1,a,f8.3,a/)')'To read data for time:',      &
           time,' secs = ',(time/60.),' mins.'

      CALL hdfrdi(sd_id,"grdflg",grdin,istat)
      CALL hdfrdi(sd_id,"basflg",basin,istat)
      CALL hdfrdi(sd_id,"varflg",varin,istat)
      CALL hdfrdi(sd_id,"mstflg",mstin,istat)
      CALL hdfrdi(sd_id,"iceflg",icein,istat)
      CALL hdfrdi(sd_id,"trbflg",trbin,istat)
      CALL hdfrdi(sd_id,"sfcflg",sfcin,istat)
      CALL hdfrdi(sd_id,"rainflg",rainin,istat)
      !Don't read landin (land data only written to base file,
      !landin flag set there).
      !CALL hdfrdi(sd_id,"landflg",landin,istat)
      CALL hdfrdi(sd_id,"totflg",totin,istat)
      CALL hdfrdi(sd_id,"tkeflg",tkein,istat)

      IF (intver >= intver530) THEN       ! New version

        CALL hdfrdi(sd_id,"nscalar",nscalarin,istat)
        DO nq = 1,nscalar
          CALL hdfrdi(sd_id,'P_'//upcase(qnames(nq)),nqscalarin(nq),istat)
        END DO

      ELSE                        ! maybe old version
        nscalarin = 0
        nqscalarin(:) = 0

        IF (mstin == 1) THEN

          IF (P_QC > 0) THEN
            nqscalarin(P_QC) = 1
            nscalarin = nscalarin + 1
          END IF
          IF (P_QR > 0) THEN
            nqscalarin(P_QR) = 2
            nscalarin = nscalarin + 1
          END IF

          IF (icein == 1) THEN
            IF (P_QI > 0) THEN
              nqscalarin(P_QI) = 3
              nscalarin = nscalarin + 1
            END IF
            IF (P_QS > 0) THEN
              nqscalarin(P_QS) = 4
              nscalarin = nscalarin + 1
            END IF
            IF (P_QH > 0) THEN
              nqscalarin(P_QH) = 5
              nscalarin = nscalarin + 1
            END IF
          END IF

        END IF

      END IF

      WRITE(6,*) 'nscalarin = ',nscalarin
      DO nq = 1,nscalar
        WRITE(6,FMT='(1x,a,I3,a)',ADVANCE='NO') TRIM(qnames(nq))//'in = ', &
                       nqscalarin(nq),','
        IF (MOD(nq,6) == 0) WRITE(6,*)
      END DO
      WRITE(6,'(/)')

    END IF

    CALL mpupdatei(grdin, 1)
    CALL mpupdatei(basin, 1)
    CALL mpupdatei(varin, 1)
    CALL mpupdatei(mstin, 1)
    CALL mpupdatei(icein, 1)
    CALL mpupdatei(trbin, 1)
    CALL mpupdatei(sfcin, 1)
    CALL mpupdatei(rainin, 1)
    CALL mpupdatei(totin, 1)
    CALL mpupdatei(tkein, 1)

    CALL mpupdatei(nscalarin,  1)
    CALL mpupdatei(nqscalarin, nscalar)

  END IF

  IF(myproc == 0) THEN
    CALL hdfrdi(sd_id,"nstyp",nstyp1,istat)

    IF ( nstyp1 < 1 ) nstyp1 = 1

    IF (nstyp1 > nstyp) THEN
      WRITE (6,'(1x,a,I2,a,/,21x,a,I2,a,I2,a)')                         &
         'HDFREADSPLIT: INFO: nstyp in file (',nstyp1,                  &
         ') is greater than that specified','in namelist (',nstyp,   &
         '), using only ',nstyp,' stypes.'
      nstypin = nstyp
    ELSE
      nstypin = nstyp1
    ENDIF

    CALL hdfrdi(sd_id,"prcflg",prcin,istat)
    CALL hdfrdi(sd_id,"radflg",radin,istat)
    CALL hdfrdi(sd_id,"flxflg",flxin,istat)
    CALL hdfrdi(sd_id,"snowflg",snowin,istat)

    CALL hdfrdi(sd_id,"month",month,istat)
    CALL hdfrdi(sd_id,"day",day,istat)
    CALL hdfrdi(sd_id,"year",year,istat)
    CALL hdfrdi(sd_id,"hour",hour,istat)
    CALL hdfrdi(sd_id,"minute",minute,istat)
    CALL hdfrdi(sd_id,"second",second,istat)

    CALL hdfrdr(sd_id,"umove",umove,istat)
    CALL hdfrdr(sd_id,"vmove",vmove,istat)
    CALL hdfrdr(sd_id,"xgrdorg",xgrdorg,istat)
    CALL hdfrdr(sd_id,"ygrdorg",ygrdorg,istat)

    CALL hdfrdi(sd_id,"mapproj",mapproj,istat)
    CALL hdfrdr(sd_id,"trulat1",trulat1,istat)
    CALL hdfrdr(sd_id,"trulat2",trulat2,istat)
    CALL hdfrdr(sd_id,"trulon",trulon,istat)
    CALL hdfrdr(sd_id,"sclfct",sclfct,istat)
    CALL hdfrdr(sd_id,"tstop",tstop,istat)
    CALL hdfrdr(sd_id,"thisdmp",thisdmp,istat)
    CALL hdfrdr(sd_id,"latitud",latitud,istat)
    CALL hdfrdr(sd_id,"ctrlat",ctrlat,istat)
    CALL hdfrdr(sd_id,"ctrlon",ctrlon,istat)

    CALL hdfrdr(sd_id,"n0rain",n0rain,istat)
    CALL hdfrdr(sd_id,"n0snow",n0snow,istat)
    CALL hdfrdr(sd_id,"n0hail",n0hail,istat)
    CALL hdfrdr(sd_id,"rhosnow",rhosnow,istat)
    CALL hdfrdr(sd_id,"rhohail",rhohail,istat)
    IF (intver >= intver530) THEN
      CALL hdfrdr(sd_id,'ntcloud',ntcloud,istat)
      CALL hdfrdr(sd_id,'n0grpl', n0grpl,istat)
      CALL hdfrdr(sd_id,'rhoice', rhoice,istat)
      CALL hdfrdr(sd_id,'rhogrpl',rhogrpl,istat)
      CALL hdfrdr(sd_id,'alpharain',alpharain,istat)
      CALL hdfrdr(sd_id,'alphaice', alphaice, istat)
      CALL hdfrdr(sd_id,'alphasnow',alphasnow,istat)
      CALL hdfrdr(sd_id,'alphagrpl',alphagrpl,istat)
      CALL hdfrdr(sd_id,'alphahail',alphahail,istat)
    END IF
  END IF

  CALL mpupdatei(nstypin, 1)
  CALL mpupdatei(nstyp1, 1)

  CALL mpupdatei(prcin,1)
  CALL mpupdatei(radin,1)
  CALL mpupdatei(flxin,1)
  CALL mpupdatei(snowin,1)

  CALL mpupdatei(month,1)
  CALL mpupdatei(day,1)
  CALL mpupdatei(year,1)
  CALL mpupdatei(hour,1)
  CALL mpupdatei(minute,1)
  CALL mpupdatei(second,1)

  CALL mpupdater(umove,1)
  CALL mpupdater(vmove,1)
  CALL mpupdater(xgrdorg,1)
  CALL mpupdater(ygrdorg,1)

  CALL mpupdatei(mapproj,1)
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

  IF (intver >= intver530) THEN
    CALL mpupdater( ntcloud, 1)
    CALL mpupdater( n0grpl,  1)
    CALL mpupdater( rhoice,  1)
    CALL mpupdater( rhogrpl, 1)
    CALL mpupdater( alpharain, 1)
    CALL mpupdater( alphaice,  1)
    CALL mpupdater( alphasnow, 1)
    CALL mpupdater( alphagrpl, 1)
    CALL mpupdater( alphahail, 1)
  END IF

  ALLOCATE (var1d(MAX(nxlg,nylg,nzlg)),stat=istat)
  CALL check_alloc_status(istat, "hdfreadsplit:var1d")

  ALLOCATE (var3d(nxlg,nylg,MAX(nzlg,nzsoillg,nstyps+1,4)),stat=istat)
  CALL check_alloc_status(istat, "hdfreadsplit:var3d")

  ALLOCATE (var3di(nxlg,nylg,nstypin),stat=istat)
  CALL check_alloc_status(istat, "hdfreadsplit:var3di")

  ALLOCATE (var4d(nxlg,nylg,nzsoillg,0:nstyps),stat=istat)
  CALL check_alloc_status(istat, "hdfreadsplit:var4d")

!-----------------------------------------------------------------------
!
!  Read in x,y and z at grid cell centers (scalar points).
!
!-----------------------------------------------------------------------

  IF( grdin == 1 .OR. grdbas == 1 ) THEN

    IF(myproc == 0) THEN
      CALL hdfrd1d(sd_id,"x",nxlg,var1d,istat)
      IF (istat /= 0) GO TO 110
    END IF
    CALL mpisplit1dx(var1d,nx,x)

    IF(myproc == 0) THEN
      CALL hdfrd1d(sd_id,"y",nylg,var1d,istat)
      IF (istat /= 0) GO TO 110
    END IF
    CALL mpisplit1dy(var1d,ny,y)

    IF(myproc == 0) THEN
      CALL hdfrd1d(sd_id,"z",nzlg,z,istat)
      IF (istat /= 0) GO TO 110
    END IF
    CALL mpupdater(z, nz)

    IF(myproc == 0) THEN
      CALL hdfrd3d(sd_id,"zp",nxlg,nylg,nzlg,var3d,istat,itmp,hmax,hmin)
      IF (istat /= 0) GO TO 110
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
        WRITE(6,*) ' Assign zpsoil. '
        WRITE(6,*) ' Assume zpsoil(,,1) is zp(,,2). '
        WRITE(6,*) ' Assume zpsoil(,,2) is zp(,,2)-1. '
      END IF

    ELSE IF (intver >= intver500) THEN
      IF(myproc == 0) THEN
        CALL hdfrd3d(sd_id,"zpsoil",nxlg,nylg,nzsoillg,var3d,istat, &
                 itmpsoil,hmaxsoil,hminsoil)
        IF (istat /= 0) GO TO 110
      END IF
      CALL mpisplit3d(var3d,nx,ny,nzsoil,zpsoil)
    END IF


  END IF  ! grdin

!-----------------------------------------------------------------------
!
!  Read in base state fields
!
!-----------------------------------------------------------------------

  IF( basin == 1 .OR. grdbas == 1 ) THEN

    IF(myproc == 0) THEN
      CALL hdfrd3d(sd_id,"ubar",nxlg,nylg,nzlg,var3d,istat,itmp,hmax,hmin)
      IF (istat /= 0) GO TO 110
    END IF
    CALL mpisplit3d(var3d,nx,ny,nz,ubar)

    IF(myproc == 0) THEN
      CALL hdfrd3d(sd_id,"vbar",nxlg,nylg,nzlg,var3d,istat,itmp,hmax,hmin)
      IF (istat /= 0) GO TO 110
    END IF
    CALL mpisplit3d(var3d,nx,ny,nz,vbar)

    IF(myproc == 0) THEN
      CALL hdfrd3d(sd_id,"wbar",nxlg,nylg,nzlg,var3d,istat,itmp,hmax,hmin)
      IF (istat /= 0) GO TO 110
    END IF
    CALL mpisplit3d(var3d,nx,ny,nz,wbar)

    IF(myproc == 0) THEN
      CALL hdfrd3d(sd_id,"ptbar",nxlg,nylg,nzlg,var3d,istat,itmp,hmax,hmin)
      IF (istat /= 0) GO TO 110
    END IF
    CALL mpisplit3d(var3d,nx,ny,nz,ptbar)

    IF(myproc == 0) THEN
      CALL hdfrd3d(sd_id,"pbar",nxlg,nylg,nzlg,var3d,istat,itmp,hmax,hmin)
      IF (istat /= 0) GO TO 110
    END IF
    CALL mpisplit3d(var3d,nx,ny,nz,pbar)

    IF( mstin == 1) THEN
      IF(myproc == 0) THEN
        CALL hdfrd3d(sd_id,"qvbar",nxlg,nylg,nzlg,var3d,istat,itmp,hmax,hmin)
        IF (istat /= 0) GO TO 110
      END IF
      CALL mpisplit3d(var3d,nx,ny,nz,qvbar)
    END IF

    IF (landin == 1) THEN

      IF(myproc == 0) THEN
        CALL hdfrd3di(sd_id,"soiltyp",nxlg,nylg,nstypin,var3di,istat)
        IF (istat /= 0) GO TO 110
      END IF
      CALL mpisplit3di(var3di,nx,ny,nstypin,soiltyp)

      IF(myproc == 0) THEN
        CALL hdfrd3d(sd_id,"stypfrct",nxlg,nylg,nstypin,                &
                     var3d,istat,itmp,hmax,hmin)
        IF (istat /= 0) GO TO 110
      END IF
      CALL mpisplit3d(var3d,nx,ny,nstypin,stypfrct)

      CALL fix_stypfrct_nstyp(nx,ny,nstyp1,nstyp,stypfrct)

      IF(myproc == 0) THEN
        CALL hdfrd2di(sd_id,"vegtyp",nxlg,nylg,var3di,istat)
        IF (istat /= 0) GO TO 110
      END IF
      CALL mpisplit3di(var3di,nx,ny,1,vegtyp)

      IF(myproc == 0) THEN
        CALL hdfrd2d(sd_id,"lai",nxlg,nylg,var3d,istat,itmp)
        IF (istat /= 0) GO TO 110
      END IF
      CALL mpisplit3d(var3d,nx,ny,1,lai)

      IF(myproc == 0) THEN
        CALL hdfrd2d(sd_id,"roufns",nxlg,nylg,var3d,istat,itmp)
        IF (istat /= 0) GO TO 110
      END IF
      CALL mpisplit3d(var3d,nx,ny,1,roufns)

      IF(myproc == 0) THEN
        CALL hdfrd2d(sd_id,"veg",nxlg,nylg,var3d,istat,itmp)
        IF (istat /= 0) GO TO 110
      END IF
      CALL mpisplit3d(var3d,nx,ny,1,veg)

    END IF

  END IF

  IF( grdbas == 1 ) GO TO 930

  IF( varin == 1 ) THEN

    IF ( totin == 0 ) THEN

!-----------------------------------------------------------------------
!
!  Read in perturbations from history dump
!
!-----------------------------------------------------------------------

      IF(myproc == 0) THEN
        CALL hdfrd3d(sd_id,"uprt",nxlg,nylg,nzlg,var3d,istat,itmp,hmax,hmin)
        IF (istat /= 0) GO TO 110
      END IF
      CALL mpisplit3d(var3d,nx,ny,nz,uprt)

      IF(myproc == 0) THEN
        CALL hdfrd3d(sd_id,"vprt",nxlg,nylg,nzlg,var3d,istat,itmp,hmax,hmin)
        IF (istat /= 0) GO TO 110
      END IF
      CALL mpisplit3d(var3d,nx,ny,nz,vprt)

      IF(myproc == 0) THEN
        CALL hdfrd3d(sd_id,"wprt",nxlg,nylg,nzlg,var3d,istat,itmp,hmax,hmin)
        IF (istat /= 0) GO TO 110
      END IF
      CALL mpisplit3d(var3d,nx,ny,nz,wprt)

!-----------------------------------------------------------------------
!
!  Read in scalars
!
!-----------------------------------------------------------------------

      IF(myproc == 0) THEN
        CALL hdfrd3d(sd_id,"ptprt",nxlg,nylg,nzlg,var3d,istat,itmp,hmax,hmin)
        IF (istat /= 0) GO TO 110
      END IF
      CALL mpisplit3d(var3d,nx,ny,nz,ptprt)

      IF(myproc == 0) THEN
        CALL hdfrd3d(sd_id,"pprt",nxlg,nylg,nzlg,var3d,istat,itmp,hmax,hmin)
        IF (istat /= 0) GO TO 110
      END IF
      CALL mpisplit3d(var3d,nx,ny,nz,pprt)

    ELSE

!-----------------------------------------------------------------------
!
!  Read in total values of variables from history dump
!
!-----------------------------------------------------------------------

      IF(myproc == 0) THEN
        CALL hdfrd3d(sd_id,"u",nxlg,nylg,nzlg,var3d,istat,itmp,hmax,hmin)
        IF (istat /= 0) GO TO 110
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
        CALL hdfrd3d(sd_id,"v",nxlg,nylg,nzlg,var3d,istat,itmp,hmax,hmin)
        IF (istat /= 0) GO TO 110
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
        CALL hdfrd3d(sd_id,"w",nxlg,nylg,nzlg,var3d,istat,itmp,hmax,hmin)
        IF (istat /= 0) GO TO 110
      END IF
      CALL mpisplit3d(var3d,nx,ny,nz,wprt)

      IF(myproc == 0) THEN
        CALL hdfrd3d(sd_id,"pt",nxlg,nylg,nzlg,var3d,istat,itmp,hmax,hmin)
        IF (istat /= 0) GO TO 110
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
        CALL hdfrd3d(sd_id,"p",nxlg,nylg,nzlg,var3d,istat,itmp,hmax,hmin)
        IF (istat /= 0) GO TO 110
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

!-----------------------------------------------------------------------
!
!  Read in moisture variables
!
!-----------------------------------------------------------------------

  IF( mstin == 1 ) THEN

    IF ( totin == 0 ) THEN

      IF(myproc == 0) THEN
        CALL hdfrd3d(sd_id,"qvprt",nxlg,nylg,nzlg,var3d,istat,itmp,hmax,hmin)
        IF (istat /= 0) GO TO 110
      END IF
      CALL mpisplit3d(var3d,nx,ny,nz,qvprt)

    ELSE

      IF(myproc == 0) THEN
        CALL hdfrd3d(sd_id,"qv",nxlg,nylg,nzlg,var3d,istat,itmp,hmax,hmin)
        IF (istat /= 0) GO TO 110
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

    IF (intver < intver530) THEN

      IF(myproc == 0) THEN
        CALL hdfrd3d(sd_id,"qc",nxlg,nylg,nzlg,var3d,istat,itmp,hmax,hmin)
        IF (istat /= 0) GO TO 110
      END IF
      IF (P_QC > 0) CALL mpisplit3d(var3d,nx,ny,nz,qscalar(:,:,:,P_QC))

      IF(myproc == 0) THEN
        CALL hdfrd3d(sd_id,"qr",nxlg,nylg,nzlg,var3d,istat,itmp,hmax,hmin)
        IF (istat /= 0) GO TO 110
      END IF
      IF (P_QR > 0) CALL mpisplit3d(var3d,nx,ny,nz,qscalar(:,:,:,P_QR))

    ELSE

      DO nq = 1,nscalar
        IF( nqscalarin(nq) > 0) THEN
          IF(myproc == 0 ) THEN
            CALL hdfrd3d(sd_id,qnames(nq),nxlg,nylg,nzlg,var3d,istat,itmp,hmax,hmin)
            IF (istat /= 0) GO TO 110
          END IF
          CALL mpisplit3d(var3d,nx,ny,nz,qscalar(:,:,:,nq))
        END IF
      END DO

    END IF

    IF( rainin == 1 ) THEN
      IF(myproc == 0) THEN
        CALL hdfrd2d(sd_id,"raing",nxlg,nylg,var3d,istat,itmp)
        IF (istat /= 0) GO TO 110
      END IF
      CALL mpisplit3d(var3d,nx,ny,1,raing)

      IF(myproc == 0) THEN
        CALL hdfrd2d(sd_id,"rainc",nxlg,nylg,var3d,istat,itmp)
        IF (istat /= 0) GO TO 110
      END IF
      CALL mpisplit3d(var3d,nx,ny,1,rainc)
    END IF

    IF( prcin == 1 ) THEN
      IF(myproc == 0) THEN
        CALL hdfrd2d(sd_id,"prcrate1",nxlg,nylg,var3d(1,1,1),istat,itmp)
        IF (istat /= 0) GO TO 110

        CALL hdfrd2d(sd_id,"prcrate2",nxlg,nylg,var3d(1,1,2),istat,itmp)
        IF (istat /= 0) GO TO 110

        CALL hdfrd2d(sd_id,"prcrate3",nxlg,nylg,var3d(1,1,3),istat,itmp)
        IF (istat /= 0) GO TO 110

        CALL hdfrd2d(sd_id,"prcrate4",nxlg,nylg,var3d(1,1,4),istat,itmp)
        IF (istat /= 0) GO TO 110
      END IF
      CALL mpisplit3d(var3d,nx,ny,4,prcrate)
    END IF

    IF( icein == 1 .AND. intver < intver530) THEN

      IF(myproc == 0) THEN
        CALL hdfrd3d(sd_id,"qi",nxlg,nylg,nzlg,var3d,istat,itmp,hmax,hmin)
        IF (istat /= 0) GO TO 110
      END IF
      IF (P_QI > 0) CALL mpisplit3d(var3d,nx,ny,nz,qscalar(:,:,:,P_QI))

      IF(myproc == 0) THEN
        CALL hdfrd3d(sd_id,"qs",nxlg,nylg,nzlg,var3d,istat,itmp,hmax,hmin)
        IF (istat /= 0) GO TO 110
      END IF
      IF (P_QS > 0) CALL mpisplit3d(var3d,nx,ny,nz,qscalar(:,:,:,P_QS))

      IF(myproc == 0) THEN
        CALL hdfrd3d(sd_id,"qh",nxlg,nylg,nzlg,var3d,istat,itmp,hmax,hmin)
        IF (istat /= 0) GO TO 110
      END IF
      IF (P_QH > 0) CALL mpisplit3d(var3d,nx,ny,nz,qscalar(:,:,:,P_QH))

    END IF
  END IF

  IF( tkein == 1 ) THEN

    IF(myproc == 0) THEN
      CALL hdfrd3d(sd_id,"tke",nxlg,nylg,nzlg,var3d,istat,itmp,hmax,hmin)
      IF (istat /= 0) GO TO 110
    END IF
    CALL mpisplit3d(var3d,nx,ny,nz,tke)

  END IF

  IF( trbin == 1 ) THEN

    IF(myproc == 0) THEN
      CALL hdfrd3d(sd_id,"kmh",nxlg,nylg,nzlg,var3d,istat,itmp,hmax,hmin)
      IF (istat /= 0) GO TO 110
    END IF
    CALL mpisplit3d(var3d,nx,ny,nz,kmh)

    IF(myproc == 0) THEN
      CALL hdfrd3d(sd_id,"kmv",nxlg,nylg,nzlg,var3d,istat,itmp,hmax,hmin)
      IF (istat /= 0) GO TO 110
    END IF
    CALL mpisplit3d(var3d,nx,ny,nz,kmv)

  END IF

  IF( sfcin == 1 ) THEN

    IF (intver <= intver410) THEN

      ALLOCATE (tem2(nx,ny,0:nstyps),stat=istat)

      tem2=0.
      IF (myproc == 0) THEN
        CALL hdfrd3d(sd_id,"tsfc",nxlg,nylg,nstypin+1,var3d,istat,    &
                     itmp,hmax,hmin)
        IF (istat /= 0) GO TO 110
      END IF
      CALL mpisplit3d(var3d, nx, ny, nstypin+1, tem2)

      DO is=0,nstypin
        DO j=1,ny
          DO i=1,nx
            tsoil(i,j,1,is)=tem2(i,j,is)
          END DO
        END DO
      END DO

      tem2=0.
      IF (myproc == 0) THEN
        CALL hdfrd3d(sd_id,"tsoil",nxlg,nylg,nstypin+1,var3d,istat,   &
                     itmp,hmax,hmin)
        IF (istat /= 0) GO TO 110
      END IF
      CALL mpisplit3d(var3d, nx, ny, nstypin+1, tem2)

      DO is=0,nstypin
        DO j=1,ny
          DO i=1,nx
            tsoil(i,j,2,is)=tem2(i,j,is)
          END DO
        END DO
      END DO

      tem2=0.
      IF (myproc == 0) THEN
        CALL hdfrd3d(sd_id,"wetsfc",nxlg,nylg,nstypin+1,var3d,       &
                     istat,itmp,hmax,hmin)
        IF (istat /= 0) GO TO 110
      END IF
      CALL mpisplit3d(var3d, nx, ny, nstypin+1, tem2)

      DO is=0,nstypin
        DO j=1,ny
          DO i=1,nx
            qsoil(i,j,1,is)=tem2(i,j,is)
          END DO
        END DO
      END DO

      tem2 = 0.0
      IF (myproc == 0) THEN
        CALL hdfrd3d(sd_id,"wetdp",nxlg,nylg,nstypin+1,var3d,       &
                     istat,itmp,hmax,hmin)
        IF (istat /= 0) GO TO 110
      END IF
      CALL mpisplit3d(var3d, nx, ny, nstypin+1, tem2)

      DO is=0,nstypin
        DO j=1,ny
          DO i=1,nx
            qsoil(i,j,2,is)=tem2(i,j,is)
          END DO
        END DO
      END DO

      DEALLOCATE (tem2,stat=istat)

    ELSE IF (intver >= intver500) THEN

      IF(myproc == 0) THEN
        CALL hdfrd4d(sd_id,"tsoil",nxlg,nylg,nzsoillg,nstypin+1,var4d,istat, &
                       itmpsoil,hmaxsoil,hminsoil)
        IF (istat /= 0) GO TO 110
      END IF
      CALL mpisplit4d(var4d,nx,ny,nzsoil,nstypin+1,tsoil)

      IF(myproc == 0) THEN
      CALL hdfrd4d(sd_id,"qsoil",nxlg,nylg,nzsoillg,nstypin+1,var4d,istat, &
                       itmpsoil,hmaxsoil,hminsoil)
      IF (istat /= 0) GO TO 110
      END IF
      CALL mpisplit4d(var4d,nx,ny,nzsoil,nstypin+1,qsoil)

    END IF

    IF(myproc == 0) THEN
      CALL hdfrd3d(sd_id,"wetcanp",nxlg,nylg,nstypin+1,var3d,istat,  &
                   itmp,hmax,hmin)
      IF (istat /= 0) GO TO 110
    END IF
    CALL mpisplit3d(var3d,nx,ny,nstypin+1,wetcanp)

    CALL fix_soil_nstyp(nx,ny,nzsoil,nstyp1,nstyp,tsoil,qsoil,wetcanp)

    IF (snowin == 1) THEN
      IF(myproc == 0) THEN
        CALL hdfrd2d(sd_id,"snowdpth",nxlg,nylg,var3d,istat,itmp)
        IF (istat /= 0) GO TO 110
      END IF
      CALL mpisplit3d(var3d,nx,ny,1,snowdpth)
    END IF

  END IF

  IF( radin == 1 ) THEN

    IF(myproc == 0) THEN
      CALL hdfrd3d(sd_id,"radfrc",nxlg,nylg,nzlg,var3d,istat,         &
                   itmp,hmax,hmin)
      IF (istat /= 0) GO TO 110
    END IF
    CALL mpisplit3d(var3d,nx,ny,nz,radfrc)

    IF(myproc == 0) THEN
      CALL hdfrd2d(sd_id,"radsw",nxlg,nylg,var3d,istat,itmp)
      IF (istat /= 0) GO TO 110
    END IF
    CALL mpisplit3d(var3d,nx,ny,1,radsw)

    IF(myproc == 0) THEN
      CALL hdfrd2d(sd_id,"rnflx",nxlg,nylg,var3d,istat,itmp)
      IF (istat /= 0) GO TO 110
    END IF
    CALL mpisplit3d(var3d,nx,ny,1,rnflx)

    IF (intver <= intver410) THEN

      radswnet=0.
      radlwin=0.

    ELSE IF (intver >= intver500) THEN

      IF(myproc == 0) THEN
        CALL hdfrd2d(sd_id,"radswnet",nxlg,nylg,var3d,istat,itmp)
        IF (istat /= 0) GO TO 110
      END IF
      CALL mpisplit3d(var3d,nx,ny,1,radswnet)

      IF(myproc == 0) THEN
        CALL hdfrd2d(sd_id,"radlwin",nxlg,nylg,var3d,istat,itmp)
        IF (istat /= 0) GO TO 110
      END IF
      CALL mpisplit3d(var3d,nx,ny,1,radlwin)

    END IF

  END IF

  IF( flxin == 1 ) THEN

    IF(myproc == 0) THEN
      CALL hdfrd2d(sd_id,"usflx",nxlg,nylg,var3d,istat,itmp)
      IF (istat /= 0) GO TO 110
    END IF
    CALL mpisplit3d(var3d,nx,ny,1,usflx)

      IF(myproc == 0) THEN
    CALL hdfrd2d(sd_id,"vsflx",nxlg,nylg,var3d,istat,itmp)
    IF (istat /= 0) GO TO 110
      END IF
      CALL mpisplit3d(var3d,nx,ny,1,vsflx)

      IF(myproc == 0) THEN
    CALL hdfrd2d(sd_id,"ptsflx",nxlg,nylg,var3d,istat,itmp)
    IF (istat /= 0) GO TO 110
      END IF
      CALL mpisplit3d(var3d,nx,ny,1,ptsflx)

    IF(myproc == 0) THEN
      CALL hdfrd2d(sd_id,"qvsflx",nxlg,nylg,var3d,istat,itmp)
      IF (istat /= 0) GO TO 110
    END IF
    CALL mpisplit3d(var3d,nx,ny,1,qvsflx)

  END IF

!-----------------------------------------------------------------------
!
!  Friendly exit message
!
!-----------------------------------------------------------------------

  930   CONTINUE

  IF(myproc == 0) &
    WRITE(6,'(/a,F8.1,a/)')                                            &
      ' Data at time=', time/60,' (min) were successfully read.'

  ireturn = 0

  GO TO 130

!-----------------------------------------------------------------------
!
!  Error during read
!
!-----------------------------------------------------------------------

  110   CONTINUE
  IF(myproc == 0) &
    WRITE(6,'(/a/)') ' Error reading data in HDFREADSPLIT'
  ireturn=1

  GO TO 130

!-----------------------------------------------------------------------
!
!  End-of-file during read
!
!-----------------------------------------------------------------------

!  120   CONTINUE
!  IF(myproc == 0) &
!  WRITE(6,'(/a/)') ' End of file reached in HDFREADSPLIT'
!  ireturn=2

  130   CONTINUE

  IF(myproc == 0) THEN
    CALL hdfclose(sd_id,istat)
    IF (ireturn == 0) THEN
      IF (istat == 0) THEN
        WRITE(*,*) "HDFREADSPLIT: Successfully read ", trim(filename)
      ELSE
        WRITE(*,*) "HDFREADSPLIT: ERROR (status=", istat, ") closing ", trim(filename)
      END IF
    END IF
  END IF

  DEALLOCATE (itmp,stat=istat)
  DEALLOCATE (hmax,stat=istat)
  DEALLOCATE (hmin,stat=istat)
  DEALLOCATE (itmpsoil,stat=istat)
  DEALLOCATE (hmaxsoil,stat=istat)
  DEALLOCATE (hminsoil,stat=istat)

  DEALLOCATE (var1d, var3d, var3di, var4d)

  RETURN
END SUBROUTINE hdfreadsplit

!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE HDFDUMP                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################

SUBROUTINE hdfdump(nx,ny,nz,nzsoil,nstyps, filename, grdbas,            &
           u,v,w,ptprt,pprt,qv,qscalar,tke,kmh,kmv,                     &
           ubar,vbar,ptbar,pbar,rhobar,qvbar,                           &
           x,y,z,zp,zpsoil,                                             &
           soiltyp,stypfrct,vegtyp,lai,roufns,veg,                      &
           tsoil,qsoil,wetcanp,snowdpth,                                &
           raing,rainc,prcrate,                                         &
           radfrc,radsw,rnflx,radswnet,radlwin,                         &
           usflx,vsflx,ptsflx,qvsflx,                                   &
           tem1)

!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Produces a history data file "filename" in the NCSA HDF4 format by
!  calling HDF library subroutines.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Gene Bassett
!  2000/03/15
!
!  MODIFICATION HISTORY:
!
!  05/15/2002  (J. Brotzge)
!  Added to allow for multiple soil schemes
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!    nzsoil   Number of grid points in the vertical
!
!    filename File name of history dump data.
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
!    tem1     work array.
!    tem2     work array.
!
!-----------------------------------------------------------------------
!
!  The following parameters are passed into this subroutine through
!  a common block in globcst.inc, and they determine which
!  variables are output.
!
!  grdout =0 or 1. If grdout=0, grid variables are not dumped.
!  basout =0 or 1. If basout=0, base state variables are not dumped.
!  varout =0 or 1. If varout=0, perturbation variables are not dumped.
!  mstout =0 or 1. If mstout=0, water variables are not dumped.
!  rainout=0 or 1. If rainout=0, rain variables are not dumped.
!  prcout =0 or 1. If prcout=0, precipitation rates are not dumped.
!  iceout =0 or 1. If iceout=0, qi, qs and qh are not dumped.
!  tkeout =0 or 1. If tkeout=0, tke is not dumped.
!  trbout =0 or 1. If trbout=0, the eddy viscosity km is not dumped.
!  radout =0 or 1. If radout=0, the radiation arrays are not dumped.
!  flxout =0 or 1. If flxout=0, the surface fluxes are not dumped.
!
!-----------------------------------------------------------------------

  USE arps_precision

  IMPLICIT NONE

!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------

  INCLUDE 'globcst.inc'
  INCLUDE 'grid.inc'        ! Grid parameters
  INCLUDE 'phycst.inc'
  INCLUDE 'mp.inc'

!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------

  INTEGER :: nx,ny,nz          ! Number of grid points in 3 directions
  INTEGER :: nzsoil            ! Number of grid points in the soil

  CHARACTER (LEN=*) :: filename
  INTEGER :: grdbas            ! If this is a grid/base state dump

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
                               ! w-point of the soil.

  INTEGER :: nstyps             ! Number of soil type

  INTEGER :: soiltyp (nx,ny,nstyps)    ! Soil type
  REAL    :: stypfrct(nx,ny,nstyps)    ! Fraction of soil types
  INTEGER :: vegtyp (nx,ny)            ! Vegetation type
  REAL :: lai    (nx,ny)    ! Leaf Area Index
  REAL :: roufns (nx,ny)    ! Surface roughness
  REAL :: veg    (nx,ny)    ! Vegetation fraction

  REAL :: tsoil (nx,ny,nzsoil,0:nstyps)   ! Soil temperature (K)
  REAL :: qsoil (nx,ny,nzsoil,0:nstyps)   ! Soil moisture (m**3/m**3)
  REAL :: wetcanp(nx,ny,0:nstyps)         ! Canopy water amount
  REAL :: snowdpth(nx,ny)                 ! Snow depth (m)

  REAL :: raing(nx,ny)         ! Grid supersaturation rain
  REAL :: rainc(nx,ny)         ! Cumulus convective rain
  REAL :: prcrate(nx,ny,4)     ! precipitation rates (kg/(m**2*s))
                               ! prcrate(1,1,1) = total precip. rate
                               ! prcrate(1,1,2) = grid scale precip. rate
                               ! prcrate(1,1,3) = cumulative precip. rate
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

  REAL :: tem1  (nx,ny,nz)                     ! Temporary work array

  INTEGER(KIND=INT16), ALLOCATABLE :: itmp(:,:,:)             ! Temporary array
  REAL(SP),            ALLOCATABLE :: hmax(:), hmin(:)
  INTEGER(KIND=INT16), ALLOCATABLE :: itmpsoil(:,:,:,:)
  REAL(SP),            ALLOCATABLE :: hmaxsoil(:), hminsoil(:)! Temporary array
!
! 06/28/2002  Zuwen He
!
! Create a tem2 which will be used when dump data in previous version.
! It may be better to pass the array from the argument list
! in the future.
!
  REAL, ALLOCATABLE :: tem2(:,:,:) ! Temporary array

  REAL :: dx_out,dy_out

!-----------------------------------------------------------------------
!
!  Parameters describing this routine
!
!-----------------------------------------------------------------------
! 06/28/2002 Zuwen He
!
! fmtver??: to label each data a version.
! intver??: an integer to allow faster comparison than fmtver??,
!           which are strings.
!
! Verion 5.00: significant change in soil variables since version 4.10.
!
  CHARACTER (LEN=40) :: fmtver,fmtver410,fmtver500,fmtver530
  INTEGER  :: intver,intver410,intver500,intver530

  PARAMETER (fmtver410='004.10 HDF4 Coded Data',intver410=410)
  PARAMETER (fmtver500='005.00 HDF4 Coded Data',intver500=500)
  PARAMETER (fmtver530='005.30 HDF4 Coded Data',intver530=530)

  CHARACTER (LEN=10) :: tmunit
  PARAMETER (tmunit='seconds   ')

!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------

  INTEGER :: i,j,k,is
  INTEGER :: nq
  INTEGER :: nstypout
  INTEGER :: istat, sd_id

  INTEGER :: hdfcomprtmp, nztmp

  CHARACTER(LEN=256) :: cmntstr

!-----------------------------------------------------------------------
!
!  Functions
!
!-----------------------------------------------------------------------

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  IF ( (nproc_x_out > nproc_x .OR. nproc_y_out > nproc_y ) .AND. splitdmp > 0) THEN
    CALL hdfsplitdump(nx,ny,nz,nzsoil,nstyps, filename , grdbas,        &
             u,v,w,ptprt,pprt,qv,qscalar,tke,kmh,kmv,                   &
             ubar,vbar,ptbar,pbar,rhobar,qvbar,                         &
             x,y,z,zp,zpsoil,                                           &
             soiltyp,stypfrct,vegtyp,lai,roufns,veg,                    &
             tsoil,qsoil,wetcanp,snowdpth,                              &
             raing,rainc,prcrate,                                       &
             radfrc,radsw,rnflx,radswnet,radlwin,                       &
             usflx,vsflx,ptsflx,qvsflx,                                 &
             tem1)
    RETURN
  END IF

!-----------------------------------------------------------------------

  nstypout = max(1,nstyps)
  IF (hdfcompr > 3) THEN
!
!   If "nz" is a small number, the arrays may not be big enough for the
!   soil family of variables.
!
    nztmp = max(nz,nstypout+1)

    ALLOCATE (itmp(nx,ny,nztmp),stat=istat)
    CALL check_alloc_status(istat, "hdfdump:itmp")

    ALLOCATE (hmax(nztmp),stat=istat)
    CALL check_alloc_status(istat, "hdfdump:hmax")

    ALLOCATE (hmin(nztmp),stat=istat)
    CALL check_alloc_status(istat, "hdfdump:hmin")

    ALLOCATE (itmpsoil(nx,ny,nzsoil,0:nstyps),stat=istat)
    CALL check_alloc_status(istat, "hdfdump:itmpsoil")

    ALLOCATE (hmaxsoil(nzsoil),stat=istat)
    CALL check_alloc_status(istat, "hdfdump:hmaxsoil")

    ALLOCATE (hminsoil(nzsoil),stat=istat)
    CALL check_alloc_status(istat, "hdfdump:hminsoil")

  END IF

  IF ( grdbas == 1 ) THEN
    IF (myproc == 0) WRITE(6,'(1x,a,a/)')                               &
      'Writing HDF4 grid and base file to ',filename
  ELSE
    IF (myproc == 0) WRITE(6,'(1x,a,f13.3,a,a/)')                       &
      'Writing HDF4 data at time=', curtim,' into file ',filename
  ENDIF

!-----------------------------------------------------------------------
!
!  Create the HDF4 file.
!
!-----------------------------------------------------------------------

  CALL hdfopen(filename,2,sd_id)
  IF (sd_id < 0) THEN
    IF(myproc == 0) &
    WRITE (6,*) "HDFDUMP: ERROR creating HDF4 file: ",                  &
                trim(filename)
    GO TO 600
  END IF

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
    IF (myproc == 0) WRITE(6,'(/1x,a,i10,a/)')                        &
      ' Data format, intver=',intver, ', not found. The Job stopped.'
    CALL arpsstop('arpstop called from HDFDUMP. ',1)
  END IF

  CALL hdfwrtc(sd_id, 40, 'fmtver',  fmtver,  istat)
  CALL hdfwrtc(sd_id, 80, 'runname', runname, istat)
  CALL hdfwrti(sd_id,     'nocmnt',  nocmnt,  istat)

  IF( nocmnt > 0 ) THEN
    IF (nocmnt*80 > 256) THEN
      WRITE(*,'(1x,a,I0,a)') "ERROR: comment array is ",80*nocmnt,    &
                           ", but comment str is only 256 characters."
      CALL arpsstop('ERROR:Comment string is too short.',1)
    END IF
    DO i = 1,nocmnt
      cmntstr((i-1)*80+1:i*80) = cmnt(i)
    END DO
    CALL hdfwrtc(sd_id, 80*nocmnt, 'cmnt', cmntstr, istat)
  END IF

  CALL hdfwrtc(sd_id, 7, 'tmunit', 'seconds', istat)
  CALL hdfwrtr(sd_id,    'time',   curtim,    istat)

  CALL hdfwrti(sd_id, 'nx', nx, istat)
  CALL hdfwrti(sd_id, 'ny', ny, istat)
  CALL hdfwrti(sd_id, 'nz', nz, istat)

  IF (intver >= intver500)  THEN
    CALL hdfwrti(sd_id, 'nzsoil', nzsoil, istat)
  END IF

  IF( grdbas == 1 ) THEN

    CALL hdfwrti(sd_id, 'grdflg', 1, istat)
    CALL hdfwrti(sd_id, 'basflg', 1, istat)
    CALL hdfwrti(sd_id, 'varflg', 0, istat)
    CALL hdfwrti(sd_id, 'mstflg', 1, istat)
    CALL hdfwrti(sd_id, 'iceflg', 0, istat)
    CALL hdfwrti(sd_id, 'trbflg', 0, istat)
    CALL hdfwrti(sd_id, 'sfcflg', 0, istat)
    CALL hdfwrti(sd_id, 'rainflg', 0, istat)
    CALL hdfwrti(sd_id, 'landflg', landout, istat)
    CALL hdfwrti(sd_id, 'totflg', 1, istat)
    CALL hdfwrti(sd_id, 'tkeflg', 0, istat)

  ELSE

    CALL hdfwrti(sd_id, 'grdflg', grdout, istat)
    CALL hdfwrti(sd_id, 'basflg', basout, istat)
    CALL hdfwrti(sd_id, 'varflg', varout, istat)
    CALL hdfwrti(sd_id, 'mstflg', mstout, istat)
    CALL hdfwrti(sd_id, 'iceflg', iceout, istat)
    CALL hdfwrti(sd_id, 'trbflg', trbout, istat)
    CALL hdfwrti(sd_id, 'sfcflg', sfcout, istat)
    CALL hdfwrti(sd_id, 'rainflg', rainout, istat)
    CALL hdfwrti(sd_id, 'landflg', landout*basout, istat)
    CALL hdfwrti(sd_id, 'totflg', totout, istat)
    CALL hdfwrti(sd_id, 'tkeflg', tkeout, istat)

    CALL hdfwrti(sd_id, 'nscalar', nscalar, istat)
    CALL hdfwrti(sd_id, 'P_QC',    P_QC, istat)
    CALL hdfwrti(sd_id, 'P_QR',    P_QR, istat)
    CALL hdfwrti(sd_id, 'P_QI',    P_QI, istat)
    CALL hdfwrti(sd_id, 'P_QS',    P_QS, istat)
    CALL hdfwrti(sd_id, 'P_QG',    P_QG, istat)
    CALL hdfwrti(sd_id, 'P_QH',    P_QH, istat)

    CALL hdfwrti(sd_id, 'P_NC',    P_NC, istat)
    CALL hdfwrti(sd_id, 'P_NR',    P_NR, istat)
    CALL hdfwrti(sd_id, 'P_NI',    P_NI, istat)
    CALL hdfwrti(sd_id, 'P_NS',    P_NS, istat)
    CALL hdfwrti(sd_id, 'P_NG',    P_NG, istat)
    CALL hdfwrti(sd_id, 'P_NH',    P_NH, istat)

    CALL hdfwrti(sd_id, 'P_ZR',    P_ZR, istat)
    CALL hdfwrti(sd_id, 'P_ZI',    P_ZI, istat)
    CALL hdfwrti(sd_id, 'P_ZS',    P_ZS, istat)
    CALL hdfwrti(sd_id, 'P_ZG',    P_ZG, istat)
    CALL hdfwrti(sd_id, 'P_ZH',    P_ZH, istat)

    CALL hdfwrti(sd_id, 'P_CC',    P_CC, istat)

  END IF

  nstypout = max(1,nstyps)
  CALL hdfwrti(sd_id, 'nstyp', nstypout, istat)
  CALL hdfwrti(sd_id, 'prcflg', prcout, istat)
  CALL hdfwrti(sd_id, 'radflg', radout, istat)
  CALL hdfwrti(sd_id, 'flxflg', flxout, istat)
  CALL hdfwrti(sd_id, 'snowflg', snowout, istat)

  CALL hdfwrti(sd_id, 'day', day, istat)
  CALL hdfwrti(sd_id, 'year', year, istat)
  CALL hdfwrti(sd_id, 'month', month, istat)
  CALL hdfwrti(sd_id, 'hour', hour, istat)
  CALL hdfwrti(sd_id, 'minute', minute, istat)
  CALL hdfwrti(sd_id, 'second', second, istat)

  CALL hdfwrtr(sd_id, 'umove', umove, istat)
  CALL hdfwrtr(sd_id, 'vmove', vmove, istat)
  CALL hdfwrtr(sd_id, 'xgrdorg', xgrdorg, istat)
  CALL hdfwrtr(sd_id, 'ygrdorg', ygrdorg, istat)

  CALL hdfwrti(sd_id, 'mapproj', mapproj, istat)
  CALL hdfwrtr(sd_id, 'trulat1', trulat1, istat)
  CALL hdfwrtr(sd_id, 'trulat2', trulat2, istat)
  CALL hdfwrtr(sd_id, 'trulon', trulon, istat)
  CALL hdfwrtr(sd_id, 'sclfct', sclfct, istat)
  CALL hdfwrtr(sd_id, 'tstop', tstop, istat)
  CALL hdfwrtr(sd_id, 'thisdmp', thisdmp, istat)
  CALL hdfwrtr(sd_id, 'latitud', latitud, istat)
  CALL hdfwrtr(sd_id, 'ctrlat', ctrlat, istat)
  CALL hdfwrtr(sd_id, 'ctrlon', ctrlon, istat)

  CALL hdfwrtr(sd_id, 'ntcloud',  ntcloud,  istat)
  CALL hdfwrtr(sd_id, 'n0rain',  n0rain,  istat)
  CALL hdfwrtr(sd_id, 'n0snow',  n0snow,  istat)
  CALL hdfwrtr(sd_id, 'n0grpl',  n0grpl,  istat)
  CALL hdfwrtr(sd_id, 'n0hail',  n0hail,  istat)
  CALL hdfwrtr(sd_id, 'rhoice', rhoice, istat)
  CALL hdfwrtr(sd_id, 'rhosnow', rhosnow, istat)
  CALL hdfwrtr(sd_id, 'rhogrpl', rhogrpl, istat)
  CALL hdfwrtr(sd_id, 'rhohail', rhohail, istat)
  CALL hdfwrtr(sd_id, 'alpharain', alpharain, istat)
  CALL hdfwrtr(sd_id, 'alphaice', alphaice, istat)
  CALL hdfwrtr(sd_id, 'alphasnow', alphasnow, istat)
  CALL hdfwrtr(sd_id, 'alphagrpl', alphagrpl, istat)
  CALL hdfwrtr(sd_id, 'alphahail', alphahail, istat)

  dx_out = x(2) - x(1)
  dy_out = y(2) - y(1)

  CALL hdfwrtr(sd_id, 'dx', dx_out, istat)
  CALL hdfwrtr(sd_id, 'dy', dy_out, istat)

!-----------------------------------------------------------------------
!
!  If grdout=1 or grdbas=1, write out grid variables
!
!-----------------------------------------------------------------------

  IF(grdout == 1 .OR. grdbas == 1 ) THEN

    CALL hdfwrt1d(x,nx,sd_id,'x','x coordinate','m')
    CALL hdfwrt1d(y,ny,sd_id,'y','y coordinate','m')
    CALL hdfwrt1d(z,nz,sd_id,'z','z coordinate','m')
    CALL hdfwrt3d(zp,nx,ny,nz,sd_id,0,hdfcompr,                         &
                  'zp','Physical height coordinate','m',                &
                   itmp,hmax,hmin)

    IF (intver >= intver500) THEN
      CALL hdfwrt3d(zpsoil,nx,ny,nzsoil,sd_id,0,hdfcompr,               &
                  'zpsoil','Physical height coordinate (soil)','m',     &
                   itmpsoil,hmaxsoil,hminsoil)
    END IF

  END IF

!-----------------------------------------------------------------------
!
!  If basout=1, write out base state variables.
!
!-----------------------------------------------------------------------

  IF(basout == 1 .OR. grdbas == 1 ) THEN

    CALL edgfill(ubar,1,nx,1,nx, 1,ny,1,ny-1, 1,nz,1,nz-1)
    CALL hdfwrt3d(ubar,nx,ny,nz,sd_id,1,hdfcompr,                       &
                  'ubar','Base state u-velocity','m/s',                 &
                   itmp,hmax,hmin)

    CALL edgfill(vbar,1,nx,1,nx-1, 1,ny,1,ny, 1,nz,1,nz-1)
    CALL hdfwrt3d(vbar,nx,ny,nz,sd_id,2,hdfcompr,                       &
                  'vbar','Base state v-velocity','m/s',                 &
                   itmp,hmax,hmin)

    DO k=1,nz
      DO j=1,ny
        DO i=1,nx
          tem1(i,j,k) = 0.0
        END DO
      END DO
    END DO
    CALL hdfwrt3d(tem1,nx,ny,nz,sd_id,3,hdfcompr,                       &
                  'wbar','Base state w-velocity','m/s',                 &
                   itmp,hmax,hmin)

    CALL edgfill(ptbar,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)
    CALL hdfwrt3d(ptbar,nx,ny,nz,sd_id,0,hdfcompr,                      &
                  'ptbar','Base state potential temperature','K',       &
                   itmp,hmax,hmin)

    CALL edgfill(pbar,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)
    CALL hdfwrt3d(pbar,nx,ny,nz,sd_id,0,hdfcompr,                       &
                  'pbar','Base state pressure','Pascal',                &
                   itmp,hmax,hmin)

    IF(mstout == 1) THEN

      CALL edgfill(qvbar,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)
      CALL hdfwrt3d(qvbar,nx,ny,nz,sd_id,0,hdfcompr,                    &
          'qvbar','Base state water vapor specific humidity','kg/kg',   &
                   itmp,hmax,hmin)

    END IF

    IF(landout == 1) THEN

      CALL iedgfill(soiltyp(1,1,1),1,nx,1,nx-1, 1,ny,1,ny-1,            &
                      1,nstypout,1,nstypout)
      CALL hdfwrt3di(soiltyp,nx,ny,nstypout,sd_id,0,0,                  &
                'soiltyp','Soil type','index')
      CALL edgfill(stypfrct(1,1,1),1,nx,1,nx-1, 1,ny,1,ny-1,            &
                     1,nstypout,1,nstypout)
      CALL hdfwrt3d(stypfrct(1,1,1),nx,ny,nstypout,sd_id,0,hdfcompr,    &
          'stypfrct','Soil type fractional coverage','fraction',        &
                 itmp,hmax,hmin)

      CALL iedgfill(vegtyp ,1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)
      CALL hdfwrt2di(vegtyp,nx,ny,sd_id,0,0,                            &
                  'vegtyp','Vegetation type','index')

      CALL edgfill(lai    ,1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)
      CALL hdfwrt2d(lai,nx,ny,sd_id,0,hdfcompr,                         &
                  'lai','Leaf Area Index','index',itmp)

      CALL edgfill(roufns ,1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)
      CALL hdfwrt2d(roufns,nx,ny,sd_id,0,hdfcompr,                      &
                  'roufns','Surface roughness','0-1',itmp)

      CALL edgfill(veg    ,1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)
      CALL hdfwrt2d(veg,nx,ny,sd_id,0,hdfcompr,                         &
                  'veg','Vegetation fraction','fraction',itmp)

    END IF

  END IF

  IF ( grdbas == 1 ) GO TO 600

!-----------------------------------------------------------------------
!
!  If varout = 1, Write out uprt, vprt, wprt, ptprt, pprt.
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!  Write out u,v and w.
!
!-----------------------------------------------------------------------

  IF(varout == 1) THEN

    IF ( totout == 0 ) THEN

!-----------------------------------------------------------------------
!
!  Write out perturbations to history dump
!
!-----------------------------------------------------------------------

      DO k=1,nz-1
        DO j=1,ny-1
          DO i=1,nx
            tem1(i,j,k)=u(i,j,k)-ubar(i,j,k)
          END DO
        END DO
      END DO
      CALL edgfill(tem1,1,nx,1,nx, 1,ny,1,ny-1, 1,nz,1,nz-1)
      CALL hdfwrt3d(tem1,nx,ny,nz,sd_id,1,hdfcompr,                     &
                    'uprt','Perturbation u-velocity','m/s',             &
                     itmp,hmax,hmin)

      DO k=1,nz-1
        DO j=1,ny
          DO i=1,nx-1
            tem1(i,j,k)=v(i,j,k)-vbar(i,j,k)
          END DO
        END DO
      END DO
      CALL edgfill(tem1,1,nx,1,nx-1, 1,ny,1,ny, 1,nz,1,nz-1)
      CALL hdfwrt3d(tem1,nx,ny,nz,sd_id,2,hdfcompr,                     &
                    'vprt','Perturbation v-velocity','m/s',             &
                     itmp,hmax,hmin)

      CALL edgfill(w,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz)
      CALL hdfwrt3d(w,nx,ny,nz,sd_id,3,hdfcompr,                        &
                    'wprt','Perturbation w-velocity','m/s',             &
                     itmp,hmax,hmin)

!-----------------------------------------------------------------------
!
!  Write out scalars
!
!-----------------------------------------------------------------------

      CALL edgfill(ptprt,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)
      CALL hdfwrt3d(ptprt,nx,ny,nz,sd_id,0,hdfcompr,                    &
                    'ptprt','Perturbation potential temperature','K',   &
                     itmp,hmax,hmin)

      CALL edgfill(pprt,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)
      CALL hdfwrt3d(pprt,nx,ny,nz,sd_id,0,hdfcompr,                     &
                    'pprt','Perturbation pressure','Pascal',            &
                     itmp,hmax,hmin)

    ELSE

!-----------------------------------------------------------------------
!
!  Write out total values to history dump
!
!-----------------------------------------------------------------------

      CALL edgfill(u,1,nx,1,nx, 1,ny,1,ny-1, 1,nz,1,nz-1)
      CALL hdfwrt3d(u,nx,ny,nz,sd_id,1,hdfcompr,                        &
                    'u','u-velocity','m/s',                             &
                     itmp,hmax,hmin)

      CALL edgfill(v,1,nx,1,nx-1, 1,ny,1,ny, 1,nz,1,nz-1)
      CALL hdfwrt3d(v,nx,ny,nz,sd_id,2,hdfcompr,                        &
                    'v','v-velocity','m/s',                             &
                     itmp,hmax,hmin)

      CALL edgfill(w,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz)
      CALL hdfwrt3d(w,nx,ny,nz,sd_id,3,hdfcompr,                        &
                    'w','w-velocity','m/s',                             &
                     itmp,hmax,hmin)

!-----------------------------------------------------------------------
!
!  Write out scalars
!
!-----------------------------------------------------------------------

      DO k=1,nz-1
        DO j=1,ny-1
          DO i=1,nx-1
            tem1(i,j,k) = ptbar(i,j,k) + ptprt(i,j,k)
          END DO
        END DO
      END DO
      CALL edgfill(tem1,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)
      CALL hdfwrt3d(tem1,nx,ny,nz,sd_id,0,hdfcompr,                     &
                    'pt','Potential temperature','K',                   &
                     itmp,hmax,hmin)

      DO k=1,nz-1
        DO j=1,ny-1
          DO i=1,nx-1
            tem1(i,j,k) = pbar(i,j,k) + pprt(i,j,k)
          END DO
        END DO
      END DO
      CALL edgfill(tem1,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)
      CALL hdfwrt3d(tem1,nx,ny,nz,sd_id,0,hdfcompr,                     &
                    'p','Pressure','Pascal',                            &
                     itmp,hmax,hmin)

    END IF

  END IF     ! varout

!-----------------------------------------------------------------------
!
!  If mstout = 1, write out moisture scalars.
!
!-----------------------------------------------------------------------

  IF(mstout == 1) THEN

    IF( totout == 0 ) THEN

!-----------------------------------------------------------------------
!
!  Write out perturbation to history dump
!
!-----------------------------------------------------------------------

      DO k=1,nz-1
        DO j=1,ny-1
          DO i=1,nx-1
            tem1(i,j,k)=qv(i,j,k)-qvbar(i,j,k)
          END DO
        END DO
      END DO

      CALL edgfill(tem1,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)
      CALL hdfwrt3d(tem1,nx,ny,nz,sd_id,0,hdfcompr,                     &
          'qvprt','Pert. water vapor specific humidity','kg/kg',        &
                     itmp,hmax,hmin)

    ELSE

!-----------------------------------------------------------------------
!
!  Write out total values to history dump
!
!-----------------------------------------------------------------------

      CALL edgfill(qv,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)
      CALL hdfwrt3d(qv,nx,ny,nz,sd_id,0,hdfcompr,                       &
                    'qv','Water vapor specific humidity','kg/kg',       &
                     itmp,hmax,hmin)

    END IF

    DO nq = 1,nscalar
      CALL edgfill(qscalar(:,:,:,nq),1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)

      ! DTD: Turn off bit-packing for Z array if hdfcompr > 3 and mphyopt == 11

      hdfcomprtmp = hdfcompr
      IF(nq >= 13 .and. mphyopt == 11 .and. hdfcompr > 3) THEN
        hdfcomprtmp = hdfcompr - 3
      END IF
      CALL hdfwrt3d(qscalar(:,:,:,nq),nx,ny,nz,sd_id,0,hdfcomprtmp,     &
                    TRIM(qnames(nq)),TRIM(qdescp(nq)),'kg/kg',          &
                    itmp,hmax,hmin)
    END DO

    IF(rainout == 1) THEN

      CALL edgfill(raing,   1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)
      CALL hdfwrt2d(raing,nx,ny,sd_id,0,hdfcompr,                       &
                  'raing','Grid supersaturation rain','mm',itmp)

      CALL edgfill(rainc,1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)
      CALL hdfwrt2d(rainc,nx,ny,sd_id,0,hdfcompr,                       &
                  'rainc','Cumulus convective rain','mm',itmp)

    END IF   !rainout

    IF ( prcout == 1 ) THEN
      CALL edgfill(prcrate,1,nx,1,nx-1, 1,ny,1,ny-1, 1,4,1,4)
      CALL hdfwrt2d(prcrate(1,1,1),nx,ny,sd_id,0,hdfcompr,              &
          'prcrate1','Total precip. rate','kg/(m**2*s)',itmp)
      CALL hdfwrt2d(prcrate(1,1,2),nx,ny,sd_id,0,hdfcompr,              &
          'prcrate2','Grid scale precip. rate','kg/(m**2*s)',itmp)
      CALL hdfwrt2d(prcrate(1,1,3),nx,ny,sd_id,0,hdfcompr,              &
          'prcrate3','Cumulative precip. rate','kg/(m**2*s)',itmp)
      CALL hdfwrt2d(prcrate(1,1,4),nx,ny,sd_id,0,hdfcompr,              &
          'prcrate4','Microphysics precip. rate','kg/(m**2*s)',itmp)
    END IF

  END IF   !mstout

!-----------------------------------------------------------------------
!
!  If tkeout = 1, write out tke.
!
!-----------------------------------------------------------------------

  IF( tkeout == 1 ) THEN

    CALL edgfill(tke,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)
    CALL hdfwrt3d(tke,nx,ny,nz,sd_id,0,hdfcompr,                        &
                  'tke','Turbulent Kinetic Energy','(m/s)**2',          &
                   itmp,hmax,hmin)

  END IF

!-----------------------------------------------------------------------
!
!  If trbout = 1, write out the turbulence parameter, km.
!
!-----------------------------------------------------------------------

  IF( trbout == 1 ) THEN

    CALL edgfill(kmh,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)
    CALL hdfwrt3d(kmh,nx,ny,nz,sd_id,0,hdfcompr,                        &
        'kmh','Hori. turb. mixing coef. for momentum','m**2/s',         &
                   itmp,hmax,hmin)

    CALL edgfill(kmv,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)
    CALL hdfwrt3d(kmv,nx,ny,nz,sd_id,0,hdfcompr,                        &
        'kmv','Vert. turb. mixing coef. for momentum','m**2/s',         &
                   itmp,hmax,hmin)

  END IF   ! trbout


!-----------------------------------------------------------------------
!
!  If sfcout = 1, write out the surface variables,
!  tsoil, qsoil, and wetcanp.
!
!-----------------------------------------------------------------------

  IF( sfcout == 1) THEN

    DO is=0,nstypout

      CALL edgfill(tsoil(1,1,1,is),  1,nx,1,nx-1, 1,ny,1,ny-1,          &
                   1,nzsoil,1,nzsoil)
      CALL edgfill(qsoil(1,1,1,is), 1,nx,1,nx-1, 1,ny,1,ny-1,           &
                   1,nzsoil,1,nzsoil)
      CALL edgfill(wetcanp(1,1,is),1,nx,1,nx-1, 1,ny,1,ny-1,            &
                   1,1,1,1)

    END DO

    IF (intver == intver410) THEN

!
! 06/28/2002 Zuwen He
!
! In version 410, tsoil is the average of tsoil from 2 to nzsoil in later
! version, and wetdp is similar.
!
      ALLOCATE (tem2(nx,ny,0:nstyps),stat=istat)

      tem2=0.
      DO is=0,nstypout
        DO j=1,ny
          DO i=1,nx
            tem2(i,j,is)=tsoil(i,j,1,is)
          END DO
        END DO
      END DO

      CALL hdfwrt3d(tem2,nx,ny,nstypout+1,sd_id,0,hdfcompr,             &
                    'tsfc','Surface ground temperature','K',            &
                    itmp,hmax,hmin)

      tem2=0.
      DO is=0,nstypout

        DO k=2,nzsoil
          DO j=1,ny
            DO i=1,nx
              tem2(i,j,is)=tem2(i,j,is)+tsoil(i,j,k,is)
            END DO
          END DO
        END DO

        DO j=1,ny
          DO i=1,nx
            tem1(i,j,is)=tem1(i,j,is)/float(nzsoil-1)
          END DO
        END DO

      END DO

      CALL hdfwrt3d(tem1,nx,ny,nstypout+1,sd_id,0,hdfcompr,             &
                 'tsoil','Deep soil temperature','K',                   &
                 itmp,hmax,hmin)


      tem2=0.
      DO is=0,nstypout
        DO j=1,ny
          DO i=1,nx
            tem2(i,j,is)=qsoil(i,j,1,is)
          END DO
        END DO
      END DO

      CALL hdfwrt3d(tem2,nx,ny,nstypout+1,sd_id,0,hdfcompr,             &
                    'wetsfc','Surface soil moisture','fraction',        &
                    itmp,hmax,hmin)

      tem2=0.
      DO is=0,nstypout

        DO k=2,nzsoil
          DO j=1,ny
            DO i=1,nx
              tem2(i,j,is)=tem2(i,j,is)+qsoil(i,j,k,is)
            END DO
          END DO
        END DO

        DO j=1,ny
          DO i=1,nx
            tem1(i,j,is)=tem1(i,j,is)/float(nzsoil-1)
          END DO
        END DO

      END DO

      CALL hdfwrt3d(tem1,nx,ny,nstypout+1,sd_id,0,hdfcompr,             &
                    'wetdp','Deep soil moisture','fraction',            &
                    itmp,hmax,hmin)

      DEALLOCATE (tem2,stat=istat)

    ELSE IF (intver >= intver500) THEN

      CALL hdfwrt4d(tsoil,nx,ny,nzsoil,nstypout+1,sd_id,0,hdfcompr,     &
                  'tsoil','Soil temperature','K',                       &
                   itmpsoil,hmaxsoil,hminsoil)
      CALL hdfwrt4d(qsoil,nx,ny,nzsoil,nstypout+1,sd_id,0,hdfcompr,     &
                  'qsoil','Soil moisture','fraction',                   &
                   itmpsoil,hmaxsoil,hminsoil)
    END IF

    CALL hdfwrt3d(wetcanp,nx,ny,nstypout+1,sd_id,0,hdfcompr,            &
                'wetcanp','Canopy water amount','fraction',             &
                 itmp,hmax,hmin)

    IF (snowout == 1) THEN

      CALL edgfill(snowdpth,1,nx,1,nx-1, 1,ny,1,ny-1,1,1,1,1)
      CALL hdfwrt2d(snowdpth,nx,ny,sd_id,0,hdfcompr,                    &
                  'snowdpth','Snow depth','m',itmp)
    END IF

  END IF

!-----------------------------------------------------------------------
!
!  If radout = 1, write out the radiation arrays
!
!-----------------------------------------------------------------------

  IF( radout == 1 ) THEN

    CALL edgfill(radfrc,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)
    CALL hdfwrt3d(radfrc,nx,ny,nz,sd_id,0,hdfcompr,                     &
                  'radfrc','Radiation forcing','K/s',                   &
                   itmp,hmax,hmin)

    CALL edgfill(radsw,1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)
    CALL hdfwrt2d(radsw,nx,ny,sd_id,0,hdfcompr,                         &
        'radsw','Solar radiation reaching the surface','W/m**2',itmp)

    CALL edgfill(rnflx,1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)
    CALL hdfwrt2d(rnflx,nx,ny,sd_id,0,hdfcompr,                         &
        'rnflx','Net radiation flux absorbed by surface','W/m**2',      &
        itmp)

    IF (intver >= intver500) THEN

      CALL edgfill(radswnet,1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)
      CALL hdfwrt2d(radswnet,nx,ny,sd_id,0,hdfcompr,                    &
          'radswnet','Net solar radiation','W/m**2',itmp)

      CALL edgfill(radlwin,1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)
      CALL hdfwrt2d(radlwin,nx,ny,sd_id,0,hdfcompr,                     &
          'radlwin','Incoming longwave radiation','W/m**2',             &
          itmp)

    END IF

  END IF   ! radout

!-----------------------------------------------------------------------
!
!  If flxout = 1, write out the surface fluxes
!
!-----------------------------------------------------------------------

  IF( flxout == 1 ) THEN

    CALL edgfill(usflx,1,nx,1,nx, 1,ny,1,ny-1, 1,1,1,1)
    CALL hdfwrt2d(usflx,nx,ny,sd_id,0,hdfcompr,                         &
        'usflx','Surface flux of u-momentum','kg/(m*s**2)',itmp)

    CALL edgfill(vsflx,1,nx,1,nx-1, 1,ny,1,ny, 1,1,1,1)
    CALL hdfwrt2d(vsflx,nx,ny,sd_id,0,hdfcompr,                         &
        'vsflx','Surface flux of v-momentum','kg/(m*s**2)',itmp)

    CALL edgfill(ptsflx,1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)
    CALL hdfwrt2d(ptsflx,nx,ny,sd_id,0,hdfcompr,                        &
        'ptsflx','Surface heat flux','K*kg/(m**2*s)',itmp)

    CALL edgfill(qvsflx,1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)
    CALL hdfwrt2d(qvsflx,nx,ny,sd_id,0,hdfcompr,                        &
        'qvsflx','Surface moisture flux','kg/(m**2*s)',itmp)

  END IF   ! flxout

  600   CONTINUE

  CALL hdfclose(sd_id,istat)
  IF (istat /= 0) THEN
    IF(myproc == 0) &
    WRITE (6,*) "HDFDUMP: ERROR on closing file ",trim(filename),       &
                " (status",istat,")"
  END IF

  IF (hdfcompr > 3) THEN
    DEALLOCATE (itmp,stat=istat)
    DEALLOCATE (hmax,stat=istat)
    DEALLOCATE (hmin,stat=istat)
    DEALLOCATE (itmpsoil,stat=istat)
    DEALLOCATE (hmaxsoil,stat=istat)
    DEALLOCATE (hminsoil,stat=istat)
  END IF

  RETURN
END SUBROUTINE hdfdump

!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE HDFJOINDUMP                ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################

SUBROUTINE hdfjoindump(nx,ny,nz,nzsoil,nstyps, filename , grdbas,       &
           u,v,w,ptprt,pprt,qv,qscalar,tke,kmh,kmv,                     &
           ubar,vbar,ptbar,pbar,rhobar,qvbar,                           &
           x,y,z,zp,zpsoil,                                             &
           soiltyp,stypfrct,vegtyp,lai,roufns,veg,                      &
           tsoil,qsoil,wetcanp,snowdpth,                                &
           raing,rainc,prcrate,                                         &
           radfrc,radsw,rnflx,radswnet,radlwin,                         &
           usflx,vsflx,ptsflx,qvsflx,                                   &
           tem1)

!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Produces a joined history data file "filename" in the NCSA HDF4 format
!  for MP mode by calling HDF library subroutines.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yunheng Wang
!  2002/08/15
!  Based on subroutine hdfdump
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
!    nzsoil   Number of grid points in the vertical
!
!    filename File name of history dump data.
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
!    tem1     work array.
!    tem2     work array.
!
!    out1d      work array.
!    out3d      work array.
!    out3di     work array.
!    outtsoil   work array.
!    outqsoil   work array.
!
!-----------------------------------------------------------------------
!
!  The following parameters are passed into this subroutine through
!  a common block in globcst.inc, and they determine which
!  variables are output.
!
!  grdout =0 or 1. If grdout=0, grid variables are not dumped.
!  basout =0 or 1. If basout=0, base state variables are not dumped.
!  varout =0 or 1. If varout=0, perturbation variables are not dumped.
!  mstout =0 or 1. If mstout=0, water variables are not dumped.
!  rainout=0 or 1. If rainout=0, rain variables are not dumped.
!  prcout =0 or 1. If prcout=0, precipitation rates are not dumped.
!  iceout =0 or 1. If iceout=0, qi, qs and qh are not dumped.
!  tkeout =0 or 1. If tkeout=0, tke is not dumped.
!  trbout =0 or 1. If trbout=0, the eddy viscosity km is not dumped.
!  radout =0 or 1. If radout=0, the radiation arrays are not dumped.
!  flxout =0 or 1. If flxout=0, the surface fluxes are not dumped.
!
!-----------------------------------------------------------------------

  USE arps_precision

  IMPLICIT NONE

!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------

  INCLUDE 'globcst.inc'
  INCLUDE 'grid.inc'        ! Grid parameters
  INCLUDE 'phycst.inc'
  INCLUDE 'mp.inc'

!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------

  INTEGER :: nx,ny,nz          ! Number of grid points in 3 directions
  INTEGER :: nzsoil            ! Number of grid points in the soil

  CHARACTER (LEN=*) :: filename
  INTEGER :: grdbas            ! If this is a grid/base state dump

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
                               ! w-point of the soil.

  INTEGER :: nstyps             ! Number of soil type

  INTEGER :: soiltyp (nx,ny,nstyps)    ! Soil type
  REAL :: stypfrct(nx,ny,nstyps)       ! Fraction of soil types
  INTEGER :: vegtyp (nx,ny)            ! Vegetation type
  REAL :: lai    (nx,ny)    ! Leaf Area Index
  REAL :: roufns (nx,ny)    ! Surface roughness
  REAL :: veg    (nx,ny)    ! Vegetation fraction

  REAL :: tsoil (nx,ny,nzsoil,0:nstyps)   ! Soil temperature (K)
  REAL :: qsoil (nx,ny,nzsoil,0:nstyps)   ! Soil moisture (m**3/m**3)
  REAL :: wetcanp(nx,ny,0:nstyps)         ! Canopy water amount
  REAL :: snowdpth(nx,ny)                 ! Snow depth (m)

  REAL :: raing(nx,ny)         ! Grid supersaturation rain
  REAL :: rainc(nx,ny)         ! Cumulus convective rain
  REAL :: prcrate(nx,ny,4)     ! precipitation rates (kg/(m**2*s))
                               ! prcrate(1,1,1) = total precip. rate
                               ! prcrate(1,1,2) = grid scale precip. rate
                               ! prcrate(1,1,3) = cumulative precip. rate
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

  REAL :: tem1  (nx,ny,nz)                     ! Temporary work array
  INTEGER(KIND=INT16), ALLOCATABLE :: itmp(:,:,:)       ! Temporary array
  REAL(SP),            ALLOCATABLE :: hmax(:), hmin(:)  ! Temporary array
  INTEGER(KIND=INT16), ALLOCATABLE :: itmpsoil(:,:,:,:) ! Temporary array
  REAL(SP),            ALLOCATABLE :: hmaxsoil(:), hminsoil(:)! Temporary array

  REAL :: dx_out,dy_out

  REAL,    ALLOCATABLE :: out1d(:)
  REAL,    ALLOCATABLE :: out3d(:,:,:)
  REAL,    ALLOCATABLE :: outtsoil(:,:,:,:)
  REAL,    ALLOCATABLE :: outqsoil(:,:,:,:)
  INTEGER, ALLOCATABLE :: out3di(:,:,:)

!-----------------------------------------------------------------------
!
!  Parameters describing this routine
!
!-----------------------------------------------------------------------
!
! fmtver??: to label each data a version.
! intver??: an integer to allow faster comparison than fmtver??,
!           which are strings.
!
! Verion 5.00: significant change in soil variables since version 4.10.
!
  CHARACTER (LEN=40) :: fmtver,fmtver410,fmtver500,fmtver530
  INTEGER  :: intver,intver410,intver500,intver530

  PARAMETER (fmtver410='004.10 HDF4 Coded Data',intver410=410)
  PARAMETER (fmtver500='005.00 HDF4 Coded Data',intver500=500)
  PARAMETER (fmtver530='005.30 HDF4 Coded Data',intver530=530)

  CHARACTER(LEN=10), PARAMETER :: tmunit = 'seconds   '

!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------

  INTEGER :: i,j,k,is
  INTEGER :: nq
  INTEGER :: nstypout
  INTEGER :: istat, sd_id
  INTEGER :: nztmp

  INTEGER :: nxlg, nylg
  INTEGER :: n3d
  INTEGER :: hdfcomprtmp

  CHARACTER(LEN=256) :: cmntstr

!-----------------------------------------------------------------------
!
!  Functions
!
!-----------------------------------------------------------------------

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  nxlg = nproc_x*(nx-3)+3
  nylg = nproc_y*(ny-3)+3

  nstypout = MAX(1,nstyps)
  n3d      = MAX(nz, nzsoil, nstypout+1, 4)   ! 3rd dimenson for out3d

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
      WRITE(6,'(/1x,a,i10,a/)')                        &
        ' Data format, intver=',intver,                &
        ', not found. The Job stopped.'
      CALL arpsstop('arpstop called from HDFDUMP. ',1)
  END IF

  IF (hdfcompr > 3) THEN

!
!   If "nz" is a small number, the arrays may not be big enough for the
!   soil family of variables.
!
    nztmp = max(nz,nstypout+1)

    ALLOCATE (itmp(nxlg,nylg,nztmp),stat=istat)
    CALL check_alloc_status(istat,'HDFJOINDUMP:itmp')

    ALLOCATE (hmax(nztmp),stat=istat)
    CALL check_alloc_status(istat,'HDFJOINDUMP:hmax')

    ALLOCATE (hmin(nztmp),stat=istat)
    CALL check_alloc_status(istat,'HDFJOINDUMP:hmin')

    ALLOCATE (itmpsoil(nxlg,nylg,nzsoil,0:nstyps),stat=istat)
    CALL check_alloc_status(istat,'HDFJOINDUMP:itmpsoil')

    ALLOCATE (hmaxsoil(nzsoil),stat=istat)
    CALL check_alloc_status(istat,'HDFJOINDUMP:hmaxsoil')

    ALLOCATE (hminsoil(nzsoil),stat=istat)
    CALL check_alloc_status(istat,'HDFJOINDUMP:hminsoil')

  END IF

  IF(myproc == 0) THEN

    WRITE(6,'(1x,a,f13.3,a,a/)')                                        &
       'Writing HDF4 data at time=', curtim,' into file ',filename

!-----------------------------------------------------------------------
!
!  Create the HDF4 file.
!
!-----------------------------------------------------------------------

    CALL hdfopen(filename,2,sd_id)
    IF (sd_id < 0) THEN
      WRITE (6,*) "HDFJOINDUMP: ERROR creating HDF4 file: ",            &
                trim(filename)
      GO TO 600
    END IF

    CALL hdfwrtc(sd_id, 40, 'fmtver', fmtver,  istat)
    CALL hdfwrtc(sd_id, 80, 'runname', runname, istat)
    CALL hdfwrti(sd_id,     'nocmnt', nocmnt,  istat)

    IF( nocmnt > 0 ) THEN
      IF (nocmnt*80 > 256) THEN
        WRITE(*,'(1x,a,I0,a)') "ERROR: comment array is ",80*nocmnt,    &
                             ", but comment str is only 256 characters."
        CALL arpsstop('ERROR:Comment string is too short.',1)
      END IF
      DO i = 1,nocmnt
        cmntstr((i-1)*80+1:i*80) = cmnt(i)
      END DO
      CALL hdfwrtc(sd_id, 80*nocmnt, 'cmnt', cmntstr, istat)
    END IF

    CALL hdfwrtc(sd_id, 7, 'tmunit', 'seconds', istat)
    CALL hdfwrtr(sd_id, 'time', curtim, istat)

    CALL hdfwrti(sd_id, 'nx', nxlg, istat)
    CALL hdfwrti(sd_id, 'ny', nylg, istat)
    CALL hdfwrti(sd_id, 'nz', nz, istat)

    IF (intver >= intver500)  THEN
      CALL hdfwrti(sd_id, 'nzsoil', nzsoil, istat)
    END IF

    IF( grdbas == 1 ) THEN

      CALL hdfwrti(sd_id, 'grdflg', 1, istat)
      CALL hdfwrti(sd_id, 'basflg', 1, istat)
      CALL hdfwrti(sd_id, 'varflg', 0, istat)
      CALL hdfwrti(sd_id, 'mstflg', 1, istat)
      CALL hdfwrti(sd_id, 'iceflg', 0, istat)
      CALL hdfwrti(sd_id, 'trbflg', 0, istat)
      CALL hdfwrti(sd_id, 'sfcflg', 0, istat)
      CALL hdfwrti(sd_id, 'rainflg', 0, istat)
      CALL hdfwrti(sd_id, 'landflg', landout, istat)
      CALL hdfwrti(sd_id, 'totflg', 1, istat)
      CALL hdfwrti(sd_id, 'tkeflg', 0, istat)

    ELSE

      CALL hdfwrti(sd_id, 'grdflg', grdout, istat)
      CALL hdfwrti(sd_id, 'basflg', basout, istat)
      CALL hdfwrti(sd_id, 'varflg', varout, istat)
      CALL hdfwrti(sd_id, 'mstflg', mstout, istat)
      CALL hdfwrti(sd_id, 'iceflg', iceout, istat)
      CALL hdfwrti(sd_id, 'trbflg', trbout, istat)
      CALL hdfwrti(sd_id, 'sfcflg', sfcout, istat)
      CALL hdfwrti(sd_id, 'rainflg', rainout, istat)
      CALL hdfwrti(sd_id, 'landflg', landout*basout, istat)
      CALL hdfwrti(sd_id, 'totflg', totout, istat)
      CALL hdfwrti(sd_id, 'tkeflg', tkeout, istat)

      CALL hdfwrti(sd_id, 'nscalar', nscalar, istat)
      CALL hdfwrti(sd_id, 'P_QC',    P_QC, istat)
      CALL hdfwrti(sd_id, 'P_QR',    P_QR, istat)
      CALL hdfwrti(sd_id, 'P_QI',    P_QI, istat)
      CALL hdfwrti(sd_id, 'P_QS',    P_QS, istat)
      CALL hdfwrti(sd_id, 'P_QG',    P_QG, istat)
      CALL hdfwrti(sd_id, 'P_QH',    P_QH, istat)

      CALL hdfwrti(sd_id, 'P_NC',    P_NC, istat)
      CALL hdfwrti(sd_id, 'P_NR',    P_NR, istat)
      CALL hdfwrti(sd_id, 'P_NI',    P_NI, istat)
      CALL hdfwrti(sd_id, 'P_NS',    P_NS, istat)
      CALL hdfwrti(sd_id, 'P_NG',    P_NG, istat)
      CALL hdfwrti(sd_id, 'P_NH',    P_NH, istat)

      CALL hdfwrti(sd_id, 'P_ZR',    P_ZR, istat)
      CALL hdfwrti(sd_id, 'P_ZI',    P_ZI, istat)
      CALL hdfwrti(sd_id, 'P_ZS',    P_ZS, istat)
      CALL hdfwrti(sd_id, 'P_ZG',    P_ZG, istat)
      CALL hdfwrti(sd_id, 'P_ZH',    P_ZH, istat)

      CALL hdfwrti(sd_id, 'P_CC',    P_CC, istat)

    END IF

    CALL hdfwrti(sd_id, 'nstyp', nstypout, istat)
    CALL hdfwrti(sd_id, 'prcflg', prcout, istat)
    CALL hdfwrti(sd_id, 'radflg', radout, istat)
    CALL hdfwrti(sd_id, 'flxflg', flxout, istat)
    CALL hdfwrti(sd_id, 'snowflg', snowout, istat)

    CALL hdfwrti(sd_id, 'day',   day,   istat)
    CALL hdfwrti(sd_id, 'year',  year,  istat)
    CALL hdfwrti(sd_id, 'month', month, istat)
    CALL hdfwrti(sd_id, 'hour',  hour,  istat)
    CALL hdfwrti(sd_id, 'minute',minute,istat)
    CALL hdfwrti(sd_id, 'second',second,istat)

    CALL hdfwrtr(sd_id, 'umove', umove, istat)
    CALL hdfwrtr(sd_id, 'vmove', vmove, istat)
    CALL hdfwrtr(sd_id, 'xgrdorg', xgrdorg, istat)
    CALL hdfwrtr(sd_id, 'ygrdorg', ygrdorg, istat)

    CALL hdfwrti(sd_id, 'mapproj', mapproj, istat)
    CALL hdfwrtr(sd_id, 'trulat1', trulat1, istat)
    CALL hdfwrtr(sd_id, 'trulat2', trulat2, istat)
    CALL hdfwrtr(sd_id, 'trulon',  trulon,  istat)
    CALL hdfwrtr(sd_id, 'sclfct',  sclfct,  istat)
    CALL hdfwrtr(sd_id, 'tstop',   tstop,   istat)
    CALL hdfwrtr(sd_id, 'thisdmp', thisdmp, istat)
    CALL hdfwrtr(sd_id, 'latitud', latitud, istat)
    CALL hdfwrtr(sd_id, 'ctrlat',  ctrlat,  istat)
    CALL hdfwrtr(sd_id, 'ctrlon',  ctrlon,  istat)

    CALL hdfwrtr(sd_id, 'ntcloud',  ntcloud,  istat)
    CALL hdfwrtr(sd_id, 'n0rain',  n0rain,  istat)
    CALL hdfwrtr(sd_id, 'n0snow',  n0snow,  istat)
    CALL hdfwrtr(sd_id, 'n0grpl',  n0grpl,  istat)
    CALL hdfwrtr(sd_id, 'n0hail',  n0hail,  istat)
    CALL hdfwrtr(sd_id, 'rhoice', rhoice, istat)
    CALL hdfwrtr(sd_id, 'rhosnow', rhosnow, istat)
    CALL hdfwrtr(sd_id, 'rhogrpl', rhogrpl, istat)
    CALL hdfwrtr(sd_id, 'rhohail', rhohail, istat)
    CALL hdfwrtr(sd_id, 'alpharain', alpharain, istat)
    CALL hdfwrtr(sd_id, 'alphaice', alphaice, istat)
    CALL hdfwrtr(sd_id, 'alphasnow', alphasnow, istat)
    CALL hdfwrtr(sd_id, 'alphagrpl', alphagrpl, istat)
    CALL hdfwrtr(sd_id, 'alphahail', alphahail, istat)

    dx_out = x(2) - x(1)
    dy_out = y(2) - y(1)
    CALL hdfwrtr(sd_id, 'dx', dx_out, istat)
    CALL hdfwrtr(sd_id, 'dy', dy_out, istat)

  END IF  ! myproc == 0

  ALLOCATE (out1d( MAX(nxlg,nylg) ),stat=istat)
  CALL check_alloc_status(istat,'HDFJOINDUMP:out1d')

  ALLOCATE (out3d( nxlg,nylg, n3d ),stat=istat)
  CALL check_alloc_status(istat,'HDFJOINDUMP:out3d')

  ALLOCATE (out3di( nxlg,nylg, nstypout ),stat=istat)
  CALL check_alloc_status(istat,'HDFJOINDUMP:out3di')

  ALLOCATE (outtsoil( nxlg,nylg, nzsoil, 0:nstypout ),stat=istat)
  CALL check_alloc_status(istat,'HDFJOINDUMP:outtsoil')

  ALLOCATE (outqsoil( nxlg,nylg, nzsoil, 0:nstypout ),stat=istat)
  CALL check_alloc_status(istat,'HDFJOINDUMP:outqsoil')

!-----------------------------------------------------------------------
!
!  If grdout=1 or grdbas=1, write out grid variables
!
!-----------------------------------------------------------------------

  IF(grdout == 1 .OR. grdbas == 1 ) THEN

    CALL mpimerge1dx(x,nx,out1d)
    IF (myproc == 0) &
      CALL hdfwrt1d(out1d,nxlg,sd_id,'x','x coordinate','m')

    CALL mpimerge1dy(y,ny,out1d)
    IF (myproc == 0) &
      CALL hdfwrt1d(out1d,nylg,sd_id,'y','y coordinate','m')

    IF (myproc == 0) &
      CALL hdfwrt1d(z,nz,sd_id,'z','z coordinate','m')

    CALL mpimerge3d(zp,nx,ny,nz,out3d)
    IF (myproc == 0) &
      CALL hdfwrt3d(out3d,nxlg,nylg,nz,sd_id,0,hdfcompr,                &
                 'zp','Physical height coordinate','m',itmp,hmax,hmin)

    IF (intver >= intver500) THEN
      CALL mpimerge3d(zpsoil,nx,ny,nzsoil,out3d)
      IF (myproc == 0) &
        CALL hdfwrt3d(out3d,nxlg,nylg,nzsoil,sd_id,0,hdfcompr,          &
                  'zpsoil','Physical height coordinate (soil)','m',     &
                   itmpsoil,hmaxsoil,hminsoil)
    END IF

  END IF

!-----------------------------------------------------------------------
!
!  If basout=1, write out base state variables.
!
!-----------------------------------------------------------------------

  IF(basout == 1 .OR. grdbas == 1 ) THEN

    CALL edgfill(ubar,1,nx,1,nx, 1,ny,1,ny-1, 1,nz,1,nz-1)
    CALL mpimerge3d(ubar,nx,ny,nz,out3d)
    IF (myproc == 0) THEN
      CALL hdfwrt3d(out3d,nxlg,nylg,nz,sd_id,1,hdfcompr,'ubar',         &
                    'Base state u-velocity','m/s',itmp,hmax,hmin)
    END IF

    CALL edgfill(vbar,1,nx,1,nx-1, 1,ny,1,ny, 1,nz,1,nz-1)
    CALL mpimerge3d(vbar,nx,ny,nz,out3d)
    IF (myproc == 0) THEN
      CALL hdfwrt3d(out3d,nxlg,nylg,nz,sd_id,2,hdfcompr,'vbar',         &
                    'Base state v-velocity','m/s',itmp,hmax,hmin)
    END IF

    IF (myproc == 0) THEN
      out3d = 0.0
      CALL hdfwrt3d(out3d,nxlg,nylg,nz,sd_id,3,hdfcompr,'wbar',         &
                    'Base state w-velocity','m/s',itmp,hmax,hmin)
    END IF

    CALL edgfill(ptbar,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)
    CALL mpimerge3d(ptbar,nx,ny,nz,out3d)
    IF (myproc == 0) THEN
      CALL hdfwrt3d(out3d,nxlg,nylg,nz,sd_id,0,hdfcompr,'ptbar',        &
                'Base state potential temperature','K',itmp,hmax,hmin)
    END IF

    CALL edgfill(pbar,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)
    CALL mpimerge3d(pbar,nx,ny,nz,out3d)
    IF (myproc == 0) THEN
      CALL hdfwrt3d(out3d,nxlg,nylg,nz,sd_id,0,hdfcompr,'pbar',         &
                    'Base state pressure','Pascal',itmp,hmax,hmin)
    END IF

    IF(mstout == 1) THEN

      CALL edgfill(qvbar,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)
      CALL mpimerge3d(qvbar,nx,ny,nz,out3d)
      IF (myproc == 0) THEN
        CALL hdfwrt3d(out3d,nxlg,nylg,nz,sd_id,0,hdfcompr,'qvbar',      &
                    'Base state water vapor specific humidity','kg/kg', &
                    itmp,hmax,hmin)
      END IF
    END IF

    IF(landout == 1) THEN

      CALL iedgfill(soiltyp(1,1,1),1,nx,1,nx-1, 1,ny,1,ny-1,            &
                    1,nstypout,1,nstypout)
      CALL mpimerge3di(soiltyp,nx,ny,nstypout,out3di)
      IF (myproc == 0) THEN
        CALL hdfwrt3di(out3di,nxlg,nylg,nstypout,sd_id,0,0,             &
                'soiltyp','Soil type','index')
      END IF

      CALL edgfill(stypfrct(1,1,1),1,nx,1,nx-1, 1,ny,1,ny-1,            &
                   1,nstypout,1,nstypout)
      CALL mpimerge3d(stypfrct,nx,ny,nstypout,out3d)
      IF (myproc == 0) THEN
        CALL hdfwrt3d(out3d,nxlg,nylg,nstypout,sd_id,0,hdfcompr,        &
                'stypfrct','Soil type fractional coverage','fraction',  &
                itmp,hmax,hmin)
      END IF

      CALL iedgfill(vegtyp ,1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)
      CALL mpimerge2di(vegtyp,nx,ny,out3di)
      IF (myproc == 0) THEN
        CALL hdfwrt2di(out3di,nxlg,nylg,sd_id,0,0,'vegtyp',             &
                       'Vegetation type','index')
      END IF

      CALL edgfill(lai    ,1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)
      CALL mpimerge2d(lai,nx,ny,out3d)
      IF (myproc == 0) THEN
        CALL hdfwrt2d(out3d,nxlg,nylg,sd_id,0,hdfcompr,'lai',           &
                      'Leaf Area Index','index',itmp)
      END IF

      CALL edgfill(roufns ,1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)
      CALL mpimerge2d(roufns,nx,ny,out3d)
      IF (myproc == 0) THEN
        CALL hdfwrt2d(out3d,nxlg,nylg,sd_id,0,hdfcompr,'roufns',        &
                      'Surface roughness','0-1',itmp)
      END IF

      CALL edgfill(veg    ,1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)
      CALL mpimerge2d(veg,nx,ny,out3d)
      IF (myproc == 0) THEN
        CALL hdfwrt2d(out3d,nxlg,nylg,sd_id,0,hdfcompr,'veg',           &
                      'Vegetation fraction','fraction',itmp)
      END IF

    END IF

  END IF

  IF ( grdbas == 1 ) GO TO 600

!-----------------------------------------------------------------------
!
!  If varout = 1, Write out uprt, vprt, wprt, ptprt, pprt.
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!  Write out u,v and w.
!
!-----------------------------------------------------------------------

  IF(varout == 1) THEN

    IF ( totout == 0 ) THEN

!-----------------------------------------------------------------------
!
!  Write out perturbations to history dump
!
!-----------------------------------------------------------------------

      DO k=1,nz-1
        DO j=1,ny-1
          DO i=1,nx
            tem1(i,j,k)=u(i,j,k)-ubar(i,j,k)
          END DO
        END DO
      END DO
      CALL edgfill(tem1,1,nx,1,nx, 1,ny,1,ny-1, 1,nz,1,nz-1)
      CALL mpimerge3d(tem1,nx,ny,nz,out3d)
      IF (myproc == 0) THEN
        CALL hdfwrt3d(out3d,nxlg,nylg,nz,sd_id,1,hdfcompr,'uprt',       &
                      'Perturbation u-velocity','m/s',itmp,hmax,hmin)
      END IF

      DO k=1,nz-1
        DO j=1,ny
          DO i=1,nx-1
            tem1(i,j,k)=v(i,j,k)-vbar(i,j,k)
          END DO
        END DO
      END DO
      CALL edgfill(tem1,1,nx,1,nx-1, 1,ny,1,ny, 1,nz,1,nz-1)
      CALL mpimerge3d(tem1,nx,ny,nz,out3d)
      IF (myproc == 0) THEN
        CALL hdfwrt3d(out3d,nxlg,nylg,nz,sd_id,2,hdfcompr,'vprt',       &
                      'Perturbation v-velocity','m/s',itmp,hmax,hmin)
      END IF

      CALL edgfill(w,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz)
      CALL mpimerge3d(w,nx,ny,nz,out3d)
      IF (myproc == 0) THEN
        CALL hdfwrt3d(out3d,nxlg,nylg,nz,sd_id,3,hdfcompr,'wprt',       &
                      'Perturbation w-velocity','m/s',itmp,hmax,hmin)
      END IF

!-----------------------------------------------------------------------
!
!  Write out scalars
!
!-----------------------------------------------------------------------

      CALL edgfill(ptprt,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)
      CALL mpimerge3d(ptprt,nx,ny,nz,out3d)
      IF (myproc == 0) THEN
        CALL hdfwrt3d(out3d,nxlg,nylg,nz,sd_id,0,hdfcompr,'ptprt',      &
                'Perturbation potential temperature','K',itmp,hmax,hmin)
      END IF

      CALL edgfill(pprt,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)
      CALL mpimerge3d(pprt,nx,ny,nz,out3d)
      IF (myproc == 0) THEN
        CALL hdfwrt3d(out3d,nxlg,nylg,nz,sd_id,0,hdfcompr,'pprt',       &
                      'Perturbation pressure','Pascal',itmp,hmax,hmin)
      END IF

    ELSE

!-----------------------------------------------------------------------
!
!  Write out total values to history dump
!
!-----------------------------------------------------------------------

      CALL edgfill(u,1,nx,1,nx, 1,ny,1,ny-1, 1,nz,1,nz-1)
      CALL mpimerge3d(u,nx,ny,nz,out3d)
      IF (myproc == 0) THEN
        CALL hdfwrt3d(out3d,nxlg,nylg,nz,sd_id,1,hdfcompr,'u',          &
                      'u-velocity','m/s',itmp,hmax,hmin)
      END IF

      CALL edgfill(v,1,nx,1,nx-1, 1,ny,1,ny, 1,nz,1,nz-1)
      CALL mpimerge3d(v,nx,ny,nz,out3d)
      IF (myproc == 0) THEN
        CALL hdfwrt3d(out3d,nxlg,nylg,nz,sd_id,2,hdfcompr,'v',          &
                      'v-velocity','m/s',itmp,hmax,hmin)
      END IF

      CALL edgfill(w,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz)
      CALL mpimerge3d(w,nx,ny,nz,out3d)
      IF (myproc == 0) THEN
        CALL hdfwrt3d(out3d,nxlg,nylg,nz,sd_id,3,hdfcompr,'w',          &
                      'w-velocity','m/s',itmp,hmax,hmin)
      END IF

!-----------------------------------------------------------------------
!
!  Write out scalars
!
!-----------------------------------------------------------------------

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
        CALL hdfwrt3d(out3d,nxlg,nylg,nz,sd_id,0,hdfcompr,'pt',         &
                      'Potential temperature','K',itmp,hmax,hmin)
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
        CALL hdfwrt3d(out3d,nxlg,nylg,nz,sd_id,0,hdfcompr,'p',          &
                      'Pressure','Pascal',itmp,hmax,hmin)
      END IF

    END IF

  END IF     ! varout

!-----------------------------------------------------------------------
!
!  If mstout = 1, write out moisture scalars.
!
!-----------------------------------------------------------------------

  IF(mstout == 1) THEN

    IF( totout == 0 ) THEN

!-----------------------------------------------------------------------
!
!  Write out perturbation to history dump
!
!-----------------------------------------------------------------------

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
        CALL hdfwrt3d(out3d,nxlg,nylg,nz,sd_id,0,hdfcompr,'qvprt',      &
           'Pert. water vapor specific humidity','kg/kg',itmp,hmax,hmin)
      END IF

    ELSE

!-----------------------------------------------------------------------
!
!  Write out total values to history dump
!
!-----------------------------------------------------------------------

      CALL edgfill(qv,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)
      CALL mpimerge3d(qv,nx,ny,nz,out3d)
      IF (myproc == 0) THEN
        CALL hdfwrt3d(out3d,nxlg,nylg,nz,sd_id,0,hdfcompr,'qv',         &
                 'Water vapor specific humidity','kg/kg',itmp,hmax,hmin)
      END IF

    END IF

    DO nq = 1,nscalar
      ! DTD: Turn off bit-packing for Z array if hdfcompr > 3 and mphyopt == 11

      hdfcomprtmp = hdfcompr
      IF(nq >= 13 .and. mphyopt == 11 .and. hdfcompr > 3) THEN
        hdfcomprtmp = hdfcompr - 3
      END IF
      CALL edgfill(qscalar(:,:,:,nq),1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)
      CALL mpimerge3d(qscalar(:,:,:,nq),nx,ny,nz,out3d)
      IF (myproc == 0) THEN
        CALL hdfwrt3d(out3d,nxlg,nylg,nz,sd_id,0,hdfcomprtmp,           &
               TRIM(qnames(nq)),TRIM(qdescp(nq)),'kg/kg',itmp,hmax,hmin)
      END IF
    END DO

    IF(rainout == 1) THEN

      CALL edgfill(raing,   1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)
      CALL mpimerge2d(raing,nx,ny,out3d)
      IF (myproc == 0) THEN
        CALL hdfwrt2d(out3d(:,:,1),nxlg,nylg,sd_id,0,hdfcompr,'raing',  &
                      'Grid supersaturation rain','mm',itmp)
      END IF

      CALL edgfill(rainc,1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)
      CALL mpimerge2d(rainc,nx,ny,out3d)
      IF (myproc == 0) THEN
        CALL hdfwrt2d(out3d(:,:,1),nxlg,nylg,sd_id,0,hdfcompr,'rainc',  &
                      'Cumulus convective rain','mm',itmp)
      END IF

    END IF   !rainout

    IF ( prcout == 1 ) THEN
      CALL edgfill(prcrate,1,nx,1,nx-1, 1,ny,1,ny-1, 1,4,1,4)
      CALL mpimerge3d(prcrate,nx,ny,4,out3d)
      IF (myproc == 0) THEN
        CALL hdfwrt2d(out3d(1,1,1),nxlg,nylg,sd_id,0,hdfcompr,          &
                'prcrate1','Total precip. rate','kg/(m**2*s)',itmp)
        CALL hdfwrt2d(out3d(1,1,2),nxlg,nylg,sd_id,0,hdfcompr,          &
                'prcrate2','Grid scale precip. rate','kg/(m**2*s)',itmp)
        CALL hdfwrt2d(out3d(1,1,3),nxlg,nylg,sd_id,0,hdfcompr,          &
                'prcrate3','Cumulative precip. rate','kg/(m**2*s)',itmp)
        CALL hdfwrt2d(out3d(1,1,4),nxlg,nylg,sd_id,0,hdfcompr,          &
                'prcrate4','Microphysics precip. rate','kg/(m**2*s)',itmp)
      END IF
    END IF

  END IF   !mstout

!-----------------------------------------------------------------------
!
!  If tkeout = 1, write out tke.
!
!-----------------------------------------------------------------------

  IF( tkeout == 1 ) THEN

    CALL edgfill(tke,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)
    CALL mpimerge3d(tke,nx,ny,nz,out3d)
    IF (myproc == 0) THEN
      CALL hdfwrt3d(out3d,nxlg,nylg,nz,sd_id,0,hdfcompr,'tke',          &
                    'Turbulent Kinetic Energy','(m/s)**2',itmp,hmax,hmin)
    END IF

  END IF

!-----------------------------------------------------------------------
!
!  If trbout = 1, write out the turbulence parameter, km.
!
!-----------------------------------------------------------------------

  IF( trbout == 1 ) THEN

    CALL edgfill(kmh,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)
    CALL mpimerge3d(kmh,nx,ny,nz,out3d)
    IF (myproc == 0) THEN
      CALL hdfwrt3d(out3d,nxlg,nylg,nz,sd_id,0,hdfcompr,'kmh',          &
               'Hori. turb. mixing coef. for momentum','m**2/s',        &
               itmp,hmax,hmin)
    END IF

    CALL edgfill(kmv,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)
    CALL mpimerge3d(kmv,nx,ny,nz,out3d)
    IF (myproc == 0) THEN
      CALL hdfwrt3d(out3d,nxlg,nylg,nz,sd_id,0,hdfcompr,'kmv',          &
               'Vert. turb. mixing coef. for momentum','m**2/s',        &
               itmp,hmax,hmin)
    END IF

  END IF   ! trbout


!-----------------------------------------------------------------------
!
!  If sfcout = 1, write out the surface variables,
!  tsoil, qsoil, and wetcanp.
!
!-----------------------------------------------------------------------

  IF( sfcout == 1) THEN

    DO is=0,nstypout

      CALL edgfill(tsoil(1,1,1,is),  1,nx,1,nx-1, 1,ny,1,ny-1,        &
                   1,nzsoil,1,nzsoil)
      CALL edgfill(qsoil(1,1,1,is), 1,nx,1,nx-1, 1,ny,1,ny-1,         &
                   1,nzsoil,1,nzsoil)
      CALL edgfill(wetcanp(1,1,is),1,nx,1,nx-1, 1,ny,1,ny-1,          &
                   1,1,1,1)
    END DO

    CALL mpimerge4d(tsoil,nx,ny,nzsoil,nstypout+1,outtsoil)
    CALL mpimerge4d(qsoil,nx,ny,nzsoil,nstypout+1,outqsoil)
    CALL mpimerge3d(wetcanp,nx,ny,nstypout+1,out3d)

    IF (myproc == 0) THEN

      CALL hdfwrt4d(outtsoil,nxlg,nylg,nzsoil,nstypout+1,sd_id,0,       &
                  hdfcompr,'tsoil','Soil temperature','K',              &
                  itmpsoil,hmaxsoil,hminsoil)
      CALL hdfwrt4d(outqsoil,nxlg,nylg,nzsoil,nstypout+1,sd_id,0,       &
                  hdfcompr,'qsoil','Soil moisture','fraction',          &
                  itmpsoil,hmaxsoil,hminsoil)

      CALL hdfwrt3d(out3d,nxlg,nylg,nstypout+1,sd_id,0,hdfcompr,        &
                 'wetcanp','Canopy water amount','fraction',            &
                 itmp,hmax,hmin)

    END IF

    IF (snowout == 1) THEN

      CALL edgfill(snowdpth,1,nx,1,nx-1, 1,ny,1,ny-1,1,1,1,1)
      CALL mpimerge2d(snowdpth,nx,ny,out3d)
      IF (myproc == 0) THEN
        CALL hdfwrt2d(out3d,nxlg,nylg,sd_id,0,hdfcompr,'snowdpth',      &
                      'Snow depth','m',itmp)
      END IF
    END IF

  END IF

!-----------------------------------------------------------------------
!
!  If radout = 1, write out the radiation arrays
!
!-----------------------------------------------------------------------

  IF( radout == 1 ) THEN

    CALL edgfill(radfrc,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)
    CALL mpimerge3d(radfrc,nx,ny,nz,out3d)
    IF (myproc == 0) THEN
      CALL hdfwrt3d(out3d,nxlg,nylg,nz,sd_id,0,hdfcompr,'radfrc',       &
                    'Radiation forcing','K/s',itmp,hmax,hmin)
    END IF

    CALL edgfill(radsw,1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)
    CALL mpimerge2d(radsw,nx,ny,out3d)
    IF (myproc == 0) THEN
      CALL hdfwrt2d(out3d,nxlg,nylg,sd_id,0,hdfcompr,'radsw',           &
                  'Solar radiation reaching the surface','W/m**2',itmp)
    END IF

    CALL edgfill(rnflx,1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)
    CALL mpimerge2d(rnflx,nx,ny,out3d)
    IF (myproc == 0) THEN
      CALL hdfwrt2d(out3d,nxlg,nylg,sd_id,0,hdfcompr,'rnflx',           &
                'Net radiation flux absorbed by surface','W/m**2',itmp)
    END IF

    IF (intver >= intver500) THEN

      CALL edgfill(radswnet,1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)
      CALL mpimerge2d(radswnet,nx,ny,out3d)
      IF (myproc == 0) THEN
        CALL hdfwrt2d(out3d,nxlg,nylg,sd_id,0,hdfcompr,'radswnet',      &
                      'Net solar radiation','W/m**2',itmp)
      END IF

      CALL edgfill(radlwin,1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)
      CALL mpimerge2d(radlwin,nx,ny,out3d)
      IF (myproc == 0) THEN
        CALL hdfwrt2d(out3d,nxlg,nylg,sd_id,0,hdfcompr,'radlwin',       &
                      'Incoming longwave radiation','W/m**2',itmp)
      END IF

    END IF

  END IF   ! radout

!-----------------------------------------------------------------------
!
!  If flxout = 1, write out the surface fluxes
!
!-----------------------------------------------------------------------

  IF( flxout == 1 ) THEN

    CALL edgfill(usflx,1,nx,1,nx, 1,ny,1,ny-1, 1,1,1,1)
    CALL mpimerge2d(usflx,nx,ny,out3d)
    IF (myproc == 0) THEN
      CALL hdfwrt2d(out3d,nxlg,nylg,sd_id,0,hdfcompr,'usflx',           &
                    'Surface flux of u-momentum','kg/(m*s**2)',itmp)
    END IF

    CALL edgfill(vsflx,1,nx,1,nx-1, 1,ny,1,ny, 1,1,1,1)
    CALL mpimerge2d(vsflx,nx,ny,out3d)
    IF (myproc == 0) THEN
      CALL hdfwrt2d(out3d,nxlg,nylg,sd_id,0,hdfcompr,'vsflx',           &
                    'Surface flux of v-momentum','kg/(m*s**2)',itmp)
    END IF

    CALL edgfill(ptsflx,1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)
    CALL mpimerge2d(ptsflx,nx,ny,out3d)
    IF (myproc == 0) THEN
      CALL hdfwrt2d(out3d,nxlg,nylg,sd_id,0,hdfcompr,'ptsflx',          &
                    'Surface heat flux','K*kg/(m**2*s)',itmp)
    END IF

    CALL edgfill(qvsflx,1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)
    CALL mpimerge2d(qvsflx,nx,ny,out3d)
    IF (myproc == 0) THEN
      CALL hdfwrt2d(out3d,nxlg,nylg,sd_id,0,hdfcompr,'qvsflx',          &
                    'Surface moisture flux','kg/(m**2*s)',itmp)
    END IF

  END IF   ! flxout

  600   CONTINUE

  IF (myproc == 0) THEN
    CALL hdfclose(sd_id,istat)
    IF (istat /= 0) WRITE (6,'(1x,3a,I3,a)')      &
      'HDFJOINDUMP: ERROR on closing file ',trim(filename),' (status',istat,')'
  END IF

  IF (hdfcompr > 3) THEN
    DEALLOCATE (itmp,stat=istat)
    DEALLOCATE (hmax,stat=istat)
    DEALLOCATE (hmin,stat=istat)
    DEALLOCATE (itmpsoil,stat=istat)
    DEALLOCATE (hmaxsoil,stat=istat)
    DEALLOCATE (hminsoil,stat=istat)
  END IF

  DEALLOCATE (out1d)
  DEALLOCATE (out3d, out3di)
  DEALLOCATE (outtsoil, outqsoil)

  RETURN
END SUBROUTINE hdfjoindump

!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE HDFSPLITDUMP               ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################

SUBROUTINE hdfsplitdump(nx,ny,nz,nzsoil,nstyps, filename , grdbas,      &
           u,v,w,ptprt,pprt,qv,qscalar,tke,kmh,kmv,                     &
           ubar,vbar,ptbar,pbar,rhobar,qvbar,                           &
           x,y,z,zp,zpsoil,                                             &
           soiltyp,stypfrct,vegtyp,lai,roufns,veg,                      &
           tsoil,qsoil,wetcanp,snowdpth,                                &
           raing,rainc,prcrate,                                         &
           radfrc,radsw,rnflx,radswnet,radlwin,                         &
           usflx,vsflx,ptsflx,qvsflx,                                   &
           tem1)

!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Produces a history data file "filename" in the NCSA HDF4 format by
!  calling HDF library subroutines.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Gene Bassett
!  2000/03/15
!
!  MODIFICATION HISTORY:
!
!  05/15/2002  (J. Brotzge)
!  Added to allow for multiple soil schemes
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!    nzsoil   Number of grid points in the vertical
!
!    filename File name of history dump data.
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
!    tem1     work array.
!    tem2     work array.
!
!-----------------------------------------------------------------------
!
!  The following parameters are passed into this subroutine through
!  a common block in globcst.inc, and they determine which
!  variables are output.
!
!  grdout =0 or 1. If grdout=0, grid variables are not dumped.
!  basout =0 or 1. If basout=0, base state variables are not dumped.
!  varout =0 or 1. If varout=0, perturbation variables are not dumped.
!  mstout =0 or 1. If mstout=0, water variables are not dumped.
!  rainout=0 or 1. If rainout=0, rain variables are not dumped.
!  prcout =0 or 1. If prcout=0, precipitation rates are not dumped.
!  iceout =0 or 1. If iceout=0, qi, qs and qh are not dumped.
!  tkeout =0 or 1. If tkeout=0, tke is not dumped.
!  trbout =0 or 1. If trbout=0, the eddy viscosity km is not dumped.
!  radout =0 or 1. If radout=0, the radiation arrays are not dumped.
!  flxout =0 or 1. If flxout=0, the surface fluxes are not dumped.
!
!-----------------------------------------------------------------------

  USE arps_precision

  IMPLICIT NONE

!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------

  INCLUDE 'globcst.inc'
  INCLUDE 'grid.inc'        ! Grid parameters
  INCLUDE 'phycst.inc'
  INCLUDE 'mp.inc'

!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------

  INTEGER, INTENT(IN) :: nx,ny,nz          ! Number of grid points in 3 directions
  INTEGER, INTENT(IN) :: nzsoil            ! Number of grid points in the soil

  CHARACTER (LEN=*), INTENT(IN) :: filename
  INTEGER, INTENT(IN) :: grdbas            ! If this is a grid/base state dump

  REAL, INTENT(IN) :: u     (nx,ny,nz)     ! Total u-velocity (m/s)
  REAL, INTENT(IN) :: v     (nx,ny,nz)     ! Total v-velocity (m/s)
  REAL, INTENT(IN) :: w     (nx,ny,nz)     ! Total w-velocity (m/s)
  REAL, INTENT(IN) :: ptprt (nx,ny,nz)     ! Perturbation potential temperature (K)
  REAL, INTENT(IN) :: pprt  (nx,ny,nz)     ! Perturbation pressure (Pascal)

  REAL, INTENT(IN) :: qv    (nx,ny,nz)     ! Water vapor specific humidity (kg/kg)
  REAL, INTENT(IN) :: qscalar(nx,ny,nz,nscalar)
  REAL, INTENT(IN) :: tke   (nx,ny,nz)     ! Turbulent Kinetic Energy ((m/s)**2)
  REAL, INTENT(IN) :: kmh   (nx,ny,nz)     ! Horizontal turb. mixing coef. for
                               ! momentum. ( m**2/s )
  REAL, INTENT(IN) :: kmv   (nx,ny,nz)     ! Vertical turb. mixing coef. for
                               ! momentum. ( m**2/s )

  REAL, INTENT(IN) :: ubar  (nx,ny,nz)     ! Base state u-velocity (m/s)
  REAL, INTENT(IN) :: vbar  (nx,ny,nz)     ! Base state v-velocity (m/s)
  REAL, INTENT(IN) :: ptbar (nx,ny,nz)     ! Base state potential temperature (K)
  REAL, INTENT(IN) :: pbar  (nx,ny,nz)     ! Base state pressure (Pascal)
  REAL, INTENT(IN) :: rhobar(nx,ny,nz)     ! Base state air density (kg/m**3)
  REAL, INTENT(IN) :: qvbar (nx,ny,nz)     ! Base state water vapor specific humidity
                               ! (kg/kg)

  REAL, INTENT(IN) :: x     (nx)           ! The x-coord. of the physical and
                               ! computational grid. Defined at u-point.
  REAL, INTENT(IN) :: y     (ny)           ! The y-coord. of the physical and
                               ! computational grid. Defined at v-point.
  REAL, INTENT(IN) :: z     (nz)           ! The z-coord. of the computational grid.
                               ! Defined at w-point on the staggered grid.
  REAL, INTENT(IN) :: zp    (nx,ny,nz)     ! The physical height coordinate defined at
                               ! w-point of the staggered grid.
  REAL, INTENT(IN) :: zpsoil(nx,ny,nzsoil) ! The physical height coordinate defined at
                               ! w-point of the soil.

  INTEGER, INTENT(IN) :: nstyps             ! Number of soil type

  INTEGER, INTENT(IN) :: soiltyp (nx,ny,nstyps)    ! Soil type
  REAL   , INTENT(IN) :: stypfrct(nx,ny,nstyps)    ! Fraction of soil types
  INTEGER, INTENT(IN) :: vegtyp (nx,ny)            ! Vegetation type
  REAL, INTENT(IN) :: lai    (nx,ny)    ! Leaf Area Index
  REAL, INTENT(IN) :: roufns (nx,ny)    ! Surface roughness
  REAL, INTENT(IN) :: veg    (nx,ny)    ! Vegetation fraction

  REAL, INTENT(IN) :: tsoil (nx,ny,nzsoil,0:nstyps)   ! Soil temperature (K)
  REAL, INTENT(IN) :: qsoil (nx,ny,nzsoil,0:nstyps)   ! Soil moisture (m**3/m**3)
  REAL, INTENT(IN) :: wetcanp(nx,ny,0:nstyps)         ! Canopy water amount
  REAL, INTENT(IN) :: snowdpth(nx,ny)                 ! Snow depth (m)

  REAL, INTENT(IN) :: raing(nx,ny)         ! Grid supersaturation rain
  REAL, INTENT(IN) :: rainc(nx,ny)         ! Cumulus convective rain
  REAL, INTENT(IN) :: prcrate(nx,ny,4)     ! precipitation rates (kg/(m**2*s))
                               ! prcrate(1,1,1) = total precip. rate
                               ! prcrate(1,1,2) = grid scale precip. rate
                               ! prcrate(1,1,3) = cumulative precip. rate
                               ! prcrate(1,1,4) = microphysics precip. rate

  REAL, INTENT(IN) :: radfrc(nx,ny,nz)     ! Radiation forcing (K/s)
  REAL, INTENT(IN) :: radsw (nx,ny)        ! Solar radiation reaching the surface
  REAL, INTENT(IN) :: rnflx (nx,ny)        ! Net radiation flux absorbed by surface
  REAL, INTENT(IN) :: radswnet(nx,ny)      ! Net shortwave radiation
  REAL, INTENT(IN) :: radlwin(nx,ny)       ! Incoming longwave radiation

  REAL, INTENT(IN) :: usflx (nx,ny)        ! Surface flux of u-momentum (kg/(m*s**2))
  REAL, INTENT(IN) :: vsflx (nx,ny)        ! Surface flux of v-momentum (kg/(m*s**2))
  REAL, INTENT(IN) :: ptsflx(nx,ny)        ! Surface heat flux (K*kg/(m**2*s))
  REAL, INTENT(IN) :: qvsflx(nx,ny)        ! Surface moisture flux (kg/(m**2*s))

  REAL :: tem1  (nx,ny,nz)                     ! Temporary work array

!-----------------------------------------------------------------------
!
!  Parameters describing this routine
!
!-----------------------------------------------------------------------
! 06/28/2002 Zuwen He
!
! fmtver??: to label each data a version.
! intver??: an integer to allow faster comparison than fmtver??,
!           which are strings.
!
! Verion 5.00: significant change in soil variables since version 4.10.
!
  CHARACTER (LEN=40) :: fmtver,fmtver410,fmtver500,fmtver530
  INTEGER  :: intver,intver410,intver500,intver530

  PARAMETER (fmtver410='004.10 HDF4 Coded Data',intver410=410)
  PARAMETER (fmtver500='005.00 HDF4 Coded Data',intver500=500)
  PARAMETER (fmtver530='005.30 HDF4 Coded Data',intver530=530)

  CHARACTER (LEN=10) :: tmunit
  PARAMETER (tmunit='seconds   ')

!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------

  INTEGER :: i,j,k,is,nq
  INTEGER :: nstypout
  INTEGER :: istat, sd_id
  INTEGER :: hdfcomprtmp, nztmp

  INTEGER :: npxout, npyout
  INTEGER :: nxout,  nyout
  INTEGER :: ipx,    jpy
  INTEGER :: ia,     ja
  LOGICAL :: OUTFLAG

  INTEGER :: innerx, innery
  INTEGER :: innerx_stag, innery_stag

  INTEGER :: lenbase
  CHARACTER(LEN=256) :: outfilename

  REAL(SP), ALLOCATABLE :: varoutr(:,:,:)
  REAL(SP), ALLOCATABLE :: varout4d(:,:,:,:)
  INTEGER,  ALLOCATABLE :: varouti(:,:,:)

  INTEGER(2), ALLOCATABLE :: itmp(:,:,:)       ! Temporary array
  REAL,       ALLOCATABLE :: hmax(:), hmin(:)  ! Temporary array
  INTEGER(2), ALLOCATABLE :: itmpsoil(:,:,:,:) ! Temporary array
  REAL,       ALLOCATABLE :: hmaxsoil(:), hminsoil(:)! Temporary array

  REAL :: dx_out,dy_out

  CHARACTER(LEN=256) :: cmntstr

!-----------------------------------------------------------------------
!
!  Functions
!
!-----------------------------------------------------------------------

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  npxout = nproc_x_out/nproc_x
  npyout = nproc_y_out/nproc_y

  nxout = (nx-3)/npxout + 3
  nyout = (ny-3)/npyout + 3

  lenbase = INDEX(filename,'_',.TRUE.)

!-----------------------------------------------------------------------

  nstypout = max(1,nstyps)

  nztmp = max(nz,nzsoil,nstypout+1)

  IF (hdfcompr > 3) THEN
!
!   If "nz" is a small number, the arrays may not be big enough for the
!   soil family of variables.
!
    ALLOCATE (itmp(nxout,nyout,nztmp),stat=istat)
    CALL check_alloc_status(istat,'HDFSPLITDUMP:itmp')
    ALLOCATE (hmax(nztmp),stat=istat)
    CALL check_alloc_status(istat,'HDFSPLITDUMP:hmax')
    ALLOCATE (hmin(nztmp),stat=istat)
    CALL check_alloc_status(istat,'HDFSPLITDUMP:hmin')
    ALLOCATE (itmpsoil(nxout,nyout,nzsoil,0:nstyps),stat=istat)
    CALL check_alloc_status(istat,'HDFSPLITDUMP:itmpsoil')
    ALLOCATE (hmaxsoil(nzsoil),stat=istat)
    CALL check_alloc_status(istat,'HDFSPLITDUMP:hmaxsoil')
    ALLOCATE (hminsoil(nzsoil),stat=istat)
    CALL check_alloc_status(istat,'HDFSPLITDUMP:hminsoil')
  END IF

  ALLOCATE(varoutr (nxout,nyout,nztmp),             STAT = istat)
  CALL check_alloc_status(istat,'HDFSPLITDUMP:varoutr')
  ALLOCATE(varouti (nxout,nyout,nstypout+1),        STAT = istat)
  CALL check_alloc_status(istat,'HDFSPLITDUMP:varouti')
  ALLOCATE(varout4d(nxout,nyout,nzsoil,nstypout+1), STAT = istat)
  CALL check_alloc_status(istat,'HDFSPLITDUMP:varout4d')

  IF (myproc == 0) WRITE(6,'(1x,a,f13.3,a,a/)')                         &
       'Writing HDF4 data at time=', curtim,' into file ',filename

!-----------------------------------------------------------------------
!
!  Create the HDF4 file.
!
!-----------------------------------------------------------------------

  DO jpy = 1,npyout
    DO ipx = 1,npxout

      OUTFLAG = (myproc == 0 .AND. ipx == 1 .AND. jpy == 1)

      ia = (ipx-1)*(nxout-3)
      ja = (jpy-1)*(nyout-3)

      innerx = nxout    ! edgefilling index
      innery = nyout
      innerx_stag = nxout
      innery_stag = nyout
      IF (loc_x == nproc_x .AND. ipx == npxout) THEN
        innerx = nxout-1
      END IF
      IF (loc_y == nproc_y .AND. jpy == npyout) THEN
        innery = nyout-1
      END IF

      CALL gtsplitfn(filename(1:lenbase-1),npxout,npyout,loc_x,loc_y,ipx,jpy,  &
                     0,0,0,lvldbg,outfilename,istat)

      CALL hdfopen(TRIM(outfilename),2,sd_id)
      IF (sd_id < 0) THEN
        IF(OUTFLAG) WRITE (6,*) "HDFSPLITDUMP: ERROR creating HDF4 file: ",  &
                    trim(outfilename)
        GO TO 600
      END IF

      intver = intver500  !  for the time being, in the future, we will
                          !  allow to dump data in the different version
                          !  intver will be assigned from input file

      fmtver=fmtver500

      CALL hdfwrtc(sd_id, 40, 'fmtver', fmtver, istat)
      CALL hdfwrtc(sd_id, 80, 'runname', runname, istat)
      CALL hdfwrti(sd_id, 'nocmnt', nocmnt, istat)

      IF( nocmnt > 0 ) THEN
        IF (nocmnt*80 > 256) THEN
          WRITE(*,'(1x,a,I0,a)') "ERROR: comment array is ",80*nocmnt,    &
                               ", but comment str is only 256 characters."
          CALL arpsstop('ERROR:Comment string is too short.',1)
        END IF
        DO i = 1,nocmnt
          cmntstr((i-1)*80+1:i*80) = cmnt(i)
        END DO
        CALL hdfwrtc(sd_id, 80*nocmnt, 'cmnt', cmntstr, istat)
      END IF

      CALL hdfwrtc(sd_id, 7, 'tmunit', 'seconds', istat)
      CALL hdfwrtr(sd_id, 'time', curtim, istat)

      CALL hdfwrti(sd_id, 'nx', nxout, istat)
      CALL hdfwrti(sd_id, 'ny', nyout, istat)
      CALL hdfwrti(sd_id, 'nz', nz, istat)

      IF (intver >= intver500)  THEN
        CALL hdfwrti(sd_id, 'nzsoil', nzsoil, istat)
      END IF

      IF( grdbas == 1 ) THEN

        CALL hdfwrti(sd_id, 'grdflg', 1, istat)
        CALL hdfwrti(sd_id, 'basflg', 1, istat)
        CALL hdfwrti(sd_id, 'varflg', 0, istat)
        CALL hdfwrti(sd_id, 'mstflg', 1, istat)
        CALL hdfwrti(sd_id, 'iceflg', 0, istat)
        CALL hdfwrti(sd_id, 'trbflg', 0, istat)
        CALL hdfwrti(sd_id, 'sfcflg', 0, istat)
        CALL hdfwrti(sd_id, 'rainflg', 0, istat)
        CALL hdfwrti(sd_id, 'landflg', landout, istat)
        CALL hdfwrti(sd_id, 'totflg', 1, istat)
        CALL hdfwrti(sd_id, 'tkeflg', 0, istat)

      ELSE

        CALL hdfwrti(sd_id, 'grdflg', grdout, istat)
        CALL hdfwrti(sd_id, 'basflg', basout, istat)
        CALL hdfwrti(sd_id, 'varflg', varout, istat)
        CALL hdfwrti(sd_id, 'mstflg', mstout, istat)
        CALL hdfwrti(sd_id, 'iceflg', iceout, istat)
        CALL hdfwrti(sd_id, 'trbflg', trbout, istat)
        CALL hdfwrti(sd_id, 'sfcflg', sfcout, istat)
        CALL hdfwrti(sd_id, 'rainflg', rainout, istat)
        CALL hdfwrti(sd_id, 'landflg', landout*basout, istat)
        CALL hdfwrti(sd_id, 'totflg', totout, istat)
        CALL hdfwrti(sd_id, 'tkeflg', tkeout, istat)

        CALL hdfwrti(sd_id, 'nscalar', nscalar, istat)
        CALL hdfwrti(sd_id, 'P_QC',    P_QC, istat)
        CALL hdfwrti(sd_id, 'P_QR',    P_QR, istat)
        CALL hdfwrti(sd_id, 'P_QI',    P_QI, istat)
        CALL hdfwrti(sd_id, 'P_QS',    P_QS, istat)
        CALL hdfwrti(sd_id, 'P_QG',    P_QG, istat)
        CALL hdfwrti(sd_id, 'P_QH',    P_QH, istat)

        CALL hdfwrti(sd_id, 'P_NC',    P_NC, istat)
        CALL hdfwrti(sd_id, 'P_NR',    P_NR, istat)
        CALL hdfwrti(sd_id, 'P_NI',    P_NI, istat)
        CALL hdfwrti(sd_id, 'P_NS',    P_NS, istat)
        CALL hdfwrti(sd_id, 'P_NG',    P_NG, istat)
        CALL hdfwrti(sd_id, 'P_NH',    P_NH, istat)

        CALL hdfwrti(sd_id, 'P_ZR',    P_ZR, istat)
        CALL hdfwrti(sd_id, 'P_ZI',    P_ZI, istat)
        CALL hdfwrti(sd_id, 'P_ZS',    P_ZS, istat)
        CALL hdfwrti(sd_id, 'P_ZG',    P_ZG, istat)
        CALL hdfwrti(sd_id, 'P_ZH',    P_ZH, istat)

        CALL hdfwrti(sd_id, 'P_CC',    P_CC, istat)

      END IF

      nstypout = max(1,nstyps)
      CALL hdfwrti(sd_id, 'nstyp', nstypout, istat)
      CALL hdfwrti(sd_id, 'prcflg', prcout, istat)
      CALL hdfwrti(sd_id, 'radflg', radout, istat)
      CALL hdfwrti(sd_id, 'flxflg', flxout, istat)
      CALL hdfwrti(sd_id, 'snowflg', snowout, istat)

      CALL hdfwrti(sd_id, 'day', day, istat)
      CALL hdfwrti(sd_id, 'year', year, istat)
      CALL hdfwrti(sd_id, 'month', month, istat)
      CALL hdfwrti(sd_id, 'hour', hour, istat)
      CALL hdfwrti(sd_id, 'minute', minute, istat)
      CALL hdfwrti(sd_id, 'second', second, istat)

      CALL hdfwrtr(sd_id, 'umove', umove, istat)
      CALL hdfwrtr(sd_id, 'vmove', vmove, istat)
      CALL hdfwrtr(sd_id, 'xgrdorg', xgrdorg, istat)
      CALL hdfwrtr(sd_id, 'ygrdorg', ygrdorg, istat)

      CALL hdfwrti(sd_id, 'mapproj', mapproj, istat)
      CALL hdfwrtr(sd_id, 'trulat1', trulat1, istat)
      CALL hdfwrtr(sd_id, 'trulat2', trulat2, istat)
      CALL hdfwrtr(sd_id, 'trulon', trulon, istat)
      CALL hdfwrtr(sd_id, 'sclfct', sclfct, istat)
      CALL hdfwrtr(sd_id, 'tstop', tstop, istat)
      CALL hdfwrtr(sd_id, 'thisdmp', thisdmp, istat)
      CALL hdfwrtr(sd_id, 'latitud', latitud, istat)
      CALL hdfwrtr(sd_id, 'ctrlat', ctrlat, istat)
      CALL hdfwrtr(sd_id, 'ctrlon', ctrlon, istat)

      CALL hdfwrtr(sd_id, 'ntcloud', ntcloud, istat)
      CALL hdfwrtr(sd_id, 'n0rain',  n0rain,  istat)
      CALL hdfwrtr(sd_id, 'n0snow',  n0snow,  istat)
      CALL hdfwrtr(sd_id, 'n0grpl',  n0grpl,  istat)
      CALL hdfwrtr(sd_id, 'n0hail',  n0hail,  istat)
      CALL hdfwrtr(sd_id, 'rhoice',  rhoice,  istat)
      CALL hdfwrtr(sd_id, 'rhosnow', rhosnow, istat)
      CALL hdfwrtr(sd_id, 'rhogrpl', rhogrpl, istat)
      CALL hdfwrtr(sd_id, 'rhohail', rhohail, istat)

      CALL hdfwrtr(sd_id, 'alpharain', alpharain, istat)
      CALL hdfwrtr(sd_id, 'alphaice',  alphaice,  istat)
      CALL hdfwrtr(sd_id, 'alphasnow', alphasnow, istat)
      CALL hdfwrtr(sd_id, 'alphagrpl', alphagrpl, istat)
      CALL hdfwrtr(sd_id, 'alphahail', alphahail, istat)

      dx_out = x(2) - x(1)
      dy_out = y(2) - y(1)
      CALL hdfwrtr(sd_id, 'dx', dx_out, istat)
      CALL hdfwrtr(sd_id, 'dy', dy_out, istat)

!-----------------------------------------------------------------------
!
!  If grdout=1 or grdbas=1, write out grid variables
!
!-----------------------------------------------------------------------

      IF(grdout == 1 .OR. grdbas == 1 ) THEN

        CALL hdfwrt1d(x(ia+1),nxout,sd_id,'x','x coordinate','m')
        CALL hdfwrt1d(y(ja+1),nyout,sd_id,'y','y coordinate','m')
        CALL hdfwrt1d(z,nz,sd_id,'z','z coordinate','m')

        DO k = 1,nz
          DO j = 1, nyout
            DO i = 1,nxout
              varoutr(i,j,k) = zp(i+ia,j+ja,k)
            END DO
          END DO
        END DO
        CALL hdfwrt3d(varoutr,nxout,nyout,nz,sd_id,0,hdfcompr,            &
                      'zp','Physical height coordinate','m',              &
                       itmp,hmax,hmin)

        DO k = 1,nzsoil
          DO j = 1, nyout
            DO i = 1,nxout
              varoutr(i,j,k) = zpsoil(i+ia,j+ja,k)
            END DO
          END DO
        END DO
        CALL hdfwrt3d(varoutr,nxout,nyout,nzsoil,sd_id,0,hdfcompr,        &
                      'zpsoil','Physical height coordinate (soil)','m',   &
                       itmpsoil,hmaxsoil,hminsoil)
      END IF

!-----------------------------------------------------------------------
!
!  If basout=1, write out base state variables.
!
!-----------------------------------------------------------------------

      IF(basout == 1 .OR. grdbas == 1 ) THEN

        DO k = 1,nz
          DO j = 1, nyout
            DO i = 1, nxout
              varoutr(i,j,k) = ubar(i+ia,j+ja,k)
            END DO
          END DO
        END DO
        CALL edgfill(varoutr,1,nxout,1,innerx_stag, 1,nyout,1,innery, 1,nz,1,nz-1)
        CALL hdfwrt3d(varoutr,nxout,nyout,nz,sd_id,1,hdfcompr,              &
                      'ubar','Base state u-velocity','m/s',                 &
                       itmp,hmax,hmin)

        DO k = 1,nz
          DO j = 1, nyout
            DO i = 1,nxout
              varoutr(i,j,k) = vbar(i+ia,j+ja,k)
            END DO
          END DO
        END DO
        CALL edgfill(varoutr,1,nxout,1,innerx, 1,nyout,1,innery_stag, 1,nz,1,nz-1)
        CALL hdfwrt3d(varoutr,nxout,nyout,nz,sd_id,2,hdfcompr,              &
                      'vbar','Base state v-velocity','m/s',                 &
                       itmp,hmax,hmin)

        DO k=1,nz
          DO j=1,nyout
            DO i=1,nxout
              varoutr(i,j,k) = 0.0
            END DO
          END DO
        END DO
        CALL hdfwrt3d(varoutr,nxout,nyout,nz,sd_id,3,hdfcompr,              &
                      'wbar','Base state w-velocity','m/s',                 &
                       itmp,hmax,hmin)

        DO k = 1,nz
          DO j = 1, nyout
            DO i = 1,nxout
              varoutr(i,j,k) = ptbar(i+ia,j+ja,k)
            END DO
          END DO
        END DO
        CALL edgfill(varoutr,1,nxout,1,innerx, 1,nyout,1,innery, 1,nz,1,nz-1)
        CALL hdfwrt3d(varoutr,nxout,nyout,nz,sd_id,0,hdfcompr,              &
                      'ptbar','Base state potential temperature','K',       &
                       itmp,hmax,hmin)

        DO k = 1,nz
          DO j = 1, nyout
            DO i = 1,nxout
              varoutr(i,j,k) = pbar(i+ia,j+ja,k)
            END DO
          END DO
        END DO
        CALL edgfill(varoutr,1,nxout,1,innerx, 1,nyout,1,innery, 1,nz,1,nz-1)
        CALL hdfwrt3d(varoutr,nxout,nyout,nz,sd_id,0,hdfcompr,              &
                      'pbar','Base state pressure','Pascal',                &
                       itmp,hmax,hmin)

        IF(mstout == 1) THEN

          DO k = 1,nz
            DO j = 1, nyout
              DO i = 1,nxout
                varoutr(i,j,k) = qvbar(i+ia,j+ja,k)
              END DO
            END DO
          END DO
          CALL edgfill(varoutr,1,nxout,1,innerx, 1,nyout,1,innery, 1,nz,1,nz-1)
          CALL hdfwrt3d(varoutr,nxout,nyout,nz,sd_id,0,hdfcompr,            &
              'qvbar','Base state water vapor specific humidity','kg/kg',   &
                       itmp,hmax,hmin)

        END IF

        IF(landout == 1) THEN

          DO k = 1,nstypout
            DO j = 1, nyout
              DO i = 1,nxout
                varouti(i,j,k) = soiltyp(i+ia,j+ja,k)
              END DO
            END DO
          END DO
          CALL iedgfill(varouti,1,nxout,1,innerx, 1,nyout,1,innery,     &
                          1,nstypout,1,nstypout)
          CALL hdfwrt3di(varouti,nxout,nyout,nstypout,sd_id,0,0,        &
                    'soiltyp','Soil type','index')

          DO k = 1,nstypout
            DO j = 1, nyout
              DO i = 1,nxout
                varoutr(i,j,k) = stypfrct(i+ia,j+ja,k)
              END DO
            END DO
          END DO
          CALL edgfill(varoutr,1,nxout,1,innerx, 1,nyout,1,innery,      &
                          1,nstypout,1,nstypout)
          CALL hdfwrt3d(varoutr,nxout,nyout,nstypout,sd_id,0,hdfcompr,  &
              'stypfrct','Soil type fractional coverage','fraction',    &
                     itmp,hmax,hmin)

          DO j = 1, nyout
            DO i = 1, nxout
              varouti(i,j,1) = vegtyp(i+ia,j+ja)
            END DO
          END DO
          CALL iedgfill(varouti,1,nxout,1,innerx, 1,nyout,1,innery,1,1,1,1)
          CALL hdfwrt2di(varouti,nxout,nyout,sd_id,0,0,                 &
                      'vegtyp','Vegetation type','index')

          DO j = 1, nyout
            DO i = 1, nxout
              varoutr(i,j,1) = lai(i+ia,j+ja)
            END DO
          END DO
          CALL edgfill(varoutr,1,nxout,1,innerx, 1,nyout,1,innery,1,1,1,1)
          CALL hdfwrt2d(varoutr,nxout,nyout,sd_id,0,hdfcompr,           &
                      'lai','Leaf Area Index','index',itmp)

          DO j = 1, nyout
            DO i = 1, nxout
              varoutr(i,j,1) = roufns(i+ia,j+ja)
            END DO
          END DO
          CALL edgfill(varoutr,1,nxout,1,innerx, 1,nyout,1,innery,1,1,1,1)
          CALL hdfwrt2d(varoutr,nxout,nyout,sd_id,0,hdfcompr,           &
                      'roufns','Surface roughness','0-1',itmp)

          DO j = 1, nyout
            DO i = 1, nxout
              varoutr(i,j,1) = veg(i+ia,j+ja)
            END DO
          END DO
          CALL edgfill(varoutr,1,nxout,1,innerx, 1,nyout,1,innery,1,1,1,1)
          CALL hdfwrt2d(varoutr,nxout,nyout,sd_id,0,hdfcompr,           &
                      'veg','Vegetation fraction','fraction',itmp)

        END IF

      END IF

      IF ( grdbas == 1 ) GO TO 600

!-----------------------------------------------------------------------
!
!  If varout = 1, Write out uprt, vprt, wprt, ptprt, pprt.
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!  Write out u,v and w.
!
!-----------------------------------------------------------------------

      IF(varout == 1) THEN

        IF ( totout == 0 ) THEN

!-----------------------------------------------------------------------
!
!  Write out perturbations to history dump
!
!-----------------------------------------------------------------------

          DO k=1,nz-1
            DO j=1,nyout
              DO i=1, nxout
                varoutr(i,j,k)=u(i+ia,j+ja,k)-ubar(i+ia,j+ja,k)
              END DO
            END DO
          END DO
          CALL edgfill(varoutr,1,nxout,1,innerx_stag, 1,nyout,1,innery, 1,nz,1,nz-1)
          CALL hdfwrt3d(varoutr,nxout,nyout,nz,sd_id,1,hdfcompr,        &
                        'uprt','Perturbation u-velocity','m/s',         &
                         itmp,hmax,hmin)

          DO k=1,nz-1
            DO j=1, nyout
              DO i=1, nxout
                varoutr(i,j,k)=v(i+ia,j+ja,k)-vbar(i+ia,j+ja,k)
              END DO
            END DO
          END DO
          CALL edgfill(varoutr,1,nxout,1,innerx, 1,nyout,1,innery_stag, 1,nz,1,nz-1)
          CALL hdfwrt3d(varoutr,nxout,nyout,nz,sd_id,2,hdfcompr,        &
                        'vprt','Perturbation v-velocity','m/s',         &
                         itmp,hmax,hmin)

          DO k = 1,nz
            DO j = 1, nyout
              DO i = 1, nxout
                varoutr(i,j,k) = w(i+ia,j+ja,k)
              END DO
            END DO
          END DO
          CALL edgfill(varoutr,1,nxout,1,innerx, 1,nyout,1,innery, 1,nz,1,nz)
          CALL hdfwrt3d(varoutr,nxout,nyout,nz,sd_id,3,hdfcompr,         &
                        'wprt','Perturbation w-velocity','m/s',         &
                         itmp,hmax,hmin)

!-----------------------------------------------------------------------
!
!  Write out scalars
!
!-----------------------------------------------------------------------

          DO k = 1,nz
            DO j = 1, nyout
              DO i = 1, nxout
                varoutr(i,j,k) = ptprt(i+ia,j+ja,k)
              END DO
            END DO
          END DO
          CALL edgfill(varoutr,1,nxout,1,innerx, 1,nyout,1,innery, 1,nz,1,nz-1)
          CALL hdfwrt3d(varoutr,nxout,nyout,nz,sd_id,0,hdfcompr,         &
                        'ptprt','Perturbation potential temperature','K',&
                         itmp,hmax,hmin)

          DO k = 1,nz
            DO j = 1, nyout
              DO i = 1, nxout
                varoutr(i,j,k) = pprt(i+ia,j+ja,k)
              END DO
            END DO
          END DO
          CALL edgfill(varoutr,1,nxout,1,innerx, 1,nyout,1,innery, 1,nz,1,nz-1)
          CALL hdfwrt3d(varoutr,nxout,nyout,nz,sd_id,0,hdfcompr,         &
                        'pprt','Perturbation pressure','Pascal',         &
                         itmp,hmax,hmin)

        ELSE

!-----------------------------------------------------------------------
!
!  Write out total values to history dump
!
!-----------------------------------------------------------------------

          DO k = 1,nz
            DO j = 1, nyout
              DO i = 1,nxout
                varoutr(i,j,k) = u(i+ia,j+ja,k)
              END DO
            END DO
          END DO
          CALL edgfill(varoutr,1,nxout,1,innerx_stag, 1,nyout,1,innery, 1,nz,1,nz-1)
          CALL hdfwrt3d(varoutr,nxout,nyout,nz,sd_id,1,hdfcompr,         &
                        'u','u-velocity','m/s',                          &
                         itmp,hmax,hmin)

          DO k = 1,nz
            DO j = 1, nyout
              DO i = 1,nxout
                varoutr(i,j,k) = v(i+ia,j+ja,k)
              END DO
            END DO
          END DO
          CALL edgfill(varoutr,1,nxout,1,innerx, 1,nyout,1,innery_stag, 1,nz,1,nz-1)
          CALL hdfwrt3d(varoutr,nxout,nyout,nz,sd_id,2,hdfcompr,         &
                        'v','v-velocity','m/s',                          &
                         itmp,hmax,hmin)

          DO k = 1,nz
            DO j = 1, nyout
              DO i = 1, nxout
                varoutr(i,j,k) = w(i+ia,j+ja,k)
              END DO
            END DO
          END DO
          CALL edgfill(varoutr,1,nxout,1,innerx, 1,nyout,1,innery, 1,nz,1,nz)
          CALL hdfwrt3d(varoutr,nxout,nyout,nz,sd_id,3,hdfcompr,         &
                        'w','w-velocity','m/s',                          &
                         itmp,hmax,hmin)

!-----------------------------------------------------------------------
!
!  Write out scalars
!
!-----------------------------------------------------------------------

          DO k=1,nz-1
            DO j=1,nyout
              DO i=1,nxout
                varoutr(i,j,k) = ptbar(i+ia,j+ja,k) + ptprt(i+ia,j+ja,k)
              END DO
            END DO
          END DO
          CALL edgfill(varoutr,1,nxout,1,innerx, 1,nyout,1,innery, 1,nz,1,nz-1)
          CALL hdfwrt3d(varoutr,nxout,nyout,nz,sd_id,0,hdfcompr,        &
                        'pt','Potential temperature','K',               &
                        itmp,hmax,hmin)

          DO k=1,nz-1
            DO j=1,nyout
              DO i=1, nxout
                varoutr(i,j,k) = pbar(i+ia,j+ja,k) + pprt(i+ia,j+ja,k)
              END DO
            END DO
          END DO
          CALL edgfill(varoutr,1,nxout,1,innerx, 1,nyout,1,innery, 1,nz,1,nz-1)
          CALL hdfwrt3d(varoutr,nxout,nyout,nz,sd_id,0,hdfcompr,        &
                        'p','Pressure','Pascal',                        &
                        itmp,hmax,hmin)
        END IF

      END IF     ! varout

!-----------------------------------------------------------------------
!
!  If mstout = 1, write out moisture scalars.
!
!-----------------------------------------------------------------------

      IF(mstout == 1) THEN

        IF( totout == 0 ) THEN

!-----------------------------------------------------------------------
!
!  Write out perturbation to history dump
!
!-----------------------------------------------------------------------

          DO k=1,nz-1
            DO j=1,nyout
              DO i=1, nxout
                varoutr(i,j,k)=qv(i+ia,j+ja,k)-qvbar(i+ia,j+ja,k)
              END DO
            END DO
          END DO

          CALL edgfill(varoutr,1,nxout,1,innerx, 1,nyout,1,innery, 1,nz,1,nz-1)
          CALL hdfwrt3d(varoutr,nxout,nyout,nz,sd_id,0,hdfcompr,        &
              'qvprt','Pert. water vapor specific humidity','kg/kg',    &
                         itmp,hmax,hmin)

        ELSE

!-----------------------------------------------------------------------
!
!  Write out total values to history dump
!
!-----------------------------------------------------------------------

          DO k = 1,nz
            DO j = 1, nyout
              DO i = 1, nxout
                varoutr(i,j,k) = qv(i+ia,j+ja,k)
              END DO
            END DO
          END DO
          CALL edgfill(varoutr,1,nxout,1,innerx, 1,nyout,1,innery, 1,nz,1,nz-1)
          CALL hdfwrt3d(varoutr,nxout,nyout,nz,sd_id,0,hdfcompr,         &
                        'qv','Water vapor specific humidity','kg/kg',    &
                         itmp,hmax,hmin)

        END IF

        DO nq = 1, nscalar
          DO k = 1,nz
            DO j = 1, nyout
              DO i = 1, nxout
                varoutr(i,j,k) = qscalar(i+ia,j+ja,k,nq)
              END DO
            END DO
          END DO
          CALL edgfill(varoutr,1,nxout,1,innerx, 1,nyout,1,innery, 1,nz,1,nz-1)

          ! Turn off bit-packing for Z array if hdfcompr > 3 and mphyopt == 11 - DTD
          hdfcomprtmp = hdfcomprtmp
          IF (nq >= 13 .AND. mphyopt == 11 .AND. hdfcompr > 3) hdfcomprtmp = hdfcompr - 3

          CALL hdfwrt3d(varoutr,nxout,nyout,nz,sd_id,0,hdfcomprtmp,      &
                      'qc','Cloud water mixing ratio','kg/kg',           &
                       itmp,hmax,hmin)
        END DO

        IF(rainout == 1) THEN

          DO j = 1, nyout
            DO i = 1,nxout
              varoutr(i,j,1) = raing(i+ia,j+ja)
            END DO
          END DO
          CALL edgfill(varoutr,1,nxout,1,innerx, 1,nyout,1,innery, 1,1,1,1)
          CALL hdfwrt2d(varoutr,nxout,nyout,sd_id,0,hdfcompr,         &
                      'raing','Grid supersaturation rain','mm',itmp)

          DO j = 1, nyout
            DO i = 1,nxout
              varoutr(i,j,1) = rainc(i+ia,j+ja)
            END DO
          END DO
          CALL edgfill(varoutr,1,nxout,1,innerx, 1,nyout,1,innery, 1,1,1,1)
          CALL hdfwrt2d(varoutr,nxout,nyout,sd_id,0,hdfcompr,         &
                      'rainc','Cumulus convective rain','mm',itmp)

        END IF   !rainout

        IF ( prcout == 1 ) THEN
          DO k = 1,4
            DO j = 1, nyout
              DO i = 1,nxout
                varoutr(i,j,k) = prcrate(i+ia,j+ja,k)
              END DO
            END DO
          END DO
          CALL edgfill(varoutr,1,nxout,1,innerx, 1,nyout,1,innery, 1,4,1,4)

          CALL hdfwrt2d(varoutr(1,1,1),nxout,nyout,sd_id,0,hdfcompr,    &
              'prcrate1','Total precip. rate','kg/(m**2*s)',itmp)
          CALL hdfwrt2d(varoutr(1,1,2),nxout,nyout,sd_id,0,hdfcompr,    &
              'prcrate2','Grid scale precip. rate','kg/(m**2*s)',itmp)
          CALL hdfwrt2d(varoutr(1,1,3),nxout,nyout,sd_id,0,hdfcompr,    &
              'prcrate3','Cumulative precip. rate','kg/(m**2*s)',itmp)
          CALL hdfwrt2d(varoutr(1,1,4),nxout,nyout,sd_id,0,hdfcompr,    &
              'prcrate4','Microphysics precip. rate','kg/(m**2*s)',itmp)
        END IF

      END IF   !mstout

!-----------------------------------------------------------------------
!
!  If tkeout = 1, write out tke.
!
!-----------------------------------------------------------------------

      IF( tkeout == 1 ) THEN

        DO k = 1,nz
          DO j = 1, nyout
            DO i = 1, nxout
              varoutr(i,j,k) = tke(i+ia,j+ja,k)
            END DO
          END DO
        END DO
        CALL edgfill(varoutr,1,nxout,1,innerx, 1,nyout,1,innery, 1,nz,1,nz-1)
        CALL hdfwrt3d(varoutr,nxout,nyout,nz,sd_id,0,hdfcompr,         &
                      'tke','Turbulent Kinetic Energy','(m/s)**2',     &
                       itmp,hmax,hmin)

      END IF

!-----------------------------------------------------------------------
!
!  If trbout = 1, write out the turbulence parameter, km.
!
!-----------------------------------------------------------------------

      IF( trbout == 1 ) THEN

        DO k = 1,nz
          DO j = 1, nyout
            DO i = 1, nxout
              varoutr(i,j,k) = kmh(i+ia,j+ja,k)
            END DO
          END DO
        END DO
        CALL edgfill(varoutr,1,nxout,1,innerx, 1,nyout,1,innery, 1,nz,1,nz-1)
        CALL hdfwrt3d(varoutr,nxout,nyout,nz,sd_id,0,hdfcompr,         &
            'kmh','Hori. turb. mixing coef. for momentum','m**2/s',    &
                       itmp,hmax,hmin)

        DO k = 1,nz
          DO j = 1, nyout
            DO i = 1, nxout
              varoutr(i,j,k) = kmv(i+ia,j+ja,k)
            END DO
          END DO
        END DO
        CALL edgfill(varoutr,1,nxout,1,innerx, 1,nyout,1,innery, 1,nz,1,nz-1)
        CALL hdfwrt3d(varoutr,nxout,nyout,nz,sd_id,0,hdfcompr,         &
            'kmv','Vert. turb. mixing coef. for momentum','m**2/s',    &
                       itmp,hmax,hmin)

      END IF   ! trbout


!-----------------------------------------------------------------------
!
!  If sfcout = 1, write out the surface variables,
!  tsoil, qsoil, and wetcanp.
!
!-----------------------------------------------------------------------

      IF( sfcout == 1) THEN

        DO is=0,nstypout
          DO k = 1,nzsoil
            DO j = 1, nyout
              DO i = 1, nxout
                varout4d(i,j,k,is+1) = tsoil(i+ia,j+ja,k,is)
              END DO
            END DO
          END DO
          CALL edgfill(varout4d(:,:,:,is+1),1,nxout,1,innerx, 1,nyout,1,innery, 1,nzsoil,1,nzsoil)
        END DO
        CALL hdfwrt4d(varout4d,nxout,nyout,nzsoil,nstypout+1,sd_id,0,hdfcompr, &
                    'tsoil','Soil temperature','K',                        &
                     itmpsoil,hmaxsoil,hminsoil)

        DO is=0,nstypout
          DO k = 1,nzsoil
            DO j = 1, nyout
              DO i = 1, nxout
                varout4d(i,j,k,is+1) = qsoil(i+ia,j+ja,k,is)
              END DO
            END DO
          END DO
          CALL edgfill(varout4d(:,:,:,is+1),1,nxout,1,innerx, 1,nyout,1,innery, 1,nzsoil,1,nzsoil)
        END DO
        CALL hdfwrt4d(varout4d,nxout,nyout,nzsoil,nstypout+1,sd_id,0,hdfcompr, &
                      'qsoil','Soil moisture','fraction',                    &
                       itmpsoil,hmaxsoil,hminsoil)

        DO is=0,nstypout
          DO j = 1, nyout
            DO i = 1, nxout
              varoutr(i,j,is+1) = wetcanp(i+ia,j+ja,is)
            END DO
          END DO
        END DO
        CALL edgfill(varoutr(:,:,:),1,nxout,1,innerx, 1,nyout,1,innery, 1,nstypout+1,1,nstypout+1)
        CALL hdfwrt3d(varoutr,nxout,nyout,nstypout+1,sd_id,0,hdfcompr,  &
                    'wetcanp','Canopy water amount','fraction',         &
                     itmp,hmax,hmin)

        IF (snowout == 1) THEN

          DO j = 1, nyout
            DO i = 1, nxout
              varoutr(i,j,1) = snowdpth(i+ia,j+ja)
            END DO
          END DO
          CALL edgfill(varoutr,1,nxout,1,innerx, 1,nyout,1,innery,1,1,1,1)
          CALL hdfwrt2d(varoutr,nxout,nyout,sd_id,0,hdfcompr,           &
                      'snowdpth','Snow depth','m',itmp)
        END IF

      END IF

!-----------------------------------------------------------------------
!
!  If radout = 1, write out the radiation arrays
!
!-----------------------------------------------------------------------

      IF( radout == 1 ) THEN

        DO k = 1,nz
          DO j = 1, nyout
            DO i = 1, nxout
              varoutr(i,j,k) = radfrc(i+ia,j+ja,k)
            END DO
          END DO
        END DO
        CALL edgfill(varoutr,1,nxout,1,innerx, 1,nyout,1,innery, 1,nz,1,nz-1)
        CALL hdfwrt3d(varoutr,nxout,nyout,nz,sd_id,0,hdfcompr,          &
                      'radfrc','Radiation forcing','K/s',               &
                       itmp,hmax,hmin)

        DO j = 1, nyout
          DO i = 1,nxout
            varoutr(i,j,1) = radsw(i+ia,j+ja)
          END DO
        END DO
        CALL edgfill(varoutr,1,nxout,1,innerx, 1,nyout,1,innery, 1,1,1,1)
        CALL hdfwrt2d(varoutr,nxout,nyout,sd_id,0,hdfcompr,             &
            'radsw','Solar radiation reaching the surface','W/m**2',itmp)

        DO j = 1, nyout
          DO i = 1,nxout
            varoutr(i,j,1) = rnflx(i+ia,j+ja)
          END DO
        END DO
        CALL edgfill(varoutr,1,nxout,1,innerx, 1,nyout,1,innery, 1,1,1,1)
        CALL hdfwrt2d(varoutr,nxout,nyout,sd_id,0,hdfcompr,             &
            'rnflx','Net radiation flux absorbed by surface','W/m**2',  &
            itmp)

        DO j = 1, nyout
          DO i = 1,nxout
            varoutr(i,j,1) = radswnet(i+ia,j+ja)
          END DO
        END DO
        CALL edgfill(varoutr,1,nxout,1,innerx, 1,nyout,1,innery, 1,1,1,1)
        CALL hdfwrt2d(varoutr,nxout,nyout,sd_id,0,hdfcompr,             &
              'radswnet','Net solar radiation','W/m**2',itmp)

        DO j = 1, nyout
          DO i = 1,nxout
            varoutr(i,j,1) = radlwin(i+ia,j+ja)
          END DO
        END DO
        CALL edgfill(varoutr,1,nxout,1,innerx, 1,nyout,1,innery, 1,1,1,1)
        CALL hdfwrt2d(varoutr,nxout,nyout,sd_id,0,hdfcompr,             &
              'radlwin','Incoming longwave radiation','W/m**2',         &
              itmp)

      END IF   ! radout

!-----------------------------------------------------------------------
!
!  If flxout = 1, write out the surface fluxes
!
!-----------------------------------------------------------------------

      IF( flxout == 1 ) THEN

        DO j = 1, nyout
          DO i = 1,nxout
            varoutr(i,j,1) = usflx(i+ia,j+ja)
          END DO
        END DO
        CALL edgfill(varoutr,1,nxout,1,innerx, 1,nyout,1,innery, 1,1,1,1)
        CALL hdfwrt2d(varoutr,nxout,nyout,sd_id,0,hdfcompr,         &
            'usflx','Surface flux of u-momentum','kg/(m*s**2)',itmp)

        DO j = 1, nyout
          DO i = 1,nxout
            varoutr(i,j,1) = vsflx(i+ia,j+ja)
          END DO
        END DO
        CALL edgfill(varoutr,1,nxout,1,innerx, 1,nyout,1,innery, 1,1,1,1)
        CALL hdfwrt2d(varoutr,nxout,nyout,sd_id,0,hdfcompr,         &
            'vsflx','Surface flux of v-momentum','kg/(m*s**2)',itmp)

        DO j = 1, nyout
          DO i = 1,nxout
            varoutr(i,j,1) = ptsflx(i+ia,j+ja)
          END DO
        END DO
        CALL edgfill(varoutr,1,nxout,1,innerx, 1,nyout,1,innery, 1,1,1,1)
        CALL hdfwrt2d(varoutr,nxout,nyout,sd_id,0,hdfcompr,         &
            'ptsflx','Surface heat flux','K*kg/(m**2*s)',itmp)

        DO j = 1, nyout
          DO i = 1,nxout
            varoutr(i,j,1) = qvsflx(i+ia,j+ja)
          END DO
        END DO
        CALL edgfill(varoutr,1,nxout,1,innerx, 1,nyout,1,innery, 1,1,1,1)
        CALL hdfwrt2d(varoutr,nxout,nyout,sd_id,0,hdfcompr,          &
            'qvsflx','Surface moisture flux','kg/(m**2*s)',itmp)

      END IF   ! flxout

      600   CONTINUE

      CALL hdfclose(sd_id,istat)
      IF (istat /= 0) THEN
        IF(OUTFLAG) WRITE (6,*) 'HDFSPLITDUMP: ERROR on closing file ',trim(filename),       &
                    ' (status = ',istat,')'
      END IF
    END DO
  END DO

  DEALLOCATE(varoutr, varout4d, varouti)

  IF (hdfcompr > 3) THEN
    DEALLOCATE (itmp,stat=istat)
    DEALLOCATE (hmax,stat=istat)
    DEALLOCATE (hmin,stat=istat)
    DEALLOCATE (itmpsoil,stat=istat)
    DEALLOCATE (hmaxsoil,stat=istat)
    DEALLOCATE (hminsoil,stat=istat)
  END IF

  RETURN
END SUBROUTINE hdfsplitdump

!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE HDFWRT3D                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE hdfwrt3d(var,nx,ny,nz,sd_id,stag_dim,hdfcompr,           &
                    dname,comment,units,itmp,hmax,hmin)

!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Write out a 3-D real array to an HDF4 file.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Gene Bassett
!  2000/03/15
!
!  MODIFICATION HISTORY:
!
!  02/26/2004 Yunheng Wang
!  Added a working array tem to make sure that the data is written as
!  32-bit floating even the default KIND of REAL is not 4. It is
!  allocated in the heap only when needed. Please note that hmax
!  and hmin were redeclared as KIND = 4 floating number in this subroutine.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    var      Array to be written to the file
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!    sd_id    HDF id of the output file
!
!    stag_dim Dimension of grid staggering (0-none, 1-x, 2-y, 3-z)
!
!    hdfcompr Compression flag (0-none)
!
!    name     Variable name
!    comment  Destriptive string
!    units    String destribing units of var
!
!    itmp     Scratch array for mapping reals to integers (used for
!             some values of hdfcompr)
!    hmax     Used to store maximum values as a function of z
!    hmin     Used to store minimum values as a function of z
!
!  OUTPUT:
!
!    None.
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------

  USE arps_precision

  IMPLICIT NONE
  INTEGER :: nx,ny,nz
  REAL(P) :: var(nx,ny,nz)
  INTEGER :: sd_id, stag_dim, hdfcompr
  CHARACTER (LEN=*) :: dname, comment, units
  INTEGER (KIND=INT16) :: itmp(nx,ny,nz)
  REAL(SP)             :: hmax(nz), hmin(nz)

!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------

  INTEGER :: istat
  INTEGER :: dims(3),start(3),stride(3)
  INTEGER :: comp_prm(1)
  INTEGER :: sds_id

  INTEGER :: i,j,k

  REAL    :: amax,amin,scalef
  INTEGER :: itmp1

  REAL(SP), ALLOCATABLE :: tem(:,:,:)
  INTEGER :: istatus

!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------

  INCLUDE 'hdf.f90'    ! HDF4 library include file
  INCLUDE 'mp.inc'

!-----------------------------------------------------------------------
!
!  Functions
!
!-----------------------------------------------------------------------

  INTEGER :: sfcreate, sfscompress, sfscatt, sfsnatt,                   &
             sfwdata, sfendacc

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!-----------------------------------------------------------------------
!
!  Initialize dimension parameters
!
!-----------------------------------------------------------------------

  dims(1) = nx
  dims(2) = ny
  dims(3) = nz
  start(1) = 0
  start(2) = 0
  start(3) = 0
  stride(1) = 1
  stride(2) = 1
  stride(3) = 1

!-----------------------------------------------------------------------
!
!  Create an entry for the variable, set compression and add attributes.
!
!-----------------------------------------------------------------------

  IF(myproc == 0) &
    WRITE (*,*) "HDFWRT3D : Writing variable <",trim(dname),">"

  IF (hdfcompr <= 3) THEN
    sds_id = sfcreate(sd_id, trim(dname), dfnt_float32, 3, dims)
  ELSE
    sds_id = sfcreate(sd_id, trim(dname), dfnt_int16, 3, dims)
  END IF

  comp_prm(1) = 0
  IF (hdfcompr == 1 .OR. hdfcompr == 5) THEN  ! quick gzip
    comp_prm(1) = 1
    istat = sfscompress(sds_id, 4, comp_prm)
    istat = sfsnatt(sds_id, 'hdf_comp_code', dfnt_int32, 1, 4)
  ELSE IF (hdfcompr == 2 .OR. hdfcompr == 6) THEN  ! high gzip
    comp_prm(1) = 6
    istat = sfscompress(sds_id, 4, comp_prm)
    istat = sfsnatt(sds_id, 'hdf_comp_code', dfnt_int32, 1, 4)
  ELSE IF (hdfcompr == 3 .OR. hdfcompr == 7) THEN  ! huffman
    comp_prm(1) = 4   ! this may be a problem on a Cray
    istat = sfscompress(sds_id, 3, comp_prm)
    istat = sfsnatt(sds_id, 'hdf_comp_code', dfnt_int32, 1, 3)
  END IF
  istat = sfsnatt(sds_id, 'hdf_comp_prm', dfnt_int32, 1, comp_prm)

  IF (len_trim(comment) > 0) THEN
    istat = sfscatt(sds_id, 'comment', dfnt_char8,                       &
              len_trim(comment), comment)
  ELSE
    istat = sfscatt(sds_id, 'comment', dfnt_char8, 1, ' ')
  END IF
  IF (len_trim(units) > 0) THEN
    istat = sfscatt(sds_id, 'units', dfnt_char8,                         &
            len_trim(units), units)
  ELSE
    istat = sfscatt(sds_id, 'units', dfnt_char8, 1, ' ')
  END IF
  istat = sfsnatt(sds_id, 'stag_dim', dfnt_int32, 1, stag_dim)

!-----------------------------------------------------------------------
!
!  If called for, map reals to 16 bit integers
!
!-----------------------------------------------------------------------

  IF (hdfcompr > 3) THEN
    DO k=1,nz

!      IF ( mp_opt > 0 .AND. joindmp <= 0 ) THEN
!        CALL a3dmax0   (var(1,1,k),1,nx,1,nx,1,ny,1,ny,1,1,1,1,amax,amin)
!      ELSE
        CALL a3dmax0lcl(var(1,1,k),1,nx,1,nx,1,ny,1,ny,1,1,1,1,amax,amin)
!      ENDIF

      hmax(k) = amax
      hmin(k) = amin

      IF (ABS(hmax(k)) < 1.0E-25) hmax(k)= 0.0  !Added by Xue and Dawson
      IF (ABS(hmin(k)) < 1.0E-25) hmin(k)= 0.0  !Added by "" and ""

      IF (ABS(hmax(k)-hmin(k)) > 1.0E-10) THEN
        scalef = 65534.0 / (hmax(k) - hmin(k))
      ELSE
        scalef = 65534.0
      END IF

      DO j=1,ny
        DO i=1,nx
          itmp1 = nint(scalef * (var(i,j,k) - hmin(k))) - 32767
          itmp(i,j,k) = itmp1
        END DO
      END DO
    END DO
    istat = sfsnatt(sds_id, 'packed16', dfnt_int32, 1, 1)
    istat = sfsnatt(sds_id, 'max', dfnt_float32, nz, hmax)
    istat = sfsnatt(sds_id, 'min', dfnt_float32, nz, hmin)
  END IF

!-----------------------------------------------------------------------
!
!  Write the data in var out to the file
!
!-----------------------------------------------------------------------

  IF (hdfcompr <= 3) THEN
    IF( P == SP) THEN
       istat = sfwdata(sds_id, start, stride, dims, var)
    ELSE
       ALLOCATE(tem(nx,ny,nz), STAT = istatus)
       tem(:,:,:) = var(:,:,:)
       istat = sfwdata(sds_id, start, stride, dims, tem)
       DEALLOCATE(tem)
    END IF
  ELSE
    istat = sfwdata(sds_id, start, stride, dims, itmp)
  END IF

  IF (istat /= 0) THEN
    WRITE (6,'(1x,3a,I3)') "HDFWRT3D : ERROR writing variable <",       &
                           trim(dname),"> with istat = ",istat
  END IF

  istat = sfendacc(sds_id)
  IF (istat /= 0) THEN
    WRITE (6,'(1x,3a,I3)') "HDFWRT3D : ERROR closing dataset for variable <", &
                           trim(dname),"> with istat = ",istat
  END IF

  RETURN
END SUBROUTINE hdfwrt3d

!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE HDFWRT4D                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE hdfwrt4d(var,nx,ny,nz,nl,sd_id,stag_dim,hdfcompr,&
                    dname,comment,units,itmp,hmax,hmin)

!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Write out a 4-D real array to an HDF4 file.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Eric Kemp, 23 May 2002, Based on HDFWRT3D
!
!  MODIFICATION HISTORY:
!
!  07/02/2002 Zuwen He
!  Clean up and make it a general 4d hdfwrt.
!
!  02/26/2004 Yunheng Wang
!  Added a working array tem to make sure that the data is written as
!  32-bit floating even the default KIND of REAL is not 4. It is
!  allocated in the heap only when needed. Please note that hmax
!  and hmin are also used as KIND = 4 floating number in this subroutine.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    var      Array to be written to the file
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!    nl       Number of grid points in the 4th dimension
!
!    sd_id    HDF id of the output file
!
!    stag_dim Dimension of grid staggering (0-none, 1-x, 2-y, 3-z)
!
!    hdfcompr Compression flag (0-none)
!
!    name     Variable name
!    comment  Destriptive string
!    units    String destribing units of var
!
!    itmp     Scratch array for mapping reals to integers (used for
!             some values of hdfcompr)
!    hmax     Used to store maximum values as a function of z
!    hmin     Used to store minimum values as a function of z
!
!  OUTPUT:
!
!    None.
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------

  USE arps_precision

  IMPLICIT NONE
  INTEGER :: nx,ny,nz,nl
  REAL(P) :: var(nx,ny,nz,nl)
  INTEGER :: sd_id, stag_dim, hdfcompr
  CHARACTER (LEN=*) :: dname, comment, units
  INTEGER (KIND=INT16) :: itmp(nx,ny,nz,nl)
  REAL(SP)             :: hmax(nz), hmin(nz)

!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------

  INTEGER :: istat
  INTEGER :: dims(4),start(4),stride(4)
  INTEGER :: comp_prm(1)
  INTEGER :: sds_id

  INTEGER :: i,j,k,l

  REAL    :: scalef
  INTEGER :: itmp1
  REAL    :: amax, amin
  REAL(SP) :: amax4, amin4

  REAL(SP), ALLOCATABLE :: tem(:,:,:,:)
  INTEGER :: istatus

!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------

  INCLUDE 'hdf.f90'    ! HDF4 library include file
  INCLUDE 'mp.inc'

!-----------------------------------------------------------------------
!
!  Functions
!
!-----------------------------------------------------------------------

  INTEGER :: sfcreate, sfscompress, sfscatt, sfsnatt,                   &
             sfwdata, sfendacc

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!-----------------------------------------------------------------------
!
!  Initialize dimension parameters
!
!-----------------------------------------------------------------------

  dims(1) = nx
  dims(2) = ny
  dims(3) = nz
  dims(4) = nl
  start(1) = 0
  start(2) = 0
  start(3) = 0
  start(4) = 0
  stride(1) = 1
  stride(2) = 1
  stride(3) = 1
  stride(4) = 1

!-----------------------------------------------------------------------
!
!  Create an entry for the variable, set compression and add attributes.
!
!-----------------------------------------------------------------------

  IF (myproc == 0) WRITE (*,*) "HDFWRT4D : Writing variable <",trim(dname),">"

  IF (hdfcompr <= 3) THEN
    sds_id = sfcreate(sd_id, trim(dname), dfnt_float32, 4, dims)
  ELSE
    sds_id = sfcreate(sd_id, trim(dname), dfnt_int16, 4, dims)
  END IF

  comp_prm(1) = 0
  IF (hdfcompr == 1 .OR. hdfcompr == 5) THEN  ! quick gzip
    comp_prm(1) = 1
    istat = sfscompress(sds_id, 4, comp_prm)
    istat = sfsnatt(sds_id, 'hdf_comp_code', dfnt_int32, 1, 4)
  ELSE IF (hdfcompr == 2 .OR. hdfcompr == 6) THEN  ! high gzip
    comp_prm(1) = 6
    istat = sfscompress(sds_id, 4, comp_prm)
    istat = sfsnatt(sds_id, 'hdf_comp_code', dfnt_int32, 1, 4)
  ELSE IF (hdfcompr == 3 .OR. hdfcompr == 7) THEN  ! huffman
    comp_prm(1) = 4   ! this may be a problem on a Cray
    istat = sfscompress(sds_id, 3, comp_prm)
    istat = sfsnatt(sds_id, 'hdf_comp_code', dfnt_int32, 1, 3)
  END IF
  istat = sfsnatt(sds_id, 'hdf_comp_prm', dfnt_int32, 1, comp_prm)

  IF (len_trim(comment) > 0) THEN
    istat = sfscatt(sds_id, 'comment', dfnt_char8,                       &
              len_trim(comment), comment)
  ELSE
    istat = sfscatt(sds_id, 'comment', dfnt_char8, 1, ' ')
  END IF
  IF (len_trim(units) > 0) THEN
    istat = sfscatt(sds_id, 'units', dfnt_char8,                         &
            len_trim(units), units)
  ELSE
    istat = sfscatt(sds_id, 'units', dfnt_char8, 1, ' ')
  END IF
  istat = sfsnatt(sds_id, 'stag_dim', dfnt_int32, 1, stag_dim)

!-----------------------------------------------------------------------
!
!  If called for, map reals to 16 bit integers
!
!-----------------------------------------------------------------------

  IF (hdfcompr > 3) THEN
    ! find hmax(nz) and hmin(nz)
    DO l=1,nl
      DO k=1,nz

!        IF ( mp_opt > 0 .AND. joindmp <= 0 ) THEN
!          CALL a3dmax0   (var(1,1,k,l),1,nx,1,nx,1,ny,1,ny,1,1,1,1,amax,amin)
!        ELSE
          CALL a3dmax0lcl(var(1,1,k,l),1,nx,1,nx,1,ny,1,ny,1,1,1,1,amax,amin)
!        ENDIF

        IF ( l == 1 ) THEN   ! initialize
          hmax(k) = amax
          hmin(k) = amin
        END IF

        amax4 = amax
        amin4 = amin
        hmax(k) = MAX(hmax(k), amax4)
        hmin(k) = MIN(hmin(k), amin4)
      END DO
    END DO

    DO l=1,nl
      DO k=1,nz

        if( abs(hmax(k)) < 1.0e-25 ) hmax(k) = 0.0  !Added by Xue and Dawson
        if( abs(hmin(k)) < 1.0e-25 ) hmin(k) = 0.0  !"" "" "" ""

        IF (ABS(hmax(k)-hmin(k)) > 1.0E-10) THEN
          scalef = 65534.0/ (hmax(k) - hmin(k))
        ELSE
          scalef = 65534.0
        END IF

    ! Do map using hmax(nz), hmin(nz)

        DO j=1,ny
          DO i=1,nx
            itmp1 = nint(scalef * (var(i,j,k,l) - hmin(k))) - 32767
            itmp(i,j,k,l) = itmp1
          END DO
        END DO
      END DO
    END DO

    ! Write attributes, hmax(nz) and hmin(nz)
    istat = sfsnatt(sds_id, 'packed16', dfnt_int32, 1, 1)
    istat = sfsnatt(sds_id, 'max', dfnt_float32, nz, hmax)
    istat = sfsnatt(sds_id, 'min', dfnt_float32, nz, hmin)
  END IF

!-----------------------------------------------------------------------
!
!  Write the data in var out to the file
!
!-----------------------------------------------------------------------

  IF (hdfcompr <= 3) THEN
    IF (P == SP) THEN
      istat = sfwdata(sds_id, start, stride, dims, var)
    ELSE
      ALLOCATE(tem(nx,ny,nz,nl), STAT = istatus)
      tem(:,:,:,:) = var(:,:,:,:)
      istat = sfwdata(sds_id, start, stride, dims, tem)
      DEALLOCATE(tem)
    END IF
  ELSE
    istat = sfwdata(sds_id, start, stride, dims, itmp)
  END IF

  IF (istat /= 0) THEN
    WRITE (6,*) "HDFWRT4D : ERROR writing variable <",trim(dname),">"
  END IF

  istat = sfendacc(sds_id)
  IF (istat /= 0) THEN
    WRITE (6,*) "HDFWRT4D : ERROR writing variable <",trim(dname),">"
  END IF

  RETURN
END SUBROUTINE hdfwrt4d

!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE HDFWRT3DI                  ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE hdfwrt3di(var,nx,ny,nz,sd_id,stag_dim,hdfcompr,              &
           name,comment,units)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Write out a 3-D integer array to an HDF4 file.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Gene Bassett
!  2000/03/15
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    var      Array to be written to the file
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!    sd_id    HDF id of the output file
!
!    stag_dim Dimension of grid staggering (0-none, 1-x, 2-y, 3-z)
!
!    hdfcompr Compression flag (0-none)
!
!    name     Variable name
!    comment  Destriptive string
!    units    String destribing units of var
!
!  OUTPUT:
!
!    None.
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------

  IMPLICIT NONE
  INTEGER :: nx,ny,nz
  INTEGER :: var(nx,ny,nz)
  INTEGER :: sd_id, stag_dim, hdfcompr
  CHARACTER (LEN=*) :: name, comment, units

!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------

  INTEGER :: istat
  INTEGER :: dims(3),start(3),stride(3)
  INTEGER :: sds_id
  INTEGER :: comp_prm(1)

!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------

  INCLUDE 'hdf.f90'    ! HDF4 library include file
  INCLUDE 'mp.inc'

!-----------------------------------------------------------------------
!
!  Functions
!
!-----------------------------------------------------------------------

  INTEGER :: sfcreate, sfscompress, sfscatt, sfsnatt,                   &
          sfwdata, sfendacc

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!-----------------------------------------------------------------------
!
!  Initialize dimension parameters
!
!-----------------------------------------------------------------------

  dims(1) = nx
  dims(2) = ny
  dims(3) = nz
  start(1) = 0
  start(2) = 0
  start(3) = 0
  stride(1) = 1
  stride(2) = 1
  stride(3) = 1

!-----------------------------------------------------------------------
!
!  Create an entry for the variable, set compression and add attributes.
!
!-----------------------------------------------------------------------

  IF (myproc == 0) &
  WRITE (*,*) "HDFWRT3DI: Writing variable <",trim(name),">"
  sds_id = sfcreate(sd_id, name, dfnt_int32, 3, dims)

  comp_prm(1) = 0
  IF (hdfcompr == 1 .OR. hdfcompr == 5) THEN  ! quick gzip
    comp_prm(1) = 1
    istat = sfscompress(sds_id, 4, comp_prm)
    istat = sfsnatt(sds_id, 'hdf_comp_code', dfnt_int32, 1, 4)
  ELSE IF (hdfcompr == 2 .OR. hdfcompr == 6) THEN  ! high gzip
    comp_prm(1) = 6
    istat = sfscompress(sds_id, 4, comp_prm)
    istat = sfsnatt(sds_id, 'hdf_comp_code', dfnt_int32, 1, 4)
  ELSE IF (hdfcompr == 3 .OR. hdfcompr == 7) THEN  ! huffman
    comp_prm(1) = 4   ! this may be a problem on a Cray
    istat = sfscompress(sds_id, 3, comp_prm)
    istat = sfsnatt(sds_id, 'hdf_comp_code', dfnt_int32, 1, 3)
  END IF
  istat = sfsnatt(sds_id, 'hdf_comp_prm', dfnt_int32, 1, comp_prm)

  IF (len_trim(comment) > 0) THEN
    istat = sfscatt(sds_id, 'comment', dfnt_char8,                       &
              len_trim(comment), comment)
  ELSE
    istat = sfscatt(sds_id, 'comment', dfnt_char8, 1, ' ')
  END IF
  IF (len_trim(units) > 0) THEN
    istat = sfscatt(sds_id, 'units', dfnt_char8,                         &
            len_trim(units), units)
  ELSE
    istat = sfscatt(sds_id, 'units', dfnt_char8, 1, ' ')
  END IF
  istat = sfsnatt(sds_id, 'stag_dim', dfnt_int32, 1, stag_dim)

!-----------------------------------------------------------------------
!
!  Write the data in var out to the file
!
!-----------------------------------------------------------------------

  istat = sfwdata(sds_id, start, stride, dims, var)

  IF (istat /= 0) THEN
    IF(myproc == 0) &
    WRITE (6,*) "HDFWRT3DI: ERROR writing variable <",trim(name),">"
  END IF
  istat = sfendacc(sds_id)
  IF (istat /= 0) THEN
    IF(myproc == 0) &
    WRITE (6,*) "HDFWRT3DI: ERROR writing variable <",trim(name),">"
  END IF

  RETURN
END SUBROUTINE hdfwrt3di

!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE HDFWRT2D                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE hdfwrt2d(var,nx,ny,sd_id,stag_dim,hdfcompr,                  &
                    dname,comment,units,itmp)

!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Write out a 2-D real array to an HDF4 file.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Gene Bassett
!  2000/03/15
!
!  MODIFICATION HISTORY:
!
!  02/26/2004 Yunheng Wang
!  Added a working array tem to make sure that the data is written as
!  32-bit floating even the default KIND of REAL is not 4. It is
!  allocated in the heap only when needed. Please note that amax
!  and amin were declared as KIND = 4 floating number.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    var      Array to be written to the file
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!
!    sd_id    HDF id of the output file
!
!    stag_dim Dimension of grid staggering (0-none, 1-x, 2-y, 3-z)
!
!    hdfcompr Compression flag (0-none)
!
!    name     Variable name
!    comment  Destriptive string
!    units    String destribing units of var
!
!    itmp     Scratch array for mapping reals to integers (used for
!             some values of hdfcompr)
!
!  OUTPUT:
!
!    None.
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------

  USE arps_precision

  IMPLICIT NONE
  INTEGER :: nx,ny
  REAL(P) :: var(nx,ny)
  INTEGER :: sd_id, stag_dim, hdfcompr
  CHARACTER (LEN=*) :: dname, comment, units
  INTEGER (KIND=INT16) :: itmp(nx,ny)

!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------

  INTEGER :: istat
  INTEGER :: dims(2),start(2),stride(2)
  INTEGER :: sds_id
  INTEGER :: comp_prm(1)
  REAL    :: amax,amin,scalef
  REAL(SP) :: amax4,amin4
  INTEGER :: i,j
  INTEGER :: itmp1

  REAL(SP), ALLOCATABLE :: tem(:,:)
  INTEGER :: istatus

!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------

  INCLUDE 'hdf.f90'    ! HDF4 library include file
  INCLUDE 'mp.inc'

!-----------------------------------------------------------------------
!
!  Functions
!
!-----------------------------------------------------------------------

  INTEGER :: sfcreate, sfscompress, sfscatt, sfsnatt,                   &
             sfwdata, sfendacc

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  IF (myproc == 0) &
    WRITE (*,*) "HDFWRT2D : Writing variable <",trim(dname),">"

!-----------------------------------------------------------------------
!
!  Initialize dimension parameters
!
!-----------------------------------------------------------------------

  dims(1) = nx
  dims(2) = ny
  start(1) = 0
  start(2) = 0
  stride(1) = 1
  stride(2) = 1

!-----------------------------------------------------------------------
!
!  Create an entry for the variable, set compression and add attributes.
!
!-----------------------------------------------------------------------

  IF (hdfcompr <= 3) THEN
    sds_id = sfcreate(sd_id, trim(dname), dfnt_float32, 2, dims)
  ELSE
    sds_id = sfcreate(sd_id, trim(dname), dfnt_int16, 2, dims)
  END IF

  comp_prm(1) = 0
  IF (hdfcompr == 1 .OR. hdfcompr == 5) THEN       ! quick gzip
    comp_prm(1) = 1
    istat = sfscompress(sds_id, 4, comp_prm)
    istat = sfsnatt(sds_id, 'hdf_comp_code', dfnt_int32, 1, 4)
  ELSE IF (hdfcompr == 2 .OR. hdfcompr == 6) THEN  ! high gzip
    comp_prm(1) = 6
    istat = sfscompress(sds_id, 4, comp_prm)
    istat = sfsnatt(sds_id, 'hdf_comp_code', dfnt_int32, 1, 4)
  ELSE IF (hdfcompr == 3 .OR. hdfcompr == 7) THEN  ! huffman
    comp_prm(1) = 4                    ! this may be a problem on a Cray
    istat = sfscompress(sds_id, 3, comp_prm)
    istat = sfsnatt(sds_id, 'hdf_comp_code', dfnt_int32, 1, 3)
  END IF
  istat = sfsnatt(sds_id, 'hdf_comp_prm', dfnt_int32, 1, comp_prm)

  IF (len_trim(comment) > 0) THEN
    istat = sfscatt(sds_id, 'comment', dfnt_char8,                       &
              len_trim(comment), comment)
  ELSE
    istat = sfscatt(sds_id, 'comment', dfnt_char8, 1, ' ')
  END IF
  IF (len_trim(units) > 0) THEN
    istat = sfscatt(sds_id, 'units', dfnt_char8,                         &
            len_trim(units), units)
  ELSE
    istat = sfscatt(sds_id, 'units', dfnt_char8, 1, ' ')
  END IF
  istat = sfsnatt(sds_id, 'stag_dim', dfnt_int32, 1, stag_dim)

!-----------------------------------------------------------------------
!
!  If called for, map reals to 16 bit integers
!
!-----------------------------------------------------------------------

  IF (hdfcompr > 3) THEN

!    IF ( mp_opt > 0 .AND. joindmp <= 0 ) THEN
!      CALL a3dmax0   (var,1,nx,1,nx,1,ny,1,ny,1,1,1,1,amax,amin)
!    ELSE
      CALL a3dmax0lcl(var,1,nx,1,nx,1,ny,1,ny,1,1,1,1,amax,amin)
!    ENDIF

    if( abs(amax) < 1.0e-25 ) amax = 0.0
    if( abs(amin) < 1.0e-25 ) amin = 0.0

    IF (ABS(amax-amin) > 1.0E-10) THEN
      scalef = 65534.0 / (amax - amin)
    ELSE
      scalef = 65534.0
    END IF

    DO j=1,ny
      DO i=1,nx
        itmp1 = nint(scalef * (var(i,j) - amin)) - 32767
        itmp(i,j) = itmp1
      END DO
    END DO
    amax4 = amax
    amin4 = amin
    istat = sfsnatt(sds_id, 'packed16', dfnt_int32, 1, 1)
    istat = sfsnatt(sds_id, 'max', dfnt_float32, 1, amax4)
    istat = sfsnatt(sds_id, 'min', dfnt_float32, 1, amin4)
  END IF

!-----------------------------------------------------------------------
!
!  Write the data in var out to the file
!
!-----------------------------------------------------------------------

  IF (hdfcompr <= 3) THEN
    IF (P == SP) THEN
      istat = sfwdata(sds_id, start, stride, dims, var)
    ELSE
      ALLOCATE(tem(nx,ny), STAT = istatus)
      tem(:,:) = var(:,:)
      istat = sfwdata(sds_id, start, stride, dims, tem)
      DEALLOCATE(tem)
    END IF
  ELSE
    istat = sfwdata(sds_id, start, stride, dims, itmp)
  END IF

  IF (istat /= 0) THEN
    WRITE (6,*) "HDFWRT2D : ERROR writing variable <",trim(dname),">"
  END IF
  istat = sfendacc(sds_id)
  IF (istat /= 0) THEN
    WRITE (6,*) "HDFWRT2D : ERROR writing variable <",trim(dname),">"
  END IF

  RETURN
END SUBROUTINE hdfwrt2d

!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE HDFWRT2DI                  ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################

SUBROUTINE hdfwrt2di(var,nx,ny,sd_id,stag_dim,hdfcompr,                 &
           dname,comment,units)

!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Write out a 2-D integer array to an HDF4 file.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Gene Bassett
!  2000/03/15
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    var      Array to be written to the file
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!
!    sd_id    HDF id of the output file
!
!    stag_dim Dimension of grid staggering (0-none, 1-x, 2-y, 3-z)
!
!    hdfcompr Compression flag (0-none)
!
!    name     Variable name
!    comment  Destriptive string
!    units    String destribing units of var
!
!  OUTPUT:
!
!    None.
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
  INTEGER :: nx,ny
  INTEGER :: var(nx,ny)
  INTEGER :: sd_id, stag_dim, hdfcompr
  CHARACTER (LEN=*) :: dname, comment, units

!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------

  INTEGER :: istat
  INTEGER :: dims(2),start(2),stride(2)
  INTEGER :: sds_id
  INTEGER :: comp_prm(1)

!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------

  INCLUDE 'hdf.f90'    ! HDF4 library include file
  INCLUDE 'mp.inc'

!-----------------------------------------------------------------------
!
!  Functions
!
!-----------------------------------------------------------------------

  INTEGER :: sfcreate, sfscompress, sfscatt, sfsnatt,                   &
             sfwdata, sfendacc

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!-----------------------------------------------------------------------
!
!  Initialize dimension parameters
!
!-----------------------------------------------------------------------

  dims(1)   = nx
  dims(2)   = ny
  start(1)  = 0
  start(2)  = 0
  stride(1) = 1
  stride(2) = 1

!-----------------------------------------------------------------------
!
!  Create an entry for the variable, set compression and add attributes.
!
!-----------------------------------------------------------------------

  IF (myproc == 0) &
  WRITE (*,*) "HDFWRT2DI: Writing variable <",trim(dname),">"

  sds_id = sfcreate(sd_id, trim(dname), dfnt_int32, 2, dims)

  comp_prm(1) = 0
  IF (hdfcompr == 1 .OR. hdfcompr == 5) THEN  ! quick gzip
    comp_prm(1) = 1
    istat = sfscompress(sds_id, 4, comp_prm)
    istat = sfsnatt(sds_id, 'hdf_comp_code', dfnt_int32, 1, 4)
  ELSE IF (hdfcompr == 2 .OR. hdfcompr == 6) THEN  ! high gzip
    comp_prm(1) = 6
    istat = sfscompress(sds_id, 4, comp_prm)
    istat = sfsnatt(sds_id, 'hdf_comp_code', dfnt_int32, 1, 4)
  ELSE IF (hdfcompr == 3 .OR. hdfcompr == 7) THEN  ! huffman
    comp_prm(1) = 4   ! this may be a problem on a Cray
    istat = sfscompress(sds_id, 3, comp_prm)
    istat = sfsnatt(sds_id, 'hdf_comp_code', dfnt_int32, 1, 3)
  END IF
  istat = sfsnatt(sds_id, 'hdf_comp_prm', dfnt_int32, 1, comp_prm)

  IF (len_trim(comment) > 0) THEN
    istat = sfscatt(sds_id, 'comment', dfnt_char8,                       &
              len_trim(comment), comment)
  ELSE
    istat = sfscatt(sds_id, 'comment', dfnt_char8, 1, ' ')
  END IF
  IF (len_trim(units) > 0) THEN
    istat = sfscatt(sds_id, 'units', dfnt_char8,                         &
            len_trim(units), units)
  ELSE
    istat = sfscatt(sds_id, 'units', dfnt_char8, 1, ' ')
  END IF
  istat = sfsnatt(sds_id, 'stag_dim', dfnt_int32, 1, stag_dim)

!-----------------------------------------------------------------------
!
!  Write the data in var out to the file
!
!-----------------------------------------------------------------------

  istat = sfwdata(sds_id, start, stride, dims, var)

  IF (istat /= 0) THEN
    IF(myproc == 0) &
    WRITE (6,*) "HDFWRT2DI: ERROR writing variable <",trim(dname),">"
  END IF
  istat = sfendacc(sds_id)
  IF (istat /= 0) THEN
    IF(myproc == 0) &
    WRITE (6,*) "HDFWRT2DI: ERROR writing variable <",trim(dname),">"
  END IF

  RETURN
END SUBROUTINE hdfwrt2di

!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE HDFWRT1D                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################

SUBROUTINE hdfwrt1d(var,num,sd_id,dname,comment,units)

!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Write out a 1-D real array to an HDF4 file.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Gene Bassett
!  2000/03/15
!
!  MODIFICATION HISTORY:
!
!  02/26/2004 Yunheng Wang
!  Added a working array tem to make sure that the data is written as
!  32-bit floating even the default KIND of REAL is not 4. It is
!  allocated in the heap only when needed.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    var      Array to be written to the file
!
!    num      Number of grid points
!
!    sd_id    HDF id of the output file
!
!    name     Variable name
!    comment  Destriptive string
!    units    String destribing units of var
!
!  OUTPUT:
!
!    None.
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------

  USE arps_precision

  IMPLICIT NONE
  INTEGER :: num
  REAL(P) :: var(num)
  INTEGER :: sd_id
  CHARACTER (LEN=*) :: dname, comment, units

!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------

  INTEGER :: istat
  INTEGER :: dims(1),start(1),stride(1)
  INTEGER :: sds_id

  REAL(SP), ALLOCATABLE :: tem(:)
  INTEGER :: istatus

!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------

  INCLUDE 'hdf.f90'    ! HDF4 library include file
  INCLUDE 'mp.inc'

!-----------------------------------------------------------------------
!
!  Functions
!
!-----------------------------------------------------------------------

  INTEGER :: sfcreate, sfscatt, sfwdata, sfendacc

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!-----------------------------------------------------------------------
!
!  Initialize dimension parameters
!
!-----------------------------------------------------------------------

  dims(1)   = num
  start(1)  = 0
  stride(1) = 1

!-----------------------------------------------------------------------
!
!  Create an entry for the variable, set compression and add attributes.
!
!-----------------------------------------------------------------------

  IF (myproc == 0) &
  WRITE (*,*) "HDFWRT1D : Writing variable <",trim(dname),">"

  sds_id = sfcreate(sd_id, trim(dname), dfnt_float32, 1, dims)

  IF (len_trim(comment) > 0) THEN
    istat = sfscatt(sds_id, 'comment', dfnt_char8,                       &
              len_trim(comment), comment)
  ELSE
    istat = sfscatt(sds_id, 'comment', dfnt_char8, 1, ' ')
  END IF
  IF (len_trim(units) > 0) THEN
    istat = sfscatt(sds_id, 'units', dfnt_char8,                         &
            len_trim(units), units)
  ELSE
    istat = sfscatt(sds_id, 'units', dfnt_char8, 1, ' ')
  END IF

!-----------------------------------------------------------------------
!
!  Write the data in var out to the file
!
!-----------------------------------------------------------------------

  IF( P == SP) THEN
    istat = sfwdata(sds_id, start, stride, dims, var)
  ELSE
    ALLOCATE(tem(num), STAT = istatus)
    tem(:) = var(:)
    istat = sfwdata(sds_id, start, stride, dims, tem)
    DEALLOCATE(tem)
  END IF

  IF (istat /= 0) THEN
    WRITE (6,*) "HDFWRT1D : ERROR writing variable <",trim(dname),">"
  END IF
  istat = sfendacc(sds_id)
  IF (istat /= 0) THEN
    WRITE (6,*) "HDFWRT1D : ERROR writing variable <",trim(dname),">"
  END IF

  RETURN
END SUBROUTINE hdfwrt1d

!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE HDFWRT1DI                  ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################

SUBROUTINE hdfwrt1di(var,num,sd_id,name,comment,units)

!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Write out a 1-D integer array to an HDF4 file.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster, based on hdfwrt1d
!  10/30/2002
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    var      Integer array to be written to the file
!
!    num      Number of grid points
!
!    sd_id    HDF id of the output file
!
!    name     Variable name
!    comment  Destriptive string
!    units    String destribing units of var
!
!  OUTPUT:
!
!    None.
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------

  IMPLICIT NONE
  INTEGER :: num
  INTEGER :: var(num)
  INTEGER :: sd_id
  CHARACTER (LEN=*) :: name, comment, units

!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------

  INTEGER :: istat
  INTEGER :: dims(1),start(1),stride(1)
  INTEGER :: sds_id

!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------

  INCLUDE 'hdf.f90'    ! HDF4 library include file
  INCLUDE 'mp.inc'

!-----------------------------------------------------------------------
!
!  Functions
!
!-----------------------------------------------------------------------

  INTEGER :: sfcreate, sfscatt, sfwdata, sfendacc

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!-----------------------------------------------------------------------
!
!  Initialize dimension parameters
!
!-----------------------------------------------------------------------

  dims(1) = num
  start(1) = 0
  stride(1) = 1

!-----------------------------------------------------------------------
!
!  Create an entry for the variable and add attributes.
!  No compression is done for 1D
!!
!-----------------------------------------------------------------------

  IF (myproc == 0) &
  WRITE (*,*) "HDFWRT1DI: Writing variable <",trim(name),">"

  sds_id = sfcreate(sd_id, trim(name), dfnt_int32, 1, dims)

  IF (len_trim(comment) > 0) THEN
    istat = sfscatt(sds_id, 'comment', dfnt_char8,                       &
              len_trim(comment), comment)
  ELSE
    istat = sfscatt(sds_id, 'comment', dfnt_char8, 1, ' ')
  END IF
  IF (len_trim(units) > 0) THEN
    istat = sfscatt(sds_id, 'units', dfnt_char8,                         &
            len_trim(units), units)
  ELSE
    istat = sfscatt(sds_id, 'units', dfnt_char8, 1, ' ')
  END IF

!-----------------------------------------------------------------------
!
!  Write the data in var out to the file
!
!-----------------------------------------------------------------------

  istat = sfwdata(sds_id, start, stride, dims, var)

  IF (istat /= 0) THEN
    IF(myproc == 0) &
    WRITE (6,*) "HDFWRT1DI: ERROR writing variable <",trim(name),">"
  END IF
  istat = sfendacc(sds_id)
  IF (istat /= 0) THEN
    WRITE (6,*) "HDFWRT1DI: ERROR writing variable <",trim(name),">"
  END IF

  RETURN
END SUBROUTINE hdfwrt1di

!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE HDFWRTR                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################

SUBROUTINE hdfwrtr(sd_id,name,val,istat)

!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Write out a real attribute
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Gene Bassett
!  2000/04/05
!
!  MODIFICATION HISTORY:
!
!  02/26/2004 Yunheng Wang
!  Added a temporary variable tem to make sure that the data is written
!  as 32-bit floating even the default KIND of REAL is not 4.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    name     Variable name
!
!    sd_id    HDF id of the file or variable containing the
!             named attribute
!
!  OUTPUT:
!
!    val      The value of the attribute
!    istat     Status of the read (0-okay, 1-write error)
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------

  IMPLICIT NONE
  REAL              :: val
  INTEGER           :: sd_id
  CHARACTER (LEN=*) :: name
  INTEGER           :: istat

!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------

  REAL(4) :: tem

!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------

  INCLUDE 'hdf.f90'    ! HDF4 library include file
  INCLUDE 'mp.inc'     ! mpi include file
!-----------------------------------------------------------------------
!
!  Functions
!
!-----------------------------------------------------------------------

  INTEGER :: sfsnatt

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  tem = val
  istat = sfsnatt(sd_id, trim(name), dfnt_float32, 1, tem)

  IF (istat == -1) THEN
    WRITE (6,*) "HDFWRTR: ERROR writing variable ",trim(name),"."
    istat = 1
  ELSE
!  IF(myproc == 0)  &
!    WRITE (*,*) "HDFWRTR: Wrote variable ",trim(name)," value:",val
  END IF

  RETURN
END SUBROUTINE hdfwrtr

!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE HDFWRTI                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################

SUBROUTINE hdfwrti(sd_id,name,val,istat)

!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Write out an integer attribute
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Gene Bassett
!  2000/04/05
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    name     Variable name
!
!    sd_id    HDF id of the file or variable containing the
!             named attribute
!
!  OUTPUT:
!
!    val      The value of the attribute
!    istat     Status of the read (0-okay, 1-write error)
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------

  IMPLICIT NONE
  INTEGER           :: val
  INTEGER           :: sd_id
  CHARACTER (LEN=*) :: name
  INTEGER           :: istat

!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------

  INCLUDE 'hdf.f90'    ! HDF4 library include file
  INCLUDE 'mp.inc'     ! mpi parameters

!-----------------------------------------------------------------------
!
!  Functions
!
!-----------------------------------------------------------------------

  INTEGER :: sfsnatt

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  istat = sfsnatt(sd_id, trim(name), dfnt_int32, 1, val)

  IF (istat == -1) THEN
    IF (myproc == 0) &
    WRITE (6,*) "HDFWRTI: ERROR writing variable ",trim(name),"."
    istat = 1
  ELSE
    !IF (myproc == 0) &
    !WRITE (*,*) "HDFWRTI: Wrote variable ",trim(name)," value:",val
  END IF

  RETURN
END SUBROUTINE hdfwrti

!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE HDFWRTC                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################

SUBROUTINE hdfwrtc(sd_id,strlen,name,string,istat)

!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Write out a string attribute
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Gene Bassett
!  2000/04/11
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    name     Variable name
!
!    strlen   Length of string to be written out (set to len(string)
!             if strlen passed in as 0)
!
!    sd_id    HDF id of the file or variable containing the
!             named attribute
!
!  OUTPUT:
!
!    string   The string to be written out
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------

  IMPLICIT NONE
  INTEGER :: sd_id
  CHARACTER (LEN=*) :: name
  CHARACTER (LEN=*) :: string
  INTEGER :: strlen, istat

!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------

 INTEGER :: tmplen     ! This replaces the use of strlen which is passed
                       ! in as an expression sometimes and a "expression
                       ! is changed by subprogram" warning will arise in
                       ! such a case.   -- WYH.

!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------

  INCLUDE 'hdf.f90'    ! HDF4 library include file
  INCLUDE 'mp.inc'     ! mpi parameters

!-----------------------------------------------------------------------
!
!  Functions
!
!-----------------------------------------------------------------------

  INTEGER :: sfscatt

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


  IF (strlen == 0) THEN
    tmplen = LEN_TRIM(string)
  ELSE
    tmplen = strlen
  END IF

!  IF (strlen ==0) strlen = LEN(string)

  istat = sfscatt(sd_id, trim(name), dfnt_char8, tmplen, string)

  IF (istat == -1) THEN
    IF (myproc == 0) &
    WRITE (6,*) "HDFWRTC: ERROR writing variable ",trim(name),"."
    istat = 1
  ELSE
    !IF (myproc == 0) &
    !WRITE (*,*) "HDFWRTC: Wrote variable ",trim(name)," value:", trim(string)
  END IF

  RETURN
END SUBROUTINE hdfwrtc

!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE HDFRD3D                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################

SUBROUTINE hdfrd3d(sd_id,dname,nx,ny,nz,var,istat,itmp,hmax,hmin)

!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Read in a 3-D real array from an HDF4 file.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Gene Bassett
!  2000/03/15
!
!  MODIFICATION HISTORY:
!
!  02/26/2004 Yunheng Wang
!  Added a working array tem to make sure that the data is read in
!  corectly even the default KIND of REAL is not 4. It is
!  allocated in the heap only when needed. Please note that hmax
!  and hmin were redeclared as KIND = 4 floating number.
!
!  06/02/2011 Yunheng Wang
!  Use module arps_precision and always assum var is single precision.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    name     Variable name
!
!    sd_id    HDF id of the output file
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!    itmp     Scratch array for mapping reals to integers (used for
!             some values of hdfcompr)
!    hmax     Used to store maximum values as a function of z
!    hmin     Used to store minimum values as a function of z
!
!  OUTPUT:
!
!    var      Array to be read in
!
!    istat     Status of read (0-okay, 1-error when reading)
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------

  USE arps_precision

  IMPLICIT NONE
  INTEGER           :: nx,ny,nz
  REAL(P)           :: var(nx,ny,nz)
  INTEGER           :: sd_id
  CHARACTER (LEN=*) :: dname
  INTEGER           :: istat
  INTEGER(INT16)    :: itmp(nx,ny,nz)
  REAL(SP)          :: hmax(nz),hmin(nz)

!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------

  INTEGER :: dims(3),start(3),stride(3)
  INTEGER :: sds_index,sds_id,attr_index
  INTEGER :: packed16
  INTEGER :: istat1,istat2,istat3
  REAL    :: scalef
  INTEGER :: i,j,k

  INTEGER :: istatus

  REAL(SP), ALLOCATABLE :: dtain(:,:,:)

!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------

  INCLUDE 'hdf.f90'    ! HDF4 library include file
  INCLUDE 'mp.inc'

!-----------------------------------------------------------------------
!
!  Functions
!
!-----------------------------------------------------------------------

  INTEGER :: sfn2index, sfselect, sfrdata, sffattr, sfrnatt, sfendacc

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!-----------------------------------------------------------------------
!
!  Initialize dimension parameters
!
!-----------------------------------------------------------------------

  dims(1) = nx
  dims(2) = ny
  dims(3) = nz
  start(1) = 0
  start(2) = 0
  start(3) = 0
  stride(1) = 1
  stride(2) = 1
  stride(3) = 1

!-----------------------------------------------------------------------
!
!  Get the SDS ID for the variable.
!
!-----------------------------------------------------------------------

  sds_index = sfn2index(sd_id, trim(dname))

  IF (sds_index == -1) THEN
    WRITE (6,*) "HDFRD3D: WARNING: variable ",                            &
                 trim(dname)," not found in file."
    istat = 1
    RETURN
  END IF

  sds_id = sfselect(sd_id, sds_index)

  attr_index = sffattr(sds_id, "packed16")
  IF (attr_index >= 0) THEN
    istat1 = sfrnatt(sds_id, attr_index, packed16)
    attr_index = sffattr(sds_id, "max")
    istat2 = sfrnatt(sds_id, attr_index, hmax)
    attr_index = sffattr(sds_id, "min")
    istat3 = sfrnatt(sds_id, attr_index, hmin)
    IF (istat1 == -1 .OR. istat2 == -1 .OR. istat3 == -1) THEN
      WRITE (6,*) "HDFRD3D: ERROR reading max/min for ",trim(dname)
      istat = 2
      RETURN
    END IF
  ELSE
    packed16 = 0
  END IF

!-----------------------------------------------------------------------
!
!  Read data into var.
!
!-----------------------------------------------------------------------

  !IF (myproc == 0) WRITE (*,*) "HDFRD3D: Reading variable ",trim(dname)
  IF (packed16 > 0) THEN
    istat = sfrdata(sds_id, start, stride, dims, itmp)
  ELSE
    IF (P == SP) THEN
      istat = sfrdata(sds_id, start, stride, dims, var)
    ELSE
      ALLOCATE(dtain(nx,ny,nz), STAT = istatus)
      istat = sfrdata(sds_id, start, stride, dims, dtain)
      var(:,:,:) = dtain(:,:,:)
      DEALLOCATE(dtain)
    END IF
  END IF

  IF (istat == -1) THEN
    WRITE (6,*) "HDFRD3D: ERROR reading variable ",trim(dname),"."
    istat = 2
    RETURN
  END IF

!-----------------------------------------------------------------------
!
!  If called for, map 16 bit integers to reals
!
!-----------------------------------------------------------------------

  IF (packed16 /= 0) THEN
    DO k=1,nz
      scalef = (hmax(k)-hmin(k)) / 65534.0
      DO j=1,ny
        DO i=1,nx
          var(i,j,k) = scalef * (itmp(i,j,k) + 32767) + hmin(k)
        END DO
      END DO
    END DO
  END IF

  istat = sfendacc(sds_id)
  IF (istat /= 0) THEN
    WRITE (6,*) "HDFRD3D: ERROR reading variable ",trim(dname)
    istat = 2
  END IF

  RETURN
END SUBROUTINE hdfrd3d

!##################################################################
!##################################################################
!######                                                      ######
!######              SUBROUTINE HDFRD4D                      ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################

SUBROUTINE hdfrd4d(sd_id,dname,nx,ny,nz,nl,var,istat,        &
                   itmp,hmax,hmin)

!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Read in a 3-D real array from an HDF4 file.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: J. Brotzge (modified from hdfrd3d)
!  05/23/2002
!
!  MODIFICATION HISTORY:
!
!  07/02/2002 Zuwen He
!  Clean up and make it a general 4d HDF writer
!
!  02/26/2004 Yunheng Wang
!  Added a working array tem to make sure that the data is read in
!  corectly even the default KIND of REAL is not 4. It is
!  allocated in the heap only when needed. Please note that hmax
!  and hmin were redeclared as KIND = 4 floating number.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    name     Variable name
!
!    sd_id    HDF id of the output file
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!    nl       Number of grid points in the 4th dimension
!
!    itmp     Scratch array for mapping reals to integers (used for
!             some values of hdfcompr)
!    hmax     Used to store maximum values as a function of z
!    hmin     Used to store minimum values as a function of z
!
!  OUTPUT:
!
!    var      Array to be read in
!
!    istat     Status of read (0-okay, 1-error when reading)
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------

  USE arps_precision
  IMPLICIT NONE
  INTEGER , INTENT(IN)  :: nx,ny,nz,nl
  INTEGER , INTENT(IN)  :: sd_id
  CHARACTER (LEN=*), INTENT(IN) :: dname
  REAL(P),  INTENT(OUT) :: var(nx,ny,nz,nl)
  INTEGER,  INTENT(OUT) :: istat
  INTEGER(INT16)        :: itmp(nx,ny,nz,nl)
  REAL(SP)              :: hmax(nz),hmin(nz)

!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------

  INTEGER :: dims(4),start(4),stride(4)
  INTEGER :: sds_index,sds_id,attr_index
  INTEGER :: packed16
  INTEGER :: istat1,istat2,istat3
  REAL(P) :: scalef
  INTEGER :: i,j,k,l

  INTEGER :: istatus

  REAL(SP), ALLOCATABLE :: dtain(:,:,:,:)

!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------

  INCLUDE 'hdf.f90'    ! HDF4 library include file
  INCLUDE 'mp.inc'

!-----------------------------------------------------------------------
!
!  Functions
!
!-----------------------------------------------------------------------

  INTEGER :: sfn2index, sfselect, sfrdata, sffattr, sfrnatt, sfendacc

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!-----------------------------------------------------------------------
!
!  Initialize dimension parameters
!
!-----------------------------------------------------------------------

  dims(1) = nx
  dims(2) = ny
  dims(3) = nz
  dims(4) = nl
  start(1) = 0
  start(2) = 0
  start(3) = 0
  start(4) = 0
  stride(1) = 1
  stride(2) = 1
  stride(3) = 1
  stride(4) = 1

!-----------------------------------------------------------------------
!
!  Get the SDS ID for the variable.
!
!-----------------------------------------------------------------------

  sds_index = sfn2index(sd_id, trim(dname))

  IF (sds_index == -1) THEN
    WRITE (6,*) "HDFRD4D: WARNING: variable ",                    &
                 trim(dname)," not found in file."
    istat = 1
    RETURN
  END IF

  sds_id = sfselect(sd_id, sds_index)

  attr_index = sffattr(sds_id, "packed16")
  IF (attr_index >= 0) THEN
    istat1 = sfrnatt(sds_id, attr_index, packed16)
    attr_index = sffattr(sds_id, "max")
    istat2 = sfrnatt(sds_id, attr_index, hmax)
    attr_index = sffattr(sds_id, "min")
    istat3 = sfrnatt(sds_id, attr_index, hmin)
    IF (istat1 == -1 .OR. istat2 == -1 .OR. istat3 == -1) THEN
      WRITE (6,*) "HDFRD4D: ERROR reading max/min for ",trim(dname)
      istat = 2
      RETURN
    END IF
  ELSE
    packed16 = 0
  END IF

!-----------------------------------------------------------------------
!
!  Read data into var.
!
!-----------------------------------------------------------------------

  !IF (myproc == 0) WRITE (*,*) "HDFRD4D: Reading variable ",trim(dname)
  IF (packed16 == 0) THEN
    IF (P == SP) THEN
      istat = sfrdata(sds_id, start, stride, dims, var)
    ELSE
      ALLOCATE(dtain(nx,ny,nz,nl), STAT = istatus)
      istat = sfrdata(sds_id, start, stride, dims, dtain)
      var(:,:,:,:) = dtain(:,:,:,:)
      DEALLOCATE(dtain)
    END IF
  ELSE
    istat = sfrdata(sds_id, start, stride, dims, itmp)
  END IF
  IF (istat == -1) THEN
    WRITE (6,*) "HDFRD4D: ERROR reading variable ",trim(dname),"."
    istat = 2
    RETURN
  END IF

!-----------------------------------------------------------------------
!
!  If called for, map 16 bit integers to reals
!
!-----------------------------------------------------------------------

  IF (packed16 /= 0) THEN
    DO l=1,nl
      DO k=1,nz
        scalef = (hmax(k)-hmin(k)) / 65534.0
        DO j=1,ny
          DO i=1,nx
            var(i,j,k,l) = scalef * (itmp(i,j,k,l) + 32767) + hmin(k)
          END DO
        END DO
      END DO
    END DO
  END IF

  istat = sfendacc(sds_id)
  IF (istat /= 0) THEN
    WRITE (6,*) "HDFRD4D: ERROR reading variable ",trim(dname)
    istat = 2
  END IF

  RETURN
END SUBROUTINE hdfrd4d

!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE HDFRD3DI                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################

SUBROUTINE hdfrd3di(sd_id,name,nx,ny,nz,var,istat)

!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Read in a 3-D integer array from an HDF4 file.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Gene Bassett
!  2000/03/15
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    name     Variable name
!
!    sd_id    HDF id of the output file
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!  OUTPUT:
!
!    var      Array to be read in
!
!    istat     Status of read (0-okay, 1-error when reading)
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------

  IMPLICIT NONE
  INTEGER :: nx,ny,nz
  INTEGER :: var(nx,ny,nz)
  INTEGER :: sd_id
  CHARACTER (LEN=*) :: name
  INTEGER :: istat

!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------

  INTEGER :: dims(3),start(3),stride(3)
  INTEGER :: sds_index,sds_id

!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------

  INCLUDE 'hdf.f90'    ! HDF4 library include file
  INCLUDE 'mp.inc'

!-----------------------------------------------------------------------
!
!  Functions
!
!-----------------------------------------------------------------------

  INTEGER :: sfn2index, sfselect, sfrdata, sfendacc

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!-----------------------------------------------------------------------
!
!  Initialize dimension parameters
!
!-----------------------------------------------------------------------

  dims(1) = nx
  dims(2) = ny
  dims(3) = nz
  start(1) = 0
  start(2) = 0
  start(3) = 0
  stride(1) = 1
  stride(2) = 1
  stride(3) = 1

!-----------------------------------------------------------------------
!
!  Get the SDS ID for the variable.
!
!-----------------------------------------------------------------------

  sds_index = sfn2index(sd_id, trim(name))

  IF (sds_index == -1) THEN
    IF (myproc == 0) &
    WRITE (6,*) "HDFRD3DI: WARNING, variable ",                           &
                 trim(name)," not found in file."
    istat = 1
    RETURN
  END IF

  sds_id = sfselect(sd_id, sds_index)

!-----------------------------------------------------------------------
!
!  Read data into var.
!
!-----------------------------------------------------------------------

! IF (myproc == 0) WRITE (*,*) "HDFRD3DI: Reading variable ",trim(name)
  istat = sfrdata(sds_id, start, stride, dims, var)
  IF (istat == -1) THEN
    IF (myproc == 0) &
    WRITE (6,*) "HDFRD3DI: ERROR reading variable ",trim(name),"."
    istat = 2
  END IF

  istat = sfendacc(sds_id)
  IF (istat /= 0) THEN
    IF (myproc == 0) &
    WRITE (6,*) "HDFRD3DI: ERROR reading variable ",trim(name)
    istat = 2
  END IF

  RETURN
END SUBROUTINE hdfrd3di

!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE HDFRD2D                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################

SUBROUTINE hdfrd2d(sd_id,name,nx,ny,var,istat,itmp)

!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Read in a 2-D real array from an HDF4 file.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Gene Bassett
!  2000/03/15
!
!  MODIFICATION HISTORY:
!
!  02/26/2004 Yunheng Wang
!  Added a working array tem to make sure that the data is read in
!  corectly even the default KIND of REAL is not 4. It is
!  allocated in the heap only when needed. Please note that amax
!  and amin were declared as KIND = 4 floating number.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    name     Variable name
!
!    sd_id    HDF id of the output file
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!
!    itmp     Scratch array for mapping reals to integers (used for
!             some values of hdfcompr)
!
!  OUTPUT:
!
!    var      Array to be read in
!
!    istat     Status of read (0-okay, 1-error when reading)
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------

  USE arps_precision

  IMPLICIT NONE
  INTEGER :: nx,ny
  REAL(P) :: var(nx,ny)
  INTEGER :: sd_id
  CHARACTER (LEN=*)    :: name
  INTEGER              :: istat
  INTEGER (KIND=INT16) :: itmp(nx,ny)

!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------

  INTEGER :: dims(2),start(2),stride(2)
  INTEGER :: sds_index,sds_id,attr_index
  REAL(SP):: amax,amin
  REAL(P) :: scalef
  INTEGER :: istat1,istat2,istat3
  INTEGER :: i,j
  INTEGER :: packed16

  REAL(SP), ALLOCATABLE :: tem(:,:)
  INTEGER :: istatus

!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------

  INCLUDE 'hdf.f90'    ! HDF4 library include file
  INCLUDE 'mp.inc'

!-----------------------------------------------------------------------
!
!  Functions
!
!-----------------------------------------------------------------------

  INTEGER :: sfn2index, sfselect, sfrdata, sffattr, sfrnatt, sfendacc

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!-----------------------------------------------------------------------
!
!  Initialize dimension parameters
!
!-----------------------------------------------------------------------

  dims(1) = nx
  dims(2) = ny
  start(1) = 0
  start(2) = 0
  stride(1) = 1
  stride(2) = 1

!-----------------------------------------------------------------------
!
!  Get the SDS ID for the variable.
!
!-----------------------------------------------------------------------

  sds_index = sfn2index(sd_id, trim(name))

  IF (sds_index == -1) THEN
    WRITE (6,*) "HDFRD2D: WARNING, variable ",                            &
                trim(name)," not found in file."
    istat = 1
    RETURN
  END IF

  sds_id = sfselect(sd_id, sds_index)

  attr_index = sffattr(sds_id, "packed16")
  IF (attr_index >= 0) THEN
    istat1 = sfrnatt(sds_id, attr_index, packed16)
    attr_index = sffattr(sds_id, "max")
    istat2 = sfrnatt(sds_id, attr_index, amax)
    attr_index = sffattr(sds_id, "min")
    istat3 = sfrnatt(sds_id, attr_index, amin)
    IF (istat1 == -1 .OR. istat2 == -1 .OR. istat3 == -1) THEN
      WRITE (6,*) "HDFRD2D: ERROR reading max/min for ",trim(name)
      istat = sfendacc(sds_id)
      istat = 2
      RETURN
    END IF
  ELSE
    packed16 = 0
  END IF

!-----------------------------------------------------------------------
!
!  Read data into var.
!
!-----------------------------------------------------------------------

!  IF (myproc == 0) WRITE (*,*) "HDFRD2D: Reading variable ",trim(name)
  IF (packed16 == 0) THEN
    IF( P == SP) THEN
      istat = sfrdata(sds_id, start, stride, dims, var)
    ELSE
      ALLOCATE(tem(nx,ny), STAT = istatus)
      istat = sfrdata(sds_id, start, stride, dims, tem)
      var(:,:) = tem(:,:)
      DEALLOCATE(tem)
    END IF
  ELSE
    istat = sfrdata(sds_id, start, stride, dims, itmp)
  END IF
  IF (istat == -1) THEN
    WRITE (6,*) "HDFRD2D: ERROR reading variable ",trim(name),"."
    istat = sfendacc(sds_id)
    istat = 2
    RETURN
  END IF

!-----------------------------------------------------------------------
!
!  If called for, map 16 bit integers to reals
!
!-----------------------------------------------------------------------

  IF (packed16 /= 0) THEN
    scalef = (amax - amin) / 65534.0
    DO j=1,ny
      DO i=1,nx
        var(i,j) = scalef * (itmp(i,j) + 32767) + amin
      END DO
    END DO
  END IF

  istat = sfendacc(sds_id)
  IF (istat /= 0) THEN
    WRITE (6,*) "HDFRD2D: ERROR reading variable ",trim(name)
    istat = 2
  END IF

  RETURN
END SUBROUTINE hdfrd2d

!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE HDFRD2DI                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################

SUBROUTINE hdfrd2di(sd_id,name,nx,ny,var,istat)

!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Read in a 2-D integer array from an HDF4 file.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Gene Bassett
!  2000/03/15
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    name     Variable name
!
!    sd_id    HDF id of the output file
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!
!  OUTPUT:
!
!    var      Array to be read in
!
!    istat     Status of read (0-okay, 1-error when reading)
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------

  IMPLICIT NONE
  INTEGER :: nx,ny
  INTEGER :: var(nx,ny)
  INTEGER :: sd_id
  CHARACTER (LEN=*) :: name
  INTEGER :: istat

!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------

  INTEGER :: dims(2),start(2),stride(2)
  INTEGER :: sds_index,sds_id

!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------

  INCLUDE 'hdf.f90'    ! HDF4 library include file
  INCLUDE 'mp.inc'

!-----------------------------------------------------------------------
!
!  Functions
!
!-----------------------------------------------------------------------

  INTEGER :: sfn2index, sfselect, sfrdata, sfendacc

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!-----------------------------------------------------------------------
!
!  Initialize dimension parameters
!
!-----------------------------------------------------------------------

  dims(1) = nx
  dims(2) = ny
  start(1) = 0
  start(2) = 0
  stride(1) = 1
  stride(2) = 1

!-----------------------------------------------------------------------
!
!  Get the SDS ID for the variable.
!
!-----------------------------------------------------------------------

  sds_index = sfn2index(sd_id, trim(name))

  IF (sds_index == -1) THEN
    IF (myproc == 0) &
    WRITE (6,*) "HDFRD2DI: WARNING, variable ",                           &
                trim(name)," not found in file."
    istat = 1
    RETURN
  END IF

  sds_id = sfselect(sd_id, sds_index)

!-----------------------------------------------------------------------
!
!  Read data into var.
!
!-----------------------------------------------------------------------

!  IF (myproc == 0) WRITE (*,*) "HDFRD2DI: Reading variable ",trim(name)
  istat = sfrdata(sds_id, start, stride, dims, var)
  IF (istat == -1) THEN
    IF (myproc == 0) &
    WRITE (6,*) "HDFRD2DI: ERROR reading variable ",trim(name),"."
    istat = sfendacc(sds_id)
    istat = 2
    RETURN
  END IF

  istat = sfendacc(sds_id)
  IF (istat /= 0) THEN
    IF (myproc == 0) &
    WRITE (6,*) "HDFRD2DI: ERROR reading variable ",trim(name)
    istat = 2
  END IF

  RETURN
END SUBROUTINE hdfrd2di

!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE HDFRD1D                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################

SUBROUTINE hdfrd1d(sd_id,name,num,var,istat)

!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Read in a 2-D real array from an HDF4 file.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Gene Bassett
!  2000/03/15
!
!  MODIFICATION HISTORY:
!
!  02/26/2004 Yunheng Wang
!  Added a working array tem to make sure that the data is read in
!  correctly even the default KIND of REAL is not 4.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    name     Variable name
!
!    sd_id    HDF id of the output file
!
!    num      Number of grid points
!
!  OUTPUT:
!
!    var      Array to be read in
!
!    istat     Status of read (0-okay, 1-error when reading)
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
  USE arps_precision

  IMPLICIT NONE
  INTEGER           :: num
  REAL(P)           :: var(num)
  INTEGER           :: sd_id
  CHARACTER (LEN=*) :: name
  INTEGER           :: istat

!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------

  INTEGER :: dims(1),start(1),stride(1)
  INTEGER :: sds_index,sds_id

  REAL(SP), ALLOCATABLE :: tem(:)
  INTEGER               :: istatus

!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------

  INCLUDE 'hdf.f90'    ! HDF4 library include file
  INCLUDE 'mp.inc'

!-----------------------------------------------------------------------
!
!  Functions
!
!-----------------------------------------------------------------------

  INTEGER :: sfn2index, sfselect, sfrdata, sfendacc

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!-----------------------------------------------------------------------
!
!  Initialize dimension parameters
!
!-----------------------------------------------------------------------

  dims(1)   = num
  start(1)  = 0
  stride(1) = 1

!-----------------------------------------------------------------------
!
!  Get the SDS ID for the variable.
!
!-----------------------------------------------------------------------

  sds_index = sfn2index(sd_id, trim(name))

  IF (sds_index == -1) THEN
    IF (myproc == 0) &
    WRITE (6,*) "HDFRD1D: WARNING, variable ",trim(name)," not found in file."
    istat = 1
    RETURN
  END IF

  sds_id = sfselect(sd_id, sds_index)

!-----------------------------------------------------------------------
!
!  Read data into var.
!
!-----------------------------------------------------------------------

!  IF (myproc == 0) WRITE (*,*) "HDFRD1D: Reading variable ",trim(name)
  IF ( P == SP ) THEN
    istat = sfrdata(sds_id, start, stride, dims, var)
  ELSE
    ALLOCATE(tem(num), STAT = istatus)
    istat = sfrdata(sds_id, start, stride, dims, tem)
    var(:) = tem(:)
    DEALLOCATE(tem)
  END IF

  IF (istat == -1) THEN
    WRITE (6,*) "HDFRD1D: ERROR reading variable ",trim(name),"."
    istat = sfendacc(sds_id)
    istat = 2
    RETURN
  END IF

  istat = sfendacc(sds_id)
  IF (istat /= 0) THEN
    WRITE (6,*) "HDFRD1D: ERROR reading variable ",trim(name)
    istat = 2
  END IF

  RETURN
END SUBROUTINE hdfrd1d
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE HDFRD1DI                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################

SUBROUTINE hdfrd1di(sd_id,name,num,var,istat)

!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Read in a 1-D integer array from an HDF4 file.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster, after hdfrd1d
!  10/29/2002
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    name     Variable name
!
!    sd_id    HDF id of the output file
!
!    num      Number of grid points
!
!  OUTPUT:
!
!    var      Integer array to be read in
!
!    istat     Status of read (0-okay, 1-error when reading)
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------

  IMPLICIT NONE
  INTEGER :: num
  INTEGER :: var(num)
  INTEGER :: sd_id
  CHARACTER (LEN=*) :: name
  INTEGER :: istat

!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------

  INTEGER :: dims(1),start(1),stride(1)
  INTEGER :: sds_index,sds_id

!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------

  INCLUDE 'hdf.f90'    ! HDF4 library include file
  INCLUDE 'mp.inc'

!-----------------------------------------------------------------------
!
!  Functions
!
!-----------------------------------------------------------------------

  INTEGER :: sfn2index, sfselect, sfrdata, sfendacc

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!-----------------------------------------------------------------------
!
!  Initialize dimension parameters
!
!-----------------------------------------------------------------------

  dims(1) = num
  start(1) = 0
  stride(1) = 1

!-----------------------------------------------------------------------
!
!  Get the SDS ID for the variable.
!
!-----------------------------------------------------------------------

  sds_index = sfn2index(sd_id, trim(name))

  IF (sds_index == -1) THEN
    IF (myproc == 0) &
    WRITE (6,*) "HDFRD1D: WARNING, variable ",                            &
                trim(name)," not found in file."
    istat = 1
    RETURN
  END IF

  sds_id = sfselect(sd_id, sds_index)

!-----------------------------------------------------------------------
!
!  Read data into var.
!
!-----------------------------------------------------------------------

!  IF (myproc == 0) WRITE (*,*) "HDFRD1D: Reading variable ",trim(name)
  istat = sfrdata(sds_id, start, stride, dims, var)
  IF (istat == -1) THEN
    IF (myproc == 0) &
    WRITE (6,*) "HDFRD1DI: ERROR reading variable ",trim(name),"."
    istat = sfendacc(sds_id)
    istat = 2
    RETURN
  END IF

  istat = sfendacc(sds_id)
  IF (istat /= 0) THEN
    IF (myproc == 0) &
    WRITE (6,*) "HDFRD1DI: ERROR reading variable ",trim(name)
    istat = 2
  END IF

  RETURN
END SUBROUTINE hdfrd1di

!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE HDFRDR                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################

SUBROUTINE hdfrdr(sd_id,dname,val,istat)

!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Read in a real attribute
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Gene Bassett
!  2000/03/15
!
!  MODIFICATION HISTORY:
!
!  02/26/2004 Yunheng Wang
!  Added a temporary variable tem to make sure that the data is read in
!  correctly even the default KIND of REAL is not 4. It is because data
!  was written as dfnt_float32 before.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    name     Variable name
!
!    sd_id    HDF id of the file or variable containing the
!             named attribute
!
!  OUTPUT:
!
!    val      The value of the attribute
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------

  USE arps_precision

  IMPLICIT NONE

  INTEGER,          INTENT(IN)  :: sd_id
  CHARACTER(LEN=*), INTENT(IN)  :: dname
  REAL(P),          INTENT(OUT) :: val
  INTEGER,          INTENT(OUT) :: istat

!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------

  INTEGER  :: attr_index
  REAL(SP) :: temval

!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------

  INCLUDE 'hdf.f90'    ! HDF4 library include file
  INCLUDE 'mp.inc'

!-----------------------------------------------------------------------
!
!  Functions
!
!-----------------------------------------------------------------------

  INTEGER :: sffattr, sfrnatt

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  attr_index = sffattr(sd_id, trim(dname))
  IF (attr_index < 0) THEN
    WRITE (6,*) "HDFRDR: WARNING, variable ",                           &
                trim(dname)," not found in file."
    val   = 0.0
    istat = 1
    RETURN
  END IF
  istat = sfrnatt(sd_id, attr_index, temval)
  val   = temval

  IF (istat == -1) THEN
    WRITE (6,*) "HDFRDR: ERROR reading variable ",trim(dname),"."
    istat = 2
  ELSE
!    IF (myproc == 0)  &
!       WRITE (*,*) "HDFRDR: Read variable ",trim(name),", value =",val
  END IF

  RETURN
END SUBROUTINE hdfrdr

!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE HDFRDI                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################

SUBROUTINE hdfrdi(sd_id,dname,val,istat)

!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Read in an integer attribute
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Gene Bassett
!  2000/03/15
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    name     Variable name
!
!    sd_id    HDF id of the file or variable containing the
!             named attribute
!
!  OUTPUT:
!
!    val      The value of the attribute
!    istat    Status indicator (0-okay, 1-variable not found, 2-read error)
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------

  IMPLICIT NONE
  INTEGER :: val
  INTEGER :: sd_id
  CHARACTER (LEN=*) :: dname

!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------

  INTEGER :: attr_index, istat

!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------

  INCLUDE 'hdf.f90'    ! HDF4 library include file
  INCLUDE 'mp.inc'

!-----------------------------------------------------------------------
!
!  Functions
!
!-----------------------------------------------------------------------

  INTEGER :: sffattr, sfrnatt

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  attr_index = sffattr(sd_id, trim(dname))
  IF (attr_index < 0) THEN
    IF (myproc == 0) WRITE (6,*)                                        &
                'HDFRDI: WARNING, variable ', trim(dname),               &
                ' not found in file. Value was reset to 0.'
    val   = 0
    istat = 1
    RETURN
  END IF
  istat = sfrnatt(sd_id, attr_index, val)

  IF (istat == -1) THEN
    IF (myproc == 0) WRITE (6,*)                                        &
                'HDFRDI: ERROR reading variable ',trim(dname),'.'
    istat = 2
    val = -1
!  ELSE
!    IF (myproc == 0) WRITE (*,*) 'HDFRDI: Read variable ',trim(name),  &
!                                 ' value: ',val
  END IF

  RETURN
END SUBROUTINE hdfrdi

!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE HDFRDC                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################

SUBROUTINE hdfrdc(sd_id,max_len,name,string,istat)

!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Read in a string attribute
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Gene Bassett
!  2000/03/15
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    name     Variable name
!
!    max_len  Maximum allowable length of string to be read in
!
!    sd_id    HDF id of the file or variable containing the
!             named attribute
!
!  OUTPUT:
!
!    string   The value of the attribute
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------

  IMPLICIT NONE
  INTEGER           :: sd_id
  CHARACTER (LEN=*) :: name
  CHARACTER (LEN=*) :: string
  INTEGER           :: max_len

!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------

  INTEGER :: attr_index, istat
  INTEGER :: data_type, n_values
  CHARACTER (LEN=256) :: attr_name
  INTEGER :: i

!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------

  INCLUDE 'hdf.f90'    ! HDF4 library include file
  INCLUDE 'mp.inc'

!-----------------------------------------------------------------------
!
!  Functions
!
!-----------------------------------------------------------------------

  INTEGER :: sffattr, sfgainfo, sfrnatt, &
             sfrcatt

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


  attr_index = sffattr(sd_id, trim(name))

  IF (attr_index < 0) THEN
    IF (myproc == 0) &
    WRITE (6,*) "HDFRDC: WARNING, variable ",trim(name),                  &
                " not located."
    istat = 1
    RETURN
  END IF

  istat = sfgainfo(sd_id, attr_index, attr_name, data_type,             &
                   n_values)
  IF (istat /= 0) THEN
    IF (myproc == 0) &
    WRITE (6,*) "HDFRDC: ERROR in looking up attributes for ",          &
                trim(name)
    istat = 2
    RETURN
  END IF

  IF (n_values <= max_len) THEN
    istat = sfrcatt(sd_id, attr_index, string)
  ELSE
    IF (myproc == 0) &
    WRITE (6,'(1x,3a,I3,a,I3)') "HDFRDC: ERROR: string length for variable ",           &
                trim(name),", ",max_len,                                &
                ", is less than string in file:",n_values,","
    IF (myproc == 0) &
    WRITE (6,'(1x,a)') "value not read in."
  END IF

  DO i=n_values+1,len(string)
    string(i:i) = " "
  END DO

  IF (istat == -1) THEN
    IF (myproc == 0) WRITE (6,*) 'HDFRDC: ERROR reading variable ',trim(name),'.'
    istat = 1
!  ELSE
!    IF (myproc == 0)  &
!       WRITE (*,*) "HDFRDC: Read variable ",trim(name)," value:", trim(string)
  END IF

  RETURN
END SUBROUTINE hdfrdc

!##################################################################
!##################################################################
!######                                                      ######
!######                  SUBROUTINE HDFOPEN                  ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################

SUBROUTINE hdfopen(filename,rdwrtflg,sd_id)

!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Open a HDF file.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Gene Bassett
!  2000/04/11
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    filename     File name
!
!    rdwrtflg Read/write flag (1-read, 2-write)
!
!  OUTPUT:
!
!    sd_id    HDF id of the file or variable containing the
!             named attribute
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------

  IMPLICIT NONE
  CHARACTER (LEN=*), INTENT(IN)  :: filename
  INTEGER,           INTENT(IN)  :: rdwrtflg
  INTEGER,           INTENT(OUT) :: sd_id

!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------

  INCLUDE 'hdf.f90'    ! HDF4 library include file
  INCLUDE 'mp.inc'     ! mpi parameters

!-----------------------------------------------------------------------
!
!  Functions
!
!-----------------------------------------------------------------------

  INTEGER :: sfstart

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  IF (rdwrtflg == 1) THEN
    sd_id = sfstart(trim(filename), dfacc_read)
  ELSE IF (rdwrtflg == 2) THEN
    sd_id = sfstart(trim(filename), dfacc_create)
  ELSE IF (rdwrtflg == 3) THEN
    sd_id = sfstart(trim(filename), dfacc_write)
  ELSE
    sd_id = -1
    IF (myproc == 0) WRITE (6,*)                                        &
      'HDFOPEN: ERROR, unsupported rdwrtflg value of ', rdwrtflg
  END IF

  IF (sd_id <= 0 .AND. myproc == 0 ) THEN
    WRITE (6,*) "HDFOPEN: ERROR opening ",trim(filename)
  END IF

  RETURN
END SUBROUTINE hdfopen

!##################################################################
!##################################################################
!######                                                      ######
!######                  SUBROUTINE HDFCLOSE                 ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################

SUBROUTINE hdfclose(sd_id,istat)

!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Close a HDF file.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Gene Bassett
!  2000/04/11
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    sd_id    HDF id of file to close
!
!  OUTPUT:
!
!    istat     Status returned by close
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------

  IMPLICIT NONE
  INTEGER :: sd_id
  INTEGER :: istat

!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------

  INCLUDE 'hdf.f90'    ! HDF4 library include file

!-----------------------------------------------------------------------
!
!  Functions
!
!-----------------------------------------------------------------------

  INTEGER :: sfend

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  istat = sfend(sd_id)

  RETURN
END SUBROUTINE hdfclose

!##################################################################
!##################################################################
!######                                                      ######
!######                  SUBROUTINE HDFINFO                  ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################

SUBROUTINE hdfinfo(sd_id,ndata,nattr,istat)

!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Get the number of data sets and attributes contained in a HDF data file.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Gene Bassett
!  2000/10/18
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    sd_id     HDF id of the data file
!
!  OUTPUT:
!
!    ndata     Number of data sets in the file
!    nattr     Number of attributes in the file
!    istat     Status returned by sffinfo
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------

  IMPLICIT NONE
  INTEGER :: sd_id
  INTEGER :: ndata,nattr,istat

!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------

  INCLUDE 'hdf.f90'    ! HDF4 library include file

!-----------------------------------------------------------------------
!
!  Functions
!
!-----------------------------------------------------------------------

  INTEGER :: sffinfo

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  istat = sffinfo(sd_id,ndata,nattr)

  RETURN
END SUBROUTINE hdfinfo

!##################################################################
!##################################################################
!######                                                      ######
!######                  SUBROUTINE HDFDINFO                 ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################

SUBROUTINE hdfdinfo(sds_id,name,rank,dims,dtype,nattr,istat)

!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Get the information about a data set
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Gene Bassett
!  2000/10/18
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    sds_id    HDF id of the data set
!
!  OUTPUT:
!
!    name      Name of the data set
!    rank      Number of dimensions in the data set
!    dims      Number of points for each dimension
!    dtype     Data type
!    nattr     Number of attributes in the data set
!    istat     Status returned by sfginfo
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------

  IMPLICIT NONE
  INTEGER :: sds_id
  CHARACTER (LEN=*) :: name
  INTEGER :: rank,dims(6),dtype,nattr,istat

!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------

  INCLUDE 'hdf.f90'    ! HDF4 library include file

!-----------------------------------------------------------------------
!
!  Functions
!
!-----------------------------------------------------------------------

  INTEGER :: sfginfo

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  istat = sfginfo(sds_id,name,rank,dims,dtype,nattr)

  RETURN
END SUBROUTINE hdfdinfo

!##################################################################
!##################################################################
!######                                                      ######
!######                  SUBROUTINE HDFAINFO                 ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################

SUBROUTINE hdfainfo(sd_id,aindex,name,dtype,nvalues,istat)

!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Get the information about an attribute.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Gene Bassett
!  2000/10/25
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    sd_id     HDF id of the data set (for HDF file or SDS data set in a file)
!    aindex    Index number for the attribute
!
!  OUTPUT:
!
!    name      Name of the data set
!    dtype     Data type
!    nvalues   Number of values (number of characters if a string)
!    istat     Status returned by sfgainfo
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------

  IMPLICIT NONE
  INTEGER :: sd_id,aindex
  CHARACTER (LEN=*) :: name
  INTEGER :: dtype,nvalues,istat

!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------

  INCLUDE 'hdf.f90'    ! HDF4 library include file

!-----------------------------------------------------------------------
!
!  Functions
!
!-----------------------------------------------------------------------

  INTEGER :: sfgainfo

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  istat = sfgainfo(sd_id,aindex,name,dtype,nvalues)

  RETURN
END SUBROUTINE hdfainfo

!##################################################################
!##################################################################
!######                                                      ######
!######              SUBROUTINE GET_DIMS_FROM_HDF            ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################

SUBROUTINE get_dims_from_hdf(filename,nx,ny,nz,nzsoil,nstyps,ireturn)

!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Get the grid dimensions from an ARPS hdf file.  (Similar to
!  get_dims_from_data.)
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Gene Bassett
!  2001/04/23
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    filename  File name
!
!  OUTPUT:
!
!    nx,ny,nz  Grid dimensions
!    nstyps    Number of soil types
!    ireturn   Return status indicator
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------

  IMPLICIT NONE
  CHARACTER (LEN=*) :: filename
!  INTEGER :: nx,ny,nz,nstyps,ireturn
  INTEGER :: nx,ny,nz,nzsoil,nstyps,ireturn


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

!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------

  INTEGER :: sd_id,istat

!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!  Functions
!
!-----------------------------------------------------------------------

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  ireturn = 0

  CALL hdfopen(filename,1,sd_id)

  IF (sd_id <= 0) THEN
    ireturn = 1
    RETURN
  ENDIF

  CALL hdfrdc(sd_id,40,"fmtver",fmtverin,istat)

  IF (fmtverin == fmtver410) THEN
    intver=intver410
  ELSE IF (fmtverin == fmtver500) THEN
    intver=intver500
  ELSE IF (fmtverin == fmtver530) THEN
    intver=intver530
  ELSE
    WRITE(6,'(/1x,a,a,a/)')                        &
        'Incoming data format, fmtverin=',fmtverin,                 &
        ', not found. The Job stopped.'
    CALL arpsstop('arpstop called from HDFREAD. ',1)
  END IF

  WRITE(6,'(/1x,a,a/)')'Incoming data format, fmtverin=',fmtverin

  CALL hdfrdi(sd_id,"nx",nx,istat)
  IF (istat > 0) ireturn = 1
  CALL hdfrdi(sd_id,"ny",ny,istat)
  IF (istat > 0) ireturn = 1
  CALL hdfrdi(sd_id,"nz",nz,istat)
  IF (istat > 0) ireturn = 1
  CALL hdfrdi(sd_id,"nstyp",nstyps,istat)
  IF (istat > 0) ireturn = 1

  IF (intver >= intver500) THEN
    CALL hdfrdi(sd_id,"nzsoil",nzsoil,istat)
    IF (istat > 0) ireturn = 1
  ELSE
    nzsoil = 2
  END IF

  CALL hdfclose(sd_id,istat)

  RETURN
END SUBROUTINE get_dims_from_hdf

SUBROUTINE get_var_attr_from_hdf(sd_id,varname,attrname,attrvalue,vallen,istatus)

!-----------------------------------------------------------------------
!
! Given a variable name, find its string attributes from HDF file
!
!      long name
!      units
!
!-----------------------------------------------------------------------

  IMPLICIT NONE
  INTEGER,      INTENT(IN)  :: sd_id
  CHARACTER(*), INTENT(IN)  :: varname
  CHARACTER(*), INTENT(IN)  :: attrname
  CHARACTER(*), INTENT(OUT) :: attrvalue
  INTEGER,      INTENT(IN)  :: vallen
  INTEGER,      INTENT(OUT) :: istatus

  INCLUDE 'hdf.f90'

  INTEGER  :: sds_index, sds_id

  INTEGER  :: sfn2index, sfselect

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Executable code below ...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  sds_index = sfn2index(sd_id, trim(varname))

  IF (sds_index == -1) THEN
    WRITE (6,*) "get_var_attr: WARNING, variable ", trim(varname)," not found in file."
    istatus = 1
    RETURN
  END IF

  sds_id = sfselect(sd_id, sds_index)

  CALL hdfrdc(sds_id,vallen,attrname,attrvalue, istatus)

  RETURN
END SUBROUTINE get_var_attr_from_hdf

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Subroutines that work with module arps_dtaread
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE HDF_READ_3D                ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################

SUBROUTINE hdf_read_3d(IAMROOT, IREAD, ISPLIT, ncompressx, ncompressy,  &
           fileid,nx,ny,nz,nzsoil,nstyps,nscalar, qnames,nxdta,nydta,   &
           time,u, v, w, pt, prs,                                       &
           qv, qscalar, tke, kmh,kmv,                                   &
           tsoil,qsoil,wetcanp,snowdpth,                                &
           raing,rainc,prcrate, lvldbg,                                 &
           ireturn)

!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Read in history data in the NCSA HDF4 format. Read in atmospheric state
!  arrays only.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yunheng Wang
!  2011/06/01
!  Based on subroutine hdfread.
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
!    nstyps
!
!  OUTPUT:
!
!    time     Time in seconds of data read from "filename"
!
!    u        x component of total velocity (m/s)
!    v        y component of total velocity (m/s)
!    w        Vertical component of total velocity in Cartesian
!             coordinates (m/s).
!    pt       Total potential temperature (K)
!    prs      Total pressure (Pascal)
!
!    qv       Perturbation water vapor mixing ratio (kg/kg)
!    qscalar  Hydrometeors
!    tke      Turbulent Kinetic Energy ((m/s)**2)
!    kmh      Horizontal turb. mixing coef. for momentum ( m**2/s )
!    kmv      Vertical turb. mixing coef. for momentum ( m**2/s )
!
!    tsoil    Soil temperature (K)
!    qsoil    Soil moisture (m**3/m**3)
!    wetcanp  Canopy water amount
!
!    raing    Grid supersaturation rain
!    rainc    Cumulus convective rain
!    prcrate  Precipitation rates
!
!    ireturn  Return status indicator
!
!
!-----------------------------------------------------------------------

  USE arps_precision
  IMPLICIT NONE

!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------

  LOGICAL, INTENT(IN) :: IAMROOT, IREAD, ISPLIT
  INTEGER, INTENT(IN) :: ncompressx, ncompressy
  INTEGER, INTENT(IN) :: fileid(ncompressx,ncompressy)

  INTEGER, INTENT(IN) :: nx,ny,nz,nscalar
  INTEGER, INTENT(IN) :: nzsoil
  INTEGER, INTENT(IN) :: nstyps                    ! Number of soil type
  INTEGER, INTENT(IN) :: nxdta, nydta

  CHARACTER(LEN=4), INTENT(IN) :: qnames(nscalar)

  REAL(P), INTENT(OUT) :: u  (nx,ny,nz)     ! u-velocity (m/s)
  REAL(P), INTENT(OUT) :: v  (nx,ny,nz)     ! v-velocity (m/s)
  REAL(P), INTENT(OUT) :: w  (nx,ny,nz)     ! w-velocity (m/s)
  REAL(P), INTENT(OUT) :: pt (nx,ny,nz)     ! potential temperature (K)
  REAL(P), INTENT(OUT) :: prs(nx,ny,nz)     ! pressure (Pascal)

  REAL(P), INTENT(OUT) :: qv (nx,ny,nz)     ! water vapor mixing ratio (kg/kg)

  REAL(P), INTENT(OUT) :: qscalar(nx,ny,nz,nscalar)

  REAL(P), INTENT(OUT) :: tke   (nx,ny,nz)     ! Turbulent Kinetic Energy ((m/s)**2)
  REAL(P), INTENT(OUT) :: kmh   (nx,ny,nz)     ! Horizontal turb. mixing coef. for
                               ! momentum. ( m**2/s )
  REAL(P), INTENT(OUT) :: kmv   (nx,ny,nz)     ! Vertical turb. mixing coef. for
                               ! momentum. ( m**2/s )
  REAL(P), INTENT(OUT) :: tsoil (nx,ny,nzsoil,0:nstyps)  ! Soil temperature (K)
  REAL(P), INTENT(OUT) :: qsoil (nx,ny,nzsoil,0:nstyps)  ! Soil moisture (m**3/m**3)
  REAL(P), INTENT(OUT) :: wetcanp(nx,ny,0:nstyps)   ! Canopy water amount
  REAL(P), INTENT(OUT) :: snowdpth(nx,ny)           ! Snow depth (m)

  REAL(P), INTENT(OUT) :: raing(nx,ny)         ! Grid supersaturation rain
  REAL(P), INTENT(OUT) :: rainc(nx,ny)         ! Cumulus convective rain
  REAL(P), INTENT(OUT) :: prcrate(nx,ny,4)     ! precipitation rates (kg/(m**2*s))
                               ! prcrate(1,1,1) = total precip. rate
                               ! prcrate(1,1,2) = grid scale precip. rate
                               ! prcrate(1,1,3) = cumulative precip. rate
                               ! prcrate(1,1,4) = microphysics precip. rate


  REAL(P), INTENT(OUT) :: time

  INTEGER, INTENT(IN)  :: lvldbg
  INTEGER, INTENT(OUT) :: ireturn

!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------


  INTEGER :: i,j,k,is,nq
  INTEGER :: ia,ja
  INTEGER :: ii, jj

  INTEGER :: istatus
  INTEGER :: totin

  REAL(P), ALLOCATABLE :: dta3d(:,:,:), dta4d(:,:,:,:)
  REAL(P), ALLOCATABLE :: tem1(:,:,:), temsoil(:,:,:,:)

  INTEGER :: nztmp

  INTEGER(INT16), ALLOCATABLE :: itmp(:,:,:)              ! Temporary array
  REAL(SP),       ALLOCATABLE :: hmax(:), hmin(:)         ! Temporary array
  REAL(SP),       ALLOCATABLE :: hmaxsoil(:), hminsoil(:) ! Temporary array
  INTEGER(INT16), ALLOCATABLE :: itmpsoil(:,:,:,:)


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  nztmp = max(nz,nstyps+1)
  ALLOCATE (itmp(nxdta,nydta,nztmp),stat=istatus)
  CALL check_alloc_status(istatus, "hdf_read_3d:itmp")

  ALLOCATE (hmax(nztmp),stat=istatus)
  CALL check_alloc_status(istatus, "hdf_read_3d:hmax")

  ALLOCATE (hmin(nztmp),stat=istatus)
  CALL check_alloc_status(istatus, "hdf_read_3d:hmin")

  ALLOCATE (dta3d(nxdta,nydta,MAX(nz,nzsoil,nstyps+1,4)),stat=istatus)
  CALL check_alloc_status(istatus, "hdf_read_3d:dta3d")

  ALLOCATE (tem1(MAX(nx,nxdta),MAX(ny,nydta),nztmp),stat=istatus)
  CALL check_alloc_status(istatus, "hdf_read_3d:tem1")


  ALLOCATE (itmpsoil(nxdta,nydta,nzsoil,0:nstyps),stat=istatus)
  CALL check_alloc_status(istatus, "hdf_read_3d:itmpsoil")

  ALLOCATE (hmaxsoil(nzsoil),stat=istatus)
  CALL check_alloc_status(istatus, "hdf_read_3d:hmaxsoil")

  ALLOCATE (hminsoil(nzsoil),stat=istatus)
  CALL check_alloc_status(istatus, "hdf_read_3d:hminsoil")

  ALLOCATE (dta4d(nxdta,nydta,nzsoil,nstyps+1),stat=istatus)
  CALL check_alloc_status(istatus, "hdf_read_3d:dta4d")

  ALLOCATE (temsoil(MAX(nx,nxdta),MAX(ny,nydta),nzsoil,0:nstyps),stat=istatus)
  CALL check_alloc_status(istatus, "hdf_read_3d:temsoil")

!-----------------------------------------------------------------------
!
!  Read in current history file time
!
!-----------------------------------------------------------------------

  IF (IREAD) THEN
    CALL hdfrdr(fileid(1,1),"time",time,istatus)
    CALL hdfrdi(fileid(1,1),"totflg",totin,istatus)
  END IF
  IF (ISPLIT) THEN
    CALL mpupdater(time,1)
    CALL mpupdatei(totin,1)
  END IF

  IF (totin /= 1) THEN
    WRITE(*,'(1x,a,I2)') 'ERROR: Expect totflag = 1, but find totin = ',totin
    CALL arpsstop('ERROR wrong flag totout',1)
  END IF

!-----------------------------------------------------------------------
!
!  Read in total values of variables from history dump
!
!-----------------------------------------------------------------------

  CALL hdfread3d(IREAD,ISPLIT,ncompressx,ncompressy,fileid,"u",         &
                 nx,ny,nz,u,                                            &
                 nxdta,nydta,dta3d,tem1,itmp,hmax,hmin,istatus)

  CALL hdfread3d(IREAD,ISPLIT,ncompressx,ncompressy,fileid,"v",         &
                 nx,ny,nz,v,                                            &
                 nxdta,nydta,dta3d,tem1,itmp,hmax,hmin,istatus)

  CALL hdfread3d(IREAD,ISPLIT,ncompressx,ncompressy,fileid,"w",         &
                 nx,ny,nz,w,                                            &
                 nxdta,nydta,dta3d,tem1,itmp,hmax,hmin,istatus)

  CALL hdfread3d(IREAD,ISPLIT,ncompressx,ncompressy,fileid,"pt",        &
                 nx,ny,nz,pt,                                           &
                 nxdta,nydta,dta3d,tem1,itmp,hmax,hmin,istatus)

  CALL hdfread3d(IREAD,ISPLIT,ncompressx,ncompressy,fileid,"p",         &
                 nx,ny,nz,prs,                                          &
                 nxdta,nydta,dta3d,tem1,itmp,hmax,hmin,istatus)

!-----------------------------------------------------------------------
!
!  Read in moisture variables
!
!-----------------------------------------------------------------------

  CALL hdfread3d(IREAD,ISPLIT,ncompressx,ncompressy,fileid,"qv",        &
                 nx,ny,nz,qv,                                           &
                 nxdta,nydta,dta3d,tem1,itmp,hmax,hmin,istatus)

  DO nq = 1,nscalar
    CALL hdfread3d(IREAD,ISPLIT,ncompressx,ncompressy,fileid,qnames(nq),&
                   nx,ny,nz,qscalar(:,:,:,nq),                          &
                   nxdta,nydta,dta3d,tem1,itmp,hmax,hmin,istatus)
  END DO

  CALL hdfread2d(IREAD,ISPLIT,ncompressx,ncompressy,fileid,"raing",     &
                 nx,ny,raing,                                           &
                 nxdta,nydta,dta3d,tem1,itmp,istatus)
  IF (istatus /=0 .AND. lvldbg > 0) WRITE(*,'(1x,a)') 'WARNING: Raing is not found.'

  CALL hdfread2d(IREAD,ISPLIT,ncompressx,ncompressy,fileid,"rainc",     &
                 nx,ny,rainc,                                           &
                 nxdta,nydta,dta3d,tem1,itmp,istatus)
  IF (istatus /=0 .AND. lvldbg > 0) WRITE(*,'(1x,a)') 'WARNING: Rainc is not found.'

  CALL hdfread2d(IREAD,ISPLIT,ncompressx,ncompressy,fileid,"prcrate1",  &
                 nx,ny,prcrate(:,:,1),                                  &
                 nxdta,nydta,dta3d,tem1,itmp,istatus)
  IF (istatus /=0 .AND. lvldbg > 0) WRITE(*,'(1x,a)') 'WARNING: prcrate1 is not found.'

  CALL hdfread2d(IREAD,ISPLIT,ncompressx,ncompressy,fileid,"prcrate2",  &
                 nx,ny,prcrate(:,:,2),                                  &
                 nxdta,nydta,dta3d,tem1,itmp,istatus)
  IF (istatus /=0 .AND. lvldbg > 0) WRITE(*,'(1x,a)') 'WARNING: prcrate2 is not found.'

  CALL hdfread2d(IREAD,ISPLIT,ncompressx,ncompressy,fileid,"prcrate3",  &
                 nx,ny,prcrate(:,:,3),                                  &
                 nxdta,nydta,dta3d,tem1,itmp,istatus)
  IF (istatus /=0 .AND. lvldbg > 0) WRITE(*,'(1x,a)') 'WARNING: prcrate3 is not found.'

  CALL hdfread2d(IREAD,ISPLIT,ncompressx,ncompressy,fileid,"prcrate4",  &
                 nx,ny,prcrate(:,:,4),                                  &
                 nxdta,nydta,dta3d,tem1,itmp,istatus)
  IF (istatus /=0 .AND. lvldbg > 0) WRITE(*,'(1x,a)') 'WARNING: prcrate4 is not found.'

  CALL hdfread3d(IREAD,ISPLIT,ncompressx,ncompressy,fileid,"tke",       &
                 nx,ny,nz,tke,                                          &
                 nxdta,nydta,dta3d,tem1,itmp,hmax,hmin,istatus)
  IF (istatus /=0 .AND. lvldbg > 0) WRITE(*,'(1x,a)') 'WARNING: tke is not found.'

  CALL hdfread3d(IREAD,ISPLIT,ncompressx,ncompressy,fileid,"kmh",       &
                 nx,ny,nz,kmh,                                          &
                 nxdta,nydta,dta3d,tem1,itmp,hmax,hmin,istatus)
  IF (istatus /=0 .AND. lvldbg > 0) WRITE(*,'(1x,a)') 'WARNING: kmh is not found.'

  CALL hdfread3d(IREAD,ISPLIT,ncompressx,ncompressy,fileid,"kmv",       &
                 nx,ny,nz,kmv,                                          &
                 nxdta,nydta,dta3d,tem1,itmp,hmax,hmin,istatus)
  IF (istatus /=0 .AND. lvldbg > 0) WRITE(*,'(1x,a)') 'WARNING: kmv is not found.'

  CALL hdfread4d(IREAD,ISPLIT,ncompressx,ncompressy,fileid,"tsoil",     &
                 nx,ny,nzsoil,nstyps+1,tsoil,                           &
                 nxdta,nydta,dta4d,temsoil,itmpsoil,hmaxsoil,hminsoil,istatus)
  IF (istatus /=0 .AND. lvldbg > 0) WRITE(*,'(1x,a)') 'WARNING: tsoil is not found.'

  CALL hdfread4d(IREAD,ISPLIT,ncompressx,ncompressy,fileid,"qsoil",     &
                 nx,ny,nzsoil,nstyps+1,qsoil,                           &
                 nxdta,nydta,dta4d,temsoil,itmpsoil,hmaxsoil,hminsoil,istatus)
  IF (istatus /=0 .AND. lvldbg > 0) WRITE(*,'(1x,a)') 'WARNING: qsoil is not found.'


  CALL hdfread3d(IREAD,ISPLIT,ncompressx,ncompressy,fileid,"wetcanp",   &
                 nx,ny,nstyps+1,wetcanp,                                &
                 nxdta,nydta,dta3d,tem1,itmp,hmax,hmin,istatus)
  IF (istatus /=0 .AND. lvldbg > 0) WRITE(*,'(1x,a)') 'WARNING: wetcanp is not found.'

  !CALL fix_soil_nstyp(nx,ny,nzsoil,nstyp1,nstyp,tsoil,qsoil,wetcanp)

  CALL hdfread2d(IREAD,ISPLIT,ncompressx,ncompressy,fileid,"snowdpth",  &
                 nx,ny,snowdpth,                                        &
                 nxdta,nydta,dta3d,tem1,itmp,istatus)
  IF (istatus /=0 .AND. lvldbg > 0) WRITE(*,'(1x,a)') 'WARNING: snowdpth is not found.'

!-----------------------------------------------------------------------
!
!  Friendly exit message
!
!-----------------------------------------------------------------------

  ireturn = 0

  DEALLOCATE (itmp,hmax,hmin, stat=istatus)
  DEALLOCATE (itmpsoil,hmaxsoil,hminsoil, stat=istatus)

  DEALLOCATE (dta3d, dta4d)
  DEALLOCATE (tem1, temsoil)

  RETURN
END SUBROUTINE hdf_read_3d

!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE HDF_READ_BASE              ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################

SUBROUTINE hdf_read_base(IAMROOT, IREAD, ISPLIT, ncompressx, ncompressy,&
                         fileid,nx,ny,nz,nxdta,nydta,                   &
                         ubar, vbar, wbar, ptbar, pbar, qvbar,ireturn)

!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Read in atmospheric base state arrays only.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yunheng Wang
!  2011/06/02
!  Based on subroutine hdfread.
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
!  OUTPUT:
!
!    ubar     Base state x velocity component (m/s)
!    vbar     Base state y velocity component (m/s)
!    wbar     Base state z velocity component (m/s)
!    ptbar    Base state potential temperature (K)
!    pbar     Base state pressure (Pascal)
!    qvbar    Base state water vapor mixing ratio (kg/kg)
!
!    ireturn  Return status indicator
!
!
!-----------------------------------------------------------------------

  USE arps_precision
  IMPLICIT NONE

!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------

  LOGICAL, INTENT(IN) :: IAMROOT, IREAD, ISPLIT
  INTEGER, INTENT(IN) :: ncompressx, ncompressy
  INTEGER, INTENT(IN) :: fileid(ncompressx,ncompressy)

  INTEGER, INTENT(IN) :: nx,ny,nz
  INTEGER, INTENT(IN) :: nxdta, nydta

  REAL(P), INTENT(OUT) :: ubar  (nx,ny,nz)     ! Base state u-velocity (m/s)
  REAL(P), INTENT(OUT) :: vbar  (nx,ny,nz)     ! Base state v-velocity (m/s)
  REAL(P), INTENT(OUT) :: wbar  (nx,ny,nz)     ! Base state w-velocity (m/s)
  REAL(P), INTENT(OUT) :: ptbar (nx,ny,nz)     ! Base state potential temperature (K)
  REAL(P), INTENT(OUT) :: pbar  (nx,ny,nz)     ! Base state pressure (Pascal).
  REAL(P), INTENT(OUT) :: qvbar (nx,ny,nz)     ! Base state water vapor mixing ratio

  INTEGER, INTENT(OUT) :: ireturn

!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------


  INTEGER :: i,j,k
  INTEGER :: ia,ja
  INTEGER :: ii, jj

  INTEGER :: istatus

  REAL(SP), ALLOCATABLE :: dta3d(:,:,:)
  REAL(P),  ALLOCATABLE :: tem1(:,:,:)

  INTEGER(INT16), ALLOCATABLE :: itmp(:,:,:)              ! Temporary array
  REAL(SP),       ALLOCATABLE :: hmax(:), hmin(:)         ! Temporary array


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  ALLOCATE (itmp(nxdta,nydta,nz),stat=istatus)
  CALL check_alloc_status(istatus, "hdf_read_base:itmp")

  ALLOCATE (hmax(nz),stat=istatus)
  CALL check_alloc_status(istatus, "hdf_read_base:hmax")

  ALLOCATE (hmin(nz),stat=istatus)
  CALL check_alloc_status(istatus, "hdf_read_base:hmin")

  ALLOCATE (dta3d(nxdta,nydta,nz),stat=istatus)
  CALL check_alloc_status(istatus, "hdf_read_base:dta3d")

  ALLOCATE (tem1(MAX(nx,nxdta),MAX(ny,nydta),nz),stat=istatus)
  CALL check_alloc_status(istatus, "hdf_read_base:tem1")

!-----------------------------------------------------------------------
!
!  Read in base state fields
!
!-----------------------------------------------------------------------

  CALL hdfread3d(IREAD,ISPLIT,ncompressx,ncompressy,fileid,"ubar",      &
                 nx,ny,nz,ubar,                                         &
                 nxdta,nydta,dta3d,tem1,itmp,hmax,hmin,istatus)

  CALL hdfread3d(IREAD,ISPLIT,ncompressx,ncompressy,fileid,"vbar",      &
                 nx,ny,nz,vbar,                                         &
                 nxdta,nydta,dta3d,tem1,itmp,hmax,hmin,istatus)

  CALL hdfread3d(IREAD,ISPLIT,ncompressx,ncompressy,fileid,"wbar",      &
                 nx,ny,nz,wbar,                                         &
                 nxdta,nydta,dta3d,tem1,itmp,hmax,hmin,istatus)

  CALL hdfread3d(IREAD,ISPLIT,ncompressx,ncompressy,fileid,"ptbar",     &
                 nx,ny,nz,ptbar,                                        &
                 nxdta,nydta,dta3d,tem1,itmp,hmax,hmin,istatus)

  CALL hdfread3d(IREAD,ISPLIT,ncompressx,ncompressy,fileid,"pbar",      &
                 nx,ny,nz,pbar,                                         &
                 nxdta,nydta,dta3d,tem1,itmp,hmax,hmin,istatus)

  CALL hdfread3d(IREAD,ISPLIT,ncompressx,ncompressy,fileid,"qvbar",     &
                 nx,ny,nz,qvbar,                                        &
                 nxdta,nydta,dta3d,tem1,itmp,hmax,hmin,istatus)

!-----------------------------------------------------------------------
!
!  Friendly exit message
!
!-----------------------------------------------------------------------

  ireturn = istatus

  DEALLOCATE (itmp,hmax,hmin, stat=istatus)
  DEALLOCATE (dta3d, tem1)

  RETURN
END SUBROUTINE hdf_read_base

!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE HDF_READ_GRID              ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################

SUBROUTINE hdf_read_grid(IAMROOT, IREAD, ISPLIT, ncompressx, ncompressy,&
                         fileid,nx,ny,nz,nzsoil,nxdta,nydta,            &
                         x,y,z,zp,zpsoil,ireturn)

!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Read in arps grid arrays only.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yunheng Wang
!  2011/06/02
!  Based on subroutine hdfread.
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
!  OUTPUT:
!
!    x        x coordinate of grid points in physical/comp. space (m)
!    y        y coordinate of grid points in physical/comp. space (m)
!    z        z coordinate of grid points in computational space (m)
!    zp       z coordinate of grid points in physical space (m)
!    zpsoil   z coordinate of grid points in the soil (m)
!
!    ireturn  Return status indicator
!
!
!-----------------------------------------------------------------------

  USE arps_precision
  IMPLICIT NONE

!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------

  LOGICAL, INTENT(IN) :: IAMROOT, IREAD, ISPLIT
  INTEGER, INTENT(IN) :: ncompressx, ncompressy
  INTEGER, INTENT(IN) :: fileid(ncompressx,ncompressy)

  INTEGER, INTENT(IN) :: nx,ny,nz,nzsoil
  INTEGER, INTENT(IN) :: nxdta, nydta

  REAL(P), INTENT(OUT) :: x     (nx)           ! x coord.
  REAL(P), INTENT(OUT) :: y     (ny)           ! y coord.
  REAL(P), INTENT(OUT) :: z     (nz)           ! z coord.
  REAL(P), INTENT(OUT) :: zp    (nx,ny,nz)     ! physical x coord.
  REAL(P), INTENT(OUT) :: zpsoil(nx,ny,nzsoil) ! physical x coord. for soil (m)

  INTEGER, INTENT(OUT) :: ireturn

!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------


  INTEGER :: i,j,k
  INTEGER :: ia,ja
  INTEGER :: ii, jj

  INTEGER :: istatus

  INTEGER :: nztmp

  REAL(SP), ALLOCATABLE :: dta3d(:,:,:)
  REAL(P),  ALLOCATABLE :: tem1(:,:,:)

  INTEGER(INT16), ALLOCATABLE :: itmp(:,:,:)              ! Temporary array
  REAL(SP),       ALLOCATABLE :: hmax(:), hmin(:)         ! Temporary array


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  nztmp = MAX(nz,nzsoil)

  ALLOCATE (itmp(nxdta,nydta,nztmp),stat=istatus)
  CALL check_alloc_status(istatus, "hdf_read_base:itmp")

  ALLOCATE (hmax(nztmp),stat=istatus)
  CALL check_alloc_status(istatus, "hdf_read_base:hmax")

  ALLOCATE (hmin(nztmp),stat=istatus)
  CALL check_alloc_status(istatus, "hdf_read_base:hmin")

  ALLOCATE (dta3d(nxdta,nydta,nztmp),stat=istatus)
  CALL check_alloc_status(istatus, "hdf_read_base:dta3d")

  ALLOCATE (tem1(MAX(nx,nxdta),MAX(ny,nydta),nztmp),stat=istatus)
  CALL check_alloc_status(istatus, "hdf_read_base:tem1")

!-----------------------------------------------------------------------
!
!  Read in base state fields
!
!-----------------------------------------------------------------------

  CALL hdfread1d(IREAD,ISPLIT,ncompressx,ncompressy,fileid,"x",1,       &
                 nx,x,nxdta,dta3d,tem1,istatus)

  CALL hdfread1d(IREAD,ISPLIT,ncompressx,ncompressy,fileid,"y",2,       &
                 ny,y,nydta,dta3d,tem1,istatus)

  CALL hdfread1d(IREAD,ISPLIT,ncompressx,ncompressy,fileid,"z",3,       &
                 nz,z,nz,dta3d,tem1,istatus)

  CALL hdfread3d(IREAD,ISPLIT,ncompressx,ncompressy,fileid,"zp",        &
                 nx,ny,nz,zp,                                           &
                 nxdta,nydta,dta3d,tem1,itmp,hmax,hmin,istatus)

  CALL hdfread3d(IREAD,ISPLIT,ncompressx,ncompressy,fileid,"zpsoil",    &
                 nx,ny,nzsoil,zpsoil,                                   &
                 nxdta,nydta,dta3d,tem1,itmp,hmax,hmin,istatus)

!-----------------------------------------------------------------------
!
!  Friendly exit message
!
!-----------------------------------------------------------------------

  ireturn = 0

  DEALLOCATE (itmp,hmax,hmin, stat=istatus)
  DEALLOCATE (dta3d, tem1)

  RETURN
END SUBROUTINE hdf_read_grid

!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE HDF_READ_STATIC            ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################

SUBROUTINE hdf_read_static(IAMROOT, IREAD, ISPLIT, ncompressx, ncompressy,&
                         fileid,nx,ny,nstyps,nxdta,nydta,               &
                         landin,soiltyp,stypfrct,vegtyp,lai,roufns,veg,ireturn)

!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Read in surface static fields only.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yunheng Wang
!  2011/06/02
!  Based on subroutine hdfread.
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nstyps   Number of soil types
!
!  OUTPUT:
!
!    soiltyp  Soil type
!    vegtyp   Vegetation type
!    lai      Leaf Area Index
!    roufns   Surface roughness
!    veg      Vegetation fraction
!
!    ireturn  Return status indicator
!
!
!-----------------------------------------------------------------------

  USE arps_precision
  IMPLICIT NONE

!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------

  LOGICAL, INTENT(IN) :: IAMROOT, IREAD, ISPLIT
  INTEGER, INTENT(IN) :: ncompressx, ncompressy
  INTEGER, INTENT(IN) :: fileid(ncompressx,ncompressy)

  INTEGER, INTENT(IN) :: nx,ny,nstyps
  INTEGER, INTENT(IN) :: nxdta, nydta

  INTEGER, INTENT(OUT) :: landin
  INTEGER :: soiltyp (nx,ny,nstyps)    ! Soil type
  REAL(P) :: stypfrct(nx,ny,nstyps)    ! Fraction of soil types
  INTEGER :: vegtyp(nx,ny)             ! Vegetation type
  REAL(P) :: lai    (nx,ny)        ! Leaf Area Index
  REAL(P) :: roufns (nx,ny)        ! Surface roughness
  REAL(P) :: veg    (nx,ny)        ! Vegetation fraction

  INTEGER, INTENT(OUT) :: ireturn

!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------


  INTEGER :: i,j,k
  INTEGER :: ia,ja
  INTEGER :: ii, jj

  INTEGER :: istatus

  REAL(SP), ALLOCATABLE :: dta3d(:,:,:)
  INTEGER,  ALLOCATABLE :: dta3di(:,:,:),item1(:,:,:)
  REAL(P),  ALLOCATABLE :: tem1(:,:,:)

  INTEGER(INT16), ALLOCATABLE :: itmp(:,:,:)              ! Temporary array
  REAL(SP),       ALLOCATABLE :: hmax(:), hmin(:)         ! Temporary array


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  ALLOCATE (itmp(nxdta,nydta,nstyps),stat=istatus)
  CALL check_alloc_status(istatus, "hdf_read_static:itmp")

  ALLOCATE (hmax(nstyps),stat=istatus)
  CALL check_alloc_status(istatus, "hdf_read_static:hmax")

  ALLOCATE (hmin(nstyps),stat=istatus)
  CALL check_alloc_status(istatus, "hdf_read_static:hmin")

  ALLOCATE (dta3d(nxdta,nydta,nstyps),stat=istatus)
  CALL check_alloc_status(istatus, "hdf_read_static:dta3d")

  ALLOCATE (dta3di(nxdta,nydta,nstyps),stat=istatus)
  CALL check_alloc_status(istatus, "hdf_read_static:dta3d")

  ALLOCATE (tem1(MAX(nx,nxdta),MAX(ny,nydta),nstyps),stat=istatus)
  CALL check_alloc_status(istatus, "hdf_read_static:tem1")

  ALLOCATE (item1(MAX(nx,nxdta),MAX(ny,nydta),nstyps),stat=istatus)
  CALL check_alloc_status(istatus, "hdf_read_static:item1")

!-----------------------------------------------------------------------
!
! Check flag
!
!-----------------------------------------------------------------------

  IF (IREAD) THEN
    CALL hdfrdi(fileid(1,1),"landflg",landin,istatus)
  END IF
  IF (ISPLIT) THEN
    CALL mpupdatei(landin,1)
  END IF

  IF (landin /= 1) THEN
    WRITE(*,'(1x,a,I2)') 'ERROR: Expect landflg = 1, but find landin = ',landin
    istatus = -1
    RETURN
    !CALL arpsstop('ERROR wrong flag landflg',1)
  END IF

!-----------------------------------------------------------------------
!
!  Read in base state fields
!
!-----------------------------------------------------------------------

  CALL hdfread3di(IREAD,ISPLIT,ncompressx,ncompressy,fileid,"soiltyp",  &
                 nx,ny,nstyps,soiltyp,                                  &
                 nxdta,nydta,dta3di,item1,istatus)

  CALL hdfread3d(IREAD,ISPLIT,ncompressx,ncompressy,fileid,"stypfrct",  &
                 nx,ny,nstyps,stypfrct,                                 &
                 nxdta,nydta,dta3d,tem1,itmp,hmax,hmin,istatus)

  CALL hdfread3di(IREAD,ISPLIT,ncompressx,ncompressy,fileid,"vegtyp",   &
                 nx,ny,1,vegtyp,                                        &
                 nxdta,nydta,dta3di,item1,istatus)

  CALL hdfread2d(IREAD,ISPLIT,ncompressx,ncompressy,fileid,"lai",       &
                 nx,ny,lai,                                             &
                 nxdta,nydta,dta3d,tem1,itmp,istatus)

  CALL hdfread2d(IREAD,ISPLIT,ncompressx,ncompressy,fileid,"roufns",    &
                 nx,ny,roufns,                                          &
                 nxdta,nydta,dta3d,tem1,itmp,istatus)

  CALL hdfread2d(IREAD,ISPLIT,ncompressx,ncompressy,fileid,"veg",       &
                 nx,ny,veg,                                             &
                 nxdta,nydta,dta3d,tem1,itmp,istatus)

!-----------------------------------------------------------------------
!
!  Friendly exit message
!
!-----------------------------------------------------------------------

  ireturn = 0

  DEALLOCATE (itmp,hmax,hmin, stat=istatus)
  DEALLOCATE (dta3d, dta3di, tem1, item1)

  RETURN
END SUBROUTINE hdf_read_static

!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE HDF_READ_OPTIONAL          ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################

SUBROUTINE hdf_read_optional(IAMROOT, IREAD, ISPLIT, ncompressx, ncompressy,&
                         fileid,nx,ny,nz,nxdta,nydta,                   &
                         radfrc,radsw,rnflx,radswnet,radlwin,           &
                         usflx,vsflx,ptsflx,qvsflx,ireturn)

!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Read in radiation arrays and flux arrays.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yunheng Wang
!  2011/06/02
!  Based on subroutine hdfread.
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
!!
!  OUTPUT:
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
!
!
!-----------------------------------------------------------------------

  USE arps_precision
  IMPLICIT NONE

!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------

  LOGICAL, INTENT(IN) :: IAMROOT, IREAD, ISPLIT
  INTEGER, INTENT(IN) :: ncompressx, ncompressy
  INTEGER, INTENT(IN) :: fileid(ncompressx,ncompressy)

  INTEGER, INTENT(IN) :: nx,ny,nz
  INTEGER, INTENT(IN) :: nxdta, nydta

  REAL(P), INTENT(OUT) :: radfrc(nx,ny,nz)  ! Radiation forcing (K/s)
  REAL(P), INTENT(OUT) :: radsw (nx,ny)     ! Solar radiation reaching the surface
  REAL(P), INTENT(OUT) :: rnflx (nx,ny)     ! Net radiation flux absorbed by surface
  REAL(P), INTENT(OUT) :: radswnet(nx,ny)   ! Net shortwave radiation
  REAL(P), INTENT(OUT) :: radlwin(nx,ny)    ! Incoming longwave radiation

  REAL(P), INTENT(OUT) :: usflx (nx,ny)     ! Surface flux of u-momentum (kg/(m*s**2))
  REAL(P), INTENT(OUT) :: vsflx (nx,ny)     ! Surface flux of v-momentum (kg/(m*s**2))
  REAL(P), INTENT(OUT) :: ptsflx(nx,ny)     ! Surface heat flux (K*kg/(m**2*s))
  REAL(P), INTENT(OUT) :: qvsflx(nx,ny)     ! Surface moisture flux (kg/(m**2*s))

  INTEGER, INTENT(OUT) :: ireturn

!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------


  INTEGER :: i,j,k
  INTEGER :: ia,ja
  INTEGER :: ii, jj

  INTEGER :: istatus
  INTEGER :: radin, flxin

  REAL(SP), ALLOCATABLE :: dta3d(:,:,:)
  REAL(P),  ALLOCATABLE :: tem1(:,:,:)

  INTEGER(INT16), ALLOCATABLE :: itmp(:,:,:)              ! Temporary array
  REAL(SP),       ALLOCATABLE :: hmax(:), hmin(:)         ! Temporary array


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  ALLOCATE (itmp(nxdta,nydta,nz),stat=istatus)
  CALL check_alloc_status(istatus, "hdf_read_static:itmp")

  ALLOCATE (hmax(nz),stat=istatus)
  CALL check_alloc_status(istatus, "hdf_read_static:hmax")

  ALLOCATE (hmin(nz),stat=istatus)
  CALL check_alloc_status(istatus, "hdf_read_static:hmin")

  ALLOCATE (dta3d(nxdta,nydta,nz),stat=istatus)
  CALL check_alloc_status(istatus, "hdf_read_static:dta3d")

  ALLOCATE (tem1(MAX(nx,nxdta),MAX(ny,nydta),nz),stat=istatus)
  CALL check_alloc_status(istatus, "hdf_read_static:tem1")

!-----------------------------------------------------------------------
!
! Check flags
!
!-----------------------------------------------------------------------

  IF (IREAD) THEN
    CALL hdfrdi(fileid(1,1),"radflg",radin,istatus)
    CALL hdfrdi(fileid(1,1),"flxflg",flxin,istatus)
  END IF
  IF (ISPLIT) THEN
    CALL mpupdatei(radin,1)
    CALL mpupdatei(flxin,1)
  END IF

  IF (radin /= 1 .AND. flxin /= 1) THEN
    WRITE(*,'(1x,2(a,I2))') 'ERROR: Expect radflg & flxflg, but find radin = ',radin,', and flxin = ',flxin
    istatus = -1
    RETURN
    !CALL arpsstop('ERROR wrong flag radflg & flxin',1)
  END IF

!-----------------------------------------------------------------------
!
!  Read in base state fields
!
!-----------------------------------------------------------------------

  IF (radin == 1) THEN
    CALL hdfread3d(IREAD,ISPLIT,ncompressx,ncompressy,fileid,"radfrc",  &
                   nx,ny,nz,radfrc,                                     &
                   nxdta,nydta,dta3d,tem1,itmp,hmax,hmin,istatus)

    CALL hdfread2d(IREAD,ISPLIT,ncompressx,ncompressy,fileid,"radsw",   &
                   nx,ny,radsw,                                         &
                   nxdta,nydta,dta3d,tem1,itmp,istatus)

    CALL hdfread2d(IREAD,ISPLIT,ncompressx,ncompressy,fileid,"rnflx",   &
                   nx,ny,rnflx,                                         &
                   nxdta,nydta,dta3d,tem1,itmp,istatus)

    CALL hdfread2d(IREAD,ISPLIT,ncompressx,ncompressy,fileid,"radswnet",&
                   nx,ny,radswnet,                                      &
                   nxdta,nydta,dta3d,tem1,itmp,istatus)

    CALL hdfread2d(IREAD,ISPLIT,ncompressx,ncompressy,fileid,"radlwin", &
                   nx,ny,radlwin,                                       &
                   nxdta,nydta,dta3d,tem1,itmp,istatus)

  ELSE
    WRITE(*,'(1x,a)') 'WARNING: Radiation variables are not available in the data file.'
  END IF

  IF (flxin == 1) THEN
    CALL hdfread2d(IREAD,ISPLIT,ncompressx,ncompressy,fileid,"usflx",   &
                   nx,ny,usflx,                                         &
                   nxdta,nydta,dta3d,tem1,itmp,istatus)

    CALL hdfread2d(IREAD,ISPLIT,ncompressx,ncompressy,fileid,"vsflx",   &
                   nx,ny,vsflx,                                         &
                   nxdta,nydta,dta3d,tem1,itmp,istatus)

    CALL hdfread2d(IREAD,ISPLIT,ncompressx,ncompressy,fileid,"ptsflx",  &
                   nx,ny,ptsflx,                                        &
                   nxdta,nydta,dta3d,tem1,itmp,istatus)

    CALL hdfread2d(IREAD,ISPLIT,ncompressx,ncompressy,fileid,"qvsflx",  &
                   nx,ny,qvsflx,                                        &
                   nxdta,nydta,dta3d,tem1,itmp,istatus)

  ELSE
    WRITE(*,'(1x,a)') 'WARNING: Flux variables are not available in the data file.'
  END IF

!-----------------------------------------------------------------------
!
!  Friendly exit message
!
!-----------------------------------------------------------------------

  ireturn = 0

  DEALLOCATE (itmp,hmax,hmin, stat=istatus)
  DEALLOCATE (dta3d, tem1)

  RETURN
END SUBROUTINE hdf_read_optional

!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE HDF_READ_VAR               ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################

SUBROUTINE hdf_read_var(IAMROOT, IREAD, ISPLIT, ncompressx, ncompressy, &
                        fileid,nx,ny,nxdta,nydta,                       &
                        varname,varrank,vardim,varjoin,varout,varsize,  &
                        ireturn)

!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Read in an arbitrary arrays.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yunheng Wang
!  2012/09/21
!  Based on subroutine hdfread.
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!
!  OUTPUT:
!
!    varout   Array to be read
!
!    ireturn  Return status indicator
!
!-----------------------------------------------------------------------

  USE arps_precision
  IMPLICIT NONE

!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------

  LOGICAL, INTENT(IN) :: IAMROOT, IREAD, ISPLIT
  INTEGER, INTENT(IN) :: ncompressx, ncompressy
  INTEGER, INTENT(IN) :: fileid(ncompressx,ncompressy)

  INTEGER, INTENT(IN) :: nx,ny
  INTEGER, INTENT(IN) :: nxdta, nydta

  CHARACTER(LEN=*), INTENT(IN) :: varname
  INTEGER,          INTENT(IN) :: varrank
  INTEGER,          INTENT(IN) :: vardim(varrank)
  INTEGER,          INTENT(IN) :: varjoin  ! 0 : no join
                                           ! 1 : Join in X direction
                                           ! 2 : Join in Y direction
                                           ! 3 : Join in Bother X and Y direction for the first 2 dimensions
  INTEGER,          INTENT(IN) :: varsize

  REAL(P), INTENT(OUT) :: varout  (varsize)

  INTEGER, INTENT(OUT) :: ireturn

!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------

  INTEGER :: i,j,k
  INTEGER :: ia,ja
  INTEGER :: ii, jj

  INTEGER :: n1in, n2in, n3in, n4in
  INTEGER :: n1lo, n2lo

  INTEGER :: istatus

  REAL(SP), ALLOCATABLE :: dtain(:,:,:,:)
  REAL(P),  ALLOCATABLE :: tem1(:,:,:,:)

  INTEGER(INT16), ALLOCATABLE :: itmp(:,:,:,:)              ! Temporary array
  REAL(SP),       ALLOCATABLE :: hmax(:), hmin(:)         ! Temporary array


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  n1in = 1;  n1lo = 1
  n2in = 1;  n2lo = 1
  n3in = 1
  n4in = 1
  SELECT CASE (varrank)
  CASE (1)

    IF (varjoin == 1) THEN
      n1in = nxdta
    ELSE IF (varjoin == 2) THEN
      n1in = nydta
    ELSE
      n1in = vardim(1)
    END IF

    n1lo = vardim(1)     ! local working array size

  CASE (2,3,4)

    IF (varjoin == 3) THEN

      n1in = nxdta       ! data size in the file
      n2in = nydta

      n1lo = MAX(vardim(1),n1in)
      n2lo = MAX(vardim(2),n2in)

    ELSE
      WRITE(*,'(1x,3a,I0,a)') 'ERROR: unsupport variable ',TRIM(varname), &
                              ' with varjoin = ',varjoin,'.'
      CALL arpsstop('ERROR: wrong argument in hdf_read_var.',1)
    END IF

  CASE default

    WRITE(*,'(1x,a,I0,a)') 'ERROR: unsupport variable rank varrank = ',varrank,'.'
    CALL arpsstop('ERROR: wrong argument in hdf_read_var.',1)

  END SELECT

  IF (varrank >= 3) n3in = vardim(3)
  IF (varrank >= 4) n4in = vardim(4)

  ALLOCATE (tem1(n1lo,n2lo,n3in,n4in),stat=istatus)
  CALL check_alloc_status(istatus, "hdf_read_var:tem1")

  ALLOCATE (dtain(n1in,n2in,n3in,n4in),stat=istatus)
  CALL check_alloc_status(istatus, "hdf_read_var:dta3d")

  IF (varrank >= 2) THEN
    ALLOCATE (itmp(n1in,n2in,n3in,n4in),stat=istatus)
    CALL check_alloc_status(istatus, "hdf_read_var:itmp")
  END IF

  IF (varrank >= 3) THEN
    ALLOCATE (hmax(n3in),stat=istatus)
    CALL check_alloc_status(istatus, "hdf_read_var:hmax")

    ALLOCATE (hmin(n3in),stat=istatus)
    CALL check_alloc_status(istatus, "hdf_read_var:hmin")
  END IF

!-----------------------------------------------------------------------
!
!  Read in the fields
!
!-----------------------------------------------------------------------

  SELECT CASE (varrank)
  CASE (1)

    CALL hdfread1d(IREAD,ISPLIT,ncompressx,ncompressy,fileid,varname,varjoin,&
                   vardim(1),varout,n1in,dtain,tem1,istatus)

  CASE (2)
    CALL hdfread2d(IREAD,ISPLIT,ncompressx,ncompressy,fileid,varname,   &
                   vardim(1),vardim(2),varout,                          &
                   n1in,n2in,dtain,tem1,itmp,istatus)

  CASE (3)
    CALL hdfread3d(IREAD,ISPLIT,ncompressx,ncompressy,fileid,varname,   &
                   vardim(1),vardim(2),vardim(3),varout,                &
                   n1in,n2in,dtain,tem1,itmp,hmax,hmin,istatus)

  CASE (4)

    CALL hdfread4d(IREAD,ISPLIT,ncompressx,ncompressy,fileid,varname,   &
                  vardim(1),vardim(2),vardim(3),vardim(4),varout,       &
                  n1in,n2in,dtain,tem1,itmp, hmax, hmin, istatus)

  END SELECT

!-----------------------------------------------------------------------
!
!  Friendly exit message
!
!-----------------------------------------------------------------------

  ireturn = 0

  IF ( ALLOCATED(itmp) ) DEALLOCATE ( itmp )
  IF ( ALLOCATED(hmax) ) DEALLOCATE ( hmax, hmin )
  DEALLOCATE (dtain, tem1)

  RETURN
END SUBROUTINE hdf_read_var

!#######################################################################

SUBROUTINE hdfread1d(IREAD,ISPLIT,ncompressx,ncompressy,fileid,fieldname,&
                     stagger,nout,var1d,nin,dta1d,tem1,istatus)

!-----------------------------------------------------------------------
!  Read one 1d field from the HDF file and determins whether
!  it should be joined or split.
!-----------------------------------------------------------------------
  use arps_precision

  IMPLICIT NONE
  LOGICAL, INTENT(IN) :: IREAD, ISPLIT
  INTEGER, INTENT(IN) :: ncompressx, ncompressy
  INTEGER, INTENT(IN) :: fileid(ncompressx,ncompressy)
  CHARACTER(LEN=*), INTENT(IN) :: fieldname
  INTEGER, INTENT(IN) :: nout
  INTEGER, INTENT(IN) :: nin
  INTEGER, INTENT(IN) :: stagger

  REAL(P),  INTENT(OUT) :: var1d(nout)

  ! Work arrays
  REAL(SP), INTENT(OUT) :: dta1d(nin)
  REAL(P),  INTENT(OUT) :: tem1(MAX(nout,nin))

  INTEGER, INTENT(OUT) :: istatus

!-----------------------------------------------------------------------

  INTEGER :: i,j
  INTEGER :: ii,jj
  INTEGER :: ia,ja

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  istatus = 0

  IF (IREAD) THEN
    IF (stagger == 1) THEN       ! X direction join

      DO ii = 1, ncompressx
        CALL hdfrd1d(fileid(ii,1),TRIM(fieldname),nin,dta1d,istatus)
        IF (istatus /= 0) THEN
          WRITE(*,'(1x,a,I3)') 'WARNING: Field '//TRIM(fieldname)//' error with hdfrd1d, istat = ',istatus
          EXIT
        END IF

        DO i = 1, nin
          ia = (ii-1)*(nin-3)+i
          tem1(ia)  = dta1d(i)   ! join and convert precision
        END DO
      END DO

    ELSE IF (stagger == 2) THEN  ! Y direction join

      DO jj = 1, ncompressy
        CALL hdfrd1d(fileid(1,jj),TRIM(fieldname),nin,dta1d,istatus)
        IF (istatus /= 0) THEN
          WRITE(*,'(1x,a,I3)') 'WARNING: Field '//TRIM(fieldname)//' error with hdfrd1d, istat = ',istatus
          EXIT
        END IF

        DO j = 1, nin
          ja = (jj-1)*(nin-3)+j
          tem1(ja)  = dta1d(j)   ! join and convert precision
        END DO
      END DO

    ELSE IF (stagger == 3 .OR. stagger == 0) THEN  ! Z direction, no join

      CALL hdfrd1d(fileid(1,1),TRIM(fieldname),nin,dta1d,istatus)
      IF (istatus /= 0) THEN
        WRITE(*,'(1x,a,I3)') 'WARNING: Field '//TRIM(fieldname)//' error with hdfrd1d, istat = ',istatus
        RETURN
      END IF

      DO j = 1, nin
        tem1(j)  = dta1d(j)   ! convert precision
      END DO

    END IF
  END IF

  IF (ISPLIT) THEN
    CALL mpupdatei(istatus,1)
    IF (istatus /= 0) RETURN

    IF (stagger == 1) THEN
      CALL mpisplit1dx(tem1,nout,var1d)
    ELSE IF (stagger == 2) THEN
      CALL mpisplit1dy(tem1,nout,var1d)
    ELSE IF (stagger == 3 .OR. stagger == 0) THEN
      CALL mpupdater(tem1,nout)
      var1d(:) = tem1(:)
    END IF
  ELSE
    IF (istatus /= 0) RETURN
    var1d(:) = tem1(:)
  END IF

  RETURN
END SUBROUTINE hdfread1d

!#######################################################################

SUBROUTINE hdfread2d(IREAD,ISPLIT,ncompressx,ncompressy,fileid,fieldname,&
                     nx,ny,var2d,                                       &
                     nxdta,nydta,dta2d,tem1,itmp,istatus)
!-----------------------------------------------------------------------
!  Read one 3d field from the HDF file and determins whether
!  it should be joined or split.
!-----------------------------------------------------------------------
  use arps_precision

  IMPLICIT NONE
  LOGICAL, INTENT(IN) :: IREAD, ISPLIT
  INTEGER, INTENT(IN) :: ncompressx, ncompressy
  INTEGER, INTENT(IN) :: fileid(ncompressx,ncompressy)
  CHARACTER(LEN=*), INTENT(IN) :: fieldname
  INTEGER, INTENT(IN) :: nx, ny
  INTEGER, INTENT(IN) :: nxdta, nydta

  REAL(P),  INTENT(OUT) :: var2d(nx,ny)

  ! Work arrays
  REAL(SP), INTENT(OUT) :: dta2d(nxdta,nydta)
  REAL(P),  INTENT(OUT) :: tem1(MAX(nx,nxdta),MAX(ny,nydta))

  INTEGER(INT16), INTENT(OUT) :: itmp(nxdta,nydta)

  INTEGER, INTENT(OUT) :: istatus

!-----------------------------------------------------------------------

  INTEGER :: i,j
  INTEGER :: ii,jj
  INTEGER :: ia,ja

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  istatus = 0

  IF (IREAD) THEN
    DO jj = 1, ncompressy
      DO ii = 1, ncompressx
        CALL hdfrd2d(fileid(ii,jj),TRIM(fieldname),nxdta,nydta,dta2d,   &
                     istatus,itmp)
        IF (istatus /= 0) THEN
          WRITE(*,'(1x,a,I3)') 'WARNING: Field '//TRIM(fieldname)//' error with hdfrd2d, istat = ',istatus
          EXIT
        END IF

        DO j = 1, nydta
          ja = (jj-1)*(nydta-3)+j
          DO i = 1, nxdta
            ia = (ii-1)*(nxdta-3)+i
            tem1(ia,ja)  = dta2d(i,j)   ! join and convert precision
          END DO
        END DO

      END DO
    END DO
  END IF
  IF (ISPLIT) THEN
    CALL mpupdatei(istatus,1)
    IF (istatus /= 0) RETURN
    CALL mpisplit3d(tem1,nx,ny,1,var2d)
  ELSE
    IF (istatus /= 0) RETURN
    var2d(:,:) = tem1(:,:)
  END IF

  RETURN
END SUBROUTINE hdfread2d

!#######################################################################

SUBROUTINE hdfread3d(IREAD,ISPLIT,ncompressx,ncompressy,fileid,fieldname,&
                     nx,ny,nz,var3d,                                    &
                     nxdta,nydta,dta3d,tem1,itmp,hmax,hmin,istatus)
!-----------------------------------------------------------------------
!  Read one 3d field from the HDF file and determins whether
!  it should be joined or split.
!-----------------------------------------------------------------------
  use arps_precision

  IMPLICIT NONE
  LOGICAL, INTENT(IN) :: IREAD, ISPLIT
  INTEGER, INTENT(IN) :: ncompressx, ncompressy
  INTEGER, INTENT(IN) :: fileid(ncompressx,ncompressy)
  CHARACTER(LEN=*), INTENT(IN) :: fieldname
  INTEGER, INTENT(IN) :: nx, ny, nz
  INTEGER, INTENT(IN) :: nxdta, nydta

  REAL(P),  INTENT(OUT) :: var3d(nx,ny,nz)

  ! Work arrays
  REAL(SP), INTENT(OUT) :: dta3d(nxdta,nydta,nz)
  REAL(P),  INTENT(OUT) :: tem1(MAX(nx,nxdta),MAX(ny,nydta),nz)

  INTEGER(INT16), INTENT(OUT) :: itmp(nxdta,nydta,nz)
  REAL(SP),       INTENT(OUT) :: hmax(nz), hmin(nz)

  INTEGER, INTENT(OUT) :: istatus

!-----------------------------------------------------------------------

  INTEGER :: i,j,k
  INTEGER :: ii,jj
  INTEGER :: ia,ja

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  istatus = 0

  IF (IREAD) THEN
    DO jj = 1, ncompressy
      DO ii = 1, ncompressx
        CALL hdfrd3d(fileid(ii,jj),TRIM(fieldname),nxdta,nydta,nz,dta3d,&
                     istatus,itmp,hmax,hmin)
        IF (istatus /= 0) THEN
          WRITE(*,'(1x,a,I3)') 'WARNING: Field '//TRIM(fieldname)//' error with hdfrd3d, istat = ',istatus
          EXIT
        END IF

        DO k = 1, nz
          DO j = 1, nydta
            ja = (jj-1)*(nydta-3)+j
            DO i = 1, nxdta
              ia = (ii-1)*(nxdta-3)+i
              tem1(ia,ja,k)  = dta3d(i,j,k)   ! join and convert precision
            END DO
          END DO
        END DO

      END DO
    END DO
  END IF
  IF (ISPLIT) THEN
    CALL mpupdatei(istatus,1)
    IF (istatus /= 0) RETURN
    CALL mpisplit3d(tem1,nx,ny,nz,var3d)
  ELSE
    IF (istatus /= 0) RETURN
    var3d(:,:,:) = tem1(:,:,:)
  END IF

  RETURN
END SUBROUTINE hdfread3d

!#######################################################################

SUBROUTINE hdfread3di(IREAD,ISPLIT,ncompressx,ncompressy,fileid,fieldname,&
                     nx,ny,nz,var3d,                                    &
                     nxdta,nydta,dta3d,tem1,istatus)
!-----------------------------------------------------------------------
!  Read one 3d field from the HDF file and determins whether
!  it should be joined or split.
!-----------------------------------------------------------------------

  IMPLICIT NONE
  LOGICAL, INTENT(IN) :: IREAD, ISPLIT
  INTEGER, INTENT(IN) :: ncompressx, ncompressy
  INTEGER, INTENT(IN) :: fileid(ncompressx,ncompressy)
  CHARACTER(LEN=*), INTENT(IN) :: fieldname
  INTEGER, INTENT(IN) :: nx, ny, nz
  INTEGER, INTENT(IN) :: nxdta, nydta

  INTEGER,  INTENT(OUT) :: var3d(nx,ny,nz)

  ! Work arrays
  INTEGER, INTENT(OUT) :: dta3d(nxdta,nydta,nz)
  INTEGER, INTENT(OUT) :: tem1(MAX(nx,nxdta),MAX(ny,nydta),nz)

  INTEGER, INTENT(OUT) :: istatus

!-----------------------------------------------------------------------

  INTEGER :: i,j,k
  INTEGER :: ii,jj
  INTEGER :: ia,ja

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  istatus = 0

  IF (IREAD) THEN
    DO jj = 1, ncompressy
      DO ii = 1, ncompressx
        CALL hdfrd3di(fileid(ii,jj),TRIM(fieldname),nxdta,nydta,nz,dta3d,istatus)
        IF (istatus /= 0) THEN
          WRITE(*,'(1x,a,I3)') 'WARNING: Field '//TRIM(fieldname)//' error with hdfrd3di, istat = ',istatus
          EXIT
        END IF

        DO k = 1, nz
          DO j = 1, nydta
            ja = (jj-1)*(nydta-3)+j
            DO i = 1, nxdta
              ia = (ii-1)*(nxdta-3)+i
              tem1(ia,ja,k)  = dta3d(i,j,k)   ! join and convert precision
            END DO
          END DO
        END DO

      END DO
    END DO
  END IF
  IF (ISPLIT) THEN
    CALL mpupdatei(istatus,1)
    IF (istatus /= 0) RETURN
    CALL mpisplit3di(tem1,nx,ny,nz,var3d)
  ELSE
    IF (istatus /= 0) RETURN
    var3d(:,:,:) = tem1(:,:,:)
  END IF

  RETURN
END SUBROUTINE hdfread3di

!#######################################################################

SUBROUTINE hdfread4d(IREAD,ISPLIT,ncompressx,ncompressy,fileid,fieldname,&
                     nx,ny,nzsoil,nstyp,var4d,                          &
                     nxdta,nydta,dta4d,tem1,itmp,hmax,hmin,istatus)
!-----------------------------------------------------------------------
!  Read one 3d field from the HDF file and determins whether
!  it should be joined or split.
!-----------------------------------------------------------------------
  use arps_precision

  IMPLICIT NONE
  LOGICAL, INTENT(IN) :: IREAD, ISPLIT
  INTEGER, INTENT(IN) :: ncompressx, ncompressy
  INTEGER, INTENT(IN) :: fileid(ncompressx,ncompressy)
  CHARACTER(LEN=*), INTENT(IN) :: fieldname
  INTEGER, INTENT(IN) :: nx, ny, nzsoil,nstyp
  INTEGER, INTENT(IN) :: nxdta, nydta

  REAL(P),  INTENT(OUT) :: var4d(nx,ny,nzsoil,nstyp)

  ! Work arrays
  REAL(SP), INTENT(OUT) :: dta4d(nxdta,nydta,nzsoil,nstyp)
  REAL(P),  INTENT(OUT) :: tem1(MAX(nx,nxdta),MAX(ny,nydta),nzsoil,nstyp)

  INTEGER(INT16), INTENT(OUT) :: itmp(nxdta,nydta,nzsoil,nstyp)
  REAL(SP),       INTENT(OUT) :: hmax(nzsoil), hmin(nzsoil)

  INTEGER, INTENT(OUT) :: istatus

!-----------------------------------------------------------------------

  INTEGER :: i,j,k,is
  INTEGER :: ii,jj
  INTEGER :: ia,ja

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  istatus = 0

  IF (IREAD) THEN
    DO jj = 1, ncompressy
      DO ii = 1, ncompressx
        CALL hdfrd4d(fileid(ii,jj),TRIM(fieldname),nxdta,nydta,nzsoil,nstyp,dta4d,&
                     istatus,itmp,hmax,hmin)
        IF (istatus /= 0) THEN
          WRITE(*,'(1x,a,I3)') 'WARNING: Field '//TRIM(fieldname)//' error with hdfrd4d, istat = ',istatus
          EXIT
        END IF

        DO is = 1,nstyp
          DO k = 1, nzsoil
            DO j = 1, nydta
              ja = (jj-1)*(nydta-3)+j
              DO i = 1, nxdta
                ia = (ii-1)*(nxdta-3)+i
                tem1(ia,ja,k,is)  = dta4d(i,j,k,is)   ! join and convert precision
              END DO
            END DO
          END DO
        END DO

      END DO
    END DO
  END IF
  IF (ISPLIT) THEN
    CALL mpupdatei(istatus,1)
    IF (istatus /= 0) RETURN
    CALL mpisplit4d(tem1,nx,ny,nzsoil,nstyp,var4d)
  ELSE
    IF (istatus /= 0) RETURN
    var4d(:,:,:,:) = tem1(:,:,:,:)
  END IF

  RETURN
END SUBROUTINE hdfread4d

SUBROUTINE hdf_read_dsd(IAMROOT, IREAD, ISPLIT, ncompressx, ncompressy, &
                        fileid,                                         &
                        ntcloud,n0rain, n0snow,  n0grpl,   n0hail,      &
                        rhoice,  rhosnow,  rhogrpl,  rhohail,           &
                        alpharain,alphaice,alphasnow,alphagrpl,alphahail, &
                        ireturn)
  USE arps_precision
  IMPLICIT NONE

!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------

  LOGICAL, INTENT(IN) :: IAMROOT, IREAD, ISPLIT
  INTEGER, INTENT(IN) :: ncompressx, ncompressy
  INTEGER, INTENT(IN) :: fileid(ncompressx,ncompressy)

  REAL(P), INTENT(OUT) :: ntcloud
  REAL(P), INTENT(OUT) :: n0rain
  REAL(P), INTENT(OUT) :: n0snow
  REAL(P), INTENT(OUT) :: n0grpl
  REAL(P), INTENT(OUT) :: n0hail
  REAL(P), INTENT(OUT) :: rhoice
  REAL(P), INTENT(OUT) :: rhosnow
  REAL(P), INTENT(OUT) :: rhogrpl
  REAL(P), INTENT(OUT) :: rhohail
  REAL(P), INTENT(OUT) :: alpharain
  REAL(P), INTENT(OUT) :: alphaice
  REAL(P), INTENT(OUT) :: alphasnow
  REAL(P), INTENT(OUT) :: alphagrpl
  REAL(P), INTENT(OUT) :: alphahail

  INTEGER, INTENT(OUT) :: ireturn

!-----------------------------------------------------------------------

  INTEGER :: istat

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  ireturn = 0

  IF (IREAD) THEN
    CALL hdfrdr(fileid(1,1),"n0rain", n0rain,istat)
    IF (istat /= 0)  WRITE(*,'(1x,a)') 'WARNING: n0rain is not available in the data file.'
    CALL hdfrdr(fileid(1,1),"n0snow", n0snow,istat)
    IF (istat /= 0)  WRITE(*,'(1x,a)') 'WARNING: n0snow is not available in the data file.'
    CALL hdfrdr(fileid(1,1),"n0hail", n0hail,istat)
    IF (istat /= 0)  WRITE(*,'(1x,a)') 'WARNING: n0hail is not available in the data file.'
    CALL hdfrdr(fileid(1,1),"rhosnow",rhosnow,istat)
    IF (istat /= 0)  WRITE(*,'(1x,a)') 'WARNING: rhosnow is not available in the data file.'
    CALL hdfrdr(fileid(1,1),"rhohail",rhohail,istat)
    IF (istat /= 0)  WRITE(*,'(1x,a)') 'WARNING: rhohail is not available in the data file.'

    CALL hdfrdr(fileid(1,1),'ntcloud',ntcloud,istat)
    IF (istat /= 0)  WRITE(*,'(1x,a)') 'WARNING: ntcloud is not available in the data file.'
    CALL hdfrdr(fileid(1,1),'n0grpl', n0grpl,istat)
    IF (istat /= 0)  WRITE(*,'(1x,a)') 'WARNING: n0grpl is not available in the data file.'
    CALL hdfrdr(fileid(1,1),'rhoice', rhoice,istat)
    IF (istat /= 0)  WRITE(*,'(1x,a)') 'WARNING: rhoice is not available in the data file.'
    CALL hdfrdr(fileid(1,1),'rhogrpl',rhogrpl,istat)
    IF (istat /= 0)  WRITE(*,'(1x,a)') 'WARNING: rhogrpl is not available in the data file.'
    CALL hdfrdr(fileid(1,1),'alpharain',alpharain,istat)
    IF (istat /= 0)  WRITE(*,'(1x,a)') 'WARNING: alpharain is not available in the data file.'
    CALL hdfrdr(fileid(1,1),'alphaice', alphaice, istat)
    IF (istat /= 0)  WRITE(*,'(1x,a)') 'WARNING: alphaice is not available in the data file.'
    CALL hdfrdr(fileid(1,1),'alphasnow',alphasnow,istat)
    IF (istat /= 0)  WRITE(*,'(1x,a)') 'WARNING: alphasnow is not available in the data file.'
    CALL hdfrdr(fileid(1,1),'alphagrpl',alphagrpl,istat)
    IF (istat /= 0)  WRITE(*,'(1x,a)') 'WARNING: alphagrpl is not available in the data file.'
    CALL hdfrdr(fileid(1,1),'alphahail',alphahail,istat)
    IF (istat /= 0)  WRITE(*,'(1x,a)') 'WARNING: alphahail is not available in the data file.'

  END IF

  IF (ISPLIT) THEN
    CALL mpupdater(ntcloud  ,1)
    CALL mpupdater(n0rain   ,1)
    CALL mpupdater(n0snow   ,1)
    CALL mpupdater(n0grpl   ,1)
    CALL mpupdater(n0hail   ,1)
    CALL mpupdater(rhoice   ,1)
    CALL mpupdater(rhosnow  ,1)
    CALL mpupdater(rhogrpl  ,1)
    CALL mpupdater(rhohail  ,1)
    CALL mpupdater(alpharain,1)
    CALL mpupdater(alphaice ,1)
    CALL mpupdater(alphasnow,1)
    CALL mpupdater(alphagrpl,1)
    CALL mpupdater(alphahail,1)
  END IF

  RETURN
END SUBROUTINE hdf_read_dsd

SUBROUTINE hdf_read_head(IAMROOT, IREAD, ISPLIT, ncompressx, ncompressy,&
                        fileid,                                         &
                        month,day,year, hour,minute,second,             &
                        umove,vmove,xgrdorg,ygrdorg,                    &
                        mapproj,trulat1,trulat2,trulon, sclfct,         &
                        tstop,thisdmp,latitud,ctrlat,ctrlon,            &
                        grdin,basin, varin,mstin,icein,trbin,totin,     &
                        sfcin, rainin,tkein, prcin, radin,flxin,snowin, &
                        ireturn)
  USE arps_precision
  IMPLICIT NONE

!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------

  LOGICAL, INTENT(IN) :: IAMROOT, IREAD, ISPLIT
  INTEGER, INTENT(IN) :: ncompressx, ncompressy
  INTEGER, INTENT(IN) :: fileid(ncompressx,ncompressy)

  INTEGER, INTENT(OUT) :: month,day,year
  INTEGER, INTENT(OUT) :: hour,minute,second

  REAL(P), INTENT(OUT) :: umove,vmove
  REAL(P), INTENT(OUT) :: xgrdorg,ygrdorg

  INTEGER, INTENT(OUT) :: mapproj
  REAL(P), INTENT(OUT) :: trulat1,trulat2,trulon
  REAL(P), INTENT(OUT) :: sclfct
  REAL(P), INTENT(OUT) :: tstop,thisdmp
  REAL(P), INTENT(OUT) :: latitud
  REAL(P), INTENT(OUT) :: ctrlat,ctrlon

  INTEGER, INTENT(OUT) :: grdin,basin, varin,mstin,icein,trbin,totin
  INTEGER, INTENT(OUT) :: sfcin, rainin,tkein, prcin, radin,flxin,snowin

  INTEGER, INTENT(OUT) :: ireturn

!-----------------------------------------------------------------------

  INTEGER :: istat

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  ireturn = 0

  IF (IREAD) THEN

    CALL hdfrdi(fileid(1,1),"month", month, istat)
    CALL hdfrdi(fileid(1,1),"day",   day,   istat)
    CALL hdfrdi(fileid(1,1),"year",  year,  istat)
    CALL hdfrdi(fileid(1,1),"hour",  hour,  istat)
    CALL hdfrdi(fileid(1,1),"minute",minute,istat)
    CALL hdfrdi(fileid(1,1),"second",second,istat)

    CALL hdfrdr(fileid(1,1),"umove",  umove,  istat)
    CALL hdfrdr(fileid(1,1),"vmove",  vmove,  istat)
    CALL hdfrdr(fileid(1,1),"xgrdorg",xgrdorg,istat)
    CALL hdfrdr(fileid(1,1),"ygrdorg",ygrdorg,istat)

    CALL hdfrdi(fileid(1,1),"mapproj",mapproj,istat)
    CALL hdfrdr(fileid(1,1),"trulat1",trulat1,istat)
    CALL hdfrdr(fileid(1,1),"trulat2",trulat2,istat)
    CALL hdfrdr(fileid(1,1),"trulon", trulon, istat)
    CALL hdfrdr(fileid(1,1),"sclfct", sclfct, istat)
    CALL hdfrdr(fileid(1,1),"tstop",  tstop,  istat)
    CALL hdfrdr(fileid(1,1),"thisdmp",thisdmp,istat)
    CALL hdfrdr(fileid(1,1),"latitud",latitud,istat)
    CALL hdfrdr(fileid(1,1),"ctrlat", ctrlat, istat)
    CALL hdfrdr(fileid(1,1),"ctrlon", ctrlon, istat)

    CALL hdfrdi(fileid(1,1),"grdflg",grdin,istat)
    CALL hdfrdi(fileid(1,1),"basflg",basin,istat)
    CALL hdfrdi(fileid(1,1),"varflg",varin,istat)
    CALL hdfrdi(fileid(1,1),"mstflg",mstin,istat)
    CALL hdfrdi(fileid(1,1),"iceflg",icein,istat)
    CALL hdfrdi(fileid(1,1),"trbflg",trbin,istat)
    CALL hdfrdi(fileid(1,1),"sfcflg",sfcin,istat)
    CALL hdfrdi(fileid(1,1),"rainflg",rainin,istat)
    CALL hdfrdi(fileid(1,1),"totflg",totin,istat)
    CALL hdfrdi(fileid(1,1),"tkeflg",tkein,istat)

    CALL hdfrdi(fileid(1,1),"prcflg",prcin,istat)
    CALL hdfrdi(fileid(1,1),"radflg",radin,istat)
    CALL hdfrdi(fileid(1,1),"flxflg",flxin,istat)
    CALL hdfrdi(fileid(1,1),"snowflg",snowin,istat)

  END IF

  IF (ISPLIT) THEN
    CALL mpupdatei(month, 1)
    CALL mpupdatei(day,   1)
    CALL mpupdatei(year,  1)
    CALL mpupdatei(hour,  1)
    CALL mpupdatei(minute,1)
    CALL mpupdatei(second,1)

    CALL mpupdater(umove,  1)
    CALL mpupdater(vmove,  1)
    CALL mpupdater(xgrdorg,1)
    CALL mpupdater(ygrdorg,1)

    CALL mpupdatei(mapproj,1)
    CALL mpupdater(trulat1,1)
    CALL mpupdater(trulat2,1)
    CALL mpupdater(trulon, 1)
    CALL mpupdater(sclfct, 1)
    CALL mpupdater(tstop,  1)
    CALL mpupdater(thisdmp,1)
    CALL mpupdater(latitud,1)
    CALL mpupdater(ctrlat, 1)
    CALL mpupdater(ctrlon, 1)

    CALL mpupdatei(grdin, 1)
    CALL mpupdatei(basin, 1)
    CALL mpupdatei(varin, 1)
    CALL mpupdatei(mstin, 1)
    CALL mpupdatei(icein, 1)
    CALL mpupdatei(trbin, 1)
    CALL mpupdatei(sfcin, 1)
    CALL mpupdatei(rainin,1)
    CALL mpupdatei(totin, 1)
    CALL mpupdatei(tkein, 1)

    CALL mpupdatei(prcin, 1)
    CALL mpupdatei(radin, 1)
    CALL mpupdatei(flxin, 1)
    CALL mpupdatei(snowin,1)

  END IF

  RETURN
END SUBROUTINE hdf_read_head
