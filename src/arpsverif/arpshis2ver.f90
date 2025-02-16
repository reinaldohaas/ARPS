!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINE HIS2VER                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE HIS2VER(nx,ny,nz,nzsoil,x,y,z,zp,zpsoil,hterain,          &
               uprt,vprt,wprt,ptprt,pprt,qvprt,                      &
               qscalar,tke,kmh,kmv,                                  &
               ubar,vbar,wbar,ptbar,pbar,rhobar,qvbar,               &
               nstyps,soiltyp,stypfrct,vegtyp,lai,roufns,            &
               veg,tsoil,qsoil,wetcanp,snowdpth,                     &
               raing,rainc,prcrate,radfrc,radsw,rnflx,               &
               radswnet,radlwin,                                     &
               usflx,vsflx,ptsflx,qvsflx,                            &
               hinfmt,fgrdbasfn,hisfile,nhisfile,                    &
               sndopt,sndlist,sndrunname,snddomlist,                 &
               proopt,prolist,prorunname,prodomlist,                 &
               sfcopt,sfclist,                                       &
               sfcobsdir,sfcpre,sfcpost,                             &
               mesoobsdir,mesopre,mesopost,                          &
               sfcrunname,blackfile,                                 &
               precveropt,preclist,precrunname,precdomlist,          &
               mosopt,moslist,mosrunname,mosdomlist,                 &
               arpsnn_opt,nnrunname,wtsdir,nnoutputfn,               &
               model_data,obsrv_data)

!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Performs extraction of sounding, profiler, and surface data
!  (similar to what is done using ARPS with veridmp.f) using only
!  the ARPS' history dump data.  The structure of this program is
!  similar to arpsextsnd.f.  In order to avoid having to make
!  alterations to the makefile, to compile this program:
!
!    1)  Copy arpsextsnd.f to arpsextsnd.f.orig (or something
!        similar) in the ARPS directory
!
!    2)  Copy arpshis2ver.f to arpsextsnd.f in the ARPS directory
!
!    3)  Copy vericst.inc to the ARPS include directory
!
!    4)  Compile with:  'makearps arpsextsnd'
!
!    5)  Move the executable 'arpsextsnd' to 'his2ver'
!
!  The executable uses only a 4-line input file (arpshis2ver.input).
!
!  NOTE:  The surface data extraction portion of the program will
!         only work if the necessary arrays are contained in the
!         history dump data.  For these arrays to be present the
!         model run must have used surface physics, terrain, etc.
!         and varout, rainout, sfcout, and landout must all be set
!         to 1.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: John Mewes
!
!  Original Coding: 09/22/97
!
!  MODIFICATION HISTORY:
!
!  09/29/97:  Added adjustable size map projection arrays to the
!             subroutine calls to snddump, prodump, and sfcdump
!             as the SUN compilers will not allow adjustable length
!             arrays to be declared within a subroutine.
!
!  06/29/98:  E. Kemp, Project COMET-Tinker
!             Modified for use in ARPS 4.3.5.  Also added seperate
!             *.prec* dumps (subroutine precdump) for STORMWAVE
!             precipitation verification.
!
!  18 May 2000:  E. Kemp
!                Updated for ARPS 4.5.1.  Added namelists for
!                input file, modified code to process multiple
!                history dumps, require user to name the .tbl
!                files (instead of hardwired), added valid time
!                and forecast information to output files, and
!                added options to verify only types of data that the
!                user is interested in.
!
!  14 June 2000: E. Kemp
!                Added code to output MOS data.
!
!  23 Mar 2002:  Nick Mirsky
!
!                1) Upgraded code to f90
!
!                2) Added sounding data for mandatory pressure levels
!                   in SUBROUTINE SFCDUMP using simple linear
!                   interpolation scheme (as in SUBROUTINE PRODUMP)
!
!                3) After #1 and #2, copied code "arpshis2ver.f90" to
!                   "arpshis2ver_sort.f90".  The new code sorts and
!                   appends data to individual station files
!
!  29 Mar 2002:  Nick Mirsky
!
!  04/04/02: Jason Levit
!  Changed original code into a subroutine, to be used by the
!  arpsverif driver program.
!
!  1 June 2002 Eric Kemp
!  Soil variable updates.
!
!  18 Nov 2002 Kevin W. Thomas
!  As written, the code processes every station in the domain.  For a CONUS
!  run of 36 hours, the run time is close to two hours!  This update allows
!  a user list to indicate which stations to process.
!
!  8 August 2005 Kevin W. Thomas.
!  Update the code to run MPI.
!
!-----------------------------------------------------------------------
!
!  DATA ARRAYS READ IN:
!
!    x        x coordinate of grid points in physical/comp. space (m)
!    y        y coordinate of grid points in physical/comp. space (m)
!    z        z coordinate of grid points in computational space (m)
!    zp       z coordinate of grid points in physical space (m)
!    zpsoil   z coordinate of soil model
!
!    uprt     Perturbation x component of velocity (m/s)
!    vprt     Perturbation y component of velocity (m/s)
!    wprt     Perturbation z component of velocity (m/s)
!
!    ptprt    Perturbation potential temperature (K)
!    pprt     Perturbation pressure (Pascal)
!
!    qc       Cloud water mixing ratio (kg/kg)
!    qr       Rainwater mixing ratio (kg/kg)
!    qi       Cloud ice mixing ratio (kg/kg)
!    qs       Snow mixing ratio (kg/kg)
!    qh       Hail mixing ratio (kg/kg)
!
!    tke      Turbulent Kinetic Energy ((m/s)**2)
!    kmh      Horizontal turb. mixing coef. for momentum ( m**2/s )
!    kmv      Vertical turb. mixing coef. for momentum ( m**2/s )
!
!    ubar     Base state x velocity component (m/s)
!    vbar     Base state y velocity component (m/s)
!    wbar     Base state z velocity component (m/s)
!    ptbar    Base state potential temperature (K)
!    pbar     Base state pressure (Pascal)
!    rhobar   Base state density (kg/m**3)
!    qvbar    Base state water vapor mixing ratio (kg/kg)
!
!    soiltyp  Soil type
!    vegtyp   Vegetation type
!    lai      Leaf Area Index
!    roufns   Surface roughness
!    veg      Vegetation fraction
!
!    tsoil    soil temperature (K)
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
!    radswnet Net shortwave radiation, SWin - SWout
!    radlwin  Incoming longwave radiation
!
!    usflx    Surface flux of u-momentum (kg/(m*s**2))
!    vsflx    Surface flux of v-momentum (kg/(m*s**2))
!    ptsflx   Surface heat flux (K*kg/(m**2 * s ))
!    qvsflx   Surface moisture flux of (kg/(m**2 * s))
!
!  CALCULATED DATA ARRAYS:
!
!    su       Sounding x component of velocity (m/s)
!    sv       Sounding y component of velocity (m/s)
!    stheta   Sounding potential temperature (K)
!    spres    Sounding pressure (mb)
!    stemp    Sounding temperature (C)
!    sdewp    Sounding dew-point (C)
!    sdrct    Sounding wind direction (degrees)
!    ssped    Sounding wind speed (m/s)
!    shght    Sounding height (m)
!
!  WORK ARRAYS:
!
!    tem1     Temporary work array.
!
!   Temporary arrays are defined and used differently by each
!   subroutine.
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
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'vericst.inc'
  INCLUDE 'mp.inc'
!
!-----------------------------------------------------------------------
!
!  Arrays to be read in:
!
!-----------------------------------------------------------------------
!
  INTEGER :: nx, ny, nz, nzsoil
  REAL :: x     (nx)     ! The x-coord. of the physical and
                                     ! computational grid. Defined at u-point.
  REAL :: y     (ny)     ! The y-coord. of the physical and
                                     ! computational grid. Defined at v-point.
  REAL :: z     (nz)     ! The z-coord. of the computational grid.
                                     ! Defined at w-point on the staggered grid.
  REAL :: zp    (nx,ny,nz) ! The physical height coordinate defined at
                                     ! w-point of the staggered grid.
  REAL :: zpsoil(nx,ny,nzsoil) ! The physical height coordinate of soil model
!
  REAL :: hterain (nx,ny)       ! Terrain height
!
  REAL :: uprt   (nx,ny,nz)    ! Perturbation u-velocity (m/s)
  REAL :: vprt   (nx,ny,nz)    ! Perturbation v-velocity (m/s)
  REAL :: wprt   (nx,ny,nz)    ! Perturbation w-velocity (m/s)
  REAL :: ptprt  (nx,ny,nz)    ! Perturbation potential temperature (K)
  REAL :: pprt   (nx,ny,nz)    ! Perturbation pressure (Pascal)
  REAL :: qvprt  (nx,ny,nz)    ! Perturbation water vapor specific
                                         ! humidity (kg/kg)
  REAL :: qscalar(nx,ny,nz,nscalar)

  REAL :: tke    (nx,ny,nz)    ! Turbulent Kinetic Energy ((m/s)**2)
  REAL :: kmh    (nx,ny,nz)    ! Horizontal turb. mixing coef. for
                                         ! momentum. ( m**2/s )
  REAL :: kmv    (nx,ny,nz)    ! Vertical turb. mixing coef. for
                                         ! momentum. ( m**2/s )

  REAL :: ubar   (nx,ny,nz)    ! Base state u-velocity (m/s)
  REAL :: vbar   (nx,ny,nz)    ! Base state v-velocity (m/s)
  REAL :: wbar   (nx,ny,nz)    ! Base state w-velocity (m/s)
  REAL :: ptbar  (nx,ny,nz)    ! Base state potential temperature (K)
  REAL :: pbar   (nx,ny,nz)    ! Base state pressure (Pascal)
  REAL :: rhobar (nx,ny,nz)    ! Base state air density (kg/m**3)
  REAL :: qvbar  (nx,ny,nz)    ! Base state water vapor specific

  INTEGER :: nstyps

  INTEGER :: soiltyp (nx,ny,nstyps)       ! Soil type
  REAL :: stypfrct(nx,ny,nstyps)          ! Soil type
  INTEGER :: vegtyp  (nx,ny)         ! Vegetation type
  REAL :: lai     (nx,ny)            ! Leaf Area Index
  REAL :: roufns  (nx,ny)            ! Surface roughness
  REAL :: veg     (nx,ny)            ! Vegetation fraction

  REAL :: tsoil  (nx,ny,nzsoil,0:nstyps) ! Soil temperature (K)
  REAL :: qsoil  (nx,ny,nzsoil,0:nstyps) ! Soil moisture (m**3/m**3)
  REAL :: wetcanp(nx,ny,0:nstyps)       ! Canopy water amount
  REAL :: snowdpth(nx,ny)        ! Snow depth (m)

  REAL :: raing(nx,ny)         ! Grid supersaturation rain
  REAL :: rainc(nx,ny)         ! Cumulus convective rain
  REAL :: prcrate(nx,ny,4)     ! precipitation rate (kg/(m**2*s))
                               !   prcrate(1,1,1) = total precip. rate
                               !   prcrate(1,1,2) = grid scale precip. rate
                               !   prcrate(1,1,3) = cumulative precip. rate
                               !   prcrate(1,1,4) = microphysics precip. rate

  REAL :: radfrc(nx,ny,nz)     ! Radiation forcing (K/s)
  REAL :: radsw (nx,ny)        ! Solar radiation reaching the surface
  REAL :: rnflx (nx,ny)        ! Net radiation flux absorbed by surface
  REAL :: radswnet(nx,ny)      ! Net shortwave radiation
  REAL :: radlwin(nx,ny)       ! Incoming longwave radiation

  REAL :: usflx (nx,ny)        ! Surface flux of u-momentum (kg/(m*s**2))
  REAL :: vsflx (nx,ny)        ! Surface flux of v-momentum (kg/(m*s**2))
  REAL :: ptsflx(nx,ny)        ! Surface heat flux (K*kg/(m*s**2))
  REAL :: qvsflx(nx,ny)        ! Surface moisture flux (kg/(m**2*s)

  INTEGER :: nhisfile

  REAL :: model_data(sfcmax,nhisfile,5)
  REAL :: obsrv_data(sfcmax,nhisfile,5)

!
!
!-----------------------------------------------------------------------
!
!  Map variables, declared here for use in subroutines...
!
!-----------------------------------------------------------------------
!
  REAL, ALLOCATABLE :: xs(:)
  REAL, ALLOCATABLE :: ys(:)
  REAL, ALLOCATABLE :: xmap(:)
  REAL, ALLOCATABLE :: ymap(:)
  REAL, ALLOCATABLE :: latgr(:,:)
  REAL, ALLOCATABLE :: longr(:,:)
!
!-----------------------------------------------------------------------
!
!  Work Arrays
!
!-----------------------------------------------------------------------
!
  REAL, ALLOCATABLE :: tem1(:,:,:)
  REAL, ALLOCATABLE :: tem2(:,:,:)
  REAL, ALLOCATABLE :: tem3(:,:,:)
  REAL, ALLOCATABLE :: tem4(:,:)
!
!-----------------------------------------------------------------------
!
!  Misc. internal variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: len1, istatus
  REAL :: time
  INTEGER :: i,j
  INTEGER :: ireturn,lengbf,lenfil,nchin,lengbf00Z,lengbf12Z
  INTEGER :: hinfmt
  INTEGER :: nhisfile_max
  PARAMETER (nhisfile_max=200)
  CHARACTER (LEN=256) :: filename
  CHARACTER (LEN=256) :: fgrdbasfn
  CHARACTER (LEN=256) :: grdbasfn
  CHARACTER (LEN=256) :: hisfile(nhisfile_max)
!
  CHARACTER (LEN=256) :: snddmpfn,prodmpfn,sfcdmpfn,precdmpfn,          &
                mosdmpfn
  LOGICAL :: needsndstns,needprostns,needsfcstns,needprecstns,          &
          needmosstns
  LOGICAL :: fexist
  INTEGER :: nf
  CHARACTER (LEN=8) :: the_date
  CHARACTER (LEN=10) :: the_time
  CHARACTER (LEN=5) :: the_zone
  INTEGER :: the_values(8)

  INTEGER :: sndopt
  CHARACTER (LEN=256) :: sndlist
  CHARACTER (LEN=256) :: sndrunname
  CHARACTER (LEN=256) :: snddomlist
  INTEGER :: proopt
  CHARACTER (LEN=256) :: prolist
  CHARACTER (LEN=256) :: prorunname
  CHARACTER (LEN=256) :: prodomlist
  INTEGER :: sfcopt
  CHARACTER (LEN=256) :: sfclist
  CHARACTER (LEN=256) :: sfcobsdir, sfcpre, sfcpost
  CHARACTER (LEN=256) :: mesoobsdir, mesopre, mesopost
  CHARACTER (LEN=256) :: sfcrunname
  CHARACTER (LEN=256) :: blackfile
  INTEGER :: precveropt
  CHARACTER (LEN=256) :: preclist
  CHARACTER (LEN=256) :: precrunname
  CHARACTER (LEN=256) :: precdomlist
  INTEGER :: mosopt
  CHARACTER (LEN=256) :: moslist
  CHARACTER (LEN=256) :: mosrunname
  CHARACTER (LEN=256) :: mosdomlist
  INTEGER :: arpsnn_opt
  CHARACTER (LEN=256) :: nnrunname
  CHARACTER (LEN=256) :: wtsdir
  CHARACTER (LEN=256) :: nnoutputfn

  CHARACTER (LEN=256) :: hdfname

  INTEGER :: sd_id,istat

!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

!  CALL DATE_AND_TIME(the_date,the_time,the_zone,the_values)

!  DO nf=1,nhisfile
!    lenfil =256
!    CALL slength( hisfile(nf), lenfil)
!    WRITE(6,'(1x,a,a)') 'History file is ',hisfile(nf)(1:lenfil)
!  END DO

  nstyp = nstyps

!  ALLOCATE
!
  ALLOCATE(xs (nx),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:xs")
  xs = 0.0
  ALLOCATE(ys (ny),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:ys")
  ys = 0.0
  ALLOCATE(tem1 (nx,ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:tem1")
  tem1 = 0.0
  ALLOCATE(tem2 (nx,ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:tem2")
  tem2 = 0.0
  ALLOCATE(tem3 (nx,ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:tem3")
  tem3 = 0.0
  ALLOCATE(tem4 (nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:tem4")
  tem4 = 0.0
  ALLOCATE(xmap (nx),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:xmap")
  xmap = 0.0
  ALLOCATE(ymap (ny),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:ymap")
  ymap = 0.0
  ALLOCATE(latgr (nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:latgr")
  latgr = 0.0
  ALLOCATE(longr (nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:longr")
  longr = 0.0

!-----------------------------------------------------------------------
!
!  Get selections for verification.
!
!-----------------------------------------------------------------------
!

!  READ(5,sndverif,ERR=100)
!  WRITE(6,*)'sndopt = ',sndopt
!  IF (sndopt == 1) THEN
!    lenfil =256
!    CALL slength( sndlist, lenfil)
!    WRITE(6,'(1x,a,a)') 'Sounding list file was ',sndlist(1:lenfil)
!  END IF
!
!  READ(5,proverif,ERR=100)
!  WRITE(6,*)'proopt = ',proopt
!  IF (proopt == 1) THEN
!    lenfil =256
!    CALL slength( prolist, lenfil)
!    WRITE(6,'(1x,a,a)') 'Sounding list file was ',prolist(1:lenfil)
!  END IF
!
!  READ(5,sfcverif,ERR=100)
!  WRITE(6,*)'sfcopt = ',sfcopt
!  IF (sfcopt == 1) THEN
!    lenfil =256
!    CALL slength( sfclist, lenfil)
!    WRITE(6,'(1x,a,a)') 'Surface ob list file was ',                    &
!                        sfclist(1:lenfil)
!  END IF
!
!  READ(5,sfcobsverif,ERR=100)
!  WRITE(6,*)'sfcveropt = ',sfcveropt
!  IF (sfcveropt == 1) THEN
!    lenfil =256
!    CALL slength( sfcobslist, lenfil)
!    WRITE(6,'(1x,a,a)') 'Surface obs verif list file was ',              &
!                        sfcobslist(1:lenfil)
!  END IF
!
!  READ(5,precverif,ERR=100)
!  WRITE(6,*)'precopt = ',precopt
!  IF (precopt == 1) THEN
!    lenfil =256
!    CALL slength(preclist, lenfil)
!    WRITE(6,'(1x,a,a)') 'Surface ob list file was ',                    &
!                        preclist(1:lenfil)
!  END IF
!
!  READ(5,mosverif,ERR=100)
!  WRITE(6,*)'mosopt = ',mosopt
!  IF (mosopt == 1) THEN
!    lenfil =256
!    CALL slength(moslist, lenfil)
!    WRITE(6,'(1x,a,a)') 'Surface ob list file was ',                    &
!                        moslist(1:lenfil)
!  END IF
!
!
!-----------------------------------------------------------------------
!
!  Loop through all history dumps
!
!-----------------------------------------------------------------------
!
  grdbasfn = fgrdbasfn
  lengbf=256
  CALL slength(grdbasfn,lengbf)
  IF ( mp_opt > 0 .AND. readsplit(FINDX_H) <= 0 ) THEN
    CALL gtsplitfn(fgrdbasfn,1,1,loc_x,loc_y,1,1,                      &
                   0,0,1,lvldbg,grdbasfn,ireturn)
    lengbf = LEN_TRIM(grdbasfn)
  END IF

  DO nf = 1, nhisfile
!
!-----------------------------------------------------------------------
!
!    Read all input data arrays
!
!-----------------------------------------------------------------------
!
    filename=hisfile(nf)
    lenfil =256
    CALL slength( hisfile(nf), lenfil)
    IF ( mp_opt > 0 .AND. readsplit(FINDX_H) <= 0 ) THEN
      CALL gtsplitfn(hisfile(nf),1,1,loc_x,loc_y,1,1,                   &
                     0,0,1,lvldbg,filename,ireturn)
      lenfil = LEN_TRIM(filename)
    END IF
    IF ( myproc == 0 )                                                  &
        WRITE(6,'(1x,a,a)') 'History file is ',filename(1:lenfil)

!
!-----------------------------------------------------------------------
!
!    Subroutine DTAREAD aborts if the file doesn't exist.  For ARPSVERIF,
!    this is very bad, as no output file will ever be written.  Trap it
!    here.
!-----------------------------------------------------------------------
!

    INQUIRE(FILE=trim(filename(1:lenfil)), EXIST = fexist )
    IF( .NOT. fexist ) THEN
      IF (myproc == 0)                                                  &
        WRITE(6,'(a)') " File doesn't exist.  Skipping."
      cycle
    END IF

    CALL dtaread(nx,ny,nz,nzsoil,nstyps,                                &
               hinfmt,nchin,grdbasfn(1:lengbf),lengbf,                  &
               filename(1:lenfil),lenfil,time,                          &
               x,y,z,zp,zpsoil,uprt,vprt,wprt,ptprt,pprt,               &
               qvprt,qscalar,tke,kmh,kmv,                               &
               ubar,vbar,wbar,ptbar,pbar,rhobar,qvbar,                  &
               soiltyp,stypfrct,vegtyp,lai,roufns,veg,                  &
               tsoil,qsoil,wetcanp,snowdpth,                            &
               raing,rainc,prcrate,                                     &
               radfrc,radsw,rnflx,radswnet,radlwin,                     &
               usflx,vsflx,ptsflx,qvsflx,                               &
               ireturn, tem1,tem2,tem3)

    DO j=1,ny
      DO i=1,nx
        hterain(i,j)=zp(i,j,2)
      END DO
    END DO
!
!-----------------------------------------------------------------------
!
!    ireturn = 0 for a successful read
!
!-----------------------------------------------------------------------
!
    IF( ireturn == 0 ) THEN   ! successful read
!
!-----------------------------------------------------------------------
!
!      Extract the sounding data at selected locations.
!
!-----------------------------------------------------------------------
!

      needsndstns=.false.
      IF (sndopt == 1) needsndstns=.true.

      IF (needsndstns) THEN
        len1 =256
        CALL slength( sndrunname, len1)
        snddmpfn=sndrunname(1:len1)//'.snd'//                           &
            filename(lenfil-5:lenfil)

        CALL snddump(nx,ny,nz,                                          &
                 uprt,vprt,                                             &
                 ptprt,pprt,qvprt,                                      &
                 ubar,vbar,ptbar,                                       &
                 pbar,rhobar,qvbar,                                     &
                 x,y,z,zp,hterain,                                      &
                 roufns,tsoil(1,1,1,0),                                 &
                 snddmpfn,needsndstns,                                  &
                 xs,ys,xmap,ymap,                                       &
                 latgr,longr,tem1,tem4,                                 &
                 sndlist,snddomlist,time)
      END IF
!
!-----------------------------------------------------------------------
!
!      Extract the profiler data corresponding to the profiler network.
!      (wind data at the *standard* levels)
!
!-----------------------------------------------------------------------
!
      needprostns=.false.
      IF (proopt == 1) needprostns=.true.

      IF (needprostns) THEN

        len1 =256
        CALL slength( prorunname, len1)
        prodmpfn=prorunname(1:len1)//'.pro'//                           &
            filename(lenfil-5:lenfil)

        CALL prodump(nx,ny,nz,                                          &
                 uprt,vprt,                                             &
                 ptprt,pprt,qvprt,                                      &
                 ubar,vbar,ptbar,                                       &
                 pbar,rhobar,qvbar,                                     &
                 x,y,z,zp,hterain,                                      &
                 roufns,tsoil(1,1,1,0),                                 &
                 prodmpfn,needprostns,                                  &
                 xs,ys,xmap,ymap,                                       &
                 latgr,longr,tem1,tem4,                                 &
                 prolist,prodomlist,time)
      END IF
!
!-----------------------------------------------------------------------
!
!      Extract the surface data corresponding to the observations.
!
!-----------------------------------------------------------------------
!

        CALL sfcdump(nx,ny,nz,                                          &
                 uprt,vprt,                                             &
                 ptprt,pprt,qvprt,                                      &
                 ubar,vbar,ptbar,                                       &
                 pbar,rhobar,qvbar,                                     &
                 x,y,z,zp,hterain,                                      &
                 roufns,tsoil(1,1,1,0),                                 &
                 sfcdmpfn,needsfcstns,nhisfile,nf,time,                 &
                 raing,rainc,                                           &
                 xs,ys,xmap,ymap,                                       &
                 latgr,longr,tem1,                                      &
                 sfclist,                                               &
                 sfcobsdir,sfcpre,sfcpost,                              &
                 mesoobsdir,mesopre,mesopost,                           &
                 blackfile,arpsnn_opt,model_data,obsrv_data)

!
!-----------------------------------------------------------------------
!
!      Extract the precipitation data corresponding to the
!      observations.
!
!-----------------------------------------------------------------------
!
      needprecstns=.false.
      IF (precveropt == 1) needprecstns=.true.

      IF (needprecstns) THEN

        len1 =256
        CALL slength( precrunname, len1)
        precdmpfn=precrunname(1:len1)//'.prec'//                        &
            filename(lenfil-5:lenfil)

        CALL precdump(nx,ny,nz,                                         &
                 uprt,vprt,                                             &
                 ptprt,pprt,qvprt,                                      &
                 ubar,vbar,ptbar,                                       &
                 pbar,rhobar,qvbar,                                     &
                 x,y,z,zp,hterain,                                      &
                 roufns,tsoil(1,1,1,0),                                 &
                 precdmpfn,needprecstns,                                &
                 raing,rainc,                                           &
                 xs,ys,xmap,ymap,                                       &
                 latgr,longr,tem1,                                      &
                 preclist,precdomlist,time)
      END IF

!
!-----------------------------------------------------------------------
!
!      Extract model data for MOS.
!
!-----------------------------------------------------------------------
!
      needmosstns=.false.
      IF (mosopt == 1) needmosstns=.true.

      IF (needmosstns) THEN

        len1 =256
        CALL slength( mosrunname, len1)
        mosdmpfn=mosrunname(1:len1)//'.mos'//                           &
            filename(lenfil-5:lenfil)

        CALL mosdump(nx,ny,nz,nzsoil,nstyps,x,y,z,zp,hterain,           &
                 uprt,ubar,vprt,vbar,wprt,wbar,                         &
                 ptprt,ptbar,pprt,pbar,qvprt,qvbar,                     &
                 qscalar,tke,kmh,kmv,rhobar,                            &
                 zpsoil,tsoil,qsoil,wetcanp,                            &
                 snowdpth,soiltyp,stypfrct,vegtyp,                      &
                 lai,roufns,veg,raing,rainc,prcrate,                    &
                 radfrc,radsw,rnflx,radswnet,radlwin,                   &
                 usflx,vsflx,ptsflx,                                    &
                 qvsflx,mosdmpfn,needmosstns,                           &
                 xs,ys,xmap,ymap,latgr,longr,tem1,tem4,                 &
                 moslist,mosdomlist,time)

      END IF
    ELSE
      WRITE(6,'(a)') ' Error reading data.  HIS2VER ends'
      STOP
    END IF
!
!-----------------------------------------------------------------------
!
!  Go to next history dump
!
!-----------------------------------------------------------------------
!
  END DO
  RETURN
!
!-----------------------------------------------------------------------
!
!  Process error from namelist read attempt.
!
!-----------------------------------------------------------------------
!

  100   CONTINUE
  WRITE(6,*)'ERROR READING NAMELIST!!!'
  STOP

END SUBROUTINE his2ver
!
!#################################################################
!#################################################################
!######                                                     ######
!######                SUBROUTINE SNDDUMP                   ######
!######                                                     ######
!######                Copyright (c) 1996                   ######
!######    Center for Analysis and Prediction of Storms     ######
!######    University of Oklahoma.  All rights reserved.    ######
!######                                                     ######
!#################################################################
!#################################################################
!

SUBROUTINE snddump(nx,ny,nz,                                            &
           uprt,vprt,ptprt,pprt,qvprt,                                  &
           ubar,vbar,ptbar,                                             &
           pbar,rhobar,qvbar,                                           &
           x,y,z,zp,hterain,                                            &
           roufns,tsfc,                                                 &
           dumpfile,needsndstns,                                        &
           xs,ys,xmap,ymap,                                             &
           latgr,longr,tem1,tem4,                                       &
           sndlist,snddomlist,time)
!
!----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Coordinate the verification sounding dumps, which are converted
!  to GEMPAK format using a post-processor.
!
!----------------------------------------------------------------------
!
!  AUTHOR:  John Mewes
!
!  MODIFICATION HISTORY:
!
!  09/15/97  Major code revamping.  Changed variable names and put
!            station information arrays into a common block in order
!            to avoid having to reread information with each new
!            data dump (adds vericst.inc to the verification
!            software).
!
!  19 May 2000:  Eric Kemp
!                File name for list of stations in the domain
!                is now passed to the subroutine, along with
!                time from initialization.  Also, valid
!                time and forecast domain info are now included
!                in all output files.
!
!  23 Mar 2002:  Nick Mirsky
!                Added sounding data for mandatory pressure levels
!                using simple linear interpolation scheme
!                as in SUBROUTINE PRODUMP
!
!-----------------------------------------------------------------------
!
!  Variables to be read in:
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
  INTEGER :: nx,ny,nz          ! Number of grid points in 3 directions
!
  INTEGER :: nzsoil            ! Number of soil layers
!
  LOGICAL :: needsndstns       ! Create a stn list of stns in grid
!
  REAL :: uprt  (nx,ny,nz)     ! Perturbation u velocity (m/s)
  REAL :: vprt  (nx,ny,nz)     ! Perturbation v velocity (m/s)
  REAL :: ptprt (nx,ny,nz)     ! Perturbation potential temperature (K)
  REAL :: pprt  (nx,ny,nz)     ! Perturbation pressure (Pascal)
  REAL :: qvprt (nx,ny,nz)     ! Perturbation wv mixing ration (kg/kg)
!
  REAL :: ubar  (nx,ny,nz)     ! Base state u-velocity (m/s)
  REAL :: vbar  (nx,ny,nz)     ! Base state v-velocity (m/s)
  REAL :: ptbar (nx,ny,nz)     ! Base state potential temperature (K)
  REAL :: pbar  (nx,ny,nz)     ! Base state pressure (Pascal)
  REAL :: rhobar(nx,ny,nz)     ! Base state air density (kg/m**3)
  REAL :: qvbar (nx,ny,nz)     ! Base state water vapor specific humidity
                               ! (kg/kg)
!
  REAL :: x     (nx)           ! The x-coord. of the physical and
                               ! computational grid. Defined at u-point.
  REAL :: y     (ny)           ! The y-coord. of the physical and
                               ! computational grid. Defined at v-point.
  REAL :: z     (nz)           ! The z-coord. of the computational grid.
                               ! Defined at w-point on the staggered grid.
  REAL :: zp    (nx,ny,nz)     ! The physical height coordinate defined at
                               ! w-point of the staggered grid.
!
  REAL :: hterain(nx,ny)       ! Terrain height.
!
  REAL :: roufns (nx,ny)       ! Surface roughness
!
  REAL :: tsfc(nx,ny)          ! Ground sfc. temperature (K)

  REAL :: time
!
!
!-----------------------------------------------------------------------
!
!  Computed variables
!
!-----------------------------------------------------------------------
!
  REAL :: xs(nx)      ! x location of scalar points
  REAL :: ys(ny)      ! y location of scalar points
!
!-----------------------------------------------------------------------
!
!  Extracted sounding variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: nlevs
  INTEGER :: nzmax
  PARAMETER (nzmax=250)
!
  REAL :: su(nzmax)
  REAL :: sv(nzmax)
  REAL :: stheta(nzmax)
  REAL :: sqv(nzmax)
  REAL :: spres(nzmax)
  REAL :: stemp(nzmax)
  REAL :: sdewp(nzmax)
  REAL :: sdrct(nzmax)
  REAL :: ssped(nzmax)
  REAL :: shght(nzmax)
  REAL :: srain
!
!-----------------------------------------------------------------------
!
!  Work Arrays
!
!-----------------------------------------------------------------------
!
  REAL :: tem1(nx,ny,nz)
  REAL :: tem4(nx,ny)
!
!-----------------------------------------------------------------------
!
!  Misc. internal variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: istnm      ! Station number of the sounding location
  INTEGER :: iselev_act ! Actual station elevation in integer format
  INTEGER :: i,j,k      ! Index variables
  INTEGER :: ii,ij      ! More index variables
  INTEGER :: ireturn    ! Flag for stations outside of the domain
!
  CHARACTER (LEN=256) :: dumpfile     ! Sounding data history dump file name
  CHARACTER (LEN=256) :: sndlist   ! File to read the sounding locations from
  CHARACTER (LEN=256) :: line         ! Temp variable to read lines from files
  CHARACTER (LEN=256) :: snddomlist
  INTEGER :: LEN
!
!-----------------------------------------------------------------------
!
!  Map projection variables
!
!-----------------------------------------------------------------------
!
  REAL :: xctr,yctr
  REAL :: xll,yll
  REAL :: latnot(2)
  REAL :: xmap(nx)
  REAL :: ymap(ny)
  REAL :: latgr(nx,ny)
  REAL :: longr(nx,ny)
!
!-----------------------------------------------------------------------
!
!  Variables used in interpolation for mandatory levels
!
!-----------------------------------------------------------------------
!
  REAL :: topwt,botwt,mandu,mandv
  REAL :: mandlev
  REAL :: mandhght,mandtemp,manddewp,manddir,mandspd
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'phycst.inc'
  INCLUDE 'globcst.inc'
  INCLUDE 'bndry.inc'
  INCLUDE 'exbc.inc'
  INCLUDE 'vericst.inc'
  INCLUDE 'grid.inc'
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  IF (nz > nzmax) THEN
    PRINT *,'Reset nzmax to greater than:',nz
    STOP
  END IF
!
!-----------------------------------------------------------------------
!
!  Calculate scalar locations and set up the map projection and
!  grid parameters
!
!-----------------------------------------------------------------------
!
  DO i=1,nx-1
    xs(i)=0.5*(x(i)+x(i+1))
  END DO
  DO j=1,ny-1
    ys(j)=0.5*(y(j)+y(j+1))
  END DO
  dx=x(2)-x(1)
  dy=y(2)-y(1)
!
  latnot(1)=trulat1
  latnot(2)=trulat2
  CALL setmapr(mapproj,sclfct,latnot,trulon)
  CALL lltoxy(1,1,ctrlat,ctrlon,xctr,yctr)
  xll=xctr-(0.5*(nx-3)*dx)
  yll=yctr-(0.5*(ny-3)*dy)

  DO i=1,nx-1
    xmap(i)=xll+xs(i)
  END DO
  xmap(nx)=2.*xmap(nx-1)-xmap(nx-2)
  DO j=1,ny-1
    ymap(j)=yll+ys(j)
  END DO
  ymap(ny)=2.*ymap(ny-1)-ymap(ny-2)
  CALL xytoll(nx,ny,xmap,ymap,latgr,longr)
!
!-----------------------------------------------------------------------
!
!  Find location in ARPS grid of all the stations in sndlist, then
!  rewrite only the ones that are in the grid to common arrays
!  so as to not make the program read all of the stations from
!  the original sounding stn list and check their locations with each
!  history dump.
!
!-----------------------------------------------------------------------
!

  OPEN(UNIT=1,FILE=sndlist,STATUS='old')
!  OPEN(UNIT=2,FILE=dumpfile,STATUS='unknown')
!
!  LEN =256
!  CALL slength( dumpfile, LEN)
!  WRITE(6,*) "Opening ",dumpfile(1:LEN)," for sounding output"

  105  FORMAT(a8,a16,a20)
!
  IF (needsndstns) THEN
!
    needsndstns = .false.
    sndstn=0
!
    LEN =80
    CALL slength(sndlist,LEN)
    WRITE(6,*) 'Reading sounding stations from file: ',                 &
                sndlist(1:LEN)
!
    OPEN(UNIT=3,FILE=snddomlist,STATUS='unknown')

    WRITE(3,990) year,'-',month,'-',day,'.',hour,':',minute,':',        &
               second,'+',INT(time)

    LEN =80
    CALL slength(runname,LEN)
    WRITE(3,'(a)') runname(1:LEN)
    WRITE(3,'(i1.1)') nocmnt
    DO i = 1,nocmnt
      LEN =80
      CALL slength(cmnt(i),LEN)
      WRITE(3,'(a)') cmnt(i)(1:LEN)
    END DO

!
    110    FORMAT(a80)
!
    READ(1,110,END=30) line ! read header lines and discard them..
    READ(1,110,END=30) line
    READ(1,110,END=30) line
!
    DO i=1,sndmax
!
      115      FORMAT(a3,2X,f5.2,2X,f7.2,1X,i4)
!
      READ(1,110,END=30) line
!
      sndstn=sndstn+1
!
      READ(line,115) sndstid(sndstn),sndlat(sndstn),sndlon(sndstn),     &
                     iselev_act
      istnm=0 !  Assumes a station number is unavailable
      sndelev_act(sndstn)=FLOAT(iselev_act)
!
      CALL lltoxy(1,1,sndlat(sndstn),sndlon(sndstn),                    &
                       sndxpt(sndstn),sndypt(sndstn))
!
      sndxpt(sndstn)=sndxpt(sndstn)-xll
      sndypt(sndstn)=sndypt(sndstn)-yll
!
      CALL findlc2(nx,ny,xs,ys,sndxpt(sndstn),sndypt(sndstn),           &
                 sndipt(sndstn),sndjpt(sndstn),ireturn)
!
      IF (ireturn < 0) THEN  ! stn is outside the grid
        sndstn=sndstn-1
      ELSE  ! stn is inside the grid
        WRITE(3,110) line ! write location data to a file
      END IF
!
    END DO
!
  END IF
!
!-----------------------------------------------------------------------
!
!  Interpolate (in the horizontal) for the whole vertical column
!  for each station, derive the sounding variables, and write the
!  extracted sounding to dumpfile...
!
!-----------------------------------------------------------------------
!
  30   CONTINUE

!  WRITE(2,990) year,'-',month,'-',day,'.',hour,':',minute,':',          &
!               second,'+',INT(time)
!  990   FORMAT(i4.4,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i6.6)

!  LEN =80
!  CALL slength(runname,LEN)
!  WRITE(2,'(a)') runname(1:LEN)
!  WRITE(2,'(i1.1)') nocmnt
!  DO i = 1,nocmnt
!    LEN =80
!    CALL slength(cmnt(i),LEN)
!    WRITE(2,'(a)') cmnt(i)(1:LEN)
!  END DO

!
  DO i=1,sndstn ! do for each station 'i' ...
!
    dumpfile = './snd_dumps/'
    LEN =256
    CALL slength( dumpfile, LEN)
    WRITE(dumpfile(LEN+1:LEN+7),'(a7)') sndstid(i)(1:3)//'_snd'
    OPEN(UNIT=2,FILE=dumpfile,STATUS='unknown',POSITION='append')

    WRITE(2,990) year,'-',month,'-',day,'.',hour,':',minute,':',        &
               second,'+',INT(time)
    990   FORMAT(i4.4,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i6.6)

    DO ii=1,nx
      DO ij=1,ny
        tem4(ii,ij)=0.
      END DO
    END DO
!
!-----------------------------------------------------------------------
!
!  Interpolate (in the horizontal) for the whole vertical column.
!
!-----------------------------------------------------------------------
!
    CALL colintb(nx,ny,nz,nzmax,                                        &
                xs,ys,zp,sndxpt(i),sndypt(i),sndipt(i),sndjpt(i),       &
                uprt, vprt, ptprt, pprt, qvprt,                         &
                ubar, vbar, ptbar, pbar, qvbar,                         &
                tem4,tem4,srain,                                        &
                su,sv,stheta,spres,shght,sqv,                           &
                tem1,nlevs)


!
!-----------------------------------------------------------------------
!
!  Convert sounding to desired units/quantities (winds: m/s  theta: K
!  temp/dewp: C  press: mb  qv: kg/kg)
!
!-----------------------------------------------------------------------
!
    CALL cnvsnd(su,sv,stheta,spres,sqv,sndlon(i),                       &
                sdrct,ssped,stemp,sdewp,nlevs)
!
!-----------------------------------------------------------------------
!
!  Write out soundings for the current model time..
!
!-----------------------------------------------------------------------
!
    WRITE(2,120) sndstid(i),(nlevs-2)
    120    FORMAT(a4,3X,i2)
!
    WRITE(2,125)
    125    FORMAT(6X,                                                   &
           'PRES     HGHT     TMPC     DWPC     DRCT     SPED')
!
!-----------------------------------------------------------------------
!
!  Perform a simple linear interpolation to retrieve data at the
!  mandatory levels.  This replaces the 1-d Barnes' analysis.  The
!  replacement was made for computational speed considerations (gets
!  rid of the exponentials) and because the possible errors intro-
!  duced through linear interpolation are likely small compared to
!  the error of observations in sounding data.
!
!-----------------------------------------------------------------------
!
    DO k=2,(nlevs-1)
!
       mandlev=0.0
       IF (spres(k) > 1000.0 .AND. spres(k+1) < 1000.0) THEN
          mandlev=1000.0
       ELSE IF (spres(k) > 925.0 .AND. spres(k+1) < 925.0) THEN
          mandlev=925
       ELSE IF (spres(k) > 850.0 .AND. spres(k+1) < 850.0) THEN
          mandlev=850.0
       ELSE IF (spres(k) > 700.0 .AND. spres(k+1) < 700.0) THEN
          mandlev=700.0
       ELSE IF (spres(k) > 500.0 .AND. spres(k+1) < 500.0) THEN
          mandlev=500.0
       ELSE IF (spres(k) > 400.0 .AND. spres(k+1) < 400.0) THEN
          mandlev=400.0
       ELSE IF (spres(k) > 300.0 .AND. spres(k+1) < 300.0) THEN
          mandlev=300.0
       ELSE IF (spres(k) > 250.0 .AND. spres(k+1) < 250.0) THEN
          mandlev=250.0
       ELSE IF (spres(k) > 200.0 .AND. spres(k+1) < 200.0) THEN
          mandlev=200.0
       ELSE IF (spres(k) > 150.0 .AND. spres(k+1) < 150.0) THEN
          mandlev=150.0
       ELSE IF (spres(k) > 100.0 .AND. spres(k+1) < 100.0) THEN
          mandlev=100.0
       ELSE IF (spres(k) > 70.0 .AND. spres(k+1) < 70.0) THEN
          mandlev=70.0
       ELSE IF (spres(k) > 50.0 .AND. spres(k+1) < 50.0) THEN
          mandlev=50.0
       ELSE IF (spres(k) > 10.0 .AND. spres(k+1) < 10.0) THEN
          mandlev=10.0
       END IF
!
       IF (mandlev/=0) THEN
!
          topwt = (mandlev-spres(k))/(spres(k+1)-spres(k))
          botwt = 1.0 - topwt
!
          mandhght = botwt*shght(k) + topwt*shght(k+1)
          mandtemp = botwt*stemp(k) + topwt*stemp(k+1)
          manddewp = botwt*sdewp(k) + topwt*sdewp(k+1)
          mandu = botwt*su(k) + topwt*su(k+1)
          mandv = botwt*sv(k) + topwt*sv(k+1)
!
          CALL get_ddff(mandu,mandv,manddir,mandspd)
!
          WRITE(2,130) spres(k),shght(k),stemp(k),sdewp(k),             &
                 sdrct(k),ssped(k)
          130      FORMAT(1X,6(f9.2))

          WRITE(2,131) mandlev,mandhght,mandtemp,manddewp,              &
                 manddir,mandspd
          131      FORMAT(1X,6(f9.2))
!
       ELSE
!
          WRITE(2,132) spres(k),shght(k),stemp(k),sdewp(k),             &
                 sdrct(k),ssped(k)
          132      FORMAT(1X,6(f9.2))
!
       END IF
!
    END DO
!
    CLOSE(2)
!
  END DO ! Move on to the next station..
!
  CLOSE(1)
!  CLOSE(2)
  CLOSE(3)
!
  RETURN
!
END SUBROUTINE snddump
!
!
!#################################################################
!#################################################################
!######                                                     ######
!######               SUBROUTINE PRODUMP                    ######
!######                                                     ######
!######                Copyright (c) 1996                   ######
!######    Center for Analysis and Prediction of Storms     ######
!######    University of Oklahoma.  All rights reserved.    ######
!######                                                     ######
!#################################################################
!#################################################################
!

SUBROUTINE prodump(nx,ny,nz,                                            &
           uprt,vprt,ptprt,pprt,qvprt,                                  &
           ubar,vbar,ptbar,                                             &
           pbar,rhobar,qvbar,                                           &
           x,y,z,zp,hterain,                                            &
           roufns,tsfc,                                                 &
           dumpfile,needprostns,                                        &
           xs,ys,xmap,ymap,                                             &
           latgr,longr,tem1,tem4,                                       &
           proflist,profdomlist,time)
!
!----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Coordinate the verification profiler dumps.
!
!----------------------------------------------------------------------
!
!  AUTHOR:  John Mewes
!
!  MODIFICATION HISTORY:
!
!  09/15/97  Major code revamping.  Changed variable names and put
!            station information arrays into a common block in order
!            to avoid having to reread information with each new
!            data dump (adds vericst.inc to the verification
!            software).
!
!  09/15/97  Removed the Barnes' analysis of winds to the standard
!            levels and replaced it with a much simpler (but also
!            much quicker) linear interpolation.
!
!  09/15/97  Changed subroutine name from PROFDUMP to PRODUMP in
!            order to bring the subroutine name into a consistent
!            format (to go with SNDDUMP and SFCDUMP)
!
!  19 May 2000:  Eric Kemp
!                File name for list of stations in the domain
!                is now passed to the subroutine, along with
!                time from initialization.  Also, valid
!                time and forecast domain info are now included
!                in all output files.
!
!-----------------------------------------------------------------------
!
!  Variables to be read in:
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
  INCLUDE 'bndry.inc'
  INCLUDE 'exbc.inc'
  INCLUDE 'vericst.inc'
  INCLUDE 'grid.inc'
!

!
  INTEGER :: nx,ny,nz          ! Number of grid points in 3 directions
!
  INTEGER :: nzsoil
!
  LOGICAL :: needprostns       ! Need to look up profiler information?
!
  REAL :: uprt  (nx,ny,nz)     ! Perturbation u velocity (m/s)
  REAL :: vprt  (nx,ny,nz)     ! Perturbation v velocity (m/s)
  REAL :: ptprt (nx,ny,nz)     ! Perturbation potential temperature (K)
  REAL :: pprt  (nx,ny,nz)     ! Perturbation pressure (Pascal)
  REAL :: qvprt (nx,ny,nz)     ! Perturbation wv mixing ration (kg/kg)
!
  REAL :: ubar  (nx,ny,nz)     ! Base state u-velocity (m/s)
  REAL :: vbar  (nx,ny,nz)     ! Base state v-velocity (m/s)
  REAL :: ptbar (nx,ny,nz)     ! Base state potential temperature (K)
  REAL :: pbar  (nx,ny,nz)     ! Base state pressure (Pascal)
  REAL :: rhobar(nx,ny,nz)     ! Base state air density (kg/m**3)
  REAL :: qvbar (nx,ny,nz)     ! Base state water vapor specific humidity
                               ! (kg/kg)
!
  REAL :: x     (nx)           ! The x-coord. of the physical and
                               ! computational grid. Defined at u-point.
  REAL :: y     (ny)           ! The y-coord. of the physical and
                               ! computational grid. Defined at v-point.
  REAL :: z     (nz)           ! The z-coord. of the computational grid.
                               ! Defined at w-point on the staggered grid.
  REAL :: zp    (nx,ny,nz)     ! The physical height coordinate defined at
                               ! w-point of the staggered grid.
!
  REAL :: hterain(nx,ny)       ! Terrain height.
!
  REAL :: roufns (nx,ny)       ! Surface roughness
!
  REAL :: tsfc(nx,ny)          ! Ground sfc. temperature (K)

  REAL :: time
!
!
!-----------------------------------------------------------------------
!
!  Computed variables
!
!-----------------------------------------------------------------------
!
  REAL :: xs(nx)      ! x location of scalar points
  REAL :: ys(ny)      ! y location of scalar points
!
!-----------------------------------------------------------------------
!
!  Extracted sounding variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: nlevs
  INTEGER :: nzmax
  PARAMETER (nzmax=250)
!
  REAL :: su(nzmax)
  REAL :: sv(nzmax)
  REAL :: stheta(nzmax)
  REAL :: sqv(nzmax)
  REAL :: spres(nzmax)
  REAL :: stemp(nzmax)
  REAL :: sdewp(nzmax)
  REAL :: sdrct(nzmax)
  REAL :: ssped(nzmax)
  REAL :: shght(nzmax)
  REAL :: srain
!
!-----------------------------------------------------------------------
!
!  Inserted for prodump routine
!
!-----------------------------------------------------------------------
!
  INTEGER :: stdlevs       ! Number of standard levels
  PARAMETER (stdlevs=30)
!
  REAL :: deltaz           ! Meters between standard levels
  PARAMETER (deltaz=500.)
!
  REAL :: stdlev(stdlevs)  ! Chosen standard levels (meters AGL)
  REAL :: stddir(stdlevs)  ! Wind direction at std levels
  REAL :: stdspd(stdlevs)  ! Wind speed at std levels (m/s)
  REAL :: stdu(stdlevs)    ! U-component at std levels (m/s)
  REAL :: stdv(stdlevs)    ! V-component at std levels (m/s)
!
!-----------------------------------------------------------------------
!
!  Variables used in estimating winds at the standard levels
!  model's winds
!
!-----------------------------------------------------------------------
!
  REAL :: topwt,botwt
!
!-----------------------------------------------------------------------
!
!  Work Arrays
!
!-----------------------------------------------------------------------
!
  REAL :: tem1(nx,ny,nz)
  REAL :: tem4(nx,ny)
!
!-----------------------------------------------------------------------
!
!  Misc. internal variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: istnm      ! Station number of the sounding location
  INTEGER :: iselev_act ! Actual station elevation in integer format
  INTEGER :: i,j        ! Index variables
  INTEGER :: ii,ij      ! More index variables
  INTEGER :: ireturn    ! Return status for whether stn is in the grid
  INTEGER :: LEN
!
  CHARACTER (LEN=256) :: dumpfile ! Profiler data history dump file name
  CHARACTER (LEN=256) :: proflist ! File to read the profiler locations from
  CHARACTER (LEN=256) :: line     ! Temporary variable to read lines from files
  CHARACTER (LEN=256) :: profdomlist
!
  REAL :: selev_mod     ! Model estimate of station elevation
!
!-----------------------------------------------------------------------
!
!  Map projection variables
!
!-----------------------------------------------------------------------
!
  REAL :: xctr,yctr
  REAL :: xll,yll
  REAL :: latnot(2)
  REAL :: xmap(nx)
  REAL :: ymap(ny)
  REAL :: latgr(nx,ny)
  REAL :: longr(nx,ny)
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  IF (nz > nzmax) THEN
    PRINT *,'Reset nzmax to greater than:',nz
    STOP
  END IF
!
!-----------------------------------------------------------------------
!
!  Calculate scalar locations and set up the map projection and
!  grid parameters
!
!-----------------------------------------------------------------------
!
  DO i=1,nx-1
    xs(i)=0.5*(x(i)+x(i+1))
  END DO
  DO j=1,ny-1
    ys(j)=0.5*(y(j)+y(j+1))
  END DO
  dx=x(2)-x(1)
  dy=y(2)-y(1)
!
  latnot(1)=trulat1
  latnot(2)=trulat2
  CALL setmapr(mapproj,sclfct,latnot,trulon)
  CALL lltoxy(1,1,ctrlat,ctrlon,xctr,yctr)
  xll=xctr-(0.5*(nx-3)*dx)
  yll=yctr-(0.5*(ny-3)*dy)

  DO i=1,nx-1
    xmap(i)=xll+xs(i)
  END DO
  xmap(nx)=2.*xmap(nx-1)-xmap(nx-2)
  DO j=1,ny-1
    ymap(j)=yll+ys(j)
  END DO
  ymap(ny)=2.*ymap(ny-1)-ymap(ny-2)
  CALL xytoll(nx,ny,xmap,ymap,latgr,longr)
!
!-----------------------------------------------------------------------
!
!  Find location in ARPS grid of all the stations in proflist, then
!  rewrite only the ones that are in the grid to common arrays
!  so as to not make the program read all of the stations from
!  the original profiler stn list and check their locations with each
!  history dump.
!
!-----------------------------------------------------------------------
!
  OPEN(UNIT=1,FILE=proflist,STATUS='old')
!  OPEN(UNIT=2,FILE=dumpfile,STATUS='unknown')
!
!  LEN =256
!  CALL slength( dumpfile, LEN)
!  WRITE(6,*) "Opening ",dumpfile(1:LEN)," for profiler output"

  105  FORMAT(a8,a16,a20)
!
  IF (needprostns) THEN
!
    needprostns=.false.
    prostn=0
!
    LEN =80
    CALL slength(proflist,LEN)
    WRITE(6,*) 'Reading profiler stations from file: ',                 &
                proflist(1:LEN)
!
!     open(unit=3,file='prostns.out',status='unknown')
    OPEN(UNIT=3,FILE=profdomlist,STATUS='unknown')


    WRITE(3,990) year,'-',month,'-',day,'.',hour,':',minute,':',        &
               second,'+',INT(time)

    LEN =80
    CALL slength(runname,LEN)
    WRITE(3,'(a)') runname(1:LEN)
    WRITE(3,'(i1.1)') nocmnt
    DO i = 1,nocmnt
      LEN =80
      CALL slength(cmnt(i),LEN)
      WRITE(3,'(a)') cmnt(i)(1:LEN)
    END DO

!
    110    FORMAT(a80)
!
    READ(1,110,END=30) line ! read header lines and discard them..
    READ(1,110,END=30) line
    READ(1,110,END=30) line
!
    DO i=1,promax
!
      115      FORMAT(a3,1X,f6.2,1X,f7.2,1X,i4)
!
      READ(1,110,END=30) line
!
      prostn=prostn+1
!
      READ(line,115) prostid(prostn),prolat(prostn),                    &
                     prolon(prostn),iselev_act
      istnm=0
      proelev_act(prostn)=FLOAT(iselev_act)
!
      CALL lltoxy(1,1,prolat(prostn),prolon(prostn),                    &
                       proxpt(prostn),proypt(prostn))
!
      proxpt(prostn)=proxpt(prostn)-xll
      proypt(prostn)=proypt(prostn)-yll
!
      CALL findlc2(nx,ny,xs,ys,proxpt(prostn),proypt(prostn),           &
                  proipt(prostn),projpt(prostn),ireturn)
!
      IF (ireturn < 0) THEN  ! stn is outside the grid
        prostn=prostn-1
      ELSE  ! stn is within the grid
        WRITE(3,110) line
      END IF
!
    END DO
!
  END IF
!
  30   CONTINUE
!
!-----------------------------------------------------------------------
!
!  Use the sounding extraction subroutine and linear interpolation
!  to extract winds at standard levels (standard levels used in
!  comparing with observed data)
!
!-----------------------------------------------------------------------
!
!  WRITE(2,990) year,'-',month,'-',day,'.',hour,':',minute,':',          &
!               second,'+',INT(time)
!  990   FORMAT(i4.4,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i6.6)
!
!  LEN =80
!  CALL slength(runname,LEN)
!  WRITE(2,'(a)') runname(1:LEN)
!  WRITE(2,'(i1.1)') nocmnt
!  DO i = 1,nocmnt
!    LEN =80
!    CALL slength(cmnt(i),LEN)
!    WRITE(2,'(a)') cmnt(i)(1:LEN)
!  END DO

  DO i=1,prostn  ! do for each station 'i'
!
    dumpfile = './pro_dumps/'
    LEN =256
    CALL slength( dumpfile, LEN)
    WRITE(dumpfile(LEN+1:LEN+7),'(a4)') sndstid(i)//'_pro'
    OPEN(UNIT=2,FILE=dumpfile,STATUS='unknown',POSITION='append')

    WRITE(2,990) year,'-',month,'-',day,'.',hour,':',minute,':',        &
               second,'+',INT(time)
    990   FORMAT(i4.4,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i6.6)

    DO ii=1,nx
      DO ij=1,ny
        tem4(ii,ij)=0.
      END DO
    END DO
!
!-----------------------------------------------------------------------
!
!  Interpolate (in the horizontal) for the whole vertical column.
!
!-----------------------------------------------------------------------
!
    CALL colintb(nx,ny,nz,nzmax,                                        &
                xs,ys,zp,proxpt(i),proypt(i),proipt(i),projpt(i),       &
                uprt, vprt, ptprt, pprt, qvprt,                         &
                ubar, vbar, ptbar, pbar, qvbar,                         &
                tem4,tem4,srain,                                        &
                su,sv,stheta,spres,shght,sqv,                           &
                tem1,nlevs)
!
!-----------------------------------------------------------------------
!
!  Convert sounding to desired units/quantities (winds: m/s  theta: K
!  temp/dewp: C  press: mb  qv: kg/kg)
!
!-----------------------------------------------------------------------
!
    CALL cnvsnd(su,sv,stheta,spres,sqv,prolon(i),                       &
                sdrct,ssped,stemp,sdewp,nlevs)
!
!-----------------------------------------------------------------------
!
!  Write header output for each station
!
!-----------------------------------------------------------------------
!
    selev_mod=(shght(1)+shght(2))/2.
!
    WRITE(2,120) prostid(i),selev_mod,prolat(i),prolon(i),'Model'
    120    FORMAT(2X,a3,3X,f5.0,3X,f5.2,3X,f7.2,3X,a8)
!
    DO ii=1,stdlevs ! set up heights of the std profiler data levels
      stdlev(ii) = ii*deltaz + selev_mod
    END DO
!
!-----------------------------------------------------------------------
!
!  Perform a simple linear interpolation to retrieve winds at the
!  'standard' levels.  This replaces the 1-d Barnes' analysis.  The
!  replacement was made for computational speed considerations (gets
!  rid of the exponentials) and because the possible errors intro-
!  duced through linear interpolation are likely small compared to
!  the error of observation in profiler data.
!
!-----------------------------------------------------------------------
!
    DO ii=1,stdlevs,1
!
      DO ij=1,nlevs
!
        IF (stdlev(ii) > shght(ij-1) .AND. stdlev(ii) <= shght(ij)) THEN
!
          topwt = (stdlev(ii)-shght(ij-1))/(shght(ij)-shght(ij-1))
          botwt = 1.0 - topwt
!
          stdu(ii) = botwt*su(ij-1) + topwt*su(ij)
          stdv(ii) = botwt*sv(ij-1) + topwt*sv(ij)
!
          CALL get_ddff(stdu(ii),stdv(ii),stddir(ii),stdspd(ii))
!
        END IF
!
      END DO
!
    END DO
!
    DO ii=1,stdlevs
      WRITE(2,125) (stdlev(ii)-selev_mod),stddir(ii),                   &
          stdspd(ii)
      125      FORMAT(5X,f6.0,5X,f4.0,5X,f5.1)
    END DO
!
    CLOSE(2)
!
  END DO ! on to the next profiler station ...
!
  CLOSE(1)
!  CLOSE(2)
  CLOSE(3)
!
  RETURN
!
END SUBROUTINE prodump
!
!
!#################################################################
!#################################################################
!######                                                     ######
!######               SUBROUTINE GET_DDFF                   ######
!######                                                     ######
!######             Copyright (c) 1995-1996                 ######
!######    Center for Analysis and Prediction of Storms     ######
!######    University of Oklahoma.  All rights reserved.    ######
!######                                                     ######
!#################################################################
!#################################################################
!

SUBROUTINE get_ddff(u,v,dd,ff)
!
!----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Convert u and v winds into direction and speed.
!
!----------------------------------------------------------------------
!
  IMPLICIT NONE
!
  REAL :: dd,ff,u,v,rad2d,spval,mis_val
!
  PARAMETER (rad2d=57.29577951, spval=9999.,                            &
             mis_val=99999.0)
!
  IF(u < spval .AND. v < spval) THEN
    ff = SQRT((u*u + v*v))
    IF(ff /= 0.) THEN
      dd = rad2d*ATAN2(u,v)
      dd = dd+180.
      IF (dd > 360.) dd=dd-360.
    ELSE
      dd=0.
    END IF
  ELSE
    dd = mis_val
    ff = mis_val
  END IF
!
  RETURN
!
END SUBROUTINE get_ddff
!
!
!#################################################################
!#################################################################
!######                                                     ######
!######                SUBROUTINE SFCDUMP                   ######
!######                                                     ######
!######             Copyright (c) 1995-1996                 ######
!######    Center for Analysis and Prediction of Storms     ######
!######    University of Oklahoma.  All rights reserved.    ######
!######                                                     ######
!#################################################################
!#################################################################
!

SUBROUTINE sfcdump(nx,ny,nz,                                            &
           uprt,vprt,ptprt,pprt,qvprt,                                  &
           ubar,vbar,ptbar,                                             &
           pbar,rhobar,qvbar,                                           &
           x,y,z,zp,hterain,                                            &
           roufns,tsfc,                                                 &
           dumpfile,needsfcstns,nhisfile,nf,time,                       &
           raing,rainc,                                                 &
           xs,ys,xmap,ymap,                                             &
           latgr,longr,tem1,                                            &
           sfclist,                                                     &
           sfcobsdir,sfcpre,sfcpost,                                    &
           mesoobsdir,mesopre,mesopost,                                 &
           blackfile,arpsnn_opt,model_data,obsrv_data)
!
!----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Coordinate the verification surface dumps, which are converted
!  from files arranged by time to files arranged by station by a
!  post processor.
!
!----------------------------------------------------------------------
!
!  AUTHOR:  John Mewes
!
!  MODIFICATION HISTORY:
!
!  09/15/97  Changed subroutine name from MESODUMP to SFCDUMP in
!            order to bring the subroutine name into a consistent
!            format (to go with SNDDUMP and PRODUMP).
!
!  09/16/97  Major code revamping.  Changed variable names and put
!            station information arrays into a common block in order
!            to avoid having to reread information with each new
!            data dump (adds vericst.inc to the verification
!            software).
!
!  19 May 2000:  Eric Kemp
!                File name for list of stations in the domain
!                is now passed to the subroutine, along with
!                time from initialization.  Also, valid
!                time and forecast domain info are now included
!                in all output files.
!
!  08 Apr 2002: Jason Levit
!               Modifications to automatically find SAO files,
!               plus some changes to the input and output
!               for addition to the arpsverif package.
!
!-----------------------------------------------------------------------
!
!  Variables to be read in:
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
  INCLUDE 'bndry.inc'
  INCLUDE 'exbc.inc'
  INCLUDE 'vericst.inc'
  INCLUDE 'grid.inc'
  INCLUDE 'adas.inc'
  INCLUDE 'mp.inc'

!
  INTEGER :: nx,ny,nz          ! Number of grid points in 3 directions
!
  INTEGER :: nzsoil
!
  LOGICAL :: needsfcstns       ! Whether not need to create station file
!
  REAL :: uprt  (nx,ny,nz)     ! Perturbation u velocity (m/s)
  REAL :: vprt  (nx,ny,nz)     ! Perturbation v velocity (m/s)
  REAL :: ptprt (nx,ny,nz)     ! Perturbation potential temperature (K)
  REAL :: pprt  (nx,ny,nz)     ! Perturbation pressure (Pascal)
  REAL :: qvprt (nx,ny,nz)     ! Perturbation wv mixing ration (kg/kg)
!
  REAL :: ubar  (nx,ny,nz)     ! Base state u-velocity (m/s)
  REAL :: vbar  (nx,ny,nz)     ! Base state v-velocity (m/s)
  REAL :: ptbar (nx,ny,nz)     ! Base state potential temperature (K)
  REAL :: pbar  (nx,ny,nz)     ! Base state pressure (Pascal)
  REAL :: rhobar(nx,ny,nz)     ! Base state air density (kg/m**3)
  REAL :: qvbar (nx,ny,nz)     ! Base state water vapor specific humidity
                               ! (kg/kg)
!
  REAL :: x     (nx)           ! The x-coord. of the physical and
                               ! computational grid. Defined at u-point.
  REAL :: y     (ny)           ! The y-coord. of the physical and
                               ! computational grid. Defined at v-point.
  REAL :: z     (nz)           ! The z-coord. of the computational grid.
                               ! Defined at w-point on the staggered grid.
  REAL :: zp    (nx,ny,nz)     ! The physical height coordinate defined at
                               ! w-point of the staggered grid.
!
  REAL :: hterain(nx,ny)       ! Terrain height.
!
  REAL :: roufns (nx,ny)       ! Surface roughness
!
  REAL :: tsfc(nx,ny)          ! Ground sfc. temperature (K)
!
  REAL :: raing(nx,ny)         ! Grid supersaturation rain (mm)
  REAL :: rainc(nx,ny)         ! Cumulus convective rain (mm)

  REAL :: time
  INTEGER :: arpsnn_opt
!
!-----------------------------------------------------------------------
!
!  Computed variables
!
!-----------------------------------------------------------------------
!
  REAL :: xs(nx)      ! x location of scalar points
  REAL :: ys(ny)      ! y location of scalar points
!
!-----------------------------------------------------------------------
!
!  Extracted sounding variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: nlevs
  INTEGER :: nzmax
  PARAMETER (nzmax=250)
!
  REAL :: su(nzmax)
  REAL :: sv(nzmax)
  REAL :: stheta(nzmax)
  REAL :: sqv(nzmax)
  REAL :: spres(nzmax)
  REAL :: stemp(nzmax)
  REAL :: sdewp(nzmax)
  REAL :: sdrct(nzmax)
  REAL :: ssped(nzmax)
  REAL :: shght(nzmax)
  REAL :: srain
!
!-----------------------------------------------------------------------
!
!  Work Arrays
!
!-----------------------------------------------------------------------
!
  REAL :: tem1(nx,ny,nz)
!
!-----------------------------------------------------------------------
!
!  Misc. internal variables
!
!------------------------------------------------------------------------
!
  INTEGER :: islat      ! Integer latitude of the sfc station
  INTEGER :: islon      ! Integer longitude of the sfc station
  INTEGER :: istnm      ! Station number of the sfc station
  INTEGER :: iselev_act ! Actual station elevation in integer format
  INTEGER :: i,j,k      ! Index variables
  INTEGER :: ireturn    ! Return status for whether stn is in the grid
  INTEGER :: LEN
!
  CHARACTER (LEN=256) :: dumpfile ! Sounding data history dump file name
  CHARACTER (LEN=256) :: sfclist  ! File to read the sfc station locations from
  CHARACTER (LEN=256) :: line     ! Temporary variable to read lines from files
  CHARACTER (LEN=256) :: sfcdomlist
  CHARACTER (LEN=256) :: sfcobsdir, sfcpre, sfcpost
  CHARACTER (LEN=256) :: mesoobsdir, mesopre, mesopost

  REAL :: press,presl       ! Pressures at sfc and level 2
  REAL :: hts,htl           ! Heights of sfc and level 2
  REAL :: hto               ! Observational height
  REAL :: temps,templ       ! Temperatures at ground skin and level 2
  REAL :: dewps,dewpl       ! Dew points at sfc and level 2
  REAL :: dirs,dirl         ! Wind direction at sfc and level 2
  REAL :: spds,spdl         ! Wind speeds at sfc and level 2
  REAL :: rough             ! Sfc roughness length
  REAL :: zeta              ! Monin-Obhuhov stability parameter
  REAL :: mixlength         ! Mixing length (constant in sfc. layer)
  REAL :: fricvel           ! Friction velocity
  REAL :: thetastar         ! Monin-Obhukov parameter
  REAL :: lapse             ! Low-level lapse rate
!
  INTEGER :: arps_drct      ! Estimated sfc wind direction
  REAL :: arps_spd          ! Estimated wind speed at 10 meters (m/s)
  REAL :: arps_temp         ! Estimated temperature at 2 meters (C)
  REAL :: arps_dewp         ! Estimated sfc dewpoint temperature (C)
  REAL :: arps_pres         ! Estimated sfc pressure (mb)
  REAL :: arps_theta        ! Estimated potential temp near 2 m (C)
  REAL :: arps_rain         ! Estimated precipitation (mm)

  INTEGER :: nhisfile
  REAL :: model_data(sfcmax,nhisfile,5)   ! Array of model obs
  REAL :: obsrv_data(sfcmax,nhisfile,5)   ! Array of real obs

!----------------------------------------------------------------------
!
!  Caren's NN variables
!
!----------------------------------------------------------------------
  REAL :: nn_input(2)       ! input for Caren's NN
  REAL :: nn_output         ! Caren's NN-corrected temp (C)
  CHARACTER (LEN=41) :: nn_loc
  LOGICAL :: is_there

!-----------------------------------------------------------------------
!
!  SAO obs. variables
!
!-----------------------------------------------------------------------
!
  REAL :: ob_lat,ob_lon,ob_elev
  REAL :: ob_t,ob_td
  REAL :: ob_dd,ob_ff
  REAL :: ob_ddg,ob_ffg
  REAL :: ob_pstn,ob_pmsl,ob_alt
  REAL :: ob_ceil,ob_lowcld,ob_cover
  REAL :: ob_vis,ob_rad,badflag
  CHARACTER (LEN=8) :: obstype
  CHARACTER (LEN=8) :: ob_wx
  CHARACTER (LEN=24) :: atime
  CHARACTER (LEN=5) :: ob_stn
  CHARACTER (LEN=256) :: obs_file
  CHARACTER (LEN=256) :: obsinputfl
  CHARACTER (LEN=8) :: ob_type

  REAL :: the_year
  INTEGER :: ob_kloud,ob_idp3
  INTEGER :: ob_time
  INTEGER :: valid_hour, valid_day, valid_month, valid_year
  INTEGER :: days_fwd
  INTEGER :: febr, next_day, month_day
  INTEGER :: the_int
  INTEGER :: ios,ii
  CHARACTER (LEN=4) :: year_ch
  CHARACTER (LEN=2) :: day_ch, hour_ch, month_ch, min_ch
  CHARACTER (LEN=100) :: the_line
  CHARACTER (LEN=100) :: the_char
  CHARACTER (LEN=14) :: arpsvarid
  CHARACTER (LEN=24) :: arpsname
  CHARACTER (LEN=13) :: obsvarid
  CHARACTER (LEN=25) :: obsname
  CHARACTER (LEN=256) :: hdfname
  CHARACTER (LEN=256) :: blackfile
  INTEGER :: sdirnam
  INTEGER :: lfnam
  INTEGER :: filepre
  INTEGER :: filepost

  INTEGER :: iyr,imon,idy,ihr,imin,isec,abstsec
  INTEGER :: filetime
  INTEGER :: sd_id
  INTEGER :: istat
  INTEGER :: nf

  REAL :: ddrot
  INTEGER :: n_meso_g,                                                  &
      n_meso_pos,n_sao_g,n_sao_pos_g,n_sao_b,n_sao_pos_b,n_obs_g,       &
      n_obs_pos_g,n_obs_b,n_obs_pos_b
  INTEGER :: nxt,iob,ivar,istatus
integer :: kwt

  REAL :: latsta(mx_sng,ntime),lonsta(mx_sng,ntime)
  REAL :: elevsta(mx_sng,ntime)
  REAL :: store_hgt(mx_sng,5,ntime)
  REAL :: obsrd(mx_sng,nvar_sng,ntime)
  CHARACTER (LEN=5) :: stn(mx_sng,ntime)
  CHARACTER (LEN=8) :: wx(mx_sng,ntime)
  CHARACTER (LEN=1) :: store_emv(mx_sng,5,ntime)
  CHARACTER (LEN=4) :: store_amt(mx_sng,5,ntime)
  INTEGER :: kloud(mx_sng,ntime),idp3(mx_sng,ntime)
  CHARACTER (LEN=8) :: chsrc(mx_sng,ntime)
  INTEGER :: obstime(mx_sng,ntime)

  integer :: keep
  integer :: nsfc

  REAL,EXTERNAL :: stprtopmsl

!
!-----------------------------------------------------------------------
!
!  Map projection variables
!
!-----------------------------------------------------------------------
!
  REAL :: xctr,yctr
  REAL :: xll,yll
  REAL :: latnot(2)
  REAL :: xmap(nx)
  REAL :: ymap(ny)
  REAL :: latgr(nx,ny)
  REAL :: longr(nx,ny)
  INTEGER :: nxlg, nylg
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  IF (nz > nzmax) THEN
    PRINT *,'Reset nzmax to greater than:',nz
    STOP
  END IF

!
!-----------------------------------------------------------------------
!
!  Calculate scalar locations and set up the map projection and
!  grid parameters
!
!-----------------------------------------------------------------------
!
  DO i=1,nx-1
    xs(i)=0.5*(x(i)+x(i+1))
  END DO
  DO j=1,ny-1
    ys(j)=0.5*(y(j)+y(j+1))
  END DO
  dx=x(2)-x(1)
  dy=y(2)-y(1)
!
  latnot(1)=trulat1
  latnot(2)=trulat2
  CALL setmapr(mapproj,sclfct,latnot,trulon)
  CALL lltoxy(1,1,ctrlat,ctrlon,xctr,yctr)
! xll=xctr-(0.5*(nx-3)*dx)
! yll=yctr-(0.5*(ny-3)*dy)
  nxlg = (nx-3)*nproc_x+3
  nylg = (ny-3)*nproc_y+3
  xll=xctr-(0.5*(nxlg-3)*dx)
  yll=yctr-(0.5*(nylg-3)*dy)

  DO i=1,nx-1
    xmap(i)=xll+xs(i)
  END DO
  xmap(nx)=2.*xmap(nx-1)-xmap(nx-2)
  DO j=1,ny-1
    ymap(j)=yll+ys(j)
  END DO
  ymap(ny)=2.*ymap(ny-1)-ymap(ny-2)
  CALL xytoll(nx,ny,xmap,ymap,latgr,longr)
!
!-----------------------------------------------------------------------
!
!  Find location in ARPS grid of all the stations in sfclist, then
!  rewrite only the ones that are in the grid to common arrays
!  so as to not make the program read all of the stations from
!  the original surface stn list and check their locations with each
!  history dump.
!
!-----------------------------------------------------------------------
!

  IF (readstns .eqv. .false. ) THEN

    IF (myproc == 0) THEN

      OPEN(UNIT=1,FILE=sfclist,STATUS='old')

      LEN =80
      CALL slength(sfclist,LEN)
      WRITE(6,*) 'Reading surface stations from file: ',                  &
                sfclist(1:LEN)


      READ(1,110,END=30) line ! read header lines and discard them..

      110    FORMAT(a80)
!     115      FORMAT(a3,3X,f5.2,2X,f7.2,1X,i4)
      115      FORMAT(a4,2X,f5.2,2X,f7.2,1X,i4)

      DO sfcstn=1,sfcmax
        READ(1,110,END=30) line
        READ(line,115) sfcstid(sfcstn),sfclat(sfcstn),sfclon(sfcstn),iselev_act
        sfcelev_act(sfcstn)=iselev_act/1.
      END DO
30    CONTINUE
      CLOSE(1)
      nsfc = sfcstn - 1

      IF (nsfc == sfcmax ) THEN
        WRITE(6,*)  &
           "WARNING:  Not all stations have been read.  Increase 'sfcmax'."
      END IF

    END IF

    CALL mpupdatei(nsfc,1)
    CALL mpupdatec(sfcstid,4*nsfc)
    CALL mpupdater(sfclat,nsfc)
    CALL mpupdater(sfclon,nsfc)
    CALL mpupdater(sfcelev_act,nsfc)

    istnm=0 !  Assumes a station number is unavailable

    sfcstn = 0
    sfcstn_master = 0

!
!   We'll print out what processors get what stations, however, due to race
!   issues with printing to STDOUT from different processors, some of the
!   output may never make it.
!

    DO i=1,nsfc
      IF (mp_opt > 0) THEN
        CALL flush(6)
        CALL mpbarrier
      END IF
      sfcstn=sfcstn+1

!  If we have a user list, use it!

      keep = 1
      IF ( nsfcuse .ne. 0 ) THEN
        keep = 0
        DO k=1,nsfcuse
          IF ( sfcstid(i) == sfcuse(k) ) keep = 1
        end do
      END IF
      IF ( keep == 0 ) THEN
        sfcstn = sfcstn -1
        CYCLE          ! this is a silly statement.  I prefer GOTOs!
      END IF

      sfcstid(sfcstn) = sfcstid(i)
      sfclat(sfcstn) = sfclat(i)
      sfclon(sfcstn) = sfclon(i)
      sfcelev_act(sfcstn) = sfcelev_act(i)

      IF ( myproc == 0 .AND. keep == 1) THEN
        sfcstn_master = sfcstn_master + 1
!       WRITE(6,*) 'sfcstn_master now ',sfcstn_master
        sfcstid_master(sfcstn_master) = sfcstid(i)
      END IF

      CALL lltoxy(1,1,sfclat(sfcstn),sfclon(sfcstn),                    &
                       sfcxpt(sfcstn),sfcypt(sfcstn))

      sfcxpt(sfcstn)=sfcxpt(sfcstn)-xll
      sfcypt(sfcstn)=sfcypt(sfcstn)-yll

      CALL findlc2(nx,ny,xs,ys,sfcxpt(sfcstn),sfcypt(sfcstn),           &
                 sfcipt(sfcstn),sfcjpt(sfcstn),ireturn)

      IF (mp_opt > 0) THEN
      END IF
      IF (ireturn < 0) THEN  ! stn is outside the grid
        sfcstn=sfcstn-1
      ELSE  ! stn is inside the grid
!   Tell the user what we know
        IF (mp_opt == 0) THEN
          WRITE(6,*) 'Station ',sfcstid(sfcstn)
        ELSE
          WRITE(6,*) 'Station ',sfcstid(sfcstn),' handled by processor ',myproc,i
!         WRITE(myproc+100,*) 'Station ',sfcstid(sfcstn),' handled by processor ',myproc
        ENDIF
      END IF

    END DO

  END IF

!
!-----------------------------------------------------------------------
!
!  Use Monin-Obhukov theory to estimate the surface layer profiles
!  of wind and temperature in order to get estimates of the model's
!  values that should correspond to those observations taken by the
!  mesonet.  If sufficient vertical resolution is available (i.e.
!  model levels surround the observation heights), perform a simpler
!  interpolation to get values of the variables.
!
!-----------------------------------------------------------------------
!

  readstns=.true.

  DO i=1,sfcstn ! do for each station 'i' ...
!
!-----------------------------------------------------------------------
!
!  Interpolate (in the horizontal) for the whole vertical column.
!
!-----------------------------------------------------------------------
!
    CALL colintb(nx,ny,nz,nzmax,                                        &
                xs,ys,zp,sfcxpt(i),sfcypt(i),sfcipt(i),sfcjpt(i),       &
                uprt, vprt, ptprt, pprt, qvprt,                         &
                ubar, vbar, ptbar, pbar, qvbar,                         &
                raing,rainc,srain,                                      &
                su,sv,stheta,spres,shght,sqv,                           &
                tem1,nlevs)
!
!-----------------------------------------------------------------------
!
!  Convert sounding to desired units/quantities (winds: m/s  theta: K
!
!-----------------------------------------------------------------------
!
    CALL cnvsnd(su,sv,stheta,spres,sqv,sfclon(i),                       &
                sdrct,ssped,stemp,sdewp,nlevs)
!
!-----------------------------------------------------------------------
!
!  Solve for the temperature and wind at observation levels using
!  Monin-Obukhov similarity theory.  Recall that levels 1 & 2
!  mirror around the surface.
!
!-----------------------------------------------------------------------
!
    press=EXP((LOG(spres(1))+LOG(spres(2)))/2.)
    presl=spres(2)

    hts=(shght(1)+shght(2))/2.
    htl=shght(2)-hts ! change heights to heights above the sfc
    hts=0.0          ! reset the height of the sfc to 0

    temps=tsfc(sfcipt(i),sfcjpt(i))-273.16
    templ=stemp(2)

    dewps=sdewp(2)   ! assume sfc dewpoint = dewpoint(level 2)
    dewpl=sdewp(2)

    dirs=(sdrct(1)+sdrct(2))/2.
    dirl=sdrct(2)

    spds=0.          ! no-slip condition
    spdl=ssped(2)
!
!-----------------------------------------------------------------------
!
!  Get the value of zeta (the M-O stability parameter: zeta=z/L)
!  using the equations from Byun, 1990.
!
!-----------------------------------------------------------------------
!
    rough=roufns(sfcipt(i),sfcjpt(i))

    CALL mosolns(press,hts,temps,dewps,dirs,spds,                       &
        presl,htl,templ,dewpl,dirl,spdl,rough,zeta)
!
!-----------------------------------------------------------------------
!
!  Solve for the M-O mixing length which in M-O theory
!  remains constant within the surface layer.  Recall
!  zeta=z/L where the z used is here set to be 1/2 of the
!  height htl (i.e. where the finite difference between
!  hts(=0.0) and htl should be most valid).  L is equivalent
!  to mixlength.
!
!-----------------------------------------------------------------------
!
    mixlength=htl/(zeta*2.0)
!
!----------------------------------------------------------------------
!
!  Solve for the friction velocity and theta* in Byun, 1990
!
!----------------------------------------------------------------------
!
    CALL ustar_thstar(zeta,spdl,htl,rough,templ,temps,                  &
        presl,press,fricvel,thetastar)
!
!----------------------------------------------------------------------
!
!  Interpolate the ARPS temps to 2 meters above the model
!  terrain and the winds to 10 meters above the model terrain.
!  Makes no correction (yet) for the differences between the
!  model terrain height and the actual station elevation.  If
!  sufficient vertical resolution is present, perform a simpler
!  linear interpolation to get temps and winds.  Set the dewpoint
!  at the observation level in the model to that at the surface,
!  which is taken to be the average of that at levels 1 & 2.  Set
!  the wind direction at the observational level to that at the
!  lowest model level.  (Either of these could be changed to
!  more accurate methods in the future.)
!
!----------------------------------------------------------------------
!
    lapse=-(templ-temps)/(htl-hts)

    CALL arps2obs(fricvel,thetastar,temps,press,presl,hts,htl,          &
        mixlength,lapse,rough,arps_spd,arps_temp)

    IF (dewps > arps_temp) THEN
      arps_dewp=arps_temp
    ELSE
      arps_dewp=dewps
    END IF

    arps_drct=nint(sdrct(2))

    arps_pres=press

    arps_rain=srain

    arps_theta=(arps_temp+273.16)*(1000./arps_pres)**0.286 - 273.16

    IF ((htl-hts) < 2.) THEN

      hto=hts+2.

      DO j=1,nlevs

        IF (shght(j) < hto .AND. shght(j+1) >= hto) THEN
          arps_temp=(hto-shght(j))/(shght(j+1)-shght(j)) *              &
              (stemp(j+1)-stemp(j)) + stemp(j)
        END IF

      END DO

    END IF

    IF ((htl-hts) < 10.) THEN

      hto=hts+10.

      DO j=1,nlevs

        IF (shght(j) < hto .AND. shght(j+1) >= hto) THEN
          arps_spd=(hto-shght(j))/(shght(j+1)-shght(j)) *            &
              (ssped(j+1)-ssped(j)) + ssped(j)
          arps_drct=sdrct(j)
        END IF

      END DO

    END IF

    hts=(shght(1)+shght(2))/2.
!
!---------------------------------------------------------------------
!
!         MAKE A CALL TO CAREN'S NN FOR TEMP. CORRECTION
!
!---------------------------------------------------------------------
!
! nn_input(1) = the forecast time in hours after init time
!   IMPORTANT: NN_INPUT(1) is assuming that the model run is 12Z!!!
! nn_input(2) = the arps forecast temp.

! first check to see if NN exists for the station...
    IF (arpsnn_opt.eq.1) THEN
    nn_loc = '/home/nmirsky/arpsverif/src_h2v/wts_K'
!
!  The state of this part of the code is unknown, as I've never used it.
!  There willprobably be consequences that I'm not following (including a
!  few lines below) that are a consequence of going to a 4 character station
!  id.

!   WRITE(nn_loc(38:40),'(a3)') sfcstid(i)(1:3)
    WRITE(nn_loc(38:41),'(a4)') sfcstid(i)(1:4)
    INQUIRE(FILE=nn_loc, EXIST=is_there)

    IF (is_there) THEN
      nn_input(1) = (INT(time)/3600.0) + 6.0
      nn_input(2) = arps_temp
      CALL arps_nn('K'//sfcstid(i)(1:3),nn_input, nn_output)
      nn_output = (anint(nn_output*10.0))/10.0
    ELSE
      nn_output = -999.0
    END IF

    END IF

!--------------------------------------------------------------------
!
!    Set up the default values if surf observatiopns are not used for
!    verification.
!
!--------------------------------------------------------------------

     ob_t = -999.0
     ob_td = -999.0
     ob_dd = -99
     ob_ff = -999.0
     ob_ddg = -999.0
     ob_ffg = -999.0
     ob_pstn = -999.0
     ob_pmsl = -999.0
     ob_alt = -999.0
     ob_kloud = -99
     ob_ceil = -999.0
     ob_lowcld = -999.0
     ob_cover = -999.0
     ob_vis = -999.0
     ob_rad = -999.0
     ob_idp3 = -999

!
! Convert the calculated ARPS observation values and place into the
! model_data array.
!

    model_data(i,nf,1)=(arps_temp*(9./5.))+32.0
    model_data(i,nf,2)=(arps_dewp*(9./5.))+32.0
    model_data(i,nf,3)=arps_drct
    model_data(i,nf,4)=arps_spd
    model_data(i,nf,5)=stprtopmsl(arps_pres,(arps_temp+273.16),         &
                       sfcelev_act(i))

  END DO

  IF (myproc == 0) THEN
!
! Read the surface observations and place into the obsrv_data array
! for future processing.
! The SAO data is read in using the read_surface_obs subroutine in
! ADAS, so the SAO data format needs to be in the standard LAPS
! format used by ADAS.
!

    CALL ctim2abss(year,month,day,hour,minute,second,abstsec)
    filetime=abstsec+int(time)
    CALL abss2ctim(filetime,iyr,imon,idy,ihr,imin,isec)
    WRITE (year_ch, 50) iyr
    WRITE (month_ch, 55) imon
    WRITE (day_ch, 55) idy
    WRITE (hour_ch, 55) ihr
    WRITE (min_ch, 55) imin
50  FORMAT (i4.4)
55  FORMAT (i2.2)

    filepre=LEN_TRIM(sfcpre)
    filepost=LEN_TRIM(sfcpost)
    WRITE(obs_file, 60) year_ch,month_ch,day_ch,hour_ch,min_ch
60  FORMAT(a4,a2,a2,a2,a2,'.')
    sdirnam=LEN_TRIM(sfcobsdir)
    lfnam=LEN_TRIM(obs_file)
    obsinputfl= sfcobsdir(1:sdirnam)//                                  &
      sfcpre(1:filepre)//obs_file(1:lfnam)//sfcpost(1:filepost)

    nxt=1
    CALL read_surface_obs(obsinputfl,blackfile,mx_sng,nxt,atime,n_meso_g, &
      n_meso_pos,n_sao_g,n_sao_pos_g,n_sao_b,n_sao_pos_b,n_obs_g,       &
      n_obs_pos_g,n_obs_b,n_obs_pos_b,stn(1,1),chsrc(1,1),              &
      latsta(1,1),lonsta(1,1),elevsta(1,1),wx(1,1),                     &
      obsrd(1,1,1),obsrd(1,2,1),obsrd(1,3,1),obsrd(1,4,1),              &
      obsrd(1,5,1),obsrd(1,6,1),obsrd(1,7,1),obsrd(1,8,1),              &
      obsrd(1,9,1),  kloud(1,1),obsrd(1,10,1),obsrd(1,11,1),            &
      obsrd(1,12,1),obsrd(1,13,1),idp3(1,1),store_emv(1,1,1),           &
      store_amt(1,1,1),store_hgt(1,1,1),obsrd(1,14,1),                  &
      obstime(1,1),istatus)

!
! Same drill, except use the mesonet file.
!

    IF (mesoobsdir /= "NULL" ) THEN
      filepre=LEN_TRIM(mesopre)
      filepost=LEN_TRIM(mesopost)
      sdirnam=LEN_TRIM(mesoobsdir)
      obsinputfl=mesoobsdir(1:sdirnam)//                                &
        mesopre(1:filepre)//obs_file(1:lfnam)//mesopost(1:filepost)

      nxt=n_obs_b

!
!   Only variables that are used will have unique names for mesonet obs.
!   If more are used in the fugure, this code will have to be modified!
!

      CALL read_surface_obs(obsinputfl,blackfile,mx_sng,nxt,atime,n_meso_g, &
        n_meso_pos,n_sao_g,n_sao_pos_g,n_sao_b,n_sao_pos_b,n_obs_g,       &
        n_obs_pos_g,n_obs_b,n_obs_pos_b,stn(1,1),chsrc(1,1),              &
        latsta(1,1),lonsta(1,1),elevsta(1,1),wx(1,1),                     &
        obsrd(1,1,1),obsrd(1,2,1),obsrd(1,3,1),obsrd(1,4,1),              &
        obsrd(1,5,1),obsrd(1,6,1),obsrd(1,7,1),obsrd(1,8,1),              &
        obsrd(1,9,1),  kloud(1,1),obsrd(1,10,1),obsrd(1,11,1),            &
        obsrd(1,12,1),obsrd(1,13,1),idp3(1,1),store_emv(1,1,1),           &
        store_amt(1,1,1),store_hgt(1,1,1),obsrd(1,14,1),                  &
        obstime(1,1),istatus)

        n_obs_b = n_obs_b + nxt
    END IF

  END IF

  CALL mpupdatei(n_obs_b,1)
  CALL mpupdatec(stn(1,1),5*n_obs_b)
  CALL mpupdater(obsrd(1,1,1),n_obs_b)
  CALL mpupdater(obsrd(1,2,1),n_obs_b)
  CALL mpupdater(obsrd(1,3,1),n_obs_b)
  CALL mpupdater(obsrd(1,4,1),n_obs_b)
  CALL mpupdater(obsrd(1,8,1),n_obs_b)

  DO j=1,sfcstn
   DO k=1,n_obs_b
   IF (sfcstid(j) == stn(k,1)) THEN
    obsrv_data(j,nf,1)=obsrd(k,1,1)
    obsrv_data(j,nf,2)=obsrd(k,2,1)
    obsrv_data(j,nf,3)=obsrd(k,3,1)
    obsrv_data(j,nf,4)=obsrd(k,4,1)
    obsrv_data(j,nf,5)=obsrd(k,8,1)
   END IF
   END DO
  END DO

!
!  Fix missing data.  LAPS files use -100, while the rest of this software
!  uses -99.9.
!
!  Also check for a bogus wind direction.  I've seen one, which was probably
!  a coding error, or network noise when the data arrived.
!

  DO k=1,5
    DO i=1,sfcstn
      IF (obsrv_data(i,nf,k) < -90.0) obsrv_data(i,nf,k) = -99.9
      IF (k == 3 ) THEN
        IF (obsrv_data(i,nf,k) < 0.0 .OR. obsrv_data(i,nf,k) > 360.0)      &
          obsrv_data(i,nf,k) = -99.9
      END IF
    END DO
  END DO

  RETURN
!
END SUBROUTINE sfcdump


!
!
!#################################################################
!#################################################################
!######                                                     ######
!######                SUBROUTINE PRECDUMP                  ######
!######                                                     ######
!######             Copyright (c) 1995-1996                 ######
!######    Center for Analysis and Prediction of Storms     ######
!######    University of Oklahoma.  All rights reserved.    ######
!######                                                     ######
!#################################################################
!#################################################################
!

SUBROUTINE precdump(nx,ny,nz,                                           &
           uprt,vprt,ptprt,pprt,qvprt,                                  &
           ubar,vbar,ptbar,                                             &
           pbar,rhobar,qvbar,                                           &
           x,y,z,zp,hterain,                                            &
           roufns,tsfc,                                                 &
           dumpfile,needsfcstns,                                        &
           raing,rainc,                                                 &
           xs,ys,xmap,ymap,                                             &
           latgr,longr,tem1,                                            &
           sfclist,sfcdomlist,time)
!
!----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Coordinate the verification precipitation dumps, which are
!  used in program PRECANAL.
!
!----------------------------------------------------------------------
!
!  AUTHOR:  John Mewes
!
!  MODIFICATION HISTORY:
!
!  09/15/97  Changed subroutine name from MESODUMP to SFCDUMP in
!            order to bring the subroutine name into a consistent
!            format (to go with SNDDUMP and PRODUMP).
!
!  09/16/97  Major code revamping.  Changed variable names and put
!            station information arrays into a common block in order
!            to avoid having to reread information with each new
!            data dump (adds vericst.inc to the verification
!            software).
!
!  06/30/98  Modified to dump only precip data and information
!            on the history dump grid (lat/lon of corners,
!            map projection, etc.) (Eric Kemp, Project COMET-Tinker)
!
!  19 May 2000:  Eric Kemp
!                File name for list of stations in the domain
!                is now passed to the subroutine, along with
!                time from initialization.  Also, valid
!                time and forecast domain info are now included
!                in all output files.
!
!-----------------------------------------------------------------------
!
!  Variables to be read in:
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
  INTEGER :: nx,ny,nz          ! Number of grid points in 3 directions
!
  INTEGER :: nzsoil
!
  LOGICAL :: needsfcstns       ! Whether not need to create station file
!
  REAL :: uprt  (nx,ny,nz)     ! Perturbation u velocity (m/s)
  REAL :: vprt  (nx,ny,nz)     ! Perturbation v velocity (m/s)
  REAL :: ptprt (nx,ny,nz)     ! Perturbation potential temperature (K)
  REAL :: pprt  (nx,ny,nz)     ! Perturbation pressure (Pascal)
  REAL :: qvprt (nx,ny,nz)     ! Perturbation wv mixing ration (kg/kg)
!
  REAL :: ubar  (nx,ny,nz)     ! Base state u-velocity (m/s)
  REAL :: vbar  (nx,ny,nz)     ! Base state v-velocity (m/s)
  REAL :: ptbar (nx,ny,nz)     ! Base state potential temperature (K)
  REAL :: pbar  (nx,ny,nz)     ! Base state pressure (Pascal)
  REAL :: rhobar(nx,ny,nz)     ! Base state air density (kg/m**3)
  REAL :: qvbar (nx,ny,nz)     ! Base state water vapor specific humidity
                               ! (kg/kg)
!
  REAL :: x     (nx)           ! The x-coord. of the physical and
                               ! computational grid. Defined at u-point.
  REAL :: y     (ny)           ! The y-coord. of the physical and
                               ! computational grid. Defined at v-point.
  REAL :: z     (nz)           ! The z-coord. of the computational grid.
                               ! Defined at w-point on the staggered grid.
  REAL :: zp    (nx,ny,nz)     ! The physical height coordinate defined at
                               ! w-point of the staggered grid.
!
  REAL :: hterain(nx,ny)       ! Terrain height.
!
  REAL :: roufns (nx,ny)       ! Surface roughness
!
  REAL :: tsfc(nx,ny)          ! Ground sfc. temperature (K)
!
  REAL :: raing(nx,ny)         ! Grid supersaturation rain (mm)
  REAL :: rainc(nx,ny)         ! Cumulus convective rain (mm)

  REAL :: time
!
!-----------------------------------------------------------------------
!
!  Computed variables
!
!-----------------------------------------------------------------------
!
  REAL :: xs(nx)      ! x location of scalar points
  REAL :: ys(ny)      ! y location of scalar points
!
!-----------------------------------------------------------------------
!
!  Extracted sounding variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: nlevs
  INTEGER :: nzmax
  PARAMETER (nzmax=250)
!
  REAL :: su(nzmax)
  REAL :: sv(nzmax)
  REAL :: stheta(nzmax)
  REAL :: sqv(nzmax)
  REAL :: spres(nzmax)
  REAL :: stemp(nzmax)
  REAL :: sdewp(nzmax)
  REAL :: sdrct(nzmax)
  REAL :: ssped(nzmax)
  REAL :: shght(nzmax)
  REAL :: srain
!
!-----------------------------------------------------------------------
!
!  Work Arrays
!
!-----------------------------------------------------------------------
!
  REAL :: tem1(nx,ny,nz)
!
!-----------------------------------------------------------------------
!
!  Misc. internal variables
!
!------------------------------------------------------------------------
!
  INTEGER :: islat      ! Integer latitude of the sfc station
  INTEGER :: islon      ! Integer longitude of the sfc station
  INTEGER :: istnm      ! Station number of the sfc station
  INTEGER :: iselev_act ! Actual station elevation in integer format
  INTEGER :: i,j        ! Index variables
  INTEGER :: ireturn    ! Return status for whether stn is in the grid
  INTEGER :: LEN
!
  CHARACTER (LEN=256) :: dumpfile ! Sounding data history dump file name
  CHARACTER (LEN=256) :: sfclist  ! File to read the sfc station locations from
  CHARACTER (LEN=256) :: line     ! Temporary variable to read lines from files
  CHARACTER (LEN=256) :: sfcdomlist
!
  REAL :: press,presl       ! Pressures at sfc and level 2
  REAL :: hts,htl           ! Heights of sfc and level 2
  REAL :: hto               ! Observational height
  REAL :: temps,templ       ! Temperatures at ground skin and level 2
  REAL :: dewps,dewpl       ! Dew points at sfc and level 2
  REAL :: dirs,dirl         ! Wind direction at sfc and level 2
  REAL :: spds,spdl         ! Wind speeds at sfc and level 2
  REAL :: rough             ! Sfc roughness length
  REAL :: zeta              ! Monin-Obhuhov stability parameter
  REAL :: mixlength         ! Mixing length (constant in sfc. layer)
  REAL :: fricvel           ! Friction velocity
  REAL :: thetastar         ! Monin-Obhukov parameter
  REAL :: lapse             ! Low-level lapse rate
!
  INTEGER :: arps_drct      ! Estimated sfc wind direction
  REAL :: arps_spd          ! Estimated wind speed at 10 meters (m/s)
  REAL :: arps_temp         ! Estimated temperature at 2 meters (C)
  REAL :: arps_dewp         ! Estimated sfc dewpoint temperature (C)
  REAL :: arps_pres         ! Estimated sfc pressure (mb)
  REAL :: arps_theta        ! Estimated potential temp near 2 m (C)
  REAL :: arps_rain         ! Estimated precipitation (mm)
!
!-----------------------------------------------------------------------
!
!  Map projection variables
!
!-----------------------------------------------------------------------
!
  REAL :: xctr,yctr
  REAL :: xll,yll
  REAL :: latnot(2)
  REAL :: xmap(nx)
  REAL :: ymap(ny)
  REAL :: latgr(nx,ny)
  REAL :: longr(nx,ny)

  REAL :: rloc1(2),rloc2(2),rloc3(2),rloc4(2)

  REAL :: xc(nx,ny,nz),yc(nx,ny,nz)
  REAL :: swx,swy
  INTEGER :: k
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'phycst.inc'
  INCLUDE 'globcst.inc'
  INCLUDE 'bndry.inc'
  INCLUDE 'exbc.inc'
  INCLUDE 'vericst.inc'
  INCLUDE 'grid.inc'
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  IF (nz > nzmax) THEN
    PRINT *,'Reset nzmax to greater than:',nz
    STOP
  END IF
!
!-----------------------------------------------------------------------
!
!  Calculate scalar locations and set up the map projection and
!  grid parameters
!
!-----------------------------------------------------------------------
!
  DO i=1,nx-1
    xs(i)=0.5*(x(i)+x(i+1))
  END DO
  DO j=1,ny-1
    ys(j)=0.5*(y(j)+y(j+1))
  END DO
  dx=x(2)-x(1)
  dy=y(2)-y(1)
!

  latnot(1)=trulat1
  latnot(2)=trulat2
  CALL setmapr(mapproj,sclfct,latnot,trulon)
  CALL lltoxy(1,1,ctrlat,ctrlon,xctr,yctr)
  xll=xctr-(0.5*(nx-3)*dx)
  yll=yctr-(0.5*(ny-3)*dy)

  DO i=1,nx-1
    xmap(i)=xll+xs(i)
  END DO
  xmap(nx)=2.*xmap(nx-1)-xmap(nx-2)
  DO j=1,ny-1
    ymap(j)=yll+ys(j)
  END DO
  ymap(ny)=2.*ymap(ny-1)-ymap(ny-2)
  CALL xytoll(nx,ny,xmap,ymap,latgr,longr)
!
!-----------------------------------------------------------------------
!
!  Find location in ARPS grid of all the stations in sfclist, then
!  rewrite only the ones that are in the grid to common arrays
!  so as to not make the program read all of the stations from
!  the original surface stn list and check their locations with each
!  history dump.
!
!-----------------------------------------------------------------------
!
  OPEN(UNIT=1,FILE=sfclist,STATUS='old')
!  OPEN(UNIT=2,FILE=dumpfile,STATUS='unknown')
!
!  LEN =256
!  CALL slength( dumpfile, LEN)
!  WRITE(6,*) "Opening ",dumpfile(1:LEN)," for surface output"


  WRITE(2,990) year,'-',month,'-',day,'.',hour,':',minute,':',          &
               second,'+',INT(time)
  990   FORMAT(i4.4,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i6.6)

  LEN =80
  CALL slength(runname,LEN)
  WRITE(2,'(a)') runname(1:LEN)
  WRITE(2,'(i1.1)') nocmnt
  DO i = 1,nocmnt
    LEN =80
    CALL slength(cmnt(i),LEN)
    WRITE(2,'(a)') cmnt(i)(1:LEN)
  END DO

  rloc1(1) = latgr(1,1)
  rloc1(2) = longr(1,1)

  rloc2(1) = latgr(nx-1,1)
  rloc2(2) = longr(nx-1,1)

  rloc3(1) = latgr(1,ny-1)
  rloc3(2) = longr(1,ny-1)

  rloc4(1) = latgr(nx-1,ny-1)
  rloc4(2) = longr(nx-1,ny-1)


  WRITE(2,*) ctrlat
  WRITE(2,*) ctrlon

  WRITE(2,*) rloc1(1)
  WRITE(2,*) rloc1(2)

  WRITE(2,*) rloc2(1)
  WRITE(2,*) rloc2(2)

  WRITE(2,*) rloc3(1)
  WRITE(2,*) rloc3(2)

  WRITE(2,*) rloc4(1)
  WRITE(2,*) rloc4(2)

  WRITE(2,*) trulat1
  WRITE(2,*) trulat2

  WRITE(2,*) trulon
  WRITE(2,*) sclfct

  WRITE(2,*) mapproj

  105  FORMAT(a8,a16,a20)
!
  IF (needsfcstns) THEN
!
    needsfcstns = .false.
    sfcstn=0
    sfcstn=0
!
    LEN =80
    CALL slength(sfclist,LEN)
    WRITE(6,*) 'Reading surface stations from file: ',                  &
                sfclist(1:LEN)
!
    OPEN(UNIT=3,FILE=sfcdomlist,STATUS='unknown')
    WRITE(3,990) year,'-',month,'-',day,'.',hour,':',minute,':',        &
               second,'+',INT(time)

    LEN =80
    CALL slength(runname,LEN)
    WRITE(3,'(a)') runname(1:LEN)
    WRITE(3,'(i1.1)') nocmnt
    DO i = 1,nocmnt
      LEN =80
      CALL slength(cmnt(i),LEN)
      WRITE(3,'(a)') cmnt(i)(1:LEN)
    END DO
!

    110    FORMAT(a80)
!
    READ(1,110,END=30) line ! read header lines and discard them..
    READ(1,110,END=30) line
    READ(1,110,END=30) line
    READ(1,110,END=30) line
!
    DO i=1,sfcmax
!
      115      FORMAT(a3,1X,f6.2,1X,f7.2,1X,i4)
!
      READ(1,110,END=30) line
!
      sfcstn=sfcstn+1
!
      READ(line,115) precstid(sfcstn),islat,islon,iselev_act

      istnm=0 !  Assumes a station number is unavailable
      sfclat(sfcstn)=islat/100000.
      sfclon(sfcstn)=islon/100000.
      sfcelev_act(sfcstn)=iselev_act/1.
!
      CALL lltoxy(1,1,sfclat(sfcstn),sfclon(sfcstn),                    &
                       sfcxpt(sfcstn),sfcypt(sfcstn))
!
      sfcxpt(sfcstn)=sfcxpt(sfcstn)-xll
      sfcypt(sfcstn)=sfcypt(sfcstn)-yll
!
      CALL findlc2(nx,ny,xs,ys,sfcxpt(sfcstn),sfcypt(sfcstn),           &
                 sfcipt(sfcstn),sfcjpt(sfcstn),ireturn)
!
      IF (ireturn < 0) THEN  ! stn is outside the grid
        sfcstn=sfcstn-1
      ELSE  ! stn is inside the grid
        WRITE(3,110) line ! write location data to a file
      END IF
!
    END DO
!
  END IF
!
  30   CONTINUE
!
  DO i=1,sfcstn ! do for each station 'i' ...

    dumpfile = './prec_dumps/'
    LEN =256
    CALL slength( dumpfile, LEN)
!   WRITE(dumpfile(LEN+1:LEN+8),'(a8)') sfcstid(i)(1:3)//'_prec'
    WRITE(dumpfile(LEN+1:LEN+8),'(a8)') sfcstid(i)(1:4)//'_prec'
    OPEN(UNIT=2,FILE=dumpfile,STATUS='unknown',POSITION='append')

    WRITE(2,991) year,'-',month,'-',day,'.',hour,':',minute,':',        &
               second,'+',INT(time)
    991   FORMAT(i4.4,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i6.6)

!
!-----------------------------------------------------------------------
!
!  Interpolate (in the horizontal) for the whole vertical column.
!
!-----------------------------------------------------------------------
!
    CALL colintb(nx,ny,nz,nzmax,                                        &
                xs,ys,zp,sfcxpt(i),sfcypt(i),sfcipt(i),sfcjpt(i),       &
                uprt, vprt, ptprt, pprt, qvprt,                         &
                ubar, vbar, ptbar, pbar, qvbar,                         &
                raing,rainc,srain,                                      &
                su,sv,stheta,spres,shght,sqv,                           &
                tem1,nlevs)
!
    arps_rain=srain


    120    FORMAT('Stn:',a16,' rain:',f6.2)
    WRITE(2,120) precstid(i),arps_rain
!
    CLOSE(2)
!
  END DO
!
  CLOSE(1)
!  CLOSE(2)
  CLOSE(3)
!
  RETURN
!
END SUBROUTINE precdump
!
!
!#################################################################
!#################################################################
!######                                                     ######
!######                  Function AINT2D                    ######
!######                                                     ######
!######                Copyright (c) 1996                   ######
!######    Center for Analysis and Prediction of Storms     ######
!######    University of Oklahoma.  All rights reserved.    ######
!######                                                     ######
!#################################################################
!#################################################################
!
!

  FUNCTION aint2dtmp(nx,ny,nz,a,im,jm,k,in,jn,w1,w2,w3,w4)
!
  IMPLICIT NONE
!
  REAL :: aint2dtmp
  INTEGER :: nx,ny,nz
  REAL :: a(nx,ny,nz)
  INTEGER :: im,jm,in,jn,k
  REAL :: w1,w2,w3,w4

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  aint2dtmp = w1*a(im,jm,k) + w2*a(in,jm,k) +                           &
              w3*a(in,jn,k) + w4*a(im,jn,k)

  RETURN

  END FUNCTION aint2dtmp
!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINE mosolns                   ######
!######                                                      ######
!######                Copyright (c) 1995-1996               ######
!######    Center for Analysis and Prediction of Storms      ######
!######    University of Oklahoma.  All rights reserved.     ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE mosolns(press,hts,temps,dewps,dirs,                          &
           spds,presl,htl,templ,dewpl,dirl,spdl,rough,zeta)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Uses the equations from Byun, 1990 (JAM) to solve for the value
!  of the stability parameter 'zeta' in Monin-Obukhov theory.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: John Mewes
!  September 1995.
!
!  MODIFICATION HISTORY:
!
!  November, 1995 (J. Mewes)
!  Cleaned up and commented.
!
!  09/16/97  Changed some variable names and recommented the code.
!            Removed function call to calculate theta (done explic-
!            itly now).
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
  REAL :: bh,bm        ! Constants (beta-h, beta-m)
  PARAMETER (bh=6.35,bm=4.7)
!
  REAL :: pr           ! Prandtl # for neutral stability
  PARAMETER (pr=0.74)
!
  REAL :: gm,gh        ! Constants (gamma-m, gamma-h)
  PARAMETER (gh=9.0,gm=15.0)
!
  REAL :: fctr         ! Conversion for using the Bulk Richardson #
  REAL :: deltau       ! Difference in wind spd between sfc and level 2
  REAL :: deltaz       ! Difference in height between sfc and level 2
  REAL :: deltatheta   ! Difference in theta between sfc and level 2
  REAL :: thetanot     ! Average theta between sfc and level 2
  REAL :: bulkrich     ! Bulk Richardson #, from Byun 1990
  REAL :: rough        ! Roughness length at station
  REAL :: zeta         ! M-O stability parameter
!
  REAL :: thetai,si,ti,qi,pi
                     ! Terms used in calculating 'zeta'
!
  REAL :: templ,dewpl,htl,presl,dirl,spdl,thetal
                     ! Values of the ARPS variables at stn at level 2
!
  REAL :: temps,dewps,hts,press,dirs,spds,thetas
                     ! Values of the ARPS variables at stn at the sfc
!
  REAL :: check        ! Used to check for stability
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  thetal=(templ+273.16)*(1000./presl)**.286
  thetas=(temps+273.16)*(1000./press)**.286
  thetanot=(thetas+thetal)/2
  deltatheta=thetal-thetas
!
  deltaz=htl-hts
!
  deltau=spdl       ! Since speed at sfc = 0.0
!
!----------------------------------------------------------------------
!
!  Calculate the Richardson # between the sfc and level 2.
!
!  Note:  In the future may want to evaluate stability based on
!         theta-e rather than theta.
!
!----------------------------------------------------------------------
!
  bulkrich = (9.81*(deltatheta*deltaz))/(thetanot*(deltau**2))
  bulkrich = AMAX1(bulkrich,-10.0) ! limit for unstable case
!
  fctr=deltaz/(deltaz-rough)*LOG(deltaz/rough)
!
!
  IF (bulkrich > 0.) THEN
!
!---------------------------------------------------------------------
!
!  A Richardson # > 0 indicates a stable atmosphere (i.e.
!  it is positive because the potential temperature is
!  greater at the higher level than at the lower level).
!
!---------------------------------------------------------------------
!
    zeta=(-(2.*bh*bulkrich-1.)-(1.+4.*(bh-bm)*bulkrich/pr)**0.5) /      &
         (2.*bh*(bm*bulkrich-1.))*fctr
!
  ELSE
!
    si=bulkrich/pr
    qi=(1./(gm**2.0)+3.*(gh/gm)*(si**2.0))/9.
    pi=(-2./(gm**3.0)+9./gm*(-gh/gm+3.)*(si**2.0))/54.
    check=qi**3.0-pi**2.0 ! Delineates the two cases for unstable
!
!--------------------------------------------------------------------
!
!  There are two cases for unstable (depending on how
!  unstable the atmosphere is):
!
!--------------------------------------------------------------------
!
    IF (check >= 0.) THEN
!
      thetai=ACOS(pi/SQRT(qi**3.0))
      zeta=(-2.*SQRT(qi)*COS(thetai/3.0)+1/(3.*gm))*fctr
!
    ELSE
!
      ti=(SQRT(-check)+ABS(pi))**(1./3.)
      zeta=(-(ti+qi/ti)+1./(3.*gm))*fctr
!
    END IF
!
  END IF
!
  RETURN
!
END SUBROUTINE mosolns
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE ustar_thstar               ######
!######                                                      ######
!######                Copyright (c) 1995-1996               ######
!######    Center for Analysis and Prediction of Storms      ######
!######    University of Oklahoma.  All rights reserved.     ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE ustar_thstar(zeta,spdl,z,zo,templ,temps,presl,press,         &
           fricvel,thetastar)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Uses the equations from Byun, 1990 (JAM) to solve for the value
!  of the friction velocity (u*) and the scaling temperature (theta*)
!  using the mixing length.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: John Mewes
!  September 1995.
!
!  MODIFICATION HISTORY:
!
!  November, 1995 (J. Mewes)
!  Cleaned up and commented.
!
!  09/16/97  Changed some variable names and recommented the code.
!            Removed function call to calculate theta (done explic-
!            itly now).  Changed subroutine name.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
  REAL :: bh,bm              ! Constants (beta-h, beta-m)
  PARAMETER (bh=6.35,bm=4.7)
!
  REAL :: pr                 ! Prandtl # for neutral stability
  PARAMETER (pr=0.74)
!
  REAL :: gm,gh,k            ! Constants (gamma-m, gamma-h)
  PARAMETER (gm=15.0,gh=9.0,k=0.375)
!
  REAL :: templ,temps        ! Model temps at level 2 and sfc
  REAL :: presl,press        ! Model pressure at level 2 and sfc
  REAL :: thetal,thetas      ! Model theta at level 2 and sfc
  REAL :: spdl               ! Model wind speed at level 2
  REAL :: z                  ! Height above sfc of level 2
  REAL :: zo                 ! Roughness length
  REAL :: fricvel            ! Friction velocity (solving for this!)
  REAL :: thetastar          ! Theta* (solving for this also!)
  REAL :: x,xo               ! Calculation terms
  REAL :: y,yo               ! Calculation terms
  REAL :: zeta,zetao         ! Stability parameters
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  zetao=zeta*zo/z
!
  thetas=(temps+273.16)*(1000./press)**.286
  thetal=(templ+273.16)*(1000./presl)**.286
!
  IF (zeta > 0.) THEN  ! zeta > 0 denotes stable conditions...
!
    fricvel=k*spdl/(LOG(z/zo)+bm*(zeta-zetao))
!
  ELSE   ! ...unstable conditions
!
    x=(1-gm*zeta)**(1./4.)
    xo=(1-gm*zo/z*zeta)**(1./4.)
    y=(1-gh*zeta)**(1./2.)
    yo=(1-gh*zo/z*zeta)**(1./2.)
    fricvel=k*spdl/(LOG(z/zo)-(2*LOG((1+x)/(1+xo))+                     &
        LOG((1+x**2)/(1+xo**2))-2*ATAN(x)+2*ATAN(xo)))
!
  END IF
!
  IF (zeta > 0.) THEN  ! again, stable conditions...
!
    thetastar=(k*(thetal-thetas)/pr)/(LOG(z/zo)+bh*                     &
        (zeta-zetao))
!
  ELSE  ! ...and unstable conditions.
!
    thetastar=(k*(thetal-thetas)/pr)/(LOG(z/zo)-2*                      &
        LOG((y+1)/(yo+1)))
!
  END IF
!
  RETURN
!
END SUBROUTINE ustar_thstar
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINE arps2obs                  ######
!######                                                      ######
!######                Copyright (c) 1995-1996               ######
!######    Center for Analysis and Prediction of Storms      ######
!######    University of Oklahoma.  All rights reserved.     ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE arps2obs(fricvel,thetastar,temps,press,presl,hts,htl,        &
           mixlength,lapse,zo,arps_spd,arps_temp)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Interpolates the ARPS data to the correct elevation above the
!  model terrain so as to correspond with the observation heights
!  above the surface.  Does not take into account the differences
!  between the model terrain and the actual elevations.  Uses Monin-
!  Obukhov theory.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: John Mewes
!  September 1995.
!
!  MODIFICATION HISTORY:
!
!  November, 1995 (J. Mewes)
!  Cleaned up and commented.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
  REAL :: bh,bm              ! Constants (beta-h, beta-m)
  PARAMETER (bh=6.35,bm=4.7)
!
  REAL :: pr                 ! Prandtl # for neutral stability
  PARAMETER (pr=0.74)
!
  REAL :: gm,gh,k            ! Constants (gamma-m, gamma-h)
  PARAMETER (gm=15.0,gh=9.0,k=0.375)
!

  REAL :: fricvel            ! Friction velocity
  REAL :: thetastar          ! Theta*
  REAL :: arps_spd           ! Wind speed at 10 meters
  REAL :: arps_temp          ! Temperature at 2 meters
  REAL :: temps              ! Surface temperature
  REAL :: hts,htl            ! Height of sfc and level 2
  REAL :: press,presl        ! Pressure at sfc and level 2
  REAL :: thetas             ! Surface potential temperature
  REAL :: dtheta             ! Change in theta from zo to obs level
  REAL :: zo                 ! Roughness length
  REAL :: x,xo               ! Calculation terms
  REAL :: y,yo               ! Calculation terms
  REAL :: mixlength          ! Mixing length
  REAL :: pres2m             ! Pressure at 2 meters above the sfc.
  REAL :: lapse              ! Average lapse rate -- sfc. to 50 meters
  REAL :: frac               ! used in calculating press at 2 meters
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  thetas=(temps+273.16)*(1000./press)**.286
  frac=2./(htl-hts)
!  pres2m=press-EXP(frac*(LOG(press)-LOG(presl)))    ! origin
  pres2m=press/EXP(frac*(LOG(press)-LOG(presl)))     ! Ting Lei's fix
!
!----------------------------------------------------------------------
!
!  'x' is used in velocity profile and 'y' in temperature profile.
!  The 10 meter wind speed is calculated first, followed by the
!  2 meter temperature.
!
!----------------------------------------------------------------------
!
  IF (mixlength > 0.) THEN  ! (stable conditions)
    arps_spd=fricvel/k*(LOG(10/zo)+bm*(10/mixlength-zo/mixlength))
  ELSE ! (unstable conditions)
    x=(1-gm*10/mixlength)**(1./4.)
    xo=(1-gm*zo/mixlength)**(1./4.)
    arps_spd=fricvel/k*(LOG(10/zo)-(2*LOG((1+x)/(1+xo))+                &
        LOG((1+x**2)/(1+xo**2))-2*ATAN(x)+2*ATAN(xo)))
  END IF
!
  IF (mixlength > 0.) THEN
    dtheta=thetastar*pr/k*(LOG(2./zo)+bh*                               &
        (2./mixlength-zo/mixlength))
  ELSE
    y=(1-gh*2./mixlength)**(1./2.)
    yo=(1-gh*zo/mixlength)**(1./2.)
    dtheta=thetastar*pr/k*(LOG(2./zo)-2*                                &
        LOG((y+1)/(yo+1)))
  END IF
!
  arps_temp=(thetas+dtheta)*(pres2m/1000.)**.286 - 273.16
!
  RETURN
!
END SUBROUTINE arps2obs
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINE FINDLC2                   ######
!######                                                      ######
!######                Copyright (c) 1992-1994               ######
!######    Center for Analysis and Prediction of Storms      ######
!######    University of Oklahoma.  All rights reserved.     ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE findlc2(nx,ny,xs,ys,xpt,ypt,ipt,jpt,ireturn)
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
!    ireturn  status indicator, 0 = good
!                              -1 = extrapolation in x
!                              -2 = extrapolation in y
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
  INTEGER :: ireturn
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
  DO i=2,nx-2
    IF(xpt < xs(i)) EXIT
  END DO
  101 CONTINUE
  ipt=i-1
! Use "nx-2" not "nx-1", as MPI will overlap, and we don't want that.
! The cost is the edges in non-mpi mode, however, that is still in the
! fake zone, so it might not be that useful.
  IF(xpt > xs(nx-2) .OR. xpt < xs(1)) ireturn=-1
  DO j=2,ny-2
    IF(ypt < ys(j)) EXIT
  END DO
  201 CONTINUE
  jpt=j-1
  IF(ypt > ys(ny-2) .OR. ypt < ys(1)) ireturn=-2

  RETURN
END SUBROUTINE findlc2
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINE COLINTB                   ######
!######                                                      ######
!######                Copyright (c) 1992-1994               ######
!######    Center for Analysis and Prediction of Storms      ######
!######    University of Oklahoma.  All rights reserved.     ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE colintb(nx,ny,nz,nzmax,                                      &
           xs,ys,zp,xpt,ypt,ipt,jpt,                                    &
           uprt, vprt, ptprt, pprt, qvprt,                              &
           ubar, vbar, ptbar, pbar, qvbar,                              &
           raing,rainc,srain,                                           &
           su,sv,stheta,spres,shght,sqv,                                &
           tem1,nlevs)
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
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx,ny,nz Dimensions of ARPS grids.
!
!    nzmax    Maximum number of vertical levels allowed
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
!    uprt     x component of perturbation velocity (m/s)
!    vprt     y component of perturbation velocity (m/s)
!
!    ptprt    Perturbation potential temperature (K)
!    pprt     Perturbation pressure (Pascal)
!
!    qvprt    Perturbation water vapor mixing ratio (kg/kg)
!
!    ubar     Base state x velocity component (m/s)
!    vbar     Base state y velocity component (m/s)
!    ptbar    Base state potential temperature (K)
!    pbar     Base state pressure (Pascal)
!    qvbar    Base state water vapor mixing ratio (kg/kg)
!
!    raing    Supersaturation rainfall
!    rainc    Cumulus convective rainfall
!
!  OUTPUT:
!
!    su       Interpolated u wind component.  (m/s)
!    sv       Interpolated v wind component.  (m/s)
!    stheta   Interpolated potential temperature (K).
!    spres    Interpolated pressure. (Pascals)
!    shght    Interpolated height (meters)
!    sqv      Interpolated water vapor mixing ratio (kg/kg).
!    srain    Interpolated accumulated rainfall (m)
!    nlevs    Number of above-ground sounding levels.
!
!  WORK ARRAYS:
!
!    tem1     Temporary work array.
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

  INTEGER :: nx,ny,nz          ! Dimensions of ARPS grids.
  INTEGER :: nzmax             ! Maximum # of vertical levels allowed
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
!  Arguments -- model data
!
!-----------------------------------------------------------------------
!
  REAL :: uprt   (nx,ny,nz)    ! Perturbation u-velocity (m/s)
  REAL :: vprt   (nx,ny,nz)    ! Perturbation v-velocity (m/s)
  REAL :: ptprt  (nx,ny,nz)    ! Perturbation potential temperature (K)
  REAL :: pprt   (nx,ny,nz)    ! Perturbation pressure (Pascal)
  REAL :: qvprt  (nx,ny,nz)    ! Perturbation water vapor specific
                               ! humidity (kg/kg)

  REAL :: ubar   (nx,ny,nz)    ! Base state u-velocity (m/s)
  REAL :: vbar   (nx,ny,nz)    ! Base state v-velocity (m/s)
  REAL :: ptbar  (nx,ny,nz)    ! Base state potential temperature (K)
  REAL :: pbar   (nx,ny,nz)    ! Base state pressure (Pascal)
  REAL :: qvbar  (nx,ny,nz)    ! Base state water vapor specific
                               ! humidity (kg/kg)
  REAL :: raing  (nx,ny)       ! Grid supersaturation rainfall
  REAL :: rainc  (nx,ny)       ! Cumulus convective rainfall
!
!-----------------------------------------------------------------------
!
!  Arguments -- Extracted sounding variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: nlevs
!
  REAL :: su(nzmax)
  REAL :: sv(nzmax)
  REAL :: stheta(nzmax)
  REAL :: sqv(nzmax)
  REAL :: spres(nzmax)
  REAL :: shght(nzmax)
  REAL :: srain
!
!-----------------------------------------------------------------------
!
!  Temporary work arrays
!
!-----------------------------------------------------------------------
!
  REAL :: tem1(nx,ny,nz)
!
!-----------------------------------------------------------------------
!
!  Functions called
!
!-----------------------------------------------------------------------
!
  REAL :: aint2dtmp
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
  INTEGER :: i,j,k,in,jn
  REAL :: delx,ddx,dely,ddy,w1,w2,w3,w4
!
  REAL :: mindist,d1,d2,d3,d4
  INTEGER :: numprec
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
  DO k=1,nz
    shght(k)=                                                           &
        aint2dtmp(nx,ny,nz,    zp,ipt,jpt,k,in,jn,w1,w2,w3,w4)
    stheta(k)=                                                          &
        aint2dtmp(nx,ny,nz, ptprt,ipt,jpt,k,in,jn,w1,w2,w3,w4)          &
        +aint2dtmp(nx,ny,nz, ptbar,ipt,jpt,k,in,jn,w1,w2,w3,w4)
    sqv(k)=                                                             &
        aint2dtmp(nx,ny,nz, qvprt,ipt,jpt,k,in,jn,w1,w2,w3,w4)          &
        +aint2dtmp(nx,ny,nz, qvbar,ipt,jpt,k,in,jn,w1,w2,w3,w4)
    spres(k)=                                                           &
        aint2dtmp(nx,ny,nz,  pprt,ipt,jpt,k,in,jn,w1,w2,w3,w4)          &
        +aint2dtmp(nx,ny,nz,  pbar,ipt,jpt,k,in,jn,w1,w2,w3,w4)
  END DO
!
!-----------------------------------------------------------------------
!
!  Interpolate accumulated rainfall
!
!-----------------------------------------------------------------------
!
!
  DO i=1,nx
    DO j=1,ny
      tem1(i,j,1)=raing(i,j)+rainc(i,j)
    END DO
  END DO
!
  numprec=0 ! # of g.p.'s surrounding stn with precip forecast
  IF (tem1(ipt,jpt,1) > 0) THEN
    numprec=numprec+1
  END IF
  IF (tem1(ipt+1,jpt,1) > 0) THEN
    numprec=numprec+1
  END IF
  IF (tem1(ipt,jpt+1,1) > 0) THEN
    numprec=numprec+1
  END IF
  IF (tem1(ipt+1,jpt+1,1) > 0) THEN
    numprec=numprec+1
  END IF
!
  IF (numprec == 0) THEN ! no g.p.'s have precip forecast
!
    srain=0.
    GO TO 210
!
  ELSE ! at least 1 g.p. has precip forecast
!
    srain=aint2dtmp(nx,ny,nz,tem1,ipt,jpt,1,in,jn,w1,w2,w3,w4)
!
    IF (numprec >= 3) THEN ! if at least 3 have precip,
                               ! interpolation is finished
      GO TO 210
!
    END IF
!
  END IF
!
!----------------------------------------------------------------------
!
!  If less than 3 surrounding grid points have precip
!  forecast, do further checking to see if the nearest
!  grid point has precip forecast, and use that to determine
!  whether the station's precip should be set to zero or
!  not.
!
!----------------------------------------------------------------------
!
! SLOW!  Rewrite!!!
! d1=((xpt-xs(ipt))**2+(ypt-ys(jpt))**2)**(1./2.)
! d2=((xpt-xs(ipt+1))**2+(ypt-ys(jpt))**2)**(1./2.)
! d3=((xpt-xs(ipt))**2+(ypt-ys(jpt+1))**2)**(1./2.)
! d4=((xpt-xs(ipt+1))**2+(ypt-ys(jpt+1))**2)**(1./2.)

  d1=sqrt( (xpt-xs(ipt))   * (xpt-xs(ipt))   + (ypt-ys(jpt)) * (ypt-ys(jpt)))
  d2=sqrt( (xpt-xs(ipt+1)) * (xpt-xs(ipt+1)) + (ypt-ys(jpt)) * (ypt-ys(jpt)))
  d3=sqrt( (xpt-xs(ipt))   * (xpt-xs(ipt))   + (ypt-ys(jpt+1)) * (ypt-ys(jpt+1)))
  d4=sqrt( (xpt-xs(ipt+1)) * (xpt-xs(ipt+1)) + (ypt-ys(jpt+1)) * (ypt-ys(jpt+1)))

  mindist=AMIN1(d1,d2)
  mindist=AMIN1(mindist,d3)
  mindist=AMIN1(mindist,d4)
!
  IF (mindist == d1) THEN
    IF (tem1(ipt,jpt,1) == 0) THEN
      srain=0.
    END IF
  ELSE IF (mindist == d2) THEN
    IF (tem1(ipt+1,jpt,1) == 0) THEN
      srain=0.
    END IF
  ELSE IF (mindist == d3) THEN
    IF (tem1(ipt,jpt+1,1) == 0) THEN
      srain=0.
    END IF
  ELSE IF (mindist == d4) THEN
    IF (tem1(ipt+1,jpt+1,1) == 0) THEN
      srain=0.
    END IF
  ELSE
    WRITE(6,*) "mindist-d's:",mindist,d1,d2,d3,d4
  END IF
!
!
!-----------------------------------------------------------------------
!
!  Get height at scalar points, since zp was defined at w points.
!
!-----------------------------------------------------------------------
!
  210  CONTINUE
  DO k=1,nz-1
    shght(k)=0.5*(shght(k+1)+shght(k))
  END DO
!
!-----------------------------------------------------------------------
!
!  Form total u wind component at scalar points
!
!-----------------------------------------------------------------------
!
  DO k=1,nz-1
    DO j=1,ny-1
      DO i=1,nx-1
        tem1(i,j,k)=0.5*(ubar(i,j,k)+ubar(i+1,j,k))+                    &
                    0.5*(uprt(i,j,k)+uprt(i+1,j,k))
      END DO
    END DO
  END DO

  DO k=1,nz-1
    su(k)=                                                              &
        aint2dtmp(nx,ny,nz,  tem1,ipt,jpt,k,in,jn,w1,w2,w3,w4)
  END DO
!
!-----------------------------------------------------------------------
!
!  Form total v wind component at scalar points
!
!-----------------------------------------------------------------------
!
  DO k=1,nz-1
    DO j=1,ny-1
      DO i=1,nx-1
        tem1(i,j,k)=0.5*(vbar(i,j,k)+vbar(i,j+1,k)) +                   &
                    0.5*(vprt(i,j,k)+vprt(i,j+1,k))
      END DO
    END DO
  END DO

  DO k=1,nz-1
    sv(k)=                                                              &
        aint2dtmp(nx,ny,nz,  tem1,ipt,jpt,k,in,jn,w1,w2,w3,w4)
  END DO

!
  RETURN
END SUBROUTINE colintb
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINE CNVSND                    ######
!######                                                      ######
!######                Copyright (c) 1992-1994               ######
!######    Center for Analysis and Prediction of Storms      ######
!######    University of Oklahoma.  All rights reserved.     ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE cnvsnd(su,sv,stheta,spres,sqv,slon,                          &
           sdrct,ssped,stemp,sdewp,nlevs)

!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Converts units of data extracted from ARPS history data to
!  those required of sounding data. Determines direction and
!  speed from u and v wind components.
!
!  Dew-point formula from Bolton, 1980, MWR pp 1046-1053.
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
!  10/28/1992 (K. Brewster)
!  Special allowance for low qv-to-dew pt
!
!  09/12/1995 (K. Brewster)
!  Changed the determination of direc and speed to account
!  for map rotation factor.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    su       Sounding u wind component.  (m/s)
!    sv       Sounding v wind component.  (m/s)
!    stheta   Sounding potential temperature (K).
!    spres    Sounding pressure. (Pascals)
!    sqv      Sounding water vapor mixing ratio (kg/kg).
!    slon     Station longitude (degrees E)
!    nlevs    Number of above-ground sounding levels.
!
!  OUTPUT:
!
!    spres    Sounding pressure. (mb)
!    sdrct    Sounding wind direction (degrees from north)
!    ssped    Sounding wind speed (m/s)
!    stemp    Sounding temperature (degrees C)
!    sdewp    Sounding dew point temperature (degrees C)
!
!-----------------------------------------------------------------------
!
!  Variable declarations
!
!-----------------------------------------------------------------------
!
!  Input arguments
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nlevs             ! Number of above-ground sounding levels
  REAL :: su    (nlevs)        ! Sounding u wind component (m/s)
  REAL :: sv    (nlevs)        ! Sounding v wind component (m/s)
  REAL :: stheta(nlevs)        ! Sounding potential temperature (K)
  REAL :: spres (nlevs)        ! Sounding pressure. (Pascals)
  REAL :: sqv   (nlevs)        ! Sounding water vapor mixing
                               ! ratio (kg/kg)
  REAL :: slon
!
!-----------------------------------------------------------------------
!
!  Output arguments
!
!-----------------------------------------------------------------------
!
  REAL :: sdrct(nlevs)         ! Sounding wind direction
                               ! (degrees from north)
  REAL :: ssped(nlevs)         ! Sounding wind speed (m/s)
  REAL :: stemp(nlevs)         ! Sounding temperature (degrees C)
  REAL :: sdewp(nlevs)         ! Sounding dew point temperature
                               ! (degrees C)
!
!-----------------------------------------------------------------------
!
!  Constants
!
!-----------------------------------------------------------------------
!
  REAL :: rd2dg
  PARAMETER (rd2dg=(180./3.141592654))
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
  INTEGER :: k
  REAL :: smix,e,bige,alge
!
!-----------------------------------------------------------------------
!
!  Convert u,v to direction and speed
!
!-----------------------------------------------------------------------
!
  DO k=1,nlevs
    CALL uvrotdd(1,1,slon,su(k),sv(k),sdrct(k),ssped(k))
!
!-----------------------------------------------------------------------
!
!  Convert pressure from Pascals to mb
!
!-----------------------------------------------------------------------
!
    spres(k)=spres(k)*0.01
!
!-----------------------------------------------------------------------
!
!  Convert theta to temperature in degrees C.
!
!-----------------------------------------------------------------------
!
    stemp(k)=stheta(k)*((spres(k)/1000.)**rddcp)
    stemp(k)=stemp(k)-273.15
!
!-----------------------------------------------------------------------
!
!  Convert QV to dew-point in degrees C.
!
!-----------------------------------------------------------------------
!
    IF( sqv(k) > 1E-05) THEN
      smix=sqv(k)
      e=(spres(k)*smix)/(0.62197 + smix)
      bige=e/( 1.001 + ( (spres(k) - 100.) / 900) * 0.0034)
      alge = ALOG(bige/6.112)
      sdewp(k)= (alge * 243.5) / (17.67 - alge)
    ELSE
      sdewp(k)= stemp(k)-30.
    END IF

  END DO
!
  RETURN
END SUBROUTINE cnvsnd


!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE SLENGTH                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE slength ( string, length )

!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Return the length of the non-blank part of a string.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!
!  MODIFICATION HISTORY:
!
!  6/09/92 (K. Brewster)
!  Added full documentation and streamlined logic
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    string   character string to be sized
!
!  INPUT/OUTPUT:
!
!    length   on input, full size of character string
!             on output, true length of string as measured by the
!                        location of the last non-blank character
!
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
  INTEGER :: length
  CHARACTER (LEN=*) :: string
!
!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: i
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  DO i = length,1,-1
    IF(string(i:i) /= ' ') EXIT
  END DO
  101   CONTINUE

  length = i

  RETURN
END SUBROUTINE slength

!
!#################################################################
!#################################################################
!######                                                     ######
!######                SUBROUTINE MOSDUMP                   ######
!######                                                     ######
!######                Copyright (c) 1996                   ######
!######    Center for Analysis and Prediction of Storms     ######
!######    University of Oklahoma.  All rights reserved.    ######
!######                                                     ######
!#################################################################
!#################################################################
!

SUBROUTINE mosdump(nx,ny,nz,nzsoil,nstyps,x,y,z,zp,hterain,             &
           uprt,ubar,vprt,vbar,wprt,wbar,                               &
           ptprt,ptbar,pprt,pbar,qvprt,qvbar,                           &
           qscalar,tke,kmh,kmv,rhobar,                                  &
           zpsoil,tsoil,qsoil,wetcanp,                                  &
           snowdpth,soiltyp,stypfrct,vegtyp,                            &
           lai,roufns,veg,raing,rainc,prcrate,                          &
           radfrc,radsw,rnflx,radswnet,radlwin,                         &
           usflx,vsflx,ptsflx,                                          &
           qvsflx,dumpfile,needmosstns,                                 &
           xs,ys,xmap,ymap,latgr,longr,tem1,tem4,                       &
           moslist,mosdomlist,time)
!
!----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Coordinate the verification MOS dumps.  Based on subroutine
!  snddump
!
!----------------------------------------------------------------------
!
!  AUTHOR:  Eric Kemp
!
!  Modifications:
!
!  1 June 2002 Eric Kemp
!  Soil variable updates.
!
!-----------------------------------------------------------------------
!
!  Variables to be read in:
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
  INCLUDE 'bndry.inc'
  INCLUDE 'exbc.inc'
  INCLUDE 'vericst.inc'
  INCLUDE 'grid.inc'

!-----------------------------------------------------------------------

  INTEGER :: nx,ny,nz          ! Number of grid points in 3 directions
  INTEGER :: nzsoil            ! Number of soil levels.
  INTEGER :: nstyps            ! Number of soil type

  REAL :: x     (nx)           ! The x-coord. of the physical and
                               ! computational grid. Defined at u-point.
  REAL :: y     (ny)           ! The y-coord. of the physical and
                               ! computational grid. Defined at v-point.
  REAL :: z     (nz)           ! The z-coord. of the computational grid.
                               ! Defined at w-point on the staggered grid.
  REAL :: zp    (nx,ny,nz)     ! The physical height coordinate defined at
                               ! w-point of the staggered grid.
  REAL :: zpsoil(nx,ny,nzsoil) ! Height of soil levels
  REAL :: hterain(nx,ny)       ! Terrain height.

  REAL :: uprt  (nx,ny,nz)     ! Perturbation u velocity (m/s)
  REAL :: ubar  (nx,ny,nz)     ! Base state u-velocity (m/s)

  REAL :: vprt  (nx,ny,nz)     ! Perturbation v velocity (m/s)
  REAL :: vbar  (nx,ny,nz)     ! Base state v-velocity (m/s)

  REAL :: wprt  (nx,ny,nz)     ! Perturbation w velocity (m/s)
  REAL :: wbar  (nx,ny,nz)     ! Base state z velocity component (m/s)

  REAL :: ptprt (nx,ny,nz)     ! Perturbation potential temperature (K)
  REAL :: ptbar (nx,ny,nz)     ! Base state potential temperature (K)

  REAL :: pprt  (nx,ny,nz)     ! Perturbation pressure (Pascal)
  REAL :: pbar  (nx,ny,nz)     ! Base state pressure (Pascal)

  REAL :: qvprt (nx,ny,nz)     ! Perturbation wv mixing ration (kg/kg)
  REAL :: qvbar (nx,ny,nz)     ! Base state water vapor specific humidity
                               ! (kg/kg)
  REAL :: qscalar(nx,ny,nz,nscalar)

  REAL :: tke  (nx,ny,nz)      ! Turbulent Kinetic Energy ((m/s)**2)
  REAL :: kmh   (nx,ny,nz)     ! Horizontal turb. mixing coef. for
                               ! momentum. ( m**2/s )
  REAL :: kmv   (nx,ny,nz)     ! Vertical turb. mixing coef. for
                               ! momentum. ( m**2/s )
  REAL :: rhobar(nx,ny,nz)     ! Base state air density (kg/m**3)

  REAL :: tsoil  (nx,ny,nzsoil,0:nstyps) ! Soil temperature (K)
  REAL :: qsoil  (nx,ny,nzsoil,0:nstyps) ! Soil moisture (m**3/m**3))
  REAL :: wetcanp(nx,ny,0:nstyps)      ! Canopy water amount
  REAL :: snowdpth(nx,ny)              ! Snow depth (m)

  INTEGER :: soiltyp (nx,ny,nstyps)    ! Soil type
  REAL :: stypfrct(nx,ny,nstyps)    ! Soil type
  INTEGER :: vegtyp (nx,ny)            ! Vegetation type
  REAL :: lai    (nx,ny)            ! Leaf Area Index
  REAL :: roufns (nx,ny)            ! Surface roughness
  REAL :: veg    (nx,ny)            ! Vegetation fraction

  REAL :: raing(nx,ny)         ! Grid supersaturation rain
  REAL :: rainc(nx,ny)         ! Cumulus convective rain
  REAL :: prcrate(nx,ny,4)     ! precipitation rate (kg/(m**2*s))
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
  REAL :: ptsflx(nx,ny)        ! Surface heat flux (K*kg/(m*s**2))
  REAL :: qvsflx(nx,ny)        ! Surface moisture flux (kg/(m**2*s))


  CHARACTER (LEN=256) :: dumpfile     ! Sounding data history dump file name
  LOGICAL :: needmosstns       ! Create a stn list of stns in grid
  CHARACTER (LEN=256) :: moslist   ! File to read the sounding locations from
  CHARACTER (LEN=256) :: mosdomlist
!
  REAL :: time
!
!-----------------------------------------------------------------------
!
!  Computed variables
!
!-----------------------------------------------------------------------
!
  REAL :: xs(nx)      ! x location of scalar points
  REAL :: ys(ny)      ! y location of scalar points
!
!-----------------------------------------------------------------------
!
!  Extracted sounding variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: nlevs
  INTEGER :: nzmax
  PARAMETER (nzmax=250)
!
  REAL :: su(nzmax)
  REAL :: sv(nzmax)
  REAL :: sw(nzmax)
  REAL :: spt(nzmax)
  REAL :: sp(nzmax)
  REAL :: sqv(nzmax)
  REAL :: sqc(nzmax)
  REAL :: sqr(nzmax)
  REAL :: sqi(nzmax)
  REAL :: sqs(nzmax)
  REAL :: sqh(nzmax)
  REAL :: stke(nzmax)
  REAL :: skmh(nzmax)
  REAL :: skmv(nzmax)
  REAL :: srhobar(nzmax)
  REAL :: stsfc(0:nstyps)
  REAL :: stsoil(0:nstyps)
  REAL :: swetsfc(0:nstyps)
  REAL :: swetdp(0:nstyps)
  REAL :: swetcanp(0:nstyps)
  REAL :: ssnowdpth
  INTEGER :: ssoiltyp(nstyps)
  REAL :: sstypfrct(nstyps)
  INTEGER :: svegtyp
  REAL :: slai
  REAL :: sroufns
  REAL :: sveg
  REAL :: sraing
  REAL :: srainc
  REAL :: sprcrate(4)
  REAL :: sradfrc(nzmax)
  REAL :: sradsw
  REAL :: srnflx
  REAL :: sradswnet
  REAL :: sradlwin
  REAL :: susflx
  REAL :: sptflx
  REAL :: svsflx
  REAL :: sptsflx
  REAL :: sqvsflx
  REAL :: shght(nzmax)
  REAL :: shterain
  REAL :: srain
!
!-----------------------------------------------------------------------
!
!  Work Arrays
!
!-----------------------------------------------------------------------
!
  REAL :: tem1(nx,ny,nz)
  REAL :: tem4(nx,ny)
!
!-----------------------------------------------------------------------
!
!  Misc. internal variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: istnm      ! Station number of the sounding location
  INTEGER :: iselev_act ! Actual station elevation in integer format
  INTEGER :: i,j,k      ! Index variables
  INTEGER :: ii,ij      ! More index variables
  INTEGER :: ireturn    ! Flag for stations outside of the domain
!
  CHARACTER (LEN=256) :: line         ! Temp variable to read lines from files
  INTEGER :: LEN
!
!-----------------------------------------------------------------------
!
!  Map projection variables
!
!-----------------------------------------------------------------------
!
  REAL :: xctr,yctr
  REAL :: xll,yll
  REAL :: latnot(2)
  REAL :: xmap(nx)
  REAL :: ymap(ny)
  REAL :: latgr(nx,ny)
  REAL :: longr(nx,ny)
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  IF (nz > nzmax) THEN
    PRINT *,'Reset nzmax to greater than:',nz
    STOP
  END IF
!
!-----------------------------------------------------------------------
!
!  Calculate scalar locations and set up the map projection and
!  grid parameters
!
!-----------------------------------------------------------------------
!
  DO i=1,nx-1
    xs(i)=0.5*(x(i)+x(i+1))
  END DO
  DO j=1,ny-1
    ys(j)=0.5*(y(j)+y(j+1))
  END DO
  dx=x(2)-x(1)
  dy=y(2)-y(1)
!
  latnot(1)=trulat1
  latnot(2)=trulat2
  CALL setmapr(mapproj,sclfct,latnot,trulon)
  CALL lltoxy(1,1,ctrlat,ctrlon,xctr,yctr)
  xll=xctr-(0.5*(nx-3)*dx)
  yll=yctr-(0.5*(ny-3)*dy)

  DO i=1,nx-1
    xmap(i)=xll+xs(i)
  END DO
  xmap(nx)=2.*xmap(nx-1)-xmap(nx-2)
  DO j=1,ny-1
    ymap(j)=yll+ys(j)
  END DO
  ymap(ny)=2.*ymap(ny-1)-ymap(ny-2)
  CALL xytoll(nx,ny,xmap,ymap,latgr,longr)
!
!-----------------------------------------------------------------------
!
!  Find location in ARPS grid of all the stations in sndlist, then
!  rewrite only the ones that are in the grid to common arrays
!  so as to not make the program read all of the stations from
!  the original sounding stn list and check their locations with each
!  history dump.
!
!-----------------------------------------------------------------------
!

  OPEN(UNIT=1,FILE=moslist,STATUS='old')
!  OPEN(UNIT=2,FILE=dumpfile,STATUS='unknown')
!
!  LEN =256
!  CALL slength( dumpfile, LEN)
!  WRITE(6,*) "Opening ",dumpfile(1:LEN)," for MOS output"

  105  FORMAT(a8,a16,a20)
!
  IF (needmosstns) THEN
!
    needmosstns = .false.
    mosstn=0
!
    LEN =80
    CALL slength(moslist,LEN)
    WRITE(6,*) 'Reading sounding stations from file: ',                 &
                moslist(1:LEN)
!
    OPEN(UNIT=3,FILE=mosdomlist,STATUS='unknown')

    WRITE(3,990) year,'-',month,'-',day,'.',hour,':',minute,':',        &
               second,'+',INT(time)

    LEN =80
    CALL slength(runname,LEN)
    WRITE(3,'(a)') runname(1:LEN)
    WRITE(3,'(i1.1)') nocmnt
    DO i = 1,nocmnt
      LEN =80
      CALL slength(cmnt(i),LEN)
      WRITE(3,'(a)') cmnt(i)(1:LEN)
    END DO
!
    110    FORMAT(a80)
!
    READ(1,110,END=30) line ! read header lines and discard them..
    READ(1,110,END=30) line
    READ(1,110,END=30) line
!
    DO i=1,mosmax
!
      115      FORMAT(a3,3X,f5.2,2X,f7.2,1X,i4)
!
      READ(1,110,END=30) line
!
      mosstn=mosstn+1
!
      READ(line,115) mosstid(mosstn),moslat(mosstn),moslon(mosstn),     &
                     iselev_act
      istnm=0 !  Assumes a station number is unavailable
      moselev_act(mosstn)=FLOAT(iselev_act)
!
      CALL lltoxy(1,1,moslat(mosstn),moslon(mosstn),                    &
                       mosxpt(mosstn),mosypt(mosstn))
!
      mosxpt(mosstn)=mosxpt(mosstn)-xll
      mosypt(mosstn)=mosypt(mosstn)-yll
!
      CALL findlc2(nx,ny,xs,ys,mosxpt(mosstn),mosypt(mosstn),           &
                 mosipt(mosstn),mosjpt(mosstn),ireturn)
!
      IF (ireturn < 0) THEN  ! stn is outside the grid
        mosstn=mosstn-1
      ELSE  ! stn is inside the grid
        WRITE(3,110) line ! write location data to a file
      END IF
!
    END DO
!
  END IF
!
!-----------------------------------------------------------------------
!
!  Interpolate (in the horizontal) for the whole vertical column
!  for each station, derive the sounding variables, and write the
!  extracted sounding to dumpfile...
!
!-----------------------------------------------------------------------
!
  30   CONTINUE
!
  DO i=1,mosstn ! do for each station 'i' ...
!
    dumpfile = './mos_dumps/'
    LEN =256
    CALL slength( dumpfile, LEN)
    WRITE(dumpfile(LEN+1:LEN+7),'(a7)') mosstid(i)(1:3)//'_mos'
    OPEN(UNIT=2,FILE=dumpfile,STATUS='unknown',POSITION='append')

    WRITE(2,991) year,'-',month,'-',day,'.',hour,':',minute,':',        &
              second,'+',INT(time)
    991   FORMAT(i4.4,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i6.6)

    DO ii=1,nx
      DO ij=1,ny
        tem4(ii,ij)=0.
      END DO
    END DO
!
!-----------------------------------------------------------------------
!
!  Interpolate (in the horizontal) for the whole vertical column.
!
!-----------------------------------------------------------------------
!

!TODO:  Interpolate all tsoil and qsoil to station.
    CALL colintmos(nx,ny,nz,nzsoil,nstyps,nzmax,                        &
             xs,ys,zp,mosxpt(i),mosypt(i),mosipt(i),mosjpt(i),          &
             uprt,ubar,vprt,vbar,wprt,wbar,                             &
             ptprt,ptbar,pprt,pbar,qvprt,qvbar,                         &
             qscalar,tke,kmh,kmv,rhobar,                                &
             zpsoil,tsoil,qsoil,wetcanp,                                &
             snowdpth,soiltyp,stypfrct,vegtyp,                          &
             lai,roufns,veg,raing,rainc,prcrate,                        &
             radfrc,radsw,rnflx,radswnet,radlwin,                       &
             usflx,vsflx,ptsflx,                                        &
             qvsflx,su,sv,sw,spt,sp,sqv,sqc,sqr,                        &
             sqi,sqs,sqh,stke,skmh,skmv,srhobar,                        &
             stsfc,stsoil,swetsfc,swetdp,                               &
             swetcanp,ssnowdpth,ssoiltyp,sstypfrct,                     &
             svegtyp,slai,sroufns,sveg,sraing,srainc,                   &
             sprcrate,sradfrc,sradsw,srnflx,sradswnet,sradlwin,         &
             susflx,svsflx,sptsflx,sqvsflx,shght,shterain,              &
             tem1,nlevs)

!
!-----------------------------------------------------------------------
!
!    Convert u and v from map coordinates to earth coordinates
!
!-----------------------------------------------------------------------
!
    DO k = 1,nz-1
      CALL uvmptoe(1,1,su(k),sv(k),moslon(i),su(k),sv(k))
    END DO
!
!-----------------------------------------------------------------------
!
!    Write out soundings for the current model time..
!
!-----------------------------------------------------------------------
!
    WRITE(line,116) mosstid(i),moslat(i),moslon(i),                     &
                    shterain, (nlevs-2), nstyps
    116    FORMAT(a4,2X,f5.2,2X,f7.2,2X,f6.1,2X,i2,2X,i2)
    WRITE(2,110) line ! write location data to a file

    WRITE(2,990) year,'-',month,'-',day,'.',hour,':',minute,':',        &
               second,'+',INT(time)
    990     FORMAT(i4.4,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i6.6)

    LEN =80
    CALL slength(runname,LEN)
    WRITE(2,'(a)') runname(1:LEN)
    WRITE(2,'(i1.1)') nocmnt
    DO j = 1,nocmnt
      LEN =80
      CALL slength(cmnt(j),LEN)
      WRITE(2,'(a)') cmnt(j)(1:LEN)
    END DO

    WRITE(2,125)
    125    FORMAT(7X,                                                   &
           'hght             u                v                ',       &
           'w                ')
    WRITE(2,126)
    126    FORMAT(7X,                                                   &
           'm                m/s              m/s              ',       &
           'm/s              ')
    DO k = 2,(nlevs-1)
      WRITE(2,230) shght(k),su(k),sv(k),sw(k)
      230      FORMAT(1X,4(e16.10,1X))
    END DO

    WRITE(2,127)
    127    FORMAT(7X,                                                   &
           'pt               p                qv               ',       &
           'qc               ')
    WRITE(2,128)
    128    FORMAT(7X,                                                   &
           'K                Pa               kg/kg            ',       &
           'kg/kg            ')
    DO k = 2,(nlevs-1)
      WRITE(2,230) spt(k),sp(k),sqv(k),sqc(k)
    END DO

    WRITE(2,129)
    129    FORMAT(7X,                                                   &
           'qr               qi               qs               ',       &
           'qh               ')
    WRITE(2,130)
    130    FORMAT(7X,                                                   &
           'kg/kg            kg/kg            kg/kg            ',       &
           'kg/kg            ')
    DO k = 2,(nlevs-1)
      WRITE(2,230) sqr(k),sqi(k),sqs(k),sqh(k)
    END DO

    WRITE(2,131)
    131    FORMAT(7X,                                                   &
           'tke              kmh              kmv              ',       &
           'rhobar           ')
    WRITE(2,132)
    132    FORMAT(7X,                                                   &
           'm^2/s^2          m^2/s            m^2/s            ',       &
           'kg/m^3           ')
    DO k = 2,(nlevs-1)
      WRITE(2,230) stke(k),skmh(k),skmv(k),srhobar(k)
    END DO

    WRITE(2,133)
    133    FORMAT(7X,                                                   &
           'radfrc           ')
    WRITE(2,134)
    134    FORMAT(7X,                                                   &
           'K/s              ')
    DO k = 2,(nlevs-1)
      WRITE(2,231) sradfrc(k)
      231      FORMAT(1X,1(e16.10,1X))
    END DO

    WRITE(2,135)
    135    FORMAT(7X,                                                   &
           'tsfc             tsoil            wetsfc          ',        &
           'wetdp            ')
    WRITE(2,136)
    136    FORMAT(7X,                                                   &
           'K                K                                ',        &
           '                 ')
    DO k = 0,nstyps
      WRITE(2,240) stsfc(k),stsoil(k),swetsfc(k),swetdp(k)
      240      FORMAT(1X,4(e16.10,1X))
    END DO

    WRITE(2,137)
    137    FORMAT(7X,                                                   &
           'wetcanp          ')
    WRITE(2,138)
    138    FORMAT(7X,                                                   &
           '                 ')
    DO k = 0,nstyps
      WRITE(2,241) swetcanp(k)
      241      FORMAT(1X,1(e16.10,1X))
    END DO


    WRITE(2,145)
    145    FORMAT(7X,                                                   &
           'ssoiltyp         sstypfrct        ')
    146    FORMAT(7X,                                                   &
           '                                  ')
    WRITE(2,146)
    DO k = 1,nstyps
      WRITE(2,150) ssoiltyp(k),sstypfrct(k)
      150      FORMAT(1X,i16,1X e16.10)
    END DO

    WRITE(2,155)
    155    FORMAT(7X,                                                   &
           'snowdpth         vegtyp           lai             ',        &
           'roufns           ')
    WRITE(2,156)
    156    FORMAT(7X,                                                   &
           'm                                                 ',        &
           '                 ')
    WRITE(2,260) ssnowdpth,svegtyp,slai,sroufns
    260    FORMAT(1X,e16.10,1X,i16,1X,2(e16.10,1X))


    WRITE(2,157)
    157    FORMAT(7X,                                                   &
           'veg              raing            rainc           ',        &
           'radsw            ')
    WRITE(2,158)
    158    FORMAT(7X,                                                   &
           '                                                  ',        &
           '                 ')
    WRITE(2,261) sveg,sraing,srainc,sradsw
    261    FORMAT(1X,4(e16.10,1X))

!    WRITE(2,159)
!    159    FORMAT(7X,                                                   &
!           'rnflx            usflx            vsflx           ',        &
!           'ptsflx           ')
!    WRITE(2,160)
!    160    FORMAT(7X,                                                   &
!           '                 kg/m/s^2         kg/m/s^2        ',        &
!           'K kg/m/s^2       ')
!    WRITE(2,261) srnflx,susflx,svsflx,sptsflx

    WRITE(2,159)
    159    FORMAT(7X,                                                   &
           'rnflx            radswnet         radlwin         ',        &
           'usflx            ')
    WRITE(2,160)
    160    FORMAT(7X,                                                   &
           '                                                  ',        &
           'K kg/m/s^2       ')
    WRITE(2,261) srnflx,sradswnet,radlwin,susflx

!    WRITE(2,161)
!    161    FORMAT(7X,                                                   &
!           'qvsflx           ')
!    WRITE(2,162)
!    162    FORMAT(7X,                                                   &
!           'kg/s/m^2         ')
!    WRITE(2,262) sqvsflx
!    262    FORMAT(1X,1(e16.10,1X))

    WRITE(2,161)
    161    FORMAT(7X,                                                   &
           'vsflx            ptsflx           qvsflx          ')
    WRITE(2,162)
    162    FORMAT(7X,                                                   &
           'kg/s/m^2         K kg/m/s^2       kg/m/s^2        ')
    WRITE(2,262) svsflx,sptflx,sqvsflx
    262    FORMAT(1X,3(e16.10,1X))


    WRITE(2,165)
    165    FORMAT(7X,                                                   &
           'prcrate         ')
    WRITE(2,166)
    166    FORMAT(7X,                                                   &
           'kg/m^2/s        ')
    DO k = 1,4
      WRITE(2,170) sprcrate(k)
      170      FORMAT(1X,e16.10,i16,11(e16.10))
    END DO

    CLOSE(2)

  END DO ! Move on to the next station..
!
  CLOSE(1)
!  CLOSE(2)
  CLOSE(3)
!
  RETURN
!
END SUBROUTINE mosdump

!##################################################################
!##################################################################
!######                                                      ######
!######               SUBROUTINE COLINTMOS                   ######
!######                                                      ######
!######    Center for Analysis and Prediction of Storms      ######
!######    University of Oklahoma.  All rights reserved.     ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE colintmos(nx,ny,nz,nzsoil,nstyps,nzmax,                      &
           xs,ys,zp,xpt,ypt,ipt,jpt,                                    &
           uprt,ubar,vprt,vbar,wprt,wbar,                               &
           ptprt,ptbar,pprt,pbar,qvprt,qvbar,                           &
           qscalar,tke,kmh,kmv,rhobar,                                  &
           zpsoil,tsoil,qsoil,wetcanp,                                  &
           snowdpth,soiltyp,stypfrct,vegtyp,                            &
           lai,roufns,veg,raing,rainc,prcrate,                          &
           radfrc,radsw,rnflx,radswnet,radlwin,                         &
           usflx,vsflx,ptsflx,                                          &
           qvsflx,su,sv,sw,spt,sp,sqv,sqc,sqr,                          &
           sqi,sqs,sqh,stke,skmh,skmv,srhobar,                          &
           stsfc,stsoil,swetsfc,swetdp,                                 &
           swetcanp,ssnowdpth,ssoiltyp,sstypfrct,                       &
           svegtyp,slai,sroufns,sveg,sraing,srainc,                     &
           sprcrate,sradfrc,sradsw,srnflx,sradswnet,sradlwin,           &
           susflx,svsflx,sptsflx,sqvsflx,shght,shterain,                &
           tem1,nlevs)

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
!  Based on subroutine COLINTB
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Eric Kemp
!  May 2000
!
!  Modified 1 June 2002 Eric Kemp
!  Updated soil variables.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx,ny,nz Dimensions of ARPS grids.
!
!    nzmax    Maximum number of vertical levels allowed
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
!    uprt     x component of perturbation velocity (m/s)
!    vprt     y component of perturbation velocity (m/s)
!
!    ptprt    Perturbation potential temperature (K)
!    pprt     Perturbation pressure (Pascal)
!
!    qvprt    Perturbation water vapor mixing ratio (kg/kg)
!
!    ubar     Base state x velocity component (m/s)
!    vbar     Base state y velocity component (m/s)
!    ptbar    Base state potential temperature (K)
!    pbar     Base state pressure (Pascal)
!    qvbar    Base state water vapor mixing ratio (kg/kg)
!
!    raing    Supersaturation rainfall
!    rainc    Cumulus convective rainfall
!
!  OUTPUT:
!
!    su       Interpolated u wind component.  (m/s)
!    sv       Interpolated v wind component.  (m/s)
!    stheta   Interpolated potential temperature (K).
!    spres    Interpolated pressure. (Pascals)
!    shght    Interpolated height (meters)
!    sqv      Interpolated water vapor mixing ratio (kg/kg).
!    srain    Interpolated accumulated rainfall (m)
!    nlevs    Number of scalar levels.
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
!  Include files
!
!-----------------------------------------------------------------------
!
  INCLUDE 'phycst.inc'
  INCLUDE 'globcst.inc'

!-----------------------------------------------------------------------
!
!  Variable Declarations:
!
!-----------------------------------------------------------------------
!
!  Arguments -- location data
!
!-----------------------------------------------------------------------

  INTEGER :: nx,ny,nz          ! Dimensions of ARPS grids.
  INTEGER :: nzsoil            ! Number of soil levels
  INTEGER :: nstyps
  INTEGER :: nzmax             ! Maximum # of vertical levels allowed
  REAL :: xs(nx)               ! x coordinate of grid points in
                               ! physical/comp. space (m)
  REAL :: ys(ny)               ! y coordinate of grid points in
                               ! physical/comp. space (m)
  REAL :: zp(nx,ny,nz)         ! z coordinate of grid points in
                               ! physical space (m)
  REAL :: zpsoil(nx,ny,nzsoil) ! Soil level height (m)
  REAL :: xpt                  ! location to find in x coordinate (m)
  REAL :: ypt                  ! location to find in y coordinate (m)
  INTEGER :: ipt               ! i index to the west of desired
                               ! location
  INTEGER :: jpt               ! j index to the south of desired
                               ! location
!
!-----------------------------------------------------------------------
!
!  Arguments -- model data
!
!-----------------------------------------------------------------------
!

  REAL :: uprt  (nx,ny,nz)     ! Perturbation u velocity (m/s)
  REAL :: ubar  (nx,ny,nz)     ! Base state u-velocity (m/s)

  REAL :: vprt  (nx,ny,nz)     ! Perturbation v velocity (m/s)
  REAL :: vbar  (nx,ny,nz)     ! Base state v-velocity (m/s)

  REAL :: wprt  (nx,ny,nz)     ! Perturbation w velocity (m/s)
  REAL :: wbar  (nx,ny,nz)     ! Base state z velocity component (m/s)

  REAL :: ptprt (nx,ny,nz)     ! Perturbation potential temperature (K)
  REAL :: ptbar (nx,ny,nz)     ! Base state potential temperature (K)

  REAL :: pprt  (nx,ny,nz)     ! Perturbation pressure (Pascal)
  REAL :: pbar  (nx,ny,nz)     ! Base state pressure (Pascal)

  REAL :: qvprt (nx,ny,nz)     ! Perturbation wv mixing ration (kg/kg)
  REAL :: qvbar (nx,ny,nz)     ! Base state water vapor specific humidity
                               ! (kg/kg)
  REAL :: qscalar(nx,ny,nz,nscalar)

  REAL :: tke  (nx,ny,nz)      ! Turbulent Kinetic Energy ((m/s)**2)
  REAL :: kmh   (nx,ny,nz)     ! Horizontal turb. mixing coef. for
                               ! momentum. ( m**2/s )
  REAL :: kmv   (nx,ny,nz)     ! Vertical turb. mixing coef. for
                               ! momentum. ( m**2/s )
  REAL :: rhobar(nx,ny,nz)     ! Base state air density (kg/m**3)

  REAL :: tsoil  (nx,ny,nzsoil,0:nstyps) ! Soil temperature (K)
  REAL :: qsoil  (nx,ny,nzsoil,0:nstyps) ! Soil moisture (m**3/m**3)
  REAL :: wetcanp(nx,ny,0:nstyps)      ! Canopy water amount
  REAL :: snowdpth(nx,ny)              ! Snow depth (m)

  INTEGER :: soiltyp (nx,ny,nstyps)    ! Soil type
  REAL :: stypfrct(nx,ny,nstyps)    ! Soil type
  INTEGER :: vegtyp (nx,ny)            ! Vegetation type
  REAL :: lai    (nx,ny)            ! Leaf Area Index
  REAL :: roufns (nx,ny)            ! Surface roughness
  REAL :: veg    (nx,ny)            ! Vegetation fraction

  REAL :: raing(nx,ny)         ! Grid supersaturation rain
  REAL :: rainc(nx,ny)         ! Cumulus convective rain
  REAL :: prcrate(nx,ny,4)     ! precipitation rate (kg/(m**2*s))
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
  REAL :: ptsflx(nx,ny)        ! Surface heat flux (K*kg/(m*s**2))
  REAL :: qvsflx(nx,ny)        ! Surface moisture flux (kg/(m**2*s))

!
!-----------------------------------------------------------------------
!
!  Arguments -- Extracted sounding variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: nlevs
!
  REAL :: su(nzmax)
  REAL :: sv(nzmax)
  REAL :: sw(nzmax)
  REAL :: spt(nzmax)
  REAL :: sp(nzmax)
  REAL :: sqv(nzmax)
  REAL :: sqc(nzmax)
  REAL :: sqr(nzmax)
  REAL :: sqi(nzmax)
  REAL :: sqs(nzmax)
  REAL :: sqh(nzmax)
  REAL :: stke(nzmax)
  REAL :: skmh(nzmax)
  REAL :: skmv(nzmax)
  REAL :: srhobar(nzmax)
  REAL :: stsfc(0:nstyps)
  REAL :: stsoil(0:nstyps)
  REAL :: swetsfc(0:nstyps)
  REAL :: swetdp(0:nstyps)
  REAL :: swetcanp(0:nstyps)
  REAL :: ssnowdpth
  INTEGER :: ssoiltyp(nstyps)
  REAL :: sstypfrct(nstyps)
  INTEGER :: svegtyp
  REAL :: slai
  REAL :: sroufns
  REAL :: sveg
  REAL :: sraing
  REAL :: srainc
  REAL :: sprcrate(4)
  REAL :: sradfrc(nzmax)
  REAL :: sradsw
  REAL :: srnflx
  REAL :: sradswnet
  REAL :: sradlwin
  REAL :: susflx
  REAL :: svsflx
  REAL :: sptsflx
  REAL :: sqvsflx
  REAL :: shght(nzmax)
  REAL :: shterain
  REAL :: srain
  REAL :: stmp

  REAL :: p(nx,ny,nz)
  REAL :: pt(nx,ny,nz)
  REAL :: qv(nx,ny,nz)
!
!-----------------------------------------------------------------------
!
!  Temporary work arrays
!
!-----------------------------------------------------------------------
!
  REAL :: tem1(nx,ny,nz)
!
!-----------------------------------------------------------------------
!
!  Functions called
!
!-----------------------------------------------------------------------
!
  REAL :: aint2dtmp
!
!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k,in,jn
  REAL :: delx,ddx,dely,ddy,w1,w2,w3,w4
!
  REAL :: mindist,d1,d2,d3,d4
  INTEGER :: numprec
!
!-----------------------------------------------------------------------
!
!  Find corner weights
!
!-----------------------------------------------------------------------
!
  DO k = 1,nz
    DO j = 1,ny
      DO i = 1,nx
        p(i,j,k) = pbar(i,j,k) + pprt(i,j,k)
        qv(i,j,k) = qvbar(i,j,k) + qvprt(i,j,k)
        pt(i,j,k) = ptbar(i,j,k) + ptprt(i,j,k)
      END DO
    END DO
  END DO

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

  w1=(1.-ddx*(1.-ddy))
  w3=ddx*ddy
  w4=(1.-ddx)*ddy
!
!-----------------------------------------------------------------------
!
!  Interpolate all variables at all levels.
!
!-----------------------------------------------------------------------
!
  sqc(:) = 0.0
  sqr(:) = 0.0
  sqi(:) = 0.0
  sqs(:) = 0.0
  sqh(:) = 0.0

  nlevs=nz-1
  DO k=1,nz
    shght(k)=                                                           &
        aint2dtmp(nx,ny,nz,    zp,ipt,jpt,k,in,jn,w1,w2,w3,w4)
    spt(k)=                                                             &
        aint2dtmp(nx,ny,nz, pt,ipt,jpt,k,in,jn,w1,w2,w3,w4)
    sqv(k)=                                                             &
        aint2dtmp(nx,ny,nz, qv,ipt,jpt,k,in,jn,w1,w2,w3,w4)
    sp(k)=                                                              &
        aint2dtmp(nx,ny,nz,  p,ipt,jpt,k,in,jn,w1,w2,w3,w4)

    IF (P_QC > 0)  sqc(k)=  aint2dtmp(nx,ny,nz,                         &
            qscalar(:,:,:,P_QC),ipt,jpt,k,in,jn,w1,w2,w3,w4)
    IF (P_QR > 0)  sqr(k)=  aint2dtmp(nx,ny,nz,                         &
            qscalar(:,:,:,P_QR),ipt,jpt,k,in,jn,w1,w2,w3,w4)
    IF (P_QI > 0)  sqi(k)=  aint2dtmp(nx,ny,nz,                         &
            qscalar(:,:,:,P_QI),ipt,jpt,k,in,jn,w1,w2,w3,w4)
    IF (P_QS > 0)  sqs(k)=  aint2dtmp(nx,ny,nz,                         &
            qscalar(:,:,:,P_QS),ipt,jpt,k,in,jn,w1,w2,w3,w4)
    IF (P_QH > 0)  sqh(k)=  aint2dtmp(nx,ny,nz,                         &
            qscalar(:,:,:,P_QH),ipt,jpt,k,in,jn,w1,w2,w3,w4)

    stke(k)=                                                            &
        aint2dtmp(nx,ny,nz,tke,ipt,jpt,k,in,jn,w1,w2,w3,w4)
    skmh(k)=                                                            &
        aint2dtmp(nx,ny,nz,kmh,ipt,jpt,k,in,jn,w1,w2,w3,w4)
    skmv(k)=                                                            &
        aint2dtmp(nx,ny,nz,kmv,ipt,jpt,k,in,jn,w1,w2,w3,w4)
    srhobar(k)=                                                         &
        aint2dtmp(nx,ny,nz,rhobar,ipt,jpt,k,in,jn,w1,w2,w3,w4)
    sradfrc(k)=                                                         &
        aint2dtmp(nx,ny,nz,radfrc,ipt,jpt,k,in,jn,w1,w2,w3,w4)
  END DO

!
!-----------------------------------------------------------------------
!
!  Bilinearly interpolate snowdpth and radsw, and use nearest
!  neighbor values for other soil variables.
!
!-----------------------------------------------------------------------
!

  CALL colintsoil2d(nx,ny,nz,tem1,snowdpth,ssnowdpth,                   &
                          ipt,jpt,in,jn,w1,w2,w3,w4)
  CALL colintsoil2d(nx,ny,nz,tem1,radsw,sradsw,                         &
                          ipt,jpt,in,jn,w1,w2,w3,w4)

  d1=((xpt-xs(ipt))**2+(ypt-ys(jpt))**2)**(1./2.)
  d2=((xpt-xs(ipt+1))**2+(ypt-ys(jpt))**2)**(1./2.)
  d3=((xpt-xs(ipt))**2+(ypt-ys(jpt+1))**2)**(1./2.)
  d4=((xpt-xs(ipt+1))**2+(ypt-ys(jpt+1))**2)**(1./2.)

  mindist=AMIN1(d1,d2)
  mindist=AMIN1(mindist,d3)
  mindist=AMIN1(mindist,d4)
!
  IF (mindist == d1) THEN
    DO k = 1,nstyps
      ssoiltyp(k) = soiltyp(ipt,jpt,k)
      sstypfrct(k) = stypfrct(ipt,jpt,k)
    END DO
    DO k = 0,nstyps
      stsfc(k) = tsoil(ipt,jpt,1,k)
      stsoil(k) = tsoil(ipt,jpt,nzsoil,k)
      swetsfc(k) = qsoil(ipt,jpt,1,k)
      swetdp(k) = qsoil(ipt,jpt,nzsoil,k)
      swetcanp(k) = wetcanp(ipt,jpt,k)
    END DO
    svegtyp = vegtyp(ipt,jpt)
    slai = lai(ipt,jpt)
    sroufns = roufns(ipt,jpt)
    sveg = veg(ipt,jpt)
    srnflx = rnflx(ipt,jpt)
    sradswnet = radswnet(ipt,jpt)
    sradlwin = radlwin(ipt,jpt)
    susflx = usflx(ipt,jpt)
    svsflx = vsflx(ipt,jpt)
    sptsflx = ptsflx(ipt,jpt)
    sqvsflx = qvsflx(ipt,jpt)
  ELSE IF (mindist == d2) THEN
    DO k = 1,nstyps
      ssoiltyp(k) = soiltyp(ipt+1,jpt,k)
      sstypfrct(k) = stypfrct(ipt+1,jpt,k)
    END DO
    DO k = 0,nstyps
      stsfc(k) = tsoil(ipt+1,jpt,1,k)
      stsoil(k) = tsoil(ipt+1,jpt,nzsoil,k)
      swetsfc(k) = qsoil(ipt+1,jpt,1,k)
      swetdp(k) = qsoil(ipt+1,jpt,nzsoil,k)
      swetcanp(k) = wetcanp(ipt+1,jpt,k)
    END DO
    svegtyp = vegtyp(ipt+1,jpt)
    slai = lai(ipt+1,jpt)
    sroufns = roufns(ipt+1,jpt)
    sveg = veg(ipt+1,jpt)
    srnflx = rnflx(ipt+1,jpt)
    sradswnet = radswnet(ipt+1,jpt)
    sradlwin = radlwin(ipt+1,jpt)
    susflx = usflx(ipt+1,jpt)
    svsflx = vsflx(ipt+1,jpt)
    sptsflx = ptsflx(ipt+1,jpt)
    sqvsflx = qvsflx(ipt+1,jpt)
  ELSE IF (mindist == d3) THEN
    DO k = 1,nstyps
      ssoiltyp(k) = soiltyp(ipt,jpt+1,k)
      sstypfrct(k) = stypfrct(ipt,jpt+1,k)
    END DO
    DO k = 0,nstyps
      stsfc(k) = tsoil(ipt,jpt+1,1,k)
      stsoil(k) = tsoil(ipt,jpt+1,nzsoil,k)
      swetsfc(k) = qsoil(ipt,jpt+1,1,k)
      swetdp(k) = qsoil(ipt,jpt+1,nzsoil,k)
      swetcanp(k) = wetcanp(ipt,jpt+1,k)
    END DO
    svegtyp = vegtyp(ipt,jpt+1)
    slai = lai(ipt,jpt+1)
    sroufns = roufns(ipt,jpt+1)
    sveg = veg(ipt,jpt+1)
    srnflx = rnflx(ipt,jpt+1)
    sradswnet = radswnet(ipt,jpt+1)
    sradlwin = radlwin(ipt,jpt+1)
    susflx = usflx(ipt,jpt+1)
    svsflx = vsflx(ipt,jpt+1)
    sptsflx = ptsflx(ipt,jpt+1)
    sqvsflx = qvsflx(ipt,jpt+1)
  ELSE IF (mindist == d4) THEN
    DO k = 1,nstyps
      ssoiltyp(k) = soiltyp(ipt+1,jpt+1,k)
      sstypfrct(k) = stypfrct(ipt+1,jpt+1,k)
    END DO
    DO k = 0,nstyps
      stsfc(k) = tsoil(ipt+1,jpt+1,1,k)
      stsoil(k) = tsoil(ipt+1,jpt+1,nzsoil,k)
      swetsfc(k) = qsoil(ipt+1,jpt+1,1,k)
      swetdp(k) = qsoil(ipt+1,jpt+1,nzsoil,k)
      swetcanp(k) = wetcanp(ipt+1,jpt+1,k)
    END DO
    svegtyp = vegtyp(ipt+1,jpt+1)
    slai = lai(ipt+1,jpt+1)
    sroufns = roufns(ipt+1,jpt+1)
    sveg = veg(ipt+1,jpt+1)
    srnflx = rnflx(ipt+1,jpt+1)
    sradswnet = radswnet(ipt+1,jpt+1)
    sradlwin = radlwin(ipt+1,jpt+1)
    susflx = usflx(ipt+1,jpt+1)
    svsflx = vsflx(ipt+1,jpt+1)
    sptsflx = ptsflx(ipt+1,jpt+1)
    sqvsflx = qvsflx(ipt+1,jpt+1)
  ELSE
    WRITE(6,*) "mindist-d's:",mindist,d1,d2,d3,d4
    DO k = 1,nstyps
      ssoiltyp(k) = soiltyp(ipt,jpt,k)
      sstypfrct(k) = stypfrct(ipt,jpt,k)
    END DO
    DO k = 0,nstyps
      stsfc(k) = tsoil(ipt,jpt,1,k)
      stsoil(k) = tsoil(ipt,jpt,nzsoil,k)
      swetsfc(k) = qsoil(ipt,jpt,1,k)
      swetdp(k) = qsoil(ipt,jpt,nzsoil,k)
      swetcanp(k) = wetcanp(ipt,jpt,k)
    END DO
    svegtyp = vegtyp(ipt,jpt)
    slai = lai(ipt,jpt)
    sroufns = roufns(ipt,jpt)
    sveg = veg(ipt,jpt)
    srnflx = rnflx(ipt,jpt)
    sradswnet = radswnet(ipt,jpt)
    sradlwin = radlwin(ipt,jpt)
    susflx = usflx(ipt,jpt)
    svsflx = vsflx(ipt,jpt)
    sptsflx = ptsflx(ipt,jpt)
    sqvsflx = qvsflx(ipt,jpt)
  END IF
!
!-----------------------------------------------------------------------
!
!  Interpolate accumulated precip. and rainrates
!
!-----------------------------------------------------------------------
!
  DO i=1,nx
    DO j=1,ny
      tem1(i,j,1)=raing(i,j)
    END DO
  END DO
  CALL colintprec(nx,ny,nz,tem1,sraing,ipt,jpt,in,jn,                   &
                        xpt,ypt,xs,ys,w1,w2,w3,w4)

  DO i=1,nx
    DO j=1,ny
      tem1(i,j,1)=rainc(i,j)
    END DO
  END DO
  CALL colintprec(nx,ny,nz,tem1,srainc,ipt,jpt,in,jn,                   &
                        xpt,ypt,xs,ys,w1,w2,w3,w4)

  srain = srainc+sraing

  DO i=1,nx
    DO j=1,ny
      tem1(i,j,1)=prcrate(i,j,1)
    END DO
  END DO
  CALL colintprec(nx,ny,nz,tem1,stmp,ipt,jpt,in,jn,                     &
                        xpt,ypt,xs,ys,w1,w2,w3,w4)
  sprcrate(1) =stmp

  DO i=1,nx
    DO j=1,ny
      tem1(i,j,1)=prcrate(i,j,2)
    END DO
  END DO
  CALL colintprec(nx,ny,nz,tem1,stmp,ipt,jpt,in,jn,                     &
                        xpt,ypt,xs,ys,w1,w2,w3,w4)
  sprcrate(2) =stmp

  DO i=1,nx
    DO j=1,ny
      tem1(i,j,1)=prcrate(i,j,3)
    END DO
  END DO
  CALL colintprec(nx,ny,nz,tem1,stmp,ipt,jpt,in,jn,                     &
                        xpt,ypt,xs,ys,w1,w2,w3,w4)
  sprcrate(3) =stmp

  DO i=1,nx
    DO j=1,ny
      tem1(i,j,1)=prcrate(i,j,4)
    END DO
  END DO
  CALL colintprec(nx,ny,nz,tem1,stmp,ipt,jpt,in,jn,                     &
                        xpt,ypt,xs,ys,w1,w2,w3,w4)
  sprcrate(4) =stmp
!
!
!-----------------------------------------------------------------------
!
!  Get height at scalar points, since zp was defined at w points.
!
!-----------------------------------------------------------------------
!
  shterain = shght(2)
  DO k=1,nz-1
    shght(k)=0.5*(shght(k+1)+shght(k))
  END DO
!
!-----------------------------------------------------------------------
!
!  Form total u wind component at scalar points
!
!-----------------------------------------------------------------------
!
  DO k=1,nz-1
    DO j=1,ny-1
      DO i=1,nx-1
        tem1(i,j,k)=0.5*(ubar(i,j,k)+ubar(i+1,j,k))+                    &
                    0.5*(uprt(i,j,k)+uprt(i+1,j,k))
      END DO
    END DO
  END DO

  DO k=1,nz-1
    su(k)=                                                              &
        aint2dtmp(nx,ny,nz,  tem1,ipt,jpt,k,in,jn,w1,w2,w3,w4)
  END DO
!
!-----------------------------------------------------------------------
!
!  Form total v wind component at scalar points
!
!-----------------------------------------------------------------------
!
  DO k=1,nz-1
    DO j=1,ny-1
      DO i=1,nx-1
        tem1(i,j,k)=0.5*(vbar(i,j,k)+vbar(i,j+1,k)) +                   &
                    0.5*(vprt(i,j,k)+vprt(i,j+1,k))
      END DO
    END DO
  END DO

  DO k=1,nz-1
    sv(k)=                                                              &
        aint2dtmp(nx,ny,nz,  tem1,ipt,jpt,k,in,jn,w1,w2,w3,w4)
  END DO
!

!-----------------------------------------------------------------------
!
!  Form total w wind component at scalar points
!
!-----------------------------------------------------------------------
!
  DO k=1,nz-1
    DO j=1,ny-1
      DO i=1,nx-1
        tem1(i,j,k)=0.5*(wbar(i,j,k)+wbar(i,j,k+1)) +                   &
                    0.5*(wprt(i,j,k)+wprt(i,j,k+1))
      END DO
    END DO
  END DO

  DO k=1,nz-1
    sw(k)=                                                              &
        aint2dtmp(nx,ny,nz,  tem1,ipt,jpt,k,in,jn,w1,w2,w3,w4)
  END DO
!
  RETURN
END SUBROUTINE colintmos


!##################################################################
!##################################################################
!######                                                      ######
!######             SUBROUTINE COLINTSOILSD                  ######
!######                                                      ######
!######    Center for Analysis and Prediction of Storms      ######
!######    University of Oklahoma.  All rights reserved.     ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE colintsoilsd(nx,ny,nstyps,tem1,soilvar,ssoilvar,             &
           ipt,jpt,in,jn,w1,w2,w3,w4)

!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Interpolates "special dimension" soil data array in the
!  horizontal to create a column of data.  Must be called by
!  subroutine COLINTMOS.
!
!  Bilinear interpolation is used.
!
!  Based on subroutine COLINTB
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Eric Kemp
!  May 2000
!
!-----------------------------------------------------------------------
!
!  Arguments
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER :: nx,ny,nstyps
  REAL :: tem1(nx,ny,nstyps)
  REAL :: soilvar(nx,ny,0:nstyps),ssoilvar(0:nstyps)

  INTEGER :: ipt               ! i index to the west of desired
                               ! location
  INTEGER :: jpt               ! j index to the south of desired
                               ! location
  INTEGER :: in,jn
  REAL :: w1,w2,w3,w4

!-----------------------------------------------------------------------
!
!  Miscellaneous variables.
!
!-----------------------------------------------------------------------

  INTEGER :: i,j,k

!-----------------------------------------------------------------------
!
!  Functions called
!
!-----------------------------------------------------------------------

  REAL :: aint2dtmp ! function

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  DO k=1,nstyps+1
    DO j=1,ny-1
      DO i=1,nx-1
        tem1(i,j,k)= soilvar(i,j,k-1)
      END DO
    END DO
  END DO

  DO k=1,nstyps+1
    ssoilvar(k-1)=                                                      &
        aint2dtmp(nx,ny,(nstyps+1),tem1,ipt,jpt,k,in,jn,w1,w2,          &
                  w3,w4)

  END DO

  RETURN
END SUBROUTINE colintsoilsd

!
!##################################################################
!##################################################################
!######                                                      ######
!######             SUBROUTINE COLINTSOIL2D                  ######
!######                                                      ######
!######    Center for Analysis and Prediction of Storms      ######
!######    University of Oklahoma.  All rights reserved.     ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE colintsoil2d(nx,ny,nz,tem1,soilvar,ssoilvar,                 &
           ipt,jpt,in,jn,w1,w2,w3,w4)

!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Interpolates two dimension soil data array in the
!  horizontal to create a column of data.  Must be called by
!  subroutine COLINTMOS.
!
!  Bilinear interpolation is used.
!
!  Based on subroutine COLINTB.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Eric Kemp
!  May 2000
!
!-----------------------------------------------------------------------
!
!  Arguments
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER :: nx,ny,nz
  REAL :: tem1(nx,ny,nz)
  REAL :: soilvar(nx,ny),ssoilvar

  INTEGER :: ipt               ! i index to the west of desired
                               ! location
  INTEGER :: jpt               ! j index to the south of desired
                               ! location
  INTEGER :: in,jn
  REAL :: w1,w2,w3,w4

!-----------------------------------------------------------------------
!
!  Miscellaneous variables.
!
!-----------------------------------------------------------------------

  INTEGER :: i,j,k

!-----------------------------------------------------------------------
!
!  Functions called
!
!-----------------------------------------------------------------------

  REAL :: aint2dtmp ! function

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  DO k=1,nz-1
    DO j=1,ny-1
      DO i=1,nx-1
        tem1(i,j,k)= soilvar(i,j)
      END DO
    END DO
  END DO

  ssoilvar=                                                             &
        aint2dtmp(nx,ny,nz,tem1,ipt,jpt,1,in,jn,w1,w2,                  &
                  w3,w4)

  RETURN
END SUBROUTINE colintsoil2d

!##################################################################
!##################################################################
!######                                                      ######
!######              SUBROUTINE COLINTPREC                   ######
!######                                                      ######
!######    Center for Analysis and Prediction of Storms      ######
!######    University of Oklahoma.  All rights reserved.     ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE colintprec(nx,ny,nz,tem1,sprec,ipt,jpt,in,jn,                &
           xpt,ypt,xs,ys,w1,w2,w3,w4)

!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Interpolates two dimension soil precip. in the
!  horizontal to create a column of data.  Must be called by
!  subroutine COLINTMOS.
!
!  Bilinear interpolation is used.
!
!  Based on subroutine COLINTB.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Eric Kemp
!  May 2000
!
!-----------------------------------------------------------------------
!
!  Arguments
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER :: nx,ny,nz
  REAL :: tem1(nx,ny,nz), sprec
  INTEGER :: in,jn
  REAL :: w1,w2,w3,w4
  INTEGER :: ipt,jpt
  REAL :: xs(nx),ys(nx)
  REAL :: xpt,ypt

!-----------------------------------------------------------------------
!
!  Misc. variables
!
!-----------------------------------------------------------------------

  INTEGER :: i,j,k,numprec
  REAL :: mindist,d1,d2,d3,d4

!-----------------------------------------------------------------------
!
!  Functions
!
!-----------------------------------------------------------------------

  REAL :: aint2dtmp ! function

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  numprec=0 ! # of g.p.'s surrounding stn with precip forecast
  IF (tem1(ipt,jpt,1) > 0) THEN
    numprec=numprec+1
  END IF
  IF (tem1(ipt+1,jpt,1) > 0) THEN
    numprec=numprec+1
  END IF
  IF (tem1(ipt,jpt+1,1) > 0) THEN
    numprec=numprec+1
  END IF
  IF (tem1(ipt+1,jpt+1,1) > 0) THEN
    numprec=numprec+1
  END IF
!
  IF (numprec == 0) THEN ! no g.p.'s have precip forecast
!
    sprec=0.
    GO TO 210
!
  ELSE ! at least 1 g.p. has precip forecast
!
    sprec=aint2dtmp(nx,ny,nz,tem1,ipt,jpt,1,in,jn,w1,w2,w3,w4)
!
    IF (numprec >= 3) THEN ! if at least 3 have precip,
                               ! interpolation is finished
      GO TO 210
!
    END IF
!
  END IF
!
!----------------------------------------------------------------------
!
!  If less than 3 surrounding grid points have precip
!  forecast, do further checking to see if the nearest
!  grid point has precip forecast, and use that to determine
!  whether the station's precip should be set to zero or
!  not.
!
!----------------------------------------------------------------------
!
  d1=((xpt-xs(ipt))**2+(ypt-ys(jpt))**2)**(1./2.)
  d2=((xpt-xs(ipt+1))**2+(ypt-ys(jpt))**2)**(1./2.)
  d3=((xpt-xs(ipt))**2+(ypt-ys(jpt+1))**2)**(1./2.)
  d4=((xpt-xs(ipt+1))**2+(ypt-ys(jpt+1))**2)**(1./2.)

  mindist=AMIN1(d1,d2)
  mindist=AMIN1(mindist,d3)
  mindist=AMIN1(mindist,d4)
!
  IF (mindist == d1) THEN
    IF (tem1(ipt,jpt,1) == 0) THEN
      sprec=0.
    END IF
  ELSE IF (mindist == d2) THEN
    IF (tem1(ipt+1,jpt,1) == 0) THEN
      sprec=0.
    END IF
  ELSE IF (mindist == d3) THEN
    IF (tem1(ipt,jpt+1,1) == 0) THEN
      sprec=0.
    END IF
  ELSE IF (mindist == d4) THEN
    IF (tem1(ipt+1,jpt+1,1) == 0) THEN
      sprec=0.
    END IF
  ELSE
    WRITE(6,*) "mindist-d's:",mindist,d1,d2,d3,d4
  END IF

  210  CONTINUE

  RETURN
END SUBROUTINE colintprec
