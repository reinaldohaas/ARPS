!########################################################################
!########################################################################
!######                                                            ######
!######                    PROGRAM ARPSVERIF                       ######
!######                                                            ######
!######                      Developed by                          ######
!######         Center for Analysis and Prediction of Storms       ######
!######                   University of Oklahoma                   ######
!######                                                            ######
!########################################################################
!########################################################################

PROGRAM arpsverif

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! Calculates verification statistics for ARPS forecasts.
!
!-----------------------------------------------------------------------
!
! HISTORY:
!   2002/04/03  First written by Jason J. Levit
!
!   2005/08/05  MPI version.  Kevin W. Thomas
!
!   2008/09/11  Kevin W. Thomas
!   Changes:
!
!  o  Many more stations in the station list including Oklahoma and West Texas
!     Mesonet stations.
!  o  Warn the user if the "sfcmax" array isn't large enough, meaning that the
!     entire station list isn't being read.
!  o  Fix "Makefile" so that things get rebuilt if "vericst.inc" is changed.
!     File "vericst.inc" contains "sfcmax".
!  o  The program no longer fails, leaving no HDF output file, if there are
!     missing history files (likely cause, model run failed to reach completion).
!  o  Array "model_data" is now initialized to -99.9 (aka "missing" value) instead
!     of 0.0.  This will allow the right things to happen if there are missing
!     history files.
!  o  Division by zero for RMS and BIAS calculations is prevented, with -99.9
!     returned for the case of no history files for a computation hour.
!  o  LAPS code used -100.0 instead of -99.9 for missing data.  This caused plot
!     issues with NCL code.  -100.0 changed to -99.9.
!  o  4-letter station id'are now work.  The current code would only look at
!     the first three letters of read in station files.
!  o  The program now can handle two surface files, which is what the realtime
!     system uses.  Standard SAO's go to one file, with Oklahoma and WTX Mesonet
!     going to the second.
!  o  Users can now change the naming scheme for the files.  The current code
!     looks for sao200812101200.lso.  The user can now change the "sao" and "lso"
!     part.  I did this because the realtime system uses names like
!     "okmeso200812101200.lso", which might not make sense to some people.  The
!     default is "meso200812101200.lso" (of course I'm talking about "meso" and
!     "lso", not the date part which gets computed).
!
!  12/01/2011 (Y. Wang)
!  Reorganized the program code especially the allocation statements.
!  Upgraded reflecitivity formula in his2verGrid for multi-moment microphysics.
!
!-----------------------------------------------------------------------
!
! Use modules
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
! Variable Declarations:
!
!-----------------------------------------------------------------------

  IMPLICIT NONE                    ! Force explicit declarations

  INCLUDE 'globcst.inc'
  INCLUDE 'phycst.inc'
  INCLUDE 'vericst.inc'
  INCLUDE 'grid.inc'
  INCLUDE 'mp.inc'

!
!-----------------------------------------------------------------------
!
!  Data arrays
!
!-----------------------------------------------------------------------
!
  REAL, ALLOCATABLE :: x     (:)     ! The x-coord. of the physical and
                                     ! computational grid. Defined at u-point.
  REAL, ALLOCATABLE :: y     (:)     ! The y-coord. of the physical and
                                     ! computational grid. Defined at v-point.
  REAL, ALLOCATABLE :: z     (:)     ! The z-coord. of the computational grid.
                                     ! Defined at w-point on the staggered grid.
  REAL, ALLOCATABLE :: zp    (:,:,:) ! The physical height coordinate defined at
                                     ! w-point of the staggered grid.
  REAL, ALLOCATABLE :: zpsoil(:,:,:) ! The physical height coordinate of soil
                                     ! model

!
  REAL, ALLOCATABLE :: hterain (:,:)       ! Terrain height
!
  REAL, ALLOCATABLE :: uprt   (:,:,:)    ! Perturbation u-velocity (m/s)
  REAL, ALLOCATABLE :: vprt   (:,:,:)    ! Perturbation v-velocity (m/s)
  REAL, ALLOCATABLE :: wprt   (:,:,:)    ! Perturbation w-velocity (m/s)
  REAL, ALLOCATABLE :: ptprt  (:,:,:)    ! Perturbation potential temperature (K)
  REAL, ALLOCATABLE :: pprt   (:,:,:)    ! Perturbation pressure (Pascal)
  REAL, ALLOCATABLE :: qvprt  (:,:,:)    ! Perturbation water vapor specific
                                         ! humidity (kg/kg)
  REAL, ALLOCATABLE :: qscalar(:,:,:,:)  ! Cloud water mixing ratio (kg/kg)

  REAL, ALLOCATABLE :: tke    (:,:,:)    ! Turbulent Kinetic Energy ((m/s)**2)
  REAL, ALLOCATABLE :: kmh    (:,:,:)    ! Horizontal turb. mixing coef. for
                                         ! momentum. ( m**2/s )
  REAL, ALLOCATABLE :: kmv    (:,:,:)    ! Vertical turb. mixing coef. for
                                         ! momentum. ( m**2/s )

  REAL, ALLOCATABLE :: ubar   (:,:,:)    ! Base state u-velocity (m/s)
  REAL, ALLOCATABLE :: vbar   (:,:,:)    ! Base state v-velocity (m/s)
  REAL, ALLOCATABLE :: wbar   (:,:,:)    ! Base state w-velocity (m/s)
  REAL, ALLOCATABLE :: ptbar  (:,:,:)    ! Base state potential temperature (K)
  REAL, ALLOCATABLE :: pbar   (:,:,:)    ! Base state pressure (Pascal)
  REAL, ALLOCATABLE :: rhobar (:,:,:)    ! Base state air density (kg/m**3)
  REAL, ALLOCATABLE :: qvbar  (:,:,:)    ! Base state water vapor specific
                                                                ! humidity (kg/kg)
  INTEGER :: nstyps

  INTEGER, ALLOCATABLE :: soiltyp (:,:,:)       ! Soil type
  REAL, ALLOCATABLE :: stypfrct(:,:,:)          ! Soil type
  INTEGER, ALLOCATABLE :: vegtyp  (:,:)         ! Vegetation type
  REAL, ALLOCATABLE :: lai     (:,:)            ! Leaf Area Index
  REAL, ALLOCATABLE :: roufns  (:,:)            ! Surface roughness
  REAL, ALLOCATABLE :: veg     (:,:)            ! Vegetation fraction

  REAL, ALLOCATABLE :: tsoil  (:,:,:,:)     ! soil temperature (K)
  REAL, ALLOCATABLE :: qsoil  (:,:,:,:)     ! soil moisture
  REAL, ALLOCATABLE :: wetcanp(:,:,:)       ! Canopy water amount
  REAL, ALLOCATABLE :: snowdpth(:,:)        ! Snow depth (m)

  REAL, ALLOCATABLE :: raing(:,:)         ! Grid supersaturation rain
  REAL, ALLOCATABLE :: rainc(:,:)         ! Cumulus convective rain
  REAL, ALLOCATABLE :: prcrate(:,:,:)     ! precipitation rate (kg/(m**2*s))
                                          ! prcrate(1,1,1) = total precip. rate
                                          ! prcrate(1,1,2) = grid scale precip. rate
                                          ! prcrate(1,1,3) = cumulative precip. rate
                                          ! prcrate(1,1,4) = microphysics precip. rate

  REAL, ALLOCATABLE :: radfrc(:,:,:)      ! Radiation forcing (K/s)
  REAL, ALLOCATABLE :: radsw (:,:)        ! Solar radiation reaching the surface
  REAL, ALLOCATABLE :: rnflx (:,:)        ! Net radiation flux absorbed by surface
  REAL, ALLOCATABLE :: radswnet(:,:)      ! Net shortwave radiation
  REAL, ALLOCATABLE :: radlwin(:,:)       ! Incoming longwave radiation

  REAL, ALLOCATABLE :: usflx (:,:)        ! Surface flux of u-momentum (kg/(m*s**2))
  REAL, ALLOCATABLE :: vsflx (:,:)        ! Surface flux of v-momentum (kg/(m*s**2))
  REAL, ALLOCATABLE :: ptsflx(:,:)        ! Surface heat flux (K*kg/(m*s**2))
  REAL, ALLOCATABLE :: qvsflx(:,:)        ! Surface moisture flux (kg/(m**2*s)

  REAL, ALLOCATABLE :: model_data(:,:,:)
  REAL, ALLOCATABLE :: obsrv_data(:,:,:)
                         !(:,:,1) = sfc temp
                         !(:,:,2) = sfc dewpoint
                         !(:,:,3) = sfc wind dir
                         !(:,:,4) = sfc wind speed
                         !(:,:,5) = sfc pressure

  REAL, ALLOCATABLE :: stats_data(:,:,:)
                         !(:,1,1) = RMSE data (temp)
                         !(:,1,2) = RMSE data (dewpoint)
                         !(:,1,3) = RMSE data (sfc wind dir)
                         !(:,1,4) = RMSE data (sfc wind speed)
                         !(:,1,5) = RMSE data (sfc pressure)
                         !(:,2,1) = Bias data (temp)
                         !(:,2,2) = Bias data (dewpoint)
                         !(:,2,3) = Bias data (sfc wind dir)
                         !(:,2,4) = Bias data (sfc wind speed)
                         !(:,2,5) = Bias data (sfc pressure)
                                          ! etc.
!
!-----------------------------------------------------------------------
!
!  Namelists
!
!-----------------------------------------------------------------------
!

  INTEGER :: verifopt
  NAMELIST /verif_opt/ verifopt

  INTEGER :: sndopt
  CHARACTER (LEN=256) :: sndlist, sndrunname, snddomlist
  INTEGER :: proopt
  CHARACTER (LEN=256) :: prolist, prorunname, prodomlist
  INTEGER :: sfcopt
  CHARACTER (LEN=256) :: sfclist, sfcobsdir, sfcpre, sfcpost
  CHARACTER (LEN=256) :: mesoobsdir, mesopre, mesopost
  CHARACTER (LEN=256) :: sfcrunname
  CHARACTER (LEN=256) :: blackfile
  INTEGER :: precveropt
  CHARACTER (LEN=256) :: preclist, precrunname, precdomlist
  INTEGER :: mosopt
  CHARACTER (LEN=256) :: moslist, mosrunname, mosdomlist
  INTEGER :: arpsnn_opt
  CHARACTER (LEN=256) :: nnrunname, wtsdir, nnoutputfn

  NAMELIST /single/ sndopt,sndlist,sndrunname,snddomlist,               &
                    proopt,prolist,prorunname,prodomlist,               &
                    sfcopt,sfclist,                                     &
                    sfcobsdir,sfcpre,sfcpost,                           &
                    mesoobsdir,mesopre,mesopost,                        &
                    sfcrunname,blackfile,nsfcuse,sfcuse,                &
                    precveropt,preclist,precrunname,precdomlist,        &
                    mosopt,moslist,mosrunname,mosdomlist,               &
                    arpsnn_opt,nnrunname,wtsdir,nnoutputfn

  INTEGER :: mRefopt
  CHARACTER (LEN=256) :: mReflist, mRefdir, mRefrunname
  INTEGER :: reffmt,refopt,refcomp,nreflvl
  INTEGER :: mxscorelvl
  REAL    :: threshold(20)
  INTEGER :: ibgn,iend,jbgn,jend
  NAMELIST /gridded/ mRefopt,mReflist,mRefdir,mRefrunname,reffmt,       &
                     ibgn,iend,jbgn,jend,                               &
                     refopt,refcomp,nreflvl,mxscorelvl,threshold

!-----------------------------------------------------------------------
!
! Variables for mpi jobs
!
!-----------------------------------------------------------------------
!
  INTEGER :: nprocx_in, nprocy_in   ! number of processors in input data
  INTEGER :: ncompressx, ncompressy ! compression in x and y direction:
                                    ! ncompressx=nprocx_in/nproc_x
                                    ! ncompressy=nprocy_in/nproc_y
  INTEGER :: nproc_node, readsplitin
  NAMELIST /message_passing/ nproc_x, nproc_y, max_fopen, nproc_node,   &
                             nprocx_in, nprocy_in,readsplitin
  COMMON /init1_mpi/ ncompressx, ncompressy, nproc_node

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

  INTEGER :: i,j,k
  INTEGER :: nx,ny,nz
  INTEGER :: nzsoil
  INTEGER :: ireturn, istatus
  INTEGER :: hinfmt
  INTEGER :: nhisfile_max,nhisfile
  PARAMETER (nhisfile_max=200)
  CHARACTER (LEN=256) :: grdbasfn
  CHARACTER (LEN=256) :: filename
  CHARACTER (LEN=256) :: hisfile(nhisfile_max)
  CHARACTER (LEN=4), PARAMETER :: tail=".hdf"

  INTEGER :: nf,lenfil

  INTEGER, PARAMETER :: max_dim=200


  INTEGER, PARAMETER :: unum=5 ! unit number for namelist read
  INTEGER :: istat,sd_id,sds_id
  CHARACTER (LEN=256) :: hdfname

  INTEGER :: dims(2)
  INTEGER :: lname

!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  nsfcuse = 0

  CALL mpinit_proc(0)

  IF (myproc == 0) THEN
    WRITE(6,'(/,8(5x,a,/),/)')                                          &
      '###############################################################',&
      '###############################################################',&
      '###                                                         ###',&
      '###               Welcome to ARPSVERIF                      ###',&
      '### A program that computes verification of ARPS forecasts. ###',&
      '###                                                         ###',&
      '###############################################################',&
      '###############################################################'

    sfcpre = 'sao'
    sfcpost = 'lso'
    mesoobsdir = 'NULL'
    mesopre = 'meso'
    mesopost = 'lso'

    READ (unum,message_passing)
    WRITE(6,'(a)')'Namelist block message_passing sucessfully read.'
    READ (unum,verif_opt)
    WRITE(6,'(a)')'Namelist block verif_opt sucessfully read.'
    CALL get_input_file_names(unum,hinfmt,grdbasfn,hisfile,nhisfile)
    READ (unum,single)
    WRITE(6,'(a)')'Namelist block single sucessfully read.'
    ibgn=0
    iend=0
    jbgn=0
    jend=0
    READ (unum,gridded)
    WRITE(6,'(a)')'Namelist block gridded sucessfully read.'


    readsplit(:) = 0
    IF (nprocx_in == 1 .AND. nprocy_in == 1) readsplit(:) = 1
    filename = hisfile(1)
    IF ( mp_opt > 0 .AND. readsplit(FINDX_H) <= 0 ) THEN
      CALL gtsplitfn(hisfile(1),1,1,loc_x,loc_y,1,1,0,0,1,0,filename,istatus)
    END IF
    CALL get_dims_from_data(hinfmt,TRIM(filename),                      &
                            nx,ny,nz,nzsoil,nstyps, ireturn)

    IF( ireturn /= 0 ) THEN
      PRINT*,'Problem occured when trying to get dimensions from data.'
      PRINT*,'Program stopped.'
      CALL arpsstop("Fail in get_dims_from_data",1)
    END IF

    WRITE(6,'(3(a,i5))') 'nx =',nx,', ny=',ny,', nz=',nz

    print*,'nstyps =', nstyps

  END IF

  CALL mpupdatei(nproc_x,1)
  CALL mpupdatei(nproc_y,1)
  CALL mpupdatei(readsplit,FINDX_NUM)

  CALL mpupdatei(nx,1)
  CALL mpupdatei(ny,1)
  CALL mpupdatei(nz,1)
  CALL mpupdatei(nzsoil,1)
  CALL mpupdatei(nstyps,1)

  CALL mpupdatei(nscalar, 1)
  CALL mpupdatei(nscalarq,1)
  CALL mpupdatei(P_QC,1);
  CALL mpupdatei(P_QR,1)
  CALL mpupdatei(P_QI,1)
  CALL mpupdatei(P_QS,1)
  CALL mpupdatei(P_QG,1)
  CALL mpupdatei(P_QH,1)
  CALL mpupdatei(P_NC,1)
  CALL mpupdatei(P_NR,1)
  CALL mpupdatei(P_NI,1)
  CALL mpupdatei(P_NS,1)
  CALL mpupdatei(P_NG,1)
  CALL mpupdatei(P_NH,1)
  CALL mpupdatei(P_ZR,1)
  CALL mpupdatei(P_ZI,1)
  CALL mpupdatei(P_ZS,1)
  CALL mpupdatei(P_ZG,1)
  CALL mpupdatei(P_ZH,1)

  CALL mpupdatec(qnames,nscalar*40)
  CALL mpupdatec(qdescp,nscalar*40)

  CALL mpupdatei(nhisfile,1)
  CALL mpupdatei(hinfmt,1)
  CALL mpupdatec(grdbasfn,256)
  CALL mpupdatec(hisfile,nhisfile*256)

  CALL mpupdatei(sndopt,1)
  CALL mpupdatei(proopt,1)

  CALL mpupdatei(sfcopt,1)               !sfcopt ignored, always enabled
  CALL mpupdatei(nsfcuse,1)
  CALL mpupdatec(sfcuse,nsfcuse*4)

  CALL mpupdatei(precveropt,1)
  CALL mpupdatei(mosopt,1)
  CALL mpupdatei(arpsnn_opt,1)

  CALL mpupdatei(verifopt,1)

  IF ( mp_opt == 0 ) THEN
    nproc_x = 1
    nproc_y = 1
    readsplit(:) = 0
  END IF

!-----------------------------------------------------------------------
!
!  call mpinit_var
!
!-----------------------------------------------------------------------

  IF ( mp_opt > 0 .AND. readsplit(FINDX_H) > 0 ) THEN
    nx = (nx - 3) / nproc_x + 3
    ny = (ny - 3) / nproc_y + 3
    IF ( myproc == 0 ) WRITE(6,*) 'Processor nx/ny now',nx,ny
  END IF

  CALL mpinit_var

!-----------------------------------------------------------------------
!
!  ALLOCATE arrays
!
!-----------------------------------------------------------------------

  ALLOCATE(x (nx),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:x")
  ALLOCATE(y (ny),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:y")
  ALLOCATE(z (nz),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:z")
  ALLOCATE(zp (nx,ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:zp")
  ALLOCATE(zpsoil (nx,ny,nzsoil),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:zpsoil")
  ALLOCATE(uprt (nx,ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:uprt")
  ALLOCATE(vprt (nx,ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:vprt")
  ALLOCATE(wprt (nx,ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:wprt")
  ALLOCATE(hterain (nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:hterain")
  ALLOCATE(ptprt (nx,ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:ptprt")
  ALLOCATE(pprt (nx,ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:pprt")
  ALLOCATE(qvprt (nx,ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:qvprt")
  ALLOCATE(qscalar (nx,ny,nz,nscalar),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:qscalar")
  ALLOCATE(tke (nx,ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:tke")
  ALLOCATE(kmh (nx,ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:kmh")
  ALLOCATE(kmv (nx,ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:kmv")
  ALLOCATE(ubar (nx,ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:ubar")
  ALLOCATE(vbar (nx,ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:vbar")
  ALLOCATE(wbar (nx,ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:wbar")
  ALLOCATE(ptbar (nx,ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:ptbar")
  ALLOCATE(pbar (nx,ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:pbar")
  ALLOCATE(rhobar (nx,ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:rhobar")
  ALLOCATE(qvbar (nx,ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:qvbar")
  ALLOCATE(soiltyp (nx,ny,nstyps),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:soiltyp")
  ALLOCATE(stypfrct (nx,ny,nstyps),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:stypfrct")
  ALLOCATE(vegtyp (nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:vegtyp")
  ALLOCATE(lai (nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:lai")
  ALLOCATE(roufns (nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:roufns")
  ALLOCATE(veg (nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:veg")
  ALLOCATE(tsoil (nx,ny,nzsoil,0:nstyps),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:tsoil")
  ALLOCATE(qsoil (nx,ny,nzsoil,0:nstyps),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:qsoil")
  ALLOCATE(wetcanp (nx,ny,0:nstyps),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:wetcanp")
  ALLOCATE(snowdpth (nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:snowdpth")
  ALLOCATE(raing (nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:raing")
  ALLOCATE(rainc (nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:rainc")
  ALLOCATE(prcrate (nx,ny,4),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:prcrate")
  ALLOCATE(radfrc (nx,ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:radfrc")
  ALLOCATE(radsw (nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:radsw")
  ALLOCATE(rnflx (nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:rnflx")
  ALLOCATE(radswnet (nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:radswnet")
  ALLOCATE(radlwin (nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:radlwin")
  ALLOCATE(usflx (nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:usflx")
  ALLOCATE(vsflx (nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:vsflx")
  ALLOCATE(ptsflx (nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:ptsflx")
  ALLOCATE(qvsflx (nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:qvsflx")
  ALLOCATE(tem1 (nx,ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:tem1")
  ALLOCATE(tem2 (nx,ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:tem2")
  ALLOCATE(model_data(sfcmax,nhisfile,5),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:model_data")
  ALLOCATE(obsrv_data(sfcmax,nhisfile,5),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:obsrv_data")
  ALLOCATE(stats_data(nhisfile,2,5),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:stats_data")
!
!-----------------------------------------------------------------------
!
!  Single-point verification
!
!-----------------------------------------------------------------------
!

  IF (verifopt == 1 .OR. verifopt == 3) THEN

    x = 0.0
    y  = 0.0
    z = 0.0
    zp = 0.0
    zpsoil = 0.0
    uprt = 0.0
    vprt  = 0.0
    wprt = 0.0
    hterain = 0.0
    ptprt  = 0.0
    pprt = 0.0
    qvprt = 0.0
    qscalar = 0.0
    tke = 0.0
    kmh  = 0.0
    kmv = 0.0
    ubar = 0.0
    vbar = 0.0
    wbar = 0.0
    ptbar  = 0.0
    pbar = 0.0
    rhobar = 0.0
    qvbar = 0.0
    soiltyp = 0.0
    stypfrct  = 0.0
    vegtyp = 0.0
    lai = 0.0
    roufns = 0.0
    veg = 0.0
    tsoil = 0.0
    qsoil = 0.0
    wetcanp  = 0.0
    snowdpth = 0.0
    raing = 0.0
    rainc = 0.0
    prcrate = 0.0
    radfrc = 0.0
    radsw = 0.0
    rnflx = 0.0
    radswnet = 0.0
    radlwin = 0.0
    usflx = 0.0
    vsflx = 0.0
    ptsflx = 0.0
    qvsflx = 0.0
    tem1 = 0.0
    tem2 = 0.0
    model_data=-99.9  ! was 0.0, however, missing files are now allowed
    obsrv_data=-99.9  ! Set all obs equal to an initial missing value
    stats_data=0.0

    readstns = .false.

    IF (myproc == 0) THEN
      IF (sndopt == 1) THEN
        write (6,*) "Sounding verification is currently unavailable."
        write (6,*) "Program will terminate."
        CALL arpsstop("Configuration Error",1)
      END IF

      IF (proopt == 1) THEN
        write (6,*) "Profiler verification is currently unavailable."
        write (6,*) "Program will terminate."
        CALL arpsstop("Configuration Error",1)
      END IF

      IF (precveropt == 1) THEN
        write (6,*) "Precip. verification is currently unavailable."
        write (6,*) "Program will terminate."
        CALL arpsstop("Configuration Error",1)
      END IF

      IF (mosopt == 1) THEN
        write (6,*) "MOS computation is currently unavailable."
        write (6,*) "Program will terminate."
        CALL arpsstop("Configuration Error",1)
      END IF

      IF (arpsnn_opt == 1) THEN
        write (6,*) "The ARPS Neural Network option is currently unavailable."
        write (6,*) "Program will terminate."
        CALL arpsstop("Configuration Error",1)
      END IF
    ENDIF

    CALL his2ver(nx,ny,nz,nzsoil,x,y,z,zp,zpsoil,hterain,               &
                 uprt,vprt,wprt,ptprt,pprt,qvprt,                       &
                 qscalar,tke,kmh,kmv,                                   &
                 ubar,vbar,wbar,ptbar,pbar,rhobar,qvbar,                &
                 nstyps,soiltyp,stypfrct,vegtyp,lai,roufns,             &
                 veg,tsoil,qsoil,wetcanp,snowdpth,                      &
                 raing,rainc,prcrate,radfrc,radsw,rnflx,                &
                 radswnet,radlwin,                                      &
                 usflx,vsflx,ptsflx,qvsflx,                             &
                 hinfmt,grdbasfn,hisfile,nhisfile,                      &
                 sndopt,sndlist,sndrunname,snddomlist,                  &
                 proopt,prolist,prorunname,prodomlist,                  &
                 sfcopt,sfclist,                                        &
                 sfcobsdir,sfcpre,sfcpost,                              &
                 mesoobsdir,mesopre,mesopost,                           &
                 sfcrunname,blackfile,                                  &
                 precveropt,preclist,precrunname,precdomlist,           &
                 mosopt,moslist,mosrunname,mosdomlist,                  &
                 arpsnn_opt,nnrunname,wtsdir,nnoutputfn,                &
                 model_data,obsrv_data)

    !
    ! Collect output if we're MPI.
    !

    IF (mp_opt > 0 ) THEN
       CALL verif_collect(model_data,obsrv_data,nhisfile)
    END IF

    IF (myproc == 0) THEN
      !
      ! Calculate statistics on surface data.
      !

      DO i=1,nhisfile

        ! RSME calculations
        tem1(1:sfcstn,1,1)=obsrv_data(1:sfcstn,i,1)
        tem2(1:sfcstn,1,1)=model_data(1:sfcstn,i,1)
        CALL root_mean_square(sfcstn,1,1,tem1,tem2,0,0,0,stats_data(i,1,1))
        tem1=0.0
        tem2=0.0
        tem1(1:sfcstn,1,1)=obsrv_data(1:sfcstn,i,2)
        tem2(1:sfcstn,1,1)=model_data(1:sfcstn,i,2)
        CALL root_mean_square(sfcstn,1,1,tem1,tem2,0,0,0,stats_data(i,1,2))
        tem1=0.0
        tem2=0.0
        tem1(1:sfcstn,1,1)=obsrv_data(1:sfcstn,i,3)
        tem2(1:sfcstn,1,1)=model_data(1:sfcstn,i,3)
        CALL root_mean_square(sfcstn,1,1,tem1,tem2,0,0,0,stats_data(i,1,3))
        tem1=0.0
        tem2=0.0
        tem1(1:sfcstn,1,1)=obsrv_data(1:sfcstn,i,4)
        tem2(1:sfcstn,1,1)=model_data(1:sfcstn,i,4)
        CALL root_mean_square(sfcstn,1,1,tem1,tem2,0,0,0,stats_data(i,1,4))
        tem1=0.0
        tem2=0.0
        tem1(1:sfcstn,1,1)=obsrv_data(1:sfcstn,i,5)
        tem2(1:sfcstn,1,1)=model_data(1:sfcstn,i,5)
        CALL root_mean_square(sfcstn,1,1,tem1,tem2,0,0,0,stats_data(i,1,5))

        ! Bias calculations

        tem1=0.0
        tem2=0.0
        tem1(1:sfcstn,1,1)=obsrv_data(1:sfcstn,i,1)
        tem2(1:sfcstn,1,1)=model_data(1:sfcstn,i,1)
        CALL bias_calc(sfcstn,1,1,tem1,tem2,0,0,0,stats_data(i,2,1))
        tem1=0.0
        tem2=0.0
        tem1(1:sfcstn,1,1)=obsrv_data(1:sfcstn,i,2)
        tem2(1:sfcstn,1,1)=model_data(1:sfcstn,i,2)
        CALL bias_calc(sfcstn,1,1,tem1,tem2,0,0,0,stats_data(i,2,2))
        tem1=0.0
        tem2=0.0
        tem1(1:sfcstn,1,1)=obsrv_data(1:sfcstn,i,3)
        tem2(1:sfcstn,1,1)=model_data(1:sfcstn,i,3)
        CALL bias_calc(sfcstn,1,1,tem1,tem2,0,0,0,stats_data(i,2,3))
        tem1=0.0
        tem2=0.0
        tem1(1:sfcstn,1,1)=obsrv_data(1:sfcstn,i,4)
        tem2(1:sfcstn,1,1)=model_data(1:sfcstn,i,4)
        CALL bias_calc(sfcstn,1,1,tem1,tem2,0,0,0,stats_data(i,2,4))
        tem1=0.0
        tem2=0.0
        tem1(1:sfcstn,1,1)=obsrv_data(1:sfcstn,i,5)
        tem2(1:sfcstn,1,1)=model_data(1:sfcstn,i,5)
        CALL bias_calc(sfcstn,1,1,tem1,tem2,0,0,0,stats_data(i,2,5))

      END DO

      !
      ! Write out all verification data into a specific HDF file.
      !

      lname=LEN_TRIM(sfcrunname)
      hdfname=sfcrunname(1:lname)//".hdf"
      CALL hdfopen(hdfname,2,sd_id)
      dims(2)=sfcstn
      dims(1)=len(sfcstid(1))
      CALL write_sds_char(sds_id,sd_id,'sfcstid',sfcstid,2,dims)
      CALL hdfwrt2d(model_data(:,:,1),sfcmax,nhisfile,sd_id,0,0,'model_temp', &
                 '2m ARPS Temperature','F',istatus)
      CALL hdfwrt2d(obsrv_data(:,:,1),sfcmax,nhisfile,sd_id,0,0,'obs_temp', &
                 'Observation Temperature','F',istatus)
      CALL hdfwrt2d(model_data(:,:,2),sfcmax,nhisfile,sd_id,0,0,'model_dewp', &
                 '2m ARPS Dewpoint','F',istatus)
      CALL hdfwrt2d(obsrv_data(:,:,2),sfcmax,nhisfile,sd_id,0,0,'obs_dewp', &
                 'Observation Dewpoint','F',istatus)
      CALL hdfwrt2d(model_data(:,:,3),sfcmax,nhisfile,sd_id,0,0,'model_wdir', &
                 'ARPS surface wind direction','m/s',istatus)
      CALL hdfwrt2d(obsrv_data(:,:,3),sfcmax,nhisfile,sd_id,0,0,'obs_wdir', &
                 'Observation surface wind direction','m/s',istatus)
      CALL hdfwrt2d(model_data(:,:,4),sfcmax,nhisfile,sd_id,0,0,'model_wspd', &
                 'ARPS surface wind speed','m/s',istatus)
      CALL hdfwrt2d(obsrv_data(:,:,4),sfcmax,nhisfile,sd_id,0,0,'obs_wspd', &
                 'Observation surface wind speed','m/s',istatus)
      CALL hdfwrt2d(model_data(:,:,5),sfcmax,nhisfile,sd_id,0,0,'model_pres', &
                 'ARPS surface pressure','mb',istatus)
      CALL hdfwrt2d(obsrv_data(:,:,5),sfcmax,nhisfile,sd_id,0,0,'obs_pres', &
                 'Observation surface pressure','mb',istatus)
      CALL hdfwrt1d(stats_data(:,1,1),nhisfile,sd_id,'rmse_temp', &
                 '2m Surface Temperature RMSE','none')
      CALL hdfwrt1d(stats_data(:,1,2),nhisfile,sd_id,'rmse_dewp', &
                 '2m Surface Dewpoint RMSE','none')
      CALL hdfwrt1d(stats_data(:,1,3),nhisfile,sd_id,'rmse_wdir', &
                 'Surface Wind Direction RMSE','none')
      CALL hdfwrt1d(stats_data(:,1,4),nhisfile,sd_id,'rmse_wspd', &
                 'Surface Wind Speed RMSE','none')
      CALL hdfwrt1d(stats_data(:,1,5),nhisfile,sd_id,'rmse_pres', &
                 'Surface Pressure RMSE','none')
      CALL hdfwrt1d(stats_data(:,2,1),nhisfile,sd_id,'bias_temp', &
                 '2m Surface Temperature Bias','none')
      CALL hdfwrt1d(stats_data(:,2,2),nhisfile,sd_id,'bias_dewp', &
                 '2m Surface Dewpoint Bias','none')
      CALL hdfwrt1d(stats_data(:,2,3),nhisfile,sd_id,'bias_wdir', &
                 'Surface Wind Direction Bias','none')
      CALL hdfwrt1d(stats_data(:,2,4),nhisfile,sd_id,'bias_wspd', &
                 'Surface Wind Speed Bias','none')
      CALL hdfwrt1d(stats_data(:,2,5),nhisfile,sd_id,'bias_pres', &
                 'Surface Pressure Bias','none')
      CALL hdfclose(sd_id,istat)
    END IF
  END IF

!-----------------------------------------------------------------------
!
! calculate verficiations of gridded data against NSSL reflectivity
!
!-----------------------------------------------------------------------

  IF (verifopt == 2 .OR. verifopt == 3) THEN

    WRITE(*,*) 'calculate verficiations against gridded data'

    x = 0.0
    y  = 0.0
    z = 0.0
    zp = 0.0
    zpsoil = 0.0
    uprt = 0.0
    vprt  = 0.0
    wprt = 0.0
    hterain = 0.0
    ptprt  = 0.0
    pprt = 0.0
    qvprt = 0.0
    qscalar = 0.0
    tke = 0.0
    kmh  = 0.0
    kmv = 0.0
    ubar = 0.0
    vbar = 0.0
    wbar = 0.0
    ptbar  = 0.0
    pbar = 0.0
    rhobar = 0.0
    qvbar = 0.0
    soiltyp = 0.0
    stypfrct  = 0.0
    vegtyp = 0.0
    lai = 0.0
    roufns = 0.0
    veg = 0.0
    tsoil = 0.0
    qsoil = 0.0
    wetcanp  = 0.0
    snowdpth = 0.0
    raing = 0.0
    rainc = 0.0
    prcrate = 0.0
    radfrc = 0.0
    radsw = 0.0
    rnflx = 0.0
    radswnet = 0.0
    radlwin = 0.0
    usflx = 0.0
    vsflx = 0.0
    ptsflx = 0.0
    qvsflx = 0.0
    tem1 = 0.0
    tem2 = 0.0
    model_data= 0.0
    obsrv_data=-99.9  ! Set all obs equal to an initial missing value
    stats_data=0.0

    readstns=.false.

    IF (myproc == 0) THEN
      IF (sndopt == 1) THEN
        WRITE (6,*) "Sounding verification is currently unavailable."
        WRITE (6,*) "Program will terminate."
        CALL arpsstop("Configuration Error",1)
      END IF

    ENDIF

    CALL his2verGrid(nx,ny,nz,nzsoil,x,y,z,zp,zpsoil,hterain,           &
               uprt,vprt,wprt,ptprt,pprt,qvprt,                         &
               qscalar,tke,kmh,kmv,                                     &
               ubar,vbar,wbar,ptbar,pbar,rhobar,qvbar,                  &
               nstyps,soiltyp,stypfrct,vegtyp,lai,roufns,               &
               veg,tsoil,qsoil,wetcanp,snowdpth,                        &
               raing,rainc,prcrate,radfrc,radsw,rnflx,                  &
               radswnet,radlwin,                                        &
               usflx,vsflx,ptsflx,qvsflx,                               &
               hinfmt,grdbasfn,hisfile,nhisfile,                        &
               ibgn,iend,jbgn,jend,                                     &
               mRefopt,mReflist,mRefdir,mRefrunname,                    &
              reffmt,refopt,refcomp,nreflvl,mxscorelvl,threshold)

  END IF

!-----------------------------------------------------------------------
!
! Ending of the Program
!
!-----------------------------------------------------------------------

  IF (mp_opt > 0 ) CALL mpexit(0)

  STOP
END PROGRAM arpsverif
