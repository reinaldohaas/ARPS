!
!##################################################################
!##################################################################
!######                                                      ######
!######               SUBROUTINE INITREMAPOPT                ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE initremapopt(nx,ny,nz,nzsoil,nstyps,                         &
                        nprocx_in_cmd,nprocy_in_cmd,iskip,              &
                        namelist_filename,istatus)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Initialize the control parameters for radar remapping that are read
!  from an input file.   This file shares some NAMELIST blocks with the
!  ARPS model and also has its own NAMELIST blocks.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster
!
!  12/29/2009
!
!  MODIFICATION HISTORY:
!  Keith Brewster (2011/05/01)
!  Added northazim and altoffset to NAMELIST input for mobile radar data
!
!  Keith Brewster (2012/09/07)
!  Added options to control processing and output of dual pol variables.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx       Number of grid points in the x-direction
!    ny       Number of grid points in the y-direction
!    nz       Number of grid points in the z-direction
!    nzsoil   Number of grid points in the soil profile
!
!  OUTPUT:
!
!    Control parameters defined in include files.
!
!-----------------------------------------------------------------------
!
!  Variable Declarations. (Local Variables)
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx,ny,nz      ! The number of grid points in 3 dimensions.
  INTEGER :: nzsoil        ! The number of grid points in the soil
  INTEGER :: nstyps        ! Maximum number of soil types per grid point.
  INTEGER :: nprocx_in_cmd ! Command line override for "nprocx_in".
  INTEGER :: nprocy_in_cmd ! Command line override for "nprocy_in".
  INTEGER :: iskip         ! When set, skip checking the external model
                           ! data files, for "nids2arps".
  CHARACTER(LEN=*), INTENT(IN) :: namelist_filename
  INTEGER, INTENT(OUT) :: istatus
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  REAL    :: wrmax            ! Maximun value of canopy moisture
  INTEGER :: i

  INTEGER :: nxin,nyin,nzin
  INTEGER :: nprocx_lw,nprocy_lw
  INTEGER :: lengbf,lenstr,lenfil
  INTEGER :: opstat
  INTEGER :: ireturn

  LOGICAL :: iexist      ! Flag set by inquire statement for file
                         ! existence
  REAL    :: temr
  REAL    :: dtsml0,dtsfc0        ! Temporary variable

  CHARACTER (LEN=19)  :: initime  ! Real time in form of 'year-mo-dy:hr:mn:ss'
  CHARACTER (LEN=256) :: strtmp

  INTEGER :: unum, inisplited, dmp_out_joined

  INTEGER :: nxdta, nydta

!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!
!  Global constants and parameters, most of them specify the
!  model run options.
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'phycst.inc'
!
!-----------------------------------------------------------------------
!
!  Grid and map parameters.
!
!-----------------------------------------------------------------------
!
  INCLUDE 'grid.inc'
!
!-----------------------------------------------------------------------
!
!  Radar remapping variables.
!
!-----------------------------------------------------------------------
!
  INCLUDE 'remapcst.inc'
!
!-----------------------------------------------------------------------
!
!  Message passing parameters.
!
!-----------------------------------------------------------------------
!
  INCLUDE 'mp.inc'

  INTEGER :: ncompressx, ncompressy ! compression in x and y direction:
                                    ! ncompressx=nprocx_in/nproc_x
                                    ! ncompressy=nprocy_in/nproc_y
  INTEGER :: nproc_node
!
!-----------------------------------------------------------------------
!
! Flags for fields readin from history format initial condition file
! e.g., tkein is used to determine it tke exists in history IC file.
!
!-----------------------------------------------------------------------

  INCLUDE 'indtflg.inc'
!
!-----------------------------------------------------------------------
!
!  namelist Declarations:
!
!-----------------------------------------------------------------------
!
  NAMELIST /grid_dims/ nx,ny,nz

  NAMELIST /message_passing/ nproc_x,nproc_y,max_fopen,nproc_node,      &
            nprocx_in,nprocy_in,nprocx_lw,nprocy_lw

  NAMELIST /initialization/ initime,initopt,inibasopt,viniopt,ubar0,    &
            vbar0,pt0opt,ptpert0,pt0radx,pt0rady,pt0radz,pt0ctrx,       &
            pt0ctry,pt0ctrz,rstinf,inifmt,inisplited,inifile,inigbf,    &
            sndfile,soilinitopt,soiltintv,timeopt,tsfcopt,              &
            ptground,pttrop,ttrop,htrop,qvmixed,rhmixed,mixtop,zshear

  NAMELIST /terrain/ ternopt,mntopt,hmount,mntwidx,mntwidy,             &
            mntctrx,mntctry,terndta,ternfmt

  NAMELIST /grid/ dx,dy,dz,strhopt,dzmin,zrefsfc,dlayer1,dlayer2,       &
            strhtune,zflat,ctrlat,ctrlon, crdorgnopt

  NAMELIST /projection/ mapproj, trulat1,trulat2,trulon, sclfct,        &
            mpfctopt, mptrmopt, maptest

  NAMELIST /radremapopt/ radname,radband,radfname,nlistfil,listfil,     &
            refvarname,velvarname,rhvvarname,spwvarname,snrvarname,     &
            zdpvarname,zdrvarname,zvvvarname,phivarname,kdpvarname,     &
            ref_time,ncdcfile,rad98opt,tradopt,tintrpopt,               &
            velprocopt,dlpprocopt,clrvcpopt,                            &
            unfdiagopt,qcopt,gclopt,fillref,medfilt,vadopt

  NAMELIST /radremapcst/ bmwidth,dazim,rngmin,rngmax,envavgr,           &
           northazim,altoffset,iordref,iordvel,                         &
           refchek,refmiss,velchek,velmiss,                             &
           rhvchek,rhvmiss,zdrchek,zdrmiss,kdpchek,kdpmiss,             &
           refmedl,velmedl,zdrmedl,kdpmedl,rhvmedl,refdazl,veldazl,     &
           spwthrrat,rhohvthr,snrthr,gcvrlim,                           &
           nfixed_nyqv,fixed_nyqv,elevlim_fnyqv,                        &
           vadradius,vadwidth,vadminknt

  NAMELIST /radremapout/ dirname,fntimopt,dmpfmt,hdf4cmpr,              &
           dmpzero,dualpout,ref2dopt,ref3dopt,vel2dopt,vel3dopt,        &
           zdr2dopt,zdr3dopt,kdp2dopt,kdp3dopt,rhv2dopt,rhv3dopt,       &
           wttiltopt,wtgrdtiltopt,grdtiltver,wtsoloopt,                 &
           runname,outtime

  NAMELIST /debug/ lvldbg

  COMMON /init1_mpi/ ncompressx,ncompressy,nproc_node,nprocx_lw,nprocy_lw

  REAL    :: dh
  INTEGER :: err_no = 0
  INTEGER :: idummy

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!-----------------------------------------------------------------------
!
!  Set up the default values for all the variables to be read in
!  using the namelist method. In case the user does not specify a
!  particular value, this value will be used.
!
!-----------------------------------------------------------------------
!
  istatus = 0

  nx = 67
  ny = 67
  nz = 35

  nproc_x = 1
  nproc_y = 1
  nprocx_in = 1
  nprocy_in = 1
  nprocx_lw = 1
  nprocy_lw = 1
  max_fopen = 1

  runmod=1

  dirname='./'
  runname='radremap'
  CALL gtlfnkey( runname, lfnkey )
  outtime=0.

  initime = '1977-05-20.21:00:00'
  timeopt = 0
  initopt = 1
  inibasopt  = 1
  inisplited = 0
  viniopt = 1
  pt0opt  = 0

  ubar0   = 0.0
  vbar0   = 0.0

  zshear = 3000.0
  pttrop = 343.0           ! Tropopause pot.temp.
  ttrop  = 213.0           ! Tropopause temp.
  ptground = 300.0         ! Groud urface pot.temp.
  htrop  = 12000.0         ! Tropopause height
  qvmixed= 0.015           ! Mixed layer mixing ratio
  rhmixed= 0.95
  mixtop = 1200.0          ! Mixed layer height

  ptpert0 = 0.0
  pt0radx = 10000.0
  pt0rady = 10000.0
  pt0radz =  1500.0
  pt0ctrx = 32000.0
  pt0ctry = 32000.0
  pt0ctrz =  1500.0

  rstinf  = 'may20.rst003600'
  inifmt  = 1
  inifile = 'may20.bin003600'
  inigbf  = 'may20.bingrdbas'
  sndfile = 'may20.snd'

  soilinitopt = 0
  soiltintv   = 0.0
  tsfcopt = 0

  ternopt = 0
  mntopt =  1
  hmount =  0.0
  mntwidx = 1.0E4
  mntwidy = 1.0E4
  mntctrx = 1.0E4
  mntctry = 1.0E4
  terndta ='arpstern.data'
  ternfmt = 1

  dx = 1000.0
  dy = 1000.0
  dz =  500.0

  strhopt  = 0
  dzmin    = 500.0
  zrefsfc  =   0.0
  dlayer1  =   0.0
  dlayer2  =   1.0E5
  strhtune =   1.0
  zflat    =   1.0E5

  ctrlat  =  35.0
  ctrlon  = -100.0
  crdorgnopt = 0

  mapproj  = 2
  trulat1  =   30.0
  trulat2  =   40.0
  trulon   = -100.0
  sclfct   =    1.0
  mpfctopt = 1
  mptrmopt = 1
  maptest  = 0

  moist    = 0
  mphyopt  = 0
  nmphystp = 1
  cnvctopt = 0
  subsatopt = 0
  kffbfct  = 0.0
  kfsubsattrig = 0
  ice      = 0
  wcldbs   = 0.005
  confrq   = 600.0
  qpfgfrq  = 120.0
  idownd   = 1
  dsdpref  = 0
  n0rain   = 8.0e6
  n0snow   = 3.0e6
  n0hail   = 4.0e4
  rhosnow  = 100.0
  rhohail  = 913.0

  impfallopt = 0
  fallopt = 1

  rhsat    = 0.80
  rhsatmin = 0.80
  dx_rhsatmin = 50000.
  dx_rhsat100 = 5000.

  sfcphy   = 0
  landwtr  = 1
  cdhwtropt= 0
  cdmlnd   = 3.0E-3
  cdmwtr   = 1.0E-3
  cdhlnd   = 3.0E-3
  cdhwtr   = 1.0E-3
  cdqlnd   = 2.1E-3
  cdqwtr   = 0.7E-3
  pbldopt  = 0
  pbldpth0 = 1400.0
  lsclpbl0 = 0.15
  sflxdis  = 0
  tqflxdis = 0
  dtqflxdis= 200.0
  smthflx  = 0
  numsmth  = 1
  sfcdiag  = 0

  sfcdat  = 1
  nstyp   = 4
  nstyps  = 1
  styp    = 3
  vtyp    = 10
  lai0    = 0.31
  roufns0 = 0.1
  veg0    = 0.0
  sfcdtfl = 'arpssfc.data'
  sfcfmt  = 1

  soilmodel_forced = 0
  sitemeso = '../../arpsdata.dir/mtsnorm.dir/Meso'
  siteflux = '../../arpsdata.dir/mtsnorm.dir/Flux'
  siternet = '../../arpsdata.dir/mtsnorm.dir/Radd'
  sitesoil = '../../arpsdata.dir/mtsnorm.dir/Soil'
  siteveg =  '../../arpsdata.dir/mtsnorm.dir/Veg'

  soilmodel_option = 1
  nzsoil     = 2
  dzsoil     = 1.00
  zrefsoil   = 0.0
  DO i=1,nzsoil
    tsoilint(i) = 273.15
    qsoilint(i) = 0.50
  END DO

  soilstrhopt = 0
  soildzmin   = 0.025
  soildlayer1 = 0.0
  soildlayer2 = 1.0
  soilstrhtune = 1.0

  soilinit = 1
  ptslnd0  = 300.16
  ptswtr0  = 288.16
  wetcanp0 = 0.0
  snowdpth0 = 0
  soilinfl = 'may20.soilinit'
  soilfmt  = 1

  dtsfc = 10.0

  prtsoilflx = 0

  basin = 0
  varin = 0
  mstin = 0
  rainin= 0
  prcin = 0
  icein = 0
  tkein = 0
  trbin = 0
  sfcin = 0
  landin= 0
  totin = 0
  radin = 0
  flxin = 0
  snowcin=0
  snowin= 0

  dmpfmt=1
  hdf4cmpr=2
  dmpzero=0
  dualpout=0
  ref2dopt=0
  ref3dopt=0
  vel2dopt=0
  vel3dopt=0
  zdr2dopt=0
  zdr3dopt=0
  kdp2dopt=0
  kdp3dopt=0
  rhv2dopt=0
  rhv3dopt=0

  tintrpopt=0
  velprocopt=1
  unfdiagopt=0
  qcopt=1
  gclopt = 1    ! ground clutter option
  tradopt=0
  fntimopt=0
  fillref=0
  medfilt=0
  clrvcpopt=2
!
! For radar variables large negative numbers signify missing data
! Set limit to how low valid data may reach.
!
  refchek=-20.
  refmiss=-99.
  velchek=-200.
  velmiss=-999.
  rhvchek=-0.1
  rhvmiss=-1.0
  zdrchek=-40.
  zdrmiss=-99.
  kdpchek=-40.
  kdpmiss=-99.

  winszazim=10.
  winszrad=1000.

  iordref = 2
  iordvel = 2

  refmedl = 20.
  velmedl = 15.
  zdrmedl = 5.
  kdpmedl = 20.
  rhvmedl = 0.30
  refdazl = 360.
  veldazl = 30.

  spwthrrat = 0.67
  rhohvthr  = 0.75
  snrthr = 0.0     ! dBZ
  gcvrlim = 1.6    ! m/s

  dazim = 1.0
  rngmin = 3000.
  rngmax = 230.E03
  envavgr= (rngmax/1.4142)
  bmwidth = 1.0

  northazim=0.
  altoffset=0.

  mgrid = 1

  vadopt=0
  vadradius=40.0E03
  vadwidth=1000.
  vadminknt=100

  refvarname='Reflectivity'
  velvarname='RawVelocity'
  spwvarname='SpectrumWidth'
  snrvarname='SignalToNioseRatio'
  rhvvarname='RhoHV'
  zdpvarname='Zdp'
  zdrvarname='Zdr'
  zvvvarname='Zvv'
  kdpvarname='Kdp'

  nfixed_nyqv=0
  fixed_nyqv=-999.
  elevlim_fnyqv=-99.
!
!-----------------------------------------------------------------------
!
!  Initialize message passing processors.
!  to Non-MPI defaults:
!
!  Since "mpinit_proc" has already been made, comment out the variables
!  that have already been set.
!
!-----------------------------------------------------------------------
!
! mp_opt = 0
! myproc = 0
  nprocs = 1
  loc_x = 1
  loc_y = 1
  readsplit(:) = 0
  joindmp(:)   = 1     ! by default, we want to join everything
!
!-----------------------------------------------------------------------
!
!      Initialize the processors for an MPI job.
!
!-----------------------------------------------------------------------
!
  unum = 0
  IF (myproc == 0)THEN
    WRITE(6,'(5(/3x,a))')                                                &
      'The radar remapper begins by reading a number of control ',      &
      'parameters that are specified in namelist format through the',   &
      'standard input stream (unit 5).  See the ARPS Users Guide ',     &
      'radar remapping supplements for information and guidance on ',   &
      'specifying these parameters.'

    IF (LEN_TRIM(namelist_filename) <= 0 .OR. namelist_filename == ' ') THEN
      unum = 5
      WRITE(6,'(2(1x,a,/))') 'Awaiting namelist from standard input ... ', &
                             '========================================'
    ELSE
      CALL getunit( unum )
      OPEN(unum,FILE=TRIM(namelist_filename),STATUS='OLD',FORM='FORMATTED', &
           IOSTAT=opstat)
      IF( opstat == 0) THEN
        WRITE(6,'(1x,3a,/,1x,a,/)') 'Reading namelist from file - ',   &
              TRIM(namelist_filename),' ... ',                         &
              '========================================'
      ELSE
        WRITE(6,'(1x,3a,/,1x,a,/)') 'Error opening namelistfile - ',   &
              TRIM(namelist_filename),' ... ',                         &
        '========================================'
        istatus=opstat
        err_no = err_no + 1
      END IF
    END IF

  END IF
!
!-----------------------------------------------------------------------
!
!  Read in grid dimensions.
!
!-----------------------------------------------------------------------
!
  IF (myproc == 0) THEN
    READ (unum,grid_dims, END=100)
    WRITE(6,'(a)')'Namelist block grid_dims sucessfully read.'

    WRITE(6,'(5x,a,i5)') "nx =",nx
    WRITE(6,'(5x,a,i5)') "ny =",ny
    WRITE(6,'(5x,a,i5)') "nz =",nz
  END IF

  CALL mpupdatei(nx,1)
  CALL mpupdatei(ny,1)
  CALL mpupdatei(nz,1)
!
!-----------------------------------------------------------------------
!
!  Read in message passing options.
!
!-----------------------------------------------------------------------
!
  IF (myproc == 0) THEN
    READ (unum,message_passing, END=100)
    WRITE(6,'(a)')'Namelist block message_passing sucessfully read.'

    IF (nprocx_in_cmd > 1 .OR. nprocy_in_cmd > 1 ) THEN
      nprocx_in = nprocx_in_cmd
      nprocy_in = nprocy_in_cmd
    END IF

  END IF

  CALL mpupdatei(nproc_x,1)
  CALL mpupdatei(nproc_y,1)
  CALL mpupdatei(max_fopen,1)
  CALL mpupdatei(nproc_node,1)
  CALL mpupdatei(nprocx_in,1)
  CALL mpupdatei(nprocy_in,1)
  CALL mpupdatei(nprocx_lw,1)
  CALL mpupdatei(nprocy_lw,1)

  IF (myproc == 0)THEN
    WRITE(6,'(5x,a,i4)')                                                &
         "Number of processors in the x-direction is:",nproc_x
    WRITE(6,'(5x,a,i4)')                                                &
        "Number of processors in the y-direction is:",nproc_y
    WRITE(6,'(5x,a,i4)')                                                &
        "Maximum number of files open:",max_fopen
  END IF


!
!-----------------------------------------------------------------------
!
! Note that for MP version nx & ny here are global values.  They will
! be reassigned to their per-processor value below.
!
!-----------------------------------------------------------------------
!
  IF (mp_opt > 0) THEN

    IF (nx /= nproc_x*int((nx-3)/nproc_x)+3) THEN
      nx = nproc_x*int((nx-3)/nproc_x+0.9999999999999) + 3
      IF (myproc == 0) THEN
        WRITE (6,*) "NOTICE: adjusting nx to fit on ",nproc_x," processors:"
        WRITE(6,'(5x,a,i5)') "   new nx =",nx
      ENDIF
    ENDIF
    IF (ny /= nproc_y*int((ny-3)/nproc_y)+3) THEN
      ny = nproc_y*int((ny-3)/nproc_y+0.9999999999999) + 3
      IF (myproc == 0) THEN
        WRITE (6,*) "NOTICE: adjusting ny to fit on ",nproc_y," processors:"
        WRITE(6,'(5x,a,i5)') "   new ny =",ny
      ENDIF
    ENDIF

  ELSE

    nproc_x = 1
    nproc_y = 1
    nprocs = 1
    max_fopen = 1

  ENDIF
!
!-----------------------------------------------------------------------
!
!  Initialize message passing variables.
!
!-----------------------------------------------------------------------
!
  CALL mpinit_var

  IF(nprocs < 2) mp_opt=0
!
!-----------------------------------------------------------------------
!
!  Read in control parameters for model initialization
!
!-----------------------------------------------------------------------
!
  IF (myproc == 0) THEN
    READ (unum,initialization,END=100)
    WRITE(6,'(a)')'Namelist block initialization sucessfully read.'
  END IF

  DO i = FINDX_NUM, 1, -1
    idummy = inisplited / 10**(i-1)     ! ith digit is odd (1,3,5..), no readsplit
    readsplit(i) = MOD(idummy+1, 2)     ! ith digit is even(0,2,4..), readsplit
  END DO
!
!
!-----------------------------------------------------------------------
!
! Allow the user to use non-split sfc/trn files even though the ARPS
! history files are split.
!
!-----------------------------------------------------------------------
!
  CALL mpupdatec(initime,19)
  CALL mpupdatei(initopt,1)
  CALL mpupdatei(inibasopt,1)
  CALL mpupdatei(viniopt,1)
  CALL mpupdater(ubar0,1)
  CALL mpupdater(vbar0,1)

  CALL mpupdater(zshear,1)
  CALL mpupdater(pttrop,1)
  CALL mpupdater(ttrop,1)
  CALL mpupdater(ptground,1)
  CALL mpupdater(htrop,1)
  CALL mpupdater(qvmixed,1)
  CALL mpupdater(rhmixed,1)
  CALL mpupdater(mixtop,1)

  CALL mpupdatei(pt0opt,1)
  CALL mpupdater(ptpert0,1)
  CALL mpupdater(pt0radx,1)
  CALL mpupdater(pt0rady,1)
  CALL mpupdater(pt0radz,1)
  CALL mpupdater(pt0ctrx,1)
  CALL mpupdater(pt0ctry,1)
  CALL mpupdater(pt0ctrz,1)
  CALL mpupdatec(rstinf,256)
  CALL mpupdatei(inifmt,1)
  CALL mpupdatei(readsplit,FINDX_NUM)
  CALL mpupdatec(inifile,256)
  CALL mpupdatec(inigbf, 256)
  CALL mpupdatec(sndfile,256)
  CALL mpupdatei(soilinitopt,1)
  CALL mpupdater(tsfcopt,1)
  CALL mpupdater(soiltintv,1)
  CALL mpupdatei(timeopt,1)

  READ (initime, '(i4.4,1x,i2.2,1x,i2.2,1x,i2.2,1x,i2.2,1x,i2.2)' )     &
       year,month,day,hour,minute,second

  IF (myproc == 0)THEN
    WRITE(6,'(a,a,i4.4,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i2.2)')            &
        '     The initial local time for this run is ',                 &
        '     year-mo-dy:hr:mn:ss = ',                                  &
        year,'-',month,'-',day,'.',hour,':',minute,':',second

  END IF
!
!-----------------------------------------------------------------------
!
!  Input data files for initialization:
!
!-----------------------------------------------------------------------
!

  IF (iskip /= 1) THEN
    IF (myproc == 0 .AND. initopt == 2) THEN
      WRITE(6,'(5x,a,a)')                                                 &
          'The two time level restart data to be read in is ',            &
          TRIM(rstinf)

      WRITE(6,'(5x,a,i4)')                                                &
          'The history dump type restart data format was ',inifmt
    END IF

    IF (myproc == 0 .AND. initopt == 3) THEN
      WRITE(6,'(5x,a,a)')                                                 &
        'The time-dependent history dump format data to be read is ',&
        TRIM(inifile)

      WRITE(6,'(5x,a,a,a)')                                               &
          'The base state/grid history dump ',                            &
          'format restart data to be read is ', TRIM(inigbf)

      IF(mp_opt > 0 .AND. readsplit(FINDX_H) <= 0) THEN
        CALL gtsplitfn(inifile,1,1,loc_x,loc_y,1,1,                      &
                       nprocx_lw-1,nprocy_lw-1,1,lvldbg,strtmp,ireturn)
      ELSE IF (mp_opt == 0 .AND. (nprocx_in > 1 .OR. nprocy_in > 1) ) THEN
        CALL gtsplitfn(inifile,1,1,loc_x,loc_y,1,1,                      &
                       nprocx_lw-1,nprocy_lw-1,1,lvldbg,strtmp,ireturn)
      ELSE
        strtmp = inifile
      END IF

      CALL get_dims_from_data(inifmt,strtmp,                    &
                              nxin,nyin,nzin,nzsoil,nstyps, ireturn)

      IF( ireturn /= 0 ) THEN
        WRITE(6,'(a)') 'Problem occured when trying to get dimensions from data.'
        WRITE(6,'(a)') 'Program stopped.'
        STOP
      END IF

      nxdta = nx
      nydta = ny

      IF (readsplit(FINDX_H) <= 0) THEN   ! no split on-the-fly, so should read split files.
        nxdta = (nx-3)/nprocx_in + 3
        nydta = (ny-3)/nprocy_in + 3
      ENDIF

      IF(nxin /= nxdta .OR. nyin /= nydta .OR. nzin /= nz) THEN
        WRITE(6,'(a)')                                                       &
         'Dimension mis-match between specified and initial background file.'
        WRITE(6,'(a,i6,a,i6)') ' Specified nx: ',nxdta,'  background: ',nxin
        WRITE(6,'(a,i6,a,i6)') ' Specified ny: ',nydta,'  background: ',nyin
        WRITE(6,'(a,i6,a,i6)') ' Specified nz: ',nz,'  background: ',nzin
  !     CALL ARPSSTOP(' dimension mis-match from get_dims_from_data',1)
      END IF

      WRITE(6,'(a,i6,a,i6)')                                                 &
        ' Found nzsoil = ',nzsoil,' and nstyps = ',nstyps
      nstyp=nstyps

    END IF
  END IF ! iskip

  CALL mpupdatei(nx,1)
  CALL mpupdatei(ny,1)
  CALL mpupdatei(nz,1)
  CALL mpupdatei(nzsoil,1)
  CALL mpupdatei(nstyp,1)
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
!
!-----------------------------------------------------------------------
!
!  Input the environmental sounding.
!
!-----------------------------------------------------------------------
!
  IF (myproc == 0 .AND. initopt == 1)THEN
    WRITE(6,'(5x,a,a)')                                                 &
         'Sounding file to be used is ', TRIM(sndfile)
  END IF
!
!-----------------------------------------------------------------------
!
!  Specify the terrain:
!
!-----------------------------------------------------------------------
!
  IF (myproc == 0)THEN
    READ (unum,terrain,END=100)
    WRITE(6,'(a)')'Namelist block terrain sucessfully read.'
  END IF

  CALL mpupdatei(ternopt,1)
  CALL mpupdatei(mntopt, 1)
  CALL mpupdater(hmount, 1)
  CALL mpupdater(mntwidx,1)
  CALL mpupdater(mntwidy,1)
  CALL mpupdater(mntctrx,1)
  CALL mpupdater(mntctry,1)
  CALL mpupdatec(terndta,256)
  CALL mpupdatei(ternfmt,1)

  IF (myproc == 0 .AND. initopt == 1) THEN

    WRITE(6,'(5x,a,i4)') 'The mountain type option was ', mntopt

    WRITE(6,'(5x,a,f10.3,a)')                                           &
        'The height of mountain is ', hmount,' (m).'

    WRITE(6,'(5x,a,2e10.3,a,/5x,a,a,/5x,a,e10.3,a,e10.3,a)')            &
        'The input half-width of bell-shaped mountain are ',            &
        mntwidx, mntwidy, ' (m)',                                       &
        'in x and y direction respectively, and the center is ',        &
        'located at','x=',mntctrx,' y=',mntctry,' (m).'
    WRITE(6,'(5x,a,i4)') 'The terrain option was ', ternopt

    WRITE(6,'(5x,a,a)')                                              &
      'The terrain data file is ',TRIM(terndta)
    WRITE(6,'(5x,a,i4)') 'The terrain data file format is ', ternfmt
  END IF
!
!-----------------------------------------------------------------------
!
!  Input horizontal grid size
!
!-----------------------------------------------------------------------
!
  IF (myproc == 0) THEN
    READ (unum,grid,END=100)
    WRITE(6,'(a)')'Namelist block grid sucessfully read.'
  END IF

  CALL mpupdater(dx,1)
  CALL mpupdater(dy,1)
  CALL mpupdater(dz,1)
  CALL mpupdatei(strhopt,1)
  CALL mpupdater(dzmin,1)
  CALL mpupdater(zrefsfc,1)
  CALL mpupdater(dlayer1,1)
  CALL mpupdater(dlayer2,1)
  CALL mpupdater(strhtune,1)
  CALL mpupdater(zflat,1)
  CALL mpupdater(ctrlat,1)
  CALL mpupdater(ctrlon,1)
  CALL mpupdatei(crdorgnopt,1)

  IF( strhopt == 0.AND.dzmin /= dz ) THEN
    dzmin = dz
    IF (myproc == 0)  WRITE(6,'(5x,a)')                               &
         'For non-stretched case, dzmin was reset to dz.'
  END IF

  IF (myproc == 0)THEN
    WRITE(6,'(5x,a,f10.3,a)')                                         &
       'Input dx was',dx,' meters'

    WRITE(6,'(5x,a,f10.3,a)')                                         &
       'Input dy was ',dy,' meters'

    WRITE(6,'(5x,a,i4)') 'The stretch option was ', strhopt

    WRITE(6,'(5x,a,f10.3,a)')                                         &
       'Input dz was ',dz,' meters'

    WRITE(6,'(5x,a,f10.3,a)')                                         &
       'Input ctrlat was ',ctrlat,' degrees North'

    WRITE(6,'(5x,a,f10.3,a)')                                         &
       'Input ctrlon was ',ctrlon,' degrees East'

    WRITE(6,'(5x,a,i4)') 'domain origin option was', crdorgnopt

    WRITE(6,'(5x,a,f10.3,a)')                                         &
       'dzmin is ',dzmin,' meters'

    WRITE(6,'(5x,a,f10.3,a)')                                         &
       'zrefsfc is ',zrefsfc ,' meters'

    WRITE(6,'(5x,a,f10.3,a)')                                         &
       'dlayer1 is ',dlayer1 ,' meters'

    WRITE(6,'(5x,a,f10.3,a)')                                         &
       'dlayer2 is ',dlayer2 ,' meters'

    WRITE(6,'(5x,a,f10.3,a)')                                         &
       'zflat   is ',zflat   ,' meters'
  END IF
!
!-----------------------------------------------------------------------
!
!  Input map projection parameters
!
!-----------------------------------------------------------------------
!
  IF (myproc == 0) THEN
    READ (unum,projection,END=100)
    WRITE(6,'(a)')'Namelist block projection sucessfully read.'
  END IF

  CALL mpupdatei(mapproj,1)
  CALL mpupdater(trulat1,1)
  CALL mpupdater(trulat2,1)
  CALL mpupdater(trulon,1)
  CALL mpupdater(sclfct,1)
  CALL mpupdatei(mpfctopt,1)
  CALL mpupdatei(mptrmopt,1)
  CALL mpupdatei(maptest,1)

  IF (myproc == 0)THEN
    WRITE(6,'(5x,a,i4)')                                            &
         'Input mapproj was ',mapproj

    WRITE(6,'(5x,a,f10.3,a)')                                       &
         'Input trulat1 was ',trulat1,' degrees North'

    WRITE(6,'(5x,a,f10.3,a)')                                       &
         'Input trulat2 was ',trulat2,' degrees North'

    WRITE(6,'(5x,a,f10.3)')                                         &
        'The latitude of the center of the model domain was ',ctrlat

    WRITE(6,'(5x,a,f10.3)')                                         &
        'The longitude of the center of the model domain was ',ctrlon

    WRITE(6,'(5x,a,f10.3,a)')                                       &
         'Input trulon was ',trulon,' degrees East'

    WRITE(6,'(5x,a,e15.5)')                                         &
         'Input sclfct was ',sclfct

    WRITE(6,'(5x,a,i5)')                                            &
         'The option for map factor was ', mpfctopt

    WRITE(6,'(5x,a,i5)')                                            &
       'The option for map factor term in u and v advection was ',  &
       mptrmopt
  END IF
!
!-----------------------------------------------------------------------
!
!  Input remap control paramaters:
!
!-----------------------------------------------------------------------
!
  rad98opt = 0
  IF (myproc == 0) THEN
    READ (unum,radremapopt,END=100)
    WRITE(6,'(a)')'Namelist block radremapopt sucessfully read.'
  END IF

  CALL mpupdatec(radname,4)
  CALL mpupdatec(radfname,256)

  CALL mpupdatec(refvarname,80)
  CALL mpupdatec(velvarname,80)
  CALL mpupdatec(spwvarname,80)
  CALL mpupdatec(snrvarname,80)
  CALL mpupdatec(rhvvarname,80)
  CALL mpupdatec(zdpvarname,80)
  CALL mpupdatec(zdrvarname,80)
  CALL mpupdatec(zvvvarname,80)
  CALL mpupdatec(kdpvarname,80)

  CALL mpupdatei(radband,1)
  CALL mpupdatei(tintrpopt,1)
  CALL mpupdatei(qcopt,1)
  CALL mpupdatei(unfdiagopt,1)
  CALL mpupdatei(gclopt,1)
  CALL mpupdatei(velprocopt,1)
  CALL mpupdatei(fillref,1)
  CALL mpupdatei(medfilt,1)
  CALL mpupdatei(tradopt,1)
  CALL mpupdatei(rad98opt,1)

  CALL mpupdatei(nfixed_nyqv,1)
  CALL mpupdater(fixed_nyqv,maxfixed_nyqv)
  CALL mpupdater(elevlim_fnyqv,maxfixed_nyqv)

  IF (myproc == 0) THEN
    READ (unum,radremapcst,END=100)
    WRITE(6,'(a)')'Namelist block remapcst sucessfully read.'
  END IF

  CALL mpupdater(bmwidth,1)
  CALL mpupdater(dazim,1)
  CALL mpupdater(rngmin,1)
  CALL mpupdater(rngmax,1)
  CALL mpupdater(envavgr,1)
  CALL mpupdater(northazim,1)
  CALL mpupdater(altoffset,1)
  CALL mpupdater(refchek,1)
  CALL mpupdater(refmiss,1)
  CALL mpupdater(velchek,1)
  CALL mpupdater(velmiss,1)
  CALL mpupdater(refmedl,1)
  CALL mpupdater(velmedl,1)
  CALL mpupdater(zdrmedl,1)
  CALL mpupdater(kdpmedl,1)
  CALL mpupdater(rhvmedl,1)
  CALL mpupdater(refdazl,1)
  CALL mpupdater(veldazl,1)
  CALL mpupdater(spwthrrat,1)
  CALL mpupdater(rhohvthr,1)
  CALL mpupdater(snrthr,1)
  CALL mpupdater(gcvrlim,1)

  IF (myproc == 0) THEN
    READ (unum,radremapout,END=100)
    WRITE(6,'(a)')'Namelist block radremapout sucessfully read.'
  END IF

  CALL mpupdatec(dirname,256)
  CALL mpupdatei(fntimopt,1)
  CALL mpupdatei(dmpfmt,1)
  CALL mpupdatei(hdf4cmpr,1)
  CALL mpupdatei(ref2dopt,1)
  CALL mpupdatei(ref3dopt,1)
  CALL mpupdatei(vel2dopt,1)
  CALL mpupdatei(vel3dopt,1)
  CALL mpupdatei(zdr2dopt,1)
  CALL mpupdatei(zdr3dopt,1)
  CALL mpupdatei(kdp2dopt,1)
  CALL mpupdatei(kdp3dopt,1)
  CALL mpupdatei(rhv2dopt,1)
  CALL mpupdatei(rhv3dopt,1)
  CALL mpupdatec(runname,80)
  CALL mpupdater(outtime,1)

  CALL gtlfnkey( runname, lfnkey )
!
!-----------------------------------------------------------------------
!
!  Input debug information print controls:
!
!-----------------------------------------------------------------------
!
  IF (myproc == 0) THEN
    READ (unum,debug,END=100)
    WRITE(6,'(a)')'Namelist block debug sucessfully read.'
  END IF
  CALL mpupdatei(lvldbg,1)

  IF (myproc == 0) &
  WRITE(6,'(5x,a,i4)') 'The debug printing level was ', lvldbg

!-----------------------------------------------------------------------
!
! Hardcode some parameters
! (they have been removed from the namelist)
!
!----------------------------------------------------------------------

  fallvalpha = 0.5
  alfcoef    = 0.5
  IF ( trbvimp == 0 ) alfcoef = 1.0

!
! Hard-code dsdpref because we do not have namelist microphysics with
! radar remapping.
!

  dsdpref = 1
!
!-----------------------------------------------------------------------
!
!  Adjust max_fopen if readsplit > 0 and/or joindmp > 0
!
!-----------------------------------------------------------------------

  IF (mp_opt > 0 ) THEN

    readstride = max_fopen
    IF (ANY(readsplit > 0) ) readstride = nprocs

    dumpstride = max_fopen
    IF (ANY(joindmp > 0) )  dumpstride = nprocs

  ELSE
    readsplit(:)  = 0
    joindmp(:)    = 0
    readstride = max_fopen
    dumpstride = max_fopen
  END IF
!
!-----------------------------------------------------------------------
!
!  Compute derived variables.
!
!-----------------------------------------------------------------------
!
  IF (mp_opt > 0) THEN  ! Convert from global to processor specific values.
    nx = (nx - 3)/nproc_x + 3
    ny = (ny - 3)/nproc_y + 3
    IF (myproc == 0) THEN
      WRITE(6,'(5x,a,i5)') "Processor nx =",nx
      WRITE(6,'(5x,a,i5)') "Processor ny =",ny
    END IF
  END IF

  IF( initopt == 2 .or.initopt == 4 ) THEN
    restrt = 1
  ELSE
    restrt = 0
  END IF

  CALL julday( year, month, day, jday )         ! Get the Julian day

  IF ( mapproj == 0 ) mpfctopt = 0
  mptrmopt = mptrmopt * mpfctopt
  latitud = ctrlat
  longitud= ctrlon
  IF ( mapproj == 0 ) THEN
    trulat1 = ctrlat
    trulat2 = ctrlat
    trulon  = ctrlon
  END IF

  dxinv = 1.0/dx
  dyinv = 1.0/dy
  dzinv = 1.0/dz
  dzsoilinv = 1.0/dzsoil
  dzsoilinv2 = dzsoilinv * dzsoilinv

  xl = (nx-3)*dx
  yl = (ny-3)*dy
  zh = (nz-3)*dz

  IF( ternopt /= 0 .OR. strhopt /= 0 ) THEN
    crdtrns = 1
  ELSE
    crdtrns = 0
  END IF

  IF( runmod == 2 ) THEN
    dh = dx
  ELSE IF( runmod == 3 ) THEN
    dh = dy
  ELSE
    dh = SQRT(dx*dy)
  END IF

  IF (unum /= 5 .AND. myproc == 0) THEN
    CLOSE(unum)
    CALL retunit(unum)
  END IF

  IF( err_no /= 0 ) THEN
    IF (myproc == 0) WRITE(6,'(5x,i4,a,/5x,a)')                         &
        err_no, ' fatal errors found with the input parameters.',       &
        'Please check the input parameters carefully.'
    CALL arpsstop('arpsstop called from INITREMAPOPT',1)
  END IF

  RETURN

  100 CONTINUE
  IF (myproc == 0) WRITE(6,'(1x,a,/1x,a)')                              &
        'End of file found while reading input parameters.',            &
        'Please check the input file carefully.'
  CALL arpsstop('arpsstop called from INITREMAPOPT',1)

  RETURN
END SUBROUTINE initremapopt
!
!##################################################################
!##################################################################
!######                                                      ######
!######              SUBROUTINE MKRADFNM                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE mkradfnm(dmpfmt,dir,ldir,radar,iyr,imo,ida,ihr,imin,isec,    &
                    fname,lfnm)
!
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
!  Input arguments
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=256) :: dir
  INTEGER :: ldir
  CHARACTER (LEN=4) :: radar
  INTEGER :: iyr
  INTEGER :: imo
  INTEGER :: ida
  INTEGER :: ihr
  INTEGER :: imin
  INTEGER :: isec
  INTEGER :: dmpfmt
!
!-----------------------------------------------------------------------
!
!  Output arguments
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=256) :: fname
  INTEGER :: lfnm
!
  INTEGER :: rdtime,jyr,jmo,jda,jhr,jmin,jsec
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  CALL ctim2abss( iyr,imo,ida,ihr,imin,isec, rdtime )
  !  Round to nearest minute
  IF(isec >= 30) rdtime=rdtime+(60-isec)
  CALL abss2ctim( rdtime, jyr, jmo, jda, jhr, jmin, jsec )
  IF(dmpfmt==1)THEN
    WRITE(fname,'(a,a4,a1,i4.4,2(i2.2),a1,2(i2.2))')                    &
          dir(1:ldir),radar,'.',                                        &
          jyr,jmo,jda,'.',jhr,jmin
    lfnm=ldir+16
  ELSE
    WRITE(fname,'(a,a4,a1,i4.4,2(i2.2),a1,2(i2.2),a5)')                 &
          dir(1:ldir),radar,'.',                                        &
          jyr,jmo,jda,'.',jhr,jmin,'.hdf4'
    lfnm=ldir+21
  END IF

  RETURN
END SUBROUTINE mkradfnm

!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE WRTRADCOL                  ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE wrtradcol(nx,ny,nxny,nz,ncolp,ncoltot,                       &
                  dmpfmt,iradfmt,hdf4cmpr,dmpzero,dualpol,              &
                  rfname,radid,radlat,radlon,radalt,                    &
                  iyr,imon,iday,ihr,imin,isec,vcpnum,isource,           &
                  refelvmin,refelvmax,refrngmin,refrngmax,              &
                  icolp,jcolp,xcolp,ycolp,zcolp,                        &
                  colref,colvel,colnyq,coltim,colrhv,colzdr,colkdp,     &
                  istatus)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Writes gridded radar data to a file as columns with
!  individual lat,lons.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster
!  06/22/95
!
!  MODIFICATION HISTORY:
!
!  2000/09/11 (Gene Bassett)
!  Use only reflectivity to accept or reject a column (thus allowing one
!  to output nids columns without processing velocity data).
!
!  04/29/02 Leilei Wang and Keith Brewster
!  Added hdf option, including two new variables in the argument list.
!
!  02/02/04 Keith Brewster
!  Modified to replace Time and Nyquist values with missing at points
!  where vel and ref are missing for hdf write.  Should save compressed
!  file space.
!
!  02/08/2010 Keith Brewster
!  New revised MPI handling.
!
!  03/07/2010 Keith Brewster
!  Revision to handle column-based data input.
!  Changed name to wrtradcol to avoid confusion with previous version.
!
!  04/25/2011 Yunheng Wang
!  Simplified the code significantly and fixed the scalarbility issue
!  with too many messages that overflows the MPI buffer when the number
!  of processes is large.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!    dmpfmt    file format (1:binary, 3:hdf)
!    iradfmt   binary format
!    hdf4cmpr  hdf4 compression level
!    rfname    radar file name (character*80)
!    radid     radar id (character*4)
!    radlat    latitude of radar (degrees N)
!    radlon    longitude of radar (degrees E)
!    radalt    elevation of radar (m MSL)
!    iyr       year
!    imon      month
!    iday      day
!    ihr       hour
!    imin      min
!    isec      sec
!    vcpnum    VCP (scan type) number
!    isource)  source number
!                1= WSR-88D raw
!                2= WSR-88D NIDS
!
!  OUTPUT:
!    data are written to file
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
  INTEGER, INTENT(IN) :: nx,ny,nxny,nz
  INTEGER, INTENT(IN) :: ncolp,ncoltot
!
  INTEGER, INTENT(IN) :: dmpfmt
  INTEGER, INTENT(IN) :: iradfmt
  INTEGER, INTENT(IN) :: hdf4cmpr
  INTEGER, INTENT(IN) :: dmpzero
  INTEGER, INTENT(IN) :: dualpol
  CHARACTER (LEN=256), INTENT(IN) :: rfname
  CHARACTER (LEN=4),   INTENT(IN) :: radid
  REAL,    INTENT(IN)  :: radlat
  REAL,    INTENT(IN)  :: radlon
  REAL,    INTENT(IN)  :: radalt
  REAL,    INTENT(IN)  :: refelvmin
  REAL,    INTENT(IN)  :: refelvmax
  REAL,    INTENT(IN)  :: refrngmin
  REAL,    INTENT(IN)  :: refrngmax
  INTEGER, INTENT(IN)  :: iyr,imon,iday,ihr,imin,isec
  INTEGER, INTENT(IN)  :: vcpnum
  INTEGER, INTENT(IN)  :: isource
!
  INTEGER, INTENT(IN)  :: icolp(nxny)
  INTEGER, INTENT(IN)  :: jcolp(nxny)
  REAL, INTENT(IN)     :: xcolp(nxny)
  REAL, INTENT(IN)     :: ycolp(nxny)
  REAL, INTENT(IN)     :: zcolp(nz,nxny)
  REAL, INTENT(IN)     :: colref(nz,nxny)
  REAL, INTENT(IN)     :: colvel(nz,nxny)
  REAL, INTENT(INOUT)  :: colnyq(nz,nxny)
  REAL, INTENT(INOUT)  :: coltim(nz,nxny)
  REAL, INTENT(IN)     :: colrhv(nz,nxny)
  REAL, INTENT(IN)     :: colzdr(nz,nxny)
  REAL, INTENT(IN)     :: colkdp(nz,nxny)

  INTEGER, INTENT(OUT) :: istatus
!
!-----------------------------------------------------------------------
!
!  Radar output variables
!
!-----------------------------------------------------------------------
!
  INTEGER, ALLOCATABLE :: outi(:)
  INTEGER, ALLOCATABLE :: outj(:)
  INTEGER, ALLOCATABLE :: outk(:,:)
  INTEGER, ALLOCATABLE :: outnumlev(:)

  REAL,    ALLOCATABLE :: outx(:)
  REAL,    ALLOCATABLE :: outy(:)
  REAL,    ALLOCATABLE :: outlat(:)
  REAL,    ALLOCATABLE :: outlon(:)
  REAL,    ALLOCATABLE :: outsfc(:)

  REAL,    ALLOCATABLE :: outhgt(:,:)
  REAL,    ALLOCATABLE :: outref(:,:)
  REAL,    ALLOCATABLE :: outvel(:,:)
  REAL,    ALLOCATABLE :: outnyq(:,:)
  REAL,    ALLOCATABLE :: outtim(:,:)
  REAL,    ALLOCATABLE :: outrhv(:,:)
  REAL,    ALLOCATABLE :: outzdr(:,:)
  REAL,    ALLOCATABLE :: outkdp(:,:)

  INTEGER, ALLOCATABLE :: outk1d(:)
  REAL,    ALLOCATABLE :: outhgt1d(:)
  REAL,    ALLOCATABLE :: outref1d(:)
  REAL,    ALLOCATABLE :: outvel1d(:)
  REAL,    ALLOCATABLE :: outnyq1d(:)
  REAL,    ALLOCATABLE :: outtim1d(:)
  REAL,    ALLOCATABLE :: outrhv1d(:)
  REAL,    ALLOCATABLE :: outzdr1d(:)
  REAL,    ALLOCATABLE :: outkdp1d(:)
!
!-----------------------------------------------------------------------
!
!  Include files
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'grid.inc'
  INCLUDE 'remap.inc'
  INCLUDE 'mp.inc'
!
!-----------------------------------------------------------------------
!
!  Radar output descriptors
!
!-----------------------------------------------------------------------
!
  INTEGER, PARAMETER :: mxradvr=10
  INTEGER, PARAMETER :: iradvr(mxradvr) = (/1,2,3,4,5,6,0,0,0,0/)
!
!-----------------------------------------------------------------------
!
!  Radar output thresholds
!
!-----------------------------------------------------------------------
!
  REAL, PARAMETER :: refmin=-30.0
  REAL, PARAMETER :: refmax=80.
  REAL, PARAMETER :: velmin=-300.
  REAL, PARAMETER :: velmax=300.
  REAL, PARAMETER :: misval=-999.0
!
!-----------------------------------------------------------------------
!
!  Misc local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: iunit,itime
  INTEGER :: idat,jdat,k,klev,kk,jproc
  INTEGER :: tagbase,itag
  INTEGER :: kntcol,kntcolprc,numcol,numcolprc,nn
  INTEGER :: nxlg,nylg,numlev
  INTEGER :: idummy,ncolj,nsize,nradvr
  INTEGER :: istat,ierr,sd_id
  INTEGER :: irngmin,irngmax
  INTEGER :: ioffset,joffset
  REAL :: xmin,xmax,ymin,ymax
  REAL :: collat,collon,elev,rdummy
  LOGICAL :: dualpdata
  INTEGER, PARAMETER :: typelev = 1
  INTEGER, PARAMETER :: tagmult = 10000
!
! Temporary transfer arrays for binary MPI
!
  INTEGER, ALLOCATABLE :: icolp_tmp(:)
  INTEGER, ALLOCATABLE :: jcolp_tmp(:)
  REAL, ALLOCATABLE :: xcolp_tmp(:)
  REAL, ALLOCATABLE :: ycolp_tmp(:)
  REAL, ALLOCATABLE :: zcolp_tmp(:,:)
  REAL, ALLOCATABLE :: colref_tmp(:,:)
  REAL, ALLOCATABLE :: colvel_tmp(:,:)
  REAL, ALLOCATABLE :: colnyq_tmp(:,:)
  REAL, ALLOCATABLE :: coltim_tmp(:,:)
  REAL, ALLOCATABLE :: colrhv_tmp(:,:)
  REAL, ALLOCATABLE :: colzdr_tmp(:,:)
  REAL, ALLOCATABLE :: colkdp_tmp(:,:)
!
! Temporary gather arrays for HDF MPI
!
  LOGICAL, ALLOCATABLE :: havdat(:)
  INTEGER, ALLOCATABLE :: itmp1dlg(:)
  INTEGER, ALLOCATABLE :: itmp2dlg(:,:)
  REAL, ALLOCATABLE :: rtmp1dlg(:)
  REAL, ALLOCATABLE :: rtmp2dlg(:,:)

  INTEGER, ALLOCATABLE :: disp1d(:), disp2d(:)       ! WYH
  INTEGER, ALLOCATABLE :: lsize1d(:), lsize2d(:)
  INTEGER  :: ncolp2d, ncoltot2d, iloc

  LOGICAL :: fileopnd

  INTEGER(KIND=selected_int_kind(4)), ALLOCATABLE  :: itmp(:,:)

  REAL :: vmax,vmin
  REAL :: rmax,rmin

!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  print *, ' inside wrtradcol, dualpol =',dualpol
  dualpdata = (dualpol > 0)
  fileopnd=.FALSE.
  kntcol=0
  numcol=0
  numcolprc=0
  idummy=-999
  rdummy=-999.
  nxlg=(nproc_x*(nx-3))+3
  nylg=(nproc_y*(ny-3))+3
  IF (lvldbg > 0) THEN
    IF ( myproc == 0 ) print *, ' nxlg = ',nxlg,'   nylg=',nylg
    print *, ' myproc = ',myproc,' mp_opt = ',mp_opt,'  nprocs = ',nprocs
  END IF

  IF(nprocs > 1 .AND. myproc == 0) THEN
    IF (dmpfmt == 1) THEN ! binary
      ALLOCATE (icolp_tmp(nxny),stat=istat)
      CALL check_alloc_status(istat, 'wrtradcol:icolp_tmp')
      ALLOCATE (jcolp_tmp(nxny),stat=istat)
      CALL check_alloc_status(istat, 'wrtradcol:jcolp_tmp')
      ALLOCATE (xcolp_tmp(nxny),stat=istat)
      CALL check_alloc_status(istat, 'wrtradcol:xcolp_tmp')
      ALLOCATE (ycolp_tmp(nxny),stat=istat)
      CALL check_alloc_status(istat, 'wrtradcol:ycolp_tmp')
      ALLOCATE (zcolp_tmp(nz,nxny),stat=istat)
      CALL check_alloc_status(istat, 'wrtradcol:zcolp_tmp')
      ALLOCATE (colref_tmp(nz,nxny),stat=istat)
      CALL check_alloc_status(istat, 'wrtradcol:colref_tmp')
      ALLOCATE (colvel_tmp(nz,nxny),stat=istat)
      CALL check_alloc_status(istat, 'wrtradcol:colvel_tmp')
      ALLOCATE (colnyq_tmp(nz,nxny),stat=istat)
      CALL check_alloc_status(istat, 'wrtradcol:colnyq_tmp')
      ALLOCATE (coltim_tmp(nz,nxny),stat=istat)
      CALL check_alloc_status(istat, 'wrtradcol:coltim_tmp')
      IF( dualpdata ) THEN
        ALLOCATE (colrhv_tmp(nz,nxny),stat=istat)
        CALL check_alloc_status(istat, 'wrtradcol:colrhv_tmp')
        ALLOCATE (colzdr_tmp(nz,nxny),stat=istat)
        CALL check_alloc_status(istat, 'wrtradcol:colzdr_tmp')
        ALLOCATE (colkdp_tmp(nz,nxny),stat=istat)
        CALL check_alloc_status(istat, 'wrtradcol:colkdp_tmp')
      END IF
    END IF
  END IF

  CALL ctim2abss(iyr,imon,iday,ihr,imin,isec,itime)

  IF(dualpdata) THEN
    nradvr=9
  ELSE
    nradvr=6
  END IF

  irngmin=nint(refrngmin)
  irngmax=nint(refrngmax)

  IF(myproc == 0) THEN
    WRITE(6,'(a,i9)')   ' irngmin: ',irngmin
    WRITE(6,'(a,i9)')   ' irngmax: ',irngmax
    WRITE(6,'(a,f7.2)') ' refelvmin: ',refelvmin
    WRITE(6,'(a,f7.2)') ' refelvmax: ',refelvmax
  END IF

  IF(myproc == 0 .AND. dmpzero > 0) THEN

    IF (dmpfmt == 1) THEN
      CALL wtradhdrbin(rfname,radid,itime,vcpnum,isource,iradfmt,       &
                       ncoltot,nz,typelev,dualpol,                      &
                       irngmin,irngmax,radlat,radlon,radalt,            &
                       refelvmin,refelvmax,mxradvr,nradvr,iradvr,       &
                       iunit,fileopnd,istatus)
    ELSE
      CALL wtradhdrhdf(rfname,radid,itime,vcpnum,isource,iradfmt,       &
                       nxlg,nylg,nz,typelev,dualpol,                    &
                       irngmin,irngmax,radlat,radlon,radalt,            &
                       refelvmin,refelvmax,mxradvr,nradvr,iradvr,       &
                       sd_id,fileopnd,istatus)
    END IF

  END IF

!
!-----------------------------------------------------------------------
!
!  For each horizontal grid point form a column of remapped
!  data containing the non-missing grid points
!
!-----------------------------------------------------------------------
!
  IF( dmpfmt==1 )THEN
    IF(myproc == 0) THEN
      ALLOCATE(outk1d(nz),stat=istat)
      CALL check_alloc_status(istat, 'wrtradcol:outk1d')
      ALLOCATE(outhgt1d(nz),stat=istat)
      CALL check_alloc_status(istat, 'wrtradcol:outhgt1d')
      ALLOCATE(outref1d(nz),stat=istat)
      CALL check_alloc_status(istat, 'wrtradcol:outref1d')
      ALLOCATE(outvel1d(nz),stat=istat)
      CALL check_alloc_status(istat, 'wrtradcol:outvel1d')
      ALLOCATE(outnyq1d(nz),stat=istat)
      CALL check_alloc_status(istat, 'wrtradcol:outnyq1d')
      ALLOCATE(outtim1d(nz),stat=istat)
      CALL check_alloc_status(istat, 'wrtradcol:outtim1d')
      IF(dualpdata) THEN
        ALLOCATE(outrhv1d(nz),stat=istat)
        CALL check_alloc_status(istat, 'wrtradcol:outrhv1d')
        ALLOCATE(outzdr1d(nz),stat=istat)
        CALL check_alloc_status(istat, 'wrtradcol:outzdr1d')
        ALLOCATE(outkdp1d(nz),stat=istat)
        CALL check_alloc_status(istat, 'wrtradcol:outkdp1d')
      END IF
!
!-----------------------------------------------------------------------
!
!  First write the data calculated within processor zero.
!  For non-MPI case these are all the data.
!
!-----------------------------------------------------------------------
!
      kntcol=0
      DO idat=1,ncolp
        outk1d=0
        outhgt1d=misval
        outref1d=misval
        outvel1d=misval
        outnyq1d=misval
        outtim1d=misval
        IF(dualpdata) THEN
          outrhv1d=misval
          outzdr1d=misval
          outkdp1d=misval
        END IF
        klev=0
        DO k=2,nz-2
          IF((colref(k,idat)>refmin .AND. colref(k,idat)<refmax) .OR.   &
             (colvel(k,idat)>velmin .AND. colvel(k,idat)<velmax))THEN
            klev=klev+1
            outk1d(klev)=k
            outhgt1d(klev)=zcolp(k,idat)
            outref1d(klev)=colref(k,idat)
            outvel1d(klev)=colvel(k,idat)
            outnyq1d(klev)=colnyq(k,idat)
            outtim1d(klev)=coltim(k,idat)
            IF(dualpdata) THEN
              outrhv1d(klev)=colrhv(k,idat)
              outzdr1d(klev)=colzdr(k,idat)
              outkdp1d(klev)=colkdp(k,idat)
            END IF

          END IF
        END DO
!
!-----------------------------------------------------------------------
!
!  If there are data in this column, write them to the file.
!
!-----------------------------------------------------------------------
!
        IF(klev > 0) THEN
          IF(.NOT. fileopnd)                                            &
            CALL wtradhdrbin(rfname,radid,itime,vcpnum,isource,iradfmt, &
                     ncoltot,nz,typelev,dualpol,                        &
                     irngmin,irngmax,radlat,radlon,radalt,              &
                     refelvmin,refelvmax,mxradvr,nradvr,iradvr,         &
                     iunit,fileopnd,istatus)
          kntcol=kntcol+1
          CALL xytoll(1,1,xcolp(idat),ycolp(idat),collat,collon)
          elev=0.5*(zcolp(1,idat)+zcolp(2,idat))
          WRITE(iunit) icolp(idat),jcolp(idat),xcolp(idat),ycolp(idat), &
                       collat,collon,elev,klev
          WRITE(iunit) (outk1d(k),k=1,klev)
          WRITE(iunit) (outhgt1d(k),k=1,klev)
          WRITE(iunit) (outref1d(k),k=1,klev)
          WRITE(iunit) (outvel1d(k),k=1,klev)
          WRITE(iunit) (outnyq1d(k),k=1,klev)
          WRITE(iunit) (outtim1d(k),k=1,klev)
          IF(dualpdata) THEN
            WRITE(iunit) (outrhv1d(k),k=1,klev)
            WRITE(iunit) (outzdr1d(k),k=1,klev)
            WRITE(iunit) (outkdp1d(k),k=1,klev)
          END IF
        END IF
      END DO
!
!-----------------------------------------------------------------------
!
!  Gather data from the other processors, and write them in the
!  same manner.
!
!-----------------------------------------------------------------------
      IF( nprocs > 1) THEN
        WRITE(6,'(a,i4)') ' mp_opt= ',mp_opt
        WRITE(6,'(a,i6)') ' nprocs= ',nprocs
        WRITE(6,'(a,i6,a,i9)') ' kntcol processor: ',myproc,' = ',kntcol

        DO jproc=1,nprocs-1
          kntcolprc=0

          colref_tmp=misval
          colvel_tmp=misval
          IF(dualpdata) THEN
            colrhv_tmp=misval
            colzdr_tmp=misval
            colkdp_tmp=misval
          END IF

!
!  Receive data sent from processor jproc
!
          WRITE(6,'(a,i6,a,i6)') ' myproc=',myproc,' getting data from jproc=',jproc
          tagbase=jproc*tagmult

          CALL mprecvi(ncolj,1,jproc,tagbase,ierr)
          IF(ierr /= 0) THEN
            WRITE(6,'(a,i6,a,i6)') 'Error status: ',ierr, &
                    ' receiving icolp sent from jproc:',jproc
            RETURN
          END IF

          CALL mprecvi(icolp_tmp,nxny,jproc,(tagbase+1),ierr)
          IF(ierr /= 0) THEN
            WRITE(6,'(a,i6,a,i6)') 'Error status: ',ierr, &
                    ' receiving icolp sent from jproc:',jproc
            RETURN
          END IF

          CALL mprecvi(jcolp_tmp,nxny,jproc,(tagbase+2),ierr)
          IF(ierr /= 0) THEN
            WRITE(6,'(a,i6,a,i6)') 'Error status: ',ierr, &
                    ' receiving jcolp sent from jproc:',jproc
            RETURN
          END IF

          CALL mprecvr(xcolp_tmp,nxny,jproc,(tagbase+3),ierr)
          IF(ierr /= 0) THEN
            WRITE(6,'(a,i6,a,i6)') 'Error status: ',ierr, &
                    ' receiving xcolp sent from jproc:',jproc
            RETURN
          END IF

          CALL mprecvr(ycolp_tmp,nxny,jproc,(tagbase+4),ierr)
          IF(ierr /= 0) THEN
            WRITE(6,'(a,i6,a,i6)') 'Error status: ',ierr, &
                    ' receiving ycolp sent from jproc:',jproc
            RETURN
          END IF

          nsize=nz*nxny
          CALL mprecvr(zcolp_tmp,nsize,jproc,(tagbase+5),ierr)
          IF(ierr /= 0) THEN
            WRITE(6,'(a,i6,a,i6)') 'Error status: ',ierr, &
                    ' receiving zcolp sent from jproc:',jproc
            RETURN
          END IF

          CALL mprecvr(colref_tmp,nsize,jproc,(tagbase+6),ierr)
          IF(ierr /= 0) THEN
            WRITE(6,'(a,i6,a,i6)') 'Error status: ',ierr, &
                    ' receiving col_ref sent from jproc:',jproc
            RETURN
          END IF

          CALL mprecvr(colvel_tmp,nsize,jproc,(tagbase+7),ierr)
          IF(ierr /= 0) THEN
            WRITE(6,'(a,i6,a,i6)') 'Error status: ',ierr, &
                    ' receiving col_vel sent from jproc:',jproc
            RETURN
          END IF

          CALL mprecvr(colnyq_tmp,nsize,jproc,(tagbase+8),ierr)
          IF(ierr /= 0) THEN
            WRITE(6,'(a,i6,a,i6)') 'Error status: ',ierr, &
                    ' receiving colnyq sent from jproc:',jproc
            RETURN
          END IF

          CALL mprecvr(coltim_tmp,nsize,jproc,(tagbase+9),ierr)
          IF(ierr /= 0) THEN
            WRITE(6,'(a,i6,a,i6)') 'Error status: ',ierr, &
                    ' receiving coltim sent from jproc:',jproc
            RETURN
          END IF

          IF(dualpdata) THEN
            CALL mprecvr(colrhv_tmp,nsize,jproc,(tagbase+10),ierr)
            IF(ierr /= 0) THEN
              WRITE(6,'(a,i6,a,i6)') 'Error status: ',ierr, &
                      ' receiving colrhv sent from jproc:',jproc
              RETURN
            END IF

            CALL mprecvr(colzdr_tmp,nsize,jproc,(tagbase+11),ierr)
            IF(ierr /= 0) THEN
              WRITE(6,'(a,i6,a,i6)') 'Error status: ',ierr, &
                      ' receiving colzdr sent from jproc:',jproc
              RETURN
            END IF

            CALL mprecvr(colkdp_tmp,nsize,jproc,(tagbase+12),ierr)
            IF(ierr /= 0) THEN
              WRITE(6,'(a,i6,a,i6)') 'Error status: ',ierr, &
                      ' receiving colkdp sent from jproc:',jproc
              RETURN
            END IF
          END IF ! dualpdata

          DO idat=1,ncolj

            outk1d=0
            outhgt1d=misval
            outref1d=misval
            outvel1d=misval
            outnyq1d=misval
            outtim1d=misval
            IF(dualpdata) THEN
              outrhv1d=misval
              outzdr1d=misval
              outkdp1d=misval
            END IF

            klev=0
            DO k=2,nz-2
              IF((colref_tmp(k,idat)>refmin .AND. colref_tmp(k,idat)<refmax) .OR.   &
                 (colvel_tmp(k,idat)>velmin .AND. colvel_tmp(k,idat)<velmax))THEN
                klev=klev+1
                outk1d(klev)=k
                outhgt1d(klev)=zcolp_tmp(k,idat)
                outref1d(klev)=colref_tmp(k,idat)
                outvel1d(klev)=colvel_tmp(k,idat)
                outnyq1d(klev)=colnyq_tmp(k,idat)
                outtim1d(klev)=coltim_tmp(k,idat)
                IF(dualpdata) THEN
                  outrhv1d(klev)=colrhv_tmp(k,idat)
                  outzdr1d(klev)=colzdr_tmp(k,idat)
                  outkdp1d(klev)=colkdp_tmp(k,idat)
                END IF

              END IF
            END DO
!
!-----------------------------------------------------------------------
!
!  If there are data in this column, write them to the file.
!
!-----------------------------------------------------------------------
!
            IF(klev > 0) THEN
              IF(.NOT.fileopnd)                                         &
                CALL wtradhdrbin(rfname,radid,itime,vcpnum,isource,     &
                     iradfmt,ncoltot,nz,typelev,dualpol,                &
                     irngmin,irngmax,radlat,radlon,radalt,              &
                     refelvmin,refelvmax,mxradvr,nradvr,iradvr,         &
                     iunit,fileopnd,istatus)
              kntcolprc=kntcolprc+1
              CALL xytoll(1,1,xcolp_tmp(idat),ycolp_tmp(idat),collat,collon)
              elev=0.5*(zcolp_tmp(1,idat)+zcolp_tmp(2,idat))
              WRITE(iunit) icolp_tmp(idat),jcolp_tmp(idat),             &
                           xcolp_tmp(idat),ycolp_tmp(idat),             &
                           collat,collon,elev,klev
              WRITE(iunit) (outk1d(k),k=1,klev)
              WRITE(iunit) (outhgt1d(k),k=1,klev)
              WRITE(iunit) (outref1d(k),k=1,klev)
              WRITE(iunit) (outvel1d(k),k=1,klev)
              WRITE(iunit) (outtim1d(k),k=1,klev)
              IF(dualpdata) THEN
                WRITE(iunit) (outrhv1d(k),k=1,klev)
                WRITE(iunit) (outzdr1d(k),k=1,klev)
                WRITE(iunit) (outkdp1d(k),k=1,klev)
              END IF

            END IF
          END DO
          WRITE(6,'(a,i9)') ' kntcol processor: ',jproc,' = ',kntcolprc
          kntcol=kntcol+kntcolprc
          WRITE(6,'(a,i9)') ' kntcol sum: ',kntcol
        END DO ! process data from each processor

        DEALLOCATE (icolp_tmp,jcolp_tmp,xcolp_tmp,ycolp_tmp,zcolp_tmp)
        DEALLOCATE (colref_tmp,colvel_tmp,colnyq_tmp,coltim_tmp)

      END IF ! nprocs > 0

      DEALLOCATE(outk1d, outhgt1d, outref1d, outvel1d, outnyq1d, outtim1d)

    ELSE !  not processor zero
!
!-----------------------------------------------------------------------
!
!  Other processors send their data to processor zero.
!
!-----------------------------------------------------------------------
!
      print *, ' myproc = ',myproc,' sending data to iproc=0'
      tagbase=myproc*tagmult

      CALL mpsendi(ncolj,1,0,tagbase,ierr)

      CALL mpsendi(icolp,nxny,0,(tagbase+1),ierr)
      IF(ierr /= 0) THEN
        WRITE(6,'(a,i6,a,i6)') 'Error status: ',ierr, &
                    ' sending icolp to proc0 from:',myproc
        RETURN
      END IF

      CALL mpsendi(jcolp,nxny,0,(tagbase+2),ierr)
      IF(ierr /= 0) THEN
        WRITE(6,'(a,i6,a,i6)') 'Error status: ',ierr, &
                ' sending jcolp to proc0 from:',myproc
        RETURN
      END IF

      CALL mpsendr(xcolp,nxny,0,(tagbase+3),ierr)
      IF(ierr /= 0) THEN
        WRITE(6,'(a,i6,a,i6)') 'Error status: ',ierr, &
                ' sending xcolp to proc0 from:',myproc
        RETURN
      END IF

      CALL mpsendr(ycolp,nxny,0,(tagbase+4),ierr)
      IF(ierr /= 0) THEN
        WRITE(6,'(a,i6,a,i6)') 'Error status: ',ierr, &
                ' sending ycolp to proc0 from:',myproc
        RETURN
      END IF

      nsize=nz*nxny
      CALL mpsendr(zcolp,nsize,0,(tagbase+5),ierr)
      IF(ierr /= 0) THEN
        WRITE(6,'(a,i6,a,i6)') 'Error status: ',ierr, &
                ' sending zcolp to proc0 from:',myproc
        RETURN
      END IF

      CALL mpsendr(colref,nsize,0,(tagbase+6),ierr)
      IF(ierr /= 0) THEN
        WRITE(6,'(a,i6,a,i6)') 'Error status: ',ierr, &
                ' sending col_ref to proc0 from:',myproc
        RETURN
      END IF

      CALL mpsendr(colvel,nsize,0,(tagbase+7),ierr)
      IF(ierr /= 0) THEN
        WRITE(6,'(a,i6,a,i6)') 'Error status: ',ierr, &
                ' sending col_vel to proc0 from:',myproc
        RETURN
      END IF

      CALL mpsendr(colnyq,nsize,0,(tagbase+8),ierr)
      IF(ierr /= 0) THEN
        WRITE(6,'(a,i6,a,i6)') 'Error status: ',ierr, &
                ' sending colnyq to proc0 from:',myproc
        RETURN
      END IF

      CALL mpsendr(coltim,nsize,0,(tagbase+9),ierr)
      IF(ierr /= 0) THEN
        WRITE(6,'(a,i6,a,i6)') 'Error status: ',ierr, &
                ' sending coltim to proc0 from:',myproc
        RETURN
      END IF

      IF(dualpdata) THEN
        CALL mpsendr(colrhv,nsize,0,(tagbase+10),ierr)
        IF(ierr /= 0) THEN
          WRITE(6,'(a,i6,a,i6)') 'Error status: ',ierr, &
                  ' sending colrhv to proc0 from:',myproc
          RETURN
        END IF

        CALL mpsendr(colzdr,nsize,0,(tagbase+11),ierr)
        IF(ierr /= 0) THEN
          WRITE(6,'(a,i6,a,i6)') 'Error status: ',ierr, &
                  ' sending colzdr to proc0 from:',myproc
          RETURN
        END IF

        CALL mpsendr(colkdp,nsize,0,(tagbase+12),ierr)
        IF(ierr /= 0) THEN
          WRITE(6,'(a,i6,a,i6)') 'Error status: ',ierr, &
                  ' sending colkdp to proc0 from:',myproc
          RETURN
        END IF

      END IF

    END IF  ! processor = 0

  ELSE IF (dmpfmt > 1) THEN   ! HDF format
    !
    ! Improve compression by setting Time and Nyquist values to misval
    ! where corresponding velocity and reflectivity values are missing.
    !
    DO idat=1,ncolp
      DO k=1,nz
        IF(colvel(k,idat) < velmin .OR. colvel(k,idat) > velmax) THEN
          colnyq(k,idat)=misval
          IF(colref(k,idat) < refmin .OR. colref(k,idat) > refmax) &
            coltim(k,idat)=misval
        END IF
      END DO
    END DO

    kntcol=0
    numlev = 0
    rmax=-999.
    rmin=999.
    vmax=-999.
    vmin=999.
    xmin =  9.E37
    ymin =  9.E37
    xmax = -9.E37
    ymax = -9.E37
    ALLOCATE(havdat(nxny),stat=istat)
    CALL check_alloc_status(istat, 'wrtradcol:havdat')
    havdat=.FALSE.

    IF(myproc == 0) THEN
      WRITE(6,'(a,2f9.1)') ' check limits, ref: ',refmin,refmax
      WRITE(6,'(a,2f9.1)') ' check limits, vel: ',velmin,velmax
    END IF

    DO idat=1,ncolp
      klev=0
      DO k=2,nz-2
        IF((colref(k,idat)>refmin .AND. colref(k,idat)<refmax) .OR.   &
           (colvel(k,idat)>velmin .AND. colvel(k,idat)<velmax))THEN
          rmin=min(rmin,colref(k,idat))
          rmax=max(rmax,colref(k,idat))
          vmin=min(vmin,colvel(k,idat))
          vmax=max(vmax,colvel(k,idat))
          klev=klev+1
        END IF
      END DO

      IF(lvldbg > 10 .AND. mod(idat,50) == 0) THEN
        WRITE(6,'(a,2i6)') ' idat,klev: ',idat,klev
        WRITE(6,'(a,2f9.1)') ' reflectvity rmin,rmax: ',rmin,rmax
        WRITE(6,'(a,2f9.1)') ' velocity vmin,vmax: ',vmin,vmax
        WRITE(6,'(a,3f9.1)') ' colvel: ',colvel(5,idat),colvel(10,idat),colvel(20,idat)
      END IF

      IF(klev > 0) THEN
        kntcol=kntcol+1
        havdat(idat)=.TRUE.
        xmin = MIN(xmin,xcolp(idat))
        xmax = MAX(xmax,xcolp(idat))
        ymin = MIN(ymin,ycolp(idat))
        ymax = MAX(ymax,ycolp(idat))
        numlev=max(numlev,klev)
      END IF

    END DO

    WRITE(6,'(a,i6,a,i12,a,i6)') '  Processor:',myproc,'  ncolp:',ncolp,'  numlev:',numlev

    CALL mptotali(kntcol)

    IF ( kntcol > 0 .OR. dmpzero == 1 ) THEN

      !WRITE(6,'(a,i6,a,i6)') '  Processor:',myproc,' numlev = ',numlev

      CALL mpmaxi(numlev)

      !WRITE(6,'(a,i6,a,i6)') '  Processor:',myproc,' new numlev:',numlev

      !
      ! Allocate HDF output arrays
      !
      ALLOCATE(outi(ncolp), STAT = istat)
      CALL check_alloc_status(istat, 'wrtradcol:outi')
      ALLOCATE(outj(ncolp), STAT = istat)
      CALL check_alloc_status(istat, 'wrtradcol:outj')
      ALLOCATE(outk(numlev,ncolp), STAT = istat)
      CALL check_alloc_status(istat, 'wrtradcol:outk')
      ALLOCATE(outnumlev(ncolp), STAT = istat)
      CALL check_alloc_status(istat, 'wrtradcol:outnumlev')

      ALLOCATE(outx(ncolp), STAT = istat)
      CALL check_alloc_status(istat, 'wrtradcol:outx')
      ALLOCATE(outy(ncolp), STAT = istat)
      CALL check_alloc_status(istat, 'wrtradcol:outy')
      ALLOCATE(outlat(ncolp), STAT = istat)
      CALL check_alloc_status(istat, 'wrtradcol:outlat')
      ALLOCATE(outlon(ncolp), STAT = istat)
      CALL check_alloc_status(istat, 'wrtradcol:outlon')
      ALLOCATE(outsfc(ncolp), STAT = istat)
      CALL check_alloc_status(istat, 'wrtradcol:outsfc')

      outi = 0
      outj = 0
      outk = 0
      outnumlev = 0

      outx = misval
      outy = misval
      outlat = misval
      outlon = misval
      outsfc = misval

      ALLOCATE(outhgt(numlev,ncolp), STAT = istat)
      CALL check_alloc_status(istat, 'wrtradcol:outsfc')
      ALLOCATE(outref(numlev,ncolp), STAT = istat)
      CALL check_alloc_status(istat, 'wrtradcol:outref')
      ALLOCATE(outvel(numlev,ncolp), STAT = istat)
      CALL check_alloc_status(istat, 'wrtradcol:outvel')
      ALLOCATE(outnyq(numlev,ncolp), STAT = istat)
      CALL check_alloc_status(istat, 'wrtradcol:outnyq')
      ALLOCATE(outtim(numlev,ncolp), STAT = istat)
      CALL check_alloc_status(istat, 'wrtradcol:outtim')


      outhgt = misval
      outref = misval
      outvel = misval
      outnyq = misval
      outtim = misval

      IF(dualpdata) THEN
        ALLOCATE(outrhv(numlev,ncolp), STAT = istat)
        CALL check_alloc_status(istat, 'wrtradcol:outrhv')
        ALLOCATE(outzdr(numlev,ncolp), STAT = istat)
        CALL check_alloc_status(istat, 'wrtradcol:outzdr')
        ALLOCATE(outkdp(numlev,ncolp), STAT = istat)
        CALL check_alloc_status(istat, 'wrtradcol:outkdp')

        outrhv = misval
        outzdr = misval
        outkdp = misval
      END IF

!
! Transfer to HDF output arrays
!
      jdat=0
      DO idat=1,ncolp
        IF(havdat(idat)) THEN
          jdat=jdat+1
          outi(jdat)=icolp(idat)
          outj(jdat)=jcolp(idat)
          CALL xytoll(1,1,xcolp(idat),ycolp(idat),collat,collon)
          outlat(jdat)=collat
          outlon(jdat)=collon
          klev=0
          DO k=2,nz-2
            IF((colref(k,idat)>refmin .AND. colref(k,idat)<refmax) .OR.   &
               (colvel(k,idat)>velmin .AND. colvel(k,idat)<velmax))THEN
              klev=klev+1
              outk(klev,jdat)=k
              outhgt(klev,jdat)=zcolp(k,idat)
              outref(klev,jdat)=colref(k,idat)
              outvel(klev,jdat)=colvel(k,idat)
              outnyq(klev,jdat)=colnyq(k,idat)
              outtim(klev,jdat)=coltim(k,idat)
              IF(dualpdata) THEN
                outrhv(klev,jdat)=colrhv(k,idat)
                outzdr(klev,jdat)=colzdr(k,idat)
                outkdp(klev,jdat)=colkdp(k,idat)
              END IF

            END IF
          END DO
          outnumlev(jdat)=klev
        END IF  ! data here?
      END DO

      IF (myproc == 0)                                                 &
        WRITE(6,'(a,i12)') ' Total column count:',ncoltot

      IF( ncoltot > 0 .OR. dmpzero == 1) THEN
        IF (myproc == 0) THEN
          WRITE(6,'(a)') ' Writing HDF variables ...'
          IF(.NOT. fileopnd) THEN
            CALL wtradhdrhdf(rfname,radid,itime,vcpnum,isource,iradfmt,&
                   nxlg,nylg,nz,typelev,dualpol,                       &
                   irngmin,irngmax,radlat,radlon,radalt,               &
                   refelvmin,refelvmax,mxradvr,nradvr,iradvr,          &
                   sd_id,fileopnd,istatus)
            print *, ' back from wtradhdrhdf '
          END IF
        END IF

        CALL mpmax0(xmax,xmin)
        CALL mpmax0(ymax,ymin)

        IF (myproc == 0) THEN
          CALL hdfwrtr(sd_id, 'xmin',   xmin,    istat)
          CALL hdfwrtr(sd_id, 'xmax',   xmax,    istat)
          CALL hdfwrtr(sd_id, 'ymin',   xmin,    istat)
          CALL hdfwrtr(sd_id, 'ymax',   ymax,    istat)

          CALL hdfwrti(sd_id, 'numradcol',ncoltot,istat)
          CALL hdfwrti(sd_id, 'nummaxelv',numlev,istat)
          CALL hdfwrti(sd_id, 'dualpol',dualpol,istat)

          ALLOCATE (itmp1dlg(ncoltot),stat=istat)
          CALL check_alloc_status(istat, 'wrtradcol:itmp1dlg')
          ALLOCATE (itmp2dlg(numlev,ncoltot),stat=istat)
          CALL check_alloc_status(istat, 'wrtradcol:itmp2dlg')
          ALLOCATE (rtmp1dlg(ncoltot),stat=istat)
          CALL check_alloc_status(istat, 'wrtradcol:rtmp1dlg')
          ALLOCATE (rtmp2dlg(numlev,ncoltot),stat=istat)
          CALL check_alloc_status(istat, 'wrtradcol:rtmp2dlg')
        ELSE
          ALLOCATE (itmp1dlg(1),stat=istat)
          CALL check_alloc_status(istat, 'wrtradcol:itmp1dlg')
          ALLOCATE (itmp2dlg(1,1),stat=istat)
          CALL check_alloc_status(istat, 'wrtradcol:itmp2dlg')
          ALLOCATE (rtmp1dlg(1),stat=istat)
          CALL check_alloc_status(istat, 'wrtradcol:rtmp1dlg')
          ALLOCATE (rtmp2dlg(1,1),stat=istat)
          CALL check_alloc_status(istat, 'wrtradcol:rtmp2dlg')
        END IF

        ALLOCATE(disp1d(nprocs), STAT = istat)
        CALL check_alloc_status(istat, 'wrtradcol:disp1d')
        ALLOCATE(disp2d(nprocs), STAT = istat)
        CALL check_alloc_status(istat, 'wrtradcol:disp2d')

        ALLOCATE(lsize1d(nprocs), STAT = istat)
        CALL check_alloc_status(istat, 'wrtradcol:lsize1d')
        ALLOCATE(lsize2d(nprocs), STAT = istat)
        CALL check_alloc_status(istat, 'wrtradcol:lsize2d')

        CALL mpgatheri(ncolp,1,lsize1d,nprocs,istatus)

        IF (myproc == 0) THEN
          disp1d(1) = 0
          DO iloc = 2,nprocs
            disp1d(iloc) = disp1d(iloc-1)+lsize1d(iloc-1)
          END DO

          DO iloc = 1,nprocs
            lsize2d(iloc) = numlev*lsize1d(iloc)
            disp2d(iloc)  = numlev*disp1d(iloc)
            !WRITE(6,'(a,5I10)') 'iloc = ',iloc, lsize1d(iloc),disp1d(iloc),lsize2d(iloc),disp2d(iloc)
            !CALL flush(6)
          END DO
        ELSE
          disp1d(:) = 999    ! any number > 0
          disp2d(:) = 999
          lsize1d(:) = 999
          lsize2d(:) = 999
        END IF

        CALL mpigatherr(outlat,ncolp,rtmp1dlg,ncoltot,                  &
                        lsize1d,disp1d,0,nprocs,istatus)

        IF(myproc == 0) CALL hdfwrt1d(rtmp1dlg,ncoltot,sd_id,           &
           'radcollat','Latitude array for the radar columns','degrees N')

        CALL mpigatherr(outlon,ncolp,rtmp1dlg,ncoltot,                  &
                        lsize1d,disp1d,0,nprocs,istatus)

        IF(myproc == 0) CALL hdfwrt1d(rtmp1dlg,ncoltot,sd_id,           &
          'radcollon','Longitude array for the radar columns','degrees E')

        CALL mpigatheri(outnumlev,ncolp,itmp1dlg,ncoltot,               &
                        lsize1d,disp1d,0,nprocs,istatus)

        IF(myproc == 0) CALL hdfwrt1di(itmp1dlg,ncoltot,sd_id,          &
           'numelev','Number of vertical levels for each column','index')

        CALL mpigatheri(outi,ncolp,itmp1dlg,ncoltot,                    &
                        lsize1d,disp1d,0,nprocs,istatus)

        IF(myproc == 0) CALL hdfwrt1di(itmp1dlg,ncoltot,sd_id,          &
           'radcoli','i-index in the ARPS grid','index')

        CALL mpigatheri(outj,ncolp,itmp1dlg,ncoltot,                    &
                        lsize1d,disp1d,0,nprocs,istatus)

        IF(myproc == 0) CALL hdfwrt1di(itmp1dlg,ncoltot,sd_id,          &
           'radcolj','j-index in the ARPS grid','index')

        ncolp2d   = numlev*ncolp
        ncoltot2d = numlev*ncoltot

        CALL mpigatheri(outk,ncolp2d,itmp2dlg,ncoltot2d,                &
                        lsize2d,disp2d,0,nprocs,istatus)

        IF(myproc == 0) CALL hdfwrt2di(itmp2dlg,numlev,ncoltot,sd_id,0,0, &
            'radcolk','k-index in the ARPS grid','index')

        IF (myproc == 0 .AND. hdf4cmpr > 3) THEN
          ALLOCATE (itmp(numlev,ncoltot),stat=istat)
          CALL check_alloc_status(istat, 'wrtradcol:itmp')
        ELSE
          ALLOCATE (itmp(1,1),stat=istat)
          CALL check_alloc_status(istat, 'wrtradcol:itmp')
        END IF

        CALL mpigatherr(outhgt,ncolp2d,rtmp2dlg,ncoltot2d,              &
                        lsize2d,disp2d,0,nprocs,istatus)

        IF(myproc == 0) CALL hdfwrt2d(rtmp2dlg,numlev,ncoltot,sd_id,    &
                 0,hdf4cmpr,'radcolhgt','height','m',itmp)

        CALL mpigatherr(outref,ncolp2d,rtmp2dlg,ncoltot2d,              &
                        lsize2d,disp2d,0,nprocs,istatus)

        IF(myproc == 0) CALL hdfwrt2d(rtmp2dlg,numlev,ncoltot,sd_id,    &
                 0,hdf4cmpr,'radcolref','reflectivity','dBZ',itmp)

        CALL mpigatherr(outvel,ncolp2d,rtmp2dlg,ncoltot2d,              &
                        lsize2d,disp2d,0,nprocs,istatus)

        IF(myproc == 0) CALL hdfwrt2d(rtmp2dlg,numlev,ncoltot,sd_id,    &
                 0,hdf4cmpr,'radcolvel','radial velocity','m/s',itmp)

        CALL mpigatherr(outnyq,ncolp2d,rtmp2dlg,ncoltot2d,              &
                        lsize2d,disp2d,0,nprocs,istatus)

        IF(myproc == 0) CALL hdfwrt2d(rtmp2dlg,numlev,ncoltot,sd_id,    &
                 0,hdf4cmpr,'radcolnyq','nyquist velocity','m/s',itmp)

        CALL mpigatherr(outtim,ncolp2d,rtmp2dlg,ncoltot2d,              &
                        lsize2d,disp2d,0,nprocs,istatus)

        IF(myproc == 0) CALL hdfwrt2d(rtmp2dlg,numlev,ncoltot,sd_id,    &
                 0,hdf4cmpr,'radcoltim','relative time','seconds',itmp)

        IF(dualpdata) THEN
          CALL mpigatherr(outrhv,ncolp2d,rtmp2dlg,ncoltot2d,            &
                        lsize2d,disp2d,0,nprocs,istatus)

          IF(myproc == 0) CALL hdfwrt2d(rtmp2dlg,numlev,ncoltot,sd_id,  &
                 0,hdf4cmpr,'radcolrhv','Rho-HV',' ',itmp)

          CALL mpigatherr(outzdr,ncolp2d,rtmp2dlg,ncoltot2d,              &
                        lsize2d,disp2d,0,nprocs,istatus)

          IF(myproc == 0) CALL hdfwrt2d(rtmp2dlg,numlev,ncoltot,sd_id,    &
                 0,hdf4cmpr,'radcolzdr','Zdr','dBZ',itmp)

          CALL mpigatherr(outkdp,ncolp2d,rtmp2dlg,ncoltot2d,              &
                        lsize2d,disp2d,0,nprocs,istatus)

          IF(myproc == 0) CALL hdfwrt2d(rtmp2dlg,numlev,ncoltot,sd_id,    &
                   0,hdf4cmpr,'radcolkdp','Kdp','degrees/km',itmp)
        END IF ! dualpdata

        !IF (myproc == 0) THEN
        !  WRITE(6,'(a)') ' MPI HDF write complete. Cleaning-up'
        !  !CALL flush(6)
        !END IF

        DEALLOCATE(itmp1dlg,itmp2dlg)
        DEALLOCATE(rtmp1dlg,rtmp2dlg)
        DEALLOCATE(havdat)

      ELSE   ! zero columns and no zero-write

        IF (myproc == 0) THEN
          WRITE(6,'(a,i6)') ' No columns to write.  Exiting.'
          CALL flush(6)
        END IF
        RETURN

      END IF ! valid columns exist

    END IF  ! nprocs

  END IF  ! output format
!
  IF(fileopnd)THEN
    IF(dmpfmt == 1) THEN
      CLOSE(iunit)
      CALL retunit(iunit)
    ELSE
      CALL hdfclose(sd_id,istat)
      IF (istat /= 0) THEN
        WRITE (6,*) 'WRTRADCOL HDF: ERROR on closing file ',trim(rfname), &
                ' (status',istat,')'
    END IF
  END IF

  IF ( dmpfmt > 1 .AND. hdf4cmpr > 3) THEN
    DEALLOCATE (itmp)
  END IF

  END IF
!
!-----------------------------------------------------------------------
!
!  Report on what data were written
!
!-----------------------------------------------------------------------
!
  IF( myproc == 0 ) THEN
    WRITE(6,'(//a,i4.4,i2.2,i2.2,a1,i2.2,a1,i2.2)')                     &
                    ' Output statistics for time ',                     &
                      iyr,imon,iday,' ',ihr,':',imin
    IF( kntcol > 0 .OR. dmpzero == 1) THEN
      WRITE(6,'(a,i6,a,/a,i6,a//)')                                     &
           ' There were ',ncoltot,' columns written ',                  &
           ' of a total ',((nxlg-3)*(nylg-3)),' possible.'
    ELSE
      WRITE(6,'(a,i6,a,/a,//)')                                         &
             ' There were ',kntcol,' columns with non-missing data',    &
           ' no file was written.'
    END IF
  END IF
!
  RETURN
END SUBROUTINE wrtradcol


SUBROUTINE refract(nx,ny,nz,ptprt,pprt,qv,ptbar,pbar,rfrct)
  IMPLICIT NONE
!
  INTEGER, INTENT(IN)  :: nx,ny,nz
  REAL,    INTENT(IN)  :: ptprt(nx,ny,nz)
  REAL,    INTENT(IN)  :: pprt(nx,ny,nz)
  REAL,    INTENT(IN)  :: qv(nx,ny,nz)
  REAL,    INTENT(IN)  :: ptbar(nx,ny,nz)
  REAL,    INTENT(IN)  :: pbar(nx,ny,nz)
  REAL,    INTENT(OUT) :: rfrct(nx,ny,nz)
!
! Constants from Doviak and Zrnic', 1st Ed, Eq 2.18.
! After Bean and Dutton, 1966.
!
  REAL, PARAMETER :: cd1=0.776   ! K/Pa
  REAL, PARAMETER :: cw1=0.716   ! K/Pa
  REAL, PARAMETER :: cw2=3700.   ! K*K/Pa
!
! Misc local variables
!
  INTEGER :: i,j,k
  REAL    :: pr,tk,pw,pd,tkinv,refr
!
  INCLUDE 'phycst.inc'
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  rfrct =0.
!
! Calculate refractivity index from ARPS pressure,
! potential temp and specific humidity.
!
! pr=total pressure, Pa
! tk=temperature, K
! pw=partial pressure of water vapor, Pa
! pd=partial pressure of dry air, Pa
!
  DO k=2,nz-1
    DO j=1,ny-1
      DO i=1,nx-1
        pr=pbar(i,j,k)+pprt(i,j,k)
        tk=(ptbar(i,j,k)+ptprt(i,j,k))*((pr/p0)**rddcp)
        pw=pr*(qv(i,j,k)/((0.378*qv(i,j,k))+0.622))
        pd=pr-pw
        tkinv=1./tk
        rfrct(i,j,k)=(cd1*pd*tkinv)+(cw1*pw*tkinv)+(cw2*pw*tkinv*tkinv)
      END DO
    END DO
  END DO
  RETURN
END SUBROUTINE refract
!
SUBROUTINE refract2(nx,ny,nz,pt,p,qv,rfrct)
  IMPLICIT NONE
!
  INTEGER, INTENT(IN)  :: nx,ny,nz
  REAL,    INTENT(IN)  :: pt(nx,ny,nz)
  REAL,    INTENT(IN)  :: p(nx,ny,nz)
  REAL,    INTENT(IN)  :: qv(nx,ny,nz)
  REAL,    INTENT(OUT) :: rfrct(nx,ny,nz)
!
! Constants from Doviak and Zrnic', 1st Ed, Eq 2.18.
! After Bean and Dutton, 1966.
!
  REAL, PARAMETER :: cd1=0.776   ! K/Pa
  REAL, PARAMETER :: cw1=0.716   ! K/Pa
  REAL, PARAMETER :: cw2=3700.   ! K*K/Pa
!
! Misc local variables
!
  INTEGER :: i,j,k
  REAL    :: tk,pw,pd,tkinv,refr
!
  INCLUDE 'phycst.inc'

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!
  rfrct =0.
!
! Calculate refractivity index from ARPS pressure,
! potential temp and specific humidity.
!
! pr=total pressure, Pa
! tk=temperature, K
! pw=partial pressure of water vapor, Pa
! pd=partial pressure of dry air, Pa
!
!  print *, ' pressure: ',(0.01*p(5,5,5))
!  print *, ' potential temp: ',pt(5,5,5)
!  print *, ' qv:',(1000.*qv(5,5,5))
  DO k=2,nz-1
    DO j=1,ny-1
      DO i=1,nx-1
        tk=pt(i,j,k)*((p(i,j,k)/p0)**rddcp)
        pw=p(i,j,k)*(qv(i,j,k)/((0.378*qv(i,j,k))+0.622))
        pd=p(i,j,k)-pw
        tkinv=1./tk
        rfrct(i,j,k)=(cd1*pd*tkinv)+(cw1*pw*tkinv)+(cw2*pw*tkinv*tkinv)
      END DO
    END DO
  END DO
  RETURN
END SUBROUTINE refract2
!
!##################################################################
!##################################################################
!######                                                      ######
!######                  SUBROUTINE VADVOL                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE vadvol(maxgate,maxazim,maxelev,                              &
                  radalt,velchek,rngvad,rngwid,minknt,                  &
                  kntvgat,kntvazm,kntvelv,                              &
                  rngvvol,azmvvol,elvvvol,velvol,                       &
                  vadhgt,vaddir,vadspd)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Compute vad from volume of radar data from 2-D tilts.
!
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster, CAPS
!  November, 2003
!
!  MODIFICATION HISTORY:
!  Kun Zhao, added RMS check
!  Keith Brewster (3/18/2011), added radar elevation to height and
!                              cleaned-up some code
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    maxgate   Maximum gates in a radial
!    maxazim   Maximum radials per tilt
!    maxelev   Maximum number of tilts
!
!    ngate     Number of gates in each radial
!    nazim     Number of radials in each tilt
!
!    radalt    Height of radar above sea level
!    velchek   Threshold value to determine good vs. flagged data
!
!    kntvgat   Number of gates in each radial (3-D)
!    kntvazm   Number of radials in each tilt (3-D)
!    kntelev   Number of elevation angles (tilts)
!
!    kntvazm   Number of radials in each tilt (3-D)
!    kntvelv   Number of elevation angles (tilts)
!
!    rngvvol    Range to gate
!    azmvvol    Azimuth angle
!    elvvvol    Elevation angle
!    velvol     Radial velocity volume.
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

  INTEGER, INTENT(IN)  :: maxgate
  INTEGER, INTENT(IN)  :: maxazim
  INTEGER, INTENT(IN)  :: maxelev

  INTEGER, INTENT(IN) :: kntvgat(maxazim,maxelev)
  INTEGER, INTENT(IN) :: kntvazm(maxelev)
  INTEGER, INTENT(IN) :: kntvelv

  REAL, INTENT(IN)    :: radalt
  REAL, INTENT(IN)    :: velchek
  REAL, INTENT(IN)    :: rngvad
  REAL, INTENT(IN)    :: rngwid
  INTEGER, INTENT(IN) :: minknt

  REAL, INTENT(IN)    :: rngvvol(maxgate,maxelev)
  REAL, INTENT(IN)    :: azmvvol(maxazim,maxelev)
  REAL, INTENT(IN)    :: elvvvol(maxazim,maxelev)
  REAL, INTENT(IN)    :: velvol(maxgate,maxazim,maxelev)

  REAL, INTENT(OUT)   :: vadhgt(maxelev)
  REAL, INTENT(OUT)   :: vaddir(maxelev)
  REAL, INTENT(OUT)   :: vadspd(maxelev)
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  REAL, PARAMETER :: misval  = -9999.
  REAL, PARAMETER :: mischek = -9990.
  REAL, PARAMETER :: rmslim = 6.0
  INTEGER :: igate,jazim,kelev
  INTEGER :: irgate,nwidth,igbgn,igend,igendj,kntpt,kntpt2
  REAL :: deg2rad,rad2deg,azmrad,sinaz,cosaz,sin2az,cos2az,delrng,vr
  REAL :: sum_q0r,sum_q5r,sum_q5i,sum_q4r,sum_q4i,sum_q3r,sum_q3i
  REAL :: flknt,twon,sum_elv,elvavg,height,sfcrng
  COMPLEX :: q0,q1,q2,q3,q4,q5,ccj_q4,int_coeff,qq_int,qq
  REAL :: cf1,cf2,cf3,dd,ff,delv,sum2,rms
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  deg2rad=acos(-1.)/180.
  rad2deg=1./deg2rad

  vadhgt=misval
  vaddir=misval
  vadspd=misval

  DO kelev = 1, kntvelv
    delrng=rngvvol(2,kelev)-rngvvol(1,kelev)
    nwidth=MIN(NINT((0.5*rngwid)/delrng),1)
    irgate=NINT((rngvad-rngvvol(1,kelev))/delrng)
    igbgn=MAX(1,(irgate-nwidth))
    igend=(irgate+nwidth)
    kntpt=0
    sum_elv=0.
    sum_q0r=0.
    sum_q5r=0.
    sum_q5i=0.
    sum_q4r=0.
    sum_q4i=0.
    sum_q3r=0.
    sum_q3i=0.
!
!   Form sums going around in azim
!
    DO jazim = 1, kntvazm(kelev)
      azmrad=deg2rad*azmvvol(jazim,kelev)
      sinaz=sin(azmrad)
      cosaz=cos(azmrad)
      sin2az=sin(2.0*azmrad)
      cos2az=cos(2.0*azmrad)
      igendj=MIN(igend,kntvgat(jazim,kelev))
      DO igate=igbgn,igendj
        IF(velvol(igate,jazim,kelev) > velchek) THEN
          vr=velvol(igate,jazim,kelev)
          kntpt=kntpt+1
          sum_elv=sum_elv+elvvvol(jazim,kelev)
          sum_q0r=sum_q0r+vr
          sum_q5r=sum_q5r+cos2az
          sum_q5i=sum_q5i+sin2az
          sum_q4r=sum_q4r+cosaz
          sum_q4i=sum_q4i+sinaz
          sum_q3r=sum_q3r+vr*cosaz
          sum_q3i=sum_q3i+vr*sinaz
        END IF
      END DO
    END DO
!
    cf1=misval
    cf2=misval
    cf3=misval
!
    IF(kntpt > minknt) THEN
      flknt=float(kntpt)
      twon=float(2*kntpt)
      elvavg=sum_elv/flknt
      q0=CMPLX(sum_q0r/flknt)
      q5=CMPLX((sum_q5r/twon),-(sum_q5i/twon))
      q4=CMPLX((sum_q4r/twon),(sum_q4i/twon))
      q3=CMPLX((sum_q3r/flknt),-(sum_q3i/flknt))
!
!   Complex conjugate of q4
!
      ccj_q4=CONJG(q4)
      qq=q4-(1./(4.*ccj_q4))
!
      IF( qq /= 0.) THEN

        q2=(ccj_q4-(q5/(2.*ccj_q4)))/qq
        q1=(q0-(q3/(2.*ccj_q4)))/qq
        qq_int=1.-(CABS(q2)*CABS(q2))

        IF( qq_int /= 0.) THEN
          int_coeff=(q1-q2*CONJG(q1))/qq_int
          cf3=IMAG(int_coeff)
          cf2=REAL(int_coeff)
          cf1=REAL(q0)-(2.*REAL(int_coeff*q4))
        END IF

      END IF

    END IF
!
!   RMS check
!
    rms = 0.
    sum2 = 0.
    kntpt2 = 0
    IF( cf1 > mischek ) THEN
     DO jazim = 1, kntvazm(kelev)
       dd=180.-rad2deg*(atan2(cf3,cf2))
       IF(dd > 360.) dd = dd-360.
       dd = deg2rad*dd
       azmrad=deg2rad*azmvvol(jazim,kelev)
       sinaz=sin(azmrad)
       cosaz=cos(azmrad)
       sin2az=sin(2.0*azmrad)
       cos2az=cos(2.0*azmrad)
       igendj=MIN(igend,kntvgat(jazim,kelev))
       DO igate=igbgn,igendj
         IF(velvol(igate,jazim,kelev) > velchek) THEN
           kntpt2 = kntpt2 + 1
           vr =  velvol(igate,jazim,kelev)
           delv = (-cos(azmrad-dd)*sqrt(cf2*cf2+cf3*cf3)+cf1-vr)
           sum2 = sum2 + delv*delv
         END IF
       END DO
     END DO
     IF( kntpt2 > 0 ) THEN
       rms = SQRT(sum2 / float(kntpt2))
       IF ( rms > rmslim ) cf1 = misval
     END IF

    END IF

    IF( cf1 > mischek ) THEN
      ff=sqrt(cf2*cf2 + cf3*cf3)/cos(deg2rad*elvavg)
      dd=180.-rad2deg*(atan2(cf3,cf2))
      IF(dd > 360.) dd = dd-360.
      CALL beamhgt(elvavg,rngvad,height,sfcrng)
      print *, ' elv:',elvavg,' hgt:',height,' dd:',dd,' ff:',ff
      vadhgt(kelev)=height+radalt
      vaddir(kelev)=dd
      vadspd(kelev)=ff
    END IF

  END DO

END SUBROUTINE vadvol
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######              SUBROUTINE MKVADFNM                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE mkvadfnm(dir,radar,iyr,imo,ida,ihr,imin,isec,fname)
!
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
!  Input arguments
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=256), INTENT(IN) :: dir
  CHARACTER (LEN=4), INTENT(IN):: radar
  INTEGER, INTENT(IN) :: iyr
  INTEGER, INTENT(IN) :: imo
  INTEGER, INTENT(IN) :: ida
  INTEGER, INTENT(IN) :: ihr
  INTEGER, INTENT(IN) :: imin
  INTEGER, INTENT(IN) :: isec
!
!-----------------------------------------------------------------------
!
!  Output arguments
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=256), INTENT(OUT) :: fname
!
  INTEGER :: rdtime,jyr,jmo,jda,jhr,jmin,jsec
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!-----------------------------------------------------------------------
!
!  Round to nearest minute
!
!-----------------------------------------------------------------------
!
  CALL ctim2abss( iyr,imo,ida,ihr,imin,isec, rdtime )
  IF(isec >= 30) rdtime=rdtime+(60-isec)

  CALL abss2ctim( rdtime, jyr, jmo, jda, jhr, jmin, jsec )
  jyr = MOD(jyr,100)
  WRITE(fname,'(a,a4,a1,3(i2.2),a1,2(i2.2),a4)')                        &
        TRIM(dir),radar,'.',                                            &
        jyr,jmo,jda,'.',jhr,jmin,'.vad'
  RETURN
END SUBROUTINE mkvadfnm
!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINE WTVADPRF                  ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE wtvadprf(maxelev,vfname,radar,vadsrc,latrad,lonrad,elvrad,  &
                    vadhgt,vaddir,vadspd)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Write VAD wind profile to file to be read by ADAS.
!  Format of wind profiler data is used.
!
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster, CAPS
!  November, 2003
!
!  MODIFICATION HISTORY:
!  08/22/2011  (Keith Brewster)
!  Added vadsrc to argument list and station info output line
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: maxelev
  CHARACTER (LEN=256), INTENT(IN) :: vfname
  CHARACTER (LEN=5), INTENT(IN) :: radar
  CHARACTER(LEN=8), INTENT(IN) :: vadsrc
  REAL, INTENT(IN) :: latrad
  REAL, INTENT(IN) :: lonrad
  REAL, INTENT(IN) :: elvrad
  REAL, INTENT(IN) :: vadhgt(maxelev)
  REAL, INTENT(IN) :: vaddir(maxelev)
  REAL, INTENT(IN) :: vadspd(maxelev)
!
! Misc local variables
!
  INTEGER :: kelev
  INTEGER :: kntvalid,numsta,iunit

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!
! Determine if there are any valid data in this VAD profile
!
  numsta=88888
  kntvalid=0
  DO kelev=1,maxelev
    IF(vaddir(kelev) <= 360. .AND. vaddir(kelev) >= 0.)  &
       kntvalid=kntvalid+1
  END DO

  WRITE(6,'(i6,a)') kntvalid,' valid heights in VAD profile'
  IF ( kntvalid > 0 ) THEN
    WRITE(6,'(a,a,a)') ' Opening ',TRIM(vfname),' for writing'
    CALL GETUNIT(iunit)
    OPEN(iunit,file=TRIM(vfname),status='unknown')
    WRITE(iunit,'(i12,i12,f11.4,f15.4,f15.0,6x,a4,1x,a8)')        &
        numsta,kntvalid,latrad,lonrad,elvrad,radar(1:4),vadsrc
    DO kelev=1,maxelev
      IF(vaddir(kelev) <= 360. .AND. vaddir(kelev) >= 0.) THEN
        WRITE(iunit,'(f10.0,f10.0,f10.1)') &
              vadhgt(kelev),vaddir(kelev),vadspd(kelev)
      END IF
    END DO
    CLOSE(iunit)
    CALL RETUNIT(iunit)
  ELSE
    WRITE(6,'(//a//)') ' No valid VAD data.  Skipping VAD write.'
  END IF
  RETURN
END SUBROUTINE wtvadprf
!
!##################################################################
!##################################################################
!######                                                      ######
!######             SUBROUTINE GET_INFILNAME                 ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################

SUBROUTINE get_infilname(iniopt,hinfmt,grdbasfn,hisfile)
  IMPLICIT NONE
  INTEGER :: iniopt
  INTEGER :: hinfmt
  CHARACTER(LEN=256) :: grdbasfn
  CHARACTER(LEN=256) :: hisfile

  INCLUDE 'globcst.inc'

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  iniopt=initopt
  hinfmt=inifmt
  grdbasfn=inigbf
  hisfile=inifile
  RETURN
END SUBROUTINE get_infilname
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE READUVT                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################

SUBROUTINE readuvt(nx,ny,nz,filenamein,filefmt,nprocx_in,nprocy_in,     &
                   time,u,v,p,pt,qv,ireturn)

!-----------------------------------------------------------------------
!  PURPOSE:
!  Read in ARPS wind data in the NCSA HDF4 format.
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  2000/04/15
!
!  MODIFICATION HISTORY:
!  Keith Brewster, CAPS
!  2005/02/06 Created hdfreaduvw, streamlined from hdfreadwind
!  2008/02/26 Created hdfreaduv
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!    filename  Character variable naming the input HDF file

!-----------------------------------------------------------------------
!  Variable Declarations.
!-----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: nx,ny,nz

  CHARACTER (LEN=256), INTENT(IN) :: filenamein
  INTEGER, INTENT(IN) :: filefmt
  INTEGER, INTENT(IN) :: nprocx_in, nprocy_in

  REAL :: time
  REAL :: u(nx,ny,nz)
  REAL :: v(nx,ny,nz)
  REAL :: p(nx,ny,nz)
  REAL :: pt(nx,ny,nz)
  REAL :: qv(nx,ny,nz)

  INTEGER, INTENT(OUT) :: ireturn

  INTEGER (KIND=selected_int_kind(4)), ALLOCATABLE :: itmp(:,:,:) ! Temporary array
  REAL, ALLOCATABLE :: hmin(:) ! Temporary array
  REAL, ALLOCATABLE :: hmax(:) ! Temporary array

  REAL, ALLOCATABLE :: invar3d(:,:,:)

  CHARACTER (LEN=10) :: tmunit

!-----------------------------------------------------------------------
!  Misc. local variables
!-----------------------------------------------------------------------

  INTEGER :: lchanl, nchanl
  PARAMETER (lchanl=6)      ! Channel number for formatted printing.

  INTEGER :: i,j,k,is
  INTEGER :: nxin,nyin,nzin,nzsoilin
  INTEGER :: nstyp1

  INTEGER :: istat, sd_id
  INTEGER :: varflg, istatus

  INTEGER :: idummy
  REAL    :: rdummy

  CHARACTER (LEN=40) :: fmtverin
  CHARACTER (LEN=12) :: label

  CHARACTER (LEN=256) :: filename
  INTEGER :: ia, ja, iloc, jloc
  INTEGER :: nxdta, nydta
  INTEGER :: nxpatch, nypatch

  LOGICAL :: success

!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------

  INCLUDE 'grid.inc'        ! Grid parameters
  INCLUDE 'indtflg.inc'
  INCLUDE 'mp.inc'

  CHARACTER(LEN=80) :: runname
  INTEGER :: nocmnt
  CHARACTER(LEN=80) :: cmnt(50)
  INTEGER :: year, month, day, hour, minute, second
  REAL    :: umove, vmove

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  IF (mp_opt > 0) THEN

    IF ( MOD(nprocx_in, nproc_x) /= 0 .OR. MOD(nprocy_in, nproc_y) /= 0) THEN
      IF (myproc == 0) WRITE(6,'(3x,a/,2(3x,2(a,I2)/))')                &
      'nprocx_in (nprocy_in) must be a multiplier of nproc_x(nproc_y)', &
      'nprocx_in = ',nprocx_in, ', nprocy_in = ',nprocy_in,             &
      'nproc_x   = ',nproc_x,   ', nproc_y   = ', nproc_y
      CALL arpsstop('ERROR: command line parameter.',1);
    END IF
    nxpatch = nprocx_in / nproc_x
    nypatch = nprocy_in / nproc_y
  ELSE
    nxpatch = nprocx_in
    nypatch = nprocy_in
  END IF

  nxdta = (nx-3)/nxpatch+3
  nydta = (ny-3)/nypatch+3

  ALLOCATE(invar3d(nxdta,nydta,nz),         STAT = istat)
  CALL check_alloc_status(istat, "READUVT:invar3d")
  IF (filefmt == 3) THEN
    ALLOCATE(itmp(nxdta,nydta,nz), STAT = istat)
    CALL check_alloc_status(istat, "READUVT:itmp")
    ALLOCATE(hmin(nz),             STAT = istat)
    CALL check_alloc_status(istat, "READUVT:hmin")
    ALLOCATE(hmax(nz),             STAT = istat)
    CALL check_alloc_status(istat, "READUVT:hmax")
  END IF

  IF (myproc == 0) WRITE(6,'(/1x,2a)') 'READUVT: Reading file: ', trim(filenamein)

!-----------------------------------------------------------------------
! Open file for reading
!-----------------------------------------------------------------------

  DO jloc = 1, nypatch
    DO iloc = 1, nxpatch

      IF ( mp_opt > 0 .OR. (nxpatch > 1 .OR. nypatch > 1) ) THEN
        CALL gtsplitfn(filenamein,nxpatch,nypatch,loc_x,loc_y,iloc,jloc,&
                       0,0,1,2,filename,istatus)
      ELSE
        WRITE(filename,'(a)') TRIM(filenamein)
      END IF

      IF (filefmt == 1) THEN

        CALL getunit( nchanl )

        OPEN(UNIT=nchanl,FILE=TRIM(filename),                           &
             STATUS='old',FORM='unformatted',ACTION='READ',IOSTAT=istat)

        READ(nchanl,ERR=110,END=120) fmtverin

        IF (myproc == 0) WRITE(6,'(1x,a,a/)')                           &
          'Incoming data format, fmtverin = ',TRIM(fmtverin)

        READ(nchanl,ERR=110,END=120) runname
        READ(nchanl,ERR=110,END=120) nocmnt
        IF( nocmnt > 0 ) THEN
          DO i=1,nocmnt
            READ(nchanl,ERR=110,END=120) cmnt(i)
          END DO
        END IF

        IF (myproc == 0) WRITE(6,'(1x,2a,/)') &
          'THE NAME OF THIS RUN IS: ', TRIM(runname)

        IF( nocmnt > 0 .AND. myproc == 0 ) THEN
          DO i=1,nocmnt
            WRITE(6,'(1x,a)') cmnt(i)
          END DO
        END IF

        READ(nchanl,ERR=110,END=120) time,tmunit

        READ(nchanl,ERR=110,END=120) nxin, nyin, nzin,nzsoilin

        IF ( nxin /= nxdta .OR. nyin /= nydta .OR. nzin /= nz ) THEN
          WRITE(6,'(1x,a)') ' Dimensions in HDFREADUVW inconsistent with data.'
          WRITE(6,'(1x,a,3I15)') ' Read were: ', nxin,nyin,nzin
          WRITE(6,'(1x,a,3I15)') ' Expected:  ', nxdta, nydta, nz
          WRITE(6,'(1x,a)') ' Program aborted in HDFREAD.'
          CALL arpsstop('arpsstop called from HDFREADUVW due to nxin...',1)
        END IF

        IF (myproc == 0) WRITE(lchanl,'(1x,a,f8.1,a,f8.3,a/)')  &
           'To read data for time:',time,' secs = ',(time/60.),' mins.'

        READ(nchanl,ERR=110,END=120)                                    &
                  grdin,basin,varin,mstin,idummy,                       &
                  idummy, idummy,idummy,landin,totin,                   &
                  idummy,idummy,idummy,mapproj, month,                  &
                  day, year, hour, minute, second

        READ(nchanl,ERR=110,END=120)                                    &
                  umove,   vmove,xgrdorg,ygrdorg, rdummy,               &
                  rdummy, rdummy, rdummy, rdummy, rdummy,               &
                  rdummy, rdummy, rdummy, rdummy, rdummy,               &
                  rdummy, rdummy, rdummy, rdummy, rdummy

        IF ( totin /= 0 ) THEN
          READ(nchanl,ERR=110,END=120)                                  &
                 nstyp1,idummy,idummy,idummy,idummy,                    &
                 idummy,idummy,idummy,idummy,idummy,                    &
                 idummy,idummy,idummy,idummy,idummy,                    &
                 idummy,idummy,idummy,idummy,idummy

          IF ( nstyp1 < 1 ) nstyp1 = 1

          READ(nchanl,ERR=110,END=120)                                  &
                 rdummy,rdummy,rdummy,rdummy,rdummy,                    &
                 rdummy,rdummy,rdummy,rdummy,rdummy,                    &
                 rdummy,rdummy,rdummy,rdummy,rdummy,                    &
                 rdummy,rdummy,rdummy,rdummy,rdummy
        END IF

        IF( grdin == 1 ) THEN
          READ(nchanl,ERR=110,END=120) label
          READ(nchanl,ERR=110,END=120)
          IF (myproc == 0) WRITE(lchanl,920) label

          READ(nchanl,ERR=110,END=120) label
          READ(nchanl,ERR=110,END=120)
          IF (myproc == 0) WRITE(lchanl,920) label

          READ(nchanl,ERR=110,END=120) label
          READ(nchanl,ERR=110,END=120)
          IF (myproc == 0) WRITE(lchanl,920) label

          READ(nchanl,ERR=110,END=120) label
          READ(nchanl,ERR=110,END=120)
          IF (myproc == 0) WRITE(lchanl,920) label

          READ(nchanl,ERR=110,END=120) label
          READ(nchanl,ERR=110,END=120)
          IF (myproc == 0) WRITE(lchanl,920) label
        END IF  ! grdin

        IF( basin == 1 ) THEN

          READ(nchanl,ERR=110,END=120) label
          READ(nchanl,ERR=110,END=120)
          IF (myproc == 0) WRITE(lchanl,920) label

          READ(nchanl,ERR=110,END=120) label
          READ(nchanl,ERR=110,END=120)
          IF (myproc == 0) WRITE(lchanl,920) label

          READ(nchanl,ERR=110,END=120) label
          READ(nchanl,ERR=110,END=120)
          IF (myproc == 0) WRITE(lchanl,920) label

          READ(nchanl,ERR=110,END=120) label
          READ(nchanl,ERR=110,END=120)
          IF (myproc == 0) WRITE(lchanl,920) label

          READ(nchanl,ERR=110,END=120) label
          READ(nchanl,ERR=110,END=120)
          IF (myproc == 0) WRITE(lchanl,920) label

          IF( mstin == 1) THEN
            READ(nchanl,ERR=110,END=120) label
            READ(nchanl,ERR=110,END=120)
            IF (myproc == 0) WRITE(lchanl,920) label
          END IF

          IF (landin == 1) THEN

            DO is=1,nstyp1
              READ(nchanl,ERR=110,END=120) label
              READ(nchanl,ERR=110,END=120)
              IF (myproc == 0) WRITE(lchanl,920) label

              READ(nchanl,ERR=110,END=120) label
              READ(nchanl,ERR=110,END=120)
              IF (myproc == 0) WRITE(lchanl,920) label
            END DO

            READ(nchanl,ERR=110,END=120) label
            READ(nchanl,ERR=110,END=120)
            IF (myproc == 0) WRITE(lchanl,920) label

            READ(nchanl,ERR=110,END=120) label
            READ(nchanl,ERR=110,END=120)
            IF (myproc == 0) WRITE(lchanl,920) label

            READ(nchanl,ERR=110,END=120) label
            READ(nchanl,ERR=110,END=120)
            IF (myproc == 0) WRITE(lchanl,920) label

            READ(nchanl,ERR=110,END=120) label
            READ(nchanl,ERR=110,END=120)
            IF (myproc == 0) WRITE(lchanl,920) label

          END IF
        END IF

        IF( varin == 1 ) THEN

          IF ( totin == 0 ) THEN
            WRITE(lchanl,'(1x,3(a/))')   &
        'We are expecting total variables in the file.',     &
        'You can generate file with total variables by set totout = 1 in the namelist.',&
        'Program stoping ...'
            CALL arpsstop('total variable expected.',1)
          END IF
!
!-----------------------------------------------------------------------
!
!  Read in total values of variables from history dump
!
!----------------------------------------------------------------------
!
          READ(nchanl,ERR=110,END=120) label
          READ(nchanl,ERR=110,END=120) invar3d
          IF (myproc == 0) WRITE(lchanl,910) label,' u.'
          DO k = 1,nz
            DO j = 1, nyin
              ja = (jloc-1)*(nyin-3)+j
              DO i = 1, nxin
                ia = (iloc-1)*(nxin-3)+i
                u(ia,ja,k) = invar3d(i,j,k)
              END DO
            END DO
          END DO

          READ(nchanl,ERR=110,END=120) label
          READ(nchanl,ERR=110,END=120) invar3d
          IF (myproc == 0) WRITE(lchanl,910) label,' v.'
          DO k = 1,nz
            DO j = 1, nyin
              ja = (jloc-1)*(nyin-3)+j
              DO i = 1, nxin
                ia = (iloc-1)*(nxin-3)+i
                v(ia,ja,k) = invar3d(i,j,k)
              END DO
            END DO
          END DO

          READ(nchanl,ERR=110,END=120) label
          READ(nchanl,ERR=110,END=120)
          IF (myproc == 0) WRITE(lchanl,920) label

          READ(nchanl,ERR=110,END=120) label
          READ(nchanl,ERR=110,END=120) invar3d
          IF (myproc == 0) WRITE(lchanl,910) label,' pt.'
          DO k = 1,nz
            DO j = 1, nyin
              ja = (jloc-1)*(nyin-3)+j
              DO i = 1, nxin
                ia = (iloc-1)*(nxin-3)+i
                pt(ia,ja,k) = invar3d(i,j,k)
              END DO
            END DO
          END DO

          READ(nchanl,ERR=110,END=120) label
          READ(nchanl,ERR=110,END=120) invar3d
          IF (myproc == 0) WRITE(lchanl,910) label,' p.'
          DO k = 1,nz
            DO j = 1, nyin
              ja = (jloc-1)*(nyin-3)+j
              DO i = 1, nxin
                ia = (iloc-1)*(nxin-3)+i
                p(ia,ja,k) = invar3d(i,j,k)
              END DO
            END DO
          END DO

        END IF
!
!-----------------------------------------------------------------------
!
!  Read in moisture variables
!
!-----------------------------------------------------------------------
!
        IF( mstin == 1 ) THEN

          READ(nchanl,ERR=110,END=120) label
          READ(nchanl,ERR=110,END=120) invar3d
          IF (myproc == 0) WRITE(lchanl,910) label,' qv.'
          DO k = 1,nz
            DO j = 1, nyin
              ja = (jloc-1)*(nyin-3)+j
              DO i = 1, nxin
                ia = (iloc-1)*(nxin-3)+i
                qv(ia,ja,k) = invar3d(i,j,k)
              END DO
            END DO
          END DO

        END IF

        GOTO 130

        110   CONTINUE
        WRITE(6,'(/a/)') ' Error reading data in READUV'
        ireturn=-1
        RETURN

        120   CONTINUE
        WRITE(6,'(/a/)') ' End of file reached in READUV'
        ireturn=-2
        RETURN

        130   CONTINUE
        CLOSE(UNIT=nchanl)
        CALL retunit( nchanl )
        WRITE(6,*)

        910   FORMAT(1X,'Field ',a12,' was read into array',a)
        920   FORMAT(1X,'Field ',a12,' was skipped.')

      ELSE IF (filefmt == 3) THEN

        CALL hdfopen(filename,1,sd_id)
        IF (sd_id < 0) THEN
          WRITE (6,*) "HDFREADUV: ERROR opening ", trim(filename)," for reading."
          ireturn=-1
          RETURN
        END IF

        IF (nxpatch == 1 .AND. nypatch == 1) THEN
          CALL hdfrdc(sd_id,80,"runname",runname,istat)
          CALL hdfrdi(sd_id,"nocmnt",nocmnt,istat)
          IF( nocmnt > 0 ) THEN
            CALL hdfrdc(sd_id,80*nocmnt,"cmnt",cmnt,istat)
          END IF

          IF (myproc == 0) WRITE(6,'(1x,2a)')     &
            ' THE NAME OF THIS RUN IS:  ', trim(runname)

!          IF( nocmnt > 0 .AND. myproc == 0) THEN
!            WRITE (6,*) "Comments:"
!            DO i=1,nocmnt
!              WRITE(6,'(1x,a)') cmnt(i)
!            END DO
!          END IF
!
!          IF (myproc == 0) WRITE (6,*) " "

          CALL hdfrdc(sd_id,10,"tmunit",tmunit,istat)
        END IF

!-----------------------------------------------------------------------
!  Get dimensions of data in binary file and check against
!  the dimensions passed to HDFREAD
!-----------------------------------------------------------------------
        CALL hdfrdr(sd_id,"time",time,istat)

        CALL hdfrdi(sd_id,"nx",nxin,istat)
        CALL hdfrdi(sd_id,"ny",nyin,istat)
        CALL hdfrdi(sd_id,"nz",nzin,istat)

        IF ( nxin /= nxdta .OR. nyin /= nydta .OR. nzin /= nz ) THEN
          IF (myproc == 0) THEN
            WRITE(6,'(1x,a)') ' Dimensions in HDFREADUVW inconsistent with data.'
            WRITE(6,'(1x,a,3I15)') ' Read were: ', nxin, nyin, nzin
            WRITE(6,'(1x,a,3I15)') ' Expected:  ', nxdta, nydta, nz
            WRITE(6,'(1x,a)') ' Program aborted in READUVT.'
          END IF
          CALL arpsstop('arpsstop called from READUVW due to nxin...',1)
        END IF

        IF (myproc == 0) WRITE(lchanl,'(1x,a,f8.1,a,f8.3,a/)')  &
           'To read data for time:', time,' secs = ',(time/60.),' mins.'

        CALL hdfrdi(sd_id,"grdflg",grdin,istat)
        CALL hdfrdi(sd_id,"basflg",basin,istat)
        CALL hdfrdi(sd_id,"varflg",varin,istat)

!        IF (myproc == 0) WRITE(6,'(a)') ' Done reading parameters'

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

        success=.true.
        IF( varin == 1 ) then

!-----------------------------------------------------------------------
!  Read in total values of variables from history dump
!-----------------------------------------------------------------------

          CALL hdfrd3d(sd_id,"u",nxin,nyin,nz,invar3d,istat,itmp,hmax,hmin)
          IF (istat == 0) THEN
            DO k = 1,nz
              DO j = 1, nyin
                ja = (jloc-1)*(nyin-3)+j
                DO i = 1, nxin
                  ia = (iloc-1)*(nxin-3)+i
                  u(ia,ja,k) = invar3d(i,j,k)
                END DO
              END DO
            END DO
          ELSE
            success=.false.
          END IF

          CALL hdfrd3d(sd_id,"v",nxin,nyin,nz,invar3d,istat,itmp,hmax,hmin)
          IF (istat == 0) THEN
            DO k = 1,nz
              DO j = 1, nyin
                ja = (jloc-1)*(nyin-3)+j
                DO i = 1, nxin
                  ia = (iloc-1)*(nxin-3)+i
                  v(ia,ja,k) = invar3d(i,j,k)
                END DO
              END DO
            END DO
          ELSE
            success=.false.
          END IF

          CALL hdfrd3d(sd_id,"p",nxin,nyin,nz,invar3d,istat,itmp,hmax,hmin)
          IF (istat == 0) THEN
            DO k = 1,nz
              DO j = 1, nyin
                ja = (jloc-1)*(nyin-3)+j
                DO i = 1, nxin
                  ia = (iloc-1)*(nxin-3)+i
                  p(ia,ja,k) = invar3d(i,j,k)
                END DO
              END DO
            END DO
          ELSE
            success=.false.
          END IF

          CALL hdfrd3d(sd_id,"pt",nxin,nyin,nz,invar3d,istat,itmp,hmax,hmin)
          IF (istat == 0) THEN
            DO k = 1,nz
              DO j = 1, nyin
                ja = (jloc-1)*(nyin-3)+j
                DO i = 1, nxin
                  ia = (iloc-1)*(nxin-3)+i
                  pt(ia,ja,k) = invar3d(i,j,k)
                END DO
              END DO
            END DO
          ELSE
            success=.false.
          END IF

          CALL hdfrd3d(sd_id,"qv",nxin,nyin,nz,invar3d,istat,itmp,hmax,hmin)
          IF (istat == 0) THEN
            DO k = 1,nz
              DO j = 1, nyin
                ja = (jloc-1)*(nyin-3)+j
                DO i = 1, nxin
                  ia = (iloc-1)*(nxin-3)+i
                  qv(ia,ja,k) = invar3d(i,j,k)
                END DO
              END DO
            END DO
          ELSE
            success=.false.
          END IF

!-----------------------------------------------------------------------
!
!  Friendly exit message
!
!-----------------------------------------------------------------------

          IF(success) THEN
            IF (myproc == 0) WRITE(6,'(/a,F8.1,a/)')                    &
          ' Data at time=', time/60,' (min) were successfully read.'

            ireturn = 0
          ELSE
            WRITE(6,'(/a,F8.1,a/)')                                     &
          ' Error reading u,v data at time=', time/60,' (min)'
            ireturn = -3
            RETURN
          END IF
          CALL hdfclose(sd_id,istat)
        ELSE
          WRITE(6,'(/a/a,a/)') ' Error reading data in READUV:',        &
        ' No time-dependent data in file',trim(filename)
          ireturn = -2
          CALL hdfclose(sd_id,istat)
          RETURN
        END IF

      ELSE
        WRITE(6,'(1x,a,I2,a,/)') 'The file formt - ',filefmt,' is still not supported.'
        ireturn = -9
        RETURN
      END IF

    END DO
  END DO

  DEALLOCATE(invar3d)
  IF (filefmt == 3) DEALLOCATE(itmp, hmin, hmax)

  RETURN
END SUBROUTINE readuvt
!
!##################################################################
!##################################################################
!######                                                      ######
!######                  SUBROUTINE MAPINIT                  ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE mapinit(nx,ny)
!
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx,ny

  REAL :: alatpro(2)
  REAL :: sclf,dxscl,dyscl
  REAL :: ctrx,ctry,swx,swy

  INCLUDE 'grid.inc'
  INCLUDE 'mp.inc'

  alatpro(1) = trulat1
  alatpro(2) = trulat2

  IF( sclfct /= 1.0) THEN
    sclf  = 1.0/sclfct
    dxscl = dx*sclf
    dyscl = dy*sclf
  ELSE
    sclf  = 1.0
    dxscl = dx
    dyscl = dy
  END IF

  CALL setmapr( mapproj,sclf,alatpro,trulon )

  CALL lltoxy( 1,1, ctrlat,ctrlon, ctrx, ctry )

!       swx = ctrx - (float(nx-3)/2.) * dxscl
!       swy = ctry - (float(ny-3)/2.) * dyscl
  swx = ctrx - (FLOAT(nproc_x*(nx-3))/2.) * dxscl
  swy = ctry - (FLOAT(nproc_y*(ny-3))/2.) * dyscl

  CALL setorig( 1, swx, swy)

  RETURN
END SUBROUTINE mapinit
!
!##################################################################
!##################################################################
!######                                                      ######
!######               SUBROUTINE WTRADHDRBIN                 ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE wtradhdrbin(rfname,radid,itime,vcpnum,isource,iradfmt,      &
                numradcol,nummaxlev,typelev,dualpol,                   &
                irngmin,irngmax,latrad,lonrad,elvrad,                  &
                refelvmin,refelvmax,mxradvr,nradvr,iradvr,             &
                iunit,fileopnd,istatus)
  IMPLICIT NONE
  CHARACTER(LEN=256), INTENT(IN) :: rfname
  CHARACTER(LEN=4),   INTENT(IN) :: radid
  INTEGER, INTENT(IN)  :: itime
  INTEGER, INTENT(IN)  :: vcpnum
  INTEGER, INTENT(IN)  :: isource
  INTEGER, INTENT(IN)  :: iradfmt
  INTEGER, INTENT(IN)  :: numradcol
  INTEGER, INTENT(IN)  :: nummaxlev
  INTEGER, INTENT(IN)  :: typelev
  INTEGER, INTENT(IN)  :: dualpol
  INTEGER, INTENT(IN)  :: irngmin
  INTEGER, INTENT(IN)  :: irngmax
  REAL,    INTENT(IN)  :: latrad
  REAL,    INTENT(IN)  :: lonrad
  REAL,    INTENT(IN)  :: elvrad
  REAL,    INTENT(IN)  :: refelvmin
  REAL,    INTENT(IN)  :: refelvmax
  INTEGER, INTENT(IN)  :: mxradvr
  INTEGER, INTENT(IN)  :: nradvr
  INTEGER, INTENT(IN)  :: iradvr(mxradvr)
  INTEGER, INTENT(OUT) :: iunit
  LOGICAL, INTENT(OUT) :: fileopnd
  INTEGER, INTENT(OUT) :: istatus

  INTEGER, PARAMETER   :: idummy = 0
  REAL,    PARAMETER   :: rdummy = 0

  INCLUDE 'globcst.inc'
  INCLUDE 'grid.inc'
  INCLUDE 'remap.inc'

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  CALL getunit(iunit)
!
!-----------------------------------------------------------------------
!
!  Open file for output
!
!-----------------------------------------------------------------------
!
  OPEN(iunit,FILE=TRIM(rfname),STATUS='UNKNOWN',FORM='UNFORMATTED')
  fileopnd=.TRUE.
!
!-----------------------------------------------------------------------
!
!  Write radar description variables
!
!-----------------------------------------------------------------------
!
  WRITE(iunit) radid
  WRITE(iunit) ireftim,itime,vcpnum,isource,idummy,                     &
               idummy,idummy,idummy,idummy,idummy
!
!-----------------------------------------------------------------------
!
!  Write grid description variables
!  This should provide enough info to verify that the
!  proper grid has been chosen.  To recreate the grid,
!  including elevation information,
!  the reading program should get a grid-base-file
!  named runname.grdbasfil
!
!-----------------------------------------------------------------------
!
  WRITE(iunit) runname
  WRITE(iunit) iradfmt,strhopt,mapproj,irngmin,irngmax,                 &
               typelev,numradcol,nummaxlev,dualpol,idummy
  WRITE(iunit) dx,dy,dz,dzmin,ctrlat,                                   &
               ctrlon,trulat1,trulat2,trulon,sclfct,                    &
               latrad,lonrad,elvrad,refelvmin,refelvmax
  WRITE(iunit) nradvr,iradvr
  istatus=0

  RETURN
END SUBROUTINE wtradhdrbin

SUBROUTINE wtradhdrhdf(rfname,radid,itime,vcpnum,isource,iradfmt,   &
                       nxlg,nylg,nz,typelev,dualpol,                &
                       irngmin,irngmax,latrad,lonrad,elvrad,        &
                       refelvmin,refelvmax,mxradvr,nradvr,iradvr,   &
                       sd_id,fileopnd,istatus)
!
!
! Write header variables to HDF radar column file
!
! Keith Brewster, CAPS
! 02/08/2010
!
  IMPLICIT NONE
  CHARACTER(LEN=256), INTENT(IN) :: rfname
  CHARACTER(LEN=4), INTENT(IN) :: radid
  INTEGER, INTENT(IN)  :: itime
  INTEGER, INTENT(IN)  :: vcpnum
  INTEGER, INTENT(IN)  :: isource
  INTEGER, INTENT(IN)  :: iradfmt
  INTEGER, INTENT(IN)  :: nxlg
  INTEGER, INTENT(IN)  :: nylg
  INTEGER, INTENT(IN)  :: nz
  INTEGER, INTENT(IN)  :: typelev
  INTEGER, INTENT(IN)  :: dualpol
  INTEGER, INTENT(IN)  :: irngmin
  INTEGER, INTENT(IN)  :: irngmax
  REAL,    INTENT(IN)  :: latrad
  REAL,    INTENT(IN)  :: lonrad
  REAL,    INTENT(IN)  :: elvrad
  REAL,    INTENT(IN)  :: refelvmin
  REAL,    INTENT(IN)  :: refelvmax
  INTEGER, INTENT(IN)  :: mxradvr
  INTEGER, INTENT(IN)  :: nradvr
  INTEGER, INTENT(IN)  :: iradvr(mxradvr)
  INTEGER, INTENT(OUT) :: sd_id
  LOGICAL, INTENT(OUT) :: fileopnd
  INTEGER, INTENT(OUT) :: istatus

  INTEGER, PARAMETER :: idummy = 0
  REAL, PARAMETER :: rdummy = 0

  INTEGER :: istat

  INCLUDE 'globcst.inc'
  INCLUDE 'grid.inc'
  INCLUDE 'remap.inc'

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!-----------------------------------------------------------------------
!
!  Open file for output
!
!-----------------------------------------------------------------------
!
    CALL hdfopen(trim(rfname), 2, sd_id)
    IF (sd_id < 0) THEN
      WRITE (6,*) "WTRADCOL: ERROR opening ",                           &
                  trim(rfname)," for writing."
      istatus = -1
      STOP
    END IF
    fileopnd=.TRUE.
!
!-----------------------------------------------------------------------
!
!  Write radar description variables
!
!-----------------------------------------------------------------------
!
    CALL hdfwrtc(sd_id, 4, 'radid', radid, istat)
    CALL hdfwrti(sd_id, 'ireftim', ireftim, istat)
    CALL hdfwrti(sd_id, 'itime', itime, istat)
    CALL hdfwrti(sd_id, 'vcpnum', vcpnum, istat)
    CALL hdfwrti(sd_id, 'isource', isource, istat)
!
!-----------------------------------------------------------------------
!
!  Write grid description variables
!  This should provide enough info to verify that the
!  proper grid has been chosen.  To recreate the grid,
!  icluding elevation information,
!  the reading program should get a grid-base-file
!  named runname.grdbasfil
!
!-----------------------------------------------------------------------
!
    CALL hdfwrtc(sd_id, 4, 'runname', runname, istat)
    CALL hdfwrti(sd_id, 'iradfmt', iradfmt, istat)
    CALL hdfwrti(sd_id, 'strhopt', strhopt, istat)
    CALL hdfwrti(sd_id, 'mapproj', mapproj, istat)
    CALL hdfwrti(sd_id, 'nx', nxlg, istat)
    CALL hdfwrti(sd_id, 'ny', nylg, istat)
    CALL hdfwrti(sd_id, 'nz', nz, istat)
    CALL hdfwrtr(sd_id, 'dx', dx, istat)
    CALL hdfwrtr(sd_id, 'dy', dy, istat)
    CALL hdfwrtr(sd_id, 'dz', dz, istat)
    CALL hdfwrtr(sd_id, 'dzmin', dzmin, istat)
    CALL hdfwrtr(sd_id, 'ctrlat', ctrlat, istat)
    CALL hdfwrtr(sd_id, 'ctrlon', ctrlon, istat)
    CALL hdfwrtr(sd_id, 'trulat1', trulat1, istat)
    CALL hdfwrtr(sd_id, 'trulat2', trulat2, istat)
    CALL hdfwrtr(sd_id, 'trulon', trulon, istat)
    CALL hdfwrtr(sd_id, 'sclfct', sclfct, istat)
    CALL hdfwrtr(sd_id, 'latrad', latrad, istat)
    CALL hdfwrtr(sd_id, 'lonrad', lonrad, istat)
    CALL hdfwrtr(sd_id, 'elvrad', elvrad, istat)
    CALL hdfwrti(sd_id, 'irngmin', irngmin, istat)
    CALL hdfwrti(sd_id, 'irngmax', irngmax, istat)
    CALL hdfwrtr(sd_id, 'refelvmin', refelvmin, istat)
    CALL hdfwrtr(sd_id, 'refelvmax', refelvmax, istat)
    CALL hdfwrti(sd_id, 'nradvr', nradvr, istat)
    CALL hdfwrt1di(iradvr,mxradvr,sd_id,'iradvr','iradvr','null')
    CALL hdfwrti(sd_id, 'typelev', typelev, istat)
    CALL hdfwrti(sd_id, 'dualpol', dualpol, istat)
    istatus=istat
  RETURN
END SUBROUTINE wtradhdrhdf

SUBROUTINE distrallcol(nxny,nz,ncolp,ncoltot,                           &
                       icolp,jcolp,xcolp,ycolp,zcolp,lvldbg,istatus)
!#######################################################################
!
! Distribute all the location information of the non-missing
! radar columns among all the processors.
!
! Keith Brewster, CAPS
! 02/28/2010
!
! Yunheng Wang (04/20/2011)
! Reformed the distribution algorithm for efficiency and less message
! communications.
!
!#######################################################################

  IMPLICIT NONE
  INTEGER, INTENT(IN)    :: nxny
  INTEGER, INTENT(IN)    :: nz
  INTEGER, INTENT(INOUT) :: ncolp
  INTEGER, INTENT(IN)    :: ncoltot
  INTEGER, INTENT(INOUT) :: icolp(nxny)
  INTEGER, INTENT(INOUT) :: jcolp(nxny)
  REAL,    INTENT(INOUT) :: xcolp(nxny)
  REAL,    INTENT(INOUT) :: ycolp(nxny)
  REAL,    INTENT(INOUT) :: zcolp(nz,nxny)
  INTEGER, INTENT(IN)    :: lvldbg
  INTEGER, INTENT(OUT)   :: istatus

!-----------------------------------------------------------------------
!
! Temporary large arrays to be used by processor zero
!
!-----------------------------------------------------------------------

  !INTEGER, ALLOCATABLE :: icollg(:)
  !INTEGER, ALLOCATABLE :: jcollg(:)
  !REAL,    ALLOCATABLE :: xcollg(:)
  !REAL,    ALLOCATABLE :: ycollg(:)
  !REAL,    ALLOCATABLE :: zcollg(:,:)

  INTEGER, ALLOCATABLE :: numcol(:)   ! number of columns to be processed
  INTEGER, ALLOCATABLE :: procsq(:)   ! process sequence number corresponding to numcol
  INTEGER, ALLOCATABLE :: procmt(:,:) ! process mate
  INTEGER, ALLOCATABLE :: nsend(:)    ! number of columns to be send to each proces
                                      ! negative means receive from each proces
!-----------------------------------------------------------------------
!
! Misc local variables
!
!-----------------------------------------------------------------------

  !INTEGER :: idxbgn,idxend,ncolj,nsize2d
  !INTEGER :: ncpr
  !INTEGER :: i,k,jproc,ierr,istat

  INTEGER :: iloc, jloc, ierr
  INTEGER :: sendpnt, recvcnt, extracol, starveproc
  INTEGER, PARAMETER :: itag = 500

  INTEGER :: ncpr   ! number of columns when evenly distributed

  INCLUDE 'mp.inc'

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!
! Collect all columns in the "lg" arrays in processor zero.
!
!  IF( myproc == 0) THEN
!
!    WRITE(6,'(a,i12)') ' distrallcol: allocating temporary arrays'
    !ALLOCATE (icollg(ncoltot+1), stat=istat)
    !CALL check_alloc_status(istat,'distrallcol:icollg')
    !ALLOCATE (jcollg(ncoltot+1), stat=istat)
    !CALL check_alloc_status(istat,'distrallcol:jcollg')
    !ALLOCATE (xcollg(ncoltot+1), stat=istat)
    !CALL check_alloc_status(istat,'distrallcol:xcollg')
    !ALLOCATE (ycollg(ncoltot+1), stat=istat)
    !CALL check_alloc_status(istat,'distrallcol:ycollg')
    !ALLOCATE (zcollg(nz,ncoltot+1), stat=istat)
    !CALL check_alloc_status(istat,'distrallcol:zcollg')

    !DO i=1,ncolp
    !  icollg(i)=icolp(i)
    !  jcollg(i)=jcolp(i)
    !  xcollg(i)=xcolp(i)
    !  ycollg(i)=ycolp(i)
    !  DO k=1,nz
    !    zcollg(k,i)=zcolp(k,i)
    !  END DO
    !END DO
    !idxbgn=ncolp+1
    !
    !WRITE(6,'(a,i12)') ' distrallcol: transfer 0 complete at ncolp:',ncolp
    !
    !DO jproc=1,(nprocs-1)
    !  WRITE(6,'(a,i12)') ' distrallcol: getting ncolj from jproc:',jproc
    !  CALL mprecvi(ncolj, 1, jproc, itag, ierr )
    !  WRITE(6,'(a,i12)') ' distrallcol: got ncolj:',ncolj
    !  nsize2d=ncolj*nz
    !  CALL mprecvi(icollg(idxbgn),ncolj,jproc,(itag+1),ierr)
    !  CALL mprecvi(jcollg(idxbgn),ncolj,jproc,(itag+2),ierr)
    !  CALL mprecvr(xcollg(idxbgn),ncolj,jproc,(itag+3),ierr)
    !  CALL mprecvr(ycollg(idxbgn),ncolj,jproc,(itag+4),ierr)
    !  CALL mprecvr(zcollg(1,idxbgn),nsize2d,jproc,(itag+5),ierr)
    !  idxbgn=idxbgn+ncolj
    !
    !  WRITE(6,'(a,i10,a,i12)') ' distrallcol: jproc:',jproc,' done knt=',(idxbgn-1)
    !
    !END DO
    !
    !WRITE(6,'(a,i12)') ' End of collection ncolumns: ',(idxbgn-1)

  !ELSE

    !nsize2d=ncolp*nz
    !CALL mpsendi(ncolp, 1,       0, itag, ierr )
    !CALL mpsendi(icolp, ncolp,   0, (itag+1), ierr )
    !CALL mpsendi(jcolp, ncolp,   0, (itag+2), ierr )
    !CALL mpsendr(xcolp, ncolp,   0, (itag+3), ierr )
    !CALL mpsendr(ycolp, ncolp,   0, (itag+4), ierr )
    !CALL mpsendr(zcolp, nsize2d, 0, (itag+5), ierr )

!  END IF

!!
!! Now redistribute the columns to the processors evenly.
!! Reset ncolp accordingly.
!!
!  ncpr=ncoltot/nprocs
!  IF(mod(ncoltot,ncpr) > 0) ncpr=ncpr+1
!  idxbgn=(myproc*ncpr)+1
!  idxend=min((idxbgn+(ncpr-1)),ncoltot)
!  ncolp=(idxend-idxbgn)+1
!  WRITE(6,'(a,i6,a,i10,a,i10,a,i10)')                                &
!        ' distrallcol myproc:',myproc,' ncolp:',ncolp,               &
!        ' idxbgn:',idxbgn,' to idxend:',idxend
!
!  IF( myproc == 0) THEN
!    DO i=idxbgn,idxend
!      icolp(i)=icollg(i)
!      jcolp(i)=jcollg(i)
!      xcolp(i)=xcollg(i)
!      ycolp(i)=ycollg(i)
!      DO k=1,nz
!        zcolp(k,i)=zcollg(k,i)
!      END DO
!    END DO
!
!    DO jproc=1,(nprocs-1)
!      idxbgn=(jproc*ncpr)+1
!      idxend=min((idxbgn+(ncpr-1)),ncoltot)
!      ncolj=(idxend-idxbgn)+1
!      nsize2d=ncolj*nz
!      CALL mpsendi(icollg(idxbgn),ncolj,jproc,(itag+6),ierr)
!      CALL mpsendi(jcollg(idxbgn),ncolj,jproc,(itag+7),ierr)
!      CALL mpsendr(xcollg(idxbgn),ncolj,jproc,(itag+8),ierr)
!      CALL mpsendr(ycollg(idxbgn),ncolj,jproc,(itag+9),ierr)
!      CALL mpsendr(zcollg(1,idxbgn),nsize2d,jproc,(itag+10),ierr)
!    END DO
!
!    WRITE(6,'(a)') ' distrallcol: Processor ZERO done SEND.'
!
!    DEALLOCATE( icollg,jcollg,xcollg,ycollg,zcollg )
!
!  ELSE
!
!    nsize2d=ncolp*nz
!    CALL mprecvi(icolp, ncolp, 0, (itag+6), ierr )
!    CALL mprecvi(jcolp, ncolp, 0, (itag+7), ierr )
!    CALL mprecvr(xcolp, ncolp, 0, (itag+8), ierr )
!    CALL mprecvr(ycolp, ncolp, 0, (itag+9), ierr )
!    CALL mprecvr(zcolp, nsize2d, 0, (itag+10), ierr )
!
!    WRITE(6,'(a,i10,a)') ' distrallcol: Processor ',myproc,' done RECV.'
!
!  END IF

!-----------------------------------------------------------------------
!
! New code, total rewrote by Y. Wang on 4/20/2011.
!
!-----------------------------------------------------------------------

  WRITE(6,'(1x,a,I4,a,i10,a)')                                          &
        'Entering distrallcol, myproc: ',myproc,' ncolp: ',ncolp,'.'
  CALL flush(6)

  IF( myproc == 0) THEN

    !WRITE(6,'(a,i12)') ' distrallcol: allocating temporary arrays'

    ALLOCATE (numcol(nprocs), stat=istatus)
    CALL check_alloc_status(istatus,'distrallcol:numcol')
    numcol(:) = 0
    ALLOCATE (procsq(nprocs), stat=istatus)
    CALL check_alloc_status(istatus,'distrallcol:procsq')
    !procsq(:) = -1
    ALLOCATE (procmt(nprocs,nprocs), stat=istatus)
    CALL check_alloc_status(istatus,'distrallcol:procmt')
    procmt(:,:) = 0       ! send and receiv nothing

  END IF

  IF (myproc == 0 .AND. lvldbg > 0) WRITE(6,'(a,i12)') ' distrallcol: beginning collection phase'

  CALL mpgatheri(ncolp,1,numcol,nprocs,istatus)

!-----------------------------------------------------------------------
!
! Redistribution metadata (only root process does the job)
!
!-----------------------------------------------------------------------

  ncpr=ncoltot/nprocs
  IF(mod(ncoltot,nprocs) > 0) ncpr=ncpr+1

  IF (myproc == 0) THEN
    DO iloc = 1, nprocs
      procsq(iloc) = iloc    ! Process sequence corresponding to numcol
    END DO

    !
    ! Sorting the arrays
    !
    IF (lvldbg > 2) THEN
      WRITE(6,'(1x,a)') 'Before quick sort: '
      DO iloc = 1, nprocs
        WRITE(6,'(3x,2(a,I10))') 'myproc = ',procsq(iloc),', ncolp = ',numcol(iloc)
      END DO
      CALL flush(6)
    END IF

    CALL qsortisq(numcol,nprocs,procsq)

    IF (lvldbg > 2) THEN
      WRITE(6,'(1x,a)') 'After quick sort: '
      DO iloc = 1, nprocs
        WRITE(6,'(3x,2(a,I10))') 'myproc = ',procsq(iloc),', ncolp = ',numcol(iloc)
      END DO
      CALL flush(6)
    END IF

    !
    ! Evaluate column passing
    !
    starveproc = 1   ! search from the tail that contains extra columns
    DO jloc = nprocs,2,-1

      DO WHILE (numcol(jloc) > ncpr)
        extracol = numcol(jloc) - ncpr

        DO iloc = starveproc, nprocs-1
          recvcnt = ncpr - numcol(iloc)
          IF (recvcnt <= 0) CYCLE

          IF ( extracol <= recvcnt) THEN  ! iloc is still starve and need to find next process with more columns
            numcol(iloc) = numcol(iloc) + extracol
            numcol(jloc) = numcol(jloc) - extracol
            procmt(procsq(jloc),procsq(iloc)) = -1*extracol       ! iloc recv from nloc
            procmt(procsq(iloc),procsq(jloc)) = extracol          ! jloc send to iloc
            IF (lvldbg > 2) THEN
              WRITE(6,'(1x,5(a,I6))') 'Proc ',jloc,' = ',procsq(jloc),' send ',extracol,' columns to ',iloc,' = ',procsq(iloc)
              CALL flush(6)
            END IF
            EXIT
          ELSE                           ! iloc is full, find next hungry one
            numcol(iloc) = numcol(iloc) + recvcnt
            numcol(jloc) = numcol(jloc) - recvcnt
            procmt(procsq(jloc),procsq(iloc)) = -1*recvcnt       ! iloc recv from jloc
            procmt(procsq(iloc),procsq(jloc)) = recvcnt          ! jloc send to iloc
            extracol = extracol - recvcnt
            starveproc = iloc+1        ! since we have filled iloc
            IF (lvldbg > 2) THEN
              WRITE(6,'(1x,5(a,I6))') 'Proc ',jloc,' = ',procsq(jloc),' send ',recvcnt,' columns to ',iloc,' = ',procsq(iloc)
              CALL flush(6)
            END IF
          END IF
        END DO

      END DO   ! while loop for iloc

    END DO

    IF (lvldbg > 1) THEN
      WRITE(6,'(1x,a)') 'Send/receive relationship, negative for receiving, positvie for sending'
      WRITE(6,'(1x,a,100I10)') '          ',(iloc-1,  iloc=1,nprocs)
      WRITE(6,'(1x,a,100A10)') '          ',(' ---------',iloc=1,nprocs)
      DO jloc = 1,nprocs
        WRITE(6,'(1x,a,I4,100I10)') 'Proc: ',jloc-1,(procmt(iloc,jloc),iloc=1,nprocs)
      END DO
      WRITE(6,'(1x,a)') '------------------------------------------------'
      CALL flush(6)
    END IF

  END IF

!-----------------------------------------------------------------------
!
! Redistribute radar columns
!
!-----------------------------------------------------------------------
  ALLOCATE (nsend(nprocs), stat=istatus)
  CALL check_alloc_status(istatus,'distrallcol:nsend')
  nsend(:) = 0

  CALL mpscatteri(procmt,nprocs*nprocs,nsend,nprocs,istatus)
  DO iloc = 1,nprocs
    IF (nsend(iloc) > 0) THEN       ! send to nloc
      sendpnt = ncolp-nsend(iloc)+1

      IF (lvldbg > 1) THEN
        WRITE(6,'(1x,3(a,I10))')                                          &
          'Process ',myproc,' sending ',nsend(iloc), ' columns to Process ',iloc-1
        CALL flush(6)
      END IF

      CALL mpsendi(icolp(sendpnt),     nsend(iloc), iloc-1, (itag+1), ierr )
      CALL mpsendi(jcolp(sendpnt),     nsend(iloc), iloc-1, (itag+2), ierr )
      CALL mpsendr(xcolp(sendpnt),     nsend(iloc), iloc-1, (itag+3), ierr )
      CALL mpsendr(ycolp(sendpnt),     nsend(iloc), iloc-1, (itag+4), ierr )
      CALL mpsendr(zcolp(1,sendpnt),nz*nsend(iloc), iloc-1, (itag+5), ierr )

      ncolp = ncolp-nsend(iloc)     ! since it has sent, its ncolp should reduce

    ELSE if (nsend(iloc) < 0) THEN  ! recv from nloc
      recvcnt = -1*nsend(iloc)

      IF (lvldbg > 1) THEN
        WRITE(6,'(1x,3(a,I10))')                                          &
          'Process ',myproc,' recving ',recvcnt, ' columns from Process ',iloc-1
        CALL flush(6)
      END IF

      CALL mprecvi(icolp(ncolp+1),      recvcnt, iloc-1, (itag+1), ierr )
      CALL mprecvi(jcolp(ncolp+1),      recvcnt, iloc-1, (itag+2), ierr )
      CALL mprecvr(xcolp(ncolp+1),      recvcnt, iloc-1, (itag+3), ierr )
      CALL mprecvr(ycolp(ncolp+1),      recvcnt, iloc-1, (itag+4), ierr )
      CALL mprecvr(zcolp(1,ncolp+1), nz*recvcnt, iloc-1, (itag+5), ierr )

      ncolp = ncolp+recvcnt         ! since it receives, its ncolp increased
    END IF
  END DO

  WRITE(6,'(1x,a,I4,a,i10,a)')                                          &
        'After distrallcol, myproc: ',myproc,' ncolp: ',ncolp,'.'
  CALL flush(6)

!-----------------------------------------------------------------------
!
! Ending the subroutine nicely
!
!-----------------------------------------------------------------------

  DEALLOCATE(nsend)
  IF (myproc == 0) DEALLOCATE(numcol, procsq, procmt)

  RETURN
END SUBROUTINE distrallcol

RECURSIVE SUBROUTINE qsortisq(iarr,nsize,indx)
!
! Sort an integer array in-place and keep its index array also.
!
  INTEGER, INTENT(IN)    :: nsize
  INTEGER, INTENT(INOUT) :: iarr(nsize), indx(nsize)

  INTEGER :: iq

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  if (nsize > 1) then
     call partition(iarr, nsize,indx, iq)
     call qsortisq(iarr(1),iq-1,indx(1))
     call qsortisq(iarr(iq),nsize-iq+1,indx(iq))
  endif

  RETURN
END SUBROUTINE qsortisq

SUBROUTINE Partition(A, nsize, indx, marker)
  INTEGER, INTENT(INOUT) :: A(nsize), indx(nsize)
  INTEGER, INTENT(OUT)   :: marker

  integer :: i, j
  INTEGER :: temv, temi
  INTEGER :: x      ! pivot point

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  x = A(1)
  i= 0
  j= nsize + 1

  do
     j = j-1
     do
        if (A(j) <= x) exit
        j = j-1
     end do
     i = i+1
     do
        if (A(i) >= x) exit
        i = i+1
     end do
     if (i < j) then
        ! exchange A(i) and A(j)
        temv = A(i)
        A(i) = A(j)
        A(j) = temv

        temi    = indx(i)
        indx(i) = indx(j)
        indx(j) = temi
     elseif (i == j) then
        marker = i+1
        return
     else
        marker = i
        return
     endif
  end do

END SUBROUTINE Partition
!
!##################################################################
!##################################################################
!######                                                      ######
!######              SUBROUTINE WTRADTILTCOL                 ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE wrtgridtilt (rfname,fntimopt,maxelev,nx,ny,nz,kntrelv,radar,   &
                latrad,lonrad,radarx,radary,elvrad,dazim,rngmin,rngmax,   &
                gridtilthighref,gridrangeref,gridslrref,gridazmref,       &
                gridtiltdata,tilttimeref,tilttimevel,itimfrst,elevmeanref,&
                gridtilthighvel,gridrangevel,gridslrvel,gridazmvel,       &
                elevmeanvel,xsc,ysc,zpsc,vcpnum,isource,iver_flg,numvar)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Writes gridded radar data to a file as columns with
!  individual lat,lons for EnKF analysis.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster
!  06/22/95
!
!  MODIFICATION HISTORY:
!
!  10/13/10 Youngsun Jung
!  Adopt wtradcol and merge with wrtgridtilt (see wrtgridtilt_old)
!  for radar data on horizontal arps grid and vertical radar elevation
!  angle (EnKF format).
!
!  06/16/2011 (Y. Wang)
!  Wrapped MPI calls so that it will also be compiled in serial mode
!  without linking MPI library.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    maxelev  Maximum number of tilts
!    kntrelv  Number of tilts
!
!    latrad   : latitude of radar location
!    lonrad   : longitude of  radar location
!    radarx   : x position of radar location
!    radary   : y position of radar location
!    elvrad   : elevation of radar
!    dazim    : average azimuth angle
!    rngmin   Minimum range (m) of data to use
!            (10 000 m or more to eliminate near field ground targets).
!    rngmax   Maximum range (m) of data to use
!
!
!   gridtilthigh(nx,ny,maxelev)     height of  observation point
!   gridrange(nx,ny,maxelev)        range of observation point
!   gridslr(nx,ny)                  surface distance of observation point
!   gridazm(nx,ny)                  azimuth angle
!
!   gridtiltval(nx,ny,maxelev)      observation value
!   gridtilttime(maxelev)           onservation time of each tilt
!
!   elevmean(maxelev)               mean elevation angle
!
!   istatus                         Status indicator
!
!  OUTPUT:
!    data are written to file
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
  INCLUDE 'globcst.inc'
  INCLUDE 'grid.inc'
  INCLUDE 'remap.inc'
  INCLUDE 'mp.inc'

  INTEGER, INTENT(IN) :: maxelev
  INTEGER, INTENT(IN) :: nx
  INTEGER, INTENT(IN) :: ny
  INTEGER, INTENT(IN) :: nz
  INTEGER, INTENT(IN) :: kntrelv

  REAL, INTENT(IN)    :: latrad
  REAL, INTENT(IN)    :: lonrad
  REAL, INTENT(IN)    :: radarx
  REAL, INTENT(IN)    :: radary
  REAL, INTENT(IN)    :: elvrad
  REAL, INTENT(IN)    :: dazim
  REAL, INTENT(IN)    :: rngmin
  REAL, INTENT(IN)    :: rngmax

  REAL, INTENT(IN)   :: gridtilthighref(nx,ny,maxelev)
  REAL, INTENT(IN)   :: gridrangeref(nx,ny,maxelev)
  REAL, INTENT(IN)   :: gridslrref(nx,ny)
  REAL, INTENT(IN)   :: gridazmref(nx,ny)

  REAL, INTENT(INOUT):: gridtiltdata(nx,ny,maxelev,numvar)
  REAL, INTENT(IN)   :: elevmeanref(maxelev)

  REAL, INTENT(IN)   :: gridtilthighvel(nx,ny,maxelev)
  REAL, INTENT(IN)   :: gridrangevel(nx,ny,maxelev)
  REAL, INTENT(IN)   :: gridslrvel(nx,ny)
  REAL, INTENT(IN)   :: gridazmvel(nx,ny)

  REAL, INTENT(IN)   :: elevmeanvel(maxelev)

  INTEGER, INTENT(IN)   :: tilttimeref(maxelev)
  INTEGER, INTENT(IN)   :: tilttimevel(maxelev)

  INTEGER, INTENT(IN) :: numvar

  REAL   :: refsngl(nx,ny)
  REAL   :: velsngl(nx,ny)
  REAL   :: elevref
  REAL   :: elevvel

  CHARACTER (LEN=4) :: radar

  CHARACTER (LEN=*)   :: rfname
  integer :: ilen
  integer :: numscan
  integer :: numelev

  REAL,    ALLOCATABLE :: elevmean_recount_vel(:),        elevmean_recount_ref(:)
  REAL,    ALLOCATABLE :: gridrangevel_recount(:,:,:),    gridrangeref_recount(:,:,:)
  REAL,    ALLOCATABLE :: gridtilthighvel_recount(:,:,:), gridtilthighref_recount(:,:,:)
  REAL,    ALLOCATABLE :: gridtiltdata_recount(:,:,:,:)
  INTEGER, ALLOCATABLE :: tilttime_recount(:)

  INTEGER :: iradfmt

  REAL :: xsc(nx)
  REAL :: ysc(ny)
  REAL :: zpsc(nx,ny,nz)

  INTEGER :: count_elev_vel, count_elev_ref

  REAL :: kcol(nx,ny)

  INTEGER :: fntimopt
  INTEGER :: itimfrst

  INTEGER :: vcpnum
  INTEGER :: isource
  INTEGER :: iver_flg
!
!-----------------------------------------------------------------------
!
!  Radar output descriptors
!
!-----------------------------------------------------------------------
!
  INTEGER, PARAMETER :: mxradvr=10
  INTEGER, PARAMETER :: nradvr=7
  INTEGER :: iradvr(mxradvr)
  DATA iradvr /1,2,3,4,5,6,7,0,0,0/
!
!-----------------------------------------------------------------------
!
!  Radar output thresholds
!
!-----------------------------------------------------------------------
!
  REAL, PARAMETER :: refmin = -5.0, refmax=100., velmin=-200., velmax=200.
  REAL, PARAMETER :: misval = -999.0
!
!-----------------------------------------------------------------------
!
!  Misc local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: iunit
  INTEGER :: i,j,k,klev,kntcol
  INTEGER :: idummy
  INTEGER :: istat,sd_id
  INTEGER :: irngmin,irngmax
  REAL    :: gridlat,gridlon
  INTEGER :: iyr,imon,iday,ihour,imin,isec,time1st
  INTEGER :: samefirstlevel,kntrelv14

  INTEGER(2), ALLOCATABLE :: itmp(:,:) ! Temporary array

  INTEGER :: idp
  CHARACTER :: c1dummy*1, c5dummy*5
!
!-----------------------------------------------------------------------
!
!  Radar output variables
!
!-----------------------------------------------------------------------
!
  INTEGER, PARAMETER :: typelev = 2

  INTEGER :: numradcol
  INTEGER :: outnumelev
  REAL    :: xmin, xmax, ymin, ymax

  REAL,    ALLOCATABLE, TARGET :: colx(:),   coly(:)
  REAL,    ALLOCATABLE, TARGET :: collat(:), collon(:)
  INTEGER, ALLOCATABLE, TARGET :: coli(:),   colj(:)
  REAL,    ALLOCATABLE, TARGET :: elevsfc(:)
  REAL,    ALLOCATABLE, TARGET :: colelev(:)
  INTEGER, ALLOCATABLE, TARGET :: colelevtim(:)

  REAL,    ALLOCATABLE, TARGET :: radcolhgt(:,:), radcolrng(:,:)
  REAL,    ALLOCATABLE, TARGET :: radcoldata(:,:,:)

! For old version
  REAL,    ALLOCATABLE, TARGET :: gridelevlg(:)
  REAL,    ALLOCATABLE, TARGET :: gridhgtlg(:,:,:)
  REAL,    ALLOCATABLE, TARGET :: gridrnglg(:,:,:)
  REAL,    ALLOCATABLE, TARGET :: griddatalg(:,:,:,:)

!-----------------------------------------------------------------------
!
! non-MPI output arrays
!
!-----------------------------------------------------------------------

  REAL,    POINTER :: outelev(:)
  INTEGER, POINTER :: outtim(:)

!-----------------------------------------------------------------------
!
! MPI output arrays
!
!-----------------------------------------------------------------------

  REAL,    POINTER :: outx(:),   outy(:)
  REAL,    POINTER :: outlat(:), outlon(:)
  INTEGER, POINTER :: outi(:),   outj(:)
  REAL,    POINTER :: outsfc(:)

  REAL,    POINTER :: outhgt(:,:), outrng(:,:)
  REAL,    POINTER :: outdata(:,:,:)

  INTEGER, PARAMETER :: ROOT = 0
  INTEGER :: ibgn, iend, jbgn,jend
  INTEGER :: nxlg, nylg
  INTEGER :: numradcol_total

! For old version
  REAL,    POINTER :: outelevmean(:)
  REAL,    POINTER :: outgridtilthgt(:,:,:)
  REAL,    POINTER :: outgridtiltrng(:,:,:)
  REAL,    POINTER :: outgridtiltdata(:,:,:,:)

  INTEGER :: numradcol2d
  INTEGER :: outsum
  INTEGER :: outnumcol(max_proc),   outdisp(max_proc)
  INTEGER :: outnumcol2d(max_proc), outdisp2d(max_proc)

  INTEGER :: istatus

!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  idummy=-999

  ibgn = 1
  iend = nx
  jbgn = 1
  jend = ny
  IF (loc_x > 1)       ibgn = 2
  IF (loc_x < nproc_x) iend = nx-2
  IF (loc_y > 1)       jbgn = 2
  IF (loc_y < nproc_y) jend = ny-2

  nxlg = (nx-3)*nproc_x+3
  nylg = (ny-3)*nproc_y+3

! Adopted from wrtgridtilt modified by Jili Dong
  samefirstlevel=0
  numscan=7

  ! Redefine tilt number,remove redundant scan
  numelev=maxelev
  DO i=1,maxelev
     if(elevmeanref(i) == 0.0) numelev=numelev-1
  END DO

  ALLOCATE( elevmean_recount_vel(maxelev), STAT = istat )
  ALLOCATE( elevmean_recount_ref(maxelev), STAT = istat )
  ALLOCATE( gridrangevel_recount(nx,ny,maxelev), STAT = istat )
  ALLOCATE( gridrangeref_recount(nx,ny,maxelev), STAT = istat )
  ALLOCATE( gridtilthighvel_recount(nx,ny,maxelev), STAT = istat )
  ALLOCATE( gridtilthighref_recount(nx,ny,maxelev), STAT = istat )
  ALLOCATE( gridtiltdata_recount(nx,ny,maxelev,numvar), STAT = istat )
  gridtiltdata_recount = -99999.9

  ALLOCATE( tilttime_recount(maxelev), STAT = istat )

  count_elev_ref=1
  elevmean_recount_ref(count_elev_ref)=elevmeanref(1)
  gridtilthighref_recount(1:nx,1:ny,count_elev_ref) = gridtilthighref(1:nx,1:ny,1)
  gridtiltdata_recount   (1:nx,1:ny,count_elev_ref,:) = gridtiltdata    (1:nx,1:ny,1,:)
  gridrangeref_recount   (1:nx,1:ny,count_elev_ref) = gridrangeref   (1:nx,1:ny,1)

  tilttime_recount(count_elev_ref)=tilttimeref(1)

  count_elev_vel=1
  elevmean_recount_vel(count_elev_vel)=elevmeanvel(1)
  gridtilthighvel_recount(1:nx,1:ny,count_elev_vel) = gridtilthighvel(1:nx,1:ny,1)
  gridrangevel_recount   (1:nx,1:ny,count_elev_vel) = gridrangevel   (1:nx,1:ny,1)

  DO i=2,numelev
     if(elevmeanref(i) >= 0.0 .and. elevmeanref(i) < 90.0) then
       if((elevmeanref(i)-elevmeanref(i-1))> 0.2) then
         count_elev_ref=count_elev_ref+1
         elevmean_recount_ref(count_elev_ref)=elevmeanref(i)
         gridtilthighref_recount(1:nx,1:ny,count_elev_ref)=gridtilthighref(1:nx,1:ny,i)
         DO idp = 2,numvar
           gridtiltdata_recount(1:nx,1:ny,count_elev_ref,idp)=gridtiltdata(1:nx,1:ny,i,idp)
         ENDDO
         gridrangeref_recount(1:nx,1:ny,count_elev_ref)=gridrangeref   (1:nx,1:ny,i)
         tilttime_recount(count_elev_ref)=tilttimeref(i)
       end if
     end if

     if(elevmeanvel(i) >= 0.0 .and. elevmeanvel(i) < 90.0) then
       if((elevmeanvel(i)-elevmeanvel(i-1))> 0.2) then
         count_elev_vel=count_elev_vel+1
         elevmean_recount_vel(count_elev_vel)=elevmeanvel(i)
         gridtilthighvel_recount(1:nx,1:ny,count_elev_vel)=gridtilthighvel(1:nx,1:ny,i)
         gridtiltdata_recount    (1:nx,1:ny,count_elev_vel,1)=gridtiltdata(1:nx,1:ny,i,1)
         gridrangevel_recount   (1:nx,1:ny,count_elev_vel)=gridrangevel(1:nx,1:ny,i)
       endif
     endif

  END DO

  outnumelev = count_elev_vel

  DO k =1, outnumelev
    IF(ABS(elevmean_recount_vel(k)-elevmean_recount_ref(k)) > 0.2) THEN
      WRITE(*,*) ' Reflectivity is not at the same tilt as radial velocity'
      WRITE(*,*) ' Z and Vr tilts are ', elevmean_recount_ref(k), elevmean_recount_vel(k)
      stop 12345
    ENDIF
  ENDDO

  ! End of removing redundant scans
  ! Set clear air missing Z to 0
  ! Set Z and Vr within min range 3km to -99999.9
  DO k =1, outnumelev
    DO j=jbgn,jend
      DO i=ibgn,iend
        if( gridtiltdata_recount(i,j,k,2) < 0.0 .and.                    &
            gridtiltdata_recount(i,j,k,2) > -78000.0 ) then
             gridtiltdata_recount(i,j,k,2)=0.0
        end if
        if( gridtiltdata_recount(i,j,k,1) < -78000.0 ) then
             gridtiltdata_recount(i,j,k,1)= -99999.9
        end if
          if( gridrangeref_recount(i,j,k) < 3000.0 .and.                 &
              gridtilthighref_recount(i,j,k) < 3000.0)                   &
              gridtiltdata_recount(i,j,k,2) = -99999.9
          if( gridrangevel_recount(i,j,k) < 3000.0 .and.                 &
              gridtilthighvel_recount(i,j,k) < 3000.0)                   &
              gridtiltdata_recount(i,j,k,1) = -99999.9
      END DO
    END DO
  END DO

  time1st = itimfrst
  CALL abss2ctim(time1st,iyr,imon,iday,ihour,imin,isec)

  IF(iver_flg == 1) THEN

    ALLOCATE(gridelevlg(outnumelev), STAT = istat)
    ALLOCATE(gridhgtlg (nx,ny,outnumelev), STAT = istat)
    ALLOCATE(gridrnglg (nx,ny,outnumelev), STAT = istat)
    ALLOCATE(griddatalg (nx,ny,outnumelev,numvar), STAT = istat)

    gridelevlg = elevmean_recount_ref(1:outnumelev)
    gridhgtlg  = gridtilthighref_recount(:,:,1:outnumelev)
    gridrnglg  = gridrangeref_recount(:,:,1:outnumelev)
    griddatalg  = gridtiltdata_recount(:,:,1:outnumelev,:)

    IF (mp_opt > 0) THEN
      IF (myproc == ROOT) THEN
        ALLOCATE(outelevmean(outnumelev), STAT = istat)
        ALLOCATE(outgridtilthgt(nxlg,nylg,outnumelev), STAT = istat)
        ALLOCATE(outgridtiltrng(nxlg,nylg,outnumelev), STAT = istat)
        ALLOCATE(outgridtiltdata(nxlg,nylg,outnumelev,numvar), STAT = istat)

      ELSE
        ALLOCATE(outelevmean(1), STAT = istat)
        ALLOCATE(outgridtilthgt(1,1,1), STAT = istat)
        ALLOCATE(outgridtiltrng(1,1,1), STAT = istat)
        ALLOCATE(outgridtiltdata(1,1,1,1), STAT = istat)
      END IF
    ELSE
      outelevmean    => gridelevlg
      outgridtilthgt => gridhgtlg
      outgridtiltrng => gridrnglg
      outgridtiltdata => griddatalg
    END IF

    outelevmean = gridelevlg
    CALL mpimerge3d(gridhgtlg,nx,ny,outnumelev,outgridtilthgt)
    CALL mpimerge3d(gridrnglg,nx,ny,outnumelev,outgridtiltrng)
    CALL mpimerge4d(griddatalg,nx,ny,outnumelev,numvar,outgridtiltdata)

    IF (myproc == ROOT) THEN

      CALL getunit(iunit)

      OPEN(iunit, FILE=trim(rfname), STATUS='unknown', FORM='unformatted')

      WRITE(iunit) time1st,iyr,imon,iday,ihour,imin,isec
      WRITE(iunit) outnumelev,nxlg,nylg
      WRITE(iunit) radar//'      '
      WRITE(iunit) latrad,lonrad,radarx,radary,elvrad
      WRITE(iunit) dazim,rngmin,rngmax
      WRITE(iunit) outelevmean
      WRITE(iunit) outgridtilthgt
      WRITE(iunit) outgridtiltrng
      DO idp = 1, numvar
        WRITE(iunit) outgridtiltdata(:,:,:,idp)
      ENDDO

      CLOSE(iunit)
      CALL retunit(iunit)

    ENDIF

    DEALLOCATE( elevmean_recount_vel, elevmean_recount_ref )
    DEALLOCATE( gridrangevel_recount, gridrangeref_recount )
    DEALLOCATE( gridtilthighvel_recount, gridtilthighref_recount )
    DEALLOCATE( gridtiltdata_recount )

    DEALLOCATE( gridelevlg )
    DEALLOCATE( gridhgtlg, gridrnglg )
    DEALLOCATE( griddatalg )

    IF (mp_opt > 0) THEN
      DEALLOCATE( outelevmean )
      DEALLOCATE( outgridtilthgt,outgridtiltrng )
      DEALLOCATE( outgridtiltdata )
    END IF

    CALL mpbarrier

    RETURN
  ENDIF

!-----------------------------------------------------------------------
!
! Go through the 3D reflectivity arrays to find the right
! number of radar columns and the number vertical levels for
! allocating working arrays.
!
!  xmin, xmax, ymin, ymax  - Valid range for this radar
!
!  numradcol      - Number of radar columns
!  numelev        - Maximun number of vertical observation among all radar columns
!  typelev        - Vertical elevation type
!                   0: unknow
!                   1: ARPS vertical levels (regular grid)
!                   2: Actual radar observation elevation (irregular elevations)
!
!-----------------------------------------------------------------------

    xmin =  9.E37
    ymin =  9.E37
    xmax = -9.E37
    ymax = -9.E37

    kcol(:,:) = 0

    numelev  = 0
    numradcol= 0
    DO j=jbgn,jend
      DO i=ibgn,iend

        klev=1        ! Note 1 not 0
        DO k=1,outnumelev
          IF((gridtiltdata_recount(i,j,k,2)>refmin .AND.                  &
              gridtiltdata_recount(i,j,k,2)<refmax) .OR.                  &
             (gridtiltdata_recount(i,j,k,1)>velmin .AND.                  &
              gridtiltdata_recount(i,j,k,1)<velmax)) THEN
            klev=klev+1
          END IF
        END DO

        IF(klev > 1) THEN  ! Note 1 not 0
          xmin = MIN(xmin,xsc(i))
          xmax = MAX(xmax,xsc(i))
          ymin = MIN(ymin,ysc(j))
          ymax = MAX(ymax,ysc(j))

          kcol(i,j) = klev
          numelev = max(numelev,klev)
          numradcol=numradcol+1
       END IF
      END DO
    END DO

    numradcol_total = numradcol

    CALL mptotali(numradcol_total)

    outdisp(:)   = 0
    outnumcol(:) = 1
    CALL mpigatheri(numradcol,1,outnumcol,nprocs,                        &
                    outdisp,outdisp,ROOT,nprocs,istatus)
    !CALL MPI_Gather(numradcol,1,MPI_INTEGER,outnumcol,1,MPI_INTEGER,    &
    !                ROOT,MPI_COMM_WORLD,imstat)

    CALL mpmaxi(numelev)

    CALL mpmax0(xmax,xmin)
    CALL mpmax0(ymax,ymin)
!    write(0,*) 'myproc = ',myproc, ' numradcol / total = ', numradcol,numradcol_total

!-----------------------------------------------------------------------
!
! Assign output arrays
!
!  radcollat(:)   - radar column latitude
!  radcollon(:)   - radar column longitude
!  radnumelev(:)  - Number of observation in each column
!
!  radcoldata(:,:,:)  - observations in each radar column at each level
!
!-----------------------------------------------------------------------
!
    ALLOCATE(colx(numradcol), STAT = istat)
    ALLOCATE(coly(numradcol), STAT = istat)
    ALLOCATE(coli(numradcol), STAT = istat)
    ALLOCATE(colj(numradcol), STAT = istat)

    ALLOCATE(collat(numradcol), STAT = istat)
    ALLOCATE(collon(numradcol), STAT = istat)

    ALLOCATE(elevsfc(numradcol), STAT = istat)
    ALLOCATE(colelev(outnumelev), STAT = istat)
    ALLOCATE(colelevtim(outnumelev), STAT = istat)

    colx(:) = misval
    coly(:) = misval
    coli(:) = 0
    colj(:) = 0
    elevsfc(:) = misval
    collat(:)  = misval
    collon(:)  = misval
    colelev(:) = 0.0
    colelevtim(:) = 0

    ALLOCATE(radcolhgt(outnumelev,numradcol), STAT = istat)
    ALLOCATE(radcolrng(outnumelev,numradcol), STAT = istat)
    ALLOCATE(radcoldata(outnumelev,numradcol,numvar), STAT = istat)

    radcolhgt(:,:) = misval
    radcolrng(:,:) = misval
    radcoldata(:,:,:) = misval


    DO k=1,outnumelev
      colelev(k) = (elevmean_recount_ref(k)+elevmean_recount_vel(k))/2.
      colelevtim(k) = tilttime_recount(k)
    END DO

    kntcol = 1
    DO j = jbgn, jend
      DO i = ibgn, iend

        IF (kcol(i,j) > 0) THEN
          DO k=1,outnumelev
            IF((gridtiltdata_recount(i,j,k,2)>refmin .AND.             &
                gridtiltdata_recount(i,j,k,2)<refmax) .OR.             &
               (gridtiltdata_recount(i,j,k,1)>velmin .AND.             &
                gridtiltdata_recount(i,j,k,1)<velmax))THEN
              radcolhgt(k,kntcol) = gridtilthighref_recount(i,j,k)
              radcolrng(k,kntcol) = gridrangeref_recount(i,j,k)
              DO idp = 1, numvar
                radcoldata(k,kntcol,idp) = gridtiltdata_recount(i,j,k,idp)
              ENDDO
            END IF
          END DO
!
!-----------------------------------------------------------------------
!
!  If there are data in this column, write them to the file.
!
!-----------------------------------------------------------------------
!

          coli(kntcol) = (loc_x-1)*(nx-3) + i      ! use in binary format only
          colj(kntcol) = (loc_y-1)*(ny-3) + j
          colx(kntcol) = xsc(i)
          coly(kntcol) = ysc(j)

          elevsfc(kntcol)=0.5*(zpsc(i,j,1)+zpsc(i,j,2))  ! use in binary format only

          CALL xytoll(1,1,xsc(i),ysc(j),gridlat,gridlon)
          collat(kntcol) = gridlat
          collon(kntcol) = gridlon

          kntcol=kntcol+1    ! next radcol NO.
        END IF
      END DO
    END DO

!-----------------------------------------------------------------------
!
! Allocate arrays for outputs
!
!-----------------------------------------------------------------------

  IF (numradcol_total <= 0) GOTO 9999   ! Nothing is valid

  IF (mp_opt > 0) THEN
    IF (myproc == ROOT) THEN

      ALLOCATE(outx(numradcol_total), STAT = istat)
      ALLOCATE(outy(numradcol_total), STAT = istat)

      ALLOCATE(outi(numradcol_total), STAT = istat)
      ALLOCATE(outj(numradcol_total), STAT = istat)

      ALLOCATE(outlat(numradcol_total), STAT = istat)
      ALLOCATE(outlon(numradcol_total), STAT = istat)

      ALLOCATE(outsfc(numradcol_total), STAT = istat)

      outx(:) = misval
      outy(:) = misval
      outi(:) = 0
      outj(:) = 0
      outsfc(:) = misval
      outlat(:) = misval
      outlon(:) = misval

      ALLOCATE(outhgt(outnumelev,numradcol_total), STAT = istat)
      ALLOCATE(outrng(outnumelev,numradcol_total), STAT = istat)
      ALLOCATE(outdata(outnumelev,numradcol_total,numvar), STAT = istat)

      outhgt(:,:) = misval
      outrng(:,:) = misval
      outdata(:,:,:) = misval

      ALLOCATE(outelev(outnumelev), STAT = istat)
      ALLOCATE(outtim(outnumelev), STAT = istat)

      outelev(:) = 0
      outtim(:)   = misval

    ELSE

      ALLOCATE(outi(1), STAT = istat)
      ALLOCATE(outj(1), STAT = istat)

      ALLOCATE(outx(1), STAT = istat)
      ALLOCATE(outy(1), STAT = istat)

      ALLOCATE(outlat(1), STAT = istat)
      ALLOCATE(outlon(1), STAT = istat)

      ALLOCATE(outsfc(1), STAT = istat)

      ALLOCATE(outhgt(1,1), STAT = istat)
      ALLOCATE(outrng(1,1), STAT = istat)
      ALLOCATE(outdata(1,1,1), STAT = istat)

      ALLOCATE(outelev(1), STAT = istat)
      ALLOCATE(outtim(1),  STAT = istat)
    END IF

    outsum = 0                 ! For 1D variables
    DO i = 1,nprocs
      outdisp(i) = outsum
      outsum = outsum + outnumcol(i)
    END DO
    CALL mpigatherr(colx,numradcol,outx,numradcol_total,                &
                    outnumcol,outdisp,ROOT,nprocs,istatus)
    CALL mpigatherr(coly,numradcol,outy,numradcol_total,                &
                    outnumcol,outdisp,ROOT,nprocs,istatus)
    CALL mpigatheri(coli,numradcol,outi,numradcol_total,                &
                    outnumcol,outdisp,ROOT,nprocs,istatus)
    CALL mpigatheri(colj,numradcol,outj,numradcol_total,                &
                    outnumcol,outdisp,ROOT,nprocs,istatus)

    CALL mpigatherr(elevsfc,numradcol,outsfc,numradcol_total,           &
                    outnumcol,outdisp,ROOT,nprocs,istatus)
    CALL mpigatherr(collat,numradcol,outlat,numradcol_total,            &
                    outnumcol,outdisp,ROOT,nprocs,istatus)
    CALL mpigatherr(collon,numradcol,outlon,numradcol_total,            &
                    outnumcol,outdisp,ROOT,nprocs,istatus)

!    CALL MPI_Gatherv(colx,numradcol,MPI_REAL,                           &
!                     outx,outnumcol,outdisp,MPI_REAL,                   &
!                     ROOT,MPI_COMM_WORLD,imstat)
!    CALL MPI_Gatherv(coly,numradcol,MPI_REAL,                           &
!                     outy,outnumcol,outdisp,MPI_REAL,                   &
!                     ROOT,MPI_COMM_WORLD,imstat)
!    CALL MPI_Gatherv(coli,numradcol,MPI_INTEGER,                        &
!                     outi,outnumcol,outdisp,MPI_INTEGER,                &
!                     ROOT,MPI_COMM_WORLD,imstat)
!    CALL MPI_Gatherv(colj,numradcol,MPI_INTEGER,                        &
!                     outj,outnumcol,outdisp,MPI_INTEGER,                &
!                     ROOT,MPI_COMM_WORLD,imstat)
!
!    CALL MPI_Gatherv(elevsfc,numradcol,MPI_REAL,                        &
!                     outsfc,outnumcol,outdisp,MPI_REAL,                 &
!                     ROOT,MPI_COMM_WORLD,imstat)
!    CALL MPI_Gatherv(collat,numradcol,MPI_REAL,                         &
!                     outlat,outnumcol,outdisp,                          &
!                     MPI_REAL,ROOT,MPI_COMM_WORLD,imstat)
!    CALL MPI_Gatherv(collon,numradcol,MPI_REAL,                         &
!                     outlon,outnumcol,outdisp,                          &
!                     MPI_REAL,ROOT,MPI_COMM_WORLD,imstat)
!
    outsum = 0                  ! For 2D variables
    outnumcol2d(:) = outnumelev*outnumcol(:)
    numradcol2d    = outnumelev*numradcol
    DO i = 1,nprocs
      outdisp2d(i) = outsum
      outsum = outsum + outnumelev*outnumcol(i)
    END DO

    outelev(:) = colelev(:)
    outtim(:)  = colelevtim(:)

    outsum =  outnumelev*numradcol_total  ! Reuse this tempeorary variable
    CALL mpigatherr(radcolhgt,numradcol2d,outhgt,outsum,                &
                    outnumcol2d,outdisp2d,ROOT,nprocs,istatus)
    CALL mpigatherr(radcolrng,numradcol2d,outrng,outsum,                &
                    outnumcol2d,outdisp2d,ROOT,nprocs,istatus)
    DO idp = 1, numvar
      CALL mpigatherr(radcoldata(:,:,idp),numradcol2d,outdata(:,:,idp), &
                      outsum,outnumcol2d,outdisp2d,ROOT,nprocs,istatus)
    ENDDO
  ELSE
    outx => colx
    outy => coly
    outi => coli
    outj => colj
    outsfc => elevsfc
    outlat => collat
    outlon => collon

    outhgt => radcolhgt
    outrng => radcolrng
    outdata => radcoldata

    outelev => colelev
    outtim  => colelevtim
  END IF

!-----------------------------------------------------------------------
!
! Binary format outputs (be compatible with previous version)
!
!-----------------------------------------------------------------------

  IF (myproc == ROOT) THEN

    CALL getunit(iunit)
    OPEN(iunit,FILE=TRIM(rfname),STATUS='UNKNOWN',FORM='UNFORMATTED')
!
!-----------------------------------------------------------------------
!
!  Write radar description variables
!
!-----------------------------------------------------------------------
!
    WRITE(iunit) radar
    WRITE(iunit) ireftim,time1st,vcpnum,isource,idummy,                &
                 idummy,idummy,nxlg,nylg,outnumelev
!
!-----------------------------------------------------------------------
!
!  Write radar and grid description variables
!  This should provide enough info to verify that the
!  proper grid has been chosen.
!
!-----------------------------------------------------------------------
!
    idummy=0
    iradfmt=1
    strhopt=0
    WRITE(iunit) runname
    WRITE(iunit) iradfmt,strhopt,mapproj,rngmin,rngmax,                 &
                 typelev,numradcol_total,outnumelev,idummy,idummy
    WRITE(iunit) dx,dy,dz,dzmin,ctrlat,ctrlon,trulat1,trulat2,trulon,   &
                 sclfct,latrad,lonrad,elvrad,outelev(1),outelev(outnumelev)
    WRITE(iunit) nradvr,iradvr
    WRITE(iunit) xmin, xmax, ymin, ymax

    WRITE(iunit) iyr,imon,iday,ihour,imin,isec
    WRITE(iunit) radarx,radary,dazim
    WRITE(iunit) outelev
    WRITE(iunit) outtim
!
!-----------------------------------------------------------------------
!
!  If there are data in this column, write them to the file.
!
!-----------------------------------------------------------------------
!
    klev = outnumelev
    DO kntcol = 1, numradcol_total
      WRITE(iunit) outi(kntcol),outj(kntcol),outx(kntcol),outy(kntcol), &
                   outlat(kntcol),outlon(kntcol),outsfc(kntcol),klev

      WRITE(iunit) (outhgt(k,kntcol),k=1,klev)
      WRITE(iunit) (outrng(k,kntcol),k=1,klev)
      DO idp = 1, numvar
        WRITE(iunit) (outdata(k,kntcol,idp),k=1,klev)
      ENDDO
    END DO

    CLOSE(iunit)
    CALL retunit(iunit)
!
!-----------------------------------------------------------------------
!
!  Report on what data were written
!
!-----------------------------------------------------------------------
!
  END IF    ! ROOT processor

  IF (mp_opt > 0) THEN
  DEALLOCATE(outx,   outy)
  DEALLOCATE(outlat, outlon)
  DEALLOCATE(outi,   outj)
  DEALLOCATE(outsfc)

  DEALLOCATE(outhgt, outrng)
  DEALLOCATE(outdata)

  DEALLOCATE(outelev, outtim)
  END IF

  9999 CONTINUE

  IF (myproc == ROOT) THEN
    WRITE(6,'(//a,i4.4,i2.2,i2.2,a1,i2.2,a1,i2.2)')                     &
                  ' Output statistics for time: ',                      &
                      iyr,imon,iday,'_',ihour,':',imin
    WRITE(6,'(2(a,i7,a,/))')                                            &
         ' There were ',numradcol_total,' columns written ',            &
         ' of a total ',(nxlg*nylg),' possible.'

    IF (numradcol_total == 0)  WRITE(6,'(1x,a)')                        &
             'Since there are no columns, no file is created.'
  END IF

  DEALLOCATE(colx,   coly)
  DEALLOCATE(collat, collon)
  DEALLOCATE(coli,   colj)
  DEALLOCATE(elevsfc)
  DEALLOCATE(colelev,colelevtim)

  DEALLOCATE(radcolhgt, radcolrng)
  DEALLOCATE(radcoldata)

  DEALLOCATE(elevmean_recount_vel,elevmean_recount_ref)
  DEALLOCATE(gridrangevel_recount,gridrangeref_recount)
  DEALLOCATE(gridtilthighvel_recount,gridtilthighref_recount)
  DEALLOCATE(gridtiltdata_recount)
  DEALLOCATE(tilttime_recount)

  RETURN
END SUBROUTINE wrtgridtilt
