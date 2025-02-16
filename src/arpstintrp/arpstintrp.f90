!
!##################################################################
!##################################################################
!######                                                      ######
!######                   PROGRAM ARPSTINTRP                 ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

PROGRAM arpstintrp
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  This program interpolates two ARPS history data on grid of the same
!  size to a time inbetween them. The output will be written into a new
!  history dump file.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!
!  2/25/1999. Written based on ARPSTINTRP.
!
!  MODIFICATION HISTORY:
!
!  1999/10/13 (Gene Bassett)
!  Corrected a roundoff error problem for the history dump reference
!  times.  Made the history dump output characterists similar
!  to arpsintrp.
!
!  2001/06/18 (Gene Bassett)
!  Corrected error with absolute time (iabstinit variables).
!
!  2002/03/19 (Keith Brewster)
!  Corrected time calculations for the case when the user chooses
!  to have curtim be relative to initime of the input file.
!
!  1 June 2002 Eric Kemp
!  Soil variable updates.
!
!  02/25/2009 (Y. Wang)
!  Added capability to read split files in MPI mode.
!
!  05/08/2012 (Y. Wang)
!  Added capability to read command line for namelist file.
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
!
!-----------------------------------------------------------------------
!
!  Dimension of the base grid (input data).
!
!-----------------------------------------------------------------------
!
  INTEGER :: nx,ny,nz
  INTEGER :: nzsoil, nstyps
!
!-----------------------------------------------------------------------
!
!  ARPS arrays. The last dimension is for the two sets of variables.
!
!-----------------------------------------------------------------------
!
  REAL, ALLOCATABLE :: u     (:,:,:,:)  ! Total u-velocity (m/s).
  REAL, ALLOCATABLE :: v     (:,:,:,:)  ! Total v-velocity (m/s).
  REAL, ALLOCATABLE :: w     (:,:,:,:)  ! Total w-velocity (m/s).
  REAL, ALLOCATABLE :: ptprt (:,:,:,:)  ! Perturbation potential temperature
                              ! from that of base state atmosphere (Kelvin).
  REAL, ALLOCATABLE :: pprt  (:,:,:,:)  ! Perturbation pressure from that
                              ! of base state atmosphere (Pascal).
  REAL, ALLOCATABLE :: qv    (:,:,:,:)  ! Water vapor mixing ratio (kg/kg).

  REAL, ALLOCATABLE :: qscalar(:,:,:,:,:)

  REAL, ALLOCATABLE :: tke   (:,:,:,:)  ! Turbulent Kinetic Energy ((m/s)**2)
  REAL, ALLOCATABLE :: kmh   (:,:,:,:)  ! Horizontal turb. mixing coef. for
                              ! momentum. ( m**2/s )
  REAL, ALLOCATABLE :: kmv   (:,:,:,:)  ! Vertical turb. mixing coef. for
                              ! momentum. ( m**2/s )

  INTEGER, ALLOCATABLE :: soiltyp (:,:,:)   ! Soil type
  REAL, ALLOCATABLE :: stypfrct(:,:,:)   ! Fraction of soil type
  INTEGER, ALLOCATABLE :: vegtyp  (:,:)          ! Vegetation type
  REAL, ALLOCATABLE :: roufns  (:,:)          ! Surface roughness
  REAL, ALLOCATABLE :: lai     (:,:)          ! Leaf Area Index
  REAL, ALLOCATABLE :: veg     (:,:)          ! Vegetation fraction

  REAL, ALLOCATABLE :: tsoil  (:,:,:,:,:)   ! Deep soil temperature (K)
  REAL, ALLOCATABLE :: qsoil  (:,:,:,:,:)   ! Deep soil temperature (K)
  REAL, ALLOCATABLE :: wetcanp(:,:,:,:)   ! Canopy water amount
  REAL, ALLOCATABLE :: snowdpth(:,:,:)           ! Snow cover

  REAL, ALLOCATABLE :: raing(:,:,:)         ! Grid supersaturation rain
  REAL, ALLOCATABLE :: rainc(:,:,:)         ! Cumulus convective rain
  REAL, ALLOCATABLE :: prcrate(:,:,:,:)     ! precipitation rate (kg/(m**2*s))
                                 ! prcrate(:,:,:,:) = total precip. rate
                                 ! prcrate(:,:,:,:) = grid scale precip. rate
                                 ! prcrate(:,:,:,:) = cumulus precip. rate
                                 ! prcrate(:,:,:,:) = microphysics precip. rate

  REAL, ALLOCATABLE :: radfrc(:,:,:,:)  ! Radiation forcing (K/s)
  REAL, ALLOCATABLE :: radsw (:,:,:)    ! Solar radiation reaching the surface
  REAL, ALLOCATABLE :: rnflx (:,:,:)    ! Net radiation flux absorbed by surface
  REAL, ALLOCATABLE :: radswnet(:,:,:)  ! Net shortwave radiation
  REAL, ALLOCATABLE :: radlwin(:,:,:)   ! Incoming longwave radiation

  REAL, ALLOCATABLE :: usflx (:,:,:)    ! Surface flux of u-momentum (kg/(m*s**2))
  REAL, ALLOCATABLE :: vsflx (:,:,:)    ! Surface flux of v-momentum (kg/(m*s**2))
  REAL, ALLOCATABLE :: ptsflx(:,:,:)    ! Surface heat flux (K*kg/(m*s**2))
  REAL, ALLOCATABLE :: qvsflx(:,:,:)    ! Surface moisture flux (kg/(m**2*s))

  REAL, ALLOCATABLE :: ubar  (:,:,:)   ! Base state u-velocity (m/s).
  REAL, ALLOCATABLE :: vbar  (:,:,:)   ! Base state v-velocity (m/s).
  REAL, ALLOCATABLE :: wbar  (:,:,:)   ! Base state w-velocity (m/s).
  REAL, ALLOCATABLE :: ptbar (:,:,:)   ! Base state potential temperature (K)
  REAL, ALLOCATABLE :: pbar  (:,:,:)   ! Base state pressure (Pascal).
  REAL, ALLOCATABLE :: rhobar(:,:,:)   ! Base state air density (kg/m**3)
  REAL, ALLOCATABLE :: qvbar (:,:,:)   ! Base state water vapor specific humidity
                                       ! (kg/kg,2).
  REAL, ALLOCATABLE :: x     (:)   ! The x-coord. of the physical and
                                   !   computational grid. Defined at u-point.
  REAL, ALLOCATABLE :: y     (:)   ! The y-coord. of the physical and
                                   !   computational grid. Defined at v-point.
  REAL, ALLOCATABLE :: z     (:)   ! The z-coord. of the computational grid.
                                   !   Defined at w-point on the staggered grid.
  REAL, ALLOCATABLE :: zp  (:,:,:) ! The physical height coordinate defined at
                                   !   w-point on the staggered grid.
  REAL, ALLOCATABLE :: zpsoil  (:,:,:)   ! The physical height coordinate defined at

  REAL, ALLOCATABLE :: uprt   (:,:,:,:)    ! Perturbation u-velocity (m/s)
  REAL, ALLOCATABLE :: vprt   (:,:,:,:)    ! Perturbation v-velocity (m/s)
  REAL, ALLOCATABLE :: qvprt  (:,:,:,:)    ! Perturbation water vapor specific humidity (kg/kg)

  REAL, ALLOCATABLE :: tem1  (:,:,:)       ! Temporary array
  REAL, ALLOCATABLE :: tem2  (:,:,:)       ! Temporary array
  REAL, ALLOCATABLE :: tem3  (:,:,:)       ! Temporary array
  REAL, ALLOCATABLE :: tem4  (:,:,:)       ! Temporary array
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i, j, k, l, nq
  REAL :: amin, amax

  CHARACTER (LEN=256) :: basdmpfn
  INTEGER :: lbasdmpf
  CHARACTER (LEN=256) :: ternfn,sfcoutfl,soiloutfl,temchar
  INTEGER :: lternfn,lfn
  INTEGER :: iss,is
  INTEGER :: isub,ksub

  REAL :: zpmax
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
  INCLUDE 'indtflg.inc'
  INCLUDE 'grid.inc'
  INCLUDE 'exbc.inc'
  INCLUDE 'mp.inc'

!-----------------------------------------------------------------------

  INTEGER :: hinfmt,houtfmt, nchin, nchout
  CHARACTER (LEN=256) :: filename
  CHARACTER (LEN=256) :: grdbasfn
  INTEGER :: lenfil,lengbf

  INTEGER :: grdbas
  INTEGER :: ireturn

  REAL :: time
  INTEGER :: gboutcnt, vroutcnt

! DATA    gboutcnt, vroutcnt /0,0/
  DATA    vroutcnt /0/

  INTEGER :: nfilemax
  PARAMETER (nfilemax=2)

  INTEGER :: nouttime
  INTEGER, PARAMETER :: nouttimemax=40

  CHARACTER (LEN=256) :: hisfile(nfilemax)
  INTEGER :: nhisfile,nd, length, lenstr
  CHARACTER (LEN=80) :: timsnd
  CHARACTER (LEN=80) :: new_runname
  INTEGER :: tmstrln

  REAL    :: times(nfilemax), outtime(nouttimemax), alpha, beta
  INTEGER :: iabstinit,iabstinit1
  INTEGER :: ioffset
  INTEGER :: year1,month1,day1,hour1,minute1,second1,ioutabst,          &
             ioutabstinit

  INTEGER :: nxlg, nylg
!
!-----------------------------------------------------------------------
!
!  namelist Declarations:
!
!-----------------------------------------------------------------------
!
  INTEGER :: use_data_t
  CHARACTER (LEN=19) :: initime  ! Real time in form of 'year-mo-dy:hr:mn:ss'

  NAMELIST /message_passing/ nproc_x, nproc_y,                          &
                             nproc_x_in,  nproc_y_in,                   &
                             nproc_x_out, nproc_y_out

  NAMELIST /INPUT/hinfmt,nhisfile,grdbasfn,hisfile

  NAMELIST /output/ runname,use_data_t,initime, nouttime,outtime,       &
            dirname,exbcdmp,exbchdfcompr,hdmpfmt,grbpkbit,hdfcompr,     &
            grdout,basout,varout,mstout,rainout,prcout,iceout,          &
            tkeout, trbout,sfcout,landout,radout,flxout,                &
            qcexout,qrexout,qiexout,qsexout,qhexout,                    &
            totout,filcmprs,sfcdmp,soildmp,ngbrz,zbrdmp,gboutcnt

  INTEGER :: unum         ! unit number for reading in namelist
  CHARACTER(LEN=256)   :: nlfile
  INTEGER :: istatus
  LOGICAL :: iexist
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  mgrid = 1
  nestgrd = 0

  CALL mpinit_proc(0)

  IF(myproc == 0) THEN
    WRITE(6,'(/9(/2x,a)/)')                                             &
     '###############################################################', &
     '###############################################################', &
     '###                                                         ###', &
     '###                Welcome to ARPSTINTRP                    ###', &
     '###                                                         ###', &
     '###############################################################', &
     '###############################################################'

    unum = COMMAND_ARGUMENT_COUNT()
    IF (unum > 0) THEN
      CALL GET_COMMAND_ARGUMENT(1, nlfile, lenstr, istatus )
      IF ( nlfile(1:1) == ' ' .OR. istatus /= 0 ) THEN  ! Use standard input to be backward-compatible
        unum = 5
      ELSE
        INQUIRE(FILE=TRIM(nlfile),EXIST=iexist)
        IF (.NOT. iexist) THEN
          WRITE(6,'(1x,3a)') 'WARNING: namelist file - ',               &
                TRIM(nlfile),' does not exist. Falling back to standard input.'
          unum = 5
        END IF
      END IF
    ELSE
      unum = 5
    END IF

    IF (unum /= 5) THEN
      CALL getunit( unum )
      OPEN(unum,FILE=TRIM(nlfile),STATUS='OLD',FORM='FORMATTED')
      WRITE(*,'(1x,3a,/,1x,a,/)') 'Reading ARPSTINTRP namelist from file - ', &
              TRIM(nlfile),' ... ','========================================'
    ELSE
      WRITE(*,'(2(1x,a,/))') 'Waiting namelist from standard input ... ', &
                             '========================================'
    END IF

  END IF

!
!-----------------------------------------------------------------------
!
!  Set the default parameters
!
!-----------------------------------------------------------------------
!
  hinfmt    = 10
  grdbasfn  = 'X'
  nhisfile  = 1
  hisfile(1) = 'X'
  use_data_t = 1
  initime ='0000-00-00:00:00:00'

  nouttime = 1
  outtime(1) = 0.0

  gboutcnt  = 0

  runname   = 'runname_not_set'

  dirname   = './'
  exbcdmp   = 1
  exbchdfcompr = 0
  hdmpfmt   = 1
  grbpkbit  = 16
  hdfcompr  = 0
  filcmprs  = 0
  basout    = 0
  grdout    = 0
  varout    = 1
  mstout    = 1
  iceout    = 1
  tkeout    = 1
  trbout    = 0
  rainout   = 0
  sfcout    = 0
  radout    = 0
  flxout    = 0
  snowout   = 0
  landout   = 0
  qcexout   = 0
  qrexout   = 0
  qiexout   = 0
  qsexout   = 0
  qhexout   = 0
  sfcdmp    = 1
  soildmp   = 1
  ngbrz     = 5
  zbrdmp    = 10000.0
!
!-----------------------------------------------------------------------
!
!  Read namelist
!
!-----------------------------------------------------------------------
!
  nproc_x_in  = 1
  nproc_y_in  = 1

  nproc_x_out = 1
  nproc_y_out = 1

  IF(myproc == 0) THEN
    READ(unum,message_passing)
    WRITE(6,'(a)')'Namelist message_passing was successfully read.'
  END IF
  CALL mpupdatei(nproc_x,1)
  CALL mpupdatei(nproc_y,1)
!  CALL mpupdatei(readsplit,1)
  CALL mpupdatei(nproc_x_in, 1)
  CALL mpupdatei(nproc_y_in, 1)
  CALL mpupdatei(nproc_x_out,1)
  CALL mpupdatei(nproc_y_out,1)

  CALL mpinit_var

  readsplit(:) = 0
  IF (mp_opt > 0) THEN
    IF (nproc_x_in == 1 .AND. nproc_y_in ==1) THEN
      readsplit(:) = 1
    ELSE
      IF (nproc_x_in /= nproc_x) nproc_x_in = nproc_x      ! still did not support join or split reading
      IF (nproc_y_in /= nproc_y) nproc_y_in = nproc_y
    END IF
  END IF

  nhisfile = 2

  IF(myproc == 0) THEN
    READ(unum,INPUT, END=100)
    WRITE(6,'(/a/)') ' Input control parameters read in are:'

    IF( hinfmt == 5) THEN
      WRITE(6,'(2(2x,a))')                                                &
          'The Savi3D not supported as an INPUT file format',             &
          'Job stopped in ARPSTINTRP.'
      CALL arpsstop('ERROR: unsuportted file format.',1)
    END IF

    WRITE(6,'(1x,a,i3)') ' hinfmt      =', hinfmt
    WRITE(6,'(1x,a,i3)') ' nhisfile    =', nhisfile

    length = LEN( grdbasfn )
    CALL strlnth(  grdbasfn, length )
    WRITE(6,'(1x,a,a)')  ' grdbasfn    =', grdbasfn(1:length)

    DO i=1,nhisfile
      length = LEN( hisfile(i) )
      CALL strlnth(  hisfile(i), length )
      WRITE(6,'(1x,a,i3,a,a)') ' hisfile(',i,')=',hisfile(i)(1:length)
    END DO
  END IF
  CALL mpupdatei(hinfmt,1)
  CALL mpupdatec(grdbasfn,256)
  CALL mpupdatec(hisfile, 256*nhisfile)

!
!-----------------------------------------------------------------------
!
!  Set the control parameters for output:
!
!-----------------------------------------------------------------------
!
  IF (myproc == 0) THEN
    WRITE(6,'(/a/)')                                                      &
        ' Reading in control parameters for the output data files..'

    READ(unum,output,END=100)

    IF ( exbchdfcompr > 4 ) rayklow = -1

    WRITE(6,'(/2x,a,a)')                                                  &
        'The run name to be used for constructing output file name is ',  &
        runname
  END IF
  CALL mpupdatec(runname,80)

  new_runname = runname

  CALL mpupdatei(use_data_t,1)
  CALL mpupdatec(initime,19)

  CALL mpupdatei(nouttime,1)
  CALL mpupdater(outtime,nouttimemax)
  CALL mpupdatei(gboutcnt,1)

  CALL mpupdatec(dirname,256)

  CALL mpupdatei(hdmpfmt,1)
  CALL mpupdatei(grbpkbit,1)
  CALL mpupdatei(hdfcompr,1)

  CALL mpupdatei(filcmprs,1)
  CALL mpupdatei(basout,1)
  CALL mpupdatei(grdout,1)
  CALL mpupdatei(varout,1)
  CALL mpupdatei(mstout,1)
  CALL mpupdatei(iceout,1)
  CALL mpupdatei(tkeout,1)
  CALL mpupdatei(trbout,1)
  CALL mpupdatei(rainout,1)
  CALL mpupdatei(sfcout,1)
  CALL mpupdatei(landout,1)
  CALL mpupdatei(exbcdmp,1)
  CALL mpupdatei(exbchdfcompr,1)

  CALL mpupdatei(qcexout,1)
  CALL mpupdatei(qrexout,1)
  CALL mpupdatei(qiexout,1)
  CALL mpupdatei(qsexout,1)
  CALL mpupdatei(qhexout,1)

  CALL mpupdatei(ngbrz,1)
  CALL mpupdater(zbrdmp,1)

  CALL mpupdatei(sfcdmp,1)
  CALL mpupdatei(soildmp,1)

  totout = 1

  IF (mp_opt > 0 .AND. readsplit(FINDX_H) <= 0) THEN
    CALL gtsplitfn(hisfile(1),1,1,loc_x,loc_y,1,1,                      &
                   0,0,1,lvldbg,filename,ireturn)
  ELSE
    WRITE(filename,'(a)') TRIM(hisfile(1))
  END IF

  IF (myproc == 0) THEN
    CALL get_dims_from_data(hinfmt,trim(filename),                      &
                            nx,ny,nz,nzsoil,nstyps, ireturn)
    IF (readsplit(FINDX_H) > 0) THEN
      nx = (nx-3) / nproc_x + 3
      ny = (ny-3) / nproc_y + 3
    END IF
  END IF
  CALL mpupdatei(nx,1)
  CALL mpupdatei(ny,1)
  CALL mpupdatei(nz,1)
  CALL mpupdatei(nzsoil,1)

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

  CALL mpupdatei(nstyps,1)
  IF (nstyps <= 0) nstyps = 1
  nstyp = nstyps
  IF (nstyps <= 0) nstyps = 1
  nstyp = nstyps

  splitdmp = 0
  IF (hdmpfmt > 100 .OR. soildmp > 100 .OR. exbcdmp > 100) THEN
    splitdmp = hdmpfmt/100
    hdmpfmt  = MOD(hdmpfmt,100)

    splitexbc = exbcdmp/100
    exbcdmp   = MOD(exbcdmp,100)

    splitsoil = soildmp/100
    soildmp   = MOD(soildmp,100)

    IF (nproc_x_out > 1 .OR. nproc_y_out > 1) THEN

      IF (MOD(nproc_x_out,nproc_x) /= 0) THEN
        WRITE(6,'(/,1x,2a,I4,a,I4,/)') 'ERROR: ',                       &
        'wrong size of nproc_x_out = ',nproc_x_out,'.',                 &
        'It must be a multipler of nproc_x = ',nproc_x
        CALL arpsstop('Wrong value of nproc_x_out.',1)
      END IF

      IF (MOD(nproc_y_out,nproc_y) /= 0) THEN
        WRITE(6,'(/,1x,2a,I4,a,I4,/)') 'ERROR: ',                       &
        'wrong size of nproc_y_out = ',nproc_y_out,'.',                 &
        'It must be a multipler of nproc_y = ',nproc_y
        CALL arpsstop('Wrong value of nproc_x_out.',1)
      END IF

      nxlg = (nx-3)*nproc_x + 3
      nylg = (ny-3)*nproc_y + 3
      IF (MOD((nxlg-3),nproc_x_out) /= 0) THEN
        WRITE(6,'(/,1x,2a,I4,a,/,8x,a,I5,a,I4,a,/)') 'ERROR: ',         &
        'wrong size of nproc_x_out = ',nproc_x_out,'.',                 &
        'The grid size is ',nxlg,                                       &
        ' and it''s physical size is not dividable by ',nproc_x_out,'.'
        CALL arpsstop('Wrong value of nproc_x_out.',1)
      END IF

      IF (MOD((nylg-3),nproc_y_out) /= 0) THEN
        WRITE(6,'(/,1x,2a,I4,a,/,8x,a,I5,a,I4,a,/)') 'ERROR: ',         &
        'wrong size of nproc_y_out = ',nproc_y_out,'.',                 &
        'The grid size is ',nylg,                                       &
        ' and it''s physical size is not dividable by ',nproc_y_out,'.'
        CALL arpsstop('Wrong value of nproc_y_out.',1)
      END IF

    END IF
  END IF

  joindmp(:) = 0
  IF (mp_opt > 0 .AND. nproc_x_out == 1 .AND. nproc_y_out == 1) THEN
    joindmp(:) = 1
  END IF

  IF (nouttime == 1) THEN
    isub = 1
    ksub = 2
  ELSE
    isub = 3
    ksub = 3
  END IF

  IF (unum /= 5 .AND. myproc == 0) THEN
    CLOSE( unum )
    CALL retunit( unum )
  END IF
!
!-----------------------------------------------------------------------
!
!  Allocate the arrays.
!
!-----------------------------------------------------------------------
!
  ALLOCATE(u(nx,ny,nz,ksub))
  ALLOCATE(v(nx,ny,nz,ksub))
  ALLOCATE(w(nx,ny,nz,ksub))
  ALLOCATE(ptprt(nx,ny,nz,ksub))
  ALLOCATE(pprt(nx,ny,nz,ksub))
  ALLOCATE(qv(nx,ny,nz,ksub))
  ALLOCATE(qscalar(nx,ny,nz,nscalar,ksub))
  ALLOCATE(tke(nx,ny,nz,ksub))
  ALLOCATE(kmh(nx,ny,nz,ksub))
  ALLOCATE(kmv(nx,ny,nz,ksub))

  ALLOCATE(soiltyp(nx,ny,nstyps))
  ALLOCATE(stypfrct(nx,ny,nstyps))
  ALLOCATE(vegtyp(nx,ny))
  ALLOCATE(roufns(nx,ny))
  ALLOCATE(lai(nx,ny))
  ALLOCATE(veg(nx,ny))

  ALLOCATE(tsoil(nx,ny,nzsoil,0:nstyps,ksub))
  ALLOCATE(qsoil(nx,ny,nzsoil,0:nstyps,ksub))
  ALLOCATE(wetcanp(nx,ny,0:nstyps,ksub))
  ALLOCATE(snowdpth(nx,ny,ksub))

  ALLOCATE(raing(nx,ny,ksub))
  ALLOCATE(rainc(nx,ny,ksub))
  ALLOCATE(prcrate(nx,ny,4,ksub))

  ALLOCATE(radfrc(nx,ny,nz,ksub))
  ALLOCATE(radsw(nx,ny,ksub))
  ALLOCATE(rnflx(nx,ny,ksub))
  ALLOCATE(radswnet(nx,ny,ksub))
  ALLOCATE(radlwin(nx,ny,ksub))

  ALLOCATE(usflx(nx,ny,ksub))
  ALLOCATE(vsflx(nx,ny,ksub))
  ALLOCATE(ptsflx(nx,ny,ksub))
  ALLOCATE(qvsflx(nx,ny,ksub))

  ALLOCATE(ubar(nx,ny,nz))
  ALLOCATE(vbar(nx,ny,nz))
  ALLOCATE(wbar(nx,ny,nz))
  ALLOCATE(ptbar(nx,ny,nz))
  ALLOCATE(pbar(nx,ny,nz))
  ALLOCATE(rhobar(nx,ny,nz))
  ALLOCATE(qvbar(nx,ny,nz))
  ALLOCATE(x(nx))
  ALLOCATE(y(ny))
  ALLOCATE(z(nz))
  ALLOCATE(zp(nx,ny,nz))
  ALLOCATE(zpsoil(nx,ny,nzsoil))

  ALLOCATE(uprt(nx,ny,nz,ksub))
  ALLOCATE(vprt(nx,ny,nz,ksub))
  ALLOCATE(qvprt(nx,ny,nz,ksub))

  ALLOCATE(tem1(nx,ny,nz))
  ALLOCATE(tem2(nx,ny,nz))
  ALLOCATE(tem3(nx,ny,nz))
  ALLOCATE(tem4(nx,ny,nz))

  u     (:,:,:,:) = 0.0
  v     (:,:,:,:) = 0.0
  uprt  (:,:,:,:) = 0.0
  vprt  (:,:,:,:) = 0.0
  w     (:,:,:,:) = 0.0
  ptprt (:,:,:,:) = 0.0
  pprt  (:,:,:,:) = 0.0
  qvprt (:,:,:,:) = 0.0
  qscalar(:,:,:,:,:) = 0.0
  tke   (:,:,:,:) = 0.0
  kmh   (:,:,:,:) = 0.0
  kmv   (:,:,:,:) = 0.0
  radfrc(:,:,:,:) = 0.0

  ubar  (:,:,:)    = 0.0
  vbar  (:,:,:)    = 0.0
  wbar  (:,:,:)    = 0.0
  ptbar (:,:,:)    = 0.0
  pbar  (:,:,:)    = 0.0
  rhobar(:,:,:)    = 0.0
  qvbar (:,:,:)    = 0.0

  x     (:)        = 0.0
  y     (:)        = 0.0
  z     (:)        = 0.0
  zp    (:,:,:)    = 0.0

  zpsoil(:,:,:)    = 0.0

  tsoil(:,:,:,:,:) = 0.0
  qsoil(:,:,:,:,:) = 0.0

  wetcanp(:,:,:,:) = 0.0

  soiltyp (:,:,:) = 0
  stypfrct(:,:,:) = 0.0

  snowdpth(:,:,:) = 0.0

  vegtyp (:,:) = 0
  lai    (:,:) = 0.0
  roufns (:,:) = 0.0
  veg    (:,:) = 0.0

  raing  (:,:,:)   = 0.0
  rainc  (:,:,:)   = 0.0
  prcrate(:,:,:,:) = 0.0
  radsw  (:,:,:)   = 0.0
  rnflx  (:,:,:)   = 0.0
  radswnet(:,:,:)  = 0.0
  radlwin(:,:,:)   = 0.0
  usflx  (:,:,:)   = 0.0
  vsflx  (:,:,:)   = 0.0
  ptsflx (:,:,:)   = 0.0
  qvsflx (:,:,:)   = 0.0


  ldirnam=LEN(dirname)
  CALL strlnth( dirname , ldirnam)

  lengbf=LEN(grdbasfn)
  CALL strlnth( grdbasfn, lengbf)

  IF (mp_opt > 0 .AND. readsplit(FINDX_H) <= 0) THEN
    CALL gtsplitfn(grdbasfn,1,1,loc_x,loc_y,1,1,                      &
                   0,0,1,lvldbg,filename,ireturn)
    lengbf = LEN_TRIM(filename)
    grdbasfn(1:lengbf) = TRIM(filename)
  END IF

  IF (myproc == 0) WRITE(6,'(/a,a)')' The grid/base name is ', grdbasfn(1:lengbf)

!
!-----------------------------------------------------------------------
!
!  Loop over data files
!
!-----------------------------------------------------------------------
!
  ireturn = 0

  DO nd=1,2

    IF (mp_opt > 0 .AND. readsplit(FINDX_H) <= 0) THEN
      CALL gtsplitfn(hisfile(nd),1,1,loc_x,loc_y,1,1,                   &
                     0,0,1,lvldbg,filename,ireturn)
    ELSE
      filename = hisfile(nd)
    END IF

    lenfil=LEN(filename)
    CALL strlnth( filename, lenfil)
    IF (myproc == 0) WRITE(6,'(/a,a,a)')                                &
        ' Data set ', filename(1:lenfil) ,' to be processed.'
!
!-----------------------------------------------------------------------
!
!  Read all input data arrays
!
!-----------------------------------------------------------------------
!
    CALL dtaread(nx,ny,nz,nzsoil,nstyps,                                &
        hinfmt,nchin,grdbasfn(1:lengbf),lengbf,                         &
        filename(1:lenfil),lenfil,time, x,y,z,zp,zpsoil,                &
        uprt(1,1,1,nd),vprt (1,1,1,nd),w  (1,1,1,nd),ptprt(1,1,1,nd),   &
        pprt(1,1,1,nd),qvprt(1,1,1,nd),qscalar(1,1,1,1,nd),             &
        tke  (1,1,1,nd),                                                &
        kmh (1,1,1,nd),kmv  (1,1,1,nd),                                 &
        ubar,vbar,wbar,ptbar,pbar,rhobar,qvbar,                         &
        soiltyp,stypfrct,vegtyp,lai,roufns,veg,                         &
        tsoil(1,1,1,0,nd),qsoil(1,1,1,0,nd),                            &
        wetcanp(1,1,0,nd),snowdpth(1,1,nd),                             &
        raing(1,1,nd),rainc(1,1,nd),prcrate(1,1,1,nd),                  &
        radfrc(1,1,1,nd),radsw(1,1,nd),rnflx(1,1,nd),                   &
        radswnet(1,1,nd),radlwin(1,1,nd),                               &
        usflx(1,1,nd),vsflx(1,1,nd),ptsflx(1,1,nd),qvsflx(1,1,nd),      &
        ireturn, tem1, tem2, tem3)

    IF (myproc == 0) WRITE(6,'(/a,/2(a,i2),a,i4,a,/3(a,i2),a,f13.3,a/)')&
        'History data read in for time: ',                              &
        'month=',month,',   day=',   day,',   year=',year,',',          &
        'hour =',hour ,',minute=',minute,', second=',second,            &
        ', time=',time,'(s)'

    CALL ctim2abss(year,month,day,hour,minute,second,iabstinit)

    IF(nd == 1) THEN  ! Save the values for data set 1
      year1  = year
      month1 = month
      day1   = day
      hour1  = hour
      minute1= minute
      second1= second
      iabstinit1 = iabstinit
    END IF

    times(nd) = time + int(iabstinit - iabstinit1)

  END DO

  IF( hinfmt == 9 .AND. ireturn == 2 ) THEN
    WRITE(6,'(/1x,a/)') 'The end of GrADS file was reached.'
    CLOSE ( nchin )
    CALL retunit( nchin )
    GO TO 9001
  END IF

  IF( ireturn /= 0 ) GO TO 9002            ! Read was unsuccessful
!

  IF( use_data_t == 1) THEN ! Use init time in input file
    ioutabstinit=iabstinit1
    year  = year1
    month = month1
    day   = day1
    hour  = hour1
    minute= minute1
    second= second1
  ELSE
    READ(initime,'(i4.4,1x,i2.2,1x,i2.2,1x,i2.2,1x,i2.2,1x,i2.2)')      &
                   year,month,day,hour,minute,second
    CALL ctim2abss(year,month,day,hour,minute,second,ioutabstinit)
  END IF

  !PRINT*,'ioutabstinit=',ioutabstinit
  !PRINT*,'iabstinit1=',iabstinit1

  ioffset = ioutabstinit - iabstinit1
  times(1) = times(1) - ioffset
  times(2) = times(2) - ioffset

  DO l=1,nouttime
    curtim = outtime(l)

    IF (myproc == 0) WRITE(6,'(/a,/2(a,i2),a,i4,a,/3(a,i2),a,f13.3,a/)')&
        'In output file, the reference time is',                        &
        'month=',month,',   day=',   day,',   year=',year,',',          &
        'hour =',hour ,',minute=',minute,', second=',second,            &
        ', & the time relative to this reference =',curtim

    IF ( curtim > MAX(times(1),times(2)) .OR.                           &
         curtim < MIN(times(1),times(2))) THEN
      WRITE (*,*) "WARNING: Performing extrapolation.  Desired time ",  &
                  "is outside the range of the reference files."
    END IF

    IF( times(2) == times(1)) THEN
      WRITE (*,*) "ERROR: times in reference files are the same, ",     &
                  "can't perform interpolation."
      CALL arpsstop('No interpolation',1)
    END IF

    alpha = (times(2)-curtim)/(times(2)-times(1))
    beta = 1.0-alpha

    IF (myproc == 0) WRITE (6,'(/,1x,2(a,F12.7))')                      &
         'Relative weights: file 1 - ',alpha,', file 2 - ',beta

!-----------------------------------------------------------------------
!
!  Calculate total fields from that for base state and perturbations
!
!-----------------------------------------------------------------------
!
    DO k=1,nz
      DO j=1,ny
        DO i=1,nx
          uprt  (i,j,k,isub)=alpha*uprt  (i,j,k,1)+beta*uprt  (i,j,k,2)
          vprt  (i,j,k,isub)=alpha*vprt  (i,j,k,1)+beta*vprt  (i,j,k,2)
          w     (i,j,k,isub)=alpha*w     (i,j,k,1)+beta*w     (i,j,k,2)
          ptprt (i,j,k,isub)=alpha*ptprt (i,j,k,1)+beta*ptprt (i,j,k,2)
          pprt  (i,j,k,isub)=alpha*pprt  (i,j,k,1)+beta*pprt  (i,j,k,2)
          qvprt (i,j,k,isub)=alpha*qvprt (i,j,k,1)+beta*qvprt (i,j,k,2)
          DO nq=1,nscalar
            qscalar(i,j,k,nq,isub)=alpha*qscalar(i,j,k,nq,1)+beta*qscalar(i,j,k,nq,2)
          END DO
          tke   (i,j,k,isub)=alpha*tke   (i,j,k,1)+beta*tke   (i,j,k,2)
          kmh   (i,j,k,isub)=alpha*kmh   (i,j,k,1)+beta*kmh   (i,j,k,2)
          kmv   (i,j,k,isub)=alpha*kmv   (i,j,k,1)+beta*kmv   (i,j,k,2)
          u     (i,j,k,isub)=uprt (i,j,k,isub)+ubar (i,j,k)
          v     (i,j,k,isub)=vprt (i,j,k,isub)+vbar (i,j,k)
          qv    (i,j,k,isub)=qvprt(i,j,k,isub)+qvbar(i,j,k)
          radfrc(i,j,k,isub)=alpha*radfrc(i,j,k,1)+beta*radfrc(i,j,k,2)
        END DO
      END DO
    END DO

    DO is=0,nstyp
      DO k=1,nzsoil
        DO j=1,ny-1
          DO i=1,nx-1
          tsoil(i,j,k,is,isub)=alpha*tsoil(i,j,k,is,1)                  &
                           +beta*tsoil(i,j,k,is,2)
          qsoil(i,j,k,is,isub)=alpha*qsoil(i,j,k,is,1)                  &
                           +beta*qsoil(i,j,k,is,2)
          END DO
        END DO
      END DO
      DO j=1,ny-1
        DO i=1,nx-1
          wetcanp(i,j,is,isub)=alpha*wetcanp(i,j,is,1)                  &
                            +beta*wetcanp(i,j,is,2)
        END DO
      END DO
    END DO

    DO j=1,ny-1
      DO i=1,nx-1
        snowdpth(i,j,isub)=alpha*snowdpth(i,j,1)+beta*snowdpth(i,j,2)
        raing  (i,j,isub)=alpha*raing  (i,j,1)+beta*raing  (i,j,2)
        rainc  (i,j,isub)=alpha*rainc  (i,j,1)+beta*rainc  (i,j,2)

        prcrate(i,j,1,isub)=alpha*prcrate(i,j,1,1)+beta*prcrate(i,j,1,2)
        prcrate(i,j,2,isub)=alpha*prcrate(i,j,2,1)+beta*prcrate(i,j,2,2)
        prcrate(i,j,3,isub)=alpha*prcrate(i,j,3,1)+beta*prcrate(i,j,3,2)
        prcrate(i,j,4,isub)=alpha*prcrate(i,j,4,1)+beta*prcrate(i,j,4,2)

        radsw  (i,j,isub)=alpha*radsw  (i,j,1)+beta*radsw  (i,j,2)
        rnflx  (i,j,isub)=alpha*rnflx  (i,j,1)+beta*rnflx  (i,j,2)
        radswnet(i,j,isub)=alpha*radswnet(i,j,1)+beta*radswnet(i,j,2)
        radlwin(i,j,isub)=alpha*radlwin(i,j,1)+beta*radlwin(i,j,2)
        usflx  (i,j,isub)=alpha*usflx  (i,j,1)+beta*usflx  (i,j,2)
        vsflx  (i,j,isub)=alpha*vsflx  (i,j,1)+beta*vsflx  (i,j,2)
        ptsflx (i,j,isub)=alpha*ptsflx (i,j,1)+beta*ptsflx (i,j,2)
        qvsflx (i,j,isub)=alpha*qvsflx (i,j,1)+beta*qvsflx (i,j,2)
      END DO
    END DO

!
!-----------------------------------------------------------------------
!
!  Print out the max/min of output variables.
!
!-----------------------------------------------------------------------
!
    IF (myproc == 0) WRITE(6,'(/1x,a/)')                                &
        'Min. and max. of data interpolated to the new time:'

    CALL a3dmax0(x,1,nx,1,nx,1,1,1,1,1,1,1,1, amax,amin)
    IF (myproc == 0) WRITE(6,'(/1x,2(a,e13.6))') 'xmin    = ', amin,',  xmax    =',amax

    CALL a3dmax0(y,1,ny,1,ny,1,1,1,1,1,1,1,1, amax,amin)
    IF (myproc == 0) WRITE(6,'(1x,2(a,e13.6))') 'ymin    = ', amin,',  ymax    =',amax

    CALL a3dmax0(z,1,nz,1,nz,1,1,1,1,1,1,1,1, amax,amin)
    IF (myproc == 0) WRITE(6,'(1x,2(a,e13.6))') 'zmin    = ', amin,',  zmax    =',amax

    CALL a3dmax0(zp,1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz,                  &
                amax,amin)
    IF (myproc == 0) WRITE(6,'(1x,2(a,e13.6))') 'zpmin   = ', amin,', zpmax    =',amax

    CALL a3dmax0(ubar,1,nx,1,nx,1,ny,1,ny-1,1,nz,1,nz-1,                &
                amax,amin)
    IF (myproc == 0) WRITE(6,'(1x,2(a,e13.6))') 'ubarmin = ', amin,',  ubarmax =',amax

    CALL a3dmax0(vbar,1,nx,1,nx-1,1,ny,1,ny,1,nz,1,nz-1,                &
                amax,amin)
    IF (myproc == 0) WRITE(6,'(1x,2(a,e13.6))') 'vbarmin = ', amin,',  vbarmax =',amax

    CALL a3dmax0(ptbar,1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz-1,             &
                amax,amin)
    IF (myproc == 0) WRITE(6,'(1x,2(a,e13.6))') 'ptbarmin= ', amin,',  ptbarmax=',amax

    CALL a3dmax0(pbar,1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz-1,              &
                amax,amin)
    IF (myproc == 0) WRITE(6,'(1x,2(a,e13.6))') 'pbarmin = ', amin,',  pbarmax =',amax

    CALL a3dmax0(rhobar,1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz-1,            &
                amax,amin)
    IF (myproc == 0) WRITE(6,'(1x,2(a,e13.6))') 'rhobarmin=', amin,', rhobarmax=',amax

    CALL a3dmax0(qvbar,1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz-1,             &
                amax,amin)
    IF (myproc == 0) WRITE(6,'(1x,2(a,e13.6))') 'qvbarmin= ', amin,',  qvbarmax=',amax

    CALL a3dmax0(uprt(1,1,1,isub),1,nx,1,nx,1,ny,1,ny-1,1,nz,1,nz-1,    &
                amax,amin)
    IF (myproc == 0) WRITE(6,'(1x,2(a,e13.6))') 'uprtmin = ', amin,',  uprtmax =',amax

    CALL a3dmax0(vprt(1,1,1,isub),1,nx,1,nx-1,1,ny,1,ny,1,nz,1,nz-1,    &
                amax,amin)
    IF (myproc == 0) WRITE(6,'(1x,2(a,e13.6))') 'vprtmin = ', amin,',  vprtmax =',amax

    CALL a3dmax0(w(1,1,1,isub),1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz,       &
                amax,amin)
    IF (myproc == 0) WRITE(6,'(1x,2(a,e13.6))') 'wmin    = ', amin,',  wmax    =',amax

    CALL a3dmax0(ptprt(1,1,1,isub),1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz-1, &
                amax,amin)
    IF (myproc == 0) WRITE(6,'(1x,2(a,e13.6))') 'ptprtmin= ', amin,',  ptprtmax=',amax

    CALL a3dmax0(pprt(1,1,1,isub),1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz-1,  &
                amax,amin)
    IF (myproc == 0) WRITE(6,'(1x,2(a,e13.6))') 'pprtmin = ', amin,',  pprtmax =',amax

    CALL a3dmax0(qvprt(1,1,1,isub),1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz-1, &
                amax,amin)
    IF (myproc == 0) WRITE(6,'(1x,2(a,e13.6))') 'qvprtmin= ', amin,',  qvprtmax=',amax

    DO nq = 1,nscalar

      CALL a3dmax0(qscalar(:,:,:,nq,isub),1,nx,1,nx-1,1,ny,1,ny-1,         &
                   1,nz,1,nz-1,amax,amin)

      IF (myproc == 0) WRITE(6,'(1x,2(a,e13.6))') TRIM(qnames(nq))//'min   = ', amin,    &
                            ',  '//TRIM(qnames(nq))//'max   =',  amax

    END DO

    CALL a3dmax0(tke(1,1,1,isub),1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz-1,   &
                amax,amin)
    IF (myproc == 0) WRITE(6,'(1x,2(a,e13.6))') 'tkemin  = ', amin,',  tkemax  =',amax

    CALL a3dmax0(kmh(1,1,1,isub),1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz-1,   &
                amax,amin)
    IF (myproc == 0) WRITE(6,'(1x,2(a,e13.6))') 'kmhmin  = ', amin,',  kmhmax  =',amax

    CALL a3dmax0(kmv(1,1,1,isub),1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz-1,   &
                amax,amin)
    IF (myproc == 0) WRITE(6,'(1x,2(a,e13.6))') 'kmvmin  = ', amin,',  kmvmax  =',amax

    CALL a3dmax0(raing(1,1,isub),1,nx,1,nx-1,1,ny,1,ny-1,1,1,1,1,amax,amin)
    IF (myproc == 0) WRITE(6,'(1x,2(a,e13.6))') 'raingmin= ', amin,',  raingmax=',amax

    CALL a3dmax0(rainc(1,1,isub),1,nx,1,nx-1,1,ny,1,ny-1,1,1,1,1,amax,amin)
    IF (myproc == 0) WRITE(6,'(1x,2(a,e13.6))') 'raincmin= ', amin,',  raincmax=',amax

    CALL a3dmax0(prcrate(1,1,1,isub),1,nx,1,nx-1,1,ny,1,ny-1,           &
                 1,1,1,1,amax,amin)
    IF (myproc == 0) WRITE(6,'(1x,2(a,e13.6))') 'prcr1min= ', amin,',  prcr1max=',amax

    CALL a3dmax0(prcrate(1,1,2,isub),1,nx,1,nx-1,1,ny,1,ny-1,           &
                 1,1,1,1,amax,amin)
    IF (myproc == 0) WRITE(6,'(1x,2(a,e13.6))') 'prcr2min= ', amin,',  prcr2max=',amax

    CALL a3dmax0(prcrate(1,1,3,isub),1,nx,1,nx-1,1,ny,1,ny-1,           &
                 1,1,1,1,amax,amin)
    IF (myproc == 0) WRITE(6,'(1x,2(a,e13.6))') 'prcr3min= ', amin,',  prcr3max=',amax

    CALL a3dmax0(prcrate(1,1,4,isub),1,nx,1,nx-1,1,ny,1,ny-1,           &
                 1,1,1,1,amax,amin)
    IF (myproc == 0) WRITE(6,'(1x,2(a,e13.6))') 'prcr4min= ', amin,',  prcr4max=',amax

    DO iss = 0, nstyp

      CALL a3dmax0(tsoil(1,1,1,iss,isub),1,nx,1,nx-1,1,ny,1,ny-1,1,1,1,1, &
                   amax,amin)
      IF (myproc == 0) WRITE(6,'(1x,2(a,e13.6),a,i3)')                    &
          'tsoil_sfcmin = ', amin,',  tsoil_sfcmax =',amax,' for soil type=',iss

      CALL a3dmax0(tsoil(1,1,2,iss,isub),1,nx,1,nx-1,1,ny,1,ny-1,1,1,1,1, &
                   amax,amin)
      IF (myproc == 0) WRITE(6,'(1x,2(a,e13.6),a,i3)')                    &
          'tsoil_dpmin= ', amin,',  tsoil_dpmax=',amax,' for soil type=',iss

      CALL a3dmax0(qsoil(1,1,1,iss,isub),1,nx,1,nx-1,1,ny,1,ny-1,1,1,1,1, &
                   amax,amin)
      IF (myproc == 0) WRITE(6,'(1x,2(a,e13.6),a,i3)')                    &
          'qsoil_sfcmin = ', amin,',  qsoil_sfcmax =',amax,' for soil type=',iss

      CALL a3dmax0(qsoil(1,1,2,iss,isub),1,nx,1,nx-1,1,ny,1,ny-1,1,1,1,1, &
                   amax,amin)
      IF (myproc == 0) WRITE(6,'(1x,2(a,e13.6),a,i3)')                    &
          'qsoil_dpmin = ', amin,',  qsoil_dpmax =',amax,' for soil type=',iss

      CALL a3dmax0(wetcanp(1,1,iss,isub),                                 &
                   1,nx,1,nx-1,1,ny,1,ny-1,1,1,1,1,amax,amin)
      IF (myproc == 0) WRITE(6,'(1x,2(a,e13.6),a,i3)')                    &
          'wetcmin = ', amin,',  wetcmax =',amax,' for soil type=',iss

    END DO

    CALL a3dmax0(roufns,1,nx,1,nx-1,1,ny,1,ny-1,1,1,1,1,                  &
                 amax,amin)
    IF (myproc == 0) WRITE(6,'(1x,2(a,e13.6))') 'roufnmin =', amin,', roufnmax =',amax

    CALL a3dmax0(veg,1,nx,1,nx-1,1,ny,1,ny-1,1,1,1,1,                     &
                 amax,amin)
    IF (myproc == 0) WRITE(6,'(1x,2(a,e13.6))') 'vegmin  = ', amin,',   vegmax =',amax

    CALL a3dmax0(radfrc,1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz-1,              &
                 amax,amin)
    IF (myproc == 0) WRITE(6,'(1x,2(a,e13.6))') 'radfnmin =', amin,', radfnmax =',amax

    CALL a3dmax0(radsw(1,1,isub),1,nx,1,nx-1,1,ny,1,ny-1,1,1,1,1,         &
                 amax,amin)
    IF (myproc == 0) WRITE(6,'(1x,2(a,e13.6))') 'radswmin =', amin,', radswmax =',amax

    CALL a3dmax0(rnflx(1,1,isub),1,nx,1,nx-1,1,ny,1,ny-1,1,1,1,1,         &
                 amax,amin)
    IF (myproc == 0) WRITE(6,'(1x,2(a,e13.6))') 'rnflxmin =', amin,', rnflxmax =',amax

    CALL a3dmax0(radswnet(1,1,isub),1,nx,1,nx-1,1,ny,1,ny-1,1,1,1,1,      &
                 amax,amin)
    IF (myproc == 0) WRITE(6,'(1x,2(a,e13.6))') 'radswnetmin =', amin,', radswnetmax =',amax

    CALL a3dmax0(radlwin(1,1,isub),1,nx,1,nx-1,1,ny,1,ny-1,1,1,1,1,       &
                 amax,amin)
    IF (myproc == 0) WRITE(6,'(1x,2(a,e13.6))') 'radlwinmin =', amin,', radlwinmax =',amax

    CALL a3dmax0(usflx(1,1,isub),1,nx,1,nx-1,1,ny,1,ny-1,1,1,1,1,         &
                 amax,amin)
    IF (myproc == 0) WRITE(6,'(1x,2(a,e13.6))') 'usflxnmin =', amin,', usflxmax =',amax

    CALL a3dmax0(vsflx(1,1,isub),1,nx,1,nx-1,1,ny,1,ny-1,1,1,1,1,         &
                 amax,amin)
    IF (myproc == 0) WRITE(6,'(1x,2(a,e13.6))') 'vsflxmin =', amin,', vsflxmax =',amax

    CALL a3dmax0(ptsflx(1,1,isub),1,nx,1,nx-1,1,ny,1,ny-1,1,1,1,1,        &
                 amax,amin)
    IF (myproc == 0) WRITE(6,'(1x,2(a,e13.6))') 'ptflxmin =', amin,', ptflxmax =',amax

    CALL a3dmax0(qvsflx(1,1,isub),1,nx,1,nx-1,1,ny,1,ny-1,1,1,1,1,        &
                 amax,amin)
    IF (myproc == 0) WRITE(6,'(1x,2(a,e13.6))') 'qvflxmin =', amin,', qvflxmax =',amax

    DO k=1,nz
      DO j=1,ny
        DO i=1,nx
          tem1 (i,j,k)= 0.0        ! To be put in place of wbar
        END DO
      END DO
    END DO
!
!-----------------------------------------------------------------------
!
!  Data dump of the model grid and base state arrays:
!
!  First find a unique name basdmpfn(1:lbasdmpf) for the grid and
!  base state array dump file
!
!  If grid/base state data has been written out once, skip
!  the following writing block. Also no need to write out
!  separate data for Savi3D dump. The same for GrADS dump.
!
!-----------------------------------------------------------------------
!
    CALL mpupdatei(nocmnt,1)
    IF (mp_opt > 0 .AND. joindmp(FINDX_H) > 0) THEN
      WRITE (cmnt(nocmnt),'(a,i4,a,i4,a,i4)')                           &
          ' nx =',(nx-3)*nproc_x+3,', ny =',(ny-3)*nproc_y+3,', nz =',nz
    ELSE
      WRITE (cmnt(nocmnt),'(a,i4,a,i4,a,i4)')                           &
          ' nx =',nx,', ny =',ny,', nz =',nz
    END IF

    runname = new_runname
    houtfmt = hdmpfmt
    grbpkbit = 16

    CALL gtlfnkey(runname, lfnkey)

    IF(houtfmt /= 9 ) THEN

      IF( gboutcnt == 1 ) GO TO 500 ! If done already, skip this part.

      CALL gtbasfn(runname(1:lfnkey),dirname,ldirnam,hdmpfmt,           &
                   1,0,basdmpfn,lbasdmpf)

      IF (myproc == 0) PRINT*,'Output grid/base state file is ', basdmpfn(1:lbasdmpf)

      grdbas = 1      ! Dump out grd and base state arrays only

      CALL dtadump(nx,ny,nz,nzsoil,nstyps,hdmpfmt,nchout,               &
                   basdmpfn(1:lbasdmpf),grdbas,filcmprs,                &
                   u(1,1,1,isub),v(1,1,1,isub),                         &
                   w(1,1,1,isub),ptprt(1,1,1,isub),                     &
                   pprt(1,1,1,isub),qv(1,1,1,isub),                     &
                   qscalar(1,1,1,1,isub),                               &
                   tke(1,1,1,isub),kmh(1,1,1,isub),kmv(1,1,1,isub),     &
                   ubar,vbar,tem1,ptbar,pbar,rhobar,qvbar,              &
                   x,y,z,zp,zpsoil,                                     &
                   soiltyp,stypfrct,vegtyp,lai,roufns,veg,              &
                   tsoil(1,1,1,0,isub),qsoil(1,1,1,0,isub),             &
                   wetcanp(1,1,0,isub),snowdpth(1,1,isub),              &
                   raing(1,1,isub),rainc(1,1,isub),prcrate(1,1,1,isub), &
                   radfrc(1,1,1,isub),radsw(1,1,isub),rnflx(1,1,isub),  &
                   radswnet(1,1,isub),radlwin(1,1,isub),                &
                   usflx(1,1,isub),vsflx(1,1,isub),                     &
                   ptsflx(1,1,isub),qvsflx(1,1,isub),                   &
                   tem2,tem3,tem4)

      gboutcnt = 1

      500     CONTINUE

    END IF
!
!-----------------------------------------------------------------------
!
!  Then the time dependent fields:
!
!-----------------------------------------------------------------------
!
    IF( .NOT. (houtfmt == 9 .AND. vroutcnt == 1) ) THEN
!
!-----------------------------------------------------------------------
!
!  Reconstruct the file name using the specified directory name
!
!-----------------------------------------------------------------------
!
      CALL gtdmpfn(runname(1:lfnkey),dirname,                           &
                   ldirnam,curtim,hdmpfmt,1,0, hdmpfn, ldmpf)

    END IF

    IF (myproc == 0) WRITE(6,'(a,a)') 'Writing t-dependent variable history dump ', &
                                      hdmpfn(1:ldmpf)
    grdbas = 0

    CALL dtadump(nx,ny,nz,nzsoil,nstyps,hdmpfmt,nchout,                 &
                 hdmpfn(1:ldmpf),grdbas,filcmprs,                       &
                 u(1,1,1,isub),v(1,1,1,isub),                           &
                 w(1,1,1,isub),ptprt(1,1,1,isub),                       &
                 pprt(1,1,1,isub),qv(1,1,1,isub),                       &
                 qscalar(1,1,1,1,isub),                                 &
                 tke(1,1,1,isub),kmh(1,1,1,isub),kmv(1,1,1,isub),       &
                 ubar,vbar,tem1,ptbar,pbar,rhobar,qvbar,                &
                 x,y,z,zp,zpsoil,                                       &
                 soiltyp,stypfrct,vegtyp,lai,roufns,veg,                &
                 tsoil(1,1,1,0,isub),qsoil(1,1,1,0,isub),               &
                 wetcanp(1,1,0,isub),snowdpth(1,1,isub),                &
                 raing(1,1,isub),rainc(1,1,isub),prcrate(1,1,1,isub),   &
                 radfrc(1,1,1,isub),radsw(1,1,isub),rnflx(1,1,isub),    &
                 radswnet(1,1,isub),radlwin(1,1,isub),                  &
                 usflx(1,1,isub),vsflx(1,1,isub),                       &
                 ptsflx(1,1,isub),qvsflx(1,1,isub),                     &
                 tem2,tem3,tem4)

!
!-----------------------------------------------------------------------
!
!  Write out soil model variable file
!
!-----------------------------------------------------------------------
!
    IF ( sfcin == 1 ) THEN

      CALL cvttsnd( curtim, timsnd, tmstrln )

      soiloutfl = runname(1:lfnkey)//".soilvar."//timsnd(1:tmstrln)
      lfn = lfnkey + 9 + tmstrln

      IF( dirname /= ' ' ) THEN
        temchar = soiloutfl
        soiloutfl = dirname(1:ldirnam)//'/'//temchar
        lfn  = lfn + ldirnam + 1
      END IF

      !CALL fnversn(soiloutfl, lfn)

      IF (soildmp > 0) THEN
        !IF (mp_opt > 0 .AND. joindmp == 0) THEN
        !  CALL gtsplitfn(soiloutfl,1,1,loc_x,loc_y,1,1,                 &
        !                 0,0,0,lvldbg,filename,ireturn)
        !  lfn = LEN_TRIM(filename)
        !  soiloutfl(1:lfn) = TRIM(filename)
        !END IF

        !IF (myproc == 0) PRINT *, 'Writing soil data to ',soiloutfl(1:lfn)

        IF(mp_opt > 0 .AND. joindmp(FINDX_S) > 0) THEN
          CALL wrtjoinsoil(nx,ny,nzsoil,nstyps, soiloutfl(1:lfn),       &
                   dx,dy,zpsoil,                                        &
                   mapproj,trulat1,trulat2,trulon,sclfct,ctrlat,ctrlon, &
                   1,1,1,1,1,                                           &
                   tsoil(1,1,1,0,isub),qsoil(1,1,1,0,isub),             &
                   wetcanp(1,1,0,isub),snowdpth(1,1,isub),soiltyp)
        ELSE
          CALL wrtsoil(nx,ny,nzsoil,nstyps,soiloutfl(1:lfn),            &
                   dx,dy,zpsoil,                                        &
                   mapproj,trulat1,trulat2,trulon,sclfct,ctrlat,ctrlon, &
                   1,1,1,1,1,                                           &
                   tsoil(1,1,1,0,isub),qsoil(1,1,1,0,isub),             &
                   wetcanp(1,1,0,isub),snowdpth(1,1,isub),soiltyp)
        END IF

        IF (soildmp == 1 .AND. (  mp_opt == 0 .OR.                      &
                                 (myproc == 0 .AND. joindmp(FINDX_S) > 0) ) )    &
          CALL soilcntl(nx,ny,nzsoil,zpsoil,soiloutfl(1:lfn),           &
                  1,1,1,1,1,x,y)
      END IF

    END IF       ! sfcin.eq.1

!-----------------------------------------------------------------------
!
!  Write out surface property data file: sfcoutfl .
!
!-----------------------------------------------------------------------
!
    IF ( landin == 1 ) THEN

      sfcoutfl = runname(1:lfnkey)//".sfcdata"
      lfn = lfnkey + 8

      IF( dirname /= ' ' ) THEN

        temchar = sfcoutfl
        sfcoutfl = dirname(1:ldirnam)//'/'//temchar
        lfn  = lfn + ldirnam + 1

      END IF

      CALL fnversn(sfcoutfl, lfn)

      IF (sfcdmp > 0) THEN

        IF (mp_opt > 0 .AND. joindmp(FINDX_T) == 0) THEN
          CALL gtsplitfn(sfcoutfl,1,1,loc_x,loc_y,1,1,                  &
                         0,0,0,lvldbg,filename,ireturn)
          lfn = LEN_TRIM(filename)
          sfcoutfl(1:lfn) = TRIM(filename)
        END IF

        IF (myproc == 0) PRINT *, 'Write surface property data in ',sfcoutfl(1:lfn)

        IF(mp_opt > 0 .AND. joindmp(FINDX_T) > 0) THEN
          CALL wrtjoinsfcdt(nx,ny,nstyps,sfcoutfl(1:lfn), dx,dy,        &
                mapproj,trulat1,trulat2,trulon,sclfct,ctrlat,ctrlon,    &
                1,1,1,1,1,0,                                            &
                soiltyp,stypfrct,vegtyp,lai,roufns,veg,veg)
        ELSE
          CALL wrtsfcdt(nx,ny,nstyps,sfcoutfl(1:lfn), dx,dy,            &
                mapproj,trulat1,trulat2,trulon,sclfct,ctrlat,ctrlon,    &
                1,1,1,1,1,0,                                            &
                soiltyp,stypfrct,vegtyp,lai,roufns,veg,veg)
        END IF

        IF (sfcdmp == 1 .AND. ( (myproc == 0 .AND. joindmp(FINDX_T) > 0) .OR. (mp_opt == 0) ) ) &
          CALL sfccntl(nx,ny, sfcoutfl(1:lfn),1,1,1,1,1,0, x,y, tem1,tem2)

      END IF

    END IF       ! landin.eq.1
  END DO       ! l=1,nouttime


  IF (myproc == 0) WRITE(6,'(a)') ' ==== Normal successful completion of ARPSTINTRP. ===='
  CALL arpsstop(' ',0)

  100   WRITE(6,'(a)')                                                  &
          'Error reading NAMELIST file. Program ARPSTINTRP stopped.'
  CALL arpsstop(' ',9104)

  9001  CONTINUE

  WRITE(6,'(/2x,a)')'For the output grid:'
  WRITE(6,'(2x,a,f12.4)')                                               &
      'The latitude  of the output grid center, ctrlat=',ctrlat
  WRITE(6,'(2x,a,f12.4/)')                                              &
      'The longitude of the output grid center, ctrlon=',ctrlon
  WRITE(6,'(2x,a/2x,a,2f15.4,a)')                                       &
      'The SW corner (i,j)=(2,2) of the grid is located at ',           &
      '(',xgrdorg,ygrdorg,') of the input grid.'

  CALL arpsstop(' ',9001)

  9002  CONTINUE
  WRITE(6,'(1x,a,i2,/1x,a)')                                            &
      'Data read was unsuccessful. ireturn =', ireturn,                 &
      'Job stopped in ARPSTINTRP.'

  CALL arpsstop(' ',9002)
END PROGRAM arpstintrp
