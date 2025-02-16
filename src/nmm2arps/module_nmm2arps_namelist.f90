MODULE module_nmm2arps_namelist
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! This module contains namelist variables and handle namelist IO.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  USE module_wrfgrid_constants
  USE module_arpsgrid_constants

  INTEGER                   :: nwrfdfil
  CHARACTER(LEN=MAXFILELEN) :: wrffiles(MAXHIS)
  INTEGER                   :: frames_per_outfile(MAXHIS)
  INTEGER                   :: io_form, grid_id,magnitude_processor
  CHARACTER(LEN=19)         :: init_time_str
  CHARACTER(LEN=19)         :: datestrs(MAXHIS)
  INTEGER                   :: initsec, abstimei, abstimee
  CHARACTER(LEN=MAXFILELEN) :: dir_extd, dir_extp

  INTEGER                   :: nx, ny, nz
  REAL                      :: dx, dy, dz
  INTEGER                   :: mapproj
  REAL                      :: trulat1, trulat2, trulon, sclfct
  REAL                      :: ctrlat, ctrlon

  INTEGER                   :: strhopt
  REAL                      :: dzmin, zrefsfc, dlayer1, dlayer2, strhtune
  REAL                      :: zflat
  INTEGER                   :: nzsoil, nstyp, soilmodel_option
  REAL                      :: dzsoil


  INTEGER                   :: intropt,intrp_method
  INTEGER                   :: nsmooth
  INTEGER                   :: ext_lbc, ext_vbc, bbc, tbc, fftopt
  INTEGER                   :: wrfexttrnopt
  CHARACTER(LEN=MAXFILELEN) :: terndta
  INTEGER                   :: ternfmt, extntmrg

  INTEGER                   :: csopt, hydradj, wndadj, obropt
  REAL                      :: csfactr, csound, obrzero

  CHARACTER(LEN=80)         :: runname
  INTEGER                   :: nocmnt
  CHARACTER(LEN=80)         :: cmnt(50)

  CHARACTER(LEN=MAXFILELEN) :: dirname
  INTEGER                   :: dmp_out_joined
  INTEGER                   :: exbcdmp, qcexout, qrexout, qiexout, qsexout, qhexout
  INTEGER                   :: soildmp, terndmp
  INTEGER                   :: hdmpfmt, hdfcompr
  INTEGER                   :: basout, grdout, varout, mstout
  INTEGER                   :: iceout, tkeout, trbout, rainout
  INTEGER                   :: sfcout, landout,prcout,radout, flxout
  INTEGER                   :: snowout, totout
  INTEGER                   :: filcmprs, readyfl

  CHARACTER(LEN=80)         :: outheader, gemoutheader
  INTEGER                   :: icape, iaccu, iascii, i2dfmt, igempak
  INTEGER                   :: ilite,iltgci,icrtm,isatid,chbgn,chend,icitm,user_emis
  INTEGER                   :: ibeg_offset,iend_offset,jbeg_offset,jend_offset

  INTEGER                   :: lvldbg

  LOGICAL :: multifile
  INTEGER :: ncompressx, ncompressy ! compression in x and y direction:
                                    ! ncompressx=nprocx_in/nproc_x
                                    ! ncompressy=nprocy_in/nproc_y

  CONTAINS

  SUBROUTINE read_namelist_params(UNTIN,UNTOUT,IAMROOT,mp_opt,          &
                      nproc_x, nproc_y,nprocx_in,nprocy_in,max_fopen,   &
                      istatus)
!-----------------------------------------------------------------------
!
! PURPOSE:
!
!   Read namelist parameters
!
!-----------------------------------------------------------------------

    IMPLICIT NONE
    INTEGER, INTENT(IN)  :: UNTIN, UNTOUT
    LOGICAL, INTENT(IN)  :: IAMROOT
    INTEGER, INTENT(IN)  :: mp_opt

    INTEGER, INTENT(OUT) :: nproc_x, nproc_y
    INTEGER, INTENT(OUT) :: nprocx_in, nprocy_in
    INTEGER, INTENT(OUT) :: max_fopen
    INTEGER, INTENT(OUT) :: istatus

!-----------------------------------------------------------------------
!
! Namelist variables (local to this subroutine)
!
!-----------------------------------------------------------------------

    CHARACTER(LEN=19) :: start_time_str, end_time_str
    CHARACTER(LEN=11) :: history_interval

!-----------------------------------------------------------------------
!
! Namelist blocks
!
!-----------------------------------------------------------------------

    NAMELIST /message_passing/ nproc_x, nproc_y, nprocx_in, nprocy_in,  &
                               max_fopen
    ! these variables will be kept in both mp.inc (arps) and wrf_parallel_module (wrf)
    NAMELIST /wrfdfile/ dir_extd,dir_extp,init_time_str,io_form,grid_id,&
                        start_time_str, history_interval,end_time_str,  &
                        magnitude_processor

    NAMELIST /arpsgrid/ nx,ny,nz,nzsoil,nstyp,dx,dy,                    &
                        strhopt,dzmin,zrefsfc,dlayer1,dlayer2,strhtune, &
                        zflat,dz,dzsoil,soilmodel_option,               &
                        mapproj,trulat1,trulat2,                        &
                        trulon,sclfct,ctrlat,ctrlon

    NAMELIST /intrp_opts/ intrp_method, intropt, nsmooth, ext_lbc, ext_vbc,   &
                          bbc, tbc, fftopt,                             &
                          wrfexttrnopt, terndta, ternfmt, extntmrg

    NAMELIST /adjust/   csopt,csfactr,csound,hydradj,wndadj,obropt,obrzero

    NAMELIST /comment_lines/ runname,nocmnt,cmnt,lvldbg

    NAMELIST /output/     dirname,dmp_out_joined,exbcdmp,soildmp,terndmp, &
                          hdmpfmt,hdfcompr,                             &
                          qcexout,qrexout,qiexout,qsexout,qhexout,      &
                          basout,grdout,varout,mstout,iceout,tkeout,    &
                          trbout,rainout,sfcout,landout,prcout,         &
                          radout,flxout,snowout,totout,filcmprs,readyfl

    NAMELIST /output_2d/  outheader,gemoutheader,icape,iaccu,iascii,    &
                          i2dfmt,igempak,ilite,iltgci,icrtm,isatid,chbgn,chend,&
                          icitm,user_emis,                              &
                          ibeg_offset,iend_offset,jbeg_offset,jend_offset

!-----------------------------------------------------------------------

    INTEGER          :: strlen
    LOGICAL          :: rewindyr
    CHARACTER(LEN=1) :: ach
    INTEGER          :: myr,year,month,day,hour,minute,second
    INTEGER          :: abstimes

    LOGICAL :: fexist
    INTEGER :: n, nf, ifile

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Begin of executable code ...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    istatus = 0
!
!-----------------------------------------------------------------------
!
!  Read in message passing options.
!
!-----------------------------------------------------------------------
!
    IF (IAMROOT) THEN
      READ (UNTIN,message_passing)
      WRITE(UNTOUT,'(a)')'  Namelist block message_passing sucessfully read.'
      IF (mp_opt == 0) THEN
      	nproc_x = 1
      	nproc_y = 1
      END IF
    END IF
    CALL mpupdatei(nproc_x,1)
    CALL mpupdatei(nproc_y,1)
    CALL mpupdatei(max_fopen,1)
    CALL mpupdatei(nprocx_in,1)
    CALL mpupdatei(nprocy_in,1)
!
!-----------------------------------------------------------------------
!
!  Read in namelist &wrfdfile
!
!-----------------------------------------------------------------------
!
    dir_extd = './'
    dir_extp = ''

    init_time_str         = '0000-00-00_00:00:00'
    start_time_str        = '0000-00-00_00:00:00'
    history_interval      = '00_00:00:00'
    end_time_str          = '0000-00-00_00:00:00'

    io_form               = 7
    grid_id               = 1
    magnitude_processor   = 4

    multifile             = .FALSE.          ! not in namelist
    IF (IAMROOT) THEN
      READ(UNTIN,wrfdfile)
      WRITE(UNTOUT,'(a)') '  Namelist wrfdfile read in successfully.'

      strlen = LEN_TRIM(dir_extd)
      IF(strlen > 0) THEN
        IF(dir_extd(strlen:strlen) /= '/' .AND. dir_extd(strlen:strlen) /= '>' ) THEN
          dir_extd(strlen+1:strlen+1) = '/'
          strlen = strlen + 1
        END IF
      ELSE
        dir_extd = './'
      END IF

      IF (LEN_TRIM(dir_extp) == 0 .OR. dir_extp(1:1) == ' ') THEN
        dir_extp = dir_extd
      ELSE IF ( dir_extp(1:5) == '<dir>') THEN
      ELSE
        strlen = LEN_TRIM(dir_extp)
        IF(strlen > 0) THEN
          IF(dir_extp(strlen:strlen) /= '/') THEN
            dir_extp(strlen+1:strlen+1) = '/'
          END IF
        ELSE
          dir_extp = './'
        END IF
      END IF

      IF (io_form > 100) THEN
        io_form = MOD(io_form,100)
        multifile = .TRUE.
      END IF

      rewindyr = .FALSE.

      READ(start_time_str,'(I4.4,5(a,I2.2))')      &
                      year,ach,month,ach,day,ach,hour,ach,minute,ach,second
      IF (year < 1960) THEN   ! maybe ideal case
        myr  =  year
        year =  1960
        rewindyr = .TRUE.
      END IF
      CALL ctim2abss(year,month,day,hour,minute,second,abstimes)

      READ(end_time_str,'(I4.4,a,I2.2,a,I2.2,a,I2.2,a,I2.2,a,I2.2)')    &
                   year,ach,month,ach,day,ach,hour,ach,minute,ach,second
      IF (rewindyr)  year = 1960
      CALL ctim2abss(year,month,day,hour,minute,second,abstimee)

      READ(history_interval,'(I2.2,a,I2.2,a,I2.2,a,I2.2)')               &
                                       day,ach,hour,ach,minute,ach,second
      abstimei = day*24*3600+hour*3600+minute*60+second

      IF (multifile) THEN
        IF (MOD(nprocx_in,nproc_x) /= 0 .OR. MOD(nprocy_in,nproc_y) /= 0) THEN
          WRITE(UNTOUT,*) 'nprocx_in (nprocy_in) must be dividable by nproc_x (nproc_y).'
          istatus = -1
        END IF

        ncompressx = nprocx_in/nproc_x
        ncompressy = nprocy_in/nproc_y

      ELSE     ! non-mpi or mpi with one file

        ncompressx = 1
        ncompressy = 1

      END IF

      frames_per_outfile(:) = 1
      nwrfdfil = 0

      CALL check_wrf_files(multifile,MAXHIS,grid_id,io_form,            &
        nprocx_in,  magnitude_processor, ncompressx,ncompressy,         &
        abstimes,abstimei,abstimee,rewindyr,myr,                        &
        dir_extd,wrffiles,nwrfdfil,frames_per_outfile,nproc_x,istatus)

!-----------------------------------------------------------------------
!
!   Set datestrs
!
!-----------------------------------------------------------------------

      ifile = abstimes
      DO WHILE (ifile <= abstimee)

        CALL abss2ctim(ifile,year,month,day,hour,minute,second)
        IF (rewindyr) year = myr

        n = (ifile - abstimes)/abstimei + 1
        WRITE(datestrs(n),'(I4.4,5(a,I2.2))')                       &
           year,'-',month,'-',day,'_',hour,':',minute,':',second

        ifile = ifile + abstimei
      END DO

      READ(init_time_str,'(I4.4,5(a,I2.2))')      &
                  year,ach,month,ach,day,ach,hour,ach,minute,ach,second
      IF (year < 1960) THEN   ! maybe ideal case
        year =  1960
      END IF
      CALL ctim2abss(year,month,day,hour,minute,second,initsec)

    END IF   ! IAMROOT
    CALL mpupdatei(istatus,1)
    IF (istatus /= 0) CALL arpsstop('',1)

    CALL mpupdatei(ncompressx,1)
    CALL mpupdatei(ncompressy,1)

    CALL mpupdatei(io_form,1)
    CALL mpupdatel(multifile,1)
    CALL mpupdatec(init_time_str,19)
    CALL mpupdatei(grid_id, 1)

    CALL mpupdatei(nwrfdfil,1)
    CALL mpupdatei(frames_per_outfile,MAXHIS)
    CALL mpupdatei(magnitude_processor,1)
    CALL mpupdatec(wrffiles,MAXFILELEN*MAXHIS)
    CALL mpupdatec(datestrs,19*MAXHIS)
    !CALL mpupdatec(dir_extd,256)
    CALL mpupdatec(dir_extp,256)

    CALL mpupdatei(initsec,  1)
    CALL mpupdatei(abstimei, 1)
    CALL mpupdatei(abstimee, 1)
!
!-----------------------------------------------------------------------
!
!  Read in namelist &arpsgrid
!
!-----------------------------------------------------------------------
!
    nx = 67
    ny = 67
    nz = 35
    nzsoil = 2
    nstyp  = 4
    dx     = 3200
    dy     = 3200
    dz     = 500
    dzsoil = 1

    strhopt  = 0
    dzmin    = 100.000
    zrefsfc  =   0.0
    dlayer1  =   0.0
    dlayer2  =   1.0e5
    strhtune =   1.0
    zflat    =   1.0e5

    soilmodel_option = 1
    mapproj  = 2
    trulat1  = 30.0
    trulat2  = 60.0
    trulon   = -100.0
    sclfct   = 1.0
    ctrlat   = 35.0
    ctrlon   = -98.0

    IF (IAMROOT) THEN
      READ(UNTIN,arpsgrid)
      WRITE(UNTOUT,'(a)') '  Namelist arpsgrid read in successfully.'

      ! Note that for MP version nx & ny here are global values.  They will
      ! be reassigned to their per-processor value below.

      IF (nx /= nproc_x*((nx-3)/nproc_x)+3) THEN
        WRITE (6,'(1x,2(a,I4),a)') 'WARNING: nx = ',nx,' does not fit within nproc_x = ',nproc_x,' processors.'
        WRITE (6,'(1x,a,/,1x,a)' ) 'nmm2arps_mpi cannot handle this case at present.', &
                    'Use no-mpi version. Program stopping ...'
        CALL arpsstop('nx does not fit.',1)
      END IF
      IF (ny /= nproc_y*int((ny-3)/nproc_y)+3) THEN
        WRITE (6,'(1x,2(a,I4),a)') 'WARNING: ny = ',ny,' does not fit within nproc_y = ',nproc_y,' processors.'
        WRITE (6,'(1x,a,/,1x,a)' ) 'nmm2arps_mpi cannot handle this case at present.', &
                    'Use no-mpi version Program stopping ...'
        CALL arpsstop('ny does not fit.',1)
      END IF

      nx = (nx - 3)/nproc_x + 3             ! fake zone for ARPS is 3
      ny = (ny - 3)/nproc_y + 3

    END IF
    CALL mpupdatei(nx,1)
    CALL mpupdatei(ny,1)
    CALL mpupdater(dx,1)
    CALL mpupdater(dy,1)
    CALL mpupdatei(mapproj,1)
    CALL mpupdater(trulat1,1)
    CALL mpupdater(trulat2,1)
    CALL mpupdater(trulon,1)
    CALL mpupdater(sclfct,1)
    CALL mpupdater(ctrlat,1)
    CALL mpupdater(ctrlon,1)
    CALL mpupdatei(strhopt,1)
    CALL mpupdater(dzmin,1)
    CALL mpupdater(dlayer1,1)
    CALL mpupdater(dlayer2,1)
    CALL mpupdater(strhtune,1)
    CALL mpupdater(zflat,1)
    CALL mpupdatei(nz,1)
    CALL mpupdater(dz,1)
    CALL mpupdatei(nzsoil,1)
    CALL mpupdater(dzsoil,1)
    CALL mpupdatei(nstyp,1)
    CALL mpupdatei(soilmodel_option,1)
!
!-----------------------------------------------------------------------
!
!  Read in namelist &intrp_opts
!
!-----------------------------------------------------------------------
!
    intrp_method   = 2
    intropt  = 1
    nsmooth  = 1
    ext_lbc  = 1
    ext_vbc  = 1
    bbc      = 1
    tbc      = 1
    fftopt   = 2
    wrfexttrnopt= 2
    extntmrg = 1
    terndta  = 'arpstern.data'
    ternfmt  = 1
    IF (IAMROOT) THEN
      READ (UNTIN,intrp_opts)
      WRITE(UNTOUT,'(a)') '  Namelist intrp_opts read in successfully.'

      IF(wrfexttrnopt == 2 .OR. wrfexttrnopt == 3) THEN
        INQUIRE(FILE=TRIM(terndta), EXIST = fexist)
        IF(.NOT. fexist) THEN
          WRITE(UNTOUT,*) 'The tern file ', TRIM(terndta),                       &
                     ' you specified does not exist.'
          CALL arpsstop( 'tern file not found',1)
        END IF
      END IF
    END IF
    CALL mpupdatei(intrp_method,1)
    CALL mpupdatei(intropt,1)
    CALL mpupdatei(nsmooth,1)
    CALL mpupdatei(ext_lbc,1)
    CALL mpupdatei(ext_vbc,1)
    CALL mpupdatei(bbc,1)
    CALL mpupdatei(tbc,1)
    CALL mpupdatei(fftopt,1)
    CALL mpupdatei(wrfexttrnopt,1)
    CALL mpupdatec(terndta,LEN(terndta))
    CALL mpupdatei(ternfmt,1)
!
!-----------------------------------------------------------------------
!
!  Read in namelist &adjust
!
!-----------------------------------------------------------------------
!
    csopt    = 1
    csfactr  = 0.5
    csound   = 150.0
    hydradj  = 0
    wndadj   = 2
    obropt   = 12
    obrzero  = 12000.

    IF (IAMROOT) THEN
      READ(UNTIN,adjust)
      WRITE(UNTOUT,'(a)') '  Namelist adjust read in successfully.'
    END IF
    CALL mpupdatei(csopt,1)
    CALL mpupdater(csfactr,1)
    CALL mpupdater(csound,1)
    CALL mpupdatei(hydradj,1)
    CALL mpupdatei(wndadj,1)
    CALL mpupdatei(obropt,1)
    CALL mpupdatei(obrzero,1)

!
!-----------------------------------------------------------------------
!
!  Read in namelist &comment_lines
!
!-----------------------------------------------------------------------
!
    runname = 'WRF2ARPS, Version 2.2'
    nocmnt  = 2
    cmnt(1) = ' '
    cmnt(2) = ' '

    IF (IAMROOT) THEN
      READ(UNTIN,comment_lines)
      WRITE(UNTOUT,'(a)') '  Namelist comment_lines read in successfully.'
    END IF
    CALL mpupdatec(runname,LEN(runname))
    CALL mpupdatei(nocmnt,1)
    CALL mpupdatec(cmnt(1),LEN(cmnt(1)))
    CALL mpupdatec(cmnt(2),LEN(cmnt(2)))
  !  CALL mpupdatei(lvldbg,1)
!
!-----------------------------------------------------------------------
!
!  Read in namelist &output
!
!-----------------------------------------------------------------------
!
    dirname  = './'
    exbcdmp  = 1
    soildmp  = 0
    terndmp  = 0
    hdmpfmt  = 1
    hdfcompr = 0

    qcexout = 0
    qrexout = 0
    qiexout = 0
    qsexout = 0
    qhexout = 0
    grdout   = 0

    basout   = 0
    varout   = 1
    mstout   = 1
    rainout  = 0
    prcout   = 0
    iceout   = 0
    totout   = 1
    tkeout   = 1
    trbout   = 0
    sfcout   = 0
    snowout  = 0
    landout  = 0
    radout   = 0
    flxout   = 0
    IF (IAMROOT) THEN
      READ(UNTIN,output)
      WRITE(UNTOUT,'(a)') '  Namelist output read in successfully.'
    END IF
    CALL mpupdatec(dirname,LEN(dirname))
    CALL mpupdatei(dmp_out_joined,1)
    CALL mpupdatei(exbcdmp,1)
    CALL mpupdatei(qcexout,1)
    CALL mpupdatei(qrexout,1)
    CALL mpupdatei(qiexout,1)
    CALL mpupdatei(qsexout,1)
    CALL mpupdatei(qhexout,1)
    CALL mpupdatei(soildmp,1)
    CALL mpupdatei(terndmp,1)
    CALL mpupdatei(hdmpfmt,1)
    CALL mpupdatei(hdfcompr,1)
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
    CALL mpupdatei(prcout,1)
    CALL mpupdatei(radout,1)
    CALL mpupdatei(flxout,1)
    CALL mpupdatei(filcmprs,1)
    CALL mpupdatei(readyfl,1)

!-----------------------------------------------------------------------
!
! output_2d
!
!-----------------------------------------------------------------------

    outheader    = 'output_2d'
    gemoutheader = 'output_2d_gem'
    icape = 1
    iaccu = 1
    icrtm = 0
    icitm = 0
    user_emis = 0
    isatid = 2
    chbgn = 3
    chend = 3
    iascii  = 0
    i2dfmt = 0
    igempak = 0
    ilite = 0
    iltgci = 0
    ibeg_offset = 0
    iend_offset = 0
    jbeg_offset = 0
    jend_offset = 0

    IF(IAMROOT) THEN
      READ(UNTIN,output_2d,END=200)
      WRITE(UNTOUT,'(a)') '  Namelist output_2d read in successfully.'
  !    WRITE(UNTOUT,output_2d)
    GO TO 20
    200   WRITE(UNTOUT,'(a)')                                           &
           'Error reading NAMELIST block output_2d.  ',                 &
           'Default values used.'
    20 CONTINUE
    END IF
    CALL mpupdatec(outheader,80)
    CALL mpupdatec(gemoutheader,80)
    CALL mpupdatei(icape,1)
    CALL mpupdatei(iaccu,1)
    CALL mpupdatei(iascii,1)
    CALL mpupdatei(i2dfmt,1)
    CALL mpupdatei(igempak,1)
    CALL mpupdatei(ilite,1)
    CALL mpupdatei(iltgci,1)
    CALL mpupdatei(icrtm,1)
    CALL mpupdatei(icitm,1)
    CALL mpupdatei(user_emis,1)
    CALL mpupdatei(isatid,1)
    CALL mpupdatei(chbgn,1)
    CALL mpupdatei(chend,1)
    CALL mpupdatei(ibeg_offset,1)
    CALL mpupdatei(iend_offset,1)
    CALL mpupdatei(jbeg_offset,1)
    CALL mpupdatei(jend_offset,1)

!-----------------------------------------------------------------------

    RETURN
  END SUBROUTINE read_namelist_params

END MODULE module_nmm2arps_namelist

!
!##################################################################
!##################################################################
!######                                                      ######
!######              SUBROUTINE check_wrf_files              ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE check_wrf_files(multifile,MAXWRFFIL,grid_id,io_form,         &
           nprocx_in, numdigits,ncompressx,ncompressy,                  &
           abstimes,abstimei,abstimee,rewindyr,myr,                     &
           dir_extd,extdname,nextdfil,frames_per_outfile,nproc_xa,istatus)

!-----------------------------------------------------------------------
!
! PURPOSE: Check the existence of WRF files to be read and return the
!          valid file number, file names for each processor, number of
!          frames in each WRF file.
!
!-----------------------------------------------------------------------
!
!  AUTHOR:
!  Yunheng Wang (03/26/2007)
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

  LOGICAL, INTENT(IN)    :: multifile
  INTEGER, INTENT(IN)    :: MAXWRFFIL
  INTEGER, INTENT(IN)    :: grid_id
  INTEGER, INTENT(IN)    :: io_form
  INTEGER, INTENT(IN)    :: nprocx_in
  INTEGER, INTENT(IN)    :: ncompressx, ncompressy, numdigits
  INTEGER, INTENT(IN)    :: abstimes, abstimei, abstimee
  LOGICAL, INTENT(IN)    :: rewindyr
  INTEGER, INTENT(IN)    :: myr
  CHARACTER(LEN=256), INTENT(IN)  :: dir_extd
  CHARACTER(LEN=256), INTENT(OUT) :: extdname(MAXWRFFIL)
  INTEGER,            INTENT(OUT) :: nextdfil
  INTEGER,            INTENT(OUT) :: frames_per_outfile(MAXWRFFIL)
  INTEGER,            INTENT(OUT) :: nproc_xa
  INTEGER,            INTENT(OUT) :: istatus

!-----------------------------------------------------------------------
!
! Misc. local variables
!
!-----------------------------------------------------------------------

  INTEGER :: iloc, jloc
  INTEGER :: iloc_x, jloc_y, loc_proc

  INTEGER :: year, month, day, hour, minute, second
  INTEGER :: ifile

  CHARACTER(LEN=256) :: tmpstr

  LOGICAL :: fexist
  LOGICAL :: First_in_process

  INTEGER :: itmp
  INTEGER :: loc_xa, loc_ya

  CHARACTER(LEN=20) :: fmtstr
  CHARACTER(LEN=40) :: tmpwrftime

!-----------------------------------------------------------------------
!
! Include files
!
!-----------------------------------------------------------------------

  INCLUDE 'mp.inc'

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code ....
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  nextdfil    = 0
  extdname(:) = ' '
  istatus     = 0

  IF (myproc == 0) WRITE(6,'(1x,a,/,1x,a)')                             &
              '============================','WRF files to be read are:'

  loc_xa = MOD(myproc, nproc_xa) + 1
  loc_ya = myproc / nproc_xa + 1

  IF (multifile) THEN

    iloc_x = (loc_xa-1)*ncompressx    ! column of processors
    jloc_y = (loc_ya-1)*ncompressy    ! rows of processors

    WRITE(fmtstr,'(a,2(I1,a))') '(2a,I',numdigits,'.',numdigits,')'

    ifile = abstimes
    nextdfil_loop1: DO WHILE(ifile <= abstimee)

      CALL abss2ctim(ifile,year,month,day,hour,minute,second)
      IF (rewindyr) year = myr

      nextdfil = nextdfil + 1
      WRITE(tmpwrftime,'(a,I2.2,a,I4.4,5(a,I2.2))')                     &
             'wrfout_d',grid_id,'_',                                    &
             year,'-',month,'-',day,'_',hour,':',minute,':',second

      IF (dir_extd == '<dir>') THEN
        WRITE(extdname(nextdfil),'(3a)') TRIM(tmpwrftime),'/',TRIM(tmpwrftime)
      ELSE
        WRITE(extdname(nextdfil),'(2a)') TRIM(dir_extd),TRIM(tmpwrftime)
      END IF

      First_in_process = .TRUE.
      DO jloc = 1,ncompressy
        DO iloc = 1, ncompressx

          loc_proc = (jloc_y+jloc-1)*nprocx_in + iloc_x+(iloc-1)
          WRITE(tmpstr,FMT=fmtstr) TRIM(extdname(nextdfil)),'_',loc_proc

          INQUIRE(FILE=TRIM(tmpstr), EXIST = fexist )
          IF(.NOT. fexist) THEN
            WRITE(6,'(1x,a)') 'WARNING: The WRF file ',TRIM(tmpstr),    &
                             ' does not exist.'
            nextdfil = nextdfil - 1
            ifile = ifile + abstimei   ! try next time level
            CYCLE nextdfil_loop1
          ELSE
            CALL get_wrf_frames_per_outfile(TRIM(tmpstr),io_form,itmp,istatus)
            IF ( First_in_process) THEN
              frames_per_outfile(nextdfil) = itmp
              First_in_process = .FALSE.
              IF(myproc == 0) WRITE(6,'(5x,a,I2.2,3a,I4,a)')            &
                 'WRF file ',nextdfil,': ', TRIM(tmpstr),               &
                 ', frames = ',frames_per_outfile(nextdfil),'.'
            ELSE
              ! Ensure all files read by this processor have same frames
              IF (itmp /= frames_per_outfile(nextdfil) ) THEN
                WRITE(6,'(1x,a,I4,2(a,I2),a,/,10x,a,I4,a,I2,a,I4,a,/)') &
                  'WARNING: unlimited dimension in processor - ',myproc,&
                  ' for file (',iloc,', ',jloc,                         &
                  ') is different with others files.',                  &
                  'We got iframes = ',itmp,                             &
                  ', but other files have frames_per_outfile(',         &
                  nextdfil,') = ',frames_per_outfile(nextdfil),'.'
                istatus = -1
                RETURN
              END IF
            END IF
          END IF
        END DO
      END DO

      !
      ! Ensure all processors have the same frames
      !
!      itmp = frames_per_outfile(nextdfil)
!      CALL mpmaxi(itmp)
!      IF ( itmp /= frames_per_outfile(nextdfil) ) THEN
!        WRITE(6,'(1x,a,/,8x,3(a,I4),/,8x,a,I4,a,/)')                    &
!          'ERROR: Each processor has got different numbers for frames in WRF files.', &
!          'Processor: ',myproc,' found ',frames_per_outfile(nextdfil),  &
!          ' frames in WRF file - ',nextdfil,                            &
!          'While other processor may found ',itmp,' frames in that file.'
!        istatus = -2
!        RETURN
!      END IF

      ifile = ifile + frames_per_outfile(nextdfil)*abstimei
    END DO nextdfil_loop1

  ELSE     ! non-mpi or mpi with one file

    ifile = abstimes
    nextdfil_loop2: DO WHILE (ifile <= abstimee)

      CALL abss2ctim(ifile,year,month,day,hour,minute,second)
      IF (rewindyr) year = myr

      nextdfil = nextdfil + 1
      WRITE(tmpwrftime,'(a,I2.2,a,I4.4,5(a,I2.2))')                     &
             'wrfout_d',grid_id,'_',                                    &
             year,'-',month,'-',day,'_',hour,':',minute,':',second

      IF (dir_extd == '<dir>') THEN
        WRITE(extdname(nextdfil),'(3a)') TRIM(tmpwrftime),'/',TRIM(tmpwrftime)
      ELSE
        WRITE(extdname(nextdfil),'(2a)') TRIM(dir_extd),TRIM(tmpwrftime)
      END IF

      INQUIRE(FILE=TRIM(extdname(nextdfil)), EXIST = fexist )
      IF(.NOT. fexist) THEN
        WRITE(6,'(1x,a)') 'WARNING: The WRF file ',                   &
                       TRIM(extdname(nextdfil)),' does not exist.'
        ifile = ifile + abstimei
        nextdfil = nextdfil - 1
        CYCLE nextdfil_loop2
      ELSE
        CALL get_wrf_frames_per_outfile(TRIM(extdname(nextdfil)),     &
                                        io_form,itmp,istatus)
        frames_per_outfile(nextdfil) = itmp
        WRITE(6,'(5x,a,I2.2,3a,I4,a)')                                &
           'WRF file ',nextdfil,': ', TRIM(extdname(nextdfil)),       &
           ', frames = ',frames_per_outfile(nextdfil),'.'
      END IF

      ifile = ifile + frames_per_outfile(nextdfil)*abstimei
    END DO nextdfil_loop2

  END IF     ! multifile

!-----------------------------------------------------------------------
!
! Validate nextdfil before return
!
!-----------------------------------------------------------------------

  IF(nextdfil < 1) THEN
    WRITE(6,'(a)') 'No input WRF file was valid. Please check the input file.'
    istatus = -3
    RETURN
  END IF

!  itmp = nextdfil
!  CALL mpmaxi(nextdfil)
!  IF (itmp /= nextdfil) THEN
!    WRITE(6,'(a)') 'Each processor may process different number of WRF files.'
!    WRITE(6,'(2(a,I4),a,/,a,I4,a)') 'Processor: ',myproc,' found ',   &
!                    itmp,' WRF files,',                               &
!                    'While other processor may found ',nextdfil,' files.'
!    istatus = -4
!    RETURN
!  END IF

  RETURN
END SUBROUTINE check_wrf_files
!
!##################################################################
!##################################################################
!######                                                      ######
!######       SUBROUTINE get_wrf_frames_per_outfile          ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE get_wrf_frames_per_outfile(filename,io_form,numframes,istatus)
!-----------------------------------------------------------------------
!
!  PURPOSE:
!    Get the size of unlimited dimension in the WRF data file
!
!-----------------------------------------------------------------------
!
!  AUTHOR:
!  Yunheng Wang (03/26/2007)
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

  CHARACTER(LEN=*), INTENT(IN)  :: filename
  INTEGER,          INTENT(IN)  :: io_form
  INTEGER,          INTENT(OUT) :: numframes
  INTEGER,          INTENT(OUT) :: istatus

!------------------------------------------------------------------
!
!  Misc. local variables
!
!------------------------------------------------------------------

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  istatus = 0
!  IF (io_form == NETCDF) THEN
    numframes = 1
!    CALL get_ncd_frames_per_outfile(TRIM(filename),numframes,istatus)
!  ELSE
!    istatus   = 1
!    numframes = 1
!    WRITE(6,'(1x,a,/,10x,a,/)')       &
!      'WARNING: Multiple frames in one file is only support on netCDF file at present', &
!      'frames_per_outfile is reset to 1. Pleace check for correctness.'
!  END IF

  RETURN
END SUBROUTINE get_wrf_frames_per_outfile
