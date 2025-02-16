MODULE module_arps4wrf_namelist
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! This module contains namelist variables and handle namelist IO.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  USE module_arpsgrid_constants
  USE module_wrfgrid_constants

  INTEGER                   :: max_dom
  CHARACTER(LEN=3)          :: wrf_core

  INTEGER                   :: nhisfile(MAXDOM)
  CHARACTER(LEN=MAXFILELEN) :: hisfile(MAXHIS,MAXDOM)
  CHARACTER(LEN=MAXFILELEN) :: soilfile(MAXHIS,MAXDOM)
  CHARACTER(LEN=MAXFILELEN) :: grdbasfn(MAXDOM)
  INTEGER                   :: hinfmt(MAXDOM)
  INTEGER                   :: soilopt(MAXDOM)
  LOGICAL                   :: extsoil(MAXDOM)
  INTEGER                   :: intrp_opt
  INTEGER                   :: wtrnopt(MAXDOM)

  CHARACTER(LEN=MAXFILELEN) :: geofile(MAXDOM)
  INTEGER                   :: io_form_input

  CHARACTER(LEN=MAXFILELEN) :: output_path
  INTEGER                   :: io_form_output
  INTEGER                   :: magnitude_processor
  INTEGER                   :: lvldbg

  CONTAINS

  SUBROUTINE read_namelist_params(UNTIN,UNTOUT,IAMROOT,nproc_x, nproc_y,&
                         nprocx_in,nprocy_in,max_fopen,istatus)
!-----------------------------------------------------------------------
!
! PURPOSE:
!
!   Read namelist parameters
!
!-----------------------------------------------------------------------

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: UNTIN, UNTOUT
    LOGICAL, INTENT(IN) :: IAMROOT

    INTEGER,          INTENT(OUT) :: nproc_x, nproc_y
    INTEGER,          INTENT(OUT) :: nprocx_in, nprocy_in
    INTEGER,          INTENT(OUT) :: max_fopen
    INTEGER,          INTENT(OUT) :: istatus

!-----------------------------------------------------------------------
!
! Namelist variables (local to this subroutine)
!
!-----------------------------------------------------------------------

    INTEGER                   :: hdmpinopt(MAXDOM)
    CHARACTER(LEN=MAXFILELEN) :: hdmpfheader(MAXDOM), hdmpftrailer(MAXDOM)
    CHARACTER(LEN=MAXFILELEN) :: hsoilfheader(MAXDOM)
    REAL                      :: tintv_dmpin(MAXDOM), tbgn_dmpin(MAXDOM), &
                                 tend_dmpin(MAXDOM)

!-----------------------------------------------------------------------
!
! Namelist blocks
!
!-----------------------------------------------------------------------

    NAMELIST /message_passing/ nproc_x, nproc_y, nprocx_in, nprocy_in, max_fopen
    ! these variables will be kept in both mp.inc (arps) and wrf_parallel_module (wrf)
    NAMELIST /domains/ wrf_core, max_dom
    NAMELIST /history_data/ hinfmt, hdmpinopt,                          &
                            hdmpfheader,hdmpftrailer,                   &
                            tintv_dmpin,tbgn_dmpin,tend_dmpin,          &
                            grdbasfn,nhisfile,hisfile
    NAMELIST /interp_opt/   soilopt,wtrnopt,intrp_opt, hsoilfheader
    NAMELIST /geo_data/ geofile, io_form_input
    NAMELIST /output/   output_path, io_form_output, magnitude_processor, lvldbg

!-----------------------------------------------------------------------

    LOGICAL :: fexist
    INTEGER :: n, nf
    INTEGER :: indx

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Begin of executable code ...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    istatus = 0
!-----------------------------------------------------------------------
!
! Default values
!
!-----------------------------------------------------------------------

    nproc_x   = 1
    nproc_y   = 1
    max_fopen = 8
    nprocx_in = 0
    nprocx_in = 0

    wrf_core = 'ARW'
    max_dom  = 1

    hinfmt = 3
    hdmpinopt = 1
    nhisfile  = 0
    soilopt   = 0
    wtrnopt   = 1
    intrp_opt = FOUR_PNT

    geofile(:)     = ' '
    io_form_input  = NETCDF

    output_path = './'
    io_form_output = NETCDF
    magnitude_processor = 4
    lvldbg  = 0

    grdbasfn(:)  = ' '
    hisfile(:,:) = ' '
    soilfile(:,:) = ' '

    extsoil(:) = .FALSE.
!-----------------------------------------------------------------------
!
! Read in namelist variables
!
!-----------------------------------------------------------------------

    IF (IAMROOT) THEN
      READ(UNTIN,NML=message_passing,ERR=999)
      WRITE(UNTOUT,'(1x,a)') 'Namelist block message_passing successfully read.'
      IF (nprocx_in == 1 .AND. nprocy_in == 1) max_fopen = nproc_x*nproc_y
                  ! To avoid dead loop when read and split

!--------------------------- DOMAINS -----------------------------------

      READ(UNTIN,NML=domains,ERR=999)
      WRITE(UNTOUT,'(1x,a)') 'Namelist block domains successfully read.'
      IF (wrf_core == 'NMM') max_dom = 1

!--------------------------- HISTORY_DATA ------------------------------

      READ(UNTIN,NML=history_data,ERR=999)
      WRITE(UNTOUT,'(1x,a)') 'Namelist block history_data successfully read.'

      DO n = 1,max_dom
        IF( hdmpinopt(n) == 1 ) THEN
          CALL gthinfns(hdmpfheader(n),hdmpftrailer(n),hinfmt(n),       &
                  tintv_dmpin(n), tbgn_dmpin(n), tend_dmpin(n),         &
                  grdbasfn(n),hisfile(:,n),nhisfile(n))
        END IF

        WRITE(UNTOUT,'(/2x,a,a)')'The grid base-state file name is ',   &
                    TRIM(grdbasfn(n))
!        INQUIRE(FILE=TRIM(grdbasfn(n)), EXIST=fexist)
!        IF ( .NOT. fexist ) THEN
!          WRITE(UNTOUT,'(1x,a)') 'Specified grid base-state file not found.'
!          CALL arpsstop('ERROR: ARPS grid base file not found.',1)
!        END IF

        WRITE(UNTOUT,'(2x,a,I2.2,a)') '=========== Domain: ',n,' ============'
        DO nf=1,nhisfile(n)
          WRITE(UNTOUT,'(2x,a,i5,a,a)')                                 &
          'History file No. ',nf,'   is ',trim(hisfile(nf,n))
!          INQUIRE(FILE=TRIM(hisfile(nf,n)), EXIST=fexist)
!          IF ( .NOT. fexist ) THEN
!            WRITE(UNTOUT,'(1x,a)') 'file not found.'
!            CALL arpsstop('ERROR: ARPS history file not found.',1)
!          END IF
        END DO
        WRITE(UNTOUT,'(a)') ' '
      END DO

!--------------------------- INTERP_OPT ----------------------------------

      READ(UNTIN,NML=interp_opt,ERR=999)
      WRITE(UNTOUT,'(1x,a)') 'Namelist block interp_opt successfully read.'

      DO n = 1,max_dom
        WRITE(UNTOUT,'(2x,a,I2.2,a)') '=========== Domain: ',n,' ============'

        IF (soilopt(n) > 100) THEN
          extsoil(n) = .TRUE.
          soilopt(n) = MOD(soilopt(n), 100)
          WRITE(UNTOUT,'(1x,a)') 'Soil variable will be fetched fron another set of ARPS history files as:'
          DO nf=1,nhisfile(n)
            indx = INDEX(hisfile(nf,n),'.',.TRUE.)
            !indx = LEN_TRIM(hdmpfheader(n))
            !WRITE(soilfile(nf,n),'(2a)') TRIM(hsoilfheader(n)),TRIM(hisfile(nf,n)(indx:))
            WRITE(soilfile(nf,n),'(2a)') TRIM(hsoilfheader(n)),'.hdf000000'
            WRITE(UNTOUT,'(2x,a,i5,a,a)')                               &
            'Soil file No. ',nf,'   is ',trim(soilfile(nf,n))
            !INQUIRE(FILE=TRIM(soilfile(nf,n)), EXIST=fexist)
            !IF ( .NOT. fexist ) THEN
            !  WRITE(UNTOUT,'(1x,a)') 'file not found.'
            !  CALL arpsstop('ERROR: ARPS history file not found.',1)
            !END IF
          END DO
          WRITE(UNTOUT,'(a)') ' '
        END IF
      END DO

!--------------------------- GEO_DATA ----------------------------------

      READ(UNTIN,NML=geo_data,ERR=999)
      WRITE(UNTOUT,'(1x,a)') 'Namelist block geo_data successfully read.'

      DO n = 1, max_dom
        INQUIRE(FILE=TRIM(geofile(n)), EXIST=fexist)
        IF ( .NOT. fexist ) THEN
          WRITE(UNTOUT,'(1x,3a)') 'Specified geo data file - ',         &
                                 TRIM(geofile(n)),' was not found'
          CALL arpsstop('ERROR: GEO file not found.',1)
        END IF
      END DO

!--------------------------- OUTPUT ------------------------------------

      READ(UNTIN,NML=output,ERR=999)
      WRITE(UNTOUT,'(1x,a)') 'Namelist block output successfully read.'

      CALL inquiredir(TRIM(output_path),fexist)

      IF( .NOT. fexist ) THEN
        !CALL unixcmd( 'mkdir -p '//TRIM(output_path) )
        CALL makedir(output_path,istatus)
        WRITE(UNTOUT,'(5x,a,2(/5x,a))') 'Specified output directory '   &
          //TRIM(output_path)//' not found.',                           &
          'It was created by the program.'
      END IF
   
      indx = LEN_TRIM(output_path)
      IF (output_path(indx:indx) /= '/') output_path(indx+1:indx+1) = '/' 

      WRITE(UNTOUT,'(5x,a)')                                            &
            'Output files will be in directory '//TRIM(output_path)//'.'

      GOTO 900
!-----------------------------------------------------------------------
!
! Error with reading namelist
!
!-----------------------------------------------------------------------
      999 CONTINUE
      istatus = -1

      900 CONTINUE
    END IF
    CALL mpupdatei(istatus,1)

    IF (istatus /= 0) THEN
      IF (IAMROOT) WRITE(UNTOUT,'(1x,2a)') 'Error reading NAMELIST file.', &
                         'Job stopped in read_namelist_params.'

      CALL arpsstop('Namelist input error',1)
    END IF

    CALL mpupdatei(max_dom,1)
    CALL mpupdatec(wrf_core,3)
    CALL mpupdatei(nhisfile,MAXDOM)
    CALL mpupdatec(hisfile,MAXFILELEN*MAXHIS*MAXDOM)
    CALL mpupdatec(grdbasfn,MAXFILELEN*MAXDOM)
    CALL mpupdatei(hinfmt,MAXDOM)
    CALL mpupdatei(soilopt,MAXDOM)
    CALL mpupdatel(extsoil,MAXDOM)
    CALL mpupdatec(soilfile,MAXFILELEN*MAXHIS*MAXDOM)

    CALL mpupdatei(wtrnopt,MAXDOM)
    CALL mpupdatei(intrp_opt,1)
    CALL mpupdatec(geofile,MAXDOM)
    CALL mpupdatei(io_form_input,1)
    CALL mpupdatec(output_path,MAXFILELEN)
    CALL mpupdatei(io_form_output,1)
    CALL mpupdatei(magnitude_processor,1)
!    CALL mpupdatei(lvldbg,1)  ! only Root output debug message

    CALL mpupdatei(nproc_x,1)
    CALL mpupdatei(nproc_y,1)
    CALL mpupdatei(nprocx_in,1)
    CALL mpupdatei(nprocy_in,1)
    CALL mpupdatei(max_fopen,1)

    RETURN
  END SUBROUTINE read_namelist_params

END MODULE module_arps4wrf_namelist

