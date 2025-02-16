PROGRAM arps4wrf
!
!##################################################################
!##################################################################
!######                                                      ######
!######                  PROGRAM ARPS4WRF                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
!-----------------------------------------------------------------------
!
! PURPOSE:
!
! Prepare ARPS data for WRF (either ARW core or NMM core)
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yunheng Wang
!  10/23/2008
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
  USE module_wrfgrid_constants
  USE module_arpsgrid_constants
  USE module_arps4wrf_namelist
  USE wrf_parallel_module

  IMPLICIT NONE

!-----------------------------------------------------------------------
!
! Misc. Local variables
!
!-----------------------------------------------------------------------

  LOGICAL :: IAMROOT

  INTEGER :: UNTIN  = STDIN
  INTEGER :: UNTOUT = STDOUT
  INTEGER :: UNTERR = STDERR

  INTEGER :: istatus
  INTEGER :: nargs

  INTEGER            :: nprocx_in, nprocy_in
  CHARACTER(LEN=256) :: NML_FILE

  INTEGER :: n

!-----------------------------------------------------------------------
!
! Including files
!
!-----------------------------------------------------------------------

  INCLUDE 'mp.inc'

  INTEGER :: iargc

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Begin of executable code
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  istatus = 0

  CALL mpinit_proc(0)

  IAMROOT = .FALSE.
  IF (myproc == ROOT) IAMROOT = .TRUE.

!-----------------------------------------------------------------------
!
! Process command line for namelist input if present,
! otherwise read STDIN.
!
!-----------------------------------------------------------------------

  IF( IAMROOT ) WRITE(6,'(9(/5x,a),/)')                                 &
      '###############################################################',&
      '###############################################################',&
      '####                                                       ####',&
      '####              Welcome to ARPS4WRF                      ####',&
      '####                                                       ####',&
      '####           Prepare ARPS data for WRF runs.             ####',&
      '####                                                       ####',&
      '###############################################################',&
      '###############################################################'

  IF ( IAMROOT ) THEN
    nargs = COMMAND_ARGUMENT_COUNT()
    IF (nargs == 1) THEN
      !CALL getarg(nargs, NML_FILE)
      CALL GET_COMMAND_ARGUMENT(1, NML_FILE)
      CALL getunit( UNTIN )
      WRITE(UNTOUT,'(1x,3a)') 'Reading namelist parameters from file "',&
                              TRIM(NML_FILE),'"'
      OPEN(UNIT=UNTIN,FILE=TRIM(NML_FILE),STATUS='OLD',ACTION='READ',   &
           FORM='FORMATTED',IOSTAT=istatus)
    ELSE
      WRITE(UNTOUT,'(1x,a)') 'Reading namelist parameters from STDIN.'
    END IF
  END IF

  CALL mpupdatei(istatus,1)
  IF (istatus /= 0) THEN
    IF (IAMROOT) WRITE(UNTOUT,'(1x,3a,I3)') 'IO ERROR with file ',      &
                 TRIM(NML_FILE),', with error code = ',istatus
    CALL arpsstop('IO ERROR with namelist file.',1)
  END IF

!-----------------------------------------------------------------------
!
! Read namelist parameters
!
!-----------------------------------------------------------------------

  CALL read_namelist_params(UNTIN,UNTOUT,IAMROOT,nproc_x,nproc_y,       &
                            nprocx_in,nprocy_in,max_fopen,istatus)

  IF ( IAMROOT .AND. UNTIN /= STDIN ) THEN
    CLOSE(UNTIN)
    CALL retunit( UNTIN )
  END IF

  CALL mpinit_var
  CALL wrf_parallel_init(myproc,nprocs,nproc_x,nproc_y)
  CALL wrf_parallel_set_halo_width(0,istatus)

!-----------------------------------------------------------------------
!
! Now the main processing work for ARW or NMM. Loop over all domains
! (note that NMM should only contain 1 domain)
!
!-----------------------------------------------------------------------

  DO n = 1,max_dom
    CALL process_domain(UNTOUT,n,IAMROOT,nproc_x,nproc_y,               &
                        nprocx_in,nprocy_in,istatus)
  END DO

!-----------------------------------------------------------------------
!
! Terminate the program nicely
!
!-----------------------------------------------------------------------

  IF (mp_opt > 0) THEN
    IF ( IAMROOT ) WRITE(UNTOUT,'(/,5x,a)') '==== ARPS4WRF_MPI terminated normally ===='
    CALL mpexit(0)
  ELSE
    WRITE(UNTOUT,'(/,5x,a)') '==== ARPS4WRF terminated normally ===='
    STOP
  END IF

END PROGRAM arps4wrf
