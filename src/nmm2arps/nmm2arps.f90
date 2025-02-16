PROGRAM nmm2arps
!
!##################################################################
!##################################################################
!######                                                      ######
!######                  PROGRAM nmm2arps                    ######
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
!   Interpolation WRF NMM history file to the ARPS grid.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yunheng Wang
!  12/18/2008
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
  USE module_nmm2arps_namelist
  USE wrf_parallel_module
  USE module_nmm2arps_nmmgrid
  USE module_output_arpsgrid
  USE module_interpolation

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
  INTEGER :: nargs, n

  CHARACTER(LEN=256) :: NML_FILE

  INTEGER :: nprocx_in, nprocy_in

  INTEGER :: fake_width

!-----------------------------------------------------------------------
!
! Wroking arrays
!
!-----------------------------------------------------------------------


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

  ! Non-MPI defaults: All others initialized in mpinit_var
  mp_opt = 0
  myproc = 0
  nprocx_in  = 1
  nprocy_in  = 1
  dumpstride = 1
  readstride = 1

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
      '####              Welcome to nmm2arps                      ####',&
      '####                                                       ####',&
      '####           Interpolation NMM runs to the ARPS grid     ####',&
      '####                                                       ####',&
      '###############################################################',&
      '######################################4########################'

  IF ( IAMROOT ) THEN
    nargs = COMMAND_ARGUMENT_COUNT()
    IF (nargs == 1) THEN
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

  CALL read_namelist_params(UNTIN,UNTOUT,IAMROOT,mp_opt,nproc_x,nproc_y,&
                            nprocx_in,nprocy_in,max_fopen,istatus)

  IF ( IAMROOT .AND. UNTIN /= STDIN ) THEN
    CLOSE(UNTIN)
    CALL retunit( UNTIN )
  END IF

  IF (.NOT. IAMROOT) lvldbg = 0

  CALL mpinit_var
  CALL wrf_parallel_init(myproc,nprocs,nproc_x,nproc_y)

!-----------------------------------------------------------------------
!
! Probe first WRF file and establish ARPS grid
!
!-----------------------------------------------------------------------
  IF (lvldbg > 2) WRITE(STDOUT,'(1x,a)') '=== Initialize NMMGRID ...'

  CALL nmmgrid_init(IAMROOT,wrffiles(1),io_form,                        &
                    multifile,ncompressx, ncompressy,magnitude_processor, &
                    lvldbg,istatus)
       ! initialize nmmgrid based on meta data from the first file

  IF (lvldbg > 2) WRITE(STDOUT,'(1x,a)') '=== Allocating ARPSGRID ...'

  CALL arpsgrid_meta(mapproj,trulat1,trulat2,trulon,sclfct,             &
                     nx,ny,nz,nstyp,dx,dy,dz,ctrlat,ctrlon,             &
                     strhopt,dzmin,zrefsfc,dlayer1,dlayer2,             &
                     strhtune,zflat,nzsoil,dzsoil,soilmodel_option,     &
                     terndta, ternfmt, runname, dirname,                &
                     nocmnt,cmnt, nmmgrid%num_moist,nmmgrid%num_qq,     &
                 nmmgrid%p_qc,nmmgrid%p_qr,nmmgrid%p_qi,                &
                 nmmgrid%p_qs,nmmgrid%p_qg,nmmgrid%p_qh,                &
                 nmmgrid%p_nc,nmmgrid%p_nr,nmmgrid%p_ni,                &
                 nmmgrid%p_ns,nmmgrid%p_ng,nmmgrid%p_nh,                &
                              nmmgrid%p_zr,nmmgrid%p_zi,                &
                 nmmgrid%p_zs,nmmgrid%p_zg,nmmgrid%p_zh,                &
                 nmmgrid%p_nn,nmmgrid%p_rim,                            &
                     dmp_out_joined, terndmp, soildmp, hdmpfmt,         &
                     basout,grdout,varout,iceout,sfcout,                &
                     landout,radout, tkeout, flxout,                    &
                     mstout, rainout, prcout, trbout,                   &
                     exbcdmp,qcexout,qrexout,qiexout,qsexout,qhexout,   &
                     hdfcompr,filcmprs,readyfl,                         &
                     init_time_str,initsec, abstimei, abstimee,         &
                     lvldbg, istatus)

  CALL arpsgrid_alloc_init(wrfexttrnopt,nx,ny,nz,istatus)

  fake_width = 0

  CALL find_nmm_halo_width(nmmgrid%ips,nmmgrid%ipe,nmmgrid%jps,nmmgrid%jpe, &
                           nx,ny,arpsgrid%xlat,arpsgrid%xlon,           &
                           arpsgrid%ylat,arpsgrid%ylon,                 &
                           arpsgrid%slat,arpsgrid%slon,                 &
                           myproc,lvldbg,fake_width,istatus)

  IF (lvldbg > 2) WRITE(STDOUT,'(1x,a,I4)') '--- NMMGRID Halo Width is set to ',fake_width

  IF (istatus /= 0) CALL arpsstop('ERROR: NMM Grid fake zone.',1)

  CALL wrf_parallel_set_halo_width(fake_width,istatus)

  IF (lvldbg > 2) WRITE(STDOUT,'(1x,a)') '=== Allocating NMMGRID ...'

  CALL nmmgrid_alloc(IAMROOT,lvldbg,istatus)

!-----------------------------------------------------------------------
!
! Now loop all time series
!
!-----------------------------------------------------------------------

  DO n = 1,nwrfdfil
    IF (lvldbg > 0) WRITE(STDOUT,'(1x,a,2(I4,a))')                      &
                    '=== Process: ',myproc, ' processing file - ',n,' ...'
    CALL process_file(UNTOUT,n,IAMROOT,istatus)
    IF (istatus /= 0) EXIT
  END DO

  CALL nmmgrid_dealloc( istatus )

  CALL arpsgrid_dealloc( istatus )

  CALL deallocate_interpolation_arrays( istatus )

!-----------------------------------------------------------------------
!
! Terminate the program nicely
!
!-----------------------------------------------------------------------

  IF (mp_opt > 0) THEN
    IF ( IAMROOT ) WRITE(UNTOUT,'(a/a,i4,a)')                          &
          ' ==== Normal succesful completion of NMM2ARPS_MPI. ====',   &
          '      Processed ',nwrfdfil,' file(s)'
    CALL mpexit(0)
  ELSE
    WRITE(UNTOUT,'(a/a,i4,a)')                                         &
          ' ==== Normal succesful completion of NMM2ARPS. ====',       &
          '      Processed ',nwrfdfil,' file(s)'
    STOP
  END IF

END PROGRAM nmm2arps
