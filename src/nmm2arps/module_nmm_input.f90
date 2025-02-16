MODULE module_nmm_input
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! This module handle WRF/NMM IO
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  USE module_wrfgrid_constants
  USE module_arpsgrid_constants
  USE wrf_parallel_module

  INTEGER, ALLOCATABLE :: fHndl(:,:)
  INTEGER, ALLOCATABLE :: ipstart(:,:), ipend(:,:), jpstart(:,:), jpend(:,:)
  INTEGER              :: nxsubf, nysubf
  LOGICAL              :: multi_file
  INTEGER              :: io_form_input
  INTEGER              :: numdigits

  LOGICAL              :: IAMIONODE

  CHARACTER(LEN=80)    :: sysdepinfo

  LOGICAL              :: input_initialized = .FALSE.

  INTEGER              :: dmrsize, dmisize
  REAL,    ALLOCATABLE :: temdmr(:)
  INTEGER, ALLOCATABLE :: temdmi(:)

  CONTAINS

  SUBROUTINE nmm_input_init(IAMROOT,io_form_in,numdigits_in,            &
                            multi_filein,ncompressx,ncompressy,istatus)

!#######################################################################
    IMPLICIT NONE

    LOGICAL, INTENT(IN)  :: IAMROOT
    LOGICAL, INTENT(IN)  :: multi_filein
    INTEGER, INTENT(IN)  :: ncompressx, ncompressy, numdigits_in
    INTEGER, INTENT(IN)  :: io_form_in

    INTEGER, INTENT(OUT) :: istatus
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    istatus = 0

    multi_file      = multi_filein
    nxsubf          = ncompressx
    nysubf          = ncompressy
    numdigits       = numdigits_in
    io_form_input   = io_form_in
    IAMIONODE       = IAMROOT

    IF (multi_file) THEN
      IAMIONODE = .TRUE.
      ALLOCATE(fHndl  (nxsubf,nysubf), STAT = istatus)
      ALLOCATE(ipstart(nxsubf,nysubf), STAT = istatus)
      ALLOCATE(ipend  (nxsubf,nysubf), STAT = istatus)
      ALLOCATE(jpstart(nxsubf,nysubf), STAT = istatus)
      ALLOCATE(jpend  (nxsubf,nysubf), STAT = istatus)
    ELSE
      ALLOCATE(fHndl(1,1), STAT = istatus)
    END IF

    sysdepinfo = 'DATASET=HISTORY'

    IF (io_form_input == BINARY) THEN
      CALL ext_int_ioinit( SysDepInfo, iStatus )
    ELSE IF (io_form_input == NETCDF) THEN
      CALL ext_ncd_ioinit( SysDepInfo, iStatus )
    ELSE
      WRITE(6,'(1x,a,I4,a)') 'Unsupported IO format - ',io_form_input,'.'
      CAlL arpsstop('Unsupported IO format.',1)
    END IF

    dmrsize = 1
    dmisize = 1
    ALLOCATE(temdmi(dmrsize), STAT = istatus)
    ALLOCATE(temdmr(dmisize), STAT = istatus)

    input_initialized = .TRUE.

    RETURN
  END SUBROUTINE nmm_input_init

  SUBROUTINE nmm_input_open(filename,dbglvl,istatus)

!#######################################################################

    IMPLICIT NONE

    CHARACTER(LEN=MAXFILELEN) :: filename
    INTEGER, INTENT(IN)       :: dbglvl
    INTEGER, INTENT(OUT)      :: istatus

!-----------------------------------------------------------------------

    LOGICAL :: fexists

    CHARACTER(LEN=MAXFILELEN), ALLOCATABLE :: tmpstr(:,:)

    INTEGER :: nprocx_in
    INTEGER :: iloc_x, jloc_y

    INTEGER :: loc_proc
    INTEGER :: iloc, jloc

    CHARACTER(LEN=20) :: fmtstr

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    istatus = 0

    IF (.NOT. input_initialized) THEN
      WRITE(6,'(1x,a)') 'Module nmm_input is still not initialized.'
      istatus = -1
      RETURN
    END IF

    IF (multi_file) THEN

      nprocx_in = nxsubf*nxprocs
      iloc_x    = my_x*nxsubf
      jloc_y    = my_y*nysubf

      WRITE(fmtstr,'(a,2(I1,a))') '(2a,I',numdigits,'.',numdigits,')'

      ALLOCATE(tmpstr(nxsubf,nysubf), STAT = istatus)
      DO jloc = 1, nysubf
        DO iloc = 1, nxsubf
          loc_proc = (jloc_y + jloc-1)*nprocx_in + iloc_x + (iloc-1)
          WRITE(tmpstr(iloc,jloc),FMT=fmtstr) TRIM(filename),'_',loc_proc

          INQUIRE(FILE = TRIM(tmpstr(iloc,jloc)), EXIST = fexists)
          IF ( .NOT. fexists ) THEN
            WRITE(6,'(3a)') 'File not found: ',tmpstr(iloc,jloc),' in nmm_input_open'
            istatus = -1
            EXIT
          END IF
        END DO
      END DO
      CALL mpmini(istatus)
      IF (istatus < 0) CALL arpsstop('WRF file not exist.',1)

      DO jloc = 1, nysubf
        DO iloc = 1, nxsubf
          IF (io_form_input == BINARY) THEN            ! initialize explicitly

            CALL ext_int_open_for_read(tmpstr(iloc,jloc), 0, 0,         &
                                 SysDepInfo, fHndl(iloc,jloc), iStatus )

          ELSE IF (io_form_input == NETCDF) THEN       ! no initialization needed

            CALL ext_ncd_open_for_read(tmpstr(iloc,jloc),0,0,           &
                                 SysDepInfo, fHndl(iloc,jloc), iStatus )

          END IF

        END DO
      END DO

      CALL mpmini(istatus)

    ELSE

      INQUIRE(FILE = TRIM(filename), EXIST = fexists)
      IF ( .NOT. fexists ) THEN
        WRITE(6,'(3a)') 'File not found: ',TRIM(filename),' in nmm_input_open.'
        CALL arpsstop('WRF file not exist.',1)
      END IF

      IF (dbglvl > 1) WRITE(6,'(3x,2a)') 'Opening file - ',TRIM(filename)

      IF (io_form_input == BINARY) THEN            ! initialize explicitly

        IF (IAMIONODE)  CALL ext_int_open_for_read(filename, 0, 0,       &
                                        SysDepInfo, fHndl(1,1), iStatus )

      ELSE IF (io_form_input == NETCDF) THEN       ! no initialization needed

        IF (IAMIONODE) CALL ext_ncd_open_for_read(filename,0,0,          &
                                        SysDepInfo, fHndl(1,1), iStatus)

      END IF

    END IF

    IF (istatus /= 0) THEN
      WRITE(6,'(1x,2a)') 'ERROR: Opening file ',TRIM(filename)
      !CALL arpsstop('Open WRF file error.',1)
      istatus = -2
      RETURN
    END IF

    RETURN
  END SUBROUTINE nmm_input_open

  SUBROUTINE nmm_input_close(istatus)
!#######################################################################
    IMPLICIT NONE
    INTEGER, INTENT(OUT)      :: istatus

!-----------------------------------------------------------------------
    INTEGER :: i,j

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    istatus = 0

    DO j = 1,nysubf
      DO i = 1,nxsubf
        IF (io_form_input == BINARY) THEN
          IF (IAMIONODE) CALL ext_int_ioclose(fHndl(i,j),iStatus)
        ELSE IF(io_form_input == NETCDF) THEN
          IF (IAMIONODE) CALL ext_ncd_ioclose(fHndl(i,j),iStatus)
        END IF
      END DO
    END DO

    IF (istatus /= 0) THEN
      WRITE(0,'(1x,a)') 'ERROR: closing file handler.'
      CALL arpsstop('Error in nmm_input_close.',1)
    END IF

    RETURN
  END SUBROUTINE nmm_input_close

  SUBROUTINE nmm_input_shutdown(istatus)
!#######################################################################
    IMPLICIT NONE
    INTEGER, INTENT(OUT) :: istatus

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    istatus = 0

    IF (ALLOCATED(temdmr)) DEALLOCATE(temdmr)
    IF (ALLOCATED(temdmi)) DEALLOCATE(temdmi)

    IF (io_form_input == BINARY) THEN
      CALL ext_int_ioexit( iStatus )
    ELSE IF (io_form_input == NETCDF) THEN
      CALL ext_ncd_ioexit( iStatus )
    END IF

    RETURN
  END SUBROUTINE nmm_input_shutdown

  SUBROUTINE nmm_input_get_meta(nx,ny,nz,nzsoil,mp_physics,             &
                            iproj,ctrlat,ctrlon,dx,dy,                  &
                            ips,ipe,jps,jpe,                            &
                            iswater, isice, isurban, isoilwater,istatus)
!#######################################################################
    IMPLICIT NONE
    INTEGER, INTENT(OUT) :: nx,ny,nz,nzsoil
    INTEGER, INTENT(OUT) :: mp_physics
    INTEGER, INTENT(OUT) :: iproj
    REAL,    INTENT(OUT) :: ctrlat,ctrlon,dx,dy
    INTEGER, INTENT(OUT) :: ips,ipe,jps,jpe
    INTEGER, INTENT(OUT) :: iswater, isice, isurban, isoilwater
    INTEGER, INTENT(OUT) :: istatus

!-----------------------------------------------------------------------
    INTEGER :: sfcphys

    INTEGER :: iloc, jloc

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    istatus = 0

    CALL nmm_input_get_dom_ti_integer('MAP_PROJ',  iproj,  istatus)
    CALL nmm_input_get_dom_ti_real   ('CEN_LAT',   ctrlat, istatus)
    CALL nmm_input_get_dom_ti_real   ('CEN_LON',   ctrlon, istatus)
    CALL nmm_input_get_dom_ti_real   ('DX',        dx,     istatus)
    CALL nmm_input_get_dom_ti_real   ('DY',        dy,     istatus)

    CALL nmm_input_get_dom_ti_integer('MP_PHYSICS',  mp_physics,  istatus)

    CALL nmm_input_get_dom_ti_integer('ISWATER',    iswater,   istatus)
    CALL nmm_input_get_dom_ti_integer('ISICE',      isice,     istatus)
    CALL nmm_input_get_dom_ti_integer('ISURBAN',    isurban,   istatus)
    CALL nmm_input_get_dom_ti_integer('ISOILWATER', isoilwater,istatus)

    CALL nmm_input_get_dom_ti_integer( 'WEST-EAST_GRID_DIMENSION',  nx, istatus)
    CALL nmm_input_get_dom_ti_integer( 'SOUTH-NORTH_GRID_DIMENSION',ny, istatus)
    CALL nmm_input_get_dom_ti_integer( 'BOTTOM-TOP_GRID_DIMENSION', nz, istatus)

    CALL nmm_input_get_dom_ti_integer( 'SF_SURFACE_PHYSICS',sfcphys,istatus)

    !
    ! Processor special
    !
    IF (multi_file) THEN
      DO jloc = 1, nysubf
        DO iloc = 1,nxsubf
          CALL nmm_get_dom_ti_integer( fHndl(iloc,jloc), 'WEST-EAST_PATCH_START_UNSTAG',  &
                                       ipstart(iloc,jloc), istatus)
          CALL nmm_get_dom_ti_integer( fHndl(iloc,jloc), 'WEST-EAST_PATCH_END_UNSTAG',    &
                                       ipend(iloc,jloc), istatus)
          CALL nmm_get_dom_ti_integer( fHndl(iloc,jloc), 'SOUTH-NORTH_PATCH_START_UNSTAG',&
                                       jpstart(iloc,jloc), istatus)
          CALL nmm_get_dom_ti_integer( fHndl(iloc,jloc), 'SOUTH-NORTH_PATCH_END_UNSTAG',  &
                                       jpend(iloc,jloc), istatus)
        END DO
      END DO
      ips = ipstart(1,1)
      ipe = ipend(nxsubf,1)
      jps = jpstart(1,1)
      jpe = jpend(1,nysubf)
    ELSE
      CALL nmm_input_get_dom_ti_integer('WEST-EAST_PATCH_START_UNSTAG',   ips, istatus)
      CALL nmm_input_get_dom_ti_integer('WEST-EAST_PATCH_END_UNSTAG',     ipe, istatus)
      CALL nmm_input_get_dom_ti_integer('SOUTH-NORTH_PATCH_START_UNSTAG', jps, istatus)
      CALL nmm_input_get_dom_ti_integer('SOUTH-NORTH_PATCH_END_UNSTAG',   jpe, istatus)
    END IF

!    CALL nmm_input_get_dom_ti_integer( 'BOTTOM-TOP_PATCH_START_STAG', kps, istatus)
!    CALL nmm_input_get_dom_ti_integer( 'BOTTOM-TOP_PATCH_END_STAG',   kpe, istatus)

!-----------------------------------------------------------------------
!
!  Determine soil layers from surface physics option
!
!-----------------------------------------------------------------------

    IF (sfcphys == 1) THEN
      nzsoil = 5
    ELSE IF (sfcphys == 2 .OR. sfcphys == 99) THEN
      nzsoil = 4
    ELSE IF (sfcphys == 3) THEN
      nzsoil = 6
    ELSE
      WRITE(6,*) '=============================================='
      WRITE(6,*) 'WARNING: unknown sf_surface_physics = ',sfcphys
      WRITE(6,*) '=============================================='
      nzsoil = 5
    END IF

    RETURN
  END SUBROUTINE nmm_input_get_meta

  SUBROUTINE nmm_input_get_dom_ti_integer( element, val, ireturn)
    IMPLICIT NONE

    CHARACTER(*), INTENT(IN)  :: element
    INTEGER,      INTENT(OUT) :: val

    INTEGER,      INTENT(OUT) :: ireturn

!-----------------------------------------------------------------------

    INTEGER :: nid
    INTEGER    :: outcount

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Begin of executable code
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    nid = fHndl(1,1)
    IF (IAMIONODE) THEN
      IF (io_form_input == NETCDF) THEN
         CALL ext_ncd_get_dom_ti_integer(nid,element,val,1,outcount,ireturn)
      ELSE IF (io_form_input == 1) THEN
         CALL ext_int_get_dom_ti_integer(nid,element,val,1,outcount,ireturn)
      END IF
    END IF
    CALL mpupdatei(val,1)

    RETURN
  END SUBROUTINE nmm_input_get_dom_ti_integer

  SUBROUTINE nmm_input_get_dom_ti_real( element, val, ireturn)
    IMPLICIT NONE

    CHARACTER(*), INTENT(IN)  :: element
    REAL,         INTENT(OUT) :: val
    INTEGER,      INTENT(OUT) :: ireturn

!-----------------------------------------------------------------------

    INTEGER   :: nid
    INTEGER   :: outcount

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Begin of executable code
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    nid = fHndl(1,1)
    IF (IAMIONODE) THEN
      IF (io_form_input == NETCDF) THEN
         CALL ext_ncd_get_dom_ti_real(nid,element,val,1,outcount,ireturn)
      ELSE IF (io_form_input == BINARY) THEN
         CALL ext_int_get_dom_ti_real(nid,element,val,1, outcount,ireturn)
      END IF
    END IF
    CALL mpupdater(val,1)

    RETURN
  END SUBROUTINE nmm_input_get_dom_ti_real
  !
  SUBROUTINE nmm_input_get_dom_ti_char( element, val, ireturn)
    IMPLICIT NONE

    CHARACTER(*), INTENT(IN)  :: element
    CHARACTER(*), INTENT(OUT) :: val
    INTEGER,      INTENT(OUT) :: ireturn

!-----------------------------------------------------------------------

    INTEGER   :: nid
    INTEGER   :: outcount

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Begin of executable code
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    nid = fHndl(1,1)
    IF (IAMIONODE) THEN
      IF (io_form_input == NETCDF) THEN
         CALL ext_ncd_get_dom_ti_char(nid,element,val,ireturn)
      ELSE IF (io_form_input == 1) THEN
         CALL ext_int_get_dom_ti_char(nid,element,val,ireturn)
      END IF
    END IF
    CALL mpupdatec(val,LEN(val))

    RETURN
  END SUBROUTINE nmm_input_get_dom_ti_char

  SUBROUTINE nmm_input_get_3dvar(datestr,varname,stagger,dimnames,      &
                        MemoryStart,MemoryEnd,PatchStart,PatchEnd,      &
                        DomainStart,DomainEnd,var3d,                    &
                        dbglvl, istatus)

!-----------------------------------------------------------------------
!
!  PURPOSE:
!    Read in a 3D array from the WRF NetCDF file. Only root processor do
!    the reads, and broadcast to all other processors if necessary.
!
!
!-----------------------------------------------------------------------
    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN)  :: datestr
    CHARACTER(LEN=*), INTENT(IN)  :: varname
    CHARACTER(LEN=*), INTENT(IN)  :: stagger

    CHARACTER(LEN=*), DIMENSION(MAXDIMENSIONS), INTENT(IN) :: dimnames
    INTEGER,          DIMENSION(MAXDIMENSIONS), INTENT(IN) :: MemoryStart
    INTEGER,          DIMENSION(MAXDIMENSIONS), INTENT(IN) :: MemoryEnd
    INTEGER,          DIMENSION(MAXDIMENSIONS), INTENT(IN) :: PatchStart
    INTEGER,          DIMENSION(MAXDIMENSIONS), INTENT(IN) :: PatchEnd
    INTEGER,          DIMENSION(MAXDIMENSIONS), INTENT(IN) :: DomainStart
    INTEGER,          DIMENSION(MAXDIMENSIONS), INTENT(IN) :: DomainEnd

    REAL,    INTENT(OUT) :: var3d(MemoryStart(1):MemoryEnd(1),          &
                                  MemoryStart(2):MemoryEnd(2),          &
                                  MemoryStart(3):MemoryEnd(3))

    INTEGER, INTENT(IN)  :: dbglvl
    INTEGER, INTENT(OUT) :: istatus

!-----------------------------------------------------------------------
!
! Misc. Local variables
!
!-----------------------------------------------------------------------
!
    INTEGER  :: iloc, jloc
    INTEGER  :: msloc(MAXDIMENSIONS),meloc(MAXDIMENSIONS)
    INTEGER  :: i, j, k
    INTEGER  :: nfid

    INTEGER  :: nx, ny, nz, kin, jin, kjiin

    CHARACTER(LEN=3)  :: mOrder
    CHARACTER(LEN=20) :: fmtstr

    INTEGER  :: mInsize

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    istatus = 0

    IF (dbglvl > 0) &
    WRITE(6,FMT='(3x,2a)',ADVANCE='NO') 'Reading 3D variable ', TRIM(varname)

    IF ( DomainEnd(3) /= DomainStart(3) ) THEN
      mOrder = 'XYZ'
    ELSE IF ( DomainEnd(2) /= DomainStart(2) ) THEN
      mOrder = 'XY '
    ELSE IF ( DomainEnd(1) /= DomainStart(1) ) THEN
      mOrder = 'Z  '
    ELSE
      mOrder = '0  '
    END IF

    IF ( multi_file .AND. mOrder(1:2) == 'XY' ) THEN

      DO jloc = 1, nysubf
        DO iloc = 1, nxsubf
          msloc(1) = ipstart(iloc,jloc); meloc(1) = ipend(iloc,jloc)
          msloc(2) = jpstart(iloc,jloc); meloc(2) = jpend(iloc,jloc)
          msloc(3) = patchStart(3);      meloc(3) = patchEnd(3)

          nx  = ipend(iloc,jloc) - ipstart(iloc,jloc) + 1
          ny  = jpend(iloc,jloc) - jpstart(iloc,jloc) + 1
          nz  = patchEnd(3) - patchStart(3) + 1

          mInsize = nx*ny*nz

          IF (mInsize > dmrsize) THEN
          	DEALLOCATE ( temdmr )
          	dmrsize = mInsize
          	ALLOCATE (temdmr(dmrsize), STAT = istatus)
          END IF

          IF (io_form_input == NETCDF) THEN
            CALL ext_ncd_read_field(fHndl(iloc,jloc), DateStr,VarName,  &
                                    temdmr,WRF_REAL,                    &
                                    0,0,1, mOrder, Stagger, DimNames,   &
                                    mSloc,mEloc,                        &
                                    mSloc,mEloc,                        &
                                    mSloc,mEloc,                        &
                                    iStatus)

          ELSE IF (io_form_input == BINARY) THEN

            CALL ext_int_read_field(fHndl(iloc,jloc), DateStr,VarName,  &
                                    temdmr,WRF_REAL,                    &
                                    0, 0, 1, mOrder, Stagger, DimNames, &
                                    DomainStart,DomainEnd,              &
                                    mSloc,mEloc,                        &
                                    mSloc,mEloc,                        &
                                    iStatus)
          END IF

          ! To merge
          DO k = patchStart(3), patchEnd(3)
            kin = (k - patchStart(3))*nx*ny
        	  DO j = jpstart(iloc,jloc), jpend(iloc,jloc)
        	    jin = (j - jpstart(iloc,jloc))*nx
        	    DO i = ipstart(iloc,jloc), ipend(iloc,jloc)
        	      kjiin = kin + jin + i - ipstart(iloc,jloc) + 1    ! relative position in memory
        	      var3d(i,j,k) = temdmr(kjiin)
        	    END DO
        	  END DO
          END DO

        END DO
      END DO

    ELSE

      nx = domainEnd(1) - domainStart(1) + 1
      ny = domainEnd(2) - domainStart(2) + 1
      nz = domainEnd(3) - domainStart(3) + 1
      mInsize = nx*ny*nz

      IF (mInsize > dmrsize) THEN
        DEALLOCATE ( temdmr )
      	dmrsize = mInsize
        ALLOCATE ( temdmr(dmrsize), STAT = istatus)
      END IF

      nfid = fHndl(1,1)
      IF (IAMIONODE) THEN
        IF (io_form_input == NETCDF) THEN
          CALL ext_ncd_read_field(nfid, DateStr,VarName,temdmr,WRF_REAL,&
                                  0,0,1, mOrder, Stagger, DimNames,     &
                                  DomainStart,DomainEnd,                &
                                  DomainStart,DomainEnd,                &
                                  DomainStart,DomainEnd,                &
                                  iStatus)

        ELSE IF (io_form_input == BINARY) THEN

          CALL ext_int_read_field(nfid, DateStr,VarName,temdmr,WRF_REAL,&
                                  0, 0, 1, mOrder, Stagger, DimNames,   &
                                  DomainStart,DomainEnd,                &
                                  DomainStart,DomainEnd,                &
                                  DomainStart,DomainEnd,                &
                                  iStatus)
        END IF
      END IF  ! IAMIONODE
      CALL parallel_bcast_int(iStatus)

      IF (iStatus == 0) THEN    ! successful read

        IF (ntprocs > 1 .AND. mOrder(1:2) == 'XY') THEN   ! To scatter
            CALL scatter_whole_field_r(var3d,                           &
                                     memoryStart(1), memoryEnd(1),      &
                                     memoryStart(2), memoryEnd(2),      &
                                     memoryStart(3), memoryEnd(3),      &
                                     patchStart(1),  patchEnd(1),       &
                                     patchStart(2),  patchEnd(2),       &
                                     patchStart(3),  patchEnd(3),       &
                                     temdmr,                            &
                                     domainstart(1), domainend(1),      &
                                     domainstart(2), domainend(2),      &
                                     domainstart(3), domainend(3) )

        ELSE
          IF (ntprocs > 1) CALL parallel_bcast_real_array(temdmr,mInsize)

          DO k = patchStart(3), patchEnd(3)
            kin = ( k - patchStart(3) )*ny*nx  ! patch and domain should be the same
        	  DO j = patchStart(2), patchEnd(2)
        	    jin = ( j - patchStart(2) )*nx
        	    DO i = patchStart(1), patchEnd(1)
        	      kjiin = kin + jin + i - patchStart(1) + 1
        	  	  var3d(i,j,k) = temdmr(kjiin)
        	    END DO
        	  END DO
    	    END DO

        END IF
      END IF
    END IF

    WRITE(fmtstr,'(a,I2.2,a)') '(',10-LEN_TRIM(varname),'x,a)'

    IF (istatus == 0) THEN

      IF (mOrder(1:2) == 'XY') THEN
        CALL exchange_halo_r(var3d,                                     &
                           memoryStart(1), memoryEnd(1),                &
                           memoryStart(2), memoryEnd(2),                &
                           memoryStart(3), memoryEnd(3),                &
                           patchStart(1),  patchEnd(1),                 &
                           patchStart(2),  patchEnd(2),                 &
                           patchStart(3),  patchEnd(3) )
      END IF

      IF (dbglvl > 0) WRITE(6,FMT=fmtstr) ' ...  DONE.'
    ELSE
      IF (dbglvl > 0) WRITE(6,FMT=fmtstr) ' ...  ERROR.'
    END IF

    RETURN
  END SUBROUTINE nmm_input_get_3dvar

  SUBROUTINE nmm_input_get_3dvari(datestr,varname,stagger,dimnames,      &
                        MemoryStart,MemoryEnd,PatchStart,PatchEnd,      &
                        DomainStart,DomainEnd,var3d,                    &
                        dbglvl, istatus)

!-----------------------------------------------------------------------
!
!  PURPOSE:
!    Read in a 3D array from the WRF NetCDF file. Only root processor do
!    the reads, and broadcast to all other processors if necessary.
!
!
!-----------------------------------------------------------------------
    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN)  :: datestr
    CHARACTER(LEN=*), INTENT(IN)  :: varname
    CHARACTER(LEN=*), INTENT(IN)  :: stagger

    CHARACTER(LEN=*), DIMENSION(MAXDIMENSIONS), INTENT(IN) :: dimnames
    INTEGER,          DIMENSION(MAXDIMENSIONS), INTENT(IN) :: MemoryStart
    INTEGER,          DIMENSION(MAXDIMENSIONS), INTENT(IN) :: MemoryEnd
    INTEGER,          DIMENSION(MAXDIMENSIONS), INTENT(IN) :: PatchStart
    INTEGER,          DIMENSION(MAXDIMENSIONS), INTENT(IN) :: PatchEnd
    INTEGER,          DIMENSION(MAXDIMENSIONS), INTENT(IN) :: DomainStart
    INTEGER,          DIMENSION(MAXDIMENSIONS), INTENT(IN) :: DomainEnd

    INTEGER, INTENT(OUT) :: var3d(MemoryStart(1):MemoryEnd(1),          &
                                  MemoryStart(2):MemoryEnd(2),          &
                                  MemoryStart(3):MemoryEnd(3))

    INTEGER, INTENT(IN)  :: dbglvl
    INTEGER, INTENT(OUT) :: istatus

!-----------------------------------------------------------------------
!
! Misc. Local variables
!
!-----------------------------------------------------------------------
!
    INTEGER  :: iloc, jloc
    INTEGER  :: msloc(MAXDIMENSIONS),meloc(MAXDIMENSIONS)
    INTEGER  :: i, j, k
    INTEGER  :: nfid

    INTEGER  :: nx, ny, nz, kin, jin, kjiin

    CHARACTER(LEN=3)  :: mOrder
    CHARACTER(LEN=20) :: fmtstr

    INTEGER  :: mInsize

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    istatus = 0

    IF (dbglvl > 0) &
    WRITE(6,FMT='(3x,2a)',ADVANCE='NO') 'Reading 3D variable ', TRIM(varname)

    IF ( DomainEnd(3) /= DomainStart(3) ) THEN
      mOrder = 'XYZ'
    ELSE IF ( DomainEnd(2) /= DomainStart(2) ) THEN
      mOrder = 'XY '
    ELSE IF ( DomainEnd(1) /= DomainStart(1) ) THEN
      mOrder = 'Z  '
    ELSE
      mOrder = '0  '
    END IF

    IF ( multi_file .AND. mOrder(1:2) == 'XY' ) THEN

      DO jloc = 1, nysubf
        DO iloc = 1, nxsubf
          msloc(1) = ipstart(iloc,jloc); meloc(1) = ipend(iloc,jloc)
          msloc(2) = jpstart(iloc,jloc); meloc(2) = jpend(iloc,jloc)
          msloc(3) = patchStart(3);      meloc(3) = patchEnd(3)

          nx  = ipend(iloc,jloc) - ipstart(iloc,jloc) + 1
          ny  = jpend(iloc,jloc) - jpstart(iloc,jloc) + 1
          nz  = patchEnd(3) - patchStart(3) + 1

          mInsize = nx*ny*nz

          IF (mInsize > dmisize) THEN
          	DEALLOCATE (temdmi)
          	dmisize = mInsize
          	ALLOCATE (temdmi(dmisize), STAT = istatus)
          END IF

          IF (io_form_input == NETCDF) THEN
            CALL ext_ncd_read_field(fHndl(iloc,jloc), DateStr,VarName,  &
                                    temdmi,WRF_INTEGER,                 &
                                    0,0,1, mOrder, Stagger, DimNames,   &
                                    mSloc,mEloc,                        &
                                    mSloc,mEloc,                        &
                                    mSloc,mEloc,                        &
                                    iStatus)

          ELSE IF (io_form_input == BINARY) THEN

            CALL ext_int_read_field(fHndl(iloc,jloc), DateStr,VarName,  &
                                    temdmi,WRF_INTEGER,                 &
                                    0, 0, 1, mOrder, Stagger, DimNames, &
                                    mSloc,mEloc,                        &
                                    mSloc,mEloc,                        &
                                    mSloc,mEloc,                        &
                                    iStatus)
          END IF

          ! To merge
          DO k = patchStart(3), patchEnd(3)
            kin = (k - patchStart(3))*nx*ny
        	  DO j = jpstart(iloc,jloc), jpend(iloc,jloc)
        	    jin = (j - jpstart(iloc,jloc))*nx
        	    DO i = ipstart(iloc,jloc), ipend(iloc,jloc)
        	      kjiin = kin + jin + i - ipstart(iloc,jloc) + 1    ! relative position in memory
        	      var3d(i,j,k) = temdmi(kjiin)
        	    END DO
        	  END DO
          END DO

        END DO
      END DO

    ELSE

      nx = domainEnd(1) - domainStart(1) + 1
      ny = domainEnd(2) - domainStart(2) + 1
      nz = domainEnd(3) - domainStart(3) + 1
      mInsize = nx*ny*nz

      IF (mInsize > dmisize) THEN
      	DEALLOCATE (temdmi )
      	dmisize = mInsize
      	ALLOCATE (temdmi(dmisize), STAT = istatus)
      END IF

      nfid = fHndl(1,1)
      IF (IAMIONODE) THEN
        IF (io_form_input == NETCDF) THEN
          CALL ext_ncd_read_field(nfid, DateStr,VarName,temdmi,WRF_INTEGER,&
                                  0,0,1, mOrder, Stagger, DimNames,     &
                                  DomainStart,DomainEnd,                &
                                  DomainStart,DomainEnd,                &
                                  DomainStart,DomainEnd,                &
                                  iStatus)

        ELSE IF (io_form_input == BINARY) THEN

          CALL ext_int_read_field(nfid, DateStr,VarName,temdmi,WRF_INTEGER,&
                                  0, 0, 1, mOrder, Stagger, DimNames,   &
                                  DomainStart,DomainEnd,                &
                                  DomainStart,DomainEnd,                &
                                  DomainStart,DomainEnd,                &
                                  iStatus)
        END IF
      END IF  ! IAMIONODE
      CALL parallel_bcast_int(iStatus)

      IF (iStatus == 0) THEN    ! successful read

        IF (ntprocs > 1 .AND. mOrder(1:2) == 'XY') THEN   ! To scatter
            CALL scatter_whole_field_i(var3d,                           &
                                     memoryStart(1), memoryEnd(1),      &
                                     memoryStart(2), memoryEnd(2),      &
                                     memoryStart(3), memoryEnd(3),      &
                                     patchStart(1),  patchEnd(1),       &
                                     patchStart(2),  patchEnd(2),       &
                                     patchStart(3),  patchEnd(3),       &
                                     temdmi,                            &
                                     domainstart(1), domainend(1),      &
                                     domainstart(2), domainend(2),      &
                                     domainstart(3), domainend(3) )

        ELSE
          !IF (ntprocs > 1) CALL parallel_bcast_real_array(temdm,mInsize)

          DO k = patchStart(3), patchEnd(3)
            kin = ( k - patchStart(3) )*ny*nx  ! patch and domain should be the same
        	  DO j = patchStart(2), patchEnd(2)
        	    jin = ( j - patchStart(2) )*nx
        	    DO i = patchStart(1), patchEnd(1)
        	      kjiin = kin + jin + i - patchStart(1) + 1
        	  	  var3d(i,j,k) = temdmi(kjiin)
        	    END DO
        	  END DO
    	    END DO

        END IF
      END IF
    END IF

    WRITE(fmtstr,'(a,I2.2,a)') '(',10-LEN_TRIM(varname),'x,a)'

    IF (istatus == 0) THEN

      IF (mOrder(1:2) == 'XY') THEN
        CALL exchange_halo_i(var3d,                                     &
                           memoryStart(1), memoryEnd(1),                &
                           memoryStart(2), memoryEnd(2),                &
                           memoryStart(3), memoryEnd(3),                &
                           patchStart(1),  patchEnd(1),                 &
                           patchStart(2),  patchEnd(2),                 &
                           patchStart(3),  patchEnd(3) )
      END IF

      IF (dbglvl > 0) WRITE(6,FMT=fmtstr) ' ...  DONE.'
    ELSE
      IF (dbglvl > 0) WRITE(6,FMT=fmtstr) ' ...  ERROR.'
    END IF

    RETURN
  END SUBROUTINE nmm_input_get_3dvari

  SUBROUTINE nmm_get_rain(filename,datestr,varname1,varname2, stagger,  &
                          dimnames, MemoryStart, MemoryEnd,             &
                          PatchStart, PatchEnd, raing, rainc,           &
                          DomainStart,DomainEnd, dbglvl,istatus )

!#######################################################################
! This is a special subprogram that read precipiation variables,
! but does not interface file handles within this Module.
!
! Actually, it is a combination of three subroutines below:
!   nmm_input_open, nmm_input_close, nmm_input_get_3dvar.
!

    IMPLICIT NONE
    CHARACTER(LEN=MAXFILELEN) :: filename

    CHARACTER(LEN=*), INTENT(IN)  :: datestr
    CHARACTER(LEN=*), INTENT(IN)  :: varname1, varname2
    CHARACTER(LEN=*), INTENT(IN)  :: stagger

    CHARACTER(LEN=*), DIMENSION(MAXDIMENSIONS), INTENT(IN) :: dimnames
    INTEGER, DIMENSION(MAXDIMENSIONS), INTENT(IN) :: MemoryStart
    INTEGER, DIMENSION(MAXDIMENSIONS), INTENT(IN) :: MemoryEnd
    INTEGER, DIMENSION(MAXDIMENSIONS), INTENT(IN) :: PatchStart
    INTEGER, DIMENSION(MAXDIMENSIONS), INTENT(IN) :: PatchEnd
    REAL,    INTENT(OUT) :: raing(MemoryStart(1):MemoryEnd(1),          &
                                  MemoryStart(2):MemoryEnd(2),          &
                                  MemoryStart(3):MemoryEnd(3))
    REAL,    INTENT(OUT) :: rainc(MemoryStart(1):MemoryEnd(1),          &
                                  MemoryStart(2):MemoryEnd(2),          &
                                  MemoryStart(3):MemoryEnd(3))
    INTEGER, DIMENSION(MAXDIMENSIONS), INTENT(IN) :: DomainStart
    INTEGER, DIMENSION(MAXDIMENSIONS), INTENT(IN) :: DomainEnd

    INTEGER, INTENT(IN)  :: dbglvl
    INTEGER, INTENT(OUT) :: istatus

!-----------------------------------------------------------------------

    INTEGER :: fHndlr(nxsubf, nysubf)
    LOGICAL :: fexists
    INTEGER :: i, j, k

    INTEGER  :: nfid
    CHARACTER(LEN=3)  :: mOrder
    CHARACTER(LEN=20) :: fmtstr

    INTEGER :: nprocx_in
    INTEGER :: iloc_x, jloc_y, loc_proc
    CHARACTER(LEN=MAXFILELEN), ALLOCATABLE :: tmpstr(:,:)

    INTEGER  :: iloc, jloc
    INTEGER  :: msloc(MAXDIMENSIONS),meloc(MAXDIMENSIONS)

    INTEGER  :: nx, ny, nz, kin, jin, kjiin
    INTEGER  :: mInsize

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    istatus = 0

    IF (.NOT. input_initialized) THEN
      WRITE(6,'(1x,a)') 'Module nmm_input is still not initialized.'
      istatus = -1
      RETURN
    END IF

    IF (multi_file) THEN
      nprocx_in = nxsubf*nxprocs
      iloc_x    = my_x*nxsubf
      jloc_y    = my_y*nysubf

      WRITE(fmtstr,'(a,2(I1,a))') '(2a,I',numdigits,'.',numdigits,')'

      ALLOCATE(tmpstr(nxsubf,nysubf), STAT = istatus)
      DO jloc = 1, nysubf
        DO iloc = 1, nxsubf
          loc_proc = (jloc_y + jloc-1)*nprocx_in + iloc_x + (iloc-1)
          WRITE(tmpstr(iloc,jloc),FMT=fmtstr) filename,'_',loc_proc

          INQUIRE(FILE = TRIM(tmpstr(iloc,jloc)), EXIST = fexists)
          IF ( .NOT. fexists ) THEN
            WRITE(6,'(3a)') 'File not found: ',tmpstr(iloc,jloc),' in nmm_input_open'
            istatus = -1
            EXIT
          END IF
        END DO
      END DO
      CALL mpmini(istatus)
      IF (istatus < 0) CALL arpsstop('WRF file not exist.',1)

      DO jloc = 1, nysubf
        DO iloc = 1, nxsubf
          IF (io_form_input == BINARY) THEN            ! initialize explicitly

            CALL ext_int_open_for_read(tmpstr(iloc,jloc), 0, 0,         &
                                 SysDepInfo, fHndlr(iloc,jloc), iStatus )

          ELSE IF (io_form_input == NETCDF) THEN       ! no initialization needed

            CALL ext_ncd_open_for_read(tmpstr(iloc,jloc),0,0,           &
                                 SysDepInfo, fHndlr(iloc,jloc), iStatus )

          END IF

        END DO
      END DO

      CALL mpmini(istatus)

    ELSE

      INQUIRE(FILE = TRIM(filename), EXIST = fexists)
      IF ( .NOT. fexists ) THEN
        WRITE(6,'(3a)') 'File not found: ',TRIM(filename),' in nmm_get_rain.'
        istatus = -2
        RETURN
      END IF

      IF (dbglvl > 1 .AND. IAMIONODE) WRITE(6,'(1x,2a)') 'Opening file - ',TRIM(filename)

      IF (io_form_input == BINARY) THEN            ! initialize explicitly

        IF (IAMIONODE)  CALL ext_int_open_for_read(filename, 0, 0,       &
                                        SysDepInfo, fHndlr(1,1), iStatus )

      ELSE IF (io_form_input == NETCDF) THEN       ! no initialization needed

        IF (IAMIONODE) CALL ext_ncd_open_for_read(filename,0,0,          &
                                        SysDepInfo, fHndlr(1,1), iStatus)

      END IF

    END IF

    IF (istatus /= 0) THEN
      WRITE(0,'(1x,2a)') 'ERROR: Opening file ',TRIM(filename)
      istatus = -2
      RETURN
    END IF

!-------------------- Reading ------------------------------------------

    IF ( DomainEnd(3) /= DomainStart(3) ) THEN
      mOrder = 'XYZ'
    ELSE IF ( DomainEnd(2) /= DomainStart(2) ) THEN
      mOrder = 'XY '
    ELSE IF ( DomainEnd(1) /= DomainStart(1) ) THEN
      mOrder = 'Z  '
    ELSE
      mOrder = '0  '
    END IF

!-------------------- Reading RAING ------------------------------------

    IF (dbglvl > 0) &
    WRITE(6,FMT='(3x,2a)',ADVANCE='NO') 'Reading 3D variable ', TRIM(varname1)

    IF ( multi_file .AND. mOrder(1:2) == 'XY' ) THEN

      DO jloc = 1, nysubf
        DO iloc = 1, nxsubf
          msloc(1) = ipstart(iloc,jloc); meloc(1) = ipend(iloc,jloc)
          msloc(2) = jpstart(iloc,jloc); meloc(2) = jpend(iloc,jloc)
          msloc(3) = patchStart(3);      meloc(3) = patchEnd(3)

          nx  = ipend(iloc,jloc) - ipstart(iloc,jloc) + 1
          ny  = jpend(iloc,jloc) - jpstart(iloc,jloc) + 1
          nz  = patchEnd(3) - patchStart(3) + 1

          mInsize = nx*ny*nz

          IF (mInsize > dmrsize) THEN
          	DEALLOCATE ( temdmr )
          	dmrsize = mInsize
          	ALLOCATE (temdmr(dmrsize), STAT = istatus)
          END IF

          IF (io_form_input == NETCDF) THEN
            CALL ext_ncd_read_field(fHndlr(iloc,jloc), DateStr,VarName1,&
                     temdmr,WRF_REAL, 0,0,1, mOrder, Stagger, DimNames, &
                     mSloc,mEloc, mSloc,mEloc, mSloc,mEloc, iStatus)

          ELSE IF (io_form_input == BINARY) THEN

            CALL ext_int_read_field(fHndlr(iloc,jloc), DateStr,VarName1,&
                     temdmr,WRF_REAL, 0,0,1, mOrder, Stagger, DimNames, &
                     mSloc,mEloc, mSloc,mEloc, mSloc,mEloc, iStatus)
          END IF

          ! To merge
          DO k = patchStart(3), patchEnd(3)
            kin = (k - patchStart(3))*nx*ny
        	  DO j = jpstart(iloc,jloc), jpend(iloc,jloc)
        	    jin = (j - jpstart(iloc,jloc))*nx
        	    DO i = ipstart(iloc,jloc), ipend(iloc,jloc)
        	      kjiin = kin + jin + i - ipstart(iloc,jloc) + 1    ! relative position in memory
        	      raing(i,j,k) = temdmr(kjiin)
        	    END DO
        	  END DO
          END DO

        END DO
      END DO

    ELSE

      nx = domainEnd(1) - domainStart(1) + 1
      ny = domainEnd(2) - domainStart(2) + 1
      nz = domainEnd(3) - domainStart(3) + 1
      mInsize = nx*ny*nz

      IF (mInsize > dmrsize) THEN
        DEALLOCATE ( temdmr )
      	dmrsize = mInsize
        ALLOCATE ( temdmr(dmrsize), STAT = istatus)
      END IF

      nfid = fHndlr(1,1)
      IF (IAMIONODE) THEN
        IF (io_form_input == NETCDF) THEN
          CALL ext_ncd_read_field(nfid, DateStr,VarName1,temdmr,WRF_REAL,&
                                  0,0,1, mOrder, Stagger, DimNames,     &
                                  DomainStart,DomainEnd,                &
                                  DomainStart,DomainEnd,                &
                                  DomainStart,DomainEnd,                &
                                  iStatus)

        ELSE IF (io_form_input == BINARY) THEN

          CALL ext_int_read_field(nfid, DateStr,VarName1,temdmr,WRF_REAL,&
                                  0, 0, 1, mOrder, Stagger, DimNames,   &
                                  DomainStart,DomainEnd,                &
                                  DomainStart,DomainEnd,                &
                                  DomainStart,DomainEnd,                &
                                  iStatus)
        END IF
      END IF  ! IAMIONODE
      CALL parallel_bcast_int(iStatus)

      IF (iStatus == 0) THEN    ! successful read
        IF (ntprocs > 1 .AND. mOrder(1:2) == 'XY') THEN
            CALL scatter_whole_field_r(raing,                           &
                                     memoryStart(1), memoryEnd(1),      &
                                     memoryStart(2), memoryEnd(2),      &
                                     memoryStart(3), memoryEnd(3),      &
                                     patchStart(1),  patchEnd(1),       &
                                     patchStart(2),  patchEnd(2),       &
                                     patchStart(3),  patchEnd(3),       &
                                     temdmr,                            &
                                     domainstart(1), domainend(1),      &
                                     domainstart(2), domainend(2),      &
                                     domainstart(3), domainend(3) )
        ELSE
          IF (ntprocs > 1) CALL parallel_bcast_real_array(temdmr,mInsize)

          DO k = patchStart(3), patchEnd(3)
            kin = ( k - patchStart(3) )*ny*nx  ! patch and domain should be the same
            DO j = patchStart(2), patchEnd(2)
              jin = ( j - patchStart(2) )*nx
              DO i = patchStart(1), patchEnd(1)
                kjiin = kin + jin + i - patchStart(1) + 1
                raing(i,j,k) = temdmr(kjiin)
              END DO
            END DO
          END DO

        END IF
      END IF
    END IF

    WRITE(fmtstr,'(a,I2.2,a)') '(',10-LEN_TRIM(varname1),'x,a)'

    IF (istatus == 0) THEN
      IF (mOrder(1:2) == 'XY') THEN
        CALL exchange_halo_r(raing,                               &
                           memoryStart(1), memoryEnd(1),          &
                           memoryStart(2), memoryEnd(2),          &
                           memoryStart(3), memoryEnd(3),          &
                           patchStart(1),  patchEnd(1),           &
                           patchStart(2),  patchEnd(2),           &
                           patchStart(3),  patchEnd(3) )
      END IF

      IF (dbglvl > 0) WRITE(6,FMT=fmtstr) ' ...  DONE.'
    ELSE
      IF (dbglvl > 0) WRITE(6,FMT=fmtstr) ' ...  ERROR.'
    END IF
!-------------------- Reading RAINC ------------------------------------

    IF (dbglvl > 0) &
    WRITE(6,FMT='(3x,2a)',ADVANCE='NO') 'Reading 3D variable ', TRIM(varname2)

    IF ( multi_file .AND. mOrder(1:2) == 'XY' ) THEN

      DO jloc = 1, nysubf
        DO iloc = 1, nxsubf
          msloc(1) = ipstart(iloc,jloc); meloc(1) = ipend(iloc,jloc)
          msloc(2) = jpstart(iloc,jloc); meloc(2) = jpend(iloc,jloc)
          msloc(3) = patchStart(3);      meloc(3) = patchEnd(3)

          nx  = ipend(iloc,jloc) - ipstart(iloc,jloc) + 1
          ny  = jpend(iloc,jloc) - jpstart(iloc,jloc) + 1
          nz  = patchEnd(3) - patchStart(3) + 1

          mInsize = nx*ny*nz

          IF (mInsize > dmrsize) THEN
            DEALLOCATE ( temdmr )
          	dmrsize = mInsize
            ALLOCATE (temdmr(dmrsize), STAT = istatus)
          END IF

          IF (io_form_input == NETCDF) THEN
            CALL ext_ncd_read_field(fHndlr(iloc,jloc), DateStr,VarName2,&
                     temdmr,WRF_REAL, 0,0,1, mOrder, Stagger, DimNames, &
                     mSloc,mEloc, mSloc,mEloc, mSloc,mEloc,iStatus)

          ELSE IF (io_form_input == BINARY) THEN

            CALL ext_int_read_field(fHndlr(iloc,jloc), DateStr,VarName2,&
                     temdmr,WRF_REAL,0, 0, 1, mOrder, Stagger, DimNames,&
                     mSloc,mEloc, mSloc,mEloc, mSloc,mEloc, iStatus)
          END IF

          ! To merge
          DO k = patchStart(3), patchEnd(3)
            kin = (k - patchStart(3))*nx*ny
        	  DO j = jpstart(iloc,jloc), jpend(iloc,jloc)
        	    jin = (j - jpstart(iloc,jloc))*nx
        	    DO i = ipstart(iloc,jloc), ipend(iloc,jloc)
        	      kjiin = kin + jin + i - ipstart(iloc,jloc) + 1    ! relative position in memory
        	      rainc(i,j,k) = temdmr(kjiin)
        	    END DO
        	  END DO
          END DO

        END DO
      END DO

    ELSE

      nx = domainEnd(1) - domainStart(1) + 1
      ny = domainEnd(2) - domainStart(2) + 1
      nz = domainEnd(3) - domainStart(3) + 1
      mInsize = nx*ny*nz

      IF (mInsize > dmrsize) THEN
        DEALLOCATE ( temdmr )
      	dmrsize = mInsize
        ALLOCATE ( temdmr(dmrsize), STAT = istatus)
      END IF

      nfid = fHndlr(1,1)
      IF (IAMIONODE) THEN
        IF (io_form_input == NETCDF) THEN
          CALL ext_ncd_read_field(nfid, DateStr,VarName2,temdmr,WRF_REAL,&
                                  0,0,1, mOrder, Stagger, DimNames,     &
                                  DomainStart,DomainEnd,                &
                                  DomainStart,DomainEnd,                &
                                  DomainStart,DomainEnd,                &
                                  iStatus)

        ELSE IF (io_form_input == BINARY) THEN

          CALL ext_int_read_field(nfid, DateStr,VarName2,temdmr,WRF_REAL,&
                                  0, 0, 1, mOrder, Stagger, DimNames,   &
                                  DomainStart,DomainEnd,                &
                                  DomainStart,DomainEnd,                &
                                  DomainStart,DomainEnd,                &
                                  iStatus)
        END IF
      END IF  ! IAMIONODE
      CALL parallel_bcast_int(iStatus)

      IF (iStatus == 0) THEN    ! successful read
        IF (ntprocs > 1 .AND. mOrder(1:2) == 'XY') THEN
            CALL scatter_whole_field_r(rainc,                           &
                                     memoryStart(1), memoryEnd(1),      &
                                     memoryStart(2), memoryEnd(2),      &
                                     memoryStart(3), memoryEnd(3),      &
                                     patchStart(1),  patchEnd(1),       &
                                     patchStart(2),  patchEnd(2),       &
                                     patchStart(3),  patchEnd(3),       &
                                     temdmr,                            &
                                     domainstart(1), domainend(1),      &
                                     domainstart(2), domainend(2),      &
                                     domainstart(3), domainend(3) )
        ELSE
          IF (ntprocs > 1) CALL parallel_bcast_real_array(temdmr,mInsize)

          DO k = patchStart(3), patchEnd(3)
            kin = ( k - patchStart(3) )*ny*nx  ! patch and domain should be the same
        	  DO j = patchStart(2), patchEnd(2)
        	    jin = ( j - patchStart(2) )*nx
        	    DO i = patchStart(1), patchEnd(1)
        	      kjiin = kin + jin + i - patchStart(1) + 1
        	  	  rainc(i,j,k) = temdmr(kjiin)
        	    END DO
        	  END DO
    	    END DO

        END IF
      END IF
    END IF

    WRITE(fmtstr,'(a,I2.2,a)') '(',10-LEN_TRIM(varname2),'x,a)'

    IF (istatus == 0) THEN
      IF (mOrder(1:2) == 'XY') THEN
        CALL exchange_halo_r(rainc,                               &
                           memoryStart(1), memoryEnd(1),          &
                           memoryStart(2), memoryEnd(2),          &
                           memoryStart(3), memoryEnd(3),          &
                           patchStart(1),  patchEnd(1),           &
                           patchStart(2),  patchEnd(2),           &
                           patchStart(3),  patchEnd(3) )
      END IF

      IF (dbglvl > 0) WRITE(6,FMT=fmtstr) ' ...  DONE.'
    ELSE
      IF (dbglvl > 0) WRITE(6,FMT=fmtstr) ' ...  ERROR.'
    END IF

!--------------------------- Close the file ----------------------------

    DO j = 1,nysubf
      DO i = 1,nxsubf
        IF (io_form_input == BINARY) THEN
          IF (IAMIONODE) CALL ext_int_ioclose(fHndlr(i,j),iStatus)
        ELSE IF(io_form_input == NETCDF) THEN
          IF (IAMIONODE) CALL ext_ncd_ioclose(fHndlr(i,j),iStatus)
        END IF
      END DO
    END DO

    IF (istatus /= 0) THEN
      WRITE(0,'(1x,a)') 'ERROR: closing file handler.'
      CALL arpsstop('Error in nmm_get_rain.',1)
    END IF

    RETURN
  END SUBROUTINE nmm_get_rain

  SUBROUTINE nmm_get_dom_ti_integer( fid, element, val, ireturn)
    IMPLICIT NONE

    INTEGER,      INTENT(IN)  :: fid
    CHARACTER(*), INTENT(IN)  :: element
    INTEGER,      INTENT(OUT) :: val

    INTEGER,      INTENT(OUT) :: ireturn

!-----------------------------------------------------------------------

    INTEGER    :: outcount

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Begin of executable code
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    IF (io_form_input == NETCDF) THEN
       CALL ext_ncd_get_dom_ti_integer(fid,element,val,1,outcount,ireturn)
    ELSE IF (io_form_input == 1) THEN
       CALL ext_int_get_dom_ti_integer(fid,element,val,1,outcount,ireturn)
    END IF

    RETURN
  END SUBROUTINE nmm_get_dom_ti_integer

END MODULE module_nmm_input
