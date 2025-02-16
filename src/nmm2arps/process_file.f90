SUBROUTINE process_file(OUNT,nf,IAMROOT,istatus)

!#######################################################################

  USE module_nmm2arps_namelist
  USE module_nmm2arps_nmmgrid
  USE module_output_arpsgrid
  USE module_nmm_input

  USE module_interpolation

  IMPLICIT NONE

  INTEGER, INTENT(IN)  :: OUNT
  INTEGER, INTENT(IN)  :: nf
  LOGICAL, INTENT(IN)  :: IAMROOT

  INTEGER, INTENT(OUT) :: istatus
!-----------------------------------------------------------------------

  INTEGER :: it, itimeb, ifile
  CHARACTER(LEN=19) :: datetimestr
  REAL    :: curtime
  LOGICAL :: irain, i2dout

  INTEGER :: lvldbg_global

  REAL, PARAMETER :: latdiag = 34.5606
  REAL, PARAMETER :: londiag = -103.0820

  REAL, ALLOCATABLE :: tem1(:,:,:)

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  istatus = 0

  i2dout = i2dfmt > 0 .OR. igempak == 1
  irain  = i2dout .AND. iaccu == 1

  IF ( .NOT. input_initialized) THEN
    CALL nmm_input_init(IAMROOT,io_form,magnitude_processor,            &
                        multifile,ncompressx,ncompressy,istatus)
  END IF

  CALL nmm_input_open(wrffiles(nf),lvldbg,istatus)

  itimeb = 0                 ! to extract the right time string
  DO ifile = 1,nf-1
    itimeb = itimeb + frames_per_outfile(ifile)
  END DO

  lvldbg_global = lvldbg
  CALL mpupdatei(lvldbg_global,1)          ! For some diagnostic outputs

  DO it = 1, frames_per_outfile(nf)

    datetimestr = datestrs(itimeb+it)

    !
    ! Get WRF data
    !
    CALL nmmgrid_getdata( datetimestr, irain, initsec, itimeb+it,       &
                          abstimei, dir_extp, grid_id, lvldbg,          &
                          istatus )

    CALL nmmgrid_adjustdata(IAMROOT, lvldbg, istatus)
    !
    ! Set up interpolation working arrays
    !

    IF ( arpsgrid%first_time ) THEN

      IF (lvldbg > 2) WRITE(6,'(1x,a)') '=== Preparing horizontal interpolation arrays ...'

      CALL set_interpolation_arrays( nmmgrid%ims,nmmgrid%ime,           &
             nmmgrid%jms,nmmgrid%jme, nmmgrid%kps, nmmgrid%kpse,        &
             arpsgrid%nx, arpsgrid%ny,                                  &
             arpsgrid%xlat, arpsgrid%xlon,                              &
             arpsgrid%ylat, arpsgrid%ylon,                              &
             arpsgrid%slat, arpsgrid%slon,                              &
             istatus )

      IF (lvldbg > 2) WRITE(6,'(1x,a)') '=== Preparing vertical interpolation arrays ...'

      CALL set_interpolation_heights( arpsgrid%nx,arpsgrid%ny,intrp_method, &
                    nmmgrid%ims,nmmgrid%ime,nmmgrid%jms, nmmgrid%jme,   &
                    nmmgrid%kps,nmmgrid%kpe,nmmgrid%kpse,               &
                    nmmgrid%kzs,nmmgrid%kzse,                           &
                    nmmgrid%zp%rdata_arr,nmmgrid%hgt%rdata_arr,         &
                    nmmgrid%sldpth%rdata_arr,lvldbg,istatus )

      IF (lvldbg > 2) WRITE(6,'(1x,a)') '=== Retouching ARPS grid ...'

      CALL arpsgrid_terrain_heights (wrfexttrnopt, extntmrg, intrp_method,  &
          nmmgrid%ims,nmmgrid%ime,nmmgrid%jms,nmmgrid%jme,              &
          nmmgrid%fis%rdata_arr, lvldbg, istatus )

    END IF

!
!-----------------------------------------------------------------------
!
!  Calculate base-state sounding (vertical profile)
!
!-----------------------------------------------------------------------
!
    IF (lvldbg > 2) WRITE(6,'(1x,a)') '=== Preparing base sounding ...'

    CALL set_interpolation_basesnd(arpsgrid%nx,arpsgrid%ny,             &
           nmmgrid%ims,nmmgrid%ime,nmmgrid%jms,nmmgrid%jme,             &
           nmmgrid%kps,nmmgrid%kpe,                                     &
           arpsgrid%xs,arpsgrid%ys,                                     &
           nmmgrid%glat%rdata_arr,nmmgrid%glon%rdata_arr,               &
           nmmgrid%hgt%rdata_arr,nmmgrid%algp%rdata_arr,                &
           nmmgrid%tt%rdata_arr, nmmgrid%qv%rdata_arr,                  &
           nmmgrid%uu%rdata_arr, nmmgrid%vv%rdata_arr, lvldbg, istatus)

    IF (lvldbg > 2) WRITE(6,'(1x,a)') '=== Doing interpolation ...'

    !
    ! Do the interpolation
    !
    CALL interpolate_from_nmm_to_arps( IAMROOT, intrp_method, intropt,  &
                    nsmooth, wndadj, i2dout, lvldbg, istatus)

    IF (lvldbg > 2) WRITE(6,'(1x,a)') '=== Post-procesing of ARPS arrays ...'

    !
    ! Make necessary adjustment base on namelist parameter
    !
    CALL arpsgrid_adjustdata(IAMROOT, .TRUE.,nmmgrid%num_qq,            &
                             csopt, csfactr, csound,                    &
                             hydradj, dz, wndadj, obropt, obrzero,      &
                             ext_lbc, ext_vbc, tbc, bbc,lvldbg, istatus)
    !
    ! Diagnostic output if desired
    !
    IF (lvldbg > 0) THEN
      CALL arpsgrid_diagout( latdiag, londiag, lvldbg, istatus )
    END IF

    IF (lvldbg > 2) WRITE(6,'(1x,a)') 'Now, it is ready for output ...'

    !
    ! Do ARPS output
    !
    curtime = nmmgrid%xtime * 60   ! change from minutes to seconds

    IF ( terndmp > 0 .AND. arpsgrid%first_time ) THEN
      CALL arpsgrid_ternout( lvldbg, istatus )
    END IF

    IF ( soildmp  > 0 ) THEN
      CALL arpsgrid_soilout( curtime, lvldbg, istatus )
    END IF

    !IF( nmmgrid%mp_physics == 5 ) THEN  ! Backup F_RIMEF for 2D write
    !  ALLOCATE(tem1(nx,ny,nz),STAT = istatus)
    !  CALL arpsgrid_saveAndResetQH(tem1,nx,ny,nz)
    !  !tem1 = arpsgrid%qscalar(:,:,:,P_QH)
    !  !arpsgrid%qscalar(:,:,:,P_QH)   = 0.0
    !END IF

    CALL arpsgrid_histout( IAMROOT, curtime, lvldbg_global, istatus )

    !IF( nmmgrid%mp_physics == 5 ) THEN  ! Restore F_RIMEF to qh for 2D write
    !  CALL arpsgrid_restoreQH(tem1,nx,ny,nz)
    !  !arpsgrid%qscalar(:,:,:,P_QH)   = tem1
    !  DEALLOCATE(tem1)
    !END IF

    IF ( i2dout )  THEN     ! must be after histout

      CALL arpsgrid_2dout( IAMROOT, curtime, nmmgrid%mp_physics+100,    &
             outheader, gemoutheader,icape,iaccu,iascii,i2dfmt,igempak, &
             ilite,iltgci,icrtm,isatid,chbgn,chend,icitm,user_emis,     &
             ibeg_offset,iend_offset,jbeg_offset,jend_offset,           &
             lvldbg, istatus )
    END IF

    IF (arpsgrid%first_time) arpsgrid%first_time = .FALSE.

  END DO

!  IF (lvldbg > 2) WRITE(6,'(1x,a,I4)') '=== Done with file No ', nf

  CALL nmm_input_close(istatus)

  RETURN
END SUBROUTINE process_file

SUBROUTINE interpolate_from_nmm_to_arps( IAMROOT,intrpmethod,intropt, &
                             nsmooth, wndadj, i2dout, dbglvl, istatus)

  USE module_nmm2arps_nmmgrid
  USE module_output_arpsgrid
  USE module_interpolation

!#######################################################################

  IMPLICIT NONE

  LOGICAL, INTENT(IN)  :: IAMROOT
  INTEGER, INTENT(IN)  :: intrpmethod, intropt
  INTEGER, INTENT(IN)  :: nsmooth, wndadj
  LOGICAL, INTENT(IN)  :: i2dout
  INTEGER, INTENT(IN)  :: dbglvl
  INTEGER, INTENT(OUT) :: istatus

!-----------------------------------------------------------------------

  INCLUDE 'mp.inc'
  INCLUDE 'globcst.inc'
  INCLUDE 'bndry.inc'

!-----------------------------------------------------------------------
!
!  Non-dimensional smoothing coefficient
!
!-----------------------------------------------------------------------
!
  REAL, PARAMETER :: smfct1 = 0.5

  INTEGER :: ims, ime, jms, jme, kms, kme, ksme
  INTEGER :: nx, ny, nz

  INTEGER :: i,j,k, nq, ksmth

  INTEGER :: iprtopt

  REAL    :: amax, amin

  REAL, ALLOCATABLE :: tsoil_tmp(:,:,:), qsoil_tmp(:,:,:)

  REAL, ALLOCATABLE :: tem1(:,:,:)
  REAL, ALLOCATABLE :: tem1_ext(:,:,:)

  REAL :: intrp_pnt_from_2d

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  istatus = 0

  ims = nmmgrid%ims
  ime = nmmgrid%ime
  jms = nmmgrid%jms
  jme = nmmgrid%jme
  kms = nmmgrid%kps
  kme = nmmgrid%kpe
  ksme = nmmgrid%kpse

  nx = arpsgrid%nx
  ny = arpsgrid%ny
  nz = arpsgrid%nz

  ALLOCATE(tem1(nx,ny,nz),     STAT = istatus)
  ALLOCATE(tem1_ext(ims:ime,jms:jme,kms:ksme), STAT = istatus)
!
!-----------------------------------------------------------------------
!
!  Process Pressure data
!
!-----------------------------------------------------------------------
!
  IF (dbglvl > 9) WRITE(6,'(3x,a)') '--- Processing pressure data ...'

  iprtopt=1  ! Produce mean and perturbation field
  CALL mkarpsvlz(ims,ime,jms,jme,kms,kme,nx,ny,nz, 'P', M,              &
                 intrpmethod,iprtopt,intropt,plsnd,                     &
                 nmmgrid%hgt%rdata_arr,arpsgrid%zps,                    &
                 nmmgrid%pp%rdata_arr,                                  &
                 arpsgrid%pbar,arpsgrid%pprt, tem1_ext,                 &
                 dbglvl,istatus)

!  IF (nsmooth > 0 .AND. mp_opt > 0) THEN      ! Ensure passed in array is mpi valid
!    CALL acct_interrupt(mp_acct)
!    CALL mpsendrecv2dew(arpsgrid%pprt,nx,ny,nz,ebc,wbc,0,tem1)
!    CALL mpsendrecv2dns(arpsgrid%pprt,nx,ny,nz,nbc,sbc,0,tem1)
!  END IF

  DO ksmth=1,nsmooth
    DO k=1,nz
      CALL smooth9p(arpsgrid%pprt(:,:,k), nx,ny, 1,nx,1,ny,0, tem1)
    END DO

    IF (nsmooth > ksmth .AND. mp_opt > 0) THEN
                          ! Ensure output is good for followed smoothing
      CALL acct_interrupt(mp_acct)
      CALL mpsendrecv2dew(arpsgrid%pprt,nx,ny,nz,ebc,wbc,0,tem1)
      CALL mpsendrecv2dns(arpsgrid%pprt,nx,ny,nz,nbc,sbc,0,tem1)
    END IF
  END DO

!
!-----------------------------------------------------------------------
!
!  Process potential temperature data
!
!-----------------------------------------------------------------------
!
  IF (dbglvl > 9) WRITE(6,'(3x,a)') '--- Processing potential temperature data ...'

  iprtopt=1
  CALL mkarpsvar(ims,ime,jms,jme,kms,kme,nx,ny,nz,'PT',M,               &
                 intrpmethod,iprtopt,intropt,ptsnd,                     &
                 nmmgrid%hgt%rdata_arr,arpsgrid%zps,                    &
                 nmmgrid%pt%rdata_arr,                                  &
                 arpsgrid%ptbar,arpsgrid%ptprt, tem1_ext, dbglvl, istatus)

!  IF (mp_opt > 0) THEN      ! Ensure passed in array is mpi valid
!    CALL acct_interrupt(mp_acct)
!    CALL mpsendrecv2dew(arpsgrid%ptprt,nx,ny,nz,ebc,wbc,0,tem1)
!    CALL mpsendrecv2dns(arpsgrid%ptprt,nx,ny,nz,nbc,sbc,0,tem1)
!  END IF

  DO ksmth=1,nsmooth
    CALL smooth3d(nx,ny,nz, 1,nx,1,ny,1,nz,0, smfct1,arpsgrid%zps,      &
                  arpsgrid%ptprt,tem1,arpsgrid%ptprt)
  END DO

!
!-----------------------------------------------------------------------
!
!  Process qv data
!
!-----------------------------------------------------------------------
!
  IF (dbglvl > 9) WRITE(6,'(3x,a)') '--- Processing RHstar(QV) data ...'
  !
  ! QV stores RHstar data now
  !

  iprtopt=0        ! produce total field instead of only perturbation
  CALL mkarpsvar(ims,ime,jms,jme,kms,kme,nx,ny,nz,'RHS',M,              &
                 intrpmethod,iprtopt,intropt,rhssnd,                    &
                 nmmgrid%hgt%rdata_arr,arpsgrid%zps,                    &
                 nmmgrid%qv%rdata_arr, arpsgrid%qvbar,arpsgrid%qv,      &
                 tem1_ext,dbglvl, istatus )

!  IF (nsmooth > 0 .AND. mp_opt > 0) THEN      ! Ensure passed in array is mpi valid
!    CALL acct_interrupt(mp_acct)
!    CALL mpsendrecv2dew(arpsgrid%qv,nx,ny,nz,ebc,wbc,0,tem1)
!    CALL mpsendrecv2dns(arpsgrid%qv,nx,ny,nz,nbc,sbc,0,tem1)
!  END IF

  DO k = 1,arpsgrid%nz
    DO j = 1, arpsgrid%ny
      DO i = 1, arpsgrid%nx
        arpsgrid%qv(i,j,k) = MAX(0.0,arpsgrid%qv(i,j,k))
      END DO
    END DO
  END DO

  DO ksmth=1,nsmooth
    CALL smooth3d(nx,ny,nz,1,nx,1,ny,1,nz,0,smfct1,arpsgrid%zps,        &
                  arpsgrid%qv,tem1,arpsgrid%qv)
  END DO

  !
  !  Convert rhstar back to qv for writing, further calculations
  !
  CALL nmmgrid_qv2rhstar( IAMROOT, -1, dbglvl, istatus )

  CALL arpsgrid_qv2rhstar( IAMROOT, -1, dbglvl, istatus )

!-----------------------------------------------------------------------
!
!  Process qc data.
!
!  If no qc info is provided...set to zero.
!
!-----------------------------------------------------------------------
!
  DO nq = 1, nscalar
    IF (dbglvl > 9) WRITE(6,'(3x,a)') '--- Processing '//TRIM(qnames(nq))//' data ...'

    iprtopt=0
    CALL mkarpsvar(ims,ime,jms,jme,kms,kme,nx,ny,nz,TRIM(qnames(nq)),M, &
                   intrpmethod,iprtopt,intropt,dumsnd,                  &
                   nmmgrid%hgt%rdata_arr,arpsgrid%zps,                  &
                   nmmgrid%qscalar(nq)%rdata_arr,tem1,arpsgrid%qscalar(:,:,:,nq), &
                   tem1_ext,dbglvl, istatus )

!    IF (nsmooth > 0 .AND. mp_opt > 0) THEN      ! Ensure passed in array is mpi valid
!      CALL acct_interrupt(mp_acct)
!      CALL mpsendrecv2dew(arpsgrid%qc,nx,ny,nz,ebc,wbc,0,tem1)
!      CALL mpsendrecv2dns(arpsgrid%qc,nx,ny,nz,nbc,sbc,0,tem1)
!    END IF

    DO k = 1,arpsgrid%nz
      DO j = 1, arpsgrid%ny
        DO i = 1, arpsgrid%nx
          arpsgrid%qscalar(i,j,k,nq) = MAX(0.0,arpsgrid%qscalar(i,j,k,nq))
        END DO
      END DO
    END DO

    DO ksmth=1,nsmooth
      CALL smooth3d(nx,ny,nz,1,nx,1,ny,1,nz,0,smfct1,arpsgrid%zps,      &
                    arpsgrid%qscalar(:,:,:,nq),tem1,arpsgrid%qscalar(:,:,:,nq))
    END DO

  END DO

!-----------------------------------------------------------------------
!
!  Process tke data.
!
!-----------------------------------------------------------------------
!
!  IF (dbglvl > 9) WRITE(6,'(3x,a)') '--- Processing TKE data ...'
!
!  iprtopt=0
!  CALL mkarpsvar(ims,ime,jms,jme,kms,kme,nx,ny,nz,'TKE',M,              &
!                 intrpmethod,iprtopt,intropt,dumsnd,                    &
!                 nmmgrid%hgt%rdata_arr,arpsgrid%zps,                    &
!                 nmmgrid%tke%rdata_arr,tem1,arpsgrid%tke,               &
!                 tem1_ext,dbglvl, istatus )
!
!!  IF (nsmooth > 0 .AND. mp_opt > 0) THEN      ! Ensure passed in array is mpi valid
!!    CALL acct_interrupt(mp_acct)
!!    CALL mpsendrecv2dew(arpsgrid%qc,nx,ny,nz,ebc,wbc,0,tem1)
!!    CALL mpsendrecv2dns(arpsgrid%qc,nx,ny,nz,nbc,sbc,0,tem1)
!!  END IF
!
!  DO ksmth=1,nsmooth
!    CALL smooth3d(nx,ny,nz,1,nx,1,ny,1,nz,0,smfct1,arpsgrid%zps,        &
!                  arpsgrid%tke,tem1,arpsgrid%tke)
!  END DO
!
!-----------------------------------------------------------------------
!
!  Process density which has been stuffed into tem5_ext array.
!  We really only need rhobar, so set iprtopt=1 and pass
!  a temporary array in place of rhoprt.
!
!-----------------------------------------------------------------------
!
  IF (dbglvl > 9) WRITE(6,'(3x,a)') '--- Processing air density ...'

  iprtopt=1
  CALL mkarpsvar(ims,ime,jms,jme,kms,kme,nx,ny,nz,'RHO',M,              &
                intrpmethod,iprtopt,intropt,rhosnd,                     &
                nmmgrid%hgt%rdata_arr,arpsgrid%zps,                     &
                nmmgrid%rho%rdata_arr, arpsgrid%rhobar,tem1,            &
                tem1_ext, dbglvl, istatus )

  DO ksmth=1,nsmooth
    CALL smooth3d(nx,ny,nz, 1,nx,1,ny,1,nz,0,smfct1,arpsgrid%zps,       &
                  arpsgrid%rhobar,tem1,arpsgrid%rhobar)
  END DO

!
!-----------------------------------------------------------------------
!
!    Process U wind data
!
!-----------------------------------------------------------------------
!
  IF (dbglvl > 9) WRITE(6,'(3x,a)') '--- Processing U wind data ...'

  iprtopt=0
  CALL mkarpsvar(ims,ime,jms,jme,kms,kme,nx,ny,nz,'U',U,                &
                 intrpmethod,iprtopt,intropt,usnd,                      &
                 nmmgrid%hgt%rdata_arr,arpsgrid%zpu,                    &
                 nmmgrid%uu%rdata_arr,arpsgrid%ubar,arpsgrid%u,         &
                 tem1_ext, dbglvl, istatus )

!  IF (mp_opt > 0) THEN      ! Ensure passed in array is mpi valid
!    CALL acct_interrupt(mp_acct)
!    CALL mpsendrecv2dew(arpsgrid%u,nx,ny,nz,ebc,wbc,1,tem1)
!    CALL mpsendrecv2dns(arpsgrid%u,nx,ny,nz,nbc,sbc,1,tem1)
!  END IF

  DO ksmth=1,nsmooth
    CALL smooth3d(nx,ny,nz,1,nx,1,ny,1,nz,1,smfct1,arpsgrid%zpu,        &
                  arpsgrid%u,tem1,arpsgrid%u)
  END DO

!  IF (dbglvl > 2) THEN
!    CALL a3dmax0(arpsgrid%u,1,nx,1,nx,1,ny,1,ny-1,1,nz,1,nz-1,amax,amin)
!    IF( IAMROOT ) WRITE(6,'(1x,2(a,e13.6))')                        &
!            'umin final= ', amin,',  umax final=',amax
!  END IF
!
!-----------------------------------------------------------------------
!
!  Process v component
!
!-----------------------------------------------------------------------
!
  IF (dbglvl > 9) WRITE(6,'(3x,a)') '--- Processing V wind data ...'

  iprtopt=0
  CALL mkarpsvar(ims,ime,jms,jme,kms,kme,nx,ny,nz,'V',V,                &
                 intrpmethod,iprtopt,intropt,vsnd,                      &
                 nmmgrid%hgt%rdata_arr,arpsgrid%zpv,                    &
                 nmmgrid%vv%rdata_arr, arpsgrid%vbar,arpsgrid%v,        &
                 tem1_ext,dbglvl, istatus )

!  IF (mp_opt > 0) THEN      ! Ensure passed in array is mpi valid
!    CALL acct_interrupt(mp_acct)
!    CALL mpsendrecv2dew(arpsgrid%v,nx,ny,nz,ebc,wbc,2,tem1)
!    CALL mpsendrecv2dns(arpsgrid%v,nx,ny,nz,nbc,sbc,2,tem1)
!  END IF

  DO ksmth=1,nsmooth
    CALL smooth3d(nx,ny,nz,1,nx,1,ny,1,nz,2,smfct1,arpsgrid%zpv,        &
                  arpsgrid%v,tem1,arpsgrid%v)
  END DO

!  IF (dbglvl > 2) THEN
!    CALL a3dmax0(arpsgrid%v,1,nx,1,nx-1,1,ny,1,ny,1,nz,1,nz-1,amax,amin)
!    IF(myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                          &
!            'vmin final= ', amin,',  vmax final=',amax
!  END IF

!-----------------------------------------------------------------------
!
!  Find vertical velocity
!
!-----------------------------------------------------------------------
  IF(wndadj == 0) THEN  ! Otherwise, derive W from U, V

    IF (dbglvl > 9) WRITE(6,'(3x,a)') '--- Processing W wind data ...'

    iprtopt=0
    CALL mkarpsvar(ims,ime,jms,jme,kms,ksme,nx,ny,nz,'W',W,             &
                   intrpmethod,iprtopt,intropt,dumsnd,                  &
                   nmmgrid%hgt%rdata_arr,arpsgrid%zp,                   &
                   nmmgrid%ww%rdata_arr, arpsgrid%wbar,arpsgrid%w,      &
                   tem1_ext,dbglvl,istatus )

!    IF (nsmooth > 0 .AND. mp_opt > 0) THEN      ! Ensure passed in array is mpi valid
!      CALL acct_interrupt(mp_acct)
!      CALL mpsendrecv2dew(w,nx,ny,nz,ebc,wbc,3,tem1)
!      CALL mpsendrecv2dns(w,nx,ny,nz,nbc,sbc,3,tem1)
!    END IF

    DO ksmth=1,nsmooth
      CALL smooth3d(nx,ny,nz, 1,nx,1,ny,1,nz,3,smfct1,arpsgrid%zp,      &
                    arpsgrid%w,tem1,arpsgrid%w)
    END DO
  END IF

!
!-----------------------------------------------------------------------
!
!  Process 4D surface fields if required (sfcout = 1)
!
!-----------------------------------------------------------------------
!
!  wetcout  = 1
!  tsoilout = 1
!  qsoilout = 1

  IF (dbglvl > 9) WRITE(6,'(3x,a)') '--- Processing soil fields ...'

  ALLOCATE(tsoil_tmp(nx,ny,nmmgrid%soil_layer_stag), STAT = istatus)
  ALLOCATE(qsoil_tmp(nx,ny,nmmgrid%soil_layer_stag), STAT = istatus)

  CALL egridfill(nmmgrid%tsoil%rdata_arr,ims,ime,jms,jme,               &
                 nmmgrid%kzs, nmmgrid%kzse,HH,istatus)

  ! Horizontally interpolations
  DO k = nmmgrid%kzs, nmmgrid%kzse
    DO j = 1,ny
      DO i = 1,nx
        ! Soil temperature
        tsoil_tmp(i,j,k) = intrp_pnt_from_2d(intrpmethod,               &
                      isr(i,j),jsr(i,j),iEms,iEme,jms,jme,              &
                      nmmgrid%kzs, nmmgrid%kzse,                        &
                      k,workarr,dbglvl,istatus)
      END DO
    END DO
  END DO

  CALL egridfill(nmmgrid%qsoil%rdata_arr,ims,ime,jms,jme,               &
                 nmmgrid%kzs, nmmgrid%kzse,HH,istatus)

  DO k = nmmgrid%kzs, nmmgrid%kzse
    DO j = 1,ny
      DO i = 1,nx
        ! Soil moisture
        qsoil_tmp(i,j,k) = intrp_pnt_from_2d(intrpmethod,               &
                      isr(i,j),jsr(i,j),iEms,iEme,jms,jme,              &
                      nmmgrid%kzs, nmmgrid%kzse,                        &
                      k,workarr,dbglvl,istatus)
      END DO
    END DO
  END DO

  CALL egridfill(nmmgrid%canpwet%rdata_arr,ims,ime,jms,jme,             &
                 1, 1, HH,istatus)

  DO j = 1,ny
    DO i = 1,nx
      ! wet canopy
      arpsgrid%wetcanp(i,j,0) = intrp_pnt_from_2d(intrpmethod,          &
                      isr(i,j),jsr(i,j),iEms,iEme,jms,jme,1,1,          &
                      1,workarr,dbglvl,istatus)
    END DO
  END DO


  ! First convert zpsoil to soil depth.
  DO k = 1,arpsgrid%nzsoil
    DO j = 1,ny
      DO i = 1,nx
        arpsgrid%zpsoil(i,j,k) = arpsgrid%hterain(i,j) - arpsgrid%zpsoil(i,j,k)
      END DO
    END DO
  END DO

  CALL vintrpsoil_wrf(nx,ny,nmmgrid%soil_layer_stag,arpsgrid%nzsoil,    &
                      soilmodel_option,zpsoil_tmp,arpsgrid%zpsoil,      &
                      tsoil_tmp,qsoil_tmp,                              &
                      arpsgrid%tsoil(:,:,:,0),arpsgrid%qsoil(:,:,:,0))

  CALL fix_soil_nstyp(nx,ny,arpsgrid%nzsoil,1,arpsgrid%nstyps,          &
                      arpsgrid%tsoil,arpsgrid%qsoil,arpsgrid%wetcanp)

  !DO is = 1,nstyp
  !  CALL edgfill(tsoil(:,:,:,is),1,nx,1,nx-1,1,ny,1,ny-1,           &
  !               1,nzsoil,1,nzsoil)
  !  CALL edgfill(qsoil(:,:,:,is),1,nx,1,nx-1,1,ny,1,ny-1,           &
  !               1,nzsoil,1,nzsoil)
  !  CALL edgfill(wetcanp(:,:,is),1,nx,1,nx-1,1,ny,1,ny-1,1,1,1,1)
  !END DO

  ! Convert zpsoil back from soil depth to MSL
  DO k = 1,arpsgrid%nzsoil
    DO j = 1,ny
      DO i = 1,nx
        arpsgrid%zpsoil(i,j,k) = arpsgrid%hterain(i,j) - arpsgrid%zpsoil(i,j,k)
      END DO
    END DO
  END DO

  DEALLOCATE( tsoil_tmp, qsoil_tmp )

!-----------------------------------------------------------------------
!
! Processing Snow depth, soil type and vegtation type etc.
!
!-----------------------------------------------------------------------

  ! Snow depth using nearest neighbor
  ! snowdout = 1
  IF (dbglvl > 9) WRITE(6,'(3x,a)') '--- Processing 2D surface fields ...'

  CALL mkarps2d(nx,ny,NEARNEIGHBOR,M,ims,ime,jms,jme,                   &
                nmmgrid%snowh%rdata_arr,arpsgrid%snowdpth,dbglvl,istatus )

  CALL mkarps2di(nx,ny,NEARNEIGHBOR1,M,ims,ime,jms,jme,                 &
                nmmgrid%isltyp,arpsgrid%soiltyp,dbglvl,istatus )

  CALL mkarps2di(nx,ny,NEARNEIGHBOR2,M,ims,ime,jms,jme,                 &
                nmmgrid%ivgtyp,arpsgrid%vegtyp,dbglvl,istatus )

  CALL mkarps2d(nx,ny,NEARNEIGHBOR,M,ims,ime,jms,jme,                   &
                nmmgrid%vegfrc%rdata_arr,arpsgrid%veg,dbglvl,istatus )

  CALL mkarps2d(nx,ny,NEARNEIGHBOR,M,ims,ime,jms,jme,                   &
                nmmgrid%z0%rdata_arr,arpsgrid%roufns,dbglvl,istatus )

  arpsgrid%stypfrct(:,:,:) = 0.0
  arpsgrid%stypfrct(:,:,1) = 1.0

!-----------------------------------------------------------------------
!
!  Processing PBL grid scale precipitation
!
!-----------------------------------------------------------------------

  IF(rainout == 1) THEN
    iprtopt = 0

    IF (dbglvl > 9) WRITE(6,'(3x,a)') '--- Processing precipitation data ...'

    CALL mkarps2d(nx,ny,intrpmethod,M,ims,ime,jms,jme,                  &
                  nmmgrid%acprec%rdata_arr,arpsgrid%raing,              &
                  dbglvl,istatus )

!-----------------------------------------------------------------------
!
!  Processing cumulus precipitation
!
!-----------------------------------------------------------------------

    CALL mkarps2d(nx,ny,intrpmethod,M,ims,ime,jms,jme,                  &
                  nmmgrid%cuprec%rdata_arr,arpsgrid%rainc,              &
                  dbglvl,istatus )

    DO j = 1,arpsgrid%ny           ! Convert WRF total precipitation to
    	 DO i = 1, arpsgrid%nx       ! PBL grid scale precipitation.
    	   arpsgrid%raing(i,j) = arpsgrid%raing(i,j) - arpsgrid%rainc(i,j)
    	 END DO
    END DO

  ELSE
    arpsgrid%raing(:,:) = 0.0
    arpsgrid%rainc(:,:) = 0.0
  END IF
!
!-----------------------------------------------------------------------
!
!  Processing t2m, th2m, qv2m, u10m, v10m
!
!-----------------------------------------------------------------------

  IF ( i2dout ) THEN
    IF (dbglvl > 9) WRITE(6,'(3x,a)') '--- Processing 2m & 10m data ...'

    CALL mkarps2d(nx,ny,intrpmethod,M,ims,ime,jms,jme,                  &
                  nmmgrid%th2%rdata_arr,arpsgrid%th2m,dbglvl,istatus )

    CALL mkarps2d(nx,ny,intrpmethod,M,ims,ime,jms,jme,                  &
                  nmmgrid%tshltr%rdata_arr,arpsgrid%t2m,dbglvl,istatus )

    CALL mkarps2d(nx,ny,intrpmethod,M,ims,ime,jms,jme,                  &
                  nmmgrid%qshltr%rdata_arr,arpsgrid%qv2m,dbglvl,istatus )

    CALL mkarps2d(nx,ny,intrpmethod,M,ims,ime,jms,jme,                  &
                  nmmgrid%u10%rdata_arr,arpsgrid%u10m,dbglvl,istatus )

    CALL mkarps2d(nx,ny,intrpmethod,M,ims,ime,jms,jme,                  &
                  nmmgrid%v10%rdata_arr,arpsgrid%v10m,dbglvl,istatus )
  END IF

!
!-----------------------------------------------------------------------
!
!  Processing SPRING2011 Specific Fields
!
!-----------------------------------------------------------------------

  IF ( nmmgrid%SPRING2011 ) THEN
    IF (dbglvl > 9) WRITE(6,'(3x,a)') '--- Processing hourly MAXIMUM fields for Spring 2010 ...'

    arpsgrid%refd_max=ARPS_MISSING
    arpsgrid%grpl_max=ARPS_MISSING

    CALL mkarps2d(nx,ny,intrpmethod,M,ims,ime,jms,jme,                  &
                  nmmgrid%wspd10max%rdata_arr,arpsgrid%wspd10max,dbglvl,istatus )

    CALL mkarps2d(nx,ny,intrpmethod,M,ims,ime,jms,jme,                  &
                  nmmgrid%w_up_max%rdata_arr,arpsgrid%w_up_max,dbglvl,istatus )

    CALL mkarps2d(nx,ny,intrpmethod,M,ims,ime,jms,jme,                  &
                  nmmgrid%w_dn_max%rdata_arr,arpsgrid%w_dn_max,dbglvl,istatus )

    IF( nmmgrid%mp_physics == 5 ) THEN
    CALL mkarps2d(nx,ny,intrpmethod,M,ims,ime,jms,jme,                  &
                  nmmgrid%refd_max%rdata_arr,arpsgrid%refd_max,dbglvl,istatus )
    END IF

    CALL mkarps2d(nx,ny,intrpmethod,M,ims,ime,jms,jme,                  &
                  nmmgrid%up_heli_max%rdata_arr,arpsgrid%up_heli_max,dbglvl,istatus )

    CALL mkarps2d(nx,ny,intrpmethod,M,ims,ime,jms,jme,                  &
                  nmmgrid%rswdn%rdata_arr,arpsgrid%raddn,dbglvl,istatus )
    ! rswdn is actually holding total sfc downward radiation flux

    !IF (dbglvl > 2) THEN
    !  CALL a3dmax0(arpsgrid%raddn,1,nx,1,nx-1,1,ny,1,ny,1,1,1,1,amax,amin)
    !  IF(myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                        &
    !          'raddn_min= ', amin,',  raddn_max=',amax
    !END IF
  END IF

!------------------ Before return --------------------------------------

  DEALLOCATE( tem1 )
  DEALLOCATE( tem1_ext )

  RETURN
END SUBROUTINE interpolate_from_nmm_to_arps
