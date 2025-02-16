!
SUBROUTINE find_nmm_halo_width(ips,ipe,jps,jpe,nx,ny,                   &
                               xlat,xlon,ylat,ylon,slat,slon,           &
                               myproc,dbglvl,fake_width,istatus)

  USE module_wrfgrid_constants
  USE wrf_llxy_module

  IMPLICIT NONE
  INTEGER, INTENT(IN)  :: ips, ipe, jps, jpe
  INTEGER, INTENT(IN)  :: nx, ny
  REAL,    INTENT(IN)  :: xlat(nx,ny), xlon(nx,ny)
  REAL,    INTENT(IN)  :: ylat(nx,ny), ylon(nx,ny)
  REAL,    INTENT(IN)  :: slat(nx,ny), slon(nx,ny)
  INTEGER, INTENT(IN)  :: myproc
  INTEGER, INTENT(IN)  :: dbglvl
  INTEGER, INTENT(OUT) :: fake_width
  INTEGER, INTENT(OUT) :: istatus

!-----------------------------------------------------------------------

  REAL :: iur, jur
  REAL :: ixmin, ixmax, jymin, jymax

  INTEGER :: ixwidth, jywidth
  INTEGER :: ipEs, ipEe

  INTEGER :: i, j

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  istatus = 0

  ipEs = 2*ips-1
  ipEe = 2*ipe

  ixmin = ipEs
  ixmax = ipEe

  jymin = jps
  jymax = jpe

!----------------------------------------------------------------------
!
! Find minimum and Maximum
!
!----------------------------------------------------------------------

  DO j = 1,ny,ny-1          ! U staggered points
    DO i = 1,nx,nx-1
      CALL wps_lltoxy(xlat(i,j),xlon(i,j),iur,jur,VV)

      IF (iur < ixmin) ixmin = iur
      IF (iur > ixmax) ixmax = iur
      IF (jur < jymin) jymin = jur
      IF (jur > jymax) jymax = jur

    END DO
  END DO

  IF (dbglvl > 10) WRITE(6,'(3x,a,I4,a,4(F5.0,a))')                      &
    'Processor -',myproc,' U point ixmin,ixmax,jymin,jymax = ',         &
    ixmin,',',ixmax,',',jymin,',',jymax,'.'

  DO j = 1,ny,ny-1           ! V staggered points
    DO i = 1,nx,nx-1
      CALL wps_lltoxy(ylat(i,j),ylon(i,j),iur,jur,VV)

      IF (iur < ixmin) ixmin = iur
      IF (iur > ixmax) ixmax = iur
      IF (jur < jymin) jymin = jur
      IF (jur > jymax) jymax = jur

    END DO
  END DO

  IF (dbglvl > 10) WRITE(6,'(3x,a,I4,a,4(F5.0,a))')                      &
    'Processor -',myproc,' V point ixmin,ixmax,jymin,jymax = ',         &
    ixmin,',',ixmax,',',jymin,',',jymax,'.'

!  DO j = 1,ny           ! Scalar points
!    DO i = 1,nx
!      CALL wps_lltoxy(slat(i,j),slon(i,j),iur,jur,HH)
!
!      IF (iur < ixmin) ixmin = iur
!      IF (iur > ixmax) ixmax = iur
!      IF (jur < jymin) jymin = jur
!      IF (jur > jymax) jymax = jur
!
!    END DO
!  END DO

  ixmin = ixmin - 1
  ixmax = ixmax + 1        ! We want to do at least 4 point interpolation
  jymin = jymin - 1
  jymax = jymax + 1

!-----------------------------------------------------------------------
!
! Find width
!
!-----------------------------------------------------------------------

  IF (dbglvl > 10) WRITE(6,'(3x,a,I4,a,4(I5,a))')                       &
    'Processor -',myproc,' Patch size ipEs,ipEe,jps,jpe    = ',         &
    ipEs,',',ipEe,',',jps,',',jpe,'.'

  ixwidth = MAX(ipEs - FLOOR(ixmin),CEILING(ixmax) - ipEe,0)
  jywidth = MAX(jps  - FLOOR(jymin),CEILING(jymax) - jpe, 0)

  ixwidth = CEILING(ixwidth/2.0)             ! it has been full E-grid

  IF (dbglvl > 0) THEN
    IF (ixwidth > (ipe-ips)) WRITE(6,'(3x,a,I4,a,2(I5,a),2(F5.0,a))')   &
      'Processor -',myproc,' Patch size ipEs,ipEe = ',                  &
      ipEs,',',ipEe,', ixmin,ixmax = ',ixmin,',',ixmax,'.'
    IF (jywidth > (jpe-jps)) WRITE(6,'(3x,a,I4,a,2(I5,a),2(F5.0,a))')   &
      'Processor -',myproc,' Patch size jps,jpe   = ',                  &
      jps,',',jpe,', jymin,jymax = ',jymin,',',jymax,'.'
  END IF

  fake_width = MAX(ixwidth, jywidth)         ! One width for both directions

  CALL mpmaxi(fake_width)                    ! Global maximum

  IF (fake_width > (ipe-ips) .OR. fake_width > (jpe-jps) ) THEN
    WRITE(6,'(1x,a,I4,a,/,8x,a,/,8x,a,/,8x,a,/)')                       &
      'ERROR: NMM grid halo width is too large (',fake_width,')',       &
      'Please note that the ARPS grid and the NMM grid must be aligned together so that', &
      'each patch does not need data more than twice of its size.',     &
      'Otherwise, you can consider program "joinwrfh" or "hdfsubdomain" for subsetting any of the grid.'
    istatus = -1
  END IF

  RETURN
END SUBROUTINE find_nmm_halo_width

!#######################################################################
!#######################################################################
!####                                                               ####
!#### Interpolation method based on EXT2ARPS & WRF2ARPS             ####
!####                                                               ####
!#######################################################################
!#######################################################################

SUBROUTINE mkarpsvlz(ims,ime,jms,jme,kms,kme,nx,ny,nz,varname,stagger,  &
                     intrpmethod,iprtopt,intropt,vlnsnd,                &
                     hgtext,za,varext, arpsbar,arpsprt,                 &
                     tem_ext, dbglvl,istatus)

  USE module_interpolation
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Take data from external file on pressure levels and interpolate in
!  the horizontal and vertical to ARPS height levels and horizontal
!  locations. Then, form the ARPS mean and perturbation quantities.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster
!  5/01/1994.
!
!  MODIFICATION HISTORY:
!  Modified from SUBROUTINE mkarpsvlz in src/ext2arps/extlib.f90 for NMM2ARPS.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx_ext   Number of grid points in the x-direction (east/west)
!             for the external grid
!    ny_ext   Number of grid points in the y-direction (north/south)
!             for the external grid
!    nz_ext   Number of grid points in the vertical
!             for the external grid
!
!    nx       Number of grid points in the x-direction (east/west)
!             for the ARPS grid
!    ny       Number of grid points in the y-direction (north/south)
!             for the ARPS grid
!    nz       Number of grid points in the vertical
!             for the ARPS grid
!
!    intrpopt Order of polynomial for interpolation (1, 2 or 3)
!    intropt  Option indicating to interpolate perturbation or total
!             variables:
!             = 1  Interpolate perturbation variables and add to base
!                  sounding (default);
!             = 2  Interploate total variables.
!    iprtopt  Flag for producing a perturbation variable
!             iprtopt = 1   Produce mean and perturbation field
!             iprtopt = 0   Produce mean and total field
!
!    iloc     x-index location of ARPS grid point in the external array
!    jloc     y-index location of ARPS grid point in the external array
!
!    x2d      x coordinate of ARPS grid point in external coordinate
!    y2d      x coordinate of ARPS grid point in external coordinate
!
!    zext     Array of heights of the external grid
!    za       Array of heights of ARPS physical heights
!    zsnd     1-D array of height representing a mean sounding
!             over domain
!    vlnsnd   1-D array of variable representing the log of the
!             mean sounding of the variable.
!
!  OUTPUT:
!
!    arpsbar  3-D array of mean field on ARPS levels
!    arpsprt  3-D array of perturbation (iprtopt=1) or
!             total (iprtopt=0) variable at ARPS grid locations
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
!  Input variables
!
!-----------------------------------------------------------------------
!
  INTEGER, INTENT(IN) :: ims,ime,jms,jme,kms,kme  ! external grid dimensions
  INTEGER, INTENT(IN) :: nx,ny,nz              ! ARPS grid dimensions
  CHARACTER(LEN=*), INTENT(IN) :: varname
  INTEGER, INTENT(IN) :: stagger
  INTEGER, INTENT(IN) :: intrpmethod           ! interpolation option
  INTEGER, INTENT(IN) :: iprtopt               ! perturbation generation flag
  INTEGER, INTENT(IN) :: intropt               ! interpolation option

  REAL,    INTENT(IN) :: vlnsnd(lvlprof)       ! 1-D array of level-means

  REAL,    INTENT(IN) :: hgtext(ims:ime,jms:jme,kms:kme)     ! heights of external levels
  REAL,    INTENT(IN) :: za(nx,ny,nz)             ! ARPS physical heights
  REAL,    INTENT(IN) :: varext(ims:ime,jms:jme,kms:kme)     ! variable to convert
!
!-----------------------------------------------------------------------
!
!  Output variables
!
!-----------------------------------------------------------------------
!
  REAL,    INTENT(OUT) :: arpsbar( nx, ny, nz)     ! 3-D array of level-means
  REAL,    INTENT(OUT) :: arpsprt( nx, ny, nz)     ! Output array, perturnbation variable
                                   ! or total variable (see iprtopt)
  INTEGER, INTENT(IN)  :: dbglvl
  INTEGER, INTENT(OUT) :: istatus
!
!-----------------------------------------------------------------------
!
!  Temporary work arrays
!
!-----------------------------------------------------------------------
!
  REAL, INTENT(INOUT) :: tem_ext(ims:ime,jms:jme,kms:kme)

  REAL, POINTER :: zext(:,:,:)      ! heights of external levels in ARPS horizontal grid

  REAL, POINTER :: iloc(:,:)        ! external x-index of ARPS grid point
  REAL, POINTER :: jloc(:,:)        ! external y-index of ARPS grid point

!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k,ia,ja,ka,kl
  REAL    :: wlow
  REAL    :: topprt,botprt,arpstop,arpsbot
  REAL    :: intrp_pnt_from_2d
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  istatus = 0

!
!-----------------------------------------------------------------------
!
!  Subtract horizontal mean from external fields.
!
!-----------------------------------------------------------------------
!
  IF (intropt == 1) THEN  ! tem_ext is now pressure perturbation
    DO j=jms,jme
      DO i=ims,ime
        DO k=kms,kme
          DO kl=2,lvlprof-1
            IF(zsnd(kl) > hgtext(i,j,k)) EXIT
          END DO
          wlow=(zsnd(kl)-hgtext(i,j,k))/(zsnd(kl)-zsnd(kl-1))
          tem_ext(i,j,k) = varext(i,j,k)-EXP(                           &
                             ((1.-wlow)*vlnsnd(kl) + wlow*vlnsnd(kl-1)) )
        END DO
      END DO
    END DO
  ELSE ! intropt=2       ! tem_ext is not ln(total pressure)
    tem_ext = LOG(varext)
  END IF

!
!-----------------------------------------------------------------------
!
!  Retrieve temporary interpolation arrays from interpolation module_interpolation
!
!-----------------------------------------------------------------------
!
  zext   => zps_tmp
  iloc   => isr
  jloc   => jsr

  CALL egridfill(tem_ext,ims,ime,jms,jme,kms,kme,HH, istatus)
!
!-----------------------------------------------------------------------
!
!  Loop through all ARPS grid points
!
!-----------------------------------------------------------------------
!
  DO ja=1,ny
    DO ia=1,nx
!
!-----------------------------------------------------------------------
!
!  Interpolate from the mean sounding to get arpsbar
!
!-----------------------------------------------------------------------
!
      DO ka=1,nz
        DO kl=2,lvlprof-1
          IF(zsnd(kl) > za(ia,ja,ka)) EXIT
        END DO
        wlow=(zsnd(kl)-za(ia,ja,ka)) / (zsnd(kl)-zsnd(kl-1))
        arpsbar(ia,ja,ka)=EXP((1.-wlow)*vlnsnd(kl) + wlow*vlnsnd(kl-1))
      END DO
!
!-----------------------------------------------------------------------
!
!    Find vertical location
!    and interpolate in vertical between two horizontal
!    interpolations on the external grid surfaces.
!
!    Extrapolation is done by assuming the perturbation
!    from mean is constant.
!
!-----------------------------------------------------------------------
!
      botprt = intrp_pnt_from_2d(intrpmethod, iloc(ia,ja),jloc(ia,ja),  &
                   iEms,iEme,jms,jme,kms,kme,kms,workarr,dbglvl,istatus)

      topprt = intrp_pnt_from_2d(intrpmethod, iloc(ia,ja),jloc(ia,ja),  &
                   iEms,iEme,jms,jme,kms,kme,kme,workarr,dbglvl,istatus)

      DO ka=1,nz
        IF(za(ia,ja,ka) < zext(ia,ja,kms)) THEN
          arpsprt(ia,ja,ka)=botprt
        ELSE IF(za(ia,ja,ka) >= zext(ia,ja,kme)) THEN
          arpsprt(ia,ja,ka)=topprt
        ELSE
          DO kl=kms+1,kme
            IF(zext(ia,ja,kl) > za(ia,ja,ka)) EXIT
          END DO
          wlow = (zext(ia,ja,kl)-   za(ia,ja,ka)) /                     &
                 (zext(ia,ja,kl)-zext(ia,ja,kl-1))

          arpstop = intrp_pnt_from_2d(intrpmethod,                      &
                      iloc(ia,ja),jloc(ia,ja),iEms,iEme,jms,jme,kms,kme,&
                      kl,workarr,dbglvl,istatus)

          arpsbot = intrp_pnt_from_2d(intrpmethod,                      &
                      iloc(ia,ja),jloc(ia,ja),iEms,iEme,jms,jme,kms,kme,&
                      kl-1,workarr,dbglvl,istatus)

          arpsprt(ia,ja,ka)=(1.-wlow)*arpstop+wlow*arpsbot
        END IF
      END DO
    END DO
  END DO

  IF(intropt == 2) THEN
    arpsprt(:,:,:) = EXP(arpsprt(:,:,:))
  END IF
!
!-----------------------------------------------------------------------
!
!  If total quantity, rather than perturbation quantity desired,
!  add the level mean to the interpolated perturbation variable
!  at each grid point.
!
!-----------------------------------------------------------------------
!
  IF(iprtopt == 0 .AND. intropt == 1) THEN
    arpsprt = arpsprt + arpsbar
  ELSE IF(iprtopt == 1 .AND. intropt == 2) THEN
    arpsprt = arpsprt - arpsbar
  END IF

  RETURN
END SUBROUTINE mkarpsvlz

SUBROUTINE mkarpsvar(ims,ime,jms,jme,kms,kme,nx,ny,nz,varname, stagger, &
                     intrpmethod,iprtopt,intropt, varsnd,               &
                     hgtext,za,varext, arpsbar,arpsprt,                 &
                     tem_ext, dbglvl, istatus )

  USE module_interpolation
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Take data from external file on pressure levels and interpolate in
!  the horizontal and vertical to ARPS height levels and horizontal
!  locations. Then, form the ARPS mean and perturbation quantities.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster
!  5/01/1994.
!
!  MODIFICATION HISTORY:
!  Modified from mkarpsvar from EXT2ARPS
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx_ext   Number of grid points in the x-direction (east/west)
!             for the external grid
!    ny_ext   Number of grid points in the y-direction (north/south)
!             for the external grid
!    nz_ext   Number of grid points in the vertical
!             for the external grid
!
!    nx       Number of grid points in the x-direction (east/west)
!             for the ARPS grid
!    ny       Number of grid points in the y-direction (north/south)
!             for the ARPS grid
!    nz       Number of grid points in the vertical
!             for the ARPS grid
!
!    iorder   order of polynomial for interpolation (1, 2 or 3)
!    intropt  Option indicating to interpolate perturbation or total
!             variables:
!             = 1  Interpolate perturbation variables and add to base
!                  sounding (default);
!             = 2  Interploate total variables.
!    iprtopt  Flag for producing a perturbation variable
!             iprtopt = 1   Produce mean and perturbation field
!             iprtopt = 0   Produce mean and total field
!
!    iloc     x-index location of ARPS grid point in the external array
!    jloc     y-index location of ARPS grid point in the external array
!
!    x2d      x coordinate of ARPS grid point in external coordinate
!    y2d      x coordinate of ARPS grid point in external coordinate
!
!    zext     Array of heights of the external grid
!    za       Array of heights of ARPS physical heights
!    zsnd     1-D array of height representing a mean sounding
!             over domain
!    varsnd   1-D array of variable representing a mean sounding
!             over domain
!
!  OUTPUT:
!
!    arpsbar  3-D array of mean field on ARPS levels
!    arpsprt  3-D array of perturbation (iprtopt=1) or
!             total (iprtopt=0) variable at ARPS grid locations
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

!
!-----------------------------------------------------------------------
!
!  Input variables
!
!-----------------------------------------------------------------------
!
  INTEGER, INTENT(IN) :: ims,ime,jms,jme,kms,kme  ! external grid dimensions
  INTEGER, INTENT(IN) :: nx,ny,nz                 ! ARPS grid dimensions
  CHARACTER(LEN=*), INTENT(IN) :: varname
  INTEGER, INTENT(IN) :: stagger
  INTEGER, INTENT(IN) :: intrpmethod              ! interpolating polynomial order
  INTEGER, INTENT(IN) :: intropt                  ! interpolation option
  INTEGER, INTENT(IN) :: iprtopt                  ! perturbation generation flag
  REAL,    INTENT(IN) :: varsnd(lvlprof)          ! 1-D array of level-means
  REAL,    INTENT(IN) :: hgtext(ims:ime,jms:jme,kms:kme)     ! heights of external levels
  REAL,    INTENT(IN) :: za(nx,ny,nz)              ! ARPS physical heights
  REAL,    INTENT(IN) :: varext(ims:ime,jms:jme,kms:kme)     ! variable to convert
!
!-----------------------------------------------------------------------
!
!  Output variables
!
!-----------------------------------------------------------------------
!
  REAL,    INTENT(OUT) :: arpsbar( nx, ny, nz)     ! 3-D array of level-means
  REAL,    INTENT(OUT) :: arpsprt( nx, ny, nz)     ! Output array, perturnbation variable
                                   ! or total variable (see iprtopt)

  INTEGER, INTENT(IN)  :: dbglvl
  INTEGER, INTENT(OUT) :: istatus
!
!-----------------------------------------------------------------------
!
!  Temporary work arrays
!
!-----------------------------------------------------------------------
!
  REAL :: tem_ext(ims:ime,jms:jme,kms:kme)

  REAL, POINTER :: iloc(:,:)         ! external x-index of ARPS grid point
  REAL, POINTER :: jloc(:,:)         ! external y-index of ARPS grid point
  REAL, POINTER :: zext(:,:,:)       ! heights of external levels
                                     ! interpolated to ARPS grid locs
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k,ia,ja,ka,kl
  REAL    :: wlow
  REAL    :: topprt,botprt,arpstop,arpsbot

  REAL    :: intrp_pnt_from_2d

  INTEGER :: extstagger
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  istatus = 0
!
!-----------------------------------------------------------------------
!
!  Subtract horizontal mean from external fields.
!
!-----------------------------------------------------------------------
!
  IF (intropt == 1) THEN
    DO j=jms,jme
      DO i=ims,ime
        DO k=kms,kme
          DO kl=2,lvlprof-1
            IF(zsnd(kl) > hgtext(i,j,k)) EXIT
          END DO
          wlow = (zsnd(kl)-hgtext(i,j,k)) / (zsnd(kl)-zsnd(kl-1))
          tem_ext(i,j,k) = varext(i,j,k)                                &
                         - ((1.-wlow)*varsnd(kl) + wlow*varsnd(kl-1))
        END DO
      END DO
    END DO
  ELSE
    tem_ext(:,:,:) = varext(:,:,:)
  END IF

!
!-----------------------------------------------------------------------
!
!  Retrieve temporary interpolation arrays from interpolation module_interpolation
!
!-----------------------------------------------------------------------
!
  SELECT CASE (stagger)
  CASE (M)
    zext => zps_tmp
    iloc => isr
    jloc => jsr
    extstagger = HH
  CASE (U)
    zext => zpv_tmp
    iloc => iur
    jloc => jur
    extstagger = VV
  CASE (V)
    zext => zpv_tmp
    iloc => ivr
    jloc => jvr
    extstagger = VV
  CASE (W)
    zext => zp_tmp
    iloc => isr
    jloc => jsr
    extstagger = HH

!  CASE (SS)
!    zext => zpsoil_tmp
!    extstagger = HH
  CASE DEFAULT
    WRITE(6,'(1x,a,I4,a)') 'ERROR: unsupported stagger - ',stagger,' within mkarpsver.'
    CALL arpsstop('Unsupported stagger.',1)
  END SELECT

  CALL egridfill(tem_ext,ims,ime,jms,jme,kms,kme,extstagger, istatus)

!
!
!-----------------------------------------------------------------------
!
!  Loop through all ARPS grid points
!
!-----------------------------------------------------------------------
!
  DO ja=1,ny
    DO ia=1,nx
!
!-----------------------------------------------------------------------
!
!  Interpolate from the mean sounding to get arpsbar
!
!-----------------------------------------------------------------------
!
      DO ka=1,nz
        DO kl=2,lvlprof-1
          IF(zsnd(kl) > za(ia,ja,ka)) EXIT
        END DO
        wlow = (zsnd(kl)-za(ia,ja,ka)) / (zsnd(kl)-zsnd(kl-1))
        arpsbar(ia,ja,ka) = (1.-wlow)*varsnd(kl) + wlow*varsnd(kl-1)
      END DO
!
!-----------------------------------------------------------------------
!
!    Find vertical location
!    and interpolate in vertical between two horizontal
!    interpolations on the external grid surfaces.
!
!    Extrapolation is done by assuming the perturbation
!    from mean is constant.
!
!-----------------------------------------------------------------------
!
      botprt = intrp_pnt_from_2d(intrpmethod,                           &
                      iloc(ia,ja),jloc(ia,ja),iEms,iEme,jms,jme,kms,kme,&
                      kms,workarr,dbglvl,istatus)

      topprt = intrp_pnt_from_2d(intrpmethod,                           &
                      iloc(ia,ja),jloc(ia,ja),iEms,iEme,jms,jme,kms,kme,&
                      kme,workarr,dbglvl,istatus)

      DO ka=1,nz
        IF(za(ia,ja,ka) < zext(ia,ja,kms)) THEN
          arpsprt(ia,ja,ka)=botprt
        ELSE IF(za(ia,ja,ka) >= zext(ia,ja,kme)) THEN
          arpsprt(ia,ja,ka)=topprt
        ELSE
          DO kl=kms+1,kme
            IF(zext(ia,ja,kl) > za(ia,ja,ka)) EXIT
          END DO

          wlow=(zext(ia,ja,kl)-   za(ia,ja,ka))/                        &
               (zext(ia,ja,kl)-zext(ia,ja,kl-1))

          arpstop = intrp_pnt_from_2d(intrpmethod,                      &
                      iloc(ia,ja),jloc(ia,ja),iEms,iEme,jms,jme,kms,kme,&
                      kl,workarr,dbglvl,istatus)

          arpsbot = intrp_pnt_from_2d(intrpmethod,                      &
                      iloc(ia,ja),jloc(ia,ja),iEms,iEme,jms,jme,kms,kme,&
                      kl-1,workarr,dbglvl,istatus)

          arpsprt(ia,ja,ka)=(1.-wlow)*arpstop+wlow*arpsbot

        END IF
      END DO
    END DO
  END DO
!
!-----------------------------------------------------------------------
!
!  If total quantity, rather than perturbation quantity desired,
!  add the level mean to the interpolated perturbation variable
!  at each grid point.
!
!-----------------------------------------------------------------------
!
  IF(iprtopt == 0 .AND. intropt == 1) THEN
    arpsprt = arpsprt + arpsbar
  ELSE IF(iprtopt == 1 .AND. intropt == 2) THEN
    arpsprt = arpsprt - arpsbar
  END IF

  RETURN
END SUBROUTINE mkarpsvar

SUBROUTINE mkarps2d(nx,ny,interp_method,stagger,ims,ime,jms,jme,        &
                    varin, varout,dbglvl,istatus)

!#######################################################################

  USE module_interpolation

  IMPLICIT NONE

  INTEGER, INTENT(IN)  :: nx, ny
  INTEGER, INTENT(IN)  :: interp_method, stagger
  INTEGER, INTENT(IN)  :: ims, ime, jms, jme
  REAL,    INTENT(IN)  :: varin(ims:ime,jms:jme)
  REAL,    INTENT(OUT) :: varout(nx,ny)
  INTEGER, INTENT(IN)  :: dbglvl
  INTEGER, INTENT(OUT) :: istatus

!-----------------------------------------------------------------------

  INTEGER :: ibgn, iend, jbgn, jend
  INTEGER :: i, j

  INTEGER :: extstagger
  REAL, POINTER  :: xi(:,:), yj(:,:)

  REAL    :: interp_4pnt, interp_16pnt, interp_bilinear
  REAL    :: interp_nearneighbor

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  istatus = 0

  ibgn = 1
  iend = nx
  jbgn = 1
  jend = ny

!---------------------- Prepare for Interpolation ------------------------

  SELECT CASE (interp_method)
  CASE (NEARNEIGHBOR, NEARNEIGHBOR1, NEARNEIGHBOR2)

    xi => nmmi
    yj => nmmj

  CASE (FOUR_PNT, SIXTEEN_PNT, BILINEAR)

    SELECT CASE (stagger)
    CASE (M)
      extstagger = HH
      xi => isr
      yj => jsr
    CASE (U)
      extstagger = VV
      xi => iur
      yj => iur
    CASE (V)
      extstagger = VV
      xi => ivr
      yj => ivr
    CASE DEFAULT
      WRITE(6,'(1x,a,I4,a)') 'ERROR: unsupported stagger - ',stagger,' within mkarps2d.'
      CALL arpsstop('Unsupported stagger.',1)
    END SELECT

    CALL egridfill(varin,ims,ime,jms,jme,1,1,extstagger, istatus)

  CASE DEFAULT
    WRITE(6,'(1x,a,I2)') 'ERROR: unknown interpolation method - ',interp_method,' in mkarps2d.'
    CALL arpsstop('ERROR: unknown interpolation method.',1)
  END SELECT

!---------------------- Do Interpolation -----------------------------

  SELECT CASE (interp_method)
  CASE (NEARNEIGHBOR, NEARNEIGHBOR1, NEARNEIGHBOR2)

    DO j = jbgn, jend
      DO i = ibgn, iend
        varout(i,j) = interp_nearneighbor( xi(i,j),yj(i,j),1,varin,     &
                                   ims,ime,jms,jme,1,1,-999.0 )
      END DO
    END DO

  CASE (FOUR_PNT)
    DO j = jbgn, jend
      DO i = ibgn, iend
        varout(i,j) = interp_4pnt( xi(i,j),yj(i,j),1,workarr,           &
                                   iEms,iEme,jms,jme,1,1,-999.0 )
      END DO
    END DO

  CASE ( SIXTEEN_PNT )

    DO j = jbgn, jend
      DO i = ibgn, iend
        varout(i,j) = interp_16pnt( xi(i,j),yj(i,j),1,workarr,          &
                                    iEms,iEme,jms,jme,1,1,-999.0)

      END DO
    END DO

  CASE ( BILINEAR )

    DO j = jbgn, jend
      DO i = ibgn, iend
        varout(i,j) = interp_bilinear( xi(i,j),yj(i,j),1,workarr,       &
                                    iEms,iEme,jms,jme,1,1,-999.0)

      END DO
    END DO

  END SELECT

  RETURN
END SUBROUTINE mkarps2d

SUBROUTINE mkarps2di(nx,ny,interp_method,stagger,ims,ime,jms,jme,       &
                     varin,varout,dbglvl,istatus)

!#######################################################################

  USE module_interpolation

  IMPLICIT NONE

  INTEGER, INTENT(IN)  :: nx, ny
  INTEGER, INTENT(IN)  :: interp_method, stagger
  INTEGER, INTENT(IN)  :: ims, ime, jms, jme
  INTEGER, INTENT(IN)  :: varin(ims:ime,jms:jme)
  INTEGER, INTENT(OUT) :: varout(nx,ny)
  INTEGER, INTENT(IN)  :: dbglvl
  INTEGER, INTENT(OUT) :: istatus

!-----------------------------------------------------------------------

  INTEGER :: ibgn, iend, jbgn, jend
  INTEGER :: i, j

  INTEGER :: item
  REAL, POINTER  :: xi(:,:), yj(:,:)

  INTEGER :: interp_4pnti
  INTEGER :: interp_nearneighbori

!
! These two tables are needed to convert WRF soil and vegetation index
! to those used in ARPS, althought they may not correspond each other
! one to one.
!
! WRF uses 16 soil categories and 24-category (USGS) vegetation
! ARPS uses 13-category soil  and 14-category (ND) vegetation
!
! The following two tables were provided by Jerry Brotzge on Oct. 20, 2003
!
  INTEGER, PARAMETER :: soil_table(17) = (/ 1, 2, 3, 4, 4,              &
                                            5, 6, 7, 8, 9,              &
                                           10,11, 6,13, 1,              &
                                            2, 2/)

! WDT RLC 2004-02-12 Changed veg_table(1)/Urban from 7/EvgrnForest to 1/Desert
! WDT RLC 2004-02-12 Changed veg_table(25)/Playa from 10/Cultiv to 1/Desert
! WDT RLC 2004-02-12 Added veg_table(26:27)
! WDT RLC 2004-02-16 Not accepting CAPS changes:
                     ! WYH - (01.16.2004) veg_table(18) = 11 =>  8
                     !                    veg_table(19) = 1  => 13
  !INTEGER, PARAMETER :: veg_table(25)  = (/ 7,10,10,10,10,              &
  !                                          5, 3,12, 4,12,              &
  !                                          6, 6, 7, 7, 6,              &
  !                                         14,11,11, 1, 2,              &
  !                                          2, 2, 2, 9,10/)
  INTEGER, PARAMETER :: veg_table(27)  = (/ 1,10,10,10,10,              &
                                            5, 3,12, 4,12,              &
                                            6, 6, 7, 7, 6,              &
                                           14,11,11, 1, 2,              &
                                            2, 2, 2, 9, 1,              &
                                            1, 1 /)

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  istatus = 0

  ibgn = 1
  iend = nx
  jbgn = 1
  jend = ny

  xi => nmmi
  yj => nmmj

!-----------------------------------------------------------------------

  SELECT CASE (interp_method)
  CASE (NEARNEIGHBOR)

    DO j = jbgn, jend
      DO i = ibgn, iend
        varout(i,j) = interp_nearneighbori( xi(i,j),yj(i,j),1,varin,    &
                                   ims,ime,jms,jme,1,1,-999 )
      END DO
    END DO

  CASE ( NEARNEIGHBOR1 )    ! Convert soil type to ARPS

    DO j = jbgn, jend
      DO i = ibgn, iend
        item = interp_nearneighbori( xi(i,j),yj(i,j),1,varin,           &
                                   ims,ime,jms,jme,1,1,-999 )
        varout(i,j) = soil_table ( item )
      END DO
    END DO

  CASE ( NEARNEIGHBOR2 )    ! Convert veg type to ARPS

    DO j = jbgn, jend
      DO i = ibgn, iend
        item = interp_nearneighbori( xi(i,j),yj(i,j),1,varin,           &
                                   ims,ime,jms,jme,1,1,-999 )
        varout(i,j) = veg_table ( item )
      END DO
    END DO

  CASE (FOUR_PNT)
!    DO j = jbgn, jend
!      DO i = ibgn, iend
!        varout(i,j) = interp_4pnt( xi(i,j),yj(i,j),1,varin,             &
!                                   ims,ime,jms,jme,1,1,-999.0 )
!      END DO
!    END DO
    WRITE(6,'(1x,a,I2)') 'ERROR: interpolation method is still not implemented - ',interp_method
    CALL arpsstop('ERROR: not implemented interpolation method.',1)

  CASE ( SIXTEEN_PNT )
    WRITE(6,'(1x,a,I2)') 'ERROR: interpolation method is still not implemented - ',interp_method
    CALL arpsstop('ERROR: not implemented interpolation method.',1)
  CASE DEFAULT
    WRITE(6,'(1x,a,I2)') 'ERROR: unknown interpolation method - ',interp_method
    CALL arpsstop('ERROR: unknown interpolation method.',1)
  END SELECT
  RETURN
END SUBROUTINE mkarps2di

SUBROUTINE vintrpsoil_wrf(nx,ny,nzsoil_ext,nzsoil,soilmodel_option,     &
                         zpsoilext,zpsoil,tsoil_in,qsoil_in,            &
                         tsoil_out,qsoil_out)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Vertically interpolate WRF soil variables which is in the same horizontal
!  grid as ARPS data to ARPS soil layers
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yunheng Wang
!  9/20/2003. Based on mkarpsvar and intrp_soil etc.
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!             for the ARPS grid
!    ny       Number of grid points in the y-direction (north/south)
!             for the ARPS grid
!
!    nzsoil_ext   Number of grid points in the vertical soil layers
!                 for the external grid
!
!    nzsoil       Number of grid points in the vertical soil layers
!                 for the ARPS grid
!
!    soilmodel_option  Flag for ARPS soil model scheme
!                     = 1   Two-layer Force-restore model (Noilhan/Planton scheme)
!                     = 2   Multi-layer 'OUSoil' scheme (Based on OSU/NCEP ETA scheme)
!
!    zpsoilext     Array of heights of the external grid
!    zpsoil        Array of heights of ARPS physical heights
!
!    tsoil_in
!    qsoil_in
!
!  OUTPUT:
!
!    tsoil_out     3D array for soil temperature in ARPS grid
!    qsoil_out     3D array for soil moisture in ARPS grid
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
!  Input variables
!
!-----------------------------------------------------------------------
!
  INTEGER, INTENT(IN)    :: nx,ny           ! ARPS grid dimensions
  INTEGER, INTENT(IN)    :: nzsoil_ext      ! external grid dimensions
  INTEGER, INTENT(INOUT) :: nzsoil       ! ARPS grid dimensions
  INTEGER, INTENT(IN)    :: soilmodel_option
  REAL,    INTENT(IN)    :: zpsoilext(nzsoil_ext)
  REAL,    INTENT(INOUT) :: zpsoil(nx,ny,nzsoil)
  REAL,    INTENT(IN)    :: tsoil_in(nx,ny,nzsoil_ext)
  REAL,    INTENT(IN)    :: qsoil_in(nx,ny,nzsoil_ext)
!
!-----------------------------------------------------------------------
!
!  Output variables
!
!-----------------------------------------------------------------------
!
  REAL,    INTENT(OUT)   :: tsoil_out(nx,ny,nzsoil)
  REAL,    INTENT(OUT)   :: qsoil_out(nx,ny, nzsoil)
!
!-----------------------------------------------------------------------
!
!  Temporary work arrays
!
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k,m
  REAL    :: w1,w2, dist
  REAL    :: w(nzsoil_ext), wtotal

!  include 'mp.inc'
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  IF (soilmodel_option == 1) THEN    ! old ARPS Force-Restore Soil model
    zpsoil(:,:,1) = 0.05
    zpsoil(:,:,2) = 0.55
    nzsoil = 2

    tsoil_out(:,:,1) = tsoil_in(:,:,1)
    qsoil_out(:,:,1) = qsoil_in(:,:,1)

    wtotal = 0.0                    ! The second layer is weighted average
    DO m = 1, nzsoil_ext
      w(m) = MAX(0.0, 1.0-abs(zpsoilext(m)-0.55))
      wtotal = wtotal + w(m)
    END DO

    DO j=1,ny
      DO i=1,nx
        tsoil_out(i,j,2) = 0.0
        qsoil_out(i,j,2) = 0.0

        DO m = 1, nzsoil_ext
          tsoil_out(i,j,2) = tsoil_out(i,j,2) + w(m)*tsoil_in(i,j,m)/wtotal
          qsoil_out(i,j,2) = qsoil_out(i,j,2) + w(m)*qsoil_in(i,j,m)/wtotal
        END DO
      END DO
    END DO
  ELSE
!
!-----------------------------------------------------------------------
!
!  Loop through all ARPS grid points
!
!-----------------------------------------------------------------------
!
    DO j=1,ny
      DO i=1,nx
        DO k = 1, nzsoil
          IF(zpsoil(i,j,k) < zpsoilext(1)) THEN            ! extrapolation

            dist = zpsoilext(1) - zpsoil(i,j,k)
            w1 = (tsoil_in(i,j,1) - tsoil_in(i,j,2)) /                  &
                                           (zpsoilext(2)- zpsoilext(1))
            w2 = (qsoil_in(i,j,1) - qsoil_in(i,j,2)) /                  &
                                           (zpsoilext(2)- zpsoilext(1))
            tsoil_out(i,j,k) = tsoil_in(i,j,1) + w1*dist
            qsoil_out(i,j,k) = qsoil_in(i,j,1) + w2*dist

          ELSE IF(zpsoil(i,j,k) >= zpsoilext(nzsoil_ext)) THEN ! constant layers

            tsoil_out(i,j,k) = tsoil_in(i,j,nzsoil_ext)
            qsoil_out(i,j,k) = qsoil_in(i,j,nzsoil_ext)

          ELSE                         ! linear interpolation

            DO m = 1, nzsoil_ext
              IF(zpsoilext(m) > zpsoil(i,j,k)) EXIT
            END DO
            w1 = (zpsoilext(m) - zpsoil(i,j,k)) /                       &
                                      (zpsoilext(m)- zpsoilext(m-1))
            w2 = (zpsoil(i,j,k) - zpsoilext(m-1) ) /                &
                                      (zpsoilext(m)- zpsoilext(m-1))
            tsoil_out(:,:,k) = tsoil_in(:,:,m-1)*w1 + tsoil_in(:,:,m)*w2
            qsoil_out(:,:,k) = qsoil_in(:,:,m-1)*w1 + qsoil_in(:,:,m)*w2
          END IF

        END DO
      END DO
    END DO

  END IF

  RETURN
END SUBROUTINE vintrpsoil_wrf

REAL FUNCTION intrp_pnt_from_2d(interp_method,xx,yy,ims,ime,jms,jme,    &
                                kms,kme,izz, varin, dbglvl,istatus)

!#######################################################################

  USE module_wrfgrid_constants

  IMPLICIT NONE

  INTEGER, INTENT(IN)  :: interp_method
  INTEGER, INTENT(IN)  :: ims, ime, jms, jme
  INTEGER, INTENT(IN)  :: kms, kme, izz
  REAL,    INTENT(IN)  :: xx,yy
  REAL,    INTENT(IN)  :: varin(ims:ime,jms:jme,kms:kme)
  INTEGER, INTENT(IN)  :: dbglvl
  INTEGER, INTENT(OUT) :: istatus

!-----------------------------------------------------------------------
  REAL :: interp_4pnt, interp_16pnt, interp_bilinear

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  istatus = 0

!-----------------------------------------------------------------------

  SELECT CASE (interp_method)
  CASE (NEARNEIGHBOR, NEARNEIGHBOR1, NEARNEIGHBOR2)
    WRITE(6,'(1x,a,I2)') 'ERROR: interpolation method is still not implemented - ',interp_method
    CALL arpsstop('ERROR: not implemented interpolation method.',1)

  CASE (FOUR_PNT)
    intrp_pnt_from_2d = interp_4pnt(xx,yy,izz,varin,                    &
                                    ims,ime,jms,jme,kms,kme,-999.0)

  CASE ( SIXTEEN_PNT )
    intrp_pnt_from_2d = interp_16pnt(xx,yy,izz,varin,                   &
                                    ims,ime,jms,jme,kms,kme,-999.0)

  CASE ( BILINEAR )
    intrp_pnt_from_2d = interp_bilinear(xx,yy,izz,varin,                &
                                    ims,ime,jms,jme,kms,kme,-999.0)

  CASE DEFAULT
    WRITE(6,'(1x,a,I2)') 'ERROR: unknown interpolation method - ',interp_method
    CALL arpsstop('ERROR: unknown interpolation method.',1)
  END SELECT
  RETURN
END FUNCTION intrp_pnt_from_2d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!                                                                  !!!
!!! Interpolation method modified based on those on WPS              !!!
!!!                                                                  !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

REAL FUNCTION interp_4pnt(xx, yy, izz, array,                           &
                          start_x,end_x, start_y,end_y,start_z,end_z,   &
                          msgval )

!#######################################################################
! Inherited from function wt_four_pt_average in WPS.

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: start_x, start_y, end_x, end_y
  INTEGER, INTENT(IN) :: start_z, end_z
  INTEGER, INTENT(IN) :: izz
  REAL,    INTENT(IN) :: xx, yy
  REAL,    INTENT(IN) :: msgval
  REAL,    INTENT(IN) :: array(start_x:end_x, start_y:end_y, start_z:end_z)

!  REAL,    INTENT(IN), OPTIONAL :: mask_array(start_x:end_x, start_y:end_y)
!  REAL,    INTENT(IN), OPTIONAL :: maskval

!----------------------------------------------------------------------

  INTEGER :: istatus
  INTEGER :: ifx, ify, icx, icy
  REAL    :: fxfy, fxcy, cxfy, cxcy

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  istatus = 0

  ifx = FLOOR(xx)
  icx = CEILING(xx)
  ify = FLOOR(yy)
  icy = CEILING(yy)

  fxfy = max(0., 1.0 - sqrt((xx-real(ifx))**2+(yy-real(ify))**2))
  fxcy = max(0., 1.0 - sqrt((xx-real(ifx))**2+(yy-real(icy))**2))
  cxfy = max(0., 1.0 - sqrt((xx-real(icx))**2+(yy-real(ify))**2))
  cxcy = max(0., 1.0 - sqrt((xx-real(icx))**2+(yy-real(icy))**2))

  ! First, make sure that the point is contained in the source array
  if (ifx < start_x .or. icx > end_x .or. &
      ify < start_y .or. icy > end_y) then

     ! But if the point is at most half a grid point out, we can
     !   still proceed with modified ifx, icx, ify, and icy.
     if (xx > real(start_x)-0.5 .and. ifx < start_x) then
        ifx = start_x
        icx = start_x
     else if (xx < real(end_x)+0.5 .and. icx > end_x) then
        ifx = end_x
        icx = end_x
     end if

     if (yy > real(start_y)-0.5 .and. ifx < start_y) then
        ify = start_y
        icy = start_y
     else if (yy < real(end_y)+0.5 .and. icy > end_y) then
        ify = end_y
        icy = end_y
     end if

     if (ifx < start_x .or. icx > end_x .or. &
         ify < start_y .or. icy > end_y) then
        interp_4pnt = msgval         ! ext2arps missing value
        WRITE(6,'(1x,a,/,8x,2(a,F12.3),/,8x,a,4(I6,a))')                &
             'ERROR: Target grid point is outside of the source array', &
             'Target point xx = ',xx,', yy = ',yy,                      &
             'Source array start_x = ',start_x,', end_x = ',end_x,      &
                        '; start_y = ',start_y,', end_y = ',end_y,'.'
        istatus = -1
        return
     end if
  end if

!  if (PRESENT(mask_array) .AND. PRESENT(maskval)) then
!    write(0,* ) 'here',PRESENT(mask_array),PRESENT(maskval)
!      if (array(ifx, ify, izz) == msgval .or. mask_array(ifx,ify) == maskval) fxfy = 0.0
!      if (array(ifx, icy, izz) == msgval .or. mask_array(ifx,icy) == maskval) fxcy = 0.0
!      if (array(icx, ify, izz) == msgval .or. mask_array(icx,ify) == maskval) cxfy = 0.0
!      if (array(icx, icy, izz) == msgval .or. mask_array(icx,icy) == maskval) cxcy = 0.0
!  else
!      if (array(ifx, ify, izz) == msgval) fxfy = 0.0
!      if (array(ifx, icy, izz) == msgval) fxcy = 0.0
!      if (array(icx, ify, izz) == msgval) cxfy = 0.0
!      if (array(icx, icy, izz) == msgval) cxcy = 0.0
!  end if

  ! If all four points are missing, try the next interpolation method in the sequence
  if (fxfy == 0.0 .and. fxcy == 0.0 .and. cxfy == 0.0 .and. cxcy == 0.0) then
!      interp_4pnt = interp_sequence(xx, yy, izz, array, start_x, end_x, start_y, end_y, &
!                           start_z, end_z, msgval, interp_list, idx, maskval, mask_array)
      WRITE(6,'(1x,a)') 'ERROR: Source array are all missing.'
      istatus = -2
      return
  else
      interp_4pnt = (fxfy * array(ifx, ify, izz) +                      &
                     fxcy * array(ifx, icy, izz) +                      &
                     cxfy * array(icx, ify, izz) +                      &
                     cxcy * array(icx, icy, izz) ) / (fxfy + fxcy + cxfy + cxcy)
  end if

  RETURN
END FUNCTION interp_4pnt

REAL FUNCTION interp_nearneighbor(xx,yy,izz, array, start_x,end_x,      &
                 start_y,end_y,start_z,end_z,msgval)

!#######################################################################
! Inherited from function wt_four_pt_average in WPS.

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: start_x, start_y, end_x, end_y
  INTEGER, INTENT(IN) :: start_z, end_z
  INTEGER, INTENT(IN) :: izz
  REAL,    INTENT(IN) :: xx, yy
  REAL,    INTENT(IN) :: msgval
  REAL,    INTENT(IN) :: array(start_x:end_x, start_y:end_y, start_z:end_z)

!  REAL,    INTENT(IN), OPTIONAL :: mask_array(start_x:end_x, start_y:end_y)
!  REAL,    INTENT(IN), OPTIONAL :: maskval
!
!----------------------------------------------------------------------

  INTEGER :: istatus
  INTEGER :: ix, iy

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  istatus = 0

  ix = nint(xx)
  iy = nint(yy)

  ! The first thing to do is to ensure that the point (xx,yy) is within the array
  if (ix < start_x .OR. ix > end_x) then
     interp_nearneighbor = msgval
     istatus = -1
     return
  end if

  if (iy < start_y .OR. iy > end_y) then
     interp_nearneighbor = msgval
     istatus = -1
     return
  end if

!  if (present(mask_array) .AND. present(maskval)) then
!     if (mask_array(ix,iy) == maskval) then
!        interp_nearneighbor = msgval
!     else
!        interp_nearneighbor = array(ix,iy,izz)
!     end if
!  else
     interp_nearneighbor = array(ix,iy,izz)
!  end if

  RETURN
END FUNCTION interp_nearneighbor

INTEGER FUNCTION interp_nearneighbori(xx,yy,izz, array, start_x,end_x,  &
                 start_y,end_y,start_z,end_z, msgval )

!#######################################################################
! Inherited from function wt_four_pt_average in WPS.

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: start_x, start_y, end_x, end_y
  INTEGER, INTENT(IN) :: start_z, end_z
  INTEGER, INTENT(IN) :: izz
  REAL,    INTENT(IN) :: xx, yy
  INTEGER, INTENT(IN) :: msgval
  INTEGER, INTENT(IN) :: array(start_x:end_x, start_y:end_y, start_z:end_z)

!  REAL,    INTENT(IN), OPTIONAL :: mask_array(start_x:end_x, start_y:end_y)
!  REAL,    INTENT(IN), OPTIONAL :: maskval

!----------------------------------------------------------------------

  INTEGER :: istatus
  INTEGER :: ix, iy

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  istatus = 0

  ix = nint(xx)
  iy = nint(yy)

  ! The first thing to do is to ensure that the point (xx,yy) is within the array
  if (ix < start_x .OR. ix > end_x) then
     interp_nearneighbori = msgval
     istatus = -1
     return
  end if

  if (iy < start_y .OR. iy > end_y) then
     interp_nearneighbori = msgval
     istatus = -1
     return
  end if

!  if (present(mask_array) .AND. present(maskval)) then
!     if (mask_array(ix,iy) == maskval) then
!        interp_nearneighbori = msgval
!     else
!        interp_nearneighbori = array(ix,iy,izz)
!     end if
!  else
     interp_nearneighbori = array(ix,iy,izz)
!  end if
!
  RETURN
END FUNCTION interp_nearneighbori

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Name: sixteen_pt
!
! Purpose: Overlapping parabolic interpolation among sixteen grid values
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL FUNCTION interp_16pnt(xx, yy, izz, array,                          &
                         start_x, end_x, start_y, end_y,start_z, end_z, &
                         msgval )

   IMPLICIT NONE

   ! Arguments
   integer, intent(in) :: izz    ! z-index of 2d-array to interpolate within
   integer, intent(in) :: start_x, start_y, start_z
   integer, intent(in) :: end_x, end_y, end_z
   real,    intent(in) :: xx , yy              ! The location to interpolate to
   real,    intent(in) :: msgval
!   real, intent(in), optional :: maskval
   real,    intent(in) :: array(start_x:end_x, start_y:end_y, start_z:end_z)

!   integer, dimension(:), intent(in) :: interp_list
!   integer, intent(in) :: idx
!   real, dimension(start_x:end_x, start_y:end_y), &
!       intent(in), optional :: mask_array

   ! Return value
!   real :: interp_16pnt

   REAL :: oned
!-----------------------------------------------------------------------

   ! Local variables
   integer :: n , i , j , k , kk , l , ll
   real    :: x , y , a , b , c , d , e , f , g , h
   real    :: stl(4,4)
   logical :: is_masked

   INTEGER :: istatus

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

   is_masked = .false.

   if (int(xx) < start_x .or. int(xx) > end_x .or. &
       int(yy) < start_y .or. int(yy) > end_y) then

      interp_16pnt = msgval         ! ext2arps missing value
      WRITE(6,'(1x,a,/,8x,2(a,F12.3),/,8x,a,4(I6,a))')                &
           'ERROR: Target grid point is outside of the source array', &
           'Target point xx = ',xx,', yy = ',yy,                      &
           'Source array start_x = ',start_x,', end_x = ',end_x,      &
                      '; start_y = ',start_y,', end_y = ',end_y,'.'
      istatus = -1
      return

   end if

   interp_16pnt = 0.0
   n = 0
   i = int(xx + 0.00001)
   j = int(yy + 0.00001)
   x = xx - i
   y = yy - j

   if ( ( abs(x) > 0.0001 ) .or. ( abs(y) > 0.0001 ) ) then

      loop_1 : do k = 1,4
         kk = i + k - 2
         if ( kk < start_x) then
            kk = start_x
         else if ( kk > end_x) then
            kk = end_x
         end if
         loop_2 : do l = 1,4
            stl(k,l) = 0.
            ll = j + l - 2
            if ( ll < start_y ) then
               ll = start_y
            else if ( ll > end_y) then
               ll = end_y
            end if
            stl(k,l) = array(kk,ll,izz)
            n = n + 1
!            if (present(mask_array) .and. present(maskval)) then
!               if (mask_array(kk,ll) == maskval) is_masked = .true.
!            end if
            if ( stl(k,l) == 0. .and. msgval /= 0.) then
               stl(k,l) = 1.E-20
            end if
         end do loop_2
      end do loop_1

      ! If we have a missing value, try the next interpolation method in the sequence
!      if (present(mask_array) .and. present(maskval)) then
!         do k=1,4
!            do l=1,4
!               if (stl(k,l) == msgval .or. is_masked) then
!                  interp_16pnt = interp_sequence(xx, yy, izz, array, start_x, end_x, start_y, end_y, &
!                                               start_z, end_z, msgval, interp_list, idx, maskval, mask_array)
!                  return
!               end if
!            end do
!         end do
!      else
!         do k=1,4
!            do l=1,4
!               if (stl(k,l) == msgval) then
!                  interp_16pnt = interp_sequence(xx, yy, izz, array, start_x, end_x, start_y, end_y, &
!                                               start_z, end_z, msgval, interp_list, idx)
!                  return
!               end if
!            end do
!         end do
!      end if

      a = oned(x,stl(1,1),stl(2,1),stl(3,1),stl(4,1))
      b = oned(x,stl(1,2),stl(2,2),stl(3,2),stl(4,2))
      c = oned(x,stl(1,3),stl(2,3),stl(3,3),stl(4,3))
      d = oned(x,stl(1,4),stl(2,4),stl(3,4),stl(4,4))
      interp_16pnt = oned(y,a,b,c,d)

      if (n /= 16) then
         e = oned(y,stl(1,1),stl(1,2),stl(1,3),stl(1,4))
         f = oned(y,stl(2,1),stl(2,2),stl(2,3),stl(2,4))
         g = oned(y,stl(3,1),stl(3,2),stl(3,3),stl(3,4))
         h = oned(y,stl(4,1),stl(4,2),stl(4,3),stl(4,4))
         interp_16pnt = (interp_16pnt+oned(x,e,f,g,h)) * 0.5
      end if

      if (interp_16pnt == 1.E-20) interp_16pnt = 0.

   else
      if (i >= start_x .and. i <= end_x .and. j >= start_y .and. j <= end_y) then
         interp_16pnt = array(i,j,izz)
      else
         interp_16pnt = msgval
      end if
   end if

END FUNCTION interp_16pnt

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Name: oned
!
! Purpose: 1-dimensional overlapping parabolic interpolation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL FUNCTION oned(x,a,b,c,d)

   implicit none

   ! Arguments
   real, intent(in) :: x,a,b,c,d

   ! Return value
!   real :: oned

   oned = 0.

   if ( x == 0. ) then
      oned = b
   else if ( x == 1. ) then
      oned = c
   end if

   if (b*c /= 0.) then
      if ( a*d == 0. ) then
         if ( ( a == 0 ) .and. ( d == 0 ) ) then
            oned = b*(1.0-x)+c*x
         else if ( a /= 0. ) then
            oned = b+x*(0.5*(c-a)+x*(0.5*(c+a)-b))
         else if ( d /= 0. ) then
            oned = c+(1.0-x)*(0.5*(b-d)+(1.0-x)*(0.5*(b+d)-c))
         end if
      else
         oned = (1.0-x)*(b+x*(0.5*(c-a)+x*(0.5*(c+a)-b)))+x*(c+(1.0-x)*(0.5*(b-d)+(1.0-x)*(0.5*(b+d)-c)))
      end if
   end if

END FUNCTION oned

REAL FUNCTION interp_bilinear(xx, yy, izz, array,                           &
                        start_x,end_x, start_y,end_y,start_z,end_z,   &
                        msgval )

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: start_x, start_y, end_x, end_y
  INTEGER, INTENT(IN) :: start_z, end_z
  INTEGER, INTENT(IN) :: izz
  REAL,    INTENT(IN) :: xx, yy
  REAL,    INTENT(IN) :: msgval
  REAL,    INTENT(IN) :: array(start_x:end_x, start_y:end_y, start_z:end_z)


!----------------------------------------------------------------------

  INTEGER :: istatus
  INTEGER :: il, ir, jt, jb, ratlr, ratbt
  REAL    :: fac1, fac2, fac3, fac4
  REAL    :: wk1, wk2, wk3, wk4

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  istatus = 0

  il=FLOOR(xx)
  ir=il+1
  jb=FLOOR(yy)
  jt=jb+1
  ratlr=xx-il
  ratbt=yy-jb

  if (il < start_x .OR. il > end_x .OR. ir < start_x .OR. ir > end_x .OR. &
      jb < start_y .OR. jb > end_y .OR. jt < start_y .OR. jt > end_y) THEN

     interp_bilinear = msgval

  else

     fac1=(1.-ratlr)*(   ratbt)
     fac2=(   ratlr)*(   ratbt)
     fac3=(1.-ratlr)*(1.-ratbt)
     fac4=(   ratlr)*(1.-ratbt)

     wk1=array(il,jt,izz)
     wk2=array(ir,jt,izz)
     wk3=array(il,jb,izz)
     wk4=array(ir,jb,izz)

     interp_bilinear = fac1*wk1+fac2*wk2+fac3*wk3+fac4*wk4

  end if

  RETURN
END FUNCTION interp_bilinear
