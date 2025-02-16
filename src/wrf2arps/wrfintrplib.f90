!
!##################################################################
!##################################################################
!######                                                      ######
!######               SUBROUTINE VINTRPVAR_WRF               ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE vintrpvar_wrf(nx,ny,nz_ext,nz,kbgn,kend,lvlprof,             &
                         iprtopt,intropt,zext,za,varext,                &
                         zsnd,varsnd,arpsbar,arpsprt,tem_ext)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Vertically interpolate WRF data which is in the same horizontal
!  grid as ARPS data to ARPS height levels. External data (WRF) is
!  valid vertically from kbgn to kend
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yunheng Wang
!  9/20/2003 Changed from mkarpsvar in src/ext2arps/extlib.f90.
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
!    nz_ext   Number of grid points in the vertical
!             for the external grid
!
!    nz       Number of grid points in the vertical
!             for the ARPS grid
!
!    intropt  Option indicating to interpolate perturbation or total
!             variables:
!             = 1  Interpolate perturbation variables and add to base
!                  sounding (default);
!             = 2  Interploate total variables.
!    iprtopt  Flag for producing a perturbation variable
!             iprtopt = 1   Produce mean and perturbation field
!             iprtopt = 0   Produce mean and total field
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
  INTEGER, INTENT(IN)  :: nx,ny                 ! ARPS grid dimensions
  INTEGER, INTENT(IN)  :: nz_ext                ! external grid dimensions
  INTEGER, INTENT(IN)  :: nz                    ! ARPS grid dimensions
  INTEGER, INTENT(IN)  :: kbgn
  INTEGER, INTENT(IN)  :: kend
  INTEGER, INTENT(IN)  :: lvlprof               ! levels in mean sounding
  INTEGER, INTENT(IN)  :: intropt               ! interpolation option
  INTEGER, INTENT(IN)  :: iprtopt               ! perturbation generation flag
  REAL,    INTENT(IN)  :: zext(nx,ny,nz_ext)    ! heights of external levels
                                                ! interpolated to ARPS grid locs
  REAL,    INTENT(IN)  :: za(nx,ny,nz)          ! ARPS physical heights
  REAL,    INTENT(IN)  :: varext(nx,ny,nz_ext)  ! variable to convert
  REAL,    INTENT(IN)  :: zsnd(lvlprof)         ! 1-D array of level heights
  REAL,    INTENT(IN)  :: varsnd(lvlprof)       ! 1-D array of level-means
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

  REAL,    INTENT(INOUT) :: tem_ext(nx,ny,nz_ext)
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k,ka,kl
  REAL    :: wlow
  REAL    :: topprt,botprt,arpstop,arpsbot
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

!
!-----------------------------------------------------------------------
!
!  Subtract horizontal mean from external fields.
!
!-----------------------------------------------------------------------
!
  IF (intropt == 1) THEN
    DO j=1,ny
      DO i=1,nx
        DO k=1,nz_ext
          DO kl=2,lvlprof-1
            IF(zsnd(kl) > zext(i,j,k)) EXIT
          END DO
          wlow=(zsnd(kl)-zext(i,j,k)) / (zsnd(kl)-zsnd(kl-1))
          tem_ext(i,j,k)= varext(i,j,k)-                                &
                          ((1.-wlow)*varsnd(kl) + wlow*varsnd(kl-1))
        END DO
      END DO
    END DO
  ELSE
    tem_ext(:,:,:) = varext(:,:,:)
  ENDIF
!
!-----------------------------------------------------------------------
!
!  Loop through all ARPS grid points
!
!-----------------------------------------------------------------------
!
  DO j = 1,ny
    DO i = 1,nx
!
!-----------------------------------------------------------------------
!
!  Interpolate from the mean sounding to get arpsbar
!
!-----------------------------------------------------------------------
!
      DO ka=1,nz
        DO kl=2,lvlprof-1
          IF(zsnd(kl) > za(i,j,ka)) EXIT
        END DO
        wlow = (zsnd(kl)-za(i,j,ka)) / (zsnd(kl)-zsnd(kl-1))
        arpsbar(i,j,ka)= (1.-wlow)*varsnd(kl) + wlow*varsnd(kl-1)
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
      botprt = tem_ext(i,j,kbgn)
      topprt = tem_ext(i,j,kend)

      DO ka=1,nz
        IF(za(i,j,ka) < zext(i,j,kbgn)) THEN
          arpsprt(i,j,ka) = botprt
        ELSE IF(za(i,j,ka) > zext(i,j,kend)) THEN
          arpsprt(i,j,ka) = topprt
        ELSE
          DO kl=kbgn+1,kend-1
            IF(zext(i,j,kl) > za(i,j,ka)) EXIT
          END DO
          wlow = (zext(i,j,kl)-  za(i,j,ka))/                        &
                 (zext(i,j,kl)-zext(i,j,kl-1))

          arpstop = tem_ext(i,j,kl)

          arpsbot = tem_ext(i,j,kl-1)

          arpsprt(i,j,ka)=(1.-wlow)*arpstop+wlow*arpsbot

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
  IF(iprtopt == 0 .and. intropt == 1) THEN
    arpsprt(:,:,:) = arpsprt(:,:,:) + arpsbar(:,:,:)
  ELSE IF(iprtopt == 1 .and. intropt == 2) THEN
    arpsprt(:,:,:) = arpsprt(:,:,:) - arpsbar(:,:,:)
  END IF

  RETURN
END SUBROUTINE vintrpvar_wrf
!
!##################################################################
!##################################################################
!######                                                      ######
!######               SUBROUTINE VINTRPVLZ_WRF               ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE vintrpvlz_wrf(nx,ny,nz_ext,nz,lvlprof,iprtopt,intropt,       &
                         zext,za,varext,zsnd,vlnsnd,arpsbar,arpsprt,    &
                         tem_ext)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Vertically interpolate WRF pressure which is in the same horizontal
!  grid as ARPS data to ARPS height levels
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yunheng Wang
!  9/20/2003. Changed from mkarpsvlz in src/ext2arps/extlib.f90
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
!    nz_ext   Number of grid points in the vertical
!             for the external grid
!
!    nz       Number of grid points in the vertical
!             for the ARPS grid
!
!    intropt  Option indicating to interpolate perturbation or total
!             variables:
!             = 1  Interpolate perturbation variables and add to base
!                  sounding (default);
!             = 2  Interploate total variables.
!    iprtopt  Flag for producing a perturbation variable
!             iprtopt = 1   Produce mean and perturbation field
!             iprtopt = 0   Produce mean and total field
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
  INTEGER, INTENT(IN)  :: nx,ny           ! ARPS grid dimensions
  INTEGER, INTENT(IN)  :: nz_ext          ! external grid dimensions
  INTEGER, INTENT(IN)  :: nz              ! ARPS grid dimensions
  INTEGER, INTENT(IN)  :: lvlprof         ! levels in mean sounding
  INTEGER, INTENT(IN)  :: intropt         ! interpolation option
  INTEGER, INTENT(IN)  :: iprtopt         ! perturbation generation flag
  REAL,    INTENT(IN)  :: zext(nx,ny,nz_ext)  ! heights of external levels
  REAL,    INTENT(IN)  :: za(nx,ny,nz)        ! ARPS physical heights
  REAL,    INTENT(IN)  :: varext(nx,ny,nz_ext)! variable to convert
  REAL,    INTENT(IN)  :: zsnd(lvlprof)       ! 1-D array of level heights
  REAL,    INTENT(IN)  :: vlnsnd(lvlprof)     ! 1-D array of level-means
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
!
!-----------------------------------------------------------------------
!
!  Temporary work arrays
!
!-----------------------------------------------------------------------
!
  REAL, INTENT(INOUT) :: tem_ext(nx,ny,nz_ext)
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k,ka,kl
  REAL :: wlow
  REAL :: topprt,botprt,arpstop,arpsbot
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

!
!-----------------------------------------------------------------------
!
!  Subtract horizontal mean from external fields.
!
!-----------------------------------------------------------------------
!
  IF (intropt == 1) THEN
    DO j=1,ny
      DO i=1,nx
        DO k=1,nz_ext
          DO kl=2,lvlprof-1
            IF(zsnd(kl) > zext(i,j,k)) EXIT
          END DO
          wlow=(zsnd(kl)-zext(i,j,k)) / (zsnd(kl)-zsnd(kl-1))
          tem_ext(i,j,k)=varext(i,j,k)-EXP(                               &
                         ((1.-wlow)*vlnsnd(kl) + wlow*vlnsnd(kl-1)))
        END DO
      END DO
    END DO
  ELSE        ! intropt=2
    tem_ext(:,:,:) = LOG(varext(:,:,:))
  ENDIF
!
!-----------------------------------------------------------------------
!
!  Loop through all ARPS grid points
!
!-----------------------------------------------------------------------
!
  DO j=1,ny
    DO i=1,nx
!
!-----------------------------------------------------------------------
!
!  Interpolate from the mean sounding to get arpsbar
!
!-----------------------------------------------------------------------
!
      DO ka=1,nz
        DO kl=2,lvlprof-1
          IF(zsnd(kl) > za(i,j,ka)) EXIT
        END DO
        wlow=(zsnd(kl)-za(i,j,ka)) / (zsnd(kl)-zsnd(kl-1))
        arpsbar(i,j,ka)=EXP((1.-wlow)*vlnsnd(kl) + wlow*vlnsnd(kl-1))
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
      botprt= tem_ext(i,j,1)
      topprt= tem_ext(i,j,nz_ext)

      DO ka=1,nz
        IF(za(i,j,ka) < zext(i,j,1)) THEN
          arpsprt(i,j,ka)=botprt
        ELSE IF(za(i,j,ka) > zext(i,j,nz_ext)) THEN
          arpsprt(i,j,ka)=topprt
        ELSE
          DO kl=2,nz_ext-1
            IF(zext(i,j,kl) > za(i,j,ka)) EXIT
          END DO
          wlow=(zext(i,j,kl)-  za(i,j,ka))/                        &
               (zext(i,j,kl)-zext(i,j,kl-1))

          arpstop= tem_ext(i,j,kl)

          arpsbot= tem_ext(i,j,kl-1)

          arpsprt(i,j,ka)=(1.-wlow)*arpstop+wlow*arpsbot
        END IF
      END DO
    END DO
  END DO

  IF(intropt == 2) THEN
    arpsprt(:,:,:) = EXP(arpsprt(:,:,:))
  ENDIF
!
!-----------------------------------------------------------------------
!
!  If total quantity, rather than perturbation quantity desired,
!  add the level mean to the interpolated perturbation variable
!  at each grid point.
!
!-----------------------------------------------------------------------
!
  IF(iprtopt == 0 .and. intropt == 1) THEN
    arpsprt(:,:,:) = arpsprt(:,:,:) + arpsbar(:,:,:)
  ELSE IF(iprtopt == 1 .and. intropt == 2) THEN
    arpsprt(:,:,:) = arpsprt(:,:,:) - arpsbar(:,:,:)
  END IF
!
  RETURN
END SUBROUTINE vintrpvlz_wrf
!
!##################################################################
!##################################################################
!######                                                      ######
!######               SUBROUTINE VINTRPSOIL_WRF              ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

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
  REAL,    INTENT(IN)    :: zpsoilext(nx,ny,nzsoil_ext)
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
  END IF
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
        IF(zpsoil(i,j,k) < zpsoilext(i,j,1)) THEN            ! extrapolation
          dist = zpsoilext(i,j,1) - zpsoil(i,j,k)
          w1 = (tsoil_in(i,j,1) - tsoil_in(i,j,2)) /                    &
               (zpsoilext(i,j,2)- zpsoilext(i,j,1))
          w2 = (qsoil_in(i,j,1) - qsoil_in(i,j,2)) /                    &
               (zpsoilext(i,j,2)- zpsoilext(i,j,1))
          tsoil_out(i,j,k) = tsoil_in(i,j,1) + w1*dist
          qsoil_out(i,j,k) = qsoil_in(i,j,1) + w2*dist
        ELSE IF(zpsoil(i,j,k) >= zpsoilext(i,j,nzsoil_ext)) THEN ! constant layers
          tsoil_out(i,j,k) = tsoil_in(i,j,nzsoil_ext)
          qsoil_out(i,j,k) = qsoil_in(i,j,nzsoil_ext)
        ELSE                         ! linear interpolation
          DO m = 1, nzsoil_ext
            IF(zpsoilext(i,j,m) > zpsoil(i,j,k)) EXIT
          END DO
          w1 = (zpsoilext(i,j,m) - zpsoil(i,j,k)) /                     &
               (zpsoilext(i,j,m)- zpsoilext(i,j,m-1))
          w2 = (zpsoil(i,j,k) - zpsoilext(i,j,m-1) ) /                  &
               (zpsoilext(i,j,m)- zpsoilext(i,j,m-1))
          tsoil_out(i,j,k) = tsoil_in(i,j,m-1)*w1 + tsoil_in(i,j,m)*w2
          qsoil_out(i,j,k) = qsoil_in(i,j,m-1)*w1 + qsoil_in(i,j,m)*w2
        END IF

      END DO
    END DO
  END DO

  RETURN
END SUBROUTINE vintrpsoil_wrf
