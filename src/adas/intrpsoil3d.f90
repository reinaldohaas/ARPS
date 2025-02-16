!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINE INTRP_SOIL                ######
!######                                                      ######
!######                    Developed by                      ######
!######            Weather Decision Technologies             ######
!######                                                      ######
!##################################################################
!##################################################################

SUBROUTINE intrp_soil(nx,ny,nx1,ny1,nstyp,nstyp1,wx,wy,ix,jy,     &
                      tsfc,tsoil,wetsfc,wetdp,wetcanp,            &
                      soiltyp,stypfrct,vegtyp,                    &
                      tsfc1,tsoil1,wetsfc1,wetdp1,wetcanp1,       &
                      soiltyp1,stypfrct1,vegtyp1)
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Interpolate soil properties onto another grid.
!
!  NOTE:  This should be used with the old ARPS Force-Restore Soil
!         Model.  When using the new OUSoil model, use intrpsoil3d_avg
!         or intrpsoil3d_pst.  (Eric Kemp, 18 June 2002).
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Gene Bassett
!  2000/11/01
!
!  05/20/2011 Y. Wang
!  Fixed an issue that assigns zero to soiltyp1. Zero value for soil
!  type will cause a problem with the ARPS forecast model.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!  nx
!  ny
!  nx1
!  ny1
!  xw           Weight factor in x-direction
!  yw           Weight factor in y-direction
!  ix           Old grid index (lower left) for new grid point
!  jy           Old grid index (lower left) for new grid point
!  nstyp        Number of soil types on original grid
!  nstyp1       Number of soil types to use on new grid
!  tsfc         Temperature at surface in data set (K)
!  tsoil        Deep soil temperature in data set (K)
!  wetsfc       Surface soil moisture in data set
!  wetdp        Deep soil moisture in data set
!  wetcanp      Canopy water amount in data set
!  soiltyp      Soil type in data set
!  stypfrct     Soil type fraction
!  vegtyp       Vegetation type in data set
!
!  OUTPUT:
!
!  tsfc1        Temperature at surface in data set (K)
!  tsoil1       Deep soil temperature in data set (K)
!  wetsfc1      Surface soil moisture in data set
!  wetdp1       Deep soil moisture in data set
!  wetcanp1     Canopy water amount in data set
!  soiltyp1     Soil type in data set
!  stypfrct1    Soil type fraction
!  vegtyp1      Vegetation type in data set
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER :: nx,ny,nx1,ny1
  INTEGER :: nstyp,nstyp1

  REAL :: wx(nx1,ny1)    ! Weight factor in x-direction
  REAL :: wy(nx1,ny1)    ! Weight factor in y-direction
  INTEGER :: ix(nx1,ny1) ! Old grid index (lower left) for new grid point
  INTEGER :: jy(nx1,ny1) ! Old grid index (lower left) for new grid point

  REAL :: tsfc   (nx,ny,0:nstyp)   ! Temperature at surface (K)
  REAL :: tsoil  (nx,ny,0:nstyp)   ! Deep soil temperature (K)
  REAL :: wetsfc (nx,ny,0:nstyp)   ! Surface soil moisture
  REAL :: wetdp  (nx,ny,0:nstyp)   ! Deep soil moisture
  REAL :: wetcanp(nx,ny,0:nstyp)   ! Canopy water amount
  INTEGER :: soiltyp(nx,ny,nstyp)  ! Soil type in model domain
  REAL :: stypfrct(nx,ny,nstyp)
  INTEGER :: vegtyp(nx,ny)

  REAL    :: tsfc1   (nx1,ny1,0:nstyp1) ! Temperature at surface (K)
  REAL    :: tsoil1  (nx1,ny1,0:nstyp1) ! Deep soil temperature (K)
  REAL    :: wetsfc1 (nx1,ny1,0:nstyp1) ! Surface soil moisture
  REAL    :: wetdp1  (nx1,ny1,0:nstyp1) ! Deep soil moisture
  REAL    :: wetcanp1(nx1,ny1,0:nstyp1) ! Canopy water amount
  INTEGER :: soiltyp1(nx1,ny1,nstyp1)   ! Soil type in model domain
  REAL    :: stypfrct1(nx1,ny1,nstyp1)
  INTEGER :: vegtyp1(nx1,ny1)

  LOGICAL :: warned=.FALSE.

!-----------------------------------------------------------------------
!
!  Include file:
!
!-----------------------------------------------------------------------

  INCLUDE 'arpssfc.inc'

!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------

  ! Arrays start at zero in case no soil type defined (i.e. soiltyp=0)
  REAL :: tsfc1sum    (0:nsoiltyp)     ! Ground sfc. temperature (K)
  REAL :: tsoil1sum   (0:nsoiltyp)     ! Deep soil temperature (K)
  REAL :: wetsfc1sum  (0:nsoiltyp)     ! Surface soil moisture
  REAL :: wetdp1sum   (0:nsoiltyp)     ! Deep soil moisture
  REAL :: wetcanp1sum (0:nsoiltyp)     ! Canopy water amount
  REAL :: stypfrct1sum(0:nsoiltyp)     ! Frction of soil type

  INTEGER :: i,j,i1,j1,is,ii
  REAL    :: weight,maxweight,frctot
  REAL    :: totweight(0:nsoiltyp)
  INTEGER :: maxtype, im
  INTEGER :: soiltype

  REAL    :: inverse

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  tsfc1  (:,:,:) = 0.
  tsoil1 (:,:,:) = 0.
  wetsfc1 (:,:,:) = 0.
  wetdp1  (:,:,:) = 0.
  wetcanp1(:,:,:) = 0.
  soiltyp1 (:,:,:) = 0
  stypfrct1(:,:,:) = 0.
  vegtyp1(:,:) = 0

  IF (nstyp == 1) THEN
    stypfrct = 1.   ! make sure stypfrct is set for nstyp=1
    ! copy level 0 to level 1 if level 1 is undefined:
    IF (tsfc(1,1,1)    <= 0) tsfc   (:,:,1) = tsfc   (:,:,0)
    IF (tsoil(1,1,1)   <= 0) tsoil  (:,:,1) = tsoil  (:,:,0)
    IF (wetsfc(1,1,1)  <  0) wetsfc (:,:,1) = wetsfc (:,:,0)
    IF (wetdp(1,1,1)   <  0) wetdp  (:,:,1) = wetdp  (:,:,0)
    IF (wetcanp(1,1,1) <  0) wetcanp(:,:,1) = wetcanp(:,:,0)
  ENDIF

  DO j1 = 1,ny1-1                 ! desired grid
    DO i1 = 1,nx1-1
      tsfc1sum   = 0.
      tsoil1sum  = 0.
      wetsfc1sum = 0.
      wetdp1sum  = 0.
      wetcanp1sum  = 0.
      stypfrct1sum = 0.
      !vegtyp1sum = 0.
      maxweight = 0.
      totweight = 0.
      vegtyp1(i1,j1) = 0
      DO j = jy(i1,j1), jy(i1,j1)+1       ! external input grid
        DO i = ix(i1,j1), ix(i1,j1)+1
          weight = wx(i1,j1)*(1.+ix(i1,j1)-i)*wy(i1,j1)*(1+jy(i1,j1)-j)     &
                 + (1.-wx(i1,j1))*(i-ix(i1,j1))*wy(i1,j1)*(1+jy(i1,j1)-j)   &
                 + wx(i1,j1)*(1.+ix(i1,j1)-i)*(1.-wy(i1,j1))*(j-jy(i1,j1))  &
                 + (1.-wx(i1,j1))*(i-ix(i1,j1))*(1.-wy(i1,j1))*(j-jy(i1,j1))
          DO is=1,nstyp
            IF (stypfrct(i,j,is) > 0) THEN
              !soiltype = soiltyp(i,j,is)         ! What is soiltype used for? Comment it out by Y.W.
              !IF (soiltype < 1 .OR. soiltype > nsoiltyp) soiltype = 0
              tsfc1sum(soiltyp(i,j,is)) = tsfc1sum(soiltyp(i,j,is))     &
                                + weight*stypfrct(i,j,is)*tsfc(i,j,is)
              tsoil1sum(soiltyp(i,j,is)) = tsoil1sum(soiltyp(i,j,is))   &
                                + weight*stypfrct(i,j,is)*tsoil(i,j,is)
              wetsfc1sum(soiltyp(i,j,is)) = wetsfc1sum(soiltyp(i,j,is)) &
                                + weight*stypfrct(i,j,is)*wetsfc(i,j,is)
              wetdp1sum(soiltyp(i,j,is)) = wetdp1sum(soiltyp(i,j,is))   &
                                + weight*stypfrct(i,j,is)*wetdp(i,j,is)
              wetcanp1sum(soiltyp(i,j,is)) = wetcanp1sum(soiltyp(i,j,is)) &
                                + weight*stypfrct(i,j,is)*wetcanp(i,j,is)
              stypfrct1sum(soiltyp(i,j,is)) = stypfrct1sum(soiltyp(i,j,is)) &
                                + weight*stypfrct(i,j,is)
              totweight(soiltyp(i,j,is)) = totweight(soiltyp(i,j,is)) + weight
            END IF
          END DO
          ! use the vegtyp of the closest point
          IF (weight > maxweight) THEN
            vegtyp1(i1,j1) = vegtyp(i,j)
            maxweight = weight
          END IF
          !!use most popular vegtyp weighted by distance
          !IF (vegtyp(i,j) > 0 .and. vegtyp(i,j) <= nvegtyp) THEN
          !  vegtyp1sum(vegtyp(i,j)) =  vegtyp1sum(vegtyp(i,j)) + weight
          !END IF
        END DO
      END DO

      !!use most popular vegtyp weighted by distance
      !maxweight = 0
      !vegtyp1(i1,j1) = 0
      !DO ii = 1,nvegtyp
      !  IF (vegtyp1sum(ii) > maxweight) THEN
      !    vegtyp1(i1,j1) = ii
      !    maxweight = vegtyp1sum(ii)
      !  END IF
      !END DO

      ! get the nstyp1 largest weighted soil types

      DO is = 1,nstyp1
        maxtype   = -1
        maxweight = 0.
        DO ii = 0,nsoiltyp
          IF (stypfrct1sum(ii) > maxweight) THEN
            maxweight = stypfrct1sum(ii)
            maxtype   = ii
          END IF
        END DO
        IF (maxtype /= -1) THEN
          soiltyp1(i1,j1,is) = maxtype
          inverse = 1./stypfrct1sum(maxtype)

          tsfc1(i1,j1,is)    = tsfc1sum(maxtype)    * inverse
          tsoil1(i1,j1,is)   = tsoil1sum(maxtype)   * inverse
          wetsfc1(i1,j1,is)  = wetsfc1sum(maxtype)  * inverse
          wetdp1(i1,j1,is)   = wetdp1sum(maxtype)   * inverse
          wetcanp1(i1,j1,is) = wetcanp1sum(maxtype) * inverse

          IF (nstyp == 1) THEN
            stypfrct1(i1,j1,is)= stypfrct1sum(maxtype)
          ELSE
            stypfrct1(i1,j1,is)= stypfrct1sum(maxtype)/totweight(maxtype)
          END IF
          stypfrct1sum(maxtype) = 0.  ! maxtype should never be reused for this grid point (i1,j1)

          !stypfrct1(i1,j1,is) = stypfrct1sum(maxtype)/totweight(maxtype)
        ELSE
          IF (is == 1) THEN
            IF (.NOT. warned) THEN
              WRITE (6,*) 'INTRP_SOIL: WARNING, no soil type found, ',  &
                          'variables not assigned!'
              warned = .TRUE.
            END IF
          ELSE
            stypfrct1(i1,j1,is) = 0.

            im = MIN(is,nstyp)  ! this was added to avoid zeros
            soiltyp1(i1,j1,is) = soiltyp(ix(i1,j1),jy(i1,j1),im)

            tsfc1(i1,j1,is)    = tsfc(ix(i1,j1),jy(i1,j1),im)
            tsoil1(i1,j1,is)   = tsoil(ix(i1,j1),jy(i1,j1),im)
            wetsfc1(i1,j1,is)  = wetsfc(ix(i1,j1),jy(i1,j1),im)
            wetdp1(i1,j1,is)   = wetdp(ix(i1,j1),jy(i1,j1),im)
            wetcanp1(i1,j1,is) = wetcanp(ix(i1,j1),jy(i1,j1),im)

          ENDIF
        ENDIF
      END DO

      ! renormalize

      frctot = 0.
      DO is = 1,nstyp1
        frctot = frctot + stypfrct1(i1,j1,is)
      END DO
      IF (frctot /= 0) THEN
        inverse = 1.0 / frctot
        DO is = 1,nstyp1
          stypfrct1(i1,j1,is) = stypfrct1(i1,j1,is) * inverse
        END DO
      ELSE
        stypfrct1(i1,j1,1) = 1.
        IF (nstyp1 .gt. 1) stypfrct1(i1,j1,2:nstyp1) = 0.
      END IF
      tsfc1(i1,j1,0)    = 0.
      tsoil1(i1,j1,0)   = 0.
      wetsfc1(i1,j1,0)  = 0.
      wetdp1(i1,j1,0)   = 0.
      wetcanp1(i1,j1,0) = 0.
      DO is = 1,nstyp1
        tsfc1(i1,j1,0)    = tsfc1(i1,j1,0)                              &
                          + stypfrct1(i1,j1,is)*tsfc1(i1,j1,is)
        tsoil1(i1,j1,0)   = tsoil1(i1,j1,0)                             &
                          + stypfrct1(i1,j1,is)*tsoil1(i1,j1,is)
        wetsfc1(i1,j1,0)  = wetsfc1(i1,j1,0)                            &
                          + stypfrct1(i1,j1,is)*wetsfc1(i1,j1,is)
        wetdp1(i1,j1,0)   = wetdp1(i1,j1,0)                             &
                          + stypfrct1(i1,j1,is)*wetdp1(i1,j1,is)
        wetcanp1(i1,j1,0) = wetcanp1(i1,j1,0)                           &
                          + stypfrct1(i1,j1,is)*wetcanp1(i1,j1,is)
      END DO

    END DO
  END DO

END SUBROUTINE intrp_soil

!##################################################################
!##################################################################
!######                                                      ######
!######              SUBROUTINE intrpsoil3d_avg              ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################

SUBROUTINE intrpsoil3d_avg(nx_ext,ny_ext,nzsoil_ext,                    &
                           vegtyp_ext,                                  &
                           tsoil_ext,qsoil_ext,wetcanp_ext,             &
                           x_ext,y_ext,soildepth_ext,                   &
                           rdxfld_ext,rdyfld_ext,rdzsoilfld_ext,        &
                           nx,ny,nzsoil,nstyp,                          &
                           vegtyp,                                      &
                           tsoil,qsoil,wetcanp,                         &
                           x2d,y2d,soildepth,                           &
                           i2d,j2d,k3d)

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! Tri-linearly/bi-linearly interpolates average tsoil,qsoil,wetcanp
! from external grid, and assigns the interpolated value to all specific
! tsoil,qsoil,wetcanp soil type values on new grid.  Nearest neighbor
! method is used for remapping vegtyp.  This subroutine is recommended
! for use with Eta or RUC data where external soil types are unknown.
!
! NOTE:  This should not be used with the old ARPS Force-Restore Soil
!        Model.  When using that model, use intrp_soil.
!
!-----------------------------------------------------------------------
!
! HISTORY:
!
! First written 6 June 2002 (Eric Kemp)
!
! Renamed zpsoil arrays to soildepth.  7 June 2002 (Eric Kemp)
!
! Bug fix for calculating distances for vegtyp remapping.  20 June 2002
! (Eric Kemp)
!
!-----------------------------------------------------------------------
!
! Force explicit declarations
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------
!
! Variables and arrays for external grid
!
!-----------------------------------------------------------------------

  INTEGER :: nx_ext,ny_ext             ! Grid dimensions
  INTEGER :: nzsoil_ext                ! Number of soil levels

  REAL :: vegtyp_ext(nx_ext,ny_ext)    ! Vegetation type
  REAL :: soildepth_ext(nx_ext,ny_ext,nzsoil_ext)
                                       ! Depth of soil level (m)
  REAL :: tsoil_ext(nx_ext,ny_ext,nzsoil_ext)
                                       ! Average soil temperature (K)
  REAL :: qsoil_ext(nx_ext,ny_ext,nzsoil_ext)
                                       ! Average soil moisture (m**3/m**3)
  REAL :: wetcanp_ext(nx_ext,ny_ext)
                                       ! Average canopy water amount
  REAL :: x_ext(nx_ext) ! x-coord of external points
  REAL :: y_ext(ny_ext) ! y-coord of external points
  REAL :: rdxfld_ext(nx_ext) ! Reciprocal of local dx for external grid.
  REAL :: rdyfld_ext(ny_ext) ! Reciprocal of local dy for external grid.
  REAL :: rdzsoilfld_ext(nx_ext,ny_ext,nzsoil_ext) ! Reciprocal of local
                                                   ! dzsoil for external
                                                   ! grid.

!-----------------------------------------------------------------------
!
! Variables and arrays for grid to interpolate to.
!
!-----------------------------------------------------------------------

  INTEGER :: nx,ny                     ! Grid dimensions
  INTEGER :: nzsoil                    ! Number of soil levels
  INTEGER :: nstyp                     ! Number of soil types per grid box
  INTEGER :: soiltyp(nx,ny,nstyp)      ! Soil type in model domain
  REAL :: stypfrct(nx,ny,nstyp)        ! Soil fraction
  REAL :: vegtyp(nx,ny)                ! Vegetation type
  REAL :: soildepth(nx,ny,nzsoil)      ! Depth of soil level (m)
  REAL :: tsoil(nx,ny,nzsoil,0:nstyp)  ! Soil temperature (K)
  REAL :: qsoil(nx,ny,nzsoil,0:nstyp)  ! Soil moisture (m**3/m**3)
  REAL :: wetcanp(nx,ny,0:nstyp)       ! Canopy water amount
  REAL :: x2d(nx,ny)        ! x-coord of interpolation point w.r.t.
                            ! external grid.
  REAL :: y2d(nx,ny)        ! y-coord of interpolation point w.r.t.
                            ! external grid.
  INTEGER :: i2d(nx,ny)        ! i-index of interpolation point w.r.t.
                               ! external grid.
  INTEGER :: j2d(nx,ny) ! j-index of interpolation point w.r.t.
                               ! external grid.
  INTEGER :: k3d(nx,ny,nzsoil) ! k-index of interpolation point w.r.t.
                               ! external grid.

!-----------------------------------------------------------------------
!
! Internal variables
!
!-----------------------------------------------------------------------

  INTEGER :: i,j,k,is
  INTEGER :: i_ext,j_ext
  REAL    :: dist11,dist12,dist21,dist22,mindist
  REAL    :: d1x, d1y, d2x, d2y
  INTEGER :: ibeg,iend
  INTEGER :: jbeg,jend
  INTEGER :: kbeg,kend

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


!-----------------------------------------------------------------------
!
! Get i,j,k external grid anchor points.
!
!-----------------------------------------------------------------------

  CALL setijkloc3d(nx_ext,ny_ext,nzsoil_ext,x_ext,y_ext,soildepth_ext, &
                   nx,ny,nzsoil,x2d,y2d,soildepth,i2d,j2d,k3d)

!-----------------------------------------------------------------------
!
! Copy vegtyp of the closest point
!
!-----------------------------------------------------------------------

  DO j = 1,ny-1
    DO i = 1,nx-1
      i_ext = i2d(i,j)
      j_ext = j2d(i,j)

      d1x = x2d(i,j) - x_ext(i_ext)
      d1x = d1x * d1x

      d1y = y2d(i,j) - y_ext(j_ext)
      d1y = d1y * d1y

      d2x = x2d(i,j) - x_ext(i_ext+1)
      d2x = d2x * d2x

      d2y = y2d(i,j) - y_ext(j_ext+1)
      d2y = d2y * d2y

      dist11 = SQRT( d1x + d1y)
      dist12 = SQRT( d1x + d2y)
      dist21 = SQRT( d2x + d1y)
      dist22 = SQRT( d2x + d2y)

      mindist = dist11
      vegtyp(i,j) = vegtyp_ext(i_ext,j_ext)
      IF (dist12 < mindist) THEN
        mindist = dist12
        vegtyp(i,j) = vegtyp_ext(i_ext,j_ext+1)
      END IF
      IF (dist21 < mindist) THEN
        mindist = dist21
        vegtyp(i,j) = vegtyp_ext(i_ext+1,j_ext)
      END IF
      IF (dist22 < mindist) THEN
        mindist = dist22
        vegtyp(i,j) = vegtyp_ext(i_ext+1,j_ext+1)
      END IF

    END DO       ! DO i = 1,nx-1
  END DO       ! DO j = 1,ny-1

!-----------------------------------------------------------------------
!
! Interpolate average tsoil, qsoil, and wetcanp to new grid.
!
!-----------------------------------------------------------------------

  ibeg = 1
  iend = nx-1
  jbeg = 1
  jend = ny-1
  kbeg = 1
  kend = nzsoil

  CALL tri_linear_intrp(nx_ext,ny_ext,nzsoil_ext,x_ext,y_ext,soildepth_ext, &
                        rdxfld_ext,rdyfld_ext,rdzsoilfld_ext,           &
                        tsoil_ext,                                      &
                        nx,ny,nzsoil, ibeg,iend, jbeg,jend, kbeg,kend,  &
                        x2d,y2d,soildepth, i2d,j2d,k3d,                 &
                        tsoil(1,1,1,0) )

  CALL tri_linear_intrp(nx_ext,ny_ext,nzsoil_ext,x_ext,y_ext,soildepth_ext, &
                        rdxfld_ext,rdyfld_ext,rdzsoilfld_ext,           &
                        qsoil_ext,                                      &
                        nx,ny,nzsoil,ibeg,iend,jbeg,jend,kbeg,kend,     &
                        x2d,y2d,soildepth, i2d,j2d,k3d,                 &
                        qsoil(1,1,1,0) )

  CALL bi_linear_intrp(nx_ext,ny_ext, x_ext,y_ext,                      &
                       rdxfld_ext,rdyfld_ext,                           &
                       wetcanp_ext,                                     &
                       nx,ny, ibeg,iend, jbeg,jend,                     &
                       x2d,y2d, i2d,j2d,                                &
                       wetcanp(1,1,0) )

!-----------------------------------------------------------------------
!
! Now copy the interpolated averages of tsoil, qsoil, and wetcanp to the
! values for each soil type.
!
!-----------------------------------------------------------------------

  DO is = 1,nstyp
    tsoil(:,:,:,is) = tsoil(:,:,:,0)
    qsoil(:,:,:,is) = qsoil(:,:,:,0)
    wetcanp(:,:,is) = wetcanp(:,:,0)
  END DO ! DO is = 1,nstyp

  RETURN
END SUBROUTINE intrpsoil3d_avg

!##################################################################
!##################################################################
!######                                                      ######
!######              SUBROUTINE intrpsoil3d_pst              ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################

SUBROUTINE intrpsoil3d_pst(nx_ext,ny_ext,nzsoil_ext, nstyp_ext, &
                           soiltyp_ext,stypfrct_ext,vegtyp_ext, &
                           tsoil_ext,qsoil_ext,wetcanp_ext, &
                           x_ext,y_ext,soildepth_ext, &
                           rdxfld_ext,rdyfld_ext,rdzsoilfld_ext, &
                           nx,ny,nzsoil,nstyp, &
                           soiltyp,stypfrct,vegtyp, &
                           tsoil,qsoil,wetcanp, &
                           x2d,y2d,soildepth, &
                           i2d,j2d,k3d)

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! Remaps tsoil, qsoil, wetcanp, vegtyp, and soil types from external
! grid, preserving the external grid soil types to new grid.  Designed
! for use in interpolated outer/coarse ARPS grid to inner/fine ARPS
! grid without using ARPSSFC for the inner grid.
!
!-----------------------------------------------------------------------
!
! HISTORY:
!
! First written June 2002.  Based on intrp_soil subroutine (Eric Kemp).
!
! 13 June 2002 (Eric Kemp)
! Added check for new soil levels outside of original grid.  Also
! added minor bug fixes with help from Yunheng Wang and Dan Weber.
!
! NOTE:  This should not be used with the old ARPS Force-Restore Soil
!        Model.  When using that model, use intrp_soil.
!
!-----------------------------------------------------------------------
!
! Force explicit declarations
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------
!
! Include files
!
!-----------------------------------------------------------------------

  INCLUDE 'arpssfc.inc'

!-----------------------------------------------------------------------
!
! Variables and arrays for external grid
!
!-----------------------------------------------------------------------

  INTEGER :: nx_ext,ny_ext             ! Grid dimensions
  INTEGER :: nzsoil_ext                ! Number of soil levels
  INTEGER :: nstyp_ext                 ! Number of soil types per grid box
  INTEGER :: soiltyp_ext(nx_ext,ny_ext,nstyp_ext)
                                       ! Soil type in model domain
  REAL :: stypfrct_ext(nx_ext,ny_ext,nstyp_ext)
                                       ! Soil fraction
  REAL :: vegtyp_ext(nx_ext,ny_ext)    ! Vegetation type
  REAL :: soildepth_ext(nx_ext,ny_ext,nzsoil_ext)
                                       ! Depth of soil level (m)
  REAL :: tsoil_ext(nx_ext,ny_ext,nzsoil_ext,0:nstyp_ext)
                                       ! Soil temperature (K)
  REAL :: qsoil_ext(nx_ext,ny_ext,nzsoil_ext,0:nstyp_ext)
                                       ! Soil moisture (m**3/m**3)
  REAL :: wetcanp_ext(nx_ext,ny_ext,0:nstyp_ext)
                                       ! Canopy water amount
  REAL :: x_ext(nx_ext) ! x-coord of external points
  REAL :: y_ext(ny_ext) ! y-coord of external points
  REAL :: rdxfld_ext(nx_ext) ! Reciprocal of local dx for external grid.
  REAL :: rdyfld_ext(ny_ext) ! Reciprocal of local dy for external grid.
  REAL :: rdzsoilfld_ext(nx_ext,ny_ext,nzsoil_ext) ! Reciprocal of local
                                                   ! dz for external grid.

!-----------------------------------------------------------------------
!
! Variables and arrays for grid to interpolate to.
!
!-----------------------------------------------------------------------

  INTEGER :: nx,ny                     ! Grid dimensions
  INTEGER :: nzsoil                    ! Number of soil levels
  INTEGER :: nstyp                     ! Number of soil types per grid box
  INTEGER :: soiltyp(nx,ny,nstyp)      ! Soil type in model domain
  REAL :: stypfrct(nx,ny,nstyp)        ! Soil fraction
  REAL :: vegtyp(nx,ny)                ! Vegetation type
  REAL :: soildepth(nx,ny,nzsoil)      ! Depth of soil level (m)
  REAL :: tsoil(nx,ny,nzsoil,0:nstyp)  ! Soil temperature (K)
  REAL :: qsoil(nx,ny,nzsoil,0:nstyp)  ! Soil moisture (m**3/m**3)
  REAL :: wetcanp(nx,ny,0:nstyp)       ! Canopy water amount
  REAL :: x2d(nx,ny)        ! x-coord of interpolation point w.r.t.
                            ! external grid.
  REAL :: y2d(nx,ny)        ! y-coord of interpolation point w.r.t.
                            ! external grid.
  INTEGER :: i2d(nx,ny)        ! i-index of interpolation point w.r.t.
                               ! external grid.
  INTEGER :: j2d(nx,ny) ! j-index of interpolation point w.r.t.
                               ! external grid.
  INTEGER :: k3d(nx,ny,nzsoil) ! k-index of interpolation point w.r.t.
                               ! external grid.

!-----------------------------------------------------------------------
!
! Internal variables
!
!-----------------------------------------------------------------------

  REAL :: tsoilsum(nsoiltyp)
  REAL :: qsoilsum(nsoiltyp)
  REAL :: wetcanpsum(nsoiltyp)
  REAL :: stypfrctsum(nsoiltyp)
  REAL :: totweight(nsoiltyp)

  INTEGER :: soiltyp_ext11,soiltyp_ext21,soiltyp_ext12,soiltyp_ext22
  REAL :: wetcanp_ext11,wetcanp_ext21,wetcanp_ext12,wetcanp_ext22
  REAL :: stypfrct_ext11,stypfrct_ext21,stypfrct_ext12,stypfrct_ext22
  REAL :: tsoil_ext111,tsoil_ext112,tsoil_ext121,tsoil_ext122, &
          tsoil_ext211,tsoil_ext212,tsoil_ext221,tsoil_ext222
  REAL :: qsoil_ext111,qsoil_ext112,qsoil_ext121,qsoil_ext122, &
          qsoil_ext211,qsoil_ext212,qsoil_ext221,qsoil_ext222
  REAL :: c1,c2,c3,c4,c5,c6
  REAL :: c11,c12,c21,c22
  REAL :: c111,c112,c121,c122,c211,c212,c221,c222
  REAL :: temx,temy,temz
  REAL :: temxext,temyext,temzext
!  REAL :: dx_ext,dy_ext,dz_ext
!  REAL :: dxinv_ext,dyinv_ext
!  REAL :: dzinv_ext(nzsoil_ext)
  INTEGER :: i,j,k,is,ii
  INTEGER :: i_ext,j_ext,k_ext
  REAL :: frctot
  INTEGER :: maxtype
  REAL :: maxweight

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!-----------------------------------------------------------------------
!
! Get i,j,k external grid anchor points.
!
!-----------------------------------------------------------------------

  CALL setijkloc3d(nx_ext,ny_ext,nzsoil_ext,x_ext,y_ext,soildepth_ext, &
                   nx,ny,nzsoil,x2d,y2d,soildepth,i2d,j2d,k3d)

!-----------------------------------------------------------------------
!
! Define external grid dx, dy, dz  (dz varies locally).
!
!-----------------------------------------------------------------------
!
!  dxinv_ext = 1./(x_ext(2) - x_ext(1))
!  dyinv_ext = 1./(y_ext(2) - y_ext(1))
!
!  dzinv_ext(:) = 0.
!  DO k = 1,nzsoil_ext-1
!    dzinv_ext(k) = 1./ &
!      (soildepth_ext(1,1,k+1)-soildepth_ext(1,1,k))
!  END DO ! DO k = 1,nzsoil_ext
!  dzinv_ext(nzsoil_ext) = dzinv_ext(nzsoil_ext-1)
!
!-----------------------------------------------------------------------
!
! Loop through each column on new grid.
!
!-----------------------------------------------------------------------

  DO k = 1,nzsoil
    DO j = 1,ny-1
      DO i = 1,nx-1

!-----------------------------------------------------------------------
!
!       Get i,j indices relative to external grid.
!
!-----------------------------------------------------------------------

        i_ext = i2d(i,j)
        i_ext = MIN(MAX(i_ext,1),(nx-1))

        j_ext = j2d(i,j)
        j_ext = MIN(MAX(j_ext,1),(ny-1))

        k_ext = k3d(i,j,k)
!        k_ext = MIN(MAX(k_ext,1),(nzsoil-1))

        IF (k_ext < 1 .OR. k_ext >= nzsoil_ext) THEN
          WRITE(6,*)'INTRPSOIL3D_PST: Soil level outside external grid. i,j,k: ',i,j,k
          IF (k > 1) THEN
            WRITE(6,*)'INTRPSOIL3D_PST: Copying tsoil/qsoil from k-1 level.'
            tsoil(i,j,k,:) = tsoil(i,j,k-1,:)
            qsoil(i,j,k,:) = qsoil(i,j,k-1,:)

!            WRITE(6,*)'EMK: tsoil(i,j,k,:) = ',tsoil(i,j,k,:)
!            WRITE(6,*)'EMK: qsoil(i,j,k,:) = ',qsoil(i,j,k,:)

            CYCLE ! Move to next i,j,k point
          END IF
        END IF

!-----------------------------------------------------------------------
!
!       Get x,y,z coordinates relative to external grid.
!
!-----------------------------------------------------------------------

        temx = x2d(i,j)
        temy = y2d(i,j)
        temz = soildepth(i,j,k)

!-----------------------------------------------------------------------
!
!       Get x,y,z coordinates of external "anchor" grid point.
!
!-----------------------------------------------------------------------

        temxext = x_ext(i_ext)
        temyext = y_ext(j_ext)
        temzext = soildepth_ext(i_ext,j_ext,k_ext)

!-----------------------------------------------------------------------
!
!       Calculate bi-linear and tri-linear interpolation weights.
!
!-----------------------------------------------------------------------

!        c2 = (temx - temxext)*dxinv_ext
        c2 = (temx - temxext)*rdxfld_ext(i_ext)
        c1 = 1.0 - c2

!        c4 = (temy - temyext)*dyinv_ext
        c4 = (temy - temyext)*rdyfld_ext(j_ext)
        c3 = 1.0 - c4

!        c6 = (temz - temzext)*dzinv_ext(k_ext)
        c6 = (temz - temzext)*rdzsoilfld_ext(i_ext,j_ext,k_ext)
        c5 = 1.0 - c6

        c11 = c1*c3
        c21 = c2*c3
        c12 = c1*c4
        c22 = c2*c4

        c111 = c1*c3*c5
        c112 = c1*c3*c6
        c121 = c1*c4*c5
        c122 = c1*c4*c6
        c211 = c2*c3*c5
        c212 = c2*c3*c6
        c221 = c2*c4*c5
        c222 = c2*c4*c6

!-----------------------------------------------------------------------
!
!       Initialize sums.
!
!-----------------------------------------------------------------------

        tsoilsum(:) = 0.
        qsoilsum(:) = 0.
        wetcanpsum(:) = 0.
        stypfrctsum(:) = 0.
        totweight(:) = 0.
        maxweight = 0.
        totweight = 0.

!-----------------------------------------------------------------------
!
!       Loop through each soil type on external grid.
!
!-----------------------------------------------------------------------

        DO is = 1,nstyp_ext

!-----------------------------------------------------------------------
!
!         Extract and save soil type, soil type fraction, soil
!         temperature and soil moisture from surrounding points.
!
!-----------------------------------------------------------------------

          soiltyp_ext11 = soiltyp_ext(i_ext  ,j_ext  ,is)
          soiltyp_ext21 = soiltyp_ext(i_ext+1,j_ext  ,is)
          soiltyp_ext12 = soiltyp_ext(i_ext  ,j_ext+1,is)
          soiltyp_ext22 = soiltyp_ext(i_ext+1,j_ext+1,is)

          stypfrct_ext11 = stypfrct_ext(i_ext  ,j_ext  ,is)
          stypfrct_ext21 = stypfrct_ext(i_ext+1,j_ext  ,is)
          stypfrct_ext12 = stypfrct_ext(i_ext  ,j_ext+1,is)
          stypfrct_ext22 = stypfrct_ext(i_ext+1,j_ext+1,is)

          tsoil_ext111 = tsoil_ext(i_ext  ,j_ext  ,k_ext  ,is)
          tsoil_ext112 = tsoil_ext(i_ext  ,j_ext  ,k_ext+1,is)
          tsoil_ext121 = tsoil_ext(i_ext  ,j_ext+1,k_ext  ,is)
          tsoil_ext122 = tsoil_ext(i_ext  ,j_ext+1,k_ext+1,is)
          tsoil_ext211 = tsoil_ext(i_ext+1,j_ext  ,k_ext  ,is)
          tsoil_ext212 = tsoil_ext(i_ext+1,j_ext  ,k_ext+1,is)
          tsoil_ext221 = tsoil_ext(i_ext+1,j_ext+1,k_ext  ,is)
          tsoil_ext222 = tsoil_ext(i_ext+1,j_ext+1,k_ext+1,is)

          qsoil_ext111 = qsoil_ext(i_ext  ,j_ext  ,k_ext  ,is)
          qsoil_ext112 = qsoil_ext(i_ext  ,j_ext  ,k_ext+1,is)
          qsoil_ext121 = qsoil_ext(i_ext  ,j_ext+1,k_ext  ,is)
          qsoil_ext122 = qsoil_ext(i_ext  ,j_ext+1,k_ext+1,is)
          qsoil_ext211 = qsoil_ext(i_ext+1,j_ext  ,k_ext  ,is)
          qsoil_ext212 = qsoil_ext(i_ext+1,j_ext  ,k_ext+1,is)
          qsoil_ext221 = qsoil_ext(i_ext+1,j_ext+1,k_ext  ,is)
          qsoil_ext222 = qsoil_ext(i_ext+1,j_ext+1,k_ext+1,is)


!-----------------------------------------------------------------------
!
!         If the first k level is being operated on, extract wet canopy.
!
!-----------------------------------------------------------------------

          IF (k == 1) THEN
            wetcanp_ext11 = wetcanp_ext(i_ext  ,j_ext  ,is)
            wetcanp_ext21 = wetcanp_ext(i_ext+1,j_ext  ,is)
            wetcanp_ext12 = wetcanp_ext(i_ext  ,j_ext+1,is)
            wetcanp_ext22 = wetcanp_ext(i_ext+1,j_ext+1,is)

!-----------------------------------------------------------------------
!
!           Calculate a weighted sum of wet canopy for each soil type as
!           a function of distance (using bi-linear interpolation
!           weight) and the soil type fraction.
!
!-----------------------------------------------------------------------

            IF (stypfrct_ext11 > 0.) THEN
              wetcanpsum(soiltyp_ext11) = wetcanpsum(soiltyp_ext11) + &
                c11*stypfrct_ext11*wetcanp_ext11
              totweight(soiltyp_ext11) = totweight(soiltyp_ext11) + c11
            END IF

            IF (stypfrct_ext21 > 0.) THEN
              wetcanpsum(soiltyp_ext21) = wetcanpsum(soiltyp_ext21) + &
                c21*stypfrct_ext21*wetcanp_ext21
              totweight(soiltyp_ext21) = totweight(soiltyp_ext21) + c21
            END IF

            IF (stypfrct_ext12 > 0.) THEN
              wetcanpsum(soiltyp_ext12) = wetcanpsum(soiltyp_ext12) + &
                c12*stypfrct_ext12*wetcanp_ext12
              totweight(soiltyp_ext12) = totweight(soiltyp_ext12) + c12
            END IF

            IF (stypfrct_ext22 > 0.) THEN
              wetcanpsum(soiltyp_ext22) = wetcanpsum(soiltyp_ext22) + &
                c22*stypfrct_ext22*wetcanp_ext22
              totweight(soiltyp_ext22) = totweight(soiltyp_ext22) + c22
            END IF

!-----------------------------------------------------------------------
!
!           Copy nearest neighbor vegtyp from external grid to new grid.
!
!-----------------------------------------------------------------------

            maxweight = c11
            vegtyp(i,j) = vegtyp_ext(i_ext  ,j_ext  )
            IF (c12 > maxweight) THEN
              maxweight = c12
              vegtyp(i,j) = vegtyp_ext(i_ext+1,j_ext )
            END IF
            IF (c21 > maxweight) THEN
              maxweight = c21
              vegtyp(i,j) = vegtyp_ext(i_ext  ,j_ext+1)
            END IF
            IF (c22 > maxweight) THEN
              maxweight = c22
              vegtyp(i,j) = vegtyp_ext(i_ext+1,j_ext+1)
            END IF
          END IF ! IF (k == 1)

!-----------------------------------------------------------------------
!
!         Calculate a weighted sum of soil temperature, soil moisture,
!         and soil type fraction for each soil type as a function of
!         distance (using bi-linear or tri-linear interpolation weight)
!         and the soil fraction.
!
!-----------------------------------------------------------------------

          IF (stypfrct_ext11 > 0.) THEN
            tsoilsum(soiltyp_ext11) = tsoilsum(soiltyp_ext11) + &
                c111*stypfrct_ext11*tsoil_ext111 + &
                c112*stypfrct_ext11*tsoil_ext112
            qsoilsum(soiltyp_ext11) = qsoilsum(soiltyp_ext11) + &
                c111*stypfrct_ext11*qsoil_ext111 + &
                c112*stypfrct_ext11*qsoil_ext112
            stypfrctsum(soiltyp_ext11) = stypfrctsum(soiltyp_ext11) + &
                c11*stypfrct_ext11
          END IF

          IF (stypfrct_ext21 > 0.) THEN
            tsoilsum(soiltyp_ext21) = tsoilsum(soiltyp_ext21) + &
                c211*stypfrct_ext21*tsoil_ext211 + &
                c212*stypfrct_ext21*tsoil_ext212
            qsoilsum(soiltyp_ext21) = qsoilsum(soiltyp_ext21) + &
                c211*stypfrct_ext21*qsoil_ext211 + &
                c212*stypfrct_ext21*qsoil_ext212
            stypfrctsum(soiltyp_ext21) = stypfrctsum(soiltyp_ext21) + &
               c21*stypfrct_ext21
          END IF

          IF (stypfrct_ext12 > 0.) THEN
            tsoilsum(soiltyp_ext12) = tsoilsum(soiltyp_ext12) + &
                c121*stypfrct_ext12*tsoil_ext121 + &
                c122*stypfrct_ext12*tsoil_ext122
            qsoilsum(soiltyp_ext12) = qsoilsum(soiltyp_ext12) + &
                c121*stypfrct_ext12*qsoil_ext121 + &
                c122*stypfrct_ext12*qsoil_ext122
            stypfrctsum(soiltyp_ext12) = stypfrctsum(soiltyp_ext12) + &
                c12*stypfrct_ext12
          END IF

          IF (stypfrct_ext22 > 0.) THEN
            tsoilsum(soiltyp_ext22) = tsoilsum(soiltyp_ext22) + &
                c221*stypfrct_ext22*tsoil_ext221 + &
                c222*stypfrct_ext22*tsoil_ext222
            qsoilsum(soiltyp_ext22) = qsoilsum(soiltyp_ext22) + &
                c221*stypfrct_ext22*qsoil_ext221 + &
                c222*stypfrct_ext22*qsoil_ext222
            stypfrctsum(soiltyp_ext22) = stypfrctsum(soiltyp_ext22) + &
                c22*stypfrct_ext22
          END IF

        END DO ! DO is = 1,nstyp_ext

!-----------------------------------------------------------------------
!
!       Loop through all soil types for new grid.
!
!-----------------------------------------------------------------------

        DO is = 1,nstyp
          maxtype = -1
          maxweight = 0.

!-----------------------------------------------------------------------
!
!         Find largest weighted soil type.  (Later on, the current
!         "largest" value is reset to zero, so we don't pick the same
!         type twice.)
!
!-----------------------------------------------------------------------

          DO ii = 1,nsoiltyp
            IF (stypfrctsum(ii) > maxweight) THEN
              maxweight = stypfrctsum(ii)
              maxtype = ii
            END IF
          END DO ! DO ii = 1,nsoiltyp

!-----------------------------------------------------------------------
!
!         If a soil type was picked and we are at the first k level,
!         save the soil type and then calculate the (averaged) wet
!         canopy and soil type fraction for that type at the new grid
!         point.  Otherwise, move on.
!
!-----------------------------------------------------------------------

          IF (k == 1) THEN
            IF (maxtype /= -1) THEN
              soiltyp(i,j,is) = maxtype
              wetcanp(i,j,is) = wetcanpsum(maxtype)/stypfrctsum(maxtype)
              stypfrct(i,j,is) = stypfrctsum(maxtype)/totweight(maxtype)
!Commented out, postponed to tsoil/qsoil calculation further down.
!              stypfrctsum(maxtype) = 0.
            ELSE
              IF (is == 1)  THEN
                WRITE(6,*) &
  "INTRPSOIL3D_PST:  WARNING, no soil type found, variables not assigned!"
              ELSE
                soiltyp(i,j,is) = soiltyp(i,j,is-1)
                stypfrct(i,j,is) = 0.
              END IF
            END IF
          END IF ! IF (k == 1) THEN

!-----------------------------------------------------------------------
!
!         If a soil type was picked, calculate the (averaged) soil
!         temperature and soil moisture for that type at the
!         new grid point.  Otherwise, move on.
!
!-----------------------------------------------------------------------

          IF (maxtype /= -1) THEN
            tsoil(i,j,k,is) = tsoilsum(maxtype)/stypfrctsum(maxtype)
            qsoil(i,j,k,is) = qsoilsum(maxtype)/stypfrctsum(maxtype)
            stypfrctsum(maxtype) = 0.
          END IF

        END DO ! DO is = 1,nstyp

!-----------------------------------------------------------------------
!
!       If we are at the first k level, renormalize the soil type
!       fractions at the new grid point to make sure they add up to
!       unity.  Then, calculate the average wet canopy for the new
!       grid point.
!
!-----------------------------------------------------------------------

        IF (k == 1) THEN
          frctot = 0.
          DO is = 1,nstyp
            frctot = frctot + stypfrct(i,j,is)
          END DO ! DO is = 1,nstyp
          IF (frctot /= 0) THEN
            DO is = 1,nstyp
              stypfrct(i,j,is) = stypfrct(i,j,is)/frctot
            END DO ! DO is = 1,nstyp
          ELSE
            stypfrct(i,j,1) = 1.
            IF (nstyp .GT. 1) stypfrct(i,j,2:nstyp) = 0.
          END IF

          wetcanp(i,j,0) = 0.
          DO is = 1,nstyp
            wetcanp(i,j,0) = wetcanp(i,j,0) + &
                           stypfrct(i,j,is)*wetcanp(i,j,is)
          END DO ! DO is = 1,nstyp
        END IF ! IF (k == 1)

!-----------------------------------------------------------------------
!
!       Finally, calculate the average soil temperature and soil
!       moisture for the new grid point.
!
!-----------------------------------------------------------------------

        tsoil(i,j,k,0) = 0.
        qsoil(i,j,k,0) = 0.
        DO is = 1,nstyp
          tsoil(i,j,k,0) = tsoil(i,j,k,0) + &
                           stypfrct(i,j,is)*tsoil(i,j,k,is)
          qsoil(i,j,k,0) = qsoil(i,j,k,0) + &
                           stypfrct(i,j,is)*qsoil(i,j,k,is)
        END DO ! DO is = 1,nstyp

        DO is = 1,nstyp
          IF (stypfrct(i,j,is) == 0.) THEN
            IF (k == 1) wetcanp(i,j,is) = wetcanp(i,j,0)
            tsoil(i,j,k,is) = tsoil(i,j,k,0)
            qsoil(i,j,k,is) = qsoil(i,j,k,0)
          END IF
        END DO

      END DO ! DO i = 1,nx-1
    END DO ! DO j = 1,ny-1
  END DO ! DO k = 1,nzsoil

  RETURN
END SUBROUTINE intrpsoil3d_pst


!##################################################################
!##################################################################
!######                                                      ######
!######              SUBROUTINE tri_linear_intrp             ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################

SUBROUTINE tri_linear_intrp(nx_ext,ny_ext,nz_ext, x_ext,y_ext,z3d_ext,  &
                            rdxfld_ext,rdyfld_ext,rdzfld_ext,           &
                            var_ext,                                    &
                            nx,ny,nz, ibeg,iend,jbeg,jend, kbeg,kend,   &
                            x2d,y2d,z3d, i2d,j2d,k3d,                   &
                            var)

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! Remaps 3D variable from external grid to internal grid using
! tri-linear interpolation.
!
!-----------------------------------------------------------------------
!
! HISTORY:
!
! First written 4 June 2002. Based on research code from Dan Weber.
! (Eric Kemp)
!
!-----------------------------------------------------------------------
!
! Force explicit declarations
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------
!
! Variables and arrays for external grid
!
!-----------------------------------------------------------------------

  INTEGER :: nx_ext,ny_ext,nz_ext       ! Dimensions of external grid
  REAL :: x_ext(nx_ext) ! x-coord of external grid points.
  REAL :: y_ext(ny_ext) ! y-coord of external grid points.
  REAL :: z3d_ext(nx_ext,ny_ext,nz_ext) ! z-coord of external grid points.
  REAL :: var_ext(nx_ext,ny_ext,nz_ext) ! 3d variable on external grid
                                        ! to be interpolated
  REAL :: rdxfld_ext(nx_ext) ! Reciprocal of local dx on external grid
  REAL :: rdyfld_ext(ny_ext) ! Reciprocal of local dy on external grid
  REAL :: rdzfld_ext(nx_ext,ny_ext,nz_ext) ! Reciprocal of local dz on
                                           ! external grid

!-----------------------------------------------------------------------
!
! Variables and arrays for grid to be interpolated to
!
!-----------------------------------------------------------------------

  INTEGER :: nx,ny,nz      ! Dimensions of grid to interpolate to.
  INTEGER :: ibeg,iend     ! Range of i indices on internal grid to loop
                           ! through.
  INTEGER :: jbeg,jend     ! Range of j indices on internal grid to loop
                           ! through.
  INTEGER :: kbeg,kend     ! Range of k indices on internal grid to loop
                           ! through.
  REAL :: x2d(nx,ny)       ! x-coord of interpolation points w.r.t.
                           ! external grid.
  REAL :: y2d(nx,ny)       ! y-coord of interpolation points w.r.t.
                           ! external grid.
  REAL :: z3d(nx,ny,nz)    ! z-coord of interpolation points w.r.t.
                           ! external grid.
  INTEGER :: i2d(nx,ny)    ! i-index of interpolation point w.r.t.
                           ! external grid.
  INTEGER :: j2d(nx,ny)    ! j-index of interpolation point w.r.t.
                           ! external grid.
  INTEGER :: k3d(nx,ny,nz) ! k-index of interpolation point w.r.t.
                           ! external grid.
  REAL :: var(nx,ny,nz)    ! Interpolated 3d variable.

!-----------------------------------------------------------------------
!
! Internal variables
!
!-----------------------------------------------------------------------

  REAL :: c1,c2,c3,c4,c5,c6 ! Interpolation weights
  REAL :: var111,var211,var121,var221,var112,var212,var122,var222
  INTEGER :: i,j,k,i_ext,j_ext,k_ext
  REAL :: temx,temy,temz
  REAL :: temxext,temyext,temzext
!  REAL :: dxinv_ext,dyinv_ext
!  REAL :: dzinv_ext(nz_ext)

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  IF ((ibeg < 1) .OR. (iend > nx) .OR.                                 &
      (jbeg < 1) .OR. (jend > ny) .OR.                                 &
      (kbeg < 1) .OR. (kend > nz)) THEN
    WRITE(6,*)'tri_linear_intrp:  ERROR -- Mismatch.'
    WRITE(6,*)'ibeg = ',ibeg,' Dimension start = ',1
    WRITE(6,*)'iend = ',iend,' Dimension end (nx) = ',nx

    WRITE(6,*)'jbeg = ',jbeg,' Dimension start = ',1
    WRITE(6,*)'jend = ',jend,' Dimension end (ny) = ',ny

    WRITE(6,*)'kbeg = ',kbeg,' Dimension start = ',1
    WRITE(6,*)'kend = ',kend,' Dimension end (nz) = ',nz

    CALL arpsstop('Aborting...',1)
  END IF

!-----------------------------------------------------------------------
!
! Define external grid dx, dy, dz.  (dz varies locally).
!
!-----------------------------------------------------------------------
!
!  dxinv_ext = 1./(x_ext(2) - x_ext(1))
!  dyinv_ext = 1./(y_ext(2) - y_ext(1))
!
!  dzinv_ext(:) = 0.
!  DO k = 1,nz_ext-1
!    dzinv_ext(k) = 1./(z3d_ext(1,1,k+1) - z3d_ext(1,1,k))
!  END DO ! DO k = 1,nz_ext
!
!-----------------------------------------------------------------------
!
! Loop through all points on grid to be interpolated to.
!
!-----------------------------------------------------------------------

  DO k = kbeg,kend
    DO j = jbeg,jend
      DO i = ibeg,iend

!-----------------------------------------------------------------------
!
!       Get x,y,z coordinates of current interpolation point w.r.t.
!       external grid.
!
!-----------------------------------------------------------------------

        temx = x2d(i,j)
        temy = y2d(i,j)
        temz = z3d(i,j,k)

!        WRITE(6,*)'temx,temy,temz: ',temx,temy,temz

!-----------------------------------------------------------------------
!
!       Get i,j,k indices of current interpolation point w.r.t.
!       external grid.
!
!-----------------------------------------------------------------------

        i_ext = i2d(i,j)
        i_ext = MIN(MAX(i_ext,1),(nx_ext-1))

        j_ext = j2d(i,j)
        j_ext = MIN(MAX(j_ext,1),(ny_ext-1))

        k_ext = k3d(i,j,k)
        !k_ext = MIN(MAX(k_ext,1),(nz_ext-1))

        IF (k_ext < 1 .OR. k_ext >= nz_ext) THEN
          !WRITE(6,*)'TRI_LINEAR_INTRP: Level outside external grid. i,j,k: ',i,j,k
          IF (k > 1) THEN
            !WRITE(6,*)'TRI_LINEAR_INTRP: Copying var from k-1 level.'
            var(i,j,k) = var(i,j,k-1)
            CYCLE ! Move to next i,j,k point
          END IF
        END IF

        !WRITE(6,*)'i_ext,j_ext,k_ext: ',i_ext,j_ext,k_ext

!-----------------------------------------------------------------------
!
!       Get x,y,z coordinates of "anchor" external grid point (southwest
!       and below interpolation point.)
!
!-----------------------------------------------------------------------

        temxext = x_ext(i_ext)
        temyext = y_ext(j_ext)
        temzext = z3d_ext(i_ext,j_ext,k_ext)

!-----------------------------------------------------------------------
!
!       Calculate weights for interpolation.
!
!-----------------------------------------------------------------------

        !c2 = (temx - temxext)*dxinv_ext
        c2 = (temx - temxext)*rdxfld_ext(i_ext)
        c1 = 1.0 - c2

        !c4 = (temy - temyext)*dyinv_ext
        c4 = (temy - temyext)*rdyfld_ext(j_ext)
        c3 = 1.0 - c4

        !c6 = (temz - temzext)*dzinv_ext(k_ext)
        c6 = (temz - temzext)*rdzfld_ext(i_ext,j_ext,k_ext)
        c5 = 1.0 - c6
        !print *, '=== ',i,j,k, k_ext,c5,c6
!EMK TEST
        IF (c1 < 0. .OR. c1 > 1. .OR. &
            c2 < 0. .OR. c2 > 1. .OR. &
            c3 < 0. .OR. c3 > 1. .OR. &
            c4 < 0. .OR. c4 > 1. .OR. &
            c5 < 0. .OR. c5 > 1. .OR. &
            c6 < 0. .OR. c6 > 1. ) THEN
          WRITE(6,*)'tri_linear_intrp3d:  ERROR in weight calculation!'
          WRITE(6,*)'c1: ',c1,' c2: ',c2
          WRITE(6,*)'c3: ',c3,' c4: ',c4
          WRITE(6,*)'c5: ',c5,' c6: ',c6

          WRITE(6,'(1x,a,3I4)')  'i,j,k:             ',i,j,k
          WRITE(6,'(1x,a,3I4)')  'i_ext,j_ext,k_ext: ',i_ext,j_ext,k_ext
          !WRITE(6,'(1x,3F12.5)') temz, temzext, rdzfld_ext(i_ext,j_ext,k_ext)
          CALL arpsstop('Aborting...',1)
        END IF
!EMK END TEST

!-----------------------------------------------------------------------
!
!       Copy surrounding eight points
!
!-----------------------------------------------------------------------

        var111 = var_ext(i_ext  ,j_ext  ,k_ext  ) ! "anchor" point
        var211 = var_ext(i_ext+1,j_ext  ,k_ext  )
        var121 = var_ext(i_ext  ,j_ext+1,k_ext  )
        var221 = var_ext(i_ext+1,j_ext+1,k_ext  )
        var112 = var_ext(i_ext  ,j_ext  ,k_ext+1)
        var212 = var_ext(i_ext+1,j_ext  ,k_ext+1)
        var122 = var_ext(i_ext,  j_ext+1,k_ext+1)
        var222 = var_ext(i_ext+1,j_ext+1,k_ext+1)

!-----------------------------------------------------------------------
!
!       Tri-linear interpolation.
!
!-----------------------------------------------------------------------

        var(i,j,k) = ((var111*c1 + var211*c2)*c3 +                    &
                      (var121*c1 + var221*c2)*c4)*c5 +                &
                     ((var112*c1 + var212*c2)*c3 +                    &
                      (var122*c1 + var222*c2)*c4)*c6

!        WRITE(6,*)'i,j,k,var: ',i,j,k,var(i,j,k)
      END DO ! DO i = ibeg,iend
    END DO ! DO j = jbeg,jend
  END DO ! DO k = kbeg,kend

  RETURN
END SUBROUTINE tri_linear_intrp

!##################################################################
!##################################################################
!######                                                      ######
!######              SUBROUTINE bi_linear_intrp              ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################

SUBROUTINE bi_linear_intrp(nx_ext,ny_ext, &
                           x_ext,y_ext, &
                           rdxfld_ext,rdyfld_ext, &
                           var_ext, &
                           nx,ny, &
                           ibeg,iend, &
                           jbeg,jend, &
                           x2d,y2d, &
                           i2d,j2d, &
                           var)

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! Remaps 2D variable from external grid to internal grid using
! bi-linear interpolation.
!
!-----------------------------------------------------------------------
!
! HISTORY:
!
! First written 6 June 2002.  Based on subroutine tri_linear_interp
! (Eric Kemp)
!
!-----------------------------------------------------------------------
!
! Force explicit declarations
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------
!
! Variables and arrays for external grid
!
!-----------------------------------------------------------------------

  INTEGER :: nx_ext,ny_ext       ! Dimensions of external grid
  REAL :: x_ext(nx_ext)          ! x-coord of external grid points.
  REAL :: y_ext(ny_ext)          ! y-coord of external grid points.
  REAL :: var_ext(nx_ext,ny_ext) ! 2d variable on external grid
                                 ! to be interpolated
  REAL :: rdxfld_ext(nx_ext) ! Reciprocal of local dx on external grid
  REAL :: rdyfld_ext(ny_ext) ! Reciprocal of local dy on external grid

!-----------------------------------------------------------------------
!
! Variables and arrays for grid to be interpolated to
!
!-----------------------------------------------------------------------

  INTEGER :: nx,ny         ! Dimensions of grid to interpolate to.
  INTEGER :: ibeg,iend     ! Range of i indices on internal grid to loop
                           ! through.
  INTEGER :: jbeg,jend     ! Range of j indices on internal grid to loop
                           ! through.
  REAL :: x2d(nx,ny)       ! x-coord of interpolation points w.r.t.
                           ! external grid.
  REAL :: y2d(nx,ny)       ! y-coord of interpolation points w.r.t.
                           ! external grid.
  INTEGER :: i2d(nx,ny)    ! i-index of interpolation point w.r.t.
                           ! external grid.
  INTEGER :: j2d(nx,ny)    ! j-index of interpolation point w.r.t.
                           ! external grid.
  REAL :: var(nx,ny)       ! Interpolated 2d variable.

!-----------------------------------------------------------------------
!
! Internal variables
!
!-----------------------------------------------------------------------

  REAL :: c1,c2,c3,c4  ! Interpolation weights
  REAL :: var11,var21,var12,var22
  INTEGER :: i,j,k,i_ext,j_ext
  REAL :: temx,temy
  REAL :: temxext,temyext
!  REAL :: dxinv_ext,dyinv_ext

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  IF ((ibeg < 1) .OR. (iend > nx) .OR.                                 &
      (jbeg < 1) .OR. (jend > ny)) THEN
    WRITE(6,*)'tri_linear_intrp:  ERROR -- Mismatch.'
    WRITE(6,*)'ibeg = ',ibeg,' Dimension start = ',1
    WRITE(6,*)'iend = ',iend,' Dimension end (nx) = ',nx

    WRITE(6,*)'jbeg = ',jbeg,' Dimension start = ',1
    WRITE(6,*)'jend = ',jend,' Dimension end (ny) = ',ny

    CALL arpsstop('Aborting...',1)
  END IF

!-----------------------------------------------------------------------
!
! Define external grid dx, dy.
!
!-----------------------------------------------------------------------
!
!  dxinv_ext = 1./(x_ext(2) - x_ext(1))
!  dyinv_ext = 1./(y_ext(2) - y_ext(1))
!
!-----------------------------------------------------------------------
!
! Loop through all points on grid to be interpolated to.
!
!-----------------------------------------------------------------------

  DO j = jbeg,jend
    DO i = ibeg,iend

!-----------------------------------------------------------------------
!
!     Get x,y coordinates of current interpolation point w.r.t.
!     external grid.
!
!-----------------------------------------------------------------------

      temx = x2d(i,j)
      temy = y2d(i,j)

!-----------------------------------------------------------------------
!
!     Get i,j indices of current interpolation point w.r.t.
!     external grid.
!
!-----------------------------------------------------------------------

      i_ext = i2d(i,j)
      i_ext = MIN(MAX(i_ext,1),(nx_ext-1))

      j_ext = j2d(i,j)
      j_ext = MIN(MAX(j_ext,1),(ny_ext-1))

!-----------------------------------------------------------------------
!
!     Get x,y coordinates of "anchor" external grid point (southwest
!     of interpolation point.)
!
!-----------------------------------------------------------------------

      temxext = x_ext(i_ext)
      temyext = y_ext(j_ext)

!-----------------------------------------------------------------------
!
!     Calculate weights for interpolation.
!
!-----------------------------------------------------------------------

!      c2 = (temx - temxext)*dxinv_ext
      c2 = (temx - temxext)*rdxfld_ext(i_ext)
      c1 = 1.0 - c2

!      c4 = (temy - temyext)*dyinv_ext
      c4 = (temy - temyext)*rdyfld_ext(j_ext)
      c3 = 1.0 - c4

!-----------------------------------------------------------------------
!
!     Copy surrounding four points
!
!-----------------------------------------------------------------------

      var11 = var_ext(i_ext  ,j_ext  ) ! "anchor" point
      var21 = var_ext(i_ext+1,j_ext  )
      var12 = var_ext(i_ext  ,j_ext+1)
      var22 = var_ext(i_ext+1,j_ext+1)

!-----------------------------------------------------------------------
!
!     Bi-linear interpolation.
!
!-----------------------------------------------------------------------

      var(i,j) = ((var11*c1 + var21*c2)*c3 +                     &
                  (var12*c1 + var22*c2)*c4)

    END DO ! DO i = ibeg,iend
  END DO ! DO j = jbeg,jend

  RETURN
END SUBROUTINE bi_linear_intrp


!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE setijkloc3d                ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################

SUBROUTINE setijkloc3d(nx_ext,ny_ext,nz_ext,x_ext,y_ext,z3d_ext, &
                       nx,ny,nz,x2d,y2d,z3d,i2d,j2d,k3d)

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! Identifies i,j,k coordinates of interior grid relative to external
! grid.
!
!-----------------------------------------------------------------------
!
! HISTORY:
!
! First written 5 June 2002 (Eric Kemp and Dan Weber).
!
!-----------------------------------------------------------------------
!
! Force explicit declarations
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------
!
! External grid variables
!
!-----------------------------------------------------------------------

  INTEGER :: nx_ext,ny_ext,nz_ext        ! External grid dimensions
  REAL :: x_ext(nx_ext)  ! x-coord of external grid
  REAL :: y_ext(ny_ext)  ! y-coord of external grid
  REAL :: z3d_ext(nx_ext,ny_ext,nz_ext)  ! z-coord of external grid

!-----------------------------------------------------------------------
!
! Variables for grid to interpolate to.
!
!-----------------------------------------------------------------------

  INTEGER :: nx,ny,nz                    ! Dimensions of grid to
                                         ! interpolate to.
  REAL :: x2d(nx,ny)                     ! x-coord of interpolation points
                                         ! w.r.t. external grid.
  REAL :: y2d(nx,ny)                     ! y-coord of interpolation points
                                         ! w.r.t. external grid.
  REAL :: z3d(nx,ny,nz)                  ! z-coord of interpolation points
                                         ! w.r.t. external grid.
  INTEGER :: i2d(nx,ny)                  ! i-index of interpolation points
                                         ! w.r.t. external grid.
  INTEGER :: j2d(nx,ny)                  ! j-index of interpolation points
                                         ! w.r.t. external grid.
  INTEGER :: k3d(nx,ny,nz)               ! k-index of interpolation points
                                         ! w.r.t. external grid.

!-----------------------------------------------------------------------
!
! Internal variables
!
!-----------------------------------------------------------------------

  INTEGER :: i,j,k
  INTEGER :: i_ext,j_ext,k_ext

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


!-----------------------------------------------------------------------
!
! Find i and j indices of interpolation points in external grid.
!
!-----------------------------------------------------------------------

  CALL setijloc(nx,ny,nx_ext,ny_ext,x2d,y2d,x_ext,y_ext,i2d,j2d)

!-----------------------------------------------------------------------
!
! Find the i_ext,j_ext,k_ext anchor points for each i,j,k point.
!
!-----------------------------------------------------------------------

  DO k = 1,nz
    DO j = 1,ny
      DO i = 1,nx

!       Brute force method

        i_ext = i2d(i,j)
        j_ext = j2d(i,j)

        DO k_ext = 1,nz_ext-1

          IF ( z3d(i,j,k) <= z3d_ext(i_ext  ,j_ext  ,k_ext+1  ) .AND.   &
               z3d(i,j,k) > z3d_ext(i_ext  ,j_ext  ,k_ext  ))THEN
            k3d(i,j,k) = k_ext

!            WRITE(6,*) 'Found ...',i,j,k,k_ext

            EXIT
          ELSE IF ( z3d(i,j,k) == z3d_ext(i_ext  ,j_ext  ,k_ext  ) .AND. &
                    k_ext == 1  )THEN
            k3d(i,j,k) = k_ext

!            WRITE(6,*) 'Found ...',i,j,k,k_ext

            EXIT
          ELSE IF ( z3d(i,j,k) == z3d_ext(i_ext  ,j_ext  ,k_ext+1  ) .AND. &
                    k_ext+1 == nz_ext  )THEN

            k3d(i,j,k) = k_ext

!            WRITE(6,*)'Found ...',i,j,k,k_ext

            EXIT
          ELSE

!            WRITE(6,*) 'Not found yet...',k_ext

          END IF

        END DO ! DO k_ext = 1,nz_ext-1

!      WRITE(6,*)'i,j,k,k3d = ',i,j,k,k3d(i,j,k)

      END DO ! DO i = 1,nx-1
    END DO ! DO j = 1,ny-1
  END DO ! DO k = 1,nz-1

  RETURN
END SUBROUTINE setijkloc3d


!##################################################################
!##################################################################
!######                                                      ######
!######           SUBROUTINE get_nstyps_from_sfcdat          ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################

SUBROUTINE get_nstyps_from_sfcdat(nstyps,sfcfile)

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! Determines number of soil types in a sfcdata file.
!
!-----------------------------------------------------------------------
!
! HISTORY:
!
! Written 8 June 2002 by Eric Kemp
!
!-----------------------------------------------------------------------
!
! Force explicit declarations
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------
!
! Include files
!
!-----------------------------------------------------------------------

  INCLUDE 'globcst.inc'
  INCLUDE 'mp.inc'            ! Message passing parameters

!-----------------------------------------------------------------------
!
! Declare arguments
!
!-----------------------------------------------------------------------

  INTEGER :: nstyps
  CHARACTER (LEN=*) :: sfcfile

!-----------------------------------------------------------------------
!
! Local variables
!
!-----------------------------------------------------------------------

  CHARACTER (LEN=256) :: savename
  INTEGER :: ierr,istat,idummy
  INTEGER :: nstyp1
  INTEGER :: sd_id

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  IF (mp_opt > 0) THEN
    savename(1:256) = sfcfile(1:256)
    CALL gtsplitfn(savename,1,1,loc_x,loc_y,1,1,0,0,1,2,sfcfile,istat)
  END IF

  WRITE (6,*) "GET_NSTYPS_FROM_SFCDAT: reading in external supplied surface", &
              "data from file ",trim(sfcfile)

!-----------------------------------------------------------------------
!
! Read in necessary header information.
!
!-----------------------------------------------------------------------

  IF (sfcfmt == 1) THEN ! Unformatted Fortran binary

    CALL getunit( sfcunit )

    CALL asnctl ('NEWLOCAL', 1, ierr)
    CALL asnfile(sfcfile, '-F f77 -N ieee', ierr)

    OPEN(UNIT=sfcunit,FILE=TRIM(sfcfile),FORM='unformatted',            &
         STATUS='old',IOSTAT=istat)

    IF( istat /= 0 ) THEN

      WRITE(6,'(/1x,a,i2,/1x,a/)')                                      &
          'Error occured when opening the surface data file '           &
          //sfcfile//' using FORTRAN unit ',sfcunit,                    &
          ' Program stopped in GET_NSTYPS_FROM_SFCDAT.'
      CALL arpsstop("arpsstop called from GET_NSTYPS_FROM_SFCDAT opening file",1)

    END IF

    WRITE(6,'(/1x,a,/1x,a,i2/)')                                        &
        'This run will start from an external supplied surface ',       &
        'data file '//sfcfile//' using FORTRAN unit ',sfcunit

    READ (sfcunit,ERR=998) idummy,idummy

    READ (sfcunit,ERR=998) idummy,idummy,idummy,idummy,idummy,          &
                         idummy, idummy,nstyp1,idummy,idummy,           &
                         idummy, idummy,idummy,idummy,idummy,           &
                         idummy, idummy,idummy,idummy,idummy

    CLOSE ( sfcunit )
    CALL retunit ( sfcunit )

    nstyps = MAX(nstyp1,1)

  ELSE IF (sfcfmt == 3) THEN ! HDF4 format

    CALL hdfopen(trim(sfcfile), 1, sd_id)
    IF (sd_id < 0) THEN
      WRITE (6,*) "GET_NSTYPS_FROM_SFCDAT: ERROR opening ",             &
                 trim(sfcfile)," for reading."
      GO TO 998
    END IF

    CALL hdfrdi(sd_id,"nstyp",nstyp1,istat)
    IF (istat > 1) GO TO 998
    IF (istat == 0) THEN
      nstyps = MAX(nstyp1,1)
    ELSE
      WRITE (6, '(a)') 'Variable nstyp is not in the data set.'
    END IF

    CALL hdfclose(sd_id, istat)

  ELSE IF (sfcfmt == 2) THEN ! Direct output (Jerry Brotzge)

    nstyp = 1

  ELSE

    ! alternate dump format ...
    WRITE(6,*) 'The supported Surface data format are ',                &
               'binary (sfcfmt=1) and HDF4 no compressed (sfcfmt = 3).'
    CALL arpsstop('Surface data format is not supported.',1)

  END IF

  IF (mp_opt > 0) sfcfile(1:256) = savename(1:256)

  RETURN

!-----------------------------------------------------------------------
!
! Error handling
!
!-----------------------------------------------------------------------

  998   WRITE (6,'(/a,i2/a)')                                           &
         'GET_NSTYPS_FROM_SFCDAT: Read error in surface data file '     &
         //sfcfile//' with the I/O unit ',sfcunit,                      &
         'The model will STOP in subroutine GET_DIMS_FROM_SFCDAT.'

  CALL arpsstop("arpsstop called from GET_NSTYPS_FROM_SFCDAT reading sfc file",1)

END SUBROUTINE get_nstyps_from_sfcdat
