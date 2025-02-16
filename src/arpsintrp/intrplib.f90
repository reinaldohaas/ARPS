!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINE INTRPX3D                  ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################

SUBROUTINE intrpx3d(ain,nx,is,ie, ny,js,je, nz,ks,ke, wgtx,ix,          &
                    intrphopt,aout,nx1,is1,ie1)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!  Perform interpolation in the first dimension
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  4/1/1999.
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER :: nx,is,ie, ny,js,je, nz,ks,ke
  INTEGER :: nx1,is1,ie1
  REAL    :: ain (nx ,ny,nz)
  REAL    :: aout(nx1,ny,nz)
  REAL    :: wgtx(nx1,3)
  INTEGER :: ix(nx1)
  INTEGER :: intrphopt
  INTEGER :: i,j,k

  IF(intrphopt == 1) THEN
    DO k=ks ,ke
      DO j=js ,je
        DO i=is1,ie1
          aout(i,j,k)=      wgtx(i,1) *ain(ix(i)  ,j,k)                 &
                      +(1.0-wgtx(i,1))*ain(ix(i)+1,j,k)
        END DO
      END DO
    END DO
  ELSE
    DO k=ks ,ke
      DO j=js ,je
        DO i=is1,ie1
          aout(i,j,k)=wgtx(i,1)*ain(ix(i)-1,j,k)                        &
                     +wgtx(i,2)*ain(ix(i)  ,j,k)                        &
                     +wgtx(i,3)*ain(ix(i)+1,j,k)
        END DO
      END DO
    END DO
  END IF

  RETURN
END SUBROUTINE intrpx3d

SUBROUTINE intrpy3d(ain,nx,is,ie, ny,js,je, nz,ks,ke, wgty,jy,          &
           intrphopt,aout,ny1,js1,je1)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!  Perform interpolation in the second dimension
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  4/1/1999.
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx,is,ie, ny,js,je, nz,ks,ke
  INTEGER :: ny1,js1,je1
  REAL :: ain (nx,ny ,nz)
  REAL :: aout(nx,ny1,nz)
  REAL :: wgty(ny1,3)
  INTEGER :: jy(ny1)
  INTEGER :: intrphopt
  INTEGER :: i,j,k

  IF(intrphopt == 1) THEN
    DO k=ks ,ke
      DO j=js1,je1
        DO i=is ,ie
          aout(i,j,k)=      wgty(j,1) *ain(i,jy(j)  ,k)                 &
                      +(1.0-wgty(j,1))*ain(i,jy(j)+1,k)
        END DO
      END DO
    END DO
  ELSE
    DO k=ks ,ke
      DO j=js1,je1
        DO i=is ,ie
          aout(i,j,k)=wgty(j,1)*ain(i,jy(j)-1,k)                        &
                     +wgty(j,2)*ain(i,jy(j)  ,k)                        &
                     +wgty(j,3)*ain(i,jy(j)+1,k)
        END DO
      END DO
    END DO
  END IF

  RETURN
END SUBROUTINE intrpy3d

SUBROUTINE intrpz3d(ain,nx,is,ie, ny,js,je, nz,ks,ke, wgtz,kz,          &
           intrpvopt,aout,nz1,ks1,ke1)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!  Perform interpolation in the third dimension
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  4/1/1999.
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx,is,ie, ny,js,je, nz,ks,ke
  INTEGER :: nz1,ks1,ke1
  REAL :: ain (nx,ny,nz)
  REAL :: aout(nx,ny,nz1)
  REAL :: wgtz(nx,ny,nz1,3)
  INTEGER :: kz(nx,ny,nz1)
  INTEGER :: intrpvopt
  INTEGER :: i,j,k

  IF(intrpvopt == 1) THEN
    DO k=ks1,ke1
      DO i=is ,ie
        DO j=js ,je
          aout(i,j,k)=      wgtz(i,j,k,1) *ain(i,j,kz(i,j,k)  )         &
                      +(1.0-wgtz(i,j,k,1))*ain(i,j,kz(i,j,k)+1)
        END DO
      END DO
    END DO
  ELSE
    DO k=ks1,ke1
      DO i=is ,ie
        DO j=js ,je
          aout(i,j,k)=wgtz(i,j,k,1)*ain(i,j,kz(i,j,k)-1)                &
                     +wgtz(i,j,k,2)*ain(i,j,kz(i,j,k)  )                &
                     +wgtz(i,j,k,3)*ain(i,j,kz(i,j,k)+1)
        END DO
      END DO
    END DO
  END IF

  RETURN
END SUBROUTINE intrpz3d

SUBROUTINE intrpxy3d(ain,nx,is,ie, ny,js,je, nz,ks,ke,                  &
           wgtx,ix,wgty,jy,intrphopt,                                   &
           aout,nx1,is1,ie1, ny1,js1,je1,                               &
           temx1yz)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!  Perform interpolation in the first and second dimensions
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  4/1/1999.
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx,is,ie, ny,js,je, nz,ks,ke
  INTEGER :: nx1,is1,ie1, ny1,js1,je1
  REAL :: ain (nx ,ny ,nz)
  REAL :: aout(nx1,ny1,nz)
  REAL :: wgtx(nx1,3),wgty(ny1,3)
  INTEGER :: ix(nx1),jy(ny1)
  INTEGER :: intrphopt

  REAL :: temx1yz(nx1,ny,nz)

  CALL intrpx3d(ain,nx,is,ie, ny,js,je, nz,ks,ke,                       &
                wgtx,ix,intrphopt, temx1yz,nx1,is1,ie1)

  CALL intrpy3d(temx1yz,nx1,is1,ie1, ny,js,je, nz,ks,ke,                &
                wgty,jy,intrphopt, aout,   ny1,js1,je1)

  RETURN
END SUBROUTINE intrpxy3d

SUBROUTINE intrpxyz3d(ain,nx,is,ie, ny,js,je, nz,ks,ke,                 &
           wgtx,ix,wgty,jy,wgtz,kz,                                     &
           intrphopt,intrpvopt,                                         &
           aout,nx1,is1,ie1, ny1,js1,je1, nz1,ks1,ke1,                  &
           temx1yz,temx1y1z)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!  Perform interpolation in all three dimensions
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  4/1/1999.
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx,is,ie, ny,js,je, nz,ks,ke
  INTEGER :: nx1,is1,ie1, ny1,js1,je1, nz1,ks1,ke1
  REAL :: ain (nx ,ny ,nz)
  REAL :: aout(nx1,ny1,nz1)
  REAL :: wgtx(nx1,3),wgty(ny1,3),wgtz(nx1,ny1,nz1,3)
  INTEGER :: ix(nx1),jy(ny1),kz(nx1,ny1,nz1)
  INTEGER :: intrphopt
  INTEGER :: intrpvopt

  REAL :: temx1yz (nx1,ny ,nz)
  REAL :: temx1y1z(nx1,ny1,nz)

  CALL intrpxy3d(ain, nx,is,ie, ny,js,je, nz,ks,ke,                     &
                 wgtx,ix,wgty,jy,intrphopt,                             &
                 temx1y1z, nx1,is1,ie1, ny1,js1,je1, temx1yz)

  CALL intrpz3d (temx1y1z, nx1,is1,ie1, ny1,js1,je1, nz,ks,ke,          &
                 wgtz,kz,intrpvopt,                                     &
                 aout, nz1,ks1,ke1)

  RETURN
END SUBROUTINE intrpxyz3d

!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINE INTRP_SOIL_Real           ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################

SUBROUTINE intrp_soil_real(nx,ny,nx1,ny1,nstyp,nstyp1,nzsoil,wx,wy,ix,jy,     &
                      tsoil,soiltyp,stypfrct,tsoil1)

!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Interpolate soil properties onto another grid.
!  real variables used in soil model: tsoil, qsoil,wetcanp
!
!  NOTE:  This should be used with the old ARPS Force-Restore Soil 
!         Model.  When using the new OUSoil model, use intrpsoil3d_avg
!         or intrpsoil3d_pst.  (Eric Kemp, 18 June 2002).
!
!  Based on intrp_soil
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Hu     
!  2006/08/29
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
!  tsoil        Deep soil temperature in data set (K)
!  soiltyp      Soil type in data set
!  stypfrct     Soil type fraction
!
!  OUTPUT:
!
!  tsoil1       Deep soil temperature in data set (K)
!  soiltyp1     Soil type in data set
!  stypfrct1    Soil type fraction
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER :: nx,ny,nx1,ny1
  INTEGER :: nstyp,nstyp1
  INTEGER :: nzsoil

  REAL :: wx(nx1)    ! Weight factor in x-direction
  REAL :: wy(ny1)    ! Weight factor in y-direction
  INTEGER :: ix(nx1) ! Old grid index (lower left) for new grid point
  INTEGER :: jy(ny1) ! Old grid index (lower left) for new grid point

  REAL :: tsoil  (nx,ny,nzsoil,0:nstyp)   ! Deep soil temperature (K)
  INTEGER :: soiltyp(nx,ny,nstyp)  ! Soil type in model domain
  REAL :: stypfrct(nx,ny,nstyp)

  REAL :: tsoil1  (nx1,ny1,nzsoil,0:nstyp1)      ! Deep soil temperature (K)
  INTEGER :: soiltyp1(nx1,ny1,nstyp1)     ! Soil type in model domain
  REAL :: stypfrct1(nx1,ny1,nstyp1)

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
  REAL :: tsoil1sum   (0:nsoiltyp)     ! Deep soil temperature (K)
  REAL :: stypfrct1sum(0:nsoiltyp)     ! Frction of soil type

  INTEGER :: i,j,i1,j1,is,ii
  REAL :: weight,maxweight,frctot
  REAL :: totweight(0:nsoiltyp)
  INTEGER :: maxtype
  INTEGER :: soiltype

  REAL    :: inverse
  INTEGER :: kk

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  tsoil1 (:,:,:,:) = 0

  DO kk=1,nzsoil
    soiltyp1 (:,:,:) = 0
    stypfrct1(:,:,:) = 0

    IF (nstyp == 1) THEN
      stypfrct = 1   ! make sure stypfrct is set for nstyp=1
      ! copy level 0 to level 1 if level 1 is undefined:
      IF (tsoil(1,1,kk,1)   <= 0) tsoil  (:,:,kk,1) = tsoil  (:,:,kk,0)
    ENDIF

    DO j1 = 1,ny1-1                 ! desired grid
      DO i1 = 1,nx1-1
        tsoil1sum = 0.
        stypfrct1sum = 0.
        maxweight = 0.
        totweight = 0.
        DO j = jy(j1), jy(j1)+1       ! external input grid
          DO i = ix(i1), ix(i1)+1
            weight =      wx(i1)*(1.+ix(i1)-i)*wy(j1)     *(1+jy(j1)-j)&
                   + (1.-wx(i1))*(i-ix(i1)   )*wy(j1)     *(1+jy(j1)-j)&
                   +      wx(i1)*(1.+ix(i1)-i)*(1.-wy(j1))*(j-jy(j1))  &
                   + (1.-wx(i1))*(i-ix(i1)   )*(1.-wy(j1))*(j-jy(j1))
            DO is=1,nstyp
              IF (stypfrct(i,j,is) > 0) THEN
                soiltype = soiltyp(i,j,is)
                IF (soiltype < 1 .or. soiltype > nsoiltyp) soiltype = 0
                stypfrct1sum(soiltyp(i,j,is)) = stypfrct1sum(soiltyp(i,j,is)) &
                                + weight*stypfrct(i,j,is)
                totweight(soiltyp(i,j,is)) = totweight(soiltyp(i,j,is)) + weight

                tsoil1sum(soiltyp(i,j,is)) = tsoil1sum(soiltyp(i,j,is))   &
                                  + weight*stypfrct(i,j,is)*tsoil(i,j,kk,is)
              END IF
            END DO
          END DO
        END DO

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
            tsoil1(i1,j1,kk,is) = tsoil1sum(maxtype) * inverse
            IF (nstyp /= 1) THEN
              stypfrct1(i1,j1,is)= stypfrct1sum(maxtype)/totweight(maxtype)
            ELSE
              stypfrct1(i1,j1,is)= stypfrct1sum(maxtype)
            END IF
          !wjm
            stypfrct1sum(maxtype) = 0.
          ELSE 
            IF (is == 1) THEN
              IF (.NOT. warned) THEN
                WRITE (6,*) 'INTRP_SOIL: WARNING, no soil type found, ',  &
                            'variables not assigned!'
                warned = .TRUE.
              END IF
            ELSE
              stypfrct1(i1,j1,is) = 0.
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
        tsoil1(i1,j1,kk,0)   = 0.
        DO is = 1,nstyp1
          tsoil1(i1,j1,kk,0) = tsoil1(i1,j1,kk,0)                               &
                           + stypfrct1(i1,j1,is)*tsoil1(i1,j1,kk,is)
        END DO
      END DO
    END DO

  ENDDO  ! kk

END SUBROUTINE intrp_soil_real

!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINE INTRP_SOIL_int            ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################

SUBROUTINE intrp_soil_int(nx,ny,nx1,ny1,nstyp,nstyp1,wx,wy,ix,jy, &
                      soiltyp,stypfrct,vegtyp,                    &
                      soiltyp1,stypfrct1,vegtyp1)
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Interpolate soil properties onto another grid. For soiltyp, stypfrct
!                        and vegtyp only
!
!  NOTE:  This should be used with the old ARPS Force-Restore Soil 
!         Model.  When using the new OUSoil model, use intrpsoil3d_avg
!         or intrpsoil3d_pst.  (Eric Kemp, 18 June 2002).
!
!  based on intrp_soil
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Hu     
!  2006/08/29
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
!  soiltyp      Soil type in data set
!  stypfrct     Soil type fraction
!  vegtyp       Vegetation type in data set
!
!  OUTPUT:
!
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

  REAL :: wx(nx1)    ! Weight factor in x-direction
  REAL :: wy(ny1)    ! Weight factor in y-direction
  INTEGER :: ix(nx1) ! Old grid index (lower left) for new grid point
  INTEGER :: jy(ny1) ! Old grid index (lower left) for new grid point

  INTEGER :: soiltyp(nx,ny,nstyp)  ! Soil type in model domain
  REAL :: stypfrct(nx,ny,nstyp)
  INTEGER :: vegtyp(nx,ny)

  INTEGER :: soiltyp1(nx1,ny1,nstyp1)     ! Soil type in model domain
  REAL :: stypfrct1(nx1,ny1,nstyp1)
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
  REAL :: stypfrct1sum(0:nsoiltyp)     ! Frction of soil type

  INTEGER :: i,j,i1,j1,is,ii
  REAL :: weight,maxweight,frctot
  REAL :: totweight(0:nsoiltyp)
  INTEGER :: maxtype
  INTEGER :: soiltype

  REAL    :: inverse

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  soiltyp1 (:,:,:) = 0
  stypfrct1(:,:,:) = 0
  vegtyp1(:,:) = 0

  DO j1 = 1,ny1-1                 ! desired grid
    DO i1 = 1,nx1-1
      stypfrct1sum = 0.
      maxweight = 0.
      totweight = 0.
      vegtyp1(i1,j1) = 0
      DO j = jy(j1), jy(j1)+1       ! external input grid
        DO i = ix(i1), ix(i1)+1
          weight =      wx(i1)*(1.+ix(i1)-i)*    wy(j1) *(1+jy(j1)-j)  &
                 + (1.-wx(i1))*(i-ix(i1)   )*    wy(j1) *(1+jy(j1)-j)  &
                 +      wx(i1)*(1.+ix(i1)-i)*(1.-wy(j1))*(j-jy(j1))    &
                 + (1.-wx(i1))*(i-ix(i1)   )*(1.-wy(j1))*(j-jy(j1))
          DO is=1,nstyp
            IF (stypfrct(i,j,is) > 0) THEN
              soiltype = soiltyp(i,j,is)
              IF (soiltype < 1 .or. soiltype > nsoiltyp) soiltype = 0
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
          IF (nstyp /= 1) THEN
            stypfrct1(i1,j1,is)= stypfrct1sum(maxtype)/totweight(maxtype)
          ELSE
            stypfrct1(i1,j1,is)= stypfrct1sum(maxtype)
          END IF
          !wjm
          stypfrct1sum(maxtype) = 0.
        ELSE 
          IF (is == 1) THEN
            IF (.NOT. warned) THEN
              WRITE (6,*) 'INTRP_SOIL: WARNING, no soil type found, ',  &
                          'variables not assigned!'
              warned = .TRUE.
            END IF
          ELSE
            stypfrct1(i1,j1,is) = 0.
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
    END DO
  END DO

END SUBROUTINE intrp_soil_int

