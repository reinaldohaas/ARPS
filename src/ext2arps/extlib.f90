!
!##################################################################
!##################################################################
!######                                                      ######
!######               SUBROUTINE MKARPSVAR                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE mkarpsvar(nx_ext,ny_ext,nz_ext,nx,ny,nz,lvlprof,             &
           iorder,iprtopt,intropt,iloc,jloc,x_ext,y_ext,                &
           hgtext,zext,x2d,y2d,za,varext,                               &
           zsnd,varsnd,arpsbar,arpsprt,                                 &
           dxfld,dyfld,rdxfld,rdyfld,                                   &
           slopey,alphay,betay,                                         &
           tem_ext)
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
!
!  11/21/1994 (KB)
!  Added full documentation.
!
!  2000/08/16 (Gene Bassett)
!  Added intropt.
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

  INCLUDE 'bndry.inc'
  INCLUDE 'mp.inc'
!
!-----------------------------------------------------------------------
!
!  Input variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: nx_ext,ny_ext,nz_ext  ! external grid dimensions
  INTEGER :: nx,ny,nz              ! ARPS grid dimensions
  INTEGER :: lvlprof               ! levels in mean sounding
  INTEGER :: iorder                ! interpolating polynomial order
  INTEGER :: intropt               ! interpolation option
  INTEGER :: iprtopt               ! perturbation generation flag
  INTEGER :: iloc(nx,ny)           ! external x-index of ARPS grid point
  INTEGER :: jloc(nx,ny)           ! external y-index of ARPS grid point
  REAL :: x_ext(nx_ext)            ! external x-coord
  REAL :: y_ext(ny_ext)            ! external y-coord
  REAL :: hgtext(nx_ext,ny_ext,nz_ext)     ! heights of external levels
  REAL :: zext(nx,ny,nz_ext)       ! heights of external levels
!                                  interpolated to ARPS grid locs
  REAL :: x2d(nx,ny)
  REAL :: y2d(nx,ny)
  REAL :: za(nx,ny,nz)             ! ARPS physical heights
  REAL :: varext(nx_ext,ny_ext,nz_ext)     ! variable to convert
  REAL :: zsnd(lvlprof)            ! 1-D array of level heights
  REAL :: varsnd(lvlprof)          ! 1-D array of level-means
!
!-----------------------------------------------------------------------
!
!  Output variables
!
!-----------------------------------------------------------------------
!
  REAL :: arpsbar( nx, ny, nz)     ! 3-D array of level-means
  REAL :: arpsprt( nx, ny, nz)     ! Output array, perturnbation variable
                                   ! or total variable (see iprtopt)
!
!-----------------------------------------------------------------------
!
!  Temporary work arrays
!
!-----------------------------------------------------------------------
!
  REAL :: dxfld(nx_ext)
  REAL :: dyfld(ny_ext)
  REAL :: rdxfld(nx_ext)
  REAL :: rdyfld(ny_ext)
  REAL :: slopey(nx_ext,ny_ext,nz_ext)
  REAL :: alphay(nx_ext,ny_ext,nz_ext)
  REAL :: betay(nx_ext,ny_ext,nz_ext)
  REAL :: tem_ext(nx_ext,ny_ext,nz_ext)
  REAL :: tem1(nx,ny,nz)
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k,ia,ja,ka,kl
  REAL :: wlow
  REAL :: topprt,botprt,arpstop,arpsbot
  REAL :: pntint2d
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
    DO j=1,ny_ext
      DO i=1,nx_ext
        DO k=1,nz_ext
          DO kl=2,lvlprof-1
            IF(zsnd(kl) > hgtext(i,j,k)) EXIT
          END DO
          wlow=(zsnd(kl)-hgtext(i,j,k))/                                &
               (zsnd(kl)-zsnd(kl-1))
          tem_ext(i,j,k)=varext(i,j,k)-                                   &
               ((1.-wlow)*varsnd(kl) + wlow*varsnd(kl-1))
        END DO
      END DO
    END DO
  ELSE
    tem_ext = varext
  END IF
!
!-----------------------------------------------------------------------
!
!  Compute derivative terms
!
!-----------------------------------------------------------------------
!
  CALL setdrvy(nx_ext,ny_ext,nz_ext,                                    &
               1,nx_ext,1,ny_ext,1,nz_ext,                              &
               dyfld,rdyfld,tem_ext,                                    &
               slopey,alphay,betay)
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
!        51     CONTINUE
        wlow=(zsnd(kl)-za(ia,ja,ka))/                                   &
             (zsnd(kl)-zsnd(kl-1))
        arpsbar(ia,ja,ka)=                                              &
              (1.-wlow)*varsnd(kl) + wlow*varsnd(kl-1)
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
      botprt=pntint2d(nx_ext,ny_ext,                                    &
               2,nx_ext-1,2,ny_ext-1,                                   &
               iorder,x_ext,y_ext,x2d(ia,ja),y2d(ia,ja),                &
               iloc(ia,ja),jloc(ia,ja),tem_ext(1,1,1),                  &
               dxfld,dyfld,rdxfld,rdyfld,                               &
               slopey(1,1,1),alphay(1,1,1),betay(1,1,1))

      topprt=pntint2d(nx_ext,ny_ext,                                    &
               2,nx_ext-1,2,ny_ext-1,                                   &
               iorder,x_ext,y_ext,x2d(ia,ja),y2d(ia,ja),                &
               iloc(ia,ja),jloc(ia,ja),tem_ext(1,1,nz_ext),             &
               dxfld,dyfld,rdxfld,rdyfld,                               &
               slopey(1,1,nz_ext),                                      &
               alphay(1,1,nz_ext),betay(1,1,nz_ext))

      DO ka=1,nz
        IF(za(ia,ja,ka) < zext(ia,ja,1)) THEN
          arpsprt(ia,ja,ka)=botprt
        ELSE IF(za(ia,ja,ka) > zext(ia,ja,nz_ext)) THEN
          arpsprt(ia,ja,ka)=topprt
        ELSE
          DO kl=2,nz_ext-1
            IF(zext(ia,ja,kl) > za(ia,ja,ka)) EXIT
          END DO
!          351       CONTINUE
          wlow=(zext(ia,ja,kl)-   za(ia,ja,ka))/                        &
               (zext(ia,ja,kl)-zext(ia,ja,kl-1))

          arpstop=pntint2d(nx_ext,ny_ext,                               &
               2,nx_ext-1,2,ny_ext-1,                                   &
               iorder,x_ext,y_ext,x2d(ia,ja),y2d(ia,ja),                &
               iloc(ia,ja),jloc(ia,ja),tem_ext(1,1,kl),                 &
               dxfld,dyfld,rdxfld,rdyfld,                               &
               slopey(1,1,kl),                                          &
               alphay(1,1,kl),betay(1,1,kl))

          arpsbot=pntint2d(nx_ext,ny_ext,                               &
               2,nx_ext-1,2,ny_ext-1,                                   &
               iorder,x_ext,y_ext,x2d(ia,ja),y2d(ia,ja),                &
               iloc(ia,ja),jloc(ia,ja),tem_ext(1,1,kl-1),               &
               dxfld,dyfld,rdxfld,rdyfld,                               &
               slopey(1,1,kl-1),                                        &
               alphay(1,1,kl-1),betay(1,1,kl-1))

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
  IF(iprtopt == 0 .and. intropt == 1) THEN
    arpsprt = arpsprt + arpsbar
  ELSE IF(iprtopt == 1 .and. intropt == 2) THEN
    arpsprt = arpsprt - arpsbar
  END IF

  CALL mpsendrecv2dew(arpsbar, nx, ny, nz, ebc, wbc, 1, tem1)
  CALL mpsendrecv2dew(arpsprt, nx, ny, nz, ebc, wbc, 1, tem1)

  CALL mpsendrecv2dns(arpsbar, nx, ny, nz, nbc, sbc, 1, tem1)
  CALL mpsendrecv2dns(arpsprt, nx, ny, nz, nbc, sbc, 1, tem1)

!
  RETURN
END SUBROUTINE mkarpsvar
!
!##################################################################
!##################################################################
!######                                                      ######
!######               SUBROUTINE MKARPSVAR                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE mkarpsvar1(nx_ext,ny_ext,nz_ext,nx,ny,nz,lvlprof,            &
           iorder,iprtopt,intropt,iloc,jloc,x_ext,y_ext,                &
           hgtext,zext,x2d,y2d,za,varext,                               &
           zsnd,varsnd,arpsbar,arpsprt,                                 &
           trn_ext,var_h0,h0,                                           &
           dxfld,dyfld,rdxfld,rdyfld,                                   &
           slopey,alphay,betay,                                         &
           tem_ext)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Take data from external file on pressure levels and interpolate in
!  the horizontal and vertical to ARPS height levels and horizontal
!  locations. Then, form the ARPS mean and perturbation quantities.
!  Near surface fields are used in low level interpolation.
!  (Written from Keith Brewster's original mkarpsvar code)
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Fanyou Kong   (Modified based on Keith Brewster's mkarpsvar)
!  10/25/2003
!
!  MODIFICATION HISTORY:
!
!  05/20/2004 (Yunheng Wang)
!  Removed the unused parameter extsfcopt.
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
!    trn_ext  External grid terrain height
!    var_h0   Near surface variable at h0 height above ground
!    h0       Near surface variable height (m)
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
  INTEGER :: nx_ext,ny_ext,nz_ext  ! external grid dimensions
  INTEGER :: nx,ny,nz              ! ARPS grid dimensions
  INTEGER :: lvlprof               ! levels in mean sounding
  INTEGER :: iorder                ! interpolating polynomial order
  INTEGER :: intropt               ! interpolation option
  INTEGER :: iprtopt               ! perturbation generation flag
  INTEGER :: iloc(nx,ny)           ! external x-index of ARPS grid point
  INTEGER :: jloc(nx,ny)           ! external y-index of ARPS grid point
  REAL :: x_ext(nx_ext)            ! external x-coord
  REAL :: y_ext(ny_ext)            ! external y-coord
  REAL :: hgtext(nx_ext,ny_ext,nz_ext)     ! heights of external levels
  REAL :: zext(nx,ny,nz_ext)       ! heights of external levels
                                   ! interpolated to ARPS grid locs
  REAL :: x2d(nx,ny)
  REAL :: y2d(nx,ny)
  REAL :: za(nx,ny,nz)             ! ARPS physical heights
  REAL :: varext(nx_ext,ny_ext,nz_ext)     ! variable to convert
  REAL :: zsnd(lvlprof)            ! 1-D array of level heights
  REAL :: varsnd(lvlprof)          ! 1-D array of level-means

  REAL :: trn_ext(nx_ext,ny_ext)   ! external terrain height (m)
  REAL :: var_h0(nx_ext,ny_ext)    ! near surface values
  REAL :: h0                       ! near surface varaiable height (m)
  REAL :: tem_var_h0(nx_ext,ny_ext)
!
!-----------------------------------------------------------------------
!
!  Output variables
!
!-----------------------------------------------------------------------
!
  REAL :: arpsbar( nx, ny, nz)     ! 3-D array of level-means
  REAL :: arpsprt( nx, ny, nz)     ! Output array, perturnbation variable
                                   ! or total variable (see iprtopt)
!
!-----------------------------------------------------------------------
!
!  Temporary work arrays
!
!-----------------------------------------------------------------------
!
  REAL :: dxfld(nx_ext)
  REAL :: dyfld(ny_ext)
  REAL :: rdxfld(nx_ext)
  REAL :: rdyfld(ny_ext)
  REAL :: slopey(nx_ext,ny_ext,nz_ext)
  REAL :: alphay(nx_ext,ny_ext,nz_ext)
  REAL :: betay(nx_ext,ny_ext,nz_ext)
  REAL :: tem_ext(nx_ext,ny_ext,nz_ext)
  REAL :: tem_1(nx_ext,ny_ext)
  REAL :: tem_2(nx_ext,ny_ext)
  REAL :: tem_3(nx_ext,ny_ext)
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k,ia,ja,ka,kl,ks
  REAL :: wlow,zsfc
  REAL :: topprt,botprt,arpstop,arpsbot
  REAL :: pntint2d
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
    DO j=1,ny_ext
      DO i=1,nx_ext
        DO k=1,nz_ext
          DO kl=2,lvlprof-1
            IF(zsnd(kl) > hgtext(i,j,k)) EXIT
          END DO
!          11     CONTINUE
          wlow=(zsnd(kl)-hgtext(i,j,k))/                                  &
               (zsnd(kl)-zsnd(kl-1))
          tem_ext(i,j,k)=varext(i,j,k)-                                   &
               ((1.-wlow)*varsnd(kl) + wlow*varsnd(kl-1))
        END DO

! processing near surface variable

        zsfc=trn_ext(i,j) + h0
        DO kl=2,lvlprof-1
          IF(zsnd(kl) > zsfc) EXIT
        END DO
        wlow=(zsnd(kl)-zsfc)/(zsnd(kl)-zsnd(kl-1))
        tem_var_h0(i,j)=var_h0(i,j)-                                &
              ((1.-wlow)*varsnd(kl) + wlow*varsnd(kl-1))

      END DO
    END DO
  ELSE
    tem_ext = varext
    tem_var_h0=var_h0
  ENDIF
!
!-----------------------------------------------------------------------
!
!  Compute derivative terms
!
!-----------------------------------------------------------------------
!
  CALL setdrvy(nx_ext,ny_ext,nz_ext,                                    &
               1,nx_ext,1,ny_ext,1,nz_ext,                              &
               dyfld,rdyfld,tem_ext,                                    &
               slopey,alphay,betay)

! processing near surface variable

  CALL setdrvy(nx_ext,ny_ext,1,                                         &
               1,nx_ext,1,ny_ext,1,1,                                   &
               dyfld,rdyfld,tem_var_h0,                                 &
               tem_1,tem_2,tem_3)
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
        wlow=(zsnd(kl)-za(ia,ja,ka))/                                   &
             (zsnd(kl)-zsnd(kl-1))
        arpsbar(ia,ja,ka)=                                              &
              (1.-wlow)*varsnd(kl) + wlow*varsnd(kl-1)
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
!      botprt=pntint2d(nx_ext,ny_ext,                                    &
!               2,nx_ext-1,2,ny_ext-1,                                   &
!               iorder,x_ext,y_ext,x2d(ia,ja),y2d(ia,ja),                &
!               iloc(ia,ja),jloc(ia,ja),tem_ext(1,1,1),                  &
!               dxfld,dyfld,rdxfld,rdyfld,                               &
!               slopey(1,1,1),alphay(1,1,1),betay(1,1,1))

      topprt=pntint2d(nx_ext,ny_ext,                                    &
               2,nx_ext-1,2,ny_ext-1,                                   &
               iorder,x_ext,y_ext,x2d(ia,ja),y2d(ia,ja),                &
               iloc(ia,ja),jloc(ia,ja),tem_ext(1,1,nz_ext),             &
               dxfld,dyfld,rdxfld,rdyfld,                               &
               slopey(1,1,nz_ext),                                      &
               alphay(1,1,nz_ext),betay(1,1,nz_ext))

! find first external layer index (ks) above zsfc

      zsfc=0.5*(za(ia,ja,1)+za(ia,ja,2)) + h0
      DO ks=1,nz_ext
        IF(zext(ia,ja,ks) > zsfc) EXIT  ! ks determined
      END DO

! vertical loop (arps domain)

      DO ka=2,nz
!        IF(za(ia,ja,ka) < zext(ia,ja,1)) THEN
!          arpsprt(ia,ja,ka)=botprt
!        ELSE IF(za(ia,ja,ka) > zext(ia,ja,nz_ext)) THEN
        IF(za(ia,ja,ka) > zext(ia,ja,nz_ext)) THEN
          arpsprt(ia,ja,ka)=topprt
        ELSE
          DO kl=1,nz_ext-1
            IF(zext(ia,ja,kl) > za(ia,ja,ka)) EXIT
          END DO
          arpstop=pntint2d(nx_ext,ny_ext,                               &
               2,nx_ext-1,2,ny_ext-1,                                   &
               iorder,x_ext,y_ext,x2d(ia,ja),y2d(ia,ja),                &
               iloc(ia,ja),jloc(ia,ja),tem_ext(1,1,kl),                 &
               dxfld,dyfld,rdxfld,rdyfld,                               &
               slopey(1,1,kl),                                          &
               alphay(1,1,kl),betay(1,1,kl))

          IF(za(ia,ja,ka) >= zext(ia,ja,ks)) THEN  ! using external layer
          wlow=(zext(ia,ja,kl)-   za(ia,ja,ka))/                        &
               (zext(ia,ja,kl)-zext(ia,ja,kl-1))

          arpsbot=pntint2d(nx_ext,ny_ext,                               &
               2,nx_ext-1,2,ny_ext-1,                                   &
               iorder,x_ext,y_ext,x2d(ia,ja),y2d(ia,ja),                &
               iloc(ia,ja),jloc(ia,ja),tem_ext(1,1,kl-1),               &
               dxfld,dyfld,rdxfld,rdyfld,                               &
               slopey(1,1,kl-1),                                        &
               alphay(1,1,kl-1),betay(1,1,kl-1))

          ELSE      ! using near surface value

            wlow=(zext(ia,ja,kl)-   za(ia,ja,ka))/                      &
                 (zext(ia,ja,kl)-zsfc)

! The line below reflects rare situation when za(ka)<zsfc (occurring for u,v)
            IF ( za(ia,ja,ka) <= zsfc ) wlow=1.0

          arpsbot=pntint2d(nx_ext,ny_ext,                               &
               2,nx_ext-1,2,ny_ext-1,                                   &
               iorder,x_ext,y_ext,x2d(ia,ja),y2d(ia,ja),                &
               iloc(ia,ja),jloc(ia,ja),tem_var_h0,                      &
               dxfld,dyfld,rdxfld,rdyfld,                               &
               tem_1,tem_2,tem_3)

          END IF
          arpsprt(ia,ja,ka)=(1.-wlow)*arpstop+wlow*arpsbot

        END IF
      END DO

! Extropolating to k=1 level
      wlow=(za(ia,ja,3)-za(ia,ja,1))/(za(ia,ja,3)-za(ia,ja,2))
      arpsprt(ia,ja,1)=wlow*arpsprt(ia,ja,2)+                   &
                       (1.-wlow)*arpsprt(ia,ja,3)

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
    arpsprt = arpsprt + arpsbar
  ELSE IF(iprtopt == 1 .and. intropt == 2) THEN
    arpsprt = arpsprt - arpsbar
  END IF
!
  RETURN
END SUBROUTINE mkarpsvar1
!
!##################################################################
!##################################################################
!######                                                      ######
!######               SUBROUTINE MKARPSVLZ                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE mkarpsvlz(nx_ext,ny_ext,nz_ext,nx,ny,nz,lvlprof,             &
           iorder,iprtopt,intropt,iloc,jloc,x_ext,y_ext,                &
           hgtext,zext,x2d,y2d,za,varext,                               &
           zsnd,vlnsnd,arpsbar,arpsprt,                                 &
           dxfld,dyfld,rdxfld,rdyfld,                                   &
           slopey,alphay,betay,                                         &
           tem_ext)
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
!
!  11/21/1994 (KB)
!  Added full documentation.
!
!  2000/08/16 (Gene Bassett)
!  Added intropt.
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
  INTEGER :: nx_ext,ny_ext,nz_ext  ! external grid dimensions
  INTEGER :: nx,ny,nz              ! ARPS grid dimensions
  INTEGER :: lvlprof               ! levels in mean sounding
  INTEGER :: iorder                ! interpolating polynomial order
  INTEGER :: intropt               ! interpolation option
  INTEGER :: iprtopt               ! perturbation generation flag
  INTEGER :: iloc(nx,ny)           ! external x-index of ARPS grid point
  INTEGER :: jloc(nx,ny)           ! external y-index of ARPS grid point
  REAL :: x_ext(nx_ext)            ! external x-coord
  REAL :: y_ext(ny_ext)            ! external y-coord
  REAL :: hgtext(nx_ext,ny_ext,nz_ext)     ! heights of external levels
  REAL :: zext(nx,ny,nz_ext)       ! heights of external levels
!                                  interpolated to ARPS grid locs
  REAL :: x2d(nx,ny)
  REAL :: y2d(nx,ny)
  REAL :: za(nx,ny,nz)             ! ARPS physical heights
  REAL :: varext(nx_ext,ny_ext,nz_ext)     ! variable to convert
  REAL :: zsnd(lvlprof)            ! 1-D array of level heights
  REAL :: vlnsnd(lvlprof)          ! 1-D array of level-means
!
!-----------------------------------------------------------------------
!
!  Output variables
!
!-----------------------------------------------------------------------
!
  REAL :: arpsbar( nx, ny, nz)     ! 3-D array of level-means
  REAL :: arpsprt( nx, ny, nz)     ! Output array, perturnbation variable
                                   ! or total variable (see iprtopt)
!
!-----------------------------------------------------------------------
!
!  Temporary work arrays
!
!-----------------------------------------------------------------------
!
  REAL :: dxfld(nx_ext)
  REAL :: dyfld(ny_ext)
  REAL :: rdxfld(nx_ext)
  REAL :: rdyfld(ny_ext)
  REAL :: slopey(nx_ext,ny_ext,nz_ext)
  REAL :: alphay(nx_ext,ny_ext,nz_ext)
  REAL :: betay(nx_ext,ny_ext,nz_ext)
  REAL :: tem_ext(nx_ext,ny_ext,nz_ext)
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k,ia,ja,ka,kl
  REAL :: wlow
  REAL :: topprt,botprt,arpstop,arpsbot
  REAL :: pntint2d
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
    DO j=1,ny_ext
      DO i=1,nx_ext
        DO k=1,nz_ext
          DO kl=2,lvlprof-1
            IF(zsnd(kl) > hgtext(i,j,k)) EXIT
          END DO
!          11     CONTINUE
          wlow=(zsnd(kl)-hgtext(i,j,k))/                                  &
               (zsnd(kl)-zsnd(kl-1))
          tem_ext(i,j,k)=varext(i,j,k)-EXP(                               &
               ((1.-wlow)*vlnsnd(kl) + wlow*vlnsnd(kl-1)))
        END DO
      END DO
    END DO
  ELSE ! intropt=2
    tem_ext = LOG(varext)
  END IF
!
!-----------------------------------------------------------------------
!
!  Compute derivative terms
!
!-----------------------------------------------------------------------
!
  CALL setdrvy(nx_ext,ny_ext,nz_ext,                                    &
               1,nx_ext,1,ny_ext,1,nz_ext,                              &
               dyfld,rdyfld,tem_ext,                                    &
               slopey,alphay,betay)
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
!        51     CONTINUE
        wlow=(zsnd(kl)-za(ia,ja,ka))/                                   &
             (zsnd(kl)-zsnd(kl-1))
        arpsbar(ia,ja,ka)=EXP(                                          &
              (1.-wlow)*vlnsnd(kl) + wlow*vlnsnd(kl-1))
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
      botprt=pntint2d(nx_ext,ny_ext,                                    &
               2,nx_ext-1,2,ny_ext-1,                                   &
               iorder,x_ext,y_ext,x2d(ia,ja),y2d(ia,ja),                &
               iloc(ia,ja),jloc(ia,ja),tem_ext(1,1,1),                  &
               dxfld,dyfld,rdxfld,rdyfld,                               &
               slopey(1,1,1),alphay(1,1,1),betay(1,1,1))

      topprt=pntint2d(nx_ext,ny_ext,                                    &
               2,nx_ext-1,2,ny_ext-1,                                   &
               iorder,x_ext,y_ext,x2d(ia,ja),y2d(ia,ja),                &
               iloc(ia,ja),jloc(ia,ja),tem_ext(1,1,nz_ext),             &
               dxfld,dyfld,rdxfld,rdyfld,                               &
               slopey(1,1,nz_ext),                                      &
               alphay(1,1,nz_ext),betay(1,1,nz_ext))

      DO ka=1,nz
        IF(za(ia,ja,ka) < zext(ia,ja,1)) THEN
          arpsprt(ia,ja,ka)=botprt
        ELSE IF(za(ia,ja,ka) > zext(ia,ja,nz_ext)) THEN
          arpsprt(ia,ja,ka)=topprt
        ELSE
          DO kl=2,nz_ext-1
            IF(zext(ia,ja,kl) > za(ia,ja,ka)) EXIT
          END DO
!          351       CONTINUE
          wlow=(zext(ia,ja,kl)-   za(ia,ja,ka))/                        &
               (zext(ia,ja,kl)-zext(ia,ja,kl-1))

          arpstop=pntint2d(nx_ext,ny_ext,                               &
               2,nx_ext-1,2,ny_ext-1,                                   &
               iorder,x_ext,y_ext,x2d(ia,ja),y2d(ia,ja),                &
               iloc(ia,ja),jloc(ia,ja),tem_ext(1,1,kl),                 &
               dxfld,dyfld,rdxfld,rdyfld,                               &
               slopey(1,1,kl),                                          &
               alphay(1,1,kl),betay(1,1,kl))

          arpsbot=pntint2d(nx_ext,ny_ext,                               &
               2,nx_ext-1,2,ny_ext-1,                                   &
               iorder,x_ext,y_ext,x2d(ia,ja),y2d(ia,ja),                &
               iloc(ia,ja),jloc(ia,ja),tem_ext(1,1,kl-1),               &
               dxfld,dyfld,rdxfld,rdyfld,                               &
               slopey(1,1,kl-1),                                        &
               alphay(1,1,kl-1),betay(1,1,kl-1))

          arpsprt(ia,ja,ka)=(1.-wlow)*arpstop+wlow*arpsbot
        END IF
      END DO
    END DO
  END DO

  IF(intropt == 2) THEN
    arpsprt = EXP(arpsprt)
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
  IF(iprtopt == 0 .and. intropt == 1) THEN
    arpsprt = arpsprt + arpsbar
  ELSE IF(iprtopt == 1 .and. intropt == 2) THEN
    arpsprt = arpsprt - arpsbar
  END IF
!
  RETURN
END SUBROUTINE mkarpsvlz
!
!##################################################################
!##################################################################
!######                                                      ######
!######               SUBROUTINE MKARPS2D                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE mkarps2d (nx_ext,ny_ext,nx,ny,                               &
           iorder,iloc,jloc,x_ext,y_ext,                                &
           x2d,y2d,varext,arpsvar,                                      &
           dxfld,dyfld,rdxfld,rdyfld,                                   &
           slopey,alphay,betay,                                         &
           tem_ext)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Take 2D surface data from external file and interpolate in
!  the horizontal to form the ARPS 2D data set
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Fanyou Kong
!  6/10/1997.
!
!  MODIFICATION HISTORY:
!
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx_ext   Number of grid points in the x-direction (east/west)
!             for the external grid
!    ny_ext   Number of grid points in the y-direction (north/south)
!             for the external grid
!
!    nx       Number of grid points in the x-direction (east/west)
!             for the ARPS grid
!    ny       Number of grid points in the y-direction (north/south)
!             for the ARPS grid
!
!    iorder   order of polynomial for interpolation (1, 2 or 3)
!
!    iloc     x-index location of ARPS grid point in the external array
!    jloc     y-index location of ARPS grid point in the external array
!
!    x2d      x coordinate of ARPS grid point in external coordinate
!    y2d      x coordinate of ARPS grid point in external coordinate
!
!    varext   external 2D data
!
!  OUTPUT:
!
!    arpsvar  2D array of variable at ARPS grid locations
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
  INTEGER :: nx_ext,ny_ext         ! external grid dimensions
  INTEGER :: nx,ny                 ! ARPS grid dimensions
  INTEGER :: iorder                ! interpolating polynomial order
  INTEGER :: iloc(nx,ny)           ! external x-index of ARPS grid point
  INTEGER :: jloc(nx,ny)           ! external y-index of ARPS grid point
  REAL :: x_ext(nx_ext)            ! external x-coord
  REAL :: y_ext(ny_ext)            ! external y-coord
!                                  interpolated to ARPS grid locs
  REAL :: x2d(nx,ny)
  REAL :: y2d(nx,ny)
  REAL :: varext(nx_ext,ny_ext)    ! 2D variable to convert
!
!-----------------------------------------------------------------------
!
!  Output variables
!
!-----------------------------------------------------------------------
!
  REAL :: arpsvar( nx, ny)         ! 2-D array of ARPS grid data
!
!-----------------------------------------------------------------------
!
!  Temporary work arrays
!
!-----------------------------------------------------------------------
!
  REAL :: dxfld(nx_ext)
  REAL :: dyfld(ny_ext)
  REAL :: rdxfld(nx_ext)
  REAL :: rdyfld(ny_ext)
  REAL :: slopey(nx_ext,ny_ext)
  REAL :: alphay(nx_ext,ny_ext)
  REAL :: betay(nx_ext,ny_ext)
  REAL :: tem_ext(nx_ext,ny_ext)
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: ia,ja
  REAL :: arpsdata
  REAL :: pntint2d
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
!  Compute derivative terms
!
!-----------------------------------------------------------------------
!
  CALL setdrvy(nx_ext,ny_ext,1,                                         &
               1,nx_ext,1,ny_ext,1,1,                                   &
               dyfld,rdyfld,varext,                                     &
               slopey,alphay,betay)
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
!    Horizontal interpolation
!
!-----------------------------------------------------------------------
!
      arpsdata=pntint2d(nx_ext,ny_ext,                                  &
               2,nx_ext-1,2,ny_ext-1,                                   &
               iorder,x_ext,y_ext,x2d(ia,ja),y2d(ia,ja),                &
               iloc(ia,ja),jloc(ia,ja),varext,                          &
               dxfld,dyfld,rdxfld,rdyfld,                               &
               slopey,alphay,betay)

      arpsvar(ia,ja)=arpsdata

    END DO
  END DO

  RETURN
END SUBROUTINE mkarps2d

!
!##################################################################
!##################################################################
!######                                                      ######
!######              SUBROUTINE INITSOILEXT                  ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE initsoilext(nx,ny,nzsoil,nzsoil_ext, nstyp,zpsoil,  &
                       zpsoil_ext,tsoil,tsoil_ext,qsoil,qsoil_ext)

!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Take soil data from external file and interpolate in
!  the vertical to form the ARPS 2D data set profile
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Jerry Brotzge
!  05/15/2002
!
!  MODIFICATION HISTORY:
!
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!             for the ARPS grid
!    ny       Number of grid points in the y-direction (north/south)
!             for the ARPS grid
!    nzsoil   Number of grid points in the soil
!
!    nzsoil_ext   Number of grid points in the soil external file
!
!    nstyp    Number of soil types
!
!    zpsoil   Physical depth of soil layers (m)
!
!    zpsoil_ext Physical depth of external model soil layers (m)
!
!    tsoil    Soil temperature (K)
!
!    qsoil    Soil moisture (m**3/m**3)
!
!  OUTPUT:
!
!    tsoil    Soil temperature profile (K)
!
!    qsoil    Soil moisture profile (m**3/m**3)
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
!
!-----------------------------------------------------------------------
!
!  Input/output variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k,kk

  INTEGER :: nx,ny                 ! ARPS grid dimensions
  INTEGER :: nzsoil                ! ARPS soil levels
  INTEGER :: nzsoil_ext            ! External file soil levels
  INTEGER :: nstyp                 ! Number of soil types
!
  REAL :: zpsoil(nx,ny,nzsoil)     ! Physical depth of ARPS soil layers
  REAL :: zpsoil_ext(nx,ny,nzsoil_ext) ! Physical depth of ext. file soil layers

  REAL :: tsoil(nx,ny,nzsoil)      ! Soil temperature (K)
  REAL :: qsoil(nx,ny,nzsoil)      ! Soil moisture (m**3/m**3)

  REAL :: tsoil_ext(nx,ny,nzsoil_ext) ! External file soil temperature (K)
  REAL :: qsoil_ext(nx,ny,nzsoil_ext) ! External file soil moisture (m**3/m**3)


!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------

  REAL :: dampdepth

!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------

  INCLUDE 'grid.inc'

!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

!-----------------------------------------------------------------------
!
!    Vertical interpolation
!
!-----------------------------------------------------------------------

!----------------------------------------------------------------
!     Set top boundary condition at level = 1
!-------------------------------------------------------------

  DO j=1,ny
    DO i=1,nx

      tsoil (i,j,1) = tsoil_ext(i,j,1)
      qsoil (i,j,1) = qsoil_ext(i,j,1)

      DO k=1,nzsoil
        zpsoil(i,j,k) = zpsoil(i,j,k) + zrefsfc
      END DO
      DO kk=1,nzsoil_ext
        zpsoil_ext(i,j,kk) = zpsoil_ext(i,j,kk) + zrefsfc
      END DO

    END DO
  END DO

!----------------------------------------------------------------------
!   Initialize soil temp and moisture profiles
!----------------------------------------------------------------------

! Note: All soil depths are negative (downward)

  DO k=2,nzsoil
    DO j=1,ny
      DO i=1,nx

        dampdepth = -0.15 - zpsoil(i,j,k) !Typical damping depth = 15 cm

        DO kk=2,nzsoil_ext

!----------------------------------------------------------------------
!   Linear fit between initialized levels below damping depth (15 cm)
!----------------------------------------------------------------------

        IF (zpsoil(i,j,k) >= zpsoil_ext(i,j,kk) .AND.   &
            zpsoil(i,j,k) <= zpsoil_ext(i,j,kk-1)) THEN

            tsoil(i,j,k) = tsoil_ext(i,j,kk)+((tsoil_ext(i,j,kk)-   &
            tsoil_ext(i,j,kk-1))*  (zpsoil(i,j,k)/zpsoil_ext(i,j,kk)) )

            qsoil(i,j,k) = qsoil_ext(i,j,kk)+((qsoil_ext(i,j,kk)-   &
            qsoil_ext(i,j,kk-1))*  (zpsoil(i,j,k)/zpsoil_ext(i,j,kk)) )


!----------------------------------------------------------------------
!   Exponential fit to initialization above damping depth (15 cm)
!----------------------------------------------------------------------

        ELSE IF (zpsoil(i,j,k) > zpsoil_ext(i,j,kk)) THEN

            IF (zpsoil_ext(i,j,kk) < dampdepth) THEN

              IF (zpsoil(i,j,k) >= dampdepth) THEN

                tsoil(i,j,k)=tsoil(i,j,k-1)+(tsoil_ext(i,j,kk)-tsoil(i,j,k-1))* &
                  EXP( - (zpsoil(i,j,k)/dampdepth) )

                qsoil(i,j,k)=qsoil(i,j,k-1)+(qsoil_ext(i,j,kk)-qsoil(i,j,k-1))* &
                  EXP( - (zpsoil(i,j,k)/dampdepth) )

              ELSE IF (zpsoil(i,j,k) < dampdepth) THEN

                 tsoil(i,j,k) = tsoil(i,j,k-1)+((tsoil_ext(i,j,kk)-             &
                       tsoil(i,j,k-1))*  (zpsoil(i,j,k)/zpsoil_ext(i,j,kk)) )

                 qsoil(i,j,k) = qsoil(i,j,k-1)+((qsoil_ext(i,j,kk)-             &
                       qsoil(i,j,k-1))*  (zpsoil(i,j,k)/zpsoil_ext(i,j,kk)) )

               END IF

            ELSE IF (zpsoil_ext(i,j,kk) >= dampdepth) THEN

                tsoil(i,j,k)=tsoil(i,j,k-1)+(tsoil_ext(i,j,kk)-tsoil(i,j,k-1))* &
                  EXP( - (zpsoil(i,j,k)/zpsoil_ext(i,j,kk)) )

                qsoil(i,j,k)=qsoil(i,j,k-1)+(qsoil_ext(i,j,kk)-qsoil(i,j,k-1))* &
                  EXP( - (zpsoil(i,j,k)/zpsoil_ext(i,j,kk)) )

            END IF


!----------------------------------------------------------------------
!   Set constant below lowest initialized level
!----------------------------------------------------------------------

        ELSE IF (zpsoil(i,j,k) < zpsoil_ext(i,j,kk)) THEN

            tsoil(i,j,k) = tsoil_ext(i,j,kk)
            qsoil(i,j,k) = qsoil_ext(i,j,kk)

        END IF
        END DO

       END DO
     END DO
   END DO

  RETURN
END SUBROUTINE initsoilext
