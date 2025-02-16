!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE RDMPUVW                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE rdmpuvw(nx,ny,nz,exbcbufsz,                                  &
           u,v,w, ubar,vbar,rhostr, zp,                                 &
           umix,vmix,wmix,                                              &
           exbcbuf,uexbc,vexbc,wexbc,                                   &
           rdmp,urdmp,vrdmp,wrdmp)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Apply Rayleigh sponge to w, and to perturbations of u and v in the
!  momentum equations. The Rayleigh damping terms are then added to
!  arrays umix, vmix and wmix.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  4/21/1992.
!
!  MODIFICATION HISTORY:
!
!  5/05/92 (M. Xue)
!  Added full documentation.
!
!  6/2/92 (M. Xue and H. Jin)
!  Further facelift.
!
!  2/10/93 (K. Droegemeier)
!  Cleaned up documentation.
!
!  9/23/93 (M. Xue)
!  Rayleigh damping was extended to include the points on the lateral
!  boundaries. These values will be used in the case of radiation
!  lateral boundary conditions, where prognostic equations are
!  integrated on the y boundaries for u, x boundaries for v and
!  all four lateral boundaries for w and scalar (except for pprt)
!  variables.
!
!  9/10/94 (D. Weber & Y. Lu)
!  Cleaned up documentation.
!
!  5/7/96  (Donghai Wang and M. Xue)
!  Added a parameter, raydmp=2, for Rayleigh damping to relax the
!  total fields toward external fields defined in the EXBC file.
!  When raydmp=1, it still damps the difference between the total
!  and the base state fields.
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical direction.
!
!    u        x component of velocity (m/s)
!    v        y component of velocity (m/s)
!    w        Vertical component of velocity in Cartesian
!             coordinates (m/s).
!
!    ubar     Base state zonal velocity component (m/s)
!    vbar     Base state meridional velocity component (m/s)
!    rhostr   Base state density times j3 (kg/m**3)
!
!    zp       Vertical coordinate of grid points in physical space(m)
!
!    umix     Array containing turbulent and computational
!             mixing on u
!    vmix     Array containing turbulent and computational
!             mixing on v
!    wmix     Array containing turbulent and computational
!             mixing on w
!
!  OUTPUT:
!
!    umix     Total mixing including, Rayleigh damping on u
!    vmix     Total mixing including, Rayleigh damping on v
!    wmix     Total mixing including, Rayleigh damping on w
!
!  WORK ARRAYS:
!
!    rdmp     Rayleigh damping coefficient.
!             A temporary array defined locally.
!    urdmp    A temporary array
!    vrdmp    A temporary array
!    wrdmp    A temporary array
!    uexbc    A temporary array
!    vexbc    A temporary array
!    wexbc    A temporary array
!
!-----------------------------------------------------------------------
!

!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx,ny,nz          ! Number of grid points in 3 directions

  REAL :: u     (nx,ny,nz)     ! Total u-velocity (m/s)
  REAL :: v     (nx,ny,nz)     ! Total v-velocity (m/s)
  REAL :: w     (nx,ny,nz)     ! Total w-velocity (m/s)

  REAL :: ubar  (nx,ny,nz)     ! Base state u-velocity (m/s)
  REAL :: vbar  (nx,ny,nz)     ! Base state v-velocity (m/s)
  REAL :: rhostr(nx,ny,nz)     ! Base state air density times j3
                               ! (kg/m**3)

  REAL :: zp    (nx,ny,nz)     ! The physical height coordinate defined
                               ! at w-point of the staggered grid.

  REAL :: umix  (nx,ny,nz)     ! Total mixing in u-eq.
  REAL :: vmix  (nx,ny,nz)     ! Total mixing in v-eq.
  REAL :: wmix  (nx,ny,nz)     ! Total mixing in w-eq.

  REAL :: rdmp  (nx,ny,nz)     ! Rayleigh damping coefficient.
                               ! Temporary array defined locally.
  REAL :: urdmp (nx,ny,nz)     ! Temporary array defined locally.
  REAL :: vrdmp (nx,ny,nz)     ! Temporary array defined locally.
  REAL :: wrdmp (nx,ny,nz)     ! Temporary array defined locally.

  INTEGER :: exbcbufsz         ! EXBC buffer size
  REAL :: exbcbuf( exbcbufsz ) ! EXBC buffer array. It carrys different
                               ! contents for nested grids from base
                               ! grid

  REAL :: uexbc (nx,ny,nz)     ! Temporary array defined locally.
  REAL :: vexbc (nx,ny,nz)     ! Temporary array defined locally.
  REAL :: wexbc (nx,ny,nz)     ! Temporary array defined locally.
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i, j, k, ijk
  REAL :: pi, tema
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'exbc.inc'
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  IF( raydmp == 0 .OR. cfrdmp == 0.0 ) RETURN

  IF( zbrdmp > zp(2,2,nz-1) ) RETURN
!
!-----------------------------------------------------------------------
!
!  Calculate the Rayleigh damping coefficient defined at a w-point.
!
!-----------------------------------------------------------------------
!
  pi = 4.0*ATAN( 1.0 )

  DO k= rayklow,nz
    DO j=1,ny-1
      DO i=1,nx-1

        rdmp(i,j,k) = cfrdmp*                                           &
            (1.0-COS(pi*MIN(1.0,(zp(i,j,k)-zbrdmp)/(zp(i,j,nz-1)-zbrdmp)) &
            ))*0.5

      END DO
    END DO
  END DO
!
!-----------------------------------------------------------------------
!
!  Rayleigh damping relaxes the total fields toward the base state
!
!-----------------------------------------------------------------------
!
  IF (raydmp == 1) THEN     ! Using the base state fields

    DO k=rayklow,nz-2
      DO j=1,ny-1
        DO i=2,nx-1

          urdmp(i,j,k) = u(i,j,k)-ubar(i,j,k)

        END DO
      END DO
    END DO

    DO k=rayklow,nz-2
      DO j=2,ny-1
        DO i=1,nx-1

          vrdmp(i,j,k) = v(i,j,k)-vbar(i,j,k)

        END DO
      END DO
    END DO

    DO k=rayklow,nz-1
      DO j=1,ny-1
        DO i=1,nx-1

          wrdmp(i,j,k) = w(i,j,k)

        END DO
      END DO
    END DO

!-----------------------------------------------------------------------
!
!  Rayleigh damping relaxes the total fields toward external fields
!  defined in the EXBC file.
!
!-----------------------------------------------------------------------


  ELSE IF (raydmp == 2) THEN      ! Using EXBC fields

    IF ( mgrid == 1 ) THEN         ! for base grid
      tema = curtim - ( abstfcst0 - abstinit )
      DO k = rayklow-1, nz-1
        DO j = 1, ny-1
          DO i = 2, nx-1
            ijk = (k-1)*nx*ny + (j-1)*nx + i
            uexbc(i,j,k) = exbcbuf(nu0exb+ijk-1)                        &
                         + exbcbuf(nudtexb+ijk-1) * tema
          END DO
        END DO
      END DO

      DO k = rayklow-1, nz-1
        DO j = 2, ny-1
          DO i = 1, nx-1
            ijk = (k-1)*nx*ny + (j-1)*nx + i
            vexbc(i,j,k) = exbcbuf(nv0exb+ijk-1)                        &
                         + exbcbuf(nvdtexb+ijk-1) * tema
          END DO
        END DO
      END DO

      DO k = rayklow-1, nz
        DO j = 1, ny-1
          DO i = 1, nx-1
            ijk = (k-1)*nx*ny + (j-1)*nx + i
            wexbc(i,j,k) = exbcbuf(nw0exb+ijk-1)                        &
                         + exbcbuf(nwdtexb+ijk-1) * tema
          END DO
        END DO
      END DO

    ELSE

      DO k = rayklow-1, nz-1
        DO j = 1, ny-1
          DO i = 2, nx-1
            ijk = (k-1)*nx*ny + (j-1)*nx + i
            uexbc(i,j,k) = exbcbuf(ijk)
          END DO
        END DO
      END DO

      DO k = rayklow-1, nz-1
        DO j = 2, ny-1
          DO i = 1, nx-1
            ijk = (k-1)*nx*ny + (j-1)*nx + i
            vexbc(i,j,k) = exbcbuf(1*nx*ny*nz+ijk)
          END DO
        END DO
      END DO

      DO k = rayklow-1, nz
        DO j = 1, ny-1
          DO i = 1, nx-1
            ijk = (k-1)*nx*ny + (j-1)*nx + i
            wexbc(i,j,k) = exbcbuf(2*nx*ny*nz+ijk)
          END DO
        END DO
      END DO

    END IF

    DO k=rayklow,nz-2
      DO j=1,ny-1
        DO i=2,nx-1

          urdmp(i,j,k) = u(i,j,k)-uexbc(i,j,k)

        END DO
      END DO
    END DO

    DO k=rayklow,nz-2
      DO j=2,ny-1
        DO i=1,nx-1

          vrdmp(i,j,k) = v(i,j,k)-vexbc(i,j,k)

        END DO
      END DO
    END DO

    DO k=rayklow,nz-1
      DO j=1,ny-1
        DO i=1,nx-1

          wrdmp(i,j,k) = w(i,j,k)-wexbc(i,j,k)

        END DO
      END DO
    END DO

  END IF

!
!-----------------------------------------------------------------------
!
!  Calculate the Rayleigh damping term in the u-momentum equation and
!  add to array umix.
!
!-----------------------------------------------------------------------
!
  DO k=rayklow,nz-2
    DO j=1,ny-1
      DO i=2,nx-1

        umix(i,j,k)=umix(i,j,k) -                                       &
            (rdmp(i,j,k)+rdmp(i-1,j,k)+rdmp(i,j,k+1)+rdmp(i-1,j,k+1))*0.25* &
            (rhostr(i,j,k)+rhostr(i-1,j,k))*0.5*                        &
            urdmp(i,j,k)

      END DO
    END DO
  END DO
!
!-----------------------------------------------------------------------
!
!  Calculate the Rayleigh damping term in the v-momentum equation and
!  add to array vmix.
!
!-----------------------------------------------------------------------
!
  DO k=rayklow,nz-2
    DO j=2,ny-1
      DO i=1,nx-1

        vmix(i,j,k)=vmix(i,j,k) -                                       &
            (rdmp(i,j,k)+rdmp(i,j-1,k)+rdmp(i,j,k+1)+rdmp(i,j-1,k+1))*0.25* &
            (rhostr(i,j,k)+rhostr(i,j-1,k))*0.5*                        &
            vrdmp(i,j,k)

      END DO
    END DO
  END DO

!
!-----------------------------------------------------------------------
!
!  Calculate the Rayleigh damping term in the w-momentum equation and
!  add to array wmix.
!
!-----------------------------------------------------------------------
!
  DO k=rayklow,nz-1
    DO j=1,ny-1
      DO i=1,nx-1

        wmix(i,j,k)=wmix(i,j,k) - rdmp(i,j,k)*                          &
                    (rhostr(i,j,k)+rhostr(i,j,k-1))*0.5*                &
                    wrdmp(i,j,k)

      END DO
    END DO
  END DO

  RETURN
END SUBROUTINE rdmpuvw

!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE RDMPPT                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE rdmppt(nx,ny,nz, exbcbufsz,                                  &
           ptprt, rhostr, zp,                                           &
           ptmix,                                                       &
           exbcbuf,ptexbc,                                              &
           rdmp,ptrdmp)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Apply Rayleigh sponge to the perturbation potential temperature,
!  and accumulate the results in array ptmix.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  4/21/1992.
!
!  MODIFICATION HISTORY:
!
!  5/05/92 (M. Xue)
!  Added full documentation.
!
!  6/2/92 (M. Xue and H. Jin)
!  Further facelift.
!
!  2/10/93 (K. Droegemeier)
!  Cleaned up documentation.
!
!  9/23/93 (MX)
!  Rayleigh damping was extended to include the points on the lateral
!  boundaries. These values will be used in the case of radiation
!  lateral boundary conditions, where prognostic equations are
!  integrated on the y boundaries for u, x boundaries for v and
!  all four lateral boundaries for w and scalar (except for pprt)
!  variables.
!
!  9/1/94 (D. Weber & Y. Lu)
!  Cleaned up documentation
!
!  5/7/96  (Donghai Wang and M. Xue)
!  Added a parameter, raydmp=2, for Rayleigh damping to relax the
!  total fields toward external fields defined in the EXBC file.
!  When raydmp=1, it still damps the difference between the total
!  and the base state fields.
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical direction.
!
!    ptprt    Perturbation potential temperature (K)
!
!    rhostr   Base state air density times j3 (kg/m**3)
!
!    zp       Vertical coordinate of grid points in physical space
!             (m)
!
!  OUTPUT:
!
!    ptmix    Total mixing in potential temperature equation.
!
!
!  WORK ARRAYS:
!
!    rdmp     Rayleigh damping coefficient.
!             A temporary array defined locally.
!    ptrdmp    A temporary array
!    ptexbc    A temporary array
!
!-----------------------------------------------------------------------
!

!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx,ny,nz          ! Number of grid points in 3 directions

  REAL :: ptprt (nx,ny,nz)     ! Perturbation potential temperature (K)
  REAL :: rhostr(nx,ny,nz)     ! Base state air density times j3
                               ! (kg/m**3)
  REAL :: zp    (nx,ny,nz)     ! The physical height coordinate defined
                               ! at w-point of the staggered grid.

  REAL :: ptmix (nx,ny,nz)     ! Total mixing on potential temperature

  REAL :: rdmp  (nx,ny,nz)     ! Rayleigh damping coefficient.
  REAL :: ptrdmp  (nx,ny,nz)   ! A temporary array.

  INTEGER :: exbcbufsz         ! EXBC buffer size
  REAL :: exbcbuf( exbcbufsz ) ! EXBC buffer array. It carrys different
                               ! contents for nested grids from base
                               ! grid

  REAL :: ptexbc  (nx,ny,nz)   ! A temporary array.
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i, j, k, ijk
  REAL :: pi, zpmax,tema
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'exbc.inc'
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  IF( raydmp == 0 .OR. cfrdmp == 0.0 ) RETURN
!
!-----------------------------------------------------------------------
!
!  Calculate the Rayleigh damping coefficient defined at a w-point.
!
!-----------------------------------------------------------------------
!
  pi = 4.0*ATAN( 1.0 )

  DO k= rayklow,nz
    DO j=1,ny-1
      DO i=1,nx-1

        rdmp(i,j,k) = cfrdmp*                                           &
            (1.0-COS(pi*MIN(1.0,(zp(i,j,k)-zbrdmp)/(zp(i,j,nz-1)-zbrdmp)) &
            ))*0.5


      END DO
    END DO
  END DO

!
!-----------------------------------------------------------------------
!
!  Rayleigh damping relaxes the total fields toward the base state
!
!-----------------------------------------------------------------------


  IF (raydmp == 1) THEN     ! Using the base state fields

    DO k=rayklow,nz-2
      DO j=1,ny-1
        DO i=1,nx-1

          ptrdmp(i,j,k) = ptprt(i,j,k)

        END DO
      END DO
    END DO

!-----------------------------------------------------------------------
!
!  Rayleigh damping relaxes the total fields toward external fields
!  defined in the EXBC file
!
!-----------------------------------------------------------------------

  ELSE IF ( raydmp == 2. ) THEN     ! Using EXBC data

    IF ( mgrid == 1 ) THEN         ! for base grid
      tema = curtim - ( abstfcst0 - abstinit )

      DO k = 1, nz-1
        DO j = 1, ny-1
          DO i = 1, nx-1
            ijk = (k-1)*nx*ny + (j-1)*nx + i
            ptexbc(i,j,k) = exbcbuf(npt0exb+ijk-1)                      &
                          + exbcbuf(nptdtexb+ijk-1) * tema
          END DO
        END DO
      END DO
    ELSE
      DO k = 1, nz-1
        DO j = 1, ny-1
          DO i = 1, nx-1
            ijk = (k-1)*nx*ny + (j-1)*nx + i
            ptexbc(i,j,k) = exbcbuf(3*nx*ny*nz+ijk)
          END DO
        END DO
      END DO
    END IF

    DO k=rayklow,nz-2
      DO j=1,ny-1
        DO i=1,nx-1

          ptrdmp(i,j,k) = ptprt(i,j,k) - ptexbc(i,j,k)

        END DO
      END DO
    END DO

  END IF
!
!-----------------------------------------------------------------------
!
!  Calculate the Rayleigh damping term for the perturbation potential
!  temperature equation and add to array ptmix.
!
!-----------------------------------------------------------------------
!
  DO k=rayklow,nz-2
    DO j=1,ny-1
      DO i=1,nx-1

        ptmix(i,j,k)=ptmix(i,j,k) -                                     &
                    (rdmp(i,j,k)+rdmp(i,j,k+1))*0.5*                    &
                     rhostr(i,j,k)*                                     &
                    ptrdmp(i,j,k)

      END DO
    END DO
  END DO

  RETURN
END SUBROUTINE rdmppt
