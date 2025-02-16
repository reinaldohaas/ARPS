!########################################################################
!########################################################################
!#########                                                      #########
!#########                SUBROUTINE adjust_refl                #########
!#########                                                      #########
!#########                     Developed by                     #########
!#########    Center for Analysis and Prediction of Storms      #########
!#########                University of Oklahoma                #########
!#########                                                      #########
!########################################################################
!########################################################################

SUBROUTINE adjust_refl(nx,ny,nz,refl,zs,et_nids,vil_nids)

!------------------------------------------------------------------------
!
! PURPOSE:
!
! Adjusts the remapped 3D reflectivity field using NIDS echo top and VIL
! data.  A linear function R = A*H + B is used as a first guess between
! the NIDS echo top product level and the top of the original remapped
! reflectivity field (R being reflectivity, A being a slope, H is height,
! and B is the intercept).  However, a quadratic function R = (A*H*H) +
! (B*H) + C is used if "too much" reflectivity is inserted by the linear
! function.
!
!------------------------------------------------------------------------
!
! AUTHOR:
!
! Eric Kemp, 6 September 2001.  Developed for WDT.
!
!------------------------------------------------------------------------
!
! Force explicit declarations.
!
!------------------------------------------------------------------------

  IMPLICIT NONE

!------------------------------------------------------------------------
!
! Declare arguments
!
!------------------------------------------------------------------------

  INTEGER, INTENT(IN) :: nx,ny,nz        ! Dimensions of grid
  REAL, INTENT(INOUT) :: refl(nx,ny,nz)  ! 3D reflectivity (dBZ)
  REAL, INTENT(IN)    :: zs(nx,ny,nz)    ! Height of scalar points (m MSL)
  REAL, INTENT(IN)    :: et_nids(nx,ny)  ! NIDS Echo top (m MSL)
  REAL, INTENT(IN)    :: vil_nids(nx,ny) ! NIDS VIL (kg m**-2)

!-----------------------------------------------------------------------
!
! Functions
!
!-----------------------------------------------------------------------

  REAL, EXTERNAL :: vil_88d

!-----------------------------------------------------------------------
!
! Internal parameters
!
!-----------------------------------------------------------------------

  REAL,    PARAMETER :: K_const = 3.44E-6
  REAL,    PARAMETER :: foursevenths = 4./7.
  REAL,    PARAMETER :: sevenfourths = 7./4.
  INTEGER, PARAMETER :: iterationmax = 40

!-----------------------------------------------------------------------
!
! Internal variables
!
!-----------------------------------------------------------------------

  REAL :: refl_adj(nx,ny,nz)
  REAL :: et_refl(nx,ny)
  INTEGER :: k_et_refl(nx,ny)
  INTEGER :: k_et_nids(nx,ny)
  REAL :: vil_refl(nx,ny)
  REAL :: vil_residual(nx,ny)

  REAL :: reflthresh

  REAL :: Z, r
  REAL :: slope, intercept
  REAL :: Mtop, Mr, M
  REAL :: reflthreshexp

  REAL :: maxresidual,minresidual
  INTEGER :: i,j,k
  INTEGER :: iteration

  REAL :: A, B, C
  REAL :: N, N1, N2, N3, N4
  REAL :: D, D1, D2, D3, D4

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!------------------------------------------------------------------------
!
! Calculate the echo top from the 3D reflectivity field.  Also,
! determine k-levels of NIDS echo tops.
!
!------------------------------------------------------------------------

  reflthresh = 18.5

  CALL calc_et_refl(nx,ny,nz,refl,reflthresh,zs,et_refl,k_et_refl)
  CALL calc_k_et(nx,ny,nz,zs,et_nids,k_et_nids)

  refl_adj(:,:,:) = refl(:,:,:)

  DO j = 1,ny-1
    DO i = 1,nx-1

      IF (et_nids(i,j) == 0.) THEN
!        WRITE(6,*)'WARNING:  No NIDS echo top.'
!        WRITE(6,*)'i: ',i,' j: ',j
        CYCLE ! Skip to next column
      ELSE IF (vil_nids(i,j) == 0.) THEN
!        WRITE(6,*)'WARNING:  No NIDS VIL.'
!        WRITE(6,*)'i: ',i,' j: ',j
        CYCLE ! Skip to next column
      ELSE IF (et_nids(i,j) == -9999. .OR. k_et_nids(i,j) == -9999) THEN
!        WRITE(6,*)'WARNING:  NIDS echo top value missing.'
!        WRITE(6,*)'i: ',i,' j: ',j
        CYCLE ! Skip to next column
      ELSE IF (vil_nids(i,j) == -9999.) THEN
!        WRITE(6,*)'WARNING:  NIDS VIL value missing.'
!        WRITE(6,*)'i: ',i,' j: ',j
        CYCLE ! Skip to next column
      ELSE IF (et_refl(i,j) == -9999. .OR. k_et_refl(i,j) == -9999) THEN
!        WRITE(6,*)'WARNING:  Estimated echo top value missing.'
!        WRITE(6,*)'i: ',i,' j: ',j
        CYCLE ! Skip to next column
      ELSE IF (vil_refl(i,j) == -9999.) THEN
!        WRITE(6,*)'WARNING:  Estimated VIL value missing.'
!        WRITE(6,*)'i: ',i,' j: ',j
        CYCLE ! Skip to next column
      ELSE IF (k_et_nids(i,j) == k_et_refl(i,j)) THEN
!        WRITE(6,*)'NOTE:  Echo top k levels are equal.  Skipping...'
!        WRITE(6,*)'i: ',i,' j: ',j
        CYCLE ! Skip to next column
      ELSE IF (k_et_nids(i,j) < k_et_refl(i,j)) THEN
!        WRITE(6,*)'WARNING:  NIDS echo top below estimated value.'
!        WRITE(6,*)'i: ',i,' j: ',j
!        WRITE(6,*)'et_nids: ',et_nids(i,j),' et_refl: ',et_refl(i,j)
        CYCLE ! Skip to next column
      END IF

      reflthresh = 18.5

      DO iteration = 1,iterationmax

        vil_refl(i,j) = vil_88d(nz,refl_adj(i,j,:),zs(i,j,:))

        IF ( (vil_refl(i,j) < (vil_nids(i,j) - 5.) ) .OR.                &
             (iteration == 1) ) THEN

!------------------------------------------------------------------------
!
!         To use the linear function, we must first calculate the slope
!         and intercept values.  With some algebra and calculus, it can
!         be shown that:
!
!         slope = (Mtop - Mr)/(et_nids - et_refl), and
!
!         intercept = Mr - (slope*et_refl)
!
!         where Mtop is M at et_nids and Mr is M at et_refl.  So, we need
!         to calculate Mtop and Mr, and then slope and intercept.
!
!------------------------------------------------------------------------

          Z = 10.**(refl(i,j,k_et_refl(i,j))*0.1)
          Mr = K_const*(Z**foursevenths)

          reflthreshexp = reflthresh*0.1
          Z = 10.**(reflthreshexp)
          Mtop = K_const*(Z**foursevenths)

          slope = (Mtop - Mr)/(et_nids(i,j) - et_refl(i,j))

          intercept = Mr - (slope*et_refl(i,j))

!------------------------------------------------------------------------
!
!         Now apply the linear function to the refl_adj field between
!         k_et_refl and k_ef_nids.
!
!------------------------------------------------------------------------

          DO k = k_et_refl(i,j),k_et_nids(i,j)
            M = (slope*zs(i,j,k)) + intercept

            Z = (M/K_const)**sevenfourths
            refl_adj(i,j,k) = 10.*ALOG10(Z)

          END DO ! k loop

          vil_refl(i,j) = vil_88d(nz,refl_adj(i,j,:),zs(i,j,:))

          reflthresh = reflthresh + 1.5

          IF (reflthresh > 55.) EXIT

        ELSE IF (vil_refl(i,j) > (vil_nids(i,j) + 5.) ) THEN

          reflthresh = 18.5

!------------------------------------------------------------------------
!
!         To use the quadratic function, we must first calculate the A,
!         B, and C coefficients.  By performing much algebra and
!         calculus, it can be shown that:
!
!         A = N/D
!
!           N = N1 + N2 + N3 + N4
!
!             N1 = 6.*r*(et_nids - et_refl)
!
!             N2 = -3.*(Mtop - Mr)*(et_nids*et_nids - et_refl*et_refl)
!
!             N3 = 6.*et_refl*(Mtop - Mr)*(et_nids - et_refl)
!
!             N4 = -6.*Mr*(et_nids - et_refl)**3.
!
!           D = D1 + D2 + D3 + D4
!
!             D1 = 2.*(et_nids**3 - et_refl**3.)*(et_nids-et_refl)
!
!             D2 = -6.*(et_refl*et_refl)*(et_nids - et_refl)**2.
!
!             D3 = -3.*(et_nids**2. - et_refl**2.)
!
!             D4 = 6.*et_refl*(et_nids**2. - et_refl**2.)*
!                  (et_nids - et_refl)
!
!         B = (Mtop - Mr - A*(et_nids**2. - et_refl**2.))/(et_nids-et_refl)
!
!         C = Mr - (A*et_refl*et_refl) - (B*et_refl)
!
!         where Mtop is M at et_nids, Mr is M at et_refl, and r is the
!         residual VIL (NIDS VIL - reflectivity field VIL).  So, we need
!         to calculate Mtop, Mr, and r, and then A, B, and C
!
!------------------------------------------------------------------------

          IF (iteration > 2) THEN
            refl_adj(i,j,k_et_refl(i,j)) =                               &
              refl_adj(i,j,k_et_refl(i,j)) - 1.5
          END IF
          Z = 10.**(refl_adj(i,j,k_et_refl(i,j))*0.1)
          Mr = K_const*(Z**foursevenths)

          reflthreshexp = reflthresh*0.1
          Z = 10.**(reflthreshexp)
          Mtop = K_const*(Z**foursevenths)

          r = vil_nids(i,j) - vil_refl(i,j)

          N1 = 6.*r*(et_nids(i,j) - et_refl(i,j))
          N2 = -3.*(Mtop - Mr)*(et_nids(i,j)*et_nids(i,j) -              &
               et_refl(i,j)*et_refl(i,j))
          N3 = 6.*et_refl(i,j)*(Mtop - Mr)*(et_nids(i,j) - et_refl(i,j))
          N4 = -6.*Mr*(et_nids(i,j) - et_refl(i,j))**3.
          N = N1 + N2 + N3 + N4

          D1 = 2.*(et_nids(i,j)**3. -                                    &
               et_refl(i,j)**3.)*(et_nids(i,j)-et_refl(i,j))
          D2 = -6.*(et_refl(i,j)*et_refl(i,j))*                          &
               (et_nids(i,j)-et_refl(i,j))**2.
          D3 = -3.*(et_nids(i,j)*et_nids(i,j) -                          &
               et_refl(i,j)*et_refl(i,j))
          D4 = 6.*et_refl(i,j)*(et_nids(i,j)*et_nids(i,j) -              &
               et_refl(i,j)*et_refl(i,j))*(et_nids(i,j) - et_refl(i,j))
          D = D1 + D2 + D3 + D4

          A = N/D

          B = (Mtop - Mr - A*(et_nids(i,j)*et_nids(i,j) -                &
               et_refl(i,j)*et_refl(i,j)))/(et_nids(i,j) - et_refl(i,j))

          C = Mr - (A*et_refl(i,j)*et_refl(i,j)) - (B*et_refl(i,j))

!------------------------------------------------------------------------
!
!         Now apply the quadratic function to the refl_adj field between
!         k_et_refl and k_ef_nids.
!
!------------------------------------------------------------------------

          DO k = k_et_refl(i,j),k_et_nids(i,j)
            M = (A*zs(i,j,k)*zs(i,j,k)) + (B*zs(i,j,k)) + C

            Z = (M/K_const)**sevenfourths
            refl_adj(i,j,k) = 10.*ALOG10(Z)

          END DO ! k loop

        END IF ! linear or quadratic
        vil_refl(i,j) = vil_88d(nz,refl_adj(i,j,:),zs(i,j,:))
      END DO ! iteration loop

!------------------------------------------------------------------------
!
!     Move on to next column.
!
!------------------------------------------------------------------------

    END DO ! i loop
  END DO ! j loop

!------------------------------------------------------------------------
!
! Calculate new VIL and residual VIL from adjusted reflectivity field.
!
!------------------------------------------------------------------------

  vil_residual(:,:) = -9999.

  maxresidual = 0.
  minresidual = 0.

  DO j = 1,ny-1
    DO i = 1,nx-1
      IF (vil_nids(i,j) /= -9999. .AND. et_nids(i,j) /= -9999. .AND.     &
          vil_refl(i,j) /= -9999. .AND. et_refl(i,j) /= -9999. .AND.     &
          k_et_nids(i,j) /= -9999 .AND. k_et_refl(i,j) /= -9999 .AND.    &
          k_et_nids(i,j) > k_et_refl(i,j)) THEN
        vil_residual(i,j) = vil_nids(i,j) - vil_refl(i,j)
        IF (vil_residual(i,j) < minresidual) THEN
          minresidual = vil_residual(i,j)
        END IF
        IF (vil_residual(i,j) > maxresidual) THEN
          maxresidual = vil_residual(i,j)
        END IF
      END IF
    END DO
  END DO

  WRITE(6,*)'Max Residual VIL: ',maxresidual
  WRITE(6,*)'Min Residual VIL: ',minresidual

!------------------------------------------------------------------------
!
! Overwrite original reflectivity and exit.
!
!------------------------------------------------------------------------

  refl(:,:,:) = refl_adj(:,:,:)

  RETURN
END SUBROUTINE adjust_refl

!########################################################################
!########################################################################
!#########                                                      #########
!#########                SUBROUTINE calc_et_refl               #########
!#########                                                      #########
!#########                     Developed by                     #########
!#########    Center for Analysis and Prediction of Storms      #########
!#########                University of Oklahoma                #########
!#########                                                      #########
!########################################################################
!########################################################################

SUBROUTINE calc_et_refl(nx,ny,nz,refl,reflthresh,zs,et,k_et)

!------------------------------------------------------------------------
!
! PURPOSE:
!
! Determines echo tops for 3D reflectivity field using a dBZe threshold
! specified as an argument.  Also returns the k-levels.
!
!------------------------------------------------------------------------
!
! AUTHOR:
!
! Eric Kemp, 31 August 2001.  Developed for WDT.
!
!------------------------------------------------------------------------
!
! Force explicit declarations.
!
!------------------------------------------------------------------------

  IMPLICIT NONE

!------------------------------------------------------------------------
!
! Declare arguments.
!
!------------------------------------------------------------------------

  INTEGER, INTENT(IN) :: nx,ny,nz        ! Dimensions of grid

  REAL, INTENT(IN) :: refl(nx,ny,nz)     ! 3D reflectivity (dBZ)
  REAL, INTENT(IN) :: reflthresh         ! Reflectivity threshold for
                                         ! echo top (dBZ)
  REAL, INTENT(IN) :: zs(nx,ny,nz)       ! Scalar heights (m MSL)
  REAL, INTENT(INOUT) :: et(nx,ny)       ! Echo top heights (m MSL)
  INTEGER, INTENT(INOUT) :: k_et(nx,ny)  ! k-level of echo top height.

!------------------------------------------------------------------------
!
! Declare internal variables
!
!------------------------------------------------------------------------

  INTEGER :: i,j,k

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  et(:,:) = -9999.
  k_et(:,:) = -9999

  DO j = 1,ny-1
    DO i = 1,nx-1
      DO k = nz-1,2,-1
         IF (refl(i,j,k) > reflthresh) THEN
           et(i,j) = zs(i,j,k)
           k_et(i,j) = k
           EXIT ! Get out of k loop
         END IF
      END DO ! k loop
    END DO ! i loop
  END DO ! j loop
  RETURN
END SUBROUTINE calc_et_refl

!########################################################################
!########################################################################
!#########                                                      #########
!#########                 SUBROUTINE calc_k_et                 #########
!#########                                                      #########
!#########                     Developed by                     #########
!#########    Center for Analysis and Prediction of Storms      #########
!#########                University of Oklahoma                #########
!#########                                                      #########
!########################################################################
!########################################################################

SUBROUTINE calc_k_et(nx,ny,nz,zs,et,k_et)

!------------------------------------------------------------------------
!
! PURPOSE:
!
! Determines the k-level of the (NIDS) echo tops.  Developed for WDT.
!
!------------------------------------------------------------------------
!
! AUTHOR:
!
! Eric Kemp, 31 August 2001
!
!------------------------------------------------------------------------
!
! Force explicit declarations.
!
!------------------------------------------------------------------------

  IMPLICIT NONE

!------------------------------------------------------------------------
!
! Declare arguments.
!
!------------------------------------------------------------------------

  INTEGER, INTENT(IN) :: nx,ny,nz       ! Dimensions of grid
  REAL, INTENT(IN) :: zs(nx,ny,nz)      ! Scalar heights (m MSL)
  REAL, INTENT(IN) :: et(nx,ny)         ! Echo top heights (m MSL)

  INTEGER, INTENT(INOUT) :: k_et(nx,ny) ! k-level of echo tops

!------------------------------------------------------------------------
!
! Declare internal variables
!
!------------------------------------------------------------------------

  INTEGER :: i,j,k

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  k_et(:,:) = -9999

  DO j = 1,ny-1
    DO i = 1,nx-1
       DO k = nz-1,2,-1
         IF (et(i,j) > zs(i,j,k)) THEN
           k_et(i,j) = k
           EXIT ! Get out of k loop
         END IF
       END DO ! k loop
    END DO ! i loop
  END DO ! j loop
  RETURN
END SUBROUTINE calc_k_et

!########################################################################
!########################################################################
!#########                                                      #########
!#########                   FUNCTION vil_88d                   #########
!#########                                                      #########
!#########                     Developed by                     #########
!#########    Center for Analysis and Prediction of Storms      #########
!#########                University of Oklahoma                #########
!#########                                                      #########
!########################################################################
!########################################################################

REAL FUNCTION vil_88d(nz,refl,zs)

!------------------------------------------------------------------------
!
! PURPOSE:
!
! Calculates vertically integrated liquid from 3D reflectivity data
! following the WSR-88D algorithm.
!
!------------------------------------------------------------------------
!
! REFERENCES:
!
! Edwards, R., and R. L. Thompson, 1998:  Nationwide comparisons of hail
!   size with WSR-88D vertically integrated liquid water and derived
!   thermodynamic sounding data.  _Wea. Forecasting_, 12, 277-285.
! Greene, D. R., and R. A. Clark, 1972:  Vertically integrated liquid
!   water -- A new analysis tool.  _Mon. Wea. Rev._, 100, 548-552.
! OSF, Undated Document:  Vertically-integrated liquid water algorithm
!   description.  NX-DR-03-006/25, 7pp.  Available on the Internet at
!   http://120.125.144.136/~swimm/wsr88d/algorithms.htm
! Witt, A., and Coauthors, 1998:  An enhanced hail detection algorithm
!   for the WSR-88D.  _Wea. Forecasting_, 13, 286-303.
!
!------------------------------------------------------------------------
!
! AUTHOR:
!
! Eric Kemp, 7 September 2001.  Developed for WDT.
!
!------------------------------------------------------------------------
!
! Force explicit declarations
!
!------------------------------------------------------------------------

  IMPLICIT NONE

!------------------------------------------------------------------------
!
! Declare arguments.
!
!------------------------------------------------------------------------

  INTEGER, INTENT(IN) :: nz    ! Dimension of column
  REAL, INTENT(IN) :: refl(nz) ! Reflectivity in column (dBZ)
  REAL, INTENT(IN) :: zs(nz)   ! Scalar heights (m MSL)

!------------------------------------------------------------------------
!
! Declare internal parameters.
!
!------------------------------------------------------------------------

  REAL, PARAMETER :: minrefl = 18.5  ! Minimum point reflectivity used
                                     ! in VIL algorithm.
  REAL, PARAMETER :: maxrefl = 55.   ! Maximum point reflectivity used
                                     ! in VIL algorithm.
  REAL, PARAMETER :: maxreflexp = 55.*0.1 ! Exponent used for cases
                                          ! with reflectivity > maxrefl

  REAL, PARAMETER :: K_const = 3.44E-6 ! Conversion factor for Z to M.

  REAL, PARAMETER :: maxvil = 80. ! Maximum VIL permitted (kg m**-2)

  REAL, PARAMETER :: foursevenths = 4./7.

!------------------------------------------------------------------------
!
! Declare internal variables.
!
!------------------------------------------------------------------------

  REAL :: Z, M, dh, VIL
  INTEGER :: k

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  VIL = 0.
  DO k = 2,nz-1

!------------------------------------------------------------------------
!
!   Calculate reflectivity factor Z in mm**6 m**-3.  Logarithmic
!   reflectivities less than 18.5 dBZ are set to zero, while
!   those values above 55 dBZ are truncated.
!
!------------------------------------------------------------------------

     IF (refl(k) < minrefl) THEN
       Z = 0.
     ELSE IF (refl(k) > maxrefl) THEN
       Z = 10.**(maxreflexp)
     ELSE
       Z = 10.**(refl(k)*0.1)
     END IF

!------------------------------------------------------------------------
!
!    Now calculate liquid water content (M) in kg m**-3.
!
!------------------------------------------------------------------------

     M = K_const*(Z**foursevenths)

!------------------------------------------------------------------------
!
!    Calculate local dh.
!
!------------------------------------------------------------------------

     dh = zs(k+1) - zs(k)

!------------------------------------------------------------------------
!
!    Add to VIL.
!
!------------------------------------------------------------------------

     VIL = VIL + (M*dh)
  END DO

!------------------------------------------------------------------------
!
! Return value and exit.
!
!------------------------------------------------------------------------

  vil_88d = MIN(VIL,maxvil)
  RETURN
END FUNCTION vil_88d
