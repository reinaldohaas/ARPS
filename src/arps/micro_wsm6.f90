!
!##################################################################
!##################################################################
!######                                                      ######
!######            SUBROUTINE MICROPH_wsm6_driver            ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE micro_wsm6_driver(mscheme,nx,ny,nz,dtbig1,zp,w,              &
              ptprt,ptbar,pprt,pbar,ppi,                                &
              qv,qvbar,qscalar,                                         &
              raing,prcrate,                                            &
              tk,p,q,qci,qrs,ww,den,delz,rainncv,mpteqnterms,N0x)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Call WRF ice microphysics parameterization scheme WSM6.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yunheng Wang 
!  03/01/2006.
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------

  USE module_mp_wsm6

  IMPLICIT NONE

  INCLUDE 'globcst.inc'

  INTEGER, INTENT(IN)    :: mscheme     ! 0 - 2
             ! 0 - WRF WSM6
             ! 1 - WSM6 with gamma distribution constraint for rain
             ! 2 - WSM6 with diagnostic N0
  INTEGER, INTENT(IN)    :: nx, ny, nz
  REAL,    INTENT(IN)    :: dtbig1
  REAL,    INTENT(IN)    :: zp(nx,ny,nz)
  REAL,    INTENT(IN)    :: w (nx,ny,nz)

  REAL,    INTENT(INOUT) :: ptprt(nx,ny,nz)
  REAL,    INTENT(INOUT) :: pprt (nx,ny,nz)
  REAL,    INTENT(INOUT) :: qv   (nx,ny,nz)
  REAL,    INTENT(IN)    :: ptbar(nx,ny,nz)
  REAL,    INTENT(IN)    :: pbar (nx,ny,nz)
  REAL,    INTENT(IN)    :: qvbar(nx,ny,nz)

  REAL,    INTENT(IN)    :: ppi  (nx,ny,nz)

  REAL,    INTENT(INOUT) :: qscalar(nx,ny,nz,nscalar)

  REAL,    INTENT(INOUT) :: raing  (nx,ny)
  REAL,    INTENT(INOUT) :: prcrate(nx,ny)

  REAL                   :: tk(nx,nz,ny)         ! ikj memory order
  REAL                   :: q (nx,nz,ny)         ! ikj memory order
  REAL                   :: p (nx,nz,ny)         ! ikj memory order

  REAL                   :: qci(nx,nz,2)
  REAL                   :: qrs(nx,nz,3)
  REAL                   :: ww(nx,nz), den(nx,nz), delz(nx,nz)
  REAL                   :: rainncv(nx)
  REAL                   :: mpteqnterms(nx,ny,nz,27)       ! microphysical heating terms
  REAL                   :: N0x(nx,ny,nz,nscalar)  ! intercept parameter
  REAL                   :: mpteqnterms_ikj(nx,nz,27,ny)       ! microphysical heating terms

!-----------------------------------------------------------------------
!
! Misc. local varibles
!
!-----------------------------------------------------------------------

  INTEGER :: i, j, k, nq
  INTEGER :: its,ite,jts,jte,kts,kte
  REAL    :: deltat

  LOGICAL, SAVE :: initialized = .FALSE.

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!-----------------------------------------------------------------------
!
! Since this scheme comtain microphysical ice process, it must at least
! define qc, qr, qi, qs, qg
!
!-----------------------------------------------------------------------

  IF (P_QC < 1 .OR. P_QR < 1 .OR. P_QI < 1 .OR. P_QS < 1 .OR. P_QG < 1) THEN

    WRITE(6,'(2a,/,5(a,I2),/,a)')                                       &
               'No enough microphysical array was defined ',            &
               'inside subroutine micro_wsm6.',                         &
               'P_QC = ',P_QC,' P_QR = ',P_QR,' P_QI = ',P_QI,          &
              ' P_QS = ',P_QS,' P_QG = ',P_QG,                          &
               'Program aborting ...'
    CALL arpsstop('Wrong size for microphysics array, qscalar.',1)

  END IF

!-----------------------------------------------------------------------
!
! Initialize if the first time
!
!-----------------------------------------------------------------------

  IF (mscheme < 0 .OR. mscheme > 2) THEN
    WRITE(6,*) 'Unsupported WSM6 microphysics scheme ',mscheme
    CALL arpsstop('Unsupported WSM6 scheme.',mscheme)
  END IF

  IF (.NOT. initialized) THEN
    CALL wsm6init(mscheme)
    initialized = .TRUE.
  END IF

!-----------------------------------------------------------------------
!
! Time steps assignment
!
!-----------------------------------------------------------------------

  IF( sadvopt /= 4 .and. tintegopt == 1) THEN                  ! Leapfrog scheme
    deltat = 2*dtbig1
  ELSE                                    ! Forward scheme
    deltat = dtbig1
  END IF

!-----------------------------------------------------------------------
!
! Get temperature from potential temperature
! Get total specific humidity (actually it should be mixing ratio,
!      we use them alternately because the difference is small)
! Get total pressure
!
!-----------------------------------------------------------------------

  DO k = 1,nz-1
    DO j = 1,ny-1
      DO i = 1,nx-1
        tk (i,k,j) = (ptprt(i,j,k) + ptbar(i,j,k))*ppi(i,j,k)
        p  (i,k,j) =  pprt (i,j,k) + pbar (i,j,k)
        q  (i,k,j) =  qv   (i,j,k)
      END DO
    END DO
  END DO

!-----------------------------------------------------------------------
!
! Call the microphysics 2D model with each x-z slab
!
!-----------------------------------------------------------------------

  its = 1               ! Should fix these index if OMP is desired.
  ite = nx-1
  jts = 1
  jts = ny-1
  kts = 1
  kte = nz-1

  DO j = 1,ny-1         ! for each x-z slab

!-----------------------------------------------------------------------
!
! Convert to WSM6 arrays from ARPS arrays
!
!-----------------------------------------------------------------------

    DO k = 1,nz-1
      DO i = 1,nx-1
        ww(i,k)    = w(i,j,k)
        qci(i,k,1) = qscalar(i,j,k,P_QC)
        qci(i,k,2) = qscalar(i,j,k,P_QI)
        qrs(i,k,1) = qscalar(i,j,k,P_QR)
        qrs(i,k,2) = qscalar(i,j,k,P_QS)
        qrs(i,k,3) = qscalar(i,j,k,P_QG)
        delz(i,k)  = zp(i,j,k+1)-zp(i,j,k)
        den (i,k)  = p(i,k,j)/(rd*tk(i,k,j))   ! air density at time tfuture
      END DO
    END DO
    
    IF (mscheme == 0) THEN       ! original WRF WSM6 scheme in WRFV2.1.2

      CALL wsm62D_WRF(tk(:,:,j),  q(:,:,j), qci, qrs, ww, den, p(:,:,j),&
                delz, raing(:,j), rainncv, deltat,                      &
                j,                                                      &
                1,   nx-1,  1, ny-1, 1,   nz-1,                         &
                1,   nx,    1, ny,   1,   nz,                           &
                its, ite, jts, jte,  kts, kte,mpteqnterms_ikj(its:ite,kts:kte,:,j))

    ELSE IF (mscheme == 1) THEN  ! Simplified Gamma distribution for rain

      CALL wsm62D_GR (tk(:,:,j),  q(:,:,j), qci, qrs, ww, den, p(:,:,j),&
                delz, raing(:,j), rainncv, deltat,                      &
                j,                                                      &
                1,   nx-1,  1, ny-1, 1,   nz-1,                         &
                1,   nx,    1, ny,   1,   nz,                           &
                its, ite, jts, jte,  kts, kte)

    ELSE IF (mscheme == 2) THEN ! Diagnostic N0

      CALL wsm62D_N0 (tk(:,:,j),  q(:,:,j), qci, qrs, ww, den, p(:,:,j),&
                delz, raing(:,j), rainncv, deltat,                      &
                j,                                                      &
                1,   nx-1,  1, ny-1, 1,   nz-1,                         &
                1,   nx,    1, ny,   1,   nz,                           &
                its, ite, jts, jte,  kts, kte,mpteqnterms_ikj(its:ite,kts:kte,:,j))
    END IF

!-----------------------------------------------------------------------
!
! Convert back to ARPS arrays
!
!-----------------------------------------------------------------------

    DO k = 1,nz-1
      DO i = 1,nx-1
        qscalar(i,j,k,P_QC) = qci(i,k,1)
        qscalar(i,j,k,P_QI) = qci(i,k,2)
        qscalar(i,j,k,P_QR) = qrs(i,k,1)
        qscalar(i,j,k,P_QS) = qrs(i,k,2)
        qscalar(i,j,k,P_QG) = qrs(i,k,3)
        DO nq=1,27
          mpteqnterms(i,j,k,nq) = mpteqnterms_ikj(i,k,nq,j)
        END DO
      END DO
    END DO

    DO i = 1,nx-1
      prcrate(i,j) = rainncv(i)*denr/(1000*dtbig)   ! from mm -> kg m-2 s-1
    END DO

  END DO

!-----------------------------------------------------------------------
!
! Beside hydrometeor arrays, the microphysics scheme changes ARPS 
! potential temperature and mixing ratio only.
!
!-----------------------------------------------------------------------

  N0x = 0.0

  DO k = 1,nz-1
    DO j = 1,ny-1
      DO i = 1,nx-1
        ptprt(i,j,k) = tk(i,k,j)/ppi(i,j,k) - ptbar(i,j,k)
        qv   (i,j,k) =  q(i,k,j)

        ! Added by DTD, 01/07.  Calculate intercept parameter for various species
        ! Cloud (not applicable = 0.0)
        N0x(i,j,k,1) = 0.0 

        ! Rain
        IF (mscheme == 0) THEN  ! Original WRF 
          N0x(i,j,k,2) = 8.0e6 
        ELSE IF (mscheme ==  2) THEN ! Diagnostic N0r (ignore SCG version for now)
          IF(qscalar(i,j,k,P_QR) >= 1.0e-9) THEN
            N0x(i,j,k,2) = 7835.5*1000*((p(i,k,j)/(rd*tk(i,k,j)))*qscalar(i,j,k,P_QR)*1000)**0.681
          ELSE
            N0x(i,j,k,2) = 0.0
          END IF
        END IF

        ! Ice
        N0x(i,j,k,3) = 0.0

        ! Snow (WSM6 scheme uses a temperature-dependent formula for N0s)
        N0x(i,j,k,4) = min(2.0e6*exp(0.12*(273.15-tk(i,k,j))),1.0e11)

        !IF(N0x(i,j,k,2) /= 0.0) THEN
        !  print*,"nonzero N0r before WSM6 exit"
        !END IF

        ! Graupel/hail
        N0x(i,j,k,5) = 4.0e6 
 
      END DO
    END DO
  END DO

  RETURN
END SUBROUTINE micro_wsm6_driver
