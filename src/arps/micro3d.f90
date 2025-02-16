!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE MICROPH                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE microph(nx,ny,nz,mphyflg,dtbig1,                             &
           ptprt,pprt,qv,qscalar,raing,prcrate,                         &
           rhostr,pbar,ptbar,ppi,j3,j3inv,                              &
           rhobar,p,lathvt,tem1,tem2,tem3,tem4,tem5,tem6,tem7,tem8)

!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate and apply the microphysical contributions to the water
!  and temperature fields. The current version implements the Kessler
!  warm rain microphysics parameterization scheme. Ice processes are
!  not included in this version.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  10/10/1991.
!
!  MODIFICATION HISTORY:
!
!  5/02/92 (M. Xue)
!  Added full documentation.
!
!  5/28/92 (K. Brewster)
!  Further facelift.
!
!  11/08/98 (K. Brewster)
!  Added calculation of total p and pi for use in revap and satadj
!
!  12/15/2005 (Y. Wang)
!  Changed microphysics arrays to qscalar. Subroutines autocac, revap,
!  qrfall, satadj in this file were not changed.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!    mphyflg  Flag to indicate whether the microphysic processes will be
!             computed and the adjustments applied.
!             = 0, Do not compute;
!             = 1, Compute.
!             Hydrometeor sedimentation will still be calculated.
!
!    dtbig1   Big time step size (s)
!    ptprt    Perturbation potential temperature at all time levels (K)
!    pprt     Perturbation pressure at all time levels (Pascal)
!    qv       Water vapor specific humidity at all time levels (kg/kg)
!    qscalar  A larger array contains former 4d arrays, qc, qr, qi, qs, qh
!    raing    Accumulated grid-scale rainfall (mm)
!
!    rhostr   Base state air density times j3 (kg/m**3)
!    pbar     Base state pressure (Pascal)
!    ptbar    Base state potential temperature (K)
!    ppi      Exner function at tfuture
!    j3       Coordinate transformation Jacobian  d(zp)/dz
!
!  OUTPUT:
!
!    ptprt    Perturbation potential temperature at time tfuture (K)
!    pprt     Perturbation pressure at time tfuture (Pascal)
!    qv       Water vapor specific humidity at time tfuture (kg/kg)
!    qc       Cloud water mixing ratio at time tfuture (kg/kg)
!    qr       Rainwater mixing ratio at time tfuture (kg/kg)
!    qi       Cloud ice mixing ratio at time tfuture (kg/kg)
!    qs       Snow mixing ratio at time tfuture (kg/kg)
!    qh       Hail mixing ratio at time tfuture (kg/kg)
!    raing    Accumulated grid-scale rainfall (mm)
!    prcrate  Precipitation rate (kg/(m**2*s))
!
!  WORK ARRAYS:
!
!    rhobar   Base state air density (kg/m**3)
!    lathvt   Temperature-dependent latent heat of evaporation
!    p        Total pressure at tfuture
!    tem1     Temporary work array.
!    tem2     Temporary work array.
!    tem3     Temporary work array.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'phycst.inc'
  INCLUDE 'globcst.inc'
  INCLUDE 'timelvls.inc'
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!

  INTEGER :: nx,ny,nz          ! Number of grid points in 3 directions
  INTEGER :: mphyflg           ! microphysic flag

  REAL :: dtbig1               ! Big time step size (s)
  REAL :: ptprt (nx,ny,nz,nt)  ! Perturbation potential temperature (K)
  REAL :: pprt  (nx,ny,nz,nt)  ! Perturbation pressure (Pascal)

  REAL :: qv    (nx,ny,nz,nt)  ! Water vapor specific humidity (kg/kg)
  REAL :: qscalar(nx,ny,nz,nt,nscalar)
  REAL :: raing (nx,ny)        ! Accumulated grid-scale rainfall (mm)
  REAL :: prcrate(nx,ny)       ! Precipitation rate (kg/(m**2*s))


  REAL :: rhostr(nx,ny,nz)     ! Base state air density times j3 (kg/m**3)
  REAL :: pbar  (nx,ny,nz)     ! Base state pressure (Pascal).
  REAL :: ptbar (nx,ny,nz)     ! Base state potential temperature (K)
  REAL :: j3    (nx,ny,nz)     ! Coordinate transformation Jacobian defined as
                               ! d( zp )/d( z ).
  REAL :: ppi   (nx,ny,nz)     ! Exner function at tfuture
  REAL :: j3inv (nx,ny,nz)     ! Coordinate transformation Jacobian defined as
                               ! d( zp )/d( z ).
!
!-----------------------------------------------------------------------
!
!  Temporary arrays
!
!-----------------------------------------------------------------------
!
  REAL :: p     (nx,ny,nz)    ! Total pressure at tfuture (Pa)
  REAL :: rhobar(nx,ny,nz)    ! Base state air density (kg/m**3)
  REAL :: lathvt(nx,ny,nz)    ! Temperature dependent latent heat
                              ! of evaporation. Local array

  REAL :: tem1(nx,ny,nz)      ! Temporary work array
  REAL :: tem2(nx,ny,nz)      ! Temporary work array
  REAL :: tem3(nx,ny,nz)      ! Temporary work array
  REAL :: tem4(nx,ny,nz)      ! Temporary work array
  REAL :: tem5(nx,ny,nz)      ! Temporary work array
  REAL :: tem6(nx,ny,nz)      ! Temporary work array
  REAL :: tem7(nx,ny,nz)      ! Temporary work array
  REAL :: tem8(nx,ny,nz)      ! Temporary work array
!
!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------
!
!
  INTEGER :: i,j,k,tindex
  REAL    :: tk
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!-----------------------------------------------------------------------
!
! Since only warm rain microphysics is computed, P_QC & P_QR must be 
! defined.
!
!-----------------------------------------------------------------------

  IF (P_QC < 1 .OR. P_QR < 1) THEN
    WRITE(6,'(a,/,2(a,I2),/,a)')                                        &
            'Either P_QC or P_QR was not defined inside subroutine microph.', &
            'P_QC = ',P_QC,' P_QR = ',P_QR,                             &
            'Program aborting ...'
    CALL arpsstop('Wrong size for microphysics array, qscalar.',1)
  END IF
!
!-----------------------------------------------------------------------
!
!  Since rhobar is used most often in the microphysics parameterizations,
!  we calculate and save rhobar for later use.
!
!-----------------------------------------------------------------------
!
  DO k = 1,nz-1
    DO j = 1,ny-1
      DO i = 1,nx-1
        rhobar(i,j,k) = rhostr(i,j,k)*j3inv(i,j,k)
      END DO
    END DO
  END DO

  IF(mphyflg == 1) THEN       ! Only compute when mphyflg = 1 
!
!-----------------------------------------------------------------------
!
!  To remove negative mixing ratios, which result from computational
!  inaccuracies in the advection process, we set all negative mixing
!  ratios to zero. This is an artificial adjustment and, as a result,
!  total water will not be conserved. The adjustment can be averted
!  by enhancing the numerical accuracy of subsequent model versions.
!
!-----------------------------------------------------------------------
!
    DO k = 1,nz
      DO j = 1,ny
        DO i = 1,nx
          qv(i,j,k,tfuture) = MAX(0.0,qv(i,j,k,tfuture))
          qscalar(i,j,k,tfuture,P_QC) = MAX(0.0,qscalar(i,j,k,tfuture,P_QC))
          qscalar(i,j,k,tfuture,P_QR) = MAX(0.0,qscalar(i,j,k,tfuture,P_QR))
        END DO
      END DO
    END DO
!
!-----------------------------------------------------------------------
!
!  Use the look up table to calculate the temperature dependent
!  latent heat of evaporation lathvt.
!  Also calculate Exner function.
!
!-----------------------------------------------------------------------
!
    DO k=1,nz-1
      DO j=1,ny-1
        DO i=1,nx-1

          p(i,j,k) = pprt(i,j,k,tfuture)+pbar(i,j,k)
          tk = (ptprt(i,j,k,tfuture)+ptbar(i,j,k)) * ppi(i,j,k)
          tindex = MAX( MIN(151, nint(tk-172.15)), 1)
          lathvt(i,j,k) = latheatv(tindex)

        END DO
      END DO
    END DO

  END IF        ! mphyflg == 1

  IF( mphyopt /= 0 .OR. cnvctopt == 2) THEN

    IF(mphyflg == 1) THEN

!-----------------------------------------------------------------------
!
!  Calculate the autoconversion of cloud water to rainwater
!  and the accretion (collection) of cloud droplets by rain drops.
!
!-----------------------------------------------------------------------
!
      CALL autocac(nx,ny,nz,dtbig1,qscalar(:,:,:,:,P_QC),qscalar(:,:,:,:,P_QR))
!
!
!-----------------------------------------------------------------------
!
!  Rainwater evaporates when the air is subsaturated. Under these
!  conditions, water vapor specific humidity increases at the expense
!  of rainwater and the air is subsequently cooled due to evaporation.
!
!  The evaporation rate is subject to certain specified limits.
!  Rainwater, water vapor and temperature fields are then adjusted.
!
!-----------------------------------------------------------------------
!

      CALL revap(nx,ny,nz, dtbig1,ptprt,qv,qscalar(:,:,:,:,P_QR),       &
                 rhobar,ptbar,p,ppi,lathvt, tem1,tem2)


    END IF     ! mphyflg == 1
!
!-----------------------------------------------------------------------
!
!                        RAINWATER SEDIMENTATION
!
!  Calculate the rain fallout rate and apply it to the rain
!  water field.
!
!  Since the rainwater fallout is an advective process,
!  negative values of rainwater can be generated by the finite
!  difference scheme. If this occurs, we set the negative water
!  content to zero.
!
!-----------------------------------------------------------------------
!
    CALL qrfall (nx,ny,nz,dtbig1,                                       &
                 qscalar(:,:,:,:,P_QR),rhobar,j3,j3inv,raing,prcrate,   &
                 tem1,tem2,tem3,tem4,tem5,tem6,tem7,tem8)

    DO k = 1,nz
      DO j = 1,ny
        DO i = 1,nx
          qscalar(i,j,k,tfuture,P_QR) = MAX(0.0,qscalar(i,j,k,tfuture,P_QR))
        END DO
      END DO
    END DO
!
!-----------------------------------------------------------------------
!
!                     SATURATION ADJUSTMENT
!             (CONDENSATION AND EVAPORATION BETWEEN QV AND QC)
!
!  Under the conditions of supersaturation, condensation occurs
!  which converts water vapor to cloud water.
!  .
!
!  If the air is subsaturated, cloud water evaporates and becomes
!  water vapor until either the air is saturated or the cloud water
!  has completely evaporated.
!
!  During these processes, the air temperature changes.
!
!  Adjust the potential temperature, water vapor and cloud water
!  accordingly.
!
!-----------------------------------------------------------------------
!
!
  END IF

  IF(mphyflg == 1) THEN

    CALL satadj(nx,ny,nz, ptprt,qv,qscalar(:,:,:,:,P_QC),ptbar,         &
                p,ppi,lathvt,tem1,tem2,tem3)

    DO k = 1,nz
      DO j = 1,ny
        DO i = 1,nx
          qv(i,j,k,tfuture)           = MAX(0.0,qv(i,j,k,tfuture))
          qscalar(i,j,k,tfuture,P_QC) = MAX(0.0,qscalar(i,j,k,tfuture,P_QC))
        END DO
      END DO
    END DO

  END IF    ! mphyflg == 1

  RETURN
END SUBROUTINE microph
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE AUTOCAC                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE autocac(nx,ny,nz, dtbig1, qc,qr)

!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate the autoconversion of cloud water to rainwater and the
!  accretion (collection) of cloud droplets by rain drops.
!  First calculate the conversion and accretion rates, then impose
!  limits on the amount of conversion and accretion. Adjust the
!  cloud and rainwater fields accordingly.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  10/10/1991.
!
!  MODIFICATION HISTORY:
!
!  5/02/92 (M. Xue)
!  Added full documentation.
!
!  5/28/92 (K. Brewster)
!  Further facelift.
!
!  6/10/94 (M Xue & A. Sathye)
!  Power calculation for the accretion rate is replaced by a
!  lookup table function.
!
!  07/10/97 (Fanyou Kong - CMRP)
!  Include MPDCD advection option (sadvopt = 5)
!
!  11/18/98 (Keith Brewster)
!  Changed pibar to ppi (full pi).
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!    qc       Cloud water mixing ratio at all time levels (kg/kg)
!    qr       Rainwater mixing ratio at all time levels (kg/kg)
!
!  OUTPUT:
!
!    qc       Cloud water mixing ratio at time tfuture (kg/kg)
!    qr       Rainwater mixing ratio at time tfuture (kg/kg)
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
!
  INTEGER :: nt           ! The no. of t-levels of t-dependent arrays.
  INTEGER :: tpast        ! Index of time level for the past time.
  INTEGER :: tpresent     ! Index of time level for the present time.
  INTEGER :: tfuture      ! Index of time level for the future time.

  PARAMETER (nt=3, tpast=1, tpresent=2, tfuture=3)
!
  INTEGER :: nx,ny,nz          ! Number of grid points in 3 directions

  REAL :: dtbig1               ! Big time step size (s)

  REAL :: qc    (nx,ny,nz,nt)  ! Cloud water mixing ratio (kg/kg)
  REAL :: qr    (nx,ny,nz,nt)  ! Rainwater mixing ratio (kg/kg)
!
!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------
!
!
  INTEGER :: i,j,k,lvlq
  REAL :: deltat,qcplus,qrplus,ar,arcrdt,cr
  REAL :: temp, f1, f2, rstep
  INTEGER :: INDEX
!
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!

  INCLUDE 'phycst.inc'
  INCLUDE 'globcst.inc'

!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  lvlq = tfuture
!
!-----------------------------------------------------------------------
!
!  When water and temperture fields are advected using leapfrog scheme:
!  (if using forward scheme, deltat = dtbig1.)
!
!-----------------------------------------------------------------------
!
!C    IF( sadvopt.ge.1.and.sadvopt.le.3.or.sadvopt.eq.5) THEN ! Leapfrog scheme
  IF( sadvopt /= 4 .and. tintegopt == 1) THEN                  ! Leapfrog scheme
    deltat = 2*dtbig1
  ELSE                                    ! Forward scheme
    deltat = dtbig1
  END IF
!
!-----------------------------------------------------------------------
!
!  Autoconversion rate of cloud water to rainwater:
!
!      AR = autort * (QC - autotr)
!
!  autort : Autoconversion rate parameter (defined in globcst.inc)
!  autotr:  Threshold value for the autoconversion process to begin.
!
!
!  Accretion (collection) rate of cloud water by rainwater
!
!      CR = accrrt * QC * (QR ** 0.875)
!
!  accrrt:  Accretion rate multiplier (defined in globcst.inc)
!
!-----------------------------------------------------------------------

  rstep = 1.0/50.0E-6

  DO k=1,nz-1
    DO j=1,ny-1
      DO i=1,nx-1

!
!-----------------------------------------------------------------------
!
!  Take only the positive values of qc and qr.
!
!-----------------------------------------------------------------------
!
        qcplus = MAX(0.0, qc(i,j,k,lvlq) )
        qrplus = MAX(0.0, qr(i,j,k,lvlq) )
!
!-----------------------------------------------------------------------
!
!  The autoconversion rate:
!
!-----------------------------------------------------------------------
!
        ar = autort *( qcplus - autotr)
!
!-----------------------------------------------------------------------
!
!  Autoconversion only occurs when qc > autotr
!  i.e., AR calculated above should be >0.0
!
!-----------------------------------------------------------------------
!
        ar = MAX(0.0, ar)

        arcrdt = MIN( ar*deltat, qcplus )
        qc(i,j,k,lvlq) = qc(i,j,k,lvlq) - arcrdt
        qr(i,j,k,lvlq) = qr(i,j,k,lvlq) + arcrdt
        qcplus = MAX(0.0, qc(i,j,k,lvlq) )
        qrplus = MAX(0.0, qr(i,j,k,lvlq) )
!
!-----------------------------------------------------------------------
!
!  Accretion rate:
!
!      CR = accrrt * QC * (QR ** 0.875)
!
!  Here a linear interpolation from a lookup table data is used to
!  evaluate power function qr ** 0.875
!  for which it is assumed that 0.0=< qr =< 0.05.
!
!-----------------------------------------------------------------------
!

        temp = MIN(0.05, qrplus)*rstep
        INDEX = INT(temp)

        f1 = pwr875(INDEX)
        f2 = pwr875(INDEX + 1)

        cr = accrrt * qcplus * (f1 + ((f2 - f1) * (temp - INDEX) ))
!
!-----------------------------------------------------------------------
!
!  Calculate the total conversion from cloud water to rainwater due
!  to autoconversion and accretion. This should not exceed the
!  amount of cloud water present:
!
!-----------------------------------------------------------------------
!
!    arcrdt = min( (ar + cr)*deltat, qcplus )

        arcrdt = MIN( cr*deltat, qcplus )
!
!-----------------------------------------------------------------------
!
!  Update qc and qr to account for loss in cloud water
!  and gain in rainwater.
!
!-----------------------------------------------------------------------

        qc(i,j,k,lvlq) = qc(i,j,k,lvlq) - arcrdt

        qr(i,j,k,lvlq) = qr(i,j,k,lvlq) + arcrdt


      END DO
    END DO
  END DO

  RETURN
END SUBROUTINE autocac
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE REVAP                      ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE revap(nx,ny,nz, dtbig1,ptprt,qv,qr,                          &
           rhobar,ptbar,p,ppi, lathvt, qvs, tem1)

!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate the rainwater evaporation rate. Apply these changes
!  to QR and QV and then adjust the potential temperature PTPRT
!  due to the evaporative cooling.
!
!  Rainwater evaporates when the air is subsaturated. Under these
!  conditions, water vapor specific humidity increases at the expense
!  of rainwater and the air is subsequently cooled due to evaporation.
!
!  The evaporation rate is subject to certain specified limits.
!  Rainwater, water vapor and temperature fields are then adjusted.
!
!-----------------------------------------------------------------------

!
!  AUTHOR: Ming Xue
!  2/29/1992.
!
!  MODIFICATION HISTORY:
!
!  5/02/92 (M. Xue)
!  Added full documentation.
!
!  5/28/92 (K. Brewster)
!  Further facelift.
!
!  6/10/94 (M Xue & A. Sathye)
!  Power calculations are replaced by lookup table functions.
!
!  11/08/98 (K. Brewster)
!  Changed calculation of qvs to be based on total pi and p, not pbar
!
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!    ptprt    Perturbation potential temperature at all time levels (K)
!    qv       Water vapor specific humidity at all time levels (kg/kg)
!    qr       Rainwater mixing ratio at all time levels (kg/kg)
!
!    ptbar    Base state potential temperature (K)
!    p        Total pressure at tfuture
!    ppi      Exner function of pressure
!    lathvt   Temperature-dependent latent heat of evaporation
!
!  OUTPUT:
!
!    ptprt    Perturbation potential temperature at time tfuture (K)
!    qv       Water vapor specific humidity at time tfuture (kg/kg)
!    qr       Rainwater mixing ratio at time tfuture (kg/kg)
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
!
  INTEGER :: nt           ! The no. of t-levels of t-dependent arrays.
  INTEGER :: tpast        ! Index of time level for the past time.
  INTEGER :: tpresent     ! Index of time level for the present time.
  INTEGER :: tfuture      ! Index of time level for the future time.

  PARAMETER (nt=3, tpast=1, tpresent=2, tfuture=3)
!
  INTEGER :: nx,ny,nz          ! Number of grid points in 3 directions

  REAL :: dtbig1               ! Big time step size (s)

  REAL :: ptprt (nx,ny,nz,nt)  ! Perturbation potential temperature (K)
  REAL :: qv    (nx,ny,nz,nt)  ! Water vapor specific humidity (kg/kg)
  REAL :: qr    (nx,ny,nz,nt)  ! Rainwater mixing ratio (kg/kg)

  REAL :: ptbar (nx,ny,nz)     ! Base state potential temperature (K)
  REAL :: rhobar(nx,ny,nz)     ! Base state air density (kg/m**3)
!
  REAL :: p     (nx,ny,nz)     ! Total pressure at tfuture.
  REAL :: ppi   (nx,ny,nz)     ! Exner function of pressure.
  REAL :: lathvt(nx,ny,nz)     ! Temperature dependent latent heat

!-----------------------------------------------------------------------
!
!  Temporary arrays
!
!-----------------------------------------------------------------------
!
  REAL :: qvs   (nx,ny,nz)      ! Saturation specific humidity
                                ! with respect to liquid water.
  REAL :: tem1  (nx,ny,nz)      ! Temperature (K)
!
!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k,lvlq,lvlt
  REAL :: deltat,qrplus,qvplus,c,er,erdt
  REAL :: interp, temp, f1, f2, rstep
  INTEGER :: INDEX
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'phycst.inc'
  INCLUDE 'globcst.inc'
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

  lvlq = tfuture
  lvlt = tfuture
!
!-----------------------------------------------------------------------
!
!  The time increment when water and temperature fields are advected
!  using the leapfrog scheme (if using forward scheme, deltat =
!  dtbig1).
!
!-----------------------------------------------------------------------
!
!C    IF( sadvopt.ge.1.and.sadvopt.le.3) THEN ! Leapfrog scheme
  IF( sadvopt /= 4 .and. tintegopt == 1) THEN                  ! Leapfrog scheme
    deltat = 2*dtbig1
  ELSE                                    ! Forward scheme
    deltat = dtbig1
  END IF
!
!-----------------------------------------------------------------------
!
!   Calculate qvs, the saturation specific humidity, using
!   enhanced Teten's formula
!
!-----------------------------------------------------------------------

  DO k = 1,nz-1
    DO j = 1,ny-1
      DO i = 1,nx-1
        tem1(i,j,k) = (ptprt(i,j,k,lvlt)+ptbar(i,j,k)) * ppi(i,j,k)
      END DO
    END DO
  END DO

  CALL getqvs(nx,ny,nz, 1,nx-1,1,ny-1,1,nz-1, p,tem1, qvs)


! write (*,*) "ZUWEN subsatopt/rhsat", subsatopt, rhsat 
  IF (subsatopt /= 0) THEN 

    DO k = 1,nz-1
      DO j = 1,ny-1
        DO i = 1,nx-1
          qvs(i,j,k) = rhsat*qvs(i,j,k)   ! adjust qvs for sub-saturation
        END DO
      END DO
    END DO

  END IF 

!
!-----------------------------------------------------------------------
!
!  Calculate the amount of rainwater evaporation and adjust
!  the rainwater, water vapor and temperature fields.
!
!-----------------------------------------------------------------------
!
  rstep = 1.0/50.0E-6

  DO k=1,nz-1
    DO j=1,ny-1
      DO i=1,nx-1

        qrplus = MAX(0.0, qr(i,j,k,lvlq) )
        qvplus = MAX(0.0, qv(i,j,k,lvlq) )
!
!-----------------------------------------------------------------------
!
!  The formula for the calculation of the coefficient for the
!  rainwater evaporation rate
!
!    c = 1.6 + 30.3922 * ( rhobar(i,j,k)*qrplus )**0.2046.
!
!  Here a linear interpolation from a lookup table data is used to
!  evaluate power function (rhobar*qr) ** 0.2046
!  for which it is assumed that 0.0=< rhobar*qr =< 0.05.
!
!-----------------------------------------------------------------------
!
        temp = MIN (0.05, qrplus * rhobar(i,j,k)) * rstep
        INDEX = INT(temp)

        f1 = pwr2046(INDEX)
        f2 = pwr2046(INDEX + 1)

        interp= f1 + ((f2 - f1) * (temp - INDEX) )

        c = 1.6 + 30.3922 * interp
!
!-----------------------------------------------------------------------
!
!  The formula for calculation of the rainwater evaporation rate
!
!     (rhobar(i,j,k)*qrplus)**0.525
!
!  Here a linear interpolation from a lookup table data is used to
!  evaluate power function (rhobar*qr) ** 0.525
!  for which it is assumed that 0.0=< rhobar*qr =< 0.05.
!
!-----------------------------------------------------------------------
!

        f1 = pwr525(INDEX)
        f2 = pwr525(INDEX + 1)

        interp= f1 + ((f2 - f1) * (temp - INDEX) )

        er = c * (1.0- qvplus/qvs(i,j,k) )*                             &
                interp/( (2.03E4 + 9.584E6/(p(i,j,k)*qvs(i,j,k)))       &
                *rhobar(i,j,k) )
!
!-----------------------------------------------------------------------
!
!  The amount of rain evaporation is limited by the available
!  rainwater:
!
!-----------------------------------------------------------------------
!
        erdt = er * deltat
        erdt = MIN( qrplus, MAX(0.0, erdt) )
!
!-----------------------------------------------------------------------
!
!  Adjust qr, qv and ptprt accordingly.
!
!-----------------------------------------------------------------------
!
        qr(i,j,k,lvlq) = qr(i,j,k,lvlq) - erdt

        qv(i,j,k,lvlq) = qv(i,j,k,lvlq) + erdt

        ptprt(i,j,k,lvlt) = ptprt(i,j,k,lvlt)                           &
              - lathvt(i,j,k)*erdt/( ppi(i,j,k)*cp )

      END DO
    END DO
  END DO


  RETURN
END SUBROUTINE revap
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE QRFALL                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE qrfall(nx,ny,nz,dtbig1,                                      &
           qr,rhobar,j3,j3inv,raing, prcrate,                           &
           draing,vtr,tem1,tem2,tem3,tem4,tem5,tem6)

!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate the fall-out rate of rainwater and apply this effect to
!  the rainwater field, qr. Since the leapfrog-centered scheme is
!  used, the rainwater flux divergence (i.e., the advection by the
!  terminal fallout speed of rain) is calculated at time ""tpresent".
!  This subroutine is called by the warm rain package.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  10/10/1991.
!
!  MODIFICATION HISTORY:
!
!  5/02/92 (M. Xue)
!  Added full documentation.
!
!  5/28/92 (K. Brewster)
!  Further facelift.
!
!  6/10/94 (M Xue & A. Sathye)
!  Power calculations are replaced by lookup table functions.
!
!  5/8/1995 (M. Xue)
!  Time-splitting upstream advection used for rainwater sedimentation
!
!  10/20/1996 (M. Xue)
!  Created a general routine for hydrometeor sedimentation called
!  QFALLOUT, which is now also called by ice microphysics routine.
!  Precipitation rate array prcrate is now correctly calculated
!  for split-step integration.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!    dtbig1   Big time step size (s)
!    qr       Rainwater mixing ratio at all time levels (kg/kg)
!
!    rhobar   Base state air density (kg/m**3)
!    j3       Coordinate transformation Jacobian  d(zp)/dz
!    raing    Accumulated grid-scale rainfall (mm)
!
!  OUTPUT:
!
!    qr       Rainwater mixing ratio at time tfuture (kg/kg)
!    raing    Accumulated grid-scale rainfall (mm)
!    prcrate  Precipitation rate
!
!  WORK ARRAYS:
!
!    draing   Gridscale rainfall (mm)
!    vtr      terminal velocity
!    tem1     Work array
!
!   (These arrays are defined and used locally (i.e. inside this
!    subroutine), they may also be passed into routines called by
!    this one. Exiting the call to this subroutine, these temporary
!    work arrays may be used for other purposes and therefore their
!    contents may be overwritten. Please examine the usage of work
!    arrays before you alter the code.)
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
!
  INTEGER :: nt           ! The no. of t-levels of t-dependent arrays.
  INTEGER :: tpast        ! Index of time level for the past time.
  INTEGER :: tpresent     ! Index of time level for the present time.
  INTEGER :: tfuture      ! Index of time level for the future time.

  PARAMETER (nt=3, tpast=1, tpresent=2, tfuture=3)
!
  INTEGER :: nx,ny,nz          ! Number of grid points in 3 directions

  REAL :: dtbig1               ! Big time step size (s)
  REAL :: qr    (nx,ny,nz,nt)  ! Rainwater mixing ratio (kg/kg)
  REAL :: rhobar(nx,ny,nz)     ! Base state air density (kg/m**3)
  REAL :: j3    (nx,ny,nz)     ! Coordinate transformation Jacobian
                               ! defined as d( zp )/d( z ).
  REAL :: j3inv (nx,ny,nz)     ! Coordinate transformation Jacobian
                               ! defined as d( zp )/d( z ).
  REAL :: raing (nx,ny)        ! Accumulated grid-scale rainfall (mm)
  REAL :: prcrate(nx,ny)       ! Precipitation rate (kg/(m**2*s))
!
!-----------------------------------------------------------------------
!
!  Temporary arrays
!
!-----------------------------------------------------------------------
!
  REAL :: draing(nx,ny)       ! Gridscale rainfall (mm)
  REAL :: vtr   (nx,ny,nz)    ! terminal velocify

  REAL :: tem1  (nx,ny,nz)    ! work array
  REAL :: tem2  (nx,ny,nz)    ! work array
  REAL :: tem3  (nx,ny,nz)    ! work array
  REAL :: tem4  (nx,ny,nz)    ! work array
  REAL :: tem5  (nx,ny,nz)    ! work array
  REAL :: tem6  (nx,ny,nz)    ! work array

!
!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k,lvlq
  REAL :: temp, interp, f1, f2, rstep
  INTEGER :: INDEX
  REAL :: denwater,tem
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'phycst.inc'

!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  lvlq = tfuture
!
!-----------------------------------------------------------------------
!
!  The formula for terminal fallout speed of rainwater:
!
!  Vtr = 36.34 * ((0.001 * rhobar * qr) ** 0.1364) * (rho0 /rhobar)**0.5
!  Vtr is positive downwards.
!
!  A linear interpolation from a lookup table data is used to
!  evaluate power function (0.001 * rhobar * qr) ** 0.1364
!  for which it is assumed that 0.0=< rhobar * qr =< 0.05.
!
!
!-----------------------------------------------------------------------
!
  rstep = 1.0/50.0E-9

  DO k=1,nz-1
    DO j=1,ny-1
      DO i=1,nx-1

        temp = MIN (0.00005, MAX(0.0,0.001 * rhobar(i,j,k)              &
              * qr(i,j,k,lvlq-1)))*rstep
        INDEX = INT(temp)

        f1 = pwr1364(INDEX)
        f2 = pwr1364(INDEX+1)

        interp= f1 + ((f2 - f1) * (temp - INDEX) )

        vtr(i,j,k) = 36.34 *  interp * SQRT(rho0/rhobar(i,j,k))

      END DO
    END DO
  END DO
!
!-----------------------------------------------------------------------
!
!  Calculate the rainwater fallout term and the precipitation rate.
!
!-----------------------------------------------------------------------
!

  CALL qfallout(nx,ny,nz,dtbig1,qr,rhobar,j3,j3inv,draing,vtr,          &
                tem1,tem2,tem3,tem4,tem5,tem6)

  denwater = 1000.0   ! Density of liquid water (kg/m**3)
  tem = denwater/(1000.0*dtbig)

  DO j=1,ny-1
    DO i=1,nx-1
      raing(i,j)=raing(i,j)+draing(i,j)
      prcrate(i,j) = draing(i,j)*tem
    END DO
  END DO

  RETURN
END SUBROUTINE qrfall
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE QFALLOUT                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE qfallout(nx,ny,nz,dtbig1,q,rhobar,j3,j3inv,draing,vtr,       &
                    rhovq,tem1,tem2,tem3,tem4,tem5)

!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate the fall-out of hydrometeor given it's fallout velocity.
!  Upstream advection with split time steps is used here.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  10/10/1991.
!
!  MODIFICATION HISTORY:
!
!  5/8/1995 (M. Xue)
!  Time-splitting upstream advection used for rainwater sedimentation
!
!  12/6/1996 (Yuhe Liu)
!  Moved the initialization of draing to the beginning of executable
!  statements. This solved the problem that the returned draing
!  became undefined when no rain sedimentation needed.
!
!  03/20/2002 (Dan Weber, M. Xue & X. Jin)
!  added an option for explcit or implicit fall velocity scheme 
!  
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!    dtbig1   Big time step size (s)
!    q        Hydrometeor mixing ratio at all time levels (kg/kg)
!
!    rhobar   Base state air density (kg/m**3)
!    j3       Coordinate transformation Jacobian  d(zp)/dz
!    vtr      terminal velocity (m/s)
!
!  OUTPUT:
!
!    q        Rainwater mixing ratio at time tfuture (kg/kg)
!    draing   Grid-scale rainfall (mm)
!
!  WORK ARRAYS:
!
!    rhovq    rhobar*terminal velocity*q
!    tem1     Temporary Work Array
!    tem2     Temporary Work Array
!    tem3     Temporary Work Array
!    tem4     Temporary Work Array
!    tem5     Temporary Work Array
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
!
  INTEGER :: nt           ! The no. of t-levels of t-dependent arrays.
  INTEGER :: tpast        ! Index of time level for the past time.
  INTEGER :: tpresent     ! Index of time level for the present time.
  INTEGER :: tfuture      ! Index of time level for the future time.

  PARAMETER (nt=3, tpast=1, tpresent=2, tfuture=3)
!
  INTEGER :: nx,ny,nz          ! Number of grid points in 3 directions

  REAL :: dtbig1               ! Big time step size (s)
  REAL :: q     (nx,ny,nz,nt)  ! Rainwater mixing ratio (kg/kg)
  REAL :: rhobar(nx,ny,nz)     ! Base state air density (kg/m**3)
  REAL :: j3    (nx,ny,nz)     ! Coordinate transformation Jacobian
                               ! defined as d( zp )/d( z ).
  REAL :: j3inv (nx,ny,nz)     ! Coordinate transformation Jacobian
                               ! defined as d( zp )/d( z ).
  REAL :: draing (nx,ny)       ! Grid-scale rainfall (mm)
  REAL :: vtr   (nx,ny,nz)     ! terminal velocify
!
!-----------------------------------------------------------------------
!
!  Temporary arrays
!
!-----------------------------------------------------------------------
!
  REAL :: rhovq (nx,ny,nz)    ! rhobar*terminal velocity*q

  REAL :: tem1(nx,ny,nz)      ! temporary work array
  REAL :: tem2(nx,ny,nz)      ! temporary work array
  REAL :: tem3(nx,ny,nz)      ! temporary work array
  REAL :: tem4(nx,ny,nz)      ! temporary work array
  REAL :: tem5(nx,ny,nz)      ! temporary work array 

!
!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k,nrnstp,irnstp
  REAL    :: denwater,tem, vtrdzmax,dtrnf,dtrnf0
  REAL    :: deltat
  REAL    :: temmin
  INTEGER :: memory_status
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'grid.inc'          ! Grid & map parameters.
  INCLUDE 'phycst.inc'
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  CALL acct_interrupt(qfall_acct)

  vtrdzmax = 0.0
  DO k=1,nz-1
    DO j=1,ny-1
      DO i=1,nx-1
        vtrdzmax=MAX(vtrdzmax,vtr(i,j,k)*(j3inv(i,j,k)*dzinv))
      END DO
    END DO
  END DO

! for bit-for-bit MP accuracy, slows execution down significantly:
  temmin = 0
  CALL mpmax0(vtrdzmax,temmin)

  IF( vtrdzmax == 0.0 ) RETURN  ! vtr=0, no rain sedimentation.

  IF( impfallopt == 0) THEN
    DO j=1,ny-1
      DO i=1,nx-1
        draing(i,j)=0.0
      END DO
    END DO

    IF( sadvopt /= 4 .and. tintegopt == 1) THEN                  ! Leapfrog scheme
      deltat = 2*dtbig1
    ELSE                                    ! Forward scheme
      deltat = dtbig1
    END IF

    dtrnf = MIN(deltat, 0.5/MAX(1.0E-20,vtrdzmax))
                         ! 0.5 is the max Courant number

    dtrnf0 = dtrnf
    nrnstp = nint ( deltat/dtrnf )
    IF(nrnstp == 0) nrnstp=1

    dtrnf  = deltat/nrnstp
    IF (dtrnf > dtrnf0) THEN
      nrnstp = nrnstp + 1
      dtrnf  = deltat/nrnstp
    END IF

    DO irnstp = 1, nrnstp
!
!-----------------------------------------------------------------------
!
!  The density weighted rainwater fallout flux:
!  rhovq = rhobar * vtr * q
!
!-----------------------------------------------------------------------
!
      DO k=1,nz-1
        DO j=1,ny-1
          DO i=1,nx-1
            rhovq(i,j,k) = rhobar(i,j,k)*vtr(i,j,k)*q(i,j,k,tfuture)
          END DO
        END DO
      END DO
!
!-----------------------------------------------------------------------
!
!  Upstream-foreward advection for the rainwater sedimentation
!
!-----------------------------------------------------------------------
!
      DO k=2,nz-2
        DO j=1,ny-1
          DO i=1,nx-1
            q(i,j,k,tfuture) = q(i,j,k,tfuture)+ dtrnf *                  &
              ((rhovq(i,j,k+1)-rhovq(i,j,k))/(dz*rhobar(i,j,k)*j3(i,j,k)))
          END DO
        END DO
      END DO
!
!-----------------------------------------------------------------------
!
!  Accumulated precipitation.  The depth of water
!  accumulated per timestep is given by:
!
!  depth (m) = terminal velocity * q * dt * (rhobar/denwater)
!
!  This equation is derived by temporally integrating the vertical
!  flux of rain through the lowest model level.  The total is given
!  for the entire computational domain.
!
!-----------------------------------------------------------------------

      denwater = 1000.0   ! Density of liquid water.


      tem = dtbig1/deltat *dtrnf/denwater *1000.0
                          ! dtbig1/deltat=0.5 for leapfrog scheme

      DO j=1,ny-1
        DO i=1,nx-1
          draing(i,j)=draing(i,j)+rhovq(i,j,2) * tem
        END DO
      END DO
    END DO

  ELSE IF ( impfallopt == 1) THEN    ! implicit scheme
!
!-----------------------------------------------------------------------
!
!  The density weighted rainwater fallout flux at k=2
!
!-----------------------------------------------------------------------
!
    DO j=1,ny-1
      DO i=1,nx-1
        rhovq(i,j,2) = rhobar(i,j,2)*vtr(i,j,2)*q(i,j,2,tfuture)
      END DO
    END DO

    CALL qfallout_im(nx,ny,nz,dtbig1,rhobar,j3,j3inv,vtr, tem1,tem2,  &
                     tem3,tem4,tem5,q)

    DO j=1,ny-1
      DO i=1,nx-1
        draing(i,j)=rhovq(i,j,2) * dtbig1
      END DO
    END DO

  END IF

  CALL acct_stop_inter

  RETURN
END SUBROUTINE qfallout
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE SATADJ                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE satadj(nx,ny,nz, ptprt,qv,qc,ptbar,                          &
           p, ppi, lathvt,  t,qvs, tem1 )

!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Perform the saturation adjustment (condensation and evaporation
!  between qv and qc). Under the conditions of supersaturation
!  (qv > qvs) water vapor, qv, is condensed into cloud water, qc.
!  When the air is subsaturated, (qv <= qvs) cloud water is then
!  evaporated. Adjustments are then made to the potential
!  temperature, water vapor and cloud water fields.

!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  2/29/1992.
!
!  MODIFICATION HISTORY:
!
!  5/02/92 (M. Xue)
!  Added full documentation.
!
!  5/28/92 (K. Brewster)
!  Further facelift.
!
!  11/05/98 (K. Brewster)
!  Added "time constant" to dqv calculation.
!
!  11/08/98 (K. Brewster)
!  Changed calculation of qvs to be based on total pi and p, not pbar
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!    ptprt    Perturbation potential temperature at all time levels(K)
!    pprt     Perturbation pressure at all time levels (Pascal)
!    qv       Water vapor specific humidity at all time levels (kg/kg)
!    qc       Cloud water mixing ratio at all time levels (kg/kg)
!
!    ptbar    Base state potential temperature (K)
!
!  OUTPUT:
!
!    ptprt    Perturbation potential temperature at time tfuture (K)
!    qv       Water vapor specific humidity at time tfuture (kg/kg)
!    qc       Cloud water mixing ratio at time tfuture (kg/kg)
!
!  WORK ARRAYS:
!
!    qvs      Saturation specific humidity (kg/kg)
!    t        Temperature (K)
!
!   (These arrays are defined and used locally (i.e. inside this
!    subroutine), they may also be passed into routines called by
!    this one. Exiting the call to this subroutine, these temporary
!    work arrays may be used for other purposes and therefore their
!    contents may be overwritten. Please examine the usage of work
!    arrays before you alter the code.)
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
!
  INTEGER :: nt           ! The no. of t-levels of t-dependent arrays.
  INTEGER :: tpast        ! Index of time level for the past time.
  INTEGER :: tpresent     ! Index of time level for the present time.
  INTEGER :: tfuture      ! Index of time level for the future time.

  PARAMETER (nt=3, tpast=1, tpresent=2, tfuture=3)
!
  INTEGER :: nx,ny,nz          ! Number of grid points in 3 directions

  REAL :: ptprt (nx,ny,nz,nt)  ! Perturbation potential temperature (K)
  REAL :: qv    (nx,ny,nz,nt)  ! Water vapor specific humidity (kg/kg)
  REAL :: qc    (nx,ny,nz,nt)  ! Cloud water mixing ratio (kg/kg)

  REAL :: ptbar (nx,ny,nz)     ! Base state potential temperature (K)

  REAL :: lathvt(nx,ny,nz)     ! Temperature dependent latent heat of
                               ! vaporization (kg/(m*s**2).
!
!-----------------------------------------------------------------------
!
!  Temporary arrays
!
!-----------------------------------------------------------------------
!
  REAL :: p     (nx,ny,nz)     ! Total pressure at tfuture (Pa)
  REAL :: ppi   (nx,ny,nz)     ! Exner function
  REAL :: t     (nx,ny,nz)     ! Temperature (K)
  REAL :: qvs   (nx,ny,nz)     ! Saturation specific humidity (kg/kg)
  REAL :: tem1  (nx,ny,nz)     ! Temporary arrays
!
!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k,lvlq,lvlt
  REAL :: dqv,qvplus,qcplus
  REAL :: dqvcst
  PARAMETER (dqvcst=1.0)
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'phycst.inc'
!
!-----------------------------------------------------------------------
!
!  Function f_desdt and inline directive for Cray PVP
!
!-----------------------------------------------------------------------
!
  REAL :: f_desdt

!fpp$ expand (f_desdt)
!!dir$ inline always f_desdt
!*$*  inline routine (f_desdt)

!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  lvlq = tfuture
  lvlt = tfuture
!
!-----------------------------------------------------------------------
!
!  Calculate the temperature.
!
!-----------------------------------------------------------------------
!
  DO k=1,nz-1
    DO j=1,ny-1
      DO i=1,nx-1
        t(i,j,k) = (ptprt(i,j,k,lvlt)+ptbar(i,j,k)) * ppi(i,j,k)
      END DO
    END DO
  END DO
!
!-----------------------------------------------------------------------
!
!  Calculate the saturation specific humidity.
!
!-----------------------------------------------------------------------
!
  CALL getqvs(nx,ny,nz, 1,nx-1,1,ny-1,1,nz-1, p, t, qvs)

! write (*,*) "ZUWEN subsatopt/rhsat", subsatopt, rhsat 
  IF (subsatopt /= 0) THEN 

    DO k=1,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          qvs(i,j,k) = rhsat*qvs(i,j,k)  ! adjust qvs for sub-saturation
        END DO
      END DO
    END DO

  END IF 

  DO k=1,nz-1
    DO j=1,ny-1
      DO i=1,nx-1
        tem1(i,j,k) = f_desdt( t(i,j,k) )     ! d(es)/dt/es
        tem1(i,j,k) = (rddrv+(1.-rddrv)*qvs(i,j,k))/rddrv * tem1(i,j,k)
                                              ! d(qvs)/dt/qvs
      END DO
    END DO
  END DO

  DO k=1,nz-1
    DO j=1,ny-1
      DO i=1,nx-1
!
!-----------------------------------------------------------------------
!
!  Calculate the amount of evaporation (when subsaturated, dqv>0)
!  or condensation (when supersaturated, dqv<0).
!
!-----------------------------------------------------------------------
!
        qvplus = MAX(0.0,qv(i,j,k,lvlq))
        qcplus = MAX(0.0,qc(i,j,k,lvlq))

        dqv = (qvs(i,j,k)-qvplus)                                       &
            / (1.0+tem1(i,j,k)*qvs(i,j,k)*lathvt(i,j,k)/cp)
!
!-----------------------------------------------------------------------
!
!  If evaporation occurs (dqv>0), the amount should not exceed
!  the available cloud water.
!
!  Make evaporation and condensation happen gradually.
!
!-----------------------------------------------------------------------
!
        dqv = dqvcst*MIN( dqv, qcplus )
!
!-----------------------------------------------------------------------
!
!  Adjust qv, qc and ptprt accordingly.
!
!-----------------------------------------------------------------------
!
        qv(i,j,k,lvlq) = qv(i,j,k,lvlq) + dqv

        qc(i,j,k,lvlq) = qc(i,j,k,lvlq) - dqv

        ptprt(i,j,k,lvlt) = ptprt(i,j,k,lvlt) -                         &
              dqv*lathvt(i,j,k)/(ppi(i,j,k)*cp)


      END DO
    END DO
  END DO

  RETURN
END SUBROUTINE satadj
 

!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINE QFALLOUT_IM               ######
!######                                                      ######
!######                     Developed by                     ######
!######                                                      ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE qfallout_im(nx,ny,nz,dtbig1,rhobar,j3,j3inv,vtr, &
           a,b,c,d,q_force,q)

!
!  PURPOSE:
!
!  Use vertically implicit scheme to solve the fall-out of hydrometeor
!  given its fallout velocity. The resulting equation is tridiagonal. 
!
!-----------------------------------------------------------------------
!
!
!  AUTHOR: Dan Weber, M. Xue, X. Jin
!  03/20/02
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    nx       Number of grid points in the x-direction
!    ny       Number of grid points in the y-direction
!    nz       Number of grid points in the z-direction
!
!    dtbig1   big time step size (s)
!    q        hydrometeor mixing ratio at all time levels (kg/kg)
!    rhobar   base state air density (kg/m**3)
!    j3       coordinate transformation Jocobian  d(zp)/dz
!    vtr      terminal velocity (m/s)
!
!  OUTPUT:
!
!    q        rainwater mixing ratio at time future (kg/kg)
!
!  WORK ARRAYS:
!
!    a        left of main diagonal, useful range [kbgn+1,kend]
!    b        main diagonal, useful range [kbgn,kend]
!    c        right of main diagonal, useful range [kbgn,kend-1]
!    d        right hand side of equations
!    q_force  a temporary array to account for contributions of 
!             turbulent/computational mixing and  advection in water 
!             substance equation
!
!-----------------------------------------------------------------------
!
!  Variable Declarations:
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE             ! Force explicit declarations
!
  INTEGER :: nx,ny,nz          ! Number of grid points in 3 directions

  REAL :: dtbig1               ! Big time step size (s)
  REAL :: q     (nx,ny,nz,3)  ! Rainwater mixing ratio (kg/kg)
  REAL :: rhobar(nx,ny,nz)     ! Base state air density (kg/m**3)
  REAL :: j3    (nx,ny,nz)     ! Coordinate transformation Jacobian
                               ! defined as d( zp )/d( z ).
  REAL :: j3inv (nx,ny,nz)     ! inverse of j3
  REAL :: vtr   (nx,ny,nz)     ! terminal velocify

!
!-----------------------------------------------------------------------
!
!  temporary arrays:
!
!-----------------------------------------------------------------------
!
  REAL :: a(nx,ny,nz)       ! Left of main diagonal
  REAL :: b(nx,ny,nz)       ! Main diagonal
  REAL :: c(nx,ny,nz)       ! Right of main diagonal
  REAL :: d(nx,ny,nz)       ! Right hand side of equations
  REAL :: q_force(nx,ny,nz)    ! values of q_force for time level n

!
!-----------------------------------------------------------------------
!
! Misc. local variable:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j, k
  REAL :: courant1, courant2, courant3, tema, temb    ! courant number
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'grid.inc'
  INCLUDE 'timelvls.inc'
  INCLUDE 'globcst.inc' 
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

!----------------------------------------------------------------------
!
! calculate q_force from the previously computed adv. and mix. terms.
!
!----------------------------------------------------------------------
!
  temb=0.5/dtbig1
  DO k=2,nz-2
    DO j=1,ny-1
      DO i=1, nx-1
        q_force(i,j,k)=rhobar(i,j,k)*j3(i,j,k)*(q(i,j,k,tfuture)-   &
                       q(i,j,k,tpast))*temb
      END DO
    END DO
  END DO

!----------------------------------------------------------------------
!
! initialize coefficients a, b, c, d
!
!----------------------------------------------------------------------
!
    
  tema = dtbig1*dzinv 
  DO k=2,nz-2
    DO j=1,ny-1
      DO i=1,nx-1
        courant1=vtr(i,j,k-1)*tema   ! courant number at k-1
        courant2=vtr(i,j,k+1)*tema    ! courant number at k+1
!       courant3=vtr(i,j,k)*tema     ! courant number at k
        a(i,j,k)=fallvalpha*rhobar(i,j,k-1)*courant1
        b(i,j,k)=rhobar(i,j,k)*j3(i,j,k)
        c(i,j,k)=-fallvalpha*rhobar(i,j,k+1)*courant2
        d(i,j,k)=2.0*q_force(i,j,k)*dtbig1+rhobar(i,j,k)*              &
                   j3(i,j,k)*q(i,j,k, tpast)+(1.0-fallvalpha)*         &
                   (rhobar(i,j,k+1)*courant2*q(i,j,k+1,tpresent)       &
                   -rhobar(i,j,k-1)*courant1*q(i,j,k-1,tpresent))  
      END DO
    END DO
  END DO
   
!
!----------------------------------------------------------------------
!
! initialize coefficients b at k=2. Bourndary condition is rigid 
!
!----------------------------------------------------------------------
!

  k=2
  DO j=1,ny-1
    DO i=1,nx-1
      b(i,j,k)=b(i,j,k)+a(i,j,k)
    END DO
  END DO

!
!----------------------------------------------------------------------
! 
!  initialize coefficients b for k=nz-2. Bourndary condition is rigid
!
!----------------------------------------------------------------------
! 
     
  k=nz-2
  DO j=1,ny-1    
    DO i=1,nx-1
      b(i,j,k)=b(i,j,k)+c(i,j,k)
    END DO
  END DO
  
!
!----------------------------------------------------------------------
!
! call tridiagonal equation solver
!
!----------------------------------------------------------------------
!
  CALL tridiag2(nx,ny,nz,1,nx-1,1,ny-1,2,nz-2,a,b,c,d)

!
!---------------------------------------------------------------------- 
!
!  assign values of array d to q(i,j,k,tfuture)
!
!----------------------------------------------------------------------
!

  DO k=2,nz-2
    DO j=1,ny-1
      DO i=1,nx-1
        q(i,j,k,tfuture)=d(i,j,k)
      END DO   
    END DO
  END DO  
  
!
!----------------------------------------------------------------------
!
! boundary values at k=1 and nz-1. Boundary is rigid 
!
!----------------------------------------------------------------------
!

  DO j=1,ny-1
    DO i=1,nx-1
      q(i,j,1,tfuture)=q(i,j,2,tfuture)
    END DO  
  END DO

  DO j=1,ny-1
    DO i=1,nx-1
      q(i,j,nz-1,tfuture)=q(i,j,nz-2,tfuture)
    END DO  
  END DO

  RETURN
END SUBROUTINE qfallout_im

