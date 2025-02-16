!##################################################################
!##################################################################
!######                                                      ######
!######              SUBROUTINE MICROPH_NEM                  ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE microph_nem(nx,ny,nz,mphyflg,dtbig1,                         &
           t,p,qv,qscalar,raing,prcrate,                                &
           rhostr,pbar,ptbar,ppi,j3,j3inv,                              &
           rho,rsat,rliq,rice,x,t1,                                     &
           change,rate,maxmelt,rp_nuc,                                  &
           needi,needr,needs,tvq,tem15,tem16,                           &
           tem1,tem2,tem3,tem4,tem5)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Paul Schultz's (FSL) Ice microphysics scheme.  Contains
!  equations for water vapor, cloud water, rain water, ice,
!  snow, and hail.  The scheme is based on:
!
!  Schultz, P., 1995: An Explicit Cloud Physics Parameterization
!  for Operational Numerical Weather Prediction.  Monthly Weather
!  Review, 123, 3331-3343.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Jason J. Levit
!  10/01/1996.
!
!  MODIFICATION HISTORY:
!
!  02/19/97 (J. Levit)
!  A major bug was found in the diffusional crystal growth
!  equation.  The conversion parameter 'c2p' should have been
!  'v2p'.
!
!  01/27/97 (J. Levit)
!  Added a variable 'delt', which is equal to 2.0*dtbig1, to fix
!  a bug associated with the temporal finite differencing in some of
!  the subroutines.  Also, uncommented out the code which allows
!  riming to occur.
!  Cleaned up some more of the code and documentation.
!
!  12/09/1996 (J. Levit)
!  Added some extra code at the end of the subroutine to
!  account for generation of negative values in the precipitation
!  fields.  The terminal velocity (tqv) was set equal to
!  zero before calculation, and a DO LOOP was added to
!  take only the positive values of the precipitation fields
!  once the fallout was calculated.
!
!  11/15/1996 (J. Levit)
!  Completed conversion from RAMS code into ARPS format.
!  Cleaned up code and documentation.
!
!  07/10/1997 (Fanyou Kong)
!  Include MPDCD advection option (sadvopt = 5)
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
!    qc       Cloud water mixing ratio at all time levels (kg/kg)
!    qr       Rainwater mixing ratio at all time levels (kg/kg)
!    qi       Cloud ice mixing ratio at all time levels (kg/kg)
!    qs       Snow mixing ratio at all time levels (kg/kg)
!    qh       Hail mixing ratio at all time levels (kg/kg)
!    raing    Accumulated grid-scale rainfall (mm)
!
!    rhostr   Base state air density times j3 (kg/m**3)
!    pbar     Base state pressure (Pascal)
!    ptbar    Base state potential temperature (K)
!    ppi      Exner function of pressure
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
!    tem1     Temporary work array.
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
  INCLUDE 'nemcst.inc'
  INCLUDE 'globcst.inc'
  INCLUDE 'phycst.inc'
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
!
  INTEGER :: nt           ! The no. of t-levels of t-dependent arrays.
  INTEGER :: tpast        ! Index of time level for the past time.
  INTEGER :: tpresent     ! Index of time level for the present time.
  INTEGER :: tfuture      ! Index of time level for the future time.
  PARAMETER (nt=3, tpast=1, tpresent=2, tfuture=3)

  INTEGER :: nx,ny,nz          ! Number of grid points in 3 directions
  INTEGER :: mphyflg           ! microphysic flag

  REAL :: dtbig1           ! Big time step size (s)
  REAL :: t (nx,ny,nz,nt)  ! Perturbation potential temperature (K)
  REAL :: p (nx,ny,nz,nt)  ! Perturbation pressure (Pascal)

  REAL :: qv    (nx,ny,nz,nt)  ! Water vapor specific humidity (kg/kg)

  REAL, TARGET :: qscalar(nx,ny,nz,nt,nscalar)

  REAL :: raing (nx,ny)        ! Accumulated grid-scale rainfall (mm)
  REAL :: prcrate(nx,ny)       ! Precipitation rate (kg/(m**2*s))


  REAL :: rhostr(nx,ny,nz)     ! Base state air density times j3 (kg/m**3)
  REAL :: pbar  (nx,ny,nz)     ! Base state pressure (Pascal).
  REAL :: ptbar (nx,ny,nz)     ! Base state potential temperature (K)
  REAL :: ppi   (nx,ny,nz)     ! Exner function.
  REAL :: j3    (nx,ny,nz)     ! Coordinate transformation Jacobian defined as
                               ! d( zp )/d( z ).
  REAL :: j3inv (nx,ny,nz)     ! Coordinate transformation Jacobian defined as
                               ! d( zp )/d( z ).
!
!-----------------------------------------------------------------------
!
!  Temporary arrays
!
!-----------------------------------------------------------------------
!
  REAL :: rho(nx,ny,nz)           ! Base state air density (kg/m**3)
  REAL :: rsat    (nx,ny,nz)      ! Temporary work array
  REAL :: rliq    (nx,ny,nz)      ! Temporary work array
  REAL :: rice    (nx,ny,nz)      ! Temporary work array
  REAL :: x       (nx,ny,nz)      ! Temporary work array
  REAL :: t1      (nx,ny,nz)      ! Temporary work array
  REAL :: change  (nx,ny,nz)      ! Temporary work array
  REAL :: rate    (nx,ny,nz)      ! Temporary work array
  REAL :: maxmelt (nx,ny,nz)      ! Temporary work array
  REAL :: rp_nuc  (nx,ny,nz)      ! Temporary work array
  REAL :: needi   (nx,ny,nz)      ! Temporary work array
  REAL :: needr   (nx,ny,nz)      ! Temporary work array
  REAL :: needs   (nx,ny,nz)      ! Temporary work array
  REAL :: tvq     (nx,ny,nz)      ! Temporary work array
  REAL :: tem15   (nx,ny,nz)      ! Temporary work array
  REAL :: tem16   (nx,ny,nz)      ! Temporary work array
  REAL :: tem1    (nx,ny,nz)      ! Temporary work array
  REAL :: tem2    (nx,ny,nz)      ! Temporary work array
  REAL :: tem3    (nx,ny,nz)      ! Temporary work array
  REAL :: tem4    (nx,ny,nz)      ! Temporary work array
  REAL :: tem5    (nx,ny,nz)      ! Temporary work array
!
!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------
!
!
  INTEGER :: i,j,k                ! Temporary variable
  REAL    :: denwater,tema,tem    ! Temporary variable
  REAL    :: delt                 ! Temporary variable
  REAL    :: eff                  ! Temporary variable
  REAL    :: err_num,err_den,dr   ! Temporary variable
  INTEGER :: lvlq                 ! Temporayr variable

  REAL, POINTER :: qc(:,:,:,:)
  REAL, POINTER :: qr(:,:,:,:)
  REAL, POINTER :: qi(:,:,:,:)
  REAL, POINTER :: qs(:,:,:,:)
  REAL, POINTER :: qh(:,:,:,:)
!
!-----------------------------------------------------------------------
!
!  Function f_qvsatl, f_qvsati, f_desdtl, and f_desdti, and inline
!  directive for Cray PVP
!
!-----------------------------------------------------------------------
!
  REAL :: f_qvsatl, f_qvsati, f_desdtl, f_desdti

!fpp$ expand (f_qvsatl)
!fpp$ expand (f_qvsati)
!fpp$ expand (f_desdtl)
!fpp$ expand (f_desdti)
!!dir$ inline always f_qvsatl, f_qvsati, f_desdtl, f_desdti
!*$*  inline routine (f_qvsatl, f_qvsati, f_desdtl, f_desdti)

!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

  lvlq=3

!-----------------------------------------------------------------------
!
! Since this scheme comtain microphysical ice process, it must at least
! define qc, qr, qi, qs, qh
!
!-----------------------------------------------------------------------

  IF (P_QC < 1 .OR. P_QR < 1 .OR. P_QI < 1 .OR. P_QS < 1 .OR. P_QH < 1) THEN

    WRITE(6,'(2a,/,5(a,I2),/,a)')                                       &
               'No enough microphysical array was defined ',            &
               'inside subroutine microph_ic.',                         &
               'P_QC = ',P_QC,' P_QR = ',P_QR,' P_QI = ',P_QI,          &
              ' P_QS = ',P_QS,' P_QH = ',P_QH,                          &
               'Program aborting ...'
    CALL arpsstop('Wrong size for microphysics array, qscalar.',1)

  END IF

  qc => qscalar(:,:,:,:,P_QC)
  qr => qscalar(:,:,:,:,P_QR)
  qi => qscalar(:,:,:,:,P_QI)
  qs => qscalar(:,:,:,:,P_QS)
  qh => qscalar(:,:,:,:,P_QH)

  DO k=1,nz-1
    DO j=1,ny-1
      DO i=1,nx-1
        qv(i,j,k,lvlq)=MAX(0.0,qv(i,j,k,lvlq))
        qc(i,j,k,lvlq)=MAX(0.0,qc(i,j,k,lvlq))
        qr(i,j,k,lvlq)=MAX(0.0,qr(i,j,k,lvlq))
        qi(i,j,k,lvlq)=MAX(0.0,qi(i,j,k,lvlq))
        qs(i,j,k,lvlq)=MAX(0.0,qs(i,j,k,lvlq))
        qh(i,j,k,lvlq)=MAX(0.0,qh(i,j,k,lvlq))
        rho(i,j,k)=rhostr(i,j,k)*j3inv(i,j,k)
      END DO
    END DO
  END DO

  IF(mphyflg == 1) THEN

!C    IF( sadvopt.ge.1.and.sadvopt.le.3) THEN ! Leapfrog scheme
    IF( sadvopt /= 4 .and. tintegopt == 1) THEN                  ! Leapfrog scheme
      delt = 2*dtbig1
    ELSE                                    ! Forward scheme
      delt = dtbig1
    END IF

    DO k=1,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          t(i,j,k,lvlq)=(t(i,j,k,lvlq)+ptbar(i,j,k))*ppi(i,j,k)
          p(i,j,k,lvlq)=(p(i,j,k,lvlq)+pbar(i,j,k))
        END DO
      END DO
    END DO
!
!-----------------------------------------------------------------------
!
!  Nothing in this routine can change the quantities til, which is
!  the temperature with corrections for latent heating by ice and
!  liquid, and rtot.  til is similar to theta-il of Tripoli and
!  Cotton (1981 MWR).  These are checked at the end.
!
!-----------------------------------------------------------------------
!

    DO k=1,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          rliq(i,j,k) = qc(i,j,k,lvlq) + qr(i,j,k,lvlq)
          rice(i,j,k) = qi(i,j,k,lvlq) + qs(i,j,k,lvlq) + qh(i,j,k,lvlq)
        END DO
      END DO
    END DO

! rtot0 = qv(i,j,k,lvlq) + rliq(i,j,k) + rice(i,j,k)
! til0 = t  - (Lvl*rliq(i,j,k)+Lvi*rice(i,j,k))/cp

!
!-----------------------------------------------------------------------
!
!  Condensation and evaporation of liquid.  Cloud liquid is assumed
!  to occur instantaneously.  Evaporation of rain is as in Dudhia
!  and Moncrieff (JAS 89).  No condensational growth of rain.
!
!  The evaporation process eats up liquid before ice, and small
!  particles before large, so the order is cloud liquid, rain,
!  pristine crystals, snow, and finally ice.
!
!-----------------------------------------------------------------------
!
! write (*,*) "ZUWEN subsatopt/rhsat", subsatopt, rhsat 

    DO k=1,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          rsat(i,j,k) = rhsat*f_qvsatl( p(i,j,k,lvlq), t(i,j,k,lvlq) )

          IF ( .NOT. (rliq(i,j,k) == 0. .AND.                             &
                      qv(i,j,k,lvlq) < rsat(i,j,k)) ) THEN

! Call SatAdjL (qv(i,j,k,lvlq), p(i,j,k,lvlq), t(i,j,k,lvlq), rsat(i,j,k), t1)

            t1(i,j,k) = t(i,j,k,lvlq)
            rsat(i,j,k) = rhsat*f_qvsatl( p(i,j,k,lvlq), t1(i,j,k) )
            dr = qv(i,j,k,lvlq) - rsat(i,j,k)
            err_num = t1(i,j,k) - t(i,j,k,lvlq) - lvl/cp * dr
            err_den = 1. + lvl/cp * rsat(i,j,k)      & ! 1+Lv1/cp*d(qvs)/dt
                     * (rddrv+(1.-rddrv)*rsat(i,j,k))                     &
                         * f_desdtl(t1(i,j,k))
  
            t1(i,j,k) = t1(i,j,k) - err_num/err_den
            rsat(i,j,k) = rhsat*f_qvsatl( p(i,j,k,lvlq), t1(i,j,k) )
            dr = qv(i,j,k,lvlq) - rsat(i,j,k)
            err_num = t1(i,j,k) - t(i,j,k,lvlq) - lvl/cp * dr
            err_den = 1. + lvl/cp * rsat(i,j,k)      & ! 1+Lv1/cp*d(qvs)/dt
                     * (rddrv+(1.-rddrv)*rsat(i,j,k))                     &
                         * f_desdtl(t1(i,j,k))
  
            t1(i,j,k) = t1(i,j,k) - err_num/err_den
            rsat(i,j,k) = rhsat*f_qvsatl( p(i,j,k,lvlq), t1(i,j,k) )
            dr = qv(i,j,k,lvlq) - rsat(i,j,k)
            err_num = t1(i,j,k) - t(i,j,k,lvlq) - lvl/cp * dr
            err_den = 1. + lvl/cp * rsat(i,j,k)      & ! 1+Lv1/cp*d(qvs)/dt
                     * (rddrv+(1.-rddrv)*rsat(i,j,k))                     &
                         * f_desdtl(t1(i,j,k))
  
            t1(i,j,k) = t1(i,j,k) - err_num/err_den
            rsat(i,j,k) = rhsat*f_qvsatl( p(i,j,k,lvlq), t1(i,j,k) )
            dr = qv(i,j,k,lvlq) - rsat(i,j,k)
  
            qc(i,j,k,lvlq) = qc(i,j,k,lvlq) + (qv(i,j,k,lvlq)-rsat(i,j,k))
  
            IF (qc(i,j,k,lvlq) > 0.) THEN
  
              qv(i,j,k,lvlq) = rsat(i,j,k)
              t(i,j,k,lvlq) = t1(i,j,k)
  
            ELSE ! evaporation of cloud liquid, then rain
  
              qv(i,j,k,lvlq) = rsat(i,j,k) + qc(i,j,k,lvlq)
              qc(i,j,k,lvlq) = 0.
  
              IF (qr(i,j,k,lvlq) > 0. .AND. qv(i,j,k,lvlq) < rsat(i,j,k)) THEN
                rate(i,j,k) = r2v * (rsat(i,j,k)-qv(i,j,k,lvlq))
                change(i,j,k) = rate(i,j,k) * delt
                IF (change(i,j,k) > (rsat(i,j,k)-qv(i,j,k,lvlq)))         &
                    change(i,j,k)=(rsat(i,j,k)-qv(i,j,k,lvlq))
                IF (change(i,j,k) > qr(i,j,k,lvlq)) THEN
                  qv(i,j,k,lvlq) = qv(i,j,k,lvlq) + qr(i,j,k,lvlq)
                  qr(i,j,k,lvlq) = 0.
                ELSE
                  qv(i,j,k,lvlq) = qv(i,j,k,lvlq) + change(i,j,k)
                  qr(i,j,k,lvlq) = qr(i,j,k,lvlq) - change(i,j,k)
                END IF
              END IF
  
              x(i,j,k) = qr(i,j,k,lvlq) - rliq(i,j,k)
              t(i,j,k,lvlq) = t(i,j,k,lvlq) + x(i,j,k)*lvl/cp
  
            END IF
  
            rliq(i,j,k) = qc(i,j,k,lvlq) + qr(i,j,k,lvlq)
  
          END IF
  
        END DO
      END DO
    END DO

!
!-----------------------------------------------------------------------
!
!  Ice growth and evaporation.  Ice nucleation is more a function of
!  supersaturation than temperature (the same is true of diffusional
!  ice crystal growth).  In the presence of cloud water, ice
!  supersaturation is greatest at -12C, and actually decreases at
!  very cold temperatures.
!
!  Ice nucleation generates ice mass much slower than diffusional
!  growth, so we don't compute that function it if there's any cloud
!  ice present.  Also, no nucleation in the absence of cloud water.
!
!  The diffusion from vapor to pristine crystals is proportional to
!  the vapor excess and is a function of the crystal mass already
!  there.  In the presence of cloud water, the excess is .17 g/kg
!   at 1000 mb and .85 g/kg at 200 mb.  This is a supersaturation of
!  12.4% at -12 C, the temperature at which the difference is
!  greatest.
!
!  Diffusional growth of snow is not allowed at this time.  The
!  water mass will get there anyway via crystal growth and
!  collection.
!
!  If the vapor transfer should be independent of pressure, the
!  equation should be rate(i,j,k) = V2P * (rho*(rv-rsat)) * (rho*rp)]
!
!-----------------------------------------------------------------------
!
    DO k=1,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          rsat(i,j,k) = rhsat*f_qvsati( p(i,j,k,lvlq), t(i,j,k,lvlq) )
          IF (.NOT.(rice(i,j,k) == 0. .AND. qv(i,j,k,lvlq) < rsat(i,j,k))) THEN
  
  ! Call SatAdjI (qv(i,j,k,lvlq), p, t, rsat(i,j,k), t1(i,j,k))
  
            t1(i,j,k) = t(i,j,k,lvlq)
            rsat(i,j,k) = rhsat*f_qvsati( p(i,j,k,lvlq), t1(i,j,k) )
            dr = qv(i,j,k,lvlq) - rsat(i,j,k)
            err_num = t1(i,j,k) - t(i,j,k,lvlq) - lvi/cp * dr
            err_den = 1. + lvl/cp * rsat(i,j,k)      & ! 1+Lv1/cp*d(qvs)/dt
                     * (rddrv+(1.-rddrv)*rsat(i,j,k))                     &
                         * f_desdti(t1(i,j,k))
  
            t1(i,j,k) = t1(i,j,k) - err_num/err_den
            rsat(i,j,k) = rhsat*f_qvsati( p(i,j,k,lvlq), t1(i,j,k) )
            dr = qv(i,j,k,lvlq) - rsat(i,j,k)
            err_num = t1(i,j,k) - t(i,j,k,lvlq) - lvi/cp * dr
            err_den = 1. + lvl/cp * rsat(i,j,k)      & ! 1+Lv1/cp*d(qvs)/dt
                     * (rddrv+(1.-rddrv)*rsat(i,j,k))                     &
                         * f_desdti(t1(i,j,k))
  
            t1(i,j,k) = t1(i,j,k) - err_num/err_den
            rsat(i,j,k) = rhsat*f_qvsati( p(i,j,k,lvlq), t1(i,j,k) )
            dr = qv(i,j,k,lvlq) - rsat(i,j,k)
            err_num = t1(i,j,k) - t(i,j,k,lvlq) - lvi/cp * dr
            err_den = 1. + lvl/cp * rsat(i,j,k)      & ! 1+Lv1/cp*d(qvs)/dt
                     * (rddrv+(1.-rddrv)*rsat(i,j,k))                     &
                         * f_desdti(t1(i,j,k))
  
            t1(i,j,k) = t1(i,j,k) - err_num/err_den
            rsat(i,j,k) = rhsat*f_qvsati( p(i,j,k,lvlq), t1(i,j,k) )
            dr = qv(i,j,k,lvlq) - rsat(i,j,k)
  
            IF (qv(i,j,k,lvlq) > rsat(i,j,k)) THEN ! ice growth
  
              IF (qc(i,j,k,lvlq) > 0. .AND. qi(i,j,k,lvlq) < 1E-6) THEN
  
  !    Call Nucleate_pristine (qv(i,j,k,lvlq), rsat(i,j,k),
  !  :       t(i,j,k,lvlq), rho(i,j,k), rp_nuc(i,j,k))
  
                rp_nuc(i,j,k) = 0.
                IF (.NOT.(t(i,j,k,lvlq) > 268.)) THEN
  
                  IF (rp_nuc(i,j,k) > qi(i,j,k,lvlq)) THEN
                    qv(i,j,k,lvlq) = qv(i,j,k,lvlq) + qi(i,j,k,lvlq)
                    qi(i,j,k,lvlq) = rp_nuc(i,j,k)
                    qv(i,j,k,lvlq) = qv(i,j,k,lvlq) - rp_nuc(i,j,k)
                  END IF
  
                ELSE
  
                  rp_nuc(i,j,k) = 1E3 * EXP(-.639+12.96*(qv(i,j,k,lvlq)   &
                                            / rsat(i,j,k)-1.))            &
                                * pmas / rho(i,j,k)
  
                  IF (rp_nuc(i,j,k) > (qv(i,j,k,lvlq)-rsat(i,j,k))/2.)    &
                      rp_nuc(i,j,k) = (qv(i,j,k,lvlq)-rsat(i,j,k))/2.
  
                  IF (rp_nuc(i,j,k) > qi(i,j,k,lvlq)) THEN
                    qv(i,j,k,lvlq) = qv(i,j,k,lvlq) + qi(i,j,k,lvlq)
                    qi(i,j,k,lvlq) = rp_nuc(i,j,k)
                    qv(i,j,k,lvlq) = qv(i,j,k,lvlq) - rp_nuc(i,j,k)
                  END IF
                END IF
              END IF
  
              IF (qi(i,j,k,lvlq) > 0. .AND. qv(i,j,k,lvlq) > rsat(i,j,k)) THEN
                rate(i,j,k) = v2p*(qv(i,j,k,lvlq)-rsat(i,j,k)) *          &
                              (qi(i,j,k,lvlq)*rho(i,j,k))
                change(i,j,k) = rate(i,j,k) * delt
  
                IF (change(i,j,k) > (qv(i,j,k,lvlq)-rsat(i,j,k)))         &
                    change(i,j,k)=(qv(i,j,k,lvlq)-rsat(i,j,k))
  
                qi(i,j,k,lvlq) = qi(i,j,k,lvlq) + change(i,j,k)
                qv(i,j,k,lvlq) = qv(i,j,k,lvlq) - change(i,j,k)
              END IF
  
            ELSE ! ice evaporation

!
!-----------------------------------------------------------------------
!
!    Pristine crystals.  Might make this instantaneous.
!
!-----------------------------------------------------------------------
!
              IF (qi(i,j,k,lvlq) > 0. .AND. qv(i,j,k,lvlq) < rsat(i,j,k)) THEN
  
                rate(i,j,k) = p2v * (rsat(i,j,k)-qv(i,j,k,lvlq))
                change(i,j,k) = rate(i,j,k) * delt
  
                IF (change(i,j,k) > (rsat(i,j,k)-qv(i,j,k,lvlq)))         &
                    change(i,j,k)=(rsat(i,j,k)-qv(i,j,k,lvlq))
  
                IF (change(i,j,k) > qi(i,j,k,lvlq)) THEN
                  qv(i,j,k,lvlq) = qv(i,j,k,lvlq) + qi(i,j,k,lvlq)
                  qi(i,j,k,lvlq) = 0.
                ELSE
                  qv(i,j,k,lvlq) = qv(i,j,k,lvlq) + change(i,j,k)
                  qi(i,j,k,lvlq) = qi(i,j,k,lvlq) - change(i,j,k)
                END IF
  
              END IF

!
!-----------------------------------------------------------------------
!
!    Then snow.
!
!-----------------------------------------------------------------------
!
              IF (qs(i,j,k,lvlq) > 0. .AND. qv(i,j,k,lvlq) < rsat(i,j,k)) THEN
  
                rate(i,j,k) = s2v * (rsat(i,j,k)-qv(i,j,k,lvlq))
                change(i,j,k) = rate(i,j,k) * delt
  
                IF (change(i,j,k) > (rsat(i,j,k)-qv(i,j,k,lvlq)))         &
                    change(i,j,k)=(rsat(i,j,k)-qv(i,j,k,lvlq))
  
                IF (change(i,j,k) > qs(i,j,k,lvlq)) THEN
                  qv(i,j,k,lvlq) = qv(i,j,k,lvlq) + qs(i,j,k,lvlq)
                  qs(i,j,k,lvlq) = 0.
                ELSE
                  qv(i,j,k,lvlq) = qv(i,j,k,lvlq) + change(i,j,k)
                  qs(i,j,k,lvlq) = qs(i,j,k,lvlq) - change(i,j,k)
                END IF
              END IF

!
!-----------------------------------------------------------------------
!
!    And finally ice.  It might be argued that graupel and hail
!    can be water-coated and thus evaporate wrt liquid
!    saturation (i.e., faster).
!
!-----------------------------------------------------------------------
!
              IF (qh(i,j,k,lvlq) > 0. .AND. qv(i,j,k,lvlq) < rsat(i,j,k)) THEN
  
                rate(i,j,k) = i2v * (rsat(i,j,k)-qv(i,j,k,lvlq))
                change(i,j,k) = rate(i,j,k) * delt
  
                IF (change(i,j,k) > (rsat(i,j,k)-qv(i,j,k,lvlq)))         &
                    change(i,j,k)=(rsat(i,j,k)-qv(i,j,k,lvlq))
  
                IF (change(i,j,k) > qh(i,j,k,lvlq)) THEN
                  qv(i,j,k,lvlq) = qv(i,j,k,lvlq) + qh(i,j,k,lvlq)
                  qh(i,j,k,lvlq) = 0.
                ELSE
                  qv(i,j,k,lvlq) = qv(i,j,k,lvlq) + change(i,j,k)
                  qh(i,j,k,lvlq) = qh(i,j,k,lvlq) - change(i,j,k)
                END IF
              END IF
  
            END IF
  
            x(i,j,k) = (qi(i,j,k,lvlq) + qs(i,j,k,lvlq) +                 &
                       qh(i,j,k,lvlq)) - rice(i,j,k)
            t(i,j,k,lvlq) = t(i,j,k,lvlq) + x(i,j,k)*lvi/cp
            rice(i,j,k) = qi(i,j,k,lvlq) +                                &
                          qs(i,j,k,lvlq) + qh(i,j,k,lvlq)
  
          END IF
  
        END DO
      END DO
    END DO

!
!-----------------------------------------------------------------------
!
!  MELTING AND FREEZING.
!  Melting first.  Cloud ice melts before snow, but they both melt
!  immediately.  Cloud ice melts into cloud liquid.  Snow melts into
!  rain.  Ice also melts into rain, but not immediately.  In all
!  cases, the amount of melting is limited to the amount it takes to
!  cool the parcel to the freezing point.  It works out to about
!  3 g/kg per centigrade degree.  Start by calculating
!  the maximum amount of melting possible in this time step.
!
!-----------------------------------------------------------------------
!

    DO k=1,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
  
          IF (t(i,j,k,lvlq) > 273.1) THEN
  
            maxmelt(i,j,k) = (t(i,j,k,lvlq)-273.1) * cp/lli
  
            IF (qi(i,j,k,lvlq) > 0.) THEN
              IF (maxmelt(i,j,k) < qi(i,j,k,lvlq)) THEN
                qc(i,j,k,lvlq) = qc(i,j,k,lvlq) + maxmelt(i,j,k)
                qi(i,j,k,lvlq) = qi(i,j,k,lvlq) - maxmelt(i,j,k)
                GO TO 501
              ELSE
                maxmelt(i,j,k) = maxmelt(i,j,k) - qi(i,j,k,lvlq)
                qc(i,j,k,lvlq) = qc(i,j,k,lvlq) + qi(i,j,k,lvlq)
                qi(i,j,k,lvlq) = 0.
              END IF
            END IF
  
            IF (qs(i,j,k,lvlq) > 0.) THEN
              IF (maxmelt(i,j,k) < qs(i,j,k,lvlq)) THEN
                qr(i,j,k,lvlq) = qr(i,j,k,lvlq) + maxmelt(i,j,k)
                qs(i,j,k,lvlq) = qs(i,j,k,lvlq) - maxmelt(i,j,k)
                GO TO 501
              ELSE
                maxmelt(i,j,k) = maxmelt(i,j,k) - qs(i,j,k,lvlq)
                qr(i,j,k,lvlq) = qr(i,j,k,lvlq) + qs(i,j,k,lvlq)
                qs(i,j,k,lvlq) = 0.
              END IF
            END IF
  
            IF (qh(i,j,k,lvlq) > 0.) THEN
              rate(i,j,k) = i2r * (t(i,j,k,lvlq)-273.1)/5.
              change(i,j,k) = rate(i,j,k)*delt
              IF (change(i,j,k) > maxmelt(i,j,k)) change(i,j,k)=maxmelt(i,j,k)
              IF (change(i,j,k) < qh(i,j,k,lvlq)) THEN
                qr(i,j,k,lvlq) = qr(i,j,k,lvlq) + change(i,j,k)
                qh(i,j,k,lvlq) = qh(i,j,k,lvlq) - change(i,j,k)
              ELSE
                qr(i,j,k,lvlq) = qr(i,j,k,lvlq) + qh(i,j,k,lvlq)
                qh(i,j,k,lvlq) = 0.
              END IF
            END IF
  
            501           CONTINUE

!
!-----------------------------------------------------------------------
!
!  Now freezing.  First the cloud liquid, then the rain.
!
!-----------------------------------------------------------------------
!
          ELSE
!
!-----------------------------------------------------------------------
!
!  Cloud water stays supercooled well below freezing, but how much?
!  This is just freezing because it's cold.  There is very little
!  cloud water below -25C, and almost none observed below -40C.
!  Ramp parabolically from -20 to -40.
!
!-----------------------------------------------------------------------
!
            IF (qc(i,j,k,lvlq) > 0. .AND. t(i,j,k,lvlq) < 253.) THEN
              rate(i,j,k) = c2p * ((253.-t(i,j,k,lvlq))/20.)**2
              change(i,j,k) = rate(i,j,k) * delt
              IF (change(i,j,k) > qc(i,j,k,lvlq)) THEN
                qi(i,j,k,lvlq) = qi(i,j,k,lvlq) + qc(i,j,k,lvlq)
                qc(i,j,k,lvlq) = 0.
              ELSE
                qi(i,j,k,lvlq) = qi(i,j,k,lvlq) + change(i,j,k)
                qc(i,j,k,lvlq) = qc(i,j,k,lvlq) - change(i,j,k)
              END IF
            END IF
!
!-----------------------------------------------------------------------
!
!  Rain freezing into ice; parabolic function similar to C2P.
!  Based loosely on Fig.1 from Cotton (MWR 72b).
!
!-----------------------------------------------------------------------
!
            IF (qr(i,j,k,lvlq) > 0. .AND. t(i,j,k,lvlq) < 267.) THEN
              rate(i,j,k) = r2i * ((267.-t(i,j,k,lvlq))/14.)**2
              change(i,j,k) = rate(i,j,k) * delt
              IF (change(i,j,k) > qr(i,j,k,lvlq)) THEN
                qh(i,j,k,lvlq) = qh(i,j,k,lvlq) + qr(i,j,k,lvlq)
                qr(i,j,k,lvlq) = 0.
              ELSE
                qh(i,j,k,lvlq) = qh(i,j,k,lvlq) + change(i,j,k)
                qr(i,j,k,lvlq) = qr(i,j,k,lvlq) - change(i,j,k)
              END IF
            END IF
  
          END IF
!
!-----------------------------------------------------------------------
!
!  Temperature after freezing or melting.
!  This should give the same answer.
!  x = (qi(i,j,k,lvlq) + qs(i,j,k,lvlq) + qh(i,j,k,lvlq)) -
!    :     rice(i,j,k)
!  t = t + x*Lli/cp
!
!-----------------------------------------------------------------------
!
          x(i,j,k) = (qc(i,j,k,lvlq) + qr(i,j,k,lvlq)) - rliq(i,j,k)
          t(i,j,k,lvlq) = t(i,j,k,lvlq) - x(i,j,k)*lli/cp
  
          rliq(i,j,k) = qc(i,j,k,lvlq) + qr(i,j,k,lvlq)
          rice(i,j,k) = qi(i,j,k,lvlq) + qs(i,j,k,lvlq) + qh(i,j,k,lvlq)
  
        END DO
      END DO
    END DO

!
!-----------------------------------------------------------------------
!
!  COLLECTION.  These processes are determined by spacing between
!  particles regardless of how much gas is also in the volume.  So
!  we'll first convert to specific contents, and then later back to
!  mixing ratios.
!
!-----------------------------------------------------------------------
!
    DO k=1,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          qc(i,j,k,lvlq) = qc(i,j,k,lvlq) * rho(i,j,k)
          qi(i,j,k,lvlq) = qi(i,j,k,lvlq) * rho(i,j,k)
          qr(i,j,k,lvlq) = qr(i,j,k,lvlq) * rho(i,j,k)
          qs(i,j,k,lvlq) = qs(i,j,k,lvlq) * rho(i,j,k)
          qh(i,j,k,lvlq) = qh(i,j,k,lvlq) * rho(i,j,k)
        END DO
      END DO
    END DO

!
!-----------------------------------------------------------------------
!
!  Autoconversion.  As soon as you build up enough cloud matter, it
!  starts converting to rain or snow.  This "nucleates" the
!  collection process, which is nonlinear.  The nucleated amount
!  determines how long before rapid collection occurs.
!  If some precipitate is already present, the autoconv procedure
!  just makes sure there's enough.
!
!-----------------------------------------------------------------------
!
    DO k=1,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          IF (qc(i,j,k,lvlq) > qcmin) THEN
            change(i,j,k) = qc(i,j,k,lvlq) - qcmin
            IF (change(i,j,k) > qr(i,j,k,lvlq)) THEN
              change(i,j,k) = change(i,j,k) - qr(i,j,k,lvlq)
              qr(i,j,k,lvlq) = qr(i,j,k,lvlq) + change(i,j,k)
              qc(i,j,k,lvlq) = qc(i,j,k,lvlq) - change(i,j,k)
            END IF
          END IF
        END DO
      END DO
    END DO
  
    DO k=1,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          IF (qi(i,j,k,lvlq) > qpmin) THEN
            change(i,j,k) = qi(i,j,k,lvlq) - qpmin
            IF (change(i,j,k) > qs(i,j,k,lvlq)) THEN
              change(i,j,k) = change(i,j,k) - qs(i,j,k,lvlq)
              qs(i,j,k,lvlq) = qs(i,j,k,lvlq) + change(i,j,k)
              qi(i,j,k,lvlq) = qi(i,j,k,lvlq) - change(i,j,k)
            END IF
          END IF
        END DO
      END DO
    END DO

!
!-----------------------------------------------------------------------
!
!  There can be a three-way competition for cloud water among the
!  rain, snow, and ice, so we distribute it into the three
!  categories.  The "need" variables are the amount of cloud water
!  that process would use up in a time step if it didn't have to
!  compete.  The collection formula for rain is very nearly the same
!  as Soong and Ogura (1973); the other functions are based on that.
!
!-----------------------------------------------------------------------
!
    DO k=1,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          IF (qc(i,j,k,lvlq) > 0.) THEN
            rate(i,j,k) = c2r * qc(i,j,k,lvlq) * qr(i,j,k,lvlq)
            needr(i,j,k) = rate(i,j,k)*delt
            rate(i,j,k) = c2s * qc(i,j,k,lvlq) * qs(i,j,k,lvlq)
            needs(i,j,k) = rate(i,j,k)*delt
            rate(i,j,k) = c2i * qc(i,j,k,lvlq) * qh(i,j,k,lvlq)
            needi(i,j,k) = rate(i,j,k)*delt
            tem15(i,j,k) = needr(i,j,k) + needs(i,j,k) + needi(i,j,k)
            IF (tem15(i,j,k) > qc(i,j,k,lvlq)) THEN
              needr(i,j,k) = needr(i,j,k) * qc(i,j,k,lvlq) / tem15(i,j,k)
              needs(i,j,k) = needs(i,j,k) * qc(i,j,k,lvlq) / tem15(i,j,k)
              needi(i,j,k) = needi(i,j,k) * qc(i,j,k,lvlq) / tem15(i,j,k)
            END IF

!
!-----------------------------------------------------------------------
!
!  The riming process nucleates a little cloud ice.  Until a better
!  number comes along we'll say 1% of the collected liquid, both for
!  snow and graupel.
!
!-----------------------------------------------------------------------
!
            qi(i,j,k,lvlq) = qi(i,j,k,lvlq) + .01*needs(i,j,k)
            needs(i,j,k) = needs(i,j,k) - .01*needs(i,j,k)
            qi(i,j,k,lvlq) = qi(i,j,k,lvlq) + .01*needi(i,j,k)
            needi(i,j,k) = needi(i,j,k) - .01*needi(i,j,k)
!
!-----------------------------------------------------------------------
!
!  For simplicity (i.e., until a better way to do it comes along),
!  we assume that the result of snow riming is to convert the
!  cloud water, but not the snow, to the graupel category.
!
!-----------------------------------------------------------------------
!
            qr(i,j,k,lvlq) = qr(i,j,k,lvlq) + needr(i,j,k)
            qh(i,j,k,lvlq) = qh(i,j,k,lvlq) + needi(i,j,k) + needs(i,j,k)
            qc(i,j,k,lvlq) = qc(i,j,k,lvlq) - needr(i,j,k) -              &
                             needs(i,j,k) - needi(i,j,k)

          END IF

        END DO
      END DO
    END DO

!
!-----------------------------------------------------------------------
!
!  Unlike the collection of cloud water, which can be complete, we
!  don't want to zero out the pristine crystals.  For one thing,
!  some are so tiny they won't get collected, but also, we want to
!  be able to produce more condensate if it still supersaturated wrt
!  ice, which won't happen if there's zero cloud ice.  We'll leave
!  behind the equivalent of 100 per liter.  At 1E-11 kg per crystal
!  and 1000 liters per m**3, that would be 1E-6 kg/m**3 (1 mg/m**3).
!  The temperature-dependent efficiency follows Lin et al (JCAM 83).
!
!-----------------------------------------------------------------------
!
    DO k=1,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          IF (qi(i,j,k,lvlq) > 0. .AND. qs(i,j,k,lvlq) > 0.) THEN
            eff = 1. - (273.1-t(i,j,k,lvlq))/50.
            rate(i,j,k) = p2s * eff * qi(i,j,k,lvlq) * qs(i,j,k,lvlq)
            change(i,j,k) = rate(i,j,k) * delt
            IF (change(i,j,k) > qi(i,j,k,lvlq)) THEN
              qs(i,j,k,lvlq) = qs(i,j,k,lvlq) + qi(i,j,k,lvlq) - 1E-6
              qi(i,j,k,lvlq) = 1E-6
            ELSE
              qs(i,j,k,lvlq) = qs(i,j,k,lvlq) + change(i,j,k)
              qi(i,j,k,lvlq) = qi(i,j,k,lvlq) - change(i,j,k)
            END IF
          END IF
        END DO
      END DO
    END DO

!
!-----------------------------------------------------------------------
!
!   Convert back to mixing ratios.
!
!-----------------------------------------------------------------------
!
    DO k=1,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          qc(i,j,k,lvlq) = qc(i,j,k,lvlq) / rho(i,j,k)
          qi(i,j,k,lvlq) = qi(i,j,k,lvlq) / rho(i,j,k)
          qr(i,j,k,lvlq) = qr(i,j,k,lvlq) / rho(i,j,k)
          qs(i,j,k,lvlq) = qs(i,j,k,lvlq) / rho(i,j,k)
          qh(i,j,k,lvlq) = qh(i,j,k,lvlq) / rho(i,j,k)
  
!
!-----------------------------------------------------------------------
!
!  Calculate temperature after phase change(i,j,k)s resulting from
!  collection.
!
!-----------------------------------------------------------------------
!
          x(i,j,k) = (qc(i,j,k,lvlq) + qr(i,j,k,lvlq)) - rliq(i,j,k)
          t(i,j,k,lvlq) = t(i,j,k,lvlq) - x(i,j,k)*lli/cp
          rliq(i,j,k) = qc(i,j,k,lvlq) + qr(i,j,k,lvlq)
          rice(i,j,k) = qi(i,j,k,lvlq) + qs(i,j,k,lvlq) + qh(i,j,k,lvlq)
        END DO
      END DO
    END DO

!
!-----------------------------------------------------------------------
!
!  Now comes the final saturation adjustment, if necessary.
!
!-----------------------------------------------------------------------
!
    DO k=1,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          rsat(i,j,k) = rhsat*f_qvsatl( p(i,j,k,lvlq), t(i,j,k,lvlq) )
  
          IF (.NOT.(qc(i,j,k,lvlq) == 0. .AND. qv(i,j,k,lvlq) < rsat(i,j,k))) THEN
  
  ! Call SatAdjL (qv(i,j,k,lvlq), p, t, rsat(i,j,k), t1)
  
            t1(i,j,k) = t(i,j,k,lvlq)
            dr = qv(i,j,k,lvlq) - rsat(i,j,k)
            err_num = t1(i,j,k) - t(i,j,k,lvlq) - lvl/cp * dr
            err_den = 1. + lvl/cp * rsat(i,j,k) * f_desdtl( t1(i,j,k) )
  
            t1(i,j,k) = t1(i,j,k) - err_num/err_den
            rsat(i,j,k) = rhsat*f_qvsatl( p(i,j,k,lvlq), t1(i,j,k) )
            dr = qv(i,j,k,lvlq) - rsat(i,j,k)
            err_num = t1(i,j,k) - t(i,j,k,lvlq) - lvl/cp * dr
            err_den = 1. + lvl/cp * rsat(i,j,k) * f_desdtl( t1(i,j,k) )
  
            t1(i,j,k) = t1(i,j,k) - err_num/err_den
            rsat(i,j,k) = rhsat*f_qvsatl( p(i,j,k,lvlq), t1(i,j,k) )
            dr = qv(i,j,k,lvlq) - rsat(i,j,k)
            err_num = t1(i,j,k) - t(i,j,k,lvlq) - lvl/cp * dr
            err_den = 1. + lvl/cp * rsat(i,j,k) * f_desdtl( t1(i,j,k) )
  
            t1(i,j,k) = t1(i,j,k) - err_num/err_den
            rsat(i,j,k) = rhsat*f_qvsatl( p(i,j,k,lvlq), t1(i,j,k) )
            dr = qv(i,j,k,lvlq) - rsat(i,j,k)
  
            qc(i,j,k,lvlq) = qc(i,j,k,lvlq) + (qv(i,j,k,lvlq)-rsat(i,j,k))
  
            IF (qc(i,j,k,lvlq) > 0.) THEN
              qv(i,j,k,lvlq) = rsat(i,j,k)
            ELSE
              qv(i,j,k,lvlq) = rsat(i,j,k) + qc(i,j,k,lvlq)
              qc(i,j,k,lvlq) = 0.
            END IF
  
            x(i,j,k) = qc(i,j,k,lvlq) + qr(i,j,k,lvlq) - rliq(i,j,k)
            t(i,j,k,lvlq) = t(i,j,k,lvlq) + x(i,j,k)*lvl/cp
            rliq(i,j,k) = qc(i,j,k,lvlq) + qr(i,j,k,lvlq)
  
          END IF
  
        END DO
      END DO
    END DO

    DO k=1,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          t(i,j,k,lvlq)=(t(i,j,k,lvlq)/ppi(i,j,k))-ptbar(i,j,k)
          p(i,j,k,lvlq)=p(i,j,k,lvlq)-pbar(i,j,k)
          qv(i,j,k,lvlq)=MAX(0.0,qv(i,j,k,lvlq))
          qc(i,j,k,lvlq)=MAX(0.0,qc(i,j,k,lvlq))
          qr(i,j,k,lvlq)=MAX(0.0,qr(i,j,k,lvlq))
          qi(i,j,k,lvlq)=MAX(0.0,qi(i,j,k,lvlq))
          qs(i,j,k,lvlq)=MAX(0.0,qs(i,j,k,lvlq))
          qh(i,j,k,lvlq)=MAX(0.0,qh(i,j,k,lvlq))
        END DO
      END DO
    END DO

  END IF  ! mphyflg == 1
!
!-----------------------------------------------------------------------
!
!  Conservation checks.
!
!-----------------------------------------------------------------------
!

! rtot1(i,j,k) = qv(i,j,k,lvlq) + rliq(i,j,k) + rice(i,j,k)
! If (abs(rtot0-rtot1(i,j,k)).gt..000001) then
!  print*, 'rtot check', rtot0, rtot1(i,j,k)
!  print*, qv(i,j,k,lvlq), qc(i,j,k,lvlq), qr(i,j,k,lvlq), qi(i,j,k,lvlq), qs(i,j,k,lvlq), qh(i,j,k,lvlq)
! End if

! til1 = t - (Lvl*rliq(i,j,k)+Lvi*rice(i,j,k))/cp
! If (abs(til1-til0).gt..0001) then
!  pqh(i,j,k,lvlq)nt*, 'til check', til0, til1
! End if

!
!-----------------------------------------------------------------------
!
!  These terminal velocity formulations are similar in form to Ogura
!  and Takahashi (1971).  The curves for rain and snow were tweaked
!  until they matched the curves on page 241 of Pielke's book (they
!  are very, very close).  The curve for ice is based on the curve
!  for rain; the only difference is in the exponent.  The effect is
!  that small values of ice, presumed to be heavily-rimed snow, fall
!  slower than rain of the same concentration, but higher values,
!  presumed to be big graupel or hail, fall faster than rain.  The
!  transition is at 1 g/kg.  The terminal velocity for pristine
!  crystals does not have a dependency on mixing ratio like the
!  others, because the others incorporate(i,j,k) the assumption that
!  higher mixing ratios imply bigger particles, which fall faster,
!  but that's not the case for pristine crystals.
!
!-----------------------------------------------------------------------
!
  DO k=1,nz-1
    DO j=1,ny-1
      DO i=1,nx-1
        tvq(i,j,k)=0.0
        IF (qi(i,j,k,lvlq) > 0.) tvq(i,j,k) = 0.5 * SQRT(1./rho(i,j,k))
      END DO
    END DO
  END DO

  CALL qfallout(nx,ny,nz,dtbig1,qi,rho,j3,j3inv,tem15(1,1,1),           &
                tvq,tem16,tem1,tem2,tem3,tem4,tem5)

  DO k=1,nz-1
    DO j=1,ny-1
      DO i=1,nx-1
        tvq(i,j,k)=0.0
        IF (qr(i,j,k,lvlq) > 0.) tvq(i,j,k) =                           &
            5.5 * (1000.*qr(i,j,k,lvlq))**.125 * SQRT(1./rho(i,j,k))
      END DO
    END DO
  END DO

  CALL qfallout(nx,ny,nz,dtbig1,qr,rho,j3,j3inv,tem15(1,1,2),           &
                tvq,tem16,tem1,tem2,tem3,tem4,tem5)

  DO k=1,nz-1
    DO j=1,ny-1
      DO i=1,nx-1
        tvq(i,j,k)=0.0
        IF (qs(i,j,k,lvlq) > 0.) tvq(i,j,k) =                           &
            2.0 * (1000.*qs(i,j,k,lvlq))**.100 * SQRT(1./rho(i,j,k))
      END DO
    END DO
  END DO

  CALL qfallout(nx,ny,nz,dtbig1,qs,rho,j3,j3inv,tem15(1,1,3),           &
                tvq,tem16,tem1,tem2,tem3,tem4,tem5)

  DO k=1,nz-1
    DO j=1,ny-1
      DO i=1,nx-1
        tvq(i,j,k)=0.0
        IF (qh(i,j,k,lvlq) > 0.) tvq(i,j,k) = 5.5 *                     &
            (1000.*qh(i,j,k,lvlq))**.333 * SQRT(1./rho(i,j,k))
      END DO
    END DO
  END DO

  CALL qfallout(nx,ny,nz,dtbig1,qh,rho,j3,j3inv,tem15(1,1,4),           &
                tvq,tem16,tem1,tem2,tem3,tem4,tem5)

  denwater = 1000.0   ! Density of liquid water (kg/m**3)
  tem = denwater/(1000.0*dtbig1)

  DO j=1,ny-1
    DO i=1,nx-1
      tema=tem15(i,j,2)+tem15(i,j,3)+tem15(i,j,4)
      raing(i,j)=raing(i,j)+tema
      prcrate(i,j) = tema*tem
    END DO
  END DO

  DO k=1,nz-1
    DO j=1,ny-1
      DO i=1,nx-1
        qr(i,j,k,lvlq)=MAX(0.0,qr(i,j,k,lvlq))
        qi(i,j,k,lvlq)=MAX(0.0,qi(i,j,k,lvlq))
        qs(i,j,k,lvlq)=MAX(0.0,qs(i,j,k,lvlq))
        qh(i,j,k,lvlq)=MAX(0.0,qh(i,j,k,lvlq))
      END DO
    END DO
  END DO

  RETURN
END SUBROUTINE microph_nem
