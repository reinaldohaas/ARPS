!########################################################################
!########################################################################
!#########                                                      #########
!#########           SUBROUTINE interface_wrf_bmjdrv            #########
!#########                                                      #########
!#########                     Developed by                     #########
!#########     Center for Analysis and Prediction of Storms     #########
!#########                University of Oklahoma                #########
!#########                                                      #########
!########################################################################
!########################################################################

SUBROUTINE interface_wrf_bmjdrv(nx,ny,nz,pprt,ptprt,qv,pbar,ptbar,zp,   &
                                ptcumsrc,qcumsrc,bmjraincv,prcrate,     &
                                cldefi,xland)

!------------------------------------------------------------------------
!
! PURPOSE:
!
! Interfaces with the WRF version of the Betts-Miller-Janjic convective
! adjustment scheme.
!
!------------------------------------------------------------------------
!
! AUTHOR:  Eric Kemp, 10 October 2001
!
! MODIFICATION HISTORY:
!
! Eric Kemp, 1 November 2001.  Fixed dimension error with array wrf_t.
!
! Eric Kemp, 12 March 2002.  Removed lowlyr array from argument list.
! Replaced with new automatic array wrf_lowlyr.
!
!------------------------------------------------------------------------
!
! Use WRF Betts-Miller-Janjic module.
!
!------------------------------------------------------------------------

  USE module_cu_bmj                             

!------------------------------------------------------------------------
!
! Force explicit declarations.
!
!------------------------------------------------------------------------

  IMPLICIT NONE

!------------------------------------------------------------------------
!
! List include files.
!
!------------------------------------------------------------------------

  INCLUDE 'globcst.inc'
  INCLUDE 'grid.inc'
  INCLUDE 'phycst.inc'

!------------------------------------------------------------------------
!
! Arguments.
! 
!------------------------------------------------------------------------

  INTEGER,INTENT(IN) :: nx,ny,nz           ! Grid dimensions.
  REAL,INTENT(IN) :: pprt(nx,ny,nz)        ! Perturbation pressure (Pa).
  REAL,INTENT(IN) :: ptprt(nx,ny,nz)       ! Perturbation potential
                                           !   temperature (K).
  REAL,INTENT(IN) :: qv(nx,ny,nz)          ! Specific humidity (kg/kg).
  REAL,INTENT(IN) :: pbar(nx,ny,nz)        ! Base-state pressure (Pa).
  REAL,INTENT(IN) :: ptbar(nx,ny,nz)       ! Base-state potential
                                           !   temperature (K).
  REAL,INTENT(IN) :: zp(nx,ny,nz)          ! Height at w-points (m).
  REAL,INTENT(INOUT) :: ptcumsrc(nx,ny,nz) ! Potential temperature
                                           !   tendency.
  REAL,INTENT(INOUT) :: qcumsrc(nx,ny,nz,5)! Moisture tendencies.
                               ! qcumsrc(1,1,1,1) for qv
                               ! qcumsrc(1,1,1,2) for qc
                               ! qcumsrc(1,1,1,3) for qr 
                               ! qcumsrc(1,1,1,4) for qi
                               ! qcumsrc(1,1,1,5) for qs

  REAL,INTENT(INOUT) :: bmjraincv(nx,ny)   ! BMJ rainfall (cm).
  REAL,INTENT(INOUT) :: prcrate(nx,ny)     ! Precipitation rate (mm/s)
  REAL,INTENT(INOUT) :: cldefi(nx,ny)      ! Cloud efficiency
                                           !   (dimensionless).
  REAL,INTENT(IN) :: xland(nx,ny)          ! Land-sea mask (1.0 for land;
                                           !   2.0 for water)

!------------------------------------------------------------------------
!
! WRF 3-D arrays (dimensioned i,k,j).
! 
!------------------------------------------------------------------------

  REAL :: wrf_rr(nx,nz,ny)                 ! Dry air density (kg/m^3)
  REAL :: wrf_rthcuten(nx,nz,ny)           ! Rho_dTheta_m tendency due to 
                                           !   cumulus scheme
                                           !   precipitation 
                                           !   (kg/m^3 . K)
  REAL :: wrf_rqvcuten(nx,nz,ny)           ! Rho_dQv tendency due to
                                           !   cumulus scheme
                                           !   precipitation 
                                           !   (kg/m^3 . kg/kg)
  REAL :: wrf_th(nx,nz,ny)                 ! Potential temperature (K)
  REAL :: wrf_t(nx,nz,ny)                  ! Temperature (K)
  REAL :: wrf_qvmix(nx,nz,ny)              ! Water vapor mixing ratio
                                           !   (kg/kg)
  REAL :: wrf_pint(nx,nz,ny)               ! Pressure at w-points (Pa)
  REAL :: wrf_pmid(nx,nz,ny)               ! Pressure (Pa)
  REAL :: wrf_pi(nx,nz,ny)                 ! Exner function
                                           !   (dimensionless)
  REAL :: wrf_rho(nx,nz,ny)                ! Density (kg/m^3)
  REAL :: wrf_dz8w(nx,nz,ny)               ! dz between full levels (m)

!------------------------------------------------------------------------
!
! Other local arrays and variables.
! 
!------------------------------------------------------------------------

  INTEGER :: wrf_lowlyr(nx,ny)             ! Index of lowest model level.

  REAL :: wrf_raincv(nx,ny)                ! Cumulus scheme 
                                           ! precipitation (mm)

  INTEGER :: ids,ide,jds,jde,kds,kde,                                    &
              ims,ime,jms,jme,kms,kme,                                   &
              its,ite,jts,jte,kts,kte
  INTEGER :: itimestep                     ! Number of timestep.
  INTEGER :: stepcu                        ! Number of fundamental
                                           !   timesteps between
                                           !   convection calls.
  REAL :: d608                             ! rvovrd - 1.

  INTEGER :: i,j,k
  REAL :: qvcumsrctmp,ptcumsrctmp

!------------------------------------------------------------------------
!
! Local parameters.
! 
!------------------------------------------------------------------------

  REAL,PARAMETER :: trfz = 273.15    ! Freezing point (273.15 K)

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  wrf_lowlyr(:,:) = 2
  wrf_raincv(:,:) = 0

!------------------------------------------------------------------------
!
! Set the WRF "domain," "memory," and "tile" dimensions based on the
! ARPS dimensions.  
!
!------------------------------------------------------------------------

  CALL interface_wrf_dims(nx,ny,nz,                                     &
                          ids,ide,jds,jde,kds,kde,                      &
                          ims,ime,jms,jme,kms,kme,                      &
                          its,ite,jts,jte,kts,kte)

!------------------------------------------------------------------------
!
! Set several scalar arguments.
!
!------------------------------------------------------------------------
  
  d608 = rvdrd - 1.

  stepcu = confrq/dtbig ! Number of actual time steps per convection
                        ! time step.
  itimestep = (curtim-tstart)/dtbig ! Current actual time step.

!------------------------------------------------------------------------
!
! Fill WRF 3-D arrays.  Note that WRF arrays are indexed i,k,j.
! 
!------------------------------------------------------------------------

  DO k = 1,nz-1
    DO j = 1,ny-1
      DO i = 1,nx-1
        wrf_pmid(i,k,j) = pprt(i,j,k) + pbar(i,j,k)
        wrf_pi(i,k,j) = (wrf_pmid(i,k,j)*1.0E-5)**rddcp
        wrf_th(i,k,j) = ptprt(i,j,k) + ptbar(i,j,k)
        wrf_t(i,k,j) = wrf_th(i,k,j)*wrf_pi(i,k,j)
        wrf_rr(i,k,j) = wrf_pmid(i,k,j)/(rd*wrf_t(i,k,j)) ! Dry air
                                                          !   density.
        wrf_qvmix(i,k,j) = MAX(0.,qv(i,j,k)/(1. - qv(i,j,k)))
        wrf_rho(i,k,j) = wrf_rr(i,k,j)/(1. + 0.608*wrf_qvmix(i,k,j))
        wrf_dz8w(i,k,j) = zp(i,j,k+1) - zp(i,j,k)
      END DO ! DO i = 1,nx-1
    END DO ! DO j = 1,ny-1
  END DO ! DO k = 1,nz-1

  wrf_pint(:,:,:) = wrf_pmid(:,:,:) ! Not actually used by WRF code,
                                    ! but passed as argument.

!------------------------------------------------------------------------
!
! Call WRF Betts-Miller-Janjic code.
! 
!------------------------------------------------------------------------

  CALL bmjdrv(dtbig,itimestep,stepcu,wrf_rr,wrf_rthcuten,wrf_rqvcuten,  &
              wrf_raincv,wrf_th,wrf_t,wrf_qvmix,wrf_pint,wrf_pmid,      &
              wrf_pi,wrf_rho,wrf_dz8w,cp,rd,rvdrd,lathv,                &
              g,trfz,d608,cldefi,wrf_lowlyr,xland,                      &
              ids,ide,jds,jde,kds,kde,                                  &
              ims,ime,jms,jme,kms,kme,                                  &
              its,ite,jts,jte,kts,kte)

!------------------------------------------------------------------------
!
! Save precipitation rate and accumulation.
! 
!------------------------------------------------------------------------

  prcrate(:,:) = wrf_raincv(:,:)/dtbig  ! mm to mm/s
  bmjraincv(:,:) = wrf_raincv(:,:)*1.E-1   ! mm to cm

!------------------------------------------------------------------------
!
! Convert rho_dthetam_dt and rho_dqvmix_dt to dtheta_dt and dqv_dt.
! 
!------------------------------------------------------------------------

  DO k = 2,nz-2
    DO j = 1,ny-1
      DO i = 1,nx-1
 
        qvcumsrctmp = wrf_rqvcuten(i,k,j)* &
                      (1. - (wrf_qvmix(i,k,j)*wrf_qvmix(i,k,j)))/ &
                      (wrf_rr(i,k,j))
        qcumsrc(i,j,k,1) = qcumsrc(i,j,k,1) + qvcumsrctmp

        ptcumsrctmp = (wrf_rthcuten(i,k,j)/wrf_rr(i,k,j)) - &
                          (rvdrd*wrf_th(i,k,j)*qvcumsrctmp)
        ptcumsrctmp = ptcumsrctmp / &
                          (1. + (rvdrd*wrf_qvmix(i,k,j)))
        ptcumsrc(i,j,k) = ptcumsrc(i,j,k) + ptcumsrctmp

      END DO ! DO i = 1,nx-1
    END DO ! DO j = 1,ny-1
  END DO ! DO k = 2,nz-2

!------------------------------------------------------------------------
!
! The end.
! 
!------------------------------------------------------------------------

  RETURN
END SUBROUTINE interface_wrf_bmjdrv

!########################################################################
!########################################################################
!#########                                                      #########
!#########           SUBROUTINE interface_wrf_bmjinit           #########
!#########                                                      #########
!#########                     Developed by                     #########
!#########     Center for Analysis and Prediction of Storms     #########
!#########                University of Oklahoma                #########
!#########                                                      #########
!########################################################################
!########################################################################

SUBROUTINE interface_wrf_bmjinit(nx,ny,nz,cldefi,restart)

!------------------------------------------------------------------------
!
! PURPOSE:
!
! Initializes variables and look-up tables used by the WRF version of
! the Betts-Miller-Janjic convective adjustment scheme.
!
!------------------------------------------------------------------------
!
! AUTHOR:  Eric Kemp, 12 October 2001
!
! MODIFICATION HISTORY:
!
! Eric Kemp, 12 March 2002.  Removed lowlyr array.
!
!------------------------------------------------------------------------
!
! Use WRF Betts-Miller-Janjic module.
!
!------------------------------------------------------------------------

  USE module_cu_bmj

!------------------------------------------------------------------------
!
! Force explicit declarations.
!
!------------------------------------------------------------------------

  IMPLICIT NONE

!------------------------------------------------------------------------
!
! List include files.
!
!------------------------------------------------------------------------

  INCLUDE 'phycst.inc'

!------------------------------------------------------------------------
!
! Declare arguments.
!
!------------------------------------------------------------------------

  INTEGER,INTENT(IN) :: nx,ny,nz        ! Grid dimensions
  LOGICAL,INTENT(IN) :: restart         ! Restart flag
  REAL,INTENT(INOUT) :: cldefi(nx,ny)   ! BMJ cloud efficiency.

!------------------------------------------------------------------------
!
! Local variables.  Note that the arrays, while passed as arguments
! to the BMJ code, are not actually used by ARPS.  Also, the 3-D
! arrays are dimensioned i,k,j.
!
!------------------------------------------------------------------------

  REAL :: wrf_rthcuten(nx,nz,ny)
  REAL :: wrf_rqvcuten(nx,nz,ny)
  REAL :: wrf_rqccuten(nx,nz,ny)
  REAL :: wrf_rqrcuten(nx,nz,ny)

  INTEGER :: wrf_lowlyr(nx,ny)

  INTEGER :: ids,ide,jds,jde,kds,kde, &
             ims,ime,jms,jme,kms,kme, &
             its,ite,jts,jte,kts,kte

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!------------------------------------------------------------------------
!
! Set the WRF "domain," "memory," and "tile" dimensions based on the
! ARPS dimensions.  
!
!------------------------------------------------------------------------

  CALL interface_wrf_dims(nx,ny,nz,                                     &
                          ids,ide,jds,jde,kds,kde,                      &
                          ims,ime,jms,jme,kms,kme,                      &
                          its,ite,jts,jte,kts,kte)

!------------------------------------------------------------------------
!
! Call WRF Betts-Miller-Janjic initializer.
!
!------------------------------------------------------------------------

  CALL bmjinit(wrf_rthcuten,wrf_rqvcuten,wrf_rqccuten,wrf_rqrcuten,     &
               cldefi,wrf_lowlyr,cp,rd,restart,                         &
               ids,ide,jds,jde,kds,kde,                                 &
               ims,ime,jms,jme,kms,kme,                                 &
               its,ite,jts,jte,kts,kte)

!------------------------------------------------------------------------
!
! The end.
!
!------------------------------------------------------------------------

  RETURN
END SUBROUTINE interface_wrf_bmjinit

!########################################################################
!########################################################################
!#########                                                      #########
!#########             SUBROUTINE interface_wrf_dims            #########
!#########                                                      #########
!#########                     Developed by                     #########
!#########     Center for Analysis and Prediction of Storms     #########
!#########                University of Oklahoma                #########
!#########                                                      #########
!########################################################################
!########################################################################

SUBROUTINE interface_wrf_dims(nx,ny,nz,                                 &
                              ids,ide,jds,jde,kds,kde,                  &
                              ims,ime,jms,jme,kms,kme,                  &
                              its,ite,jts,jte,kts,kte)

!------------------------------------------------------------------------
!
! PURPOSE:
!
! Assigns WRF "domain," "memory," and "tile" dimensions based on ARPS
! dimensions.
!
!------------------------------------------------------------------------
!
! AUTHOR:  Eric Kemp, 11 October 2001
!
!------------------------------------------------------------------------
!
! Force explicit declarations.
!
!------------------------------------------------------------------------

  IMPLICIT NONE

!------------------------------------------------------------------------
!
! Arguments
!
!------------------------------------------------------------------------
  
  INTEGER,INTENT(IN) :: nx,ny,nz                 ! ARPS grid dimensions
  INTEGER,INTENT(OUT) :: ids,ide,jds,jde,kds,kde ! WRF "domain" dims.
  INTEGER,INTENT(OUT) :: ims,ime,jms,jme,kms,kme ! WRF "memory" dims.
  INTEGER,INTENT(OUT) :: its,ite,jts,jte,kts,kte ! WRF "tile" dims.

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!------------------------------------------------------------------------
!
! The "domain" dimensions, while passed as arguments, are not actually
! used by the WRF BMJ scheme.
!
!------------------------------------------------------------------------

  ids = 1
  ide = nx
  jds = 1
  jde = ny
  kds = 1
  kde = nz-2

!------------------------------------------------------------------------
!
! The "memory" dimensions are used to allocate the 3-D and 2-D argument
! arrays in the WRF BMJ code.
!
!------------------------------------------------------------------------

  ims = 1
  ime = nx
  jms = 1
  jme = ny
  kms = 1
  kme = nz

!------------------------------------------------------------------------
!
! The "tile" dimensions are used to allocate the 1-D arrays and are
! used as constraints for the DO loops in the WRF BMJ code.  
! Also, KTE + 1 = KME.
!
!------------------------------------------------------------------------

  its = 1
  ite = nx-1
  jts = 1
  jte = ny-1
  kts = 1
  kte = nz-3
!  kts = 2
!  kte = nz-1

!------------------------------------------------------------------------
!
! The end.
!
!------------------------------------------------------------------------

  RETURN
END SUBROUTINE interface_wrf_dims
