!########################################################################
!########################################################################
!#########                                                      #########
!#########           SUBROUTINE interface_wrf_kfetadrv          #########
!#########                                                      #########
!#########                     Developed by                     #########
!#########     Center for Analysis and Prediction of Storms     #########
!#########                University of Oklahoma                #########
!#########                                                      #########
!########################################################################
!########################################################################

SUBROUTINE interface_wrf_kfetadrv(nx,ny,nz,u,v,w0avg,                 &
                                pprt,ptprt,qv,pbar,ptbar,zp,          &
                                ptcumsrc,qcumsrc,prcrate,             &
                                inca,raincv)

!------------------------------------------------------------------------
!
! PURPOSE:
!
! Interfaces with the WRF version of the new Kain-Fritsch cumulus
! parameterization scheme (April 2002 version).
!
!------------------------------------------------------------------------
!
! AUTHOR:  Fanyou Kong, April 2002
!
! MODIFICATION HISTORY:
!
!------------------------------------------------------------------------
!
! Use WRF KF_ETA module.
!
!------------------------------------------------------------------------

  USE module_cu_kfeta

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
  REAL,INTENT(IN) :: u   (nx,ny,nz)        ! x component of velocity (m/s)
  REAL,INTENT(IN) :: v   (nx,ny,nz)        ! y component of velocity (m/s)
  REAL,INTENT(IN) :: w0avg(nx,ny,nz)       ! Running average w (m/s)
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

  INTEGER,INTENT(INOUT) :: inca(nx,ny) ! Counter of the cloud relaxation
                                       !   time in KF scheme (integer)
  REAL,INTENT(INOUT) :: prcrate(nx,ny) ! Precipitation rate (mm/s)
  REAL,INTENT(INOUT) :: raincv(nx,ny)  ! KF Precipitation (cm)

!------------------------------------------------------------------------
!
! WRF 1-D arrays
! 
!------------------------------------------------------------------------

  REAL, DIMENSION(nz-3) :: u1d,v1d,t1d,dz1d,qv1d,p1d,rh1d,w0avg1d,pi1d
  REAL, DIMENSION(nz-3) :: dqdt,dqidt,dqcdt,dqrdt,dqsdt,dtdt

!------------------------------------------------------------------------
!
! Other local arrays and variables.
! 
!------------------------------------------------------------------------

  REAL :: nca(nx,ny)                       ! Real type of inca(nx,ny)

  INTEGER :: ids,ide,jds,jde,kds,kde,                                    &
              ims,ime,jms,jme,kms,kme,                                   &
              its,ite,jts,jte,kts,kte
  INTEGER :: stepcu                        ! Number of fundamental
                                           !   timesteps between
                                           !   convection calls.
  LOGICAL :: warm_rain                     ! Warm rain flag
  LOGICAL :: cu_act_flag(nx,ny)            ! used in WRF model

  REAL :: dxsq
  REAL :: d608                             ! rv/rd - 1.
  REAL :: xlv0,xlv1,xls0,xls1,svp1,svp2,svp3,svpt0

!  INTEGER :: p_qc,p_qr,p_qi,p_qs,p_first_scalar
  INTEGER :: p_first_scalar
  INTEGER :: i,j,k

  DATA xlv0,xlv1,xls0,xls1/3.15e6,2370.0,2.905e6,259.532/
  DATA svp1,svp2,svp3,svpt0/0.6112,17.67,29.65,273.15/

  DATA p_first_scalar /1/ 
!  DATA p_qc,p_qr,p_qi,p_qs,p_first_scalar/1,1,1,1,2/ ! initial values of
                                           ! microphy. index used in WRF 
                                           ! (p_qx >= p_first_scalar to 
                                           ! be in effect)
!------------------------------------------------------------------------
!
! Local parameters.
! 
!------------------------------------------------------------------------

  warm_rain = .false.
  IF ( mphyopt == 0 .OR. mphyopt == 1 ) THEN
    warm_rain = .true.
!    p_qc = 3            ! they should be set inside initpara3d.f90
!    p_qr = 4
!  ELSE IF ( mphyopt == 2 ) THEN
  ELSE        ! Added mphyopt == 4 then how about mphyopt == 3? Ask Nate.
!    p_qc = 3            ! they should be set inside initpara3d.f90
!    p_qr = 4
!    p_qi = 5
!    p_qs = 6
  ENDIF

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
! Set several scalar arguments.
!
!------------------------------------------------------------------------
  
  dxsq=dx*dx
  d608 = rvdrd - 1.

  stepcu = confrq/dtbig ! Number of actual time steps per convection
                        ! time step.

  DO i=its,ite
    DO j=jts,jte
      cu_act_flag(i,j) = .true.       ! WRF array, not used in ARPS model
      nca(i,j)=float(inca(i,j))
    ENDDO
  ENDDO

!------------------------------------------------------------------------
! I,J loop to call KF_ETA
!------------------------------------------------------------------------

  DO i=its,ite
    DO j=jts,jte

     IF ( inca(i,j) .gt. 0 ) THEN
       cu_act_flag(i,j) = .false.
     ELSE

       DO k=kts,kte
          dqdt(k)=0.
          dqidt(k)=0.
          dqcdt(k)=0.
          dqrdt(k)=0.
          dqsdt(k)=0.
          dtdt(k)=0.
       ENDDO
       raincv(i,j)=0.

       DO k=kts,kte
         u1d(k)=0.5*(u(i,j,k+1)+u(i+1,j,k+1))
         v1d(k)=0.5*(v(i,j,k+1)+v(i,j+1,k+1))
         qv1d(k)=qv(i,j,k+1)
         p1d(k)=pbar(i,j,k+1)+pprt(i,j,k+1)
         w0avg1d(k)=0.5*(w0avg(i,j,k+1)+w0avg(i,j,k+2))
         dz1d(k)=zp(i,j,k+2)-zp(i,j,k+1)
         pi1d(k)=((pbar(i,j,k+1)+pprt(i,j,k+1))/p0)**rddcp
       ENDDO
       DO k=kts,kte
         t1d(k)=(ptbar(i,j,k+1)+ptprt(i,j,k+1))*pi1d(k)
       ENDDO
       DO k=kts,kte
         rh1d(k)=p1d(k)/(rd*t1d(k)*(1.+d608*qv1d(k)))
       ENDDO

!if(i == 5 .and. j == 5) then
!write(17,*) 'u1d,v1d,qv1d,p1d,w0avg1d,dz1d,pi1d,t1d,rh1d:'
!write(17,'(9e10.3)') (u1d(k),v1d(k),qv1d(k),p1d(k),w0avg1d(k),dz1d(k), &
!                             pi1d(k),t1d(k),rh1d(k),k=kts,kte)
!endif

!------------------------------------------------------------------------
!
! Call WRF KF_ETA code
! 
!------------------------------------------------------------------------

       CALL KF_eta_PARA(i, j,                  &
            u1d,v1d,t1d,qv1d,p1d,dz1d,         &
            w0avg1d,dtbig,dx,dxsq,rh1d,        &
            xlv0,xlv1,xls0,xls1,cp,rd,g,       &
            rddrv,svp1,svp2,svp3,svpt0,        &
            dqdt,dqidt,dqcdt,dqrdt,dqsdt,dtdt, &
            raincv,nca,stepcu,kffbfct,kfsubsattrig,  &  ! modified to pass
                                                        ! kffbfct and
                                                        ! kfsubsattrig
            p_qi,p_qs,p_first_scalar,warm_rain,&
            ids,ide, jds,jde, kds,kde,         &
            ims,ime, jms,jme, kms,kme,         &
            its,ite, jts,jte, kts,kte)

       DO k=kts,kte
         ptcumsrc(i,j,k+1)=ptcumsrc(i,j,k+1)+dtdt(k)/pi1d(k)
         qcumsrc(i,j,k+1,1)=qcumsrc(i,j,k+1,1)+dqdt(k)
       ENDDO

       IF ( p_qc >= p_first_scalar ) THEN
         DO k=kts,kte
           qcumsrc(i,j,k+1,2)=dqcdt(k)
         ENDDO
       ENDIF
       IF ( p_qr >= p_first_scalar ) THEN
         DO k=kts,kte
           qcumsrc(i,j,k+1,3)=dqrdt(k)
         ENDDO
       ENDIF

!......     QSTEN STORES GRAUPEL TENDENCY IF IT EXISTS, OTHERISE SNOW (V2)
       IF ( p_qi >= p_first_scalar ) THEN
         DO k=kts,kte
           qcumsrc(i,j,k+1,4)=dqidt(k)
         ENDDO
       ENDIF

       IF ( p_qs >= p_first_scalar ) THEN
         DO k=kts,kte
           qcumsrc(i,j,k+1,5)=dqsdt(k)
         ENDDO
       ENDIF

       raincv(i,j)=0.1*raincv(i,j)    ! new KF_ETA output raincv in (mm)

     ENDIF

    END DO
  END DO       ! End I,J loop

!------------------------------------------------------------------------
!
! Save precipitation rate and restore integer nca
! 
!------------------------------------------------------------------------

  DO i=its,ite
    DO j=jts,jte
      prcrate(i,j) = 10.0*raincv(i,j)/dtbig
    END DO
  END DO

  DO i=its,ite
    DO j=jts,jte
      inca(i,j)=nint(nca(i,j))
    END DO
  END DO

!------------------------------------------------------------------------
!
! The end.
! 
!------------------------------------------------------------------------

  RETURN
END SUBROUTINE interface_wrf_kfetadrv

!########################################################################
!########################################################################
!#########                                                      #########
!#########           SUBROUTINE interface_wrf_kfinit            #########
!#########                                                      #########
!#########                     Developed by                     #########
!#########     Center for Analysis and Prediction of Storms     #########
!#########                University of Oklahoma                #########
!#########                                                      #########
!########################################################################
!########################################################################

SUBROUTINE interface_wrf_kfinit(nx,ny,nz,inca,restart)

!------------------------------------------------------------------------
!
! PURPOSE:
!
! Initializes variables and look-up tables used by the WRF version of
! the new Kain-Fritsch convective adjustment scheme.
!
!------------------------------------------------------------------------
!
! AUTHOR:  Fanyou Kong, April 2002
!
! MODIFICATION HISTORY:
!
!------------------------------------------------------------------------
!
! Use WRF KF_ETA module.
!
!------------------------------------------------------------------------

  USE module_cu_kfeta

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

!  INCLUDE 'phycst.inc'

!------------------------------------------------------------------------
!
! Declare arguments.
!
!------------------------------------------------------------------------

  INTEGER,INTENT(IN) :: nx,ny,nz       ! Grid dimensions
  INTEGER,INTENT(INOUT) :: inca(nx,ny) ! Counter of the cloud relaxation
                                       !   time in KF scheme (integer)
  LOGICAL,INTENT(IN) :: restart        ! Restart flag

!------------------------------------------------------------------------
!
! Local variables.  Note that the arrays, while passed as arguments
! to the BMJ code, are not actually used by ARPS.  Also, the 3-D
! arrays are dimensioned i,k,j.
!
!------------------------------------------------------------------------

  REAL :: nca(nx,ny)
  REAL :: wrf_w0avg(nx,nz,ny)
  REAL :: wrf_rthcuten(nx,nz,ny)
  REAL :: wrf_rqvcuten(nx,nz,ny)
  REAL :: wrf_rqccuten(nx,nz,ny)
  REAL :: wrf_rqrcuten(nx,nz,ny)
  REAL :: wrf_rqicuten(nx,nz,ny)
  REAL :: wrf_rqscuten(nx,nz,ny)

  INTEGER :: ids,ide,jds,jde,kds,kde, &
             ims,ime,jms,jme,kms,kme, &
             its,ite,jts,jte,kts,kte

  REAL :: svp1,svp2,svp3,svpt0

  INTEGER :: p_qi,p_qs,p_first_scalar
  INTEGER :: i,j

  DATA svp1,svp2,svp3,svpt0 /0.6112,17.67,29.65,273.15/
  DATA p_qi,p_qs,p_first_scalar /3,4,1/

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
! Call WRF KF_ETA initializer.
!
!------------------------------------------------------------------------

  CALL kf_eta_init(wrf_rthcuten,wrf_rqvcuten,wrf_rqccuten,wrf_rqrcuten, &
               wrf_rqicuten,wrf_rqscuten,nca,wrf_w0avg,p_qi,p_qs,       &
               svp1,svp2,svp3,svpt0,                                    &
               p_first_scalar,restart,                                  &
               ids,ide,jds,jde,kds,kde,                                 &
               ims,ime,jms,jme,kms,kme,                                 &
               its,ite,jts,jte,kts,kte)

  DO i=its,ite
    DO j=jts,jte
      inca(i,j)=nint(nca(i,j))
    ENDDO
  ENDDO

!------------------------------------------------------------------------
!
! The end.
!
!------------------------------------------------------------------------

  RETURN
END SUBROUTINE interface_wrf_kfinit
