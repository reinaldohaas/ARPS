!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE QPFGRID                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE qpfgrid(nx,ny,nz,prcrate,pprt,ptprt,qv,pbar,ptbar,           &
           rhostr,zp,j3,j3inv,raing,temxy1,temxy2,pi,dqv)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate the source/sink terms in temperature and moisture
!  equations as a results of large-scale (as opposed to subgird scale)
!  condensation. Apply the source terms to ptprt and qv.
!
!  The surface gridscale rainfall is accumulated in raing.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: V. Wong, L. Zhao and X. Song
!  03/12/95
!
!  MODIFICATION HISTORY:
!
!  05/10/95 (L. Zhao)
!  Correction made to the calculation of accumulated rainfall
!  for restart runs.
!
!  8/12/95 (M. Xue)
!  Rearranged the argument list to follow ARPS convention.
!
!  01/31/1996 (V. Wong and X. Song)
!  Related the accumulated grid-scale precipitation to the parameter,
!  qpfgfrq, that controls the the frequency of calling QPFGRID.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!    pprt     Perturbation pressure (Pascal)
!    ptprt    Perturbation potential temperature at a given time
!             level (K)
!    qv       Water vapor specific humidity at a given time level
!             (kg/kg)
!    ptbar    Base state potential temperature (K)
!    pbar     Base state pressure (Pascal)
!
!  OUTPUT:
!
!    raing    Accumulated gridscale rainfall (mm).
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx,ny,nz          ! The number grid points in 3 directions

  REAL :: pprt   (nx,ny,nz) ! Perturbation pressure (Pascal)
  REAL :: ptprt  (nx,ny,nz) ! Perturbation potential temperature (K)
  REAL :: qv     (nx,ny,nz) ! Water vapor specific humidity (kg/kg)
  REAL :: pbar   (nx,ny,nz) ! Base state pressure (Pascal)
  REAL :: ptbar  (nx,ny,nz) ! Base state potential temperature (K)
  REAL :: rhostr (nx,ny,nz) ! Base state density rhobar times j3.

  REAL :: zp     (nx,ny,nz) ! The height of the terrain.
  REAL :: j3     (nx,ny,nz) ! Coordinate transformation
                            ! Jacobian  d(zp)/d(z)
  REAL :: j3inv  (nx,ny,nz) ! Coordinate transformation
                            ! Jacobian  d(zp)/d(z)


  REAL :: raing  (nx,ny)    ! Accumulated grid precipitation (mm)
  REAL :: prcrate(nx,ny)    ! precipitation rate (kg/(m**2*s))

  REAL :: temxy1 (nx,ny)    ! 2-D work array
  REAL :: temxy2 (nx,ny)    ! 2-D work array
  REAL :: pi     (nx,ny)    ! 2-D work array
  REAL :: dqv    (nx,ny)    ! 2-D work array
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i, j, k

  REAL :: rhmin,rhmax,al2ovcr,alovcp
  DATA rhmin/0.05/
!  DATA rhmax/1.0/
  DATA al2ovcr/1.346946E+07/              !Lv**2/(CP*Rv)
  DATA alovcp/2.487458E+03/               !Lv/Cp

  REAL :: es,temp2df,pres,tdf,theta
  REAL :: p0inv,pterm,tema,temb
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
  p0inv=1./p0
!
!-----------------------------------------------------------------------
!
!  Initializing the working array for each (i,j) point at beginning.
!
!-----------------------------------------------------------------------
!
  DO i=1,nx-1
    DO j=1,ny-1
      temxy1(i,j)=0.0
    END DO
  END DO
!
!-----------------------------------------------------------------------
!

  DO k=nz-2,2,-1
!
!-----------------------------------------------------------------------
!
    DO i=1,nx-1
      DO j=1,ny-1
!
!-----------------------------------------------------------------------
!
!  Computer needed variables from model inputs. p0 and rddcp are
!  included in 'physct.inc'
!
!-----------------------------------------------------------------------
!
        pres=pprt(i,j,k) + pbar(i,j,k)              ! p
        pi(i,j)=(pres*p0inv)**rddcp                    ! pi
        tdf=(ptbar(i,j,k)+ptprt(i,j,k))*pi(i,j)     ! temperature
!
!-----------------------------------------------------------------------
!
!  Calculate saturation vapor pressure(ES) at (P,T).
!
!-----------------------------------------------------------------------
!
        pterm = 0.5 + SIGN(0.5, tdf-260.0)

        tema = pterm*17.2694 + (1.-pterm)*21.87456
        temb = pterm*35.86   + (1.-pterm)*7.66

        es = 610.78*EXP(tema*(tdf-273.16)/(tdf-temb))
!
!-----------------------------------------------------------------------
!
!  Computer the saturation specific humidity.
!
!-----------------------------------------------------------------------
!
        dqv(i,j)=0.622*es/(pres-0.378*es)
!
!-----------------------------------------------------------------------
!
!  Set a lower limit to 'qv field' to avoid negative values.
!
!-----------------------------------------------------------------------
!
        qv(i,j,k) = MAX( qv(i,j,k), rhmin*dqv(i,j) )
!
!-----------------------------------------------------------------------
!
!  Set the upper limit of 'qv field' such that the RHMAX is equal to
!  or less than 1.0.
!
!  Modified by Zuwen He (04/2002) 
!  adopt new naming convention: rhsat
!
!  original code: dqv(i,j)=rhmax*dqv(i,j)
!  new code:      dqv(i,j)=rhsat*dqv(i,j)
!
!-----------------------------------------------------------------------
!
        dqv(i,j)=rhsat*dqv(i,j)
!
!-----------------------------------------------------------------------
!
!  Compute Dq/Dt term.
!
!-----------------------------------------------------------------------
!
        temp2df=1.0+al2ovcr*dqv(i,j)/(tdf*tdf)
!
!-----------------------------------------------------------------------
!
!    Check for condensation.
!
!-----------------------------------------------------------------------
!
        dqv(i,j)=(qv(i,j,k)-dqv(i,j))/temp2df
        IF (dqv(i,j) >= 0.0) temxy1(i,j)=temxy1(i,j)+dqv(i,j)           &
             *rhostr(i,j,k)*(zp(i,j,k+1)-zp(i,j,k))*j3inv(i,j,k)
!
!-----------------------------------------------------------------------
!
!  Check for evaporation.
!
!-----------------------------------------------------------------------
!
        temp2df=temxy1(i,j)+dqv(i,j)                                    &
             *rhostr(i,j,k)*(zp(i,j,k+1)-zp(i,j,k))*j3inv(i,j,k)
        IF(temp2df < 0) THEN
          temp2df=0.0
          dqv(i,j)=-temxy1(i,j)                                         &
              /(rhostr(i,j,k)*(zp(i,j,k+1)-zp(i,j,k))*j3inv(i,j,k))
        END IF
!
        IF(dqv(i,j) < 0)temxy1(i,j)=temp2df

        temxy2(i,j)=tdf

        temxy1(i,j)=MAX(0.0,temxy1(i,j))      !keep it positive

      END DO
    END DO
!
!-----------------------------------------------------------------------
!
!    Modify T and q fields at the Kth level.
!
!-----------------------------------------------------------------------
!
    IF (curtim > dtbig) THEN

      DO i=1,nx-1
        DO j=1,ny-1
          qv(i,j,k)=MAX(0.0,qv(i,j,k)-dqv(i,j))
          tdf=temxy2(i,j) + alovcp*dqv(i,j)
          theta=tdf/pi(i,j)
          ptprt(i,j,k)=theta-ptbar(i,j,k)
        END DO
      END DO
    END IF

  END DO

!
!-----------------------------------------------------------------------
!
!  Add in condensation to stable preiciptation array (accumulated
!  grid precipitation). The value is cumulated until the cumulus
!  parameterization is turn on in the mode, which is control by
!  'confrg'. The unit for rainfall is 'mm'.
!
!-----------------------------------------------------------------------
!
!   IF (mod(curtim+0.001,confrq).lt.(0.5*dtbig)) THEN
  DO i=1,nx-1
    DO j=1,ny-1
      tema=temxy1(i,j)/rhow
!       raing(i,j)=raing(i,j) +
!  :               1000.0* temxy1(i,j) * ( pprt(i,j,2)+pbar(i,j,2) -
!  :               pprt(i,j,3)-pbar(i,j,3) ) /(rhow*g)
      raing(i,j)=raing(i,j) + 1000.0*tema
      prcrate(i,j)=tema*rhow/qpfgfrq
    END DO
  END DO
!   ENDIF

  RETURN
END SUBROUTINE qpfgrid
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE QPFCUMS                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE qpfcums(nx,ny,nz,prcrate,qvsflx,                             &
           u,v,w,pprt,ptprt,qv,pbar,ptbar,qvbar,rhostr,zp,j3,           &
           ptcumsrc,qcumsrc,rainc,nca,kfraincv,                         &
           cldefi,xland,bmjraincv,                                      &
           tem1,tem2,tem3,tem4,tem5,tem6,tem7,tem8,tem9,tem10)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate the source/sink terms in temperature and moisture
!  equations as a results of subgrid scale (as opposed to grid scale)
!  cumulus convection.
!
!  Kuo, Kain-Fritsch, and WRF BMJ schemes are used here.
!
!  The surface convective rainfall is accumulated in rainc.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: V. Wong, L. Zhao and X. Song
!  3/06/95
!
!  MODIFICATION HISTORY:
!
!  8/12/95 (M. Xue)
!  Rearranged the argument list to follow ARPS conventions.
!
!  08/01/97 (Zonghui Huo)
!  Added Kain-fritsch cumulus parameterization scheme.
!
!  13 March 2002 (Eric Kemp)
!  Added WRF BMJ cumulus parameterization scheme.
!
!  April 2002 (Fanyou Kong)
!  Added WRF new Kain-Fritsch (April 2002 version: KF_ETA) scheme
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!    u        x component of velocity at a given time level (m/s)
!    v        y component of velocity at a given time level (m/s)
!    w        * z component of velocity at a given time level (m/s)
!               for KUO scheme
!             * time-average vertical velocity (m/s)
!               for Kain Fritsch scheme
!
!    pprt     Perturbation pressure (Pascal)
!    ptprt    Perturbation potential temperature (K)
!    qv       Water vapor specific humidity (kg/kg)
!    pbar     Base state pressure (Pascal)
!    ptbar    Base state potential temperature (K)
!    qvbar    Base state water vapor specific humidity (kg/kg)
!    rhostr   Base state density rhobar times j3 (kg/m**3)
!
!    zp       The physical height coordinate defined at w-point
!    j3       Coordinate transformation Jacobian, d(zp)/d(z)
!
!  OUTPUT:
!
!    ptcumsrc Potential temperature source term.
!    qcumsrc  Water specific humidity source term
!    rainc    Accumulated precipitation by cumulus convection (mm)
!    prcrate  precipitatioon rate (mm/s)
!    kfraincv K-F convective rainfall (cm)
!    nca      K-F counter for CAPE release
!    cldefi   BMJ cloud efficiency
!    xland    BMJ land/sea mask
!    bmjraincv   BMJ convective rainfall (cm)
!
!  WORK ARRAYS:
!
!    tem1     Temporary work array
!    tem2     Temporary work array
!    tem3     Temporary work array
!    tem4     Temporary work array
!    tem5     Temporary work array
!    tem6     Temporary work array
!    tem7     Temporary work array
!    tem8     Temporary work array
!    tem9     Temporary work array
!    tem10    Temporary work array
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx,ny,nz          ! The number grid points in 3 directions
  REAL :: u      (nx,ny,nz) ! Total u-velocity (m/s)
  REAL :: v      (nx,ny,nz) ! Total v-velocity (m/s)
  REAL :: w      (nx,ny,nz) ! Total w-velocity (m/s)

  REAL :: pprt   (nx,ny,nz) ! Perturbation pressure (Pascal)
  REAL :: ptprt  (nx,ny,nz) ! Perturbation potential temperature (K)
  REAL :: qv     (nx,ny,nz) ! Water vapor specific humidity (kg/kg)
  REAL :: pbar   (nx,ny,nz) ! Base state pressure (Pascal)
  REAL :: ptbar  (nx,ny,nz) ! Base state potential temperature (K)
  REAL :: qvbar  (nx,ny,nz) ! Base state water vapor specific
  REAL :: rhostr (nx,ny,nz) ! Base state density rhobar times j3.

  REAL :: zp     (nx,ny,nz) ! The height of the terrain.
  REAL :: j3     (nx,ny,nz) ! Coordinate transformation
                            ! Jacobian  d(zp)/d(z)

  REAL :: ptcumsrc(nx,ny,nz) ! Source term in pt equation.

  REAL :: qcumsrc(nx,ny,nz,5)  ! Source term in water equations due
                               ! to cumulus parameterization:
                               ! qcumsrc(1,1,1,1) for qv equation
                               ! qcumsrc(1,1,1,2) for qc equation
                               ! qcumsrc(1,1,1,3) for qr equation
                               ! qcumsrc(1,1,1,4) for qi equation
                               ! qcumsrc(1,1,1,5) for qs equation
  REAL :: rainc(nx,ny)       ! Accumulated precipitation by convection
  REAL :: prcrate(nx,ny)     ! precipitation rate (kg/(m**2*s))
  REAL :: qvsflx(nx,ny)      ! Surface moisture flux (kg/(m**2*s))

  REAL :: kfraincv    (nx,ny)   ! K-F convective rainfall (cm)
  INTEGER :: nca       (nx,ny)  ! K-F counter for CAPE release

!EMK BMJ
  REAL,INTENT(INOUT) :: cldefi(nx,ny)    ! BMJ cloud efficiency
  REAL,INTENT(IN) :: xland(nx,ny)        ! BMJ land mask
                                         ! (1.0 = land, 2.0 = sea)
  REAL,INTENT(INOUT) :: bmjraincv(nx,ny) ! BMJ convective rainfall (cm)
!EMK END

  REAL :: tem1 (nx,ny,nz)    ! Temporary work array
  REAL :: tem2 (nx,ny,nz)    ! Temporary work array
  REAL :: tem3 (nx,ny,nz)    ! Temporary work array
  REAL :: tem4 (nx,ny,nz)    ! Temporary work array
  REAL :: tem5 (nx,ny,nz)    ! Temporary work array
  REAL :: tem6 (nx,ny,nz)    ! Temporary work array
  REAL :: tem7 (nx,ny,nz)    ! Temporary work array
  REAL :: tem8 (nx,ny,nz)    ! Temporary work array
  REAL :: tem9 (nx,ny,nz)    ! Temporary work array
  REAL :: tem10(nx,ny,nz)    ! Temporary work array
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i, j, k
  REAL :: p0inv,dthmax
  REAL :: kqmax,jqmax,iqmax
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'phycst.inc'
  INCLUDE 'cumucst.inc'
  INCLUDE 'mp.inc'
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  p0inv=1./p0

  CALL set_acct(cum_acct)

  IF(cnvctopt == 3) THEN


    IF (myproc == 0) &
      WRITE(6,'(1x,a)') 'Use Kain-Fritcsh Scheme,Enter kfinterfc'

    CALL  kfinterfc(nx,ny,nz,u,v,w,pprt,ptprt,qv,pbar,ptbar,zp,         &
                    ptcumsrc,qcumsrc,rainc,prcrate,                     &
                    nca,kfraincv,                                       &
                    tem1(1,1,1),tem1(1,1,2),tem1(1,1,3),                &
                    tem2,tem3,tem4,tem5,tem6,tem7,tem8,tem9,tem10)

  ELSE IF (cnvctopt == 5) THEN

    IF (myproc == 0) &
      WRITE(6,'(1x,a)') 'Use WRF Kain-Fritcsh Scheme'

    CALL interface_wrf_kfetadrv(nx,ny,nz,u,v,w,                         &
                                pprt,ptprt,qv,pbar,ptbar,zp,            &
                                ptcumsrc,qcumsrc,prcrate,               &
                                nca,kfraincv)

  ELSE IF (cnvctopt == 4) THEN
  
    IF (myproc == 0) &
      WRITE(6,'(1x,a)') 'Use WRF Betts-Miller-Janjic Scheme'

    CALL interface_wrf_bmjdrv(nx,ny,nz,pprt,ptprt,qv,pbar,ptbar,zp,     &
                              ptcumsrc,qcumsrc,bmjraincv,prcrate,cldefi,&
                              xland)
    
  ELSE         ! use Kuo-scheme

    IF (myproc == 0) &
      WRITE(6,'(1x,a)') 'Use KUO cumulus scheme'
!
!-----------------------------------------------------------------------
!
!    For each grid point, add the convective parameterization.
!
!-----------------------------------------------------------------------
!
    dthmax=0.0
    DO j=1,ny-1
      DO i=1,nx-1
        msflux=-qvsflx(i,j)
        DO k=1,nz-1
          ucon(k)=(u(i,j,k)+u(i+1,j,k))/2.0
          vcon(k)=(v(i,j,k)+v(i,j+1,k))/2.0
          wcon(k)=w(i,j,k)
          picon(k)=cp*((pprt(i,j,k) + pbar(i,j,k))*p0inv)**rddcp
          thtcon(k)=ptbar(i,j,k)+ptprt(i,j,k)
          tmpcon(k)=thtcon(k)*picon(k)/cp
!
! DNCON = rhobar
!
          dncon(k)=rhostr(i,j,k)
          qvcon(k)=qv(i,j,k)

          ptcumsrc(i,j,k)=0.    !initializing forcing terms
          qcumsrc(i,j,k,1)=0.
!
          zzcon(k) =zp(i,j,k)
          zcon(k)=0.5*(zp(i,j,k)+zp(i,j,k+1))

        END DO
!
!-----------------------------------------------------------------------
!
!  Calculate enviromental parameters.
!
!-----------------------------------------------------------------------
!
        CALL environc(nz-1)
!
!-----------------------------------------------------------------------
!
!  IF there is convective heating, Kuo's cumulus parameterization
!  scheme is used.
!
!-----------------------------------------------------------------------
!
        IF(igo /= 0) THEN
          CALL kuocp
        END IF
!
!-----------------------------------------------------------------------
!
!    If there is convective heating, adjust and output the forcnig
!    terms into the model.
!
!-----------------------------------------------------------------------
!
        IF(igo /= 0) THEN
          CALL cp2mod (nz-1)
!
          DO k=2,nz-2
            ptcumsrc(i,j,k)=ftcon(k)
            qcumsrc(i,j,k,1)=frcon(k)
          END DO

!       IF (mod(curtim+0.001,confrq).lt.(0.5*dtbig))THEN
          rainc(i,j)=rainc(i,j)+1000.0*cprecip*confrq/rhow
          prcrate(i,j)=cprecip
!       ENDIF
!
!-----------------------------------------------------------------------
!
!   Find the location of maximum convective heating rate.
!
!-----------------------------------------------------------------------
!
          DO k=2,nz-2
            IF(ptcumsrc(i,j,k) > dthmax) THEN
              dthmax=ptcumsrc(i,j,k)
              iqmax=i
              jqmax=j
              kqmax=k
            END IF
          END DO

        END IF

      END DO
    END DO

  END IF           ! End of cumulus schemes

  RETURN
END SUBROUTINE qpfcums
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE KUOCP                      ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE kuocp
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Performs the Kuo's convective cumulus parameterization.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: V. Wong, L. Zhao and X. Song
!  3/06/95
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: k
  INTEGER :: kdiv,klfs,kover,kdet
  REAL :: supply
!  real SUPPLYW

  REAL :: zdetr,bkuo,dzdet,wtgnd,wtdiv,preff,envshr,vhint,overmax,      &
       factr,vdint,vmint,avtdiff,avgmin,wtlcl,                          &
       anegl,dzdiv,dddt,wtlfs,anegh,                                    &
       apos
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'phycst.inc'
  INCLUDE 'cumucst.inc'
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!-----------------------------------------------------------------------
!
!      Downdraft flag - 0 - no downdrafts
!                       1 - simple downdraft model
!
!-----------------------------------------------------------------------
!
!    Initializing the data array for theta and qv forcing terms.
!
!-----------------------------------------------------------------------
!
  DO k=1,nkp
    ftcon(k)=0.
    frcon(k)=0.
  END DO
!
!-----------------------------------------------------------------------
!
!    Compute vertical moisture convergence into the cloud layer.
!    Neglect the horizontal advection of specific humidity since
!    it can be shown to be relatively small. Vertical flux out cloud
!    top is also assumed small.
!
!-----------------------------------------------------------------------
!
!  SUPPLYW=RHOE(KLCL)*RVE(KLCL)*(WPE(KLCL)+WPE(KLCL-1))*.5
!
!  SUPPLY=MSFLUX-0.125*(RHOE(2)+RHOE(1))
!    :                   *(WPE(2)+WPE(1))*(RVE(2)-RVE(1))

  supply=msflux-0.125*(rhoe(klcl+1)+rhoe(klcl))                         &
                     *(wpe(klcl+1)+wpe(klcl))                           &
                     *(rve(klcl+1)-rve(klcl))

!  DO 105 K=3,KMT-1
  DO k=klcl+2,kmt-1
    supply=supply-0.25*(rhoe(k)+rhoe(k-1))                              &
                      *(wpe(k)+wpe(k-1))*(rve(k)-rve(k-1))
  END DO

  supply=supply-0.125*(rhoe(kmt)+rhoe(kmt-1))                           &
                     *(wpe(kmt)+wpe(kmt-1))*(rve(kmt)-rve(kmt-1))

!  SUPPLY=SUPPLYW
!
!-----------------------------------------------------------------------
!
!    If water supply is not large enough, no convective is allowed.
!
!-----------------------------------------------------------------------
!
  IF (supply <= 1.0E-3) THEN
    igo=0
    RETURN
  END IF
!
!-----------------------------------------------------------------------
!
!    This is the cloud model.  Updraft is constant THETA e and
!    saturated with respect to water.  There is no ice.
!    Cloud top is one level above ETL.
!
!    THETA e of the updraft
!
!-----------------------------------------------------------------------
!
  theu(klcl)=(3.376/tlcl-0.00254)*1000*rve(kcon)                        &
      *(1+0.81*rve(kcon))
  theu(klcl)=the(kcon)*EXP(theu(klcl))
  PRINT*,'THEU(KLCL)=',theu(klcl),'  KLCL=',klcl
!
!-----------------------------------------------------------------------
!
!    Equilibrium Temperature Level of the source level air.
!
!-----------------------------------------------------------------------
!
  igo=0                            ! index variable for convection
  DO k=klcl,kmt
    CALL the2t(theu(klcl),pe(k),thu(k),tu(k),rsu(k))
    IF(thu(k) > the(k).AND.igo == 0) THEN
      igo=1
      klfc=k                       ! level of free convection
    END IF
    IF(thu(k) <= the(k).AND.igo == 1) GO TO 500           !cloud top
  END DO

  IF (igo == 0) RETURN

  PRINT*,' Convection beyond model top - THup, THenv ',thu(kmt)         &
        ,the(kmt)
  k=kmt-1

  500   CONTINUE

!
!  KETL=MIN(K,KMT)
!  KCT=MIN(KETL+1,KMT)
!
  ketl=MIN(k-1,kmt)
  kct=MIN(ketl,kmt)

  CALL the2t(theu(klcl),pe(kct),thu(kct),tu(kct),rsu(kct))

  DO k=1,klfc-1
    thu(k)=the(k)
  END DO
!
!-----------------------------------------------------------------------
!
!   If the cloud is not at least CDZMIN deep or cloud top is
!   under 500 mb, no convection is allowed.
!
!-----------------------------------------------------------------------
!
  IF (ze(ketl)-ze(klfc) < cdzmin.OR.pe(kct) > 50000.) THEN
    igo=0
    RETURN
  END IF
!
!-----------------------------------------------------------------------
!
!  Require the positive area be 50% greater than the negative
!  area below the LFC and 5% greater in total.
!
!-----------------------------------------------------------------------
!
  anegl=0.                        ! negative area
  DO k=klcl,klfc-1
    anegl=anegl+(thu(k)-the(k))*(zc(k)-zc(k-1))
  END DO
!
  apos=0.                         ! positive area
  DO k=klfc,ketl-1
    apos=apos+(thu(k)-the(k))*(zc(k)-zc(k-1))
  END DO
!
  anegh=0.
  DO k=ketl,kct
    anegh=anegh+(thu(k)-the(k))*(zc(k)-zc(k-1))
  END DO
!
  IF(apos < ABS(anegl)*1.5.OR.apos < ABS(anegl+anegh)*1.05) THEN
    igo=0
    RETURN
  END IF
!
!-----------------------------------------------------------------------
!
!   The downdraft model - starts at THETA e minimum (LFS).
!          Downdraft is 2 degrees colder than
!          environment at cloud base increasing to 5 degrees
!          colder at the ground.
!
!   Find LFS as THETA e minimum
!
!-----------------------------------------------------------------------
!
  IF (idownd == 1) THEN
!
    DO k=kct,2,-1
      IF (thee(k) < thee(k+1).AND.thee(k) < thee(k-1)) GO TO 510
    END DO
    k=2
    510     CONTINUE
    klfs=k
    IF(klfs <= klcl)klfs=klcl+1
    thd(klfs)=the(klfs)
!
!-----------------------------------------------------------------------
!
!     Limit dd deficit at the ground to the maximum of positive
!       temperature difference of updraft if less than 2.5 degrees.
!
!-----------------------------------------------------------------------
!
    dddt=0.
    DO k=klcl,kct
      dddt=MAX(dddt,thu(k)-the(k))
    END DO
    IF(dddt > 2.5) dddt=5.
!
    thd(2)=the(2)-dddt
    thd(klcl)=the(klcl)-dddt*.2
    DO k=klcl,klfs
      thd(k)=thd(klcl)+(thd(klfs)-thd(klcl))/(ze(klfs)-ze(klcl))        &
             *(ze(k)-ze(klcl))
    END DO
    DO k=3,klcl-1
      thd(k)=thd(2)+(thd(klcl)-thd(2))/(ze(klcl)-ze(2))                 &
              *(ze(k)-ze(2))
    END DO
!
!-----------------------------------------------------------------------
!
!    Now we need to weight the downdraft relative to the updraft.
!    Assume that the dd weight is zero at the LFS, 1/2 of
!    updraft at cloud base, and equal to the updraft at cloud
!    base at the ground.
!
!-----------------------------------------------------------------------
!
    dzdiv=1E20
    DO k=1,kmt
      IF(ABS(ze(k)-800.) < dzdiv)THEN
        kdiv=k
        dzdiv=ABS(ze(k)-800.)
      END IF
    END DO
    kdiv=MAX(MIN(klcl,kdiv),2)
    IF(kdiv == klcl) kdiv=klcl-1

    DO k=1,nkp
      wtd(k)=0.
    END DO

    wtlfs=0.
    wtlcl=.1
    wtdiv=.2
    wtgnd=1.
    DO k=klcl+1,klfs
      wtd(k)=wtlcl+(wtlfs-wtlcl)/(ze(klfs)-ze(klcl))                    &
          *(ze(k)-ze(klcl))
    END DO

    DO k=kdiv,klcl
      wtd(k)=wtdiv+(wtlcl-wtdiv)/(ze(klcl)-ze(kdiv))                    &
          *(ze(k)-ze(kdiv))
    END DO

    DO k=2,kdiv-1
      wtd(k)=wtgnd+(wtdiv-wtgnd)/(ze(kdiv)-ze(2))                       &
             *(ze(k)-ze(2))
    END DO

  ELSE                           ! for idownd=0 option

    DO k=1,nkp               ! set the downdraft weights to zero
      wtd(k)=0.
    END DO

    DO k=2,klcl-1            ! below cloud base, no difference
      thu(k)=the(k)
    END DO

  END IF
!
!-----------------------------------------------------------------------
!
!    Compute the b parameter. Use Fritsch/Chappell's precipitation
!    efficiency. Notice that the unit of (DV/DZ) is 10e-3 s-1.
!
!-----------------------------------------------------------------------
!
  envshr=SQRT((upe(kct)-upe(klfc))**2                                   &
             +(vpe(kct)-vpe(klfc))**2)                                  &
             /(ze(kct)-ze(klfc))*1E3
  IF (envshr > 1.35) THEN
    preff=1.591-.639*envshr+.0953*envshr**2-.00496*envshr**3
  ELSE
    preff=.9                   !parameter (1-b)
  END IF

  preff=MAX(0.0,preff)         !keep it positive

  bkuo=1.-preff                !parameter b
!
!-----------------------------------------------------------------------
!
!    Initializing the vertical profiles of convective heating and
!    moistening
!
!-----------------------------------------------------------------------
!
  DO k=2,kmt
    vheat(k)=0.
    vmois(k)=0.
    vmdry(k)=0.
  END DO
!
!-----------------------------------------------------------------------
!
!   Find the weighted THETA to use for the convection. Where WTD is the
!   weights for downdraft. When WTD(k)=0, THCON(k)=THU(k).
!
!-----------------------------------------------------------------------
!
  DO k=2,kct
    thcon(k)=wtd(k)*thd(k)+(1.-wtd(k))*thu(k)
  END DO
!
!-----------------------------------------------------------------------
!
!    Heating profile is difference between convective THETAs and
!    environment.
!
!-----------------------------------------------------------------------
!
  DO k=2,kct
    vheat(k)=thcon(k)-the(k)                !delta theta
  END DO
!
!-----------------------------------------------------------------------
!
!    Moisture profile is difference between vapor's of updraft and
!    environment in the cloud layer.  Below cloud base, air is
!    dried by SUPPLY.  Downdrafts are assumed to have no effect
!    on this.
!
!-----------------------------------------------------------------------
!
  zdetr=.66667*ze(kct)
  dzdet=1000000.
  DO k=klcl,kct
    IF(ABS(ze(k)-zdetr) < dzdet)THEN
      dzdet=ABS(ze(k)-zdetr)
      kdet=k
    END IF
  END DO

  DO k=kdet,kct
    vmois(k)=1.
  END DO

  DO k=klcl,kdet-1
    vmois(k)=rsu(k)-rve(k)
  END DO

  DO k=2,klcl-1
    vmdry(k)=rve(k)
  END DO

  vhint=0.
  vmint=0.
  vdint=0.
  DO k=2,kmt
    vhint=vhint+vheat(k)*(zc(k)-zc(k-1)) ! sum (dtheta*dz)
    vmint=vmint+vmois(k)*(zc(k)-zc(k-1)) ! sum(dq*dz) for moist case
    vdint=vdint+vmdry(k)*(zc(k)-zc(k-1)) ! sum(dq*dz) for dry case
  END DO
!
!-----------------------------------------------------------------------
!
!   If VHINT is less than 0, there is more negative area than
!   positive area.  No convection allowed.
!
!-----------------------------------------------------------------------
!
  IF (vhint <= 0.) THEN
    igo=0
    RETURN
  END IF
!
!-----------------------------------------------------------------------
!
!   Also require that there is a minimum average
!   temperature difference between the updraft and environment
!   from the LFC to the ETL.  This eliminates the cases where
!   VHINT is very small and the heating and cooling rates get
!   astronomically large.
!
!-----------------------------------------------------------------------
!
  avgmin=.10
  avtdiff=0.
  DO k=klfc,ketl-1
    avtdiff=avtdiff+(thcon(k)-the(k))
  END DO
  avtdiff=avtdiff/MAX(1,ketl-klfc)
  IF (avtdiff < avgmin) THEN
    igo=0
    RETURN
  END IF
!
!-----------------------------------------------------------------------
!
!    Heating and moistening rates
!
!-----------------------------------------------------------------------
!
  3100  CONTINUE
  DO k=2,kmt
    ftcon(k)=alvl*preff*supply*vheat(k)                                 &
            /(pke(k)*rhoe(k)*vhint)
  END DO
  DO k=klcl,kct
    frcon(k)=bkuo*supply*vmois(k)/(rhoe(k)*vmint)
  END DO
  DO k=2,klcl-1
!     FRCON(K)=-SUPPLY*VMDRY(K)/(RHOE(K)*VDINT)
    frcon(k)=0.0
  END DO
!
!-----------------------------------------------------------------------
!
  DO k=klfc,ketl-1
    qvct1(k)=the(k)+confrq*ftcon(k)
  END DO

  overmax=0.
  DO k=klfc,ketl-1
    IF(qvct1(k)-thu(k) > overmax)THEN
      overmax=(qvct1(k)-thu(k))/(ftcon(k)*confrq)
      kover=k
    END IF
  END DO

  IF(overmax > 0.) THEN
    factr=1.-overmax
    supply=factr*supply
    GO TO 3100
  END IF

  cprecip=preff*supply

  RETURN
END SUBROUTINE kuocp
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE CP2MOD                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE cp2mod(nzz)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Convert the cumulus heating and moistening in the convective grid
!  to model grid.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: V. Wong, L. Zhao and X. Song
!  3/06/95
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nzz
!    See the include file cumucst.inc
!
!  OUTPUT:
!
!    See the include file cumucst.inc
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nzz
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: k
  INTEGER :: nz1,l
  REAL :: tftm,frres,ftres,tfrm,tftc,s_sum,dzlft,tfrc
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'phycst.inc'
  INCLUDE 'cumucst.inc'

  REAL :: vctr5(nkp),vctr6(nkp)
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!-----------------------------------------------------------------------
!
!   Compute integrated heating and moistening energy tendencies
!
!-----------------------------------------------------------------------
!

  DO k=2,kmt
    qvct1(k)=rhoe(k)*ftcon(k)*pke(k)
    qvct2(k)=rhoe(k)*alvl*frcon(k)
    qvct3(k)=(zc(k)-zc(k-1))*qvct1(k)      !delta (ZE(k))
    qvct4(k)=(zc(k)-zc(k-1))*qvct2(k)
  END DO

  tftc=s_sum(kmt-1,qvct3(2),1)             !vertical integration
  tfrc=s_sum(kmt-1,qvct4(2),1)
!
!-----------------------------------------------------------------------
!
!   Transfer energy tendencies to model grids
!
!-----------------------------------------------------------------------
!
  nz1=nzz-1
  DO k=1,nzz                           !initializing data arrays
    vctr5(k)=0.
    vctr6(k)=0.
  END DO

  dzlft=0.
  l=2
  DO k=2,nzz-1
    IF(dzlft /= 0.) THEN
      vctr5(k)=vctr5(k)+qvct1(l)*dzlft
      vctr6(k)=vctr6(k)+qvct2(l)*dzlft
      l=l+1
    END IF
    500     CONTINUE
    IF(zc(l) <= zzcon(k)) THEN
      vctr5(k)=vctr5(k)+qvct1(l)*(zc(l)-zc(l-1))
      vctr6(k)=vctr6(k)+qvct2(l)*(zc(l)-zc(l-1))
      l=l+1
      dzlft=0.
      GO TO 500
    ELSE
      vctr5(k)=vctr5(k)+qvct1(l)*(zzcon(k)-zc(l-1))
      vctr6(k)=vctr6(k)+qvct2(l)*(zzcon(k)-zc(l-1))
      dzlft=zc(l)-zzcon(k)
    END IF
  END DO

  tftm=s_sum(nz1,vctr5(2),1)
  tfrm=s_sum(nz1,vctr6(2),1)
!
!-----------------------------------------------------------------------
!
!  Make sure the transfer from the convective grid to the model
!  grid happened correctly.
!
!-----------------------------------------------------------------------
!
  ftres=tftm-tftc
  frres=tfrm-tfrc
  IF(ABS(ftres) > .01*ABS(tftc)) THEN
    PRINT*,' Energy error in grid tranfser in convective param.'
    PRINT*,' TFTM,TFTC ',tftm,tftc
  END IF
!
!-----------------------------------------------------------------------
!
!  Change energy tendencies to temperature and mixing ratio
!  tendencies.
!
!-----------------------------------------------------------------------
!
  DO k=2,nzz-1
    ftcon(k)=vctr5(k)/((zzcon(k+1)-zzcon(k))*dncon(k)*picon(k))
    frcon(k)=vctr6(k)/((zzcon(k+1)-zzcon(k))*dncon(k)*alvl)
  END DO

  RETURN
END SUBROUTINE cp2mod
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE ENVIRONC                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE environc(nzz)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!     Stability check and determination of cloud base and top
!
!-----------------------------------------------------------------------
!
!  AUTHOR: V. Wong, L. Zhao and X. Song
!  3/06/95
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nzz
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: k,nkmid
  REAL :: rlll,zlll,plll,themax,tlll,dzlll,tdu,rdsu,thdu,dzdd,          &
       abe,znz,wcpmax
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'phycst.inc'
  INCLUDE 'cumucst.inc'
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  REAL :: hz(nkp)
!
!-----------------------------------------------------------------------
!
! Basic constants
!
!-----------------------------------------------------------------------
!
  cpr=3.4965                ! Cp/Rd
  alvl=2.68E6               ! Binbin Zhuo, July 24, 1997
  aliv=2.837E6              ! Lf
  dzlow=200.                ! higher vertical increment
  dzhigh=500.               ! lower vertical increment
  zmid=3000.                ! height of cloud base
  cdzmin=3000.              ! depth of cloud
!
!-----------------------------------------------------------------------
!
!   Compute moist static energy profile. Cp*T+G*Z+Lv*q
!
!-----------------------------------------------------------------------
!
  DO k=1,nzz
    hz(k)=cp*tmpcon(k)+g*zcon(k)+alvl*qvcon(k)
  END DO
!
!-----------------------------------------------------------------------
!
!   Check for conditional instability. If no instability exist, there
!   will be no convection.
!
!-----------------------------------------------------------------------
!
  igo=0                        !index variable for convection
  DO k=2,nzz
    IF(hz(k) > hz(k+1))THEN
      igo=1
      EXIT
    END IF
  END DO

!  500   CONTINUE

  IF(igo == 0) RETURN
!
!-----------------------------------------------------------------------
!
!   Check the upward motion. Only it is greater than some value, like
!   WCONMIN, under ZMID, the convection is possible.
!
!-----------------------------------------------------------------------
!
  igo=0
  wcpmax=-1.e10
  DO k=2,nzz
    IF(zcon(k) > zmid)EXIT
    wcpmax=MAX(wcpmax,wcon(k))        !w in cloud base
  END DO

!  510   CONTINUE

  IF (wcpmax > 0.0.AND.wcpmax > wcldbs ) igo=1  ! WCLDBS is input parm.
  IF (igo == 0) RETURN
!
!-----------------------------------------------------------------------
!
!    If the above two conditions are satisfied, set the convective analysis
!    grids. There are two type of resolutions. Lower resolution for upper
!    levels (above ZMID), higher resolution for lower levels (below ZMID)
!
!-----------------------------------------------------------------------
!
  nkmid=zmid/dzlow+1        !# of grid points below cloud base

  zc(1)=0.0
  DO k=2,nkmid
    zc(k)=zc(1)+(k-1)*dzlow
  END DO

  DO k=nkmid+1,nkp
    zc(k)=zc(nkmid)+(k-nkmid)*dzhigh
  END DO
!
!-----------------------------------------------------------------------
!
!    Computer the height at PI grid points.
!
!-----------------------------------------------------------------------
!
  ze(1)=0.
  DO k=2,nkp
    ze(k)=(zc(k)+zc(k-1))*.5
  END DO
!
!-----------------------------------------------------------------------
!
!    Find the model top on the convective analysis grids.
!
!-----------------------------------------------------------------------
!
  znz=zcon(nzz)      ! the height of the model top
  DO k=nkp,1,-1  ! to find the cloud grid close to the model top
    IF(ze(k) < znz)GO TO 520
  END DO

  CALL arpsstop('arpsstop called from environc',1)
         ! means cloud is out of the analysis domain

  520   CONTINUE

  kmt=k                  ! the level of cloud top in cloud grids
!
!-----------------------------------------------------------------------
!
!    Interpolate model(environment) variables to the convective analysis
!    grids.
!
!-----------------------------------------------------------------------
!
  CALL htint(nzz,ucon,zcon,kmt,upe,ze)        ! u
  CALL htint(nzz,vcon,zcon,kmt,vpe,ze)        ! v
  CALL htint(nzz,wcon,zzcon,kmt,wpe,ze)       ! w
  CALL htint(nzz,thtcon,zcon,kmt,the,ze)      ! theta
  CALL htint(nzz,qvcon,zcon,kmt,rve,ze)       ! q
!
!-----------------------------------------------------------------------
!
!   Set the minimum value of qv.
!
!-----------------------------------------------------------------------
!
  DO k=1,kmt
    rve(k)=MAX(rve(k),1E-8)          ! set the limited qv
  END DO
!
!-----------------------------------------------------------------------
!
!    Computer theta V, theta E, and get PI pressure profile.
!
!-----------------------------------------------------------------------
!
  DO k=1,kmt
    thve(k)=the(k)*(1.+.61*rve(k))   ! virtual potential temperature
  END DO

  pke(1)=picon(1)
  DO k=2,kmt
    pke(k)=pke(k-1)-g*2.*(ze(k)-ze(k-1))                                &
          /(thve(k)+thve(k-1))
  END DO

  DO k=1,kmt
    te(k)=the(k)*pke(k)/cp                          ! temp
    pe(k)=(pke(k)/cp)**cpr*p0                       ! pressure
    rhoe(k)=pe(k)/(rd*te(k)*(1.+.61*rve(k)))        ! density
  END DO
!
!-----------------------------------------------------------------------
!
!   Computer theta E.
!
!-----------------------------------------------------------------------
!
  DO k=1,kmt
    CALL thetae(pe(k),te(k),rve(k),thee(k))
  END DO
!
!-----------------------------------------------------------------------
!
!    Find the main source level of the updraft. First, test if any
!    inversion below 1.2 km.
!
!-----------------------------------------------------------------------
!
  DO k=3,nkmid                !check the levels below ZMID
    IF(te(k) > te(k-1).AND.te(k) > te(k+1)                              &
                         .AND.ze(k) <= 1200.)THEN
      kcon=k
      GO TO 530
    END IF
  END DO
!
!-----------------------------------------------------------------------
!
!    If there isn't an inversion, use the level of highest theta E.
!
!-----------------------------------------------------------------------
!
  themax=0.
  DO k=2,nkmid
    IF(thee(k) > themax)THEN
      themax=thee(k)
      kcon=k
    END IF
  END DO
!
!-----------------------------------------------------------------------
!
!    Find the LCL of a layer average around the source level.
!
!-----------------------------------------------------------------------
!

  530   CONTINUE

  tlll=(te(kcon)+te(kcon+1)+te(kcon-1))/3.      ! averaged
  plll=pe(kcon)
  rlll=(rve(kcon)+rve(kcon+1)+rve(kcon-1))/3.   ! averaged
  zlll=ze(kcon)
!
!-----------------------------------------------------------------------
!
!   Find the (p,t) at the condesation level and the height of LCL relative
!   to the source level.
!
!-----------------------------------------------------------------------
!
  CALL lcl(tlll,plll,rlll,tlcl,plcl,dzlcl)      ! condesation level
!
!-----------------------------------------------------------------------
!
!    Find the closest level on the convective grid near the LCL.
!
!-----------------------------------------------------------------------
!
  dzlll=1E20
  DO k=1,kmt
    dzdd=ABS(ze(k)-(zlll+dzlcl))
    IF(dzdd < dzlll)THEN
      dzlll=dzdd
      klcl=k
    END IF
  END DO
!
!-----------------------------------------------------------------------
!
!    If there is not upward motion at the LCL, no convection (must be
!    greater than wconmin)
!
!-----------------------------------------------------------------------
!
  IF(wpe(klcl) < 0.0.OR.wpe(klcl) < wcldbs)THEN
    igo=0
    RETURN
  END IF
!
!-----------------------------------------------------------------------
!
!    Locate equilibrium temperature level of an unentrained parcel.
!
!-----------------------------------------------------------------------
!
  theu(klcl)=(3.376/tlcl-0.00254)*1000*rve(kcon)                        &
      *(1+0.81*rve(kcon))
  theu(klcl)=the(kcon)*EXP(theu(klcl))

  DO k=klcl,kmt
    IF(theu(klcl) <= thve(k))THEN
      ketl=k
      GO TO 540
    END IF
  END DO

  CALL arpsstop("arpsstop called from environc at 540",1)

  540   CONTINUE
!
!-----------------------------------------------------------------------
!
!   If the cloud depth is less than a critical value (CDZMID), then
!   convection is not allowed.
!
!-----------------------------------------------------------------------
!
  IF(ze(ketl)-ze(klcl) < cdzmin)THEN
    igo=0
    RETURN
  END IF
!
!-----------------------------------------------------------------------
!
!    Computer initial ABE (Avaliable Buoyant Energy). Convection is allowed
!    only if ABE is positive.
!
!-----------------------------------------------------------------------
!
  abe=0.
  DO k=klcl,ketl
    CALL the2t(theu(klcl),pe(k),thdu,tdu,rdsu)
    abe=abe+(thdu*(1.+.61*rdsu)-thve(k))/thve(k)*(zc(k)-zc(k-1))
  END DO
  IF(abe <= 0.)THEN
    igo=0
    RETURN
  END IF
!
  RETURN
END SUBROUTINE environc
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE LCL                        ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE lcl(t0,pp0,r0,tlcl,plcl,dzlcl)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Determine the LCL.
!
!-----------------------------------------------------------------------
!
!  AUTHOR:
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!  OUTPUT:
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  REAL :: t0,pp0,r0,tlcl,plcl,dzlcl

  REAL :: cpg,cpr,p00k
  DATA cpg/102.45/,cpr/3.4965/,p00k/26.91535/
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: nitt

  REAL :: dz,td,ttth0,ttd,rvs,r_rs,ti,pki,pppi,p0k,pi0i
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'phycst.inc'
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  plcl=pp0
  tlcl=t0
  p0k=pp0**rddcp

  pi0i=p0k/p00k*cp
  ttth0=t0*p00k/p0k
  ttd=td(pp0,r0)
  dz=cpg*(t0-ttd)

  IF(dz <= 0.)THEN
    dzlcl=0.
    RETURN
  END IF
  DO nitt=1,50
    pki=pi0i-g*dz/(ttth0*(1.+.61*r0))
    pppi=(pki/cp)**cpr*p0
    ti=ttth0*pki/cp
    rvs=r_rs(pppi,ti)
    IF (ABS(rvs-r0) < .00001) GO TO 110
    ttd=td(pppi,r0)
    dz=dz+cpg*(ti-ttd)
  END DO

  CALL arpsstop("arpstop called from lcl at 100",1)

  110   CONTINUE
  plcl=pppi
  tlcl=ti
  dzlcl=dz

  RETURN
END SUBROUTINE lcl
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE THE2T                      ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE the2t(the,p,th,t,r)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Convert potential temperature to temperature.
!
!-----------------------------------------------------------------------
!
!  AUTHOR:
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!  OUTPUT:
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  REAL :: the,p,th,t,r

!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: itter

  REAL :: TO,pi,tn,r_rs

!  REAL :: cp,alvl
!  DATA cp/1004./,alvl/2.68E6/
                                   ! new, suggested by G S Bhat
                                   ! and modified by Binbin Zhou
  pi=(p*1E-5)**.286

  TO=the/EXP((3.376/295.-0.00254)*1000*.012*(1+0.81*.012))*pi
  DO itter=1,50
    r=r_rs(p,TO)

    th=the/EXP((3.376/TO-0.00254)*1000*r*(1+0.81*r))
    tn=th*pi
    IF(ABS(TO-tn) < 0.005)GO TO 12
    TO=TO+(tn-TO)*.3
  END DO

  WRITE(6,1) the,p,TO,tn,th,r
  1     FORMAT(' STOP IN ROUTINE THE2T '/' THE,P,TO,TN,TH,R',6E15.6)

  CALL arpsstop("arpstop called from the2t",1)

  12    CONTINUE

  t=tn

  RETURN
END SUBROUTINE the2t
!
!##################################################################
!##################################################################
!######                                                      ######
!######                FUNCTION R_RS                         ######
!######                                                      ######
!##################################################################
!##################################################################
!

  FUNCTION r_rs(p,t)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Compute saturation mixing ratio.
!
!-----------------------------------------------------------------------
!
!  AUTHOR:
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!  OUTPUT:
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  REAL :: r_rs,p,t
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  REAL :: es
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  es=610.78*EXP(17.269*(t-273.16)/(t-35.86))
  r_rs=.622*es/(p-es)
  RETURN
  END FUNCTION r_rs
!
!##################################################################
!##################################################################
!######                                                      ######
!######                FUNCTION TD                           ######
!######                                                      ######
!##################################################################
!##################################################################
!

  FUNCTION td(p,rs)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!-----------------------------------------------------------------------
!
!  AUTHOR:
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!  OUTPUT:
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  REAL :: td,p,rs
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  REAL :: rr,es,esln
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  rr=rs+1E-8
  es=p*rr/(.622+rr)
  esln=LOG(es)
  td=(35.86*esln-4947.2325)/(esln-23.6837)
  RETURN
  END FUNCTION td
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE HTINT                      ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE htint(nzz1,vctra,eleva,nzz2,vctrb,elevb)
!
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!-----------------------------------------------------------------------
!
!  AUTHOR:
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!  OUTPUT:
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nzz1,nzz2

  REAL :: vctra(nzz1),vctrb(nzz2),eleva(nzz1),elevb(nzz2)
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: l,k
  REAL :: wt
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! stop SGI pfa compiler from parallelize the code

!*$*NOCONCURRENTIZE

  l=1
  DO k=1,nzz2
    30      CONTINUE

    IF (elevb(k) < eleva(1)) GO TO 35
    IF (elevb(k) >= eleva(l).AND.elevb(k) <= eleva(l+1)) GO TO 35
    IF (elevb(k) > eleva(nzz1)) GO TO 36

    l=l+1
    IF (l == nzz1) CALL arpsstop('arpstop called from htint',1)

    GO TO 30

    35      CONTINUE

    wt=(elevb(k)-eleva(l))/(eleva(l+1)-eleva(l))
    vctrb(k)=vctra(l)+(vctra(l+1)-vctra(l))*wt

    CYCLE

    36      CONTINUE
    wt=(elevb(k)-eleva(nzz1))/(eleva(nzz1-1)-eleva(nzz1))
    vctrb(k)=vctra(nzz1)+(vctra(nzz1-1)-vctra(nzz1))*wt

  END DO

! SGI pfa compiler directive:
!*$*CONCURRENTIZE

  RETURN
END SUBROUTINE htint
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE THETAE                     ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE thetae(p,t,rv,the)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Compute equivalent potential temperature.
!
!-----------------------------------------------------------------------
!
!  AUTHOR:
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!  OUTPUT:
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  REAL :: p,t,rv,the
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: itter

  REAL :: ttd,td,pit,tupo,tupn,tmn,dz,cpg

  REAL :: cp,g,r
  DATA cp/1004./,g/9.8/,r/287./    ! modified by Binbin Zhou
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  cpg=cp/g
  pit=p
  tupo=t
  ttd=td(p,rv)
  dz=cpg*(t-ttd)
  IF(dz <= 0.)GO TO 20
  DO itter=1,50
    tupn=t-g/cp*dz
    tmn=(tupn+t)*.5*(1.+.61*rv)
    pit=p*EXP(-g*dz/(r*tmn))
    IF(ABS(tupn-tupo) < 0.001)GO TO 20
    ttd=td(pit,rv)
    tupo=tupn
    dz=dz+cpg*(tupn-ttd)

  END DO

  CALL arpsstop("arpstop called from thetae",1)

  20    CONTINUE

  the=(3.376/tupo-0.00254)*1000*rv*(1+0.81*rv)
  the=tupo*(1E5/pit)**.286*EXP(the)

  RETURN
END SUBROUTINE thetae
!
!##################################################################
!##################################################################
!######                                                      ######
!######                FUNCTION S_SUM                        ######
!######                                                      ######
!##################################################################
!##################################################################
!

  FUNCTION s_sum(nn,vctr,inc)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!-----------------------------------------------------------------------
!
!  AUTHOR:
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!  OUTPUT:
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nn,inc

  REAL :: vctr(nn)
  REAL :: s_sum
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: k
  REAL :: sum
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  sum=0.
  DO k=1,nn,inc
    sum=sum+vctr(k)
  END DO
  s_sum=sum

  RETURN
  END FUNCTION s_sum
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE ES_CAL                     ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE es_cal(t,ess)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  This function computer vapor pressure over water (in Pascal).
!
!-----------------------------------------------------------------------
!
!  AUTHOR:
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!  OUTPUT:
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  REAL :: t,ess

!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  REAL :: a1,a2,b1,b2,tt,t1,t2,tice
  DATA a1/17.2694/a2/21.87456/
  DATA b1/35.86/b2/7.66/
  DATA tice/260.0/
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  tt=t
  IF (tt <= tice) GO TO 10
  t1=a1*(tt-273.16)
  t2=tt-b1
  GO TO 11
  10    t1=a2*(tt-273.16)
  t2=tt-b2
  11    ess=610.78*EXP(t1/t2)

  RETURN
END SUBROUTINE es_cal
