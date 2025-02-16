SUBROUTINE CALREF9s(nx_p,ny_p,nztn_p,qr,qs,qg,p,tb,refl)

! SUBPROGRAM:    CALREFL      CALCULATE Reflectivity
!
!   PRGMMR:  BENJAMIN, STAN ORG: NOAA/ESRL     DATE: 2007-02-10
!
! ABSTRACT:
!
!   Calculate reflectivity with both resolved
!     hydrometeors and latest 1h sub-grid-scale
!     (parameterized) precip.
!   This routine is a modified version of RUC calvis.f, which
!    uses many of the same input parameters.
!
!   This routine computes max 2-d reflectivity (in dBz)
!    from qr, qs, and qg.
!   qr--rain water mixing ratio  (kg/kg)
!   qs--snow mixing ratio        (kg/kg)
!   qg--graupel mixing ratio     (kg/kg)
!   p  -- pressure (Pa)
!   tb -- temperature (K)
!   refl -- reflectivity (dbZ) output
!   nx_p -- x dimension
!   ny_p -- y dimension
!   nztn_p -- z dimension
!
! HISTORY
! PROGRAM HISTORY LOG:
!     2007-02-12      S. Benjamin      Original version
!     2007-02-24      John M. Brown, S. Benjamin
!                       More detailed calculations
!			from John taken from NCAR/Thompson
!			microphysics.   John's original code
!                       was written in Apr 2005.
!
!     2007-04-04      Ming Hu
!                       convert to Fortran 90 format
!                       insert meteorology constant to make subroutine
!                               stand alone
!                       take all constant out of loop
!     2007-05-25      Fanyou Kong
!                       Rearrange arguement list
!
!------------------------------------------------------------------
!
   IMPLICIT   NONE
!      include 'IMPLICIT'
!      include 'METCON'
!      include 'PMICRPH'
!
!----------------------------------------------------
!   meteorological constant from METCON
!----------------------------------------------------
!
   REAL        PI_P
   REAL        R_P

   PARAMETER ( PI_P    = 3.14159265  )
   PARAMETER ( R_P     =  8.31451     )
!C--------------------------------------------------------------------
!C DRY AIR CONSTANTS
!C--------------------------------------------------------------------

   REAL        MD_P
   REAL        RD_P
   PARAMETER ( MD_P     =  0.0289645          )
   PARAMETER ( RD_P     =  R_P/MD_P           )
!
!----------------------------------------------------
!   microphysics variables from PMICRPH that computered by paramr.f
!----------------------------------------------------
   REAL  PI
   REAL  R1
   REAL  RON, RON2, SON, GON

   REAL RON_min,const1r,const2r,qr0,delqr0
   REAL DRAIN,DSNOW,DGRAUPEL,TOPR,TOPS,TOPG
   REAL constd,constgb,const_ng1,const_ng2
!
!
!----------------------------------------------------
!

   integer nx_p,ny_p,nztn_p,i,j,k

   REAL     &
       QR(NX_P,NY_P,nztn_p),  &
       QS(NX_P,NY_P,nztn_p),  &
       QG(NX_P,NY_P,nztn_p)   &
      , p (NX_P,NY_P,nztn_p)  &
      , tb(NX_P,NY_P,nztn_p)  &
      , refl(NX_P,NY_P,nztn_p)       &
      , ref_conv

   real refl_wat, refl_sno
   real reflx (nztn_p)

   real qrain, qsnow, qgraupel

   REAL rain,ronv,slor,snow,rhoqs,temp_c,sonv,slos
   REAL graupel,rhoqg,gonv,slog
   real alpha, rhod, bb
   real ze_s, ze_r, ze_g, ze_max, ze_nc, ze_conv, ze_sum

   real ze_smax, ze_rmax,ze_gmax

   real Arain, Asnow, Ahail
   real Brain, Bsnow, Bhail
   real Hail_min, Hail_max, Hail_cfa, Hail_cfb

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!------------------------------------------------------------------
!------------------------------------------------------------------
!      call paramr

   PI = PI_P
!  Min value for hydrometeor mixing ratios
   R1 = 1.E-15

! SLOPE INTERCEPT FOR RAIN, SNOW, AND GRAUPEL

!jmb--Roy R. suggests a larger value for the slope-intercept for rain.
!      This will slow down the fall speed.--16dec98
      RON=8.E6          ! Original M-P value.
!gt   RON2=1.E10        ! GREG T.  MWR Part1 vrbl intercept
      RON2=1.E9         ! GREG T.  changed 01 Dec 2004
!     SON=2.E6          ! Original M-P value.
      SON=2.E7
!jmb--According to Roy Rasmussens data (from a QJRMS paper he was reviewing)
!      the value of the M-P slope intercept can be as large as 3.E7 for
!      graupel.  The value GON = 4.E6 as an upper bound on the intercept value
!      appears too small.  Use same value as for snow.--17oct96
!     GON=4.E6          ! Original M-P value.
!gt   GON=5.e7          ! Roy R., summer 1998, 19 Jan 00, Oct 00
      GON=4.E6          ! Original M-P value.  GREG T.  changed 01 Dec 2004

!C DENSITY OF RAIN, SNOW, AND GRAUPEL

      DRAIN=1000.
      DSNOW=100.
      DGRAUPEL=400.

!C TOP OF SLOPE FOR RAIN, SNOW, AND GRAUPEL
!cjmb--By top of slope is meant numerator Marshall-Palmer slope parameter
!      [See (4) in Reisner et al 1998, QJRMS].

      TOPR=PI*DRAIN*RON
      TOPS=PI*DSNOW*SON
      TOPG=PI*DGRAUPEL*GON


! CONSTANTS FOR VARIABLE RON

!jmb  qr0 is center value of rain mixing ratio for transition from
!       M-P slope-intercept for drizzle formed by a collision-coalescence
!       process to M-P slope-intercept for traditional rain.--nov00
!jmb  delqr0 governs sharpness of transition: small delt_qr0 makes the
!      transition sharper:
!      if the rate of change of zero intercept wrt rain mixing ratio
!      were linear, with the slope at QR0 given by present tanh formula,
!      the transition would occur between qr0-delqr0 and qr0+delqr0.--nov00
!Cgt   RON_min = RON
      RON_min = 2.e7
!Cgt   qr0 = 0.0001  !  Roy R. 26oct00
      qr0 = 0.0002  !  GREG T.  01 Dec 2004
!Cgt   delqr0 = 0.25*qr0
      delqr0 = 0.5*qr0    !  GREG T.  01 Dec 2004
      const1r=(ron2-ron_min)*0.5
      const2r=(ron2+ron_min)*0.5

!CONSTANTS FOR VARIABLE GON
!     Based on Roy R s formulation, Jun 96

      constd=1./0.52
      constgb=12./13.
      const_ng1=(1.57**constd)*((pi*dgraupel)**constgb)
      const_ng2=-constgb
      Hail_max=(const_ng1/1.e4)**(1.0/constgb)
      Hail_min=(const_ng1/gon)**(1.0/constgb)
      Hail_cfa=constgb*0.75+1.75
      Hail_cfb=constgb*0.75-0.25

!
!-----------------------------------------------------------

      bb = 0.       !  bright band effect - yes or no (0)

      alpha = 0.224 ! = (1000kg/m^3/917kg/m^3)**2)*(0.176/0.930)
!                      1000kg/m^3 is density of liquid water
!                       917kg/m^3 is density of solid ice
!                      0.176 = dielectric factor of ice
!                      0.930 = dielectric factor of liquid water

      ze_smax = -1.E30
      ze_rmax = -1.E30
      ze_gmax = -1.E30
!
!  CONSTANT for RAIN reflectivity factor
!
  Arain=720/(pi**1.75 * Drain**1.75)
  Asnow=720*alpha*sqrt(sqrt(DSNOW))/(pi**1.75 * Drain*Drain)
  Ahail=720*alpha/( (1.57**(constd*0.75)) * pi**Hail_cfa *     &
           Drain*Drain * DGRAUPEL**Hail_cfb)
!      Ahail=720*alpha/(1.916164*pi**2.44 * Drain*Drain * DGRAUPEL**0.44)

  DO J=1,ny_p
    DO I=1,nx_p
      do k = 1,nztn_p

        rhod = p(i,j,k) / (rd_p * tb(i,j,k) )
        Brain=Arain*rhod**1.75
        Bsnow=Asnow*rhod**1.75
        Bhail=Ahail*rhod**Hail_cfa

!      Particle size distributions and reflectivity
!      ---------------------------------------------
!       Much of this code borrowed from EXMOISG loop 20 to get particle size
!       distributions

!jmb--    Note that SLOR, SLOS and SLOG are inverse slopes!!!! Also,
!          RONV,SONV,GONV, M-P zero intercept values, normalized by
!          max allowable values.

!         Have to set min values of hydrometeors (r1) large enough to avoid
!          underflow problems with log later on.

!   -- rain
        ze_r = 1.e-35
        if (qr(i,j,k).ge.1.e-6) then
          rain = max(r1,qr(i,j,k))
          ronv = (const1r*tanh((qr0 - rain)/delqr0) + const2r)
          ze_r=Brain * rain**1.75 / (ronv**0.75)
        endif

!   -- snow
        ze_s = 1.e-35
        if (qs(i,j,k).ge.1.e-6) then
          snow = max(r1,qs(i,j,k))
          temp_C = min(-0.001, Tb(i,j,k)-273.15)
          sonv = (min(2.0E8, 2.0E6*exp(-0.12*temp_C)))
          ze_s = Bsnow*snow**1.75/(sonv**0.75)

!         For bright band, increase reflectivity by factor of 5.28,
!          which is ratio of dielectric factors for water/ice (.930/.176)

          IF (tb(i,j,k) .gt. 273.15)  ze_s = ze_s*(1. + 4.28*bb)

        endif

!   -- graupel
        ze_g = 1.e-35
        if (qg(i,j,k).ge.1.e-6) then
            graupel = max(r1,qg(i,j,k))
! graupel = max(Hail_min/RHOd, min(graupel,Hail_max/RHOd))
          ze_g = Bhail * graupel**Hail_cfa

!         For bright band
          IF (tb(i,j,k) .gt. 273.15) ze_g = ze_g*(1. + 4.28*bb)

        endif

!   -- total grid scale
        ze_nc =  ze_r + ze_s + ze_g
        if( ze_nc > 1.e-30) then
          refl(i,j,k) = 10.*LOG10(ze_nc * 1.E18)
        else
          refl(i,j,k) =0.0
        endif

      end do
    end do
  end do
!            write(*,*)ze_r,ze_s,ze_g
!            write(*,*) refl(1,1,1),ze_nc

  RETURN

END SUBROUTINE CALREF9s

SUBROUTINE reflec_wrf(nx,ny,nz,Q,QQR,QQS,QQG,IICE,PMID,T,DBZ)

! Calculate reflectivity for microphysics schemes other than
! Ferrier and Thompson, based on WRF_POST algorithm
!
! Fanyou Kong
! 5/25/2007
!
! INPUT:
! nx,ny,nz dimension
! IICE    flag for warm-rain (0) or mixed-phase (1)
! Q       water vapor mixing ratio (kg/kg)
! QQR     rain water          (kg/kg)
! QQS     snow water          (kg/kg)
! QQG     grapeul/hail water  (kg/kg)
! PMID    air pressure (Pa)
! T       air temperature (K)
!
! OUTPUT:
! DBZ     reflectivity (dBZ)

  IMPLICIT   NONE
  INTEGER :: i,j,k,nx,ny,nz
  INTEGER :: IICE
  REAL :: Q  (nx,ny,nz)
  REAL :: QQR(nx,ny,nz)
  REAL :: QQS(nx,ny,nz)
  REAL :: QQG(nx,ny,nz)
  REAL :: PMID(nx,ny,nz)
  REAL :: T  (nx,ny,nz)
  REAL :: DBZ(nx,ny,nz)
  REAL :: DBZR(nx,ny,nz)
  REAL :: DBZI(nx,ny,nz)

  REAL, PARAMETER :: RD=287.04,TFRZ=273.16,D608=0.608,DBZmin=-20.
  REAL:: DENS

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  DBZI = 0.0

  DO k=1,nz-1
    DO J=1,ny-1
      DO I=1,nx-1
        IF(T(I,J,k) .LT. 1.0E-3) print*,'ZERO T'
        IF(T(I,J,k) .gt. 1.0E-3) DENS=PMID(I,J,k)/             &
           (RD*T(I,J,k)*(Q(I,J,k)*D608+1.0))      ! DENSITY

        IF(QQR(I,J,k).LT. 0.0)QQR(I,J,k)=0.0
        IF (IICE.EQ.0) THEN
          IF (T(I,J,k) .GE. TFRZ) THEN
             DBZ(I,J,k)=((QQR(I,J,k)*DENS)**1.75)*            &
                3.630803E-9 * 1.E18                  ! Z FOR RAIN
             DBZR(I,J,k)=DBZ(I,J,k)
          ELSE
             DBZ(I,J,k)=((QQR(I,J,k)*DENS)**1.75)*            &
                2.18500E-10 * 1.E18                  ! Z FOR SNOW
             DBZI(I,J,k)=DBZ(I,J,k)
          ENDIF
        ELSEIF (IICE.EQ.1) THEN
          IF(QQS(I,J,k).LT. 0.0)QQS(I,J,k)=0.0
          IF(QQG(I,J,k).LT. 0.0)QQG(I,J,k)=0.0
          DBZR(I,J,k)=((QQR(I,J,k)*DENS)**1.75)*              &
                3.630803E-9 * 1.E18                  ! Z FOR RAIN
          DBZI(I,J,k)= DBZI(I,J,k)+ &
             ((QQS(I,J,k)*DENS)**1.75)* &
                2.18500E-10 * 1.E18                  ! Z FOR SNOW
          DBZI(I,J,k)= DBZI(I,J,k)+ &
             ((QQG(I,J,k)*DENS)**1.75)* &
                1.033267E-9 * 1.E18                  ! Z FOR GRAUP
          DBZ(I,J,k)=DBZR(I,J,k)+DBZI(I,J,k)
        ELSE
          print *,'IICE=',IICE,' is invalid, STOP!'
          STOP
        ENDIF
        IF (DBZ(I,J,k).GT.0.) &
          DBZ(I,J,k)=10.0*LOG10(DBZ(I,J,k)) ! DBZ
!         IF (DBZR(I,J,k).GT.0.) &
!            DBZR(I,J,k)=10.0*LOG10(DBZR(I,J,k)) ! DBZ
!         IF (DBZI(I,J,k).GT.0.) &
!            DBZI(I,J,k)=10.0*LOG10(DBZI(I,J,k)) ! DBZ

          DBZ(I,J,k)=MAX(DBZmin, DBZ(I,J,k))
!         DBZR(I,J,k)=MAX(DBZmin, DBZR(I,J,k))
!         DBZI(I,J,k)=MAX(DBZmin, DBZI(I,J,k))
      ENDDO
    ENDDO
  ENDDO

  RETURN
END SUBROUTINE reflec_wrf

SUBROUTINE reflec_ferrier_wrf(nx,ny,nz,Q,QQC,QQR,QQS,PMID,T,DBZ,QQL,FS1D)

! Calculate reflectivity for Ferrier microphysics scheme, based on
! WRF_POST algorithm (CALMICT)
!
! Fanyou Kong
! 5/25/2007
!
! 3/28/2009 (F KONG)
! Modified to store and return large ice mixing ratio (FS1D) to approximate
! as graupel/hail category in Ferrier scheme
!
! 3/26/2011 (F KONG)
! - Modified to use WRF output F_RIMEF_PHY (as input in QQL) in determining
! graupel/hail from total ice
! - set new threshold F_RIMEF_PHY >=5 for graupel/hail (Ferrier personal
!   communication)
!
! INPUT:
! nx,ny,nz dimension
! Q       water vapor mixing ratio (kg/kg)
! QQC     cloud water         (kg/kg)
! QQR     rain water          (kg/kg)
! QQS     total ice water     (kg/kg)
! FS1D    F_RIMEF_PHY
! PMID    air pressure (Pa)
! T       air temperature (K)
!
! OUTPUT:
! DBZ     reflectivity (dBZ)
! QQL     graupel/hail        (kg/kg)

  IMPLICIT   NONE

  include "grid.inc"
  INCLUDE 'mp.inc'

  INTEGER :: i,j,k,nx,ny,nz

  REAL :: Q  (nx,ny,nz)
  REAL :: QQC(nx,ny,nz)
  REAL :: QQR(nx,ny,nz)
  REAL :: QQS(nx,ny,nz)
  REAL :: PMID(nx,ny,nz)
  REAL :: T  (nx,ny,nz)
  REAL :: DBZ(nx,ny,nz)
  REAL :: DBZR(nx,ny,nz)
  REAL :: DBZI(nx,ny,nz)

  REAL :: QQL(nx,ny,nz)
  REAL :: FS1D(nx,ny,nz)

  REAL :: MASSR(50:450),MASSI(50:1000)

  REAL, PARAMETER :: RD=287.04,TFRZ=273.16,D608=0.608,DBZmin=-20.
  REAL, PARAMETER :: EPSQ=1.E-12,RHOL=1000.0,ONEPS=1.-18.015/28.964
  REAL, PARAMETER :: DMRmin=.05E-3, DMRmax=.45E-3, N0r0=8.E6,N0rmin=1.e4
  REAL, PARAMETER :: XMRmin=1.E6*DMRmin, XMRmax=1.E6*DMRmax
  INTEGER, PARAMETER :: MDRmin=XMRmin, MDRmax=XMRmax
  REAL, PARAMETER :: DMImin=.05E-3, DMImax=1E-3
  REAL, PARAMETER :: XMImin=1.E6*DMImin, XMImax=1.E6*DMImax
  INTEGER, PARAMETER :: MDImin=50, MDImax=1000
  REAL, PARAMETER :: Cice=1.634e13
  REAL, PARAMETER :: NLImin=1.E3, NLImax=5.E3

  REAL:: f_qvsat
  REAL:: QSAT,RHO,RRHO,TC,PI,C_N0r0,CN0r0,CN0r_DMRmin,CN0r_DMRmax
  REAL:: RQR,N0r,DRmm,Zrain,Zice,Ztot
  REAL:: QICE,QSIgrd,WVQW,RQR_DRmin,RQR_DRmax,RHgrd
  REAL:: Zmin,AX,XLI,DUM,RimeF,QLICE,NLICE,DLI
  REAL:: FLARGE,FSMALL,XSIMASS,XLIMASS,FLIMASS
  INTEGER :: INDEXR,INDEXS

  INTEGER :: istatus

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  DBZ = DBZmin
  DBZR= DBZmin
  DBZI= DBZmin

  Zmin=10.**(0.1*DBZmin)
  PI=ACOS(-1.)
  C_N0r0=PI*RHOL*N0r0
  CN0r0=1.E6/C_N0r0**.25
  CN0r_DMRmin=1./(PI*RHOL*DMRmin**4)
  CN0r_DMRmax=1./(PI*RHOL*DMRmax**4)
  !RQR_DRmin=N0r0*MASSR(MDRmin)    ! Rain content for mean drop diameter of .05 mm
  !RQR_DRmax=N0r0*MASSR(MDRmax)    ! Rain content for mean drop diameter of .45 mm
  RQR_DRmin=N0r0*1.9274059E-14    ! Rain content for mean drop diameter of .05 mm
  RQR_DRmax=N0r0*1.2882169E-10    ! Rain content for mean drop diameter of .45 mm

  !AX=111.*(DYVAL/1000.**2+DXVAL/1000.**2)**.5
  AX=111.*(dy/1000.**2+dx/1000.**2)**.5
  AX=MIN(100., MAX(5., AX) )
  RHgrd=0.90+.08*((100.-AX)/95.)**.5

!  OPEN (UNIT=1,FILE="/usr/users/5/fkong/input0/eta_micro_lookup.dat", &
!  OPEN (UNIT=1,FILE="eta_micro_lookup.dat", FORM="UNFORMATTED")
!!  OPEN (UNIT=1,FILE="ETAMPNEW_DATA", FORM="UNFORMATTED")
!!
!!  DO I=1,3
!!    READ(1)
!!  ENDDO
!!  READ(1) MASSR
!!  DO I=1,5
!!    READ(1)
!!  ENDDO
!!  READ(1) MASSI
!!  CLOSE(1)

  CALL init_massr_massi(massr,50,450,massi,50,1000,istatus)

  IF (myproc == 0) print *,'MASSI(MDImin),MASSI(MDImax):',MASSI(50),MASSI(1000)

  !QQL = QQC + QQR + QQS
  QQL(:,:,:) = 0.0

  DO k=1,nz
    DO J=1,ny
      DO I=1,nx

        Zrain=0.            !--- Radar reflectivity from rain
        Zice=0.             !--- Radar reflectivity from ice
        TC=T(I,J,k)-TFRZ

        IF(QQC(I,J,k)+QQR(I,J,k)+QQS(I,J,k) .LE. EPSQ) GO TO 10


!       ESAT=1000.*FPVS(T1D(I,J,k))
!       QSAT=EPS*ESAT/(P1D(I,J,k)-ESAT)
        QSAT=f_qvsat(PMID(I,J,k), T(I,J,k))         ! using arps f_qvsat
        RHO=PMID(I,J,k)/(RD*T(I,J,k)*(Q(I,J,k)*D608+1.0))      ! DENSITY
        RRHO=1./RHO
!
!--- Based on code from GSMCOLUMN in model to determine reflectivity from rain
!
        IF (QQR(I,J,k) .GT. EPSQ) THEN
          RQR=RHO*QQR(I,J,k)
          IF (RQR .LE. RQR_DRmin) THEN
            N0r=MAX(N0rmin, CN0r_DMRmin*RQR)
            INDEXR=MDRmin
          ELSE IF (RQR .GE. RQR_DRmax) THEN
            N0r=CN0r_DMRmax*RQR
            INDEXR=MDRmax
          ELSE
            N0r=N0r0
            INDEXR=MAX( XMRmin, MIN(CN0r0*RQR**.25, XMRmax) )
          ENDIF
!
!--- INEXR is the mean drop size in microns; convert to mm
!
          DRmm=1.e-3*REAL(INDEXR)
          Zrain=0.72*N0r*DRmm*DRmm*DRmm*DRmm*DRmm*DRmm*DRmm
        ENDIF        !--- End IF (QQR(I,J,k) .GT. EPSQ) block
!
!--- Based on code from GSMCOLUMN in model to determine partition of
!    total ice into cloud ice & snow (precipitation ice)
        IF (QQS(I,J,k) .GT. EPSQ) THEN
          QICE=QQS(I,J,k)
          RHO=PMID(I,J,k)/ &
             (RD*T(I,J,k)*(Q(I,J,k)*ONEPS+1.0))      ! DENSITY
          RRHO=1./RHO
          QSIgrd=RHgrd*QSAT
          WVQW=QQC(I,J,k)+Q(I,J,k)/(1.-Q(I,J,k))
!
! * FLARGE  - ratio of number of large ice to total (large & small) ice
! * FSMALL  - ratio of number of small ice crystals to large ice particles
!  ->  Small ice particles are assumed to have a mean diameter of 50 microns.
!  * XSIMASS - used for calculating small ice mixing ratio
!  * XLIMASS - used for calculating large ice mixing ratio
!  * INDEXS  - mean size of snow to the nearest micron (units of microns)
!  * RimeF   - Rime Factor, which is the mass ratio of total (unrimed &
!              rimed) ice mass to the unrimed ice mass (>=1)
!  * FLIMASS - mass fraction of large ice
!  * QTICE   - time-averaged mixing ratio of total ice
!  * QLICE   - time-averaged mixing ratio of large ice
!  * NLICE   - time-averaged number concentration of large ice
!
          IF (TC.GE.0. .OR. WVQW.LT.QSIgrd) THEN
            FLARGE=1.
          ELSE
            FLARGE=.2
            IF (TC.GE.-8. .AND. TC.LE.-3.) FLARGE=.5*FLARGE
          ENDIF
          FSMALL=(1.-FLARGE)/FLARGE
          XSIMASS=RRHO*MASSI(MDImin)*FSMALL
          DUM=XMImax*EXP(.0536*TC)
          INDEXS=MIN(MDImax, MAX(MDImin, INT(DUM) ) )
          RimeF=AMAX1(1., FS1D(I,J,k) )
!          RimeF=2.0
          XLIMASS=RRHO*RimeF*MASSI(INDEXS)
          FLIMASS=XLIMASS/(XLIMASS+XSIMASS)
          QLICE=FLIMASS*QICE
          NLICE=QLICE/XLIMASS
          IF (NLICE.LT.NLImin .OR. NLICE.GT.NLImax) THEN
!
!--- Force NLICE to be between NLImin and NLImax
!
             DUM=MAX(NLImin, MIN(NLImax, NLICE) )
             XLI=RHO*(QICE/DUM-XSIMASS)/RimeF
             IF (XLI .LE. MASSI(MDImin) ) THEN
               INDEXS=MDImin
             ELSE IF (XLI .LE. MASSI(450) ) THEN
               DLI=9.5885E5*XLI**.42066         ! DLI in microns
               INDEXS=MIN(MDImax, MAX(MDImin, INT(DLI) ) )
             ELSE IF (XLI .LE. MASSI(MDImax) ) THEN
               DLI=3.9751E6*XLI**.49870         ! DLI in microns
               INDEXS=MIN(MDImax, MAX(MDImin, INT(DLI) ) )
             ELSE
               INDEXS=MDImax
!
!--- 8/22/01: Increase density of large ice if maximum limits
!    are reached for number concentration (NLImax) and mean size
!    (MDImax).  Done to increase fall out of ice.
!
               IF (DUM .GE. NLImax)                                     &
                 RimeF=RHO*(QICE/NLImax-XSIMASS)/MASSI(INDEXS)
             END IF             ! End IF (XLI .LE. MASSI(MDImin) )
             XLIMASS=RRHO*RimeF*MASSI(INDEXS)
             FLIMASS=XLIMASS/(XLIMASS+XSIMASS)
             QLICE=FLIMASS*QICE
!            NLICE=QLICE/XLIMASS
          ENDIF               ! End IF (NLICE.LT.NLImin ...
!            QS1(I,J,k)=AMIN1(QI1(I,J,k), QLICE)
!            QI1(I,J,k)=AMAX1(0., QI1(I,J,k)-QS1(I,J,k))
   !
   !--- Equation (C.8) in Ferrier (1994, JAS, p. 272), which when
   !    converted from cgs units to mks units results in the same
   !    value for Cice, which is equal to the {} term below:
   !
   !    Zi={.224*720*(10**18)/[(PI*RHOL)**2]}*(RHO*QLICE)**2/NLICE,
   !    where RHOL=1000 kg/m**3 is the density of liquid water
   !
   !--- Valid only for exponential ice distributions
   !
!         Zice=Cice*RHO*RHO*QLICE*QLICE/NLICE
          Zice=Cice*RHO*RHO*QLICE*XLIMASS

          IF(RimeF .GE. 5) QQL(I,J,k) = QLICE

        ENDIF                 ! End IF (QI1(I,J,k) .GT. 0.) THEN
!
!---  Calculate total (grid-scale) radar reflectivity
        10 Ztot=Zrain+Zice
        IF (Ztot .GT. Zmin)  DBZ(I,J,k)= 10.*ALOG10(Ztot)
!          IF (Zrain .GT. Zmin) DBZR(I,J,k)=10.*ALOG10(Zrain)
!          IF (Zice .GT. Zmin)  DBZI(I,J,k)=10.*ALOG10(Zice)

      ENDDO
    ENDDO
  ENDDO
  RETURN
END SUBROUTINE reflec_ferrier_wrf

SUBROUTINE init_massr_massi(massr,rbgn,rend,massi,ibgn,iend,istatus)
!#######################################################################
!
! Hard-coded data read in from file ETAMPNEW_DATA in WRFV3.3.1 with code:
!
!!  OPEN (UNIT=1,FILE="ETAMPNEW_DATA", FORM="UNFORMATTED")
!!
!!  DO I=1,3
!!    READ(1)
!!  ENDDO
!!  READ(1) MASSR
!!  DO I=1,5
!!    READ(1)
!!  ENDDO
!!  READ(1) MASSI
!!  CLOSE(1)
!
!#######################################################################

  IMPLICIT NONE

  INTEGER, INTENT(IN)  :: rbgn, rend
  REAL,    INTENT(OUT) :: massr(rbgn:rend)
  INTEGER, INTENT(IN)  :: ibgn, iend
  REAL,    INTENT(OUT) :: massi(ibgn:iend)
  INTEGER, INTENT(OUT) :: istatus

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  !           401          951

  istatus = 0

! ... massr
   massr(  50 ) =    1.9274059E-14
   massr(  51 ) =    2.0887141E-14
   massr(  52 ) =    2.2598456E-14
   massr(  53 ) =    2.4411794E-14
   massr(  54 ) =    2.6331197E-14
   massr(  55 ) =    2.8360623E-14
   massr(  56 ) =    3.0504204E-14
   massr(  57 ) =    3.2766094E-14
   massr(  58 ) =    3.5150665E-14
   massr(  59 ) =    3.7662144E-14
   massr(  60 ) =    4.0304945E-14
   massr(  61 ) =    4.3083524E-14
   massr(  62 ) =    4.6002552E-14
   massr(  63 ) =    4.9066664E-14
   massr(  64 ) =    5.2280259E-14
   massr(  65 ) =    5.5648537E-14
   massr(  66 ) =    5.9176046E-14
   massr(  67 ) =    6.2867992E-14
   massr(  68 ) =    6.6729228E-14
   massr(  69 ) =    7.0764798E-14
   massr(  70 ) =    7.4980000E-14
   massr(  71 ) =    7.9380086E-14
   massr(  72 ) =    8.3970184E-14
   massr(  73 ) =    8.8755826E-14
   massr(  74 ) =    9.3742485E-14
   massr(  75 ) =    9.8935339E-14
   massr(  76 ) =    1.0434071E-13
   massr(  77 ) =    1.0996378E-13
   massr(  78 ) =    1.1581040E-13
   massr(  79 ) =    1.2188632E-13
   massr(  80 ) =    1.2819770E-13
   massr(  81 ) =    1.3475030E-13
   massr(  82 ) =    1.4155048E-13
   massr(  83 ) =    1.4860425E-13
   massr(  84 ) =    1.5591758E-13
   massr(  85 ) =    1.6349736E-13
   massr(  86 ) =    1.7134905E-13
   massr(  87 ) =    1.7947998E-13
   massr(  88 ) =    1.8789610E-13
   massr(  89 ) =    1.9660433E-13
   massr(  90 ) =    2.0561118E-13
   massr(  91 ) =    2.1492357E-13
   massr(  92 ) =    2.2454796E-13
   massr(  93 ) =    2.3449141E-13
   massr(  94 ) =    2.4476122E-13
   massr(  95 ) =    2.5536414E-13
   massr(  96 ) =    2.6630691E-13
   massr(  97 ) =    2.7759739E-13
   massr(  98 ) =    2.8924264E-13
   massr(  99 ) =    3.0125064E-13
   massr( 100 ) =    3.1362719E-13
   massr( 101 ) =    3.2638204E-13
   massr( 102 ) =    3.3952027E-13
   massr( 103 ) =    3.5305117E-13
   massr( 104 ) =    3.6698219E-13
   massr( 105 ) =    3.8132095E-13
   massr( 106 ) =    3.9607578E-13
   massr( 107 ) =    4.1125317E-13
   massr( 108 ) =    4.2686248E-13
   massr( 109 ) =    4.4291204E-13
   massr( 110 ) =    4.5941002E-13
   massr( 111 ) =    4.7636222E-13
   massr( 112 ) =    4.9378069E-13
   massr( 113 ) =    5.1167149E-13
   massr( 114 ) =    5.3004367E-13
   massr( 115 ) =    5.4890624E-13
   massr( 116 ) =    5.6826787E-13
   massr( 117 ) =    5.8813664E-13
   massr( 118 ) =    6.0851964E-13
   massr( 119 ) =    6.2942923E-13
   massr( 120 ) =    6.5087145E-13
   massr( 121 ) =    6.7285825E-13
   massr( 122 ) =    6.9539584E-13
   massr( 123 ) =    7.1849449E-13
   massr( 124 ) =    7.4216458E-13
   massr( 125 ) =    7.6641417E-13
   massr( 126 ) =    7.9125308E-13
   massr( 127 ) =    8.1668944E-13
   massr( 128 ) =    8.4273371E-13
   massr( 129 ) =    8.6939706E-13
   massr( 130 ) =    8.9668670E-13
   massr( 131 ) =    9.2461314E-13
   massr( 132 ) =    9.5318458E-13
   massr( 133 ) =    9.8241565E-13
   massr( 134 ) =    1.0123140E-12
   massr( 135 ) =    1.0428884E-12
   massr( 136 ) =    1.0741509E-12
   massr( 137 ) =    1.1061098E-12
   massr( 138 ) =    1.1387767E-12
   massr( 139 ) =    1.1721621E-12
   massr( 140 ) =    1.2062764E-12
   massr( 141 ) =    1.2411270E-12
   massr( 142 ) =    1.2767308E-12
   massr( 143 ) =    1.3130921E-12
   massr( 144 ) =    1.3502225E-12
   massr( 145 ) =    1.3881367E-12
   massr( 146 ) =    1.4268456E-12
   massr( 147 ) =    1.4663549E-12
   massr( 148 ) =    1.5066804E-12
   massr( 149 ) =    1.5478326E-12
   massr( 150 ) =    1.5898213E-12
   massr( 151 ) =    1.6326557E-12
   massr( 152 ) =    1.6763520E-12
   massr( 153 ) =    1.7209174E-12
   massr( 154 ) =    1.7663683E-12
   massr( 155 ) =    1.8127135E-12
   massr( 156 ) =    1.8599625E-12
   massr( 157 ) =    1.9081243E-12
   massr( 158 ) =    1.9572163E-12
   massr( 159 ) =    2.0072557E-12
   massr( 160 ) =    2.0582457E-12
   massr( 161 ) =    2.1102007E-12
   massr( 162 ) =    2.1631353E-12
   massr( 163 ) =    2.2170557E-12
   massr( 164 ) =    2.2719795E-12
   massr( 165 ) =    2.3279200E-12
   massr( 166 ) =    2.3848801E-12
   massr( 167 ) =    2.4428849E-12
   massr( 168 ) =    2.5019327E-12
   massr( 169 ) =    2.5620545E-12
   massr( 170 ) =    2.6232506E-12
   massr( 171 ) =    2.6855373E-12
   massr( 172 ) =    2.7489250E-12
   massr( 173 ) =    2.8134292E-12
   massr( 174 ) =    2.8790542E-12
   massr( 175 ) =    2.9458283E-12
   massr( 176 ) =    3.0137550E-12
   massr( 177 ) =    3.0828478E-12
   massr( 178 ) =    3.1531225E-12
   massr( 179 ) =    3.2245929E-12
   massr( 180 ) =    3.2972704E-12
   massr( 181 ) =    3.3711703E-12
   massr( 182 ) =    3.4463040E-12
   massr( 183 ) =    3.5226895E-12
   massr( 184 ) =    3.6003383E-12
   massr( 185 ) =    3.6792540E-12
   massr( 186 ) =    3.7594524E-12
   massr( 187 ) =    3.8409714E-12
   massr( 188 ) =    3.9237979E-12
   massr( 189 ) =    4.0079680E-12
   massr( 190 ) =    4.0934786E-12
   massr( 191 ) =    4.1803574E-12
   massr( 192 ) =    4.2686080E-12
   massr( 193 ) =    4.3582447E-12
   massr( 194 ) =    4.4492903E-12
   massr( 195 ) =    4.5417524E-12
   massr( 196 ) =    4.6356417E-12
   massr( 197 ) =    4.7309916E-12
   massr( 198 ) =    4.8277970E-12
   massr( 199 ) =    4.9260873E-12
   massr( 200 ) =    5.0258612E-12
   massr( 201 ) =    5.1271483E-12
   massr( 202 ) =    5.2299593E-12
   massr( 203 ) =    5.3343168E-12
   massr( 204 ) =    5.4402229E-12
   massr( 205 ) =    5.5476851E-12
   massr( 206 ) =    5.6567472E-12
   massr( 207 ) =    5.7673900E-12
   massr( 208 ) =    5.8796661E-12
   massr( 209 ) =    5.9935707E-12
   massr( 210 ) =    6.1091176E-12
   massr( 211 ) =    6.2263293E-12
   massr( 212 ) =    6.3452082E-12
   massr( 213 ) =    6.4657919E-12
   massr( 214 ) =    6.5880899E-12
   massr( 215 ) =    6.7121109E-12
   massr( 216 ) =    6.8378753E-12
   massr( 217 ) =    6.9653918E-12
   massr( 218 ) =    7.0946864E-12
   massr( 219 ) =    7.2257816E-12
   massr( 220 ) =    7.3586779E-12
   massr( 221 ) =    7.4933514E-12
   massr( 222 ) =    7.6299262E-12
   massr( 223 ) =    7.7683302E-12
   massr( 224 ) =    7.9086338E-12
   massr( 225 ) =    8.0508187E-12
   massr( 226 ) =    8.1949005E-12
   massr( 227 ) =    8.3409243E-12
   massr( 228 ) =    8.4889084E-12
   massr( 229 ) =    8.6388197E-12
   massr( 230 ) =    8.7907138E-12
   massr( 231 ) =    8.9446115E-12
   massr( 232 ) =    9.1005172E-12
   massr( 233 ) =    9.2584464E-12
   massr( 234 ) =    9.4184539E-12
   massr( 235 ) =    9.5804926E-12
   massr( 236 ) =    9.7446218E-12
   massr( 237 ) =    9.9108395E-12
   massr( 238 ) =    1.0079187E-11
   massr( 239 ) =    1.0249670E-11
   massr( 240 ) =    1.0422278E-11
   massr( 241 ) =    1.0597109E-11
   massr( 242 ) =    1.0774064E-11
   massr( 243 ) =    1.0953300E-11
   massr( 244 ) =    1.1134731E-11
   massr( 245 ) =    1.1318405E-11
   massr( 246 ) =    1.1504339E-11
   massr( 247 ) =    1.1692545E-11
   massr( 248 ) =    1.1883086E-11
   massr( 249 ) =    1.2075890E-11
   massr( 250 ) =    1.2271101E-11
   massr( 251 ) =    1.2468622E-11
   massr( 252 ) =    1.2668517E-11
   massr( 253 ) =    1.2870819E-11
   massr( 254 ) =    1.3075559E-11
   massr( 255 ) =    1.3282673E-11
   massr( 256 ) =    1.3492272E-11
   massr( 257 ) =    1.3704334E-11
   massr( 258 ) =    1.3918898E-11
   massr( 259 ) =    1.4135964E-11
   massr( 260 ) =    1.4355555E-11
   massr( 261 ) =    1.4577702E-11
   massr( 262 ) =    1.4802397E-11
   massr( 263 ) =    1.5029596E-11
   massr( 264 ) =    1.5259544E-11
   massr( 265 ) =    1.5492024E-11
   massr( 266 ) =    1.5727232E-11
   massr( 267 ) =    1.5965097E-11
   massr( 268 ) =    1.6205676E-11
   massr( 269 ) =    1.6448830E-11
   massr( 270 ) =    1.6694823E-11
   massr( 271 ) =    1.6943547E-11
   massr( 272 ) =    1.7195068E-11
   massr( 273 ) =    1.7449306E-11
   massr( 274 ) =    1.7706418E-11
   massr( 275 ) =    1.7966319E-11
   massr( 276 ) =    1.8229081E-11
   massr( 277 ) =    1.8494737E-11
   massr( 278 ) =    1.8763255E-11
   massr( 279 ) =    1.9034678E-11
   massr( 280 ) =    1.9309115E-11
   massr( 281 ) =    1.9586452E-11
   massr( 282 ) =    1.9866797E-11
   massr( 283 ) =    2.0150099E-11
   massr( 284 ) =    2.0436392E-11
   massr( 285 ) =    2.0725772E-11
   massr( 286 ) =    2.1018198E-11
   massr( 287 ) =    2.1313729E-11
   massr( 288 ) =    2.1612363E-11
   massr( 289 ) =    2.1914085E-11
   massr( 290 ) =    2.2218982E-11
   massr( 291 ) =    2.2527038E-11
   massr( 292 ) =    2.2838316E-11
   massr( 293 ) =    2.3152815E-11
   massr( 294 ) =    2.3470507E-11
   massr( 295 ) =    2.3791467E-11
   massr( 296 ) =    2.4115753E-11
   massr( 297 ) =    2.4443253E-11
   massr( 298 ) =    2.4774151E-11
   massr( 299 ) =    2.5108424E-11
   massr( 300 ) =    2.5445993E-11
   massr( 301 ) =    2.5786984E-11
   massr( 302 ) =    2.6131352E-11
   massr( 303 ) =    2.6479267E-11
   massr( 304 ) =    2.6830550E-11
   massr( 305 ) =    2.7185348E-11
   massr( 306 ) =    2.7543568E-11
   massr( 307 ) =    2.7905459E-11
   massr( 308 ) =    2.8270809E-11
   massr( 309 ) =    2.8639801E-11
   massr( 310 ) =    2.9012350E-11
   massr( 311 ) =    2.9388582E-11
   massr( 312 ) =    2.9768157E-11
   massr( 313 ) =    3.0151621E-11
   massr( 314 ) =    3.0538818E-11
   massr( 315 ) =    3.0929672E-11
   massr( 316 ) =    3.1324294E-11
   massr( 317 ) =    3.1722833E-11
   massr( 318 ) =    3.2124973E-11
   massr( 319 ) =    3.2530974E-11
   massr( 320 ) =    3.2940827E-11
   massr( 321 ) =    3.3354517E-11
   massr( 322 ) =    3.3772159E-11
   massr( 323 ) =    3.4193669E-11
   massr( 324 ) =    3.4618936E-11
   massr( 325 ) =    3.5048492E-11
   massr( 326 ) =    3.5481791E-11
   massr( 327 ) =    3.5919209E-11
   massr( 328 ) =    3.6360644E-11
   massr( 329 ) =    3.6806027E-11
   massr( 330 ) =    3.7255615E-11
   massr( 331 ) =    3.7709239E-11
   massr( 332 ) =    3.8167077E-11
   massr( 333 ) =    3.8628920E-11
   massr( 334 ) =    3.9095123E-11
   massr( 335 ) =    3.9565455E-11
   massr( 336 ) =    4.0040062E-11
   massr( 337 ) =    4.0518863E-11
   massr( 338 ) =    4.1001941E-11
   massr( 339 ) =    4.1489250E-11
   massr( 340 ) =    4.1981089E-11
   massr( 341 ) =    4.2477067E-11
   massr( 342 ) =    4.2977601E-11
   massr( 343 ) =    4.3482534E-11
   massr( 344 ) =    4.3991779E-11
   massr( 345 ) =    4.4505552E-11
   massr( 346 ) =    4.5023780E-11
   massr( 347 ) =    4.5546771E-11
   massr( 348 ) =    4.6073922E-11
   massr( 349 ) =    4.6605865E-11
   massr( 350 ) =    4.7142266E-11
   massr( 351 ) =    4.7683465E-11
   massr( 352 ) =    4.8229223E-11
   massr( 353 ) =    4.8779602E-11
   massr( 354 ) =    4.9334713E-11
   massr( 355 ) =    4.9894509E-11
   massr( 356 ) =    5.0459179E-11
   massr( 357 ) =    5.1028466E-11
   massr( 358 ) =    5.1602594E-11
   massr( 359 ) =    5.2181631E-11
   massr( 360 ) =    5.2765532E-11
   massr( 361 ) =    5.3354245E-11
   massr( 362 ) =    5.3947923E-11
   massr( 363 ) =    5.4546527E-11
   massr( 364 ) =    5.5149933E-11
   massr( 365 ) =    5.5758585E-11
   massr( 366 ) =    5.6372167E-11
   massr( 367 ) =    5.6990825E-11
   massr( 368 ) =    5.7614472E-11
   massr( 369 ) =    5.8243285E-11
   massr( 370 ) =    5.8877271E-11
   massr( 371 ) =    5.9516336E-11
   massr( 372 ) =    6.0160044E-11
   massr( 373 ) =    6.0809642E-11
   massr( 374 ) =    6.1464417E-11
   massr( 375 ) =    6.2124320E-11
   massr( 376 ) =    6.2789711E-11
   massr( 377 ) =    6.3460459E-11
   massr( 378 ) =    6.4136405E-11
   massr( 379 ) =    6.4817859E-11
   massr( 380 ) =    6.5504650E-11
   massr( 381 ) =    6.6196937E-11
   massr( 382 ) =    6.6894650E-11
   massr( 383 ) =    6.7597976E-11
   massr( 384 ) =    6.8306590E-11
   massr( 385 ) =    6.9021067E-11
   massr( 386 ) =    6.9740998E-11
   massr( 387 ) =    7.0466549E-11
   massr( 388 ) =    7.1197707E-11
   massr( 389 ) =    7.1934535E-11
   massr( 390 ) =    7.2677107E-11
   massr( 391 ) =    7.3425356E-11
   massr( 392 ) =    7.4179343E-11
   massr( 393 ) =    7.4939402E-11
   massr( 394 ) =    7.5704831E-11
   massr( 395 ) =    7.6476568E-11
   massr( 396 ) =    7.7253794E-11
   massr( 397 ) =    7.8037209E-11
   massr( 398 ) =    7.8826425E-11
   massr( 399 ) =    7.9621601E-11
   massr( 400 ) =    8.0422898E-11
   massr( 401 ) =    8.1230099E-11
   massr( 402 ) =    8.2043587E-11
   massr( 403 ) =    8.2863098E-11
   massr( 404 ) =    8.3688570E-11
   massr( 405 ) =    8.4520085E-11
   massr( 406 ) =    8.5357998E-11
   massr( 407 ) =    8.6202032E-11
   massr( 408 ) =    8.7052358E-11
   massr( 409 ) =    8.7908951E-11
   massr( 410 ) =    8.8771809E-11
   massr( 411 ) =    8.9641281E-11
   massr( 412 ) =    9.0516955E-11
   massr( 413 ) =    9.1398743E-11
   massr( 414 ) =    9.2287296E-11
   massr( 415 ) =    9.3182295E-11
   massr( 416 ) =    9.4083616E-11
   massr( 417 ) =    9.4991626E-11
   massr( 418 ) =    9.5906130E-11
   massr( 419 ) =    9.6827220E-11
   massr( 420 ) =    9.7754908E-11
   massr( 421 ) =    9.8689064E-11
   massr( 422 ) =    9.9630394E-11
   massr( 423 ) =    1.0057801E-10
   massr( 424 ) =    1.0153237E-10
   massr( 425 ) =    1.0249372E-10
   massr( 426 ) =    1.0346197E-10
   massr( 427 ) =    1.0443691E-10
   massr( 428 ) =    1.0541857E-10
   massr( 429 ) =    1.0640708E-10
   massr( 430 ) =    1.0740271E-10
   massr( 431 ) =    1.0840540E-10
   massr( 432 ) =    1.0941510E-10
   massr( 433 ) =    1.1043177E-10
   massr( 434 ) =    1.1145548E-10
   massr( 435 ) =    1.1248615E-10
   massr( 436 ) =    1.1352422E-10
   massr( 437 ) =    1.1456933E-10
   massr( 438 ) =    1.1562159E-10
   massr( 439 ) =    1.1668112E-10
   massr( 440 ) =    1.1774788E-10
   massr( 441 ) =    1.1882213E-10
   massr( 442 ) =    1.1990231E-10
   massr( 443 ) =    1.2099076E-10
   massr( 444 ) =    1.2208716E-10
   massr( 445 ) =    1.2319085E-10
   massr( 446 ) =    1.2430205E-10
   massr( 447 ) =    1.2542048E-10
   massr( 448 ) =    1.2654673E-10
   massr( 449 ) =    1.2768003E-10
   massr( 450 ) =    1.2882169E-10

! ... massi
   massi(   50 ) =    7.5019720E-11
   massi(   51 ) =    7.8376493E-11
   massi(   52 ) =    8.1811120E-11
   massi(   53 ) =    8.5323797E-11
   massi(   54 ) =    8.8915125E-11
   massi(   55 ) =    9.2585772E-11
   massi(   56 ) =    9.6336057E-11
   massi(   57 ) =    1.0016658E-10
   massi(   58 ) =    1.0407813E-10
   massi(   59 ) =    1.0807087E-10
   massi(   60 ) =    1.1214574E-10
   massi(   61 ) =    1.1630324E-10
   massi(   62 ) =    1.2054374E-10
   massi(   63 ) =    1.2486814E-10
   massi(   64 ) =    1.2927692E-10
   massi(   65 ) =    1.3377034E-10
   massi(   66 ) =    1.3834967E-10
   massi(   67 ) =    1.4301516E-10
   massi(   68 ) =    1.4776722E-10
   massi(   69 ) =    1.5260679E-10
   massi(   70 ) =    1.5753424E-10
   massi(   71 ) =    1.6255045E-10
   massi(   72 ) =    1.6765551E-10
   massi(   73 ) =    1.7285053E-10
   massi(   74 ) =    1.7813617E-10
   massi(   75 ) =    1.8351269E-10
   massi(   76 ) =    1.8898057E-10
   massi(   77 ) =    1.9454142E-10
   massi(   78 ) =    2.0019507E-10
   massi(   79 ) =    2.0594212E-10
   massi(   80 ) =    2.1178335E-10
   massi(   81 ) =    2.1771963E-10
   massi(   82 ) =    2.2375099E-10
   massi(   83 ) =    2.2987881E-10
   massi(   84 ) =    2.3610280E-10
   massi(   85 ) =    2.4242464E-10
   massi(   86 ) =    2.4884478E-10
   massi(   87 ) =    2.5536279E-10
   massi(   88 ) =    2.6198002E-10
   massi(   89 ) =    2.6869759E-10
   massi(   90 ) =    2.7551539E-10
   massi(   91 ) =    2.8243441E-10
   massi(   92 ) =    2.8945515E-10
   massi(   93 ) =    2.9657837E-10
   massi(   94 ) =    3.0380432E-10
   massi(   95 ) =    3.1113437E-10
   massi(   96 ) =    3.1856787E-10
   massi(   97 ) =    3.2610720E-10
   massi(   98 ) =    3.3375289E-10
   massi(   99 ) =    3.4150402E-10
   massi(  100 ) =    3.4936148E-10
   massi(  101 ) =    3.5732653E-10
   massi(  102 ) =    3.6540065E-10
   massi(  103 ) =    3.7358239E-10
   massi(  104 ) =    3.8187434E-10
   massi(  105 ) =    3.9027678E-10
   massi(  106 ) =    3.9878875E-10
   massi(  107 ) =    4.0741316E-10
   massi(  108 ) =    4.1614817E-10
   massi(  109 ) =    4.2499662E-10
   massi(  110 ) =    4.3395895E-10
   massi(  111 ) =    4.4303372E-10
   massi(  112 ) =    4.5222379E-10
   massi(  113 ) =    4.6152926E-10
   massi(  114 ) =    4.7095000E-10
   massi(  115 ) =    4.8048776E-10
   massi(  116 ) =    4.9013987E-10
   massi(  117 ) =    4.9991222E-10
   massi(  118 ) =    5.0980231E-10
   massi(  119 ) =    5.1981069E-10
   massi(  120 ) =    5.2994015E-10
   massi(  121 ) =    5.4018801E-10
   massi(  122 ) =    5.5055799E-10
   massi(  123 ) =    5.6104976E-10
   massi(  124 ) =    5.7166288E-10
   massi(  125 ) =    5.8239774E-10
   massi(  126 ) =    5.9325711E-10
   massi(  127 ) =    6.0423982E-10
   massi(  128 ) =    6.1534794E-10
   massi(  129 ) =    6.2657912E-10
   massi(  130 ) =    6.3793859E-10
   massi(  131 ) =    6.4942246E-10
   massi(  132 ) =    6.6103467E-10
   massi(  133 ) =    6.7277417E-10
   massi(  134 ) =    6.8464051E-10
   massi(  135 ) =    6.9663603E-10
   massi(  136 ) =    7.0876210E-10
   massi(  137 ) =    7.2101747E-10
   massi(  138 ) =    7.3340334E-10
   massi(  139 ) =    7.4592033E-10
   massi(  140 ) =    7.5856926E-10
   massi(  141 ) =    7.7135076E-10
   massi(  142 ) =    7.8426488E-10
   massi(  143 ) =    7.9730872E-10
   massi(  144 ) =    8.1049095E-10
   massi(  145 ) =    8.2380774E-10
   massi(  146 ) =    8.3725943E-10
   massi(  147 ) =    8.5084795E-10
   massi(  148 ) =    8.6457125E-10
   massi(  149 ) =    8.7843144E-10
   massi(  150 ) =    8.9242980E-10
   massi(  151 ) =    9.0656632E-10
   massi(  152 ) =    9.2084140E-10
   massi(  153 ) =    9.3525554E-10
   massi(  154 ) =    9.4981034E-10
   massi(  155 ) =    9.6450326E-10
   massi(  156 ) =    9.7933739E-10
   massi(  157 ) =    9.9431530E-10
   massi(  158 ) =    1.0094341E-09
   massi(  159 ) =    1.0246936E-09
   massi(  160 ) =    1.0400991E-09
   massi(  161 ) =    1.0556451E-09
   massi(  162 ) =    1.0713370E-09
   massi(  163 ) =    1.0871724E-09
   massi(  164 ) =    1.1031525E-09
   massi(  165 ) =    1.1192796E-09
   massi(  166 ) =    1.1355504E-09
   massi(  167 ) =    1.1519697E-09
   massi(  168 ) =    1.1685330E-09
   massi(  169 ) =    1.1852470E-09
   massi(  170 ) =    1.2021063E-09
   massi(  171 ) =    1.2191166E-09
   massi(  172 ) =    1.2362751E-09
   massi(  173 ) =    1.2535828E-09
   massi(  174 ) =    1.2710402E-09
   massi(  175 ) =    1.2886465E-09
   massi(  176 ) =    1.3063992E-09
   massi(  177 ) =    1.3243081E-09
   massi(  178 ) =    1.3423700E-09
   massi(  179 ) =    1.3605822E-09
   massi(  180 ) =    1.3789452E-09
   massi(  181 ) =    1.3974645E-09
   massi(  182 ) =    1.4161350E-09
   massi(  183 ) =    1.4349588E-09
   massi(  184 ) =    1.4539381E-09
   massi(  185 ) =    1.4730707E-09
   massi(  186 ) =    1.4923568E-09
   massi(  187 ) =    1.5118036E-09
   massi(  188 ) =    1.5314028E-09
   massi(  189 ) =    1.5511580E-09
   massi(  190 ) =    1.5710688E-09
   massi(  191 ) =    1.5911356E-09
   massi(  192 ) =    1.6113614E-09
   massi(  193 ) =    1.6317455E-09
   massi(  194 ) =    1.6522842E-09
   massi(  195 ) =    1.6729830E-09
   massi(  196 ) =    1.6938416E-09
   massi(  197 ) =    1.7148520E-09
   massi(  198 ) =    1.7360287E-09
   massi(  199 ) =    1.7573615E-09
   massi(  200 ) =    1.7788534E-09
   massi(  201 ) =    1.8005094E-09
   massi(  202 ) =    1.8223267E-09
   massi(  203 ) =    1.8442978E-09
   massi(  204 ) =    1.8664337E-09
   massi(  205 ) =    1.8887292E-09
   massi(  206 ) =    1.9111872E-09
   massi(  207 ) =    1.9338116E-09
   massi(  208 ) =    1.9565913E-09
   massi(  209 ) =    1.9795368E-09
   massi(  210 ) =    2.0026427E-09
   massi(  211 ) =    2.0259090E-09
   massi(  212 ) =    2.0493456E-09
   massi(  213 ) =    2.0729400E-09
   massi(  214 ) =    2.0966995E-09
   massi(  215 ) =    2.1206115E-09
   massi(  216 ) =    2.1446969E-09
   massi(  217 ) =    2.1689508E-09
   massi(  218 ) =    2.1933630E-09
   massi(  219 ) =    2.2179414E-09
   massi(  220 ) =    2.2426845E-09
   massi(  221 ) =    2.2676006E-09
   massi(  222 ) =    2.2926758E-09
   massi(  223 ) =    2.3179134E-09
   massi(  224 ) =    2.3433171E-09
   massi(  225 ) =    2.3688924E-09
   massi(  226 ) =    2.3946234E-09
   massi(  227 ) =    2.4205282E-09
   massi(  228 ) =    2.4465969E-09
   massi(  229 ) =    2.4728337E-09
   massi(  230 ) =    2.4992328E-09
   massi(  231 ) =    2.5258022E-09
   massi(  232 ) =    2.5525311E-09
   massi(  233 ) =    2.5794373E-09
   massi(  234 ) =    2.6065001E-09
   massi(  235 ) =    2.6337388E-09
   massi(  236 ) =    2.6611380E-09
   massi(  237 ) =    2.6887059E-09
   massi(  238 ) =    2.7164433E-09
   massi(  239 ) =    2.7443501E-09
   massi(  240 ) =    2.7724221E-09
   massi(  241 ) =    2.8006535E-09
   massi(  242 ) =    2.8290623E-09
   massi(  243 ) =    2.8576388E-09
   massi(  244 ) =    2.8863758E-09
   massi(  245 ) =    2.9152893E-09
   massi(  246 ) =    2.9443670E-09
   massi(  247 ) =    2.9736134E-09
   massi(  248 ) =    3.0030232E-09
   massi(  249 ) =    3.0325995E-09
   massi(  250 ) =    3.0623544E-09
   massi(  251 ) =    3.0922733E-09
   massi(  252 ) =    3.1223562E-09
   massi(  253 ) =    3.1526077E-09
   massi(  254 ) =    3.1830307E-09
   massi(  255 ) =    3.2136194E-09
   massi(  256 ) =    3.2443799E-09
   massi(  257 ) =    3.2753011E-09
   massi(  258 ) =    3.3063969E-09
   massi(  259 ) =    3.3376579E-09
   massi(  260 ) =    3.3690855E-09
   massi(  261 ) =    3.4006844E-09
   massi(  262 ) =    3.4324537E-09
   massi(  263 ) =    3.4643723E-09
   massi(  264 ) =    3.4964711E-09
   massi(  265 ) =    3.5287444E-09
   massi(  266 ) =    3.5611829E-09
   massi(  267 ) =    3.5937875E-09
   massi(  268 ) =    3.6265611E-09
   massi(  269 ) =    3.6594989E-09
   massi(  270 ) =    3.6926149E-09
   massi(  271 ) =    3.7258889E-09
   massi(  272 ) =    3.7593408E-09
   massi(  273 ) =    3.7929424E-09
   massi(  274 ) =    3.8267309E-09
   massi(  275 ) =    3.8606793E-09
   massi(  276 ) =    3.8947987E-09
   massi(  277 ) =    3.9290731E-09
   massi(  278 ) =    3.9635277E-09
   massi(  279 ) =    3.9981427E-09
   massi(  280 ) =    4.0329273E-09
   massi(  281 ) =    4.0678771E-09
   massi(  282 ) =    4.1029971E-09
   massi(  283 ) =    4.1382848E-09
   massi(  284 ) =    4.1737382E-09
   massi(  285 ) =    4.2093546E-09
   massi(  286 ) =    4.2451402E-09
   massi(  287 ) =    4.2810995E-09
   massi(  288 ) =    4.3172119E-09
   massi(  289 ) =    4.3534931E-09
   massi(  290 ) =    4.3899480E-09
   massi(  291 ) =    4.4265640E-09
   massi(  292 ) =    4.4633439E-09
   massi(  293 ) =    4.5002939E-09
   massi(  294 ) =    4.5374082E-09
   massi(  295 ) =    4.5746908E-09
   massi(  296 ) =    4.6121333E-09
   massi(  297 ) =    4.6497344E-09
   massi(  298 ) =    4.6875095E-09
   massi(  299 ) =    4.7254618E-09
   massi(  300 ) =    4.7635629E-09
   massi(  301 ) =    4.8018274E-09
   massi(  302 ) =    4.8402602E-09
   massi(  303 ) =    4.8788684E-09
   massi(  304 ) =    4.9176294E-09
   massi(  305 ) =    4.9565476E-09
   massi(  306 ) =    4.9956337E-09
   massi(  307 ) =    5.0348921E-09
   massi(  308 ) =    5.0743063E-09
   massi(  309 ) =    5.1138804E-09
   massi(  310 ) =    5.1536397E-09
   massi(  311 ) =    5.1935363E-09
   massi(  312 ) =    5.2336007E-09
   massi(  313 ) =    5.2738387E-09
   massi(  314 ) =    5.3142339E-09
   massi(  315 ) =    5.3548002E-09
   massi(  316 ) =    5.3955049E-09
   massi(  317 ) =    5.4363940E-09
   massi(  318 ) =    5.4774270E-09
   massi(  319 ) =    5.5186393E-09
   massi(  320 ) =    5.5600049E-09
   massi(  321 ) =    5.6015259E-09
   massi(  322 ) =    5.6432272E-09
   massi(  323 ) =    5.6850227E-09
   massi(  324 ) =    5.7270331E-09
   massi(  325 ) =    5.7692198E-09
   massi(  326 ) =    5.8115357E-09
   massi(  327 ) =    5.8540235E-09
   massi(  328 ) =    5.8966783E-09
   massi(  329 ) =    5.9394893E-09
   massi(  330 ) =    5.9824532E-09
   massi(  331 ) =    6.0255747E-09
   massi(  332 ) =    6.0688699E-09
   massi(  333 ) =    6.1123182E-09
   massi(  334 ) =    6.1559180E-09
   massi(  335 ) =    6.1996812E-09
   massi(  336 ) =    6.2436087E-09
   massi(  337 ) =    6.2876917E-09
   massi(  338 ) =    6.3319296E-09
   massi(  339 ) =    6.3763128E-09
   massi(  340 ) =    6.4208536E-09
   massi(  341 ) =    6.4655952E-09
   massi(  342 ) =    6.5104504E-09
   massi(  343 ) =    6.5554544E-09
   massi(  344 ) =    6.6006365E-09
   massi(  345 ) =    6.6459704E-09
   massi(  346 ) =    6.6914629E-09
   massi(  347 ) =    6.7371224E-09
   massi(  348 ) =    6.7829058E-09
   massi(  349 ) =    6.8288744E-09
   massi(  350 ) =    6.8749744E-09
   massi(  351 ) =    6.9212551E-09
   massi(  352 ) =    6.9676709E-09
   massi(  353 ) =    7.0142381E-09
   massi(  354 ) =    7.0609771E-09
   massi(  355 ) =    7.1078690E-09
   massi(  356 ) =    7.1548998E-09
   massi(  357 ) =    7.2020985E-09
   massi(  358 ) =    7.2494415E-09
   massi(  359 ) =    7.2969359E-09
   massi(  360 ) =    7.3446023E-09
   massi(  361 ) =    7.3923911E-09
   massi(  362 ) =    7.4403363E-09
   massi(  363 ) =    7.4884419E-09
   massi(  364 ) =    7.5367064E-09
   massi(  365 ) =    7.5851219E-09
   massi(  366 ) =    7.6336759E-09
   massi(  367 ) =    7.6823969E-09
   massi(  368 ) =    7.7312654E-09
   massi(  369 ) =    7.7802698E-09
   massi(  370 ) =    7.8294358E-09
   massi(  371 ) =    7.8787448E-09
   massi(  372 ) =    7.9282128E-09
   massi(  373 ) =    7.9778379E-09
   massi(  374 ) =    8.0276044E-09
   massi(  375 ) =    8.0775200E-09
   massi(  376 ) =    8.1275733E-09
   massi(  377 ) =    8.1777802E-09
   massi(  378 ) =    8.2281435E-09
   massi(  379 ) =    8.2786498E-09
   massi(  380 ) =    8.3292990E-09
   massi(  381 ) =    8.3801019E-09
   massi(  382 ) =    8.4310496E-09
   massi(  383 ) =    8.4821501E-09
   massi(  384 ) =    8.5333847E-09
   massi(  385 ) =    8.5847729E-09
   massi(  386 ) =    8.6363237E-09
   massi(  387 ) =    8.6880041E-09
   massi(  388 ) =    8.7398284E-09
   massi(  389 ) =    8.7918091E-09
   massi(  390 ) =    8.8439229E-09
   massi(  391 ) =    8.8961878E-09
   massi(  392 ) =    8.9485876E-09
   massi(  393 ) =    9.0011474E-09
   massi(  394 ) =    9.0538430E-09
   massi(  395 ) =    9.1066923E-09
   massi(  396 ) =    9.1596757E-09
   massi(  397 ) =    9.2127284E-09
   massi(  398 ) =    9.2660235E-09
   massi(  399 ) =    9.3194270E-09
   massi(  400 ) =    9.3729824E-09
   massi(  401 ) =    9.4266746E-09
   massi(  402 ) =    9.4805532E-09
   massi(  403 ) =    9.5345243E-09
   massi(  404 ) =    9.5886437E-09
   massi(  405 ) =    9.6429211E-09
   massi(  406 ) =    9.6973389E-09
   massi(  407 ) =    9.7518802E-09
   massi(  408 ) =    9.8065893E-09
   massi(  409 ) =    9.8614290E-09
   massi(  410 ) =    9.9163771E-09
   massi(  411 ) =    9.9715001E-09
   massi(  412 ) =    1.0026770E-08
   massi(  413 ) =    1.0082156E-08
   massi(  414 ) =    1.0137732E-08
   massi(  415 ) =    1.0193389E-08
   massi(  416 ) =    1.0249217E-08
   massi(  417 ) =    1.0305138E-08
   massi(  418 ) =    1.0361234E-08
   massi(  419 ) =    1.0417470E-08
   massi(  420 ) =    1.0473836E-08
   massi(  421 ) =    1.0530332E-08
   massi(  422 ) =    1.0586974E-08
   massi(  423 ) =    1.0643731E-08
   massi(  424 ) =    1.0700678E-08
   massi(  425 ) =    1.0757740E-08
   massi(  426 ) =    1.0814914E-08
   massi(  427 ) =    1.0872250E-08
   massi(  428 ) =    1.0929705E-08
   massi(  429 ) =    1.0987320E-08
   massi(  430 ) =    1.1045032E-08
   massi(  431 ) =    1.1102897E-08
   massi(  432 ) =    1.1160925E-08
   massi(  433 ) =    1.1219035E-08
   massi(  434 ) =    1.1277319E-08
   massi(  435 ) =    1.1335743E-08
   massi(  436 ) =    1.1394302E-08
   massi(  437 ) =    1.1452972E-08
   massi(  438 ) =    1.1511775E-08
   massi(  439 ) =    1.1570734E-08
   massi(  440 ) =    1.1629776E-08
   massi(  441 ) =    1.1688992E-08
   massi(  442 ) =    1.1748331E-08
   massi(  443 ) =    1.1807789E-08
   massi(  444 ) =    1.1867404E-08
   massi(  445 ) =    1.1927166E-08
   massi(  446 ) =    1.1987026E-08
   massi(  447 ) =    1.2047035E-08
   massi(  448 ) =    1.2107157E-08
   massi(  449 ) =    1.2167435E-08
   massi(  450 ) =    1.2227797E-08
   massi(  451 ) =    1.2288330E-08
   massi(  452 ) =    1.2348966E-08
   massi(  453 ) =    1.2409811E-08
   massi(  454 ) =    1.2470704E-08
   massi(  455 ) =    1.2531734E-08
   massi(  456 ) =    1.2592886E-08
   massi(  457 ) =    1.2654191E-08
   massi(  458 ) =    1.2715618E-08
   massi(  459 ) =    1.2777200E-08
   massi(  460 ) =    1.2838903E-08
   massi(  461 ) =    1.2900722E-08
   massi(  462 ) =    1.2962626E-08
   massi(  463 ) =    1.3024688E-08
   massi(  464 ) =    1.3086929E-08
   massi(  465 ) =    1.3149250E-08
   massi(  466 ) =    1.3211729E-08
   massi(  467 ) =    1.3274268E-08
   massi(  468 ) =    1.3337013E-08
   massi(  469 ) =    1.3399824E-08
   massi(  470 ) =    1.3462774E-08
   massi(  471 ) =    1.3525849E-08
   massi(  472 ) =    1.3589078E-08
   massi(  473 ) =    1.3652425E-08
   massi(  474 ) =    1.3715887E-08
   massi(  475 ) =    1.3779463E-08
   massi(  476 ) =    1.3843160E-08
   massi(  477 ) =    1.3907028E-08
   massi(  478 ) =    1.3970989E-08
   massi(  479 ) =    1.4035066E-08
   massi(  480 ) =    1.4099279E-08
   massi(  481 ) =    1.4163609E-08
   massi(  482 ) =    1.4228044E-08
   massi(  483 ) =    1.4292671E-08
   massi(  484 ) =    1.4357329E-08
   massi(  485 ) =    1.4422177E-08
   massi(  486 ) =    1.4487101E-08
   massi(  487 ) =    1.4552209E-08
   massi(  488 ) =    1.4617378E-08
   massi(  489 ) =    1.4682700E-08
   massi(  490 ) =    1.4748133E-08
   massi(  491 ) =    1.4813669E-08
   massi(  492 ) =    1.4879225E-08
   massi(  493 ) =    1.4944977E-08
   massi(  494 ) =    1.5010889E-08
   massi(  495 ) =    1.5076957E-08
   massi(  496 ) =    1.5143138E-08
   massi(  497 ) =    1.5209405E-08
   massi(  498 ) =    1.5275784E-08
   massi(  499 ) =    1.5342318E-08
   massi(  500 ) =    1.5408999E-08
   massi(  501 ) =    1.5475692E-08
   massi(  502 ) =    1.5542575E-08
   massi(  503 ) =    1.5609567E-08
   massi(  504 ) =    1.5676706E-08
   massi(  505 ) =    1.5743975E-08
   massi(  506 ) =    1.5811302E-08
   massi(  507 ) =    1.5878754E-08
   massi(  508 ) =    1.5946384E-08
   massi(  509 ) =    1.6014077E-08
   massi(  510 ) =    1.6081893E-08
   massi(  511 ) =    1.6149848E-08
   massi(  512 ) =    1.6217907E-08
   massi(  513 ) =    1.6286052E-08
   massi(  514 ) =    1.6354392E-08
   massi(  515 ) =    1.6422742E-08
   massi(  516 ) =    1.6491290E-08
   massi(  517 ) =    1.6559977E-08
   massi(  518 ) =    1.6628723E-08
   massi(  519 ) =    1.6697577E-08
   massi(  520 ) =    1.6766595E-08
   massi(  521 ) =    1.6835671E-08
   massi(  522 ) =    1.6904906E-08
   massi(  523 ) =    1.6974205E-08
   massi(  524 ) =    1.7043652E-08
   massi(  525 ) =    1.7113285E-08
   massi(  526 ) =    1.7182890E-08
   massi(  527 ) =    1.7252669E-08
   massi(  528 ) =    1.7322598E-08
   massi(  529 ) =    1.7392617E-08
   massi(  530 ) =    1.7462758E-08
   massi(  531 ) =    1.7533024E-08
   massi(  532 ) =    1.7603353E-08
   massi(  533 ) =    1.7673875E-08
   massi(  534 ) =    1.7744417E-08
   massi(  535 ) =    1.7815156E-08
   massi(  536 ) =    1.7885904E-08
   massi(  537 ) =    1.7956822E-08
   massi(  538 ) =    1.8027855E-08
   massi(  539 ) =    1.8098987E-08
   massi(  540 ) =    1.8170303E-08
   massi(  541 ) =    1.8241652E-08
   massi(  542 ) =    1.8313093E-08
   massi(  543 ) =    1.8384727E-08
   massi(  544 ) =    1.8456356E-08
   massi(  545 ) =    1.8528258E-08
   massi(  546 ) =    1.8600177E-08
   massi(  547 ) =    1.8672164E-08
   massi(  548 ) =    1.8744331E-08
   massi(  549 ) =    1.8816577E-08
   massi(  550 ) =    1.8888937E-08
   massi(  551 ) =    1.8961401E-08
   massi(  552 ) =    1.9034021E-08
   massi(  553 ) =    1.9106697E-08
   massi(  554 ) =    1.9179451E-08
   massi(  555 ) =    1.9252376E-08
   massi(  556 ) =    1.9325418E-08
   massi(  557 ) =    1.9398568E-08
   massi(  558 ) =    1.9471759E-08
   massi(  559 ) =    1.9545130E-08
   massi(  560 ) =    1.9618584E-08
   massi(  561 ) =    1.9692115E-08
   massi(  562 ) =    1.9765771E-08
   massi(  563 ) =    1.9839536E-08
   massi(  564 ) =    1.9913474E-08
   massi(  565 ) =    1.9987478E-08
   massi(  566 ) =    2.0061542E-08
   massi(  567 ) =    2.0135758E-08
   massi(  568 ) =    2.0210075E-08
   massi(  569 ) =    2.0284464E-08
   massi(  570 ) =    2.0358991E-08
   massi(  571 ) =    2.0433639E-08
   massi(  572 ) =    2.0508352E-08
   massi(  573 ) =    2.0583231E-08
   massi(  574 ) =    2.0658169E-08
   massi(  575 ) =    2.0733193E-08
   massi(  576 ) =    2.0808310E-08
   massi(  577 ) =    2.0883606E-08
   massi(  578 ) =    2.0958950E-08
   massi(  579 ) =    2.1034467E-08
   massi(  580 ) =    2.1110036E-08
   massi(  581 ) =    2.1185658E-08
   massi(  582 ) =    2.1261467E-08
   massi(  583 ) =    2.1337327E-08
   massi(  584 ) =    2.1413365E-08
   massi(  585 ) =    2.1489400E-08
   massi(  586 ) =    2.1565665E-08
   massi(  587 ) =    2.1641922E-08
   massi(  588 ) =    2.1718378E-08
   massi(  589 ) =    2.1794865E-08
   massi(  590 ) =    2.1871474E-08
   massi(  591 ) =    2.1948180E-08
   massi(  592 ) =    2.2025006E-08
   massi(  593 ) =    2.2101931E-08
   massi(  594 ) =    2.2179002E-08
   massi(  595 ) =    2.2256090E-08
   massi(  596 ) =    2.2333300E-08
   massi(  597 ) =    2.2410630E-08
   massi(  598 ) =    2.2488042E-08
   massi(  599 ) =    2.2565578E-08
   massi(  600 ) =    2.2643235E-08
   massi(  601 ) =    2.2720933E-08
   massi(  602 ) =    2.2798814E-08
   massi(  603 ) =    2.2876751E-08
   massi(  604 ) =    2.2954758E-08
   massi(  605 ) =    2.3032911E-08
   massi(  606 ) =    2.3111163E-08
   massi(  607 ) =    2.3189520E-08
   massi(  608 ) =    2.3267940E-08
   massi(  609 ) =    2.3346507E-08
   massi(  610 ) =    2.3425113E-08
   massi(  611 ) =    2.3503832E-08
   massi(  612 ) =    2.3582677E-08
   massi(  613 ) =    2.3661615E-08
   massi(  614 ) =    2.3740339E-08
   massi(  615 ) =    2.3819487E-08
   massi(  616 ) =    2.3898714E-08
   massi(  617 ) =    2.3978041E-08
   massi(  618 ) =    2.4057455E-08
   massi(  619 ) =    2.4136989E-08
   massi(  620 ) =    2.4216659E-08
   massi(  621 ) =    2.4296360E-08
   massi(  622 ) =    2.4376206E-08
   massi(  623 ) =    2.4456114E-08
   massi(  624 ) =    2.4536099E-08
   massi(  625 ) =    2.4616270E-08
   massi(  626 ) =    2.4696469E-08
   massi(  627 ) =    2.4776815E-08
   massi(  628 ) =    2.4857247E-08
   massi(  629 ) =    2.4937739E-08
   massi(  630 ) =    2.5018336E-08
   massi(  631 ) =    2.5099034E-08
   massi(  632 ) =    2.5179883E-08
   massi(  633 ) =    2.5260769E-08
   massi(  634 ) =    2.5341711E-08
   massi(  635 ) =    2.5422869E-08
   massi(  636 ) =    2.5504017E-08
   massi(  637 ) =    2.5585289E-08
   massi(  638 ) =    2.5666717E-08
   massi(  639 ) =    2.5748111E-08
   massi(  640 ) =    2.5829779E-08
   massi(  641 ) =    2.5911408E-08
   massi(  642 ) =    2.5993153E-08
   massi(  643 ) =    2.6075009E-08
   massi(  644 ) =    2.6156988E-08
   massi(  645 ) =    2.6239023E-08
   massi(  646 ) =    2.6321157E-08
   massi(  647 ) =    2.6403361E-08
   massi(  648 ) =    2.6485687E-08
   massi(  649 ) =    2.6568213E-08
   massi(  650 ) =    2.6650666E-08
   massi(  651 ) =    2.6733256E-08
   massi(  652 ) =    2.6816013E-08
   massi(  653 ) =    2.6898817E-08
   massi(  654 ) =    2.6981731E-08
   massi(  655 ) =    2.7064740E-08
   massi(  656 ) =    2.7147749E-08
   massi(  657 ) =    2.7230971E-08
   massi(  658 ) =    2.7314231E-08
   massi(  659 ) =    2.7397585E-08
   massi(  660 ) =    2.7481116E-08
   massi(  661 ) =    2.7564630E-08
   massi(  662 ) =    2.7648234E-08
   massi(  663 ) =    2.7731994E-08
   massi(  664 ) =    2.7815881E-08
   massi(  665 ) =    2.7899720E-08
   massi(  666 ) =    2.7983793E-08
   massi(  667 ) =    2.8067882E-08
   massi(  668 ) =    2.8152073E-08
   massi(  669 ) =    2.8236350E-08
   massi(  670 ) =    2.8320731E-08
   massi(  671 ) =    2.8405193E-08
   massi(  672 ) =    2.8489806E-08
   massi(  673 ) =    2.8574528E-08
   massi(  674 ) =    2.8659143E-08
   massi(  675 ) =    2.8744079E-08
   massi(  676 ) =    2.8829035E-08
   massi(  677 ) =    2.8914005E-08
   massi(  678 ) =    2.8999102E-08
   massi(  679 ) =    2.9084379E-08
   massi(  680 ) =    2.9169605E-08
   massi(  681 ) =    2.9254977E-08
   massi(  682 ) =    2.9340496E-08
   massi(  683 ) =    2.9426065E-08
   massi(  684 ) =    2.9511698E-08
   massi(  685 ) =    2.9597485E-08
   massi(  686 ) =    2.9683305E-08
   massi(  687 ) =    2.9769200E-08
   massi(  688 ) =    2.9855254E-08
   massi(  689 ) =    2.9941404E-08
   massi(  690 ) =    3.0027529E-08
   massi(  691 ) =    3.0113839E-08
   massi(  692 ) =    3.0200237E-08
   massi(  693 ) =    3.0286696E-08
   massi(  694 ) =    3.0373219E-08
   massi(  695 ) =    3.0459926E-08
   massi(  696 ) =    3.0546630E-08
   massi(  697 ) =    3.0633437E-08
   massi(  698 ) =    3.0720333E-08
   massi(  699 ) =    3.0807396E-08
   massi(  700 ) =    3.0894391E-08
   massi(  701 ) =    3.0981713E-08
   massi(  702 ) =    3.1068854E-08
   massi(  703 ) =    3.1156297E-08
   massi(  704 ) =    3.1243633E-08
   massi(  705 ) =    3.1331258E-08
   massi(  706 ) =    3.1418832E-08
   massi(  707 ) =    3.1506545E-08
   massi(  708 ) =    3.1594261E-08
   massi(  709 ) =    3.1682259E-08
   massi(  710 ) =    3.1770188E-08
   massi(  711 ) =    3.1858221E-08
   massi(  712 ) =    3.1946438E-08
   massi(  713 ) =    3.2034585E-08
   massi(  714 ) =    3.2122934E-08
   massi(  715 ) =    3.2211297E-08
   massi(  716 ) =    3.2299827E-08
   massi(  717 ) =    3.2388439E-08
   massi(  718 ) =    3.2477054E-08
   massi(  719 ) =    3.2565783E-08
   massi(  720 ) =    3.2654611E-08
   massi(  721 ) =    3.2743579E-08
   massi(  722 ) =    3.2832610E-08
   massi(  723 ) =    3.2921644E-08
   massi(  724 ) =    3.3010888E-08
   massi(  725 ) =    3.3100171E-08
   massi(  726 ) =    3.3189504E-08
   massi(  727 ) =    3.3278919E-08
   massi(  728 ) =    3.3368490E-08
   massi(  729 ) =    3.3458115E-08
   massi(  730 ) =    3.3547710E-08
   massi(  731 ) =    3.3637512E-08
   massi(  732 ) =    3.3727410E-08
   massi(  733 ) =    3.3817393E-08
   massi(  734 ) =    3.3907394E-08
   massi(  735 ) =    3.3997544E-08
   massi(  736 ) =    3.4087705E-08
   massi(  737 ) =    3.4178090E-08
   massi(  738 ) =    3.4268414E-08
   massi(  739 ) =    3.4358852E-08
   massi(  740 ) =    3.4449407E-08
   massi(  741 ) =    3.4540093E-08
   massi(  742 ) =    3.4630812E-08
   massi(  743 ) =    3.4721509E-08
   massi(  744 ) =    3.4812444E-08
   massi(  745 ) =    3.4903405E-08
   massi(  746 ) =    3.4994464E-08
   massi(  747 ) =    3.5085584E-08
   massi(  748 ) =    3.5176850E-08
   massi(  749 ) =    3.5268087E-08
   massi(  750 ) =    3.5359584E-08
   massi(  751 ) =    3.5451045E-08
   massi(  752 ) =    3.5542509E-08
   massi(  753 ) =    3.5634145E-08
   massi(  754 ) =    3.5725900E-08
   massi(  755 ) =    3.5817688E-08
   massi(  756 ) =    3.5909526E-08
   massi(  757 ) =    3.6001563E-08
   massi(  758 ) =    3.6093589E-08
   massi(  759 ) =    3.6185750E-08
   massi(  760 ) =    3.6278028E-08
   massi(  761 ) =    3.6370274E-08
   massi(  762 ) =    3.6462648E-08
   massi(  763 ) =    3.6555161E-08
   massi(  764 ) =    3.6647716E-08
   massi(  765 ) =    3.6740303E-08
   massi(  766 ) =    3.6832965E-08
   massi(  767 ) =    3.6925840E-08
   massi(  768 ) =    3.7018676E-08
   massi(  769 ) =    3.7111629E-08
   massi(  770 ) =    3.7204661E-08
   massi(  771 ) =    3.7297720E-08
   massi(  772 ) =    3.7391004E-08
   massi(  773 ) =    3.7484245E-08
   massi(  774 ) =    3.7577681E-08
   massi(  775 ) =    3.7670368E-08
   massi(  776 ) =    3.7763957E-08
   massi(  777 ) =    3.7857561E-08
   massi(  778 ) =    3.7951256E-08
   massi(  779 ) =    3.8045041E-08
   massi(  780 ) =    3.8138875E-08
   massi(  781 ) =    3.8232873E-08
   massi(  782 ) =    3.8326753E-08
   massi(  783 ) =    3.8420851E-08
   massi(  784 ) =    3.8515118E-08
   massi(  785 ) =    3.8609436E-08
   massi(  786 ) =    3.8703657E-08
   massi(  787 ) =    3.8798092E-08
   massi(  788 ) =    3.8892601E-08
   massi(  789 ) =    3.8987238E-08
   massi(  790 ) =    3.9081879E-08
   massi(  791 ) =    3.9176619E-08
   massi(  792 ) =    3.9271423E-08
   massi(  793 ) =    3.9366419E-08
   massi(  794 ) =    3.9461320E-08
   massi(  795 ) =    3.9556372E-08
   massi(  796 ) =    3.9651514E-08
   massi(  797 ) =    3.9746766E-08
   massi(  798 ) =    3.9842064E-08
   massi(  799 ) =    3.9937422E-08
   massi(  800 ) =    4.0032862E-08
   massi(  801 ) =    4.0128427E-08
   massi(  802 ) =    4.0224098E-08
   massi(  803 ) =    4.0319762E-08
   massi(  804 ) =    4.0415635E-08
   massi(  805 ) =    4.0511473E-08
   massi(  806 ) =    4.0607375E-08
   massi(  807 ) =    4.0703323E-08
   massi(  808 ) =    4.0799382E-08
   massi(  809 ) =    4.0895547E-08
   massi(  810 ) =    4.0991765E-08
   massi(  811 ) =    4.1088192E-08
   massi(  812 ) =    4.1184663E-08
   massi(  813 ) =    4.1281055E-08
   massi(  814 ) =    4.1377678E-08
   massi(  815 ) =    4.1474276E-08
   massi(  816 ) =    4.1570974E-08
   massi(  817 ) =    4.1667885E-08
   massi(  818 ) =    4.1764700E-08
   massi(  819 ) =    4.1861750E-08
   massi(  820 ) =    4.1958721E-08
   massi(  821 ) =    4.2055721E-08
   massi(  822 ) =    4.2152905E-08
   massi(  823 ) =    4.2250129E-08
   massi(  824 ) =    4.2347544E-08
   massi(  825 ) =    4.2444999E-08
   massi(  826 ) =    4.2542474E-08
   massi(  827 ) =    4.2640089E-08
   massi(  828 ) =    4.2737732E-08
   massi(  829 ) =    4.2835463E-08
   massi(  830 ) =    4.2933227E-08
   massi(  831 ) =    4.3031175E-08
   massi(  832 ) =    4.3129280E-08
   massi(  833 ) =    4.3227384E-08
   massi(  834 ) =    4.3325382E-08
   massi(  835 ) =    4.3423409E-08
   massi(  836 ) =    4.3521691E-08
   massi(  837 ) =    4.3620062E-08
   massi(  838 ) =    4.3718508E-08
   massi(  839 ) =    4.3817053E-08
   massi(  840 ) =    4.3915634E-08
   massi(  841 ) =    4.4014314E-08
   massi(  842 ) =    4.4113079E-08
   massi(  843 ) =    4.4211848E-08
   massi(  844 ) =    4.4310752E-08
   massi(  845 ) =    4.4409703E-08
   massi(  846 ) =    4.4508749E-08
   massi(  847 ) =    4.4607816E-08
   massi(  848 ) =    4.4707075E-08
   massi(  849 ) =    4.4806420E-08
   massi(  850 ) =    4.4905740E-08
   massi(  851 ) =    4.5005137E-08
   massi(  852 ) =    4.5104553E-08
   massi(  853 ) =    4.5204047E-08
   massi(  854 ) =    4.5303715E-08
   massi(  855 ) =    4.5403485E-08
   massi(  856 ) =    4.5503281E-08
   massi(  857 ) =    4.5603151E-08
   massi(  858 ) =    4.5703104E-08
   massi(  859 ) =    4.5803123E-08
   massi(  860 ) =    4.5903267E-08
   massi(  861 ) =    4.6003400E-08
   massi(  862 ) =    4.6103686E-08
   massi(  863 ) =    4.6204004E-08
   massi(  864 ) =    4.6304375E-08
   massi(  865 ) =    4.6404846E-08
   massi(  866 ) =    4.6505445E-08
   massi(  867 ) =    4.6606079E-08
   massi(  868 ) =    4.6706791E-08
   massi(  869 ) =    4.6807560E-08
   massi(  870 ) =    4.6908447E-08
   massi(  871 ) =    4.7009340E-08
   massi(  872 ) =    4.7110394E-08
   massi(  873 ) =    4.7211440E-08
   massi(  874 ) =    4.7312692E-08
   massi(  875 ) =    4.7413739E-08
   massi(  876 ) =    4.7515059E-08
   massi(  877 ) =    4.7616425E-08
   massi(  878 ) =    4.7717819E-08
   massi(  879 ) =    4.7819395E-08
   massi(  880 ) =    4.7920931E-08
   massi(  881 ) =    4.8022681E-08
   massi(  882 ) =    4.8124367E-08
   massi(  883 ) =    4.8226166E-08
   massi(  884 ) =    4.8328111E-08
   massi(  885 ) =    4.8430067E-08
   massi(  886 ) =    4.8532030E-08
   massi(  887 ) =    4.8634167E-08
   massi(  888 ) =    4.8736300E-08
   massi(  889 ) =    4.8838565E-08
   massi(  890 ) =    4.8940954E-08
   massi(  891 ) =    4.9043372E-08
   massi(  892 ) =    4.9145868E-08
   massi(  893 ) =    4.9248335E-08
   massi(  894 ) =    4.9351023E-08
   massi(  895 ) =    4.9453732E-08
   massi(  896 ) =    4.9556434E-08
   massi(  897 ) =    4.9659420E-08
   massi(  898 ) =    4.9762278E-08
   massi(  899 ) =    4.9865207E-08
   massi(  900 ) =    4.9968360E-08
   massi(  901 ) =    5.0071407E-08
   massi(  902 ) =    5.0174634E-08
   massi(  903 ) =    5.0277976E-08
   massi(  904 ) =    5.0381271E-08
   massi(  905 ) =    5.0484765E-08
   massi(  906 ) =    5.0588191E-08
   massi(  907 ) =    5.0691781E-08
   massi(  908 ) =    5.0795474E-08
   massi(  909 ) =    5.0899231E-08
   massi(  910 ) =    5.1002988E-08
   massi(  911 ) =    5.1106870E-08
   massi(  912 ) =    5.1210783E-08
   massi(  913 ) =    5.1314874E-08
   massi(  914 ) =    5.1418819E-08
   massi(  915 ) =    5.1523028E-08
   massi(  916 ) =    5.1627371E-08
   massi(  917 ) =    5.1731600E-08
   massi(  918 ) =    5.1836018E-08
   massi(  919 ) =    5.1940436E-08
   massi(  920 ) =    5.2044960E-08
   massi(  921 ) =    5.2149552E-08
   massi(  922 ) =    5.2254158E-08
   massi(  923 ) =    5.2358999E-08
   massi(  924 ) =    5.2463790E-08
   massi(  925 ) =    5.2568737E-08
   massi(  926 ) =    5.2673649E-08
   massi(  927 ) =    5.2778603E-08
   massi(  928 ) =    5.2883745E-08
   massi(  929 ) =    5.2988870E-08
   massi(  930 ) =    5.3094134E-08
   massi(  931 ) =    5.3199418E-08
   massi(  932 ) =    5.3304742E-08
   massi(  933 ) =    5.3410169E-08
   massi(  934 ) =    5.3515766E-08
   massi(  935 ) =    5.3621271E-08
   massi(  936 ) =    5.3726943E-08
   massi(  937 ) =    5.3832682E-08
   massi(  938 ) =    5.3938493E-08
   massi(  939 ) =    5.4044428E-08
   massi(  940 ) =    5.4150284E-08
   massi(  941 ) =    5.4256400E-08
   massi(  942 ) =    5.4362435E-08
   massi(  943 ) =    5.4468597E-08
   massi(  944 ) =    5.4574869E-08
   massi(  945 ) =    5.4681184E-08
   massi(  946 ) =    5.4787471E-08
   massi(  947 ) =    5.4893892E-08
   massi(  948 ) =    5.5000399E-08
   massi(  949 ) =    5.5107083E-08
   massi(  950 ) =    5.5213903E-08
   massi(  951 ) =    5.5320509E-08
   massi(  952 ) =    5.5427336E-08
   massi(  953 ) =    5.5534056E-08
   massi(  954 ) =    5.5641081E-08
   massi(  955 ) =    5.5748060E-08
   massi(  956 ) =    5.5855192E-08
   massi(  957 ) =    5.5962253E-08
   massi(  958 ) =    5.6069453E-08
   massi(  959 ) =    5.6176731E-08
   massi(  960 ) =    5.6284080E-08
   massi(  961 ) =    5.6391539E-08
   massi(  962 ) =    5.6499022E-08
   massi(  963 ) =    5.6606719E-08
   massi(  964 ) =    5.6714217E-08
   massi(  965 ) =    5.6822028E-08
   massi(  966 ) =    5.6929739E-08
   massi(  967 ) =    5.7037724E-08
   massi(  968 ) =    5.7145549E-08
   massi(  969 ) =    5.7253594E-08
   massi(  970 ) =    5.7361625E-08
   massi(  971 ) =    5.7469709E-08
   massi(  972 ) =    5.7577946E-08
   massi(  973 ) =    5.7686151E-08
   massi(  974 ) =    5.7794505E-08
   massi(  975 ) =    5.7902895E-08
   massi(  976 ) =    5.8011320E-08
   massi(  977 ) =    5.8119884E-08
   massi(  978 ) =    5.8228455E-08
   massi(  979 ) =    5.8337179E-08
   massi(  980 ) =    5.8445963E-08
   massi(  981 ) =    5.8554871E-08
   massi(  982 ) =    5.8663716E-08
   massi(  983 ) =    5.8772603E-08
   massi(  984 ) =    5.8881664E-08
   massi(  985 ) =    5.8990707E-08
   massi(  986 ) =    5.9098593E-08
   massi(  987 ) =    5.9207657E-08
   massi(  988 ) =    5.9317003E-08
   massi(  989 ) =    5.9426387E-08
   massi(  990 ) =    5.9535818E-08
   massi(  991 ) =    5.9645508E-08
   massi(  992 ) =    5.9755060E-08
   massi(  993 ) =    5.9864696E-08
   massi(  994 ) =    5.9974383E-08
   massi(  995 ) =    6.0084169E-08
   massi(  996 ) =    6.0194033E-08
   massi(  997 ) =    6.0303890E-08
   massi(  998 ) =    6.0413910E-08
   massi(  999 ) =    6.0523995E-08
   massi( 1000 ) =    6.0634157E-08

  RETURN
END SUBROUTINE init_massr_massi

!#######################################################################
!
! Borrowed and modified from COAMPS diagnostic package.
!
! by Yunheng Wang on 05/11/2011.
!
!#######################################################################
!
SUBROUTINE coamps_rrf(tz,press,qr,qs,qi,qc,qg,qv,m,n,kk,radr) !,radc)
!
!********0*********0*********0*********0*********0*********0*********
!   Now, calculate radar reflectivity based on the following
!        information:
!        The formulas used here are from Louis J. Battan (1981)
!            Z(cloud) = 4.8*10**(-2)*M**(2.0)        (Atlas, 1954)
!            Z(ice)   = 3.8*10**(-1)*M**(2.2)
!            Z(rain)  = 2.4*10**(+4)*M**(1.82)       (Douglas, 1964)
!            Z(snow)  = 3.8*10**(+4)*M**(2.2)        (Douglas, 1964)
!            Z(graupel,wet,r=3cm)
!                     = 8.0*10**(+5)*M**(0.98)       (Douglas, 1964)
!            Z(graupel,dry,r=3cm)
!                     = 2.6*10**(+5)*M**(1.03)       (Douglas, 1964)
!        where Z is the radar reflectivity factor and in (mm**6/m**3)
!              M is the water content and in (g/m**3)
!        then radar reflectivity ref = 10*log(Z)  is of unit of (dbZ)
!********0*********0*********0*********0*********0*********0*********
!

      IMPLICIT NONE

      integer, intent(in) :: m, n, kk

      real,    intent(in) :: tz(m,n,kk)
      real,    intent(in) :: press(m,n,kk)
      real,    intent(in) :: qr(m,n,kk)
      real,    intent(in) :: qs(m,n,kk)
      real,    intent(in) :: qi(m,n,kk)
      real,    intent(in) :: qc(m,n,kk)
      real,    intent(in) :: qg(m,n,kk)
      real,    intent(in) :: qv(m,n,kk)

      real,   intent(out) :: radr(m,n,kk)
      !real,   intent(out) :: radc(m,n)

!-----------------------------------------------------------------------

      integer :: i,j,k

      real    :: qcm, qrm, qim, qsm, qgm

      real    :: vovermd

      real    :: tt, tv
      real    :: rhoair, rhowat, rhoice

      !include 'constant.h'
      real, parameter :: p00  = 1.e5
      real, parameter :: rgas = 287.04
      real, parameter :: cp   = 1004.64
      real, parameter :: rocp = rgas/cp

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  !
  !  compute radar reflectivity
  !
  rhowat=1000.0
  rhoice=917.0

  do k=1,kk
    do j=1,n
      do i=1,m
        !tt=th(i,1,k)*((press(i,j,k)/p00)**rocp)
        tt = tz(i,j,k)
        tv=tt*(1.0+0.61*qv(i,j,k))
        rhoair=press(i,j,k)/(rgas*tv)

        vovermd=( (1.+qv(i,j,k))/rhoair + (qc(i,j,k)+qr(i,j,k))/rhowat  &
                 +(qi(i,j,k)+qs(i,j,k)+qg(i,j,k))/rhoice)*0.001

        qcm = qc(i,j,k)/vovermd
        qrm = qr(i,j,k)/vovermd
        qim = qi(i,j,k)/vovermd
        qsm = qs(i,j,k)/vovermd
        qgm = qg(i,j,k)/vovermd
        if(tt >= 273.0) then
           radr(i,j,k)=0.04800*qcm**2.0 +0.38000*qim**2.2               &
                      +24000.0*qrm**1.82+38000.0*qsm**2.2               &
                      +24000.0*qgm**1.82
        else
           radr(i,j,k)=0.04800*qcm**2.0 +0.38000*qim**2.2               &
                      +24000.0*qrm**1.82+38000.0*qsm**2.2               &
                      +38000.0*qgm**2.2
        endif
        if(radr(i,j,k) >= 1.0e-25) then
           radr(i,j,k)=10.0*alog10(radr(i,j,k))
        else
           radr(i,j,k)=0.0
        endif
      enddo
    enddo
  enddo
  !
  ! get composite scan
  !
  !do i=1,m
  !  do j=1,n
  !    radc(i,j)=0.0
  !    do k=1,kk
  !      if (radr(i,j,k) >= radc(i,j)) radc(i,j)=radr(i,j,k)
  !    enddo
  !  enddo
  !enddo

  RETURN
END SUBROUTINE coamps_rrf
