!
!##################################################################
!##################################################################
!######                                                      ######
!######                   SUBROUTINE KFPARA                  ######
!######                                                      ######
!######                       Adapted by                     ######
!######          Coastal Meteorology Research Program        ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE kfpara (nx,ny,nz,j,dt,dx,mphyopt,kffbfct,               & 
           kfsubsattrig,                                           & 
           ptop,psb,sigma,a,ub,vb,w0avg,tb,qvb,                    & 
           dtdt,dqdt,dqldt,dqrdt,dqidt,dqsdt,                      &
           nca,raincv,ncuyes,icuyes,lsb)
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                                                                   C
!     THIS SUBROUTINE COMPUTES THE EFFECTS OF DEEP CONVECTION USING C
!   THE KAIN-FRITSCH CONVECTIVE PARAMETERIZATION SCHEME.            C
!                                              MKS UNITS.           C
!                                                                   C
!  INPUT:   TEMPERTURE (T0, K) ;    SPECIFIC HUMIDITY (Q0, KG/KG) ; C
!           HORIZONTAL WIND SPEED (U0 AND V0, M/S) ;                C
!           PRESSURE (P0, PASCAL) ;  VERTICAL MOTION (W0AVG, M/S).  C
!  OUTPUT:  CONVECTIVE TEMPERATURE (DTDT), WATER VAPOR (DQDT),      C
!           CLOUD WATER (DQLDT), CLOUD ICE (DQIDT),                 C
!           RAIN WATER (DQRDT),  SNOW (DQSDT), AND                  C
!           RAINFALL (RAINCV) TENDENCIES.                           C
!                                                                   C
!                                DOCUMENTED BY JACK KAIN            C
!                                JANUARY 1995                       C
!                                UPDATED SEPTEMBER 1995             C
!                                        NOVEMBER 1996              C
!                                                                   C
!  REFERENCES:                                                      C
!                                                                   C
!      KAIN AND FRITSCH (1993):  "CONVECTIVE PARAMETERIZATION IN    C
!      MESOSCALE MODELS:  THE KAIN-FRITSCH SCHEME" IN THE REPRESEN- C
!      TATION OF CUMULUS CONVECTION IN NUMERICAL MODELS, A.M.S.     C
!      MONOGRAPH, K.A. EMANUEL AND D.J. RAYMOND, EDS., 165-170.     C
!                                                                   C
!      FRITSCH AND KAIN (1993):  "CONVECTIVE PARAMETERIZATION IN    C
!      MESOSCALE MODELS:  THE FRITSCH-CHAPPELL SCHEME" IN THE REP-  C
!      RESENTATION OF CUMULUS CONVECTION IN NUMERICAL MODELS, A.M.S.C
!      MONOGRAPH, K.A. EMANUEL AND D.J. RAYMOND, EDS., 165-170.     C
!                                                                   C
!      FRITSCH AND CHAPPELL (1980), J. ATMOS. SCI., 1722-1733.      C
!                                                                   C
! NOTE:  THE KAIN-FRITSCH SCHEME IS UNDER CONTINUED DEVELOPMENT.    C
!     IF YOU WOULD LIKE TO BE NOTIFIED WHEN THERE ARE UPDATES OR    C
!     CHANGES TO THE SCHEME, OR IF YOU HAVE PROBLEMS, QUESTIONS, OR C
!     SUGGESTIONS PLEASE NOTIFY ME AT                               C
!                                                                   C
!                               KAIN@ESSC.PSU.EDU                   C
!                                                                   C
!     JACK KAIN                                                     C
!     JANUARY 1995                                                  C
!     November 1996                                                 C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!    This subrooutine is implemented in the ARPS by Zonghui Huo,
!    08/01/97
!
!    11/21/2001  (Yunheng Wang)
!    Corrected error u0(k) v0(k) should be  column value and no average
!    needed again in this subroutine.
!
!    04/18/2002  (Zuwen He) 
! 
!    Allow sub-saturation in KF cumulus parameterization, see below 
!    for detail. Add a pass-in argument, kfsubsattrig, 
!    from the arps.input. If kfsubsattrig=1, the trigger is on. 
!
!    7/26/2002 (Fanyou Kong)
!    Rename subroutine dtfrznew, envirtht, condload, and prof5 to
!    dtfrznew_old, envirtht_old, condload_old, and prof5_old,
!    respectively, to avoid ambiguity with the new WRF KF code 
!
!-----------------------------------------------------------------------
!
!  Variable Declarations:
!
!-----------------------------------------------------------------------
!
  PARAMETER (kmx=100)  ! increase if vertical levels >100

  DIMENSION ptop(nx,ny)
  DIMENSION sigma(nx,ny,nz), a(nx,ny,nz)

  DIMENSION ub(nx,ny,nz),vb(nx,ny,nz),w0avg(nx,ny,nz)
  DIMENSION tb(nx,ny,nz),qvb(nx,ny,nz)
  DIMENSION psb(nx,ny)

  DIMENSION dtdt(nx,ny,nz),  dqdt(nx,ny,nz)
!  DIMENSION DUDT(NX,NY,NZ),  DVDT(NX,NY,NZ)
  DIMENSION dqldt(nx,ny,nz), dqrdt(nx,ny,nz)
  DIMENSION dqidt(nx,ny,nz), dqsdt(nx,ny,nz)
  DIMENSION raincv(nx,ny)
  DIMENSION nca(nx,ny)
  DIMENSION icuyes(nx),lsb(nx)

  REAL :: kffbfct, dt
  INTEGER :: mphyopt
!
! To allow sub-saturation
!
  INTEGER :: kfsubsattrig 

  REAL :: rhlcl 
  REAL :: dqsdt_rh
  REAL :: dtrh 
!
!...DEFINE LOCAL VARIABLES...
!
  DIMENSION p0(kmx),z0(kmx),t0(kmx),tv0(kmx),q0(kmx),                   &
      u0(kmx),v0(kmx),tu(kmx),tvu(kmx),qu(kmx),tz(kmx),                 &
      tvd(kmx),qd(kmx),qes(kmx),thtes(kmx),tg(kmx),tvg(kmx),            &
      qg(kmx),wu(kmx),wd(kmx),w0(kmx),ems(kmx),emsd(kmx),               &
      umf(kmx),uer(kmx),udr(kmx),dmf(kmx),der(kmx),ddr(kmx),            &
      dzq(kmx),umf2(kmx),uer2(kmx),udr2(kmx),dmf2(kmx),der2(kmx),       &
      ddr2(kmx),dza(kmx),thta0(kmx),thetee(kmx),thtau(kmx),             &
      theteu(kmx),thtad(kmx),theted(kmx),qliq(kmx),qice(kmx),           &
      qlqout(kmx),qicout(kmx),pptliq(kmx),pptice(kmx),detlq(kmx),       &
      detic(kmx),detlq2(kmx),detic2(kmx),ratio(kmx),ratio2(kmx)
  DIMENSION domgdp(kmx),exn(kmx),rhoe(kmx),tvqu(kmx),dp(kmx),           &
      rh(kmx),eqfrc(kmx),wspd(kmx),qdt(kmx),fxm(kmx),thtag(kmx),        &
      thtesg(kmx),thpa(kmx),thfxin(kmx),thfxout(kmx),qpa(kmx),          &
      qfxin(kmx),qfxout(kmx),qlpa(kmx),qlfxin(kmx),qlfxout(kmx),        &
      qipa(kmx),qifxin(kmx),qifxout(kmx),qrpa(kmx),qrfxin(kmx),         &
      qrfxout(kmx),qspa(kmx),qsfxin(kmx),qsfxout(kmx),ql0(kmx),         &
      qlg(kmx),qi0(kmx),qig(kmx),qr0(kmx),qrg(kmx),qs0(kmx),qsg(kmx)

  DIMENSION omg(kmx+1)
  DIMENSION rainfb(kmx),snowfb(kmx)
!
  DATA p00,t00/1.e5,273.16/
  DATA cv,b61,rlf/717.,0.608,3.339E5/
  DATA rhic,rhbc/1.,0.90/
  DATA pie,ttfrz,tbfrz,c5/3.141592654,268.16,248.16,1.0723E-3/
  DATA rate,rad/0.01,1500./
  DATA r,g,cp/287.04,9.81,1005.7/
  DATA xlv0,xlv1,xls0,xls1/3.147E6,2369.0,2.905E6,259.532/
  DATA aliq,bliq,cliq,dliq/613.3,17.502,4780.8,32.19/
  DATA aice,bice,cice,dice/613.2,22.452,6133.0,0.61/

  ix=nx-1
  kx=nz-3
  kxp1=kx+1
  kl=kx
  klm=kx-1

  gdry = -g/cp
  rovg = r/g
  dxsq = dx*dx
  dt2 = 2*dt

  xtime = 10.0
!
!*****************************************************************************
!...ADDITIONAL PARAMETERS THAT WILL NEED TO BE DEFINED SOMEWHERE...
!...IN MM4, THEY ARE DEFINED IN subroutine param AND STORED IN
!...VARIOUS COMMON BLOCKS
!
!   DX=25000.             ! HORIZONTAL GRID LENGTH (m)
!   DXSQ=DX*DX            ! GRID AREA (M**2)
!   DT2=80.               ! 2*DT (s)
!
!...DEFINE CONSTANTS FOR CALCULATION OF LATENT HEATING...
!
!    XLV0=3.147E6         ! LV = XLV0 + TMP*XLV1
!    XLV1=2369.
!    XLS0=2.905E6         ! LS = XLS0 + TMP*XLS1
!    XLS1=259.532
!
!...DEFINE CONSTANTS FOR CALCULATION OF SATURATION VAPOR PRESSURE
!...ACCORDING TO BUCK (J. APPL. METEO., DECEMBER, 1981)...
!
!    ALIQ = 613.3   !ES(LIQ)=ALIQ*EXP((BLIQ*TMP-CLIQ)/(TMP-DLIQ))
!    BLIQ = 17.502
!    CLIQ = 4780.8
!    DLIQ = 32.19
!    AICE = 613.2   !ES(ICE)=AICE*EXP((BICE*TMP-CICE)/(TMP-DICE))
!    BICE = 22.452
!    CICE = 6133.0
!    DICE = 0.61
!
!****************************************************************************
                                                    ! PPT FB MODS
!...OPTION TO FEED CONVECTIVELY GENERATED RAINWATER    ! PPT FB MODS
!...INTO GRID-RESOLVED RAINWATER (OR SNOW/GRAUPEL)     ! PPT FB MODS
!...FIELD.  'FBFRC' IS THE FRACTION OF AVAILABLE       ! PPT FB MODS
!...PRECIPITATION TO BE FED BACK (0.0 - 1.0)...        ! PPT FB MODS

  fbfrc=kffbfct         !from ARPS input file arps.input
!
!...SCHEME IS CALLED ONCE  ON EACH NORTH-SOUTH SLICE, THE LOOP BELOW
!...CHECKS FOR THE POSSIBILITY OF INITIATING PARAMETERIZED
!...CONVECTION AT EACH POINT WITHIN THE SLICE
!
!   DO 325 I=1,IX
  DO nc = 1,ncuyes
    i = icuyes(nc)
!
!...SEE IF IT IS NECESSARY TO CHECK FOR CONVECTIVE TRIGGERING AT THIS
!...GRID POINT...IF NCA>0, CONVECTION IS ALREADY ACTIVE AT THIS POINT,
!...JUST FEED BACK THE TENDENCIES SAVED FROM THE TIME WHEN CONVECTION
!...WAS INITIATED.  IF NCA<0, CONVECTION IS NOT ACTIVE
!...AND YOU MAY WANT TO CHECK TO SEE IF IT CAN BE ACTIVATED FOR THE
!...CURRENT CONDITIONS.  IN PREVIOUS APLICATIONS OF THIS SCHEME,
!...THE VARIABLE ICLDCK WAS USED BELOW TO SAVE TIME BY ONLY CHECKING
!...FOR THE POSSIBILITY OF CONVECTIVE INITIATION EVERY 5 OR 10
!...MINUTES...
!
    IF(nca(i,j) >= 1) CYCLE

!     IF (ICLDCK.NE.0) GO TO 325
!
    p300=1000.*(psb(i,j)*a(i,j,kl)+ptop(i,j)-30.)
!
!...INPUT A VERTICAL SOUNDING ...
!
!...*********** NOTE THAT MODEL LAYERS ARE NUMBERED ***************
!...*********** FROM THE BOTTOM-UP IN THE KF SCHEME ***************
!
    ml=0
    cell=ptop(i,j)/psb(i,j)
    DO k=1,kx
      nk=kx-k+1
      p0(k)=1.e3*(a(i,j,nk)*psb(i,j)+ptop(i,j))
      t0(k)=tb(i,j,nk)
      q0(k)=qvb(i,j,nk)
      q0(k)=AMAX1(q0(k),1.0E-10)
!
!...IF Q0 IS ABOVE SATURATION VALUE, REDUCE IT TO SATURATION LEVEL...
!
      es=aliq*EXP((bliq*t0(k)-cliq)/(t0(k)-dliq))
      qes(k)=0.622*es/(p0(k)-es)
      q0(k)=AMIN1(qes(k),q0(k))
      ql0(k)=0.
      qi0(k)=0.
      qr0(k)=0.
      qs0(k)=0.
!      u0(k)=.25*(ub(i,j,nk)+ub(i+1,j,nk)+ub(i,j+1,nk)+ub(i+1,j+1,nk))
!      v0(k)=.25*(vb(i,j,nk)+vb(i+1,j,nk)+vb(i,j+1,nk)+vb(i+1,j+1,nk))
      u0(k) = ub(i,j,nk)
      v0(k) = vb(i,j,nk)
      tv0(k)=t0(k)*(1.+b61*q0(k))
      rhoe(k)=p0(k)/(r*tv0(k))
!        W0(K) = -101.9368*SCR9(I,NK)/RHOE(K)
      dp(k)=(sigma(i,j,nk+1)-sigma(i,j,nk))*psb(i,j)*1.e3
      dzq(k)=rovg*tv0(k)*ALOG((sigma(i,j,nk+1)+cell)/                   &
             (sigma(i,j,nk)+cell))
!
!...DZQ IS DZ BETWEEN SIGMA SURFACES, DZA IS DZ BETWEEN MODEL HALF LEVEL
!...DP IS THE PRESSURE INTERVAL BETWEEN FULL SIGMA LEVELS...
!...
!
      IF(p0(k) >= 500E2)l5=k
      IF(p0(k) >= 400E2)l4=k
      IF(p0(k) >= p300)llfc=k
      IF(t0(k) > t00)ml=k
    END DO
    z0(1)=.5*dzq(1)
    DO k=2,kl
      z0(k)=z0(k-1)+.5*(dzq(k)+dzq(k-1))
      dza(k-1)=z0(k)-z0(k-1)
    END DO
    dza(kl)=0.
!     KMIX=1
    kmix=lsb(i)
    25   low=kmix
    IF(low > llfc)CYCLE
    lc=low
    mxlayr=0
!
!...ASSUME THAT IN ORDER TO SUPPORT A DEEP UPDRAFT YOU NEED A LAYER OF
!...UNSTABLE AIR 50 TO 100 mb DEEP...TO APPROXIMATE THIS, ISOLATE A
!...GROUP OF ADJACENT INDIVIDUAL MODEL LAYERS, WITH THE BASE AT LEVEL
!...LC, SUCH THAT THE COMBINED DEPTH OF THESE LAYERS IS AT LEAST 60 mb..
!
    nlayrs=0
    dpthmx=0.
    DO nk=lc,kx
      dpthmx=dpthmx+dp(nk)
      nlayrs=nlayrs+1
      IF(dpthmx > 6.e3)GO TO 64
    END DO
    CYCLE
    64   kpbl=lc+nlayrs-1
!     KMIX=LC+1
!
!...go ahead and determine what level to start with for the
!...next mixture in case the current mixture, with base at
!...level LC, is not buoyant...
!...instead of checking mixtures using every single layer,
!...move up in increments of at least 25 mb...
!...!!! make that 20 mb !!!!!! 9/15/97...
!...!!! make that 15 mb !!!!!! 10/97...
!      KMIX=LC+1
!     PM25 = P0(LC)-25.E2
!     PM20 = P0(LC)-20.E2
    pm15 = p0(lc)-15.e2
    DO nk = lc+1,kl
!        IF(P0(NK).LT.PM25)THEN
!        IF(P0(NK).LT.PM20)THEN
      IF(p0(nk) < pm15)THEN
        kmix = nk
        GO TO 67
      END IF
    END DO
    CYCLE
    67      CONTINUE
    thmix=0.
    qmix=0.
    zmix=0.
    pmix=0.
    dpthmx=0.
!
!...FIND THE THERMODYNAMIC CHARACTERISTICS OF THE LAYER BY
!...MASS-WEIGHTING THE CHARACTERISTICS OF THE INDIVIDUAL MODEL
!...LAYERS...
!
    DO nk=lc,kpbl
      dpthmx=dpthmx+dp(nk)
      rocpq=0.2854*(1.-0.28*q0(nk))
      thmix=thmix+dp(nk)*t0(nk)*(p00/p0(nk))**rocpq
      qmix=qmix+dp(nk)*q0(nk)
      zmix=zmix+dp(nk)*z0(nk)
      pmix=pmix+dp(nk)*p0(nk)
    END DO
    thmix=thmix/dpthmx
    qmix=qmix/dpthmx
    zmix=zmix/dpthmx
    pmix=pmix/dpthmx
    rocpq=0.2854*(1.-0.28*qmix)
    tmix=thmix*(pmix/p00)**rocpq
    emix=qmix*pmix/(0.622+qmix)
!
!...FIND THE TEMPERATURE OF THE MIXTURE AT ITS LCL, PRESSURE
!...LEVEL OF LCL...
!
    tlog=ALOG(emix/aliq)
    tdpt=(cliq-dliq*tlog)/(bliq-tlog)
    tlcl=tdpt-(.212+1.571E-3*(tdpt-t00)-4.36E-4*(tmix-t00))*(tmix-      &
         tdpt)
    tlcl=AMIN1(tlcl,tmix)
    tvlcl=tlcl*(1.+0.608*qmix)
    cporq=1./rocpq
    plcl=p00*(tlcl/thmix)**cporq
    DO nk=lc,kl
      klcl=nk
!      IF(PLCL.GE.P0(NK))GOTO 35
! test
      IF(plcl >= p0(nk) .AND. klcl >= 2) GO TO 35
    END DO
    CYCLE
    35   k=klcl-1
    dlp=ALOG(plcl/p0(k))/ALOG(p0(klcl)/p0(k))
!
!...ESTIMATE ENVIRONMENTAL TEMPERATURE AND MIXING RATIO AT THE LCL...
!
    tenv=t0(k)+(t0(klcl)-t0(k))*dlp
    qenv=q0(k)+(q0(klcl)-q0(k))*dlp
    tven=tenv*(1.+0.608*qenv)
    tvbar=0.5*(tv0(k)+tven)
!     ZLCL=Z0(K)+R*TVBAR*ALOG(P0(K)/PLCL)/G
    zlcl=z0(k)+(z0(klcl)-z0(k))*dlp
!
!...CHECK TO SEE IF CLOUD IS BUOYANT USING FRITSCH-CHAPPELL TRIGGER
!...FUNCTION DESCRIBED IN KAIN AND FRITSCH (1992)...W0AVG IS AN
!...APROXIMATE VALUE FOR THE RUNNING-MEAN GRID-SCALE VERTICAL
!...VELOCITY, WHICH GIVES SMOOTHER FIELDS OF CONVECTIVE INITIATION
!...THAN THE INSTANTANEOUS VALUE...FORMULA RELATING TEMPERATURE
!...PERTURBATION TO VERTICAL VELOCITY HAS BEEN USED WITH THE MOST
!...SUCCESS AT GRID LENGTHS NEAR 25 km.  FOR DIFFERENT GRID-LENGTHS,
!...ADJUST VERTICAL VELOCITY TO EQUIVALENT VALUE FOR 25 KM GRID
!...LENGTH, ASSUMING LINEAR DEPENDENCE OF W ON GRID LENGTH...
!
!    WKLCL=0.02*ZLCL/2.5E3
!    WKLCL=0.02*ZLCL/2.0E3
    wklcl=0.02*zlcl/2.0E3
    wkl=(w0avg(i,j,k)+(w0avg(i,j,klcl)-w0avg(i,j,k))*dlp)*dx/25.e3-     &
        wklcl
    wabs=ABS(wkl)+1.e-10
    wsigne=wkl/wabs
    dtlcl=4.64*wsigne*wabs**0.33
!
! TO ALLOW SUB-SATURATION IN THE CUMULUS PARAMETERIZATION SCHEME
! 
! Original Author: Jack Kain
! Slightly Modified By: Zuwen He, 04/2002
!
! Original code from an email of Jack on 04/18/2002, which is his 
! version for the ETA model to allow subsaturation. 
! 
! Quote: 
!
! An additional term to the trigger function to enhance the 
! likelihood of convective initiation when the environment 
! at the LCL is near saturation.   --- Jack 
!
! End of quote
!
! rhsat is the threshold value for condensation to begin.  This code 
! is designed to give a maximum value of DTRH at a relative humidity 
! of 95%, decreasing rapidly to zero for higher values of RH.  DTRH 
! is then added to DTLCL to determine whether a parcel might be 
! buoyant at the LCL. 
!
! BEGIN OF CODE OF ALLOWING SUB-SATURATION
!
    IF (kfsubsattrig == 1) THEN 
!      rhsat=0.75 
!      IF(rhsat < 1.)THEN 
        rhlcl = qenv/(qes(k)+(qes(klcl)-qes(k))*dlp)
        dqsdt_rh = qmix*(cliq-bliq*dliq)/((tlcl-dliq)*(tlcl-dliq)) 
        if(rhlcl >= 0.75 .AND. rhlcl <= 0.95)THEN 
          dtrh = 0.25*(rhlcl-0.75)*qmix/dqsdt_rh
        ELSEIF(rhlcl > 0.95)THEN 
          dtrh = (1./rhlcl-1.)*qmix/dqsdt_rh
        ELSE 
          dtrh = 0. 
        ENDIF 
!      ENDIF 
      dtlcl=dtlcl+dtrh 
    END IF       ! kfsubsattrig==1
!
! END OF CODE OF ALLOWING SUB-SATURATION 
!
    gdt=g*dtlcl*(zlcl-z0(lc))/(tv0(lc)+tven)
    wlcl=1.+.5*wsigne*SQRT(ABS(gdt)+1.e-10)
    IF(tlcl+dtlcl > tenv)GO TO 45
    IF(kpbl >= llfc)CYCLE
    GO TO 25
!
!...CONVECTIVE TRIGGERING CRITERIA HAS BEEN SATISFIED...COMPUTE
!...EQUIVALENT POTENTIAL TEMPERATURE
!...(THETEU) AND VERTICAL VELOCITY OF THE RISING PARCEL AT THE LCL...
!
    45   CONTINUE
    theteu(k)=tmix*(1.e5/pmix)**(0.2854*(1.-0.28*qmix))*                &
              EXP((3374.6525/tlcl-2.5403)*qmix*(1.+0.81*qmix))
    es=aliq*EXP((tenv*bliq-cliq)/(tenv-dliq))
    tvavg=0.5*(tv0(klcl)+tenv*(1.+0.608*qenv))
    plcl=p0(klcl)*EXP(g/(r*tvavg)*(z0(klcl)-zlcl))
    qese=0.622*es/(plcl-es)
    gdt=g*dtlcl*(zlcl-z0(lc))/(tv0(lc)+tven)
    wlcl=1.+.5*wsigne*SQRT(ABS(gdt)+1.e-10)
    thtes(k)=tenv*(1.e5/plcl)**(0.2854*(1.-0.28*qese))*                 &
             EXP((3374.6525/tenv-2.5403)*qese*(1.+0.81*qese))
    wtw=wlcl*wlcl
    IF(wlcl < 0.)GO TO 25
    tvlcl=tlcl*(1.+0.608*qmix)
    rholcl=plcl/(r*tvlcl)
!
    lcl=klcl
    let=lcl
!
!*******************************************************************
!                                                               *
!              COMPUTE UPDRAFT PROPERTIES                       *
!                                                               *
!*******************************************************************
!
!
!...ESTIMATE INITIAL UPDRAFT MASS FLUX (UMF(K))...
!
    wu(k)=wlcl
    au0=pie*rad*rad
    umf(k)=rholcl*au0
    vmflcl=umf(k)
    upold=vmflcl
    upnew=upold
!
!...RATIO2 IS THE DEGREE OF GLACIATION IN THE CLOUD (0 TO 1),
!...UER IS THE ENVIR ENTRAINMENT RATE, ABE IS AVAILABLE
!...BUOYANT ENERGY, TRPPT IS THE TOTAL RATE OF PRECIPITATION
!...PRODUCTION...
!
    ratio2(k)=0.
    uer(k)=0.
    abe=0.
    trppt=0.
    tu(k)=tlcl
    tvu(k)=tvlcl
    qu(k)=qmix
    eqfrc(k)=1.
    qliq(k)=0.
    qice(k)=0.
    qlqout(k)=0.
    qicout(k)=0.
    detlq(k)=0.
    detic(k)=0.
    pptliq(k)=0.
    pptice(k)=0.
    iflag=0
    kfrz=lc
!
!...THE AMOUNT OF CONV AVAIL POT ENERGY (CAPE) IS CALCULATED WITH
!...RESPECT TO UNDILUTE PARCEL ASCENT; EQ POT TEMP OF UNDILUTE
!...PARCEL IS THTUDL, UNDILUTE TEMPERATURE IS GIVEN BY TUDL...
!
    thtudl=theteu(k)
    tudl=tlcl
!
!...TTEMP IS USED DURING CALCULATION OF THE LINEAR GLACIATION
!...PROCESS; IT IS INITIALLY SET TO THE TEMPERATURE AT WHICH
!...FREEZING IS SPECIFIED TO BEGIN.  WITHIN THE GLACIATION
!...INTERVAL, IT IS SET EQUAL TO THE UPDRAFT TEMP AT THE
!...PREVIOUS MODEL LEVEL...
!
    ttemp=ttfrz
!
!...ENTER THE LOOP FOR UPDRAFT CALCULATIONS...CALCULATE UPDRAFT TEMP,
!...MIXING RATIO, VERTICAL MASS FLUX, LATERAL DETRAINMENT OF MASS AND
!...MOISTURE, PRECIPITATION RATES AT EACH MODEL LEVEL...
!
    DO nk=k,klm
      nk1=nk+1
      ratio2(nk1)=ratio2(nk)
!
!...UPDATE UPDRAFT PROPERTIES AT THE NEXT MODEL LVL TO REFLECT
!...ENTRAINMENT OF ENVIRONMENTAL AIR...
!
      frc1=0.
      tu(nk1)=t0(nk1)
      theteu(nk1)=theteu(nk)
      qu(nk1)=qu(nk)
      qliq(nk1)=qliq(nk)
      qice(nk1)=qice(nk)
      CALL tpmix(p0(nk1),theteu(nk1),tu(nk1),qu(nk1),qliq(nk1),         &
           qice(nk1),qnewlq,qnewic,ratio2(nk1),rl,xlv0,xlv1,            &
           xls0,xls1,aliq,bliq,cliq,dliq,aice,bice,cice,dice)
      tvu(nk1)=tu(nk1)*(1.+0.608*qu(nk1))
!
!...CHECK TO SEE IF UPDRAFT TEMP IS WITHIN THE FREEZING INTERVAL,
!...IF IT IS, CALCULATE THE FRACTIONAL CONVERSION TO GLACIATION
!...AND ADJUST QNEWLQ TO REFLECT THE GRADUAL CHANGE IN THETAU
!...SINCE THE LAST MODEL LEVEL...THE GLACIATION EFFECTS WILL BE
!...DETERMINED AFTER THE AMOUNT OF CONDENSATE AVAILABLE AFTER
!...PRECIP FALLOUT IS DETERMINED...TTFRZ IS THE TEMP AT WHICH
!...GLACIATION BEGINS, TBFRZ THE TEMP AT WHICH IT ENDS...
!
      IF(tu(nk1) <= ttfrz.AND.iflag < 1)THEN
        IF(tu(nk1) > tbfrz)THEN
          IF(ttemp > ttfrz)ttemp=ttfrz
          frc1=(ttemp-tu(nk1))/(ttfrz-tbfrz)
          r1=(ttemp-tu(nk1))/(ttemp-tbfrz)
        ELSE
          frc1=(ttemp-tbfrz)/(ttfrz-tbfrz)
          r1=1.
          iflag=1
        END IF
        qnwfrz=qnewlq
        qnewic=qnewic+qnewlq*r1*0.5
        qnewlq=qnewlq-qnewlq*r1*0.5
        effq=(ttfrz-tbfrz)/(ttemp-tbfrz)
        ttemp=tu(nk1)
      END IF
!
!  CALCULATE UPDRAFT VERTICAL VELOCITY AND PRECIPITATION FALLOUT...
!
      IF(nk == k)THEN
        be=(tvlcl+tvu(nk1))/(tven+tv0(nk1))-1.
        boterm=2.*(z0(nk1)-zlcl)*g*be/1.5
        enterm=0.
        dzz=z0(nk1)-zlcl
      ELSE
        be=(tvu(nk)+tvu(nk1))/(tv0(nk)+tv0(nk1))-1.
        boterm=2.*dza(nk)*g*be/1.5
        enterm=2.*uer(nk)*wtw/upold
        dzz=dza(nk)
      END IF
      wsq=wtw
      CALL condload_old(qliq(nk1),qice(nk1),wtw,dzz,boterm,enterm,rate,     &
           qnewlq,qnewic,qlqout(nk1),qicout(nk1))

!
!...IF VERT VELOCITY IS LESS THAN ZERO, EXIT THE UPDRAFT LOOP AND,
!...IF CLOUD IS TALL ENOUGH, FINALIZE UPDRAFT CALCULATIONS...
!
     IF(wtw <= 0.) GOTO 65   !cfb

      wabs=SQRT(ABS(wtw))
      wu(nk1)=wtw/wabs

!cfb     IF(wu(nk1) < 0.)EXIT
!
!  UPDATE THE ABE FOR UNDILUTE ASCENT...
!
      thtes(nk1)=t0(nk1)*                                               &
                 (1.e5/p0(nk1))**(0.2854*(1.-0.28*qes(nk1)))*           &
          EXP((3374.6525/t0(nk1)-2.5403)*qes(nk1)*(1.+0.81*qes(nk1)))
      udlbe=((2.*thtudl)/(thtes(nk)+thtes(nk1))-1.)*dzz
      IF(udlbe > 0.)abe=abe+udlbe*g
!
!  DETERMINE THE EFFECTS OF CLOUD GLACIATION IF WITHIN THE SPECIFIED
!  TEMP INTERVAL...
!
      IF(frc1 > 1.e-6)THEN
        CALL dtfrznew_old(tu(nk1),p0(nk1),theteu(nk1),qu(nk1),qliq(nk1),    &
             qice(nk1),ratio2(nk1),ttfrz,tbfrz,qnwfrz,rl,frc1,effq,     &
             iflag,xlv0,xlv1,xls0,xls1,aliq,bliq,cliq,dliq,aice,bice    &
             ,cice,dice)
      END IF

!
!  CALL SUBROUTINE TO CALCULATE ENVIRONMENTAL EQUIVALENT POTENTIAL TEMP.
!  WITHIN GLACIATION INTERVAL, THETAE MUST BE CALCULATED WITH RESPECT TO
!  SAME DEGREE OF GLACIATION FOR ALL ENTRAINING AIR...
!
      CALL envirtht_old(p0(nk1),t0(nk1),q0(nk1),thetee(nk1),ratio2(nk1),    &
           rl,aliq,bliq,cliq,dliq,aice,bice,cice,dice)

!
!...REI IS THE RATE OF ENVIRONMENTAL INFLOW...
!
      rei=vmflcl*dp(nk1)*0.03/rad
      tvqu(nk1)=tu(nk1)*(1.+0.608*qu(nk1)-qliq(nk1)-qice(nk1))
!
!...IF CLOUD PARCELS ARE VIRTUALLY COLDER THAN THE ENVIRONMENT, NO
!   ENTRAINMENT IS ALLOWED AT THIS LEVEL...
!
      IF(tvqu(nk1) <= tv0(nk1))THEN
        uer(nk1)=0.0
        udr(nk1)=rei
        ee2=0.
        ud2=1.
        eqfrc(nk1)=0.
        GO TO 55
      END IF
      let=nk1
      ttmp=tvqu(nk1)
!
!...DETERMINE THE CRITICAL MIXED FRACTION OF UPDRAFT AND ENVIRONMENTAL
!...FOR ESTIMATION OF ENTRAINMENT AND DETRAINMENT RATES...
!
      f1=0.95
      f2=1.-f1
      thttmp=f1*thetee(nk1)+f2*theteu(nk1)
      qtmp=f1*q0(nk1)+f2*qu(nk1)
      tmpliq=f2*qliq(nk1)
      tmpice=f2*qice(nk1)
      CALL tpmix(p0(nk1),thttmp,ttmp,qtmp,tmpliq,tmpice,qnewlq,         &
           qnewic,ratio2(nk1),rl,xlv0,xlv1,xls0,xls1,aliq,bliq,cliq,    &
           dliq,aice,bice,cice,dice)
      tu95=ttmp*(1.+0.608*qtmp-tmpliq-tmpice)
      IF(tu95 > tv0(nk1))THEN
        ee2=1.
        ud2=0.
        eqfrc(nk1)=1.0
        GO TO 50
      END IF
      f1=0.10
      f2=1.-f1
      thttmp=f1*thetee(nk1)+f2*theteu(nk1)
      qtmp=f1*q0(nk1)+f2*qu(nk1)
      tmpliq=f2*qliq(nk1)
      tmpice=f2*qice(nk1)
      CALL tpmix(p0(nk1),thttmp,ttmp,qtmp,tmpliq,tmpice,qnewlq,         &
           qnewic,ratio2(nk1),rl,xlv0,xlv1,xls0,xls1,aliq,bliq,cliq,    &
           dliq,aice,bice,cice,dice)
      tu10=ttmp*(1.+0.608*qtmp-tmpliq-tmpice)
      IF(tu10 == tvqu(nk1))THEN
        ee2=1.
        ud2=0.
        eqfrc(nk1)=1.0
        GO TO 50
      END IF
      eqfrc(nk1)=(tv0(nk1)-tvqu(nk1))*f1/(tu10-tvqu(nk1))
      eqfrc(nk1)=AMAX1(0.0,eqfrc(nk1))
      eqfrc(nk1)=AMIN1(1.0,eqfrc(nk1))
      IF(eqfrc(nk1) == 1)THEN
        ee2=1.
        ud2=0.
        GO TO 50
      ELSE IF(eqfrc(nk1) == 0.)THEN
        ee2=0.
        ud2=1.
        GO TO 50
      ELSE
!
!...SUBROUTINE PROF5 INTEGRATES OVER THE GAUSSIAN DIST TO DETERMINE THE
!   FRACTIONAL ENTRAINMENT AND DETRAINMENT RATES...
!
        CALL prof5_old(eqfrc(nk1),ee2,ud2)
      END IF
!
      50     IF(nk == k)THEN
        ee1=1.
        ud1=0.
      END IF
!
!...NET ENTRAINMENT AND DETRAINMENT RATES ARE GIVEN BY THE AVERAGE FRACT
!   VALUES IN THE LAYER...
!
      uer(nk1)=0.5*rei*(ee1+ee2)
      udr(nk1)=0.5*rei*(ud1+ud2)
!
!...IF THE CALCULATED UPDRAFT DETRAINMENT RATE IS GREATER THAN THE TOTAL
!   UPDRAFT MASS FLUX, ALL CLOUD MASS DETRAINS, EXIT UPDRAFT CALCULATION
!
      55     IF(umf(nk)-udr(nk1) < 10.)THEN
!
!...IF THE CALCULATED DETRAINED MASS FLUX IS GREATER THAN THE TOTAL UPD
!   FLUX, IMPOSE TOTAL DETRAINMENT OF UPDRAFT MASS AT THE PREVIOUS MODEL
!
        IF(udlbe > 0.)abe=abe-udlbe*g
        let=nk
!      WRITE(98,1015)P0(NK1)/100.
        EXIT
      END IF
      ee1=ee2
      ud1=ud2
      upold=umf(nk)-udr(nk1)
      upnew=upold+uer(nk1)
      umf(nk1)=upnew
!
!...DETLQ AND DETIC ARE THE RATES OF DETRAINMENT OF LIQUID AND ICE IN TH
!   DETRAINING UPDRAFT MASS...
!
      detlq(nk1)=qliq(nk1)*udr(nk1)
      detic(nk1)=qice(nk1)*udr(nk1)
      qdt(nk1)=qu(nk1)
      qu(nk1)=(upold*qu(nk1)+uer(nk1)*q0(nk1))/upnew
      theteu(nk1)=(theteu(nk1)*upold+thetee(nk1)*uer(nk1))/upnew
      qliq(nk1)=qliq(nk1)*upold/upnew
      qice(nk1)=qice(nk1)*upold/upnew
!
!...KFRZ IS THE HIGHEST MODEL LEVEL AT WHICH LIQUID CONDENSATE IS GENERA
!   PPTLIQ IS THE RATE OF GENERATION (FALLOUT) OF LIQUID PRECIP AT A GIV
!   MODEL LVL, PPTICE THE SAME FOR ICE, TRPPT IS THE TOTAL RATE OF PRODU
!   OF PRECIP UP TO THE CURRENT MODEL LEVEL...
!
      IF(ABS(ratio2(nk1)-1.) > 1.e-6)kfrz=nk1
      pptliq(nk1)=qlqout(nk1)*(umf(nk)-udr(nk1))
      pptice(nk1)=qicout(nk1)*(umf(nk)-udr(nk1))
      trppt=trppt+pptliq(nk1)+pptice(nk1)
      IF(nk1 <= kpbl)uer(nk1)=uer(nk1)+vmflcl*dp(nk1)/dpthmx
    END DO
!
!...CHECK CLOUD DEPTH...IF CLOUD IS TALL ENOUGH, ESTIMATE THE EQUILIBRIU
!   TEMPERATURE LEVEL (LET) AND ADJUST MASS FLUX PROFILE AT CLOUD TOP SO
!   THAT MASS FLUX DECREASES TO ZERO AS A LINEAR FUNCTION OF PRESSURE BE
!   THE LET AND CLOUD TOP...
!
!...LTOP IS THE MODEL LEVEL JUST BELOW THE LEVEL AT WHICH VERTICAL VELOC
!   FIRST BECOMES NEGATIVE...
!

65    ltop=nk
    cldhgt=z0(ltop)-zlcl
!
!...IF CLOUD TOP HGT IS LESS THAN SPECIFIED MINIMUM HEIGHT, GO BACK AND
!   THE NEXT HIGHEST 60MB LAYER TO SEE IF A BIGGER CLOUD CAN BE OBTAINED
!   THAT SOURCE AIR...
!
!     IF(CLDHGT.LT.4.E3.OR.ABE.LT.1.)THEN            ! JSK MODS
    IF(cldhgt < 3.e3.OR.abe < 1.)THEN             ! JSK MODS
      DO nk=k,ltop
        umf(nk)=0.
        udr(nk)=0.
        uer(nk)=0.
        detlq(nk)=0.
        detic(nk)=0.
        pptliq(nk)=0.
        pptice(nk)=0.
      END DO
      GO TO 25
    END IF
!
!...IF THE LET AND LTOP ARE THE SAME, DETRAIN ALL OF THE UPDRAFT MASS FL
!   THIS LEVEL...
!
    IF(let == ltop)THEN
      udr(ltop)=umf(ltop)+udr(ltop)-uer(ltop)
      detlq(ltop)=qliq(ltop)*udr(ltop)*upnew/upold
      detic(ltop)=qice(ltop)*udr(ltop)*upnew/upold
      trppt=trppt-(pptliq(ltop)+pptice(ltop))
      uer(ltop)=0.
      umf(ltop)=0.
      GO TO 85
    END IF
!
!   BEGIN TOTAL DETRAINMENT AT THE LEVEL ABOVE THE LET...
!
    dptt=0.
    DO nj=let+1,ltop
      dptt=dptt+dp(nj)
    END DO
    dumfdp=umf(let)/dptt
!
!...ADJUST MASS FLUX PROFILES, DETRAINMENT RATES, AND PRECIPITATION FALL
!   RATES TO REFLECT THE LINEAR DECREASE IN MASS FLX BETWEEN THE LET AND
!
    DO nk=let+1,ltop
      udr(nk)=dp(nk)*dumfdp
      umf(nk)=umf(nk-1)-udr(nk)
      detlq(nk)=qliq(nk)*udr(nk)
      detic(nk)=qice(nk)*udr(nk)
      trppt=trppt-pptliq(nk)-pptice(nk)
      pptliq(nk)=(umf(nk-1)-udr(nk))*qlqout(nk)
      pptice(nk)=(umf(nk-1)-udr(nk))*qicout(nk)
      trppt=trppt+pptliq(nk)+pptice(nk)
    END DO
!
!...SEND UPDRAFT CHARACTERISTICS TO OUTPUT FILES...
!
    85   CONTINUE
!  WRITE(98,3001)
!    3001  FORMAT(' ')
!  WRITE(98,3002)XTIME,INEST,I,J
!    3002  FORMAT('XTIME =',f9.2,', INEST =',i2,', GRID POINT: I =',i2,  &
!           ' J =',i2)
!   IF(MXLAYR.EQ.1)THEN
!    PRINT *,'HFX, QFX, PBL, THGV =',HFX(I,J),QFX(I,J),PBL(I,J),THGV
!  PRINT *,'RHOE(1), WPTHP, WPQP, WPTHVP, WSTR =',RHOE(1),WPTHP,WPQP,
!  *WPTHVP,WSTR
!   ENDIF
!  PRINT *,'PRS AT BASE OF UPDRAFT SOURCE LAYER =',P0(LC)/100.
!   PRINT *,'DEPTH OF UPDRAFT SOURCE MIXED LAYER =',DPTHMX/100.
!   PRINT 1000,'  P  ','  VMFU  ','  TU  ',' EQFRC',' CLDWAT ',
!  *' CLDICE ',' PPTLIQ ',' PPTICE ',' DETLQ ',' DETIC '
!   DO 1020 NK=KLCL-1,LTOP
!      IF(NK.EQ.KLCL-1)THEN
!        PRS=PLCL/100.
!      ELSE
!        PRS=P0(NK)/100.
!      ENDIF
!C         PRS=CVMGT(PLCL/100.,P0(NK)/100.,NK.EQ.KLCL-1)
!      PRINT 1005,PRS,UMF(NK)/VMFLCL,TU(NK)-273.16,EQFRC(NK),
!  *   QLIQ(NK)*1.E3,QICE(NK)*1.E3,QLQOUT(NK)*1.E3,QICOUT(NK)*1.E3,
!  *   DETLQ(NK)/1.E3,DETIC(NK)/1.E3
!1020  CONTINUE
!   PRINT 1010,P0(NK1)/100.
!
!...EXTEND THE UPDRAFT MASS FLUX PROFILE DOWN TO THE SOURCE LAYER FOR TH
!   UPDRAFT AIR...ALSO, DEFINE THETAE FOR LEVELS BELOW THE LCL...
!
    DO nk=1,k
      IF(nk >= lc)THEN
        IF(nk == lc)THEN
          umf(nk)=vmflcl*dp(nk)/dpthmx
          uer(nk)=vmflcl*dp(nk)/dpthmx
        ELSE IF(nk <= kpbl)THEN
          uer(nk)=vmflcl*dp(nk)/dpthmx
          umf(nk)=umf(nk-1)+uer(nk)
        ELSE
          umf(nk)=vmflcl
          uer(nk)=0.
        END IF
        tu(nk)=tmix+(z0(nk)-zmix)*gdry
        qu(nk)=qmix
        wu(nk)=wlcl
      ELSE
        tu(nk)=0.
        qu(nk)=0.
        umf(nk)=0.
        wu(nk)=0.
        uer(nk)=0.
      END IF
      udr(nk)=0.
      qdt(nk)=0.
      qliq(nk)=0.
      qice(nk)=0.
      qlqout(nk)=0.
      qicout(nk)=0.
      pptliq(nk)=0.
      pptice(nk)=0.
      detlq(nk)=0.
      detic(nk)=0.
      ratio2(nk)=0.
      ee=q0(nk)*p0(nk)/(0.622+q0(nk))
      tlog=ALOG(ee/aliq)
      tdpt=(cliq-dliq*tlog)/(bliq-tlog)
      tsat=tdpt-(.212+1.571E-3*(tdpt-t00)-4.36E-4*(t0(nk)-t00))*(       &
           t0(nk)-tdpt)
      thta=t0(nk)*(1.e5/p0(nk))**(0.2854*(1.-0.28*q0(nk)))
      thetee(nk)=thta*                                                  &
                 EXP((3374.6525/tsat-2.5403)*q0(nk)*(1.+0.81*q0(nk))    &
                 )
      thtes(nk)=thta*                                                   &
                EXP((3374.6525/t0(nk)-2.5403)*qes(nk)*(1.+0.81*         &
                qes(nk)))
      eqfrc(nk)=1.0
    END DO
!
    ltop1=ltop+1
    ltopm1=ltop-1
!
!...DEFINE VARIABLES ABOVE CLOUD TOP...
!
    DO nk=ltop1,kx
      umf(nk)=0.
      udr(nk)=0.
      uer(nk)=0.
      qdt(nk)=0.
      qliq(nk)=0.
      qice(nk)=0.
      qlqout(nk)=0.
      qicout(nk)=0.
      detlq(nk)=0.
      detic(nk)=0.
      pptliq(nk)=0.
      pptice(nk)=0.
      IF(nk > ltop1)THEN
        tu(nk)=0.
        qu(nk)=0.
        wu(nk)=0.
      END IF
      thta0(nk)=0.
      thtau(nk)=0.
      ems(nk)=dp(nk)*dxsq/g
      emsd(nk)=1./ems(nk)
      tg(nk)=t0(nk)
      qg(nk)=q0(nk)
      qlg(nk)=0.
      qig(nk)=0.
      qrg(nk)=0.
      qsg(nk)=0.
      omg(nk)=0.
    END DO
    omg(kxp1)=0.
    p150=p0(klcl)-1.50E4
    DO nk=1,ltop
      ems(nk)=dp(nk)*dxsq/g
      emsd(nk)=1./ems(nk)
!
!...INITIALIZE SOME VARIABLES TO BE USED LATER IN THE VERT ADVECTION SCH
!
      exn(nk)=(p00/p0(nk))**(0.2854*(1.-0.28*qdt(nk)))
      thtau(nk)=tu(nk)*exn(nk)
      exn(nk)=(p00/p0(nk))**(0.2854*(1.-0.28*q0(nk)))
      thta0(nk)=t0(nk)*exn(nk)
!
!...LVF IS THE LEVEL AT WHICH MOISTURE FLUX IS ESTIMATED AS THE BASIS FO
!...PRECIPITATION EFFICIENCY CALCULATIONS...
!
      IF(p0(nk) > p150)lvf=nk
      omg(nk)=0.
    END DO
    lvf=MIN0(lvf,let)                                  ! JSK MODS
    usr=umf(lvf+1)*(qu(lvf+1)+qliq(lvf+1)+qice(lvf+1))
    usr=AMIN1(usr,trppt)
!
!  WRITE(98,1025)KLCL,ZLCL,DTLCL,LTOP,P0(LTOP),IFLAG,
!    * TMIX-T00,PMIX,QMIX,ABE
!  WRITE(98,1030)P0(LET)/100.,P0(LTOP)/100.,VMFLCL,PLCL/100.,
!    * WLCL,CLDHGT
!
!...COMPUTE CONVECTIVE TIME SCALE(TIMEC). THE MEAN WIND AT THE LCL
!...AND MIDTROPOSPHERE IS USED.
!
    wspd(klcl)=SQRT(u0(klcl)*u0(klcl)+v0(klcl)*v0(klcl))
    wspd(l5)=SQRT(u0(l5)*u0(l5)+v0(l5)*v0(l5))
    wspd(ltop)=SQRT(u0(ltop)*u0(ltop)+v0(ltop)*v0(ltop))
    vconv=.5*(wspd(klcl)+wspd(l5))
    timec=dx/vconv
    tadvec=timec
    timec=AMAX1(1800.,timec)
    timec=AMIN1(3600.,timec)
    nic=nint(timec/(.5*dt2))
    timec=REAL(nic)*.5*dt2
!
!...COMPUTE WIND SHEAR AND PRECIPITATION EFFICIENCY.
!
!     SHSIGN = CVMGT(1.,-1.,WSPD(LTOP).GT.WSPD(KLCL))
    IF(wspd(ltop) > wspd(klcl))THEN
      shsign=1.
    ELSE
      shsign=-1.
    END IF
    vws=(u0(ltop)-u0(klcl))*(u0(ltop)-u0(klcl))+(v0(ltop)-v0(klcl))*    &
        (v0(ltop)-v0(klcl))
    vws=1.e3*shsign*SQRT(vws)/(z0(ltop)-z0(lcl))
    pef=1.591+vws*(-.639+vws*(9.53E-2-vws*4.96E-3))
    pef=AMAX1(pef,.2)
    pef=AMIN1(pef,.9)
!
!...PRECIPITATION EFFICIENCY IS A FUNCTION OF THE HEIGHT OF CLOUD BASE.
!
    cbh=(zlcl-z0(1))*3.281E-3
    IF(cbh < 3.)THEN
      rcbh=.02
    ELSE
      rcbh=.96729352+cbh*(-.70034167+cbh*(.162179896+cbh*(-             &
           1.2569798E-2+cbh*(4.2772E-4-cbh*5.44E-6))))
    END IF
    IF(cbh > 25)rcbh=2.4
    pefcbh=1./(1.+rcbh)
    pefcbh=AMIN1(pefcbh,.9)
!
!... MEAN PEF. IS USED TO COMPUTE RAINFALL.
!
    peff=.5*(pef+pefcbh)
    peff2 = peff                                ! JSK MODS
!
!     WRITE(98,1035)PEF,PEFCBH,LC,LET,WKL,VWS
!
!*****************************************************************
!                                                             *
!               COMPUTE DOWNDRAFT PROPERTIES                  *
!                                                             *
!*****************************************************************
!
!...LET DOWNDRAFT ORIGINATE AT THE LEVEL OF MINIMUM SATURATION EQUIVALEN
!...POTENTIAL TEMPERATURE (SEQT) IN THE CLOUD LAYER, EXTEND DOWNWARD TO
!...SURFACE, OR TO THE LAYER BELOW CLOUD BASE AT WHICH ENVIR SEQT IS LES
!...THAN MIN SEQT IN THE CLOUD LAYER...LET DOWNDRAFT DETRAIN OVER A LAYE
!...OF SPECIFIED PRESSURE-DEPTH (DPDD)...
!
    tder=0.
    kstart=MAX0(kpbl,klcl)
    thtmin=thtes(kstart+1)
    kmin=kstart+1
    DO nk=kstart+2,ltop-1
      thtmin=AMIN1(thtmin,thtes(nk))
      IF(thtmin == thtes(nk))kmin=nk
    END DO
    lfs=kmin
    IF(ratio2(lfs) > 0.)CALL envirtht_old(p0(lfs),t0(lfs),q0(lfs),          &
        thetee(lfs),0.,rl,aliq,bliq,cliq,dliq,aice,bice,cice,dice)
    eqfrc(lfs)=(thtes(lfs)-theteu(lfs))/(thetee(lfs)-theteu(lfs))
    eqfrc(lfs)=AMAX1(eqfrc(lfs),0.)
    eqfrc(lfs)=AMIN1(eqfrc(lfs),1.)
    theted(lfs)=thtes(lfs)
!
!...ESTIMATE THE EFFECT OF MELTING PRECIPITATION IN THE DOWNDRAFT...
!
    IF(ml > 0)THEN
      dtmltd=0.5*(qu(klcl)-qu(ltop))*rlf/cp
    ELSE
      dtmltd=0.
    END IF
    tz(lfs)=t0(lfs)-dtmltd
    es=aliq*EXP((tz(lfs)*bliq-cliq)/(tz(lfs)-dliq))
    qs=0.622*es/(p0(lfs)-es)
    qd(lfs)=eqfrc(lfs)*q0(lfs)+(1.-eqfrc(lfs))*qu(lfs)
    thtad(lfs)=tz(lfs)*(p00/p0(lfs))**(0.2854*(1.-0.28*qd(lfs)))
    IF(qd(lfs) >= qs)THEN
      theted(lfs)=thtad(lfs)*                                           &
                  EXP((3374.6525/tz(lfs)-2.5403)*qs*(1.+0.81*qs))
    ELSE
      CALL envirtht_old(p0(lfs),tz(lfs),qd(lfs),theted(lfs),0.,rl,aliq,     &
           bliq,cliq,dliq,aice,bice,cice,dice)
    END IF
    DO nk=1,lfs
      nd=lfs-nk
      IF(theted(lfs) > thtes(nd).OR.nd == 1)THEN
        ldb=nd
                                                          ! JSK MODS
!...IF DOWNDRAFT NEVER BECOMES NEGATIVELY BUOYANT OR IF IT   ! JSK MODS
!...IS SHALLOWER 50 mb, DON'T ALLOW IT TO OCCUR AT ALL...    ! JSK MODS
                                                          ! JSK MODS
        IF(nk == 1.OR.(p0(ldb)-p0(lfs)) < 50.e2)GO TO 141  ! JSK MODS
        EXIT
      END IF
    END DO
!
!...ALLOW DOWNDRAFT TO DETRAIN IN A SINGLE LAYER, BUT WITH DOWNDRAFT AIR
!...TYPICALLY FLUSHED UP INTO HIGHER LAYERS AS ALLOWED IN THE TOTAL
!...VERTICAL ADVECTION CALCULATIONS FARTHER DOWN IN THE CODE...
!
    dpdd=dp(ldb)
    ldt=ldb
    frc=1.
    dpt=0.
!   DO 115 NK=LDB,LFS
!     DPT=DPT+DP(NK)
!     IF(DPT.GT.DPDD)THEN
!       LDT=NK
!       FRC=(DPDD+DP(NK)-DPT)/DP(NK)
!       GOTO 120
!     ENDIF
!     IF(NK.EQ.LFS-1)THEN
!      LDT=NK
!     FRC=1.
!     DPDD=DPT
!     GOTO 120
!     ENDIF
!115   CONTINUE
!    120   CONTINUE
!
!...TAKE A FIRST GUESS AT THE INITIAL DOWNDRAFT MASS FLUX...
!
    tvd(lfs)=t0(lfs)*(1.+0.608*qes(lfs))
    rdd=p0(lfs)/(r*tvd(lfs))
    a1=(1.-peff)*au0
    dmf(lfs)=-a1*rdd
    der(lfs)=eqfrc(lfs)*dmf(lfs)
    ddr(lfs)=0.
    DO nd=lfs-1,ldb,-1
      nd1=nd+1
      IF(nd <= ldt)THEN
        der(nd)=0.
        ddr(nd)=-dmf(ldt+1)*dp(nd)*frc/dpdd
        dmf(nd)=dmf(nd1)+ddr(nd)
        frc=1.
        theted(nd)=theted(nd1)
        qd(nd)=qd(nd1)
      ELSE
        der(nd)=dmf(lfs)*0.03*dp(nd)/rad
        ddr(nd)=0.
        dmf(nd)=dmf(nd1)+der(nd)
        IF(ratio2(nd) > 0.)CALL envirtht_old(p0(nd),t0(nd),q0(nd),          &
            thetee(nd),0.,rl,aliq,bliq,cliq,dliq,aice,bice,cice,dice)
        theted(nd)=(theted(nd1)*dmf(nd1)+thetee(nd)*der(nd))/dmf(nd)
        qd(nd)=(qd(nd1)*dmf(nd1)+q0(nd)*der(nd))/dmf(nd)
      END IF
    END DO
    tder=0.
!
!...CALCULATION AN EVAPORATION RATE FOR GIVEN MASS FLUX...
!
    DO nd=ldb,ldt
      tz(nd)=                                                           &
             tpdd(p0(nd),theted(ldt),t0(nd),qs,qd(nd),1.0,xlv0,xlv1,    &
             aliq,bliq,cliq,dliq,aice,bice,cice,dice)
      es=aliq*EXP((tz(nd)*bliq-cliq)/(tz(nd)-dliq))
      qs=0.622*es/(p0(nd)-es)
      dssdt=(cliq-bliq*dliq)/((tz(nd)-dliq)*(tz(nd)-dliq))
      rl=xlv0-xlv1*tz(nd)
      dtmp=rl*qs*(1.-rhbc)/(cp+rl*rhbc*qs*dssdt)
      t1rh=tz(nd)+dtmp
      es=rhbc*aliq*EXP((bliq*t1rh-cliq)/(t1rh-dliq))
      qsrh=0.622*es/(p0(nd)-es)
!
!...CHECK TO SEE IF MIXING RATIO AT SPECIFIED RH IS LESS THAN ACTUAL
!...MIXING RATIO...IF SO, ADJUST TO GIVE ZERO EVAPORATION...
!
      IF(qsrh < qd(nd))THEN
        qsrh=qd(nd)
!       T1RH=T1+(QS-QSRH)*RL/CP
        t1rh=tz(nd)
      END IF
      tz(nd)=t1rh
      qs=qsrh
      tder=tder+(qs-qd(nd))*ddr(nd)
      qd(nd)=qs
      thtad(nd)=tz(nd)*(p00/p0(nd))**(0.2854*(1.-0.28*qd(nd)))
    END DO
!
!...IF DOWNDRAFT DOES NOT EVAPORATE ANY WATER FOR SPECIFIED RELATIVE
!...HUMIDITY, NO DOWNDRAFT IS ALLOWED...
!
    141   IF(tder < 1.)THEN
!       WRITE(98,3004)I,J
!      3004       FORMAT(' ','I=',i3,2X,'J=',i3)
      pptflx=trppt
      cpr=trppt
      tder=0.
      cndtnf=0.
      updinc=1.
      ldb=lfs
      DO ndk=1,ltop
        dmf(ndk)=0.
        der(ndk)=0.
        ddr(ndk)=0.
        thtad(ndk)=0.
        wd(ndk)=0.
        tz(ndk)=0.
        qd(ndk)=0.
      END DO
      aincm2=100.
      GO TO 165
    END IF
!
!...ADJUST DOWNDRAFT MASS FLUX SO THAT EVAPORATION RATE IN DOWNDRAFT IS
!...CONSISTENT WITH PRECIPITATION EFFICIENCY RELATIONSHIP...
!
    devdmf=tder/dmf(lfs)
    ppr=0.
    pptflx=peff*usr
    rced=trppt-pptflx
!
!...PPR IS THE TOTAL AMOUNT OF PRECIPITATION THAT FALLS  OUT OF THE UPDR
!...FROM CLOUD BASE TO THE LFS...UPDRAFT MASS FLUX WILL BE INCREASED UP
!...THE LFS TO ACCOUNT FOR UPDRAFT AIR MIXING WITH ENVIRONMENTAL AIR TO
!...THE UPDRAFT, SO PPR WILL INCREASE PROPORTIONATELY...
!
    DO nm=klcl,lfs
      ppr=ppr+pptliq(nm)+pptice(nm)
    END DO
    IF(lfs >= klcl)THEN
      dpptdf=(1.-peff)*ppr*(1.-eqfrc(lfs))/umf(lfs)
    ELSE
      dpptdf=0.
    END IF
!
!...CNDTNF IS THE AMOUNT OF CONDENSATE TRANSFERRED ALONG WITH UPDRAFT MA
!...THE DOWNDRAFT AT THE LFS...
!
    cndtnf=(qliq(lfs)+qice(lfs))*(1.-eqfrc(lfs))
    dmflfs=rced/(devdmf+dpptdf+cndtnf)
    IF(dmflfs > 0.)THEN
      tder=0.
      GO TO 141
    END IF
!
!...DDINC IS THE FACTOR BY WHICH TO INCREASE THE FIRST-GUESS DOWNDRAFT M
!...FLUX TO SATISFY THE PRECIP EFFICIENCY RELATIONSHIP, UPDINC IS THE FA
!...WHICH TO INCREASE THE UPDRAFT MASS FLUX BELOW THE LFS TO ACCOUNT FOR
!...TRANSFER OF MASS FROM UPDRAFT TO DOWNDRAFT...
!
!     DDINC=DMFLFS/DMF(LFS)                              ! JSK MODS

    IF(lfs >= klcl)THEN
      updinc=(umf(lfs)-(1.-eqfrc(lfs))*dmflfs)/umf(lfs)
                                                         ! JSK MODS
!...LIMIT UPDINC TO LESS THAN OR EQUAL TO 1.5...            ! JSK MODS
                                                         ! JSK MODS
      IF(updinc > 1.5)THEN                             ! JSK MODS
        updinc = 1.5                                    ! JSK MODS
        dmflfs2 = umf(lfs)*(updinc-1.)/(eqfrc(lfs)-1.)  ! JSK MODS
        rced2 = dmflfs2*(devdmf+dpptdf+cndtnf)          ! JSK MODS
        pptflx = pptflx+(rced-rced2)                    ! JSK MODS
        peff2 = pptflx/usr                              ! JSK MODS
        rced = rced2                                    ! JSK MODS
        dmflfs = dmflfs2                                ! JSK MODS
      END IF                                             ! JSK MODS
    ELSE
      updinc=1.
    END IF
    ddinc=dmflfs/dmf(lfs)                     ! JSK MODS
    DO nk=ldb,lfs
      dmf(nk)=dmf(nk)*ddinc
      der(nk)=der(nk)*ddinc
      ddr(nk)=ddr(nk)*ddinc
    END DO
    cpr=trppt+ppr*(updinc-1.)
    pptflx=pptflx+peff*ppr*(updinc-1.)
    peff=peff2                                ! JSK MODS

    tder=tder*ddinc
!
!...ADJUST UPDRAFT MASS FLUX, MASS DETRAINMENT RATE, AND LIQUID WATER AN
!   DETRAINMENT RATES TO BE CONSISTENT WITH THE TRANSFER OF THE ESTIMATE
!   FROM THE UPDRAFT TO THE DOWNDRAFT AT THE LFS...
!
    DO nk=lc,lfs
      umf(nk)=umf(nk)*updinc
      udr(nk)=udr(nk)*updinc
      uer(nk)=uer(nk)*updinc
      pptliq(nk)=pptliq(nk)*updinc
      pptice(nk)=pptice(nk)*updinc
      detlq(nk)=detlq(nk)*updinc
      detic(nk)=detic(nk)*updinc
    END DO
!
!...ZERO OUT THE ARRAYS FOR DOWNDRAFT DATA AT LEVELS ABOVE AND BELOW THE
!...DOWNDRAFT...
!
    IF(ldb > 1)THEN
      DO nk=1,ldb-1
        dmf(nk)=0.
        der(nk)=0.
        ddr(nk)=0.
        wd(nk)=0.
        tz(nk)=0.
        qd(nk)=0.
        thtad(nk)=0.
      END DO
    END IF
    DO nk=lfs+1,kx
      dmf(nk)=0.
      der(nk)=0.
      ddr(nk)=0.
      wd(nk)=0.
      tz(nk)=0.
      qd(nk)=0.
      thtad(nk)=0.
    END DO
    DO nk=ldt+1,lfs-1
      tz(nk)=0.
      qd(nk)=0.
      thtad(nk)=0.                                   ! JSK MODS
    END DO
!
!...SET LIMITS ON THE UPDRAFT AND DOWNDRAFT MASS FLUXES SO THAT THE INFL
!   INTO CONVECTIVE DRAFTS FROM A GIVEN LAYER IS NO MORE THAN IS AVAILAB
!   IN THAT LAYER INITIALLY...
!
    165   aincmx=1000.
    lmax=MAX0(klcl,lfs)
    DO nk=lc,lmax
      IF((uer(nk)-der(nk)) > 0.)aincm1=ems(nk)/((uer(nk)-der(nk))* timec)
      aincmx=AMIN1(aincmx,aincm1)
    END DO
    ainc=1.
    IF(aincmx < ainc)ainc=aincmx
!
!...SAVE THE RELEVENT VARIABLES FOR A UNIT UPDRFT AND DOWNDRFT...THEY WI
!...ITERATIVELY ADJUSTED BY THE FACTOR AINC TO SATISFY THE STABILIZATION
!...CLOSURE...
!
    ncount=0
    tder2=tder
    pptfl2=pptflx
    DO nk=1,ltop
      detlq2(nk)=detlq(nk)
      detic2(nk)=detic(nk)
      udr2(nk)=udr(nk)
      uer2(nk)=uer(nk)
      ddr2(nk)=ddr(nk)
      der2(nk)=der(nk)
      umf2(nk)=umf(nk)
      dmf2(nk)=dmf(nk)
    END DO
    fabe=1.
    stab=0.95
!   XNIN=AMIN0(I,J,ILX-I,JLX-J)
!   IF(XNIN.LT.5)STAB=STAB*(XNIN-1.)*0.25
    noitr=0
    IF(ainc/aincmx > 0.999)THEN
      ncount=0
      GO TO 255
    END IF
    istop=0
    175   ncount=ncount+1
!
!*****************************************************************
!                                                             *
!        COMPUTE PROPERTIES FOR COMPENSATIONAL SUBSIDENCE     *
!                                                             *
!*****************************************************************
!
!...DETERMINE OMEGA VALUE NECESSARY AT TOP AND BOTTOM OF EACH LAYER TO
!...SATISFY MASS CONTINUITY...
!
!    185   CONTINUE
    dtt=timec
    DO nk=1,ltop
      domgdp(nk)=-(uer(nk)-der(nk)-udr(nk)-ddr(nk))*emsd(nk)
      IF(nk > 1)THEN
        omg(nk)=omg(nk-1)-dp(nk-1)*domgdp(nk-1)
        dtt1=0.75*dp(nk-1)/(ABS(omg(nk))+1.e-10)
        dtt=AMIN1(dtt,dtt1)
      END IF
    END DO
    DO nk=1,ltop
      thpa(nk)=thta0(nk)
      qpa(nk)=q0(nk)
      nstep=nint(timec/dtt+1)
      dtime=timec/FLOAT(nstep)
      fxm(nk)=omg(nk)*dxsq/g
    END DO

!
!...DO AN UPSTREAM/FORWARD-IN-TIME ADVECTION OF THETA, QV...
!
    DO ntc=1,nstep
!
!...ASSIGN THETA AND Q VALUES AT THE TOP AND BOTTOM OF EACH LAYER BASED
!...SIGN OF OMEGA...
!
      DO nk=1,ltop
        thfxin(nk)=0.
        thfxout(nk)=0.
        qfxin(nk)=0.
        qfxout(nk)=0.
      END DO
      DO nk=2,ltop
        IF(omg(nk) <= 0.)THEN
          thfxin(nk)=-fxm(nk)*thpa(nk-1)
          qfxin(nk)=-fxm(nk)*qpa(nk-1)
          thfxout(nk-1)=thfxout(nk-1)+thfxin(nk)
          qfxout(nk-1)=qfxout(nk-1)+qfxin(nk)
        ELSE
          thfxout(nk)=fxm(nk)*thpa(nk)
          qfxout(nk)=fxm(nk)*qpa(nk)
          thfxin(nk-1)=thfxin(nk-1)+thfxout(nk)
          qfxin(nk-1)=qfxin(nk-1)+qfxout(nk)
        END IF
      END DO
!
!...UPDATE THE THETA AND QV VALUES AT EACH LEVEL...
!
      DO nk=1,ltop
        thpa(nk)=thpa(nk)+(thfxin(nk)+udr(nk)*thtau(nk)+ddr(nk)*        &
                 thtad(nk)-thfxout(nk)-(uer(nk)-der(nk))*thta0(nk))*    &
                 dtime*emsd(nk)
        qpa(nk)=qpa(nk)+(qfxin(nk)+udr(nk)*qdt(nk)+ddr(nk)*qd(nk)-      &
                qfxout(nk)-(uer(nk)-der(nk))*q0(nk))*dtime*emsd(nk)
      END DO
    END DO
    DO nk=1,ltop
      thtag(nk)=thpa(nk)
      qg(nk)=qpa(nk)
    END DO
!
!...CHECK TO SEE IF MIXING RATIO DIPS BELOW ZERO ANYWHERE;  IF SO, BORRO
!...MOISTURE FROM ADJACENT LAYERS TO BRING IT BACK UP ABOVE ZERO...
!
    DO nk=1,ltop
      IF(qg(nk) < 0.)THEN
        IF(nk == 1)THEN                             ! JSK MODS
          PRINT *,'!!!!! PROBLEM WITH KF SCHEME:  ' ! JSK MODS
          PRINT *,'QG = 0 AT THE SURFACE!!!!!!!'    ! JSK MODS
          CALL arpsstop('QG',1)                     ! JSK MODS
        END IF                                       ! JSK MODS
        nk1=nk+1
        IF(nk == ltop)nk1=klcl
        tma=qg(nk1)*ems(nk1)
        tmb=qg(nk-1)*ems(nk-1)
        tmm=(qg(nk)-1.e-9)*ems(nk  )
        bcoeff=-tmm/((tma*tma)/tmb+tmb)
        acoeff=bcoeff*tma/tmb
        tmb=tmb*(1.-bcoeff)
        tma=tma*(1.-acoeff)
        IF(nk == ltop)THEN
          qvdiff=(qg(nk1)-tma*emsd(nk1))*100./qg(nk1)
          IF(ABS(qvdiff) > 1.)THEN
          PRINT *,'& !!!WARNING!!! CLOUD BASE WATER VAPOR CHANGES BY ', &
             qvdiff,'% WHEN MOISTURE IS BORROWED TO PREVENT NEGATIVE ', &
             'VALUES IN KAIN-FRITSCH'
          END IF
        END IF
        qg(nk)=1.e-9
        qg(nk1)=tma*emsd(nk1)
        qg(nk-1)=tmb*emsd(nk-1)
      END IF
    END DO
    topomg=(udr(ltop)-uer(ltop))*dp(ltop)*emsd(ltop)
    IF(ABS(topomg-omg(ltop)) > 1.e-3)THEN
!    WRITE(98,*)'ERROR:  MASS DOES NOT BALANCE IN KF SCHEME;
!    * TOPOMG, OMG =',TOPOMG,OMG(LTOP)
      WRITE(6,*)                                                        & 
      'ERROR:  MASS DOES NOT BALANCE IN KF SCHEME; topomg, omg =',      &
      TOPOMG,OMG(LTOP)
      istop=1
      GO TO 265
    END IF
!
!...CONVERT THETA TO T...
!
    DO nk=1,ltop
      exn(nk)=(p00/p0(nk))**(0.2854*(1.-0.28*qg(nk)))
      tg(nk)=thtag(nk)/exn(nk)
      tvg(nk)=tg(nk)*(1.+0.608*qg(nk))
    END DO
!
!*******************************************************************
!                                                               *
!  COMPUTE NEW CLOUD AND CHANGE IN AVAILABLE BUOYANT ENERGY.    *
!                                                               *
!*******************************************************************
!
!...THE FOLLOWING COMPUTATIONS ARE SIMILAR TO THAT FOR UPDRAFT
!
    thmix=0.
    qmix=0.
    pmix=0.
    DO nk=lc,kpbl
      rocpq=0.2854*(1.-0.28*qg(nk))
      thmix=thmix+dp(nk)*tg(nk)*(p00/p0(nk))**rocpq
      qmix=qmix+dp(nk)*qg(nk)
      pmix=pmix+dp(nk)*p0(nk)
    END DO
    thmix=thmix/dpthmx
    qmix=qmix/dpthmx
    pmix=pmix/dpthmx
    rocpq=0.2854*(1.-0.28*qmix)
    tmix=thmix*(pmix/p00)**rocpq
    es=aliq*EXP((tmix*bliq-cliq)/(tmix-dliq))
    qs=0.622*es/(pmix-es)
!
!...REMOVE SUPERSATURATION FOR DIAGNOSTIC PURPOSES, IF NECESSARY...
!
    IF(qmix > qs)THEN
      rl=xlv0-xlv1*tmix
      cpm=cp*(1.+0.887*qmix)
      dssdt=qs*(cliq-bliq*dliq)/((tmix-dliq)*(tmix-dliq))
      dq=(qmix-qs)/(1.+rl*dssdt/cpm)
      tmix=tmix+rl/cp*dq
      qmix=qmix-dq
      rocpq=0.2854*(1.-0.28*qmix)
      thmix=tmix*(p00/pmix)**rocpq
      tlcl=tmix
      plcl=pmix
    ELSE
      qmix=AMAX1(qmix,0.)
      emix=qmix*pmix/(0.622+qmix)
      tlog=ALOG(emix/aliq)
      tdpt=(cliq-dliq*tlog)/(bliq-tlog)
      tlcl=tdpt-(.212+1.571E-3*(tdpt-t00)-4.36E-4*(tmix-t00))*(tmix-    &
           tdpt)
      tlcl=AMIN1(tlcl,tmix)
      cporq=1./rocpq
      plcl=p00*(tlcl/thmix)**cporq
    END IF
    tvlcl=tlcl*(1.+0.608*qmix)
    DO nk=lc,kl
      klcl=nk
      IF(plcl >= p0(nk))EXIT
    END DO
    k=klcl-1
    dlp=ALOG(plcl/p0(k))/ALOG(p0(klcl)/p0(k))
!
!...ESTIMATE ENVIRONMENTAL TEMPERATURE AND MIXING RATIO AT THE LCL...
!
    tenv=tg(k)+(tg(klcl)-tg(k))*dlp
    qenv=qg(k)+(qg(klcl)-qg(k))*dlp
    tven=tenv*(1.+0.608*qenv)
    tvbar=0.5*(tvg(k)+tven)
!     ZLCL=Z0(K)+R*TVBAR*ALOG(P0(K)/PLCL)/G
    zlcl=z0(k)+(z0(klcl)-z0(k))*dlp
    tvavg=0.5*(tven+tg(klcl)*(1.+0.608*qg(klcl)))
    plcl=p0(klcl)*EXP(g/(r*tvavg)*(z0(klcl)-zlcl))
    theteu(k)=tmix*(1.e5/pmix)**(0.2854*(1.-0.28*qmix))*                &
              EXP((3374.6525/tlcl-2.5403)*qmix*(1.+0.81*qmix))
    es=aliq*EXP((tenv*bliq-cliq)/(tenv-dliq))
    qese=0.622*es/(plcl-es)
    thtesg(k)=tenv*(1.e5/plcl)**(0.2854*(1.-0.28*qese))*                &
              EXP((3374.6525/tenv-2.5403)*qese*(1.+0.81*qese))
!
!...COMPUTE ADJUSTED ABE(ABEG).
!
    abeg=0.
    thtudl=theteu(k)
    DO nk=k,ltopm1
      nk1=nk+1
      es=aliq*EXP((tg(nk1)*bliq-cliq)/(tg(nk1)-dliq))
      qese=0.622*es/(p0(nk1)-es)
      thtesg(nk1)=tg(nk1)*(1.e5/p0(nk1))**(0.2854*(1.-0.28*qese))*      &
                  EXP((3374.6525/tg(nk1)-2.5403)*qese*(1.+0.81*qese)    &
                  )
!      DZZ=CVMGT(Z0(KLCL)-ZLCL,DZA(NK),NK.EQ.K)
      IF(nk == k)THEN
        dzz=z0(klcl)-zlcl
      ELSE
        dzz=dza(nk)
      END IF
      be=((2.*thtudl)/(thtesg(nk1)+thtesg(nk))-1.)*dzz
      IF(be > 0.)abeg=abeg+be*g
    END DO
!
!...ASSUME AT LEAST 90% OF CAPE (ABE) IS REMOVED BY CONVECTION DURING
!...THE PERIOD TIMEC...
!
    IF(noitr == 1)THEN
!     WRITE(98,1060)FABE
      GO TO 265
    END IF
    dabe=AMAX1(abe-abeg,0.1*abe)
    fabe=abeg/(abe+1.e-8)
    IF(fabe > 1.)THEN
!    WRITE(98,*)'UPDRAFT/DOWNDRAFT COUPLET INCREASES CAPE AT THIS
!    *GRID POINT; NO CONVECTION ALLOWED!'
      CYCLE
    END IF
    IF(ncount /= 1)THEN
      dfda=(fabe-fabeold)/(ainc-aincold)
      IF(dfda > 0.)THEN
        noitr=1
        ainc=aincold
        GO TO 255
      END IF
    END IF
    aincold=ainc
    fabeold=fabe
    IF(ainc/aincmx > 0.999.AND.fabe > 1.05-stab)THEN
!   WRITE(98,1055)FABE
      GO TO 265
    END IF
    IF(fabe <= 1.05-stab.AND.fabe >= 0.95-stab)GO TO 265
    IF(ncount > 10)THEN
!    WRITE(98,1060)FABE
      GO TO 265
    END IF
!
!...IF MORE THAN 10% OF THE ORIGINAL CAPE REMAINS, INCREASE THE CONVECTI
!...MASS FLUX BY THE FACTOR AINC:
!
    IF(fabe == 0.)THEN
      ainc=ainc*0.5
    ELSE
      ainc=ainc*stab*abe/(dabe+1.e-8)
    END IF
    255   ainc=AMIN1(aincmx,ainc)
!...IF AINC BECOMES VERY SMALL, EFFECTS OF CONVECTION ! JSK MODS
!...WILL BE MINIMAL SO JUST IGNORE IT...          ! JSK MODS
    IF(ainc < 0.05)CYCLE
!     AINC=AMAX1(AINC,0.05)                    ! JSK MODS
    tder=tder2*ainc
    pptflx=pptfl2*ainc
!  WRITE(98,1080)LFS,LDB,LDT,TIMEC,NSTEP,NCOUNT,FABEOLD,AINCOLD
    DO nk=1,ltop
      umf(nk)=umf2(nk)*ainc
      dmf(nk)=dmf2(nk)*ainc
      detlq(nk)=detlq2(nk)*ainc
      detic(nk)=detic2(nk)*ainc
      udr(nk)=udr2(nk)*ainc
      uer(nk)=uer2(nk)*ainc
      der(nk)=der2(nk)*ainc
      ddr(nk)=ddr2(nk)*ainc
    END DO
!
!...GO BACK UP FOR ANOTHER ITERATION...
!
    GO TO 175
    265   CONTINUE
!
!...CLEAN THINGS UP, CALCULATE CONVECTIVE FEEDBACK TENDENCIES FOR THIS
!...GRID POINT...
!
!...COMPUTE HYDROMETEOR TENDENCIES AS IS DONE FOR T, QV...
!
!...FRC2 IS THE FRACTION OF TOTAL CONDENSATE      !  PPT FB MODS
!...GENERATED THAT GOES INTO PRECIPITIATION       !  PPT FB MODS
    frc2=pptflx/(cpr*ainc)                    !  PPT FB MODS
    DO nk=1,ltop
      qlpa(nk)=ql0(nk)
      qipa(nk)=qi0(nk)
      qrpa(nk)=qr0(nk)
      qspa(nk)=qs0(nk)
      rainfb(nk)=pptliq(nk)*ainc*fbfrc*frc2   !  PPT FB MODS
      snowfb(nk)=pptice(nk)*ainc*fbfrc*frc2   !  PPT FB MODS
    END DO
    DO ntc=1,nstep
!
!...ASSIGN HYDROMETEORS CONCENTRATIONS AT THE TOP AND BOTTOM OF EACH LAY
!...BASED ON THE SIGN OF OMEGA...
!
      DO nk=1,ltop
        qlfxin(nk)=0.
        qlfxout(nk)=0.
        qifxin(nk)=0.
        qifxout(nk)=0.
        qrfxin(nk)=0.
        qrfxout(nk)=0.
        qsfxin(nk)=0.
        qsfxout(nk)=0.
      END DO
      DO nk=2,ltop
        IF(omg(nk) <= 0.)THEN
          qlfxin(nk)=-fxm(nk)*qlpa(nk-1)
          qifxin(nk)=-fxm(nk)*qipa(nk-1)
          qrfxin(nk)=-fxm(nk)*qrpa(nk-1)
          qsfxin(nk)=-fxm(nk)*qspa(nk-1)
          qlfxout(nk-1)=qlfxout(nk-1)+qlfxin(nk)
          qifxout(nk-1)=qifxout(nk-1)+qifxin(nk)
          qrfxout(nk-1)=qrfxout(nk-1)+qrfxin(nk)
          qsfxout(nk-1)=qsfxout(nk-1)+qsfxin(nk)
        ELSE
          qlfxout(nk)=fxm(nk)*qlpa(nk)
          qifxout(nk)=fxm(nk)*qipa(nk)
          qrfxout(nk)=fxm(nk)*qrpa(nk)
          qsfxout(nk)=fxm(nk)*qspa(nk)
          qlfxin(nk-1)=qlfxin(nk-1)+qlfxout(nk)
          qifxin(nk-1)=qifxin(nk-1)+qifxout(nk)
          qrfxin(nk-1)=qrfxin(nk-1)+qrfxout(nk)
          qsfxin(nk-1)=qsfxin(nk-1)+qsfxout(nk)
        END IF
      END DO
!
!...UPDATE THE HYDROMETEOR CONCENTRATION VALUES AT EACH LEVEL...
!
      DO nk=1,ltop
        qlpa(nk)=qlpa(nk)+(qlfxin(nk)+detlq(nk)-qlfxout(nk))*dtime*     &
                 emsd(nk)
        qipa(nk)=qipa(nk)+(qifxin(nk)+detic(nk)-qifxout(nk))*dtime*     &
                 emsd(nk)
        qrpa(nk)=qrpa(nk)+(qrfxin(nk)+qlqout(nk)*udr(nk)-qrfxout(nk)    &
                 +rainfb(nk))*dtime*emsd(nk)         !  PPT FB MODS
!  +               )*DTIME*EMSD(NK)                   !  PPT FB MODS

        qspa(nk)=qspa(nk)+(qsfxin(nk)+qicout(nk)*udr(nk)-qsfxout(nk)    &
                 +snowfb(nk))*dtime*emsd(nk)         !  PPT FB MODS

!  +               )*DTIME*EMSD(NK)                   !  PPT FB MODS

      END DO
    END DO
    DO nk=1,ltop
      qlg(nk)=qlpa(nk)
      qig(nk)=qipa(nk)
      qrg(nk)=qrpa(nk)
      qsg(nk)=qspa(nk)
    END DO

!  WRITE(98,1080)LFS,LDB,LDT,TIMEC,NSTEP,NCOUNT,FABE,AINC
!
!...SEND FINAL PARAMETERIZED VALUES TO OUTPUT FILES...
!
    IF(xtime < 10..OR.istop == 1)THEN
!   WRITE(98,1070)'  P  ','   DP ',' DT K/D ',' DR K/D ','   OMG  ',
!    *' DOMGDP ','   UMF  ','   UER  ','   UDR  ','   DMF  ','   DER  '
!    *,'   DDR  ','   EMS  ','    W0  ','  DETLQ ',' DETIC '
!       DO 300 NK=1,LTOP
!         K=LTOP-NK+1
!         DTT=(TG(K)-T0(K))*86400./TIMEC
!         RL=XLV0-XLV1*TG(K)
!         DR=-(QG(K)-Q0(K))*RL*86400./(TIMEC*CP)
!         UDFRC=UDR(K)*TIMEC*EMSD(K)
!         UEFRC=UER(K)*TIMEC*EMSD(K)
!         DDFRC=DDR(K)*TIMEC*EMSD(K)
!         DEFRC=-DER(K)*TIMEC*EMSD(K)
!   WRITE(98,1075)P0(K)/100.,DP(K)/100.,DTT,DR,OMG(K),DOMGDP(K)*1.E4,
!    *  UMF(K)/1.E6,UEFRC,UDFRC,DMF(K)/1.E6,DEFRC,DDFRC,EMS(K)/1.E11,
!    *  W0AVG(I,J,K)*1.E2,DETLQ(K)*TIMEC*EMSD(K)*1.E3,DETIC(K)*
!    *  TIMEC*EMSD(K)*1.E3
!  300     CONTINUE
!  WRITE(98,1085)'K','P','Z','T0','TG','DT','TU','TD','Q0','QG',
!    *            'DQ','QU','QD','QLG','QIG','QRG','QSG','RH0','RHG'
!       DO 305 NK=1,KX
!         K=KX-NK+1
!         DTT=TG(K)-T0(K)
!         TUC=TU(K)-T00
!         IF(K.LT.LC.OR.K.GT.LTOP)TUC=0.
!         TDC=TZ(K)-T00
!         IF((K.LT.LDB.OR.K.GT.LDT).AND.K.NE.LFS)TDC=0.
!         ES=ALIQ*EXP((BLIQ*TG(K)-CLIQ)/(TG(K)-DLIQ))
!         QGS=ES*0.622/(P0(K)-ES)
!         RH0=Q0(K)/QES(K)
!         RHG=QG(K)/QGS
!  WRITE(98,1090)K,P0(K)/100.,Z0(K),T0(K)-T00,TG(K)-T00,DTT,TUC,
!    *TDC,Q0(K)*1000.,QG(K)*1000.,(QG(K)-Q0(K))*1000.,QU(K)*
!    *1000.,QD(K)*1000.,QLG(K)*1000.,QIG(K)*1000.,QRG(K)*1000.,
!    *QSG(K)*1000.,RH0,RHG
!  305     CONTINUE
!
!...IF CALCULATIONS ABOVE SHOW AN ERROR IN THE MASS BUDGET, PRINT OUT A
!...TO BE USED LATER FOR DIAGNOSTIC PURPOSES, THEN ABORT RUN...
!
      IF(istop == 1)THEN
        DO k = 1,kx
          WRITE(98,1115)z0(k),p0(k)/100.,t0(k)-273.16,q0(k)*1000.,      &
                   u0(k),v0(k),dp(k)/100.,w0avg(i,j,k)
        END DO
        CALL arpsstop('KAIN-FRITSCH',1)
      END IF
    END IF
    cndtnf=(1.-eqfrc(lfs))*(qliq(lfs)+qice(lfs))*dmf(lfs)
!  WRITE(98,1095)CPR*AINC,TDER+PPTFLX+CNDTNF
!
!  EVALUATE MOISTURE BUDGET...
!
    qinit=0.
    qfnl=0.
    dpt=0.
    DO nk=1,ltop
      dpt=dpt+dp(nk)
      qinit=qinit+q0(nk)*ems(nk)
      qfnl=qfnl+qg(nk)*ems(nk)
      qfnl=qfnl+(qlg(nk)+qig(nk)+qrg(nk)+qsg(nk))*ems(nk)
    END DO
    qfnl=qfnl+pptflx*timec*(1.-fbfrc)       !  PPT FB MODS
!     QFNL=QFNL+PPTFLX*TIMEC                 !  PPT FB MODS
    err2=(qfnl-qinit)*100./qinit
!  WRITE(98,1110)QINIT,QFNL,ERR2
!wang   IF(ABS(ERR2).GT.0.05)CALL arpsstop('QVERR',1)
    relerr=err2*qinit/(pptflx*timec+1.e-10)
!  WRITE(98,1120)RELERR
!  WRITE(98,*)'TDER, CPR, USR, TRPPT =',
!    *TDER,CPR*AINC,USR*AINC,TRPPT*AINC
!
!...FEEDBACK TO RESOLVABLE SCALE TENDENCIES.
!
!...IF THE ADVECTIVE TIME PERIOD (TADVEC) IS LESS THAN SPECIFIED MINIMUM
!...TIMEC, ALLOW FEEDBACK TO OCCUR ONLY DURING TADVEC...
!
    IF(tadvec < timec)nic=nint(tadvec/(0.5*dt2))
    nca(i,j)=nic
    DO k=1,kx
      nk=kx-k+1
!      IF(IMOIST(INEST).NE.2)THEN
      IF(mphyopt == 0.0 ) THEN  !no microphy
!
!...IF HYDROMETEORS ARE NOT ALLOWED, THEY MUST BE EVAPORATED OR SUBLIMAT
!...AND FED BACK AS VAPOR, ALONG WITH ASSOCIATED CHANGES IN TEMPERATURE.
!...NOTE:  THIS WILL INTRODUCE CHANGES IN THE CONVECTIVE TEMPERATURE AND
!...WATER VAPOR FEEDBACK TENDENCIES AND MAY LEAD TO SUPERSATURATED VALUE
!...OF QG...
!
        rlc=xlv0-xlv1*tg(k)
        rls=xls0-xls1*tg(k)
        cpm=cp*(1.+0.887*qg(k))
        tg(k)=tg(k)-(rlc*(qlg(k)+qrg(k))+rls*(qig(k)+qsg(k)))/cpm
        qg(k)=qg(k)+(qlg(k)+qrg(k)+qig(k)+qsg(k))
        dqldt(i,j,nk)=0.
        dqidt(i,j,nk)=0.
        dqrdt(i,j,nk)=0.
        dqsdt(i,j,nk)=0.
      ELSE
!     IF(IEXICE.NE.1 .AND. IICE.NE.1) THEN
!        IF(IMPHYS(INEST).EQ.3)THEN
        IF(mphyopt == 1) THEN  !warmrain microphy
!
!...IF ICE PHASE IS NOT ALLOWED, MELT ALL FROZEN HYDROMETEORS...
!
          cpm=cp*(1.+0.887*qg(k))
          tg(k)=tg(k)-(qig(k)+qsg(k))*rlf/cpm
          dqldt(i,j,nk)=(qlg(k)+qig(k)-ql0(k)-qi0(k))/timec
          dqidt(i,j,nk)=0.
          dqrdt(i,j,nk)=(qrg(k)+qsg(k)-qr0(k)-qs0(k))/timec
          dqsdt(i,j,nk)=0.
!     ELSEIF(IEXICE.EQ.1 .AND. IICE.EQ.0)THEN
!         ELSEIF(IMPHYS(INEST).EQ.4)THEN
!
!...IF ICE PHASE IS ALLOWED, BUT MIXED PHASE IS NOT, MELT FROZEN HYDROME
!...BELOW THE MELTING LEVEL, FREEZE LIQUID WATER ABOVE THE MELTING LEVEL
!
!          CPM=CP*(1.+0.887*QG(K))
!          IF(K.LE.ML)THEN
!            TG(K)=TG(K)-(QIG(K)+QSG(K))*RLF/CPM
!          ELSEIF(K.GT.ML)THEN
!            TG(K)=TG(K)+(QLG(K)+QRG(K))*RLF/CPM
!          ENDIF
!          DQLDT(I,J,NK)=(QLG(K)+QIG(K)-QL0(K)-QI0(K))/TIMEC
!          DQIDT(I,J,NK)=0.
!          DQRDT(I,J,NK)=(QRG(K)+QSG(K)-QR0(K)-QS0(K))/TIMEC
!          DQSDT(I,J,NK)=0.
!     ELSEIF(IICE.EQ.1 .AND. IEXICE.EQ.0)THEN
!        ELSEIF(IMPHYS(INEST).GE.5)THEN
        ELSE IF(mphyopt == 2 .OR. mphyopt == 3) THEN !ice
!
!...IF MIXED PHASE HYDROMETEORS ARE ALLOWED, FEED BACK CONVECTIVE TENDEN
!...OF HYDROMETEORS DIRECTLY...
!
          dqldt(i,j,nk)=(qlg(k)-ql0(k))/timec
          dqidt(i,j,nk)=(qig(k)-qi0(k))/timec
          dqrdt(i,j,nk)=(qrg(k)-qr0(k))/timec
          dqsdt(i,j,nk)=(qsg(k)-qs0(k))/timec
!        ELSE
!      PRINT *,'THIS COMBINATION OF IMOIST, IEXICE, IICE NOT ALLOWED!'
!           CALL arpsstop('KAIN-FRITSCH',1)
        END IF
      END IF
      dtdt(i,j,nk)=(tg(k)-t0(k))/timec
      dqdt(i,j,nk)=(qg(k)-q0(k))/timec
    END DO
    raincv(i,j)=.1*.5*dt2*pptflx*(1.-fbfrc)/dxsq     !  PPT FB MODS
!     RAINCV(I,J)=.1*.5*DT2*PPTFLX/DXSQ               !  PPT FB MODS
!      RNC=0.1*TIMEC*PPTFLX/DXSQ
    rnc=raincv(i,j)*nic
!     WRITE(98,909)RNC
!    909     FORMAT(' CONVECTIVE RAINFALL =',f8.4,' CM')

  END DO

!  1000  FORMAT(' ',10A8)
!  1005  FORMAT(' ',f6.0,2X,f6.4,2X,f7.3,1X,f6.4,2X,4(f6.3,2X),2(f7.3,1X))
!  1010  FORMAT(' ',' VERTICAL VELOCITY IS NEGATIVE AT ',f4.0,' MB')
!  1015   FORMAT(' ','ALL REMAINING MASS DETRAINS BELOW ',f4.0,' MB')
!  1025   FORMAT(5X,' KLCL=',i2,' ZLCL=',f7.1,'M',                       &
!          ' DTLCL=',f5.2,' LTOP=',i2,' P0(LTOP)=',-2PF5.1,'MB FRZ LV=', &
!          i2,' TMIX=',0PF4.1,1X,'PMIX=',-2PF6.1,' QMIX=',3PF5.1,        &
!          ' CAPE=',0PF7.1)
!  1030   FORMAT(' ',' P0(LET) = ',f6.1,' P0(LTOP) = ',f6.1,' VMFLCL =', &
!        e12.3,' PLCL =',f6.1,' WLCL =',f6.3,' CLDHGT =',                &
!        f8.1)
!  1035  FORMAT(1X,'PEF(WS)=',f4.2,'(CB)=',f4.2,'LC,LET=',2I3,'WKL='     &
!        ,f6.3,'VWS=',f5.2)
!  1040  FORMAT(' ','PRECIP EFF = 100%, ENVIR CANNOT SUPPORT DOWND afts!')
!  1045  FORMAT('NUMBER OF DOWNDRAFT ITERATIONS EXCEEDS 10...PPTFLX &
!       & is different from that given by precip eff relation!')
!  1050  FORMAT(' ','LCOUNT= ',i3,' PPTFLX/CPR, PEFF= ',f5.3,1X,f5.3,    &
!      & 'DMF(LFS)/UMF(LCL)= ',f5.3)
!  1055     FORMAT(/'*** DEGREE OF STABILIZATION =',f5.3,', NO MORE MASS F &
!      & ux is allowed!')
!  1060     FORMAT(/' ITERATION DOES NOT CONVERGE TO GIVE THE SPECIFIED  &
!      & degree of stabilization!  FABE= ',F6.4)
!  1070 FORMAT (16A8)
!  1075 FORMAT (f8.2,3(f8.2),2(f8.3),f8.2,2F8.3,f8.2,6F8.3)
!  1080   FORMAT(2X,'LFS,LDB,LDT =',3I3,' TIMEC, NSTEP=',f5.0,i3,        &
!      & 'NCOUNT, FABE, AINC=',i2,1X,f5.3,f6.2)
!  1085 FORMAT (a3,16A7,2A8)
!  1090 FORMAT (i3,f7.2,f7.0,10F7.2,4F7.3,2F8.3)
!  1095   FORMAT(' ','  PPT PRODUCTION RATE= ',f10.0,' TOTAL EVAP+PPT= ', &
!      &  f10.0)
!  1105   FORMAT(' ','NET LATENT HEAT RELEASE =',e12.5,' ACTUAL HEATING =', &
!      &  e12.5,' J/KG-S, DIFFERENCE = ',f9.3,'%')
!  1110   FORMAT(' ','INITIAL WATER =',e12.5,' FINAL WATER =',e12.5,     &
!      &  ' TOTAL WATER CHANGE =',f8.2,'%')
  1115 FORMAT(2X,f6.0,2X,f7.2,2X,f5.1,2X,f6.3,2(2X,f5.1),2X,f7.2,2X,f7.4 )
!  1120   FORMAT(' ','MOISTURE ERROR AS FUNCTION OF TOTAL PPT =',f9.3,'%')
  RETURN
END SUBROUTINE kfpara
!**********************************************************************
!   BLOCK DATA
!    COMMON/VAPPRS/ALIQ,BLIQ,CLIQ,DLIQ,AICE,BICE,CICE,DICE,XLS0,XLS1
!    DATA ALIQ,BLIQ,CLIQ,DLIQ/613.3,17.502,4780.8,32.19/
!    DATA AICE,BICE,CICE,DICE/613.2,22.452,6133.0,0.61/
!    DATA XLS0,XLS1/2.905E6,259.532/
!   END
!*********************************************************************

SUBROUTINE tpmix(p,thtu,tu,qu,qliq,qice,qnewlq,qnewic,ratio2,rl,        &
           xlv0,xlv1,xls0,xls1,                                         &
           aliq,bliq,cliq,dliq,aice,bice,cice,dice)
!
!...THIS SUBROUTINE ITERATIVELY EXTRACTS WET-BULB TEMPERATURE FROM EQUIV
!   POTENTIAL TEMPERATURE, THEN CHECKS TO SEE IF SUFFICIENT MOISTURE IS
!   AVAILABLE TO ACHIEVE SATURATION...IF NOT, TEMPERATURE IS ADJUSTED
!   ACCORDINGLY, IF SO, THE RESIDUAL LIQUID WATER/ICE CONCENTRATION IS
!   DETERMINED...
!
  c5=1.0723E-3
  rv=461.51
!
!   ITERATE TO FIND WET BULB TEMPERATURE AS A FUNCTION OF EQUIVALENT POT
!   TEMP AND PRS, ASSUMING SATURATION VAPOR PRESSURE...RATIO2 IS THE DEG
!   OF GLACIATION...
!
  IF(ratio2 < 1.e-6)THEN
    es=aliq*EXP((bliq*tu-cliq)/(tu-dliq))
    qs=0.622*es/(p-es)
    pi=(1.e5/p)**(0.2854*(1.-0.28*qs))
    thtgs=tu*pi*EXP((3374.6525/tu-2.5403)*qs*(1.+0.81*qs))
  ELSE IF(ABS(ratio2-1.) < 1.e-6)THEN
    es=aice*EXP((bice*tu-cice)/(tu-dice))
    qs=0.622*es/(p-es)
    pi=(1.e5/p)**(0.2854*(1.-0.28*qs))
    thtgs=tu*pi*EXP((3114.834/tu-0.278296)*qs*(1.+0.81*qs))
  ELSE
    esliq=aliq*EXP((bliq*tu-cliq)/(tu-dliq))
    esice=aice*EXP((bice*tu-cice)/(tu-dice))
    es=(1.-ratio2)*esliq+ratio2*esice
    qs=0.622*es/(p-es)
    pi=(1.e5/p)**(0.2854*(1.-0.28*qs))
    thtgs=tu*pi*EXP(rl*qs*c5/tu*(1.+0.81*qs))
  END IF
  f0=thtgs-thtu
  t1=tu-0.5*f0
  t0=tu
  itcnt=0
  90 IF(ratio2 < 1.e-6)THEN
    es=aliq*EXP((bliq*t1-cliq)/(t1-dliq))
    qs=0.622*es/(p-es)
    pi=(1.e5/p)**(0.2854*(1.-0.28*qs))
    thtgs=t1*pi*EXP((3374.6525/t1-2.5403)*qs*(1.+0.81*qs))
  ELSE IF(ABS(ratio2-1.) < 1.e-6)THEN
    es=aice*EXP((bice*t1-cice)/(t1-dice))
    qs=0.622*es/(p-es)
    pi=(1.e5/p)**(0.2854*(1.-0.28*qs))
    thtgs=t1*pi*EXP((3114.834/t1-0.278296)*qs*(1.+0.81*qs))
  ELSE
    esliq=aliq*EXP((bliq*t1-cliq)/(t1-dliq))
    esice=aice*EXP((bice*t1-cice)/(t1-dice))
    es=(1.-ratio2)*esliq+ratio2*esice
    qs=0.622*es/(p-es)
    pi=(1.e5/p)**(0.2854*(1.-0.28*qs))
    thtgs=t1*pi*EXP(rl*qs*c5/t1*(1.+0.81*qs))
  END IF
  f1=thtgs-thtu
  IF(ABS(f1) < 0.01)GO TO 50
  itcnt=itcnt+1
  IF(itcnt > 10)GO TO 50
  dt=f1*(t1-t0)/(f1-f0)
  t0=t1
  f0=f1
  t1=t1-dt
  GO TO 90
!
!   IF THE PARCEL IS SUPERSATURATED, CALCULATE CONCENTRATION OF FRESH
!   CONDENSATE...
!
  50 IF(qs <= qu)THEN
    qnew=qu-qs
    qu=qs
    GO TO 96
  END IF
!
!   IF THE PARCEL IS SUBSATURATED, TEMPERATURE AND MIXING RATIO MUST BE
!   ADJUSTED...IF LIQUID WATER OR ICE IS PRESENT, IT IS ALLOWED TO EVAPO
!   SUBLIMATE.
!
  qnew=0.
  dq=qs-qu
  qtot=qliq+qice
!
!   IF THERE IS ENOUGH LIQUID OR ICE TO SATURATE THE PARCEL, TEMP STAYS
!   WET BULB VALUE, VAPOR MIXING RATIO IS AT SATURATED LEVEL, AND THE MI
!   RATIOS OF LIQUID AND ICE ARE ADJUSTED TO MAKE UP THE ORIGINAL SATURA
!   DEFICIT... OTHERWISE, ANY AVAILABLE LIQ OR ICE VAPORIZES AND APPROPR
!   ADJUSTMENTS TO PARCEL TEMP; VAPOR, LIQUID, AND ICE MIXING RATIOS ARE
!
!...NOTE THAT THE LIQ AND ICE MAY BE PRESENT IN PROPORTIONS SLIGHTLY DIF
!   THAN SUGGESTED BY THE VALUE OF RATIO2...CHECK TO MAKE SURE THAT LIQ
!   ICE CONCENTRATIONS ARE NOT REDUCED TO BELOW ZERO WHEN EVAPORATION/
!   SUBLIMATION OCCURS...
!
  IF(qtot >= dq)THEN
    dqice=0.0
    dqliq=0.0
    qliq=qliq-(1.-ratio2)*dq
    IF(qliq < 0.)THEN
      dqice=0.0-qliq
      qliq=0.0
    END IF
    qice=qice-ratio2*dq+dqice
    IF(qice < 0.)THEN
      dqliq=0.0-qice
      qice=0.0
    END IF
    qliq=qliq+dqliq
    qu=qs
    GO TO 96
  ELSE
    IF(ratio2 < 1.e-6)THEN
      rll=xlv0-xlv1*t1
    ELSE IF(ABS(ratio2-1.) < 1.e-6)THEN
      rll=xls0-xls1*t1
    ELSE
      rll=rl
    END IF
    cp=1005.7*(1.+0.89*qu)
    IF(qtot < 1.e-10)THEN
!
!...IF NO LIQUID WATER OR ICE IS AVAILABLE, TEMPERATURE IS GIVEN BY:
      t1=t1+rll*(dq/(1.+dq))/cp
      GO TO 96
    ELSE
!
!...IF SOME LIQ WATER/ICE IS AVAILABLE, BUT NOT ENOUGH TO ACHIEVE SATURA
!   THE TEMPERATURE IS GIVEN BY:
      t1=t1+rll*((dq-qtot)/(1+dq-qtot))/cp
      qu=qu+qtot
      qtot=0.
    END IF
    qliq=0
    qice=0.
  END IF
  96 tu=t1
  qnewlq=(1.-ratio2)*qnew
  qnewic=ratio2*qnew
  IF(itcnt > 10)PRINT*,'***** NUMBER OF ITERATIONS IN TPMIX =', itcnt
  RETURN
END SUBROUTINE tpmix
!************************* TPDD.FOR ************************************
!   THIS SUBROUTINE ITERATIVELY EXTRACTS TEMPERATURE FROM EQUIVALENT   *
!   POTENTIAL TEMP.  IT IS DESIGNED FOR USE WITH DOWNDRAFT CALCULATIONS.
!   IF RELATIVE HUMIDITY IS SPECIFIED TO BE LESS THAN 100%, PARCEL     *
!   TEMP, SPECIFIC HUMIDITY, AND LIQUID WATER CONTENT ARE ITERATIVELY  *
!   CALCULATED.                                                        *
!***********************************************************************

  FUNCTION tpdd(p,thted,tgs,rs,rd,rh,xlv0,xlv1,                         &
      aliq,bliq,cliq,dliq,aice,bice,cice,dice)
  es=aliq*EXP((bliq*tgs-cliq)/(tgs-dliq))
  rs=0.622*es/(p-es)
  pi=(1.e5/p)**(0.2854*(1.-0.28*rs))
  thtgs=tgs*pi*EXP((3374.6525/tgs-2.5403)*rs*(1.+0.81*rs))
  f0=thtgs-thted
  t1=tgs-0.5*f0
  t0=tgs
  cp=1005.7
!
!...ITERATE TO FIND WET-BULB TEMPERATURE...
!
  itcnt=0
  90 es=aliq*EXP((bliq*t1-cliq)/(t1-dliq))
  rs=0.622*es/(p-es)
  pi=(1.e5/p)**(0.2854*(1.-0.28*rs))
  thtgs=t1*pi*EXP((3374.6525/t1-2.5403)*rs*(1.+0.81*rs))
  f1=thtgs-thted
  IF(ABS(f1) < 0.05)GO TO 50
  itcnt=itcnt+1
  IF(itcnt > 10)GO TO 50
  dt=f1*(t1-t0)/(f1-f0)
  t0=t1
  f0=f1
  t1=t1-dt
  GO TO 90
  50 rl=xlv0-xlv1*t1
!
!...IF RELATIVE HUMIDITY IS SPECIFIED TO BE LESS THAN 100%, ESTIMATE THE
!   TEMPERATURE AND MIXING RATIO WHICH WILL YIELD THE APPROPRIATE VALUE.
!
  IF(rh == 1.)GO TO 110
  dssdt=(cliq-bliq*dliq)/((t1-dliq)*(t1-dliq))
  dt=rl*rs*(1.-rh)/(cp+rl*rh*rs*dssdt)
  t1rh=t1+dt
  es=rh*aliq*EXP((bliq*t1rh-cliq)/(t1rh-dliq))
  rsrh=0.622*es/(p-es)
!
!...CHECK TO SEE IF MIXING RATIO AT SPECIFIED RH IS LESS THAN ACTUAL
!...MIXING RATIO...IF SO, ADJUST TO GIVE ZERO EVAPORATION...
!
  IF(rsrh < rd)THEN
    rsrh=rd
    t1rh=t1+(rs-rsrh)*rl/cp
  END IF
  t1=t1rh
  rs=rsrh
  110 tpdd=t1
  IF(itcnt > 10)PRINT*,'***** NUMBER OF ITERATIONS IN TPDD = ', itcnt
  RETURN
  END FUNCTION tpdd
!***********************************************************************

SUBROUTINE dtfrznew_old(tu,p,thteu,qvap,qliq,qice,ratio2,ttfrz,tbfrz,       &
           qnwfrz,rl,frc1,effq,iflag,xlv0,xlv1,xls0,xls1,               &
           aliq,bliq,cliq,dliq,aice,bice,cice,dice)
!
!...ALLOW GLACIATION OF THE UPDRAFT TO OCCUR AS AN APPROXIMATELY LINEAR
!   FUNCTION OF TEMPERATURE IN THE TEMPERATURE RANGE TTFRZ TO TBFRZ...
!
  rv=461.51
  c5=1.0723E-3
!
!...ADJUST THE LIQUID WATER CONCENTRATIONS FROM FRESH CONDENSATE AND THA
!   BROUGHT UP FROM LOWER LEVELS TO AN AMOUNT THAT WOULD BE PRESENT IF N
!   LIQUID WATER HAD FROZEN THUS FAR...THIS IS NECESSARY BECAUSE THE
!   EXPRESSION FOR TEMP CHANGE IS MULTIPLIED BY THE FRACTION EQUAL TO TH
!   PARCEL TEMP DECREASE SINCE THE LAST MODEL LEVEL DIVIDED BY THE TOTAL
!   GLACIATION INTERVAL, SO THAT EFFECTIVELY THIS APPROXIMATELY ALLOWS A
!   AMOUNT OF LIQUID WATER TO FREEZE WHICH IS EQUAL TO THIS SAME FRACTIO
!   OF THE LIQUID WATER THAT WAS PRESENT BEFORE THE GLACIATION PROCESS W
!   INITIATED...ALSO, TO ALLOW THETAU TO CONVERT APPROXIMATELY LINEARLY
!   ITS VALUE WITH RESPECT TO ICE, WE NEED TO ALLOW A PORTION OF THE FRE
!   CONDENSATE TO CONTRIBUTE TO THE GLACIATION PROCESS; THE FRACTIONAL
!   AMOUNT THAT APPLIES TO THIS PORTION IS 1/2 OF THE FRACTIONAL AMOUNT
!   FROZEN OF THE "OLD" CONDENSATE BECAUSE THIS FRESH CONDENSATE IS ONLY
!   PRODUCED GRADUALLY OVER THE LAYER...NOTE THAT IN TERMS OF THE DYNAMI
!   OF THE PRECIPITATION PROCESS, IE. PRECIPITATION FALLOUT, THIS FRACTI
!   AMNT OF FRESH CONDENSATE HAS ALREADY BEEN INCLUDED IN THE ICE CATEGO
!
  qlqfrz=qliq*effq
  qnew=qnwfrz*effq*0.5
  esliq=aliq*EXP((bliq*tu-cliq)/(tu-dliq))
  esice=aice*EXP((bice*tu-cice)/(tu-dice))
  rlc=2.5E6-2369.276*(tu-273.16)
  rls=2833922.-259.532*(tu-273.16)
  rlf=rls-rlc
  cp=1005.7*(1.+0.89*qvap)
!
!  A = D(ES)/DT IS THAT CALCULATED FROM BUCK'S (1981) EMPERICAL FORMULAS
!  FOR SATURATION VAPOR PRESSURE...
!
  a=(cice-bice*dice)/((tu-dice)*(tu-dice))
  b=rls*0.622/p
  c=a*b*esice/cp
  dqvap=b*(esliq-esice)/(rls+rls*c)-rlf*(qlqfrz+qnew)/(rls+rls/c)
  dtfrz=(rlf*(qlqfrz+qnew)+b*(esliq-esice))/(cp+a*b*esice)
  tu1=tu
  qvap1=qvap
  tu=tu+frc1*dtfrz
  qvap=qvap-frc1*dqvap
  es=qvap*p/(0.622+qvap)
  esliq=aliq*EXP((bliq*tu-cliq)/(tu-dliq))
  esice=aice*EXP((bice*tu-cice)/(tu-dice))
  ratio2=(esliq-es)/(esliq-esice)
!
!  TYPICALLY, RATIO2 IS VERY CLOSE TO (TTFRZ-TU)/(TTFRZ-TBFRZ), USUALLY
!  WITHIN 1% (USING TU BEFORE GALCIATION EFFECTS ARE APPLIED);  IF THE
!  INITIAL UPDRAFT TEMP IS BELOW TBFRZ AND RATIO2 IS STILL LESS THAN 1,
!  AN ADJUSTMENT TO FRC1 AND RATIO2 IS INTRODUCED SO THAT GLACIATION
!  EFFECTS ARE NOT UNDERESTIMATED; CONVERSELY, IF RATIO2 IS GREATER THAN
!  FRC1 IS ADJUSTED SO THAT GLACIATION EFFECTS ARE NOT OVERESTIMATED...
!
  IF(iflag > 0.AND.ratio2 < 1)THEN
    frc1=frc1+(1.-ratio2)
    tu=tu1+frc1*dtfrz
    qvap=qvap1-frc1*dqvap
    ratio2=1.
    iflag=1
    GO TO 20
  END IF
  IF(ratio2 > 1.)THEN
    frc1=frc1-(ratio2-1.)
    frc1=AMAX1(0.0,frc1)
    tu=tu1+frc1*dtfrz
    qvap=qvap1-frc1*dqvap
    ratio2=1.
    iflag=1
  END IF
!
!  CALCULATE A HYBRID VALUE OF THETAU, ASSUMING THAT THE LATENT HEAT OF
!  VAPORIZATION/SUBLIMATION CAN BE ESTIMATED USING THE SAME WEIGHTING
!  FUNCTION AS THAT USED TO CALCULATE SATURATION VAPOR PRESSURE, CALCU-
!  LATE NEW LIQUID WATER AND ICE CONCENTRATIONS...
!
  20 rlc=xlv0-xlv1*tu
  rls=xls0-xls1*tu
  rl=ratio2*rls+(1.-ratio2)*rlc
  pi=(1.e5/p)**(0.2854*(1.-0.28*qvap))
  thteu=tu*pi*EXP(rl*qvap*c5/tu*(1.+0.81*qvap))
  IF(iflag == 1)THEN
    qice=qice+frc1*dqvap+qliq
    qliq=0.
  ELSE
    qice=qice+frc1*(dqvap+qlqfrz)
    qliq=qliq-frc1*qlqfrz
  END IF
  qnwfrz=0.
  RETURN
END SUBROUTINE dtfrznew_old
!**********************************************************************

SUBROUTINE envirtht_old(p1,t1,q1,tht1,r1,rl,                                &
           aliq,bliq,cliq,dliq,aice,bice,cice,dice)
  DATA t00,p00,c1,c2,c3,c4,c5/273.16,1.e5,3374.6525,2.5403,3114.834,    &
       0.278296,1.0723E-3/
!
!  CALCULATE ENVIRONMENTAL EQUIVALENT POTENTIAL TEMPERATURE...
!
  IF(r1 < 1.e-6)THEN
    ee=q1*p1/(0.622+q1)
    tlog=ALOG(ee/aliq)
    tdpt=(cliq-dliq*tlog)/(bliq-tlog)
    tsat=tdpt-(.212+1.571E-3*(tdpt-t00)-4.36E-4*(t1-t00))*(t1-tdpt)
    tht=t1*(p00/p1)**(0.2854*(1.-0.28*q1))
    tht1=tht*EXP((c1/tsat-c2)*q1*(1.+0.81*q1))
  ELSE IF(ABS(r1-1.) < 1.e-6)THEN
    ee=q1*p1/(0.622+q1)
    tlog=ALOG(ee/aice)
    tfpt=(cice-dice*tlog)/(bice-tlog)
    tht=t1*(p00/p1)**(0.2854*(1.-0.28*q1))
    tsat=tfpt-(.182+1.13E-3*(tfpt-t00)-3.58E-4*(t1-t00))*(t1-tfpt)
    tht1=tht*EXP((c3/tsat-c4)*q1*(1.+0.81*q1))
  ELSE
    ee=q1*p1/(0.622+q1)
    tlog=ALOG(ee/aliq)
    tdpt=(cliq-dliq*tlog)/(bliq-tlog)
    tlogic=ALOG(ee/aice)
    tfpt=(cice-dice*tlogic)/(bice-tlogic)
    tht=t1*(p00/p1)**(0.2854*(1.-0.28*q1))
    tsatlq=tdpt-(.212+1.571E-3*(tdpt-t00)-4.36E-4*(t1-t00))*(t1-tdpt    &
           )
    tsatic=tfpt-(.182+1.13E-3*(tfpt-t00)-3.58E-4*(t1-t00))*(t1-tfpt)
    tsat=r1*tsatic+(1.-r1)*tsatlq
    tht1=tht*EXP(rl*q1*c5/tsat*(1.+0.81*q1))
  END IF
  RETURN
END SUBROUTINE envirtht_old
!*********************************************************************

SUBROUTINE condload_old(qliq,qice,wtw,dz,boterm,enterm,rate,qnewlq,         &
           qnewic,qlqout,qicout)
!  9/18/88...THIS PRECIPITATION FALLOUT SCHEME IS BASED ON THE SCHEME US
!  BY OGURA AND CHO (1973).  LIQUID WATER FALLOUT FROM A PARCEL IS CAL-
!  CULATED USING THE EQUATION DQ=-RATE*Q*DT, BUT TO SIMULATE A QUASI-
!  CONTINUOUS PROCESS, AND TO ELIMINATE A DEPENDENCY ON VERTICAL
!  RESOLUTION THIS IS EXPRESSED AS Q=Q*EXP(-RATE*DZ).
  DATA g/9.81/
  qtot=qliq+qice
  qnew=qnewlq+qnewic
!
!  ESTIMATE THE VERTICAL VELOCITY SO THAT AN AVERAGE VERTICAL VELOCITY C
!  BE CALCULATED TO ESTIMATE THE TIME REQUIRED FOR ASCENT BETWEEN MODEL
!  LEVELS...
!
  qest=0.5*(qtot+qnew)
  g1=wtw+boterm-enterm-2.*g*dz*qest/1.5
  IF(g1 < 0.0)g1=0.
  wavg=(SQRT(wtw)+SQRT(g1))/2.
  conv=rate*dz/wavg
!
!  RATIO3 IS THE FRACTION OF LIQUID WATER IN FRESH CONDENSATE, RATIO4 IS
!  THE FRACTION OF LIQUID WATER IN THE TOTAL AMOUNT OF CONDENSATE INVOLV
!  IN THE PRECIPITATION PROCESS - NOTE THAT ONLY 60% OF THE FRESH CONDEN
!  SATE IS IS ALLOWED TO PARTICIPATE IN THE CONVERSION PROCESS...
!
  ratio3=qnewlq/(qnew+1.e-10)
!  OLDQ=QTOT
  qtot=qtot+0.6*qnew
  oldq=qtot
  ratio4=(0.6*qnewlq+qliq)/(qtot+1.e-10)
  qtot=qtot*EXP(-conv)
!
!  DETERMINE THE AMOUNT OF PRECIPITATION THAT FALLS OUT OF THE UPDRAFT
!  PARCEL AT THIS LEVEL...
!
  dq=oldq-qtot
  qlqout=ratio4*dq
  qicout=(1.-ratio4)*dq
!
!  ESTIMATE THE MEAN LOAD OF CONDENSATE ON THE UPDRAFT IN THE LAYER, CAL
!  LATE VERTICAL VELOCITY
!
  pptdrg=0.5*(oldq+qtot-0.2*qnew)
  wtw=wtw+boterm-enterm-2.*g*dz*pptdrg/1.5
!
!  DETERMINE THE NEW LIQUID WATER AND ICE CONCENTRATIONS INCLUDING LOSSE
!  DUE TO PRECIPITATION AND GAINS FROM CONDENSATION...
!
  qliq=ratio4*qtot+ratio3*0.4*qnew
  qice=(1.-ratio4)*qtot+(1.-ratio3)*0.4*qnew
  qnewlq=0.
  qnewic=0.
  RETURN
END SUBROUTINE condload_old
!***********************************************************************
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!  THIS SUBROUTINE INTEGRATES THE AREA UNDER THE CURVE IN THE GAUSSIAN
!  DISTRIBUTION...THE NUMERICAL APPROXIMATION TO THE INTEGRAL IS TAKEN F
!  "HANDBOOK OF MATHEMATICAL FUNCTIONS WITH FORMULAS, GRAPHS AND MATHEMA
!  TABLES" ED. BY ABRAMOWITZ AND STEGUN, NAT'L BUREAU OF STANDARDS APPLI
!  MATHEMATICS SERIES.  JUNE, 1964., MAY, 1968.
!                                  JACK KAIN
!                                  7/6/89
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!*********************************************************************
!*****    GAUSSIAN TYPE MIXING PROFILE....******************************

SUBROUTINE prof5_old(EQ,ee,ud)
  DATA sqrt2p,a1,a2,a3,p,sigma,fe/2.506628,0.4361836,-0.1201676,        &
      0.9372980,0.33267,0.166666667,0.202765151/
  x=(EQ-0.5)/sigma
  y=6.*EQ-3.
  ey=EXP(y*y/(-2))
  e45=EXP(-4.5)
  t2=1./(1.+p*ABS(y))
  t1=0.500498
  c1=a1*t1+a2*t1*t1+a3*t1*t1*t1
  c2=a1*t2+a2*t2*t2+a3*t2*t2*t2
  IF(y >= 0.)THEN
    ee=sigma*(0.5*(sqrt2p-e45*c1-ey*c2)+sigma*(e45-ey))-e45*EQ*EQ/2.
    ud=sigma*(0.5*(ey*c2-e45*c1)+sigma*(e45-ey))-e45*(0.5+EQ*EQ/2.-     &
        EQ)
  ELSE
    ee=sigma*(0.5*(ey*c2-e45*c1)+sigma*(e45-ey))-e45*EQ*EQ/2.
    ud=sigma*(0.5*(sqrt2p-e45*c1-ey*c2)+sigma*(e45-ey))-e45*(0.5+EQ*    &
        EQ/2.-EQ)
  END IF
  ee=ee/fe
  ud=ud/fe
  RETURN
END SUBROUTINE prof5_old
