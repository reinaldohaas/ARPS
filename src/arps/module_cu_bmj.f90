!########################################################################
!########################################################################
!#########                                                      #########
!#########     WRF:MODEL_LAYER:PHYSICS:MODULE module_cu_bmj     #########
!#########                                                      #########
!#########                      Adapted by                      #########
!#########     Center for Analysis and Prediction of Storms     #########
!#########                University of Oklahoma                #########
!#########                                                      #########
!########################################################################
!########################################################################

MODULE module_cu_bmj

!------------------------------------------------------------------------
!
! PURPOSE:
!
! Contains the routines and parameters for the Betts-Miller-Janjic 
! deep and shallow convective adjustment scheme.
!
! REFERENCES:
!
! Janjic, Z. I., 1994:  The step-mountain Eta coordinate model:  Further
!   developments of the convection, viscous sublayer, and turbulence
!   closure schemes.  _Mon. Wea. Rev._, 122, 927-945.
! Betts, A. K., and M. J. Miller, 1993:  The Betts-Miller Scheme.  _The
!   Representation of Cumulus Convection in Numerical Models_.  _Meteor.
!   Mono._, 24, 107-121.
! Janjic, Z. I., 1990:  The step-mountain coordinate:  Physical package.
!   _Mon. Wea. Rev._, 118, 1429-1443.
! Betts, A. K., 1986:  A new convective adjustment scheme.  Part I:  
!   Observational and theroetical basis.  _Quart. J. Roy. Meteor. Soc._,
!   112, 677-691.
! Betts, A. K., and M. J. Miller, 1986:  A new convective adjustment
!   scheme.  Part II:  Single column tests using GATE wave, BOMEX, ATEX
!   and arctic air-mass data sets.  _Quart. J. Roy. Meteor. Soc._, 112,
!   693-709.
! Betts, A. K., 1985:  Mixing line analysis of clouds and cloudy 
!   boundary layers.  _J. Atmos. Sci._, 42, 2751-2763.
! Betts, A. K., 1984:  Boundary layer thermodynamics of a High Plains
!   severe storm.  _Mon. Wea. Rev._, 112, 2199-2211.
! Betts, A. K., 1983a:  Thermodynamics of mixed stratocumulus layers:  
!   Saturation point budgets.  _J. Atmos. Sci._, 40, 2655-2670.
! Betts, A. K., 1983b:  Atmospheric convective structure and a convection
!   scheme based on saturation point adjustment.  Workshop on convection
!   in large-scale models, 20 Nov. to 1 Dec., ECMWF.
! Betts, A. K., 1982a:  Saturation point analysis of moist convective
!   overturning.  _J. Atmos. Sci._, 39, 1484-1505.
! Betts, A. K., 1982b:  Cloud thermodynamic models in saturation point
!   coordinates.  _J. Atmos. Sci._, 39, 2182-2191.
!
!------------------------------------------------------------------------
!
! AUTHOR: ???
!
! MODIFICATIONS:
!
! Eric Kemp, 21 September 2001
! Reformatted code.
!
!------------------------------------------------------------------------
!
! Declare module parameters and variables.
!
!------------------------------------------------------------------------

  LOGICAL :: UNIS=.TRUE.,UNIL=.FALSE.

  REAL,PARAMETER :: A2=17.2693882,A3=273.16,A4=35.86                     &
                       ,DSPC=-3000.                                      &
                       ,DTTOP=0.,EFIFC=5.0,EFIMN=.20,EFMNT=.70           &
                       ,ELIVW=2.72E6                                     &
                       ,EPSDN=1.05,EPSDT=0.,EPSNTP=1.E2                  &
                       ,EPSP=1.E-7,EPSQ=2.E-12,EPSUP=1.0                 &
                       ,FCC=.50,FSL=.850,FSS=.85                         &
                       ,PBM=13000.,PFRZ=15000.,PNO=1000.                 &
                       ,PONE=2500.,PQM=20000.,PQ0=379.90516              &
                       ,PSH=20000.,PSHU=45000.,RHF=0.10,ROW=1.E3         &
                       ,STABD=.90,STABFC=1.0                             &
                       ,STABS=1.0,STRESH=1.10                            &
                       ,T1=274.16,TREL=2400.

  REAL,PARAMETER :: DSPBFL=-4843.75,DSP0FL=-7050.0,DSPTFL=-2250.0        &
                       ,DSPBFS=-3875.,DSP0FS=-5875.,DSPTFS=-1875.

  REAL,PARAMETER :: PL=2500.,PLQ=70000.,PH=105000.                       &
                       ,THL=210.,THH=365.,THHQ=325.

  INTEGER,PARAMETER :: ITB=76,JTB=134,ITBQ=152,JTBQ=440

!------------------------------------------------------------------------
!
! Arrays for lookup tables.
!
!------------------------------------------------------------------------

  REAL,DIMENSION(ITB),PRIVATE,SAVE :: STHE,THE0
  REAL,DIMENSION(JTB),PRIVATE,SAVE :: QS0,SQS
  REAL,DIMENSION(ITBQ),PRIVATE,SAVE :: STHEQ,THE0Q
  REAL,DIMENSION(ITB,JTB),PRIVATE,SAVE :: PTBL
  REAL,DIMENSION(JTB,ITB),PRIVATE,SAVE :: TTBL
  REAL,DIMENSION(JTBQ,ITBQ),PRIVATE,SAVE :: TTBLQ

  REAL,PARAMETER :: RDP=(ITB-1.)/(PH-PL),RDPQ=(ITBQ-1.)/(PH-PLQ)        &
                       ,RDQ=ITB-1,RDTH=(JTB-1.)/(THH-THL)               &
                       ,RDTHE=JTB-1.,RDTHEQ=JTBQ-1.

CONTAINS ! Subroutine declarations

!########################################################################
!########################################################################
!#########                                                      #########
!#########                  SUBROUTINE bmjdrv                   #########
!#########                                                      #########
!#########                      Adapted by                      #########
!#########     Center for Analysis and Prediction of Storms     #########
!#########                University of Oklahoma                #########
!#########                                                      #########
!########################################################################
!########################################################################

SUBROUTINE BMJDRV(DT,ITIMESTEP,STEPCU                                   &
                       ,RR,RTHCUTEN,RQVCUTEN                            &
                       ,RAINCV                                          &
                       ,TH,T,QV,PINT,PMID,PI,RHO,DZ8W                   &
                       ,CP,R,RVOVRD,ELWV,G,TFRZ,D608                    &
                       ,CLDEFI,LOWLYR,XLAND                             &
                       ,IDS,IDE,JDS,JDE,KDS,KDE                         &
                       ,IMS,IME,JMS,JME,KMS,KME                         &
                       ,ITS,ITE,JTS,JTE,KTS,KTE)

!------------------------------------------------------------------------
!
! PURPOSE:
!
! Driver subroutine for Betts-Miller-Janjic deep and shallow convective
! adjustment scheme.
!
!------------------------------------------------------------------------
!
! AUTHOR:  ???
!
! MODIFICATIONS:
!
! Eric Kemp, 21 September 2001
! Reformatted code.
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

!  INTEGER,INTENT(IN) :: IDS,IDE,JDS,JDE,KDS,KDE                        &
!                           ,IMS,IME,JMS,JME,KMS,KME                     & 
!                           ,ITS,ITE,JTS,JTE,KTS,KTE
  INTEGER,INTENT(IN) :: IDS, & ! start index for i in domain
                        IDE, & ! end index for i in domain
                        JDS, & ! start index for j in domain
                        JDE, & ! end index for j in domain
                        KDS, & ! start index for k in domain
                        KDE, & ! end index for k in domain
                        IMS, & ! start index for i in memory
                        IME, & ! end index for i in memory
                        JMS, & ! start index for j in memory
                        JME, & ! end index for j in memory
                        KMS, & ! start index for k in memory
                        KME, & ! end index for k in memory
                        ITS, & ! start index for i in tile
                        ITE, & ! end index for i in tile
                        JTS, & ! start index for j in tile
                        JTE, & ! end index for j in tile
                        KTS, & ! start index for k in tile
                        KTE    ! end index for k in tile

  INTEGER,INTENT(IN) :: ITIMESTEP, &  ! Number of time step (integer)
                        STEPCU        ! Number of fundamental timesteps 
                                      !   between convection calls.

  INTEGER,DIMENSION(IMS:IME,JMS:JME),INTENT(IN) :: &
                     LOWLYR ! Index of lowest model layer above the ground

!  REAL,INTENT(IN) :: CP,DT,ELWV,G,R,RVOVRD,TFRZ,D608
  REAL,INTENT(IN) :: CP, &      ! specific heat at constant pressure 
                                !   (1004 J/k/kg)
                     DT, &      ! Time step (sec)
                     ELWV, &    ! Latent heat of vaporization at 0 deg C
                                !   (2.5e6 J kg**-1)
                     G,  &      ! Acceleration due to gravity 
                                !   (9.81 m s**-2)
                     R,  &      ! Dry gas constant (287 J K**-1 kg**-1)
                     RVOVRD, &  ! Ratio of gas constant for water vapor
                                !   over dry gas constant
                                !   (461.6/287 J K**-1 kg**-1
                     TFRZ, &    ! Freezing point (273.15 K)
                     D608       ! rvovrd - 1.

  REAL,DIMENSION(IMS:IME,JMS:JME),INTENT(IN) :: &
                     XLAND      ! Land-sea mask (1.0 for land; 2.0 for 
                                !   water)

  REAL,DIMENSION(IMS:IME,KMS:KME,JMS:JME),INTENT(IN) :: &
                     TH, &      ! Potential temperature (K)
                     T, &       ! Temperature (K)
                     QV, &      ! Water vapor mixing ratio (kg/kg)
                     RR, &      ! Dry air density (kg/m^3)
                     PINT, &    ! Pressure at full levels (Pa)
                     PMID, &    ! Pressure (Pa)
                     PI, &      ! Exner function (dimensionless)
                     RHO, &     ! Density (kg/m^3)
                     DZ8W       ! dz between full levels (m)

  REAL,DIMENSION(IMS:IME,KMS:KME,JMS:JME),INTENT(INOUT) ::               &
                     RQVCUTEN,& ! Rho_dTheta_m tendency due to
                                !   cumulus scheme precipitation 
                                !   (kg/m^3 . K)
                     RTHCUTEN   ! Rho_dQv tendency due to
                                !   cumulus scheme precipitation 
                                !   (kg/m^3 . kg/kg)

  REAL,DIMENSION(IMS:IME,JMS:JME),INTENT(INOUT) :: &
                     RAINCV, &  ! Cumulus scheme precipitation (mm)
                     CLDEFI     ! Precipitation efficiency (for BMJ 
                                !   scheme) (dimensionless)

!------------------------------------------------------------------------
!
! Local variables.
!
!------------------------------------------------------------------------

  REAL :: CUBOT,CUTOP,DTCNVC,LANDMASK,PCPCOL,PSFC,PTOP

  REAL,DIMENSION(KTS:KTE) :: DPCOL,DQDT,DTDT,PCOL,QCOL,TCOL

  INTEGER :: I,J,K,KFLIP,ICLDCK,LBOT,LMH,LTOP

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! 
! Beginning of executable code...
! 
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!------------------------------------------------------------------------
!
! Check to see if this is a convection timestep.
!
!------------------------------------------------------------------------

  ICLDCK=MOD(ITIMESTEP,STEPCU)                                              

!------------------------------------------------------------------------
!
! Compute convection every stepcu*dt/60.0 minutes.
!
!------------------------------------------------------------------------

  IF (ICLDCK.EQ.0) THEN

    DTCNVC=DT*STEPCU

    DO J=JTS,JTE  
      DO I=ITS,ITE
        DO K=KTS,KTE
          DQDT(K)=0. 
          DTDT(K)=0.
        END DO ! DO k = kts,kte

        RAINCV(I,J)=0.
        PCPCOL=0.
        PSFC=PINT(I,LOWLYR(I,J),J)  ! Surface pressure (Pa)
        PTOP=PINT(I,KTE+1,J)        ! Top level pressure (Pa)  KTE+1=KME

!------------------------------------------------------------------------
!
!       Convert to BMJ land mask.  (1.0 for sea, 0.0 for land)
!
!------------------------------------------------------------------------

        LANDMASK=XLAND(I,J)-1.

!------------------------------------------------------------------------
!
!       Fill 1-D vertical arrays and flip direction since BMJ scheme
!       counts downward from the domain's top.
!
!------------------------------------------------------------------------

        DO K=KTS,KTE
          KFLIP=KTE+1-K
          TCOL(K)=T(I,KFLIP,J)  ! Local temperature (K)

! change 10/18/00 
! from mixing ratio to specific humidity (prevent -ve 4/23/01)

          QCOL(K)=MAX(0.,QV(I,KFLIP,J)/(1.+QV(I,KFLIP,J))) ! Spec. Humidity
          PCOL(K)=PMID(I,KFLIP,J)                      ! Pressure
          DPCOL(K)=PINT(I,KFLIP,J)-PINT(I,KFLIP+1,J)
          DPCOL(K)=RHO(I,KFLIP,J)*G*DZ8W(I,KFLIP,J)   ! Pressure thickness
                                                      ! of layer.
        END DO ! DO k=kts,kte

!------------------------------------------------------------------------
!
!       Lowest layer above ground must also be flipped.
!
!------------------------------------------------------------------------

        LMH=KTE+1-LOWLYR(I,J)

!------------------------------------------------------------------------
!
!       Run BMJ scheme.
!
!------------------------------------------------------------------------

        CALL BMJ(I,J,DTCNVC,LMH,LANDMASK,CLDEFI(I,J)                     &
                  ,DPCOL,PCOL,QCOL,TCOL,PSFC,PTOP                        &
                  ,DQDT,DTDT,PCPCOL,LBOT,LTOP                            &
                  ,CP,R,ELWV,G,TFRZ,D608                                 &
                  ,IDS,IDE,JDS,JDE,KDS,KDE                               &
                  ,IMS,IME,JMS,JME,KMS,KME                               &
                  ,ITS,ITE,JTS,JTE,KTS,KTE)

!------------------------------------------------------------------------
!
!       Compute heating and moistening tendencies.
!
!------------------------------------------------------------------------

        DO K=KTS,KTE
          KFLIP=KTE+1-K

          RTHCUTEN(I,K,J)=RR(I,K,J)*                                     &
                          ((1.+RVOVRD*QV(I,K,J))*DTDT(KFLIP)             &
                                                /PI(I,K,J)               &
                             +RVOVRD*TH(I,K,J)*DQDT(KFLIP))

! change 10/18/00 
! from specific humidity back to mixing ration

          RQVCUTEN(I,K,J)=RR(I,K,J)*DQDT(KFLIP)/(1.-QCOL(KFLIP)**2)
        END DO ! DO k = kts,kte

!------------------------------------------------------------------------
!
!       All units in BMJ scheme are MKS, so convert precip from meters
!       to millimeters per step.
!
!------------------------------------------------------------------------

        RAINCV(I,J)=PCPCOL*1.E3/STEPCU

!------------------------------------------------------------------------
!
!       Variables cubot and cutop are the bottom and top layers, 
!       respectively, of the convective cloud.  If need be, they can be 
!       used in the radiation computations. (Note:  They are not
!       used at this time! EMK)
!
!------------------------------------------------------------------------

        IF (CUBOT.NE.0) CUBOT=REAL(KTE+1-LBOT)
        IF (CUTOP.NE.0) CUTOP=REAL(KTE+1-LTOP)

      END DO ! DO i = its,ite
    END DO ! DO j = jts,jte

  END IF ! IF (ICLDCK.EQ.0) THEN

END SUBROUTINE BMJDRV

!########################################################################
!########################################################################
!#########                                                      #########
!#########                    SUBROUTINE bmj                    #########
!#########                                                      #########
!#########                      Adapted by                      #########
!#########     Center for Analysis and Prediction of Storms     #########
!#########                University of Oklahoma                #########
!#########                                                      #########
!########################################################################
!########################################################################

SUBROUTINE BMJ(I,J,DTCNVC,LMH,SM,CLDEFI                                 &
              ,DPRS,PRSMID,Q,T,PSFC,PT                                  &
              ,DQDT,DTDT,PCPCOL,LBOT,LTOP                               &
              ,CP,R,ELWV,G,TFRZ,D608                                    &
              ,IDS,IDE,JDS,JDE,KDS,KDE                                  &
              ,IMS,IME,JMS,JME,KMS,KME                                  &
              ,ITS,ITE,JTS,JTE,KTS,KTE)

!------------------------------------------------------------------------
!
! PURPOSE:
!
! Runs the Betts-Miller-Janjic deep and shallow convective adjustment
! scheme.
!
!------------------------------------------------------------------------
!
! AUTHOR: ???
!
! MODIFICATION:
!
! Eric Kemp, 21 September 2001
! Reformatted code.
!
!------------------------------------------------------------------------
!
! Force explicit declarations
!
!------------------------------------------------------------------------

  IMPLICIT NONE

!------------------------------------------------------------------------
!
! Declare arguments
!
!------------------------------------------------------------------------

  INTEGER,INTENT(IN) :: IDS, & ! start index for i in domain
                        IDE, & ! end index for i in domain
                        JDS, & ! start index for j in domain
                        JDE, & ! end index for j in domain
                        KDS, & ! start index for k in domain
                        KDE, & ! end index for k in domain
                        IMS, & ! start index for i in memory
                        IME, & ! end index for i in memory
                        JMS, & ! start index for j in memory
                        JME, & ! end index for j in memory
                        KMS, & ! start index for k in memory
                        KME, & ! end index for k in memory
                        ITS, & ! start index for i in tile
                        ITE, & ! end index for i in tile
                        JTS, & ! start index for j in tile
                        JTE, & ! end index for j in tile
                        KTS, & ! start index for k in tile
                        KTE, & ! end index for k in tile
                        I,   & ! I index of column
                        J      ! J index of column

  INTEGER,INTENT(IN) :: LMH    ! Lowest k index above ground

  INTEGER,INTENT(OUT) :: LBOT,&! K level of bottom of convection
                         LTOP  ! K level of top of convection

  REAL,INTENT(IN) :: CP, &     ! Specific heatat constant pressure
                               !   (1004 J/K/kg)
                     D608, &   ! (Rv/Rd) - 1.
                     DTCNVC, & ! Convection time step (stepcu*dt)
                     ELWV, &   ! Latent heat of vaporization at 0 deg C
                               !   (2.5e6 J/kg)
                     G, &      ! Gravitational acceleration (9.81 m/s/s)
                     PSFC, &   ! Surface (lowest layer) pressure (Pa)
                     PT, &     ! Model top pressure (Pa)
                     R, &      ! Dry gas constant (287 J/K/kg)
                     SM, &     ! Land mask (1.0 for sea, 0.0 for land)
                     TFRZ      ! Freezing point (273.15 K)

  REAL,DIMENSION(KTS:KTE),INTENT(IN) :: &
                     DPRS, &   ! Delta pressure across layer (Pa)
                     PRSMID, & ! Pressure at layer (Pa)
                     Q, &      ! Specific humidity (kg/kg)
                     T         ! Temperature (K)

  REAL,INTENT(INOUT) :: &
                     CLDEFI, & ! Precip. efficiency (dimensionless)
                     PCPCOL    ! Convective precip. (m)

  REAL,DIMENSION(KTS:KTE),INTENT(INOUT) :: &
                     DQDT, &   ! Specific humidity time tendency (kg/kg/s)
                     DTDT      ! Temperature time tendency (K/s)

!------------------------------------------------------------------------
!
! Declare local variable arrays.
!
!------------------------------------------------------------------------

  REAL,DIMENSION(KTS:KTE) :: APEK,APESK,FPK                              &
                            ,PK,PSK,QK,QREFK,QSATK                       &
                            ,THERK,THEVRF,THSK                           &
                            ,THVMOD,THVREF,TK,TREFK

  REAL,DIMENSION(KTS:KTE) :: APE,DIFQ,DIFT,TREF

!------------------------------------------------------------------------
!
! Declare local scalar variables.
!
!------------------------------------------------------------------------

 LOGICAL :: DEEP,SHALLOW

 REAL :: APEKL,APEKXX,APEKXY,APESP,APESTS                                &
             ,AVRGT,AVRGTL,BQ,BQK,BQS00K,BQS10K                          &
             ,CAPA,CTHRS,DEN,DENTPY,DEPMIN,DEPTH                         &
             ,DEPWL,DHDT,DIFQL,DIFTL,DPKL,DPMIX                          &
             ,DQREF,DRHDP,DRHEAT,DSP                                     &
             ,DSP0,DSP0K,DSPB,DSPBK,DSPT,DSPTK                           &
             ,DST,DSTQ,DTHEM,DTDP,EFI                                    &
             ,FEFI,FPTK,HCORR,OTSUM,P,P00K,P01K,P10K,P11K                &
             ,PART1,PART2,PART3,PBOT,PBTK                                &
             ,PK0,PKB,PKL,PKT,PKXXXX,PKXXXY                              &
             ,PLMH,POTSUM,PP1,PPK,PRECK                                  &
             ,PRESK,PSP,PSUM,PTHRS,PTOP,PTPK                             &
             ,QBT,QKL,QNEW,QOTSUM,QQ1,QQK,QRFKL,QRFTP,QSUM,RDP0T         &
             ,RDPSUM,RDTCNVC,RHH,RHL,RHMAX,ROTSUM,RTBAR                  &
             ,SMIX,SQ,SQK,SQS00K,SQS10K,STABDL,SUMDE,SUMDP               &
             ,SUMDT,TAUK,TCORR,THBT,THERKX,THERKY                        &
             ,THESP,THSKL,THTPK,THVMKL,TKL                               &
             ,TPSP,TQ,TQK,TREFKX,TRFKL,TSKL,TTH,TTHBT,TTHES,TTHK

  INTEGER :: IQ,IQTB,IT,ITER,ITTB,ITTBK,KB,KNUMH,KNUML                   &
                ,L,L0,L0M1,LB,LBM1,LCOR,LPT1                             &
                ,LQM,LSHU,LTP1,LTP2,LTSH

!------------------------------------------------------------------------
!
! Declare local parameters.
!
!------------------------------------------------------------------------

  REAL,PARAMETER :: DSPBSL=DSPBFL*FSL,DSP0SL=DSP0FL*FSL                  &
                       ,DSPTSL=DSPTFL*FSL                                &
                       ,DSPBSS=DSPBFS*FSS,DSP0SS=DSP0FS*FSS              &
                       ,DSPTSS=DSPTFS*FSS

  REAL,PARAMETER :: AVGEFI=(EFIMN+1.)*.5,STEFI=1.

  REAL,PARAMETER :: SLOPBL=(DSPBFL-DSPBSL)/(1.-EFIMN)                    &
                       ,SLOP0L=(DSP0FL-DSP0SL)/(1.-EFIMN)                &
                       ,SLOPTL=(DSPTFL-DSPTSL)/(1.-EFIMN)                &
                       ,SLOPBS=(DSPBFS-DSPBSS)/(1.-EFIMN)                &
                       ,SLOP0S=(DSP0FS-DSP0SS)/(1.-EFIMN)                &
                       ,SLOPTS=(DSPTFS-DSPTSS)/(1.-EFIMN)                &
                       ,SLOPE=(1.-EFMNT)/(1.-EFIMN)

  REAL :: A23M4L,CPRLG,ELOCP,FCB,RCP

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! 
! Beginning of executable code...
! 
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  CAPA=R/CP
  CPRLG=CP/(ROW*G*ELWV)
  ELOCP=ELIVW/CP
  RCP=1./CP
  A23M4L=A2*(A3-A4)*ELWV
  FCB=1.-FCC
  RDTCNVC=1./DTCNVC

  DEEP=.FALSE.
  SHALLOW=.FALSE.

  DSP0=0.
  DSPB=0.
  DSPT=0.
  PSP=0.

  TAUK=DTCNVC/TREL
  CTHRS=(0.006350/86400.)*DTCNVC/CPRLG

!------------------------------------------------------------------------
!
! Preparations. (???)
!
!------------------------------------------------------------------------

  LBOT=LMH
  THESP=0.
  PCPCOL=0.
  TREF(1)=T(1)

  DO L=1,LMH
    APESTS=PRSMID(L)
    APE(L)=(1.E5/APESTS)**CAPA
  END DO ! DO l = 1,lmh

!------------------------------------------------------------------------
!
! Search for maximum buoyancy level.
!
!------------------------------------------------------------------------

  DO 170 KB=1,LMH

!------------------------------------------------------------------------
!
!   Trial maximum buoyancy level variables.
!
!------------------------------------------------------------------------

    PKL=PRSMID(KB)
    PLMH=PRSMID(LMH)

!------------------------------------------------------------------------
!
!   Search over a scaled depth (between bottom pressure level and 80%
!   of the bottom pressure level) to find the parcel with the maximum
!   equivalent potential temperature.
!
!------------------------------------------------------------------------

    IF (PKL.GE.0.80*PLMH) THEN
      QBT=Q(KB)
      TTHBT=T(KB)*APE(KB)   ! T to Theta
      TTH=(TTHBT-THL)*RDTH   
      QQ1=TTH-AINT(TTH)     ! Remainder used with (bi)linear 
                            ! interpolation.
      ITTB=INT(TTH)+1       ! Look-up table index.

!------------------------------------------------------------------------
!
!     Keep indices within the table. (A look-up table is used for this
!     code --EMK)
!
!------------------------------------------------------------------------

      IF (ITTB.LT.1) THEN
        ITTB=1
        QQ1=0.
      END IF ! IF (ITTB.LT.1) THEN

      IF (ITTB.GE.JTB) THEN
        ITTB=JTB-1
        QQ1=0.
      END IF ! IF (ITTB.GE.JTB) THEN

!------------------------------------------------------------------------
!
!     Base and scaling factor for saturation specific humidity.
!
!------------------------------------------------------------------------

      ITTBK=ITTB
      BQS00K=QS0(ITTBK)
      SQS00K=SQS(ITTBK)
      BQS10K=QS0(ITTBK+1)
      SQS10K=SQS(ITTBK+1)

!------------------------------------------------------------------------
!
!     Scaling specific humidity and table index using linear 
!     interpolation.
!
!------------------------------------------------------------------------

      BQ=(BQS10K-BQS00K)*QQ1+BQS00K
      SQ=(SQS10K-SQS00K)*QQ1+SQS00K
      TQ=(QBT-BQ)/SQ*RDQ
      PP1=TQ-AINT(TQ)      ! Remainder used with bilinear interpolation.
      IQTB=INT(TQ)+1       ! Index in look-up table

!------------------------------------------------------------------------
!
!     Keep indices within the table. (Refers to look-up table -- EMK)
!
!------------------------------------------------------------------------

      IF (IQTB.LT.1) THEN
        IQTB=1
        PP1=0.
      END IF ! IF (IQTB.LT.1) THEN

      IF (IQTB.GE.ITB) THEN
        IQTB=ITB-1
        PP1=0.
      END IF ! IF (IQTB.GE.ITB) THEN

!------------------------------------------------------------------------
!
!     Saturation pressure at four surrounding table points. (Refers to
!     look-up table --EMK)
!
!------------------------------------------------------------------------

      IQ=IQTB
      IT=ITTB
      P00K=PTBL(IQ  ,IT  )
      P10K=PTBL(IQ+1,IT  )
      P01K=PTBL(IQ  ,IT+1)
      P11K=PTBL(IQ+1,IT+1)

!------------------------------------------------------------------------
!
!     Saturation point variables at the bottom.
!
!------------------------------------------------------------------------

      TPSP=P00K+(P10K-P00K)*PP1+(P01K-P00K)*QQ1                          &
               +(P00K-P10K-P01K+P11K)*PP1*QQ1 ! Saturation pressure
      APESP=(1.E5/TPSP)**CAPA                 
      TTHES=TTHBT*EXP(ELOCP*QBT*APESP/TTHBT)  ! Theta-e.

!------------------------------------------------------------------------
!
!     Check for maximum buoyancy.
!
!------------------------------------------------------------------------

      IF (TTHES.GT.THESP) THEN
        PSP  =TPSP
        THBT =TTHBT
        THESP=TTHES
      END IF ! IF (TTHES.GT.THESP) THEN

    END IF ! IF (ITTB.LT.1) THEN

  170 CONTINUE ! DO 170 KB=1,LMH

!------------------------------------------------------------------------
!
! Choose cloud base as model level just below psp.
!
!------------------------------------------------------------------------

  DO L=1,LMH-1
    P=PRSMID(L)
    IF (P.LT.PSP.AND.P.GE.PQM) LBOT=L+1 ! K level of bottom of convection
  END DO ! DO L=1,LMH-1

!------------------------------------------------------------------------
!
! WARNING:  LBOT must not be greater than LMH-1 in shallow convection
! scheme.  Make sure the cloud base is at least PONE above the surface. 
! PONE is currently set to 2500 Pa.
!
!------------------------------------------------------------------------

  PBOT=PRSMID(LBOT) ! Pressure at bottom of convection.
  PLMH=PRSMID(LMH)  ! Bottom level pressure.

  IF (PBOT.GE.PLMH-PONE.OR.LBOT.GE.LMH) THEN
    DO L=1,LMH-1
      P=PRSMID(L)
      IF(P.LT.PLMH-PONE) LBOT=L
    END DO ! DO L=1,LMH-1

    PBOT=PRSMID(LBOT)
  END IF ! IF (PBOT.GE.PLMH-PONE.OR.LBOT.GE.LMH) THEN

!------------------------------------------------------------------------
!
! Cloud top computation.
!
!------------------------------------------------------------------------

  LTOP=LBOT ! Initialize k-index of top of convection
  PTOP=PBOT ! Initialize pressure of top of convection

!------------------------------------------------------------------------
!
! Compute parcel temperature along moist adiabat, scaling pressure
! and TT table index.  (???)
!
!------------------------------------------------------------------------

  DO L=LMH,1,-1

    PRESK=PRSMID(L)

    IF (PRESK.LT.PLQ) THEN ! PLQ = 70000 Pa
      CALL TTBLEX(ITB,JTB,PL,PRESK,RDP,RDTHE                             &
                 ,STHE,THE0,THESP,TTBL,TREF(L))
    ELSE
      CALL TTBLEX(ITBQ,JTBQ,PLQ,PRESK,RDPQ,RDTHEQ                        &
                 ,STHEQ,THE0Q,THESP,TTBLQ,TREF(L))
    END IF ! IF (PRESK.LT.PLQ) THEN

  END DO ! DO L=LMH,1,-1

!------------------------------------------------------------------------
!
! Buoyancy check.
!
!------------------------------------------------------------------------

  DO L=LMH,1,-1
    IF (TREF(L).GT.T(L)-DTTOP) LTOP=L ! K-index of top of cloud
  END DO ! DO L=LMH,1,-1

!------------------------------------------------------------------------
!
! Cloud top pressure.
!
!------------------------------------------------------------------------

  PTOP=PRSMID(LTOP)

!------------------------------------------------------------------------
!
! Define and smooth dsps and cldefi.  Unified or separate land/sea
! conv 
!
! DSPB is the saturation pressure difference at cloud base.
! DSP0 is the saturation pressure difference at freezing level.
! DSPT is the saturation pressure difference at cloud top.
!
!------------------------------------------------------------------------

  EFI=CLDEFI

  IF (UNIS) THEN ! Currently, UNIS is set to true.

!------------------------------------------------------------------------
!
!   All profiles are sea profiles (unified sea) (???)
! 
!------------------------------------------------------------------------

    DSPB=(EFI-EFIMN)*SLOPBS+DSPBSS
    DSP0=(EFI-EFIMN)*SLOP0S+DSP0SS
    DSPT=(EFI-EFIMN)*SLOPTS+DSPTSS

  ELSE IF (.NOT.UNIL) THEN ! Currently, UNIL is set to false

!------------------------------------------------------------------------
! 
!   All profiles are not land profiles (not unified land) (???)
!
!------------------------------------------------------------------------

    DSPB=((EFI-EFIMN)*SLOPBS+DSPBSS)*SM                                  &
            +((EFI-EFIMN)*SLOPBL+DSPBSL)*(1.-SM)
    DSP0=((EFI-EFIMN)*SLOP0S+DSP0SS)*SM                                  &
            +((EFI-EFIMN)*SLOP0L+DSP0SL)*(1.-SM)
    DSPT=((EFI-EFIMN)*SLOPTS+DSPTSS)*SM                                  &
                 +((EFI-EFIMN)*SLOPTL+DSPTSL)*(1.-SM)
  ELSE

!------------------------------------------------------------------------
! 
!   All profiles are land profiles (unified land) (???)
!
!------------------------------------------------------------------------

    DSPB=((EFI-EFIMN)*SLOPBL+DSPBSL)
    DSP0=((EFI-EFIMN)*SLOP0L+DSP0SL)
    DSPT=((EFI-EFIMN)*SLOPTL+DSPTSL) 
  END IF ! IF (UNIS) ... ELSE IF (.NOT. UNIL) ...

!------------------------------------------------------------------------
! 
! Clean up and gather deep convection points.
!
!------------------------------------------------------------------------

  IF (LTOP.GE.LBOT) THEN
    LBOT=0
    LTOP=LBOT
    PTOP=PBOT
  END IF

  IF (PTOP.GT.PBOT-PNO.OR.LTOP.GT.LBOT-2)                                &
       CLDEFI=AVGEFI*SM+STEFI*(1.-SM) ! Recalculate precip. efficiency

!------------------------------------------------------------------------
! 
! Depth of cloud required to make the point a deep convection point is
! a scaled value of psfc.
!
!------------------------------------------------------------------------

  DEPMIN=PSH*PSFC*1.E-5
  DEPTH=PBOT-PTOP

  IF (DEPTH.GE.DEPMIN) THEN
    DEEP=.TRUE.
  END IF ! IF (DEPTH.GE.DEPMIN) THEN

!------------------------------------------------------------------------
!
! Deep convection.
!
!------------------------------------------------------------------------

  IF (.NOT.DEEP) GO TO 600 ! Skip to shallow convection

  LB   =LBOT
  EFI  =CLDEFI
  DSPBK=DSPB 
  DSP0K=DSP0
  DSPTK=DSPT

!------------------------------------------------------------------------
! 
! Initialize variables in the convective column.  One should note that
! the values assigned to the array trefk in the following loop are 
! really only relevant in anchoring the reference temperature profile
! at level lb.  When building the reference profile from cloud base,
! then assigning the ambient temperature to trefk is acceptable.  
! However, when building the reference profile from some other level
! (such as one level above the ground), then trefk should be filled with
! the temperatures in tref(l) which are the temperatures of the moist
! adiabat through cloud base.  By the time the line numbered 450 has 
! been reached, trefk actually does hold the reference temperature
! profile.
!
!------------------------------------------------------------------------

  DO L=1,LMH
    DIFT  (L)=0.
    DIFQ  (L)=0.
    TKL      =T(L)
    TK    (L)=TKL           ! Temperature
    TREFK (L)=TKL           ! Temperature
    QKL      =Q(L)
    QK    (L)=QKL           ! Specific humidity
    QREFK (L)=QKL           ! Specific humidity
    PKL      =PRSMID(L)
    PK    (L)=PKL           ! Pressure
    PSK   (L)=PKL           ! Pressure
    APEKL    =APE(L)
    APEK  (L)=APEKL
    THERK (L)=TREF(L)*APEKL  ! Temperature to Theta
  END DO ! DO L=1,LMH

!------------------------------------------------------------------------
!
! Deep convection reference temperature profile.
!
!------------------------------------------------------------------------
  
  LTP1=LTOP+1    ! K-index of one level below cloud top
  LBM1=LB-1      ! K-index of one level above cloud bottom
  PKB=PK(LB)     ! Pressure at cloud bottom
  PKT=PK(LTOP)   ! Pressure at cloud top

!------------------------------------------------------------------------
!
! Temperature reference profile below freezing level.
!
!------------------------------------------------------------------------

  L0=LB
  PK0=PK(LB)
  TREFKX=TREFK(LB)     ! Temperature at cloud bottom
  THERKX=THERK(LB)     ! Theta at cloud bottom
  APEKXX=APEK(LB)
  THERKY=THERK(LBM1)   ! Theta at one level above cloud bottom
  APEKXY=APEK(LBM1)

  DO L=LBM1,LTOP,-1
    IF (T(L+1).LT.TFRZ) GO TO 430 ! Skip to code for temperature profile
                                  ! above freezing level.
    STABDL=STABD
    TREFKX=((THERKY-THERKX)*STABDL                                       &
                +TREFKX*APEKXX)/APEKXY  ! Reference profile temperature
    TREFK(L)=TREFKX                     ! Reference profile temperature

    APEKXX=APEKXY
    THERKX=THERKY                       ! Change lower level theta to
                                        ! upper level theta.
    APEKXY=APEK(L-1)
    THERKY=THERK(L-1)                   ! Change upper level theta to
                                        ! theta from next highest level.
    L0=L
    PK0=PK(L0)                          ! Pressure
  END DO ! DO L=LBM1,LTOP,-1

!------------------------------------------------------------------------
!
! Freezing level at or above the cloud top.
!
!------------------------------------------------------------------------

  L0M1=L0-1 
  GO TO 450 ! Skip to deep convection reference humidity profile.

!------------------------------------------------------------------------
!
! Temperature reference profile above freezing level.
!
!------------------------------------------------------------------------

  430 L0M1=L0-1
  RDP0T=1./(PK0-PKT) ! Reciprocal of difference in pressure between 
                     ! current level and cloud top.
  DTHEM=THERK(L0)-TREFK(L0)*APEK(L0) ! Difference between theta and
                                     ! reference theta.

  DO L=LTOP,L0M1
    TREFK(L)=(THERK(L)-(PK(L)-PKT)*DTHEM*RDP0T)/APEK(L) ! Reference 
                                                        ! profile
                                                        ! temperature.
  END DO ! DO L=LTOP,L0M1

!------------------------------------------------------------------------
!
! Deep convection reference humidity profile.
!
! DEPWL is the pressure difference between cloud base and the freezing
! level.
!
!------------------------------------------------------------------------

  450 DEPWL=PKB-PK0
  DEPTH=PFRZ*PSFC*1.E-5 ! 15000 Pa * lowest layer pressure * 1.E-5

!------------------------------------------------------------------------
!
! Saturation pressure difference.
!
!------------------------------------------------------------------------

  DO 460 L=LTOP,LB ! Cloud top to cloud base.

    IF (DEPWL.GE.DEPTH) THEN
      IF (L.LT.L0) THEN ! L0 is k-index of freezing level
        DSP=((PK0-PK(L))*DSPTK+(PK(L)-PKT)*DSP0K)/(PK0-PKT)
      ELSE
        DSP=((PKB-PK(L))*DSP0K+(PK(L)-PK0)*DSPBK)/(PKB-PK0)
      END IF ! IF (L.LT.L0) THEN
    ELSE
      DSP=DSP0K
      IF (L.LT.L0)                                                       &
        DSP=((PK0-PK(L))*DSPTK+(PK(L)-PKT)*DSP0K)/(PK0-PKT)
    END IF ! IF (DEPWL.GE.DEPTH) THEN

!------------------------------------------------------------------------
!
!   Humidity profile.
!
!------------------------------------------------------------------------

    IF (PK(L).GT.PQM) THEN ! PQM = 20000 Pa
      PSK(L)=PK(L)+DSP
      APESK(L)=(1.E5/PSK(L))**CAPA
      THSK(L)=TREFK(L)*APEK(L)     ! Reference profile theta
      QREFK(L)=PQ0/PSK(L)*EXP(A2*(THSK(L)-A3*APESK(L))                   &
                                  /(THSK(L)-A4*APESK(L))) ! Reference
                                                          ! profile
                                                          ! sp. humidity.
    ELSE
      QREFK(L)=Q(L) ! Reference profile specific humidity.
    END IF ! IF (PK(L).GT.PQM) THEN

  460 CONTINUE ! DO 460 L=LTOP,LB

!------------------------------------------------------------------------
!
! Enthalpy conservation integral.
!
!------------------------------------------------------------------------

  LQM=0

  DO 520 ITER=1,2
    SUMDE=0. ! Integral of enthalpy multiplied by delta pressure
             ! across each layer.
    SUMDP=0. ! Integral of pressure from cloud base to cloud top. 
    DO L=LTOP,LB ! Cloud top to cloud base
      SUMDE=((TK(L)-TREFK(L))*CP+(QK(L)-QREFK(L))*ELWV)*DPRS(L)          &
              +SUMDE
      SUMDP=SUMDP+DPRS(L)
    END DO ! DO L=LTOP,LB

    HCORR=SUMDE/(SUMDP-DPRS(LTOP)) ! Integrated enthalpy.
    LCOR=LTOP+1 ! One level below cloud top height.

!------------------------------------------------------------------------
! 
!   Find LQM.
!       
!------------------------------------------------------------------------

    DO L=1,LB
      IF (PK(L).LE.PQM) LQM=L ! PQM = 20000 Pa
    END DO ! DO L=1,LB

!------------------------------------------------------------------------
! 
!   Below lqm (above 200 mb level) correct temperature only.
!       
!------------------------------------------------------------------------

    IF (LCOR.LE.LQM) THEN
      DO L=LCOR,LQM
        TREFK(L)=TREFK(L)+HCORR*RCP
      END DO ! DO L=LCOR,LQM
      LCOR=LQM+1
    END IF ! IF (LCOR.LE.LQM) THEN

!------------------------------------------------------------------------
! 
!   Above lqm (below 200 mb level) correct both temperature and moisture.
!       
!------------------------------------------------------------------------

    DO L=LCOR,LB ! LQM + 1 to cloud bottom               
      TSKL=TREFK(L)*APEK(L)/APESK(L)
      DHDT=QREFK(L)*A23M4L/(TSKL-A4)**2+CP
      TREFK(L)=HCORR/DHDT+TREFK(L)
      THSKL=TREFK(L)*APEK(L)
      QREFK(L)=PQ0/PSK(L)*EXP(A2*(THSKL-A3*APESK(L))                     &
                                  /(THSKL-A4*APESK(L)))
    END DO ! DO L=LCOR,LB

  520 CONTINUE ! DO 520 ITER=1,2

!------------------------------------------------------------------------
! 
! Heating, moistening, precipitation.  (???)
!       
!------------------------------------------------------------------------

  DENTPY=0. ! Integral of difference in entropy between actual and
            ! reference thermodynamic profile.
  AVRGT =0.
  PRECK =0. ! Condensate/precipitation

  DO L=LTOP,LB ! Cloud top to cloud base
    TKL    =TK(L)
    DIFTL  =(TREFK(L)-TKL  )*TAUK ! TAUK=DTCNVC/TREL
    DIFQL  =(QREFK(L)-QK(L))*TAUK ! TAUK=DTCNVC/TREL
    AVRGTL =(TKL+TKL+DIFTL)
    DENTPY =(DIFTL*CP+DIFQL*ELWV)*DPRS(L)/AVRGTL+DENTPY
    AVRGT  =AVRGTL*DPRS(L)+AVRGT
    PRECK  =DPRS(L)*DIFTL+PRECK
    DIFT(L)=DIFTL
    DIFQ(L)=DIFQL
  END DO ! DO L=LTOP,LB

  DENTPY=DENTPY+DENTPY
  AVRGT =AVRGT/(SUMDP+SUMDP)

!------------------------------------------------------------------------
!
! Swap if entropy and/or precip is less than zero.
!
!------------------------------------------------------------------------

  IF (DENTPY.LT.EPSNTP.OR.PRECK.LT.0.) THEN ! EPSNTP = 1.E2
    CLDEFI=EFIMN*SM+STEFI*(1.-SM) ! New precip. efficiency

!------------------------------------------------------------------------
!
!   Search for shallow cloud top.
!
!------------------------------------------------------------------------

    LTSH=LBOT                 ! K-index of cloud bottom
    LBM1=LBOT-1               ! K-index of level above cloud bottom
    PBTK=PK(LBOT)             ! Pressure at cloud bottom
    DEPMIN=PSH*PSFC*1.E-5
    PTPK=PBTK-DEPMIN          ! Est. pressure at shallow cloud top

!------------------------------------------------------------------------
!
!   Cloud top is the level just below pbtk-psh.
!
!------------------------------------------------------------------------

    DO L=1,LMH
      IF(PK(L).LE.PTPK)LTOP=L+1 ! K-index of shallow cloud top
    END DO ! DO L=1,LMH

    PTPK=PK(LTOP) ! Pressure of shallow cloud top

!------------------------------------------------------------------------
!
!   Highest level allowed is level just below pshu.
!
!------------------------------------------------------------------------

    IF (PTPK.LE.PSHU) THEN ! PSHU = 45000 Pa

      DO L=1,LMH
        IF (PK(L).LE.PSHU) LSHU=L+1
      END DO ! DO L=1,LMH

      LTOP=LSHU            ! K-index of shallow cloud top
      PTPK=PK(LTOP)        ! Pressure of shallow cloud top
    END IF ! IF (PTPK.LE.PSHU) THEN

    LTP1=LTOP+1    ! K-index of first level below shallow cloud top
    LTP2=LTOP+2    ! K-index of second level below shallow cloud top

    DO L=LTOP,LBOT
      QSATK(L)=PQ0/PK(L)*EXP(A2*(TK(L)-A3)/(TK(L)-A4)) ! Saturation
                                                       ! sp. humidity.
    END DO ! DO L=LTOP,LBOT

    RHH=QK(LTOP)/QSATK(LTOP) ! Relative humidity at shallow cloud top
    RHMAX=0.

    DO L=LTP1,LBM1 ! Cloud top to cloud bottom
      RHL=QK(L)/QSATK(L)     ! Relative humidity
      DRHDP=(RHH-RHL)/(PK(L-1)-PK(L)) ! Lapse rate of relative humidity

      IF (DRHDP.GT.RHMAX) THEN
        LTSH=L-1     ! K-index of maximum relative humidity
        RHMAX=DRHDP
      END IF ! IF (DRHDP.GT.RHMAX) THEN

      RHH=RHL
    END DO ! DO L=LTP1,LBM1

    LTOP=LTSH ! K-index of shallow cloud top

!------------------------------------------------------------------------
!
!   Cloud must be at least two layers thick.
!
!------------------------------------------------------------------------

    IF (LBOT-LTOP.LT.2) LTOP=LBOT-2

    PTOP=PK(LTOP) ! Pressure of shallow cloud top
    GO TO 600 ! Go to shallow convection
  END IF ! IF (DENTPY.LT.EPSNTP.OR.PRECK.LT.0.) THEN

!------------------------------------------------------------------------
!
! Deep convection otherwise.
!
! Keep the land value of efi equal to 1 until precip surpasses
! a threshold value, currently set to 0.25 in per 24 hrs.
!
!------------------------------------------------------------------------

  PTHRS=CTHRS
  DRHEAT=(PRECK*SM+AMAX1(EPSP,PRECK-PTHRS)*(1.-SM))                      &
             *CP/AVRGT ! Latent heat release
  EFI=EFIFC*DENTPY/DRHEAT ! Precip. efficiency

!------------------------------------------------------------------------
!
! Unified or separate land/sea conv. (???)
!
!------------------------------------------------------------------------

  EFI=CLDEFI*FCB+EFI*FCC ! Precip. efficiency, between 0.20 and 1.

  IF(EFI.GT.1.   )EFI=1.
  IF(EFI.LT.EFIMN)EFI=EFIMN
  IF(PRECK.EQ.0.) EFI=1.
  CLDEFI=EFI

  FEFI=EFMNT+SLOPE*(EFI-EFIMN)
!  FEFI=AMAX1(EFI,EFMNT)
!
  PRECK=PRECK*FEFI

!------------------------------------------------------------------------
!
! Precipitation, temperature and moisture changes.
!
!------------------------------------------------------------------------

  PCPCOL=PRECK*CPRLG ! Precipitation accumuluation.

  DO L=LTOP,LB
    DTDT(L)=DIFT(L)*FEFI*RDTCNVC  ! Temperature tendency
    DQDT(L)=DIFQ(L)*FEFI*RDTCNVC  ! Moisture tendency
  END DO ! DO L=LTOP,LB

  600 CONTINUE ! End of deep convection

!------------------------------------------------------------------------
!
! Gather shallow convection points.
!
!------------------------------------------------------------------------

  IF (PTOP.LE.PBOT-PNO.AND.LTOP.LE.LBOT-2) THEN
    DEPMIN=PSH*PSFC*1.E-5 ! Minimum depth requirement

    IF (PTOP+1..GE.PBOT-DEPMIN) THEN
      SHALLOW=.TRUE.
    END IF ! IF (PTOP+1..GE.PBOT-DEPMIN) THEN

  END IF ! IF (PTOP.LE.PBOT-PNO.AND.LTOP.LE.LBOT-2) THEN

  IF (.NOT.SHALLOW) GO TO 800 ! Skip shallow convection

!------------------------------------------------------------------------
!
! Shallow convection.
!
!------------------------------------------------------------------------

  DO L=1,LMH
    TKL      =T(L)  
    TK   (L) =TKL           ! Temperature
    TREFK(L) =TKL           ! Temperature
    QKL      =Q(L)          
    QK   (L) =QKL           ! Specific humidity
    QREFK(L) =QKL           ! Specific humidity
    PKL      =PRSMID(L)
    PK   (L) =PKL           ! Pressure
    QSATK(L) =PQ0/PK(L)*EXP(A2*(TK(L)-A3)/(TK(L)-A4))
                            ! Saturation specific humidity
    APEKL    =APE(L)
    APEK (L) =APEKL
    THVMKL   =TKL*APEKL*(QKL*D608+1.) ! Theta-v
    THVREF(L)=THVMKL                  ! Theta-v
  END DO ! DO L=1,LMH

!------------------------------------------------------------------------
!
! Shallow cloud top.
!
!------------------------------------------------------------------------

  LBM1=LBOT-1     ! One k-level above cloud bottom
  PTPK=PTOP       ! Pressure at cloud top
  LTP1=LTOP-1     ! One k-level above cloud top

! NOTE:  This appears to be redundant...EMK
  IF (PTOP.GT.PBOT-PNO.OR.LTOP.GT.LBOT-2) THEN
    LBOT=0
    LTOP=LBOT
    PTOP=PBOT
    GO TO 800 ! Skip shallow convection
  END IF ! IF (PTOP.GT.PBOT-PNO.OR.LTOP.GT.LBOT-2) THEN

!------------------------------------------------------------------------
!
! Scaling potential temperature and table index at top.
!
!------------------------------------------------------------------------

  THTPK=T(LTP1)*APE(LTP1) ! Theta at one level above cloud top

  TTHK=(THTPK-THL)*RDTH   
  QQK =TTHK-AINT(TTHK)    ! Remainder for linear interpolation.
  IT  =INT(TTHK)+1        ! Index for look-up table

  IF (IT.LT.1) THEN
    IT=1
    QQK=0.
  END IF ! IF (IT.LT.1) THEN

  IF (IT.GE.JTB) THEN
    IT=JTB-1
    QQK=0.
  END IF ! IF (IT.GE.JTB) THEN

!------------------------------------------------------------------------
!
! Base and scaling factor for saturation specific humidity at top.
!
!------------------------------------------------------------------------

  BQS00K=QS0(IT)
  SQS00K=SQS(IT)
  BQS10K=QS0(IT+1)
  SQS10K=SQS(IT+1)

!------------------------------------------------------------------------
!
! Scaling saturation specific humidity and table index at top.
!
!------------------------------------------------------------------------

  BQK=(BQS10K-BQS00K)*QQK+BQS00K
  SQK=(SQS10K-SQS00K)*QQK+SQS00K

!  TQK=(Q(LTOP)-BQK)/SQK*RDQ
  TQK=(Q(LTP1)-BQK)/SQK*RDQ

  PPK=TQK-AINT(TQK)             ! Remainder for linear interpolation.
  IQ =INT(TQK)+1                ! Index for look-up table.

  IF (IQ.LT.1) THEN
    IQ=1
    PPK=0.
  END IF ! IF (IQ.LT.1) THEN

  IF (IQ.GE.ITB) THEN
    IQ=ITB-1
    PPK=0.
  END IF ! IF (IQ.GE.ITB) THEN

!------------------------------------------------------------------------
!
! Cloud top saturation point pressure.
!
!------------------------------------------------------------------------

  PART1=(PTBL(IQ+1,IT)-PTBL(IQ,IT))*PPK
  PART2=(PTBL(IQ,IT+1)-PTBL(IQ,IT))*QQK
  PART3=(PTBL(IQ  ,IT  )-PTBL(IQ+1,IT  )                                 &
         -PTBL(IQ  ,IT+1)+PTBL(IQ+1,IT+1))*PPK*QQK
  PTPK=PTBL(IQ,IT)+PART1+PART2+PART3 ! Saturation point pressure


  DPMIX=PTPK-PSP ! Difference in saturation pressure along mixing line 
                 ! between cloud base and cloud top.
  IF (ABS(DPMIX).LT.3000.) DPMIX=-3000.

!------------------------------------------------------------------------
!
! Temperature profile slope.
!
!------------------------------------------------------------------------

  SMIX=(THTPK-THBT)/DPMIX*STABS

  TREFKX=TREFK(LBOT+1)  ! Reference profile temperature one level below
                        ! cloud bottom.
  PKXXXX=PK(LBOT+1)     ! Pressure one level below cloud bottom.
  PKXXXY=PK(LBOT)       ! Pressure at cloud bottom.
  APEKXX=APEK(LBOT+1)
  APEKXY=APEK(LBOT)

  DO L=LBOT,LTOP,-1
        TREFKX=((PKXXXY-PKXXXX)*SMIX                                     &
                +TREFKX*APEKXX)/APEKXY
    TREFK(L)=TREFKX ! Reference profile temperature.
    APEKXX=APEKXY
    PKXXXX=PKXXXY   ! Switch lower-level pressure with upper-level
                    ! pressure.
    APEKXY=APEK(L-1)
    PKXXXY=PK(L-1)  ! Switch upper-level pressure with pressure at
                    ! next highest level.
  END DO ! DO L=LBOT,LTOP,-1

!------------------------------------------------------------------------
!
! Temperature reference profile correction.
!
!------------------------------------------------------------------------

  SUMDT=0. ! Integral of temperature difference times pressure thickness
           ! of layers.
  SUMDP=0. ! Integral of thickness of pressure of layers.

  DO L=LTOP,LBOT
    SUMDT=(TK(L)-TREFK(L))*DPRS(L)+SUMDT
    SUMDP=SUMDP+DPRS(L)
  END DO ! DO L=LTOP,LBOT

  RDPSUM=1./SUMDP
  FPK(LBOT)=TREFK(LBOT)

  TCORR=SUMDT*RDPSUM ! Integral of temperature difference.

  DO L=LTOP,LBOT
    TRFKL   =TREFK(L)+TCORR
    TREFK(L)=TRFKL   ! Reference profile temperature.
    FPK  (L)=TRFKL
  END DO ! DO L=LTOP,LBOT

!------------------------------------------------------------------------
!
! Humidity profile equations.
!
!------------------------------------------------------------------------

  PSUM  =0. ! Integral of pressure differences 1 (pressure of layer and
            ! pressure of cloud top) times pressure difference 2 
            ! (thickness of layers).
  QSUM  =0. ! Integral of specific humidity times pressure thickness
            ! of layers.
  POTSUM=0.
  QOTSUM=0.
  OTSUM =0.
  DST   =0. ! Change in entropy.
  FPTK  =FPK(LTOP) ! Reference profile temperature at cloud top.

  DO L=LTOP,LBOT ! Cloud top to cloud base
    DPKL  =FPK(L)-FPTK
    PSUM  =DPKL *DPRS(L)+PSUM
    QSUM  =QK(L)*DPRS(L)+QSUM
    RTBAR =2./(TREFK(L)+TK(L))
    OTSUM =DPRS(L)*RTBAR+OTSUM
    POTSUM=DPKL    *RTBAR*DPRS(L)+POTSUM
    QOTSUM=QK(L)   *RTBAR*DPRS(L)+QOTSUM
    DST   =(TREFK(L)-TK(L))*RTBAR*DPRS(L)+DST
  END DO ! DO L=LTOP,LBOT

  PSUM  =PSUM*RDPSUM
  QSUM  =QSUM*RDPSUM
  ROTSUM=1./OTSUM
  POTSUM=POTSUM*ROTSUM
  QOTSUM=QOTSUM*ROTSUM
  DST   =DST*ROTSUM*CP/ELWV ! Change in entropy

!------------------------------------------------------------------------
!
! Ensure positive entropy change.
!
!------------------------------------------------------------------------

  IF (DST.GT.0.) THEN

!    DSTQ=DST*EPSUP
    LBOT=0
    LTOP=LBOT
    PTOP=PBOT
    GO TO 800 ! Skip shallow convection

  ELSE
    DSTQ=DST*EPSDN ! EPSDN = 1.05
  END IF ! IF (DST.GT.0.) THEN

!------------------------------------------------------------------------
!
! Check for isothermal atmosphere.
!
!------------------------------------------------------------------------

  DEN=POTSUM-PSUM

  IF (-DEN/PSUM.LT.5.E-5) THEN
    LBOT=0
    LTOP=LBOT
    PTOP=PBOT
    GO TO 800 ! Skip shallow convection

  ELSE

!------------------------------------------------------------------------
!
! Slope of the reference humidity profile.
!
!------------------------------------------------------------------------

    DQREF=(QOTSUM-DSTQ-QSUM)/DEN
  END IF ! IF (-DEN/PSUM.LT.5.E-5) THEN

!------------------------------------------------------------------------
!
! Humidity does not increase with height.
!
!------------------------------------------------------------------------

  IF (DQREF.LT.0.) THEN
    LBOT=0
    LTOP=LBOT
    PTOP=PBOT
    GO TO 800 ! Skip shallow convection
  END IF ! IF (DQREF.LT.0.) THEN

!------------------------------------------------------------------------
!
! Humidity at the cloud top.
!
!------------------------------------------------------------------------

  QRFTP=QSUM-DQREF*PSUM

!------------------------------------------------------------------------
!
! Humidity profile.
!
!------------------------------------------------------------------------

  DO L=LTOP,LBOT ! Cloud top to cloud base
    QRFKL=(FPK(L)-FPTK)*DQREF+QRFTP

!------------------------------------------------------------------------
!
!   Supersaturation or negative q not allowed.
!
!------------------------------------------------------------------------

    QNEW=(QRFKL-QK(L))*TAUK+QK(L) ! New specific humidity.

    IF (QNEW.GT.QSATK(L)*STRESH .OR. QNEW.LT.0.) THEN
      LBOT=0
      LTOP=LBOT
      PTOP=PBOT
      GO TO 800 ! Skip shallow convection
    END IF ! IF (QNEW.GT.QSATK(L)*STRESH .OR. QNEW.LT.0.) THEN

    THVREF(L)=TREFK(L)*APEK(L)*(QRFKL*D608+1.) ! New theta-v.
    QREFK(L)=QRFKL  ! New reference profile specific humidity.
  END DO ! DO L=LTOP,LBOT

!------------------------------------------------------------------------
!
! Eliminate impossible slopes.  (Betts, dtheta/dq)
!
!------------------------------------------------------------------------

  DO L=LTOP,LBOT ! Cloud top to cloud base.
    DTDP=(THVREF(L-1)-THVREF(L))/(PRSMID(L)-PRSMID(L-1)) ! Dthetav/DP

    IF (DTDP.LT.EPSDT) THEN ! EPSDT = 0.
      LBOT=0
      LTOP=LBOT
      PTOP=PBOT
      GO TO 800 ! Skip shallow convection
    END IF ! IF (DTDP.LT.EPSDT) THEN

  END DO ! DO L=LTOP,LBOT

  DENTPY=0.

  DO L=LTOP,LBOT ! Cloud top to cloud base.
    DENTPY=((TREFK(L)-TK(L))*CP+(QREFK(L)-QK(L))*ELWV)                   &
               /(TK(L)+TREFK(L))*DPRS(L)+DENTPY  ! Change in entropy
  END DO ! DO L=LTOP,LBOT 

!------------------------------------------------------------------------
!
! Relaxation toward reference profiles.
!
!------------------------------------------------------------------------

  DO L=LTOP,LBOT ! Cloud top to cloud base
    DTDT(L)=(TREFK(L)-TK(L))*TAUK*RDTCNVC  ! Temperature tendency
    DQDT(L)=(QREFK(L)-QK(L))*TAUK*RDTCNVC  ! Specific humidity tendency
  END DO ! DO L=LTOP,LBOT

!------------------------------------------------------------------------
!
! End of shallow convection.
!
!------------------------------------------------------------------------

  800 CONTINUE
END SUBROUTINE BMJ

!########################################################################
!########################################################################
!#########                                                      #########
!#########                  SUBROUTINE ttblex                   #########
!#########                                                      #########
!#########                      Adapted by                      #########
!#########     Center for Analysis and Prediction of Storms     #########
!#########                University of Oklahoma                #########
!#########                                                      #########
!########################################################################
!########################################################################

SUBROUTINE TTBLEX(ITBX,JTBX,PLX,PRSMID,RDPX,RDTHEX,STHE,THE0,THESP,TTBL &
                  ,TREF)

!------------------------------------------------------------------------
!
! PURPOSE:
!
! Extract temperature of the moist adiabat from the appropriate ttbl.
!
!------------------------------------------------------------------------
!
! AUTHOR: ???
!
! MODIFICATIONS:
!
! Eric Kemp, 24 September 2001
! Reformatted code.
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

  INTEGER,INTENT(IN) :: ITBX,JTBX  ! Dimensions of look-up tables

  REAL,INTENT(IN) :: PLX, &    ! Pressure reference used with look-up
                               ! table calculation.
                     PRSMID, & ! Grid point pressure (Pa).
                     RDPX, &   ! Remainder used with look-up table
                               ! calculation. 
                     RDTHEX, & ! Remainder used with look-up table
                               ! calculation.
                     THESP     ! Saturation point potential temperature 
                               ! (K).

  REAL,DIMENSION(ITBX),INTENT(IN) :: STHE, & ! Look-up table.
                                     THE0    ! Look-up table.

  REAL,DIMENSION(JTBX,ITBX),INTENT(IN) :: TTBL ! Look-up table.

  REAL,INTENT(OUT) :: TREF ! Parcel temperature (K)

!------------------------------------------------------------------------
!
! Internal variables.
!
!------------------------------------------------------------------------

  REAL :: BTHE00K,BTHE10K,BTHK,PK,PP,QQ,STHE00K,STHE10K,STHK             &
          ,T00K,T01K,T10K,T11K,TPK,TTHK

  INTEGER :: IPTB,ITHTB

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! 
! Beginning of executable code...
! 
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!------------------------------------------------------------------------
!
! Scaling pressure and TT table index.
!
! Determine pressure index for look-up table and remainder for (bi)linear
! interpolation.
!
!------------------------------------------------------------------------

  PK=PRSMID
  TPK=(PK-PLX)*RDPX
  QQ=TPK-AINT(TPK) ! Remainder for linear interpolation.
  IPTB=INT(TPK)+1  ! Index for look-up table.

!------------------------------------------------------------------------
!
! Keeping indices within the table.
!
!------------------------------------------------------------------------

  IF (IPTB.LT.1) THEN
    IPTB=1  ! Index for look-up table.
    QQ=0.   ! Remainder for linear interpolation.
  END IF ! IF (IPTB.LT.1) THEN

  IF (IPTB.GE.ITBX) THEN
    IPTB=ITBX-1  ! Index for look-up table.
    QQ=0.        ! Remainder for linear interpolation.
  END IF ! IF (IPTB.GE.ITBX) THEN

!------------------------------------------------------------------------
!
! Base and scaling factor for thetae.
!
!------------------------------------------------------------------------

  BTHE00K=THE0(IPTB)
  STHE00K=STHE(IPTB)
  BTHE10K=THE0(IPTB+1)
  STHE10K=STHE(IPTB+1)

!------------------------------------------------------------------------
!
! Scaling the and TT table index. (???)
!
! Use linear interpolation to determine the scale and base factors for
! thetae, then use them to determine theta-e index for look-up table
! and remainder for bilinear interpolation.
!
!------------------------------------------------------------------------

  BTHK=(BTHE10K-BTHE00K)*QQ+BTHE00K
  STHK=(STHE10K-STHE00K)*QQ+STHE00K
  TTHK=(THESP-BTHK)/STHK*RDTHEX
  PP=TTHK-AINT(TTHK)   ! Remainder for linear interpolation.
  ITHTB=INT(TTHK)+1    ! Index for look-up table.

!------------------------------------------------------------------------
!
! Keeping indices within the table.
!
!------------------------------------------------------------------------

  IF (ITHTB.LT.1) THEN
    ITHTB=1            ! Index for look-up table.
    PP=0.              ! Remainder for linear interpolation.
  END IF ! IF (ITHTB.LT.1) THEN

  IF (ITHTB.GE.JTBX) THEN
    ITHTB=JTBX-1       ! Index for look-up table.
    PP=0.              ! Remainder for linear interpolation.
  END IF ! IF (ITHTB.GE.JTBX) THEN

!------------------------------------------------------------------------
!
! Temperature at four surrounding TT table points.
!
!------------------------------------------------------------------------

  T00K=TTBL(ITHTB,IPTB)
  T10K=TTBL(ITHTB+1,IPTB)
  T01K=TTBL(ITHTB,IPTB+1)
  T11K=TTBL(ITHTB+1,IPTB+1)

!------------------------------------------------------------------------
!
! Parcel temperature, using bilinear interpolation.
!
!------------------------------------------------------------------------

  TREF=(T00K+(T10K-T00K)*PP+(T01K-T00K)*QQ                               &
       +(T00K-T10K-T01K+T11K)*PP*QQ)

END SUBROUTINE TTBLEX

!########################################################################
!########################################################################
!#########                                                      #########
!#########                  SUBROUTINE bmjinit                  #########
!#########                                                      #########
!#########                      Adapted by                      #########
!#########     Center for Analysis and Prediction of Storms     #########
!#########                University of Oklahoma                #########
!#########                                                      #########
!########################################################################
!########################################################################

SUBROUTINE BMJINIT(RTHCUTEN,RQVCUTEN,RQCCUTEN,RQRCUTEN                  &
                        ,CLDEFI,LOWLYR,CP,RD,restart                    &
                        ,IDS,IDE,JDS,JDE,KDS,KDE                        &
                        ,IMS,IME,JMS,JME,KMS,KME                        &
                        ,ITS,ITE,JTS,JTE,KTS,KTE)

!------------------------------------------------------------------------
!
! PURPOSE:
!
! Initialize variables and look-up tables used by Betts-Miller-Janjic 
! convective scheme.
!
!------------------------------------------------------------------------
!
! AUTHOR: ???
!
! MODIFICATIONS:
!
! Eric Kemp, 24 September 2001
! Reformatted code.
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

  LOGICAL , INTENT(IN) :: restart ! Restart flag.

  INTEGER,INTENT(IN) :: IDS, & ! start index for i in domain
                        IDE, & ! end index for i in domain
                        JDS, & ! start index for j in domain
                        JDE, & ! end index for j in domain
                        KDS, & ! start index for k in domain
                        KDE, & ! end index for k in domain
                        IMS, & ! start index for i in memory
                        IME, & ! end index for i in memory
                        JMS, & ! start index for j in memory
                        JME, & ! end index for j in memory
                        KMS, & ! start index for k in memory
                        KME, & ! end index for k in memory
                        ITS, & ! start index for i in tile
                        ITE, & ! end index for i in tile
                        JTS, & ! start index for j in tile
                        JTE, & ! end index for j in tile
                        KTS, & ! start index for k in tile
                        KTE    ! end index for k in tile

  REAL,INTENT(IN) :: CP, & ! specific heat at constant pressure
                           !   (1004 J/k/kg)
                     RD    ! Dry air gas constant.

  REAL,DIMENSION(IMS:,KMS:,JMS:),INTENT(OUT) ::                          &
            RTHCUTEN, & ! Rho_dTheta_m tendency due to  
                        !   cumulus scheme precipitation
                        !   (kg/m^3 . K)
            RQVCUTEN, & ! Rho_dQv tendency due to
                        !   cumulus scheme precipitation
                        !   (kg/m^3 . kg/kg)
            RQCCUTEN, & ! Rho_dQc tendency due to
                        !   cumulus scheme precipitation
                        !   (kg/m^3 . kg/kg)
            RQRCUTEN    ! Rho_dQr tendency due to
                        !   cumulus scheme precipitation
                        !   (kg/m^3 . kg/kg)

  REAL,DIMENSION(IMS:,JMS:),INTENT(OUT) :: &
            CLDEFI      ! Precipitation efficiency (for BMJ scheme) 
                        !   (dimensionless)

  INTEGER,DIMENSION(IMS:,JMS:),INTENT(OUT) :: &
            LOWLYR      ! Index of lowest model layer above the ground

!------------------------------------------------------------------------
!
! Declare internal parameters.
!
!------------------------------------------------------------------------

  REAL,PARAMETER :: ELIWV=2.683E6,EPS=1.E-9

!------------------------------------------------------------------------
!
! Declare internal variables.
!
!------------------------------------------------------------------------

  REAL, DIMENSION(JTB) :: APP,APT,AQP,AQT,PNEW,POLD,QSNEW,QSOLD          &
                             ,THENEW,THEOLD,TNEW,TOLD,Y2P,Y2T

  REAL,DIMENSION(JTBQ) :: APTQ,AQTQ,THENEWQ,THEOLDQ                      &
                             ,TNEWQ,TOLDQ,Y2TQ

  INTEGER :: I,J,K,ITF,JTF,KTF
  INTEGER :: KTH,KTHM,KTHM1,KP,KPM,KPM1

  REAL :: APE,DP,DQS,DTH,DTHE,P,QS,QS0K,SQSK,STHEK                       &
             ,TH,THE0K

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! 
! Beginning of executable code...
! 
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
 
  JTF=MIN0(JTE,JDE-1)
  KTF=MIN0(KTE,KDE-1)
  ITF=MIN0(ITE,IDE-1)

  IF (.not.restart) THEN
    DO J=JTS,JTF
      DO K=KTS,KTF
        DO I=ITS,ITF
          RTHCUTEN(I,K,J)=0.
          RQVCUTEN(I,K,J)=0.
          RQCCUTEN(I,K,J)=0.
          RQRCUTEN(I,K,J)=0.
        END DO ! DO I=ITS,ITF
      END DO ! DO K=KTS,KTF
    END DO ! DO J=JTS,JTF

    DO J=JTS,JTF
      DO I=ITS,ITF
        CLDEFI(I,J)=1.
      END DO ! DO I=ITS,ITF
    END DO ! DO J=JTS,JTF
  END IF ! IF (.not.restart) THEN

!------------------------------------------------------------------------
!
! For now, assume sigma mode for lowest model layer.
!
!------------------------------------------------------------------------

  DO J=JTS,JTF
    DO I=ITS,ITF
      LOWLYR(I,J)=1
    END DO ! DO I=ITS,ITF
  END DO ! DO J=JTS,JTF

!------------------------------------------------------------------------
!
! Coarse look-up table for saturation point.
!
!------------------------------------------------------------------------

  KTHM=JTB
  KPM=ITB
  KTHM1=KTHM-1
  KPM1=KPM-1

  DTH=(THH-THL)/REAL(KTHM-1)
  DP =(PH -PL )/REAL(KPM -1)

  TH=THL-DTH

  DO 100 KTH=1,KTHM

    TH=TH+DTH
    P=PL-DP

    DO KP=1,KPM
      P=P+DP
      APE=(100000./P)**(RD/CP)
      QSOLD(KP)=PQ0/P*EXP(A2*(TH-A3*APE)/(TH-A4*APE)) ! Saturation 
                                                      ! specific humidity.
      POLD(KP)=P
    END DO ! DO KP=1,KPM

    QS0K=QSOLD(1)
    SQSK=QSOLD(KPM)-QSOLD(1) ! Scaling factor.
    QSOLD(1  )=0.
    QSOLD(KPM)=1.

    DO KP=2,KPM1
      QSOLD(KP)=(QSOLD(KP)-QS0K)/SQSK ! Base look-up table
      IF ((QSOLD(KP)-QSOLD(KP-1)).LT.EPS) QSOLD(KP)=QSOLD(KP-1)+EPS
    END DO ! DO KP=2,KPM1

    QS0(KTH)=QS0K  ! Base factor look-up table.
    SQS(KTH)=SQSK  ! Scaling factor look-up table.

    QSNEW(1  )=0.
    QSNEW(KPM)=1.
    DQS=1./REAL(KPM-1)

    DO KP=2,KPM1
      QSNEW(KP)=QSNEW(KP-1)+DQS
    END DO ! DO KP=2,KPM1

    Y2P(1   )=0.
    Y2P(KPM )=0.

    CALL SPLINE(JTB,KPM,QSOLD,POLD,Y2P,KPM,QSNEW,PNEW,APP,AQP)

    DO KP=1,KPM
      PTBL(KP,KTH)=PNEW(KP)
    END DO ! DO KP=1,KPM

  100 CONTINUE ! DO 100 KTH=1,KTHM

!------------------------------------------------------------------------
!
! Coarse look-up table for T(P) from constant THE.  (???)
!
!------------------------------------------------------------------------

  P=PL-DP

  DO 200 KP=1,KPM

    P=P+DP
    TH=THL-DTH

    DO KTH=1,KTHM
      TH=TH+DTH
      APE=(1.E5/P)**(RD/CP)
      QS=PQ0/P*EXP(A2*(TH-A3*APE)/(TH-A4*APE))
      TOLD(KTH)=TH/APE
      THEOLD(KTH)=TH*EXP(ELIWV*QS/(CP*TOLD(KTH)))
    END DO ! DO KTH=1,KTHM

    THE0K=THEOLD(1)
    STHEK=THEOLD(KTHM)-THEOLD(1)
    THEOLD(1   )=0.
    THEOLD(KTHM)=1.

    DO KTH=2,KTHM1
      THEOLD(KTH)=(THEOLD(KTH)-THE0K)/STHEK
      IF((THEOLD(KTH)-THEOLD(KTH-1)).LT.EPS)                             &
            THEOLD(KTH)=THEOLD(KTH-1)  +  EPS
    END DO ! DO KTH=2,KTHM1

    THE0(KP)=THE0K
    STHE(KP)=STHEK

    THENEW(1  )=0.
    THENEW(KTHM)=1.
    DTHE=1./REAL(KTHM-1)

    DO KTH=2,KTHM1
      THENEW(KTH)=THENEW(KTH-1)+DTHE
    END DO ! DO KTH=2,KTHM1

    Y2T(1   )=0.
    Y2T(KTHM)=0.

    CALL SPLINE(JTB,KTHM,THEOLD,TOLD,Y2T,KTHM,THENEW,TNEW,APT,AQT)

    DO KTH=1,KTHM
      TTBL(KTH,KP)=TNEW(KTH)
    END DO ! DO KTH=1,KTHM

  200 CONTINUE ! DO 200 KP=1,KPM

!------------------------------------------------------------------------
!
! Fine look-up table for saturation point.
!
!------------------------------------------------------------------------

  KTHM=JTBQ
  KPM=ITBQ
  KTHM1=KTHM-1
  KPM1=KPM-1

  DTH=(THHQ-THL)/REAL(KTHM-1)
  DP=(PH-PLQ)/REAL(KPM-1)

  TH=THL-DTH
  P=PLQ-DP

!------------------------------------------------------------------------
!
! Fine look-up table for T(P) from constant THE. (???)
!
!------------------------------------------------------------------------

  DO 300 KP=1,KPM

    P=P+DP
    TH=THL-DTH

    DO KTH=1,KTHM
      TH=TH+DTH
      APE=(1.E5/P)**(RD/CP)
      QS=PQ0/P*EXP(A2*(TH-A3*APE)/(TH-A4*APE))
      TOLDQ(KTH)=TH/APE
      THEOLDQ(KTH)=TH*EXP(ELIWV*QS/(CP*TOLDQ(KTH)))
    END DO ! DO KTH=1,KTHM

    THE0K=THEOLDQ(1)
    STHEK=THEOLDQ(KTHM)-THEOLDQ(1)
    THEOLDQ(1   )=0.
    THEOLDQ(KTHM)=1.

    DO KTH=2,KTHM1
      THEOLDQ(KTH)=(THEOLDQ(KTH)-THE0K)/STHEK
      IF ((THEOLDQ(KTH)-THEOLDQ(KTH-1)).LT.EPS)                          &
            THEOLDQ(KTH)=THEOLDQ(KTH-1)  +  EPS
    END DO ! DO KTH=2,KTHM1

    THE0Q(KP)=THE0K
    STHEQ(KP)=STHEK

    THENEWQ(1  )=0.
    THENEWQ(KTHM)=1.
    DTHE=1./REAL(KTHM-1)

    DO KTH=2,KTHM1
      THENEWQ(KTH)=THENEWQ(KTH-1)+DTHE
    END DO ! DO KTH=2,KTHM1

    Y2TQ(1   )=0.
    Y2TQ(KTHM)=0.

    CALL SPLINE(JTBQ,KTHM,THEOLDQ,TOLDQ,Y2TQ,KTHM                        &
                 ,THENEWQ,TNEWQ,APTQ,AQTQ)

    DO KTH=1,KTHM
      TTBLQ(KTH,KP)=TNEWQ(KTH)
    END DO ! DO KTH=1,KTHM

  300 CONTINUE

END SUBROUTINE BMJINIT

!########################################################################
!########################################################################
!#########                                                      #########
!#########                  SUBROUTINE spline                   #########
!#########                                                      #########
!#########                      Adapted by                      #########
!#########     Center for Analysis and Prediction of Storms     #########
!#########                University of Oklahoma                #########
!#########                                                      #########
!########################################################################
!########################################################################

SUBROUTINE SPLINE(JTBX,NOLD,XOLD,YOLD,Y2,NNEW,XNEW,YNEW,P,Q)

!------------------------------------------------------------------------
!
! PURPOSE:
!
! 1-D cubic spline fitting subroutine.
!
!------------------------------------------------------------------------
!
! AUTHOR: Z. Janjic
!
! MODIFICATIONS:
!
! Eric Kemp, 24 September 2001
! Reformatted code.  However, the overall "spagetti" structure of
! the code has been left undisturbed.
!
!------------------------------------------------------------------------
!
! Original documentation:
!
!   ******************************************************************
!   *                                                                *
!   *  THIS IS A ONE-DIMENSIONAL CUBIC SPLINE FITTING ROUTINE        *
!   *  PROGRAMED FOR A SMALL SCALAR MACHINE.                         *
!   *                                                                *
!   *  PROGRAMER Z. JANJIC                                           *
!   *                                                                *
!   *  NOLD - NUMBER OF GIVEN VALUES OF THE FUNCTION.  MUST BE GE 3. *
!   *  XOLD - LOCATIONS OF THE POINTS AT WHICH THE VALUES OF THE     *
!   *         FUNCTION ARE GIVEN.  MUST BE IN ASCENDING ORDER.       *
!   *  YOLD - THE GIVEN VALUES OF THE FUNCTION AT THE POINTS XOLD.   *
!   *  Y2   - THE SECOND DERIVATIVES AT THE POINTS XOLD.  IF NATURAL *
!   *         SPLINE IS FITTED Y2(1)=0. AND Y2(NOLD)=0. MUST BE      *
!   *         SPECIFIED.                                             *
!   *  NNEW - NUMBER OF VALUES OF THE FUNCTION TO BE CALCULATED.     *
!   *  XNEW - LOCATIONS OF THE POINTS AT WHICH THE VALUES OF THE     *
!   *         FUNCTION ARE CALCULATED.  XNEW(K) MUST BE GE XOLD(1)   *
!   *         AND LE XOLD(NOLD).                                     *
!   *  YNEW - THE VALUES OF THE FUNCTION TO BE CALCULATED.           *
!   *  P, Q - AUXILIARY VECTORS OF THE LENGTH NOLD-2.                *
!   *                                                                *
!   ******************************************************************
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

  INTEGER,INTENT(IN) :: JTBX,NNEW,NOLD
  REAL,DIMENSION(JTBX),INTENT(IN) :: XNEW,XOLD,YOLD
  REAL,DIMENSION(JTBX),INTENT(INOUT) :: P,Q,Y2
  REAL,DIMENSION(JTBX),INTENT(OUT) :: YNEW

!------------------------------------------------------------------------
!
! Declare internal variables.
!
!------------------------------------------------------------------------

  INTEGER :: K,K1,K2,KOLD,NOLDM1
  REAL :: AK,BK,CK,DEN,DX,DXC,DXL,DXR,DYDXL,DYDXR                        &
             ,RDX,RTDXC,X,XK,XSQ,Y2K,Y2KP1

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! 
! Beginning of executable code...
! 
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  NOLDM1=NOLD-1

  DXL=XOLD(2)-XOLD(1)
  DXR=XOLD(3)-XOLD(2)
  DYDXL=(YOLD(2)-YOLD(1))/DXL
  DYDXR=(YOLD(3)-YOLD(2))/DXR
  RTDXC=0.5/(DXL+DXR)

  P(1)= RTDXC*(6.*(DYDXR-DYDXL)-DXL*Y2(1))
  Q(1)=-RTDXC*DXR

  IF (NOLD.EQ.3) GO TO 150 ! Skip down

  K=3

  100 DXL=DXR
  DYDXL=DYDXR
  DXR=XOLD(K+1)-XOLD(K)
  DYDXR=(YOLD(K+1)-YOLD(K))/DXR
  DXC=DXL+DXR
  DEN=1./(DXL*Q(K-2)+DXC+DXC)

  P(K-1)= DEN*(6.*(DYDXR-DYDXL)-DXL*P(K-2))
  Q(K-1)=-DEN*DXR

  K=K+1
  IF (K.LT.NOLD) GO TO 100 ! Jump back up

  150 K=NOLDM1

  200 Y2(K)=P(K-1)+Q(K-1)*Y2(K+1)

  K=K-1
  IF (K.GT.1) GO TO 200 ! Jump back up

  K1=1

  300 XK=XNEW(K1)

  DO 400 K2=2,NOLD

    IF (XOLD(K2).GT.XK) THEN
      KOLD=K2-1
      GO TO 450 ! Exit loop, jump down
    END IF ! IF (XOLD(K2).GT.XK) THEN

  400 CONTINUE ! DO 400 K2=2,NOLD

  YNEW(K1)=YOLD(NOLD)
  GO TO 600 ! Jump down

  450 IF (K1.EQ.1) GO TO 500 ! Jump down
  IF (K.EQ.KOLD) GO TO 550 ! Jump down

  500 K=KOLD

  Y2K=Y2(K)
  Y2KP1=Y2(K+1)
  DX=XOLD(K+1)-XOLD(K)
  RDX=1./DX

  AK=.1666667*RDX*(Y2KP1-Y2K)
  BK=0.5*Y2K
  CK=RDX*(YOLD(K+1)-YOLD(K))-.1666667*DX*(Y2KP1+Y2K+Y2K)

  550 X=XK-XOLD(K)
  XSQ=X*X

  YNEW(K1)=AK*XSQ*X+BK*XSQ+CK*X+YOLD(K)

  600 K1=K1+1
  IF (K1.LE.NNEW) GO TO 300 ! Jump up

END SUBROUTINE SPLINE

!------------------------------------------------------------------------
!
! End subroutine declarations.
!
!------------------------------------------------------------------------

END MODULE MODULE_CU_BMJ

