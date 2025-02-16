!
!###############################################################################
!
!   /////////////////////         BEGIN            \\\\\\\\\\\\\\\\\\\\
!   \\\\\\\\\\\\\\\\\\\\\   SUBROUTINE LFO_DRIVER  ////////////////////
!
! LFO_ICE_DRIVE interfaces COMMAS model with Straka's 3-class LFO ice scheme
!
!###############################################################################

SUBROUTINE lfo_ice_drive(tb, qb,             & ! tbar, qvbar  (pib here is pb in ncommas)
              ptbar,ptprt,qv,qscalar,        & ! perturbation temp. and pres., water species
              qvbar,pprt,pbar,               & ! qvbar,pbar,pprt
              dt,dzc,rhostr,j3inv,rhobar,    & ! constants, rhostr, j3inverse, rhobar
              nx, ny, nz, ng)                  ! grid dimensions

!-------------------------------------------------------------------------------
!  PURPOSE
!
!  calculate and apply the microphysical contributions to the water,
!  ice and temperature fields, using an ice microphysics parameterization
!  scheme.
!-------------------------------------------------------------------------------
!
!  MODIFICATION HISTORY:
!
!  4/18/05 (Nate Snook)
!  Updated original code to use standard arps variables and perform
!  calculations in columns rather than x-z slices.
!
!  8/8/05 (Nate Snook)
!  Fixed error in density calculations--code now correctly uses rhobar
!  instead of rhostr.
!
!-------------------------------------------------------------------------------
!
! INPUT
!
!    nx       number of grid points in the x-direction (east/west)
!    ny       number of grid points in the y-direction (north/south)
!    nz       number of grid points in the vertical
!    ng       number of buffer points
!
!    dt       the large time step size for this call.
!
!    pib      base state Exner function        (no longer used - ~NS 4/18/05)
!    pin      perturb Exner function           (no longer used - ~NS 4/18/05)
!    pb       base state pressure function     (replaced by pbar - ~NS 4/18/05)
!    pn       perturbation pressure function   (replaced by pprt - ~NS 4/18/05)
!
!    tb       base state potential temperature (K)
!    qb       base state mixing ratio (kg/kg)
!    t        potential temperature  (K)
!    qv       water vapor specific humidity at all time levels (kg/kg)
!    qc       cloud water mixing ratio at all time levels (kg/kg)
!    qr       rainwater mixing ratio at all time levels (kg/kg)
!    qi       cloud ice mixing ratio at all time levels (kg/kg)
!    qs       snow mixing ratio at all time levels (kg/kg)
!    qg       hail mixing ratio at all time levels (kg/kg)
!    precip   accumulated grid-scale rainfall (mm)
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

  INCLUDE 'globcst.inc'

!-------------------------------------------------------------------------------
!
!  variable declarations.
!
!-------------------------------------------------------------------------------

  INTEGER :: nx,ny,nz       ! number of grid points in 3 directions
  INTEGER :: ng             ! number of buffer points

  REAL :: dt                ! the large time step size for this call.
  REAL :: dzc(nz), cqfac


  REAL :: tb (nz)           ! base state potential temperature (k)
  REAL :: qb (nz)           ! base state water vapor specific humidity (kg/kg)
  REAL :: qscalar(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng,nscalar)
  REAL :: pbar(nx,ny,nz)         !MSG base state pressure
  REAL :: pprt(nx,ny,nz)         !MSG perturbation pressure
  REAL :: ppi(nx,ny,nz)          !MSG exner function
  REAL :: rhostr(nx,ny,nz)       !MSG base state density (multiplied by j3)
  REAL :: j3inv(nx,ny,nz)        !MSG j3 inverse
  REAL :: rhobar(nx,ny,nz)       !Base state density (from rhostr and j3inv)
  REAL :: ptbar(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)   !MSG base state potential temperature
  REAL :: ptprt(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)   !MSG perturbation potential temperature
  REAL :: qv(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)      !MSG vapor mixing ratioy
  REAL :: qvbar(nx,ny,nz)                              !MSG base state vapor mixing ratio

!-------------------------------------------------------------------------------
!
!  temporary arrays
!
!-------------------------------------------------------------------------------
  INTERFACE
    SUBROUTINE lfo_ice(nx,ny,nz,nor,dtp,dzc,ptbar,ptprt,qv,qscalar,  &
                       qvbar,rhobar,pprt,pbar)

      IMPLICIT NONE

      INCLUDE 'globcst.inc'


      INTEGER, INTENT(IN) :: nx,ny,nz,nor
      REAL,    INTENT(IN) :: dtp
      REAL,    INTENT(IN) :: dzc(nz)    ! MSG weren't defined before

      REAL :: ptbar(-nor+1:nx+nor,-nor+1:ny+nor,-nor+1:nz+nor)
                                    ! MSG ptbar (base state potential temp)
      REAL :: ptprt(-nor+1:nx+nor,-nor+1:ny+nor,-nor+1:nz+nor)
                                    ! MSG ptprt (perturbation potential temp)
      REAL :: qv(-nor+1:nx+nor,-nor+1:ny+nor,-nor+1:nz+nor)
                                    ! MSG qv    (vapor mixing ratio)
      REAL :: qvbar(nx,ny,nz)       ! MSG qvbar (base state vapor mixing ratio)

      REAL, TARGET :: qscalar(-nor+1:nx+nor,-nor+1:ny+nor,-nor+1:nz+nor,nscalar)

      REAL :: pprt(-nor+1:nx+nor,-nor+1:ny+nor,-nor+1:nz+nor)
                                    ! MSG pprt  (perturbation pressure)
      REAL :: pbar(nx,ny,nz)        ! MSG pbar  (base state pressure)
      REAL :: rhobar(nx,ny,nz)      ! MSG rhobar (base state density = rhostr*j3inv)

    END SUBROUTINE lfo_ice
  END INTERFACE
!-------------------------------------------------------------------------------
!
!  misc. local variables
!
!-------------------------------------------------------------------------------

  INTEGER :: i,j,k

  INTEGER, SAVE      :: constset = 0
  INTEGER, PARAMETER :: hole_fill = 0

  REAL*8 cqtotn, cqvtotn, cqctotn, cqrtotn, cqitotn, cqstotn, cqhtotn
  REAL*8 cqtotp, cqvtotp, cqctotp, cqrtotp, cqitotp, cqstotp, cqhtotp
  REAL*8 cqvfac,  cqcfac,  cqrfac,  cqifac,  cqsfac,  cqhfac

  INCLUDE 'mp.inc'

!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code below ...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!-----------------------------------------------------------------------
!
! Since this scheme comtain microphysical ice process, it must at least
! define qc, qr, qi, qs, qh
!
!-----------------------------------------------------------------------

  IF (P_QC < 1 .OR. P_QR < 1 .OR. P_QI < 1 .OR. P_QS < 1 .OR. P_QH < 1) THEN

    WRITE(6,'(2a,/,5(a,I2),/,a)')                                       &
               'No enough microphysical array was defined ',            &
               'inside subroutine lfo_ice_drive.',                      &
               'P_QC = ',P_QC,' P_QR = ',P_QR,' P_QI = ',P_QI,          &
              ' P_QS = ',P_QS,' P_QH = ',P_QH,                          &
               'Program aborting ...'
    CALL arpsstop('Wrong size for microphysics array, qscalar.',1)

  END IF


! Convert from rhostr to rhobar using j3inv ~NS (8/8/05)

  DO k = 1, nz-1
    DO j = 1, ny-1
      DO i = 1, nx-1
        rhobar(i,j,k) = rhostr(i,j,k) *j3inv(i,j,k)
      END DO
    END DO
  END DO

! Print microphysics header

  IF (constset == 0) THEN

    IF (myproc == 0) THEN
      WRITE(6,*)
      WRITE(6,*) '2004 STRAKA(GILMORE) SAM microphysics'
      WRITE(6,*)
      IF( hole_fill == 0 ) WRITE(6,*) 'hole_filling OFF '
      IF( hole_fill == 1 ) WRITE(6,*) 'global HOLE FILLING TURN ON'
      IF( hole_fill == 2 ) WRITE(6,*) 'species HOLE FILLING TURN ON'
      WRITE(6,*)
    END IF

    constset = 1

  END IF

  CALL lfo_ice (nx,ny,nz,ng,dt,dzc,ptbar,ptprt,qv,qscalar,     &
                qvbar,rhobar,pprt,pbar)                ! DO MICROPHYSICS

!-----------------------------------------------------------------------
!
! HOLE FILLING SCHEME:  WATER MASS IS CONSERVED
!
!-----------------------------------------------------------------------

  IF( hole_fill == 1 ) THEN       ! GLOBAL METHOD ACROSS ALL SPECIES

    cqtotn = 0.0
    cqtotp = 0.0

    DO k =  1,nz-1
      DO j =  1,ny-1
        DO i =  1,nx-1

          cqtotn = cqtotn + rhobar(i,j,k)*(qv(i,j,k)+qscalar(i,j,k,P_QC) &
                 + qscalar(i,j,k,P_QR) + qscalar(i,j,k,P_QI)             &
                 + qscalar(i,j,k,P_QS) + qscalar(i,j,k,P_QH))

          qv(i,j,k) = MAX(qv(i,j,k),0.0)
          qscalar(i,j,k,P_QC) = MAX(qscalar(i,j,k,P_QC),0.0)
          qscalar(i,j,k,P_QR) = MAX(qscalar(i,j,k,P_QR),0.0)
          qscalar(i,j,k,P_QI) = MAX(qscalar(i,j,k,P_QI),0.0)
          qscalar(i,j,k,P_QS) = MAX(qscalar(i,j,k,P_QS),0.0)
          qscalar(i,j,k,P_QH) = MAX(qscalar(i,j,k,P_QH),0.0)

          cqtotp = cqtotp + rhobar(i,j,k)*(qv(i,j,k)+qscalar(i,j,k,P_QC) &
                            + qscalar(i,j,k,P_QR)+qscalar(i,j,k,P_QI)    &
                            + qscalar(i,j,k,P_QS)+qscalar(i,j,k,P_QH))

        END DO
      END DO
    END DO

! CREATE ratio

    cqfac = (cqtotn+1.0E-20)/(cqtotp+1.0E-20)

! ADJUST fields

    DO k =  1,nz-1
      DO j =  1,ny-1
        DO i =  1,nx-1

          qv(i,j,k) = qv(i,j,k) * cqfac
          qscalar(i,j,k,P_QC) = qscalar(i,j,k,P_QC) * cqfac
          qscalar(i,j,k,P_QR) = qscalar(i,j,k,P_QR) * cqfac
          qscalar(i,j,k,P_QI) = qscalar(i,j,k,P_QI) * cqfac
          qscalar(i,j,k,P_QS) = qscalar(i,j,k,P_QS) * cqfac
          qscalar(i,j,k,P_QH) = qscalar(i,j,k,P_QH) * cqfac

        END DO
      END DO
    END DO

  END IF

!-------------------------------------------------------------------------------
! HOLE_FILL == 2 --> LOCAL FILL METHOD FROM EACH SPECIES

  IF( hole_fill == 2 ) THEN

    cqvtotn = 0.0
    cqvtotp = 0.0
    cqctotn = 0.0
    cqctotp = 0.0
    cqrtotn = 0.0
    cqrtotp = 0.0
    cqitotn = 0.0
    cqitotp = 0.0
    cqstotn = 0.0
    cqstotp = 0.0
    cqhtotn = 0.0
    cqhtotp = 0.0

    DO i =  1,nx-1
      DO k =  1,nz-1
        DO j =  1,ny-1

          cqvtotn = cqvtotn + rhobar(i,j,k)*qv(i,j,k)
          cqctotn = cqctotn + rhobar(i,j,k)*qscalar(i,j,k,P_QC)
          cqrtotn = cqrtotn + rhobar(i,j,k)*qscalar(i,j,k,P_QR)
          cqitotn = cqitotn + rhobar(i,j,k)*qscalar(i,j,k,P_QI)
          cqstotn = cqstotn + rhobar(i,j,k)*qscalar(i,j,k,P_QS)
          cqhtotn = cqhtotn + rhobar(i,j,k)*qscalar(i,j,k,P_QH)

          qv(i,j,k) = MAX(qv(i,j,k),0.0)
          qscalar(i,j,k,P_QC) = MAX(qscalar(i,j,k,P_QC),0.0)
          qscalar(i,j,k,P_QR) = MAX(qscalar(i,j,k,P_QR),0.0)
          qscalar(i,j,k,P_QI) = MAX(qscalar(i,j,k,P_QI),0.0)
          qscalar(i,j,k,P_QS) = MAX(qscalar(i,j,k,P_QS),0.0)
          qscalar(i,j,k,P_QH) = MAX(qscalar(i,j,k,P_QH),0.0)

          cqvtotp = cqvtotp + rhobar(i,j,k)*qv(i,j,k)
          cqctotp = cqctotp + rhobar(i,j,k)*qscalar(i,j,k,P_QC)
          cqrtotp = cqrtotp + rhobar(i,j,k)*qscalar(i,j,k,P_QR)
          cqitotp = cqitotp + rhobar(i,j,k)*qscalar(i,j,k,P_QI)
          cqstotp = cqstotp + rhobar(i,j,k)*qscalar(i,j,k,P_QS)
          cqhtotp = cqhtotp + rhobar(i,j,k)*qscalar(i,j,k,P_QH)

        END DO
      END DO
    END DO

! CREATE ratio

    cqvfac = (cqvtotn+1.0E-20)/(cqvtotp+1.0E-20)
    cqcfac = (cqctotn+1.0E-20)/(cqctotp+1.0E-20)
    cqrfac = (cqrtotn+1.0E-20)/(cqrtotp+1.0E-20)
    cqifac = (cqitotn+1.0E-20)/(cqitotp+1.0E-20)
    cqsfac = (cqstotn+1.0E-20)/(cqstotp+1.0E-20)
    cqhfac = (cqhtotn+1.0E-20)/(cqhtotp+1.0E-20)

! ADJUST fields

    DO i =  1,nx-1
      DO k =  1,nz-1
        DO j =  1,ny-1

          qv(i,j,k) = qv(i,j,k) * cqvfac
          qscalar(i,j,k,P_QC) = qscalar(i,j,k,P_QC) * cqcfac
          qscalar(i,j,k,P_QR) = qscalar(i,j,k,P_QR) * cqrfac
          qscalar(i,j,k,P_QI) = qscalar(i,j,k,P_QI) * cqifac
          qscalar(i,j,k,P_QS) = qscalar(i,j,k,P_QS) * cqsfac
          qscalar(i,j,k,P_QH) = qscalar(i,j,k,P_QH) * cqhfac

        END DO
      END DO
    END DO

  END IF

  RETURN
END SUBROUTINE lfo_ice_drive
!
!--------------------------------------------------------------------------
!
!
!  3-ICE MICROPHYSICS
!
!  VERSION: 1.3 with Rates (2/10/05)
!
!  LIN-FARELY-ORVILLE-like "Simple Ice and Liquid Microphysics Scheme"
!   based on Gilmore et al. (2004) Mon. Wea. Rev.
!
!   Copyright (C) <2004>  <Jerry Straka and Matthew Gilmore>
!
!   This library is free software; you can redistribute it and/or
!   modify it under the terms of the GNU Lesser General Public
!   License as published by the Free Software Foundation; either
!   version 2.1 of the License, or (at your option) any later version.
!
!   This library is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!   Lesser General Public License for more details.
!
!   If you find this code useful and publish results using it, please reference:
!
!  Gilmore M. S., J. M. Straka, and E. N. Rasmussen,
!        Monthly Weather Review: Vol. 132, No. 8, pp. 1897-1916.
!
!--------------------------------------------------------------------------

SUBROUTINE lfo_ice(nx,ny,nz,nor,dtp,dzc,ptbar,ptprt,qv,qscalar,  &
                   qvbar,rhobar,pprt,pbar)

!--------------------------------------------------------------------------
! RELEASE NOTES
!
! MSG - 2/10/05 Ported the SAM gather-scatter version to NCOMMAS
!    This is most similar to that actually used in the published
!    Gilmore et al. manuscripts. It also runs faster than the
!    column-based version used in NCOMMAS.
!
!    This version also provides rate output (with labels consistent
!    with Gilmore et al) for post-analysis.
!
!    Known omissions/inconsistencies with this version:
!       1) min q criteria only applied to some processes
!       2) internal rate names and signs differ compared to Gilmore et al. (2004)
!          however, rate output/labels are consistent with manuscripts
!       3) pre-micro-q fallout applied to post-micro-q values (faster
!          but not as accurate)
!       4) time splitting for fallout omitted (might blow up for fine dz)
!
!   Users may want to fix 3 and 4 before porting to their own models.
!   Users will also want to compile with an extend-line option since
!   goes past 72 columns.
!
!
!--------------------------------------------------------------------------
!
  IMPLICIT NONE

  INCLUDE 'globcst.inc'


  INTEGER, INTENT(IN) :: nx,ny,nz,nor
  REAL,    INTENT(IN) :: dtp
  REAL,    INTENT(IN) :: dzc(nz)    ! MSG weren't defined before

  REAL :: ptbar(-nor+1:nx+nor,-nor+1:ny+nor,-nor+1:nz+nor)
                                ! MSG ptbar (base state potential temp)
  REAL :: ptprt(-nor+1:nx+nor,-nor+1:ny+nor,-nor+1:nz+nor)
                                ! MSG ptprt (perturbation potential temp)
  REAL :: qv(-nor+1:nx+nor,-nor+1:ny+nor,-nor+1:nz+nor)
                                ! MSG qv    (vapor mixing ratio)
  REAL :: qvbar(nx,ny,nz)       ! MSG qvbar (base state vapor mixing ratio)

  REAL, TARGET :: qscalar(-nor+1:nx+nor,-nor+1:ny+nor,-nor+1:nz+nor,nscalar)

  REAL :: pprt(-nor+1:nx+nor,-nor+1:ny+nor,-nor+1:nz+nor)
                                ! MSG pprt  (perturbation pressure)
  REAL :: pbar(nx,ny,nz)        ! MSG pbar  (base state pressure)
  REAL :: rhobar(nx,ny,nz)      ! MSG rhobar (base state density = rhostr*j3inv)

  REAL :: pn(-nor+1:ny+nor,-nor+1:nz+nor,-nor+1:nx+nor)    !MSG perturb pressure
  REAL :: pb(nz), db(nz)                                   !MSG base state pressure, density

!-----------------------------------------------------------------------
!
!  general declarations
!
!-----------------------------------------------------------------------

  INTEGER, PARAMETER :: istag=1, jstag=1, kstag=1   !MSG
  INTEGER, PARAMETER :: lt=1, lv=2, lc=3, lr=4, li=5, ls=6, lh=7
                             !MSG changed index numbers (no longer used) ~NS
!  INTEGER, PARAMETER :: nxm=260, nym=260, nzm=68

!
!  Uncomment next 6 lines if not using "param.h" include file.
!
  REAL,   PARAMETER :: g=9.8, cp=1004., cv=717., rd=287.04, rw=461.5, rcp=rd/cp
  REAL,   PARAMETER :: cwdiap=20.e-6, qcmincwrn = 2.0E-3, cwdisp = 0.15
  INTEGER,PARAMETER :: autoconversion = 0
!  REAL,   PARAMETER :: rho_qr = 1000., cnor = 8.0E6 ! rain params
!  REAL,   PARAMETER :: rho_qs =  100., cnos = 3.0E6 ! snow params
!  REAL,   PARAMETER :: rho_qh =  900., cnoh = 4.0E4 ! hail params

  REAL,   PARAMETER :: rho_qr = 1000.
  REAL :: cnor
  REAL :: rho_qs, cnos
  REAL :: rho_qh, cnoh

  INTEGER :: ix,jy,kz
!
!  declarations microphyscs and for gather/scatter
!
  INTEGER, PARAMETER ::  ngs=2000
  INTEGER :: jgs,mgs,numgs,inumgs,nxmpb,nzmpb,nxz      !MSG added innumgs
  INTEGER :: ngscnt,igs(ngs),kgs(ngs)

! Air and particle temperatures (equivalent to each other in LFO version)

  REAL :: temp(nx,ny,ngs)
  REAL :: temg(ngs),temcg(ngs),theta(ngs),thetap(ngs),theta0(ngs)

! Air pressure, density, & Exner function

  REAL :: pbz(ngs),pres(ngs),presp(ngs),pres0(ngs),poo
  REAL :: rho0x(ngs),dnz(ngs),pi0(ngs),piz(ngs), rho00, dnz00

! Accretions

  REAL :: qiacr(ngs)
  REAL :: qracw(ngs), qraci(ngs), qracs(ngs)
  REAL :: qsacw(ngs), qsaci(ngs),            qsacr(ngs)
  REAL :: qhacw(ngs), qhaci(ngs), qhacs(ngs),qhacr(ngs)
  REAL :: eri(ngs),erw(ngs),ers(ngs)
  REAL :: esw(ngs),esi(ngs)
  REAL :: ehw(ngs),ehr(ngs),ehi(ngs),ehs(ngs)
  REAL :: xracwi(ngs),xiacr(ngs)
  REAL :: xsacwi(ngs),xcnos(ngs)
  REAL :: xhacwi(ngs),xcnoh(ngs)
  REAL :: xhacx(ngs)

! Biggs Freezing

  REAL :: qrfrz(ngs), xrfrz(ngs), arz,brz

! Bergeron Process

  REAL :: qsfw(ngs),qsfi(ngs), eic(ngs), bsfw
  REAL :: cs10(32),cs11(32),cs9
  REAL :: bfa1(32),bfa2(32),cbtim(32)
  REAL :: cmn, cmi40, cmi50, ri50, vti50, a, cm50a, cm40b, cm50b
  INTEGER :: ib
!   Koenig (1971) was corrected at T=-14 C (fifteenth index of bfa1) from 0.1725e-6 to 0.1725e-4 after a reviewer noted a typo.
  DATA bfa1/0., 0.7939E-7, 0.7841E-6, 0.3369E-5, 0.4336E-5, 0.5285E-5, 0.3728E-5, 0.1852E-5, 0.2991E-6, 0.4248E-6, &
      0.7434E-6, 0.1812E-5, 0.4394E-5, 0.9145E-5, 0.1725E-4, 0.3348E-4, 0.1725E-4, 0.9175E-5, 0.4412E-5, 0.2252E-5, &
      0.9115E-6, 0.4876E-6, 0.3473E-6, 0.4758E-6, 0.6306E-6, 0.8573E-6, 0.7868E-6, 0.7192E-6, 0.6513E-6, 0.5956E-6, &
      0.5333E-6, 0.4834E-6/
  DATA bfa2/0., 0.4006,    0.4831,    0.5320,    0.5307,    0.5319,    0.5249,    0.4888,    0.3894,    0.4047, &
      0.4318,    0.4771,    0.5183,    0.5463,    0.5651,    0.5813,    0.5655,    0.5478,    0.5203,    0.4906, &
      0.4447,    0.4126,    0.3960,    0.4149,    0.4320,    0.4506,    0.4483,    0.4460,    0.4433,    0.4413, &
      0.4382,    0.4361/

! Conversions

  REAL :: qhcns(ngs), qrcnw(ngs), qscni(ngs), qdiff,argrcnw
  REAL :: ehscnv(ngs), esicnv(ngs)

! Evaporation/Deposition/Sublimation

  REAL :: qhdsv(ngs),qhdpv(ngs),qhsbv(ngs)
  REAL :: qsdsv(ngs),qsdpv(ngs),qssbv(ngs)
  REAL :: qrcev(ngs)
  REAL :: xav(ngs),xbv(ngs),xas(ngs),xbs(ngs)
  REAL :: xrcev1(ngs),xrcev2(ngs)
  REAL :: xce,xrv,xds
  REAL :: xxdsv(ngs),xxcev(ngs)

! Initiation of cloud ice

  REAL :: dqisdt(ngs),qiint(ngs),xiint(ngs),cnnt

! Melting/Wet growth of hail, Melting snow

  REAL :: qhmlr(ngs), qsmlr(ngs)
  REAL :: xvth3(ngs)
  REAL :: xmlt1(ngs),xmlt2(ngs),xmlt3(ngs)
  REAL :: xhmlt1(ngs),xhmlt2(ngs)
  REAL :: xsmlt1(ngs),xsmlt2(ngs)
  REAL :: xsv, xhv, xhsw

! Wet/dry growth, shedding

  REAL :: qhdry(ngs),qhwet(ngs)
  REAL :: xhwet1(ngs),xhwet2(ngs), xcwt,xwt1
  REAL :: qhacip(ngs),qhacsp(ngs),qhshr(ngs)

! Water Budgets

  REAL :: ptotal(ngs),ptotsat(ngs)
  REAL :: pqcwi(ngs),pqcii(ngs),pqrwi(ngs)
  REAL :: pqswi(ngs),pqhwi(ngs),pqwvi(ngs)
  REAL :: pqcwd(ngs),pqcid(ngs),pqrwd(ngs)
  REAL :: pqswd(ngs),pqhwd(ngs),pqwvd(ngs)
  INTEGER :: il2(ngs),il3(ngs),il5(ngs)

! Latent Heating Computation

  REAL :: psub(ngs),pvap(ngs),pfrz(ngs),ptem(ngs)

! Maximum Depletion Tendencies

  REAL :: qc5dt,qi5dt,qr5dt,qs5dt,qh5dt

! Functions
  REAL, EXTERNAL :: gamma_lfo

! Flags
  INTEGER :: imake = 0               !prints intercept/density to output file
  INTEGER :: ndebug= 0, nrates = 1   !prints debug stuff to output file, prints rates TO out FILE

! Species threshold,  mass,    mass constraints

  REAL :: qcmin,qhmin,qimin,qrmin,qsmin
  REAL :: qccrit, qscrit, qicrit

! Fallout, fall velocity parameters/vars

  REAL :: dtz1
  REAL :: cwflx(ngs), cflux(nx,nz)
  REAL :: piflx(ngs), pflux(nx,nz)
  REAL :: rwflx(ngs), rflux(nx,nz)
  REAL :: swflx(ngs), sflux(nx,nz)
  REAL :: hwflx(ngs), hflux(nx,nz)
  REAL :: xvtr(ngs),xvts(ngs),xvth1(ngs),xvth2(ngs)

!  Distribution parameters, Fallout, Mean diameters, mass, and mixing ratios

  REAL :: xrslop(ngs),xsslop(ngs),xhslop(ngs)
!
  REAL :: xcnor(ngs)
!
  REAL :: qwv(ngs),qcw(ngs),qci(ngs),qrw(ngs),qsw(ngs),qhw(ngs)
  REAL :: ccw(ngs),cci(ngs),crw(ngs),csw(ngs),chw(ngs)
!
  REAL :: vtwbar(ngs),cwmas(ngs),cwdn(ngs)
  REAL :: vtibar(ngs),cimas(ngs),cidn(ngs)
  REAL :: vtrbar(ngs),           rwdn(ngs)
  REAL :: vtsbar(ngs),           swdn(ngs)
  REAL :: vthbar(ngs),           hwdn(ngs)
!
  REAL :: cwdia(ngs),cwdia2(ngs)
  REAL :: cidia(ngs),cidia2(ngs)
  REAL :: rwdia(ngs),rwdia2(ngs)
  REAL :: swdia(ngs),swdia2(ngs)
  REAL :: hwdia(ngs),hwdia2(ngs)

! Saturation Adjustment

  REAL :: qvap(ngs)
  REAL :: dqvcnd(ngs),dqwv(ngs),dqcw(ngs),dqci(ngs)
  REAL :: cqv1(ngs),cqv2(ngs)
  REAL :: gamss,denom1,denom2
  REAL :: gamw(ngs),gams(ngs)
  REAL :: qwvp(ngs), qv0(ngs)
  REAL :: fraci(ngs),fracl(ngs)
  REAL :: cdw, cdi
  INTEGER :: itertd
  REAL :: qwfzi(ngs),qimlw(ngs)
  REAL :: qidep(ngs),qisub(ngs)
  REAL :: qcevp(ngs), qccnd(ngs)
  REAL :: qcevpcnd(ngs), qisubdep(ngs)

! Saturation lookup table, vapor pressures, ratios

  INTEGER :: nqsat, ltemq, l
  PARAMETER(nqsat=20001)
  REAL :: tabqvs(nqsat),tabqis(nqsat)
  REAL :: fqsat, temq
  PARAMETER (fqsat=0.01)
  REAL :: cai, caw
  REAL :: qvs(ngs),qis(ngs),qss(ngs),pqs(ngs)
  REAL :: qss0(ngs)
  REAL :: tsqr(ngs),ssi(ngs),ssw(ngs)

! Misc constants

  REAL :: advisc0,advisc1
  REAL :: elv(ngs),elf(ngs),els(ngs)
  REAL :: ar,br,bta1,cs,ds
  REAL :: cbi,cbw,cpi,cnit,c1f3,cwc1(ngs)
  REAL :: gf6,gf5,gf4,gf3,gf2,gf1,gf4ds,gf4p5,gf4br
  REAL :: pi,pid4
  REAL :: dragh
  REAL :: tfr,thnuc,tfrcbw,tfrcbi,tka0
  REAL :: gamma_2, gamma_275, gamma_5ds, gamma_5br

! Misc Variables

  REAL :: advisc(ngs),schm(ngs)
  REAL :: wvdf(ngs),akvisc(ngs),ci(ngs),cw(ngs),tka(ngs)
  REAL :: cc3(ngs),cc4(ngs),cc5(ngs)

! Initiation

  REAL :: cwmasn,cwmasx
  REAL :: cimasn,cimasx

!  Rate Output (Domain-total g/m^3)

  REAL :: hfrz, hdep, hcnd, cevap, cmelt, csub          !MSG added total cooling/heating vars
  REAL :: tqva,  tqia,  tqca
  REAL :: tqvb,  tqib,  tqcb
  REAL :: tqvap, tqsap, tqiap, tqrap, tqcap, tqhap
  REAL :: tqvbp, tqsbp, tqibp, tqrbp, tqcbp, tqhbp
  REAL :: tvsum, tssum, tisum, trsum, tcsum, thsum, tqc,tqi,tqv,tsumall
  REAL :: suma,  sumb,  psum
  REAL :: trqsacw, trqhacr, trqhshr, trqhmlr, trqsmlr, tiqiint, tiqidep, trqhacs
  REAL :: trqhacw, trqhaci, tvqhsbv, tvqcevp, tvqssbv, tvqrcev, trqrcnw, trqracw
  REAL :: tcqcmli, tvqisub, tcqccnd, tiqifzc, thqhacs, thqsacr, thqhaci, thqhacr
  REAL :: thqhacw, thqhdpv, thqrfrz, thqiacr, thqracs, thqraci, tsqsacw, tsqsacr
  REAL :: tsqscni,  tsqsfi,  tsqsfw, tsqsdpv, thqhcns, tsqsaci, tsqraci, tsqiacr, thqhwet

  REAL, POINTER :: qc(:,:,:)
  REAL, POINTER :: qr(:,:,:)
  REAL, POINTER :: qi(:,:,:)
  REAL, POINTER :: qs(:,:,:)
  REAL, POINTER :: qh(:,:,:)
  INTEGER :: dsdpref                  ! preference for the newly added parameters below
  REAL :: n0rain   ! Intercept parameter for rainwater DSD
  REAL :: n0snow   ! Intercept parameter for snow DSD
  REAL :: n0hail   ! Intercept parameter for hail DSD
  REAL :: rhosnow  ! Snow density ( kg/m**3 )
  REAL :: rhohail  ! Hail density ( kg/m**3 )

  COMMON /DSD_paras/ dsdpref,n0rain,n0snow,n0hail,rhosnow,rhohail

!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code below ...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  qc => qscalar(:,:,:,P_QC)
  qr => qscalar(:,:,:,P_QR)
  qi => qscalar(:,:,:,P_QI)
  qs => qscalar(:,:,:,P_QS)
  qh => qscalar(:,:,:,P_QH)

!
!
!  read in constants from 'inmicro.jmslfo'
!  [deleted this part - MSG]
!
! modified by mtong
!
  cnor   = 8.0E6
  rho_qs = 100.
  cnos   = 3.0E6
  rho_qh = 900.
  cnoh   = 4.0E4

!  write(*,*)'dsdpref', dsdpref
  IF(dsdpref == 1)THEN
    cnor   = n0rain
    rho_qs = rhosnow
    cnos   = n0snow
    rho_qh = rhohail
    cnoh   = n0hail
    !write(*,*)'cnor', cnor, 'rho_qs', rho_qs, 'cnos', cnos, 'rho_qh', rho_qh, 'cnoh', cnoh
  ENDIF

  IF ( ndebug > 1 ) WRITE(6,*) 'just entered micro...'

!
!  array size tests
!
!  IF ( nx > nxm )  THEN
!    WRITE(6,*) 'STOP, NXM= ',nxm,' must be => NX= ',nx
!    STOP
!  END IF
!
!  IF ( nz > nzm )  THEN
!    WRITE(6,*) 'STOP, NZM= ',nzm,' must be => NZ= ',nz
!    STOP
!  END IF

!
!  ZERO
!
!
!  totals for source / sink terms
!
!
!  vapor
!
  tvqrcev = 0.0
  tvqssbv = 0.0
  tvqhsbv = 0.0
  tvqcevp = 0.0
  tvqisub = 0.0
!
!  cloud water
!
  tcqccnd = 0.0
  tcqcmli = 0.0
!
!  rain
!
  trqrcnw = 0.0
  trqracw = 0.0
  trqhmlr = 0.0
  trqsmlr = 0.0
  trqhshr = 0.0
  trqsacw = 0.0
  trqhacr = 0.0
  trqhacw = 0.0
  trqhaci = 0.0
  trqhacs = 0.0

!
!  cloud ice
!
  tiqiint = 0.0
  tiqidep = 0.0
  tiqifzc = 0.0
!
!  snow
!
  tsqsfi  = 0.0
  tsqsfw  = 0.0
  tsqscni = 0.0
  tsqsacw = 0.0
  tsqsacr = 0.0
  tsqraci = 0.0
  tsqiacr = 0.0
  tsqsaci = 0.0
  tsqsdpv = 0.0
!
!  hail
!
  thqhcns = 0.0
  thqhacr = 0.0
  thqhacw = 0.0
  thqhaci = 0.0
  thqhacs = 0.0
  thqsacr = 0.0
  thqracs = 0.0
  thqraci = 0.0
  thqiacr = 0.0
  thqhdpv = 0.0
  thqrfrz = 0.0
  thqhwet = 0.0
!
!  total heating and cooling rates
!
  hfrz    = 0.0
  hdep    = 0.0
  hcnd    = 0.0
  cevap   = 0.0
  cmelt   = 0.0
  csub    = 0.0
!
!
!  various rate budgets
!
  tqvap     = 0.0
  tqcap     = 0.0
  tqiap     = 0.0
  tqrap     = 0.0
  tqsap     = 0.0
  tqhap     = 0.0
!
  tqvbp     = 0.0
  tqcbp     = 0.0
  tqibp     = 0.0
  tqrbp     = 0.0
  tqsbp     = 0.0
  tqhbp     = 0.0
!
  tqva     = 0.0
  tqca     = 0.0
  tqia     = 0.0
!
  tqvb     = 0.0
  tqcb     = 0.0
  tqib     = 0.0
!
  suma      = 0.0
  psum      = 0.0
!
  sumb      = 0.0
!
  tvsum    = 0.0
  tcsum    = 0.0
  tisum    = 0.0
  trsum    = 0.0
  tssum    = 0.0
  thsum    = 0.0
!
  tsumall   = 0.0
  tqv       = 0.0
  tqc       = 0.0
  tqi       = 0.0
!
!
!  end of totals
!
!
!  other constants
!


  IF ( ndebug > 1 ) PRINT*,'dbg = 0a'
!
!  constants
!
  poo = 1.0E+05
  ar = 841.99666
  br = 0.8
  bta1 = 0.6
  cnit = 1.0E-02
  dnz00 = 1.225
  rho00 = 1.225
  cs = 4.83607122
  ds = 0.25
  pi = 4.0*ATAN(1.0)
  pid4 = pi/4.0
  qccrit = 2.0E-03
  qscrit = 6.0E-04
  qicrit = 1.0E-03
  gf1 = gamma_lfo(1.0)
  gf2 = gamma_lfo(2.0)
  gf3 = gamma_lfo(3.0)
  gf4 = gamma_lfo(4.0)
  gf5 = gamma_lfo(5.0)
  gf6 = gamma_lfo(6.0)
  gf4br = gamma_lfo(4.0+br)
  gf4ds = gamma_lfo(4.0+ds)
  gf4p5 = gamma_lfo(4.0+0.5)
!  rho_qr = 1000.0
!  rho_qs = 100.0
!  rho_qh = 900.0
!  cnor = 8.0e+06
!  cnos = 3.0e+06
!  cnoh = 4.0e+04
  IF ( rho_qh < 600.0 ) THEN
    dragh = 1.00
  ELSE
    dragh = 0.60
  END IF

  IF ( imake == 0 ) THEN

    IF (ndebug > 1) THEN
      WRITE(6,*) '-----------------------------------------------------------------------'
      WRITE(6,*) '3-ICE MICROPHYSICAL CONSTANTS'
      101    FORMAT(1X,a,6(g10.5,2X))
      WRITE(6,101) 'CNOR/DENR:        ',cnor, rho_qr
      WRITE(6,101) 'CNOS/DENS:        ',cnos, rho_qs
      WRITE(6,101) 'CNOH/DENH/DRAGH:  ',cnoh, rho_qh,dragh
      WRITE(6,*)   'TIME STEP:        ', dtp
      IF( autoconversion == 0 ) THEN
        WRITE(6,*) 'Berry (1968) Critical qc g/g for autoconversion= ',qcmincwrn
      ELSE
        WRITE(6,*) 'Berry (1968) Critical qc diam for autoconversion= ',cwdiap*1E6,' microns'
      END IF
    END IF

    imake = 1
  END IF

  IF (ndebug > 1) WRITE(6,*) '-----------------------------------------------------------------------'

!  constants
!
  c1f3 = 1.0/3.0
!
!  general constants for microphysics
!
  brz = 100.0
  arz = 0.66
  cai = 21.87455
  caw = 17.2693882
  cbi = 7.66
  cbw = 35.86
  qcmin = 1.0E-09
  qimin = 1.0E-12
  qrmin = 1.0E-07
  qsmin = 1.0E-07
  qhmin = 1.0E-07

  tfr = 273.16
  thnuc = 233.16
  advisc0 = 1.832E-05
  advisc1 = 1.718E-05
  tka0 = 2.43E-02
  cpi = 1.0/cp
  tfrcbw = tfr - cbw
  tfrcbi = tfr - cbi
!
  IF ( ndebug > 1 ) PRINT*,'dbg = 0b'

  DO kz = 1,nz-kstag

!   IF (ieee_is_nan(piz(kz))) STOP 'A'
!   IF (ieee_is_nan(dnz(kz))) STOP 'B'
!   IF (ieee_is_nan(pbz(kz))) STOP 'C'
!   IF (ieee_is_nan(temp(kz))) STOP 'D'
    elv(kz) = 2500300.
    elf(kz) = 335717.
    els(kz) = elv(kz) + elf(kz)

  END DO
!
  DO l = 1,nqsat
    temq = 163.15 + (l-1)*fqsat
    tabqvs(l) = EXP(caw*(temq-273.15)/(temq-cbw))
    tabqis(l) = EXP(cai*(temq-273.15)/(temq-cbi))
!   IF (ieee_is_nan(tabqvs(l))) STOP 'L'
!   IF (ieee_is_nan(tabqis(l))) STOP 'M'
  END DO
!
!  cw constants in mks units
!
  cwmasn = 4.25E-15
  cwmasx = 5.25E-10
  DO kz = 1,nz-kstag
    cwdn(kz) = 1000.0
    cwc1(kz) = 6.0/(pi*cwdn(kz))
  END DO
!
!  ci constants in mks units
!
  cimasx = 3.23E-8
  cimasn = 4.25E-15
  DO kz = 1,nz-kstag
    cidn(kz) = 570.0
  END DO

!
!  constants for paramerization
!
  gamma_2 = 0.78*gamma_lfo(2.0)
  DO kz = 1,nz-kstag
    xcnor(kz) = cnor
    xcnos(kz) = cnos
    xcnoh(kz) = cnoh
    xhacx(kz) = pid4/gf4
    xracwi(kz) = pid4
    xsacwi(kz) = pid4
    xhacwi(kz) = pid4
    xiacr(kz) = pid4/gf4
    xvth1(kz) = (gf4p5/(6.0))
    cw(kz) = 4218.0

    xmlt3(kz) = -cw(kz)/elf(kz)
    xhmlt1(kz) = gamma_2
    xsmlt1(kz) = gamma_2
    xrcev1(kz) = gamma_2
    xiint(kz) = (elv(kz)**2)/(cp*rw)
  END DO

!
!  constants for bergeron process, note in cgs units
!
  cmn = 1.05E-15
  cmi40 = 2.4546E-07
  cmi50 = 4.8E-07
  ri50 = 5.0E-03
  vti50 = 100.0
  bsfw = (ri50**2)*vti50*pi
  cbtim(1) = 1.0E+30
  cs10(1)  = 0.0
  cs11(1)  = 0.0
  DO ib = 2,32
    cm50a = cmi50**bfa2(ib)
    a     = 1.0-bfa2(ib)
    cm40b = cmi40**a
    cm50b = cmi50**a
    cbtim(ib) = (cm50b-cm40b)/(bfa1(ib)*a)
    cs10(ib) = bfa1(ib)*cm50a
    cs11(ib) = dtp/(cbtim(ib)*cmi50)
!    write(6,*) 'ib, cbtim(ib)=',ib, cbtim(ib)
  END DO

  IF (ndebug > 1 ) PRINT*,'dbg = 1'

  IF (ndebug > 1 ) PRINT*,'dbg = 2'
!
!  start jy loop  (for doing XZ slabs)
!
!
  gamma_275 = 0.308*gamma_lfo(2.75)
  gamma_5ds = 0.308*(gamma_lfo((5.0+ds)/2.0))
  gamma_5br = 0.308*(gamma_lfo((5.0+br)/2.0))
  DO jy = 1,ny-jstag
!
!  VERY IMPORTANT:  SET jgs
!
    jgs = jy
!
!  zero precip flux arrays
!
    IF (ndebug > 1 ) PRINT*,'dbg = 3'
    DO kz = 1,nz-kstag
      DO ix = 1,nx-istag
        hflux(ix,kz) = 0.0
        cflux(ix,kz) = 0.0
        pflux(ix,kz) = 0.0
        rflux(ix,kz) = 0.0
        sflux(ix,kz) = 0.0
      END DO
    END DO
!
!..gather microphysics
!

    IF (ndebug > 1 ) PRINT*,'dbg = 4'
    nxmpb = 1
    nzmpb = 1
    nxz = nx*nz
    numgs = nxz/ngs + 1
    DO inumgs = 1,numgs
      ngscnt = 0
      DO kz = nzmpb,nz-kstag-1
        DO ix = nxmpb,nx-istag
!--------The variables in this segment have been moved to within the ix and jy loops----(NS 3/30, 4/08, 4/15)-----!
          xrfrz(kz) = (20.0)*(pi**2)*brz*(cwdn(kz)/rhobar(ix,jy,kz))
          xvtr(kz) = (ar*gf4br*(dnz00/rhobar(ix,jy,kz))**(0.5))/(6.0)
          xvts(kz) = (cs*gf4ds*(dnz00/rhobar(ix,jy,kz))**(0.5))/(6.0)
          xvth2(kz) = (4.0*g/(3.0*dragh*rhobar(ix,jy,kz)))**(0.5)
          xvth3(kz) = (4.0*g/(3.0*dragh*rhobar(ix,jy,kz)))**(0.25)
          xrslop(kz) = (rhobar(ix,jy,kz)/(pi))**(0.25)
          xsslop(kz) = (rhobar(ix,jy,kz)/(pi))**(0.25)
          xhslop(kz) = (rhobar(ix,jy,kz)/(pi))**(0.25)
          xxcev(kz) = (2.0*pi/rhobar(ix,jy,kz))
          xxdsv(kz) =   2.0*pi/rhobar(ix,jy,kz)
          xmlt1(kz) = -2.0*pi/(elf(kz)*rhobar(ix,jy,kz))
          xhwet1(kz) = 2.0*pi/rhobar(ix,jy,kz)

          piz(kz) = (pbar(ix,jy,kz)/poo)**rcp
          pbz(kz) = pbar(ix,jy,kz)
          gamw(kz) = elv(kz)*cpi/piz(kz)
          gams(kz) = els(kz)*cpi/piz(kz)
          cqv1(kz) = 4098.0258*piz(kz)*gamw(kz)
          cqv2(kz) = 5807.6953*piz(kz)*gams(kz)
          cc3(kz) = cpi*elf(kz)/piz(kz)
          cc4(kz) = cpi*elv(kz)/piz(kz)
          cc5(kz) = cpi*els(kz)/piz(kz)

          temp(ix,jy,kz) = piz(kz)*ptbar(ix,jy,kz)  !replaced ab(kz,LT) with ptbar(ix,jy,kz)
          wvdf(kz) = (2.11E-05)*((temp(ix,jy,kz)/tfr)**1.94)*(101325.0/(pbz(kz)))
          advisc(kz) = advisc0*(416.16/(temp(ix,jy,kz)+120.0))*(temp(ix,jy,kz)/296.0)**(1.5)
          akvisc(kz) = advisc(kz)/rhobar(ix,jy,kz)
          tka(kz) = tka0*advisc(kz)/advisc1

          xav(kz) = (rhobar(ix,jy,kz)**2)/(tka(kz)*rw)
          xas(kz) = (els(kz)**2)/(tka(kz)*rw)

          schm(kz) = (akvisc(kz)/wvdf(kz))
          ci(kz) = (2.118636 + 0.007371*(temp(ix,jy,kz)-tfr))*(1.0E+03)
          xbv(kz) = (1.0/(rhobar(ix,jy,kz)*wvdf(kz)))
          xbs(kz) = (1.0/(rhobar(ix,jy,kz)*wvdf(kz)))
          xmlt2(kz) = wvdf(kz)*elv(kz)*rhobar(ix,jy,kz)
          xhmlt2(kz) = gamma_275*(schm(kz)**(1./3.))*(akvisc(kz)**(-0.5))
          xsmlt2(kz) = gamma_5ds*(schm(kz)**(1./3.))        &
               *(akvisc(kz)**(-0.5))*(cs**(0.5))*((dnz00/rhobar(ix,jy,kz))**(0.25))
          xrcev2(kz) = gamma_5br*(schm(kz)**(1./3.))        &
               *(akvisc(kz)**(-0.5))*(ar**(0.5))*((dnz00/rhobar(ix,jy,kz))**(0.25))
          xhwet2(kz) = rhobar(ix,jy,kz)*elv(kz)*wvdf(kz)
!-----------------------------------------------------------------------------------------------------!
          pqs(kz) = 380.0/(pprt(ix,jy,kz)+pbar(ix,jy,kz))
          qwv(kz) = MAX(qv(ix,jy,kz)  ,0.0)
          qcw(kz) = MAX(qc(ix,jy,kz)  ,0.0)
          qci(kz) = MAX(qi(ix,jy,kz)  ,0.0)
          theta(kz) = ptprt(ix,jy,kz) + ptbar(ix,jy,kz)
          temg(kz) = theta(kz)*( (pprt(ix,jy,kz)+pbar(ix,jy,kz)) / poo ) ** rcp
          temcg(kz) = temg(kz) - tfr
          ltemq = nint((temg(kz)-163.15)/fqsat+1.5)
          ltemq = MIN(MAX(ltemq,1),nqsat)
          qvs(kz) = pqs(kz)*tabqvs(ltemq)
          qis(kz) = pqs(kz)*tabqis(ltemq)

          IF ( temg(kz) < tfr ) THEN
            IF( qcw(kz) >= 0.0 .AND. qci(kz) == 0.0 ) qss(kz) = qvs(kz)
            IF( qcw(kz) == 0.0 .AND. qci(kz) > 0.0)  qss(kz) = qis(kz)
            IF( qcw(kz) > 0.0 .AND. qci(kz) > 0.0)  qss(kz) = (qcw(kz)*qvs(kz) + qci(kz)*qis(kz)) /(qcw(kz) + qci(kz))
          ELSE
            qss(kz) = qvs(kz)
          END IF
!
          IF ( qv(ix,jy,kz) > qss(kz) .OR.  qc(ix,jy,kz) > qcmin  .OR. &
                 qi(ix,jy,kz) > qimin .OR.                           &
                 qr(ix,jy,kz) > qrmin .OR.                           &
                 qs(ix,jy,kz) > qsmin .OR.                           &
                 qh(ix,jy,kz) > qhmin ) THEN
            ngscnt = ngscnt + 1
            igs(ngscnt) = ix
            kgs(ngscnt) = kz
            IF ( ngscnt == ngs ) GO TO 1100
          END IF
        END DO     !MSG - i loop
        nxmpb = 1
      END DO     !MSG - k loop
      1100 CONTINUE
      IF ( ngscnt == 0 ) GO TO 9998
      IF ( ndebug > 1 ) PRINT*,'dbg = 5'
!
!  define temporaries to be used in calculations
!
!write(*,*)'ngscnt', ngscnt
      DO mgs = 1,ngscnt
        theta0(mgs) = ptbar(igs(mgs),jy,kgs(mgs))
        thetap(mgs) = ptprt(igs(mgs),jy,kgs(mgs))
        theta(mgs) = ptbar(igs(mgs),jy,kgs(mgs)) + ptprt(igs(mgs),jy,kgs(mgs))
        pres0(mgs) = pbar(igs(mgs),jy,kgs(mgs))
        presp(mgs) = pprt(igs(mgs),jy,kgs(mgs))
        pres(mgs) = presp(mgs) + pres0(mgs)
        rho0x(mgs) = rhobar(igs(mgs),jy,kgs(mgs))
        pi0(mgs) = piz(kgs(mgs))
        temg(mgs) = theta(mgs)*( pres(mgs) / poo ) ** rcp

        temcg(mgs) = temg(mgs) - tfr
        qss0(mgs) = (380.0)/(pres(mgs))
        pqs(mgs) = (380.0)/(pres(mgs))
        ltemq = nint((temg(mgs)-163.15)/fqsat+1.5)
        ltemq = MIN(MAX(ltemq,1),nqsat)
        qvs(mgs) = pqs(mgs)*tabqvs(ltemq)
        qis(mgs) = pqs(mgs)*tabqis(ltemq)
        qv0(mgs) = qvbar(igs(mgs),jy,kgs(mgs))
!write(*,*)igs(mgs), jy, kgs(mgs)
!write(*,*)qv(igs(mgs),jy,kgs(mgs)), qvbar(igs(mgs),jy,kgs(mgs))
        qwvp(mgs) = qv(igs(mgs),jy,kgs(mgs)) - qvbar(igs(mgs),jy,kgs(mgs))
        qwv(mgs) = qv(igs(mgs),jy,kgs(mgs))
        qcw(mgs) = MAX(qc(igs(mgs),jy,kgs(mgs)), 0.0)
        qci(mgs) = MAX(qi(igs(mgs),jy,kgs(mgs)), 0.0)
        qrw(mgs) = MAX(qr(igs(mgs),jy,kgs(mgs)), 0.0)
        qsw(mgs) = MAX(qs(igs(mgs),jy,kgs(mgs)), 0.0)
        qhw(mgs) = MAX(qh(igs(mgs),jy,kgs(mgs)), 0.0)
        il2(mgs) = 0
        il3(mgs) = 0
        il5(mgs) = 0
        IF ( temg(mgs) < tfr ) THEN
          il5(mgs) = 1
          IF ( qrw(mgs) < 1.0E-04 .AND. qsw(mgs) < 1.0E-04 ) il2(mgs) = 1
          IF ( qrw(mgs) < 1.0E-04 ) il3(mgs) = 1
        END IF
      END DO

      IF (ndebug > 1 ) PRINT*,'dbg = 6'
!
! cloud water variables
!
      DO mgs = 1,ngscnt
        ccw(mgs) = 1.e9    !LFO default (Western Plains)
        !ccw(mgs) = .6e9    !Central plains CCN value
        !ccw(mgs) = .3e9    !Maritime CCN value
        cwmas(mgs) = MIN( MAX(qcw(mgs)*rho0x(mgs)/ccw(mgs),cwmasn),cwmasx )
        cwdia(mgs) = (cwmas(mgs)*cwc1(kgs(mgs)))**c1f3
        cwdia2(mgs) = cwdia(mgs)**2
        vtwbar(mgs) = (ar*(cwdia(mgs)**br))*(rho00/rho0x(mgs))**0.5
      END DO
!
! cloud ice variables
!
      DO mgs = 1,ngscnt

        cimasx     = 3.23E-8
        cci(mgs)   = MAX(MIN(cnit*EXP(-temcg(mgs)*bta1),1.e+09),1.0)       !Fletcher's formula
        cimas(mgs) = MIN( MAX(qci(mgs)*rho0x(mgs)/cci(mgs),cimasn),cimasx )
        IF ( temcg(mgs) > 0 ) THEN
          cidia(mgs) = 0.0
        ELSE
          cidia(mgs) = 16.7*(cimas(mgs)**(0.5))
          cidia(mgs) = MAX(cidia(mgs), 1.e-5)
          cidia(mgs) = MIN(cidia(mgs), 3.e-3)
        END IF
        cidia2(mgs) = cidia(mgs)**2
        vtibar(mgs) = (cs*(cidia(mgs)**ds))*(rho00/rho0x(mgs))**0.5

      END DO
!
!  mp-distribution information for rain, snow agg's, and graupel/hail
!
!  definitions for marshall palmer distribution variables
!  (rain, snow, hail) when mixing ratio only is predicted
!
      IF (ndebug > 1 ) PRINT*,'dbg = 7a'

      DO mgs = 1,ngscnt

        rwdn(mgs) = rho_qr
        swdn(mgs) = rho_qs
        hwdn(mgs) = rho_qh
        rwdia(mgs) = 1.e-20
        swdia(mgs) = 1.e-20
        hwdia(mgs) = 1.e-20

        IF ( qrw(mgs) > 1.0E-10 ) rwdia(mgs) =     &
          xrslop(kgs(mgs))*(qrw(mgs)/(rwdn(mgs)*xcnor(kgs(mgs))))**(0.25)
        IF ( qsw(mgs) > 1.0E-10 ) swdia(mgs) =     &
          xsslop(kgs(mgs))*(qsw(mgs)/(swdn(mgs)*xcnos(kgs(mgs))))**(0.25)
        IF ( qhw(mgs) > 1.0E-10 ) hwdia(mgs) =     &
          xhslop(kgs(mgs))*(qhw(mgs)/(hwdn(mgs)*xcnoh(kgs(mgs))))**(0.25)

        rwdia2(mgs) = rwdia(mgs)**2
        swdia2(mgs) = swdia(mgs)**2
        hwdia2(mgs) = hwdia(mgs)**2

        vtrbar(mgs) = xvtr(kgs(mgs))*(rwdia(mgs)**br)
        vtsbar(mgs) = xvts(kgs(mgs))*(swdia(mgs)**ds)
        vthbar(mgs) = xvth1(kgs(mgs))*xvth2(kgs(mgs))*((hwdn(mgs)*hwdia(mgs))**(0.5))

        crw(mgs) = xcnor(kgs(mgs))*rwdia(mgs)
        csw(mgs) = xcnos(kgs(mgs))*swdia(mgs)
        chw(mgs) = xcnoh(kgs(mgs))*hwdia(mgs)

      END DO

      DO mgs = 1,ngscnt
        !
        !  maximum depletion tendency by any one source
        !
        qc5dt = 0.20*qcw(mgs)/dtp
        qi5dt = 0.20*qci(mgs)/dtp
        qr5dt = 0.20*qrw(mgs)/dtp
        qs5dt = 0.20*qsw(mgs)/dtp
        qh5dt = 0.20*qhw(mgs)/dtp

        !
        !  collection efficiencies
        !
        eic(mgs) = 1.0
        eri(mgs) = 1.0
        erw(mgs) = 1.0
        esi(mgs) = EXP(0.025*MIN(temcg(mgs),0.0))
        esicnv(mgs) = esi(mgs)
        esw(mgs) = 1.0
        ers(mgs) = 1.0
        ehw(mgs) = 1.0
        ehr(mgs) = 1.0
        ehs(mgs) = EXP(0.09*MIN(temcg(mgs),0.0))
        ehscnv(mgs) = ehs(mgs)
        IF ( temcg(mgs) > 0.0 ) ehs(mgs) = 1.0
        ehi(mgs) = 0.1

        IF ( qcw(mgs) < qcmin .OR. qci(mgs) < qimin ) eic(mgs) = 0.0
        IF ( qrw(mgs) < qrmin .OR. qci(mgs) < qimin ) eri(mgs) = 0.0
        IF ( qrw(mgs) < qrmin .OR. qcw(mgs) < qcmin ) erw(mgs) = 0.0
        IF ( qsw(mgs) < qsmin .OR. qci(mgs) < qimin ) esi(mgs) = 0.0
        IF ( qsw(mgs) < qsmin .OR. qcw(mgs) < qcmin ) esw(mgs) = 0.0
        IF ( qsw(mgs) < qsmin .OR. qrw(mgs) < qrmin ) ers(mgs) = 0.0
        IF ( qhw(mgs) < qhmin .OR. qcw(mgs) < qcmin ) ehw(mgs) = 0.0
        IF ( qhw(mgs) < qhmin .OR. qrw(mgs) < qrmin ) ehr(mgs) = 0.0
        IF ( qhw(mgs) < qhmin .OR. qsw(mgs) < qsmin ) ehs(mgs) = 0.0
        IF ( qhw(mgs) < qhmin .OR. qci(mgs) < qimin ) ehi(mgs) = 0.0
!
!  accretions:
!    marshall-palmer size distribution collection
!    of constant size distribution
!   1)  sink for constant size distribution
!   2)  source for marshall-palmer size distribution
!
        qracw(mgs) =                                                    &
            MIN(erw(mgs)*qcw(mgs)*xracwi(kgs(mgs))*crw(mgs)*ABS(vtrbar(mgs)-vtwbar(mgs)) &
            *(  gf3*rwdia2(mgs) + 2.0*gf2*rwdia(mgs)*cwdia(mgs) + gf1*cwdia2(mgs) )      , qc5dt)
        qsacw(mgs) =                                                    &
            MIN(esw(mgs)*qcw(mgs)*xsacwi(kgs(mgs))*csw(mgs)*ABS(vtsbar(mgs)-vtwbar(mgs)) &
            *(  gf3*swdia2(mgs) + 2.0*gf2*swdia(mgs)*cwdia(mgs) + gf1*cwdia2(mgs) )      , qc5dt)
        qhacw(mgs) =                                                    &
            MIN(ehw(mgs)*qcw(mgs)*xhacwi(kgs(mgs))*chw(mgs)*ABS(vthbar(mgs)-vtwbar(mgs)) &
            *(  gf3*hwdia2(mgs) + 2.0*gf2*hwdia(mgs)*cwdia(mgs) + gf1*cwdia2(mgs) )      , qc5dt)
        qraci(mgs) =                                                    &
            MIN(eri(mgs)*qci(mgs)*xracwi(kgs(mgs))*crw(mgs)*ABS(vtrbar(mgs)-vtibar(mgs)) &
            *(  gf3*rwdia2(mgs) + 2.0*gf2*rwdia(mgs)*cidia(mgs) + gf1*cidia2(mgs) )      , qi5dt)
        qsaci(mgs) =                                                    &
            MIN(esi(mgs)*qci(mgs)*xsacwi(kgs(mgs))*csw(mgs)*ABS(vtsbar(mgs)-vtibar(mgs)) &
            *(  gf3*swdia2(mgs) + 2.0*gf2*swdia(mgs)*cidia(mgs) + gf1*cidia2(mgs) )      , qi5dt)
        qhaci(mgs) =                                                    &
            MIN(ehi(mgs)*qci(mgs)*xhacwi(kgs(mgs))*chw(mgs)*ABS(vthbar(mgs)-vtibar(mgs)) &
            *(  gf3*hwdia2(mgs) + 2.0*gf2*hwdia(mgs)*cidia(mgs) + gf1*cidia2(mgs) )      , qi5dt)
        qiacr(mgs) =                                                    &
            MIN(eri(mgs)*qrw(mgs)* xiacr(kgs(mgs))*cci(mgs)*ABS(vtrbar(mgs)-vtibar(mgs)) &
            *(  gf6*rwdia2(mgs) + 2.0*gf5*rwdia(mgs)*cidia(mgs) + gf4*cidia2(mgs) )      , qr5dt)
!
!  accretions:
!    marshall-palmer size distribution collecting marshall-palmer size
!    distribution
!
        qhacs(mgs) =                                                    &
            MIN( xhacx(kgs(mgs))*ABS(vthbar(mgs)-vtsbar(mgs))*ehs(mgs)*qsw(mgs)*chw(mgs) &
            *(  gf6*gf1*swdia2(mgs) + 2.0*gf5*gf2*swdia(mgs)*hwdia(mgs) + gf4*gf3*hwdia2(mgs) ) , qs5dt)
        qhacr(mgs) =                                                    &
            MIN( xhacx(kgs(mgs))*ABS(vthbar(mgs)-vtrbar(mgs))*ehr(mgs)*qrw(mgs)*chw(mgs) &
            *(  gf6*gf1*rwdia2(mgs) + 2.0*gf5*gf2*rwdia(mgs)*hwdia(mgs) + gf4*gf3*hwdia2(mgs) ) , qr5dt)
        qracs(mgs) =                                                    &
            MIN( xhacx(kgs(mgs))*ABS(vtrbar(mgs)-vtsbar(mgs))*ers(mgs)*qsw(mgs)*crw(mgs) &
            *(  gf6*gf1*swdia2(mgs) + 2.0*gf5*gf2*swdia(mgs)*rwdia(mgs) + gf4*gf3*rwdia2(mgs) ) , qs5dt)
        qsacr(mgs) =                                                    &
            MIN( xhacx(kgs(mgs))*ABS(vtrbar(mgs)-vtsbar(mgs))*ers(mgs)*qrw(mgs)*csw(mgs) &
            *(  gf6*gf1*rwdia2(mgs) + 2.0*gf5*gf2*rwdia(mgs)*swdia(mgs) + gf4*gf3*swdia2(mgs) ) , qr5dt)
!
!  bergeron process for snow
!
        ib = MIN(MAX(1,INT(-temcg(mgs))),32)
        cs9     = bsfw*rhobar(igs(mgs),jy,kgs(mgs))*(0.001)
        qsfw(mgs) = qci(mgs)*cs11(ib)*(cs10(ib) + eic(mgs)*cs9*qcw(mgs))
        qsfw(mgs) = MIN(qsfw(mgs),qc5dt)
        qsfi(mgs) = qci(mgs)/cbtim(ib)
        qsfi(mgs) = MIN(qsfi(mgs),qi5dt)
!
!  conversions
!
        qscni(mgs) = 0.001*esicnv(mgs)*MAX((qci(mgs)-qicrit),0.0)
        qscni(mgs) = MIN(qscni(mgs),qi5dt)
!
        qhcns(mgs) = 0.001*ehscnv(mgs)*MAX((qsw(mgs)-qscrit),0.0)
        qhcns(mgs) = MIN(qhcns(mgs),qs5dt)
!
        qrfrz(mgs) = MIN(xrfrz(kgs(mgs))*crw(mgs)*(rwdia(mgs)**6)*(EXP(MAX(-arz*temcg(mgs), 0.0))-1.0), qr5dt)

! Berry (1968) Autoconversion == 0 (critical qc) or ==1 (critical diameter from Ferrier 1994)

        IF( autoconversion == 0 ) THEN
          qdiff  = MAX((qcw(mgs)-qcmincwrn),0.)
        ELSE
          qccrit = (pi/6.)*((ccw(mgs)*cwdiap**3)*cwdn(kgs(mgs)))/rho0x(mgs)
          qdiff  = MAX((qcw(mgs)-qccrit),0.)
        END IF

        qrcnw(mgs) =  0.0
        IF ( qdiff > 0.0 ) THEN
          argrcnw = ((1.2E-4)+(1.596E-12)*ccw(mgs)*(1E-6)/(cwdisp*qdiff))
          qrcnw(mgs) = rho0x(mgs)*1E-3*(qdiff**2)/argrcnw
          qrcnw(mgs) = (MAX(qrcnw(mgs),0.0))
        END IF
        qrcnw(mgs) = MIN(qrcnw(mgs),qc5dt)
!
!  constants for hydrometeor-vapor interactions
!
        ssi(mgs) = qwv(mgs)/qis(mgs)
        ssw(mgs) = qwv(mgs)/qvs(mgs)
        tsqr(mgs) = temg(mgs)**2
!
!  melting of snow and hail
!
        xsv = (xsmlt1(kgs(mgs))*(swdia(mgs)) + xsmlt2(kgs(mgs))*(swdia(mgs)**((3.0+ds)/2.0)))
        xhv = (xhmlt1(kgs(mgs))*(hwdia(mgs)) + (hwdn(mgs)**(0.25))*xvth3(kgs(mgs))*xhmlt2(kgs(mgs))*(hwdia(mgs)**(1.75)))
        xhsw =(tka(kgs(mgs))*temcg(mgs) + xmlt2(kgs(mgs))*(qwv(mgs)-qss0(mgs)))
        qsmlr(mgs) = MIN( (xmlt1(kgs(mgs))*csw(mgs)*xsv*xhsw + temcg(mgs)*xmlt3(kgs(mgs))*(qsacr(mgs)+qsacw(mgs)) ) , 0.0 )
        qhmlr(mgs) = MIN( (xmlt1(kgs(mgs))*chw(mgs)*xhv*xhsw + temcg(mgs)*xmlt3(kgs(mgs))*(qhacr(mgs)+qhacw(mgs)) ) , 0.0 )
        qsmlr(mgs) = MAX( qsmlr(mgs), -qs5dt )
        qhmlr(mgs) = MAX( qhmlr(mgs), -qh5dt )
!
!  deposition/sublimation of snow and hail
!
        xds = xxdsv(kgs(mgs))*(ssi(mgs)-1.0)*(1.0/(xas(kgs(mgs))/tsqr(mgs)+xbs(kgs(mgs))/qis(mgs)))
        qsdsv(mgs) =   xds*csw(mgs)*xsv
        qhdsv(mgs) =   xds*chw(mgs)*xhv
        qhsbv(mgs) = MAX( MIN(qhdsv(mgs), 0.0), -qh5dt )
        qhdpv(mgs) = MAX( qhdsv(mgs), 0.0 )
        qssbv(mgs) = MAX( MIN(qsdsv(mgs), 0.0), -qs5dt )
        qsdpv(mgs) = MAX( qsdsv(mgs), 0.0 )

!
! SHEDDING CALCULATION
! New version by MSG closer to JMS original  - Last modified 4/6/03
!
!
!  compute dry growth rate of hail regardless of location
!
        qhdry(mgs) = qhacr(mgs) + qhacw(mgs) + qhaci(mgs) + qhacs(mgs)
!
!  compute wet growth rate of hail regardless of location
!
        qhacip(mgs)= qhaci(mgs)          !ehi=0 case
        qhacsp(mgs)= qhacs(mgs)          !ehs=0 case
        IF ( ehi(mgs) > 0.0 ) qhacip(mgs) = MIN(qhaci(mgs)/ehi(mgs),qi5dt)
        IF ( ehs(mgs) > 0.0 ) qhacsp(mgs) = MIN(qhacs(mgs)/ehs(mgs),qs5dt)

        xcwt = 1.0/( elf(kgs(mgs)) +cw(kgs(mgs))*temcg(mgs) )

        xwt1 = xhwet2(kgs(mgs))*(qss0(mgs)-qwv(mgs)) -tka(kgs(mgs))*temcg(mgs)

        qhwet(mgs) =  MAX( 0.0, ( xhv*chw(mgs)*xwt1*xhwet1(kgs(mgs))*xcwt &
                                  + ( 1.0 -ci(kgs(mgs))*temcg(mgs)*xcwt )* ( qhacip(mgs)+qhacsp(mgs) )  )  )

!
!  evaluate shedding rate (effective range is 243 < T < 273 due to other "if" checks below)
!
        qhshr(mgs) = 0.0
        IF ( qhwet(mgs) < qhdry(mgs) .AND. qhwet(mgs) > 0.0 ) THEN
          qhdry(mgs) = 0.0                                        ! Wet growth
          qhshr(mgs) = qhwet(mgs) -(qhacw(mgs) +qhacr(mgs))
        ELSE                                         ! Dry growth (defaults here if qhwet<0)
          qhwet(mgs) = 0.0
          qhshr(mgs) = 0.0
        END IF
!
!  Special shedding case when warmer than freezing
!
        IF ( temg(mgs) > tfr ) THEN
          qhwet(mgs) = 0.0
          qhdry(mgs) = 0.0
          qhshr(mgs) =  -qhacr(mgs) -qhacw(mgs)-qhacip(mgs)-qhacsp(mgs)
        END IF
!
!  Special no-shedding (dry) case when T<243....
!
        IF ( temg(mgs) < 243.15 ) THEN
          qhwet(mgs) = 0.0
          qhshr(mgs) = 0.0
        END IF
!
!  Reset some vars if wet particle surface due to shedding....
!
        IF ( qhshr(mgs) < 0.0 ) THEN
          qhaci(mgs) = qhacip(mgs)
          qhacs(mgs) = qhacsp(mgs)
          qhdpv(mgs) = 0.0
          qhsbv(mgs) = 0.0
        END IF

!
!  evaporation/condensation on wet snow and hail (NOT USED)
!

!
!  evaporation of rain
!
        xce = xxcev(kgs(mgs))*(ssw(mgs)-1.0)*(1.0/(xav(kgs(mgs))/tsqr(mgs)+xbv(kgs(mgs))/qvs(mgs)))
        xrv = (xrcev1(kgs(mgs))*(rwdia(mgs)) +  xrcev2(kgs(mgs))*(rwdia(mgs)**((3.0+br)/2.0)))
        qrcev(mgs) = MAX(MIN(xce*crw(mgs)*xrv, 0.0), -qr5dt)
!
!  vapor to pristine ice crystals
!
        qiint(mgs) = 0.0
        IF ( ssi(mgs) > 1.0 ) THEN
          dqisdt(mgs)= (qwv(mgs)-qis(mgs))/ (1.0 + xiint(kgs(mgs))*qis(mgs)/tsqr(mgs))
          cnnt       = cci(mgs)
          qiint(mgs) = (1.0/dtp) *MIN((1.0E-12)*cnnt/rho0x(mgs), 0.50*dqisdt(mgs))
        END IF
!
!  Domain totals for source terms
!
!
!  vapor
!
        tvqssbv = tvqssbv - il5(mgs)*qssbv(mgs)*rhobar(igs(mgs),jy,kgs(mgs))
        tvqhsbv = tvqhsbv - il5(mgs)*qhsbv(mgs)*rhobar(igs(mgs),jy,kgs(mgs))
        tvqrcev = tvqrcev - qrcev(mgs)*rhobar(igs(mgs),jy,kgs(mgs))
        tvqcevp = tvqcevp + 0.0
!  tvqisub = tvqisub + 0.0
!
!  cloud water
!
        tcqccnd = tcqccnd     + 0.0
        tcqcmli = tcqcmli     + 0.0
!
!  rain
!
        trqracw = trqracw + qracw(mgs)*rhobar(igs(mgs),jy,kgs(mgs))
        trqrcnw = trqrcnw + qrcnw(mgs)*rhobar(igs(mgs),jy,kgs(mgs))
        trqsacw = trqsacw + (1-il5(mgs))*qsacw(mgs)*rhobar(igs(mgs),jy,kgs(mgs))
        trqsmlr = trqsmlr - (1-il5(mgs))*qsmlr(mgs)*rhobar(igs(mgs),jy,kgs(mgs))
        trqhmlr = trqhmlr - (1-il5(mgs))*qhmlr(mgs)*rhobar(igs(mgs),jy,kgs(mgs))
        trqhshr = trqhshr - qhshr(mgs)*rhobar(igs(mgs),jy,kgs(mgs))
!
!  cloud ice
!
        tiqiint = tiqiint + il5(mgs)*qiint(mgs)*rhobar(igs(mgs),jy,kgs(mgs))
        tiqidep = tiqidep + 0.0
        tiqifzc = tiqifzc + 0.0
!
!  snow
!
        tsqsacw = tsqsacw + il5(mgs)*qsacw(mgs)*rhobar(igs(mgs),jy,kgs(mgs))
        tsqscni = tsqscni + il5(mgs)*qscni(mgs)*rhobar(igs(mgs),jy,kgs(mgs))
        tsqsaci = tsqsaci + il5(mgs)*qsaci(mgs)*rhobar(igs(mgs),jy,kgs(mgs))
        tsqsfi  = tsqsfi  + il5(mgs)*qsfi(mgs) *rhobar(igs(mgs),jy,kgs(mgs))
        tsqsfw  = tsqsfw  + il5(mgs)*qsfw(mgs) *rhobar(igs(mgs),jy,kgs(mgs))
        tsqraci = tsqraci + il5(mgs)*il3(mgs)*qraci(mgs)*rhobar(igs(mgs),jy,kgs(mgs))
        tsqiacr = tsqiacr + il5(mgs)*il3(mgs)*qiacr(mgs)*rhobar(igs(mgs),jy,kgs(mgs))
        tsqsacr = tsqsacr + il5(mgs)*il2(mgs)*qsacr(mgs)*rhobar(igs(mgs),jy,kgs(mgs))
        tsqsdpv = tsqsdpv + il5(mgs)*qsdpv(mgs)*rhobar(igs(mgs),jy,kgs(mgs))
!
!  hail/graupel
!
        thqhcns = thqhcns + qhcns(mgs)*rhobar(igs(mgs),jy,kgs(mgs))
        thqiacr = thqiacr + il5(mgs)*(1-il3(mgs))*qiacr(mgs)*rhobar(igs(mgs),jy,kgs(mgs))
        thqraci = thqraci + il5(mgs)*(1-il3(mgs))*qraci(mgs)*rhobar(igs(mgs),jy,kgs(mgs))
        thqracs = thqracs + il5(mgs)*(1-il2(mgs))*qracs(mgs)*rhobar(igs(mgs),jy,kgs(mgs))
        thqsacr = thqsacr + il5(mgs)*(1-il2(mgs))*qsacr(mgs)*rhobar(igs(mgs),jy,kgs(mgs))
        thqhdpv = thqhdpv + il5(mgs)*qhdpv(mgs)*rhobar(igs(mgs),jy,kgs(mgs))
        thqrfrz = thqrfrz + il5(mgs)*qrfrz(mgs)*rhobar(igs(mgs),jy,kgs(mgs))
!  thqhacr = thqhacr + qhacr(mgs)*rhobar(igs(mgs),jy,kgs(mgs))           ! MSG see below instead 2/6/04
!  thqhacw = thqhacw + qhacw(mgs)*rhobar(igs(mgs),jy,kgs(mgs))           ! MSG see below instead 2/6/04
!  thqhacs = thqhacs + qhacs(mgs)*rhobar(igs(mgs),jy,kgs(mgs))           ! MSG see below instead 2/6/04
!  thqhaci = thqhaci + qhaci(mgs)*rhobar(igs(mgs),jy,kgs(mgs))           ! MSG see below instead 2/6/04

!
!  hail/graupel and rain (based upon wet growth budget)  !MSG added on 2/6/04
!
!--
        IF ( temg(mgs) >= 273.15 ) THEN
          trqhaci = trqhaci + qhaci(mgs)*rhobar(igs(mgs),jy,kgs(mgs))
          trqhacs = trqhacs + qhacs(mgs)*rhobar(igs(mgs),jy,kgs(mgs))
        ELSE
          thqhaci = thqhaci + qhaci(mgs)*rhobar(igs(mgs),jy,kgs(mgs))
          thqhacs = thqhacs + qhacs(mgs)*rhobar(igs(mgs),jy,kgs(mgs))
        END IF

        IF ((qhwet(mgs) > 0.0).OR.( temg(mgs) >= 273.15)) THEN
          trqhacw = trqhacw + qhacw(mgs)*rhobar(igs(mgs),jy,kgs(mgs))
          trqhacr = trqhacr + qhacr(mgs)*rhobar(igs(mgs),jy,kgs(mgs))
        ELSE
          thqhacw = thqhacw + qhacw(mgs)*rhobar(igs(mgs),jy,kgs(mgs))
          thqhacr = thqhacr + qhacr(mgs)*rhobar(igs(mgs),jy,kgs(mgs))
        END IF

        thqhwet  = thqhwet + qhwet(mgs)*rhobar(igs(mgs),jy,kgs(mgs))    !MSG qhwet is positive or zero here

!--
!  end of totals
!
      END DO
!
      5002 CONTINUE
!
!
      IF (ndebug > 1 ) PRINT*,'dbg = 8'


!  rain, snow, hail fluxes due to gravity    !MSG an() was ad() in SAM version
!
      DO mgs = 1,ngscnt
        hwflx(mgs) = rho0x(mgs)*qh(igs(mgs),jgs,kgs(mgs))*vthbar(mgs)  !an() --> qh()
        piflx(mgs) = rho0x(mgs)*qi(igs(mgs),jgs,kgs(mgs))*vtibar(mgs)  !an() --> qi()
        cwflx(mgs) = rho0x(mgs)*qc(igs(mgs),jgs,kgs(mgs))*vtwbar(mgs)  !an() --> qc()
        rwflx(mgs) = rho0x(mgs)*qr(igs(mgs),jgs,kgs(mgs))*vtrbar(mgs)  !an() --> qr()
        swflx(mgs) = rho0x(mgs)*qs(igs(mgs),jgs,kgs(mgs))*vtsbar(mgs)  !an() --> qs()
      END DO

!
!  Compute total-domain content (g/m^3) before production rates
!
      DO mgs = 1,ngscnt                        !MSG domain total of each species before microphysics
        tqvbp = tqvbp + qwvp(mgs)*rhobar(igs(mgs),jy,kgs(mgs))
        tqcbp = tqcbp + qcw(mgs)*rhobar(igs(mgs),jy,kgs(mgs))
        tqibp = tqibp + qci(mgs)*rhobar(igs(mgs),jy,kgs(mgs))
        tqrbp = tqrbp + qrw(mgs)*rhobar(igs(mgs),jy,kgs(mgs))
        tqsbp = tqsbp + qsw(mgs)*rhobar(igs(mgs),jy,kgs(mgs))
        tqhbp = tqhbp + qhw(mgs)*rhobar(igs(mgs),jy,kgs(mgs))
        sumb = sumb+tqvbp+tqcbp+tqibp+tqrbp+tqsbp+tqhbp
      END DO

!
! CALCULATE RATE TOTALS
!
      DO mgs = 1,ngscnt

!
        pqwvi(mgs) =  il5(mgs)*( -qhsbv(mgs) -qssbv(mgs)              )           - qrcev(mgs)
        pqwvd(mgs) =  il5(mgs)*( -qhdpv(mgs) -qsdpv(mgs) - qiint(mgs) )
!
        pqcii(mgs) =  il5(mgs)*qiint(mgs)
        pqcid(mgs) =  il5(mgs)*( -qscni(mgs) -qsaci(mgs) -qraci(mgs) -qsfi(mgs))  - qhaci(mgs)
!
        pqcwi(mgs) =  0.0
        pqcwd(mgs) =  (-il5(mgs)*qsfw(mgs))  -qracw(mgs) -qsacw(mgs) -qrcnw(mgs) -qhacw(mgs)
!
        pqrwi(mgs) =  qracw(mgs) +qrcnw(mgs) +(1-il5(mgs))*(qsacw(mgs)-qhmlr(mgs) -qsmlr(mgs)) -qhshr(mgs)
        pqrwd(mgs) =  il5(mgs)*(-qiacr(mgs) -qrfrz(mgs)  -qsacr(mgs)) +qrcev(mgs) -qhacr(mgs)
!
        pqswi(mgs) =  il5(mgs)*( qsacw(mgs) +qscni(mgs) +qsaci(mgs) + qsfi(mgs) +qsfw(mgs) &
                                 +il3(mgs)*(qraci(mgs) +qiacr(mgs)) +il2(mgs)*qsacr(mgs) +qsdpv(mgs)     )
        pqswd(mgs) = -qhcns(mgs) -qhacs(mgs) +(1-il5(mgs))*qsmlr(mgs)  +il5(mgs)*(qssbv(mgs) -(1-il2(mgs))*qracs(mgs))
!
        pqhwi(mgs) =  qhcns(mgs) +qhacr(mgs) +qhacw(mgs) +qhacs(mgs) +qhaci(mgs) &
            + il5(mgs)*( (1-il3(mgs))*(qraci(mgs)+qiacr(mgs)) +(1-il2(mgs))*(qsacr(mgs)+qracs(mgs))+qhdpv(mgs)+qrfrz(mgs))
        pqhwd(mgs) =  qhshr(mgs) +(1-il5(mgs))*qhmlr(mgs) + il5(mgs)*qhsbv(mgs)
!
        ptotal(mgs) = pqwvi(mgs) +pqwvd(mgs) + pqcwi(mgs) +pqcwd(mgs) + pqcii(mgs) +pqcid(mgs) + &
                      pqrwi(mgs) +pqrwd(mgs) + pqswi(mgs) +pqswd(mgs) + pqhwi(mgs) +pqhwd(mgs)
!
        IF(ABS(ptotal(mgs)) > 1.e-7)THEN                             !on both these lines...
          PRINT*,'NOTICE:PTOTAL>1e-7 ', mgs, kgs(mgs), ptotal(mgs)   !...I changed 1e-9 to 1e-7
        END IF
!
        psum = psum + ptotal(mgs)
!
      END DO
!
!
!  latent heating from phase changes (except qcw, qci cond, and evap)
!   (22 processes involve phase changes, 10 do not)
!
      DO mgs = 1,ngscnt
        pfrz(mgs) = (1.-il5(mgs))*(qhmlr(mgs) + qsmlr(mgs))             &
                  + ( il5(mgs)  )*(qiacr(mgs)+qsacr(mgs)+ qsfw(mgs)+qrfrz(mgs)+qsacw(mgs)+qhacw(mgs)+qhacr(mgs)+qhshr(mgs))
        psub(mgs) = ( il5(mgs)  )*(qhdpv(mgs)+qhsbv(mgs)+qiint(mgs)+qsdpv(mgs)+qssbv(mgs))
        pvap(mgs) = qrcev(mgs)
        ptem(mgs) = cc3(kgs(mgs))*pfrz(mgs) + cc5(kgs(mgs))*psub(mgs) + cc4(kgs(mgs))*pvap(mgs)
        thetap(mgs) = thetap(mgs) + dtp*ptem(mgs)

!
!  partitioned domain-total heating and cooling rates
!  (all are adjusted again later within saturation adjustment)
!

        hfrz  = hfrz  + ( qiacr(mgs) + qsacr(mgs) + qsfw(mgs) +qsacw(mgs) &
                         +qhacw(mgs) + qhacr(mgs) + qhshr(mgs)+qrfrz(mgs))*il5(mgs) &
                                                                   *cc3(kgs(mgs))
        hdep  = hdep  + (qiint(mgs)+qhdpv(mgs)+qsdpv(mgs))*il5(mgs)*cc5(kgs(mgs))
        hcnd  = hcnd  + 0.0
        cevap = cevap +  qrcev(mgs)                                *cc4(kgs(mgs))
        cmelt = cmelt + (qhmlr(mgs)+qsmlr(mgs))*(1.-il5(mgs))      *cc3(kgs(mgs))
        csub  = csub  + (qhsbv(mgs)+qssbv(mgs))*il5(mgs)           *cc5(kgs(mgs))

      END DO
      9004 CONTINUE
!
!  sum the sources and sinks for qwvp, qcw, qci, qrw, qsw
!
      DO mgs = 1,ngscnt
        qwvp(mgs)= qwvp(mgs) + dtp*(pqwvi(mgs)+pqwvd(mgs))      !initial qwvp is being adjusted by all source/sink
        qcw(mgs) = qcw(mgs) +  dtp*(pqcwi(mgs)+pqcwd(mgs))
        qci(mgs) = qci(mgs) +  dtp*(pqcii(mgs)+pqcid(mgs))
        qrw(mgs) = qrw(mgs) +  dtp*(pqrwi(mgs)+pqrwd(mgs))
        qsw(mgs) = qsw(mgs) +  dtp*(pqswi(mgs)+pqswd(mgs))
        qhw(mgs) = qhw(mgs) +  dtp*(pqhwi(mgs)+pqhwd(mgs))

      END DO
!
!  domain-total content (g/m^3) before saturation adjustment
!
      DO mgs = 1,ngscnt
        tqvap = tqvap + qwvp(mgs)*rhobar(igs(mgs),jy,kgs(mgs))
        tqcap = tqcap + qcw(mgs)*rhobar(igs(mgs),jy,kgs(mgs))
        tqiap = tqiap + qci(mgs)*rhobar(igs(mgs),jy,kgs(mgs))
        tqrap = tqrap + qrw(mgs)*rhobar(igs(mgs),jy,kgs(mgs))
        tqsap = tqsap + qsw(mgs)*rhobar(igs(mgs),jy,kgs(mgs))
        tqhap = tqhap + qhw(mgs)*rhobar(igs(mgs),jy,kgs(mgs))
        suma = suma+tqvap+tqcap+tqiap+tqrap+tqsap+tqhap
      END DO


!
      IF (ndebug > 1 ) PRINT*,'dbg = 10a'

!
!  set up temperature and vapor arrays
!
      IF ( ndebug > 1 ) PRINT*,'dbg = 10.1'
!
      DO mgs = 1,ngscnt
        pqs(mgs) = (380.0)/(pres(mgs))
        theta(mgs) = thetap(mgs) + theta0(mgs)
        qvap(mgs) = MAX( (qwvp(mgs) + qv0(mgs)), 0.0 )         !MSG Current total qwv
        temg(mgs) = theta(mgs)*( pres(mgs) / poo ) ** rcp
      END DO
!
!  melting of cloud ice
!
      IF ( ndebug > 1 ) PRINT*,'dbg = 10.2'
!
      DO mgs = 1,ngscnt
        IF( temg(mgs) > tfr .AND. qci(mgs) > 0.0 ) THEN
          qimlw(mgs) = - qci(mgs)/dtp                         !MSG Rate of cloudice melting
          tcqcmli = tcqcmli - qimlw(mgs)*rhobar(igs(mgs),jy,kgs(mgs))        !MSG updated 11/22/03 (domain-total rate of qc increase)
          thetap(mgs) = thetap(mgs) - cc3(kgs(mgs))*qci(mgs)  !MSG heat decrease
          cmelt   = cmelt   + cc3(kgs(mgs))*qimlw(mgs)        !MGS cooling rate
          qcw(mgs) = qcw(mgs) + qci(mgs)
          qci(mgs) = 0.0
        END IF
      END DO
!
!
!  homogeneous freezing of cloud water
!
      IF ( ndebug > 1 ) PRINT*,'dbg = 10.3'
!
      DO mgs = 1,ngscnt
        IF( temg(mgs) < thnuc .AND. qcw(mgs) > 0.0 ) THEN
          qwfzi(mgs)  = -qcw(mgs)/dtp                         ! MSG Rate of clouwater freezing
          tiqifzc = tiqifzc - qwfzi(mgs)*rhobar(igs(mgs),jy,kgs(mgs))        ! MSG updated 11/22/03 (domain-total rate of qi increase)
          thetap(mgs) = thetap(mgs) + cc3(kgs(mgs))*qcw(mgs)  ! MSG heat increase
          hfrz    = hfrz    + cc3(kgs(mgs))*(-qwfzi(mgs))     ! MSG heating rate
          qci(mgs) = qci(mgs) + qcw(mgs)
          qcw(mgs) = 0.0
        END IF
      END DO

!
!  Saturation adjustment iteration procedure
!
!  Modified Straka adjustment (nearly identical to Tao et al. 1989 MWR)
!
!

!
!  reset temporaries for cloud particles and vapor
!
      IF ( ndebug > 1 ) PRINT*,'dbg = 10.4'
      DO mgs = 1,ngscnt
        qwv(mgs) = MAX( 0.0, qvap(mgs) )
        qcw(mgs) = MAX( 0.0, qcw(mgs) )
        qci(mgs) = MAX( 0.0, qci(mgs) )
        ptotsat(mgs) = qvap(mgs)+qci(mgs)+qcw(mgs)       !MSG updated just before sat adj. (qwv+qci+qcw)
        qcevpcnd(mgs) = 0.0
        qisubdep(mgs) = 0.0
      END DO
!
      tqvb = tqvap              !MSG domain-total vapor perturb. prior to sat adj.
      DO mgs = 1,ngscnt         !MSG domain-total qcw and qci prior to sat adj.
        tqcb = tqcb + qcw(mgs)*rhobar(igs(mgs),jy,kgs(mgs))
        tqib = tqib + qci(mgs)*rhobar(igs(mgs),jy,kgs(mgs))
      END DO

      DO mgs = 1,ngscnt
        theta(mgs) = thetap(mgs) + theta0(mgs)
        temg(mgs) = theta(mgs)*( pres(mgs) / poo ) ** rcp
        temcg(mgs) = temg(mgs) - tfr
        ltemq = nint((temg(mgs)-163.15)/fqsat+1.5)
        ltemq = MIN(MAX(ltemq,1),nqsat)
        qvs(mgs) = pqs(mgs)*tabqvs(ltemq)
        qis(mgs) = pqs(mgs)*tabqis(ltemq)
        IF ( temg(mgs) < tfr ) THEN
          IF( qcw(mgs) >= 0.0 .AND. qci(mgs) == 0.0 ) qss(mgs) = qvs(mgs)
          IF( qcw(mgs) == 0.0 .AND. qci(mgs) > 0.0)  qss(mgs) = qis(mgs)
          IF( qcw(mgs) > 0.0 .AND. qci(mgs) > 0.0)  qss(mgs) = (qcw(mgs)*qvs(mgs) + qci(mgs)*qis(mgs)) &
                                                                   / (qcw(mgs) + qci(mgs))
        ELSE
          qss(mgs) = qvs(mgs)
        END IF
      END DO
!
!  iterate  adjustment
!
      IF ( ndebug > 1 ) PRINT*,'dbg = 10.5'
      DO itertd = 1,2
!
        DO mgs = 1,ngscnt
!
!  calculate super/sub-saturation
!
          dqcw(mgs) = 0.0
          dqci(mgs) = 0.0
          dqwv(mgs) = ( qwv(mgs) - qss(mgs) )
!
!  evaporation and sublimation adjustment
!
          IF( dqwv(mgs) < 0. ) THEN
            IF( qcw(mgs) > -dqwv(mgs) ) THEN          !evap some of qc
              dqcw(mgs) = dqwv(mgs)
              dqwv(mgs) = 0.
            ELSE                          !Evap all of qc
              dqcw(mgs) = -qcw(mgs)
              dqwv(mgs) = dqwv(mgs) + qcw(mgs)
            END IF
!
            IF( qci(mgs) > -dqwv(mgs) ) THEN          !sublimate some of qi
              dqci(mgs) = dqwv(mgs)
              dqwv(mgs) = 0.
            ELSE                              !Sublimate all of qi
              dqci(mgs) = -qci(mgs)
              dqwv(mgs) = dqwv(mgs) + qci(mgs)
            END IF
!
            qwvp(mgs) = qwvp(mgs) - ( dqcw(mgs) + dqci(mgs) )     !Increase vapor

            qcw(mgs) = qcw(mgs) + dqcw(mgs)               !Decrease cloudwater (dqcw<0)
            qci(mgs) = qci(mgs) + dqci(mgs)                !Decrease cloudice   (dqci<0)
            thetap(mgs) = thetap(mgs) + cpi/pi0(mgs)*(elv(kgs(mgs))*dqcw(mgs) +els(kgs(mgs))*dqci(mgs))
            qcevpcnd(mgs) = qcevpcnd(mgs) + (dqcw(mgs)/dtp)
            qisubdep(mgs) = qisubdep(mgs) + (dqci(mgs)/dtp)
          END IF
!
! condensation/deposition
!
          IF( dqwv(mgs) >= 0. ) THEN
!
            fracl(mgs) = 1.0
            fraci(mgs) = 0.0
            IF ( temg(mgs) < tfr .AND. temg(mgs) > thnuc ) THEN
              fracl(mgs) = MAX(MIN(1.,(temg(mgs)-233.15)/(20.)),0.0)
              fraci(mgs) = 1.0-fracl(mgs)
            END IF
            IF ( temg(mgs) <= thnuc ) THEN
              fraci(mgs) = 1.0
              fracl(mgs) = 0.0
            END IF
            fraci(mgs) = 1.0-fracl(mgs)
!
            gamss = (elv(kgs(mgs))*fracl(mgs) + els(kgs(mgs))*fraci(mgs))/ (pi0(mgs)*cp)
!
            IF ( temg(mgs) < tfr ) THEN
              IF (qcw(mgs) >= 0.0 .AND. qci(mgs) <= 0.0 ) THEN
                dqvcnd(mgs) = dqwv(mgs)/(1. + cqv1(kgs(mgs))*qss(mgs)/((temg(mgs)-cbw)**2))
              END IF
              IF( qcw(mgs) == 0.0 .AND. qci(mgs) > 0.0 ) THEN
                dqvcnd(mgs) = dqwv(mgs)/(1. + cqv2(kgs(mgs))*qss(mgs)/((temg(mgs)-cbi)**2))
              END IF
              IF ( qcw(mgs) > 0.0 .AND. qci(mgs) > 0.0 ) THEN
                cdw = caw*pi0(mgs)*tfrcbw/((temg(mgs)-cbw)**2)
                cdi = cai*pi0(mgs)*tfrcbi/((temg(mgs)-cbi)**2)
                denom1 = qcw(mgs) + qci(mgs)
                denom2 = 1.0 + gamss*(qcw(mgs)*qvs(mgs)*cdw + qci(mgs)*qis(mgs)*cdi) / denom1
                dqvcnd(mgs) =  dqwv(mgs) / denom2
              END IF
            END IF
!
            IF ( temg(mgs) >= tfr ) THEN
              dqvcnd(mgs) = dqwv(mgs)/(1. + cqv1(kgs(mgs))*qss(mgs)/ ((temg(mgs)-cbw)**2))
            END IF
!
            dqcw(mgs) = dqvcnd(mgs)*fracl(mgs)
            dqci(mgs) = dqvcnd(mgs)*fraci(mgs)
!
            thetap(mgs) = thetap(mgs) + (elv(kgs(mgs))*dqcw(mgs) + els(kgs(mgs))*dqci(mgs))/ (pi0(mgs)*cp)
            qwvp(mgs) = qwvp(mgs) - ( dqvcnd(mgs) )          !Decrease vapor
            qcw(mgs) = qcw(mgs) + dqcw(mgs)               !Increase cloudwater (dqcw>0)
            qci(mgs) = qci(mgs) + dqci(mgs)               !Increase cloudice   (dqci>0)
            qcevpcnd(mgs) = qcevpcnd(mgs) + (dqcw(mgs)/dtp)
            qisubdep(mgs) = qisubdep(mgs) + (dqci(mgs)/dtp)

!
          END IF
        END DO
!
        IF ( ndebug > 1 ) PRINT*,'dbg = 10.51'
        DO mgs = 1,ngscnt
          theta(mgs) = thetap(mgs) + theta0(mgs)
          temg(mgs) = theta(mgs)*( pres(mgs) / poo ) ** rcp
          qvap(mgs) =MAX((qwvp(mgs) + qv0(mgs)), 0.0)
          temcg(mgs) = temg(mgs) - tfr

          ltemq = nint((temg(mgs)-163.15)/fqsat+1.5)
          ltemq = MIN(MAX(ltemq,1),nqsat)
          qvs(mgs) = pqs(mgs)*tabqvs(ltemq)
          qis(mgs) = pqs(mgs)*tabqis(ltemq)
          qcw(mgs) = MAX( 0.0, qcw(mgs) )
          qwv(mgs) = MAX( 0.0, qvap(mgs))
          qci(mgs) = MAX( 0.0, qci(mgs) )

          IF ( temg(mgs) < tfr ) THEN
            IF( qcw(mgs) >= 0.0 .AND. qci(mgs) == 0.0 ) qss(mgs) = qvs(mgs)
            IF( qcw(mgs) == 0.0 .AND. qci(mgs) > 0.0)  qss(mgs) = qis(mgs)
            IF( qcw(mgs) > 0.0 .AND. qci(mgs) > 0.0)  qss(mgs) = (qcw(mgs)*qvs(mgs) + qci(mgs)*qis(mgs)) &
                                                                   / (qcw(mgs) + qci(mgs))
          ELSE
            qss(mgs) = qvs(mgs)
          END IF
        END DO
        IF ( ndebug > 1 ) PRINT*,'dbg = 10.52'
!
!  end the saturation adjustment iteration loop
!
      END DO
      IF ( ndebug > 1 ) PRINT*,'dbg = 10.6'

      DO mgs = 1,ngscnt                     !MSG net at each gpt after all sat adj. iterations are finished
        qcevp(mgs) = MIN(qcevpcnd(mgs),0.)  ! qcevp <=0
        qisub(mgs) = MIN(qisubdep(mgs),0.)  ! qisub <=0
        qccnd(mgs) = MAX(qcevpcnd(mgs),0.)  ! qccnd >=0
        qidep(mgs) = MAX(qisubdep(mgs),0.)  ! qidep >=0
      END DO

      DO mgs = 1,ngscnt
        tcqccnd = tcqccnd + qccnd(mgs)*rhobar(igs(mgs),jy,kgs(mgs))    ! MSG updated 2/12/05 domain-total condensation
        tvqcevp = tvqcevp - qcevp(mgs)*rhobar(igs(mgs),jy,kgs(mgs))    ! MSG updated 2/12/05 domain-total qc evaporation
        tvqisub = tvqisub - qisub(mgs)*rhobar(igs(mgs),jy,kgs(mgs))    ! MSG updated 2/12/05 domain-total sublimation
        tiqidep = tiqidep + qidep(mgs)*rhobar(igs(mgs),jy,kgs(mgs))    ! MSG updated 2/12/05 domain-total deposition

        hcnd = hcnd    + qccnd(mgs)*cc4(kgs(mgs))    ! MSG Update domain-total heating rate via qc condensation
        cevap = cevap   + qcevp(mgs)*cc4(kgs(mgs))    ! MSG Update domain-total cooling rate via qc evap
        csub = csub    + qisub(mgs)*cc5(kgs(mgs))    ! MSG Update domain-total cooling rate via qi sublim
        hdep = hdep    + qidep(mgs)*cc5(kgs(mgs))    ! MSG Update domain-total heating rate via qi deposition
      END DO

!
!  Compute vapor, ice, and cloud totals after saturation adjustment.
!
      DO mgs = 1,ngscnt
        IF(ABS(ptotsat(mgs)-qwv(mgs)-qci(mgs)-qcw(mgs)) > 1.e-7)THEN   !on both these lines...
          PRINT*,'NOTICE:PTOTSAT>1e-7 ', mgs, kgs(mgs), ptotsat(mgs)   !...I changed from e-9 to e-7. ~NS (4/15)
        END IF
      END DO

      IF (ndebug > 1 ) PRINT*,'dbg = 10b'
!
      DO mgs = 1,ngscnt
        tqva = tqva + qwvp(mgs)*rhobar(igs(mgs),jy,kgs(mgs))
        tqca = tqca +  qcw(mgs)*rhobar(igs(mgs),jy,kgs(mgs))
        tqia = tqia +  qci(mgs)*rhobar(igs(mgs),jy,kgs(mgs))
      END DO

      tqv=tqva-tqvb          ! Change in vapor due to sat adj. (Should equal +tvqcevp+tvqisub-tiqidep-tcqccnd  )
      tqc=tqca-tqcb          ! Change in cloud due to sat adj. (Should equal +tcqccnd+tcqcmli-tvqcevp-tiqifzc  )
      tqi=tqia-tqib          ! Change in ice due to sat adj.   (Should equal +tiqidep+tiqifzc-tvqisub-tcqcmli  )
!
!
!
!  end of saturation adjustment
!
!  scatter precipitation fluxes, and thetap, and hydrometeors
!
!DIR$ IVDEP
      DO mgs = 1,ngscnt
        pflux(igs(mgs),kgs(mgs)) = piflx(mgs)
        cflux(igs(mgs),kgs(mgs)) = cwflx(mgs)
        rflux(igs(mgs),kgs(mgs)) = rwflx(mgs)
        sflux(igs(mgs),kgs(mgs)) = swflx(mgs)
        hflux(igs(mgs),kgs(mgs)) = hwflx(mgs)
        ptprt(igs(mgs),jy,kgs(mgs)) = thetap(mgs)
        qv(igs(mgs),jy,kgs(mgs)) = qvbar(igs(mgs),jy,kgs(mgs)) + qwvp(mgs)
        qc(igs(mgs),jy,kgs(mgs)) = qcw(mgs)  + MIN( qc(igs(mgs),jy,kgs(mgs)), 0.0 ) !msg putting any neg gibbs values back
        qi(igs(mgs),jy,kgs(mgs)) = qci(mgs)  + MIN( qi(igs(mgs),jy,kgs(mgs)), 0.0 )
        qr(igs(mgs),jy,kgs(mgs)) = qrw(mgs)  + MIN( qr(igs(mgs),jy,kgs(mgs)), 0.0 )
        qs(igs(mgs),jy,kgs(mgs)) = qsw(mgs)  + MIN( qs(igs(mgs),jy,kgs(mgs)), 0.0 )
        qh(igs(mgs),jy,kgs(mgs)) = qhw(mgs)  + MIN( qh(igs(mgs),jy,kgs(mgs)), 0.0 )
      END DO
!
!

!
!
      9998 CONTINUE
!
      IF ( kz > nz-kstag-1 .AND. ix > nx-istag ) THEN
        GO TO 1200
      ELSE
        nzmpb = kz
      END IF

      IF ( ix+1 > nx ) THEN
        nxmpb = 1
      ELSE
        nxmpb = ix+1
      END IF

    END DO
    1200 CONTINUE
!
!  end of gather scatter
!
!  precipitation fallout contributions
!
!
    IF( sadvopt /= 4) THEN                  ! Leapfrog scheme
      dtz1 = 2*dzc(kz)*dtp/1.0        !MSG exchanged 1/dzc for dz
    ELSE
      dtz1 = dzc(kz)*dtp/1.0        !MSG exchanged 1/dzc for dz
    ENDIF
!
!  MSG Technically we should really be re-computing the fluxes based upon new mixing ratio's but we don't.
    IF (ndebug > 1 ) PRINT*,'dbg = 10g'
    DO kz = 1,nz-kstag-1
      DO ix = 1,nx-istag
        qi(ix,jy,kz) = qi(ix,jy,kz) + dtz1*(pflux(ix,kz+1)-pflux(ix,kz))/rhobar(ix,jy,kz)
        qc(ix,jy,kz) = qc(ix,jy,kz) + dtz1*(cflux(ix,kz+1)-cflux(ix,kz))/rhobar(ix,jy,kz)
        qr(ix,jy,kz) = qr(ix,jy,kz) + dtz1*(rflux(ix,kz+1)-rflux(ix,kz))/rhobar(ix,jy,kz)
        qs(ix,jy,kz) = qs(ix,jy,kz) + dtz1*(sflux(ix,kz+1)-sflux(ix,kz))/rhobar(ix,jy,kz)
        qh(ix,jy,kz) = qh(ix,jy,kz) + dtz1*(hflux(ix,kz+1)-hflux(ix,kz))/rhobar(ix,jy,kz)
      END DO
    END DO
!
!
!  end of jy loop
!
  END DO
!
  IF (ndebug > 1 ) PRINT*,'dbg = 10h'
!
!
!  WRITE totals for source / sink terms
  IF (nrates == 1 ) THEN
!
    tvsum  = tvqrcev -tiqiint -thqhdpv -tsqsdpv +tvqhsbv +tvqssbv
    tcsum  = -tsqsfw -thqhacw -trqrcnw -trqracw -tsqsacw -trqsacw
    trsum  = -tvqrcev -thqhacr +trqrcnw -thqrfrz +trqhshr +trqhmlr +trqsmlr -tsqiacr -thqiacr -tsqsacr -thqsacr &
             +trqracw +trqsacw
    tisum  = -tsqsfi -thqhaci -tsqscni +tiqiint -tsqraci -thqraci  -tsqsaci
    tssum  = tsqsfw +tsqsfi -thqhacs +tsqscni -thqhcns -trqsmlr +tsqiacr +tsqraci +tsqsacr +tsqsdpv -tvqssbv &
            +tsqsaci +tsqsacw -thqracs
    thsum  = thqhacw +thqhacr +thqhacs +thqhaci +thqhcns +thqrfrz -trqhshr -trqhmlr +thqiacr +thqraci +thqsacr +thqhdpv &
             +thqracs -tvqhsbv
    tsumall =  tvsum + tcsum + tisum + trsum + tssum + thsum        !MSG gives ~1e-9 (machine PRECISION)
!
    IF (ndebug > 1) THEN    ! MSG only print details if debugging on
!
      WRITE(6,*) 'Sum species source/sink domain totals for all but Sat. Adj. (kg s^-1 m^-3)'
      WRITE(6,*) 'qv', tvsum            ! Total rate for qv (not including sat adj.)
      WRITE(6,*) 'qc', tcsum
      WRITE(6,*) 'qi', tisum
      WRITE(6,*) 'qr', trsum
      WRITE(6,*) 'qs', tssum
      WRITE(6,*) 'qh', thsum
!
      WRITE(6,*) 'Sum only rates for saturation adjustment (kg s^-1 m^-3)'
      WRITE(6,*) 'qv', tqv              !MSG total change in vapor only due to sat adj.
      WRITE(6,*) 'qc', tqc              !MSG change in cloud only due to sat adj.
      WRITE(6,*) 'qi', tqi              !MSG change in ice only due to sat adj.
!
      WRITE(6,*) 'Sum all species rates but not sat. adj. (kg s^-1 m^-3)'
      WRITE(6,*) 'sum', tsumall         !MSG - all sources/sink rates except for sat adj.(comes from tsum rates)
!    write(6,*) 'psum', psum           !MSG - all sources/sink rates except for sat adj.(comes from model rates)
!    write(6,*) 'sumchange', suma-sumb !MSG - change in content (g/m^3) over timestep (not incl. sat adj.)
!
      WRITE(6,*) 'Sum of all rates (kg s^-1 m^-3)'
      WRITE(6,*) 'tsa', tsumall+tqv+tqc+tqi  !MSG all source/sink rates including sat adj.
!    write(6,*) 'psa', psum+tqv+tqc+tqi     !MSG all source/sink rates including sat adj.
!
!    write(6,*) 'tqvp',tqvap-tqvbp          !Changes due to all (except sat adj or fallout)
!    write(6,*) 'tqcp',tqcap-tqcbp
!    write(6,*) 'tqip',tqiap-tqibp
!    write(6,*) 'tqrp',tqrap-tqrbp
!    write(6,*) 'tqsp',tqsap-tqsbp
!    write(6,*) 'tqhp',tqhap-tqhbp
!
    END IF
!
!
    IF (ndebug >= 1) THEN
      WRITE(6,*) 'Individual domain total rates (kg s^-1 m^-3)'
      WRITE(6,*) 'Using Gilmore et al. (2004b) terminology'
      WRITE(6,*) 'vapor sources'
!
      WRITE(6,*) 'qvevr ' , tvqrcev
      WRITE(6,*) 'qvevw ' , tvqcevp
      WRITE(6,*) 'qvsbi ' , tvqisub
      WRITE(6,*) 'qvsbs ' , tvqssbv
      WRITE(6,*) 'qvsbh ' , tvqhsbv

      WRITE(6,*) 'cloud water sources'

      WRITE(6,*) 'qwcdv ' , tcqccnd
      WRITE(6,*) 'qwmli ' , tcqcmli

      WRITE(6,*) 'rain sources'

      WRITE(6,*) 'qrhacr' , trqhacr   ! Added 2/6/04  (rain accreted
                                    ! and shed during same timestep)
      WRITE(6,*) 'qrmlh ' , trqhmlr
      WRITE(6,*) 'qracw ' , trqracw
      WRITE(6,*) 'qrcnw ' , trqrcnw
  !    write(6,*) 'qrshh ' , trqhshr  ! MSG ambiguous since sign can switch.
                                      ! Instead, use 4 individual terms (2/6/04)
      WRITE(6,*) 'qrhacw' , trqhacw   ! Added 2/6/04  (from shedding)
      WRITE(6,*) 'qrhacs' , trqhacs   ! Added 2/6/04  (from shedding)
      WRITE(6,*) 'qrmls ' , trqsmlr
      WRITE(6,*) 'qrsacw' , trqsacw
      WRITE(6,*) 'qrhaci' , trqhaci   ! Added 2/6/04  (from shedding)

      WRITE(6,*) 'cloud ice sources'

      WRITE(6,*) 'qidpv ' , tiqidep
      WRITE(6,*) 'qiint ' , tiqiint
      WRITE(6,*) 'qifzw ' , tiqifzc

      WRITE(6,*) 'snow sources'

      WRITE(6,*) 'qsfi  ' , tsqsfi
      WRITE(6,*) 'qsacw ' , tsqsacw
      WRITE(6,*) 'qsaci ' , tsqsaci
      WRITE(6,*) 'qsdpv ' , tsqsdpv
      WRITE(6,*) 'qsiacr' , tsqiacr
      WRITE(6,*) 'qsacr ' , tsqsacr
      WRITE(6,*) 'qsfw  ' , tsqsfw
      WRITE(6,*) 'qsraci' , tsqraci
      WRITE(6,*) 'qscni ' , tsqscni

      WRITE(6,*) 'hail sources'

      WRITE(6,*) 'qhacw ' , thqhacw   ! 2/6/04 (only that which is not shed as rain)
      WRITE(6,*) 'qhacr ' , thqhacr   ! 2/6/04 (only that which is not shed as rain)
      WRITE(6,*) 'qhdpv ' , thqhdpv
      WRITE(6,*) 'qhsacr' , thqsacr
      WRITE(6,*) 'qhwtr ' , thqhwet   ! 2/6/04 (only that which is not shed as rain)
      WRITE(6,*) 'qhacs ' , thqhacs   ! 2/6/04 (only that which is not shed as rain)
      WRITE(6,*) 'qhfzr ' , thqrfrz
      WRITE(6,*) 'qhaci ' , thqhaci   ! 2/6/04 (only that which is not shed as rain)
      WRITE(6,*) 'qhracs' , thqracs
      WRITE(6,*) 'qhcns ' , thqhcns
      WRITE(6,*) 'qhiacr' , thqiacr
      WRITE(6,*) 'qhraci' , thqraci

      WRITE(6,*) 'total heating and cooling rates'

      WRITE(6,*) 'hfrz '   , hfrz
      WRITE(6,*) 'hdep '   , hdep
      WRITE(6,*) 'hcnd '   , hcnd
      WRITE(6,*) 'cmelt'   , cmelt
      WRITE(6,*) 'cevap'   , cevap
      WRITE(6,*) 'csub '   , csub

    END IF

  END IF

  IF (ndebug > 1 ) PRINT*,'dbg = 11a'

  RETURN
END SUBROUTINE lfo_ice
!
!##########################################################################
!
! Routine from Numerical Recipes to replace other gamma function
! using 32-bit reals, this is accurate to 6th decimal place.
!
!#######################################################################

REAL FUNCTION gamma_lfo(xx)

  IMPLICIT NONE
  REAL, INTENT(IN) :: xx

! Double precision ser,stp,tmp,x,y,cof(6)

!  REAL*8 ser,stp,tmp,x,y,cof(6)
  DOUBLE PRECISION :: ser, tmp, x, y
  DOUBLE PRECISION, SAVE :: stp    = 2.5066282746310005D0
  DOUBLE PRECISION, SAVE :: cof(6) = (/  76.18009172947146D+0,          &
                                        -86.50532032941677D0,           &
                                         24.01409824083091D0,           &
                                        -1.231739572450155D0,           &
                                        0.1208650973866179D-2,          &
                                          -0.5395239384953D-5  /)
  INTEGER :: j

  x = xx
  y = x
  tmp = x + 5.5D0
  tmp = (x + 0.5D0)*LOG(tmp) - tmp
  ser = 1.000000000190015D0
  DO j=1,6
    y = y + 1.0D0
    ser = ser + cof(j)/y
  END DO
  gamma_lfo = EXP(tmp + LOG(stp*ser/x))

  RETURN
END FUNCTION gamma_lfo
