!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINE ICECVT                    ######
!######                                                      ######
!######                     Developed by                     ######
!######                                                      ######
!######    Goddard Cumulus Ensemble Modeling Group, NASA     ######
!######                                                      ######
!######     Center for Analysis and Prediction of Storms     ######
!######               University of Oklahoma                 ######
!######                                                      ######
!##################################################################
!##################################################################
!
!

SUBROUTINE icecvt(nx,ny,nz, d2t,                                        &
           ptprt,qv,qscalar,rhobar,ptbar,qvbar,pres,ppi,                &
           wgacr, scv, tca, dwv, zr, vr, zs, vs,                        &
           zg, vg, psaut, psaci,psacw, qsacw,praci,piacr,               &
           praut, pracw, psfw, psfi, dgacs, dgacw, dgaci,               &
           dgacr, pgacs, wgacs, qgacw,wgaci,qgacr,pgaut,                &
           pracs, psacr, qsacr, pgfr, tair, tairc,                      &
           tair3, tairc3,                                               &
           pr, pg, ps, dlt2, dlt3, rtair, tem2d1, tem2d2, mpteqnterms)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate and apply the microphysical contributions to the water/ice
!  and temperature fields. It uses Y.L. Lin's ice microphysics
!  parameterization scheme and the code is adopted from W. G. Tao's
!  ice package.
!
!-----------------------------------------------------------------------
!
!  REFERENCES:
!
!  Lin et. al. (1983)
!  TAo and Simpson (1989, JAS)
!
!-----------------------------------------------------------------------
!
!  AUTHOR: W.G. Tao, Goddard Cumulus Ensemble Modeling Group, NASA
!
!  MODIFICATION HISTORY:
!
!  12/29/93 (A. Sathye)
!  Adopted the original code and converted 2-D into 3-D.
!
!  01/05/94 (Y. Liu)
!  Formatting the code in accordance with the ARPS coding format.
!  Continue the conversion of 2-D to 3-D.
!
!  8/9/94 (MX)
!  A statement using unset variable r4f corrected.
!
!  10/14/94 (V. Wong, A. Sathye, Y. Liu, & X. Song)
!  Walked through the code and fixed a few bugs.
!
!  03/01/96 (Vince Wong)
!  Modified to prevent qs divided by zero.
!
!  10/23/1996 (M. Xue)
!  Altered codes in loops 150 and 200 to avoid trunction errors that
!  were causing the program to bomb.
!
!  1/23/1997 (J. Zong, M. Xue, V. Wong)
!  Corrected bugs in loop 150 and 200 introduced in modification
!  on 10/23/1996.
!
!  2/11/1997 (J. Zong)
!  Set a positive lower limit to qi when calculating deposition of
!  qi to avoid arithmetic exception.
!
!  2/26/1997 (J. Zong, M. Xue and Yuhe Liu)
!  Fractional power calculations are replaced by lookup table functions.
!
!  3/05/1997 (Fanyou Kong)
!  Modify the time levels used to calculate microphysical source and
!  sink terms, to make the scheme more physically rational;
!
!  07/10/97 (Fanyou Kong - CMRP)
!  Change sqrt(sqrt(zr**11)) to zr**2.75 to allow the code
!  executable on DEC Alpha system
!
!  5/19/1998 (M. Xue)
!  Changed reference density at surface from rhobar(1,1,2) to rhobar0.
!
!  11/18/98 (Keith Brewster)
!  Changed pibar to ppi (full pi).
!
!  July 2001
!  Correction to rainwater evaporation term.  See ARPSSUPPORT e-mails
!  by Vince Wong and Eric Kemp.
!
!  Fall 2002
!  Added fallopt option to specify which hydrometeor fall velocity
!  formulations -- original Lin-Tao or Ferrier -- should be used.
!  ***NOTE*** Additional work is required for some of the transfer
!  terms when the Ferrier fall speeds are used.
!
!  13 May 2003 (Eric Kemp)
!  Bug fix to Praut term.  Added rhobar to denominator following
!  Orville and Kopp, JAS, 1977, equation 17.  The exclusion of the
!  rhobar in the Lin et al., JAM, 1983 manuscript is a typo.  This
!  error was detected by Jianzhong Wang.
!
!-----------------------------------------------------------------------
!
!  ORIGINAL COMMENTS:
!
!  lin et al (83) ice phase microphysical processes
!  modified and coded by goddard cumulus ensemble modeling group
!  tao, simpson and mccumber's saturation technique (mwr, 1989)
!  tao and simpson (jas, 1989; tao, 1993)
!
!  d2t - leapfrog time step
!  dpt - deviation potential temperature field
!  dqv - deviation water vapor field
!  qcl - cloud water field
!  qrn - rain field
!  qci - cloud ice field
!  qcs - snow field
!  qcg - hail field
!  rho - air density
!  ta1 - base air potential temperature at the level
!  qa1 - base water vapor at the level
!  p0 - air pressure
!  ppi - exner function
!
!  nx: dimension in x-direction
!  nz: dimension in z-direction
!  note: physical domain extends from k=2 to k=nz-1 and i=2 to i=nx-1
!
!  c.g.s units
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

  INTEGER :: nx, ny, nz

  REAL :: ptprt (nx,ny,nz,nt)  ! Perturbation potential temperature (K)
  REAL :: qv    (nx,ny,nz,nt)  ! Water vapor specific humidity (g/g)

  REAL, TARGET :: qscalar(nx,ny,nz,nt,nscalar)

  REAL :: rhobar(nx,ny,nz)  ! Base state air density (g/cm**3)
  REAL :: pres  (nx,ny,nz)  ! Pressure ((g cm/s**2)/cm**2)
  REAL :: ptbar (nx,ny,nz)  ! Base state potential temperature (K)
  REAL :: qvbar (nx,ny,nz)  ! Base state air density (g/g)
  REAL :: ppi   (nx,ny,nz)  ! Exner function
!
!-----------------------------------------------------------------------
!
!  Local work arrays for production terms etc.
!
!-----------------------------------------------------------------------
!

  REAL :: wgacr(nx,ny)
  REAL :: psaut(nx,ny)
  REAL :: psaci(nx,ny)
  REAL :: psacw(nx,ny)
  REAL :: qsacw(nx,ny)
  REAL :: praci(nx,ny)
  REAL :: piacr(nx,ny)
  REAL :: praut(nx,ny)
  REAL :: pracw(nx,ny)
  REAL :: psfw (nx,ny)
  REAL :: psfi (nx,ny)
  REAL :: dgacs(nx,ny)
  REAL :: dgacw(nx,ny)
  REAL :: dgaci(nx,ny)
  REAL :: dgacr(nx,ny)
  REAL :: pgacs(nx,ny)
  REAL :: wgacs(nx,ny)
  REAL :: qgacw(nx,ny)
  REAL :: wgaci(nx,ny)
  REAL :: qgacr(nx,ny)
  REAL :: pgaut(nx,ny)
  REAL :: pracs(nx,ny)
  REAL :: psacr(nx,ny)
  REAL :: qsacr(nx,ny)
  REAL :: pgfr (nx,ny)

  REAL :: tair (nx,ny)
  REAL :: tairc(nx,ny)
  REAL :: tair3 (nx,ny)
  REAL :: tairc3(nx,ny)
  REAL :: pr   (nx,ny)
  REAL :: pg   (nx,ny)
  REAL :: ps   (nx,ny)
  REAL :: dlt2 (nx,ny)
  REAL :: dlt3 (nx,ny)
  REAL :: rtair(nx,ny)

  REAL :: scv  (nx,ny)
  REAL :: tca  (nx,ny)
  REAL :: dwv  (nx,ny)
  REAL :: zr   (nx,ny)
  REAL :: vr   (nx,ny)
  REAL :: zs   (nx,ny)
  REAL :: vs   (nx,ny)
  REAL :: zg   (nx,ny)
  REAL :: vg   (nx,ny)

  REAL :: tem2d1(nx,ny)
  REAL :: tem2d2(nx,ny)
!
!-----------------------------------------------------------------------
!
!  Following constants are defined and set in subroutine STCSTICE
!
!  (Suggest to move the following COMMON block to a inlcude file.)
!
!-----------------------------------------------------------------------
!
  REAL :: tnw,tns,tng,roqr,roqs,roqg
  REAL :: c76,c358,c172,c409,c218,c580,c610,c149,c879,c141
  REAL :: zrc,zgc,zsc,vrc,vgc,vsc
  REAL :: ag,bg,as,bs,aww,bww,bgh,bgq,bsh,bsq,bwh,bwq
  REAL :: alv,alf,als,t0,t00,avc,afc,asc,rn1,bnd1,rn2,bnd2,             &
      rn3,rn4,rn5,rn6,rn7,rn8,rn9,rn10,rn101,rn10a,rn11,rn11a,          &
      rn12,rn12a(31),rn12b(31),rn13(31),rn14,rn15,rn15a,rn16,rn17,      &
      rn17a,rn17b,rn17c,rn18,rn18a,rn19,rn19a,rn19b,rn20,rn20a,rn20b,   &
      bnd3,rn21,rn22,rn23,rn23a,rn23b,rn25,rn25a(31),rn30a,rn30b,       &
      rn30c,rn31,beta,rn32

  COMMON/size/ tnw,tns,tng,roqr,roqs,roqg
  COMMON/cont/ c76,c358,c172,c409,c218,c580,c610,c149,c879,c141
  COMMON/bterv/ zrc,zgc,zsc,vrc,vgc,vsc
  COMMON/b3cs/ ag,bg,as,bs,aww,bww,bgh,bgq,bsh,bsq,bwh,bwq
  COMMON/bsnw/ alv,alf,als,t0,t00,avc,afc,asc,rn1,bnd1,rn2,bnd2,        &
      rn3,rn4,rn5,rn6,rn7,rn8,rn9,rn10,rn101,rn10a,rn11,rn11a,          &
      rn12,rn12a,rn12b,rn13,rn14,rn15,rn15a,rn16,rn17,                  &
      rn17a,rn17b,rn17c,rn18,rn18a,rn19,rn19a,rn19b,rn20,rn20a,rn20b,   &
      bnd3,rn21,rn22,rn23,rn23a,rn23b,rn25,rn25a,rn30a,rn30b,           &
      rn30c,rn31,beta,rn32

  REAL :: r19bt,r20t,r23t,betah,r10t,r11at,ee2,a1,a2,ee1,bs6,x3,        &
       cpcgs,rt0,d2t,x1,x2,bw3,bsh5,bgh5,bw6,bs3,bg3,bwh5
!
!-----------------------------------------------------------------------
!
!  Above constants defined in subroutine STCSTICE
!
!  The followings are locally used variables
!
!-----------------------------------------------------------------------
!
  REAL :: dep
  REAL :: dd
  REAL :: dd1
  REAL :: dm
  REAL :: ern
  REAL :: rsub1
  REAL :: cnd
  REAL :: pgwet
  REAL :: egs
  REAL :: esi
  REAL :: qsi
  REAL :: ssi
  REAL :: qsw
  REAL :: pihom
  REAL :: ssw
  REAL :: pidw
  REAL :: pimlt
  REAL :: psmlt
  REAL :: pgmlt
  REAL :: psdep
  REAL :: pssub
  REAL :: pgsub
  REAL :: pint
  REAL :: pidep
  REAL :: prn
  REAL :: psn
  REAL :: dlt1

  REAL :: y1, y2, y3, y4, y5, y6

  INTEGER :: it, i,j,k
  REAL    :: t1, t2, t3, t4
  REAL    :: tem

  REAL    :: temp0, temp, interp, f1, f2, rstep, scvstep
  INTEGER :: nindex
  REAL    :: rhobar0

  REAL, POINTER :: qc(:,:,:,:)
  REAL, POINTER :: qr(:,:,:,:)
  REAL, POINTER :: qi(:,:,:,:)
  REAL, POINTER :: qs(:,:,:,:)
  REAL, POINTER :: qh(:,:,:,:)
!
!-----------------------------------------------------------------------
!
!  Define an inline function cvmgp(x1,x2,x3).
!
!  cvmgp=x1 for x3>=0.0
!  cvmgp=x2 for x3< 0.0
!
!-----------------------------------------------------------------------
!
  REAL :: eps
  REAL :: cvmgp
  cvmgp(x1,x2,x3)= x1*(0.5+SIGN(0.5,x3))+x2*(0.5-SIGN(0.5,x3))

! DTD: evaporation and melting rate array

  REAL :: mpteqnterms(nx,ny,nz,28)

  mpteqnterms = 0.0

!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

!-----------------------------------------------------------------------
!
! Since this subroutine is only called by microph_ice and it is also
! confirmed inside microph_ice that P_QC, P_QR, P_QI, P_QS, P_QH are
! defined, we do not have to make any changes.
!
!-----------------------------------------------------------------------

  qc => qscalar(:,:,:,:,P_QC)
  qr => qscalar(:,:,:,:,P_QR)
  qi => qscalar(:,:,:,:,P_QI)
  qs => qscalar(:,:,:,:,P_QS)
  qh => qscalar(:,:,:,:,P_QH)

!
!-----------------------------------------------------------------------
!
!  Two water and three classes of ice-phase.
!
!-----------------------------------------------------------------------
!
  it=1

  cpcgs=1.003E7    ! Missing in Tao's original code
                   ! Since ARPS has defined cp, better to use that.

  rt0=1./(t0-t00)

  bw3=bww+3.
  bs3=bs+3.
  bg3=bg+3.
  bwh5=2.5+bwh
  bsh5=2.5+bsh
  bgh5=2.5+bgh
  bw6=bww+6.
  bs6=bs+6.
  betah=.5*beta
  r10t=rn10*d2t
  r11at=rn11a*d2t
  r19bt=rn19b*d2t
  r20t=-rn20*d2t
  r23t=-rn23*d2t

  rhobar0 = 1.275*0.001     ! Air density in cgs at 0 C and 1000 mb
  t1 = SQRT(rhobar0)

  t2 = 1.e5 * zrc
  t3 = 1.e5 * zsc
  t4 = 1.e5 * zgc
!
!-----------------------------------------------------------------------
!
!  Begin k loop for vertical dimension
!
!-----------------------------------------------------------------------
!

!
!-----------------------------------------------------------------------
!
!  The inverse of the increment of lookup tables for rhobar*q is held
!  by rstep, and that for (rhobar/mu/psi**2) by scvstep, where mu is
!  dynamic viscosity of air and psi is diffusivity of water vapor in
!  air.
!
!-----------------------------------------------------------------------
!
  rstep   = 1.0 / 50.0E-10
  scvstep = 1.0 / 0.3

  DO k=2,nz-2

    DO j = 1, ny-1
      DO i = 1, nx-1
        tem2d1(i,j) = SQRT(rhobar(i,j,k))
        tem2d2(i,j) = t1 / tem2d1(i,j)
      END DO
    END DO

    DO j = 1,ny-1
      DO i = 1,nx-1
        ! The following clipping is already done just before the call to this subroutine,
        ! so the following should have no effect.
        qc(i,j,k,2) = MAX(0.0, qc(i,j,k,2))
        qr(i,j,k,2) = MAX(0.0, qr(i,j,k,2))
        qi(i,j,k,2) = MAX(0.0, qi(i,j,k,2))
        qs(i,j,k,2) = MAX(0.0, qs(i,j,k,2))
        qh(i,j,k,2) = MAX(0.0, qh(i,j,k,2))

        tair(i,j) = (ptprt(i,j,k,2)+ptbar(i,j,k))*ppi(i,j,k)
        tairc(i,j) = tair(i,j)-t0
!
!-----------------------------------------------------------------------
!
!  Compute zr,zs,zg,vr,vs,vg
!
!-----------------------------------------------------------------------
!
        zr(i,j) = t2/SQRT(tem2d1(i,j))
        zs(i,j) = t3/SQRT(tem2d1(i,j))
        zg(i,j) = t4/SQRT(tem2d1(i,j))
        vr(i,j) = 0.0
        vs(i,j) = 0.0
        vg(i,j) = 0.0

!
!      for rain we have a more complicated modification
!       vtr =  4854*sqrt(rho0/rho)* [gamma(5)*lambda**(4)/
!                                   gamma(4)*(lambda+f)**(5)
!

!        IF (qr(i,j,k,2) > 1.e-20) THEN
        IF (qr(i,j,k,2) > 1.e-16) THEN ! EMK
          dd = rhobar(i,j,k)*qr(i,j,k,2)
          y1 = SQRT(SQRT(dd))
          zr(i,j) = zrc/y1

          IF (fallopt == 2) THEN  ! Ferrier with updated coeff.

!            vr(i,j)=vrc*tem2d2(i,j)*(zr(i,j)**4)/((zr(i,j)+1.95)**5)

            temp  = min( 50.0e-6, max(0.0,dd) ) * rstep
            nindex = int(temp)
            f1 = pwrlam195ratio(nindex)
            f2 = pwrlam195ratio(nindex+1)
           vr(i,j)=vrc*tem2d2(i,j)*( f1 + (f2-f1)*(temp-nindex) )

          ELSE            ! original Lin scheme...

            temp  = min( 50.0e-6, max(0.0,dd) ) * rstep
            nindex = int(temp)
            f1 = pwr2(nindex)
            f2 = pwr2(nindex+1)
            vr(i,j) = vrc * tem2d2(i,j) * ( f1 + (f2-f1)*(temp-nindex) )

          END IF

        END IF

!        IF (qs(i,j,k,2) > 1.e-20) THEN
        IF (qs(i,j,k,2) > 1.e-16) THEN ! EMK
          dd = rhobar(i,j,k)*qs(i,j,k,2)
          y1 = SQRT(SQRT(dd))
          zs(i,j) = zsc/y1
          IF (fallopt == 2) THEN  ! Ferrier scheme with new coeff.
!            vs(i,j) = vsc*tem2d2(i,j)*(y1**.42)
!            vs(i,j) = vsc*tem2d2(i,j)*(dd**0.105) ! EMK

            temp  = min( 50.0e-6, max(0.0,dd) ) * rstep
            nindex = int(temp)
            f1 = pwr105(nindex)
            f2 = pwr105(nindex+1)
            vs(i,j) = vsc*tem2d2(i,j)*( f1 + (f2-f1)*(temp-nindex) ) ! EMK

          ELSE                  !  Lin scheme..
!            vs(i,j) = vsc*tem2d2(i,j)*sqrt(sqrt(y1))

            ! EMK...Adapted code from subroutine TERV
            temp  = min( 50.0e-6, max(0.0,dd) ) * rstep
            nindex = int(temp)
            f1 = pwr0625(nindex)
            f2 = pwr0625(nindex+1)
            vs(i,j) = vsc * tem2d2(i,j) * ( f1 + (f2-f1)*(temp-nindex) )

          END IF
        END IF

!  note zs is lambda for snow and vs is mass weighted mean fall
!  velocity for snow.

!        IF (qh(i,j,k,2) > 1.e-20) THEN
        IF (qh(i,j,k,2) > 1.e-16) THEN
          dd = rhobar(i,j,k)*qh(i,j,k,2)
          y1 = SQRT(SQRT(dd))
          zg(i,j) = zgc/y1

          IF (fallopt == 2) THEN  ! Ferrier scheme with new coeff.
!            vg(i,j) = vgc*tem2d2(i,j)*(y1**.6384)
!            vg(i,j) = vgc*tem2d2(i,j)*(dd**.1596) ! EMK

            temp  = min( 50.0e-6, max(0.0,dd) ) * rstep
            nindex = int(temp)
            f1 = pwr1596(nindex)
            f2 = pwr1596(nindex+1)
            vg(i,j) = vgc*tem2d2(i,j)*( f1 + (f2-f1)*(temp-nindex) )

          ELSE                  !  Lin scheme..
!            vg(i,j) = vgc/tem2d1(i,j)*sqrt(y1)

            ! EMK...Adapted code from subroutine TERV
            temp  = min( 50.0e-6, max(0.0,dd) ) * rstep
            nindex = int(temp)
            f1 = pwr0625(nindex)
            f2 = pwr0625(nindex+1)
            interp = f1 + (f2 - f1) * (temp - nindex)
            vg(i,j) = vgc / sqrt(rhobar(i,j,k)) * interp * interp

          END IF
        END IF
!
!-----------------------------------------------------------------------
!
!  y1  : dynamic viscosity of air (u)
!  dwv : diffusivity of water vapor in air (pi)
!  tca : thermal conductivity of air (ka)
!  y2  : kinetic viscosity (v)
!
!-----------------------------------------------------------------------
!
        y1 = c149*SQRT(tair(i,j)**3)/(tair(i,j)+120.)

        temp  = MIN( 150.0, MAX( 0.0, (tair(i,j)-173.15) ) )
        nindex = INT(temp)
        f1 = pwr81(nindex)
        f2 = pwr81(nindex+1)
        interp = f1 + (f2-f1)*(temp-nindex)

        dwv(i,j) = c879/pres(i,j,k)*tair(i,j) * interp
        tca(i,j) = c141*y1

        temp = MIN( 3000.0, rhobar(i,j,k)/(y1*dwv(i,j)**2) )            &
                * scvstep
        nindex = INT( temp )
        f1 = pwr1666(nindex)
        f2 = pwr1666(nindex+1)

        scv(i,j) = f1 + (f2-f1)*(temp-nindex)

      END DO
    END DO
!
!-----------------------------------------------------------------------
!
!   1 * psaut : autoconversion of qi to qs
!   3 * psaci : accretion of qi to qs
!   4 * psacw : accretion of qc by qs (riming) (qsacw for psmlt)
!   5 * praci : accretion of qi by qr
!   6 * piacr : accretion of qr or qg by qi
!
!-----------------------------------------------------------------------
!
    DO j = 1,ny-1
      DO i = 1,nx-1

        psaut(i,j) = 0.0
        psaci(i,j) = 0.0
        praci(i,j) = 0.0
        piacr(i,j) = 0.0
        psacw(i,j) = 0.0
        qsacw(i,j) = 0.0

        tem = 1.0/zs(i,j)
        dd = tem**3 * SQRT(SQRT(tem))  ! dd= (1/zs)**3.25

        IF(qr(i,j,k,2) > 1.0E-20) THEN
          temp0 = rhobar(i,j,k)*qr(i,j,k,2)
        ELSE
          temp0 = rhobar(i,j,k)*1.0E-20
        END IF
        temp  = MIN( 50.0E-6, temp0 ) * rstep
        nindex = INT(temp)
        f1 = pwr2(nindex)
        f2 = pwr2(nindex+1)
        interp = f1 + ( (f2 - f1) * (temp - nindex) )

        IF (tair(i,j) < t0) THEN
          esi = EXP(.025*tairc(i,j))
          psaut(i,j) = AMAX1(rn1*esi*(qi(i,j,k,2)-bnd1) ,0.0)
          psaci(i,j) = rn3*tem2d2(i,j)*esi*qi(i,j,k,2)*dd
          psacw(i,j) = rn4*tem2d2(i,j)*qc(i,j,k,2)*dd
          praci(i,j) = rn5*tem2d2(i,j)*qi(i,j,k,2)*interp/zr(i,j)**3
          piacr(i,j) = rn6*tem2d2(i,j)*qi(i,j,k,2)*interp               &
                          *(zr(i,j))**(-6)
!                          *(1.0/zr(i,j))**6
        ELSE
          qsacw(i,j) = rn4*tem2d2(i,j)*qc(i,j,k,2)*dd
        END IF
!
!-----------------------------------------------------------------------
!
!   praut : autoconversion of qc to qr                         (21)
!   pracw : accretion of qc by qr                              (22)
!
!-----------------------------------------------------------------------
!
        pracw(i,j) = rn22*tem2d2(i,j)*qc(i,j,k,2)*interp/zr(i,j)**3
        praut(i,j) = 0.0
        y1 = qc(i,j,k,2)-bnd3

! Original
!        IF (y1 > 0.0) praut(i,j) = rhobar(i,j,k)*y1*y1/(1.2E-4+rn21/y1)
! EMK 12 May 2003 Density correction...see equation 27 in
! Orville and Kopp, JAS, 1977.
        IF (y1 > 0.0) praut(i,j) = &
          rhobar(i,j,k)*y1*y1/(1.2E-4+rn21/(rhobar(i,j,k)*y1))
!
!-----------------------------------------------------------------------
!
!   psfw : bergeron processes for qs (koening, 1971)        (12)
!   psfi : bergeron processes for qs                        (13)
!
!-----------------------------------------------------------------------
!
        psfw(i,j) = 0.0
        psfi(i,j) = 0.0

        IF (tair(i,j) < t0) THEN
          y1 = AMAX1( AMIN1(tairc(i,j), -1.), -31.)
          it = INT(ABS(y1))
          y1 = rn12a(it)
          y2 = rn12b(it)
          psfw(i,j) = AMAX1(d2t*y1*(y2+rn12*rhobar(i,j,k)               &
                     *qc(i,j,k,2))*qi(i,j,k,2),0.0)
          psfi(i,j) = rn13(it)*qi(i,j,k,2)
        END IF
!
!-----------------------------------------------------------------------
!
!  qg = qg+min(pgdry,pgwet)
!  pgacs : accretion of qs by qg (dgacs,wgacs: dry and wet)    (9)
!  dgacw : accretion of qc by qg (qgacw for pgmlt)            (14)
!  dgacr : accretion of qr to qg (qgacr for pgmlt)            (16)
!
!-----------------------------------------------------------------------
!
        ee1 = 1.
        ee2 = 0.09
        egs = ee1*EXP(ee2*tairc(i,j))

        IF (tair(i,j) >= t0) egs = 1.0
        y1 = ABS(vg(i,j)-vs(i,j))
        y2 = zs(i,j)*zg(i,j)
        y3 = 5./y2
        y4 = .08*y3*y3
        y5 = .05*y3*y4
        dd = y1*(y3/zs(i,j)**5+y4/zs(i,j)**3+y5/zs(i,j))
        pgacs(i,j) = rn9/rhobar(i,j,k)*egs*dd
        dgacs(i,j) = pgacs(i,j)
        wgacs(i,j) = rn9/rhobar(i,j,k)*dd
        y1 = 1.0 / ( zg(i,j)**3 * SQRT(zg(i,j)) )
        dgacw(i,j) = AMAX1(rn14/tem2d1(i,j)*qc(i,j,k,2)*y1, 0.0)
        qgacw(i,j) = dgacw(i,j)
        y1 = ABS(vg(i,j)-vr(i,j))
        y2 = zr(i,j)*zg(i,j)
        y3 = 5./y2
        y4 = .08*y3*y3
        y5 = .05*y3*y4
        dd = rn16/rhobar(i,j,k)                                         &
                *y1*(y3/zr(i,j)**5+y4/zr(i,j)**3                        &
                +y5/zr(i,j))
        dgacr(i,j) = AMAX1(dd, 0.0)
        qgacr(i,j) = dgacr(i,j)

        IF (tair(i,j) >= t0) THEN
          dgacs(i,j) = 0.0
          wgacs(i,j) = 0.0
          dgacw(i,j) = 0.0
          dgacr(i,j) = 0.0
        ELSE
          pgacs(i,j) = 0.0
          qgacw(i,j) = 0.0
          qgacr(i,j) = 0.0
        END IF
!
!-----------------------------------------------------------------------
!
!  dgaci : accretion of qi by qg (wgaci for wet growth)       (15)
!  pgwet : wet growth of qg                                   (17)
!
!-----------------------------------------------------------------------
!
        dgaci(i,j) = 0.0
        wgaci(i,j) = 0.0
        pgwet = 0.0

        IF (tair(i,j) < t0) THEN

          y1 = qi(i,j,k,2)/( zg(i,j)**3 * SQRT(zg(i,j)) )
          dgaci(i,j) = rn15/tem2d1(i,j)*y1
          wgaci(i,j) = rn15a/tem2d1(i,j)*y1

!         y1 = 1./(alf+rn17c*tairc(i,j))

          IF(qh(i,j,k,2) > 1.0E-20) THEN
            temp0 = rhobar(i,j,k)*qh(i,j,k,2)
          ELSE
            temp0 = rhobar(i,j,k)*1.0E-20
          END IF
          temp  = MIN( 50.0E-6, temp0 ) * rstep
          nindex = INT(temp)
          f1 = pwr0625(nindex)
          f2 = pwr0625(nindex+1)
          interp = f1 + (f2 - f1) * (temp - nindex)

          y3 = .78/zg(i,j)**2+rn17a/SQRT(tem2d1(i,j))                   &
               *scv(i,j)/( zg(i,j)**3 * interp )
          y4 = rhobar(i,j,k)*alv*dwv(i,j)                               &
               *(3.799052E3/pres(i,j,k)-(qv(i,j,k,2)))                  &
               -tca(i,j)*tairc(i,j)
!mx
          y6 = alf+rn17c*tairc(i,j)
          y6 = sign(1.0, y6) * max(1.0e-20, abs(y6))

          dd = (rn17/rhobar(i,j,k)*y4*y3+(wgaci(i,j)+wgacs(i,j))*  &
               (alf+rn17b*tairc(i,j)))/y6
          pgwet = max(0.0, dd)

        END IF
!
!-----------------------------------------------------------------------
!
!  Shed process (wgacr = pgwet-dgacw-wgaci-wgacs)
!
!-----------------------------------------------------------------------
!

        wgacr(i,j) = pgwet-dgacw(i,j)-wgaci(i,j)-wgacs(i,j)
        y2 = dgacw(i,j)+dgaci(i,j)+dgacr(i,j)+dgacs(i,j)

        IF (pgwet >= y2) THEN
          wgacr(i,j) = 0.0
          wgaci(i,j) = 0.0
          wgacs(i,j) = 0.0
        ELSE
          dgacr(i,j) = 0.0
          dgaci(i,j) = 0.0
          dgacs(i,j) = 0.0
        END IF

      END DO
    END DO
!
!-----------------------------------------------------------------------
!
!  Handling the negative cloud water (qc)
!  Handling the negative cloud ice (qi)
!
!-----------------------------------------------------------------------
!
    eps = 1.0E-30

    DO j = 1,ny-1
      DO i = 1,nx-1
        y1 = (psacw(i,j)+praut(i,j)+pracw(i,j)+psfw(i,j)+               &
              dgacw(i,j)+qsacw(i,j)+qgacw(i,j))*d2t

        IF (qc(i,j,k,3) < y1) THEN
          a1 = qc(i,j,k,3)/(y1+eps)
          psacw(i,j) = psacw(i,j)*a1
          praut(i,j) = praut(i,j)*a1
          pracw(i,j) = pracw(i,j)*a1
          psfw (i,j) = psfw (i,j)*a1
          dgacw(i,j) = dgacw(i,j)*a1
          qsacw(i,j) = qsacw(i,j)*a1
          qgacw(i,j) = qgacw(i,j)*a1
          qc(i,j,k,3)  = 0.0
        ELSE
          qc(i,j,k,3)  = qc(i,j,k,3)-y1
        END IF

        y2 = (psaut(i,j)+psaci(i,j)+praci(i,j)+psfi(i,j)+               &
              dgaci(i,j)+wgaci(i,j))*d2t

        IF (y2 > qi(i,j,k,3)) THEN
          a2 = qi(i,j,k,3)/(y2+eps)
          psaut(i,j) = psaut(i,j)*a2
          psaci(i,j) = psaci(i,j)*a2
          praci(i,j) = praci(i,j)*a2
          psfi(i,j) = psfi(i,j)*a2
          dgaci(i,j) = dgaci(i,j)*a2
          wgaci(i,j) = wgaci(i,j)*a2
          qi(i,j,k,3) = 0.0
        ELSE
          qi(i,j,k,3) = qi(i,j,k,3)-y2
        END IF

        dlt3(i,j) = 0.0
        dlt2(i,j) = 0.0

        IF (tair(i,j) < t0) THEN

          IF (qr(i,j,k,2) < 1.e-4) THEN
            dlt3(i,j) = 1.0
            dlt2(i,j) = 1.0
          END IF

          IF (qs(i,j,k,2) >= 1.e-4) dlt2(i,j) = 0.0

        END IF

        pr(i,j) = (qsacw(i,j)+praut(i,j)+pracw(i,j)+                    &
                   qgacw(i,j))*d2t
        ps(i,j) = (psaut(i,j)+psaci(i,j)+psacw(i,j)+                    &
                   psfw(i,j)+psfi(i,j)+dlt3(i,j)*praci(i,j))*d2t
        pg(i,j) = ((1.-dlt3(i,j))*praci(i,j)+dgaci(i,j)+                &
                   wgaci(i,j)+dgacw(i,j))*d2t

      END DO
    END DO


!
!-----------------------------------------------------------------------
!
!  pracs : accretion of qs by qr                               (7)
!  psacr : accretion of qr by qs (qsacr for psmlt)             (8)
!
!-----------------------------------------------------------------------
!

    DO j = 1,ny-1
      DO i = 1,nx-1

        y1 = ABS(vr(i,j)-vs(i,j))
        y2 = zr(i,j)*zs(i,j)
        y3 = 5./y2
        y4 = .08*y3*y3
        y5 = .05*y3*y4
        pracs(i,j) = rn7/rhobar(i,j,k)                                  &
                   *y1*(y3/zs(i,j)**5+y4/zs(i,j)**3                     &
                   +y5/zs(i,j))
        psacr(i,j) = rn8/rhobar(i,j,k)                                  &
                   *y1*(y3/zr(i,j)**5+y4/zr(i,j)**3                     &
                   +y5/zr(i,j))
        qsacr(i,j) = psacr(i,j)

        IF (tair(i,j) >= t0) THEN
          pracs(i,j) = 0.0
          psacr(i,j) = 0.0
        ELSE
          qsacr(i,j) = 0.0
        END IF
!
!-----------------------------------------------------------------------
!
!  pgaut : autoconversion of qs to qg                          (2)
!  pgfr  : freezing of qr to qg                               (18)
!
!-----------------------------------------------------------------------
!
        pgaut(i,j) = 0.0
        pgfr (i,j) = 0.0

        IF (tair(i,j) < t0) THEN
!          y1 = EXP(.09*tairc(i,j))
!          pgaut(i,j) = AMAX1(rn2*y1*(qs(i,j,k,2)-bnd2),0.0)
          y2 = EXP(rn18a*(t0-tair(i,j)))
!          pgfr(i,j) = AMAX1(rn18/rhobar(i,j,k)*(y2-1.)                  &
!                       *(zr(i,j))**(-7), 0.0)
          IF (zr(i,j) > 4.0E05) THEN   ! It will be float underflow
            pgfr(i,j) = 0.0
          ELSE
            pgfr(i,j) = AMAX1(rn18/rhobar(i,j,k)*(y2-1.)                 &
                              *(zr(i,j))**(-7), 0.0)
          END IF

        END IF

      END DO
    END DO

!-----------------------------------------------------------------------
!
!  Handling the negative rain water (qr)
!  Handling the negative snow (qs)
!
!-----------------------------------------------------------------------

    DO j = 1,ny-1
      DO i = 1,nx-1

        y1 = (piacr(i,j)+dgacr(i,j)+wgacr(i,j)+psacr(i,j)+              &
              pgfr(i,j))*d2t

        tem = qr(i,j,k,3)+pr(i,j)
        IF (tem < y1) THEN
          a1 = tem/(y1+eps)
          piacr(i,j) = piacr(i,j)*a1
          dgacr(i,j) = dgacr(i,j)*a1
          wgacr(i,j) = wgacr(i,j)*a1
          pgfr (i,j) = pgfr (i,j)*a1
          psacr(i,j) = psacr(i,j)*a1
          qr(i,j,k,3)  = 0.0
        ELSE
          qr(i,j,k,3)  = qr(i,j,k,3)+pr(i,j)-y1
        END IF

        prn = d2t*((1.-dlt3(i,j))*piacr(i,j)+dgacr(i,j)+                &
              wgacr(i,j)+(1.-dlt2(i,j))*psacr(i,j)+pgfr(i,j))
        ps(i,j) = ps(i,j)+d2t*(dlt3(i,j)*piacr(i,j)+                    &
               dlt2(i,j)* psacr(i,j))
        pracs(i,j) = (1.-dlt2(i,j))*pracs(i,j)
        psn = d2t*(pgacs(i,j)+dgacs(i,j)+wgacs(i,j)+pgaut(i,j)+         &
               pracs(i,j))

        tem = qs(i,j,k,3)+ps(i,j)
        qs(i,j,k,3) = qs(i,j,k,3)+ps(i,j)-psn

        IF (qs(i,j,k,3) < 0.0) THEN

          IF ( psn > 0.0 ) THEN
            a2 = tem/psn
            pgacs(i,j) = pgacs(i,j)*a2
            dgacs(i,j) = dgacs(i,j)*a2
            wgacs(i,j) = wgacs(i,j)*a2
            pgaut(i,j) = pgaut(i,j)*a2
            pracs(i,j) = pracs(i,j)*a2
            psn = psn*a2
          ELSE
            psn = psn + qs(i,j,k,3)
          END IF

          qs(i,j,k,3) = 0.0

        END IF

        y2 = d2t*(psacw(i,j)+psfw(i,j)+dgacw(i,j)+piacr(i,j)+           &
              dgacr(i,j)+wgacr(i,j)+psacr(i,j)+pgfr(i,j))
        ptprt(i,j,k,3) = ptprt(i,j,k,3)+afc/ppi(i,j,k)*y2
        qh   (i,j,k,3) = qh(i,j,k,3)+pg(i,j)+prn+psn

        ! DTD

        mpteqnterms(i,j,k,20) = psacw(i,j)
        mpteqnterms(i,j,k,24) = piacr(i,j)
        mpteqnterms(i,j,k,22) = dgacw(i,j)
        mpteqnterms(i,j,k,27) = wgacr(i,j) + dgacr(i,j)  ! One of these will always be zero
        mpteqnterms(i,j,k,25) = psacr(i,j)
        mpteqnterms(i,j,k,23) = pgfr(i,j)

      END DO
    END DO


!
!-----------------------------------------------------------------------
!
!  psmlt : melting of qs                                      (11)
!  pgmlt : melting of qg to qr                                (19)
!
!-----------------------------------------------------------------------
!

!   write (*,*) "ZUWEN subsatopt/rhsat", subsatopt, rhsat

    DO j = 1,ny-1
      DO i = 1,nx-1

        psmlt = 0.0
        pgmlt = 0.0
!C        tair(i,j) = (ptprt(i,j,k,2)+ptbar(i,j,k))*ppi(i,j,k)

        IF (tair(i,j) >= t0) THEN

!C          tairc(i,j) = tair(i,j)-t0
          y1 = tca(i,j)*tairc(i,j)-rhobar(i,j,k)*alv                    &
               *dwv(i,j)*(3.799052E3/pres(i,j,k)-(qv(i,j,k,2)))

          IF(qs(i,j,k,2) > 1.0E-20) THEN
            temp0 = rhobar(i,j,k)*qs(i,j,k,2)
          ELSE
            temp0 = rhobar(i,j,k)*1.0E-20
          END IF
          temp  = MIN( 50.0E-6, temp0 ) * rstep
          nindex = INT(temp)
          f1 = pwr15625(nindex)
          f2 = pwr15625(nindex+1)
          interp = f1 + (f2 - f1) * (temp - nindex)

          y2 = .78/zs(i,j)**2+rn101*SQRT(tem2d2(i,j))                   &
               *scv(i,j)*interp/zs(i,j)**2
          dd = rn11/rhobar(i,j,k)*d2t*y1*y2                             &
               +r11at*tairc(i,j)*(qsacw(i,j)+qsacr(i,j))
          psmlt = AMAX1(0.0, AMIN1(dd, qs(i,j,k,3)))

          IF(qh(i,j,k,2) > 1.0E-20) THEN
            temp0 = rhobar(i,j,k)*qh(i,j,k,2)
          ELSE
            temp0 = rhobar(i,j,k)*1.0E-20
          END IF
          temp  = MIN( 50.0E-6, temp0 ) * rstep
          nindex = INT(temp)
          f1 = pwr0625(nindex)
          f2 = pwr0625(nindex+1)
          interp = f1 + (f2 - f1) * (temp - nindex)

          y3  = .78/zg(i,j)**2+rn19a/SQRT(tem2d1(i,j))                  &
               *scv(i,j)/( zg(i,j)**3 * interp )
          dd1 = rn19/rhobar(i,j,k)*d2t*y1*y3                            &
                +r19bt*tairc(i,j)*(qgacw(i,j)+qgacr(i,j))
          pgmlt = AMAX1(0.0, AMIN1(dd1, qh(i,j,k,3)))
          ptprt(i,j,k,3) = ptprt(i,j,k,3)-afc/ppi(i,j,k)*(psmlt+pgmlt)
          qr(i,j,k,3) = qr(i,j,k,3)+psmlt+pgmlt
          qs(i,j,k,3) = qs(i,j,k,3)-psmlt
          qh(i,j,k,3) = qh(i,j,k,3)-pgmlt

          ! DTD: store snow and graupel melting rates
          mpteqnterms(i,j,k,8) = psmlt/d2t
          mpteqnterms(i,j,k,10) = pgmlt/d2t

        END IF
!
!-----------------------------------------------------------------------
!
!  pihom : homogeneous freezing of qc to qi (t < t00)         (24)
!  pidw : deposition growth of qc to qi ( t0 < t < =  t00)    (25)
!  pimlt : melting of qi to qc (t > =  t0)                    (26)
!
!-----------------------------------------------------------------------

!C        IF (qc(i,j,k,2).le.1.e-20) qc(i,j,k,2) = 0.0
!C        IF (qi(i,j,k,2).le.1.e-20) qi(i,j,k,2) = 0.0

        tair3(i,j) = (ptprt(i,j,k,3)+ptbar(i,j,k))*ppi(i,j,k)
        pihom = cvmgp(0.,qc(i,j,k,3),tair3(i,j)-t00)
        pimlt = cvmgp(qi(i,j,k,3),0.,tair3(i,j)-t0)
        pidw = 0.0

        IF (tair(i,j) < t0 .AND. tair(i,j) > t00) THEN
!C          tairc(i,j) = tair(i,j)-t0
          y1 = AMAX1( AMIN1(tairc(i,j), -1.), -31.)
          it = INT(ABS(y1))
          y2 = rn25a(it)
          y3 = EXP(.5*ABS(tairc(i,j)))
          pidw = AMIN1(rn25/rhobar(i,j,k)*d2t*y2*                       &
               y3, qc(i,j,k,3))
        END IF

        y1 = pihom-pimlt+pidw
        ptprt(i,j,k,3) = ptprt(i,j,k,3)+afc/ppi(i,j,k)*y1
        qc(i,j,k,3) = qc(i,j,k,3)-y1
        qi(i,j,k,3) = qi(i,j,k,3)+y1

        ! DTD: store melting rate of qi to qc
        mpteqnterms(i,j,k,7) = pimlt/d2t
        ! store freezing rate of cloud to ice
        ! store depositional growth of ice at expense of cloud
        mpteqnterms(i,j,k,18) = (pihom+pidw)/d2t

!
!-----------------------------------------------------------------------
!
!  pint  : initiation of qi                                   (31)
!  pidep : deposition of qi                                   (32)
!
!  Note:   rhsat is the minimum RH for condensation to occur.
!          It is used to adjust saturated specific humidity,
!          qsw and qsi, in the following code. Zuwen, 04/2002
!
!-----------------------------------------------------------------------
!
        pint = 0.0
!C        tair(i,j) = (ptprt(i,j,k,2)+ptbar(i,j,k))*ppi(i,j,k)

        IF (tair(i,j) < t0) THEN
!C          tairc(i,j) = tair(i,j)-t0
          dd = rn31/rhobar(i,j,k)*EXP(beta*tairc(i,j))
          rtair(i,j) = 1./(tair(i,j)-c76)
          y2 = EXP(c218-c580*rtair(i,j))
          qsi = rhsat*(3.799052E3/pres(i,j,k))*y2
          esi = c610*y2
          ssi = (qv(i,j,k,2))/qsi-1.
          dm = AMAX1( (qv(i,j,k,3)-qsi), 0.)
          rsub1 = c580*asc*qsi*rtair(i,j)*rtair(i,j)
          pint = AMIN1(dd,dm)
          y1 = 1./tair(i,j)
          y2 = EXP(betah*tairc(i,j))
          IF(qi(i,j,k,2) > 1.0E-30) THEN
            y3 = SQRT(qi(i,j,k,2))
          ELSE
            y3 = 1.0E-15
          END IF
          dd = y1*(rn30a*y1-rn30b)+rn30c*                               &
               tair(i,j)/esi
          pidep = AMAX1(rn32*d2t/tem2d1(i,j)                            &
               *ssi*y2*y3/dd, 0.0)
          pint = pint+pidep
          dep = dm/(1.+rsub1)
          pint = AMIN1(pint,dep)
          ptprt(i,j,k,3) = ptprt(i,j,k,3)+asc/ppi(i,j,k)*pint
          qv(i,j,k,3) = qv(i,j,k,3)-pint
          qi(i,j,k,3) = qi(i,j,k,3)+pint

          ! DTD store initiation rate (includes initial deposition) of ice
          mpteqnterms(i,j,k,13) = pint/d2t

        END IF

      END DO
    END DO
!
!-----------------------------------------------------------------------
!
!  Tao et al (1989) saturation technique
!
!-----------------------------------------------------------------------
!
    DO j = 1,ny-1
      DO i = 1,nx-1

        tair3(i,j) = (ptprt(i,j,k,3)+ptbar(i,j,k))*ppi(i,j,k)
        cnd = rt0*(tair3(i,j)-t00)
        dep = rt0*(t0-tair3(i,j))
        y1 = 1./(tair3(i,j)-c358)
        y2 = 1./(tair3(i,j)-c76)
        qsw = rhsat*(3.799052E3/pres(i,j,k))*EXP(c172-c409*y1)
        qsi = rhsat*(3.799052E3/pres(i,j,k))*EXP(c218-c580*y2)
        dd = c409*ppi(i,j,k)*y1*y1
        dd1 = c580*ppi(i,j,k)*y2*y2

        IF (qc(i,j,k,3) <= 1.e-20) qc(i,j,k,3) = 1.e-20
        IF (qi(i,j,k,3) <= 1.e-20) qi(i,j,k,3) = 1.e-20

        IF (tair3(i,j) >= t0) THEN
          dep = 0.0
          cnd = 1.0
          qi(i,j,k,3) = 0.0
        END IF

        IF (tair3(i,j) <= t00) THEN
          cnd = 0.0
          dep = 1.0
          qc(i,j,k,3) = 0.0
        END IF

        y5 = avc/ppi(i,j,k)*cnd+asc/ppi(i,j,k)*dep
        y1 = qc(i,j,k,3)*qsw/(qc(i,j,k,3)+qi(i,j,k,3))
        y2 = qi(i,j,k,3)*qsi/(qc(i,j,k,3)+qi(i,j,k,3))
        y4 = dd*y1+dd1*y2
        dm = qv(i,j,k,3)-y1-y2
        rsub1 = dm/(1.+y4*y5)
        cnd = cnd*rsub1
        dep = dep*rsub1

        IF (qc(i,j,k,3) <= 1.e-20) qc(i,j,k,3) = 0.
        IF (qi(i,j,k,3) <= 1.e-20) qi(i,j,k,3) = 0.
!
!-----------------------------------------------------------------------
!
!  Condensation or evaporation of qc
!
!-----------------------------------------------------------------------
!
        cnd = AMAX1(-qc(i,j,k,3),cnd)
!
!-----------------------------------------------------------------------
!
!  Deposition or sublimation of qi
!
!-----------------------------------------------------------------------
!
        dep = AMAX1(-qi(i,j,k,3),dep)
        ptprt(i,j,k,3) = ptprt(i,j,k,3)+avc/ppi(i,j,k)*                 &
                     cnd+asc/ppi(i,j,k)*dep
        qv(i,j,k,3) = qv(i,j,k,3)-cnd-dep
        qc(i,j,k,3) = qc(i,j,k,3)+cnd
        qi(i,j,k,3) = qi(i,j,k,3)+dep

        !DTD: store evaporation of qc and sublimation of qi
        IF(cnd < 0.0) THEN
          mpteqnterms(i,j,k,1) = -cnd/d2t
        ELSE
          mpteqnterms(i,j,k,11) = cnd/d2t
        END IF

        IF(dep < 0.0) THEN
          mpteqnterms(i,j,k,3) = -dep/d2t
        ELSE
          mpteqnterms(i,j,k,14) = dep/d2t
        END IF

      END DO
    END DO
!
!-----------------------------------------------------------------------
!
!  psdep : deposition or sublimation of qs                    (10)
!  pgsub : sublimation of qg                                  (20)
!
!-----------------------------------------------------------------------
!
    DO j = 1,ny-1
      DO i = 1,nx-1

        psdep = 0.0
        pssub = 0.0
        pgsub = 0.0
!C        tair(i,j) = (ptprt(i,j,k,2)+ptbar(i,j,k))*ppi(i,j,k)

        IF (tair(i,j) < t0) THEN

          IF (qs(i,j,k,2) < 1.e-20) qs(i,j,k,2) = 0.0
          IF (qh(i,j,k,2) < 1.e-20) qh(i,j,k,2) = 0.0

          rtair(i,j) = 1./(tair(i,j)-c76)
          qsi = rhsat*(3.799052E3/pres(i,j,k))*EXP(c218-c580*rtair(i,j))
          ssi = (qv(i,j,k,2))/qsi-1.
          y1 = rn10a*rhobar(i,j,k)/(tca(i,j)*tair(i,j)**2)              &
                +1./(dwv(i,j)*qsi)

          IF(qs(i,j,k,2) > 1.0E-20) THEN
            temp0 = rhobar(i,j,k)*qs(i,j,k,2)
          ELSE
            temp0 = rhobar(i,j,k)*1.0E-20
          END IF
          temp  = MIN( 50.0E-6, temp0 ) * rstep
          nindex = INT(temp)
          f1 = pwr15625(nindex)
          f2 = pwr15625(nindex+1)
          interp = f1 + (f2 - f1) * (temp - nindex)

          y2 = .78/zs(i,j)**2+rn101*SQRT(tem2d2(i,j))                   &
               *scv(i,j)*interp/zs(i,j)**2
          psdep = r10t*ssi*y2/y1
          pssub = psdep
          psdep = AMAX1(psdep, 0.)
          pssub = AMAX1(-qs(i,j,k,3), AMIN1(pssub, 0.))

          IF(qh(i,j,k,2) > 1.0E-20) THEN
            temp0 = rhobar(i,j,k)*qh(i,j,k,2)
          ELSE
            temp0 = rhobar(i,j,k)*1.0E-20
          END IF
          temp  = MIN( 50.0E-6, temp0 ) * rstep
          nindex = INT(temp)
          f1 = pwr0625(nindex)
          f2 = pwr0625(nindex+1)
          interp = f1 + (f2 - f1) * (temp - nindex)

          y2 = .78/zg(i,j)**2+rn20b/SQRT(tem2d1(i,j))                   &
                *scv(i,j)/( zg(i,j)**3 * interp )
          pgsub = r20t*ssi*y2/y1
          dm = qv(i,j,k,2)-qsi
          rsub1 = c580*asc*qsi*rtair(i,j)*rtair(i,j)
!
!-----------------------------------------------------------------------
!
!  Deposition or sublimation of qs
!
!-----------------------------------------------------------------------
!
          y1 = dm/(1.+rsub1)
          psdep = AMIN1(psdep,AMAX1(y1,0.))
          y2 = AMIN1(y1,0.)
          pssub = AMAX1(pssub,y2)
!
!-----------------------------------------------------------------------
!
!  Sublimation of qg
!
!-----------------------------------------------------------------------
!
          dd = AMAX1((-y2-qs(i,j,k,3)), 0.)
          pgsub = AMIN1(dd, qh(i,j,k,3), AMAX1(pgsub,0.))
          dlt1 = cvmgp(1.,0.,qc(i,j,k,2)+qi(i,j,k,2)-1.e-8)
          psdep = dlt1*psdep
          pssub = (1.-dlt1)*pssub
          pgsub = (1.-dlt1)*pgsub
          ptprt(i,j,k,3) = ptprt(i,j,k,3)+asc/ppi(i,j,k)                &
                *(psdep+ pssub-pgsub)
          qv(i,j,k,3) = qv(i,j,k,3)+pgsub-pssub-psdep
          qs(i,j,k,3) = qs(i,j,k,3)+psdep+pssub
          qh(i,j,k,3) = qh(i,j,k,3)-pgsub

          ! DTD: store deposition and sublimation rate of qs and qh
          mpteqnterms(i,j,k,4) = -pssub/d2t
          IF(pgsub > 0.0) THEN
            mpteqnterms(i,j,k,6) = pgsub/d2t
          ELSE
            mpteqnterms(i,j,k,17) = -pgsub/d2t
          END IF

          mpteqnterms(i,j,k,15) = psdep/d2t

        END IF

!
!-----------------------------------------------------------------------
!
!* 23 * ern : evaporation of qr (subsaturation)                    (23**
!
!-----------------------------------------------------------------------
!
        ern = 0.0

        IF (qr(i,j,k,2) > 1.e-20) THEN

!C          tair(i,j) = (ptprt(i,j,k,2)+ptbar(i,j,k))*ppi(i,j,k)
          rtair(i,j) = 1./(tair(i,j)-c358)
          qsw = rhsat*(3.799052E3/pres(i,j,k))*EXP(c172-c409*rtair(i,j))
          ssw = (qv(i,j,k,2))/qsw-1.0
          dm = qv(i,j,k,3)-qsw
          rsub1 = c409*avc*qsw*rtair(i,j)*rtair(i,j)
          dd1 = AMAX1(-dm/(1.+rsub1), 0.0)
!wdt update: (see arpssupport emails by Vince Wong & Eric Kemp
!"Re: micro_ice3d.f (ARPS Version 4.5.2.)" 9 Jul 2001)
!DTD update (06/02/2009), added switch for calculation of y1 (bulk
!ventilation coefficient for the two different fall speed
!formulations. Previously, in the final term, the Lin fall coefficients
!were being used even when the Ferrier option (fallopt == 2) was chosen.

          IF(fallopt == 2) THEN  ! Ferrier fall speed option
            y1 = .78/zr(i,j)**2+rn23a*SQRT(tem2d2(i,j))                   &
                    *scv(i,j)/(zr(i,j)+1.95)**3.0
          ELSE
          y1 = .78/zr(i,j)**2+rn23a*SQRT(tem2d2(i,j))                   &
                  *scv(i,j)/zr(i,j)**2.9
          END IF

          y2 = rn23b*rhobar(i,j,k)/(tca(i,j)*tair(i,j)**2)              &
                +1./(dwv(i,j)*qsw)
          ern = r23t*ssw*y1/y2
          ern = AMIN1(dd1,qr(i,j,k,3),AMAX1(ern,0.))
          ptprt(i,j,k,3) = ptprt(i,j,k,3)-avc/ppi(i,j,k)*ern
          qv(i,j,k,3) = qv(i,j,k,3)+ern
          qr(i,j,k,3) = qr(i,j,k,3)-ern

          ! DTD: store evaporation rate of rain
          IF(ern > 0.0) THEN
            mpteqnterms(i,j,k,2) = ern/d2t
          ELSE
            mpteqnterms(i,j,k,12) = -ern/d2t
          END IF

        END IF

        qc(i,j,k,3) = MAX(0., qc(i,j,k,3))
        qs(i,j,k,3) = MAX(0., qs(i,j,k,3))
        qi(i,j,k,3) = MAX(0., qi(i,j,k,3))
        qh(i,j,k,3) = MAX(0., qh(i,j,k,3))
        qr(i,j,k,3) = MAX(0., qr(i,j,k,3))

      END DO
    END DO

  END DO

  RETURN
END SUBROUTINE icecvt

!##################################################################
!##################################################################
!######                                                      ######
!######              SUBROUTINE MICROPH_ICE                  ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE microph_ice(nx,ny,nz,mphyflg,dtbig1,                         &
           ptprt,pprt,qv,qscalar,raing,prcrate,                         &
           rhostr,pbar,ptbar,qvbar,ppi,j3,j3inv,                        &
           rhobar,cgsrhobar,cgspres,tem1,tem2,tem3,tem4,tem5,           &
           tem6,tem7,tem8,tem9,tem10,tem11,tem12,mpteqnterms)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate and apply the microphysical contributions to the water,
!  ice and temperature fields, using an ice microphysics parameterization
!  scheme.
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
!  9/10/92 (M. Xue)
!  Negative water contents are not zeroed out in this version.
!
!  09/15/94 (M. Xue)
!  Modified to adapt Tao's ice code.
!
!  10/20/94 (Yuhe Liu)
!  Fixed a bug in the argument list of calling subroutine ICECVT.
!
!  03/05/97 (Fanyou Kong -- CMRP)
!  Modify the code to apply all three time levels in calculate
!  production terms in Tao scheme (icecvt)
!
!  07/10/97 (Fanyou Kong - CMRP)
!  Include MPDCD advection option (sadvopt = 5)
!
!  11/18/98 (Keith Brewster)
!  Changed pibar to ppi (full pi).
!
!  08/27/02 (Dan Weber)
!  Added option for using Ferrier (1994) fall velocity coefficients
!  and replaces the Lin constants when fallopt=2.
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
!    dtbig1   The large time step size for this call.
!
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
!    j3       Coordinate transformation Jacobian  d(zp)/dz
!    exner
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
!    cgsrhobarBase state density in cgs unit.
!    cgspres  Base state pressure in cgs unit.
!    tem1     Temporary work array.
!    tem2     Temporary work array.
!    tem3     Temporary work array.
!    tem4     Temporary work array.
!    tem5     Temporary work array.
!    tem6     Temporary work array.
!    tem7     Temporary work array.
!    tem8     Temporary work array.
!    tem9     Temporary work array.
!    tem10    Temporary work array.
!    tem11    Temporary work array.
!    tem12    Temporary work array.
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
  INTEGER :: irsh         ! Flag identifying hydrometer fields
                          ! 0 = rain; 1 = snow; 2 = hail

  INTEGER :: nx,ny,nz          ! Number of grid points in 3 directions
  INTEGER :: mphyflg           ! microphysic flag

  REAL :: dtbig1               ! The large time step size for this call.

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
  REAL :: j3inv (nx,ny,nz)     ! Coordinate transformation Jacobian defined as
                               ! d( zp )/d( z ).
  REAL :: qvbar (nx,ny,nz)     ! Base state water vapor specific humidity
                               ! (kg/kg)
  REAL :: ppi   (nx,ny,nz)     ! Exner function.

  ! DTD: added array to store terms in temperature equation

  REAL :: mpteqnterms(nx,ny,nz,28)

!
!-----------------------------------------------------------------------
!
!  Temporary arrays
!
!-----------------------------------------------------------------------
!
  REAL :: rhobar(nx,ny,nz)     ! Temporary work array
  REAL :: cgsrhobar(nx,ny,nz)  ! Temporary work array
  REAL :: cgspres(nx,ny,nz)    ! Temporary work array

  REAL :: tem1  (nx,ny,nz)     ! Temporary work array
  REAL :: tem2  (nx,ny,nz)     ! Temporary work array
  REAL :: tem3  (nx,ny,nz)     ! Temporary work array
  REAL :: tem4  (nx,ny,nz)     ! Temporary work array
  REAL :: tem5  (nx,ny,nz)     ! Temporary work array
  REAL :: tem6  (nx,ny,nz)     ! Temporary work array
  REAL :: tem7  (nx,ny,nz)     ! Temporary work array
  REAL :: tem8  (nx,ny,nz)     ! Temporary work array
  REAL :: tem9  (nx,ny,nz)     ! Temporary work array
  REAL :: tem10 (nx,ny,nz)     ! Temporary work array
  REAL :: tem11 (nx,ny,nz)     ! Temporary work array
  REAL :: tem12 (nx,ny,nz)     ! Temporary work array
!
!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k,lvlq
!
!-----------------------------------------------------------------------
!
!  Following constants are defined in subroutine STCSTICE
!
!  (Suggest to move the following COMMON block to a inlcude file.)
!
!-----------------------------------------------------------------------
!
  REAL :: tnw,tns,tng,roqr,roqs,roqg
  REAL :: c76,c358,c172,c409,c218,c580,c610,c149,c879,c141
  REAL :: zrc,zgc,zsc,vrc,vgc,vsc
  REAL :: ag,bg,as,bs,aww,bww,bgh,bgq,bsh,bsq,bwh,bwq
  REAL :: alv,alf,als,t0,t00,avc,afc,asc,rn1,bnd1,rn2,bnd2,             &
          rn3,rn4,rn5,rn6,rn7,rn8,rn9,rn10,rn101,rn10a,rn11,rn11a,      &
          rn12,rn12a(31),rn12b(31),rn13(31),rn14,rn15,rn15a,rn16,rn17,  &
          rn17a,rn17b,rn17c,rn18,rn18a,rn19,rn19a,rn19b,rn20,rn20a,     &
          rn20b,bnd3,rn21,rn22,rn23,rn23a,rn23b,rn25,rn25a(31),rn30a,   &
          rn30b,rn30c,rn31,beta,rn32

  COMMON/size/ tnw,tns,tng,roqr,roqs,roqg
  COMMON/cont/ c76,c358,c172,c409,c218,c580,c610,c149,c879,c141
  COMMON/bterv/ zrc,zgc,zsc,vrc,vgc,vsc
  COMMON/b3cs/ ag,bg,as,bs,aww,bww,bgh,bgq,bsh,bsq,bwh,bwq
  COMMON/bsnw/ alv,alf,als,t0,t00,avc,afc,asc,rn1,bnd1,rn2,bnd2,        &
       rn3,rn4,rn5,rn6,rn7,rn8,rn9,rn10,rn101,rn10a,rn11,rn11a,         &
       rn12,rn12a,rn12b,rn13,rn14,rn15,rn15a,rn16,rn17,                 &
       rn17a,rn17b,rn17c,rn18,rn18a,rn19,rn19a,rn19b,rn20,rn20a,        &
       rn20b,bnd3,rn21,rn22,rn23,rn23a,rn23b,rn25,rn25a,rn30a,          &
       rn30b,rn30c,rn31,beta,rn32

  INTEGER :: constset
  SAVE constset
  DATA constset /0/
  REAL :: denwater,tem,tema, deltat

  INTERFACE
    SUBROUTINE icecvt(nx,ny,nz, d2t,                                    &
               ptprt,qv,qscalar,rhobar,ptbar,qvbar,pres,ppi,            &
               wgacr, scv, tca, dwv, zr, vr, zs, vs,                    &
               zg, vg, psaut, psaci,psacw, qsacw,praci,piacr,           &
               praut, pracw, psfw, psfi, dgacs, dgacw, dgaci,           &
               dgacr, pgacs, wgacs, qgacw,wgaci,qgacr,pgaut,            &
               pracs, psacr, qsacr, pgfr, tair, tairc,                  &
               tair3, tairc3,                                           &
               pr, pg, ps, dlt2, dlt3, rtair, tem2d1, tem2d2, mpteqnterms)

      IMPLICIT NONE

      INCLUDE 'phycst.inc'
      INCLUDE 'globcst.inc'
      INCLUDE 'timelvls.inc'

      INTEGER :: nx, ny, nz
      REAL    :: d2t

      REAL :: ptprt (nx,ny,nz,nt)  ! Perturbation potential temperature (K)
      REAL :: qv    (nx,ny,nz,nt)  ! Water vapor specific humidity (g/g)

      REAL, TARGET :: qscalar(nx,ny,nz,nt,nscalar)

      REAL :: rhobar(nx,ny,nz)  ! Base state air density (g/cm**3)
      REAL :: pres  (nx,ny,nz)  ! Pressure ((g cm/s**2)/cm**2)
      REAL :: ptbar (nx,ny,nz)  ! Base state potential temperature (K)
      REAL :: qvbar (nx,ny,nz)  ! Base state air density (g/g)
      REAL :: ppi   (nx,ny,nz)  ! Exner function
    !-----------------------------------------------------------------------
    !

      REAL :: wgacr(nx,ny)
      REAL :: psaut(nx,ny)
      REAL :: psaci(nx,ny)
      REAL :: psacw(nx,ny)
      REAL :: qsacw(nx,ny)
      REAL :: praci(nx,ny)
      REAL :: piacr(nx,ny)
      REAL :: praut(nx,ny)
      REAL :: pracw(nx,ny)
      REAL :: psfw (nx,ny)
      REAL :: psfi (nx,ny)
      REAL :: dgacs(nx,ny)
      REAL :: dgacw(nx,ny)
      REAL :: dgaci(nx,ny)
      REAL :: dgacr(nx,ny)
      REAL :: pgacs(nx,ny)
      REAL :: wgacs(nx,ny)
      REAL :: qgacw(nx,ny)
      REAL :: wgaci(nx,ny)
      REAL :: qgacr(nx,ny)
      REAL :: pgaut(nx,ny)
      REAL :: pracs(nx,ny)
      REAL :: psacr(nx,ny)
      REAL :: qsacr(nx,ny)
      REAL :: pgfr (nx,ny)

      REAL :: tair (nx,ny)
      REAL :: tairc(nx,ny)
      REAL :: tair3 (nx,ny)
      REAL :: tairc3(nx,ny)
      REAL :: pr   (nx,ny)
      REAL :: pg   (nx,ny)
      REAL :: ps   (nx,ny)
      REAL :: dlt2 (nx,ny)
      REAL :: dlt3 (nx,ny)
      REAL :: rtair(nx,ny)

      REAL :: scv  (nx,ny)
      REAL :: tca  (nx,ny)
      REAL :: dwv  (nx,ny)
      REAL :: zr   (nx,ny)
      REAL :: vr   (nx,ny)
      REAL :: zs   (nx,ny)
      REAL :: vs   (nx,ny)
      REAL :: zg   (nx,ny)
      REAL :: vg   (nx,ny)

      REAL :: tem2d1(nx,ny)
      REAL :: tem2d2(nx,ny)

      REAL :: mpteqnterms(nx,ny,nz,28)

    END SUBROUTINE icecvt
  END INTERFACE
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

  ! DTD
  mpteqnterms = 0.0

  lvlq = tfuture

!-----------------------------------------------------------------------
!
! Since this scheme comtain microphysical ice process, it must at least
! define qc, qr, qi, qs, qh
!
!-----------------------------------------------------------------------

  IF (P_QC < 1 .OR. P_QR < 1 .OR. P_QI < 1 .OR. P_QS < 1 .OR. P_QH < 1) THEN

    WRITE(6,'(2a,/,5(a,I2),/,a)')                                       &
               'No enough microphysical array was defined ',            &
               'inside subroutine microph_ice.',                        &
               'P_QC = ',P_QC,' P_QR = ',P_QR,' P_QI = ',P_QI,          &
              ' P_QS = ',P_QS,' P_QH = ',P_QH,                          &
               'Program aborting ...'
    CALL arpsstop('Wrong size for microphysics array, qscalar.',1)

  END IF
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
  IF (constset == 0) THEN

    CALL stcstice
    constset = 1

  END IF

  DO k = 1,nz-1
    DO j = 1,ny-1
      DO i = 1,nx-1
        qv(i,j,k,tfuture) = MAX(0.0,qv(i,j,k,tfuture))
        qscalar(i,j,k,tfuture,P_QC) = MAX(0.0,qscalar(i,j,k,tfuture,P_QC))
        qscalar(i,j,k,tfuture,P_QR) = MAX(0.0,qscalar(i,j,k,tfuture,P_QR))
        qscalar(i,j,k,tfuture,P_QI) = MAX(0.0,qscalar(i,j,k,tfuture,P_QI))
        qscalar(i,j,k,tfuture,P_QS) = MAX(0.0,qscalar(i,j,k,tfuture,P_QS))
        qscalar(i,j,k,tfuture,P_QH) = MAX(0.0,qscalar(i,j,k,tfuture,P_QH))
      END DO
    END DO
  END DO

  DO k = 1, nz-1
    DO j = 1, ny-1
      DO i = 1, nx-1
        rhobar(i,j,k) = rhostr(i,j,k) *j3inv(i,j,k)
      END DO
    END DO
  END DO

  IF(mphyflg == 1) THEN

    DO k = 1, nz-1
      DO j = 1, ny-1
        DO i = 1, nx-1
          cgsrhobar(i,j,k) = rhobar(i,j,k) * 0.001
          cgspres  (i,j,k) = (pprt(i,j,k,tfuture)+pbar(i,j,k))*10.0
        END DO
      END DO
    END DO
!
!-----------------------------------------------------------------------
!
!  Call ice parameterization routine that does conversions among
!  water and ice categories.
!
!-----------------------------------------------------------------------
!
!C    IF( sadvopt.ge.1.and.sadvopt.le.3) THEN ! Leapfrog scheme
    IF( sadvopt /= 4 .and. tintegopt == 1) THEN                  ! Leapfrog scheme
      deltat = 2*dtbig1
    ELSE                                    ! Forward scheme
      deltat = dtbig1
    END IF

    CALL icecvt(nx, ny, nz, deltat,                                     &
              ptprt,qv,qscalar,cgsrhobar, ptbar,qvbar,cgspres, ppi,     &
              tem1(1,1,1), tem1(1,1,2), tem1(1,1,3), tem1(1,1,4),       &
              tem2(1,1,1), tem2(1,1,2), tem2(1,1,3), tem2(1,1,4),       &
              tem3(1,1,1), tem3(1,1,2), tem3(1,1,3), tem3(1,1,4),       &
              tem4(1,1,1), tem4(1,1,2), tem4(1,1,3), tem4(1,1,4),       &
              tem5(1,1,1), tem5(1,1,2), tem5(1,1,3), tem5(1,1,4),       &
              tem6(1,1,1), tem6(1,1,2), tem6(1,1,3), tem6(1,1,4),       &
              tem7(1,1,1), tem7(1,1,2), tem7(1,1,3), tem7(1,1,4),       &
              tem8(1,1,1), tem8(1,1,2), tem8(1,1,3), tem8(1,1,4),       &
              tem9(1,1,1), tem9(1,1,2), tem9(1,1,3), tem9(1,1,4),       &
              tem10(1,1,1),tem10(1,1,2),tem10(1,1,3),tem10(1,1,4),      &
              tem11(1,1,1),tem11(1,1,2),tem11(1,1,3),tem11(1,1,4),      &
              tem12(1,1,1),tem12(1,1,2),mpteqnterms)

    DO k = 1,nz-1
      DO j = 1,ny-1
        DO i = 1,nx-1
          qv(i,j,k,tfuture) = MAX(0.0,qv(i,j,k,tfuture))
          qscalar(i,j,k,tfuture,P_QC) = MAX(0.0,qscalar(i,j,k,tfuture,P_QC))
          qscalar(i,j,k,tfuture,P_QR) = MAX(0.0,qscalar(i,j,k,tfuture,P_QR))
          qscalar(i,j,k,tfuture,P_QI) = MAX(0.0,qscalar(i,j,k,tfuture,P_QI))
          qscalar(i,j,k,tfuture,P_QS) = MAX(0.0,qscalar(i,j,k,tfuture,P_QS))
          qscalar(i,j,k,tfuture,P_QH) = MAX(0.0,qscalar(i,j,k,tfuture,P_QH))
        END DO
      END DO
    END DO

  END IF   ! mphyflg == 1
!
!-----------------------------------------------------------------------
!
!  Hydrometeor sedimentation
!
!-----------------------------------------------------------------------
!

  irsh = 0  ! for rainwater
  CALL qhfall (irsh,dtbig1,nx,ny,nz, qscalar(:,:,:,:,P_QR),             &
               rhobar,j3,j3inv,                                         &
               tem4(1,1,1), tem1,tem2,tem5,tem6,tem7,tem8,tem9)

  irsh = 1  ! for snow
  CALL qhfall (irsh,dtbig1,nx,ny,nz, qscalar(:,:,:,:,P_QS),             &
               rhobar,j3,j3inv,                                         &
               tem4(1,1,2), tem1,tem2,tem5,tem6,tem7,tem8,tem9)

  irsh = 2  ! for hail
  CALL qhfall (irsh,dtbig1,nx,ny,nz, qscalar(:,:,:,:,P_QH),             &
               rhobar,j3,j3inv,                                         &
               tem4(1,1,3), tem1,tem2,tem5,tem6,tem7,tem8,tem9)

  denwater = 1000.0   ! Density of liquid water (kg/m**3)
  tem = denwater/(1000.0*dtbig)

  DO j=1,ny-1
    DO i=1,nx-1
      tema = tem4(i,j,1)+tem4(i,j,2)+tem4(i,j,3)
      raing  (i,j) = raing(i,j)+tema
      prcrate(i,j) = tema*tem
    END DO
  END DO

  DO k = 1,nz-1
    DO j = 1,ny-1
      DO i = 1,nx-1
        qscalar(i,j,k,tfuture,P_QR) = MAX(0.0,qscalar(i,j,k,tfuture,P_QR))
        qscalar(i,j,k,tfuture,P_QS) = MAX(0.0,qscalar(i,j,k,tfuture,P_QS))
        qscalar(i,j,k,tfuture,P_QH) = MAX(0.0,qscalar(i,j,k,tfuture,P_QH))
      END DO
    END DO
  END DO

  RETURN
END SUBROUTINE microph_ice

!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE QHFALL                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE qhfall(irsh,dtbig1,nx,ny,nz,q,rhobar,j3,j3inv,               &
           draing, vtr3d,tem1,tem2,tem3,tem4,tem5,tem6)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate the fall-out of hydrometers.
!  This subroutine is called by the ice microphysics package.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  10/10/1991. Written based on QRFALL used by warmrain microphysics.
!
!  MODIFICATION HISTORY:
!
!  10/20/1996 (M. Xue)
!  This routine now calls a standard routine QFALLOUT.
!  Precipitation rate array prcrate is now correctly calculated
!  for split-step integration.

!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    irsh     Flag identifying hydrometer fields
!    dtbig1   The large time step size for this call.
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!    q        Hydrometer mixing ratio at all time levels (kg/kg)
!
!    rhobar   Base state air density (kg/m**3)
!    j3       Coordinate transformation Jacobian  d(zp)/dz
!    j3inv    1/j3.
!
!  OUTPUT:
!
!    q        Hydrometer mixing ratio at time tfuture (kg/kg)
!    draing   Grid-scale rainfall (mm)
!
!  WORK ARRAYS:
!
!    vtr3d    Work array
!    tem1     Work array
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
  INCLUDE 'timelvls.inc'
!
  INTEGER :: irsh         ! Flag identifying hydrometer fields
                          ! 0 = rain, 1 = snow, 2 = hail
  REAL :: dtbig1          ! The large time step size for this call.
  INTEGER :: nx,ny,nz     ! Number of grid points in 3 directions

  REAL :: q     (nx,ny,nz,nt)  ! Rainwater mixing ratio (kg/kg)
  REAL :: rhobar(nx,ny,nz)     ! Base state air density (kg/m**3)
  REAL :: j3    (nx,ny,nz)     ! Coordinate transformation Jacobian
                               ! defined as d( zp )/d( z ).
  REAL :: j3inv (nx,ny,nz)     ! 1/j3.
  REAL :: draing (nx,ny)       ! Grid-scale rainfall (mm)
!
!-----------------------------------------------------------------------
!
!  Temporary arrays
!
!-----------------------------------------------------------------------
!
  REAL :: vtr3d(nx,ny,nz)
  REAL :: tem1 (nx,ny,nz)
  REAL :: tem2 (nx,ny,nz)
  REAL :: tem3 (nx,ny,nz)
  REAL :: tem4 (nx,ny,nz)
  REAL :: tem5 (nx,ny,nz)
  REAL :: tem6 (nx,ny,nz)
!
!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  DO k = 1, nz-1
    DO j = 1, ny-1
      DO i = 1, nx-1
        tem1(i,j,k) = rhobar(i,j,k) * 0.001
      END DO
    END DO
  END DO
!
!-----------------------------------------------------------------------
!
!  Calculate terminal falling speed for q, and store in tem1
!
!-----------------------------------------------------------------------
!
  CALL terv(nx,ny,nz,irsh, tem1, q(1,1,1,tpresent), vtr3d)

!-----------------------------------------------------------------------
!
!  vtr3d is defined at the w point and is in cgs unit.
!
!-----------------------------------------------------------------------

  DO k=1,nz-1
    DO j=1,ny-1
      DO i=1,nx-1
        vtr3d(i,j,k)=vtr3d(i,j,k)*0.01
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

  CALL qfallout(nx,ny,nz,dtbig1,q, rhobar,j3,j3inv, draing,vtr3d,       &
                tem1,tem2,tem3,tem4,tem5,tem6)

  RETURN
END SUBROUTINE qhfall
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINE TERV                      ######
!######                                                      ######
!######                     Developed by                     ######
!######                                                      ######
!######    Goddard Cumulus Ensemble Modeling Group, NASA     ######
!######                                                      ######
!######     Center for Analysis and Prediction of Storms     ######
!######               University of Oklahoma                 ######
!######                                                      ######
!##################################################################
!##################################################################
!
!

SUBROUTINE terv(nx,ny,nz,irsg,rhobar,q,vtr)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Compute the terminal velocities (vtr) of rain water qr, snow qs and
!  hail qg.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Goddard Cumulus Ensemble Modeling Group, NASA
!
!  MODIFICATION HISTORY:
!
!  10/20/1996 (M. Xue)
!  Rewritten for ARPS.
!
!  02/24/1997 (J. Zong, M. Xue and Yuhe Liu)
!  Power calculations are replaced by lookup table functions for
!  terminal velocity.
!
!  08/27/2002 (D. Weber)
!  Added fallopt option to provide a choice of fall velocity
!  computations.
!
!-----------------------------------------------------------------------
!

!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!  irsg    Flag identifying hydrometeor field as rain, snow or hail
!          =  0 for rain;  qpt is rain
!          =  1 for snow;  qpt is snow
!          =  2 for hail;  qpt is hail
!  q       Hydrometeor field defined at the scalar point
!  rhobar  Air density defined at scalar point (g/cm**3)
!
!  OUTPUT:
!
!  vtr     Vertical velocity defined at the scalar point (cm/s)
!
!-----------------------------------------------------------------------
!
  INTEGER :: nx,ny,nz
  REAL :: q     (nx,ny,nz)
  REAL :: rhobar(nx,ny,nz)
  REAL :: vtr   (nx,ny,nz)
  REAL :: tema, temb, temc, temd, rho0cgs
  REAL :: temp, interp, f1, f2, rstep
  INTEGER :: nindex
!
!-----------------------------------------------------------------------
!
!  Common variables. Defined in subroutine STCSTICE
!
!-----------------------------------------------------------------------
!
  COMMON/bterv/ zrc,zgc,zsc,vrc,vgc,vsc
  COMMON/b3cs/ ag,bg,as,bs,aww,bww,bgh,bgq,bsh,bsq,bwh,bwq

  INCLUDE 'globcst.inc'
  INCLUDE 'phycst.inc'

!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  rho0cgs = rho0*0.001

  rstep = 1.0/50.0E-10

  temc = 0.0

  IF (irsg == 0) THEN    ! irsg = 0    for rain water (qr)

    DO k = 1,nz-1
      DO j = 1,ny-1
        DO i = 1,nx-1
          tema=SQRT(rho0cgs/rhobar(i,j,k))
          temb=rhobar(i,j,k)*q(i,j,k)

          IF (temb > 1.e-16) THEN

            IF (fallopt == 2) THEN  ! ferrier numbers.....

               temd=zrc/(sqrt(sqrt(temb)))
!              vtr(i,j,k)=vrc*tema*(temd**4)/((temd+1.95)**5)

              temp  = min( 50.0e-6, max(0.0,temb) ) * rstep
              nindex = int(temp)
              f1 = pwrlam195ratio(nindex)
              f2 = pwrlam195ratio(nindex+1)

              vtr(i,j,k)=vrc*tema*(f1 + (f2-f1)*(temp-nindex))

            ELSE   !  original Lin scheme

              temp  = min( 50.0e-6, max(0.0,temb) ) * rstep
              nindex = int(temp)
              f1 = pwr2(nindex)
              f2 = pwr2(nindex+1)

              vtr(i,j,k) = vrc * tema * (f1 + (f2-f1)*(temp-nindex))

            END IF

          ELSE
            vtr(i,j,k) = 0.0
          END IF
          temc = max(temc,vtr(i,j,k))
        END DO
      END DO
    END DO

!   print *,'max rain fall speed is=',temc,irsg

  ELSE IF( irsg == 1) THEN  !irsg = 1 for snow (qs)

    DO k = 1,nz-1
      DO j = 1,ny-1
        DO i = 1,nx-1
          tema=SQRT(rho0cgs/rhobar(i,j,k))
          temb=rhobar(i,j,k)*q(i,j,k)

          IF (temb > 1.e-16) THEN

            IF (fallopt == 2) THEN  ! Ferrier numbers.....
!              vtr(i,j,k) = vsc * tema * (temb**.105)
              vtr(i,j,k) = vsc * tema * (temb**.105)

              temp  = min( 50.0e-6, max(0.0,temb) ) * rstep
              nindex = int(temp)
              f1 = pwr105(nindex)
              f2 = pwr105(nindex+1)
              vtr(i,j,k) = vsc * tema * ( f1 + (f2-f1)*(temp-nindex) )

            ELSE          ! original Lin scheme.......
              temp  = min( 50.0e-6, max(0.0,temb) ) * rstep
              nindex = int(temp)
              f1 = pwr0625(nindex)
              f2 = pwr0625(nindex+1)
              vtr(i,j,k) = vsc * tema * ( f1 + (f2-f1)*(temp-nindex) )
            END IF

          ELSE
            vtr(i,j,k) = 0.0
          END IF

          temc = max(temc,vtr(i,j,k))
        END DO
      END DO
    END DO

!   print *,'max snow fall speed is=',temc,irsg

  ELSE IF( irsg == 2) THEN  ! irsg = 2   for graupel (qg)

    DO k = 1,nz-1
      DO j = 1,ny-1
        DO i = 1,nx-1
          tema=sqrt(rho0cgs/rhobar(i,j,k))
          temb=rhobar(i,j,k)*q(i,j,k)
          IF (temb > 1.e-16) THEN

            IF (fallopt == 2) THEN  ! new method

!              vtr(i,j,k) = vgc * tema * (temb**.1596)

              temp  = min( 50.0e-6, max(0.0,temb) ) * rstep
              nindex = int(temp)
              f1 = pwr1596(nindex)
              f2 = pwr1596(nindex+1)
              vtr(i,j,k) = vgc * tema * ( f1 + (f2-f1)*(temp-nindex) )

            ELSE    ! Original method
              temp  = min( 50.0e-6, max(0.0,temb) ) * rstep
              nindex = int(temp)
              f1 = pwr0625(nindex)
              f2 = pwr0625(nindex+1)
              interp = f1 + (f2 - f1) * (temp - nindex)
              vtr(i,j,k) = vgc / sqrt(rhobar(i,j,k)) * interp * interp
            END IF
            temc = max(temc,vtr(i,j,k))

          ELSE
            vtr(i,j,k) = 0.0
          END IF
            temc = max(temc,vtr(i,j,k))
        END DO
      END DO
    END DO

!   print *,'max hail fall speed is=',temc,irsg

  END IF

  RETURN
END SUBROUTINE terv
