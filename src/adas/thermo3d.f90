!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE ARPS_BE                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE arps_be(nx,ny,nz,                                            &
           pres,hgt,tk,wmr,                                             &
           lcl,lfc,el,twdf,li,cape,mcape,cin,tcap,                      &
           p1d,ht1d,t1d,tv1d,td1d,wmr1d,                                &
           partem,buoy,wload,mbuoy,pbesnd,mbesnd,tem1)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate the lifting condensation level (lcl), level of free
!  convection (lfc), equilibrium level (el), max wet-bulb potential
!  temperature difference (twdf), lifted index (LI), Convective
!  Available Potential Energy (CAPE), Moist Convective Potential
!  Energy (MCAPE, includes water loading), convective inhibition
!  (CIN) and lid strength (tcap)  over the ARPS grid.
!
!-----------------------------------------------------------------------
!
!
!  AUTHOR: Keith Brewster
!  April, 1995   Based on OLAPS, hence LAPS, version of same.
!
!  MODIFICATION HISTORY:
!  3/11/1996  (Keith Brewster)
!  Cleaned-up, removed OLAPS artifacts.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
  INTEGER :: nx,ny,nz
!
!-----------------------------------------------------------------------
!
!  3-D input variables
!
!-----------------------------------------------------------------------
!
  REAL :: pres(nx,ny,nz)
  REAL :: hgt(nx,ny,nz)
  REAL :: tk(nx,ny,nz)
  REAL :: wmr(nx,ny,nz)
!
!-----------------------------------------------------------------------
!
!  Output variables (2d)
!
!-----------------------------------------------------------------------
!
  REAL :: lcl(nx,ny)
  REAL :: lfc(nx,ny)
  REAL :: el(nx,ny)
  REAL :: twdf(nx,ny)
  REAL :: li(nx,ny)
  REAL :: cape(nx,ny)
  REAL :: mcape(nx,ny)
  REAL :: cin(nx,ny)
  REAL :: tcap(nx,ny)
!
!-----------------------------------------------------------------------
!
!  Scratch space for calculations
!
!-----------------------------------------------------------------------
!
  REAL :: p1d(nz),ht1d(nz)
  REAL :: t1d(nz),tv1d(nz),td1d(nz),wmr1d(nz)
  REAL :: partem(nz),buoy(nz),wload(nz)
  REAL :: mbuoy(nz),pbesnd(nz),mbesnd(nz)
  REAL :: tem1(nx,ny)

!-----------------------------------------------------------------------
!
! Inlucde files
!
!-----------------------------------------------------------------------
  INCLUDE 'mp.inc'
!
!-----------------------------------------------------------------------
!
!  Functions
!
!-----------------------------------------------------------------------
!
  REAL :: wmr2td,tctotv
!
!-----------------------------------------------------------------------
!
!  Misc internal variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k,n,imid,jmid,nlevel
  INTEGER :: imidproc, jmidproc
  INTEGER :: nxlg, nylg
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  nxlg = (nx-3)*nproc_x  + 3
  nylg = (ny-3)*nproc_y  + 3

  imid = nxlg/2     ! Added by Y. Wang for proper diagnostic outputs
  jmid = nylg/2
  imidproc = (imid-2) / (nx-3) + 1
  jmidproc = (jmid-2) / (ny-3) + 1

  IF (loc_x == imidproc) THEN
    imid = MOD((imid-2),(nx-3)) + 2   ! Local index for global middle point
  ELSE
    imid = -999
  END IF

  IF (loc_y == jmidproc) THEN
    jmid = MOD((jmid-2),(ny-3)) + 2   ! Local index for global middle point
  ELSE
    jmid = -999
  END IF

!
!-----------------------------------------------------------------------
!
!  Loop over all columns
!
!-----------------------------------------------------------------------
!
  DO j=1,ny-1
    DO i=1,nx-1
!
      DO k=2,nz-1
        n = k-1
        p1d(n)  = 0.01 * pres(i,j,k)
        ht1d(n) = hgt(i,j,k)*1000.
        wmr1d(n)= 1000. * wmr(i,j,k)
        t1d(n)  = tk(i,j,k) - 273.15
        td1d(n) = wmr2td(p1d(n),wmr1d(n))
        tv1d(n) = tctotv(t1d(n),wmr1d(n))
      END DO
!
      nlevel = nz-2
!
      CALL sindex(nz,nlevel,p1d,ht1d,t1d,tv1d,td1d,wmr1d,               &
                  partem,buoy,wload,mbuoy,pbesnd,mbesnd,                &
                  lcl(i,j),lfc(i,j),el(i,j),twdf(i,j),                  &
                  li(i,j),cape(i,j),mcape(i,j),                         &
                  cin(i,j),tcap(i,j))


      IF(i == imid .AND. j == jmid) THEN
        WRITE(6,'(1x,2(a,I4),a)')                                       &
          '=== Diagnostic outputs from arps_be for imid = ',imid,       &
          ', jmid = ',jmid,' ==='
        WRITE(6,*)                                                      &
            ' n    p         t      tpar      buoy      wmr       be'
        DO n = 1,nlevel
          WRITE(6,801) n,p1d(n),t1d(n),partem(n),buoy(n),               &
                      wmr1d(n),pbesnd(n)
          801       FORMAT(1X,i2,f8.1,f10.2,f10.2,f10.4,f10.4,f10.1)
        END DO
      END IF

    END DO
  END DO
!
!-----------------------------------------------------------------------
!
!  Smooth the output arrays
!
!-----------------------------------------------------------------------
!
  CALL smooth9p(  lcl,nx,ny,1,nx-1,1,ny-1,0,tem1)
  CALL smooth9p(  lfc,nx,ny,1,nx-1,1,ny-1,0,tem1)
  CALL smooth9p(   el,nx,ny,1,nx-1,1,ny-1,0,tem1)
  CALL smooth9p( twdf,nx,ny,1,nx-1,1,ny-1,0,tem1)
  CALL smooth9p(   li,nx,ny,1,nx-1,1,ny-1,0,tem1)
  CALL smooth9p( cape,nx,ny,1,nx-1,1,ny-1,0,tem1)
  CALL smooth9p(mcape,nx,ny,1,nx-1,1,ny-1,0,tem1)
  CALL smooth9p(  cin,nx,ny,1,nx-1,1,ny-1,0,tem1)
  CALL smooth9p( tcap,nx,ny,1,nx-1,1,ny-1,0,tem1)
!
  RETURN
END SUBROUTINE arps_be
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE SINDEX                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE sindex(maxlev,nlevel,p,ht,t,tv,td,w,                         &
           partem,buoy,wload,mbuoy,pbesnd,mbesnd,                       &
           lcl_pbe,lfc,el,twdf,li,cape,mcape,cin,tcap)
!
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate Convective Available Potential Energy (CAPE)
!  in each column of ARPS grid.
!
!-----------------------------------------------------------------------
!
!
!  AUTHOR: Keith Brewster
!  April, 1995   Based on OLAPS, hence LAPS, version of same.
!
!  MODIFICATION HISTORY:
!  3/11/1996  Cleaned-up, removed OLAPS artifacts.
!  5/13/1996  Added cap strength.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
  INTEGER :: maxlev,nlevel
!
  REAL :: p(maxlev),ht(maxlev),t(maxlev),tv(maxlev),                    &
       td(maxlev),w(maxlev)
  REAL :: partem(maxlev),buoy(maxlev),wload(maxlev),mbuoy(maxlev)
  REAL :: pbesnd(maxlev),mbesnd(maxlev)
!
!  Returned from sindx
!
  REAL :: lfc,el,twdf,li,cape,mcape,cin,mcin,tcap
!
!  Potbe variables
!
  REAL :: plcl_pbe,tlcl_pbe,lcl_pbe,thepcl
  REAL :: velneg,mvelneg
!
!  Functions
!
  REAL :: oe
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!
  thepcl=oe(t(1),td(1),p(1))
!
!  print *, ' theta-e of parcel: ',thepcl
!
  CALL ptlcl(p(1),t(1),td(1),plcl_pbe,tlcl_pbe)
!
!  Find height of LCL
!
  CALL intrpr(maxlev,nlevel,p,ht,plcl_pbe,lcl_pbe)
!
!  Calculate the CAPE and such
!
  CALL potbe(maxlev,nlevel,p(1),t(1),w(1),                              &
             thepcl,plcl_pbe,tlcl_pbe,lcl_pbe,                          &
             p,ht,t,tv,td,w,                                            &
             partem,buoy,wload,mbuoy,pbesnd,mbesnd,                     &
             cin,velneg,cape,mcin,mvelneg,mcape,lfc,el,tcap)
!
!  Calculate Lifted Index
!
  CALL calcli(maxlev,nlevel,thepcl,p,t,li)
!
!  Calculate max and min wet bulb potential temperature
!
  CALL thwxn(maxlev,nlevel,p,ht,t,td,ht(1),twdf)
!
  RETURN
END SUBROUTINE sindex
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE POTBE                      ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
!

SUBROUTINE potbe(maxlev,nlevel,pmean,tmean,wmean,                       &
           blthte,plcl,tlcl,lcl,                                        &
           p,ht,t,tv,td,w,                                              &
           partem,buoy,wload,mbuoy,pbesnd,mbesnd,                       &
           pbeneg,velneg,pos_max,                                       &
           mbeneg,mvelneg,mpos_max,lfc,el,tcap)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate Convective Available Potential Energy (CAPE)
!  in each column of ARPS grid.
!
!-----------------------------------------------------------------------
!
!
!  AUTHOR: Keith Brewster
!  February, 1994  Based on OLAPS, hence LAPS, version of same.
!                  from FSL, by Steve Albers 1991
!
!  MODIFICATION HISTORY:
!  3/11/1996  Cleaned-up, removed OLAPS artifacts.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
  INTEGER :: maxlev,nlevel
  REAL :: pmean,tmean,wmean,blthte,plcl,tlcl,lcl
  REAL :: p(maxlev),ht(maxlev)
  REAL :: t(maxlev),tv(maxlev),td(maxlev),w(maxlev)
  REAL :: partem(maxlev),buoy(maxlev),wload(maxlev),mbuoy(maxlev)
  REAL :: pbesnd(maxlev),mbesnd(maxlev)
  REAL :: pbeneg,velneg,pos_max
  REAL :: mbeneg,mvelneg,mpos_max,lfc,el,tcap
!
!  Parameters
!
  REAL :: g,gamma
  PARAMETER (g=9.80665,                                                 &
             gamma = .009760)   ! Dry Adiabatic Lapse Rate Deg/m
!
!  Functions
!
  REAL :: tsa_fast,tctotv
!
!  Misc internal variables
!
  INTEGER :: n,nel
  REAL :: deltah,delta_ht_dry,delta_ht_wet
  REAL :: sntlcl,buoy_lcl,wsat,partv
  REAL :: nbe_min,pbe_wet,pbe_dry,pos_area
  REAL :: wlow,htzero,adjeng
!
!-----------------------------------------------------------------------
!
!  Function f_mrsat and inline directive for Cray PVP
!
!-----------------------------------------------------------------------
!
  REAL :: f_mrsat

!fpp$ expand (f_mrsat)
!!dir$ inline always f_mrsat
!*$*  inline routine (f_mrsat)

!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!
!  Reset output variables.
!
!  These should be the same as what is assigned
!  no positive area found.  They are limited in
!  range to allow for contouring when positive
!  areas exist in some columns of a domain and
!  not in others.
!
  pbeneg=-400.
  velneg=20.
  pos_max=0.
  mbeneg=-400.
  mvelneg=20.
  mpos_max=0.
  lfc=10000.
  el=0.
  tcap=0.
!
!  Inxtialize parcel path arrays
!
  partem(1) = t(1)
  buoy(1) = 0.
  wload(1) = 0.
  mbuoy(1) = 0.
  pbesnd(1) = 0.
  mbesnd(1) = 0.

!  WRITE(6,810)pmean,tmean,wmean,plcl,tlcl,lcl
! 810 format(' pmean,tmean,wmean,plcl,tlcl,lcl',2F10.2,F10.5,2F10.2
!    +   ,F5.1)

  DO n=2,nlevel
    deltah = ht(n) - ht(n-1)
    IF(plcl < p(n-1))THEN ! lower level is below LCL
      IF(plcl < p(n))THEN ! upper level is below LCL
!        WRITE(6,*)' DRY CASE'
        partem(n)=partem(n-1)-gamma*deltah
        partv=tctotv(partem(n),w(1))
        buoy(n)=(partv-tv(n))/tv(n)
        pbesnd(n)=pbesnd(n-1)+g*0.5*(buoy(n)+buoy(n-1))*deltah
        wload(n)=0.
        mbuoy(n)=buoy(n)
        mbesnd(n)=pbesnd(n)
        IF((p(1)-p(n)) < 300.) tcap=AMAX1(tcap,(tv(n)-partv))

      ELSE ! Upper level is above LCL
!
!  BRACKETING CASE AROUND lcl - DRY ADIABATIC PART
!
!        WRITE(6,*)' DRY ADIABATIC PART'
        delta_ht_dry = lcl - ht(n-1)
!        WRITE(6,307)tlcl
!307        format(' PARCEL TEMP AT lcl= ',F10.3)
        CALL intrpr(maxlev,nlevel,p,tv,plcl,sntlcl)
        partv=tctotv(tlcl,w(1))
        buoy_lcl=(partv-sntlcl)/sntlcl
        pbe_dry=g*0.5*(buoy_lcl+buoy(n-1))*delta_ht_dry
        IF((p(1)-plcl) < 300.) tcap=AMAX1(tcap,(sntlcl-partv))
!        WRITE(6,777)N,P(N),tlcl,sntlcl,buoy_lcl
!#          ,buoy(N-1),delta_ht_dry,HT(N),pbesnd(N-1)+pbe_dry
!
!        MOIST ADIABATIC PART
!
!        WRITE(6,*)' MOIST ADIABATIC PART'
        delta_ht_wet=deltah-delta_ht_dry

        partem(n) = tsa_fast(blthte,p(n))
        wsat=1000.*f_mrsat( p(n)*100., partem(n)+273.15 )
        partv=tctotv(partem(n),wsat)
        buoy(n)=(partv-tv(n))/tv(n)
        pbe_wet = g*0.5*(buoy(n)+buoy_lcl)*delta_ht_wet
        pbesnd(n)=pbesnd(n-1) + pbe_dry + pbe_wet
!
        wload(n)=0.001*(w(1)-wsat)
        mbuoy(n)=buoy(n) - wload(n)
        pbe_wet = g*0.5*(mbuoy(n)+buoy_lcl)*delta_ht_wet
        mbesnd(n)=mbesnd(n-1) + pbe_dry + pbe_wet
        IF((p(1)-plcl) < 300.) tcap=AMAX1(tcap,(tv(n)-partv))

      END IF ! Upper level below LCL (Dry or bracket)
    ELSE ! Lower Level is above LCL
!      WRITE(6,*)' GETTING PARCEL TEMPERATURE FOR MOIST CASE'
      partem(n) = tsa_fast(blthte,p(n))

      wsat=1000.*f_mrsat( p(n)*100., partem(n)+273.15 )
      partv=tctotv(partem(n),wsat)
      buoy(n)=(partv-tv(n))/tv(n)
      pbesnd(n)=pbesnd(n-1)+g*0.5*(buoy(n)+buoy(n-1))*deltah
!
      wload(n)=0.001*(w(1)-wsat)
      mbuoy(n)=buoy(n) - wload(n)
      mbesnd(n)=mbesnd(n-1)+g*0.5*(mbuoy(n)+mbuoy(n-1))*deltah
      IF((p(1)-p(n)) < 300.) tcap=AMAX1(tcap,(tv(n)-partv))

    END IF

!    WRITE(6,777)N,P(N),partem(N),T(N),(buoy(n)*1000.),pbesnd(n)
!777    format(' PBE: P,partem,t,b,pbe=',I3,F6.1,4F8.2)
  END DO
!
!  DETERMINE ENERGY EXTREMA
!  Find heights with nuetral buoyancy
!
  pos_area=0.
  nbe_min=0.
  DO n=2,nlevel
!    WRITE(6,940)N
!940    format(
!    :' LOOKING FOR NEUTRAL BUOYANCY - ENERGY EXTREMUM, LEVEL',I3)

    IF((buoy(n)*buoy(n-1)) < 0.)THEN
      wlow=buoy(n)/(buoy(n)-buoy(n-1))
      htzero=ht(n)*(1.-wlow) + wlow*ht(n-1)
      deltah=htzero-ht(n-1)
      adjeng=pbesnd(n-1)+g*0.5*buoy(n-1)*deltah
!
      IF (p(n) >= 500.)  THEN
        nbe_min=AMIN1(adjeng,nbe_min)
      END IF
!
      pos_area=adjeng-nbe_min
      pos_max=AMAX1(pos_area,pos_max)
    END IF
  END DO

!  WRITE(6,464)ICP,ICT,N1,NLEVEL
!464  format(' ICP,ICT,N1,NLEVEL',4I5)
!
!  Case when equlibrium level is above top of domain
!
  pos_area=pbesnd(nlevel)-nbe_min
  pos_max=AMAX1(pos_area,pos_max)
!
!  At least one region of positive area in sounding
!  Make sure there is at least 1 J/kg to avoid some
!  round-off errors esp near LCL.
!
  IF(pos_max > 1.0)THEN
    pbeneg=AMAX1(nbe_min,-400.)
    velneg=SQRT(2.0*ABS(pbeneg))
    velneg=AMIN1(velneg,20.)
  ELSE ! Case when no positive area exists anywhere in sounding
    pos_max=0.0
    pbeneg =-400.
    velneg = 20.
  END IF
!  WRITE(6,485)pos_max,PBENEG,VELNEG
!485  format(' pos_max',F10.1,' PBENEG',F10.1,' VELNEG',F10.1)

!
!  DETERMINE ENERGY EXTREMA FOR MOIST BUOYANCY
!  Find heights with nuetral buoyancy
!
  pos_area=0.
  nbe_min=0.
  DO n=2,nlevel
!    WRITE(6,940)N

    IF((mbuoy(n)*mbuoy(n-1)) < 0.)THEN
      wlow=mbuoy(n)/(mbuoy(n)-mbuoy(n-1))
      htzero=ht(n)*(1.-wlow) + wlow*ht(n-1)
      deltah=htzero-ht(n-1)
      adjeng=mbesnd(n-1)+g*0.5*mbuoy(n-1)*deltah
!
      IF (p(n) >= 500.)  THEN
        nbe_min=AMIN1(adjeng,nbe_min)
      END IF
!
      pos_area=adjeng-nbe_min
      mpos_max=AMAX1(pos_area,mpos_max)
    END IF
  END DO

!  WRITE(6,464)ICP,ICT,N1,NLEVEL
!
!  Case when equlibrium level is above top of domain
!
  pos_area=mbesnd(nlevel)-nbe_min
  mpos_max=AMAX1(pos_area,mpos_max)
!
!  At least one region of positive area in sounding
!  Make sure there is at least 1 J/kg to
!  spurious pos energy due to round off.
!
  IF(mpos_max > 1.0)THEN
    mbeneg=AMAX1(nbe_min,-400.)
    mvelneg=SQRT(2.0*ABS(pbeneg))
    mvelneg=AMIN1(mvelneg,20.)
  ELSE ! Case when no positive area exists anywhere in sounding
    mpos_max=0.0
    mbeneg =-400.
    mvelneg = 20.
  END IF
!  WRITE(6,486)mpos_max,PBENEG,VELNEG
!486  format(' Mpos_max',F10.1,' MBENEG',F10.1,' mVELNEG',F10.1)
!
!    Case when equlibrium level is above top of domain
!
  mpos_max = MAX(mpos_max,(mbesnd(nlevel) - nbe_min))
!
!  Find EL and LFC
!  Unxts are set to km ASL
!
  IF(pos_max > 1.0) THEN
    IF(buoy(nlevel) > 0.) THEN
      nel=nlevel
      el=0.001*ht(nlevel)
    ELSE
      DO  n=nlevel-1,2,-1
        IF(buoy(n) > 0.) EXIT
      END DO
!      1201     CONTINUE
      nel=n
      wlow=buoy(n+1)/(buoy(n+1)-buoy(n))
      el=0.001 * (ht(n+1)*(1.-wlow) + ht(n)*wlow)
    END IF
!
    DO n=nel,1,-1
      IF(buoy(n) < 0.) EXIT
    END DO
!    1301   CONTINUE
    IF(n > 0) THEN
      wlow=buoy(n+1)/(buoy(n+1)-buoy(n))
      lfc=ht(n+1)*(1.-wlow) + ht(n)*wlow
    ELSE
      lfc=ht(1)
    END IF
  ELSE
    el=0.
    lfc=10000.
  END IF
  lfc=AMIN1(lfc,10000.)
  RETURN
END SUBROUTINE potbe
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                 FUNCTION TCTOTV                      ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
!

  FUNCTION tctotv(tt,ww)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Virtual Temperature
!
!  Given T in Celcius and mixing ratio in g/kg
!  find the virtual temperature.
!
!-----------------------------------------------------------------------
!
!
!  AUTHOR: Keith Brewster
!  April, 1995   Based on OLAPS, hence LAPS, version of same.
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
  REAL :: tctotv,tt,ww
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  tctotv=(tt+273.15)*(1.+0.0006*ww)
  RETURN
  END FUNCTION tctotv
!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINE THWXN                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
!

SUBROUTINE thwxn(maxlev,nlevel,p,ht,t,td,elev,twdf)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!
!-----------------------------------------------------------------------
!
!
!  AUTHOR: Keith Brewster
!  April, 1995   Based on OLAPS, hence LAPS, version of same.
!
!  MODIFICATION HISTORY:
!  3/11/1996  Cleaned-up, removed OLAPS artifacts.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
!  Input variables
!
  INTEGER :: maxlev,nlevel
  REAL :: p(maxlev),ht(maxlev),t(maxlev),td(maxlev)
  REAL :: elev
!
!  Output variables
!
  REAL :: twdf
!
!  Functions
!
  REAL :: oe,tsa_fast
!
!  Misc internal variables
!
  INTEGER :: n
  REAL :: h3km,thaec,thw,twx,twn
!
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  twx=-999.
  twn=999.
  h3km=elev+3000.
  DO n=1,nlevel
    IF(ht(n) >= elev) THEN
      IF(ht(n) > h3km) EXIT
      thaec=oe(t(n),td(n),p(n))
      thw=tsa_fast(thaec,1000.)
      twx=AMAX1(twx,thw)
      twn=AMIN1(twn,thw)
    END IF
  END DO
!  101 CONTINUE
!
!  Find difference between max and min
!
  twdf=twx-twn
  RETURN
END SUBROUTINE thwxn
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE CALCLI                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE calcli(maxlev,nlevel,thepcl,p,t,li)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate Convective Available Potential Energy (CAPE)
!  in each column of ARPS grid.
!
!-----------------------------------------------------------------------
!
!
!  AUTHOR: Keith Brewster
!  April, 1995   Based on OLAPS, hence LAPS, version of same.
!
!  MODIFICATION HISTORY:
!  3/11/1996  Cleaned-up, removed OLAPS artifacts.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
!  Input variables
!
  INTEGER :: maxlev,nlevel
  REAL :: thepcl
  REAL :: p(maxlev),t(maxlev)
!
!  Output variable
!
  REAL :: li
!
!  Functions
!
  REAL :: tsa_fast
!
!  Misc internal variables
!
  INTEGER :: n
  REAL :: dp,wlow,t500,par500
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!
  DO n=2,nlevel-1
    IF( p(n) <= 500.) EXIT
  END DO
!  101 CONTINUE
  dp=ALOG(p(n-1)/p(n))
  wlow=ALOG(500./p(n))/dp
  t500=t(n)*(1.-wlow) + t(n-1)*wlow
  par500=tsa_fast(thepcl,500.)
!
  li=t500-par500
!
  RETURN
END SUBROUTINE calcli
!
!##################################################################
!##################################################################
!######                                                      ######
!######                 FUNCTION WMR2TD                      ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
!

  FUNCTION wmr2td(pres,wmr)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate Convective Available Potential Energy (CAPE)
!  in each column of ARPS grid.
!
!-----------------------------------------------------------------------
!
!
!  AUTHOR: Keith Brewster
!  April, 1995   Based on GEMPAK routine of same name.
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
  REAL :: wmr2td
  REAL :: pres
  REAL :: wmr
  REAL :: wkgkg,e,evap
!
  wkgkg = 0.001 * wmr
  wkgkg = AMAX1(wkgkg,0.00005)
  e= (pres*wkgkg) / (0.62197 + wkgkg)
  evap = e /(1.001 + (( pres - 100.) /900.) * 0.0034)
  wmr2td = ALOG(evap/6.112) * 243.5 /( 17.67 - ALOG (evap/6.112))

  RETURN
  END FUNCTION wmr2td
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE INTRPR                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
!

SUBROUTINE intrpr(nz,nlev,p,var,plvl,varatp)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Interpolate variable "var" linearly in log-pressure (p).
!
!-----------------------------------------------------------------------
!
!
!  AUTHOR: Keith Brewster
!  April, 1995   Based on OLAPS, hence LAPS, version of same.
!
!  MODIFICATION HISTORY:
!  3/11/1996  Cleaned-up, removed OLAPS artifacts.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
  INTEGER :: nz,nlev
  REAL :: p(nz),var(nz)
  REAL :: plvl
  REAL :: varatp
!
  INTEGER :: k
  REAL :: w1
!
  DO k=2,nlev-1
    IF(p(k) < plvl) EXIT
  END DO
!  101 CONTINUE
!
  w1=ALOG(p(k)/plvl)/ALOG(p(k)/p(k-1))
  varatp = w1*var(k-1) + (1.-w1)*var(k)
!
  RETURN
END SUBROUTINE intrpr

!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE CALCSHR                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE calcshr(nx,ny,nz,x,y,zp,sigma,                               &
           p_pa,t3d,u3d,v3d,cape,                                       &
           shr37,ustrm,vstrm,srlfl,srmfl,helicity,brn,brnu,blcon,       &
           tem1,tem2,tem3,elev,hgt3d)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate various wind shear parameters useful for gauging
!  the potential for severe storms.
!
!-----------------------------------------------------------------------
!
!
!  AUTHOR: Keith Brewster
!  February, 1994   Based on OLAPS, hence LAPS, version of same.
!
!  MODIFICATION HISTORY:
!  3/11/1996  Keith Brewster
!  Added storm-relative flows, general clean-up to
!  meet ARPS coding standards.
!
!  3/22/1996  (Keith Brewster)
!  Fixed some bugs, added smoothing at the end.
!
!  06/20/2000 (Eric Kemp and Keith Brewster)
!  Changed BRN Shear to be the denominator of BRN, instead of wind
!  speed, now has units of speed squared.
!
!  11/24/2003 (Keith Brewster)
!  Modified corner calculations for blcon.
!
!-----------------------------------------------------------------------
!
!  Calculates some of the shear related variables from the
!  ARPS 3D wind field.
!
!  shr37         Magnitude of wind shear between 3 and 7 km AGL
!  ustrm,vstrm   Estimated storm motion (modified from Bob Johns)
!  srlfl         Low-level storm-relative wind
!  srmfl         Mid-level storm-relative wind
!  helicity      Helicity, storm relative
!  brn           Bulk Richardson Number (Weisman and Klemp)
!  brnu          Shear parameter of BRN, "U"
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
!  Input variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: nx,ny,nz
  REAL :: x(nx)
  REAL :: y(ny)
  REAL :: zp(nx,ny,nz)
  REAL :: sigma(nx,ny)
  REAL :: p_pa(nx,ny,nz)
  REAL :: t3d(nx,ny,nz)
  REAL :: u3d(nx,ny,nz)
  REAL :: v3d(nx,ny,nz)
  REAL :: cape(nx,ny)
!
!-----------------------------------------------------------------------
!
!  Output variables
!
!-----------------------------------------------------------------------
!
  REAL :: shr37(nx,ny)      ! 7km - 3km wind shear
  REAL :: ustrm(nx,ny)
  REAL :: vstrm(nx,ny)  ! Estimated storm motion (Bob Johns)
  REAL :: srlfl(nx,ny)
  REAL :: srmfl(nx,ny)
  REAL :: helicity(nx,ny)   ! Helicity, storm relative
  REAL :: brn(nx,ny)        ! Bulk Richardson Number (Weisman and Klemp)
  REAL :: brnu(nx,ny)       ! Shear parameter of BRN, "U"
  REAL :: blcon(nx,ny)
!
!-----------------------------------------------------------------------
!
!  Temporary variables
!
!-----------------------------------------------------------------------
!
  REAL :: tem1(nx,ny)
  REAL :: tem2(nx,ny)
  REAL :: tem3(nx,ny)
  REAL :: elev(nx,ny)
  REAL :: hgt3d(nx,ny,nz)
!
!-----------------------------------------------------------------------
!
!  Misc internal variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k,ksfc,k2km,k3km
  REAL :: h3km,u3km,v3km,h7km,u7km,v7km,u2,v2
  REAL :: p2km,t2km,h2km,u2km,v2km,p9km,t9km,h9km,u9km,v9km
  REAL :: p500m,t500m,h500m,u500m,v500m,p6km,t6km,h6km,u6km,v6km
  REAL :: sumu,sumv,sump,sumh,wlow,whigh,dx,dy,dz,dp,dx2,dy2
  REAL :: rhohi,rholo,rhoinv
  REAL :: dirmean,spmean,ushr,vshr,ddir,perc
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  tem1=0.
  tem2=0.
  tem3=0.
!
  DO j=1,ny-1
    DO i=1,nx-1
      elev(i,j)=zp(i,j,2)
    END DO
  END DO

  DO k=1,nz-1
    DO j=1,ny-1
      DO i=1,nx-1
        hgt3d(i,j,k)=0.5*(zp(i,j,k)+zp(i,j,k+1))
      END DO
    END DO
  END DO
!
!-----------------------------------------------------------------------
!
!  Main loop over all horizontal scalar points
!
!-----------------------------------------------------------------------
!
  ksfc=2
  DO j=1,ny-1
    DO i=1,nx-1
!
!-----------------------------------------------------------------------
!
!  Find mass weighted mean wind in first 500m
!
!-----------------------------------------------------------------------
!
      sumu=0.
      sumv=0.
      sump=0.
      h500m=elev(i,j)+500.
      DO k=ksfc+1,nz-1
        IF( hgt3d(i,j,k) < h500m ) THEN
          dp=p_pa(i,j,k-1)-p_pa(i,j,k)
          rhohi=p_pa(i,j,k)/t3d(i,j,k)
          rholo=p_pa(i,j,k-1)/t3d(i,j,k-1)
          rhoinv=1./(rhohi+rholo)
          sumu=sumu+dp*rhoinv*(rholo*u3d(i,j,k-1)+rhohi*u3d(i,j,k))
          sumv=sumv+dp*rhoinv*(rholo*v3d(i,j,k-1)+rhohi*v3d(i,j,k))
          sump=sump+dp
        ELSE
          dz=hgt3d(i,j,k)-hgt3d(i,j,k-1)
          wlow=(hgt3d(i,j,k)-h500m)/dz
          u500m=u3d(i,j,k)*(1.-wlow) + u3d(i,j,k-1)*wlow
          v500m=v3d(i,j,k)*(1.-wlow) + v3d(i,j,k-1)*wlow
          p500m=p_pa(i,j,k)*(1.-wlow) + p_pa(i,j,k-1)*wlow
          t500m=t3d(i,j,k)*(1.-wlow) + t3d(i,j,k-1)*wlow
          dp=p_pa(i,j,k-1)-p500m
          rhohi=p500m/t500m
          rholo=p_pa(i,j,k-1)/t3d(i,j,k-1)
          rhoinv=1./(rhohi+rholo)
!        print *, ' p500,t500 = ',p500m,t500m
!        print *, ' u500,v500 = ',u500m,v500m
          sumu=sumu+dp*rhoinv*(rholo*u3d(i,j,k-1)+rhohi*u500m)
          sumv=sumv+dp*rhoinv*(rholo*v3d(i,j,k-1)+rhohi*v500m)
          sump=sump+dp
          EXIT
        END IF
      END DO
!      121    CONTINUE
      u500m=sumu/sump
      v500m=sumv/sump
!
!-----------------------------------------------------------------------
!
!  Find mass weighted mean wind sfc-2km AGL
!
!-----------------------------------------------------------------------
!
      sumu=0.
      sumv=0.
      sump=0.
      h2km=elev(i,j)+2000.
      DO k=ksfc+1,nz-1
        IF( hgt3d(i,j,k) < h2km ) THEN
          dp=p_pa(i,j,k-1)-p_pa(i,j,k)
          rhohi=p_pa(i,j,k)/t3d(i,j,k)
          rholo=p_pa(i,j,k-1)/t3d(i,j,k-1)
          rhoinv=1./(rhohi+rholo)
          sumu=sumu+dp*rhoinv*(rholo*u3d(i,j,k-1)+rhohi*u3d(i,j,k))
          sumv=sumv+dp*rhoinv*(rholo*v3d(i,j,k-1)+rhohi*v3d(i,j,k))
          sump=sump+dp
        ELSE
          dz=hgt3d(i,j,k)-hgt3d(i,j,k-1)
          wlow=(hgt3d(i,j,k)-h2km)/dz
          u2km=u3d(i,j,k)*(1.-wlow) + u3d(i,j,k-1)*wlow
          v2km=v3d(i,j,k)*(1.-wlow) + v3d(i,j,k-1)*wlow
          p2km=p_pa(i,j,k)*(1.-wlow) + p_pa(i,j,k-1)*wlow
          t2km=t3d(i,j,k)*(1.-wlow) + t3d(i,j,k-1)*wlow
          dp=p_pa(i,j,k-1)-p2km
          rhohi=p2km/t2km
          rholo=p_pa(i,j,k-1)/t3d(i,j,k-1)
          rhoinv=1./(rhohi+rholo)
          sumu=sumu+dp*rhoinv*(rholo*u3d(i,j,k-1)+rhohi*u2km)
          sumv=sumv+dp*rhoinv*(rholo*v3d(i,j,k-1)+rhohi*v2km)
          sump=sump+dp
          EXIT
        END IF
      END DO
!      141    CONTINUE
      u2km=sumu/sump
      v2km=sumv/sump
      tem2(i,j)=u2km
      tem3(i,j)=v2km
!
!-----------------------------------------------------------------------
!
!  Find mass weighted mean wind sfc-6km AGL
!
!-----------------------------------------------------------------------
!
      sumu=0.
      sumv=0.
      sump=0.
      h6km=elev(i,j)+6000.
      DO k=ksfc+1,nz-1
        IF( hgt3d(i,j,k) < h6km ) THEN
          dp=p_pa(i,j,k-1)-p_pa(i,j,k)
          rhohi=p_pa(i,j,k)/t3d(i,j,k)
          rholo=p_pa(i,j,k-1)/t3d(i,j,k-1)
          rhoinv=1./(rhohi+rholo)
          sumu=sumu+dp*rhoinv*(rholo*u3d(i,j,k-1)+rhohi*u3d(i,j,k))
          sumv=sumv+dp*rhoinv*(rholo*v3d(i,j,k-1)+rhohi*v3d(i,j,k))
          sump=sump+dp
        ELSE
          dz=hgt3d(i,j,k)-hgt3d(i,j,k-1)
          wlow=(hgt3d(i,j,k)-h6km)/dz
          u6km=u3d(i,j,k)*(1.-wlow) + u3d(i,j,k-1)*wlow
          v6km=v3d(i,j,k)*(1.-wlow) + v3d(i,j,k-1)*wlow
          p6km=p_pa(i,j,k)*(1.-wlow) + p_pa(i,j,k-1)*wlow
          t6km=t3d(i,j,k)*(1.-wlow) + t3d(i,j,k-1)*wlow
          dp=p_pa(i,j,k-1)-p6km
          rhohi=p6km/t6km
          rholo=p_pa(i,j,k-1)/t3d(i,j,k-1)
          rhoinv=1./(rhohi+rholo)
          sumu=sumu+dp*rhoinv*(rholo*u3d(i,j,k-1)+rhohi*u6km)
          sumv=sumv+dp*rhoinv*(rholo*v3d(i,j,k-1)+rhohi*v6km)
          sump=sump+dp
          EXIT
        END IF
      END DO
      u6km=sumu/sump
      v6km=sumv/sump
!
!-----------------------------------------------------------------------
!
!  Find mass weighted mean wind 2km-9km AGL
!
!-----------------------------------------------------------------------
!
      sumu=0.
      sumv=0.
      sump=0.
      h9km=elev(i,j)+9000.
      DO k=ksfc+1,nz-2
        IF( hgt3d(i,j,k) > h2km ) EXIT
      END DO
!      181   CONTINUE
      k2km=k
      dz=hgt3d(i,j,k)-hgt3d(i,j,k-1)
      wlow=(hgt3d(i,j,k)-h2km)/dz
      u2=u3d(i,j,k)*(1.-wlow) + u3d(i,j,k-1)*wlow
      v2=v3d(i,j,k)*(1.-wlow) + v3d(i,j,k-1)*wlow
      p2km=p_pa(i,j,k)*(1.-wlow) + p_pa(i,j,k-1)*wlow
      t2km=t3d(i,j,k)*(1.-wlow) + t3d(i,j,k-1)*wlow
      dp=p2km-p_pa(i,j,k)
      rholo=p2km/t2km
      rhohi=p_pa(i,j,k)/t3d(i,j,k)
      rhoinv=1./(rhohi+rholo)
      sumu=sumu+dp*rhoinv*(rholo*u2+rhohi*u3d(i,j,k))
      sumv=sumv+dp*rhoinv*(rholo*v2+rhohi*v3d(i,j,k))
      sump=sump+dp
      DO k=k2km+1,nz-1
        IF( hgt3d(i,j,k) < h9km ) THEN
          dp=p_pa(i,j,k-1)-p_pa(i,j,k)
          rhohi=p_pa(i,j,k)/t3d(i,j,k)
          rholo=p_pa(i,j,k-1)/t3d(i,j,k-1)
          rhoinv=1./(rhohi+rholo)
          sumu=sumu+dp*rhoinv*(rholo*u3d(i,j,k-1)+rhohi*u3d(i,j,k))
          sumv=sumv+dp*rhoinv*(rholo*v3d(i,j,k-1)+rhohi*v3d(i,j,k))
          sump=sump+dp
        ELSE
          dz=hgt3d(i,j,k)-hgt3d(i,j,k-1)
          wlow=(hgt3d(i,j,k)-h9km)/dz
          u9km=u3d(i,j,k)*(1.-wlow) + u3d(i,j,k-1)*wlow
          v9km=v3d(i,j,k)*(1.-wlow) + v3d(i,j,k-1)*wlow
          p9km=p_pa(i,j,k)*(1.-wlow) + p_pa(i,j,k-1)*wlow
          t9km=t3d(i,j,k)*(1.-wlow) + t3d(i,j,k-1)*wlow
          dp=p_pa(i,j,k-1)-p9km
          rhohi=p9km/t9km
          rholo=p_pa(i,j,k-1)/t3d(i,j,k-1)
          rhoinv=1./(rhohi+rholo)
          sumu=sumu+dp*rhoinv*(rholo*u3d(i,j,k-1)+rhohi*u9km)
          sumv=sumv+dp*rhoinv*(rholo*v3d(i,j,k-1)+rhohi*v9km)
          sump=sump+dp
          EXIT
        END IF
      END DO
!      191    CONTINUE
      u9km=sumu/sump
      v9km=sumv/sump
!
!-----------------------------------------------------------------------
!
!  Storm motion estimation
!  From Davies and Johns, 1993
!  "Some wind and instability parameters associated With
!  strong and violent tornadoes."
!  AGU Monograph 79, The Tornado...(Page 575)
!
!  Becuase of the discontinuity produced by that method
!  at the 15.5 m/s cutoff, their rules have been modified
!  to provide a gradual transition, and accomodate all the
!  data they mention in the article.
!
!-----------------------------------------------------------------------
!
      CALL uv2ddff(u6km,v6km,dirmean,spmean)
      IF(spmean >= 20.0) THEN
        dirmean=dirmean+18.
        IF(dirmean > 360.) dirmean=dirmean-360.
        spmean=spmean*0.89
      ELSE IF (spmean > 8.0) THEN
        whigh=(spmean - 8.0)/12.
        wlow =1.-whigh
        ddir=wlow*32.0 + whigh*18.0
        perc=wlow*0.75 + whigh*0.89
        dirmean=dirmean+ddir
        IF(dirmean > 360.) dirmean=dirmean-360.
        spmean=spmean*perc
      ELSE
        dirmean=dirmean+32.
        IF(dirmean > 360.) dirmean=dirmean-360.
        spmean=spmean*0.75
      END IF
      CALL ddff2uv(dirmean,spmean,ustrm(i,j),vstrm(i,j))
!
!-----------------------------------------------------------------------
!
!  Storm-relative low-level flow
!
!-----------------------------------------------------------------------
!
      srlfl(i,j)=SQRT((ustrm(i,j)-u2km)*(ustrm(i,j)-u2km) +             &
                      (vstrm(i,j)-v2km)*(vstrm(i,j)-v2km))
!
!-----------------------------------------------------------------------
!
!  Storm relative mid-level flow
!
!-----------------------------------------------------------------------
!
      srmfl(i,j)=SQRT((ustrm(i,j)-u9km)*(ustrm(i,j)-u9km) +             &
                      (vstrm(i,j)-v9km)*(vstrm(i,j)-v9km))
!
!-----------------------------------------------------------------------
!
!  Shear parameter for Bulk Richardson number
!
!-----------------------------------------------------------------------
!
!    print *, ' density-weight mean 0-500 m ',u500m,v500m
!    print *, ' density-weight mean 0-6  km ',u6km,v6km
!
      brnu(i,j)=0.5*( (u6km-u500m)*(u6km-u500m) +                       &
                      (v6km-v500m)*(v6km-v500m) )
!
!-----------------------------------------------------------------------
!
!  Bulk Richardson number
!  A limit of 200 is imposed, since this could
!  go to inifinity.
!
!-----------------------------------------------------------------------
!
      IF(brnu(i,j) > 0.) THEN
        brn(i,j)=cape(i,j)/brnu(i,j)
        brn(i,j)=AMIN1(brn(i,j),200.)
      ELSE
        brn(i,j)=200.
      END IF
!
    END DO
  END DO
!
!-----------------------------------------------------------------------
!
!  Calculate Helicity and 3km to 7km shear
!  since both involve the 3km wind.
!
!  For more efficient computation the Helicity is
!  computed for zero storm motion and the storm
!  motion is accounted for by adding a term at the end.
!  This is mathematically equivalent to accounting
!  for the storm motion at each level.
!
!-----------------------------------------------------------------------
!
  DO j=1,ny-1
    DO i=1,nx-1
      h3km=elev(i,j)+3000.
      h7km=elev(i,j)+7000.
!
!-----------------------------------------------------------------------
!
!  Find level just above 3 km AGL
!  Note, it is assumed here that there is at least
!  one level between the sfc and 3 km.
!
!-----------------------------------------------------------------------
!
      sumh=0.
      DO k=ksfc+1,nz-2
        IF(hgt3d(i,j,k) > h3km) EXIT
        sumh=sumh +                                                     &
            ( u3d(i,j,k)*v3d(i,j,k-1) ) -                               &
            ( v3d(i,j,k)*u3d(i,j,k-1) )
      END DO
!      240    CONTINUE
      k3km=k
      dz=hgt3d(i,j,k)-hgt3d(i,j,k-1)
      wlow=(hgt3d(i,j,k)-h3km)/dz
      u3km=u3d(i,j,k)*(1.-wlow) + u3d(i,j,k-1)*wlow
      v3km=v3d(i,j,k)*(1.-wlow) + v3d(i,j,k-1)*wlow
      sumh=sumh +                                                       &
            ( u3km*v3d(i,j,k-1) ) -                                     &
            ( v3km*u3d(i,j,k-1) )
      ushr=u3km-u3d(i,j,2)
      vshr=v3km-v3d(i,j,2)
      helicity(i,j)=sumh + vshr*ustrm(i,j) - ushr*vstrm(i,j)
!
!-----------------------------------------------------------------------
!
!  Now Find 7km wind for 3-to-7km shear
!
!-----------------------------------------------------------------------
!
      DO k=k3km,nz-1
        IF(hgt3d(i,j,k) > h7km) EXIT
      END DO
!      260   CONTINUE
      dz=hgt3d(i,j,k)-hgt3d(i,j,k-1)
      wlow=(hgt3d(i,j,k)-h7km)/dz
      u7km=u3d(i,j,k)*(1.-wlow) + u3d(i,j,k-1)*wlow
      v7km=v3d(i,j,k)*(1.-wlow) + v3d(i,j,k-1)*wlow
      shr37(i,j)=(SQRT( (u7km-u3km)*(u7km-u3km) +                       &
                        (v7km-v3km)*(v7km-v3km) ))/4000.
    END DO
  END DO

!-----------------------------------------------------------------------
!
!  Calculate the low-level convergence
!
!-----------------------------------------------------------------------
!
  dx=(x(2)-x(1))
  dy=(y(2)-y(1))
  dx2=2.*dx
  dy2=2.*dy
  DO j=2,ny-2
    DO i=2,nx-2
      blcon(i,j)=-1000.*(                                               &
                 (tem2(i+1,j)-tem2(i-1,j))/(sigma(i,j)*dx2)+            &
                 (tem3(i,j+1)-tem3(i,j-1))/(sigma(i,j)*dy2) )
    END DO
  END DO

  DO j=2,ny-2
    blcon(1,j)=-1000.*(                                                 &
                (tem2(2,j)-tem2(1,j))/(sigma(1,j)*dx)+                  &
                (tem3(1,j+1)-tem3(1,j-1))/(sigma(1,j)*dy2) )
    blcon(nx-1,j)=-1000.*(                                              &
           (tem2(nx-1,j)-tem2(nx-2,j))/(sigma(nx-1,j)*dx)+              &
           (tem3(nx-1,j+1)-tem3(nx-1,j-1))/(sigma(nx-1,j)*dy2) )
  END DO
  DO i=2,nx-2
    blcon(i,1)=-1000.*(                                                 &
                (tem2(i+1,1)-tem2(i-1,1))/(sigma(i,1)*dx2)+             &
                (tem3(i,2)-tem3(i,1))/(sigma(i,1)*dy) )
    blcon(i,ny-1)=-1000.*(                                              &
            (tem2(i+1,ny-1)-tem2(i-1,ny-1))/(sigma(i,ny-1)*dx2)+        &
            (tem3(i,ny-1)-tem3(i,ny-2))/(sigma(i,ny-1)*dy) )
  END DO
  blcon(1,1)=-1000.*(                                                   &
                (tem2(2,1)-tem2(1,1))/(sigma(1,1)*dx)+                  &
                (tem3(1,2)-tem3(1,1))/(sigma(1,1)*dy) )
  blcon(nx-1,1)=-1000.*(                                                &
                (tem2(nx-1,1)-tem2(nx-2,1))/(sigma(nx-1,1)*dx)+         &
                (tem3(nx-1,2)-tem3(nx-1,1))/(sigma(nx-1,1)*dy) )
  blcon(1,ny-1)=-1000.*(                                                &
           (tem2(2,ny-1)-tem2(1,ny-1))/(sigma(1,ny-1)*dx)+              &
           (tem3(1,ny-1)-tem3(1,ny-2))/(sigma(1,ny-1)*dy) )
  blcon(nx-1,ny-1)=-1000.*(                                             &
           (tem2(nx-1,ny-1)-tem2(nx-2,ny-1))/(sigma(nx-1,ny-1)*dx)+     &
           (tem3(nx-1,ny-1)-tem3(nx-1,ny-2))/(sigma(nx-1,ny-1)*dy) )
!
!-----------------------------------------------------------------------
!
!  Smooth the output arrays
!
!-----------------------------------------------------------------------
!
  CALL smooth9p(   shr37,nx,ny,1,nx-1,1,ny-1,0,tem1)
  CALL smooth9p(   ustrm,nx,ny,1,nx-1,1,ny-1,0,tem1)
  CALL smooth9p(   vstrm,nx,ny,1,nx-1,1,ny-1,0,tem1)
  CALL smooth9p(   srlfl,nx,ny,1,nx-1,1,ny-1,0,tem1)
  CALL smooth9p(   srmfl,nx,ny,1,nx-1,1,ny-1,0,tem1)
  CALL smooth9p(helicity,nx,ny,1,nx-1,1,ny-1,0,tem1)
  CALL smooth9p(     brn,nx,ny,1,nx-1,1,ny-1,0,tem1)
  CALL smooth9p(    brnu,nx,ny,1,nx-1,1,ny-1,0,tem1)
  CALL smooth9p(   blcon,nx,ny,1,nx-1,1,ny-1,0,tem1)
!
  RETURN
END SUBROUTINE calcshr
!
!##################################################################
!##################################################################
!######                                                      ######
!######                FUNCTION TSA_FAST                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

  FUNCTION tsa_fast(os,p)
!
!   THIS FUNCTION RETURNS THE TEMPERATURE TSA (CELSIUS) ON A SATURATION
!   ADIABAT AT PRESSURE P (MILLIBARS). OS IS THE EQUIVALENT POTENTIAL
!   TEMPERATURE OF THE PARCEL (CELSIUS). SIGN(A,B) REPLACES THE
!   ALGEBRAIC SIGN OF A WITH THAT OF B.
!
!    BAKER,SCHLATTER 17-MAY-1982     Original version
!    Modification for better convergence, Keith Brewster, Feb 1994.
!
!   B IS AN EMPIRICAL CONSTANT APPROXIMATELY EQUAL TO THE LATENT HEAT
!   OF VAPORIZATION FOR WATER DIVIDED BY THE SPECIFIC HEAT AT CONSTANT
!   PRESSURE FOR DRY AIR.
!
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate Convective Available Potential Energy (CAPE)
!  in each column of ARPS grid.
!
!-----------------------------------------------------------------------
!
!
!  AUTHOR: Keith Brewster
!  April, 1995   Based on OLAPS, hence LAPS, version of same.
!
!  MODIFICATION HISTORY:
!  3/11/1996  Cleaned-up, removed OLAPS artifacts.
!
!-----------------------------------------------------------------------
!
  REAL :: b
!  PARAMETER (B=2.6518986)
  PARAMETER (b=2651.8986)
  a= os+273.15
!
!   Above 200 mb figure all the moisture is wrung-out, so
!   the temperature is that which has potential temp of theta-e.
!   Otherwise iterate to find combo of moisture and temp corresponding
!   to thetae.
!
  IF( p < 200.) THEN
    tq=a*((p/1000.)**.286)
  ELSE
!   D IS AN INITIAL VALUE USED IN THE ITERATION BELOW.
    d= 120.
!   TQ IS THE FIRST GUESS FOR TSA.
    tq= 253.15
    x = 0.
!
!   ITERATE TO OBTAIN SUFFICIENT ACCURACY....SEE TABLE 1, P.8
!   OF STIPANUK (1973) FOR EQUATION USED IN ITERATION.
    DO i= 1,25
      d= 0.5*d

      x_last = x

      x= a*EXP(-b*f_mrsat(p*100.,tq)/tq)-tq*((1000./p)**.286)

      IF (ABS(x) < 1E-3) GO TO 2 
!
      IF (x_last * x < 0.) THEN
        slope = (x-x_last) / (tq - tq_last)
        delta = - x / slope
        ad = AMIN1(ABS(delta),d)
        tq_last = tq
        tq = tq + SIGN(ad,delta)
      ELSE
        tq_last = tq
        tq= tq+SIGN(d,x)
      END IF

    END DO
  END IF
2 tsa_fast = tq-273.15
  RETURN
  END FUNCTION tsa_fast
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE UV2DDFF                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE uv2ddff(u,v,dd,ff)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate direction and speed from u and v.
!
!-----------------------------------------------------------------------
!
!
!  AUTHOR: Keith Brewster
!  3/11/1996
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
  REAL :: u,v,dd,ff
  REAL :: dlon
  REAL :: r2deg
  PARAMETER (r2deg=180./3.141592654)
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!
  ff = SQRT(u*u + v*v)

  IF(v > 0.) THEN
    dlon=r2deg*ATAN(u/v)
  ELSE IF(v < 0.) THEN
    dlon=180. + r2deg*ATAN(u/v)
  ELSE IF(u >= 0.) THEN
    dlon=90.
  ELSE
    dlon=-90.
  END IF

  dd= dlon + 180.
  dd= dd-360.*(nint(dd)/360)
  RETURN
END SUBROUTINE uv2ddff
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE DDFF2UV                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE ddff2uv(dd,ff,u,v)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate u and v wind components from direction and speed.
!
!-----------------------------------------------------------------------
!
!
!  AUTHOR: Keith Brewster
!  3/11/1996
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
  REAL :: u,v,dd,ff
  REAL :: arg
  REAL :: d2rad
  PARAMETER (d2rad=3.141592654/180.)
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  arg = (dd * d2rad)
  u = -ff * SIN(arg)
  v = -ff * COS(arg)
  RETURN
END SUBROUTINE ddff2uv
