c
c
c		    THERMODYNAMIC SUBROUTINES AND FUNCTIONS
c
c				Richard Carpenter
c
c			    Univ. of Oklahoma, 1992
c
c
c
c    ANALYZSND	- Analyze a sounding
c    ANALYZSND2	- Analyze a sounding (does not recompute quantities)
c    GETCAPE	- Computes CAPE given T(v) of sounding, parcel
c    PARCEL	- Compute T,LWC of lifted parcel.
c    MODEL2RAW	- Convert (z,Th,Qv) sounding to (p,T,Td).
c    SATVAPOR	- es(T) using Bolton's (1980) best fit.
c    THQ2T	- Compute T given ThQ.
c    THE2T	- Compute T given ThE.
c
c
c
c
c
c                   ########################################
c                   ########################################
c                   ########################################
c                   ########                        ########
c                   ########       ANALYZSND2       ########
c                   ########                        ########
c                   ########################################
c                   ########################################
c                   ########################################
c
c
      Subroutine AnalyzSnd2 (n,pres,temp,tdew, irev, 
     > z,qv,qs,rh,tv,theta,
     > presLCL,tempLCL,zLCL,presLNB,zLNB,
     > alwcP,tempP,tvP,tvwlP, 
     &  CINtvwl,CAPEtvwl, CINtvnowl,CAPEtvnowl, CINt,CAPEt, ifCB)
c
!  12/08/93  Added 2 more CIN/CAPE parameters

      Implicit None
      Include 'thermo.consts'
c
      Integer	k,n, irev, maxp, ifCB, k1
      Parameter (maxp=1000)
      Real	pres(n), temp(n), tdew(n), z(n), alwcP(n), 
     >  qv(n), qs(n), rh(n), tv(n), theta(n), tempP(n),
     >  tvP(n), tvwlP(n), 
     >  dCAPE(maxp)
      Real  zLCL,presLCL, presLNB,zLNB, tempLCL, 
     >  e,es,qP,qvP,th,CINtvwl,CAPEtvwl, 
     >  CINt,CAPEt, presLNBt,zLNBt, CINtvnowl, CAPEtvnowl
      Include 'thermo.stfunc'
c
c
c  Input:
c     n		= number of points
c     pres	= pressure [Pa]
c     temp,tdew = temperature, dew point [K]
c     irev	= 0 for pseudoadabatic ascent, 1 for reversible
c     ifCB	= 0: surface-based parcel
c		  1: cloud base parcel. presLCL,tempLCL are inputs.
c     z		= height
c     qv 	= water vapor [kg/kg]
c     qs 	= saturation water vapor [kg/kg]
c     rh	= relative humidity
c     tv	= virtual temp [K]
c     theta	= potential temp
c
c  Note for ifCB=1 (cloud base parcel): The cloud base lies between k=1,2. 
c
c  Output:
c     alwcP	= parcel adiabatic liquid water content [kg/kg]
c     tempP	= parcel temperature [K]
c     tvP	= parcel virt temp [K]
c     tvwlP	= parcel virt temp, incl water loading [K]
c     CINtvwl	= negative CAPE [J/kg]
c     CAPEp	= positive CAPE [J/kg]
c     tempLCL	= parcel temperature at LCL
c     presLCL	= pressure at LCL
c
c
c
      If (n.GT.maxp) Then
	Print *, '% AnalyzSnd2: maxp too small. maxp,n: ', maxp,n
	Stop
      End If
c
c  Find LCL
c
      If (ifCB.EQ.0) Then
        th	= Ftheta(pres(1),temp(1))
        e	= Fsvpres(tdew(1)-tfrz)
        qP	= Fmixrat(pres(1),e)
        k1	= 1
        tempLCL	= Ftlcl(temp(1),tdew(1))
        presLCL	= Fplcl(pres(1),temp(1),tempLCL,qP)
        zLCL	= Fzlcl(temp(1),tempLCL,qP)
      Print '(1x,a,3f10.2)', '% AnalyzSnd2: LCL: pres,temp,z: ', 
     >  presLCL*1.e-2,tempLCL-tfrz,zLCL
      Else
        th	= Ftheta(presLCL,tempLCL)
        e	= Fsvpres(tempLCL-tfrz)
        qP	= Fmixrat(presLCL,e)
        k1	= 2
	zLCL	= -99.
      End If
c
c     Print *, '% AnalyzSnd2: ifCB:', ifCB
c     Print *, '% AnalyzSnd2: temp(1),tempLCL:', temp(1),tempLCL
c
c     print*,'% AnalyzSnd2:pres(1),temp(1),qP:',
c    >  pres(1)/100.,temp(1)-tfrz,qP*1000.
c
c
c  Compute Lifted Parcel quantities
c  For cloud base parcel, only want to recompute parcel properties: 
c  alwcP,tempP(2:n),tvP,tvwlP.
c
      Do k=k1,n
      es	= Fsvpres(temp(k)-tfrz)
      e		= Fsvpres(tdew(k)-tfrz)
c
c  ...parcel stuff
c
      If (pres(k).LT.presLCL) Then
        Call Parcel (presLCL,tempLCL,pres(k),irev, tempP(k),alwcP(k))
      Else
        tempP(k)= th * (pres(k)/p00)**rcp
	alwcP(k)= 0.
      End If
      qvP	= qP - alwcP(k)
      tvP(k)	= Ftvirtnowl(tempP(k),qvP)
      tvwlP(k)	= Ftvirtwl(tempP(k),qvP,alwcP(k))
c     If (ifCB.NE.0) Print *, k,pres(k)/100.,tempP(k)-tfrz,
c    >   qvP*1000.,alwcP(k)*1000.
c
      End Do
c
c
c  Supply values for k=1 for cloud base parcel
c
      If (ifCB.EQ.1) Then
      k		= 1
      qvP	= qP
      tempP(k)	= tempLCL
      tvP(k)	= Ftvirtnowl(tempP(k),qvP)
      tvwlP(k)	= Ftvirtwl(tempP(k),qvP,0.)
      End If
c
c
c
c
c  compute heights
c
      If (ifCB.NE.0) Then
      zLCL	= z(1) + (z(2)-z(1)) * 
     >         (Log(pres(1))-Log(presLCL)) / (Log(pres(1))-Log(pres(2)))
      Print *, '% AnalyzSnd2: zLCL: ', zLCL
      End If
c
c
c
c  Compute CAPE 
!	1. Using Tv(wl).  2. Using Tv.  3. Using T.
c
!	
c
      Call GetCAPE (n,maxp,ifCB,tv,tvP,z,pres,		!Tv
     >  dCAPE,CAPEtvnowl,CINtvnowl,presLNB,zLNB,zLCL)
c
      Call GetCAPE (n,maxp,ifCB,tv,tvwlP,z,pres,	!Tv(wl)
     >  dCAPE,CAPEtvwl,CINtvwl,presLNB,zLNB,zLCL)
c
      Call GetCAPE (n,maxp,ifCB,temp,tempP,z,pres,	!T
     >  dCAPE,CAPEt,CINt,presLNBt,zLNBt,zLCL)
c
 9901 Format (1x,a,3f8.0)
      Print 9901, '% AnalyzSnd2: CAPE (Tvwl,Tv,T):',
     &  CAPEtvwl, CAPEtvnowl, CAPEt
      Print 9901, '% AnalyzSnd2: CIN (Tvwl,Tv,T): ',
     &  CINtvwl, CINtvnowl, CINt
      Print 9901, '% AnalyzSnd2: pLNB:', presLNB*1.e-2, presLNBt*1.e-2
      Print 9901, '% AnalyzSnd2: zLNB:', zLNB, zLNBt
c
      End
c
c
c
c
c                   ########################################
c                   ########################################
c                   ########################################
c                   ########                        ########
c                   ########        ANALYZSND       ########
c                   ########                        ########
c                   ########################################
c                   ########################################
c                   ########################################
c
c
      Subroutine AnalyzSnd (n,pres,temp,tdew, irev, 
     > z,qv,qs,rh,tv,theta,
     > presLCL,tempLCL,zLCL,presLNB,zLNB,
     > alwcP,tempP,tvP,tvwlP,capem,capep, capemt,capept, ifCB)
c
      Implicit None
      Include 'thermo.consts'
c
      Integer	k,n, irev, maxp, ifCB, k1
      Parameter (maxp=1000)
      Real	pres(n), temp(n), tdew(n), z(n), alwcP(n), 
     >  qv(n), qs(n), rh(n), tv(n), theta(n), tempP(n),
     >  tvP(n), tvwlP(n), dcape(maxp)
      Real  dz,tvavg,zLCL,presLCL, presLNB,zLNB,tempLCL, 
     >  e,es,qP,qvP,th,capem,capep, dtvPavg, 
     >  capemt,capept, presLNBt,zLNBt
      Include 'thermo.stfunc'
c
c
c  Input:
c     n		= number of points
c     pres	= pressure [Pa]
c     temp,tdew = temperature, dew point [K]
c     irev	= 0 for pseudoadabatic ascent, 1 for reversible
c     ifCB	= 0: surface-based parcel
c		  1: cloud base parcel. presLCL,tempLCL are inputs.
c
c  Note for ifCB=1 (cloud base parcel): The cloud base lies between k=1,2. 
c
c  Output:
c     z		= height
c     qv 	= water vapor [kg/kg]
c     qs 	= saturation water vapor [kg/kg]
c     rh	= relative humidity
c     tv	= virtual temp [K]
c     theta	= potential temp
c     alwcP	= parcel adiabatic liquid water content [kg/kg]
c     tempP	= parcel temperature [K]
c     tvP	= parcel virt temp [K]
c     tvwlP	= parcel virt temp, incl water loading [K]
c     capem	= negative cape [J/kg]
c     capep	= positive cape [J/kg]
c     tempLCL	= parcel temperature at LCL
c     presLCL	= pressure at LCL
c
c
c
      If (n.GT.maxp) Then
	Print *, '% AnalyzSnd: maxp too small. maxp,n: ', maxp,n
	Stop
      End If
c
c  Find LCL
c
      If (ifCB.EQ.0) Then
        th	= Ftheta(pres(1),temp(1))
        e	= Fsvpres(tdew(1)-tfrz)
        qP	= Fmixrat(pres(1),e)
        k1	= 1
        tempLCL	= Ftlcl(temp(1),tdew(1))
        presLCL	= Fplcl(pres(1),temp(1),tempLCL,qP)
        zLCL	= Fzlcl(temp(1),tempLCL,qP)
      Print '(1x,a,3f10.2)', '% AnalyzSnd: LCL: pres,temp,z: ', 
     >  presLCL*1.e-2,tempLCL-tfrz,zLCL
      Else
        th	= Ftheta(presLCL,tempLCL)
        e	= Fsvpres(tempLCL-tfrz)
        qP	= Fmixrat(presLCL,e)
        k1	= 2
	zLCL	= -99.
      End If
c
c     Print *, '% AnalyzSnd: ifCB:', ifCB
c     Print *, '% AnalyzSnd: temp(1),tempLCL:', temp(1),tempLCL
c
c     print*,'% AnalyzSnd:pres(1),temp(1),qP:',
c    >  pres(1)/100.,temp(1)-tfrz,qP*1000.
c
c
c  For cloud base parcel, only want to recompute parcel properties: 
c  alwcP,tempP(2:n),tvP,tvwlP.
c
      Do k=k1,n
      es	= Fsvpres(temp(k)-tfrz)
      e		= Fsvpres(tdew(k)-tfrz)
      If (ifCB.EQ.0) Then
      theta(k)	= temp(k) * (p00/pres(k)) ** rcp
      qv(k)	= Fmixrat(pres(k),e)
      qs(k)	= Fmixrat(pres(k),es)
      rh(k)	= Frh_q(qv(k),qs(k))
      tv(k)	= Ftvirtnowl(temp(k),qv(k))
      End If
c
c  ...parcel stuff
c
      If (pres(k).LT.presLCL) Then
        Call Parcel (presLCL,tempLCL,pres(k),irev, tempP(k),alwcP(k))
      Else
        tempP(k)= th * (pres(k)/p00)**rcp
	alwcP(k)= 0.
      End If
      qvP	= qP - alwcP(k)
      tvP(k)	= Ftvirtnowl(tempP(k),qvP)
      tvwlP(k)	= Ftvirtwl(tempP(k),qvP,alwcP(k))
c     If (ifCB.NE.0) Print *, k,pres(k),tempP(k)-tfrz,
c    >   qvP*1000.,alwcP(k)*1000.
c
      End Do
c
c
c  Supply values for k=1 for cloud base parcel
c
      If (ifCB.EQ.1) Then
      k		= 1
      qvP	= qP
      tempP(k)	= tempLCL
      tvP(k)	= Ftvirtnowl(tempP(k),qvP)
      tvwlP(k)	= Ftvirtwl(tempP(k),qvP,0.)
      End If
c
c
c
c
c  compute heights
c
      If (ifCB.EQ.0) Then
      z(1)	= 0.
c
      Do k=2,n
	tvavg	= .5 * (tv(k-1) + tv(k))
	dz	= rd/grav*tvavg * Log(pres(k-1)/pres(k))
	z(k)	= z(k-1) + dz
      End Do
c
      Else
      zLCL	= z(1) + (z(2)-z(1)) * 
     >         (Log(pres(1))-Log(presLCL)) / (Log(pres(1))-Log(pres(2)))
      Print *, '% AnalyzSnd: zLCL: ', zLCL
      End If
c
c
c
c  Compute CAPE 
c
c
      Call GetCAPE (n,maxp,ifCB,tv,tvwlP,z,pres,
     >  dcape,capep,capem,presLNB,zLNB,zLCL)
c
      Call GetCAPE (n,maxp,ifCB,temp,tempP,z,pres,
     >  dcape,capept,capemt,presLNBt,zLNBt,zLCL)
c
 9901 Format (1x,a,2f8.0)
      Print 9901, '% AnalyzSnd: capep:', capep, capept
      Print 9901, '% AnalyzSnd: capem:', capem, capemt
      Print 9901, '% AnalyzSnd: pLNB:', presLNB*1.e-2, presLNBt*1.e-2
      Print 9901, '% AnalyzSnd: zLNB:', zLNB, zLNBt
c
      End
c
c
c
c
c                   ########################################
c                   ########################################
c                   ########################################
c                   ########                        ########
c                   ########        GETCAPE         ########
c                   ########                        ########
c                   ########################################
c                   ########################################
c                   ########################################
c
c
      Subroutine GetCAPE (n,maxp,ifCB,t,tP,z,pres,
     >  dcape,capep,capem,presLNB,zLNB,zLCL)
      Implicit None
      Include 'thermo.consts'
c
c
c  INPUT: (all SI units)
c     n		= number of points in sounding
c     maxp	= size of arrays
c     ifCB	= Compute relative to cloud base (given by zLCL)
c     t		= array of (virtual) temperatures from sounding
c     tP	= array of (virtual) temperatures of lifted parcel
c     z,pres	= array of heights,pressures
c     zLCL	= Hgt of cloud base (only if ifCB=1)
c
c  OUTPUT:
c     dcape	= array of CAPE for each level
c     zLNB,presLNB = hgt and pres of LNB
c     capep	= CAPE from LFC to LNB (positive area)
c     capem	= CAPE from sfc to LFC (negative area)
c
c
      Integer	n, maxp, k, ifCB
      Real	t(maxp), tP(maxp), z(maxp), pres(maxp), dcape(maxp)
      Real	t0avg,t1avg,dz,zLNB,capep,capem, presLNB, zLCL
      Logical	lnbflag
c
      capem	= 0.
      capep	= 0.
      lnbflag	= .False.
      presLNB	= pres(1)
      zLNB	= z(1)
c
c
c  Compute CAPE in each level
c  ...account for the special case of (k=1 && ifCB=1).
c
      Do k=1,n-1
      If (ifCB.EQ.1 .AND. k.EQ.1) Then
        t0avg	= t(k+1)
        t1avg	= .5 * (tP(k+1)-t(k+1))
        dz	= z(k+1) - zLCL
      Else
        t0avg	= .5 * (t(k) + t(k+1))
        t1avg	= .5 * (tP(k)-t(k) + tP(k+1)-t(k+1))
        dz	= z(k+1) - z(k)
      End If
      dcape(k)	= dz * grav * t1avg / t0avg
      End Do
c
c
c  Find LNB, total CAPE
c
      lnbflag	= .False.
      Do k=n-1,1,-1
c
c  ...find LNB
c
      If (.NOT.lnbflag .AND. dcape(k).GT.0.) Then
	presLNB	= pres(k+1)
	zLNB	= z(k+1)
	lnbflag	= .True.
      End If
c
c  ...accumulate CAPE below LNB 
c
      If (lnbflag) Then
      If (dcape(k).LT.0.) Then
	capem	= capem + dcape(k)
      Else
	capep	= capep + dcape(k)
      End If
      End If
c
      End Do
c
      End
c
c
c
c
c                   ########################################
c                   ########################################
c                   ########################################
c                   ########                        ########
c                   ########         PARCEL         ########
c                   ########                        ########
c                   ########################################
c                   ########################################
c                   ########################################
c
c
      Subroutine Parcel (presLCL, tempLCL, pres, irev, temp, alwc)
      Implicit None
      Include 'thermo.consts'
      Integer	irev
      Real	presLCL, tempLCL, pres, temp, alwc, thetae, qsLCL, qv,es
      Include 'thermo.stfunc'
c
c
c  Input:
c     presLCL	= pressure at LCL [mb]
c     tempLCL	= temperature at LCL [C]
c     pres	= pressure level to compute parcel at [mb]
c     irev	= 0 for pseudoadabatic ascent, 1 for reversible
c  Output:
c     temp	= temperature of parcel [C]
c     alwc	= adiabatic liquid water content of parcel [g/kg]
c
c
c  I don't know if this works below the LCL (that is, for unsat parcels)
c
c  We must find not only the temp but also the mix ratio
c
c
c  Reference thetae/q
c
c
      es	= Fsvpres(tempLCL-tfrz)
      qsLCL	= Fmixrat(presLCL,es)
      If (irev.EQ.0) Then
        thetae	= Fthetae(presLCL,tempLCL,tempLCL,qsLCL)
      Else
        thetae	= Fthetaq(presLCL,tempLCL,qsLCL,Fthq_cwcptrm(qsLCL))
      End If
c
c	Print *, 'Parcel: presLCL,tempLCL,qsLCL: ',presLCL,tempLCL,qsLCL
c	Print *, 'Parcel: thetae/q: ', thetae
c  Find T
c
      If (irev.EQ.0) Then
        Call ThE2T (temp, thetae, pres)
      Else
        Call ThQ2T (temp, thetae, pres, qsLCL, 0)
      End If
c
c  get ALWC
c
      es	= Fsvpres(temp-tfrz)
      qv	= Fmixrat(pres,es)
      alwc	= qsLCL - qv
c
c
      End
c
c
c
c
c                   ########################################
c                   ########################################
c                   ########################################
c                   ########                        ########
c                   ########        MODEL2RAW       ########
c                   ########                        ########
c                   ########################################
c                   ########################################
c                   ########################################
c
c
      Subroutine Model2Raw (presmb,tempc,tdewc, n,dz,zsfc,psfc)
      Implicit None
      Include 'thermo.consts'
      Integer	k,n
      Real	presmb(n), tempc(n), tdewc(n), 
     >  dz,zsfc,psfc,pres,theta,qv,tvavg,pkp1,
     >  tempkp1,thetakp1,qvkp1,tvkp1
      Include 'thermo.stfunc'
c
c
c  Input:
c     tempc	= theta
c     tdewc	= qv [kg/kg]
c     n		= number of points
c     dz	= grid spacing
c     zsfc	= station elevation (not used)
c     psfc	= surface pres
c
c  Output:
c     presmb	= pressure [mb]
c     tempc,tdewc = temperature, dew point [deg C]
c
c
      pkp1	= psfc
      Do k=1,n
      pres	= pkp1
      presmb(k)	= pres / 100.
      theta	= tempc(k)
      qv	= tdewc(k)
      tempc(k)	= theta * (pres / p00) ** rcp - tfrz
      tdewc(k)  = Ftdewc(Fvpres(pres,qv))
      tvavg	= (tempc(k)+tfrz) * (1. + cp608*qv)	!should be layer avg
      pkp1	= pres * Exp (-grav*dz*rdi/tvavg)
c
c  compute Tv at k+1 to get more accurate integration of pres.
c
      thetakp1	= tempc(k+1)
      qvkp1	= tdewc(k+1)
      tempkp1	= theta * (pkp1 / p00) ** rcp 
      tvkp1	= Ftvirtnowl(tempkp1,qvkp1)
      tvavg	= .5 * (tvavg + tvkp1)
      pkp1	= pres * Exp (-grav*dz*rdi/tvavg)
c
      End Do
c
      End
c
c
c
c
c                   ########################################
c                   ########################################
c                   ########################################
c                   ########                        ########
c                   ########        SATVAPOR        ########
c                   ########                        ########
c                   ########################################
c                   ########################################
c                   ########################################
c
c
c  Computes saturation vapor pressure
c  SI units
c  See Bolton (1980)
c
      Real Function SatVapor (temp)
      Implicit None
      Real	temp, g(0:7)
      Integer	i
      Data g /-2.9912729e3, -6.0170128e3,  1.887643854e1,-2.8354721e-2,
     >         1.7838301e-5,-8.4150417e-10,4.4412543e-13,  2.858487e0 /
c
      SatVapor = g(7) * Log (temp)
      Do i=0,6
      SatVapor = SatVapor + g(i) * temp ** (i-2)
      End Do
      SatVapor = Exp (SatVapor)
      End
c
c
c
c
c                   ########################################
c                   ########################################
c                   ########################################
c                   ########                        ########
c                   ########          THQ2T         ########
c                   ########                        ########
c                   ########################################
c                   ########################################
c                   ########################################
c
c
      Subroutine ThQ2T (temp, thetaq, pres, q, iqflag)
      Implicit None
      Include 'thermo.consts'
c
      Integer	iter, iqflag
      Real	temp, thetaq, pres, tempc, es, qs, thetaq1, err, 
     >  delta, tempold, q, qv, qtotal,
     >  deltaold, stepfrac
      Include 'thermo.stfunc'
c
c  NOTE: THIS ASSUMES THE PARCEL IS JUST SATURATED.
c
c  Given ThetaQ and Pres, find Temp
c  Note that, for fixed pres, ThetaQ varies with both temp and q.
c
c  Input: 
c     temp	= temperature (K) (first guess)
c     thetaq	= ThetaQ
c     pres	= pressure (Pa)
c     q		= Mixing ratio (kg/kg) (iqflag=0 only)
c  Output:
c     temp	= temperature (K)
c     qtotal	= Mixing ratio (kg/kg)
c     iqflag	= 0: Hold q const.  1: Set q to satd value.
c
c
      If (temp.LT.200.) temp = 300.
      stepfrac	= 0.5
c     print *, 'ThQ2T: temp,thetaq,pres: ', temp,thetaq,pres
c
      Do iter=1,40
      tempc	= temp - tfrz
      es	= Fsvpres(tempc)
      qs	= Fmixrat(pres,es)
      If (iqflag.EQ.0) Then
	qv	= Min(qs,q)
	qtotal	= q
      Else
	qv	= qs
	qtotal	= qs
      End If
      thetaq1	= Fthetaq(pres,temp,qv,Fthq_cwcptrm(qtotal))
c     print *, 'ThQ2T: tempc,es,qs,thetaq1: ', tempc,es,qs,thetaq1
      err	= thetaq1 - thetaq
      If (Abs(err).LT.0.001) Go To 9000
c
c  limit the step to 20 K or 1/2 the error in ThE
c
      deltaold	= delta
      delta	= Min (40., stepfrac*Abs(err))
      delta	= Sign (delta, err)
      If (iter.GT.2 .AND. Sign(1.,delta).NE.Sign(1.,deltaold)) Then
	stepfrac = stepfrac * .5
c	Print *, 'ThQ2T: reducing stepfrac to :', stepfrac
        delta	= Min (40., stepfrac*Abs(err))
        delta	= Sign (delta, err)
      End If
      tempold	= temp
      temp	= temp - delta
c     Print 9999,'ThQ2T: iter,temp,err,Err:',iter,temp,err,temp-tempold
      End Do
c
      Print *, 'ThQ2T: did not converge. iter,temp: ', iter,temp
c
 9000 Continue
c     Print *, 'ThQ2T: temp,err: ', temp,temp-tempold
      End
c
c
c
c
c                   ########################################
c                   ########################################
c                   ########################################
c                   ########                        ########
c                   ########          THE2T         ########
c                   ########                        ########
c                   ########################################
c                   ########################################
c                   ########################################
c
c
      Subroutine ThE2T (temp, thetae, pres)
      Implicit None
      Include 'thermo.consts'
c
      Integer	iter
      Real	temp, thetae, pres, tempc, es, qs, thetae1, err, 
     >  delta, deltaold, stepfrac, tempold, tlcl
      Include 'thermo.stfunc'
c
c  NOTE: THIS ASSUMES THE PARCEL IS JUST SATURATED.
c
c  Given ThetaE and Pres, find Temp
c  Note that, for fixed pres, ThetaE varies with both temp and q.
c
c  Input: 
c     temp	= temperature (K) (first guess)
c     thetae	= ThetaE
c     pres	= pressure (Pa)
c  Output:
c     temp	= temperature (K)
c
c
      If (temp.LT.200.) temp = 300.
      stepfrac	= 0.5
c     print *, 'ThE2T: temp,thetae,pres: ', temp,thetae,pres
c
      Do iter=1,30
      tempc	= temp - tfrz
      es	= Fsvpres(tempc)
      qs	= Fmixrat(pres,es)
      tlcl	= Ftlcl2(temp,es)
      thetae1	= Fthetae(pres,temp,tlcl,qs)
c     print *, 'ThE2T: tempc,es,qs,thetae1: ', tempc,es,qs,thetae1
      err	= thetae1 - thetae
      If (Abs(err).LT.0.001) Go To 9000
c
c  limit the step to 20 K or 1/2 the error in ThE
c
      deltaold	= delta
      delta	= Min (40., stepfrac*Abs(err))
      delta	= Sign (delta, err)
      If (iter.GT.2 .AND. Sign(1.,delta).NE.Sign(1.,deltaold)) Then
	stepfrac = stepfrac * .5
c	Print *, 'ThE2T: reducing stepfrac to :', stepfrac
        delta	= Min (40., stepfrac*Abs(err))
        delta	= Sign (delta, err)
      End If
      tempold	= temp
      temp	= temp - delta
c     Print 9999,'ThE2T: iter,temp,err,Err:',iter,temp,err,temp-tempold
      End Do
c
      Print 9999, 'ThE2T: did not converge. iter,temp,err: ', 
     >  iter,temp,err,temp-tempold
c
 9999 Format (1x,a,i5,8f10.3)
 9000 Continue
c     Print *, 'ThE2T: temp,err: ', temp,temp-tempold
      End
