c
c
c                   ########################################
c                   ########################################
c                   ########################################
c                   ########                        ########
c                   ########        Interpress      ########
c                   ########       JJM smr '95      ########
c                   ########                        ########
c                   ########################################
c                   ########################################
c                   ########################################
c
      Subroutine Interpress (n,pres_S,tempc_S,tdewc_S,spd_S,
     > zz_S,dir_S,pres850,ht850,temp850,tdewp850,spd850,
     > dir850,pres700,ht700,temp700,tdewp700,spd700,dir700,
     > pres500,ht500,temp500,tdewp500,spd500,dir500,nmax)
c
      Integer n,nmax 
      Real pres_S(nmax),tempc_S(nmax),tdewc_S(nmax)
      Real spd_S(nmax),dir_S(nmax),zz_S(nmax)
      Real pres850,ht850,temp850,tdewp850,spd850,dir850
      Real pres700,ht700,temp700,tdewp700,spd700,dir700
      Real pres500,ht500,temp500,tdewp500,spd500,dir500
      Real inc
c
      pres850=0.
      pres700=0.
      pres500=0.
c
c      Do i=1,n,1    ! TEST TO FIX SOUNDINGS STARTING ABOVE 850.0 mb 
      Do i=2,n,1
        If (pres_S(i) .LT. 85000.0) Then
          If (pres_S(i-1) .GT. 85000.0) Then
            inc=ALOG(pres_S(i-1))-ALOG(pres_S(i))
            inc=(ALOG(pres_S(i-1))-ALOG(85000.0))/inc
c
            pres850=ALOG(pres_S(i-1))-inc*(ALOG(pres_S(i-1))-
     >       ALOG(pres_S(i)))
            pres850=EXP(pres850)/100.
c
            ht850=zz_S(i-1)-inc*(zz_S(i-1)-
     >       zz_S(i))
c
            temp850=tempc_S(i-1)-inc*
     >       (tempc_S(i-1)-tempc_S(i))
c
            tdewp850=tdewc_S(i-1)-inc*
     >       (tdewc_S(i-1)-tdewc_S(i))
c
            if ((spd_S(i-1).eq.(0.)).or.(spd_S(i).eq.(0.)))
     +      then ! EMK 1/31/98
              spd850 = 0.
            else
              spd850=ALOG(spd_S(i-1))-inc*(ALOG(spd_S(i-1))-
     >        ALOG(spd_S(i)))
              spd850=EXP(spd850)
            end if ! EMK
c
            dir850=dir_S(i)
          End If
        Else If (pres_S(i) .EQ. 85000.0) Then
          pres850=pres_S(i)/100.
          ht850=zz_S(i)	
          temp850=tempc_S(i)
          tdewp850=tdewc_S(i)
          spd850=spd_S(i)
          dir850=dir_S(i)
        End If
c
c
        If (pres_S(i) .LT. 70000.0) Then
          If (pres_S(i-1) .GT. 70000.0) Then
            inc=ALOG(pres_S(i-1))-ALOG(pres_S(i))
            inc=(ALOG(pres_S(i-1))-ALOG(70000.0))/inc
c
            pres700=ALOG(pres_S(i-1))-inc*(ALOG(pres_S(i-1))-
     >       ALOG(pres_S(i)))
            pres700=EXP(pres700)/100.
c
            ht700=zz_S(i-1)-inc*(zz_S(i-1)-
     >       zz_S(i))
c
            temp700=tempc_S(i-1)-inc*
     >       (tempc_S(i-1)-tempc_S(i))
c
            tdewp700=tdewc_S(i-1)-inc*
     >       (tdewc_S(i-1)-tdewc_S(i))
c

            if ((spd_s(i-1).eq.(0.)).or.(spd_s(i).eq.(0.)))
     +      then ! EMK 1/31/98
              spd700 = 0.
            else
              spd700=ALOG(spd_S(i-1))-inc*(ALOG(spd_S(i-1))-
     >        ALOG(spd_S(i)))
              spd700=EXP(spd700)
            end if ! EMK
c
            dir700=dir_S(i)
          End If
        Else If (pres_S(i) .EQ. 70000.0) Then
          pres700=pres_S(i)/100.
          ht700=zz_S(i)	
          temp700=tempc_S(i)
          tdewp700=tdewc_S(i)
          spd700=spd_S(i)
          dir700=dir_S(i)
        End If
c
c
        If (pres_S(i) .LT. 50000.0) Then
          If (pres_S(i-1) .GT. 50000.0) Then
            inc=ALOG(pres_S(i-1))-ALOG(pres_S(i))
            inc=(ALOG(pres_S(i-1))-ALOG(50000.0))/inc
c
            pres500=ALOG(pres_S(i-1))-inc*(ALOG(pres_S(i-1))-
     >       ALOG(pres_S(i)))
            pres500=EXP(pres500)/100.
c
            ht500=zz_S(i-1)-inc*(zz_S(i-1)-
     >       zz_S(i))
c
            temp500=tempc_S(i-1)-inc*
     >       (tempc_S(i-1)-tempc_S(i))
c
            tdewp500=tdewc_S(i-1)-inc*
     >       (tdewc_S(i-1)-tdewc_S(i))
c
            if ((spd_S(i-1).eq.(0.)).or.(spd_S(i).eq.(0.)))
     +      then ! EMK 1/31/98
              spd500 = 0.
            else
              spd500=ALOG(spd_S(i-1))-inc*(ALOG(spd_S(i-1))-
     >        ALOG(spd_S(i)))
              spd500=EXP(spd500)
            end if ! EMK
c
            dir500=dir_S(i)
          End If
        Else If (pres_S(i) .EQ. 50000.0) Then
          pres500=pres_S(i)/100.
          ht500=zz_S(i)
          temp500=tempc_S(i)
          tdewp500=tdewc_S(i)
          spd500=spd_S(i)
          dir500=dir_S(i)
        End If
      End Do
c
      Return
c
      End
c
c
c                   ########################################
c                   ########################################
c                   ########################################
c                   ########                        ########
c                   ########         INDICES        ########
c                   ########       JJM smr '95      ########
c                   ########                        ########
c                   ########################################
c                   ########################################
c                   ########################################
c
      Subroutine Indices (pres850,ht850,temp850,tdewp850,
     > spd850,dir850,pres700,ht700,temp700,tdewp700,spd700,
     > dir700,pres500,ht500,temp500,tdewp500,spd500,dir500,
     > showalter,kindex,liftedindex,totaltotals,sweat,
     > presLCL,zLCL,tempLCL,mixrat,prestopBL,pres_S,tempc_S,
     > zz_S,tdewc_S,n,nmax,convtemp,precipwat,htwbz,pc,zLFC,
     > presLFC,zLNB,uu_S,vv_S,BRNshear,lidstrength,dirmean,
     > spmean,helicity)
c
c
      Real pres850,ht850,temp850,tdewp850,spd850,dir850
      Real pres700,ht700,temp700,tdewp700,spd700,dir700
      Real pres500,ht500,temp500,tdewp500,spd500,dir500
      Real showalter,kindex,liftedindex,totaltotals
      Real sweat,presLCL,zLCL,tempLCL,shearterm
      Real wobf,thw,tcon,powt,satlft,p,drydist,showt
      Real t850parc_LCL,mixrat,prestopBL,convtemp,ct
      Real precipwat,precpw,pccl,inc,htwbz,pc
      Real wettemp,wettemplast,zz_S(nmax),u6km,v6km 
      Real pres_S(nmax),tempc_S(nmax),tdewc_S(nmax)
      Real zLFC,presLFC,equivpotTLCL,zLNB,tempclast,preslast
      Real mstTatP,mstTatPlast,tempc,pres,incT,incP
      Real brnsh,BRNshear,uu_S(nmax),vv_S(nmax)
      Real dewpterm,tt_term,lidstrength,this,max,diff
      Real theta_swl,avg_theta_w,thispowt,totpowt
      Real lidstrengthA,lidstrengthB,avg_theta_sw
      Real theta_sw_tot,theta_sw_max,theta_sw
      Real dirmean,spmean,ustrm,vstrm,whigh,wlow,perc,ddir
      Real sumh,h3km,dz,ushr,vshr,helicity,dirtemp,sptemp
      Real diff500,diff700,diff850
      Integer n,nmax,i,j,k,l,k3km
      Logical LFC,LID,x8,x7,x5
c
c
c     Check to see that each of the levels used in the
c     index computations are valid (i.e. the sounding
c     doesn't truncate below the given level)
c
c
      x8=.true.
      x7=.true.
      x5=.true.
c
      diff850=abs(pres850-850.)
      diff700=abs(pres700-700.)
      diff500=abs(pres500-500.)
c
      If (diff850 .gt. 5.) Then
        x8=.false.
        x7=.false.
        x5=.false.
      End If
c
      If (diff700 .gt. 5.) Then
        x7=.false.
        x5=.false.
      End If
c
      If (diff500 .gt. 5.) Then
        x5=.false.
      End If
c
c
c     Compute the K-Index..
c
c
      kindex=0.
      If (x8 .and. x7) Then
        kindex=(temp850-temp500)+tdewp850-(temp700-tdewp700)
      End If
c
c
c     Compute the Totaltotals Index..
c
c
      totaltotals=0.
      If (x8 .and. x5) Then
        totaltotals=temp850+tdewp850-2*temp500 
      End If
c
c
c	Compute SWEAT index..
c
c
      dewpterm=12*tdewp850
      tt_term=20*(totaltotals-49)
      shearterm=125*(sin((dir500-dir850)/360*2*3.14)+.2)
c
      If (tt_term .LT. 0) Then
        tt_term=0.
      End If
c
      If (dewpterm .LT. 0) Then
        dewpterm=0.
      End If
c
      If (dir850 .LT. 130.0) Then
        shearterm=0
      Else If (dir850 .GT. 250.0) Then
        shearterm=0
      Else If (dir500 .LT. 210.0) Then
        shearterm=0
      Else If (dir500 .GT. 310.0) Then
        shearterm=0
      Else If ((dir500-dir850) .LT. 0) Then
        shearterm=0
      Else If ((spd850/.5144) .LE. 15.0) Then
        shearterm=0
      Else If ((spd500/.5144) .LE. 15.0) Then
        shearterm=0
      End If
c
      sweat=0.
      If (x8 .and. x5) Then
        sweat=(dewpterm)+(tt_term)+(2*(spd850/.5144))+ ! EMK 1/31/98
     >   (spd500/.5144)+(shearterm)                    ! EMK 1/31/98

*        sweat=dewpterm+tt_term+2*(spd850/.5144)+ ! Original
*     >   spd500/.5144+shearterm                  ! Original
      End If

! RLC 1998/05/20
! Test for bad pressures (call to POWT)

      IF (presLCL .LE. 0.0) THEN
	PRINT *, 'INDICES: WARNING: Bad presLCL: ', presLCL, 
     &    ', resetting to 1000 mb'
	presLCL = 1000.0E2
      END IF

      IF (pres850 .LE. 0.0) THEN
	PRINT *, 'INDICES: WARNING: Bad pres850: ', pres850, 
     &    ', resetting to 850 mb'
	pres850 = 850.0E2
      END IF

c
c
c     Compute the Lifted Index (surface based) ..
c
c
      thw=powt((tempLCL-273.16),(presLCL/100),(tempLCL-273.16))
      p=500.0
c
      liftedindex=0.
      If (x5) Then
        liftedindex=temp500-satlft(thw,p)
      End If
c
c
c     Compute the Showalter Index..
c
c
      If (presLCL/100. .LT. pres850) Then
        drydist=(zLCL-ht850)/1000
        t850parc_LCL=temp850-drydist*9.77
        showt=satlft(powt(t850parc_LCL,presLCL/100.,
     >   t850parc_LCL),p)
      Else
        showt=satlft(powt(temp850,pres850,temp850),p)
      End If
c
      showalter=0.
      If (x8 .and. x5) Then
        showalter=temp500-showt
      End If
c
c
c     Compute the Convective Temperatue..
c
c
      do i=1,n
        pres_S(i)=pres_S(i)/100.
      end do
c
      pc=pccl(prestopBL/100.,pres_S,tempc_S,tdewc_S,mixrat*1000,
     > nmax,n)
      convtemp=ct(mixrat*1000.,pc,pres_S(1))
c
      do i=1,n
        pres_S(i)=pres_S(i)*100.
      end do
c
c
c     Compute Precipitable Water
c
c
      precipwat=precpw(tdewc_S,pres_S,n)
c
c
c     Compute the Wet Bulb Zero Height..
c
c

      htwbz = 0.
      wettemplast=tw(tempc_S(1),tdewc_S(1),(pres_S(1)/100.))
      Do i=1,n,1
        wettemp=tw(tempc_S(i),tdewc_S(i),(pres_S(i)/100.))
        If (wettemp .LT. 0.) Then
          If (wettemplast .GT. 0.) Then
            inc=(wettemplast)/(wettemplast-wettemp)
            htwbz=zz_S(i-1)+inc*(zz_S(i)-zz_S(i-1)) 
          End If
        End If
        wettemplast=wettemp
      End Do
c
c
c     Compute the LFC..  (seems to work, I couldn't find any
c     program to model this after, so its my try.. JJM)
c
c
      LFC = .FALSE.
      zLFC=0.
      presLFC=0.
      equivpotTLCL=os((tempLCL-273.16),(presLCL/100.))
      Do i=2,n,1
        If (pres_S(i) .LE. presLCL) Then
          mstTatP=tmlaps(equivpotTLCL,(pres_S(i)/100.))
          mstTatPlast=tmlaps(equivpotTLCL,(pres_S(i-1)/100.))
          If (tempc_S(i) .LE. mstTatP) Then
            If (tempc_S(i-1) .GT. mstTatPlast) Then
              incT=(tempc_S(i-1)-tempc_S(i))/100.
              incP=(pres_S(i-1)-pres_S(i))/100.
              tempc=tempc_S(i-1)
              pres=pres_S(i-1)
              Do j=1,100,1
                tempclast=tempc
                preslast=pres	
                tempc=tempc-incT
                pres=pres-incP
                mstTatP=tmlaps(equivpotTLCL,(pres/100.))
                If (tempc .LE. mstTatP) Then
	          If (tempclast .GT. mstTatPlast) Then
                    presLfC=(pres+preslast)/2
                    LFC = .TRUE.	
                  End If
                End If
                mstTatPlast=mstTatP
              End Do
              inc=(ALOG(pres_S(i-1))-ALOG(presLFC))/(ALOG(pres_S(i-1))
     >         -ALOG(pres_S(i)))
              zLFC=zz_S(i-1)+inc*(zz_S(i)-zz_S(i-1))
            End If
          End If 
        End If
      End Do
c
c
c     Compute the Lid Strength Index.. (This is also written
c     purely by me.  There seems to be no nice way to numer-
c     ically find the lid, so my method is arbitrary and seems
c     to work pretty fairly.  JJM smr '95)  (see Carlson, Yr:?
c     for more info on LSI)
c
c
      LID = .false. ! EMK 1/31/98
      lidstrengthA=0.
      lidstrengthB=0.
c
      If (LFC) Then
        max=-100.
        j=0
        this=powt((tempLCL-273.16),(presLCL/100.),
     >   (tempLCL-273.16))
        Do i=2,n,1
          dz=(zz_S(i)-zz_S(i-1))/1000.
          dT=tempc_S(i)-tempc_S(i-1)
          If (dT/dz .GE. -6.0) Then
            diff=tempc_S(i)-satlft(this,(pres_S(i)/100.))
            If (zz_S(i) .LE. (ht700+ht500)/2) Then
              max=amax1(diff,max)
              If (max .EQ. diff) Then
                j=i
              End If
            End If
          End If
        End Do
        If (j .NE. 0) Then
          Print *,'Lid found at:',zz_S(j),'m, ',pres_S(j)/100.,'mb'
          LID = .TRUE.
        Else
          Print *,'No lid detected!??'
        End If
      Else
        Print *,'No LFC, therefore no lid computed!'
      End If
c
      k=0
      thispowt=0.
      totpowt=0.
      avg_theta_w=0.
      theta_swl=0.
      dp=0.
      sumdp=0.
      Do i=1,n,1
        If (pres_S(i) .GE. (pres_S(1)- 5000.)) Then
          dp=pres_S(i)-pres_S(i+1)
          sumdp=sumdp+dp	
          IF (pres_S(i) .LE. 0.0) THEN
	    PRINT *, 'INDICES: WARNING: Bad i, pres_S(i): ', i,pres_S(i)
	    thispowt = 273.15
	  ELSE
          thispowt=powt(tempc_S(i),(pres_S(i)/100.),tdewc_S(i))
          END IF
          totpowt=totpowt+thispowt*dp
        End If
      End Do
      avg_theta_w=totpowt/sumdp
      If (LID) Then 
        IF (pres_S(j) .LE. 0.0) THEN
	  PRINT *, 'INDICES: WARNING: Bad j, pres_S(i): ', i,pres_S(j)
	  theta_swl = 273.15
	ELSE
          theta_swl=powt(tempc_S(j),(pres_S(j)/100.),tempc_S(j))
        END IF
        lidstrengthB=theta_swl-avg_theta_w
      Else
        lidstrengthB=0.
      End If
c
      k=0
      theta_sw_max=-100
      Do i=1,n,1
        If (pres_S(i) .GT. 50000.) Then 
          theta_sw=powt(tempc_S(i),(pres_S(i)/100.),tempc_S(i))
          theta_sw_max=amax1(theta_sw_max,theta_sw)
          If (theta_sw_max .EQ. theta_sw) Then
            k=i
          End If
        End If
      End Do
c
      l=0
      theta_sw_tot=0.
      sumdp=0.
      dp=0.
      Do i=k,n,1
        If (pres_S(i) .GE. 50000.) Then
          dp=pres_S(i)-pres_S(i+1)
          sumdp=sumdp+dp
          theta_sw=powt(tempc_S(i),(pres_S(i)/100.),tempc_S(i))
          theta_sw_tot=theta_sw_tot+theta_sw*dp
        End If
      End Do
      avg_theta_sw=theta_sw_tot/sumdp
c
      lidstrengthA=avg_theta_w-avg_theta_sw
c
      lidstrength=lidstrengthA-lidstrengthB
c
c
c     Compute the Bulk Richardson # Shear..
c
c
      BRNshear=brnsh(pres_S,zz_S,tempc_S,uu_S,vv_S,n,nmax,u6km,v6km) 
c
c
c     Compute the Storm Motion Vector and the SREH..
c
c     Storm motion estimation
c     From Davies and Johns, 1993
c     "Some wind and instability parameters associated With
c     strong and violent tornadoes."
c     AGU Monograph 79, The Tornado...(Page 575)
c
c     (modified by JJM)
c
c     Becuase of the discontinuity produced by that method
c     at the 15.5 m/s cutoff, their rules have been modified
c     to provide a gradual transition, and accomodate all the
c     data they mention in the article.
c
      call get_ddff(u6km,v6km,dirmean,spmean)
c
      IF(spmean .ge. 20.0) THEN
        dirmean=dirmean+18.
        IF(dirmean.gt.360.) dirmean=dirmean-360.
        spmean=spmean*0.89
      ELSE IF (spmean .gt. 8.0) THEN
        whigh=(spmean - 8.0)/12.
        wlow =1.-whigh
        ddir=wlow*32.0 + whigh*18.0
        perc=wlow*0.75 + whigh*0.89
        dirmean=dirmean+ddir
        IF(dirmean.gt.360.) dirmean=dirmean-360.
        spmean=spmean*perc
      ELSE
        dirmean=dirmean+32.
        IF(dirmean.gt.360.) dirmean=dirmean-360.
        spmean=spmean*0.75
      END IF
c
      dirtemp=dirmean
      sptemp=spmean
      call get_uv(dirmean,spmean,ustrm,vstrm)
      dirmean=dirtemp
      spmean=sptemp
c
c
c     Calculate Helicity 
c
c     For more efficient computation the Helicity is
c     computed for zero storm motion and the storm
c     motion is accounted for by adding a term at the end.
c     This is mathematically equivalent to accounting
c     for the storm motion at each level.
c
        h3km=3000.
c
c     Find level just above 3 km AGL
c     Note, it is assumed here that there is at least
c     one level between the sfc and 3 km.
c
      sumh=0.
      DO k=2,n
        IF(zz_S(k).gt.h3km) GO TO 240
        sumh=sumh +
     :    ( uu_S(k)*vv_S(k-1) ) -
     :    ( vv_S(k)*uu_S(k-1) )
      END DO
 240  CONTINUE
      k3km=min(k,n)
      dz=zz_S(k3km)-zz_S(k3km-1)
      wlow=(zz_S(k3km)-h3km)/dz
      u3km=uu_S(k3km)*(1.-wlow) + uu_S(k3km-1)*wlow
      v3km=vv_S(k3km)*(1.-wlow) + vv_S(k3km-1)*wlow
      sumh=sumh +
     :  (u3km*vv_S(k3km-1)) -
     :  (v3km*uu_S(k3km-1))
      ushr=u3km-uu_S(1)
      vshr=v3km-vv_S(1)
      helicity=sumh + vshr*ustrm - ushr*vstrm
c
      Return
c
      End
c
c	==============================================================
c	Function brnshear
c	==============================================================
c
      function brnsh(pres_S,zz_S,tempc_S,uu_S,vv_S,n,nmax,u6km,v6km)
c
c
      Real pres_S(nmax),zz_S(nmax),tempc_S(nmax)
      Real uu_S(nmax),vv_S(nmax)
      Integer n,nmax
c
      Real sumu,sumv,sump,tempk,tempklast
      Real rhohi,rholo,rhoinv,dp
      Real u500m,v500m,p500m,t500m
      Real u6km,v6km,p6km,t6km,wlow
c
c     Find the mass weighted mean wind in the first 500 meters
c
      sumu=0
      sumv=0
      sump=0
      tempklast=tempc_S(1)+273.16
      Do i=2,n,1
        tempk=tempc_S(i)+273.16
        If (zz_S(i) .LT. 500.) Then
          dp=pres_S(i-1)-pres_S(i)
          rhohi=pres_S(i)/tempk
          rholo=pres_S(i-1)/tempklast
          rhoinv=1./(rhohi+rholo)
          sumu=sumu+dp*(rhoinv)*(rholo*uu_S(i-1)+rhohi*uu_S(i))
          sumv=sumv+dp*(rhoinv)*(rholo*vv_S(i-1)+rhohi*vv_S(i))
          sump=sump+dp
        Else
          dz=zz_S(i)-zz_S(i-1)
          wlow=(zz_S(i)-500.)/dz
          u500m=uu_S(i)*(1.-wlow)+uu_S(i-1)*wlow
          v500m=vv_S(i)*(1.-wlow)+vv_S(i-1)*wlow
          p500m=pres_S(i)*(1.-wlow)+pres_S(i-1)*wlow
          t500m=tempk*(1.-wlow)+tempklast*wlow
          dp=pres_S(i-1)-p500m
          rhohi=p500m/t500m
          rholo=pres_S(i-1)/tempklast
          rhoinv=1./(rhohi+rholo)
          sumu=sumu+dp*rhoinv*(rholo*uu_S(i-1)+rhohi*u500m)
          sumv=sumv+dp*rhoinv*(rholo*vv_S(i-1)+rhohi*v500m)
          sump=sump+dp
          GO TO 150
        End If
      tempklast=tempk
      End Do
 150  Continue
      u500m=sumu/sump
      v500m=sumv/sump
c
c     Find the mass weighted mean wind in the first 6 km
c
      sumu=0
      sumv=0
      sump=0
      tempklast=tempc_S(1)+273.16
      Do i=2,n,1
        tempk=tempc_S(i)+273.16
        If (zz_S(i) .LT. 6000.) Then
          dp=pres_S(i-1)-pres_S(i)
          rhohi=pres_S(i)/tempk
          rholo=pres_S(i-1)/tempklast
          rhoinv=1./(rhohi+rholo)
          sumu=sumu+dp*(rhoinv)*(rholo*uu_S(i-1)+rhohi*uu_S(i))
          sumv=sumv+dp*(rhoinv)*(rholo*vv_S(i-1)+rhohi*vv_S(i))
          sump=sump+dp
        Else
          dz=zz_S(i)-zz_S(i-1)
          wlow=(zz_S(i)-6000.)/dz
          u6km=uu_S(i)*(1.-wlow)+uu_S(i-1)*wlow
          v6km=vv_S(i)*(1.-wlow)+vv_S(i-1)*wlow
          p6km=pres_S(i)*(1.-wlow)+pres_S(i-1)*wlow
          t6km=tempk*(1.-wlow)+tempklast*wlow
          dp=pres_S(i-1)-p6km
          rhohi=p6km/t6km*wlow
          dp=pres_S(i-1)-p6km
          rhohi=p6km/t6km
          rholo=pres_S(i-1)/tempklast
          rhoinv=1./(rhohi+rholo)
          sumu=sumu+dp*rhoinv*(rholo*uu_S(i-1)+rhohi*u6km)
          sumv=sumv+dp*rhoinv*(rholo*vv_S(i-1)+rhohi*v6km)
          sump=sump+dp
          GO TO 180
        End If
      tempklast=tempk
      End Do
 180  Continue
      u6km=sumu/sump
      v6km=sumv/sump
c
      brnsh=sqrt((u6km-u500m)**2+(v6km-v500m)**2)
c
	end
c
c	==============================================================
c	Function wobf(t)
c	==============================================================
c
      function wobf(t)
c
c   this function calculates the difference of the wet bulb potential
c   temperatures for saturated and dry air given the temperature.
c
c     baker,schlatter   17-may-1982 original version
c
c        let wbpts = wet-bulb potential temperature for saturated
c   air at temperature t (celsius). let wbptd = wet-bulb potential
c   temperature for completely dry air at the same temperature t.
c   the wobus function wobf (in degrees celsius) is defined by
c                      wobf(t) = wbpts-wbptd.
c   although wbpts and wbptd are functions of both pressure and
c   temperature, their difference is a function of temperature only.
c
c        to understand why, consider a parcel of dry air at tempera-
c   ture t and pressure p. the thermodynamic state of the parcel is
c   represented by a point on a pseudoadiabatic chart. the wet-bulb
c   potential temperature curve (moist adiabat) passing through this
c   point is wbpts. now t is the equivalent temperature for another
c   parcel saturated at some lower temperature tw, but at the same
c   pressure p.  to find tw, ascend along the dry adiabat through
c   (t,p). at a great height, the dry adiabat and some moist
c   adiabat will nearly coincide. descend along this moist adiabat
c   back to p. the parcel temperature is now tw. the wet-bulb
c   potential temperature curve (moist adiabat) through (tw,p) is wbptd.
c   the difference (wbpts-wbptd) is proportional to the heat imparted
c   to a parcel saturated at temperature tw if all its water vapor
c   were condensed. since the amount of water vapor a parcel can
c   hold depends upon temperature alone, (wbptd-wbpts) must depend
c   on temperature alone.
c
c        the wobus function is useful for evaluating several thermo-
c   dynamic quantities.  by definition:
c               wobf(t) = wbpts-wbptd.               (1)
c   if t is at 1000 mb, then t is a potential temperature pt and
c   wbpts = pt. thus
c               wobf(pt) = pt-wbptd.                 (2)
c   if t is at the condensation level, then t is the condensation
c   temperature tc and wbpts is the wet-bulb potential temperature
c   wbpt. thus
c               wobf(tc) = wbpt-wbptd.               (3)
c   if wbptd is eliminated from (2) and (3), there results
c               wbpt = pt-wobf(pt)+wobf(tc).
c   if wbptd is eliminated from (1) and (2), there results
c               wbpts = pt-wobf(pt)+wobf(t).
c
c        if t is an equivalent potential temperature ept (implying
c   that the air at 1000 mb is completely dry), then wbpts = ept
c   and wbptd = wbpt. thus
c               wobf(ept) = ept-wbpt.
c   this form is the basis for a polynomial approximation to wobf.
c   in table 78 on pp.319-322 of the smithsonian meteorological
c   tables by roland list (6th revised edition), one finds wet-bulb
c   potential temperatures and the corresponding equivalent potential
c   temperatures listed together. herman wobus, a mathematician for-
c   merly at the navy weather research facility, norfolk, virginia,
c   and now retired, computed the coefficients for the polynomial
c   approximation from numbers in this table.
c
c                                    notes by t.w. schlatter
c                                    noaa/erl/profs program office
c                                    august 1981
c
      x = t-20.
      if (x.gt.0.) go to 10
      pol = 1.                 +x*(-8.8416605e-03
     1           +x*( 1.4714143e-04  +x*(-9.6719890e-07
     2       +x*(-3.2607217e-08  +x*(-3.8598073e-10)))))
      wobf = 15.130/pol**4
      return
   10 continue
      pol = 1.                 +x*( 3.6182989e-03
     1       +x*(-1.3603273e-05  +x*( 4.9618922e-07
     2       +x*(-6.1059365e-09  +x*( 3.9401551e-11
     3       +x*(-1.2588129e-13  +x*( 1.6688280e-16)))))))
      wobf = 29.930/pol**4+0.96*x-14.8
      return
      end
c
c
c	==========================================================
c	Function satlft
c	==========================================================
c
      function satlft(thw,p)
c
c   input:  thw = wet-bulb potential temperature (celsius).
c             thw defines a moist adiabat.
c         p = pressure (millibars)
c   output: satlft = temperature (celsius) where the moist adiabat
c             crosses p
c
c     baker,schlatter   17-may-1982 original version
c
      data cta,akap/273.16,0.28541/
c   akap = difference between kelvin and celsius temperatures
c   akap = (gas constant for dry air) / (specific heat at constant
c         pressure for dry air)
c
c        the algorithm below can best be understood by referring to a
c   skew-t/log p chart.  it was devised by herman wobus, a mathemati-
c   cian formerly at the navy weather research facility but now retired.
c   the value returned by satlft can be checked by referring to table
c   78, pp.319-322, smithsonian meteorological tables, by roland list
c   (6th revised edition).
c
      if (p.ne.1000.) go to 5
      satlft = thw
      return
    5 continue
c   compute tone, the temperature where the dry adiabat with value thw
c   (celsius) crosses p.
      pwrp = (p/1000.)**akap
      tone = (thw+cta)*pwrp-cta
c   consider the moist adiabat ew1 through tone at p.  using the defini-
c   tion of the wobus function (see documentation on wobf), it can be
c   shown that eone = ew1-thw.
      eone = wobf(tone)-wobf(thw)
      rate = 1.
      go to 15
c   in the loop below, the estimate of satlft is iteratively improved.
   10 continue
c   rate is the ratio of a change in t to the corresponding change in
c   e.  its initial value was set to 1 above.
      rate = (ttwo-tone)/(etwo-eone)
      tone = ttwo
      eone = etwo
   15 continue
c   ttwo is an improved estimate of satlft.
      ttwo = tone-eone*rate
c   pt is the potential temperature (celsius) corresponding to ttwo at p
      pt = (ttwo+cta)/pwrp-cta
c   consider the moist adiabat ew2 through ttwo at p. using the defini-
c   tion of the wobus function, it can be shown that etwo = ew2-thw.
      etwo = pt+wobf(ttwo)-wobf(pt)-thw
c   dlt is the correction to be subtracted from ttwo.
      dlt = etwo*rate
      if (abs(dlt).gt.0.1) go to 10
      satlft = ttwo-dlt
      return
      end
c
c
c	=======================================================
c	Function powt
c	=======================================================
c
      function powt(t,p,td)
c
c   this function yields wet-bulb potential temperature powt
c   (celsius), given the following input:
c        t = temperature (celsius)
c        p = pressure (millibars)
c        td = dew point (celsius)
c
c     baker,schlatter   17-may-1982 original version
c
      data cta,akap/273.16,0.28541/
c        cta = difference between kelvin and celsius temperatures
c        akap = (gas constant for dry air) / (specific heat at
c              constant pressure for dry air)
c   compute the potential temperature (celsius)
      IF (p .LE. 0.0) PRINT *, 'POWT: ERROR: Bad p :', p
      pt = (t+cta)*(1000./p)**akap-cta
c   compute the lifting condensation level (lcl).
      tc = tcon(t,td)
c   for the origin of the following approximation, see the documen-
c   tation for the wobus function.
      powt = pt-wobf(pt)+wobf(tc)
      return
      end
c
c
c	=============================================================
c	Function tcon(t,d)
c	=============================================================
c
      function tcon(t,d)
c
c   this function returns the temperature tcon (celsius) at the lifting
c   condensation level, given the temperature t (celsius) and the
c   dew point d (celsius).
c
c     baker,schlatter   17-may-1982 original version
c
c   compute the dew point depression s.
      s = t-d
c   the approximation below, a third order polynomial in s and t,
c   is due to herman wobus. the source of data for fitting the
c   polynomial is unknown.
      dlt = s*(1.2185+1.278e-03*t+
     1        s*(-2.19e-03+1.173e-05*s-5.2e-06*t))
      tcon = t-dlt
      return
      end
c
c
c	===========================================================
c	Function XMXRAT
c	===========================================================
c
	       FUNCTION XMXRAT(PRES,DEWP)
c   COMPUTE MIXING RATIO (GM/GM) GIVEN DEW POINT TEMP
c   AND THE PRESSURE (MB)
       RATMIX=EXP(21.16-5415.0/DEWP)
       RATMIX=RATMIX/PRES
       IF(RATMIX.LT.(5.0E-05)) RATMIX=5.0E-05
       XMXRAT=RATMIX
       RETURN
       END
c
c	===========================================================
c	Function pccl
c	===========================================================
c
      function pccl(pm,p,t,td,wbar,nmax,n)
c
c   this function returns the pressure at the convective condensation
c   level given the appropriate sounding data.
c
c     baker,schlatter   17-may-1982 original version
c
c   on input:
c     p = pressure (millibars). note that p(i).gt.p(i+1).
c     t = temperature (celsius)
c     td = dew point (celsius)
c     n = number of levels in the sounding and the dimension of
c         p, t and td
c     pm = pressure (millibars) at upper boundary of the layer for
c          computing the mean mixing ratio. p(1) is the lower
c          boundary.
c   on output:
c     pccl = pressure (millibars) at the convective condensation level
c     wbar = mean mixing ratio (g/kg) in the layer bounded by
c            pressures p(1) at the bottom and pm at the top
c   the algorithm is decribed on p.17 of stipanuk, g.s.,1973:
c   "algorithms for generating a skew-t log p diagram and computing
c   selected meteorological quantities," atmospheric sciences labora-
c   tory, u.s. army electronics command, white sands missile range, new
c   mexico 88002.
      dimension t(nmax),p(nmax),td(nmax)
      if (pm.ne.p(1)) go to 5
      wbar= w(td(1),p(1))
      pc= pm
      if (abs(t(1)-td(1)).lt.0.05) go to 45
      go to 25
    5 continue
      wbar= 0.
      k= 0
   10 continue
      k = k+1
      if (p(k).gt.pm) go to 10
      k= k-1
      j= k-1
      if(j.lt.1) go to 20
c   compute the average mixing ratio....alog = natural log
      do 15 i= 1,j
         l= i+1
   15      wbar= (w(td(i),p(i))+w(td(l),p(l)))*alog(p(i)/p(l))
     *          +wbar
   20 continue
      l= k+1
c   estimate the dew point at pressure pm.
      tq= td(k)+(td(l)-td(k))*(alog(pm/p(k)))/(alog(p(l)/p(k)))
      wbar= wbar+(w(td(k),p(k))+w(tq,pm))*alog(p(k)/pm)
      wbar= wbar/(2.*alog(p(1)/pm))
c   find level at which the mixing ratio line wbar crosses the
c   environmental temperature profile.
   25 continue
      do 30 j= 1,n
         i= n-j+1
         if (p(i).lt.300.) go to 30
c   tmr = temperature (celsius) at pressure p (mb) along a mixing
c       ratio line given by wbar (g/kg)
         x= tmr(wbar,p(i))-t(i)
         if (x.le.0.) go to 35
   30 continue
      pccl= 0.0
      return
c  set up bisection routine
   35 l = i
      i= i+1
      del= p(l)-p(i)
      pc= p(i)+.5*del
      a= (t(i)-t(l))/alog(p(l)/p(i))
      do 40 j= 1,10
         del= del/2.
         x= tmr(wbar,pc)-t(l)-a*(alog(p(l)/pc))
c   the sign function replaces the sign of the first argument
c   with that of the second.
   40 pc= pc+sign(del,x)
   45 pccl = pc
      return
      end
c
c	=============================================================
c	Function w
c	=============================================================
c
      function w(t,p)
c
c   this function returns the mixing ratio (grams of water vapor per
c   kilogram of dry air) given the temperature t (celsius) and pressure
c   (millibars). the formula is quoted in most meteorological texts.
c
c     baker,schlatter   17-may-1982 original version
c
      x= esat(t)
      w= 622.*x/(p-x)
      return
      end
c
c	===============================================================
c	Function tmr
c	===============================================================
c
      function tmr(w,p)
c
c   this function returns the temperature (celsius) on a mixing
c   ratio line w (g/kg) at pressure p (mb). the formula is given in
c   table 1 on page 7 of stipanuk (1973).
c
c     baker,schlatter   17-may-1982 original version
c
c   initialize constants
      data c1/.0498646455/,c2/2.4082965/,c3/7.07475/
      data c4/38.9114/,c5/.0915/,c6/1.2035/
c
      x= alog10(w*p/(622.+w))
      tmrk= 10.**(c1*x+c2)-c3+c4*((10.**(c5*x)-c6)**2.)
      tmr= tmrk-273.16
      return
      end
c
c	==============================================================
c	Function o 
c	==============================================================
c
      function o(t,p)
c
c   this function returns potential temperature (celsius) given
c   temperature t (celsius) and pressure p (mb) by solving the poisson
c   equation.
c
c     baker,schlatter   17-may-1982 original version
c
      tk= t+273.16
      ok= tk*((1000./p)**.286)
      o= ok-273.16
      return
      end
c
c	==============================================================
c	Function tda
c	==============================================================
c
      function tda(o,p)
c
c   this function returns the temperature tda (celsius) on a dry adiabat
c   at pressure p (millibars).
c
c     baker,schlatter   17-may-1982 original version
c
c   the dry adiabat is given by
c   potential temperature o (celsius). the computation is based on
c   poisson's equation.
c
      ok= o+273.16
      tdak= ok*((p*.001)**.286)
      tda= tdak-273.16
      return
      end
c
c	===============================================================
c	Function ct
c	===============================================================
c
      function ct(wbar,pc,ps)
c
c   this function returns the convective temperature ct (celsius)
c   given the mean mixing ratio wbar (g/kg) in the surface layer,
c   the pressure pc (mb) at the convective condensation level (ccl)
c   and the surface pressure ps (mb).
c
c     baker,schlatter   17-may-1982 original version
c
c   compute the temperature (celsius) at the ccl.
      tc= tmr(wbar,pc)
c   compute the potential temperature (celsius), i.e., the dry
c   adiabat ao through the ccl.
      ao= o(tc,pc)
c        compute the surface temperature on the same dry adiabat ao.
      ct= tda(ao,ps)
      return
      end
c
c	======================================================
c	Function esat
c	======================================================
c
      function esat(t)
c
c   this function returns the saturation vapor pressure over
c   water (mb) given the temperature (celsius).
c
c     baker,schlatter   17-may-1982 original version
c
c   the algorithm is due to nordquist, w.s.,1973: "numerical approxima-
c   tions of selected meteorlolgical parameters for cloud physics prob-
c   lems," ecom-5475, atmospheric sciences laboratory, u.s. army
c   electronics command, white sands missile range, new mexico 88002.
c
      tk = t+273.16
      p1 = 11.344-0.0303998*tk
      p2 = 3.49149-1302.8844/tk
      c1 = 23.832241-5.02808*alog10(tk)
      esat = 10.**(c1-1.3816e-7*10.**p1+8.1328e-3*10.**p2-2949.076/tk)
      return
      end
c
c	========================================================
c	Function precpw
c	========================================================
c
      function precpw(td,p,n)
c
c   this function computes total precipitable water precpw (cm) in a
c   vertical column of air based upon sounding data at n levels:
c        td = dew point (celsius)
c        p = pressure (millibars)
c
c     baker,schlatter   17-may-1982 original version
c
c   calculations are done in cgs units.
      dimension td(n),p(n)
c   g = acceleration due to the earth's gravity (cm/s**2)
      data g/980.616/
c   initialize value of precipitable water
      pw = 0.
      nl = n-1
c   calculate the mixing ratio at the lowest level.
      wbot = wmr(p(1),td(1))
      do 5 i=1,nl
      wtop = wmr(p(i+1),td(i+1))
c   calculate the layer-mean mixing ratio (g/kg).
      w = 0.5*(wtop+wbot)
c   make the mixing ratio dimensionless.
      wl = .001*w
c   calculate the specific humidity.
      ql = wl/(wl+1.)
c   the factor of 1000. below converts from millibars to dynes/cm**2.
      dp = 1000.*(p(i)-p(i+1))
      pw = pw+(ql/g)*dp
      wbot = wtop
    5 continue
      precpw = pw
      return
      end
c
c	==========================================================
c	Function wmr
c	==========================================================
c
      function wmr(p,t)
c
c   this function approximates the mixing ratio wmr (grams of water
c   vapor per kilogram of dry air) given the pressure p (mb) and the
c   temperature t (celsius).
c
c     baker,schlatter   17-may-1982 original version
c
c     the formula used is given on p. 302 of the
c   smithsonian meteorological tables by roland list (6th edition).
c
c   eps = ratio of the mean molecular weight of water (18.016 g/mole)
c         to that of dry air (28.966 g/mole)
      data eps/0.62197/
c   the next two lines contain a formula by herman wobus for the
c   correction factor wfw for the departure of the mixture of air
c   and water vapor from the ideal gas law. the formula fits values
c   in table 89, p. 340 of the smithsonian meteorological tables,
c   but only for temperatures and pressures normally encountered in
c   in the atmosphere.
      x = 0.02*(t-12.5+7500./p)
      wfw = 1.+4.5e-06*p+1.4e-03*x*x
      fwesw = wfw*esw(t)
      r = eps*fwesw/(p-fwesw)
c   convert r from a dimensionless ratio to grams/kilogram.
      wmr = 1000.*r
      return
      end
c
c	==============================================================
c	Function esw
c	==============================================================
c
      function esw(t)
c
c   this function returns the saturation vapor pressure esw (millibars)
c   over liquid water given the temperature t (celsius).
c
c     baker,schlatter   17-may-1982 original version
c
c   the polynomial approximation below is due to herman wobus, a mathematician
c   who worked at the navy weather research facility, norfolk, virginia,
c   but who is now retired. the coefficients of the polynomial were
c   chosen to fit the values in table 94 on pp. 351-353 of the smith-
c   sonian meteorological tables by roland list (6th edition). the
c   approximation is valid for -50 < t < 100c.
c
c   es0 = saturation vapor ressure over liquid water at 0c
      data es0/6.1078/
      pol = 0.99999683       + t*(-0.90826951e-02 +
     1         t*(0.78736169e-04   + t*(-0.61117958e-06 +
     2     t*(0.43884187e-08   + t*(-0.29883885e-10 +
     3     t*(0.21874425e-12   + t*(-0.17892321e-14 +
     4     t*(0.11112018e-16   + t*(-0.30994571e-19)))))))))
      esw = es0/pol**8
      return
      end
c
c	================================================================
c	Function tw
c	================================================================
c
      function tw(t,td,p)
c
c   this function returns the wet-bulb temperature tw (celsius)
c   given the temperature t (celsius), dew point td (celsius)
c   and pressure p (mb).  see p.13 in stipanuk (1973), referenced
c   above, for a description of the technique.
c
c     baker,schlatter   17-may-1982 original version
c
c   determine the mixing ratio line thru td and p.
      aw = w(td,p)
c
c   determine the dry adiabat thru t and p.
      ao = o(t,p)
      pi = p
c
c   iterate to locate pressure pi at the intersection of the two
c   curves .  pi has been set to p for the initial guess.
      do 4 i= 1,10
         x= .02*(tmr(aw,pi)-tda(ao,pi))
         if (abs(x).lt.0.01) go to 5
 4       pi= pi*(2.**(x))
c   find the temperature on the dry adiabat ao at pressure pi.
 5    ti= tda(ao,pi)
c
c   the intersection has been located...now, find a saturation
c   adiabat thru this point. function os returns the equivalent
c   potential temperature (c) of a parcel saturated at temperature
c   ti and pressure pi.
      aos= os(ti,pi)
c   function tsa returns the wet-bulb temperature (c) of a parcel at
c   pressure p whose equivalent potential temperature is aos.
      tw = tsa(aos,p)
      return
      end
c
c	=========================================================
c	Function os
c	=========================================================
c
      function os(t,p)
c
c   this function returns the equivalent potential temperature os
c   (celsius) for a parcel of air saturated at temperature t (celsius)
c   and pressure p (millibars).
c
c     baker,schlatter   17-may-1982 original version
c
      data b/2.6518986/
c   b is an empirical constant approximately equal to the latent heat
c   of vaporization for water divided by the specific heat at constant
c   pressure for dry air.
      tk = t+273.16
      osk= tk*((1000./p)**.286)*(exp(b*w(t,p)/tk))
      os= osk-273.16
      return
      end
c
c	===========================================================
c	Function tsa
c	===========================================================
c
      function tsa(os,p)
c
c   this function returns the temperature tsa (celsius) on a saturation
c   adiabat at pressure p (millibars). os is the equivalent potential
c   temperature of the parcel (celsius). sign(a,b) replaces the
c   algebraic sign of a with that of b.
c
c     baker,schlatter   17-may-1982 original version
c
c   b is an empirical constant approximately equal to the latent heat
c   of vaporization for water divided by the specific heat at constant
c   pressure for dry air.
c
      data b/2.6518986/
      a= os+273.16
c   tq is the first guess for tsa.
      tq= 253.16
c   d is an initial value used in the iteration below.
      d= 120.
c
c   iterate to obtain sufficient accuracy....see table 1, p.8
c   of stipanuk (1973) for equation used in iteration.
      do 1 i= 1,12
         tqk= tq-273.16
         d= d/2.
         x= a*exp(-b*w(tqk,p)/tq)-tq*((1000./p)**.286)
         if (abs(x).lt.1e-7) go to 2
         tq= tq+sign(d,x)
 1    continue
 2    tsa= tq-273.16
      return
      end
c
c	==============================================================
c	Function oe
c	==============================================================
c
      function oe(t,td,p)
c
c   this function returns equivalent potential temperature oe (celsius)
c   of a parcel of air given its temperature t (celsius), dew point
c   td (celsius) and pressure p (millibars).
c
c     baker,schlatter   17-may-1982 original version
c
c   find the wet bulb temperature of the parcel.
      atw = tw(t,td,p)
c   find the equivalent potential temperature.
      oe = os(atw,p)
      return
      end
c
c	=======================================================================
c	Function tmlaps
c	=======================================================================
c
      function tmlaps(thetae,p)
c
c   this function returns the temperature tmlaps (celsius) at pressure
c   p (millibars) along the moist adiabat corresponding to an equivalent
c   potential temperature thetae (celsius).
c
c     baker,schlatter   17-may-1982 original version
c
c   the algorithm was written by eric smith at colorado state
c   university.
      data crit/0.1/
c   cta = difference between kelvin and celsius temperatures.
c   crit = convergence criterion (degrees kelvin)
      eq0 = thetae
c   initial guess for solution
      tlev = 25.
c   compute the saturation equivalent potential temperature correspon-
c   ding to temperature tlev and pressure p.
      eq1 = ept(tlev,tlev,p)
      dif = abs(eq1-eq0)
      if (dif.lt.crit) go to 3
      if (eq1.gt.eq0) go to 1
c   dt is the initial stepping increment.
      dt = 10.
      i = -1
      go to 2
    1 dt = -10.
      i = 1
    2 tlev = tlev+dt
      eq1 = ept(tlev,tlev,p)
      dif = abs(eq1-eq0)
      if (dif.lt.crit) go to 3
      j = -1
      if (eq1.gt.eq0) j=1
      if (i.eq.j) go to 2
c   the solution has been passed. reverse the direction of search
c   and decrease the stepping increment.
      tlev = tlev-dt
      dt = dt/10.
      go to 2
    3 tmlaps = tlev
      return
      end
c
c	=================================================================
c	Function ept
c	=================================================================
c
      function ept(t,td,p)
c
c   this function returns the equivalent potential temperature ept
c   (celsius) for a parcel of air initially at temperature t (celsius),
c   dew point td (celsius) and pressure p (millibars).
c
c     baker,schlatter   17-may-1982 original version
c
c     the formula used
c   is eq.(43) in bolton, david, 1980: "the computation of equivalent
c   potential temperature," monthly weather review, vol. 108, noc   potential temperature," monthly weather review, vol. 108, no. 7
c   (july), pp. 1046-1053. the maximum error in ept in 0.3c.  in most
c   cases the error is less than 0.1c.
c
c   compute the mixing ratio (grams of water vapor per kilogram of
c   dry air).
      w = wmr(p,td)
c   compute the temperature (celsius) atr).
      w = wmr(p,td)
c   compute the temperature (celsius) at the lifting condensation level.
      tlcl = tcon(t,td)
      tk = t+273.16
      tl = tlcl+273.16
      pt = tk*(1000./p)**(0.2854*(1.-0.00028*w))
      eptk = pt*exp((3.376/tl-0.00254)*w*(1.+0.00081*w))
      ept= eptk-273.16
      return
      end
