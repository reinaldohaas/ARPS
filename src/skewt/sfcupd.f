      SUBROUTINE SFCUPD(pres,hgt,temp,dewpt,q,theta,
     +                  u,v,sfctf,sfcdpf,orig,nlev,maxlev, 
     >  kmod,pres_mod,orgtf,orgdpf,weightq)
      IMPLICIT NONE
      Logical weightq
C
C  Mixes out lower portion of thermodynamic profile based on
C  updated surface conditions provided by user.
C
C  Important note: nlev can be changed by this routine.
C                  make sure all variables are updated so 
C                  they match reported heights/press.
C
C  Keith Brewster, March 1992     
C  OU SoM
c
c  Modified Nov 1992, Richard Carpenter, Univ. of Okla.
c    weightq:  If true, weight modified mixing ratio. If false, set it
c		   const in the modified layer.
C
C  Arguments
C
      INTEGER maxlev
      REAL pres(maxlev),hgt(maxlev),temp(maxlev),dewpt(maxlev),
     +     q(maxlev),u(maxlev),v(maxlev),theta(maxlev)
      INTEGER nlev, kmod
      LOGICAL orig
C
C  Constants
C
      REAL RCp
      PARAMETER (RCp = 0.286)
C
C  Functions
C
      REAL DWPTOQ,PR_DWPT
C
C  Misc internal variables
C
      REAL orgtf,sfctf,sfctk,sfcth, pres_mod
      REAL wgt,qsum,wgtsum,dpavg,qavg,qmix
      REAL orgdpf,sfcdpf,sfcdpc
      REAL w1,w2,dth,elen
      INTEGER nn,k,km1,kp1,ktop
C
C  Get temperature update from user
C
      orgtf=(temp(1)*9./5.) + 32.
      WRITE(6,805) orgtf,' F', temp(1),' C'
 805  FORMAT(' Sounding sfc temp:', 2(f6.1,a),
     +     /,' Enter new sfc temp (F), -99 for no change: ')
      read(5,*) sfctf
      IF(sfctf .GT. -98.) THEN     
        orig=.false.
C
C  Convert F to Kelvin
C
        sfctk=((sfctf - 32.) * (5./9.)) + 273.15
C
C  Find sfc theta
C      
        sfcth=sfctk*((1000./pres(1))**RCp)
        IF(sfcth.GT.theta(1)) THEN
C
C  Search through sounding to found where sfc theta intersects sounding
C
c  kmod=ktop is the top of the modified layer. (RLC 11/13/92)
c
          nn=nlev
          DO k=2,nlev
            IF(theta(k).GT.sfcth) THEN
C
C  Make room for new level...move what's above this
C       level up one index
C
              CALL MOVLEV(pres,k,nlev,maxlev)
              CALL MOVLEV(hgt,k,nlev,maxlev)
              CALL MOVLEV(temp,k,nlev,maxlev)
              CALL MOVLEV(theta,k,nlev,maxlev)
              CALL MOVLEV(dewpt,k,nlev,maxlev)
              CALL MOVLEV(q,k,nlev,maxlev)
              CALL MOVLEV(u,k,nlev,maxlev)
              CALL MOVLEV(v,k,nlev,maxlev)
C
C  Note through the action of MOVLEV, what was at level k 
C       is now at kp1
C
              ktop=k
              km1=k-1
              kp1=k+1
              dth=theta(kp1)-theta(k-1)
              w1=(sfcth-theta(k-1))/dth
              w2=1.-w1
              theta(k)=sfcth
              pres(k )=w2*pres(km1)  + w1*pres(kp1)
              hgt(k ) =-9999.
              q(k)    =w2*q(km1)     + w1*q(kp1)
              dewpt(k)=PR_DWPT(q(k),pres(k))
              u(k)    =w2*u(km1)     + w1*u(kp1)
              v(k)    =w2*v(km1)     + w1*v(kp1)
              temp(k)=(theta(k)*((pres(k)/1000.)**RCp)) - 273.15
              IF(dewpt(k).GT.temp(k)) THEN
                dewpt(k)=temp(k)
                q(k) = DWPTOQ(pres(k),temp(k),dewpt(k))
              END IF
              nn=nlev+1
              GO TO 25
            END IF
          END DO
 25       CONTINUE
          nlev=nn
c
c  Update modified levels (do loop moved RLC 01/93)
c
	    Do k=2,ktop-1
              theta(k)=sfcth
              temp(k)=(theta(k)*((pres(k)/1000.)**RCp)) - 273.15
              hgt(k)=-9999.
              IF(dewpt(k).GT.temp(k)) THEN
                dewpt(k)=temp(k)
                q(k) = DWPTOQ(pres(k),temp(k),dewpt(k))
              END IF
	    End Do
        ELSE
          ktop=1
        END IF
C
C  Update surface temperature itself
C
        theta(1)=sfcth
        temp(1)=sfctk-273.15
        IF(dewpt(1).GT.temp(1)) THEN
          dewpt(1)=temp(1)
          q(1) = DWPTOQ(pres(1),temp(1),dewpt(1))
        END IF
      ELSE  ! no temperature change at all
        ktop=1
      END IF
C
C  Next, attend to the moisture update.
C  Find mean Q in mixed layer.
C
      IF(ktop.GT.1) THEN
        wgt=0.5*(pres(1)-pres(2))
        qsum=wgt*q(1)
        wgtsum=wgt
	k	= 1
        DO k=2,(ktop-1)
          wgt=0.5*(pres(k-1)-pres(k+1))
          qsum=qsum+wgt*q(k)
          wgtsum=wgtsum+wgt
        END DO
        wgt=0.5*(pres(ktop-1)-pres(ktop))	!RLC 01/93
        qsum=qsum+wgt*q(ktop)
        wgtsum=wgtsum+wgt
	k	= ktop
        IF(wgtsum.NE. 0.) THEN
          qavg=qsum/wgtsum
        ELSE
          qavg=0.
        END IF
      ELSE
        qavg=q(1)
      END IF
C
C  report average depwt in mixed layer as a surface dew-point
C  in degrees F, and get user-updated surface dewpt
C
      pres_mod	= pres(ktop)
      dpavg=PR_DWPT(qavg,pres(1))
      orgdpf=(dpavg*9./5.) + 32.
      WRITE(6,810) pres(ktop),orgdpf,' F',dpavg,' C'
  810 FORMAT('  Mixed layer extends to',F8.1,' mb',/
     +'  Sfc dew-point from avg mixing ratio in mixed layer:',2(f6.1,a))

      If (ktop.NE.1) Then
        Print 811
  811 Format ('  Enter dew-point representative of mixed layer (F),',
     +/'  Use -99. for no change in moisture.')
      Else
        Print *, 'No mixed layer'
        Return
      End If
C
      read(5,*) sfcdpf
      IF(sfcdpf .GT. -98.) THEN
        orig=.false.
C
C  Find mixing ratio corresponding to that dewpt
C
        sfcdpc=(sfcdpf - 32.) * (5./9.)
        qmix= DWPTOQ(pres(1),temp(1),sfcdpc)
C
C  Set all mixing ratios in mixed layer to be equal to that q
C
	If (.NOT. weightq) Then
	Print *, '% Sfcupd: Mixing ratios in mixed layer set constant'
        DO k=1,ktop
          q(k)=qmix
        END DO
C
C  Instead, change mixing ratio to be a weighted mean of the
C  specified sfc mixing ratio and the previously observed value.
C  The e-length for the weighting is a function of the depth of the
C  mixed layer for deep mixed layers (> 300. mb depth) mixed layer.
C
      Else
	Print *, '% Sfcupd: Mixing ratios in mixed layer are weighted'
      elen=MAX(300.,(pres(1)-pres(ktop)))
      elen=0.5*elen
      q(1)=qmix
      DO k=2,ktop
        wgt=exp((pres(k)-pres(1))/elen)
        q(k)=(wgt*qmix) + ((1.-wgt)*q(k))
      END DO
      End If
      
C
C  Convert new mixing ratios to dewpts
C
        CALL CALDEWP(pres,q,dewpt,ktop,maxlev)
      END IF
C
C  Check for consistency with temperature profile
C
      DO k=1,ktop
        IF(dewpt(k).GT.temp(k)) THEN
          dewpt(k) = (temp(k) - 0.2)   ! not quite saturated
          q(k) = DWPTOQ(pres(k),temp(k),dewpt(k))
        END IF
      END DO
C
C  Done
c
      kmod	= ktop
      If (sfctf.EQ.-99.) sfctf = orgtf
      If (sfcdpf.EQ.-99.) sfcdpf = orgdpf
C
      RETURN
      END
C
      SUBROUTINE MOVLEV(Z,k,nlev,maxlev)
      IMPLICIT NONE
C
C  Arguments
C
      INTEGER k,nlev,maxlev
      REAL Z(maxlev)
C
C  Misc internal variables
C
      INTEGER i
C
      DO i=nlev,k,-1
        Z(i+1)=Z(i)
      END DO
      RETURN
      END
