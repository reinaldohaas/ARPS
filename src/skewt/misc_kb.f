      SUBROUTINE CALTHET(pres,temp,theta,nlev,maxlev)
      IMPLICIT NONE
C
C  Arguments
C
      INTEGER maxlev,nlev
      REAL pres(maxlev),temp(maxlev),theta(maxlev)
C
C Misc internal variables
C
      REAL RCp
      PARAMETER (RCp=0.286)
      INTEGER k
C
      DO k=1,nlev
        theta(k)=(temp(k)+273.15)*((1000./pres(k))**RCp)
      END DO
      RETURN
      END
C
      SUBROUTINE CALUV(drct,sped,u,v,nlev,maxlev)
      IMPLICIT NONE
C
C  Arguments
C
      INTEGER maxlev,nlev
      REAL drct(maxlev),sped(maxlev),u(maxlev),v(maxlev)
C
C Misc internal variables
C
      INTEGER k
C
      DO k=1,nlev
        CALL GET_UV(drct(k),sped(k),u(k),v(k))
      END DO
      RETURN
      END
C
      SUBROUTINE CALDDFF(u,v,drct,sped,nlev,maxlev)
      IMPLICIT NONE
C
C  Arguments
C
      INTEGER maxlev,nlev
      REAL drct(maxlev),sped(maxlev),u(maxlev),v(maxlev)
C
C Misc internal variables
C
      INTEGER k
C
      DO k=1,nlev
        CALL GET_DDFF(u(k),v(k),drct(k),sped(k))
      END DO
      RETURN
      END
C
      SUBROUTINE CALQ(pres,temp,dewp,q,nlev,maxlev)
      IMPLICIT NONE
C
C  Arguments
C
      INTEGER maxlev,nlev
      REAL pres(maxlev),temp(maxlev),dewp(maxlev),q(maxlev)
C
C Misc internal variables
C
      INTEGER k
      REAL DWPTOQ
C
      DO k=1,nlev
        q(k)=DWPTOQ(pres(k),temp(k),dewp(k))
      END DO
      RETURN
      END
 
       SUBROUTINE CALDEWP(pres,q,dewp,nlev,maxlev)
       IMPLICIT NONE
C
C  Computes dew-point (C) from mixing ratio (g/kg) and
C  pressure (mb) for all sounding levels.
C
C  Uses gempak function PR_DWPT for the conversion.
C
C  Arguments
C
       INTEGER maxlev,nlev
       REAL pres(maxlev),dewp(maxlev),q(maxlev)
C
C Misc internal variables
C
       INTEGER k
       REAL PR_DWPT
 
       DO k=1,nlev
         IF(pres(k).GT.0. .AND. q(k).GT.0.) THEN
           dewp(k)=PR_DWPT(q(k),pres(k))
         ELSE
           dewp(k)=-9999.
         END IF
       END DO
       RETURN
       END
C
      SUBROUTINE GET_UV(DD,FF,U,V)
      IMPLICIT NONE
      REAL DD,FF,U,V,DEG2R,SPVAL,MIS_VAL
      PARAMETER (DEG2R=0.0174532925, SPVAL=-9998., 
     +           MIS_VAL=-9999.0)
      IF(DD.GT.SPVAL .AND. FF.GT.SPVAL) THEN
        U=-FF*SIN(DEG2R*DD)
        V=-FF*COS(DEG2R*DD)
      ELSE
        U=MIS_VAL
        V=MIS_VAL
      END IF
      RETURN
      END
C
      SUBROUTINE GET_DDFF(U,V,DD,FF)
      IMPLICIT NONE
      REAL DD,FF,U,V,RAD2D,SPVAL,MIS_VAL
      PARAMETER (RAD2D=57.29577951, SPVAL=-9998., 
     +           MIS_VAL=-9999.0)
      IF(U.GT.SPVAL .AND. V.GT.SPVAL) THEN
        FF = SQRT((U*U + V*V))
        IF(FF.NE.0.) THEN
          DD = RAD2D*ATAN2(U,V)
          DD = DD+180.
          IF (DD.GT.360.) DD=DD-360.
        ELSE
          DD=0.
        END IF
      ELSE
        DD = MIS_VAL
        FF = MIS_VAL
      END IF
      RETURN
      END
C
      FUNCTION DWPTOQ(PRESS,T,DWPT)
C
C  Calculates mixing ratio (g/kg) given pressure (mb) and
C  dewpoint (C).
C  
C  If dew point is missing (.LT.-90.) or the dew point depression
C  is reported as 30 degrees, then the QV is calculated as 20 percent
C  of the saturation QV.
C  
C  If temp and dew point are missing, QV=0.
C
C  Version for surface data Keith Brewster April, 1991     OU SoM
C  ARPS version             Keith Brewster February, 1992
C  Modsnd version           Keith Brewster March, 1992
C  Based on routines documented in GEMPAK manual
C
C
      IMPLICIT NONE
C
C  Function
C
      REAL DWPTOQ
      REAL PRESS,T,DWPT
C
C  Misc internal variables
C
      REAL DEPRS,VAPR,E,QVSAT
C
      DEPRS=T-DWPT
      IF(DWPT.GT.-90. .AND. DEPRS.LT.29.999) THEN
        vapr=6.112 * exp((17.67 * dwpT)/ (dwpT + 243.5))
        e = vapr *( 1.001 + (PRESS-100.)/ 900. * .0034)
        DWPTOQ = .62197 *( e/(PRESS - e)) * 1000.
      ELSE IF(T.GT.-90.) THEN
        vapr=6.112 * exp((17.67 * T)/ (T + 243.5))
        e = vapr *( 1.001 + (PRESS-100.)/ 900. * .0034)
        QVSAT = .62197 *( e/(PRESS - e)) * 1000.
        DWPTOQ=0.2*QVSAT
      ELSE
        DWPTOQ=0.
      END IF
      DWPTOQ=AMAX1(DWPTOQ,0.)
C      print *, '  PRESS, TEMP, DEWPT = ',PRESS,T,DWPT
C      print *, '      Yer mixing ratio is: ',DWPTOQ
      RETURN
      END
C
c     SUBROUTINE CALHGT(pres,hght,temp,dewp,
c    +                  selev,nlev,maxlev)
c
C  fills in missing heights by integrating the hydrostatic eqn
c  uses GEMPAK subroutines where possible
c
c     IMPLICIT NONE
c
C  Arguments
C
c     INTEGER nlev,maxlev
c     REAL pres(maxlev),hght(maxlev),temp(maxlev),dewp(maxlev)
c     REAL selev
C
C  GEMPAK functions
C
c     REAL PR_SCLH,PR_MHGT
c     External PR_SCLH,PR_MHGT
C
C  Misc internal variables
C
c     INTEGER k
c     REAL tdb,tdt,sclh
C
c    
c     hght(1)=selev
c     DO k=2,nlev
c       IF(hght(k).LT.selev) THEN
c         tdb=amax1(dewp(k-1),(temp(k-1)-30.))
c         tdt=amax1(dewp(k),(temp(k)-30.))
c         SCLH=PR_SCLH(temp(k-1),temp(k),tdb,tdt,pres(k-1),pres(k))
c         hght(k)=PR_MHGT(hght(k-1),pres(k-1),pres(k),sclh)
c       END IF
c     END DO
c     RETURN
c     END
c
c
c
c
      Real Function PR_DWPT (qvgkg,presmb)
      Implicit None
      Include 'thermo.consts'
      Real	qvgkg,presmb
      Include 'thermo.stfunc'
c
c  Convert qv,p -> Td. Meteorological units.
c
      PR_DWPT = Ftdewc(Fvpres(presmb*100.,qvgkg*1.e-3))
      End

