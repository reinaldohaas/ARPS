!
!#######################################################################
!#######################################################################
!####                                                               ####
!####                   READSOUND                                   ####
!####                                                               ####
!#######################################################################
!#######################################################################
!
      Subroutine ReadSound (pres_S, zz_S, theta_S, qv_S, uu_S, vv_S, 
     >  snd_pres, snd_zz, pmbCB, tcCB, dz, stnelev,psfc,
     >  sndfile, nmax, n, ifCB, plot_wind, gempak_format, snd_wind,
     >  stid,stnm,slat,slon, iyr,imon,iday,ihr,imin, istatus)

!-----------------------------------------------------------------------
!
! 1994/11/14  Added code to keep input qv >= 0.
! 1995/02/15  New GEMPAK 5.0 has extra line. Both fmts supported here
! 1997/03/06  Change GEMPAK format istnm (INTEGER) to stnm (CHAR).
! 1997/04/22  Heights are now MSL, not AGL.
! 1999/10/11  Added slat,slon. [RLC]
! 2007/09/20  Added extra check for wmo > in FSL file type detection. [KB]
!
!-----------------------------------------------------------------------
!
!  Input: 
!     sndfile
!
!  Output:
!     pres_S, zz_S, theta_S, qv_S = Pres, height, theta, mixing ratio
!     snd_pres                    = True if sndfile has pressure
!     snd_zz                      = True if sndfile has heights
!     ifCB                        = 1 if sndfile contains cloud base information
!     pmbCB, tcCB                 = Pres (mb) and T (C) of cloud base (if ifCB = 1)
!     dz                   = Spacing between levels, if known and const
!     stnelev              = Station elevation
!     psfc                 = Pres of surface
!     nmax                 = Max size of arrays
!     n                    = number of points in sndfile
!     plot_wind            = True if sndfile has winds and winds are non-zero.
!     snd_wind             = True if sndfile has winds
!     gempak_format        = True if sndfile is in Gempak format
!     stid,stnm,iyr...imin = Read in for Gempak format
!     slat,slno            = stn lat and lon (GEMPAK only)
!
!
!  Sounding File Formats
!
!  A. Non-Spiffy formats
!     1. 'SAM' format (theta, qv)
!     2a. 'Raw' format (p, T, Td)
!     2b. 'Gempak' format (p, T, Td)
!     2c. 'FSL' format (p, T, Td)
!
!  B. Spiffy Formats
!     01. p, T, Td
!     02. z, theta, qv
!
!-----------------------------------------------------------------------
!
      IMPLICIT NONE

      INCLUDE 'thermo.consts'

      INTEGER            :: nmax, ifCB, n, k, iyr,imon,iday,ihr,imin,   &
     &                      nwords, j,j1,j2,jformat, k1,k2,kinc, iter

      INTEGER, PARAMETER :: nwordmax = 4, nwrk = 1000
      REAL               :: pmbCB, tcCB, dz, stnelev,psfc, xx, temp, rh,&
     &                      tv,tvkm1,tvavg,tdew, slat,slon
      REAL               :: pres_S(nmax), zz_S(nmax), theta_S(nmax),    &
     &                      qv_S(nmax), uu_S(nmax), vv_S(nmax),         &
     &                      tempc_S(nwrk), tdewc_S(nwrk)
      REAL               :: wspd,wdir,spcvt,dtr
      REAL,    PARAMETER :: kts2ms = 0.514444
      LOGICAL            :: snd_pres, snd_zz, plot_wind, gempak_format, &
     &                      snd_wind, minus_z
      Logical   fsl_format
      Integer idata(6)
      Integer iunit,ios
      Integer lintyp,iline,sonde,wmo,wban,ielev,ipsfc,itime
      Integer hydro,mxwd,tropl,lines,nlines,tindex,source
      Integer ktop,kbot,kk,kn,nlevs
      Real    dtdz,dtddz,dudz,dvdz
      Logical foundt,foundw

      Character*4 monlist(12)
      Data monlist /' JAN',' FEB',' MAR',' APR',' MAY',' JUN',
     :              ' JUL',' AUG',' SEP',' OCT',' NOV',' DEC'/

      Character*(*) sndfile,stid,stnm
      Character*49  line2
      Character*64  line,word(nwordmax)
      Character*24  cheight,ctemp,cmoist
      Character*4   chmon
      Character*3   staid
      Character*2   wsunits
      Character*1   ns,ew

      INTEGER, INTENT(OUT) :: istatus

      INCLUDE  'thermo.stfunc'
!
!----------------------------------------------------------------------=
!
      dtr           = acos(-1.)/180.
      ifCB          = 0
      jformat       = -1
      snd_pres      = .False.
      snd_zz        = .False.
      gempak_format = .False.
      fsl_format    = .False.
      stnelev       = 0.
!
      CALL ChkExist (sndfile,'*** ')
      PRINT *
      PRINT *, 'Reading from: ', TRIM(sndfile)
      iunit=31
      OPEN (iunit, File=sndfile, Status='Old')
!-----------------------------------------------------------------------
!
!...Read the first line, Check for my format
!
!-----------------------------------------------------------------------
!
      READ (iunit,'(a64)') line
      PRINT *, 'Reading first line - ',TRIM(line)
!
!-----------------------------------------------------------------------
!
!  Parsing the file header
!
!-----------------------------------------------------------------------
      !
      ! Format begins with '%'
      !
      IF (line(1:1).EQ.'%') Then
        PRINT '(a)', line
        j1        = 1
        j2        = Index(line(2:),'%') + j1
        Read (line(2:),'(I2)') jformat
        PRINT '(a,i2.2)', '% Sounding format: ', jformat
        nwords        = 0
        IF (j2.GT.4) CALL Parse (line(4:j2-1),word,nwords,nwordmax)
        PRINT *, nwords, line(4:j2-1)
        DO j=1,nwords
          PRINT *, j, word(j)
        END Do
        DO j=1,nwords
          CALL Uprc (word(j))
          PRINT *, 'Sounding flag: ', word(j)
          IF (word(j).EQ.'W') Then
            snd_wind        = .True.
          ELSE IF (word(j).EQ.'-Z') Then
            minus_z        = .True.
          ELSE IF (word(j).EQ.'CB') Then
            ifCB        = 1
          ELSE 
            PRINT *, 'Unrecognized flag: ', word(j)
            Stop
          END IF
        END DO

      ELSE IF (line(1:1).EQ.'&') Then
      !
      ! Format begins with '&' -- ARPS input format
      !
        PRINT *, 'ARPS input format'
        jformat        = -99
        snd_wind       = .True.
        minus_z        = .True.
        READ (iunit,'(a)') line
        PRINT *, line
        READ (iunit,'(a)') line
        PRINT *, line
        READ (iunit,'(a)') line
        PRINT *, line
        READ (iunit,'(a)') line
        PRINT *, line
        READ (iunit,*) cheight,ctemp,cmoist
        CALL UPRC (cheight)
        CALL UPRC (ctemp)
        CALL UPRC (cmoist)

        IF (cheight(1:1).EQ.'H') Then
          snd_zz = .True.
        ELSE IF (cheight(1:1).EQ.'P') Then
          snd_pres = .True.
        ELSE
          PRINT *, 'Unknown vertical coord: ', cheight
          STOP
        END IF
        PRINT *, 'Vertical coord: ', cheight
        !
        !  1=Pot.temp., 2=Temp
        !
        IF (ctemp(1:1).NE.'P' .AND. ctemp(1:1).NE.'T') Then
          PRINT *, 'Unknown temperature variable: ', ctemp
          Stop
        END IF
        PRINT *, 'Temperature variable: ', ctemp
        !
        !  1=Mixing ratio, 2=RH, 3=Td.
        !
        IF (cmoist(1:1).NE.'M' .AND. cmoist(1:1).NE.'R'
     >                         .AND. cmoist(1:1).NE.'D') Then
          PRINT *, 'Unknown moisture variable: ', cmoist
          Stop
        END IF
        PRINT *, 'Moisture variable: ', cmoist

        Read (iunit,*) stnelev, psfc
        Read (iunit,*) n
        Read (iunit,'(a)') line
      END IF

!-----------------------------------------------------------------------
!
!  Read soundings,   So far, the first two lines have been read
!
!  jformat > 0    First character is "%"
!  jformat = -99  First character is "&"
!  jformat = -1   Others
!
!  Older or Odd-ball formats
!
!----------------------------------------------------------------------
      IF (jformat == -99) Then
        ! ARPS Input Format

        PRINT *, 'Reading ARPS format sounding'
        k1        = 1
        k2        = n
        kinc      = 1
        IF (minus_z) Then
          k1   = n
          k2   = 1
          kinc = -1
        END IF
        Do k=k1,k2,kinc
          Read(iunit,*) zz_S(k),theta_S(k),qv_S(k),uu_S(k),vv_S(k)
          qv_S(k) = MAX( 0.0, qv_S(k) )
          PRINT *, k,zz_s(k),theta_S(k),qv_S(k)
        END Do
        !
        ! Compute 'first guess' pressure if given Z (ignoring moisture)
        !
        IF (snd_zz) Then
          PRINT *, 'Computing first guess pressure'
          k           = 1
          pres_S(k)   = psfc
          temp        = theta_S(k)
          IF (ctemp(1:1).EQ.'P') temp=theta_S(k)*(pres_S(k)*p00inv)**rcp
          tv          = temp

          Do k=2,n
            tvkm1     = tv
            dz        = zz_S(k) - zz_S(k-1)
            pres_S(k) = pres_S(k-1) * Exp (-grav*dz*rdi/tvkm1)
            !
            !  estimate Tv at k+1/2 to get more accurate integration of pres.
            !
            DO iter=1,2
              temp      = theta_S(k)
              IF (ctemp(1:1).EQ.'P')
     >          temp    = theta_S(k)*(pres_S(k)*p00inv)**rcp
              tv        = temp
              tvavg     = .5 * (tv + tvkm1)
              pres_S(k) = pres_S(k-1) * Exp (-grav*dz*rdi/tvavg)
            END DO

          END DO
        END IF
        !
        ! Compute theta if given temp
        !
        IF (ctemp(1:1).EQ.'T') Then
          PRINT *, 'Computing theta from ', ctemp
          Do k=1,n
            temp        = theta_S(k)
            theta_S(k)  = temp * (p00/pres_S(k)) ** rcp
          END Do
        END IF
        !
        ! Compute Qv if given RH or Td
        !
        IF (cmoist(1:1).EQ.'R') Then
          PRINT *, 'Computing mixing ratio from ', cmoist
          Do k=1,n
            rh      = qv_S(k)
            temp    = theta_S(k) * (pres_S(k)/p00) ** rcp
            qv_S(k) = rh * Fmixrat(pres_S(k),Fsvpres(temp-tfrz))
            qv_S(k) = MAX( 0.0, qv_S(k) )
          END Do
        ELSE IF (cmoist(1:1).EQ.'D') Then
          PRINT *, 'Computing mixing ratio from ', cmoist
          Do k=1,n
            tdew    = qv_S(k)
            temp    = theta_S(k) * (pres_S(k)/p00) ** rcp
            qv_S(k) = Fsvpres(tdew-tfrz)
            qv_S(k) = MAX( 0.0, qv_S(k) )
          END Do
        END IF

        IF (.NOT.snd_pres) Then
          CALL HydroPres (n, psfc, pres_S, zz_S, theta_S, qv_S)
        END IF

        stid = 'XXX'    ! We do not know
        istatus = 0

      ELSE IF (jformat <= 0) Then            ! still do not know the format
      !
      !  Non-spiffy sounding format
      !

        PRINT *, 'Non-spiffy sounding format'

        IF (sndfile.EQ.'arps.sound') Then    ! Special file name "arps.sound"
          PRINT *, 'Special format for file ', sndfile
          Read(iunit,*)
          DO j=1,100
            Read(iunit,*,END=9800) k, zz_S(k), pres_S(k), theta_S(k),
     &      xx,  qv_S(k), uu_S(k), vv_S(k)
            qv_S(k) = MAX( 0.0, qv_S(k) )
          END dO
 9800     CONTINUE
          n         = j - 1

          stnelev   = zz_S(1)
          psfc      = pres_S(1)
          plot_wind = .True.
          snd_pres= .True.
          snd_zz        = .True.
          !Do k=1,n
          !  zz_S(k)        = zz_S(k) - stnelev
          !END Do

          istatus = 0

          Return

        END IF
!
!-----------------------------------------------------------------------
!
!  Read the second line, check for GEMPAK format or FSL format
!  (First line does not contain either '%' or '&')
!
!-----------------------------------------------------------------------
!
      PRINT *,' Checking contents of second line'
      READ(iunit,'(a)') line2
      PRINT *, ' line2: ',line2

      IF (line2(:9).EQ.' SNPARM =') Then
        PRINT *, '% Sounding file is GEMPAK format.'
        ifCB  = -99
        pmbCB = 0.
        tcCB  = 0.
      ELSE
        PRINT *,' Checking for FSL format'
        read(line2,'(16x,i5,f7.2,1x,f6.2)',iostat=ios) wmo,slat,slon
        IF(ios .eq. 0) THEN
          PRINT *, ' wmo,slat,slon:',wmo,slat,slon
          IF( wmo .gt. 0 .and. 
     >        slat .gt. -90.1 .and. slat .lt. 90.1 .and.
     >        slon .gt. -180.1 .and. slon .lt. 180.1 ) THEN
            PRINT *, '% Sounding file is FSL format.'
            ifCB = -199 
            pmbCB = 0.
            tcCB = 0.
          END IF
        END IF 
      END IF
!
!-----------------------------------------------------------------------
!
!   Read the third line if not GEMPAK or FSL
!
!-----------------------------------------------------------------------

      IF (ifCB .GT. -99) READ(iunit,*) ifCB, pmbCB, tcCB

      IF (ifCB.GT.2) Then
      !
      !  Model format
      !
        snd_zz = .True.
        ifCB   = 0
        pmbCB  = 0.
        tcCB   = 0.
        Backspace (1)
        Read(iunit,*) n,dz,stnelev,psfc, ifCB, pmbCB, tcCB
        n                = n * 2
        dz        = dz * .5
        PRINT '(1x,a,i3,3f9.2)', 
     >    'Model format sounding: n,dz,stnelev,psfc:', n,dz,stnelev,psfc
        Read(iunit,*) (theta_S(k),k=1,n)
        Read(iunit,*) (qv_S(k),k=1,n)
        Read(iunit,*) (uu_S(k),k=1,n)
        Read(iunit,*) (vv_S(k),k=1,n)
        Close (1)
        Do k=1,n
          zz_S(k) = (k-1) * dz
          qv_S(k) = MAX( 0.0, qv_S(k) )
          IF (uu_S(k).NE.0. .OR. vv_S(k).NE.0.) plot_wind = .True.
        END Do

        CALL HydroPres (n, psfc, pres_S, zz_S, theta_S, qv_S)

      ELSE   ! may be GEMPAK foramt, FSL format, or new snd_wind
      !
      !  Raw format
      !
        snd_pres  = .True.
        snd_zz    = .False.
        plot_wind = .False.
        IF (snd_wind) plot_wind = .TRUE.
        gempak_format = .False.
        !
        ! Gempak format
        ! New GEMPAK 5.0 has extra line. Both fmts supported here (1995/02/15)
        !
        IF (ifCB.EQ.-99) Then
          PRINT *, ' Reading GEMPAK sounding'
          snd_zz        = .True.
          plot_wind        = .True.
          gempak_format = .True.
          Read(iunit,*)
!STID = OUN        STNM =    72357   TIME = 921119/ 0 0
!23456789012345678901234567890123456789012345678901234567890123456789012
          Read(iunit,9987) stid,stnm,iyr,imon,iday,ihr,imin
 9987     Format (8x,a3,17x,A6,10x,3i2,1x,2i2)
            Read(iunit,1972) slat,slon 
 1972     Format (11x,f5.2,11x,f7.2)
          Read(iunit,*)
          Read(iunit,'(A)') line
          IF (line(7:10).NE.'PRES') READ (1,*) line
          iyr = iyr + 1900
          IF (iyr .LT. 1950) iyr = iyr + 100
        ELSE if (ifCB .EQ. -199) THEN
          PRINT *, ' Reading FSL-formatted sounding'
          snd_zz        = .True.
          plot_wind        = .True.
          fsl_format = .True.
          Read(line,'(3i7,6x,a4,i7)') lintyp,ihr,iday,chmon,iyr
          Read(line2,'(3i7,f7.2,a1,f6.2,a1,i6,i7)')
     >     lintyp,wban,wmo,slat,ns,slon,ew,ielev,itime
          Read(iunit,'(7i7)') 
     >         lintyp,hydro,mxwd,tropl,lines,tindex,source
          Read(iunit,'(i7,11x,a3,14x,i7,5x,a2)') 
     >         lintyp,staid,sonde,wsunits
          stid=staid
          DO imon=1,11
             IF(chmon .eq. monlist(imon)) GO TO 91
          END Do
  91      CONTINUE
          PRINT *, ' chmon, imon =',chmon,imon

          IF(wsunits.eq.'kt') THEN
            spcvt=kts2ms
          ELSE
            spcvt=1.0
          END IF
          IF(ns .eq. 'S') slat=-slat
          IF(ew .eq. 'W') slon=-slon
        END IF
!
!
!   Read T(C),Td(C), convert to theta, qv.
!   For Gempak, also convert DDSS to U,V
!
        IF (gempak_format) THEN

          DO k=1,nmax
            READ(iunit,*,end=1000)
     >          pres_S(k),zz_S(k),tempc_S(k),tdewc_S(k),uu_S(k),vv_S(k)
            IF (uu_S(k) < 0. .OR. vv_S(k) <= 0) Then
              uu_S(k) = 0.
              vv_S(k) = 0.
            END IF
            pres_S(k) = pres_S(k) * 100.
            IF (tdewc_S(k).LT.-199.) tdewc_S(k) = -199.
          END DO

 1000     CONTINUE

        ELSE IF (fsl_format) Then
          k=0
          ipsfc=999999
          nlines=lines-4
          DO iline=1,nlines
            Read(iunit,'(7i7)',end=1001) lintyp,(idata(j),j=1,6)
            PRINT *, lintyp,(idata(j),j=1,6)
            IF(idata(2).ge.ielev .AND. idata(2).lt.99990 .AND.
     :        idata(1).lt.ipsfc ) THEN
              k=k+1
              IF(idata(1).lt.99990) THEN
                pres_S(k)=10.0*idata(1)
              ELSE
                pres_S(k)=-199.
              END IF
              IF(k .eq. 1) ipsfc=idata(1)
              zz_S(k)=float(idata(2))
              IF(idata(3).lt.99990) THEN
                tempc_S(k)=idata(3)*0.1
              ELSE
                tempc_S(k)=-199.
              END IF
              IF(idata(4).lt.9990) THEN
                tdewc_S(k)=idata(4)*0.1
              ELSE
                tdewc_S(k)=-199.
              END IF
              IF(idata(5).lt.99990 .and. idata(6).lt.200) THEN
                wdir=dtr*float(idata(5))
                wspd=spcvt*idata(6)
                uu_S(k)=-wspd*sin(wdir)
                vv_S(k)=-wspd*cos(wdir)
              ELSE
                uu_S(k)=-199.
                vv_S(k)=-199.
              END IF
            END IF
          END Do
 1001     CONTINUE
          nlevs=k
!
! Interpolate missing temperatures
!
          DO k=2,nlevs
            IF(tempc_S(k) .lt. -190.) THEN
              kbot=k-1
              foundt=.false.
              DO kn=k+1,nlevs
                IF(tempc_S(kn) .gt. -190.) THEN
                  foundt=.true.
                  EXIT
                END IF
              END DO
              IF(foundt) THEN
                ktop=kn
                dtdz=(tempc_S(ktop)-tempc_S(kbot))/     
     >               (zz_S(ktop)-zz_S(kbot))
                DO kk=k,ktop-1
                  tempc_S(kk)=tempc_S(kbot)+dtdz*(zz_S(kk)-zz_S(kbot))
                END DO
              END IF
            END IF
          END DO
!
! Interpolate missing dew points
!
          DO k=2,nlevs
            IF(tdewc_S(k) .lt. -190.) THEN
              kbot=k-1
              foundt=.false.
              DO kn=k+1,nlevs
                IF(tdewc_S(kn) .gt. -190.) THEN
                  foundt=.true.
                  EXIT
                END IF
              END DO
              IF(foundt) THEN
                ktop=kn
                dtddz=(tdewc_S(ktop)-tdewc_S(kbot))/
     >                (zz_S(ktop)-zz_S(kbot))
                DO kk=k,ktop-1
                  tdewc_S(kk)=tdewc_S(kbot)+dtddz*(zz_S(kk)-zz_S(kbot))
                END DO
              END IF
            END IF
          END DO
!
! Interpolate missing winds
!
          DO k=2,nlevs
            IF(uu_S(k) .lt. -190.) THEN
              kbot=k-1
              foundw=.false.
              DO kn=k+1,nlevs
                IF(uu_S(kn) .gt. -190.) THEN
                  foundw=.true.
                  EXIT
                END IF
              END DO
              IF(foundw) THEN
                ktop=kn
                dudz=(uu_S(ktop)-uu_S(kbot)) / (zz_S(ktop)-zz_S(kbot))
                dvdz=(vv_S(ktop)-vv_S(kbot)) / (zz_S(ktop)-zz_S(kbot))
                DO kk=k,ktop-1
                  uu_S(kk)=uu_S(kbot)+dudz*(zz_S(kk)-zz_S(kbot))
                  vv_S(kk)=vv_S(kbot)+dvdz*(zz_S(kk)-zz_S(kbot))
                END DO
              END IF
            END IF
          END DO
      
        ELSE IF (snd_wind) Then  ! neither GEMPAK format nor FSL format
          DO k=1,nmax
            READ(iunit,*,END=1002) 
     >           pres_S(k),tempc_S(k),tdewc_S(k),uu_S(k),vv_S(k)
            pres_S(k) = pres_S(k) * 100.
            IF (tdewc_S(k).LT.-199.) tdewc_S(k) = -199.

          END DO
 1002     CONTINUE
        ELSE
          DO k=1,nmax
            Read(iunit,*,END=1003) pres_S(k), tempc_S(k), tdewc_S(k)
            pres_S(k) = pres_S(k) * 100.
            IF (tdewc_S(k).LT.-199.) tdewc_S(k) = -199.
          END DO
 1003     CONTINUE
        END IF

        CLOSE (1)
        n = k - 1
        IF (gempak_format) Then
          stnelev  = zz_S(1)
          psfc     = pres_S(1)
          !Do k=1,n
          !  zz_S(k) = zz_S(k) - stnelev
          !END Do
          CALL DDSS2UV (n, uu_S, vv_S)
        END IF
!
        CALL ConvertSound (n, pres_S, zz_S, theta_S, qv_S, 
     >    tempc_S, tdewc_S, psfc, snd_pres, snd_zz)
!
!        CALL HydroPres (n, psfc, pres_S, zz_S, theta_S, qv_S) !*** test ***
!
!9977 Format (i3,f8.2,f7.0,3f7.2)
!       Do k=1,n
!          PRINT 9977, k,pres_S(k)/100., zz_S(k), theta_S(k), qv_S(k)*1000.
!        end do
      END IF    ! 

      ELSE IF (jformat.NE.-99) Then   ! ifCB > 2
      !
      !     Spiffy New Format
      !
      !
      !...Format 01: P, T, Td
      !
      IF (jformat.EQ.01) Then
        Read(iunit,*)
        Read(iunit,*)
        Read(iunit,*)
        Read(iunit,*) n, dz, stnelev, psfc
        IF (ifCB.GT.0) Read(iunit,*) ifCB, pmbCB, tcCB
        Do k=1,nmax
          IF (snd_wind) Then
            Read(iunit,*,END=1111) zz_S(k), tempc_S(k), tdewc_S(k),
     >        uu_S(k), vv_S(k)
          ELSE 
            Read(iunit,*,END=1111) zz_S(k), tempc_S(k), tdewc_S(k)
          END IF
        END Do
 1111   Continue
        n        = k - 1
        !Do k=1,n
        !  zz_S(k) = zz_S(k) - stnelev
        !END Do
!
!
!
!...Format 02: Z, Theta, Qv, U, V
!
      ELSE IF (jformat.LE.02) Then
        Read(iunit,*)
        Read(iunit,*)
        Read(iunit,*)
        Read(iunit,*) n, dz, stnelev, psfc
        IF (ifCB.GT.0) Read(iunit,*) ifCB, pmbCB, tcCB
        k1        = 1
        k2        = n
        kinc        = 1
        IF (minus_z) Then
          k1        = n
          k2        = 1
          kinc        = -1
        END IF
        Do k=k1,k2,kinc
          IF (snd_wind) Then
            Read(iunit,*) zz_S(k),theta_S(k),qv_S(k),uu_S(k),vv_S(k)
          ELSE 
            Read(iunit,*) zz_S(k),theta_S(k),qv_S(k)
          END IF
          qv_S(k) = MAX( 0.0, qv_S(k) )
          PRINT *, k, zz_s(k), theta_S(k), qv_S(k)
        END Do
        !Do k=1,n
        !  zz_S(k)        = zz_S(k) - stnelev
        !END Do

        CALL HydroPres (n, psfc, pres_S, zz_S, theta_S, qv_S)

      ELSE
        !
        !...Unrecognized Format 
        !
        PRINT *, 'Unrecognized sounding format: ', jformat
        STOP
      END IF  ! jformat

      END IF  ! ifCB
!----------------------------------------------------------------------=
!
!  Done reading the sounding file
!
!-----------------------------------------------------------------------

      IF (snd_wind) plot_wind = .TRUE.

      PRINT '(i5,2a)', n ,' points read in from ', TRIM(sndfile)
      istatus = 0

      END

!
!########################################
!########################################
!########                        ########
!########        DDSS2UV         ########
!########                        ########
!########################################
!########################################
!
      Subroutine DDSS2UV (n, uu, vv)
      Implicit None
      Integer        k,n
      Real        uu(n), vv(n), dir, spd, deg2rad
      Parameter (deg2rad=.0174532925)
!
!  On input, uu is dir, vv is spd
!  On output, uu, vv are Cartesian speeds
!
      Do k=1,n
        dir        = uu(k)
        spd        = vv(k)
        uu(k)        = - spd * Sin(deg2rad*dir)
        vv(k)        = - spd * Cos(deg2rad*dir)
      END Do
!
      END
!
!
!
!
!                   ########################################
!                   ########################################
!                   ########################################
!                   ########                        ########
!                   ########      CONVERTSOUND      ########
!                   ########                        ########
!                   ########################################
!                   ########################################
!                   ########################################
!
!
      Subroutine ConvertSound (n, pres, zz, theta, qv, tempc, tdewc, 
     >    psfc, snd_pres, snd_zz)
      Implicit None
      Include 'thermo.consts'
      Integer        k, n
      Real        pres(n), zz(n), theta(n), qv(n), tempc(n), tdewc(n)
      Real        psfc, dz, tvavg, tvkp1, tvkm1, tv, temp
      Logical        snd_pres, snd_zz
      Include 'thermo.stfunc'
!
!
!  Convert sounding containing (p|z,T,Td) to a standard set of variables
!  (p,z,theta,qv).
!
!
!...Has Z but not P -- compute hydrostatically
!
      IF (.NOT.snd_pres) Then
      PRINT *, '% ConvertSound: computing (p,theta,qv) from (z,T,Td)'
!
      k                = 1
      pres(k)        = psfc
      temp        = tempc(k) + tfrz
      theta(k)        = Ftheta(pres(k),temp)
      qv(k)        = Fmixrat(pres(k),Fsvpres(tdewc(k)))
      tv        = Ftvirtnowl(temp,qv(k))
!
      Do k=1,n-1
        dz        = zz(k+1) - zz(k)
        pres(k+1) = pres(k) * Exp (-grav*dz*rdi/tv)
!
!  estimate Tv at k+1 to get more accurate integration of pres.
!
        temp        = tempc(k+1) + tfrz
        theta(k+1) = Ftheta(pres(k+1),temp)
        qv(k+1)        = Fmixrat(pres(k+1),Fsvpres(tdewc(k+1)))
        tvkp1        = Ftvirtnowl(temp,qv(k+1))
        tvavg        = .5 * (tv + tvkp1)
        pres(k+1) = pres(k) * Exp (-grav*dz*rdi/tvavg)
!
!  recompute Th,Qv
!
        theta(k+1) = Ftheta(pres(k+1),temp)
        qv(k+1)        = Fmixrat(pres(k+1),Fsvpres(tdewc(k+1)))
      END Do
!
!
!...Has pressure and heights
!
      ELSE IF (snd_zz) Then
      PRINT *, '% ConvertSound: computing (Th,Qv) from (p,z,T,Td)'
!
      Do k=1,n
        theta(k)= Ftheta(pres(k),tempc(k)+tfrz)
        qv(k)        = Fmixrat(pres(k),Fsvpres(tdewc(k)))
      END Do
!
!
!...Has pressure but not heights
!
      ELSE IF (.NOT.snd_zz) Then
      PRINT *, '% ConvertSound: computing (z,Th,Qv) from (p,T,Td)'
!
      Do k=1,n
        theta(k)= Ftheta(pres(k),tempc(k)+tfrz)
        qv(k)        = Fmixrat(pres(k),Fsvpres(tdewc(k)))
      END Do
!
      zz(1)        = 0.
      tvkm1        = Ftvirtnowl(tempc(1)+tfrz,qv(1))
!
      Do k=2,n
        tv        = Ftvirtnowl(tempc(k)+tfrz,qv(k))
        tvavg        = .5 * (tvkm1 + tv)
        dz        = rd/grav*tvavg * Log(pres(k-1)/pres(k))
        zz(k)        = zz(k-1) + dz
        tvkm1        = tv
      END Do
!
!
      END IF
!
!
      END
!
!
!
!
!                   ########################################
!                   ########################################
!                   ########################################
!                   ########                        ########
!                   ########        HYDROPRES       ########
!                   ########                        ########
!                   ########################################
!                   ########################################
!                   ########################################
!
!
      Subroutine HydroPres (n, psfc, pres, zz, theta, qv)
      Implicit None
      Include 'thermo.consts'
      Integer        k,n, iter
      Real        pres(n), zz(n), theta(n), qv(n)
      Real        dz, tv, tvkm1, tvavg, temp, psfc
      Include 'thermo.stfunc'
!
!
!
      PRINT *, '% HydroPres: compute p from (z,Th,Qv)'
!
      k                = 1
      pres(1)        = psfc
      temp        = theta(k) * (pres(k) * p00inv) ** rcp
      tv        = Ftvirtnowl(temp,qv(k))
!
      Do k=2,n
        tvkm1        = tv
        dz        = zz(k) - zz(k-1)
        pres(k) = pres(k-1) * Exp (-grav*dz*rdi/tvkm1)
!        PRINT '(2i3,f7.2)', 0, k, pres(k)/100.
!
!  estimate Tv at k+1/2 to get more accurate integration of pres.
!
        Do iter=1,2
          temp        = theta(k) * (pres(k) * p00inv) ** rcp
          tv        = Ftvirtnowl(temp,qv(k))
          tvavg        = .5 * (tv + tvkm1)
          pres(k) = pres(k-1) * Exp (-grav*dz*rdi/tvavg)
!          PRINT '(2i3,f7.2)', iter, k, pres(k)/100.
        END Do
!
      END Do

!     do k=1,n
!     tempc        = theta(k) * (pres(k) * p00inv) ** rcp - tfrz
!     tdewc        = Ftdewc(Fvpres(pres(k),qv(k)))
!     PRINT 9999, k, zz(k), pres(k)/100., tempc, tdewc,
!    >  theta(k), qv(k)*1000.
!     end do
!9999 Format (i4,f7.0,f8.2, 4f7.2)
!
      END
!
!#######################################################################
!#######################################################################
!########                                                       ########
!########      WRITEARPSSND                                     ########
!########                                                       ########
!#######################################################################
!#######################################################################
!
      SUBROUTINE WriteARPSSnd(filename, zrefsfc, psfc, n,               
     &                        pres, z, theta, temp, qv, rh, tdew, u, v)
      IMPLICIT NONE

!     RLC 1997/04/22
!     All units MKS

      INTEGER n
      REAL zrefsfc, psfc
      REAL pres(n), z(n), theta(n), temp(n), qv(n), rh(n), tdew(n),
     &     u(n), v(n)
      CHARACTER*(*) filename

      INTEGER k
      CHARACTER*16 hgtstr, tempstr, mststr, windstr
      DATA hgtstr /"'height'"/,        tempstr /"'temp.'"/,
     &     mststr /"'rel. humidity'"/, windstr /"'uv'"/

!----------------------------------------------------------------------=

      PRINT *, '% WriteARPSSnd: ', filename
      OPEN (UNIT=1, FILE=filename)

      WRITE (1,9000) 'Line 1'
      WRITE (1,9000) 'Line 2'
      WRITE (1,9000) 'Line 3'
      WRITE (1,9000) 'Line 4'
      WRITE (1,9000) 'Line 5'
      WRITE (1,'(12(A," "))') hgtstr, tempstr, mststr, windstr
      WRITE (1,'(2F12.1)')    zrefsfc, psfc
      WRITE (1,*) n
      WRITE (1,*) '-----------------------------------------'

      DO k=n,1,-1
        WRITE (1,9010) z(k), temp(k), rh(k), u(k), v(k)
      END DO

      CLOSE (1)

 9000 FORMAT (16A)
 9010 FORMAT (2F12.2, F10.4, 2F12.3)
      END
