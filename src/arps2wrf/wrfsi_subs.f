c
cdis    Forecast Systems Laboratory
cdis    NOAA/OAR/ERL/FSL
cdis    325 Broadway
cdis    Boulder, CO     80303
cdis
cdis    Forecast Research Division
cdis    Local Analysis and Prediction Branch
cdis    LAPS   
cdis
cdis    This software and its documentation are in the public domain and
cdis    are furnished "as is."  The United States government, its
cdis    instrumentalities, officers, employees, and agents make no
cdis    warranty, express or implied, as to the usefulness of the software
cdis    and documentation for any purpose.  They assume no responsibility
cdis    (1) for the use of the software and documentation; or (2) to provide
cdis     technical support to users.
cdis
cdis    Permission to use, copy, modify, and distribute this software is
cdis    hereby granted, provided that the entire disclaimer notice appears
cdis    in all copies.  All modifications to this software must be clearly
cdis    documented, and are solely the responsibility of the agent making
cdis    the modifications.  If significant modifications or enhancements
cdis    are made to this software, the FSL Software Policy Manager
cdis    (softwaremgr@fsl.noaa.gov) should be notified.
cdis

C
C     ******************************************************************
C
      SUBROUTINE SFCOPQR(NO,MOF,NP,NIQ,NJQ,N2,N3,lcat
     +          ,XT,YT,RLAT,WLON1,ERAD,RWOFF,RSOFF
     +          ,DELTALLO,DELTAXP,DELTAYP,DELTAXQ,DELTAYQ
     +          ,IBLKSIZO,ISBEGO,IWBEGO,DATR,DATS,DATLN,DATLT
     +          ,OFN,WVLN,SILWT,dem_data,maxdatacat,istat_files)

c JS: removed dato array from subroutine argument list
c JS: added RWOFF/RSOFF - West and South offset of tile data

      real,  allocatable ::  dato(:,:,:,:)    !dato(no,no,mof,lcat)

      real,  allocatable ::  DATP(:,:,:)
      real,  allocatable ::  DATQ(:,:)
      real,  allocatable ::  DATQS(:,:,:)
      real,  allocatable ::  DATSM(:,:)
      real,  allocatable ::  DATSMX(:,:)
      real,  allocatable ::  DATSLN(:,:)
      real,  allocatable ::  DATSLT(:,:)

      real,  intent(inout) ::  DATLN(N2,N3)
      real,  intent(out)   ::  DATLT(N2,N3)
      real,  intent(out)   ::  DATR(N2,N3)
      real,  intent(out)   ::  DATS(N2,N3,maxdatacat)

      real ISO(MOF),IWO(MOF),XT(N2),YT(N3),rlat,wlon1,
     +     erad,deltallo,deltaxp,deltayp,deltaxq,deltayq,
     +     wvln,silwt,xq,yq,xp,yp,xcentr,ycentr,glatp,               ! pla,plo,
     +     glonp,rio,rjo,wio2,wio1,wjo2,wjo1,xq1,yq1,
     +     rwoff,rsoff

      real r_missing_data
      real xr,yr,rval,sh,sha,rh,rha,rhn,rht,shn,sht
      real shln,shlt,rhln,rhlt
      real delta_ln(np,np),delta_lt(np,np)

c     real xpmn,xpmx,ypmn,ypmx
      real xp1,xp2,yp1,yp2
      real xpcentr,ypcentr

      real  pctcat(maxdatacat)

      integer lp
      integer ixr,iyr
      integer lent

      CHARACTER*180 OFN
      CHARACTER*180 TITLE3,TITLE3_last_read,TITLE3_last_inquire
      CHARACTER*3   TITLE1
      CHARACTER*4   TITLE2
      CHARACTER*10  cdatatype

      LOGICAL L1,L2,dem_data,l_string_contains

      data icnt/0/
      save icnt
C
      print *,'no,mof,np,niq,njq=',no,mof,np,niq,njq

      istat_files = 1

      NONO=NO*NO
      XCENTR=0.5*(XT(1)+XT(N2))
      YCENTR=0.5*(YT(1)+YT(N3))
      print *,xt(1),xt(n2),xcentr
      print *,'deltaxp=',deltaxp
      NOFR=0
      DO 11 IOF=1,MOF
         ISO(IOF)=0
         IWO(IOF)=0
  11  continue

      TITLE3_last_read    = '/dev/null'
      TITLE3_last_inquire = '/dev/null'

      lcat=1
      len=index(ofn,' ')
      if(ofn(len-1:len-1).eq.'V'.or.
     &   ofn(len-1:len-1).eq.'I')then
         icnt = 0
         cdatatype='landuse'
         if(ofn(len-1:len-1).eq.'I')cdatatype='islope'
      elseif(ofn(len-1:len-1).eq.'O')then
         icnt = 0
         cdatatype='soiltype'
      elseif(ofn(len-1:len-1).eq.'U' .or.
     &       ofn(len-1:len-1).eq.'H')then
         cdatatype='topography'
      endif

      print*,'SFCOPQR: cdatatype = ',cdatatype

!      call s_len(cdatatype,lent)
      lent = LEN_TRIM(cdatatype)

      allocate(dato(no,no,mof,lcat))

      allocate (DATP(NP,NP,lcat),
     &          DATQ(NIQ,NJQ),
     &          DATSM(NIQ,NJQ),
     &          DATSMX(NIQ,NJQ),
     &          DATSLN(NIQ,NJQ),
     &          DATSLT(NIQ,NJQ),
     &          DATQS(NIQ,NJQ,maxdatacat))

c      call get_r_missing_data(r_missing_data,istatus)
      r_missing_data = 1e+37
      istatus = 1

      if(istatus.ne.1)then
         print*,'failed to get r_missing_data'
         return
      endif

      DO 15 JQ=1,NJQ
         print *,'jq,njq,niq,nofr=',jq,njq,niq,nofr
         DO 16 IQ=1,NIQ

            XQ=(FLOAT(IQ)-0.5*FLOAT(NIQ+1))*DELTAXQ+XCENTR
            YQ=(FLOAT(JQ)-0.5*FLOAT(NJQ+1))*DELTAYQ+YCENTR

            xpmn=1.0e30
            ypmn=1.0e30
c           xpmx=-1.0e30
c           ypmx=-1.0e30

            DO 17 JP=1,NP
               DO 18 IP=1,NP

                  XP=XQ+(FLOAT(IP)-0.5*FLOAT(NP+1))*DELTAXP
                  YP=YQ+(FLOAT(JP)-0.5*FLOAT(NP+1))*DELTAYP

!                  call xy_to_latlon(XP,YP,erad,GLATP,GLONP) 
                  CALL xytoll(1,1,xp,yp,GLATP,GLONP)
!         print *, 'GLATP = ',glatp,' GLONP = ',glonp

                  glatp = max(-89.9999,min(89.9999,glatp - rsoff))
                  glonp = glonp - rwoff
                  if(glonp.ge.180.) glonp = glonp - 360.
                  if(glonp.le.-180.) glonp = glonp + 360.

c                 print *,'rlat,wlon1=',rlat,wlon1

                  ISOC=(INT((GLATP-FLOAT(ISBEGO))/FLOAT(IBLKSIZO)
     &          +200.)-200)*IBLKSIZO+ISBEGO
            IWOC=(INT((GLONP-FLOAT(IWBEGO))/FLOAT(IBLKSIZO)
     &          +400.)-400)*IBLKSIZO+IWBEGO

                  DO 19 IOFR=1,NOFR
                     JOFR=IOFR
                     IF(ISO(IOFR).EQ.ISOC.AND.IWO(IOFR).EQ.IWOC)GO TO 10
 19                 continue
                  ISOCPT=ABS(ISOC)/10
                  ISOCPO=ABS(ISOC)-ISOCPT*10
                  IWOCPH=ABS(IWOC)/100
                  IWOCPT=(ABS(IWOC)-IWOCPH*100)/10
                  IWOCPO=ABS(IWOC)-IWOCPH*100-IWOCPT*10
                  IF(ISOC.GE.0)THEN
                     WRITE(TITLE1,'(2I1,A1)')ISOCPT,ISOCPO,'N'
                  ELSE
                     WRITE(TITLE1,'(2I1,A1)')ISOCPT,ISOCPO,'S'
                  ENDIF

                  IF(IWOC.GE.0 
     1               .and. IWOC .ne. 180                    ! 1998 Steve Albers
     1                                      )THEN
                     WRITE(TITLE2,'(3I1,A1)')IWOCPH,IWOCPT,IWOCPO,'E'
                  ELSE
                     WRITE(TITLE2,'(3I1,A1)')IWOCPH,IWOCPT,IWOCPO,'W'
                  ENDIF

                  LB=LEN_TRIM(ofn)
                  TITLE3=OFN(1:LB)//TITLE1//TITLE2
                  LB=INDEX(TITLE3,' ')-1

!        print *, 'Title3 = ',title3
!        print *, 'Title3_last_inquire = ',title3_last_inquire
                  if(TITLE3 .ne. TITLE3_last_inquire)then
                     INQUIRE(FILE=TITLE3(1:LB),EXIST=L1,OPENED=L2)
                     TITLE3_last_inquire = TITLE3
                  endif

                  IF(.NOT.L1)THEN
                     iwrite = 0

                     if(icnt .le. 100)then ! Reduce the output
                         iwrite=1

                     elseif(icnt .le. 1000)then
                         if(icnt .eq. (icnt/100)*100)iwrite=1

                     elseif(icnt .le. 10000)then
                         if(icnt .eq. (icnt/1000)*1000)iwrite=1

                     elseif(icnt .le. 100000)then
                         if(icnt .eq. (icnt/10000)*10000)iwrite=1

                     else
                         if(icnt .eq. (icnt/100000)*100000)iwrite=1

                     endif

                     if(iwrite .eq. 1)then
                        if(l_string_contains(TITLE3(1:LB),
     1                                       'world_topo_30s',
     1                                       istatus)             )then       
                           PRINT*, ' ERROR: ',TITLE3(1:LB)
     1                            ,' DOES NOT EXIST ',icnt

                        elseif(l_string_contains(TITLE3(1:LB),
     1                                           'topo_30s',
     1                                       istatus)             )then       
                             PRINT*, ' topo_30s file ',TITLE3(1:LB)
     1                            ,' does not exist, using topo_10m '
     1                            ,icnt

                        else ! Generic warning message
                           PRINT*, ' WARNING: ',TITLE3(1:LB)
     1                            ,' DOES NOT EXIST ',icnt

                        endif

                     endif ! iwrite

                     icnt = icnt + 1

c initialize these arrays as they may have some garbage in them
c if we dont actually read in any data.
c
                     DATP(IP,JP,:) = 0.
                     DELTA_LN(IP,JP) = 0.
                     DELTA_LT(IP,JP) = 0.
                     istat_files = 0
                     GO TO 20

                  ENDIF

                  IF(NOFR.GE.MOF)THEN
                     DO 21 IOF=1,MOF
                        ISO(IOF)=0
                        IWO(IOF)=0
21                    continue
                     NOFR=0
                  ENDIF
                  NOFR=NOFR+1
                  JOFR=NOFR

!                 Read the tile
                  if(TITLE3 .ne. TITLE3_last_read)then
                    if( (ofn(len-1:len).eq.'U').and.(no.eq.1200).or.
     .                   no.eq.1201 )then
                         if(no.eq.1201)then
                            print*,'Reading ', title3(1:lb)
                            CALL READ_DEM(29,TITLE3(1:LB),no,no,4,4, ! topo_3s experimental
     .                              DATO(1,1,NOFR,1),istat)
                         else
                            CALL READ_DEM(29,TITLE3(1:LB),no,no,2,2, ! world topo_30s
     .                              DATO(1,1,NOFR,1),istat)
                         endif
                      dem_data=.true.
                    elseif( (ofn(len-1:len).eq.'O') )then      ! soiltype top and bot layer
                      CALL READ_DEM(29,TITLE3(1:LB),no,no,1,4,
     .                              DATO(1,1,NOFR,1),istat)
                      dem_data=.true.
                    elseif( (ofn(len-1:len).eq.'V') )then      ! world USGS 30s landuse
                      CALL READ_DEM(29,TITLE3(1:LB),no,no,1,4,
     .                              DATO(1,1,NOFR,1),istat)
                      dem_data=.true.
                    elseif((ofn(len-1:len).eq.'I'))then      ! only islope in this code section
                      CALL READ_DEM_G(29,TITLE3(1:LB),no,no,1,lcat
     .                     ,nofr, 1,4, DATO,istat)
                      dem_data=.true.
                    elseif( (ofn(len-1:len).eq.'T') )then      ! soiltemp - obsolete in this code
                      CALL READ_DEM(29,TITLE3(1:LB),no,no,2,2,
     .                              DATO(1,1,NOFR,1),istat)
                      dem_data=.true.
                    else                                       ! other
!                      CALL JCLGET(29,TITLE3(1:LB),'FORMATTED',0,istatus)      
             OPEN(29,STATUS='OLD',FILE=TITLE3(1:LB),FORM='FORMATTED')
                      CALL VFIREC(29,DATO(1,1,NOFR,1),NONO,'LIN')
                      if ((ofn(len-1:len).eq.'U').and.(no.eq.121)) then
                        dem_data=.false.                       ! topo_30s
                      endif
                    endif

                    if(istat.ne.0)then
                       print*,'Error returned: SFCOPQR: READ_DEM'
                       return
                    endif


                    TITLE3_last_read = TITLE3

c                   print *,'nofr,dato=',nofr,dato(1,1,nofr)
                    CLOSE(29)

                  else
                    write(6,*)' We have made the code more efficient'

                  endif ! Is this a new file we haven't read yet?

                  ISO(NOFR)=ISOC
                  IWO(NOFR)=IWOC
10		  continue

                  RIO=(GLONP-FLOAT(IWOC))/DELTALLO+1.
                  RJO=(GLATP-FLOAT(ISOC))/DELTALLO+1.

!                 Prevent Bounds Error (Steve Albers)
                  if(RIO .lt. 1.0)then
                      if(RIO .gt. 0.98)then
                          write(6,*)' Reset RIO for Machine Epsilon'      
                          RIO = 1.0
                      elseif(RIO .lt. 0.5)then
                          write(6,*)' ERROR: RIO out of bounds',RIO
                          stop
                      endif
                  endif

                  if(RJO .lt. 1.0)then
                      if(RJO .gt. 0.98)then
                          write(6,*)' Reset RJO for Machine Epsilon'      
                          write(6,*)JQ,IQ,
     1                          IP,JP,IO1,JO1,JOFR,RIO,RJO,GLATP,ISOC
                          RJO = 1.0
                      elseif(RJO .lt. 0.5)then
                          write(6,*)' ERROR: RJO out of bounds',RJO
                          write(6,*)JQ,IQ,
     1                          IP,JP,IO1,JO1,JOFR,RIO,RJO,GLATP,ISOC
                          stop
                      endif
                  endif

C Interp OK for continuous data such as topography

                  if(cdatatype.eq.'topography')then

                   IO1=INT(RIO)
                   JO1=INT(RJO)
                   IO2=IO1+1
                   JO2=JO1+1
                   WIO2=RIO-FLOAT(IO1)
                   WJO2=RJO-FLOAT(JO1)
                   WIO1=1.0-WIO2
                   WJO1=1.0-WJO2

                   do LP = 1,lcat

                   DATP(IP,JP,LP)=WIO1*(WJO1*DATO(IO1,JO1,JOFR,LP)
     +                                 +WJO2*DATO(IO1,JO2,JOFR,LP))
     +                           +WIO2*(WJO1*DATO(IO2,JO1,JOFR,LP)
     +                                 +WJO2*DATO(IO2,JO2,JOFR,LP))

!S & W-facing slopes > 0.
                   DELTA_LN(IP,JP)=
     .           ((DATO(IO2,JO1,JOFR,LP)-DATO(IO1,JO1,JOFR,LP))+
     .            (DATO(IO2,JO2,JOFR,LP)-DATO(IO1,JO2,JOFR,LP)))*.5

                   DELTA_LT(IP,JP)=
     .           ((DATO(IO1,JO2,JOFR,LP)-DATO(IO1,JO1,JOFR,LP))+
     .            (DATO(IO2,JO2,JOFR,LP)-DATO(IO2,JO1,JOFR,LP)))*.5

                   enddo !LP = 1,lcat

                  else

C Nearest grid point for landuse and soiltype

                   IO1=NINT(RIO)
                   JO1=NINT(RJO)
                   do LP = 1,lcat
                    DATP(IP,JP,LP)= DATO(IO1,JO1,JOFR,LP)
                   enddo

                  endif ! cdatatype eq topography.
                   
20               CONTINUE
18             continue ! IP
17           continue ! JP


!           print*,'xpmx/xpmn//ypmx/ypmn/ ',xpmx,xpmn,ypmx,ypmn

! Calculate average and silhouette terrain, then apply SILWT weight

            if(cdatatype(1:lent).eq.'topography')then

             SHA=0.
             RHA=0.
             RHLN=0.
             RHLT=0.
             shmax=0.

             DO 22 JP=1,NP
               SH=0.
               RH=0.
               RHN=0.
               RHT=0.
               DO 23 IP=1,NP
!                 Test for missing - then go to 16?
                  SH=max(SH,DATP(IP,JP,1)) 
                  RH=RH+DATP(IP,JP,1)
                  RHN=RHN+DELTA_LN(IP,JP)
                  RHT=RHT+DELTA_LT(IP,JP)
23             continue ! IP
               SHA=SHA+SH/(2.*FLOAT(NP))
               RHA=RHA+RH
               RHLN=RHLN+RHN
               RHLT=RHLT+RHT
               SHMAX=max(SHMAX,SH)
22           continue ! JP
 
             RHA=RHA/FLOAT(NP*NP)
             RMS=0.0
             DO 24 IP=1,NP 
               SH=0.
               DO 25 JP=1,NP
                  SH=max(SH,DATP(IP,JP,1))
                  RMS=RMS+((DATP(IP,JP,1)-RHA)*(DATP(IP,JP,1)-RHA))
25             continue ! JP
               SHA=SHA+SH/(2.*FLOAT(NP))
24           continue ! IP

             DATQS(IQ,JQ,1)=SQRT(RMS/FLOAT(NP*NP))
             DATQ(IQ,JQ)=SHA*SILWT+RHA*(1.-SILWT)
             DATSM(IQ,JQ)=RHA                           !mean value of points used for IQ,JQ
             DATSMX(IQ,JQ)=SHMAX                        !max value from points used for IQ,JQ
             DATSLN(IQ,JQ)=RHLN/FLOAT(NP*NP)/DELTAXP
             DATSLT(IQ,JQ)=RHLT/FLOAT(NP*NP)/DELTAYP

c            print *,'datq=',datq(iq,jq)

            elseif(cdatatype(1:lent).eq.'islope'    .or.
     &             cdatatype(1:lent).eq.'landuse'   .or.
     &             cdatatype(1:lent).eq.'soiltype'    )then

             call compute_categories(cdatatype,np*np,DATP(1,1,1)
     &               ,maxdatacat,domcat,pctcat)
             datq(iq,jq)=domcat 
             datqs(iq,jq,:)=pctcat(:)

            elseif(cdatatype(1:lent).eq.'greenfrac'   )then

c dominant greenness fraction for each month

             do lp=1,lcat
              call compute_categories(cdatatype,np*np,DATP(1,1,lp)
     &                               ,1,domcat,pctcat)
              datqs(iq,jq,lp)=domcat
             enddo

            endif

16       continue ! IQ
15    continue ! JQ

      print *,'after 15'
 
      XQ1=(1.-0.5*FLOAT(NIQ+1))*DELTAXQ+XCENTR
      YQ1=(1.-0.5*FLOAT(NJQ+1))*DELTAYQ+YCENTR

      if(cdatatype(1:lent).eq.'topography')then

        print*
        print*,'Before GDTOST2'
        print*,'--------------'
        print*,'datq(1,1)/(niq,njq)= ',datq(1,1),datq(niq,njq)
        print*,'datqs(1,1,1)/(niq,njq)= ',datqs(1,1,1),datqs(niq,njq,1)
        print*,'datsln(1,1)/(niq,njq)= ',datsln(1,1),datsln(niq,njq)
        print*,'datslt(1,1)/(niq,njq)= ',datslt(1,1),datslt(niq,njq)
        print*,'Mean/Max topo at IQ,JQ (1,1)/(niq,njq): '
     +,datsm(1,1),datsmx(1,1),datsm(niq,njq),datsmx(niq,njq)

        DO 28 JR=1,N3
         DO 29 IR=1,N2
           XR=(XT(IR)-XQ1)/DELTAXQ+1.
           YR=(YT(JR)-YQ1)/DELTAYQ+1.

           CALL GDTOST2(DATQ,NIQ,NJQ,XR,YR,RVAL)
           DATR(IR,JR)=max(0.,RVAL)
           if( DATR(IR,JR).gt.30000. )then
               print*,'Warning: value out of bounds'
           endif    

           CALL GDTOST2(DATQS,NIQ,NJQ,XR,YR,RVAL)
           DATS(IR,JR,1)=max(0.,RVAL)
           CALL GDTOST2(DATSLN,NIQ,NJQ,XR,YR,RVAL)
           DATLN(IR,JR)=RVAL
           CALL GDTOST2(DATSLT,NIQ,NJQ,XR,YR,RVAL)
           DATLT(IR,JR)=RVAL

 29      CONTINUE
 28     CONTINUE

        print*,'After GDTOST2'
        print*,'-------------'
        print*,'datr(1,1)/(n2,n3)= ',datr(1,1),datr(N2,N3)
        print*,'dats(1,1,1)/(n2,n3)= ',dats(1,1,1),dats(n2,n3,1)
        print*,'datln(1,1)/(n2,n3)= ',datln(1,1),datln(n2,n3)
        print*,'datlt(1,1)/(n2,n3)= ',datlt(1,1),datlt(n2,n3)
 
      elseif(cdatatype(1:lent).eq.'landuse'.or.
     +       cdatatype(1:lent).eq.'islope' .or.
     +       cdatatype(1:lent).eq.'soiltype')then

        DO 38 JR=1,N3
         DO 39 IR=1,N2
            IXR=NINT((XT(IR)-XQ1)/DELTAXQ)+1.
            IYR=NINT((YT(JR)-YQ1)/DELTAYQ)+1.
            if(ixr.lt.1)ixr=1
            if(iyr.lt.1)iyr=1
            if(ixr.gt.n2)ixr=niq
            if(iyr.gt.n3)iyr=njq

            datr(ir,jr)=  datq(ixr,iyr)     !dominant category
            dats(ir,jr,:)=datqs(ixr,iyr,:)  !percent dist for ea category

 39      CONTINUE
 38     CONTINUE

      endif

      deallocate(dato)
      deallocate(DATP,
     &           DATQ,
     &           DATQS,
     &           DATSM, 
     &           DATSMX,
     &           DATSLN,
     &           DATSLT)

      RETURN
      END

      SUBROUTINE GDTOST2(A,IX,IY,STAX,STAY,STAVAL)
*  SUBROUTINE TO RETURN STATIONS BACK-INTERPOLATED VALUES(STAVAL)
*  FROM UNIFORM GRID POINTS USING OVERLAPPING-QUADRATICS.
*  GRIDDED VALUES OF INPUT ARRAY A DIMENSIONED A(IX,IY),WHERE
*  IX=GRID POINTS IN X, IY = GRID POINTS IN Y .  STATION
*  LOCATION GIVEN IN TERMS OF GRID RELATIVE STATION X (STAX)
*  AND STATION COLUMN.
*  VALUES GREATER THAN 1.0E30 INDICATE MISSING DATA.
*
      real A(IX,IY),R(4),SCR(4),stax,stay,staval
     +  ,fixm2,fiym2,yy,xx
      IY1=INT(STAY)-1
      IY2=IY1+3
      IX1=INT(STAX)-1
      IX2=IX1+3
      STAVAL=1E30
      FIYM2=FLOAT(IY1)-1
      FIXM2=FLOAT(IX1)-1
      II=0
      DO 100 I=IX1,IX2
      II=II+1
      IF(I.GE.1.AND.I.LE.IX) GO TO 101
      SCR(II)=1E30
      GO TO 100
101   JJ=0
      DO 111 J=IY1,IY2
      JJ=JJ+1
      IF(J.GE.1.AND.J.LE.IY) GO TO 112
      R(JJ)=1E30
      GO TO 111
112   R(JJ)=A(I,J)
111   CONTINUE
      YY=STAY-FIYM2
      CALL BINOM(1.,2.,3.,4.,R(1),R(2),R(3),R(4),YY,SCR(II))
100   CONTINUE
      XX=STAX-FIXM2
      CALL BINOM(1.,2.,3.,4.,SCR(1),SCR(2),SCR(3),SCR(4),XX,STAVAL)
      RETURN
      END
c
cc ------------------------------------------------------------------
c
      subroutine vfirec(iunit,a,n,type)
      character*1 vc
      character*(*) type
      common/vform/vc(0:63)
      character line*80, cs*1
      dimension a(*)

      if(vc(0).ne.'0') call vfinit

      ich0=ichar('0')
      ich9=ichar('9')
      ichcz=ichar('Z')
      ichlz=ichar('z')
      ichca=ichar('A')
      ichla=ichar('a')
      
      read(iunit,10)nn,nbits,bias,fact
 10   format(2i8,2e20.10)
      if(nn.ne.n) then
         print*,' Word count mismatch on vfirec record '
         print*,' Words on record - ',nn
         print*,' Words expected  - ',n
         stop 'vfirec'
      endif

      nvalline=(78*6)/nbits
      nchs=nbits/6
      do 20 i=1,n,nvalline
         read(iunit,'(a78)') line
         ic=0
         do 30 ii=i,i+nvalline-1
            isval=0
            if(ii.gt.n) go to 20
            do 40 iii=1,nchs
               ic=ic+1
               cs=line(ic:ic)
               ics=ichar(cs)
               if(ics.le.ich9)then
                  nc=ics-ich0
               elseif(ics.le.ichcz) then
                  nc=ics-ichca+10
               else
                  nc=ics-ichla+36
               endif
               isval=ior(ishft(nc,6*(nchs-iii)),isval)
 40         continue
            a(ii)=isval
 30      continue
 20   continue

      facti=1./fact
      if(type.eq.'LIN') then
         do 48 i=1,n
            a(i)=a(i)*facti-bias
 48      continue
      elseif(type.eq.'LOG') then
         scfct=2.**(nbits-1)
         do 55 i=1,n
            a(i)=sign(1.,a(i)-scfct)
     +           *(10.**(abs(20.*(a(i)/scfct-1.))-10.))
 55      continue
      endif

      return
      end
c
cc ------------------------------------------------------------------
c
      subroutine vfinit                                                  
      character*1vc,vcscr(0:63)                                         
      common/vform/vc(0:63)                                             
      data vcscr/'0','1','2','3','4','5','6','7','8','9'                 
     +,'A','B','C','D','E','F','G','H','I','J'                          
     +,'K','L','M','N','O','P','Q','R','S','T'                          
     +,'U','V','W','X','Y','Z','a','b','c','d'                          
     +,'e','f','g','h','i','j','k','l','m','n'                          
     +,'o','p','q','r','s','t','u','v','w','x'                          
     +,'y','z','{','|'/                                                 
                                                                        
      do10n=0,63                                                        
      vc(n)=vcscr(n)                                                    
  10  continue                                                          
                                                                        
      return                                                            
      end
C +------------------------------------------------------------------+
      FUNCTION INTLSHFT(IWORD,NSHFT)
C
C       This function shifts IWORD to the left NSHFT bits in a
C         circular manner.
C
      INTLSHFT=ISHFT(IWORD,NSHFT)
      RETURN
      END
C +------------------------------------------------------------------+
c
c determine dominant category 1-05-01 JS
c

      subroutine compute_categories(ctype,nnp,data,nlcat,domcat
     +,pctcat)

      implicit none

      character*(*) ctype

      integer nnp
      integer igc
      integer i,j,k
      integer nlcat
      integer maxcat
      integer lcat(nlcat)
      real    pctcat(nlcat)
      real    data(nnp)
      real    domcat
      real    sum_g

c categorical data types
      if(ctype.eq.'landuse'   .or.
     +   ctype.eq.'soiltype'  .or.
     +   ctype.eq.'islope'  )then

         do k=1,nlcat
            lcat(k)=0
            pctcat(k)=0.
         enddo

         do i=1,nnp
         do k=1,nlcat
            if(nint(data(i)).eq.k)then
               lcat(k)=lcat(k)+1
            endif
         enddo
         enddo

         maxcat=-1
         do k=1,nlcat
            pctcat(k)=lcat(k)/float(nnp)
            if(lcat(k).gt.maxcat)then
               maxcat=lcat(k)
               domcat=float(k)
            endif
         enddo
         if(ctype.eq.'landuse')then
          if(pctcat(16).ge.0.5) then                            !!JBresch
	    domcat = 16                                         !!JBresch
          else if(domcat.eq.16.and.pctcat(16).lt.0.5)then       !!JBresch
            maxcat=-1
            do k=1,nlcat
               if(k.ne.16)then
                  if(lcat(k).gt.maxcat)then
                     maxcat=lcat(k)
                     domcat=float(k)
                  endif
               endif
            enddo
          endif
         endif

c quantitative data types
      elseif(ctype.eq.'greenfrac'.or.ctype.eq.'soiltemp')then

         sum_g=0.
         igc=0
         do i=1,nnp
          if(data(i).gt.0.0)then
             sum_g=sum_g+data(i)
             igc = igc+1
          endif
         enddo
         if(igc.gt.0)then
            domcat=sum_g/float(igc)
         else
            domcat=0.0
         endif

      endif

      return
      end


c ********************************************************************

	subroutine read_dem(unit_no,unit_name,nn1,nn2,i1,i2
     &,data,istat)

	implicit none

	integer countx,county,unit_no,nn1,nn2

        character  cdata(nn1,nn2)*2
        integer,   allocatable :: idata(:,:)
	real       data(nn1,nn2)
	integer len, lend, i1, i2, ia, istat
        real multiplier
c       logical l1,l2
	character*(*) unit_name
        character*1   ctiletype

C	open(unit_no,file=unit_name,status='old',access='direct',
C	. recl=nn2*nn1*2)
C	inquire(unit_no,exist=l1,opened=l2)
C	read(unit_no,rec=1) idata

!	call s_len(unit_name,len)
        len = LEN_TRIM(unit_name)

        call get_directory_length(unit_name,lend)
        ctiletype=unit_name(lend+1:lend+1)

        if(.not.allocated(idata))then
           allocate (idata(nn1,nn2),stat=istat)
           if(istat.ne.0)then
              print*,'unable to allocate idata array: read_dem'
              print*,'nn1/nn2/istat: ',nn1,nn2,istat
              return
           endif
        endif


        multiplier=1.0
	if(ctiletype.eq.'T'.or.ctiletype.eq.'U')then
           if(nn1.eq.1201)then
              open(unit_no,file=unit_name,status='old',
     .form='unformatted')
              read(unit_no)data
              close(unit_no)
           else
              call read_binary_field(cdata,i1,i2,nn1*nn2,unit_name,len)

              do county=1,nn2
              do countx=1,nn1
               idata(countx,county) = ia (cdata(countx,county),2,0)
              enddo
              enddo
              if(ctiletype.eq.'T')multiplier=.01  !(for T data these are temps * 100)
           endif
        else
           call read_binary_field(idata,i1,i2,nn1*nn2,unit_name,len)
        endif

        if(nn1.le.1200)then
	do county=1,nn2
	do countx=1,nn1
	 if(idata(countx,county).eq.-9999) idata(countx,county)=0
	  data(countx,county)=float(idata(countx,nn2-county+1))
     &*multiplier
c SG97 initial data (DEM format) starts in the lower-left corner;
c SG97 this format is wrapped around to have upper-left corner as its start.
c
c JS00 some machines do not account for signed integers
	   if(data(countx,county).ge.15535.0)
     &data(countx,county)=data(countx,county)-65535

	enddo
	enddo
        endif
 
        if(allocated (idata))deallocate(idata)

ccc	 close(unit_no)
	return
	end

c ********************************************************************

        subroutine read_dem_g(unit_no,unit_name,nn1,nn2,nn3,nn4
     &,nofr,i1,i2,data,istat)

        implicit none
        integer  countx,county,countz
        integer  unit_no,nn1,nn2,nn3,nn4,nofr
        integer  len, lend, i1, i2, i
        integer  istat

        real     data(nn1,nn2,nn3,nn4)
        integer, allocatable ::  idata(:,:,:)

c       logical  l1,l2

        character*(*) unit_name
        character*1   ctype

C       open(unit_no,file=unit_name,status='old',access='direct',
C       . recl=nn2*nn1*2)
C       inquire(unit_no,exist=l1,opened=l2)
C       read(unit_no,rec=1) idata

        if(.not.allocated(idata))then
           print*,'allocate idata in read_dem_g'
           allocate (idata(nn4,nn1,nn2),stat=istat)
           if(istat.ne.0)then
              print*,'unable to allocate idata array: read_dem_g'
              print*,'nn1/nn2/nn4/istat: ',nn1,nn2,nn4,istat
              return
           endif
        endif

!        call s_len(unit_name,len)
        len = LEN_TRIM(unit_name)

        call get_directory_length(unit_name,lend)
        ctype=unit_name(lend+1:lend+1)

        print*,'read_dem_g: tile type = ',ctype

        call read_binary_field(idata,i1,i2,nn1*nn2*nn4,unit_name,len)

c       if(nn1.ne.1250 .and. nn2.ne.1250)then
        if(ctype.ne.'A' .and. 
     &     ctype.ne.'G' .and.
     &     ctype.ne.'I' .and.
     &     ctype.ne.'M')then

           do county=1,nn2
           do countx=1,nn1
           do countz=1,nn4

              if(idata(countz,countx,county).eq.-9999)
     &idata(countz,countx,county)=0

              data(countx,county,nofr,countz)=
     &float(idata(countz,countx,nn2-county+1))

c SG97 initial data (DEM format) starts in the lower-left corner;
c SG97 this format is wrapped around to have upper-left corner as its start.


           enddo
           enddo
           enddo

        else   !new greenfrac data starts at 90S

           do county=1,nn2
           do countx=1,nn1
           do countz=1,nn4

              data(countx,county,nofr,countz)=
     &float(idata(countz,countx,county))

           enddo
           enddo
           enddo

        endif

c we'll resurrect the actual categories later
c but for now we want categories 1-9 for
c terrain slope index.
        if(ctype == 'I')then
           do county=1,nn2
           do countx=1,nn1
              if(data(countx,county,1,1).eq.13)then
                 data(countx,county,1,1)=8
              elseif(data(countx,county,1,1).eq.0)then
                 data(countx,county,1,1)=9
              endif
           enddo
           enddo
c          where(data .eq. 13)data = 8
c          where(data .eq. 0)data = 9
        endif

        if(allocated(idata))deallocate(idata)

ccc      close(unit_no)
        return
        end
C
C********************************************************************
C
        subroutine get_directory_length(c_fname,lenf)

cdoc    This routine takes as input a full path filename and returns an
cdoc    index that points at the end of the directory portion of the pathname.
C       Simple-minded algorithm just searches backwards for the first
C       occurance of a `/' (for UNIX) or a `]' (for VMS).
C
C       Input/Output:
C
C       Name            Type    I/O     Description
C       ----            ---     --      -----------
C       c_fnames        char    I       file name.
C       lenf            I       O       index to end of directory
C
C********************************************************************
C
        character c_fname*(*)
        integer*4 lenf

        integer*4 i, strlen
C
C****************************
C
        strlen = len(c_fname)

        i = strlen
        do while (i .gt. 0)
        if( (c_fname(i:i) .ne. ']')
     1.and. (c_fname(i:i) .ne. '/') )then
           i = i-1
        else
           goto 100
        endif
        enddo

100     lenf = i

        return
        end

C*****************************************************************

        function l_string_contains(string,substring,istatus)
 
cdoc    Returns boolean on whether 'string' contains 'substring'
cdoc    The 'string' and 'substring' should be free of blanks.
 
        logical l_string_contains
 
        character*(*) string
        character*(*) substring
 
        l_string_contains = .false.
 
!        call s_len(string,len1)
!        call s_len(substring,len2)
        len1 = LEN_TRIM(string)
        len2 = LEN_TRIM(substring)

        if(len1 .ge. len2)then
            i_search_end = len1 - len2 + 1
 
            do i = 1,i_search_end
                if(string(i:i+len2-1) .eq. substring(1:len2))then
                    l_string_contains = .true.
                endif
            enddo ! i
 
            istatus = 1
            return
 
        else
            istatus = 0
            return
 
        endif
 
        end

c
c===============================================================================
c
      subroutine binom(x1,x2,x3,x4,y1,y2,y3,y4,xxx,yyy)
c
!      include 'bgdata.inc'
!      yyy=missingflag
      yyy = 1.0e37

      if (x2 .gt. 1.e19 .or. x3 .gt. 1.e19 .or.
     .    y2 .gt. 1.e19 .or. y3 .gt. 1.e19) return
c
      wt1=(xxx-x3)/(x2-x3)
      wt2=1.0-wt1
c
      if (y4 .lt. 1.e19 .and. x4 .lt. 1.e19) then
c        yz22=(xxx-x3)*(xxx-x4)/((x2-x3)*(x2-x4))
         yz22=wt1*(xxx-x4)/(x2-x4)
c        yz23=(xxx-x2)*(xxx-x4)/((x3-x2)*(x3-x4))
         yz23=wt2*(xxx-x4)/(x3-x4)
         yz24=(xxx-x2)*(xxx-x3)/((x4-x2)*(x4-x3))
      else
         yz22=wt1
         yz23=wt2
         yz24=0.0
      endif
c
      if (y1 .lt. 1.e19 .and. x1 .lt. 1.e19) then
         yz11=(xxx-x2)*(xxx-x3)/((x1-x2)*(x1-x3))
c        yz12=(xxx-x1)*(xxx-x3)/((x2-x1)*(x2-x3))
         yz12=wt1*(xxx-x1)/(x2-x1)
c        yz13=(xxx-x1)*(xxx-x2)/((x3-x1)*(x3-x2))
         yz13=wt2*(xxx-x1)/(x3-x1)
      else
         yz11=0.0
         yz12=wt1
         yz13=wt2
      endif
c
      if (yz11 .eq. 0. .and. yz24 .eq. 0.) then
         yyy=wt1*y2+wt2*y3
      else
         yyy=wt1*(yz11*y1+yz12*y2+yz13*y3)+wt2*(yz22*y2+yz23*y3+yz24*y4)
      endif
c
      return
      end

!----------------------------------------------------------------------
        FUNCTION IA(CHR,N,ISPVAL)
!                                                              
!  PURPOSE: TO CONVERT A N-BYTES CHARACTER (CHR) TO INTEGER IA. 
!        ** THE INTEGER DATA FILE IS SAVED AS A N-BYTE CHARACTER
!           DATA FILE. THIS FUNCTION IS USED TO RECOVER THE    
!           CHARACTER DATA TO THE INTEGER DATA.               
!                                                            
!  N      --- THE NUMBER OF BYTES IN CHR                    
!  ISPVAL --- DEFAULT VALUE FOR THE NEGATIVE INTEGER.      
!                                                       
        CHARACTER*(*) :: CHR
        integer  N, II1, II2, JJ, ISN, M, NBIT, MSHFT, IA2, ispval
        INTEGER  BIT_1, BIT_2
!                                                    
        BIT_1 = O'200'     ! BINARY '10000000'        
        BIT_2 = O'377'     ! BINARY '11111111'       
        IA    = 0
!                                                
        II1 = ICHAR(CHR(1:1))
        if(II1 < 0) II1=II1+256
 
! .. GET THE SIGN -- ISN=0 POSITIVE, ISN=1 NEGATIVE:
        JJ  = IAND(II1,BIT_1)
        ISN = ISHFT(JJ,-7)
!                                                
! .. FOR NEGATIVE NUMBER:
!    BECAUSE THE NEGATIVE INTEGERS ARE REPRESENTED BY THE SUPPLEMENTARY
!    BINARY CODE INSIDE MACHINE.
!                              
        IF (ISN.EQ.1) THEN
          DO M = N+1,4
           NBIT = (M-1)*8
           JJ = ISHFT(BIT_2,NBIT)
           IA = IEOR(JJ,IA)
          END DO
        ENDIF
!                              
!   .. GET THE BYTE FROM CHR: 
        DO M = 1,N
         II2 = ICHAR(CHR(M:M))
         if(II2 < 0) II2=II2+256
         MSHFT = (N-M)*8
         IA2   = ISHFT(II2,MSHFT)
!   .. THE ABS(INTEGER):          
         IA = IEOR(IA,IA2)
        END DO
!                              
        IF (IA.LT.0) IA = ISPVAL
!                            
        RETURN
        END

! See WRFSI src/grid/proc_geodat.f
!
      subroutine proc_geodat(nx_dom,ny_dom,ncat
     1,path_to_tile_data,dom_lats_in,dom_lons_in,lmask_out
     1,geodat,istatus)
!    1,cat_pct)
!
! J. Smart (NOAA/FSL) : 2002            original version
! T. Hume/T. Simmers Meteorological Service of New Zealand
!                     : 2002
!                    Corrected problems with crossing the dateline
!                    and constrained albedo to one tile.
! J. Smart (NOAA/FSL) : 2003
!                    Further refined processing of tiles to accomodate
!                    dateline, Greenwich mean, and Poles. Added smoothing
!                    of temp field (subroutine one_two_one).

!      use horiz_interp

      implicit none 
      integer maxtiles 
      integer nx_dom 
      integer ny_dom
      integer ncat
      integer ntn,nt
      integer bgmodel
      integer itilesize_d
      integer lentd,lenp,lenf
      integer iblksizo
      integer isbego,iwbego
      integer no
      integer icurNS,icurEW
      integer itile
      integer itilex
      integer itiley
      integer,allocatable:: itile_lat(:)
      integer,allocatable:: itile_lon(:)
      integer itilelon
      integer it,jt
      integer is,js
      integer ie,je
      integer itx,ity
      integer i,j,k,ii,jj
      integer ix,iy
      integer icnta
      integer tile_n_lat
      integer tile_s_lat
      integer tile_w_lon
      integer tile_e_lon
      integer istatus
      integer istat

      integer nx_tile,ny_tile
      integer dom_i,dom_j
      integer point_value
      integer ilon,ilat
      integer method

      real    dom_lats_in(nx_dom,ny_dom)
      real    dom_lons_in(nx_dom,ny_dom)
      real    grx(nx_dom,ny_dom)
      real    gry(nx_dom,ny_dom)

!Define Array to keep count of number of raw tile data points found for each
!model grid point

!     real    cat_pct(nx_dom,ny_dom,ncat)
!     integer, allocatable :: num_raw_points_total (:,:)
!     integer, allocatable :: num_raw_points_cat(:,:,:)

      real,    allocatable :: raw_data(:,:,:)
      real,    allocatable :: data_proc(:,:,:)
      real,    allocatable :: lmask_src(:,:)
      real,    allocatable,   dimension(:,:)::dom_lats,dom_lons

      real    geodat(nx_dom,ny_dom,ncat)
      real    lmask_out(nx_dom,ny_dom)
      real    amean(ncat)
      real    asum(ncat)

      real    dlat_tile
      real    dlon_tile
      real    sw(2),ne(2)
      real    deltallo
      real    min_lat
      real    max_lat
      real    min_lon
      real    max_lon
      real    minlat
      real    maxlat
      real    minlon
      real    maxlon
      real    rwoff,rsoff
      real    rlondif
      real    offset
      real    ri,rj
      real    rmult
      real    min_val,max_val,def_val,val_mask
      real    r_missing_data
 
      logical  dem_data
      logical  make_srcmask
      logical  lexist,lopen

      character*(*) path_to_tile_data
      character*255 title
      character*255 cfname
      character(len=8),allocatable:: ctile_name_list(:)
      character*8   ctilename
      character*1   curNS,curEW
      character*1   ctiletype

      character*1   cgrddef
      integer nxst,nyst,nz
      real    lat0,lon0,dlat,dlon
!      common /llgrid/nxst,nyst,nz,lat0,lon0,dlat,dlon,cgrddef

c original create from Brent Shaw pseudocode

      print*,'Start proc_geodat'

      if(.not. allocated(dom_lats))then
          allocate(dom_lats(nx_dom,ny_dom),dom_lons(nx_dom,ny_dom))
      endif
      dom_lats=dom_lats_in
      dom_lons=dom_lons_in

! TH: 8 Aug 2002 Now we may need to adjust the longitude values in 
! dom_lons so that they are always monotonically increasing, even if 
! we have the bad fortune to cross the date line.
      DO jj=1,ny_dom-1
         DO ii=1,nx_dom
            IF ((dom_lons(ii,jj+1) - dom_lons(ii,jj)) < -180.) THEN
               dom_lons(ii,jj+1) = dom_lons(ii,jj+1) + 360.
            ELSE IF ((dom_lons(ii,jj+1) - dom_lons(ii,jj)) > 180.) THEN
               dom_lons(ii,jj+1) = dom_lons(ii,jj+1) - 360.
            END IF
         END DO
      END DO
! TH: 8 Aug 2002 We also have to cope with spastic grids that have
! "horizontal" date lines (where we might cross the date line by
! moving up and down the grid rather than left and right).
      DO jj=1,ny_dom
         DO ii=1,nx_dom-1
            IF ((dom_lons(ii+1,jj) - dom_lons(ii,jj)) < -180.) THEN
               dom_lons(ii+1,jj) = dom_lons(ii+1,jj) + 360.
            ELSE IF ((dom_lons(ii+1,jj) - dom_lons(ii,jj)) > 180.) THEN
               dom_lons(ii+1,jj) = dom_lons(ii+1,jj) - 360.
            END IF
         END DO
      END DO
! TH: 8 Aug 2002 Now we want to make sure our longitudes don't fall
! outside the range (-360,360).
      IF (MAXVAL(dom_lons(:,:)) > 360.) dom_lons(:,:) = dom_lons(:,:)
     1   - 360
      max_lat = MAXVAL(dom_lats)
      IF (MINVAL(dom_lons(:,:)) < -360.) dom_lons(:,:) = dom_lons(:,:)
     1   + 360

! nx_dom = number of East-West points in anal/model domain   I
! ny_dom = number of North-south points      "               I
! dom_lats(nx_dom,ny_dom)  ! Latitude array                  I
! dom_lons(nx_dom,ny_dom)  ! Longitude array                 I

!c      call s_len(path_to_tile_data,lenp)
      lenp = LEN_TRIM(path_to_tile_data)

      ctiletype=path_to_tile_data(lenp:lenp)
      if(ctiletype.eq.'T'.or.
     1   ctiletype.eq.'A'.or.
     1   ctiletype.eq.'G'.or.
     1   ctiletype.eq.'M' )then
         print*,'tile type to process = ',ctiletype
      else
         print*,'Unknown tile type in proc_geodat_tiles'
         print*,'path_to_tile_data: ',path_to_tile_data(1:lenp)
         stop
      endif

      bgmodel=6

      TITLE=path_to_tile_data(1:lenp)//'HEADER'
      lentd=INDEX(TITLE,' ')-1
!c     CALL JCLGET(29,TITLE(1:lentd),'FORMATTED',1,istatus)
      OPEN(29,FILE=title(1:lentd),FORM='FORMATTED',STATUS='old')

      if(istatus .ne. 1)then
         write(6,*)'Warning: proc_geodat_tiles opening HEADER: check'
     1            ,'geog paths and HEADER file'
         return
      endif

      READ(29,*)IBLKSIZO,NO,ISBEGO,IWBEGO,RWOFF,RSOFF
      print *,'title=',title(1:lentd)
      print *,'RWOFF,RSOFF = ',RWOFF,RSOFF
      print *,'isbego,iwbego=',isbego,iwbego
      print *,'iblksizo,no=',iblksizo,no
      CLOSE(29)

      if(ctiletype.eq.'T'.and.iblksizo.eq.10)then
       print*
       print*,'Error: ****************************************'
       print*,'Error: ********** TERMINATING *****************'
       print*,'Error: Old Soil Temp database!!!'
       print*,'Error: You need to update these'
       print*,'Error: data: ',TRIM(path_to_tile_data)
       print*,'Error: by downloading new soil temp data'
       print*,'Error: ********** TERMINATING *****************'
       print*,'Error: ****************************************'
       print*
      endif

      maxtiles=(360/iblksizo)*(180/iblksizo)
      print*,'Maxtiles = ',maxtiles

      if(NO .GT. 0)then
         DELTALLO=FLOAT(IBLKSIZO)/FLOAT(NO)
      else
         print *,'ERROR:  NO <= 0 ... '
         stop
      endif

      itilesize_d = iblksizo

! Define number of points in a raw tile and what the lat/lon increment
! between the points in a tile are

      nx_tile=no
      ny_tile=no
      dlat_tile=deltallo
      dlon_tile=dlat_tile

! ncat ! Number of categories (1 for topo, 24 for veg, etc.)
!Define Array to keep count of number of raw tile data points found for each
!model grid point

! Find min/max atitude and longitude so we can compute which tiles to 
! read

      minlat = MINVAL(dom_lats)
      maxlat = MAXVAL(dom_lats)
      minlon = MINVAL(dom_lons)
      maxlon = MAXVAL(dom_lons)

! no offsets needed here since these are the tiles, not the points
! in the tiles. However, since domains may need data from tiles with
! offset considered, add it in.
      min_lat = max(-89.9999, min(89.9999, minlat - abs(rsoff)))
      max_lat = max(-89.9999, min(89.9999, maxlat + abs(rsoff)))
      min_lon = max(-359.9999,min(359.9999,minlon - abs(rwoff)))
      max_lon = max(-359.9999,min(359.9999,maxlon + abs(rwoff)))

c     if(ctiletype .eq. 'G'.or.ctiletype.eq.'A')then
c        if(min_lon.lt.-180)min_lon=360.+min_lon
c     endif

      print*,'max/min lats: ',max_lat,min_lat
      print*,'max/min lons: ',max_lon,min_lon

      deallocate(dom_lats,dom_lons)

!  Compute a list of tiles needed to fulfill lat/lon range just computed

      print*,'allocate ctile_name_list lat/lon arrays'
      print*,'with maxtiles = ',maxtiles

      allocate (ctile_name_list(maxtiles))

      call get_tile_list(min_lat,max_lat,min_lon,max_lon
     1,maxtiles,isbego,iwbego,itilesize_d,ctiletype
     1,ntn,ctile_name_list,tile_w_lon,tile_e_lon
     1,tile_s_lat,tile_n_lat,istatus)

      if(istatus.ne.1)then
         print*,'ERROR:  returned from get_tile_list'
         return
      endif

      allocate (raw_data(nx_tile,ny_tile,ncat))
      allocate (itile_lat(ntn))
      allocate (itile_lon(ntn))

!     ALLOCATE(num_raw_points_total(nx_dom,ny_dom))
!     ALLOCATE(num_raw_points_cat(nx_dom,ny_dom,ncat))
!     num_raw_points_total(:,:)=0
!     num_raw_points_cat(:,:,:)=0

      do itile=1,ntn
         read(ctile_name_list(itile)(1:2),'(i2.2)')itile_lat(itile)
         read(ctile_name_list(itile)(4:6),'(i3.3)')itile_lon(itile)
         if(ctile_name_list(itile)(7:7).eq.'W')then
            if(itile_lon(itile).lt.180)then
               itile_lon(itile)=360-itile_lon(itile)
            endif
         endif
         if(ctile_name_list(itile)(3:3).eq.'S')then
            itile_lat(itile)=-1.0*itile_lat(itile)
         endif
      enddo

      print*,'S/N tile lat points = ',tile_s_lat,tile_n_lat
      print*,'W/E tile lon points = ',tile_w_lon,tile_e_lon

c determine the x/y size dimensions and allocate/initialize the super tile

      rlondif=abs(tile_e_lon-tile_w_lon)
      itx=nx_tile*nint((rlondif+itilesize_d)/float(itilesize_d))
      ity=ny_tile*nint((abs(tile_n_lat-tile_s_lat)+itilesize_d)
     ./float(itilesize_d))
      if(ctiletype.eq.'T' .or. ctiletype.eq.'M' .and.
     &itx.gt.360)itx=360
      if(ctiletype.eq.'G' .or. ctiletype.eq.'A' .and. 
     &itx.gt.2500)itx=2500
      print*,'allocate data_proc: nx/ny/ncat ',itx,ity,ncat
      allocate(data_proc(itx,ity,ncat))
!      call get_r_missing_data(r_missing_data,istatus)
      r_missing_data = 1.e+37
      data_proc = r_missing_data

c for the supertile lat-lon "grid" definition:
c the LL setup must consider the offsets since this
c is relevant to the actual data points within the tile.
      nxst=itx
      nyst=ity
      dlat=dlat_tile
      dlon=dlat
      lat0=tile_s_lat+rsoff
      if(tile_w_lon.gt.180)then
         lon0=tile_w_lon-360+rwoff 
      elseif(tile_w_lon.lt.-180)then
         lon0=360+tile_w_lon+rwoff
      else
         lon0=tile_w_lon+rwoff
      endif
      sw(1)=lat0
      sw(2)=lon0
      ne(1)=tile_n_lat+itilesize_d+rsoff
      ne(2)=tile_e_lon+itilesize_d+rwoff
      cgrddef='S'

      print*,'generate supertile for domain'
      print*,'number of small tiles needed= ',ntn
      print*

      DO itile = 1, ntn  !number of tiles needed
 
       ctilename=ctile_name_list(itile)
       cfname = path_to_tile_data(1:lenp)//ctilename(1:3)
       read(ctilename(4:6),'(i3.3)')icurEW
       read(ctilename(7:7),'(a1)')curEW
       IF (icurEW > 180)THEN    ! .and. icurEW < 360) THEN
          icurEW = 360 - icurEW
       ENDIF
       IF (icurEW == 180.and.curEW.ne.'W') THEN
          curEW = 'W' 
       END IF 

       write(cfname(lenp+4:lenp+6),'(i3.3)')icurEW
       write(cfname(lenp+7:lenp+7),'(a1)')curEW

!       call s_len(cfname,lenf)
       lenf = LEN_TRIM(cfname)

       print*,'processing tile: ',cfname(1:lenf)

       inquire(file=cfname,exist=lexist,opened=lopen)
       if(.not.lexist)then
          print*,'Error: file does not exist: ',cfname(1:lenp)
          goto 1000
       endif
       if(lopen)then
          print*,'Error: file already open: ',cfname(1:lenp)
          goto 1000
       endif


! Open the tile and read the points
       if( (ctiletype.eq.'U').and.(no.eq.1200) )then
           CALL READ_DEM(29,cfname,no,no,2,2,raw_data,istat) ! world topo_30s
           dem_data=.true.
       elseif( ctiletype.eq.'O' )then      ! soiltype top and bot layer
           CALL READ_DEM(29,cfname,no,no,1,4,raw_data,istat)
           dem_data=.true.
       elseif( ctiletype.eq.'V' )then      ! world USGS 30s landuse
           CALL READ_DEM(29,cfname,no,no,1,4,raw_data,istat)
           dem_data=.true.
       elseif( (ctiletype.eq.'G')
     1     .or.(ctiletype.eq.'A')
     1     .or.(ctiletype.eq.'M') )then      ! greenfrac/albedo/maxsnowalb
           CALL READ_DEM_G(29,cfname,no,no,1,ncat,1,1,4,raw_data,istat)
           dem_data=.true.
       elseif( ctiletype.eq.'T')then ! .or. ctiletype.eq.'M')then      ! soiltemp
           CALL READ_DEM(29,cfname,no,no,2,2,raw_data,istat)
           dem_data=.true.
       else                                ! other  like albedo
!           CALL JCLGET(29,cfname,'FORMATTED',0,istatus)
           OPEN(29,FILE=cfname,FORM='FORMATTED',STATUS='OLD')
           CALL VFIREC(29,raw_data,no*no,'LIN')
           if ((ctiletype.eq.'U').and.(no.eq.121)) then
                dem_data=.false.           ! topo_30s
           endif
       endif

       if(istat.ne.0)then
          print*,'Error returned: proc_geodat: READ_DEM'
          istatus=0
          return
       endif
c
c make "super tile" from all smaller tiles
c
       read(ctilename(1:2),'(i2.2)')icurNS
       read(ctilename(3:3),'(a1)')curNS
       read(ctilename(4:6),'(i3.3)')icurEW
       read(ctilename(7:7),'(a1)')curEW

c      print*,'compute itiley/itilex'
       if(curNS .eq. 'S') icurNS=-1.0*icurNS
       itiley=nint(1.0+(float(icurNS)-tile_s_lat)/float(itilesize_d))
c
       rlondif=float(itile_lon(itile))-tile_w_lon
       if(rlondif.lt.0)then
          rlondif=rlondif+360.
       elseif(rlondif.ge.360)then
          rlondif=rlondif-360
       endif

       itilex=nint(1.0+(rlondif/float(itilesize_d)))

       itx=1+(itilex-1)*nx_tile
       ity=1+(itiley-1)*ny_tile

       print*,'itilex/itiley ',itilex,itiley
       print*,'itx/ity ',itx,ity

       jj=0
       do iy=ity,ny_tile*itiley
          jj=jj+1
          ii=0
          do ix=itx,nx_tile*itilex
             ii=ii+1
             data_proc(ix,iy,:)=raw_data(ii,jj,:)
          enddo
       enddo

      ENDDO

      deallocate (raw_data,itile_lat,itile_lon,ctile_name_list)

      print*,'initializing hinterp grx/gry arrays'

      print*,'nxst= ',nxst
      print*,'nyst= ',nyst
      print*,'dlat= ',dlat
      print*,'dlon= ',dlon
      print*,'lat0= ',lat0
      print*,'lon0= ',lon0

      call init_hinterp(nxst,nyst,nx_dom,ny_dom,'LL',
     .dom_lats_in,dom_lons_in,grx,gry,bgmodel,
     &lat0,lon0,dlat,dlon,cgrddef)

      print*,'grid rx/ry computed'
      print*,'SW: grx/gry    1,1: ',grx(1,1),gry(1,1)
      print*,'SE: grx/gry   nx,1: ',grx(nx_dom,1),gry(nx_dom,1)
      print*,'NW: grx/gry   1,ny: ',grx(1,ny_dom),gry(1,ny_dom)
      print*,'NE: grx/gry  nx,ny: ',grx(nx_dom,ny_dom)
     .,gry(nx_dom,ny_dom)

c compute mean value to use as def_value

      is=nint(MINVAL(grx))
      ie=nint(MAXVAL(grx))
      js=nint(MINVAL(gry))
      je=nint(MAXVAL(gry))

      print*,'is/ie/js/je ',is,ie,js,je

      do k=1,ncat
         asum(k)=0.0
         icnta=0
         do j=js,je
         do i=is,ie
            if(data_proc(i,j,k).ne.0.0)then
              icnta=icnta+1
              asum(k)=asum(k)+data_proc(i,j,k)
            endif
         enddo
         enddo
         if(icnta.gt.0)then
            amean(k)=asum(k)/icnta
         else
            amean(k)=r_missing_data
         endif
      enddo

      allocate (lmask_src(nxst,nyst))
      make_srcmask=.true.
      val_mask = 1
      nt=1

      if(ctiletype.eq.'T')then

         min_val = 239.73
         max_val = 305.09
         method  = 1
         geodat(:,:,:)=r_missing_data

      elseif(ctiletype.eq.'A')then

         min_val = 2.0
         max_val = 100.
         method  = 2
         where(amean .eq. r_missing_data)amean=18.

      elseif(ctiletype.eq.'G')then

         min_val = 1.0
         max_val = 100.
         method  = 2
         where(amean .eq. r_missing_data)amean=65.

      elseif(ctiletype.eq.'M')then

         min_val = 1.0
         max_val = 100.
         method  = 1
         geodat(:,:,:)=65.

      endif

      do ii=1,ncat

         def_val = amean(ii)

         call interpolate_masked_val(nxst, nyst
     ., lmask_src, data_proc(1,1,ii), nx_dom, ny_dom, lmask_out
     ., geodat(1,1,ii), grx, gry, make_srcmask, min_val, max_val
     ., def_val, val_mask, method)

         if(ctiletype.eq.'A'.or.ctiletype.eq.'M')then
            geodat(:,:,ii)=geodat(:,:,ii)/100.
         endif

         if(ctiletype.eq.'T')then  !.or.ctiletype.eq.'M')then
            print*,'Filering Deep Soil Temp with 1-2-1'
            call one_two_one(nx_dom,ny_dom,nt,geodat(1,1,ii))
         endif

      enddo

      deallocate (data_proc,lmask_src)

! --------------------------------------------------------
! ---- we are not yet processing categorical data (like
! ---- landuse or soil type, so this section is commented
! ---- for now.
! --------------------------------------------------------
! Loop over the raw data points in the tile

!       DO jt = 1, ny_tile
!       DO it = 1, nx_tile

! Compute lat/lon of this raw data point using lat/lon from
! file name, the offset of the data, and the increment

!        raw_lat = float(ilat) + RSOFF + (jt-1)*dlat_tile
!        raw_lon = float(ilon) + RWOFF + (it-1)*dlon_tile

! Use lat/lon to ij to find wheter this tile point is in our domain

!        call latlon_to_rlapsgrid(raw_lat,raw_lon,dom_lats_in
!    &                ,dom_lons_in,nx_dom,ny_dom,ri,rj,istatus)
!        dom_i = NINT(ri)
!        dom_j = NINT(rj)
!        IF( (dom_i .GE. 1).AND.(dom_i .LE. nx_dom).AND.
!    &   (dom_j .GE. 1).AND.(dom_j .LE. ny_dom) ) THEN

!         if(ctiletype.eq.'V')then

!            point_value=int(raw_data(it,jt,1))
!            num_raw_points_cat(dom_i,dom_j,point_value) =
!    &       num_raw_points_cat(dom_i,dom_j,point_value) + 1

!            num_raw_points_total(dom_i, dom_j) =
!    &       num_raw_points_total(dom_i, dom_j)+ 1

!         elseif(ctiletype.eq.'T')then

!            num_raw_points_cat(dom_i,dom_j,1) =
!    &       num_raw_points_cat(dom_i,dom_j,1)+1

!            num_raw_points_total(dom_i, dom_j) =
!    &       num_raw_points_total(dom_i, dom_j)+raw_data(it,jt,1)

!         endif
         
!         ! The point is in our domain, so increment the counters.  The 
          ! value of the raw data point is equal to the category ID for
          ! things like veg-type and soil-type.  For topography, we would
          ! have to have additional arrays to keep the sum of the terrain
          ! values to compute the average value at the end. This example
          ! code only does category data

!        ENDIF
!       ENDDO
!       ENDDO 


! Now you can compute percentage of each category like so:
!     if(ctiletype.eq.'V')then
!        DO icat = 1, ncat
!           cat_pct(:,:,icat)=num_raw_points_cat(:,:,icat) / 
!    &                         num_raw_points_total(:,:)
!        ENDDO
!     elseif(ctiletype.eq.'T')then
!           cat_pct(:,:,1)= num_raw_points_total(:,:)/
!    &                       num_raw_points_cat(:,:,1)
!     endif

      return

1000  print*,'returning to main. no data processed'
      return

      end
c
c-----------------------------------------------------------
c
      subroutine one_two_one(nx,ny,nt,data)

      implicit none
      integer i,j,ii,jj,icnt,n
      integer nx,ny,nt
      integer istatus
      real  r_missing_data
      real  fact,sum
      real  data(nx,ny)
      real, allocatable :: temp(:,:)

      allocate(temp(nx,ny))

!      call get_r_missing_data(r_missing_data,istatus)
      r_missing_data = 1e+37

      temp=data

      do n=1,nt

         do j=2,ny-1
         do i=2,nx-1

            sum=0.0
            icnt=0
            do jj=j-1,j+1
            do ii=i-1,i+1

               fact=1
               if(jj==j.and.ii==i)fact=2
               if(data(ii,jj).lt.r_missing_data)then
                  icnt=icnt+fact
                  sum=sum+data(ii,jj)*fact
               endif

            enddo
            enddo

            if(icnt.gt.0)then
               temp(i,j)=sum/float(icnt)
            else
               temp(i,j)=data(i,j)
            endif

         enddo
         enddo

         data=temp

      enddo

      deallocate (temp)

      return
      end

      subroutine get_tile_list(min_lat,max_lat,min_lon,max_lon
     1,maxtiles,isbego,iwbego,itilesize,ctiletype
     1,num_tiles_needed,ctile_name_list,iwoc1,iwoc2
     1,isoc1,isoc2,istatus)

      implicit  none

      real      max_lat,min_lat
      real      max_lon,min_lon

      integer   num_tiles_needed
      integer   isbego,iwbego
      integer   itilesize
      integer   itile_ns
      integer   iwoc0
      integer   iwoc,iwoc1,iwoc2
      integer   isoc,isoc1,isoc2
      integer   isocpt,isocpo
      integer   iwocpt,iwocpo,iwocph
      integer   iwocdif
      integer   istatus
      integer   maxtiles


      character*1 ctiletype
      character*3 nstitle
      character*4 ewtitle
      character*8 ctilenamelist(500)
      character*(*) ctile_name_list(maxtiles)

      double precision r8term

      istatus = 1

      r8term=(min_lat-float(isbego))/float(itilesize)
      ISOC1=(INT(r8term+200.)-200)*itilesize+ISBEGO
      r8term=(min_lon-float(iwbego))/float(itilesize)
      IWOC1=(INT(r8term+400.)-400)*itilesize+IWBEGO
      r8term=(max_lat-float(isbego))/float(itilesize)
      ISOC2=(INT(r8term+200.)-200)*itilesize+ISBEGO
      r8term=(max_lon-float(iwbego))/float(itilesize)
      IWOC2=(INT(r8term+400.)-400)*itilesize+IWBEGO

      num_tiles_needed=0

      if(IWOC1.lt.-180)then
         if(itilesize.eq.180)then
            iwocdif=iwoc2-iwoc1
            IWOC1=IWOC1+IWOCDIF
            IWOC2=IWOC2+abs(IWOCDIF)
c        else
c           IWOC1=360+IWOC1
         endif
c     elseif(IWOC1.gt.180)then
c        IWOC1=360-IWOC1
      endif
      print*,'Noddy IWOC1, IWOC2 ',IWOC1,IWOC2
      do IWOC = IWOC1,IWOC2,itilesize

         IWOC0 = IWOC
         IF(IWOC.LT.-180)IWOC0=360+IWOC0
c        IF(IWOC.GT.+180)IWOC0=360-IWOC0

         IWOCPH=ABS(IWOC0)/100
         IWOCPT=(ABS(IWOC0)-IWOCPH*100)/10
         IWOCPO=ABS(IWOC0)-IWOCPH*100-IWOCPT*10
!
! TH: 8 Aug 2002 We now allow 180E longitudes (and greater). The only 
! time we want to assign W is when the longitude is less than 0.
!
         IF(IWOC0.GE.0.and.IWOC0.LT.180) THEN
            WRITE(EWTITLE,'(3I1,A1)')IWOCPH,IWOCPT,IWOCPO,'E'
         ELSE
            WRITE(EWTITLE,'(3I1,A1)')IWOCPH,IWOCPT,IWOCPO,'W'
         ENDIF

         if(ewtitle(1:1).eq.' ')ewtitle(1:1)='0'
         if(ewtitle(2:2).eq.' ')ewtitle(2:2)='0'

c        ewtitle=ewtitle2//ewtitle1

         do ISOC = ISOC1,ISOC2,itilesize

            ISOCPT=ABS(ISOC)/10
            ISOCPO=ABS(ISOC)-ISOCPT*10

            IF(ISOC.GE.0)THEN
              WRITE(NSTITLE,'(2I1,A1)')ISOCPT,ISOCPO,'N'
            ELSE
              WRITE(NSTITLE,'(2I1,A1)')ISOCPT,ISOCPO,'S'
            ENDIF

            num_tiles_needed=num_tiles_needed+1
            ctilenamelist(num_tiles_needed)=nstitle//ewtitle

         enddo
      enddo

      if(num_tiles_needed.eq.3.and.maxtiles.eq.2)then
         num_tiles_needed = 2
      endif

      if(num_tiles_needed.le.maxtiles)then
         do itile_ns=1,num_tiles_needed
            ctile_name_list(itile_ns)=ctilenamelist(itile_ns)
         enddo
      else
         print*,'more tiles than array allocation'
         istatus = 0
      endif

      return
      end


      subroutine adjust_geog(nnxp,nnyp,ncat,ctype
     &,istat_dat,lat,topt_out           !istattmp,istatslp,lat,topt_out,path_to_soiltemp
     &,landmask,geog_data,istatus)      !soiltemp_1deg,greenfrac,islope,istatus)
c
c This routine uses landmask to make the course resolution
c soil temp and green fraction data conform to water-land mask.
c It also fills small islands or isolated land
c bodies with appropriate geog values as necessary
c
c J. Smart NOAA/FSL
c     "    06-22-01 original
c     "    12-01-03 put subroutine into gridgen_utils, make arrangements for soil
c                   type data, and allow separate calls for each geog type to reduce
c                   memory requirement.
c
      implicit none

      integer nnxp,nnyp
      integer i,j,l,ii,jj
      integer is,js
      integer istatus
      integer istat_dat    !Input. =0 indicates input geog data was NOT processed properly; =1 otherwise.
      integer ijsthresh    !search distance to look for representative vlues to fill inconsistency.
      integer isum
      integer l1,l2
      integer lent
      integer ncat
      integer ifixw,ifixl
      integer ifixws,ifixwt,ifixls

      integer isc
      integer ic(ncat) 

      character*132 path_to_soiltemp
      character*(*) ctype

      logical endsearch

      real,intent(inout)  ::    lat(nnxp,nnyp)
      real,intent(inout)  ::    landmask(nnxp,nnyp)
      real,intent(inout)  ::    topt_out(nnxp,nnyp)
      real,intent(inout)  ::    geog_data(nnxp,nnyp,ncat)

c     real,intent(inout)  ::    soiltemp_1deg(nnxp,nnyp)
c     real,intent(inout)  ::    greenfrac(nnxp,nnyp,ncat)
c     real,intent(inout)  ::    islope(nnxp,nnyp)

      real,   allocatable ::    geog_tmp      (:,:,:) !temporary holder of input geog data.

c     real,   allocatable ::    grnfrctmp     (:,:,:)
c     real,   allocatable ::    soiltmp       (:,:)
c     real,   allocatable ::    islopetmp     (:,:)

      real,   allocatable ::    rmeanlattemp  (:)

      real    avgtmp
      real    avggrn(ncat)
      real    avgcat
      real    tatlmx,tatlmn
      real    rlatmx,rlatmn
      real    r_missing_data
      real    rmngrn
      real    islp
      real    sumt
      real    sumg
      real    sum(ncat)
      real    tslp

      istatus=0
      if(istat_dat.eq.0)then  !.and.istattmp.eq.0.and.istatslp.eq.0)then
	 print*,'Unable to process geog data in adjust_geog ...'
     &,' processsing of data failed prior to subroutine call.'
	 return
      endif


c use moist adiabatic laps rate (6.5 deg/km) to get new temp
 
      allocate (geog_tmp(nnxp,nnyp,ncat))

!      call get_r_missing_data(r_missing_data,istatus)
      r_missing_data = 1e+37

      geog_tmp=geog_data

      where(geog_tmp.eq.0.0)geog_tmp=r_missing_data

      if(ctype.eq.'greenfrac')then
         rmngrn=minval(geog_tmp(1:nnxp,1:nnyp,1:1))
         print*,'minimum green fraction = ',rmngrn
      endif

c determine average soiltemp and greenfrac in domain
      sumt=0.0
      sum=0.0
      isc=0
      ic=0
      islp=0
      isum=0
      if(ctype.eq.'soiltemp')then
         do j = 1,nnyp
         do i = 1,nnxp
            if(geog_tmp(i,j,1).ne. r_missing_data)then
               sumt=geog_tmp(i,j,1)+sumt
               isc=isc+1
            endif
         enddo
         enddo
      elseif(ctype.eq.'greenfrac')then
c greenfrac is assumed to be continuous globally; only use land points
         do j = 1,nnyp
         do i = 1,nnxp
            do l=1,ncat
               if(landmask(i,j) .ne. 0 .and.
     &geog_tmp(i,j,l).lt.r_missing_data)then
                  sum(l)=geog_tmp(i,j,l)+sum(l)
                  ic(l)=ic(l)+1
               endif
            enddo
         enddo
         enddo
      elseif(ctype.eq.'islope'.or.ctype.eq.'soiltype')then
         do j = 1,nnyp
         do i = 1,nnxp
            if(geog_tmp(i,j,1).ne. r_missing_data)then
               isum=geog_tmp(i,j,1)+isum
               islp=islp+1
            endif
         enddo
         enddo
      endif
c
c get average information for filling as necessary
c ----------------------------------------------------
c
c 1. deep soil temp.
c
      if(isc.gt.0)then
         avgtmp=sumt/float(isc)
         print*,'Domain average annual mean temp = ',avgtmp
      elseif(ctype.eq.'soiltemp')then
c
c this section uses mean latitudinally averaged temps derived from 
c the raw 1 deg temp data. File in raw geog annual mean deep soil
c temp directory (see wrfsi.nl variable soiltemp_1deg).
c
         call get_path_to_soiltemp_1deg(path_to_soiltemp,istatus)
         sumt=0.0
         print*,'*** Using mean lat temps -> LATMEANTEMP ***'
         allocate (rmeanlattemp(180))
         call  get_directory_length(path_to_soiltemp,lent)
         call  get_meanlattemp(path_to_soiltemp(1:lent-1)
     &,rmeanlattemp,istatus)
         if(istatus.ne.1)then
            print*,'Error returned: get_meanlattemp'
            return
         endif
         l1=90-nint(minval(lat(:,1)))+1
         l2=90-nint(maxval(lat(:,nnyp)))+1
         if(l1.gt.180)l1=180
         if(l2.gt.180)l2=180

         rlatmn=minval(lat(:,1))
         rlatmx=maxval(lat(:,nnyp))
         tatlmx=rmeanlattemp(l2)
         tatlmn=rmeanlattemp(l1)
         if(tatlmx.ne.0.and.tatlmn.ne.0.and.
     .      tatlmx.ne.tatlmn)then
            tslp=(rlatmx-rlatmn)/(tatlmx-tatlmn)
            do j=1,nnyp
            do i=1,nnxp
               geog_tmp(i,j,1)=tatlmn+(rlatmx-lat(i,j))/tslp
               sumt=geog_tmp(i,j,1)+sumt
            enddo
            enddo
            avgtmp=sumt/(nnyp*nnxp)
         elseif(tatlmx.gt.0)then
            avgtmp=tatlmx
         elseif(tatlmn.gt.0)then
            avgtmp=tatlmn
         else
            print*,'Unusual condition with meanlattemp'
            istatus = 0
            return
         endif
         deallocate (rmeanlattemp)
         print*,'Domain average annual mean temp = ',avgtmp
      endif
c
c 2. green fraction.
c
      if(ctype.eq.'greenfrac')then
         avggrn=r_missing_data
         do l=1,ncat
            if(ic(l).gt.0)then
               avggrn(l)=sum(l)/float(ic(l))
            else
               avggrn(l)=0.0
            endif
            print*,'Domain average greenfrac = ',l,avggrn(l)
         enddo
c
c 3. Terrain slope index or soil type category..
c
      elseif(ctype.eq.'islope'.or.ctype.eq.'soiltype')then
        avgcat=r_missing_data
        if(islp.gt.0)then
           avgcat=float(isum)/islp
           print*,'Average category: ',avgcat
        else
           print*,'Could not compute average category ?'
           print*,'Average category: ',avgcat
        endif
      endif

c extend search to a fraction of the domain size. Could improve this for
c ratio geog-data-res/domain-res (possibly) to avoid unreasonable
c search distance.

      ijsthresh = int(nnxp/2)
      ifixw=0
      ifixl=0

      if(ctype.eq.'soiltemp')then

         do j = 1,nnyp
         do i = 1,nnxp

            if(landmask(i,j).eq.1)then                   !a land point
               if(geog_tmp(i,j,1).eq.r_missing_data)then !inconsistent because it is water
 
                  is=1
                  js=1
                  endsearch = .false.

                  sumt=0.0
                  isc=0

                  do while (.not.endsearch)

                   do jj=j-js,j+js
                   do ii=i-is,i+is

                      if((ii.ge.1) .and. (ii.le.nnxp) .and.
     &                   (jj.ge.1) .and. (jj.le.nnyp)) then

                       if(geog_tmp(ii,jj,1).ne.r_missing_data)then
                          sumt=sumt+geog_tmp(ii,jj,1)
                          isc=isc+1
                       endif

                      endif

                   enddo
                   enddo

                   if(isc.gt.0)then
                      geog_data(i,j,1)=-0.0065*topt_out(i,j)+sumt/isc
                      ifixw=ifixw+1    !count the # of inconsistent water points fixed with rep land value
                      endsearch=.true.
                   else
                      is=is+1
                      js=js+1
                      if(is.gt.ijsthresh)endsearch=.true.
                   endif

                  enddo

               else  !no inconsistency

                  geog_data(i,j,1)=-0.0065*topt_out(i,j)
     &+geog_data(i,j,1)

               endif

            else     !landmask says this is a water point

               if(geog_tmp(i,j,1).ne.r_missing_data)then
                  geog_data(i,j,1)=r_missing_data
                  ifixl=ifixl+1            !count the # of fixed land points
               endif

            endif

         enddo
         enddo

      elseif(ctype.eq.'greenfrac')then

         do j = 1,nnyp
         do i = 1,nnxp

            if(landmask(i,j).eq.1)then                  !a land point

               if(geog_tmp(i,j,1).eq.r_missing_data
     &.and.rmngrn.lt.r_missing_data)then  !inconsistency and there is valid data to search for.

                  endsearch = .false.

                  sum=0.0
                  ic=0
                  is=1
                  js=1
        
                  do while (.not.endsearch)

                   do ii=i-is,i+is
                   do jj=j-js,j+js

                    if( (ii.ge.1) .and. (ii.le.nnxp)
     &             .and.(jj.ge.1) .and. (jj.le.nnyp)) then

                       if(landmask(ii,jj).eq.1.and.
     &         geog_tmp(ii,jj,1).lt.r_missing_data)then

                          do l=1,12
                             sum(l)=sum(l)+geog_tmp(ii,jj,l)
                             ic(l)=ic(l)+1
                          enddo
                       endif

                    endif

                   enddo
                   enddo

                   if(ic(1).gt.0)then
                      do l=1,12
                         geog_data(i,j,l)=sum(l)/float(ic(l))
                      enddo
                      ifixw=ifixw+1
                      endsearch=.true.
                   else
                      is=is+1
                      js=js+1
                      if(is.gt.ijsthresh.or.js.gt.ijsthresh)
     &                   endsearch=.true.
                   endif

                  enddo

               else

                  geog_data(i,j,:)=geog_tmp(i,j,:)

               endif

            else  !this is a water point

c              do l=1,12
c                 if(geog_tmp(i,j,l).ne. 0.0 .and.
c    &               geog_tmp(i,j,l).lt.r_missing_data)then
c                    geog_data(i,j,l)=0.0
c                    ifixl=ifixl+1  !count # of land points that should be water
c                 endif
c              enddo

               if(geog_tmp(i,j,1).ne. 0.0 .and.
     &               geog_tmp(i,j,1).lt.r_missing_data)then
                     geog_data(i,j,1:12)=0.0
                     ifixl=ifixl+1  !count # of land points that should be water
               endif

            endif

         enddo
         enddo
c
c -------------------------------------------------------
c this section for categories (terrain slope or soiltype)
c
      elseif(ctype.eq.'islope'.or.ctype.eq.'soiltype')then

         do j = 1,nnyp
         do i = 1,nnxp

            if(landmask(i,j).eq.1)then                  !a land point

               if(geog_tmp(i,j,1).eq.r_missing_data)then !an inconsistency

                  endsearch = .false.

                  sumt=0.0
                  isc=0
                  islp=-9
                  is=1
                  js=1

                  do while (.not.endsearch)
                     do ii=i-is,i+is
                     do jj=j-js,j+js

                        if( (ii.ge.1) .and. (ii.le.nnxp)
     &                 .and.(jj.ge.1) .and. (jj.le.nnyp)) then

                           if(landmask(ii,jj).eq.1.and.
     &                        geog_tmp(ii,jj,1).lt.r_missing_data)then
                              islp=geog_tmp(ii,jj,1)
                              sumt=sumt+islp
                              isc=isc+1
                           endif

                        endif

                     enddo
                     enddo
 
                     if(ctype.eq.'islope')then
                        if(isc.gt.0)then
                           if(islp.gt.7.0)islp=13.0 !force mean to glacial ice as necessary
                           if(sumt.lt.isc)then
                              islp=1.0
                           else
                              islp=float(nint(sumt/float(isc)))
                           endif
                           geog_data(i,j,1)=islp
                           endsearch=.true.
                           ifixw=ifixw+1  !water point fixed to be land point
                        else
                           is=is+1
                           js=js+1
                           if(is.gt.ijsthresh.or.js.gt.ijsthresh)
     &                        endsearch=.true.
                        endif
                     else
                        if(isc.gt.0)then
                           if(sumt.lt.isc)then
                              islp=1.0
                           else
                              islp=float(nint(sumt/float(isc)))
                           endif
                           geog_data(i,j,1)=islp
                           endsearch=.true.
                           ifixw=ifixw+1  !water point fixed to be land point
                        else
                           is=is+1
                           js=js+1
                           if(is.gt.ijsthresh.or.js.gt.ijsthresh)
     &                        endsearch=.true.
                        endif
                     endif

                  enddo

               else

                  geog_data(i,j,1)=geog_tmp(i,j,1)

               endif

            else     !this is a water point

               if(geog_tmp(i,j,1).ne.r_missing_data)then
                  geog_data(i,j,1)=r_missing_data
                  if(ctype.eq.'islope')
     &            ifixl=ifixl+1   !land point fixed to be water point
               endif

            endif

         enddo
         enddo

      else

         print*,'Unknown type input to adjust_geog!'
         return

      endif ! ctype switch

      deallocate (geog_tmp)

c if the above search failed to find nearby soil temp or greenness frac
c then use average value

      if(ctype.eq.'soiltemp')then
         do j = 1,nnyp
         do i = 1,nnxp
            if(landmask(i,j).eq.1)then  !a land point
               if(geog_data(i,j,1).eq.r_missing_data)then
                  geog_data(i,j,1) = avgtmp-0.0065*topt_out(i,j)
                  ifixw=ifixw+1
               endif
            endif
         enddo
         enddo
         print*
         print*,'Soil Temp stats'

      elseif(ctype.eq.'greenfrac')then

         do j = 1,nnyp
         do i = 1,nnxp
            if(landmask(i,j).eq.1)then  !a land point
               do l = 1,12
                  if(geog_data(i,j,l).eq.0.0)then
                     geog_data(i,j,l)=float(nint(avggrn(l)))
                     if(l==1)ifixw=ifixw+1
                  endif
               enddo
            endif
         enddo
         enddo
         print*
         print*,'Green fraction stats'

      elseif(ctype.eq.'islope'.or.ctype.eq.'soiltype')then

         do j = 1,nnyp
         do i = 1,nnxp

            if(landmask(i,j).eq.1)then  !a land point
               if(geog_data(i,j,1).eq.r_missing_data)then
                  geog_data(i,j,1)=float(nint(avgcat))
                  if(ctype.eq.'islope')then
                     ifixw=ifixw+1
                  else
                     ifixl=ifixl+1
                  endif
               endif
            endif

        enddo
        enddo
        print*
        if(ctype.eq.'islope')then
           print*,'Terrain Slope Index stats'
        else
           print*,'Soil Type Category stats'
        endif

        where(geog_data.eq.r_missing_data)geog_data=0.0  !put water category back to original value

      else

        print*,'unknown geog data from ctype ',ctype

      endif

      print*,'--------------------------'
      if(ctype.ne.'islope'.and.ctype.ne.'soiltype')then
         if(ctype.eq.'greenfrac')then
            print*,'Fixed ',ifixw,'  points with rep land value'
            print*,'Fixed ',ifixl,'  points with missing value = water'
            print*
            do j = 1,nnyp
            do i = 1,nnxp
            do l = 1,ncat
               if(landmask(i,j) .ne. 0 .and.
     &geog_data(i,j,l).gt.0.0)then
                  sum(l)=geog_data(i,j,l)+sum(l)
                  ic(l)=ic(l)+1
               endif
            enddo
            enddo
            enddo
            do l=1,ncat
               if(ic(l).gt.0)then
                  avggrn(l)=sum(l)/float(ic(l))
               else
                  avggrn(l)=0.0
               endif
               print*,'Domain average greenfrac = ',l,avggrn(l)
            enddo
         else
            print*,'Fixed ',ifixw,'  points with rep land value'
            print*,'Fixed ',ifixl,'  points with missing value = water'
         endif
      elseif(ctype.eq.'islope')then
         print*,'Fixed ',ifixw,'  points with rep land value'
         print*,'Fixed ',ifixl,'  points with missing value = water'
      else
         print*,'Fixed ',ifixw,'  points with rep land value'
         print*,'Fixed ',ifixl,'  points with missing value = water'
      endif

      istatus=1

      return
      end

      subroutine filter_2dx(field,ix,iy,iz,smth)
c
c *** Subprogram:  smooth - Smooth a meteorological field.
c     Author:  Stan Benjamin 
c     Date  :  90-06-15
c
c *** Abstract:  Shapiro smoother. 
c 
c *** Program history log: 
c        85-12-09  S. Benjamin - Original version
c        96-06-16  J. Snook    - Modified to do 3d RAMS fields
c                              - hold array is dynamically allocated
c 
c *** Usage:  call smooth(field,ix,iy,iz,smth) 
c
c *** Input argument list: 
c        field    - real array  field(ix,iy,iz)
c                               Meteorological field
c        ix       - integer     x coordinates of field
c        iy       - integer     y coordinates of field
c        iz       - integer     z coordinates of field
c        smth     - real      
c
c *** Output argument list:   
c        field    - real array  field(ix,iy,iz)
c                               Smoothed meteorological field
c 
c *** Remarks:  Reference:  Shapiro, 1970: "Smoothing, filtering, and
c        boundary effects", Rev. Geophys. Sp. Phys., 359-387.
c
c     This filter is of the type 
c        z(i) = (1-s)z(i) + s(z(i+1)+z(i-1))/2
c     for a filter which is supposed to damp 2dx waves completely
c     but leave 4dx and longer with little damping,
c     it should be run with 2 passes using smth (or s) of 0.5
c     and -0.5.
c_______________________________________________________________________________
c   
      
c
      integer ix,iy,iz,i,j,k,i1,i2,it
c
      real field(ix,iy,iz),
     .     hold(ix,2),
     .     smth,smth1,smth2,smth3,smth4,smth5,
     .     sum1,sum2
c_______________________________________________________________________________
c
      smth1=0.25*smth*smth
      smth2=0.50*smth*(1.-smth)
      smth3=(1.-smth)*(1.-smth)
      smth4=(1.-smth)
      smth5=0.5*smth
c
      do k=1,iz
c
         do j=1,2
         do i=1,ix
            hold(i,j)=0.
         enddo
         enddo
c
         i1=2
         i2=1
         do j=2,iy-1
            it=i1
            i1=i2
            i2=it
            do i=2,ix-1
               sum1=field(i-1,j+1,k)+field(i-1,j-1,k)
     .             +field(i+1,j+1,k)+field(i+1,j-1,k)
               sum2=field(i  ,j+1,k)+field(i+1,j  ,k)
     .             +field(i  ,j-1,k)+field(i-1,j  ,k)
               hold(i,i1)=smth1*sum1+smth2*sum2+smth3*field(i,j,k)
            enddo
            if (j .eq. 2) goto 200
            do i=2,ix-1
               field(i,j-1,k)=hold(i,i2)
            enddo
200         continue
         enddo
c
         do i=2,ix-1
            field(i,iy-1,k)=hold(i,i1)
         enddo
c
         do i=2,ix-1
            field(i,1,k)=smth4*field(i,1,k) 
     .                  +smth5*(field(i-1,1,k)+field(i+1,1,k))
            field(i,iy,k)=smth4*field(i,iy,k) 
     .                   +smth5*(field(i-1,iy,k)+field(i+1,iy,k))
         enddo
c
         do j=2,iy-1
            field(1,j,k)=smth4*field(1,j,k) 
     .                  +smth5*(field(1,j-1,k)+field(1,j+1,k))
            field(ix,j,k)=smth4*field(ix,j,k) 
     .                   +smth5*(field(ix,j-1,k)+field(ix,j+1,k))
         enddo
c
      enddo
c
      return
      end

c
c ---------------------------------------------------------------
c
      subroutine get_meanlattemp(path_to_tiles,temp,istat)
c
      implicit      none

      integer       istat
      integer       n,i
      integer       ldir
      character*(*) path_to_tiles

      real          temp(180)

      istat=0
!      call s_len(path_to_tiles,ldir)
      ldir = LEN_TRIM(path_to_tiles)

      open(22,file=path_to_tiles(1:ldir)//'/LATMEANTEMP.DAT'
     &,form='formatted',status='old',iostat=istat)
      if(istat.ne.0) then
	write(6,*) 'insert bogus temp of 280'
	temp=280.0
	istat=1
	return
!	goto 3
      endif
      do i=1,180
         read(22,222,err=4)temp(i)
      enddo
 
      close(22)

      istat=1

      return

c     print*,'rmeantemp 1/2/3/4/5; ',rmeantemp(1),rmeantemp(45),rmeantemp(90)&
c          ,rmeantemp(135),rmeantemp(180)

222   format(1x,f6.2)

  3   print*,'Error: opening LATMEANTEMP file '
      print*,'path_to_tiles: ',path_to_tiles(1:ldir+3),ldir
      return
  4   print*,'Error: reading LATMEANTEMP file '

      return
      end

c
c===============================================================================
c
      subroutine init_hinterp(nx_bg,ny_bg,nx_laps,ny_laps,gproj,
     .     lat,lon,grx,gry,bgmodel,
     &     lat0,lon0,dlat,dlon,cgrddef)

! only part in subourine init_hinterp of WRFSI file lib/gridconv.f was copied.
!
      implicit none

      integer nx_bg,ny_bg,nx_laps,ny_laps,bgmodel

      real*4 lat(nx_laps,ny_laps),lon(nx_laps,ny_laps),
     .       grx(nx_laps,ny_laps),gry(nx_laps,ny_laps)

      character*1  cgrddef 
      character*2  gproj

      real   dlat, dlon
      real   lat0, lon0

      integer i,j,k
      integer istatus

      integer lenc

      integer nxc,nyc,nzc
      integer nx,ny
      logical lprintmessage
      real sw(2),ne(2),rota
      real nw(2),se(2),rlatc,rlonc
      real tolx,toly
      real grxdifsum1,grxdifsum2
      real grydifsum1,grydifsum2
      real r_missing_data
c________________________________________________________________________________
c
!      call get_r_missing_data(r_missing_data,istatus)

      r_missing_data = 1.e+37
      istatus = 1

      if(istatus.ne. 1)then
         print*,'Error getting r_missing_data - init_hinterp'
         return
      endif
c
c *** Determine location of LAPS grid point in background data i,j space.
c
      if (gproj .eq. 'LL') then
         call latlon_2_llij(nx_laps*ny_laps,lat,lon,grx,gry,
     &                      lat0,lon0,dlat,dlon,cgrddef)
      else
         print*,"Error: Unknown gproj spec in gridconv ",gproj
      endif
c
c *** Check that all LAPS grid points are within the background data coverage.
c

c set tolerance based upon the grid spacing as a function of grx/gry
      grxdifsum1=0.0
      grxdifsum2=0.0
      grydifsum1=0.0
      grydifsum2=0.0
      do j=1,ny_laps
         grxdifsum1=grxdifsum1+(grx(2,j)-grx(1,j))
         grxdifsum2=grxdifsum2+(grx(nx_laps-1,j)-grx(nx_laps,j))
      enddo
      do i=1,nx_laps
         grydifsum1=grydifsum1+(gry(i,2)-gry(i,1))
         grydifsum2=grydifsum2+(gry(i,ny_laps-1)-gry(i,ny_laps))
      enddo

      tolx=(abs(grxdifsum1)/ny_laps+abs(grxdifsum2)/ny_laps)*0.5
      toly=(abs(grydifsum1)/nx_laps+abs(grydifsum2)/nx_laps)*0.5

      print*,'horiz mapping tolerance x/y: ',tolx,toly
      lprintmessage=.true.
c
c *** First, check for wrapping if a global data set.
c
      if ( bgmodel .eq. 6 .or. 
     .     bgmodel .eq. 8) then
         do j=1,ny_laps
            do i=1,nx_laps
               if (grx(i,j) .lt. 1.) grx(i,j)=grx(i,j)+float(nx_bg)
               if (grx(i,j) .gt. nx_bg) grx(i,j)=grx(i,j)-float(nx_bg)
               if (gry(i,j) .lt. 1.) then
                  gry(i,j)=2.-gry(i,j)
                  grx(i,j)=grx(i,j)-float(nx_bg/2)
                  if (grx(i,j) .lt. 1.) grx(i,j)=grx(i,j)+float(nx_bg)
                  if (grx(i,j).gt.nx_bg) grx(i,j)=grx(i,j)-float(nx_bg)
               endif
               if (gry(i,j) .gt. ny_bg) then
                  gry(i,j)=float(2*ny_bg)-gry(i,j)
                  grx(i,j)=grx(i,j)-float(nx_bg/2)
                  if (grx(i,j) .lt. 1.) grx(i,j)=grx(i,j)+float(nx_bg)
                  if (grx(i,j).gt.nx_bg) grx(i,j)=grx(i,j)-float(nx_bg)
               endif
               if (grx(i,j) .lt. 1) then
                  grx(i,j)=grx(i,j)+1.
               endif
            enddo
         enddo
c
c ****** If not a global data set, then check that LAPS domain is fully
c           within background domain.
c
      else
         do j=1,ny_laps
            do i=1,nx_laps
c
c LAPS must fit into model grid which must also fit into LAPS grid thus we
c introduce a small fudge factor on the grid boundaries.
c               

               if(grx(i,j).gt.1.-tolx) grx(i,j) = max(1.,grx(i,j))
               if(gry(i,j).gt.1.-toly) gry(i,j) = max(1.,gry(i,j))

               if(grx(i,j).lt.float(nx_bg)+tolx) 
     +              grx(i,j) = min(float(nx_bg)-tolx,grx(i,j))
               if(gry(i,j).lt.float(ny_bg)+toly) 
     +              gry(i,j) = min(float(ny_bg)-toly,gry(i,j))

               if (grx(i,j) .lt. 1 .or. grx(i,j) .gt. nx_bg .or.
     .             gry(i,j) .lt. 1 .or. gry(i,j) .gt. ny_bg) then

                  grx(i,j) = r_missing_data
                  gry(i,j) = r_missing_data

           if(lprintmessage)then
              print*,'Domain gridpt outside of bkgd data coverage.'
              print*,'   data i,j,lat,lon - ',i,j,lat(i,j),lon(i,j)
              print*,'   grx, gry:',grx(i,j),gry(i,j)
              lprintmessage=.false.
c                 stop 'init_hinterp'
           endif

               endif
            enddo
         enddo
      endif
c
      return
      end
