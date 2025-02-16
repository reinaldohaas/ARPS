      SUBROUTINE XTPNUP(X,Y)
C position pen at point (x,y) defined in maths space
      COMMON /XPEN11/ XPEN,YPEN,FLEN,BLEN,NPD,XMPEN,YMPEN
      CALL PPENUP( X,Y)
      XPEN= X
      YPEN= Y
      RETURN
      END
      SUBROUTINE XTPNDN (X,Y)
C Join point (x,y) defined in maths space with current line thickness
      COMMON /XPEN11/ XPEN,YPEN,FLEN,BLEN,NPD,XMPEN,YMPEN
      COMMON /XLPN13/ HF1,HB1,HF2,HB2,LFULL,lfull0,LTHICK,DTHICK
      XP1=XPEN
      YP1=YPEN
      XP2=X
      YP2=Y
c     IF( LTHICK.EQ.1) THEN
        CALL PPENDN( XP2, YP2)
c     ELSE
c       DX=XP2-XP1
c       DY=YP2-YP1
c       IF( ABS(DX).LT.1.0E-6) THEN
c         DPX=DTHICK
c         DPY=0.0
c       ELSE
c         AK=DY/DX
c         A=SQRT(1.0+AK*AK)
c         DPX=DTHICK*AK/A
c         DPY=DTHICK/A
c       ENDIF
c       CALL PPENDN( XP2,YP2)
c       CALL PPENUP( XP2+DPX, YP2-DPY)
c       CALL PPENDN( XP1+DPX, YP1-DPY)
c       CALL PPENUP( XP1-DPX, YP1+DPY)
c       CALL PPENDN( XP2-DPX, YP2+DPY)
c       CALL PPENUP( XP2,YP2)
c     ENDIF
      XPEN=XP2
      YPEN=YP2
      RETURN
      END

      SUBROUTINE XBROKN(IF1,IB1,IF2,IB2)
C Set line atribute as broken line segments.
C IF1 IF2 set length of  line segments plotted in unit of 0.001
C that of the total vertical ND-space range.
C IB1 IB2 set length of blanks between line segments in unit of 0.001
C that of the total vertical ND-space range.
      COMMON /XLPN13/ HF1,HB1,HF2,HB2,LFULL,lfull0,LTHICK,DTHICK
      H=0.001
      HF1=H*IF1
      HB1=H*IB1
      HF2=H*IF2
      HB2=H*IB2
      lfull=0   ! Is the current line setting 'full'?
      lfull0=0  ! To use own dash line plotting algorithm? Used in XPENDN.
                ! Set it to lfull if want to use own dash line plotting algorithm
                !   and remove call to gsln
c     call gsln(2)
      RETURN

      ENTRY XBROKN0
      RETURN

      ENTRY XDASH
C Set line atribute as dash line.
      H=0.001
      HF1=H*10
      HB1=H*5
      HF2=H*10
      HB2=H*5
      LFULL=0
      LFULL0=0
c     call gsln(2)
      RETURN

      ENTRY XDOT
C Set line atribute as dash line.
      H=0.001
      HF1=0.0
      HB1=H*6
      HF2=0.0
      HB2=H*6
      LFULL=0
      LFULL0=0
c     call gsln(3)
      RETURN

      ENTRY XQBRKN(KF1,KB1,KF2,KB2)
      H=0.001
      KF1=HF1/H
      KB1=HB1/H
      KF2=HF2/H
      KB2=HB2/H
      RETURN

      ENTRY XFULL

C Set line atribute as solid (full) line.
      LFULL =1
      LFULL0=1
c     call gsln(1)
      RETURN

      ENTRY XQFULL(KFULL)
      KFULL=LFULL
      RETURN
      END

      SUBROUTINE XTHICK(ITHICK)
      COMMON /XLPN13/ HF1,HB1,HF2,HB2,LFULL,lfull0,LTHICK, DTHICK
      integer alnmag,blnmag

C Set thickness of lines. ITHICK=1 or 2.
c     LTHICK=min(ITHICK,2)
      LTHICK=ITHICK
      call plotif(0.0, 0.0, 2)
      call gslwsc(1.0+0.5*(ithick-1))
      RETURN

      entry xlnmag(alnmag)
c
c     Magnify the line thickness by a factor of nlnmag.
c
      return

      entry xqlnmag(blnmag)

      return

      entry xqthik(kthick)
C Enquiry routine for line thickness.
      kthick=lthick
      return
      END

      SUBROUTINE XCHARL(XO,YO,STRING)
      CHARACTER*(*) STRING
      call xqthik(lthick)
      call plotif(0.0, 0.0, 2)
      call gslwsc(1.0)
      call gqln(ierr,lstyle)
      call gsln(1)
      CALL XLETER(XO,YO,STRING, -1)
      call plotif(0.0, 0.0, 2)
      call gslwsc(1.0+0.5*(lthick-1))
      if (ierr.eq.0) call gsln(lstyle)
      RETURN
      END

      SUBROUTINE XCHARR(XO,YO,STRING)
      CHARACTER*(*) STRING
      call xqthik(lthick)
      call plotif(0.0, 0.0, 2)
      call gslwsc(1.0)
      call gqln(ierr,lstyle)
      call gsln(1)
      CALL XLETER(XO,YO,STRING,  1)
      call plotif(0.0, 0.0, 2)
      call gslwsc(1.0+0.5*(lthick-1))
      if (ierr.eq.0) call gsln(lstyle)
      RETURN
      END

      SUBROUTINE XCHARC(XO,YO,STRING)
      CHARACTER*(*) STRING
      call xqthik(lthick)
      call plotif(0.0, 0.0, 2)
      call gslwsc(1.0)
      call gqln(ierr,lstyle)
      call gsln(1)
      CALL XLETER(XO,YO,STRING,  0)
      call plotif(0.0, 0.0, 2)
      call gslwsc(1.0+0.5*(lthick-1))
      if (ierr.eq.0) call gsln(lstyle)
      RETURN
      END

      SUBROUTINE XCHMAG(H)
      COMMON /XCHA20/ HCTR,SCTR,CRATIO, KFONT,NUNDLN
      COMMON /XFTR06/ XFACTR,YFACTR
      HCTR=ABS(H)
      SCTR=HCTR/YFACTR
      RETURN

      ENTRY XQCHMG(HH )
      HH=HCTR
      RETURN

      ENTRY XCHSIZ( HC )
      SCTR= abs(HC)
      HCTR= SCTR*YFACTR
      RETURN

      ENTRY XQCHSZ( CS1 )
      CS1= HCTR/YFACTR
      RETURN
      END

      SUBROUTINE XCFONT( IFONT )
      COMMON /XCHA20/ HCTR,SCTR,CRATIO, KFONT,NUNDLN
      COMMON /XCHR31/ CHDATA
      COMMON /XCHR32/ ICDATA
      CHARACTER CHDATA(127)*300
      common /xoutch/ nch
      INTEGER ICDATA (0:150 , 32:127)
      KF=KFONT
      KFONT=IFONT
      IF ((KFONT .LT. 1) .OR. (KFONT .GT. 4)) THEN
         WRITE(NCH,'(1x,a,i2,a,/1x,a)')'Font number ',ifont,
     :   ' not available with the NCAR graphics version of ZXPLOT. ',
     :   'Default font number 2 will be used.'
         kfont = 2
      ENDIF
      IF( KF.EQ. KFONT) RETURN
      GOTO ( 501, 502, 503, 504 ) KFONT
 501  CALL XCSETB(CHDATA)
      GOTO 505
 502  CALL XCSETC(CHDATA)
      GOTO 505
 503  CALL XCSETA(CHDATA)
      GOTO 505
 504  CALL XCSETD(CHDATA)
 505  CONTINUE
      RETURN

      ENTRY XQCFNT ( NFONT )
      NFONT=KFONT
      RETURN
      END

      SUBROUTINE XARROW(U,V,X0,Y0, XLENG,UUNIT)
C Plot vector (U,V) at (X0,Y0). Unit X-component UUNIT is plotted with
C length XLENG in mapped maths coordinate. Length of vector in
C Y-direction is scaled according to the mapping.
c     PARAMETER(PI=3.14159,ANGLE1=(15./180+1)*PI,ANGLE2=(-15./180+1)*PI)
c     PARAMETER(SINA1=-.25882,COSA1=-.96593,SINA2=-SINA1,COSA2=COSA1)

      parameter(pi=3.14159,angle1=(20./180+1)*pi,angle2=(-20./180+1)*pi)
      PARAMETER(SINA1=-.342,COSA1=-.93969,SINA2=-SINA1,COSA2=COSA1)

      COMMON /XART36/ KARTYP,KVMODE,VSC

      REAL :: vlenfactor
      COMMON /VFACTOR/ vlenfactor
C
      IF( ABS(U)+ABS(V).EQ.0.0) RETURN
      ALPHA= XLENG/UUNIT *0.5
      XC0=0.0
      YC0=0.0
      XPLENG =XPNTSD( XC0,YC0,XC0+XLENG,YC0 )

      DX=U*ALPHA

      IF( KARTYP==3 .or. KARTYP==4 ) then
        call XQPSPC(pxl, pxr, pyb, pyt)
        call XQMAP (xl,xr,yb,yt)
        DY=V*ALPHA* (yt-yb)/(xr-xl) * (pxr-pxl)/(pyt-pyb)
      ELSE
        DY=V*ALPHA*vlenfactor
      ENDIF

C TO PLOT ARROW IN ABSOLUTE SPACE (FOR CONFORMALITY)
      PX0=X0
      PY0=Y0
      PX1=X0+DX
      PY1=Y0+DY
      CALL XTRANS(PX0,PY0)
      CALL XTRANS(PX1,PY1)
      DPX=PX1-PX0
      DPY=PY1-PY0

      DPXY=SQRT( DPX*DPX+DPY*DPY)
      IF(DPXY.GT.1.0E-30) THEN
        SINTA=DPY/DPXY
        COSTA=DPX/DPXY
        ARROW=0.40* MIN(XPLENG,2*DPXY)
        IF( KARTYP.EQ.2 .or. KARTYP.EQ.4 ) ARROW=0.40* DPXY*2
        DX1=ARROW*(COSTA*COSA1-SINTA*SINA1)
        DY1=ARROW*(SINTA*COSA1+COSTA*SINA1)
        DX2=ARROW*(COSTA*COSA2-SINTA*SINA2)
        DY2=ARROW*(SINTA*COSA2+COSTA*SINA2)
        CALL XPENUP(X0-DX, Y0-DY)
        CALL XPENDN(X0+DX, Y0+DY)
!       CALL XTPNUP(PX1    , PY1    )
!       CALL XTPNDN(PX1+DX1, PY1+DY1)
!       CALL XTPNUP(PX1    , PY1    )
!       CALL XTPNDN(PX1+DX2, PY1+DY2)

        CALL XTPNUP(PX1+DX1, PY1+DY1)
        CALL XTPNDN(PX1    , PY1    )
        CALL XTPNDN(PX1+DX2, PY1+DY2)
      ENDIF
      RETURN

      ENTRY XARTYP(KTYPE)
      KARTYP=KTYPE
      RETURN
      END SUBROUTINE xarrow

      SUBROUTINE XVECTK_old(X0,Y0, XLENG, UUNIT, KEY )
C Plot unit vectors  starting at (X0,Y0)
C KEY=-1, 0, 1, 2, for none,in both X and Y-direction,X only, Y only
      PARAMETER(PI=3.14159,ANGLE1=(10./180+1)*PI,ANGLE2=(-10./180+1)*PI)
      PARAMETER(SINA1=-.17365,COSA1=-.98481,SINA2=-SINA1,COSA2=COSA1)
      CHARACTER CH*20
      CALL XQRANG(XRG,YRG)
      CALL XQMAP(XL,XR,YB,YT)
      XSCALE=XR-XL
      YSCALE=YT-YB
      YF=0.4+( MIN(XRG,YRG)-0.4)*0.5
      CTRSIZ=3.0*YF*0.01
      CALL XCHMAG(CTRSIZ)

      VUNIT=UUNIT
      DX=XLENG
      DY=XLENG
      PXO=X0
      PYO=Y0
      PX1=X0+DX
      PY1=Y0
      PX2=X0
      PY2=Y0+DY
      CALL XTRANS(PXO,PYO)
      CALL XTRANS(PX1,PY1)
      CALL XTRANS(PX2,PY2)
      DPH=SQRT( (PX1-PXO)**2+(PY1-PYO)**2  )
      DPV=SQRT( (PX2-PXO)**2+(PY2-PYO)**2  )
 5    IF( DPV.GT.1.5*DPH ) THEN
          DPV=DPV*0.5
          DY=DY*0.5
          VUNIT=VUNIT*0.5
          GOTO 5
      ENDIF
 6    IF( DPV.LT.0.75*DPH ) THEN
          DPV=DPV*2
          DY=DY*2
          VUNIT=VUNIT*2
          GOTO 6
      ENDIF
      PX2=X0
      PY2=Y0+DY
      CALL XTRANS(PX2,PY2)

      IF( KEY.EQ.0.OR.KEY.EQ.1) THEN
        ARROW=0.30*DPH
        COSTA=(PX1-PXO)/DPH
        SINTA=(PY1-PYO)/DPH
        DX1=ARROW*(COSTA*COSA1-SINTA*SINA1)
        DY1=ARROW*(SINTA*COSA1+COSTA*SINA1)
        DX2=ARROW*(COSTA*COSA2-SINTA*SINA2)
        DY2=ARROW*(SINTA*COSA2+COSTA*SINA2)
        CALL XTPNUP(PXO,PYO)
        CALL XTPNDN(PX1    , PY1)
C       CALL XTPNUP(PX1    , PY1    )
        CALL XTPNDN(PX1+DX1, PY1+DY1)
        CALL XTPNUP(PX1    , PY1    )
        CALL XTPNDN(PX1+DX2, PY1+DY2)
      WRITE(CH,'(F6.1,'' m/s'')') UUNIT
      LCH=10
      CALL XCHLJ (CH,LCH)
      CALL XCHARL(X0+DX   +0.01*XSCALE,Y0,CH(1:LCH) )
      ENDIF
      IF(KEY.EQ.0.OR.KEY.EQ.2)  THEN
        ARROW=0.20*DPH
        COSTA=(PX2-PXO)/DPV
        SINTA=(PY2-PYO)/DPV
        DX1=ARROW*(COSTA*COSA1-SINTA*SINA1)
        DY1=ARROW*(SINTA*COSA1+COSTA*SINA1)
        DX2=ARROW*(COSTA*COSA2-SINTA*SINA2)
        DY2=ARROW*(SINTA*COSA2+COSTA*SINA2)
        CALL XTPNUP(PXO,PYO)
        CALL XTPNDN(PX2    , PY2)
C       CALL XTPNUP(PX2    , PY2    )
        CALL XTPNDN(PX2+DX1, PY2+DY1)
        CALL XTPNUP(PX2    , PY2    )
        CALL XTPNDN(PX2+DX2, PY2+DY2)
        CALL XQOBAG( XANG, YANG )
        CALL XQCHOR( ASYM )
        CALL XCHORI(90.0+ YANG- XANG)
        WRITE(CH,'(F6.1,'' m/s'')') VUNIT
        LCH=10
        CALL XCHLJ (CH,LCH)
        CALL XCHARL(X0-.02*XSCALE ,Y0 ,CH(1:LCH) )
        CALL XCHORI( ASYM )
      ENDIF
      RETURN
      END SUBROUTINE XVECTK_old
