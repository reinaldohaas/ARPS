c
c#######################################################################
c
c Postscript version of ZXplot routines using postscript character writing,
c Line pattern generating facilities.
c Written by Ming Xue at CIMMS/CAPS, Feb. 1990.
c
c#######################################################################
c
      SUBROUTINE xtpnup(x,y)
C position pen at point (x,y) defined in maths space
      common /xpen11/ xpen,ypen,flen,blen,npd,xmpen,ympen
      call ppenup( x,y)
      xpen= x
      ypen= y
      return
      END SUBROUTINE xtpnup

      SUBROUTINE xtpndn (x,y)
C Join point (x,y) defined in maths space with current line thickness
      common /xpen11/ xpen,ypen,flen,blen,npd,xmpen,ympen
      common /xlpn13/ hf1,hb1,hf2,hb2,lfull,lfull0,lthick, dthick
      xp1=xpen
      yp1=ypen
      xp2=x
      yp2=y
      call ppendn( xp2, yp2)
      xpen=xp2
      ypen=yp2
      return
      END SUBROUTINE xtpndn

      SUBROUTINE xbrokn(if1,ib1,if2,ib2)
C Set broken line patten in the one thousandth of the total vetical ND-space
C range unit.
      common /xlpn13/ hf1,hb1,hf2,hb2,lfull,lfull0,lthick, dthick
      common/psdef/io
      character char_io*132
      common /psscal/ p1,p2,p3,p4, xa,xb,ya,yb, xsca,ysca
      COMMON /XPHY01/PL,PR,PB,PT,XRANGE,YRANGE
      common /xfctr1/ fctr

      fctr = sqrt( abs(xrange*yrange) )
      h=0.0013*fctr
      hf1=h*if1
      hb1=h*ib1
      hf2=h*if2
      hb2=h*ib2
      p=p4-p3
      write(char_io,'(a)') 'S '
      call write_ps(char_io)
      write(char_io,100) max(1,nint(hf1*p)),max(1,nint(hb1*p)),
     :              max(1,nint(hf1*p)),max(1,nint(hb2*p))
      call write_ps(char_io)
 100  format(' [',4i3,' ] 0 d')
      lfull=0   ! Is the current line setting 'full'?
      lfull0=1  ! To use own dash line plotting algorithm? Used in XPENDN.

      entry xbrokn0
c
c To be called by XPSPAC to reset the dash line lengths when
c the size of plotting space as measured by xrange or yange change
c
      if( lfull.eq.1)  return
      fctr1 = sqrt( abs(xrange*yrange) )
      if( abs(fctr-fctr1).gt. 0.001) then
        p=p4-p3
        hf1=hf1*fctr1/fctr
        hb1=hb1*fctr1/fctr
        hf2=hf2*fctr1/fctr
        hb2=hb2*fctr1/fctr
        write(char_io,'(a)') 'S '
        call write_ps(char_io)
        write(char_io,100) max(1,nint(hf1*p)),max(1,nint(hb1*p)),
     :                max(1,nint(hf1*p)),max(1,nint(hb2*p))
        call write_ps(char_io)
        fctr = fctr1
      endif

      RETURN

      ENTRY xdash
C Set line atribute as dash line.
      fctr = sqrt( abs(xrange*yrange) )
      h=0.0013*fctr
      p =p4-p3
      hf1=h*10
      hb1=h*5
      hf2=h*10
      hb2=h*5
      write(char_io,'(a)') 'S '
      call write_ps(char_io)
      write(char_io,100) max(1,nint(hf1*p)),max(1,nint(hb1*p)),
     :              max(1,nint(hf1*p)),max(1,nint(hb2*p))
      call write_ps(char_io)
      lfull=0
      lfull0=1
      return

      ENTRY xdot
C Set line atribute as dash line.
      fctr = sqrt( abs(xrange*yrange) )
      h=0.0013*fctr
      p = p4-p3
      hf1=h*1
      hb1=h*6
      hf2=h*1
      hb2=h*6
      write(char_io,'(a)') 'S '
      call write_ps(char_io)
      write(char_io,100) max(1,nint(hf1*p)),max(1,nint(hb1*p)),
     :              max(1,nint(hf1*p)),max(1,nint(hb2*p))
      call write_ps(char_io)
      lfull=0
      lfull0=1
      return

      entry xqbrkn(kf1,kb1,kf2,kb2)
      h=0.0013*fctr
      kf1=hf1/h
      kb1=hb1/h
      kf2=hf2/h
      kb2=hb2/h
      return

      entry xfull
C Set line atribute as solid (full) line.
      lfull =1
      lfull0=1
      write(char_io,'(a)') 'S '
      call write_ps(char_io)
      write(char_io,'(a)') ' [] 0 d'
      call write_ps(char_io)
      return

      entry xqfull(kfull)
      kfull=lfull
      RETURN
      END SUBROUTINE xbrokn 

      SUBROUTINE xthick(ithick)
C Set thickness of lines.
      COMMON /XPHY01/PL,PR,PB,PT,XRANGE,YRANGE
      common /xlpn13/ hf1,hb1,hf2,hb2,lfull,lfull0,lthick, dthick
      integer lnmag,alnmag,blnmag
      save lnmag
      data lnmag /1/

      fctr = sqrt( abs(xrange*yrange) )
c     print*,'in xthick, lnmag,ithick,fctr=',lnmag,ithick,fctr

      lthick=ithick

c     print*,'ithick,lnmag,fctr=',ithick,lnmag,fctr

      call PSlnwd(0.25*ithick*lnmag*fctr)
      return

      entry xlnmag(alnmag)
c
c     Magnify the line thickness by a factor of nlnmag.
c
      lnmag = alnmag

      fctr = sqrt( abs(xrange*yrange) )
c     print*,'in xlnmag, lnmag,ithick,fctr=',lnmag,lthick,fctr
      call PSlnwd(0.5*lthick*lnmag*fctr)

      return

      entry xqlnmag(blnmag)
C Enquiry routine for line magnification factor
      blnmag = lnmag
      return

      entry xqthik(kthick)
C Enquiry routine for line thickness.
      kthick=lthick
      return
      END SUBROUTINE xthick

      SUBROUTINE Xcharl(x,y,ch)
      character*(*) ch
      COMMON /XCHP21/ XCHPEN, YCHPEN ,XCHMO,YCHMO,XCHPO,YCHPO
      integer wrtch
c
c#######################################################################
c
c Save the postion of last character string plotting. Strictly it should be
c the last pen postion at the end of string plotting, which is not calculated.
c
c#######################################################################
c
      call xtstchwrt(x,y,wrtch)
      if(wrtch.eq.0) return

      xchpen=x
      ychpen=y
      x1=x
      y1=y
      CALL xtrans(x1,y1)
      CALL PSstrg(x1,y1,ch,-1)
      RETURN
      END SUBROUTINE Xcharl
c
      SUBROUTINE Xcharr(x,y,ch)
      character*(*) ch
      COMMON /XCHP21/ XCHPEN, YCHPEN ,XCHMO,YCHMO,XCHPO,YCHPO
      integer wrtch

      call xtstchwrt(x,y,wrtch)
      if(wrtch.eq.0) return

      xchpen=x
      ychpen=y
      x1=x
      y1=y
      CALL xtrans(x1,y1)
      CALL PSstrg(x1,y1,ch,+1)
      RETURN
      END SUBROUTINE Xcharr

      SUBROUTINE Xcharc(x,y,ch)
      character*(*) ch
      COMMON /XCHP21/ XCHPEN, YCHPEN ,XCHMO,YCHMO,XCHPO,YCHPO
      integer wrtch

      call xtstchwrt(x,y,wrtch)
      if(wrtch.eq.0) return

      xchpen=x
      ychpen=y
      x1=x
      y1=y
      CALL xtrans(x1,y1)
      CALL PSstrg(x1,y1,ch,0)
      RETURN
      END SUBROUTINE Xcharc

      SUBROUTINE xchmag(h)
      common /xpsize/ psize1
      common /xcha20/ hctr,sctr,cratio, kfont,nundln
      common /xftr06/ xfactr,yfactr
      common /psscal/ p1,p2,p3,p4, xa,xb,ya,yb, xsca,ysca
      common /xpsd01/ xside, yside

      sctr=abs(h)/yfactr
      hctr = h

      if(abs(yside-1.0).lt.0.001) then ! yside=1.0
        fntsiz= abs(h)*psize1*550 ! when yside =1.0
      else
        fntsiz= abs(h)*psize1*460 ! when yside =1.5
      endif

      call PSftsz(nint(fntsiz*1.2))

      return

      ENTRY xqchmg(hh )
      hh=hctr
      return

      ENTRY xchsiz( hc )

      sctr= abs(hc)
      hctr = hc*yfactr

      if(abs(yside-1.0).lt.0.001) then ! yside=1.0
        fntsiz= abs(hc)*yfactr*psize1*550 ! when yside =1.0
      else
        fntsiz= abs(hc)*yfactr*psize1*460 ! when yside =1.5
      endif
      call PSftsz(nint(fntsiz*1.2))

      RETURN

      ENTRY xqchsz( cs1 )
      cs1= hctr/yfactr
      RETURN
      END SUBROUTINE xchmag

      SUBROUTINE xcfont( ifont )
      common /xcha20/ hctr,sctr,cratio, kfont,nundln

      if( ifont.eq. kfont) return
      call PSfont(ifont)
      kfont=ifont
      RETURN

      ENTRY xqcfnt ( nfont )
      nfont=kfont
      RETURN
      END SUBROUTINE xcfont

      SUBROUTINE xarrow(u,v,x0,y0, xleng,uunit)
C
C Plot vector (U,V) at (X0,Y0). by making use of PostScript
C arrow procedure. Unit X-component UUNIT is plotted with
C length XLENG in mapped maths coordinate. Length of vector in
C Y-direction is scaled according to the mapping.
C
c     parameter(pi=3.14159,angle1=(15./180+1)*pi,angle2=(-15./180+1)*pi)
c     PARAMETER(SINA1=-.25882,COSA1=-.96593,SINA2=-SINA1,COSA2=COSA1)

      parameter(pi=3.14159,angle1=(20./180+1)*pi,angle2=(-20./180+1)*pi)
      PARAMETER(SINA1=-.342,COSA1=-.93969,SINA2=-SINA1,COSA2=COSA1)

      common /xart36/ kartyp,kvmode,vsc
      
      REAL :: vlenfactor
      COMMON /VFACTOR/ vlenfactor
            
      common/psdef/io
      character char_io*132
c
      if( abs(u)+abs(v).eq.0.0) return
      alpha= xleng/uunit *0.5
      xc0=0.0
      yc0=0.0
      xpleng =xpntsd( xc0,yc0,xc0+xleng,yc0 )
      dx=u*alpha

      IF( KARTYP==3 .or. KARTYP==4 ) then
        call XQPSPC(pxl, pxr, pyb, pyt)
        call XQMAP (xl,xr,yb,yt)
        DY=V*ALPHA* (yt-yb)/(xr-xl) * (pxr-pxl)/(pyt-pyb)
      ELSE
        DY=V*ALPHA*vlenfactor
      ENDIF

c to plot arrow in absolute space (for conformality)
      px0=x0
      py0=y0
      px1=x0+dx
      py1=y0+dy
      call xtrans(px0,py0)
      call xtrans(px1,py1)
      dpx=px1-px0
      dpy=py1-py0
      pxa=px0-dpx
      pya=py0-dpy

      dpxy=sqrt( dpx*dpx+dpy*dpy)
      if(dpxy.gt.1.0e-30) then
        sinta=dpy/dpxy
        costa=dpx/dpxy
        arrow=0.40* min(xpleng,2*dpxy)
        if( kartyp.eq.2 ) arrow=0.40* dpxy*2
        px3=px1+arrow*(costa*cosa1-sinta*sina1)
        py3=py1+arrow*(sinta*cosa1+costa*sina1)
        px4=px1+arrow*(costa*cosa2-sinta*sina2)
        py4=py1+arrow*(sinta*cosa2+costa*sina2)
	call PStran(pxa,pya)
	call PStran(px1,py1)
	call PStran(px3,py3)
	call PStran(px4,py4)
C       write(char_io,100) pxa,pya,px1,py1,px3,py3,px4,py4
C       call write_ps(char_io)
C Or without using arw procedure:
C100  format(1x,8f9.4,' arw')
        write(char_io,200) pxa,pya,px1,py1,px3
        call write_ps(char_io)
        write(char_io,201) py3,px4,py4,px1,py1
        call write_ps(char_io)
      endif
200   format(1x,'S',2f9.4,' m',2f9.4,' l',f9.4)
201   format(f9.4,' l S',2f9.4,' m',2f9.4,' l S')
      return

      entry xartyp(ktype)
      kartyp=ktype
      RETURN
      END SUBROUTINE xarrow

      SUBROUTINE xvectk_old(x0,y0, xleng, uunit, key )
C Plot unit vectors  starting at (X0,Y0), by making use of PostScript
C arrow procedure.
C KEY=-1, 0, 1, 2, for none,in both X and Y-direction,X only, Y only
c
c Corrected an error with the vector key plot.
c
      parameter(pi=3.14159,angle1=(10./180+1)*pi,angle2=(-10./180+1)*pi)
      parameter(sina1=-.17365,cosa1=-.98481,sina2=-sina1,cosa2=cosa1)
      character ch*20
      common/psdef/io
      character char_io*132
      call xqrang(xrg,yrg)
      call xqmap(xl,xr,yb,yt)
      xscale=xr-xl
      yscale=yt-yb
      yf=0.4+( min(xrg,yrg)    -0.4)*0.5
      ctrsiz=3.0*yf*0.01
      call xchmag(ctrsiz)

      vunit=uunit
      dx=xleng
      dy=xleng
      pxo=x0
      pyo=y0
      px1=x0+dx
      py1=y0
      px2=x0
      py2=y0+dy
      call xtrans(pxo,pyo)
      call xtrans(px1,py1)
      call xtrans(px2,py2)
      dph=sqrt( (px1-pxo)**2+(py1-pyo)**2  )
      dpv=sqrt( (px2-pxo)**2+(py2-pyo)**2  )
 5    if( dpv.gt.1.5*dph ) then
          dpv=dpv*0.5
          dy=dy*0.5
          vunit=vunit*0.5
          goto 5
      endif
 6    if( dpv.lt.0.75*dph ) then
          dpv=dpv*2
          dy=dy*2
          vunit=vunit*2
          goto 6
      endif
      px2=x0
      py2=y0+dy
      call xtrans(px2,py2)

      if( key.eq.0.or.key.eq.1) then
        arrow=0.30*dph
        costa=(px1-pxo)/dph
        sinta=(py1-pyo)/dph
      endif
      if(key.eq.0.or.key.eq.2)  then
        arrow=0.30*dph
        costa=(px2-pxo)/dpv
        sinta=(py2-pyo)/dpv
      endif
        px3=px1+arrow*(costa*cosa1-sinta*sina1)
        py3=py1+arrow*(sinta*cosa1+costa*sina1)
        px4=px1+arrow*(costa*cosa2-sinta*sina2)
        py4=py1+arrow*(sinta*cosa2+costa*sina2)
	call PStran(pxo,pyo)
	call PStran(px1,py1)
	call PStran(px3,py3)
	call PStran(px4,py4)
C       write(char_io,100) pxo,pyo,px1,py1,px3,py3,px4,py4
c       call write_ps(char_io)
C Or without using arw procedure:
      if(  key.eq.0.or.key.eq.1) then
        write(char_io,200) pxo,pyo,px1,py1,px3
        call write_ps(char_io)
        write(char_io,201) py3,px4,py4,px1,py1
        call write_ps(char_io)
        write(ch,'(f6.1,'' m/s'')') uunit
        lch=10
        call xchlj (ch,lch)
        call xcharl(x0+dx   +0.01*xscale,y0,ch(1:lch) )
      endif
      if(key.eq.0.or.key.eq.2)  then
        write(char_io,200) pxo,pyo,px2,py2,px3
        call write_ps(char_io)
        write(char_io,201) py3,px4,py4,px2,py2
        call write_ps(char_io)
        call xqobag( xang, yang )
        call xqchor( asym )
        call xchori(90.0+ yang- xang)
        write(ch,'(f6.1,'' m/s'')') vunit
        lch=10
        call xchlj (ch,lch)
        call xcharl(x0-.02*xscale ,y0 ,ch(1:lch) )
        call xchori( asym )
      endif
 100  format(1x,8f7.2,' arw')
 200  format(1x,'S',2f7.2,' m',2f7.2,' l',f7.2)
 201  format(f7.2,' l S',2f7.2,' m',2f7.2,' l S')
      RETURN
      END SUBROUTINE xvectk_old
      
