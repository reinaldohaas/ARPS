c
c
c
c
c                   ########################################
c                   ########################################
c                   ########################################
c                   ########                        ########
c                   ########        XSCATTER        ########
c                   ########                        ########
c                   ########################################
c                   ########################################
c                   ########################################
c
c
      Subroutine Xscatter (x, y, n, itype, diameter)
      Implicit None
      Integer	n, itype, i
      Real	x(n), y(n), dx, dy
      Real	pi,size,deg2rad,radius,diameter, c22,s22,
     >  xx,x1,x2,x3,x4, yy,y1,y2,y3,y4
      Parameter (pi=3.14159265, deg2rad=pi/180.)
c
      radius	= .5 * diameter
      Call Xfull
      Call Plt2Wrld (0.,0.,x1,y1)
      Call Plt2Wrld (radius,radius,x2,y2)
      dx	= Abs (x2 - x1)
      dy	= Abs (y2 - y1)
c  
c...dot
c
      If (itype.EQ.0) Then
c       Call Xpntsz (0.002)
        Do i=1,n
        Call Xpoint (x(i), y(i))
        End Do
c       Call Xpntsz (0.0005)
c  
c...X
c
      Else If (itype.EQ.1) Then
        Do i=1,n
        Call Xpenup (x(i)+dx, y(i)+dy)
        Call Xpendn (x(i)-dx, y(i)-dy)
        Call Xpenup (x(i)-dx, y(i)+dy)
        Call Xpendn (x(i)+dx, y(i)-dy)
        End Do
c  
c...crosshair
c
      Else If (itype.EQ.2) Then
        Do i=1,n
        Call Xpenup (x(i)+dx, y(i))
        Call Xpendn (x(i)-dx, y(i))
        Call Xpenup (x(i), y(i)+dy)
        Call Xpendn (x(i), y(i)-dy)
        End Do
c  
c...triangle
c
      Else If (itype.EQ.3) Then
	xx	= dx * Cos(30.*deg2rad)
	yy	= dy * Cos(60.*deg2rad)
        Do i=1,n
        Call Xpenup (x(i), y(i)+dy)
        Call Xpendn (x(i)+xx, y(i)-yy)
        Call Xpendn (x(i)-xx, y(i)-yy)
        Call Xpendn (x(i), y(i)+dy)
        End Do
c  
c...square
c
      Else If (itype.EQ.4) Then
        Do i=1,n
        Call Xpenup (x(i)+dx, y(i)+dy)
        Call Xpendn (x(i)+dx, y(i)-dy)
        Call Xpendn (x(i)-dx, y(i)-dy)
        Call Xpendn (x(i)-dx, y(i)+dy)
        Call Xpendn (x(i)+dx, y(i)+dy)
        End Do
c  
c...diamond
c
      Else If (itype.EQ.5) Then
        Do i=1,n
        Call Xpenup (x(i)+dx, y(i))
        Call Xpendn (x(i), y(i)-dy)
        Call Xpendn (x(i)-dx, y(i))
        Call Xpendn (x(i), y(i)+dy)
        Call Xpendn (x(i)+dx, y(i))
        End Do
c  
c...octagon
c
      Else If (itype.EQ.8) Then
        c22	= Cos (22.5*deg2rad)
        s22	= Sin (22.5*deg2rad)
        xx	= x2 - x1
        yy	= y2 - y1
c
        Do i=1,n
          x4	= x(i) + xx * c22
          x3	= x(i) + xx * s22
          x2	= x(i) - xx * s22
          x1	= x(i) - xx * c22
          y4	= y(i) + yy * c22
          y3	= y(i) + yy * s22
          y2	= y(i) - yy * s22
          y1	= y(i) - yy * c22
c
          Call Xpenup (x3,y4)
          Call Xpendn (x4,y3)
          Call Xpendn (x4,y2)
          Call Xpendn (x3,y1)
          Call Xpendn (x2,y1)
          Call Xpendn (x1,y2)
          Call Xpendn (x1,y3)
          Call Xpendn (x2,y4)
          Call Xpendn (x3,y4)
        End Do
c  
      Else 
	Print *, 'Xscatter: no such line type: ', itype
      End If
c  
      End
c
c
c
c                   ########################################
c                   ########################################
c                   ########################################
c                   ########                        ########
c                   ########        PLT2WRLD        ########
c                   ########        WRLD2PLT        ########
c                   ########                        ########
c                   ########################################
c                   ########################################
c                   ########################################
c
c
      Subroutine Plt2wrld (px,py,x,y)
      Implicit None
      Real x,y,px,py
      x = px
      y = py
      Call Xlinvt (x,y)
      End
c
c
      Subroutine Wrld2Plt (x,y,px,py)
      Implicit None
      Real x,y,px,py
      px        = x
      py        = y
      Call Xtrans (px,py)       !input:x,y in world coords. output: plot coords
      End

c
c
c
c
c                   ########################################
c                   ########################################
c                   ########################################
c                   ########                        ########
c		    ########         BOXER 	    ########
c                   ########                        ########
c                   ########################################
c                   ########################################
c                   ########################################
c
c
      Subroutine Boxer (nbracket,xl,xr,yl,yr,dx,dy)
      Implicit None
      Integer	nbracket
      Real	xl,yl,xr,yr,dx,dy
c
      If (nbracket.EQ.1) Then
	Call Xpenup (xl,yl+dy)
	Call Xpendn (xl,yl)
	Call Xpendn (xl+dx,yl)
c
	Call Xpenup (xr-dx,yl)
	Call Xpendn (xr,yl)
	Call Xpendn (xr,yl+dy)
c
	Call Xpenup (xr,yr-dy)
	Call Xpendn (xr,yr)
	Call Xpendn (xr-dx,yr)
c
	Call Xpenup (xl+dx,yr)
	Call Xpendn (xl,yr)
	Call Xpendn (xl,yr-dy)
c
      Else
	Call Xbox   (xl,xr,yl,yr)
      End If
      End
c
c
c
c
c
c		    ########################################
c		    ########################################
c		    ########################################
c		    ########			    ########
c		    ########         MYLINE 	    ########
c		    ########			    ########
c		    ########################################
c		    ########################################
c		    ########################################
c
c
      Subroutine MyLine (x1,y1,x2,y2)
      Implicit None
      Real      x1,y1,x2,y2
      Call Xpenup (x1,y1)
      Call Xpendn (x2,y2)
      End
c
c
c
c
c                   ########################################
c                   ########################################
c                   ########################################
c                   ########                        ########
c                   ########         PCURVES        ########
c                   ########                        ########
c                   ########################################
c                   ########################################
c                   ########################################
c
c
      Subroutine Pcurves (px1,px2,py1,py2, xmin,xmax,ymin,ymax, ifclip,
     >  labelx,xlabel,fmtxax,xinc, labely,ylabel,fmtyax,yinc,
     >  n,ncurv,x1,y1,y2,y3,y4, l1,l2,l3,l4, lc1,lc2,lc3,lc4)
      Implicit None
c
c  labelx:
c    0: plot/label B, plot T
c    1: plot/label B
c    2: plot B
c    3: plot T
c    4: plot T,B
c
c  labely:
c    0: no y axis
c   -1: plot/label on L
c   -2: plot/label on R
c    1: also plot on R
c    2: also plot on L
c
c  xinc,yinc	= Axis tick intervals. Use 0 for default.
c
c
      Integer	n,ncurv, labelx,labely,ifclip,
     >  l1,l2,l3,l4, n1,n2, i,nn, lc1,lc2,lc3,lc4
      Real	px1,px2,py1,py2, pxc,pyc, 
     >  xmin,xmax,xinc, ymin,ymax,yinc, xx,yy, 
     >  x1(n), y1(n), y2(n), y3(n), y4(n) 
      Character*(*) xlabel,fmtxax,ylabel,fmtyax
c
      Call Xwdwof
      Call Xaxdef
      Call Xaxfmt ('*')
c     Call Xaxtik (1,1)
c     Call Xaxant (-1,-1)
      Call Xthick (1)
      Call Xfull
c
c     Print 9001, 'Pcurves: px1,px2,py1,py2: ', px1,px2,py1,py2
c     Print 9001, 'Pcurves: xmin,xmax,ymin,ymax: ', xmin,xmax,ymin,ymax
 9001 Format (1x,a,4g14.4)
      Call Xpspac (px1,px2,py1,py2)
      Call Xmap   (xmin,xmax,ymin,ymax)
c
c  x axis
c
      If ((labelx.GE.0 .AND. labelx.LE.2) .OR. labelx.EQ.4) Then
	Call Xaxdef
        Call Xaxant (0,1)
        Call Xaxfmt (fmtxax)
	If (labelx.LE.1) Then
          pxc	= px1 + .5 * (px2-px1)
	  pyc	= Max (0., py1-0.06)
          Call Xaxant (-1,1)
          Call Plt2Wrld (pxc,pyc,xx,yy)
          Call Xcharc (xx,yy,xlabel)
	End If
        Call Xaxisx (xmin, ymin, xinc)	!label bot
      End If
c
      If (labelx.EQ.0 .OR. labelx.EQ.3 .OR. labelx.EQ.4) Then
        Call Xaxtik (-1,1)
        Call Xaxant (0,1)
        Call Xaxisx (xmin, ymax, xinc)	!label top
      End If
c
c  y axis
c
      Call Xaxdef
c
c     yinc	= 0.
c     If (py2-py1.LE.0.15) Then
cIf (ymax-ymin.LE.5.01) Then
c  yinc = 1.
cElse If (ymax-ymin.LE.10.01) Then
c  yinc = 2.
c       End If
c     End If
c     If (yinc.NE.0.) Print *, '% Pcurves: yaxis inc: ', yinc
c
      If (Abs(labely).EQ.1) Then
        Call Xaxfmt (fmtyax)
        Call Xaxisy (xmin, ymin, yinc)
        pyc	= py1 + .5 * (py2-py1)
	pxc	= Max (0., px1-0.06)
        Call Plt2Wrld (px1-0.06,pyc,xx,yy)
        Call Xchori (90.)
        Call Xcharc (xx,yy,ylabel)
        Call Xchori (0.)
c
	If (labely.GT.0) Then
        Call Xaxtik (1,-1)
        Call Xaxant (1,0)
        Call Xaxisy (xmax, ymin, yinc)
	End If
      Else If (Abs(labely).EQ.2) Then
        pyc	= py1 + .5 * (py2-py1)
	pxc	= Min (1., px2+0.06)
        Call Xaxfmt (fmtyax)
        Call Xaxtik (1,-1)
        Call Xaxant (1,1)
        Call Xaxisy (xmax, ymin, yinc)
        pyc	= py1 + .5 * (py2-py1)
        Call Plt2Wrld (pxc,pyc,xx,yy)
        Call Xchori (-90.)
        Call Xcharc (xx,yy,ylabel)
        Call Xchori (0.)
c
	If (labely.GT.0) Then
        Call Xaxtik (1,1)
        Call Xaxant (1,0)
        Call Xaxisy (xmin, ymin, yinc)
	End If
      End If
c
c  clip
c  always clip the x axis
c
      If (ifclip.NE.0) Call Xwindw (xmin,xmax,ymin,ymax)
      n1	= 1
      n2	= n
      Do i=1,n
      If (xmin.LE.x1(i)) Then 
      n1	= i
      Go To 1001
      End If
      End Do
 1001 Continue
c
      Do i=n,1,-1
      If (xmax.GE.x1(i)) Then
      n2	= i
      Go To 1002
      End If
      End Do
 1002 Continue
c
      nn	= n2 - n1 + 1
      If (nn.NE.n) 
     >  Print '(1p,2(1x,a,i4,g12.4))', 
     >  'Pcurves: plotting points ', n1,x1(n1), 'through ', n2,x1(n2)
c
c  Plot curves
c
      If (ncurv.LT.1) Go To 1000
c     print *, '...curve 1'
      Call Lintyp (l1)
      Call Color (lc1)
      Call Xcurve (x1(n1), y1(n1), nn, 0)
c 
      If (ncurv.LT.2) Go To 1000
c     print *, '...curve 2'
      Call Lintyp (l2)
      Call Color (lc2)
      Call Xcurve (x1(n1), y2(n1), nn, 0)
c 
      If (ncurv.LT.3) Go To 1000
c     print *, '...curve 3'
      Call Lintyp (l3)
      Call Color (lc3)
      Call Xcurve (x1(n1), y3(n1), nn, 0)
c 
      If (ncurv.LT.4) Go To 1000
c     print *, '...curve 4'
      Call Lintyp (l4)
      Call Color (lc4)
      Call Xcurve (x1(n1), y4(n1), nn, 0)
c
 1000 Continue
      Call Xwdwof
      Call Color (0)
      Call Lintyp (0)
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
c                   ########         LINTYP         ########
c                   ########                        ########
c                   ########################################
c                   ########################################
c                   ########################################
c
c
      Subroutine Lintyp (l)
      Implicit None
      Integer l, ll
c
c
c  LINE TYPES:
c    (negative: bold)
c    0 default
c    1 full
c    2 dash
c    3 dot
c    4 my dot
c    5 my dash
c    6 my dash-dot
c
c
      If (l.GE.0) Then
	Call Xthick (1)
      Else
	Call Xthick (2)
      End If
c
      ll	= Abs (l)
      If (ll.LE.1) Then
        Call Xfull
      Else If (ll.EQ.2) Then
        Call Xdash
      Else If (ll.EQ.3) Then
        Call Xdot
      Else If (ll.EQ.4) Then
        Call Xbrokn (1,6,1,6)
      Else If (ll.EQ.5) Then
        Call Xbrokn (4,4,4,4)
      Else If (ll.EQ.6) Then
        Call Xbrokn (1,6,4,6)
      Else
	Call Xfull
      End If
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
c                   ########          INC           ########
c                   ########                        ########
c                   ########################################
c                   ########################################
c                   ########################################
c
c  Useful for setting max/min for plotting.
c
      Subroutine Inc (a, b, d, n)
c
      If (d.GT.0.) Then
        Do i=1,n
        If (b.GT.a) a = a + d
        End Do
      Else 
        Do i=1,n
        If (b.LT.a) a = a + d
        End Do
      End If
      End
c
c
c
c
c                   ########################################
c                   ########################################
c                   ########################################
c                   ########                        ########
c                   ########         LEGEND         ########
c                   ########                        ########
c                   ########################################
c                   ########################################
c                   ########################################
c
c
c
      Subroutine Legend (nitems,strings,ltypes,lcolors,
     >  px1,py2,ibox,siz)
      Implicit None
c
      Integer nitems, ltypes(nitems), lcolors(nitems),
     >  ibox, i, lnb, Lnblnk
      Real	py1,py2,py3,py4, px1,px2,px3,px4, yy,y1,y2, dpx,dpy, 
     >  py,px, xx,x1,x2,x3,x4,x5, siz0,siz, pyl,yyl
      Character*(*) strings(nitems)
c
c  px1 = left edge, px2 = start of line, px3 = end of line, px4 = start of str
c
c     siz	= 0.02
      Call Xqchsz(siz0) !siz0 in mathematical space
      Call Xchmag(siz)	!siz in plotter space
c
      px	= px1
      dpy	= 1.2 * siz
      dpx	= 1.5 * siz
      px2	= px1 + dpx
      px3	= px2 + 3. * dpx
      px4	= px3 + 0.5 * dpx
      Call Plt2Wrld (px1,yy,x1,yy)
      Call Plt2Wrld (px2,yy,x2,yy)
      Call Plt2Wrld (px3,yy,x3,yy)
      Call Plt2Wrld (px4,yy,x4,yy)
      py	= py2 - 0.5 * dpy
c
      Do i=1,nitems
      py	= py - dpy
      pyl	= py + 0.3 * siz
      Call Lintyp (ltypes(i))
      Call Color (lcolors(i))
      Call Plt2Wrld (xx,py,xx,yy)
      Call Plt2Wrld (xx,pyl,xx,yyl)
      Call MyLine (x2,yyl,x3,yyl)
      lnb	= Lnblnk(strings(i))
      Call Xcharl (x4,yy,strings(i)(:lnb))
      End Do
c
      Call Xchsiz(siz0)	!siz0 in mathematical space
      Call Color (0)
      End
c
c
c
c
c                   ########################################
c                   ########################################
c                   ########################################
c                   ########                        ########
c                   ########        LABELER         ########
c                   ########                        ########
c                   ########################################
c                   ########################################
c                   ########################################
c
c
      Subroutine Labeler (xx,yy,zz,n, lftrgt,fmt,siz,ang)
      Implicit None
      Integer	n,k,Lnblnk,lftrgt,i_format
      Real	xx(n),yy(n),zz(n), siz,ang
      Character*(*) fmt
      Character*8 string
c
      Call Xchmag (siz)
      Call Xchori (ang)
c
c...Plot as integer or real, depending on format string
c
      i_format	= Index(fmt,'I') + Index(fmt,'i')
c
      Do k=1,n
      If (i_format.EQ.0) Then
        Write (string,fmt) zz(k)
      Else
        Write (string,fmt) Int(zz(k))
      End If
c
      If (lftrgt.EQ.0) Then
        Call Xcharr (xx(k),yy(k),string(:Lnblnk(string)+1))
      Else 
        Call Xcharl (xx(k),yy(k),string(:Lnblnk(string)))
      End If
      End Do
c
      Call Xchori (0.)
      End
