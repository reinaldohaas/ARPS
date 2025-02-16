c
c
c  Richard Carpenter
c  April 1993
c
c
c
c                   ########################################
c                   ########################################
c                   ########################################
c                   ########                        ########
c                   ########          BARB          ########
c                   ########                        ########
c                   ########################################
c                   ########################################
c                   ########################################
c
c
      Subroutine Barb (x1,y1,dir,spd, junits,shaftlen)
      Implicit None
      Real	uu,vv,dpx,dpy,x1,y1,x2,y2,dir,spd,deg2rad,spd0,
     >  px1,px2,px3,py1,py2,py3,dir0, angbarb,angbarb0, x3,y3,
     >  featherlen, featherspace, shaftlen, dist
      Integer	nfull,nfifty,nbarb,nhalf,ibarb,jbarb,ififty,junits
      Data deg2rad /0.017453/!, shaftlen/0.035/
c
c
      featherlen= shaftlen * .5
      featherspace = shaftlen * .12
c
      uu        = spd * Sin(deg2rad*dir)	!really -U
      vv        = spd * Cos(deg2rad*dir)	!really -V
      If (junits.EQ.0) Then
        nhalf	= Nint (spd/2.5)	!for m/s
      Else If (junits.EQ.1) Then
        nhalf	= Nint (spd/5.0)	!for mph/knots
      End If
      nfifty	= nhalf / 10
      nhalf	= nhalf - nfifty * 10
      nfull	= nhalf / 2
      nhalf	= nhalf - nfull * 2
      nbarb	= nfifty+nfull+nhalf
      If (nfifty.EQ.0 .AND. nfull.EQ.0 .AND. nhalf.EQ.1) nbarb = 2 !fudge 5 kt
c     Print *, spd, uu,vv, nfifty,nfull,nhalf
c
      spd0	= .001 * spd/10.
      spd0	= .001
      dir0	= 90. - dir
      angbarb	= dir + 60.
      angbarb0	= 90. - angbarb

c
      Call Wrld2Plt (x1,y1,px1,py1)

c  calm

      If (spd.LE.0.0) RETURN
c
c  shaft
c
      dpx	= shaftlen * Cos(deg2rad*dir0)
      dpy	= shaftlen * Sin(deg2rad*dir0)
      px2	= px1 + dpx
      py2	= py1 + dpy
      Call Plt2Wrld (px2,py2,x2,y2)
      Call MyLine (x1,y1,x2,y2)
c
c  feathers
c
      Do ibarb=1,nbarb
      dpx	= featherlen * Cos(deg2rad*angbarb0)
      dpy	= featherlen * Sin(deg2rad*angbarb0)
c
      If (ibarb.EQ.nbarb .AND. nhalf.GE.1) Then
        dpx	= dpx * .5
        dpy	= dpy * .5
      End If
c
      px3	= px2 + dpx
      py3	= py2 + dpy
      Call Plt2Wrld (px3,py3,x3,y3)
c
c..fudge for 5 knot wind
c
      If (ibarb.EQ.1 .AND. nfifty.EQ.0 .AND. nfull.EQ.0) Then
      Else
	Call MyLine (x2,y2,x3,y3)
      End If
c
c  ...reset posn for next barb
c
      ififty	= Min (ibarb,nfifty)
      dist	= shaftlen - (ibarb + ififty) * featherspace
c     If (nfifty.GT.0 .AND. ibarb.EQ.nfifty+1)
c    >  dist	= shaftlen - (ibarb + nfifty+1) * featherspace
      dpx	= dist * Cos(deg2rad*dir0)
      dpy	= dist * Sin(deg2rad*dir0)
      px2	= px1 + dpx
      py2	= py1 + dpy
      Call Plt2Wrld (px2,py2,x2,y2)
c
c  ...fifty barb
c
      If (ibarb.LE.nfifty) Then
        Call MyLine (x2,y2,x3,y3)
	jbarb	= ibarb+Min(ibarb,nfifty)
        dpx	= (shaftlen-jbarb*featherspace) * Cos(deg2rad*dir0)
        dpy	= (shaftlen-jbarb*featherspace) * Sin(deg2rad*dir0)
        px2	= px1 + dpx
        py2	= py1 + dpy
        Call Plt2Wrld (px2,py2,x2,y2)
      End If
c
      End Do
c
      End
c
c
