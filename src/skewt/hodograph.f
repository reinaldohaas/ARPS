c
c                               Richard Carpenter
c                                Univ. of Oklahoma
c                                  August 1992
c
c
c
c
c                   ########################################
c                   ########################################
c                   ########################################
c                   ########                        ########
c                   ########       HODOGRAPH        ########
c                   ########                        ########
c                   ########################################
c                   ########################################
c                   ########################################
c
c
      Subroutine Hodograph (pres, zz, theta, qv, uu, vv, 
     &                      n, isnd,sndfile, helcontrs, hodo_denmwind,
     &                      jdefcolor, jgray, jhodoclr,
     &                      px1, px2, py1, py2)

      IMPLICIT NONE
      Include 'thermo.consts'
c
c
c  UNITS: All quantities are SI units unless labeled otherwise (e.g., 
c  presmb, tempc). Mixing ratios are kg/kg.
c
c  m/s        = mph / 2.237 = kts / 1.944
c  mph        = kts / 1.151
c
c----------------------------------------------------------------------=
c
      Integer        nmax  , i     , k     , n     , j     , isnd, 
     >               ncl   , mode  , mhodo1, nhodo , 
     >               jdefcolor, jgray, jhodoclr
      Logical        exist , helcontrs, hodo_denmwind

      Parameter (nmax=1000, nhodo=41)

      Real            
     >  px1   , px2   , py1   , py2   , 
     >  deg2rad, pi   , rhodo , uuold , vvold , duu ,
     >  cl(100),u06   , v06   , csqrt2,circleint,radius,
     >  px    , py    , stormu, stormv, pxc   , pyc   , umax  , umin  ,
     >  vmax  , vmin  , yy    , xx    ,circlemax,temp,  tmp

      Real                
     >  zzkm(nmax)    ,
     >  zz(n)    ,
     >  theta(n) , rho(nmax)   ,
     >  pres(n)  , qv(n)    ,
     >  uu(n)    , vv(n)    ,
     >  uu_grid(nhodo,nhodo), vv_grid(nhodo,nhodo),
     >  hel3km(nhodo,nhodo) , iwork(nhodo*nhodo)  

c     REAL uuold, vvold
      Character*(*) sndfile
      Character*79 string, stormfile
      Data        stormfile/'storm.track'/

      INCLUDE 'thermo.stfunc'
c
c
c----------------------------------------------------------------------=
c                        READ THE SOUNDING
c----------------------------------------------------------------------=
c

      pi        = 4. * Atan (1.)
      deg2rad        = pi/180.
c
c
      Do k=1,n
      zzkm(k)        = zz(k) * 1.e-3
      temp        = theta(k) * (pres(k) * p00inv) ** rcp
      If (rho(k).EQ.0.) rho(k) = pres(k) * rdi /
     >                               Ftvirtnowl(temp,qv(k))
c     print '(i3,f8.0,4f8.2)', k, zz(k), uu(k), vv(k), rho(k)
      End Do
c
c----------------------------------------------------------------------=
c                                PLOTTING
c----------------------------------------------------------------------=
c
      Print *, '...Drawing Hodograph'
      pxc        = px1 + .5 * (px2-px1)
      pyc        = py1 + .5 * (py2-py1)
c
c...Find radius of plot based on max winds.
c
      DO i=1,n
        IF (ABS(uu(i)) .LT. 200.0) THEN
          umax = MAX(umax,uu(i))
          umin = MIN(umin,uu(i))
        END IF
      END DO
!     Call Mxmn1 (uu, umax,umin, 1,n, 1,n)
!     Call Mxmn1 (vv, vmax,vmin, 1,n, 1,n)
      umax        = Max (Abs(umax),Abs(umin),Abs(vmax),Abs(vmin)) - 4.
      umax        = Int(umax+10.0)/10 * 10.0 + 2.0
      umax        = MIN( umax, 62.0 )
      umax        = MAX( umax, 42.0 )
      vmax        = umax
      umin        = -umax
      vmin        = -vmax
c
      Call Xmap (umin, umax, vmin, vmax)
      CALL XPSPAC (px1, px2, py1, py2) 
c
      Call Color (jgray)
      Call Lintyp (01)
c
c...Draw circles & axes
c
      IF (isnd.EQ.1) THEN

      Call Xchmag (0.015)
      circlemax = Int (umax/10.) * 10.0
      circleint = 10.0
c
      Do i=Int(circleint),Int(circlemax),Int(circleint)
        radius  = i
        Call Xcircle (radius)
        Write (string,'(I2)') i
        tmp        = radius / SQRT(2.0)
        Call Xcharc (-tmp,-tmp,TRIM(string))
      End Do
c
      Call MyLine (-umax,-umax, umax,-umax)
      Call MyLine (umax,-umax, umax,umax)
      Call MyLine (-circlemax,0., circlemax,0.)
      Call MyLine (0.,-circlemax, 0.,circlemax)
      csqrt2    = circlemax/Sqrt(2.)
      Call MyLine (-csqrt2,-csqrt2, csqrt2,csqrt2)
      Call MyLine (-csqrt2,csqrt2, csqrt2,-csqrt2)
c
      Call Plt2Wrld (px1+0.01,py1+0.01,xx,yy)
      Write (string,'(2a)') 'm/s'
      Call Xcharl (xx,yy,TRIM(string))

      END IF
c
c...Plot the sounding
c
      Call Color (jhodoclr)
      uuold = uu(1)
      vvold = vv(1)
      DO i=1,n
        IF (ABS(uu(i)) .LT. 200.0 .AND. ABS(vv(i)) .LT. 200.0) THEN
          CALL MyLine(uuold,vvold, uu(i),vv(i))
          uuold = uu(i)
          vvold = vv(i)
          CALL Xscatter (uu(i),vv(i),1, 0,0.005)
        END IF
      END DO
!     Call Xcurve (uu, vv, n, 0)
!     Call Color (jhodoclr)
!     Call Xscatter (uu,vv,n, 0,0.005)

      IF (isnd.GT.1) RETURN
c
c...Label points every km
c
c     Call Labeler (uu,vv,zzkm,1, 0,'(F4.1)',0.012,0.)
      zzkm(n+1)        = 9999.
c
      duu        = 5.0
      uuold        = -99.0
      vvold        = -99.0
      DO j=2,n        ! RLC 1994/03/31
        IF ( (ABS(uuold-uu(j)).GT.duu) .AND.
     &       (ABS(vvold-vv(j)).GT.duu) ) THEN
            uuold        = uu(j)
            vvold        = vv(j)
            Call Labeler (uu(j),vv(j),zzkm(j),1, 0,'(F4.1)',0.014,0.)
        END IF
      END DO
c
c
c...Read in Mean Storm Motion
c
      Call Xchmag (0.02)
      Print *, 'Storm motion file: ', stormfile
      Inquire (File=stormfile, Exist=exist)
      If (exist) Then
        Open (Unit=1, File=stormfile, Status='Old')
        Read (1,*) stormu, stormv
        Call Color (11)
        Call Xscatter (stormu,stormv,1, 3,0.020)
      Else
        Print *, 'Storm motion file not found, ignoring'
      End If
c
c
c...Legend
c
      px        = px1 + 0.05
      py        = py2 - 0.01

!  legend for storm motion

      IF (exist) THEN
      Call Plt2Wrld (px1+0.01,py,xx,yy)
      Call Color (jdefcolor)
      Call Xscatter (xx,yy,1, 3,0.020)
      Call Plt2Wrld (px,py,xx,yy)
      Call Xcharl (xx,yy,'Mean storm motion (model)')
      END IF

!  legend for density wtgd mean wind

      IF (hodo_denmwind) THEN
      py        = py - 0.04
      Call Plt2Wrld (px1+0.01,py,xx,yy)
      Call Color (jdefcolor)
      Call Xscatter (xx,yy,1, 1,0.020)
      Call Plt2Wrld (px,py,xx,yy)
      Call Xcharl (xx,yy,'0-6 km mean wind (sounding)')
      END IF
c
c
c...Helicity contours
c
      IF (helcontrs) THEN

      rhodo        = Int (umax/10.)
      mhodo1        = rhodo + 1.
      rhodo        = nhodo - 1.
      mhodo1        = nhodo
      Print *, rhodo, mhodo1, nhodo
      Call Helcont (uu, vv, zz, n, hel3km, rhodo, mhodo1)

      cl(1)        = 0.
      cl(2)        = 200.
      mode        = 1
c
      Do i=1,nhodo
      Do j=1,nhodo
        uu_grid(i,j) = -rhodo + (i-1)*2.
        vv_grid(i,j) = -rhodo + (j-1)*2.
      End Do
      End Do
c
      Call Color (25)
c
c...Donnot contour the whole array
c
      Call Xwindw (umin,umax,0.85*vmin,0.85*vmax)  !turn window clipping on
      Call Xchmag (0.012)
      Call Xclfmt ('(I5)')
      Call Xnctrs (6,16)
      Call Xhilit (0)
c     write(0,*) 'before xconta1'
      Call Xconta (hel3km,uu_grid(1,1),vv_grid(1,1), iwork,
     >  nhodo,nhodo,nhodo, cl,ncl,mode)
c     Call Xconta (hel3km,uu_grid,vv_grid, iwork, 61,61,61, cl,ncl,mode)
      Call Xwdwof                !turn window clipping off
c
      Call Xchmag (0.015)
      Call Plt2Wrld (px1+0.02,py1-0.01,xx,yy)
      Write (string,'(a,I3)')
     >  'Contours of 0-3 km storm-relative helicity, intvl=',
     >  Int(cl(2)-cl(1))
      Call Xcharl (xx,yy,string)

      END IF
c
c...Density-weighted mean wind
c
      IF (hodo_denmwind) THEN
      Call Denmw (uu, vv, rho, zz, n, u06, v06)
      Call Color (04)
      Call Xscatter (u06,v06,1, 1,0.02)
      END IF
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
c                   ########         XCIRCLE        ########
c                   ########                        ########
c                   ########################################
c                   ########################################
c                   ########################################
c
c
      Subroutine Xcircle (radius)
      Implicit None
      Integer        i, inc
      Real        xx,yy, radius, deg2rad, pi, angle
      Parameter (pi=3.14159265, deg2rad=pi/180.)
c
c...assumes aspect ratio of one.
c
      Call Xpenup (radius,0.)
c
      inc = 10
      Do i=inc,360,inc
      angle        = deg2rad * i
      xx        = radius * Cos (angle)
      yy        = radius * Sin (angle)
      Call Xpendn (xx,yy)
      End Do
c
      End

