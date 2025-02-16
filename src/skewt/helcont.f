c
c     ##################################################################
c     ##################################################################
c     ######                                                      ######
c     ######                SUBROUTINE HELCONT                    ######
c     ######                                                      ######
c     ######                Copyright (c) 1993                    ######
c     ######    Center for Analysis and Prediction of Storms      ######
c     ######    University of Oklahoma.  All rights reserved.     ######
c     ######                                                      ######
c     ##################################################################
c     ##################################################################
c
c
       SUBROUTINE HELCONT(u,v,z,npts,hel3km,rhodo,mhodo1)
c
c#######################################################################
c
c     Purpose:
c
c     Calculate the storm-relative 0-3 km integrated helicity, H:
c
c               H = INTEGRAL [ (v-Cy)du/dz - (u-Cx)dv/dz ]dz,
c
c     where     Cx = x-component of the storm motion vector
c               Cy = y-component of the storm motion vector
c            du/dz = x-component of the shear vector
c            dv/dz = y-component of the shear vector
c
c     See Davies-Jones (1990) Severe Storms Conference Preprint  
c
c
c#######################################################################
c
c     Author: Steve Lazarus
c             4/24/93
c
c     Modification History:
c		RLC 04/27/93
c
c
c#######################################################################
c
c     Input :
c
c       u        x component of velocity from sounding (m/s)
c       v        y component of velocity from sounding (m/s)
c       z        Sounding height   (m)
c
c       npts     Number of levels in sounding
c	rhodo	Radius of hodograph (m/s)
c	mhodo1	Dimension of hel3km
c
c  Note: rhodo should be divisible by 10. mhodo1 should equal rhodo + 1.
c
c     Output:
c
c       hel3km   0-3 km integrated helicity as a function of storm
c                motion f(Cx,Cy)
c
c#######################################################################
c

c
c#######################################################################
c
c     Variable declarations.
c
c#######################################################################
c
      implicit none 

      integer npts    ! Number of points in sounding.

      real u(npts)    ! x component of velocity from sounding (m/s)
      real v(npts)    ! y component of velocity from sounding (m/s)
      
      real z(npts)      ! Sounding height (m)

c
c#######################################################################
c
c     Misc. local variables: 
c
c#######################################################################
c
      integer i,j,k     ! Indices for do-loops.      

      real c,c3
      real dudz,dvdz
      real deltz,dz     ! Interpolation intervals
      real u3km,v3km    ! Interpolated u,v 3 km velocity components
      real cx,cy        ! Storm motion in x and y directions

      real rhodo        ! Hodograph radius
      integer mhodo1

      real hel3km(mhodo1,mhodo1)
c
c#######################################################################
c
c     Initialize Variables
c
c#######################################################################
c

c
c@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
c
c     Beginning of executable code...
c
c@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
c

c
c#######################################################################
c
c     Set the storm motion components Cx and Cy as a function of 
c     hodograph radius, where:
c
c           rhodo = (hodograph radius) 
c   
c#######################################################################
c
	DO 400 i = 1,mhodo1

          cx = -rhodo + float(i-1)*2.0

        DO 300 j = 1,mhodo1

          cy = -rhodo + float(j-1)*2.0
        
          hel3km(i,j)=0.0
c
c#######################################################################
c
c     Calculate the 0-3 km Helicity as a function of Cx and Cy
c     If z(k) does not equal 3000. then linearly interpolate the
c     winds (u,v) to this level.
c
c#######################################################################
c
          DO 200 k=2,npts-1
            c = (u(k+1) - cx) * (v(k) - cy)
            c = c - (u(k) - cx) * (v(k+1) - cy)
            c3= c
            IF (z(k).lt.3000.0.and.z(k+1).gt.3000.0) THEN

              dz   = z(k+1) - z(k)
              dudz = (u(k+1)-u(k))/dz
              dvdz = (v(k+1)-v(k))/dz
              deltz= 3000. - z(k)
              u3km = u(k) + dudz*deltz
              v3km = v(k) + dvdz*deltz
              c3   = (u3km-cx) * (v(k)-cy) - (u(k)-cx) * (v3km - cy)
              
            END IF
            
            IF (z(k).lt.3000.0) hel3km(i,j) = hel3km(i,j) + c3

            c  = 0.0
            c3 = 0.0

200      CONTINUE

300      CONTINUE

400    CONTINUE
       RETURN
       END
