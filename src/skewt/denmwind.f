c
c     ##################################################################
c     ##################################################################
c     ######                                                      ######
c     ######                  Subroutine DENMW                    ######
c     ######                                                      ######
c     ######                Copyright (c) 1993                    ######
c     ######    Center for Analysis and Prediction of Storms      ######
c     ######    University of Oklahoma.  All rights reserved.     ######
c     ######                                                      ######
c     ##################################################################
c     ##################################################################
c
      Subroutine Denmw (u, v, rhobar, z, npts, u06, v06)
c
c#######################################################################
c
c     Purpose:
c
c     Calculate the density-weighted mean wind for the following layers:
c
c                              0-3  km
c                              0-6  km
c                              0-12 km
c                              3-6  km
c                              3-9  km
c                              3-12 km
c
c#######################################################################
c
c     Author: Steve Lazarus
c             4/15/93
c
c     Modification History:
c
c
c#######################################################################
c
c     Input :
c
c       u        x component of velocity from sounding (m/s)
c       v        y component of velocity from sounding (m/s)
c       z        Sounding height   (m)
c       rhobar   Sounding density (kg/m**3)
c
c       npts     Number of levels in sounding
c
c     Output:
c
c       u500     Density-weighted u-velocity over 0-500 m.
c       v500     Density-weighted v-velocity over 0-500 m.
c       u03      Density-weighted u-velocity over 0-3 km.
c       v03      Density-weighted v-velocity over 0-3 km.
c       u06      Density-weighted u-velocity over 0-6 km.
c       v06      Density-weighted v-velocity over 0-6 km.
c       u12      Density-weighted u-velocity over 0-12 km.
c       v12      Density-weighted v-velocity over 0-12 km.
c       u36      Density-weighted u-velocity over 3-6 km.
c       v36      Density-weighted v-velocity over 3-6 km.
c       u39      Density-weighted u-velocity over 3-9 km.
c       v39      Density-weighted v-velocity over 3-9 km.
c       u312     Density-weighted u-velocity over 3-12 km.
c       v312     Density-weighted v-velocity over 3-12 km.
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

      real u500,v500  ! 0-500 m mean winds.
      real u03,v03    ! 0-3 km mean winds.
      real u06,v06    ! 0-6 km mean winds.
      real u12,v12    ! 0-12 km mean winds.
      real u36,v36    ! 3-6 km mean winds.
      real u39,v39    ! 3-9 km mean winds.
      real u312,v312  ! 3-12 km mean winds.

      real u(npts)    ! x component of velocity from sounding (m/s)
      real v(npts)    ! y component of velocity from sounding (m/s)
      
      real z(npts)      ! Sounding height (m)
      real rhobar(npts) ! Sounding density (kg/m**3)

c
c#######################################################################
c
c     Misc. local variables: 
c
c#######################################################################
c
      integer k         ! K-index for do-loop.      

      real denu500,denv500,denu03,denv03,denu06,denv06
      real denu12,denv12,denu312,denv312,denu36,denv36
      real denu39,denv39
      real zsum,zsum500,zsum03,zsum06,zsum36,zsum39,zsum312
      real dz 
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
c     Determine the 0-3,0-6,0-12,3-6,3-9,& 3-12 km mean winds from 
c     input sounding.  
c     
c     NOTE: Assumes that the winds have been interpolated to an evenly-
c           spaced (dz=500 m) grid colocated with the density, i.e.
c
c                   rhobar,u,v --|-- 1250 m   k=4
c                                |
c                   rhobar,u,v --|--  750 m   k=3   Z
c                                |
c                   rhobar,u,v --|--  250 m   k=2
c
c            If the grid differs from that shown above, changes in the
c            following conditionals must be made!
c
c
c#######################################################################
c
      DO 300 k=2,npts
        IF (z(k).le.11750.) THEN

          dz     = z(k+1) - z(k)
          denu12 = denu12 + u(k)*rhobar(k)*dz
          denv12 = denv12 + v(k)*rhobar(k)*dz
          zsum   = zsum + rhobar(k)*dz

           IF (z(k).le.250.) THEN
         
             denu500 = denu12
             denv500 = denv12
             zsum500 = zsum

           END IF

           IF (z(k).le.2750.) THEN

             denu03 = denu12
             denv03 = denv12
             zsum03 = zsum

           END IF

           IF (z(k).ge.3250.) THEN

             denu312 = denu312 + u(k)*rhobar(k)*dz
             denv312 = denv312 + v(k)*rhobar(k)*dz
             zsum312 = zsum312 + rhobar(k)*dz

               IF (z(k).le.5750.) THEN

                 denu36 = denu36 + u(k)*rhobar(k)*dz
                 denv36 = denv36 + v(k)*rhobar(k)*dz
                 zsum36 = zsum312

               END IF

               IF (z(k).le.8750.) THEN

                 denu39 = denu312
                 denv39 = denv312
                 zsum39 = zsum312

               END IF

           END IF

           IF (z(k).le.5750.) THEN

             denu06 = denu12
             denv06 = denv12
             zsum06 = zsum

           END IF

        END IF 

300   CONTINUE

      u500 = denu500/zsum500
      v500 = denv500/zsum500
      u03  = denu03/zsum03
      v03  = denv03/zsum03
      u06  = denu06/zsum06
      v06  = denv06/zsum06
      u12  = denu12/zsum
      v12  = denv12/zsum
      u36  = denu36/zsum36
      v36  = denv36/zsum36
      u39  = denu39/zsum39
      v39  = denv39/zsum39
      u312 = denu312/zsum312
      v312 = denv312/zsum312

c
c#######################################################################
c
c     Write out the density-weighted winds
c
c#######################################################################
c       
      write(6,400)
400   format(45x,'  U      V')
      write(6,425) u500,v500
425   format(2x,'Density-weighted mean wind 0-500  m     ',2f7.1)
      write(6,450) u03,v03
450   format(2x,'Density-weighted mean wind 0-3   km     ',2f7.1)
      write(6,475) u06,v06
475   format(2x,'Density-weighted mean wind 0-6   km     ',2f7.1)
      write(6,500) u12,v12
500   format(2x,'Density-weighted mean wind 0-12  km     ',2f7.1)
      write(6,525) u36,v36
525   format(2x,'Density-weighted mean wind 3-6   km     ',2f7.1)
      write(6,550) u39,v39
550   format(2x,'Density-weighted mean wind 3-9   km     ',2f7.1)
      write(6,575) u312,v312
575   format(2x,'Density-weighted mean wind 3-12  km     ',2f7.1)

c
      end   
