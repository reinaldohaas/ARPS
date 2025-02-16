c  Last file modification: 4/16/98

c
c
c
c
c                   ########################################
c                   ########################################
c                   ########################################
c                   ########                        ########
c                   ########       PLOTPRECIP       ########
c                   ########                        ########
c                   ########################################
c                   ########################################
c                   ########################################
c
c
      SUBROUTINE PlotPrecip (presmb, n, pspacing, jscheme,
     &  label,prectype_1d)
      IMPLICIT NONE   
      INTEGER   k,n, jscheme
      INTEGER prectype_1d(n) ! EMK
      REAL      presmb(n), pspacing
      REAL      xx,yy,yykm1,xmin,xmax,ymin,ymax, dlnp, px1,px2,py1,py2
      REAL      px, siz0
      DATA      px/0.0/, siz0/0.02/
      CHARACTER*(*) label
      CHARACTER*1 symbol ! EMK
      SAVE      px
c
c  spacing between printed levels
c
      CALL XQMAP (xmin,xmax,ymin,ymax)
      dlnp      = Abs(ymax-ymin)/pspacing
      yykm1     = 10000.   
c
c  location of precip. symbols
c
      IF (px .EQ. 0.0) THEN
        CALL XQPSPC (px1,px2,py1,py2)
        px = px2 + 0.075 ! EMK
*        px = px2 + 0.11 ! Original setting
        CALL Plt2Wrld (px,0.05,xx,yy)
      ELSE
        px = px + 0.15 ! EMK
*        px = px + 0.05 ! Original setting
        CALL Plt2Wrld (px,0.0,xx,yy)
      END IF  
c
      call XCHMAG(1.5*siz0)
*      call XCHMAG(1.2*siz0) 

      DO k=1,n-1   
*      DO k=1,n   
      yy        = LOG (presmb(k))
      IF ( yykm1-yy.GT.dlnp ) THEN
        yykm1   = yy
!       IF (jscheme.NE.2) THEN
!        IF (spd(k).GT.117.5/2.) THEN
!         CALL Color (04)
!        Else IF (spd(k).GT.77.5/2.) THEN
!         CALL Color (24)
!        Else IF (spd(k).GT.37.5/2.) THEN
!         CALL Color (23)
!        Else
!         CALL Color (01)
!       END IF
!       END IF
        if (prectype_1d(k).eq.1) then
         symbol = '*'
        else if (prectype_1d(k).eq.2) then
         symbol = 'P'
        else if (prectype_1d(k).eq.3) then
         symbol = '~'
        else if (prectype_1d(k).eq.4) then
         symbol = '.'
c         symbol = ','
        else
         symbol = ' '
        end if
        if (symbol.eq.',' .or. symbol.eq.',' ) then
         call xcharc(xx-0.5,yy,symbol)
        else
         Call xcharc(xx-0.5,yy+0.0425,symbol)         
        end if
      END IF
      END DO

      CALL XCHMAG (0.6*siz0)
      
      CALL Plt2Wrld (px,py1-0.02,xx,yy)
      CALL XCHARC (xx,yy,TRIM(label))
     
c
      IF (jscheme.NE.2) CALL Color (01)

      END


c
c
c     ##################################################################
c     ##################################################################
c     ######                                                      ######
c     ######           SUBROUTINE PRECTYPE_LAPS_SKEWT             ######
c     ######                                                      ######
c     ######                Copyright (c) 1998                    ######
c     ######    Center for Analysis and Prediction of Storms      ######
c     ######    University of Oklahoma.  All rights reserved.     ######
c     ######                                                      ######
c     ##################################################################
c     ##################################################################
c
 
      SUBROUTINE PRECTYPE_LAPS_SKEWT (nmax,nz,temp_1d,p_pa_1d,tw_1d,
     :               cldpcp_type_1d)
 
c
c#######################################################################
c
c     PURPOSE:
c     This routine returns a 1-D precipitation type array using the
c     LAPS algorithm.
c
c#######################################################################
c       
c     AUTHOR:  Jian Zhang
c     05/1996  Based on the LAPS cloud analysis code developed by
c              Steve Albers. 
c
c     This program modifies the most significant 4 bits of the integer
c     array by inserting multiples of 16.
c
c     MODIFICATION HISTORY:   
c
c     05/16/96 (J. Zhang)
c              Modified for ADAS format. Added full documentation.
c
c     01/16/98 (E. Kemp, Project COMET-Tinker)
c              Modified for use in ARPSSKEWT (creating subroutine
c              PRECTYPE_LAPS_SKEWT from subroutine PCP_TYPE_3D). Added
c              some new documentation.c
c
c     01/23/98 (E. Kemp, Project COMET-Tinker)
c              Corrected bugs involving the retesting of sleet for
c              freezing in a sub-zero wet bulb temperature layer.
c
c     01/31/98 (E. Kemp, Project COMET-Tinker)
c              Restrict analyzing precip. type to below 500mb.
c
c     04/16/98 (E. Kemp, Project COMET-Tinker)
c              Added check for supercooled liquid (freezing rain) at
c              the top of a reflectivity column, per changes to the
c              LAPS code at FSL.
c
c     05/13/98 (E. Kemp, Project COMET-Tinker)
c              Corrected setting of melting flag for freezing rain.
c              Also corrected generation of rain at top of column for
c              Tw between 0 and 1.3 C.
c
c#######################################################################
c
c#######################################################################
c
c     Variable Declarations.
c
c#######################################################################
c
      implicit none
c
c#######################################################################
c     
c     INPUT:
      integer nmax
      integer nz                  ! Model grid size
      real temp_1d(nmax)            ! temperature (K)
      real tw_1d(nmax)              ! Wet bulb temperature (K) (EMK)
      real p_pa_1d(nmax)            ! pressure (Pascal)
c
c     OUTPUT:
      integer istatus
      integer cldpcp_type_1d(nmax)! precip type
      integer*4 itype                   ! precip type index
c     
c     LOCAL functions:

c       
c#######################################################################
c
c     
c     Misc local variables   
c
c#######################################################################
c     
      integer k,k_upper
      real t_c,t_wb_c,temp_lower_c,temp_upper_c,tbar_c
     :     ,thickns,frac_below_zero
      integer iprecip_type,iprecip_type_last,iflag_melt
     :        ,iflag_refreez
      real zero_c,rlayer_refreez_max,rlayer_refreez
      integer n_zr,n_sl,n_last
      real tmelt_c
c              
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
C
C     Beginning of executable code...
C
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
c
c#######################################################################
c
      istatus=0
c
c
c#######################################################################
c
c     Stuff precip type into cloud type array
c     No Precip = 0
c     Rain = 4
c     Snow = 1
c     Freezing Rain = 3
c     Sleet = 2
c     Hail = 5
c
c#######################################################################
c
      zero_c = 273.15
      rlayer_refreez_max = 0.0
      
      n_zr = 0
      n_sl = 0
      n_last = 0
        
 
        iflag_melt = 0
        iflag_refreez = 0
        rlayer_refreez = 0.0
      
        iprecip_type_last = 0
     
      DO k = nz-1,1,-1
     
c
c#######################################################################
c              
c     Set refreezing flag
c
c#######################################################################
c
c        if (p_pa_1d(k).lt.50000.) then ! EMK 1/31/98
c          iprecip_type = 0
c          goto 99
c        end if
          
        t_c  = temp_1d(k)
        t_wb_c = tw_1d(k) - zero_c
      
        tmelt_c = 1.3 ! Radar reflectivity not available
 
        if(t_wb_c .lt. 0.) then
          if(iflag_melt .eq. 1) then
c     
c#######################################################################
c     
c     Integrate below freezing temperature times column
c     thickness - ONLY for portion of layer below freezing
c     
c#######################################################################
c
            temp_lower_c = t_wb_c
            k_upper = min(k+1,nz-1)
c
c#######################################################################
c     
c     For simplicity and efficiency, the assumption is
c     here made that the wet bulb depression is constant
c     throughout the level.
c       
c#######################################################################
c
c
            temp_upper_c = t_wb_c + ( temp_1d(k_upper)
     :                                  - temp_1d(k))
            if(temp_upper_c .le. 0.) then
              frac_below_zero = 1.0
              tbar_c = 0.5 * (temp_lower_c + temp_upper_c)
      
            else ! Layer straddles the freezing level
              frac_below_zero = temp_lower_c
     :                                / (temp_lower_c - temp_upper_c)
              tbar_c = 0.5 * temp_lower_c
 
            endif
               
            thickns = p_pa_1d(k_upper) - p_pa_1d(k)
            rlayer_refreez = rlayer_refreez
     :             + abs(tbar_c * thickns * frac_below_zero)
 
            if(rlayer_refreez .ge. 25000.) then
              iflag_refreez = 1
            endif
          
            rlayer_refreez_max =
     :                        max(rlayer_refreez_max,rlayer_refreez)
 
          endif ! iflag_melt = 1
            
        else ! Temp > 0C
          iflag_refreez = 0
          rlayer_refreez = 0.0
 
        endif ! T < 0.0c, Temp is below freezing
c     
c#######################################################################
c
c     Set melting flag
c
c#######################################################################
c
        if(t_wb_c .ge. tmelt_c) then 
          iflag_melt = 1
        endif
 
        if(t_wb_c .ge. tmelt_c) then  ! Melted to Rain
          iprecip_type = 4
        else ! Check if below zero_c (Refrozen Precip or Snow)
          if(t_wb_c .lt. 0.0) then
            if (iprecip_type_last .eq. 0)then  ! Generating lyr
              if (t_wb_c .ge. -6.)then       ! Supercooled pcp
*              if (.false.)then               ! Supercooled pcp ! EMK
                iflag_melt = 1 ! EMK
*                iflag_melt = 0 ! LAPS
                iprecip_type = 3          ! Freezing Rain
              else
                iprecip_type = 1          ! Snow
              endif
            else if (iflag_melt .eq. 1) then
              if ((iprecip_type_last.eq.4).or.
     +            (iprecip_type_last.eq.3)) then ! Test RN or FZ RN for freezing
                if(iflag_refreez .eq. 0) then ! Freezing Rain
                  iprecip_type = 3
                else  ! (iflag_refreez = 1)  ! Sleet
                  n_sl = n_sl + 1
                  iprecip_type = 2
                endif ! iflag_refreez .eq. 0
              else ! Precip. above is sleet, remains unchanged
                iprecip_type = iprecip_type_last  
                n_last = n_last + 1
                if(n_last .lt. 5) then
                  write(6,*)'Unchanged Precip',k,t_wb_c
                endif
              endif ! liquid precip. above level being tested?
               
            else    ! iflag_melt =0        ! Snow
              iprecip_type = 1
            endif   ! iflag_melt = 1
              
          else ! t_wb_c >= 0c, and t_wb_c < tmelt_c
     
            if (iprecip_type_last.eq.0) then        !   1/20/98
              iprecip_type = 4    ! rain:at echo top and 0<Tw<1.3C
              iflag_melt = 1
            else
              iprecip_type = iprecip_type_last
              n_last = n_last + 1
              if(n_last .lt. 5) then
                write(6,*)'Unchanged Precip',k,t_wb_c
              endif
            end if

          endif ! t_wb_c < 0c
        endif  ! t_wb_c >= tmelt_c
c
c#######################################################################
c
c     Insert most sig 4 bits into array
c
c#######################################################################
c

 99     CONTINUE
          
        itype = 0
        itype = itype - (itype/16)*16     ! Initialize the 4 bits
        itype = itype + iprecip_type      ! Add in the new value
              
        cldpcp_type_1d(k) = itype

        iprecip_type_last = iprecip_type
     
      ENDDO ! k
                  
      write(6,*)' rlayer_refreez_max = ',rlayer_refreez_max
      write(6,*)' n_frz_rain/n_sleet = ',n_zr,n_sl
      istatus=1
 
      RETURN
      END
          

      
c 
c  
 
c                   ########################################
c                   ########                        ########
c                   ########    MEPRECTYPESKEWT     ########
c                   ########                        ########
c                   ########################################
c                   ########################################
c          
c 
 
      function meprectypeskewt(nmax,nz,press,zp,tz,td,wbt)
 
c Function for calculating precipitation type from a sounding (assuming
c precipitation occurs)
c
c Based on NCEP MesoEta Precipitation Type Algorithm outlined in:
c  
c Baldwin, M. E., and S. P. Contorno, 1993:  Development of a Weather-Type
c Prediction System for NMC's Mesoscale Eta Model.  Preprints, 13th
c Conf. on Weather Analysis and Forecasting, Vienna, VA, Amer. Meteor.
c Soc., 86-87.
c  
c Updated algorithm described through personal communication with
c Mike Baldwin of the General Sciences Corporation.
c          
c Eric Kemp
c Project COMET-Tinker
c 12/4/97
c --
 
 
c  Description of Variables and Arrays:
c 
c  nmax == parameter for maximum possible number of levels in
c          sounding (defined in VERSKEWT).
c  nz == actual number of sounding levels (variable n in
c        VERSKEWT)
c  tz(nmax) == temperature (C) at each sounding level
c  td(nmax) == dew point (C) at each sounding level
c  press(nmax) == pressure (Pa) at each sounding level
c  Low150pres == pressure level (Pa) at the top of the lowest 150 mb in the sounding 
c                (the surface being the bottom of the layer)
c  delLow150lnp == change in natural logarithm pressure between Low150pres and the pressure of
c                  the sounding level immediately below it.
c  delLow150T == difference in temperature (C) between Low150pres and the
c                sounding level immediately below it.
c  delLow150z == difference in height (m) between Low150pres and the sounding level immediately
c                below it.
c  delLow150Td == difference in dew point (C) between Low150pres and the sounding level
c                 immediately below it.
c  Low150T == temperature (C) at Low150pres
c  Low150z == height (m) at Low150pres
c  Low150Td == dew point (C) at Low150pres
c  Low150WBT == wet bulb temperature at Low150pres
c  delT == difference in temperature (C) between the sounding levels immediately below and
c          immediately above Low150pres
c  delTd == difference in temperature (C) between the sounding levels immediately below and
c           immediately abovt Low150pres 
c  tw() == function for calculating wet bulb temperature, found in
c          u.skewt.f
c  wbt(nmax) == wet bulb temperature (K) at each sounding level (a
c               miscellaneous array is read in from VERSKEWT for
c               this variable array)
c  coldsatT == coldest temperature (K) of saturated layer for sounding
c              below 700 mb.
c  coldsatk == k-index of the bottom of the coldest saturated layer
c              for sounding below 700 mb.
c  LayerT == average temperature (K) of a layer (between two adjacent
c            levels in the sounding)
c  LayerTd == average dew point temperature (K) of a layer (between two
c             adjacent levels in the sounding)
c  Layerdewdep == average dew point depression (K) of a layer (between
c                 two adjacent levels in the sounding)
c  AWBT == area wet bulb temperature (K m)
c  zp(nmax) == height of sounding levels (m)
c  delWBT == change in wet bulb temperature (K) between two adjacent
c            vertical grid points
c  dellnp == change in natural logarithms of pressure between two adjacent
c            vertical grid points
c  delsnWBT == change in wet bulb temperature (K) between a level with wet
c              bulb temperature of 269 K and the grid point immediately
c              below it
c  delsnlnp == change in natural logarithms of pressure between a level
c              with wet bulb temperature of 269 K and the grid point
c              immediately below it
c  delz == change in height (m) between two adjacent sounding levels
c  snheight == height (m) of a level with wet bulb temperature of 269 K
c  snpress == pressure (Pa) of a level with wet bulb temperature of 269 K
c  AWBT1 == Area Wet Bulb Temperatures (K m) above 269 K from
c                  the surface to the coldest saturated layer for the
c                  sounding
c  frWBT == freezing wet bulb temperature (273.15 K)
c  kL150mb == k-index of the level 150mb less than the pressure at the
c            surface.
c  delfrWBT == change in wet bulb temperature (K) between a level at frWBT
c              and the sounding level immediately below it
c  delfrlnp == change in natural logarithms of pressure between a level   
c              at frWBT and the sounding level immediately below it
c  frheight == height (m) of a level at frWBT
c  frpress == pressure (Pa) of a level at frWBT
c  AWBT2 == net area wet bulb temperature (K) with respect to freezing
c                  in lowest 150mb 
c  AWBT3 == area (K m) of surface based wet bulb temperature above  
c                  freezing for the sounding.
c  AWBT4 == area (K m) of surface based wet bulb temperature below
c                  freezing for the sounding.
c  LowLayerT == average temperature (K) of lowest layer (between
c                      k=1 and k=2 sounding levels)
 
      implicit none

c  Function declaration
 
      integer meprectypeskewt
 
c  Input arrays and variables
 
      integer nz,nmax
      
      real press(nmax),zp(nmax),tz(nmax),td(nmax),wbt(nmax)
       
c  Internal arrays and variables used for calculations
 
      real coldsatT,LayerT,LayerTd,
     +     Layerdewdep,AWBT,delWBT,dellnp,delsnWBT,
     +     delsnlnp,delz,snheight,snpress,AWBT1,delfrWBT,
     +     delfrlnp, frheight,frpress,AWBT2,AWBT3,
     +     AWBT4,frWBT,lowlayerT,
     +     Low150pres,delLow150lnp,delLow150T,delLow150z,
     +     delLow150Td,Low150T,Low150z,Low150Td, Low150WBT,
     +     delT,delTd
      
      parameter (frWBT = 273.15)
 
      integer coldsatk, k
      
      real lowunsatdewdep ! new (EMK)

*  External function called in this function.
       
      real tw
            
*  Initialize ColdSatT so that, if no saturated layer exists, the program
*  will still check for all types of precip.

      ColdSatT = 9999.
*      ColdSatT = -9999.
      
*  Calculate wet bulb temperature at each sounding level
      
c      print*,'Calculating wet bulb temperatures at each sounding
c     +level...'

*      k = 1
*      do while (k.le.nz)
*       wbt(k) = tw((tz(k)),(td(k)),(press(k)/100.))
*       wbt(k) = wbt(k) + 273.15
*       print*,'k = ',k
*       print*,'Wet Bulb Temperature (K) = ',wbt(k)
*       k = k + 1
*      end do
      
*  Find the temperature of the coldest saturated layer (coldsatT), and 
*  the k-index of the top of the coldest saturated layer (coldsatk).

c      print*,'Finding coldest saturated layer...'

      Lowunsatdewdep = 9999. ! new variable (EMK)
      k = 1
      coldsatk = k
      do while ((k+1).le.nz)
       LayerT = (TZ(k)+TZ(k+1))/2.
       LayerT = LayerT + 273.15
       LayerTd = (Td(k)+Td(k+1))/2.
       LayerTd = LayerTd + 273.15
       Layerdewdep = LayerT - LayerTd
       If (Layerdewdep.le.(2.)) then
         lowunsatdewdep = 2.
         if (LayerT.le.coldsatT) then
           print*,'Saturated Layer w/ LayerT = ',LayerT
           print*,'ColdSatT = ',coldsatT
           coldsatT = LayerT
           print*,'ColdsatT changed to ',coldsatT
           coldsatk = k 
         end if         
         print*
       else if ((layerdewdep.le.lowunsatdewdep)) then
         print*,'Unsaturated layer w/ Layerdewdep = ',layerdewdep
         print*,'LayerT = ',layerT
         print*,'ColdSatT = ',coldsatT
         print*,'Lowunsatdewdep = ',lowunsatdewdep
         coldsatT = LayerT
         print*,'ColdsatT changed to ',coldsatT
         coldsatk = k
         lowunsatdewdep = layerdewdep
         print*,'Lowunsatdewdep changed to ',lowunsatdewdep
         print*
       end if
       k = k + 1
      end do
        
       
*  Calculate Area Wet Bulb Temperature Above 269 K (AWBT1) from the
*  surface to the coldest saturated layer.

c     print*,'Calculating Area Wet Bulb Temperature Above 269K...'      

      k = 1
      AWBT = 0.
      Do while ((k).le.(coldsatk))
       If ((WBT(k).le.(269.)).and.(WBT(k+1).le.(269.))) then
        k = k + 1
       else if ((WBT(k).ge.(269.)).and.(WBT(k+1).ge.(269.))) then
        AWBT = AWBT +
     +          ((WBT(k+1)+WBT(k)-(2.*269.))*
     +          (0.5)*(zp(k+1)-zp(k)))
        k = k + 1
       else
        delWBT=WBT(k+1)-WBT(k)
        dellnp=Alog(press(k+1)) - Alog(press(k))
        delsnWBT = 269. - WBT(k)
        delsnlnp = (delsnWBT)*(dellnp)/(delWBT)
        delz = zp(k+1) - zp(k)
        snheight=zp(k) + (delz*delsnlnp/dellnp)
        snpress=alog(Press(k)) + delsnlnp
        snpress = exp(snpress)
        if ((WBT(k).gt.(269.)).and.(WBT(k+1).lt.(269.))) then
         AWBT=AWBT +
     +         (((269.)+WBT(k)-(2.*269.))*
     +         (0.5)*(snheight-zp(k)))
        else
         AWBT=AWBT +
     +         ((WBT(k+1)+269.-(2.*269.))*
     +         (0.5)*(zp(k+1)-snheight))
        end if
        k = k + 1
       end if
      end do
      AWBT1 = AWBT
        
*  Calculate area wet bulb temperature with respect to 0 degrees C in lowest
*  150 mb (AWBT2)

c      print*,'Determining lowest 150mb of sounding...'     

      k = 1
      AWBT = 0.
      do while (((press(1)-press(k)).le.(15000.)).and.((k+1).le.nz))
       k = k + 1
      end do
      Low150pres = press(1) - 15000.
      delT = TZ(k+1) - TZ(k)
      dellnp = Alog(press(k+1)) - Alog(press(k))
      delLow150lnp = Alog(Low150pres) - Alog(press(k))
      delLow150T = (delLow150lnP)*(delT)/(dellnp)
      delz = zp(k+1) - zp(k)
      delLow150z = (delLow150lnP)*(delz)/(dellnp)
      delTd = Td(k+1) - Td(k)
      delLow150Td = (delLow150lnP)*(delTd)/(dellnp)
      Low150T = TZ(k) + delLow150T
      Low150z = zp(k) + delLow150z
      Low150Td = Td(K) + delLow150Td
      Low150WBT = tw(Low150T,Low150Td,(Low150pres/100.))
     
c      print*,'Calculating Net Area Wet Bulb Temperature with respect to
c     +freezing...'

      k = 1
      do while (press(k+1).ge.Low150pres)
       if ((WBT(k).le.(273.15)).and.(WBT(k+1).le.(273.15))) then
        AWBT = AWBT +
     +         ((WBT(k+1)+WBT(k)-(2.*frWBT))*
     +         (0.5)*(zp(k+1)-zp(k)))
        k = k + 1
       else if ((WBT(k).ge.(273.15)).and.
     +         (WBT(k+1).ge.(273.15))) then
        AWBT = AWBT +
     +         (((WBT(k+1)+WBT(k)-(2.*frWBT))*(0.5)*
     +         (zp(k+1)-zp(k))))
        k = k + 1
       else
        delWBT = WBT(k+1) - WBT(k)
        dellnp = Alog(press(k+1)) - Alog(press(k))
        delfrWBT = 273.15 - WBT(k)
        delfrlnp = (delfrWBT)*(dellnp)/delWBT
        delz = zp(k+1) - zp(k)
        frheight = zp(k) + (delz*delfrlnp/dellnp)
        frpress = alog(press(k)) + delfrlnp  
        frpress = exp(frpress)
        AWBT = AWBT +
     +         ((frWBT+WBT(k)-(2.*frWBT))*
     +         (0.5)*(frheight-zp(k)))
        AWBT = AWBT +
     +         ((WBT(k+1)+frWBT-(2.*frWBT))*
     +         (0.5)*(zp(k+1)-frheight))
        k = k + 1
       end if
      end do
      if ((WBT(k).le.(273.15)).and.
     +      (Low150WBT.le.(273.15))) then
       AWBT = AWBT +
     +        ((Low150WBT+WBT(k)-(2.*frWBT))*
     +        (0.5)*(Low150z-zp(k)))
       k = k + 1
      else if ((WBT(k).ge.(273.15)).and.
     +        (Low150WBT.ge.(273.15))) then
       AWBT = AWBT +
     +        (((Low150WBT+WBT(k)-(2.*frWBT))*(0.5)*
     +        (Low150z-zp(k))))
       k = k + 1
      else
       delWBT = Low150WBT - WBT(k)
       dellnp = Alog(Low150pres) - Alog(press(k))
       delfrWBT = 273.15 - WBT(k)
       delfrlnp = (delfrWBT)*(dellnp)/delWBT
       delz = Low150z - zp(k)
       frheight = zp(k) + (delz*delfrlnp/dellnp)
       frpress = alog(press(k)) + delfrlnp  
       frpress = exp(frpress)
       AWBT = AWBT +
     +        ((frWBT+WBT(k)-(2.*frWBT))*
     +        (0.5)*(frheight-zp(k)))
       AWBT = AWBT +
     +        ((Low150WBT+frWBT-(2.*frWBT))*
     +        (0.5)*(Low150z-frheight))
       k = k + 1
      end if
      AWBT2 = AWBT
        
*  Calculate area of surface based wet bulb temperature above 273.15 K
*  (AWBT3) and below 273.15 K (AWBT4)

c      print*,'Calculating area of surface based wet bulb temperature 
c     +above 273.15K...'        

      k = 1
      AWBT = 0.
      do while ((WBT(k).ge.(273.15)).and.((k+1).le.nz))
       if ((WBT(k).ge.(273.15)).and.(WBT(k+1).ge.(273.15))) then
        AWBT=AWBT +
     +       ((WBT(k+1)+WBT(k)-(2.*frWBT))*
     +       (0.5)*(zp(k+1)-zp(k)))
        k = k + 1
       else if ((WBT(k).gt.(273.15)).and.
     +          (WBT(k+1).lt.(273.15))) then
        delWBT=WBT(k+1) - WBT(k)
        dellnp=alog(press(k+1)) - alog(press(k))
        delfrWBT=273.15-WBT(k)
        delfrlnp=(delfrWBT)*(dellnp)/(delWBT)
        delz = zp(k+1) - zp(k)
        frheight = zp(k) + (delz*delfrlnp/dellnp)
        frpress=alog(Press(k)) + delfrlnp
        frpress = exp(frpress)
        AWBT=AWBT + ((frWBT+WBT(k)-(2.*frWBT))*
     +              (0.5)*(frheight-zp(k)))
        k = k + 1  
       end if
      end do
      AWBT3 = AWBT
       
c      print*,'Calculating area of surface based wet bulb temperature 
c     +below 273.15K...'

      k = 1
      AWBT = 0.
      do while ((WBT(k).le.(273.15)).and.(press(k+1).ge.Low150pres))
       if ((WBT(k).le.(273.15)).and.(WBT(k+1).le.(273.15))) then
        AWBT = AWBT +
     +         (((WBT(k+1)+WBT(k)-(2.*273.15)))*
     +         (0.5)*(zp(k+1)-zp(k)))
        k = k + 1
       else if ((WBT(k).lt.(273.15)).and.
     +          (WBT(k+1).gt.(273.15))) then   
        delWBT = WBT(k+1)-WBT(k)
        dellnp = alog(press(k+1)) - alog(press(k))
        delfrWBT = 273.15 - WBT(k)
        delfrlnp = (delfrWBT)*(dellnp)/(delWBT)
        delz = zp(k+1) - zp(k)
        frheight = zp(k) + (delz*delfrlnp/dellnp)
        frpress = alog(press(k))+ delfrlnp
        frpress=exp(frpress)
        AWBT = AWBT +
     +         ((frWBT+WBT(k)-(2.*273.15))*(0.5)*
     +         (frheight-zp(k)))
        k = k + 1
       end if
      end do
      if ((WBT(k).le.(273.15)).and.(press(k+1).lt.Low150pres).and.
     +      (Low150WBT.le.(273.15))) then
       AWBT = AWBT +
     +        (((Low150WBT+WBT(k)-(2.*273.15)))*
     +        (0.5)*(Low150z-zp(k)))
       k = k + 1
      else if ((WBT(k).lt.(273.15)).and.(press(k+1).lt.Low150pres).and.
     +          (Low150WBT.gt.(273.15))) then   
       delWBT = Low150WBT-WBT(k)
       dellnp = alog(Low150pres) - alog(press(k))
       delfrWBT = 273.15 - WBT(k)
       delfrlnp = (delfrWBT)*(dellnp)/(delWBT)
       delz = Low150z - zp(k)
       frheight = zp(k) + (delz*delfrlnp/dellnp)
       frpress = alog(press(k))+ delfrlnp
       frpress=exp(frpress)
       AWBT = AWBT +
     +        ((frWBT+WBT(k)-(2.*273.15))*(0.5)*
     +        (frheight-zp(k)))
       k = k + 1
      end if
      AWBT4 = AWBT
 
*  Calculate Lowest Layer Temperature (LowLayerT)
                
c      Print*,'Calculating Lowest Layer Temperature...'

      k = 1
      LowLayerT = 0.5*(TZ(k)+TZ(k+1))
      LowLayerT = LowLayerT + 273.15
       
*  Determine Precipitation type for each model column
C
C  Precipitation type (assuming precipitation occurs):
C                       1 -- Snow (SN)
C                       2 -- Ice Pellets (IP)
C                       3 -- Freezing Rain (FZ RN)
C                       4 -- Rain (RN)
      
       print*,'ColdSatT (K) = ',ColdSatT
       print*,'AWBT1 (K m) = ',AWBT1
       print*,'AWBT2 (K m) = ',AWBT2
       print*,'AWBT3 (K m) = ',AWBT3
       print*,'AWBT4 (K m) = ',AWBT4
       print*,'LowLayerT = ',lowlayert
c      print*,'Determining precipitation type...'
  
c      print*,'ColdsatT = ', coldsatT

      if (coldsatT.gt.(269.)) then
       if (LowLayerT.lt.(273.)) then
        meprectypeskewt = 3
       else
        meprectypeskewt = 4
       end if
      else if (AWBT1.lt.(3000.)) then
       meprectypeskewt = 1
      else if ((AWBT2.lt.(-3000.)).and.
     + (AWBT3.lt.(50.))) then
       meprectypeskewt = 2
      else if (AWBT4.lt.(-3000.))  then
       meprectypeskewt = 2
      else if (LowLayerT.lt.(273.)) then
       meprectypeskewt = 3
      else
       meprectypeskewt = 4
      end if
 
c      print*,'Cond. Precip. Type = ',meprectypeskewt

      return
      end  


