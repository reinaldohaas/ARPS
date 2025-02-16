!
!     ##################################################################
!     ##################################################################
!     ######                                                      ######
!     ######                 SUBROUTINE RDARPSQPF                 ######
!     ######                                                      ######
!     ######                     Developed by                     ######
!     ######     Center for Analysis and Prediction of Storms     ######
!     ######                University of Oklahoma                ######
!     ######                                                      ######
!     ##################################################################
!     ##################################################################

      SUBROUTINE RDARPSQPF(nx,ny,nz,tpcp,tpcpold,foundpcp,addpcp,nfile, &
     &                kgds,ibi,kbms,iyear,imon,iday,ihr,ifhr1, &
     &                ifhr2,jmaxin,ibufsize, &
     &                dir_extd,extdname,extdfcst,extdfmt,iimin,iimax, &
     &                jjmin,jjmax)

!#######################################################################
!
!     PURPOSE:
!     
!     Reads in ARPS history dump and calculates accumulated 
!     precipitation.
!     
!     AUTHOR:  Eric Kemp, February 2000.
!     Based on program ARPSENSCV.
!     
!#######################################################################
!     
!     Variable Declarations:
!
!#######################################################################

      IMPLICIT NONE

!#######################################################################
!
!     Include files:
!
!#######################################################################

      include 'indtflg.inc'
      include 'globcst.inc'
      include 'phycst.inc'
      include 'grid.inc'

!#######################################################################
!
!     External variables
!
!#######################################################################

      INTEGER nx,ny,nz
      INTEGER nzsoil
      INTEGER JMAXIN,IBUFSIZE
      REAL tpcp(JMAXIN)
      REAL tpcpold(JMAXIN)
      INTEGER foundpcp,addpcp   
      CHARACTER*132 NFILE
      INTEGER IYEAR,IMON,IDAY,IHR,IFHR1,IFHR2
      INTEGER KGDS(*)
      
      INTEGER IBI
      LOGICAL*1 KBMS(JMAXIN)
   
      CHARACTER*(*) dir_extd,extdname
      CHARACTER*9 extdfcst
      INTEGER extdfmt
      INTEGER iimin,iimax,jjmin,jjmax

!#######################################################################
!
!     Arrays to be read in:
!
!#######################################################################

      real x     (nx)        ! The x-coord. of the physical and
                             ! computational grid. Defined at u-point.  
      real y     (ny)        ! The y-coord. of the physical and
                             ! computational grid. Defined at v-point.
      real z     (nz)        ! The z-coord. of the computational grid.
                             ! Defined at w-point on the staggered grid.
      real zp    (nx,ny,nz)  ! The physical height coordinate defined at
                             ! w-point of the staggered grid.
      real zpsoil(nx,ny,nz)  ! The physical height coordinate of soil model
     
      real uprt  (nx,ny,nz)  ! Perturbation u-velocity (m/s)
      real vprt  (nx,ny,nz)  ! Perturbation v-velocity (m/s)
      real wprt  (nx,ny,nz)  ! Perturbation w-velocity (m/s)
      real pt    (nx,ny,nz)  ! Total poten
      real qvprt (nx,ny,nz)

      real u     (nx,ny,nz)  ! Total u-velocity (m/s)
      real v     (nx,ny,nz)  ! Total v-velocity (m/s)
      real w     (nx,ny,nz)  ! Total w-velocity (m/s)
      real ptprt (nx,ny,nz)  ! Perturbation potential temperature
                             ! from that of base state atmosphere (K)
      real pprt  (nx,ny,nz)  ! Perturbation pressure from that
                             ! of base state atmosphere (Pascal)
      real qv    (nx,ny,nz)  ! Water vapor specific humidity (kg/kg)
      real qc    (nx,ny,nz)  ! Cloud water mixing ratio (kg/kg)
      real qr    (nx,ny,nz)  ! Rain water mixing ratio (kg/kg)
      real qi    (nx,ny,nz)  ! Cloud ice mixing ratio (kg/kg)  
      real qs    (nx,ny,nz)  ! Snow mixing ratio (kg/kg)
      real qh    (nx,ny,nz)  ! Hail mixing ratio (kg/kg)
      real tke   (nx,ny,nz)  ! Turbulent Kinetic Energy ((m/s)**2)
      real kmh   (nx,ny,nz)  ! Horizontal turb. mixing coef. for
                             ! momentum. ( m**2/s )
      real kmv   (nx,ny,nz)  ! Vertical turb. mixing coef. for
                             ! momentum. ( m**2/s )  
      real ubar  (nx,ny,nz)  ! Base state u-velocity (m/s)
      real vbar  (nx,ny,nz)  ! Base state v-velocity (m/s)
      real wbar  (nx,ny,nz)  ! Base state w-velocity (m/s)
      real ptbar (nx,ny,nz)  ! Base state potential temperature (K)  
      real pbar  (nx,ny,nz)  ! Base state pressure (Pascal)   
      real rhobar(nx,ny,nz)  ! Base state density rhobar
      real qvbar (nx,ny,nz)  ! Base state water vapor specific humidity
                             ! (kg/kg)
      integer nstyps             ! Number of soil type
      parameter ( nstyps = 4 )
      
      integer soiltyp(nx,ny,nstyps)     ! Soil type
      real    stypfrct(nx,ny,nstyps)    ! Soil type fraction
      integer vegtyp (nx,ny) ! Vegetation type
      real    lai    (nx,ny) ! Leaf Area Index
      real    roufns (nx,ny) ! Surface roughness
      real    veg    (nx,ny) ! Vegetation fraction   
      
      real tsoil  (nx,ny,nzsoil,0:nstyps)    ! Deep soil temperature (K)
      real qsoil  (nx,ny,nzsoil,0:nstyps)    ! Deep soil temperature (K)
      real wetcanp(nx,ny,0:nstyps)    ! Canopy water amount
      real snowdpth(nx,ny)            ! Snow depth (m)
                             
      real raing(nx,ny)      ! Grid supersaturation rain
      real rainc(nx,ny)      ! Cumulus convective rain
      real raint(nx,ny)      ! Total rain (rainc+raing)
      
      real prcrate(nx,ny,4)     ! precipitation rate (kg/(m**2*s))
                                ! prcrate(1,1,1) = total precip. rate
                                ! prcrate(1,1,2) = grid scale precip. rate
                                ! prcrate(1,1,3) = cumulative precip. rate
                                ! prcrate(1,1,4) = microphysics precip.
                                !                  rate
      
      real radfrc(nx,ny,nz)     ! Radiation forcing (K/s)
      real radsw (nx,ny)        ! Solar radiation reaching the surface
      real rnflx (nx,ny)        ! Net radiation flux absorbed by surface
      REAL radswnet(nx,ny)   ! Net shortwave radiation
      REAL radlwin(nx,ny)    ! Incoming longwave radiation
      
      real usflx (nx,ny)        ! Surface flux of u-momentum (kg/(m*s**2))
      real vsflx (nx,ny)        ! Surface flux of v-momentum (kg/(m*s**2))
      real ptsflx(nx,ny)        ! Surface heat flux (K*kg/(m*s**2))
      real qvsflx(nx,ny)        ! Surface moisture flux (kg/(m**2*s))
      
!#######################################################################
!
!     Temporary work arrays
!
!#######################################################################

      real tem1(nx,ny,nz)
      real tem2(nx,ny,nz)
      real tem3(nx,ny,nz)

!#######################################################################
!
!     Other variables
!
!#######################################################################

      REAL accppt(nx-1,ny-1),tpptold(nx-1,ny-1)

      INTEGER nchin,lengbf,lenfil
      REAL time

      INTEGER i,j,ni,nj,nf,ireturn

      CHARACTER*80 timsnd
      INTEGER tmstrln
      REAL time_ext
      INTEGER ihr,imin,isec
      CHARACTER*3 fmtn
      INTEGER lenrun,ldir
      INTEGER hinfmt
      CHARACTER grdbasfn*132,filename*132

      INTEGER ICOMP
      CHARACTER*1 GDS(400),PDS(400)
      INTEGER KPTR(200),KPDS(200),IPDS(200),IGDS(200)
      DATA KPTR/200*0/,KPDS/200*0/
      INTEGER LENGDS,NPTS,IRET,KRET,k

      REAL TMPLAT

!#######################################################################
!
!     Temporary code
!
!#######################################################################

      REAL maxtpcp
      integer iproj
      real scale,trlon
      real latnot(2)
      real x0,y0
      real lat(nx,ny)
      real lon(nx,ny)
      real xctr,yctr,dx,dy
!     real x(nx)
!     real y(ny)
      real xsc(nx)
      real ysc(ny)


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!     Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

      maxtpcp = REAL(0)

      ni=nx-1
      nj=ny-1

      k = 1
      DO j = 1,ny
      DO i = 1,nx
        IF (i.lt.nx.AND.j.lt.ny) THEN 
          tpptold(i,j) = tpcpold(k) 
        ENDIF
        k = k + 1
      END DO
      END DO

      DO k = 1, JMAXIN
        tpcp(k) = REAL(0)
        tpcpold(k) = REAL(0)
      END DO

!#######################################################################
!
!     Build file names
!
!#######################################################################

      IF ( extdfcst .EQ. '         ') extdfcst='000:00:00'

      lenrun=LEN(dir_extd)   
      ldir=lenrun   
      CALL strlnth( dir_extd, ldir )
    
      IF ( ldir .EQ. 0 .OR. dir_extd(1:ldir) .EQ. ' ' ) THEN
        dir_extd = '.'
        ldir = 1
      END IF
 
      IF( dir_extd(ldir:ldir) .NE. '/' .AND.  ldir .LT. lenrun ) THEN
        ldir = ldir + 1
        dir_extd(ldir:ldir) = '/'
      END IF

      lenrun = LEN( extdname )
      CALL strlnth( extdname, lenrun )

      IF( extdfmt .EQ. 1 ) THEN
        fmtn = 'bin'
      ELSE IF ( extdfmt .EQ. 2 ) THEN 
        fmtn = 'asc'
      ELSE IF ( extdfmt .EQ. 3 ) THEN
        fmtn = 'hdf'  
      ELSE IF ( extdfmt .EQ. 4 ) THEN
        fmtn = 'pak'
      ELSE IF ( extdfmt .EQ. 6 ) THEN
        fmtn = 'bn2'
      ELSE IF ( extdfmt .EQ. 7 ) THEN
        fmtn = 'net'
      ELSE IF ( extdfmt .EQ. 8 ) THEN
        fmtn = 'npk'
      ELSE IF ( extdfmt .EQ. 9 ) THEN
        fmtn = 'gad'
      ELSE IF ( extdfmt .EQ. 10 ) THEN
        fmtn = 'grb'
      ELSE
        WRITE(6,'(a,a,a)') &
     & 'Unknown format, ', extdfmt, '. Program stopped in RDARPSQPF.'
        STOP
      END IF

      READ(extdfcst,'(i3,1x,i2,1x,i2)') ihr,imin,isec
 
      time_ext = FLOAT( (ihr*3600)+(imin*60)+isec )
      CALL cvttsnd( time_ext, timsnd, tmstrln )
    
      grdbasfn = dir_extd(1:ldir)//extdname(1:lenrun) &
     &               //'.'//fmtn//'grdbas'
      lenfil = ldir + lenrun + 10
    
      filename = dir_extd(1:ldir)//extdname(1:lenrun) &
     &             //'.'//fmtn//timsnd(1:tmstrln)
      lengbf = ldir + lenrun + 4 + tmstrln
    
      WRITE(6,*) 'The external grid and base file, grdbasfn = ', &
     &           grdbasfn(1:lengbf)

      WRITE(6,*) 'The external time dependent file, filename = ', &
     &           filename(1:lenfil)

!#######################################################################
!
!     Read ARPS history dump.
!
!#######################################################################

      lengbf = 80
      CALL strlnth( grdbasfn, lengbf)
      lenfil = 80
      CALL strlnth( filename, lenfil)
      hinfmt = extdfmt

      CALL dtaread(nx,ny,nz,nzsoil, nstyps,  &
     &            hinfmt, nchin,grdbasfn(1:lengbf),lengbf, &
     &            filename(1:lenfil),lenfil,time, &
     &            x,y,z,zp,zpsoil,uprt ,vprt ,wprt ,ptprt, pprt , &
     &            qvprt, qc, qr, qi, qs, qh, tke,kmh,kmv, &
     &            ubar, vbar, wbar, ptbar, pbar, rhobar, qvbar, &
     &            soiltyp,stypfrct,vegtyp,lai,roufns,veg, &
     &            tsoil,qsoil,wetcanp,snowdpth, &
     &            raing,rainc,prcrate, &
     &            radfrc,radsw,rnflx, &
     &            radswnet,radlwin, &
     &            usflx,vsflx,ptsflx,qvsflx, &
     &            ireturn, tem1,tem2, tem3)

      IYEAR = year
      IMON = month
      IDAY = day
      IHR = hour
      IFHR1 = 0
      IFHR2 = time/3600

      IF (trulat1.gt.trulat2) THEN
        tmplat = trulat1
        trulat1 = trulat2
        trulat2 = tmplat
      ENDIF

!#######################################################################
!
!     Create decoded Grid Description Section data (array KGDS) for 
!     subroutine ipolates.
!
!#######################################################################

      CALL mkigds(nx,ny,nz,0,igds)
      ICOMP=MOD(IGDS(8)/8,2)                      
      CALL W3FI74(IGDS,ICOMP,GDS,LENGDS,NPTS,IRET) 
      IF(IRET.NE.0) THEN
        WRITE(6,*)'RDARPSQPF: ERROR -- Could not create GDS.'
        WRITE(6,*)'IRET = ',IRET
        WRITE(6,*)'Aborting...'
        STOP
      ENDIF
      CALL FI633(GDS,KPTR,KGDS,IRET)
      IF(IRET.NE.0) THEN
        WRITE(6,*)'RDARPSQPF: ERROR -- Could not create KGDS.'
        WRITE(6,*)'IRET = ',IRET
        WRITE(6,*)'Aborting...'
        STOP
      ENDIF
      KGDS(14) = 1 ! Flag for subroutine gdswiz03
      KGDS(15) = 0
      KGDS(19) = 0 

!#######################################################################
!
!     Create bitmap for subroutine ipolates
!
!#######################################################################

      IBI = 1 ! Flag indicating that bitmap will be used.

      DO i = 1,JMAXIN
        KBMS(i) = .FALSE.
      END DO
     
      k = 1
      DO j = 1,ny
      DO i = 1,nx
        IF (i.ge.iimin .AND. i.le.iimax .AND. &
     &      j.ge.jjmin .AND. j.le.jjmax) THEN
          KBMS(k) = .TRUE.
        ENDIF
        k = k + 1
      END DO
      END DO
           
!#######################################################################
!
!     Convert the rainfall arrays and store the old ones.
!
!#######################################################################

      DO 300 j=1,ny-1
      DO 300 i=1,nx-1
        raint(i,j)=raing(i,j)+rainc(i,j)
300   CONTINUE
      CALL ARY2DCV(nx,ny,raint,ni,nj,accppt)

!#######################################################################
!
!     Convert accppt from 0h-now accul. rain to (Tnf-1->Tnf) accul. rain
!     (by subtracting tpptold) and store the original value in tpptold
!
!#######################################################################

      CALL pptsto(ni,nj,accppt,tpptold)

      k = 1
      DO j = 1,ny
      DO i = 1,nx
        IF (i.lt.nx.AND.j.lt.ny) THEN 
          tpcp(k) = accppt(i,j)
          tpcpold(k) = tpptold(i,j)
          maxtpcp = MAX(maxtpcp,tpcp(k))
        ENDIF
        k = k + 1
      END DO
      END DO

!#######################################################################
!
!     Set foundpcp to one (indicating that QPF has been found) and 
!     return.
!
!#######################################################################

      foundpcp = 1

      RETURN
      END

!
!     ##################################################################
!     ##################################################################
!     ######                                                      ######
!     ######                  SUBROUTINE ARY2DCV                  ######
!     ######                                                      ######
!     ######                     Developed by                     ######
!     ######     Center for Analysis and Prediction of Storms     ######
!     ######                University of Oklahoma                ######
!     ######                                                      ######
!     ##################################################################
!     ##################################################################

      SUBROUTINE ARY2DCV(nx,ny,a,ni,nj,b)   

!#######################################################################
!
!     PURPOSE:
!     
!     Copies the contents of one 2-D array into another.
!     
!     AUTHOR:  ?????
!
!     MODIFICATION HISTORY:
!     Eric Kemp, February 2000
!     Added Documentation.
!     
!#######################################################################
!     
!     Variable Declarations:
!
!#######################################################################

      IMPLICIT none

      INTEGER nx,ny,ni,nj
      REAL a(nx,ny),b(ni,nj)
      INTEGER i,j

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!     Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!     print *,a(35,35),a(78,68)
      DO j=1,nj
      DO i=1,ni  
        b(i,j)=a(i,j)
      ENDDO
      ENDDO

      RETURN  
      END  

