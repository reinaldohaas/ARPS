
!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINE WRTTILTS                  ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
  SUBROUTINE wrttilts(rfname,maxgatevel,maxgateref,maxazim,maxelev,       &
                      rngvvol,azmvvol,elvvvol,velvol,                     &
                      rngrvol,azmrvol,elvrvol,refvol,                     &
                      radar_alt,radar_lat,radar_lon,                      &
                      iyear,imon,iday,ihour,imin,isec )
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  write out radar observation in  radar coordinate
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Hu, CAPS
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    maxgate   Maximum gates in a radial
!    maxazim   Maximum radials per tilt
!    maxelev   Maximum number of tilts
!
  IMPLICIT NONE
!
  INTEGER, INTENT(IN) :: maxgatevel
  INTEGER, INTENT(IN) :: maxgateref
  INTEGER, INTENT(IN) :: maxazim
  INTEGER, INTENT(IN) :: maxelev

  REAL, INTENT(IN)    :: rngvvol(maxgatevel,maxelev)
  REAL, INTENT(IN)    :: azmvvol(maxazim,maxelev)
  REAL, INTENT(IN)    :: elvvvol(maxazim,maxelev)
  REAL, INTENT(IN)    :: velvol(maxgatevel,maxazim,maxelev)

  REAL, INTENT(IN)    :: rngrvol(maxgateref,maxelev)
  REAL, INTENT(IN)    :: azmrvol(maxazim,maxelev)
  REAL, INTENT(IN)    :: elvrvol(maxazim,maxelev)
  REAL, INTENT(IN)    :: refvol(maxgateref,maxazim,maxelev)

  INTEGER, INTENT(IN) :: iyear,imon,iday,ihour,imin,isec
  REAL, INTENT(IN)    :: radar_alt
  REAL, INTENT(IN)    :: radar_lat
  REAL, INTENT(IN)    :: radar_lon

  CHARACTER(LEN=256)  :: rfname
  CHARACTER(LEN=256)  :: rfnametilt
  INTEGER             :: ilen

!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code.
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  ilen=len(trim(rfname))
  rfnametilt=rfname(1:ilen-1)//'.tilts'
  WRITE(6,'(1x,2a)') 'Tilt data file name: ',TRIM(rfnametilt)
  OPEN(13,FILE=trim(rfnametilt),STATUS='unknown',                          &
       FORM='unformatted')
  write(13) maxgatevel,maxgateref,maxazim,maxelev
  write(13) radar_alt,radar_lat,radar_lon
  write(13) iyear,imon,iday,ihour,imin,isec
  write(13) rngrvol
  write(13) azmrvol
  write(13) elvrvol
  write(13) refvol

  write(13) rngvvol
  write(13) azmvvol
  write(13) elvvvol
  write(13) velvol

  close(13)

  RETURN
  END SUBROUTINE wrttilts

!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINE REMAP2DCTS                ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE remap2DCts(maxgate,maxazim,maxelev,nx,ny,nz,nzsnd,           &
                    kntazim,kntelev,                                    &
                    rdrlat,rdrlon,radarx,radary,rdralt,dazim,           &
                    rngmin,rngmax,time1st,                              &
                    rngvol,azmvol,elvvol,                               &
                    varvol,timevol,                                     &
                    xs,ys,zps,zsnd,rfrsnd,                              &
                    gridtilthigh,gridrange,gridslr,gridazm,             &
                    gridtiltval,gridtilttime,elevmean,n_gate_elev,istatus)
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Take data from a radar volume in polar coordinates and map it to
!  ARPS Cartesian terrain-following grid, but don't do vertical interoplation
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Hu, CAPS
!
!      Oct. 13, 2004
!
!  MODIFICATION HISTORY:
!     how to calculate distance between two points on earth
!     how to get radar azimuth angle, range and height
!     how to interpolation on radar tilt
!  22 June 2011 Keith Brewster
!               Corrected problem that occurs if there is an azimuth gap
!               in the data.  Program will not try to interpolate across
!               a large gap (> parameter azmgapmax).
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    maxgate   Maximum gates in a radial
!    maxazim   Maximum radials per tilt
!    maxelev   Maximum number of tilts
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!    varid    Radar variable ID for diagnostic file writing
!    varname  Name of radar variable for diagnostic file writing
!    varunit  Units of radar variable for diagnostic file writing
!
!    varchek  Threshold for checking data, good vs. flagged
!    varmiss  Value to assign to data for missing
!    vmedlim  Threshold limit for median check
!    dazlim   Maximum value of azimuth difference (grid vs data) to accept
!             Generally should be 30 degrees or less for velocity, 360 for refl
!    rngmin   Minimum range (m) of data to use
!            (10 000 m or more to eliminate near field ground targets).
!    rngmax   Maximum range (m) of data to use
!
!    rngvvol  Range to gate in velocity 3-D volume
!    azmvvol  Azimuth angle in velocity 3-D volume
!    elvvvol  Elevation angle in velocity 3-D volume
!    varvol   Radar data 3-D volume
!
!    xs       x coordinate of scalar grid points in physical/comp. space (m)
!    ys       y coordinate of scalar grid points in physical/comp. space (m)
!    zps      Vertical coordinate of scalar grid points in physical space(m)
!
!  OUTPUT:
!
!  REAL, INTENT(OUT)   :: gridtilthigh(nx,ny,maxelev)
!  REAL, INTENT(OUT)   :: gridrange(nx,ny,maxelev)
!  REAL, INTENT(OUT)   :: gridslr(nx,ny)
!  REAL, INTENT(OUT)   :: gridazm(nx,ny)
!
!  REAL, INTENT(OUT)   :: gridtiltval(nx,ny,maxelev)
!  INTEGER, INTENT(OUT)   :: gridtilttime(maxelev)
!
!  REAL, INTENT(OUT) :: elevmean(maxelev)
!
!    istatus  Status indicator
!
!-----------------------------------------------------------------------
!

!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
  INTEGER, INTENT(IN) :: maxgate
  INTEGER, INTENT(IN) :: maxazim
  INTEGER, INTENT(IN) :: maxelev
  INTEGER, INTENT(IN) :: nx
  INTEGER, INTENT(IN) :: ny
  INTEGER, INTENT(IN) :: nz
  INTEGER, INTENT(IN) :: nzsnd

!  REAL, INTENT(IN)    :: varchek
!  REAL, INTENT(IN)    :: varmiss
!  REAL, INTENT(IN)    :: vmedlim
!  REAL, INTENT(IN)    :: dazlim
!  INTEGER, INTENT(IN) :: iorder

  INTEGER, INTENT(IN) :: kntazim(maxelev)
  INTEGER, INTENT(IN) :: kntelev

  REAL, INTENT(IN)    :: rdrlat
  REAL, INTENT(IN)    :: rdrlon
  REAL, INTENT(IN)    :: radarx
  REAL, INTENT(IN)    :: radary
  REAL, INTENT(IN)    :: rdralt
  REAL, INTENT(IN)    :: dazim
  REAL, INTENT(IN)    :: rngmin
  REAL, INTENT(IN)    :: rngmax

  INTEGER, INTENT(IN) :: time1st
  INTEGER, INTENT(IN) :: timevol(maxazim,maxelev)
  REAL, INTENT(IN)    :: rngvol(maxgate,maxelev)
  REAL, INTENT(IN)    :: azmvol(maxazim,maxelev)
  REAL, INTENT(IN)    :: elvvol(maxazim,maxelev)
  REAL, INTENT(IN)    :: varvol(maxgate,maxazim,maxelev)

  REAL, INTENT(IN)    :: xs(nx)
  REAL, INTENT(IN)    :: ys(ny)
  REAL, INTENT(IN)    :: zps(nx,ny,nz)
  REAL, INTENT(IN)    :: zsnd(nzsnd)
  REAL, INTENT(IN)    :: rfrsnd(nzsnd)

  REAL, INTENT(OUT)   :: gridtilthigh(nx,ny,maxelev)
  REAL, INTENT(OUT)   :: gridrange(nx,ny,maxelev)
  REAL, INTENT(OUT)   :: gridslr(nx,ny)
  REAL, INTENT(OUT)   :: gridazm(nx,ny)

  REAL, INTENT(OUT)   :: gridtiltval(nx,ny,maxelev)
  INTEGER, INTENT(OUT):: gridtilttime(maxelev)

  REAL, INTENT(OUT) :: elevmean(maxelev)
  INTEGER, INTENT(IN) :: n_gate_elev(30)
  INTEGER, INTENT(OUT) :: istatus

!
!-----------------------------------------------------------------------
!
! Misc. Local Variables
!
!-----------------------------------------------------------------------
!
  REAL, PARAMETER :: azmgapmax = 5.0
  INTEGER :: iyr,imon,iday,ihour,imin,isec
  REAL    :: xsll(nx,ny),ysll(nx,ny)

  INTEGER :: i,j,k,knt,kinbox
  INTEGER :: kok,isort,jsort,mid
  INTEGER :: kbgn,kend
  INTEGER :: igate,jazim,kelev,jazmin,jmirror,jend
  INTEGER :: istatal,istatwrt

  INTEGER :: nweight
  REAL :: val, weight

  REAL :: deg2rad,rad2deg
  REAL :: elevmin,elevmax
  REAL :: delx,dely,delz,dazimr,daz,azdiff
  REAL :: ddx,ddxy,ddx2,ddy,ddy2,ddz,dxthr,dxthr0
  REAL :: mapfct,xcomp,ycomp,azmrot,sfcr,zagl
  REAL :: sum,sum2,sdev,thresh,slrange,elijk,azimijk,time
  REAL :: varmax,varmin,varavg,varmean,varmed

  REAL :: kea, skea,elevE
!  REAL :: eradius=6371000
  REAL :: missvalue = -99999.9
  REAL :: azimj1,azimj2,delazim,dazimj,drange,dgate
  INTEGER :: ii,jj
  REAL :: valrd,valru,valld,vallu
  REAL :: rngminMulti(20),rngmink
  REAL :: xylat,xylon
  REAL :: latnot(2)
  REAL :: swx,swy, ctrx, ctry
  INTEGER :: n_gate_elev_new(30),gate_count

  INTEGER :: nxlg, nylg
!
!-----------------------------------------------------------------------
!
! Include files:
!
!-----------------------------------------------------------------------

  INCLUDE 'globcst.inc'
  INCLUDE 'grid.inc'
  INCLUDE 'mp.inc'

!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  rngminMulti=10
  rngmink=rngmin
  rngmink=10

  ! Reorganize n_gate_elev
  gate_count=0
  DO k=1,30
    if (n_gate_elev(k) > 1) then
        gate_count=gate_count+1
        n_gate_elev_new(gate_count)=n_gate_elev(k)
    endif
  ENDDO

  deg2rad = atan(1.)/45.
  rad2deg = 1./deg2rad
  dazimr = dazim*deg2rad
  dxthr0=0.6*max((xs(3)-xs(2)),(ys(3)-ys(2)))

  CALL lattomf(1,1,rdrlat,mapfct)

  time=0.
  elevmax=-999.
  elevmin=999.
  DO kelev=1,kntelev
    sum=0.
    DO jazim=1,kntazim(kelev)
      sum=sum+elvvol(jazim,kelev)
      elevmin=min(elevmin,elvvol(jazim,kelev))
      elevmax=max(elevmax,elvvol(jazim,kelev))
    END DO
    elevmean(kelev)=sum/float(kntazim(kelev))
  END DO

   DO kelev=1,kntelev
      jazmin=kntazim(kelev)/2
      gridtilttime(kelev)=timevol(jazmin,kelev)
      if( gridtilttime(kelev) < 0 ) Then
         write(*,*) ' WRONG: radar observation time has problem !!!'
         stop 222
      endif
!      write(*,*) 'time if observation for tilt',kelev,' is ',gridtilttime(kelev)
   ENDDO
!
! get the latitude and longitude of each grid point
!
  latnot(1)=trulat1
  latnot(2)=trulat2
  CALL setmapr(mapproj,sclfct,latnot,trulon)

  CALL lltoxy( 1,1, ctrlat,ctrlon, ctrx, ctry )
  nxlg = (nx-3)*nproc_x+3
  nylg = (ny-3)*nproc_y+3
  swx = ctrx - (FLOAT((nxlg-3))/2.) * dx
  swy = ctry - (FLOAT((nylg-3))/2.) * dy
  CALL setorig( 1, swx, swy)

!  CALL xytoll(1,1,radarx,radary,radarlat,radarlon)
!  CALL lltoxy(1,1,rdrlat,rdrlon,radarx,radary)
  gridslr = missvalue
  gridazm = missvalue
  DO j=1,ny-1
    DO i=1,nx-1
      CALL xytoll(1,1,xs(i),ys(j),xylat,xylon)
      CALL disthead(rdrlat,rdrlon,xylat,xylon,azimijk,slrange)
      IF( slrange >= rngmink .and. slrange<=rngmax ) THEN
          gridslr(i,j) = slrange
          gridazm(i,j) = azimijk
      ENDIF  ! slrange >= 10000.0 .and. slrange<=230000.0
    ENDDO
  ENDDO
!
  IF (myproc == 0) write(*,'(1x,a,F12.3)') ' Radar altitude = ',rdralt
  kea=4.0/3.0*eradius
  gridtilthigh = missvalue
  gridrange   = missvalue
  DO j=1,ny-1
    DO i=1,nx-1
      !
      ! get gridtilthigh and gridrange
      !
      slrange=gridslr(i,j)
      IF( slrange >= 0.0 ) THEN
        skea=slrange/kea
        DO k=1,kntelev
           rngmink=rngminMulti(k)
           elevE=elevmean(k)*deg2rad

           ! according to 2.28a
           gridtilthigh(i,j,k)=kea*(cos(elevE)/cos(elevE+skea)-1.) + rdralt
           ! according to 2.28c
           gridrange(i,j,k)=(kea+gridtilthigh(i,j,k))*sin(skea)/cos(elevE)

           IF( gridtilthigh(i,j,k) > 20000.0 .or.                                &
               gridrange(i,j,k) < rngmink .or. gridrange(i,j,k) > rngmax ) THEN
             gridtilthigh(i,j,k)=  missvalue
             gridrange(i,j,k) =  missvalue
           END IF
        END DO
      END IF ! slrange >= 0.0
    END DO
  END DO

!
!  interpolation
!
  gridtiltval=missvalue
  DO k=1,kntelev
    dgate=rngvol(3,k)-rngvol(2,k)
    !write(*,*) k
    DO j=1,ny-1
      DO i=1,nx-1
         azimijk=gridazm(i,j)
         slrange=gridrange(i,j,k)

         if( azimijk > -88888.0 .AND. slrange > -88888.0 ) THEN
           ii=-99999
           jj=-99999
           DO jazim=1,kntazim(k)-1
             if( jazim == 1 ) then
                azimj1=azmvol(jazim,k)
                if(azimj1 > azimijk ) azimijk=azimijk+360.0
             else
                azimj1=azimj2
             endif
             azimj2=azmvol(jazim+1,k)
             if ( azimj1 > azimj2 ) azimj2 = azimj2 + 360.0
             if( azimijk >= azimj1 .and. azimijk < azimj2 ) then
               delazim = azimj2 - azimj1
               IF( delazim < azmgapmax ) THEN
                 ii=jazim
                 dazimj=(azimijk-azimj1)/delazim
               END IF
               EXIT
             end if
           END DO

           IF( ii > -8888 ) THEN
             jj=int(slrange/dgate)+1
             drange=mod(slrange,dgate)/dgate

             IF( ii<=maxazim-1 .and. jj<=n_gate_elev_new(k)-1 ) THEN

                valld=varvol(jj,ii,k)
                vallu=varvol(jj+1,ii,k)
                valrd=varvol(jj,ii+1,k)
                valru=varvol(jj+1,ii+1,k)
                if(valld > -488.0 .and. vallu > -488.0 .and.                &
                   valrd > -488.0 .and. valru > -488.0 ) then
                   gridtiltval(i,j,k)=valld*(1.-dazimj)*(1.-drange) +         &
                        vallu*(1.-dazimj)*drange + valrd*dazimj*(1.-drange)+  &
                        valru*dazimj*drange
                else
                  nweight=0
                  val=0
                  weight=0
                  if(valld > -488.0) then
                     nweight = nweight + 1
                     val= val + valld * (1.-dazimj)*(1.-drange)
                     weight = weight + (1.-dazimj)*(1.-drange)
                  endif
                  if(vallu > -488.0) then
                     nweight = nweight + 1
                     val= val + vallu * (1.-dazimj)*drange
                     weight = weight + (1.-dazimj)*drange
                  endif
                  if(valrd > -488.0) then
                     nweight = nweight + 1
                     val= val + valrd * dazimj*(1.-drange)
                     weight = weight + dazimj*(1.-drange)
                  endif
                  if(valru > -488.0) then
                     nweight = nweight + 1
                     val= val + valru * dazimj*drange
                     weight = weight + dazimj*drange
                  endif

                  if (nweight > 1 ) then
                    IF (abs(weight) > 0) THEN
                      gridtiltval(i,j,k)=val/weight
                    ELSE
                      gridtiltval(i,j,k) = val
                    END IF
                  else
                    gridtiltval(i,j,k)=-77777.7
                  endif
                END IF
             END IF
           END IF
         END IF
      END DO
    END DO
  END DO

  RETURN
  END SUBROUTINE remap2DCts
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE WRTGRIDTILT                ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE wrtgridtilt_old(rfname,maxelev,nx,ny,kntrelv,radar,              &
                rdrlat,rdrlon,radarx,radary,rdralt,dazim,rngmin,rngmax, &
                gridtilthighref,gridrangeref,gridslrref,gridazmref,     &
                gridtiltref,tilttimeref,tilttimevel,elevmeanref,        &
                gridtilthighvel,gridrangevel,gridslrvel,gridazmvel,     &
                gridtiltvel,elevmeanvel)

!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!       Write out data in radar elevation level and horizontal grid column
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Hu, CAPS
!
!  MODIFICATION HISTORY:
!
!  Modifid by Jili Dong to output EnKF format (2008/11/04)
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    maxelev   Maximum number of tilts
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    maxelev   Maximum number of tilts
!    kntrelv  Number of tilts
!
!    rdrlat   : latitude of radar location
!    rdrlon   : longitude of  radar location
!    radarx   : x position of radar location
!    radary   : y position of radar location
!    rdralt   : elevation of radar
!    dazim    : average azimuth angle
!    rngmin   Minimum range (m) of data to use
!            (10 000 m or more to eliminate near field ground targets).
!    rngmax   Maximum range (m) of data to use
!
!
!   gridtilthigh(nx,ny,maxelev)     high of  observation point
!   gridrange(nx,ny,maxelev)        range of observation point
!   gridslr(nx,ny)                  surface distance of observation point
!   gridazm(nx,ny)                  azimuth angle
!
!   gridtiltval(nx,ny,maxelev)      observation value
!   gridtilttime(maxelev)     onservation time of each tilt
!
!   elevmean(maxelev)               mean elevation angle
!
!    istatus  Status indicator
!
!-----------------------------------------------------------------------
!
!
!  Note: tilttimeref is uncorrect for time is not readed in at all for reflectivity.
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
  INTEGER, INTENT(IN) :: maxelev
  INTEGER, INTENT(IN) :: nx
  INTEGER, INTENT(IN) :: ny
  INTEGER, INTENT(IN) :: kntrelv

  REAL, INTENT(IN)    :: rdrlat
  REAL, INTENT(IN)    :: rdrlon
  REAL, INTENT(IN)    :: radarx
  REAL, INTENT(IN)    :: radary
  REAL, INTENT(IN)    :: rdralt
  REAL, INTENT(IN)    :: dazim
  REAL, INTENT(IN)    :: rngmin
  REAL, INTENT(IN)    :: rngmax

  REAL, INTENT(IN)   :: gridtilthighref(nx,ny,maxelev)
  REAL, INTENT(IN)   :: gridrangeref(nx,ny,maxelev)
  REAL, INTENT(IN)   :: gridslrref(nx,ny)
  REAL, INTENT(IN)   :: gridazmref(nx,ny)

  REAL, INTENT(INOUT):: gridtiltref(nx,ny,maxelev)
  REAL, INTENT(IN)   :: elevmeanref(maxelev)

  REAL, INTENT(IN)   :: gridtilthighvel(nx,ny,maxelev)
  REAL, INTENT(IN)   :: gridrangevel(nx,ny,maxelev)
  REAL, INTENT(IN)   :: gridslrvel(nx,ny)
  REAL, INTENT(IN)   :: gridazmvel(nx,ny)

  REAL, INTENT(INOUT):: gridtiltvel(nx,ny,maxelev)
  REAL, INTENT(IN)   :: elevmeanvel(maxelev)

  INTEGER, INTENT(IN)   :: tilttimeref(maxelev)
  INTEGER, INTENT(IN)   :: tilttimevel(maxelev)

  REAL   :: refsngl(nx,ny)
  REAL   :: velsngl(nx,ny)
  REAL   :: elevref
  REAL   :: elevvel

  CHARACTER (LEN=4) :: radar

  CHARACTER (LEN=100) :: rfname
  CHARACTER (LEN=100) :: rfnametilt
  integer :: ilen
  integer :: numscan
  integer :: numelev

  REAL, ALLOCATABLE :: elevmean_recount_vel(:), elevmean_recount_ref(:)
  REAL, ALLOCATABLE :: gridrangevel_recount(:,:,:), gridrangeref_recount(:,:,:)
  REAL, ALLOCATABLE :: gridtilthighvel_recount(:,:,:), gridtilthighref_recount(:,:,:)
  REAL, ALLOCATABLE :: gridtiltvel_recount(:,:,:), gridtiltref_recount(:,:,:)

  INTEGER :: count_elev_vel, count_elev_ref

!
!-----------------------------------------------------------------------
!
! Misc. Local Variables
!
!-----------------------------------------------------------------------
!

  INTEGER :: i,j,k
  INTEGER :: iyr,imon,iday,ihour,imin,isec,time1st
  INTEGER :: samefirstlevel,kntrelv14
  INTEGER :: istat

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  radar='KTLX'
  samefirstlevel=0
  numscan=7

  IF( 1 == 2) then
    write(*,*) elevmeanvel
    write(*,*) elevmeanref
    write(*,*) ' Mean elevation angle and time of each tilt'
    if(abs(elevmeanref(1)-elevmeanvel(1)) < 0.2 ) then
       samefirstlevel=1
       kntrelv14 = kntrelv
    else
       kntrelv14 = kntrelv+1
    end if
    DO k=1, kntrelv14
        refsngl=-9999.9
        velsngl=-9999.9
        elevref=-9999.9
        elevvel=-9999.9
        time1st=tilttimevel(k)
        call abss2ctim(time1st,iyr,imon,iday,ihour,imin,isec)
        write(*,*) iyr,imon,iday,ihour,imin,isec

  !      WRITE(rfnametilt,'(a4,a1,i4.4,2(i2.2),a1,3(i2.2),a9)')           &
  !        radar,'.', iyr,imon,iday,'.',ihour,imin,isec,'.gridtilt'

        WRITE(rfnametilt,'(a4,a1,i4.4,2(i2.2),a1,2(i2.2),a9)')          &
          radar,'.', iyr,imon,iday,'.',numscan,k,'.gridtilt'

         elevvel=elevmeanvel(k)
         DO j=1,ny
         DO i=1,nx
           velsngl(i,j)=gridtiltvel(i,j,k)
         END DO
         END DO

         if(samefirstlevel == 1 ) then
             elevref=elevmeanref(k)
             DO j=1,ny
             DO i=1,nx
               refsngl(i,j)=gridtiltref(i,j,k)
             ENDDO
             ENDDO
          else
           if(k>=2) then
             elevref=elevmeanref(k-1)
             DO j=1,ny
             DO i=1,nx
               refsngl(i,j)=gridtiltref(i,j,k-1)
             ENDDO
             ENDDO
           endif
         endif
         write(*,*) 'write out tilt=',k,elevref,elevvel,time1st
         if( (elevref > -8888.8 .and. elevvel > -8888.8 ) .and.         &
               abs(elevref-elevvel) > 0.2 ) then
            write(*,*) ' Reflectivity is not at the same tilt as radial velocity'
            stop 12345
         endif

        OPEN(55,FILE=trim(rfnametilt),STATUS='unknown',                 &
                FORM='unformatted')
            write(55) time1st,iyr,imon,iday,ihour,imin,isec
        write(55) k,nx,ny
        write(55) rdrlat,rdrlon,radarx,radary,rdralt
        write(55) dazim,rngmin,rngmax

        write(55) elevvel
        write(55) ((gridtilthighvel(i,j,k),i=1,nx),j=1,ny)
        write(55) ((gridrangevel(i,j,k),i=1,nx),j=1,ny)
        write(55) gridslrvel
        write(55) gridazmvel
        write(55) velsngl

            write(55) elevref
            write(55) refsngl
       CLOSE(55)
    END DO
  END IF

! Modified by Jili Dong
! The following is modified to write the data for EnKF without doing time interpolation
! If you want to use time interpolation, pls use the backup f90

  IF( .TRUE. ) then
    time1st=tilttimevel(1)
    CALL abss2ctim(time1st,iyr,imon,iday,ihour,imin,isec)
    WRITE(*,*) iyr,imon,iday,ihour,imin,isec

    WRITE(rfnametilt,'(a4,a1,i4.4,2(i2.2),a1,3(i2.2),a9)')           &
       radar,'.', iyr,imon,iday,'.',ihour,imin,isec,'.gridtilt'

    ! Redefine tilt number,remove redundant scan
    numelev=maxelev
    DO i=1,maxelev
       if(elevmeanref(i)==0.0) numelev=numelev-1
    END DO

    ALLOCATE( elevmean_recount_vel(maxelev), STAT = istat )
    ALLOCATE( elevmean_recount_ref(maxelev), STAT = istat )
    ALLOCATE( gridrangevel_recount(nx,ny,maxelev), STAT = istat )
    ALLOCATE( gridrangeref_recount(nx,ny,maxelev), STAT = istat )
    ALLOCATE( gridtilthighvel_recount(nx,ny,maxelev), STAT = istat )
    ALLOCATE( gridtilthighref_recount(nx,ny,maxelev), STAT = istat )
    ALLOCATE( gridtiltvel_recount(nx,ny,maxelev), STAT = istat )
    ALLOCATE( gridtiltref_recount(nx,ny,maxelev), STAT = istat )

    count_elev_ref=1
    count_elev_vel=1
    elevmean_recount_ref(count_elev_ref)=elevmeanref(1)
    gridtilthighref_recount(:,:,count_elev_ref)=gridtilthighref(:,:,1)
    gridtiltref_recount(:,:,count_elev_ref)=gridtiltref(:,:,1)
    gridrangeref_recount(:,:,count_elev_ref)=gridrangeref(:,:,1)

    elevmean_recount_vel(count_elev_vel)=elevmeanvel(1)
    gridtilthighvel_recount(:,:,count_elev_vel)=gridtilthighvel(:,:,1)
    gridtiltvel_recount(:,:,count_elev_vel)=gridtiltvel(:,:,1)
    gridrangevel_recount(:,:,count_elev_vel)=gridrangevel(:,:,1)

    DO i=2,numelev

       IF((elevmeanref(i)-elevmeanref(i-1))> 0.2) THEN
         count_elev_ref=count_elev_ref+1
         elevmean_recount_ref(count_elev_ref)=elevmeanref(i)
         gridtilthighref_recount(:,:,count_elev_ref)=gridtilthighref(:,:,i)
         gridtiltref_recount(:,:,count_elev_ref)=gridtiltref(:,:,i)
         gridrangeref_recount(:,:,count_elev_ref)=gridrangeref(:,:,i)
       END IF

       IF((elevmeanvel(i)-elevmeanvel(i-1))> 0.2) THEN
         count_elev_vel=count_elev_vel+1
         elevmean_recount_vel(count_elev_vel)=elevmeanvel(i)
         gridtilthighvel_recount(:,:,count_elev_vel)=gridtilthighvel(:,:,i)
         gridtiltvel_recount(:,:,count_elev_vel)=gridtiltvel(:,:,i)
         gridrangevel_recount(:,:,count_elev_vel)=gridrangevel(:,:,i)
       END IF

    END DO
    ! End of removing redundant scans

    ! Set clear air missing Z to 0
    ! Set Z and Vr within min range 3km to -99999.9
    DO k =1, count_elev_vel
      DO j=1,ny
        DO i=1,nx
          IF( gridtiltref_recount(i,j,k) < 0.0 .AND. gridtiltref_recount(i,j,k) > -80000.0 ) THEN
              gridtiltref_recount(i,j,k)=0.0
          END if
          IF( gridtiltvel_recount(i,j,k) < -70000.0 ) THEN
              gridtiltvel_recount(i,j,k)= -99999.9
          END IF
          IF( gridrangevel_recount(i,j,k) < 3000.0 .AND. gridtilthighref_recount(i,j,k) < 3000.0) &
              gridtiltvel_recount(i,j,k) = -99999.9
          IF( gridrangevel_recount(i,j,k) < 3000.0 .AND. gridtilthighref_recount(i,j,k) < 3000.0) &
              gridtiltref_recount(i,j,k) = -99999.9
        END DO
      END DO
    END DO

    !nxlg = (nx-3)*nproc_x + 3
    !nylg = (ny-3)*nproc_y + 3

    OPEN(55, FILE=trim(rfnametilt), STATUS='unknown', FORM='unformatted')

    WRITE(55) time1st,iyr,imon,iday,ihour,imin,isec
    WRITE(55) count_elev_vel,nx,ny
    WRITE(55) radar//'      '
    WRITE(55) rdrlat,rdrlon,radarx,radary,rdralt
    WRITE(55) dazim,rngmin,rngmax
    WRITE(55) elevmean_recount_ref(1:count_elev_vel)
    WRITE(55) gridtilthighref_recount(:,:,1:count_elev_vel)
    WRITE(55) gridrangeref_recount(:,:,1:count_elev_vel)
    WRITE(55) gridtiltvel_recount(:,:,1:count_elev_vel)
    WRITE(55) gridtiltref_recount(:,:,1:count_elev_vel)

    CLOSE(55)

    DEALLOCATE( elevmean_recount_vel, elevmean_recount_ref )
    DEALLOCATE( gridrangevel_recount, gridrangeref_recount )
    DEALLOCATE( gridtilthighvel_recount, gridtilthighref_recount )
    DEALLOCATE( gridtiltvel_recount, gridtiltref_recount )

  END IF

  RETURN
END SUBROUTINE wrtgridtilt_old

!##################################################################
!##################################################################
!######                                                      ######
!######               SUBROUTINE REMAP_ARPSPPI               ######
!######                                                      ######
!######                     Developed by                     ######
!######         School of Meteorology / CASA / CAPS          ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE remap_arpsppi(maxgate,maxazim,maxelev,nx,ny,                  &
                           kntazim,kntelev,rdrlat,rdrlon,radarx,           &
                           radary,rdralt,dazim,rngmin,rngmax,              &
                           gatesp,azmvol,elvvol,varvol,time_arpsppi,xs,ys, &
                           gridtilthigh,gridrange,gridslr,gridazm,         &
                           gridtiltval,gridtilttime,elevmean)
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Take data from a radar volume in polar coordinates and map it to
!  ARPS Cartesian terrain-following grid, but don't do vertical interoplation
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Hu, CAPS
!
!      Oct. 13, 2004
!
!  MODIFICATION HISTORY:
!     how to calculate distance between two points on earth
!     how to get radar azimuth angle, range and height
!     how to interpolation on radar tilt
!
!  Nate Snook, SoM/CASA  --  (Jan. 8, 2008)
!     Renamed subroutine 'remap_arpsppi'.
!     Removed unused variables from the subroutine call.
!     Changed implementation for processing of CASA data.
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------
! Input Variables
!-----------------------------------------------------------------------

  INTEGER, INTENT(IN) :: maxgate    !Number of gates per radial.  Equivalent to n_gates in 'casa2arpsppi'.
  INTEGER, INTENT(IN) :: maxazim    !Number of radials per elevation.  Equivalent to n_radials in 'casa2arpsppi'.
  INTEGER, INTENT(IN) :: maxelev    !Number of elevations present.  Equivalent to n_elevs in 'casa2arpsppi'.
  INTEGER, INTENT(IN) :: nx         !Number of ARPS gridpoints in the x-direction.
  INTEGER, INTENT(IN) :: ny         !Number of ARPS gridpoints in the y-direction.

  INTEGER, INTENT(IN) :: kntazim(maxelev)   !REPLACE THIS WITH n_radials FROM CASA2ARPSPPI
  INTEGER, INTENT(IN) :: kntelev            !Number of elevation angles.  In 'casa2arpsppi': n_elevs

  REAL, INTENT(IN)    :: rdrlat     !Radar latitude.  In 'casa2arpsppi': radar_lat
  REAL, INTENT(IN)    :: rdrlon     !Radar longitude.  In 'casa2arpsppi': radar_lon
  REAL, INTENT(IN)    :: radarx     !Radar x-location on ARPS grid.  In 'casa2arpsppi': radar_x
  REAL, INTENT(IN)    :: radary     !Radar y-location on ARPS grid.  In 'casa2arpsppi': radar_y
  REAL, INTENT(IN)    :: rdralt     !Altitude (elevation) of radar.  In 'casa2arpsppi': radar_altitude
  REAL, INTENT(IN)    :: dazim      !Azimuth angle increment.  In 'casa2arpsppi': delta_azim
  REAL, INTENT(IN)    :: rngmin     !Minimum range to use data.  In 'casa2arpsppi': range_min
  REAL, INTENT(IN)    :: rngmax     !Maximum range to use data.  In 'casa2arpsppi': range_max

  INTEGER, INTENT(IN) :: time_arpsppi  !Assimilation time, in absolute format (seconds since 1960).
  REAL, INTENT(IN)    :: gatesp     !Range gate spacing.
  REAL, INTENT(IN)    :: azmvol(maxazim,maxelev)   !REPLACE WITH SOMETHING ELSE OR CALCULATE IN 'casa2arpsppi'
  REAL, INTENT(IN)    :: elvvol(maxazim,maxelev)   !REPLACE WITH elevs FROM 'casa2arpsppi'
  REAL, INTENT(IN)    :: varvol(maxgate,maxazim,maxelev)   !3D volume of a variable.  In 'casa2arpsppi': Z_3d_valid or Vr_3d_valid

  REAL, INTENT(IN)    :: xs(nx)     !x-location of the scalar points.
  REAL, INTENT(IN)    :: ys(ny)     !y-location of the scalar points.

!-----------------------------------------------------------------------
! Output Variables
!-----------------------------------------------------------------------

  REAL, INTENT(OUT)   :: gridtilthigh(nx,ny,maxelev)
  REAL, INTENT(OUT)   :: gridrange(nx,ny,maxelev)
  REAL, INTENT(OUT)   :: gridslr(nx,ny)
  REAL, INTENT(OUT)   :: gridazm(nx,ny)

  REAL, INTENT(OUT)   :: gridtiltval(nx,ny,maxelev)
  INTEGER, INTENT(OUT)   :: gridtilttime(maxelev)

  REAL, INTENT(OUT) :: elevmean(maxelev)

!-----------------------------------------------------------------------
! Local Variables
!-----------------------------------------------------------------------

  INTEGER :: iyr,imon,iday,ihour,imin,isec
  REAL    :: xsll(nx,ny),ysll(nx,ny)

  INTEGER :: i,j,k,knt,kinbox
  INTEGER :: kok,isort,jsort,mid
  INTEGER :: kbgn,kend
  INTEGER :: igate,jazim,kelev,jazmin,jmirror,jend
  INTEGER :: istatal,istatwrt

  INTEGER :: nweight
  REAL :: val, weight

  REAL :: deg2rad,rad2deg
  REAL :: elevmin,elevmax
  REAL :: delx,dely,delz,dazimr,daz,azdiff
  REAL :: ddx,ddxy,ddx2,ddy,ddy2,ddz,dxthr,dxthr0
  REAL :: mapfct,xcomp,ycomp,azmrot,sfcr,zagl
  REAL :: sum,sum2,sdev,thresh,slrange,elijk,azimijk,time
  REAL :: varmax,varmin,varavg,varmean,varmed

  REAL :: kea, skea,elevE
!  REAL :: eradius=6371000     !It never hurts to know how big the planet is.
  REAL :: missvalue = -99999.9
  REAL :: azimj1,azimj2,dazimj,drange,dgate
  INTEGER :: ii,jj
  REAL :: valrd,valru,valld,vallu
  REAL :: rngminMulti(20),rngmink
  REAL :: xylat,xylon
  REAL :: latnot(2)
  REAL :: swx,swy, ctrx, ctry

!-----------------------------------------------------------------------
! Include files
!-----------------------------------------------------------------------

  INCLUDE 'globcst.inc'
  INCLUDE 'grid.inc'

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!  Beginning of executable code...
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  rngminMulti=10
  rngmink=rngmin
  rngmink=10

  deg2rad = atan(1.)/45.
  rad2deg = 1./deg2rad
  dazimr = dazim*deg2rad
  dxthr0=0.6*max((xs(3)-xs(2)),(ys(3)-ys(2)))

  CALL lattomf(1,1,rdrlat,mapfct)

  time=0.
  elevmax=-999.
  elevmin=999.
  DO kelev=1,kntelev
     sum=0.
     DO jazim=1,kntazim(kelev)
        sum=sum+elvvol(jazim,kelev)
        elevmin=min(elevmin,elvvol(jazim,kelev))
        elevmax=max(elevmax,elvvol(jazim,kelev))
     END DO
     elevmean(kelev)=sum/float(kntazim(kelev))
  END DO

  DO kelev=1,kntelev
     jazmin=kntazim(kelev)/2
     gridtilttime(kelev)=time_arpsppi
     if( gridtilttime(kelev) < 0 ) Then
        write(*,*) ' WRONG: radar observation time has problem !!!'
        stop 222
     endif
  ENDDO

!-----------------------------------------------------------------------
! get the latitude and longitude of each grid point
!-----------------------------------------------------------------------

  latnot(1)=trulat1
  latnot(2)=trulat2
  CALL setmapr(mapproj,sclfct,latnot,trulon)

  CALL lltoxy( 1,1, ctrlat,ctrlon, ctrx, ctry )
  swx = ctrx - (FLOAT((nx-3))/2.) * dx
  swy = ctry - (FLOAT((ny-3))/2.) * dy
  CALL setorig( 1, swx, swy)

  gridslr = missvalue
  gridazm = missvalue
  DO j=1,ny-1
     DO i=1,nx-1
        CALL xytoll(1,1,xs(i),ys(j),xylat,xylon)
        CALL disthead(rdrlat,rdrlon,xylat,xylon,azimijk,slrange)
        IF( slrange >= rngmink .and. slrange<=rngmax ) THEN
           gridslr(i,j) = slrange
           gridazm(i,j) = azimijk
        ENDIF  ! slrange >= 10000.0 .and. slrange<=230000.0
     ENDDO
  ENDDO

  write(*,*) ' Radar altitude ==',rdralt
  kea=4.0/3.0*eradius
  gridtilthigh = missvalue
  gridrange   = missvalue
  DO j=1,ny-1
     DO i=1,nx-1

!--------------------------------------------------------------------
! get gridtilthigh and gridrange
!--------------------------------------------------------------------

        slrange=gridslr(i,j)
        IF( slrange >= 0.0 ) THEN
           skea=slrange/kea
           DO k=1,kntelev
              rngmink=rngminMulti(k)
              elevE=elevmean(k)*deg2rad

! according to 2.28a
              gridtilthigh(i,j,k)=kea*(cos(elevE)/cos(elevE+skea)-1.) + rdralt
! according to 2.28c
              gridrange(i,j,k)=(kea+gridtilthigh(i,j,k))*sin(skea)/cos(elevE)

              IF( gridtilthigh(i,j,k) > 20000.0 .or.                                &
                 gridrange(i,j,k) < rngmink .or. gridrange(i,j,k) > rngmax ) THEN
                 gridtilthigh(i,j,k)=  missvalue
                 gridrange(i,j,k) =  missvalue
              ENDIF
           END DO
        ENDIF ! slrange >= 0.0
     ENDDO
  ENDDO

!--------------------------------------------------------------------
!  interpolation
!--------------------------------------------------------------------

  gridtiltval=missvalue
  DO k=1,kntelev
     dgate=gatesp
!    write(*,*) k
     DO j=1,ny-1
        DO i=1,nx-1
           azimijk=gridazm(i,j)
           slrange=gridrange(i,j,k)
           if( azimijk > -88888.0 .and. slrange > -88888.0 ) THEN
              ii=-99999
              jj=-99999
              DO jazim=1,kntazim(k)-1
                 if( jazim == 1 ) then
                    azimj1=azmvol(jazim,k)
                    if(azimj1 > azimijk ) azimijk=azimijk+360.0
                 else
                    azimj1=azimj2
                 endif
                 azimj2=azmvol(jazim+1,k)
                 if ( azimj1 > azimj2 ) azimj2 = azimj2 + 360.0
                 if( azimijk >= azimj1 .and. azimijk < azimj2 ) then
                    ii=jazim
                    dazimj=(azimijk-azimj1)/(azimj2-azimj1)
                    cycle
                 endif
              ENDDO
              if( ii > -8888 ) THEN
                 jj=int(slrange/dgate)+1
                 drange=mod(slrange,dgate)/dgate
                 IF( ii<=maxazim-1 .and. jj<=maxgate-1 ) THEN
                    valld=varvol(jj,ii,k)
                    vallu=varvol(jj+1,ii,k)
                    valrd=varvol(jj,ii+1,k)
                    valru=varvol(jj+1,ii+1,k)
                    if(valld > -488.0 .and. vallu > -488.0 .and.                &
                         valrd > -488.0 .and. valru > -488.0 ) then
                       gridtiltval(i,j,k)=valld*(1.-dazimj)*(1.-drange) +         &
                            vallu*(1.-dazimj)*drange + valrd*dazimj*(1.-drange)+  &
                            valru*dazimj*drange
                    else
                       nweight=0
                       val=0
                       weight=0
                       if(valld > -488.0) then
                          nweight = nweight + 1
                          val= valld * (1.-dazimj)*(1.-drange)
                          weight = weight + (1.-dazimj)*(1.-drange)
                       endif
                       if(vallu > -488.0) then
                          nweight = nweight + 1
                          val= vallu * (1.-dazimj)*drange
                          weight = weight + (1.-dazimj)*drange
                       endif
                       if(valrd > -488.0) then
                          nweight = nweight + 1
                          val= valrd * dazimj*(1.-drange)
                          weight = weight + dazimj*(1.-drange)
                       endif
                       if(valru > -488.0) then
                          nweight = nweight + 1
                          val= valru * dazimj*drange
                          weight = weight + dazimj*drange
                       endif

                       if (nweight > 1 ) then
                          gridtiltval(i,j,k)=val/weight
                       else
                          gridtiltval(i,j,k)=-77777.
                       endif
                    END IF
                 END IF
              END IF
           END IF
        END DO
     END DO
  END DO

  RETURN
END SUBROUTINE remap_arpsppi

!##################################################################
!##################################################################
!######                                                      ######
!######               SUBROUTINE WRITE_ARPSPPI               ######
!######                                                      ######
!######                     Developed by                     ######
!######         School of Meteorology / CASA / CAPS          ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE write_arpsppi(time_arpsppi,iyr,imon,iday,ihr,imin,isec,cur_assim,  &
                           n_elevs,nx,ny,radar_name,radar_lat,radar_lon,        &
                           radar_x,radar_y,radar_altitude,delta_azim,           &
                           range_min,range_max,elevs,                           &
                           gridtilthigh_vel,gridrange_vel,                      &
                           gridtiltval_vel,gridtiltval_ref)
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Take gridtilt data for a CASA radar and output it in a form suitable
!  for use by 'ppiplt' and data assimilation by arpsenkf
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Nate Snook, SoM/CASA
!
!      Jan. 30, 2008
!
!  MODIFICATION HISTORY:
!
!      NONE YET
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------
! Input Variables
!-----------------------------------------------------------------------

  INTEGER, INTENT(IN) :: time_arpsppi    !As time_arpsppi in 'casa2arpsppi'
  INTEGER, INTENT(IN) :: iyr, imon, iday, ihr, imin, isec  !As in 'casa2arpsppi'
  INTEGER, INTENT(IN) :: cur_assim                         !As in 'casa2arpsppi'
  INTEGER, INTENT(IN) :: n_elevs, nx, ny                   !As in 'casa2arpsppi'
  CHARACTER (LEN = 4), INTENT(in) :: radar_name            !As in 'casa2arpsppi'
  REAL, INTENT(IN) :: radar_lat, radar_lon                 !As in 'casa2arpsppi'
  REAL, INTENT(IN) :: radar_x, radar_y                     !As in 'casa2arpsppi'
  REAL, INTENT(IN) :: radar_altitude                       !As in 'casa2arpsppi'
  REAL, INTENT(IN) :: delta_azim                           !As in 'casa2arpsppi'
  REAL, INTENT(IN) :: range_min, range_max                 !As in 'casa2arpsppi'
  REAL, INTENT(IN) :: elevs(n_elevs)                       !As in 'casa2arpsppi', but mapped from INTEGER to REAL
  REAL, INTENT(IN) :: gridtilthigh_vel(nx,ny,n_elevs)      !As in 'casa2arpsppi'
  REAL, INTENT(IN) :: gridrange_vel(nx,ny,n_elevs)         !As in 'casa2arpsppi'
  REAL, INTENT(IN) :: gridtiltval_vel(nx,ny,n_elevs)       !As in 'casa2arpsppi'
  REAL, INTENT(IN) :: gridtiltval_ref(nx,ny,n_elevs)       !As in 'casa2arpsppi'

!-----------------------------------------------------------------------
! Local Variables
!-----------------------------------------------------------------------

  CHARACTER (LEN = 10) :: radar_name_long    !a 10 character long version of 'radar_name'
  CHARACTER (LEN = 64) :: arpsppi_filename   !Name to use for the output file.
  INTEGER :: totalref, totalvel              !Counters for quality check/control of output data.
  INTEGER :: i, j, k                         !Counter variables, nothing more.

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!  Beginning of executable code...
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  WRITE(arpsppi_filename,'(a4,a1,i4.4,2(i2.2),a1,i6.6)')                   &
        radar_name,'.', iyr,imon,iday,'.',cur_assim    !Set a unique rfname for the output file

  radar_name_long = ADJUSTL(radar_name)

  OPEN(55,FILE=arpsppi_filename,STATUS='unknown',FORM='unformatted')
    write(55) time_arpsppi,iyr,imon,iday,ihr,imin,isec
    write(55) n_elevs,nx,ny
    write(55) radar_name_long
    write(55) radar_lat,radar_lon,radar_x,radar_y,radar_altitude
    write(55) delta_azim,range_min,range_max
    write(55) elevs
    write(55) gridtilthigh_vel
    write(55) gridrange_vel
    write(55) gridtiltval_vel
    write(55) gridtiltval_ref
  close(55)

  DO k =1, n_elevs
    totalref = 0
    totalvel = 0
    DO j=1,ny
    DO i=1,nx
      IF(abs(gridtiltval_vel(i,j,k)) < 77777.0 ) totalvel = totalvel + 1
      IF(abs(gridtiltval_ref(i,j,k)) < 77777.0 ) totalref = totalref + 1
    ENDDO
    ENDDO
    write(*,'(a,i2,a,i7,a,i7)') 'tilt=',k,' The number of reflectivity is:',totalref, &
                ' and velocity is:',totalvel
  ENDDO

END SUBROUTINE write_arpsppi
