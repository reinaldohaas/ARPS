!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE PREPRADAR                  ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE prepradar(nx,ny,nz,nz_tab,nvar_anx,nvar_radin,               &
           nvar_rad,mx_rad,nsrc_rad,nz_rdr,mx_colrad,mx_pass,           &
           raduvobs,radrhobs,radistride,radkstride,cloudopt,            &
           iuserad,npass,refrh,rhradobs,                                &
           xs,ys,zs,hterain,latgr,longr,anx,qback,hqback,nlvqback,      &
           nradfil,radfname,srcrad,isrcrad,qsrcrad,qcthrrad,            &
           stnrad,latrad,lonrad,elvrad,                                 &
           latradc,lonradc,xradc,yradc,irad,nlevrad,                    &
           distrad,uazmrad,vazmrad,hgtradc,theradc,trnradc,dsdr,dhdr,   &
           obsrad,oanxrad,odifrad,qobsrad,qualrad,                      &
           ncolrad,ncolrad_mpi,indexrad,oindexrad,                      &
           refmax_mos_3d,ref_mos_3d,                                    &
           rhv_mos_3d,zdr_mos_3d,kdp_mos_3d,dist_mos_3d,                &
           iprocv,jprocv,iradvel,                                       &
           istatus,istat_radar,tem1,tem2,tem3,tem4,tem5,tem6,tem7,tem8)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Read radar data stored as columns on a remapped grid
!  and prepare data for use in the analysis.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster
!  August, 1995
!
!  MODIFICATION HISTORY:
!
!  07/15/97 (K. Brewster)
!  Corrected a bug in interpolation weights.
!
!  Nov 1997 (KB)
!  Added rhmax,refsat and rhrad to argument list.
!
!  Feb, 1998 (KB)
!  Added addition option switches to argument list allowing
!  selection of either radar wind obs and/or radar RH pseudo-obs.
!
!  Mar, 1998 (KB)
!  Added uazmlim to protect against division by zero when
!  azimuth is close to perpendicular to u or v component.
!
!  Dec, 1998 (KB)
!  Made changes for RHstar to qv as moisture variable.
!
!  August, 2001 (KB)
!  Corrected call to dhdr which now returns dhdr not local
!  elevation angle.  Modified calculation of odifrad to speed
!  convergence to vr.
!
!  September, 2002 (KB)
!  Added code to pass along trnradc.
!
!  November, 2007 (KWT)
!  Added code to support MPI of radial velocity.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!    nvar_anx   number of analysis variables
!    nvar_radin number of variables in the obsrad array
!    nvar_rad   number of variables in the odifrad array
!    mx_rad     maximum number of radars
!    nsrc_rad   number of radar sources
!    nz_rdr     maximum number of levels in a radar column
!    mx_colrad  maximum number of radar columns
!
!    refsat     Reflectivity threshold for assigning moisture ob
!    rhrad      Value of relative humidity to assign
!
!    xs         x location of scalar pts (m)
!    ys         y location of scalar pts (m)
!    zs         z location of scalar pts (m)
!    hterain    height of terrain (m MSL)
!    anx        analyzed variables (or first guess)
!    qback      background (first guess) error
!
!    nradfil    number of radar files
!    radfname   file name for radar dataset
!    srcrad     name of radar sources
!    qsrcrad    radar observation error
!
!  OUTPUT:
!
!    isrcrad  index of radar source
!    stnrad   radar site name    character*5
!    latrad   latitude of radar  (degrees N)
!    lonrad   longitude of radar (degrees E)
!    elvrad   elevation of feed horn of radar (m MSL)
!    latradc  latitude of radar column   (degrees N)
!    lonradc  longitude of radar column  (degrees E)
!    xradc    x location of radar column
!    yradc    y location of radar column
!    irad     radar number
!    nlevrad  number of levels of radar data in each column
!    distrad  distance of radar column from source radar
!    uazmrad  u-component of radar beam for each column
!    vazmrad  v-component of radar beam for each column
!    hgtradc  height (m MSL) of radar observations
!    obsrad   radar observations
!    oanxrad  analysis (first guess) value at radar data location
!    odifrad  difference between radar observation and analysis
!    qobsrad  normalized observation error
!    qualrad  radar data quality indicator
!    ncolrad  number of radar columns read-in, local for MPI
!    ncolrad_mpi  number of radar columns read-in, MPI only
!    indexrad index of data point processor owners (MPI only)
!    oindexrad 2nd index of data point processor owners (MPI only, iradvel>0)
!    iprocv   max proc radius of influence x-direction for radar vel (MPI only)
!    jprocv   max proc radius of influence y-direction for radar vel (MPI only)
!    iradvel  use faster radial velocity processing (less memory also)
!    istatus  status indicator
!    istat_radar ref_mos_3d status indicator
!
!-----------------------------------------------------------------------
!
!  Variable Declarations:
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
  INTEGER :: nx,ny,nz,nz_tab,nvar_anx
  INTEGER :: nvar_radin,nvar_rad,mx_rad,nsrc_rad,nz_rdr,mx_colrad
  INTEGER :: mx_pass
  INTEGER :: raduvobs
  INTEGER :: radrhobs
  INTEGER :: radistride
  INTEGER :: radkstride
  INTEGER :: cloudopt
  INTEGER :: iuserad(0:nsrc_rad,mx_pass)
  INTEGER :: npass
  INTEGER :: np
  REAL :: refrh
  REAL :: rhradobs
!
!-----------------------------------------------------------------------
!
!  Grid variables
!
!-----------------------------------------------------------------------
!
  REAL :: xs(nx)
  REAL :: ys(ny)
  REAL :: zs(nx,ny,nz)
  REAL :: hterain(nx,ny)
  REAL :: latgr(nx,ny)
  REAL :: longr(nx,ny)
  REAL :: anx(nx,ny,nz,nvar_anx)
  REAL :: qback(nvar_anx,nz_tab)
  REAL :: hqback(nz_tab)
  INTEGER :: nlvqback

  INTEGER :: nradfil
  CHARACTER (LEN=256) :: radfname(mx_rad)
!
!-----------------------------------------------------------------------
!
!  Radar site variables
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=8) :: srcrad(nsrc_rad)
  INTEGER :: isrcrad(0:mx_rad)
  REAL :: qsrcrad(nvar_radin,nsrc_rad)
  REAL :: qcthrrad(nvar_radin,nsrc_rad)
  REAL :: latrad(mx_rad),lonrad(mx_rad)
  REAL :: elvrad(mx_rad)
  CHARACTER (LEN=5) :: stnrad(mx_rad)
!
!-----------------------------------------------------------------------
!
!  Radar observation variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: irad(mx_colrad)
  INTEGER :: nlevrad(mx_colrad)
  REAL :: distrad(mx_colrad),uazmrad(mx_colrad),vazmrad(mx_colrad)
  REAL :: latradc(mx_colrad),lonradc(mx_colrad)
  REAL :: xradc(mx_colrad),yradc(mx_colrad)
  REAL :: hgtradc(nz_rdr,mx_colrad)
  REAL :: theradc(nz_rdr,mx_colrad)
  REAL :: trnradc(mx_colrad)
  REAL ::    dsdr(nz_rdr,mx_colrad)
  REAL ::    dhdr(nz_rdr,mx_colrad)

  REAL :: obsrad(nvar_radin,nz_rdr,mx_colrad)
  REAL :: oanxrad(nvar_rad,nz_rdr,mx_colrad)
  REAL :: odifrad(nvar_rad,nz_rdr,mx_colrad)
  REAL :: qobsrad(nvar_rad,nz_rdr,mx_colrad)
  INTEGER :: qualrad(nvar_rad,nz_rdr,mx_colrad)
  REAL :: refmax_mos_3d(nx,ny,nz)
  REAL :: ref_mos_3d(nx,ny,nz)
  REAL :: rhv_mos_3d(nx,ny,nz)
  REAL :: zdr_mos_3d(nx,ny,nz)
  REAL :: kdp_mos_3d(nx,ny,nz)
  REAL :: dist_mos_3d(nx,ny,nz)
  INTEGER :: ncolrad, ncolrad_mpi
  INTEGER :: istatus
  INTEGER :: istat_radar
  INTEGER :: iradvel
!
!-----------------------------------------------------------------------
!
!  MPI variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: indexrad(mx_colrad)
  INTEGER :: oindexrad(mx_colrad)
  INTEGER :: iprocv,jprocv
!
!-----------------------------------------------------------------------
!
  REAL :: tem1(nx,ny,nz)
  REAL :: tem2(nx,ny,nz)
  REAL :: tem3(nx,ny,nz)
  REAL :: tem4(nx,ny,nz)
  REAL :: tem5(nx,ny,nz)
  REAL :: tem6(nx,ny,nz)
  REAL :: tem7(nx,ny,nz)
  REAL :: tem8(nx,ny,nz)
!
!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: iusechk
  INTEGER :: ierror
  INTEGER :: k

  INCLUDE 'mp.inc'
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!-----------------------------------------------------------------------
!
!  Read remapped radar data
!
!-----------------------------------------------------------------------
!
  isrcrad(0)=0
  ncolrad = 0
  ncolrad_mpi = 0
  IF(raduvobs > 0 .OR. radrhobs > 0 .OR. cloudopt > 0) THEN
    iusechk=1
    CALL rdradcol(nx,ny,nz,nsrc_rad,nvar_radin,                         &
             mx_rad,nz_rdr,mx_colrad,mx_pass,                           &
             raduvobs,radrhobs,radistride,radkstride,cloudopt,          &
             xs,ys,zs,latgr,longr,                                      &
             iuserad,iusechk,npass,nradfil,radfname,                    &
             srcrad,isrcrad(1),stnrad,latrad,lonrad,elvrad,             &
             latradc,lonradc,irad,nlevrad,hgtradc,obsrad,               &
             iprocv,jprocv,iradvel,indexrad,oindexrad,                  &
             refmax_mos_3d,ref_mos_3d,                                  &
             rhv_mos_3d,zdr_mos_3d,kdp_mos_3d,dist_mos_3d,              &
             ncolrad,ncolrad_mpi,istatus,istat_radar,                   &
             tem1,tem2,tem3,tem4,tem5,tem6,tem7,tem8)

    print *, ' Back from rdradcol, istat_radar =',istat_radar

    CALL mpbarrier
    IF (mp_opt == 0) THEN
      WRITE(6,'(a,i7/)') ' Prepradar: ncolrad = ',ncolrad
    ELSE
      IF (myproc == 0) WRITE(6,'(a,i7/)') ' Prepradar: ncolrad = ',ncolrad_mpi
    END IF

    IF(ncolrad > 0 .OR. ncolrad_mpi > 0) THEN
!
!-----------------------------------------------------------------------
!
!  Remove precipitation terminal velocity from the radial velocity.
!
!-----------------------------------------------------------------------
!
      CALL rmvterm(nvar_radin,mx_rad,nz_rdr,mx_colrad,                  &
               latrad,lonrad,elvrad,                                    &
               latradc,lonradc,irad,nlevrad,hgtradc,obsrad,             &
               ncolrad,istatus)
!
!-----------------------------------------------------------------------
!
!  Unfold and QC radar data
!
!-----------------------------------------------------------------------
!
      CALL prcssrad(nx,ny,nz,nz_tab,nvar_anx,nvar_radin,                &
               nvar_rad,mx_rad,nsrc_rad,nz_rdr,mx_colrad,               &
               raduvobs,radrhobs,refrh,rhradobs,                        &
               xs,ys,zs,hterain,anx,qback,hqback,nlvqback,              &
               srcrad,qsrcrad,qcthrrad,                                 &
               isrcrad,stnrad,latrad,lonrad,elvrad,                     &
               latradc,lonradc,xradc,yradc,irad,nlevrad,                &
               hgtradc,obsrad,oanxrad,odifrad,                          &
               distrad,uazmrad,vazmrad,                                 &
               theradc,trnradc,qobsrad,qualrad,dsdr,dhdr,               &
               ncolrad,indexrad,istatus,tem1)
    END IF
  ELSE
    ncolrad=0
  END IF

  RETURN
END SUBROUTINE prepradar
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE PRCSSRAD                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE prcssrad(nx,ny,nz,nz_tab,nvar_anx,nvar_radin,                &
           nvar_rad,mx_rad,nsrc_rad,nz_rdr,mx_colrad,                   &
           raduvobs,radrhobs,refrh,rhradobs,                            &
           xs,ys,zs,hterain,anx,qback,hqback,nlvqback,                  &
           srcrad,qsrcrad,qcthrrad,                                     &
           isrcrad,stnrad,latrad,lonrad,elvrad,                         &
           latradc,lonradc,xradc,yradc,irad,nlevrad,                    &
           hgtradc,obsrad,oanxrad,odifrad,                              &
           distrad,uazmrad,vazmrad,                                     &
           theradc,trnradc,qobsrad,qualrad,dsdr,dhdr,                   &
           ncolrad,indexrad,istatus,tem1)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Unfolds and QC's radar data stored as columns on a remapped grid.
!  For MPI, data point owners also get assigned.
!
!  AUTHOR: Keith Brewster
!  August, 1995
!
!  MODIFICATION HISTORY:
!
!  September, 2002 (KB)
!  Added code to calculate and store trnradc.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!    nvar_anx   number of analysis variables
!    nvar_radin number of variables in the obsrad array
!    nvar_rad   number of variables in the odifrad array
!    mx_rad     maximum number of radars
!    nsrc_rad   number of radar sources
!    nz_rdr     maximum number of levels in a radar column
!    mx_colrad  maximum number of radar columns
!
!    xs         x location of scalar pts (m)
!    ys         y location of scalar pts (m)
!    zs         z location of scalar pts (m)
!    hterain    height of terrain (m MSL)
!    anx        analyzed variables (or first guess)
!    qback      background (first guess) error
!    srcrad   name of radar sources
!    qsrcrad  radar observation error
!    isrcrad  index of radar source
!    stnrad   radar site name    character*5
!    latrad   latitude of radar  (degrees N)
!    lonrad   longitude of radar (degrees E)
!    elvrad   elevation of feed horn of radar (m MSL)
!    latradc  latitude of radar column   (degrees N)
!    lonradc  longitude of radar column  (degrees E)
!    xradc    x location of radar column (m)
!    yradc    y location of radar column (m)
!    irad     radar number
!    nlevrad  number of levels of radar data in each column
!    distrad  distance of radar column from source radar
!    uazmrad  u-component of radar beam for each column
!    vazmrad  v-component of radar beam for each column
!    hgtradc  height (m MSL) of radar observations
!    obsrad   radar observations
!
!  OUTPUT:
!
!    oanxrad  analysis (first guess) value at radar data location
!    odifrad  difference between radar observation and analysis
!    theradc  theta (potential temperature K) at each radar ob
!    trnradc  terrain height at each radar column
!    qobsrad  normalized observation error
!    qualrad  radar data quality indicator
!    ncolrad  number of radar columns read-in
!    indexrad index of data point processor owners (MPI only)
!    istatus  status indicator
!
!-----------------------------------------------------------------------
!
!  Variable Declarations:
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
  INTEGER :: nx,ny,nz,nz_tab,nvar_anx
  INTEGER :: nvar_radin,nvar_rad,mx_rad,nsrc_rad,nz_rdr,mx_colrad
  INTEGER :: raduvobs
  INTEGER :: radrhobs
  REAL :: refrh
  REAL :: rhradobs
!
!-----------------------------------------------------------------------
!
!  Grid variables
!
!-----------------------------------------------------------------------
!
  REAL :: xs(nx)
  REAL :: ys(ny)
  REAL :: zs(nx,ny,nz)
  REAL :: hterain(nx,ny)
  REAL :: anx(nx,ny,nz,nvar_anx)
  REAL :: qback(nvar_anx,nz_tab)
  REAL :: hqback(nz_tab)
  INTEGER :: nlvqback
!
!-----------------------------------------------------------------------
!
!  Radar site variables
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=8) :: srcrad(nsrc_rad)
  INTEGER :: isrcrad(0:mx_rad)
  REAL :: qsrcrad(nvar_radin,nsrc_rad)
  REAL :: qcthrrad(nvar_radin,nsrc_rad)
  REAL :: latrad(mx_rad),lonrad(mx_rad)
  REAL :: elvrad(mx_rad)
  CHARACTER (LEN=5) :: stnrad(mx_rad)
!
!-----------------------------------------------------------------------
!
!  Radar observation variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: irad(mx_colrad)
  INTEGER :: nlevrad(mx_colrad)
  REAL :: latradc(mx_colrad),lonradc(mx_colrad)
  REAL :: xradc(mx_colrad),yradc(mx_colrad)
  REAL :: hgtradc(nz_rdr,mx_colrad)
  REAL :: theradc(nz_rdr,mx_colrad)
  REAL :: trnradc(mx_colrad)
  REAL ::    dsdr(nz_rdr,mx_colrad)
  REAL ::    dhdr(nz_rdr,mx_colrad)

  REAL :: obsrad(nvar_radin,nz_rdr,mx_colrad)
  REAL :: oanxrad(nvar_rad,nz_rdr,mx_colrad)
  REAL :: odifrad(nvar_rad,nz_rdr,mx_colrad)
  REAL :: qobsrad(nvar_rad,nz_rdr,mx_colrad)
  REAL :: distrad(mx_colrad),uazmrad(mx_colrad),vazmrad(mx_colrad)
  INTEGER :: qualrad(nvar_rad,nz_rdr,mx_colrad)
  INTEGER :: ncolrad
  INTEGER :: istatus
!
!-----------------------------------------------------------------------
!
!  MPI variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: indexrad(mx_colrad)
!
!-----------------------------------------------------------------------
!
!  Temporary work array
!
!-----------------------------------------------------------------------
!
  REAL :: tem1(nx*ny*nz)
!
!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: icol,jrad,kr
  INTEGER :: i,imid,iloc
  INTEGER :: j,jmid,jloc
  INTEGER :: k,kmid,kloc,ktab
  INTEGER :: kfold
  REAL :: xmid,ymid,zmid
  REAL :: wx,wy,w11,w12,w21,w22,whi,wlo,wqhi,wqlo
  REAL :: ddrot,ugrid,vgrid,prgrid,tgrid,qvgrid,qvsgrid,rhgrid,wgrid
  REAL :: vr_miss,qv_miss
  REAL :: azim,dist,vr,vrdiff,avrdif,dz,eleva,range,dsdrinv
!
!-----------------------------------------------------------------------
!
!  Include file
!
!-----------------------------------------------------------------------
!
  INCLUDE 'phycst.inc'
  INCLUDE 'mp.inc'
!
!-----------------------------------------------------------------------
!
!  QC and conversion parameters
!
!  vrqcthr   maximum radial velocity diff from background to use
!  vrcnvthr  lower limit of the dot product between the radial
!            direction and the analysis wind component to use data
!
!-----------------------------------------------------------------------
!
  REAL :: vrmax,vrqcthr,uazmlim
  PARAMETER (vrmax=70.,                                                 &
             vrqcthr=15.,   & ! maximum vrdiff from background to use
         uazmlim=0.087) ! sin(5 degrees)
  REAL :: refmax, refmin1, refmin2, refmin
  PARAMETER (refmax = 80.)  ! maximum non-missing reflectivity
  PARAMETER (refmin1= 25.)  ! maximin non-missing reflectivity
  PARAMETER (refmin2= 15.)  ! maximin non-missing reflectivity
!
!-----------------------------------------------------------------------
!
!  Function f_qvsat and inline directive for Cray PVP
!
!-----------------------------------------------------------------------
!
  REAL :: f_qvsat

!fpp$ expand (f_qvsat)
!!dir$ inline always f_qvsat
!*$*  inline routine (f_qvsat)

!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!-----------------------------------------------------------------------
!
!  Some initializations
!
!-----------------------------------------------------------------------
!
  vr_miss=-999.
  qv_miss=-999.
!
!-----------------------------------------------------------------------
!
!  Get x and y locations of each radar ob
!
!-----------------------------------------------------------------------
!
  CALL lltoxy(ncolrad,1,latradc,lonradc,xradc,yradc)
!
  imid=nx/2
  xmid=xs(imid)
  jmid=ny/2
  ymid=ys(jmid)
  kmid=nz/2

  DO icol=1,ncolrad

    IF(mp_opt > 0 .AND. indexrad(icol) .NE. myproc) CYCLE

    IF( irad(icol) > 0 .AND. nlevrad(icol) > 0) THEN
      IF( xradc(icol) >= xs(1) .AND. xradc(icol) <= xs(nx-1) .AND.      &
          yradc(icol) >= ys(1) .AND. yradc(icol) <= ys(ny-1)) THEN
!
!-----------------------------------------------------------------------
!
!  Find column in grid
!
!-----------------------------------------------------------------------
!
        IF(xradc(icol) < xmid) THEN
          DO i=imid,2,-1
            IF(xs(i) < xradc(icol)) EXIT
          END DO
          iloc=i
        ELSE
          DO i=imid,nx-1
            IF(xs(i) > xradc(icol)) EXIT
          END DO
          iloc=i-1
        END IF
        wx =  (xradc(icol)-xs(iloc))/                                   &
              (xs(iloc+1)-xs(iloc))

        IF(yradc(icol) < ymid) THEN
          DO j=jmid,2,-1
            IF(ys(j) < yradc(icol)) EXIT
          END DO
          jloc=j
        ELSE
          DO j=jmid,ny-1
            IF(ys(j) > yradc(icol)) EXIT
          END DO
          jloc=j-1
        END IF
        wy = (yradc(icol)-ys(jloc))/                                    &
              (ys(jloc+1)-ys(jloc))
!
!-----------------------------------------------------------------------
!
!  Determine bilinear interpolation weights
!
!-----------------------------------------------------------------------
!
        w22 = wx*wy
        w12 = (1.0 - wx) * wy
        w21 = wx * (1.0 - wy)
        w11 = (1.0 - wx) * (1.0 - wy)
!
!-----------------------------------------------------------------------
!
!  Find heading and distance from radar along ground
!  Store correlation of azimuth with u and v drections
!  in uazmrad and vazmrad, respectively.
!
!-----------------------------------------------------------------------
!
        jrad=irad(icol)
        CALL disthead(latradc(icol),lonradc(icol),                      &
                      latrad(jrad),lonrad(jrad),                        &
                      azim,dist)
        distrad(icol)=dist

        CALL ddrotuv(1,lonradc(icol),azim,1.,ddrot,                     &
                       uazmrad(icol),vazmrad(icol))
!
!    in between added by GAO !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     CALL lltoxy(1,1,latradc(icol),lonradc(icol),xx1,yy1)
!     CALL lltoxy(1,1,latrad(jrad),lonrad(jrad),xx2,yy2)
!      uazmrad(icol)=xx1-xx2
!      vazmrad(icol)=yy1-yy2
!
!    in between added by GAO !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!-----------------------------------------------------------------------
!
!  Interpolate grid heights in horizontal
!
!-----------------------------------------------------------------------
!
        DO k=1,nz-1
          tem1(k)=w11*zs(  iloc,  jloc,k) +                             &
                  w21*zs(iloc+1,  jloc,k) +                             &
                  w12*zs(  iloc,jloc+1,k) +                             &
                  w22*zs(iloc+1,jloc+1,k)
        END DO
        trnradc(icol)=0.5*(tem1(1)+tem1(2))
!
!-----------------------------------------------------------------------
!
!  Loop in height
!
!-----------------------------------------------------------------------
!
        DO kr=1,nlevrad(icol)

          IF ( hgtradc(kr,icol) >= 1500.0) THEN
            refmin =refmin2
          ELSE
            refmin =refmin1
          END IF

          IF( ( ABS(obsrad(2,kr,icol)) < vrmax  .AND.                   &
                ABS(obsrad(1,kr,icol)) < refmax .AND.                   &
                ABS(obsrad(1,kr,icol)) > refmin )      .AND.            &
              ( hgtradc(kr,icol) >= tem1(1)     .AND.                   &
                hgtradc(kr,icol) <= tem1(nz-1) )                        &
            ) THEN

            dz=hgtradc(kr,icol)-elvrad(irad(icol))
            CALL beamelv(dz,dist,eleva,range)
            
!g          CALL dhdrange(eleva,range,dhdr)
!g          dsdr=SQRT(MAX(0.,(1.-dhdr*dhdr)))

            CALL dhdrange( eleva,range,dhdr(kr,icol) )
            dsdr(kr,icol)=SQRT(AMAX1(0.,(1.-dhdr(kr,icol)*dhdr(kr,icol))))

!           IF(dsdr /= 0.) THEN
!             dsdrinv=1./dsdr
!           ELSE
!             dsdrinv=0.
!           END IF

            IF(dsdr(kr,icol) /= 0.) THEN
              dsdrinv=1./dsdr(kr,icol)
            ELSE
              dsdrinv=0.
            END IF
!
!-----------------------------------------------------------------------
!
!  Find z location in grid
!
!-----------------------------------------------------------------------
!
            zmid=tem1(kmid)
            IF(hgtradc(kr,icol) < zmid) THEN
              DO k=kmid,2,-1
                IF(tem1(k) < hgtradc(kr,icol)) EXIT
              END DO
              kloc=k
            ELSE
              DO k=kmid,nz-1
                IF(tem1(k) > hgtradc(kr,icol)) EXIT
              END DO
              kloc=k-1
            END IF
!
!-----------------------------------------------------------------------
!
!  Set z weights
!
!-----------------------------------------------------------------------
!
            whi = (hgtradc(kr,icol)-tem1(kloc))/                          &
                  (tem1(kloc+1)-tem1(kloc))
            wlo = 1.0 - whi
!
!-----------------------------------------------------------------------
!
!  Set z weights for error data
!
!-----------------------------------------------------------------------
!
            DO k=2,nlvqback-1
              IF(hqback(k) > hgtradc(kr,icol)) EXIT
            END DO
            ktab=k-1
            wqhi = (hgtradc(kr,icol)-hqback(ktab))/                       &
                   (hqback(ktab+1)-hqback(ktab))
            wqlo = 1.0 - wqhi
!
!-----------------------------------------------------------------------
!
!  Interpolate background field to radar location
!
!  Note this is dependent on indices set so that
!  u=1,v=2,theta=4
!  This can be made more general in a future version.
!
!-----------------------------------------------------------------------
!
            ugrid=                                                        &
                wlo* (w11*anx(  iloc,  jloc,  kloc,1) +                   &
                    w21*anx(iloc+1,  jloc,  kloc,1) +                     &
                    w12*anx(  iloc,jloc+1,  kloc,1) +                     &
                    w22*anx(iloc+1,jloc+1,  kloc,1))                      &
                + whi* (w11*anx(  iloc,  jloc,kloc+1,1) +                 &
                    w21*anx(iloc+1,  jloc,kloc+1,1) +                     &
                    w12*anx(  iloc,jloc+1,kloc+1,1) +                     &
                    w22*anx(iloc+1,jloc+1,kloc+1,1))
            vgrid=                                                        &
                wlo* (w11*anx(  iloc,  jloc,  kloc,2) +                   &
                    w21*anx(iloc+1,  jloc,  kloc,2) +                     &
                    w12*anx(  iloc,jloc+1,  kloc,2) +                     &
                    w22*anx(iloc+1,jloc+1,  kloc,2))                      &
                + whi* (w11*anx(  iloc,  jloc,kloc+1,2) +                 &
                  w21*anx(iloc+1,  jloc,kloc+1,2) +                       &
                  w12*anx(  iloc,jloc+1,kloc+1,2) +                       &
                  w22*anx(iloc+1,jloc+1,kloc+1,2))
            prgrid=                                                       &
                wlo* (w11*anx(  iloc,  jloc,  kloc,3) +                   &
                    w21*anx(iloc+1,  jloc,  kloc,3) +                     &
                    w12*anx(  iloc,jloc+1,  kloc,3) +                     &
                    w22*anx(iloc+1,jloc+1,  kloc,3))                      &
                + whi* (w11*anx(  iloc,  jloc,kloc+1,3) +                 &
                    w21*anx(iloc+1,  jloc,kloc+1,3) +                     &
                    w12*anx(  iloc,jloc+1,kloc+1,3) +                     &
                    w22*anx(iloc+1,jloc+1,kloc+1,3))
            theradc(kr,icol)=                                             &
                wlo* (w11*anx(  iloc,  jloc,  kloc,4) +                   &
                    w21*anx(iloc+1,  jloc,  kloc,4) +                     &
                    w12*anx(  iloc,jloc+1,  kloc,4) +                     &
                    w22*anx(iloc+1,jloc+1,  kloc,4))                      &
                + whi* (w11*anx(  iloc,  jloc,kloc+1,4) +                 &
                    w21*anx(iloc+1,  jloc,kloc+1,4) +                     &
                    w12*anx(  iloc,jloc+1,kloc+1,4) +                     &
                    w22*anx(iloc+1,jloc+1,kloc+1,4))
            qvgrid=                                                       &
                wlo* (w11*anx(  iloc,  jloc,  kloc,5) +                   &
                    w21*anx(iloc+1,  jloc,  kloc,5) +                     &
                    w12*anx(  iloc,jloc+1,  kloc,5) +                     &
                    w22*anx(iloc+1,jloc+1,  kloc,5))                      &
                + whi* (w11*anx(  iloc,  jloc,kloc+1,5) +                 &
                    w21*anx(iloc+1,  jloc,kloc+1,5) +                     &
                    w12*anx(  iloc,jloc+1,kloc+1,5) +                     &
                    w22*anx(iloc+1,jloc+1,kloc+1,5))
#ifdef arps3dvar
            wgrid=                                                        &
                wlo* (w11*anx(  iloc,  jloc,  kloc,6) +                   &
                    w21*anx(iloc+1,  jloc,  kloc,6) +                     &
                    w12*anx(  iloc,jloc+1,  kloc,6) +                     &
                    w22*anx(iloc+1,jloc+1,  kloc,6))                      &
                + whi* (w11*anx(  iloc,  jloc,kloc+1,6) +                 &
                  w21*anx(iloc+1,  jloc,kloc+1,6) +                       &
                  w12*anx(  iloc,jloc+1,kloc+1,6) +                       &
                  w22*anx(iloc+1,jloc+1,kloc+1,6))
#endif

            IF(ABS(obsrad(2,kr,icol)) < vrmax .AND. raduvobs > 0) THEN
!
!-----------------------------------------------------------------------
!
!  Find earth-relative u,v
!
!-----------------------------------------------------------------------
!
#ifdef adas
              oanxrad(1,kr,icol)=ugrid
              oanxrad(2,kr,icol)=vgrid
#endif
!
!-----------------------------------------------------------------------
!
!  Find radial velocity
!  This assumes motions are quasi-horizontal.
!
!-----------------------------------------------------------------------
!
#ifdef adas
              vr=(uazmrad(icol)*ugrid + vazmrad(icol)*vgrid) * dsdr(kr,icol)
#else         
              vr=(uazmrad(icol)*ugrid + vazmrad(icol)*vgrid) * dsdr(kr,icol) &
                                                      + wgrid*dhdr(kr,icol)
              oanxrad(2,kr,icol)=vr
#endif
!
!-----------------------------------------------------------------------
!
!  Check for folding
!
!-----------------------------------------------------------------------
!
              vrdiff=obsrad(2,kr,icol)-vr
              IF(obsrad(3,kr,icol) > 0.0) THEN ! Check that Nyquist velocity is sane
                kfold=nint(vrdiff/(2.*obsrad(3,kr,icol)))
                obsrad(2,kr,icol)=obsrad(2,kr,icol)-                        &
                                  kfold*2.*obsrad(3,kr,icol)
                vrdiff=obsrad(2,kr,icol)-vr
              END IF
              avrdif=ABS(vrdiff)
!
!-----------------------------------------------------------------------
!
!  QC check and calculation of increments.
!  QC should be more sophisticated at some point.
!
!-----------------------------------------------------------------------
!
              IF( avrdif < qcthrrad(2,isrcrad(irad(icol))) ) THEN
!
!-----------------------------------------------------------------------
!
!  Find increment in analysis u and v needed to make vrdiff zero
!
!-----------------------------------------------------------------------
!
#ifdef adas
                odifrad(1,kr,icol)=dsdrinv*(uazmrad(icol)*vrdiff)
                odifrad(2,kr,icol)=dsdrinv*(vazmrad(icol)*vrdiff)
#else
                odifrad(2,kr,icol)=obsrad(2,kr,icol)-vr
#endif
!
!-----------------------------------------------------------------------
!
!  Assign obs error to u according to observation angle
!
!-----------------------------------------------------------------------
!
#ifdef adas
                IF(ABS(uazmrad(icol)) > uazmlim) THEN
                  qobsrad(1,kr,icol)=qsrcrad(2,isrcrad(irad(icol)))/      &
                       ((uazmrad(icol)*uazmrad(icol))*                    &
                        (wqlo*qback(1,ktab)+wqhi*qback(1,ktab+1)))
                  qualrad(1,kr,icol)=10
                ELSE
                  qobsrad(1,kr,icol)=vr_miss
                  qualrad(1,kr,icol)=-90
                END IF
#endif
!
!-----------------------------------------------------------------------
!
!  Assign obs error to v according to observation angle
!
!-----------------------------------------------------------------------
!
#ifdef adas
                IF(ABS(vazmrad(icol)) > uazmlim) THEN
                  qobsrad(2,kr,icol)=qsrcrad(2,isrcrad(irad(icol)))/      &
                      ((vazmrad(icol)*vazmrad(icol))*                     &
                      (wqlo*qback(2,ktab)+wqhi*qback(2,ktab+1)))
                  qualrad(2,kr,icol)=10
                ELSE
                  qobsrad(2,kr,icol)=vr_miss
                  qualrad(2,kr,icol)=-90
                END IF
#else
                qobsrad(2,kr,icol)=qsrcrad(2,isrcrad(irad(icol)))
                qualrad(2,kr,icol)=10
#endif

              ELSE

                odifrad(1,kr,icol)=vr_miss
                odifrad(2,kr,icol)=vr_miss
                qobsrad(1,kr,icol)=vr_miss
                qobsrad(2,kr,icol)=vr_miss
                qualrad(1,kr,icol)=-99
                qualrad(2,kr,icol)=-99

              END IF

            ELSE

              odifrad(1,kr,icol)=vr_miss
              odifrad(2,kr,icol)=vr_miss
              qobsrad(1,kr,icol)=vr_miss
              qobsrad(2,kr,icol)=vr_miss
              qualrad(1,kr,icol)=-99
              qualrad(2,kr,icol)=-99

            END IF
!
!-----------------------------------------------------------------------
!
!    Where reflectivity is greater than that
!    assumed for "saturation" assign an rhs corresponding to
!    rh = rhind, but not if the background rh is greater
!    than rhradobs (background more saturated).
!
!-----------------------------------------------------------------------
!
            IF(radrhobs > 0 .AND. obsrad(1,kr,icol) > refrh .AND.         &
                  obsrad(1,kr,icol) < refmax ) THEN
              tgrid=theradc(kr,icol)*(prgrid/p0)**rddcp
              qvsgrid=f_qvsat(prgrid,tgrid)
              rhgrid=qvgrid/qvsgrid
              IF( rhgrid < rhradobs ) THEN
#ifdef adas
                obsrad(1,kr,icol) = rhradobs*qvsgrid
                odifrad(3,kr,icol) = obsrad(1,kr,icol)-qvgrid
                oanxrad(3,kr,icol) = qvgrid
                qualrad(3,kr,icol) = 10
                qobsrad(3,kr,icol) = qsrcrad(1,isrcrad(irad(icol)))*      &
                                     qvsgrid
#else
                odifrad(1,kr,icol) = obsrad(1,kr,icol)-qvgrid
                oanxrad(1,kr,icol) = qvgrid
                qualrad(1,kr,icol) = 10
                qobsrad(1,kr,icol) = qsrcrad(1,isrcrad(irad(icol)))*      &
                                     qvsgrid
#endif
              ELSE
#ifdef adas
                obsrad(1,kr,icol) = qv_miss
                odifrad(3,kr,icol) = qv_miss
                oanxrad(3,kr,icol) = qv_miss
                qobsrad(3,kr,icol) = qv_miss
                qualrad(3,kr,icol) =-99
#else
                odifrad(1,kr,icol) = qv_miss
                oanxrad(1,kr,icol) = qv_miss
                qobsrad(1,kr,icol) = qv_miss
                qualrad(1,kr,icol) =-99
#endif
              END IF
            ELSE
#ifdef adas
              obsrad(1,kr,icol) = qv_miss
              odifrad(3,kr,icol) = qv_miss
              oanxrad(3,kr,icol) = qv_miss
              qobsrad(3,kr,icol) = qv_miss
              qualrad(3,kr,icol) =-99
#else
              odifrad(1,kr,icol) = qv_miss
              oanxrad(1,kr,icol) = qv_miss
              qobsrad(1,kr,icol) = qv_miss
              qualrad(1,kr,icol) =-99
#endif
            END IF

          ELSE

            odifrad(1,kr,icol)=vr_miss
            odifrad(2,kr,icol)=vr_miss
            qobsrad(1,kr,icol)=vr_miss
            qobsrad(2,kr,icol)=vr_miss
            qualrad(1,kr,icol)=-99
            qualrad(2,kr,icol)=-99
#ifdef adas
            obsrad(1,kr,icol) = qv_miss
            odifrad(3,kr,icol) = qv_miss
            oanxrad(3,kr,icol) = qv_miss
            qobsrad(3,kr,icol) = qv_miss
            qualrad(3,kr,icol) =-99
#else
            odifrad(1,kr,icol) = qv_miss
            oanxrad(1,kr,icol) = qv_miss
            qobsrad(1,kr,icol) = qv_miss
            qualrad(1,kr,icol) =-99
#endif

          END IF
  
        END DO
      ELSE
        irad(icol)=0
      END IF

    END IF
  END DO

  RETURN
END SUBROUTINE prcssrad
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE RADMCRO                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE radmcro(nx,ny,nz,ref_mos_3d,                                 &
           radqvopt,radqcopt,radqropt,radptopt,                         &
           refsat,rhradobs,                                             &
           refcld,cldrad,ceilopt,ceilmin,dzfill,                        &
           refrain,radsetrat,radreflim,radptgain,                       &
           xs,ys,zs,zp,j3,pr,pt,qv,                                     &
           ptbar,qvbar,rhostr,                                          &
           qscalar,ceillim )
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Assign values to microphysical variables based on observed
!  radar data alone.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster
!  November, 1995
!
!  MODIFICATION HISTORY:
!  December, 1997
!  February, 1998
!  September, 1998 Jingsing Zong
!  Corrected logic after DO 150 to continue to next column.
!
!  April, 2010
!  Entirely rewritten to use the ref_mos_3d array instead of the individual
!  radar columns.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!    xs         x location of scalar pts (m)
!    ys         y location of scalar pts (m)
!    zs         z location of scalar pts (m)
!
!  OUTPUT:
!
!    qc         cloud water mixing ratio
!
!  WORK ARRAYS
!
!    ceillim    Ceiling limit (m MSL)
!    tem1       Temporary array
!
!-----------------------------------------------------------------------
!
!  Variable Declarations:
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INCLUDE 'globcst.inc'

  INTEGER :: nx,ny,nz
!
!-----------------------------------------------------------------------
!
!  Control parameters
!
!-----------------------------------------------------------------------
!
  INTEGER :: radqvopt,radqcopt,radqropt,radptopt
  REAL :: refsat,rhradobs,refcld,cldrad
  INTEGER :: ceilopt
  REAL :: ceilmin,dzfill
  REAL :: refrain
  REAL :: radsetrat,radreflim,radptgain
!
!-----------------------------------------------------------------------
!
!  Radar observation variables
!
!-----------------------------------------------------------------------
!
  REAL :: ref_mos_3d(nx,ny,nz)
!
  REAL :: xs(nx)
  REAL :: ys(ny)
  REAL :: zs(nx,ny,nz)
  REAL :: zp(nx,ny,nz)
  REAL :: j3(nx,ny,nz)
  REAL :: pr(nx,ny,nz)
  REAL :: pt(nx,ny,nz)
  REAL :: qv(nx,ny,nz)
  REAL :: ptbar(nx,ny,nz)
  REAL :: qvbar(nx,ny,nz)
  REAL :: rhostr(nx,ny,nz)
  REAL :: qscalar(nx,ny,nz,nscalar)
  REAL :: ceillim(nx,ny)
!
  REAL, PARAMETER :: reflim=90.
  REAL, PARAMETER :: refchek=-10.

!
!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k
  REAL    :: ref,reflow, qr, qc, qi, qs, qh
  INTEGER :: istatus
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!-----------------------------------------------------------------------
!
!  Set minimum ceiling, stored in ceillim
!
!-----------------------------------------------------------------------
!
  IF(ceilopt == 1) THEN

    WRITE(6,'(a,f6.0,a)')                                               &
        ' Setting minimun ceiling to ',ceilmin,' m AGL'
    DO j=1,ny-1
      DO i=1,nx-1
        ceillim(i,j)=ceilmin+zp(i,j,2)
      END DO
    END DO

  ELSE IF(ceilopt == 2) THEN

    WRITE(6,'(a)')                                                      &
        ' Setting minimun ceiling based on surface LCL'
    CALL adaslcl(nx,ny,nz,zs,pr,pt,qv,ceillim)

  ELSE

    WRITE(6,'(a)')                                                      &
        ' Setting zero AGL minimun ceiling'
    DO j=1,ny-1
      DO i=1,nx-1
        ceillim(i,j)=zp(i,j,2)
      END DO
    END DO

  END IF

  DO j=1,ny-1
    DO i=1,nx-1
      DO k=2,nz-1
        IF(ref_mos_3d(i,j,k) > refchek) THEN
           CALL cldset(pr(i,j,k),pt(i,j,k),                             &
                  j3(i,j,k),ptbar(i,j,k),qvbar(i,j,k),rhostr(i,j,k),    &
                  radqvopt,radqcopt,radqropt,radptopt,                  &
                  refsat,rhradobs,refcld,cldrad,                        &
                  refrain,radsetrat,radreflim,radptgain,                &
                  zs(i,j,k),ceillim(i,j),ref_mos_3d(i,j,k),             &
                  qc,qr,qi,qs,qh)

           IF (P_QC > 0) qscalar(i,j,k,P_QC) = qc
           IF (P_QR > 0) qscalar(i,j,k,P_QR) = qr
           IF (P_QI > 0) qscalar(i,j,k,P_QI) = qi
           IF (P_QS > 0) qscalar(i,j,k,P_QS) = qs
           IF (P_QH > 0) qscalar(i,j,k,P_QH) = qh
        END IF
      END DO
    END DO
  END DO
  RETURN
END SUBROUTINE radmcro
!
!##################################################################
!##################################################################
!######                                                      ######
!######                  SUBROUTINE CLDSET                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE cldset(pr,pt,                                                &
           j3,ptbar,qvbar,rhostr,                                       &
           radqvopt,radqcopt,radqropt,radptopt,                         &
           refsat,rhradobs,refcld,cldrad,                               &
           refrain,radsetrat,radreflim,radptgain,                       &
           hgt,ceillim,ref,                                             &
           qc,qr,qi,qs,qh)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Set RH and microphysical variables at a specific point
!  that is above the minimum ceiling and has a reflectivity
!  greater than the threshold.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Dingchen Hou, CAPS
!  May, 1997
!
!  MODIFICATION HISTORY:
!
!  02/02/98 (K. Brewster)
!  Addition of refcldopt and old cldrad fixed cloud setting.
!  Addition of radsetrat to control amount of cloud added.
!  Addition of radreflim to allow limiting the reflectivity
!     used to derive cloud amounts.
!  Slight streamling of code and addition of documentation.
!
!  04/30/2010 (K. Brewster)
!  Modification for new use of ref_mos_3d here and upstream
!  routines.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
!  INPUT:
!
  REAL :: pr,pt,j3,ptbar,qvbar,rhostr
  INTEGER :: radqvopt,radqcopt,radqropt,radptopt
  REAL :: refsat,rhradobs,refcld,cldrad
  REAL :: refrain,radsetrat,radreflim,radptgain
  REAL :: hgt,ceillim,ref
!
!  OUTPUT:
!
  REAL :: qc,qr,qi,qs,qh
!
!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------
!
  REAL :: pratio,t,qv,qvsat,reflect
  REAL :: ref1,ref2,qc1,qc2
  REAL :: dqr,dqc,dpt,dpt1
  REAL :: aaa,bbb,pttem
  INTEGER :: iter
  PARAMETER (ref1=62.96,                                                &
             ref2=7.38,                                                 &
             qc1=2.0E-3,                                                &
             qc2=0.1E-3)
!
!-----------------------------------------------------------------------
!
!  Function f_qvsat and inline directive for Cray PVP
!
!-----------------------------------------------------------------------
!
  REAL :: f_qvsat

!fpp$ expand (f_qvsat)
!!dir$ inline always f_qvsat
!*$*  inline routine (f_qvsat)

!
!-----------------------------------------------------------------------
!
!  Include files
!
!-----------------------------------------------------------------------
!
  INCLUDE 'phycst.inc'
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

  dqc=0.
  dqr=0.
  dpt=0.
!
!
!-----------------------------------------------------------------------
!
!  Limit reflectivity according to radreflim
!
!-----------------------------------------------------------------------
!
  reflect=AMIN1(ref,radreflim)

  dqc=0.0

  IF( ref > refcld .AND. hgt > ceillim) THEN

    IF( radqcopt == 1 ) THEN

    dqc=AMAX1((cldrad-(qc+qi)),0.0)

    ELSE IF( radqcopt == 2 ) THEN

    dqc=(qc1-qc2)/(LOG(ref1)-LOG(ref2))*                                &
        (LOG(reflect)-LOG(ref2))+qc2
    dqc=AMAX1(((radsetrat*dqc)-(qc+qi)),0.0)

  END IF

  END IF
!
!-----------------------------------------------------------------------
!
!  Adjustment of pt due to qc and qr
!
!-----------------------------------------------------------------------
!
  IF( radqropt > 0 .AND. ref > refrain ) THEN
    dqr=(10.0**(reflect/10.0)/17300.0)**(4.0/7.0)                       &
        /(rhostr/j3*1000.0)
    dqr=AMAX1(((radsetrat*dqr)-(qr+qs+qh)),0.0)
  END IF
!
  IF(radptopt == 1 .OR. radptopt == 3) dpt=radptgain*(dqr+dqc)*ptbar/(1.0+qvbar)

  pttem=pt+dpt

!
!-----------------------------------------------------------------------
!
!  Adjustment of water vapor
!
!-----------------------------------------------------------------------
!
  IF(radqvopt > 0 .AND. ref > refsat) THEN
!
!-----------------------------------------------------------------------
!
!  Adjustment of pt due to water vapor change
!  bbb is the bouyancy multiplied by ptbar/g, to be conserved.
!
!-----------------------------------------------------------------------
!
    IF (radptopt == 2 .OR. radptopt == 3) THEN
      pratio=(pr/p0) ** rddcp
      t =pt*pratio
      qvsat=f_qvsat( pr, t )
      qv=qvsat*rhradobs
      aaa=ptbar*(1-0.622)/((0.622+qvbar)*(1.0+qvbar))
      bbb=aaa*qv+pt

      DO iter=1,10
        t=pttem*pratio
        qvsat=f_qvsat( pr, t )
        dpt1=bbb-aaa*(qvsat*rhradobs)-pttem
!        print *,iter,pttem,t,p,qvsat,aaa,bbb,dpt1
        IF (ABS(dpt1) < 1.0E-3) EXIT
        pttem=pttem+dpt1
      END DO
!      101     CONTINUE
      t=pttem*pratio
      qvsat=f_qvsat( pr, t )
      qv=qvsat*rhradobs
    END IF
  END IF

  qc=qc+dqc
  qr=qr+dqr
  pt=pttem

  RETURN
END SUBROUTINE cldset
!
!##################################################################
!##################################################################
!######                                                      ######
!######                   FUNCTION ADASLCL                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
!

SUBROUTINE adaslcl(nx,ny,nz,zs,pr,pt,qv,hgtlcl)

  IMPLICIT NONE
  INTEGER :: nx,ny,nz
  REAL :: zs(nx,ny,nz)
  REAL :: pr(nx,ny,nz)
  REAL :: pt(nx,ny,nz)
  REAL :: qv(nx,ny,nz)
  REAL :: hgtlcl(nx,ny)
!
!-----------------------------------------------------------------------
!
!  Misc local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k
  INTEGER :: imid,jmid

  REAL :: pmb,tc,tk,td,wmr,thepcl,plcl,tlcl
  REAL :: plnhi,plnlo,plnlcl,whi,wlo

  REAL :: oe,wmr2td

!
!-----------------------------------------------------------------------
!
!  Include files
!
!-----------------------------------------------------------------------
!
  INCLUDE 'phycst.inc'
!
!-----------------------------------------------------------------------
!
!  Find pressure of lcl
!  Using pressure of lcl, find height of lcl
!
!-----------------------------------------------------------------------
!
  imid=nx/2
  jmid=ny/2

  DO j=1,ny-1
    DO i=1,nx-1
      pmb=0.01*pr(i,j,2)
      tk=pt(i,j,2)*((1000./pmb)**rddcp)
      tc=tk-273.15
      wmr=1000.*(qv(i,j,2)/(1.-qv(i,j,2)))
      td=wmr2td(pmb,wmr)
      thepcl=oe(tc,td,pmb)
      CALL ptlcl(pmb,tc,td,plcl,tlcl)
      plcl=plcl*100.
      DO k=3,nz-2
        IF(pr(i,j,k) < plcl) EXIT
      END DO
!      81   CONTINUE
      plnhi=LOG(pr(i,j,k))
      plnlo=LOG(pr(i,j,k-1))
      plnlcl=LOG(plcl)
      whi=(plnlo-plnlcl)/(plnlo-plnhi)
      wlo=1.-whi
      hgtlcl(i,j)=whi*zs(i,j,k)+wlo*zs(i,j,k-1)
      IF(i == imid .AND. j == jmid) THEN
        PRINT *, 'arpslcl: ',pmb,tc,td,wmr,plcl,hgtlcl(i,j)
      END IF
    END DO
  END DO
!
  RETURN
END SUBROUTINE adaslcl
