!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE CELTRK                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE celtrk(nx,ny,nz,ibeg,iend,jbeg,jend,kbeg,kend,               &
           w, qr, rhobar, x,y,zp, ntrcell, xcw, ycw,                    &
           ireturn,                                                     &
           ihead,itail,jseg,kcell,icmin,icmax,                          &
           xs,ys,tem1,tem2)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Driver for cell tracking subroutines.
!  Part of the ARPS cell tracking subsystem.
!
!  AUTHOR: Keith Brewster
!   22 April 1993
!
!-----------------------------------------------------------------------
!
!  INPUT VARIABLES
!
!    x        x coordinate of grid points in physical/comp. space (m)
!    y        y coordinate of grid points in physical/comp. space (m)
!    zp       z coordinate of grid points in physical space (m)
!
!    w        vertical component of velocity in Cartesian
!             coordinates (m/s).
!
!    qr       Rain water mixing ratio (kg/kg)
!
!    rhobar   Mean density (kg/m^3)
!
!  OUTPUT:
!
!  Writes tracking information to a file, and returns the following
!  variables:
!
!  ntrcell    Number of cells found at this time
!  xcw        Weighted central x position of all cells at this time
!  ycw        Weighted central y position of all cells at this time
!
!  ireturn    return status
!             =0  center of "mass" successfully determined
!             =-1 no cells found so center of "mass" only estimated based
!                 on previous location
!
!  WORK ARRAYS:
!
!    tem1     Temporary work array.
!    tem2     Temporary work array.
!
!  (These arrays are defined and used locally (i.e. inside this
!   subroutine), they may also be passed into routines called by
!   this one. Exiting the call to this subroutine, these temporary
!   work arrays may be used for other purposes, and therefore their
!   contents may be overwritten. Please examine the usage of work
!   arrays before you make any change to the code.)
!
!-----------------------------------------------------------------------
!
!  Variable Declarations:
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'grid.inc'          ! Grid parameters
  INCLUDE 'mp.inc'
!
!-----------------------------------------------------------------------
!
!  Arguments
!
!-----------------------------------------------------------------------
!
  INTEGER :: nx,ny,nz
  INTEGER :: ibeg,iend,jbeg,jend,kbeg,kend

  REAL :: x(nx)
  REAL :: y(ny)
  REAL :: zp(nx,ny,nz)

  REAL :: w     (nx,ny,nz)
  REAL :: qr    (nx,ny,nz)
  REAL :: rhobar(nx,ny,nz)

  INTEGER :: ntrcell
  REAL    :: xcw,ycw
  INTEGER :: ireturn
!
!-----------------------------------------------------------------------
!
!  Work Arrays
!
!-----------------------------------------------------------------------
!
  REAL :: tem1(nx,ny,nz)
  REAL :: tem2(nx,ny,nz)
  REAL :: xs(nx),ys(ny)
!
!-----------------------------------------------------------------------
!
!  Array size parameters.
!
!  maxcell: maximum number of cells
!           not all may be present at one time
!
!  maxtime: number of cell locations over time to save
!
!-----------------------------------------------------------------------
!
  INTEGER :: maxcell,maxtime
  PARAMETER (maxcell=26, maxtime=15)
!
!-----------------------------------------------------------------------
!
!  Important adjustable parameters.
!
!  maxhgt:  maximum height to look for cells.
!  thresh:  minimum value of vertical velocity used to define a cell
!     eps:  small margin of error to allow to pass as meeting threshold
!           if thresh is not exactly met
!  maxblw:  maximum number of points meeting thresh-eps threshold but
!           not thresh threshold in a segment
!  minlen:  minimum length of a segment (meters)
!    minarea:  minimum area of a horizontal cell (m*m)
!  minvol:  minimum volume of a 3-d cell (m*m*m)
!   minqr:  minimum rainwater maximum within a 3-d cell (kg/kg)
!  disthr:  maximum distance from an extrapolated cell location
!           to match with a current cell location (m)
!
!
!-----------------------------------------------------------------------
!
  INTEGER :: maxblw
  REAL    :: thresh,eps,maxhgt,minlen,minqr,minarea,minvol,disthr
  PARAMETER (maxhgt=12000.,    & ! m
             thresh=5.0,       & ! m/s
                eps=0.5,       & ! m/s
             maxblw=2,         & ! grid points
             minlen=2000.,     & ! m
            minarea=9.0E06,    & ! m*m
             minvol=50.0E09,   & ! m*m*m
              minqr=0.001,     & ! kg/kg maximum in cell
             disthr=4000.)       ! m
!
!
!-----------------------------------------------------------------------
!
!  Variables to determine cell location at a given time.
!
!-----------------------------------------------------------------------
!

  INTEGER :: ihead(nx*ny),itail(nx*ny),jseg(nx*ny),kcell(nx*ny)

  REAL :: xcent(maxcell,2)
  REAL :: ycent(maxcell,2)
  REAL :: zcent(maxcell,2)
  REAL :: qrmax(maxcell,2)
  REAL :: wgt(maxcell,2)
  REAL :: area(maxcell,2)
  INTEGER :: icmin(maxcell,ny,2)
  INTEGER :: icmax(maxcell,ny,2)
  INTEGER :: ncell
!
!-----------------------------------------------------------------------
!
!  Location-to-motion variables
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=1) :: cname(maxcell)
  INTEGER :: ktime,ncalls
  INTEGER :: nctim(maxtime)
  REAL    :: alltime(maxtime)
  REAL    :: xpos(maxtime,maxcell)
  REAL    :: ypos(maxtime,maxcell)
  REAL    :: zpos(maxtime,maxcell)
  REAL    :: xnot(maxcell)
  REAL    :: ynot(maxcell)
  REAL    :: allvol(maxtime,maxcell)
  REAL    :: allwgt(maxtime,maxcell)
  REAL    :: umotion(maxtime,maxcell)
  REAL    :: vmotion(maxtime,maxcell)
  REAL    :: wmotion(maxtime,maxcell)
  REAL    :: growth(maxtime,maxcell)
!
!-----------------------------------------------------------------------
!
!  Save cname and all variables that are a function of time.
!
!-----------------------------------------------------------------------
!
  SAVE cname
  SAVE ktime,ncalls
  SAVE nctim
  SAVE alltime,xpos,ypos,zpos,xnot,ynot
  SAVE allvol,allwgt,umotion,vmotion,wmotion,growth

  DATA cname /'A','B','C','D','E','F','G','H','I','J',                  &
              'K','L','M','N','O','P','Q','R','S','T',                  &
              'U','V','W','X','Y','Z'/
  DATA ktime,ncalls /0,0/
!
!
!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------
!
  REAL :: prtmis
  PARAMETER (prtmis= 999.)
  CHARACTER (LEN=1) :: comma
  PARAMETER (comma=',')
  REAL :: wgtsum,xsum,ysum
  INTEGER :: maxseg
  INTEGER :: istr,ifin,jstr,jfin,kstr,kfin
  INTEGER :: i,j,kc,mtime,mpast
!
  CHARACTER (LEN=256) :: trkfn  ! Name of max./min. file
  INTEGER :: ltrkfn
  INTEGER :: istat, nchtrk
  SAVE nchtrk
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  IF (mp_opt > 0) THEN
    WRITE(6,'(1x,a,/,7x,a)')                                            &
       'Error: Cell-tracking is still not implemented for MPI mode.',   &
       'Program is aborting ...'
    CALL arpsstop('MPI mode is still not implemented for cell-tracking algorithm.',1)
  END IF

  maxseg = nx*ny

  IF( ((nx*nz)/2) < maxcell) THEN
    WRITE(6,'(1x,a,/1x,a,I4,a,/1x,a)')                                  &
    'The size of temporary arrays passed to CELTRK not big enough& !',  &
    '(nx*nz)/2 must be greater than ',maxcell,                          &
    ' for CELTRK to work properly.',                                    &
    ' No cell tracking was done.'
    RETURN
  END IF

  IF(ncalls == 0) THEN

    ncalls=1

    IF (myproc == 0) THEN
      trkfn = runname(1:lfnkey)//'.track'
      ltrkfn = lfnkey + 6

      CALL getunit(nchtrk)

      WRITE(6,'(1x,a,a,a/,1x,a)')                                       &
        'Check to see if file ',trkfn(1:ltrkfn),' already exists.',     &
        'If so, append a version number to the filename.'

      CALL fnversn( trkfn, ltrkfn )

      OPEN(nchtrk,FORM='formatted',STATUS='new',                        &
                  FILE=trkfn(1:ltrkfn),IOSTAT=istat)

      IF( istat /= 0) THEN

        WRITE(6,'(/a,i2,/a/)')                                          &
          ' Error occured when opening file '//trkfn(1:ltrkfn)//        &
          ' using FORTRAN unit ',nchtrk,' Program stopped in trkdriv.'
        CALL arpsstop ("arpstop called from celtrk",1)

      END IF

      WRITE(nchtrk,'(a)') ''''//runname//''''

      WRITE(nchtrk,'(/a,a)')                                            &
        '   time  name  xloc   yloc   zloc',                            &
        '    u      v      w     vol       wgt  '
    END IF
  END IF
!
!-----------------------------------------------------------------------
!
!    Check the limits that were passed in.
!    For bookkeeping reasons (avoid array out of bounds and some
!    storage tricks), as well as physical sense, the range of i,j,k
!    should not involve the artificial boundary points.
!
!    Use istr,ifin,jstr,jfin,kstr,kfin in place of ibeg,iend, etc.
!    for the call to loccell.
!
!-----------------------------------------------------------------------
!
  istr=MAX(ibeg,2)
  ifin=MIN(iend,nx-1)
  jstr=MAX(jbeg,2)
  jfin=MIN(jend,ny-1)
  kstr=MAX(kbeg,2)
  kfin=MIN(kend,nz-1)
!
!-----------------------------------------------------------------------
!
!  Locate cells in grid.
!
!-----------------------------------------------------------------------
!
  DO i=1,nx-1
    xs(i)=(x(i)+x(i+1))*0.5
  END DO
  xs(nx)=2.*xs(nx-1)-xs(nx-2)
  DO j=1,ny-1
    ys(j)=(y(j)+y(j+1))*0.5
  END DO
  ys(ny)=2.*ys(ny-1)-ys(ny-2)

  CALL loccell(nx,ny,nz, maxseg,maxcell,ncell,                          &
              istr,ifin,jstr,jfin,kstr,kfin,                            &
              w,qr, rhobar, xs,ys,zp,                                   &
              ihead,itail,kcell,jseg,                                   &
              xcent,ycent,zcent,icmin,icmax,                            &
              qrmax,area,wgt,maxhgt,                                    &
              thresh,eps,maxblw,minlen,minqr,minarea,minvol,            &
              ireturn, tem1,tem2)
!
!-----------------------------------------------------------------------
!
!  Add movement of grid to the cell locations.
!  Compute a weighted mean location of all cells at this time.
!
!-----------------------------------------------------------------------
!
  xsum=0.
  ysum=0.
  wgtsum=0.
  DO kc=1,ncell
    xcent(kc,1)=xcent(kc,1)+xgrdorg
    ycent(kc,1)=ycent(kc,1)+ygrdorg
    wgtsum=wgtsum+wgt(kc,1)
    xsum=xsum+xcent(kc,1)*wgt(kc,1)
    ysum=ysum+ycent(kc,1)*wgt(kc,1)
  END DO
!
  ntrcell=ncell
!
  IF(wgtsum > 0.) THEN
    ireturn=0
    xcw=xsum/wgtsum
    ycw=ysum/wgtsum
  ELSE
    ireturn=-1
    IF( ktime > 1) THEN
      xcw=xcw+umove*tceltrk
      ycw=ycw+vmove*tceltrk
    ELSE
      xcw=(xs(1)+xs(nx))*0.5
      ycw=(ys(1)+ys(ny))*0.5
    END IF
  END IF
!
!-----------------------------------------------------------------------
!
!  Spit out some diagnostics to the list file.
!
!-----------------------------------------------------------------------
!
  IF (myproc == 0) THEN
    WRITE(6,'(//a,i4,/a,a)')                                            &
      ' Cell positions for whole vol  Total cells: ',ncell,             &
      ' kc   xcabs   ycabs    xcrel    ycrel  ',                        &
      '  zcent     vol       wgt     qrmax'
    DO kc=1,ncell
      WRITE(6,'(1x,i4,2(F9.2,F9.2),F9.2,F9.0,F12.0)')                   &
           kc,(xcent(kc,1)*0.001),(ycent(kc,1)*0.001)                   &
             ,((xcent(kc,1)-xgrdorg)*0.001)                             &
             ,((ycent(kc,1)-ygrdorg)*0.001)                             &
             ,(zcent(kc,1)*0.001)                                       &
             ,( area(kc,1)*1.0E-09),wgt(kc,1)                           &
             ,(qrmax(kc,1)*1000.)
    END DO
  END IF
!
!-----------------------------------------------------------------------
!
!  Match current cell volumes with previous.
!
!-----------------------------------------------------------------------
!
!  print *, ' Calling matctim'
  CALL matctim(maxcell,maxtime,ktime,mtime,mpast,                       &
              disthr,curtim,umove,vmove,                                &
              ncell,xcent(1,1),ycent(1,1),zcent(1,1),                   &
              area(1,1),wgt(1,1),                                       &
              cname,nctim,alltime,                                      &
              xpos,ypos,zpos,xnot,ynot,                                 &
              allvol,allwgt,                                            &
              umotion,vmotion,wmotion)

  IF (myproc == 0) THEN
    WRITE(6,'(//a,i4,/a)')                                              &
      ' Cells combined in time arrays: ',nctim(mtime),                  &
      '   jc   xpos    ypos   zpos   volume     wgt'
    DO kc=1,nctim(mtime)
      IF(xpos(mtime,kc) > -1.e30) THEN
        WRITE(6,'(3x,a1,F9.2,F9.2,F9.2,F9.0,F12.0)')                    &
           cname(kc),(xpos(mtime,kc)*0.001),                            &
              (ypos(mtime,kc)*0.001),                                   &
              (zpos(mtime,kc)*0.001),                                   &
            (allvol(mtime,kc)*1.0E-09),                                 &
             allwgt(mtime,kc)
      END IF
    END DO
  END IF
!
!-----------------------------------------------------------------------
!
!  Calculate cell motion vector and growth rate.
!
!-----------------------------------------------------------------------
!
  CALL getvec(maxcell,maxtime,ktime,mtime,                              &
              cname,nctim,alltime,                                      &
              xpos,ypos,zpos,xnot,ynot,                                 &
              allvol,allwgt,                                            &
              umotion,vmotion,wmotion,growth)
!
!-----------------------------------------------------------------------
!
!  Spit out motion vector results.
!
!-----------------------------------------------------------------------
!
  IF (myproc == 0) THEN
    IF(nctim(mtime) > 0) THEN
      DO kc=1,nctim(mtime)
        IF(xpos(mtime,kc) > -1.0E30) THEN
          IF(umotion(mtime,kc) > -1.0E30) THEN
            WRITE(nchtrk,'(1x,F7.0,2x,a1,1x,f7.2,f7.2,f7.2,f7.2,f7.2,f7.2,f7.0,f10.0)') &
              alltime(mtime),cname(kc),                                 &
              (xpos(mtime,kc)*0.001),                                   &
              (ypos(mtime,kc)*0.001),                                   &
              (zpos(mtime,kc)*0.001),                                   &
              umotion(mtime,kc),                                        &
              vmotion(mtime,kc),                                        &
              wmotion(mtime,kc),                                        &
              (allvol(mtime,kc)*1.0E-09),                               &
              allwgt(mtime,kc)
          ELSE
            WRITE(nchtrk,'(1x,F7.0,2x,a1,1x,f7.2,f7.2,f7.2,f7.2,f7.2,f7.2,f7.0,f10.0)') &
              alltime(mtime),cname(kc),                                 &
              (xpos(mtime,kc)*0.001),                                   &
              (ypos(mtime,kc)*0.001),                                   &
              (zpos(mtime,kc)*0.001),                                   &
                  prtmis,                                               &
                  prtmis,                                               &
                  prtmis,                                               &
              (allvol(mtime,kc)*1.0E-09),                               &
              allwgt(mtime,kc)
          END IF
        END IF
      END DO
    END IF
    WRITE(6,'(/a,i8,f9.0,f9.0//)') ' End of tracking for time:',        &
             mtime,alltime(mtime),curtim
  END IF

  RETURN
END SUBROUTINE celtrk
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE LOCCELL                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE loccell(nx,ny,nz, maxseg,maxcell,ncell,                      &
           ibeg,iend,jbeg,jend,kbeg,kend,                               &
           tvar,qr,rhobar,x,y,zp,                                       &
           ihead,itail,kcell,jseg,                                      &
           xcent,ycent,zcent,icmin,icmax,                               &
           qrmax,area,wgt,maxhgt,                                       &
           thresh,eps,maxblw,minlen,minqr,minarea,minvol,               &
           ireturn, tem1,tem2)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Builds cells from tracking variable (tvar, i.e. w) field and
!  establishes a central postion for each which is used for
!  storm tracking.
!
!  A module of the ARPS cell tracking subsystem.
!
!  AUTHOR: Keith Brewster
!   21 Mar 1993
!
!-----------------------------------------------------------------------
!
!  DATA ARRAYS:
!
!    x        x coordinate of grid points in physical/comp. space (m)
!    y        y coordinate of grid points in physical/comp. space (m)
!    zp       z coordinate of grid points in physical space (m)
!
!    w        z component of velocity (m/s)
!
!    qr       Rainwater mixing ratio (kg/kg)
!
!    rhobar   Mean density (kg/m^3)
!
!  WORK ARRAYS:
!
!    tem1     Temporary work array.
!    tem2     Temporary work array.
!
!  (These arrays are defined and used locally (i.e. inside this
!   subroutine), they may also be passed into routines called by
!   this one. Exiting the call to this subroutine, these temporary
!   work arrays may be used for other purposes, and therefore their
!   contents may be overwritten. Please examine the usage of work
!   arrays before you make any change to the code.)
!
!-----------------------------------------------------------------------
!
!  Variable Declarations:
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
  INTEGER :: nx,ny,nz
  INTEGER :: maxseg
  INTEGER :: maxcell
  INTEGER :: ncell
  INTEGER :: ibeg,iend,jbeg,jend,kbeg,kend
  REAL :: x(nx)
  REAL :: y(ny)
  REAL :: zp(nx,ny,nz)
  REAL :: tvar(nx,ny,nz)
  REAL :: qr(nx,ny,nz)
  REAL :: rhobar(nx,ny,nz)
  REAL :: thresh
  INTEGER :: ihead(maxseg)
  INTEGER :: itail(maxseg)
  INTEGER :: jseg(maxseg)
  INTEGER :: kcell(maxseg)
  REAL :: xcent(maxcell,2)
  REAL :: ycent(maxcell,2)
  REAL :: zcent(maxcell,2)
  REAL :: wgt(maxcell,2)
  REAL :: qrmax(maxcell,2)
  REAL :: area(maxcell,2)
  INTEGER :: icmin(maxcell,ny,2)
  INTEGER :: icmax(maxcell,ny,2)
  INTEGER :: maxblw
  REAL :: eps,maxhgt,minlen,minqr,minarea,minvol
  REAL :: tem1(nx,ny,nz)
  REAL :: tem2(nx,ny,nz)
  INTEGER :: ireturn
!
!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: nseg,onvf
  INTEGER :: i,j,k,kc,ncellk
  REAL :: tmin,tmax
  INTEGER :: ilmax,jlmax,klmax,ilmin,jlmin,klmin
  INTEGER :: icengd,jcengd,iseg,ihd,itl,kntblw
  REAL :: xlen,thmin,zsfc
  LOGICAL :: onseg
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!
!-----------------------------------------------------------------------
!
!  Average qr in the vertical to get qr at w points.
!  This is stored in tem1.
!
!-----------------------------------------------------------------------
!
  onvf=1
  CALL avgz(qr, onvf, nx,ny,nz,                                         &
             1,nx-1,1,ny-1,2,nz-1, tem1)

!
!-----------------------------------------------------------------------
!
!  Average rhobar in the vertical to get rhobar at w points.
!  This is stored in tem2.
!
!-----------------------------------------------------------------------
!
  onvf=1
  CALL avgz(rhobar, onvf, nx,ny,nz,                                     &
             1,nx-1,1,ny-1,2,nz-1, tem2)
!
!-----------------------------------------------------------------------
!
!  Find the w max and report it.
!  For now this information is only used for diagnostics,
!  but in the future the cell matching could begin at the level
!  of max w and work up and down from there.
!  For now, it begins at the top and works down.
!
!-----------------------------------------------------------------------
!
  CALL a3dmax(tvar, 1,nx,ibeg,iend, 1,ny,jbeg,jend,                     &
                    1,nz,kbeg,kend, tmax,tmin,                          &
                    ilmax,jlmax,klmax,ilmin,jlmin,klmin)
  PRINT *, ' w max = ',tmax
!  print *, ' imax,jmax,kmax = ',ilmax,jlmax,klmax
!  print *, ' thresh= ',thresh
!  print *, ' minlen= ',minlen
!
!-----------------------------------------------------------------------
!
!  Search for segments, strings of points meeting the threshold.
!  The variable onseg is a flag indicating whether we are on a
!  segment looking for its end (onseg=true) or looking for the
!  beginning of the next segment (onseg=false).
!
!  When the end of a segment is reached matchsg is called which
!  matches the current segment to any previously found segments
!  on this level.
!
!  Variables ihead and itail describe the beginning and end of
!  the each segment.
!
!  At the end of each level, matchsg is called in reverse
!  order to see if any further matching of cells can be made
!  (since the matching process can be directional dependent).
!  In the following example segments A and B would not be matched
!  to each other going from 1-to-3 but on the reverse matching
!  (3-to-1) B would be matched to C which had originally been
!  matched to A.
!
!  3      DDDD     EEEE
!  2        CCCCCCCCCCCCC
!  1         AAA      BBBB
!
!-----------------------------------------------------------------------
!
  icengd=((nx-3)/2) + 1
  jcengd=((ny-3)/2) + 1
  zsfc=0.5*(zp(icengd,jcengd,1)+zp(icengd,jcengd,2))
  thmin=thresh-eps
  ncell=0
  DO k=kend,kbeg,-1
!    print *, ' k = ',k
    IF((zp(icengd,jcengd,k)-zsfc) <= maxhgt) THEN
      nseg=0
      ncellk=0
      IF(tmax > thresh) THEN
        DO j=jbeg,jend
          onseg=.false.
          DO i=ibeg,iend
            IF(onseg) THEN
              IF(tvar(i,j,k) < thresh) THEN
                kntblw=kntblw+1
                IF(tvar(i,j,k) >= thmin) THEN
                  IF(kntblw > maxblw) THEN
                    itl=(i-kntblw)
                    xlen=ABS(x(itl)-x(ihd-1))
                    IF(xlen >= minlen) THEN
                      IF(nseg < maxseg) THEN
                        nseg=nseg+1
                        ihead(nseg)=ihd
                        itail(nseg)=itl
                        jseg(nseg)=j
                        kcell(nseg)=-99
                        CALL matchsg(maxseg,maxcell,                    &
                            nseg, 1,(nseg-1),1,ihd,itl,j,               &
                            ihead,itail,jseg,kcell,ncellk)
                      END IF
                    END IF
                    onseg=.false.
                  END IF
                ELSE  ! LT.thresh and LT.thmin
                  itl=(i-kntblw)
                  xlen=ABS(x(itl)-x(ihd-1))
                  IF(xlen >= minlen) THEN
                    IF(nseg < maxseg) THEN
                      nseg=nseg+1
                      ihead(nseg)=ihd
                      itail(nseg)=itl
                      jseg(nseg)=j
                      kcell(nseg)=-99
                      CALL matchsg(maxseg,maxcell,                      &
                          nseg, 1,(nseg-1),1,ihd,itl,j,                 &
                          ihead,itail,jseg,kcell,ncellk)
                    END IF
                  END IF
                  onseg=.false.
                END IF
!
!-----------------------------------------------------------------------
!
!  Onseg and tvar is greater than thresh
!  Reset kntblw to zero,
!  then resume search for the end of the segment.
!
!-----------------------------------------------------------------------
!
              ELSE
                kntblw=0
              END IF
!
!-----------------------------------------------------------------------
!
!  Not "onseg"
!  See if the threshold is exceeded, if so, initialize
!  segment pointers, ihd and kntblw, and set onseg to true.
!
!-----------------------------------------------------------------------
!
            ELSE
              IF(tvar(i,j,k) >= thresh) THEN
                ihd=i
                kntblw=0
                onseg=.true.
              END IF
            END IF
          END DO
!
!-----------------------------------------------------------------------
!
!  If you've reached the end of the i indices and your are
!  still "onseg", then close this segment.
!
!-----------------------------------------------------------------------
!
          IF(onseg) THEN
            itl=iend-kntblw
            xlen=ABS(x(itl)-x(ihd-1))
            IF(xlen >= minlen) THEN
              IF(nseg < maxseg) THEN
                nseg=nseg+1
                ihead(nseg)=ihd
                itail(nseg)=itl
                jseg(nseg)=j
                kcell(nseg)=-99
                CALL matchsg(maxseg,maxcell,                            &
                    nseg, 1,(nseg-1),1,ihd,itl,j,                       &
                    ihead,itail,jseg,kcell,ncellk)
              END IF
            END IF
            onseg=.false.
          END IF
        END DO
      END IF     ! tmax .gt. thresh
!
!-----------------------------------------------------------------------
!
!  Do backwards matching of segments
!
!-----------------------------------------------------------------------
!
!    print *, '  before backward matching ncellk = ',ncellk
      DO iseg=(nseg-1),1,-1
        CALL matchsg(maxseg,maxcell,iseg, nseg,(iseg+1),-1,             &
                     ihead(iseg),itail(iseg),jseg(iseg),                &
                     ihead,itail,jseg,kcell,ncellk)
      END DO
!    print *, '  after backward matching ncell = ',ncellk
!
!-----------------------------------------------------------------------
!
!  Calculate cell areas and centroids.
!
!-----------------------------------------------------------------------
!
      IF(ncellk > 0) THEN
        CALL ctrwgt(nx,ny,maxseg,maxcell,nseg,ncellk,                   &
                  tvar(1,1,k),tem1(1,1,k),tem2(1,1,k),                  &
                  x,y,zp(1,1,k),                                        &
                  ihead,itail,jseg,kcell,                               &
                  xcent(1,2),ycent(1,2),zcent(1,2),                     &
                  qrmax(1,2),area(1,2),wgt(1,2))
!
!-----------------------------------------------------------------------
!
!  Create array describing min and max of i for each j column
!  in each cell.
!
!  This is required for efficient vertical matching of horizontal
!  cells to create 3d cells.
!
!-----------------------------------------------------------------------
!
        DO j=1,ny
          DO kc=1,maxcell
            icmin(kc,j,2)=nx+1
            icmax(kc,j,2)=-99
          END DO
        END DO
        DO iseg=1,nseg
          kc=kcell(iseg)
          icmin(kc,jseg(iseg),2)=                                       &
                         MIN(icmin(kc,jseg(iseg),2),ihead(iseg))
          icmax(kc,jseg(iseg),2)=                                       &
                         MAX(icmax(kc,jseg(iseg),2),itail(iseg))
        END DO
!
!-----------------------------------------------------------------------
!
!  Clean up cell list, eliminate those with small area.
!  Or those not meeting the minimum rainwater (qr) requirement.
!
!-----------------------------------------------------------------------
!
!      print *, ' inside loccell minqr,minarea= ',minqr,minarea
!      print *, ' Before sieve2d, ncell= ',ncellk
        CALL sieve2d(maxcell,ny,ncellk,                                 &
                    xcent(1,2),ycent(1,2),zcent(1,2),                   &
                    qrmax(1,2),area(1,2),wgt(1,2),                      &
                    icmin(1,1,2),icmax(1,1,2),                          &
                    minarea)
!      print *, ' After sieving, ncell= ',ncellk
!
!  tell us about icmin and icmax
!      do 256 kc=1,ncellk
!       print *, ' k, kc = ',k,kc
!      do 255 j=1,ny
!        print *,' j= ',j,' imin= ',icmin(kc,j,2),
!    :                       ' imax= ',icmax(kc,j,2)
!255       continue
!256       continue
      END IF
!
!    IF( ncellk.gt.0 ) THEN
!      write(6,'(//a,i4,a,i4,/a)')
!    : '  Cells on level: ',k,'  Total cells: ',ncellk,
!    : '   kc       xcent       ycent      zcent      area        wgt'
!      DO 300 kc=1,ncellk
!        write(6,'(1x,i4,F9.2,F9.2,F9.2,F9.0,F12.0)')
!    :      kc,xcent(kc,2),ycent(kc,2),zcent(kc,2),area(kc,2),wgt(kc,2)
!00       CONTINUE
!    END IF
!
!-----------------------------------------------------------------------
!
!  Combine cells  with cells found at previous levels.
!  The combined cell info is put in the level one index.
!
!-----------------------------------------------------------------------
!
      CALL linkcell(maxcell,nx,ny,nz,k,ncell,ncellk,                    &
                  x,y,zp,xcent,ycent,zcent,icmin,icmax,                 &
                  qrmax,area,wgt)
    END IF    ! below maxhgt?
  END DO
!
!-----------------------------------------------------------------------
!
!  Compute the central x,y,z
!  for each cell volume from the info collected on plane 1.
!
!-----------------------------------------------------------------------
!
  DO kc=1,ncell
    IF(wgt(kc,1) /= 0.0) THEN
      xcent(kc,1)=xcent(kc,1)/wgt(kc,1)
      ycent(kc,1)=ycent(kc,1)/wgt(kc,1)
      zcent(kc,1)=zcent(kc,1)/wgt(kc,1)
    END IF
  END DO
!
!-----------------------------------------------------------------------
!
!  Eliminate cells which are too small.
!
!-----------------------------------------------------------------------
  PRINT *, ' before sieve3d, ncell= ',ncell
  CALL sieve3d(maxcell,ncell,                                           &
               xcent(1,1),ycent(1,1),zcent(1,1),                        &
               qrmax(1,1),area(1,1),wgt(1,1),                           &
               minqr,minvol)
  PRINT *, ' after sieve3d, ncell= ',ncell
!
  RETURN
END SUBROUTINE loccell
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE MATCHSG                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE matchsg(maxseg,maxcell,                                      &
           kseg, mbeg,mend,mdir,ihd,itl,j,                              &
           ihead,itail,jseg,kcell,ncell)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Matches line segments to form horizontal cells.
!  A module of the ARPS cell tracking subsystem.
!
!  AUTHOR: Keith Brewster
!   21 Mar 1993
!
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!
!  Variable Declarations:
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
  INTEGER :: maxseg,maxcell,kseg
  INTEGER :: ihd,itl,j
  INTEGER :: mbeg,mend,mdir
  INTEGER :: ihead(maxseg),itail(maxseg),jseg(maxseg),                  &
          kcell(maxseg)
  INTEGER :: ncell
!
!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: mseg
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!
!-----------------------------------------------------------------------
!
!  Form cells from line segments
!
!-----------------------------------------------------------------------
!
  DO mseg=mbeg,mend,mdir
    IF(ABS(jseg(mseg)-j) == 1) THEN
      IF(ihd >= ihead(mseg) .AND. ihd <= itail(mseg)) GO TO 120
      IF(itl >= ihead(mseg) .AND. itl <= itail(mseg)) GO TO 120
      IF(itl >= itail(mseg) .AND. ihd <= ihead(mseg)) GO TO 120
    END IF
  END DO
!
!-----------------------------------------------------------------------
!
!  No cell found so this may be a new one.
!
!-----------------------------------------------------------------------
!
  IF(kcell(kseg) < 0) THEN
    IF(ncell < maxcell) THEN
      ncell=ncell+1
      kcell(kseg)=ncell
!      print *, ' kseg j,ihd,itl: ',kseg,j,ihd,itl
!      print *, ' new cell: ',ncell,kcell(kseg)
    ELSE
      PRINT *, ' WARNING: ran out of space for cells.'
    END IF
  END IF
  RETURN
!
!-----------------------------------------------------------------------
!
!  Cell found, it belongs with kcell(mseg)
!
!-----------------------------------------------------------------------
!
  120   CONTINUE
  IF(kcell(mseg) > 0) THEN
    kcell(kseg)=kcell(mseg)
  ELSE
    IF(ncell < maxcell) THEN
      ncell=ncell+1
      kcell(kseg)=ncell
      kcell(mseg)=ncell
    END IF
  END IF
!  print *, ' kseg j,ihd,itl: ',kseg,j,ihd,itl
!  print *, ' matched with cell: ',kcell(mseg),j,
!    :             ihead(kcell(mseg)),itail(kcell(mseg))
  RETURN
END SUBROUTINE matchsg
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE CTRWGT                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE ctrwgt(nx,ny,maxseg,maxcell,nseg,ncell,                      &
           tvar,qratw,rhoatw,x,y,zp,                                    &
           ihead,itail,jseg,kcell,                                      &
           xcent,ycent,zcent,qrmax,area,wgt)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Compute a central location for all cells formed by matchsg.
!  A module of the ARPS cell tracking subsystem.
!
!  AUTHOR: Keith Brewster
!   21 Mar 1993
!
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!
!  Variable Declarations:
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
  INTEGER :: nx,ny
  INTEGER :: maxseg
  INTEGER :: maxcell
  REAL :: tvar(nx,ny)
  REAL :: qratw(nx,ny)
  REAL :: rhoatw(nx,ny)
  REAL :: x(nx)
  REAL :: y(ny)
  REAL :: zp(nx,ny)
  INTEGER :: ihead(maxseg)
  INTEGER :: itail(maxseg)
  INTEGER :: jseg(maxseg)
  INTEGER :: kcell(maxseg)
  REAL :: xcent(maxcell)
  REAL :: ycent(maxcell)
  REAL :: zcent(maxcell)
  REAL :: qrmax(maxcell)
  REAL :: area(maxcell)
  REAL :: wgt(maxcell)
  INTEGER :: nseg
  INTEGER :: ncell
!
!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------
!
  REAL :: ddx,ddy,ddxy,wayt
  INTEGER :: i,j,m,kc
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  ddx=x(2)-x(1)
  ddy=y(2)-y(1)
  ddxy=ABS(ddx*ddy)
!
!-----------------------------------------------------------------------
!
!  Initializations.
!
!-----------------------------------------------------------------------
!
  DO kc=1,ncell
    wgt(kc)=0.
    qrmax(kc)=-999.
    area(kc)=0.
    xcent(kc)=0.
    ycent(kc)=0.
    zcent(kc)=0.
  END DO
!
!-----------------------------------------------------------------------
!
!  Go through all segments.  Build sums for each cell.
!
!-----------------------------------------------------------------------
!
  DO m=1,nseg
    IF(kcell(m) > 0) THEN
      j=jseg(m)
      kc=kcell(m)
      DO i=ihead(m),itail(m)
        wayt=rhoatw(i,j)*tvar(i,j)*tvar(i,j)
!         wayt=tvar(i,j)*tvar(i,j)
        wgt(kc)=wgt(kc)+wayt
        xcent(kc)=xcent(kc)+wayt*x(i)
        ycent(kc)=ycent(kc)+wayt*y(j)
        zcent(kc)=zcent(kc)+wayt*zp(i,j)
        qrmax(kc)=AMAX1(qrmax(kc),qratw(i,j))
        area(kc)=area(kc)+ddxy
      END DO
    ELSE
      PRINT *, ' Warning: segment ',m,' not assigned to a cell'
    END IF
  END DO
!
!-----------------------------------------------------------------------
!
!  Calculated cell centroid, a weight average of the locations
!  of all the grid points in the cell.
!
!-----------------------------------------------------------------------
!
  DO kc=1,ncell
    IF(wgt(kc) > 0.) THEN
      xcent(kc)=xcent(kc)/wgt(kc)
      ycent(kc)=ycent(kc)/wgt(kc)
      zcent(kc)=zcent(kc)/wgt(kc)
    END IF
  END DO
!
  RETURN
END SUBROUTINE ctrwgt
!
!##################################################################
!##################################################################
!######                                                      ######
!######               SUBROUTINE SIEVE2D                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE sieve2d(maxcell,ny,ncell,                                    &
           xcent,ycent,zcent,qrmax,area,wgt,                            &
           icmin,icmax,                                                 &
           minarea)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Discard cells less than a given areal size.
!  A module of the ARPS cell tracking subsystem.
!
!  AUTHOR: Keith Brewster
!   21 Mar 1993
!
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!
!  Variable Declarations:
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
  INTEGER :: maxcell,ny,ncell
  REAL :: xcent(maxcell)
  REAL :: ycent(maxcell)
  REAL :: zcent(maxcell)
  REAL :: qrmax(maxcell)
  REAL :: area(maxcell)
  REAL :: wgt(maxcell)
  INTEGER :: icmin(maxcell,ny)
  INTEGER :: icmax(maxcell,ny)
  REAL :: minarea
!
!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: ic,j,kc,nn
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  nn=ncell
  DO ic=1,ncell
    20      CONTINUE
    IF(ic > nn) EXIT
!    print *, ' ic area, qrmax= ',ic,area(ic),qrmax(ic)
    IF(area(ic) < minarea) THEN
      nn=nn-1
      DO kc=ic,nn
        wgt(kc)=wgt(kc+1)
        area(kc)=area(kc+1)
        qrmax(kc)=qrmax(kc+1)
        xcent(kc)=xcent(kc+1)
        ycent(kc)=ycent(kc+1)
        zcent(kc)=zcent(kc+1)
      END DO
      DO j=1,ny
        DO kc=ic,nn
          icmin(kc,j)=icmin(kc+1,j)
          icmax(kc,j)=icmax(kc+1,j)
        END DO
      END DO
      GO TO 20
    END IF
  END DO
  ncell=nn
  RETURN
END SUBROUTINE sieve2d
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######               SUBROUTINE SIEVE3D                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE sieve3d(maxcell,ncell,                                       &
           xcent,ycent,zcent,qrmax,area,wgt,                            &
           minqr,minvol)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Discard 3-dimensional cells less than a given volume size.
!  A module of the ARPS cell tracking subsystem.
!
!  AUTHOR: Keith Brewster
!   21 Mar 1993
!
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!
!  Variable Declarations:
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
  INTEGER :: maxcell
  INTEGER :: ncell
  REAL :: xcent(maxcell)
  REAL :: ycent(maxcell)
  REAL :: zcent(maxcell)
  REAL :: qrmax(maxcell)
  REAL :: area(maxcell)
  REAL :: wgt(maxcell)
  REAL :: minqr,minvol
!
!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: ic,kc,nn
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  nn=ncell
  DO ic=1,ncell
    20      CONTINUE
    IF(ic > nn) EXIT
    IF(area(ic) < minvol .OR. qrmax(ic) < minqr) THEN
      nn=nn-1
      DO kc=ic,nn
        wgt(kc)=wgt(kc+1)
        area(kc)=area(kc+1)
        qrmax(kc)=qrmax(kc+1)
        xcent(kc)=xcent(kc+1)
        ycent(kc)=ycent(kc+1)
        zcent(kc)=zcent(kc+1)
      END DO
      GO TO 20
    END IF
  END DO
  ncell=nn
  RETURN
END SUBROUTINE sieve3d
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE LINKCELL                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE linkcell(maxcell,nx,ny,nz,k,ncell,ncellk,                    &
           x,y,zp,xcent,ycent,zcent,icmin,icmax,                        &
           qrmax,area,wgt)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Link cells on adjacent horizontal surfaces to form 3-d cell.
!  Combine location info to get 3D location info.
!  A module of the ARPS cell tracking subsystem.
!
!  AUTHOR: Keith Brewster
!   21 Mar 1993
!
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!
!  Variable Declarations:
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
  INTEGER :: maxcell,nx,ny,nz
  INTEGER :: k,ncell,ncellk
  REAL :: x(nx)
  REAL :: y(ny)
  REAL :: zp(nx,ny,nz)
  REAL :: xcent(maxcell,2)
  REAL :: ycent(maxcell,2)
  REAL :: zcent(maxcell,2)
  REAL :: wgt(maxcell,2)
  REAL :: qrmax(maxcell,2)
  REAL :: area(maxcell,2)
  INTEGER :: icmin(maxcell,ny,2)
  INTEGER :: icmax(maxcell,ny,2)
!
!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: ic,jc,j,iloc,jloc
  REAL :: dx,dy
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  NOTE:
!  3-d info is in 3rd index=1,
!  2-d info for a single level is in 3rd index=2.
!
!  Find a cell that has an overlap with this cell.
!  For this nested loop, 2 is the plane being checked
!                     1 is the plane being searched for possible matches
!                    ic is the cell being checked
!                    jc is the cell being compared on plane 1 (3-d info)
!
!-----------------------------------------------------------------------
!
  dx=x(2)-x(1)
  dy=y(2)-y(1)
  DO ic=1,ncellk
    DO jc=1,ncell
!
!-----------------------------------------------------------------------
!
!  Check for overlap between these two cells
!
!-----------------------------------------------------------------------
!
      DO j=1,ny
        IF(icmax(ic,j,2) > 0 .AND. icmax(jc,j,1) > 0) THEN
          IF((icmax(ic,j,2) >= icmin(jc,j,1)) .AND.                     &
              (icmin(ic,j,2) <= icmax(jc,j,1)))                         &
              GO TO 250
        END IF
      END DO
    END DO
!
!-----------------------------------------------------------------------
!
!  No match found, must be a new cell
!
!-----------------------------------------------------------------------
!
    ncell=ncell+1
    jc=ncell
!
    wgt(jc,1)=wgt(ic,2)
    xcent(jc,1)=wgt(ic,2)*xcent(ic,2)
    ycent(jc,1)=wgt(ic,2)*ycent(ic,2)
    zcent(jc,1)=wgt(ic,2)*zcent(ic,2)
    qrmax(jc,1)=qrmax(ic,2)
    iloc=nint((xcent(ic,2)/dx) + 1.5)
    jloc=nint((ycent(ic,2)/dy) + 1.5)
    area(jc,1)=area(ic,2)*                                              &
               0.5*(zp(iloc,jloc,k+1)-zp(iloc,jloc,k-1))
!
!-----------------------------------------------------------------------
!
!    Set icmin, icmax
!
!    1 index is added to the max and 1 index is subtracted from
!    the min so that the comparisons in do loop 180 accept cells
!    that touch rather than just overlap.  Adding now saves time
!    of adding inside the IF statement.
!
!-----------------------------------------------------------------------
!
    DO j=1,ny
      icmin(jc,j,1)=icmin(ic,j,2)-1
      icmax(jc,j,1)=icmax(ic,j,2)+1
    END DO
!
    CYCLE
!
!-----------------------------------------------------------------------
!
!    Match found.
!    Update weighted sums
!
!-----------------------------------------------------------------------
!
    250     CONTINUE
    wgt(jc,1)=wgt(jc,1)+wgt(ic,2)
    xcent(jc,1)=xcent(jc,1)+wgt(ic,2)*xcent(ic,2)
    ycent(jc,1)=ycent(jc,1)+wgt(ic,2)*ycent(ic,2)
    zcent(jc,1)=zcent(jc,1)+wgt(ic,2)*zcent(ic,2)
    qrmax(jc,1)=AMAX1(qrmax(jc,1),qrmax(ic,2))
    iloc=nint((xcent(ic,2)/dx) + 1.5)
    jloc=nint((ycent(ic,2)/dy) + 1.5)
    area(jc,1)=area(jc,1)+area(ic,2)*                                   &
               0.5*(zp(iloc,jloc,k+1)-zp(iloc,jloc,k-1))
!
!-----------------------------------------------------------------------
!
!    Match found.
!    Update icmin,icmax
!
!-----------------------------------------------------------------------
!
    DO j=1,ny
      icmin(jc,j,1)=MIN(icmin(jc,j,1),(icmin(ic,j,2)-1))
      icmax(jc,j,1)=MAX(icmax(jc,j,1),(icmax(ic,j,2)+1))
    END DO
  END DO
!
  RETURN
END SUBROUTINE linkcell
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######               SUBROUTINE MATCTIM                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE matctim(maxcell,maxtime,ktime,mtime,mpast,                   &
           disthr,time,umove,vmove,                                     &
           ncell,xcent,ycent,zcent,                                     &
           area,wgt,                                                    &
           cname,nctim,alltime,                                         &
           xpos,ypos,zpos,xnot,ynot,                                    &
           allvol,allwgt,                                               &
           umotion,vmotion,wmotion)
!
  IMPLICIT NONE
  INTEGER :: maxcell,maxtime,ktime,mtime
  REAL :: time,disthr,umove,vmove
  INTEGER :: ncell
  REAL :: xcent(maxcell),ycent(maxcell),zcent(maxcell),                 &
       wgt(maxcell),area(maxcell)
  CHARACTER (LEN=2) :: cname(maxcell)
  INTEGER :: nctim(maxtime)
  REAL :: alltime(maxtime)
  REAL :: xpos(maxtime,maxcell),ypos(maxtime,maxcell),                  &
       zpos(maxtime,maxcell)
  REAL :: xnot(maxcell),ynot(maxcell)
  REAL :: allvol(maxtime,maxcell),allwgt(maxtime,maxcell)
  REAL :: umotion(maxtime,maxcell),vmotion(maxtime,maxcell),            &
       wmotion(maxtime,maxcell)
!
!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: icell,jcell,inear,mpast
  REAL :: xxtrap,yxtrap
  REAL :: d2thr,deltat,ddx,ddy,dist2,distmin,d2rat,drmin
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!-----------------------------------------------------------------------
!
!  Index control:   ktime is count of times processed here
!                   mtime is current index in arrays saving locations
!                   mpast is past index in arrays saving locations
!
!-----------------------------------------------------------------------
!
  ktime=ktime+1
  mtime=MOD(ktime,maxtime)
  IF(mtime == 0) mtime=maxtime
  mpast=MOD((ktime-1),maxtime)
  IF(mpast == 0) mpast=maxtime
  PRINT *, ' inside matctim '
  PRINT *, ' ncell,ktime,mtime = ',ncell,ktime,mtime
  PRINT *, ' mpast,nctim(mpast)= ',mpast,nctim(mpast)
!
  d2thr=disthr*disthr
  alltime(mtime)=time
!
  DO icell=1,maxcell
    nctim(mtime)=0
    xpos(mtime,icell)=-1.0E30
    ypos(mtime,icell)=-1.0E30
    zpos(mtime,icell)=-1.0E30
    allvol(mtime,icell)=-1.0E30
    allwgt(mtime,icell)=-1.0E30
    umotion(mtime,icell)=-1.0E30
    vmotion(mtime,icell)=-1.0E30
    wmotion(mtime,icell)=-1.0E30
  END DO
!
  IF(ktime == 1) THEN
    nctim(mtime)=ncell
    DO icell=1,ncell
      xpos(mtime,icell)=xcent(icell)
      ypos(mtime,icell)=ycent(icell)
      zpos(mtime,icell)=zcent(icell)
      allvol(mtime,icell)=area(icell)
      allwgt(mtime,icell)=wgt(icell)
    END DO
  ELSE      ! not first time, need to match up with past cells
!
!-----------------------------------------------------------------------
!
!    Try to match location of "current" cells with extrapolated
!    positions of past cells.  It is assumed throughout that we are
!    dealing with absolute (not grid relative locations).  This
!    is necessary for the linear fit to work right under the condition
!    that the grid motion (umove,vmove) are changing with time.
!
!
!    For extrapolation use, (in the following order of priority,
!     according to the available data for each cell):
!
!   1) Least-squares line from previous u,v fitting
!      (umotion,vmotion and xnot,ynot)
!
!   2) Last umotion and vmotion and last location
!
!   3) Guess motion (zero relative to grid, e.g.)
!      and last location
!
!-----------------------------------------------------------------------
!
    nctim(mtime)=nctim(mpast)
    deltat=time-alltime(mpast)
    DO jcell=1,nctim(mpast)
      IF(umotion(mpast,jcell) /= -1.0E30 .AND.                          &
            vmotion(mpast,jcell) /= -1.0E30 ) THEN
        PRINT *, ' Past abs umotion,vmotion (m/s):',                    &
                  umotion(mpast,jcell),                                 &
                  vmotion(mpast,jcell)
        IF(xnot(jcell) /= -1.0E30 .AND. ynot(jcell) /= -1.0E30 ) THEN
          xxtrap=xnot(jcell)+time*umotion(mpast,jcell)
          yxtrap=ynot(jcell)+time*vmotion(mpast,jcell)
        ELSE
          xxtrap=xpos(mpast,jcell)+deltat*umotion(mpast,jcell)
          yxtrap=ypos(mpast,jcell)+deltat*vmotion(mpast,jcell)
        END IF
      ELSE
        xxtrap=xpos(mpast,jcell)+deltat*umove
        yxtrap=ypos(mpast,jcell)+deltat*vmove
      END IF
      PRINT *, ' past x,xxtrap = ',xpos(mpast,jcell),xxtrap
      PRINT *, ' past y,yxtrap = ',ypos(mpast,jcell),yxtrap
!
!-----------------------------------------------------------------------
!
!    Find the location among the current cells (icell) that is the
!    closest to this past cell (jcell).  When a cell has been used
!    already, its wgt is set to -1.
!
!-----------------------------------------------------------------------
!
      inear=0
      drmin=1.0E30
      DO icell=1,ncell
        IF(wgt(icell) > 0.) THEN
          ddx=xcent(icell)-xxtrap
          ddy=ycent(icell)-yxtrap
          dist2=ddx*ddx + ddy*ddy
          IF(dist2 < d2thr) THEN
            d2rat=dist2/area(icell)
!            d2rat=dist2/wgt(icell)
            IF(d2rat < drmin) THEN
              inear=icell
              drmin=d2rat
            END IF
          END IF
        END IF
      END DO
      IF(inear > 0) THEN
        PRINT *, ' time match jcell,inear= ',jcell,inear
        distmin=SQRT(drmin*area(inear))
        PRINT *, ' min distance (m) = ',distmin
        xpos(mtime,jcell)=xcent(inear)
        ypos(mtime,jcell)=ycent(inear)
        zpos(mtime,jcell)=zcent(inear)
        allvol(mtime,jcell)=area(inear)
        allwgt(mtime,jcell)=wgt(inear)
        wgt(inear)=-1.
      END IF
    END DO
!
!-----------------------------------------------------------------------
!
!    Now make new locations for any cells not matched in the
!    above process.
!
!-----------------------------------------------------------------------
!
    DO icell=1,ncell
      IF(wgt(icell) > 0.) THEN
        PRINT *, ' New cell in time/motion matching: ', icell
        IF(nctim(mtime) < maxcell) THEN
          nctim(mtime)=nctim(mtime)+1
          PRINT *, ' New count is ',nctim(mtime)
          jcell=nctim(mtime)
          xpos(mtime,jcell)=xcent(icell)
          ypos(mtime,jcell)=ycent(icell)
          zpos(mtime,jcell)=zcent(icell)
          allvol(mtime,jcell)=area(icell)
          allwgt(mtime,jcell)=wgt(icell)
        ELSE
          PRINT *, ' WARNING ran out space for cells'
          PRINT *, '    in subroutine matctim. '
        END IF
      END IF
    END DO
  END IF

!
  PRINT *, ' Leaving matctim, nctim = ',nctim(mtime)
  RETURN
END SUBROUTINE matctim
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINE GETVEC                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE getvec(maxcell,maxtime,ktime,mtime,                          &
           cname,nctim,alltime,                                         &
           xpos,ypos,zpos,xnot,ynot,                                    &
           allvol,allwgt,                                               &
           umotion,vmotion,wmotion,growth)
!
  IMPLICIT NONE
  INTEGER :: maxcell,maxtime,ktime,mtime
  CHARACTER (LEN=2) :: cname(maxcell)
  INTEGER :: nctim(maxtime)
  REAL :: alltime(maxtime)
  REAL :: xpos(maxtime,maxcell),ypos(maxtime,maxcell),                  &
       zpos(maxtime,maxcell)
  REAL :: xnot(maxcell),ynot(maxcell)
  REAL :: allvol(maxtime,maxcell),allwgt(maxtime,maxcell)
  REAL :: umotion(maxtime,maxcell),vmotion(maxtime,maxcell),            &
       wmotion(maxtime,maxcell),growth(maxtime,maxcell)
!
!-----------------------------------------------------------------------
!
!  Parameters and variables for least squares fitting.
!
!-----------------------------------------------------------------------
!
  REAL :: twparm
  PARAMETER (twparm=900.*900.)    ! secs*secs
  REAL :: tvalid,ufit,vfit
  INTEGER :: ireturn
!
!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: icell,mpast,itime
  REAL :: deltat,znot,dummy
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  itime=MIN(ktime,maxtime)
  mpast=MOD((ktime-1),maxtime)
  IF(mpast == 0) mpast=maxtime
  PRINT *, '  Inside getvec: ktime= ',ktime
  PRINT *, '     mtime ',mtime,'  nctim(mtime) = ',nctim(mtime)
  PRINT *, '     mpast ',mpast,'  nctim(mpast) = ',nctim(mpast)
!
!-----------------------------------------------------------------------
!
!  Initialize all motions to a wacky missing value.
!
!-----------------------------------------------------------------------
!
  DO icell=1,maxcell
    xnot(icell)=-1.0E30
    ynot(icell)=-1.0E30
    umotion(mtime,icell)=-1.0E30
    vmotion(mtime,icell)=-1.0E30
    wmotion(mtime,icell)=-1.0E30
    growth(mtime,icell)=-1.0E30
  END DO
!
!-----------------------------------------------------------------------
!
!   Use the least squares fit of a line to "itime" number of
!   points to get the velocity.
!
!   In the event a solution could not be found, ireturn.ne.0,
!   use the current position less the past position to get the
!   velocity.
!
!-----------------------------------------------------------------------
!
  IF(ktime > 1) THEN
    deltat=alltime(mtime)-alltime(mpast)
    tvalid=0.5*(alltime(mtime)+alltime(mpast))
    DO icell=1,nctim(mtime)
      IF(xpos(mtime,icell) /= -1.0E30) THEN
        CALL leastsq(xpos(1,icell),alltime,                             &
            tvalid,twparm,itime,xnot(icell),ufit,ireturn)
        IF(ireturn == 0) THEN
          umotion(mtime,icell)=ufit
        ELSE IF(xpos(mpast,icell) /= -1.0E30) THEN
          umotion(mtime,icell)=                                         &
              (xpos(mtime,icell)-xpos(mpast,icell))/deltat
        END IF
        CALL leastsq(ypos(1,icell),alltime,                             &
            tvalid,twparm,itime,ynot(icell),vfit,ireturn)
        IF(ireturn == 0) THEN
          vmotion(mtime,icell)=vfit
        ELSE IF(ypos(mpast,icell) /= -1.0E30) THEN
          vmotion(mtime,icell)=                                         &
              (ypos(mtime,icell)-ypos(mpast,icell))/deltat
        END IF
        CALL leastsq(zpos(1,icell),alltime,                             &
                 tvalid,twparm,itime,                                   &
                 znot,wmotion(mtime,icell),ireturn)
        IF(ireturn /= 0 .AND. zpos(mpast,icell) /= -1.0E30) THEN
          wmotion(mtime,icell)=                                         &
              (zpos(mtime,icell)-zpos(mpast,icell))/deltat
        END IF
        CALL leastsq(allvol(1,icell),alltime,                           &
                 tvalid,twparm,itime,                                   &
                 dummy,growth(mtime,icell),ireturn)
        IF(ireturn /= 0 .AND. allvol(mpast,icell) /= -1.0E30) THEN
          growth(mtime,icell)=                                          &
              (allvol(mtime,icell)-allvol(mpast,icell))/deltat
        END IF
      END IF
    END DO
  END IF
  RETURN
END SUBROUTINE getvec
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE LEASTSQ                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE leastsq(xo,TO,tvalid,twparm,npts,xnot,uo,ireturn)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Modified least-squares estimate for the following equations:
!
!                  X = xnot + Uo*t
!
!  Modification is that what is minimized is not the X deviation
!  squared, rather X deviation squared times a weight specified
!  in the time weighting function tweight.  The idea is to
!  make the line fit the observations near the time specified as
!  the "valid time", tvalid.
!
!-----------------------------------------------------------------------
!
!  AUTHOR:  Steven Lazarus
!           4/12/93
!
!  MODIFICATIONS:
!    4/13/93  Keith Brewster
!             Added weighting based on time.
!             Changed to do one variables at a time.
!    4/19/93  Keith Brewster
!             Added ireturn variable.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    npts     Number of points in the Sample
!    xo       Observed storm position (predictand) for X = xnot + Uo*t
!    to       Observed time for xo(to),yo(to) (t=to)
!    tvalid   "Valid" time of estimation, used for time weighting
!             Closer fit to data near vaild time.
!    twparm   Time weighting parameter used as denominator in exponential
!             time weighting function  (secs*secs)
!
!   OUTPUT:
!
!    xnot   Origin of regression line    X =Xnot + Uo*t, (f(xo,to))
!    uo     Slope of the regression line X =Xnot + Uo*t, (f(xo,to))
!    ireturn  Return status indicator.
!
!-----------------------------------------------------------------------
!
!   Variable Declarations
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: npts       ! Number of points in the sample

  REAL :: xo  (npts)    ! Observed independent variable
  REAL :: TO  (npts)    ! Observed time for xo,yo
  REAL :: tvalid        ! valid time
  REAL :: twparm        ! time weighting parameter (secs*secs)
  REAL :: xnot          ! Least-Squares fit for x=xnot+uo*t
  REAL :: uo            ! Least-Squares fit for x=xnot+uo*t
  INTEGER :: ireturn    ! Return status indicator
                        ! =0  succesful completion
                        ! =1  not enough data
                        ! =2  divide by zero avoided, nothing computed

!
!-----------------------------------------------------------------------
!
!   Misc. Local Variables
!
!-----------------------------------------------------------------------
!
  REAL :: sumx,sumt,sumw,sumt2
  REAL :: sumxt
  REAL :: dt,twgt
  REAL :: denomx,denomu

  INTEGER :: i,knt
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!   Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  knt    = 0
  sumw   = 0.0
  sumx   = 0.0
  sumt   = 0.0
  sumxt  = 0.0
  sumt2  = 0.0

  DO i=1,npts
    IF(xo(i) /= -1.0E30) THEN
      knt=knt+1
      dt=tvalid-TO(i)
      twgt=EXP(-dt*dt/twparm)
      sumw  = sumw  + twgt
      sumx  = sumx  + xo(i)*twgt
      sumt  = sumt  + TO(i)*twgt
      sumt2 = sumt2 + TO(i)*TO(i)*twgt
      sumxt = sumxt + xo(i)*TO(i)*twgt
    END IF
  END DO

  IF( knt > 2) THEN
    denomx = sumt*sumt - sumt2*sumw
    denomu = sumw*sumt2 - sumt*sumt
    IF( denomx /= 0.  .AND.  denomu /= 0.) THEN
      ireturn=0
      xnot  = (sumt*sumxt - sumt2*sumx) / denomx
      uo    = (sumw*sumxt - sumx*sumt)  / denomu
    ELSE
      ireturn=2
      xnot  = -1.e30
      uo    = -1.e30
    END IF
  ELSE
    ireturn=1
    xnot  = -1.e30
    uo    = -1.e30
  END IF

  RETURN
END SUBROUTINE leastsq
