SUBROUTINE barnes_cntrl(nx,ny,mxstn,nvar,maxrng,nobs,x,y,               &
           obs,odiff,obanx,xsta,ysta,isrc,qobs,rely,                    &
           nm_var,npass,kappa,irngsel,range,rngobs,                     &
           wlim,klim,knt,rngwgt,wgtsum,zsum,istat,                      &
           anx,istatus)
!
!  Routine to control Barnes analysis looping.
!  For each pass it calls routines to analyze at grid points,
!  analyze at obs points, computes new obs differences and
!  reports rms statistics.
!
!  Keith Brewster, CAPS
!  March, 1994
!
  IMPLICIT NONE
  INTEGER :: nx,ny,mxstn,nvar,maxrng
  INTEGER :: istat
!
!  Grid arrays
!
  REAL :: x(nx),y(ny)
  REAL :: anx(nx,ny,nvar)
!
!  Observation arrays
!
  INTEGER :: nobs
  REAL :: obs(mxstn,nvar)
  REAL :: odiff(mxstn,nvar)
  REAL :: obanx(mxstn,nvar)
  REAL :: xsta(mxstn),ysta(mxstn)
  REAL :: qobs(mxstn,nvar)
  REAL :: rngobs(mxstn,maxrng)
  INTEGER :: isrc(mxstn)
  INTEGER :: rely(mxstn,nvar)
  CHARACTER (LEN=6) :: nm_var(nvar)
!
!  Barnes control variables
!
  INTEGER :: npass
  INTEGER :: irngsel(nvar)
  REAL :: kappa(npass)
  REAL :: range(nx,ny,maxrng)
  REAL :: wlim
  INTEGER :: klim(npass)
!
!  Barnes scratch space
!
  INTEGER :: knt(nvar)
  REAL :: rngwgt(maxrng)
  REAL :: wgtsum(nvar)
  REAL :: zsum(nvar)
!
  INTEGER :: istatus
!
!  Misc internal variables
!
  INTEGER :: i,j,k,ipass,ista,irng
  INTEGER :: iloc,jloc
  REAL :: delx,dely,rpass
!
!  Initialize rngobs, the appropriate analysis scale range at
!  each observation point.  It is important that this be
!  consistant with the ranges assigned at the grid points.
!  Bilinear interpolation is used rather than an independent
!  calculation.  Obs outside the domain use the nearest grid
!  point in the domain.
!
  delx=x(2)-x(1)
  dely=y(2)-y(1)
  DO ista=1,nobs
    iloc=INT((xsta(ista)-x(1))/delx) + 1
    iloc=MIN(iloc,nx-1)
    iloc=MAX(iloc,1)
    jloc=INT((ysta(ista)-y(1))/dely) + 1
    jloc=MIN(jloc,ny-1)
    jloc=MAX(jloc,1)
    DO irng=1,maxrng
      CALL bilin1(nx,ny,1,xsta(ista),ysta(ista),iloc,jloc,              &
                x,y,range(1,1,irng),rngobs(ista,irng))
    END DO
  END DO
!
!  Establish initial condition.   Analysis is zero and
!  ob differences are equal to obs.
!
  DO k=1,nvar
    DO j=1,ny
      DO i=1,nx
        anx(i,j,k)=0.
      END DO
    END DO
  END DO
!
  DO k=1,nvar
    DO i=1,mxstn
      obanx(i,k)=0.
      odiff(i,k)=obs(i,k)
    END DO
  END DO
!
  DO ipass=1,npass
!
    rpass=kappa(ipass)*kappa(ipass)
!
    CALL bargrid(nx,ny,mxstn,nvar,maxrng,x,y,                           &
                   odiff,xsta,ysta,qobs,isrc,nobs,                      &
                   knt,rngwgt,wgtsum,zsum,                              &
                   rpass,irngsel,range,wlim,klim(ipass),                &
                   anx,istatus)
!
    CALL barobs(mxstn,nvar,maxrng,                                      &
                  odiff,xsta,ysta,qobs,isrc,nobs,                       &
                  knt,rngwgt,wgtsum,zsum,                               &
                  rpass,irngsel,rngobs,wlim,klim(ipass),                &
                  obanx,istatus)
!
    CALL diffnstat(mxstn,nvar,nobs,obs,isrc,obanx,odiff,                &
                       knt,wgtsum,zsum)
!
    WRITE(6,820) ipass
    820   FORMAT('  Statistics for Barnes pass ',i4,/                   &
                 '  ivar  name    knt      bias        rms')
    DO k=1,nvar
      WRITE(6,822) k,nm_var(k),knt(k),wgtsum(k),zsum(k)
      822     FORMAT(i6,2X,a6,i8,f10.3,f10.3)
    END DO
    IF(istat > 0) THEN
      WRITE(istat,820) ipass
      DO k=1,nvar
        WRITE(istat,822) k,nm_var(k),knt(k),wgtsum(k),zsum(k)
      END DO
    END IF
!
  END DO
  istatus=0
  RETURN
END SUBROUTINE barnes_cntrl
!

SUBROUTINE diffnstat(mxstn,nvar,nobs,obs,isrc,obanx,odiff,              &
           knt,dmean,drms)
!
!  Routine to calculate new observation differences
!  and report statistics on the differences.
!
!  Keith Brewster, CAPS
!  March, 1994
!
  IMPLICIT NONE
  INTEGER :: mxstn,nvar
!
  REAL :: obs(mxstn,nvar)
  INTEGER :: isrc(mxstn)
  REAL :: obanx(mxstn,nvar)
  REAL :: odiff(mxstn,nvar)
  INTEGER :: nobs
  INTEGER :: knt(nvar)
  REAL :: dmean(nvar)
  REAL :: drms(nvar)
!
!  Misc internal variables
!
  INTEGER :: i,k
  REAL :: flknt
!
  DO k=1,nvar
    knt(k)=0
    dmean(k)=0.
    drms(k)=0.
  END DO
!
  DO i=1,nobs
    IF(isrc(i) /= 0) THEN
      DO k=1,nvar
        IF(obs(i,k) > -90.) THEN
          odiff(i,k)=obs(i,k)-obanx(i,k)
          knt(k)=knt(k)+1
          dmean(k)=dmean(k)+odiff(i,k)
          drms(k)=drms(k)+odiff(i,k)*odiff(i,k)
        END IF
      END DO
    END IF
  END DO
!
!  Calculate stats from the sums formed above
!
  DO k=1,nvar
    IF(knt(k) > 0) THEN
      flknt=FLOAT(knt(k))
      dmean(k)=dmean(k)/flknt
      drms(k)=SQRT(drms(k)/flknt)
    ELSE
      dmean(k)=0.
      drms(k)=0.
    END IF
  END DO
!
  RETURN
END SUBROUTINE diffnstat

SUBROUTINE barqc(maxsta,nvar,nsrc,nsta,z,zdiff,                         &
           obstime,xsta,ysta,isrc,                                      &
           knt,wgtsum,zsum,sqsum,                                       &
           range,wlim,klim,qclim,varnam,snam,                           &
           iqclist,iqcsave,                                             &
           flagz,flagd)
!
!   Routine to objectively quality control data using a
!   Barnes analysis of nearby stations.
!
!   If an observations is found to be bad it is set to missing flag
!   flagz, and the corresponding difference variable is set to flagz.
!
!   Keith Brewster, April, 1991
!   Based on Barnes analysis    Keith Brewster, March, 1989
!
  IMPLICIT NONE
!
!  Observation Arguments
!
  INTEGER :: maxsta,nvar,nsrc
  REAL :: z(maxsta,nvar),     & ! original variables
   zdiff(maxsta,nvar)      ! variable to analyze
  INTEGER :: obstime(maxsta)
  REAL :: xsta(maxsta),ysta(maxsta)  ! x and y locations of stations
  INTEGER :: isrc(maxsta)
  INTEGER :: nsta
!
!  Scratch space
!
  INTEGER :: knt(nvar)
  REAL :: wgtsum(nvar),zsum(nvar),sqsum(nvar)
!
!  File pointers
!
  INTEGER :: iqclist,iqcsave
!
!  Analysis specification arguments
!
!  INTEGER iwgt
  REAL :: range,    & ! e-folding range km**2 of barnes weight
   wlim      ! limit of weight to set max range
  INTEGER :: klim   ! minimum # of stations to influence grid pt
  REAL :: qclim(nsrc,nvar)   ! QC cutoffs
!
!  Diagnostics and other argument variables
!
  CHARACTER (LEN=6) :: varnam(nvar)
  CHARACTER (LEN=5) :: snam(maxsta)
  REAL :: flagz,flagd
!
!  Misc internal variables
!
  REAL :: rlimsq,dxsta,dysta,dist,                                      &
       wgt,zanx,dob,sqanx,chklim
  INTEGER :: jsta,ksta,ivar
  LOGICAL :: listit,saveit
!
!  Set printing logicals
!
  listit=(iqclist > 0)
  saveit=(iqcsave > 0)
!
!  Based on the minimum weight to consider, set
!  a maximum range
!
  rlimsq=-range*ALOG(wlim)
!
!  print *, ' Analyzing...'
!  print *, '    rlim= ', rlim
!  print *, '    range = ',range
!
!   Uses Barnes weighting function.  see variable wgt.
!
!   Ksta is location that is being examined
!
  DO ksta=nsta,1,-1
    DO ivar=1,nvar
      zsum(ivar)=0.
      sqsum(ivar)=0.
      wgtsum(ivar)=0.
      knt(ivar)=0
    END DO
    DO jsta=1,nsta
      IF(jsta /= ksta) THEN
        dxsta=xsta(ksta)-xsta(jsta)
        dysta=ysta(ksta)-ysta(jsta)
        dist=dxsta*dxsta + dysta*dysta
        IF(dist < rlimsq) THEN
          wgt=EXP(-dist/range)
          DO ivar=1,nvar
            IF(zdiff(jsta,ivar) > flagd) THEN
              knt(ivar)=knt(ivar)+1
              wgtsum(ivar)=wgtsum(ivar)+wgt
              zsum(ivar)=zsum(ivar)+wgt*zdiff(jsta,ivar)
              sqsum(ivar)=sqsum(ivar)+                                  &
                   wgt*zdiff(jsta,ivar)*zdiff(jsta,ivar)
            END IF
          END DO
        END IF
      END IF  ! dont use current station
    END DO
    DO ivar=1,nvar
      IF(zdiff(ksta,ivar) > flagd) THEN
        IF(knt(ivar) > klim) THEN
          sqanx=sqsum(ivar)-                                            &
                (zsum(ivar)*zsum(ivar)/wgtsum(ivar))
          sqanx=AMAX1(sqanx,0.)
          sqanx=3.0*SQRT(sqanx/wgtsum(ivar))
          zanx=zsum(ivar)/wgtsum(ivar)
          dob=zdiff(ksta,ivar)-zanx
          chklim=AMAX1(sqanx,qclim(isrc(ksta),ivar))
          IF(ABS(dob) > chklim) THEN
            PRINT *, ' QC flagging data...sta=',snam(ksta)
            PRINT *, '     var, ob =',                                  &
                        varnam(ivar),z(ksta,ivar)
            PRINT *, '     delta ob =',dob
            PRINT *, '  3.0*stdv,qclimit=',sqanx,                       &
                                 qclim(isrc(ksta),ivar)
            IF(listit) WRITE(iqclist,810)                               &
                snam(ksta),obstime(ksta),ivar,varnam(ivar),             &
                z(ksta,ivar),dob,                                       &
                sqanx,qclim(isrc(ksta),ivar)
            810           FORMAT(2X,a5,i5,2X,i3,1X,a6,f11.2,f11.2,f11.2,f11.2)
            z(ksta,ivar)=flagz
            zdiff(ksta,ivar)=flagd
          END IF
          zsum(ivar)=dob
        ELSE
!          print *, ' Not enuf data to check ivar= ',
!    +                   ivar,' for ',snam(ksta)
          zsum(ivar)=-999.
        END IF
      ELSE
        zsum(ivar)=-999.
      END IF  ! data missing, dont bother checking
    END DO
    IF(saveit) WRITE(iqcsave,820) snam(ksta),obstime(ksta),             &
              (zdiff(ksta,ivar),ivar=1,nvar)
    820   FORMAT(2X,a4,i6,2X,15(f7.2))
  END DO
  RETURN
END SUBROUTINE barqc

SUBROUTINE barobs(mxstn,nvar,maxrng,                                    &
           odiff,xsta,ysta,qobs,isrc,nobs,                              &
           knt,rngwgt,wgtsum,zsum,                                      &
           rpass,irngsel,rngobs,wlim,klim,                              &
           obanx,istatus)
!
!  Routine to objectively analyze data at the observation locations
!  Purpose is to create an analysis parallel to the gridded analysis
!  for obtaining a new observation difference for the next pass.
!  Barnes weight function is employed.
!
!  Keith Brewster, March, 1989
!  Streamlined to remove unneeded options, kb, April, 1990
!  Modified for use in OLAPS  Keith Brewster, March, 1994
!
  IMPLICIT NONE
!
!  Observation Arguments
!
  INTEGER :: mxstn,nvar,maxrng
  REAL :: odiff(mxstn,nvar),             & ! variable to analyse
   xsta(mxstn),ysta(mxstn),  & ! x and y locations of stations
   qobs(mxstn,nvar)
  INTEGER :: isrc(mxstn)
  INTEGER :: nobs
!
!    Scratch space
!
  INTEGER :: knt(nvar)
  REAL :: rngwgt(maxrng)
  REAL :: wgtsum(nvar),zsum(nvar)
!
!  Input/Output analysis at observation locations.Grid arguments
!
  REAL :: obanx(mxstn,nvar)    ! analysed data --  Input and OUTPUT
  INTEGER :: istatus
!
!  Analysis specification arguments
!
!  integer iwgt
  INTEGER :: irngsel(nvar)
  REAL :: rpass,                                                        &
       rngobs(mxstn,maxrng),    & ! e-folding range km**2 of barnes weight
   wlim      ! limit in km of station influence
  INTEGER :: klim   ! minimum # of stations to influence grid pt
!
!  Misc internal variables
!
  INTEGER :: ista,jsta,ivar,irng
  REAL :: denom,rlimsq,dxsta,dysta,dist,wgt,rmax
!  real wvar
!
!  print *, ' Analyzing...'
!  print *, '    rlim= ', rlim
!  print *, '    range = ',range
!
!   Uses Barnes weighting function.  see variable wgt.
!
  DO ista=1,nobs
    IF(isrc(ista) /= 0) THEN
      rmax=rngobs(ista,1)
      DO irng=1,maxrng
        rmax=AMAX1(rmax,rngobs(ista,irng))
      END DO
      denom=rpass*rmax
      rlimsq=-denom*ALOG(wlim)
      DO ivar=1,nvar
        zsum(ivar)=0.
        wgtsum(ivar)=0.
        knt(ivar)=0
      END DO
      DO jsta=1,nobs
        IF(isrc(jsta) /= 0) THEN
          dxsta=xsta(ista)-xsta(jsta)
          dysta=ysta(ista)-ysta(jsta)
          dist=dxsta*dxsta + dysta*dysta
          IF(dist < rlimsq) THEN
            DO irng=1,maxrng
              rngwgt(irng)=                                             &
                     EXP(-dist/(rpass*rngobs(ista,irng)))
            END DO
            DO ivar=1,nvar
              wgt=rngwgt(irngsel(ivar))
              IF(odiff(jsta,ivar) > -90.) THEN
                knt(ivar)=knt(ivar)+1
!                wvar=wgt/qobs(jsta,ivar)
                wgtsum(ivar)=wgtsum(ivar)+wgt
                zsum(ivar)=zsum(ivar)+wgt*odiff(jsta,ivar)
              END IF
            END DO
          END IF
        END IF
      END DO
      DO ivar=1,nvar
        IF(knt(ivar) > klim)                                            &
            obanx(ista,ivar)=obanx(ista,ivar)+zsum(ivar)/wgtsum(ivar)
      END DO
    END IF  ! valid station?
  END DO
  istatus=0
  RETURN
END SUBROUTINE barobs

SUBROUTINE bargrid(nx,ny,mxstn,nvar,maxrng,x,y,                         &
           obs,xsta,ysta,qobs,isrc,nobs,                                &
           knt,rngwgt,wgtsum,zsum,                                      &
           rpass,irngsel,range,wlim,klim,                               &
           anx,istatus)
!
!   Routine to objectively analyze data on a regular grid specified
!   by x and y.  Barnes weight function is employed.
!
!   Keith Brewster, March, 1989
!   Streamlined to remove unneeded options, kb, April, 1990
!   Modified for use in OLAPS  Keith Brewster, March, 1994
!
  IMPLICIT NONE
!
!  Observation Arguments
!
  INTEGER :: mxstn,nvar,maxrng
  REAL :: obs(mxstn,nvar),             & ! variable to analyse
   xsta(mxstn),ysta(mxstn),  & ! x and y locations of stations
   qobs(mxstn,nvar)
  INTEGER :: isrc(mxstn)
  INTEGER :: nobs
!
!  Scratch space
!
  INTEGER :: knt(nvar)
  REAL :: rngwgt(maxrng)
  REAL :: wgtsum(nvar),zsum(nvar)
!
!  Grid arguments
!
  INTEGER :: nx,ny
  REAL :: x(nx),y(ny)        ! x and y locations of grid points
  REAL :: anx(nx,ny,nvar)    ! analysed grid --  Input and OUTPUT
  INTEGER :: istatus
!
!  Analysis specification arguments
!
!  integer iwgt
  INTEGER :: irngsel(nvar)
  REAL :: rpass,                                                        &
       range(nx,ny,maxrng),  & ! e-folding range km**2 of barnes weight
   wlim                  ! limit in km of station influence
  INTEGER :: klim   ! minimum # of stations to influence grid pt
!
!  Misc internal variables
!
  REAL :: rlimsq,denom,dxsta,dysta,dist,wgt,rmax
  INTEGER :: i,j,jsta,ivar,irng
!
!   Uses Barnes weighting function.  see variable wgt.
!
  DO j=1,ny
    DO i=1,nx
      rmax=range(i,j,1)
      DO irng=1,maxrng
        rmax=AMAX1(rmax,range(i,j,irng))
      END DO
      denom=rpass*rmax
      rlimsq=-denom*ALOG(wlim)
      DO ivar=1,nvar
        zsum(ivar)=0.
        wgtsum(ivar)=0.
        knt(ivar)=0
      END DO
      DO jsta=1,nobs
        IF(isrc(jsta) /= 0) THEN
          dxsta=x(i)-xsta(jsta)
          dysta=y(j)-ysta(jsta)
          dist=dxsta*dxsta + dysta*dysta
          IF(dist < rlimsq) THEN
            DO irng=1,maxrng
              rngwgt(irng)=EXP(-dist/(rpass*range(i,j,irng)))
            END DO
            DO ivar=1,nvar
              wgt=rngwgt(irngsel(ivar))
              IF(obs(jsta,ivar) > -90.) THEN
                knt(ivar)=knt(ivar)+1
                wgtsum(ivar)=wgtsum(ivar)+wgt
                zsum(ivar)=zsum(ivar)+wgt*obs(jsta,ivar)
              END IF
            END DO
          END IF
        END IF
      END DO
      DO ivar=1,nvar
        IF(knt(ivar) > klim) anx(i,j,ivar)=anx(i,j,ivar)+zsum(ivar)/wgtsum(ivar)
      END DO
    END DO
  END DO
  istatus=0
  RETURN
END SUBROUTINE bargrid

SUBROUTINE bilin1(nx,ny,nvar,xpos,ypos,i,j,x,y,var,varint)
!
!  Does bilinear interpolation of variables (nvar variables)
!  to position (xpos,ypos) given a grid of variables "var", their
!  x and y positions of the grid (x,y) and the index of the
!  grid point (i,j) which is to the "lower-left" of xpos and ypos.
!  That means (xpos,ypos) is between (i,j) and (i+1,j+1).
!
!  A saftey feature in case xpos and ypos are outside the grid
!  is built-in.
!
  IMPLICIT NONE
!
!  Arguments, input
!
  INTEGER :: nx,ny,nvar
  REAL :: var(nx,ny,nvar)
  REAL :: x(nx),y(ny),xpos,ypos
  INTEGER :: i,j
!
!  Arguments output
!
  REAL :: varint(nvar)
!
!  Misc Internal Variables
!
  REAL :: xtem,ytem,dxpos,dypos,dx,dy
  REAL :: w11,w12,w21,w22
  INTEGER :: ivar,ip1,jp1
!
  xtem=AMAX1(xpos,x(1))
  xtem=AMIN1(xtem,x(nx))
  ytem=AMAX1(ypos,y(1))
  ytem=AMIN1(ytem,y(ny))
  ip1=i+1
  jp1=j+1
  dx=x(ip1)-x(i)
  dy=y(jp1)-y(j)
  dxpos=xtem-x(i)
  dypos=ytem-y(j)
  IF((dx /= 0.) .AND. (dy /= 0.))THEN
    dx=dxpos/dx
    dy=dypos/dy
    w11=(1.-dx)*(1.-dy)
    w21=dx*(1.-dy)
    w12=dy*(1.-dx)
    w22=dx*dy
  ELSE IF(dx /= 0.) THEN
    dx=dxpos/dx
    w11=1.-dx
    w21=dx
    w12=0.
    w22=0.
  ELSE
    dy=dypos/dy
    w11=1.-dy
    w21=0.
    w12=dy
    w22=0.
  END IF
!
  DO ivar=1,nvar
    varint(ivar)=var(  i,  j,ivar)*w11                                  &
                +var(ip1,  j,ivar)*w21                                  &
                +var(  i,jp1,ivar)*w12                                  &
                +var(ip1,jp1,ivar)*w22
  END DO
  RETURN
END SUBROUTINE bilin1
