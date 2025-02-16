!
!##################################################################
!##################################################################
!######                                                      ######
!######               SUBROUTINE SET_LBCOPT                  ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE set_lbcopt(lbcnew)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  An interface to access bndry.inc commons from c routines.
!  Allows program to reset lbcopt to save memory in call
!  to initgrdvar when it is used for purposes other than
!  initializing the model.
!
!  AUTHOR:
!  Keith Brewster, CAPS
!  Jan 30, 2002
!
!-----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: lbcnew

  INCLUDE 'bndry.inc'

  lbcopt=lbcnew

END SUBROUTINE set_lbcopt

!
!##################################################################
!##################################################################
!######                                                      ######
!######               SUBROUTINE EXTENVPRF                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE extenvprf(nx,ny,nz,lvlprof,                                  &
                     x,y,zp,xs,ys,zps,                                  &
                     u,v,ptprt,pprt,qv,ptbar,pbar,usc,vsc,rfrct,        &
                     radarx,radary,radarz,radius,rngmax,                &
                     zsnd,ktsnd,usnd,vsnd,rfrctsnd,ienvstat)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Finds environmental profile around radar location from gridded data.
!
!  AUTHOR:
!  Keith Brewster, CAPS
!  Jan 30, 2002
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: nx,ny,nz
  INTEGER, INTENT(IN) :: lvlprof
  REAL, INTENT(IN)    :: x(nx)
  REAL, INTENT(IN)    :: y(ny)
  REAL, INTENT(IN)    :: zp(nx,ny,nz)
  REAL, INTENT(OUT)   :: xs(nx)
  REAL, INTENT(OUT)   :: ys(ny)
  REAL, INTENT(OUT)   :: zps(nx,ny,nz)
  REAL, INTENT(IN)    :: u(nx,ny,nz)
  REAL, INTENT(IN)    :: v(nx,ny,nz)
  REAL, INTENT(IN)    :: ptprt(nx,ny,nz)
  REAL, INTENT(IN)    :: pprt(nx,ny,nz)
  REAL, INTENT(IN)    :: qv(nx,ny,nz)
  REAL, INTENT(IN)    :: ptbar(nx,ny,nz)
  REAL, INTENT(IN)    :: pbar(nx,ny,nz)
  REAL, INTENT(OUT)   :: usc(nx,ny,nz)
  REAL, INTENT(OUT)   :: vsc(nx,ny,nz)
  REAL, INTENT(OUT)   :: rfrct(nx,ny,nz)
  REAL, INTENT(IN)    :: radarx
  REAL, INTENT(IN)    :: radary
  REAL, INTENT(IN)    :: radarz
  REAL, INTENT(IN)    :: radius
  REAL, INTENT(IN)    :: rngmax
  REAL, INTENT(IN)    :: zsnd(lvlprof)
  REAL, INTENT(OUT)   :: ktsnd(lvlprof)
  REAL, INTENT(OUT)   :: usnd(lvlprof)
  REAL, INTENT(OUT)   :: vsnd(lvlprof)
  REAL, INTENT(OUT)   :: rfrctsnd(lvlprof)
  INTEGER, INTENT(OUT):: ienvstat

  INTEGER, PARAMETER :: nptsmin = 3

  INTEGER :: i,j,k,kbot,ktop,kext
  INTEGER :: iradar,jradar
  INTEGER :: ibeg,iend,jbeg,jend
  INTEGER :: knt,knttot
  REAL :: dx,dy
  REAL :: rngavg,radius2,dist2,distmin,whigh,wlow,accept,rfrgrad

  INCLUDE 'mp.inc'

  INTEGER :: ioutdomain

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  ienvstat=0

  dx=x(2)-x(1)
  dy=y(2)-y(1)

  DO i=1,nx-1
    xs(i)=0.5*(x(i)+x(i+1))
  END DO
  xs(nx)=xs(nx-1)+dx

  DO j=1,ny-1
    ys(j)=0.5*(y(j)+y(j+1))
  END DO
  ys(ny)=ys(ny-1)+dy

  DO k=1,nz-1
    DO j=1,ny-1
      DO i=1,nx-1
        zps(i,j,k)=0.5*(zp(i,j,k)+zp(i,j,k+1))
      END DO
    END DO
  END DO

  DO j=1,ny-1
    DO i=1,nx-1
      zps(i,j,nz)=2.0*zp(i,j,nz-1)-zp(i,j,nz-2)
    END DO
  END DO

!-----------------------------------------------------------------------
!
!  Bring u and v to a common grid, the scalar points.
!
!-----------------------------------------------------------------------

  CALL avgx(u, 0, nx,ny,nz, 1,nx-1,1,ny-1,1,nz-1, usc)
  CALL avgy(v, 0, nx,ny,nz, 1,nx-1,1,ny-1,1,nz-1, vsc)
!
!-----------------------------------------------------------------------
!
!  Find the index of scalar grid location that is closest to the radar.
!  Set bounds of searching to be within a box of
!  radarx-radius < xs < radarx+radius
!  radary-radius < ys < radary+radius
!
!-----------------------------------------------------------------------
!
  rngavg=radius
  radius2=rngavg*rngavg
  iradar=1+nint((radarx-xs(1))/dx)
  jradar=1+nint((radary-ys(1))/dy)
  IF(myproc == 0) print *, ' iradar=',iradar,', jradar=',jradar

  ibeg=iradar-(int(radius/dx)+1)
  ibeg=min(max(ibeg,2),nx-2)

  iend=iradar+(int(radius/dx)+1)
  iend=min(max(iend,2),nx-2)

  jbeg=jradar-(int(radius/dy)+1)
  jbeg=min(max(jbeg,2),ny-2)

  jend=jradar+(int(radius/dy)+1)
  jend=min(max(jend,2),ny-2)

  IF(myproc == 0) THEN
    print *, ' ibeg = ',ibeg,', iend= ',iend
    print *, ' jbeg = ',jbeg,', jend= ',jend
  END IF
!
! Special check if radar location is outside subdomain
!
  ioutdomain = 0
  IF( iradar > (nx-2) .OR. iradar < 2 .OR.                              &
      jradar > (ny-2) .OR. jradar < 2 ) THEN
    ioutdomain = 1
  END IF
  CALL mptotali(ioutdomain)

  IF (ioutdomain == nprocs) THEN   ! It is outside of all subdomains
    knt=0
    distmin=1.0E32
    DO j=jbeg,jend
      DO i=ibeg,iend
        dist2 = (xs(i)-radarx)*(xs(i)-radarx)                           &
               +(ys(j)-radary)*(ys(j)-radary)
        IF(dist2 < radius2) knt=knt+1
        distmin=min(distmin,dist2)
      END DO
    END DO
    CALL mpminr(distmin)     ! minimum in global domain
    distmin=sqrt(distmin)
    CALL mptotali(knt)      ! all should count.
    IF( knt < nptsmin ) THEN
      IF( distmin < rngmax ) THEN
        IF (myproc == 0) THEN
        WRITE(6,'(a,f9.0,a)')                                           &
        ' Too few points within radar radius ',(0.001*radius),' km.'
        WRITE(6,'(a,f9.0,a)')                                           &
        ' Expanding averaging radius to rngmax ',(0.001*rngmax),' km.'
        END IF
        rngavg=rngmax
        radius2=rngavg*rngavg

        ibeg=iradar-(int(rngmax/dx)+1)
        ibeg=min(max(ibeg,2),nx-2)

        iend=iradar+(int(rngmax/dx)+1)
        iend=min(max(iend,2),nx-2)

        jbeg=jradar-(int(rngmax/dy)+1)
        jbeg=min(max(jbeg,2),ny-2)

        jend=jradar+(int(rngmax/dy)+1)
        jend=min(max(jend,2),ny-2)

      ELSE

        IF (myproc == 0) THEN
        WRITE(6,'(a,f12.1,a)')                                         &
        ' extenvprf: Minimum distance ',(0.001*distmin),' km.'
        WRITE(6,'(a,f9.0,a)')                                          &
        ' No points in the domain within radius ',(0.001*rngmax),' km.'
        WRITE(6,'(a)') ' Radar range outside of domain.  Error exit'
        END IF

        CALL flush(6)
        ienvstat=-1
        RETURN

      END IF

    END IF

  END IF
!
!-----------------------------------------------------------------------
!
!  Initialize sums to zero
!
!-----------------------------------------------------------------------
!
  DO k=1,lvlprof
    ktsnd(k)=0.
    usnd(k)=0.
    vsnd(k)=0.
    rfrctsnd(k)=0.
  END DO
!
  CALL refract(nx,ny,nz,ptprt,pprt,qv,ptbar,pbar,rfrct)
!
  knt=0
  knttot=0
  DO j=jbeg,jend
    DO i=ibeg,iend
!
!-----------------------------------------------------------------------
!
!  Is this point within the ARPS domain?
!  Since the ARPS grid is Cartesian, need only compare
!  the external grid coordinates to the x and y limits
!  of the ARPS grid.
!
!-----------------------------------------------------------------------
!
      knttot=knttot+1
      dist2 = (xs(i)-radarx)*(xs(i)-radarx)                             &
             +(ys(j)-radary)*(ys(j)-radary)
      IF(dist2 < radius2 ) THEN
        knt=knt+1
!
!-----------------------------------------------------------------------
!
!  Interpolate external data in vertical onto profile
!  arrays.
!
!-----------------------------------------------------------------------
!
        DO k=1,lvlprof
          IF(zps(i,j,2) <= zsnd(k)) THEN
            DO kext=3,nz-1
              IF(zps(i,j,kext) >= zsnd(k)) EXIT
            END DO
            IF(kext > nz-1) EXIT
            whigh=(zsnd(k)-zps(i,j,kext-1))/                            &
                  (zps(i,j,kext)-zps(i,j,kext-1))
            wlow=1.-whigh
            ktsnd(k)=ktsnd(k)+1.
            usnd(k)=usnd(k)+                                            &
                  whigh*usc(i,j,kext)+wlow*usc(i,j,kext-1)
            vsnd(k)=vsnd(k)+                                            &
                  whigh*vsc(i,j,kext)+wlow*vsc(i,j,kext-1)
            rfrctsnd(k)=rfrctsnd(k)+                                    &
                  whigh*rfrct(i,j,kext)+wlow*rfrct(i,j,kext-1)
          END IF
        END DO
      END IF
    END DO
  END DO

  CALL mptotali(knt)
  CALL mptotali(knttot)
  CALL mpsumr(ktsnd,lvlprof)
  CALL mpsumr(usnd,lvlprof)
  CALL mpsumr(vsnd,lvlprof)
  CALL mpsumr(rfrctsnd,lvlprof)

  IF (myproc == 0) WRITE(6,'(/a,i6,a,/a,i6,a/)')                        &
       '  extenvprf: found ',knt,' points within radius',               &
       '  of ',knttot,' checked.'

  IF(knt < nptsmin) THEN
    IF (myproc == 0) THEN
    WRITE(6,'(a,f9.0,a)')                                               &
      ' Too few domain points within radar radius ',(0.001*rngavg),' km.'
    WRITE(6,'(a)') ' Radar range outside of domain.  Error exit'
    END IF
    CALL flush(6)
    ienvstat=-2
    RETURN
  END IF

  accept=0.3*float(knt)
!
!-----------------------------------------------------------------------
!
!  Find lowest height with data
!
!-----------------------------------------------------------------------
!
  IF(myproc == 0) WRITE(6,'(a)') ' Finding range of mean profile data ...'
  DO k=1,lvlprof-1
    IF(ktsnd(k) > accept) EXIT
  END DO
  kbot=k

!
! Catastrophic failure.  Abort now.
!

  IF (kbot > lvlprof) THEN
    IF(myproc == 0) WRITE(6,'(a)') ' No acceptable data in EXTENVPRF.  Abort.'
    ienvstat=-1
    RETURN
  END IF

!
!-----------------------------------------------------------------------
!
!  Find highest height with data
!
!-----------------------------------------------------------------------
!
  DO k=lvlprof,2,-1
    IF(myproc == 0) WRITE(6,'(a,f10.2,a,f6.0,a,f10.0)') ' z = ',zsnd(k),&
                  ' knt = ',ktsnd(k),' accept = ',accept
    IF(ktsnd(k) > accept) EXIT
  END DO
  ktop=k
!
  IF(myproc == 0) WRITE(6,'(a,f10.2,a,f10.2,a)')                       &
               ' Height of data for wind profile spans from ',         &
                 zsnd(kbot),' to ',zsnd(ktop),' meters.'

  IF((ktop - kbot) < 2 ) THEN
    IF (myproc == 0) THEN
      WRITE(6,'(a,i6,a,i6)')                                            &
         ' Problem creating evironmental sounding, kbot= ',kbot,        &
         ' ktop=',ktop
      WRITE(6,'(a)') ' Radar range may be outside of domain.  Error exit'
    END IF
    ienvstat=-2
    RETURN
  END IF
!
!-----------------------------------------------------------------------
!
!  Divide through to find average
!
!-----------------------------------------------------------------------
!
  DO k=kbot,ktop
    usnd(k)=usnd(k)/ktsnd(k)
    vsnd(k)=vsnd(k)/ktsnd(k)
    rfrctsnd(k)=rfrctsnd(k)/ktsnd(k)
  END DO
!
!-----------------------------------------------------------------------
!
!  Set variables "below-ground"
!  Zero gradiant assumed in usnd, vsnd.
!  Constant gradiant assumed in rfrsnd.
!
!-----------------------------------------------------------------------
!
  rfrgrad=(rfrctsnd(kbot+1)-rfrctsnd(kbot))/(zsnd(kbot+1)-zsnd(kbot))
  DO k=kbot-1,1,-1
    usnd(k)=usnd(kbot)
    vsnd(k)=vsnd(kbot)
    rfrctsnd(k)=rfrctsnd(kbot)+(zsnd(k)-zsnd(kbot))*rfrgrad
  END DO
!
!-----------------------------------------------------------------------
!
!  Set variables "above-top"
!  Zero gradiant assumed.
!
!-----------------------------------------------------------------------
!
  rfrgrad=(rfrctsnd(ktop)-rfrctsnd(ktop-1))/(zsnd(ktop)-zsnd(ktop-1))
  DO k=ktop+1,lvlprof
    usnd(k)=usnd(ktop)
    vsnd(k)=vsnd(ktop)
    rfrctsnd(k)=rfrctsnd(ktop)+(zsnd(k)-zsnd(ktop))*rfrgrad
  END DO

  RETURN
END SUBROUTINE extenvprf
!
!##################################################################
!##################################################################
!######                                                      ######
!######               SUBROUTINE EXTENVPRF                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE extenvprf2(nx,ny,nz,lvlprof,                                 &
                     x,y,zp,xs,ys,zps,                                  &
                     u,v,pt,p,qv,usc,vsc,rfrct,                         &
                     radarx,radary,radarz,radius,rngmax,                &
                     zsnd,ktsnd,usnd,vsnd,rfrctsnd,ienvstat)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Finds environmental profile around radar location from gridded data.
!
!  AUTHOR:
!  Keith Brewster, CAPS
!  Jan 30, 2002
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: nx,ny,nz
  INTEGER, INTENT(IN) :: lvlprof
  REAL, INTENT(IN)    :: x(nx)
  REAL, INTENT(IN)    :: y(ny)
  REAL, INTENT(IN)    :: zp(nx,ny,nz)
  REAL, INTENT(OUT)   :: xs(nx)
  REAL, INTENT(OUT)   :: ys(ny)
  REAL, INTENT(OUT)   :: zps(nx,ny,nz)
  REAL, INTENT(IN)    :: u(nx,ny,nz)
  REAL, INTENT(IN)    :: v(nx,ny,nz)
  REAL, INTENT(IN)    :: pt(nx,ny,nz)
  REAL, INTENT(IN)    :: p(nx,ny,nz)
  REAL, INTENT(IN)    :: qv(nx,ny,nz)
  REAL, INTENT(OUT)   :: usc(nx,ny,nz)
  REAL, INTENT(OUT)   :: vsc(nx,ny,nz)
  REAL, INTENT(OUT)   :: rfrct(nx,ny,nz)
  REAL, INTENT(IN)    :: radarx
  REAL, INTENT(IN)    :: radary
  REAL, INTENT(IN)    :: radarz
  REAL, INTENT(IN)    :: radius
  REAL, INTENT(IN)    :: rngmax
  REAL, INTENT(IN)    :: zsnd(lvlprof)
  REAL, INTENT(OUT)   :: ktsnd(lvlprof)
  REAL, INTENT(OUT)   :: usnd(lvlprof)
  REAL, INTENT(OUT)   :: vsnd(lvlprof)
  REAL, INTENT(OUT)   :: rfrctsnd(lvlprof)
  INTEGER, INTENT(OUT):: ienvstat

  INTEGER, PARAMETER :: nptsmin = 3

  INTEGER :: i,j,k,kbot,ktop,kext
  INTEGER :: iradar,jradar
  INTEGER :: ibeg,iend,jbeg,jend
  INTEGER :: knt,knttot
  REAL :: dx,dy
  REAL :: rngavg,radius2,dist2,distmin,whigh,wlow,accept,rfrgrad

  INCLUDE 'mp.inc'

  INTEGER :: ioutdomain

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  ienvstat=0

  dx=x(2)-x(1)
  dy=y(2)-y(1)

  DO i=1,nx-1
    xs(i)=0.5*(x(i)+x(i+1))
  END DO
  xs(nx)=xs(nx-1)+dx

  DO j=1,ny-1
    ys(j)=0.5*(y(j)+y(j+1))
  END DO
  ys(ny)=ys(ny-1)+dy

  DO k=1,nz-1
    DO j=1,ny-1
      DO i=1,nx-1
        zps(i,j,k)=0.5*(zp(i,j,k)+zp(i,j,k+1))
      END DO
    END DO
  END DO

  DO j=1,ny-1
    DO i=1,nx-1
      zps(i,j,nz)=2.0*zp(i,j,nz-1)-zp(i,j,nz-2)
    END DO
  END DO

!-----------------------------------------------------------------------
!
!  Bring u and v to a common grid, the scalar points.
!
!-----------------------------------------------------------------------

  CALL avgx(u, 0, nx,ny,nz, 1,nx-1,1,ny-1,1,nz-1, usc)
  CALL avgy(v, 0, nx,ny,nz, 1,nx-1,1,ny-1,1,nz-1, vsc)
!
!-----------------------------------------------------------------------
!
!  Find the index of scalar grid location that is closest to the radar.
!  Set bounds of searching to be within a box of
!  radarx-radius < xs < radarx+radius
!  radary-radius < ys < radary+radius
!
!-----------------------------------------------------------------------
!
  rngavg=radius
  radius2=rngavg*rngavg
  iradar=1+nint((radarx-xs(1))/dx)
  jradar=1+nint((radary-ys(1))/dy)
  IF(myproc == 0) print *, ' iradar=',iradar,', jradar=',jradar

  ibeg=iradar-(int(radius/dx)+1)
  ibeg=min(max(ibeg,2),nx-2)

  iend=iradar+(int(radius/dx)+1)
  iend=min(max(iend,2),nx-2)

  jbeg=jradar-(int(radius/dy)+1)
  jbeg=min(max(jbeg,2),ny-2)

  jend=jradar+(int(radius/dy)+1)
  jend=min(max(jend,2),ny-2)

  IF(myproc == 0) print *, ' ibeg = ',ibeg,', iend= ',iend
  IF(myproc == 0) print *, ' jbeg = ',jbeg,', jend= ',jend
!
! Special check if radar location is outside subdomain
!
  ioutdomain = 0
  IF( iradar > (nx-2) .OR. iradar < 2 .OR.                              &
      jradar > (ny-2) .OR. jradar < 2 ) THEN
    ioutdomain = 1
  END IF
  CALL mptotali(ioutdomain)

  IF (ioutdomain == nprocs) THEN   ! It is outside of all subdomains
    knt=0
    distmin=1.0E32
    DO j=jbeg,jend
      DO i=ibeg,iend
        dist2 = (xs(i)-radarx)*(xs(i)-radarx)                           &
               +(ys(j)-radary)*(ys(j)-radary)
        IF(dist2 < radius2) knt=knt+1
        distmin=min(distmin,dist2)
      END DO
    END DO
    CALL mpminr(distmin)     ! minimum in global domain
    distmin=sqrt(distmin)
    CALL mptotali(knt)      ! all should count.
    IF( knt < nptsmin ) THEN
      IF( distmin < rngmax ) THEN
        IF (myproc == 0) THEN
        WRITE(6,'(a,f9.0,a)')                                           &
        ' Too few points within radar radius ',(0.001*radius),' km.'
        WRITE(6,'(a,f9.0,a)')                                           &
        ' Expanding averaging radius to rngmax ',(0.001*rngmax),' km.'
        END IF
        rngavg=rngmax
        radius2=rngavg*rngavg

        ibeg=iradar-(int(rngmax/dx)+1)
        ibeg=min(max(ibeg,2),nx-2)

        iend=iradar+(int(rngmax/dx)+1)
        iend=min(max(iend,2),nx-2)

        jbeg=jradar-(int(rngmax/dy)+1)
        jbeg=min(max(jbeg,2),ny-2)

        jend=jradar+(int(rngmax/dy)+1)
        jend=min(max(jend,2),ny-2)

      ELSE

        IF (myproc == 0) THEN
        WRITE(6,'(a,f12.1,a)')                                         &
        ' extenvprf: Minimum distance ',(0.001*distmin),' km.'
        WRITE(6,'(a,f9.0,a)')                                          &
        ' No points in the domain within radius ',(0.001*rngmax),' km.'
        WRITE(6,'(a)') ' Radar range outside of domain.  Error exit'
        END IF

        CALL flush(6)
        ienvstat=-1
        RETURN

      END IF

    END IF

  END IF
!
!-----------------------------------------------------------------------
!
!  Initialize sums to zero
!
!-----------------------------------------------------------------------
!
  DO k=1,lvlprof
    ktsnd(k)=0.
    usnd(k)=0.
    vsnd(k)=0.
    rfrctsnd(k)=0.
  END DO
!
  CALL refract2(nx,ny,nz,pt,p,qv,rfrct)
!
  knt=0
  knttot=0
  DO j=jbeg,jend
    DO i=ibeg,iend
!
!-----------------------------------------------------------------------
!
!  Is this point within the ARPS domain?
!  Since the ARPS grid is Cartesian, need only compare
!  the external grid coordinates to the x and y limits
!  of the ARPS grid.
!
!-----------------------------------------------------------------------
!
      knttot=knttot+1
      dist2 = (xs(i)-radarx)*(xs(i)-radarx)                             &
             +(ys(j)-radary)*(ys(j)-radary)
      IF(dist2 < radius2 ) THEN
        knt=knt+1
!
!-----------------------------------------------------------------------
!
!  Interpolate external data in vertical onto profile
!  arrays.
!
!-----------------------------------------------------------------------
!
        DO k=1,lvlprof
          IF(zps(i,j,2) <= zsnd(k)) THEN
            DO kext=3,nz-1
              IF(zps(i,j,kext) >= zsnd(k)) EXIT
            END DO
            IF(kext > nz-1) EXIT
            whigh=(zsnd(k)-zps(i,j,kext-1))/                            &
                  (zps(i,j,kext)-zps(i,j,kext-1))
            wlow=1.-whigh
            ktsnd(k)=ktsnd(k)+1.
            usnd(k)=usnd(k)+                                            &
                  whigh*usc(i,j,kext)+wlow*usc(i,j,kext-1)
            vsnd(k)=vsnd(k)+                                            &
                  whigh*vsc(i,j,kext)+wlow*vsc(i,j,kext-1)
            rfrctsnd(k)=rfrctsnd(k)+                                    &
                  whigh*rfrct(i,j,kext)+wlow*rfrct(i,j,kext-1)
          END IF
        END DO
      END IF
    END DO
  END DO

  CALL mptotali(knt)
  CALL mptotali(knttot)
  CALL mpsumr(ktsnd,lvlprof)
  CALL mpsumr(usnd,lvlprof)
  CALL mpsumr(vsnd,lvlprof)
  CALL mpsumr(rfrctsnd,lvlprof)

  IF (myproc == 0) WRITE(6,'(/a,i6,a,/a,i6,a/)')                        &
       '  extenvprf: found ',knt,' points within radius',               &
       '  of ',knttot,' checked.'

  IF(knt < nptsmin) THEN
    IF (myproc == 0) THEN
    WRITE(6,'(a,f9.0,a)')                                               &
      ' Too few domain points within radar radius ',(0.001*rngavg),' km.'
    WRITE(6,'(a)') ' Radar range outside of domain.  Error exit'
    END IF
    CALL flush(6)
    ienvstat=-2
    RETURN
  END IF

  accept=0.3*float(knt)
!
!-----------------------------------------------------------------------
!
!  Find lowest height with data
!
!-----------------------------------------------------------------------
!
  IF(myproc == 0) WRITE(6,'(a)') ' Finding range of mean profile data ...'
  DO k=1,lvlprof-1
    IF(ktsnd(k) > accept) EXIT
  END DO
  kbot=k

!
! Catastrophic failure.  Abort now.
!

  IF (kbot > lvlprof) THEN
    IF(myproc == 0) WRITE(6,'(a)') ' No acceptable data in EXTENVPRF.  Abort.'
    ienvstat=-1
    RETURN
  END IF

!
!-----------------------------------------------------------------------
!
!  Find highest height with data
!
!-----------------------------------------------------------------------
!
  DO k=lvlprof,2,-1
    IF(myproc == 0) WRITE(6,'(a,f10.2,a,f6.0,a,f10.0)') ' z = ',zsnd(k),&
                  ' knt = ',ktsnd(k),' accept = ',accept
    IF(ktsnd(k) > accept) EXIT
  END DO
  ktop=k
!
  IF(myproc == 0) WRITE(6,'(a,f10.2,a,f10.2,a)')                       &
               ' Height of data for wind profile spans from ',         &
                 zsnd(kbot),' to ',zsnd(ktop),' meters.'

  IF((ktop - kbot) < 2 ) THEN
    IF (myproc == 0) THEN
      WRITE(6,'(a,i6,a,i6)')                                            &
         ' Problem creating evironmental sounding, kbot= ',kbot,        &
         ' ktop=',ktop
      WRITE(6,'(a)') ' Radar range may be outside of domain.  Error exit'
    END IF
    ienvstat=-2
    RETURN
  END IF
!
!-----------------------------------------------------------------------
!
!  Divide through to find average
!
!-----------------------------------------------------------------------
!
  DO k=kbot,ktop
    usnd(k)=usnd(k)/ktsnd(k)
    vsnd(k)=vsnd(k)/ktsnd(k)
    rfrctsnd(k)=rfrctsnd(k)/ktsnd(k)
  END DO
!
!-----------------------------------------------------------------------
!
!  Set variables "below-ground"
!  Zero gradiant assumed in usnd, vsnd.
!  Constant gradiant assumed in rfrsnd.
!
!-----------------------------------------------------------------------
!
  rfrgrad=(rfrctsnd(kbot+1)-rfrctsnd(kbot))/(zsnd(kbot+1)-zsnd(kbot))
  DO k=kbot-1,1,-1
    usnd(k)=usnd(kbot)
    vsnd(k)=vsnd(kbot)
    rfrctsnd(k)=rfrctsnd(kbot)+(zsnd(k)-zsnd(kbot))*rfrgrad
  END DO
!
!-----------------------------------------------------------------------
!
!  Set variables "above-top"
!  Zero gradiant assumed.
!
!-----------------------------------------------------------------------
!
  rfrgrad=(rfrctsnd(ktop)-rfrctsnd(ktop-1))/(zsnd(ktop)-zsnd(ktop-1))
  DO k=ktop+1,lvlprof
    usnd(k)=usnd(ktop)
    vsnd(k)=vsnd(ktop)
    rfrctsnd(k)=rfrctsnd(ktop)+(zsnd(k)-zsnd(ktop))*rfrgrad
  END DO

  RETURN
END SUBROUTINE extenvprf2

!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE RDENVPRF                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
  SUBROUTINE rdenvprf(nzsnd,sndfile,                              &
                      radarid,radlat,radlon,radelv,               &
                      zsnd,usnd,vsnd,rfrsnd,istatus)
  IMPLICIT NONE
  INTEGER :: nzsnd
  CHARACTER(LEN=256) :: sndfile
  CHARACTER(LEN=5) :: radarid
  REAL :: radlat
  REAL :: radlon
  REAL :: radelv
  REAL :: zsnd(nzsnd)
  REAL :: usnd(nzsnd)
  REAL :: vsnd(nzsnd)
  REAL :: rfrsnd(nzsnd)
  INTEGER :: istatus

  INTEGER, PARAMETER :: iunit=61
  REAL, PARAMETER :: epslatlon=0.001
  REAL, PARAMETER :: epselv=1.0
  CHARACTER(LEN=8) :: radarin
  CHARACTER(LEN=1) :: dummy
  INTEGER :: istat,k
  REAL :: radlatin,radlonin,radelvin

  istatus = 0
  radarin = 'NULL'
  radlatin = -999.
  radlonin = -999.
  radelvin = -999.
  OPEN(iunit,file=TRIM(sndfile),form='formatted',   &
       status='old',iostat=istatus)
  IF(istatus == 0) THEN
    READ(iunit,'(7x,a8)',iostat=istat) radarin
    IF( radarin(1:4) == radarid(1:4)) THEN
      READ(iunit,'(6x,f10.4,8x,f10.4,8x,f7.1)',iostat=istat) &
           radlatin,radlonin,radelvin
      IF(abs(radlatin-radlat) < epslatlon .AND. &
         abs(radlonin-radlon) < epslatlon .AND. &
         abs(radelvin-radelv) < epselv) THEN
        READ(iunit,'(a1)',iostat=istat) dummy
        READ(iunit,'(a1)',iostat=istat) dummy
        DO k=1,nzsnd
          READ(iunit,*,iostat=istat) &
            zsnd(k),usnd(k),vsnd(k),rfrsnd(k)
          IF(istat /= 0) THEN
            WRITE(6,'(a,i4)') ' Error reading sounding level: ',k
            istatus=-3
            EXIT
          END IF
        END DO
        WRITE(6,'(a,i4,a)')  &
         ' RDENVPRF: Success reading ',nzsnd,' levels.'
      ELSE
        WRITE(6,'(2a)') ' Mismatch of radar location:'
        WRITE(6,'(a,2f10.4)') ' Latitudes:  ',radlatin,radlat
        WRITE(6,'(a,2f10.4)') ' Longitudes: ',radlonin,radlon
        WRITE(6,'(a,2f10.4)') ' Elevations: ',radelvin,radelv
        istatus=-2
      END IF
    ELSE
      WRITE(6,'(4a)') ' Mismatch of radar names: ',TRIM(radarin),' ',TRIM(radarid)
      istatus=-1
    END IF
  ELSE
    WRITE(6,'(a,a)') ' Error opening sounding file: ',TRIM(sndfile)
  END IF
  CLOSE(iunit)
  RETURN
  END SUBROUTINE rdenvprf

!##################################################################
!##################################################################
!######                                                      ######
!######                   SUBROUTINE UV2VR                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
  SUBROUTINE uv2vr(rlen,elev, azimu,            &
                 rlat_rad,rlon_rad,ralt_rad,ugrid,vgrid,vr)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!  Project u,v to radial direction
!
!-----------------------------------------------------------------------
!
!  INPUT:
!    rlen       distance to radar location
!    elev       elevation angle
!    azimu      azimuth angle
!    rlat_rad   latitude of radar
!    rlon_rad   latitude of radar
!    ralt_rad   altitude of radar
!    ugrid      u
!    vgrid      v
!
!  OUTPUT:
!    vr         radial wind
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
!  INCLUDE FILES
!
!-----------------------------------------------------------------------
!
  INCLUDE 'phycst.inc'
!
!-----------------------------------------------------------------------
!
!  Variable Declarations:
!
!-----------------------------------------------------------------------
!
  REAL :: rlen, elev, azimu
  REAL :: rlat_rad,rlon_rad,ralt_rad
  REAL :: ugrid,vgrid
  REAL :: vr
!
!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------
!
  REAL :: azim,dist,dz,eleva,range,dhdr,dsdr,ddrot
  REAL :: xrad, yrad, r_horiz, z_verti
  REAL :: xgrid,ygrid,zgrid
  REAL :: rlat_grid,rlon_grid
  REAL :: uazmrad,vazmrad
  REAL :: DPI
  PARAMETER ( DPI=3.1415926535898/180.0 )

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
!
!-----------------------------------------------------------------------
!
!  Get x and y locations of each radar ob
!
!-----------------------------------------------------------------------
!
  CALL lltoxy(1,1,rlat_rad,rlon_rad, xrad, yrad)

  r_horiz = rlen*COS(elev*DPI)
  z_verti = rlen*SIN(elev*DPI)

  xgrid   = r_horiz*SIN( azimu*DPI ) + xrad
  ygrid   = r_horiz*COS( azimu*DPI ) + yrad
  zgrid   = z_verti

  CALL  xytoll(1,1,xgrid, ygrid, rlat_grid,rlon_grid)
!
!-----------------------------------------------------------------------
!
!  Find heading and distance from radar along ground
!  Store correlation of azimuth with u and v drections
!  in uazmrad and vazmrad, respectively.
!
!-----------------------------------------------------------------------
!
  CALL disthead(rlat_grid,rlon_grid,                                &
                rlat_rad,rlon_rad,                                  &
                azim,dist)
!
  CALL ddrotuv(1,rlon_grid, azim,1.,ddrot,                          &
               uazmrad,vazmrad)
!
!-----------------------------------------------------------------------
!
!  Loop in height
!
!-----------------------------------------------------------------------
!
  dz= zgrid - ralt_rad

  CALL beamelv(dz,dist,eleva,range)
  CALL dhdrange(eleva,range,dhdr)
  dsdr=SQRT(AMAX1(0.,(1.-dhdr*dhdr)))
!
!-----------------------------------------------------------------------
!
!  Project u-v to radial direction to get radial velocity
!
!-----------------------------------------------------------------------
!
  vr=(uazmrad*ugrid + vazmrad*vgrid) * dsdr
!
  RETURN
  END SUBROUTINE uv2vr

!##################################################################
!##################################################################
!######                                                      ######
!######                  SUBROUTINE UV2VRN                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
  SUBROUTINE uv2vrn(nzsnd,zsnd,rfrsnd,rfropt,rlen,elev,azimu,           &
                 rlat_rad,rlon_rad,ralt_rad,ugrid,vgrid,vr)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!  Project u,v to radial direction
!
!-----------------------------------------------------------------------
!
!  INPUT:
!    rlen       distance to radar location
!    elev       elevation angle
!    azimu      azimuth angle
!    rlat_rad   latitude of radar
!    rlon_rad   latitude of radar
!    ralt_rad   altitude of radar
!    ugrid      u
!    vgrid      v
!
!  OUTPUT:
!    vr         radial wind
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
!  INCLUDE FILES
!
!-----------------------------------------------------------------------
!
  INCLUDE 'phycst.inc'
!
!-----------------------------------------------------------------------
!
!  Variable Declarations:
!
!-----------------------------------------------------------------------
!
  INTEGER, INTENT(IN) :: nzsnd
  REAL, INTENT(IN)    :: zsnd(nzsnd)
  REAL, INTENT(IN)    :: rfrsnd(nzsnd)
  INTEGER, INTENT(IN) :: rfropt
  REAL, INTENT(IN)    :: rlen, elev, azimu
  REAL, INTENT(IN)    :: rlat_rad,rlon_rad,ralt_rad
  REAL, INTENT(IN)    :: ugrid,vgrid
  REAL, INTENT(OUT)   :: vr
!
!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------
!
  REAL :: azim,dist,dz,eleva,range,dhdr,dsdr,ddrot
  REAL :: xrad, yrad, r_horiz, z_verti
  REAL :: xgrid,ygrid,zgrid
  REAL :: rlat_grid,rlon_grid
  REAL :: uazmrad,vazmrad
  REAL :: DPI
  PARAMETER ( DPI=3.1415926535898/180.0 )

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
!
!-----------------------------------------------------------------------
!
!  Get x and y locations of each radar ob
!
!-----------------------------------------------------------------------
!
  CALL lltoxy(1,1,rlat_rad,rlon_rad, xrad, yrad)

  r_horiz = rlen*COS(elev*DPI)
  z_verti = rlen*SIN(elev*DPI)

  xgrid   = r_horiz*SIN( azimu*DPI ) + xrad
  ygrid   = r_horiz*COS( azimu*DPI ) + yrad
  zgrid   = z_verti

  CALL  xytoll(1,1,xgrid, ygrid, rlat_grid,rlon_grid)
!
!-----------------------------------------------------------------------
!
!  Find heading and distance from radar along ground
!  Store correlation of azimuth with u and v drections
!  in uazmrad and vazmrad, respectively.
!
!-----------------------------------------------------------------------
!
  CALL disthead(rlat_grid,rlon_grid,                                &
                rlat_rad,rlon_rad,                                  &
                azim,dist)
!
  CALL ddrotuv(1,rlon_grid, azim,1.,ddrot,                          &
               uazmrad,vazmrad)
!
!-----------------------------------------------------------------------
!
!  Loop in height
!
!-----------------------------------------------------------------------
!
  dz= zgrid - ralt_rad

  CALL beamelvn(nzsnd,zsnd,rfrsnd,ralt_rad,rfropt,dz,dist,eleva,range)
  CALL dhdrangn(nzsnd,zsnd,rfrsnd,ralt_rad,rfropt,eleva,range,dhdr)
  dsdr=SQRT(AMAX1(0.,(1.-dhdr*dhdr)))
!
!-----------------------------------------------------------------------
!
!  Project u-v to radial direction to get radial velocity
!
!-----------------------------------------------------------------------
!
  vr=(uazmrad*ugrid + vazmrad*vgrid) * dsdr
!
  RETURN
  END SUBROUTINE uv2vrn
!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINE QUADUNF                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
  SUBROUTINE quadunf(maxgate,maxazim,maxelev,nzsnd,nsort,              &
                     velchek,vmedlim,dazlim,iorder,rfropt,             &
                     sortmin,dsort,                                    &
                     kntgate,kntazim,kntelev,                          &
                     radarx,radary,rdralt,dazim,                       &
                     rngmin,rngmax,zsnd,rfrsnd,                        &
                     rngvol,azmvol,elvvol,elvmnvol,vnyqvol,            &
                     kntbin,rxvol,ryvol,rzvol,velvol,                  &
                     istatus)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Uses a least squares local quadratic fit to refill any data region
!  rejected during unfolding process
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster, CAPS
!
!  MODIFICATION HISTORY:
!
!  14-Jan-2004  Keith Brewster
!               Modified logic to avoid possibility of overflow of varsort.
!               Updated azimuth searching to match design of remapvol.
!
!  08-Sep-2008  Keith Brewster
!               Replaced median computation with a binned median method.
!
!  30-Aug-2010  Keith Brewster
!               Elvmean is now pre-computed
!               Modified logic to retain original velocity to check here
!                 for a velocity folding issue.
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    maxgate   Maximum gates in a radial
!    maxazim   Maximum radials per tilt
!    maxelev   Maximum number of tilts
!    nzsnd     Number of levels in sounding arrays
!    nsort     Number of bins in sorting vector for computing median
!
!    velchek  Threshold for checking data, good vs. flagged
!    vmedlim  Threshold limit for median check
!    dazlim   Maximum value of azimuth difference (grid vs data) to accept
!             Generally should be 30 degrees or less for velocity, 360 for refl
!    iorder   Order of polynomial to fit (1: linear, 2: quadratic)
!    rfropt   Rafractivity option (1: std atmos lapse, 2: avg actual lapse rate)
!    kntgate  Number of gates
!    kntazim  Number of azimuths for each tilt
!    kntelev  Number of elevation angles (tilts) of data
!
!    radarx   x-location of radar
!    radary   y-location of radar
!    rdralt   Elevation (m MSL) of radar
!    dazim    Approximate range resolution of data (degrees of azimuth)
!    rngmin   Minimum range (m) of data to use
!            (10 000 m or more to eliminate near field ground targets).
!    rngmax   Maximum range (m) of data to use
!
!    zsnd     Heights of levels in refractivity profile
!    rfrsnd   Refractivity profile
!
!    rngvvol  Range to gate in velocity 3-D volume
!    azmvvol  Azimuth angle in velocity 3-D volume
!    elvvvol  Elevation angle in velocity 3-D volume
!    velvol   Radar data 3-D volume
!
!    rxvol    x coordinate of radar volume elements
!    ryvol    y coordinate of radar volume elements
!    rzvol    z coordinate of radar volume elements
!
!
!  OUTPUT:
!
!    kntbin   Temporary array used for computing the median
!    velvol   Radar data 3-D volume
!    istatus  Status indicator
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
  INTEGER, INTENT(IN) :: nzsnd
  INTEGER, INTENT(IN) :: nsort

  REAL, INTENT(IN)    :: velchek
  REAL, INTENT(IN)    :: vmedlim
  REAL, INTENT(IN)    :: dazlim
  INTEGER, INTENT(IN) :: iorder
  INTEGER, INTENT(IN) :: rfropt

  REAL, INTENT(IN)    :: sortmin
  REAL, INTENT(IN)    :: dsort

  INTEGER, INTENT(IN) :: kntgate(maxazim,maxelev)
  INTEGER, INTENT(IN) :: kntazim(maxelev)
  INTEGER, INTENT(IN) :: kntelev

  REAL, INTENT(IN)    :: radarx
  REAL, INTENT(IN)    :: radary
  REAL, INTENT(IN)    :: rdralt
  REAL, INTENT(IN)    :: dazim
  REAL, INTENT(IN)    :: rngmin
  REAL, INTENT(IN)    :: rngmax

  REAL, INTENT(IN)    :: zsnd(nzsnd)
  REAL, INTENT(IN)    :: rfrsnd(nzsnd)

  REAL, INTENT(IN)    :: rngvol(maxgate,maxelev)
  REAL, INTENT(IN)    :: azmvol(maxazim,maxelev)
  REAL, INTENT(IN)    :: elvvol(maxazim,maxelev)
  REAL, INTENT(IN)    :: elvmnvol(maxelev)
  REAL, INTENT(IN)    :: vnyqvol(maxazim,maxelev)

  INTEGER, INTENT(OUT)   :: kntbin(nsort)
  REAL, INTENT(OUT)      :: rxvol(maxgate,maxazim,maxelev)
  REAL, INTENT(OUT)      :: ryvol(maxgate,maxazim,maxelev)
  REAL, INTENT(OUT)      :: rzvol(maxgate,maxazim,maxelev)

  REAL, INTENT(INOUT)    :: velvol(maxgate,maxazim,maxelev)

  INTEGER, INTENT(OUT)   :: istatus
!
!-----------------------------------------------------------------------
!
! Misc. Local Variables
!
!-----------------------------------------------------------------------
!
  INTEGER, PARAMETER :: n = 6

  REAL, PARAMETER :: eps = 1.0E-25

  REAL :: avar(n,n)
  REAL :: rhsvar(n)
  REAL :: avel(n,n)
  REAL :: rhsvel(n)
  REAL :: sol(n)
  REAL :: work(n,n+1)
  REAL :: work1d(n+1)

  REAL :: array(4,4)
  REAL :: rhsv(4)
  REAL :: solv(4)

  INTEGER :: ii,jj,kk,i,j,k,knt,kinbox,kntall,kntdaz
  INTEGER :: kok,isort,jsort,mid
  INTEGER :: kbgn,kend
  INTEGER :: igate,jazim,kelev,jazmin,jmirror,jend,jknt
  INTEGER :: iigate,jjazim,kkelev
  INTEGER :: istatal,istatwrt
  INTEGER :: nbeam,nrang
  INTEGER :: kntchk,kntfold,kntflag
!
! jsrchmn is the minimum number of radials to search to find data
! that are within the distance limits of the grid center.
!
  INTEGER, PARAMETER :: jsrchmn = 2
  REAL, PARAMETER :: dxthr0 = 1000.    ! minimum dx processing radius

  REAL :: deg2rad,rad2deg
  REAL :: twonyq,inv2nyq
  REAL :: sortscale
  REAL :: delx,dely,delz,dazimr,azdiff
  REAL :: ddx,ddxy,ddx2,ddy,ddy2
  REAL :: cosaz,sinaz,sfcr,zagl
  REAL :: sum,sum2,sdev,thresh,slrange,elijk,azimijk,time
  REAL :: varmax,varmin,varavg,varmean,varmed
  REAL :: daz,sumdaz,azspc,dxthr
  REAL :: xs,ys,zps

  REAL :: origvel,tstdev,fitvel,tstvel,tstdif,vardif
!
!-----------------------------------------------------------------------
!
! Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'

!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  deg2rad = atan(1.)/45.
  rad2deg = 1./deg2rad
  dazimr = dazim*deg2rad
  sortscale=1./dsort
!
  time=0.
!
  DO kelev=1,kntelev
    DO jazim=1,kntazim(kelev)
      cosaz=cos(deg2rad*azmvol(jazim,kelev))
      sinaz=sin(deg2rad*azmvol(jazim,kelev))
      DO igate=1,kntgate(jazim,kelev)
        CALL beamhgtn(nzsnd,zsnd,rfrsnd,rdralt,rfropt,                      &
             elvvol(jazim,kelev),rngvol(igate,kelev),zagl,sfcr)
        rxvol(igate,jazim,kelev)=radarx+sinaz*sfcr
        ryvol(igate,jazim,kelev)=radary+cosaz*sfcr
        rzvol(igate,jazim,kelev)=rdralt+zagl
      END DO
    END DO
  END DO
!
  DO kkelev = 1,kntelev
    kntdaz=0
    sumdaz=0.
    DO jjazim = 2,kntazim(kkelev)
      daz=abs(azmvol(jjazim,kkelev)-azmvol(jjazim-1,kkelev))
      daz=min(daz,(360.-daz))
      kntdaz=kntdaz+1
      sumdaz=sumdaz+daz
    END DO
    IF( kntdaz > 0 ) THEN
      azspc=sumdaz/float(kntdaz)
      print *, ' Azimuth spacing: ',azspc,' degrees'
      azspc=deg2rad*azspc
    ELSE
      azspc=1.0
    END IF
    kntchk=0
    kntfold=0
    kntflag=0
    DO jjazim = 1,kntazim(kkelev)
      twonyq=2.0*vnyqvol(jjazim,kkelev)
      inv2nyq=1./twonyq
      DO iigate= 1,kntgate(jjazim,kkelev)
        IF(velvol(iigate,jjazim,kkelev) < -1500.0) THEN
          kntchk=kntchk+1
          origvel=velvol(iigate,jjazim,kkelev)+2000.

          xs = rxvol(iigate,jjazim,kkelev)
          ys = ryvol(iigate,jjazim,kkelev)
          zps = rzvol(iigate,jjazim,kkelev)

          kok=0
          sum=0.
          sum2=0.
          sdev=0.
          varavg=999999.
          kntbin=0
          delx=xs-radarx
          dely=ys-radary
          delz=zps-rdralt
!
          slrange=rngvol(iigate,kkelev)
          elijk=elvvol(jjazim,kkelev)
          azimijk=azmvol(jjazim,kkelev)

          dxthr=max(dxthr0,(2.1*azspc*rngvol(iigate,kkelev)))
          print *, ' dxthr =',(0.001*dxthr),' km'
!
          varmax=-999.
          varmin=999.
          DO jj=1,n
            DO ii=1,n
              avar(ii,jj)=0.
            END DO
          END DO
!
          DO ii=1,n
            rhsvar(ii)=0.
          END DO

          kbgn = kkelev
          DO k=kkelev-1,1,-1
            IF((elvmnvol(kkelev)-elvmnvol(k)) > 0.1) THEN
              kbgn=k
              EXIT
            END IF
          END DO
          !kbgn = MAX(k,1)    ! Temporary fix by Y. Wang. It is an algorithm bug because
                             ! it did not consider boundary values when kkelev = 1.
          kend=kkelev
          DO k=kkelev+1,kntelev
            IF((elvmnvol(k)-elvmnvol(kkelev)) > 0.1) THEN
              kend=k
              EXIT
            END IF
          END DO
          print *, ' Using levels',kbgn,' to ',kend,' at klevel:',kkelev
          print *, '    elvmnvol(kkelev):',elvmnvol(kkelev)
!
!-----------------------------------------------------------------------
!
!  First pass, find min,max,mean,median.
!
!-----------------------------------------------------------------------
!
          DO kk=kbgn,kend
!
!-----------------------------------------------------------------------
!
!  Find nearest azimuth at this level
!
!-----------------------------------------------------------------------
!
            IF(kk == kkelev) THEN
              jazmin=jjazim
            ELSE
              azdiff=181.
              jazmin=1
              DO jazim=1,kntazim(kk)
                daz=azmvol(jazim,kk)-azimijk
                IF(daz > 180.) daz=daz-360.
                IF(daz < -180.) daz=daz+360.
                daz=abs(daz)
                IF(daz < azdiff) THEN
                  azdiff=daz
                  jazmin=jazim
                END IF
              END DO
            END IF

            jmirror=jazmin+(kntazim(kk)/2)
            IF(jmirror > kntazim(kk)) jmirror=jmirror-kntazim(kk)
!
!-----------------------------------------------------------------------
!
!  Loop forward from jazmin
!
!-----------------------------------------------------------------------
!
            jend=kntazim(kk)
            IF(jmirror > jazmin) jend=jmirror-1
            jknt=0
            DO jazim=jazmin,jend
              kinbox=0
              jknt=jknt+1
              daz=azmvol(jazim,kk)-azimijk
              IF(daz > 180.) daz=daz-360.
              IF(daz < -180.) daz=daz+360.
              IF(abs(daz) > dazlim) EXIT
              DO igate=1,kntgate(jazim,kk)
                ddx=rxvol(igate,jazim,kk)-xs
                ddy=ryvol(igate,jazim,kk)-ys
!
                IF( rngvol(igate,kk) > rngmin .AND.                    &
                    rngvol(igate,kk) < rngmax .AND.                    &
                      abs(ddx) < dxthr .AND. abs(ddy) < dxthr ) THEN
                  kinbox=kinbox+1
                  IF(velvol(igate,jazim,kk) > velchek ) THEN
                    isort=1+NINT((velvol(igate,jazim,kk)-sortmin)*sortscale)
                    isort=max(min(isort,nsort),1)
                    kntbin(isort)=kntbin(isort)+1
                    sum=sum+velvol(igate,jazim,kk)
                    sum2=sum2+(velvol(igate,jazim,kk)*velvol(igate,jazim,kk))
                    kok=kok+1
                  END IF   ! data ok
                END IF  ! inside box
              END DO ! igate
              IF(kinbox == 0 .AND. jknt > jsrchmn) EXIT
            END DO ! jazim
!
!-----------------------------------------------------------------------
!
!  IF kinbox > 0 continue from jazim=1
!
!-----------------------------------------------------------------------
!
            IF((kinbox > 0 .OR. jknt <= jsrchmn) .and. jend==kntazim(kk)) THEN
              DO jazim=1,jmirror-1
                kinbox=0
                jknt=jknt+1
                daz=azmvol(jazim,kk)-azimijk
                IF(daz > 180.) daz=daz-360.
                IF(daz < -180.) daz=daz+360.
                IF(abs(daz) > dazlim) EXIT
                DO igate=1,kntgate(jazim,kk)
!
                  ddx=rxvol(igate,jazim,kk)-xs
                  ddy=ryvol(igate,jazim,kk)-ys

                  IF( rngvol(igate,kk) > rngmin .AND.                    &
                      rngvol(igate,kk) < rngmax .AND.                    &
                      abs(ddx) < dxthr .AND. abs(ddy) < dxthr ) THEN
                    kinbox=kinbox+1
                    IF(velvol(igate,jazim,kk) > velchek ) THEN
                      isort=1+NINT((velvol(igate,jazim,kk)-sortmin)*sortscale)
                      isort=max(min(isort,nsort),1)
                      kntbin(isort)=kntbin(isort)+1
                      sum=sum+velvol(igate,jazim,kk)
                      sum2=sum2+(velvol(igate,jazim,kk)*velvol(igate,jazim,kk))
                      kok=kok+1
                    END IF
                  END IF
                END DO
                IF(kinbox == 0 .AND. jknt > jsrchmn) EXIT
              END DO
            END IF
!
!-----------------------------------------------------------------------
!
! Loop backward from jazmin
!
!-----------------------------------------------------------------------
!
            jend= 1
            IF(jmirror < jazmin) jend=jmirror
            jknt=0
            DO jazim=jazmin-1,jend,-1
              kinbox=0
              jknt=jknt+1
              daz=azmvol(jazim,kk)-azimijk
              IF(daz > 180.) daz=daz-360.
              IF(daz < -180.) daz=daz+360.
              IF(abs(daz) > dazlim) EXIT
              DO igate=1,kntgate(jazim,kk)
!
                ddx=rxvol(igate,jazim,kk)-xs
                ddy=ryvol(igate,jazim,kk)-ys

                IF( rngvol(igate,kk) > rngmin .AND.                    &
                    rngvol(igate,kk) < rngmax .AND.                    &
                    abs(ddx) < dxthr .AND. abs(ddy) < dxthr ) THEN

                  kinbox=kinbox+1

                  IF(velvol(igate,jazim,kk) > velchek ) THEN
                    isort=1+NINT((velvol(igate,jazim,kk)-sortmin)*sortscale)
                    isort=max(min(isort,nsort),1)
                    kntbin(isort)=kntbin(isort)+1
                    sum=sum+velvol(igate,jazim,kk)
                    sum2=sum2+(velvol(igate,jazim,kk)*velvol(igate,jazim,kk))
                    kok=kok+1
                  END IF
                END IF
              END DO
              IF(kinbox == 0 .AND. jknt > jsrchmn) EXIT
            END DO
!
!-----------------------------------------------------------------------
!
! If not yet outside box, continue from last radial.
!
!-----------------------------------------------------------------------
!
!           IF(i==10 .and. j==10 .and. k==5) &
!               print *, ' OUT3 jazim,kinbox= ',jazim,kinbox
            IF((kinbox > 0 .OR. jknt <= jsrchmn) .and. jend==1 ) THEN
              DO jazim=kntazim(kk),jmirror,-1
                kinbox=0
                jknt=jknt+1
                daz=azmvol(jazim,kk)-azimijk
                IF(daz > 180.) daz=daz-360.
                IF(daz < -180.) daz=daz+360.
                IF(abs(daz) > dazlim) EXIT
                DO igate=1,kntgate(jazim,kk)
!
                  ddx=rxvol(igate,jazim,kk)-xs
                  ddy=ryvol(igate,jazim,kk)-ys
!
                  IF( rngvol(igate,kk) > rngmin .AND.                    &
                      rngvol(igate,kk) < rngmax .AND.                    &
                      abs(ddx) < dxthr .AND. abs(ddy) < dxthr ) THEN
                    kinbox=kinbox+1
                    IF(velvol(igate,jazim,kk) > velchek ) THEN
                      isort=1+NINT((velvol(igate,jazim,kk)-sortmin)*sortscale)
                      isort=max(min(isort,nsort),1)
                      kntbin(isort)=kntbin(isort)+1
                      sum=sum+velvol(igate,jazim,kk)
                      sum2=sum2+(velvol(igate,jazim,kk)*velvol(igate,jazim,kk))
                      kok=kok+1
                    END IF
                  END IF
                END DO  ! igate
                IF(kinbox == 0 .AND. jknt > jsrchmn) EXIT
              END DO ! jazim
            END IF
          END DO ! kk
!
          IF( kok == 0 ) THEN
            kntflag=kntflag+1
            write(6,'(a,3i7,a,f10.2)') ' No neighbors near',           &
                   iigate,jjazim,kkelev,' Setting missing. Orig vel=', &
                   origvel
            velvol(iigate,jjazim,kkelev)=-777.
            CYCLE
          ELSE
            varavg=sum/float(kok)
            mid=(kok/2)+1
            kntall=0
            DO isort=1,nsort-1
              kntall=kntall+kntbin(isort)
              IF(kntall >= mid) EXIT
            END DO
            varmed=sortmin+((isort-1)*dsort)
            IF ( kok > 1 ) THEN
              sdev=sqrt((sum2-(sum*sum/float(kok)))/float(kok-1))
              thresh=max((2.*sdev),vmedlim)
            ELSE
              thresh=vmedlim
            END IF
!           print *, ' Velocity difference threshold:',thresh
          END IF
!
!-----------------------------------------------------------------------
!
!  Process data for local quadratic fit
!
!-----------------------------------------------------------------------
!
          DO kk=kend,kend
!
!-----------------------------------------------------------------------
!
!  Find nearest azimuth at this level
!
!-----------------------------------------------------------------------
!
            IF(kk == kkelev) THEN
              jazmin=jjazim
            ELSE
              azdiff=181.
              jazmin=1
              DO jazim=1,kntazim(kk)
                daz=azmvol(jazim,kk)-azimijk
                IF(daz > 180.) daz=daz-360.
                IF(daz < -180.) daz=daz+360.
                daz=abs(daz)
                IF(daz < azdiff) THEN
                  azdiff=daz
                  jazmin=jazim
                END IF
              END DO
            END IF

            jmirror=jazmin+(kntazim(kk)/2)
            IF(jmirror > kntazim(kk)) jmirror=jmirror-kntazim(kk)
!
!-----------------------------------------------------------------------
!
!  Loop forward from jazmin
!
!-----------------------------------------------------------------------
!
            jend=kntazim(kk)
            IF(jmirror > jazmin) jend=jmirror-1
            jknt=0
            DO jazim=jazmin,jend
              kinbox=0
              jknt=jknt+1
              daz=azmvol(jazim,kk)-azimijk
              IF(daz > 180.) daz=daz-360.
              IF(daz < -180.) daz=daz+360.
              IF(abs(daz) > dazlim) EXIT
              DO igate=1,kntgate(jazim,kk)
!
                ddx=rxvol(igate,jazim,kk)-xs
                ddy=ryvol(igate,jazim,kk)-ys
!
                IF( rngvol(igate,kk) > rngmin .AND.                    &
                    rngvol(igate,kk) < rngmax .AND.                    &
                    abs(ddx) < dxthr .AND. abs(ddy) < dxthr ) THEN

                  kinbox=kinbox+1
                  ddxy=ddx*ddy
                  ddx2=ddx*ddx
                  ddy2=ddy*ddy

                  IF(velvol(igate,jazim,kk) > velchek .AND.  &
                     abs(velvol(igate,jazim,kk)-varmed) < thresh ) THEN
!
                    varmax=max(varmax,velvol(igate,jazim,kk))
                    varmin=min(varmin,velvol(igate,jazim,kk))
!
                    rhsvar(1)=rhsvar(1)+velvol(igate,jazim,kk)
                    rhsvar(2)=rhsvar(2)+velvol(igate,jazim,kk)*ddx
                    rhsvar(3)=rhsvar(3)+velvol(igate,jazim,kk)*ddy
                    rhsvar(4)=rhsvar(4)+velvol(igate,jazim,kk)*ddxy
                    rhsvar(5)=rhsvar(5)+velvol(igate,jazim,kk)*ddx2
                    rhsvar(6)=rhsvar(6)+velvol(igate,jazim,kk)*ddy2
!
                    avar(1,1)=avar(1,1)+1.
                    avar(1,2)=avar(1,2)+ddx
                    avar(1,3)=avar(1,3)+ddy
                    avar(1,4)=avar(1,4)+ddxy
                    avar(1,5)=avar(1,5)+ddx2
                    avar(1,6)=avar(1,6)+ddy2
!
                    avar(2,1)=avar(2,1)+ddx
                    avar(2,2)=avar(2,2)+ddx2
                    avar(2,3)=avar(2,3)+ddx*ddy
                    avar(2,4)=avar(2,4)+ddx*ddxy
                    avar(2,5)=avar(2,5)+ddx*ddx2
                    avar(2,6)=avar(2,6)+ddx*ddy2
!
                    avar(3,1)=avar(3,1)+ddy
                    avar(3,2)=avar(3,2)+ddy*ddx
                    avar(3,3)=avar(3,3)+ddy2
                    avar(3,4)=avar(3,4)+ddy*ddx2
                    avar(3,5)=avar(3,5)+ddy*ddx2
                    avar(3,6)=avar(3,6)+ddy*ddy2
!
                    avar(4,1)=avar(4,1)+ddxy
                    avar(4,2)=avar(4,2)+ddxy*ddx
                    avar(4,3)=avar(4,3)+ddxy*ddy
                    avar(4,4)=avar(4,4)+ddxy*ddxy
                    avar(4,5)=avar(4,5)+ddxy*ddx2
                    avar(4,6)=avar(4,6)+ddxy*ddy2
!
                    avar(5,1)=avar(5,1)+ddx2
                    avar(5,2)=avar(5,2)+ddx2*ddx
                    avar(5,3)=avar(5,3)+ddx2*ddy
                    avar(5,4)=avar(5,4)+ddx2*ddxy
                    avar(5,5)=avar(5,5)+ddx2*ddx2
                    avar(5,6)=avar(5,6)+ddx2*ddy2
!
                    avar(6,1)=avar(6,1)+ddy2
                    avar(6,2)=avar(6,2)+ddy2*ddx
                    avar(6,3)=avar(6,3)+ddy2*ddy
                    avar(6,4)=avar(6,4)+ddy2*ddxy
                    avar(6,5)=avar(6,5)+ddy2*ddx2
                    avar(6,6)=avar(6,6)+ddy2*ddy2

                  END IF
!
                END IF
              END DO  ! igate
              IF(kinbox == 0 .AND. jknt > jsrchmn) EXIT
            END DO ! jazim
!
!-----------------------------------------------------------------------
!
!  IF kinbox > 0 continue from jazim=1
!
!-----------------------------------------------------------------------
!
            IF((kinbox > 0 .OR. jknt <= jsrchmn) .and. jend==kntazim(kk)) THEN
              DO jazim=1,jmirror-1
                kinbox=0
                jknt=jknt+1
                daz=azmvol(jazim,kk)-azimijk
                IF(daz > 180.) daz=daz-360.
                IF(daz < -180.) daz=daz+360.
                IF(abs(daz) > dazlim) EXIT
                DO igate=1,kntgate(jazim,kk)
!
                  ddx=rxvol(igate,jazim,kk)-xs
                  ddy=ryvol(igate,jazim,kk)-ys
!
                  IF( rngvol(igate,kk) > rngmin .AND.                    &
                      rngvol(igate,kk) < rngmax .AND.                    &
                      abs(ddx) < dxthr .AND. abs(ddy) < dxthr ) THEN

                    kinbox=kinbox+1
                    ddxy=ddx*ddy
                    ddx2=ddx*ddx
                    ddy2=ddy*ddy

                    IF(velvol(igate,jazim,kk) > velchek .AND.             &
                       abs(velvol(igate,jazim,kk)-varmed) < thresh ) THEN
!
                      varmax=max(varmax,velvol(igate,jazim,kk))
                      varmin=min(varmin,velvol(igate,jazim,kk))
!
                      rhsvar(1)=rhsvar(1)+velvol(igate,jazim,kk)
                      rhsvar(2)=rhsvar(2)+velvol(igate,jazim,kk)*ddx
                      rhsvar(3)=rhsvar(3)+velvol(igate,jazim,kk)*ddy
                      rhsvar(4)=rhsvar(4)+velvol(igate,jazim,kk)*ddxy
                      rhsvar(5)=rhsvar(5)+velvol(igate,jazim,kk)*ddx2
                      rhsvar(6)=rhsvar(6)+velvol(igate,jazim,kk)*ddy2
!
                      avar(1,1)=avar(1,1)+1.
                      avar(1,2)=avar(1,2)+ddx
                      avar(1,3)=avar(1,3)+ddy
                      avar(1,4)=avar(1,4)+ddxy
                      avar(1,5)=avar(1,5)+ddx2
                      avar(1,6)=avar(1,6)+ddy2
!
                      avar(2,1)=avar(2,1)+ddx
                      avar(2,2)=avar(2,2)+ddx2
                      avar(2,3)=avar(2,3)+ddx*ddy
                      avar(2,4)=avar(2,4)+ddx*ddxy
                      avar(2,5)=avar(2,5)+ddx*ddx2
                      avar(2,6)=avar(2,6)+ddx*ddy2

                      avar(3,1)=avar(3,1)+ddy
                      avar(3,2)=avar(3,2)+ddy*ddx
                      avar(3,3)=avar(3,3)+ddy2
                      avar(3,4)=avar(3,4)+ddy*ddxy
                      avar(3,5)=avar(3,5)+ddy*ddx2
                      avar(3,6)=avar(3,6)+ddy*ddy2
!
                      avar(4,1)=avar(4,1)+ddxy
                      avar(4,2)=avar(4,2)+ddxy*ddx
                      avar(4,3)=avar(4,3)+ddxy*ddy
                      avar(4,4)=avar(4,4)+ddxy*ddxy
                      avar(4,5)=avar(4,5)+ddxy*ddx2
                      avar(4,6)=avar(4,6)+ddxy*ddy2

                      avar(5,1)=avar(5,1)+ddx2
                      avar(5,2)=avar(5,2)+ddx2*ddx
                      avar(5,3)=avar(5,3)+ddx2*ddy
                      avar(5,4)=avar(5,4)+ddx2*ddxy
                      avar(5,5)=avar(5,5)+ddx2*ddx2
                      avar(5,6)=avar(5,6)+ddx2*ddy2

                      avar(6,1)=avar(6,1)+ddy2
                      avar(6,2)=avar(6,2)+ddy2*ddx
                      avar(6,3)=avar(6,3)+ddy2*ddy
                      avar(6,4)=avar(6,4)+ddy2*ddxy
                      avar(6,5)=avar(6,5)+ddy2*ddx2
                      avar(6,6)=avar(6,6)+ddy2*ddy2
!
                    END IF

                  END IF
                END DO  ! igate
                IF(kinbox == 0 .AND. jknt > jsrchmn) EXIT
              END DO ! jazim
            END IF
!
!-----------------------------------------------------------------------
!
! Loop backward from jazmin
!
!-----------------------------------------------------------------------
!
            jend= 1
            IF(jmirror < jazmin) jend=jmirror
            jknt=0
            DO jazim=jazmin-1,jend,-1
              kinbox=0
              jknt=jknt+1
              daz=azmvol(jazim,kk)-azimijk
              IF(daz > 180.) daz=daz-360.
              IF(daz < -180.) daz=daz+360.
              IF(abs(daz) > dazlim) EXIT
              DO igate=1,kntgate(jazim,kk)
!
                ddx=rxvol(igate,jazim,kk)-xs
                ddy=ryvol(igate,jazim,kk)-ys
!
                IF( rngvol(igate,kk) > rngmin .AND.                    &
                    rngvol(igate,kk) < rngmax .AND.                    &
                    abs(ddx) < dxthr .AND. abs(ddy) < dxthr ) THEN

                  kinbox=kinbox+1
                  ddxy=ddx*ddy
                  ddx2=ddx*ddx
                  ddy2=ddy*ddy

                  IF(velvol(igate,jazim,kk) > velchek .AND.             &
                     abs(velvol(igate,jazim,kk)-varmed) < thresh ) THEN
!
                    varmax=max(varmax,velvol(igate,jazim,kk))
                    varmin=min(varmin,velvol(igate,jazim,kk))
!
                    rhsvar(1)=rhsvar(1)+velvol(igate,jazim,kk)
                    rhsvar(2)=rhsvar(2)+velvol(igate,jazim,kk)*ddx
                    rhsvar(3)=rhsvar(3)+velvol(igate,jazim,kk)*ddy
                    rhsvar(4)=rhsvar(4)+velvol(igate,jazim,kk)*ddxy
                    rhsvar(5)=rhsvar(5)+velvol(igate,jazim,kk)*ddx2
                    rhsvar(6)=rhsvar(6)+velvol(igate,jazim,kk)*ddy2
!
                    avar(1,1)=avar(1,1)+1.
                    avar(1,2)=avar(1,2)+ddx
                    avar(1,3)=avar(1,3)+ddy
                    avar(1,4)=avar(1,4)+ddxy
                    avar(1,5)=avar(1,5)+ddx2
                    avar(1,6)=avar(1,6)+ddy2
!
                    avar(2,1)=avar(2,1)+ddx
                    avar(2,2)=avar(2,2)+ddx2
                    avar(2,3)=avar(2,3)+ddx*ddy
                    avar(2,4)=avar(2,4)+ddx*ddxy
                    avar(2,5)=avar(2,5)+ddx*ddx2
                    avar(2,6)=avar(2,6)+ddx*ddy2

                    avar(3,1)=avar(3,1)+ddy
                    avar(3,2)=avar(3,2)+ddy*ddx
                    avar(3,3)=avar(3,3)+ddy2
                    avar(3,4)=avar(3,4)+ddy*ddxy
                    avar(3,5)=avar(3,5)+ddy*ddx2
                    avar(3,6)=avar(3,6)+ddy*ddy2

                    avar(4,1)=avar(4,1)+ddxy
                    avar(4,2)=avar(4,2)+ddxy*ddx
                    avar(4,3)=avar(4,3)+ddxy*ddy
                    avar(4,4)=avar(4,4)+ddxy*ddxy
                    avar(4,5)=avar(4,5)+ddxy*ddx2
                    avar(4,6)=avar(4,6)+ddxy*ddy2
!
                    avar(5,1)=avar(5,1)+ddx2
                    avar(5,2)=avar(5,2)+ddx2*ddx
                    avar(5,3)=avar(5,3)+ddx2*ddy
                    avar(5,4)=avar(5,4)+ddx2*ddxy
                    avar(5,5)=avar(5,5)+ddx2*ddx2
                    avar(5,6)=avar(5,6)+ddx2*ddy2
!
                    avar(6,1)=avar(6,1)+ddy2
                    avar(6,2)=avar(6,2)+ddy2*ddx
                    avar(6,3)=avar(6,3)+ddy2*ddy
                    avar(6,4)=avar(6,4)+ddy2*ddxy
                    avar(6,5)=avar(6,5)+ddy2*ddx2
                    avar(6,6)=avar(6,6)+ddy2*ddy2
!
                  END IF
!
                END IF
              END DO  ! igate
              IF(kinbox == 0 .AND. jknt > jsrchmn) EXIT
            END DO ! jazim
!
!-----------------------------------------------------------------------
!
! If not yet outside box, continue from last radial.
!
!-----------------------------------------------------------------------
!
            IF((kinbox > 0 .OR. jknt <= jsrchmn) .and. jend==1 ) THEN
              DO jazim=kntazim(kk),jmirror,-1
                kinbox=0
                jknt=jknt+1
                daz=azmvol(jazim,kk)-azimijk
                IF(daz > 180.) daz=daz-360.
                IF(daz < -180.) daz=daz+360.
                IF(abs(daz) > dazlim) EXIT
                DO igate=1,kntgate(jazim,kk)
!
                  ddx=rxvol(igate,jazim,kk)-xs
                  ddy=ryvol(igate,jazim,kk)-ys
!
                  IF( rngvol(igate,kk) > rngmin .AND.                    &
                      rngvol(igate,kk) < rngmax .AND.                    &
                      abs(ddx) < dxthr .AND. abs(ddy) < dxthr ) THEN

                    kinbox=kinbox+1
                    ddxy=ddx*ddy
                    ddx2=ddx*ddx
                    ddy2=ddy*ddy

                    IF(velvol(igate,jazim,kk) > velchek .AND.             &
                       abs(velvol(igate,jazim,kk)-varmed) < thresh ) THEN
!
                      varmax=max(varmax,velvol(igate,jazim,kk))
                      varmin=min(varmin,velvol(igate,jazim,kk))
!
                      rhsvar(1)=rhsvar(1)+velvol(igate,jazim,kk)
                      rhsvar(2)=rhsvar(2)+velvol(igate,jazim,kk)*ddx
                      rhsvar(3)=rhsvar(3)+velvol(igate,jazim,kk)*ddy
                      rhsvar(4)=rhsvar(4)+velvol(igate,jazim,kk)*ddxy
                      rhsvar(5)=rhsvar(5)+velvol(igate,jazim,kk)*ddx2
                      rhsvar(6)=rhsvar(6)+velvol(igate,jazim,kk)*ddy2
!
                      avar(1,1)=avar(1,1)+1.
                      avar(1,2)=avar(1,2)+ddx
                      avar(1,3)=avar(1,3)+ddy
                      avar(1,4)=avar(1,4)+ddxy
                      avar(1,5)=avar(1,5)+ddx2
                      avar(1,6)=avar(1,6)+ddy2
!
                      avar(2,1)=avar(2,1)+ddx
                      avar(2,2)=avar(2,2)+ddx2
                      avar(2,3)=avar(2,3)+ddx*ddy
                      avar(2,4)=avar(2,4)+ddx*ddxy
                      avar(2,5)=avar(2,5)+ddx*ddx2
                      avar(2,6)=avar(2,6)+ddx*ddy2
!
                      avar(3,1)=avar(3,1)+ddy
                      avar(3,2)=avar(3,2)+ddy*ddx
                      avar(3,3)=avar(3,3)+ddy2
                      avar(3,4)=avar(3,4)+ddy*ddxy
                      avar(3,5)=avar(3,5)+ddy*ddx2
                      avar(3,6)=avar(3,6)+ddy*ddy2
!
                      avar(4,1)=avar(4,1)+ddx2
                      avar(4,2)=avar(4,2)+ddx2*ddx
                      avar(4,3)=avar(4,3)+ddx2*ddy
                      avar(4,4)=avar(4,4)+ddx2*ddxy
                      avar(4,5)=avar(4,5)+ddx2*ddx2
                      avar(4,6)=avar(4,6)+ddx2*ddy2
!
                      avar(5,1)=avar(5,1)+ddx2
                      avar(5,2)=avar(5,2)+ddx2*ddx
                      avar(5,3)=avar(5,3)+ddx2*ddy
                      avar(5,4)=avar(5,4)+ddx2*ddxy
                      avar(5,5)=avar(5,5)+ddx2*ddx2
                      avar(5,6)=avar(5,6)+ddx2*ddy2
!
                      avar(6,1)=avar(6,1)+ddy2
                      avar(6,2)=avar(6,2)+ddy2*ddx
                      avar(6,3)=avar(6,3)+ddy2*ddy
                      avar(6,4)=avar(6,4)+ddy2*ddxy
                      avar(6,5)=avar(6,5)+ddy2*ddx2
                      avar(6,6)=avar(6,6)+ddy2*ddy2
!
                    END IF

                  END IF
                END DO  ! igate
                IF(kinbox == 0 .AND. jknt > jsrchmn) EXIT
              END DO ! jazim
            END IF
!
          END DO ! kk
!
!-----------------------------------------------------------------------
!
!   Solve for variable at grid point
!
!-----------------------------------------------------------------------
!

          knt=nint(avar(1,1))
          IF ( iorder > 1 .and. knt > 6 ) THEN
            varmean=rhsvar(1)/avar(1,1)
            CALL GJELIM(n,avar,rhsvar,sol,work,work1d,eps,istatus)
            fitvel=min(varmax,max(varmin,sol(1)))
            write(6,'(3f7.1,a,i8,4f7.1)') azmvol(jjazim,kkelev),       &
                (0.001*rngvol(iigate,kkelev)),elvmnvol(kkelev),        &
                ' Qf Analysis1:',knt,  &
                varmin,varmean,varmax,fitvel
          ELSE IF ( iorder > 0 .and. knt > 5 ) THEN
            DO jj=1,4
              DO ii=1,4
                array(ii,jj)=avar(ii,jj)
              END DO
            END DO
            DO ii=1,4
              rhsv(ii)=rhsvar(ii)
            END DO
            CALL GJELIM(4,array,rhsv,solv,work,work1d,eps,istatus)
            fitvel=min(varmax,max(varmin,solv(1)))
            write(6,'(3f7.1,a,i8,4f7.1)') azmvol(jjazim,kkelev),      &
                (0.001*rngvol(iigate,kkelev)),elvmnvol(kkelev),       &
                ' Qf Analysis2:',knt,  &
                varmin,varmean,varmax,fitvel
          ELSE IF ( knt > 0 ) THEN
            varmean=rhsvar(1)/avar(1,1)
            fitvel=varmean
            write(6,'(3f7.1,a,i8,4f7.1)') azmvol(jjazim,kkelev),      &
                (0.001*rngvol(iigate,kkelev)),elvmnvol(kkelev),       &
                ' Qf Analysis3:',knt,  &
                varmin,varmean,varmax,fitvel
          END IF
          tstdev=twonyq*NINT((fitvel-origvel)*inv2nyq)
          tstvel=origvel+tstdev
          vardif=abs(tstvel-fitvel)

          IF(abs(vardif) < thresh) THEN
            kntfold=kntfold+1
            velvol(iigate,jjazim,kkelev) = tstvel
            write(6,'(a,f10.1,a,f10.1,a,f10.1)')            &
               ' Qf Unfold: Meas=',origvel,' Qfitvel=',     &
                              fitvel,'  New=',tstvel
          ELSE
            kntflag=kntflag+1
            velvol(iigate,jjazim,kkelev) = -777.0
               write(6,'(a,f10.1,a,f10.1,a)')               &
               ' Qf Marked bad: Meas=',origvel,' Qfitvel=', &
                              fitvel,'  New=-777.0'
          END IF
        END IF
      END DO
    END DO

    WRITE(6,'(/a/,a/,3i12)')                                     &
          '       After Least-Squares Quad Fit Folding Check',   &
          '       Gates Checked   Unfolded   Flagged: ',         &
                  kntchk,kntfold,kntflag
  END DO  ! level loop

  RETURN
  END SUBROUTINE quadunf
!
!##################################################################
!##################################################################
!######                                                      ######
!######                  SUBROUTINE WRTVEL                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
  SUBROUTINE wrtvel(i_angle,i_tilt,varid,                         &
                    iyear,imon,iday,ihr,imin,isec,                &
                    gsp_vel,rfrst_vel,vnyquist,                   &
                    radarid,radar_lat,radar_lon,radar_elv,        &
                    maxgate,maxazim,ngate,nazim,                  &
                    azim,elev,rvel);
!-----------------------------------------------------------------------
!
!  PURPOSE:
!  Write out a tilt of radar data in radar coordinates for plotting
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    maxgate  Maximum number of gates in a radial
!    maxazim  Maximum number of radials in a tilt
!    ngate    Number of gates in radial
!    nazim    Number of radials
!    i_angle  elevation angle
!    i_tilt   tilt number
!    varid    variable id
!    iyear,imon,iday,ihr,imin,isec    Data time
!    gsp_vel  gate spacing
!    rfrst_vel distance of first gate to radar
!    azim     azimuthal angles for each radial
!    elev     elevations for each radial
!    rvel     Doppler radial velocity
!
!  OUTPUT:
!    None
!
!
!  WORK ARRAYS:
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
  INTEGER :: maxgate
  INTEGER :: maxazim
  INTEGER :: ngate
  INTEGER :: nazim

  INTEGER :: i, j
  INTEGER :: i_angle,i_tilt
  INTEGER :: iyear,imon,iday,ihr,imin,isec
  INTEGER :: gsp_vel,rfrst_vel
  REAL :: vnyquist

  CHARACTER(LEN=4) :: radarid

  REAL :: radar_lat,radar_lon,radar_elv

  REAL :: rvel(maxgate,maxazim)
  REAL :: elev(maxazim)
  REAL :: azim(maxazim)
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=6) :: varid

  CHARACTER (LEN=256) :: tmpdirnam
  CHARACTER (LEN=256) :: vfnam
  INTEGER :: nunit
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
   INCLUDE 'globcst.inc'

!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  CALL gtlfnkey( runname, lfnkey )

  tmpdirnam = dirname
  IF ( LEN(tmpdirnam) == 0 .OR. tmpdirnam == ' ' ) THEN
    tmpdirnam = '.'
  END IF

  WRITE (vfnam,'(3a,i4.4,2i2.2,a,2i2.2,2a,i6.6)') TRIM(tmpdirnam),     &
   '/',TRIM(radarid),iyear,imon,iday,'_',ihr,imin,                     &
   '_',TRIM(varid),i_angle

  CALL getunit( nunit)
  write(6,'(a,a)') ' wrtvel filename: ',vfnam

  open(nunit,file=TRIM(vfnam),form='unformatted',status='unknown')
  write(nunit) radarid
  write(nunit) iyear,imon,iday,ihr,imin,isec

  write(nunit) i_angle,i_tilt,rfrst_vel,gsp_vel
  write(nunit) ngate,nazim
  write(nunit) radar_lat,radar_lon,radar_elv

  write(nunit) vnyquist
  write(nunit) (elev(i),i=1,nazim)
  write(nunit) (azim(i),i=1,nazim)

  write(nunit) ((rvel(i,j),i=1,ngate),j=1,nazim)

  close(nunit)
  CALL retunit(nunit)

  END SUBROUTINE wrtvel
