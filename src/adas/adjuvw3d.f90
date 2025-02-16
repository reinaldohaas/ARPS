!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE ADJUVW                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE adjuvw( nx,ny,nz, u,v,w,wcont,ptprt,ptbar,                   &
                   zp,j1,j2,j3,aj3z,mapfct,rhostr,wcloud,               &
                   wndadj,obropt,obrzero,cldwopt,                       &
                   rhostru,rhostrv,rhostrw,rhs,divh,tem1,tem2,phi)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Diagnose w and wcont from u and v fields. Adjust u and v to ensure
!  anelastic mass conservation at each grid point.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  3/27/1995
!
!  MODIFICATION HISTORY:
!
!  4/19/1995 (M. Xue)
!  Redefined phi as a 2-D array, and added it into the argument list.
!
!  9/25/1996 (K. Brewster)
!  Added optional processing to set wcont=0 at bottom of
!  Rayleigh damping layer.
!
!  2/16/1998 (K. Brewster)
!  Added processing of wcloud, an estimate of w from the cloud
!  analysis.
!
!  7/27/1998 (M. Xue)
!  Map factor properly included. WCTOW and WCONTRA are now called
!  to convert between w and wcont. WCONTRA required that crdtrns
!  be defined, which is added to INITADAS.
!
!  2001-06-01 (Gene Bassett)
!  Fixed use of obropt > 9.  If subroutine was called more than once,
!  subsequent calls would use obropt-10 instead of obropt.
!
!  1/17/2005 (Kevin W. Thomas)
!  Add MPI support.
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!    u        x component of velocity at a given time level (m/s)
!    v        y component of velocity at a given time level (m/s)
!    w        z component of velocity at a given time level (m/s)
!    ptprt    perturbation potential temperature (K)
!    ptbar    base state potential temperature (K)
!    zp       model physical height
!    j1       Coordinate transformation Jacobian -d(zp)/dx
!    j2       Coordinate transformation Jacobian -d(zp)/dy
!    j3       Coordinate transformation Jacobian  d(zp)/dz
!    mapfct   map factor at scalar, u and v points.
!    rhostr   Base state density rhobar times j3 (kg/m**3)
!    wcloud   W from cloud analysis, if available
!
!    wndadj   Wind adjustment option.
!    obropt   O'Brien adjustment option.  Determines
!             distribution of mean divergence error used to
!             enforce upper boundary condition of w=0
!             = 1  Linearly in computational z. (default)
!             = 2  Linearly in physical z.
!             = 3  Linearly in potential temperature.
!    cldwopt  Option to replace w with w from cloud analysis
!
!  OUTPUT:
!
!    u        x component of velocity at a given time level (m/s)
!    v        y component of velocity at a given time level (m/s)
!    w        Vertical component of velocity in Cartesian
!             coordinates at a given time level (m/s)
!    wcont    Contravariant vertical velocity (m/s)
!
!  WORK ARRAYS:
!
!    rhostru  Temporary work array.
!    rhostrv  Temporary work array.
!    rhostrw  Temporary work array.
!    rhs      Temporary work array.
!    divh     Temporary work array storing horizontal mass divergence.
!    tem1     Temporary work array.
!    tem2     Temporary work array.
!    phi      Temporary work array.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx,ny,nz
  INTEGER :: nxlg,nylg

  INTEGER :: wndadj
  INTEGER :: obropt
  REAL    :: obrzero
  INTEGER :: cldwopt

  REAL :: u     (nx,ny,nz)     ! Total u-velocity (m/s)
  REAL :: v     (nx,ny,nz)     ! Total v-velocity (m/s)
  REAL :: w     (nx,ny,nz)     ! Total w-velocity (m/s)
  REAL :: ptprt (nx,ny,nz)     ! Perturbation potential temperature (K)
  REAL :: wcont (nx,ny,nz)     ! Contravariant vertical velocity (m/s)
  REAL :: ptbar (nx,ny,nz)     ! Base state potential temperature (K)

  REAL :: mapfct(nx,ny,8)
  REAL :: rhostr(nx,ny,nz)     ! Base state density rhobar times j3.
  REAL :: wcloud(nx,ny,nz)

  REAL :: zp(nx,ny,nz)

  REAL :: j1    (nx,ny,nz)     ! Coordinate transformation Jacobian defined as
                               ! - d( zp )/d( x ).
  REAL :: j2    (nx,ny,nz)     ! Coordinate transformation Jacobian defined as
                               ! - d( zp )/d( y ).
  REAL :: j3    (nx,ny,nz)     ! Coordinate transformation Jacobian defined as
                               ! d( zp )/d( z ).
  REAL :: aj3z  (nx,ny,nz)     ! Coordinate transformation Jacobian defined
                               ! as d( zp )/d( z ) AVERAGED IN THE Z-DIR.

  REAL :: rhostru(nx,ny,nz)    ! Work array
  REAL :: rhostrv(nx,ny,nz)    ! Work array
  REAL :: rhostrw(nx,ny,nz)    ! Work array
  REAL :: rhs    (nx,ny,nz)    ! Work array
  REAL :: divh   (nx,ny,nz)    ! Work array
  REAL :: tem1   (nx,ny,nz)    ! Work array
  REAL :: tem2   (nx,ny,nz)    ! Work array

  REAL :: phi  (0:nx,0:ny)     ! Work array
!
!-----------------------------------------------------------------------
!
!  Control parameters
!
!-----------------------------------------------------------------------
!
  LOGICAL :: chkmas
  PARAMETER (chkmas=.true.)
!
!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k,ktop
  INTEGER :: iterations,imax,imin,jmax,jmin,kmax,kmin
  REAL    :: dxx,dyy,absphi,delphi,rij,dphi
  REAL    :: divmax, divmin, wcmax, wcmin, dz2, dth2, zk, zk2, depth2
  REAL    :: pi, alpha, rdxx, rdyy, tem
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'grid.inc'        ! Grid parameters
  INCLUDE 'mp.inc'
  INCLUDE 'bndry.inc'
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  CALL rhouvw(nx,ny,nz,rhostr,rhostru,rhostrv,rhostrw)


  IF(wndadj == 1) THEN

    DO k=1,nz
      DO j=1,ny-1
        DO i=1,nx-1
          wcont(i,j,k)=0.
        END DO
      END DO
    END DO

    CALL wctow(nx,ny,nz,u,v,wcont,mapfct,                               &
               j1,j2,j3,aj3z,rhostr,rhostru,rhostrv,rhostrw,1,1,        &
               w,                                                       &
               tem1,tem2)

  ELSE IF (wndadj == 2 .OR. wndadj == 3 ) THEN

    w = 0.0
    wcont = 0.0
 
!-----------------------------------------------------------------------
!
!  Storing m**2 difx(ustr/m)+dify(vstr/m) in divh.
!
!-----------------------------------------------------------------------
!
    DO k=1,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          divh(i,j,k)=((u(i+1,j,k)*rhostru(i+1,j,k)*mapfct(i+1,j,5)     &
                       -u(i  ,j,k)*rhostru(i  ,j,k)*mapfct(i  ,j,5))    &
                       *dxinv                                           &
                      +(v(i,j+1,k)*rhostrv(i,j+1,k)*mapfct(i,j+1,6)     &
                       -v(i,j  ,k)*rhostrv(i,j  ,k)*mapfct(i,j  ,6))    &
                       *dyinv) *mapfct(i,j,1)*mapfct(i,j,1)
          tem1(i,j,k)=zp(i,j,k)
        END DO
      END DO
    END DO

    IF (mp_opt > 0) THEN
      CALL acct_interrupt(mp_acct)
      CALL mpsendrecv2dew(divh,nx,ny,nz,ebc,wbc,0,tem2)
      CALL mpsendrecv2dns(divh,nx,ny,nz,nbc,sbc,0,tem2)
    END IF

    DO k=3,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          tem2(i,j,k)=0.5*(ptprt(i,j,k)  +ptbar(i,j,k) +                &
                           ptprt(i,j,k-1)+ptbar(i,j,k-1))
        END DO
      END DO
    END DO

    DO j=1,ny-1
      DO i=1,nx-1
        wcont(i,j,2)=0.
      END DO
    END DO

    DO k=3,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          wcont(i,j,k)= wcont(i,j,k-1) - dz*divh(i,j,k-1)
        END DO
      END DO
    END DO
!
!-----------------------------------------------------------------------
!
!  wcont is now carrying wcont * rhostrw
!
!-----------------------------------------------------------------------
!

    IF (myproc == 0) WRITE(6,'(a,i4)') 'obropt = ',obropt

!
!-----------------------------------------------------------------------
!
!  For the MPI case, make sure each point gets used exactly once!
!
!-----------------------------------------------------------------------
!
    IF(obropt > 10) THEN
      ktop=0
      DO j=2,ny-1
        IF ((mp_opt > 0) .AND. (j == ny-1) .AND. (loc_y .NE. nproc_y)) EXIT
        DO i=2,nx-1
          IF ((mp_opt > 0) .AND. (i == nx-1) .AND. (loc_x .NE. nproc_x)) EXIT
          DO k=nz-2,3,-1
            IF(zp(i,j,k) < obrzero) EXIT
          END DO
          ktop=ktop+k
        END DO
      END DO

!
!-----------------------------------------------------------------------
!
!  For MPI, collect all the "ktop" computations.
!
!-----------------------------------------------------------------------
!

      IF (mp_opt > 0 ) THEN
        CALL mpsumi(ktop,1)
        nxlg = (nx - 3) * nproc_x + 3
        nylg = (ny - 3) * nproc_y + 3
        ktop=nint(FLOAT(ktop)/FLOAT((nxlg-2)*(nylg-2)))
        CALL mpupdatei(ktop,1)
      ELSE
        ktop=nint(FLOAT(ktop)/FLOAT((nx-2)*(ny-2)))
      END IF

      ktop=MIN((ktop+1),nz-1)

      IF (myproc == 0 ) WRITE(6,'(a,i5,a,f9.0,a)')                      &
      'Found mean k level at which w is zero for O''Brian adjustment:', &
      ktop,' obrzero=',obrzero,' m.'

      DO j=1,ny-1
        DO i=1,nx-1
          DO k=ktop+1,nz
            wcont(i,j,k)=0.
          END DO
        END DO
      END DO


    ELSE

      ktop=nz-1

    END IF

    DO j=1,ny-1
      DO i=1,nx-1
        tem1(i,j,1)=1./((tem1(i,j,ktop)-tem1(i,j,2))*                   &
                        (tem1(i,j,ktop)-tem1(i,j,2)))
        tem2(i,j,2)=1.5*tem2(i,j,3)-0.5*tem2(i,j,4)
        tem2(i,j,1)=1./((tem2(i,j,ktop)-tem2(i,j,2))*                   &
                      (tem2(i,j,ktop)-tem2(i,j,2)))
      END DO
    END DO

    IF (myproc == 0) WRITE(6,'(a)')                                     &
        ' Before w adjust, rho*wcont(top):'

    CALL a3dmax(wcont,1,nx,1,nx-1,                                      &
              1,ny,1,ny-1,1,nz,ktop,ktop,                               &
              wcmax,wcmin, imax,jmax,kmax, imin,jmin,kmin)

    IF (myproc == 0) WRITE(6,'(2(1x,a,g12.4,3(a,i3)))')                 &
        'wcmin =',wcmin,' at i=',imin,', j=',jmin,', k=',kmin,          &
        'wcmax =',wcmax,' at i=',imax,', j=',jmax,', k=',kmax

    IF( obropt == 1 .or. obropt == 11) THEN
!
!-----------------------------------------------------------------------
!
!  For obropt = 1,
!  we apply correction to wcont * rhostrw by assuming that
!  difz ( delta(wcont * rhostr) ) = A k, where A is
!  a constant. Therefore  delta( wcont*rhostr ) = A k**2/2.
!
!-----------------------------------------------------------------------
!
      IF (myproc == 0) WRITE(6,'(a)')                                   &
          'Horizontal divergence correction distributed linearly in k.'

      depth2 = 1./((ktop-2)*(ktop-2)*dz*dz)

      DO k=3,ktop
        zk2 = ((k-2)*dz)*((k-2)*dz)
        DO j=1,ny-1
          DO i=1,nx-1
            wcont(i,j,k)=wcont(i,j,k)-                                  &
                         wcont(i,j,ktop)*(zk2*depth2)
          END DO
        END DO
      END DO

    ELSE IF( obropt == 2 .or. obropt == 12) THEN
!
!-----------------------------------------------------------------------
!
!  For obropt = 2,
!  we apply correction to wcont * rhostrw by assuming that
!  difz ( delta(wcont * rhostr) ) = A zs, where A is
!  a constant. Therefore  delta( wcont*rhostr ) = A zp**2/2.
!
!-----------------------------------------------------------------------
!
      IF (myproc == 0) WRITE(6,'(a)')                                   &
          'Horizontal divergence correction distributed linearly in z.'

      DO k=3,ktop
        DO j=1,ny-1
          DO i=1,nx-1
            dz2=(tem1(i,j,k)-tem1(i,j,2))*(tem1(i,j,k)-tem1(i,j,2))
            wcont(i,j,k)=wcont(i,j,k)-                                  &
                         wcont(i,j,ktop)*(dz2*tem1(i,j,1))
          END DO
        END DO
      END DO

    ELSE IF( obropt == 3 .or. obropt == 13) THEN
!
!-----------------------------------------------------------------------
!
!  For obropt = 3,
!  We apply correction to wcont * rhostrw by assuming that
!  difz ( delta(wcont * rhostr) ) = A theta, where A is
!  a constant. Therefore  delta( wcont*rhostr ) = A theta**2/2.
!
!-----------------------------------------------------------------------
!
      IF (myproc == 0) WRITE(6,'(a)')                                   &
          'Horizontal div. correction distributed linearly in theta.'

      DO k=3,ktop
        DO j=1,ny-1
          DO i=1,nx-1
            dth2=(tem2(i,j,k)-tem2(i,j,2))*(tem2(i,j,k)-tem2(i,j,2))
            wcont(i,j,k)=wcont(i,j,k)-                                  &
                       wcont(i,j,ktop)*(dth2*tem2(i,j,1))
          END DO
        END DO
      END DO
    ELSE
      WRITE(6,'(a)')                                                    &
          'No OBrien vertical velocity/mean divergence adjustment'
    END IF

    IF (myproc == 0) WRITE(6,'(a)')                                     &
        ' After wcont adjustment, rho*wcont(top):'

    CALL a3dmax(wcont,1,nx,1,nx-1,                                      &
              1,ny,1,ny-1,1,nz,ktop,ktop,                               &
              wcmax,wcmin, imax,jmax,kmax, imin,jmin,kmin)

    IF (myproc == 0) WRITE(6,'(2(1x,a,g12.4,3(a,i3)))')                 &
        'wcmin =',wcmin,' at i=',imin,', j=',jmin,', k=',kmin,          &
        'wcmax =',wcmax,' at i=',imax,', j=',jmax,', k=',kmax
!
!-----------------------------------------------------------------------
!
!  Obtain wcont from wcont*rhostrw.
!
!-----------------------------------------------------------------------
!
    DO k=2,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          wcont(i,j,k)=wcont(i,j,k)/rhostrw(i,j,k)
        END DO
      END DO
    END DO
    DO j=1,ny-1
      DO i=1,nx-1
        wcont(i,j,1 )=-wcont(i,j,   3)
        wcont(i,j,nz)=-wcont(i,j,nz-2)
      END DO
    END DO

!
!-----------------------------------------------------------------------
!
!  Apply w from cloud analysis, if desired.
!  Reset wcont
!
!-----------------------------------------------------------------------
!
    IF(cldwopt > 0) THEN

      DO k=1,nz-1
        DO j=1,ny-1
          DO i=1,nx-1
            w(i,j,k)=AMAX1(w(i,j,k),wcloud(i,j,k))
          END DO
        END DO
      END DO

      CALL wcontra(nx,ny,nz,u,v,w,mapfct,j1,j2,j3,aj3z,                 &
           rhostr,rhostru,rhostrv,rhostrw,wcont,tem1,tem2)


    END IF

    IF(wndadj == 3) THEN

      IF( nx*ny*nz < (nx+1)*(ny+1) ) THEN
        PRINT*,'Work array PHI in ADJUVW too small. Job stopped.'
        CALL arpsstop("work array problem",1)
      END IF

      DO k=2,nz-2
        DO j=1,ny-1
          DO i=1,nx-1
            tem1(i,j,k) = divh(i,j,k) +                                 &
                       (wcont(i,j,k+1)*rhostrw(i,j,k+1)-                &
                        wcont(i,j,k  )*rhostrw(i,j,k  ))*dzinv
          END DO
        END DO
      END DO

      IF (myproc == 0) PRINT*,' Div max/min before u-v adjustment.'

      CALL a3dmax(tem1,1,nx,1,nx-1,                                     &
              1,ny,1,ny-1,1,nz,2,nz-2,                                  &
              divmax,divmin, imax,jmax,kmax, imin,jmin,kmin)

      IF (myproc == 0) WRITE(6,'(2(1x,a,g25.12,3(a,i3)))')              &
          'divmin =',divmin,' at i=',imin,', j=',jmin,', k=',kmin,      &
          'divmax =',divmax,' at i=',imax,', j=',jmax,', k=',kmax


      dxx = dx*dx
      dyy = dy*dy

      DO k=2,nz-2
        DO j=1,ny-1
          DO i=1,nx-1
            rhs(i,j,k) = - tem1(i,j,k)/(mapfct(i,j,1)**2)
          END DO
        END DO
      END DO

      pi = 4.0*ATAN(1.0)
      alpha = 4*(0.5 - pi/(2*SQRT(2.0))*                                &
                       SQRT(1.0/(nx+1)**2+1.0/(ny+1)**2) )

      rdxx = 1.0/dxx
      rdyy = 1.0/dyy
      tem = alpha*dxx*dyy/(2*(dxx+dyy))

      IF (myproc == 0) PRINT*,'Over-relaxation coeff = ', alpha

      DO j=0,ny
        DO i=0,nx
          phi(i,j) = 0.0
        END DO
      END DO

      DO k=2,nz-2

        zk = (k-1.5)*dz
        IF(k > 2) THEN
          DO j=1,ny-1
            DO i=1,nx-1
              phi(i,j) = phi(i,j)*zk/(zk-dz)
            END DO
          END DO
        END IF
!
!-----------------------------------------------------------------------
!
!  Solve Possion equation for the potential of the corrective
!  horizontal velocity. SOR method is used.
!  phi = zero is assumed at i=0 and nx, j=0 and ny.
!
!-----------------------------------------------------------------------
!
        DO iterations = 1,1000

          absphi = 0.0
          delphi = 0.0

          DO j=1,ny-1
            DO i=1,nx-1
              rij = (phi(i+1,j)-2*phi(i,j)+phi(i-1,j))*rdxx+            &
                    (phi(i,j+1)-2*phi(i,j)+phi(i,j-1))*rdyy-            &
                    rhs(i,j,k)
              dphi = tem*rij

              phi(i,j) = phi(i,j) + dphi
              absphi = ABS(phi(i,j))+absphi
              delphi = ABS(  dphi  )+delphi
            END DO
          END DO

          IF( delphi/(absphi+1.0E-30 ) < 1.0E-5 ) EXIT

        END DO

        IF (myproc == 0) PRINT*,'Iterations ended at', iterations,      &
               ' for k=,',k,' error=', delphi/(absphi+1.0E-30 )

!
!-----------------------------------------------------------------------
!
!  Apply the adjustment to u and v.
!
!-----------------------------------------------------------------------
!
        DO j=1,ny-1
          DO i=1,nx
            u(i,j,k)=u(i,j,k)+(phi(i,j)-phi(i-1,j))*mapfct(i,j,2)       &
                     /(rhostru(i,j,k)*dx)
          END DO
        END DO

        DO j=1,ny
          DO i=1,nx-1
            v(i,j,k)=v(i,j,k)+(phi(i,j)-phi(i,j-1))*mapfct(i,j,3)       &
                     /(rhostrv(i,j,k)*dy)
          END DO
        END DO

      END DO

      IF( chkmas ) THEN  ! Check mass continuity after adjustment

        DO k=2,nz-2
          DO j=1,ny-1
            DO i=1,nx-1
              tem1(i,j,k)=((u(i+1,j,k)*rhostru(i+1,j,k)*mapfct(i+1,j,5) &
                           -u(i  ,j,k)*rhostru(i  ,j,k)*mapfct(i  ,j,5)) &
                           *dxinv                                       &
                          +(v(i,j+1,k)*rhostrv(i,j+1,k)*mapfct(i,j+1,6) &
                           -v(i,j  ,k)*rhostrv(i,j  ,k)*mapfct(i,j  ,6)) &
                           *dyinv) *mapfct(i,j,1)*mapfct(i,j,1)         &
                           +(wcont(i,j,k+1)*rhostrw(i,j,k+1)            &
                            -wcont(i,j,k  )*rhostrw(i,j,k  ))*dzinv
            END DO
          END DO
        END DO

        IF (myproc == 0) PRINT*,' Div max/min after u-v adjustment.'

        CALL a3dmax(tem1,1,nx,1,nx-1,                                   &
                1,ny,1,ny-1,1,nz,2,nz-2,                                &
                divmax,divmin, imax,jmax,kmax, imin,jmin,kmin)

        IF (myproc == 0) WRITE(6,'(2(1x,a,g25.12,3(a,i3)))')            &
            'divmin =',divmin,' at i=',imin,', j=',jmin,', k=',kmin,    &
            'divmax =',divmax,' at i=',imax,', j=',jmax,', k=',kmax

      END IF
    END IF

    CALL wctow(nx,ny,nz,u,v,wcont,mapfct,                               &
               j1,j2,j3,aj3z,rhostr,rhostru,rhostrv,rhostrw,1,1,        &
               w,                                                       &
               tem1,tem2)

  END IF

  RETURN
END SUBROUTINE adjuvw
