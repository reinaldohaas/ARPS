!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE ENERGY                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE energy(nx,ny,nz,                                             &
           u,v,w,ptprt,pprt,qv,qscalar,rhobar,                          &
           x,y,z,zp,hterain, j3,                                        &
           rhostr,ustr,vstr,wstr,tem4)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Compute the density-weighted domain average of grid-scale kinetic
!  energy, momentum, potential temperature and potential temperature
!  variance.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  11/13/91.
!
!  MODIFICATION HISTORY:
!
!  6/06/92 (M. Xue)
!  Added full documentation.
!
!  This routine need to be checked for the terrain version.
!  5/6/93, Ming Xue
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!    u        x component of velocity (m/s)
!    v        y component of velocity (m/s)
!    w        Vertical component of Cartesian velocity (m/s)
!    ptprt    Perturbation potential temperature (K)
!    pprt     Perturbation pressure (Pascal)
!    qv       Water vapor specific humidity (kg/kg)
!    qc       Cloud water mixing ratio (kg/kg)
!    qr       Rainwater mixing ratio (kg/kg)
!    qi       Cloud ice mixing ratio (kg/kg)
!    qs       Snow mixing ratio (kg/kg)
!    qh       Hail mixing ratio (kg/kg)
!
!    rhobar   Base state density (kg/m**3)
!
!    x        x coordinate of grid points in physical/comp. space (m)
!    y        y coordinate of grid points in physical/comp. space (m)
!    z        z coordinate of grid points in computational space (m)
!    zp       Vertical coordinate of grid points in physical space (m)
!    hterain  Terrain height (m)
!
!    j3       Coordinate transformation Jacobian  d(zp)/dz
!
!  OUTPUT:
!
!    None.
!
!  WORK ARRAY:
!
!    rhostr   j3 times base state density rhobar, a work array.
!    ustr     rhobar*u, work array.
!    vstr     rhobar*v, work array.
!    wstr     rhobar*w, work array.
!    tem4     Temporary work array.
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
  INCLUDE 'phycst.inc'
  INCLUDE 'globcst.inc'
  INCLUDE 'mp.inc'            ! Message passing parameters.
  INCLUDE 'timelvls.inc'
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  INTEGER :: nx,ny,nz          ! Number of grid points in 3 directions

  REAL :: u     (nx,ny,nz,nt)  ! Total u-velocity (m/s)
  REAL :: v     (nx,ny,nz,nt)  ! Total v-velocity (m/s)
  REAL :: w     (nx,ny,nz,nt)  ! Total w-velocity (m/s)
  REAL :: ptprt (nx,ny,nz,nt)  ! Perturbation potential temperature (K)
  REAL :: pprt  (nx,ny,nz,nt)  ! Perturbation pressure (Pascal)

  REAL :: qv    (nx,ny,nz,nt)  ! Water vapor specific humidity (kg/kg)

  REAL :: qscalar(nx,ny,nz,nt,nscalar)

  REAL :: rhobar(nx,ny,nz)     ! Base state air density (kg/m**3)

  REAL :: x     (nx)           ! The x-coord. of the physical and
                               ! computational grid. Defined at u-point.
  REAL :: y     (ny)           ! The y-coord. of the physical and
                               ! computational grid. Defined at v-point.
  REAL :: z     (nz)           ! The z-coord. of the computational grid.
                               ! Defined at w-point on the staggered grid.
  REAL :: zp    (nx,ny,nz)     ! The physical height coordinate defined at
                               ! w-point of the staggered grid.

  REAL :: hterain(nx,ny)       ! Terrain height.

  REAL :: j3    (nx,ny,nz)     ! Coordinate transformation Jacobian defined
                               ! as d( zp )/d( z )

  REAL :: rhostr(nx,ny,nz)     ! Work array
  REAL :: ustr  (nx,ny,nz)     ! rhostr*u, work array
  REAL :: vstr  (nx,ny,nz)     ! rhostr*v, work array
  REAL :: wstr  (nx,ny,nz)     ! rhostr*w, work array
  REAL :: tem4  (nx,ny,nz)     ! Temporary work array
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: istat
  INTEGER :: tlevel            ! Time level at which data are printed.
  INTEGER :: i,j,k

  CHARACTER (LEN=256) :: engfn  ! File name of the energy statistics
                                     ! output
  INTEGER :: lengfn            ! String length of the file name

  REAL :: tmass                ! Total mass of the air in model domain (kg)
  REAL :: keu                  ! Contribution of u to the
                               ! total kinetic energy
  REAL :: kev                  ! Contribution of v to the
                               ! total kinetic energy
  REAL :: kew                  ! Contribution of w to the
                               ! total kinetic energy
  REAL :: ke                   ! Domain average kinetic energy per
                               ! unit mass (m/s)**2
  REAL :: ptvari               ! Domain average variance of ptprt
                               ! per unit mass (K**2)

  REAL :: momntu               ! Domain average x-component of momentum
                               ! per unit mass (m/s)
  REAL :: momntv               ! Domain average y-component of momentum
                               ! per unit mass (m/s)
  REAL :: momntw               ! Domain average z-component of momentum
                               ! per unit mass (m/s)
  REAL :: ptavg                ! Domain average potential temperature
                               ! perturbation (K)

  REAL :: dxdz05,dydz05,dxdy05

  INTEGER :: ncalls
  SAVE ncalls
  DATA ncalls /0/
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
!  Dump restart data files every nrstout number of time steps.
!
!-----------------------------------------------------------------------
!
  tlevel=tpresent
!
!-----------------------------------------------------------------------
!
!  Calculate ustr=rhostr*u, vstr=rhostr*v, wstr=rhostr*w
!
!-----------------------------------------------------------------------
!
  CALL aamult(rhobar,j3,nx,ny,nz,1,nx-1,1,ny-1,1,nz-1, rhostr)

  CALL rhouvw(nx,ny,nz,rhostr, ustr, vstr, wstr )

  DO k=1,nz-1
    DO j=1,ny-1
      DO i=1,nx
        ustr(i,j,k)=u(i,j,k,tlevel)*ustr(i,j,k)
      END DO
    END DO
  END DO

  DO k=1,nz-1
    DO j=1,ny
      DO i=1,nx-1
        vstr(i,j,k)=v(i,j,k,tlevel)*vstr(i,j,k)
      END DO
    END DO
  END DO

  DO k=1,nz
    DO j=1,ny-1
      DO i=1,nx-1
        wstr(i,j,k)=w(i,j,k,tlevel)*wstr(i,j,k)
      END DO
    END DO
  END DO
!
!-----------------------------------------------------------------------
!
!  Calculate the total mass of air inside the physical boundaries
!  based on the base state air density.
!
!-----------------------------------------------------------------------
!
  tmass = 0.0

  DO k=2,nz-2
    DO j=2,ny-2
      DO i=2,nx-2
        tmass = tmass + rhostr(i,j,k)                                   &
                      *(x(i+1)-x(i))*(y(j+1)-y(j))*(z(k+1)-z(k))
      END DO
    END DO
  END DO

  IF (mp_opt > 0) THEN
    CALL mptotal(tmass)
  END IF
!
!-----------------------------------------------------------------------
!
!  Determine the contribution of u to the total kinetic energy:
!
!-----------------------------------------------------------------------
!

  keu = 0.0

  DO k=2,nz-2
    DO j=2,ny-2

      dydz05 = (y(j+1)-y(j))*(z(k+1)-z(k))*0.5

      keu = keu + 0.5*ustr(  2, j,k)*u(  2, j,k,tlevel)                 &
                    *dydz05*(x(3)-x(1))

      keu = keu + 0.5*ustr(nx-1,j,k)*u(nx-1,j,k,tlevel)                 &
                    *dydz05*(x(nx)-x(nx-2))

      DO i=3,nx-2
        keu = keu + ustr(i,j,k)*u(i,j,k,tlevel)                         &
                    *dydz05*(x(i+1)-x(i-1))
      END DO

    END DO
  END DO

  IF (mp_opt > 0) THEN
    CALL mptotal(keu)
  END IF
!
!-----------------------------------------------------------------------
!
!  Determine the contribution of v to the total kinetic energy:
!
!-----------------------------------------------------------------------
!

  kev = 0.0

  DO k=2,nz-2
    DO i=2,nx-2

      dxdz05 = (x(i+1)-x(i))*(z(k+1)-z(k))*0.5

      kev = kev + 0.5*vstr(i, 2,  k)*v(i, 2,  k,tlevel)                 &
                    *dxdz05*(y(3)-y(1))

      kev = kev + 0.5*vstr(i,ny-1,k)*v(i,ny-1,k,tlevel)                 &
                    *dxdz05*(y(ny)-y(ny-2))

      DO j=3,ny-2
        kev = kev + vstr(i,j,k)*v(i,j,k,tlevel)                         &
                    *dxdz05*(y(j+1)-y(j-1))
      END DO

    END DO
  END DO

  IF (mp_opt > 0) THEN
    CALL mptotal(kev)
  END IF
!
!-----------------------------------------------------------------------
!
!  Determine the contribution of w to the total kinetic energy:
!
!-----------------------------------------------------------------------
!
  kew = 0.0

  DO i=2,nx-2
    DO j=2,ny-2

      dxdy05 = (x(i+1)-x(i))*(y(j+1)-y(j))*0.5

      kew = kew + 0.5*wstr(i,j, 2  )*w(i,j, 2  ,tlevel)                 &
                    *dxdy05*(z(3)-z(1))
      kew = kew + 0.5*wstr(i,j,nz-1)*w(i,j,nz-1,tlevel)                 &
                    *dxdy05*(z(nz)-z(nz-2))

      DO k=3,nz-2
        kew = kew + wstr(i,j,k)*w(i,j,k,tlevel)                         &
                    *dxdy05*(z(k+1)-z(k-1))
      END DO

    END DO
  END DO

  IF (mp_opt > 0) THEN
    CALL mptotal(kew)
  END IF
!
!-----------------------------------------------------------------------
!
!  The domain average kinetic energy per unit mass:
!
!-----------------------------------------------------------------------
!
  ke = 0.5*(keu+kev+kew)/tmass
!
!-----------------------------------------------------------------------
!
!  The domain average potential temperature variance:
!
!-----------------------------------------------------------------------
!
  ptvari = 0.0

  DO k=2,nz-2
    DO j=2,ny-2
      DO i=2,nx-2

        ptvari = ptvari + rhostr(i,j,k)*ptprt(i,j,k,tlevel)**2          &
                        *(x(i+1)-x(i))*(y(j+1)-y(j))*(z(k+1)-z(k))

      END DO
    END DO
  END DO

  IF (mp_opt > 0) THEN
    CALL mptotal(ptvari)
  END IF

  ptvari = ptvari/tmass

!
!-----------------------------------------------------------------------
!
!  Calculate the domain average momentum in the x direction.
!
!-----------------------------------------------------------------------
!

  momntu = 0.0

  DO k=2,nz-2
    DO j=2,ny-2

      dydz05 = (y(j+1)-y(j))*(z(k+1)-z(k))*0.5

      momntu = momntu + 0.5*ustr(  2, j,k)*dydz05*(x(3)-x(1))

      momntu = momntu + 0.5*ustr(nx-1,j,k)*dydz05*(x(nx)-x(nx-2))

      DO i=3,nx-2
        momntu = momntu + ustr(i,j,k)*dydz05*(x(i+1)-x(i-1))
      END DO

    END DO
  END DO

  IF (mp_opt > 0) THEN
    CALL mptotal(momntu)
  END IF

  momntu = momntu/tmass
!
!-----------------------------------------------------------------------
!
!  Calculate the domain average momentum in the y direction.
!
!-----------------------------------------------------------------------
!

  momntv = 0.0

  DO k=2,nz-2
    DO i=2,nx-2

      dxdz05 = (x(i+1)-x(i))*(z(k+1)-z(k))*0.5

      momntv = momntv + 0.5*vstr(i, 2,  k)*dxdz05*(y(3)-y(1))

      momntv = momntv + 0.5*vstr(i,ny-1,k)*dxdz05*(y(ny)-y(ny-2))

      DO j=3,ny-2
        momntv = momntv + vstr(i,j,k)*dxdz05*(y(j+1)-y(j-1))
      END DO

    END DO
  END DO

  IF (mp_opt > 0) THEN
    CALL mptotal(momntv)
  END IF

  momntv = momntv/tmass
!
!-----------------------------------------------------------------------
!
!  Calculate the domain average momentum in the z direction.
!
!-----------------------------------------------------------------------
!
  momntw = 0.0

  DO i=2,nx-2
    DO j=2,ny-2

      dxdy05 = (x(i+1)-x(i))*(y(j+1)-y(j))*0.5

      momntw = momntw + 0.5*wstr(i,j, 2  )*dxdy05*(z(3)-z(1))
      momntw = momntw + 0.5*wstr(i,j,nz-1)*dxdy05*(z(nz)-z(nz-2))

      DO k=3,nz-2
        momntw = momntw + wstr(i,j,k)*dxdy05*(z(k+1)-z(k-1))
      END DO

    END DO
  END DO

  IF (mp_opt > 0) THEN
    CALL mptotal(momntw)
  END IF

  momntw = momntw/tmass

!
!-----------------------------------------------------------------------
!
!  Calculate the mass-weighted domain average perturbation potential
!  temperature (K)
!
!-----------------------------------------------------------------------
!
  ptavg = 0.0

  DO k=2,nz-2
    DO j=2,ny-2
      DO i=2,nx-2

        ptavg = ptavg + rhostr(i,j,k)*ptprt(i,j,k,tlevel)               &
                        *(x(i+1)-x(i))*(y(j+1)-y(j))*(z(k+1)-z(k))

      END DO
    END DO
  END DO

  IF (mp_opt > 0) THEN
    CALL mptotal(ptavg)
  END IF

  ptavg = ptavg/tmass
!
!-----------------------------------------------------------------------
!
!  Write the results to a file
!
!-----------------------------------------------------------------------
!
  IF (myproc == 0) THEN

    IF(ncalls == 0) THEN

      IF( dirname /= ' ' ) THEN

        engfn  = dirname(1:ldirnam)//'/'//runname(1:lfnkey)//'.eng'
        lengfn = 6 + lfnkey + ldirnam + 1

      ELSE

        engfn  = runname(1:lfnkey)//'.eng'
        lengfn = 6 + lfnkey

      END IF

      CALL getunit( ncheng )

      OPEN(ncheng,FORM='formatted',STATUS='unknown',                    &
                FILE=engfn(1:lengfn),IOSTAT=istat)

      IF( istat /= 0) THEN

        WRITE(6,'(/a,i2)')                                              &
            ' Error occured when opening file '//runname(1:lfnkey)      &
            //'.eng'//                                                  &
            ' using FORTRAN unit ',ncheng
        CALL arpsstop('arpsstop called from energy',1)

      END IF

      WRITE(ncheng,'(a)') ''''//runname//''''

      WRITE(ncheng,'(t4,a,t15,a,t30,a,t45,a,t60,a,t75,a,t90,a)')        &
          '''TIME','    KE','  U-MOMENTUM','  V-MOMENTUM',              &
           '  w-momentum','   ptvari','     ptavg'''

    END IF

    WRITE(ncheng,'(f9.3,6f15.7)')                                       &
        curtim,ke,momntu,momntv,momntw,ptvari,ptavg

  END IF

  ncalls = 1

  RETURN
END SUBROUTINE energy
