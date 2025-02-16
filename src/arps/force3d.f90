!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE FRCUVW                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE frcuvw( nx,ny,nz,nzsoil,exbcbufsz,dtbig1,                    &
           u,v,w,wcont,ptprt,pprt,qv,qscalar,tke,pbldpth,               &
           ubar,vbar,ptbar,pbar,ptbari,pbari,rhostr,qvbar,              &
           usflx,vsflx, x,y,z,zp,zpsoil, mapfct,                        &
           j1,j2,j3,j3soil,aj3x,aj3y,aj3z,j3inv, sinlat,ptsfc,          &
           uforce,vforce,wforce,kmh,kmv,rprntl,lenscl,defsq,            &
           exbcbuf, bcrlx,rhofct,phydro,                                &
           tem1,tem2,                                                   &
           tem3,tem4,tem5,tem6,tem7,tem8,tem9,                          &
           tem10,tem11,tem12,tem13,tem14,tem15,tem16 )
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate the inactive acoustic forcing terms in the momentum
!  equations. These forcing terms include the advection, mixing, Coriolis
!  force and buoyancy terms, and are accumulated into one array for
!  each equation, i.e.,
!      uforce = - uadv + umix + ucorio.
!      vforce = - vadv + vmix + vcorio.
!      wforce = - wadv + wmix + wbuoy + wcorio.
!  These forcing terms are held fixed during the small acoustic wave
!  integration time steps. The turbulent mixing coefficient for
!  momentum km is an output which will be used to calculate the
!  turbulent mixing of heat and water substances.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  3/21/92.
!
!  MODIFICATION HISTORY:
!
!  5/05/92 (M. Xue)
!  Added full documentation.
!
!  4/10/93 (M. Xue & Hao Jin)
!  Add the terrain.
!
!  9/10/94 (D. Weber & Y. Lu)
!  Cleaned up documentation.
!
!  01/28/95 (G. Bassett)
!  Option added to turn off buoyancy terms (for buoyopt=0).
!
!  01/23/96 (Donghai Wang, Ming Xue and Yuhe Liu)
!  Added the map factor to forcing terms.
!
!  4/1/96 (Donghai Wang, X. Song and M. Xue)
!  Added the implicit treatment for the vertical mixing.
!
!  10/15/97 (Donghai Wang)
!  Using total density (rho) in the calculation of the pressure
!  gradient force terms (if rhofctopt=1).
!
!  11/05/97 (D. Weber)
!  Added phydro array for use in the bottom boundary condition for
!  pertrubation pressure (hydrostatic).
!
!  11/06/97 (D. Weber)
!  Added three additional levels to the mapfct array.  The three
!  levels (4,5,6) represent the inverse of the first three in order.
!  The inverse map factors are computed to improve efficiency.
!
!  11/06/97 (D. Weber)
!  Added tem12 for use in turbulent mixing (optimizing tmix3d.f).
!
!  9/18/98 (D. Weber)
!  Removed single operation do loops and placed the calculations
!  into the subroutines that generate the forcing terms.
!  Added arrays storing the average of j3 in the x,y, and z
!  directions.
!
!  1999/11/24 (Limin Zhao and Alan Shapiro)
!  Added the forward data assimilation package, which renews the
!  the wind, temperature and pressure fields before the forcing terms
!  are calculated when the assimilation mode is turned on
!  (assimopt=1).
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!    nzsoil   Number of soil levels
!  dtbig1     Local large time step size.
!
!    u        x component of velocity (m/s)
!    v        y component of velocity (m/s)
!    w        Vertical component of velocity in Cartesian
!             coordinates (m/s).
!    wcont    Contravariant vertical velocity (m/s)
!    ptprt    Perturbation potential temperature (K)
!    pprt     Perturbation pressure (Pascal)
!    qv       Water vapor specific humidity (kg/kg)
!
!    qscalar  Water/ice mixing ratios (kg/kg)
!
!    tke      Turbulent Kinetic Energy ((m/s)**2)
!    pbldpth  Planetary boundary layer depth (m)
!
!    ubar     Base state x velocity component (m/s)
!    vbar     Base state y velocity component (m/s)
!    ptbar    Base state potential temperature (K)
!    pbar     Base state pressure (Pascal)
!    ptbari   Inverse Base state potential temperature (K)
!    pbari    Inverse Base state pressure (Pascal)
!    rhostr   Base state density rhobar times j3 (kg/m**3)
!    qvbar    Base state water vapor specific humidity (kg/kg)
!
!    usflx    Surface flux of u-momentum
!    vsflx    Surface flux of v-momentum
!
!    x        x coordinate of grid points in physical/comp. space (m)
!    y        y coordinate of grid points in physical/comp. space (m)
!    z        z coordinate of grid points in computational space (m)
!    zp       Vertical coordinate of grid points in physical space(m)
!    zpsoil   Vertical coordinate of soil levels
!    mapfct   Map factors at scalar, u and v points
!
!    j1       Coordinate transformation Jacobian -d(zp)/dx
!    j2       Coordinate transformation Jacobian -d(zp)/dy
!    j3       Coordinate transformation Jacobian  d(zp)/dz
!    aj3x     Avgx of the coordinate transformation Jacobian d(zp)/dz
!    aj3y     Avgy of the coordinate transformation Jacobian d(zp)/dz
!    aj3z     Avgz of the coordinate transformation Jacobian d(zp)/dz
!    j3inv    Inverse of the coordinate transformation j3 d(zp)/dz
!
!    sinlat   Sin of latitude at each grid point
!    ptsfc    Potential temperature at the ground level (K)
!
!  OUTPUT:
!
!    uforce   Acoustically inactive forcing terms in u-momentum
!             equation (kg/(m*s)**2). uforce= -uadv + umix + ucorio
!    vforce   Acoustically inactive forcing terms in v-momentum
!             equation (kg/(m*s)**2). vforce= -vadv + vmix + vcorio
!    wforce   Acoustically inactive forcing terms in w-momentum
!             equation (kg/(m*s)**2).
!             wforce= -wadv + wmix + wbuoy + wcorio
!
!    phydro   Big time step forcing term for use in computing the
!             hydrostatic pressure at k=1.
!
!    kmh      Horizontal turb. mixing coef. for momentum ( m**2/s )
!    kmv      Vertical turb. mixing coef. for momentum ( m**2/s )
!    rprntl   Reciprocal of Prandtl number
!    lenscl   Turbulent mixing length scale (m)
!    defsq    Deformation squared (1/s**2)
!
!  WORK ARRAYS:
!
!    tem1     Temporary work array.
!    tem2     Temporary work array.
!    tem3     Temporary work array.
!    tem4     Temporary work array.
!    tem5     Temporary work array.
!    tem6     Temporary work array.
!    tem7     Temporary work array.
!    tem8     Temporary work array.
!    tem9     Temporary work array.
!    tem10    Temporary work array.
!    tem11    Temporary work array.
!    tem12    Temporary work array.
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
  INCLUDE 'globcst.inc'     ! Global constants that control model
                            ! execution
  INCLUDE 'bndry.inc'
  INCLUDE 'exbc.inc'
  INCLUDE 'assim.inc'
  INCLUDE 'timelvls.inc'

!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  INTEGER :: nx,ny,nz          ! Number of grid points in 3 directions
  INTEGER :: nzsoil            ! Number of soil levels

  REAL :: dtbig1               ! Local large time step size.

  REAL :: u     (nx,ny,nz,nt)  ! Total u-velocity (m/s)
  REAL :: v     (nx,ny,nz,nt)  ! Total v-velocity (m/s)
  REAL :: w     (nx,ny,nz,nt)  ! Total w-velocity (m/s)
  REAL :: wcont (nx,ny,nz)     ! Contravariant vertical velocity (m/s)
  REAL :: ptprt (nx,ny,nz,nt)  ! Perturbation potential temperature (K)
  REAL :: pprt  (nx,ny,nz,nt)  ! Perturbation pressure (Pascal)
  REAL :: qv    (nx,ny,nz,nt)  ! Water vapor specific humidity (kg/kg)

  REAL :: qscalar(nx,ny,nz,nt,nscalar)

  REAL :: tke   (nx,ny,nz,nt)  ! Turbulent kinetic energy ((m/s)**2)
  REAL :: pbldpth(nx,ny,nt)    ! Planetary boundary layer depth (m)

  REAL :: ubar  (nx,ny,nz)     ! Base state u-velocity (m/s)
  REAL :: vbar  (nx,ny,nz)     ! Base state v-velocity (m/s)
  REAL :: ptbar (nx,ny,nz)     ! Base state potential temperature (K)
  REAL :: pbar  (nx,ny,nz)     ! Base state pressure (Pascal).
  REAL :: ptbari(nx,ny,nz)     ! Inverse Base state pot. temperature (K)
  REAL :: pbari (nx,ny,nz)     ! Inverse Base state pressure (Pascal).
  REAL :: rhostr(nx,ny,nz)     ! Base state density rhobar times j3.
  REAL :: qvbar (nx,ny,nz)     ! Base state water vapor specific
                               ! humidity(kg/kg)

  REAL :: x     (nx)           ! x-coord. of the physical and compu-
                               ! tational grid. Defined at u-point(m).
  REAL :: y     (ny)           ! y-coord. of the physical and compu-
                               ! tational grid. Defined at v-point(m).
  REAL :: z     (nz)           ! z-coord. of the computational grid.
                               ! Defined at w-point on the staggered
                               ! grid(m).
  REAL :: zp    (nx,ny,nz)     ! Physical height coordinate defined at
                               ! w-point of the staggered grid(m).
  REAL :: zpsoil(nx,ny,nzsoil) ! Height of soil levels (m)
  REAL :: mapfct(nx,ny,8)      ! Map factors at scalar, u and v points

  REAL :: j1    (nx,ny,nz)     ! Coordinate transformation Jacobian
                               ! -d(zp)/d(x)
  REAL :: j2    (nx,ny,nz)     ! Coordinate transformation Jacobian
                               ! -d(zp)/d(y)
  REAL :: j3    (nx,ny,nz)     ! Coordinate transformation Jacobian
                               ! d(zp)/d(z)
  REAL :: j3soil(nx,ny,nzsoil) ! Coordinate transformation Jacobian
                               ! d(zpsoil)/d(z)
  REAL :: aj3x  (nx,ny,nz)     ! Coordinate transformation Jacobian defined
                               ! as d( zp )/d( z ) AVERAGED IN THE X-DIR.
  REAL :: aj3y  (nx,ny,nz)     ! Coordinate transformation Jacobian defined
                               ! as d( zp )/d( z ) AVERAGED IN THE Y-DIR.
  REAL :: aj3z  (nx,ny,nz)     ! Coordinate transformation Jacobian defined
                               ! as d( zp )/d( z ) AVERAGED IN THE Z-DIR.
  REAL :: j3inv (nx,ny,nz)     ! Coordinate transformation Jacobian
                               ! d(zp)/d(z)

  REAL :: sinlat(nx,ny)        ! Sin of latitude at each grid point
  REAL :: ptsfc (nx,ny)        ! Potential temperature at the ground level (K)

  REAL :: usflx (nx,ny)        ! Surface flux of u-momentum
                               ! (kg/(m*s**2))
  REAL :: vsflx (nx,ny)        ! Surface flux of v-momentum
                               ! (kg/(m*s**2))


  REAL :: uforce(nx,ny,nz)     ! Acoustically inactive forcing terms
                               ! in u-momentum equation (kg/(m*s)**2)
                               ! uforce= -uadv + umix + ucorio

  REAL :: vforce(nx,ny,nz)     ! Acoustically inactive forcing terms
                               ! in v-momentum equation (kg/(m*s)**2)
                               ! vforce= -vadv + vmix + vcorio

  REAL :: wforce(nx,ny,nz)     ! Acoustically inactive forcing terms
                               ! in w-momentum equation (kg/(m*s)**2)
                               ! wforce= -wadv + wmix + wbuoy

  REAL :: phydro(nx,ny)        ! Big time step forcing for computing
                               ! hydrostatic pprt at k=1.

  REAL :: kmh   (nx,ny,nz)     ! Horizontal turb. mixing coef. for
                               ! momentum. ( m**2/s )
  REAL :: kmv   (nx,ny,nz)     ! Vertical turb. mixing coef. for
                               ! momentum. ( m**2/s )
  REAL :: rprntl(nx,ny,nz)     ! Reciprocal of Prandtl number
  REAL :: lenscl(nx,ny,nz)     ! Turbulent mixing length scale (m)
  REAL :: defsq (nx,ny,nz)     ! Deformation squared (1/s**2)

  INTEGER :: exbcbufsz         ! EXBC buffer size
  REAL :: exbcbuf( exbcbufsz ) ! EXBC buffer array
  REAL :: bcrlx(nx,ny)         ! EXBC relaxation coefficients
  REAL :: rhofct(nx,ny,nz)     ! rho-factor:rhobar/rho

  REAL :: tem1  (nx,ny,nz)     ! Temporary work array.
  REAL :: tem2  (nx,ny,nz)     ! Temporary work array.
  REAL :: tem3  (nx,ny,nz)     ! Temporary work array.
  REAL :: tem4  (nx,ny,nz)     ! Temporary work array.
  REAL :: tem5  (nx,ny,nz)     ! Temporary work array.
  REAL :: tem6  (nx,ny,nz)     ! Temporary work array.
  REAL :: tem7  (nx,ny,nz)     ! Temporary work array.
  REAL :: tem8  (nx,ny,nz)     ! Temporary work array.
  REAL :: tem9  (nx,ny,nz)     ! Temporary work array.
  REAL :: tem10 (nx,ny,nz)     ! Temporary work array.
  REAL :: tem11 (nx,ny,nz)     ! Temporary work array.
  REAL :: tem12 (nx,ny,nz)     ! Temporary work array.
  REAL :: tem13 (nx,ny,nz)     ! Temporary work array.
  REAL :: tem14 (nx,ny,nz)     ! Temporary work array.
  REAL :: tem15 (nx,ny,nz)     ! Temporary work array.
  REAL :: tem16 (nx,ny,nz)     ! Temporary work array.
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k
  INTEGER :: tlevel
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  assimopt = 0     ! It is illegal to initilize variable in assim.inc
                   ! so move the initilization here by WYH.
!
!-----------------------------------------------------------------------
!
!  Calculate the total mixing (which includes subgrid scale
!  turbulent mixing and computational mixing) terms for u, v, w
!  equations as well as the mixing coefficient kmh and kmv, and
!  the inverse Prandtl number.
!
!  These mixing terms are accumulated in the arrays
!  uforce, vforce and wforce.
!
!  Since all mixing terms are evaluated at the past time level,
!  we pass only the variable fields at time tpast into the routine.
!
!-----------------------------------------------------------------------
!
  IF ( assimopt == 1 ) THEN
!
!-----------------------------------------------------------------------
!
!    ASSIMCON is a driver for the whole data assimilation system, which
!    determines the control parameters and switches in the assimilation
!    system.
!
!------------------------------------------------------------------------
!
    CALL assimcon(nx,ny,nz,x,y,z,zp,                                    &
                  ubar,vbar,pbar,ptbar,rhostr,qvbar,                    &
                  u,v,w,ptprt,pprt,qv,qscalar,                          &
                  j1,j2,j3,                                             &
                  tem1,tem2,tem3,tem4,tem5,tem6,tem7,tem8,tem9,         &
                  tem10,tem11)
!
!------------------------------------------------------------------------
!
!    ASSIMDRIV is a driver for the variational velocity adjustment,which
!    controls the ways to perform the varaitional adjustment and the
!    velocity blending.
!
!    Note: the arrays uforce, vforce, wforce, defsq are used in ASSIMDRIV
!       as working array. This will not effect any model results since
!       these array will be specified after the data assimilation.
!
!------------------------------------------------------------------------
!
    CALL assimdriv(nx,ny,nz,nzsoil,x,y,z,zp,zpsoil,                     &
                   u,v,w,wcont,ptprt,pprt,qv,qscalar,                   &
                   tke,kmh,kmv,mapfct,                                  &
                   ubar,vbar,ptbar,pbar,rhostr,qvbar,j1,j2,j3,j3soil,   &
                   tem1,tem2,tem3,tem4,tem5,tem6,tem7,tem8,tem9,        &
                   tem10,tem11,tem12,                                   &
                   tem13,tem14,tem15,tem16)
  END IF
!
!------------------------------------------------------------------------
!
!  If the assim mode is on, the modified model variables will be used
!  to calculate the force terms in the following subroutines.
!
!
!  Calculate the total mixing (which includes subgrid scale
!  turbulent mixing and computational mixing) terms for u, v, w
!  equations as well as the mixing coefficient kmh and kmv, and
!  the inverse Prandtl number.
!
!  These mixing terms are accumulated in the arrays
!  uforce, vforce and wforce.
!
!  Since all mixing terms are evaluated at the past time level,
!  we pass only the variable fields at time tpast into the routine.
!
!-----------------------------------------------------------------------
!

  IF(tintegopt == 1) THEN

    tlevel = tpast

  ELSE IF(tintegopt == 2 .or. tintegopt == 3) THEN

    tlevel = tpresent

  END IF
  CALL mixuvw(nx,ny,nz, exbcbufsz,                                      &
              u    (1,1,1,tlevel),v   (1,1,1,tlevel),                   &
              w    (1,1,1,tlevel),                                      &
              ptprt(1,1,1,tlevel),pprt(1,1,1,tlevel),                   &
              qv   (1,1,1,tlevel),qscalar(:,:,:,tlevel,:),              &
              tke, pbldpth(1,1,tpresent),                               &
              ubar,vbar,ptbar,                                          &
              pbar,rhostr,qvbar,                                        &
              usflx,vsflx, x,y,z,zp,mapfct,                             &
              j1,j2,j3,aj3x,aj3y,aj3z,j3inv, ptsfc,                     &
              uforce,vforce,wforce,kmh,kmv,rprntl,lenscl,defsq,         &
              exbcbuf,                                                  &
              tem1,tem2,tem3,tem4,tem5,tem6,tem7,tem8,tem9,tem10,       &
              tem11,tem12)

  IF( (tmixopt /= 0 .OR. cmix2nd /= 0 .OR. cmix4th /= 0 )               &
        .AND. lvldbg >= 4 ) THEN
    CALL set_acct(misc_acct)
    CALL checkuhx(uforce, nx,ny,nz,2,nx-1,1,ny-1,2,nz-2,                &
                  'umixx', tem1)
    CALL checkuhy(uforce, nx,ny,nz,2,nx-1,1,ny-1,2,nz-2,                &
                  'umixy', tem1)
    CALL checkvhx(vforce, nx,ny,nz,1,nx-1,2,ny-1,2,nz-2,                &
                  'vmixx', tem1)
    CALL checkvhy(vforce, nx,ny,nz,1,nx-1,2,ny-1,2,nz-2,                &
                  'vmixy', tem1)
    CALL checkwhx(wforce, nx,ny,nz,1,nx-1,1,ny-1,2,nz-1,                &
                  'wmixx', tem1)
    CALL checkwhy(wforce, nx,ny,nz,1,nx-1,1,ny-1,2,nz-1,                &
                  'wmixy', tem1)
  END IF
!
!
!-----------------------------------------------------------------------
!
!  Calculate the advection terms in the momentum equations, using
!  the equivalent advective formulation.
!
!  On exit of advuvw, the advection terms are included in the
!  forcing arrays which contain the turbulent mixing terms.
!
!-----------------------------------------------------------------------
!
!
  IF(tintegopt == 1 .or. tintegopt == 3) THEN

    CALL set_acct(advuvw_acct)

    CALL advuvw(nx,ny,nz,u,v,w,wcont,rhostr,ubar,vbar, mapfct,          &
                uforce,vforce,wforce,                                   &
                tem1,tem2,tem3,                                         &
                tem4,tem5,tem6,tem7,tem8,tem9)

    IF( lvldbg >= 4 ) THEN
      CALL checkuhx(uforce, nx,ny,nz,2,nx-1,1,ny-1,2,nz-2,              &
                    'uforce after advu', tem9)
      CALL checkuhy(uforce, nx,ny,nz,2,nx-1,1,ny-1,2,nz-2,              &
                    'uforce after advu', tem9)
      CALL checkvhx(vforce, nx,ny,nz,1,nx-1,2,ny-1,2,nz-2,              &
                    'vforce after advv', tem9)
      CALL checkvhy(vforce, nx,ny,nz,1,nx-1,2,ny-1,2,nz-2,              &
                    'vforce after advv', tem9)
      CALL checkwhx(wforce, nx,ny,nz,1,nx-1,1,ny-1,2,nz-1,              &
                    'wforce after advw', tem9)
      CALL checkwhy(wforce, nx,ny,nz,1,nx-1,1,ny-1,2,nz-1,              &
                    'wforce after advw', tem9)
    END IF

  END IF

!-----------------------------------------------------------------------
!
!  Calculate the advection terms from inhomogeneous map factor for
!  u and v, and store them into tem1 and tem2.
!
!  Ustr and vstr were calculated by advuvw and stored in tem4 and
!  tem5.
!
!-----------------------------------------------------------------------

  IF ( mptrmopt /= 0 ) THEN

    CALL maptrm( nx,ny,nz, u,v, tem4,tem5, mapfct,                      &
                 uforce, vforce, tem7,tem8,tem9 )

    IF ( lvldbg >= 5 ) THEN
      CALL checkuhx(uforce, nx,ny,nz,2,nx-1,1,ny-1,2,nz-2,              &
                  'uforce after maptrm', tem9)
      CALL checkuhy(uforce, nx,ny,nz,2,nx-1,1,ny-1,2,nz-2,              &
                  'uforce after maptrm', tem9)
      CALL checkvhx(vforce, nx,ny,nz,1,nx-1,2,ny-1,2,nz-2,              &
                  'vforce after maptrm', tem9)
      CALL checkvhy(vforce, nx,ny,nz,1,nx-1,2,ny-1,2,nz-2,              &
                  'vforce after maptrm', tem9)
    END IF

  END IF
!-----------------------------------------------------------------------
!
!  Calculate the Coriolis terms in u, v and w equations and add
!  to the forcing terms INSIDE the call to coriol.
!
!-----------------------------------------------------------------------

  IF( coriopt /= 0 ) THEN

    CALL set_acct(coriol_acct)

    tlevel = tpresent

    CALL coriol(nx,ny,nz,                                               &
                u(1,1,1,tlevel),v(1,1,1,tlevel),w(1,1,1,tlevel),        &
                ubar,vbar,rhostr, sinlat,mapfct,                        &
                uforce, vforce, wforce, tem4,tem5,tem6,tem7,tem8,tem9)

    IF( lvldbg >= 4 ) THEN
      CALL checkuhx(uforce, nx,ny,nz,2,nx-1,1,ny-1,2,nz-2,              &
                    'uforce after ucorx', tem9)
      CALL checkuhy(uforce, nx,ny,nz,2,nx-1,1,ny-1,2,nz-2,              &
                    'uforce after ucory', tem9)
      CALL checkvhx(vforce, nx,ny,nz,1,nx-1,2,ny-1,2,nz-2,              &
                    'vforce after vcorx', tem9)
      CALL checkvhy(vforce, nx,ny,nz,1,nx-1,2,ny-1,2,nz-2,              &
                    'vforce after vcory', tem9)
      CALL checkwhx(wforce, nx,ny,nz,1,nx-1,1,ny-1,2,nz-1,              &
                    'wforce after wcorx', tem9)
      CALL checkwhy(wforce, nx,ny,nz,1,nx-1,1,ny-1,2,nz-1,              &
                    'wforce after wcory', tem9)
    END IF

  END IF
!
!-----------------------------------------------------------------------
!
!    The temperature and pressure fields are modified by Gal-Chen's
!    thermodynamic recovery method if the assim mode is turned on.
!
!-----------------------------------------------------------------------
!
  IF ( assimopt == 1 ) THEN
    CALL set_acct(tinteg_acct)
    CALL assimptpr(nx,ny,nz,nzsoil,x,y,z,zp,zpsoil,                     &
                   u,v,w,ptprt,pprt,qv,qscalar,                         &
                   tke,kmh,kmv,                                         &
                   ubar,vbar,ptbar,pbar,qvbar,rhostr,                   &
                   uforce,vforce,wforce,j1,j2,j3,j3soil)

    WRITE(6,*)'back from assimptpr'
  END IF
!
!-----------------------------------------------------------------------
!
!  Calculate the rho-factor: rhobar/rho, to correct the calculation
!  for the pressure gradient terms
!
!-----------------------------------------------------------------------
!
  CALL set_acct(buoy_acct)

  IF (rhofctopt /= 0) THEN

    tlevel =  tpresent

    CALL rhofactor(nx,ny,nz,ptprt(1,1,1,tlevel),pprt(1,1,1,tlevel),     &
                 qv(1,1,1,tlevel),qscalar(:,:,:,tlevel,:),              &
                 ptbar,pbar,ptbari,pbari,qvbar,rhofct,tem1)

  ELSE

    DO k=1,nz
      DO j=1,ny
        DO i=1,nx
          rhofct(i,j,k) = 1.0
        END DO
      END DO
    END DO

  END IF
!
!-----------------------------------------------------------------------
!
!  Average rhofct at w points, stored in tem3
!
!-----------------------------------------------------------------------
!
  CALL avgsw(rhofct,nx,ny,nz, 1,nx-1, 1,ny-1,tem3)
!
!-----------------------------------------------------------------------
!
!  Calculate the buoyancy term for the w-equation.
!  The buoyancy term is stored in tem1 then added to wforce.
!
!  If buoyopt = 0 then buoyancy is turned off.
!
!-----------------------------------------------------------------------
!
  tlevel =  tpresent

  IF ( buoyopt /= 0 ) THEN

    CALL buoycy(nx,ny,nz, ptprt(1,1,1,tlevel),pprt(1,1,1,tlevel),       &
                qv(1,1,1,tlevel),qscalar(:,:,:,tlevel,:),               &
                ptbar,pbar,ptbari,pbari,rhostr,qvbar,                   &
                tem1,                                 & ! wbuoy
                tem2)

    DO k=2,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          wforce(i,j,k)=wforce(i,j,k)+tem3(i,j,k)*tem1(i,j,k)
        END DO
      END DO
    END DO

    DO j=1,ny-1  ! store buoyancy at k=2 into phydro.
      DO i=1,nx-1  ! for use in pprt(i,j,1) boundary condition
        phydro(i,j)=tem1(i,j,2)
      END DO
    END DO

  END IF

  IF( lvldbg >= 4 ) THEN
    CALL checkwhx(tem1, nx,ny,nz,1,nx-1,1,ny-1,2,nz-1,                  &
                  'wbuox', tem9)
    CALL checkwhy(tem1, nx,ny,nz,1,nx-1,1,ny-1,2,nz-1,                  &
                  'wbuoy', tem9)
  END IF
!
!-----------------------------------------------------------------------
!
!  To calculate additional boundary relaxation and mixing terms for
!  the wind fields in the case of externally forced boundary condition.
!
!-----------------------------------------------------------------------
!
  IF ( lbcopt == 2 .AND. mgrid == 1 ) THEN

    CALL set_acct(cmix_acct)

    CALL brlxuvw( nx,ny,nz, dtbig1,                                     &
                  u(1,1,1,1),v(1,1,1,1),w(1,1,1,1),rhostr,              &
                  uforce,vforce,wforce,                                 &
                  exbcbuf(nu0exb), exbcbuf(nv0exb),                     &
                  exbcbuf(nw0exb), exbcbuf(nudtexb),                    &
                  exbcbuf(nvdtexb),exbcbuf(nwdtexb),bcrlx,              &
                  tem1,tem2,tem3,tem4,tem5,tem6 )

  END IF

  IF ( lbcopt == 1 .AND. mgrid == 1 ) THEN
    CALL brlxuvw_rbc( nx,ny,nz, u,v,w,ubar,vbar,rhostr,               &
           uforce,vforce,wforce, bcrlx)
  END IF
!
!-----------------------------------------------------------------------
!
!  To calculate the vertically implicit mixing terms,
!  and add to uforce, vforce, and wforce.
!
!  This is only called here if tintegopt == 1 (the original leapfrog
!  time-stepping) or tintegopt == 3 (RK3 with all terms updated).
!  For tintegopt = 2 (RK3 case with only advection terms updated),
!  it will be calculated outside this subroutine.
!
!-----------------------------------------------------------------------
!
  IF (trbvimp == 1 .and. (tintegopt == 1 .or. tintegopt == 3)) THEN
    ! Vertical implicit application

    CALL set_acct(tmix_acct)

    CALL vmiximpuvw(nx,ny,nz,dtbig1,u(1,1,1,1),v(1,1,1,1),w(1,1,1,1),   &
                    rhostr,kmv,j1,j2,j3inv,                             &
                    uforce,vforce,wforce,                               &
                    tem1,tem2,tem3,tem4,                                &
                    tem5,tem6,tem7,tem8,tem9,tem10,tem11)

  END IF

  RETURN
END SUBROUTINE frcuvw
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE FRCP                       ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE frcp( nx,ny,nz,exbcbufsz, dtbig1,                            &
           u,v,w,wcont,ptprt,pprt,qv,qscalar,                           &
           ptbar,pbar,rhostr,qvbar,mapfct,j1,j2,j3,aj3x,aj3y,aj3z,      &
           pforce,                                                      &
           exbcbuf, bcrlx,                                              &
           padv,tem1,tem2,tem3,tem4,tem5,tem6,tem7,mp_tem)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate the total acoustically inactive forcing terms in
!  the pressure equation (advection and other source/sink terms).
!
!  These terms are invariant during the small time-step integration.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  3/21/92.
!
!  MODIFICATION HISTORY:
!
!  5/05/92 (M. Xue)
!  Added full documentation.
!
!  4/10/93 (M. Xue & Hao Jin)
!  Add the terrain.
!
!  9/10/94 (D. Weber & Y. Lu)
!  Cleaned up documentation.
!
!  01/23/96 (Donghai Wang, Yuhe Liu and Ming Xue)
!  Added the map factor to forcing terms.
!
!  9/18/98 (D. Weber)
!  Added arrays aj3x,aj3y,aj3z.
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical direction.
!    dtbig1   Local large time step size.
!
!    u        x component of velocity (m/s)
!    v        y component of velocity (m/s)
!    w        Vertical component of velocity in Cartesian
!             coordinates (m/s).
!    wcont    Contravariant vertical velocity (m/s)
!    ptprt    Perturbation potential temperature (K)
!    pprt     Perturbation pressure (Pascal)
!    qv       Water vapor specific humidity (kg/kg)
!    qc       Cloud water mixing ratio (kg/kg)
!    qr       Rainwater mixing ratio (kg/kg)
!    qi       Cloud ice mixing ratio (kg/kg)
!    qs       Snow mixing ratio (kg/kg)
!    qh       Hail mixing ratio (kg/kg)
!
!    ptbar    Base state potential temperature (K)
!    pbar     Base state pressure (Pascal)
!    rhostr   Base state density rhobar times j3 (kg/m**3)
!    qvbar    Base state water vapor specific humidity (kg/kg)
!
!    mapfct   Map factors at scalar points
!
!    j1       Coordinate transformation Jacobian -d(zp)/dx
!    j2       Coordinate transformation Jacobian -d(zp)/dy
!    j3       Coordinate transformation Jacobian  d(zp)/dz
!    aj3x     Avgx of the coordinate transformation Jacobian  d(zp)/dz
!    aj3y     Avgy of the coordinate transformation Jacobian  d(zp)/dz
!    aj3z     Avgz of the coordinate transformation Jacobian  d(zp)/dz
!
!  OUTPUT:
!
!    pforce   Acoustically inactive forcing terms in pressure
!             equation (Pascal/s). pforce= -padv
!
!  WORK ARRAYS:
!
!    padv     Advection of pressure
!    tem1     Temporary work array.
!    tem2     Temporary work array.
!    tem3     Temporary work array.
!    tem4     Temporary work array.
!    tem5     Temporary work array.
!    tem6     Temporary work array.
!    tem7     Temporary work array.
!    mp_tem   Temporary work array.
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
  INCLUDE 'globcst.inc'     ! Global constants that control model
                            ! execution
  INCLUDE 'bndry.inc'
  INCLUDE 'exbc.inc'
  INCLUDE 'timelvls.inc'

!
!-----------------------------------------------------------------------
!
!  Variable Declarations
!
!-----------------------------------------------------------------------
!
  INTEGER :: nx,ny,nz          ! Number of grid points in 3 directions
  REAL :: dtbig1               ! Local large time step size.

  REAL :: u     (nx,ny,nz,nt)  ! Total u-velocity (m/s)
  REAL :: v     (nx,ny,nz,nt)  ! Total v-velocity (m/s)
  REAL :: w     (nx,ny,nz,nt)  ! Total w-velocity (m/s)
  REAL :: wcont (nx,ny,nz)     ! Contravariant vertical velocity (m/s)
  REAL :: ptprt (nx,ny,nz,nt)  ! Perturbation potential temperature
                               ! from that of base state atmosphere
                               ! (Kelvin)
  REAL :: pprt  (nx,ny,nz,nt)  ! Perturbation pressure from that
                               ! of base state atmosphere (Pascal)
  REAL :: qv    (nx,ny,nz,nt)  ! Water vapor specific humidity (kg/kg)

  REAL :: qscalar(nx,ny,nz,nt,nscalar)

  REAL :: ptbar (nx,ny,nz)     ! Base state potential temperature (K)
  REAL :: pbar  (nx,ny,nz)     ! Base state pressure (Pascal).
  REAL :: rhostr(nx,ny,nz)     ! Base state density rhobar times j3.
  REAL :: qvbar (nx,ny,nz)     ! Base state water vapor specific
                               ! humidity(kg/kg)

  REAL :: mapfct(nx,ny)        ! Map factors at scalar points

  REAL :: j1    (nx,ny,nz)     ! Coordinate transformation Jacobian
                               ! -d(zp)/dx.
  REAL :: j2    (nx,ny,nz)     ! Coordinate transformation Jacobian
                               ! -d(zp)/dy.
  REAL :: j3    (nx,ny,nz)     ! Coordinate transformation Jacobian
                               ! d(zp)/dz.
  REAL :: aj3x  (nx,ny,nz)     ! Coordinate transformation Jacobian defined
                               ! as d( zp )/d( z ) AVERAGED IN THE X-DIR.
  REAL :: aj3y  (nx,ny,nz)     ! Coordinate transformation Jacobian defined
                               ! as d( zp )/d( z ) AVERAGED IN THE Y-DIR.
  REAL :: aj3z  (nx,ny,nz)     ! Coordinate transformation Jacobian defined
                               ! as d( zp )/d( z ) AVERAGED IN THE Z-DIR.

  REAL :: pforce(nx,ny,nz)     ! Sound-wave independent forcing terms
                               ! in pressure equation (Pascal/s).
                               ! pforce = -padv

  REAL :: padv  (nx,ny,nz)     ! The advection term of pressure eq.

  INTEGER :: exbcbufsz         ! EXBC buffer size
  REAL :: exbcbuf( exbcbufsz ) ! EXBC buffer array
  REAL :: bcrlx(nx,ny)         ! EXBC relaxation coefficients

  REAL :: tem1  (nx,ny,nz)     ! Temporary work array.
  REAL :: tem2  (nx,ny,nz)     ! Temporary work array.
  REAL :: tem3  (nx,ny,nz)     ! Temporary work array.
  REAL :: tem4  (nx,ny,nz)     ! Temporary work array.
  REAL :: tem5  (nx,ny,nz)     ! Temporary work array.
  REAL :: tem6  (nx,ny,nz)     ! Temporary work array.
  REAL :: tem7  (nx,ny,nz)     ! Temporary work array.

  REAL :: mp_tem(MAX(nx+1,ny+1)*(nz+1))  ! Temporary message passing array.
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k             ! local varaibles.
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!
!-----------------------------------------------------------------------
!
!  Calculate the advection term in the pressure equation
!  using the equivalent advective formulation
!
!-----------------------------------------------------------------------
!
!
  CALL set_acct(advs_acct)
  IF( peqopt == 1 ) THEN

    CALL advp(nx,ny,nz,pprt,u,v,w,wcont,rhostr,mapfct,                  &
              j3,aj3x,aj3y,aj3z,                                        &
              padv, tem1,tem2,tem3,tem4,tem5,tem6,tem7,mp_tem)

    IF( lvldbg >= 4) THEN
      CALL checkshx(padv, nx,ny,nz,1,nx-1,1,ny-1,2,nz-2,                &
                    'padvx', tem3)
      CALL checkshy(padv, nx,ny,nz,1,nx-1,1,ny-1,2,nz-2,                &
                    'padvy', tem3)
    END IF
!
!-----------------------------------------------------------------------
!
!  Store the total forcing in array pforce.
!
!-----------------------------------------------------------------------
!
!
    DO k=2,nz-2
      DO j=1,ny-1
        DO i=1,nx-1
          pforce(i,j,k)= -padv(i,j,k)
        END DO
      END DO
    END DO

  ELSE

    DO k=2,nz-2
      DO j=1,ny-1
        DO i=1,nx-1
          pforce(i,j,k)= 0.0
        END DO
      END DO
    END DO

  END IF

!
!-----------------------------------------------------------------------
!
!  To calculate additional boundary relaxation and mixing terms for
!  pressure in the case of externally forced boundary condition.
!
!-----------------------------------------------------------------------
!
  IF ( lbcopt == 2 .AND. mgrid == 1 ) THEN

    CALL set_acct(bc_acct)

    CALL brlxp(nx,ny,nz, dtbig1, pprt(1,1,1,1),rhostr, pforce,          &
               exbcbuf(npr0exb),exbcbuf(nprdtexb), bcrlx,               &
               tem1,tem2,tem3,tem4)

  END IF

  IF ( lbcopt == 1 .AND. mgrid == 1 ) THEN
    CALL brlxs_rbc( nx,ny,nz,pprt(1,1,1,1),rhostr, pforce, bcrlx )
  END IF

  RETURN
END SUBROUTINE frcp
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE FRCPT                      ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE frcpt(nx,ny,nz, exbcbufsz, dtbig1,ptprt,u,v,w,wcont,         &
           ptbar,rhostr,rhostri,kmh,kmv,rprntl,                         &
           usflx,vsflx,ptsflx,pbldpth,                                  &
           x,y,z,zp,mapfct,j1,j2,j3,aj3x,aj3y,j3inv,ptsfc,              &
           ptforce,                                                     &
           exbcbuf, bcrlx,                                              &
           tem1,tem2,tem3,tem4,tem5,tem6,tem7,                          &
           tem8,tem9,tem10,tem11,                                       &
           tem1_0,tem2_0,tem3_0,mp_tem)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate gravity wave or inactive acoustic wave  terms in the
!  potential temperature equation.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  9/17/94
!
!  MODIFICATION HISTORY:
!
!  8/15/95 (Ming Xue)
!  Corrected a bug related the calculation of addtional boundary zone
!  relaxation. The ptmix term was effectively added twice to
!  ptforce when lbcopt=2.
!
!  01/23/96 (Donghai Wang, Yuhe Liu and Ming Xue)
!  Added the map factor to forcing terms.
!
!  4/1/96 (Donghai Wang, X. Song and M. Xue)
!  Added the implicit treatment for the vertical mixing.
!
!  7/10/1997 (Fanyou Kong -- CMRP)
!  Added the positive definite advection option (sadvopt = 5)
!
!  7/17/1998 (M. Xue)
!  Changed call to ADVPTFCT.
!
!  9/18/1998 (D. Weber)
!  Added arrays aj3x,y.
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!    dtbig1   The big time step size (s)
!
!    ptprt    Perturbation potential temperature at times tpast and
!             tpresent (K)
!
!    u        x component of velocity at all time levels (m/s)
!    v        y component of velocity at all time levels (m/s)
!    w        Vertical component of Cartesian velocity at times
!             tpast and tpresent (m/s)
!    wcont    Contravariant vertical velocity (m/s)
!
!    ptbar    Base state potential temperature (K)
!    rhostr   Base state density rhobar times j3 (kg/m**3)
!    rhostri  Inv. base state density rhobar times j3 (kg/m**3)
!
!    kmh      Horizontal turb. mixing coef. for momentum ( m**2/s )
!    kmv      Vertical turb. mixing coef. for momentum ( m**2/s )
!    rprntl   Reciprocal of Prandtl number
!
!    ptsflx   Surface flux of heat (K*kg/(m**2*s)).
!    pbldpth  Planetary boundary layer depth (m)
!
!    x        x coordinate of grid points in physical/comp. space (m)
!    y        y coordinate of grid points in physical/comp. space (m)
!    z        z coordinate of grid points in computational space (m)
!    zp       Vertical coordinate of grid points in physical space(m)
!
!    mapfct   Map factors at scalar points
!
!    j1       Coordinate transformation Jacobian -d(zp)/dx
!    j2       Coordinate transformation Jacobian -d(zp)/dy
!    j3       Coordinate transformation Jacobian  d(zp)/dz
!    aj3x     Avgx of the coordinate transformation Jacobian  d(zp)/dz
!    aj3y     Avgy of the coordinate transformation Jacobian  d(zp)/dz
!
!  OUTPUT:
!
!    ptforce  Gravity wave inactive forcing terms in pt-eq.
!             (K*kg/(m**3*s))
!
!  WORK ARRAYS:
!
!    tem1     Temporary work array.
!    tem2     Temporary work array.
!    tem3     Temporary work array.
!    tem4     Temporary work array.
!    tem5     Temporary work array.
!    tem6     Temporary work array.
!    tem7     Temporary work array.
!    tem8     Temporary work array.
!    tem9     Temporary work array.
!    tem10    Temporary work array.
!    tem11    Temporary work array.
!
!    tem1_0   Temporary work array.
!    tem2_0   Temporary work array.
!    tem3_0   Temporary work array.
!
!    mp_tem   Temporary work array.
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
  IMPLICIT NONE             ! Force explicit declarations

  INCLUDE 'timelvls.inc'

  INTEGER :: nx, ny, nz        ! Number of grid points in 3 directions

  REAL :: dtbig1               ! Local big time step size

  REAL :: ptprt (nx,ny,nz,nt)  ! Perturbation potential temperature (K)

  REAL :: u     (nx,ny,nz,nt)  ! Total u-velocity (m/s)
  REAL :: v     (nx,ny,nz,nt)  ! Total v-velocity (m/s)
  REAL :: w     (nx,ny,nz,nt)  ! Total w-velocity (m/s)
  REAL :: wcont (nx,ny,nz)     ! Contravariant vertical velocity (m/s)

  REAL :: ptbar (nx,ny,nz)     ! Base state potential temperature (K)
  REAL :: rhostr(nx,ny,nz)     ! Base state density rhobar times j3.
  REAL :: rhostri(nx,ny,nz)    ! Inv. base state density rhobar times j3.

  REAL :: kmh   (nx,ny,nz)     ! Horizontal turb. mixing coef. for
                               ! momentum. ( m**2/s )
  REAL :: kmv   (nx,ny,nz)     ! Vertical turb. mixing coef. for
                               ! momentum. ( m**2/s )
  REAL :: rprntl(nx,ny,nz)     ! Reciprocal of Prandtl number

  REAL :: usflx (nx,ny)        ! Surface flux of u-momentum
                               ! (kg/(m*s**2))
  REAL :: vsflx (nx,ny)        ! Surface flux of v-momentum
                               ! (kg/(m*s**2))
  REAL :: ptsflx(nx,ny)        ! Surface flux of heat (K*kg/(m**2*s))
  REAL :: pbldpth(nx,ny,nt)    ! Planetary boundary layer depth (m)

  REAL :: x     (nx)           ! x-coord. of the physical and compu-
                               ! tational grid. Defined at u-point.
  REAL :: y     (ny)           ! y-coord. of the physical and compu-
                               ! tational grid. Defined at v-point.
  REAL :: z     (nz)           ! z-coord. of the computational grid.
                               ! Defined at w-point on the staggered
                               ! grid.
  REAL :: zp    (nx,ny,nz)     ! Physical height coordinate defined at
                               ! w-point of the staggered grid.

  REAL :: mapfct(nx,ny,8)      ! Map factors at scalar points

  REAL :: j1    (nx,ny,nz)     ! Coordinate transformation Jacobian
                               ! defined as - d( zp )/d( x ).
  REAL :: j2    (nx,ny,nz)     ! Coordinate transformation Jacobian
                               ! defined as - d( zp )/d( y ).
  REAL :: j3    (nx,ny,nz)     ! Coordinate transformation Jacobian
                               ! defined as d( zp )/d( z ).
  REAL :: aj3x  (nx,ny,nz)     ! Coordinate transformation Jacobian defined
                               ! as d( zp )/d( z ) AVERAGED IN THE X-DIR.
  REAL :: aj3y  (nx,ny,nz)     ! Coordinate transformation Jacobian defined
                               ! as d( zp )/d( z ) AVERAGED IN THE Y-DIR.
  REAL :: j3inv (nx,ny,nz)     ! Coordinate transformation Jacobian
                               ! defined as d( zp )/d( z ).
  REAL :: ptsfc  (nx,ny)       ! Ground surface potential temperature (K)

  REAL :: ptforce(nx,ny,nz)    ! Gravity wave inactive forcing terms
                               ! in pressure equation (Pascal/s)

  INTEGER :: exbcbufsz         ! EXBC buffer size
  REAL :: exbcbuf( exbcbufsz ) ! EXBC buffer array
  REAL :: bcrlx(nx,ny)         ! EXBC relaxation coefficients

  REAL :: tem1  (nx,ny,nz)     ! Temporary work array
  REAL :: tem2  (nx,ny,nz)     ! Temporary work array
  REAL :: tem3  (nx,ny,nz)     ! Temporary work array
  REAL :: tem4  (nx,ny,nz)     ! Temporary work array
  REAL :: tem5  (nx,ny,nz)     ! Temporary work array
  REAL :: tem6  (nx,ny,nz)     ! Temporary work array
  REAL :: tem7  (nx,ny,nz)     ! Temporary work array
  REAL :: tem8  (nx,ny,nz)     ! Temporary work array
  REAL :: tem9  (nx,ny,nz)     ! Temporary work array
  REAL :: tem10 (nx,ny,nz)     ! Temporary work array
  REAL :: tem11 (nx,ny,nz)     ! Temporary work array

  REAL :: tem1_0(0:nx,0:ny,0:nz)     ! Temporary work array.
  REAL :: tem2_0(0:nx,0:ny,0:nz)     ! Temporary work array.
  REAL :: tem3_0(0:nx,0:ny,0:nz)     ! Temporary work array.

  REAL :: mp_tem(MAX(nx+1,ny+1)*(nz+1))  ! Temporary message passing array.
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i, j, k, tstrtlvl
  REAL :: deltat
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'bndry.inc'
  INCLUDE 'exbc.inc'
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
!  Compute the advection term of the potential temperature
!  equation and store it in array ptadv.
!
!-----------------------------------------------------------------------
!
  CALL set_acct(advs_acct)

  IF( sadvopt == 4 .and. (tintegopt == 1 .or. tintegopt == 3)) THEN  ! FCT advection

    CALL advptfct(nx,ny,nz,dtbig1,ptprt,u,v,w,wcont,                    &
                  rhostr,rhostri,ptbar,mapfct,j3,j3inv,                 &
                  ptforce,                                              &
                  tem1,tem2,tem3,tem4,tem5,tem6,tem7,tem8,              &
                  tem9,tem10,tem11,tem1_0,tem2_0,tem3_0,mp_tem)

    deltat = dtbig1
    tstrtlvl = tpresent

  ELSE

    IF(tintegopt == 1 .or. tintegopt == 3) THEN

      CALL advpt(nx,ny,nz,ptprt,u,v,w,wcont, rhostr,ptbar,mapfct,         &
                 j3,j3inv,ptforce,                                        &
                 tem1,tem2,tem3,tem4,tem5,tem6,tem7,tem8)

      IF(tintegopt == 1) THEN
        deltat = dtbig1*2
        tstrtlvl = tpast
      ELSE
        deltat = dtbig1
        tstrtlvl = tpresent
      END IF

    ELSE IF(tintegopt == 2) THEN

      deltat = dtbig1
      tstrtlvl = tpresent

    END IF

  END IF

!
!-----------------------------------------------------------------------
!
!  Compute the mixing terms in the potential temperature equation.
!  This includes both physical and computational mixing.
!  Store in array ptmix.
!
!-----------------------------------------------------------------------
!
  CALL set_acct(cmix_acct)
  CALL mixpt(nx,ny,nz, exbcbufsz,                                       &
             ptprt(1,1,1,tstrtlvl),ptbar,rhostr,                        &
             kmh,kmv,rprntl,                                            &
             usflx,vsflx,ptsflx,pbldpth(1,1,tpresent),                  &
             x,y,z,zp,mapfct,j1,j2,j3,aj3x,aj3y,j3inv,ptsfc,            &
             tem7, exbcbuf,                                             &
             tem1,tem2,tem3,tem4,tem5,tem6)

  DO k=1,nz-1
    DO j=1,ny-1
      DO i=1,nx-1
        IF(tintegopt == 1 .or. tintegopt == 3) THEN
          ptforce(i,j,k)=-ptforce(i,j,k)+tem7(i,j,k)
        ELSE IF(tintegopt == 2) THEN
          ptforce(i,j,k)=tem7(i,j,k)
        END IF
        ! DTD: store mixing term (tem7) in tem11 for later dumping
        tem10(i,j,k) = ptforce(i,j,k)
        tem11(i,j,k) = tem7(i,j,k)
      END DO
    END DO
  END DO
!
!-----------------------------------------------------------------------
!
!  Calculate additional relaxation and mixing on ptprt for the
!  base grid when external boundary forcing is used.
!
!-----------------------------------------------------------------------
!
  IF ( lbcopt == 2 .AND. mgrid == 1 ) THEN

    CALL set_acct(bc_acct)

    CALL brlxpt(nx,ny,nz,deltat*0.5,ptprt(1,1,1,tstrtlvl),              &
                rhostr,ptforce,                                         &
                exbcbuf(npt0exb),exbcbuf(nptdtexb),bcrlx,               &
                tem1,tem2,tem3,tem4)
  END IF

  IF ( lbcopt == 1 .AND. mgrid == 1 ) THEN
    CALL brlxs_rbc( nx,ny,nz,ptprt(1,1,1,tstrtlvl),rhostr, ptforce, bcrlx )
  END IF

!
!-----------------------------------------------------------------------
!
!  Treat the vertically implicit mixing
!
!  This is only called here if tintegopt == 1 (the original leapfrog
!  time-stepping) or tintegopt == 3 (RK3 with all terms updated).
!  For tintegopt = 2 (RK3 case with only advection terms updated),
!  it will be calculated outside this subroutine.
!
!-----------------------------------------------------------------------
!

  IF (trbvimp == 1 .and. (tintegopt == 1 .or. tintegopt == 3)) THEN     ! Vertical implicit application

    CALL set_acct(tmix_acct)

    DO k=1,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          tem1(i,j,k)=ptbar(i,j,k)+ptprt(i,j,k,tstrtlvl)
        END DO
      END DO
    END DO

    DO k=1,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          tem2(i,j,k)=rhostr(i,j,k)*kmv(i,j,k)*rprntl(i,j,k)            &
                      *j3inv(i,j,k)*j3inv(i,j,k)
        END DO
      END DO
    END DO

    CALL vmiximps(nx,ny,nz,deltat*0.5,tem1,rhostr,tem2,                 &
                  ptforce,tem3,tem4,tem5,tem6)

  END IF

  !DTD: Retrieve the completed mixing term and store in tem11 for later dumping

  DO k=1,nz-1
    DO j=1,ny-1
      DO i=1,nx-1
        tem10(i,j,k) = ptforce(i,j,k)-tem10(i,j,k) ! tem10 now contains the vertically implicit mixing term
        tem11(i,j,k) = tem11(i,j,k) + tem10(i,j,k) ! tem11 now contains the complete mixing term
        tem11(i,j,k) = tem11(i,j,k)*rhostri(i,j,k)
      END DO
    END DO
  END DO

  RETURN
END SUBROUTINE frcpt
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE CORIOL                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE coriol(nx,ny,nz,                                             &
           u,v,w,ubar,vbar,rhostr, sinlat, mapfct,                      &
           uforce, vforce, wforce,                                      &
           fcorio1, fcorio2, tem1,tem2,tem3,tem4)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate the Coriolis force terms in the u, v and w equations.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  4/20/92.
!
!  MODIFICATION HISTORY:
!
!  5/05/92 (M. Xue)
!  Added full documentation.
!
!  8/28/92 (M. Xue)
!  Included the domain translation effect into the Coriolis force
!  terms.
!
!  9/10/94 (D. Weber & Y. Lu)
!  Cleaned up documentation.
!
!  9/14/98 (D. Weber)
!  Removed operators and merged loops and added coriolis terms
!  to the forcing terms inside this subroutine.
!
!  10/23/2003 (Ming Xue & Y. Wang)
!  Added effects of spatial gradient of map factor and of earth
!  curvative
!
!  NOTE: fcoro1 is changed to a 3D array now.
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical direction.
!
!    u        x component of velocity at a given time level (m/s)
!    v        y component of velocity at a given time level (m/s)
!    w        z component of velocity at a given time level (m/s)
!    ubar     Base state x velocity component (m/s)
!    vbar     Base state y velocity component (m/s)
!    rhostr   Base state density rhobar times j3 (kg/m**3)
!    sinlat   Sin of latitude at each grid point
!    mapfct   Map factor at scalar, u and v point
!
!  OUTPUT:
!
!    uforce   Acoustically inactive forcing terms in u-momentum
!             equation (kg/(m*s)**2). uforce= -uadv + umix + ucorio
!    vforce   Acoustically inactive forcing terms in v-momentum
!             equation (kg/(m*s)**2). vforce= -vadv + vmix + vcorio
!    wforce   Acoustically inactive forcing terms in w-momentum
!             equation (kg/(m*s)**2).
!
!  WORK ARRAYS:
!
!    fcorio1  Temporary work array.
!    fcorio2  Temporary work array.
!    tem1     Temporary work array.
!    tem2     Temporary work array.
!
!-----------------------------------------------------------------------
!

!
!-----------------------------------------------------------------------
!
!  Variable Declarations
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx,ny,nz          ! Number of grid points in 3 directions

  REAL :: u     (nx,ny,nz)     ! Total u-velocity at a time level (m/s)
  REAL :: v     (nx,ny,nz)     ! Total v-velocity at a time level (m/s)
  REAL :: w     (nx,ny,nz)     ! Total w-velocity at a time level (m/s)

  REAL :: ubar  (nx,ny,nz)     ! Base state u-velocity (m/s)
  REAL :: vbar  (nx,ny,nz)     ! Base state v-velocity (m/s)
  REAL :: rhostr(nx,ny,nz)     ! Base state density rhobar times j3.
  REAL :: sinlat(nx,ny)        ! Sin of latitude at each grid point
  REAL :: mapfct(nx,ny,8)      ! Map factor at scalar (1), U(2) and V(3)
                               ! point and the corresponding reverse (4-6)
                               ! and the square (7), quarter (8)

  REAL :: uforce(nx,ny,nz)     ! Acoustically inactive forcing terms
                               ! in u-momentum equation (kg/(m*s)**2)
                               ! uforce= -uadv + umix + ucorio

  REAL :: vforce(nx,ny,nz)     ! Acoustically inactive forcing terms
                               ! in v-momentum equation (kg/(m*s)**2)
                               ! vforce= -vadv + vmix + vcorio

  REAL :: wforce(nx,ny,nz)     ! Acoustically inactive forcing terms
                               ! in w-momentum equation (kg/(m*s)**2)
                               ! wforce= -wadv + wmix + wbuoy

  REAL :: fcorio1  (nx,ny,nz)  ! Temporary work array.
  REAL :: fcorio2  (nx,ny)     ! Temporary work array.
  REAL :: tem1  (nx,ny,nz)     ! Temporary work array.
  REAL :: tem2  (nx,ny,nz)     ! Temporary work array.
  REAL :: tem3  (nx,ny,nz)     ! Temporary work array.
  REAL :: tem4  (nx,ny,nz)     ! Temporary work array.

!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k
  REAL :: omega2,sinclat,cosclat
  REAL :: gumove,gvmove
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'grid.inc'          ! Grid parameters
  INCLUDE 'phycst.inc'
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  IF ( grdtrns == 0 ) THEN
    gumove = 0.0
    gvmove = 0.0
  ELSE
    gumove = umove
    gvmove = vmove
  END IF

  omega2 = 2.0* omega

  ! Added the effects of spatial gradient of map factor on the coriolis force
  !
  IF( coriopt == 3 .OR. coriopt == 4) THEN
    DO j = 1,ny
      DO i = 1,nx
        tem3(i,j,1) = sinlat(i,j)/SQRT(1-sinlat(i,j)**2)  ! tan(lat)
      END DO
    END DO

    DO k = 1,nz
      DO j = 1,ny-1
        DO i = 1,nx-1
          !
          ! fm = U*My - V*Mx + U*TAN(lat)/eradius
          !
          tem4(i,j,k) = 0.5*( (u(i,j,k)+u(i+1,j,k))*                    &
                              ((mapfct(i,j+1,3)-mapfct(i,j,3))/dy)      &
                             -(v(i,j,k)+v(i,j+1,k))*                    &
                              ((mapfct(i+1,j,2)-mapfct(i,j,2))/dx)      &
                             +((u(i,j,k)+u(i+1,j,k))*tem3(i,j,1))/      &
                              eradius)
        END DO
      END DO
    END DO

  END IF

  IF( coriopt == 1 ) THEN

    sinclat = SIN( ATAN(1.0)/45.0 * ctrlat )
    DO k = 1,nz
      DO j=1,ny
        DO i=1,nx
          fcorio1(i,j,k) = omega2* sinclat
          fcorio2(i,j)   = 0.0
        END DO
      END DO
    END DO

  ELSE IF( coriopt == 2 ) THEN

    sinclat = SIN( ATAN(1.0)/45.0 * ctrlat )
    cosclat = SQRT( 1.0 - sinclat*sinclat  )
    DO k = 1,nz
      DO j=1,ny
        DO i=1,nx
          fcorio1(i,j,k) = omega2* sinclat
          fcorio2(i,j) = omega2* cosclat
        END DO
      END DO
    END DO

  ELSE IF( coriopt == 3 ) THEN

    DO k = 1,nz
      DO j=1,ny-1
        DO i=1,nx-1
          fcorio1(i,j,k) = omega2*sinlat(i,j) + tem4(i,j,k)
          fcorio2(i,j) = 0.0
        END DO
      END DO
    END DO

  ELSE IF( coriopt == 4 ) THEN

    DO k = 1,nz
      DO j=1,ny-1
        DO i=1,nx-1
          fcorio1(i,j,k) = omega2* sinlat(i,j) + tem4(i,j,k)
          fcorio2(i,j) = omega2* SQRT( 1.0 - sinlat(i,j)**2 )
        END DO
      END DO
    END DO

  END IF

!-----------------------------------------------------------------------
!
!  Coriolis terms in u-eq. if coriotrm=1
!
!    ucorio = avgx( rhostr * (fcorio1*avgy(v)-fcorio2*avgz(w)) )
!
!  or if coriotrm=2,
!
!    ucorio = avgx( rhostr * (fcorio1*avgy(v-vbar)-fcorio2*avgz(w)) )
!
!  where fcorio1 = 2*omega*sinlat.
!  and   fcorio2 = 2*omega* sqrt(1-sinlat**2).
!
!-----------------------------------------------------------------------

  IF( coriotrm == 1 ) THEN
    DO k=1,nz-1
      DO j=1,ny
        DO i=1,nx-1
          tem2(i,j,k) = v(i,j,k) + gvmove
        END DO
      END DO
    END DO
  ELSE
    DO k=1,nz-1
      DO j=1,ny
        DO i=1,nx-1
          tem2(i,j,k) = v(i,j,k) - vbar(i,j,k)
          ! dropped gvmove, detected by Ernani
        END DO
      END DO
    END DO
  END IF

  DO k=1,nz-1
    DO j=1,ny-1
      DO i=1,nx-1
        tem1(i,j,k) = 0.5*rhostr(i,j,k)*                                &
                        (fcorio1(i,j,k)*(tem2(i,j+1,k)+tem2(i,j,k))     &
                        -fcorio2(i,j)  *(w(i,j,k+1)+w(i,j,k)))
      END DO
    END DO
  END DO

  !
  ! earth curvature term - U*W/eradius
  !
  IF(earth_curvature == 1 .AND. (coriopt == 3 .OR. coriopt == 4) ) THEN
    DO k = 1,nz-1
      DO j = 1,ny-1
        DO i = 1,nx-1
           tem1(i,j,k) = tem1(i,j,k) - 0.25*rhostr(i,j,k)*              &
                         (u(i,j,k)+u(i+1,j,k))*(w(i,j,k)+w(i,j,k+1))/   &
                         eradius
        END DO
      END DO
    END DO
  END IF

  DO k=1,nz-1
    DO j=1,ny-1
      DO i=2,nx-1
        uforce(i,j,k) = uforce(i,j,k) + 0.5*(tem1(i-1,j,k)+tem1(i,j,k))
      END DO
    END DO
  END DO

!-----------------------------------------------------------------------
!
!  Coriolis terms in v-eq. and w-eq. if coriotrm=1:
!
!    vcorio = - avgy( rhostr * fcorio1 * avgx(u) )
!    wcorio = avgz( rhostr * fcorio2 * avgx(u) )
!
!  or if coriotrm=2:
!
!    vcorio = - avgy( rhostr * fcorio1 * avgx(u-ubar) )
!    wcorio = avgz( rhostr * fcorio2 * avgx(u) )
!
!
!-----------------------------------------------------------------------

  IF( coriotrm == 1 ) THEN
    DO k=1,nz-1
      DO j=1,ny-1
        DO i=1,nx
          tem2(i,j,k) = u(i,j,k) + gumove
        END DO
      END DO
    END DO
  ELSE
    DO k=1,nz-1
      DO j=1,ny-1
        DO i=1,nx
          tem2(i,j,k) = u(i,j,k) - ubar(i,j,k)
          ! dropped gumove, detected by Ernani
        END DO
      END DO
    END DO
  END IF

  DO k=1,nz-1
    DO j=1,ny-1
      DO i=1,nx-1
        tem1(i,j,k) = - rhostr(i,j,k)*fcorio1(i,j,k)*                   &
                        0.5*(tem2(i+1,j,k)+tem2(i,j,k))
      END DO
    END DO
  END DO

  !
  ! earth curvature term - V*W/eradius
  !
  IF(earth_curvature == 1 .AND. (coriopt == 3 .OR. coriopt == 4) ) THEN
    DO k = 1,nz-1
      DO j = 1,ny-1
        DO i = 1,nx-1
           tem1(i,j,k) = tem1(i,j,k) - 0.25*rhostr(i,j,k)*              &
                         (v(i,j,k)+v(i,j+1,k))*(w(i,j,k)+w(i,j,k+1))/   &
                         eradius
        END DO
      END DO
    END DO
  END IF

  DO k=1,nz-1
    DO j=2,ny-1
      DO i=1,nx-1
        vforce(i,j,k) = vforce(i,j,k)+0.5*(tem1(i,j-1,k)+tem1(i,j,k))
      END DO
    END DO
  END DO

!-----------------------------------------------------------------------
!
! Compute wforce below
!
!-----------------------------------------------------------------------

  IF( coriopt == 2 .OR. coriopt == 4 ) THEN ! compute the w term.

    DO k=1,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          tem1(i,j,k) = fcorio2(i,j) * rhostr(i,j,k)*                   &
                          ( 0.5*(u(i+1,j,k)+u(i,j,k)) + gumove )
        END DO
      END DO
    END DO

    !
    ! earth curvature term - (U**2+V**2)/a
    !
    IF( earth_curvature == 1 .AND. coriopt == 4 ) THEN
      DO k = 1,nz-1
        DO j = 1,ny-1
          DO i = 1,nx-1
            tem1(i,j,k) = tem1(i,j,k) + 0.25*rhostr(i,j,k)*             &
                          ( (u(i,j,k)+u(i+1,j,k))**2+                   &
                            (v(i,j,k)+v(i,j+1,k))**2 ) / eradius
          END DO
        END DO
      END DO
    END IF

    DO k=2,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          wforce(i,j,k) = wforce(i,j,k)+0.5*(tem1(i,j,k-1)+tem1(i,j,k))
        END DO
      END DO
    END DO

  END IF

  RETURN
END SUBROUTINE coriol

!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE BUOYCY                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE buoycy(nx,ny,nz,ptprt,pprt,qv,qscalar,                &
           ptbar,pbar,ptbari,pbari,rhostr,qvbar, wbuoy, tem1)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate the total buoyancy including liquid and solid water
!  loading.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  10/10/91.
!
!  MODIFICATION HISTORY:
!
!  5/05/92 (M. Xue)
!  Added full documentation.
!
!  3/10/93 (M. Xue)
!  The buoyancy term is reformulated. The previous formula was
!  in error. The water loading was calculated wrong, resulting in
!  a value of the water loading that is typically an order of
!  magnitude too small.
!
!  3/25/94 (G. Bassett & M. Xue)
!  The buoyancy terms are reformulated for better numerical accuracy.
!  Instead of storing numbers which had the form (1+eps)*(1+eps1)
!  (eps << 1 and eps1 <<1), terms were expanded out, and most of the
!  high order terms neglected, except for the second order terms
!  in ptprt, pprt and qvbar.
!
!  9/10/94 (D. Weber & Y. Lu)
!  Cleaned up documentation.
!
!  6/21/95 (Alan Shapiro)
!  Fixed bug involving missing qvpert term in buoyancy formulation.
!
!  10/15/97 (Donghai Wang)
!  Added a new option for including the second order terms.
!
!  11/05/97 (D. Weber)
!  Changed lower loop bounds in DO LOOP 400 for computing the
!  buoyancy term from k=3,nz-2 to k=2,nz-1.  Level k=2 data will be
!  used in the hydrostatic pprt lower boundary condition (removed
!  DO LOOP 410 used to set wbuoy = 0.0 for k= 2 and nz-1).
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical direction.
!
!    ptprt    Perturbation potential temperature at a time level (K)
!    pprt     Perturbation pressure at a given time level (Pascal)
!    qv       Water vapor specific humidity at a given time level
!             (kg/kg)
!    qc       Cloud water mixing ratio at a given time level (kg/kg)
!    qr       Rainwater mixing ratio at a given time level (kg/kg)
!    qi       Cloud ice mixing ratio at a given time level (kg/kg)
!    qs       Snow mixing ratio at a given time level (kg/kg)
!    qh       Hail mixing ratio at a given time level (kg/kg)
!
!    ptbar    Base state potential temperature (K)
!    pbar     Base state pressure (Pascal)
!    ptbari   Inverse base state potontial temperature (K)
!    pbari    Inverse base state pressure (Pascal)
!    rhostr   Base state density rhobar times j3 (kg/m**3)
!    qvbar    Base state water vapor specific humidity (kg/kg)
!
!  OUTPUT:
!
!    wbuoy    The total buoyancy force (kg/(m*s)**2)
!
!  WORK ARRAYS:
!
!    tem1     Temporary work array.
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
  INCLUDE 'globcst.inc'     ! Global model control parameters
  INCLUDE 'phycst.inc'      ! Physical constants

!
!-----------------------------------------------------------------------
!
!  Variable Declarations
!
!-----------------------------------------------------------------------
!
  INTEGER :: nx,ny,nz          ! Number of grid points in 3 directions

  REAL :: ptprt (nx,ny,nz)     ! Perturbation potential temperature
                               ! at a given time level (K)
  REAL :: pprt  (nx,ny,nz)     ! Perturbation pressure at a given time
                               ! level (Pascal)
  REAL :: qv    (nx,ny,nz)     ! Water vapor specific humidity (kg/kg)

  REAL :: qscalar    (nx,ny,nz,nscalar)

  REAL :: ptbar (nx,ny,nz)     ! Base state potential temperature (K)
  REAL :: pbar  (nx,ny,nz)     ! Base state pressure (Pascal).
  REAL :: ptbari(nx,ny,nz)     ! Inverse base state pot. temperature (K)
  REAL :: pbari (nx,ny,nz)     ! Inverse base state pressure (Pascal).
  REAL :: rhostr(nx,ny,nz)     ! Base state density rhobar times j3.
  REAL :: qvbar (nx,ny,nz)     ! Base state water vapor specific
                               ! humidity(kg/kg)

  REAL :: wbuoy(nx,ny,nz)      ! Total buoyancy in w-eq. (kg/(m*s)**2)

  REAL :: tem1  (nx,ny,nz)     ! Temporary work array.

!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k
  INTEGER :: nq
  REAL :: g5
  REAL :: pttem,tema
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
!  The total buoyancy
!
!    wbuoy = rhostr*g ( ptprt/ptbar-pprt/(rhobar*csndsq)+
!    qvprt/(rddrv+qvbar)-(qvprt+qc+qr+qs+qi+qh)/(1+qvbar)
!    -(ptprt*ptprt)/(ptbar*ptbar)                        !2nd-order
!    +0.5*(ptprt*pprt)/(cpdcv*ptbar*pbar))               !2nd-order
!
!  and rddrv=rd/rv, cp, cv, rd and rv are defined in phycst.inc.
!
!  Here, the contribution from pprt (i.e., term pprt/(rhobar*csndsq))
!  is evaluated inside the small time steps, therefore wbuoy
!  does not include this part.
!
!  The contribution from ptprt is calculated inside the small time
!  steps if the potential temperature equation is solved inside
!  small time steps, i.e., if ptsmlstp=1.
!
!-----------------------------------------------------------------------
!
  IF( ptsmlstp == 1 ) THEN

    DO k=1,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          tem1(i,j,k) = 0.0
        END DO
      END DO
    END DO

  ELSE
    IF (buoy2nd == 0) THEN  !1st-order

      DO k=1,nz-1
        DO j=1,ny-1
          DO i=1,nx-1
            tem1(i,j,k) = ptprt(i,j,k)*ptbari(i,j,k)
          END DO
        END DO
      END DO

    ELSE                          !2nd-order
      IF ( bsnesq == 1 ) THEN

        DO k=1,nz-1
          DO j=1,ny-1
            DO i=1,nx-1
              pttem = ptprt(i,j,k)*ptbari(i,j,k)
              tem1(i,j,k) = pttem-pttem*pttem
            END DO
          END DO
        END DO


      ELSE

        tema = 1.0/cpdcv
        DO k=1,nz-1
          DO j=1,ny-1
            DO i=1,nx-1
              pttem = ptprt(i,j,k)*ptbari(i,j,k)
              tem1(i,j,k) = pttem*                                      &
                  (1.0-pttem+0.5*pprt(i,j,k)*(tema*pbari(i,j,k)))
            END DO
          END DO
        END DO

      END IF
    END IF
  END IF

!
!-----------------------------------------------------------------------
!
!  Add on the contributions to the buoyancy from the water vapor
!  content and the liquid and ice water loading.
!
!-----------------------------------------------------------------------
!

!  IF( moist == 1 .AND. ice == 0 ) THEN  ! Moist case, no ice.
!
!    DO k=1,nz-1
!      DO j=1,ny-1
!        DO i=1,nx-1
!          tem1(i,j,k) = tem1(i,j,k)                                     &
!              + (qv(i,j,k) - qvbar(i,j,k))/(rddrv + qvbar(i,j,k))       &
!              - (qv(i,j,k) - qvbar(i,j,k) + qc(i,j,k) + qr(i,j,k))      &
!              /(1 + qvbar(i,j,k))
!        END DO
!      END DO
!    END DO
!
!  ELSE IF(moist == 1 .AND. ice == 1) THEN
!                              ! Full microphysics case, loading
!! of liquid and ice water.
!    DO k=1,nz-1
!      DO j=1,ny-1
!        DO i=1,nx-1
!          tem1(i,j,k) = tem1(i,j,k)                                     &
!              + (qv(i,j,k) - qvbar(i,j,k))/(rddrv + qvbar(i,j,k))       &
!              - (qv(i,j,k) - qvbar(i,j,k) + qc(i,j,k) + qr(i,j,k) +     &
!              qs(i,j,k) + qi(i,j,k) + qh(i,j,k))/(1 + qvbar(i,j,k))
!        END DO
!      END DO
!    END DO
!
!  END IF

  IF (moist == 1) THEN

    DO k=1,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          tem1(i,j,k) = tem1(i,j,k)                                     &
              + (qv(i,j,k) - qvbar(i,j,k))/(rddrv + qvbar(i,j,k))       &
              - (qv(i,j,k) - qvbar(i,j,k))/(1 + qvbar(i,j,k))
        END DO
      END DO
    END DO

    DO k=1,nz-1
    DO j=1,ny-1
    DO i=1,nx-1
      DO nq = 1,nscalarq
        tem1(i,j,k) = tem1(i,j,k) - (qscalar(i,j,k,nq))/(1+qvbar(i,j,k))
      END DO
    END DO
    END DO
    END DO

  END IF
!
!-----------------------------------------------------------------------
!
!  Then the total buoyancy:
!
!    wbuoy = tem1 * rhostr * g
!
!  averged to the w-point on the staggered grid.
!
!-----------------------------------------------------------------------
!
  g5 = g*0.5

  DO k=2,nz-1
    DO j=1,ny-1
      DO i=1,nx-1
        wbuoy(i,j,k)= (tem1(i,j, k )*rhostr(i,j, k )                    &
                      +tem1(i,j,k-1)*rhostr(i,j,k-1))*g5
      END DO
    END DO
  END DO

  RETURN
END SUBROUTINE buoycy
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE PGRAD                      ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE pgrad(nx,ny,nz, pprt, div,                                   &
           j1,j2,j3, upgrad,vpgrad,wpgrad,tem1,tem2)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate the pressure gradient terms in the momentum equations.
!  These terms are evaluated every small time step.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  10/10/91.
!
!  MODIFICATION HISTORY:
!
!  5/05/92 (M. Xue)
!  Added full documentation.
!
!  4/10/93 (M. Xue & Hao Jin)
!  Add the terrain.
!
!  5/25/93 (M. Xue & K. Brewster)
!  Fixed and error in the vertical pressure gradient force term.
!  The error is present in the finite differenced w equation
!  of ARPS 3.0 User's guide. The equation in continuous form is
!  correct.
!
!  9/10/94 (D. Weber & Y. Lu)
!  Cleaned up documentation.
!
!  01/23/96 (Donghai Wang, Yuhe Liu and Ming Xue)
!  Added the map factor.
!
!  07/31/96 (Ming Xue and Yuhe Liu)
!  Added the isotropic option for divergence damping
!
!  9/14/98 (D. Weber)
!  Removed operators and merged loops.
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical direction.
!
!    pprt     Perturbation pressure at a given time level (Pascal)
!
!    j1       Coordinate transformation Jacobian -d(zp)/dx
!    j2       Coordinate transformation Jacobian -d(zp)/dy
!    j3       Coordinate transformation Jacobian  d(zp)/dz
!
!  OUTPUT:
!
!    upgrad   Pressure gradient force in u-eq. (kg/(m*s)**2)
!    vpgrad   Pressure gradient force in v-eq. (kg/(m*s)**2)
!    wpgrad   Pressure gradient force in w-eq. (kg/(m*s)**2)
!
!  WORK ARRAYS:
!
!    tem1     Temporary work array.
!    tem2     Temporary work array.
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

  INTEGER :: nx,ny,nz          ! Number of grid points in 3 directions

  REAL :: pprt  (nx,ny,nz)     ! Perturbation pressure at a given time
                               ! level (Pascal).
  REAL :: div   (nx,ny,nz)     ! Mass weighted divergence

  REAL :: j1    (nx,ny,nz)     ! Coordinate transformation Jacobian
                               ! -d(zp)/dx.
  REAL :: j2    (nx,ny,nz)     ! Coordinate transformation Jacobian
                               ! -d(zp)/dy.
  REAL :: j3    (nx,ny,nz)     ! Coordinate transformation Jacobian
                               ! d(zp)/dz.

  REAL :: upgrad(nx,ny,nz)     ! Pressure gradient force in u-eq.
                               ! (kg/(m*s)**2)
  REAL :: vpgrad(nx,ny,nz)     ! Pressure gradient force in v-eq.
                               ! (kg/(m*s)**2)
  REAL :: wpgrad(nx,ny,nz)     ! Pressure gradient force in w-eq.

  REAL :: tem1  (nx,ny,nz)     ! Temporary work array.
  REAL :: tem2  (nx,ny,nz)     ! Temporary work array.
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'grid.inc'          ! Grid parameters
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
!
!-----------------------------------------------------------------------
!
!  pprt - cdvdmpv*div(i,j,k) in vertical direction
!
!-----------------------------------------------------------------------

  DO k=1,nz-1
    DO j=1,ny-1
      DO i=1,nx-1
        tem1(i,j,k)=pprt(i,j,k)-cdvdmpv*div(i,j,k)
      END DO
    END DO
  END DO

!-----------------------------------------------------------------------
!
!  d(pprt-cdvdmpv*div)/dz at w point - p grad. force in z-dir
!
!-----------------------------------------------------------------------

  DO k=2,nz-1
    DO j=1,ny-1
      DO i=1,nx-1
        wpgrad(i,j,k)=dzinv*(tem1(i,j,k)-tem1(i,j,k-1))
      END DO
    END DO
  END DO

!-----------------------------------------------------------------------
!
!  Array div may be used only for horizontal now
!
!  div = pprt - cdvdmph*div(i,j,k) in horizontal direction
!
!-----------------------------------------------------------------------

  DO k=1,nz-1
    DO j=1,ny-1
      DO i=1,nx-1
        div(i,j,k)=pprt(i,j,k)-cdvdmph*div(i,j,k)
        tem1(i,j,k)=div(i,j,k)*j3(i,j,k)
      END DO
    END DO
  END DO

!-----------------------------------------------------------------------
!
!  d(j3*pprt)/dx at u point - 1st component of p grad. force in x-dir
!
!-----------------------------------------------------------------------

  DO k=1,nz-1
    DO j=1,ny-1
      DO i=2,nx-1
        upgrad(i,j,k)=dxinv*(tem1(i,j,k)-tem1(i-1,j,k))
      END DO
    END DO
  END DO
!-----------------------------------------------------------------------
!
!  d(j3*pprt)/dy at v point - 1st component of p grad. force in y-dir
!
!-----------------------------------------------------------------------

  DO k=1,nz-1
    DO j=2,ny-1
      DO i=1,nx-1
        vpgrad(i,j,k)=dyinv*(tem1(i,j,k)-tem1(i,j-1,k))
      END DO
    END DO
  END DO

!-----------------------------------------------------------------------
!
!  If there is no terrain, i.e. the ground is flat, skip the
!  following calculations.
!
!-----------------------------------------------------------------------

  IF( ternopt /= 0 ) THEN

!-----------------------------------------------------------------------
!
!  d(j1*pprt)/dz at u point - 2nd component of p grad. force in x-dir
!  due to terrain
!
!-----------------------------------------------------------------------

    DO k=2,nz-1
      DO j=1,ny-1
        DO i=2,nx-1
          tem2(i,j,k)= j1(i,j,k)*                                       &
                       0.25*((div(i-1,j,k-1)+div(i,j,k-1))              &
                            +(div(i-1,j,k)  +div(i,j,k)))
        END DO
      END DO
    END DO

    DO k=2,nz-2
      DO j=1,ny-1
        DO i=2,nx-1
          upgrad(i,j,k)=upgrad(i,j,k)+                                  &
                   dzinv*(tem2(i,j,k+1)-tem2(i,j,k))
        END DO
      END DO
    END DO

!-----------------------------------------------------------------------
!
!  d(j2*pprt)/dz at v point - 2nd component of p grad. force in y-dir
!  due to terrain
!
!-----------------------------------------------------------------------

    DO k=2,nz-1
      DO j=2,ny-1
        DO i=1,nx-1
          tem2(i,j,k)= j2(i,j,k)*                                       &
                       0.25*((div(i,j-1,k-1)+div(i,j,k-1))              &
                            +(div(i,j-1,k)  +div(i,j,k)))
        END DO
      END DO
    END DO

    DO k=2,nz-2
      DO j=2,ny-1
        DO i=1,nx-1
          vpgrad(i,j,k)=vpgrad(i,j,k)+                                  &
                   dzinv*(tem2(i,j,k+1)-tem2(i,j,k))
        END DO
      END DO
    END DO

  END IF

  IF( lvldbg >= 5) THEN
    CALL checkuhx(upgrad, nx,ny,nz,2,nx-1,1,ny-1,2,nz-2,                &
                  'upgrdx', tem2)
    CALL checkuhy(upgrad, nx,ny,nz,2,nx-1,1,ny-1,2,nz-2,                &
                  'upgrdy', tem2)
    CALL checkvhx(vpgrad, nx,ny,nz,1,nx-1,2,ny-1,2,nz-2,                &
                  'vpgrdx', tem2)
    CALL checkvhy(vpgrad, nx,ny,nz,1,nx-1,2,ny-1,2,nz-2,                &
                  'vpgrdy', tem2)
  END IF

  RETURN
END SUBROUTINE pgrad
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE UVWRHO                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE uvwrho(nx,ny,nz,u,v,wcont,rhostr,                            &
           ustr,vstr,wstr)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Compute ustr=u*rhostr, vstr=v*rhostr, wstr=wcont*rhostr.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  10/10/91.
!
!  MODIFICATION HISTORY:
!
!  5/05/92 (M. Xue)
!  Added full documentation.
!
!  9/10/94 (D. Weber & Y. Lu)
!  Cleaned up documentation.
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
!    wcont    Contravariant vertical velocity (m/s)
!    rhostr   Base state density rhobar times j3 (kg/m**3)
!
!  OUTPUT:
!
!    ustr     u * rhostr at u-point
!    vstr     v * rhostr at v-point
!    wstr     wcont * rhostr at w-point
!
!  WORK ARRAYS:
!
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

  INTEGER :: nx,ny,nz          ! Number of grid points in 3 directions

  REAL :: u     (nx,ny,nz)     ! Total u-velocity (m/s)
  REAL :: v     (nx,ny,nz)     ! Total v-velocity (m/s)
  REAL :: wcont (nx,ny,nz)     ! Contravariant vertical velocity (m/s)
  REAL :: rhostr(nx,ny,nz)     ! Base state density rhobar times j3.

  REAL :: ustr  (nx,ny,nz)     ! u * rhostr
  REAL :: vstr  (nx,ny,nz)     ! v * rhostr
  REAL :: wstr  (nx,ny,nz)     ! wcont * rhostr

  INTEGER :: i,j,k
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

!
!-----------------------------------------------------------------------
!
!  Calculate ustr=rhostr*u, vstr=rhostr*v, wstr=rhostr*wcont
!
!-----------------------------------------------------------------------
!
  CALL rhouvw(nx,ny,nz,rhostr,ustr,vstr,wstr)

  DO k=1,nz-1
    DO j=1,ny-1
      DO i=1,nx
        ustr(i,j,k)=u(i,j,k)*ustr(i,j,k)
      END DO
    END DO
  END DO

  DO k=1,nz-1
    DO j=1,ny
      DO i=1,nx-1
        vstr(i,j,k)=v(i,j,k)*vstr(i,j,k)
      END DO
    END DO
  END DO

  DO k=1,nz
    DO j=1,ny-1
      DO i=1,nx-1
        wstr(i,j,k)=wcont(i,j,k)*wstr(i,j,k)
      END DO
    END DO
  END DO

  RETURN
END SUBROUTINE uvwrho
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE MAPTRM                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE maptrm( nx,ny,nz, u, v, ustr,vstr, mapfct,                   &
           uforce, vforce, tem1,tem2,tem3 )
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate the term caused by inhomogeneous map factor in the
!  advection terms of u and v.
!
!  mpuadv = avgx( avgy(v*)
!               * (avgy(v)*difx(mapfct_u)-avgx(u)*dify(mapfct_v)) )
!
!  mpvadv = avgy( avgx(u*)
!               * (avgy(v)*difx(mapfct_u)-avgx(u)*dify(mapfct_v)) )
!
!  where mapfct_u is the mapfct at u points and mapfct_v at v points.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Donghai Wang, Yuhe Liu and Ming Xue
!  01/30/96
!
!  MODIFICATION HISTORY:
!
!  11/06/97 (D. Weber)
!  Added three additional levels to the mapfct array.  The three
!  levels (4,5,6) represent the inverse of the first three in order.
!  The inverse map factors are computed to improve efficiency.
!
!  9/14/98 (D. Weber)
!  Removed operators and merged loops and included the results into
!  the forcing arrays.
!
!c#######################################################################
!
!  INPUT :
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!    mapfct   Map factors at scalar, u and v points
!
!    u        x component of velocity (m/s)
!    v        y component of velocity (m/s)
!    ustr     u*rhostr
!    vstr     v*rhostr
!
!  OUTPUT:
!
!    uforce   Acoustically inactive forcing terms in u-momentum
!             equation (kg/(m*s)**2).
!    vforce   Acoustically inactive forcing terms in v-momentum
!             equation (kg/(m*s)**2).
!
!  WORK ARRAYS:
!
!    tem1     Temporary work array.
!    tem2     Temporary work array.
!    tem3     Temporary work array.
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
  INTEGER :: nx,ny,nz          ! Number of grid points in 3 directions

  REAL :: mapfct(nx,ny,8)      ! Map factors at scalar, u and v points

  REAL :: u     (nx,ny,nz)     ! Total u-velocity (m/s)
  REAL :: v     (nx,ny,nz)     ! Total v-velocity (m/s)

  REAL :: ustr  (nx,ny,nz)     ! u * rhostr
  REAL :: vstr  (nx,ny,nz)     ! v * rhostr

  REAL :: uforce(nx,ny,nz)     ! Acoustically inactive forcing terms
                               ! in u-momentum equation (kg/(m*s)**2)

  REAL :: vforce(nx,ny,nz)     ! Acoustically inactive forcing terms
                               ! in v-momentum equation (kg/(m*s)**2)

  REAL :: tem1  (nx,ny,nz)     ! Temporary work array
  REAL :: tem2  (nx,ny,nz)     ! Temporary work array
  REAL :: tem3  (nx,ny,nz)     ! Temporary work array
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k
  INTEGER :: onvf
  REAL :: tema
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'grid.inc'          ! Grid parameters
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!-----------------------------------------------------------------------
!
!  Calculate avgy(v)*difx(mapfct_u) and store it to tem3
!
!-----------------------------------------------------------------------

  tema = 0.5*dxinv
  DO k=2,nz-2
    DO j=1,ny-1
      DO i=1,nx-1
        tem3(i,j,k) = tema*(v(i,j+1,k)+v(i,j,k))*                       &
                      (mapfct(i+1,j,2)-mapfct(i,j,2))
      END DO
    END DO
  END DO

!-----------------------------------------------------------------------
!
!   At this point tem3 contains avgy(v)*difx(mapfct_u)
!
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!
!  Calculate avgx(u)*dify(mapfct_v) and substract it from tem3
!
!-----------------------------------------------------------------------

  tema = 0.5*dyinv
  DO k=2,nz-2
    DO j=1,ny-1
      DO i=1,nx-1
        tem3(i,j,k) = tem3(i,j,k) - tema*(u(i+1,j,k)+u(i,j,k))*         &
                      (mapfct(i,j+1,3)-mapfct(i,j,3))
      END DO
    END DO
  END DO

!-----------------------------------------------------------------------
!
!  At this point tem3 contains
!
!      avgy(v)*difx(mapfct_u) - avgx(u)*dify(mapfct_v)
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!  Calculate avgsu(avgy(v*))*tem3 and subtract from uforce.
!
!-----------------------------------------------------------------------

  DO k=2,nz-2
    DO j=1,ny-1
      DO i=1,nx-1
        tem2(i,j,k) = tem3(i,j,k)*0.5*(vstr(i,j+1,k)+vstr(i,j,k))
      END DO
    END DO
  END DO

  DO k=2,nz-2
    DO j=1,ny-1
      DO i=2,nx-1
        uforce(i,j,k) = uforce(i,j,k) - 0.5*(tem2(i,j,k)+tem2(i-1,j,k))
      END DO
    END DO
  END DO

!-----------------------------------------------------------------------
!
!  Calculate avgsv(avgx(u*))*tem3 and add to vforce.
!
!-----------------------------------------------------------------------

  DO k=2,nz-2
    DO j=1,ny-1
      DO i=1,nx-1
        tem2(i,j,k) = tem3(i,j,k)*0.5*(ustr(i+1,j,k)+ustr(i,j,k))
      END DO
    END DO
  END DO

  DO k=2,nz-2
    DO j=2,ny-1
      DO i=1,nx-1
        vforce(i,j,k) = vforce(i,j,k)+0.5*(tem2(i,j-1,k)+tem2(i,j,k))
      END DO
    END DO
  END DO

  RETURN
END SUBROUTINE maptrm


!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE RHOFACTOR                  ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE rhofactor(nx,ny,nz,ptprt,pprt,qv,qscalar,             &
           ptbar,pbar,ptbari,pbari,qvbar,rhofct,tem1)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate the rho-factor term: rhobar/rho to correct the
!  calculation of the pressure gradient force terms
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Donghai Wang
!  05/26/97.
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical direction.
!
!    ptprt    Perturbation potential temperature at a time level (K)
!    pprt     Perturbation pressure at a given time level (Pascal)
!    qv       Water vapor specific humidity at a given time level
!             (kg/kg)
!    qc       Cloud water mixing ratio at a given time level (kg/kg)
!    qr       Rainwater mixing ratio at a given time level (kg/kg)
!    qi       Cloud ice mixing ratio at a given time level (kg/kg)
!    qs       Snow mixing ratio at a given time level (kg/kg)
!    qh       Hail mixing ratio at a given time level (kg/kg)
!
!    ptbar    Base state potential temperature (K)
!    pbar     Base state pressure (Pascal)
!    ptbari   Inverse Base state pot. temperature (K)
!    pbari    Inverse Base state pressure (Pascal)
!    qvbar    Base state water vapor specific humidity (kg/kg)
!
!  OUTPUT:
!
!    rhofct   the factor at scalar point
!
!  WORK ARRAYS:
!
!    tem1     Temporary work array.
!    rhofct     Temporary work array.
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
  INCLUDE 'globcst.inc'     ! Global model control parameters
  INCLUDE 'phycst.inc'      ! Physical constants
!
!-----------------------------------------------------------------------
!
!  Variable Declarations
!
!-----------------------------------------------------------------------
!
  INTEGER :: nx,ny,nz          ! Number of grid points in 3 directions

  REAL :: ptprt (nx,ny,nz)     ! Perturbation potential temperature
                               ! at a given time level (K)
  REAL :: pprt  (nx,ny,nz)     ! Perturbation pressure at a given time
                               ! level (Pascal)
  REAL :: qv    (nx,ny,nz)     ! Water vapor specific humidity (kg/kg)

  REAL :: qscalar    (nx,ny,nz,nscalar)

  REAL :: ptbar (nx,ny,nz)     ! Base state potential temperature (K)
  REAL :: pbar  (nx,ny,nz)     ! Base state pressure (Pascal).
  REAL :: ptbari(nx,ny,nz)     ! Inverse Base state pot. temperature (K)
  REAL :: pbari (nx,ny,nz)     ! Inverse Base state pressure (Pascal).
  REAL :: qvbar (nx,ny,nz)     ! Base state water vapor specific
                               ! humidity(kg/kg)

  REAL :: rhofct(nx,ny,nz)      !

  REAL :: tem1  (nx,ny,nz)     ! Temporary work array.

!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  REAL    :: pttem,ptem,tema
  INTEGER :: i,j,k
  INTEGER :: nq
!
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
!
!       rhofct = rhobar/rho
!              = 1.0/(1.0-B/g)
!              = 1.0/(1.0-tem1)
!
!  where   B/g = ptprt/ptbar-pprt/(cpdcv*pbar)+            !1st-order
! qvprt/(rddrv+qvbar)-(qvprt+qc+qr+qs+qi+qh)/(1+qvbar(i,j,k)) !1st-order
! -(ptprt*ptprt)/(ptbar*ptbar)-1.0/(2.0*cpdcv)*(1.0/cpdcv-1.0)!2nd-order
! *(pprt*pprt)/(pbar*pbar)+(ptprt*pprt)/(2.0*cpdcv*ptbar*pbar)!2nd-order
!
!  B is the buoyancy term,
!  and rddrv=rd/rv,cpdcv=cp/cv, cp, cv, rd and rv are defined
!  in phycst.inc.
!
!-----------------------------------------------------------------------
!
  tema = 1.0/cpdcv
  DO k=1,nz-1
    DO j=1,ny-1
      DO i=1,nx-1
        pttem = ptprt(i,j,k)*ptbari(i,j,k)
        ptem  = pprt(i,j,k)*tema*pbari(i,j,k)
        tem1(i,j,k) = pttem-pttem*pttem                                 &
                      -ptem-0.5*(1.0-cpdcv)*ptem*ptem                   &
                      +0.5*pttem*ptem
      END DO
    END DO
  END DO

!
!-----------------------------------------------------------------------
!
!  Add on the contributions to the term from the water vapor
!  content and the liquid and ice water loading.
!
!-----------------------------------------------------------------------
!

!  IF( moist == 1 .AND. ice == 0 ) THEN  ! Moist case, no ice.
!
!    DO k=1,nz-1
!      DO j=1,ny-1
!        DO i=1,nx-1
!          tem1(i,j,k) = tem1(i,j,k)                                     &
!              + (qv(i,j,k) - qvbar(i,j,k))/(rddrv + qvbar(i,j,k))       &
!              - (qv(i,j,k) - qvbar(i,j,k) + qc(i,j,k) + qr(i,j,k))      &
!              /(1 + qvbar(i,j,k))
!        END DO
!      END DO
!    END DO
!
!  ELSE IF(moist == 1 .AND. ice == 1) THEN
!    ! Full microphysics case, loading of liquid and ice water.
!
!    DO k=1,nz-1
!      DO j=1,ny-1
!        DO i=1,nx-1
!          tem1(i,j,k) = tem1(i,j,k)                                     &
!              + (qv(i,j,k) - qvbar(i,j,k))/(rddrv + qvbar(i,j,k))       &
!              - (qv(i,j,k) - qvbar(i,j,k) + qc(i,j,k) + qr(i,j,k) +     &
!              qs(i,j,k) + qi(i,j,k) + qh(i,j,k))/(1 + qvbar(i,j,k))
!        END DO
!      END DO
!    END DO
!
!  END IF

  IF (moist == 1) THEN

    DO k=1,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          tem1(i,j,k) = tem1(i,j,k)                                     &
              + (qv(i,j,k) - qvbar(i,j,k))/(rddrv + qvbar(i,j,k))       &
              - (qv(i,j,k) - qvbar(i,j,k))/(1 + qvbar(i,j,k))
        END DO
      END DO
    END DO

    DO k=1,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          DO nq = 1,nscalarq
            tem1(i,j,k) = tem1(i,j,k) - qscalar(i,j,k,nq)/(1+qvbar(i,j,k))
          END DO
        END DO
      END DO
    END DO

  END IF
!
!-----------------------------------------------------------------------
!
!
!    rhofct = 1.0/(1.0-B/g), B is the buoyancy, tem1=B/g
!
!
!-----------------------------------------------------------------------
!
  DO k=1,nz-1
    DO j=1,ny-1
      DO i=1,nx-1
        rhofct(i,j,k)=1.0/(1.0-tem1(i,j,k))
      END DO
    END DO
  END DO

  RETURN
END SUBROUTINE rhofactor
