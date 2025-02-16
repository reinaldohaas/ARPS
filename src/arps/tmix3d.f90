!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE MIXUVW                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE mixuvw(nx,ny,nz, exbcbufsz,                                  &
           u,v,w,ptprt,pprt,qv,qscalar,tke,pbldpth,                     &
           ubar,vbar,ptbar,pbar,rhostr,qvbar,                           &
           usflx,vsflx, x,y,z,zp,mapfct,                                &
           j1,j2,j3,aj3x,aj3y,aj3z,j3inv,ptsfc,                         &
           umix,vmix,wmix,kmh,kmv,rprntl,lenscl,defsq,                  &
           exbcbuf,                                                     &
           tem1,tem2,tem3,tem4,tem5,tem6,                               &
           tem7,tem8,tem9,tem10,tem11,tem12)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate the total mixing (turbulent mixing and the externally
!  imposed computational mixing) for the momentum equations. This
!  subroutine also calculates the turbulent mixing coefficient ,km,
!  which is used to calculate mixing terms for temperature and
!  water quantities. The mixing coefficient is based on Smagorinsky's
!  formulation. For the computational mixing, there are two options:
!  second order and fourth order mixing. The mixing applies only to
!  the perturbations quantities.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  10/10/91
!
!  MODIFICATION HISTORY:
!
!  5/05/92 (M. Xue)
!  Added full documentation.
!
!  6/1/92 (M. Xue and H. Jin)
!  Further facelift.
!
!  7/6/92 (M. Xue and D. Weber)
!  Terrain included.
!
!  1/23/96 (Donghai Wang & Yuhe Liu)
!  Added the map projection factor.
!
!  4/22/96 (M. Xue)
!  Modified the calls to computational mixing routines to pass
!  in rhostr instead of rhobar. Effectively, J3 is included in the
!  formulation as was in ARPS version 3.0.
!
!  11/06/97 (D. Weber)
!  Added three additional levels to the mapfct array.  The three
!  levels (4,5,6) represent the inverse of the first three in order.
!  The inverse map factors are computed to improve efficiency.
!
!  9/18/98 (D. Weber)
!  Modified do loop structure to improve code efficiency via:
!  -merging existing loops together
!  -removing and merging operators
!  -switched to DO  ENDDO structure for f90 conversion
!  -added tem12
!  -added aj3x,y,z
!  -added mapfct(i,j,7) = mapfct(i,j,1)*mapfct(i,j,1)
!  -added mapfct(i,j,8) = 0.25*mapfct(i,j,1)
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
!    w        Vertical component of velocity in Cartesian
!             coordinates at a given time level (m/s)
!    ptprt    Perturbation potential temperature at a given time level (K)
!    pprt     Perturbation pressure at a given time level (Pascal)
!    qv       Water vapor specific humidity at a given time level (kg/kg)
!    qc       Cloud water mixing ratio at a given time level (kg/kg)
!    qr       Rainwater mixing ratio at a given time level (kg/kg)
!    qi       Cloud ice mixing ratio at a given time level (kg/kg)
!    qs       Snow mixing ratio at a given time level (kg/kg)
!    qh       Hail mixing ratio at a given time level (kg/kg)
!    tke      Turbulent Kinetic Energy ((m/s)**2)
!    pbldpth  Planetary boundary layer depth (m)
!
!    ubar     Base state x velocity component (m/s)
!    vbar     Base state y velocity component (m/s)
!    ptbar    Base state potential temperature (K)
!    pbar     Base state pressure (Pascal)
!    rhostr   Base state density rhobar times j3 (kg/m**3)
!    qvbar    Base state water vapor specific humidity (kg/kg)
!
!    usflx    Surface flux of u-momentum
!    vsflx    Surface flux of v-momentum
!
!    x        x coordinate of grid points in physical/comp. space (m)
!    y        y coordinate of grid points in physical/comp. space (m)
!    z        z coordinate of grid points in computational space (m)
!    zp       Vertical coordinate of grid points in physical space (m)
!
!    mapfct   Map factors at scalar, u and v points
!
!    j1       Coordinate transformation Jacobian -d(zp)/dx
!    j2       Coordinate transformation Jacobian -d(zp)/dy
!    j3       Coordinate transformation Jacobian  d(zp)/dz
!    aj3x     Avgx of the coordinate transformation Jacobian  d(zp)/dz
!    aj3y     Avgy of the coordinate transformation Jacobian  d(zp)/dz
!    aj3z     Avgz of the coordinate transformation Jacobian  d(zp)/dz
!    j3inv    Inverse of j3
!
!  OUTPUT:
!
!    umix     Total mixing in u equation (kg/(m*s)**2)
!    vmix     Total mixing in v equation (kg/(m*s)**2)
!    wmix     Total mixing in w equation (kg/(m*s)**2)
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
!   (These arrays are defined and used locally (i.e. inside this
!    subroutine), they may also be passed into routines called by
!    this one. Exiting the call to this subroutine, these temporary
!    work arrays may be used for other purposes and therefore their
!    contents may be overwritten. Please examine the usage of work
!    arrays before you alter the code.)
!
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
  INCLUDE 'timelvls.inc'
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  INTEGER :: nx,ny,nz          ! Number of grid points in 3 directions

  REAL :: u     (nx,ny,nz)     ! Total u-velocity (m/s)
  REAL :: v     (nx,ny,nz)     ! Total v-velocity (m/s)
  REAL :: w     (nx,ny,nz)     ! Total w-velocity (m/s)
  REAL :: ptprt (nx,ny,nz)     ! Perturbation potential temperature (K)
  REAL :: pprt  (nx,ny,nz)     ! Perturbation pressure (Pascal)

  REAL :: qv    (nx,ny,nz)     ! Water vapor specific humidity (kg/kg)

  REAL :: qscalar(nx,ny,nz,nscalar)

  REAL :: tke   (nx,ny,nz,nt)  ! Turbulent Kinetic Energy ((m/s)**2)
  REAL :: pbldpth(nx,ny)       ! Planetary boundary layer depth (m)

  REAL :: ubar  (nx,ny,nz)     ! Base state u-velocity (m/s)
  REAL :: vbar  (nx,ny,nz)     ! Base state v-velocity (m/s)
  REAL :: ptbar (nx,ny,nz)     ! Base state potential temperature (K)
  REAL :: pbar  (nx,ny,nz)     ! Base state pressure (Pascal).
  REAL :: rhostr(nx,ny,nz)     ! Base state density rhobar times j3.
  REAL :: qvbar (nx,ny,nz)     ! Base state water vapor specific humidity
                               ! (kg/kg)

  REAL :: usflx (nx,ny)        ! surface flux of u-momentum (kg/(m*s**2))
  REAL :: vsflx (nx,ny)        ! surface flux of v-momentum (kg/(m*s**2))
  REAL :: x     (nx)           ! The x-coord. of the physical and
                               ! computational grid. Defined at u-point.
  REAL :: y     (ny)           ! The y-coord. of the physical and
                               ! computational grid. Defined at v-point.
  REAL :: z     (nz)           ! The z-coord. of the computational grid.

                                ! Defined at w-point on the staggered grid.
  REAL :: zp    (nx,ny,nz)     ! The physical height coordinate defined at
                               ! w-point of the staggered grid.

  REAL :: mapfct(nx,ny,8)      ! Map factors at scalar, u and v points

  REAL :: j1    (nx,ny,nz)     ! Coordinate transformation Jacobian defined as
                               ! - d( zp )/d( x ).
  REAL :: j2    (nx,ny,nz)     ! Coordinate transformation Jacobian defined as
                               ! - d( zp )/d( y ).
  REAL :: j3    (nx,ny,nz)     ! Coordinate transformation Jacobian defined as
                               ! d( zp )/d( z ).
  REAL :: aj3x  (nx,ny,nz)     ! Coordinate transformation Jacobian defined
                               ! as d( zp )/d( z ) AVERAGED IN THE X-DIR.
  REAL :: aj3y  (nx,ny,nz)     ! Coordinate transformation Jacobian defined
                               ! as d( zp )/d( z ) AVERAGED IN THE Y-DIR.
  REAL :: aj3z  (nx,ny,nz)     ! Coordinate transformation Jacobian defined
                               ! as d( zp )/d( z ) AVERAGED IN THE Z-DIR.

  REAL :: j3inv (nx,ny,nz)     ! Inverse of j3
  REAL :: ptsfc  (nx,ny)       ! Ground surface potential temperature (K)

  REAL :: umix  (nx,ny,nz)     ! Turbulent mixing on u-momentum.
  REAL :: vmix  (nx,ny,nz)     ! Turbulent mixing on v-momentum.
  REAL :: wmix  (nx,ny,nz)     ! Turbulent mixing on w-momentum.

  REAL :: kmh   (nx,ny,nz)     ! Horizontal turb. mixing coef. for
                               ! momentum. ( m**2/s )
  REAL :: kmv   (nx,ny,nz)     ! Vertical turb. mixing coef. for
                               ! momentum. ( m**2/s )
  REAL :: rprntl(nx,ny,nz)     ! Reciprocal of Prandtl number
  REAL :: lenscl(nx,ny,nz)     ! Turbulent mixing length scale (m)
  REAL :: defsq (nx,ny,nz)     ! Deformation squared (1/s**2)

  INTEGER :: exbcbufsz            ! EXBC buffer size
  REAL    :: exbcbuf( exbcbufsz ) ! EXBC buffer array

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
  REAL :: tem12 (nx,ny,nz)     ! Temporary work array
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: nxyz
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
!  Control parameters for subgrid scalar turbulent mixing are passed
!  in globcst.inc, these include:
!
!  tmixopt: Tubulent mixing option
!           = 0, zero turbulent mixing.
!           = 1, constant mixing coefficient.
!           = 2, Smagorinsky mixing coefficient.
!           = 3, Smagorinsky + constant coefficient mixing.
!           = 4, 1.5 order TKE turbulence
!  tmixcst: Constant mixing coefficient. (Options 1 and 3)
!
!-----------------------------------------------------------------------
!
  CALL set_acct(tmix_acct)

  IF( tmixopt == 0 .AND. sfcphy == 0 ) THEN

    nxyz = nx*ny*nz
    CALL flzero(umix,nxyz)
    CALL flzero(vmix,nxyz)
    CALL flzero(wmix,nxyz)

  ELSE

    CALL tmixuvw(nx,ny,nz,                                              &
         u,v,w,ptprt,pprt,qv,qscalar,tke,pbldpth,                       &
         ubar,vbar,ptbar,pbar,rhostr,qvbar,                             &
         usflx,vsflx, x,y,z,zp,mapfct,                                  &
         j1,j2,j3,aj3x,aj3y,aj3z,j3inv,ptsfc,                           &
         umix,vmix,wmix,kmh,kmv,rprntl,lenscl,defsq,                    &
         tem1,tem2,tem3,tem4,tem5,tem6,tem7,                            &
         tem8,tem9,tem10,tem11,tem12)

  END IF

  CALL set_acct(cmix_acct)

  IF( cmix2nd == 1) THEN
!
!-----------------------------------------------------------------------
!
!  Calculate the second order computational mixing terms and add
!  them to the turbulent mixing terms umix, vmix, wmix.
!
!-----------------------------------------------------------------------
!
    CALL cmix2uvw(nx,ny,nz,                                             &
                  u,v,w, ubar,vbar,rhostr,                              &
                  umix,vmix,wmix,                                       &
                  tem1,tem2,tem3)

  END IF

  IF( cmix4th == 1) THEN
!
!-----------------------------------------------------------------------
!
!  Calculate the fourth order computational mixing terms and add
!  them to the turbulent mixing terms umix, vmix, wmix.
!
!-----------------------------------------------------------------------
!
    CALL cmix4uvw(nx,ny,nz,                                             &
                  u,v,w, ubar,vbar,rhostr,                              &
                  umix,vmix,wmix,                                       &
                  tem1,tem2,tem3,tem4,tem5)
  END IF

  IF( raydmp /= 0) THEN
    CALL set_acct(raydmp_acct)
!
!-----------------------------------------------------------------------
!
!    Calculate the upper level Rayleigh damping on u, v and w perturbations.
!    The terms are accumulated in arrays umix, vmix and wmix.
!
!-----------------------------------------------------------------------
!
    CALL rdmpuvw(nx,ny,nz,exbcbufsz,                                    &
                 u,v,w, ubar,vbar,rhostr, zp,                           &
                 umix,vmix,wmix,                                        &
                 exbcbuf,tem1,tem2,tem3,                                &
                 tem4,tem5,tem6,tem7)

  END IF

  RETURN
END SUBROUTINE mixuvw
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE MIXPT                      ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE mixpt(nx,ny,nz, exbcbufsz,                                   &
           ptprt,ptbar,rhostr,kmh,kmv,rprntl,                           &
           usflx,vsflx,ptsflx,pbldpth,                                  &
           x,y,z,zp,mapfct,j1,j2,j3,aj3x,aj3y,j3inv,ptsfc,              &
           ptmix, exbcbuf,                                              &
           tem1,tem2,tem3,tem4,tem5,tem6)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate the total mixing (turbulent mixing and an externally imposed
!  computational mixing) for the potential temperature equation.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  10/10/91
!
!  MODIFICATION HISTORY:
!
!  5/05/92 (M. Xue)
!  Added full documentation.
!
!  6/1/92 (M. Xue and H. Jin)
!  Further facelift.
!
!  7/6/92 (M. Xue and D. Weber)
!  Terrain included.
!
!  1/23/96 (Donghai Wang & Yuhe Liu)
!  Added the map projection factor.
!
!  4/22/96 (M. Xue)
!  Modified the calls to computational mixing routines to pass
!  in rhostr instead of rhobar. Effectively, J3 is included in the
!  formulation as was in ARPS version 3.0.
!
!  9/18/98 (D. Weber)
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
!    ptprt    Perturbation potential temperature at a given time level (K)
!    ptbar    Base state potential temperature (K)
!    rhostr   Base state density rhobar times j3 (kg/m**3)
!
!    kmh      Horizontal turb. mixing coef. for momentum ( m**2/s )
!    kmv      Vertical turb. mixing coef. for momentum ( m**2/s )
!    rprntl   Reciprocal of Prandtl number
!
!    usflx    Surface u-momentum flux
!    vsflx    Surface v-momentum flux
!    ptsflx   Surface heat flux (K*kg/(m**2*s))
!    pbldpth  Planetary boundary layer depth (m)
!
!    x        x coordinate of grid points in physical/comp. space (m)
!    y        y coordinate of grid points in physical/comp. space (m)
!    z        z coordinate of grid points in computational space (m)
!    zp       Vertical coordinate of grid points in physical space (m)
!
!    mapfct   Map factors at scalar points
!
!    j1       Coordinate transformation Jacobian -d(zp)/dx
!    j2       Coordinate transformation Jacobian -d(zp)/dy
!    j3       Coordinate transformation Jacobian  d(zp)/dz
!    aj3x     Avgx of the coordinate transformation Jacobian  d(zp)/dz
!    aj3y     Avgy of the coordinate transformation Jacobian  d(zp)/dz
!    j3inv    Inverse of j3
!    ptsfc    Ground surface potential temperature (K)
!
!  OUTPUT:
!
!    ptmix    Total mixing in potential temperature equation
!             (K*kg/(m**3 * s))
!
!  WORK ARRAYS:
!
!    tem1     Temporary work array.
!    tem2     Temporary work array.
!    tem3     Temporary work array.
!    tem4     Temporary work array.
!    tem5     Temporary work array.
!    tem6     Temporary work array.
!
!   (These arrays are defined and used locally (i.e. inside this
!    subroutine), they may also be passed into routines called by
!    this one. Exiting the call to this subroutine, these temporary
!    work arrays may be used for other purposes and therefore their
!    contents may be overwritten. Please examine the usage of work
!    arrays before you alter the code.)
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
!
  INTEGER :: nx,ny,nz          ! Number of grid points in 3 directions
!
  REAL :: ptprt (nx,ny,nz)     ! Perturbation potential temperature (K)
  REAL :: ptbar (nx,ny,nz)     ! Base state potential temperature (K)
  REAL :: rhostr(nx,ny,nz)     ! Base state density rhobar times j3.

  REAL :: kmh   (nx,ny,nz)     ! Horizontal turb. mixing coef. for
                               ! momentum. ( m**2/s )
  REAL :: kmv   (nx,ny,nz)     ! Vertical turb. mixing coef. for
                               ! momentum. ( m**2/s )
  REAL :: rprntl(nx,ny,nz)     ! Reciprocal of Prandtl number

  REAL :: usflx(nx,ny)         ! surface flux of u-momentum
  REAL :: vsflx(nx,ny)         ! surface flux of v-momentum
  REAL :: ptsflx(nx,ny)        ! surface flux of heat (K*kg/(m**2*s))
  REAL :: pbldpth(nx,ny)       ! Planetary boundary layer depth (m)

  REAL :: x     (nx)           ! The x-coord. of the physical and
                               ! computational grid. Defined at u-point.
  REAL :: y     (ny)           ! The y-coord. of the physical and
                               ! computational grid. Defined at v-point.
  REAL :: z     (nz)           ! The z-coord. of the computational grid.
                               ! Defined at w-point on the staggered grid.
  REAL :: zp    (nx,ny,nz)     ! The physical height coordinate defined at
                               ! w-point of the staggered grid.

  REAL :: mapfct(nx,ny,8)      ! Map factors at scalar points

  REAL :: j1    (nx,ny,nz)     ! Coordinate transformation Jacobian defined as
                               ! - d( zp )/d( x ).
  REAL :: j2    (nx,ny,nz)     ! Coordinate transformation Jacobian defined as
                               ! - d( zp )/d( y ).
  REAL :: j3    (nx,ny,nz)     ! Coordinate transformation Jacobian defined as
  REAL :: aj3x  (nx,ny,nz)     ! Coordinate transformation Jacobian defined
                               ! as d( zp )/d( z ) AVERAGED IN THE X-DIR.
  REAL :: aj3y  (nx,ny,nz)     ! Coordinate transformation Jacobian defined
                               ! as d( zp )/d( z ) AVERAGED IN THE Y-DIR.
  REAL :: j3inv (nx,ny,nz)     ! Inverse of j3
  REAL :: ptsfc  (nx,ny)       ! Ground surface potential temperature (K)

  REAL :: ptmix(nx,ny,nz)      ! Mixing on potential
                               ! temperature (K*kg/(m**3 * s))

  INTEGER :: exbcbufsz         ! EXBC buffer size
  REAL :: exbcbuf( exbcbufsz ) ! EXBC buffer array

  REAL :: tem1  (nx,ny,nz)     ! Temporary work array
  REAL :: tem2  (nx,ny,nz)     ! Temporary work array
  REAL :: tem3  (nx,ny,nz)     ! Temporary work array
  REAL :: tem4  (nx,ny,nz)     ! Temporary work array
  REAL :: tem5  (nx,ny,nz)     ! Temporary work array
  REAL :: tem6  (nx,ny,nz)     ! Temporary work array

!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: nxyz
!
!-----------------------------------------------------------------------
!
!  Include files:
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
  CALL set_acct(tmix_acct)

  IF( tmixopt == 0 .AND. sfcphy == 0 ) THEN

    nxyz = nx*ny*nz
    CALL flzero(ptmix,nxyz)

  ELSE

    CALL tmixpt(nx,ny,nz, ptprt,ptbar,rhostr,kmh,kmv,rprntl,            &
               usflx,vsflx,ptsflx,pbldpth,                              &
               x,y,z,zp,mapfct, j1,j2,j3,aj3x,aj3y,j3inv,ptsfc,         &
               ptmix,                                                   &
               tem1,tem2,tem3,tem4,tem5,tem6)

  END IF

  CALL set_acct(cmix_acct)

  IF( cmix2nd == 1) THEN
!
!-----------------------------------------------------------------------
!
!  Calculate the second order computational mixing term, ptmix, for
!  the potential temperature equation.
!
!-----------------------------------------------------------------------
!
    CALL cmix2s(nx,ny,nz, ptprt,rhostr, ptmix, tem1,tem2,tem3)

  END IF

  IF( cmix4th == 1) THEN
!
!-----------------------------------------------------------------------
!
!  Calculate the fourth order computational mixing term, ptmix, for
!  the potential temperature equation.
!
!-----------------------------------------------------------------------
!

    CALL cmix4s(nx,ny,nz, ptprt,rhostr, ptmix, tem1,tem2,tem3,tem4)

  END IF

  IF( raydmp /= 0) THEN
!
!-----------------------------------------------------------------------
!
!    Calculate the upper level Rayleigh damping on the potential temperature
!    perturbations. The term is accumulated in array ptmix.
!
!-----------------------------------------------------------------------
!
    CALL set_acct(raydmp_acct)

    CALL rdmppt(nx,ny,nz, exbcbufsz,                                    &
                ptprt ,rhostr, zp,                                      &
                ptmix,                                                  &
                exbcbuf,tem1,                                           &
                tem2,tem3)

  END IF


  RETURN
END SUBROUTINE mixpt
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE MIXQV                      ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE mixqv(nx,ny,nz,                                              &
           qv,rhostr,qvbar,kmh,kmv,rprntl,qvsflx,pbldpth,               &
           x,y,z,zp,mapfct,j1,j2,j3,aj3x,aj3y,j3inv,                    &
           usflx,vsflx,ptsflx,ptsfc,qvsfc,ptbar,ptprt,                  &
           qvmix,                                                       &
           tem1,tem2,tem3,tem4,tem5,tem6)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate the total mixing (externally imposed computational mixing
!  and turbulent mixing) for the water vapor mixing ratio equation.
!  This mixing term is different from those of other water quantities
!  since QV has a base state.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  10/10/91
!
!  MODIFICATION HISTORY:
!
!  5/05/92 (M. Xue)
!  Added full documentation.
!
!  6/1/92 (M. Xue and H. Jin)
!  Further facelift.
!
!  7/6/92 (M. Xue and D. Weber)
!  Terrain included.
!
!  1/23/96 (Donghai Wang & Yuhe Liu)
!  Added the map projection factor.
!
!  4/22/96 (M. Xue)
!  Modified the calls to computational mixing routines to pass
!  in rhostr instead of rhobar. Effectively, J3 is included in the
!  formulation as was in ARPS version 3.0.
!
!  9/18/98 (D. Weber)
!  Modified do loop structure to improve code efficiency via:
!  -merging existing loops together
!  -removing and merging operators
!  -switched to DO  ENDDO structure for f90 conversion
!  -added arrays aj3x,y.
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!    qv       Water vapor specific humidity at a given time level (kg/kg)
!    qvbar    Base state water vapor specific humidity (kg/kg)
!    rhostr   Base state density rhobar times j3 (kg/m**3)
!
!    kmh      Horizontal turb. mixing coef. for momentum ( m**2/s )
!    kmv      Vertical turb. mixing coef. for momentum ( m**2/s )
!    rprntl   Reciprocal of Prandtl number
!
!    qvsflx   Surface moisture flux (K*kg/(m**2*s))
!
!    x        x coordinate of grid points in physical/comp. space (m)
!    y        y coordinate of grid points in physical/comp. space (m)
!    z        z coordinate of grid points in computational space (m)
!    zp       Vertical coordinate of grid points in physical space (m)
!
!    mapfct   Map factors at scalar points
!
!    j1       Coordinate transformation Jacobian -d(zp)/dx
!    j2       Coordinate transformation Jacobian -d(zp)/dy
!    j3       Coordinate transformation Jacobian  d(zp)/dz
!    aj3x     Avgx of the coordinate transformation Jacobian  d(zp)/dz
!    aj3y     Avgy of the coordinate transformation Jacobian  d(zp)/dz
!    j3inv    Inverse of j3
!
!  OUTPUT:
!
!    qvmix    Total mixing in water vapor equation (kg/(m**3 * s))
!
!  WORK ARRAYS:
!
!    tem1     Temporary work array.
!    tem2     Temporary work array.
!    tem3     Temporary work array.
!    tem4     Temporary work array.
!    tem5     Temporary work array.
!    tem6     Temporary work array.
!
!   (These arrays are defined and used locally (i.e. inside this
!    subroutine), they may also be passed into routines called by
!    this one. Exiting the call to this subroutine, these temporary
!    work arrays may be used for other purposes and therefore their
!    contents may be overwritten. Please examine the usage of work
!    arrays before you alter the code.)
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
  INTEGER :: nx,ny,nz          ! Number of grid points in 3 directions
!
  REAL :: qv    (nx,ny,nz)     ! Water vapor specific humidity (kg/kg)
  REAL :: qvbar (nx,ny,nz)     ! Base state water vapor specific humidity
                               ! (kg/kg)
  REAL :: rhostr(nx,ny,nz)     ! Base state density rhobar times j3.
  REAL :: kmh   (nx,ny,nz)     ! Horizontal turb. mixing coef. for
                               ! momentum. ( m**2/s )
  REAL :: kmv   (nx,ny,nz)     ! Vertical turb. mixing coef. for
                               ! momentum. ( m**2/s )
  REAL :: rprntl(nx,ny,nz)     ! Reciprocal of Prandtl number

  REAL :: qvsflx(nx,ny)        ! surface flux of moisture (kg/(m**2*s))
  REAL :: pbldpth(nx,ny)       ! Planetary boundary layer depth (m)

  REAL :: x     (nx)           ! The x-coord. of the physical and
                               ! computational grid. Defined at u-point.
  REAL :: y     (ny)           ! The y-coord. of the physical and
                               ! computational grid. Defined at v-point.
  REAL :: z     (nz)           ! The z-coord. of the computational grid.
                               ! Defined at w-point on the staggered grid.
  REAL :: zp    (nx,ny,nz)     ! The physical height coordinate defined at
                               ! w-point of the staggered grid.

  REAL :: mapfct(nx,ny,8)      ! Map factors at scalar points

  REAL :: j1    (nx,ny,nz)     ! Coordinate transformation Jacobian defined as
                               ! - d( zp )/d( x ).
  REAL :: j2    (nx,ny,nz)     ! Coordinate transformation Jacobian defined as
                               ! - d( zp )/d( y ).
  REAL :: j3    (nx,ny,nz)     ! Coordinate transformation Jacobian defined as
                               ! d( zp )/d( z ).
  REAL :: aj3x  (nx,ny,nz)     ! Coordinate transformation Jacobian defined
                               ! as d( zp )/d( z ) AVERAGED IN THE X-DIR.
  REAL :: aj3y  (nx,ny,nz)     ! Coordinate transformation Jacobian defined
                               ! as d( zp )/d( z ) AVERAGED IN THE Y-DIR.
  REAL :: j3inv (nx,ny,nz)     ! Inverse of j3

  REAL :: ptbar (nx,ny,nz)     ! Base state potential temperature (K)
  REAL :: ptprt (nx,ny,nz)     ! Perturbation potential temperature
  REAL :: ptsfc(nx,ny)         ! Temperature at ground (K) (in top 1 cm layer)
  REAL :: qvsfc(nx,ny)         ! Effective qv at the surface (kg/kg)

  REAL :: usflx(nx,ny)         ! Surface u-momentum flux
  REAL :: vsflx(nx,ny)         ! Surface v-momentum flux
  REAL :: ptsflx(nx,ny)        ! Surface heat flux (K*kg/(m**2*s))

!
  REAL :: qvmix(nx,ny,nz)      ! Mixing on water vapor specific humidity
                               ! (kg/(m**3 *s))

  REAL :: tem1  (nx,ny,nz)     ! Temporary work array
  REAL :: tem2  (nx,ny,nz)     ! Temporary work array
  REAL :: tem3  (nx,ny,nz)     ! Temporary work array
  REAL :: tem4  (nx,ny,nz)     ! Temporary work array
  REAL :: tem5  (nx,ny,nz)     ! Temporary work array
  REAL :: tem6  (nx,ny,nz)     ! Temporary work array

!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i, j, k
  INTEGER :: nxyz
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
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
  CALL set_acct(tmix_acct)

  IF( tmixopt == 0 .AND. sfcphy == 0 ) THEN

    nxyz = nx*ny*nz
    CALL flzero(qvmix,nxyz)

  ELSE

    CALL tmixqv(nx,ny,nz, qv,rhostr,kmh,kmv,rprntl,qvsflx,pbldpth,      &
               x,y,z,zp,mapfct, j1,j2,j3,aj3x,aj3y,j3inv,               &
               usflx,vsflx,ptsflx,ptsfc,qvsfc,ptbar,ptprt,              &
               qvmix,                                                   &
               tem1,tem2,tem3,tem4,tem5,tem6)


  END IF

  CALL set_acct(cmix_acct)

  IF( cmix2nd == 1) THEN
!
!-----------------------------------------------------------------------
!
!  Calculate the second order computational mixing term, qvmix,
!  for the qv equation.
!
!-----------------------------------------------------------------------
!
    DO k=1,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          tem5(i,j,k)=qv(i,j,k)-qvbar(i,j,k)
        END DO
      END DO
    END DO

    CALL cmix2s(nx,ny,nz, tem5 ,rhostr, qvmix, tem1,tem2,tem3)

  END IF

  IF( cmix4th == 1) THEN
!
!-----------------------------------------------------------------------
!
!  Calculate the fourth order computational mixing term, qvmix,
!  for qv equation.
!
!-----------------------------------------------------------------------
!
    DO k=1,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          tem5(i,j,k)=qv(i,j,k)-qvbar(i,j,k)
        END DO
      END DO
    END DO

    CALL cmix4s(nx,ny,nz, tem5 ,rhostr, qvmix, tem1,tem2,tem3,tem4)

  END IF

  RETURN
END SUBROUTINE mixqv
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE MIXQ                       ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE mixq(nx,ny,nz,                                               &
           q,rhostr,kmh,kmv,rprntl,x,y,z,zp,mapfct,j1,j2,j3,            &
           aj3x,aj3y,j3inv,                                             &
           qmix,                                                        &
           tem1,tem2,tem3,tem4,tem5,tem6)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate the total mixing (which includes the turbulent mixing
!  as well as externally imposed computational mixing) for all water
!  substance equations except water vapor.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  10/10/91
!
!  MODIFICATION HISTORY:
!
!  5/05/92 (M. Xue)
!  Added full documentation.
!
!  6/1/92 (M. Xue and H. Jin)
!  Further facelift.
!
!  7/6/92 (D. Weber)
!  Terrain included.
!
!  1/23/96 (Donghai Wang & Yuhe Liu)
!  Added the map projection factor.
!
!  4/22/96 (M. Xue)
!  Modified the calls to computational mixing routines to pass
!  in rhostr instead of rhobar. Effectively, J3 is included in the
!  formulation as was in ARPS version 3.0.
!
!  9/18/98 (D. Weber)
!  Added arrays aj3x and aj3y.
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!    q        water/ice mixing ratio (kg/kg)
!    kmh      Horizontal turb. mixing coef. for momentum ( m**2/s )
!    kmv      Vertical turb. mixing coef. for momentum ( m**2/s )
!    rprntl   Reciprocal of Prandtl number
!    rhostr   Base state density rhobar times j3 (kg/m**3)
!
!    x        x coordinate of grid points in physical/comp. space (m)
!    y        y coordinate of grid points in physical/comp. space (m)
!    z        z coordinate of grid points in computational space (m)
!    zp       Vertical coordinate of grid points in physical space (m)
!
!    mapfct   Map factors at scalar points
!
!    j1       Coordinate transformation Jacobian -d(zp)/dx
!    j2       Coordinate transformation Jacobian -d(zp)/dy
!    j3       Coordinate transformation Jacobian  d(zp)/dz
!    aj3x     Avgx of the coordinate transformation Jacobian  d(zp)/dz
!    aj3y     Avgy of the coordinate transformation Jacobian  d(zp)/dz
!    j3inv    Inverse of j3
!
!  OUTPUT:
!
!    qmix     Total mixing in water/ice substance
!             equation (kg/(m**3 * s))
!
!  WORK ARRAYS:
!
!    tem1     Temporary work array.
!    tem2     Temporary work array.
!    tem3     Temporary work array.
!    tem4     Temporary work array.
!    tem5     Temporary work array.
!    tem6     Temporary work array.
!
!   (These arrays are defined and used locally (i.e. inside this
!    subroutine), they may also be passed into routines called by
!    this one. Exiting the call to this subroutine, these temporary
!    work arrays may be used for other purposes and therefore their
!    contents may be overwritten. Please examine the usage of work
!    arrays before you alter the code.)
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
  INTEGER :: nx,ny,nz          ! Number of grid points in 3 directions
!
  REAL :: q     (nx,ny,nz)     ! Water/ice quantity (kg/kg)
  REAL :: kmh   (nx,ny,nz)     ! Horizontal turb. mixing coef. for
                               ! momentum. ( m**2/s )
  REAL :: kmv   (nx,ny,nz)     ! Vertical turb. mixing coef. for
                               ! momentum. ( m**2/s )
  REAL :: rprntl(nx,ny,nz)     ! Reciprocal of Prandtl number

  REAL :: rhostr(nx,ny,nz)     ! Base state density rhobar times j3.

  REAL :: x     (nx)           ! The x-coord. of the physical and
                               ! computational grid. Defined at u-point.
  REAL :: y     (ny)           ! The y-coord. of the physical and
                               ! computational grid. Defined at v-point.
  REAL :: z     (nz)           ! The z-coord. of the computational grid.
                               ! Defined at w-point on the staggered grid.
  REAL :: zp    (nx,ny,nz)     ! The physical height coordinate defined at
                               ! w-point of the staggered grid.
  REAL :: mapfct(nx,ny,8)      ! Map factors at scalar points

  REAL :: j1    (nx,ny,nz)     ! Coordinate transformation Jacobian defined as
                               ! - d( zp )/d( x ).
  REAL :: j2    (nx,ny,nz)     ! Coordinate transformation Jacobian defined as
                               ! - d( zp )/d( y ).
  REAL :: j3    (nx,ny,nz)     ! Coordinate transformation Jacobian defined as
                               ! d( zp )/d( z ).
  REAL :: aj3x  (nx,ny,nz)     ! Coordinate transformation Jacobian defined
                               ! as d( zp )/d( z ) AVERAGED IN THE X-DIR.
  REAL :: aj3y  (nx,ny,nz)     ! Coordinate transformation Jacobian defined
                               ! as d( zp )/d( z ) AVERAGED IN THE Y-DIR.
  REAL :: j3inv (nx,ny,nz)     ! Inverse of j3
!
  REAL :: qmix(nx,ny,nz)       ! Water/ice mixing
!
  REAL :: tem1  (nx,ny,nz)     ! Temporary work array
  REAL :: tem2  (nx,ny,nz)     ! Temporary work array
  REAL :: tem3  (nx,ny,nz)     ! Temporary work array
  REAL :: tem4  (nx,ny,nz)     ! Temporary work array
  REAL :: tem5  (nx,ny,nz)     ! Temporary work array
  REAL :: tem6  (nx,ny,nz)     ! Temporary work array

!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!

!
!-----------------------------------------------------------------------
!
!  Include files:
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
  CALL set_acct(tmix_acct)

  IF( tmixopt == 0) THEN

    qmix = 0.0

  ELSE

    CALL tmixq(nx,ny,nz, q,rhostr,kmh,kmv,rprntl,                       &
               x,y,z,zp, mapfct, j1,j2,j3,aj3x,aj3y,j3inv,              &
               qmix,                                                    &
               tem1,tem2,tem3,tem4,tem5,tem6)

  END IF

  CALL set_acct(cmix_acct)

  IF( cmix2nd == 1) THEN
!
!-----------------------------------------------------------------------
!
!  Calculate the second order computational mixing term, qmix,
!  for the q equation.
!
!-----------------------------------------------------------------------
!
    CALL cmix2s(nx,ny,nz, q ,rhostr, qmix, tem1,tem2,tem3)

  END IF

  IF( cmix4th == 1) THEN
!
!-----------------------------------------------------------------------
!
!  Calculate the fourth order computational mixing term, qmix,
!  for the q equation.
!
!-----------------------------------------------------------------------

    CALL cmix4s(nx,ny,nz, q ,rhostr, qmix, tem1,tem2,tem3,tem4)

  END IF

  RETURN
END SUBROUTINE mixq
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE MIXTKE                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE mixtke(nx,ny,nz,                                             &
           tke,rhostr,kmh,kmv,x,y,z,zp,mapfct,                          &
           j1,j2,j3,aj3x,aj3y,j3inv,                                    &
           tkemix, tem1,tem2,tem3,tem4,tem5,tem6,tem7)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate the mixing term in TKE equation.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  8/15/95
!
!  MODIFICATION HISTORY:
!
!  1/23/96 (Donghai Wang & Yuhe Liu)
!  Added the map projection factor.
!
!  4/22/96 (M. Xue)
!  Modified the calls to computational mixing routines to pass
!  in rhostr instead of rhobar. Effectively, J3 is included in the
!  formulation as was in ARPS version 3.0.
!
!  07/10/97 (Fanyou Kong - CMRP)
!  Change the upper limit in DO LOOP 210 to nx-1, ny-1, and nz-1
!  to avoid possible floating overflow
!
!  8/27/98 (D. Weber)
!  Modified do loop structure to improve code efficiency via:
!  -merging existing loops together
!  -removing and merging operators
!  -switched to DO  ENDDO structure for f90 conversion
!  -added aj3x,y arrays
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!    tke      Turbulent kinetic energy
!    kmh      Horizontal turb. mixing coef. for momentum ( m**2/s )
!    kmv      Vertical turb. mixing coef. for momentum ( m**2/s )
!    rhostr   Base state density rhobar times j3 (kg/m**3)
!
!    x        x coordinate of grid points in physical/comp. space (m)
!    y        y coordinate of grid points in physical/comp. space (m)
!    z        z coordinate of grid points in computational space (m)
!    zp       Vertical coordinate of grid points in physical space (m)
!
!    mapfct   Map factors at scalar points
!
!    j1       Coordinate transformation Jacobian -d(zp)/dx
!    j2       Coordinate transformation Jacobian -d(zp)/dy
!    j3       Coordinate transformation Jacobian  d(zp)/dz
!    aj3x     Avgx of the coordinate transformation Jacobian  d(zp)/dz
!    aj3y     Avgy of the coordinate transformation Jacobian  d(zp)/dz
!    j3inv    Inverse of j3
!
!  OUTPUT:
!
!    tkemix   Mixing term in TKE equation
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
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
  INTEGER :: nx,ny,nz          ! Number of grid points in 3 directions
!
  REAL :: tke   (nx,ny,nz)     ! Turbulent kinetic energy
  REAL :: kmh   (nx,ny,nz)     ! Horizontal turb. mixing coef. for
                               ! momentum. ( m**2/s )
  REAL :: kmv   (nx,ny,nz)     ! Vertical turb. mixing coef. for
                               ! momentum. ( m**2/s )
  REAL :: rhostr(nx,ny,nz)     ! Base state density rhobar times j3.

  REAL :: x     (nx)           ! The x-coord. of the physical and
                               ! computational grid. Defined at u-point.
  REAL :: y     (ny)           ! The y-coord. of the physical and
                               ! computational grid. Defined at v-point.
  REAL :: z     (nz)           ! The z-coord. of the computational grid.
                               ! Defined at w-point on the staggered grid.
  REAL :: zp    (nx,ny,nz)     ! The physical height coordinate defined at
                               ! w-point of the staggered grid.
  REAL :: mapfct(nx,ny,8)      ! Map factors at scalar points

  REAL :: j1    (nx,ny,nz)     ! Coordinate transformation Jacobian defined as
                               ! - d( zp )/d( x ).
  REAL :: j2    (nx,ny,nz)     ! Coordinate transformation Jacobian defined as
                               ! - d( zp )/d( y ).
  REAL :: j3    (nx,ny,nz)     ! Coordinate transformation Jacobian defined as
                               ! d( zp )/d( z ).
  REAL :: aj3x  (nx,ny,nz)     ! Coordinate transformation Jacobian defined
                               ! as d( zp )/d( z ) AVERAGED IN THE X-DIR.
  REAL :: aj3y  (nx,ny,nz)     ! Coordinate transformation Jacobian defined
                               ! as d( zp )/d( z ) AVERAGED IN THE Y-DIR.
  REAL :: j3inv (nx,ny,nz)     ! Inverse of j3

  REAL :: tkemix(nx,ny,nz)     ! Mixing term in TKE equation.

  REAL :: tem1  (nx,ny,nz)     ! Temporary work array
  REAL :: tem2  (nx,ny,nz)     ! Temporary work array
  REAL :: tem3  (nx,ny,nz)     ! Temporary work array
  REAL :: tem4  (nx,ny,nz)     ! Temporary work array
  REAL :: tem5  (nx,ny,nz)     ! Temporary work array
  REAL :: tem6  (nx,ny,nz)     ! Temporary work array
  REAL :: tem7  (nx,ny,nz)     ! Temporary work array

!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: nxyz,i,j,k
!
!-----------------------------------------------------------------------
!
!  Include files:
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
  CALL set_acct(tmix_acct)

  IF( tmixopt == 0) THEN

    nxyz = nx*ny*nz
    CALL flzero(tkemix,nxyz)

  ELSE

    DO k=1,nz
      DO j=1,ny
        DO i=1,nx
          tem7(i,j,k)=1.0
        END DO
      END DO
    END DO

    CALL tmixq(nx,ny,nz,tke,rhostr,kmh,kmv, tem7,                       &
               x,y,z,zp, mapfct, j1,j2,j3,aj3x,aj3y,j3inv,              &
               tkemix,                                                  &
               tem1,tem2,tem3,tem4,tem5,tem6)


    DO k=1,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          tkemix(i,j,k)= 2.0*tkemix(i,j,k)
        END DO
      END DO
    END DO

  END IF

  CALL set_acct(cmix_acct)

  IF( cmix2nd == 1) THEN
!
!-----------------------------------------------------------------------
!
!  Calculate the second order computational mixing term, tkemix,
!  for the TKE equation.
!
!-----------------------------------------------------------------------
!
    CALL cmix2s(nx,ny,nz, tke ,rhostr, tkemix, tem1,tem2,tem3)

  END IF

  IF( cmix4th == 1) THEN
!
!-----------------------------------------------------------------------
!
!  Calculate the fourth order computational mixing term, tkemix,
!  for the tke equation.
!
!-----------------------------------------------------------------------
!
    CALL cmix4s(nx,ny,nz, tke ,rhostr, tkemix, tem1,tem2,tem3,tem4)

  END IF

  RETURN
END SUBROUTINE mixtke
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE TMIXUVW                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE tmixuvw(nx,ny,nz,                                            &
           u,v,w,ptprt,pprt,qv,qscalar,tke,pbldpth,                     &
           ubar,vbar,ptbar,pbar,rhostr,qvbar,                           &
           usflx,vsflx, x,y,z,zp, mapfct,                               &
           j1,j2,j3,aj3x,aj3y,aj3z,j3inv,ptsfc,                         &
           umix,vmix,wmix,kmh,kmv,rprntl,lenscl,defsq,                  &
           nsqed,tau11,tau12,tau13,tau22,tau23,tau33,                   &
           tem1,tem2,tem3,tem4,tem5)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate the turbulent mixing terms for the momentum equations. These
!  terms are expressed in the form of a stress tensor.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  10/10/91
!
!  MODIFICATION HISTORY:
!
!  5/05/92 (M. Xue)
!  Added full documentation.
!
!  6/1/92 (M. Xue and H. Jin)
!  Further facelift.
!
!  6/25/92 (M. Xue and D. Weber)
!  Terrain included.
!
!  1/23/96 (Donghai Wang & Yuhe Liu)
!  Added the map projection factor.
!
!  3/16/96 (Ming Xue)
!  The model can now include the surface fluxes with the turbulent
!  mxing option off (tmixopt=0).
!
!  6/5/96 (Donghai Wang, M. Xue and X. Song)
!  Fixed a minor bug related to the vertical implicit treatment of
!  turbulence mixing. Should not affect the results much.
!
!  11/06/97 (D. Weber)
!  Added three additional levels to the mapfct array.  The three
!  levels (4,5,6) represent the inverse of the first three in order.
!  The inverse map factors are computed to improve efficiency.
!
!  9/18/98 (D. Weber)
!  Modified do loop structure to improve code efficiency via:
!  -merging existing loops together
!  -removing and merging operators
!  -switched to DO  ENDDO structure for f90 conversion
!  -added tem5
!  -added aj3x,y,z
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
!    w        Vertical component of velocity in Cartesian
!             coordinates at a given time level (m/s)
!    ptprt    Perturbation potential temperature at a given time level (K)
!    pprt     Perturbation pressure at a given time level (Pascal)
!    qv       Water vapor specific humidity at a given time level (kg/kg)
!    qc       Cloud water mixing ratio at a given time level (kg/kg)
!    qr       Rainwater mixing ratio at a given time level (kg/kg)
!    qi       Cloud ice mixing ratio at a given time level (kg/kg)
!    qs       Snow mixing ratio at a given time level (kg/kg)
!    qh       Hail mixing ratio at a given time level (kg/kg)
!    tke      Turbulent Kinetic Energy ((m/s)**2)
!    pbldpth  Planetary boundary layer depth (m)
!
!    ubar     Base state x velocity component (m/s)
!    vbar     Base state y velocity component (m/s)
!    ptbar    Base state potential temperature (K)
!    pbar     Base state pressure (Pascal)
!    rhostr   Base state density rhobar times j3 (kg/m**3)
!    qvbar    Base state water vapor specific humidity (kg/kg)
!
!    usflx    Surface flux of u-momentum
!    vsflx    Surface flux of v-momentum
!
!    x        x coordinate of grid points in physical/comp. space (m)
!    y        y coordinate of grid points in physical/comp. space (m)
!    z        z coordinate of grid points in computational space (m)
!    zp       Vertical coordinate of grid points in physical space (m)
!
!    mapfct   Map factors at scalar, u and v points
!
!    j1       Coordinate transformation Jacobian -d(zp)/dx
!    j2       Coordinate transformation Jacobian -d(zp)/dy
!    j3       Coordinate transformation Jacobian  d(zp)/dz
!    aj3x     Avgx of the coordinate transformation Jacobian  d(zp)/dz
!    aj3y     Avgy of the coordinate transformation Jacobian  d(zp)/dz
!    aj3z     Avgz of the coordinate transformation Jacobian  d(zp)/dz
!    j3inv    Inverse of j3
!
!  OUTPUT:
!
!    umix     Total mixing in u equation (kg/(m*s)**2)
!    vmix     Total mixing in v equation (kg/(m*s)**2)
!    wmix     Total mixing in w equation (kg/(m*s)**2)
!
!    kmh      Horizontal turb. mixing coef. for momentum ( m**2/s )
!    kmv      Vertical turb. mixing coef. for momentum ( m**2/s )
!    rprntl   Reciprocal of Prandtl number
!    lenscl   Turbulent mixing length scale (m)
!    defsq    Deformation squared (1/s**2)
!
!  WORK ARRAYS:
!
!    real nsqed  Temporary array containing the Brunt Vaisala frquency
!    real tau11  Temporary work array
!    real tau12  Temporary work array
!    real tau13  Temporary work array
!    real tau22  Temporary work array
!    real tau23  Temporary work array
!    real tau33  Temporary work array
!    real tem1   Temporary work array
!    real tem2   Temporary work array
!    real tem3   Temporary work array
!    real tem4   Temporary work array
!    real tem5   Temporary work array
!
!   (These arrays are defined and used locally (i.e. inside this
!    subroutine), they may also be passed into routines called by
!    this one. Exiting the call to this subroutine, these temporary
!    work arrays may be used for other purposes and therefore their
!    contents may be overwritten. Please examine the usage of work
!    arrays before you alter the code.)
!
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
  INCLUDE 'grid.inc'          ! Grid & map parameters.
  INCLUDE 'phycst.inc'

  INCLUDE 'timelvls.inc'
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  INTEGER :: nx,ny,nz          ! Number of grid points in 3 directions

  REAL :: u     (nx,ny,nz)     ! Total u-velocity (m/s)
  REAL :: v     (nx,ny,nz)     ! Total v-velocity (m/s)
  REAL :: w     (nx,ny,nz)     ! Total w-velocity (m/s)
  REAL :: ptprt (nx,ny,nz)     ! Perturbation potential temperature (K)
  REAL :: pprt  (nx,ny,nz)     ! Perturbation pressure (Pascal)
!
  REAL :: qv    (nx,ny,nz)     ! Water vapor specific humidity (kg/kg)

  REAL :: qscalar    (nx,ny,nz,nscalar)

  REAL :: tke   (nx,ny,nz,nt)  ! Turbulent Kinetic Energy ((m/s)**2)
  REAL :: pbldpth(nx,ny)       ! Planetary boundary layer depth (m)
!
  REAL :: ubar  (nx,ny,nz)     ! Base state u-velocity (m/s)
  REAL :: vbar  (nx,ny,nz)     ! Base state v-velocity (m/s)
  REAL :: ptbar (nx,ny,nz)     ! Base state potential temperature (K)
  REAL :: pbar  (nx,ny,nz)     ! Base state pressure (Pascal).
  REAL :: rhostr(nx,ny,nz)     ! Base state density rhobar times j3.
  REAL :: qvbar (nx,ny,nz)     ! Base state water vapor specific humidity
                               ! (kg/kg)

  REAL :: usflx (nx,ny)        ! surface flux of u-momentum (kg/(m*s**2))
  REAL :: vsflx (nx,ny)        ! surface flux of v-momentum (kg/(m*s**2))
  REAL :: x     (nx)           ! The x-coord. of the physical and
                               ! computational grid. Defined at u-point.
  REAL :: y     (ny)           ! The y-coord. of the physical and
                               ! computational grid. Defined at v-point.
  REAL :: z     (nz)           ! The z-coord. of the computational grid.
                               ! Defined at w-point on the staggered grid.
  REAL :: zp    (nx,ny,nz)     ! The physical height coordinate defined at
                               ! w-point of the staggered grid.

  REAL :: mapfct(nx,ny,8)      ! Map factors at scalar, u and v points

  REAL :: j1    (nx,ny,nz)     ! Coordinate transformation Jacobian defined as
                               ! - d( zp )/d( x ).
  REAL :: j2    (nx,ny,nz)     ! Coordinate transformation Jacobian defined as
                               ! - d( zp )/d( y ).
  REAL :: j3    (nx,ny,nz)     ! Coordinate transformation Jacobian defined as
                               ! d( zp )/d( z ).
  REAL :: aj3x  (nx,ny,nz)     ! Coordinate transformation Jacobian defined
                               ! as d( zp )/d( z ) AVERAGED IN THE X-DIR.
  REAL :: aj3y  (nx,ny,nz)     ! Coordinate transformation Jacobian defined
                               ! as d( zp )/d( z ) AVERAGED IN THE Y-DIR.
  REAL :: aj3z  (nx,ny,nz)     ! Coordinate transformation Jacobian defined
                               ! as d( zp )/d( z ) AVERAGED IN THE Z-DIR.
  REAL :: j3inv (nx,ny,nz)     ! Inverse of j3

  REAL :: ptsfc  (nx,ny)       ! Ground surface potential temperature (K)
!
  REAL :: umix  (nx,ny,nz)     ! Turbulent mixing on u-momentum.
  REAL :: vmix  (nx,ny,nz)     ! Turbulent mixing on v-momentum.
  REAL :: wmix  (nx,ny,nz)     ! Turbulent mixing on w-momentum.

  REAL :: kmh   (nx,ny,nz)     ! Horizontal turb. mixing coef. for
                               ! momentum. ( m**2/s )
  REAL :: kmv   (nx,ny,nz)     ! Vertical turb. mixing coef. for
                               ! momentum. ( m**2/s )
  REAL :: rprntl(nx,ny,nz)     ! Reciprocal of Prandtl number
  REAL :: lenscl(nx,ny,nz)     ! Turbulent mixing length scale (m)
  REAL :: defsq (nx,ny,nz)     ! Deformation squared (1/s**2)
!
  REAL :: nsqed (nx,ny,nz)     ! Temporary work array
  REAL :: tau11 (nx,ny,nz)     ! Temporary work array
  REAL :: tau12 (nx,ny,nz)     ! Temporary work array
  REAL :: tau13 (nx,ny,nz)     ! Temporary work array
  REAL :: tau22 (nx,ny,nz)     ! Temporary work array
  REAL :: tau23 (nx,ny,nz)     ! Temporary work array
  REAL :: tau33 (nx,ny,nz)     ! Temporary work array
  REAL :: tem1  (nx,ny,nz)     ! Temporary work array
  REAL :: tem2  (nx,ny,nz)     ! Temporary work array
  REAL :: tem3  (nx,ny,nz)     ! Temporary work array
  REAL :: tem4  (nx,ny,nz)     ! Temporary work array
  REAL :: tem5  (nx,ny,nz)     ! Temporary work array
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i, j, k
  REAL    :: zpup,zplow
  INTEGER :: cdiszero
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  IF( tmixopt == 0 )  THEN

    DO k=1,nz
      DO j=1,ny
        DO i=1,nx
          wmix(i,j,k)=0.0
          tau13(i,j,k)=0.0
          tau23(i,j,k)=0.0
          kmh(i,j,k)=0.0
          kmv(i,j,k)=0.0
        END DO
      END DO
    END DO

  ELSE
!
!-----------------------------------------------------------------------
!
!  Calculate the static stability parameter N squared (Brunt-Vaisala
!  frequency). The arrays tem1,tem2,tau11,tau12 are used as work
!  arrays for the following call.
!
!-----------------------------------------------------------------------
!
    CALL stabnsq(nx,ny,nz,                                              &
                 ptprt,pprt,qv,qscalar,ptbar,qvbar,pbar,                &
                 j3,j3inv,nsqed,tem1,tem2,tau11,tau12)
!
!-----------------------------------------------------------------------
!
!  Calculate the deformation tensor components Dij, which are stored
!  in arrays Tauij.
!
!  Please note umix and vmix are used here to temporarily store
!  d31 and d32.
!
!-----------------------------------------------------------------------
!
    CALL deform(nx,ny,nz,u,v,w,mapfct,j1,j2,j3,aj3x,aj3y,aj3z,j3inv,    &
                tau11,tau12,tau13,tau22,tau23,tau33,umix,vmix,defsq,    &
                tem1,tem2,tem3,tem4,tem5)
!
!-----------------------------------------------------------------------
!
!  Calculate the turbulent mixing coefficient, km, using the modified
!  Smagorinsky and TKE formulations.
!
!-----------------------------------------------------------------------
!
!WRITE(*,*) '==3==',MAXVAL(kmv)

    CALL cftmix(nx,ny,nz,                                               &
                nsqed,zp,tke,pbldpth,defsq, kmh,kmv,rprntl,lenscl,      &
                tem1,tem2)
!WRITE(*,*) '==4==',MAXVAL(kmv)

    !IF (tmixopt >= 5) THEN
    !  ! keep kmv and rprntl not changed, but set kmh inside cftmix
    !END IF
!
!-----------------------------------------------------------------------
!
!  Calculate the stress tensor Tauij = km * (rhostr/j3) * Dij.
!
!  Please note umix and vmix are used here to temporarily store
!  tau31 and tau32, which are used immediatedly by wmixtrm.
!  umix and vmix are freed afterwards.
!
!
!-----------------------------------------------------------------------
!
    CALL stress(nx,ny,nz,j3,j3inv,                                      &
                kmh,kmv,rhostr,tau11,tau12,tau13,tau22,tau23,tau33,     &
                umix,vmix,                                              &
                tem1,tem2,tem3,tem4,tem5)

    CALL wmixtrm(nx,ny,nz, mapfct(1,1,1), j1,j2,j3,aj3x,aj3y,           &
                 umix,vmix,tau33,                                       &
                 wmix, tem1, tem2)

  END IF
!
!-----------------------------------------------------------------------
!
!  Set the surface momentum fluxes:
!
!  When sfcphy = 0, the internally calculated (from SGS turbulence
!  parameterization) surface fluxes are used.
!  Otherwise, pre-calculated fluxes usflx and vsflx are used.
!
!-----------------------------------------------------------------------
!

  IF( (landwtr == 0 .AND. cdmlnd == 0.0) .OR.                           &
        (landwtr == 1 .AND. cdmlnd == 0.0 .AND. cdmwtr == 0.0) ) THEN
    cdiszero = 1
  ELSE
    cdiszero = 0
  END IF

  IF( (sfcphy == 1 .OR. sfcphy == 3) .AND. cdiszero == 1) GO TO 111

  IF( sfcphy /= 0 ) THEN
    IF ( (sflxdis == 0) .OR. (sflxdis == 1 .AND.                        &
           (sfcphy == 1 .OR. sfcphy == 3).AND.cdiszero == 1) .OR.       &
           (sflxdis == 2) .OR. (sflxdis == 3) ) THEN
!
!-----------------------------------------------------------------------
!
!  tau13 at the surface = usflx.
!
!-----------------------------------------------------------------------
!
      DO j=2,ny-2
        DO i=2,nx-1
          tau13(i,j,2) = usflx(i,j)
        END DO
      END DO
!
!-----------------------------------------------------------------------
!
!  tau23 at the surface = vsflx.
!
!-----------------------------------------------------------------------
!
      DO j=2,ny-1
        DO i=2,nx-2
          tau23(i,j,2) = vsflx(i,j)
        END DO
      END DO

    ELSE IF (sflxdis == 1)  THEN

      DO j=2,ny-2
        DO i=2,nx-1

          tau13(i,j,2) = usflx(i,j)
          zplow=0.5*(zp(i,j,2)+zp(i-1,j,2))

          IF(ptsfc(i,j)+ptsfc(i-1,j) > ptprt(i,j,2)                     &
                +ptbar(i,j,2)+ptprt(i-1,j,2)+ptbar(i-1,j,2))THEN

            DO k=3,nz-1
              zpup =0.5*(zp(i,j,k)+zp(i-1,j,k))
              IF ( (zpup-zplow) <= pbldpth(i,j) ) tau13(i,j,k)=usflx(i,j)* &
                  (1.0-(zpup-zplow)/pbldpth(i,j))**2
            END DO

          END IF
        END DO
      END DO

      DO j=2,ny-1
        DO i=2,nx-2

          tau23(i,j,2) = vsflx(i,j)
          zplow=0.5*(zp(i,j,2)+zp(i,j-1,2))

          IF(ptsfc(i,j)+ptsfc(i,j-1) > ptprt(i,j,2)+ptbar(i,j,2)        &
                +ptprt(i,j-1,2)+ptbar(i,j-1,2))THEN

            DO k=3,nz-1
              zpup =0.5*(zp(i,j,k)+zp(i,j-1,k))
              IF ( (zpup-zplow) <= pbldpth(i,j) ) tau23(i,j,k)=vsflx(i,j)* &
                  (1.0-(zpup-zplow)/pbldpth(i,j))**2
            END DO

          END IF
        END DO
      END DO

    END IF

  END IF

  111   CONTINUE

  IF(tmixopt == 0) THEN

    DO k=2,nz-2
      DO j=1,ny-1
        DO i=2,nx-1
          umix(i,j,k)=(tau13(i,j,k+1)-tau13(i,j,k))*dzinv
        END DO
      END DO
    END DO

    DO k=2,nz-2
      DO j=2,ny-1
        DO i=1,nx-1
          vmix(i,j,k)=(tau23(i,j,k+1)-tau23(i,j,k))*dzinv
        END DO
      END DO
    END DO

  ELSE
!
!-----------------------------------------------------------------------
!
!  Calculate the divergence of the stresses. This is the turbulent
!  mixing on u and v momentum (umix and vmix).
!
!-----------------------------------------------------------------------
!
    CALL umixtrm(nx,ny,nz, mapfct(1,1,2),j1,j2,j3,aj3y,                 &
                 tau11,tau12,tau13,                                     &
                 umix, tem1, tem2)

    CALL vmixtrm(nx,ny,nz, mapfct(1,1,3),j1,j2,j3,aj3x,                 &
                 tau12,tau22,tau23,                                     &
                 vmix, tem1, tem2)

  END IF
!
!-----------------------------------------------------------------------
!
!  Set the normal gradient of the mixing terms on the lateral
!  boundaries to zero
!
!  The umix and vmix terms are set for completness, they are not
!  used in any calculations.
!
!-----------------------------------------------------------------------

!  DO 51 k=2,nz-2
!  DO 51 j=1,ny-1
!    umix(1,j,k)=umix(2,j,k)
!    umix(nx,j,k)=umix(nx-1,j,k)
! 51    CONTINUE

!  DO 52 k=2,nz-2
!  DO 52 i=1,nx-1
!    vmix(i,1,k)=vmix(i,2,k)
!    vmix(i,ny,k)=vmix(i,ny-1,k)
!  52    CONTINUE


!
!-----------------------------------------------------------------------
!
!  Set the mixing terms at the lateral boundaries equal to those
!  at the neighboring interior points.  This is done to essure
!  the boundary points get a similar amount of mixing as the interior
!  points. The boundary values will be used only in the case of
!  radiation lateral boundary conditions.
!
!-----------------------------------------------------------------------
!
  DO k=1,nz-1
    DO i=1,nx
      umix(i,   1,k)=umix(i,   2,k)
      umix(i,ny-1,k)=umix(i,ny-2,k)
    END DO
  END DO

  DO k=1,nz-1
    DO j=1,ny-1
      umix( 1,j,k)=umix(   2,j,k)
      umix(nx,j,k)=umix(nx-1,j,k)
    END DO
  END DO

  DO k=1,nz-1
    DO i=1,nx-1
      vmix(i, 1,k)=vmix(i,   2,k)
      vmix(i,ny,k)=vmix(i,ny-1,k)
    END DO
  END DO

  DO k=1,nz-1
    DO j=1,ny
      vmix(   1,j,k)=vmix(   2,j,k)
      vmix(nx-1,j,k)=vmix(nx-2,j,k)
    END DO
  END DO

  DO k=1,nz-1
    DO i=1,nx-1
      wmix(i,   1,k)=wmix(i,   2,k)
      wmix(i,ny-1,k)=wmix(i,ny-2,k)
    END DO
  END DO

  DO k=1,nz-1
    DO j=1,ny-1
      wmix(   1,j,k)=wmix(   2,j,k)
      wmix(nx-1,j,k)=wmix(nx-2,j,k)
    END DO
  END DO

  RETURN
END SUBROUTINE tmixuvw
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE STABNSQ                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE stabnsq(nx,ny,nz,                                            &
           ptprt,pprt,qv,qscalar,ptbar,qvbar,pbar,j3,j3inv,             &
           nsqed, tmprtr,qvs,qw,pt)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate the static stability parameter N-squared (Brunt-Vaisala
!  frequency squared).
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  9/15/1992
!
!  MODIFICATION HISTORY:
!
!  3/25/94 (G. Bassett)
!  The condition for using saturated or nonsaturated formulation for
!  nsqed is changed from qc > 0 to qc >= 1e-06 in order to avoid
!  fluctuations in Km when very small values of qc are present.
!
!  6/7/94 (M. Xue)
!  Correction to the calculation of nsqed for the dry area inside
!  loop 3100. The error was introduced when modifications were made
!  on 3/25/94. Variable pt in that line was mistakenly written as ptbar.
!
!  3/4/96 (M. Xue)
!  Restructured so that no moisture related calculation is done
!  when moist=0.
!
!  3/14/96 (M. Xue)
!  Corrected a bug introduced on 3/4/96. dzinvd2 was not set for
!  moist=0 case.
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!    ptprt    Perturbation potential temperature (K)
!    pprt     Perturbation pressure (Pascal)
!    ptbar    Base state potential temperature (K)
!    qv       Water vapor specific humidity (kg/kg)
!    qc       Cloud water mixing ratio (kg/kg)
!    qr       Rainwater mixing ratio (kg/kg)
!    qi       Cloud ice mixing ratio (kg/kg)
!    qs       Snow mixing ratio (kg/kg)
!    qh       Hail mixing ratio (kg/kg)
!    qvbar    Base state water vapor specific humidity (kg/kg)
!    pbar     Base state pressure (Pascal)
!    j3       Coordinate transformation Jacobian  d(zp)/dz
!    j3inv    Inverse of j3
!
!  OUTPUT:
!
!    nsqed    Brunt-Vaisala frequency squared (1/s**2), a local array
!
!  WORK ARRAYS:
!
!    tmprtr   Work array for temperature
!    qvs      work array for saturation specific humidity
!    qw       Work array for total water mixing ratio
!    pt       Work arrays for total potential temperature
!
!   (These arrays are defined and used locally (i.e. inside this
!    subroutine), they may also be passed into routines called by
!    this one. Exiting the call to this subroutine, these temporary
!    work arrays may be used for other purposes and therefore their
!    contents may be overwritten. Please examine the usage of work
!    arrays before you alter the code.)
!
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
  INCLUDE 'grid.inc'          ! Grid & map parameters.
  INCLUDE 'phycst.inc'

!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  INTEGER :: nx,ny,nz          ! Number of grid points in 3 directions

  REAL :: ptprt (nx,ny,nz)     ! Perturbation potential temperature (K)
  REAL :: pprt  (nx,ny,nz)     ! Perturbation pressure (Pascal)
  REAL :: ptbar (nx,ny,nz)     ! Base state potential temperature (K)
  REAL :: qv    (nx,ny,nz)     ! Water vapor specific humidity (kg/kg)

  REAL :: qscalar(nx,ny,nz,nscalar)

  REAL :: qvbar (nx,ny,nz)     ! Base state water vapor specific humidity
  REAL :: pbar  (nx,ny,nz)     ! Base state pressure (Pascal)
  REAL :: j3    (nx,ny,nz)     ! Coordinate transformation Jacobian defined as
                               ! d( zp )/d( z ).
  REAL :: j3inv (nx,ny,nz)     ! Inverse of j3

  REAL :: nsqed (nx,ny,nz)     ! Brunt-Vaisala frequency squared (1/s**2)

  REAL :: tmprtr(nx,ny,nz)     ! Temporary work array
  REAL :: qvs   (nx,ny,nz)     ! Temporary work array
  REAL :: qw    (nx,ny,nz)     ! Temporary work array
  REAL :: pt    (nx,ny,nz)     ! Temporary work array
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i, j, k

  REAL :: ppi, dzinvd2, p0inv

  REAL :: eps   ! A small value of qc above which cloudwater is
                ! regarded as present.
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  p0inv = 1.0/p0
!
!-----------------------------------------------------------------------
!
!  Calculate the static stability N**2.
!  The moist static stability follows that of Durran and Klemp (1982)
!
!-----------------------------------------------------------------------
!
  DO k=1,nz-1
    DO j=1,ny-1
      DO i=1,nx-1

        pt(i,j,k) =  ptbar(i,j,k)+ptprt(i,j,k)

      END DO
    END DO
  END DO

  dzinvd2=1.0/(2*dz)

  IF( moist == 0 ) THEN

    DO k=2,nz-2
      DO j=1,ny-1
        DO i=1,nx-1
          nsqed(i,j,k)=g*(pt(i,j,k+1)-pt(i,j,k-1))*dzinvd2              &
                       /(pt(i,j,k)*j3(i,j,k))
        END DO
      END DO
    END DO

  ELSE

    DO k=1,nz-1
      DO j=1,ny-1
        DO i=1,nx-1

          ppi   = ((pbar(i,j,k)+pprt(i,j,k))*p0inv) ** rddcp
          tmprtr(i,j,k) = (ptbar(i,j,k)+ptprt(i,j,k)) * ppi

        END DO
      END DO
    END DO


    CALL getqvs(nx,ny,nz, 1,nx-1,1,ny-1,1,nz-1, pbar, tmprtr,qvs)

    DO k=1,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          qw(i,j,k) = qvs(i,j,k)
        END DO
      END DO
    END DO

    IF (P_QC > 0) THEN

      DO k=1,nz-1
        DO j=1,ny-1
          DO i=1,nx-1
            qw(i,j,k) = qw(i,j,k)+qscalar(i,j,k,P_QC)
          END DO
        END DO
      END DO

    END IF

    IF (P_QR > 0) THEN

      DO k=1,nz-1
        DO j=1,ny-1
          DO i=1,nx-1
            qw(i,j,k) = qw(i,j,k)+qscalar(i,j,k,P_QR)
          END DO
        END DO
      END DO

    END IF

    eps = 1.0E-6
    DO k=2,nz-2
      DO j=1,ny-1
        DO i=1,nx-1

          IF( (P_QC > 0 .OR. P_QR > 0) .AND. (qscalar(i,j,k,P_QC) > eps)) THEN

            nsqed(i,j,k)= g * j3inv(i,j,k) *                            &
                ((1+lathv*qvs(i,j,k)/(rd*tmprtr(i,j,k)))/               &
                (1+rd/rv*lathv*lathv*qvs(i,j,k)/(cp*rd*tmprtr(i,j,k)**2)) &
                *((pt(i,j,k+1)-pt(i,j,k-1))/pt(i,j,k)                   &
                +lathv/(cp*tmprtr(i,j,k))*(qvs(i,j,k+1)-qvs(i,j,k-1)) ) &
                -(qw(i,j,k+1)-qw(i,j,k-1)) ) *dzinvd2


          ELSE

            nsqed(i,j,k)=g*(pt(i,j,k+1)-pt(i,j,k-1))*dzinvd2            &
                         /(pt(i,j,k)*j3(i,j,k))

          END IF

        END DO
      END DO
    END DO

  END IF

  DO j=1,ny-1
    DO i=1,nx-1
      nsqed(i,j,  1 )=nsqed(i,j,2)
      nsqed(i,j,nz-1)=nsqed(i,j,nz-2)
    END DO
  END DO

  RETURN
END SUBROUTINE stabnsq
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE CFTMIX                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE cftmix(nx,ny,nz,                                             &
           nsqed,zp,tke,pbldpth,defsq, kmh,kmv,rprntl,lenscl,           &
           grdscl,tem1)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate the turbulent mixing coefficient, km, with the modified
!  Smagorinsky and TKE formulations.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  10/10/91
!
!  MODIFICATION HISTORY:
!
!  5/05/92 (M. Xue)
!  Added full documentation.
!
!  6/1/92 (M. Xue and H. Jin)
!  Further facelift.
!
!  6/1/93 (M. Xue)
!  Calculation of delta from dzp rather than dz.
!
!  6/16/93 (M. Xue)
!  The values of KM at the lateral boundaries are explicitly set.
!  The calculated values are overwritten.
!
!  4/20/1995 (M. Xue)
!  Added an upper limit on KM in ensure numerical stability.
!
!  3/8/96 (M. Xue, X. Song and V. Wong)
!  Add parameter tkeopt for three versions of 1.5 order TKE schemes.
!  They differ mainly in the specification of turbulent mixing length.
!
!  tkeopt = 1, Wyngaard formulation;
!         = 2, Deardroff/Moeng formulation;
!         = 3, Sun and Chang formulation.
!
!  3/11/96 (M. Xue)
!  Rewrote most of this subroutine.
!  Introduced kmh, kmv and rprntl arrays
!
!  3/21/96 (M. Xue)
!  Moved the calculation of magnitude of deformation (defsq) to
!  DEFORM.
!
!  8/5/96 (M. Xue and J. Zong)
!  Set an upper limit of 3.0 on the inverse Prandtl number for the
!  Sun and Chang formulation (tkeopt=3).
!
!  2/2/1998 (M. Xue and D. Weber)
!  Fixed a problem with mixing length calculation at the PBL top
!  when Sun and Chang scheme is used.
!
!  4/19/1999 (Pengfei Zhang)
!  Changed the calculation of Kmh and Kmv as frstep=1 and
!  tmixopt=4 by using Kmh=0.1*tke**2*lh and Kmv=0.1*tke**2*lv
!  in order to be consistent with future values.
!
!  6/1/2008 (Ming Xue)
!  Check if tkein is 1, i.e., tke is read in from the IC file.
!  If so, use it instead of doing special treatment for kmh,kmv etc.
!  for the first time step. Now using special_tmix instead of frstep
!  as the flag for the special treatment.
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!    nsqed    Brunt-Vaisala frequency (1/s**2), a local array
!    zp       Vertical coordinate of grid points in physical space (m)
!
!    tke      Turbulent Kinetic Energy ((m/s)**2)
!    pbldpth  Planetary boundary layer depth (m)
!    defsq    Deformation squared (1/s**2)
!
!  OUTPUT:
!
!    kmh      Horizontal turb. mixing coef. for momentum ( m**2/s )
!    kmv      Vertical turb. mixing coef. for momentum ( m**2/s )
!    rprntl   Reciprocal of Prandtl number
!    lenscl   Turbulent mixing length scale (m)
!
!  WORK ARRAYS:
!
!    grdscl   Work array to store grid scale.
!    tem1     Temporary array
!
!   (These arrays are defined and used locally (i.e. inside this
!    subroutine), they may also be passed into routines called by
!    this one. Exiting the call to this subroutine, these temporary
!    work arrays may be used for other purposes and therefore their
!    contents may be overwritten. Please examine the usage of work
!    arrays before you alter the code.)
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
!
  INTEGER :: nx,ny,nz          ! Number of grid points in 3 directions

  INCLUDE 'timelvls.inc'
!
  REAL :: nsqed (nx,ny,nz)     ! Brunt-Vaisala frequency (1/s**2)
  REAL :: zp    (nx,ny,nz)     ! The physical height coordinate defined at
                               ! w-point of the staggered grid.
!
  REAL :: tke   (nx,ny,nz,nt)  ! Turbulent Kinetic Energy ((m/s)**2)
  REAL :: pbldpth(nx,ny)       ! Planetary boundary layer depth (m)
  REAL :: defsq (nx,ny,nz)     ! Deformation squared (1/s**2)
!
  REAL :: kmh   (nx,ny,nz)     ! Horizontal turb. mixing coef. for
                               ! momentum. ( m**2/s )
  REAL :: kmv   (nx,ny,nz)     ! Vertical turb. mixing coef. for
                               ! momentum. ( m**2/s )
  REAL :: rprntl(nx,ny,nz)     ! Reciprocal of Prandtl number
  REAL :: lenscl(nx,ny,nz)     ! Turbulent mixing length scale (m)

  REAL :: grdscl(nx,ny,nz)     ! Array to store grid scale

  REAL :: tem1(nx,ny,nz)       ! Temporary array
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i, j, k
!
  REAL :: tema, temb
  REAL :: prinv, one3rd, km_nd_min,km_nd_max,km_nd_max2
  REAL :: deltah,deltah2, lpbltf, zs
  INTEGER :: special_tmix
  INTEGER :: tstrtlvl
  REAL    :: rprnum
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'grid.inc'          ! Grid & map parameters.
  INCLUDE 'phycst.inc'
  INCLUDE 'bndry.inc'
  INCLUDE 'mp.inc'            ! Message passing parameters.
  INCLUDE 'indtflg.inc'       ! Need tkein flag
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
!  Set the mixing coefficient, km, to a constant value TMIXCST.
!
!-----------------------------------------------------------------------
!
  one3rd= 1.0/3.0
  prinv = 1.0/prantl
  km_nd_max = kmlimit * 0.125/dtbig

  IF(tmixopt == 0) THEN

    DO k=1,nz
      DO j=1,ny
        DO i=1,nx
          kmh(i,j,k)= 0.0
          kmv(i,j,k)= 0.0
          rprntl(i,j,k)=prinv
          lenscl(i,j,k)= dx     ! not used for this option
        END DO
      END DO
    END DO

    RETURN

  ELSE IF(tmixopt == 1) THEN

    DO k=1,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          kmh(i,j,k)=tmixcst
          kmv(i,j,k)=tmixcst    ! Always assume isotropic
          rprntl(i,j,k)=prinv
          lenscl(i,j,k)= dx     ! not used for this option
        END DO
      END DO
    END DO

    GO TO 3000

  END IF

  deltah = SQRT( dx*dy )
  deltah2= dx*dy

  CALL stgrdscl(nx,ny,nz, zp, grdscl)
!
!-----------------------------------------------------------------------
!
!  Calculate the turbulent mixing coefficient, km, using the
!  modified Smagorinsky formulation.
!
!-----------------------------------------------------------------------
!

  IF((tmixopt==4 .OR. tmixopt==5).AND.( (ABS(curtim-tstart) <= 1.0E-10) .AND. (restrt /= 1) &
      .AND.(.NOT.(initopt==3.AND.tkein==1)) ) ) THEN
    special_tmix =1    ! Indicate that this is the initial step of
                       ! model integration.
  ELSE                 ! For non-first step or restart run
    special_tmix =0
  END IF

  IF((tmixopt == 3).OR.(tmixopt == 2).OR. special_tmix==1 ) THEN

    IF(tmixopt /= 3) tmixcst=0.0

    IF (trbisotp == 1) THEN

      DO k=1,nz-1
        DO j=1,ny-1
          DO i=1,nx-1
!
!  An upper bound is imposed on the KM value for numerical stability:
!
            tema = grdscl(i,j,k)*grdscl(i,j,k)

            kmh(i,j,k)= MIN(km_nd_max*tema, 0.0441*tema*                &
                SQRT(MAX(defsq(i,j,k)-nsqed(i,j,k)*prinv,0.0))          &
                +tmixcst)

            IF (tmixopt < 5) kmv(i,j,k)= kmh(i,j,k)

!
!  Set the initial value of tke (at tpast and tpresent) to be
!  consistent with Smagorinsky value through Kolmogorov-Prandtl relation
!
!  This could be a potential problem if users try to use tke as a
!  unused array when turn off the tke option.
!
            IF(tmixopt == 4) THEN
              tke(i,j,k,tpresent)=(kmh(i,j,k)*10./grdscl(i,j,k))**2
              tke(i,j,k,tpast   )=(kmh(i,j,k)*10./grdscl(i,j,k))**2
            END IF

            IF (tmixopt < 5) rprntl(i,j,k)=prinv
            lenscl(i,j,k)=grdscl(i,j,k) ! not used for this option

          END DO
        END DO
      END DO

    ELSE

      temb = 1.0/deltah2

      DO k=1,nz-1
        DO j=1,ny-1
          DO i=1,nx-1
!
!  An upper bound is imposed on the KM value for numerical stability:
!

            kmh(i,j,k)= MIN(km_nd_max*deltah2, 0.0441*deltah2*          &
                SQRT(MAX(defsq(i,j,k)-nsqed(i,j,k)*prinv,0.0))+tmixcst)

            tema = grdscl(i,j,k)*grdscl(i,j,k)

            IF (tmixopt < 5) kmv(i,j,k)= MIN(km_nd_max*tema, 0.0441*tema* &
                SQRT(MAX(defsq(i,j,k)-nsqed(i,j,k)*prinv,0.0))          &
                +tmixcst*tema*temb )

            IF(tmixopt == 4) THEN
              tke(i,j,k,tpresent)=(kmv(i,j,k)*10./grdscl(i,j,k))**2
              tke(i,j,k,tpast   )=(kmv(i,j,k)*10./grdscl(i,j,k))**2

              kmh(i,j,k)= 0.1*SQRT(tke(i,j,k,tpresent)*deltah2)
            END IF

            IF (tmixopt < 5) rprntl(i,j,k)=prinv
            lenscl(i,j,k)=grdscl(i,j,k) ! not used for this option

          END DO
        END DO
      END DO

    END IF

  END IF
!
!-----------------------------------------------------------------------
!
!  1.5 order TKE option:
!
!-----------------------------------------------------------------------
!
  IF(sadvopt == 4 .or. tintegopt == 2 .or. tintegopt == 3) THEN ! Forward-based FCT scheme
    tstrtlvl = tpresent
  ELSE
    tstrtlvl = tpast
  END IF
!WRITE(*,*) '==5==',MAXVAL(kmv),trbisotp

  IF( (tmixopt == 4 .OR. tmixopt==5).AND. special_tmix == 0 ) THEN

    DO k=1,nz-1
      DO j=1,ny-1
        DO i=1,nx-1

          IF( nsqed(i,j,k) > 0.0 ) THEN
            lenscl(i,j,k)=MIN(0.76*SQRT(tke(i,j,k,tstrtlvl)             &
                          /nsqed(i,j,k)),grdscl(i,j,k))
          ELSE
            lenscl(i,j,k)=grdscl(i,j,k)
          END IF

        END DO
      END DO
    END DO

    IF (tkeopt == 3) THEN
!
!-----------------------------------------------------------------------
!
!  For tkeopt=3, set length scale inside and immediately above the
!  PBL using a profile after Sun and Chang (1986).
!
!-----------------------------------------------------------------------
!
      lpbltf = lsclpbl0 * 1.8*(1.0-EXP(-4.0)-0.0003*EXP(8.0))

      DO j=1,ny-1
        DO i=1,nx-1
!
!-----------------------------------------------------------------------
!
!  When the PBL is convectively stable, the diagnosed PBL depth should
!  have been set to the thickness of the layer below first scalar
!  point. In this case, the profile will not be used.
!
!-----------------------------------------------------------------------
!
          IF( pbldpth(i,j)+zp(i,j,2) > 0.51*(zp(i,j,3)+zp(i,j,2)) ) THEN

            DO k=1,nz-1


              zs = (zp(i,j,k)+zp(i,j,k+1))*0.5
              tema=(zs-zp(i,j,2))/pbldpth(i,j)

              IF (tema <= 1.0)THEN  ! Below PBL top, use exp expression

                lenscl(i,j,k)=lsclpbl0*(1.8*pbldpth(i,j)*               &
                            (1.0-EXP(-4.0*tema)-0.0003*EXP(8.0*tema)))

              ELSE  ! Using linear function
!
!-----------------------------------------------------------------------
!
!  lenscl linearly decrease to zero from 1.0 pbldpth to 1.2 pbldpth.
!
!-----------------------------------------------------------------------
!
                temb=(lpbltf * pbldpth(i,j))                            &
                    *(1.2*pbldpth(i,j)-(zs-zp(i,j,2)))/(0.2*pbldpth(i,j))

                lenscl(i,j,k)= MAX( lenscl(i,j,k), temb )

              END IF
!
!-----------------------------------------------------------------------
!
!  Force lenscl to be smaller than the larger of grid scale and 600m.
!
!-----------------------------------------------------------------------
!
              lenscl(i,j,k)= MIN( lenscl(i,j,k),MAX(grdscl(i,j,k),600.0) )

            END DO
          END IF
        END DO
      END DO

    END IF
!
    km_nd_min=(1.e-6)
!WRITE(*,*) '==7==',MAXVAL(kmv),trbisotp

    IF (trbisotp == 1) THEN

      DO k=1,nz-1
        DO j=1,ny-1
          DO i=1,nx-1

            tema = grdscl(i,j,k)*grdscl(i,j,k)

            rprnum = 1.+2.*lenscl(i,j,k)/grdscl(i,j,k)

            rprnum = MIN( rprnum, 3.0 )

            kmh(i,j,k)=0.1*SQRT(tke(i,j,k,tstrtlvl))*lenscl(i,j,k)

            IF( rprnum*nsqed(i,j,k) < defsq(i,j,k))                     &
                kmh(i,j,k)=MAX(km_nd_min*tema,kmh(i,j,k))

            kmh(i,j,k)=MIN(km_nd_max*tema,kmh(i,j,k))

            IF (tmixopt == 4) THEN
              rprntl(i,j,k) = rprnum
              kmv(i,j,k)    = kmh(i,j,k)
            END IF

          END DO
        END DO
      END DO

    ELSE

      deltah  = SQRT(dx*dy)
      deltah2 = dx*dy

      km_nd_max2 = km_nd_max
      IF (trbvimp == 1) km_nd_max2 = 1000.0

      DO k=1,nz-1
        DO j=1,ny-1
          DO i=1,nx-1

            tema = grdscl(i,j,k)*grdscl(i,j,k)

            rprnum = 1.+2.*lenscl(i,j,k)/grdscl(i,j,k)

            rprnum = MIN( rprnum, 3.0 )

            kmh(i,j,k)=0.1*SQRT(tke(i,j,k,tstrtlvl))*deltah
            IF (tmixopt == 4) kmv(i,j,k)=0.1*SQRT(tke(i,j,k,tstrtlvl))*lenscl(i,j,k)

            IF( rprnum*nsqed(i,j,k) < defsq(i,j,k)) THEN
              kmh(i,j,k)=MAX(km_nd_min*deltah2,kmh(i,j,k))
              IF (tmixopt == 4) kmv(i,j,k)=MAX(km_nd_min*tema   ,kmv(i,j,k))
            END IF

            kmh(i,j,k)=MIN(km_nd_max*deltah2, kmh(i,j,k))

            IF (tmixopt == 4) THEN
              rprntl(i,j,k) = rprnum
              tema=km_nd_max2*tema
              kmv(i,j,k)=MIN(tema   , kmv(i,j,k))
            END IF

          END DO
        END DO
      END DO

    END IF

  END IF
!WRITE(*,*) '==6==',MAXVAL(kmv)

  3000  CONTINUE
!
!-----------------------------------------------------------------------
!
!  Reset the lateral boundary values of km and rprntl
!
!-----------------------------------------------------------------------
!

  CALL acct_interrupt(bc_acct)
  CALL bckmkh(nx,ny,nz,kmh)
  CALL bckmkh(nx,ny,nz,kmv)
  CALL bckmkh(nx,ny,nz,rprntl)
  CALL acct_stop_inter

!WRITE(*,*) '==9==',MAXVAL(kmv)
  IF (mp_opt > 0) THEN
    CALL acct_interrupt(mp_acct)
    CALL mpsendrecv2dew(kmh,nx,ny,nz,ebc,wbc,0,tem1)
    CALL mpsendrecv2dns(kmh,nx,ny,nz,nbc,sbc,0,tem1)

    CALL mpsendrecv2dew(kmv,nx,ny,nz,ebc,wbc,0,tem1)
    CALL mpsendrecv2dns(kmv,nx,ny,nz,nbc,sbc,0,tem1)

    CALL mpsendrecv2dew(rprntl,nx,ny,nz,ebc,wbc,0,tem1)
    CALL mpsendrecv2dns(rprntl,nx,ny,nz,nbc,sbc,0,tem1)
  END IF
!WRITE(*,*) '==8==',MAXVAL(kmv)

  RETURN
END SUBROUTINE cftmix
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE STGRDSCL                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE stgrdscl(nx,ny,nz, zp, grdscl )
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Set the grid scale to be used in turbulence mixing calculations.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  3/11/96
!
!  MODIFICATION HISTORY:
!
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!    zp       Vertical coordinate of grid points in physical space (m)
!
!  OUTPUT:
!
!    grdscl   The grid scale using turbulence mixing calculations (m)
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
  INTEGER :: nx,ny,nz          ! Number of grid points in 3 directions

  REAL :: zp    (nx,ny,nz)     ! The physical height coordinate defined at
                               ! w-point of the staggered grid.

  REAL :: grdscl(nx,ny,nz)     ! Grid scale (m)
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i, j, k
!
  REAL :: deltah, one3rd, dzp, dxdy3rd
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'grid.inc'          ! Grid & map parameters.
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
!  Set the mixing coefficient, km, to a constant value TMIXCST.
!
!-----------------------------------------------------------------------
!
  one3rd = 1.0/3.0
  deltah = SQRT( dx*dy )
  dxdy3rd =(dx*dy)**one3rd

  IF (trbisotp == 1) THEN  ! isotropic case

    DO k=1,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          dzp = zp(i,j,k+1)-zp(i,j,k)
          IF( runmod == 1 ) THEN
            grdscl(i,j,k) = dxdy3rd * dzp**one3rd
          ELSE IF( runmod == 2 ) THEN
            grdscl(i,j,k) = SQRT(dx*dzp)
          ELSE IF( runmod == 3 ) THEN
            grdscl(i,j,k) = SQRT(dy*dzp)
          ELSE IF( runmod == 4 ) THEN
            grdscl(i,j,k) = dzp
          END IF
        END DO
      END DO
    END DO

  ELSE                       ! anisotropic case

    DO k=1,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          grdscl(i,j,k) = zp(i,j,k+1)-zp(i,j,k)
        END DO
      END DO
    END DO

  END IF

  RETURN
END SUBROUTINE stgrdscl
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE DEFORM                     ######
!######                                                      ######
!######                     Developed by                     ######
!######  Center for the Analysis and Prediction of Storms    ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
!

SUBROUTINE deform(nx,ny,nz,u,v,w,mapfct,                                &
           j1,j2,j3,aj3x,aj3y,aj3z,j3inv,                               &
           d11,d12,d13,d22,d23,d33,d31,d32,defsq,                       &
           tem1,tem2,tem3,tem4,mp_tem)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate the deformation tensor components Dij.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  10/10/91
!
!  MODIFICATION HISTORY:
!
!  5/05/92 (M. Xue)
!  Added full documentation.
!
!  6/1/92 (M. Xue and H. Jin)
!  Further facelift.
!
!  6/25/92 (M. Xue and D. Weber)
!  Terrain included.
!
!  1/23/96 (Donghai Wang & Yuhe Liu)
!  Added the map projection factor.
!
!  3/21/96 (M. Xue)
!  Moved the calculation of magnitude of deformation (defsq) to
!  from CFTMIX to this routine.
!
!  4/1/96 (Donghai Wang, X. Song and M. Xue)
!  Added a time average weighting coefficient for implicit treatment
!  of vertical mixing.
!
!  6/5/96 (Donghai Wang, M. Xue and X. Song)
!  Fixed a minor bug related to the vertical implicit treatment of
!  turbulence mixing. Should not affect the results much.
!
!  11/06/97 (D. Weber)
!  Added three additional levels to the mapfct array.  The three
!  levels (4,5,6) represent the inverse of the first three in order.
!  The inverse map factors are computed to improve efficiency.
!
!  9/18/98 (D. Weber)
!  Modified do loop structure to improve code efficiency via:
!  -merging existing loops together
!  -removing and merging operators
!  -switched to DO  ENDDO structure for f90 conversion
!  Added tmixvert option for vertical tmixing only (for large aspect
!  ratio simulations only dx,dy>>dz).
!  -added aj3x,y,z arrays.
!
!  2/15/2002 (Yunheng Wang)
!  Fixed a typo in the first call of mpsend2dew and mprecv2dew
!  respectively (sbc replaced with wbc).
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
!    w        Vertical component of velocity in Cartesian
!             coordinates at a given time level (m/s).
!
!    mapfct   Map factors at scalar, u and v points
!
!    j1       Coordinate transform Jacobian -d(zp)/dx
!    j2       Coordinate transform Jacobian -d(zp)/dy
!    j3       Coordinate transform Jacobian  d(zp)/dz
!    aj3x     Avgx of the coordinate transformation Jacobian  d(zp)/dz
!    aj3y     Avgy of the coordinate transformation Jacobian  d(zp)/dz
!    aj3z     Avgz of the coordinate transformation Jacobian  d(zp)/dz
!    j3inv    Inverse of j3
!
!  OUTPUT:
!
!    d11      Deformation tensor component (1/s)
!    d12      Deformation tensor component (1/s)
!    d13      Deformation tensor component (1/s)
!    d22      Deformation tensor component (1/s)
!    d23      Deformation tensor component (1/s)
!    d31      Deformation tensor component (1/s)
!    d32      Deformation tensor component (1/s)
!    d33      Deformation tensor component (1/s)
!    defsq    Deformation squared (1/s**2)
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
  INTEGER :: nx,ny,nz          ! Number of grid points in 3 directions

  REAL :: u     (nx,ny,nz)     ! Total u-velocity (m/s)
  REAL :: v     (nx,ny,nz)     ! Total v-velocity (m/s)
  REAL :: w     (nx,ny,nz)     ! Total w-velocity (m/s)

  REAL :: mapfct(nx,ny,8)      ! Map factors at scalar, u and v points

  REAL :: j1    (nx,ny,nz)     ! Coordinate transform Jacobian defined as
                               ! - d( zp )/d( x ).
  REAL :: j2    (nx,ny,nz)     ! Coordinate transform Jacobian defined as
                               ! - d( zp )/d( y ).
  REAL :: j3    (nx,ny,nz)     ! Coordinate transform Jacobian defined as
                               ! d( zp )/d( z ).
  REAL :: aj3x  (nx,ny,nz)     ! Coordinate transformation Jacobian defined
                               ! as d( zp )/d( z ) AVERAGED IN THE X-DIR.
  REAL :: aj3y  (nx,ny,nz)     ! Coordinate transformation Jacobian defined
                               ! as d( zp )/d( z ) AVERAGED IN THE Y-DIR.
  REAL :: aj3z  (nx,ny,nz)     ! Coordinate transformation Jacobian defined
                               ! as d( zp )/d( z ) AVERAGED IN THE Z-DIR.
  REAL :: j3inv (nx,ny,nz)     ! Inverse of j3
!
  REAL :: d11   (nx,ny,nz)     ! Deformation tensor component (1/s)
  REAL :: d12   (nx,ny,nz)     ! Deformation tensor component (1/s)
  REAL :: d13   (nx,ny,nz)     ! Deformation tensor component (1/s)
  REAL :: d22   (nx,ny,nz)     ! Deformation tensor component (1/s)
  REAL :: d23   (nx,ny,nz)     ! Deformation tensor component (1/s)
  REAL :: d31   (nx,ny,nz)     ! Deformation tensor component (1/s)
  REAL :: d32   (nx,ny,nz)     ! Deformation tensor component (1/s)
  REAL :: d33   (nx,ny,nz)     ! Deformation tensor component (1/s)
  REAL :: defsq (nx,ny,nz)     ! Deformation squared (1/s**2)
!
  REAL :: tem1(nx,ny,nz)       ! Temproary working array
  REAL :: tem2(nx,ny,nz)       ! Temproary working array
  REAL :: tem3(nx,ny,nz)       ! Temproary working array
  REAL :: tem4(nx,ny,nz)       ! Temproary working array
  REAL :: mp_tem(nx,ny,nz)     ! Message passing work array
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k
  REAL :: tema,temb,temc
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'grid.inc'          ! Grid & map parameters.
  INCLUDE 'phycst.inc'
  INCLUDE 'bndry.inc'
  INCLUDE 'mp.inc'            ! Message passing parameters.

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!
!-----------------------------------------------------------------------
!
!  Calculate the deformation tensor components Dij, which are stored
!  in arrays dij.
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!
!  Calculate the stress tensor component d33
!  d33 = 2./j3 * difz(w)
!
!-----------------------------------------------------------------------

  tema = 2.0*dzinv
  DO k=1,nz-1
    DO j=1,ny-1
      DO i=1,nx-1
        d33(i,j,k)=tema*j3inv(i,j,k) * (w(i,j,k+1)-w(i,j,k))
      END DO
    END DO
  END DO

!-----------------------------------------------------------------------
!
!  Calculate the deformation tensor component D11:
!
!   d11 = (2./j3) * m * m * (difx(avgsu(j3) * (u/m)) +
!                           (difz(avgx(j1 * avgsw(u/m)))))
!
!-----------------------------------------------------------------------
!

!-----------------------------------------------------------------------
!
!  Calculate difx(avgsu(j3) * (u/m))
!
!-----------------------------------------------------------------------

  IF(tmixvert == 0)THEN   !  compute the horizontal terms also....

    DO k=1,nz-1     ! difx (avgx(j3)*u*mapfct)
      DO j=1,ny-1
        DO i=1,nx-1
          d11(i,j,k) = dxinv*(aj3x(i+1,j,k)*u(i+1,j,k)*mapfct(i+1,j,5)  &
                            - aj3x(i,j,k)  *u(i,j,k)  *mapfct(i,j,5) )
        END DO
      END DO
    END DO

!
!  At this point, D11 = difx(avgsu(j3) * (u/m))
!

    IF( ternopt /= 0 ) THEN  ! add in the terrain term.......

!-----------------------------------------------------------------------
!
!  Calculate the second term of d11
!  difz(avgx(j1 * avgsw(u/m)))  ! note 1/m comes outside difz...
!  NOTE:avgsw(u) is stored in d32 and will be used in the
!       third term of d12 later...
!
!-----------------------------------------------------------------------

      DO k=2,nz-1    ! average u in z-dir.
        DO j=1,ny-1
          DO i=1,nx
            d32(i,j,k) = 0.5*(u(i,j,k)+u(i,j,k-1))
          END DO
        END DO
      END DO

      CALL acct_interrupt(bc_acct)
      CALL bcsw(nx,ny,nz,1,nx,1,ny-1,tbc,bbc,d32)
      CALL acct_stop_inter

      DO k=1,nz      ! average j1*d32 in x-dir.
        DO j=1,ny-1
          DO i=1,nx-1
            tem2(i,j,k) = 0.5*(j1(i+1,j,k)*d32(i+1,j,k)                 &
                             + j1(i  ,j,k)*d32(i  ,j,k))
          END DO
        END DO
      END DO

!-----------------------------------------------------------------------
!
!  Combine the terms, ( 2 /j3)* m * (m * d11 + difz(tem2)) to form
!  the final D11 stress tensor
!
!-----------------------------------------------------------------------

      tema = 2.0*dzinv
      DO k=1,nz-1      ! difz(tem2).
        DO j=1,ny-1
          DO i=1,nx-1
            d11(i,j,k) = 2.0*j3inv(i,j,k)*mapfct(i,j,1)*                &
                (mapfct(i,j,1)*d11(i,j,k)+dzinv*(tem2(i,j,k+1)-tem2(i,j,k)))
          END DO
        END DO
      END DO   ! d11 with terrain is complete...

    ELSE

      DO k=1,nz-1
        DO j=1,ny-1
          DO i=1,nx-1
            d11(i,j,k) = 2.0*j3inv(i,j,k)*mapfct(i,j,7)*d11(i,j,k)
          END DO
        END DO
      END DO   ! d11 w/o terrain is complete...


    END IF   ! end of terrain if block for d11.....

!-----------------------------------------------------------------------
!
!  Calculate the deformation tensor component d22:
!
!   d22 = (2./j3) * m * m *(dify(avgsv(j3) * (v/m)) +
!                           difz(avgy(j2) * avgsw(v/m)))
!
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!  Calculate the first term in D22
!  dify(avgsv(j3) * (v/m))
!  NOTE:avgsv(j3) is stored in d23 and will be used in the second
!       term of d12 later...
!
!-----------------------------------------------------------------------

    DO k=1,nz-1        ! dify(d23*v*mapfct)
      DO j=1,ny-1
        DO i=1,nx-1
          d22(i,j,k) = dyinv*(aj3y(i,j+1,k)*v(i,j+1,k)*mapfct(i,j+1,6)  &
                           -  aj3y(i,j,k)  *v(i,j,k)  *mapfct(i,j,6) )
        END DO
      END DO
    END DO

!
!  At this point, D22 = dify(avgsv(j3) * v/m)
!

    IF( ternopt /= 0 ) THEN

!-----------------------------------------------------------------------
!
!  Calculate the second term of d22
!  difz(avgy(j2 * avgsw(v/m)))  Note: 1/m comes outside difz..
!  NOTE:avgsw(v) is stored in d31 and will be used in the
!       third term of d12 later...
!
!-----------------------------------------------------------------------

      DO k=2,nz-1    ! average v in z-dir.
        DO j=1,ny
          DO i=1,nx-1
            d31(i,j,k) = 0.5*(v(i,j,k)+v(i,j,k-1))
          END DO
        END DO
      END DO

      CALL acct_interrupt(bc_acct)
      CALL bcsw(nx,ny,nz,1,nx-1,1,ny,tbc,bbc,d31)
      CALL acct_stop_inter

      DO k=1,nz      ! average j2*d31 in y-dir.
        DO j=1,ny-1
          DO i=1,nx-1
            tem2(i,j,k) = 0.5*(j2(i,j+1,k)*d31(i,j+1,k)                 &
                             + j2(i,j,k)  *d31(i,j,k))
          END DO
        END DO
      END DO

!-----------------------------------------------------------------------
!
!  Combine terms 2./j3 * m * (m*D22 + difz(tem2)) to obtain the final
!  D22 deformation tensor
!
!-----------------------------------------------------------------------

      DO k=1,nz-1
        DO j=1,ny-1
          DO i=1,nx-1
            d22(i,j,k)= 2.0*j3inv(i,j,k)*mapfct(i,j,1)*                 &
                          (mapfct(i,j,1)*d22(i,j,k)                     &
                   +dzinv*(tem2(i,j,k+1)-tem2(i,j,k)))
          END DO
        END DO
      END DO    ! d22 with terrain is complete

    ELSE

      DO k=1,nz-1
        DO j=1,ny-1
          DO i=1,nx-1
            d22(i,j,k)= 2.0*j3inv(i,j,k)*mapfct(i,j,7)*d22(i,j,k)
          END DO
        END DO
      END DO    ! d22 w/o terrain is complete

    END IF   ! end of d22 terrain if block.......


!-----------------------------------------------------------------------
!
!  Calculate the deformation tensor d12
!
!  d12 = m*m/(avgsv(avgsu(j3))) *
!         (dify(avgsu(j3) * u/m) + difx(avgsv(j3) * v/m)
!        + difz(avgsu(j2) * avgsv(avgsw(u/m))
!        +      avgsv(j1) * avgsu(avgsw(v/m))))
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!  Note: In order to reduce memory usage, order dependence has been
!  introduced into the calculation of D12, which reduces the possible
!  term-wise parallellism.  This is due to the constraint of only two
!  temporary arrays being avaliable for passage into this subroutine.
!
!-----------------------------------------------------------------------

    IF( ternopt /= 0 ) THEN

!-----------------------------------------------------------------------
!
!  Calculate the first sub-term of the third term of D12
!  avgsu(j2) * avgsv(avgsw(u/m)), this term is zero when ternopt=0
!  and
!  Calculate the second sub-term in the third term of d12
!  avgsv(j1) * avgsu(avgsw(v/m)), this term is zero when ternopt=0
!  Note: from above avgsw(u) = d32
!  Note: from above avgsw(v) = d31
!  NOTE: ORDER IS VERY IMPORTANT!!
!
!-----------------------------------------------------------------------

      DO k=1,nz      ! average j2 in x-dir.
        DO j=1,ny
          DO i=2,nx-1
            tem1(i,j,k) = 0.5*(j2(i,j,k)+j2(i-1,j,k))
          END DO
        END DO
      END DO

      IF (mp_opt > 0) THEN
        CALL acct_interrupt(mp_acct)
        CALL mpsendrecv2dew(tem1,nx,ny,nz,ebc,wbc,1,mp_tem)
      END IF
      CALL acct_interrupt(bc_acct)
      CALL bcsu(nx,ny,nz,1,ny,1,nz,ebc,wbc,tem1)
      CALL acct_stop_inter

      DO k=1,nz    ! average u in y-dir.
        DO j=2,ny-1
          DO i=1,nx
            tem2(i,j,k) = 0.5*(d32(i,j,k)*mapfct(i,j,5)                 &
                              +d32(i,j-1,k)*mapfct(i,j-1,5))
          END DO
        END DO
      END DO   ! d32 is available for use...

      IF (mp_opt > 0) THEN
        CALL acct_interrupt(mp_acct)
        CALL mpsendrecv2dns(tem2,nx,ny,nz,nbc,sbc,2,mp_tem)
      END IF
      CALL acct_interrupt(bc_acct)
      CALL bcsv(nx,ny,nz,1,nx,1,nz,nbc,sbc,tem2)
      CALL acct_stop_inter

      DO k=1,nz    ! average v in x-dir.
        DO j=1,ny
          DO i=2,nx-1
            tem4(i,j,k) = 0.5*(d31(i,j,k)* mapfct(i,j,6)                &
                              +d31(i-1,j,k)* mapfct(i-1,j,6))
          END DO
        END DO
      END DO   ! d31 is available for use...

      IF (mp_opt > 0) THEN
        CALL acct_interrupt(mp_acct)
        CALL mpsendrecv2dew(tem4,nx,ny,nz,ebc,wbc,1,mp_tem)
      END IF
      CALL acct_interrupt(bc_acct)
      CALL bcsu(nx,ny,nz,1,ny,1,nz,ebc,wbc,tem4)
      CALL acct_stop_inter

      DO k=1,nz    ! average j1 in y-dir.
        DO j=2,ny-1
          DO i=1,nx
            tem3(i,j,k) = 0.5*(j1(i,j,k)+j1(i,j-1,k))
          END DO
        END DO
      END DO

      IF (mp_opt > 0) THEN
        CALL acct_interrupt(mp_acct)
        CALL mpsendrecv2dns(tem3,nx,ny,nz,nbc,sbc,2,mp_tem)
      END IF
      CALL acct_interrupt(bc_acct)
      CALL bcsv(nx,ny,nz,1,nx,1,nz,nbc,sbc,tem3)
      CALL acct_stop_inter

      DO k=1,nz-1    ! difz of the third term...
        DO j=1,ny
          DO i=1,nx
            d12(i,j,k) = dzinv*((tem1(i,j,k+1)*tem2(i,j,k+1)            &
                                +tem3(i,j,k+1)*tem4(i,j,k+1))           &
                               -(tem1(i,j,k)*tem2(i,j,k)                &
                                +tem3(i,j,k)*tem4(i,j,k)))
          END DO
        END DO
      END DO

    ELSE

      DO k=1,nz-1
        DO j=1,ny
          DO i=1,nx
            d12(i,j,k)=0.0
          END DO
        END DO
      END DO

    END IF    ! end of the terrain if block for d12....
!
!  At this point, d12 = the third term of d12
!

!-----------------------------------------------------------------------
!
!  Calculate the first and second terms in d12
!    (dify(avgsu(j3) * u/m) + difx(avgsv(j3) * v/m)
!  NOTE: D13 contains avgsu(j3) and
!        D23 contains avgsv(j3)
!
!-----------------------------------------------------------------------

    DO k=1,nz-1    ! dify(avgsu(j3) * u/m)
      DO j=2,ny-1
        DO i=1,nx
          tem2(i,j,k) = dyinv*(aj3x(i,j,k)*u(i,j,k)*mapfct(i,j,5)       &
                         - aj3x(i,j-1,k)*u(i,j-1,k)*mapfct(i,j-1,5))
        END DO
      END DO
    END DO    ! d13 is available for use

    IF (mp_opt > 0) THEN
      CALL acct_interrupt(mp_acct)
      CALL mpsendrecv2dns(tem2,nx,ny,nz,nbc,sbc,2,mp_tem)
    END IF
    CALL acct_interrupt(bc_acct)
    CALL boundv (tem2,nx,ny,nz, 1,nx, 1,nz-1)
    CALL acct_stop_inter

    DO k=1,nz-1    ! difx(avgsv(j3) * v/m)
      DO j=1,ny
        DO i=2,nx-1
          tem1(i,j,k) = dxinv*(aj3y(i,j,k)*v(i,j,k)*mapfct(i,j,6)       &
                         - aj3y(i-1,j,k)*v(i-1,j,k)*mapfct(i-1,j,6))
        END DO
      END DO
    END DO    ! d23 is available for use

    IF (mp_opt > 0) THEN
      CALL acct_interrupt(mp_acct)
      CALL mpsendrecv2dew(tem1,nx,ny,nz,ebc,wbc,1,mp_tem)
    END IF
    CALL acct_interrupt(bc_acct)
    CALL boundu (tem1,nx,ny,nz, 1,ny, 1,nz-1)
    CALL acct_stop_inter

!-----------------------------------------------------------------------
!
!  Combine the d12 terms
!
!-----------------------------------------------------------------------

    DO k=1,nz-1
      DO j=1,ny
        DO i=1,nx
          d12(i,j,k)= d12(i,j,k)+tem1(i,j,k)+tem2(i,j,k)
        END DO
      END DO
    END DO

!
!  At this point, d12 = first + second + third terms of d12
!

!-------------------------------------------------------------------
!
!  Calculate the coefficient of d12
!  avgsv(avgsu(j3*mapfct**2))
!
!-------------------------------------------------------------------

    DO k=1,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          tem2(i,j,k) = j3inv(i,j,k) * mapfct(i,j,7)
        END DO
      END DO
    END DO

    CALL avgsu(tem2,nx,ny,nz, 1,ny-1, 1,nz-1, tem1, tem4)
    CALL avgsv(tem1,nx,ny,nz, 1,nx,   1,nz-1, tem3, tem4)

!  Compute the final d12 = d12 * tem2

    DO k=1,nz-1
      DO j=1,ny
        DO i=1,nx
          d12(i,j,k)=d12(i,j,k)*tem3(i,j,k)
        END DO
      END DO
    END DO

  ELSE   ! zero the horizontal mixing terms....
         ! needed for the deformation computation below....

    DO k=1,nz-1
      DO j=1,ny
        DO i=1,nx
          d11(i,j,k)=0.0
          d22(i,j,k)=0.0
          d12(i,j,k)=0.0
        END DO
      END DO
    END DO

  END IF  !  end if tmixvert if block for d12....

!-----------------------------------------------------------------------
!
!  Calculate the deformation tensor d13
!  d13 = m/ avgsw(avgsu(j3)) * (difx(avgsw(j3) * w) +
!        difz(u/m + avgz(j1 * (avgsu(w)))))
!
!-----------------------------------------------------------------------

!-------------------------------------------------------------------
!
!  Compute the first term of d13
!  difx(avgsw(j3) * w)
!
!-------------------------------------------------------------------

  DO k=1,nz   ! difx (avgz(j3)*w)
    DO j=1,ny-1
      DO i=2,nx-1
        d13(i,j,k)=dxinv*(aj3z(i,j,k)*w(i,j,k)-                         &
                          aj3z(i-1,j,k)*w(i-1,j,k))
      END DO
    END DO
  END DO

  IF (mp_opt > 0) THEN
    CALL acct_interrupt(mp_acct)
    CALL mpsendrecv2dew(d13,nx,ny,nz,ebc,wbc,1,mp_tem)
  END IF
  CALL acct_interrupt(bc_acct)
  CALL boundu(d13,nx,ny,nz, 1,ny-1, 1,nz)
  CALL acct_stop_inter

!
!  At this point D23 = difx(avgsw(j3) * w)  for term d13
!

!-----------------------------------------------------------------------
!
!  Compute the second term of d13
!  difz(u + avgz(j1 * (avgsu(w))))
!
!-----------------------------------------------------------------------

  IF( ternopt /= 0 ) THEN

    DO k=1,nz    ! average w in x-dir.
      DO j=1,ny-1
        DO i=2,nx-1
          tem1(i,j,k) = 0.5*(w(i,j,k)+w(i-1,j,k))
        END DO
      END DO
    END DO

    IF (mp_opt > 0) THEN
      CALL acct_interrupt(mp_acct)
      CALL mpsendrecv2dew(tem1,nx,ny,nz,ebc,wbc,1,mp_tem)
    END IF
    CALL acct_interrupt(bc_acct)
    CALL bcsu(nx,ny,nz,1,ny-1,1,nz,ebc,wbc,tem1)
    CALL acct_stop_inter

    DO k=1,nz-1
      DO j=1,ny-1
        DO i=1,nx
          tem2(i,j,k) = 0.5*(tem1(i,j,k+1)*j1(i,j,k+1)                  &
                           + tem1(i,j,k)*j1(i,j,k) )
        END DO
      END DO
    END DO

!-----------------------------------------------------------------------
!
!  Compute u + tem2
!
!-----------------------------------------------------------------------

    DO k=1,nz-1
      DO j=1,ny-1
        DO i=1,nx
          tema = u(i,j,k)*mapfct(i,j,5)
          tem3(i,j,k)=tema+tem2(i,j,k)
          tem2(i,j,k)=alfcoef*tema+tem2(i,j,k)
        END DO
      END DO
    END DO

  ELSE   ! If ternopt=0, tem2=0

    DO k=1,nz-1
      DO j=1,ny-1
        DO i=1,nx
          tem3(i,j,k)=u(i,j,k)*mapfct(i,j,5)
          tem2(i,j,k)=alfcoef* tem3(i,j,k)
        END DO
      END DO
    END DO

  END IF

  DO k=2,nz-1
    DO j=1,ny-1
      DO i=1,nx
        tem1(i,j,k)=dzinv*(tem2(i,j,k)-tem2(i,j,k-1))
        tem4(i,j,k)=dzinv*(tem3(i,j,k)-tem3(i,j,k-1))
      END DO
    END DO
  END DO

  CALL acct_interrupt(bc_acct)
  CALL boundw(tem1,nx,ny,nz, 1,nx, 1,ny-1)
  CALL boundw(tem4,nx,ny,nz, 1,nx, 1,ny-1)
  CALL acct_stop_inter

!-----------------------------------------------------------------------
!
!  Compute d13 + tem1
!
!-----------------------------------------------------------------------

  DO k=1,nz
    DO j=1,ny-1
      DO i=1,nx
        d31(i,j,k)=d13(i,j,k)+tem4(i,j,k)
        d13(i,j,k)=d13(i,j,k)+tem1(i,j,k)
      END DO
    END DO
  END DO

!
!  At this point, D13 = first and second terms of d13
!

!-----------------------------------------------------------------------
!
!  Compute the coefficient of d13
!  avgsw(avgsu(j3))
!
!-----------------------------------------------------------------------

  DO k=2,nz-1
    DO j=1,ny-1
      DO i=1,nx
        tem2(i,j,k)=0.5*(aj3x(i,j,k)+aj3x(i,j,k-1))
      END DO
    END DO
  END DO

  CALL acct_interrupt(bc_acct)
  CALL bcsw(nx,ny,nz,1,nx,1,ny,tbc,bbc,tem2)
  CALL acct_stop_inter


!-----------------------------------------------------------------------
!
!  Compute the final d13 tensor = d13 / tem2
!  Note: the limits have been changed from 1,nz to 2,nz-1....
!
!-----------------------------------------------------------------------

  DO k=2,nz-1
    DO j=1,ny-1
      DO i=1,nx
        tema = mapfct(i,j,2)/tem2(i,j,k)
        d13(i,j,k)=d13(i,j,k)*tema
        d31(i,j,k)=d31(i,j,k)*tema
      END DO
    END DO
  END DO

!-----------------------------------------------------------------------
!
!  Calculate the deformation tensor d23
!  d23 = m/ avgsw(avgsv(j3)) * (dify(avgsw(j3) * w) +
!            difz(v/m + avgz(j2 * (avgsv(w)))))
!
!-----------------------------------------------------------------------

!-------------------------------------------------------------------
!
!  Compute the first term of d23
!  dify(avgsw(j3) * w)
!  NOTE: d32 = avgsw(j3)
!
!-------------------------------------------------------------------

  DO k=1,nz   ! dify (avgsw(j3)*w)
    DO j=2,ny-1
      DO i=1,nx-1
        d23(i,j,k)=dyinv*(aj3z(i,j,k)*w(i,j,k)-                         &
                          aj3z(i,j-1,k)*w(i,j-1,k))
      END DO
    END DO
  END DO

  IF (mp_opt > 0) THEN
    CALL acct_interrupt(mp_acct)
    CALL mpsendrecv2dns(d23,nx,ny,nz,nbc,sbc,2,mp_tem)
  END IF
  CALL acct_interrupt(bc_acct)
  CALL boundv(d23,nx,ny,nz, 1,nx-1, 1,nz)
  CALL acct_stop_inter

!
!  At this point d23 = dify(avgsw(j3) * w)
!

!-----------------------------------------------------------------------
!
!  Compute the second term of D23
!  difz(v + avgz(j2 * (avgsv(w))))
!
!-----------------------------------------------------------------------

  IF( ternopt /= 0 ) THEN

    DO k=1,nz
      DO j=2,ny-1
        DO i=1,nx-1
          tem1(i,j,k)=0.5*(w(i,j,k)+w(i,j-1,k))
        END DO
      END DO
    END DO

    IF (mp_opt > 0) THEN
      CALL acct_interrupt(mp_acct)
      CALL mpsendrecv2dns(tem1,nx,ny,nz,nbc,sbc,2,mp_tem)
    END IF
    CALL acct_interrupt(bc_acct)
    CALL bcsv(nx,ny,nz,1,nx-1,1,nz,nbc,sbc,tem1)
    CALL acct_stop_inter

    DO k=1,nz-1   ! avgz (j2*avgsv(w))
      DO j=1,ny
        DO i=1,nx-1
          tem2(i,j,k)=0.5*(j2(i,j,k)  *tem1(i,j,k)                      &
                         + j2(i,j,k+1)*tem1(i,j,k+1))
        END DO
      END DO
    END DO

    DO k=1,nz-1
      DO j=1,ny
        DO i=1,nx-1
          tema = v(i,j,k)*mapfct(i,j,6)
          tem3(i,j,k)=tema +tem2(i,j,k)
          tem2(i,j,k)=alfcoef * tema +tem2(i,j,k)
        END DO
      END DO
    END DO


  ELSE   ! When ternopt.eq.0, tem2=0

    DO k=1,nz-1
      DO j=1,ny
        DO i=1,nx-1
          tem3(i,j,k)=v(i,j,k)*mapfct(i,j,6)
          tem2(i,j,k)=alfcoef * tem3(i,j,k)
        END DO
      END DO
    END DO

  END IF

  DO k=2,nz-1
    DO j=1,ny
      DO i=1,nx-1
        tem1(i,j,k)=dzinv*(tem2(i,j,k)-tem2(i,j,k-1))
        tem4(i,j,k)=dzinv*(tem3(i,j,k)-tem3(i,j,k-1))
      END DO
    END DO
  END DO

  CALL acct_interrupt(bc_acct)
  CALL boundw(tem1,nx,ny,nz, 1,nx-1, 1,ny)
  CALL boundw(tem4,nx,ny,nz, 1,nx-1, 1,ny)
  CALL acct_stop_inter

!-----------------------------------------------------------------------
!
!    Compute d23 + tem1
!
!-----------------------------------------------------------------------

  DO k=1,nz
    DO j=1,ny
      DO i=1,nx-1
        d32(i,j,k)=d23(i,j,k)+tem4(i,j,k)
        d23(i,j,k)=d23(i,j,k)+tem1(i,j,k)
      END DO
    END DO
  END DO

!
!  At this point, d23 = first + second terms of d23
!

!-----------------------------------------------------------------------
!
!  Compute the coefficient of d23
!  avgsw(avgsv(j3))
!
!-----------------------------------------------------------------------

  DO k=2,nz-1
    DO j=1,ny
      DO i=1,nx-1
        tem2(i,j,k)= 0.5*(aj3y(i,j,k)+aj3y(i,j,k-1))
      END DO
    END DO
  END DO

  CALL acct_interrupt(bc_acct)
  CALL bcsw(nx,ny,nz,1,nx-1,1,ny,tbc,bbc,tem2)
  CALL acct_stop_inter

!-----------------------------------------------------------------------
!
!  Compute the final d23 deformation tensor d23 = d23 / tem2
!
!-----------------------------------------------------------------------

  DO k=1,nz
    DO j=1,ny
      DO i=1,nx-1
        tema = mapfct(i,j,3)/tem2(i,j,k)
        d23(i,j,k)=d23(i,j,k)*tema
        d32(i,j,k)=d32(i,j,k)*tema
      END DO
    END DO
  END DO

!-----------------------------------------------------------------------
!
!  Calculate the magnitude of deformation squared, Def**2.
!
!-----------------------------------------------------------------------

  DO k=1,nz-1
    DO j=1,ny-1
      DO i=1,nx-1
        tema = (d12(i,j,k)+d12(i,j+1,k))                                &
              +(d12(i+1,j,k)+d12(i+1,j+1,k))
        temb = (d31(i,j,k)+d31(i,j,k+1))                                &
              +(d31(i+1,j,k)+d31(i+1,j,k+1))
        temc = (d32(i,j,k)+d32(i,j+1,k))                                &
              +(d32(i,j,k+1)+d32(i,j+1,k+1))

        defsq(i,j,k) =(d11(i,j,k)*d11(i,j,k)                            &
                      +d22(i,j,k)*d22(i,j,k)                            &
                      +d33(i,j,k)*d33(i,j,k))*0.5                       &
                      +(tema*tema+temb*temb+temc*temc)*0.0625
      END DO
    END DO
  END DO

  DO k=1,nz-1
    DO j=1,ny-1
      DO i=1,nx-1
        d33(i,j,k)=alfcoef * d33(i,j,k)
      END DO
    END DO
  END DO

  RETURN
END SUBROUTINE deform


!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE STRESS                     ######
!######                                                      ######
!######                     Developed by                     ######
!######  Center for the Analysis and Prediction of Storms    ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE stress(nx,ny,nz,j3,j3inv,                                    &
           kmh,kmv,rhostr,tau11,tau12,tau13,tau22,tau23,tau33,          &
           tau31,tau32,                                                 &
           tem1, tem2, tem3, rjkmh, rjkmv)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate the stress tensor tauij from deformation tensor
!  Dij (input as tauij).
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  10/10/91
!
!  MODIFICATION HISTORY:
!
!  5/05/92 (M. Xue)
!  Added full documentation.
!
!  6/1/92 (M. Xue and H. Jin)
!  Further facelift.
!
!  6/26/92 (M. Xue and D. Weber)
!  Terrain included.
!
!  3/21/96 (Ming Xue)
!  Added work arrays rjkmh and rjkmv. Simplified the code.
!
!  6/5/96 (Donghai Wang, M. Xue and X. Song)
!  Fixed a minor bug related to the vertical implicit treatment of
!  turbulence mixing. Should not affect the results much.
!
!  8/27/98 (D. Weber)
!  Modified do loop structure to improve code efficiency via:
!  -merging existing loops together
!  -removing and merging operators
!  -switched to DO  ENDDO structure for f90 conversion
!  -added tem3
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!    j3       Coordinate transformation Jacobian  d(zp)/dz
!    j3inv    Inverse of j3
!
!    km       Turbulent mixing coefficient for momentum ( m**2/s ).
!    rhostr   Base state air density times j3 (kg/m**3).
!
!    tau11    Deformation tensor component (1/s)
!    tau12    Deformation tensor component (1/s)
!    tau13    Deformation tensor component (1/s)
!    tau22    Deformation tensor component (1/s)
!    tau23    Deformation tensor component (1/s)
!    tau33    Deformation tensor component (1/s)
!
!
!  OUTPUT:
!
!    tau11    Stress tensor component (kg/(m*s**2))
!    tau12    Stress tensor component (kg/(m*s**2))
!    tau13    Stress tensor component (kg/(m*s**2))
!    tau22    Stress tensor component (kg/(m*s**2))
!    tau23    Stress tensor component (kg/(m*s**2))
!    tau33    Stress tensor component (kg/(m*s**2))
!
!  WORKING ARRAYS:
!
!    tem1     Temporary working array
!    tem2     Temporary working array
!    tem3     Temporary working array
!    rjkmh    rhostr/j3*kmh
!    rjkmv    rhostr/j3*kmv
!
!  (These arrays are defined and used locally (i.e. inside this
!   subroutine), they may also be passed into routines called by
!   this one. Exiting the call to this subroutine, these temporary
!   working arrays may be used for other purposes therefore their
!   contents overwritten. Please examine the usage of working arrays
!   before you alter the code.)
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
  INTEGER :: nx,ny,nz          ! Number of grid points in 3 directions
!
  REAL :: kmh   (nx,ny,nz)     ! Horizontal turb. mixing coef. for
                               ! momentum. ( m**2/s )
  REAL :: kmv   (nx,ny,nz)     ! Vertical turb. mixing coef. for
                               ! momentum. ( m**2/s )

  REAL :: rhostr(nx,ny,nz)     ! Base state air density times j3(kg/m**3)
  REAL :: j3    (nx,ny,nz)     ! Coordinate transformation Jacobian defined as
                               ! d( zp )/d( z ).
  REAL :: j3inv (nx,ny,nz)     ! Inverse of j3
!
  REAL :: tau11 (nx,ny,nz)     ! Deformation and stress tensor
                               ! component (kg/(m*s**2))
  REAL :: tau12 (nx,ny,nz)     ! Deformation and stress tensor
                               ! component (kg/(m*s**2))
  REAL :: tau13 (nx,ny,nz)     ! Deformation and stress tensor
                               ! component (kg/(m*s**2))
  REAL :: tau22 (nx,ny,nz)     ! Deformation and stress tensor
                               ! component (kg/(m*s**2))
  REAL :: tau23 (nx,ny,nz)     ! Deformation and stress tensor
                               ! component (kg/(m*s**2))
  REAL :: tau33 (nx,ny,nz)     ! Deformation and stress tensor
                               ! component (kg/(m*s**2))
  REAL :: tau31 (nx,ny,nz)     ! Deformation and stress tensor
                               ! component (kg/(m*s**2))
  REAL :: tau32 (nx,ny,nz)     ! Deformation and stress tensor
                               ! component (kg/(m*s**2))

  REAL :: tem1  (nx,ny,nz)     ! Temporary working array
  REAL :: tem2  (nx,ny,nz)     ! Temporary working array
  REAL :: tem3  (nx,ny,nz)     ! Temporary working array
  REAL :: rjkmh (nx,ny,nz)     ! rhostr/j3*kmh
  REAL :: rjkmv (nx,ny,nz)     ! rhostr/j3*kmv
!
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
  INCLUDE 'phycst.inc'
  INCLUDE 'bndry.inc'
  INCLUDE 'mp.inc'            ! Message passing parameters.

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!
!-----------------------------------------------------------------------
!
!  Calculate the stress tensor Tau11, Tau22 and Tau33.
!
!  tau11 = rhostr/j3*kmh*d11
!  tau22 = rhostr/j3*kmh*d22
!  tau33 = rhostr/j3*kmv*d33
!
!-----------------------------------------------------------------------

  DO k=1,nz-1
    DO j=1,ny-1
      DO i=1,nx-1
        rjkmv(i,j,k)=rhostr(i,j,k)*j3inv(i,j,k)*kmv(i,j,k)
        tau33(i,j,k)=tau33(i,j,k)*rjkmv(i,j,k)
      END DO
    END DO
  END DO

  IF( tmixvert == 0)THEN   ! compute the horizontal terms...

    DO k=1,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          rjkmh(i,j,k)=rhostr(i,j,k)*j3inv(i,j,k)*kmh(i,j,k)
          tau11(i,j,k)=tau11(i,j,k)*rjkmh(i,j,k)
          tau22(i,j,k)=tau22(i,j,k)*rjkmh(i,j,k)
        END DO
      END DO
    END DO

!-----------------------------------------------------------------------
!
!  Calculate the stress tensor Tau12.
!  tau12 = avgsv(avgsu(rhostr/j3*kmh)) * d12
!
!-----------------------------------------------------------------------

    DO k=1,nz-1   ! note: do not overwrite tem3...
      DO j=1,ny-1
        DO i=2,nx-1
          tem3(i,j,k)=0.5*(rjkmh(i,j,k)+rjkmh(i-1,j,k))
        END DO
      END DO
    END DO        ! tem3 is used below...

    DO k=1,nz-1
      DO j=2,ny-1
        DO i=2,nx-1
          tem2(i,j,k)=0.5*(tem3(i,j,k)+tem3(i,j-1,k))
        END DO
      END DO
    END DO

    DO k=1,nz-1
      DO j=2,ny-1
        DO i=2,nx-1
          tau12(i,j,k)=tem2(i,j,k)*tau12(i,j,k)
        END DO
      END DO
    END DO

!-----------------------------------------------------------------------
!
!  Calculate the stress tensor Tau31.
!  tau31 = avgsw(avgsu(rhostr/j3*kmh)) * d31
!
!  Note:  avgsu(rhostr/j3*kmh) is stored in tem3 from above...
!
!-----------------------------------------------------------------------

    DO k=2,nz-1
      DO j=2,ny-2
        DO i=2,nx-1
          tem2(i,j,k)=0.5*(tem3(i,j,k)+tem3(i,j,k-1))
        END DO
      END DO
    END DO   ! tem3 is used below... -- ??? can't find it [GMB]

    CALL acct_interrupt(bc_acct)
    CALL bcsw(nx,ny,nz,2,nx-1,2,ny-2,tbc,bbc,tem2)
    CALL acct_stop_inter

    DO k=1,nz
      DO j=2,ny-2
        DO i=2,nx-1
          tau31(i,j,k)=tem2(i,j,k)*tau31(i,j,k)
        END DO
      END DO
    END DO

!-----------------------------------------------------------------------
!
!  Calculate the stress tensor Tau32.
!  tau32 = avgsw(avgsv(rhostr/j3*kmh)) * d32
!
!-----------------------------------------------------------------------

    DO k=2,nz-1
      DO j=2,ny-1
        DO i=2,nx-2
          tem1(i,j,k)=0.25*((rjkmh(i,j,k)  +rjkmh(i,j-1,k))             &
                           +(rjkmh(i,j,k-1)+rjkmh(i,j-1,k-1)))
        END DO
      END DO
    END DO

    CALL acct_interrupt(bc_acct)
    CALL bcsw(nx,ny,nz,2,nx-2,2,ny-1,tbc,bbc,tem1)
    CALL acct_stop_inter

    DO k=1,nz
      DO j=2,ny-1
        DO i=2,nx-2
          tau32(i,j,k)=tau32(i,j,k)*tem1(i,j,k)
        END DO
      END DO
    END DO

  END IF     ! end of tmixvert if block...

!-----------------------------------------------------------------------
!
!  Calculate the stress tensor Tau13.
!  tau13 = avgsw(avgsu(rhostr/j3*kmv)) * d13
!
!-----------------------------------------------------------------------

  DO k=2,nz-1
    DO j=2,ny-2
      DO i=2,nx-1
        tem2(i,j,k)=0.25*((rjkmv(i,j,k)+rjkmv(i-1,j,k))+                &
                          (rjkmv(i,j,k-1)+rjkmv(i-1,j,k-1)))
      END DO
    END DO
  END DO

  DO k=2,nz-1
    DO j=2,ny-2
      DO i=2,nx-1
        tau13(i,j,k)=tau13(i,j,k)*tem2(i,j,k)
      END DO
    END DO
  END DO

!-----------------------------------------------------------------------
!
!  Calculate the stress tensor Tau23.
!  tau23 = avgsw(avgsv(rhostr/j3*kmv)) * d23
!
!-----------------------------------------------------------------------

  DO k=2,nz-1
    DO j=2,ny-1
      DO i=2,nx-2
        tem1(i,j,k)=0.25*((rjkmv(i,j,k)+rjkmv(i,j-1,k))+                &
                          (rjkmv(i,j,k-1)+rjkmv(i,j-1,k-1)))
      END DO
    END DO
  END DO

  DO k=2,nz-1
    DO j=2,ny-1
      DO i=2,nx-2
        tau23(i,j,k)=tau23(i,j,k)*tem1(i,j,k)
      END DO
    END DO
  END DO

  RETURN
END SUBROUTINE stress


!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE TMIXPT                     ######
!######                                                      ######
!######                     Developed by                     ######
!######  Center for the Analysis and Prediction of Storms    ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE tmixpt (nx,ny,nz,ptprt,ptbar,rhostr,kmh,kmv,rprntl,          &
           usflx,vsflx,ptsflx,pbldpth,                                  &
           x,y,z,zp, mapfct, j1,j2,j3,aj3x,aj3y,j3inv,ptsfc,            &
           ptmix,                                                       &
           h1,h2,h3, tem1,tem2,tem3)
!
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate the turbulent mixing term for the potential
!  temperature equation. The term is expressed in terms of turbulent
!  heat fluxes.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  10/10/91
!
!  MODIFICATION HISTORY:
!
!  5/05/92 (M. Xue)
!
!  Added full documentation.
!
!  6/1/92 (M. Xue and H. Jin)
!  Further facelift.
!
!  7/20/92 (M. Xue and D. Weber)
!  Terrain included.
!
!  6/16/93 (MX)
!  Break down into subroutine calls.
!
!  1/23/96 (Donghai Wang & Yuhe Liu)
!  Added the map projection factor.
!
!  3/16/96 (Ming Xue)
!  The model can now include the surface fluxes with the turbulent
!  mxing option off (tmixopt=0).
!
!  3/27/1998 (M. Xue)
!  Added option tqflxdis for quadratic distribution of heat
!  and moisture fluxes in depth dtqflxdis, which is typically
!  set to 200 m.
!
!  9/18/98 (D. Weber)
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
!    ptprt    Perturbation potential temperature at a given time level (K)
!    ptbar    Base state potential temperature (K)
!    rhostr   Base state air density times j3 (kg/m**3)
!
!    usflx    Surface flux of u-momentum
!    vsflx    Surface flux of v-momentum
!
!    x        x coordinate of grid points in physical/comp. space (m)
!    y        y coordinate of grid points in physical/comp. space (m)
!    z        z coordinate of grid points in computational space (m)
!    zp       Vertical coordinate of grid points in physical space (m)
!
!    mapfct   Map factors at scalar points
!
!    j1       Coordinate transform Jacobian -d(zp)/dx
!    j2       Coordinate transform Jacobian -d(zp)/dy
!    j3       Coordinate transform Jacobian  d(zp)/dz
!    aj3x     Avgx of the coordinate transformation Jacobian  d(zp)/dz
!    aj3y     Avgy of the coordinate transformation Jacobian  d(zp)/dz
!    j3inv    Inverse of j3
!
!    kmh      Horizontal turb. mixing coef. for momentum ( m**2/s )
!    kmv      Vertical turb. mixing coef. for momentum ( m**2/s )
!    rprntl   Reciprocal of Prandtl number
!
!    ptsflx   Surface flux of heat (K*kg/(m**2*s))
!    pbldpth  Planetary boundary layer depth (m)
!
!  OUTPUT:
!
!    ptmix    Total mixing in potential temperature
!             equation (K*kg/(m**3*s)).
!
!  WORK ARRAYS:
!
!    h1       Turbulent heat flux, a local array.
!    h2       Turbulent heat flux, a local array.
!    h3       Turbulent heat flux, a local array.
!    tem1     Temporary work array
!    tem2     Temporary work array
!    tem3     Temporary work array
!
!   (These arrays are defined and used locally (i.e. inside this
!    subroutine), they may also be passed into routines called by
!    this one. Exiting the call to this subroutine, these temporary
!    work arrays may be used for other purposes, and therefore their
!    contents may be overwritten. Please examine the usage of work
!    arrays before you alter the code.)
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
  INTEGER :: nx,ny,nz          ! Number of grid points in 3 directions
!
  REAL :: ptprt (nx,ny,nz)     ! Perturbation potential temperature (K)
  REAL :: ptbar (nx,ny,nz)     ! Base state potential temperature (K)
  REAL :: rhostr(nx,ny,nz)     ! Base state air density times j3 (kg/m**3)
  REAL :: kmh   (nx,ny,nz)     ! Horizontal turb. mixing coef. for
                               ! momentum. ( m**2/s )
  REAL :: kmv   (nx,ny,nz)     ! Vertical turb. mixing coef. for
                               ! momentum. ( m**2/s )
  REAL :: rprntl(nx,ny,nz)     ! Reciprocal of Prandtl number

  REAL :: usflx(nx,ny)         ! Surface flux of u momentum
  REAL :: vsflx(nx,ny)         ! surface flux of v momentum
  REAL :: ptsflx(nx,ny)        ! surface flux of heat (K*kg/(m**2*s))
  REAL :: pbldpth(nx,ny)       ! Planetary boundary layer depth (m)

  REAL :: x     (nx)           ! The x-coord. of the physical and
                               ! computational grid. Defined at u-point.
  REAL :: y     (ny)           ! The y-coord. of the physical and
                               ! computational grid. Defined at v-point.
  REAL :: z     (nz)           ! The z-coord. of the computational grid.
                               ! Defined at w-point on the staggered grid.
  REAL :: zp    (nx,ny,nz)     ! The physical height coordinate defined at
                               ! w-point of the staggered grid.

  REAL :: mapfct(nx,ny,8)      ! Map factors at scalar points

  REAL :: j1    (nx,ny,nz)     ! Coordinate transform Jacobian defined as
                               ! - d( zp )/d( x ).
  REAL :: j2    (nx,ny,nz)     ! Coordinate transform Jacobian defined as
                               ! - d( zp )/d( y ).
  REAL :: j3    (nx,ny,nz)     ! Coordinate transform Jacobian defined as
                               ! d( zp )/d( z ).
  REAL :: aj3x  (nx,ny,nz)     ! Coordinate transformation Jacobian defined
                               ! as d( zp )/d( z ) AVERAGED IN THE X-DIR.
  REAL :: aj3y  (nx,ny,nz)     ! Coordinate transformation Jacobian defined
                               ! as d( zp )/d( z ) AVERAGED IN THE Y-DIR.
  REAL :: j3inv (nx,ny,nz)     ! Inverse of j3
  REAL :: ptsfc  (nx,ny)       ! Ground surface potential temperature (K)
!
  REAL :: ptmix (nx,ny,nz)     ! Turbulent mixing on potential
                               ! temperature (K*kg/(m**3 *s))
!
  REAL :: h1    (nx,ny,nz)     ! Turbulent heat flux
  REAL :: h2    (nx,ny,nz)     ! Turbulent heat flux
  REAL :: h3    (nx,ny,nz)     ! Turbulent heat flux
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
  INTEGER :: i,j,k,kt
  REAL    :: mdpth                ! Match layer depth (m)
  REAL    :: tem,temb,dtqflxdis1
  INTEGER :: cdiszero
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'grid.inc'          ! Grid & map parameters.
  INCLUDE 'phycst.inc'
  INCLUDE 'bndry.inc'
!
!-----------------------------------------------------------------------


!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


  IF(tmixopt == 0) THEN

    DO k=1,nz
      DO j=1,ny
        DO i=1,nx
          h3(i,j,k)=0.0
        END DO
      END DO
    END DO

  ELSE

!-----------------------------------------------------------------------
!
!  Compute turbulent heat fluxes H1, H2 and H3.
!
!-----------------------------------------------------------------------

    DO k=1,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          tem1(i,j,k)=ptbar(i,j,k)+ptprt(i,j,k)
        END DO
      END DO
    END DO

    CALL trbflxs(nx,ny,nz,tem1,rhostr,kmh,kmv,rprntl,                   &
                 x,y,z,zp, j1,j2,j3,j3inv,mapfct,                       &
                 h1,h2,h3, tem2,tem3)

  END IF

!-----------------------------------------------------------------------
!
!  The surface heat flux depends on the option chosen:
!
!  When sfcphy = 0, the internally calculated (from SGS turbulence
!  parameterization) surface fluxes are used.
!  Otherwise, pre-calculated fluxes ptsflx is used.
!
!-----------------------------------------------------------------------
!
  IF( (landwtr == 0.AND.cdhlnd == 0.0).OR.                              &
        (landwtr == 1.AND.cdhlnd == 0.0.AND.cdhwtr == 0.0) ) THEN
    cdiszero = 1
  ELSE
    cdiszero = 0
  END IF

  IF( (sfcphy == 1.OR.sfcphy == 3) .AND. cdiszero == 1) GO TO 111

  IF( sfcphy /= 0 ) THEN
!
    IF ( (sflxdis == 0) .OR. (sflxdis == 1 .AND.                        &
           (sfcphy == 1.OR.sfcphy == 3).AND.cdiszero == 1) ) THEN
!
!-----------------------------------------------------------------------
!
!  Heat flux h3 at the surface = ptsflx.
!
!-----------------------------------------------------------------------
!
      DO j=2,ny-2
        DO i=2,nx-2
          h3(i,j,2) = ptsflx(i,j)
        END DO
      END DO

    ELSE IF (sflxdis == 1 .OR. sflxdis == 2)  THEN

      DO j=2,ny-2
        DO i=2,nx-2
          tem=0.833*pbldpth(i,j)
          temb = ptsflx(i,j)/tem
          h3(i,j,2) = ptsflx(i,j)
          IF(ptsfc(i,j) > ptprt(i,j,2)+ptbar(i,j,2))THEN
            DO k=3,nz-1
              IF ((zp(i,j,k)-zp(i,j,2)) <= tem)                         &
                  h3(i,j,k)=ptsflx(i,j)-(zp(i,j,k)-zp(i,j,2))*temb
            END DO
          END IF
        END DO
      END DO

    END IF

    IF (tqflxdis == 1)  THEN

      DO j=2,ny-2
        DO i=2,nx-2
          IF(pbldpth(i,j)+zp(i,j,2) > 0.51*(zp(i,j,3)+zp(i,j,2)))THEN
            dtqflxdis1 = MIN(dtqflxdis,pbldpth(i,j))
          ELSE
            dtqflxdis1 = dtqflxdis
          END IF
          tem = 1.0/dtqflxdis1

          temb = ptsflx(i,j)*tem
          h3(i,j,2) = ptsflx(i,j)
          DO k=3,nz-1
            IF ((zp(i,j,k)-zp(i,j,2)) <= dtqflxdis1)                    &
                h3(i,j,k)= h3(i,j,k) + ptsflx(i,j) *                    &
                       ( 1. - (zp(i,j,k)-zp(i,j,2)) *tem ) ** 2
          END DO
        END DO
      END DO

    ELSE IF (tqflxdis == 2)  THEN

      DO j=2,ny-2
        DO i=2,nx-2
          tem = usflx(i,j)*usflx(i,j)+vsflx(i,j)*vsflx(i,j)
          tem = (SQRT(tem))**3
          temb = ABS(ptsflx(i,j))
          mdpth = 0.254842*tem*(ptprt(i,j,2)                            &
                         +ptbar(i,j,2))/(temb+1.e-20)
          mdpth = mdpth*2.5
          tem = MAX (pbldpth(i,j), 100.)
          mdpth = MIN(mdpth,tem)
          mdpth = MAX(mdpth,(zp(i,j,3)-zp(i,j,2)))
          kt=1
          temb =1./(0.833*tem)
          IF(ptsfc(i,j) <= ptprt(i,j,2)+ptbar(i,j,2)) THEN
            kt = 2
            temb = 1./mdpth
          END IF
          h3(i,j,2) = ptsflx(i,j)
          DO k=3,nz-1
            IF ((zp(i,j,k)-zp(i,j,2)) <= mdpth)                         &
                h3(i,j,k)=ptsflx(i,j)*(1.-(zp(i,j,k)-zp(i,j,2))*temb)   &
                      **kt+h3(i,j,k)
          END DO
        END DO
      END DO

    END IF

  END IF

  111   CONTINUE
!
!-----------------------------------------------------------------------
!
!  Calculate the turbulent mixing term for the potential temperature
!
!  ptmix = difx(avgx(j3) * h1) + dify(avgy(j3) * h2) + difz(h3 +
!          avgx(avgz(h1) * j1) + avgy(avgz(h2) * j2))
!
!-----------------------------------------------------------------------

  IF( tmixopt == 0 ) THEN

    DO k=2,nz-2
      DO j=2,ny-2
        DO i=2,nx-2
          ptmix(i,j,k)=(h3(i,j,k+1)-h3(i,j,k))*dzinv
        END DO
      END DO
    END DO

  ELSE

    CALL smixtrm(nx,ny,nz,h1,h2,h3,rhostr,x,y,z,zp,mapfct,              &
                 j1,j2,j3,aj3x,aj3y,                                    &
                 ptmix,                                                 &
                 tem1,tem2,tem3)

  END IF
!
!-----------------------------------------------------------------------
!
!  Set the boundary conditions for the mixing term
!
!-----------------------------------------------------------------------
!
  DO k=1,nz-1
    DO i=1,nx-1
      ptmix(i,   1,k)=ptmix(i,2,k)
      ptmix(i,ny-1,k)=ptmix(i,ny-2,k)
    END DO
  END DO

  DO k=1,nz-1
    DO j=1,ny-1
      ptmix(   1,j,k)=ptmix(2,j,k)
      ptmix(nx-1,j,k)=ptmix(nx-2,j,k)
    END DO
  END DO

  RETURN
END SUBROUTINE tmixpt

!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE TMIXQV                     ######
!######                                                      ######
!######                     Developed by                     ######
!######  Center for the Analysis and Prediction of Storms    ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################

SUBROUTINE tmixqv(nx,ny,nz,qv,rhostr,kmh,kmv,rprntl,                    &
           qvsflx,pbldpth,                                              &
           x,y,z,zp,mapfct, j1,j2,j3,aj3x,aj3y,j3inv,                   &
           usflx,vsflx,ptsflx,ptsfc,qvsfc,ptbar,ptprt,                  &
           qvmix,                                                       &
           h1,h2,h3, tem1,tem2,tem3)

!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate the turbulent mixing term for the water vapor specific
!  humidity equation. This term is expressed in the form of a turbulent
!  flux of the quantity in question.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  10/10/91
!
!  MODIFICATION HISTORY:
!
!  5/05/92 (M. Xue)
!  Added full documentation.
!
!  6/1/92 (M. Xue and H. Jin)
!  Further facelift.
!
!  7/20/92 (M. Xue and D. Weber)
!  Terrain included.
!
!  6/16/93 (MX)
!  Break down into subroutine calls.
!
!  1/23/96 (Donghai Wang & Yuhe Liu)
!  Added the map projection factor.
!
!  3/16/96 (Ming Xue)
!  The model can now include the surface fluxes with the turbulent
!  mxing option off (tmixopt=0).
!
!  3/27/1998 (M. Xue)
!  Added option tqflxdis for quadratic distribution of heat
!  and moisture fluxes in depth dtqflxdis, which is typically
!  set to 200 m.
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
!    qv       Water vapor specific humidity at a given time level (kg/kg)
!    pbldpth  Planetary boundary layer depth (m)
!    rhostr   Base state air density times j3 (kg/m**3)
!
!    kmh      Horizontal turb. mixing coef. for momentum ( m**2/s )
!    kmv      Vertical turb. mixing coef. for momentum ( m**2/s )
!    rprntl   Reciprocal of Prandtl number
!
!    qvsflx  Surface flux of moisture (K*kg/(m**2*s))
!
!    x        x coordinate of grid points in physical/comp. space (m)
!    y        y coordinate of grid points in physical/comp. space (m)
!    z        z coordinate of grid points in computational space (m)
!    zp       Vertical coordinate of grid points in physical space (m)
!
!    mapfct   Map factors at scalar points
!
!    j1       Coordinate transform Jacobian -d(zp)/dx
!    j2       Coordinate transform Jacobian -d(zp)/dy
!    j3       Coordinate transform Jacobian  d(zp)/dz
!    aj3x     Avgx of the coordinate transformation Jacobian  d(zp)/dz
!    aj3y     Avgy of the coordinate transformation Jacobian  d(zp)/dz
!    j3inv    Inverse of j3
!
!  OUTPUT:
!
!    qvmix    Total mixing in water vapor equation (kg/(m**3 * s))
!
!  WORK ARRAYS:
!
!    h1       Turbulent flux of moisture. A local array.
!    h2       Turbulent flux of moisture. A local array.
!    h3       Turbulent flux of moisture. A local array.
!    tem1     Temporary work array.
!    tem2     Temporary work array.
!    tem3     Temporary work array.
!
!   (These arrays are defined and used locally (i.e. inside this
!    subroutine), they may also be passed into routines called by
!    this one. Exiting the call to this subroutine, these temporary
!    work arrays may be used for other purposes, and therefore their
!    contents may be overwritten. Please examine the usage of work
!    arrays before you alter the code.)
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
  INTEGER :: nx,ny,nz          ! Number of grid points in 3 directions
!
  REAL :: qv    (nx,ny,nz)     ! Water vapor specific humidity (kg/kg)
  REAL :: pbldpth(nx,ny)        ! Planetary boundary layer depth (m)
  REAL :: rhostr(nx,ny,nz)     ! Base state air density times j3 (kg/m**3)
  REAL :: kmh   (nx,ny,nz)     ! Horizontal turb. mixing coef. for
                               ! momentum. ( m**2/s )
  REAL :: kmv   (nx,ny,nz)     ! Vertical turb. mixing coef. for
                               ! momentum. ( m**2/s )
  REAL :: rprntl(nx,ny,nz)     ! Reciprocal of Prandtl number

  REAL :: qvsflx(nx,ny)        ! surface flux of moisture (kg/(m**2*s))

  REAL :: x     (nx)           ! The x-coord. of the physical and
                               ! computational grid. Defined at u-point.
  REAL :: y     (ny)           ! The y-coord. of the physical and
                               ! computational grid. Defined at v-point.
  REAL :: z     (nz)           ! The z-coord. of the computational grid.
                               ! Defined at w-point on the staggered grid.
  REAL :: zp    (nx,ny,nz)     ! The physical height coordinate defined at
                               ! w-point of the staggered grid.

  REAL :: mapfct(nx,ny,8)      ! Map factors at scalar points

  REAL :: j1    (nx,ny,nz)     ! Coordinate transform Jacobian defined as
                               ! - d( zp )/d( x ).
  REAL :: j2    (nx,ny,nz)     ! Coordinate transform Jacobian defined as
                               ! - d( zp )/d( y ).
  REAL :: j3    (nx,ny,nz)     ! Coordinate transform Jacobian defined as
                               ! d( zp )/d( z ).
  REAL :: aj3x  (nx,ny,nz)     ! Coordinate transformation Jacobian defined
                               ! as d( zp )/d( z ) AVERAGED IN THE X-DIR.
  REAL :: aj3y  (nx,ny,nz)     ! Coordinate transformation Jacobian defined
                               ! as d( zp )/d( z ) AVERAGED IN THE Y-DIR.
  REAL :: j3inv (nx,ny,nz)     ! Inverse of j3

  REAL :: ptbar (nx,ny,nz)     ! Base state potential temperature (K)
  REAL :: ptprt (nx,ny,nz)     ! Perturbation potential temperature
  REAL :: ptsfc(nx,ny)         ! Temperature at ground (K) (in top 1 cm layer)
  REAL :: qvsfc(nx,ny)         ! Effective qv at the surface (kg/kg)

  REAL :: usflx(nx,ny)         ! Surface flux of u momentum
  REAL :: vsflx(nx,ny)         ! surface flux of v momentum
  REAL :: ptsflx(nx,ny)        ! Surface heat flux (K*kg/(m**2*s))

  REAL :: qvmix (nx,ny,nz)     ! Turbulent mixing on water substance
                               ! kg/(m**3 *s)
!
  REAL :: h1    (nx,ny,nz)     ! Turbulent moisture flux.
  REAL :: h2    (nx,ny,nz)     ! Turbulent moisture flux.
  REAL :: h3    (nx,ny,nz)     ! Turbulent moisture flux.
  REAL :: tem1  (nx,ny,nz)     ! Temporary work array
  REAL :: tem2  (nx,ny,nz)     ! Temporary work array
  REAL :: tem3  (nx,ny,nz)     ! Temporary work array


!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k,kt
  REAL    :: mdpth                ! Match layer depth (m)
  REAL    :: tem,temb,dtqflxdis1
  INTEGER :: cdiszero
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'grid.inc'          ! Grid & map parameters.
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!-----------------------------------------------------------------------
!
!  Compute turbulent moisture fluxes H1, H2 and H3.
!
!-----------------------------------------------------------------------

  IF(tmixopt == 0) THEN

    DO k=1,nz
      DO j=1,ny
        DO i=1,nx
          h3(i,j,k)=0.0
        END DO
      END DO
    END DO

  ELSE

    CALL trbflxs(nx,ny,nz,qv,rhostr,kmh,kmv,rprntl,                     &
                 x,y,z,zp, j1,j2,j3,j3inv,mapfct,                       &
                 h1,h2,h3, tem1,tem2)

  END IF

!-----------------------------------------------------------------------
!
!  Set the surface moisture fluxes:
!
!  When sfcphy = 0, the internally calculated (from SGS turbulence
!  parameterization) surface fluxes are used.
!  Otherwise, pre-calculated fluxes qvsflx is used.
!
!-----------------------------------------------------------------------
!
  IF( (landwtr == 0.AND.cdqlnd == 0.0).OR.                              &
        (landwtr == 1.AND.cdqlnd == 0.0.AND.cdqwtr == 0.0) ) THEN
    cdiszero = 1
  ELSE
    cdiszero = 0
  END IF

  IF( (sfcphy == 1.OR.sfcphy == 3) .AND. cdiszero == 1) GO TO 111

  IF( sfcphy /= 0 ) THEN

    IF ( (sflxdis == 0) .OR. (sflxdis == 1 .AND.                        &
           (sfcphy == 1.OR.sfcphy == 3).AND.cdiszero == 1) ) THEN
!
!-----------------------------------------------------------------------
!
!  Moisture flux h3 at the surface = qvsflx.
!
!-----------------------------------------------------------------------
!
      DO j=2,ny-2
        DO i=2,nx-2
          h3(i,j,2) = qvsflx(i,j)
        END DO
      END DO

    ELSE IF (sflxdis == 1 .OR. sflxdis == 2)  THEN

      DO j=2,ny-2
        DO i=2,nx-2
          tem=0.833*pbldpth(i,j)
          temb = qvsflx(i,j)/tem
          h3(i,j,2) = qvsflx(i,j)
          IF(ptsfc(i,j) > ptprt(i,j,2)+ptbar(i,j,2))THEN
            DO k=3,nz-1
              IF ((zp(i,j,k)-zp(i,j,2)) <= tem)                         &
                  h3(i,j,k)=qvsflx(i,j)-(zp(i,j,k)-zp(i,j,2))*temb
            END DO
          END IF
        END DO
      END DO

    END IF

    IF (tqflxdis == 1)  THEN

      DO j=2,ny-2
        DO i=2,nx-2
          IF(pbldpth(i,j)+zp(i,j,2) > 0.51*(zp(i,j,3)+zp(i,j,2)))THEN
            dtqflxdis1 = MIN(dtqflxdis,pbldpth(i,j))
          ELSE
            dtqflxdis1 = dtqflxdis
          END IF
          tem = 1.0/dtqflxdis1

          temb = qvsflx(i,j)*tem
          h3(i,j,2) = qvsflx(i,j)
          DO k=3,nz-1
            IF((zp(i,j,k)-zp(i,j,2)) <= dtqflxdis1)                     &
                h3(i,j,k)= h3(i,j,k) + qvsflx(i,j) *                    &
                       ( 1. - (zp(i,j,k)-zp(i,j,2)) * tem ) ** 2
          END DO
        END DO
      END DO

    ELSE IF (tqflxdis == 2)  THEN

      DO j=2,ny-2
        DO i=2,nx-2
          tem = usflx(i,j)*usflx(i,j)+vsflx(i,j)*vsflx(i,j)
          tem = (SQRT(tem))**3
          temb = ABS(ptsflx(i,j))
          mdpth = 0.254842*tem*(ptprt(i,j,2)                            &
                         +ptbar(i,j,2))/(temb+1.e-20)

          mdpth = mdpth*2.5
          tem = MAX (pbldpth(i,j), 100.)
          mdpth = MIN(mdpth,tem)

          mdpth = MAX(mdpth,(zp(i,j,3)-zp(i,j,2)))
          kt=1
          temb =1./tem
          IF(ptsfc(i,j) <= ptprt(i,j,2)+ptbar(i,j,2)) THEN
            kt = 2
            temb = 1./mdpth
          END IF
          h3(i,j,2) = qvsflx(i,j)
          DO k=3,nz-1
            IF ((zp(i,j,k)-zp(i,j,2)) <= mdpth)                         &
                h3(i,j,k)=qvsflx(i,j)*(1.-(zp(i,j,k)-zp(i,j,2))*temb)   &
                      **kt+h3(i,j,k)
          END DO
        END DO
      END DO

    END IF

  END IF

  111   CONTINUE
!
!-----------------------------------------------------------------------
!
!  Calculate the mixing term for the water vapor mixing ratio (qv)
!
!  qvmix = difx(avgx(j3) * h1) + dify(avgy(j3) * h2) + difz(h3 +
!          avgx(avgz(h1) * j1) + avgy(avgz(h2) * j2))
!
!-----------------------------------------------------------------------

  IF( tmixopt == 0 ) THEN

    DO k=2,nz-2
      DO j=2,ny-2
        DO i=2,nx-2
          qvmix(i,j,k)=(h3(i,j,k+1)-h3(i,j,k))*dzinv
        END DO
      END DO
    END DO

  ELSE

    CALL smixtrm(nx,ny,nz,h1,h2,h3,rhostr,x,y,z,zp,mapfct,              &
                 j1,j2,j3,aj3x,aj3y,                                    &
                 qvmix,                                                 &
                 tem1,tem2,tem3)

  END IF
!
!-----------------------------------------------------------------------
!
!  Set the qvmix boundary values to the inside points.
!
!-----------------------------------------------------------------------
!

  DO k=1,nz-1
    DO i=1,nx-1
      qvmix(i,   1,k)=qvmix(i,2,k)
      qvmix(i,ny-1,k)=qvmix(i,ny-2,k)
    END DO
  END DO

  DO k=1,nz-1
    DO j=1,ny-1
      qvmix(   1,j,k)=qvmix(2,j,k)
      qvmix(nx-1,j,k)=qvmix(nx-2,j,k)
    END DO
  END DO

  RETURN
END SUBROUTINE tmixqv

!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE TMIXQ                      ######
!######                                                      ######
!######                     Developed by                     ######
!######  Center for the Analysis and Prediction of Storms    ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE tmixq(nx,ny,nz, q,rhostr,kmh,kmv,rprntl,                     &
           x,y,z,zp,mapfct, j1,j2,j3,aj3x,aj3y,j3inv,                   &
           qmix,                                                        &
           h1,h2,h3, tem1,tem2,tem3)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate the turbulent mixing term for a water substance equation.
!  This term is expressed in the form of a turbulent flux of the water
!  quantity in question.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  10/10/91
!
!  MODIFICATION HISTORY:
!
!  5/05/92 (M. Xue)
!  Added full documentation.
!
!  6/1/92 (M. Xue and H. Jin)
!  Further facelift.
!
!  7/20/92 (M. Xue and D. Weber)
!  Terrain included.
!
!  6/16/93 (MX)
!  Break down into subroutine calls.
!
!  1/23/96 (Donghai Wang & Yuhe Liu)
!  Added the map projection factor.
!
!  11/06/97 (D. Weber)
!  Added three additional levels to the mapfct array.  The three
!  levels (4,5,6) represent the inverse of the first three in order.
!  The inverse map factors are computed to improve efficiency.
!
!  9/18/98 (D. Weber)
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
!    q        mixing ratio of one of the water/ice quantities at
!             a given time level (kg/kg)
!
!    rhostr   Base state air density times j3 (kg/m**3)
!
!    kmh      Horizontal turb. mixing coef. for momentum ( m**2/s )
!    kmv      Vertical turb. mixing coef. for momentum ( m**2/s )
!    rprntl   Reciprocal of Prandtl number
!
!    x        x coordinate of grid points in physical/comp. space (m)
!    y        y coordinate of grid points in physical/comp. space (m)
!    z        z coordinate of grid points in computational space (m)
!    zp       Vertical coordinate of grid points in physical space (m)
!
!    mapfct   Map factors at scalar, u, and v points
!
!    j1       Coordinate transform Jacobian -d(zp)/dx
!    j2       Coordinate transform Jacobian -d(zp)/dy
!    j3       Coordinate transform Jacobian  d(zp)/dz
!    aj3x     Avgx of the coordinate transformation Jacobian  d(zp)/dz
!    aj3y     Avgy of the coordinate transformation Jacobian  d(zp)/dz
!    j3inv    Inverse of j3
!
!  OUTPUT:
!
!    qmix     Total mixing in water/ice substance
!             equation (kg/(m**3 * s))
!
!  WORK ARRAYS:
!
!    h1       Turbulent flux a water quantity, a local array.
!    h2       Turbulent flux a water quantity, a local array.
!    h3       Turbulent flux a water quantity, a local array.
!    tem1     Temporary work array.
!    tem2     Temporary work array.
!    tem3     Temporary work array.
!
!   (These arrays are defined and used locally (i.e. inside this
!    subroutine), they may also be passed into routines called by
!    this one. Exiting the call to this subroutine, these temporary
!    work arrays may be used for other purposes, and therefore their
!    contents may be overwritten. Please examine the usage of work
!    arrays before you alter the code.)
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
  INTEGER :: nx,ny,nz          ! Number of grid points in 3 directions
!
  REAL :: q     (nx,ny,nz)     ! Water/ice mixing ratio (kg/kg)
  REAL :: rhostr(nx,ny,nz)     ! Base state air density times j3 (kg/m**3)
  REAL :: kmh   (nx,ny,nz)     ! Horizontal turb. mixing coef. for
                               ! momentum. ( m**2/s )
  REAL :: kmv   (nx,ny,nz)     ! Vertical turb. mixing coef. for
                               ! momentum. ( m**2/s )
  REAL :: rprntl(nx,ny,nz)     ! Reciprocal of Prandtl number

  REAL :: x     (nx)           ! The x-coord. of the physical and
                               ! computational grid. Defined at u-point.
  REAL :: y     (ny)           ! The y-coord. of the physical and
                               ! computational grid. Defined at v-point.
  REAL :: z     (nz)           ! The z-coord. of the computational grid.
                               ! Defined at w-point on the staggered grid.
  REAL :: zp    (nx,ny,nz)     ! The physical height coordinate defined at
                               ! w-point of the staggered grid.
  REAL :: mapfct(nx,ny,8)      ! Map factors at scalar, u, and v points

  REAL :: j1    (nx,ny,nz)     ! Coordinate transform Jacobian defined as
                               ! - d( zp )/d( x ).
  REAL :: j2    (nx,ny,nz)     ! Coordinate transform Jacobian defined as
                               ! - d( zp )/d( y ).
  REAL :: j3    (nx,ny,nz)     ! Coordinate transform Jacobian defined as
                               ! d( zp )/d( z ).
  REAL :: aj3x  (nx,ny,nz)     ! Coordinate transformation Jacobian defined
                               ! as d( zp )/d( z ) AVERAGED IN THE X-DIR.
  REAL :: aj3y  (nx,ny,nz)     ! Coordinate transformation Jacobian defined
                               ! as d( zp )/d( z ) AVERAGED IN THE Y-DIR.
  REAL :: j3inv (nx,ny,nz)     ! Inverse of j3
!
  REAL :: qmix  (nx,ny,nz)     ! Turbulent mixing on water/ice
                               ! substance (kg/(m**2 *s))
!
  REAL :: h1    (nx,ny,nz)     ! Turbulent flux of a water quantity
  REAL :: h2    (nx,ny,nz)     ! Turbulent flux of a water quantity
  REAL :: h3    (nx,ny,nz)     ! Turbulent flux of a water quantity
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
!
!-----------------------------------------------------------------------
!
!  Include files:
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


!-----------------------------------------------------------------------
!
!  Compute turbulent water fluxes H1, H2 and H3.
!
!-----------------------------------------------------------------------

    CALL trbflxs(nx,ny,nz,q,rhostr,kmh,kmv,rprntl,                      &
                 x,y,z,zp, j1,j2,j3,j3inv,mapfct,                       &
                 h1,h2,h3, tem1,tem2)

!-----------------------------------------------------------------------
!
!  Calculate the mixing term for the water mixing ratio (q)
!
!  qmix = difx(avgx(j3) * h1) + dify(avgy(j3) * h2) + difz(h3 +
!         avgx(avgz(h1) * j1) + avgy(avgz(h2) * j2))
!
!-----------------------------------------------------------------------


    CALL smixtrm(nx,ny,nz,h1,h2,h3,rhostr,x,y,z,zp,mapfct,              &
                 j1,j2,j3,aj3x,aj3y,                                    &
                 qmix,                                                  &
                 tem1,tem2,tem3)

!
!-----------------------------------------------------------------------
!
!  Set the boundary conditions for the mixing term
!
!-----------------------------------------------------------------------
!
  DO k=1,nz-1
    DO i=1,nx-1
      qmix(i,   1,k)=qmix(i,2,k)
      qmix(i,ny-1,k)=qmix(i,ny-2,k)
    END DO
  END DO

  DO k=1,nz-1
    DO j=1,ny-1
      qmix(   1,j,k)=qmix(2,j,k)
      qmix(nx-1,j,k)=qmix(nx-2,j,k)
    END DO
  END DO

  RETURN
END SUBROUTINE tmixq
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE TRBFLXS                    ######
!######                                                      ######
!######                     Developed by                     ######
!######  Center for the Analysis and Prediction of Storms    ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE trbflxs(nx,ny,nz,s,rhostr,kmh,kmv,rprntl,x,y,z,zp,           &
           j1,j2,j3,j3inv,mapfct,                                       &
           h1,h2,h3, tem1,tem2)
!
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate the turbulent fluxes in x, y and z direction for a
!  scalar quantity 's'
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  6/16/1993
!
!  MODIFICATION HISTORY:
!
!  4/1/96 (Donghai Wang, X. Song and M. Xue)
!  Added a time average weighting coefficient for implicit treatment
!  of vertical mixing.
!
!  8/27/98 (D. Weber)
!  Modified do loop structure to improve code efficiency via:
!  -merging existing loops together
!  -removing and merging operators
!  -switched to DO-ENDDO structure for f90 conversion
!  -added mapfct to h1,h2 computation (previously performed in
!   smixtrm)
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!    s        A scalar variable.
!    rhostr   Base state air density times j3 (kg/m**3)
!
!    kmh      Horizontal turb. mixing coef. for momentum ( m**2/s )
!    kmv      Vertical turb. mixing coef. for momentum ( m**2/s )
!    rprntl   Reciprocal of Prandtl number
!
!    x        x coordinate of grid points in physical/comp. space (m)
!    y        y coordinate of grid points in physical/comp. space (m)
!    z        z coordinate of grid points in computational space (m)
!    zp       Vertical coordinate of grid points in physical space (m)
!
!    j1       Coordinate transform Jacobian -d(zp)/dx
!    j2       Coordinate transform Jacobian -d(zp)/dy
!    j3       Coordinate transform Jacobian  d(zp)/dz
!    j3inv    Inverse of j3
!
!    mapfct   Map factors at scalar, u, and v points
!
!  OUTPUT:
!
!    h1       Turbulent flux in x direction.
!    h2       Turbulent flux in y direction.
!    h3       Turbulent flux in z direction.
!
!  WORK ARRAYS:
!
!    tem1     Temporary work array
!    tem2     Temporary work array
!
!   (These arrays are defined and used locally (i.e. inside this
!    subroutine), they may also be passed into routines called by
!    this one. Exiting the call to this subroutine, these temporary
!    work arrays may be used for other purposes, and therefore their
!    contents may be overwritten. Please examine the usage of work
!    arrays before you alter the code.)
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
  INTEGER :: nx,ny,nz          ! Number of grid points in 3 directions
!
  REAL :: s     (nx,ny,nz)     ! A scalar variable.
  REAL :: rhostr(nx,ny,nz)     ! Base state air density times j3 (kg/m**3)
  REAL :: kmh   (nx,ny,nz)     ! Horizontal turb. mixing coef. for
                               ! momentum. ( m**2/s )
  REAL :: kmv   (nx,ny,nz)     ! Vertical turb. mixing coef. for
                               ! momentum. ( m**2/s )
  REAL :: rprntl(nx,ny,nz)     ! Reciprocal of Prandtl number

  REAL :: x     (nx)           ! The x-coord. of the physical and
                               ! computational grid. Defined at u-point.
  REAL :: y     (ny)           ! The y-coord. of the physical and
                               ! computational grid. Defined at v-point.
  REAL :: z     (nz)           ! The z-coord. of the computational grid.
                               ! Defined at w-point on the staggered grid.
  REAL :: zp    (nx,ny,nz)     ! The physical height coordinate defined at
                               ! w-point of the staggered grid.

  REAL :: j1    (nx,ny,nz)     ! Coordinate transform Jacobian defined as
                               ! - d( zp )/d( x ).
  REAL :: j2    (nx,ny,nz)     ! Coordinate transform Jacobian defined as
                               ! - d( zp )/d( y ).
  REAL :: j3    (nx,ny,nz)     ! Coordinate transform Jacobian defined as
                               ! d( zp )/d( z ).
  REAL :: j3inv (nx,ny,nz)     ! Inverse of j3

  REAL :: mapfct(nx,ny,8)      ! Map factors at scalar, u, and v points

!
  REAL :: h1    (nx,ny,nz)     ! Turbulent flux in x direction
  REAL :: h2    (nx,ny,nz)     ! Turbulent flux in y direction
  REAL :: h3    (nx,ny,nz)     ! Turbulent flux in z direction

  REAL :: tem1  (nx,ny,nz)     ! Temporary work array
  REAL :: tem2  (nx,ny,nz)     ! Temporary work array
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k
  REAL :: del2inv
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'grid.inc'          ! Grid & map parameters.
  INCLUDE 'bndry.inc'
  INCLUDE 'mp.inc'            ! Message passing parameters.
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  del2inv = 1./(dx*dy)

!-----------------------------------------------------------------------
!
!  Compute turbulent flux in x direction (H1)
!
!  H1 = avgx(rhobar * Km / j3) * (difx(j3 * s ) +
!       difz(j1 * avgsw(avgx( s ))))
!
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!
!  Calculate the first term in the turbulent flux (H1)
!  difx(j3 * s)
!
!-----------------------------------------------------------------------

  IF( tmixvert == 0)THEN  ! compute the horizontal terms...

    DO k=1,nz-1    ! difx(j3 * s)
      DO j=2,ny-2
        DO i=2,nx-1
          h1(i,j,k) = dxinv*(j3(i,j,k)*s(i,j,k)                         &
                           - j3(i-1,j,k)*s(i-1,j,k))
        END DO
      END DO
    END DO

!-----------------------------------------------------------------------
!
!  Calculate the second term in the turbulent flux (H1)
!  difz(j1 * avgsw(avgx(s))). When ternopt=0, this term is zero.
!
!-----------------------------------------------------------------------

    IF( ternopt /= 0 ) THEN

      DO k=2,nz-1    ! avgsw(avgx(s))
        DO j=2,ny-2
          DO i=2,nx-1
            tem1(i,j,k) = 0.25*(s(i,j,k-1)+s(i-1,j,k-1)                 &
                               +s(i,j,k)  +s(i-1,j,k))
          END DO
        END DO
      END DO

      CALL acct_interrupt(bc_acct)
      CALL bcsw(nx,ny,nz,2,nx-1,2,ny-2,tbc,bbc,tem1)
      CALL acct_stop_inter

      DO k=1,nz-1
        DO j=2,ny-2
          DO i=2,nx-1
            h1(i,j,k) = h1(i,j,k) + dzinv*(j1(i,j,k+1)*tem1(i,j,k+1)    &
                                       -   j1(i,j,k)  *tem1(i,j,k))
          END DO
        END DO
      END DO

    END IF        ! end of terrain if block for H1....

!-----------------------------------------------------------------------
!
!  Calculate turbulent flux in Y direction (H2)
!
!  H2= avgy(rhobar*kh/j3)*(dify(j3 * s) + difz(j2 * avgsw(avgy(s))))
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!  Calculate the first term in the turbulent flux (H2)
!  dify(j3 * s)
!
!-----------------------------------------------------------------------

    DO k=1,nz-1    ! dify(j3 * s)
      DO j=2,ny-1
        DO i=1,nx-1
          h2(i,j,k) = dyinv*(j3(i,j,k)*s(i,j,k)                         &
                           - j3(i,j-1,k)*s(i,j-1,k))
        END DO
      END DO
    END DO

!-----------------------------------------------------------------------
!
!  Calculate the second term in the turbulent flux (H2)
!  difz(j2 * avgsw(avgy(s))). When ternopt=0, this term is zero.
!
!-----------------------------------------------------------------------

    IF( ternopt /= 0 ) THEN

      DO k=2,nz-1    ! avgsw(avgy(s))
        DO j=2,ny-1
          DO i=2,nx-2
            tem1(i,j,k) = 0.25*(s(i,j,k-1)+s(i,j-1,k-1)                 &
                               +s(i,j,k)  +s(i,j-1,k))
          END DO
        END DO
      END DO

      CALL acct_interrupt(bc_acct)
      CALL bcsw(nx,ny,nz,2,nx-2,2,ny-1,tbc,bbc,tem1)
      CALL acct_stop_inter

      DO k=1,nz-1
        DO j=2,ny-1
          DO i=2,nx-2
            h2(i,j,k) = h2(i,j,k) + dzinv*(j2(i,j,k+1)*tem1(i,j,k+1)    &
                                       -   j2(i,j,k)  *tem1(i,j,k))
          END DO
        END DO
      END DO

    END IF     ! end of terrain if block for H2....


!-----------------------------------------------------------------------
!
!  Calculate the coefficient for H1 and H2...
!  for H1  avgx(rhobar*kh/j3) = avgx(rhostr*kh/(j3*j3))
!  for H2  avgy(rhobar*kh/j3) = avgy(rhostr*kh/(j3*j3))
!
!-----------------------------------------------------------------------

    DO k=1,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          tem1(i,j,k) = rhostr(i,j,k)*kmh(i,j,k)*rprntl(i,j,k)          &
                      * (j3inv(i,j,k)*j3inv(i,j,k))
        END DO
      END DO
    END DO

    DO k=1,nz-1
      DO j=2,ny-2
        DO i=2,nx-1
          h1(i,j,k)  = 0.5*(tem1(i,j,k)+tem1(i-1,j,k))                  &
                             *h1(i,j,k)*mapfct(i,j,2)
        END DO
      END DO
    END DO     ! H1 is complete....

    DO k=1,nz-1
      DO j=2,ny-1
        DO i=2,nx-2
          h2(i,j,k) = 0.5*(tem1(i,j,k)+tem1(i,j-1,k))                   &
                          *h2(i,j,k)*mapfct(i,j,3)
        END DO
      END DO
    END DO     ! H2 is complete....

  END IF       !  end of tmixvert if block....

!-----------------------------------------------------------------------
!
!  Compute turbulent flux in z direction (H3).
!  H3= avgz (rhobar*kh/j3) * difz(s)
!
!  First calculate H3 = difz(s)
!
!-----------------------------------------------------------------------

  DO k=2,nz-1
    DO j=2,ny-2
      DO i=2,nx-2
        h3(i,j,k) = dzinv*(s(i,j,k)-s(i,j,k-1))
      END DO
    END DO
  END DO

!-----------------------------------------------------------------------
!
!  Calculate the coefficient for H3
!
!  avgz(rhobar*kh/j3)
!
!-----------------------------------------------------------------------

  DO k=1,nz-1
    DO j=2,ny-2
      DO i=2,nx-2
        tem1(i,j,k) = alfcoef * rhostr(i,j,k)*kmv(i,j,k)*rprntl(i,j,k)  &
                    * (j3inv(i,j,k)*j3inv(i,j,k))
      END DO
    END DO
  END DO

  DO k=2,nz-1
    DO j=2,ny-2
      DO i=2,nx-2
        h3(i,j,k) = 0.5*(tem1(i,j,k)+tem1(i,j,k-1))*h3(i,j,k)
      END DO
    END DO
  END DO

  RETURN
END SUBROUTINE trbflxs
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE SMIXTRM                    ######
!######                                                      ######
!######                     Developed by                     ######
!######  Center for the Analysis and Prediction of Storms    ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE smixtrm(nx,ny,nz,h1,h2,h3,rhostr,                            &
           x,y,z,zp, mapfct, j1,j2,j3,aj3x,aj3y,                        &
           smix,                                                        &
           tem1,tem2,tem3)
!
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate the turbulent mixing term for a scalar from the
!  turbulent fluxes h1, h2 and h3.
!
!  smix = difx(avgx(j3) * h1) + dify(avgy(j3) * h2) +
!         difz(h3 + avgx(avgz(h1) * j1) + avgy(avgz(h2) * j2))
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue & Dan Weber
!  6/16/93
!
!  MODIFICATION HISTORY:
!
!  1/23/96 (Donghai Wang & Yuhe Liu)
!  Added the map projection factor.
!
!  11/06/97 (D. Weber)
!  Added three additional levels to the mapfct array.  The three
!  levels (4,5,6) represent the inverse of the first three in order.
!  The inverse map factors are computed to improve efficiency.
!
!  9/18/98 (D. Weber)
!  Modified do loop structure to improve code efficiency via:
!  -merging existing loops together
!  -removing and merging operators
!  -switched to DO  ENDDO structure for f90 conversion
!  -added tmixvert option for computing the vertical turbulent
!   mixing terms only.
!  -added arrays aj3x,y.
!
!  2/24/1999 (M.Xue)
!  Set contributions of H1 and H2 to fluxes through lower boundary
!  to zero. They are zero anyway in the absence of terrain.
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!    h1       Turbulent heat flux, a local array.
!    h2       Turbulent heat flux, a local array.
!    h3       Turbulent heat flux, a local array.
!
!    rhostr   Base state air density times j3 (kg/m**3)
!
!    x        x coordinate of grid points in physical/comp. space (m)
!    y        y coordinate of grid points in physical/comp. space (m)
!    z        z coordinate of grid points in computational space (m)
!    zp       Vertical coordinate of grid points in physical space (m)
!
!    mapfct   Map factors at scalar, u, and v points
!
!    j1       Coordinate transform Jacobian -d(zp)/dx
!    j2       Coordinate transform Jacobian -d(zp)/dy
!    j3       Coordinate transform Jacobian  d(zp)/dz
!    aj3x     Avgx of the coordinate transformation Jacobian  d(zp)/dz
!    aj3y     Avgy of the coordinate transformation Jacobian  d(zp)/dz
!
!  OUTPUT:
!
!    smix     The mixing term for a scalar calc. from its turbulent fluxes
!
!  WORK ARRAYS:
!
!    tem1     Temporary work array
!    tem2     Temporary work array
!    tem3     Temporary work array
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
  INTEGER :: nx,ny,nz          ! Number of grid points in 3 directions
!
  REAL :: h1    (nx,ny,nz)     ! Turbulent flux in x direction
  REAL :: h2    (nx,ny,nz)     ! Turbulent flux in y direction
  REAL :: h3    (nx,ny,nz)     ! Turbulent flux in z direction

  REAL :: rhostr(nx,ny,nz)     ! Base state air density times j3 (kg/m**3)
                               ! momentum. ( m**2/s )

  REAL :: x     (nx)           ! The x-coord. of the physical and
                               ! computational grid. Defined at u-point.
  REAL :: y     (ny)           ! The y-coord. of the physical and
                               ! computational grid. Defined at v-point.
  REAL :: z     (nz)           ! The z-coord. of the computational grid.
                               ! Defined at w-point on the staggered grid.
  REAL :: zp    (nx,ny,nz)     ! The physical height coordinate defined at
                               ! w-point of the staggered grid.

  REAL :: mapfct(nx,ny,8)      ! Map factors at scalar, u, and v points

  REAL :: j1    (nx,ny,nz)     ! Coordinate transform Jacobian defined as
                               ! - d( zp )/d( x ).
  REAL :: j2    (nx,ny,nz)     ! Coordinate transform Jacobian defined as
                               ! - d( zp )/d( y ).
  REAL :: j3    (nx,ny,nz)     ! Coordinate transform Jacobian defined as
                               ! d( zp )/d( z ).
  REAL :: aj3x  (nx,ny,nz)     ! Coordinate transformation Jacobian defined
                               ! as d( zp )/d( z ) AVERAGED IN THE X-DIR.
  REAL :: aj3y  (nx,ny,nz)     ! Coordinate transformation Jacobian defined
                               ! as d( zp )/d( z ) AVERAGED IN THE Y-DIR.
!
  REAL :: smix  (nx,ny,nz)     ! Turbulent mixing for a scalar
!
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
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'grid.inc'          ! Grid & map parameters.
!
!-----------------------------------------------------------------------

!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


!-----------------------------------------------------------------------
!
!  Calculate the turbulent mixing term for a scalar as the divergence
!  of its turbluent fluxes:
!
!  smix = difx(avgx(j3) * h1) + dify(avgy(j3) * h2) + difz(h3 +
!         avgx(avgz(h1) * j1) + avgy(avgz(h2) * j2))
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!  Computing the vertical term first - difz(h3).
!
!-----------------------------------------------------------------------

  DO k=2,nz-2   ! add in the vertical h3 term.....
    DO j=2,ny-2
      DO i=2,nx-2
        smix(i,j,k) = dzinv*(h3(i,j,k+1)-h3(i,j,k))
      END DO
    END DO
  END DO

  IF( tmixvert == 0)THEN    ! compute the horizontal terms.

!-----------------------------------------------------------------------
!
!  Calculate term difx(avgx(j3) * h1) and dify(avgy(j3) * h2)
!  and difz(h3)  (non-terrain terms)
!
!-----------------------------------------------------------------------

    DO k=2,nz-2   ! add the two terms to smix...
      DO j=2,ny-2
        DO i=2,nx-2
          smix(i,j,k) =smix(i,j,k) + mapfct(i,j,1)*                     &
                        (dxinv*( h1(i+1,j,k)*aj3x(i+1,j,k)-             &
                                 h1(i,j,k)*aj3x(i,j,k))                 &
                        +dyinv*( h2(i,j+1,k)*aj3y(i,j+1,k)              &
                                -h2(i,j,k)*aj3y(i,j,k)))
        END DO
      END DO
    END DO     ! smix w/o terrain is complete.

!
!  At this point, smix = difx(avgx(j3) * h1) + dify(avgy(j3) * h2)
!

!-----------------------------------------------------------------------
!
!  Calculate the third term in the smix relation.
!  difz(h3 + avgx(avgz(h1) * j1) + avgy(avgz(h2) * j2)) ! old
!  difz( avgx(avgz(h1) * j1) + avgy(avgz(h2) * j2))   ! new
!
!-----------------------------------------------------------------------

    IF( ternopt /= 0 ) THEN

!-----------------------------------------------------------------------
!
!  Calculate the second and third parts of the third term.
!  mapfct( avgx(avgz(h1) * j1) + avgy(avgz(h2) * j2)) )
!
!-----------------------------------------------------------------------

      DO k=3,nz-1   ! avgz(h1) * j1 and avgz(h2) * j2
        DO j=2,ny-1
          DO i=2,nx-1
            tem1(i,j,k) = 0.5*j1(i,j,k)*(h1(i,j,k)+h1(i,j,k-1))
            tem2(i,j,k) = 0.5*j2(i,j,k)*(h2(i,j,k)+h2(i,j,k-1))
          END DO
        END DO
      END DO

      DO k=3,nz-1
        DO j=2,ny-2
          DO i=2,nx-2
            tem3(i,j,k) =   mapfct(i,j,1)*0.5*                          &
                               ((tem1(i+1,j,k)+tem1(i,j,k))             &
                               +(tem2(i,j+1,k)+tem2(i,j,k)))
          END DO
        END DO
      END DO

      DO j=2,ny-2
        DO i=2,nx-2
          tem3(i,j,2) = 0.0  ! This term should not contribute more to the
                             ! flux through the lower boundary, which is
                             ! stored in H3.
        END DO
      END DO

      DO k=2,nz-2    ! add the terrain terms to smix...
        DO j=2,ny-2
          DO i=2,nx-2
            smix(i,j,k) =smix(i,j,k)+dzinv*(tem3(i,j,k+1)-tem3(i,j,k))
          END DO
        END DO
      END DO         ! smix is complete.....


    END IF        ! end of terrain if block....

  END IF      ! end of tmixvert if block....

  RETURN
END SUBROUTINE smixtrm

!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE UMIXTRM                    ######
!######                                                      ######
!######                     Developed by                     ######
!######  Center for the Analysis and Prediction of Storms    ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE umixtrm(nx,ny,nz,mapfct,j1,j2,j3,aj3y,                       &
           tau11,tau12,tau13,                                           &
           umix,tem1,tem2)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate the turbulent mixing term in u equation as the divergence
!  of the turbulent momentum fluxes.
!
!  j3 * DIV U = difx (j3 * tau11) + dify (avgx (avgsv(j3)) * tau12)
!               + difz (tau13 + (j1 * avgz (avgx (tau11))) + (avgy
!               ((avgx (j2)) * (avgz (tau12))))
!
!-----------------------------------------------------------------------
!
!  AUTHOR: D. Weber & M. Xue
!  6/29/92.
!
!  MODIFICATION HISTORY:
!
!  1/23/96 (Donghai Wang & Yuhe Liu)
!  Added the map projection factor.
!
!  8/18/98 (D. Weber)
!  Modified do loop structure to improve code efficiency via:
!  -merging existing loops together
!  -removing and merging operators
!  -switched to DO  ENDDO structure for f90 conversion
!  -added tmixvert option to compute only vertical t-mixing terms.
!  -added array aj3y.
!
!  2/24/1999 (M.Xue)
!  Set contributions of tau11 and tau12 to momentum flux through
!  lower boundary to zero. They are zero anyway in the absence of terrain.
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!    mapfct   Map factors at u points
!
!    j1       Coordinate transform Jacobian -d(zp)/dx
!    j2       Coordinate transform Jacobian -d(zp)/dy
!    j3       Coordinate transform Jacobian  d(zp)/dz
!    aj3y     Avgy of the coordinate transformation Jacobian  d(zp)/dz
!
!    tau11    Deformation tensor component (1/s)
!    tau12    Deformation tensor component (1/s)
!    tau13    Deformation tensor component (1/s)
!
!
!  OUTPUT:
!
!    umix     The divergence of u mixing term
!
!  WORKING ARRAYS:
!
!    tem1     Temporary working array.
!    tem2     Temporary working array.
!
!   (These arrays are defined and used locally (i.e. inside this
!    subroutine), they may also be passed into routines called by
!    this one. Exiting the call to this subroutine, these temporary
!    working arrays may be used for other purposes therefore their
!    contents overwritten. Please examine the usage of working arrays
!    before you alter the code.)
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!
!  Variable Declarations
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
  INTEGER :: nx,ny,nz          ! Number of grid points in 3 directions

  REAL :: mapfct(nx,ny)        ! Map factors at u points

  REAL :: j1    (nx,ny,nz)     ! Coordinate transformation Jacobian
                               ! - d(zp)/dx
  REAL :: j2    (nx,ny,nz)     ! Coordinate transformation Jacobian
                               ! - d(zp)/dy
  REAL :: j3    (nx,ny,nz)     ! Coordinate transformation Jacobian
                               !   d(zp)/dz
  REAL :: aj3y  (nx,ny,nz)     ! Coordinate transformation Jacobian defined
                               ! as d( zp )/d( z ) AVERAGED IN THE Y-DIR.
!
  REAL :: tau11 (nx,ny,nz)     ! Deformation stress component.
  REAL :: tau12 (nx,ny,nz)     ! Deformation stress component.
  REAL :: tau13 (nx,ny,nz)     ! Deformation stress component.
!
  REAL :: umix  (nx,ny,nz)     ! Divergence of u mixing term
!
  REAL :: tem1  (nx,ny,nz)     ! Temporary working array.
  REAL :: tem2  (nx,ny,nz)     ! Temporary working array.

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
  INCLUDE 'globcst.inc'     ! Global constants that control model
  INCLUDE 'grid.inc'          ! Grid & map parameters.
  INCLUDE 'bndry.inc'
!
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!-----------------------------------------------------------------------
!
!  NOTE: In order to reduce memory usage, order dependance has been
!        introduced in the calculation of the third term in the
!        U mixing divergence term, which reduces the possible term-
!        wise parallelism.  This is due to the constraint that only two
!        temporary arrays being avaliable for passage into this subroutine.
!
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!  Calculate the vertical term first of the u mixing divergence term.
!  difz(Tau13)
!  Note:  tau13 will be used as a temporary array later...
!
!-----------------------------------------------------------------------

  DO k=2,nz-2
    DO j=2,ny-2
      DO i=2,nx-1
        umix(i,j,k) = dzinv*(tau13(i,j,k+1)-tau13(i,j,k))
      END DO
    END DO
  END DO   ! note tau13 is available for use as a temporary array..

!-----------------------------------------------------------------------
!
!  Calculate the first term of the u mixing divergence term.
!  difx(j3 * Tau11)
!
!-----------------------------------------------------------------------

  IF( tmixvert == 0)THEN   ! compute the horizontal terms.....

    DO k=2,nz-2    !  difx(j3 * Tau11)
      DO j=2,ny-2
        DO i=2,nx-1
          tau13(i,j,k) = dxinv*(j3(i,j,k)*tau11(i,j,k)-                 &
                                j3(i-1,j,k)*tau11(i-1,j,k))
        END DO
      END DO
    END DO

!-----------------------------------------------------------------------
!
!  Calculate the second term of the u mixing divergence term.
!  dify(avgx(avgsv(j3)) * Tau12)
!
!-----------------------------------------------------------------------

    DO k=2,nz-2    !  avgx(aj3y) * Tau12
      DO j=2,ny-1
        DO i=2,nx-1
          tem2(i,j,k) = 0.5*(aj3y(i,j,k)+aj3y(i-1,j,k))*tau12(i,j,k)
        END DO
      END DO
    END DO

    DO k=2,nz-2    !  mapfct*(dify(tem2) + umix)
      DO j=2,ny-2    !  tau13 is the difx term...
        DO i=2,nx-1
          umix(i,j,k) = umix(i,j,k) + mapfct(i,j) *                     &
                   (tau13(i,j,k) + dyinv*(tem2(i,j+1,k)-tem2(i,j,k)))
        END DO
      END DO
    END DO

!-----------------------------------------------------------------------
!
!  Calculate the third part of the third term first.
!  avgy (avgx (j2) * avgz (tau12)). This term is zero when ternopt=0.
!
!-----------------------------------------------------------------------

    IF( ternopt /= 0) THEN

      DO k=3,nz-1    !  avgx (j2) * avgz (tau12)....
        DO j=2,ny-1
          DO i=2,nx-1
            tem1(i,j,k) = 0.25*(j2(i,j,k)+j2(i-1,j,k))*                 &
                          (tau12(i,j,k)+tau12(i,j,k-1))
          END DO
        END DO
      END DO

!
!  At this point, tem1 =  avgx (j2) * avgz (tau12).
!

!-----------------------------------------------------------------------
!
!  Calculate the second part of the third term.
!  j1 * (avgz (avgx (tau11))). This term is zero when ternopt=0.
!
!-----------------------------------------------------------------------

      DO k=3,nz-1    !  j1 * (avgz (avgx (tau11)))....
        DO j=2,ny-2
          DO i=2,nx-1
            tem2(i,j,k)=j1(i,j,k)*0.25*((tau11(i-1,j,k)+tau11(i,j,k))   &
                                     +(tau11(i-1,j,k-1)+tau11(i,j,k-1)))
          END DO
        END DO
      END DO

      DO k=3,nz-1    !  mapfct*(avgy(tem1)+tem2)....
        DO j=2,ny-2
          DO i=2,nx-1
            tem2(i,j,k)=mapfct(i,j)*(tem2(i,j,k)+                       &
                                0.5*(tem1(i,j+1,k)+tem1(i,j,k)))
          END DO
        END DO
      END DO

      DO j=2,ny-2
        DO i=2,nx-1
          tem2(i,j,2)=0.0
        END DO
      END DO

      DO k=2,nz-2    !  add in terrain term.
        DO j=2,ny-2
          DO i=2,nx-1
            umix(i,j,k)=umix(i,j,k)+dzinv*(tem2(i,j,k+1)-tem2(i,j,k))
          END DO
        END DO
      END DO         ! umix is complete.......

    END IF       ! end of terrain if block.........

  END IF  !  end of tmixvert if block....

  RETURN
END SUBROUTINE umixtrm


!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE VMIXTRM                    ######
!######                                                      ######
!######                     Developed by                     ######
!######  Center for the Analysis and Prediction of Storms    ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE vmixtrm(nx,ny,nz,mapfct,j1,j2,j3,aj3x,                       &
           tau12,tau22,tau23,                                           &
           vmix,tem1,tem2)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate the turbulent mixing term in v equation as the
!  divergence of the turbulent momentum fluxes.
!
!  j3 * DIV V = difx (avgy(avgsu(j3)) * tau12) + dify (j3 * tau22)
!               + difz (tau23 + avgx(avgy(j1) * avgz (tau12)) +
!               (j2 * (avgz (avgy(tau22)))))
!
!-----------------------------------------------------------------------
!
!  AUTHOR: D. Weber & M. Xue
!  6/29/92.
!
!  MODIFICATION HISTORY:
!
!  1/23/96 (Donghai Wang & Yuhe Liu)
!  Added the map projection factor.
!
!  9/18/98 (D. Weber)
!  Modified do loop structure to improve code efficiency via:
!  -merging existing loops together
!  -removing and merging operators
!  -switched to DO  ENDDO structure for f90 conversion
!  -added tmixvert option for computing only the vertical
!   turbulent mixing terms.
!  -added array aj3x.
!
!  2/24/1999 (M.Xue)
!  Set contributions of tau21 and tau22 to momentum flux through
!  lower boundary to zero. They are zero anyway in absence of terrain.
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!    mapfct   Map factors at v points
!
!    j1       Coordinate transform Jacobian -d(zp)/dx
!    j2       Coordinate transform Jacobian -d(zp)/dy
!    j3       Coordinate transform Jacobian  d(zp)/dz
!    aj3x     Avgx of the coordinate transformation Jacobian  d(zp)/dz
!
!    tau12    Deformation tensor component (1/s)
!    tau22    Deformation tensor component (1/s)
!    tau23    Deformation tensor component (1/s)
!
!
!  OUTPUT:
!
!    vmix     The divergence of v mixing term
!
!  WORKING ARRAYS:
!
!    tem1     Temporary working array.
!    tem2     Temporary working array.
!
!   (These arrays are defined and used locally (i.e. inside this
!    subroutine), they may also be passed into routines called by
!    this one. Exiting the call to this subroutine, these temporary
!    working arrays may be used for other purposes therefore their
!    contents overwritten. Please examine the usage of working arrays
!    before you alter the code.)
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
!
  INTEGER :: nx,ny,nz          ! Number of grid points in 3 directions

  REAL :: mapfct(nx,ny)        ! Map factors at v points

  REAL :: j1    (nx,ny,nz)     ! Coordinate transformation Jacobian
                               ! - d(zp)/dx
  REAL :: j2    (nx,ny,nz)     ! Coordinate transformation Jacobian
                               ! - d(zp)/dy
  REAL :: j3    (nx,ny,nz)     ! Coordinate transformation Jacobian
                               !   d(zp)/dz
  REAL :: aj3x  (nx,ny,nz)     ! Coordinate transformation Jacobian defined
                               ! as d( zp )/d( z ) AVERAGED IN THE Y-DIR.
!
  REAL :: tau12 (nx,ny,nz)     ! Deformation stress component.
  REAL :: tau22 (nx,ny,nz)     ! Deformation stress component.
  REAL :: tau23 (nx,ny,nz)     ! Deformation stress component.
!
  REAL :: vmix  (nx,ny,nz)     ! Divergence of v mixing term
!
  REAL :: tem1  (nx,ny,nz)     ! Temporary working array.
  REAL :: tem2  (nx,ny,nz)     ! Temporary working array.

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
  INCLUDE 'globcst.inc'     ! Global constants that control model
  INCLUDE 'grid.inc'          ! Grid & map parameters.
  INCLUDE 'bndry.inc'
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
!  Calculate the third term of the V mixing divergence term first.
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!  NOTE: In order to reduce memory usage, order dependance has been
!        introduced in the calculation of the third term in the
!        v mixing divergence term, which reduces the possible term-
!        wise parallelism.  This is due to the constraint that only two
!        temporary arrays being avaliable for passage into this subroutine.
!
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!  Calculate the vetical mixing term of the V mixing divergence term.
!  difz(Tau23)
!
!-----------------------------------------------------------------------

  DO k=2,nz-2    !  difz(Tau23)
    DO j=2,ny-1
      DO i=2,nx-2
        vmix(i,j,k) = dzinv*(tau23(i,j,k+1)-tau23(i,j,k))
      END DO
    END DO         ! if needed tau23 can be used as a temp array.
  END DO         ! vmix for tmixvert=1 is complete.


  IF( tmixvert == 0)THEN   ! compute the horizontal mixing terms.

!-----------------------------------------------------------------------
!
!  Calculate the second term of the V mixing divergence term.
!  dify(j3 * Tau22)
!
!-----------------------------------------------------------------------

    DO k=2,nz-2    !  dify(j3 * Tau22)
      DO j=2,ny-1
        DO i=2,nx-2
          tau23(i,j,k) = dyinv*(j3(i,j,k)*tau22(i,j,k)-                 &
                                j3(i,j-1,k)*tau22(i,j-1,k))
        END DO
      END DO
    END DO

!-----------------------------------------------------------------------
!
!  Calculate the first term of the V mixing divergence term.
!  difx (avgy(avgsu(j3)) * Tau12)
!
!-----------------------------------------------------------------------

    DO k=2,nz-2    !  avgy(aj3x) * Tau12
      DO j=2,ny-1
        DO i=2,nx-1
          tem2(i,j,k) = 0.5*(aj3x(i,j,k)+aj3x(i,j-1,k))*tau12(i,j,k)
        END DO
      END DO
    END DO

    DO k=2,nz-2    !  vmix+ mapfct*(difx(tem2) + tau23)
      DO j=2,ny-1    !                             temp. array
        DO i=2,nx-2
          vmix(i,j,k) = vmix(i,j,k)+mapfct(i,j)*                        &
                      (tau23(i,j,k)+dxinv*(tem2(i+1,j,k)-tem2(i,j,k)))
        END DO
      END DO
    END DO         ! non-terrain portion of vmix is complete.

!-----------------------------------------------------------------------
!
!  Now add in the terrain term......
!  Calculate the second part of the third term (the x-terrain term).
!  avgx(avgy(j1) * avgz(Tau12)). This term is zero when ternopt=0.
!
!-----------------------------------------------------------------------

    IF( ternopt /= 0) THEN

      DO k=3,nz-1    !  avgy (j1) * avgz (tau12)....
        DO j=2,ny-1
          DO i=2,nx-1
            tem1(i,j,k) = 0.25*(j1(i,j,k)+j1(i,j-1,k))*                 &
                            (tau12(i,j,k)+tau12(i,j,k-1))
          END DO
        END DO
      END DO

!-----------------------------------------------------------------------
!
!  Calculate the third part of the third term (the y-terrain term).
!  j2 * avgz(avgy(tau22)). This term is zero when ternopt=0.
!
!-----------------------------------------------------------------------

      DO k=3,nz-1    !  j2 * (avgz (avgy (tau22)))....
        DO j=2,ny-1
          DO i=2,nx-2
            tem2(i,j,k)=j2(i,j,k)*0.25*((tau22(i,j-1,k)+tau22(i,j,k))+  &
                                      (tau22(i,j-1,k-1)+tau22(i,j,k-1)))
          END DO
        END DO
      END DO

      DO k=3,nz-1    !  mapfct*(avgx(tem1)+tem2)....
        DO j=2,ny-1
          DO i=2,nx-2
            tem2(i,j,k)=mapfct(i,j)*(tem2(i,j,k)+                       &
                                0.5*(tem1(i+1,j,k)+tem1(i,j,k)))
          END DO
        END DO
      END DO

      DO j=2,ny-1
        DO i=2,nx-2
          tem2(i,j,2)= 0.0
        END DO
      END DO

      DO k=2,nz-2    !  add in terrain term.
        DO j=2,ny-1
          DO i=2,nx-2
            vmix(i,j,k)=vmix(i,j,k)+dzinv*(tem2(i,j,k+1)-tem2(i,j,k))
          END DO
        END DO
      END DO         ! vmix is complete.......

    END IF           ! end of terrain if block.....

  END IF             ! end of tmixvert if block....

  RETURN
END SUBROUTINE vmixtrm


!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE WMIXTRM                    ######
!######                                                      ######
!######                     Developed by                     ######
!######  Center for the Analysis and Prediction of Storms    ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE wmixtrm(nx,ny,nz,mapfct,j1,j2,j3,aj3x,aj3y,                  &
           tau13,tau23,tau33,                                           &
           wmix,tem1,tem2)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate the turbulent mixing term in w equation as the
!  divergence of the turbulent momentum fluxes.
!
!  j3 * DIV W = difx (avgz(avgsu(j3)) * tau13) + dify(avgz(avgsv(j3))
!               * tau23) + difz (tau33 + avgx(avgz(j1) * avgz(tau13))
!               + (avgy(avgz(j2) * avgz(tau23))))
!
!-----------------------------------------------------------------------
!
!  AUTHOR: D. Weber & M. Xue
!  6/29/92.
!
!  10/17/93  (D. Weber)
!  Vertical loop bounds in loop 300 and 305 corrected.
!
!  1/23/96 (Donghai Wang & Yuhe Liu)
!  Added the map projection factor.
!
!  8/27/98 (D. Weber)
!  Modified code to improve code efficiency via:
!  -merging existing loops together
!  -removing and merging operators
!  -switched to DO  ENDDO structure for f90 conversion
!  -added tmixvert option for computing only the horizontal
!   turbulent mixing terms.
!  -added aj3x,y arrays.
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!    mapfct   Map factors at scalar points
!
!    j1       Coordinate transform Jacobian -d(zp)/dx
!    j2       Coordinate transform Jacobian -d(zp)/dy
!    j3       Coordinate transform Jacobian  d(zp)/dz
!    aj3x     Avgx of the coordinate transformation Jacobian  d(zp)/dz
!    aj3y     Avgy of the coordinate transformation Jacobian  d(zp)/dz
!
!    tau13    Deformation tensor component (1/s)
!    tau23    Deformation tensor component (1/s)
!    tau33    Deformation tensor component (1/s)
!
!
!  OUTPUT:
!
!    wmix     The divergence of w mixing term
!
!  WORKING ARRAYS:
!
!    tem1     Temporary working array.
!    tem2     Temporary working array.
!
!   (These arrays are defined and used locally (i.e. inside this
!    subroutine), they may also be passed into routines called by
!    this one. Exiting the call to this subroutine, these temporary
!    working arrays may be used for other purposes therefore their
!    contents overwritten. Please examine the usage of working arrays
!    before you alter the code.)
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
!
  INTEGER :: nx,ny,nz          ! Number of grid points in 3 directions

  REAL :: mapfct(nx,ny)        ! Map factors at scalar points

  REAL :: j1    (nx,ny,nz)     ! Coordinate transformation Jacobian
                               ! - d(zp)/dx
  REAL :: j2    (nx,ny,nz)     ! Coordinate transformation Jacobian
                               ! - d(zp)/dy
  REAL :: j3    (nx,ny,nz)     ! Coordinate transformation Jacobian
                               !   d(zp)/dz
  REAL :: aj3x  (nx,ny,nz)     ! Coordinate transformation Jacobian defined
                               ! as d( zp )/d( z ) AVERAGED IN THE X-DIR.
  REAL :: aj3y  (nx,ny,nz)     ! Coordinate transformation Jacobian defined
                               ! as d( zp )/d( z ) AVERAGED IN THE Y-DIR.
!
  REAL :: tau13 (nx,ny,nz)     ! Deformation stress component.
  REAL :: tau23 (nx,ny,nz)     ! Deformation stress component.
  REAL :: tau33 (nx,ny,nz)     ! Deformation stress component.
!
  REAL :: wmix  (nx,ny,nz)     ! Divergence of w mixing term
!
  REAL :: tem1  (nx,ny,nz)     ! Temporary working array.
  REAL :: tem2  (nx,ny,nz)     ! Temporary working array.

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
  INCLUDE 'globcst.inc'     ! Global constants that control model
  INCLUDE 'grid.inc'          ! Grid & map parameters.
  INCLUDE 'bndry.inc'       ! Boundry constants
!
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!-----------------------------------------------------------------------
!
!  NOTE: In order to reduce memory usage, order dependance has been
!        introduced in the calculation of the third term in the
!        w mixing divergence term, which reduces the possible term-
!        wise parallelism.  This is due to the constraint that only
!        two temporary arrays being avaliable for passage into this
!        subroutine.
!
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!  Calculate the first term of the third W mixing divergence term.
!  difz(Tau33).
!
!-----------------------------------------------------------------------

  DO k=2,nz-1    ! difz(tau33)
    DO j=2,ny-2    ! note tau33 is used as a tem array after this
      DO i=2,nx-2    ! loop.........
        wmix(i,j,k) = dzinv*(tau33(i,j,k)-tau33(i,j,k-1))
      END DO
    END DO
  END DO         ! vertical mixing for wmix is complete.

!-----------------------------------------------------------------------
!
!  NOTE: tau33 is no longer needed and will be used as a temporary
!        array...
!
!-----------------------------------------------------------------------

  IF( tmixvert == 0)THEN   ! compute the horizontal terms.

!-----------------------------------------------------------------------
!
!  Calculate the first term of the W mixing divergence term.
!  difx (avgz(avgsu(j3)) * Tau13)
!
!-----------------------------------------------------------------------

    DO k=2,nz-1    !  avgz(avgsu(j3)) * Tau13
      DO j=2,ny-2
        DO i=2,nx-1
          tem1(i,j,k) = 0.5*(aj3x(i,j,k)+aj3x(i,j,k-1))*tau13(i,j,k)
        END DO
      END DO
    END DO

!-----------------------------------------------------------------------
!
!  Calculate the second term of the W mixing divergence term.
!  dify(avgz(avgsv(j3)) * Tau23)
!
!-----------------------------------------------------------------------

    DO k=2,nz-1    !  avgz(aj3y) * Tau23
      DO j=2,ny-1
        DO i=2,nx-2
          tem2(i,j,k) = 0.5*(aj3y(i,j,k)+aj3y(i,j,k-1))*tau23(i,j,k)
        END DO
      END DO
    END DO

!-----------------------------------------------------------------------
!
!  Combine the non terrain terms with the vertical term...
!
!-----------------------------------------------------------------------

    DO k=2,nz-1    !  wmix + mapfct*(difx(tem1)+dify(tem2))
      DO j=2,ny-2
        DO i=2,nx-2
          wmix(i,j,k) = wmix(i,j,k)+                                    &
                        mapfct(i,j)*(dxinv*(tem1(i+1,j,k)-tem1(i,j,k))+ &
                                     dyinv*(tem2(i,j+1,k)-tem2(i,j,k)))
        END DO
      END DO
    END DO         ! wmix without terrain is complete.


    IF( ternopt /= 0 ) THEN

!-----------------------------------------------------------------------
!
!  Calculate the second part of the third term.
!  avgx(avgz(j1) * avgz(tau13)). This term is zero when ternopt=0.
!
!-----------------------------------------------------------------------

      DO k=1,nz-1    !  avgz(j1) * avgz(tau13)....
        DO j=2,ny-2
          DO i=2,nx-1
            tem1(i,j,k)=0.25*(j1(i,j,k+1)+   j1(i,j,k))*                &
                          (tau13(i,j,k+1)+tau13(i,j,k))
          END DO
        END DO
      END DO

!-----------------------------------------------------------------------
!
!  Calculate the third part of the third term.
!  avgy(avgz(j2) * avgz(tau23)). This term is zero when ternopt=0.
!
!-----------------------------------------------------------------------

      DO k=1,nz-1    !  avgz(j2) * avgz(tau23)....
        DO j=2,ny-1
          DO i=2,nx-2
            tem2(i,j,k)=0.25*(j2(i,j,k+1)+   j2(i,j,k))*                &
                          (tau23(i,j,k+1)+tau23(i,j,k))
          END DO
        END DO
      END DO


      DO k=1,nz-1    !  mapfct*(avgx(tem1)+avgy(tem2))....
        DO j=2,ny-2
          DO i=2,nx-2
            tau33(i,j,k)=mapfct(i,j)*0.5*((tem1(i+1,j,k)+tem1(i,j,k))+  &
                                          (tem2(i,j+1,k)+tem2(i,j,k)))
          END DO
        END DO
      END DO

      DO k=2,nz-1    !  wmix+difz(tau33)...
        DO j=2,ny-2
          DO i=2,nx-2
            wmix(i,j,k)=wmix(i,j,k)+dzinv*(tau33(i,j,k)-tau33(i,j,k-1))
          END DO
        END DO
      END DO         !  wmix is complete with terrain terms...

    END IF           !  end of terrain if block...

  END IF             !  end of tmixvert if block...

  RETURN
END SUBROUTINE wmixtrm

!
!
!##################################################################
!##################################################################
!######                                                      ######
!######               SUBROUTINE VMIXIMPUVW                  ######
!######                                                      ######
!######                     Developed by                     ######
!######                                                      ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE vmiximpuvw(nx,ny,nz,dtbig1,u,v,w,                            &
           rhostr,kmv,j1,j2,j3inv,                                      &
           uforce,vforce,wforce,                                        &
           rkj3inv2,u_temp,v_temp,w_temp,                               &
           rhostru,rhostrv,rhostrw,a,b,c,d)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Coordinate the vertical mixing terms by using implicit scheme
!  in the u,v,and w momentum equations.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Donghai Wang, M. Xue and X. Song
!  4/19/96
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    nx       Number of grid points in the x-direction
!    ny       Number of grid points in the y-direction
!    nz       Number of grid points in the vertical
!
!    dtbig1   The big time step size
!    u        x component of velocity at a given time level (m/s)
!    v        y component of velocity at a given time level (m/s)
!    w        Vertical component of velocity in Cartesian
!             coordinates at a given time level (m/s)
!    rhostr    Base state density rhobar times j3 (kg/m**3)
!    kmv      Vertical turb. mixing coef. for momentum ( m**2/s )
!    j1       Coordinate transformation Jacobian -d(zp)/dx
!    j2       Coordinate transformation Jacobian -d(zp)/dy
!    j3inv    Inverse of J3
!    uforce   Acoustically inactive forcing terms in u-momentum
!             equation (kg/(m*s)**2). uforce= -uadv + umix + ucorio
!    vforce   Acoustically inactive forcing terms in v-momentum
!             equation (kg/(m*s)**2). vforce= -vadv + vmix + vcorio
!    wforce   Acoustically inactive forcing terms in w-momentum
!             equation (kg/(m*s)**2).
!
!  OUTPUT:
!
!    uforce    Updated uforce.
!    vforce    Updated vforce.
!    wforce    Updated wforce.
!
!  WORK ARRAYS:
!
!    rkj3inv2  Array for calculating a, b, c
!              defined as rhostr*kmv*j3inv*j3inv
!    u_temp    Updated total u-velocity (m/s)
!    v_temp    Updated total v-velocity (m/s)
!    w_temp    Intermediate value of w (m/s)
!    rhostru  Average rhostr at u points (kg/m**3).
!    rhostrv  Average rhostr at v points (kg/m**3).
!    rhostrw  Average rhostr at w points (kg/m**3).
!    a         lhs coefficient of tridigonal eqations
!    b         lhs coefficient of tridigonal eqations
!    c         lhs coefficient of tridigonal eqations
!    d         rhs coefficient of tridigonal eqations
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

  INTEGER :: nx, ny, nz        ! Number of grid points in 3 directions
  REAL :: dtbig1            ! the big time step size
  REAL :: u     (nx,ny,nz)     ! Total u-velocity (m/s)
  REAL :: v     (nx,ny,nz)     ! Total v-velocity (m/s)
  REAL :: w     (nx,ny,nz)     ! Total w-velocity (m/s)
  REAL :: rhostr(nx,ny,nz)     ! Base state density rhobar times j3.
  REAL :: kmv   (nx,ny,nz)     ! Vertical turb. mixing coef. for
                               ! momentum. ( m**2/s )
  REAL :: j1    (nx,ny,nz)     ! Coordinate transformation Jacobian
  REAL :: j2    (nx,ny,nz)     ! Coordinate transformation Jacobian
  REAL :: j3inv (nx,ny,nz)     ! Inverse of J3

  REAL :: uforce(nx,ny,nz)     ! Acoustically inactive forcing terms
                               ! in u-momentum equation (kg/(m*s)**2)
  REAL :: vforce(nx,ny,nz)     ! Acoustically inactive forcing terms
                               ! in v-momentum equation (kg/(m*s)**2)
  REAL :: wforce(nx,ny,nz)     ! Acoustically inactive forcing terms
                               ! in w-momentum equation (kg/(m*s)**2)

  REAL :: rkj3inv2(nx,ny,nz)   ! Array for calculating a, b, c
                               ! defined as rhostr*kmv*j3inv*j3inv
  REAL :: u_temp(nx,ny,nz)     ! Updated total u-velocity (m/s)
  REAL :: v_temp(nx,ny,nz)     ! Updated total v-velocity (m/s)
  REAL :: w_temp(nx,ny,nz)     ! Intermediate value of w
  REAL :: rhostru(nx,ny,nz)    ! Average rhostr at u points (kg/m**3)
  REAL :: rhostrv(nx,ny,nz)    ! Average rhostr at v points (kg/m**3)
  REAL :: rhostrw(nx,ny,nz)    ! Average rhostr at w points (kg/m**3)
  REAL :: a     (nx,ny,nz)     ! lhs coefficient of tridigonal eqations
  REAL :: b     (nx,ny,nz)     ! lhs coefficient of tridigonal eqations
  REAL :: c     (nx,ny,nz)     ! lhs coefficient of tridigonal eqations
  REAL :: d     (nx,ny,nz)     ! rhs coefficient of tridigonal eqations
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'phycst.inc'      ! Physical constants
  INCLUDE 'globcst.inc'     ! Global constants that control model
                            ! execution
  INCLUDE 'grid.inc'          ! Grid & map parameters.
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

  DO k=1,nz-1
    DO j=1,ny-1
      DO i=1,nx-1
        rkj3inv2(i,j,k)=rhostr(i,j,k)*kmv(i,j,k)                        &
                        *j3inv(i,j,k)*j3inv(i,j,k)
      END DO
    END DO
  END DO

  CALL rhouvw(nx,ny,nz,rhostr,rhostru,rhostrv,rhostrw)

  CALL vmiximpu(nx,ny,nz,dtbig1,u,rhostru,rkj3inv2,                     &
                uforce,u_temp,                                          &
                a,b,c,d,w_temp)

  CALL vmiximpv(nx,ny,nz,dtbig1,v,rhostrv,rkj3inv2,                     &
                vforce,v_temp,                                          &
                a,b,c,d,w_temp)

  IF( ternopt == 0 ) THEN

    DO k=2,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          w_temp(i,j,k)=0.0
        END DO
      END DO
    END DO

  ELSE


    DO j=1,ny-1
      DO i=1,nx-1
        w_temp(i,j,nz-1)=0.0
        w_temp(i,j,2) =-((u_temp(i  ,j,2)+u_temp(i  ,j,1))*j1(i,j,2)    &
                        +(u_temp(i+1,j,2)+u_temp(i+1,j,1))*j1(i+1,j,2)  &
                        +(v_temp(i,j  ,2)+v_temp(i,j  ,1))*j2(i,j,2)    &
                        +(v_temp(i,j+1,2)+v_temp(i,j+1,1))*j2(i,j+1,2)  &
                         )*0.25
      END DO
    END DO

  END IF

  CALL vmiximpw(nx,ny,nz,dtbig1,w,rhostrw,w_temp,rkj3inv2,              &
                wforce,                                                 &
                a,b,c,d)


  RETURN
END SUBROUTINE vmiximpuvw
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE VMIXIMPU                   ######
!######                                                      ######
!######                     Developed by                     ######
!######                                                      ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE vmiximpu(nx,ny,nz,dtbig1,u,rhostru,rkj3inv2,                 &
           uforce,u_temp,                                               &
           a,b,c,d,tem2)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate the vertical mixing term using implicit scheme in the u
!  momentum equation.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Donghai Wang, X. Song and M. Xue
!  4/2/96
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    nx        Number of grid points in the x-direction
!    ny        Number of grid points in the y-direction
!    nz        Number of grid points in the vertical
!
!    dtbig1    The big time step size
!    u         Total u-velocity (m/s)
!    rhostru  Average rhostr at u points (kg/m**3).
!    rkj3inv2  Array for calculating a, b, c
!              defined as rhostr*kmv*j3inv*j3inv
!    uforce    Acoustically inactive forcing terms in u-eq.
!              (kg/(m*s)**2)
!
!  OUTPUT:
!
!    uforce    Updated uforce.
!    u_temp    Updated total u-velocity (m/s)
!
!  WORK ARRAYS:
!
!    a         lhs coefficient of tridigonal eqations
!    b         lhs coefficient of tridigonal eqations
!    c         lhs coefficient of tridigonal eqations
!    d         rhs coefficient of tridigonal eqations
!    tem2      Work array
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

  INTEGER :: nx, ny, nz        ! Number of grid points in 3 directions
  REAL :: dtbig1            !the big time step size
  REAL :: u     (nx,ny,nz)     ! Total u-velocity (m/s)
  REAL :: rhostru(nx,ny,nz)    ! Average rhostr at u points (kg/m**3)
  REAL :: rkj3inv2(nx,ny,nz)   ! Array for calculating a, b, c
                               ! defined as rhostr*kmv*j3inv*j3inv
  REAL :: uforce(nx,ny,nz)     ! Acoustically inactive forcing terms
                               ! in u-eq.(kg/(m*s)**2)
  REAL :: u_temp(nx,ny,nz)     ! Updated total u-velocity (m/s)
  REAL :: dt2inv         !Defined as 1/(2*dtbig1),
  REAL :: a     (nx,ny,nz)     ! lhs coefficient of tridigonal eqations
  REAL :: b     (nx,ny,nz)     ! lhs coefficient of tridigonal eqations
  REAL :: c     (nx,ny,nz)     ! lhs coefficient of tridigonal eqations
  REAL :: d     (nx,ny,nz)     ! rhs coefficient of tridigonal eqations
                               ! The contains solution u_temp on exit.

  REAL :: tem2  (nx,ny,nz)     ! Temporary work array

!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'phycst.inc'      ! Physical constants
  INCLUDE 'globcst.inc'     ! Global constants that control model
                            ! execution
  INCLUDE 'grid.inc'          ! Grid & map parameters.
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k
  REAL :: tema,accoef
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  CALL avgsu(rkj3inv2, nx,ny,nz, 1,ny-1, 1,nz,  tem2, a) ! a used as temp
!
!-----------------------------------------------------------------------
!
!  Calculate coefficients of the tridigonal eqation.
!
!-----------------------------------------------------------------------
!

  accoef=0.5*(1.-alfcoef)*dzinv*dzinv
  dt2inv = 0.5/dtbig1

  DO k=2,nz-2
    DO j=1,ny-1
      DO i=1,nx-1
        a(i,j,k)=-accoef*(tem2(i,j,k-1)+tem2(i,j,k  ))
        c(i,j,k)=-accoef*(tem2(i,j,k  )+tem2(i,j,k+1))
        tema    =rhostru(i,j,k)*dt2inv
        b(i,j,k)=-(a(i,j,k)+c(i,j,k))+tema
        d(i,j,k)=uforce(i,j,k)+tema*u(i,j,k)
      END DO
    END DO
  END DO

  DO j=1,ny-1
    DO i=1,nx-1
      b(i,j,2)   =b(i,j,2)   +a(i,j,2)
      b(i,j,nz-2)=b(i,j,nz-2)+c(i,j,nz-2)
    END DO
  END DO
!
!-----------------------------------------------------------------------
!
!  Call the tridiagonal solver.
!
!-----------------------------------------------------------------------
!

  CALL tridiag2(nx,ny,nz,2,nx-1,1,ny-1,2,nz-2,a,b,c,d)

!
!-----------------------------------------------------------------------
!
!  Set d(=u_temp) on the boundaries using zero gradient boundary
!  condition.
!
!-----------------------------------------------------------------------
!

  DO j=1,ny-1
    DO i=2,nx-1
      d(i,j,1   )=d(i,j,2   )
      d(i,j,nz-1)=d(i,j,nz-2)
    END DO
  END DO

!
!-----------------------------------------------------------------------
!
!  Updated uforce and u(=u_temp).
!
!-----------------------------------------------------------------------
!

  DO k=2,nz-2
    DO j=1,ny-1
      DO i=2,nx-1
        uforce(i,j,k)=rhostru(i,j,k)*(d(i,j,k)-u(i,j,k))*dt2inv
      END DO
    END DO
  END DO

  DO k=1,nz-1
    DO j=1,ny-1
      DO i=2,nx-1
        u_temp(i,j,k)=d(i,j,k)
      END DO
    END DO
  END DO

  RETURN
END SUBROUTINE vmiximpu

!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE VMIXIMPV                   ######
!######                                                      ######
!######                     Developed by                     ######
!######                                                      ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE vmiximpv(nx,ny,nz,dtbig1,v,rhostrv,rkj3inv2,                 &
           vforce,v_temp,                                               &
           a,b,c,d,tem2)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate the vertical mixing term using implicit scheme in the v
!  momentum equation.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Donghai Wang, X. Song and M. Xue
!  4/2/96
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    nx        Number of grid points in the x-direction
!    ny        Number of grid points in the y-direction
!    nz        Number of grid points in the vertical
!
!    dtbig1    the big time step size
!    v         Total v-velocity (m/s)
!    rhostrv  Average rhostr at v points (kg/m**3).
!    rkj3inv2  Array for calculating a, b, c
!              defined as rhostr*kmv*j3inv*j3inv
!    vforce    Acoustically inactive forcing terms in v-eq.
!              (kg/(m*s)**2)
!
!  OUTPUT:
!
!    vforce    Updated vforce.
!    v_temp    Updated total v-velocity (m/s)
!
!  WORK ARRAYS:
!
!    a         lhs coefficient of tridigonal eqations
!    b         lhs coefficient of tridigonal eqations
!    c         lhs coefficient of tridigonal eqations
!    d         rhs coefficient of tridigonal eqations
!    rhostrv   Work array
!    tem2      Work array
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

  INTEGER :: nx, ny, nz        ! Number of grid points in 3 directions
  REAL :: dtbig1            !the big time step size
  REAL :: v     (nx,ny,nz)     ! Total v-velocity (m/s)
  REAL :: rhostrv(nx,ny,nz)    ! Average rhostr at v points (kg/m**3)
  REAL :: rkj3inv2(nx,ny,nz)   ! Array for calculating a, b, c
                               ! defined as rhostr*kmv*j3inv*j3inv

  REAL :: vforce(nx,ny,nz)     ! Acoustically inactive forcing terms
                               ! in v-eq.(kg/(m*s)**2)
  REAL :: v_temp(nx,ny,nz)     ! Updated total v-velocity (m/s)

  REAL :: dt2inv            !Defined as 1/(2*dtbig1),
  REAL :: a     (nx,ny,nz)     ! lhs coefficient of tridigonal eqations
  REAL :: b     (nx,ny,nz)     ! lhs coefficient of tridigonal eqations
  REAL :: c     (nx,ny,nz)     ! lhs coefficient of tridigonal eqations
  REAL :: d     (nx,ny,nz)     ! rhs coefficient of tridigonal eqations
                               ! The contains solution v_temp on exit.

  REAL :: tem2  (nx,ny,nz)     ! Temporary work array
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'phycst.inc'      ! Physical constants
  INCLUDE 'globcst.inc'     ! Global constants that control model
                            ! execution
  INCLUDE 'grid.inc'          ! Grid & map parameters.
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k
  REAL :: tema,accoef
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  CALL avgsv(rkj3inv2, nx,ny,nz, 1,nx-1, 1,nz, tem2, a) ! a used as temp
!
!-----------------------------------------------------------------------
!
!  Calculate coefficients of the tridigonal eqation.
!
!-----------------------------------------------------------------------
!

  accoef=0.5*(1.-alfcoef)*dzinv*dzinv
  dt2inv = 0.5/dtbig1

  DO k=2,nz-2
    DO j=1,ny-1
      DO i=1,nx-1
        a(i,j,k)=-accoef*(tem2(i,j,k-1)+tem2(i,j,k  ))
        c(i,j,k)=-accoef*(tem2(i,j,k  )+tem2(i,j,k+1))
        tema    =rhostrv(i,j,k)*dt2inv
        b(i,j,k)=-(a(i,j,k)+c(i,j,k))+tema
        d(i,j,k)=vforce(i,j,k)+tema*v(i,j,k)
      END DO
    END DO
  END DO

  DO j=1,ny-1
    DO i=1,nx-1
      b(i,j,2)   =b(i,j,2)   +a(i,j,2)
      b(i,j,nz-2)=b(i,j,nz-2)+c(i,j,nz-2)
    END DO
  END DO

!
!-----------------------------------------------------------------------
!
!  Call the tridiagonal solver.
!
!-----------------------------------------------------------------------
!

  CALL tridiag2(nx,ny,nz,1,nx-1,2,ny-1,2,nz-2,a,b,c,d)

!
!-----------------------------------------------------------------------
!
!  Set d(=v_temp) on the boundaries using zero gradient boundary
!  condition.
!
!-----------------------------------------------------------------------
!

  DO j=2,ny-1
    DO i=1,nx-1
      d(i,j,1   )=d(i,j,2   )
      d(i,j,nz-1)=d(i,j,nz-2)
    END DO
  END DO
!
!-----------------------------------------------------------------------
!
!  Updated vforce and v(=v_temp).
!
!-----------------------------------------------------------------------
!

  DO k=2,nz-2
    DO j=2,ny-1
      DO i=1,nx-1
        vforce(i,j,k)=rhostrv(i,j,k)*(d(i,j,k)-v(i,j,k))*dt2inv
      END DO
    END DO
  END DO

  DO k=1,nz-1
    DO j=2,ny-1
      DO i=1,nx-1
        v_temp(i,j,k)=d(i,j,k)
      END DO
    END DO
  END DO

  RETURN
END SUBROUTINE vmiximpv


!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE VMIXIMPW                   ######
!######                                                      ######
!######                     Developed by                     ######
!######                                                      ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE vmiximpw(nx,ny,nz,dtbig1,w,rhostrw,w_temp,rkj3inv2,          &
           wforce,                                                      &
           a,b,c,d)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate the vertical mixing term using implicit scheme in the w
!  momentum equation.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Donghai Wang, X. Song and M. Xue
!  4/2/96
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    nx       Number of grid points in the x-direction
!    ny       Number of grid points in the y-direction
!    nz       Number of grid points in the vertical
!
!    dtbig1   The big time step size
!    w        Total w-velocity (m/s)
!    rhostrw  Average rhostr at w points (kg/m**3).
!    w_temp   Intermediate value of w
!    rkj3inv2 Array for calculating a, b, c
!             defined as rhostr*kmv*j3inv*j3inv
!    wforce   Acoustically inactive forcing terms in w-eq.
!             (kg/(m*s)**2)
!
!  OUTPUT:
!
!    wforce   Updated wforce.
!
!  WORK ARRAYS:
!
!    dt2inv   Defined as 1/(2*dtbig1)
!    a        lhs coefficient of tridigonal eqations
!    b        lhs coefficient of tridigonal eqations
!    c        lhs coefficient of tridigonal eqations
!    d        rhs coefficient of tridigonal eqations
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

  INTEGER :: nx, ny, nz        ! Number of grid points in 3 directions
  REAL :: dtbig1            ! the big time step size
  REAL :: w     (nx,ny,nz)     ! Total w-velocity (m/s)
  REAL :: rhostrw(nx,ny,nz)    ! Average rhostr at w points (kg/m**3)
  REAL :: w_temp(nx,ny,nz)     ! Intermediate value of w
  REAL :: rkj3inv2(nx,ny,nz)   ! Array for calculating a, b, c
                               ! defined as rhostr*kmv*j3inv*j3inv

  REAL :: wforce(nx,ny,nz)     ! Acoustically inactive forcing terms
                               ! in w-eq.(kg/(m*s)**2)

  REAL :: dt2inv            ! Defined as 1/(2*dtbig1),
  REAL :: a     (nx,ny,nz)     ! lhs coefficient of tridigonal eqations
  REAL :: b     (nx,ny,nz)     ! lhs coefficient of tridigonal eqations
  REAL :: c     (nx,ny,nz)     ! lhs coefficient of tridigonal eqations
  REAL :: d     (nx,ny,nz)     ! rhs coefficient of tridigonal eqations
                               ! The contains solution updated-w on exit

!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'phycst.inc'      ! Physical constants
  INCLUDE 'globcst.inc'     ! Global constants that control model
                            ! execution
  INCLUDE 'grid.inc'          ! Grid & map parameters.
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k
  REAL :: tema,accoef
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
!  Calculate coefficients of the tridigonal eqation.
!
!-----------------------------------------------------------------------
!

  accoef=2.0*(1.-alfcoef)*dzinv*dzinv
  dt2inv = 0.5/dtbig1

  DO k=2,nz-1
    DO j=1,ny-1
      DO i=1,nx-1
        a(i,j,k)=-accoef*rkj3inv2(i,j,k  )
        c(i,j,k)=-accoef*rkj3inv2(i,j,k+1)
        tema    =rhostrw(i,j,k)*dt2inv
        b(i,j,k)=-(a(i,j,k)+c(i,j,k))+tema
        d(i,j,k)=wforce(i,j,k)+tema*w(i,j,k)
      END DO
    END DO
  END DO

  DO j=1,ny-1
    DO i=1,nx-1
      d(i,j,3)   =d(i,j,3)   - w_temp(i,j,2)*a(i,j,3)
      d(i,j,nz-2)=d(i,j,nz-2)
    END DO
  END DO

!
!-----------------------------------------------------------------------
!
!  Call the tridiagonal solver.
!
!-----------------------------------------------------------------------
!

  CALL tridiag2(nx,ny,nz,1,nx-1,1,ny-1,3,nz-2,a,b,c,d)

!
!-----------------------------------------------------------------------
!
!  Set d on the boundaries using zero gradient boundary
!  condition.
!
!-----------------------------------------------------------------------
!

  DO j=1,ny-1
    DO i=1,nx-1
      d(i,j,2 )= w_temp(i,j,2)
      d(i,j,nz-1)=0.0
    END DO
  END DO

!
!-----------------------------------------------------------------------
!
!  Updated wforce.
!
!-----------------------------------------------------------------------
!
  DO k=2,nz-1
    DO j=1,ny-1
      DO i=1,nx-1
        wforce(i,j,k)=rhostrw(i,j,k)*(d(i,j,k)-w(i,j,k))*dt2inv
      END DO
    END DO
  END DO

  RETURN
END SUBROUTINE vmiximpw

!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE VMIXIMPS                   ######
!######                                                      ######
!######                     Developed by                     ######
!######                                                      ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE vmiximps(nx,ny,nz,dtbig1,s,rhostr,rkj3inv2,                  &
           sforce,a,b,c,d)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate the vertical mixing term using implicit scheme in the
!  scalar equations, including potential temperature, presure, water
!  substances and TKE.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Donghai Wang, X. Song and M. Xue
!  4/2/96
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    nx  Number of grid points in the x-direction
!    ny  Number of grid points in the y-direction
!    nz  Number of grid points in the vertical
!
!    dtbig1    the big time step size
!    s        A scalar varable
!    rhostr   Base state density rhobar times j3 (kg/m**3)
!    rkj3inv2 Array for calculating a, b, c
!             defined as rhostr*kmv*j3inv*j3inv
!    sforce   Acoustically inactive forcing terms in a scalar-eq.
!
!  OUTPUT:
!
!    sforce   Updated sforce.
!
!  WORK ARRAYS:
!
!    a        lhs coefficient of tridigonal eqations
!    b        lhs coefficient of tridigonal eqations
!    c        lhs coefficient of tridigonal eqations
!    d        rhs coefficient of tridigonal eqations
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

  INTEGER :: nx, ny, nz        ! Number of grid points in 3 directions
  REAL :: dtbig1            ! the big time step size
  REAL :: s     (nx,ny,nz)     ! A total scalar varable
  REAL :: rhostr(nx,ny,nz)     ! Base state density rhobar times j3.
  REAL :: rkj3inv2(nx,ny,nz)   ! Array for calculating a, b, c
                               ! defined as rhostr*kmv*j3inv*j3inv

  REAL :: sforce(nx,ny,nz)     ! Acoustically inactive forcing terms
                               ! in a scalar-eq

  REAL :: dt2inv            ! Defined as 1/(2*dtbig1),
  REAL :: a     (nx,ny,nz)     ! lhs coefficient of tridigonal eqations
  REAL :: b     (nx,ny,nz)     ! lhs coefficient of tridigonal eqations
  REAL :: c     (nx,ny,nz)     ! lhs coefficient of tridigonal eqations
  REAL :: d     (nx,ny,nz)     ! rhs coefficient of tridigonal eqations
                               ! The contains solution ustar on exit.

!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'phycst.inc'      ! Physical constants
  INCLUDE 'globcst.inc'     ! Global constants that control model
                            ! execution
  INCLUDE 'grid.inc'          ! Grid & map parameters.
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k
  REAL :: tema,accoef
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
!  Calculate coefficients of the tridigonal eqation.
!
!-----------------------------------------------------------------------
!

  accoef=0.5*(1.-alfcoef)*dzinv*dzinv
  dt2inv = 0.5/dtbig1

  DO k=2,nz-2
    DO j=1,ny-1
      DO i=1,nx-1
        a(i,j,k)=-accoef*(rkj3inv2(i,j,k-1)+rkj3inv2(i,j,k  ))
        c(i,j,k)=-accoef*(rkj3inv2(i,j,k  )+rkj3inv2(i,j,k+1))
        tema    =rhostr(i,j,k)*dt2inv
        b(i,j,k)=-(a(i,j,k)+c(i,j,k))+tema
        d(i,j,k)=sforce(i,j,k)+tema*s(i,j,k)
      END DO
    END DO
  END DO

  DO j=1,ny-1
    DO i=1,nx-1
      b(i,j,2)   =b(i,j,2)   +a(i,j,2)
      b(i,j,nz-2)=b(i,j,nz-2)+c(i,j,nz-2)
    END DO
  END DO

!
!-----------------------------------------------------------------------
!
!  Call the tridiagonal solver.
!
!-----------------------------------------------------------------------
!

  CALL tridiag2(nx,ny,nz,1,nx-1,1,ny-1,2,nz-2,a,b,c,d)

!
!-----------------------------------------------------------------------
!
!  Set d on the boundaries using zero gradient boundary
!  condition.
!
!-----------------------------------------------------------------------
!

  DO j=1,ny-1
    DO i=1,nx-1
      d(i,j,1   )=d(i,j,2   )
      d(i,j,nz-1)=d(i,j,nz-2)
    END DO
  END DO

!
!-----------------------------------------------------------------------
!
!  Updated sforce.
!
!-----------------------------------------------------------------------
!

  DO k=1,nz-1
    DO j=1,ny-1
      DO i=1,nx-1
        sforce(i,j,k)=rhostr(i,j,k)*(d(i,j,k)-s(i,j,k))*dt2inv
      END DO
    END DO
  END DO

  RETURN
END SUBROUTINE vmiximps
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE TRIDIAG2                   ######
!######                                                      ######
!######                     Developed by                     ######
!######                                                      ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE tridiag2(nx,ny,nz,ibgn,iend,jbgn,jend,                       &
           kbgn,kend,a,b,c,d)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate the solution for a tridiagonal equations.
!  Note:
!       a: left of main diagonal, useful range [kbgn+1,kend];
!       b: main diagonal, useful range [kbgn,kend];
!       c: right of main diagonal, useful range [kbgn,kend-1];
!       d: right hand side of equations.
!
!-----------------------------------------------------------------------
!
!
!  AUTHOR: X. Song, Donghai Wang and M. Xue
!  4/2/96
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    nx       Number of grid points in the x-direction
!    ny       Number of grid points in the y-direction
!    nz       Number of grid points in the vertical
!
!    ibgn     i-index where operation begins.
!    iend     i-index where operation ends.
!    jbgn     j-index where operation begins.
!    jend     j-index where operation ends.
!    kbgn     k-index where operation begins.
!    kend     k-index where operation ends.
!
!    a        left of main diagonal, useful range [kbgn+1,kend]
!    b        main diagonal, useful range [kbgn,kend]
!    c        right of main diagonal, useful range [kbgn,kend-1]
!    d        right hand side of equations
!
!  OUTPUT:
!
!    d        The solution
!
!-----------------------------------------------------------------------
!
!  Variable Declarations:
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE             ! Force explicit declarations

  INTEGER :: nx, ny, nz        ! Number of grid points in 3 directions

  REAL :: c(nx,ny,nz)       ! Right of main diagonal
  REAL :: a(nx,ny,nz)       ! Left of main diagonal
  REAL :: d(nx,ny,nz)       ! Right hand side of equations
  REAL :: b(nx,ny,nz)       ! Main diagonal

  INTEGER :: ibgn,iend,jbgn,jend,kbgn,kend
  INTEGER :: i,j,k
  REAL :: r
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  DO k=kbgn+1,kend
    DO j=jbgn,jend
      DO i=ibgn,iend
        r=a(i,j,k)/b(i,j,k-1)
        b(i,j,k)=b(i,j,k)-r*c(i,j,k-1)
        d(i,j,k)=d(i,j,k)-r*d(i,j,k-1)
      END DO
    END DO
  END DO

  DO j=jbgn,jend
    DO i=ibgn,iend
      d(i,j,kend)=d(i,j,kend)/b(i,j,kend)
    END DO
  END DO

  DO k=kend-1,kbgn,-1
    DO j=jbgn,jend
      DO i=ibgn,iend
        d(i,j,k)=(d(i,j,k)-c(i,j,k)*d(i,j,k+1))/b(i,j,k)
      END DO
    END DO
  END DO

  RETURN
END SUBROUTINE tridiag2


