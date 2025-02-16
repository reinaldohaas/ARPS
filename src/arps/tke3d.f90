!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE SOLVTKE                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE solvtke(RK3step,nx,ny,nz,dtbig1,u,v,wcont,                   &
           ptprt,pprt,qv,qscalar,tke,                                   &
           ubar,vbar,ptbar,pbar,ptbari,rhostr,rhostri,qvbar,            &
           x,y,z,zp, mapfct, j1,j2,j3,aj3x,aj3y,j3inv,                  &
           kmh,kmv,rprntl,lenscl,defsq,                                 &
           ptsflx,qvsflx,                                               &
           tem1,tem2,tem3,tem4,tem5,tem6,tem7,tem8,                     &
           tem9,tem10,tem11,tkeadv,tkeforce,tem1_0,tem2_0,tem3_0)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Compute the forcing terms in the turbulence kinetic energy (TKE)
!  equation and integrate the TKE equation one time step forward.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: V.Wong, Y.Tang and X. Song
!  10/1993
!
!  MODIFICATION HISTORY:
!
!  04/1994
!  V.Wong, M.Xue and X. Song
!
!  9/1/94 (Y. Lu)
!  Cleaned up documentation.
!
!  1/23/96 (Donghai Wang & Yuhe Liu)
!  Added the map projection factor.
!
!  2/26/96 (M.Xue, X. Song and V. Wong)
!  Reganized the call to SOLVTKE. TKE equation is now solved at the
!  same time level as the other scalar equations.
!
!  3/8/96 (M. Xue, X. Song and V. Wong)
!  Add parameter tkeopt for three versions of 1.5 order TKE schemes.
!  The differ mainly in the specification of turbulence mixing length.
!
!  3/11/96 (M. Xue)
!  Corrected time level error for tke mixing and dissipation
!  terms. They should be evaluated at tpast for loapfrog step.
!  Also all three time levels of scalars are now passed into
!  this subroutine.
!
!  4/1/96 (Donghai Wang, X. Song and M. Xue)
!  Added the implicit treatment of the vertical diffusion term
!  in TKE equation.
!
!  7/10/1997 (Fanyou Kong - CMRP)
!  Fixed a bug in FCT advection mode with the tem1_0(),tem2_0(),
!  and tem3_0() corrected
!  Added MPDCD positive advection option (sadvopt = 5)
!
!  7/17/1998 (M. Xue)
!  Changed call to ADVQFCT.
!
!  9/18/1998 (D. Weber)
!  Added aj3x,y rhostri arrays.
!
!  1/18/1999 (W. Martin and M. Xue)
!  Changed the coefficient in the dissipation term from 3.9
!  to 0.41 for the lowest two levels from 3.9 for Sun and Chang
!  scheme (tkeopt=3).
!
!  10/19/2009 (D. Dawson)
!  Updated to handle the RK3 time integration option, both for the case
!  in which only the advection terms are updated each RK3 step, and for
!  the case in which all forcing terms are updated each step.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    RK3step  Which RK3 step are we on (1,2,3)?
!             (only used for tintegopt == 2)
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!    dtbig1   Time step. If frstep=1 ,dtbig1=dtbig/2, otherwise,
!             dtbig1=dtbig.
!
!    u        x component of velocity at time tpast (m/s)
!    v        y component of velocity at time tpast (m/s)
!    wcont    Vertical component of Cartesian velocity at time
!             tpast (m/s)
!             computational coordinates (m/s)
!    ptprt    Perturbation potential temperature at times tpast and
!             tpresent (K)
!    pprt     Perturbation pressure at times tpast and tpresent
!             (Pascal)
!    qv       Water vapor specific humidity at times tpast and
!             tpresent (kg/kg)
!
!    qscalar  Cloud water mixing ratio at times tpast and tpresent
!             (kg/kg)
!
!    tke      Turbulent kinetic energy at times tpast and tpresent.
!
!    ubar     Base state zonal velocity component (m/s)
!    vbar     Base state meridional velocity component (m/s)
!    ptbar    Base state potential temperature (K)
!    pbar     Base state pressure (Pascal)
!    ptbari   Inverse base state potential temperature (K)
!    rhostr   Base state density rhobar times j3 (kg/m**3)
!    rhostri  Inverse base state density rhobar times j3 (kg/m**3)
!    qvbar    Base state water vapor specific humidity (kg/kg)
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
!    kmh      Horizontal turb. mixing coef. for momentum ( m**2/s )
!    kmv      Vertical turb. mixing coef. for momentum ( m**2/s )
!    rprntl   Reciprocal of Prandtl number
!    lenscl   Turbulent mixing length scale (m)
!    defsq    Deformation squared (1/s**2)
!    ptsflx   Surface heat flux
!    qvsflx   Surface moisture flux
!
!  OUTPUT:
!
!    tke      Updated (by dtbig) Turbulent Kinetic Energy at time
!             tfuture.
!
!  WORK ARRAYS:
!
!    tem1     Temporary work array
!    tem2     Temporary work array
!    tem3     Temporary work array
!    tem4     Temporary work array
!    tem5     Temporary work array
!    tem6     Temporary work array
!    tem7     Temporary work array
!    tem8     Temporary work array
!    tem9     Temporary work array
!    tkeforce Temporary work array
!    tketot   Temporary work array
!    tem1_0   Temporary work array
!    tem2_0   Temporary work array
!    tem3_0   Temporary work array
!
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE             ! Force explicit declarations
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'phycst.inc'
  INCLUDE 'globcst.inc'
  INCLUDE 'bndry.inc'       ! Boundary condition control parameters
  INCLUDE 'mp.inc'            ! Message passing parameters.
  INCLUDE 'timelvls.inc'
!
!-----------------------------------------------------------------------
!
!  Variable Declarations:
!
!-----------------------------------------------------------------------
!
  REAL :: dtbig1               ! Local value of big timestep

  INTEGER :: RK3step
  INTEGER :: nx, ny, nz        ! Number of grid points in 3 directions

  REAL :: u     (nx,ny,nz,nt)  ! Total u-velocity at time tpast (m/s)
  REAL :: v     (nx,ny,nz,nt)  ! Total v-velocity at time tpast (m/s)
  REAL :: wcont (nx,ny,nz)     ! Vertical component of Cartesian
                               ! velocity at
                               ! time tpast (m/s)

  REAL :: ptprt (nx,ny,nz,nt)  ! Perturbation potential temperature (K)
  REAL :: pprt  (nx,ny,nz,nt)  ! Perturbation pressure (Pascal)

  REAL :: qv    (nx,ny,nz,nt)  ! Water vapor specific humidity (kg/kg)

  REAL :: qscalar(nx,ny,nz,nt,nscalar)

  REAL :: tke   (nx,ny,nz,nt)  ! Turbulent Kinetic Energy ((m/s)**2)

  REAL :: ubar  (nx,ny,nz)     ! Base state u-velocity (m/s)
  REAL :: vbar  (nx,ny,nz)     ! Base state v-velocity (m/s)
  REAL :: ptbar (nx,ny,nz)     ! Base state potential temperature (K)
  REAL :: pbar  (nx,ny,nz)     ! Base state pressure (Pascal).
  REAL :: ptbari(nx,ny,nz)     ! Inverse base state pot. temperature (K)
  REAL :: rhostr(nx,ny,nz)     ! Base state density rhobar times
                               ! j3 (kg/m**3)
  REAL :: rhostri(nx,ny,nz)    ! Inverse base state density rhobar times
                               ! j3 (kg/m**3)
  REAL :: qvbar (nx,ny,nz)     ! Base state water vapor specific
                               ! humidity
                               ! (kg/kg)

  REAL :: x     (nx)           ! x-coord. of the physical and
                               ! computational grid.
                               ! Defined at u-point.
  REAL :: y     (ny)           ! y-coord. of the physical and
                               ! computational grid.
                               ! Defined at v-point.
  REAL :: z     (nz)           ! z-coord. of the
                               ! computational grid.
                               ! Defined at w-point on the
                               ! staggered grid.
  REAL :: zp    (nx,ny,nz)     ! Physical height coordinate
                               ! defined at
                               ! w-point of the staggered grid.

  REAL :: mapfct(nx,ny,8)      ! Map factors at scalar points

  REAL :: j1    (nx,ny,nz)     ! Coordinate transformation
                               ! Jacobian defined as
                               ! - d( zp )/d( x ).
  REAL :: j2    (nx,ny,nz)     ! Coordinate transformation
                               ! Jacobian defined as
                               ! - d( zp )/d( y ).
  REAL :: j3    (nx,ny,nz)     ! Coordinate transformation
                               ! Jacobian defined as
                               ! d( zp )/d( z ).
  REAL :: aj3x  (nx,ny,nz)     ! Coordinate transformation Jacobian defined
                               ! as d( zp )/d( z ) AVERAGED IN THE X-DIR.
  REAL :: aj3y  (nx,ny,nz)     ! Coordinate transformation Jacobian defined
                               ! as d( zp )/d( z ) AVERAGED IN THE Y-DIR.
  REAL :: j3inv (nx,ny,nz)     ! Inverse of j3

  REAL :: kmh   (nx,ny,nz)     ! Horizontal turb. mixing coef. for
                               ! momentum. ( m**2/s )
  REAL :: kmv   (nx,ny,nz)     ! Vertical turb. mixing coef. for
                               ! momentum. ( m**2/s )
  REAL :: rprntl(nx,ny,nz)     ! Reciprocal of Prandtl number
  REAL :: lenscl(nx,ny,nz)     ! Turbulent mixing length scale (m)
  REAL :: defsq (nx,ny,nz)     ! Deformation squared (1/s**2)

  REAL :: ptsflx(nx,ny)        ! Surface heat flux (K*kg/(m*s**2))
  REAL :: qvsflx(nx,ny)        ! Surface moisture flux (kg/(m**2*s))

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
  REAL :: tkeadv(nx,ny,nz)     ! Work array for the advection forcing terms.
  REAL :: tkeforce(nx,ny,nz)   ! Work array for all th other the forcing terms.

  REAL :: tem1_0(0:nx,0:ny,0:nz)     ! Temporary work array
  REAL :: tem2_0(0:nx,0:ny,0:nz)     ! Temporary work array
  REAL :: tem3_0(0:nx,0:ny,0:nz)     ! Temporary work array

!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k, it
  REAL    :: eps
  REAL    :: dt2
  INTEGER :: tstrtlvl
  REAL    :: deltat
  LOGICAL :: RK3sadvflag      ! Flag to indicate if sadvopt == 4 and RK3
                              ! time-stepping is being used.  It's used
                              ! to allow the normal advection to be called
                              ! for the first 2 RK3 steps, and the full
                              ! FCT only for the last one.

!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  eps=1.e-20
!
!-----------------------------------------------------------------------
!
!  Compute the advection term for TKE and store the result in tem1
!
!-----------------------------------------------------------------------
!
  CALL set_acct(advs_acct)

  IF(sadvopt == 4 .or. tintegopt == 2 .or. tintegopt == 3) THEN
    ! Forward-based FCT scheme or RK3 time-stepping
    deltat = dtbig1
    tstrtlvl = tpresent
  ELSE
    deltat = dtbig1*2
    tstrtlvl = tpast
  END IF

  ! If sadvopt == 4, but we are on one of the first 2 RK3 steps,
  ! we need to temporarily set sadvopt to 3
  IF(sadvopt == 4 .and. (RK3step == 1 .or. RK3step == 2)) THEN
    sadvopt = 3
    RK3sadvflag = .true.
  ELSE
    RK3sadvflag = .false.
  END IF

  IF (sadvopt == 1 .OR. sadvopt == 2 .OR. sadvopt == 3) THEN
                                          ! 2nd or 4th order advection
    CALL rhouvw(nx,ny,nz,rhostr,tem1,tem2,tem3)

    DO k=1,nz-1
      DO j=1,ny-1
        DO i=1,nx
          tem4(i,j,k)=u(i,j,k,2)*tem1(i,j,k)
        END DO
      END DO
    END DO

    DO k=1,nz-1
      DO j=1,ny
        DO i=1,nx-1
          tem5(i,j,k)=v(i,j,k,2)*tem2(i,j,k)
        END DO
      END DO
    END DO

    DO k=1,nz
      DO j=1,ny-1
        DO i=1,nx-1
          tem6(i,j,k)=wcont(i,j,k)*tem3(i,j,k)
        END DO
      END DO
    END DO

    DO k=1,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          tem2(i,j,k)= 0.0
        END DO
      END DO
    END DO

    CALL advq(nx,ny,nz,tke,u,v,wcont,tem4,tem5,tem6,                    &
              rhostr,tem2,mapfct,                                       &
              tkeadv,                                                   &
              tem3,tem7,tem8,tem9)

  ELSE IF( sadvopt == 4 .OR. sadvopt == 5) THEN   ! FCT advection

    DO k=1,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          tem3_0(i,j,k)=rhostr(i,j,k)
        END DO
      END DO
    END DO

    CALL extndsbc(tem3_0,nx,ny,nz,0,ebc,wbc,nbc,sbc,tbc,bbc)

    DO k=1,nz-1
      DO j=1,ny-1
        DO i=1,nx
          tem9(i,j,k)=u(i,j,k,2)*(tem3_0(i-1,j,k)+tem3_0(i,j,k))        &
                    *mapfct(i,j,5)*0.5
        END DO
      END DO
    END DO

    DO k=1,nz-1
      DO j=1,ny
        DO i=1,nx-1
          tem10(i,j,k)=v(i,j,k,2)*(tem3_0(i,j-1,k)+tem3_0(i,j,k))       &
                    *mapfct(i,j,6)*0.5
        END DO
      END DO
    END DO

    DO k=1,nz
      DO j=1,ny-1
        DO i=1,nx-1
          tem11(i,j,k)=wcont(i,j,k)                                     &
                    *(tem3_0(i,j,k-1)+tem3_0(i,j,k))*0.5
        END DO
      END DO
    END DO

    CALL advqfct(nx,ny,nz,dtbig1,tke,u,v,wcont,tem9,tem10,tem11,        &
                 rhostr,rhostri,mapfct,j3,tkeadv,                       &
                 tem1,tem2,tem3,tem4,tem5,tem6,tem7,tem8,               &
                 tem1_0,tem2_0,tem3_0)

  END IF

  IF(RK3sadvflag) sadvopt = 4   ! Reset sadvopt back to 4 if needed

  DO k=1,nz-1
    DO j=1,ny-1
      DO i=1,nx-1
        tkeadv(i,j,k)=-tkeadv(i,j,k)
      END DO
    END DO
  END DO

  ! If this is the first RK3 step, or if using the old leapfrog method
  ! Then calculate all the additional forcing terms.
  ! In the case of RK3, there are two options.  In the first (tintegopt == 2)
  ! Only the advection terms (calculated above) are recalcuated each RK3
  ! step, and added to the final tkeforce for the time integration.
  ! In the second option (tintegopt == 3), all the forcing terms are
  ! recalculated each RK3 step (more expensive).  The following block
  ! of code handles each of these possbilities.

  IF( tintegopt == 1 .OR. tintegopt == 3 .OR.                           &
      (tintegopt == 2 .and. RK3step == 1) ) THEN
!
!-----------------------------------------------------------------------
!
!  Calculate the TKE diffusion term and store the result in the array
!  tkeforce.
!
!  2*[d( rhobar*km*d(tke) /dx )/dx +d( rhobar*km*d(tke) /dy )/dy +
!    d( rhobar*km*d(tke) /dz )/dz) ]
!  + computational mixing
!
!-----------------------------------------------------------------------
!
  CALL mixtke(nx,ny,nz,tke(1,1,1,tstrtlvl),rhostr,kmh,kmv,              &
              x,y,z,zp,mapfct,j1,j2,j3,aj3x,aj3y,j3inv,tem8,            &
              tem1,tem2,tem3,tem4,tem5,tem6,tem7)

  CALL set_acct(tkesrc_acct)
  DO k=1,nz-1
    DO j=1,ny-1
      DO i=1,nx-1
        tkeforce(i,j,k)=tem8(i,j,k)
      END DO
    END DO
  END DO
!
!-----------------------------------------------------------------------
!
!  Calculate the TKE dissipation term and add the result to array
!  tkeforce.
!
!-----------------------------------------------------------------------
!
  CALL stgrdscl(nx,ny,nz, zp, tem1)
!
  IF (tkeopt == 1) THEN      ! Wyngaard formulation

    DO k=3,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          tem6(i,j,k)=0.93/tem1(i,j,k)
        END DO
      END DO
    END DO

    DO k=1,2
      DO j=1,ny-1
        DO i=1,nx-1
          tem6(i,j,k)=3.9/tem1(i,j,k)
        END DO
      END DO
    END DO

  ELSE IF (tkeopt == 2) THEN  ! Deardroff formulation

    DO k=3,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          tem6(i,j,k)=(0.19+0.51*lenscl(i,j,k)/tem1(i,j,k))             &
                      /MAX(lenscl(i,j,k),0.1*tem1(i,j,k))
        END DO
      END DO
    END DO

    DO k=1,2
      DO j=1,ny-1
        DO i=1,nx-1
          tem6(i,j,k)=3.9/MAX(lenscl(i,j,k),0.1*tem1(i,j,k))
        END DO
      END DO
    END DO

  ELSE IF (tkeopt == 3) THEN  ! Sun and Chang formulation

    DO k=1,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          tem6(i,j,k)=0.41/MAX(lenscl(i,j,k),0.1*tem1(i,j,k))
        END DO
      END DO
    END DO

  END IF
!
!-----------------------------------------------------------------------
!
!  Added the TKE dissipation term Ce*tke**(3/2)/lenscl to tkeforce.
!
!-----------------------------------------------------------------------
!
  DO k=1,nz-1
    DO j=1,ny-1
      DO i=1,nx-1
        tkeforce(i,j,k)=tkeforce(i,j,k) -rhostr(i,j,k)*tem6(i,j,k)*     &
                        tke(i,j,k,tstrtlvl)*SQRT(tke(i,j,k,tstrtlvl))
      END DO
    END DO
  END DO
!
!-----------------------------------------------------------------------
!
!  Add the shear production term into tkeforce. Note that the
!  divergence term is ignored.
!
!-----------------------------------------------------------------------
!
  DO k=1,nz-1
    DO j=1,ny-1
      DO i=1,nx-1
        tkeforce(i,j,k)=tkeforce(i,j,k) +                               &
                        rhostr(i,j,k)*kmv(i,j,k)*defsq(i,j,k)
      END DO
    END DO
  END DO
!
!-----------------------------------------------------------------------
!
!  Calculate the buoyancy production term in the TKE equation and
!  store the result in tem1, which is then added to tkeforce.
!
!  Note: The following lendel (the length scale normalized by
!          (dx*dy*dz)**(1/3)) is at time, tpresent.
!
!-----------------------------------------------------------------------
!
  it = tpresent

  CALL buoytke(nx,ny,nz,ptprt(1,1,1,it),pprt(1,1,1,it),                 &
               qv(1,1,1,it),qscalar(:,:,:,it,:),                        &
               ptbar,pbar,ptbari,rhostr,zp,j3,j3inv,kmv,                &
               rprntl,ptsflx,qvsflx,                                    &
               tem1,                                                    &
               tem2,tem3,tem4,tem5,tem6)

  DO k=1,nz-1
    DO j=1,ny-1
      DO i=1,nx-1
        tkeforce(i,j,k)=tkeforce(i,j,k) + tem1(i,j,k)
      END DO
    END DO
  END DO

  END IF ! Calculating additional forcing terms

  DO k=1,nz-1
    DO j=1,ny-1
      DO i=1,nx-1
        ! Make a copy of the total forcing in tem8, which is used in the calculation of vertically-implicit
        ! mixing.  tkeforce itself (all forcing terms except for advection) is not modified, so that it
        ! can be passed back into the subroutine for the last two RK steps (i.e. only advection/vertically-implicit mixing
        ! are updated for each RK3 step for tintegopt = 2)
        tem8(i,j,k) = tkeforce(i,j,k)+tkeadv(i,j,k)
      END DO
    END DO
  END DO

  dt2 = dtbig1*2.0  ! Not sure why this is here: it isn't used anywhere

!
!-----------------------------------------------------------------------
!
!  Treat the vertically implicit diffusion term
!  Note that the vertically-implicit mixing is updated for each RK3 step
!  since the advective forcings are updated, while the explicit mixing
!  is still held constant.
!
!-----------------------------------------------------------------------
!
  IF (trbvimp == 1) THEN     ! Vertical implicit application

    CALL set_acct(tmix_acct)

    DO k=1,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          tem1(i,j,k)=2.0*rhostr(i,j,k)*kmv(i,j,k)                      &
                      *j3inv(i,j,k)*j3inv(i,j,k)
        END DO
      END DO
    END DO

    CALL vmiximps(nx,ny,nz,deltat*0.5,tke(1,1,1,tstrtlvl),rhostr,       &
                  tem1,tem8,tem2,tem3,tem4,tem5)

  END IF

  ! Again, note that in the following, tkeforce contains all the forcing
  ! terms except advection, whereas tem8 contains the *total* forcing.
  ! tkeforce is either calculated above or
  ! passed into the subroutine from a prior call to solvtke, depending
  ! on the choice of tintegopt.

!
!-----------------------------------------------------------------------
!
!  Integrate the TKE equation forward by one timestep,
!  yielding TKE at  time = tfuture.
!
!-----------------------------------------------------------------------
!
  CALL set_acct(misc_acct)

  DO k=2,nz-2
    DO j=2,ny-2
      DO i=2,nx-2
        tke(i,j,k,tfuture)=tke(i,j,k,tstrtlvl) +                        &
                           deltat*tem8(i,j,k)*rhostri(i,j,k)
      END DO
    END DO
  END DO
!
!-----------------------------------------------------------------------
!
!  Update B.C.'s for TKE.
!
!-----------------------------------------------------------------------
!
  IF (mp_opt > 0) THEN
    CALL acct_interrupt(mp_acct)
    CALL mpsendrecv2dew(tke(1,1,1,tfuture),nx,ny,nz,ebc,wbc,0,tem1)
    CALL mpsendrecv2dns(tke(1,1,1,tfuture),nx,ny,nz,nbc,sbc,0,tem1)
  END IF
  CALL acct_interrupt(bc_acct)
  CALL bckmkh(nx,ny,nz,tke(1,1,1,tfuture))
  CALL acct_stop_inter
!
!-----------------------------------------------------------------------
!
!  Cut off the potentially negative values of TKE.
!
!-----------------------------------------------------------------------
!
  DO k=1,nz-1
    DO j=1,ny-1
      DO i=1,nx-1
        tke(i,j,k,tfuture)=MAX(0.0,tke(i,j,k,tfuture))
      END DO
    END DO
  END DO

  RETURN
END SUBROUTINE solvtke
!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINE BUOYTKE                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE buoytke(nx,ny,nz,ptprt,pprt,qv,qscalar,                      &
           ptbar,pbar,ptbari,rhostr,zp,j3,j3inv,kmv,rprntl,             &
           ptsflx,qvsflx,                                               &
           tkebuoy,                                                     &
           tem1,tem2,tem3,tem4,tem5)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Compute the buoyancy production term in the TKE equation.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: V.Wong, Y. Tang and X. Song
!  10/1993
!
!  MODIFICATION HISTORY:
!
!  9/1/94 (Y. Lu)
!  Cleaned up documentation.
!
!  11/17/1995 (Ming Xue)
!  Fixed an important bug in the loop 10, where J3 should be
!  divided not multiplied.
!
!  05/29/1999 (Ming Xue)
!  Re-written to incooperate the surface heat and moisture fluxes.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!    ptprt    Perturbation potential temperature at times tpast and
!             tpresent (K)
!    pprt     Perturbation pressure at times tpast and
!             tpresent (Pascal)
!
!    qv       Water vapor specific humidity at times tpast
!             and tpresent (kg /kg)
!
!    qscalar
!
!    ptbar    Base state potential temperature (K)
!    pbar     Base state pressure (Pascal)
!    ptbari   Inverse base state potential temperature (K)
!
!    rhostr   Base state density rhobar times j3 (kg/m**3)
!    zp       Vertical coordinate of grid points
!             in physical space (m)
!    j3       Coordinate transformation Jacobian  d(zp)/dz
!    j3inv    Inverse of j3
!
!    kmv      Vertical turb. mixing coef. for momentum ( m**2/s )
!    rprntl   Reciprocal of Prandtl number
!    ptsflx   Surface heat flux
!    qvsflx   Surface moisture flux
!
!  OUTPUT:
!
!    tkebuoy  Temporary work array for buoyancy production term
!
!  WORK ARRAYS:
!
!    tem1     Temporary work array
!    tem2     Temporary work array
!    tem3     Temporary work array
!    tem4     Temporary work array
!    tem5     Temporary work array
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE             ! Force explicit declarations
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'phycst.inc'
  INCLUDE 'globcst.inc'
  INCLUDE 'grid.inc'          ! Grid & map parameters.
!
!-----------------------------------------------------------------------
!
!  Variable Declarations:
!
!-----------------------------------------------------------------------
!
  INTEGER :: nx, ny, nz        ! Number of grid points in 3 directions

  REAL :: ptprt (nx,ny,nz)     ! Perturbation potential temperature (K)
  REAL :: pprt  (nx,ny,nz)     ! Perturbation pressure (Pascal)

  REAL :: qv    (nx,ny,nz)     ! Water vapor specific humidity (kg/kg)

  REAL :: qscalar(nx,ny,nz,nscalar)

  REAL :: ptbar (nx,ny,nz)     ! Base state potential temperature (K)
  REAL :: pbar  (nx,ny,nz)     ! Base state pressure (Pascal).
  REAL :: ptbari(nx,ny,nz)     ! Inverse base state pot. temperature (K)
  REAL :: rhostr(nx,ny,nz)     ! Base state density rhobar times j3.

  REAL :: zp    (nx,ny,nz)     ! Physical height
                               ! coordinate defined at
                               ! w-point of the staggered grid.

  REAL :: j3    (nx,ny,nz)     ! Coordinate transformation Jacobian
                               ! defined d( zp )/d( z ).
  REAL :: j3inv (nx,ny,nz)     ! Inverse of j3

  REAL :: kmv   (nx,ny,nz)     ! Vertical turb. mixing coef. for
                               ! momentum. ( m**2/s )
  REAL :: rprntl(nx,ny,nz)     ! Reciprocal of Prandtl number
  REAL :: ptsflx(nx,ny)        ! Surface heat flux (K*kg/(m*s**2))
  REAL :: qvsflx(nx,ny)        ! Surface moisture flux (kg/(m**2*s))

  REAL :: tkebuoy (nx,ny,nz)   ! Temporary work array for buoyancy
                               ! production term

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
  INTEGER :: i,j,k
  INTEGER :: nq
  INTEGER :: onvf
  REAL    :: lcp,eps,tema

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  eps=1.e-6
!
!-----------------------------------------------------------------------
!
!  Calculate avgz(rhobar*kh/j3)
!
!-----------------------------------------------------------------------
!
  DO k=1,nz-1
    DO j=1,ny-1
      DO i=1,nx-1
        tem5(i,j,k) = rhostr(i,j,k)*kmv(i,j,k)*rprntl(i,j,k)            &
                    * (j3inv(i,j,k)*j3inv(i,j,k))
      END DO
    END DO
  END DO

!-----------------------------------------------------------------------
!
!  Calculate rhobar*Khv/J3 *  difz( pt ).
!
!-----------------------------------------------------------------------
!
  DO k=2,nz-1
    DO j=1,ny-1
      DO i=1,nx-1
        tem1(i,j,k)=dzinv*((ptprt(i,j,k)+ptbar(i,j,k))-                 &
                           (ptprt(i,j,k-1)+ptbar(i,j,k-1)))
      END DO
    END DO
  END DO

  DO k=2,nz-1
    DO j=1,ny-1
      DO i=1,nx-1
        tem1(i,j,k)=tem1(i,j,k) *(tem5(i,j,k)+tem5(i,j,k-1))*0.5
      END DO
    END DO
  END DO

!-----------------------------------------------------------------------
!
!  Calculate rhobar*Khv/J3 * difz( qv ).
!
!-----------------------------------------------------------------------

  tema = 0.5*dzinv
  DO k=2,nz-1
    DO j=1,ny-1
      DO i=1,nx-1
        tem2(i,j,k)=tema*(qv(i,j,k)-qv(i,j,k-1))                        &
                        *(tem5(i,j,k)+tem5(i,j,k-1))
      END DO
    END DO
  END DO

  IF( sfcphy /= 0 ) THEN
    DO j=1,ny-1
      DO i=1,nx-1
        tem1(i,j,2)=ptsflx(i,j)
        tem2(i,j,2)=qvsflx(i,j)
      END DO
    END DO
  END IF

  DO k=2,nz-2
    DO j=1,ny-1
      DO i=1,nx-1
        tem3(i,j,k)=0.5*(tem1(i,j,k+1)+tem1(i,j,k))
        tem4(i,j,k)=0.5*(tem2(i,j,k+1)+tem2(i,j,k))
      END DO
    END DO
  END DO

!-----------------------------------------------------------------------
!
!    Calculate the buoyancy production term for the unsaturated case
!    in the TKE equation.
!
!-----------------------------------------------------------------------

  DO k=2,nz-2
    DO j=1,ny-1
      DO i=1,nx-1
        tkebuoy(i,j,k)=-g*j3(i,j,k)*                                    &
                        (tem3(i,j,k)*ptbari(i,j,k)+0.61*tem4(i,j,k))
      END DO
    END DO
  END DO

  IF( moist /= 0) THEN
!
!-----------------------------------------------------------------------
!
!  Compute the buoyancy production term in the saturated case where
!  Qv >= Qs (or Qc>0).
!
!-----------------------------------------------------------------------


!
!-----------------------------------------------------------------------
!
!  Calculate tem1 = d(eqivalent potential temperature)/dz
!  tem2 = total temperature,
!  tem3=(0.622*L*qv)/(R*T) is the equivalent potential temperature.
!
!-----------------------------------------------------------------------
!
    lcp=lathv/cp
    tema = 1.0/p0

    DO k=1,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          tem2(i,j,k)=(ptbar(i,j,k)+ptprt(i,j,k))*                      &
                      ((pbar(i,j,k)+pprt(i,j,k))*tema)**rddcp
        END DO
      END DO
    END DO

    DO k=1,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          tem3(i,j,k)=(3.376/tem2(i,j,k)-0.00254)*1000*qv(i,j,k)*       &
              (1+0.81*qv(i,j,k))
          tem3(i,j,k)=(ptprt(i,j,k)+ptbar(i,j,k))*EXP(tem3(i,j,k))
        END DO
      END DO
    END DO

    onvf = 1
    CALL difz(tem3, onvf,                                               &
              nx,ny,nz, 1,nx-1, 1,ny-1, 2,nz-1, dz, tem1)

!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!

    DO k=1,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          tem3(i,j,k)=0.622*lathv*qv(i,j,k)/(rd*tem2(i,j,k))
        END DO
      END DO
    END DO

    DO k=1,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          tem4(i,j,k)=(1.0+lcp*tem3(i,j,k)/tem2(i,j,k))
        END DO
      END DO
    END DO

    DO k=2,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          tem1(i,j,k)=tem1(i,j,k)*(tem5(i,j,k)+tem5(i,j,k-1))           &
                      /(tem4(i,j,k-1)+tem4(i,j,k))
        END DO
      END DO
    END DO

    IF( sfcphy /= 0 ) THEN
      DO j=1,ny-1
        DO i=1,nx-1
          tem1(i,j,2)=ptsflx(i,j)
        END DO
      END DO
    END IF

    DO k=1,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          tem2(i,j,k)=(1.0+0.61*tem3(i,j,k))*ptbari(i,j,k) ! testing
        END DO
      END DO
    END DO

    DO k=2,nz-2
      DO j=1,ny-1
        DO i=1,nx-1
          tem2(i,j,k)=tem2(i,j,k)*(tem1(i,j,k)+tem1(i,j,k+1))*0.5
        END DO
      END DO
    END DO

!-----------------------------------------------------------------------
!
!  Calculate qls=(qv+qc+qr+qi+qs+qg), and store the result in tem2,
!  which is the sum of all liquid and solid water. Then compute
!  the saturated part of the buoyancy production term.
!
!-----------------------------------------------------------------------
!
    tem3(:,:,:) = 0.0

    DO nq = 1,nscalarq

      DO k=1,nz-1
        DO j=1,ny-1
          DO i=1,nx-1
            tem3(i,j,k)= tem3(i,j,k) + qscalar(i,j,k,nq)
          END DO
        END DO
      END DO

    END DO
!
!-----------------------------------------------------------------------
!
!    Reset the buoyancy production term for saturated points.
!
!-----------------------------------------------------------------------
!
    tema = 0.5*dzinv

    DO k=2,nz-2
      DO j=1,ny-1
        DO i=1,nx-1
          IF( qscalar(i,j,k,P_QC) > eps) THEN   ! because it is inside if block moist /=0
                                                ! so P_QC > 0
            tkebuoy(i,j,k)=-g*j3(i,j,k)*                                &
                (tem2(i,j,k)-tem5(i,j,k)*tema*(tem3(i,j,k+1)-tem3(i,j,k-1)))
          END IF
        END DO
      END DO
    END DO

  END IF

  RETURN
END SUBROUTINE buoytke
!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINE BCKMKH                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################

SUBROUTINE bckmkh(nx,ny,nz,s)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Set the boundary conditions for scalar, s (which may be TKE). Zero
!  gradient is assumed except for the periodic boundary condition
!  case.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Vince Wong and X. Song
!  010/08/93
!
!  MODIFICATION HISTORY:
!
!  9/1/94 (Y. Lu)
!  Cleaned up documentation.
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
!    s        Array of any scalar.
!
!  OUTPUT:
!
!    s       Including boundary values
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx,ny,nz          ! Number of grid points in 3 directions

  REAL :: s   (nx,ny,nz)       ! Scalar array.
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
  INCLUDE 'bndry.inc'
  INCLUDE 'mp.inc'

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  IF (wbc == 2) THEN
    IF (mp_opt == 0) THEN
      DO k=2,nz-2
      DO j=1,ny-1
        s(1,j,k)=s(nx-2,j,k)
      END DO
      END DO
    END IF
  ELSE IF (wbc /= 0) THEN
    DO k=2,nz-2
      DO j=1,ny-1
        s(1,j,k)=s(2,j,k)
      END DO
    END DO
  END IF

  IF (ebc == 2) THEN
    IF (mp_opt == 0) THEN
      DO k=2,nz-2
      DO j=1,ny-1
        s(nx-1,j,k)=s(2,j,k)
      END DO
      END DO
    END IF
  ELSE IF (ebc /= 0) THEN
    DO k=2,nz-2
      DO j=1,ny-1
        s(nx-1,j,k)=s(nx-2,j,k)
      END DO
    END DO
  END IF

  IF (nbc == 2) THEN
    IF (mp_opt == 0) THEN
      DO k=2,nz-2
      DO i=1,nx-1
        s(i,ny-1,k)=s(i,2,k)
      END DO
      END DO
    END IF
  ELSE IF (nbc /= 0) THEN
    DO k=2,nz-2
      DO i=1,nx-1
        s(i,ny-1,k)=s(i,ny-2,k)
      END DO
    END DO
  END IF

  IF (sbc == 2) THEN
    IF (mp_opt == 0) THEN
      DO k=2,nz-2
      DO i=1,nx-1
        s(i,1,k)=s(i,ny-2,k)
      END DO
      END DO
    END IF
  ELSE IF (sbc /= 0) THEN
    DO k=2,nz-2
      DO i=1,nx-1
        s(i,1,k)=s(i,2,k)
      END DO
    END DO
  END IF

  IF (bbc == 2) THEN
    DO i=1,nx-1
      DO j=1,ny-1
        s(i,j,1)=s(i,j,nz-2)
      END DO
    END DO
  ELSE
    DO i=1,nx-1
      DO j=1,ny-1
        s(i,j,1)=s(i,j,2)
      END DO
    END DO
  END IF

  IF (tbc == 2) THEN
    DO i=1,nx-1
      DO j=1,ny-1
        s(i,j,nz-1)=s(i,j,2)
      END DO
    END DO
  ELSE
    DO i=1,nx-1
      DO j=1,ny-1
        s(i,j,nz-1)=s(i,j,nz-2)
      END DO
    END DO
  END IF

  RETURN
END SUBROUTINE bckmkh
