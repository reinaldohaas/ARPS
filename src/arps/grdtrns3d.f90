!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINE GRDTRAN                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE grdtran(nx,ny,nz,ubar,vbar,u,v,w,ptprt,pprt,                 &
           qv,qvbar,qscalar,rhostr,x,y,zp,j3,j3inv,                     &
           tem1,tem2,tem3,tem4,tem5,tem6,tem7,                          &
           tem8,tem9,tem10,tem11)
!
!--------------------------------------------------------------------------
!
!  PURPOSE:
!
!  Coordinate the calculation and adjustment of domain translation
!  speed using either cell-tracking or optimal pattern translation
!  algorithm.
!
!---------------------------------------------------------------------------
!
!  AUTHORS: Ming Xue.
!  3/29/1995.
!
!  MODIFICATION HISTORY:
!
!  04/18/1995 (Y. Liu)
!  Set the default values for variable uchange and vchange.
!
!  06/24/1995 (Alan Shapiro)
!  Documentation clean-up.
!
!  06/27/95 (M. Xue)
!  Added qvbar to argument list of GRDTRAN and ADJUVWV.
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the z-direction (vertical)
!
!    ubar     Model's base state x velocity component (m/s)
!    vbar     Model's base state y velocity component (m/s)
!
!    u        Total u-velocity (m/s)
!    v        Total v-velocity (m/s)
!    w        Total w-velocity (m/s)
!    ptprt    Perturbation potential temperature (K)
!    pprt     Perturbation pressure (Pascal)
!    qv       Water vapor specific humidity (kg/kg)
!    qc       Cloud water mixing ratio (kg/kg)
!    qr       Rain water mixing ratio (kg/kg)
!    qi       Cloud ice mixing ratio (kg/kg)
!    qs       Snow mixing ratio (kg/kg)
!    qh       Hail mixing ratio (kg/kg)
!
!    qvbar    Base-state water vapor mixing ratio (kg/kg)
!    rhostr   Base state density rhobar times j3 (kg/m**3)
!    x        x coordinate of grid points in physical/comp. space (m)
!    y        y coordinate of grid points in physical/comp. space (m)
!    zp       Vertical coordinate of grid points in physical space (m)
!    j3       Coordinate transformation Jacobian  d(zp)/dz
!
!  OUTPUT:
!
!    ubar     Redefined base state x velocity component (m/s)
!    vbar     Redefined base state y velocity component (m/s)
!
!    u        Redefined total u-velocity (m/s)
!    v        Redefined total v-velocity (m/s)
!    w        Redefined total w-velocity (m/s)
!    ptprt    Redefined perturbation potential temperature at tpast (K)
!    pprt     Redefined perturbation pressure at tpast (Pascal)
!    qv       Redefined water vapor specific humidity at tpast (kg/kg)
!    qc       Redefined cloud water mixing ratio at tpast (kg/kg)
!    qr       Redefined rain water mixing ratio at tpast (kg/kg)
!    qi       Redefined cloud ice mixing ratio at tpast (kg/kg)
!    qs       Redefined snow mixing ratio at tpast (kg/kg)
!    qh       Redefined hail mixing ratio at tpast (kg/kg)
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
  INCLUDE 'grid.inc'        ! Grid parameters
  INCLUDE 'timelvls.inc'
!
!-----------------------------------------------------------------------
!
!  Variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: nx,ny,nz         ! Number of grid points in x, y, z directions

  REAL :: ubar (nx,ny,nz)     ! Base state x velocity (m/s)
  REAL :: vbar (nx,ny,nz)     ! Base state y velocity (m/s)

  REAL :: u    (nx,ny,nz,nt)  ! Total u-velocity (m/s)
  REAL :: v    (nx,ny,nz,nt)  ! Total v-velocity (m/s)
  REAL :: w    (nx,ny,nz,nt)  ! Total w-velocity (m/s)
  REAL :: ptprt(nx,ny,nz,nt)  ! Perturbation potential temperature (K)
  REAL :: pprt (nx,ny,nz,nt)  ! Perturbation pressure (Pascal)
  REAL :: qv   (nx,ny,nz,nt)  ! Water vapor specific humidity (kg/kg)
  REAL :: qvbar(nx,ny,nz)     ! Base-state water vapor mixing ratio (kg/kg)

  REAL :: qscalar(nx,ny,nz,nt,nscalar)

  REAL :: rhostr(nx,ny,nz)    ! Base state density rhobar times j3.
  REAL :: x    (nx)           ! The x-coord. of the physical and
                              ! computational grid. Defined at u-point.
  REAL :: y    (ny)           ! The y-coord. of the physical and
                              ! computational grid. Defined at v-point.
  REAL :: zp    (nx,ny,nz)    ! The physical height coordinate defined at
                              ! w-point of staggered grid.

  REAL :: j3    (nx,ny,nz)    ! Coordinate transformation Jacobian defined
                              ! as d( zp )/d( z ).
  REAL :: j3inv (nx,ny,nz)    ! Coordinate transformation Jacobian defined
                              ! as d( zp )/d( z ).

  REAL :: tem1 (nx,ny,nz)     ! Temporary work array.
  REAL :: tem2 (nx,ny,nz)     ! Temporary work array.
  REAL :: tem3 (nx,ny,nz)     ! Temporary work array.
  REAL :: tem4 (nx,ny,nz)     ! Temporary work array.
  REAL :: tem5 (nx,ny,nz)     ! Temporary work array.
  REAL :: tem6 (nx,ny,nz)     ! Temporary work array.
  REAL :: tem7 (nx,ny,nz)     ! Temporary work array.
  REAL :: tem8 (nx,ny,nz)     ! Temporary work array.
  REAL :: tem9 (nx,ny,nz)     ! Temporary work array.
  REAL :: tem10(nx,ny,nz)     ! Temporary work array.
  REAL :: tem11(nx,ny,nz)     ! Temporary work array.

!-----------------------------------------------------------------------
!
! Misc. local variables
!
!-----------------------------------------------------------------------

  INTEGER :: i,j,k,ntrcell,ireturn
  REAL :: uvmvmag, scal ,ucelcntr,vcelcntr
  REAL :: uchange,vchange,xcctr_old,ycctr_old

  INTEGER, SAVE :: called
  REAL,    SAVE :: xcelcntr, ycelcntr
  DATA called /0/
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  uchange = 0.0
  vchange = 0.0

  IF(MOD(nstep,nceltrk) == 0 .AND.(cltkopt == 1.OR.grdtrns == 2))THEN
!
!-----------------------------------------------------------------------
!
!  Current model grid origin in absolute coord.
!
!-----------------------------------------------------------------------
!
    xgrdorg = xgrdorg + tceltrk * umove
    ygrdorg = ygrdorg + tceltrk * vmove

    DO k=1,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          tem11(i,j,k) = rhostr(i,j,k)*j3inv(i,j,k)
        END DO
      END DO
    END DO

    IF( called /= 0 ) THEN
      xcctr_old = xcelcntr
      ycctr_old = ycelcntr
    END IF
!
!-----------------------------------------------------------------------
!
!  Find new cell center.
!
!-----------------------------------------------------------------------
!
    IF (P_QR <= 0) THEN     ! make sure QR is defined
      WRITE(6,'(2a)') 'QR is not defined before calling celtrk. ',      &
                      'Please check mphyopt and try again.'
      CALL arpsstop('QR is not defined in grdtrns3d.f90.',1)
    END IF

    CALL celtrk(nx,ny,nz,2,nx-2,2,ny-2,2,nz-2,                          &
         w(1,1,1,tpresent), qscalar(:,:,:,tpresent,P_QR),               &
         tem11, x,y,zp, ntrcell, xcelcntr, ycelcntr,                    &
         ireturn,                                                       &
         tem1,tem2,tem3,tem4,tem5,tem6,tem7,tem8,tem9,tem10)
!
!-----------------------------------------------------------------------
!
!  Cell movement speed between the past two celtrk calls.
!
!-----------------------------------------------------------------------
!
    IF( ntrcell > 0 ) THEN
      IF( called /= 0 ) THEN
        ucelcntr = (xcelcntr-xcctr_old)/tceltrk
        vcelcntr = (ycelcntr-ycctr_old)/tceltrk
      ELSE
        ucelcntr = 0.0
        vcelcntr = 0.0
      END IF
      called = 1
    END IF

  END IF

  IF (grdtrns == 2 .AND. MOD(nstep,nceltrk) == 0) THEN
!
!-----------------------------------------------------------------------
!
!  Calculate adjustment to umove and vmove (uchange,vchange) which
!  tries to bring the center of mass of cells back
!  to the center of model domain during time period tcrestr
!  assuming the cell center moves at past speed.
!
!-----------------------------------------------------------------------
!
    IF( ireturn == 0 .AND. called /= 0 ) THEN
      !WYH??? x(1), x(nx), y(1), y(ny) should use global values?
      uchange = ( xcelcntr+ucelcntr*tcrestr                             &
              -((x(1)+x(nx))*0.5+xgrdorg+umove*tcrestr))/tcrestr
      vchange = ( ycelcntr+vcelcntr*tcrestr                             &
              -((y(1)+y(ny))*0.5+ygrdorg+vmove*tcrestr))/tcrestr

      uvmvmag = SQRT( uchange*uchange + vchange*vchange )

      IF( uvmvmag > 10.0 ) THEN
        scal = 10.0 /uvmvmag
        uchange = uchange * scal
        vchange = vchange * scal
      END IF

    ELSE

      uchange = 0.0
      vchange = 0.0

    END IF

  ELSE IF (grdtrns == 3) THEN
!
!-----------------------------------------------------------------------
!
!  Automatic domain translation.
!
!  Call AUTOTRANS to compute running mean of the optimal pattern
!  translation umove and vmove, and return the appropriate increments
!  (uchange, vchange) in umove and vmove at the end of time window.
!
!  The average domain translation speed in a given time window
!  (twindow) is calculated from the movement speed of traced properties
!  during each time step (autotrans is called every time step).
!  Only at the end of this time window are non-zero values of umove
!  and vmove sent out. At other times, they are zero, indicating
!  no change should be made to umove or vmove.
!
!-----------------------------------------------------------------------
!
    CALL autotrans(nx,ny,nz,                                            &
                   ubar,vbar,u,v,w,ptprt,pprt,qv,qscalar,               &
                   uchange, vchange, tem1)


  END IF
!
!-----------------------------------------------------------------------
!
!  Reset the model u and v velocity values using the new domain
!  translation speed.
!  (ADJUVMV should only be called when it's time to use it.)
!
!-----------------------------------------------------------------------
!

  IF( uchange /= 0. .OR. vchange /= 0.) THEN

    WRITE(6,'(3(/1x,a)/)')                                              &
        'ATTENTION: UMOVE or VMOVE has just been changed. ',            &
        'Subroutine ADJUVMV will now be called to adjust the time-',    &
        'dependent variables.'

    CALL adjuvmv(nx,ny,nz,                                              &
         ubar,vbar,u,v,w,ptprt,pprt,qv,qvbar,qscalar,                   &
         uchange, vchange, tem1, tem11)

  END IF


  RETURN

END SUBROUTINE grdtran

!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINE ADJUVMV                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE adjuvmv(nx,ny,nz,                                            &
                  ubar,vbar,u,v,w,ptprt,pprt,qv,qvbar,qscalar,          &
                  uchange,vchange , tem1, tem2)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Adjust the model variables per change in the domain translation
!  speed.
!
!-----------------------------------------------------------------------
!
!  AUTHORS: Ming Xue
!  Based on Alan Shapiro and Steve Lazarus' GRIDTRANS.
!  9/7/1993
!
!  MODIFICATION HISTORY:
!
!  06/27/95 (M. Xue)
!  Added qvbar to the argument list. bcqv called for qv.
!
!  11/20/97 (fanyou Kong - CMRP)
!  Add a do loop to remove possible negative Qx values after
!  calling galilei subroutine
!
!  7/16/98 (M. Xue and Y. Richardson)
!  Fixed a minor bug in call to bcp. It was introduced
!  when subroutine BCP was modified.
!
!  2000/04/21 (Gene Bassett)
!  Adjusted BC calls so that tfuture is update for all variables (rather
!  than u, v, w & pprt for tpast).
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the z-direction (vertical)
!
!    ubar     Model's base state x velocity component (m/s)
!    vbar     Model's base state y velocity component (m/s)
!
!    u        Total u-velocity (m/s)
!    v        Total v-velocity (m/s)
!    w        Total w-velocity (m/s)
!    ptprt    Perturbation potential temperature (K)
!    pprt     Perturbation pressure (Pascal)
!    qv       Water vapor specific humidity (kg/kg)
!    qc       Cloud water mixing ratio (kg/kg)
!    qr       Rain water mixing ratio (kg/kg)
!    qi       Cloud ice mixing ratio (kg/kg)
!    qs       Snow mixing ratio (kg/kg)
!    qh       Hail mixing ratio (kg/kg)
!    qvbar    Base-state water vapor mixing ratio (kg/kg)
!
!    uchange  Change in the domain translation speed in x-dir.
!    vchange  Change in the domain translation speed in y-dir.
!
!  OUTPUT:
!
!    ubar     Redefined base state x velocity component (m/s)
!    vbar     Redefined base state y velocity component (m/s)
!
!    u        Redefined total u-velocity (m/s)
!    v        Redefined total v-velocity (m/s)
!    w        Redefined total w-velocity (m/s)
!    ptprt    Redefined perturbation potential temperature at tpast (K)
!    pprt     Redefined perturbation pressure at tpast (Pascal)
!    qv       Redefined water vapor specific humidity at tpast (kg/kg)
!    qc       Redefined cloud water mixing ratio at tpast (kg/kg)
!    qr       Redefined rain water mixing ratio at tpast (kg/kg)
!    qi       Redefined cloud ice mixing ratio at tpast (kg/kg)
!    qs       Redefined snow mixing ratio at tpast (kg/kg)
!    qh       Redefined hail mixing ratio at tpast (kg/kg)
!
!  WORK ARRAYS:
!
!    tem1     Temporary work array.
!    tem2     Temporary work array.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE            ! Force explicit declarations
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'    ! Global constants that control model
  INCLUDE 'bndry.inc'      ! Boundary condition parameters
  INCLUDE 'mp.inc'         ! Message passing parameters.
  INCLUDE 'timelvls.inc'
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!

  INTEGER :: nx,ny,nz         ! Number of grid points in x, y, z directions

  REAL :: ubar (nx,ny,nz)     ! Base state x velocity (m/s)
  REAL :: vbar (nx,ny,nz)     ! Base state y velocity (m/s)

  REAL :: u    (nx,ny,nz,nt)  ! Total u-velocity (m/s)
  REAL :: v    (nx,ny,nz,nt)  ! Total v-velocity (m/s)
  REAL :: w    (nx,ny,nz,nt)  ! Total w-velocity (m/s)
  REAL :: ptprt(nx,ny,nz,nt)  ! Perturbation potential temperature (K)
  REAL :: pprt (nx,ny,nz,nt)  ! Perturbation pressure (Pascal)
  REAL :: qv   (nx,ny,nz,nt)  ! Water vapor specific humidity (kg/kg)
  REAL :: qvbar(nx,ny,nz)     ! Base-state water vapor mixing ratio (kg/kg)
  REAL :: qscalar(nx,ny,nz,nt,nscalar)

  REAL :: uchange             ! Change in domain translation speed in x-dir.
  REAL :: vchange             ! Change in domain translation speed in y-dir.

  REAL :: tem1 (nx,ny,nz)     ! Temporary work array.
  REAL :: tem2 (nx,ny,nz)     ! Temporary work array.
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!

  INTEGER :: i, j, k
  INTEGER :: nq

  INTEGER :: newebc           ! East boundary condition parameter.
  INTEGER :: newwbc           ! West boundary condition parameter.
  INTEGER :: newnbc           ! North boundary condition parameter.
  INTEGER :: newsbc           ! South boundary condition parameter.

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
!  Reset the domain translation speed if umove or vmove given in the
!  input file is different from those in the restart data.
!
!-----------------------------------------------------------------------
!
  WRITE(6,'((/1x,a,f10.3)/)')                                           &
      'ATTENTION: the domain translation speed was changed at t=',      &
      curtim

  WRITE(6,'(2(/1x,2(a,f10.4))/)')                                       &
      'umove_new=',umove   ,' vmove_new=',vmove  ,                      &
      'umove_old=',umove-uchange,' vmove_old=',vmove-vchange

  DO k=1,nz
    DO j=1,ny
      DO i=1,nx
        ubar(i,j,k)=ubar(i,j,k)-uchange
        vbar(i,j,k)=vbar(i,j,k)-vchange
        u(i,j,k,tpast   )=u(i,j,k,tpast   )-uchange
        u(i,j,k,tpresent)=u(i,j,k,tpresent)-uchange
        u(i,j,k,tfuture )=u(i,j,k,tfuture )-uchange
        v(i,j,k,tpast   )=v(i,j,k,tpast   )-vchange
        v(i,j,k,tpresent)=v(i,j,k,tpresent)-vchange
        v(i,j,k,tfuture )=v(i,j,k,tfuture )-vchange
      END DO
    END DO
  END DO
!
!-----------------------------------------------------------------------
!
!  Adjust tpast values of all the time dependent variables
!  such that their local time derivatives in the new moving reference
!  frame satisfies the Galilean transformation:
!
!  time derivative of a variable in the old system =
!
!  time derivative of the variable in the new system
!                    - uchange *d(var)/dx - vchange*d(var)/dy
!
!  where the time derivative is defined as
!
!  ( var( tpresent ) - var( tpast ) )/dt.
!
!-----------------------------------------------------------------------
!

  CALL galilei(nx,ny,nz, u,    2,nx-1,2,ny-2,1,nz-1,                    &
               uchange,vchange,tem1)
  CALL galilei(nx,ny,nz, v,    2,nx-2,2,ny-1,1,nz-1,                    &
               uchange,vchange,tem1)
  CALL galilei(nx,ny,nz, w,    2,nx-2,2,ny-2,1,nz  ,                    &
               uchange,vchange,tem1)
  CALL galilei(nx,ny,nz, ptprt,2,nx-2,2,ny-2,1,nz-1,                    &
               uchange,vchange,tem1)
  CALL galilei(nx,ny,nz, pprt, 2,nx-2,2,ny-2,1,nz-1,                    &
               uchange,vchange,tem1)
  CALL galilei(nx,ny,nz, qv,   2,nx-2,2,ny-2,1,nz-1,                    &
               uchange,vchange,tem1)
  DO nq = 1,nscalar
    CALL galilei(nx,ny,nz, qscalar(:,:,:,:,nq),2,nx-2,2,ny-2,1,nz-1,    &
                 uchange,vchange,tem1)
  END DO
!
!****************************************************************
!  Remove possible negative values for Q-fields resulted from
!  the calling of 'galilei' subroutine
!****************************************************************
!
  DO i=2,nx-2
    DO j=2,ny-2
      DO k=1,nz-1
        qv(i,j,k,1) = AMAX1(qv(i,j,k,1), 0.0)
        DO nq = 1,nscalar
          qscalar(i,j,k,1,nq) = AMAX1(qscalar(i,j,k,1,nq), 0.0)
        END DO
      END DO
    END DO
  END DO
!
!-----------------------------------------------------------------------
!
!  Apply boundary conditions to redefined dynamical, thermodynamical
!  and microphysical variables at tpast. Note: there is no
!  provision for radiation lateral boundary conditions; if radiation
!  conditions are chosen, zero-gradient conditions are applied instead.
!
!-----------------------------------------------------------------------
!

  IF ((ebc /= 1).AND.(ebc /= 2).AND.(ebc /= 3).AND.(ebc /= 0)) THEN
    newebc = 3
  ELSE
    newebc = ebc
  END IF

  IF ((wbc /= 1).AND.(wbc /= 2).AND.(wbc /= 3).AND.(wbc /= 0)) THEN
    newwbc = 3
  ELSE
    newwbc = wbc
  END IF

  IF ((nbc /= 1).AND.(nbc /= 2).AND.(nbc /= 3).AND.(nbc /= 0)) THEN
    newnbc = 3
  ELSE
    newnbc = nbc
  END IF

  IF ((sbc /= 1).AND.(sbc /= 2).AND.(sbc /= 3).AND.(sbc /= 0)) THEN
    newsbc = 3
  ELSE
    newsbc = sbc
  END IF

  IF (mp_opt > 0) THEN
    CALL acct_interrupt(mp_acct)
    CALL mpsendrecv2dew(u,nx,ny,nz,newebc,newwbc,1,tem2)
    CALL mpsendrecv2dns(u,nx,ny,nz,newnbc,newsbc,1,tem2)
  END IF
  CALL acct_interrupt(bc_acct)
  CALL bcu(nx,ny,nz,0., u(1,1,1,tfuture),                               &
           0.,0.,0.,0.,newebc,newwbc,newnbc,newsbc,tbc,bbc,             &
           ebc_global,wbc_global,nbc_global,sbc_global)
  CALL acct_stop_inter

  IF (mp_opt > 0) THEN
    CALL acct_interrupt(mp_acct)
    CALL mpsendrecv2dew(v,nx,ny,nz,newebc,newwbc,2,tem2)
    CALL mpsendrecv2dns(v,nx,ny,nz,newnbc,newsbc,2,tem2)
  END IF
  CALL acct_interrupt(bc_acct)
  CALL bcv(nx,ny,nz,0., v(1,1,1,tfuture),                               &
           0.,0.,0.,0.,newebc,newwbc,newnbc,newsbc,tbc,bbc,             &
           ebc_global,wbc_global,nbc_global,sbc_global)
  CALL acct_stop_inter

  IF (mp_opt > 0) THEN
    CALL acct_interrupt(mp_acct)
    CALL mpsendrecv2dew(tem1,nx,ny,nz,newebc,newwbc,3,tem2)
    CALL mpsendrecv2dns(tem1,nx,ny,nz,newnbc,newsbc,3,tem2)
  END IF
  CALL acct_interrupt(bc_acct)
  CALL lbcw(nx,ny,nz,0.0,w,tem1,                                        &
            0.,0.,0.,0.,newebc,newwbc,newnbc,newsbc,                    &
            ebc_global,wbc_global,nbc_global,sbc_global)
  CALL acct_stop_inter

  IF (mp_opt > 0) THEN
    CALL acct_interrupt(mp_acct)
    CALL mpsendrecv2dew(ptprt,nx,ny,nz,newebc,newwbc,0,tem2)
    CALL mpsendrecv2dns(ptprt,nx,ny,nz,newnbc,newsbc,0,tem2)
  END IF
  CALL acct_interrupt(bc_acct)
  CALL bcpt(nx,ny,nz,0.,ptprt,                                          &
            0.,0.,0.,0.,newebc,newwbc,newnbc,newsbc,tbc,bbc,            &
            ebc_global,wbc_global,nbc_global,sbc_global)
  CALL acct_stop_inter

  DO j=1,ny-1
    DO i=1,nx-1
      tem1(i,j,1)=pprt(i,j,2,1)
    END DO
  END DO
  IF (mp_opt > 0) THEN
    CALL acct_interrupt(mp_acct)
    CALL mpsendrecv2dew(pprt,nx,ny,nz,newebc,newwbc,0,tem2)
    CALL mpsendrecv2dns(pprt,nx,ny,nz,newnbc,newsbc,0,tem2)
  END IF
  CALL acct_interrupt(bc_acct)
  CALL bcp(nx,ny,nz,0.,                                                 &
           pprt(1,1,1,tfuture),0.,0.,0.,0.,tem1(1,1,1),                 &
           newebc,newwbc,newnbc,newsbc,tbc,bbc,                         &
           ebc_global,wbc_global,nbc_global,sbc_global)
  CALL acct_stop_inter

  IF (mp_opt > 0) THEN
    CALL acct_interrupt(mp_acct)
    CALL mpsendrecv2dew(qv,nx,ny,nz,newebc,newwbc,0,tem2)
    CALL mpsendrecv2dns(qv,nx,ny,nz,newnbc,newsbc,0,tem2)
  END IF
  CALL acct_interrupt(bc_acct)
  CALL bcqv(nx,ny,nz,0., qv,qvbar,                                      &
           0.,0.,0.,0.,newebc,newwbc,newnbc,newsbc,tbc,bbc,             &
           ebc_global,wbc_global,nbc_global,sbc_global)
  CALL acct_stop_inter

  DO nq = 1,nscalar

    IF (mp_opt > 0) THEN
      CALL acct_interrupt(mp_acct)
      CALL mpsendrecv2dew(qscalar(:,:,:,1,nq),nx,ny,nz,newebc,newwbc,0,tem2)
      CALL mpsendrecv2dns(qscalar(:,:,:,1,nq),nx,ny,nz,newnbc,newsbc,0,tem2)
    END IF
    CALL acct_interrupt(bc_acct)
    CALL bcq(nx,ny,nz,0., qscalar(:,:,:,1,nq),                          &
             0.,0.,0.,0.,newebc,newwbc,newnbc,newsbc,tbc,bbc,           &
             ebc_global,wbc_global,nbc_global,sbc_global)
    CALL acct_stop_inter

  END DO

  RETURN
END SUBROUTINE adjuvmv

!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINE GALILEI                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE galilei(nx,ny,nz, var, ibgn,iend,jbgn,jend,kbgn,kend,        &
           uchange,vchange, tem1)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Adjust tpast values of all the time dependent variables
!  such that their local time derivatives in the new moving reference
!  frame satisfies the Galilean transformation:
!
!  time derivative of a variable in the old system =
!
!  time derivative of the variable in the new system
!                    - vchange*d(var)/dx - vchange*d(var)/dy
!
!  where the time derivative is defined as
!
!  ( var( tpresent ) - var( tpast ) )/dt.
!
!-----------------------------------------------------------------------
!

!
!-----------------------------------------------------------------------
!
!
!  AUTHORS: Alan Shapiro and Steve Lazarus
!  5/23/93.
!
!  MODIFICATION HISTORY:
!  9/7/93 (MX)
!  Code cleanup with modifications.
!
!  3/25/94 (MX)
!  A bug fixed to remove the dependency of var(tpart,i) on
!  var(tpast,n-i) in loop 100.
!
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the z-direction (vertical)
!
!    var      Any one of the four-dimensional variable that is advected
!             using leapfrog-centered scheme.
!    uchange  New-old domain translation speed in x direction
!    vchange  New-old domain translation speed in y direction
!    ibgn     i-index where multiplication begins.
!    iend     i-index where multiplication ends.
!    jbgn     j-index where multiplication begins.
!    jend     j-index where multiplication ends.
!    kbgn     k-index where multiplication begins.
!    kend     k-index where multiplication ends.
!
!  OUTPUT:
!
!    var      var at tpast.
!
!  WORK ARRAY:
!
!    tem1     Work array.
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
  IMPLICIT NONE         ! Force explicit declarations

  INCLUDE 'timelvls.inc'

  INTEGER :: nx,ny,nz      ! Number of grid points in x, y, z directions

  REAL :: var(nx,ny,nz,nt) ! Any one of the advect variables.

  INTEGER :: ibgn,iend,jbgn,jend,kbgn,kend

  REAL :: uchange          ! New-old domain translation speed in x-dir
  REAL :: vchange          ! New-old domain translation speed in y-dir

  REAL :: tem1(nx,ny,nz)   ! Local work array.
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i, j, k
  REAL :: dvardx, dvardy
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'     ! Global constants that control model
  INCLUDE 'grid.inc'          ! Grid parameters

!-----------------------------------------------------------------------
!
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

  DO k=kbgn,kend
    DO j=jbgn,jend
      DO i=ibgn,iend

        dvardx = 0.25*dxinv*                                            &
                    (var(i+1,j,k,tpast) + var(i+1,j,k,tpresent)         &
                   - var(i-1,j,k,tpast) - var(i-1,j,k,tpresent))

        dvardy = 0.25*dyinv*                                            &
                    (var(i,j+1,k,tpast) + var(i,j+1,k,tpresent)         &
                   - var(i,j-1,k,tpast) - var(i,j-1,k,tpresent))

        tem1(i,j,k) = var(i,j,k,tpast)                                  &
                   - dtbig*(uchange*dvardx+vchange*dvardy)

      END DO
    END DO
  END DO

  DO k=kbgn,kend
    DO j=jbgn,jend
      DO i=ibgn,iend
        var(i,j,k,tpast) = tem1(i,j,k)
      END DO
    END DO
  END DO

  RETURN
END SUBROUTINE galilei
!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINE AUTOTRANS                 ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE autotrans(nx,ny,nz,                                          &
           ubar,vbar,u,v,w,ptprt,pprt,qv,qscalar,                       &
           uchange, vchange, tem1)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  This subroutine computes the optimum (least squares) pattern translation
!  speed based on vertical velocity data. This technique is described in:
!
!  T. Gal-Chen, 1982: "Errors in Fixed and Moving Frame of References:
!          Applications for Conventional and Doppler Radar Analysis",
!          J. Atmos. Sci., Vol. 39, 1982. 2279-2300.
!
!
!  Special note:  Many of the arrays in the AUTOTRANS subroutine call
!  are not actually used in AUTOTRANS.  However, these arrays will
!  likely be used in future tests of this subroutine.
!
!-----------------------------------------------------------------------
!
!  AUTHORS: Alan Shapiro and Steve Lazarus
!  5/02/93.
!
!  MODIFICATION HISTORY:
!
!  4/18/94  Alan Shapiro and Robb Randall
!  Changed weighting so that stronger storm contributes more to
!  the calculation of bigu and bigv.
!  Also, only regions of updraft (w>0) are used in the calculation.
!
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the z-direction (vertical)
!
!    ubar     Model's base state x velocity component (m/s)
!    vbar     Model's base state y velocity component (m/s)
!
!    u        Total u-velocity (m/s)
!    v        Total v-velocity (m/s)
!    w        Total w-velocity (m/s)
!    ptprt    Perturbation potential temperature (K)
!    pprt     Perturbation pressure (Pascal)
!    qv       Water vapor specific humidity (kg/kg)
!    qc       Cloud water mixing ratio (kg/kg)
!    qr       Rain water mixing ratio (kg/kg)
!    qi       Cloud ice mixing ratio (kg/kg)
!    qs       Snow mixing ratio (kg/kg)
!    qh       Hail mixing ratio (kg/kg)
!
!
!  OUTPUT:
!
!    uchange  Change in umove (=0 if it's not time to change).
!    vchange  Change in vmove (=0 if it's not time to change).
!
!
!  WORK ARRAYS:
!
!    tem1     Temporary work array.
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
  INCLUDE 'globcst.inc'     ! Global constants that control model
  INCLUDE 'grid.inc'          ! Grid parameters
  INCLUDE 'timelvls.inc'
  INCLUDE 'mp.inc'
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  INTEGER :: nx,ny,nz          ! Number of grid points in x, y, z directions

  REAL :: ubar  (nx,ny,nz)     ! Base state x velocity (m/s)
  REAL :: vbar  (nx,ny,nz)     ! Base state y velocity (m/s)

  REAL :: u     (nx,ny,nz,nt)  ! Total u-velocity (m/s)
  REAL :: v     (nx,ny,nz,nt)  ! Total v-velocity (m/s)
  REAL :: w     (nx,ny,nz,nt)  ! Total w-velocity (m/s)
  REAL :: ptprt (nx,ny,nz,nt)  ! Perturbation potential temperature (K)
  REAL :: pprt  (nx,ny,nz,nt)  ! Perturbation pressure (Pascal)
  REAL :: qv    (nx,ny,nz,nt)  ! Water vapor specific humidity (kg/kg)
  REAL :: qscalar(nx,ny,nz,nt,nscalar)

  REAL :: tem1  (nx,ny,nz)     ! Temporary work array.
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!

  REAL :: w1111, w1110, w2111,                                          &
          w0111, w1211, w1011  ! Vertical velocity.

  INTEGER :: i, j, k
  INTEGER :: kmin              ! Minimum k-level for calculation of optimal
                               ! grid-translation.
  INTEGER :: kmax              ! Maximum k-level for calculation of optimal
                               ! grid-translation.

  REAL :: dwdt, dwdx, dwdy     ! Spatial and temporal derivatives of w.

  REAL :: bigu                 ! Optimum perturbation translation speed in
                               ! x direction (computed each time step).

  REAL :: bigv                 ! Optimum perturbation translation speed in
                               ! y direction (computed each time step).

  REAL :: speed                ! speed = sqrt(bigu**2 + bigv**2)

  REAL :: a, b, c, d, e        ! Coefficients in solution for bigu and bigv.

  REAL :: dtinv                ! Reciprocal of dtbig
  REAL :: denom                ! Denominator of formula for bigu, bigv

  INTEGER :: ntlev             ! Number of big time steps in every time
                               ! window (rounded down)

  REAL :: velmax               ! Maximum allowable perturbation translation
                               ! speed.

  REAL :: uchange              ! Optimum perturbation translation speed in
                               ! x direction (averaged over a time window).

  REAL :: vchange              ! Optimum perturbation translation speed in
                               ! y direction (averaged over a time window).

  REAL :: sumu                 ! Running sum of bigu
  REAL :: sumv                 ! Running sum of bigv
  SAVE sumu, sumv

  INTEGER :: ichnge            ! Parameter that signifies whether grid
                               ! translation will be changed this time step.
  SAVE ichnge

  DATA sumu /0./
  DATA sumv /0./
  DATA ichnge /0/
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
!  Initialize uchange and vchange to zero.
!
!-----------------------------------------------------------------------
!
  uchange = 0.
  vchange = 0.
!
!-----------------------------------------------------------------------
!
!  Define maximum value of grid translation in terms of CFL criterion.
!
!-----------------------------------------------------------------------
!
  dtinv = 1./dtbig
  velmax = dtinv*SQRT(dx*dy)

!-----------------------------------------------------------------------
!
!  Compute optimal grid translation components based on low-level
!  vertical velocity data.
!
!-----------------------------------------------------------------------
!
  a = 0.
  b = 0.
  c = 0.
  d = 0.
  e = 0.

  kmin = 3
  kmax = 2 + chkdpth/dz

  DO k=kmin, kmax
    DO j=2,ny-2
      DO i=2,nx-2
        w1111 = MAX(0.0, w(i,j,k,tpresent))
        w1110 = MAX(0.0, w(i,j,k,tpast))
        w2111 = MAX(0.0, w(i+1,j,k,tpresent))
        w0111 = MAX(0.0, w(i-1,j,k,tpresent))
        w1211 = MAX(0.0, w(i,j+1,k,tpresent))
        w1011 = MAX(0.0, w(i,j-1,k,tpresent))

        dwdt = dtinv*(w1111**4 - w1110**4)
        dwdx = 0.5*dxinv*(w2111**4 - w0111**4)
        dwdy = 0.5*dyinv*(w1211**4 - w1011**4)

        a = a + dwdt*dwdx
        b = b + dwdx*dwdx
        c = c + dwdx*dwdy
        d = d + dwdt*dwdy
        e = e + dwdy*dwdy
      END DO
    END DO
  END DO

  CALL mpsumr(a,1)
  CALL mpsumr(b,1)
  CALL mpsumr(c,1)
  CALL mpsumr(d,1)
  CALL mpsumr(e,1)

  denom = c*c - b*e

  IF (denom == 0.) THEN
    PRINT *, 'Warning: denom = 0. bigu, bigv set to 0'
    bigu = 0.
    bigv = 0.
  ELSE
    bigu = (a*e - c*d)/denom
    bigv = (b*d - c*a)/denom
    PRINT *, 'bigu, bigv = ', bigu, bigv
  END IF

!
!-----------------------------------------------------------------------
!
!  Check that bigu and bigv are not too big.  If the speed of
!  translation exceeds velmax (based on CFL criterion) then set
!  bigu and bigv to zero.
!
!-----------------------------------------------------------------------
!
  speed = SQRT(bigu*bigu + bigv*bigv)

  IF (speed > velmax) THEN
    bigu = 0.
    bigv = 0.
    PRINT *, 'Warning! speed= ', speed, ' exceeds velmax= ', velmax
    PRINT *, 'Set bigu and bigv to zero.'
    PRINT *, 'bigu, bigv = ', bigu, bigv
  END IF

!
!-----------------------------------------------------------------------
!
!  Sum bigu and bigv each time step.
!
!-----------------------------------------------------------------------
!

  sumu = sumu + bigu
  sumv = sumv + bigv

!
!-----------------------------------------------------------------------
!
!  Update umove and vmove at the end of every time window.
!
!-----------------------------------------------------------------------
!

  ntlev = IFIX(twindow/dtbig)
  ichnge = MOD(nstep, ntlev)

  IF (ichnge == 0) THEN

    uchange = sumu/ntlev
    vchange = sumv/ntlev

    PRINT *, 'umove and vmove are being updated at time =', curtim

    PRINT *, 'Old umove and vmove = ', umove, vmove
    PRINT *, 'ntlev, uchange, vchange = ', ntlev, uchange, vchange

    umove = umove + uchange
    vmove = vmove + vchange

    PRINT *, 'New umove and vmove will be  = ', umove, vmove

!
!-----------------------------------------------------------------------
!
!    Set sumu and sumv to zero for the beginning of next time window.
!
!-----------------------------------------------------------------------
!

    sumu = 0.
    sumv = 0.

  END IF

  RETURN
END SUBROUTINE autotrans
