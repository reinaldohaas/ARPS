!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE SOLVUV                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE solvuv(mptr, nx,ny,nz, exbcbufsz, dtsml1, curtsml,           &
           u,v,wcont,pprt,                                              &
           udteb,udtwb,udtnb,udtsb,                                     &
           vdteb,vdtwb,vdtnb,vdtsb,                                     &
           rhostr, uforce,vforce,ubar,vbar,                             &
           x,y,z,zp, mapfct, j1,j2,j3,j3inv,                            &
           exbcbuf,rhofct,                                              &
           rhostru,rhostrv,rhostrw,dtmrxrsu,dtmryrsv,wpgrad,            &
           upgrad,vpgrad,tem1,tem2,tem3)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Coordinate the time stepping of u and v momentum equations and the
!  setting of boundary conditions.
!
!-----------------------------------------------------------------------
!
!
!  AUTHOR: Ming Xue
!  10/10/92
!
!  MODIFICATION HISTORY:
!
!  5/20/92 (K. Droegemeier and M. Xue)
!  Added full documentation.
!
!  6/07/92 ( M. Xue)
!  Now only tfuture fields of u, v, w and pprt are passed in.
!
!  2/10/93 (K. Droegemeier)
!  Cleaned up documentation.
!
!  4/10/93 (M. Xue & Hao Jin)
!  Add the terrain.
!
!  9/1/94 (Y. Lu)
!  Cleaned up documentation.
!
!  1/25/96 (Donghai Wang & Yuhe Liu)
!  Added the map projection factor to ARPS governing equations.
!
!  10/16/97 (Donghai Wang)
!  Using total density (rho) for the calculation of the pressure
!  gradient force terms.
!
!  11/06/97 (Dan Weber)
!  Added three additional levels to the mapfct array.  The three
!  levels (4,5,6) represent the inverse of the first three in order.
!  The inverse map factors are computed to improve efficiency.
!
!  9/28/98 (Dan Weber)
!  Reformulated do-loops to improve efficiency, brought in
!  pre-computed variables.  Note u,vforce contains dtsml and
!  1/rhostru,v.
!
!  07/23/2001 (K. Brewster)
!  Added mptr to argument list.
!  Added optional writing of total divergence as a diagnostic noise
!  parameter.
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!    dtsml1   Local value of small time step size
!   curtsml   Local curttim at a small time step
!
!    u        x component of velocity at tfuture (m/s)
!    v        y component of velocity at tfuture (m/s)
!    wcont    Contravariant vertical velocity (m/s)
!
!    pprt     Perturbation pressure at tfuture (Pascal)
!
!    udteb    Time tendency of the u field at the east boundary
!    udtwb    Time tendency of the u field at the west boundary
!    udtnb    Time tendency of the u field at the north boundary
!    udtsb    Time tendency of the u field at the south boundary
!
!    vdteb    Time tendency of the v field at the east boundary
!    vdtwb    Time tendency of the v field at the west boundary
!    vdtnb    Time tendency of the v field at the north boundary
!    vdtsb    Time tendency of the v field at the south boundary
!
!    rhostr   Base state density rhobar times j3 (kg/m**3)
!
!    uforce   Acoustically inactive forcing terms in u-eq.
!             (kg/(m*s)**2)
!    vforce   Acoustically inactive forcing terms in v-eq.
!             (kg/(m*s)**2)
!
!    ubar     Base state x velocity component (m/s)
!    vbar     Base state y velocity component (m/s)
!
!    x        x coordinate of grid points in physical/comp. space (m)
!    y        y coordinate of grid points in physical/comp. space (m)
!    z        z coordinate of grid points in computational space (m)
!    zp       Vertical coordinate of grid points in physical space(m)
!
!    mapfct   Map factors at scalar, u and v points
!
!    j1       Coordinate transformation Jacobian -d(zp)/dx
!    j2       Coordinate transformation Jacobian -d(zp)/dy
!    j3       Coordinate transformation Jacobian  d(zp)/dz
!    j3inv    Inverse of j3
!
!    rhostru  Average rhostr at u points (kg/m**3).
!    rhostrv  Average rhostr at v points (kg/m**3).
!    rhostrw  Average rhostr at w points (kg/m**3).
!
!  OUTPUT:
!
!    u        x component of velocity at time tfuture updated
!             in time by dtsml1 (m/s)
!    v        y component of velocity at time tfuture updated
!             in time by dtsml1 (m/s)
!    w        Vertical component of Cartesian velocity at tfuture
!             updated in time by dtsml1 (m/s)
!    wpgrad   Pressure gradient term in w-eq. (kg/(m*s)**2)
!
!  WORK ARRAYS:
!
!    upgrad   Pressure gradient term in u-eq. (kg/(m*s)**2)
!    vpgrad   Pressure gradient term in v-eq. (kg/(m*s)**2)
!    tem1     Temporary work array
!    tem2     Temporary work array
!    tem3     Temporary work array
!    dtmrxrsu dtsml*mapfct*avgx(rhofct)/rhostru
!    dtmryrsv dtsml*mapfct*avgy(rhofct)/rhostrv
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

  INTEGER :: mptr

  INTEGER :: nx, ny, nz        ! Number of grid points in 3 directions

  REAL :: u     (nx,ny,nz)     ! u-velocity at tfuture (m/s)
  REAL :: v     (nx,ny,nz)     ! v-velocity at tfuture (m/s)
  REAL :: wcont (nx,ny,nz)     ! Contravariant vertical velocity (m/s)

  REAL :: pprt  (nx,ny,nz)     ! Perturbation pressure at tfuture
                               ! (Pascal)

  REAL :: udteb (ny,nz)        ! Time tendency of u at east boundary
  REAL :: udtwb (ny,nz)        ! Time tendency of u at west boundary
  REAL :: udtnb (nx,nz)        ! Time tendency of u at north boundary
  REAL :: udtsb (nx,nz)        ! Time tendency of u at south boundary

  REAL :: vdteb (ny,nz)        ! Time tendency of v at east boundary
  REAL :: vdtwb (ny,nz)        ! Time tendency of v at west boundary
  REAL :: vdtnb (nx,nz)        ! Time tendency of v at north boundary
  REAL :: vdtsb (nx,nz)        ! Time tendency of v at south boundary

  REAL :: rhostr(nx,ny,nz)     ! Base state density rhobar times j3.
  REAL :: uforce(nx,ny,nz)     ! Acoustically inactive forcing terms
                               ! in u-momentum equation (kg/(m*s)**2)
  REAL :: vforce(nx,ny,nz)     ! Acoustically inactive forcing terms
                               ! in v-momentum equation (kg/(m*s)**2)
  REAL :: ubar(nx,ny,nz)
  REAL :: vbar(nx,ny,nz)

  REAL :: x     (nx)           ! x-coord. of the physical and compu-
                               ! tational grid. Defined at u-point.
  REAL :: y     (ny)           ! y-coord. of the physical and compui-
                               ! tational grid. Defined at v-point.
  REAL :: z     (nz)           ! z-coord. of the computational grid.
                               ! Defined at w-point on the staggered
                               ! grid.
  REAL :: zp    (nx,ny,nz)     ! Physical height coordinate defined at
                               ! w-point of the staggered grid.

  REAL :: mapfct(nx,ny,8)      ! Map factors at scalar, u and v points

  REAL :: j1    (nx,ny,nz)     ! Coordinate transformation Jacobian
                               ! defined as - d( zp )/d( x ).
  REAL :: j2    (nx,ny,nz)     ! Coordinate transformation Jacobian
                               ! defined as - d( zp )/d( y ).
  REAL :: j3    (nx,ny,nz)     ! Coordinate transformation Jacobian
                               ! defined as d( zp )/d( z ).
  REAL :: j3inv (nx,ny,nz)     ! Inverse of j3

  INTEGER :: exbcbufsz         ! EXBC buffer size
  REAL :: exbcbuf( exbcbufsz ) ! EXBC buffer array

  REAL :: rhofct(nx,ny,nz)     ! rho-factor: rhobar/rho

  REAL :: rhostru(nx,ny,nz)    ! Averaged rhostr at u points (kg/m**3).
  REAL :: rhostrv(nx,ny,nz)    ! Averaged rhostr at v points (kg/m**3).
  REAL :: rhostrw(nx,ny,nz)    ! Averaged rhostr at w points (kg/m**3).
  REAL :: wpgrad(nx,ny,nz)     ! Pressure gradient term in w-eq.
!
!-----------------------------------------------------------------------
!
!  Temporary WORK ARRAYS:
!
!-----------------------------------------------------------------------
!
  REAL :: upgrad(nx,ny,nz)     ! Pressure gradient term in u-eq.
  REAL :: vpgrad(nx,ny,nz)     ! Pressure gradient term in v-eq.
  REAL :: tem1  (nx,ny,nz)     ! Temporary work array
  REAL :: tem2  (nx,ny,nz)     ! Temporary work array
  REAL :: tem3  (nx,ny,nz)     ! Temporary work array
  REAL :: dtmrxrsu(nx,ny,nz)   ! dtsml*map*avgx(rhofct)/rhostru
  REAL :: dtmryrsv(nx,ny,nz)   ! dtsml*map*avgy(rhofct)/rhostrv
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  REAL :: curtsml           ! Local curttim at a small time step
  REAL :: dtsml1            ! Local small time step
  REAL :: sumdiv
  INTEGER :: i, j, k
  INTEGER :: istart,iend,jstart,jend
  INTEGER :: ebc1,wbc1,nbc1,sbc1
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
  INCLUDE 'exbc.inc'
  INCLUDE 'mp.inc'            ! Message passing parameters.

  INTEGER :: noisewrt
  INTEGER :: ncalls(mgrdmax), nchdiv1(mgrdmax)
  CHARACTER (LEN=256 ) :: divfn
  INTEGER :: ldivfn,istat

  SAVE ncalls,nchdiv1

  DATA ncalls /mgrdmax*0/

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
  noisewrt=0

  ncalls(mptr)=ncalls(mptr) + 1

  IF(noisewrt == 1 .AND. ncalls(mptr) == 1 ) THEN

      divfn  = runname(1:lfnkey)//'.absdiv'
      ldivfn = 7 + lfnkey

      IF(nestgrd == 1) THEN
        WRITE(divfn((ldivfn+1):(ldivfn+4)),'(a,i2.2)')'.g',mptr
        ldivfn =  ldivfn + 4
      END IF

      WRITE(6,'(1x,a,a,a/,1x,a)')                                       &
          'Check to see if file ',divfn(1:ldivfn),' already exists.',   &
          'If so, append a version number to the filename.'

      CALL fnversn( divfn, ldivfn )

      CALL getunit ( nchdiv1(mptr) )

      OPEN(nchdiv1(mptr),FORM='formatted',STATUS='new',                 &
                FILE=divfn(1:ldivfn),IOSTAT=istat)

      IF( istat /= 0) THEN

        WRITE(6,'(/a,i2,/a/)')                                          &
            ' Error occured when opening file '//divfn(1:ldivfn)//      &
            ' using FORTRAN unit ',nchdiv1(mptr),                       &
            ' Program stopped in MAXMIN.'
        CALL arpsstop(' ',1)

      END IF

      WRITE(nchdiv1(mptr),'(a)') runname
      WRITE(nchdiv1(mptr),'(t2,a,t15,a)') 'Time_(s)','MnAbsDiv'

  END IF

!-----------------------------------------------------------------------
!
!  Divergence damping is activated if divdmp = 1.
!  Otherwise, this portion of the code is skipped to save
!  computations.
!
!-----------------------------------------------------------------------

  IF( divdmp == 1 .OR. divdmp == 2 ) THEN
!
!-----------------------------------------------------------------------
!
!  The divergence damping terms are defined as
!
!    d(cdvdmph*div)/dx for u eq.,
!    d(cdvdmph*div)/dy for v eq.,
!    d(cdvdmpv*div)/dz for w eq.,
!
!  where
!
!    div = (difx(u*rhostr) + dify(v*rhostr) + difz(wcont*rhostr))/j3
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!  Compute the momentum divergence term, defined as
!
!  div = d(u*rhostr)/dx + d(v*rhostr)/dy + d(wcont*rhostr)/dz.
!
!  and combine the pressure and divergence damping into one
!  array so that the pressure gradient force and divergence damping
!  can be calculated in a single step.
!
!-----------------------------------------------------------------------

    DO k=1,nz-1
      DO j=1,ny-1
        DO i=1,nx
          upgrad(i,j,k)=u(i,j,k)*rhostru(i,j,k)*mapfct(i,j,5)
        END DO
      END DO
    END DO

    DO k=1,nz-1
      DO j=1,ny
        DO i=1,nx-1
          vpgrad(i,j,k)=v(i,j,k)*rhostrv(i,j,k)*mapfct(i,j,6)
        END DO
      END DO
    END DO

    DO k=1,nz
      DO j=1,ny-1
        DO i=1,nx-1
          tem1(i,j,k)=wcont(i,j,k)*rhostrw(i,j,k)
        END DO
      END DO
    END DO

    IF(noisewrt == 1 .AND. MOD(ncalls(mptr),2) == 0) THEN
      IF(lbcopt == 2) THEN
        istart=MAX(2,(ngbrz+1))
        iend=MIN((nx-2),(nx-ngbrz-1))
        jstart=MAX(2,(ngbrz+1))
        jend=MIN((ny-2),(ny-ngbrz-1))
      ELSE
        istart=2
        iend=nx-2
        jstart=2
        jend=ny-2
      END IF
      sumdiv=0.
      DO k=2,nz-2
        DO j=jstart,jend
          DO i=istart,iend
            sumdiv = sumdiv+abs(j3inv(i,j,k)                          &
                      * ( mapfct(i,j,7)                                 &
                        * ((upgrad(i+1,j,k)-upgrad(i,j,k))*dxinv        &
                          +(vpgrad(i,j+1,k)-vpgrad(i,j,k))*dyinv)       &
                        + (tem1(i,j,k+1)-tem1(i,j,k))*dzinv ))
          END DO
        END DO
      END DO
      sumdiv=1.0E06*sumdiv/                                           &
             float((iend-istart+1)*(jend-jstart+1)*(nz-3))
      write(nchdiv1(mptr),'(f10.2,f12.4)') (curtim+2*dtsml1),sumdiv
    END IF

    DO k=1,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          tem3(i,j,k) = j3inv(i,j,k)                                    &
                      * ( mapfct(i,j,7)                                 &
                        * ((upgrad(i+1,j,k)-upgrad(i,j,k))*dxinv        &
                          +(vpgrad(i,j+1,k)-vpgrad(i,j,k))*dyinv)       &
                        + (tem1(i,j,k+1)-tem1(i,j,k))*dzinv )
        END DO
      END DO
    END DO

  ELSE

    DO k=1,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          tem3(i,j,k)= 0.0
        END DO
      END DO
    END DO

  END IF

!-----------------------------------------------------------------------
!
!  Compute the pressure gradient force terms for the three
!  momentum equations and store the values in upgrad, vpgrad,
!  and wpgrad.
!
!  When divergence damping is activated (divdmp=1), these arrays also
!  contain the divergence damping terms.
!
!-----------------------------------------------------------------------

  CALL pgrad(nx,ny,nz, pprt, tem3,                                      &
             j1,j2,j3, upgrad,vpgrad,wpgrad,tem1,tem2)

!-----------------------------------------------------------------------
!
!  Integrate the u, and v equations forward by one small
!  timestep dtsml1 using a forward-in-time integration scheme
!  (that is, forward relative to the pressure gradient terms).
!
!-----------------------------------------------------------------------

  IF( nestgrd == 1 .AND. mgrid /= 1 ) THEN

    DO k=2,nz-2
      DO j=2,ny-2
        DO i=2,nx-1
          u(i,j,k)=u(i,j,k) + uforce(i,j,k) -                           &
                              dtmrxrsu(i,j,k)*upgrad(i,j,k)
        END DO
      END DO
    END DO

    DO k=2,nz-2
      DO j=2,ny-1
        DO i=2,nx-2
          v(i,j,k)=v(i,j,k) + vforce(i,j,k) -                           &
                              dtmryrsv(i,j,k)*vpgrad(i,j,k)
        END DO
      END DO
    END DO

!-----------------------------------------------------------------------
!
!  Set the lateral boundary conditions for u and v given the
!  time tendencies in the case of nested inner grid.
!
!-----------------------------------------------------------------------

    DO k=2,nz-2
      DO j=1,ny-1
        u( 1,j,k)=u( 1,j,k)+dtsml1 * udtwb(j,k)
        u(nx,j,k)=u(nx,j,k)+dtsml1 * udteb(j,k)
      END DO
    END DO
    DO k=2,nz-2
      DO i=2,nx-1
        u(i,   1,k)=u(i,   1,k)+dtsml1 * udtsb(i,k)
        u(i,ny-1,k)=u(i,ny-1,k)+dtsml1 * udtnb(i,k)
      END DO
    END DO

    DO k=2,nz-2
      DO j=1,ny
        v(   1,j,k)=v(   1,j,k)+dtsml1 * vdtwb(j,k)
        v(nx-1,j,k)=v(nx-1,j,k)+dtsml1 * vdteb(j,k)
      END DO
    END DO
    DO k=2,nz-2
      DO i=2,nx-2
        v(i, 1,k)=v(i, 1,k)+dtsml1 * vdtsb(i,k)
        v(i,ny,k)=v(i,ny,k)+dtsml1 * vdtnb(i,k)
      END DO
    END DO

!-----------------------------------------------------------------------
!
!  Set the top and bottom boundary conditions for u and v.
!
!-----------------------------------------------------------------------


! bc's for mp not implemented for nested grids

    CALL acct_interrupt(bc_acct)
    CALL bcu(nx,ny,nz, dtsml1, u, udteb,udtwb,udtnb,udtsb,              &
             0,0,0,0 ,tbc,bbc,                                          &
             ebc_global,wbc_global,nbc_global,sbc_global)
    CALL acct_stop_inter

    CALL acct_interrupt(bc_acct)
    CALL bcv(nx,ny,nz, dtsml1, v, vdteb,vdtwb,vdtnb,vdtsb,              &
             0,0,0,0 ,tbc,bbc,                                          &
             ebc_global,wbc_global,nbc_global,sbc_global)
    CALL acct_stop_inter

  ELSE

    DO k=2,nz-2
      DO j=1,ny-1
        DO i=2,nx-1
          u(i,j,k)=u(i,j,k) + uforce(i,j,k) -                           &
                              dtmrxrsu(i,j,k)*upgrad(i,j,k)
        END DO
      END DO
    END DO

    DO k=2,nz-2
      DO j=2,ny-1
        DO i=1,nx-1
          v(i,j,k)=v(i,j,k) + vforce(i,j,k) -                           &
                              dtmryrsv(i,j,k)*vpgrad(i,j,k)
        END DO
      END DO
    END DO

!-----------------------------------------------------------------------
!
!  Set the boundary conditions for u and v.
!
!-----------------------------------------------------------------------

    IF ( lbcopt == 1 ) THEN

      ebc1=ebc
      wbc1=wbc
      sbc1=sbc
      nbc1=nbc
      IF( ebc == 4 )  ebc1=0
      IF( wbc == 4 )  wbc1=0
      IF( nbc == 4 )  nbc1=0
      IF( sbc == 4 )  sbc1=0

      IF (mp_opt > 0) THEN
        CALL acct_interrupt(mp_acct)
        CALL mpsendrecv2dew(u,nx,ny,nz,ebc,wbc,1,tem1)
        CALL mpsendrecv2dns(u,nx,ny,nz,nbc1,sbc1,1,tem1)
      END IF
      CALL acct_interrupt(bc_acct)
      CALL bcu(nx,ny,nz, dtsml1, u, udteb,udtwb,udtnb,udtsb,            &
               ebc,wbc,nbc1,sbc1,tbc,bbc,                               &
               ebc_global,wbc_global,nbc_global,sbc_global)
      CALL acct_stop_inter

      IF (mp_opt > 0) THEN
        CALL acct_interrupt(mp_acct)
        CALL mpsendrecv2dew(v,nx,ny,nz,ebc1,wbc1,2,tem1)
        CALL mpsendrecv2dns(v,nx,ny,nz,nbc,sbc,2,tem1)
      END IF
      CALL acct_interrupt(bc_acct)
      CALL bcv(nx,ny,nz, dtsml1, v, vdteb,vdtwb,vdtnb,vdtsb,            &
               ebc1,wbc1,nbc,sbc,tbc,bbc,                               &
               ebc_global,wbc_global,nbc_global,sbc_global)
      CALL acct_stop_inter

    ELSE

      IF (mp_opt > 0) THEN
        CALL acct_interrupt(mp_acct)
        CALL mpsendrecv2dew(u,nx,ny,nz,0,0,1,tem1)
        CALL mpsendrecv2dns(u,nx,ny,nz,0,0,1,tem1)
      END IF
      CALL acct_interrupt(bc_acct)
      CALL bcu(nx,ny,nz, dtsml1, u, udteb,udtwb,udtnb,udtsb,            &
               0,0,0,0,tbc,bbc,                                         &
               ebc_global,wbc_global,nbc_global,sbc_global)
      CALL acct_stop_inter

      IF (mp_opt > 0) THEN
        CALL acct_interrupt(mp_acct)
        CALL mpsendrecv2dew(v,nx,ny,nz,0,0,2,tem1)
        CALL mpsendrecv2dns(v,nx,ny,nz,0,0,2,tem1)
      END IF
      CALL acct_interrupt(mp_acct)
      CALL bcv(nx,ny,nz, dtsml1, v, vdteb,vdtwb,vdtnb,vdtsb,            &
               0,0,0,0,tbc,bbc,                                         &
               ebc_global,wbc_global,nbc_global,sbc_global)

      CALL exbcuv(nx,ny,nz, curtsml, u,v,                               &
                  exbcbuf(nu0exb), exbcbuf(nv0exb),                     &
                  exbcbuf(nudtexb),exbcbuf(nvdtexb))
      CALL acct_stop_inter
    END IF

  END IF

  RETURN
END SUBROUTINE solvuv
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE SOLVWPEX                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE solvwpex(mptr, nx,ny,nz,exbcbufsz, dtsml1, curtsml,          &
           u,v,w,wcont,ptprt,pprt,phydro,                               &
           wdteb,wdtwb,wdtnb,wdtsb,pdteb,pdtwb,pdtnb,pdtsb,             &
           rhostr,ptbar,ptbari,pbari,csndsq,                            &
           wforce,wpgrad,pforce,                                        &
           x,y,z,zp, mapfct, j1,j2,j3,aj3x,aj3y,aj3z,j3inv,             &
           rhostru,rhostrv,rhostrw,                                     &
           exbcbuf,rhofct,                                              &
           div,pdiv,tem1,tem2,tem3,rstpbi,rstptbi,                      &
           dtrzrstw,dxj3xm,dyj3ym,rstj3i)

!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Coordinate the time stepping of the w momentum and pressure equations
!  using explicit time integration scheme. Divergence damping and the
!  acoustically inactive forcing terms in w and pressure equations
!  have been evaluated prior to this subroutine and are passed into
!  this routine.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  10/10/91
!
!  MODIFICATION HISTORY:
!
!  5/20/92 (K. Droegemeier and M. Xue)
!  Added full documentation.
!
!  6/07/92 ( M. Xue)
!  Now only tfuture fields of pprt, u, v and w are passed in.
!
!  2/10/93 (K. Droegemeier)
!  Cleaned up documentation.
!
!  4/10/93 (M. Xue & Hao Jin)
!  Add the terrain.
!
!  9/6/94 (M.Xue)
!  Pressure perturbation buoyancy and base state pressure
!  advection terms moved into the small time steps.
!
!  9/1/94 (Y. Lu)
!  Cleaned up documentation.
!
!  01/28/95 (G. Bassett)
!  Option added to turn off buoyancy terms (for buoyopt=0).
!
!  1/25/96 (Donghai Wang & Yuhe Liu)
!  Added the map projection factor to ARPS governing equations.
!
!  4/27/1996 (M. Xue)
!  Added code for peqopt=2 case.
!
!  5/6/96 (M. Xue)
!  Replaced csndsq in pressure buoyancy term by cpdcv*pbar/rhobar.
!
!  10/16/97 (Donghai Wang)
!  Using total density (rho) for the calculation of the pressure
!  gradient force terms.
!
!  10/16/97 (Donghai Wang)
!  Added the second order terms in the linerized buoyancy terms.
!
!  11/05/97 (D. Weber)
!  Added phydro array for use in the bottom boundary condition for
!  perturbation pressure (hydrostatic).
!
!  11/06/97 (D. Weber)
!  Added three additional levels to the mapfct array.  The three
!  levels (4,5,6) represent the inverse of the first three in order.
!  The inverse map factors are computed to improve efficiency.
!
!  9/18/98 (D. Weber)
!  Added arrays aj3x,y,z,pbari,ptbari,rstpbi,rstptbi,dtrzrstw,dxj3xm,
!  dyj3ym,rstj3i.
!  do not alter arrays!!!
!
!  07/23/2001 (K. Brewster)
!  Added mptr to argument list.
!  Added optional writing of mean abs(dp/dt) as a diagnostic noise
!  parameter.
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!    dtsml1   Local value of small time step size
!   curtsml   Local curttim at a small time step
!
!    u        x component of velocity at tfuture (m/s)
!    v        y component of velocity at tfuture (m/s)
!    w        Vertical component of Cartesian velocity at tfuture
!             (m/s)
!    wcont    Contravariant vertical velocity (m/s)
!    ptprt    Perturbation potential temperature at all time levels
!             (K)
!    pprt     Perturbation pressure at tfuture (Pascal)
!
!    phydro   Big time step forcing term for use in computing the
!             hydrostatic pressure at k=1.
!
!    wdteb    Time tendency of the w field at the east boundary
!    wdtwb    Time tendency of the w field at the west boundary
!    wdtnb    Time tendency of the w field at the north boundary
!    wdtsb    Time tendency of the w field at the south boundary
!
!    pdteb    Time tendency of the pressure field at the east
!             boundary
!    pdtwb    Time tendency of the pressure field at the west
!             boundary
!    pdtnb    Time tendency of the pressure field at the north
!             boundary
!    pdtsb    Time tendency of the pressure field at the south
!             boundary
!
!    rhostr   Base state density rhobar times j3 (kg/m**3)
!    ptbar    Base state potential temperature (K)
!    ptbari   Inverse base state potential temperature (K)
!    pbari    Inverse base state pressure (Pascal)
!    csndsq   Sound wave speed squared.
!
!    wforce   Acoustically inactive forcing terms in w-eq.
!             (kg/(m*s)**2)
!    pforce   Acoustically inactive terms in pressure eq. (Pascal/s)
!
!    x        x coordinate of grid points in physical/comp. space (m)
!    y        y coordinate of grid points in physical/comp. space (m)
!    z        z coordinate of grid points in computational space (m)
!    zp       Vertical coordinate of grid points in physical space(m)
!
!    mapfct   Map factors at scalar, u and v points
!
!    j1       Coordinate transformation Jacobian -d(zp)/dx
!    j2       Coordinate transformation Jacobian -d(zp)/dy
!    j3       Coordinate transformation Jacobian  d(zp)/dz
!    aj3x     Avgx of the coordinate transformation Jacobian  d(zp)/dz
!    aj3y     Avgy of the coordinate transformation Jacobian  d(zp)/dz
!    aj3z     Avgz of the coordinate transformation Jacobian  d(zp)/dz
!
!    rhostru  Average rhostr at u points (kg/m**3).
!    rhostrv  Average rhostr at v points (kg/m**3).
!    rhostrw  Average rhostr at w points (kg/m**3).
!
!  OUTPUT:
!
!    pprt     Perturbation pressure at tfuture updated in time
!             by stsml1 (Pascal)
!    w        Vertical component of Cartesian velocity at tfuture
!             updated in time by dtsml1 (m/s)
!
!  WORK ARRAYS:
!
!    wpgrad   Pressure gradient term in w-eq. (kg/(m*s)**2)
!    div      Velocity divergence (1/s), a local array
!    pdiv     Divergence term in pressure eq.
!             (Pascal/s), a local array
!    tem1     Temporary work array
!    tem2     Temporary work array
!    tem3     Temporary work array
!    rstpbi   rhostr*pbari
!    rstptbi  rhostr*ptbari
!    dtrzrstw dtsml*avgz(rhofct)/rhostrw
!    dxj3xm   dxinv*aj3x*mapfct(i,j,5)
!    dyj3ym   dyinv*aj3y*mapfct(i,j,6)
!    rstj3i   rhostr*j3inv
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

  INTEGER :: mptr

  INTEGER :: nx, ny, nz        ! Number of grid points in 3 directions

  REAL :: curtsml              ! Local curttim at a small time step
  REAL :: dtsml1               ! Local small time step size (s)


  REAL :: u     (nx,ny,nz)     ! u-velocity at tfuture (m/s)
  REAL :: v     (nx,ny,nz)     ! v-velocity at tfuture (m/s)
  REAL :: w     (nx,ny,nz)     ! w-velocity at tfuture (m/s)
  REAL :: wcont (nx,ny,nz)     ! Contravariant vertical velocity (m/s)
  REAL :: ptprt (nx,ny,nz)     ! Perturbation potential temperature (K)
  REAL :: pprt  (nx,ny,nz)     ! Perturbation pressure at tfuture
                               ! (Pascal)

  REAL :: phydro(nx,ny)        ! Big time step forcing for computing
                               ! hydrostatic pprt at k=1.

  REAL :: wdteb (ny,nz)        ! Time tendency of w at east boundary
  REAL :: wdtwb (ny,nz)        ! Time tendency of w at west boundary
  REAL :: wdtnb (nx,nz)        ! Time tendency of w at north boundary
  REAL :: wdtsb (nx,nz)        ! Time tendency of w at south boundary

  REAL :: pdteb (ny,nz)        ! Time tendency of pressure at east
                               ! boundary
  REAL :: pdtwb (ny,nz)        ! Time tendency of pressure at west
                               ! boundary
  REAL :: pdtnb (nx,nz)        ! Time tendency of pressure at north
                               ! boundary
  REAL :: pdtsb (nx,nz)        ! Time tendency of pressure at south
                               ! boundary

  REAL :: ptbar (nx,ny,nz)     ! Base state potential temperature (K)
  REAL :: ptbari(nx,ny,nz)     ! Inverse base state pot. temperature (K)
  REAL :: pbari (nx,ny,nz)     ! Inverse base state pressure (Pascal).
  REAL :: csndsq(nx,ny,nz)     ! Sound wave speed squared.
  REAL :: rhostr(nx,ny,nz)     ! Base state density rhobar times j3.
  REAL :: wforce(nx,ny,nz)     ! Acoustically inactive forcing terms
                               ! in w-momentum equation (kg/(m*s)**2)
  REAL :: pforce(nx,ny,nz)     ! Acoustically inactive forcing terms
                               ! in the pressure equation (Pascal/s)

  REAL :: x     (nx)           ! x-coord. of the physical and compu-
                               ! tational grid. Defined at u-point.
  REAL :: y     (ny)           ! y-coord. of the physical and compu-
                               ! tational grid. Defined at v-point.
  REAL :: z     (nz)           ! z-coord. of the computational grid.
                               ! Defined at w-point on the staggered
                               ! grid.
  REAL :: zp    (nx,ny,nz)     ! Physical height coordinate defined at
                               ! w-point of the staggered grid.

  REAL :: mapfct(nx,ny,8)      ! Map factors at scalar, u and v points

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
  REAL :: aj3z  (nx,ny,nz)     ! Coordinate transformation Jacobian defined
                               ! as d( zp )/d( z ) AVERAGED IN THE Z-DIR.
  REAL :: j3inv (nx,ny,nz)     ! Coordinate transformation Jacobian
                               ! defined as d( zp )/d( z ).

  INTEGER :: exbcbufsz         ! EXBC buffer size
  REAL :: exbcbuf( exbcbufsz ) ! EXBC buffer array

  REAL :: rhofct(nx,ny,nz)    ! rho-factor:rhobar/rho

  REAL :: rhostru(nx,ny,nz)    ! Averaged rhostr at u points (kg/m**3).
  REAL :: rhostrv(nx,ny,nz)    ! Averaged rhostr at v points (kg/m**3).
  REAL :: rhostrw(nx,ny,nz)    ! Averaged rhostr at w points (kg/m**3).
!
!-----------------------------------------------------------------------
!
!  Temporary work arrays:
!
!-----------------------------------------------------------------------
!
  REAL :: wpgrad(nx,ny,nz)     ! Pressure gradient term in w-eq.
  REAL :: div   (nx,ny,nz)     ! Velocity divergence (1/s)
  REAL :: pdiv  (nx,ny,nz)     ! Divergence term in pressure eq.
                               ! (Pascal/s)
  REAL :: tem1  (nx,ny,nz)     ! Temporary work array
  REAL :: tem2  (nx,ny,nz)     ! Temporary work array
  REAL :: tem3  (nx,ny,nz)     ! Temporary work array
  REAL :: rstpbi(nx,ny,nz)     ! rhostr*pbari
  REAL :: rstptbi(nx,ny,nz)    ! rhostr*ptbari
  REAL :: dtrzrstw(nx,ny,nz)   ! dtsml1*avgz(rhofct)/rhostrw
  REAL :: dxj3xm(nx,ny,nz)     ! dxinv*aj3x*mapfct(i,j,5)
  REAL :: dyj3ym(nx,ny,nz)     ! dyinv*aj3y*mapfct(i,j,6)
  REAL :: rstj3i(nx,ny,nz)     ! rhostr*j3inv
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'bndry.inc'
  INCLUDE 'exbc.inc'
  INCLUDE 'globcst.inc'
  INCLUDE 'grid.inc'          ! Grid & map parameters.
  INCLUDE 'phycst.inc'
  INCLUDE 'mp.inc'            ! Message passing parameters.
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i, j, k
  INTEGER :: istart,iend,jstart,jend
  INTEGER :: ebc1,wbc1,nbc1,sbc1
  INTEGER :: wpprt ! Switch for pprt/(rhobar*csndsq) term in w-eq.
  INTEGER :: prgw  ! Switch for rhobar*g*w term in w-eq.
  INTEGER :: noisewrt
  REAL :: prgwg5, g05
  REAL :: pttem,ptem,tem,tema,temb,temc
  REAL :: sumdpdt

  CHARACTER (LEN=256) :: dpdtfn
  INTEGER :: ldpdtfn,istat

  INTEGER :: ncalls(mgrdmax), nchdpdt1(mgrdmax)
  SAVE ncalls,nchdpdt1
  DATA ncalls /mgrdmax*0/
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

  noisewrt=0

  ncalls(mptr)=ncalls(mptr) + 1

  IF (noisewrt == 1 .AND. ncalls(mptr) == 1) THEN
      dpdtfn  = runname(1:lfnkey)//'.dpdt'
      ldpdtfn = 5 + lfnkey

      IF(nestgrd == 1) THEN
        WRITE(dpdtfn((ldpdtfn+1):(ldpdtfn+4)),'(a,i2.2)')'.g',mptr
        ldpdtfn =  ldpdtfn + 4
      END IF

      WRITE(6,'(1x,a,a,a/,1x,a)')                                       &
          'Check to see if file ',dpdtfn(1:ldpdtfn),' already exists.', &
          'If so, append a version number to the filename.'

      CALL fnversn( dpdtfn, ldpdtfn )

      CALL getunit ( nchdpdt1(mptr) )

      OPEN(nchdpdt1(mptr),FORM='formatted',STATUS='new',                &
                FILE=dpdtfn(1:ldpdtfn),IOSTAT=istat)

      IF( istat /= 0) THEN

        WRITE(6,'(/a,i2,/a/)')                                          &
            ' Error occured when opening file '//dpdtfn(1:ldpdtfn)//    &
            ' using FORTRAN unit ',nchdpdt1(mptr),                      &
            ' Program stopped in SOLVWPEX.'
        CALL arpsstop('arpsstop called from SOLVWPEX error opening '//  &
             'file',1)

      END IF

      WRITE(nchdpdt1(mptr),'(a)') runname
      WRITE(nchdpdt1(mptr),'(t2,a,t15,a)') 'Time (s)','MnAbsdpdt'

  END IF

  g05 = g*0.5

  IF( bsnesq == 1 ) THEN
    wpprt = 0   ! Switch for pprt/(rhobar*csndsq) term in w-eq.
    prgw  = 0   ! Switch for rhobar*g*w term in w-eq.
  ELSE
    wpprt = 1   ! Switch for pprt/(rhobar*csndsq) term in w-eq.
    prgw  = 1   ! Switch for rhobar*g*w term in w-eq.
  END IF

!
!-----------------------------------------------------------------------
!
!  The contribution from ptprt to the buoyancy term is calculated
!  inside the small time steps if the potential temperature equation
!  is solved inside small time steps, i.e., if ptsmlstp=1.
!
!  The contribution from pressure perturbation to the buoyancy is
!  always calculated inside the small time steps for better
!  computational stability
!
!  If buoyopt = 0 then turn off all buoyancy.
!
!-----------------------------------------------------------------------
!
  tem = 0.5*(1.0-cpdcv)
  tema = 1.0/cpdcv

  IF ( buoyopt /= 0 ) THEN
    CALL acct_interrupt(buoy_acct)

    IF ( ptsmlstp == 1 ) THEN

      IF ( buoy2nd == 0) THEN  !1st-order

        DO k=1,nz-1
          DO j=1,ny-1
            DO i=1,nx-1
              tem2(i,j,k)=rstptbi(i,j,k)*ptprt(i,j,k)                   &
                         -wpprt*pprt(i,j,k)*tema*rstpbi(i,j,k)
            END DO
          END DO
        END DO

      ELSE         !2nd-order

        DO k=1,nz-1
          DO j=1,ny-1
            DO i=1,nx-1
              pttem = ptprt(i,j,k)*ptbari(i,j,k)
              ptem  = wpprt*pprt(i,j,k)*(tema*pbari(i,j,k))
              tem2(i,j,k)=rhostr(i,j,k)*(                               &
                          pttem-pttem*pttem-ptem-tem*ptem*ptem          &
                          + 0.5*pttem*ptem)
            END DO
          END DO
        END DO

      END IF

    ELSE

      IF (buoy2nd == 0) THEN  !1st-order

        DO k=1,nz-1
          DO j=1,ny-1
            DO i=1,nx-1
              tem2(i,j,k)=                                              &
                  -wpprt*pprt(i,j,k)*rstpbi(i,j,k)*tema
            END DO
          END DO
        END DO

      ELSE     !2nd-order

        DO k=1,nz-1
          DO j=1,ny-1
            DO i=1,nx-1
              ptem  = wpprt*pprt(i,j,k)*(tema*pbari(i,j,k))
              tem2(i,j,k)=-rhostr(i,j,k)*(ptem+tem*ptem*ptem )
            END DO
          END DO
        END DO
      END IF

    END IF

    DO k=2,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          tem1(i,j,k)=(tem2(i,j,k)+tem2(i,j,k-1))*g05
        END DO
      END DO
    END DO

    CALL acct_stop_inter
  ELSE

    DO k=1,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          tem1(i,j,k) = 0.0
        END DO
      END DO
    END DO

  END IF
!
!-----------------------------------------------------------------------
!
!  To assume mirror symmetry about the top and bottom boundary,
!  the density/temperature has zero gradient across the boundary
!  but g is assumed to change sign, so that the buoyancy is zero
!  at the boundary. c.f. subroutine BUOYCY.
!
!-----------------------------------------------------------------------
!
  DO j=1,ny-1
    DO i=1,nx-1
      tem1(i,j,2)=0.0
      tem1(i,j,nz-1)=0.0
    END DO
  END DO

!-----------------------------------------------------------------------
!
!  Integrate w equations forward by one small timestep dtsml1 using
!  a forward-in-time integration scheme (that is, forward relative
!  to the pressure gradient terms).
!
!-----------------------------------------------------------------------

  IF( nestgrd == 1 .AND. mgrid /= 1 ) THEN

    DO k=2,nz-1
      DO j=2,ny-2
        DO i=2,nx-2
          w(i,j,k)=w(i,j,k) + wforce(i,j,k) -                           &
                   (wpgrad(i,j,k)-tem1(i,j,k))*dtrzrstw(i,j,k)
        END DO
      END DO
    END DO
!
!-----------------------------------------------------------------------
!
!  Set the lateral boundary conditions for w given the time tendencies
!  in the case of nested inner grid.
!
!-----------------------------------------------------------------------
!
    DO k=3,nz-2
      DO j=1,ny-1
        w(   1,j,k)=w(   1,j,k)+dtsml1 * wdtwb(j,k)
        w(nx-1,j,k)=w(nx-1,j,k)+dtsml1 * wdteb(j,k)
      END DO
    END DO
    DO k=3,nz-2
      DO i=2,nx-2
        w(i,   1,k)=w(i,   1,k)+dtsml1 * wdtsb(i,k)
        w(i,ny-1,k)=w(i,ny-1,k)+dtsml1 * wdtnb(i,k)
      END DO
    END DO

  ELSE

    DO k=2,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          w(i,j,k)=w(i,j,k) + wforce(i,j,k) -                           &
                   (wpgrad(i,j,k)-tem1(i,j,k))*dtrzrstw(i,j,k)
        END DO
      END DO
    END DO
!
!-----------------------------------------------------------------------
!
!  Apply the lateral boundary conditions for w.
!
!  For the open boundary case, w at the lateral boundarues are solved
!  from w equation directly therefore does not need to be reset.
!
!-----------------------------------------------------------------------
!

    IF( lbcopt == 1 ) THEN ! Internal boundary conditions

      ebc1=ebc
      wbc1=wbc
      sbc1=sbc
      nbc1=nbc
      IF( ebc == 4 )  ebc1=0
      IF( wbc == 4 )  wbc1=0
      IF( nbc == 4 )  nbc1=0
      IF( sbc == 4 )  sbc1=0

      IF (mp_opt > 0) THEN
        CALL acct_interrupt(mp_acct)
        CALL mpsendrecv2dew(w,nx,ny,nz,ebc1,wbc1,3,tem3)
        CALL mpsendrecv2dns(w,nx,ny,nz,nbc1,sbc1,3,tem3)
      END IF
      CALL acct_interrupt(mp_acct)
      CALL lbcw(nx,ny,nz,dtsml,w,wcont,wdteb,wdtwb,wdtnb,wdtsb,         &
                ebc1,wbc1,nbc1,sbc1,                                    &
                ebc_global,wbc_global,nbc_global,sbc_global)
      CALL acct_stop_inter

    ELSE ! Apply zero gradient condition

      ebc1 = 3
      wbc1 = 3
      sbc1 = 3
      nbc1 = 3
      IF( ebc == 0 )  ebc1 = 0
      IF( wbc == 0 )  wbc1 = 0
      IF( nbc == 0 )  nbc1 = 0
      IF( sbc == 0 )  sbc1 = 0

      IF (mp_opt > 0) THEN
        CALL acct_interrupt(mp_acct)
        CALL mpsendrecv2dew(w,nx,ny,nz,ebc1,wbc1,3,tem3)
        CALL mpsendrecv2dns(w,nx,ny,nz,nbc1,sbc1,3,tem3)
      END IF
      CALL acct_interrupt(mp_acct)
      CALL lbcw(nx,ny,nz,dtsml,w,wcont,wdteb,wdtwb,wdtnb,wdtsb,         &
                ebc1,wbc1,nbc1,sbc1,                                    &
                ebc_global,wbc_global,nbc_global,sbc_global)
      CALL acct_stop_inter

    END IF

  END IF
!
!-----------------------------------------------------------------------
!
!  Calculate wcont at time tfuture, including the boundary
!  points. Wcont at the lateral boundaries is calculated
!  from boundary u, v and w. Wcont at the top and bottom
!  depends on the boundary condition option.
!
!-----------------------------------------------------------------------
!
  CALL wcontra(nx,ny,nz,u,v,w,mapfct,j1,j2,j3,aj3z,                     &
               rhostr,rhostru,rhostrv,rhostrw,wcont,tem1,tem2)
!
!-----------------------------------------------------------------------
!
!  Set the top and bottom boundary conditions for w based on u, v and
!  wcont at the top and bottom boundaries.
!
!-----------------------------------------------------------------------
!
  CALL acct_interrupt(bc_acct)
  CALL vbcw(nx,ny,nz,w,wcont,tbc,bbc,u,v,                               &
            rhostr,rhostru,rhostrv,rhostrw,                             &
            j1,j2,j3)
  CALL acct_stop_inter

  IF( peqopt == 1) THEN

!-----------------------------------------------------------------------
!
!  Calculate the velocity divergence using newly updated velocity.
!
!-----------------------------------------------------------------------

    tema = dzinv*dtsml1
    DO k=2,nz-2
      DO j=1,ny-1
        DO i=1,nx-1
          div(i,j,k) = mapfct(i,j,7)                                    &
              *( u(i+1,j,k)*dxj3xm(i+1,j,k)- u(i  ,j,k)*dxj3xm(i  ,j,k) &
              +  v(i,j+1,k)*dyj3ym(i,j+1,k)- v(i,j  ,k)*dyj3ym(i,j  ,k)) &
                +( wcont(i,j,k+1)*aj3z(i,j,k+1)                         &
                 - wcont(i,j,k)*aj3z(i,j,k) ) * tema
        END DO
      END DO
    END DO

!-----------------------------------------------------------------------
!
!  Compute the divergence term and the base state pressure advection
!  in the pressure equation.
!
!  The pressure divergence term is: -rhostr*csndsq*div/j3
!  The base state pressure advection is -j3*w*d(pbar)/dzp = w*rhostr*g.
!
!  These two terms are saved in pdiv.
!
!-----------------------------------------------------------------------

    prgwg5=prgw * g05

    tema = dtsml1*prgw * g05
    DO k=2,nz-2
      DO j=1,ny-1
        DO i=1,nx-1
          pdiv(i,j,k)=                                                  &
                    tema*(w(i,j,k)+w(i,j,k+1))*rhostr(i,j,k)            &
                   -csndsq(i,j,k)*rstj3i(i,j,k)*div(i,j,k)
        END DO
      END DO
    END DO

  ELSE

!-----------------------------------------------------------------------
!
!  Calculate the density weighted mass divergence using newly
!  updated velocity. Then calculate the pressure divergence term:
!  -csndsq*div
!
!-----------------------------------------------------------------------
!
    DO k= 2,nz-2
      DO j= 1,ny-1
        DO i= 1,nx
          tem1(i,j,k)=u(i,j,k)*rhostru(i,j,k)*mapfct(i,j,5)
        END DO
      END DO
    END DO

    DO k= 2,nz-2
      DO j= 1,ny
        DO i= 1,nx-1
          tem2(i,j,k)=v(i,j,k)*rhostrv(i,j,k)*mapfct(i,j,6)
        END DO
      END DO
    END DO

    DO k= 2,nz-1
      DO j= 1,ny-1
        DO i= 1,nx-1
          tem3(i,j,k)=wcont(i,j,k)*rhostrw(i,j,k)
        END DO
      END DO
    END DO

    tema =  dxinv*dtsml1
    temb =  dyinv*dtsml1
    temc =  dzinv*dtsml1
    DO k=2,nz-2
      DO j=1,ny-1
        DO i=1,nx-1
          pdiv(i,j,k)= -csndsq(i,j,k)*( mapfct(i,j,7)                   &
                       *((tem1(i+1,j,k)-tem1(i,j,k))*tema               &
                       + (tem2(i,j+1,k)-tem2(i,j,k))*temb )             &
                       + (tem3(i,j,k+1)-tem3(i,j,k))*temc )
        END DO
      END DO
    END DO

  END IF
!
!-----------------------------------------------------------------------
!
!  Integrate forward the pressure equation by one small timestep,
!  using a backward (relative to the divergence term) integration
!  scheme.
!
!  And set lateral boundary conditions for the pressure equation.
!
!-----------------------------------------------------------------------
!
  IF( nestgrd == 1 .AND. mgrid /= 1 ) THEN

    DO k=2,nz-2
      DO j=2,ny-2
        DO i=2,nx-2
          pprt(i,j,k)=pprt(i,j,k)                                       &
                     +pforce(i,j,k)+pdiv(i,j,k)*j3inv(i,j,k)
        END DO
      END DO
    END DO

    DO k=2,nz-2
      DO j=1,ny-1
        pprt(   1,j,k)=pprt(   1,j,k)+dtsml1 * pdtwb(j,k)
        pprt(nx-1,j,k)=pprt(nx-1,j,k)+dtsml1 * pdteb(j,k)
      END DO
    END DO
    DO k=2,nz-2
      DO i=2,nx-2
        pprt(i,   1,k)=pprt(i,   1,k)+dtsml1 * pdtsb(i,k)
        pprt(i,ny-1,k)=pprt(i,ny-1,k)+dtsml1 * pdtnb(i,k)
      END DO
    END DO

!
!  Call the pprt bottom boundary condition subroutine to
!  compute the hydrostatic pprt at k=1.
!
    CALL acct_interrupt(bc_acct)
    CALL pprtbbc(nx,ny,nz,g05,buoy2nd,rhostr,pprt,ptprt,                &
                 pbari,ptbari,phydro,                                   &
                 tem1,tem2)          ! tem1 = new pprt at k=1.
    CALL acct_stop_inter

!
!-----------------------------------------------------------------------
!
!  And top and bottom boundary conditions for the pressure equation.
!
!-----------------------------------------------------------------------
!
! bc's for mp not implemented for nested grids

    CALL acct_interrupt(mp_acct)
    CALL bcp(nx,ny,nz, dtsml1, pprt, pdteb, pdtwb, pdtnb, pdtsb,        &
             tem1(1,1,1),0,0,0,0 ,tbc,bbc,                              &
             ebc_global,wbc_global,nbc_global,sbc_global)
    CALL acct_stop_inter

  ELSE ! Not nested grid

    DO k=2,nz-2
      DO j=1,ny-1
        DO i=1,nx-1
          pprt(i,j,k)=pprt(i,j,k)                                       &
                     +pforce(i,j,k)+pdiv(i,j,k)*j3inv(i,j,k)
        END DO
      END DO
    END DO

    IF(noisewrt == 1 .AND. MOD(ncalls(mptr),2) == 0 ) THEN
      IF(lbcopt == 2) THEN
        istart=MAX(2,(ngbrz+1))
        iend=MIN((nx-2),(nx-ngbrz-1))
        jstart=MAX(2,(ngbrz+1))
        jend=MIN((ny-2),(ny-ngbrz-1))
      ELSE
        istart=2
        iend=nx-2
        jstart=2
        jend=ny-2
      END IF
      sumdpdt=0.
      DO k=2,nz-2
        DO j=jstart,jend
          DO i=istart,iend
            sumdpdt=sumdpdt+abs(pforce(i,j,k)+pdiv(i,j,k)*j3inv(i,j,k))
          END DO
        END DO
      END DO
      sumdpdt=sumdpdt/                                                 &
              (float((iend-istart+1)*(jend-jstart+1)*(nz-3))*dtsml1)
      write(nchdpdt1(mptr),'(f10.2,f12.6)') (curtim+2*dtsml1),sumdpdt
    END IF
!
!  Call the pprt bottom boundary condition subroutine to
!  compute the hydrostatic pprt at k=1.
!
    CALL acct_interrupt(bc_acct)
    CALL pprtbbc(nx,ny,nz,g05,buoy2nd,rhostr,pprt,ptprt,                &
                 pbari,ptbari,phydro,                                   &
                 tem1,tem2)          ! tem1 = new pprt at k=1.
    CALL acct_stop_inter

!
!-----------------------------------------------------------------------
!
!  Apply the boundary conditions for the pressure equation.
!
!  For the open boundary case, the boundary value is solved from the
!  equation.
!
!-----------------------------------------------------------------------
!
    IF( lbcopt == 1 ) THEN  ! Internal boundary conditions
      ebc1=ebc
      wbc1=wbc
      sbc1=sbc
      nbc1=nbc

      IF( rbc_plbc == 1 ) then

        IF( ebc == 4 )  ebc1=0
        IF( wbc == 4 )  wbc1=0
        IF( nbc == 4 )  nbc1=0
        IF( sbc == 4 )  sbc1=0

      ELSE

        IF( ebc == 4 )  ebc1=3
        IF( wbc == 4 )  wbc1=3
        IF( nbc == 4 )  nbc1=3
        IF( sbc == 4 )  sbc1=3

      ENDIF

      IF (mp_opt > 0) THEN
        CALL acct_interrupt(mp_acct)
        CALL mpsendrecv2dew(pprt,nx,ny,nz,ebc1,wbc1,0,tem2)
        CALL mpsendrecv2dns(pprt,nx,ny,nz,nbc1,sbc1,0,tem2)
      END IF
      CALL acct_interrupt(mp_acct)
      CALL bcp(nx,ny,nz, dtsml1, pprt, pdteb, pdtwb, pdtnb, pdtsb,      &
               tem1(1,1,1), ebc1,wbc1,nbc1,sbc1,tbc,bbc,                &
               ebc_global,wbc_global,nbc_global,sbc_global)
      CALL acct_stop_inter

    ELSE  ! External boundary condition

      ebc1 = 3
      wbc1 = 3
      sbc1 = 3
      nbc1 = 3
      IF( ebc == 0 )  ebc1 = 0
      IF( wbc == 0 )  wbc1 = 0
      IF( nbc == 0 )  nbc1 = 0
      IF( sbc == 0 )  sbc1 = 0

      IF (mp_opt > 0) THEN
        CALL acct_interrupt(mp_acct)
        CALL mpsendrecv2dew(pprt,nx,ny,nz,ebc1,wbc1,0,tem2)
        CALL mpsendrecv2dns(pprt,nx,ny,nz,nbc1,sbc1,0,tem2)
      END IF
      CALL acct_interrupt(mp_acct)
      CALL bcp(nx,ny,nz, dtsml1, pprt, pdteb, pdtwb, pdtnb, pdtsb,      &
               tem1(1,1,1), 0,0,0,0,tbc,bbc,                            &
               ebc_global,wbc_global,nbc_global,sbc_global)
      CALL exbcp(nx,ny,nz, curtsml, pprt,                               &
                 exbcbuf(npr0exb),exbcbuf(nprdtexb))
      CALL acct_stop_inter

    END IF

  END IF

  RETURN
END SUBROUTINE solvwpex
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE SOLVWPIM                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE solvwpim(mptr, nx,ny,nz,exbcbufsz, dtsml1, curtsml,          &
           u,v,w,wcont,ptprt,pprt,phydro,                               &
           wdteb,wdtwb,wdtnb,wdtsb,pdteb,pdtwb,pdtnb,pdtsb,             &
           rhostr,ptbar,ptbari,pbari,csndsq,                            &
           wforce,wpgrad,pforce,                                        &
           x,y,z,zp, mapfct, j1,j2,j3,aj3x,aj3y,aj3z,j3inv,             &
           trigs1,trigs2,ifax1,ifax2,                                   &
           wsave1,wsave2,vwork1,vwork2,                                 &
           rhostru,rhostrv,rhostrw,                                     &
           exbcbuf,rhofct,                                              &
           fw,fp,wcontuv,tem1,tem2,tem3,tem4,j3irst,                    &
           csj32irst,rstpbi,rstwi,                                      &
           dxij3xm,dyij3ym,pbzi)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Perform the time integration of w and pressure equations using
!  implicit schem in vertical direction. The acoustically inactive
!  forcing terms in these equations have been evaluated prior to this
!  subroutine and are stored in wforce and pforce.
!
!-----------------------------------------------------------------------
!
!
!  AUTHOR: M. Xue & H. Jin
!  8/20/93.
!
!  9/6/94 (M.Xue)
!  A new version of vertically implicit small time step solver,
!  The perturbation pressure contribution to buoyancy and the
!  base state pressure advection terms are also treated implicitly
!  in the vertical direction inside the small time steps.
!
!  9/1/94 (Y. Lu)
!  Cleaned up documentation.
!
!  3/2/1995 (M. Xue and A. Shapiro)
!  Significant bug fix in loop 600. rhostru there was mistakenly
!  written as rhostrw.
!
!  10/31/95 (D. Weber)
!  Added linear hydrostatic upper radiation condition for w-p.
!  References are Klemp and Durran (JAS, 1983) and Chen (1991).
!  Includes the use of trigs1,trigs2,ifax1,ifax2.
!
!  1/25/96 (Donghai Wang & Yuhe Liu)
!  Added the map projection factor to ARPS governing equations.
!
!  4/27/1996 (M. Xue)
!  Corrected the formulations of the coefficients of the tridiagonal
!  equations. The oringal formulation was slightly inacurate.
!
!  5/6/96 (M. Xue)
!  Replaced csndsq in pressure buoyancy term by cpdcv*pbar/rhobar.
!
!  3/11/1997 (M. Xue and D. Weber)
!  Corrected a minor bug related to the radiation top boundary
!  condition.
!
!  07/22/97 (D. Weber)
!  Added even fft for linear hydrostatic w-p upper radiation condition.
!  References are Klemp and Durran (JAS, 1983) and Chen (1991).
!  Includes the use of wsave1,wsave2,vwork1,vwork2 (fftopt=2).
!
!  10/21/97 (Donghai Wang)
!  Using total density (rho) in the calculation of the pressure
!  gradient force terms.
!
!  10/21/97 (Donghai Wang)
!  Added the second order terms in the linerized buoyancy terms.
!
!  11/05/97 (D. Weber)
!  Added phydro array for use in the bottom boundary condition for
!  perturbation pressure (hydrostatic).
!
!  11/06/97 (D. Weber)
!  Added three additional levels to the mapfct array.  The three
!  levels (4,5,6) represent the inverse of the first three in order.
!  The inverse map factors are computed to improve efficiency.
!  Computed constants outside loops in selected loops.
!
!  9/18/98 (D. Weber)
!  Added arrays aj3x,y,z,ptbari, and pbari to improve the speed
!  of the code.  Also added tem5-12 for pre-computing groups of
!  commonly used constants.
!
!  3/18/99 (D. Weber)
!  Bug fix to second order buoyancy term (do loop 350).
!
!  07/23/2001 (K. Brewster)
!  Added mptr to argument list.
!  Added optional writing of mean abs(dp/dt) as a diagnostic noise
!  parameter.
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!    dtsml1   Local value of small time step size
!    curtsml  Local curttim at a small time step
!
!    u        x component of velocity at tfuture (m/s)
!    v        y component of velocity at tfuture (m/s)
!    w        Vertical component of Cartesian velocity at tfuture
!             (m/s)
!    wcont    Contravariant vertical velocity (m/s)
!    ptprt    Perturbation potential temperature at time tpresent (K)
!    pprt     Perturbation pressure at time tpresent (Pascal)
!
!    wdteb    Time tendency of the w field at the east boundary
!    wdtwb    Time tendency of the w field at the west boundary
!    wdtnb    Time tendency of the w field at the north boundary
!    wdtsb    Time tendency of the w field at the south boundary
!
!    pdteb    Time tendency of the pressure field at the east
!             boundary
!    pdtwb    Time tendency of the pressure field at the west
!             boundary
!    pdtnb    Time tendency of the pressure field at the north
!             boundary
!    pdtsb    Time tendency of the pressure field at the south
!             boundary
!
!    phydro   Big time step forcing term for use in computing the
!             hydrostatic pressure at k=1.
!
!    rhostr   Base state density rhobar times j3 (kg/m**3)
!    ptbar    Base state potential temperature (K)
!    ptbari   Inverse base state potential temperature (K)
!    pbari    Inverse base state pressure (Pascal)
!    csndsq   Sound wave speed squared.
!
!    wforce   Acoustically inactive forcing terms in w-eq.
!             (kg/(m*s)**2)
!    wpgrad   Pressure gradient term in w-eq. (kg/(m*s)**2)
!    pforce   Acoustically inactive terms in pressure eq. (Pascal/s)
!
!    x        x coordinate of grid points in physical/comp. space (m)
!    y        y coordinate of grid points in physical/comp. space (m)
!    z        z coordinate of grid points in computational space (m)
!    zp       Vertical coordinate of grid points in physical space(m)
!
!    mapfct   Map factors at scalar, u and v points
!
!    j1       Coordinate transformation Jacobian -d(zp)/dx
!    j2       Coordinate transformation Jacobian -d(zp)/dy
!    j3       Coordinate transformation Jacobian  d(zp)/dz
!    aj3x     Avgx of the coordinate transformation Jacobian  d(zp)/dz
!    aj3y     Avgy of the coordinate transformation Jacobian  d(zp)/dz
!    aj3z     Avgz of the coordinate transformation Jacobian  d(zp)/dz
!    trigs1   Array containing pre-computed trig function for
!             fftopt=1.
!    trigs2   Array containing pre-computed trig function for
!             fftopt=1.
!    ifax1    Array containing the factors of nx for fftopt=1.
!    ifax2    Array containing the factors of ny for fftopt=1.
!
!    vwork1   2-D work array for fftopt = 2.
!    vwork2   2-D work array for fftopt = 2.
!    wsave1   Work array for fftopt = 2.
!    wsave2   Work array for fftopt = 2.
!
!    rhostru  Average rhostr at u points (kg/m**3).
!    rhostrv  Average rhostr at v points (kg/m**3).
!    rhostrw  Average rhostr at w points (kg/m**3).
!
!  OUTPUT:
!
!    pprt     Perturbation pressure at tfuture updated in time
!             by dtsml1 (Pascal)
!    w        Vertical component of Cartesian velocity at tfuture
!             updated in time by dtsml1 (m/s)
!
!  WORK ARRAYS:
!
!    fw       Work array to carry force terms in w-eqation.
!    fp       Work array to carry force terms in p-eqation.
!    wcontuv  Work array to carry the contributions of u & v to wcont
!    tem1     Work array
!    tem2     Work array
!    tem3     Work array
!    tem4     Work array
!    j3irst   j3inv*rhostr
!    csj32irst csndsq*j3inv*j3inv*rhostr
!    rstpbi   rhostr*pbari
!    rstwi    1/rhostrw
!    dxij3xm  dxinv*aj3x*mapfct(i,j,5)
!    dyij3ym  dyinv*aj3y*mapfct(i,j,6)
!    pbzi     1/avgz(pbar)
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

  INTEGER :: mptr

  INTEGER :: nx, ny, nz        ! Number of grid points in 3 directions

  REAL :: dtsml1               ! Local small time step size (s)
  REAL :: curtsml              ! Local curttim at a small time step

  REAL :: u     (nx,ny,nz)     ! u-velocity at tfuture (m/s)
  REAL :: v     (nx,ny,nz)     ! v-velocity at tfuture (m/s)
  REAL :: w     (nx,ny,nz)     ! w-velocity at tfuture (m/s)
  REAL :: wcont (nx,ny,nz)     ! Contravariant vertical velocity (m/s)
  REAL :: ptprt (nx,ny,nz)     ! Perturbation potential temperature (K)
  REAL :: pprt  (nx,ny,nz)     ! Perturbation pressure at tfuture
                               ! (Pascal)

  REAL :: wdteb (ny,nz)        ! Time tendency of w at east boundary
  REAL :: wdtwb (ny,nz)        ! Time tendency of w at west boundary
  REAL :: wdtnb (nx,nz)        ! Time tendency of w at north boundary
  REAL :: wdtsb (nx,nz)        ! Time tendency of w at south boundary

  REAL :: pdteb (ny,nz)        ! Time tendency of pressure at east
                               ! boundary
  REAL :: pdtwb (ny,nz)        ! Time tendency of pressure at west
                               ! boundary
  REAL :: pdtnb (nx,nz)        ! Time tendency of pressure at north
                               ! boundary
  REAL :: pdtsb (nx,nz)        ! Time tendency of pressure at south
                               ! boundary

  REAL :: phydro(nx,ny)        ! Big time step forcing for computing
                               ! hydrostatic pprt at k=1.

  REAL :: rhostr(nx,ny,nz)     ! Base state density rhobar times j3.
  REAL :: ptbar (nx,ny,nz)     ! Base state potential temperature (K)
  REAL :: ptbari(nx,ny,nz)     ! Inverse base state pot. temperature (K)
  REAL :: pbari (nx,ny,nz)     ! Inverse base state pressure (Pascal).
  REAL :: csndsq(nx,ny,nz)     ! Sound wave speed squared.

  REAL :: wforce(nx,ny,nz)     ! Acoustically inactive forcing terms
                               ! in w-momentum equation (kg/(m*s)**2)
  REAL :: wpgrad(nx,ny,nz)     ! Pressure gradient term in w-eq.
                               ! (kg/(m*s)**2)
  REAL :: pforce(nx,ny,nz)     ! Acoustically inactive forcing terms
                               ! in the pressure equation (Pascal/s)

  REAL :: x     (nx)           ! x-coord. of the physical and compu-
                               ! tational grid. Defined at u-point.
  REAL :: y     (ny)           ! y-coord. of the physical and compu-
                               ! tational grid. Defined at v-point.
  REAL :: z     (nz)           ! z-coord. of the computational grid.
                               ! Defined at w-point on the staggered
                               ! grid.
  REAL :: zp    (nx,ny,nz)     ! Physical height coordinate defined at
                               ! w-point of the staggered grid.

  REAL :: mapfct(nx,ny,8)      ! Map factors at scalar, u and v points

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
  REAL :: aj3z  (nx,ny,nz)     ! Coordinate transformation Jacobian defined
                               ! as d( zp )/d( z ) AVERAGED IN THE Z-DIR.
  REAL :: j3inv (nx,ny,nz)     ! Coordinate transformation Jacobian
                               ! defined as d( zp )/d( z ).

  INTEGER :: exbcbufsz         ! EXBC buffer size
  REAL :: exbcbuf( exbcbufsz ) ! EXBC buffer array

  REAL :: rhofct(nx,ny,nz)     ! rho-factor: rhobar/rho

  REAL :: trigs1(3*(nx-1)/2+1) ! Array containing pre-computed trig
                               ! function for fftopt=1.
  REAL :: trigs2(3*(ny-1)/2+1) ! Array containing pre-computed trig
                               ! function for fftopt=1.
  INTEGER :: ifax1(13)         ! Array containing the factors of nx for
                               ! fftopt=1.
  INTEGER :: ifax2(13)         ! Array containing the factors of ny for
                               ! fftopt=1.

  REAL :: vwork1 (nx+1,ny+1)   ! 2-D work array for fftopt = 2.
  REAL :: vwork2 (ny,nx+1)     ! 2-D work array for fftopt = 2.
  REAL :: wsave1 (3*(ny-1)+15) ! Work array for fftopt = 2.
  REAL :: wsave2 (3*(nx-1)+15) ! Work array for fftopt = 2.

  REAL :: rhostru(nx,ny,nz)    ! Averaged rhostr at u points (kg/m**3).
  REAL :: rhostrv(nx,ny,nz)    ! Averaged rhostr at v points (kg/m**3).
  REAL :: rhostrw(nx,ny,nz)    ! Averaged rhostr at w points (kg/m**3).
!
!-----------------------------------------------------------------------
!
!  Temporary WORK ARRAYS:
!
!-----------------------------------------------------------------------
!
  REAL :: fp    (nx,ny,nz)     ! Pressure gradient term in w-eq.
  REAL :: fw    (nx,ny,nz)     ! Force in w equation.(kg/(m*s)**2)
  REAL :: wcontuv(nx,ny,nz)    ! Contributions of u & v to wcont
  REAL :: tem1  (nx,ny,nz)     ! Temporary work array
  REAL :: tem2  (nx,ny,nz)     ! Temporary work array
  REAL :: tem3  (nx,ny,nz)     ! Temporary work array
  REAL :: tem4  (nx,ny,nz)     ! Temporary work array
  REAL :: j3irst(nx,ny,nz)     ! j3inv*rhostr
  REAL :: csj32irst (nx,ny,nz) ! csndsq*j3inv*j3inv*rhostr
  REAL :: rstpbi(nx,ny,nz)     ! rhostr*pbari
  REAL :: rstwi (nx,ny,nz)     ! 1/rhostrw
  REAL :: dxij3xm(nx,ny,nz)    ! dxinv*aj3x*mapfct(i,j,5)
  REAL :: dyij3ym(nx,ny,nz)    ! dyinv*aj3y*mapfct(i,j,6)
  REAL :: pbzi  (nx,ny,nz)     ! 1/avgz(pbar)
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i, j, k
  INTEGER :: istart,iend,jstart,jend
  INTEGER :: ebc1,wbc1,nbc1,sbc1
  REAL :: g05,pk,nk,g05wp,g05pr

  INTEGER :: wpprt, prgw, itema
  REAL :: tema,temb,w1,w2,nrho
  INTEGER :: buoy2swt !Switch for 1st-order or 2nd-order in buoyancy
  REAL :: ptemk,ptemk1,pttemk,pttemk1

  INTEGER :: noisewrt
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
  INCLUDE 'exbc.inc'
  INCLUDE 'phycst.inc'      ! Physical constants
  INCLUDE 'mp.inc'            ! Message passing parameters.

  REAL :: sumdpdt

  CHARACTER (LEN=256) :: dpdtfn
  INTEGER :: ldpdtfn,istat

  INTEGER :: ncalls(mgrdmax), nchdpdt1(mgrdmax)
  SAVE ncalls,nchdpdt1
  DATA ncalls /mgrdmax*0/

!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  IF( bsnesq == 1 ) THEN
    wpprt = 0   ! Switch for pprt/(rhobar*csndsq) term in w-eq.
    prgw  = 0   ! Switch for rhobar*g*w term in w-eq.
  ELSE
    wpprt = 1   ! Switch for pprt/(rhobar*csndsq) term in w-eq.
    prgw  = 1   ! Switch for rhobar*g*w term in w-eq.
  END IF

  g05 = g*0.5
  g05wp = g05*wpprt
  g05pr = g05*prgw

  noisewrt = 0

  ncalls(mptr) = ncalls(mptr) + 1

  IF (noisewrt == 1 .AND. ncalls(mptr) == 1) THEN
      dpdtfn  = runname(1:lfnkey)//'.dpdt'
      ldpdtfn = 5 + lfnkey

      IF(nestgrd == 1) THEN
        WRITE(dpdtfn((ldpdtfn+1):(ldpdtfn+4)),'(a,i2.2)')'.g',mptr
        ldpdtfn =  ldpdtfn + 4
      END IF

      WRITE(6,'(1x,a,a,a/,1x,a)')                                       &
          'Check to see if file ',dpdtfn(1:ldpdtfn),' already exists.', &
          'If so, append a version number to the filename.'

      CALL fnversn( dpdtfn, ldpdtfn )

      CALL getunit ( nchdpdt1(mptr) )

      OPEN(nchdpdt1(mptr),FORM='formatted',STATUS='new',                &
                FILE=dpdtfn(1:ldpdtfn),IOSTAT=istat)

      IF( istat /= 0) THEN

        WRITE(6,'(/a,i2,/a/)')                                          &
            ' Error occured when opening file '//dpdtfn(1:ldpdtfn)//    &
            ' using FORTRAN unit ',nchdpdt1(mptr),                      &
            ' Program stopped in SOLVWPIM.'
        CALL arpsstop('arpsstop called from SOLVWPIM error opening '//  &
             'file',1)

      END IF

      WRITE(nchdpdt1(mptr),'(a)') runname
      WRITE(nchdpdt1(mptr),'(t2,a,t15,a)') 'Time (s)','MnAbsdpdt'

  END IF

!-----------------------------------------------------------------------
!
!  Calculate the horizontal velocity divergence using newly updated
!  u and v velocity plus half vertical divergence from wcont, and
!  store the result in tem3.
!
!  Namely, tem3 = difx(u*avgx(j3))+dify(v*avgy(j3))+
!                 (1-tacoef)*difz(wcont*avgz(j3))
!
!    note tem9=dtsml*aj3x*mapfct(5)*dxinv
!    note tem10=dtsml*aj3y*mapfct(6)*dyinv
!-----------------------------------------------------------------------

  tema = dtsml1*(1.0-tacoef)*dzinv
  DO k=2,nz-2
    DO j=1,ny-1
      DO i=1,nx-1
        tem3(i,j,k) = mapfct(i,j,7)                                     &
                    * ( (u(i+1,j,k)*dxij3xm(i+1,j,k)                    &
                        -u(i,j,k)*dxij3xm(i,j,k))                       &
                      + (v(i,j+1,k)*dyij3ym(i,j+1,k)                    &
                        -v(i,j,k)*dyij3ym(i,j,k)) )                     &
            +(wcont(i,j,k+1)*aj3z(i,j,k+1)-wcont(i,j,k)*aj3z(i,j,k))    &
            *tema
      END DO
    END DO
  END DO

!
!-----------------------------------------------------------------------
!
!  Compute the forcing terms in pressure equation
!
!  fp = dtsml/j3 *(pforce-rhobar*csndsq*tem3+(1-tacoef)*rhostr*g*w)
!
!  note pforce includes dtsml1 and j3inv...
!
!  rhostr*g*w is the base state pressure advection term.
!
!-----------------------------------------------------------------------
!

  tema = dtsml1*g05pr*(1.0-tacoef)
  DO k=2,nz-2
    DO j=1,ny-1
      DO i=1,nx-1
        fp(i,j,k)=pforce(i,j,k)                                         &
            -csj32irst(i,j,k)*tem3(i,j,k)                               &
            +tema*(w(i,j,k)+w(i,j,k+1))*j3irst(i,j,k)
      END DO
    END DO
  END DO
!
!-----------------------------------------------------------------------
!
!  Compute the first two terms in the contravariant vertical
!  velocity formulation:
!
!  wcontuv = (avgx(avgz(ustr)*j1)+avgy(avgz(vstr)*j2))/rhostrw
!
!-----------------------------------------------------------------------
!
  IF( ternopt == 0 ) THEN

    DO k=2,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          wcontuv(i,j,k)=0.0
        END DO
      END DO
    END DO

  ELSE

    DO k= 1,nz-1
      DO j= 1,ny-1
        DO i= 1,nx
          tem1(i,j,k)=u(i,j,k)*rhostru(i,j,k)
        END DO
      END DO
    END DO

    DO k= 1,nz-1
      DO j= 1,ny
        DO i= 1,nx-1
          tem2(i,j,k)=v(i,j,k)*rhostrv(i,j,k)
        END DO
      END DO
    END DO

    DO k= 2,nz-1
      DO j= 1,ny-1
        DO i= 1,nx-1
          wcontuv(i,j,k)=((tem1(i  ,j,k)+tem1(i  ,j,k-1))*j1(i  ,j,k)   &
                         +(tem1(i+1,j,k)+tem1(i+1,j,k-1))*j1(i+1,j,k)   &
                         +(tem2(i  ,j,k)+tem2(i  ,j,k-1))*j2(i  ,j,k)   &
                         +(tem2(i,j+1,k)+tem2(i,j+1,k-1))*j2(i,j+1,k))  &
                         *mapfct(i,j,8) * rstwi(i,j,k)
        END DO
      END DO
    END DO

    DO j=1,ny-1
      DO i=1,nx-1
        wcontuv(i,j,nz-1)=0.0
        wcontuv(i,j,2) = ((u(i  ,j,2)+u(i  ,j,1))*j1(i,j,2)             &
                         +(u(i+1,j,2)+u(i+1,j,1))*j1(i+1,j,2)           &
                         +(v(i,j  ,2)+v(i,j  ,1))*j2(i,j,2)             &
                         +(v(i,j+1,2)+v(i,j+1,1))*j2(i,j+1,2))          &
                         *mapfct(i,j,8)
      END DO
    END DO

  END IF
!
!-----------------------------------------------------------------------
!
!  Average rhofct at w points, stored in tem4
!
!-----------------------------------------------------------------------
!
  CALL avgsw(rhofct,nx,ny,nz, 1,nx-1, 1,ny-1,tem4)
!
!-----------------------------------------------------------------------
!
!  Compute the right-hand-side forcing term in tridiagonal linear
!  equation. Array w is used to store this forcing term.
!
!  First, we add the contribution of pressure perturbation to the
!  buoyancy to wforce and wpgrad.
!
!  This term has a form of - rhostr*g*pprt/(cpdcv*pbar),
!
!-----------------------------------------------------------------------
!

  tema = g05wp/cpdcv

  IF ( buoy2nd == 0) THEN  !1st-order
    buoy2swt = 0             !Switch for 1st-order or 2nd-order
  ELSE                       !2nd-order
    buoy2swt = 1
  END IF


  temb = 0.5*buoy2swt*(1.0-cpdcv)/cpdcv
  DO k=3,nz-2
    DO j=1,ny-1
      DO i=1,nx-1
        ptemk = pprt(i,j,k  )*pbari(i,j,k  )    ! new code...
        ptemk1= pprt(i,j,k-1)*pbari(i,j,k-1)
        fw(i,j,k) = wforce(i,j,k)-tem4(i,j,k)*(wpgrad(i,j,k)            &
                   +rhostrw(i,j,k)*(ptemk+ptemk1+                       &
                     temb*(ptemk*ptemk+ptemk1*ptemk1))                  &
                   *tema)
      END DO
    END DO
  END DO
!
!-----------------------------------------------------------------------
!
!  When potential temperature equation is solved within small time
!  steps, the contribution from ptprt to buoyancy term is calculated
!  here.
!
!-----------------------------------------------------------------------
!
  IF( ptsmlstp == 1 ) THEN
    tema = 1.0/(2.0*cpdcv)
    DO k=3,nz-2
      DO j=1,ny-1
        DO i=1,nx-1
          ptemk = pprt(i,j,k  )*pbari(i,j,k  )
          ptemk1= pprt(i,j,k-1)*pbari(i,j,k-1)
          pttemk = ptprt(i,j,k  )*ptbari(i,j,k  )
          pttemk1= ptprt(i,j,k-1)*ptbari(i,j,k-1)
          fw(i,j,k)=fw(i,j,k)+                                          &
                    rhostrw(i,j,k)*(pttemk+pttemk1                      &
                    -buoy2swt*(pttemk*pttemk+pttemk1*pttemk1            &
                    +tema*(ptemk*pttemk + ptemk1*pttemk1)))             &
                    *g05*tem4(i,j,k)
        END DO
      END DO
    END DO

  END IF
!
!-----------------------------------------------------------------------
!
!  tem3=fp-dtsml*rhostr*csndsq*tacoef* d(wcontuv)/dz /(j3*j3)
!
!-----------------------------------------------------------------------
!
  tema = tacoef*dtsml1*dzinv
  DO k=2,nz-2
    DO j=1,ny-1
      DO i=1,nx-1
        tem3(i,j,k)=fp(i,j,k)-tema*csj32irst(i,j,k)                     &
                              *(wcontuv(i,j,k+1)-wcontuv(i,j,k))
      END DO
    END DO
  END DO
!
!-----------------------------------------------------------------------
!
!  fw=fw-tacoef*d(tem3)/dz-g*tacoef*avgz(rhostr*tem3/(cpdcv*pbar)))
!
!-----------------------------------------------------------------------
!
  tema = g05wp/cpdcv

  DO k=3,nz-2
    DO j=1,ny-1
      DO i=1,nx-1
        fw(i,j,k)=fw(i,j,k)-tacoef*tem4(i,j,k)*(                        &
                  (tem3(i,j,k)-tem3(i,j,k-1))*dzinv +                   &
                  (rstpbi(i,j,k  )*tem3(i,j,k  )                        &
                  +rstpbi(i,j,k-1)*tem3(i,j,k-1))*tema)
      END DO
    END DO
  END DO

!
!-----------------------------------------------------------------------
!
!  fw=fw*dtsml/rhostr + w
!
!-----------------------------------------------------------------------
!
  DO k=3,nz-2
    DO j=1,ny-1
      DO i=1,nx-1
        fw(i,j,k)=dtsml1*(fw(i,j,k)*rstwi(i,j,k)) + w(i,j,k)
      END DO
    END DO
  END DO

!
!-----------------------------------------------------------------------
!
!  Prepare the left-hand-side coefficents for tridiagonal
!  equation set:
!
!     A(k)*w(k-1)+B(k)*w(k)+C(k)*w(k+1)=D(k)   for k=3,nz-2
!
!  w(2) and w(nz-1) are used as the boundary conditions.
!
!  Due to the lack of work arrays, we are storing the coefficients
!  A in rhostru, B in rhostrv, and C in rhostrw. D is stored in fw.
!
!  rhostru, rhostrv and rhostrw are re-calculated from rhostr after
!  they have been used.
!
!-----------------------------------------------------------------------
!
  DO k=2,nz-2
    DO j=1,ny-1
      DO i=1,nx-1
        tem1(i,j,k)= g05pr*j3irst(i,j,k)
        tem2(i,j,k)= dzinv*csj32irst(i,j,k)
      END DO
    END DO
  END DO

  tema = (dtsml1*tacoef)**2 * wpprt*g/cpdcv
  temb = (dtsml1*tacoef)**2 * dzinv

  DO k=3,nz-2
    DO j=1,ny-1
      DO i=1,nx-1
        pk = tem4(i,j,k)*tema*pbzi(i,j,k)
        nk = tem4(i,j,k)*temb*rstwi(i,j,k)

        rhostru(i,j,k)= (-nk+pk)*(tem1(i,j,k-1)+tem2(i,j,k-1))
        rhostrw(i,j,k)= ( nk+pk)*(tem1(i,j,k  )-tem2(i,j,k  ))
        rhostrv(i,j,k)= 1                                               &
            +nk*(tem1(i,j,k)+tem2(i,j,k)-tem1(i,j,k-1)+tem2(i,j,k-1))   &
            +pk*(tem1(i,j,k)+tem2(i,j,k)+tem1(i,j,k-1)-tem2(i,j,k-1))
      END DO
    END DO
  END DO

!
!-----------------------------------------------------------------------
!
!  Reset fw on the boundaries using the top and bottom boundary
!  conditions of w.
!
!  w=0.0 at k=nz-1. At the lower boundary, wcont=0.0, therefore
!  w=-wcontuv, where wcontuv was calculated in loop 240.
!
!-----------------------------------------------------------------------
!
  DO j=1,ny-1
    DO i=1,nx-1
      fw(i,j,3)   =fw(i,j,3)+wcontuv(i,j,2)*rhostru(i,j,3)
      fw(i,j,nz-2)=fw(i,j,nz-2)+0.0*rhostrw(i,j,nz-2)
    END DO
  END DO


!
!-----------------------------------------------------------------------
!
!  Call the tridiagonal solver for either a rigid (tbc=1) or open upper
!  boundary condition for w (tbc=4).
!
!  For tbc=4, we have a choice of fft's for the upper boundary:
!  fftopt =1, periodic fft for linearized hydrostatic radiation condition.
!  fftopt =2, even fft for linearized hydrostatic radiation condition.
!
!  The references are Klemp and Durran Jas (1983) and Chen (MWR) 1991.
!
!   w1 is the coef for w(nz-1) and w2 is the coef. for w(nz-2) in the
!   pressure equation.  The pressure equation is solved at p(nz-2)
!   for w(nz-1).   tem2(i,j,nz-2) is the summation of the known terms
!   in the pressure equation.  The only unknowns are p(nz-2), w(nz-1)
!   and w(nz-1) at the new time level.  In the tridiagonal solver
!   the elimination step is performed and w(nz-2) is given in termsx
!   of w(nz-1) and known terms.  THe next step is to use the rad.
!   condition (in fourier space):
!
!
!        p = N * rhobar * w(nz-1) / abs(kx+ky)
!
!   The end result is a relation for w(nz-1) given known quantities.
!   Trigs1 and trigs2 are the predetermined trig functions used in the
!   fft program in tridiag.
!
!-----------------------------------------------------------------------

  IF(tbc == 4)THEN ! apply linear radiation condition to w at nz-1.

    tema = rhostr(nx/2,ny/2,nz-2)*csndsq(nx/2,ny/2,nz-2)*               &
                                      j3inv(nx/2,ny/2,nz-2)
    temb = dtsml1*j3inv(nx/2,ny/2,nz-2)
    w2=temb*tacoef*(tema*dzinv+g05*prgw*rhostr(nx/2,ny/2,nz-2))
    w1=temb*tacoef*(-tema*dzinv+g05*prgw*rhostr(nx/2,ny/2,nz-2))
    DO j=1,ny-1
      DO i=1,nx-1
        tem2(i,j,2)=pprt(i,j,nz-2)+tem3(i,j,nz-2)
      END DO
    END DO

    nrho=SQRT(g/ptbar(nx/2,ny/2,nz-2)*(ptbar(nx/2,ny/2,nz-1)            &
              -ptbar(nx/2,ny/2,nz-3))*dzinv*0.5)                        &
        *rhostr(nx/2,ny/2,nz-2)*j3inv(nx/2,ny/2,nz-2)

  END IF

!
!-----------------------------------------------------------------------
!
!  NOTE:  tem2(1,1,4) is passed into subroutine tridiag and becomes
!         a 2-D array of size (nx+1,ny+1).
!
!-----------------------------------------------------------------------
!
  CALL tridiag(nx,ny,nz,rhostru,rhostrv,rhostrw,fw,tem1,                &
      tem2(1,1,1),tem2(1,1,2),w1,w2,nrho,tem2(1,1,3),trigs1,trigs2,     &
      ifax1,ifax2,wsave1,wsave2,vwork1,vwork2,tem2(1,1,4))

!
!-----------------------------------------------------------------------
!
!  Restore rhostru, rhostrv and rhostrw after they have been used
!  are work arrays.
!
!-----------------------------------------------------------------------
!
  CALL rhouvw(nx,ny,nz,rhostr,rhostru,rhostrv,rhostrw)
!
!-----------------------------------------------------------------------
!
!  On exit of tridiag, the interior solution of w is stored in fw.
!
!-----------------------------------------------------------------------

  IF(tbc == 4)THEN   ! set itema = nz-1
    itema = nz-1
  ELSE               ! set itema = nz-2
    itema = nz-2
  END IF

  IF( nestgrd /= 1 .OR. mgrid == 1 ) THEN
                                     ! For non-nesting case or the
                                     ! base-grid of the nested case

    DO k=3,itema
      DO j=1,ny-1
        DO i=1,nx-1
          w(i,j,k) = fw(i,j,k)
        END DO
      END DO
    END DO

    IF ( lbcopt == 1 ) THEN ! Internal boundary condition
      ebc1=ebc
      wbc1=wbc
      sbc1=sbc
      nbc1=nbc
      IF( ebc == 4 )  ebc1=0
      IF( wbc == 4 )  wbc1=0
      IF( nbc == 4 )  nbc1=0
      IF( sbc == 4 )  sbc1=0

      IF (mp_opt > 0) THEN
        CALL acct_interrupt(mp_acct)
        CALL mpsendrecv2dew(w,nx,ny,nz,ebc1,wbc1,3,tem3)
        CALL mpsendrecv2dns(w,nx,ny,nz,nbc1,sbc1,3,tem3)
      END IF
      CALL acct_interrupt(mp_acct)
      CALL lbcw(nx,ny,nz,dtsml1,w,wcont,wdteb,wdtwb,wdtnb,wdtsb,        &
                ebc1,wbc1,nbc1,sbc1,                                    &
                ebc_global,wbc_global,nbc_global,sbc_global)
      CALL acct_stop_inter

    ELSE  ! External boundary condition

      ebc1 = 3
      wbc1 = 3
      sbc1 = 3
      nbc1 = 3
      IF( ebc == 0 )  ebc1 = 0
      IF( wbc == 0 )  wbc1 = 0
      IF( nbc == 0 )  nbc1 = 0
      IF( sbc == 0 )  sbc1 = 0

      IF (mp_opt > 0) THEN
        CALL acct_interrupt(mp_acct)
        CALL mpsendrecv2dew(w,nx,ny,nz,ebc1,wbc1,3,tem3)
        CALL mpsendrecv2dns(w,nx,ny,nz,nbc1,sbc1,3,tem3)
      END IF
      CALL acct_interrupt(mp_acct)
      CALL lbcw(nx,ny,nz,dtsml1,w,wcont,wdteb,wdtwb,wdtnb,wdtsb,        &
                ebc1,wbc1,nbc1,sbc1,                                    &
                ebc_global,wbc_global,nbc_global,sbc_global)
      CALL acct_stop_inter

    END IF

  ELSE    ! For nested interior grid

    DO k=3,itema
      DO j=2,ny-2
        DO i=2,nx-2
          w(i,j,k) = fw(i,j,k)
        END DO
      END DO
    END DO

    DO k=3,nz-2
      DO j=1,ny-1
        w(   1,j,k)=w(   1,j,k)+dtsml1 * wdtwb(j,k)
        w(nx-1,j,k)=w(nx-1,j,k)+dtsml1 * wdteb(j,k)
      END DO
    END DO
    DO k=3,nz-2
      DO i=2,nx-2
        w(i,   1,k)=w(i,   1,k)+dtsml1 * wdtsb(i,k)
        w(i,ny-1,k)=w(i,ny-1,k)+dtsml1 * wdtnb(i,k)
      END DO
    END DO

  END IF
!
!-----------------------------------------------------------------------
!
!  Calculate wcont at time tfuture, including the boundary
!  points. Wcont at the lateral boundaries is calculated
!  from boundary u, v and w. Wcont at the top and bottom
!  depends on the boundary condition option.
!
!-----------------------------------------------------------------------
!
  CALL wcontra(nx,ny,nz,u,v,w,mapfct,j1,j2,j3,aj3z,                     &
               rhostr,rhostru,rhostrv,rhostrw,wcont,tem1,tem2)
!
!-----------------------------------------------------------------------
!
!  Set the top and bottom boundary conditions for w based on u, v and
!  wcont at the top and bottom boundaries.
!
!-----------------------------------------------------------------------
!
  CALL acct_interrupt(bc_acct)
  CALL vbcw(nx,ny,nz,w,wcont,tbc,bbc,u,v,                               &
            rhostr,rhostru,rhostrv,rhostrw,                             &
            j1,j2,j3)
  CALL acct_stop_inter
!
!-----------------------------------------------------------------------
!
!  Calculate the new pressure
!
!  pprt(new) = pprt(old)+fp+dtsml*tacoef/j3*
!             (rhostr*g*avg(w)-rhostr/j3*csndsq*difz(j3*wcont))
!
!-----------------------------------------------------------------------

  tema = tacoef*dtsml1*g05pr
  temb = tacoef*dtsml1*dzinv

  IF( nestgrd == 1 .AND. mgrid /= 1 ) THEN
                                        ! For nested interior grid

    DO k=2,nz-2
      DO j=2,ny-2
        DO i=2,nx-2
          pprt(i,j,k)=pprt(i,j,k) + fp(i,j,k)                           &
              +tema*j3irst(i,j,k)*(w(i,j,k+1)+w(i,j,k))                 &
              -temb*csj32irst(i,j,k)*                                   &
              (aj3z(i,j,k+1)*wcont(i,j,k+1)-aj3z(i,j,k)*wcont(i,j,k))
        END DO
      END DO
    END DO

    DO k=2,nz-2
      DO j=1,ny-1
        pprt(   1,j,k)=pprt(   1,j,k)+dtsml1 * pdtwb(j,k)
        pprt(nx-1,j,k)=pprt(nx-1,j,k)+dtsml1 * pdteb(j,k)
      END DO
    END DO
    DO k=2,nz-2
      DO i=2,nx-2
        pprt(i,   1,k)=pprt(i,   1,k)+dtsml1 * pdtsb(i,k)
        pprt(i,ny-1,k)=pprt(i,ny-1,k)+dtsml1 * pdtnb(i,k)
      END DO
    END DO

!
!  Call the pprt bottom boundary condition subroutine to
!  compute the hydrostatic pprt at k=1.
!
    CALL acct_interrupt(bc_acct)
    CALL pprtbbc(nx,ny,nz,g05,buoy2swt,rhostr,pprt,ptprt,               &
                 pbari,ptbari,phydro,                                   &
                 tem1,tem2)         ! tem1 = new pprt at k=1.
    CALL acct_stop_inter

! bc's for mp not implemented for nested grids
    CALL acct_interrupt(mp_acct)
    CALL bcp(nx,ny,nz, dtsml1, pprt, pdteb, pdtwb, pdtnb, pdtsb,        &
             tem1(1,1,1), 0,0,0,0,tbc,bbc,                              &
             ebc_global,wbc_global,nbc_global,sbc_global)
    CALL acct_stop_inter

  ELSE                              ! For non-nesting case or the
                                    ! base-grid of the nested case
    DO k=2,nz-2
      DO j=1,ny-1
        DO i=1,nx-1
          pprt(i,j,k)=pprt(i,j,k) + fp(i,j,k)                           &
              +tema*j3irst(i,j,k)*(w(i,j,k+1)+w(i,j,k))                 &
              -temb*csj32irst(i,j,k)*                                   &
              (aj3z(i,j,k+1)*wcont(i,j,k+1)-aj3z(i,j,k)*wcont(i,j,k))
        END DO
      END DO
    END DO

    IF (noisewrt == 1 .AND. MOD(ncalls(mptr),2) == 0 ) THEN
      IF(lbcopt == 2) THEN
        istart=MAX(2,(ngbrz+1))
        iend=MIN((nx-2),(nx-ngbrz-1))
        jstart=MAX(2,(ngbrz+1))
        jend=MIN((ny-2),(ny-ngbrz-1))
      ELSE
        istart=2
        iend=nx-2
        jstart=2
        jend=ny-2
      END IF
      sumdpdt=0.
      DO k=2,nz-2
        DO j=jstart,jend
          DO i=istart,iend
            sumdpdt=sumdpdt+abs(fp(i,j,k)                                 &
                +tema*j3irst(i,j,k)*(w(i,j,k+1)+w(i,j,k))                 &
                -temb*csj32irst(i,j,k)*                                   &
                (aj3z(i,j,k+1)*wcont(i,j,k+1)-aj3z(i,j,k)*wcont(i,j,k)))
          END DO
        END DO
      END DO
      sumdpdt=sumdpdt/                                                    &
              (float((iend-istart+1)*(jend-jstart+1)*(nz-3))*dtsml1)
      write(nchdpdt1(mptr),'(f10.2,f12.6)') (curtim+2*dtsml1),sumdpdt
    END IF
!
!  Call the pprt bottom boundary condition subroutine to
!  compute the hydrostatic pprt at k=1.
!
    CALL acct_interrupt(bc_acct)
    CALL pprtbbc(nx,ny,nz,g05,buoy2swt,rhostr,pprt,ptprt,               &
                 pbari,ptbari,phydro,                                   &
                 tem1,tem2)           ! tem1 = new pprt at k=1.
    CALL acct_stop_inter

!
!-----------------------------------------------------------------------
!
!  Apply the boundary conditions for the pressure equation.
!
!  For the open boundary case, the boundary value of pprt is
!  predicted by the pressure equation.
!
!-----------------------------------------------------------------------
!
    IF( lbcopt == 1 ) THEN  ! Internal boundary conditions
      ebc1=ebc
      wbc1=wbc
      sbc1=sbc
      nbc1=nbc

      IF( rbc_plbc == 1 ) THEN

        IF( ebc == 4 )  ebc1=0
        IF( wbc == 4 )  wbc1=0
        IF( nbc == 4 )  nbc1=0
        IF( sbc == 4 )  sbc1=0

      ELSE

        IF( ebc == 4 )  ebc1=3
        IF( wbc == 4 )  wbc1=3
        IF( nbc == 4 )  nbc1=3
        IF( sbc == 4 )  sbc1=3

      ENDIF

      IF (mp_opt > 0) THEN
        CALL acct_interrupt(mp_acct)
        CALL mpsendrecv2dew(pprt,nx,ny,nz,ebc1,wbc1,0,tem2)
        CALL mpsendrecv2dns(pprt,nx,ny,nz,nbc1,sbc1,0,tem2)
      END IF
      CALL acct_interrupt(mp_acct)
      CALL bcp(nx,ny,nz, dtsml1, pprt, pdteb, pdtwb, pdtnb, pdtsb,      &
               tem1(1,1,1), ebc1,wbc1,nbc1,sbc1,tbc,bbc,                &
               ebc_global,wbc_global,nbc_global,sbc_global)
      CALL acct_stop_inter

    ELSE  ! External boundary condition

      IF (mp_opt > 0) THEN
        CALL acct_interrupt(mp_acct)
        CALL mpsendrecv2dew(pprt,nx,ny,nz,0,0,0,tem2)
        CALL mpsendrecv2dns(pprt,nx,ny,nz,0,0,0,tem2)
      END IF
      CALL acct_interrupt(mp_acct)
      CALL bcp(nx,ny,nz, dtsml1, pprt, pdteb, pdtwb, pdtnb, pdtsb,      &
               tem1(1,1,1), 0,0,0,0,tbc,bbc,                            &
               ebc_global,wbc_global,nbc_global,sbc_global)
      CALL exbcp(nx,ny,nz, curtsml, pprt,                               &
                 exbcbuf(npr0exb),exbcbuf(nprdtexb))
      CALL acct_stop_inter

    END IF

  END IF

  RETURN
END SUBROUTINE solvwpim

!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE SOLVWPIM1                  ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE solvwpim1(nx,ny,nz,exbcbufsz, dtsml1, curtsml,               &
           u,v,w,wcont,ptprt,pprt,phydro,                               &
           wdteb,wdtwb,wdtnb,wdtsb,pdteb,pdtwb,pdtnb,pdtsb,             &
           rhostr,ptbar,ptbari,pbari,csndsq,                            &
           wforce,wpgrad,pforce,                                        &
           x,y,z,zp, mapfct, j1,j2,j3,aj3z,j3inv,                       &
           trigs1,trigs2,ifax1,ifax2,                                   &
           wsave1,wsave2,vwork1,vwork2,                                 &
           rhostru,rhostrv,rhostrw,                                     &
           exbcbuf,rhofct,                                              &
           fw,fp,tem1,tem2,tem3,tem4,tem5,                              &
           rstpbi,rstptbi,rstwi,pbzi,rstj3i)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Perform the time integration of w and pressure equations using
!  implicit scheme in vertical direction. The acoustically inactive
!  forcing terms in these equations have been evaluated prior to this
!  subroutine and are stored in wforce and pforce.
!
!  This routine is for peqopt=2.
!
!-----------------------------------------------------------------------
!
!
!  AUTHOR: M. Xue
!  4/26/1996.
!  Written for peqopt=2 based on SOLVWPIM.
!
!  MODIFICATION HISTORY:
!
!  07/22/97 (D. Weber)
!  Added wsave1,wsave2,vwork1,vwork2 for use in the even fft version
!  of the upper w-p radiation condition (fftopt=2).
!
!  10/21/97 (Donghai Wang)
!  Using total density (rho) in the calculation of the pressure
!  gradient force terms.
!
!  10/21/97 (Donghai Wang)
!  Added the second order terms in the linerized buoyancy terms.
!
!  11/05/97 (D. Weber)
!  Added phydro array for use in the bottom boundary condition for
!  perturbation pressure (hydrostatic).
!
!  11/06/97 (D. Weber)
!  Added three additional levels to the mapfct array.  The three
!  levels (4,5,6) represent the inverse of the first three in order.
!  The inverse map factors are computed to improve efficiency.
!  Removed multiplication of constants from loops.
!
!  9/18/98 (D. Weber)
!  Added arrays aj3z ptbari,pbari,rstpbi,rstptbi,rstwi,pbzi,
!  rstj3i to improve code efficiency.
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!    dtsml1   Local value of small time step size
!    curtsml  Local curttim at a small time step
!
!    u        x component of velocity at tfuture (m/s)
!    v        y component of velocity at tfuture (m/s)
!    w        Vertical component of Cartesian velocity at tfuture
!             (m/s)
!    wcont    Contravariant vertical velocity (m/s)
!    ptprt    Perturbation potential temperature at time tpresent (K)
!    pprt     Perturbation pressure at time tpresent (Pascal)
!
!    wdteb    Time tendency of the w field at the east boundary
!    wdtwb    Time tendency of the w field at the west boundary
!    wdtnb    Time tendency of the w field at the north boundary
!    wdtsb    Time tendency of the w field at the south boundary
!
!    pdteb    Time tendency of the pressure field at the east
!             boundary
!    pdtwb    Time tendency of the pressure field at the west
!             boundary
!    pdtnb    Time tendency of the pressure field at the north
!             boundary
!    pdtsb    Time tendency of the pressure field at the south
!             boundary
!
!    phydro   Big time step forcing term for use in computing the
!             hydrostatic pressure at k=1.
!
!    rhostr   Base state density rhobar times j3 (kg/m**3)
!    ptbar    Base state potential temperature (K)
!    ptbari   Inverse base state potential temperature (K)
!    pbari    Inverse base state pressure (Pascal)
!    csndsq   Sound wave speed squared.
!
!    wforce   Acoustically inactive forcing terms in w-eq.
!             (kg/(m*s)**2)
!    wpgrad   Pressure gradient term in w-eq. (kg/(m*s)**2)
!    pforce   Acoustically inactive terms in pressure eq. (Pascal/s)
!
!    x        x coordinate of grid points in physical/comp. space (m)
!    y        y coordinate of grid points in physical/comp. space (m)
!    z        z coordinate of grid points in computational space (m)
!    zp       Vertical coordinate of grid points in physical space(m)
!
!    mapfct   Map factors at scalar, u and v points
!
!    j1       Coordinate transformation Jacobian -d(zp)/dx
!    j2       Coordinate transformation Jacobian -d(zp)/dy
!    j3       Coordinate transformation Jacobian  d(zp)/dz
!    aj3z     Avgz of the coordinate transformation Jacobian  d(zp)/dz
!    trigs1   Array containing pre-computed trig function for
!             fftopt=1.
!    trigs2   Array containing pre-computed trig function for
!             fftopt=1.
!    ifax1    Array containing the factors of nx for fftopt=1.
!    ifax2    Array containing the factors of ny for fftopt=1.
!
!    vwork1   2-D work array for fftopt=2.
!    vwork2   2-D work array for fftopt=2.
!    wsave1   Work array for fftopt=2.
!    wsave2   Work array for fftopt=2.
!
!    rhostru  Average rhostr at u points (kg/m**3).
!    rhostrv  Average rhostr at v points (kg/m**3).
!    rhostrw  Average rhostr at w points (kg/m**3).
!
!  OUTPUT:
!
!    pprt     Perturbation pressure at tfuture updated in time
!             by dtsml1 (Pascal)
!    w        Vertical component of Cartesian velocity at tfuture
!             updated in time by dtsml1 (m/s)
!
!  WORK ARRAYS:
!
!    fw       Work array to carry force terms in w-eqation.
!    fp       Work array to carry force terms in p-eqation.
!    tem1     Work array
!    tem2     Work array
!    tem3     Work array
!    tem4     Work array
!    tem5     Work array
!    rstpbi   pbari*rhostr
!    rstptbi  ptbari*rhostr
!    rstwi    1/rhostrw
!    pbzi     2/avgz(pbar)
!    rstj3i   rhostr*j3inv
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

  REAL :: dtsml1               ! Local small time step size (s)
  REAL :: curtsml              ! Local curttim at a small time step

  REAL :: u     (nx,ny,nz)     ! u-velocity at tfuture (m/s)
  REAL :: v     (nx,ny,nz)     ! v-velocity at tfuture (m/s)
  REAL :: w     (nx,ny,nz)     ! w-velocity at tfuture (m/s)
  REAL :: wcont (nx,ny,nz)     ! Contravariant vertical velocity (m/s)
  REAL :: ptprt (nx,ny,nz)     ! Perturbation potential temperature (K)
  REAL :: pprt  (nx,ny,nz)     ! Perturbation pressure at tfuture
                               ! (Pascal)

  REAL :: wdteb (ny,nz)        ! Time tendency of w at east boundary
  REAL :: wdtwb (ny,nz)        ! Time tendency of w at west boundary
  REAL :: wdtnb (nx,nz)        ! Time tendency of w at north boundary
  REAL :: wdtsb (nx,nz)        ! Time tendency of w at south boundary

  REAL :: pdteb (ny,nz)        ! Time tendency of pressure at east
                               ! boundary
  REAL :: pdtwb (ny,nz)        ! Time tendency of pressure at west
                               ! boundary
  REAL :: pdtnb (nx,nz)        ! Time tendency of pressure at north
                               ! boundary
  REAL :: pdtsb (nx,nz)        ! Time tendency of pressure at south
                               ! boundary

  REAL :: phydro(nx,ny)        ! Big time step forcing for computing
                               ! hydrostatic pprt at k=1.

  REAL :: rhostr(nx,ny,nz)     ! Base state density rhobar times j3.
  REAL :: ptbar (nx,ny,nz)     ! Base state potential temperature (K)
  REAL :: ptbari(nx,ny,nz)     ! Inverse base state pot. temperature (K)
  REAL :: pbari (nx,ny,nz)     ! Inverse base state pressure (Pascal).
  REAL :: csndsq(nx,ny,nz)     ! Sound wave speed squared.

  REAL :: wforce(nx,ny,nz)     ! Acoustically inactive forcing terms
                               ! in w-momentum equation (kg/(m*s)**2)
  REAL :: wpgrad(nx,ny,nz)     ! Pressure gradient term in w-eq.
                               ! (kg/(m*s)**2)
  REAL :: pforce(nx,ny,nz)     ! Acoustically inactive forcing terms
                               ! in the pressure equation (Pascal/s)

  REAL :: x     (nx)           ! x-coord. of the physical and compu-
                               ! tational grid. Defined at u-point.
  REAL :: y     (ny)           ! y-coord. of the physical and compu-
                               ! tational grid. Defined at v-point.
  REAL :: z     (nz)           ! z-coord. of the computational grid.
                               ! Defined at w-point on the staggered
                               ! grid.
  REAL :: zp    (nx,ny,nz)     ! Physical height coordinate defined at
                               ! w-point of the staggered grid.

  REAL :: mapfct(nx,ny,8)      ! Map factors at scalar, u and v points

  REAL :: j1    (nx,ny,nz)     ! Coordinate transformation Jacobian
                               ! defined as - d( zp )/d( x ).
  REAL :: j2    (nx,ny,nz)     ! Coordinate transformation Jacobian
                               ! defined as - d( zp )/d( y ).
  REAL :: j3    (nx,ny,nz)     ! Coordinate transformation Jacobian
                               ! defined as d( zp )/d( z ).
  REAL :: aj3z  (nx,ny,nz)     ! Coordinate transformation Jacobian defined
                               ! as d( zp )/d( z ) AVERAGED IN THE Z-DIR.
  REAL :: j3inv (nx,ny,nz)     ! Coordinate transformation Jacobian
                               ! defined as d( zp )/d( z ).

  INTEGER :: exbcbufsz         ! EXBC buffer size
  REAL :: exbcbuf( exbcbufsz ) ! EXBC buffer array

  REAL :: rhofct(nx,ny,nz)     ! rho-factor: rhobar/rho

  REAL :: trigs1(3*(nx-1)/2+1) ! Array containing pre-computed trig
                               ! function for fftopt=1.
  REAL :: trigs2(3*(ny-1)/2+1) ! Array containing pre-computed trig
                               ! function for fftopt=1.
  INTEGER :: ifax1(13)         ! Array containing the factors of nx
                               ! for fftopt=1.
  INTEGER :: ifax2(13)         ! Array containing the factors of ny
                               ! for fftopt=1.

  REAL :: vwork1 (nx+1,ny+1)   ! 2-D work array for fftopt=2 option.
  REAL :: vwork2 (ny,nx+1)     ! 2-D work array for fftopt=2 option.
  REAL :: wsave1 (3*(ny-1)+15) ! Work array for fftopt=2 option.
  REAL :: wsave2 (3*(nx-1)+15) ! Work array for fftopt=2 option.

  REAL :: rhostru(nx,ny,nz)    ! Averaged rhostr at u points (kg/m**3).
  REAL :: rhostrv(nx,ny,nz)    ! Averaged rhostr at v points (kg/m**3).
  REAL :: rhostrw(nx,ny,nz)    ! Averaged rhostr at w points (kg/m**3).
!
!-----------------------------------------------------------------------
!
!  Temporary WORK ARRAYS:
!
!-----------------------------------------------------------------------
!
  REAL :: fp    (nx,ny,nz)     ! Pressure gradient term in w-eq.
  REAL :: fw    (nx,ny,nz)     ! Force in w equation.(kg/(m*s)**2)
  REAL :: tem1  (nx,ny,nz)     ! Temporary work array
  REAL :: tem2  (nx,ny,nz)     ! Temporary work array
  REAL :: tem3  (nx,ny,nz)     ! Temporary work array
  REAL :: tem4  (nx,ny,nz)     ! Temporary work array
  REAL :: tem5  (nx,ny,nz)     ! Temporary work array
  REAL :: rstpbi(nx,ny,nz)     ! pbari*rhostr
  REAL :: rstptbi(nx,ny,nz)    ! ptbari*rhostr
  REAL :: rstwi (nx,ny,nz)     ! 1/rhostrw
  REAL :: pbzi (nx,ny,nz)      ! 2/avgz(pbar)
  REAL :: rstj3i(nx,ny,nz)     ! rhostr*j3inv
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i, j, k
  INTEGER :: ebc1,wbc1,nbc1,sbc1
  REAL :: g05,pk,nk,g05wp

  INTEGER :: wpprt, itema
  REAL :: tema,temb,w1,w2,nrho
  INTEGER :: prgw
  INTEGER :: buoy2swt !Switch for 1st-order or 2nd-order in buoyancy
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
  INCLUDE 'exbc.inc'
  INCLUDE 'phycst.inc'      ! Physical constants
  INCLUDE 'mp.inc'            ! Message passing parameters.
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  IF( bsnesq == 1 ) THEN
    wpprt = 0   ! Switch for pprt/(rhobar*csndsq) term in w-eq.
  ELSE
    wpprt = 1   ! Switch for pprt/(rhobar*csndsq) term in w-eq.
  END IF

  prgw = 0

  g05 = g*0.5
  g05wp = g05*wpprt
!
!-----------------------------------------------------------------------
!
!  Calculate the horizontal velocity divergence using newly updated
!  u and v velocity plus half vertical divergence from wcont, and
!  store the result in tem3.
!
!  Namely, tem3 = difx(u*avgx(j3))+dify(v*avgy(j3))+
!                 (1-tacoef)*difz(wcont*avgz(j3))
!
!-----------------------------------------------------------------------
!
  DO k= 2,nz-2
    DO j= 1,ny-1
      DO i= 1,nx
        tem1(i,j,k)=u(i,j,k)*rhostru(i,j,k)*mapfct(i,j,5)
      END DO
    END DO
  END DO

  DO k= 2,nz-2
    DO j= 1,ny
      DO i= 1,nx-1
        tem2(i,j,k)=v(i,j,k)*rhostrv(i,j,k)*mapfct(i,j,6)
      END DO
    END DO
  END DO

  DO k= 2,nz-1
    DO j= 1,ny-1
      DO i= 1,nx-1
        tem3(i,j,k)=wcont(i,j,k)*rhostrw(i,j,k)
      END DO
    END DO
  END DO

  DO k=2,nz-2
    DO j=1,ny-1
      DO i=1,nx-1
        tem4(i,j,k)=mapfct(i,j,7)                                       &
                   *((tem1(i+1,j,k)-tem1(i,j,k))*dxinv                  &
                    +(tem2(i,j+1,k)-tem2(i,j,k))*dyinv)                 &
                    +(tem3(i,j,k+1)-tem3(i,j,k))*dzinv*(1.0-tacoef)
      END DO
    END DO
  END DO
!
!-----------------------------------------------------------------------
!
!  Compute the forcing terms in pressure equation
!
!  fp = dtsml/j3 *(pforce-csndsq*tem4)
!
!  rhostr*g*w is the base state pressure advection term.
!
!-----------------------------------------------------------------------
!
  DO k=2,nz-2
    DO j=1,ny-1
      DO i=1,nx-1
!    fp(i,j,k)=dtsml1*j3inv(i,j,k) *
!    :            (pforce(i,j,k)-csndsq(i,j,k)*tem4(i,j,k))
        fp(i,j,k)=pforce(i,j,k) -dtsml1*j3inv(i,j,k)*                   &
                                       csndsq(i,j,k)*tem4(i,j,k)
      END DO
    END DO
  END DO

!-----------------------------------------------------------------------
!
!  Compute two more terms in fp related to coordinate transformation:
!
!  tem3 = (avgx(avgz(ustr)*j1)+avgy(avgz(vstr)*j2))/avgz(j3)
!
!-----------------------------------------------------------------------
!
  IF( ternopt /= 0 ) THEN

    DO k= 1,nz-1
      DO j= 1,ny-1
        DO i= 1,nx
          tem1(i,j,k)=u(i,j,k)*rhostru(i,j,k)
        END DO
      END DO
    END DO

    DO k= 1,nz-1
      DO j= 1,ny
        DO i= 1,nx-1
          tem2(i,j,k)=v(i,j,k)*rhostrv(i,j,k)
        END DO
      END DO
    END DO

    DO k= 2,nz-1
      DO j= 1,ny-1
        DO i= 1,nx-1
          tem3(i,j,k)=((tem1(i  ,j,k)+tem1(i  ,j,k-1))*j1(i  ,j,k)      &
                      +(tem1(i+1,j,k)+tem1(i+1,j,k-1))*j1(i+1,j,k)      &
                      +(tem2(i  ,j,k)+tem2(i  ,j,k-1))*j2(i  ,j,k)      &
                      +(tem2(i,j+1,k)+tem2(i,j+1,k-1))*j2(i,j+1,k))     &
                      *mapfct(i,j,8)/aj3z(i,j,k)
        END DO
      END DO
    END DO

!
!-----------------------------------------------------------------------
!
!  fp=fp-(dtsml/j3)*csndsq*tacoef* d(tem3)/dz
!
!-----------------------------------------------------------------------
!
    DO k=2,nz-2
      DO j=1,ny-1
        DO i=1,nx-1
          fp(i,j,k)=fp(i,j,k)                                           &
                     -dtsml1*j3inv(i,j,k)*csndsq(i,j,k)                 &
                     *tacoef*(tem3(i,j,k+1)-tem3(i,j,k))*dzinv
        END DO
      END DO
    END DO

  END IF

!
!-----------------------------------------------------------------------
!
!  Average rhofct at w points, stored in tem5
!
!-----------------------------------------------------------------------
!
  CALL avgsw(rhofct,nx,ny,nz, 1,nx-1, 1,ny-1,tem5)
!
!-----------------------------------------------------------------------
!
!  Compute the right-hand-side forcing term in tridiagonal linear
!  equation. Array w is used to store this forcing term.
!
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!
!  fw=wforce-wpgrad-tacoef*d(fp)/dz
!           -g*avgz(rhostr*(pprt+taceof*pfrc)/(gamma*pbar))
!
!-----------------------------------------------------------------------
!

  IF ( buoy2nd == 0) THEN  !1st-order
    buoy2swt = 0             !Switch for 1st-order or 2nd-order
  ELSE                       !2nd-order
    buoy2swt = 1
  END IF

  tema = 0.5*buoy2swt*(1.0-cpdcv)/cpdcv
  temb = 1.0/cpdcv
  DO k=2,nz-2
    DO j=1,ny-1
      DO i=1,nx-1
!    tem1(i,j,k)=tem6(i,j,k)*((pprt(i,j,k)+tacoef*fp(i,j,k))
        tem1(i,j,k)=rstpbi(i,j,k)*((pprt(i,j,k)+tacoef*fp(i,j,k))       &
                   +tema*pprt(i,j,k)*pprt(i,j,k)                        &
                   *pbari(i,j,k))*temb

      END DO
    END DO
  END DO

!
!-----------------------------------------------------------------------
!
!  When potential temperature equation is solved within small time
!  steps, the contribution from ptprt to buoyancy term is calculated
!  here.
!
!-----------------------------------------------------------------------
!
  IF( ptsmlstp == 1 ) THEN

    DO k=2,nz-2
      DO j=1,ny-1
        DO i=1,nx-1
!      tem2(i,j,k)= ptprt(i,j,k)*tem7(i,j,k)
          tem2(i,j,k)= ptprt(i,j,k)*rstptbi(i,j,k)                      &
                       *(1.0-buoy2swt*ptprt(i,j,k)*ptbari(i,j,k))
        END DO
      END DO
    END DO

  ELSE

    DO k=2,nz-2
      DO j=1,ny-1
        DO i=1,nx-1
          tem2(i,j,k)= 0.0
        END DO
      END DO
    END DO

  END IF

  tema = tacoef*dzinv
  DO k=3,nz-2
    DO j=1,ny-1
      DO i=1,nx-1
        fw(i,j,k) = wforce(i,j,k)-tem5(i,j,k)*(wpgrad(i,j,k)            &
                  +tema*(fp(i,j,k)-fp(i,j,k-1))                         &
                  +g05wp*(tem1(i,j,k)+tem1(i,j,k-1))                    &
                  -g05  *(tem2(i,j,k)+tem2(i,j,k-1)))
      END DO
    END DO
  END DO
!
!-----------------------------------------------------------------------
!
!  fw=fw*dtsml/rhostr + w
!
!-----------------------------------------------------------------------
!
  DO k=3,nz-2
    DO j=1,ny-1
      DO i=1,nx-1
!    fw(i,j,k)=fw(i,j,k)*dtsml1*tem8(i,j,k) + w(i,j,k)
        fw(i,j,k)=fw(i,j,k)*dtsml1*rstwi(i,j,k) + w(i,j,k)
      END DO
    END DO
  END DO

!
!-----------------------------------------------------------------------
!
!  Prepare the left-hand-side coefficents for tridiagonal
!  equation set:
!
!     A(k)*w(k-1)+B(k)*w(k)+C(k)*w(k+1)=D(k)   for k=3,nz-2
!
!  w(2) and w(nz-1) are used as the boundary conditions.
!
!  Due to the lack of work arrays, we are storing the coefficients
!  A in rhostru, B in rhostrv, and C in rhostrw. D is stored in fw.
!
!  rhostru, rhostrv and rhostrw are re-calculated from rhostr after
!  they have been used.
!
!-----------------------------------------------------------------------
!
  DO k=2,nz-1
    DO j=1,ny-1
      DO i=1,nx-1
!    tem4(i,j,k)=0.5*(tem10(i,j,k) +tem10(i,j,k-1))
        tem4(i,j,k)=0.5*(rstj3i(i,j,k)+ rstj3i(i,j,k-1))
      END DO
    END DO
  END DO

  tema = (dtsml1*tacoef)**2 * wpprt*g*dzinv /cpdcv
  temb = (dtsml1*tacoef*dzinv)**2

  DO k=3,nz-2
    DO j=1,ny-1
      DO i=1,nx-1
!    pk = tem5(i,j,k)*tema*tem9(i,j,k)
        pk = tem5(i,j,k)*tema*pbzi(i,j,k)
        nk = tem5(i,j,k)*temb*rstwi(i,j,k)
!    nk = tem5(i,j,k)*temb*tem8(i,j,k)

        tem1(i,j,k)= (-nk+pk)*csndsq(i,j,k-1)*j3inv(i,j,k-1)
        tem3(i,j,k)=-( nk+pk)*csndsq(i,j,k  )*j3inv(i,j,k  )
        tem2(i,j,k)= 1-tem4(i,j,k)*(tem1(i,j,k)+tem3(i,j,k))

        tem1(i,j,k)=tem1(i,j,k)*tem4(i,j,k-1)
        tem3(i,j,k)=tem3(i,j,k)*tem4(i,j,k+1)
      END DO
    END DO
  END DO

!
!-----------------------------------------------------------------------
!
!  Reset fw on the boundaries using the top and bottom boundary
!  conditions of w.
!
!-----------------------------------------------------------------------
!
  DO j=1,ny-1
    DO i=1,nx-1
      w(i,j,nz-1)=0.0
      w(i,j,2) =-((u(i  ,j,2)+u(i  ,j,1))*j1(i,j,2)                     &
                 +(u(i+1,j,2)+u(i+1,j,1))*j1(i+1,j,2)                   &
                 +(v(i,j  ,2)+v(i,j  ,1))*j2(i,j,2)                     &
                 +(v(i,j+1,2)+v(i,j+1,1))*j2(i,j+1,2))                  &
                 *mapfct(i,j,8)

    END DO
  END DO

  DO j=1,ny-1
    DO i=1,nx-1
      fw(i,j,3)   =fw(i,j,   3)-w(i,j,   2)*tem1(i,j,   3)
      fw(i,j,nz-2)=fw(i,j,nz-2)-w(i,j,nz-1)*tem3(i,j,nz-2)
    END DO
  END DO
!
!-----------------------------------------------------------------------
!
!  Call the tridiagonal solver for either a rigid or open upper
!  boundary condition for w.
!
!  tbc = 4, periodic fft for linearized hydrostatic radiation condition.
!  tbc = 6, even fft for linearized hydrostatic radiation condition.
!  The references are Klemp and Durran Jas (1983) and Chen (MWR) 1991.
!
!   w1 is the coef for w(nz-1) and w2 is the coef. for w(nz-2 in the
!   pressure equation.  The pressure equation is solved at p(nz-2
!   for w(nz-1).   tem4(i,j,nz-2) is the summation of the known terms
!   in the pressure equation.  The only unknowns are p(nz-2), w(nz-1)
!   and w(nz-1) at the new time level.  In the tridiagonal solver
!   the elimination step is performed and w(nz-2) is given in termsx
!   of w(nz-1) and known terms.  THe next step is to use the rad.
!   condition (in fourier space):
!
!
!        p = N * rhobar * w(nz-1) / abs(kx+ky)
!
!   The end result is a relation for w(nz-1) given known quantities.
!   Trigs1 and trigs2 are the predetermined trig functions used in the
!   fft program in tridiag.
!
!-----------------------------------------------------------------------

  IF(tbc == 4)THEN ! apply linear radiation condition to w at nz-1.

    tema = rhostr(nx/2,ny/2,nz-2)*csndsq(nx/2,ny/2,nz-2)*               &
                                      j3inv(nx/2,ny/2,nz-2)
    temb = dtsml1*j3inv(nx/2,ny/2,nz-2)
    w2=temb*tacoef*(tema*dzinv+g05*prgw*rhostr(nx/2,ny/2,nz-2))
    w1=temb*tacoef*(-tema*dzinv+g05*prgw*rhostr(nx/2,ny/2,nz-2))
    DO j=1,ny-1
      DO i=1,nx-1
        tem4(i,j,2)=pprt(i,j,nz-2)+fp(i,j,nz-2)
      END DO
    END DO

    nrho=SQRT(g/ptbar(nx/2,ny/2,nz-2)*(ptbar(nx/2,ny/2,nz-1)            &
               -ptbar(nx/2,ny/2,nz-3))*dzinv*0.5)                       &
         *rhostr(nx/2,ny/2,nz-2)*j3inv(nx/2,ny/2,nz-2)

  END IF

!
!-----------------------------------------------------------------------
!
!  NOTE:  tem4(1,1,4) is passed into subroutine tridiag and becomes
!         a 2-D array of size (nx+1,ny+1).
!
!  rhostrw is used as a work array in tridiag.
!
!-----------------------------------------------------------------------
!
  CALL tridiag(nx,ny,nz,tem1,tem2,tem3,fw,rhostrw,                      &
       tem4(1,1,1),tem4(1,1,2),w1,w2,nrho,tem4(1,1,3),                  &
       trigs1,trigs2,ifax1,ifax2,wsave1,wsave2,vwork1,vwork2,           &
       tem4(1,1,4))

  CALL avgsw(rhostr,nx,ny,nz, 1,nx-1, 1,ny-1, rhostrw)


!  call TRIDIAG_old(nx,ny,nz,tem1,tem2,tem3,fw, tem4)

!
!-----------------------------------------------------------------------
!
!  On exit of tridiag, the interior solution of w is stored in fw.
!
!-----------------------------------------------------------------------

  IF(tbc == 4)THEN   ! set itema = nz-1
    itema = nz-1
  ELSE               ! set itema = nz-2
    itema = nz-2
  END IF

  IF( nestgrd /= 1 .OR. mgrid == 1 ) THEN
                                     ! For non-nesting case or the
! base-grid of the nested case

    DO k=3,itema
      DO j=1,ny-1
        DO i=1,nx-1
          w(i,j,k) = fw(i,j,k)
        END DO
      END DO
    END DO

    IF ( lbcopt == 1 ) THEN ! Internal boundary condition
      ebc1=ebc
      wbc1=wbc
      sbc1=sbc
      nbc1=nbc
      IF( ebc == 4 )  ebc1=0
      IF( wbc == 4 )  wbc1=0
      IF( nbc == 4 )  nbc1=0
      IF( sbc == 4 )  sbc1=0

      IF (mp_opt > 0) THEN
        CALL acct_interrupt(mp_acct)
        CALL mpsendrecv2dew(w,nx,ny,nz,ebc1,wbc1,3,tem3)
        CALL mpsendrecv2dns(w,nx,ny,nz,nbc1,sbc1,3,tem3)
      END IF
      CALL acct_interrupt(mp_acct)
      CALL lbcw(nx,ny,nz,dtsml1,w,wcont,wdteb,wdtwb,wdtnb,wdtsb,        &
                ebc1,wbc1,nbc1,sbc1,                                    &
                ebc_global,wbc_global,nbc_global,sbc_global)
      CALL acct_stop_inter

    ELSE  ! External boundary condition

      ebc1 = 3
      wbc1 = 3
      sbc1 = 3
      nbc1 = 3
      IF( ebc == 0 )  ebc1 = 0
      IF( wbc == 0 )  wbc1 = 0
      IF( nbc == 0 )  nbc1 = 0
      IF( sbc == 0 )  sbc1 = 0

      IF (mp_opt > 0) THEN
        CALL acct_interrupt(mp_acct)
        CALL mpsendrecv2dew(w,nx,ny,nz,ebc1,wbc1,3,tem3)
        CALL mpsendrecv2dns(w,nx,ny,nz,nbc1,sbc1,3,tem3)
      END IF
      CALL acct_interrupt(mp_acct)
      CALL lbcw(nx,ny,nz,dtsml1,w,wcont,wdteb,wdtwb,wdtnb,wdtsb,        &
                ebc1,wbc1,nbc1,sbc1,                                    &
                ebc_global,wbc_global,nbc_global,sbc_global)
      CALL acct_stop_inter

    END IF

  ELSE    ! For nested interior grid

    DO k=3,itema
      DO j=2,ny-2
        DO i=2,nx-2
          w(i,j,k) = fw(i,j,k)
        END DO
      END DO
    END DO

    DO k=3,nz-2
      DO j=1,ny-1
        w(   1,j,k)=w(   1,j,k)+dtsml1 * wdtwb(j,k)
        w(nx-1,j,k)=w(nx-1,j,k)+dtsml1 * wdteb(j,k)
      END DO
    END DO
    DO k=3,nz-2
      DO i=2,nx-2
        w(i,   1,k)=w(i,   1,k)+dtsml1 * wdtsb(i,k)
        w(i,ny-1,k)=w(i,ny-1,k)+dtsml1 * wdtnb(i,k)
      END DO
    END DO

  END IF
!
!-----------------------------------------------------------------------
!
!  Calculate wcont at time tfuture, including the boundary
!  points. Wcont at the lateral boundaries is calculated
!  from boundary u, v and w. Wcont at the top and bottom
!  depends on the boundary condition option.
!
!-----------------------------------------------------------------------
!
  CALL wcontra(nx,ny,nz,u,v,w,mapfct,j1,j2,j3,aj3z,                     &
               rhostr,rhostru,rhostrv,rhostrw,wcont,tem1,tem2)
!
!-----------------------------------------------------------------------
!
!  Set the top and bottom boundary conditions for w based on u, v and
!  wcont at the top and bottom boundaries.
!
!-----------------------------------------------------------------------
!
  CALL acct_interrupt(bc_acct)
  CALL vbcw(nx,ny,nz,w,wcont,tbc,bbc,u,v,                               &
            rhostr,rhostru,rhostrv,rhostrw,                             &
            j1,j2,j3)
  CALL acct_stop_inter
!
!-----------------------------------------------------------------------
!
!  Calculate the new pressure
!
!  pprt(new)=pprt(old)+fp-dtsml/j3*tacoef*csndsq*difz(rhobar*w))
!
!-----------------------------------------------------------------------
!
  tema = tacoef*dtsml1*dzinv

  DO k=2,nz-1
    DO j=1,ny-1
      DO i=1,nx-1
!    tem1(i,j,k)= 0.5*(tem10(i,j,k  )+tem10(i,j,k-1))
        tem1(i,j,k)= 0.5*(rstj3i(i,j,k  )+rstj3i(i,j,k-1))
      END DO
    END DO
  END DO

  IF( nestgrd == 1 .AND. mgrid /= 1 ) THEN
                                        ! For nested interior grid

    DO k=2,nz-2
      DO j=2,ny-2
        DO i=2,nx-2
          pprt(i,j,k) = pprt(i,j,k) + fp(i,j,k)                         &
                      - tema*j3inv(i,j,k)*csndsq(i,j,k)                 &
                      * (tem1(i,j,k+1)*w(i,j,k+1)-tem1(i,j,k)*w(i,j,k))
        END DO
      END DO
    END DO

    DO k=2,nz-2
      DO j=1,ny-1
        pprt(   1,j,k)=pprt(   1,j,k)+dtsml1 * pdtwb(j,k)
        pprt(nx-1,j,k)=pprt(nx-1,j,k)+dtsml1 * pdteb(j,k)
      END DO
    END DO
    DO k=2,nz-2
      DO i=2,nx-2
        pprt(i,   1,k)=pprt(i,   1,k)+dtsml1 * pdtsb(i,k)
        pprt(i,ny-1,k)=pprt(i,ny-1,k)+dtsml1 * pdtnb(i,k)
      END DO
    END DO

!
!  Call the pprt bottom boundary condition subroutine to
!  compute the hydrostatic pprt at k=1.
!
    CALL acct_interrupt(bc_acct)
    CALL pprtbbc(nx,ny,nz,g05,buoy2swt,rhostr,pprt,ptprt,               &
                 pbari,ptbari,phydro,                                   &
                 tem1,tem2)          ! tem1 = new pprt at k=1.
    CALL acct_stop_inter

! bc's for mp not implemented for nested grids
    CALL acct_interrupt(mp_acct)
    CALL bcp(nx,ny,nz, dtsml1, pprt, pdteb, pdtwb, pdtnb, pdtsb,        &
             tem1(1,1,1), 0,0,0,0 ,tbc,bbc,                             &
             ebc_global,wbc_global,nbc_global,sbc_global)
    CALL acct_stop_inter

  ELSE                              ! For non-nesting case or the
                                    ! base-grid of the nested case

    DO k=2,nz-2
      DO j=1,ny-1
        DO i=1,nx-1
          pprt(i,j,k) = pprt(i,j,k) + fp(i,j,k)                         &
                      - tema*j3inv(i,j,k)*csndsq(i,j,k)                 &
                      * (tem1(i,j,k+1)*w(i,j,k+1)-tem1(i,j,k)*w(i,j,k))
        END DO
      END DO
    END DO

!
!  Call the pprt bottom boundary condition subroutine to
!  compute the hydrostatic pprt at k=1.
!
    CALL acct_interrupt(bc_acct)
    CALL pprtbbc(nx,ny,nz,g05,buoy2swt,rhostr,pprt,ptprt,               &
                 pbari,ptbari,phydro,                                   &
                 tem1,tem2)          ! tem1 = new pprt at k=1.
    CALL acct_stop_inter

!
!-----------------------------------------------------------------------
!
!  Apply the boundary conditions for the pressure equation.
!
!  For the open boundary case, the boundary value of pprt is
!  predicted by the pressure equation.
!
!-----------------------------------------------------------------------
!
    IF( lbcopt == 1 ) THEN  ! Internal boundary conditions
      ebc1=ebc
      wbc1=wbc
      sbc1=sbc
      nbc1=nbc

      IF( rbc_plbc == 1 ) THEN

        IF( ebc == 4 )  ebc1=0
        IF( wbc == 4 )  wbc1=0
        IF( nbc == 4 )  nbc1=0
        IF( sbc == 4 )  sbc1=0

      ELSE

        IF( ebc == 4 )  ebc1=3
        IF( wbc == 4 )  wbc1=3
        IF( nbc == 4 )  nbc1=3
        IF( sbc == 4 )  sbc1=3

      ENDIF

      IF (mp_opt > 0) THEN
        CALL acct_interrupt(mp_acct)
        CALL mpsendrecv2dew(pprt,nx,ny,nz,ebc1,wbc1,0,tem2)
        CALL mpsendrecv2dns(pprt,nx,ny,nz,nbc1,sbc1,0,tem2)
      END IF
      CALL acct_interrupt(mp_acct)
      CALL bcp(nx,ny,nz, dtsml1, pprt, pdteb, pdtwb, pdtnb, pdtsb,      &
               tem1(1,1,1), ebc1,wbc1,nbc1,sbc1,tbc,bbc,                &
               ebc_global,wbc_global,nbc_global,sbc_global)
      CALL acct_stop_inter

    ELSE  ! External boundary condition

      IF (mp_opt > 0) THEN
        CALL acct_interrupt(mp_acct)
        CALL mpsendrecv2dew(pprt,nx,ny,nz,0,0,0,tem2)
        CALL mpsendrecv2dns(pprt,nx,ny,nz,0,0,0,tem2)
      END IF
      CALL acct_interrupt(mp_acct)
      CALL bcp(nx,ny,nz, dtsml1, pprt, pdteb, pdtwb, pdtnb, pdtsb,      &
               tem1(1,1,1), 0,0,0,0,tbc,bbc,                            &
               ebc_global,wbc_global,nbc_global,sbc_global)
      CALL exbcp(nx,ny,nz, curtsml, pprt,                               &
                 exbcbuf(npr0exb),exbcbuf(nprdtexb))
      CALL acct_stop_inter

    END IF

  END IF

  RETURN
END SUBROUTINE solvwpim1
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE TRIDIAG                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE tridiag(nx,ny,nz,a,b,c,d, tem2, tem1, temxy2,w1,w2,          &
           nrho,work,trigs1,trigs2,ifax1,ifax2,                         &
           wsave1,wsave2,vwork1,vwork2,                                 &
           temxy1)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Solve a tridiagonal linear equation.
!
!  The tridiagonal equation set to be solved is
!
!     a(k)*w(k-1)+b(k)*w(k)+c(k)*w(k+1)=d(k)   for k=3,nz-2
!
!  given w(2) and w(nz-1) as the boundary conditions.
!
!  The solution w(k) is stored directly into d(k).
!
!  Reference: Numerical Recipes: FORTRAN Version, 1989. Page 40.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  9/9/94
!
!  MODIFICATION HISTORY:
!
!  10/31/95 (D. Weber)
!  Added linear hydrostatic upper radiation condition for w-p.
!  References are Klemp and Durran (JAS, 1983) and Chen (1991).
!  Includes the use of trigs1,trigs2,ifax1,ifax2.
!
!  07/22/97 (D. Weber)
!  Added linear hydrostatic upper radiation condition for w-p using
!  an even fft (fftopt=2).
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    nx       Number of grid points in the x-direction
!    ny       Number of grid points in the y-direction
!    nz       Number of grid points in the vertical
!
!    a        lhs coefficient of tridigonal eqations
!    b        lhs coefficient of tridigonal eqations
!    c        lhs coefficient of tridigonal eqations
!    d        rhs coefficient of tridigonal eqations
!    temxy2   Pforce array at nz-2.
!    trigs1   Array containing pre-computed trig function for
!             fftopt=1.
!    trigs2   Array containing pre-computed trig function for
!             fftopt=1.
!    ifax1    Array containing the factors of nx for fftopt=1.
!    ifax2    Array containing the factors of ny for fftopt=1.
!
!    vwork1   2-D work array for fftopt=2 option.
!    vwork2   2-D work array for fftopt=2 option.
!    wsave1   Work array for fftopt=2.
!    wsave2   Work array for fftopt=2.
!
!  OUTPUT:
!
!    d        The solution w.
!
!  WORK ARRAYS:
!
!    temxy1   Work array at nz-2.
!    tem1     Work array
!    tem2     Work array
!    work     Work array for fft.
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

!
!  Include files:
!
  INCLUDE 'bndry.inc'
  INCLUDE 'globcst.inc'

  INTEGER :: nx, ny, nz        ! Number of grid points in 3 directions

  REAL :: a     (nx,ny,nz)     ! lhs coefficient of tridigonal eqations
  REAL :: b     (nx,ny,nz)     ! lhs coefficient of tridigonal eqations
  REAL :: c     (nx,ny,nz)     ! lhs coefficient of tridigonal eqations
  REAL :: d     (nx,ny,nz)     ! rhs coefficient of tridigonal eqations
                               ! The contains solution w on exit.
  REAL :: temxy2(nx,ny)        ! Pforce array at nz-2.

  REAL :: trigs1(3*(nx-1)/2+1) ! Array containing pre-computed trig
                               ! function for fftopt=1.
  REAL :: trigs2(3*(ny-1)/2+1) ! Array containing pre-computed trig
                               ! function for fftopt=1.
  INTEGER :: ifax1(13)         ! Array containing the factors of nx for
                               ! fftopt=1.
  INTEGER :: ifax2(13)         ! Array containing the factors of ny for
                               ! fftopt=1.

  REAL :: vwork1 (nx+1,ny+1)   ! 2-D work array for fftopt=2.
  REAL :: vwork2 (ny,nx+1)     ! 2-D work array for fftopt=2.
  REAL :: wsave1 (3*(ny-1)+15) ! Work array for fftopt=2.
  REAL :: wsave2 (3*(nx-1)+15) ! Work array for fftopt=2.

  REAL :: temxy1(nx+1,ny+1)    ! Work array at nz-2.
  REAL :: tem1  (nx,ny)        ! 2-D work array.
  REAL :: tem2  (nx,ny,nz)     ! 3-D work array.
  REAL :: work  (ny,nx)        ! 2-D work array for fft code.

!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i, j, k
  REAL :: nrho
  REAL :: w1,w2,wc
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  DO j=1,ny-1
    DO i=1,nx-1
      d(i,j,3)=d(i,j,3)/b(i,j,3)
      tem1(i,j)=b(i,j,3)
    END DO
  END DO

  DO k=4,nz-2
    DO j=1,ny-1
      DO i=1,nx-1
        tem2(i,j,k) = c(i,j,k-1)/tem1(i,j)
        tem1(i,j)=b(i,j,k)-a(i,j,k)*tem2(i,j,k)
        d(i,j,k)=(d(i,j,k)-a(i,j,k)*d(i,j,k-1))/tem1(i,j)
      END DO
    END DO
  END DO


  IF(tbc == 4)THEN   !  apply upper radiation boundary condition.

!
!  Computing the new c(nz-2) (which is stored in nz-1)....
!

    DO j=1,ny-1
      DO i=1,nx-1
        tem2(i,j,nz-1) = c(i,j,nz-2)/tem1(i,j)
      END DO
    END DO

    wc = tem2(nx/2,ny/2,nz-1)  ! coef. for w(nz-1) in tri. elim. phase

!
!  Combining the pforce array and d(nz-2) from the elimination phase.
!  Note: these arrays are a function of x and/or y.  All others variables
!  are determined from base state which are not a function
!  of x,y.  (note: zflat must be at or below the height of scalar nz-2)
!

    DO j=1,ny-1
      DO i=1,nx-1
        temxy1(i,j) = temxy2(i,j) + w2*d(i,j,nz-2)
      END DO
    END DO

!-----------------------------------------------------------------------
!
!  Call the upper radiation subroutine to update w at the top (nz-2)
!
!-----------------------------------------------------------------------

    IF(fftopt == 2)THEN     !  call the even vfftpack fft...

      CALL uprad3(nx,ny,w1,w2,wc,nrho,wsave1,wsave2,vwork1,vwork2,      &
                  temxy1,work)

    ELSE IF(fftopt == 1)THEN !  call the periodic fft....

      CALL uprad1(nx,ny,w1,w2,wc,nrho,ifax1,ifax2,trigs1,trigs2,        &
                  temxy1,work)

    END IF  ! end of fftopt....

    DO j=1,ny-1
      DO i=1,nx-1
        d(i,j,nz-1) = temxy1(i,j)
      END DO
    END DO

!-----------------------------------------------------------------------
!
!  Perform the back substitution phase of the tridiagonal solver.
!
!-----------------------------------------------------------------------

    DO k=nz-2,3,-1
      DO j=1,ny-1
        DO i=1,nx-1
          d(i,j,k)=d(i,j,k)-tem2(i,j,k+1)*d(i,j,k+1)
        END DO
      END DO
    END DO

  ELSE IF(tbc /= 4)THEN     ! perform backsubtitution for other bc's.

    DO k=nz-3,3,-1
      DO j=1,ny-1
        DO i=1,nx-1
          d(i,j,k)=d(i,j,k)-tem2(i,j,k+1)*d(i,j,k+1)
        END DO
      END DO
    END DO


  END IF  ! end of tbc loop...


  RETURN
END SUBROUTINE tridiag
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE SOLVPT_LRG                 ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE solvpt_lrg(nx,ny,nz, exbcbufsz, dtbig1,                      &
           ptprt, ptdteb,ptdtwb,ptdtnb,ptdtsb,                          &
           rhostr,rhostri,ptforce, exbcbuf, tem1)

!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Coordinate the time integration of the potential temperature
!  equation in large time steps.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  10/10/91
!
!  MODIFICATION HISTORY:
!
!  9/17/94 (M. Xue)
!  Rewritten for small time step integration of ptprt equation.
!
!  11/6/1995 (M. Xue)
!  Added option for fourth order vertical advection for ptbar.
!
!  4/17/96 (M. Xue)
!  Removed the block for 4th order advection of ptbar.
!
!  4/24/1997 (M. Xue)
!  Rewrote as SOLVPT_LRG based on SOLVPT.
!
!  10/5/1998 (Dan Weber)
!  Added rhostri for efficiency.
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!    dtbig1   The big time step size for this call (s)
!
!    ptprt    Perturbation potential temperature at times tpast and
!             tpresent (K)
!
!    ptdteb   Time tendency of the ptprt field at the east boundary
!    ptdtwb   Time tendency of the ptprt field at the west boundary
!    ptdtnb   Time tendency of the ptprt field at the north boundary
!    ptdtsb   Time tendency of the ptprt field at the south boundary
!
!    rhostr   Base state density rhobar times j3 (kg/m**3)
!    rhostri  Inverse base state density rhobar times j3 (kg/m**3)
!
!    ptforce  Gravity wave inactive forcing terms in pt-eq.
!             (K*kg/(m**3*s))
!
!  OUTPUT:
!
!    ptprt    Perturbation potential temperature at time tfuture (K)
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

  REAL :: dtbig1               ! The big time step size for this call
                               ! (s)

  REAL :: ptprt (nx,ny,nz,nt)  ! Perturbation potential temperature (K)

  REAL :: ptdteb(ny,nz)        ! Time tendency of ptprt at east
                               ! boundary
  REAL :: ptdtwb(ny,nz)        ! Time tendency of ptprt at west
                               ! boundary
  REAL :: ptdtnb(nx,nz)        ! Time tendency of ptprt at north
                               ! boundary
  REAL :: ptdtsb(nx,nz)        ! Time tendency of ptprt at south
                               ! boundary

  REAL :: rhostr(nx,ny,nz)     ! Base state density rhobar times j3.
  REAL :: rhostri(nx,ny,nz)    ! Inverse base state density rhobar times j3.

  REAL :: ptforce(nx,ny,nz)    ! Gravity wave inactive forcing terms
                               ! in potential temperature eq.
                               ! (K*kg/(m**3*s))

  INTEGER :: exbcbufsz         ! EXBC buffer size
  REAL :: exbcbuf( exbcbufsz ) ! EXBC buffer array
  REAL :: tem1(nx,ny,nz)       ! Work array
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i, j, k, tstrtlvl
  INTEGER :: ebc1,wbc1,nbc1,sbc1
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
  INCLUDE 'mp.inc'            ! Message passing parameters.
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
!  Integrate forward by one timestep the potential temperature
!  equation. When PT-eq. is integrated inside small time steps,
!  it is stepped forward by a small time step, otherwise, it is
!  stepped forward by a large time step using leapfrog scheme
!  (i.e. 2*dtbig1 from ptprt at tpast).
!
!-----------------------------------------------------------------------
!
  IF(sadvopt == 4 .OR. tintegopt == 2 .OR. tintegopt == 3) THEN ! Forward-based FCT scheme
    deltat = dtbig1
    tstrtlvl = tpresent
  ELSE
    deltat = dtbig1*2
    tstrtlvl = tpast
  END IF

  IF( nestgrd == 1 .AND. mgrid /= 1 ) THEN

    DO k=2,nz-2
      DO j=2,ny-2
        DO i=2,nx-2
          ptprt(i,j,k,tfuture)=ptprt(i,j,k,tstrtlvl)                    &
                              +deltat*ptforce(i,j,k)*rhostri(i,j,k)
        END DO
      END DO
    END DO
!
!-----------------------------------------------------------------------
!
!  Integrate PT equation and the boundary conditions for a nested grid.
!
!-----------------------------------------------------------------------
!
    DO k=2,nz-2
      DO i=1,nx-1
        ptprt(i,   1,k,tfuture)=                                        &
            2*ptprt(i,   1,k,tpresent)-ptprt(i,   1,k,tpast)
        ptprt(i,ny-1,k,tfuture)=                                        &
            2*ptprt(i,ny-1,k,tpresent)-ptprt(i,ny-1,k,tpast)
      END DO
    END DO

    DO k=2,nz-2
      DO j=2,ny-2
        ptprt(   1,j,k,tfuture)=                                        &
            2*ptprt(   1,j,k,tpresent)-ptprt(   1,j,k,tpast)
        ptprt(nx-1,j,k,tfuture)=                                        &
            2*ptprt(nx-1,j,k,tpresent)-ptprt(nx-1,j,k,tpast)
      END DO
    END DO

! bc's for mp not implemented for nested grids
    CALL acct_interrupt(mp_acct)
    CALL bcsclr(nx,ny,nz, deltat*0.5 ,                                  &
                ptprt(1,1,1,tstrtlvl),ptprt(1,1,1,tpresent),            &
                ptprt(1,1,1,tfuture),ptdteb,ptdtwb,ptdtnb,ptdtsb,       &
                0,0,0,0 ,tbc,bbc,                                       &
                ebc_global,wbc_global,nbc_global,sbc_global)
    CALL acct_stop_inter

  ELSE
!
!-----------------------------------------------------------------------
!
!  Integrate PT equation and the boundary conditions for a base grid.
!
!-----------------------------------------------------------------------
!
    DO k=2,nz-2
      DO j=1,ny-1
        DO i=1,nx-1
          ptprt(i,j,k,tfuture)=ptprt(i,j,k,tstrtlvl)                    &
                              +deltat*ptforce(i,j,k)*rhostri(i,j,k)
        END DO
      END DO
    END DO

    IF ( lbcopt == 1 ) THEN ! Internal boundary conditions

      ebc1=ebc
      wbc1=wbc
      sbc1=sbc
      nbc1=nbc
      IF( ebc == 4 )  ebc1=0
      IF( wbc == 4 )  wbc1=0
      IF( nbc == 4 )  nbc1=0
      IF( sbc == 4 )  sbc1=0

      IF (mp_opt > 0) THEN
        CALL acct_interrupt(mp_acct)
        CALL mpsendrecv2dew(ptprt(1,1,1,tfuture),nx,ny,nz,ebc1,wbc1,0,tem1)
        CALL mpsendrecv2dns(ptprt(1,1,1,tfuture),nx,ny,nz,nbc1,sbc1,0,tem1)
      END IF
      CALL acct_interrupt(mp_acct)
      CALL bcsclr(nx,ny,nz, deltat*0.5 ,                                &
                  ptprt(1,1,1,tstrtlvl),ptprt(1,1,1,tpresent),          &
                  ptprt(1,1,1,tfuture),ptdteb,ptdtwb,ptdtnb,ptdtsb,     &
                  ebc1,wbc1,nbc1,sbc1,tbc,bbc,                          &
                  ebc_global,wbc_global,nbc_global,sbc_global)
      CALL acct_stop_inter

    ELSE  ! External boundary condition

      IF (mp_opt > 0) THEN
        CALL acct_interrupt(mp_acct)
        CALL mpsendrecv2dew(ptprt(1,1,1,tfuture),nx,ny,nz,0,0,0,tem1)
        CALL mpsendrecv2dns(ptprt(1,1,1,tfuture),nx,ny,nz,0,0,0,tem1)
      END IF
      CALL acct_interrupt(mp_acct)
      CALL bcsclr(nx,ny,nz, deltat*0.5 ,                                &
                  ptprt(1,1,1,tstrtlvl),ptprt(1,1,1,tpresent),          &
                  ptprt(1,1,1,tfuture),ptdteb,ptdtwb,ptdtnb,ptdtsb,     &
                  0,0,0,0,tbc,bbc,                                      &
                  ebc_global,wbc_global,nbc_global,sbc_global)

      CALL exbcpt(nx,ny,nz, curtim+dtbig, ptprt(1,1,1,tfuture),         &
                  exbcbuf(npt0exb),exbcbuf(nptdtexb))
      CALL acct_stop_inter

    END IF

  END IF

  RETURN
END SUBROUTINE solvpt_lrg
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE SOLVPT_SML                 ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE solvpt_sml(nx,ny,nz, exbcbufsz, dtbig1,dtsml1,curtsml,       &
           ptprt,w, ptdteb,ptdtwb,ptdtnb,ptdtsb,                        &
           ptbar,rhostr,rhostri,rhostrw,j3,aj3z,ptforce,                &
           exbcbuf,                                                     &
           ptadv,tem1)

!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Coordinate the time integration of the potential temperature
!  equation in small time steps.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  10/10/91
!
!  MODIFICATION HISTORY:
!
!  9/17/94 (M. Xue)
!  Rewritten for small time step integration of ptprt equation.
!
!  11/6/1995 (M. Xue)
!  Added option for fourth order vertical advection for ptbar.
!
!  4/17/96 (M. Xue)
!  Removed the block for 4th order advection of ptbar.
!
!  4/24/1997 (M. Xue)
!  Rewrote as SOLVPT_SML based on SOLVPT.
!
!  9/18/1998 (D. Weber)
!  Added array aj3z and rhostri,w for efficiency.
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!    dtbig1   The big time step size for this call (s)
!    dtsml1   The big time step size for this call (s)
!    curtsml  Current time during small time step integration.
!
!    ptprt    Perturbation potential temperature at times tpast and
!             tpresent (K)
!    w        Vertical component of Cartesian velocity
!
!    ptdteb   Time tendency of the ptprt field at the east boundary
!    ptdtwb   Time tendency of the ptprt field at the west boundary
!    ptdtnb   Time tendency of the ptprt field at the north boundary
!    ptdtsb   Time tendency of the ptprt field at the south boundary
!
!    ptbar    Base state potential temperature (K)
!    rhostr   Base state density rhobar times j3 (kg/m**3)
!    rhostri  Inverse base state density rhobar times j3 (kg/m**3)
!    rhostrw  rhostr averaged to w-point.
!
!    j3       Coordinate transformation Jacobian  d(zp)/dz
!    aj3z     Avgz of the coordinate transformation Jacobian  d(zp)/dz
!
!    ptforce  Gravity wave inactive forcing terms in pt-eq.
!             (K*kg/(m**3*s))
!
!  OUTPUT:
!
!    ptprt    Perturbation potential temperature at time tfuture (K)
!
!  WORK ARRAYS:
!
!    ptadv    Advection of base state potential temperature
!    tem1     Temporary work array.
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

  REAL :: dtbig1               ! The big time step size for this call
                               ! (s)
  REAL :: dtsml1               ! The big time step size for this call
                               ! (s)
  REAL :: curtsml              ! Current time during small time step
                               ! integration.

  REAL :: ptprt (nx,ny,nz,nt)  ! Perturbation potential temperature (K)
  REAL :: w     (nx,ny,nz,nt)  ! Total w-velocity (m/s)

  REAL :: ptdteb(ny,nz)        ! Time tendency of ptprt at east
                               ! boundary
  REAL :: ptdtwb(ny,nz)        ! Time tendency of ptprt at west
                               ! boundary
  REAL :: ptdtnb(nx,nz)        ! Time tendency of ptprt at north
                               ! boundary
  REAL :: ptdtsb(nx,nz)        ! Time tendency of ptprt at south
                               ! boundary

  REAL :: ptbar (nx,ny,nz)     ! Base state potential temperature (K)
  REAL :: rhostr(nx,ny,nz)     ! Base state density rhobar times j3.
  REAL :: rhostri(nx,ny,nz)    ! Inverse base state density rhobar times j3.
  REAL :: rhostrw(nx,ny,nz)    ! rhostr averaged to w-point.

  REAL :: j3    (nx,ny,nz)     ! Coordinate transformation Jacobian
                               ! defined as d( zp )/d( z ).
  REAL :: aj3z  (nx,ny,nz)     ! Coordinate transformation Jacobian defined
                               ! as d( zp )/d( z ) AVERAGED IN THE Z-DIR.
  REAL :: ptforce(nx,ny,nz)    ! Gravity wave inactive forcing terms
                               ! in potential temperature eq.
                               ! (K*kg/(m**3*s))

  INTEGER :: exbcbufsz         ! EXBC buffer size
  REAL :: exbcbuf( exbcbufsz ) ! EXBC buffer array

  REAL :: ptadv (nx,ny,nz)     ! Temporary array to store base state
                               ! potential temperature advection.
  REAL :: tem1  (nx,ny,nz)     ! Temporary work array
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i, j, k, onvf
  INTEGER :: ebc1,wbc1,nbc1,sbc1
  REAL :: deltat, dz05 , time

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
  INCLUDE 'exbc.inc'
  INCLUDE 'mp.inc'            ! Message passing parameters.
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
!  Integrate forward by one timestep the potential temperature
!  equation. When PT-eq. is integrated inside small time steps,
!  it is stepped forward by a small time step, otherwise, it is
!  stepped forward by a large time step using leapfrog scheme
!  (i.e. 2*dtbig1 from ptprt at tpast).
!
!-----------------------------------------------------------------------
!

  deltat = dtsml1
  time = curtsml

!
!-----------------------------------------------------------------------
!
!  Base state potential temperature advection.  This term is added to
!  the array ptadv to yield the total potential temperature advection.
!
!  ptbar is assumed to be independent of physical x and y, therefore
!  d(ptbar)/dx and d(ptbar)/dy for constant z are zero, the base state
!  advection is -w*d(ptbar)/dzp = -w/j3*d(ptbar)/dz.
!
!  This term is responsible for internal gravity waves, when potential
!  temperature equation is solved inside the small time steps, this
!  term is evaluated inside the small time steps.
!
!-----------------------------------------------------------------------
!
!  dz05 = dz*0.5

  DO k=2,nz-1
    DO j=1,ny-1
      DO i=1,nx-1
        tem1(i,j,k)=rhostrw(i,j,k)*w(i,j,k,tfuture)                     &
                   *(ptbar(i,j,k)-ptbar(i,j,k-1))                       &
                   *dzinv/aj3z(i,j,k)
      END DO
    END DO
  END DO

  onvf = 0
  CALL avgz(tem1 , onvf,                                                &
       nx,ny,nz, 1,nx-1, 1,ny-1, 2,nz-2, ptadv)

  IF( nestgrd == 1 .AND. mgrid /= 1 ) THEN

    DO k=2,nz-2
      DO j=2,ny-2
        DO i=2,nx-2
          ptprt(i,j,k,tfuture)=ptprt(i,j,k,tfuture)                     &
              +deltat*(ptforce(i,j,k)-ptadv(i,j,k))*rhostri(i,j,k)
        END DO
      END DO
    END DO
!
!-----------------------------------------------------------------------
!
!  Integrate PT equation and the boundary conditions for a nested grid.
!
!-----------------------------------------------------------------------
!
    DO k=2,nz-2
      DO j=1,ny-1
        ptprt(   1,j,k,tfuture)=ptprt(   1,j,k,tfuture)                 &
                               +deltat* ptdtwb(j,k)
        ptprt(nx-1,j,k,tfuture)=ptprt(nx-1,j,k,tfuture)                 &
                               +deltat* ptdteb(j,k)
      END DO
    END DO

    DO k=2,nz-2
      DO i=2,nx-2
        ptprt(i,   1,k,tfuture)=ptprt(i,   1,k,tfuture)                 &
                               +deltat* ptdtsb(i,k)
        ptprt(i,ny-1,k,tfuture)=ptprt(i,ny-1,k,tfuture)                 &
                               +deltat* ptdtnb(i,k)
      END DO
    END DO

    ! bc's for mp not implemented for nested grids
    CALL acct_interrupt(mp_acct)
    CALL bcsclr(nx,ny,nz, deltat*0.5 ,                                  &
                ptprt(1,1,1,tfuture),ptprt(1,1,1,tpresent),             &
                ptprt(1,1,1,tfuture),ptdteb,ptdtwb,ptdtnb,ptdtsb,       &
                0,0,0,0 ,tbc,bbc,                                       &
                ebc_global,wbc_global,nbc_global,sbc_global)
    CALL acct_stop_inter

  ELSE
!
!-----------------------------------------------------------------------
!
!  Integrate PT equation and the boundary conditions for a base grid.
!
!-----------------------------------------------------------------------
!
    DO k=2,nz-2
      DO j=1,ny-1
        DO i=1,nx-1
          ptprt(i,j,k,tfuture)=ptprt(i,j,k,tfuture)                     &
              +deltat*(ptforce(i,j,k)-ptadv(i,j,k))*rhostri(i,j,k)
        END DO
      END DO
    END DO

    IF ( lbcopt == 1 ) THEN ! Internal boundary conditions

      ebc1=ebc
      wbc1=wbc
      sbc1=sbc
      nbc1=nbc
      IF( ebc == 4 )  ebc1=0
      IF( wbc == 4 )  wbc1=0
      IF( nbc == 4 )  nbc1=0
      IF( sbc == 4 )  sbc1=0

      IF (mp_opt > 0) THEN
        CALL acct_interrupt(mp_acct)
        CALL mpsendrecv2dew(ptprt(1,1,1,tfuture),nx,ny,nz,ebc1,wbc1,0,tem1)
        CALL mpsendrecv2dns(ptprt(1,1,1,tfuture),nx,ny,nz,nbc1,sbc1,0,tem1)
      END IF
      CALL acct_interrupt(mp_acct)
      CALL bcsclr(nx,ny,nz, deltat*0.5 ,                                &
                  ptprt(1,1,1,tfuture),ptprt(1,1,1,tpresent),           &
                  ptprt(1,1,1,tfuture),ptdteb,ptdtwb,ptdtnb,ptdtsb,     &
                  ebc1,wbc1,nbc1,sbc1,tbc,bbc,                          &
                  ebc_global,wbc_global,nbc_global,sbc_global)
      CALL acct_stop_inter

    ELSE  ! External boundary condition

      IF (mp_opt > 0) THEN
        CALL acct_interrupt(mp_acct)
        CALL mpsendrecv2dew(ptprt(1,1,1,tfuture),nx,ny,nz,0,0,0,tem1)
        CALL mpsendrecv2dns(ptprt(1,1,1,tfuture),nx,ny,nz,0,0,0,tem1)
      END IF
      CALL acct_interrupt(mp_acct)
      CALL bcsclr(nx,ny,nz, deltat*0.5 ,                                &
                  ptprt(1,1,1,tfuture),ptprt(1,1,1,tpresent),           &
                  ptprt(1,1,1,tfuture),ptdteb,ptdtwb,ptdtnb,ptdtsb,     &
                  0,0,0,0,tbc,bbc,                                      &
                  ebc_global,wbc_global,nbc_global,sbc_global)

      CALL exbcpt(nx,ny,nz, time , ptprt(1,1,1,tfuture),                &
                  exbcbuf(npt0exb),exbcbuf(nptdtexb))
      CALL acct_stop_inter

    END IF

  END IF

  RETURN
END SUBROUTINE solvpt_sml
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE SOLVQV                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
!

SUBROUTINE solvqv(RK3step,nx,ny,nz, exbcbufsz, dtbig1,                  &
         qv,u,v,wcont, ustr,vstr,wstr,                                  &
         qvdteb,qvdtwb,qvdtnb,qvdtsb,                                   &
         rhostr,rhostri,qvbar,kmh,kmv,rprntl,qvsflx,pbldpth,            &
         x,y,z,zp,mapfct,j1,j2,j3,aj3x,aj3y,j3inv,qvcumsrc,             &
         usflx,vsflx,ptsflx,ptsfc,qvsfc,ptbar,ptprt,                    &
         exbcbuf,bcrlx,                                                 &
         qadv,qmix,                                                     &
         tem1,tem2,tem3,tem4,tem5,tem6,tem7,tem8,                       &
         tem1_0,tem2_0,tem3_0)


!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Coordinate the time integration of the equation for water vapor
!  specific humidity qv, and the setting of boundary conditions.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  10/10/91
!
!  MODIFICATION HISTORY:
!
!  5/20/92 (M. Xue)
!  Added full documentation.
!
!  2/10/93 (K. Droegemeier)
!  Cleaned up documentation.
!
!  4/10/93 (M. Xue & Hao Jin)
!  Add the terrain.
!
!  8/22/95 (M. Xue)
!  Added ptcumsrc term to the right hand side of qv equation.
!  It was omitted before.
!
!  08/30/95 (M. Xue and Y. Liu)
!  Changed EXBC for water variables to zero-gradient when the
!  variables are missing in the boundary file.
!
!  08/30/95 (Yuhe Liu)
!  Fixed a bug which double computes the qv advection and mixing.
!
!  1/25/96 (Donghai Wang & Yuhe Liu)
!  Added the map projection factor to ARPS governing equations.
!
!  4/1/96 (Donghai Wang, X. Song and M. Xue)
!  Added the implicit treatment for the vertical mixing.
!
!  3/24/1997 (M. Xue)
!  Code to handle the case of forward time integration added.
!
!  7/10/1997 (Fanyou Kong -- CMRP)
!  Added the positive definite advection scheme option (sadvopt = 5)
!
!  7/17/1998 (M. Xue)
!  Changed call to ADVQFCT.
!
!  9/18/1998 (D. Weber)
!  Added arrays aj3x,y, u,v,wstr,rhostri do not alter!
!
! 10/19/2009 (D. Dawson)
!  Updated to accomodate RK3 time integration
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    RK3step  Which RK3 step are we on (1, 2, 3)?
!             (only used for tintegopt == 2)
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!    dtbig1   The big time step size (s)
!
!    qv       Water vapor specific humidity at tpast and tpresent
!             (kg/kg)
!
!    u        x component of velocity at all time levels (m/s)
!    v        y component of velocity at all time levels (m/s)
!    wcont    Contravariant vertical velocity (m/s)
!
!    ustr     Appropriate ustr for advection.
!    vstr     Appropriate vstr for advection.
!    wstr     Appropriate wstr for advection.
!
!    qvdteb   Time tendency of the qv field at the east boundary
!    qvdtwb   Time tendency of the qv field at the west boundary
!    qvdtnb   Time tendency of the qv field at the north boundary
!    qvdtsb   Time tendency of the qv field at the south boundary
!
!    rhostr   Base state density rhobar times j3 (kg/m**3)
!    rhostri  Inverse base state density rhobar times j3 (kg/m**3)
!    qvbar    Base state water vapor specific humidity (kg/kg)
!
!    kmh      Horizontal turb. mixing coef. for momentum ( m**2/s )
!    kmv      Vertical turb. mixing coef. for momentum ( m**2/s )
!    rprntl   Reciprocal of Prandtl number
!
!    qvsflx   Surface flux of moisture (kg/(m**2*s)).
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
!    qvcumsrc Source term in qv-equation due to cumulus parameterization
!
!  OUTPUT:
!
!    qv       Water vapor specific humidity at time tfuture (kg/kg)
!
!  WORK ARRAYS:
!
!    qadv     Advection term of water vapor eq. (kg/(m**3*s)).
!             A local array.
!    qmix     Total mixing in water vapor eq. (kg/(m**3*s)).
!             A local array.
!    tem1     Temporary work array.
!    tem2     Temporary work array.
!    tem3     Temporary work array.
!    tem4     Temporary work array.
!    tem5     Temporary work array.
!    tem6     Temporary work array.
!    tem6     Temporary work array.
!    tem7     Temporary work array.
!    tem8     Temporary work array.
!
!    tem1_0   Temporary work array.
!    tem2_0   Temporary work array.
!    tem3_0   Temporary work array.
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

  INTEGER :: RK3step
  INTEGER :: nx, ny, nz        ! Number of grid points in 3 directions
  REAL :: dtbig1               ! Local big time step

  REAL :: qv    (nx,ny,nz,nt)  ! Water vapor specific humidity (kg/kg)
  REAL :: u     (nx,ny,nz,nt)  ! Total u-velocity (m/s)
  REAL :: v     (nx,ny,nz,nt)  ! Total v-velocity (m/s)
  REAL :: wcont (nx,ny,nz)     ! Contravariant vertical velocity (m/s)
  REAL :: ustr  (nx,ny,nz)     ! Appropriate ustr for advection.
  REAL :: vstr  (nx,ny,nz)     ! Appropriate vstr for advection.
  REAL :: wstr  (nx,ny,nz)     ! Appropriate wstr for advection.

  REAL :: qvdteb(ny,nz)        ! Time tendency of qv at east boundary
  REAL :: qvdtwb(ny,nz)        ! Time tendency of qv at west boundary
  REAL :: qvdtnb(nx,nz)        ! Time tendency of qv at north boundary
  REAL :: qvdtsb(nx,nz)        ! Time tendency of qv at south boundary

  REAL :: rhostr(nx,ny,nz)     ! Base state density rhobar times j3.
  REAL :: rhostri(nx,ny,nz)    ! Inverse base state density rhobar times j3.
  REAL :: qvbar (nx,ny,nz)     ! Base state water vapor specific
                               ! humidity (kg/kg)

  REAL :: kmh   (nx,ny,nz)     ! Horizontal turb. mixing coef. for
                               ! momentum. ( m**2/s )
  REAL :: kmv   (nx,ny,nz)     ! Vertical turb. mixing coef. for
                               ! momentum. ( m**2/s )
  REAL :: rprntl(nx,ny,nz)     ! Reciprocal of Prandtl number

  REAL :: usflx(nx,ny)         ! Surface flux of u-momentum
  REAL :: vsflx(nx,ny)         ! Surface flux of v-momentum
  REAL :: ptsflx(nx,ny)        ! Surface heat flux (K*kg/(m*s**2))
  REAL :: qvsflx(nx,ny)        ! Surface flux of moisture (kg/(m**2*s))
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

  REAL :: qvcumsrc(nx,ny,nz)   ! Source in qv-equation due to cumulus
                               ! parameterization

  REAL :: ptbar (nx,ny,nz)     ! Base state potential temperature (K)
  REAL :: ptprt (nx,ny,nz,nt)  ! Perturbation potential temperature

  REAL :: ptsfc(nx,ny)         ! Potential temperature at the ground level (K)
  REAL :: qvsfc(nx,ny)         ! Effective qv at the surface (kg/kg)

  INTEGER :: exbcbufsz         ! EXBC buffer size
  REAL :: exbcbuf( exbcbufsz ) ! EXBC buffer array
  REAL :: bcrlx (nx,ny)        ! EXBC relaxation coefficients

  REAL :: qadv(nx,ny,nz)       ! Advection of water vapor (kg/(m**3*s))
                               ! A local array.
  REAL :: qmix(nx,ny,nz)       ! Total mixing of water vapor
                               ! (kg/(m**3*s))
                               ! A local array.
  REAL :: tem1  (nx,ny,nz)     ! Temporary work array
  REAL :: tem2  (nx,ny,nz)     ! Temporary work array
  REAL :: tem3  (nx,ny,nz)     ! Temporary work array
  REAL :: tem4  (nx,ny,nz)     ! Temporary work array
  REAL :: tem5  (nx,ny,nz)     ! Temporary work array
  REAL :: tem6  (nx,ny,nz)     ! Temporary work array
  REAL :: tem7  (nx,ny,nz)     ! Temporary work array
  REAL :: tem8  (nx,ny,nz)     ! Temporary work array

  REAL :: tem1_0(0:nx,0:ny,0:nz)     ! Temporary work array.
  REAL :: tem2_0(0:nx,0:ny,0:nz)     ! Temporary work array.
  REAL :: tem3_0(0:nx,0:ny,0:nz)     ! Temporary work array.
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k
  INTEGER :: tstrtlvl  ! Starting level of time integration
  REAL :: deltat
  INTEGER :: ebc1,wbc1,nbc1,sbc1, qflag
  LOGICAL :: RK3sadvflag      ! Flag to indicate if sadvopt == 4 and RK3
                              ! time-stepping is being used.  It's used
                              ! to allow the normal advection to be called
                              ! for the first 2 RK3 steps, and the full
                              ! FCT only for the last one.
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'bndry.inc'
  INCLUDE 'exbc.inc'
  INCLUDE 'globcst.inc'
  INCLUDE 'mp.inc'            ! Message passing parameters.
!
!-----------------------------------------------------------------------
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
!  Compute the advection term of the water vapor specific humidity
!  equation and store the result in array qadv.
!
!-----------------------------------------------------------------------
!

  CALL set_acct(advs_acct)

  IF(sadvopt == 4 .or. tintegopt == 2 .or. tintegopt == 3) THEN ! Forward-based FCT scheme
    deltat = dtbig1
    tstrtlvl = tpresent
  ELSE
    deltat = dtbig1*2
    tstrtlvl = tpast
  END IF

  ! If sadvopt == 4, but we are on one of the first 2 RK3 steps, we need to temporarily set sadvopt to 3
  IF(sadvopt == 4 .and. (RK3step == 1 .or. RK3step == 2)) THEN
    sadvopt = 3
    RK3sadvflag = .true.
  ELSE
    RK3sadvflag = .false.
  END IF

  IF (sadvopt == 1 .OR. sadvopt == 2 .OR. sadvopt == 3) THEN
                            ! 2nd or 4th order advection

    CALL advq(nx,ny,nz,qv,u,v,wcont,ustr,vstr,wstr,                     &
              rhostr,qvbar, mapfct,                                     &
              qadv,                                                     &
              tem1,tem2,tem3,tem4)

  ELSE IF( sadvopt == 4.OR.sadvopt == 5 ) THEN  ! FCT advection

    CALL advqfct(nx,ny,nz,dtbig1,qv,u,v,wcont, ustr,vstr,wstr,          &
                 rhostr,rhostri,mapfct,j3,qadv,                         &
                 tem1,tem2,tem3,tem4,tem5,tem6,tem7,tem8,               &
                 tem1_0,tem2_0,tem3_0)

  END IF

  ! Reset sadvopt back to 4 if needed
  IF(RK3sadvflag) THEN
    sadvopt = 4
  END IF

  IF(tintegopt == 1 .or. tintegopt == 3 .or. (tintegopt == 2 .and. RK3step == 1)) THEN

  CALL set_acct(cmix_acct)
!
!-----------------------------------------------------------------------
!
!  Compute the mixing terms for the water vapor specific humidity
!  equation. This includes both physical and computational
!  mixing.  Store the result in array qmix.
!
!-----------------------------------------------------------------------
!
  CALL mixqv(nx,ny,nz,                                                  &
             qv(1,1,1,tstrtlvl),rhostr,qvbar,kmh,kmv,rprntl,qvsflx,     &
             pbldpth(1,1,tpresent),                                     &
             x,y,z,zp, mapfct, j1,j2,j3,aj3x,aj3y,j3inv,                &
             usflx,vsflx,ptsflx,ptsfc,qvsfc,ptbar,ptprt,                &
             qmix,                                                      &
             tem1,tem2,tem3,tem4,tem5,tem6)
!
!-----------------------------------------------------------------------
!
!  Call BRLXQ to calculate the boundary relaxation and computation
!  mixing for qv
!
!  qmix = qmix + qvexbc_term
!
!-----------------------------------------------------------------------
!
  IF (  lbcopt == 2 .AND. mgrid == 1 ) THEN

    CALL acct_interrupt(bc_acct)
    qflag = 0
    CALL brlxq(nx,ny,nz, deltat*0.5, qflag, qv,rhostr, qmix,            &
               exbcbuf(nqv0exb), exbcbuf(nqscalar0exb(qflag)),          &
               exbcbuf(nqvdtexb),exbcbuf(nqscalardtexb(qflag)),         &
               bcrlx,tem1,tem2,tem3,tem4)

    CALL acct_stop_inter

  END IF

  IF ( lbcopt == 1 .AND. mgrid == 1 ) THEN
    tem1(:,:,:) = qv(:,:,:,1) - qvbar(:,:,:)
    CALL brlxs_rbc( nx,ny,nz, tem1 ,rhostr, qmix, bcrlx )
  END IF

  END IF ! Calculating mixing
!
!-----------------------------------------------------------------------
!
!  Calculate qvforce, store the result in tem1.
!  Note here that qmix is either calculated above (for tintegopt 1 and 3)
!  or is only calculated once at the first RK3 step and passed back in
!  as input for tintegopt = 2.
!-----------------------------------------------------------------------
!

  DO k=1,nz-1
    DO j=1,ny-1
      DO i=1,nx-1
        tem1(i,j,k)=-qadv(i,j,k)+rhostr(i,j,k)*qvcumsrc(i,j,k)          &
                     +qmix(i,j,k)
        ! Store the forcing (in tem8) before vertically-implicit treatment so that we
        ! can retrieve it later
        tem8(i,j,k)=tem1(i,j,k)
      END DO
    END DO
  END DO
!
!-----------------------------------------------------------------------
!
!  Treat the vertically implicit mixing term.
!
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
          tem2(i,j,k)=rhostr(i,j,k)*kmv(i,j,k)*rprntl(i,j,k)            &
                      *j3inv(i,j,k)*j3inv(i,j,k)
        END DO
      END DO
    END DO

    CALL vmiximps(nx,ny,nz,deltat*0.5,qv(1,1,1,tstrtlvl),rhostr,        &
                  tem2,tem1,tem3,tem4,tem5,tem6)

    ! DTD: Retrieve the completed mixing term and store in tem8 for later dumping

    DO k=1,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          tem8(i,j,k) = tem1(i,j,k)-tem8(i,j,k)   ! tem8 now contains the vertically implicit mixing term
          tem8(i,j,k) = tem8(i,j,k)+qmix(i,j,k)   ! tem8 now contains the total mixing
        END DO
      END DO
    END DO

  ELSE ! No vertically-implicit mixing

    ! DTD: Retrieve the completed mixing term and store in tem8 for later dumping

    DO k=1,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          tem8(i,j,k) = qmix(i,j,k) ! tem8 now contains the total mixing term
        END DO
      END DO
    END DO

  END IF

  ! Multiply total mixing by rhostri to get actual total mixing for later dumping and diagnostics
  DO k=1,nz-1
    DO j=1,ny-1
      DO i=1,nx-1
        tem8(i,j,k)=tem8(i,j,k)*rhostri(i,j,k)
      END DO
    END DO
  END DO

!
!-----------------------------------------------------------------------
!
!  Integrate forward by one timestep the water vapor specific humidity
!  equation, yielding qv at time = tfuture.
!
!-----------------------------------------------------------------------
!
  IF( nestgrd == 1 .AND. mgrid /= 1 ) THEN

    CALL set_acct(tinteg_acct)

    DO k=2,nz-2
      DO j=2,ny-2
        DO i=2,nx-2
          qv(i,j,k,tfuture)=qv(i,j,k,tstrtlvl)                          &
                            +deltat*tem1(i,j,k)*rhostri(i,j,k)
        END DO
      END DO
    END DO
!
!-----------------------------------------------------------------------
!
!  Set the boundary conditions on qv for an adaptive (nested)
!  grid run.  If using only one grid, this portion of the code is
!  skipped....proceed to next comment block.
!
!-----------------------------------------------------------------------
!

    DO k=2,nz-2
      DO i=1,nx-1
        qv(i,   1,k,tfuture)=                                           &
            2*qv(i,   1,k,tpresent)-qv(i,   1,k,tpast)
        qv(i,ny-1,k,tfuture)=                                           &
            2*qv(i,ny-1,k,tpresent)-qv(i,ny-1,k,tpast)
      END DO
    END DO

    DO k=2,nz-2
      DO j=2,ny-2
        qv(   1,j,k,tfuture)=                                           &
            2*qv(   1,j,k,tpresent)-qv(   1,j,k,tpast)
        qv(nx-1,j,k,tfuture)=                                           &
            2*qv(nx-1,j,k,tpresent)-qv(nx-1,j,k,tpast)
      END DO
    END DO

! bc's for mp not implemented for nested grids
    CALL acct_interrupt(bc_acct)
    CALL bcqv(nx,ny,nz, deltat*0.5,                                     &
              qv,qvbar, qvdteb,qvdtwb,qvdtnb,qvdtsb,                    &
              0,0,0,0,tbc,bbc,                                          &
              ebc_global,wbc_global,nbc_global,sbc_global)
    CALL acct_stop_inter

  ELSE

    DO k=2,nz-2
      DO j=1,ny-1
        DO i=1,nx-1
          qv(i,j,k,tfuture)=qv(i,j,k,tstrtlvl)                          &
                            +deltat*tem1(i,j,k)*rhostri(i,j,k)
        END DO
      END DO
    END DO
!
!-----------------------------------------------------------------------
!
!  Set the boundary conditions on qv for a NON-adaptive (uniform)
!  grid run.
!
!-----------------------------------------------------------------------
!
    CALL set_acct(bc_acct)

    IF ( lbcopt == 1 ) THEN

      ebc1=ebc
      wbc1=wbc
      sbc1=sbc
      nbc1=nbc

      IF( ebc == 4 )  ebc1=0
      IF( wbc == 4 )  wbc1=0
      IF( nbc == 4 )  nbc1=0
      IF( sbc == 4 )  sbc1=0

      IF (mp_opt > 0) THEN
        CALL acct_interrupt(mp_acct)
        CALL mpsendrecv2dew(qv(1,1,1,tfuture),nx,ny,nz,ebc1,wbc1,0,tem2)
        CALL mpsendrecv2dns(qv(1,1,1,tfuture),nx,ny,nz,nbc1,sbc1,0,tem2)
      END IF
      CALL acct_interrupt(mp_acct)
      CALL bcqv(nx,ny,nz, dtbig1,                                       &
                qv,qvbar, qvdteb,qvdtwb,qvdtnb,qvdtsb,                  &
                ebc1,wbc1,nbc1,sbc1,tbc,bbc,                            &
                ebc_global,wbc_global,nbc_global,sbc_global)
      CALL acct_stop_inter

    ELSE

      ebc1 = 3
      wbc1 = 3
      sbc1 = 3
      nbc1 = 3
      IF( ebc == 0 )  ebc1 = 0
      IF( wbc == 0 )  wbc1 = 0
      IF( nbc == 0 )  nbc1 = 0
      IF( sbc == 0 )  sbc1 = 0

      IF (mp_opt > 0) THEN
        CALL acct_interrupt(mp_acct)
        CALL mpsendrecv2dew(qv(1,1,1,tfuture),nx,ny,nz,ebc1,wbc1,0,tem2)
        CALL mpsendrecv2dns(qv(1,1,1,tfuture),nx,ny,nz,nbc1,sbc1,0,tem2)
      END IF
      CALL acct_interrupt(mp_acct)
      CALL bcqv(nx,ny,nz, dtbig1,                                       &
                qv,qvbar, qvdteb,qvdtwb,qvdtnb,qvdtsb,                  &
                ebc1,wbc1,nbc1,sbc1,tbc,bbc,                            &
                ebc_global,wbc_global,nbc_global,sbc_global)

      qflag = 0
      CALL exbcq(nx,ny,nz, qflag, curtim+dtbig,qv(1,1,1,tfuture),       &
                 exbcbuf(nqv0exb), exbcbuf(nqvdtexb))
      CALL acct_stop_inter

    END IF

  END IF

  RETURN
END SUBROUTINE solvqv
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE SOLVQ                      ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE solvq(RK3step,nx,ny,nz, exbcbufsz, dtbig1, qflag,            &
         q,u,v,wcont, ustr,vstr,wstr,                                   &
         qdteb,qdtwb,qdtnb,qdtsb,                                       &
         rhostr,rhostri,kmh,kmv,rprntl,                                 &
         x,y,z,zp, mapfct, j1,j2,j3,aj3x,aj3y,j3inv,                    &
         qccumsrc,qrcumsrc,qicumsrc,qscumsrc,                           &
         exbcbuf,bcrlx,                                                 &
         qadv,qmix,                                                     &
         tem1,tem2,tem3,tem4,tem5,tem6,tem7,tem8,                       &
         tem1_0,tem2_0,tem3_0)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Coordinate the time integration of the equations for the water
!  substance quantities qc,qr,qi,qs and qh, i.e. all water variables
!  except the water vapor specific humidity qv.
!
!-----------------------------------------------------------------------
!
!
!  AUTHOR: Ming Xue
!  10/10/91
!
!  MODIFICATION HISTORY:
!
!  5/20/92 (K. Droegemeier and M. Xue)
!  Added full documentation.
!
!  2/10/93 (K. Droegemeier)
!  Cleaned up documentation.
!
!  8/12/95 (M. Xue)
!  Added flag qflag.
!
!  1/25/96 (Donghai Wang & Yuhe Liu)
!  Added the map projection factor to ARPS governing equations.
!
!  4/1/96 (Donghai Wang, X. Song and M. Xue)
!  Added the implicit treatment for the vertical mixing.
!
!  3/24/1997 (M. Xue)
!  Code to handle the case of forward time integration added.
!
!  7/10/1997 (Fanyou Kong -- CMRP)
!  Added the positive definite advection scheme option (sadvopt = 5)
!
!  4/15/1998 (Donghai Wang)
!  Added the source terms to the right hand terms of the qc,qr,qi,qs
!  equations due to the K-F cumulus parameterization.
!
!  7/17/1998 (M. Xue)
!  Changed call to ADVQFCT.
!
!  9/18/1998 (D. Weber)
!  Added arrays aj3x,y, u,v,wstr,rhostri, do not alter!
!
!  10/19/2009 (D. Dawson)
!  Updated to accomodate RK3 time integration
!
!  04/25/2012 (Y. Wang)
!  Added polluant concentration integration.
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    RK3step  Which RK3 step are we on (1,2,3)?
!             (only used for tintegopt == 2)
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!    dtbig1   The big time step size (s)
!    qflag    Indicator for water/ice species
!
!    q        One of the liquid or ice variables at tpast and
!             tpresent (kg/kg)
!
!    u        x component of velocity at all time levels (m/s)
!    v        y component of velocity at all time levels (m/s)
!    wcont    Contravariant vertical velocity (m/s)
!    ustr     Appropriate ustr for advection.
!    vstr     Appropriate vstr for advection.
!    wstr     Appropriate wstr for advection.
!
!    qdteb    Time tendency of liquid/ice variable q at east boundary
!    qdtwb    Time tendency of liquid/ice variable q at west boundary
!    qdtnb    Time tendency of liquid/ice variable q at north
!             boundary
!    qdtsb    Time tendency of liquid/ice variable q at south
!             boundary
!
!    rhostr   Base state density rhobar times j3 (kg/m**3)
!    rhostri  Inverse base state density rhobar times j3 (kg/m**3)
!
!    kmh      Horizontal turb. mixing coef. for momentum ( m**2/s )
!    kmv      Vertical turb. mixing coef. for momentum ( m**2/s )
!    rprntl   Reciprocal of Prandtl number
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
!    qccumsrc Source term in qc-equation due to cumulus parameterization
!    qrcumsrc Source term in qr-equation due to cumulus parameterization
!    qicumsrc Source term in qi-equation due to cumulus parameterization
!    qscumsrc Source term in qs-equation due to cumulus parameterization
!
!
!  OUTPUT:
!
!    q        Liquid/ice variable q at tfuture (kg/kg)
!
!  WORK ARRAYS:
!
!    qadv     Advection term of water variable eq. (kg/(m**3*s)).
!             A local array.
!    qmix     Total mixing in water variable eq. (kg/(m**3*s)).
!             A local array.
!    tem1     Temporary work array.
!    tem2     Temporary work array.
!    tem3     Temporary work array.
!    tem4     Temporary work array.
!    tem5     Temporary work array.
!    tem6     Temporary work array.
!    tem7     Temporary work array.
!    tem8     Temporary work array.
!
!    tem1_0   Temporary work array.
!    tem2_0   Temporary work array.
!    tem3_0   Temporary work array.
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

  INTEGER :: RK3step

  INTEGER :: nx, ny, nz        ! Number of grid points in 3 directions
  REAL    :: dtbig1            ! Local big time step
  INTEGER :: qflag             ! Indicator for water/ice species

  REAL :: q     (nx,ny,nz,nt)  ! One of the water/ice variables (kg/kg)
  REAL :: u     (nx,ny,nz,nt)  ! Total u-velocity (m/s)
  REAL :: v     (nx,ny,nz,nt)  ! Total v-velocity (m/s)
  REAL :: wcont (nx,ny,nz)     ! Contravariant vertical velocity (m/s)
  REAL :: ustr  (nx,ny,nz)     ! Appropriate ustr for advection (m/s)
  REAL :: vstr  (nx,ny,nz)     ! Appropriate vstr for advection (m/s)
  REAL :: wstr  (nx,ny,nz)     ! Appropriate wstr for advection (m/s)

  REAL :: qdteb(ny,nz)         ! Time tendency of liquid/ice at
                               ! e-boundary
  REAL :: qdtwb(ny,nz)         ! Time tendency of liquid/ice at
                               ! w-boundary
  REAL :: qdtnb(nx,nz)         ! Time tendency of liquid/ice at
                               ! n-boundary
  REAL :: qdtsb(nx,nz)         ! Time tendency of liquid/ice at
                               ! s-boundary

  REAL :: rhostr(nx,ny,nz)     ! Base state density rhobar times j3.
  REAL :: rhostri(nx,ny,nz)    ! Inverse base state density rhobar times j3.

  REAL :: kmh   (nx,ny,nz)     ! Horizontal turb. mixing coef. for
                               ! momentum. ( m**2/s )
  REAL :: kmv   (nx,ny,nz)     ! Vertical turb. mixing coef. for
                               ! momentum. ( m**2/s )
  REAL :: rprntl(nx,ny,nz)     ! Reciprocal of Prandtl number

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
  REAL :: qccumsrc(nx,ny,nz)   ! Source in qc-equation due to cumulus
                               ! parameterization
  REAL :: qrcumsrc(nx,ny,nz)   ! Source in qr-equation due to cumulus
                               ! parameterization
  REAL :: qicumsrc(nx,ny,nz)   ! Source in qi-equation due to cumulus
                               ! parameterization
  REAL :: qscumsrc(nx,ny,nz)   ! Source in qs-equation due to cumulus
                               ! parameterization

  INTEGER :: exbcbufsz         ! EXBC buffer size
  REAL :: exbcbuf( exbcbufsz ) ! EXBC buffer array
  REAL :: bcrlx (nx,ny)        ! EXBC relaxation coefficients

  REAL :: qadv(nx,ny,nz)       ! Advection of water/ice substance
                               ! (kg/(m**3*s)). A local array.
  REAL :: qmix(nx,ny,nz)       ! Total mixing of water/ice substance
                               ! (kg/(m**3*s)). A local array.
  REAL :: tem1  (nx,ny,nz)     ! Temporary work array
  REAL :: tem2  (nx,ny,nz)     ! Temporary work array
  REAL :: tem3  (nx,ny,nz)     ! Temporary work array
  REAL :: tem4  (nx,ny,nz)     ! Temporary work array
  REAL :: tem5  (nx,ny,nz)     ! Temporary work array
  REAL :: tem6  (nx,ny,nz)     ! Temporary work array
  REAL :: tem7  (nx,ny,nz)     ! Temporary work array
  REAL :: tem8  (nx,ny,nz)     ! Temporary work array

  REAL :: tem1_0(0:nx,0:ny,0:nz) ! automatic work array
  REAL :: tem2_0(0:nx,0:ny,0:nz) ! automatic work array
  REAL :: tem3_0(0:nx,0:ny,0:nz) ! automatic work array
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k
  INTEGER :: tstrtlvl
  REAL    :: deltat
  INTEGER :: ebc1,wbc1,nbc1,sbc1
  INTEGER :: nq0exb,nqdtexb
  LOGICAL :: RK3sadvflag      ! Flag to indicate if sadvopt == 4 and RK3
                              ! time-stepping is being used.  It's used
                              ! to allow the normal advection to be called
                              ! for the first 2 RK3 steps, and the full
                              ! FCT only for the last one.
  INTEGER :: fzone,  ic, jc, ica, jca, l
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'bndry.inc'
  INCLUDE 'exbc.inc'
  INCLUDE 'globcst.inc'
  INCLUDE 'mp.inc'            ! Message passing parameters.

  INCLUDE 'grid.inc'          ! to use dxinv, dyinv, dzinv for cc
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
!  Compute the advection term for a general water substance
!  q and store the result in qadv.
!
!-----------------------------------------------------------------------
!
  CALL set_acct(advs_acct)

  IF(sadvopt == 4 .or. tintegopt == 2 .or. tintegopt == 3) THEN ! Forward-based FCT scheme
    deltat = dtbig1
    tstrtlvl = tpresent
  ELSE
    deltat = dtbig1*2
    tstrtlvl = tpast
  END IF

  ! If sadvopt == 4, but we are on one of the first 2 RK3 steps, we need to temporarily set sadvopt to 3
  IF(sadvopt == 4 .and. (RK3step == 1 .or. RK3step == 2)) THEN
    sadvopt = 3
    RK3sadvflag = .true.
  ELSE
    RK3sadvflag = .false.
  END IF

  IF (sadvopt == 1 .OR. sadvopt == 2 .OR. sadvopt == 3) THEN
                                        ! 2nd or 4th order centerd sc

    tem7(:,:,:) = 0.0

    CALL advq(nx,ny,nz,q,u,v,wcont,ustr,vstr,wstr,                      &
              rhostr, tem7, mapfct,                                     &
              qadv,                                                     &
              tem1,tem2,tem3,tem4)

  ELSE IF( sadvopt == 4 .OR. sadvopt == 5) THEN  ! FCT advection

    CALL advqfct(nx,ny,nz,dtbig1,q,u,v,wcont, ustr,vstr,wstr,           &
                 rhostr,rhostri,mapfct,j3,qadv,                         &
                 tem1,tem2,tem3,tem4,tem5,tem6,tem7,tem8,               &
                 tem1_0,tem2_0,tem3_0)

  END IF

  ! Reset sadvopt back to 4 if needed
  IF(RK3sadvflag) THEN
    sadvopt = 4
  END IF

  IF(tintegopt == 1 .or. tintegopt == 3 .or. (tintegopt == 2 .and. RK3step == 1)) THEN
!
!-----------------------------------------------------------------------
!
!  Compute the mixing terms for the general water substance (q)
!  equation, including both physical and computational mixing.
!  Store the result in array qmix.
!
!-----------------------------------------------------------------------
!
  CALL set_acct(cmix_acct)
  CALL mixq(nx,ny,nz,                                                   &
            q(1,1,1,tstrtlvl),rhostr,kmh,kmv,rprntl,                    &
            x,y,z,zp,mapfct, j1,j2,j3,aj3x,aj3y,j3inv,                  &
            qmix,                                                       &
            tem1,tem2,tem3,tem4,tem5,tem6)
!
!-----------------------------------------------------------------------
!
!  Call BRLXQ to added to qmix the additional boundary relaxation and
!  spatial mixing in the boundary zone
!
!  qmix = qmix + qexbc_term
!
!-----------------------------------------------------------------------
!
  IF (  lbcopt == 2 .AND. mgrid == 1 ) THEN

    CALL acct_interrupt(bc_acct)

    CALL brlxq(nx,ny,nz, deltat*0.5,qflag, q,rhostr, qmix,              &
               exbcbuf(nqv0exb), exbcbuf(nqscalar0exb(qflag)),          &
               exbcbuf(nqvdtexb),exbcbuf(nqscalardtexb(qflag)),         &
               bcrlx,tem1,tem2,tem3,tem4)

    CALL acct_stop_inter

  END IF

  IF ( lbcopt == 1 .AND. mgrid == 1 ) THEN
    CALL brlxs_rbc( nx,ny,nz, q,rhostr, qmix, bcrlx )
  END IF

  END IF ! Calculating mixing terms
!
!-----------------------------------------------------------------------
!
!  Calculate qforce, store the result in tem1.
!
!-----------------------------------------------------------------------
!
  IF (qflag == P_CC) THEN
!michi
!
    tem7 = 0.0
    fzone = 3

    DO l=1,cpoint

      if(curtim.ge.ccstart(l).and.curtim.le.ccend(l)) then

        do jc=1,nproc_y
          do ic=1,nproc_x

            if(myproc == (ic+(jc-1)*nproc_x-1)) then

              do k = 1, nz-1
                do j = 1, ny-1
                  do i = 1, nx-1

                  ica = i + (ic-1)*(nx-fzone)
                  jca = j + (jc-1)*(ny-fzone)

                  if(ica.eq.icc(l) .and. jca.eq.jcc(l) .and. k.eq.kcc(l)) then
                !TINA divide by dx,dy,dz - I think this is needed for consistency
                !        tem9(i,j,k) = ccemit(l)
                    tem7(i,j,k) = ccemit(l)*dxinv*dyinv*dzinv
                !  write(*,*) ' emitted point '
                !  write(*,*) ' ica,jca,k = ',ica,jca,k
                !  write(*,*) ' i,j,k     = ',i,j,k
                !  write(*,*) ' myproc = ',myproc
                  endif
                !
                  enddo
                enddo
              enddo

            endif

          enddo
        enddo

      ENDIF

    END DO
!
!michi
!
!TINA - testing 7/6/07
!mass conservation - sum up total mass emitted by source
    !summass = 0.D0
    !DO k=2,nz-2
    !  DO j=2,ny-2
    !    DO i=2,nx-2
    !      !add up mass emited/s/m3 times each cell volume times deltat
    !      summass = summass + tem9(i,j,k)/(dxinv*dyinv*dzinv)*deltat
    !    END DO
    !  END DO
    !END DO
!      print *, 'total emitted = ', summass, ' kg, at time = ', curtim
!      write(999,*) curtim, summass
!TINA added special option above
  END IF

  DO k=1,nz-1
    DO j=1,ny-1
      DO i=1,nx-1
        tem1(i,j,k)=-qadv(i,j,k)+qmix(i,j,k)
        IF(qflag == P_QC) THEN
          tem1(i,j,k)= tem1(i,j,k)+rhostr(i,j,k)*qccumsrc(i,j,k)
        ELSE IF(qflag == P_QR) THEN
          tem1(i,j,k)= tem1(i,j,k)+rhostr(i,j,k)*qrcumsrc(i,j,k)
        ELSE IF(qflag == P_QI) THEN
          tem1(i,j,k)= tem1(i,j,k)+rhostr(i,j,k)*qicumsrc(i,j,k)
        ELSE IF(qflag == P_QS) THEN
          tem1(i,j,k)= tem1(i,j,k)+rhostr(i,j,k)*qscumsrc(i,j,k)
        ELSE IF (qflag == P_CC) THEN
          tem1(i,j,k) = tem1(i,j,k)+rhostr(i,j,k)*tem7(i,j,k)           &
                      - q(i,j,k,tpresent)*0.0005*rhostr(i,j,k)
        ELSE
        END IF
        ! DTD: make a copy of tem1 in tem8 so that we can retrieve
        ! the vertically-implicit mixing term later
        tem8(i,j,k) = tem1(i,j,k)
      END DO
    END DO
  END DO

!
!-----------------------------------------------------------------------
!
!  Treat the vertically implicit mixing term.
!
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
          tem2(i,j,k)=rhostr(i,j,k)*kmv(i,j,k)*rprntl(i,j,k)            &
                      *j3inv(i,j,k)*j3inv(i,j,k)
        END DO
      END DO
    END DO

    CALL vmiximps(nx,ny,nz,deltat*0.5,q(1,1,1,tstrtlvl),rhostr,         &
                  tem2,tem1,tem3,tem4,tem5,tem6)

    ! DTD: Retrieve the completed mixing term and store in tem8 for later dumping

    DO k=1,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          tem8(i,j,k) = tem1(i,j,k)-tem8(i,j,k)   ! tem8 now contains the vertically implicit mixing term
          tem8(i,j,k) = tem8(i,j,k)+qmix(i,j,k)   ! tem8 now contains the total mixing term
        END DO
      END DO
    END DO

  ELSE ! No vertically-implicit mixing

    ! DTD: Retrieve the completed mixing term and store in tem8 for later dumping

    DO k=1,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          tem8(i,j,k) = qmix(i,j,k) ! tem8 now contains the total mixing term
        END DO
      END DO
    END DO

  END IF

  ! Multiply total mixing by rhostri to get actual total mixing for later dumping and diagnostics
  DO k=1,nz-1
    DO j=1,ny-1
      DO i=1,nx-1
        tem8(i,j,k)=tem8(i,j,k)*rhostri(i,j,k)
      END DO
    END DO
  END DO

!
!-----------------------------------------------------------------------
!
!  Integrate forward by one timestep the general water substance
!  (q) equation, yielding q at time = tfuture.
!
!-----------------------------------------------------------------------
!
  CALL set_acct(tinteg_acct)

  IF( nestgrd == 1 .AND. mgrid /= 1 ) THEN

    DO k=2,nz-2
      DO j=2,ny-2
        DO i=2,nx-2
          q(i,j,k,tfuture)=q(i,j,k,tstrtlvl)                            &
                           +deltat*tem1(i,j,k)*rhostri(i,j,k)
        END DO
      END DO
    END DO
!
!-----------------------------------------------------------------------
!
!  Set the boundary conditions on q for an adaptive (nested)
!  grid run.  If using only one grid, this portion of the code is
!  skipped....proceed to next comment block.
!
!-----------------------------------------------------------------------
!
    DO k=2,nz-2
      DO i=1,nx-1
        q(i,   1,k,tfuture)=2*q(i,   1,k,tpresent)-q(i,   1,k,tpast)
        q(i,ny-1,k,tfuture)=2*q(i,ny-1,k,tpresent)-q(i,ny-1,k,tpast)
      END DO
    END DO

    DO k=2,nz-2
      DO j=2,ny-2
        q(   1,j,k,tfuture)=2*q(   1,j,k,tpresent)-q(   1,j,k,tpast)
        q(nx-1,j,k,tfuture)=2*q(nx-1,j,k,tpresent)-q(nx-1,j,k,tpast)
      END DO
    END DO

! bc's for mp not implemented for nested grids
    CALL acct_interrupt(bc_acct)
    CALL bcq(nx,ny,nz, dtbig1, q, qdteb,qdtwb,qdtnb,qdtsb,              &
             0,0,0,0,tbc,bbc,                                           &
             ebc_global,wbc_global,nbc_global,sbc_global)
    CALL acct_stop_inter

  ELSE

    DO k=2,nz-2
      DO j=1,ny-1
        DO i=1,nx-1
          q(i,j,k,tfuture)=q(i,j,k,tstrtlvl)                            &
                          +deltat*tem1(i,j,k)*rhostri(i,j,k)
        END DO
      END DO
    END DO
!
!-----------------------------------------------------------------------
!
!  Set the boundary conditions on q for a NON-adaptive (uniform)
!  grid run.
!
!-----------------------------------------------------------------------
!
    IF ( lbcopt == 1 ) THEN

      ebc1=ebc
      wbc1=wbc
      sbc1=sbc
      nbc1=nbc

      IF( ebc == 4 )  ebc1=0
      IF( wbc == 4 )  wbc1=0
      IF( nbc == 4 )  nbc1=0
      IF( sbc == 4 )  sbc1=0

      IF (mp_opt > 0) THEN
        CALL acct_interrupt(mp_acct)
        CALL mpsendrecv2dew(q(1,1,1,tfuture),nx,ny,nz,ebc1,wbc1,0,tem2)
        CALL mpsendrecv2dns(q(1,1,1,tfuture),nx,ny,nz,nbc1,sbc1,0,tem2)
      END IF
      CALL acct_interrupt(bc_acct)
      CALL bcq(nx,ny,nz, dtbig1, q, qdteb,qdtwb,qdtnb,qdtsb,            &
               ebc1,wbc1,nbc1,sbc1,tbc,bbc,                             &
               ebc_global,wbc_global,nbc_global,sbc_global)
      CALL acct_stop_inter

    ELSE

      ebc1 = 3
      wbc1 = 3
      sbc1 = 3
      nbc1 = 3
      IF( ebc == 0 )  ebc1 = 0
      IF( wbc == 0 )  wbc1 = 0
      IF( nbc == 0 )  nbc1 = 0
      IF( sbc == 0 )  sbc1 = 0

      IF (mp_opt > 0) THEN
        CALL acct_interrupt(mp_acct)
        CALL mpsendrecv2dew(q(1,1,1,tfuture),nx,ny,nz,ebc1,wbc1,0,tem2)
        CALL mpsendrecv2dns(q(1,1,1,tfuture),nx,ny,nz,nbc1,sbc1,0,tem2)
      END IF
      CALL acct_interrupt(bc_acct)
      CALL bcq(nx,ny,nz, dtbig1, q, qdteb,qdtwb,qdtnb,qdtsb,            &
               ebc1,wbc1,nbc1,sbc1,tbc,bbc,                             &
               ebc_global,wbc_global,nbc_global,sbc_global)
          ! Zero-gradient condition will be
          ! reset if exbc data is available for q

!      IF ( qflag == 1 ) THEN
!        nq0exb = nqv0exb
!        nqdtexb = nqvdtexb
!      ELSE IF ( qflag == 2 ) THEN
!        nq0exb = nqc0exb
!        nqdtexb = nqcdtexb
!      ELSE IF ( qflag == 3 ) THEN
!        nq0exb = nqr0exb
!        nqdtexb = nqrdtexb
!      ELSE IF ( qflag == 4 ) THEN
!        nq0exb = nqi0exb
!        nqdtexb = nqidtexb
!      ELSE IF ( qflag == 5 ) THEN
!        nq0exb = nqs0exb
!        nqdtexb = nqsdtexb
!      ELSE IF ( qflag == 6 ) THEN
!        nq0exb = nqh0exb
!        nqdtexb = nqhdtexb
!      END IF

      nq0exb  = nqscalar0exb(qflag)
      nqdtexb = nqscalar0exb(qflag)

      CALL exbcq(nx,ny,nz, qflag, curtim+dtbig,q(1,1,1,tfuture),        &
                 exbcbuf(nq0exb),exbcbuf(nqdtexb))
      CALL acct_stop_inter

    END IF

  END IF

  RETURN
END SUBROUTINE solvq
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                  SUBROUTINE UPRAD1                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE uprad1(nx,ny,w1,w2,wc,nrho,ifax1,ifax2,trigs1,trigs2,        &
           temxy1,work)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  To apply the upper radiation condition using a periodic 1 or 2-D
!  fft.  The variable temxy1 contains input information and is
!  overwritten with the final w (at nz-1) on exit.
!
!  In the case of the open boundary, a linear hydrostatic relation
!  between vertical velocity and pressure is used in fourier space.
!  The use of fast fourier transform is required for this condition.
!
!  The fft routine used is FFT991 obtained from NCAR. It is a
!  version of the one used by ECMWF. Before fft991 is called, set99
!  is called to setup the variables needed for the transform routine.
!  (init3d.f)
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Dan Weber
!  07/22/1997.
!
!  Modification History:
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx       Number of grid points in the x-direction.
!    ny       Number of grid points in the y-direction.
!
!    w1       Pre-fft w(nz-1) coefficient.
!    w2       Pre-fft w(nz-2) coefficient.
!    wc       Pre-fft w(nz-1) coefficient from the reduction
!             phase of the tridiaginal solver.
!
!    trigs1   Array containing pre-computed trig function for
!             fftopt=1.
!    trigs2   Array containing pre-computed trig function for
!             fftopt=1.
!    ifax1    Array containing the factors of nx for fftopt=1.
!    ifax2    Array containing the factors of ny for fftopt=1.
!
!  OUTPUT:
!
!
!    temxy1   On input known forcing term, on output final w at nz-1.
!
!
!  WORK ARRAYS:
!
!
!    work     2-D work array for fft code.
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

  INTEGER :: nx                ! Number of grid points in the x-direction.
  INTEGER :: ny                ! Number of grid points in the y-direction.

  REAL :: w1                   ! pre-fft w(nz-1) coefficient.
  REAL :: w2                   ! pre-fft w(nz-2) coefficient.
  REAL :: wc                   ! pre-fft w(nz-1) coefficient from the
                               ! reduction phase of the tridiagonal solver.

  REAL :: trigs1(3*(nx-1)/2+1) ! Array containing pre-computed trig
                               ! function for fftopt=1.
  REAL :: trigs2(3*(ny-1)/2+1) ! Array containing pre-computed trig
                               ! function for fftopt=1.
  INTEGER :: ifax1(13)         ! Array containing the factors of nx for
                               ! fftopt=1.
  INTEGER :: ifax2(13)         ! Array containing the factors of ny for
                               ! fftopt=1.


  REAL :: temxy1 (nx+1,ny+1)   ! On input known forcing term,
                               ! on output final w at nz-1.


  REAL :: work  (ny,nx)        ! 2-D work array for fft code.

!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!

  INTEGER :: i,j,itema,itemb
  REAL :: eps,kx,ky,pi,nrho,tema,temb,temc,temd,teme,temf,temg

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

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!-----------------------------------------------------------------------
!
!  Perform the fft on temxy1
!
!-----------------------------------------------------------------------

  IF(ny == 4)THEN  ! 1-d transform in x, note nx must be ODD!

    itema = 1
    itemb = nx

    CALL fft991(temxy1,work,trigs1,ifax1,1,nx+1,nx-1,1,-1)

  ELSE IF(nx == 4)THEN ! 1-d transform in y, note ny must be ODD!

    itema = ny
    itemb = 1

    CALL fft991(temxy1,work,trigs2,ifax2,nx+1,1,ny-1,1,-1)

  ELSE       ! Do 2-d transform, note nx,ny must be ODD!

    itema = ny
    itemb = nx

    CALL fft991(temxy1,work,trigs1,ifax1,1,nx+1,nx-1,ny-1,-1)
    CALL fft991(temxy1,work,trigs2,ifax2,nx+1,1,ny-1,nx,-1)

  END IF   !  end of run type if block.

!-----------------------------------------------------------------------
!
!  Compute the wave space w at nz-1.
!
!-----------------------------------------------------------------------

  eps= 1.0E-8
  pi   = 4.0*ATAN(1.0)
  tema = pi/(nx-2)
  temb = pi/(ny-2)
  temc = 2.0*dxinv
  temd = 2.0*dyinv
  temf = w1-w2*wc
  temg = eps - nrho

  DO j=1,itema ! note: kx and ky are in finite difference form.
    ky = temd*SIN(INT((j-1)*0.5+0.001)*temb)
    ky = ky*ky
    DO i=1,itemb
      kx = temc*SIN(INT((i-1)*0.5+0.001)*tema)
      teme =SQRT(kx*kx+ky)
      temxy1(i,j)=-teme*temxy1(i,j)/(teme*temf+temg)
    END DO
  END DO

!-----------------------------------------------------------------------
!
!  Perform the reverse fft on temxy1 to obtain w(nz-1) in real space.
!
!-----------------------------------------------------------------------

  IF(ny == 4)THEN  ! 1-d transform in x, note nx must be ODD!

    CALL fft991(temxy1,work,trigs1,ifax1,1,nx+1,nx-1,1,1)

    DO j=2,ny-1
      DO i=1,nx-1
        temxy1(i,j) = temxy1(i,1)
      END DO
    END DO

  ELSE IF(nx == 4)THEN ! 1-d transform in y, note ny must be ODD!

    CALL fft991(temxy1,work,trigs2,ifax2,nx+1,1,ny-1,1,1)

    DO j=1,ny-1
      DO i=2,nx-1
        temxy1(i,j) = temxy1(1,j)
      END DO
    END DO

  ELSE              ! Do 2-d transform, note nx,ny must be ODD!

    CALL fft991(temxy1,work,trigs2,ifax2,nx+1,1,ny-1,nx, 1)
    CALL fft991(temxy1,work,trigs1,ifax1,1,nx+1,nx-1,ny-1, 1)

  END IF    !  end of run type if block.

  RETURN
END SUBROUTINE uprad1

!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                  SUBROUTINE UPRAD3                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE uprad3(nx,ny,w1,w2,wc,nrho,wsave1,wsave2,vwork1,vwork2,      &
           temxy1,work)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  To apply the upper radiation condition using an even (vfftpack)
!  sequence fft.  The variable temxy1 contains input information and
!  is overwritten with the updated w (at nz-1) on exit.
!
!  In the case of the open boundary, a linear, hydrostatic relation
!  between vertical velocity to pressure is used in fourier space.
!  The use of fast fourier transform is NOT required for this
!  option but is recommended.  If nx and ny are not of special
!  character, then a basic fourier transform is used.
!
!  The fft routine used is vcost from the vfftpack software suite
!  available at PSC.
!
!  Before vcost is called, the arrays wsave1,wsave2,vwork1,vwork2
!  must be initialized (init3d.f).
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Dan Weber
!  07/22/1997.
!
!
!  Modification History:
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx       Number of grid points in the x-direction.
!    ny       Number of grid points in the y-direction.
!
!    w1       Pre-fft w(nz-1) coefficient.
!    w2       Pre-fft w(nz-2) coefficient.
!    wc       Pre-fft w(nz-1) coefficient from the reduction
!             phase of the tridiaginal solver.
!
!    vwork1   2-D work array for fftopt=2 option.
!    vwork2   2-D work array for fftopt=2 option.
!    wsave1   Work array for fftopt=2 option.
!    wsave2   Work array for fftopt=2 option.
!
!  OUTPUT:
!
!
!    temxy1   On input known forcing term, on output final w at nz-1.
!
!
!  WORK ARRAYS:
!
!
!    work     2-D work array for fft code.
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

  IMPLICIT NONE

  INTEGER :: nx                ! Number of grid points in the x
                               ! direction.
  INTEGER :: ny                ! Number of grid points in the y
                               ! direction.

  REAL :: w1                   ! Pre-fft w(nz-1) coefficient.
  REAL :: w2                   ! Pre-fft w(nz-2) coefficient.
  REAL :: wc                   ! Pre-fft w(nz-1) coefficient from the
                               ! reduction phase of the tridiagonal
                               ! solver.

  REAL :: vwork1 (nx+1,ny+1)   ! 2-D work array for fftopt=2.
  REAL :: vwork2 (ny,nx+1)     ! 2-D work array for fftopt=2.
  REAL :: wsave1 (3*(ny-1)+15) ! Work array for fftopt=2.
  REAL :: wsave2 (3*(nx-1)+15) ! Work array for fftopt=2.

  REAL :: temxy1 (nx+1,ny+1)   ! On input known forcing term,
                               ! on output final w at nz-1.

  REAL :: work   (ny,nx)       ! 2-D work array (global)

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
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j
  REAL :: kx,ky,pi,tema,temb,temc,temd,teme,temf
  REAL :: eps,nrho
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


!-----------------------------------------------------------------------
!
!  Perform the fft on temxy1.
!
!-----------------------------------------------------------------------

  IF(ny == 4)THEN  ! 1-d transform in x.

    DO i=1,nx-1    ! swap temxy1(i,j) array into work(j,i)
      work(1,i) = temxy1(i,1)
    END DO

    CALL vcost(1,nx-1,work,vwork2,ny,wsave2)

    DO i=1,nx-1    ! swap back for computations....
      temxy1(i,1) = work(1,i)
    END DO

  ELSE IF(nx == 4)THEN ! 1-d transform in y.

    CALL vcost(1,ny-1,temxy1,vwork1,nx+1,wsave1)

  ELSE                ! Do 2-d transform.

    CALL vcost(nx-1,ny-1,temxy1,vwork1,nx+1,wsave1)

    DO j=1,ny-1    ! swapping arrays.
      DO i=1,nx-1
        work(j,i) = temxy1(i,j)
      END DO
    END DO

    CALL vcost(ny-1,nx-1,work,vwork2,ny,wsave2)

    DO j=1,ny-1    ! swapping arrays for computations.
      DO i=1,nx-1
        temxy1(i,j) = work(j,i)
      END DO
    END DO

  END IF              ! end of run type if block.

!-----------------------------------------------------------------------
!
!  Compute the wave space w at nz-1.
!
!-----------------------------------------------------------------------

  eps  = 1.0E-13
  pi   = 4.0*ATAN(1.0)
  tema = 0.5*pi/(nx-2)   ! finite difference form.
  temb = 0.5*pi/(ny-2)   ! finite difference form.
  temc = 2.0*dxinv
  temd = 2.0*dyinv
  teme = w1-w2*wc
  temf = eps - nrho

  temxy1(1,1) = -temxy1(1,1)/(teme+temf)  ! first part of 2-d case.

  IF(ny == 4)THEN        ! perform computations in x direction.

    DO i=2,nx-1   ! compute from i=2,nx-1 and j=1...
      kx = temc*SIN(INT(i-1)*tema)
      kx = SQRT(kx*kx)
      temxy1(i,1)=-kx*temxy1(i,1)/(kx*teme+temf)
    END DO

  ELSE IF(nx == 4)THEN    ! perform computations in y direction.

    DO j=2,ny-1   ! compute from j=2,ny-1 and i=1.
      ky = temd*SIN(INT(j-1)*temb)
      ky = SQRT(ky*ky)
      temxy1(1,j)=-ky*temxy1(1,j)/(ky*teme+temf)
    END DO

  ELSE                   ! compute the 2-D version.

    DO i=2,nx-1   ! compute from i=2,nx-1 and j=1.
      kx = temc*SIN(INT(i-1)*tema)
      kx = SQRT(kx*kx)
      temxy1(i,1)=-kx*temxy1(i,1)/(kx*teme+temf)
    END DO

    DO j=2,ny-1   ! compute from i=1,j=2,ny-1.
      ky = temd*SIN(INT(j-1)*temb)
      ky = SQRT(ky*ky)
      temxy1(1,j)=-ky*temxy1(1,j)/(ky*teme+temf)
    END DO

    DO j=2,ny-1  ! compute from i=2,nx-1,j=2,ny-1.
      ky = temd*SIN(INT(j-1)*temb)
      ky = ky*ky
      DO i=2,nx-1
        kx = temc*SIN(INT(i-1)*tema)
        kx =SQRT(kx*kx+ky)
        temxy1(i,j)=-kx*temxy1(i,j)/(kx*teme+temf)
      END DO
    END DO

  END IF              ! end of the run type if block.

!-----------------------------------------------------------------------
!
!  Perform the reverse fft on temxy1 to obtain w in real space.
!
!-----------------------------------------------------------------------

  IF(ny == 4)THEN  ! 1-d transform in x, note nx must be EVEN!

    DO i=1,nx-1
      work(1,i) = temxy1(i,1)
    END DO

    CALL vcost(1,nx-1,work,vwork2,ny,wsave2)

    DO j=1,ny-1
      DO i=1,nx-1
        temxy1(i,j) = work(1,i)
      END DO
    END DO

  ELSE IF(nx == 4)THEN  ! 1-d transform in y, note ny must be EVEN!

    CALL vcost(1,ny-1,temxy1,vwork1,nx+1,wsave1)

    DO j=1,ny-1
      DO i=1,nx-1
        temxy1(i,j) = temxy1(1,j)
      END DO
    END DO

  ELSE                 ! Do 2-d transform, note nx,ny must be EVEN!

    DO j=1,ny-1    ! swapping arrays.
      DO i=1,nx-1
        work(j,i) = temxy1(i,j)
      END DO
    END DO

    CALL vcost(ny-1,nx-1,work,vwork2,ny,wsave2)

    DO j=1,ny-1    ! swapping arrays.
      DO i=1,nx-1
        temxy1(i,j) = work(j,i)
      END DO
    END DO

    CALL vcost(nx-1,ny-1,temxy1,vwork1,nx+1,wsave1)

  END IF           ! end of runopt if block, w (nz-1) is updated.

  RETURN
END SUBROUTINE uprad3
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                  SUBROUTINE PPRTBBC                  ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE pprtbbc(nx,ny,nz,g05,buoy2swt,rhostr,pprt,ptprt,             &
           pbari,ptbari,phydro,                                         &
           pprt1,tem1)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  To compute the perturbation pressure below the ground (scalar k=1)
!  The technique uses the perturbation hydrostatic relation obtained
!  from the model vertical momentum equation.  The technique is
!  summarized here.
!
!  The calculation of the hydrostatic pressure bottom boundary
!  condition involves two steps.  The first step is used to
!  approximate the second order pprt term in the buoyancy term
!  and the second step obtains the final pprt at k=1.
!
!  Step A: pprt(1)=D=(pprt(2) - F - dz*phydro + 2.0*A*B)/(1-C)
!
!  where,      A = dz*g*avgz(rhobar*J3)/(2*cpdcv)
!              B = pprt(i,j,2)/pbar(i,j,2)
!              C = A/pbar(i,j,1)
!              F = alpha*div*(1)-alpha*div*(2)  not included....
!  note:
!              D = the intermediate result pprt(1) for use in the
!                  second order pprt term
!
!  Step B:  pprt(1) = D + A*(1/cpdcv -1)*(B*B + E*E)/(1.0-C)
!
!  where,      E = D/pbar(i,j,1)
!
!  The result is stored in phydro(i,j) and passed into bcp to set
!  pprt(i,j,1).
!
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Dan Weber
!  11/06/1997.
!
!
!  Modification History:
!
!  3/18/99 (D. Weber)
!  Bug fix to second order terms.

!-----------------------------------------------------------------------
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
!    ptprt    Perturbation potential temperature at time tpresent (K)
!    pprt     Perturbation pressure at time tpresent (Pascal)
!
!    phydro   Big time step forcing term for use in computing the
!             hydrostatic pressure at k=1.
!
!    rhostr   Base state density rhobar times j3 (kg/m**3)
!    ptbari   Inverse base state potential temperature (K)
!    pbari    Inverse base state pressure (Pascal)
!
!  OUTPUT:
!
!    pprt1    Holds pprt(i,j,1) via hydrostatic balance.
!
!  WORK ARRAYS:
!
!    tem2     Work array.
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE             ! Force explicit declarations

  INTEGER :: nx, ny, nz        ! Number of grid points in 3 directions

  REAL :: ptprt (nx,ny,nz)     ! Perturbation potential temperature (K)
  REAL :: pprt  (nx,ny,nz)     ! Perturbation pressure at tfuture
                               ! (Pascal)

  REAL :: phydro(nx,ny)        ! Big time step forcing for computing
                               ! hydrostatic pprt at k=1.

  REAL :: rhostr(nx,ny,nz)     ! Base state density rhobar times j3.
  REAL :: ptbari(nx,ny,nz)     ! Inverse base state pot. temperature (K)
  REAL :: pbari (nx,ny,nz)     ! Inverse base state pressure (Pascal).
  REAL :: pprt1 (nx,ny)        ! pprt1

!
!-----------------------------------------------------------------------
!
!  Temporary WORK ARRAYS:
!
!-----------------------------------------------------------------------
!
  REAL :: tem1  (nx,ny,nz)     ! Temporary work array
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i, j
  REAL :: g05

  REAL :: tema,temb,a,b,c,d,e
  INTEGER :: buoy2swt !Switch for 1st-order or 2nd-order in buoyancy
  REAL :: ptemk,ptemk1,pttemk,pttemk1

!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'bndry.inc'
  INCLUDE 'globcst.inc'
  INCLUDE 'grid.inc'          ! Grid & map parameters.
  INCLUDE 'phycst.inc'      ! Physical constants
  INCLUDE 'mp.inc'            ! Message passing parameters.
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

  IF( ptsmlstp == 1 ) THEN  ! add in the 1st and 2nd order pt terms.
    tema = 0.5/cpdcv
    DO j=1,ny-1
      DO i=1,nx-1
        ptemk = pprt(i,j,2)*pbari(i,j,2)   ! new code.....
        ptemk1= pprt(i,j,1)*pbari(i,j,1)
        pttemk = ptprt(i,j,2)*ptbari(i,j,2)
        pttemk1= ptprt(i,j,1)*ptbari(i,j,1)
        temb = 0.5*(rhostr(i,j,1)+rhostr(i,j,2))
        pprt1(i,j)=phydro(i,j)+temb*g05*(pttemk+pttemk1                 &
                   -buoy2swt*(tema*(ptemk*pttemk+ptemk1*pttemk1)        &
                   +pttemk*pttemk+pttemk1*pttemk1))
      END DO
    END DO

  ELSE

    DO j=1,ny-1
      DO i=1,nx-1
        pprt1(i,j)= phydro(i,j)
      END DO
    END DO

  END IF

  tema =buoy2swt*0.5*(1.0/cpdcv -1.0)
  temb = dz*g*0.25/cpdcv
  DO j=1,ny-1
    DO i=1,nx-1
!  compute the intermediate result...
      a = temb*(rhostr(i,j,2)+rhostr(i,j,1))
      b = pprt(i,j,2)*pbari(i,j,2)
      c = a*pbari(i,j,1)
      d = (pprt(i,j,2) - dz*pprt1(i,j) + a*b) / (1.0-c)

!  compute the final pprt(i,j,1)....
      e = d*pbari(i,j,1)
      pprt1(i,j) = d + a*tema*(b*b + e*e) / (1.0-c)
    END DO
  END DO

  IF (mp_opt > 0) THEN
    CALL acct_interrupt(mp_acct)
    CALL mpsendrecv1dew(pprt1,nx,ny,ebc,wbc,0,tem1)
    CALL mpsendrecv1dns(pprt1,nx,ny,nbc,sbc,0,tem1)
  END IF
  CALL acct_interrupt(bc_acct)
  CALL bcs2d(nx,ny,pprt1, ebc,wbc,nbc,sbc)
  CALL acct_stop_inter

  RETURN
END SUBROUTINE pprtbbc

